#!/usr/bin/env python3
"""
Host-side peak-RSS sampler for a Cell Ranger `multi` run.

WHY THIS EXISTS
  Cell Ranger's Martian runtime (`mrp`) forks stage workers that double-fork and
  REPARENT away from mrp (to PID 1, or to the launching shell). A tree-walk rooted at
  the container-init / mrp pid — what lib/collect_metrics.py does — therefore finds only
  mrp itself (~0.04 GB) and misses every heavy worker. Observed live (readlevel CR):
  mrp VmHWM=40 MB while ~20 `cr_lib` workers each held ~1.3 GB, cwd'd deep under the
  pipestance dir. So the standard sampler logs peak_rss = null / near-zero for CR.

WHAT THIS DOES
  Every `--interval` seconds, enumerate ALL host processes and select those belonging to
  the CR run by a REPARENT-PROOF union of signals:
    (a) /proc/<pid>/cwd resolves to (or under) the pipestance dir, OR
    (b) the run-id or pipestance path appears in /proc/<pid>/cmdline, OR
    (c) [container mode] the pid shares the container init's PID namespace
        (/proc/<pid>/ns/pid == container-init's ns/pid).
  For each selected LIVE pid read VmHWM (peak high-water) and VmRSS (current) from
  /proc/<pid>/status. peak_rss = max over samples of (Σ VmHWM of live members) — the
  IDENTICAL metric semantics as collect_metrics.py's host:proc_tree_vmhwm_sum and the
  alevin-fry paper's Nextflow peak_rss. Σ VmRSS peak is recorded as a secondary field.

  In container mode it ALSO tracks collect_metrics.py's proc-tree-from-init value as a
  side-by-side sanity key (peak_rss_proctree_init_gb) to make the reparenting effect
  explicit in the small-dataset validation.

TERMINATION
  --duration N     : sample for N seconds, then exit (used for the read-only live
                     validation against an already-running CR).
  --container NAME : sample until `docker inspect` reports the container not Running.
  --stop-file PATH : sample until PATH exists.
  At least one terminator must be given.

Stdlib only. Pure /proc reads + optional `docker inspect`; root-free, negligible load.
"""
import argparse, json, os, subprocess, sys, time

DEVNULL = subprocess.DEVNULL


def read_cmdline(pid):
    try:
        with open("/proc/%d/cmdline" % pid, "rb") as f:
            return f.read().replace(b"\0", b" ").decode("utf-8", "replace")
    except Exception:
        return ""


def read_cwd(pid):
    try:
        return os.readlink("/proc/%d/cwd" % pid)
    except Exception:
        return ""


def read_hwm_rss(pid):
    """(VmHWM bytes or None, VmRSS bytes or None) from /proc/<pid>/status."""
    hwm = rss = None
    try:
        with open("/proc/%d/status" % pid) as f:
            for line in f:
                if line.startswith("VmHWM:"):
                    hwm = int(line.split()[1]) * 1024
                elif line.startswith("VmRSS:"):
                    rss = int(line.split()[1]) * 1024
                if hwm is not None and rss is not None:
                    break
    except Exception:
        pass
    return hwm, rss


def pid_ns(pid):
    try:
        return os.readlink("/proc/%d/ns/pid" % pid)
    except Exception:
        return None


def docker_init_pid(name):
    try:
        out = subprocess.check_output(
            ["docker", "inspect", "-f", "{{.State.Pid}}", name],
            text=True, stderr=DEVNULL).strip()
        p = int(out)
        return p if p > 0 else None
    except Exception:
        return None


def container_running(name):
    try:
        out = subprocess.check_output(
            ["docker", "inspect", "-f", "{{.State.Running}}", name],
            text=True, stderr=DEVNULL).strip()
        return out == "true"
    except Exception:
        return False


def proc_tree_from(root):
    """collect_metrics.py-style descendant set from a root pid (for side-by-side)."""
    children = {}
    for e in os.listdir("/proc"):
        if not e.isdigit():
            continue
        try:
            with open("/proc/%s/stat" % e) as f:
                after = f.read().rsplit(") ", 1)[1].split()
            ppid = int(after[1])
        except Exception:
            continue
        children.setdefault(ppid, []).append(int(e))
    seen, stack = set(), [root]
    while stack:
        p = stack.pop()
        if p in seen:
            continue
        seen.add(p)
        stack.extend(children.get(p, []))
    return seen


def all_pids():
    return [int(e) for e in os.listdir("/proc") if e.isdigit()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pipestance", required=True,
                    help="CR --id output dir (workers cwd under here)")
    ap.add_argument("--run-id", default="",
                    help="unique run id / substring to match in cmdline")
    ap.add_argument("--out", required=True)
    ap.add_argument("--series")
    ap.add_argument("--interval", type=float, default=1.0)
    ap.add_argument("--container", default="",
                    help="docker container name: adds PID-ns membership + proc-tree side-by-side, and terminates when it stops")
    ap.add_argument("--duration", type=float, default=0.0,
                    help="sample this many seconds then exit (read-only validation)")
    ap.add_argument("--stop-file", default="")
    ap.add_argument("--wait-container-s", type=float, default=0.0,
                    help="in container mode, wait up to this long for the container to appear before sampling")
    ap.add_argument("--meta", nargs="*", default=[])
    args = ap.parse_args()

    if not (args.duration > 0 or args.container or args.stop_file):
        sys.exit("need at least one terminator: --duration, --container, or --stop-file")

    meta = {}
    for kv in args.meta:
        if "=" in kv:
            k, v = kv.split("=", 1)
            meta[k] = v

    pipe = os.path.realpath(args.pipestance)
    pipe_prefix = pipe + "/"
    run_id = args.run_id

    # container mode: optionally wait for it to appear, capture its init pid + pid-ns
    init_pid = None
    init_ns = None
    container_seen = False
    if args.container and args.wait_container_s > 0:
        deadline = time.monotonic() + args.wait_container_s
        while time.monotonic() < deadline:
            if container_running(args.container):
                container_seen = True
                break
            time.sleep(1.0)
    if args.container:
        init_pid = docker_init_pid(args.container)
        if init_pid is not None:
            container_seen = True
            init_ns = pid_ns(init_pid)
    # If the container never appeared within the wait window there is nothing to sample:
    # record the miss and exit immediately, so the caller's `wait` never hangs (RULE #5:
    # a container that never launched is a real failure, not something to silently poll on).
    if args.container and not container_seen:
        rec = {**meta, "pipestance": pipe, "run_id": run_id,
               "container": args.container, "mem_source": "host:cwd_cmdline_vmhwm_sum",
               "peak_rss_bytes": None, "peak_rss_gb": None, "peak_vmrss_gb": None,
               "n_pids_at_peak": 0, "t_at_peak_s": 0.0, "n_samples": 0, "wall_sampled_s": 0.0,
               "peak_rss_proctree_init_gb": None,
               "error": "container %r never appeared within %.0fs" % (args.container, args.wait_container_s)}
        with open(args.out, "w") as f:
            json.dump(rec, f, indent=2)
        print(json.dumps(rec, indent=2))
        sys.exit(3)

    def selected(pid):
        # (a) cwd under pipestance
        cwd = read_cwd(pid)
        if cwd and (cwd == pipe or cwd.startswith(pipe_prefix)):
            return True
        # (b) cmdline mentions run-id or pipestance path
        cl = read_cmdline(pid)
        if run_id and run_id in cl:
            return True
        if pipe in cl:
            return True
        # (c) container PID-namespace membership
        if init_ns is not None:
            ns = pid_ns(pid)
            if ns is not None and ns == init_ns:
                return True
        return False

    series = open(args.series, "w") if args.series else None
    if series:
        series.write("t_s,n_pids,sum_vmhwm_bytes,sum_vmrss_bytes,proctree_init_vmhwm_bytes\n")

    peak_hwm = 0
    peak_vmrss = 0
    n_at_peak = 0
    t_at_peak = 0.0
    peak_proctree = 0
    n_samples = 0
    t0 = time.monotonic()

    while True:
        now = time.monotonic() - t0
        sum_hwm = 0
        sum_rss = 0
        npids = 0
        for pid in all_pids():
            try:
                if not selected(pid):
                    continue
                hwm, rss = read_hwm_rss(pid)
                if hwm is not None:
                    sum_hwm += hwm
                    npids += 1
                if rss is not None:
                    sum_rss += rss
            except Exception:
                # per-pid failure (race: pid exited mid-scan) must never abort a sample
                continue

        # side-by-side: collect_metrics.py tree-walk from container init (container mode)
        proctree_hwm = 0
        if init_pid is not None:
            try:
                for pid in proc_tree_from(init_pid):
                    hwm, _ = read_hwm_rss(pid)
                    if hwm is not None:
                        proctree_hwm += hwm
            except Exception:
                proctree_hwm = 0
            peak_proctree = max(peak_proctree, proctree_hwm)

        if sum_hwm > peak_hwm:
            peak_hwm = sum_hwm
            n_at_peak = npids
            t_at_peak = now
        peak_vmrss = max(peak_vmrss, sum_rss)
        n_samples += 1

        if series:
            series.write("%.2f,%d,%d,%d,%d\n" % (now, npids, sum_hwm, sum_rss, proctree_hwm))
            series.flush()

        # termination checks
        if args.duration > 0 and now >= args.duration:
            break
        if args.stop_file and os.path.exists(args.stop_file):
            break
        if args.container:
            # if it appeared and is now gone, stop. Guard against init_pid never resolving:
            # the authoritative liveness signal is container_running, checked regardless.
            if init_pid is None:
                init_pid = docker_init_pid(args.container)
                if init_pid is not None:
                    init_ns = pid_ns(init_pid)
            if not container_running(args.container):
                break
        time.sleep(args.interval)

    wall = time.monotonic() - t0
    rec = {
        **meta,
        "pipestance": pipe,
        "run_id": run_id,
        "container": args.container or None,
        "mem_source": "host:cwd_cmdline_vmhwm_sum",
        "peak_rss_bytes": peak_hwm if peak_hwm > 0 else None,
        "peak_rss_gb": round(peak_hwm / 2**30, 3) if peak_hwm > 0 else None,
        "peak_vmrss_gb": round(peak_vmrss / 2**30, 3) if peak_vmrss > 0 else None,
        "n_pids_at_peak": n_at_peak,
        "t_at_peak_s": round(t_at_peak, 2),
        "n_samples": n_samples,
        "wall_sampled_s": round(wall, 2),
        "peak_rss_proctree_init_gb": (round(peak_proctree / 2**30, 3)
                                      if init_pid is not None else None),
    }
    if series:
        series.close()
    with open(args.out, "w") as f:
        json.dump(rec, f, indent=2)
    print(json.dumps(rec, indent=2))


if __name__ == "__main__":
    main()
