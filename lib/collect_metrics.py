#!/usr/bin/env python3
"""
Poll a running Docker container's resource usage at a fixed interval and emit an
exact record: peak RSS, CPU-seconds (-> mean %CPU), and block I/O.

Portable across cgroup layouts (the source used is recorded in the JSON as mem_source/cpu_source):
  * cgroup v2 (preferred): memory.peak / memory.current / cpu.stat / io.stat read from
    INSIDE the container via `docker exec` (nested-container-safe; host-side docker-<id>.scope
    files read 0 when we run nested).
  * cgroup v1: memory.max_usage_in_bytes / cpuacct.usage (read inside via docker exec).
  * NO memory controller (where /sys/fs/cgroup is a bare tmpfs and
    only cpuacct is populated): HOST-side process-tree sampler -- sum VmHWM (peak RSS high-water)
    over the container's PID subtree (docker inspect .State.Pid -> /proc/<pid>) and track the peak.
    This is correct for multi-process tools (Cell Ranger / Martian fork many workers). CPU then
    comes from v1 cpuacct.usage if present, else summed /proc utime+stime over the subtree.

Run as a sidecar to `docker run`: samples until the container exits, then writes a JSON summary
(+ optional time-series CSV).

Usage:
  collect_metrics.py --container NAME --out run.json [--series run_series.csv]
                     [--interval 1.0] [--meta key=val ...]
Stdlib only.
"""
import argparse, json, subprocess, sys, time, os

DEVNULL = subprocess.DEVNULL
CLK_TCK = os.sysconf("SC_CLK_TCK")


def docker_inspect(name):
    try:
        out = subprocess.check_output(["docker", "inspect", name], text=True, stderr=DEVNULL)
        return json.loads(out)[0]
    except Exception:
        return None


def container_pid(name):
    """Host PID of the container's init process (0/None if not running)."""
    try:
        out = subprocess.check_output(["docker", "inspect", "-f", "{{.State.Pid}}", name],
                                      text=True, stderr=DEVNULL).strip()
        p = int(out)
        return p if p > 0 else None
    except Exception:
        return None


def cg_sample(name):
    """One `docker exec` reading cgroup stats from inside the container (v2 preferred, v1 fallback).
    Returns dict or None if the exec itself failed."""
    script = (
        "cat /sys/fs/cgroup/memory.peak 2>/dev/null; echo ===;"
        "cat /sys/fs/cgroup/memory.current 2>/dev/null; echo ===;"
        "grep usage_usec /sys/fs/cgroup/cpu.stat 2>/dev/null; echo ===;"
        "cat /sys/fs/cgroup/io.stat 2>/dev/null; echo ===;"
        "grep -E '^(anon|file) ' /sys/fs/cgroup/memory.stat 2>/dev/null; echo ===;"
        "cat /sys/fs/cgroup/memory/memory.max_usage_in_bytes 2>/dev/null; echo ===;"
        "cat /sys/fs/cgroup/memory/memory.usage_in_bytes 2>/dev/null; echo ===;"
        "cat /sys/fs/cgroup/cpu,cpuacct/cpuacct.usage 2>/dev/null"
    )
    try:
        out = subprocess.check_output(["docker", "exec", name, "sh", "-c", script],
                                      text=True, stderr=DEVNULL)
    except Exception:
        return None
    parts = out.split("===")

    def fi(i):
        try:
            return int(parts[i].split()[0])
        except Exception:
            return None

    d = {"mem_peak_v2": fi(0), "mem_cur_v2": fi(1)}
    cpu = None
    try:
        toks = parts[2].split()
        if len(toks) >= 2 and toks[0] == "usage_usec":
            cpu = int(toks[1])
    except Exception:
        pass
    d["cpu_usec_v2"] = cpu
    r = w = 0
    try:
        for tok in parts[3].split():
            if tok.startswith("rbytes="):
                r += int(tok[7:])
            elif tok.startswith("wbytes="):
                w += int(tok[7:])
    except Exception:
        pass
    d["io_r"], d["io_w"] = r, w
    anon = file = None
    try:
        toks = parts[4].split()
        for i in range(0, len(toks) - 1, 2):
            if toks[i] == "anon":
                anon = int(toks[i + 1])
            elif toks[i] == "file":
                file = int(toks[i + 1])
    except Exception:
        pass
    d["anon"], d["file"] = anon, file
    d["mem_max_v1"] = fi(5)
    d["mem_cur_v1"] = fi(6)
    d["cpuacct_ns_v1"] = fi(7)
    return d


def proc_tree(root):
    """Set of host PIDs = root + all descendants (via /proc/*/stat ppid)."""
    children = {}
    for e in os.listdir("/proc"):
        if not e.isdigit():
            continue
        try:
            with open("/proc/%s/stat" % e) as f:
                after = f.read().rsplit(") ", 1)[1].split()  # fields after "comm)"; robust to spaces in comm
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


def proc_rss_ticks(pids):
    """Sum VmHWM (bytes) and CPU ticks (utime+stime) over the given host PIDs.
    VmHWM = 'peak resident set size (high-water mark)' from /proc/<pid>/status — the kernel's
    EXPLICIT peak-RSS field, so the peak is never missed between samples (unlike VmRSS, which is
    the *current* resident). This is exactly the field Nextflow's trace `peak_rss` reads (the
    alevin-fry paper's memory metric). Returns (peak_rss_bytes or None, ticks)."""
    rss = 0
    ticks = 0
    ok = False
    for p in pids:
        try:
            with open("/proc/%d/status" % p) as f:
                for line in f:
                    if line.startswith("VmHWM:"):
                        rss += int(line.split()[1]) * 1024
                        ok = True
                        break
        except Exception:
            pass
        try:
            with open("/proc/%d/stat" % p) as f:
                after = f.read().rsplit(") ", 1)[1].split()
            ticks += int(after[11]) + int(after[12])  # utime + stime
        except Exception:
            pass
    return (rss if ok else None), ticks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--container", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--series")
    ap.add_argument("--interval", type=float, default=1.0)
    ap.add_argument("--meta", nargs="*", default=[])
    args = ap.parse_args()

    meta = {}
    for kv in args.meta:
        if "=" in kv:
            k, v = kv.split("=", 1)
            meta[k] = v

    # wait for the container to be running AND observable (either cgroup exec works, or we can
    # see its host PID for the /proc sampler). Either is enough to measure.
    pid = None
    t_deadline = time.monotonic() + 300
    while time.monotonic() < t_deadline:
        info = docker_inspect(args.container)
        running = bool(info and info.get("State", {}).get("Running"))
        if running:
            pid = container_pid(args.container)
            if pid is not None or cg_sample(args.container) is not None:
                break
        time.sleep(0.2)
    else:
        json.dump({"error": "container not observable", **meta}, open(args.out, "w"), indent=2)
        sys.exit("could not observe container " + args.container)

    series = open(args.series, "w") if args.series else None
    if series:
        series.write("t_s,rss_bytes,cpu_usec,io_r_bytes,io_w_bytes,mem_source\n")

    t0 = time.monotonic()
    peak_v2 = peak_v1 = peak_proc = 0          # peak RSS candidates (bytes)
    # CPU baselines (captured from the first successful read of each source)
    cpu_usec0 = cpuacct0 = proc_ticks0 = None
    last_cpu_usec = last_cpuacct = last_proc_ticks = None
    ior0 = iow0 = None
    last_r = last_w = 0

    while True:
        info = docker_inspect(args.container)
        running = bool(info and info.get("State", {}).get("Running"))
        cg = cg_sample(args.container) if running else None
        # host-side process-tree RSS (works with no memory controller)
        prss = None; pticks = None
        if running and pid is not None:
            tree = proc_tree(pid)
            prss, pticks = proc_rss_ticks(tree)

        now = time.monotonic() - t0
        rss_now = None; msrc = None
        if cg:
            if cg.get("mem_peak_v2"):
                peak_v2 = max(peak_v2, cg["mem_peak_v2"]); rss_now = cg.get("mem_cur_v2"); msrc = "cgroup_v2"
            if cg.get("mem_max_v1"):
                peak_v1 = max(peak_v1, cg["mem_max_v1"])
                if rss_now is None: rss_now = cg.get("mem_cur_v1"); msrc = "cgroup_v1"
            # CPU accumulators from cgroup
            if cg.get("cpu_usec_v2") is not None:
                if cpu_usec0 is None: cpu_usec0 = cg["cpu_usec_v2"]
                last_cpu_usec = cg["cpu_usec_v2"]
            if cg.get("cpuacct_ns_v1") is not None:
                if cpuacct0 is None: cpuacct0 = cg["cpuacct_ns_v1"]
                last_cpuacct = cg["cpuacct_ns_v1"]
            if cg.get("io_r") is not None:
                if ior0 is None: ior0, iow0 = cg["io_r"], cg["io_w"]
                last_r, last_w = max(last_r, cg["io_r"]), max(last_w, cg["io_w"])
        if prss is not None:
            peak_proc = max(peak_proc, prss)
            if rss_now is None: rss_now = prss; msrc = "proc_rss"
        if pticks is not None:
            if proc_ticks0 is None: proc_ticks0 = pticks
            last_proc_ticks = pticks

        if series and rss_now is not None:
            cu = last_cpu_usec if last_cpu_usec is not None else (
                int(last_cpuacct / 1000) if last_cpuacct is not None else "")
            series.write(f"{now:.2f},{rss_now},{cu},{last_r},{last_w},{msrc}\n"); series.flush()

        if not running:
            break
        time.sleep(args.interval)

    wall = time.monotonic() - t0
    info = docker_inspect(args.container) or {}
    state = info.get("State", {})

    # PEAK RSS = peak resident memory over the run — every source below is an EXPLICIT peak field:
    # cgroup memory.peak (v2) / max_usage_in_bytes (v1), else summed VmHWM (peak high-water) over
    # the container's process tree. Same peak-RSS metric Nextflow's trace `peak_rss` reports —
    # the alevin-fry paper's method (notebooks/plot-time-rss.ipynb reads the trace `peak_rss`).
    if peak_v2 > 0:
        peak_rss = peak_v2; mem_source = "cgroup_v2:memory.peak"
    elif peak_v1 > 0:
        peak_rss = peak_v1; mem_source = "cgroup_v1:memory.max_usage_in_bytes"
    else:
        peak_rss = peak_proc; mem_source = "host:proc_tree_vmhwm_sum"

    # choose CPU seconds: prefer cgroup v2 cpu.stat, then v1 cpuacct, then /proc ticks
    cpu_sec = None; cpu_source = None
    if cpu_usec0 is not None and last_cpu_usec is not None:
        cpu_sec = (last_cpu_usec - cpu_usec0) / 1e6; cpu_source = "cgroup_v2:cpu.stat"
    elif cpuacct0 is not None and last_cpuacct is not None:
        cpu_sec = (last_cpuacct - cpuacct0) / 1e9; cpu_source = "cgroup_v1:cpuacct.usage"
    elif proc_ticks0 is not None and last_proc_ticks is not None:
        cpu_sec = (last_proc_ticks - proc_ticks0) / CLK_TCK; cpu_source = "host:proc_utime_stime"

    rec = {
        **meta,
        "container": args.container,
        "exit_code": state.get("ExitCode"),
        "wall_s": round(wall, 3),
        "peak_rss_bytes": peak_rss,
        "peak_rss_gb": round(peak_rss / 2**30, 3),
        "mem_source": mem_source,
        "cpu_seconds": None if cpu_sec is None else round(cpu_sec, 1),
        "cpu_source": cpu_source,
        "mean_cpu_pct": None if (cpu_sec is None or wall == 0) else round(100 * cpu_sec / wall, 1),
        "io_read_bytes": (last_r - ior0) if ior0 is not None else None,
        "io_write_bytes": (last_w - iow0) if iow0 is not None else None,
        "io_read_gb": round((last_r - ior0) / 2**30, 3) if ior0 is not None else None,
        "io_write_gb": round((last_w - iow0) / 2**30, 3) if iow0 is not None else None,
        "docker_started_at": state.get("StartedAt"),
        "docker_finished_at": state.get("FinishedAt"),
    }
    if series:
        series.close()
    json.dump(rec, open(args.out, "w"), indent=2)
    print(json.dumps(rec, indent=2))


if __name__ == "__main__":
    main()
