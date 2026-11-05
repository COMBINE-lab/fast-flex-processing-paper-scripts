#!/usr/bin/env bash
# 00: local NVMe scratch, docker access, analysis env, host facts.
#   - requires a fast local filesystem mounted at the parent of $SCR
#   - checks docker-group access (needed by the timed container runs)
#   - checks the scbp env imports (python + R packages used by the downstream analysis)
#   - (re)builds the venv $EXT_ENV (scbp + celltypist/harmonypy) from
#     envs/scbp_ext/{freeze.txt,scbp_constraints.txt} via lib/ext_env.sh; only verifies it when complete
#   - the same for the venv $PROBE_ENV (scbp + pydeseq2 / celltypist; $PY_PROBE) from envs/probe/freeze.txt
#   - records host facts to $DUR/00_scratch_and_env/host.txt
set -euo pipefail
source "$(dirname "$0")/../config.env"
source "$REPRO/lib/common.sh"
source "$REPRO/lib/ext_env.sh"

MP="$(dirname "$SCR")"
mountpoint -q "$MP" || fail "$MP is not a mount point: mount a fast local filesystem there (the timings were measured on striped instance-store NVMe)"
mkdir -p "$SCR" && touch "$SCR/.w" && rm -f "$SCR/.w" || fail "$SCR not writable"
say "scratch OK: $(findmnt -n -o SOURCE,FSTYPE "$MP") $(df -h --output=avail "$MP" | tail -1) free"

id -nG | tr ' ' '\n' | grep -qx docker || fail "user is not in the docker group"
docker info --format '{{.ServerVersion}}' >/dev/null || fail "docker daemon not reachable"
say "docker OK: server $(docker info --format '{{.ServerVersion}}')"

"$PY" -c "import scanpy, anndata, igraph, leidenalg, numpy, pandas; print('python', scanpy.__version__, anndata.__version__)" \
  > "$LOGD/env_python.txt" || fail "scbp python imports failed ($PY)"
"$(dirname "$PY")/Rscript" -e 'for (p in c("SoupX","scDblFinder","scran","scry")) { suppressPackageStartupMessages(library(p, character.only=TRUE)); cat(p, as.character(packageVersion(p)), "\n") }' \
  > "$LOGD/env_R.txt" 2>&1 || { cat "$LOGD/env_R.txt"; fail "scbp R packages failed"; }
say "env OK: $(cat "$LOGD/env_python.txt"); R: $(paste -sd';' "$LOGD/env_R.txt")"
build_ext_env "$EXT_ENV" "$EXT_ENV/freeze.txt" "$EXT_ENV/scbp_constraints.txt" > "$LOGD/ext_env.log" 2>&1 \
  || { cat "$LOGD/ext_env.log"; fail "scbp_ext venv build/verify failed ($EXT_ENV)"; }
say "ext env OK: $EXT_ENV = freeze.txt ($(wc -l < "$EXT_ENV/freeze.txt") packages; $(cat "$EXT_ENV/.done"))"
build_probe_env "$PROBE_ENV" "$PROBE_ENV/freeze.txt" "$EXT_ENV/scbp_constraints.txt" > "$LOGD/probe_env.log" 2>&1 \
  || { cat "$LOGD/probe_env.log"; fail "probe venv build/verify failed ($PROBE_ENV)"; }
say "probe env OK: $PROBE_ENV = freeze.txt ($(wc -l < "$PROBE_ENV/freeze.txt") packages; $(cat "$PROBE_ENV/.done"))"

{
  echo "date_utc	$(date -u +%FT%TZ)"
  echo "hostname	$(hostname)"
  echo "cpu_model	$(lscpu | awk -F: '/^Model name/{gsub(/^ +/,"",$2); print $2}')"
  echo "sockets	$(lscpu | awk -F: '/^Socket\(s\)/{gsub(/ /,"",$2); print $2}')"
  echo "threads_per_core	$(lscpu | awk -F: '/^Thread\(s\) per core/{gsub(/ /,"",$2); print $2}')"
  echo "nproc	$(nproc)"
  echo "mem_total_gib	$(awk '/^MemTotal/{printf "%.1f", $2/1048576}' /proc/meminfo)"
  echo "kernel	$(uname -r)"
  echo "os	$(. /etc/os-release && echo "$PRETTY_NAME")"
  echo "scratch	$(findmnt -n -o SOURCE,FSTYPE,SIZE "$MP")"
  echo "docker	$(docker info --format '{{.ServerVersion}}')"
  echo "rustc	$(rustc --version)"
  echo "python_env	$(cat "$LOGD/env_python.txt")"
  echo "ext_env	$EXT_ENV freeze.txt sha256 $(sha256sum "$EXT_ENV/freeze.txt" | cut -c1-16)"
  echo "probe_env	$PROBE_ENV freeze.txt sha256 $(sha256sum "$PROBE_ENV/freeze.txt" | cut -c1-16)"
  sed 's/^/R	/' "$LOGD/env_R.txt"
} > "$DD/host.txt"
cat "$DD/host.txt"
mkdir -p "$SD"; date -u +%FT%TZ > "$SD/.done"   # scratch-side marker: run_all reruns 00 after an instance stop
write_done "host.txt	$(sha256sum "$DD/host.txt" | cut -c1-64)"
