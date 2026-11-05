#!/usr/bin/env bash
# Run a command inside a resource-capped Docker container while a host-side sidecar
# records exact cgroup-v2 metrics (peak RSS, CPU-seconds, I/O). This is THE fairness
# primitive: every benchmarked invocation goes through here with identical caps.
#
# Usage:
#   run_in_container.sh \
#     --name run_xyz --metrics-out m.json [--series-out m.csv] \
#     [--image ubuntu:24.04] [--cpus 30] [--cpuset 0-29] [--mem 500g] \
#     [--mount HOST:CONT[:ro]]... [--env K=V]... [--workdir DIR] \
#     [--meta K=V]... -- <command and args...>
#
# Notes:
#  * --memory-swap == --memory disables swap (true RAM ceiling).
#  * Same image/caps for Cell Ranger and simpleaf -> identical kernel + cgroup constraints.
#  * The command runs as the invoking uid:gid so outputs are owned by you.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

NAME=""; IMAGE="ubuntu:24.04"; CPUS="30"; CPUSET="0-29"; MEM="256g"
MOUNTS=(); ENVS=(); METRICS=""; SERIES=""; WORKDIR=""; META=()
while [[ $# -gt 0 ]]; do case "$1" in
  --name) NAME="$2"; shift 2;;
  --image) IMAGE="$2"; shift 2;;
  --cpus) CPUS="$2"; shift 2;;
  --cpuset) CPUSET="$2"; shift 2;;
  --mem) MEM="$2"; shift 2;;
  --mount) MOUNTS+=( -v "$2" ); shift 2;;
  --env) ENVS+=( -e "$2" ); shift 2;;
  --workdir) WORKDIR="$2"; shift 2;;
  --metrics-out) METRICS="$2"; shift 2;;
  --series-out) SERIES="$2"; shift 2;;
  --meta) META+=( "$2" ); shift 2;;
  --) shift; break;;
  *) echo "unknown arg: $1" >&2; exit 2;;
esac; done
[[ -z "$NAME" || -z "$METRICS" || $# -eq 0 ]] && { echo "missing --name/--metrics-out/command" >&2; exit 2; }
mkdir -p "$(dirname "$METRICS")"

docker rm -f "$NAME" >/dev/null 2>&1 || true

DRUN=( docker run --name "$NAME" --rm=false
       --cpus "$CPUS" --cpuset-cpus "$CPUSET" --memory "$MEM" --memory-swap "$MEM"
       --user "$(id -u):$(id -g)" --ulimit nofile=8192:8192 )
[[ -n "$WORKDIR" ]] && DRUN+=( -w "$WORKDIR" )
DRUN+=( "${MOUNTS[@]}" "${ENVS[@]}" "$IMAGE" "$@" )

# start metrics sidecar (waits for the container to appear, then polls its cgroup)
SERIES_ARG=(); [[ -n "$SERIES" ]] && SERIES_ARG=( --series "$SERIES" )
python3 "$HERE/collect_metrics.py" --container "$NAME" --out "$METRICS" \
        "${SERIES_ARG[@]}" --meta "${META[@]}" name="$NAME" image="$IMAGE" \
        cpus="$CPUS" cpuset="$CPUSET" mem="$MEM" cmd="$*" &
COLLECTOR=$!

set +e
"${DRUN[@]}"
RC=$?
set -e

wait "$COLLECTOR" 2>/dev/null || true
docker logs "$NAME" > "${METRICS%.json}.container.log" 2>&1 || true
docker rm -f "$NAME" >/dev/null 2>&1 || true
echo "run_in_container: $NAME exit=$RC metrics=$METRICS"
exit "$RC"
