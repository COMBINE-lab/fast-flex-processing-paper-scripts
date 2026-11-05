# Shared helpers for steps/*.sh. Source AFTER config.env. Defines STEP, say/fail, write_done, and for a
# step (steps/NN[x]_<slug>.sh, steps/run_all.sh, or an explicit STEP_NAME of that form) SD (scratch dir),
# DD (durable dir) and LOGD, creating $DD and $LOGD.
# Any other caller (a test or a shell sourcing this for its helpers) gets the helpers only: no directory
# is created under $DUR / $SCR and SD / DD / LOGD stay unset, so using them fails under set -u.
STEP="${STEP_NAME:-$(basename "$0" .sh)}"
if [[ "$STEP" =~ ^([0-9]{2}[a-z]?_[a-z0-9_]+|run_all)$ ]]; then
  SD="$SCR/$STEP"; DD="$DUR/$STEP"; LOGD="$DD/logs"
  mkdir -p "$DD" "$LOGD"
else
  unset SD DD LOGD
fi

say(){ echo "[$(date -u +%FT%TZ)] [$STEP] $*"; }
fail(){ say "FATAL: $*" >&2; exit 1; }

# write_done "name<TAB>checksum" ... -> $DD/DONE with date, pinned commits and one line per deliverable
write_done(){
  {
    echo "date	$(date -u +%FT%TZ)"
    echo "simpleaf	$SIMPLEAF_COMMIT"; echo "alevin-fry	$ALEVINFRY_COMMIT"
    echo "piscem	$PISCEM_COMMIT (piscem-rs $PISCEM_RS_VERSION)"
    echo "full	${FULL:-0}"   # run_all --full: every optional analysis was produced
    printf '%s\n' "$@"
  } > "$DD/DONE"
  say "DONE -> $DD/DONE"
}
