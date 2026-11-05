# Marker validity rules shared by the steps and run_all.sh (source AFTER config.env).
# export_marker_current MARKER DS: 0 iff MARKER is a valid step-07 export completion under the CURRENT
#   config. A "skipped (not in PROBE_EXPORT_DATASETS)" marker is valid only while DS is still outside
#   $PROBE_EXPORT_DATASETS; once DS is added, the export is due and the marker is stale.
export_marker_current(){
  [[ -f "$1" ]] || return 1
  if [[ "$(head -1 "$1")" == "skipped (not in PROBE_EXPORT_DATASETS)" ]]; then
    [[ " $PROBE_EXPORT_DATASETS " != *" $2 "* ]]
  fi
}
