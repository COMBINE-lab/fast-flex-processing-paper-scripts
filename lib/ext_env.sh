# build_venv DIR FREEZE CONSTRAINTS WANT CHECK: (re)create a venv over the scbp env reproducibly from
# envs/<name>/freeze.txt and envs/scbp_ext/scbp_constraints.txt. Source AFTER config.env.
#   - the scbp env itself must still equal CONSTRAINTS (`pip list --format=freeze`), else the freeze is stale
#   - venv over scbp (--system-site-packages); the venv-local packages are exactly the `==` pins of FREEZE
#     that scbp does not already provide at that version, installed --no-deps (no resolver drift)
#   - verified: `pip freeze` == FREEZE byte for byte, no venv package shadows an scbp one, and the python
#     snippet CHECK runs; then DIR/.done = WANT (the string the steps check)
# A complete DIR (.done == WANT + bin/python) is only re-verified, never rebuilt.
#   build_ext_env DIR FREEZE CONSTRAINTS   = envs/scbp_ext (celltypist + harmonypy; steps 10/11)
#   build_probe_env DIR FREEZE CONSTRAINTS = envs/probe (pydeseq2 + celltypist; $PY_PROBE, steps 12/14/15)
build_venv(){
  local dir="$1" freeze="$2" cons="$3" want="$4" check="$5" pins shadow
  [[ -s "$freeze" && -s "$cons" ]] || { echo "build_venv: missing $freeze or $cons" >&2; return 1; }
  diff -q <("$PY" -m pip list --format=freeze 2>/dev/null) "$cons" > /dev/null \
    || { echo "build_venv: the scbp env ($PY) no longer equals $cons" >&2; return 1; }
  if [[ ! -x "$dir/bin/python" || "$(cat "$dir/.done" 2>/dev/null)" != "$want" ]]; then
    pins=$(grep -E '^[A-Za-z0-9_.-]+==' "$freeze" | grep -vxFf "$cons")
    [[ -n "$pins" ]] || { echo "build_venv: no venv-local pins in $freeze" >&2; return 1; }
    rm -f "$dir/.done"; rm -rf "$dir/bin" "$dir/lib" "$dir/lib64" "$dir/include" "$dir/pyvenv.cfg"
    "$PY" -m venv --system-site-packages "$dir"
    # shellcheck disable=SC2086
    "$dir/bin/python" -m pip install -q --disable-pip-version-check --no-deps $pins
  fi
  diff <("$dir/bin/python" -m pip freeze 2>/dev/null) "$freeze" \
    || { echo "build_venv: pip freeze of $dir differs from $freeze (diff above)" >&2; return 1; }
  shadow=$(comm -12 <(cut -d= -f1 "$cons" | tr 'A-Z_' 'a-z-' | sort -u) \
    <("$dir/bin/python" -m pip list --local --format=freeze 2>/dev/null | cut -d= -f1 | tr 'A-Z_' 'a-z-' | grep -vx 'pip\|setuptools' | sort -u))
  [[ -z "$shadow" ]] || { echo "build_venv: venv shadows scbp packages: $shadow" >&2; return 1; }
  "$dir/bin/python" -c "$check" || { echo "build_venv: import check failed ($dir)" >&2; return 1; }
  echo "$want" > "$dir/.done"
}
build_ext_env(){
  build_venv "$1" "$2" "$3" "celltypist==$CELLTYPIST_VERSION harmonypy==$HARMONYPY_VERSION" "import celltypist, harmonypy, scanpy, openpyxl
from importlib.metadata import version as v
assert celltypist.__version__ == '$CELLTYPIST_VERSION', celltypist.__version__
assert v('harmonypy') == '$HARMONYPY_VERSION', v('harmonypy')"
}
build_probe_env(){
  build_venv "$1" "$2" "$3" "pydeseq2==$PYDESEQ2_VERSION celltypist==$CELLTYPIST_VERSION" "import pydeseq2, statsmodels, celltypist, scanpy
assert pydeseq2.__version__ == '$PYDESEQ2_VERSION', pydeseq2.__version__
assert celltypist.__version__ == '$CELLTYPIST_VERSION', celltypist.__version__"
}
