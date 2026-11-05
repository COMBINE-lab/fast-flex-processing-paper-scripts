#!/usr/bin/env bash
# 01: build the pinned toolchain from fresh clones and register it with simpleaf.
#   - clones simpleaf / alevin-fry / piscem into $SRC at the commits in config.env
#   - cargo build --release --locked with RUSTFLAGS=$RUSTFLAGS_BUILD (per-tool .done marker)
#   - self-verifies commits, clean worktrees, --version strings, Cargo.lock pins
#     (piscem -> piscem-rs $PISCEM_RS_VERSION, alevin-fry -> libradicl $LIBRADICL_MINOR.x),
#     the linked piscem-rs crate, af `quant --umi-edit-dist`, and the multiplex-quant flags
#   - runs `simpleaf set-paths` into $AF_HOME with the built alevin-fry and piscem (multiplex-quant forwards
#     --umi-edit-dist itself since simpleaf PR #224); copies simpleaf_info.json to $DUR
#   - clones cyto ($CYTO_REPO, tag $CYTO_TAG = $CYTO_COMMIT) into $SRC/cyto and builds it --locked with the
#     $CYTO_LOCK (cyto ships no Cargo.lock; ours pins the dependency versions of the cyto 0.4.6
#     build); records cyto's uv tools (cell-filter, pycyto: cyto installs them itself)
#   - records the pre-installed Cell Ranger $CELLRANGER_VERSION binary (version + sha256; README "Host requirements")
# Builds take ~10-20 min in total: run detached (setsid nohup) and poll the log.
set -euo pipefail
source "$(dirname "$0")/../config.env"
source "$REPRO/lib/common.sh"
mkdir -p "$SRC"

SA="$SRC/simpleaf/target/release/simpleaf"
AF="$SRC/alevin-fry/target/release/alevin-fry"
PISCEM="$SRC/piscem/target/release/piscem"

build(){ # name repo commit [Cargo.lock to install before the build]
  local name="$1" repo="$2" commit="$3" lock="${4:-}" dir="$SRC/$1"
  if [[ -f "$dir/.build.done" && "$(cat "$dir/.build.done")" == "$(git -C "$dir" rev-parse HEAD)" \
        && -x "$dir/target/release/$name" ]] && { [[ -z "$lock" ]] || cmp -s "$lock" "$dir/Cargo.lock"; }; then
    say "$name @ $commit already built -> skip"; return 0
  fi
  if [[ ! -d "$dir/.git" ]]; then
    say "clone $repo -> $dir"; git clone -q "$repo" "$dir"
  fi
  git -C "$dir" fetch -q origin
  git -C "$dir" checkout -q --detach "$commit"
  if [[ -n "$lock" ]]; then cp "$lock" "$dir/Cargo.lock"; fi
  say "build $name @ $(git -C "$dir" rev-parse --short HEAD) (RUSTFLAGS=$RUSTFLAGS_BUILD)"
  ( cd "$dir" && RUSTFLAGS="$RUSTFLAGS_BUILD" cargo build --release --locked -j "$CARGO_JOBS" ) > "$LOGD/build_$name.log" 2>&1 \
    || { tail -30 "$LOGD/build_$name.log"; fail "cargo build $name failed (see $LOGD/build_$name.log)"; }
  git -C "$dir" rev-parse HEAD > "$dir/.build.done"
}
build piscem     "$PISCEM_REPO"    "$PISCEM_COMMIT"
build alevin-fry "$ALEVINFRY_REPO" "$ALEVINFRY_COMMIT"
build simpleaf   "$SIMPLEAF_REPO"  "$SIMPLEAF_COMMIT"
build cyto       "$CYTO_REPO"      "$CYTO_COMMIT" "$CYTO_LOCK"

# ---------------- self-verify ----------------
for t in simpleaf:$SIMPLEAF_COMMIT alevin-fry:$ALEVINFRY_COMMIT piscem:$PISCEM_COMMIT cyto:$CYTO_COMMIT; do
  d="$SRC/${t%%:*}"; h=$(git -C "$d" rev-parse HEAD)
  [[ "$h" == "${t#*:}"* ]] || fail "${t%%:*} HEAD $h != ${t#*:}"
  [[ -z "$(git -C "$d" status --porcelain --untracked-files=no)" ]] || fail "${t%%:*} worktree is dirty"
done
[[ "$(git -C "$SRC/cyto" describe --tags --exact-match HEAD)" == "$CYTO_TAG" ]] || fail "cyto HEAD is not tag $CYTO_TAG"
cmp -s "$CYTO_LOCK" "$SRC/cyto/Cargo.lock" || fail "cyto Cargo.lock differs from $CYTO_LOCK"
mkdir -p "$AF_HOME"
SAVER=$(ALEVIN_FRY_HOME="$AF_HOME" "$SA" --version 2>&1 | tail -1)
[[ "$SAVER" == "simpleaf $SIMPLEAF_VERSION" ]] || fail "simpleaf --version is '$SAVER'"
[[ "$("$AF" --version 2>&1)" == "alevin-fry $ALEVINFRY_VERSION" ]] || fail "alevin-fry --version is '$("$AF" --version 2>&1)'"
[[ "$("$PISCEM" --version 2>&1 | head -1)" == "piscem $PISCEM_VERSION" ]] || fail "piscem --version is '$("$PISCEM" --version 2>&1 | head -1)'"
PRS=$(grep -A1 '^name = "piscem-rs"$' "$SRC/piscem/Cargo.lock" | grep -oP '(?<=version = ")[^"]+')
[[ "$PRS" == "$PISCEM_RS_VERSION" ]] || fail "piscem Cargo.lock pins piscem-rs '$PRS' != $PISCEM_RS_VERSION"
grep -q "piscem-rs-$PISCEM_RS_VERSION" "$SRC"/piscem/target/release/deps/piscem_rs-*.d \
  || fail "piscem binary not built against piscem-rs $PISCEM_RS_VERSION"
LRC=$(grep -A1 '^name = "libradicl"$' "$SRC/alevin-fry/Cargo.lock" | grep -oP '(?<=version = ")[^"]+')
[[ "$LRC" == "$LIBRADICL_MINOR".* ]] || fail "alevin-fry Cargo.lock pins libradicl '$LRC' (want $LIBRADICL_MINOR.x)"
git -C "$SRC/alevin-fry" merge-base --is-ancestor "$AF_PR196_ANCESTOR" HEAD || fail "alevin-fry PR #196 ($AF_PR196_ANCESTOR) not an ancestor"
# --umi-edit-dist is hidden from `quant --help`, so probe the parser: a known-bogus flag
# must be rejected as "unexpected argument" (control), --umi-edit-dist must not be.
# (outputs captured first: under pipefail, `af ... | grep -q` returns af's non-zero exit even on a match)
CTL=$("$AF" quant --no-such-flag-xyz 1 2>&1 || true)
grep -q "unexpected argument '--no-such-flag-xyz'" <<<"$CTL" || fail "alevin-fry quant does not reject an unknown flag (probe control broken): $CTL"
UED=$("$AF" quant --umi-edit-dist "$UMI_EDIT_DIST" 2>&1 || true)
grep -qE "unexpected argument|unrecognized" <<<"$UED" && fail "alevin-fry quant rejects --umi-edit-dist: $UED"
grep -q "required arguments were not provided" <<<"$UED" || fail "unexpected af quant --umi-edit-dist probe output: $UED"
say "alevin-fry quant accepts --umi-edit-dist (probe: '$(head -1 <<<"$UED")')"
MQ=$(ALEVIN_FRY_HOME="$AF_HOME" "$SA" multiplex-quant --help 2>&1)
for f in --cell-bc-correction --sample-bc-correction --excluded-probes --kmer-length --compress --decoder --thread-policy --umi-edit-dist; do
  grep -q -- "$f" <<<"$MQ" || fail "multiplex-quant --help lacks $f"
done
IX=$(ALEVIN_FRY_HOME="$AF_HOME" "$SA" index --help 2>&1)
for f in --probe-csv --excluded-probes --minimizer-length --kmer-length; do
  grep -q -- "$f" <<<"$IX" || fail "simpleaf index --help lacks $f"
done
say "self-verify OK: $SAVER ; alevin-fry $ALEVINFRY_VERSION (libradicl $LRC) ; piscem $PISCEM_VERSION (piscem-rs $PRS)"

# ---------------- set-paths ----------------
ALEVIN_FRY_HOME="$AF_HOME" "$SA" set-paths -p "$PISCEM" -a "$AF" > "$LOGD/set_paths.log" 2>&1
INFO="$AF_HOME/simpleaf_info.json"
read -r RECAF RECPI < <("$PY" -c "import json,sys;p=json.load(open(sys.argv[1]))['prog_info'];print(p['alevin_fry']['exe_path'],p['piscem']['exe_path'])" "$INFO")
[[ "$(readlink -f "$RECAF")" == "$(readlink -f "$AF")" ]] || fail "set-paths recorded alevin-fry '$RECAF'"
[[ "$(readlink -f "$RECPI")" == "$(readlink -f "$PISCEM")" ]] || fail "set-paths recorded piscem '$RECPI'"
cp -f "$INFO" "$DD/"
say "set-paths OK: ALEVIN_FRY_HOME=$AF_HOME ; alevin-fry -> $AF ; piscem -> $PISCEM"

# ---------------- cyto (built above) and Cell Ranger (pre-installed) ----------------
CYV=$("$CYTO_BIN" --version 2>&1 | head -1); CRV=$("$CELLRANGER_BIN" --version 2>&1 | head -1)
[[ "$CYV" == "cyto $CYTO_VERSION" ]] || fail "cyto is '$CYV' ($CYTO_BIN)"
[[ "$CRV" == "cellranger $CELLRANGER_VERSION" ]] || fail "cellranger is '$CRV' ($CELLRANGER_BIN)"
{ uv --version; uv tool list; } > "$DD/uv_tools.txt" 2>&1 || fail "uv (needed by cyto's cell-filter / pycyto) not usable"

{
  printf 'tool\tcommit\tversion\tbinary\tsha256\n'
  for t in simpleaf alevin-fry piscem; do
    b="$SRC/$t/target/release/$t"
    printf '%s\t%s\t%s\t%s\t%s\n' "$t" "$(git -C "$SRC/$t" rev-parse HEAD)" "$(ALEVIN_FRY_HOME="$AF_HOME" "$b" --version 2>&1 | tail -1)" "$b" "$(sha256sum "$b" | cut -c1-64)"
  done
  printf 'cyto\t%s\t%s\t%s\t%s\n' "$(git -C "$SRC/cyto" rev-parse HEAD)" "$CYV" "$CYTO_BIN" "$(sha256sum "$CYTO_BIN" | cut -c1-64)"
  printf 'cellranger\t-\t%s\t%s\t%s\n' "$CRV" "$CELLRANGER_BIN" "$(sha256sum "$CELLRANGER_BIN" | cut -c1-64)"
  printf 'piscem-rs\t-\t%s\t(Cargo.lock)\t-\nlibradicl\t-\t%s\t(Cargo.lock)\t-\n' "$PRS" "$LRC"
  printf 'rustc\t-\t%s\tRUSTFLAGS=%s\t-\n' "$(rustc --version)" "$RUSTFLAGS_BUILD"
} > "$DD/toolchain.tsv"
column -t -s$'\t' "$DD/toolchain.tsv"
write_done "toolchain.tsv	$(sha256sum "$DD/toolchain.tsv" | cut -c1-64)" \
  "simpleaf_info.json	$(sha256sum "$DD/simpleaf_info.json" | cut -c1-64)" \
  "uv_tools.txt	$(sha256sum "$DD/uv_tools.txt" | cut -c1-64)"
