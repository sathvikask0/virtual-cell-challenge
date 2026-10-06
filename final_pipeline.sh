#!/bin/zsh
# Final round (D/E/F controls drop Oct 22; deadline Nov 5 23:59 UTC). Run from the repo root, step by step.
#   ./final_pipeline.sh DIR "D,E,F" NAME [recipe options...]
# DIR = unzipped final controls (context_X.h5ad, gene_names.csv, pert_counts.csv).
# Download first (needs the user's OK): .venv/bin/vcc datasets list ; .venv/bin/vcc datasets download <id>
# Dry run on the validation set: ./final_pipeline.sh data/vcc/controls "A,B,C" dry_rmpc5 <rmpc5 options>
set -e
DIR=$1; CTX=$2; NAME=$3; shift 3
export VCC_CTRL_DIR=$(cd $DIR && pwd) VCC_CONTEXTS=$CTX ATLAS_CACHE=$PWD/data/atlas_shift_$NAME
PY=$PWD/.venv/bin/python
mkdir -p $ATLAS_CACHE
cd src
NAMES=(); for c in ${(s:,:)CTX}; do NAMES+=ctx_$c; done
echo "== 1. context control profiles"; for n in $NAMES; do [ -f ../data/lines/$n.npz ] || $PY lines.py $n; done
echo "== 2. identity (DepMap)"; $PY line_identity.py 2000 $NAMES
echo "== 3. source caches for this target panel"; $PY atlas_shift.py build
$PY -c "import atlas_cd4 as C; C.build()"
echo "== 4. build + validate $NAME"; $PY agree_alloc.py submission $NAME "$@"
$PY validate_atlas_submission.py $NAME && rm -f ../data/submissions/$NAME.h5ad
echo "== done: data/submissions/$NAME.vcc (NOT uploaded; upload needs the user's explicit OK)"
