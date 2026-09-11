#!/usr/bin/env bash
#
# Publish this working copy to SotoLF/Epitope-Evaluator as a branch, so the
# change can be reviewed as a pull request before main moves.
#
#   tools/push_to_github.sh                 # dry run against the default branch
#   tools/push_to_github.sh --push          # push to the v2 branch (reviewable)
#   tools/push_to_github.sh --push --main   # commit straight onto main
#   tools/push_to_github.sh --push --branch v2-rewrite
#
# The default target is a branch, so the change can be reviewed as a pull
# request before main moves -- this repository is cited in the 2022 paper.
# --main skips that. It is still an ordinary commit on top of main, not a
# force-push, so v1 stays in history and `git revert` undoes it.
set -euo pipefail

REMOTE="git@github.com:SotoLF/Epitope-Evaluator.git"
BRANCH="v2"
DO_PUSH=0
SRC="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${TMPDIR:-/tmp}/ee-push-$$"

while [ $# -gt 0 ]; do
  case "$1" in
    --push)   DO_PUSH=1 ;;
    --main)   BRANCH="main" ;;
    --branch) BRANCH="$2"; shift ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
  shift
done

# --- 1. Refuse to run if a credential could reach the commit ---------------
# .Renviron holds the shinyapps.io token. It is gitignored, but "gitignored"
# is worth verifying rather than trusting, because a leaked deploy token for a
# shared lab account is not something a later commit can undo.
if [ -f "$SRC/.Renviron" ]; then
  TOK=$(sed -n 's/^SHINYAPPS_TOKEN=//p'  "$SRC/.Renviron" | head -1)
  SEC=$(sed -n 's/^SHINYAPPS_SECRET=//p' "$SRC/.Renviron" | head -1)
else
  TOK=""; SEC=""
fi

if [ "$BRANCH" = "main" ]; then
  echo "== Target: main (published branch -- an ordinary commit, revertible)"
else
  echo "== Target: branch $BRANCH (main untouched)"
fi
echo "== Preparing $BRANCH from $SRC"

# --- 2. Fresh full clone: the existing checkout is a --depth 1 shallow copy,
#        which cannot push a reviewable history.
rm -rf "$WORK"
git clone --quiet "$REMOTE" "$WORK"
cd "$WORK"
git checkout --quiet -B "$BRANCH" origin/main
echo "   branched from origin/main @ $(git rev-parse --short HEAD)"

# --- 3. Sync the tree.
#        --delete so files v1 had and v2 does not are removed, with two
#        exceptions that keep published paths resolving on a cited repository:
#        Biological_Applications_Data/ (the paper's datasets, which the app now
#        reads directly) and Images/ (13 figures the v1 README embedded by
#        absolute GitHub URL, so external deep links to them still work).
rsync -a --delete \
  --exclude '.git/' \
  --exclude '.Renviron' \
  --exclude '.Rhistory' \
  --exclude '.RData' \
  --exclude 'rsconnect/' \
  --exclude 'Biological_Applications_Data/' \
  --exclude 'Images/' \
  --exclude '*.log' \
  "$SRC"/ ./

# --- 4. Verify no credential made it into the tree, then stage.
git add -A
if [ -n "$TOK" ] && git grep -qF -- "$TOK" -- $(git diff --cached --name-only) 2>/dev/null; then
  echo "!! ABORTING: the shinyapps.io token appears in a staged file." >&2; exit 1
fi
if [ -n "$SEC" ] && git grep -qF -- "$SEC" -- $(git diff --cached --name-only) 2>/dev/null; then
  echo "!! ABORTING: the shinyapps.io secret appears in a staged file." >&2; exit 1
fi
if git ls-files --cached | grep -qx '.Renviron'; then
  echo "!! ABORTING: .Renviron is staged." >&2; exit 1
fi
echo "   credential check: clean"

# --- 5. Report.
echo
echo "== Change against main"
git diff --cached --stat origin/main | tail -20
echo
ADD=$(git diff --cached --diff-filter=A --name-only origin/main | wc -l)
DEL=$(git diff --cached --diff-filter=D --name-only origin/main | wc -l)
MOD=$(git diff --cached --diff-filter=M --name-only origin/main | wc -l)
printf "   %s added, %s removed, %s modified\n" "$ADD" "$DEL" "$MOD"
echo "   largest tracked files:"
git ls-files --cached -z | xargs -0 du -h 2>/dev/null | sort -rh | head -4 | sed 's/^/     /'

if [ "$DO_PUSH" -eq 0 ]; then
  echo
  echo "Dry run. Nothing was pushed and main was not touched."
  echo "Re-run with --push to create the branch on GitHub."
  echo "Staged tree left at: $WORK"
  exit 0
fi

# --- 6. Commit and push the branch.
git commit --quiet -F - <<'MSG'
Rewrite: faster engine, corrected parsing, native plotly figures

Same six tools and the same scientific purpose; the parsing, analysis
engine, figures and error handling are new. See CHANGELOG.md for the
full account, including the results that differ from v1.

Corrections that change numbers:
- Protein length is the true amino-acid length from the FASTA. v1 used
  nchar(sequence) - 8 for every predictor, which was the denominator of
  every epitope density.
- Protein IDs are matched to the FASTA by name with a length cross-check.
  v1 matched by position, which silently misassigned every ID whenever a
  protein was too short to yield a peptide.
- NetMHCpan/NetMHCIIpan "Score" reads the right column. v1's column
  stride landed on a text field, so those analyses were all-NA.
- Promiscuity counts alleles with <=, matching its own documentation and
  heatmap; v1's table used a strict <.
- Missing scores (MHCflurry, IEDB) are excluded from the AND/OR
  combination instead of poisoning the whole peptide.

Defects that stopped v1 running on current R:
- UpSet plots called dplyr::filter_(), defunct since dplyr 1.1.
- Input validation evaluated if() on a length-n logical, an error since
  R 4.2, so it failed on every multi-FASTA.
- Class-aware default cutoffs were applied from the wrong module session
  and never reached the tools.

Performance, measured on a synthetic 731 MB NetMHCpan table
(10^6 peptides x 20 alleles): parsing 4.0 s, all six tools 8.0 s.
Viewer lane packing went from 119 s to 0.004 s for 800 epitopes by
replacing the pairwise overlap search with a first-fit interval sweep.
UpSet no longer enumerates 2^k combinations.

Dependencies reduced from 22 to 7 (shiny, bslib, plotly, DT,
data.table, matrixStats, stringi). sf/ggVennDiagram are gone, so
deployment no longer needs GDAL/GEOS/PROJ on the host; shinyBS,
reshape and reshape2 are archived on CRAN and no longer used.

modules/example/ is deleted: 2,149 lines duplicating the seven modules
with an e_ prefix. Shiny namespaces already allow one module set to
serve both the Analyse and Run example tabs.

Adds tests/test_core.R (155 assertions, no Shiny required),
tests/test_app.R (226 assertions, drives every module server),
tests/benchmark.R, and tools/ for icons, screenshots and deployment.
MSG

PREV=$(git rev-parse --short origin/main)
git push --quiet -u origin "$BRANCH"
NOW=$(git rev-parse --short HEAD)
echo
echo "== Pushed $BRANCH  ($PREV -> $NOW)"
if [ "$BRANCH" = "main" ]; then
  echo "   https://github.com/SotoLF/Epitope-Evaluator"
  echo
  echo "   v1 remains at $PREV. To undo:"
  echo "     git clone git@github.com:SotoLF/Epitope-Evaluator.git && cd Epitope-Evaluator"
  echo "     git revert --no-commit $NOW && git commit -m 'Revert to v1' && git push"
else
  echo "   Open a pull request:"
  echo "   https://github.com/SotoLF/Epitope-Evaluator/compare/main...$BRANCH?expand=1"
fi
echo "   Working clone left at: $WORK"
