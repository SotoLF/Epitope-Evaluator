# Epitope-Evaluator

An interactive web application to study predicted T-cell epitopes.

Give it a T-cell epitope predictor's output plus the FASTA it was run on, and it
turns a table of scores into six interactive analyses. A rewrite of
[Epitope-Evaluator](https://github.com/SotoLF/Epitope-Evaluator) — same six
tools and the same scientific purpose, rebuilt for speed, correctness and
robustness. See [CHANGELOG.md](CHANGELOG.md) for what changed and why,
including the results that differ from v1.

## Quick start

```r
install.packages(c("shiny", "bslib", "plotly", "DT",
                   "data.table", "matrixStats", "stringi"))
shiny::runApp(".")
```

Then open the **Run example** tab: a prediction for the 17-protein SARS-CoV-2
proteome loads on arrival and every tool is immediately usable, with no upload.
The dropdown also carries the datasets behind the 2022 paper, including Spike
from six SARS-CoV-2 variants for the Conservation tool.

## The six tools

| Tool | Question |
|---|---|
| **Distribution** | How do epitopes spread across the score range, per allele or across allele combinations? |
| **Intersection** | Which epitopes are shared between allele combinations, and which are private to one? |
| **Density** | Which proteins are unusually epitope-rich for their length? |
| **Viewer** | Where along a protein do the epitopes sit? |
| **Promiscuity** | Which epitopes bind the most alleles at once? |
| **Conservation** | Which epitopes survive across proteins, strains or variants? |

## Input

Two files:

1. **A prediction file** from NetMHC, NetMHCpan, NetMHCIIpan, MHCflurry or IEDB
   consensus — exactly as downloaded. Do not re-save it through a spreadsheet:
   NetMHC's `.xls` files are tab-separated text and Excel rewrites the two-row
   header that names the alleles.
2. **The FASTA** submitted to the predictor. It supplies the full protein
   identifiers (predictors truncate them) and the true protein lengths.

The format is detected from the file header. MHC class is inferred from the
allele names and sets the default cutoffs: 2 %rank for class I, 10 for class II.

Any other table works too — pick **Other** and supply peptide, start position,
protein ID, protein length, then one column per allele. See the in-app
Documentation tab for each format's exact layout.

## Structure

```
app.R                      UI shell, theme, and the tabset reused by Analyse and Run example
global.R                   options, dependencies, source order
utils/
  constants.R              palette, default cutoffs, render limits
  parsing_functions.R      FASTA + five predictor formats -> one canonical dataset
  core_functions.R         the analysis engine (pure R, no Shiny)
  plot_functions.R         native plotly figures
  ui_helpers.R             error handling, cards, download handlers
modules/                   one Shiny module per tool
ui/                        About / Documentation / Tutorial
tests/
  test_core.R              engine tests -- runs without Shiny
  test_app.R               UI + reactive tests, drives every module server
  benchmark.R              scale check with a synthetic dataset
tools/screenshots.py       re-capture the documentation screenshots
tools/make_icons.R         regenerate the logo and favicon set
deploy.R                   publish to shinyapps.io (staging / production)
data/                      bundled examples (five formats + the paper's datasets)
www/styles.css             app styling
www/favicon.svg            the mark (also used as the navbar logo)
```

`utils/parsing_functions.R` and `utils/core_functions.R` have no Shiny
dependency, so the engine can be used directly:

```r
library(data.table); library(stringi); library(matrixStats)
source("utils/constants.R")
source("utils/parsing_functions.R")
source("utils/core_functions.R")

ds <- ee_parse("data/example_NetMHCPAN.xls", "data/example.fasta",
               predictor = "auto", score_type = "Rank")
ds
#> <ee_dataset> NetMHCpan / Rank
#>   14,303 rows, 9,910 unique peptides, 17 proteins, 12 alleles (MHC class I)

ee_density_table(ds, ds$alleles, mode = "Union", cutoff = 2)
ee_promiscuity(ds, ds$alleles, strong = 0.5, weak = 2, min_alleles = 8)$table
```

The canonical object is a `data.table` of peptides plus a numeric
peptide × allele score matrix:

```r
str(ds[c("peptides", "alleles", "mhc_class", "score_type")], max.level = 1)
```

## Tests

```
Rscript tests/test_core.R      # 155 assertions, needs only data.table/stringi/matrixStats
Rscript tests/test_app.R       # 226 assertions, needs the full dependency set
Rscript tests/benchmark.R      # synthetic 10^6-peptide scale check (~1 min)
```

`tests/benchmark.R [n_peptides] [n_alleles]` takes its size from the command
line, so `Rscript tests/benchmark.R 200000 8` is a quick smoke run.



## Citation

Soto, L. F., Requena, D., & Fuxman Bass, J. I. (2022). Epitope-Evaluator: An
interactive web application to study predicted T-cell epitopes. *PLoS ONE*,
17(8), e0273577. [PMC9417011](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9417011/)

## Contact

Juan Fuxman Bass — fuxman@bu.edu
Luis F. Soto — lufesu98@gmail.com

MIT licensed. See [LICENSE](LICENSE).
