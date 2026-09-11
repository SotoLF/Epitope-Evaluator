# Example data

Predicted MHC binding for the SARS-CoV-2 proteome (17 proteins, UniProt), used
by the **Run example** tab. `EE_EXAMPLES` in `modules/data_input.R` maps each
label in the dropdown to one of these files.

| File | Predictor | Class | Alleles |
|---|---|---|---|
| `example_NetMHCPAN.xls`     | NetMHCpan 4.1   | I  | 12 |
| `example_NetMHCIIPAN.xls`   | NetMHCIIpan 4.0 | II | 8  |
| `example_NetMHC.xls`        | NetMHC 4.0      | I  | 12 |
| `example_MHCFlurry.txt`     | MHCflurry 2.0   | I  | 11 |
| `example_IEDB_consensus.txt`| IEDB consensus  | I  | 7  |

`example.fasta` is the multi-FASTA submitted to every predictor, and is
required alongside any of the prediction files.

These files are also the fixtures for `tests/test_core.R` and
`tests/test_app.R`, so removing one will fail the test suite. If you need a
smaller deployment bundle, drop an entry from `EE_EXAMPLES` and delete the
matching file and its test cases together.
