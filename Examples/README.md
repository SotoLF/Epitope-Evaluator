## Example data

### example.fasta
The Fasta file contains the whole proteome of SARS-CoV2 in a fasta format obtained from Uniprot.

### example.xls
The output prediction file obtained using [NetMHCPan4.1](https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1) with the following parameters:
* Fasta file: example.fasta
* Peptide length: 9mer peptides
* Alleles: HLA-A01:01,HLA-A02:01,HLA-A03:01,HLA-A24:02,HLA-A26:01,HLA-B07:02,HLA-B08:01,HLA-B15:01,HLA-B27:05,HLA-B39:01,HLA-B40:01,HLA-B58:01
* Strong Binder cutoff: 0.5
* Weak Binder cutoff: 2
* Save prediction to XLS file: yes

### example_IEDB_conensus.txt
The output prediction file obtained using [IEDB Consensus](http://tools.iedb.org/mhci/) with the following parameters:
* Fasta file: example.fasta
* Prediciton Method: Consensus
* Alleles: HLA-A01:01,HLA-A02:01,HLA-A02:06,HLA-A03:01,HLA-A11:01,HLA-A23:01,HLA-A32:01

### example_MHCFlurry.txt
The output prediction file obtained using [MHCFlurry 2.0](https://https://github.com/openvax/mhcflurry) with the following command line:
* mhcflurry-predict-scan example.fasta --alleles HLA-A*01:01 HLA-A*02:01 HLA-A*03:01 HLA-A*24:02 HLA-A*26:01 HLA-B*07:02 HLA-B*08:01 HLA-B*27:05 HLA-B*39:01 HLA-B*40:01 HLA-B*15:01 --peptide-lengths 9 --no-throw --results-all --out example_MHCFlurry.txt.

### example_NetMHC.xls
The output prediction file obtained using [NetMHC4.0](https://services.healthtech.dtu.dk/service.php?NetMHC-4.0) with the following parameters:
* Fasta file: example.fasta
* Peptide length: 9mer peptides
* Alleles: HLA-A01:01,HLA-A02:01,HLA-A03:01,HLA-A24:02,HLA-A26:01,HLA-B07:02,HLA-B08:01,HLA-B15:01,HLA-B27:05,HLA-B39:01,HLA-B40:01,HLA-B58:01
* Strong Binder cutoff: 0.5
* Weak Binder cutoff: 2
* Save prediction to XLS file: yes

### example_NetMHCIIPAN.xls
The output prediction file obtained using [NetMHCIIPan4.0](https://services.healthtech.dtu.dk/service.php?NetMHCIIpan-4.0) with the following parameters:
* Fasta file: example.fasta
* Peptide length: 15mer peptides
* Alleles: DRB1_0101, DRB1_0103, DRB1_0301, DRB1_0101,  DRB1_0401, DRB1_0402, DRB1_0403, DRB1_0404, DRB1_0405     
* Strong Binder cutoff: 1
* Weak Binder cutoff: 5
* Save prediction to XLS file: yes
