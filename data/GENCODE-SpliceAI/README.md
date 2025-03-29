# data generation

We will require bedtools for grabing sequences to create the dataset that we can refer [here](https://bedtools.readthedocs.io/en/latest/content/installation.html) for installing it.

Firstly, Update the constants.py file:
- ref_genome: path of the genome.fa file (hg19/GRCh37) or (hg38.fa)
- splice_table: path for reference splicing sequences (canonical_dataset.txt for hg19 and hg38V46_splice_table.txt for hg38)
- sequence: for sequence name
- version: is used for identifying the genome version and naming files

__Download fasta file and prepare datasets__ Finally, use the following commands for data preprocessing:

- Modify ```data/constants.py``` to point to the correct version of the genome.

```sh
cd data/
/bin/bash get_transcriptome.sh

python create_datafile.py train all
python create_datafile.py test 0

python create_dataset.py train all
python create_dataset.py test 0
```