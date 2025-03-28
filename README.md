# ml-models-splicing

A repository for a collection of ML models that identify and/or quantify slice sites. All the models are stored as submodules and have tehyr own ```environment.yml``` file for setting up a conda environment using

```bash
conda env create --file environment.yml
```

The folder structure is at follows. Each model has its own README.md where the installation and/or details for loading models can be found.

```bash
.
├── hyena-dna
│   ├── assets
│   ├── configs
│   ├── csrc
│   ├── evals
│   ├── flash-attention
│   ├── src
│   ├── Dockerfile
│   ├── environment.yml
│   ├── huggingface.py
│   ├── LICENSE
│   ├── README.md
│   ├── standalone_hyenadna.py
│   └── train.py
├── nucleotide-transformer
│   ├── examples
│   ├── imgs
│   ├── nucleotide_transformer
│   ├── environment.yml
│   ├── LICENSE.md
│   ├── README.md
│   └── setup.py
├── Pangolin
│   ├── examples
│   ├── pangolin
│   ├── scripts
│   ├── brca_pangolin.vcf
│   ├── environment.yml
│   ├── LICENSE
│   ├── PangolinColab.ipynb
│   ├── README.md
│   └── setup.py
├── RiNALMo
│   ├── ft_schedules
│   ├── imgs
│   ├── rinalmo
│   ├── weights
│   ├── environment.yml
│   ├── LICENSE
│   ├── pyproject.toml
│   ├── README.md
│   ├── train_expression_level.py
│   ├── train_ncrna_classification.py
│   ├── train_ribosome_loading.py
│   ├── train_sec_struct_prediction.py
│   ├── train_splice_site_prediction.py
│   └── train_translation_efficiency.py
├── SpliceTransformer
│   ├── data
│   ├── model
│   ├── custom_usage.py
│   ├── environment.yml
│   ├── LICENSE
│   ├── motif_analyze.py
│   ├── old_annotator_pytorch_2.py
│   ├── old_annotator_pytorch.py
│   ├── README.md
│   ├── sptransformer.py
│   └── tasks_annotate_mutations.py
├── LICENSE
└── README.md
```