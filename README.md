![Tests](https://github.com/pha4ge/primaschema/actions/workflows/test.yml/badge.svg) [![PyPI version](https://badge.fury.io/py/primaschema.svg)](https://pypi.org/project/primaschema)
# Primaschema

**Currently under active development and unstable**. Generates PHA4GE-compatible tiling amplicon scheme bundles from either 6 or 7 column Primal Scheme-like bed files and metadata contained in a YAML file.

## Install (Python 3.10+)
```
pip install primaschema  # Latest release
pip install https://github.com/pha4ge/primascheman  # Main branch
```



## Usage

```
% primaschema --help
usage: primaschema [-h] [--version] {hash-ref,hash-bed,validate,build,6to7} ...

positional arguments:
  {hash-ref,hash-bed,validate,build,6to7}
    hash-ref            Generate reference sequence checksum
    hash-bed            Generate a bed file checksum
    validate            Validate a primer scheme bundle containing info.yaml, primer.bed and reference.fasta
    build               Build a primer scheme bundle containing info.yaml, primer.bed and reference.fasta
    6to7                Convert a 6 column scheme.bed file to a 7 column primer.bed file using a reference sequence

options:
  -h, --help            show this help message and exit
  --version             show program's version number and exit


% primaschema build tests/data/primer-schemes/eden/v1
INFO: Scheme bed file has the expected number of columns (6)
INFO: Writing info.yaml with checksums
INFO: Generating primer.bed from scheme.bed and reference.fasta
```



## `info.yaml` example

```yaml
schema_version: 1-0-0
name: eden-v1
organism: SARS-CoV-2
organism_aliases:
  - nCoV-2019
aliases:
  - sydney
developers:
  - name: John-Sebastian Eden
    url: https://orcid.org/0000-0003-1374-3551
  - name: Eby Sim
    url: https://orcid.org/0000-0002-3716-7344
vendors: []
amplicon_size: 2500
repository_url: https://github.com/pha4ge/primer-schemes/tree/main/sars-cov-2/eden/v1
citations:
  - https://www.protocols.io/view/sars-cov-2-genome-sequencing-using-long-pooled-amp-kxygxeob4v8j/v1
notes:
  - Protocol includes addendum for ONT sequencing
primer_checksum: "primaschema:9e4c6a3b84cbd76cb3e38b893d0322b5799ecafe28d8cf7bf347ce6dcc5ee8cb"
reference_checksum: "primaschema:7d5621cd3b3e498d0c27fcca9d3d3c5168c7f3d3f9776f3005c7011bd90068ca"
```



## Development

```shell
# Inside a clean Python 3.10+ environment
git clone https://github.com/pha4ge/primaschema.git
cd primaschema
pip install --editable ./
pytest
```
