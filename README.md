# Primaschema

**Currently under active development and unstable**. Generates PHA4GE-compatible tiling amplicon scheme bundles from either 6 or 7 column Primal Scheme-like bed files and metadata contained in a YAML file.

## Install (Python 3.10+)
```
pip install https://github.com/bede/primaschema
```



## Usage

```
% primaschema --help
usage: primaschema [-h] [--version] {hash-ref,hash-bed,validate,build} ...

positional arguments:
  {hash-ref,hash-bed,validate,build}
    hash-ref
    hash-bed
    validate
    build

options:
  -h, --help            show this help message and exit
  --version             show program's version number and exit


% ls tests/data/primer-schemes/eden/v1
reference.fasta	scheme.bed	info.yaml


% primaschema build tests/data/primer-schemes/eden/v1 --force
INFO: Scheme bed file has the expected number of columns (6)
INFO: Writing info.yaml with checksums
INFO: Generating primer.bed from scheme.bed and reference.fasta


% ls eden-v1
primer.bed	info.yaml

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
notes: Protocol includes addendum for ONT sequencing
primer_checksum: "primaschema:9e4c6a3b84cbd76cb3e38b893d0322b5799ecafe28d8cf7bf347ce6dcc5ee8cb"
reference_checksum: "primaschema:7d5621cd3b3e498d0c27fcca9d3d3c5168c7f3d3f9776f3005c7011bd90068ca"
```



## Test

```
$ pytest
```
