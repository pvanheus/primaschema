[![Tests](https://github.com/pha4ge/primaschema/actions/workflows/test.yml/badge.svg)](https://github.com/pha4ge/primaschema/actions/workflows/test.yml) [![PyPI version](https://badge.fury.io/py/primaschema.svg)](https://pypi.org/project/primaschema)

# Primaschema

The toolkit for validating and building tiling amplicon PCR primer scheme definitions for inclusion in the [PHA4GE primer-schemes repository, using 6 or 7 column Primal Scheme-like BED files and scheme metadata contained in a YAML file.



## Install (Python 3.10+)

```
# Latest stable release
pip install primaschema

# From main branch
git clone https://github.com/pha4ge/primaschema
pip install ./primaschema

# Development
git clone https://github.com/pha4ge/primaschema.git
cd primaschema
pip install --editable ./
pytest
```

Some Primaschema commands use components from the [primer-schemes](https://github.com/pha4ge/primer-schemes) repository. To show Primaschema where to find these, create the environment variable `PRIMER_SCHEMES_PATH` pointing to the location of the primer-schemes directory on your machine:

```
git clone https://github.com/pha4ge/primer-schemes.git
export PRIMER_SCHEMES_PATH="/path/to/primer-schemes"
```



## Usage

```
% primaschema --help
usage: primaschema [-h] [--version]
                   {hash-ref,hash-bed,validate,validate-recursive,build,build-recursive,build-manifest,diff,6to7,7to6,show-non-ref-alts}
                   ...

positional arguments:
  {hash-ref,hash-bed,validate,validate-recursive,build,build-recursive,build-manifest,diff,6to7,7to6,show-non-ref-alts}
    hash-ref            Generate reference sequence checksum
    hash-bed            Generate a bed file checksum
    validate            Validate a primer scheme bundle containing info.yml, primer.bed and reference.fasta
    validate-recursive  Recursively validate primer scheme bundles in the specified directory
    build               Build a primer scheme bundle containing info.yml, primer.bed and reference.fasta
    build-recursive     Recursively build primer scheme bundles in the specified directory
    build-manifest      Build a complete manifest of schemes contained in the specified directory
    diff                Show the symmetric difference of records in two bed files
    6to7                Convert a 6 column scheme.bed file to a 7 column primer.bed file using a reference sequence
    7to6                Convert a 7 column primer.bed file to a 6 column scheme.bed file by droppign a column
    show-non-ref-alts   Show primer records with sequences not matching the reference sequence

options:
  -h, --help            show this help message and exit
  --version             show program's version number and exit


% primaschema build test/data/primer-schemes/eden/v1
INFO: Scheme bed file has the expected number of columns (6)
INFO: Writing info.yml with checksums
INFO: Generating primer.bed from scheme.bed and reference.fasta
```
