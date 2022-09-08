import sys
import hashlib
from pathlib import Path

import defopt

from Bio import SeqIO


import pandas as pd


def hash_bed(bed_path: Path):
    bed_fields = ['chrom', 'chromStart', 'chromEnd', 'name', 'poolName', 'orientation']
    fields_to_hash = ['chrom', 'chromStart', 'chromEnd', 'orientation']

    df = pd.read_csv(
        bed_path,
        sep='\t',
        names=bed_fields,
        dtype=dict(
            chrom=str,
            chromStart=int,
            chromEnd=int,
            name=str,
            poolName=str,
            orientation=str,
        )
    )

    # Normalisation
    df_norm = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)  # Strip trailing and leading whitespace
    df_norm = df_norm[[*fields_to_hash]].sort_values('chromStart')  # Sort by start pos
    df_norm_text = df_norm.to_csv(sep='\t', header=False, index=False)  

    hex_digest = hashlib.md5(df_norm_text.encode()).hexdigest()

    print('BED checksum:', file=sys.stderr)
    print(hex_digest)
    return hex_digest
    # print('Hash function input:', df_norm_text.strip(), sep='\n', file=sys.stderr)


def hash_ref(ref_path: Path):
    record = SeqIO.read(ref_path, "fasta")
    hex_digest = hashlib.md5(str(record.seq).upper().encode()).hexdigest()
    print('Reference checksum:', file=sys.stderr)
    print(hex_digest)
    return hex_digest


def build():
    pass


def main():
    defopt.run(
        {
            "hash-bed": hash_bed,
            "hash-ref": hash_ref,
            "build": build,
        },
        no_negated_flags=True,
        strict_kwonly=False,
        short={},
    )


if __name__ == "__main__":
    main()