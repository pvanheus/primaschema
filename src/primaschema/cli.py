import sys
import hashlib
from pathlib import Path

import defopt


import pandas as pd


def build(bed_path: Path):
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

    hexdigest = hashlib.md5(df_norm_text.encode()).hexdigest()

    print('Primer scheme hash:', file=sys.stderr)
    print(hexdigest[:16])  # Use first 64 bits only
    print('Hash function input:', df_norm_text.strip(), sep='\n', file=sys.stderr)




def main():
    defopt.run(
        {
            "build": build,
        },
        no_negated_flags=True,
        strict_kwonly=False,
        short={},
    )


if __name__ == "__main__":
    main()