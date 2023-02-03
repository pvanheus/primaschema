#!/usr/bin/env python3

import argparse
import os
import sys

import yaml

def is_known(schemes, scheme, organism):
    for scheme_dict in schemes:
        scheme_organism = scheme_dict['organism'].lower()
        scheme_name = scheme_dict['name'].lower()
        if scheme == scheme_name and organism == scheme_organism:
            return True
    else:
        return False
        
if __name__ == '__main__':
    organism_dict = {
        'sars-cov-2': 'SARS-CoV-2',
        'mpxv': 'MPXV'
    }

    parser = argparse.ArgumentParser(description='Script to read the primers-schemes repo and YAML-ize contents for the index.yml')
    parser.add_argument('--index_filename', default='index.yml')
    parser.add_argument('--primer_schemes_dir', default='.')
    parser.add_argument('output_file', type=argparse.FileType('w'), nargs='?', default=sys.stdout)
    args = parser.parse_args()

    base_url = 'https://raw.githubusercontent.com/PHA4GE/primer-schemes/main'
    scheme_info = yaml.safe_load(open(args.index_filename))
    schemes = scheme_info['schemes']
    new_schemes = {}
    for root, dirs, files in os.walk(args.primer_schemes_dir):
        for filename in files:
            if filename.endswith('info.yaml'):
                (_, organism, scheme, version) = root.split('/')
                if not is_known(schemes, scheme, organism):
                    key = scheme+organism
                    if key in new_schemes:
                        record = new_schemes[key]
                        url = f'{base_url}/{organism}/{scheme}/{version}'
                        new_version = dict(version=version, url=url)
                        record['versions'].append(new_version)
                    else:
                        url = f'{base_url}/{organism}/{scheme}/{version}'
                        record = dict(
                            name=scheme,
                            organism=organism_dict[organism],
                            versions=[dict(version=version, url=url)])
                        new_schemes[key] = record

    schemes_added = []
    for key in new_schemes:
        schemes_added.append(new_schemes[key])
    if not schemes_added:
        print("no new schemes found", file=sys.stderr)
    else:
        yaml.dump(schemes_added, args.output_file, sort_keys=False)