#!/usr/bin/env python3
"""Reconstruct the pinned mouse TRD supplement from user-acquired sources.

No data is downloaded and no terms are accepted by this script. Outputs are local.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

EXPECTED = {
    'trdv': '41956333b48e440e2adeaf7b9546a2fe1f60ba7853b2027fb4e1b5067729c49c',
    'trdj': '6276c6ea88f4f248fa5c37886a6194a5d6b1b47f1f0927828b4b497e8401ca34',
    'base': '1d7fcac37ad6d955b3a3f693526e36623186464aef2c528cf8e0322dc89706cd',
}
OUTPUT_HASH = '482280a9dfacda0d45acf5243b68e088e7649819a05bc24b90a21f776536f649'
CODONS = {''.join(c): a for c, a in zip(itertools.product('TCAG', repeat=3),
    'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')}
AA = set('ARNDCQEGHILKMFPSTWYV.')


def checked(path, expected, label):
    try:
        data = path.read_bytes()
    except OSError as e:
        raise SystemExit(f'{label}: cannot read {path}: {e}') from e
    if hashlib.sha256(data).hexdigest() != expected:
        raise SystemExit(f'{label}: SHA-256 mismatch; the exact archived input is required')
    return data


def fasta(data):
    records = []
    for line in data.decode('utf-8').splitlines():
        if line.startswith('>'):
            records.append([line[1:].split('|'), ''])
        elif line.strip():
            if not records:
                raise SystemExit('invalid FASTA: sequence before header')
            records[-1][1] += line.strip().upper()
    if not records or len({r[0][1] for r in records}) != len(records):
        raise SystemExit('invalid FASTA: missing or duplicate records')
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ('trdv', 'trdj', 'base', 'output'):
        parser.add_argument('--' + option, required=True, type=Path)
    args = parser.parse_args()
    data = {name: checked(getattr(args, name), expected, name) for name, expected in EXPECTED.items()}
    if args.output.exists():
        raise SystemExit('output already exists; choose a new path')
    metadata = json.loads(Path(__file__).with_name('mouse_trd_provenance.json').read_text())
    accepted, excluded = [], []
    for fields, nt in fasta(data['trdv']):
        if not fields[2].startswith('Mus musculus') or fields[4] != 'V-REGION' or fields[7] != '1':
            raise SystemExit('unexpected TRDV record')
        protein = ''.join('.' if nt[i:i+3] == '...' else CODONS.get(nt[i:i+3], '?') for i in range(0, len(nt)-2, 3))
        regions = [protein[a-1:b] for a, b in [(28,39), (57,66), (82,88)]]
        reason = None
        if fields[3] not in ('F', '(F)', '[F]'):
            reason = 'nonfunctional_or_unqualified_functionality'
        elif len(protein) < 88 or any(not s.replace('.', '') or set(s) - AA for s in regions):
            reason = 'missing_or_ambiguous_scored_region'
        gene = fields[1]
        kind = 'dedicated' if gene.startswith('TRDV') else 'shared'
        if kind == 'shared' and '/DV' not in gene:
            raise SystemExit('unexpected TRDV allele namespace')
        if reason:
            excluded.append({'id': gene, 'reason': reason})
        else:
            accepted.append({'id': gene, 'kind': kind, 'regions': regions})
    j_alleles = sorted(fields[1] for fields, _ in fasta(data['trdj']))
    if len(accepted) != 51 or len(excluded) != 2 or len(j_alleles) != 3:
        raise SystemExit('unexpected admission counts')
    contract = {**metadata, 'entries': accepted, 'j_alleles': j_alleles, 'excluded': excluded}
    output = (json.dumps(contract, indent=2) + '\n').encode()
    if hashlib.sha256(output).hexdigest() != OUTPUT_HASH:
        raise SystemExit('constructed supplement has unexpected identity; no output written')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('xb') as out:
        out.write(output)
    print('Prepared pinned mouse TRD supplement; SHA-256:', OUTPUT_HASH)


if __name__ == '__main__':
    main()
