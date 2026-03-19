#!/usr/bin/env python3
"""
Convert ProteinMPNN output PDBs to Boltz2 YAML input files.

This script reads antibody-antigen PDB files (HLT format) and generates
YAML files compatible with `boltz predict`. Each YAML specifies the
protein chains and their sequences.

Usage:
    python scripts/prepare_boltz2_input.py -i mpnn_designs/ -o boltz2_input/

Runs in the RFantibody conda environment (only needs standard deps).
"""

import argparse
import sys
from pathlib import Path

import yaml


# ---------------------------------------------------------------------------
# PDB parsing (lightweight — avoids BioPython dependency)
# ---------------------------------------------------------------------------

AA_3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'MSE': 'M',  # Selenomethionine → Met
}


def extract_chains_from_pdb(pdb_path: str) -> dict[str, str]:
    """Extract chain sequences from a PDB file.

    Reads CA atoms to get one residue per position, preserving chain order.

    Returns:
        Dict of {chain_id: sequence_string}
    """
    chains: dict[str, list[str]] = {}
    seen: set[tuple[str, int]] = set()

    with open(pdb_path) as f:
        for line in f:
            if not (line.startswith('ATOM') or line.startswith('HETATM')):
                continue
            atom_name = line[12:16].strip()
            if atom_name != 'CA':
                continue

            chain_id = line[21].strip()
            resnum = int(line[22:26].strip())
            resname = line[17:20].strip()

            key = (chain_id, resnum)
            if key in seen:
                continue
            seen.add(key)

            aa = AA_3TO1.get(resname, 'X')
            chains.setdefault(chain_id, []).append(aa)

    return {ch: ''.join(residues) for ch, residues in chains.items()}


# ---------------------------------------------------------------------------
# YAML generation
# ---------------------------------------------------------------------------

def generate_boltz_yaml(
    complex_name: str,
    chain_sequences: dict[str, str],
    output_dir: Path,
    remap_chains: bool = False,
    target_chains: set[str] | None = None,
) -> Path:
    """Generate a Boltz2-compatible YAML file for a protein complex.

    When remap_chains is True, chains are reordered and renamed:
      - Target chains come first (A, B, C, ...)
      - Binder chains follow (next letters)
    Target chains are identified by target_chains set (default: {'T'}).

    Args:
        complex_name: Name for the output file (without extension).
        chain_sequences: Dict of {chain_id: sequence}.
        output_dir: Directory to write the YAML file.
        remap_chains: If True, reorder target-first and rename to A/B/C/...
        target_chains: Set of chain IDs that are target (default: {'T'}).

    Returns:
        Path to the generated YAML file.
    """
    if target_chains is None:
        target_chains = {'T'}

    if remap_chains:
        # Reorder: target chains first, then binder chains
        target_seqs = [(ch, seq) for ch, seq in chain_sequences.items() if ch in target_chains]
        binder_seqs = [(ch, seq) for ch, seq in chain_sequences.items() if ch not in target_chains]
        ordered = target_seqs + binder_seqs

        # Assign new chain IDs: A, B, C, D, ...
        new_ids = [chr(ord('A') + i) for i in range(len(ordered))]
        sequences = []
        for (old_ch, seq), new_id in zip(ordered, new_ids):
            cleaned = ''.join(c if c in 'ACDEFGHIKLMNPQRSTVWY' else 'X' for c in seq)
            sequences.append({
                'protein': {
                    'id': new_id,
                    'sequence': cleaned,
                }
            })

        n_target = len(target_seqs)
        n_binder = len(binder_seqs)
        target_ids = new_ids[:n_target]
        binder_ids = new_ids[n_target:]
    else:
        sequences = []
        for chain_id, seq in chain_sequences.items():
            cleaned = ''.join(c if c in 'ACDEFGHIKLMNPQRSTVWY' else 'X' for c in seq)
            sequences.append({
                'protein': {
                    'id': chain_id,
                    'sequence': cleaned,
                }
            })

    yaml_data = {
        'version': 1,
        'sequences': sequences,
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    yaml_path = output_dir / f'{complex_name}.yaml'
    with open(yaml_path, 'w') as f:
        yaml.dump(yaml_data, f, default_flow_style=False, sort_keys=False)

    return yaml_path


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Convert ProteinMPNN output PDBs to Boltz2 YAML inputs'
    )
    parser.add_argument(
        '-i', '--input-dir', type=str, required=True,
        help='Directory containing PDB files from ProteinMPNN'
    )
    parser.add_argument(
        '-o', '--output-dir', type=str, required=True,
        help='Directory to write Boltz2 YAML files'
    )
    parser.add_argument(
        '--remap-chains', action='store_true', default=False,
        help='Reorder and rename chains: target first (A,B,...), then binder'
    )
    parser.add_argument(
        '--target-chains', type=str, default='T',
        help='Comma-separated chain IDs that are target (default: T)'
    )
    parser.add_argument(
        '--min-chain-length', type=int, default=5,
        help='Skip chains shorter than this (default: 5)'
    )
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)

    pdb_files = sorted(input_dir.glob('*.pdb'))
    if not pdb_files:
        print(f'Error: No PDB files found in {input_dir}')
        sys.exit(1)

    target_chain_set = set(args.target_chains.split(','))

    print(f'Preparing Boltz2 inputs for {len(pdb_files)} PDB files...')
    if args.remap_chains:
        print(f'  Chain ordering: target ({",".join(target_chain_set)}) → A,B,... then binder → next letters')
    generated = 0
    skipped = 0

    for pdb_path in pdb_files:
        try:
            chains = extract_chains_from_pdb(str(pdb_path))

            # Filter short chains
            chains = {
                ch: seq for ch, seq in chains.items()
                if len(seq) >= args.min_chain_length
            }

            if len(chains) < 2:
                print(f'  Skipping {pdb_path.name}: fewer than 2 chains after filtering')
                skipped += 1
                continue

            name = pdb_path.stem
            yaml_path = generate_boltz_yaml(
                name, chains, output_dir,
                remap_chains=args.remap_chains,
                target_chains=target_chain_set,
            )
            generated += 1

        except Exception as e:
            print(f'  Error processing {pdb_path.name}: {e}')
            skipped += 1

    print(f'\nDone: {generated} YAML files generated, {skipped} skipped')
    print(f'Output: {output_dir}/')


if __name__ == '__main__':
    main()
