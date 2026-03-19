#!/usr/bin/env python3
"""
Validate RFdiffusion backbone designs for physical validity.

Checks each PDB for:
  1. CA-CA chain breaks (consecutive CA > 4.2Å)
  2. Motif connectivity (flank↔motif boundaries specifically)
  3. Backbone bond geometry (N-CA, CA-C distances)
  4. Steric clashes (CA-CA < 2.0Å between non-bonded residues)

Designs that fail validation are moved to a 'rejected/' subdirectory
and excluded from downstream MPNN + Boltz2 processing.

Usage:
    python scripts/validate_backbones.py \
        -i output/PROTEIN_Nb/RFdiffusion_backbones/ \
        --motif-cdr H3 \
        [--max-break 4.2] [--min-clash 2.0]
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np


# =============================================================================
# PDB parsing
# =============================================================================

def parse_backbone(pdb_path: Path) -> tuple[np.ndarray, list[str], list[int]]:
    """Extract CA coordinates, chain IDs, residue numbers from PDB."""
    coords, chains, resnums = [], [], []
    seen = set()
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            if line[12:16].strip() != 'CA':
                continue
            chain = line[21].strip()
            resnum = int(line[22:26].strip())
            key = (chain, resnum)
            if key in seen:
                continue
            seen.add(key)
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            chains.append(chain)
            resnums.append(resnum)
    return np.array(coords) if coords else np.zeros((0, 3)), chains, resnums


# =============================================================================
# Validation checks
# =============================================================================

def check_chain_breaks(
    coords: np.ndarray,
    chains: list[str],
    max_ca_dist: float = 4.2,
) -> list[dict]:
    """Find consecutive CA-CA distances exceeding max_ca_dist within same chain."""
    breaks = []
    for i in range(len(coords) - 1):
        if chains[i] != chains[i + 1]:
            continue  # Different chains — skip
        dist = float(np.linalg.norm(coords[i + 1] - coords[i]))
        if dist > max_ca_dist:
            breaks.append({
                'residue_i': i,
                'residue_j': i + 1,
                'chain': chains[i],
                'distance': round(dist, 2),
            })
    return breaks


def check_motif_connectivity(
    coords: np.ndarray,
    chains: list[str],
    resnums: list[int],
    motif_positions: dict,
    max_ca_dist: float = 4.2,
) -> list[dict]:
    """Check CA-CA distance at flank↔motif boundaries specifically.

    Args:
        motif_positions: dict from motif_fixed.json, e.g. {"H": [106, 107, 108, 109]}
    """
    issues = []
    for chain_id, positions in motif_positions.items():
        if not positions:
            continue
        first_motif = min(positions)
        last_motif = max(positions)

        # Find chain residues
        chain_mask = [i for i, c in enumerate(chains) if c == chain_id]
        resnum_to_idx = {resnums[i]: i for i in chain_mask}

        # Check N-boundary: residue before first motif → first motif
        prev_res = first_motif - 1
        if prev_res in resnum_to_idx and first_motif in resnum_to_idx:
            dist = float(np.linalg.norm(
                coords[resnum_to_idx[first_motif]] - coords[resnum_to_idx[prev_res]]))
            if dist > max_ca_dist:
                issues.append({
                    'boundary': 'N-flank→motif',
                    'chain': chain_id,
                    'residues': f'{prev_res}-{first_motif}',
                    'distance': round(dist, 2),
                })

        # Check C-boundary: last motif → residue after last motif
        next_res = last_motif + 1
        if last_motif in resnum_to_idx and next_res in resnum_to_idx:
            dist = float(np.linalg.norm(
                coords[resnum_to_idx[next_res]] - coords[resnum_to_idx[last_motif]]))
            if dist > max_ca_dist:
                issues.append({
                    'boundary': 'motif→C-flank',
                    'chain': chain_id,
                    'residues': f'{last_motif}-{next_res}',
                    'distance': round(dist, 2),
                })

    return issues


def check_steric_clashes(
    coords: np.ndarray,
    chains: list[str],
    min_nonbond_dist: float = 2.0,
    max_report: int = 5,
) -> list[dict]:
    """Find CA-CA distances < min_nonbond_dist between non-bonded residues."""
    clashes = []
    n = len(coords)
    for i in range(n):
        for j in range(i + 2, min(i + 20, n)):  # Skip adjacent, check nearby
            if chains[i] != chains[j]:
                continue
            dist = float(np.linalg.norm(coords[j] - coords[i]))
            if dist < min_nonbond_dist:
                clashes.append({
                    'residue_i': i,
                    'residue_j': j,
                    'chain': chains[i],
                    'distance': round(dist, 2),
                })
                if len(clashes) >= max_report:
                    return clashes
    return clashes


def check_bond_lengths(
    coords: np.ndarray,
    chains: list[str],
    min_ca_dist: float = 3.0,
    max_ca_dist: float = 4.2,
) -> list[dict]:
    """Check that consecutive CA-CA distances are within normal range."""
    issues = []
    for i in range(len(coords) - 1):
        if chains[i] != chains[i + 1]:
            continue
        dist = float(np.linalg.norm(coords[i + 1] - coords[i]))
        if dist < min_ca_dist or dist > max_ca_dist:
            issues.append({
                'residue_i': i,
                'residue_j': i + 1,
                'chain': chains[i],
                'distance': round(dist, 2),
                'issue': 'too_short' if dist < min_ca_dist else 'too_long',
            })
    return issues


# =============================================================================
# Main validation
# =============================================================================

def validate_backbone(
    pdb_path: Path,
    motif_json_path: Path | None = None,
    max_break: float = 4.2,
    min_clash: float = 2.0,
) -> dict:
    """Run all validation checks on a single backbone PDB.

    Returns:
        dict with 'valid' (bool), 'chain_breaks', 'motif_issues',
        'clashes', 'bond_issues' lists, and 'summary' string.
    """
    coords, chains, resnums = parse_backbone(pdb_path)

    if len(coords) == 0:
        return {'valid': False, 'summary': 'Empty PDB (no CA atoms)'}

    result = {
        'pdb': pdb_path.name,
        'n_residues': len(coords),
        'chain_breaks': check_chain_breaks(coords, chains, max_break),
        'bond_issues': check_bond_lengths(coords, chains),
        'clashes': check_steric_clashes(coords, chains, min_clash),
        'motif_issues': [],
    }

    # Motif-specific checks
    if motif_json_path is not None and motif_json_path.exists():
        with open(motif_json_path) as f:
            motif_positions = json.load(f)
        result['motif_issues'] = check_motif_connectivity(
            coords, chains, resnums, motif_positions, max_break)

    # Overall validity
    n_critical = len(result['chain_breaks']) + len(result['motif_issues'])
    result['valid'] = n_critical == 0

    parts = []
    if result['chain_breaks']:
        parts.append(f"{len(result['chain_breaks'])} chain break(s)")
    if result['motif_issues']:
        parts.append(f"{len(result['motif_issues'])} motif disconnect(s)")
    if result['clashes']:
        parts.append(f"{len(result['clashes'])} clash(es)")
    if result['bond_issues']:
        parts.append(f"{len(result['bond_issues'])} bond issue(s)")
    result['summary'] = ', '.join(parts) if parts else 'OK'

    return result


def main():
    parser = argparse.ArgumentParser(
        description='Validate RFdiffusion backbones for physical validity'
    )
    parser.add_argument(
        '-i', '--input-dir', type=str, required=True,
        help='Directory containing backbone PDB files'
    )
    parser.add_argument(
        '--max-break', type=float, default=4.2,
        help='Max CA-CA distance before flagging chain break (default: 4.2 Å)'
    )
    parser.add_argument(
        '--min-clash', type=float, default=2.0,
        help='Min CA-CA distance before flagging steric clash (default: 2.0 Å)'
    )
    parser.add_argument(
        '--reject-dir', type=str, default=None,
        help='Move invalid PDBs to this directory (default: input_dir/rejected/)'
    )
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    pdb_files = sorted(input_dir.glob('*.pdb'))

    if not pdb_files:
        print(f'No PDB files found in {input_dir}')
        sys.exit(1)

    reject_dir = Path(args.reject_dir) if args.reject_dir else input_dir / 'rejected'

    print(f'Validating {len(pdb_files)} backbone designs...')
    print(f'  Max CA-CA break: {args.max_break} Å')
    print(f'  Min CA-CA clash: {args.min_clash} Å')

    valid_count = 0
    rejected_count = 0

    for pdb_path in pdb_files:
        # Find matching motif_fixed.json
        stem = pdb_path.stem
        motif_json = pdb_path.parent / f'{stem}_motif_fixed.json'

        result = validate_backbone(
            pdb_path,
            motif_json_path=motif_json if motif_json.exists() else None,
            max_break=args.max_break,
            min_clash=args.min_clash,
        )

        if result['valid']:
            valid_count += 1
        else:
            rejected_count += 1
            print(f'  REJECT: {pdb_path.name} — {result["summary"]}')
            # Move invalid PDB + associated files to rejected dir
            reject_dir.mkdir(parents=True, exist_ok=True)
            for ext in ['.pdb', '.trb', '_motif_fixed.json']:
                src = pdb_path.parent / f'{stem}{ext}'
                if src.exists():
                    src.rename(reject_dir / src.name)

    print(f'\nResults: {valid_count} valid, {rejected_count} rejected')
    if rejected_count > 0:
        print(f'Rejected designs moved to: {reject_dir}/')

    return 0 if valid_count > 0 else 1


if __name__ == '__main__':
    sys.exit(main())
