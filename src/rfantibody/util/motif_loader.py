"""
Motif loader for aromatic peptide motif scaffolding.

Loads and validates combined PDB files containing a target protein and an
aromatic peptide motif (2-6+ residues with PHE/TYR/TRP anchors). The motif
is the last chain in the PDB file; all preceding chains are the target.
"""

import numpy as np

from rfantibody.rfdiffusion.chemical import aa2long, aa2num, aa_321, seq2chars


# Aromatic residue types (integer indices in aa2num)
AROMATIC_AA = {'PHE', 'TYR', 'TRP'}
AROMATIC_IDX = {aa2num[aa] for aa in AROMATIC_AA}

# Pre-compute (aa_idx, atom_name) → atom_slot_index for O(1) lookup
# Replaces O(27) linear scan per atom line in _parse_chain_atoms()
_ATOM_SLOT_MAP = {}
for _aa_idx, _atoms in enumerate(aa2long):
    for _slot, _atm in enumerate(_atoms):
        if _atm is not None:
            _ATOM_SLOT_MAP[(_aa_idx, _atm.strip())] = _slot


def load_motif_combined_pdb(pdb_path: str) -> dict:
    """
    Parse a combined PDB file containing target chain(s) and a motif chain.

    The motif chain is the LAST chain in the file. All earlier chains are
    treated as the target. Handles Rosetta PDB3 hydrogen naming and
    skips energy table lines appended after TER.

    Returns:
        dict with keys:
            motif_xyz:    np.ndarray [M, 27, 3] — atom coordinates for motif residues
            motif_mask:   np.ndarray [M, 27]    — atom presence mask
            motif_seq:    np.ndarray [M]         — integer AA sequence
            motif_seq_str: str                   — 1-letter AA sequence
            motif_pdb_idx: list of (chain, resnum) tuples
            aromatic_indices: list[int]          — indices of PHE/TYR/TRP in motif
            motif_chain_id: str                  — chain letter of motif
            target_xyz:   np.ndarray [T, 27, 3]
            target_mask:  np.ndarray [T, 27]
            target_seq:   np.ndarray [T]
            target_pdb_idx: list of (chain, resnum) tuples
    """
    with open(pdb_path, 'r') as f:
        lines = f.readlines()

    # Collect all ATOM lines grouped by chain
    atom_lines_by_chain = {}
    chain_order = []
    for line in lines:
        if line[:4] != "ATOM":
            continue
        chain = line[21:22].strip()
        if chain not in atom_lines_by_chain:
            atom_lines_by_chain[chain] = []
            chain_order.append(chain)
        atom_lines_by_chain[chain].append(line)

    if len(chain_order) < 2:
        raise ValueError(
            f"Combined motif PDB must have at least 2 chains (target + motif), "
            f"found {len(chain_order)}: {chain_order}"
        )

    motif_chain_id = chain_order[-1]
    target_chain_ids = chain_order[:-1]

    # Parse each group
    motif_result = _parse_chain_atoms(atom_lines_by_chain[motif_chain_id])
    target_lines = []
    for tc in target_chain_ids:
        target_lines.extend(atom_lines_by_chain[tc])
    target_result = _parse_chain_atoms(target_lines)

    # Identify aromatic residues in motif
    aromatic_indices = [
        i for i, aa_idx in enumerate(motif_result['seq'])
        if aa_idx in AROMATIC_IDX
    ]

    # Build 1-letter sequence string using canonical mapping from chemical.py
    motif_seq_str = seq2chars(motif_result['seq'])

    return {
        'motif_xyz': motif_result['xyz'],
        'motif_mask': motif_result['mask'],
        'motif_seq': motif_result['seq'],
        'motif_seq_str': motif_seq_str,
        'motif_pdb_idx': motif_result['pdb_idx'],
        'aromatic_indices': aromatic_indices,
        'motif_chain_id': motif_chain_id,
        'target_xyz': target_result['xyz'],
        'target_mask': target_result['mask'],
        'target_seq': target_result['seq'],
        'target_pdb_idx': target_result['pdb_idx'],
    }


def validate_motif(motif_data: dict, min_length: int = 2) -> None:
    """
    Validate that the loaded motif is suitable for scaffolding.

    Raises ValueError if validation fails.
    """
    motif_len = len(motif_data['motif_seq'])

    if motif_len < min_length:
        raise ValueError(
            f"Motif too short ({motif_len} residues), minimum is {min_length}"
        )

    if not motif_data['aromatic_indices']:
        raise ValueError(
            "Motif must contain at least one aromatic residue (PHE/TYR/TRP). "
            f"Found sequence: {motif_data['motif_seq_str']}"
        )

    # Check backbone atoms (N, CA, C, O = indices 0-3) are present
    bb_mask = motif_data['motif_mask'][:, :4]
    if not bb_mask.all():
        missing = np.where(~bb_mask.all(axis=1))[0]
        raise ValueError(
            f"Motif residues {missing.tolist()} are missing backbone atoms"
        )

    # Check for NaN coordinates in backbone
    bb_xyz = motif_data['motif_xyz'][:, :4, :]
    if np.isnan(bb_xyz).any():
        raise ValueError("Motif contains NaN backbone coordinates")


def _parse_chain_atoms(atom_lines: list) -> dict:
    """
    Parse a list of ATOM lines into xyz, mask, seq, pdb_idx arrays.

    Uses the same 27-atom representation and aa2long mapping as
    rfantibody.rfdiffusion.parsers.HLT_pdb_parser.
    """
    # Identify unique residues by (chain, resnum) from CA atoms
    pdb_idx = []
    seq = []
    for line in atom_lines:
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        chain = line[21:22].strip()
        resnum = line[22:26].strip()
        aa3 = line[17:20].strip()
        pdb_idx.append((chain, resnum))
        seq.append(aa2num.get(aa3, 20))  # 20 = unknown

    n_res = len(pdb_idx)
    xyz = np.full((n_res, 27, 3), np.nan, dtype=np.float32)

    # Build O(1) lookup dict for residue indices
    pdb_idx_map = {key: i for i, key in enumerate(pdb_idx)}

    for line in atom_lines:
        chain = line[21:22].strip()
        resnum = line[22:26].strip()
        atom_name = ' ' + line[12:16].strip().ljust(3)
        aa3 = line[17:20].strip()

        key = (chain, resnum)
        idx = pdb_idx_map.get(key)
        if idx is None:
            continue

        aa_idx = aa2num.get(aa3, 20)
        if aa_idx >= len(aa2long):
            continue

        slot = _ATOM_SLOT_MAP.get((aa_idx, atom_name.strip()))
        if slot is not None:
            xyz[idx, slot, :] = [
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            ]

    mask = np.logical_not(np.isnan(xyz[..., 0]))
    xyz[np.isnan(xyz[..., 0])] = 0.0

    return {
        'xyz': xyz,
        'mask': mask,
        'seq': np.array(seq),
        'pdb_idx': pdb_idx,
    }
