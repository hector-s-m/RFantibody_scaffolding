#!/usr/bin/env python3
"""
Extract confidence and interface metrics from Boltz2 prediction outputs.

Metrics computed per design:
  From Boltz2 JSON:  ipTM, iPAE (complex_ipae), pLDDT, iPLDDT
  From PAE NPZ:      ipSAE, pDockQ, pDockQ2, LIS
  From structures:   binder RMSD (designed vs Boltz2-predicted)

Usage:
    python scripts/extract_boltz2_metrics.py \
        -i boltz2_output/ \
        -d mpnn_designs/ \
        -o boltz2_metrics.csv

Runs in the RFantibody conda environment (needs numpy + pandas only).

References:
  ipSAE:   Dunbrack. https://www.biorxiv.org/content/10.1101/2025.02.10.637595v2
  pDockQ:  Bryant, Pozotti, Elofsson. https://doi.org/10.1038/s41467-022-28865-w
  pDockQ2: Zhu, Shenoy, Kundrotas, Elofsson. https://doi.org/10.1093/bioinformatics/btad424
  LIS:     Kim, Hu, Comjean, et al. https://doi.org/10.1101/2024.02.19.580970
"""

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# =============================================================================
# Boltz2 output discovery
# =============================================================================

def find_predictions(output_dir: Path) -> list[tuple[str, Path]]:
    """Find Boltz2 prediction directories.

    Returns list of (stem_name, prediction_directory).
    """
    predictions = []

    for pattern in [
        'boltz_results_*/predictions/*',
        'predictions/*',
    ]:
        for pred_dir in sorted(output_dir.glob(pattern)):
            if pred_dir.is_dir():
                predictions.append((pred_dir.name, pred_dir))
        if predictions:
            break

    # Fallback: directories containing confidence JSON
    if not predictions:
        for json_file in sorted(output_dir.glob('*/confidence_*.json')):
            pred_dir = json_file.parent
            predictions.append((pred_dir.name, pred_dir))

    return predictions


# =============================================================================
# File loaders
# =============================================================================

def load_confidence_json(pred_dir: Path, stem: str) -> dict | None:
    """Load Boltz2 confidence JSON."""
    candidates = [pred_dir / f'confidence_{stem}_model_0.json']
    candidates.extend(pred_dir.glob('confidence_*_model_0.json'))
    for path in candidates:
        if path.exists():
            with open(path) as f:
                return json.load(f)
    return None


def load_npz(pred_dir: Path, prefix: str, stem: str, key: str) -> np.ndarray | None:
    """Load a Boltz2 NPZ file."""
    candidates = [pred_dir / f'{prefix}_{stem}_model_0.npz']
    candidates.extend(pred_dir.glob(f'{prefix}_*_model_0.npz'))
    for path in candidates:
        if path.exists():
            data = np.load(str(path))
            if key in data:
                return data[key]
    return None


def find_predicted_structure(pred_dir: Path, stem: str) -> Path | None:
    """Find the Boltz2 predicted structure (PDB or mmCIF)."""
    for ext in ['.pdb', '.cif']:
        path = pred_dir / f'{stem}_model_0{ext}'
        if path.exists():
            return path
    # Glob fallback
    for ext in ['*.pdb', '*.cif']:
        matches = list(pred_dir.glob(f'*_model_0{ext[1:]}'))
        if matches:
            return matches[0]
    return None


# =============================================================================
# PDB parsing utilities
# =============================================================================

def parse_ca_coords(pdb_path: Path) -> tuple[np.ndarray, list[str], list[int]]:
    """Extract CA coordinates, chain IDs, and residue numbers from PDB.

    Returns:
        (coords [N,3], chain_ids [N], resnums [N])
    """
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
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords.append([x, y, z])
            chains.append(chain)
            resnums.append(resnum)
    return np.array(coords), chains, resnums


def parse_cb_coords(pdb_path: Path) -> tuple[np.ndarray, list[str]]:
    """Extract CB coordinates (CA for GLY) and chain IDs.

    Returns:
        (coords [N,3], chain_ids [N])
    """
    residue_coords = {}
    residue_chains = {}
    order = []
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            atom = line[12:16].strip()
            chain = line[21].strip()
            resnum = int(line[22:26].strip())
            key = (chain, resnum)
            coord = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            if key not in residue_coords:
                order.append(key)
                residue_chains[key] = chain
            if atom == 'CB':
                residue_coords[key] = coord
            elif atom == 'CA' and key not in residue_coords:
                residue_coords[key] = coord
    coords = np.array([residue_coords[k] for k in order])
    chains = [residue_chains[k] for k in order]
    return coords, chains


def parse_cif_ca_coords(cif_path: Path) -> tuple[np.ndarray, list[str], list[int]]:
    """Extract CA coordinates from mmCIF file.

    Returns:
        (coords [N,3], chain_ids [N], resnums [N])
    """
    coords, chains, resnums = [], [], []
    seen = set()
    with open(cif_path) as f:
        for line in f:
            if not line.startswith('ATOM') and not line.startswith('HETATM'):
                continue
            parts = line.split()
            if len(parts) < 15:
                continue
            # _atom_site columns: group_PDB id type_symbol label_atom_id ...
            # Typical mmCIF: ATOM 1 N N . ALA A 1 1 ? x y z ...
            atom_name = parts[3]
            if atom_name != 'CA':
                continue
            chain = parts[6]
            resnum = int(parts[8]) if parts[8].isdigit() else int(parts[7])
            key = (chain, resnum)
            if key in seen:
                continue
            seen.add(key)
            x = float(parts[10])
            y = float(parts[11])
            z = float(parts[12])
            coords.append([x, y, z])
            chains.append(chain)
            resnums.append(resnum)
    return np.array(coords) if coords else np.zeros((0, 3)), chains, resnums


# =============================================================================
# Binder RMSD
# =============================================================================

def compute_binder_rmsd(
    designed_pdb: Path,
    predicted_pdb: Path,
) -> float | None:
    """Compute CA RMSD of binder chains between designed and predicted structures.

    Auto-detects binder chains: uses all chains EXCEPT the last chain
    (assumed to be the target). Handles both HLT and ABC naming.

    Args:
        designed_pdb: PDB from AntiBMPNN (input to Boltz2).
        predicted_pdb: Structure predicted by Boltz2.

    Returns:
        RMSD in Angstroms, or None if structures can't be aligned.
    """
    # Parse designed structure
    if str(designed_pdb).endswith('.cif'):
        des_coords, des_chains, _ = parse_cif_ca_coords(designed_pdb)
    else:
        des_coords, des_chains, _ = parse_ca_coords(designed_pdb)

    # Parse predicted structure
    if str(predicted_pdb).endswith('.cif'):
        pred_coords, pred_chains, _ = parse_cif_ca_coords(predicted_pdb)
    else:
        pred_coords, pred_chains, _ = parse_ca_coords(predicted_pdb)

    if len(des_coords) == 0 or len(pred_coords) == 0:
        return None

    # Auto-detect binder chains: first chain(s) = target, remaining = binder
    # Convention: target chains listed first (A,B,...), binder chains after
    # For designed PDB (HLT): T=target, H/L=binder
    # For predicted PDB (ABC): A=target, B/C=binder
    des_unique = list(dict.fromkeys(des_chains))    # preserve order
    pred_unique = list(dict.fromkeys(pred_chains))

    # Designed PDB: last chain = target (HLT convention), rest = binder
    des_binder_chains = set(des_unique[:-1]) if len(des_unique) > 1 else set()
    # Predicted PDB: first chain = target (remapped convention), rest = binder
    pred_binder_chains = set(pred_unique[1:]) if len(pred_unique) > 1 else set()

    if not des_binder_chains or not pred_binder_chains:
        return None  # Cannot compute binder RMSD without binder chains

    des_mask = np.array([c in des_binder_chains for c in des_chains])
    pred_mask = np.array([c in pred_binder_chains for c in pred_chains])

    des_binder = des_coords[des_mask]
    pred_binder = pred_coords[pred_mask]

    if len(des_binder) == 0 or len(pred_binder) == 0:
        return None

    # Align lengths (take minimum)
    n = min(len(des_binder), len(pred_binder))
    des_binder = des_binder[:n]
    pred_binder = pred_binder[:n]

    # RMSD (no superposition — direct positional RMSD)
    diff = des_binder - pred_binder
    rmsd = float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))
    return rmsd


# =============================================================================
# Motif RMSD (designed vs Boltz2-predicted motif positions)
# =============================================================================

def compute_motif_rmsd(
    designed_pdb: Path,
    predicted_pdb: Path,
    motif_fixed_json: Path,
) -> float | None:
    """Compute CA RMSD of motif residues between designed and predicted structures.

    Uses the motif_fixed.json (from RFdiffusion) to identify which residues
    are motif positions. Computes RMSD only for those positions.

    Args:
        designed_pdb: PDB from AntiBMPNN (input to Boltz2).
        predicted_pdb: Structure predicted by Boltz2.
        motif_fixed_json: JSON file with motif fixed positions
            (e.g., {"H": [45, 46, 47, 48]}).

    Returns:
        RMSD in Angstroms, or None if computation fails.
    """
    with open(motif_fixed_json) as f:
        motif_positions = json.load(f)

    # Parse designed structure
    if str(designed_pdb).endswith('.cif'):
        des_coords, des_chains, des_resnums = parse_cif_ca_coords(designed_pdb)
    else:
        des_coords, des_chains, des_resnums = parse_ca_coords(designed_pdb)

    # Parse predicted structure
    if str(predicted_pdb).endswith('.cif'):
        pred_coords, pred_chains, pred_resnums = parse_cif_ca_coords(predicted_pdb)
    else:
        pred_coords, pred_chains, pred_resnums = parse_ca_coords(predicted_pdb)

    # Auto-detect chain mapping between designed and predicted structures.
    # Designed PDB: HLT order (binder H,L first, target T last)
    # Predicted PDB: target first (A), binder after (B,C,...)
    # Map by role: designed binder chains → predicted binder chains,
    #              designed target chains → predicted target chains
    des_unique = list(dict.fromkeys(des_chains))
    pred_unique = list(dict.fromkeys(pred_chains))

    # Designed: last chain(s) = target, rest = binder
    des_binder = des_unique[:-1] if len(des_unique) > 1 else des_unique
    des_target = des_unique[-1:] if len(des_unique) > 1 else []
    # Predicted: first chain(s) = target, rest = binder
    pred_target = pred_unique[:len(des_target)]
    pred_binder = pred_unique[len(des_target):]

    chain_remap = {}
    for old, new in zip(des_target, pred_target):
        chain_remap[old] = new
    for old, new in zip(des_binder, pred_binder):
        chain_remap[old] = new

    des_motif_coords = []
    pred_motif_coords = []

    for chain, positions in motif_positions.items():
        # Find motif residues in designed structure (uses original chain IDs)
        for pos in positions:
            for i, (c, r) in enumerate(zip(des_chains, des_resnums)):
                if c == chain and r == pos:
                    des_motif_coords.append(des_coords[i])
                    break

        # Find motif residues in predicted structure (use auto-detected mapping)
        pred_chain = chain_remap.get(chain, chain)
        for pos in positions:
            for i, (c, r) in enumerate(zip(pred_chains, pred_resnums)):
                if c == pred_chain and r == pos:
                    pred_motif_coords.append(pred_coords[i])
                    break

    total_expected = sum(len(pos) for pos in motif_positions.values())
    if len(des_motif_coords) == 0 or len(pred_motif_coords) == 0:
        print(f'  Warning: Motif RMSD skipped — found {len(des_motif_coords)}/{total_expected} '
              f'designed, {len(pred_motif_coords)}/{total_expected} predicted motif residues')
        return None

    if len(des_motif_coords) != len(pred_motif_coords):
        print(f'  Warning: Motif residue count mismatch — designed={len(des_motif_coords)}, '
              f'predicted={len(pred_motif_coords)} (expected {total_expected})')

    n = min(len(des_motif_coords), len(pred_motif_coords))
    des_arr = np.array(des_motif_coords[:n])
    pred_arr = np.array(pred_motif_coords[:n])

    diff = des_arr - pred_arr
    rmsd = float(np.sqrt(np.mean(np.sum(diff ** 2, axis=1))))
    return rmsd


# =============================================================================
# ipSAE calculation (Dunbrack method)
# =============================================================================

def _ptm_func(x, d0):
    """PTM-style score: 1 / (1 + (x/d0)^2)."""
    return 1.0 / (1.0 + (x / d0) ** 2.0)


def _calc_d0(L, pair_type='protein'):
    """d0 from Yang & Skolnick (2004)."""
    min_val = 2.0 if pair_type == 'nucleic_acid' else 1.0
    if L <= 27:
        return min_val
    return max(min_val, 1.24 * (float(L) - 15) ** (1.0 / 3.0) - 1.8)


def calculate_ipsae(pae_matrix, chain_ids, pae_cutoff=10.0):
    """Calculate ipSAE for all chain pairs.

    Args:
        pae_matrix: (N, N) PAE matrix from Boltz2.
        chain_ids: list of chain IDs per token (length N).
        pae_cutoff: PAE cutoff for valid interactions.

    Returns:
        Dict of {f'{ch1}_{ch2}': score, ...} plus 'ipSAE' (max over all pairs).
    """
    chain_ids = np.array(chain_ids)
    unique_chains = list(dict.fromkeys(chain_ids))  # preserve order
    results = {}

    for ch1 in unique_chains:
        for ch2 in unique_chains:
            if ch1 == ch2:
                continue
            mask1 = chain_ids == ch1
            mask2 = chain_ids == ch2
            sub_pae = pae_matrix[np.ix_(mask1, mask2)]
            valid = sub_pae < pae_cutoff

            scores = []
            for i in range(sub_pae.shape[0]):
                v = valid[i]
                if v.any():
                    d0 = _calc_d0(int(v.sum()))
                    scores.append(float(_ptm_func(sub_pae[i][v], d0).mean()))
                else:
                    scores.append(0.0)
            results[f'{ch1}_{ch2}'] = float(max(scores)) if scores else 0.0

    results['ipSAE'] = max(results.values()) if results else 0.0
    return results


# =============================================================================
# pDockQ, pDockQ2, LIS
# =============================================================================

def calculate_contact_scores(
    pdb_path: Path,
    pae_matrix: np.ndarray,
    plddt_array: np.ndarray,
    chain_ids: list[str],
    binder_chains: set[str] = None,
    target_chains: set[str] = None,
    contact_dist: float = 8.0,
) -> dict:
    """Calculate pDockQ, pDockQ2, and LIS.

    Args:
        pdb_path: Path to predicted structure (PDB format).
        pae_matrix: Full PAE matrix (N, N).
        plddt_array: Per-residue pLDDT (N,), 0-1 scale.
        chain_ids: Chain ID per token (length N).
        binder_chains: Binder chain IDs (default: H, L).
        target_chains: Target chain IDs (default: T, C).
        contact_dist: CB distance cutoff for contacts (default: 8 A).

    Returns:
        Dict with 'pDockQ', 'pDockQ2', 'LIS'.
    """
    chain_ids = np.array(chain_ids)
    # Auto-detect binder/target from predicted structure chain order:
    # Convention: first chain(s) = target, remaining = binder
    unique_chains = list(dict.fromkeys(chain_ids))
    if binder_chains is None:
        binder_chains = set(unique_chains[1:]) if len(unique_chains) > 1 else set()
    if target_chains is None:
        target_chains = {unique_chains[0]} if unique_chains else set()

    binder_mask = np.array([c in binder_chains for c in chain_ids])
    target_mask = np.array([c in target_chains for c in chain_ids])

    # Need at least some residues on each side
    if not binder_mask.any() or not target_mask.any():
        return {'pDockQ': 0.0, 'pDockQ2': 0.0, 'LIS': 0.0}

    binder_idx = np.where(binder_mask)[0]
    target_idx = np.where(target_mask)[0]

    # --- LIS (PAE-only) ---
    pae_bt = pae_matrix[np.ix_(binder_mask, target_mask)].flatten()
    pae_tb = pae_matrix[np.ix_(target_mask, binder_mask)].flatten()
    valid_bt = pae_bt[pae_bt < 12.0]
    valid_tb = pae_tb[pae_tb < 12.0]
    lis_bt = float(np.mean((12.0 - valid_bt) / 12.0)) if valid_bt.size > 0 else 0.0
    lis_tb = float(np.mean((12.0 - valid_tb) / 12.0)) if valid_tb.size > 0 else 0.0
    if valid_bt.size > 0 and valid_tb.size > 0:
        lis = (lis_bt + lis_tb) / 2.0
    else:
        lis = lis_bt + lis_tb  # one is 0.0

    # --- Contact-based scores (need CB coords) ---
    cb_coords, cb_chains = parse_cb_coords(pdb_path)
    cb_chains = np.array(cb_chains)

    cb_binder = cb_coords[np.array([c in binder_chains for c in cb_chains])]
    cb_target = cb_coords[np.array([c in target_chains for c in cb_chains])]

    if len(cb_binder) == 0 or len(cb_target) == 0:
        return {'pDockQ': 0.0, 'pDockQ2': 0.0, 'LIS': round(lis, 4)}

    # Distance matrix (target x binder)
    diff = cb_target[:, None, :] - cb_binder[None, :, :]
    dist_matrix = np.sqrt((diff ** 2).sum(axis=2))
    contacts = dist_matrix < contact_dist
    npairs = int(contacts.sum())

    if npairs == 0:
        return {'pDockQ': 0.0, 'pDockQ2': 0.0, 'LIS': round(lis, 4)}

    # Interface residues
    target_contact_local = np.where(contacts.any(axis=1))[0]
    binder_contact_local = np.where(contacts.any(axis=0))[0]

    # Map local indices back to global for pLDDT lookup
    target_global = target_idx[:len(cb_target)]
    binder_global = binder_idx[:len(cb_binder)]
    interface_global = np.concatenate([
        target_global[target_contact_local],
        binder_global[binder_contact_local],
    ])

    plddt_100 = plddt_array * 100.0 if plddt_array.max() <= 1.0 else plddt_array
    mean_plddt = float(plddt_100[interface_global].mean())

    # pDockQ
    x_pdockq = mean_plddt * math.log10(max(npairs, 1))
    pdockq = 0.724 / (1.0 + math.exp(-0.052 * (x_pdockq - 152.611))) + 0.018

    # pDockQ2 (max of both directions)
    def _pdockq2_dir(row_global, col_global, contact_mat):
        ri, ci = np.where(contact_mat)
        if len(ri) == 0:
            return 0.0
        pae_vals = pae_matrix[row_global[ri], col_global[ci]]
        mean_ptm = float(_ptm_func(pae_vals, 10.0).mean())
        x = mean_plddt * mean_ptm
        return 1.31 / (1.0 + math.exp(-0.075 * (x - 84.733))) + 0.005

    pdockq2_fwd = _pdockq2_dir(target_global[:len(cb_target)],
                                binder_global[:len(cb_binder)], contacts)
    pdockq2_rev = _pdockq2_dir(binder_global[:len(cb_binder)],
                                target_global[:len(cb_target)], contacts.T)
    pdockq2 = max(pdockq2_fwd, pdockq2_rev)

    return {
        'pDockQ': round(pdockq, 4),
        'pDockQ2': round(pdockq2, 4),
        'LIS': round(lis, 4),
    }


# =============================================================================
# Per-prediction metric extraction
# =============================================================================

def get_chain_ids_from_pae(pae_matrix: np.ndarray, pred_dir: Path, stem: str) -> list[str] | None:
    """Get chain IDs by parsing the predicted structure."""
    struct_path = find_predicted_structure(pred_dir, stem)
    if struct_path is None:
        return None
    if str(struct_path).endswith('.cif'):
        _, chains, _ = parse_cif_ca_coords(struct_path)
    else:
        _, chains, _ = parse_ca_coords(struct_path)
    return chains


def extract_metrics_for_prediction(
    stem: str,
    pred_dir: Path,
    designed_pdb_dir: Path | None = None,
) -> dict | None:
    """Extract all metrics for one Boltz2 prediction."""

    confidence = load_confidence_json(pred_dir, stem)
    if confidence is None:
        print(f'  Warning: No confidence JSON for {stem}')
        return None

    result = {'design': stem}

    # --- Boltz2 JSON metrics ---
    result['ipTM'] = confidence.get('iptm')
    result['pTM'] = confidence.get('ptm')
    result['pLDDT'] = confidence.get('complex_plddt')
    result['iPLDDT'] = confidence.get('complex_iplddt')
    result['confidence_score'] = confidence.get('confidence_score')

    # iPAE from JSON (if available) or compute from NPZ
    result['iPAE'] = confidence.get('complex_ipae')

    # --- PAE / pLDDT from NPZ ---
    pae = load_npz(pred_dir, 'pae', stem, 'pae')
    plddt = load_npz(pred_dir, 'plddt', stem, 'plddt')

    # Get chain IDs from predicted structure
    struct_path = find_predicted_structure(pred_dir, stem)
    chain_ids = None
    if struct_path is not None:
        if str(struct_path).endswith('.cif'):
            _, chain_ids, _ = parse_cif_ca_coords(struct_path)
        else:
            _, chain_ids, _ = parse_ca_coords(struct_path)

    # --- ipSAE ---
    if pae is not None and chain_ids is not None and len(chain_ids) == pae.shape[0]:
        try:
            ipsae_results = calculate_ipsae(pae, chain_ids)
            result['ipSAE'] = round(ipsae_results.get('ipSAE', 0.0), 4)
        except Exception as e:
            print(f'  Warning: ipSAE failed for {stem}: {e}')
            result['ipSAE'] = None
    else:
        result['ipSAE'] = None

    # --- pDockQ, pDockQ2, LIS ---
    if (pae is not None and plddt is not None
            and struct_path is not None and chain_ids is not None
            and len(chain_ids) == pae.shape[0]):
        try:
            # Truncate plddt if needed
            plddt_trunc = plddt[:len(chain_ids)]
            contact_scores = calculate_contact_scores(
                struct_path, pae, plddt_trunc, chain_ids,
            )
            result['pDockQ'] = contact_scores['pDockQ']
            result['pDockQ2'] = contact_scores['pDockQ2']
            result['LIS'] = contact_scores['LIS']
        except Exception as e:
            print(f'  Warning: Contact scores failed for {stem}: {e}')
            result['pDockQ'] = None
            result['pDockQ2'] = None
            result['LIS'] = None
    else:
        result['pDockQ'] = None
        result['pDockQ2'] = None
        result['LIS'] = None

    # --- Binder RMSD (designed vs predicted) ---
    result['binder_RMSD'] = None
    result['motif_RMSD'] = None
    if designed_pdb_dir is not None and struct_path is not None:
        # Try to find matching designed PDB
        designed_pdb = designed_pdb_dir / f'{stem}.pdb'
        if not designed_pdb.exists():
            # Try without _model_0 suffix or other variants
            for candidate in designed_pdb_dir.glob(f'{stem}*.pdb'):
                designed_pdb = candidate
                break
        if designed_pdb.exists():
            try:
                rmsd = compute_binder_rmsd(designed_pdb, struct_path)
                if rmsd is not None:
                    result['binder_RMSD'] = round(rmsd, 3)
            except Exception as e:
                print(f'  Warning: RMSD failed for {stem}: {e}')

            # --- Motif RMSD ---
            # Look for motif_fixed.json from RFdiffusion output
            # New naming: PREFIX_RFN_mpnnM → RF file is PREFIX_RFN
            # Old naming: ab_des_N_dldesign_M → ab_des_N
            # Try new naming first, then fall back to old
            import re
            motif_json = None
            rfdiff_dir = designed_pdb_dir.parent / 'RFdiffusion_backbones'
            designs_dir = designed_pdb_dir.parent / 'designs'

            # New convention: PREFIX_RFN_mpnnM → PREFIX_RFN_motif_fixed.json
            m = re.match(r'^(.+_RF\d+)_mpnn\d+$', stem)
            if m and rfdiff_dir.exists():
                rf_base = m.group(1)  # e.g. PROTEIN_Nb_RF0
                candidate = rfdiff_dir / f'{rf_base}_motif_fixed.json'
                if candidate.exists():
                    motif_json = candidate

            # Old convention fallback: ab_des_N_dldesign_M → designs/ab_des_N_motif_fixed.json
            if motif_json is None:
                base_design = stem.split('_dldesign_')[0] if '_dldesign_' in stem else stem
                if designs_dir.exists():
                    candidate = designs_dir / f'{base_design}_motif_fixed.json'
                    if candidate.exists():
                        motif_json = candidate
            if motif_json is not None:
                try:
                    mrmsd = compute_motif_rmsd(designed_pdb, struct_path, motif_json)
                    if mrmsd is not None:
                        result['motif_RMSD'] = round(mrmsd, 3)
                except Exception as e:
                    print(f'  Warning: Motif RMSD failed for {stem}: {e}')

    return result


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Extract Boltz2 metrics: ipTM, iPAE, pLDDT, iPLDDT, '
                    'ipSAE, pDockQ, pDockQ2, LIS, binder RMSD, motif RMSD'
    )
    parser.add_argument(
        '-i', '--input-dir', type=str, required=True,
        help='Boltz2 output directory (containing boltz_results_*/)'
    )
    parser.add_argument(
        '-d', '--designed-pdb-dir', type=str, default=None,
        help='Directory with designed PDBs (from AntiBMPNN) for RMSD calculation'
    )
    parser.add_argument(
        '-o', '--output-csv', type=str, default='boltz2_metrics.csv',
        help='Output CSV path (default: boltz2_metrics.csv)'
    )
    parser.add_argument(
        '--rank-by', type=str, default='ipTM',
        help='Metric to rank by, descending (default: ipTM)'
    )
    args = parser.parse_args()

    output_dir = Path(args.input_dir)
    designed_dir = Path(args.designed_pdb_dir) if args.designed_pdb_dir else None

    predictions = find_predictions(output_dir)
    if not predictions:
        print(f'Error: No Boltz2 predictions found in {output_dir}')
        sys.exit(1)

    print(f'Found {len(predictions)} predictions in {output_dir}')
    if designed_dir:
        print(f'Designed PDBs for RMSD: {designed_dir}')

    results = []
    failed = 0
    for stem, pred_dir in predictions:
        try:
            metrics = extract_metrics_for_prediction(stem, pred_dir, designed_dir)
            if metrics is not None:
                results.append(metrics)
            else:
                failed += 1
        except Exception as e:
            print(f'  Error processing {stem}: {e}')
            failed += 1

    if not results:
        print('Error: No metrics extracted.')
        sys.exit(1)

    df = pd.DataFrame(results)

    # Rank
    if args.rank_by in df.columns:
        df = df.sort_values(args.rank_by, ascending=False, na_position='last')

    csv_path = Path(args.output_csv)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(csv_path, index=False, float_format='%.4f')

    print(f'\nMetrics written to {csv_path}')
    print(f'  Successful: {len(results)}, Failed: {failed}')

    # Print top 10
    display_cols = ['design']
    for col in ['ipTM', 'ipSAE', 'pDockQ', 'pDockQ2', 'LIS', 'pLDDT', 'iPLDDT', 'iPAE', 'binder_RMSD', 'motif_RMSD']:
        if col in df.columns:
            display_cols.append(col)

    if len(df) > 0:
        print(f'\nTop 10 designs (ranked by {args.rank_by}):')
        print(df[display_cols].head(10).to_string(index=False))


if __name__ == '__main__':
    main()
