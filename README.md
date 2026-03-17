# RFantibody Scaffolding

Structure-based *de novo* antibody design with aromatic motif scaffolding.

Fork of [RFantibody](https://github.com/RosettaCommons/RFantibody) adding support for fixing aromatic peptide motifs (PHE/TYR/TRP) into CDR loops during antibody design.

## Pipeline

1. **RFdiffusion** — Design antibody backbone with fixed motif coordinates in CDR loop
2. **AntiBMPNN** — Antibody-finetuned ProteinMPNN for CDR sequence design (motif residues stay fixed)
3. **RF2** — Predict/refine final structures

## Requirements

- NVIDIA GPU with CUDA 11.8+
- Linux (Ubuntu 22.04 recommended)

## Setup

```bash
# 1. Clone
git clone https://github.com/hector-s-m/RFantibody_scaffolding.git
cd RFantibody_scaffolding

# 2. Install uv (if not already installed)
curl -LsSf https://astral.sh/uv/install.sh | sh

# 3. Download model weights (RFdiffusion, RF2)
bash include/download_weights.sh

# 4. Download AntiBMPNN weights (antibody-finetuned ProteinMPNN)
bash scripts/setup_antibmpnn.sh

# 5. Install dependencies
uv sync

# 6. Verify
uv run rfdiffusion --help
```

## Quick Start

### Standard Antibody Design

```bash
uv run rfdiffusion \
    -t target.pdb \
    -f framework.pdb \
    -o designs/ab \
    -n 100 \
    -l "H1:7,H2:6,H3:5-13" \
    -h "B146,B170,B177"
```

### Motif Scaffolding

Provide a combined PDB with target (chain A) + aromatic motif (last chain). The target is extracted automatically — no separate `--target` needed.

```bash
# Step 1: RFdiffusion — design backbones with motif fixed in H3
uv run rfdiffusion \
    -f framework.pdb \
    -m motif_combined.pdb \
    -o designs/motif_ab \
    -n 100 \
    -l "H1:,H2:,H3:10-16,L1:,L2:,L3:" \
    --motif-cdr H3

# Step 2: AntiBMPNN — design sequences (motif residues stay fixed)
# Uses AntiBMPNN weights automatically if installed, else falls back to vanilla ProteinMPNN
uv run proteinmpnn \
    -q designs.qv \
    --output-quiver sequences.qv \
    -n 4 -t 0.2

# Step 3: RF2 — predict structures
uv run rf2 \
    -q sequences.qv \
    --output-quiver predictions.qv \
    -r 10
```

A complete example script is at `scripts/examples/motif_scaffolding_pipeline.sh`.

## Input Files

### Framework PDB

An HLT-formatted antibody or nanobody scaffold. Chains labeled H (heavy), L (light), T (target). CDR loop positions annotated via REMARK lines.

Convert from Chothia format:
```bash
python scripts/util/chothia_to_HLT.py -inpdb mychothia.pdb -outpdb myHLT.pdb
```

Provided frameworks:
- Nanobody: `scripts/examples/example_inputs/h-NbBCII10.pdb`
- ScFv: `scripts/examples/example_inputs/hu-4D5-8_Fv.pdb`

### Combined Motif PDB (for motif scaffolding)

A PDB file containing:
- **Chain A** (or earlier chains): Target protein
- **Last chain**: Aromatic peptide motif (2-6+ residues with PHE/TYR/TRP)

The motif chain must have full backbone coordinates (N, CA, C, O). Sidechain atoms are preserved. Rosetta PDB3 hydrogen naming is handled automatically.

## CLI Reference

### RFdiffusion

| Flag | Description |
|------|-------------|
| `-t, --target` | Target PDB (optional if `--motif` provided) |
| `-f, --framework` | Framework PDB (required) |
| `-m, --motif` | Combined motif+target PDB |
| `--motif-cdr` | CDR loop for motif (default: H3) |
| `-o, --output` | Output prefix |
| `-q, --output-quiver` | Output Quiver file |
| `-n, --num-designs` | Number of designs |
| `-l, --design-loops` | Loop lengths, e.g., `"H3:10-16"` |
| `-h, --hotspots` | Target hotspot residues |
| `--deterministic` | Reproducible results |

### AntiBMPNN / ProteinMPNN

Uses AntiBMPNN weights by default if installed (via `scripts/setup_antibmpnn.sh`). Falls back to vanilla ProteinMPNN otherwise. Override with `-w` to use specific weights.

| Flag | Description |
|------|-------------|
| `-i, --input-dir` | Input PDB directory |
| `-q, --input-quiver` | Input Quiver file |
| `-o, --output-dir` | Output PDB directory |
| `--output-quiver` | Output Quiver file |
| `-n, --seqs-per-struct` | Sequences per structure |
| `-t, --temperature` | Sampling temperature |
| `-l, --loops` | Loops to design (default: all CDRs) |
| `-w, --weights` | Override model weights (e.g., vanilla ProteinMPNN) |

### RF2

| Flag | Description |
|------|-------------|
| `-p, --input-pdb` | Single input PDB |
| `-i, --input-dir` | Input PDB directory |
| `-q, --input-quiver` | Input Quiver file |
| `-o, --output-dir` | Output PDB directory |
| `--output-quiver` | Output Quiver file |
| `-r, --num-recycles` | Recycling iterations (default: 10) |

All commands support `--help` for full options.

## Design Loop Syntax

```
-l "H1:7,H2:6,H3:10-16,L1:,L2:,L3:"
```

- `H3:10-16` — sample H3 length uniformly from 10 to 16
- `H1:7` — fix H1 at length 7
- `L1:` — design L1 but keep framework length
- Omitted loops are not designed (sequence + structure fixed)

For motif scaffolding, ensure the motif CDR range can fit `motif_length + 2` (minimum 1 flank each side).

## Quiver Files

Efficient storage for large design campaigns. Key commands:

```bash
qvls my.qv                    # List designs
qvextract my.qv -o pdbs/      # Extract all PDBs
qvscorefile my.qv              # Export scores to TSV
qvfrompdbs *.pdb > my.qv       # Create from PDBs
```

## How Motif Scaffolding Works

1. The combined PDB is parsed: last chain = motif, earlier chains = target
2. CDR loop length is sampled from `--design-loops` range (e.g., H3:10-16)
3. Flanking residues are computed: `flanks = loop_length - motif_length`, split evenly with random ±1 offset
4. Motif backbone coordinates are placed at the computed position within the CDR loop
5. During diffusion, motif residues are masked (not denoised) — coordinates stay fixed
6. During AntiBMPNN/ProteinMPNN, motif residue identities are fixed — only flanking positions are redesigned
7. RF2 predicts the final structure; motif positions can be verified by RMSD

## Filtering

Recommended criteria:
- RF2 pAE < 10
- RMSD (design vs RF2 prediction) < 2 A
- Optional: Rosetta ddG < -20

## License

MIT License. See LICENSE file.

Based on [RFantibody](https://www.biorxiv.org/content/10.1101/2024.03.14.585103v1) by Nate Bennett, Joe Watson, and the RFantibody Team.
