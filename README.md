# RFantibody Scaffolding

Structure-based *de novo* antibody design with aromatic motif scaffolding.

Fork of [RFantibody](https://github.com/RosettaCommons/RFantibody) adding support for fixing aromatic peptide motifs (PHE/TYR/TRP) into CDR loops during antibody design.

## Pipeline

1. **RFdiffusion** — Design antibody backbone with fixed motif coordinates in CDR loop
2. **AntiBMPNN** — Antibody-finetuned ProteinMPNN for CDR sequence design (motif residues stay fixed)
3. **Boltz2** — Predict structures + confidence metrics (ipTM, pLDDT, PAE) for design ranking

## Requirements

- NVIDIA GPU with CUDA 11.8+ (L40S / A100 / H100 recommended)
- Linux (Ubuntu 22.04 recommended)
- Two conda environments: `RFantibody` (Steps 1-2) and `boltz_2.2.1` (Step 3)

## Setup

```bash
# 1. Clone
git clone https://github.com/hector-s-m/RFantibody_scaffolding.git
cd RFantibody_scaffolding

# 2. Install uv (if not already installed)
curl -LsSf https://astral.sh/uv/install.sh | sh

# 3. Download model weights (RFdiffusion)
bash include/download_weights.sh

# 4. Download AntiBMPNN weights (antibody-finetuned ProteinMPNN)
bash scripts/setup_antibmpnn.sh

# 5. Install RFantibody dependencies
uv sync

# 6. Set up Boltz2 environment (separate due to dependency conflicts)
conda create -n boltz_2.2.1 python=3.11 -y
conda activate boltz_2.2.1
pip install boltz[cuda]
conda deactivate

# 7. Verify
uv run rfdiffusion --help
conda run -n boltz_2.2.1 boltz --help
```

## Quick Start

### Full Pipeline (Recommended)

Run the entire pipeline (RFdiffusion → AntiBMPNN → Boltz2) in one command. The script handles conda environment switching automatically.

```bash
# Nanobody design
bash scripts/examples/motif_scaffolding_pipeline.sh \
    --Nb -m "inputs/target+motif.pdb"

# scFv design
bash scripts/examples/motif_scaffolding_pipeline.sh \
    --scFv -m "inputs/target+motif.pdb"
```

Common options:
```bash
bash scripts/examples/motif_scaffolding_pipeline.sh \
    --Nb \
    -m "inputs/target+motif.pdb" \
    -n 100 \                          # Number of backbone designs
    --design-loops "H1:,H2:,H3:10-16" \
    --motif-cdr H3 \                  # CDR loop for motif placement
    --hotspots "A100,A105" \          # Target hotspot residues
    --num-seqs 4 \                    # Sequences per backbone
    --temperature 0.2 \               # Sampling temperature
    --diffusion-samples 1 \           # Boltz2 diffusion samples
    -o output_Nb                      # Output directory
```

Outputs:
```
output_Nb/
├── designs/          # RFdiffusion backbone PDBs
├── mpnn_designs/     # AntiBMPNN sequence-designed PDBs
├── boltz2_input/     # Boltz2 YAML inputs
├── boltz2_output/    # Boltz2 predicted structures + confidence
└── boltz2_metrics.csv  # Ranked metrics (ipTM, ipSAE, pDockQ, ...)
```

### Standard Antibody Design (without motif)

```bash
uv run rfdiffusion \
    -t target.pdb \
    -f inputs/Nb.pdb \
    -o designs/ab \
    -n 100 \
    -l "H1:7,H2:6,H3:5-13" \
    -h "B146,B170,B177"
```

### Step-by-Step Motif Scaffolding

For more control, run each step individually:

```bash
# Step 1: RFdiffusion — design backbones with motif fixed in H3
uv run rfdiffusion \
    -f inputs/Nb.pdb \
    -m "inputs/target+motif.pdb" \
    -o designs/motif_ab \
    -n 100 \
    -l "H1:,H2:,H3:10-16" \
    --motif-cdr H3

# Step 2: AntiBMPNN — design sequences (motif residues stay fixed)
uv run proteinmpnn \
    -i designs/ \
    -o mpnn_designs/ \
    -n 4 -t 0.2

# Step 3: Boltz2 — predict structures + score (separate conda env)
python scripts/prepare_boltz2_input.py -i mpnn_designs/ -o boltz2_input/

conda activate boltz_2.2.1
boltz predict boltz2_input/ --out_dir boltz2_output/ \
    --use_msa_server --use_potentials --write_full_pae --write_full_pde
conda deactivate

python scripts/extract_boltz2_metrics.py \
    -i boltz2_output/ -d mpnn_designs/ -o boltz2_metrics.csv
```

## Input Files

All input files go in the `inputs/` directory.

| File | Description | Used by |
|------|-------------|---------|
| `inputs/Nb.pdb` | Nanobody framework (HLT format) | `--Nb` flag |
| `inputs/scFv.pdb` | scFv framework (HLT format) | `--scFv` flag |
| `inputs/target+motif.pdb` | Combined target + aromatic motif | `-m` flag |

### Framework PDB (`Nb.pdb`, `scFv.pdb`)

HLT-formatted antibody or nanobody scaffold:
- **Chain H**: Heavy chain
- **Chain L**: Light chain (scFv only)
- **Chain T**: Target protein
- Sequential 1-indexed residue numbering
- CDR annotations via `REMARK PDBinfo-LABEL` lines

Convert from Chothia format:
```bash
python scripts/util/chothia_to_HLT.py -inpdb mychothia.pdb -outpdb inputs/Nb.pdb
```

Example frameworks (for reference):
- Nanobody: `scripts/examples/example_inputs/h-NbBCII10.pdb`
- ScFv: `scripts/examples/example_inputs/hu-4D5-8_Fv.pdb`

### Combined Motif PDB (`target+motif.pdb`)

A single PDB file containing:
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

### Boltz2 (Structure Prediction + Scoring)

Runs in a separate conda environment (`boltz_2.2.1`). Three sub-steps:

```bash
# 3a. Convert PDBs → Boltz2 YAML (RFantibody env)
python scripts/prepare_boltz2_input.py -i mpnn_designs/ -o boltz2_input/

# 3b. Predict structures (boltz_2.2.1 env)
conda activate boltz_2.2.1
boltz predict boltz2_input/ --out_dir boltz2_output/ --use_msa_server --use_potentials --write_full_pae --write_full_pde

# 3c. Extract metrics (RFantibody env)
python scripts/extract_boltz2_metrics.py -i boltz2_output/ -d mpnn_designs/ -o metrics.csv
```

Or use the wrapper: `bash scripts/run_boltz2_predict.sh -i boltz2_input/ -o boltz2_output/`

**Metrics per design** (output CSV):

| Metric | Source | Description |
|--------|--------|-------------|
| `ipTM` | Boltz2 JSON | Interface predicted TM-score |
| `iPAE` | Boltz2 JSON | Interface Predicted Aligned Error |
| `pLDDT` | Boltz2 JSON/NPZ | Predicted Local Distance Difference Test |
| `iPLDDT` | Boltz2 JSON | Interface pLDDT |
| `ipSAE` | PAE NPZ | Interaction prediction Score from Aligned Errors (Dunbrack) |
| `pDockQ` | PAE + structure | Predicted DockQ score (Bryant et al.) |
| `pDockQ2` | PAE + structure | Improved pDockQ (Zhu et al.) |
| `LIS` | PAE NPZ | Local Interaction Score (Kim et al.) |
| `binder_RMSD` | Designed vs predicted | CA RMSD of binder chains (designed → Boltz2) |

### RF2 (Legacy)

RF2 is still available for lightweight prediction without Boltz2:

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
7. Boltz2 predicts the final complex structure and provides confidence metrics for ranking

## Filtering

Recommended criteria for design selection:
- ipTM > 0.7 (good interface prediction confidence)
- ipSAE > 0.5 (strong interface interaction score)
- pDockQ > 0.23 (acceptable docking quality)
- binder RMSD < 2 Å (designed structure matches prediction)
- pLDDT > 70 (overall structural confidence)

## License

MIT License. See LICENSE file.

Based on [RFantibody](https://www.biorxiv.org/content/10.1101/2024.03.14.585103v1) by Nate Bennett, Joe Watson, and the RFantibody Team.
