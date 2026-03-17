#!/bin/bash

# ============================================================================
# Aromatic Motif Scaffolding Pipeline
# ============================================================================
# Full end-to-end pipeline:
#   1. RFdiffusion  - Design antibody backbone with fixed motif in CDR loop
#   2. AntiBMPNN    - Design sequences (motif residues stay fixed)
#   3. Boltz2       - Predict structures + confidence metrics
#
# Conda environments:
#   - Steps 1 & 2 run in the "RFantibody" conda environment
#   - Step 3 runs in the "boltz_2.2.1" conda environment
#
# Usage:
#   bash scripts/examples/motif_scaffolding_pipeline.sh --scFv -m inputs/motif_combined.pdb
#   bash scripts/examples/motif_scaffolding_pipeline.sh --Nb   -m inputs/motif_combined.pdb
#
# Required:
#   --scFv  OR  --Nb        Select framework type
#   -m, --motif PATH        Combined motif+target PDB
#
# Optional:
#   -o, --output DIR        Output directory (default: output_scFv/ or output_Nb/)
#   -n, --num-designs N     Number of backbone designs (default: 100)
#   --design-loops STR      Loop length specification
#   --motif-cdr CDR         CDR loop for motif (default: H3)
#   --hotspots STR          Target hotspot residues
#   --num-seqs N            Sequences per backbone (default: 4)
#   --temperature FLOAT     Sampling temperature (default: 0.2)
#   --diffusion-samples N   Boltz2 diffusion samples (default: 1)
#   --boltz-cache PATH      Boltz2 cache directory
#   --diffuser-t N          Diffusion timesteps (default: 50)
# ============================================================================

set -e  # Exit on error

# Initialize conda for this shell (required for conda activate in scripts)
eval "$(conda shell.bash hook)"

# ============================================================================
# DEFAULTS
# ============================================================================

FRAMEWORK_TYPE=""               # Set by --scFv or --Nb
MOTIF_COMBINED_PDB=""           # Required: -m / --motif
OUTPUT_DIR=""                   # Auto-set based on framework if empty

# Conda environments
RFANTIBODY_ENV="RFantibody"
BOLTZ2_ENV="boltz_2.2.1"

# RFdiffusion parameters
NUM_DESIGNS=100
MOTIF_CDR="H3"
DESIGN_LOOPS=""                 # Auto-set based on framework if empty
HOTSPOTS=""
DIFFUSER_T=50

# AntiBMPNN parameters
NUM_SEQS=4
SAMPLING_TEMP=0.2

# Boltz2 parameters
DIFFUSION_SAMPLES=1
MSA_SERVER_URL="http://a3m-2023.mmseqs.com"
BOLTZ_CACHE=""

# ============================================================================
# PARSE ARGUMENTS
# ============================================================================

print_usage() {
    cat <<EOF
Usage: $(basename "$0") (--scFv | --Nb) -m MOTIF_PDB [OPTIONS]

Framework (required, pick one):
  --scFv                  Use scFv framework (inputs/scFv.pdb)
  --Nb                    Use nanobody framework (inputs/Nb.pdb)

Required:
  -m, --motif PATH        Combined motif+target PDB

Optional:
  -o, --output DIR        Output directory (default: output_scFv/ or output_Nb/)
  -n, --num-designs N     Number of backbone designs (default: $NUM_DESIGNS)
  --design-loops STR      Loop lengths, e.g. "H1:,H2:,H3:10-16,L1:,L2:,L3:"
  --motif-cdr CDR         CDR loop for motif (default: $MOTIF_CDR)
  --hotspots STR          Target hotspot residues (e.g. "A100,A105")
  --num-seqs N            Sequences per backbone (default: $NUM_SEQS)
  --temperature FLOAT     Sampling temperature (default: $SAMPLING_TEMP)
  --diffusion-samples N   Boltz2 diffusion samples (default: $DIFFUSION_SAMPLES)
  --boltz-cache PATH      Boltz2 cache directory
  --diffuser-t N          Diffusion timesteps (default: $DIFFUSER_T)
  -h, --help              Show this help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --scFv)              FRAMEWORK_TYPE="scFv"; shift ;;
        --Nb)                FRAMEWORK_TYPE="Nb"; shift ;;
        -m|--motif)          MOTIF_COMBINED_PDB="$2"; shift 2 ;;
        -o|--output)         OUTPUT_DIR="$2"; shift 2 ;;
        -n|--num-designs)    NUM_DESIGNS="$2"; shift 2 ;;
        --design-loops)      DESIGN_LOOPS="$2"; shift 2 ;;
        --motif-cdr)         MOTIF_CDR="$2"; shift 2 ;;
        --hotspots)          HOTSPOTS="$2"; shift 2 ;;
        --num-seqs)          NUM_SEQS="$2"; shift 2 ;;
        --temperature)       SAMPLING_TEMP="$2"; shift 2 ;;
        --diffusion-samples) DIFFUSION_SAMPLES="$2"; shift 2 ;;
        --boltz-cache)       BOLTZ_CACHE="$2"; shift 2 ;;
        --diffuser-t)        DIFFUSER_T="$2"; shift 2 ;;
        -h|--help)           print_usage; exit 0 ;;
        *)                   echo "Unknown option: $1"; print_usage; exit 1 ;;
    esac
done

# ============================================================================
# VALIDATE & CONFIGURE
# ============================================================================

if [ -z "$FRAMEWORK_TYPE" ]; then
    echo "Error: Must specify --scFv or --Nb"
    print_usage
    exit 1
fi

if [ -z "$MOTIF_COMBINED_PDB" ]; then
    echo "Error: Must specify --motif / -m with path to combined motif+target PDB"
    print_usage
    exit 1
fi

if [ ! -f "$MOTIF_COMBINED_PDB" ]; then
    echo "Error: Motif PDB not found: $MOTIF_COMBINED_PDB"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$PROJECT_ROOT"

# Framework-specific configuration
case "$FRAMEWORK_TYPE" in
    scFv)
        FRAMEWORK_PDB="inputs/scFv.pdb"
        LOOP_STRING="H1,H2,H3,L1,L2,L3"
        [ -z "$DESIGN_LOOPS" ] && DESIGN_LOOPS="H1:,H2:,H3:10-16,L1:,L2:,L3:"
        [ -z "$OUTPUT_DIR" ]   && OUTPUT_DIR="output_scFv"
        ;;
    Nb)
        FRAMEWORK_PDB="inputs/Nb.pdb"
        LOOP_STRING="H1,H2,H3"
        [ -z "$DESIGN_LOOPS" ] && DESIGN_LOOPS="H1:,H2:,H3:10-16"
        [ -z "$OUTPUT_DIR" ]   && OUTPUT_DIR="output_Nb"
        ;;
esac

if [ ! -f "$FRAMEWORK_PDB" ]; then
    echo "Error: Framework PDB not found: $FRAMEWORK_PDB"
    echo "Place your framework in inputs/"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Initialize conda
if command -v conda &>/dev/null; then
    CONDA_BASE=$(conda info --base 2>/dev/null)
    source "$CONDA_BASE/etc/profile.d/conda.sh"
else
    echo "Error: conda not found. Both environments require conda."
    exit 1
fi

echo "=============================================="
echo " Aromatic Motif Scaffolding Pipeline"
echo "=============================================="
echo "  Mode:        $FRAMEWORK_TYPE"
echo "  Framework:   $FRAMEWORK_PDB"
echo "  Motif:       $MOTIF_COMBINED_PDB -> CDR $MOTIF_CDR"
echo "  Designs:     $NUM_DESIGNS backbones x $NUM_SEQS seqs = $((NUM_DESIGNS * NUM_SEQS)) total"
echo "  Loops:       $DESIGN_LOOPS"
echo "  Output:      $OUTPUT_DIR"
echo "  Envs:        $RFANTIBODY_ENV (Steps 1-2) | $BOLTZ2_ENV (Step 3)"
echo "=============================================="

START_PIPELINE=$(date +%s)

# ============================================================================
# STEP 1: RFdiffusion with motif scaffolding (RFantibody env)
# ============================================================================

echo ""
echo "[Step 1/3] Running RFdiffusion with motif scaffolding..."
echo "  - Designing $NUM_DESIGNS backbones"
echo "  - Motif scaffolded into $MOTIF_CDR"
echo "  - Loop lengths: $DESIGN_LOOPS"

conda activate "$RFANTIBODY_ENV"

RFDIFF_CMD="rfdiffusion \
    --framework \"$FRAMEWORK_PDB\" \
    --output \"$OUTPUT_DIR/designs/ab_des\" \
    --num-designs $NUM_DESIGNS \
    --design-loops \"$DESIGN_LOOPS\" \
    --motif \"$MOTIF_COMBINED_PDB\" \
    --motif-cdr \"$MOTIF_CDR\" \
    --diffuser-t $DIFFUSER_T"

if [ -n "$HOTSPOTS" ]; then
    RFDIFF_CMD="$RFDIFF_CMD --hotspots \"$HOTSPOTS\""
fi

eval $RFDIFF_CMD

echo "[Step 1/3] RFdiffusion complete."

# ============================================================================
# STEP 2: AntiBMPNN with motif fixed positions (RFantibody env)
# ============================================================================

echo ""
echo "[Step 2/3] Running AntiBMPNN / ProteinMPNN..."
echo "  - Generating $NUM_SEQS sequences per backbone"
echo "  - Motif residues will remain fixed"
echo "  - Designing loops: $LOOP_STRING"

for design_pdb in "$OUTPUT_DIR"/designs/ab_des_*.pdb; do
    base=$(basename "$design_pdb" .pdb)
    motif_json="$OUTPUT_DIR/designs/${base}_motif_fixed.json"

    MPNN_ARGS="-pdbdir $OUTPUT_DIR/designs -outpdbdir $OUTPUT_DIR/mpnn_designs \
        -loop_string $LOOP_STRING \
        -seqs_per_struct $NUM_SEQS \
        -temperature $SAMPLING_TEMP"

    if [ -f "$motif_json" ]; then
        MPNN_ARGS="$MPNN_ARGS -motif_fixed_positions $motif_json"
    fi

    python scripts/proteinmpnn_interface_design.py $MPNN_ARGS
done

echo "[Step 2/3] AntiBMPNN complete."

conda deactivate

# ============================================================================
# STEP 3: Boltz2 — Structure prediction + scoring (boltz_2.2.1 env)
# ============================================================================

echo ""
echo "[Step 3/3] Running Boltz2 structure prediction + scoring..."
echo "  - Converting PDBs to Boltz2 YAML format"
echo "  - Predicting with $DIFFUSION_SAMPLES diffusion sample(s)"

# 3a. Prepare Boltz2 YAML inputs (RFantibody env — only needs pyyaml)
conda activate "$RFANTIBODY_ENV"

python scripts/prepare_boltz2_input.py \
    -i "$OUTPUT_DIR/mpnn_designs" \
    -o "$OUTPUT_DIR/boltz2_input"

conda deactivate

# 3b. Run Boltz2 prediction (boltz_2.2.1 env)
conda activate "$BOLTZ2_ENV"

BOLTZ_CMD="boltz predict \"$OUTPUT_DIR/boltz2_input\" \
    --out_dir \"$OUTPUT_DIR/boltz2_output\" \
    --diffusion_samples $DIFFUSION_SAMPLES \
    --use_msa_server \
    --msa_server_url=$MSA_SERVER_URL \
    --use_potentials \
    --write_full_pae \
    --write_full_pde"

if [ -n "$BOLTZ_CACHE" ]; then
    BOLTZ_CMD="$BOLTZ_CMD --cache \"$BOLTZ_CACHE\""
fi

echo "Running: $BOLTZ_CMD"
eval $BOLTZ_CMD

conda deactivate

# 3c. Extract metrics (RFantibody env — needs numpy + pandas)
conda activate "$RFANTIBODY_ENV"

python scripts/extract_boltz2_metrics.py \
    -i "$OUTPUT_DIR/boltz2_output" \
    -d "$OUTPUT_DIR/mpnn_designs" \
    -o "$OUTPUT_DIR/boltz2_metrics.csv" \
    --rank-by ipTM

conda deactivate

echo "[Step 3/3] Boltz2 complete."

# ============================================================================
# SUMMARY
# ============================================================================

END_PIPELINE=$(date +%s)
ELAPSED=$(( (END_PIPELINE - START_PIPELINE) / 60 ))

echo ""
echo "=============================================="
echo " Pipeline Complete! ($FRAMEWORK_TYPE, ${ELAPSED} min)"
echo "=============================================="
echo "Outputs:"
echo "  1. RFdiffusion backbones: $OUTPUT_DIR/designs/"
echo "  2. AntiBMPNN sequences:   $OUTPUT_DIR/mpnn_designs/"
echo "  3. Boltz2 structures:     $OUTPUT_DIR/boltz2_output/"
echo "  4. Ranked metrics CSV:    $OUTPUT_DIR/boltz2_metrics.csv"
echo ""
echo "Metrics: ipTM, ipSAE, pDockQ, pDockQ2, LIS, pLDDT, iPLDDT, iPAE, binder_RMSD"
echo ""
echo "Review boltz2_metrics.csv to identify top candidates."
echo "=============================================="
