#!/bin/bash

# ============================================================================
# Aromatic Motif Scaffolding Pipeline
# ============================================================================
# Full end-to-end pipeline:
#   1. RFdiffusion  - Design antibody backbone with fixed motif in CDR loop
#   2. AntiBMPNN    - Design sequences (motif residues stay fixed)
#   3. Boltz2       - Predict structures + confidence metrics
#
# Configuration:
#   Edit pipeline_parameters.json to set all parameters, or override via CLI.
#
# Usage:
#   bash motif_scaffolding_pipeline.sh                        # uses pipeline_parameters.json
#   bash motif_scaffolding_pipeline.sh --config my_run.json   # custom config
#   bash motif_scaffolding_pipeline.sh --Nb --diffuser-t 50   # CLI overrides
# ============================================================================

set -e  # Exit on error

# Track current step for error diagnostics
CURRENT_STEP="init"
trap 'echo ""; echo "ERROR: Pipeline failed during [$CURRENT_STEP] (line $LINENO)" >&2; exit 1' ERR

# Initialize conda for this shell (required for conda activate in scripts)
eval "$(conda shell.bash hook)"

# ============================================================================
# LOCATE PROJECT ROOT & CONFIG
# ============================================================================

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_ROOT"

CONFIG_FILE="$PROJECT_ROOT/pipeline_parameters.json"

# ============================================================================
# HELPER: Read a value from JSON config using python (no jq dependency)
# ============================================================================

json_get() {
    # Usage: json_get KEY [NESTED_KEY]
    # Examples:
    #   json_get "framework_type"         -> "Nb"
    #   json_get "rfdiffusion" "num_designs"  -> 100
    local file="$CONFIG_FILE"
    if [ $# -eq 1 ]; then
        python3 -c "import json; d=json.load(open('$file')); print(d.get('$1', ''))"
    elif [ $# -eq 2 ]; then
        python3 -c "import json; d=json.load(open('$file')); print(d.get('$1', {}).get('$2', ''))"
    fi
}

# ============================================================================
# DEFAULTS (from JSON config, then CLI overrides)
# ============================================================================

# Parse --config first if provided (before loading defaults)
for arg in "$@"; do
    if [[ "$arg" == "--config" ]]; then
        shift
        CONFIG_FILE="$1"
        shift
        break
    fi
done

if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Config file not found: $CONFIG_FILE"
    echo "Create pipeline_parameters.json or specify --config PATH"
    exit 1
fi

echo "Loading config from: $CONFIG_FILE"

# Load from JSON
FRAMEWORK_TYPE=$(json_get "framework_type")
FRAMEWORK_PDB=$(json_get "framework_pdb")
MOTIF_COMBINED_PDB=$(json_get "motif_pdb")
OUTPUT_DIR=$(json_get "output_dir")

RFANTIBODY_ENV=$(json_get "conda_env_rfantibody")
BOLTZ2_ENV=$(json_get "conda_env_boltz2")

NUM_DESIGNS=$(json_get "rfdiffusion" "num_designs")
MOTIF_CDR=$(json_get "rfdiffusion" "motif_cdr")
DESIGN_LOOPS=$(json_get "rfdiffusion" "design_loops")
HOTSPOTS=$(json_get "rfdiffusion" "hotspots")
DIFFUSER_T=$(json_get "rfdiffusion" "diffuser_T")

NUM_SEQS=$(json_get "mpnn" "num_seqs")
SAMPLING_TEMP=$(json_get "mpnn" "sampling_temp")

DIFFUSION_SAMPLES=$(json_get "boltz2" "diffusion_samples")
BOLTZ_BATCH_SIZE=$(json_get "boltz2" "batch_size")
MSA_SERVER_URL=$(json_get "boltz2" "msa_server_url")
BOLTZ_CACHE=$(json_get "boltz2" "cache")

# ============================================================================
# CLI OVERRIDES (take precedence over JSON)
# ============================================================================

print_usage() {
    cat <<EOF
Usage: $(basename "$0") [--config CONFIG_JSON] [OPTIONS]

Config:
  --config PATH             JSON config file (default: pipeline_parameters.json)

Framework override (sets framework_type + framework_pdb):
  --scFv                    Use scFv framework (inputs/scFv.pdb)
  --Nb                      Use nanobody framework (inputs/Nb.pdb)

Parameter overrides (override JSON values):
  -m, --motif PATH          Combined motif+target PDB
  -o, --output DIR          Output directory
  -n, --num-designs N       Number of backbone designs
  --design-loops STR        Loop lengths, e.g. "H1:,H2:,H3:10-16"
  --motif-cdr CDR           CDR loop for motif
  --hotspots STR            Target hotspot residues (e.g. "T5,T12")
  --num-seqs N              Sequences per backbone
  --temperature FLOAT       Sampling temperature
  --diffuser-t N            Diffusion timesteps (200=best quality, 50=fast)
  --diffusion-samples N     Boltz2 diffusion samples per design
  --boltz-batch-size N      Boltz2 parallel predictions per batch
  --boltz-cache PATH        Boltz2 cache directory
  -h, --help                Show this help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)            CONFIG_FILE="$2"; shift 2 ;;  # Already handled above, skip
        --scFv)              FRAMEWORK_TYPE="scFv"; FRAMEWORK_PDB="inputs/scFv.pdb"; shift ;;
        --Nb)                FRAMEWORK_TYPE="Nb"; FRAMEWORK_PDB="inputs/Nb.pdb"; shift ;;
        -m|--motif)          MOTIF_COMBINED_PDB="$2"; shift 2 ;;
        -o|--output)         OUTPUT_DIR="$2"; shift 2 ;;
        -n|--num-designs)    NUM_DESIGNS="$2"; shift 2 ;;
        --design-loops)      DESIGN_LOOPS="$2"; shift 2 ;;
        --motif-cdr)         MOTIF_CDR="$2"; shift 2 ;;
        --hotspots)          HOTSPOTS="$2"; shift 2 ;;
        --num-seqs)          NUM_SEQS="$2"; shift 2 ;;
        --temperature)       SAMPLING_TEMP="$2"; shift 2 ;;
        --diffusion-samples) DIFFUSION_SAMPLES="$2"; shift 2 ;;
        --boltz-batch-size)  BOLTZ_BATCH_SIZE="$2"; shift 2 ;;
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
    echo "Error: framework_type not set in JSON or via --scFv/--Nb"
    exit 1
fi

if [ -z "$MOTIF_COMBINED_PDB" ]; then
    echo "Error: motif_pdb not set in JSON or via --motif / -m"
    exit 1
fi

if [ ! -f "$MOTIF_COMBINED_PDB" ]; then
    echo "Error: Motif PDB not found: $MOTIF_COMBINED_PDB"
    exit 1
fi

# Framework-specific LOOP_STRING for MPNN (which loops to redesign)
case "$FRAMEWORK_TYPE" in
    scFv)
        LOOP_STRING="H1,H2,H3,L1,L2,L3"
        [ -z "$DESIGN_LOOPS" ] && DESIGN_LOOPS="H1:,H2:,H3:10-16,L1:,L2:,L3:"
        [ -z "$OUTPUT_DIR" ]   && OUTPUT_DIR="output_scFv"
        ;;
    Nb)
        LOOP_STRING="H1,H2,H3"
        [ -z "$DESIGN_LOOPS" ] && DESIGN_LOOPS="H1:,H2:,H3:10-16"
        [ -z "$OUTPUT_DIR" ]   && OUTPUT_DIR="output_Nb"
        ;;
    *)
        echo "Error: Unknown framework_type '$FRAMEWORK_TYPE'. Use 'scFv' or 'Nb'."
        exit 1
        ;;
esac

if [ ! -f "$FRAMEWORK_PDB" ]; then
    echo "Error: Framework PDB not found: $FRAMEWORK_PDB"
    echo "Set framework_pdb in pipeline_parameters.json or place in inputs/"
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

# Copy config to output directory for reproducibility
cp "$CONFIG_FILE" "$OUTPUT_DIR/pipeline_parameters.json"

echo "=============================================="
echo " Aromatic Motif Scaffolding Pipeline"
echo "=============================================="
echo "  Config:      $CONFIG_FILE"
echo "  Mode:        $FRAMEWORK_TYPE"
echo "  Framework:   $FRAMEWORK_PDB"
echo "  Motif:       $MOTIF_COMBINED_PDB -> CDR $MOTIF_CDR"
echo "  Designs:     $NUM_DESIGNS backbones x $NUM_SEQS seqs = $((NUM_DESIGNS * NUM_SEQS)) total"
echo "  Loops:       $DESIGN_LOOPS"
echo "  Diffuser T:  $DIFFUSER_T"
echo "  Output:      $OUTPUT_DIR"
echo "  Envs:        $RFANTIBODY_ENV (Steps 1-2) | $BOLTZ2_ENV (Step 3)"
echo "=============================================="

START_PIPELINE=$(date +%s)

# ============================================================================
# STEP 1: RFdiffusion with motif scaffolding (RFantibody env)
# ============================================================================

CURRENT_STEP="Step 1/3 — RFdiffusion backbone design"

echo ""
echo "[Step 1/3] Running RFdiffusion with motif scaffolding..."
echo "  - Designing $NUM_DESIGNS backbones"
echo "  - Motif scaffolded into $MOTIF_CDR"
echo "  - Loop lengths: $DESIGN_LOOPS"

conda activate "$RFANTIBODY_ENV"

# Resolve paths to absolute (Hydra changes cwd)
ABS_FRAMEWORK=$(realpath "$FRAMEWORK_PDB")
ABS_OUTPUT_PREFIX=$(mkdir -p "$OUTPUT_DIR/designs" && realpath "$OUTPUT_DIR/designs")/ab_des
ABS_MOTIF=$(realpath "$MOTIF_COMBINED_PDB")

# Parse design loops into Hydra format
IFS=',' read -ra LOOP_ITEMS <<< "$DESIGN_LOOPS"
LOOP_STR=$(printf ",%s" "${LOOP_ITEMS[@]}")
LOOP_STR="${LOOP_STR:1}"  # Remove leading comma

# Build RFdiffusion command (calling script directly — no entry point dependency)
RFDIFF_CMD=(python scripts/rfdiffusion_inference.py --config-name antibody
    "antibody.framework_pdb=$ABS_FRAMEWORK"
    "inference.output_prefix=$ABS_OUTPUT_PREFIX"
    "inference.num_designs=$NUM_DESIGNS"
    "diffuser.T=$DIFFUSER_T"
    "antibody.motif_pdb=$ABS_MOTIF"
    "antibody.motif_cdr_loop=$MOTIF_CDR"
    "antibody.design_loops=[$LOOP_STR]"
)

if [ -n "$HOTSPOTS" ]; then
    IFS=',' read -ra HS_ITEMS <<< "$HOTSPOTS"
    HS_STR=$(printf ",%s" "${HS_ITEMS[@]}")
    HS_STR="${HS_STR:1}"
    RFDIFF_CMD+=("ppi.hotspot_res=[$HS_STR]")
fi

# Auto-detect weights
if [ -f "weights/RFdiffusion_Ab.pt" ]; then
    RFDIFF_CMD+=("inference.ckpt_override_path=$(realpath weights/RFdiffusion_Ab.pt)")
fi

"${RFDIFF_CMD[@]}"

echo "[Step 1/3] RFdiffusion complete."

# ============================================================================
# STEP 2: AntiBMPNN with motif fixed positions (RFantibody env)
# ============================================================================

CURRENT_STEP="Step 2/3 — AntiBMPNN sequence design"

echo ""
echo "[Step 2/3] Running AntiBMPNN / ProteinMPNN..."
echo "  - Generating $NUM_SEQS sequences per backbone"
echo "  - Motif residues will remain fixed"
echo "  - Designing loops: $LOOP_STRING"

DESIGN_COUNT=0
for design_pdb in "$OUTPUT_DIR"/designs/ab_des_*.pdb; do
    DESIGN_COUNT=$((DESIGN_COUNT + 1))
done
echo "  - Processing $DESIGN_COUNT backbone designs"

MPNN_IDX=0
for design_pdb in "$OUTPUT_DIR"/designs/ab_des_*.pdb; do
    MPNN_IDX=$((MPNN_IDX + 1))
    base=$(basename "$design_pdb" .pdb)
    motif_json="$OUTPUT_DIR/designs/${base}_motif_fixed.json"

    # Create a runlist with just this one design (so -pdbdir processes only it)
    RUNLIST_FILE="$OUTPUT_DIR/designs/_runlist_tmp.txt"
    echo "$base" > "$RUNLIST_FILE"

    MPNN_ARGS="-pdbdir $OUTPUT_DIR/designs -outpdbdir $OUTPUT_DIR/mpnn_designs \
        -runlist $RUNLIST_FILE \
        -loop_string $LOOP_STRING \
        -seqs_per_struct $NUM_SEQS \
        -temperature $SAMPLING_TEMP"

    if [ -f "$motif_json" ]; then
        MPNN_ARGS="$MPNN_ARGS -motif_fixed_positions $motif_json"
    fi

    echo "  [$MPNN_IDX/$DESIGN_COUNT] $base"
    python scripts/proteinmpnn_interface_design.py $MPNN_ARGS
done
rm -f "$OUTPUT_DIR/designs/_runlist_tmp.txt"

MPNN_COUNT=$(find "$OUTPUT_DIR/mpnn_designs" -name "*.pdb" 2>/dev/null | wc -l)
echo "[Step 2/3] AntiBMPNN complete. Generated $MPNN_COUNT sequence designs."

conda deactivate

# ============================================================================
# STEP 3: Boltz2 — Structure prediction + scoring (boltz_2.2.1 env)
# ============================================================================

CURRENT_STEP="Step 3a/3 — Prepare Boltz2 YAML inputs"

echo ""
echo "[Step 3/3] Running Boltz2 structure prediction + scoring..."
echo "  - Converting PDBs to Boltz2 YAML format"
echo "  - Predicting with $DIFFUSION_SAMPLES diffusion sample(s) per design (batch=$BOLTZ_BATCH_SIZE parallel)"

# 3a. Prepare Boltz2 YAML inputs (RFantibody env — only needs pyyaml)
conda deactivate 2>/dev/null || true
conda activate "$RFANTIBODY_ENV"

python scripts/prepare_boltz2_input.py \
    -i "$OUTPUT_DIR/mpnn_designs" \
    -o "$OUTPUT_DIR/boltz2_input" \
    --remap-chains

YAML_COUNT=$(find "$OUTPUT_DIR/boltz2_input" -name "*.yaml" 2>/dev/null | wc -l)
echo "  Generated $YAML_COUNT Boltz2 YAML input files"

conda deactivate 2>/dev/null || true

# 3b. Run Boltz2 prediction (boltz_2.2.1 env)
CURRENT_STEP="Step 3b/3 — Boltz2 structure prediction"

conda activate "$BOLTZ2_ENV"
echo "  Activated environment: $BOLTZ2_ENV"
echo "  boltz location: $(which boltz 2>/dev/null || echo 'NOT FOUND — is boltz installed in $BOLTZ2_ENV?')"
echo "  Parallel batch size: $BOLTZ_BATCH_SIZE"

BOLTZ_PARALLEL_CMD=(python scripts/run_boltz2_parallel.py
    -i "$OUTPUT_DIR/boltz2_input"
    -o "$OUTPUT_DIR/boltz2_output"
    --samples "$DIFFUSION_SAMPLES"
    --batch-size "$BOLTZ_BATCH_SIZE"
    --msa-server-url "$MSA_SERVER_URL"
)

if [ -n "$BOLTZ_CACHE" ]; then
    BOLTZ_PARALLEL_CMD+=(--cache "$BOLTZ_CACHE")
fi

BOLTZ_START=$(date +%s)

"${BOLTZ_PARALLEL_CMD[@]}"

BOLTZ_END=$(date +%s)
BOLTZ_ELAPSED=$(( (BOLTZ_END - BOLTZ_START) / 60 ))
echo "  Boltz2 prediction complete (${BOLTZ_ELAPSED} min)"

conda deactivate 2>/dev/null || true

# 3c. Extract metrics (RFantibody env — needs numpy + pandas)
CURRENT_STEP="Step 3c/3 — Extract Boltz2 metrics"

conda activate "$RFANTIBODY_ENV"

python scripts/extract_boltz2_metrics.py \
    -i "$OUTPUT_DIR/boltz2_output" \
    -d "$OUTPUT_DIR/mpnn_designs" \
    -o "$OUTPUT_DIR/boltz2_metrics.csv" \
    --rank-by ipTM

conda deactivate 2>/dev/null || true

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
echo "  5. Config snapshot:       $OUTPUT_DIR/pipeline_parameters.json"
echo ""
echo "Metrics: ipTM, ipSAE, pDockQ, pDockQ2, LIS, pLDDT, iPLDDT, iPAE, binder_RMSD, motif_RMSD"
echo ""
echo "Review boltz2_metrics.csv to identify top candidates."
echo "=============================================="
