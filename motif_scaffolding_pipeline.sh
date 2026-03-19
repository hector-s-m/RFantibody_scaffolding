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
# Naming convention (example: target_name=PROTEIN, framework_type=Nb):
#   output/PROTEIN_Nb/
#     RFdiffusion_backbones/  PROTEIN_Nb_RF0.pdb
#     AntiBMPNN_sequences/    PROTEIN_Nb_RF0_mpnn0.pdb
#     Boltz-2_predictions/    PROTEIN_Nb_RF0_mpnn0_0.pdb  (converted from .cif)
#     designs_registry.csv
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
TARGET_NAME=$(json_get "target_name")
FRAMEWORK_TYPE=$(json_get "framework_type")
FRAMEWORK_PDB=$(json_get "framework_pdb")
MOTIF_COMBINED_PDB=$(json_get "motif_pdb")

RFANTIBODY_ENV=$(json_get "conda_env_rfantibody")
BOLTZ2_ENV=$(json_get "conda_env_boltz2")

NUM_DESIGNS=$(json_get "rfdiffusion" "num_designs")
MOTIF_CDR=$(json_get "rfdiffusion" "motif_cdr")
DESIGN_LOOPS=$(json_get "rfdiffusion" "design_loops")
HOTSPOTS=$(json_get "rfdiffusion" "hotspots")
DIFFUSER_T=$(json_get "rfdiffusion" "diffuser_T")

NUM_SEQS=$(json_get "antibmpnn" "num_seqs")
SAMPLING_TEMP=$(json_get "antibmpnn" "sampling_temp")

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
  --target-name NAME        Target protein name (used in output naming)
  -m, --motif PATH          Combined motif+target PDB
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
        --target-name)       TARGET_NAME="$2"; shift 2 ;;
        -m|--motif)          MOTIF_COMBINED_PDB="$2"; shift 2 ;;
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

if [ -z "$TARGET_NAME" ]; then
    echo "Error: target_name not set in JSON or via --target-name"
    exit 1
fi

if [ -z "$TARGET_NAME" ]; then
    echo "Error: target_name not set in JSON or via --target-name"
    exit 1
fi

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

# Validate numeric parameters early (before expensive steps)
if [ "$NUM_DESIGNS" -lt 1 ] 2>/dev/null; then
    echo "Error: num_designs must be >= 1 (got: $NUM_DESIGNS)"
    exit 1
fi
if [ "$NUM_SEQS" -lt 1 ] 2>/dev/null; then
    echo "Error: num_seqs must be >= 1 (got: $NUM_SEQS)"
    exit 1
fi

# Naming prefix: TARGET_TYPE (e.g., PROTEIN_Nb)
PREFIX="${TARGET_NAME}_${FRAMEWORK_TYPE}"

# Output directory: output/PREFIX
OUTPUT_DIR="output/${PREFIX}"

# Framework-specific LOOP_STRING for MPNN (which loops to redesign)
case "$FRAMEWORK_TYPE" in
    scFv)
        LOOP_STRING="H1,H2,H3,L1,L2,L3"
        [ -z "$DESIGN_LOOPS" ] && DESIGN_LOOPS="H1:,H2:,H3:10-16,L1:,L2:,L3:"
        ;;
    Nb)
        LOOP_STRING="H1,H2,H3"
        [ -z "$DESIGN_LOOPS" ] && DESIGN_LOOPS="H1:,H2:,H3:10-16"
        ;;
    *)
        echo "Error: Unknown framework_type '$FRAMEWORK_TYPE'. Use 'scFv' or 'Nb'."
        exit 1
        ;;
esac

if [ ! -f "$FRAMEWORK_PDB" ]; then
    echo "Error: Framework PDB not found: $FRAMEWORK_PDB"
    exit 1
fi

# Create output subdirectories
RFDIFF_DIR="$OUTPUT_DIR/RFdiffusion_backbones"
MPNN_DIR="$OUTPUT_DIR/AntiBMPNN_sequences"
BOLTZ_PRED_DIR="$OUTPUT_DIR/Boltz-2_predictions"
REGISTRY_CSV="$OUTPUT_DIR/designs_registry.csv"

mkdir -p "$RFDIFF_DIR" "$MPNN_DIR" "$BOLTZ_PRED_DIR"

# Internal working directories (not user-facing) — clean at start for idempotent re-runs
MPNN_RAW_DIR="$OUTPUT_DIR/_mpnn_raw"
BOLTZ_YAML_DIR="$OUTPUT_DIR/_boltz_yaml"
BOLTZ_RAW_DIR="$OUTPUT_DIR/_boltz_raw"
rm -rf "$MPNN_RAW_DIR" "$BOLTZ_YAML_DIR" "$BOLTZ_RAW_DIR"
mkdir -p "$MPNN_RAW_DIR" "$BOLTZ_YAML_DIR" "$BOLTZ_RAW_DIR"

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
echo "  Target:      $TARGET_NAME"
echo "  Mode:        $FRAMEWORK_TYPE"
echo "  Prefix:      $PREFIX"
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
ABS_OUTPUT_PREFIX=$(mkdir -p "$RFDIFF_DIR" && realpath "$RFDIFF_DIR")/${PREFIX}_RF
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

# Rename RFdiffusion outputs: PREFIX_RF_N → PREFIX_RFN (no underscore before number)
# Uses bash parameter expansion (not sed) to avoid regex metacharacter issues in PREFIX
OLD_PATTERN="${PREFIX}_RF_"
NEW_PATTERN="${PREFIX}_RF"
for f in "$RFDIFF_DIR"/${PREFIX}_RF_*.pdb "$RFDIFF_DIR"/${PREFIX}_RF_*.trb "$RFDIFF_DIR"/${PREFIX}_RF_*_motif_fixed.json; do
    [ -f "$f" ] || continue
    dir=$(dirname "$f")
    base=$(basename "$f")
    new_base="${NEW_PATTERN}${base#${OLD_PATTERN}}"
    [ "$base" != "$new_base" ] && mv "$f" "$dir/$new_base"
done

echo "[Step 1/3] RFdiffusion complete."
echo "  Output: $RFDIFF_DIR/${PREFIX}_RF*.pdb"

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
for design_pdb in "$RFDIFF_DIR"/${PREFIX}_RF*.pdb; do
    DESIGN_COUNT=$((DESIGN_COUNT + 1))
done
echo "  - Processing $DESIGN_COUNT backbone designs"

MPNN_IDX=0
for design_pdb in "$RFDIFF_DIR"/${PREFIX}_RF*.pdb; do
    MPNN_IDX=$((MPNN_IDX + 1))
    base=$(basename "$design_pdb" .pdb)
    motif_json="$RFDIFF_DIR/${base}_motif_fixed.json"

    # Create a runlist with just this one design
    RUNLIST_FILE="$RFDIFF_DIR/_runlist_tmp.txt"
    echo "$base" > "$RUNLIST_FILE"

    MPNN_ARGS="-pdbdir $RFDIFF_DIR -outpdbdir $MPNN_RAW_DIR \
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
rm -f "$RFDIFF_DIR/_runlist_tmp.txt"

# --- Rename MPNN outputs to final naming convention ---
# MPNN generates: {PREFIX}_RF{N}_dldesign_{M}.pdb
# Rename to:      {PREFIX}_RF{N}_mpnn{M}.pdb
CURRENT_STEP="Step 2/3 — Renaming MPNN outputs"
echo "  Renaming MPNN outputs..."

for f in "$MPNN_RAW_DIR"/${PREFIX}_RF*_dldesign_*.pdb; do
    [ -f "$f" ] || continue
    fname=$(basename "$f" .pdb)
    # Extract RFN and M from {PREFIX}_RF{N}_dldesign_{M}
    # Remove PREFIX_ prefix to get RF{N}_dldesign_{M}
    suffix="${fname#${PREFIX}_}"
    # suffix is now e.g. RF0_dldesign_2
    RF_PART="${suffix%%_dldesign_*}"    # RF0
    M="${suffix##*_dldesign_}"          # 2
    new_name="${PREFIX}_${RF_PART}_mpnn${M}.pdb"
    cp "$f" "$MPNN_DIR/$new_name"
done

MPNN_COUNT=$(find "$MPNN_DIR" -name "*.pdb" 2>/dev/null | wc -l)
echo "[Step 2/3] AntiBMPNN complete. Generated $MPNN_COUNT sequence designs."
echo "  Output: $MPNN_DIR/${PREFIX}_RF*_mpnn*.pdb"

conda deactivate

# ============================================================================
# STEP 3: Boltz2 — Structure prediction + scoring (boltz_2.2.1 env)
# ============================================================================

CURRENT_STEP="Step 3a/3 — Prepare Boltz2 YAML inputs"

echo ""
echo "[Step 3/3] Running Boltz2 structure prediction + scoring..."
echo "  - Converting PDBs to Boltz2 YAML format"
echo "  - Predicting with $DIFFUSION_SAMPLES diffusion sample(s) per design (batch=$BOLTZ_BATCH_SIZE parallel)"

# 3a. Prepare Boltz2 YAML inputs from the renamed MPNN outputs
conda deactivate 2>/dev/null || true
conda activate "$RFANTIBODY_ENV"

python scripts/prepare_boltz2_input.py \
    -i "$MPNN_DIR" \
    -o "$BOLTZ_YAML_DIR" \
    --remap-chains \
    --target-chains T

YAML_COUNT=$(find "$BOLTZ_YAML_DIR" -name "*.yaml" 2>/dev/null | wc -l)
echo "  Generated $YAML_COUNT Boltz2 YAML input files"

conda deactivate 2>/dev/null || true

# 3b. Run Boltz2 prediction (boltz_2.2.1 env) — mmcif output (Boltz2 native)
CURRENT_STEP="Step 3b/3 — Boltz2 structure prediction"

conda activate "$BOLTZ2_ENV"
echo "  Activated environment: $BOLTZ2_ENV"
echo "  boltz location: $(which boltz 2>/dev/null || echo 'NOT FOUND')"
echo "  Parallel batch size: $BOLTZ_BATCH_SIZE"

BOLTZ_PARALLEL_CMD=(python scripts/run_boltz2_parallel.py
    -i "$BOLTZ_YAML_DIR"
    -o "$BOLTZ_RAW_DIR"
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

# 3c. Convert CIF → PDB and rename to final convention
CURRENT_STEP="Step 3c/3 — Convert CIF to PDB + rename"

conda activate "$RFANTIBODY_ENV"

echo "  Converting CIF to PDB and organizing outputs..."

# Find all Boltz2 prediction CIF files and convert+rename
# Boltz2 creates: boltz_results_*/predictions/{stem}/{stem}_model_{S}.cif
# Rename to: {PREFIX}_{N}_mpnn{M}_{S}.pdb
for cif_file in $(find "$BOLTZ_RAW_DIR" -name "*.cif" -type f | sort); do
    cif_name=$(basename "$cif_file" .cif)
    # Extract model number: {stem}_model_{S} → S
    if [[ "$cif_name" =~ _model_([0-9]+)$ ]]; then
        model_num="${BASH_REMATCH[1]}"
        stem="${cif_name%_model_*}"
        # stem is e.g. PROTEIN_Nb_0_mpnn0 (already the right format from YAML naming)
        pdb_name="${stem}_${model_num}.pdb"
        python scripts/convert_cif_to_pdb.py "$cif_file" "$BOLTZ_PRED_DIR/$pdb_name"
    fi
done

PRED_COUNT=$(find "$BOLTZ_PRED_DIR" -name "*.pdb" 2>/dev/null | wc -l)
echo "  Converted $PRED_COUNT structures to PDB format"

# 3d. Extract metrics (reads from _boltz_raw for JSON/NPZ, MPNN_DIR for designed PDBs)
CURRENT_STEP="Step 3d/3 — Extract Boltz2 metrics"

python scripts/extract_boltz2_metrics.py \
    -i "$BOLTZ_RAW_DIR" \
    -d "$MPNN_DIR" \
    -o "$REGISTRY_CSV" \
    --rank-by ipTM

conda deactivate 2>/dev/null || true

echo "[Step 3/3] Boltz2 complete."

# ============================================================================
# CLEANUP: Remove internal working directories
# ============================================================================

rm -rf "$MPNN_RAW_DIR" "$BOLTZ_YAML_DIR"
# Keep BOLTZ_RAW_DIR for now (contains NPZ/JSON needed for re-analysis)

# ============================================================================
# SUMMARY
# ============================================================================

END_PIPELINE=$(date +%s)
ELAPSED=$(( (END_PIPELINE - START_PIPELINE) / 60 ))

echo ""
echo "=============================================="
echo " Pipeline Complete! ($PREFIX, ${ELAPSED} min)"
echo "=============================================="
echo "Outputs:"
echo "  $OUTPUT_DIR/"
echo "    RFdiffusion_backbones/  ${PREFIX}_RF*.pdb"
echo "    AntiBMPNN_sequences/    ${PREFIX}_RF*_mpnn*.pdb"
echo "    Boltz-2_predictions/    ${PREFIX}_RF*_mpnn*_*.pdb"
echo "    designs_registry.csv    (ranked by ipTM)"
echo "    pipeline_parameters.json"
echo ""
echo "Metrics: ipTM, ipSAE, pDockQ, pDockQ2, LIS, pLDDT, iPLDDT, iPAE, binder_RMSD, motif_RMSD"
echo ""
echo "Review designs_registry.csv to identify top candidates."
echo "=============================================="
