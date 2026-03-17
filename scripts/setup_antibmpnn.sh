#!/usr/bin/env bash
# =============================================================================
# Download AntiBMPNN weights for antibody-specific CDR sequence design.
#
# AntiBMPNN is a ProteinMPNN model finetuned on antibody 3D structures,
# achieving >80% sequence recovery and sub-nM experimental binders.
# It is a drop-in replacement for vanilla ProteinMPNN weights.
#
# Source: https://github.com/zeysun/AntiBMPNN
# License: MIT
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
WEIGHTS_DIR="${PROJECT_ROOT}/weights"

ZENODO_URL="https://zenodo.org/records/13387792/files/model_weights.zip"
ZIP_FILE="${WEIGHTS_DIR}/antibmpnn_weights.zip"
ANTIBMPNN_DIR="${WEIGHTS_DIR}/antibmpnn"

# Target weight file (the main model used by default)
TARGET_WEIGHT="${ANTIBMPNN_DIR}/AntiBMPNN_v48_noise_0.2.pt"

echo "============================================"
echo " AntiBMPNN Weight Setup"
echo "============================================"
echo ""

# Check if weights already exist
if [ -f "$TARGET_WEIGHT" ]; then
    echo "AntiBMPNN weights already exist at:"
    echo "  $TARGET_WEIGHT"
    echo ""
    echo "To re-download, delete the file and run this script again."
    exit 0
fi

# Create weights directory
mkdir -p "$ANTIBMPNN_DIR"

# Download
echo "Downloading AntiBMPNN weights from Zenodo..."
echo "  URL: $ZENODO_URL"
echo ""

if command -v wget &>/dev/null; then
    wget -O "$ZIP_FILE" "$ZENODO_URL"
elif command -v curl &>/dev/null; then
    curl -L -o "$ZIP_FILE" "$ZENODO_URL"
else
    echo "Error: Neither wget nor curl found. Please install one and try again."
    exit 1
fi

# Extract
echo ""
echo "Extracting weights..."
unzip -o "$ZIP_FILE" -d "$ANTIBMPNN_DIR"

# Find the actual .pt file(s) and normalize naming
# AntiBMPNN ships weights with varying names — find and symlink to our expected name
echo ""
echo "Locating model checkpoint..."

PT_FILES=$(find "$ANTIBMPNN_DIR" -name "*.pt" -type f 2>/dev/null)
PT_COUNT=$(echo "$PT_FILES" | wc -l)

if [ "$PT_COUNT" -eq 0 ]; then
    echo "Error: No .pt files found in extracted archive."
    echo "Contents of $ANTIBMPNN_DIR:"
    ls -la "$ANTIBMPNN_DIR"
    rm -f "$ZIP_FILE"
    exit 1
fi

# Symlink the first .pt file found to the expected target name
FIRST_PT=$(echo "$PT_FILES" | head -1)
echo "Found checkpoint: $FIRST_PT"

if [ "$FIRST_PT" != "$TARGET_WEIGHT" ]; then
    ln -sf "$FIRST_PT" "$TARGET_WEIGHT"
    echo "Linked to: $TARGET_WEIGHT"
fi

# Cleanup zip
rm -f "$ZIP_FILE"

echo ""
echo "============================================"
echo " Setup complete!"
echo "============================================"
echo ""
echo "AntiBMPNN weights installed at:"
echo "  $ANTIBMPNN_DIR"
echo ""
echo "Available weight files:"
find "$ANTIBMPNN_DIR" -name "*.pt" -type f -o -name "*.pt" -type l | while read f; do
    echo "  $(basename "$f")"
done
echo ""
echo "To use with RFantibody:"
echo "  proteinmpnn -i structures/ -o designed/   # Uses AntiBMPNN by default"
echo "  proteinmpnn -i structures/ -o designed/ -w weights/ProteinMPNN_v48_noise_0.2.pt  # Vanilla fallback"
echo ""
