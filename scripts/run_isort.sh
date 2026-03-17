#!/usr/bin/env bash
# Script to run isort on the codebase

set -e

echo "Running isort to organize imports..."

isort src/
isort scripts/
isort test/

echo "Import formatting complete!"
