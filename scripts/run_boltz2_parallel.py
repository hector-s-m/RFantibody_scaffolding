#!/usr/bin/env python3
"""
Run Boltz2 predictions in parallel batches on a single GPU.

Spawns up to BATCH_SIZE concurrent `boltz predict` processes, each
processing a single YAML file. This overlaps MSA server I/O with GPU
computation, significantly reducing total wall-clock time.

Usage:
    python scripts/run_boltz2_parallel.py \
        -i boltz2_input/ \
        -o boltz2_output/ \
        --samples 3 \
        --batch-size 6

Designed to be called from the motif scaffolding pipeline script.
Requires the boltz conda environment to already be active.
"""

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path


def spawn_boltz_predict(
    yaml_path: str,
    out_dir: str,
    diffusion_samples: int,
    msa_server_url: str,
    cache: str = "",
    extra_args: list = None,
) -> subprocess.Popen:
    """Spawn a non-blocking `boltz predict` process for a single YAML file."""
    cmd = [
        "boltz", "predict", yaml_path,
        "--out_dir", out_dir,
        "--diffusion_samples", str(diffusion_samples),
        "--output_format", "pdb",
        "--use_msa_server",
        "--msa_server_url", msa_server_url,
        "--use_potentials",
        "--write_full_pae",
        "--write_full_pde",
    ]
    if cache:
        cmd.extend(["--cache", cache])
    if extra_args:
        cmd.extend(extra_args)

    return subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )


def main():
    parser = argparse.ArgumentParser(
        description="Run Boltz2 predictions in parallel batches"
    )
    parser.add_argument(
        "-i", "--input-dir", type=str, required=True,
        help="Directory containing Boltz2 YAML files"
    )
    parser.add_argument(
        "-o", "--output-dir", type=str, required=True,
        help="Output directory for predictions"
    )
    parser.add_argument(
        "--samples", type=int, default=3,
        help="Diffusion samples per design (default: 3)"
    )
    parser.add_argument(
        "--batch-size", type=int, default=6,
        help="Number of parallel predictions (default: 6)"
    )
    parser.add_argument(
        "--msa-server-url", type=str, default="http://a3m-2023.mmseqs.com",
        help="MSA server URL"
    )
    parser.add_argument(
        "--cache", type=str, default="",
        help="Boltz2 cache directory"
    )
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    yaml_files = sorted(input_dir.glob("*.yaml"))
    if not yaml_files:
        print(f"Error: No YAML files found in {input_dir}")
        sys.exit(1)

    total = len(yaml_files)
    batch_size = args.batch_size
    print(f"Running {total} Boltz2 predictions in batches of {batch_size}")
    print(f"  Input:   {input_dir}")
    print(f"  Output:  {output_dir}")
    print(f"  Samples: {args.samples} per design")
    print()

    completed = 0
    failed = 0
    start_time = time.time()

    for batch_start in range(0, total, batch_size):
        batch = yaml_files[batch_start:batch_start + batch_size]
        batch_num = batch_start // batch_size + 1
        total_batches = (total + batch_size - 1) // batch_size

        print(f"[Batch {batch_num}/{total_batches}] Spawning {len(batch)} predictions...")

        # Spawn all processes in this batch
        procs = {}
        for yaml_path in batch:
            stem = yaml_path.stem
            proc = spawn_boltz_predict(
                yaml_path=str(yaml_path),
                out_dir=str(output_dir),
                diffusion_samples=args.samples,
                msa_server_url=args.msa_server_url,
                cache=args.cache,
            )
            procs[stem] = proc

        # Wait for all processes in this batch
        for stem, proc in procs.items():
            rc = proc.wait()
            if rc == 0:
                completed += 1
                print(f"  OK: {stem} ({completed}/{total})")
            else:
                failed += 1
                # Capture output for debugging
                stdout = proc.stdout.read().decode() if proc.stdout else ""
                print(f"  FAILED: {stem} (exit code {rc})")
                if stdout:
                    # Print last few lines of output
                    lines = stdout.strip().split("\n")
                    for line in lines[-5:]:
                        print(f"    {line}")

        elapsed = time.time() - start_time
        rate = completed / elapsed if elapsed > 0 else 0
        remaining = (total - completed - failed) / rate if rate > 0 else 0
        print(f"  Elapsed: {elapsed/60:.1f}min, "
              f"~{remaining/60:.1f}min remaining")
        print()

    elapsed = time.time() - start_time
    print(f"Boltz2 parallel predictions complete:")
    print(f"  Completed: {completed}/{total}")
    print(f"  Failed:    {failed}/{total}")
    print(f"  Time:      {elapsed/60:.1f} min")

    if failed > 0:
        print(f"\nWARNING: {failed} predictions failed. Check output above.")
        sys.exit(1)


if __name__ == "__main__":
    main()
