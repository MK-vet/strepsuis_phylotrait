#!/usr/bin/env python3
"""
Phylogenetic_clustering_2025_03_21.py - Canonical Entry Point
==============================================================

StrepSuis-PhyloTrait: Phylogenetic clustering and trait analysis

This is the canonical entry point for workflow orchestration (Nextflow/Snakemake).
It wraps the strepsuis_phylotrait package with standardized I/O following the
4-layer architecture.

Usage:
    python Phylogenetic_clustering_2025_03_21.py --config config.yaml
    python Phylogenetic_clustering_2025_03_21.py --data-dir input/raw_data --tree phylogeny.newick --output out/run_20260131

Architecture Compliance:
    Input:  input/raw_data/*.csv + phylogeny.newick + config.yaml
    Output: out/run_<ID>/
            ├── manifest.json
            ├── summary.json
            ├── results/*.parquet
            ├── figures/*.png
            ├── exports/*.csv
            ├── report.pdf
            └── site/

Module ID: StrepSuis-PhyloTrait
Canonical Name: Phylogenetic_clustering_2025_03_21.py
Date: 2025-03-21
"""

import sys
import argparse
from pathlib import Path

# Import the actual implementation from the package
from strepsuis_phylotrait.cli import main as cli_main

def main():
    """
    Canonical entry point with architecture-compliant defaults.
    """
    parser = argparse.ArgumentParser(
        description='StrepSuis-PhyloTrait: Phylogenetic clustering and trait analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument(
        '--config',
        type=Path,
        help='Path to config.yaml (default: input/config.yaml)'
    )

    parser.add_argument(
        '--data-dir',
        type=Path,
        help='Input data directory (default: input/raw_data)'
    )

    parser.add_argument(
        '--tree',
        type=Path,
        help='Phylogenetic tree in Newick format (default: input/raw_data/phylogeny.newick)'
    )

    parser.add_argument(
        '--output',
        type=Path,
        help='Output directory (default: out/run_<timestamp>)'
    )

    parser.add_argument(
        '--run-id',
        type=str,
        help='Run identifier for output directory (default: timestamp)'
    )

    # Pass through to underlying CLI
    sys.exit(cli_main())


if __name__ == '__main__':
    main()
