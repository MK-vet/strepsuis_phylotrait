#!/usr/bin/env python3
"""
strepsuis_phylotrait.py - Canonical Entry Point
================================================

Phylogenetic clustering and trait analysis

Usage:
    python strepsuis_phylotrait.py --data-dir input/raw_data --tree phylogeny.newick --output out/run_20260131
    python strepsuis_phylotrait.py --config config.yaml

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
"""

import sys
from strepsuis_phylotrait.cli import main

if __name__ == '__main__':
    sys.exit(main())
