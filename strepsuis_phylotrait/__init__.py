"""
StrepSuis-PhyloTrait: Phylogenetic Clustering and Binary Trait Analysis
========================================================================

A bioinformatics tool for phylogenetic clustering and binary trait analysis
in bacterial genomics.

Features:
    - Tree-aware clustering respecting phylogenetic structure
    - Phylogenetic diversity metrics (Faith's PD, MPD, MNTD)
    - Phylogenetic signal detection (D-statistic, permutation tests)
    - Trait correlation analysis with phylogenetic correction
    - UMAP visualization using phylogenetic distances
    - Association rule mining for trait co-occurrence

Example:
    >>> from strepsuis_phylotrait import PhyloTraitAnalyzer, Config
    >>> config = Config(
    ...     tree_file="tree.newick",
    ...     data_dir="./data",
    ...     output_dir="./output"
    ... )
    >>> analyzer = PhyloTraitAnalyzer(config)
    >>> results = analyzer.run()

Author: MK-vet
License: MIT
"""

__version__ = "1.0.0"
__author__ = "MK-vet"
__license__ = "MIT"

from .analyzer import PhyloTraitAnalyzer
from .config import Config

# Advanced statistical features from shared module
try:
    from shared.advanced_statistics import (
        multiview_concordance,
        outlier_consistency_index,
        bootstrap_stability_matrix,
        rare_pattern_detector,
        entropy_weighted_importance,
    )
    _HAS_ADVANCED_STATS = True
except ImportError:
    _HAS_ADVANCED_STATS = False

__all__ = ["PhyloTraitAnalyzer", "Config", "__version__"]

# Add advanced statistics if available
if _HAS_ADVANCED_STATS:
    __all__.extend([
        "multiview_concordance",
        "outlier_consistency_index",
        "bootstrap_stability_matrix",
        "rare_pattern_detector",
        "entropy_weighted_importance",
    ])
