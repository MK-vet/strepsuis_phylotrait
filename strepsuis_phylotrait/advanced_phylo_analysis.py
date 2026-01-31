#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Advanced Phylogenetic Analysis Module for StrepSuis-PhyloTrait
=============================================================

This module integrates advanced statistical methods specifically for
phylogenetic and clustering analysis.

Features:
    - Multi-view Concordance: Compare phylogenetic vs clustering groupings
    - Outlier Consistency: Identify outlier samples using multiple methods
    - Bootstrap Stability: Test robustness of phylogenetic clusters
    - Rare Pattern Detection: Find rare but consistent trait combinations
    - Entropy-weighted Importance: Adjust trait importance by information content

Author: MK-vet
Version: 1.0.0
License: MIT
"""

from typing import Dict, List, Optional
import pandas as pd
import numpy as np
import logging

try:
    from shared.advanced_statistics import (
        multiview_concordance,
        outlier_consistency_index,
        bootstrap_stability_matrix,
        rare_pattern_detector,
        entropy_weighted_importance,
    )
    HAS_ADVANCED_STATS = True
except ImportError:
    HAS_ADVANCED_STATS = False
    logging.warning(
        "Advanced statistics module not available. "
        "Install with: pip install -e ../shared"
    )

logger = logging.getLogger(__name__)


def phylogenetic_multiview_analysis(
    phylo_clusters: np.ndarray,
    trait_clusters: np.ndarray,
    network_communities: Optional[np.ndarray] = None,
) -> Dict:
    """
    Compare agreement between phylogenetic clustering and other analysis views.

    This function uses Normalized Mutual Information (NMI) and Adjusted Rand Index (ARI)
    to quantify how well different analysis methods agree with each other.

    Parameters
    ----------
    phylo_clusters : np.ndarray
        Cluster assignments from phylogenetic tree-aware clustering
    trait_clusters : np.ndarray
        Cluster assignments from trait-based clustering (e.g., K-Modes)
    network_communities : np.ndarray, optional
        Community assignments from network analysis

    Returns
    -------
    dict
        Dictionary with pairwise concordance metrics and interpretation

    Example
    -------
    >>> from strepsuis_phylotrait.advanced_phylo_analysis import phylogenetic_multiview_analysis
    >>>
    >>> # After running phylogenetic clustering
    >>> phylo_labels = tree_aware_clustering(tree, traits)
    >>> trait_labels = kmodes_clustering(traits)
    >>>
    >>> concordance = phylogenetic_multiview_analysis(
    ...     phylo_clusters=phylo_labels,
    ...     trait_clusters=trait_labels
    ... )
    >>>
    >>> print(f"Agreement NMI: {concordance['overall']['mean_nmi']:.3f}")
    >>> print(f"Interpretation: {concordance['overall']['interpretation']}")
    """
    if not HAS_ADVANCED_STATS:
        raise ImportError(
            "Advanced statistics module is required. "
            "Install shared module: pip install -e ../shared"
        )

    logger.info("Computing multi-view concordance...")
    results = multiview_concordance(
        clustering_labels=trait_clusters,
        network_communities=phylo_clusters,
        association_groups=network_communities
    )

    return results


def detect_phylogenetic_outliers(
    trait_data: pd.DataFrame,
    contamination: float = 0.1,
    density_percentile: float = 10
) -> Dict:
    """
    Identify outlier samples using multiple detection methods.

    Combines Isolation Forest and density-based outlier detection to find
    samples that are outliers by both methods (high consistency).

    Parameters
    ----------
    trait_data : pd.DataFrame
        Binary trait matrix (samples × traits)
    contamination : float, default=0.1
        Expected proportion of outliers
    density_percentile : float, default=10
        Percentile threshold for low density

    Returns
    -------
    dict
        Dictionary with outlier analysis results:
        - 'outlier_df': DataFrame with outlier scores and flags
        - 'consistent_outliers': List of samples flagged by both methods
        - 'n_outliers': Number of consistent outliers

    Example
    -------
    >>> outliers = detect_phylogenetic_outliers(trait_data, contamination=0.05)
    >>> consistent = outliers['consistent_outliers']
    >>> print(f"Found {len(consistent)} consistent outliers")
    """
    if not HAS_ADVANCED_STATS:
        raise ImportError(
            "Advanced statistics module is required. "
            "Install shared module: pip install -e ../shared"
        )

    logger.info("Detecting outliers with multi-method consistency...")
    outlier_df = outlier_consistency_index(
        data=trait_data,
        contamination=contamination,
        density_percentile=density_percentile
    )

    consistent_outliers = outlier_df[outlier_df['both_methods_outlier']].index.tolist()

    return {
        'outlier_df': outlier_df,
        'consistent_outliers': consistent_outliers,
        'n_outliers': len(consistent_outliers),
        'isolation_forest_only': outlier_df[
            outlier_df['isolation_forest_outlier'] & ~outlier_df['density_outlier']
        ].index.tolist(),
        'density_only': outlier_df[
            ~outlier_df['isolation_forest_outlier'] & outlier_df['density_outlier']
        ].index.tolist(),
    }


def phylogenetic_trait_stability(
    trait_data: pd.DataFrame,
    tree,
    clustering_func,
    n_bootstrap: int = 500,
    n_jobs: int = -1
) -> pd.DataFrame:
    """
    Test stability of trait-cluster associations via bootstrapping.

    Parameters
    ----------
    trait_data : pd.DataFrame
        Binary trait matrix
    tree
        Phylogenetic tree object (Bio.Phylo tree)
    clustering_func : callable
        Function that takes (tree, trait_data) and returns cluster labels
    n_bootstrap : int, default=500
        Number of bootstrap iterations
    n_jobs : int, default=-1
        Number of parallel jobs

    Returns
    -------
    pd.DataFrame
        Stability percentages for each trait

    Example
    -------
    >>> def my_clustering(tree, data):
    ...     # Your tree-aware clustering
    ...     return cluster_labels
    >>>
    >>> stability = phylogenetic_trait_stability(
    ...     trait_data=traits,
    ...     tree=phylo_tree,
    ...     clustering_func=my_clustering,
    ...     n_bootstrap=500
    ... )
    >>> print(stability.sort_values(ascending=False).head(10))
    """
    if not HAS_ADVANCED_STATS:
        raise ImportError(
            "Advanced statistics module is required. "
            "Install shared module: pip install -e ../shared"
        )

    logger.info(f"Computing trait stability across {n_bootstrap} bootstrap samples...")

    def analyze_significance(boot_data):
        """Identify significant traits in bootstrap sample."""
        try:
            # Run clustering on bootstrap sample
            boot_labels = clustering_func(tree, boot_data)

            # Simple chi-square test for trait-cluster association
            from scipy.stats import chi2_contingency
            significant_traits = []

            for trait in boot_data.columns:
                try:
                    contingency = pd.crosstab(boot_data[trait], boot_labels)
                    _, p_value, _, _ = chi2_contingency(contingency)
                    if p_value < 0.05:
                        significant_traits.append(trait)
                except (ValueError, TypeError, np.linalg.LinAlgError) as e:

                    logger.warning(f"Operation failed: {e}")
                    pass

            return significant_traits
        except (ValueError, TypeError, np.linalg.LinAlgError) as e:

            logger.warning(f"Operation failed: {e}")
            return []

    stability = bootstrap_stability_matrix(
        data=trait_data,
        analysis_func=analyze_significance,
        n_bootstrap=n_bootstrap,
        n_jobs=n_jobs
    )

    return stability


def find_rare_phylogenetic_patterns(
    trait_data: pd.DataFrame,
    min_support: float = 0.02,
    max_support: float = 0.05,
    min_confidence: float = 0.8,
    min_samples: int = 3
) -> List[Dict]:
    """
    Detect rare but consistent trait combinations (micro-phenotypes).

    These patterns have low support (few samples) but high confidence,
    representing potential strain-specific or rare phenotypic signatures.

    Parameters
    ----------
    trait_data : pd.DataFrame
        Binary trait matrix
    min_support : float, default=0.02
        Minimum support threshold (2%)
    max_support : float, default=0.05
        Maximum support for rare patterns (5%)
    min_confidence : float, default=0.8
        Minimum confidence threshold (80%)
    min_samples : int, default=3
        Minimum absolute sample count

    Returns
    -------
    list of dict
        List of rare patterns with metrics

    Example
    -------
    >>> rare_patterns = find_rare_phylogenetic_patterns(
    ...     trait_data=traits,
    ...     min_support=0.02,
    ...     max_support=0.05
    ... )
    >>> for pattern in rare_patterns[:5]:
    ...     print(pattern['interpretation'])
    """
    if not HAS_ADVANCED_STATS:
        raise ImportError(
            "Advanced statistics module is required. "
            "Install shared module: pip install -e ../shared"
        )

    logger.info("Searching for rare phylogenetic patterns...")
    patterns = rare_pattern_detector(
        data=trait_data,
        min_support=min_support,
        max_support=max_support,
        min_confidence=min_confidence,
        min_samples=min_samples
    )

    logger.info(f"Found {len(patterns)} rare patterns")
    return patterns


def adjust_trait_importance_by_entropy(
    trait_data: pd.DataFrame,
    base_importance: pd.Series,
    entropy_penalty: float = 0.5
) -> pd.DataFrame:
    """
    Adjust trait importance scores by their entropy (information content).

    Traits that are nearly constant (low entropy) are down-weighted, as they
    provide less discriminatory power.

    Parameters
    ----------
    trait_data : pd.DataFrame
        Binary trait matrix
    base_importance : pd.Series
        Initial importance scores (e.g., from Random Forest)
    entropy_penalty : float, default=0.5
        How much to penalize low-entropy traits (0-1)

    Returns
    -------
    pd.DataFrame
        DataFrame with original and entropy-adjusted importance

    Example
    -------
    >>> from sklearn.ensemble import RandomForestClassifier
    >>>
    >>> # Get base importance from Random Forest
    >>> rf = RandomForestClassifier(n_estimators=100, random_state=42)
    >>> rf.fit(traits, labels)
    >>> base_imp = pd.Series(rf.feature_importances_, index=traits.columns)
    >>>
    >>> # Adjust by entropy
    >>> adjusted = adjust_trait_importance_by_entropy(
    ...     trait_data=traits,
    ...     base_importance=base_imp
    ... )
    >>> print(adjusted.head(10))
    """
    if not HAS_ADVANCED_STATS:
        raise ImportError(
            "Advanced statistics module is required. "
            "Install shared module: pip install -e ../shared"
        )

    logger.info("Adjusting trait importance by entropy...")
    adjusted = entropy_weighted_importance(
        data=trait_data,
        base_importance=base_importance,
        entropy_penalty=entropy_penalty
    )

    return adjusted


def generate_phylo_advanced_report(
    multiview_results: Optional[Dict] = None,
    outlier_results: Optional[Dict] = None,
    stability_results: Optional[pd.DataFrame] = None,
    rare_patterns: Optional[List[Dict]] = None,
    output_file: str = "phylo_advanced_report.html"
) -> None:
    """
    Generate HTML report for advanced phylogenetic analysis.

    Parameters
    ----------
    multiview_results : dict, optional
        Results from phylogenetic_multiview_analysis()
    outlier_results : dict, optional
        Results from detect_phylogenetic_outliers()
    stability_results : pd.DataFrame, optional
        Results from phylogenetic_trait_stability()
    rare_patterns : list, optional
        Results from find_rare_phylogenetic_patterns()
    output_file : str, default='phylo_advanced_report.html'
        Output file path
    """
    html = ["<html><head><title>Phylogenetic Advanced Analysis</title>"]
    html.append("<style>")
    html.append("body { font-family: Arial, sans-serif; margin: 20px; }")
    html.append("table { border-collapse: collapse; width: 100%; margin: 20px 0; }")
    html.append("th, td { border: 1px solid #ddd; padding: 8px; }")
    html.append("th { background-color: #2E7D32; color: white; }")
    html.append("h1 { color: #1B5E20; }")
    html.append("h2 { color: #388E3C; }")
    html.append(".metric { background-color: #E8F5E9; padding: 10px; margin: 10px 0; border-radius: 5px; }")
    html.append("</style></head><body>")

    html.append("<h1>Advanced Phylogenetic Analysis Report</h1>")

    # Multi-view concordance
    if multiview_results:
        html.append("<h2>Multi-View Concordance</h2>")
        overall = multiview_results['overall']
        html.append(f"<div class='metric'>")
        html.append(f"<strong>Mean NMI:</strong> {overall['mean_nmi']:.3f}<br>")
        html.append(f"<strong>Mean ARI:</strong> {overall['mean_ari']:.3f}<br>")
        html.append(f"<strong>Interpretation:</strong> {overall['interpretation']}")
        html.append(f"</div>")

    # Outliers
    if outlier_results:
        html.append("<h2>Outlier Detection</h2>")
        html.append(f"<p>Detected <strong>{outlier_results['n_outliers']}</strong> consistent outliers</p>")
        if outlier_results['n_outliers'] > 0:
            html.append(f"<p>Outlier samples: {outlier_results['consistent_outliers'][:20]}</p>")

    # Stability
    if stability_results is not None:
        html.append("<h2>Trait Stability (Top 15)</h2>")
        top_stable = stability_results.sort_values(ascending=False).head(15)
        html.append(top_stable.to_frame('Stability (%)').to_html())

    # Rare patterns
    if rare_patterns:
        html.append("<h2>Rare Phylogenetic Patterns</h2>")
        html.append(f"<p>Found <strong>{len(rare_patterns)}</strong> rare but consistent patterns</p>")
        html.append("<ul>")
        for pattern in rare_patterns[:10]:
            html.append(f"<li>{pattern['interpretation']}</li>")
        html.append("</ul>")

    html.append("</body></html>")

    with open(output_file, 'w') as f:
        f.write('\n'.join(html))

    logger.info(f"Phylogenetic advanced analysis report saved to: {output_file}")
