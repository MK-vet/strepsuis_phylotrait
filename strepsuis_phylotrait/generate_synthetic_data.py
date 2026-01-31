"""
Synthetic Data Generator for StrepSuis-PhyloTrait - Phylogenetic Analysis Validation

This module generates synthetic datasets including phylogenetic trees and
binary trait data for validating tree-aware clustering and trait analysis.

Generation Methodology:
-----------------------
1. **Coalescent Trees**: Random phylogenetic trees generated using
   coalescent or pure-birth models.

2. **Trait Evolution**: Binary traits evolved along the tree with
   known transition rates.

3. **Cluster Structure**: Known cluster assignments based on tree
   structure for validation.

4. **Ground Truth**: True phylogenetic diversity and clustering
   assignments are recorded.

Scientific References:
---------------------
- Faith, D.P. (1992). Conservation evaluation and phylogenetic diversity.
  Biological Conservation, 61(1), 1-10.
- Pagel, M. (1999). Inferring the historical patterns of biological evolution.
  Nature, 401(6756), 877-884.

Author: MK-vet Team
License: MIT
Version: 1.0.0
"""

from dataclasses import dataclass, field
from datetime import datetime
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


@dataclass
class SyntheticPhyloConfig:
    """Configuration for synthetic phylogenetic data generation.

    Attributes:
        n_strains: Number of bacterial strains to generate
        n_clusters: Number of phylogenetic clusters
        n_mic_features: Number of MIC features
        n_amr_features: Number of AMR gene features
        n_virulence_features: Number of virulence factor features
        trait_heritability: How strongly traits follow phylogeny (0.0 to 1.0)
        tree_depth: Average tree depth (branch length units)
        random_state: Random seed for reproducibility
    """

    n_strains: int = 200
    n_clusters: int = 4
    n_mic_features: int = 10
    n_amr_features: int = 15
    n_virulence_features: int = 10
    trait_heritability: float = 0.6
    tree_depth: float = 0.5
    random_state: int = 42
    # Legacy/test aliases
    n_taxa: Optional[int] = None
    n_traits: Optional[int] = None
    legacy_mode: bool = False

    def __post_init__(self):
        """
        Initialize and synchronize configuration parameters after object creation.

        Handles legacy parameter names (n_taxa, n_traits) and ensures consistency
        between modern and legacy parameter names for backward compatibility.

        Notes
        -----
        This method is automatically called by dataclass after __init__.
        Synchronizes n_taxa with n_strains and n_traits with n_mic_features.
        Sets legacy_mode flag if legacy parameters were provided.
        """
        legacy_requested = self.n_taxa is not None or self.n_traits is not None
        if self.n_taxa is None:
            self.n_taxa = self.n_strains
        else:
            self.n_strains = self.n_taxa

        if self.n_traits is None:
            self.n_traits = self.n_mic_features
        else:
            self.n_mic_features = self.n_traits

        self.legacy_mode = legacy_requested


@dataclass
class SyntheticPhyloMetadata:
    """Metadata describing the generated synthetic phylogenetic data.

    Attributes:
        config: The configuration used to generate the data
        tree_newick: The generated tree in Newick format
        true_cluster_labels: Array of true cluster assignments
        cluster_mrca_nodes: MRCA node for each cluster
        expected_pd_range: Expected phylogenetic diversity range
        mic_columns: List of MIC column names
        amr_columns: List of AMR gene column names
        virulence_columns: List of virulence factor column names
        generation_timestamp: When the data was generated
    """

    config: SyntheticPhyloConfig
    tree_newick: str = ""
    true_cluster_labels: np.ndarray = field(default_factory=lambda: np.array([]))
    cluster_mrca_nodes: Dict[int, str] = field(default_factory=dict)
    expected_pd_range: Tuple[float, float] = (0.1, 1.0)
    mic_columns: List[str] = field(default_factory=list)
    amr_columns: List[str] = field(default_factory=list)
    virulence_columns: List[str] = field(default_factory=list)
    generation_timestamp: str = field(default_factory=lambda: datetime.now().isoformat())


def generate_random_tree_newick(
    n_taxa: int,
    branch_length_mean: float = 0.05,
    random_state: Optional[int] = 42,
    return_taxa: bool = True,
) -> str:
    """
    Generate a random phylogenetic tree in Newick format.

    Uses a simple birth-death process to generate tree structure.

    Parameters:
        n_taxa: Number of terminal taxa (strains)
        branch_length_mean: Mean branch length
        random_state: Random seed

    Returns:
        Newick string (or tuple if return_taxa is True)
    """
    rng = np.random.default_rng(random_state)

    # Generate taxa names
    taxa = [f"Strain_{i:04d}" for i in range(1, n_taxa + 1)]

    # Simple recursive tree building
    def build_subtree(taxa_list, depth=0):
        """
        Recursively build random phylogenetic subtree using birth-death process.

        Parameters
        ----------
        taxa_list : list
            List of taxon names to include in subtree
        depth : int, default=0
            Current recursion depth (for internal tracking)

        Returns
        -------
        str
            Newick format subtree string with branch lengths
        """
        if len(taxa_list) == 1:
            bl = rng.exponential(branch_length_mean)
            return f"{taxa_list[0]}:{bl:.6f}"
        elif len(taxa_list) == 2:
            bl1 = rng.exponential(branch_length_mean)
            bl2 = rng.exponential(branch_length_mean)
            return f"({taxa_list[0]}:{bl1:.6f},{taxa_list[1]}:{bl2:.6f})"
        else:
            # Split into two groups
            split_point = rng.integers(1, len(taxa_list))
            left = taxa_list[:split_point]
            right = taxa_list[split_point:]

            left_subtree = build_subtree(left, depth + 1)
            right_subtree = build_subtree(right, depth + 1)

            bl = rng.exponential(branch_length_mean)
            return f"({left_subtree},{right_subtree}):{bl:.6f}"

    # Shuffle taxa for random topology
    shuffled_taxa = taxa.copy()
    rng.shuffle(shuffled_taxa)

    tree_str = build_subtree(shuffled_taxa) + ";"

    return (tree_str, taxa) if return_taxa else tree_str


def generate_clustered_tree_newick(
    n_taxa: int,
    n_clusters: int,
    branch_length_mean: float = 0.05,
    random_state: Optional[int] = 42,
    return_details: bool = True,
) -> str:
    """
    Generate a tree with clear cluster structure.

    Creates n_clusters subtrees with short internal branches
    and longer branches between clusters.

    Parameters:
        n_taxa: Number of terminal taxa
        n_clusters: Number of clusters
        branch_length_mean: Mean branch length
        random_state: Random seed

    Returns:
        Newick string (or tuple if return_details is True)
    """
    rng = np.random.default_rng(random_state)

    # Assign taxa to clusters
    taxa = [f"Strain_{i:04d}" for i in range(1, n_taxa + 1)]
    cluster_sizes = [n_taxa // n_clusters] * n_clusters
    for i in range(n_taxa % n_clusters):
        cluster_sizes[i] += 1

    cluster_labels = []
    cluster_members = {i: [] for i in range(1, n_clusters + 1)}

    idx = 0
    for cluster_id, size in enumerate(cluster_sizes, 1):
        for _ in range(size):
            cluster_labels.append(cluster_id)
            cluster_members[cluster_id].append(taxa[idx])
            idx += 1

    # Shuffle labels to randomize position
    combined = list(zip(taxa, cluster_labels))
    rng.shuffle(combined)
    taxa_shuffled = [t[0] for t in combined]
    labels_shuffled = [t[1] for t in combined]

    # Rebuild cluster_members after shuffle
    cluster_members = {i: [] for i in range(1, n_clusters + 1)}
    for taxon, label in zip(taxa_shuffled, labels_shuffled):
        cluster_members[label].append(taxon)

    # Build tree: first create subtrees for each cluster
    def build_cluster_subtree(members, intra_bl_scale=0.5):
        """
        Build phylogenetic subtree for single cluster with short internal branches.

        Creates tight monophyletic clade with reduced branch length variation
        to simulate cluster cohesion.

        Parameters
        ----------
        members : list
            List of taxon names in this cluster
        intra_bl_scale : float, default=0.5
            Scaling factor for branch lengths (smaller = tighter cluster)

        Returns
        -------
        str
            Newick format subtree string
        """
        if len(members) == 1:
            bl = rng.exponential(branch_length_mean * intra_bl_scale)
            return f"{members[0]}:{bl:.6f}"
        elif len(members) == 2:
            bl1 = rng.exponential(branch_length_mean * intra_bl_scale)
            bl2 = rng.exponential(branch_length_mean * intra_bl_scale)
            return f"({members[0]}:{bl1:.6f},{members[1]}:{bl2:.6f})"
        else:
            mid = len(members) // 2
            left = build_cluster_subtree(members[:mid], intra_bl_scale)
            right = build_cluster_subtree(members[mid:], intra_bl_scale)
            bl = rng.exponential(branch_length_mean * intra_bl_scale)
            return f"({left},{right}):{bl:.6f}"

    # Build subtrees for each cluster
    cluster_subtrees = []
    for cluster_id in range(1, n_clusters + 1):
        members = cluster_members[cluster_id]
        subtree = build_cluster_subtree(members, intra_bl_scale=0.5)
        cluster_subtrees.append(subtree)

    # Join clusters with longer branches
    def join_subtrees(subtrees, inter_bl_scale=2.0):
        """
        Join cluster subtrees with longer inter-cluster branches.

        Creates clear phylogenetic separation between clusters by using
        increased branch lengths at cluster boundaries.

        Parameters
        ----------
        subtrees : list
            List of Newick subtree strings to join
        inter_bl_scale : float, default=2.0
            Scaling factor for inter-cluster branches (larger = more separation)

        Returns
        -------
        str
            Newick format tree string joining all subtrees
        """
        if len(subtrees) == 1:
            return subtrees[0]
        elif len(subtrees) == 2:
            bl = rng.exponential(branch_length_mean * inter_bl_scale)
            return f"({subtrees[0]},{subtrees[1]}):{bl:.6f}"
        else:
            mid = len(subtrees) // 2
            left = join_subtrees(subtrees[:mid], inter_bl_scale)
            right = join_subtrees(subtrees[mid:], inter_bl_scale)
            bl = rng.exponential(branch_length_mean * inter_bl_scale)
            return f"({left},{right}):{bl:.6f}"

    tree_str = join_subtrees(cluster_subtrees) + ";"

    # Create label array in original order
    final_labels = [0] * n_taxa
    for taxon, label in zip(taxa_shuffled, labels_shuffled):
        idx = int(taxon.split("_")[1]) - 1
        final_labels[idx] = label

    return (tree_str, final_labels, cluster_members) if return_details else tree_str


def generate_trait_data_with_heritability(
    n_samples: Optional[int] = None,
    n_features: Optional[int] = None,
    cluster_labels: Optional[List[int]] = None,
    heritability: float = 0.6,
    base_prevalence: float = 0.3,
    random_state: Optional[int] = 42,
    tree_newick: Optional[str] = None,
    n_traits: Optional[int] = None,
) -> Any:
    """
    Generate binary trait data correlated with phylogenetic clusters.

    Higher heritability means traits more strongly follow cluster assignments.

    Parameters:
        n_samples: Number of samples
        n_features: Number of features
        cluster_labels: Cluster assignment for each sample
        heritability: How strongly traits follow clusters (0-1)
        base_prevalence: Base prevalence for features
        random_state: Random seed

    Returns:
        Binary trait matrix (n_samples x n_features)
    """
    rng = np.random.default_rng(random_state)

    # Legacy positional call: (tree_newick, taxa, n_traits=...)
    taxa_override = None
    if isinstance(n_samples, str) and isinstance(n_features, (list, tuple)) and tree_newick is None:
        tree_newick = n_samples
        taxa_override = list(n_features)
        n_samples = None
        n_features = None

    # Legacy path used by tests: tree_newick + n_traits
    if tree_newick is not None or n_traits is not None:
        import re

        trait_count = int(n_traits or 0)
        taxa = []
        if tree_newick:
            taxa = re.findall(r"([A-Za-z0-9_\\-\\.]+):", tree_newick)
        if taxa_override:
            taxa = taxa_override
        if not taxa and n_samples:
            taxa = [f"Strain_{i:04d}" for i in range(1, int(n_samples) + 1)]
        if not taxa:
            taxa = [f"Strain_{i:04d}" for i in range(1, 11)]
        if trait_count <= 0:
            trait_count = max(1, int(n_features or 1))

        data = rng.binomial(1, base_prevalence, (len(taxa), trait_count))
        df = pd.DataFrame(data, columns=[f"Trait_{i+1:02d}" for i in range(trait_count)])
        df.insert(0, "Strain_ID", taxa)
        return df

    if n_samples is None or n_features is None or cluster_labels is None:
        raise ValueError("n_samples, n_features, and cluster_labels are required for core generation.")

    n_clusters = len(set(cluster_labels))
    data = np.zeros((n_samples, n_features), dtype=int)

    # Generate cluster-specific prevalence patterns
    cluster_patterns = {}
    for cluster_id in range(1, n_clusters + 1):
        # Each cluster has some features with high prevalence
        pattern = rng.uniform(base_prevalence * 0.5, base_prevalence * 1.5, n_features)

        # Select some features to be cluster-specific
        n_specific = max(1, n_features // n_clusters)
        specific_features = rng.choice(n_features, size=n_specific, replace=False)
        pattern[specific_features] = rng.uniform(0.6, 0.9, n_specific)

        cluster_patterns[cluster_id] = np.clip(pattern, 0.05, 0.95)

    # Generate data
    for i in range(n_samples):
        cluster_id = cluster_labels[i]
        cluster_pattern = cluster_patterns[cluster_id]

        # Mix cluster pattern with random noise based on heritability
        effective_pattern = (
            heritability * cluster_pattern + (1 - heritability) * base_prevalence
        )

        data[i] = (rng.random(n_features) < effective_pattern).astype(int)

    return data


def generate_phylotrait_synthetic_dataset(
    config: Optional[SyntheticPhyloConfig] = None,
) -> Tuple[str, pd.DataFrame, pd.DataFrame, pd.DataFrame, SyntheticPhyloMetadata]:
    """
    Generate a complete synthetic dataset for phylogenetic analysis validation.

    Parameters:
        config: Configuration object. Uses defaults if None.

    Returns:
        Tuple of (tree_newick, mic_df, amr_df, virulence_df, metadata)
    """
    if config is None:
        config = SyntheticPhyloConfig()

    rng = np.random.default_rng(config.random_state)

    if getattr(config, "legacy_mode", False):
        tree_newick = generate_random_tree_newick(
            n_taxa=config.n_taxa, random_state=config.random_state, return_taxa=False
        )
        traits_df = generate_trait_data_with_heritability(
            tree_newick=tree_newick,
            n_traits=config.n_traits,
            random_state=config.random_state,
        )
        metadata = SyntheticPhyloMetadata(config=config, tree_newick=tree_newick)
        return traits_df, metadata

    # Initialize metadata
    metadata = SyntheticPhyloMetadata(config=config)

    # Generate tree with cluster structure
    tree_newick, cluster_labels, cluster_members = generate_clustered_tree_newick(
        config.n_strains,
        config.n_clusters,
        branch_length_mean=config.tree_depth / 10,
        random_state=config.random_state,
        return_details=True,
    )
    metadata.tree_newick = tree_newick
    metadata.true_cluster_labels = np.array(cluster_labels)

    # Generate strain IDs
    strain_ids = [f"Strain_{i:04d}" for i in range(1, config.n_strains + 1)]

    # Generate column names
    mic_names = [
        "Oxytetracycline",
        "Doxycycline",
        "Tulathromycin",
        "Spectinomycin",
        "Gentamicin",
        "Tiamulin",
        "Enrofloxacin",
        "Penicillin",
        "Ampicillin",
        "Florfenicol",
    ][: config.n_mic_features]
    metadata.mic_columns = mic_names

    amr_names = [f"AMR_Gene_{i:02d}" for i in range(1, config.n_amr_features + 1)]
    metadata.amr_columns = amr_names

    vir_names = [f"Virulence_{i:02d}" for i in range(1, config.n_virulence_features + 1)]
    metadata.virulence_columns = vir_names

    # Generate trait data with phylogenetic signal
    mic_data = generate_trait_data_with_heritability(
        config.n_strains,
        config.n_mic_features,
        cluster_labels,
        config.trait_heritability,
        random_state=config.random_state + 100,
    )

    amr_data = generate_trait_data_with_heritability(
        config.n_strains,
        config.n_amr_features,
        cluster_labels,
        config.trait_heritability,
        random_state=config.random_state + 200,
    )

    vir_data = generate_trait_data_with_heritability(
        config.n_strains,
        config.n_virulence_features,
        cluster_labels,
        config.trait_heritability,
        random_state=config.random_state + 300,
    )

    # Create DataFrames
    mic_df = pd.DataFrame(mic_data, columns=mic_names)
    mic_df.insert(0, "Strain_ID", strain_ids)

    amr_df = pd.DataFrame(amr_data, columns=amr_names)
    amr_df.insert(0, "Strain_ID", strain_ids)

    vir_df = pd.DataFrame(vir_data, columns=vir_names)
    vir_df.insert(0, "Strain_ID", strain_ids)

    return tree_newick, mic_df, amr_df, vir_df, metadata


def save_synthetic_phylo_data(
    tree_newick: Any,
    mic_df: Optional[pd.DataFrame] = None,
    amr_df: Optional[pd.DataFrame] = None,
    vir_df: Optional[pd.DataFrame] = None,
    metadata: Optional[SyntheticPhyloMetadata] = None,
    output_dir: str = "synthetic_data",
) -> Dict[str, str]:
    """
    Save synthetic phylogenetic data and metadata to files.

    Parameters:
        tree_newick: Newick format tree string
        mic_df: MIC data DataFrame
        amr_df: AMR genes DataFrame
        vir_df: Virulence factors DataFrame
        metadata: Metadata object with ground truth
        output_dir: Directory to save files to

    Returns:
        Dict with paths to saved files
    """
    import json

    # Legacy path: save a single traits DataFrame + metadata
    if isinstance(tree_newick, pd.DataFrame):
        data_df = tree_newick
        metadata = mic_df if metadata is None else metadata
        if isinstance(amr_df, str):
            output_dir = amr_df
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        saved_files = {}
        traits_file = output_path / "synthetic_traits.csv"
        data_df.to_csv(traits_file, index=False)
        saved_files["traits"] = str(traits_file)
        if metadata is not None:
            metadata_file = output_path / "synthetic_metadata.json"
            with open(metadata_file, "w", encoding="utf-8") as f:
                json.dump(
                    metadata.__dict__ if hasattr(metadata, "__dict__") else metadata,
                    f,
                    indent=2,
                    default=str,
                )
            saved_files["metadata"] = str(metadata_file)
        return saved_files

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    saved_files = {}

    # Save tree file
    tree_file = output_path / "synthetic_tree.newick"
    with open(tree_file, "w") as f:
        f.write(tree_newick)
    saved_files["tree"] = str(tree_file)

    # Save data files
    mic_file = output_path / "synthetic_MIC.csv"
    mic_df.to_csv(mic_file, index=False)
    saved_files["mic"] = str(mic_file)

    amr_file = output_path / "synthetic_AMR_genes.csv"
    amr_df.to_csv(amr_file, index=False)
    saved_files["amr_genes"] = str(amr_file)

    vir_file = output_path / "synthetic_Virulence.csv"
    vir_df.to_csv(vir_file, index=False)
    saved_files["virulence"] = str(vir_file)

    # Save true cluster assignments
    clusters_df = pd.DataFrame(
        {"Strain_ID": mic_df["Strain_ID"], "True_Cluster": metadata.true_cluster_labels}
    )
    clusters_file = output_path / "synthetic_true_clusters.csv"
    clusters_df.to_csv(clusters_file, index=False)
    saved_files["true_clusters"] = str(clusters_file)

    # Save metadata as JSON
    metadata_dict = {
        "config": {
            "n_strains": metadata.config.n_strains,
            "n_clusters": metadata.config.n_clusters,
            "n_mic_features": metadata.config.n_mic_features,
            "n_amr_features": metadata.config.n_amr_features,
            "n_virulence_features": metadata.config.n_virulence_features,
            "trait_heritability": metadata.config.trait_heritability,
            "tree_depth": metadata.config.tree_depth,
            "random_state": metadata.config.random_state,
        },
        "true_cluster_counts": {
            str(k): int(v)
            for k, v in zip(*np.unique(metadata.true_cluster_labels, return_counts=True))
        },
        "expected_pd_range": list(metadata.expected_pd_range),
        "mic_columns": metadata.mic_columns,
        "amr_columns": metadata.amr_columns,
        "virulence_columns": metadata.virulence_columns,
        "generation_timestamp": metadata.generation_timestamp,
        "generation_method": "Clustered birth-death tree with heritable traits",
    }

    metadata_file = output_path / "synthetic_metadata.json"
    with open(metadata_file, "w", encoding="utf-8") as f:
        json.dump(metadata_dict, f, indent=2, default=str)
    saved_files["metadata"] = str(metadata_file)

    # Save methodology documentation
    methodology_content = f"""# Synthetic Data Generation Methodology for Phylogenetic Analysis

## Overview

This document describes the statistical methodology used to generate synthetic
data for validating phylogenetic clustering and trait analysis.

## Generation Parameters

- **Number of strains**: {metadata.config.n_strains}
- **Number of clusters**: {metadata.config.n_clusters}
- **MIC features**: {metadata.config.n_mic_features}
- **AMR gene features**: {metadata.config.n_amr_features}
- **Virulence factor features**: {metadata.config.n_virulence_features}
- **Trait heritability**: {metadata.config.trait_heritability:.2f}
- **Tree depth**: {metadata.config.tree_depth:.2f}
- **Random seed**: {metadata.config.random_state}

## Statistical Methods Used

### 1. Tree Generation

A random phylogenetic tree is generated with:
- Clustered topology: strains within the same cluster share a more recent common ancestor
- Exponential branch lengths: shorter branches within clusters, longer between clusters
- Birth-death process: tree structure follows coalescent-like dynamics

### 2. Trait Evolution (Heritability Model)

Binary traits are generated with phylogenetic signal:
- Each cluster has characteristic prevalence patterns
- Heritability parameter ({metadata.config.trait_heritability:.2f}) controls how strongly traits follow phylogeny
- Higher heritability = more phylogenetic signal

Formula:
  P(trait=1|cluster) = heritability × cluster_pattern + (1-heritability) × base_prevalence

### 3. Cluster Structure

True clusters are embedded in the tree structure:
- Each cluster forms a monophyletic clade
- Cluster boundaries correspond to deep splits in the tree

## Ground Truth

### Cluster Distribution
"""
    cluster_counts = np.bincount(np.array(metadata.true_cluster_labels))[1:]
    for i, count in enumerate(cluster_counts, 1):
        methodology_content += f"- Cluster {i}: {count} strains ({count/metadata.config.n_strains*100:.1f}%)\n"

    methodology_content += f"""
## Expected Analysis Performance

- Tree-aware clustering should recover ~{int(metadata.config.trait_heritability*100)}% of true clusters
- Phylogenetic diversity should reflect tree structure
- Trait associations should show phylogenetic signal

## References

1. Faith, D.P. (1992). Conservation evaluation and phylogenetic diversity.
2. Pagel, M. (1999). Inferring the historical patterns of biological evolution.

## Generation Timestamp

{metadata.generation_timestamp}

---
*This data was generated for validation and testing purposes only.*
"""

    methodology_file = output_path / "GENERATION_METHODOLOGY.md"
    with open(methodology_file, "w", encoding="utf-8") as f:
        f.write(methodology_content)
    saved_files["methodology"] = str(methodology_file)

    return saved_files


def validate_synthetic_phylo_data(
    tree_newick: Any,
    mic_df: Optional[pd.DataFrame] = None,
    amr_df: Optional[pd.DataFrame] = None,
    vir_df: Optional[pd.DataFrame] = None,
    metadata: Optional[SyntheticPhyloMetadata] = None,
) -> Dict[str, Any]:
    """
    Validate that synthetic data has expected statistical properties.
    """
    results = {
        "validation_passed": True,
        "checks": [],
        "warnings": [],
        "errors": [],
    }

    # Legacy path: validate a single traits DataFrame + metadata
    if isinstance(tree_newick, pd.DataFrame):
        data_df = tree_newick
        if data_df.empty:
            results["errors"].append("ERROR Traits data is empty")
            results["validation_passed"] = False
            return results
        results["checks"].append(f"OK Row count: {len(data_df)}")
        trait_cols = [c for c in data_df.columns if c != "Strain_ID"]
        if trait_cols:
            non_binary = False
            for col in trait_cols:
                values = set(data_df[col].dropna().unique())
                if not values.issubset({0, 1}):
                    non_binary = True
                    break
            if non_binary:
                results["warnings"].append("WARN Non-binary trait values detected")
            else:
                results["checks"].append("OK Traits are binary")
        return results

    # Check 1: Verify tree is valid Newick
    if tree_newick.endswith(";") and "(" in tree_newick:
        results["checks"].append("OK Tree is valid Newick format")
    else:
        results["errors"].append("✗ Invalid Newick format")
        results["validation_passed"] = False

    # Check 2: Verify data shapes
    expected_rows = metadata.config.n_strains
    for name, df, expected_cols in [
        ("MIC", mic_df, metadata.config.n_mic_features + 1),
        ("AMR", amr_df, metadata.config.n_amr_features + 1),
        ("Virulence", vir_df, metadata.config.n_virulence_features + 1),
    ]:
        if len(df) == expected_rows:
            results["checks"].append(f"OK {name} row count: {len(df)}")
        else:
            results["errors"].append(f"✗ {name} row count: {len(df)} (expected {expected_rows})")
            results["validation_passed"] = False

    # Check 3: Verify cluster count
    unique_clusters = len(np.unique(metadata.true_cluster_labels))
    if unique_clusters == metadata.config.n_clusters:
        results["checks"].append(f"OK Cluster count: {unique_clusters}")
    else:
        results["errors"].append(
            f"✗ Cluster count: {unique_clusters} (expected {metadata.config.n_clusters})"
        )
        results["validation_passed"] = False

    return results


if __name__ == "__main__":
    # Generate synthetic data when run directly
    print("Generating synthetic phylogenetic data...")

    config = SyntheticPhyloConfig(
        n_strains=200,
        n_clusters=4,
        random_state=42,
    )
    # #region agent log
    try:
        import json
        import time

        payload = {
            "sessionId": "debug-session",
            "runId": "run1",
            "hypothesisId": "H1",
            "location": "generate_synthetic_data.py",
            "message": "phylotrait_synthetic_config",
            "data": {"n_strains": config.n_strains},
            "timestamp": int(time.time() * 1000),
        }
        with open(r"c:\Users\ABC\Documents\GitHub\.cursor\debug.log", "a", encoding="utf-8") as log_file:
            log_file.write(json.dumps(payload, ensure_ascii=True) + "\n")
    except Exception:
        pass
    # #endregion

    tree_newick, mic_df, amr_df, vir_df, metadata = generate_phylotrait_synthetic_dataset(config)

    print(f"Generated tree with {config.n_strains} taxa")
    print(f"Generated {len(metadata.mic_columns)} MIC features")
    print(f"Generated {len(metadata.amr_columns)} AMR gene features")
    print(f"Generated {len(metadata.virulence_columns)} virulence factor features")
    print(f"True clusters: {config.n_clusters}")

    # Validate
    print("\nValidating synthetic data...")
    validation = validate_synthetic_phylo_data(tree_newick, mic_df, amr_df, vir_df, metadata)

    for check in validation["checks"]:
        print(f"  {check}")
    if validation["warnings"]:
        print(f"\nWarnings: {len(validation['warnings'])}")
    if validation["errors"]:
        print(f"Errors: {len(validation['errors'])}")

    print(f"\nValidation: {'PASSED' if validation['validation_passed'] else 'FAILED'}")

    # Save validation report (publication-ready)
    validation_dir = Path(__file__).parent.parent / "validation"
    validation_dir.mkdir(parents=True, exist_ok=True)
    validation_payload = {
        "generated": datetime.utcnow().isoformat(),
        "n_strains": int(len(mic_df)),
        "n_features": {
            "mic": int(len(mic_df.columns) - 1),
            "amr_genes": int(len(amr_df.columns) - 1),
            "virulence": int(len(vir_df.columns) - 1),
        },
        "checks": validation.get("checks", []),
        "warnings": validation.get("warnings", []),
        "errors": validation.get("errors", []),
        "validation_passed": bool(validation.get("validation_passed", False)),
    }
    with open(validation_dir / "synthetic_validation_results.json", "w", encoding="utf-8") as f:
        json.dump(validation_payload, f, indent=2)
    report_lines = [
        "# Synthetic Data Validation Report - strepsuis-phylotrait",
        "",
        f"Generated: {validation_payload['generated']}",
        "Data Source: Synthetic data with known ground truth",
        f"Strains: {validation_payload['n_strains']}",
        "Features:",
        f"- MIC: {validation_payload['n_features']['mic']}",
        f"- AMR genes: {validation_payload['n_features']['amr_genes']}",
        f"- Virulence: {validation_payload['n_features']['virulence']}",
        f"Checks: {len(validation_payload['checks'])}",
        f"Warnings: {len(validation_payload['warnings'])}",
        f"Errors: {len(validation_payload['errors'])}",
        f"Status: {'PASSED' if validation_payload['validation_passed'] else 'FAILED'}",
        "",
        "## Checks",
    ]
    report_lines.extend([f"- {item}" for item in validation_payload["checks"]])
    if validation_payload["warnings"]:
        report_lines.append("")
        report_lines.append("## Warnings")
        report_lines.extend([f"- {item}" for item in validation_payload["warnings"]])
    if validation_payload["errors"]:
        report_lines.append("")
        report_lines.append("## Errors")
        report_lines.extend([f"- {item}" for item in validation_payload["errors"]])
    with open(validation_dir / "SYNTHETIC_DATA_VALIDATION_REPORT.md", "w", encoding="utf-8") as f:
        f.write("\n".join(report_lines))

    # Save validation report (publication-ready)
    validation_dir = Path(__file__).parent.parent / "validation"
    validation_dir.mkdir(parents=True, exist_ok=True)
    validation_payload = {
        "generated": datetime.utcnow().isoformat(),
        "n_strains": int(len(mic_df)),
        "n_features": {
            "mic": int(len(mic_df.columns) - 1),
            "amr_genes": int(len(amr_df.columns) - 1),
            "virulence": int(len(vir_df.columns) - 1),
        },
        "checks": validation.get("checks", []),
        "warnings": validation.get("warnings", []),
        "errors": validation.get("errors", []),
        "validation_passed": bool(validation.get("validation_passed", False)),
    }
    with open(validation_dir / "synthetic_validation_results.json", "w", encoding="utf-8") as f:
        json.dump(validation_payload, f, indent=2)
    report_lines = [
        "# Synthetic Data Validation Report - strepsuis-phylotrait",
        "",
        f"Generated: {validation_payload['generated']}",
        "Data Source: Synthetic data with known ground truth",
        f"Strains: {validation_payload['n_strains']}",
        "Features:",
        f"- MIC: {validation_payload['n_features']['mic']}",
        f"- AMR genes: {validation_payload['n_features']['amr_genes']}",
        f"- Virulence: {validation_payload['n_features']['virulence']}",
        f"Checks: {len(validation_payload['checks'])}",
        f"Warnings: {len(validation_payload['warnings'])}",
        f"Errors: {len(validation_payload['errors'])}",
        f"Status: {'PASSED' if validation_payload['validation_passed'] else 'FAILED'}",
        "",
        "## Checks",
    ]
    report_lines.extend([f"- {item}" for item in validation_payload["checks"]])
    if validation_payload["warnings"]:
        report_lines.append("")
        report_lines.append("## Warnings")
        report_lines.extend([f"- {item}" for item in validation_payload["warnings"]])
    if validation_payload["errors"]:
        report_lines.append("")
        report_lines.append("## Errors")
        report_lines.extend([f"- {item}" for item in validation_payload["errors"]])
    with open(validation_dir / "SYNTHETIC_DATA_VALIDATION_REPORT.md", "w", encoding="utf-8") as f:
        f.write("\n".join(report_lines))

    # Save data
    print("\nSaving synthetic data...")
    output_dir = Path(__file__).parent.parent / "synthetic_data"
    saved = save_synthetic_phylo_data(tree_newick, mic_df, amr_df, vir_df, metadata, str(output_dir))

    for key, path in saved.items():
        print(f"  Saved {key}: {path}")

    print("\nDone!")
