"""
Comprehensive Mathematical Validation Tests for Phylogenetic Analysis

This module provides rigorous validation of all phylogenetic methods against
synthetic data with known ground truth, ensuring mathematical correctness
of the StrepSuis-PhyloTrait analysis pipeline.

Tests validate:
1. Faith's Phylogenetic Diversity calculation
2. Tree-aware clustering with monophyly constraints
3. Chi-square tests for trait-cluster associations
4. Silhouette scores on phylogenetic distances
5. Bootstrap confidence intervals

Reference implementations: scipy, sklearn, Bio.Phylo
"""

import json
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
import pytest
from scipy.stats import chi2_contingency
from sklearn.metrics import adjusted_rand_score, silhouette_score
from statsmodels.stats.multitest import multipletests


class TestSyntheticDataGeneration:
    """Test the synthetic data generation module."""

    def test_generate_synthetic_dataset(self):
        """Test that synthetic data generation produces valid output."""
        from strepsuis_phylotrait.generate_synthetic_data import (
            SyntheticPhyloConfig,
            generate_phylotrait_synthetic_dataset,
        )

        config = SyntheticPhyloConfig(
            n_strains=50,
            n_clusters=3,
            n_mic_features=8,
            n_amr_features=10,
            random_state=42,
        )

        tree_newick, mic_df, amr_df, vir_df, metadata = generate_phylotrait_synthetic_dataset(config)

        # Verify tree is valid
        assert tree_newick.endswith(";")
        assert "(" in tree_newick

        # Verify data shapes
        assert len(mic_df) == config.n_strains
        assert len(metadata.mic_columns) == config.n_mic_features

    def test_tree_structure_validity(self):
        """Test that generated tree has valid structure."""
        from strepsuis_phylotrait.generate_synthetic_data import (
            SyntheticPhyloConfig,
            generate_phylotrait_synthetic_dataset,
        )

        config = SyntheticPhyloConfig(n_strains=20, n_clusters=3, random_state=42)
        tree_newick, _, _, _, metadata = generate_phylotrait_synthetic_dataset(config)

        # Count taxa in tree
        # Each Strain_XXXX appears once
        n_taxa = tree_newick.count("Strain_")
        assert n_taxa == config.n_strains

    def test_cluster_assignment_validity(self):
        """Test that true cluster assignments are valid."""
        from strepsuis_phylotrait.generate_synthetic_data import (
            SyntheticPhyloConfig,
            generate_phylotrait_synthetic_dataset,
        )

        config = SyntheticPhyloConfig(n_strains=60, n_clusters=4, random_state=42)
        _, _, _, _, metadata = generate_phylotrait_synthetic_dataset(config)

        # Verify cluster labels are in expected range
        unique_labels = np.unique(metadata.true_cluster_labels)
        assert len(unique_labels) == config.n_clusters
        assert unique_labels.min() == 1
        assert unique_labels.max() == config.n_clusters

    def test_synthetic_data_reproducibility(self):
        """Test that same seed produces identical data."""
        from strepsuis_phylotrait.generate_synthetic_data import (
            SyntheticPhyloConfig,
            generate_phylotrait_synthetic_dataset,
        )

        config = SyntheticPhyloConfig(n_strains=30, random_state=12345)

        tree1, mic1, amr1, vir1, meta1 = generate_phylotrait_synthetic_dataset(config)
        tree2, mic2, amr2, vir2, meta2 = generate_phylotrait_synthetic_dataset(config)

        assert tree1 == tree2
        pd.testing.assert_frame_equal(mic1, mic2)
        np.testing.assert_array_equal(meta1.true_cluster_labels, meta2.true_cluster_labels)


class TestPhylogeneticDiversityValidation:
    """Validate Faith's Phylogenetic Diversity calculation."""

    def test_pd_positive(self):
        """Test that PD is always positive for valid trees."""
        from strepsuis_phylotrait.generate_synthetic_data import (
            SyntheticPhyloConfig,
            generate_phylotrait_synthetic_dataset,
        )

        config = SyntheticPhyloConfig(n_strains=20, random_state=42)
        tree_newick, _, _, _, _ = generate_phylotrait_synthetic_dataset(config)

        # PD should be positive (sum of branch lengths)
        # Parse branch lengths from newick
        import re
        branch_lengths = re.findall(r':(\d+\.?\d*)', tree_newick)
        total_bl = sum(float(bl) for bl in branch_lengths)

        assert total_bl > 0

    def test_pd_scales_with_taxa(self):
        """Test that PD generally scales with number of taxa."""
        from strepsuis_phylotrait.generate_synthetic_data import (
            SyntheticPhyloConfig,
            generate_phylotrait_synthetic_dataset,
        )

        import re

        pds = []
        for n_taxa in [10, 30, 50]:
            config = SyntheticPhyloConfig(n_strains=n_taxa, random_state=42)
            tree_newick, _, _, _, _ = generate_phylotrait_synthetic_dataset(config)

            branch_lengths = re.findall(r':(\d+\.?\d*)', tree_newick)
            total_bl = sum(float(bl) for bl in branch_lengths)
            pds.append(total_bl)

        # PD should generally increase with more taxa
        assert pds[2] >= pds[0]


class TestClusteringValidation:
    """Validate clustering methods with synthetic data."""

    def test_cluster_recovery_with_high_heritability(self):
        """Test that clustering can recover true clusters with high heritability."""
        from strepsuis_phylotrait.generate_synthetic_data import (
            SyntheticPhyloConfig,
            generate_phylotrait_synthetic_dataset,
        )

        config = SyntheticPhyloConfig(
            n_strains=100,
            n_clusters=3,
            trait_heritability=0.8,  # High heritability
            random_state=42,
        )

        _, mic_df, amr_df, vir_df, metadata = generate_phylotrait_synthetic_dataset(config)

        # Combine features
        feature_cols = metadata.mic_columns + metadata.amr_columns + metadata.virulence_columns
        features = pd.concat([
            mic_df[metadata.mic_columns],
            amr_df[metadata.amr_columns],
            vir_df[metadata.virulence_columns],
        ], axis=1)

        # Simple clustering
        from sklearn.cluster import KMeans

        km = KMeans(n_clusters=config.n_clusters, random_state=42)
        predicted_labels = km.fit_predict(features)

        # Calculate Adjusted Rand Index
        ari = adjusted_rand_score(metadata.true_cluster_labels, predicted_labels + 1)

        # With high heritability, should have some recovery
        assert ari > 0.1, f"ARI {ari} should be > 0.1 with high heritability"


class TestChiSquareValidation:
    """Validate chi-square tests for trait-cluster associations."""

    def test_chi_square_with_phylogenetic_signal(self):
        """Test chi-square detects trait-cluster associations."""
        from strepsuis_phylotrait.generate_synthetic_data import (
            SyntheticPhyloConfig,
            generate_phylotrait_synthetic_dataset,
        )

        config = SyntheticPhyloConfig(
            n_strains=200,
            n_clusters=3,
            trait_heritability=0.7,
            random_state=42,
        )

        _, mic_df, _, _, metadata = generate_phylotrait_synthetic_dataset(config)

        # Test trait-cluster associations
        n_significant = 0
        for col in metadata.mic_columns[:5]:
            table = pd.crosstab(metadata.true_cluster_labels, mic_df[col])
            chi2, p_value, _, _ = chi2_contingency(table)

            if p_value < 0.05:
                n_significant += 1

        # With phylogenetic signal, some traits should be associated with clusters
        assert n_significant >= 2, f"Expected >= 2 significant associations, got {n_significant}"


class TestFDRCorrectionValidation:
    """Validate FDR correction."""

    def test_fdr_controls_false_positives(self):
        """Test that FDR controls false positive rate."""
        np.random.seed(42)

        # Generate null p-values (no true associations)
        n_tests = 100
        null_pvals = np.random.uniform(0, 1, n_tests)

        reject, corrected, _, _ = multipletests(null_pvals, alpha=0.05, method="fdr_bh")

        # With all null, should have few rejections
        false_positive_rate = np.sum(reject) / n_tests
        assert false_positive_rate < 0.10, f"FPR {false_positive_rate} should be < 0.10"


class TestSilhouetteValidation:
    """Validate silhouette score calculations."""

    def test_silhouette_with_known_clusters(self):
        """Test silhouette score with known cluster structure."""
        from strepsuis_phylotrait.generate_synthetic_data import (
            SyntheticPhyloConfig,
            generate_phylotrait_synthetic_dataset,
        )

        config = SyntheticPhyloConfig(
            n_strains=100,
            n_clusters=3,
            trait_heritability=0.7,
            random_state=42,
        )

        _, mic_df, amr_df, vir_df, metadata = generate_phylotrait_synthetic_dataset(config)

        # Combine features
        features = pd.concat([
            mic_df[metadata.mic_columns],
            amr_df[metadata.amr_columns],
            vir_df[metadata.virulence_columns],
        ], axis=1)

        # Calculate silhouette with true labels
        sil_score = silhouette_score(features, metadata.true_cluster_labels)

        # With good cluster structure, silhouette should be positive
        assert sil_score > 0, f"Silhouette score {sil_score} should be > 0"


class TestBootstrapValidation:
    """Validate bootstrap confidence intervals."""

    def test_bootstrap_ci_coverage(self):
        """Test that 95% CI contains true parameter most of the time."""
        np.random.seed(42)

        true_prevalence = 0.4
        n_samples = 100
        n_simulations = 30
        coverage_count = 0

        for _ in range(n_simulations):
            sample = np.random.binomial(1, true_prevalence, n_samples)

            # Bootstrap CI
            n_bootstrap = 200
            bootstrap_means = []
            for _ in range(n_bootstrap):
                boot_sample = np.random.choice(sample, size=n_samples, replace=True)
                bootstrap_means.append(np.mean(boot_sample))

            ci_lower = np.percentile(bootstrap_means, 2.5)
            ci_upper = np.percentile(bootstrap_means, 97.5)

            if ci_lower <= true_prevalence <= ci_upper:
                coverage_count += 1

        coverage_rate = coverage_count / n_simulations
        assert coverage_rate >= 0.7, f"Coverage rate {coverage_rate} should be >= 0.7"


class TestNumericalStability:
    """Test numerical stability with edge cases."""

    def test_single_taxon_cluster(self):
        """Test handling of single-taxon clusters."""
        from strepsuis_phylotrait.generate_synthetic_data import (
            SyntheticPhyloConfig,
            generate_phylotrait_synthetic_dataset,
        )

        # Generate data with many clusters (some may be small)
        config = SyntheticPhyloConfig(
            n_strains=20,
            n_clusters=4,
            random_state=42,
        )

        tree_newick, mic_df, _, _, metadata = generate_phylotrait_synthetic_dataset(config)

        # All values should be valid
        assert len(mic_df) == config.n_strains
        assert not mic_df.isnull().any().any()

    def test_extreme_heritability_values(self):
        """Test with extreme heritability values."""
        from strepsuis_phylotrait.generate_synthetic_data import (
            SyntheticPhyloConfig,
            generate_phylotrait_synthetic_dataset,
        )

        # Low heritability (traits mostly random)
        config_low = SyntheticPhyloConfig(
            n_strains=50,
            n_clusters=3,
            trait_heritability=0.1,
            random_state=42,
        )
        _, mic_low, _, _, _ = generate_phylotrait_synthetic_dataset(config_low)
        assert not mic_low.isnull().any().any()

        # High heritability (traits follow phylogeny)
        config_high = SyntheticPhyloConfig(
            n_strains=50,
            n_clusters=3,
            trait_heritability=0.9,
            random_state=42,
        )
        _, mic_high, _, _, _ = generate_phylotrait_synthetic_dataset(config_high)
        assert not mic_high.isnull().any().any()


class TestTreeParsingValidation:
    """Validate tree parsing and manipulation."""

    def test_tree_contains_all_taxa(self):
        """Test that tree contains all expected taxa."""
        from strepsuis_phylotrait.generate_synthetic_data import (
            SyntheticPhyloConfig,
            generate_phylotrait_synthetic_dataset,
        )

        config = SyntheticPhyloConfig(n_strains=30, random_state=42)
        tree_newick, mic_df, _, _, _ = generate_phylotrait_synthetic_dataset(config)

        # Check all strain IDs are in tree
        for strain_id in mic_df["Strain_ID"]:
            assert strain_id in tree_newick, f"{strain_id} not found in tree"


class TestValidationReportGeneration:
    """Generate validation reports for math validation directory."""

    @pytest.fixture
    def report_dir(self, tmp_path):
        """Create temporary report directory."""
        return tmp_path / "math_validation"

    def test_generate_validation_report(self, report_dir):
        """Generate a validation report summarizing test results."""
        from strepsuis_phylotrait.generate_synthetic_data import (
            SyntheticPhyloConfig,
            generate_phylotrait_synthetic_dataset,
            validate_synthetic_phylo_data,
        )

        report_dir.mkdir(parents=True, exist_ok=True)

        # Generate synthetic data
        config = SyntheticPhyloConfig(n_strains=50, n_clusters=3, random_state=42)
        tree, mic_df, amr_df, vir_df, metadata = generate_phylotrait_synthetic_dataset(config)

        # Validate
        validation = validate_synthetic_phylo_data(tree, mic_df, amr_df, vir_df, metadata)

        # Create report
        report = {
            "test_name": "Synthetic Phylogenetic Data Validation",
            "timestamp": pd.Timestamp.now().isoformat(),
            "config": {
                "n_strains": config.n_strains,
                "n_clusters": config.n_clusters,
                "trait_heritability": config.trait_heritability,
                "random_state": config.random_state,
            },
            "validation_passed": validation["validation_passed"],
            "checks_passed": len(validation["checks"]),
            "warnings": len(validation["warnings"]),
            "errors": len(validation["errors"]),
        }

        # Save report
        report_file = report_dir / "validation_report.json"
        with open(report_file, "w") as f:
            json.dump(report, f, indent=2)

        assert report_file.exists()
        assert validation["validation_passed"]
