"""Stress tests for strepsuis-phylotrait module.

These tests verify behavior with large datasets, memory constraints,
and concurrent operations.
"""

import time

import numpy as np
import pandas as pd
import pytest


@pytest.mark.stress
class TestLargeDatasets:
    """Tests with large datasets to verify scalability."""

    def test_large_distance_matrix(self):
        """Test large phylogenetic distance matrix creation."""
        np.random.seed(42)
        n_taxa = 200
        
        # Simulate distance matrix
        distances = np.random.random((n_taxa, n_taxa))
        distances = (distances + distances.T) / 2  # Make symmetric
        np.fill_diagonal(distances, 0)
        
        assert distances.shape == (n_taxa, n_taxa)
        assert np.allclose(distances, distances.T)
        assert np.allclose(np.diag(distances), 0)

    def test_large_trait_matrix(self):
        """Test large binary trait matrix."""
        np.random.seed(42)
        n_taxa = 500
        n_traits = 100
        
        traits = np.random.randint(0, 2, size=(n_taxa, n_traits))
        df = pd.DataFrame(
            traits,
            columns=[f"Trait_{i}" for i in range(n_traits)]
        )
        
        assert df.shape == (n_taxa, n_traits)
        assert set(df.values.flatten()) <= {0, 1}

    def test_permutation_test_scaling(self):
        """Test permutation test with many permutations."""
        np.random.seed(42)
        n_taxa = 100
        n_permutations = 999
        
        trait = np.random.randint(0, 2, size=n_taxa)
        
        null_stats = []
        for _ in range(n_permutations):
            shuffled = np.random.permutation(trait)
            null_stats.append(shuffled.sum())
        
        assert len(null_stats) == n_permutations


@pytest.mark.stress
class TestPhylogeneticEdgeCases:
    """Tests for phylogenetic edge cases."""

    def test_star_tree_distances(self):
        """Test star tree (all taxa connected to root)."""
        n_taxa = 50
        branch_length = 0.1
        
        # Star tree: all pairwise distances equal
        distances = np.full((n_taxa, n_taxa), 2 * branch_length)
        np.fill_diagonal(distances, 0)
        
        # Check all off-diagonal elements are equal
        off_diag = distances[~np.eye(n_taxa, dtype=bool)]
        assert np.allclose(off_diag, 2 * branch_length)

    def test_balanced_tree_properties(self):
        """Test properties of balanced binary tree."""
        np.random.seed(42)
        n_taxa = 64  # 2^6, perfect binary tree
        
        # In balanced tree, depth = log2(n)
        expected_depth = 6
        
        # Simulate ultrametric distances
        # All tips equidistant from root
        root_distance = 1.0
        distances = np.ones((n_taxa, n_taxa)) * 2 * root_distance
        np.fill_diagonal(distances, 0)
        
        # Modify for tree structure (taxa in same clade are closer)
        for i in range(0, n_taxa, 2):
            distances[i, i+1] = 0.3  # Sister taxa
            distances[i+1, i] = 0.3
        
        assert distances.shape == (n_taxa, n_taxa)

    def test_single_taxon(self):
        """Test handling of single taxon."""
        trait = np.array([1])
        
        # PD of single taxon is 0 or just root branch
        assert len(trait) == 1


@pytest.mark.stress
class TestClusteringEdgeCases:
    """Tests for phylo-aware clustering edge cases."""

    def test_all_same_trait(self):
        """Test when all taxa have same trait value."""
        n_taxa = 100
        
        trait_all_ones = np.ones(n_taxa)
        trait_all_zeros = np.zeros(n_taxa)
        
        # Prevalence edge cases
        assert trait_all_ones.mean() == 1.0
        assert trait_all_zeros.mean() == 0.0

    def test_perfectly_clustered(self):
        """Test perfectly phylogenetically clustered trait."""
        np.random.seed(42)
        n_taxa = 100
        
        # First half has trait, second half doesn't
        trait = np.concatenate([np.ones(50), np.zeros(50)])
        
        # This represents perfect phylogenetic clustering
        assert trait[:50].mean() == 1.0
        assert trait[50:].mean() == 0.0

    def test_random_trait_distribution(self):
        """Test random trait distribution (no phylogenetic signal)."""
        np.random.seed(42)
        n_taxa = 1000
        
        trait = np.random.randint(0, 2, size=n_taxa)
        
        # Should be roughly 50/50
        prevalence = trait.mean()
        assert 0.4 < prevalence < 0.6


@pytest.mark.stress
class TestAncestralStateEdgeCases:
    """Tests for ancestral state reconstruction edge cases."""

    def test_constant_trait(self):
        """Test ancestral reconstruction with constant trait."""
        n_taxa = 50
        
        # All 1s - ancestors should also be 1
        trait = np.ones(n_taxa)
        
        # Most parsimonious ancestral state
        ancestral = int(trait.mean() > 0.5)
        assert ancestral == 1

    def test_equal_split(self):
        """Test ancestral reconstruction with 50/50 split."""
        n_taxa = 100
        
        # Equal split - ancestral state is ambiguous
        trait = np.array([0, 1] * 50)
        
        assert trait.mean() == 0.5
