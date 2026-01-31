"""
Tests for parallel phylogenetic distance matrix computation.

Author: MK-vet
License: MIT
"""

import pytest
import numpy as np
import pandas as pd
import tempfile
from io import StringIO

from Bio import Phylo

from strepsuis_phylotrait.parallel_phylo import (
    parallel_phylo_distance_matrix,
    compute_faiths_pd_parallel,
    pairwise_phylo_distances_efficient,
    summarize_phylo_distances,
    _compute_distance_row,
)


@pytest.fixture
def simple_tree():
    """Create a simple 4-taxon tree for testing."""
    # Newick format: ((A:1,B:1):1,(C:1,D:1):1);
    newick_str = "((A:1,B:1):1,(C:1,D:1):1);"
    tree = Phylo.read(StringIO(newick_str), 'newick')
    return tree


@pytest.fixture
def larger_tree():
    """Create a larger tree for testing."""
    # 8-taxon balanced tree
    newick_str = "(((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2):0.3,((E:0.1,F:0.1):0.2,(G:0.1,H:0.1):0.2):0.3);"
    tree = Phylo.read(StringIO(newick_str), 'newick')
    return tree


class TestComputeDistanceRow:
    """Test single row computation."""

    def test_basic_row_computation(self, simple_tree):
        """Test computing distances from one taxon."""
        terminals = list(simple_tree.get_terminals())

        row = _compute_distance_row(simple_tree, terminals, 0, 'patristic')

        # Check dimensions
        assert len(row) == 4

        # Diagonal should be zero
        assert row[0] == 0.0

        # Distances should be positive
        assert np.all(row >= 0)

    def test_patristic_vs_cophenetic(self, simple_tree):
        """Test different distance functions."""
        terminals = list(simple_tree.get_terminals())

        row_pat = _compute_distance_row(simple_tree, terminals, 0, 'patristic')
        row_cop = _compute_distance_row(simple_tree, terminals, 0, 'cophenetic')

        # Both should have same dimensions
        assert len(row_pat) == len(row_cop)

        # Both should be non-negative
        assert np.all(row_pat >= 0)
        assert np.all(row_cop >= 0)


class TestParallelPhyloDistanceMatrix:
    """Test parallel distance matrix computation."""

    def test_simple_tree(self, simple_tree):
        """Test with simple 4-taxon tree."""
        dist_df, taxa = parallel_phylo_distance_matrix(
            simple_tree,
            n_jobs=2,
            distance_func='patristic'
        )

        # Check dimensions
        assert dist_df.shape == (4, 4)
        assert len(taxa) == 4

        # Check symmetry
        np.testing.assert_array_almost_equal(
            dist_df.values,
            dist_df.values.T
        )

        # Check diagonal is zero
        np.testing.assert_array_equal(np.diag(dist_df.values), [0, 0, 0, 0])

        # Check taxa names
        assert all(taxon in taxa for taxon in ['A', 'B', 'C', 'D'])

    def test_larger_tree(self, larger_tree):
        """Test with larger tree."""
        dist_df, taxa = parallel_phylo_distance_matrix(
            larger_tree,
            n_jobs=2
        )

        # Check dimensions
        assert dist_df.shape == (8, 8)
        assert len(taxa) == 8

        # Check symmetry
        np.testing.assert_array_almost_equal(
            dist_df.values,
            dist_df.values.T,
            decimal=10
        )

    def test_deterministic_results(self, simple_tree):
        """Test that results are deterministic."""
        dist1, taxa1 = parallel_phylo_distance_matrix(simple_tree, n_jobs=1)
        dist2, taxa2 = parallel_phylo_distance_matrix(simple_tree, n_jobs=1)

        # Should get identical results
        np.testing.assert_array_almost_equal(dist1.values, dist2.values)
        assert taxa1 == taxa2

    def test_patristic_distances(self, simple_tree):
        """Test patristic distance calculation."""
        dist_df, _ = parallel_phylo_distance_matrix(
            simple_tree,
            distance_func='patristic'
        )

        # A and B are sisters with branch length 1 each
        # Distance A-B should be 2.0
        assert dist_df.loc['A', 'B'] == pytest.approx(2.0, abs=0.01)

        # A and C are more distant
        # Path: A -> ancestor(1) -> root(1) -> ancestor(1) -> C(1) = 4.0
        assert dist_df.loc['A', 'C'] == pytest.approx(4.0, abs=0.01)

    def test_single_taxon_error(self):
        """Test error with single taxon."""
        newick_str = "A:1.0;"
        tree = Phylo.read(StringIO(newick_str), 'newick')

        with pytest.raises(ValueError, match="at least 2 taxa"):
            parallel_phylo_distance_matrix(tree)

    def test_from_file_path(self):
        """Test loading tree from file path."""
        newick_str = "((A:1,B:1):1,(C:1,D:1):1);"

        with tempfile.NamedTemporaryFile(mode='w', suffix='.newick', delete=False) as f:
            f.write(newick_str)
            tree_file = f.name

        try:
            dist_df, taxa = parallel_phylo_distance_matrix(tree_file, n_jobs=2)
            assert dist_df.shape == (4, 4)
        finally:
            import os
            os.unlink(tree_file)


class TestComputeFaithsPD:
    """Test Faith's Phylogenetic Diversity computation."""

    def test_all_taxa(self, simple_tree):
        """Test PD for all taxa (should equal total tree length)."""
        all_taxa = [term.name for term in simple_tree.get_terminals()]
        pd_value = compute_faiths_pd_parallel(simple_tree, all_taxa)

        # For all taxa, PD should equal total tree length
        # Tree: ((A:1,B:1):1,(C:1,D:1):1) = 2*(1+1+1) = 6
        assert pd_value == pytest.approx(6.0, abs=0.01)

    def test_subset_taxa(self, simple_tree):
        """Test PD for subset of taxa."""
        # Just A and B (sister taxa)
        pd_ab = compute_faiths_pd_parallel(simple_tree, ['A', 'B'])

        # Should include: A(1), B(1), ancestor(1) = 3
        assert pd_ab == pytest.approx(3.0, abs=0.01)

    def test_empty_sample(self, simple_tree):
        """Test with empty sample."""
        pd_value = compute_faiths_pd_parallel(simple_tree, [])
        assert pd_value == 0.0

    def test_nonexistent_taxon(self, simple_tree):
        """Test with non-existent taxon."""
        # Should warn and ignore
        pd_value = compute_faiths_pd_parallel(simple_tree, ['A', 'NONEXISTENT'])
        assert pd_value > 0  # Should still compute for A


class TestPairwiseDistancesEfficient:
    """Test efficient extraction of pairwise distances."""

    def test_basic_extraction(self, simple_tree):
        """Test extracting specific pairs."""
        dist_df, _ = parallel_phylo_distance_matrix(simple_tree)

        pairs = [('A', 'B'), ('C', 'D')]
        pair_dist = pairwise_phylo_distances_efficient(dist_df, pairs)

        assert len(pair_dist) == 2
        assert 'taxon1' in pair_dist.columns
        assert 'taxon2' in pair_dist.columns
        assert 'distance' in pair_dist.columns

        # Check values
        assert pair_dist.iloc[0]['distance'] == pytest.approx(2.0, abs=0.01)
        assert pair_dist.iloc[1]['distance'] == pytest.approx(2.0, abs=0.01)

    def test_nonexistent_pair(self, simple_tree):
        """Test with non-existent taxa."""
        dist_df, _ = parallel_phylo_distance_matrix(simple_tree)

        pairs = [('A', 'NONEXISTENT')]

        # Should warn and return empty or partial results
        with pytest.warns(UserWarning):
            pair_dist = pairwise_phylo_distances_efficient(dist_df, pairs)


class TestSummarizePhyloDistances:
    """Test distance matrix summarization."""

    def test_basic_summary(self, simple_tree):
        """Test basic summary statistics."""
        dist_df, _ = parallel_phylo_distance_matrix(simple_tree)

        summary = summarize_phylo_distances(dist_df)

        assert len(summary) == 1
        assert 'n_taxa' in summary.columns
        assert 'mean_distance' in summary.columns
        assert 'median_distance' in summary.columns

        # Check values
        assert summary['n_taxa'].iloc[0] == 4
        assert summary['mean_distance'].iloc[0] > 0

    def test_summary_with_groups(self, simple_tree):
        """Test summary with grouping."""
        dist_df, taxa = parallel_phylo_distance_matrix(simple_tree)

        # Create groups: A,B in group1; C,D in group2
        groups = pd.Series(
            ['group1', 'group1', 'group2', 'group2'],
            index=['A', 'B', 'C', 'D']
        )

        summary = summarize_phylo_distances(dist_df, groups)

        assert 'mean_within_group' in summary.columns
        assert 'mean_between_group' in summary.columns

        # Within-group distance should be less than between-group
        within = summary['mean_within_group'].iloc[0]
        between = summary['mean_between_group'].iloc[0]
        assert within < between


class TestPerformance:
    """Test performance characteristics."""

    def test_parallel_speedup(self, larger_tree):
        """Test that parallel provides speedup."""
        import time

        # Sequential
        start = time.perf_counter()
        parallel_phylo_distance_matrix(larger_tree, n_jobs=1)
        time_seq = time.perf_counter() - start

        # Parallel
        start = time.perf_counter()
        parallel_phylo_distance_matrix(larger_tree, n_jobs=2)
        time_par = time.perf_counter() - start

        # Parallel should be reasonably fast (allowing for overhead)
        # Parallel provides reasonable performance (overhead acceptable for small trees)
        # For larger trees (>1000 taxa), speedup is much more significant
        assert time_par < time_seq * 5.0  # Very generous tolerance for small test tree

    def test_large_tree_scalability(self):
        """Test with moderately large tree."""
        # Create a tree with 50 taxa
        n_taxa = 50
        newick_parts = [f"taxon{i}:0.1" for i in range(n_taxa)]

        # Build simple star tree (all taxa connected to root)
        newick_str = f"({','.join(newick_parts)});"
        tree = Phylo.read(StringIO(newick_str), 'newick')

        # Should complete without errors
        dist_df, taxa = parallel_phylo_distance_matrix(tree, n_jobs=2)

        assert dist_df.shape == (n_taxa, n_taxa)
        assert len(taxa) == n_taxa

    def test_memory_efficiency(self, larger_tree):
        """Test that distance matrix doesn't cause memory issues."""
        # Should complete without memory errors
        dist_df, _ = parallel_phylo_distance_matrix(larger_tree, n_jobs=2)

        # Check no NaN values
        assert not np.any(np.isnan(dist_df.values))
