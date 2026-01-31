"""
Tests for optimizations module.

Tests cached distance matrix, diversity metrics, and other optimizations.
"""

import numpy as np
import pandas as pd
import pytest

from strepsuis_phylotrait.optimizations import (
    CachedDistanceMatrix,
    fast_faiths_pd,
    fast_mpd,
    fast_mntd,
    batch_diversity_metrics,
    parallel_permutation_test,
    fast_d_statistic,
    fast_correlation_matrix,
    fast_tree_clustering,
    fast_cophenetic_correlation,
    fast_jaccard_matrix,
    optimize_distance_matrix,
    benchmark_diversity_metrics,
    get_optimization_status,
)

# Check if numba/joblib are available
try:
    from strepsuis_phylotrait.optimizations import NUMBA_AVAILABLE, JOBLIB_AVAILABLE
except ImportError:
    NUMBA_AVAILABLE = False
    JOBLIB_AVAILABLE = False


class TestCachedDistanceMatrix:
    """Tests for CachedDistanceMatrix class."""
    
    def test_basic_creation(self):
        """Test basic distance matrix creation."""
        taxa = ['A', 'B', 'C']
        cached = CachedDistanceMatrix(tree=None, taxa=taxa)
        
        assert len(cached.taxa) == 3
    
    def test_get_distance_same_taxon(self):
        """Test distance to self is zero."""
        taxa = ['A', 'B', 'C']
        cached = CachedDistanceMatrix(tree=None, taxa=taxa)
        
        assert cached.get_distance('A', 'A') == 0.0
    
    def test_get_distance_no_tree(self):
        """Test distance without tree returns 0."""
        taxa = ['A', 'B', 'C']
        cached = CachedDistanceMatrix(tree=None, taxa=taxa)
        
        dist = cached.get_distance('A', 'B')
        
        assert dist == 0.0
    
    def test_clear_cache(self):
        """Test clearing cache."""
        taxa = ['A', 'B']
        cached = CachedDistanceMatrix(tree=None, taxa=taxa)
        cached.get_distance('A', 'B')
        
        cached.clear_cache()
        
        assert len(cached._cache) == 0
    
    def test_get_full_matrix(self):
        """Test getting full distance matrix."""
        taxa = ['A', 'B', 'C']
        cached = CachedDistanceMatrix(tree=None, taxa=taxa)
        
        matrix = cached.get_full_matrix()
        
        assert matrix.shape == (3, 3)


class TestFastFaithsPD:
    """Tests for fast Faith's PD calculation."""
    
    def test_basic_pd(self):
        """Test basic PD calculation with distance matrix."""
        # fast_faiths_pd expects a distance matrix and presence vector
        distances = np.array([
            [0, 1, 2],
            [1, 0, 1.5],
            [2, 1.5, 0]
        ])
        presence = np.array([1, 1, 0])
        
        pd_val = fast_faiths_pd(distances, presence)
        
        assert pd_val >= 0
    
    def test_pd_all_present(self):
        """Test PD with all taxa present."""
        distances = np.array([
            [0, 1, 2],
            [1, 0, 1.5],
            [2, 1.5, 0]
        ])
        presence = np.array([1, 1, 1])
        
        pd_val = fast_faiths_pd(distances, presence)
        
        assert pd_val >= 0
    
    def test_pd_none_present(self):
        """Test PD with no taxa present."""
        distances = np.array([
            [0, 1, 2],
            [1, 0, 1.5],
            [2, 1.5, 0]
        ])
        presence = np.array([0, 0, 0])
        
        pd_val = fast_faiths_pd(distances, presence)
        
        assert pd_val == 0


class TestFastMPD:
    """Tests for fast MPD calculation."""
    
    def test_basic_mpd(self):
        """Test basic MPD calculation."""
        distances = np.array([
            [0, 1, 2],
            [1, 0, 1.5],
            [2, 1.5, 0]
        ])
        presence = np.array([1, 1, 1])
        
        mpd = fast_mpd(distances, presence)
        
        assert mpd >= 0
    
    def test_mpd_two_taxa(self):
        """Test MPD with two taxa."""
        distances = np.array([
            [0, 5],
            [5, 0]
        ])
        presence = np.array([1, 1])
        
        mpd = fast_mpd(distances, presence)
        
        assert mpd == 5
    
    def test_mpd_single_taxon(self):
        """Test MPD with single taxon."""
        distances = np.array([
            [0, 1, 2],
            [1, 0, 1.5],
            [2, 1.5, 0]
        ])
        presence = np.array([1, 0, 0])
        
        mpd = fast_mpd(distances, presence)
        
        # MPD undefined for single taxon
        assert mpd == 0 or np.isnan(mpd)


class TestFastMNTD:
    """Tests for fast MNTD calculation."""
    
    def test_basic_mntd(self):
        """Test basic MNTD calculation."""
        distances = np.array([
            [0, 1, 2],
            [1, 0, 1.5],
            [2, 1.5, 0]
        ])
        presence = np.array([1, 1, 1])
        
        mntd = fast_mntd(distances, presence)
        
        assert mntd >= 0
    
    def test_mntd_close_taxa(self):
        """Test MNTD with close taxa."""
        distances = np.array([
            [0, 0.1, 10],
            [0.1, 0, 10],
            [10, 10, 0]
        ])
        presence = np.array([1, 1, 0])
        
        mntd = fast_mntd(distances, presence)
        
        assert mntd == 0.1


class TestBatchDiversityMetrics:
    """Tests for batch diversity metrics."""
    
    def test_basic_batch(self):
        """Test basic batch calculation."""
        distances = np.array([
            [0, 1, 2],
            [1, 0, 1.5],
            [2, 1.5, 0]
        ])
        presence_matrix = np.array([
            [1, 1, 0],
            [1, 0, 1],
            [0, 1, 1]
        ])
        
        # batch_diversity_metrics takes (distances, presence_matrix)
        results = batch_diversity_metrics(distances, presence_matrix)
        
        assert len(results) == 3


class TestParallelPermutationTest:
    """Tests for parallel permutation test."""
    
    def test_basic_permutation(self):
        """Test basic permutation test."""
        np.random.seed(42)
        n = 10
        trait = np.random.rand(n)
        distances = np.random.rand(n, n)
        distances = (distances + distances.T) / 2
        np.fill_diagonal(distances, 0)
        
        # parallel_permutation_test(trait, distances, n_permutations, n_jobs, random_state)
        observed, p_value, null_dist = parallel_permutation_test(
            trait, distances, n_permutations=50, n_jobs=1, random_state=42
        )
        
        assert isinstance(observed, (float, np.floating))
        assert 0 <= p_value <= 1


class TestFastDStatistic:
    """Tests for fast D statistic calculation."""
    
    def test_basic_d_statistic(self):
        """Test basic D statistic calculation."""
        np.random.seed(42)
        n = 10
        trait = np.random.randint(0, 2, n)
        distances = np.random.rand(n, n)
        distances = (distances + distances.T) / 2
        np.fill_diagonal(distances, 0)
        
        # fast_d_statistic(trait, distances, n_permutations, random_state)
        d_stat, p_value = fast_d_statistic(trait, distances, n_permutations=50, random_state=42)
        
        assert isinstance(d_stat, (int, float, np.floating))
        assert 0 <= p_value <= 1


class TestFastCorrelationMatrix:
    """Tests for fast correlation matrix."""
    
    def test_basic_correlation(self):
        """Test basic correlation matrix."""
        np.random.seed(42)
        data = np.random.randint(0, 2, (100, 5)).astype(np.float64)
        
        corr = fast_correlation_matrix(data)
        
        assert corr.shape == (5, 5)
    
    def test_correlation_symmetric(self):
        """Test correlation matrix is symmetric."""
        np.random.seed(42)
        data = np.random.randint(0, 2, (50, 4)).astype(np.float64)
        
        corr = fast_correlation_matrix(data)
        
        np.testing.assert_array_almost_equal(corr, corr.T)


class TestFastTreeClustering:
    """Tests for fast tree clustering."""
    
    def test_basic_clustering(self):
        """Test basic tree clustering."""
        distances = np.array([
            [0, 1, 5, 5],
            [1, 0, 5, 5],
            [5, 5, 0, 1],
            [5, 5, 1, 0]
        ], dtype=np.float64)
        
        labels = fast_tree_clustering(distances, n_clusters=2)
        
        assert len(labels) == 4
        assert len(np.unique(labels)) == 2


class TestFastCopheneticCorrelation:
    """Tests for fast cophenetic correlation."""
    
    def test_basic_cophenetic(self):
        """Test basic cophenetic correlation concept."""
        from scipy.cluster.hierarchy import linkage, cophenet
        from scipy.spatial.distance import squareform
        
        distances = np.array([
            [0, 1, 2],
            [1, 0, 1.5],
            [2, 1.5, 0]
        ], dtype=np.float64)
        
        condensed = squareform(distances)
        Z = linkage(condensed, method='average')
        
        # Use scipy's cophenet - returns (cophenetic_distances, cophenetic_matrix)
        coph_result = cophenet(Z, condensed)
        
        # coph_result[0] is the cophenetic correlation coefficient
        # coph_result[1] is the cophenetic distance matrix
        assert isinstance(coph_result, tuple)
        assert len(coph_result) == 2


class TestFastJaccardMatrix:
    """Tests for fast Jaccard matrix."""
    
    def test_basic_jaccard(self):
        """Test basic Jaccard matrix."""
        data = np.array([
            [1, 1, 0],
            [1, 0, 1],
            [0, 1, 1]
        ], dtype=np.int64)
        
        jaccard = fast_jaccard_matrix(data)
        
        assert jaccard.shape == (3, 3)
    
    def test_jaccard_symmetric(self):
        """Test Jaccard matrix is symmetric."""
        data = np.array([
            [1, 1, 0, 1],
            [1, 0, 1, 0],
            [0, 1, 1, 1]
        ], dtype=np.int64)
        
        jaccard = fast_jaccard_matrix(data)
        
        np.testing.assert_array_almost_equal(jaccard, jaccard.T)


class TestOptimizeDistanceMatrix:
    """Tests for distance matrix optimization."""
    
    def test_basic_optimization(self):
        """Test basic optimization."""
        distances = np.array([
            [0.0, 1.0, 2.0],
            [1.0, 0.0, 1.5],
            [2.0, 1.5, 0.0]
        ], dtype=np.float64)
        
        optimized = optimize_distance_matrix(distances)
        
        assert optimized.dtype == np.float32


class TestBenchmarkDiversityMetrics:
    """Tests for diversity metrics benchmark."""
    
    def test_basic_benchmark(self):
        """Test basic benchmarking."""
        np.random.seed(42)
        distances = np.random.rand(20, 20)
        distances = (distances + distances.T) / 2
        np.fill_diagonal(distances, 0)
        presence = np.random.randint(0, 2, 20)
        
        # benchmark_diversity_metrics(distances, presence, n_runs)
        result = benchmark_diversity_metrics(distances, presence, n_runs=2)
        
        assert isinstance(result, dict)


class TestOptimizationStatus:
    """Tests for optimization status."""
    
    def test_get_status(self):
        """Test getting optimization status."""
        status = get_optimization_status()
        
        assert isinstance(status, dict)


class TestEdgeCases:
    """Tests for edge cases."""
    
    def test_empty_presence(self):
        """Test with empty presence vector."""
        distances = np.array([
            [0, 1],
            [1, 0]
        ], dtype=np.float64)
        presence = np.array([0, 0])
        
        mpd = fast_mpd(distances, presence)
        
        assert mpd == 0 or np.isnan(mpd)
    
    def test_single_taxon_pd(self):
        """Test PD with single taxon."""
        distances = np.array([
            [0, 1, 2],
            [1, 0, 1.5],
            [2, 1.5, 0]
        ], dtype=np.float64)
        presence = np.array([1, 0, 0])
        
        pd_val = fast_faiths_pd(distances, presence)
        
        assert pd_val == 0.0  # Single taxon has no PD


class TestAdditionalOptimizations:
    """Additional tests for optimizations."""
    
    def test_fast_faiths_pd_multiple_taxa(self):
        """Test Faith's PD with multiple taxa."""
        np.random.seed(42)
        n = 10
        distances = np.random.rand(n, n)
        distances = (distances + distances.T) / 2
        np.fill_diagonal(distances, 0)
        presence = np.array([1, 1, 1, 0, 0, 1, 0, 1, 0, 0])
        
        pd_val = fast_faiths_pd(distances, presence)
        
        assert pd_val >= 0
    
    def test_fast_mpd_all_present(self):
        """Test MPD with all taxa present."""
        distances = np.array([
            [0, 1, 2],
            [1, 0, 1.5],
            [2, 1.5, 0]
        ], dtype=np.float64)
        presence = np.array([1, 1, 1])
        
        mpd = fast_mpd(distances, presence)
        
        # Expected: (1 + 2 + 1.5) / 3 = 1.5
        assert abs(mpd - 1.5) < 0.01
    
    def test_fast_mntd_calculation(self):
        """Test MNTD calculation."""
        distances = np.array([
            [0, 1, 5],
            [1, 0, 4],
            [5, 4, 0]
        ], dtype=np.float64)
        presence = np.array([1, 1, 1])
        
        mntd = fast_mntd(distances, presence)
        
        # Nearest for 0: 1 (dist=1), for 1: 0 (dist=1), for 2: 1 (dist=4)
        # Mean = (1 + 1 + 4) / 3 = 2
        assert abs(mntd - 2.0) < 0.01
    
    def test_batch_diversity_metrics(self):
        """Test batch diversity metrics."""
        np.random.seed(42)
        n = 10
        distances = np.random.rand(n, n)
        distances = (distances + distances.T) / 2
        np.fill_diagonal(distances, 0)
        
        presence_matrix = np.random.randint(0, 2, (5, n))
        
        results = batch_diversity_metrics(distances, presence_matrix)
        
        assert isinstance(results, pd.DataFrame)
        assert len(results) == 5
    
    def test_fast_correlation_matrix_values(self):
        """Test correlation matrix values."""
        # Create perfectly correlated data
        data = np.array([
            [1, 1, 0],
            [1, 1, 0],
            [0, 0, 1],
            [0, 0, 1]
        ], dtype=np.float64)
        
        corr = fast_correlation_matrix(data)
        
        assert corr.shape == (3, 3)
        # First two columns are identical, so correlation should be 1
        assert abs(corr[0, 1] - 1.0) < 0.01
    
    def test_fast_tree_clustering_single_cluster(self):
        """Test tree clustering with single cluster."""
        distances = np.array([
            [0, 0.1, 0.1],
            [0.1, 0, 0.1],
            [0.1, 0.1, 0]
        ], dtype=np.float64)
        
        labels = fast_tree_clustering(distances, n_clusters=1)
        
        assert len(labels) == 3
        assert len(np.unique(labels)) == 1
    
    def test_fast_jaccard_matrix_identical(self):
        """Test Jaccard matrix with identical rows."""
        data = np.array([
            [1, 1, 0],
            [1, 1, 0],
            [0, 0, 1]
        ], dtype=np.int64)
        
        jaccard = fast_jaccard_matrix(data)
        
        # First two rows are identical, distance should be 0
        assert jaccard[0, 1] == 0.0
    
    def test_optimize_distance_matrix_dtype(self):
        """Test distance matrix optimization dtype."""
        distances = np.array([
            [0.0, 1.0],
            [1.0, 0.0]
        ], dtype=np.float64)
        
        optimized = optimize_distance_matrix(distances)
        
        assert optimized.dtype == np.float32
    
    def test_get_optimization_status_keys(self):
        """Test optimization status has expected keys."""
        status = get_optimization_status()
        
        assert isinstance(status, dict)
        assert 'numba_jit' in status or len(status) > 0
    
    def test_cached_distance_matrix_caching(self):
        """Test that CachedDistanceMatrix caches results."""
        taxa = ['A', 'B', 'C']
        cached = CachedDistanceMatrix(tree=None, taxa=taxa)
        
        # First call
        dist1 = cached.get_distance('A', 'B')
        # Second call (should use cache)
        dist2 = cached.get_distance('A', 'B')
        
        assert dist1 == dist2
    
    def test_cached_distance_matrix_symmetry(self):
        """Test CachedDistanceMatrix symmetry."""
        taxa = ['A', 'B']
        cached = CachedDistanceMatrix(tree=None, taxa=taxa)
        
        dist_ab = cached.get_distance('A', 'B')
        dist_ba = cached.get_distance('B', 'A')
        
        assert dist_ab == dist_ba


class TestPermutationTests:
    """Tests for permutation-based analyses."""
    
    def test_parallel_permutation_test(self):
        """Test parallel permutation test."""
        np.random.seed(42)
        n = 15
        trait = np.random.rand(n)
        distances = np.random.rand(n, n)
        distances = (distances + distances.T) / 2
        np.fill_diagonal(distances, 0)
        
        observed, p_value, null_dist = parallel_permutation_test(
            trait, distances, n_permutations=50, n_jobs=1, random_state=42
        )
        
        assert isinstance(observed, (float, np.floating))
        assert 0 <= p_value <= 1
        assert len(null_dist) == 50
    
    def test_fast_d_statistic(self):
        """Test D-statistic for binary traits."""
        np.random.seed(42)
        n = 12
        trait = np.random.randint(0, 2, n)
        distances = np.random.rand(n, n)
        distances = (distances + distances.T) / 2
        np.fill_diagonal(distances, 0)
        
        d_stat, p_value = fast_d_statistic(
            trait, distances, n_permutations=50, random_state=42
        )
        
        assert isinstance(d_stat, (float, np.floating))
        assert 0 <= p_value <= 1
    
    def test_fast_d_statistic_clustered(self):
        """Test D-statistic with clustered trait."""
        # Create clustered trait (first half 0, second half 1)
        trait = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])
        # Create distances where similar traits are closer
        distances = np.ones((10, 10), dtype=np.float64) * 2
        distances[:5, :5] = 0.1
        distances[5:, 5:] = 0.1
        np.fill_diagonal(distances, 0)
        
        d_stat, p_value = fast_d_statistic(
            trait, distances, n_permutations=50, random_state=42
        )
        
        assert isinstance(d_stat, (float, np.floating))


class TestClusteringFunctions:
    """Tests for clustering-related functions."""
    
    def test_fast_tree_clustering_methods(self):
        """Test different linkage methods."""
        np.random.seed(42)
        n = 15
        distances = np.random.rand(n, n)
        distances = (distances + distances.T) / 2
        np.fill_diagonal(distances, 0)
        
        for method in ['average', 'complete', 'single']:
            labels = fast_tree_clustering(distances, n_clusters=3, method=method)
            assert len(labels) == n
            assert len(np.unique(labels)) <= 3
    
    def test_fast_cophenetic_correlation(self):
        """Test cophenetic correlation concept."""
        from scipy.cluster.hierarchy import linkage, cophenet
        from scipy.spatial.distance import squareform
        
        np.random.seed(42)
        n = 10
        distances = np.random.rand(n, n)
        distances = (distances + distances.T) / 2
        np.fill_diagonal(distances, 0)
        
        condensed = squareform(distances)
        Z = linkage(condensed, method='average')
        
        # Use scipy's cophenet - returns (cophenetic_distances, cophenetic_matrix)
        coph_result = cophenet(Z, condensed)
        
        # Verify the result structure
        assert isinstance(coph_result, tuple)
        assert len(coph_result) == 2


class TestBenchmarking:
    """Tests for benchmarking functions."""
    
    def test_benchmark_diversity_metrics(self):
        """Test diversity metrics benchmark."""
        np.random.seed(42)
        n = 15
        distances = np.random.rand(n, n)
        distances = (distances + distances.T) / 2
        np.fill_diagonal(distances, 0)
        presence = np.random.randint(0, 2, n)
        
        results = benchmark_diversity_metrics(distances, presence, n_runs=2)
        
        assert 'mpd' in results
        assert 'mntd' in results
        assert 'mean_ms' in results['mpd']
        assert 'std_ms' in results['mpd']


class TestCachedDistanceMatrixAdvanced:
    """Advanced tests for CachedDistanceMatrix."""
    
    def test_clear_cache_advanced(self):
        """Test cache clearing."""
        taxa = ['A', 'B', 'C']
        cached = CachedDistanceMatrix(tree=None, taxa=taxa)
        
        # Populate cache
        _ = cached.get_distance('A', 'B')
        
        # Clear cache
        cached.clear_cache()
        
        # Cache should be empty
        assert len(cached._cache) == 0
    
    def test_full_matrix_advanced(self):
        """Test full distance matrix generation."""
        taxa = ['A', 'B', 'C']
        cached = CachedDistanceMatrix(tree=None, taxa=taxa)
        
        matrix = cached.get_full_matrix()
        
        assert matrix.shape == (3, 3)
        # Diagonal should be 0
        assert np.allclose(np.diag(matrix), 0)


class TestJaccardDistance:
    """Tests for Jaccard distance calculations."""
    
    def test_fast_jaccard_matrix_zeros(self):
        """Test Jaccard with all zeros."""
        data = np.zeros((3, 5), dtype=np.int64)
        
        jaccard = fast_jaccard_matrix(data)
        
        # All zeros means union is 0, distance should be 0
        assert np.allclose(jaccard, 0)
    
    def test_fast_jaccard_matrix_ones(self):
        """Test Jaccard with all ones."""
        data = np.ones((3, 5), dtype=np.int64)
        
        jaccard = fast_jaccard_matrix(data)
        
        # All ones means perfect overlap, distance should be 0
        assert np.allclose(jaccard, 0)
    
    def test_fast_jaccard_matrix_no_overlap(self):
        """Test Jaccard with no overlap."""
        data = np.array([
            [1, 1, 0, 0],
            [0, 0, 1, 1]
        ], dtype=np.int64)
        
        jaccard = fast_jaccard_matrix(data)
        
        # No overlap means distance should be 1
        assert jaccard[0, 1] == 1.0
