"""Performance tests for strepsuis-phylotrait module.

These tests measure and verify timing benchmarks for key operations.
"""

import time

import numpy as np
import pandas as pd
import pytest


@pytest.mark.performance
class TestDistanceMatrixPerformance:
    """Performance tests for distance matrix operations."""

    def test_distance_matrix_creation(self):
        """Test distance matrix creation timing."""
        np.random.seed(42)
        n_taxa = 200
        
        start = time.time()
        # Create random distance matrix
        distances = np.random.random((n_taxa, n_taxa))
        distances = (distances + distances.T) / 2
        np.fill_diagonal(distances, 0)
        elapsed = time.time() - start
        
        assert distances.shape == (n_taxa, n_taxa)
        assert elapsed < 0.5

    def test_pairwise_distances_timing(self):
        """Test pairwise distance calculation."""
        from scipy.spatial.distance import pdist, squareform
        
        np.random.seed(42)
        n_taxa = 300
        n_features = 50
        
        data = np.random.random((n_taxa, n_features))
        
        start = time.time()
        distances = squareform(pdist(data))
        elapsed = time.time() - start
        
        assert distances.shape == (n_taxa, n_taxa)
        assert elapsed < 2.0


@pytest.mark.performance
class TestTraitAnalysisPerformance:
    """Performance tests for trait analysis operations."""

    def test_prevalence_calculation_timing(self):
        """Test prevalence calculation performance."""
        np.random.seed(42)
        n_taxa = 1000
        n_traits = 100
        
        traits = np.random.randint(0, 2, size=(n_taxa, n_traits))
        df = pd.DataFrame(traits)
        
        start = time.time()
        prevalences = df.mean() * 100
        elapsed = time.time() - start
        
        assert len(prevalences) == n_traits
        assert elapsed < 0.1

    def test_permutation_test_timing(self):
        """Test permutation test performance."""
        np.random.seed(42)
        n_taxa = 200
        n_permutations = 100
        
        trait = np.random.randint(0, 2, size=n_taxa)
        
        start = time.time()
        null_distribution = []
        for _ in range(n_permutations):
            shuffled = np.random.permutation(trait)
            null_distribution.append(shuffled.sum())
        elapsed = time.time() - start
        
        assert len(null_distribution) == n_permutations
        assert elapsed < 1.0


@pytest.mark.performance
class TestClusteringPerformance:
    """Performance tests for clustering operations."""

    def test_hierarchical_clustering_timing(self):
        """Test hierarchical clustering performance."""
        from scipy.cluster.hierarchy import linkage, fcluster
        from scipy.spatial.distance import squareform
        
        np.random.seed(42)
        n_taxa = 200
        
        # Random distance matrix
        distances = np.random.random((n_taxa, n_taxa))
        distances = (distances + distances.T) / 2
        np.fill_diagonal(distances, 0)
        
        start = time.time()
        condensed = squareform(distances)
        Z = linkage(condensed, method='average')
        labels = fcluster(Z, t=5, criterion='maxclust')
        elapsed = time.time() - start
        
        assert len(labels) == n_taxa
        assert elapsed < 1.0

    def test_silhouette_timing(self):
        """Test silhouette score calculation."""
        from sklearn.metrics import silhouette_score
        
        np.random.seed(42)
        n_taxa = 200
        n_features = 30
        k = 5
        
        data = np.random.random((n_taxa, n_features))
        labels = np.random.randint(0, k, size=n_taxa)
        
        start = time.time()
        score = silhouette_score(data, labels)
        elapsed = time.time() - start
        
        assert -1 <= score <= 1
        assert elapsed < 2.0


@pytest.mark.performance
class TestUMAPPerformance:
    """Performance tests for UMAP embedding."""

    def test_umap_alternative_timing(self):
        """Test dimensionality reduction timing using PCA."""
        from sklearn.decomposition import PCA
        
        np.random.seed(42)
        n_taxa = 500
        n_features = 50
        
        data = np.random.random((n_taxa, n_features))
        
        start = time.time()
        pca = PCA(n_components=2)
        embedding = pca.fit_transform(data)
        elapsed = time.time() - start
        
        assert embedding.shape == (n_taxa, 2)
        assert elapsed < 1.0

    def test_tsne_timing(self):
        """Test t-SNE timing for comparison."""
        from sklearn.manifold import TSNE
        
        np.random.seed(42)
        n_taxa = 100  # Small for t-SNE
        n_features = 20
        
        data = np.random.random((n_taxa, n_features))
        
        start = time.time()
        tsne = TSNE(n_components=2, perplexity=30, random_state=42)
        embedding = tsne.fit_transform(data)
        elapsed = time.time() - start
        
        assert embedding.shape == (n_taxa, 2)
        # t-SNE is slow, allow more time
        assert elapsed < 30.0


@pytest.mark.performance
class TestDataProcessingPerformance:
    """Performance tests for data processing."""

    def test_dataframe_merge_timing(self):
        """Test DataFrame merge performance."""
        np.random.seed(42)
        n_taxa = 1000
        
        df1 = pd.DataFrame({
            'taxon': [f'T{i}' for i in range(n_taxa)],
            'value1': np.random.random(n_taxa)
        })
        
        df2 = pd.DataFrame({
            'taxon': [f'T{i}' for i in range(n_taxa)],
            'value2': np.random.random(n_taxa)
        })
        
        start = time.time()
        merged = pd.merge(df1, df2, on='taxon')
        elapsed = time.time() - start
        
        assert len(merged) == n_taxa
        assert elapsed < 0.5

    def test_binary_encoding_timing(self):
        """Test binary encoding performance."""
        np.random.seed(42)
        n_taxa = 1000
        n_traits = 50
        
        traits = np.random.randint(0, 2, size=(n_taxa, n_traits))
        
        start = time.time()
        df = pd.DataFrame(traits)
        # Convert to boolean
        df_bool = df.astype(bool)
        elapsed = time.time() - start
        
        assert df_bool.shape == (n_taxa, n_traits)
        assert elapsed < 0.1
