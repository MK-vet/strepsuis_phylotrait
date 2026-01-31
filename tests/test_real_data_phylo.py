#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests using real S. suis data (91 strains) for phylogenetic analysis.
"""

import os
import pytest
import pandas as pd
import numpy as np

try:
    from Bio import Phylo
    from io import StringIO
    BIO_AVAILABLE = True
except ImportError:
    BIO_AVAILABLE = False

# Path to real data
REAL_DATA_PATH = os.path.join(os.path.dirname(__file__), '..', '..', '..', 'data')

try:
    from strepsuis_phylotrait.phylo_analysis_core import (
        setup_logging,
        print_memory_usage,
        print_section_header,
        print_step,
        PhylogeneticCore,
        ClusteringModule,
        EvolutionaryAnalysis,
        DataLoader,
        Visualizer,
        TraitAnalyzer,
        MCAAnalyzer,
        Config,
    )
    CORE_AVAILABLE = True
except (ImportError, OSError):
    CORE_AVAILABLE = False


@pytest.fixture
def real_tree():
    """Load real SNP tree."""
    tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
    if os.path.exists(tree_path) and BIO_AVAILABLE:
        return Phylo.read(tree_path, 'newick')
    pytest.skip("Snp_tree.newick not found or Bio not available")


@pytest.fixture
def amr_data():
    """Load real AMR genes data."""
    path = os.path.join(REAL_DATA_PATH, 'AMR_genes.csv')
    if os.path.exists(path):
        df = pd.read_csv(path)
        df = df.set_index('Strain_ID')
        return df
    pytest.skip("AMR_genes.csv not found")


@pytest.fixture
def virulence_data():
    """Load real virulence data."""
    path = os.path.join(REAL_DATA_PATH, 'Virulence.csv')
    if os.path.exists(path):
        df = pd.read_csv(path)
        df = df.set_index('Strain_ID')
        return df
    pytest.skip("Virulence.csv not found")


@pytest.fixture
def mic_data():
    """Load real MIC data."""
    path = os.path.join(REAL_DATA_PATH, 'MIC.csv')
    if os.path.exists(path):
        df = pd.read_csv(path)
        df = df.set_index('Strain_ID')
        return df
    pytest.skip("MIC.csv not found")


@pytest.mark.skipif(not BIO_AVAILABLE, reason="Biopython not available")
class TestRealTreeLoading:
    """Test loading real phylogenetic tree."""
    
    def test_load_real_tree(self, real_tree):
        """Test loading real SNP tree."""
        assert real_tree is not None
        
        terminals = real_tree.get_terminals()
        assert len(terminals) > 50  # Should have many strains
    
    def test_real_tree_terminal_names(self, real_tree):
        """Test terminal names in real tree."""
        terminals = real_tree.get_terminals()
        names = [t.name for t in terminals]
        
        # Should have strain IDs
        assert all(name is not None for name in names)
    
    def test_real_tree_distances(self, real_tree):
        """Test calculating distances in real tree."""
        terminals = real_tree.get_terminals()
        
        if len(terminals) >= 2:
            dist = real_tree.distance(terminals[0], terminals[1])
            assert isinstance(dist, float)
            assert dist >= 0


@pytest.mark.skipif(not CORE_AVAILABLE or not BIO_AVAILABLE, reason="Core or Bio not available")
class TestPhylogeneticCoreRealData:
    """Test PhylogeneticCore with real data."""
    
    def test_load_real_tree_via_core(self):
        """Test loading real tree via PhylogeneticCore."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        if not os.path.exists(tree_path):
            pytest.skip("Tree file not found")
        
        tree = PhylogeneticCore.load_tree(tree_path)
        
        assert tree is not None
        terminals = tree.get_terminals()
        assert len(terminals) > 50


@pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
class TestTraitAnalyzerRealData:
    """Test TraitAnalyzer with real data."""
    
    def test_trait_analyzer_amr(self, amr_data):
        """Test TraitAnalyzer with real AMR data."""
        analyzer = TraitAnalyzer(amr_data)
        
        assert analyzer is not None
    
    def test_trait_analyzer_virulence(self, virulence_data):
        """Test TraitAnalyzer with real virulence data."""
        analyzer = TraitAnalyzer(virulence_data)
        
        assert analyzer is not None


@pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
class TestMCAAnalyzerRealData:
    """Test MCAAnalyzer with real data."""
    
    def test_mca_analyzer_amr(self, amr_data):
        """Test MCAAnalyzer with real AMR data."""
        mca = MCAAnalyzer(amr_data)
        
        assert mca is not None
    
    def test_mca_analyzer_combined(self, amr_data, virulence_data):
        """Test MCAAnalyzer with combined data."""
        # Combine AMR and virulence data
        combined = amr_data.join(virulence_data, how='inner', rsuffix='_vir')
        
        mca = MCAAnalyzer(combined)
        
        assert mca is not None


@pytest.mark.skipif(not CORE_AVAILABLE or not BIO_AVAILABLE, reason="Core or Bio not available")
class TestEvolutionaryAnalysisRealData:
    """Test EvolutionaryAnalysis with real data."""
    
    def test_evolutionary_cluster_analysis_real(self, real_tree, amr_data):
        """Test evolutionary cluster analysis with real data."""
        # Get common strains
        tree_terminals = [t.name for t in real_tree.get_terminals()]
        common_strains = [s for s in amr_data.index.astype(str) if s in tree_terminals]
        
        if len(common_strains) < 10:
            pytest.skip("Not enough common strains")
        
        # Create simple cluster labels
        n_strains = len(common_strains)
        labels = np.array([i % 3 for i in range(n_strains)])
        mask = np.array([True] * n_strains)
        
        result = EvolutionaryAnalysis.evolutionary_cluster_analysis(
            real_tree, labels, common_strains, mask
        )
        
        assert isinstance(result, pd.DataFrame)
    
    def test_calculate_beta_diversity_real(self, real_tree, amr_data):
        """Test beta diversity calculation with real data."""
        tree_terminals = [t.name for t in real_tree.get_terminals()]
        common_strains = [s for s in amr_data.index.astype(str) if s in tree_terminals]
        
        if len(common_strains) < 10:
            pytest.skip("Not enough common strains")
        
        n_strains = len(common_strains)
        labels = np.array([i % 3 for i in range(n_strains)])
        mask = np.array([True] * n_strains)
        
        result = EvolutionaryAnalysis.calculate_beta_diversity(
            real_tree, labels, common_strains, mask
        )
        
        assert isinstance(result, pd.DataFrame)


class TestRealDataStatistics:
    """Test statistical operations on real data."""
    
    def test_amr_prevalence(self, amr_data):
        """Test AMR gene prevalence calculation."""
        prevalence = amr_data.mean()
        
        assert isinstance(prevalence, pd.Series)
        assert all(0 <= p <= 1 for p in prevalence)
    
    def test_virulence_prevalence(self, virulence_data):
        """Test virulence factor prevalence."""
        prevalence = virulence_data.mean()
        
        assert isinstance(prevalence, pd.Series)
        assert all(0 <= p <= 1 for p in prevalence)
    
    def test_amr_cooccurrence(self, amr_data):
        """Test AMR gene co-occurrence."""
        cooccurrence = amr_data.T.dot(amr_data)
        
        assert cooccurrence.shape[0] == cooccurrence.shape[1]
        assert cooccurrence.shape[0] == len(amr_data.columns)
    
    def test_amr_correlation(self, amr_data):
        """Test AMR gene correlation."""
        # Filter genes with variance
        var_genes = amr_data.columns[amr_data.var() > 0]
        
        if len(var_genes) >= 2:
            correlation = amr_data[var_genes].corr()
            
            assert correlation.shape[0] == len(var_genes)
            assert np.allclose(np.diag(correlation.values), 1.0)


class TestRealDataClustering:
    """Test clustering operations on real data."""
    
    def test_hamming_distance_amr(self, amr_data):
        """Test Hamming distance calculation on AMR data."""
        from scipy.spatial.distance import pdist, squareform
        
        dist_matrix = squareform(pdist(amr_data.values, metric='hamming'))
        
        assert dist_matrix.shape == (len(amr_data), len(amr_data))
        assert np.allclose(np.diag(dist_matrix), 0)
    
    def test_kmodes_clustering_amr(self, amr_data):
        """Test K-Modes clustering on AMR data."""
        try:
            from kmodes.kmodes import KModes
            
            km = KModes(n_clusters=3, random_state=42, n_init=5)
            labels = km.fit_predict(amr_data.values)
            
            assert len(labels) == len(amr_data)
            assert len(np.unique(labels)) <= 3
        except ImportError:
            pytest.skip("kmodes not available")


@pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
class TestConfigRealData:
    """Test Config with real data paths."""
    
    def test_config_with_real_paths(self):
        """Test Config with real data paths."""
        config = Config()
        
        config.data_folder = REAL_DATA_PATH
        config.tree_file = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        config.n_clusters = 5
        config.random_state = 42
        
        assert config.data_folder == REAL_DATA_PATH
        assert os.path.exists(config.data_folder)
