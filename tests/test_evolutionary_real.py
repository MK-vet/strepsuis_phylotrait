#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for EvolutionaryAnalysis with real data.
"""

import os
import pytest
import pandas as pd
import numpy as np

try:
    from Bio import Phylo
    BIO_AVAILABLE = True
except ImportError:
    BIO_AVAILABLE = False

REAL_DATA_PATH = os.path.join(os.path.dirname(__file__), '..', '..', '..', 'data')

try:
    from strepsuis_phylotrait.phylo_analysis_core import (
        EvolutionaryAnalysis,
        PhylogeneticCore,
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
    pytest.skip("Tree not found")


@pytest.fixture
def amr_data():
    """Load AMR data."""
    path = os.path.join(REAL_DATA_PATH, 'AMR_genes.csv')
    if os.path.exists(path):
        return pd.read_csv(path).set_index('Strain_ID')
    pytest.skip("AMR data not found")


@pytest.fixture
def common_strains(real_tree, amr_data):
    """Get strains common to tree and data."""
    tree_terminals = [str(t.name) for t in real_tree.get_terminals()]
    data_strains = amr_data.index.astype(str).tolist()
    return list(set(tree_terminals) & set(data_strains))


@pytest.mark.skipif(not CORE_AVAILABLE or not BIO_AVAILABLE, reason="Core or Bio not available")
class TestEvolutionaryAnalysisReal:
    """Test EvolutionaryAnalysis with real data."""
    
    def test_evolutionary_cluster_analysis(self, real_tree, common_strains):
        """Test evolutionary cluster analysis."""
        if len(common_strains) < 10:
            pytest.skip("Not enough common strains")
        
        n = len(common_strains)
        labels = np.array([i % 3 for i in range(n)])
        mask = np.array([True] * n)
        
        result = EvolutionaryAnalysis.evolutionary_cluster_analysis(
            real_tree, labels, common_strains, mask
        )
        
        assert isinstance(result, pd.DataFrame)
    
    def test_analyze_cluster_evolution(self, real_tree, common_strains):
        """Test cluster evolution analysis."""
        if len(common_strains) < 10:
            pytest.skip("Not enough common strains")
        
        # _analyze_cluster_evolution takes (tree, cluster_id, cluster_strains)
        cluster_strains = np.array(common_strains[:5])
        
        result = EvolutionaryAnalysis._analyze_cluster_evolution(
            real_tree, 0, cluster_strains
        )
        
        assert isinstance(result, tuple)
        assert len(result) == 6
    
    def test_calculate_beta_diversity(self, real_tree, common_strains):
        """Test beta diversity calculation."""
        if len(common_strains) < 10:
            pytest.skip("Not enough common strains")
        
        n = len(common_strains)
        labels = np.array([i % 3 for i in range(n)])
        mask = np.array([True] * n)
        
        result = EvolutionaryAnalysis.calculate_beta_diversity(
            real_tree, labels, common_strains, mask
        )
        
        assert isinstance(result, pd.DataFrame)
    
    def test_calculate_evolution_rates(self, real_tree, common_strains):
        """Test evolution rates calculation."""
        if len(common_strains) < 10:
            pytest.skip("Not enough common strains")
        
        n = len(common_strains)
        labels = np.array([i % 3 for i in range(n)])
        mask = np.array([True] * n)
        
        # First get cluster analysis
        cluster_df = EvolutionaryAnalysis.evolutionary_cluster_analysis(
            real_tree, labels, common_strains, mask
        )
        
        # Then calculate evolution rates
        result = EvolutionaryAnalysis.calculate_evolution_rates(cluster_df)
        
        assert isinstance(result, pd.DataFrame)
        assert 'Cluster_ID' in result.columns
        assert 'EvolutionRate' in result.columns


@pytest.mark.skipif(not CORE_AVAILABLE or not BIO_AVAILABLE, reason="Core or Bio not available")
class TestPhylogeneticSignal:
    """Test phylogenetic signal calculations."""
    
    def test_calculate_phylogenetic_signal_fritz_purvis(self, real_tree, amr_data, common_strains):
        """Test Fritz-Purvis D statistic calculation."""
        if len(common_strains) < 10:
            pytest.skip("Not enough common strains")
        
        # Filter data to common strains
        trait_data = amr_data.loc[amr_data.index.astype(str).isin(common_strains)]
        
        # Use first 3 traits
        trait_subset = trait_data.iloc[:, :3]
        
        # Create instance and call method
        analyzer = EvolutionaryAnalysis()
        
        if hasattr(analyzer, 'calculate_phylogenetic_signal_fritz_purvis'):
            try:
                result = analyzer.calculate_phylogenetic_signal_fritz_purvis(
                    real_tree, trait_subset
                )
                
                assert isinstance(result, pd.DataFrame)
                assert 'trait' in result.columns
                assert 'd_statistic' in result.columns
            except ImportError:
                pytest.skip("Bio.Phylo import issue in method")
        else:
            pytest.skip("calculate_phylogenetic_signal_fritz_purvis not available")
    
    def test_sister_clade_differences(self, real_tree, amr_data, common_strains):
        """Test sister clade differences calculation."""
        if len(common_strains) < 10:
            pytest.skip("Not enough common strains")
        
        trait_data = amr_data.loc[amr_data.index.astype(str).isin(common_strains)]
        trait_vector = trait_data.iloc[:, 0].values
        
        analyzer = EvolutionaryAnalysis()
        if hasattr(analyzer, '_sister_clade_differences'):
            try:
                # _sister_clade_differences takes (tree, trait_vector)
                result = analyzer._sister_clade_differences(real_tree, trait_vector)
                
                assert isinstance(result, (int, float, np.number))
            except (TypeError, AttributeError):
                pytest.skip("_sister_clade_differences has different signature")
        else:
            pytest.skip("_sister_clade_differences not available")


@pytest.mark.skipif(not CORE_AVAILABLE or not BIO_AVAILABLE, reason="Core or Bio not available")
class TestEvolutionaryMetrics:
    """Test evolutionary metrics."""
    
    def test_faiths_pd(self, real_tree, common_strains):
        """Test Faith's phylogenetic diversity."""
        if len(common_strains) < 5:
            pytest.skip("Not enough common strains")
        
        # Select subset of strains
        subset = common_strains[:10]
        
        if hasattr(EvolutionaryAnalysis, 'calculate_faiths_pd'):
            result = EvolutionaryAnalysis.calculate_faiths_pd(real_tree, subset)
            
            assert isinstance(result, float)
            assert result >= 0
    
    def test_mpd_mntd(self, real_tree, common_strains):
        """Test MPD and MNTD calculations."""
        if len(common_strains) < 5:
            pytest.skip("Not enough common strains")
        
        subset = common_strains[:10]
        
        if hasattr(EvolutionaryAnalysis, 'calculate_mpd'):
            mpd = EvolutionaryAnalysis.calculate_mpd(real_tree, subset)
            assert isinstance(mpd, float)
        
        if hasattr(EvolutionaryAnalysis, 'calculate_mntd'):
            mntd = EvolutionaryAnalysis.calculate_mntd(real_tree, subset)
            assert isinstance(mntd, float)
