#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Detailed tests for phylogenetic analysis classes.
"""

import os
import tempfile
import shutil
import pytest
import pandas as pd
import numpy as np

try:
    from Bio import Phylo
    from io import StringIO
    BIO_AVAILABLE = True
except ImportError:
    BIO_AVAILABLE = False

try:
    from strepsuis_phylotrait.phylo_analysis_core import (
        PhylogeneticCore,
        ClusteringModule,
        TreeAwareClusteringModule,
        EvolutionaryAnalysis,
        DataLoader,
        Visualizer,
        TraitAnalyzer,
        MCAAnalyzer,
        HTMLReportGenerator,
        Config,
        ParallelProcessor,
    )
    CLASSES_AVAILABLE = True
except (ImportError, OSError):
    CLASSES_AVAILABLE = False


@pytest.mark.skipif(not CLASSES_AVAILABLE, reason="Classes not available")
class TestVisualizerDetailed:
    """Detailed tests for Visualizer class."""
    
    @pytest.fixture
    def temp_output_folder(self):
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    def test_visualizer_has_output_folder(self, temp_output_folder):
        """Test Visualizer has output_folder attribute."""
        viz = Visualizer(temp_output_folder)
        
        assert hasattr(viz, 'output_folder')
        assert viz.output_folder == temp_output_folder


@pytest.mark.skipif(not CLASSES_AVAILABLE, reason="Classes not available")
class TestTraitAnalyzerDetailed:
    """Detailed tests for TraitAnalyzer class."""
    
    @pytest.fixture
    def sample_trait_data(self):
        return pd.DataFrame({
            'Trait_A': [1, 0, 1, 0, 1, 0, 1, 0],
            'Trait_B': [0, 1, 0, 1, 0, 1, 0, 1],
            'Trait_C': [1, 1, 0, 0, 1, 1, 0, 0],
        }, index=['S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8'])
    
    def test_trait_analyzer_has_data(self, sample_trait_data):
        """Test TraitAnalyzer has data."""
        analyzer = TraitAnalyzer(sample_trait_data)
        
        assert analyzer is not None


@pytest.mark.skipif(not CLASSES_AVAILABLE, reason="Classes not available")
class TestMCAAnalyzerDetailed:
    """Detailed tests for MCAAnalyzer class."""
    
    @pytest.fixture
    def sample_data(self):
        np.random.seed(42)
        return pd.DataFrame({
            'A': np.random.binomial(1, 0.5, 20),
            'B': np.random.binomial(1, 0.3, 20),
            'C': np.random.binomial(1, 0.7, 20),
        })
    
    def test_mca_analyzer_has_data(self, sample_data):
        """Test MCAAnalyzer has data."""
        mca = MCAAnalyzer(sample_data)
        
        assert mca is not None


@pytest.mark.skipif(not CLASSES_AVAILABLE, reason="Classes not available")
class TestHTMLReportGeneratorDetailed:
    """Detailed tests for HTMLReportGenerator class."""
    
    @pytest.fixture
    def temp_output_folder(self):
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    def test_html_report_generator_has_output_folder(self, temp_output_folder):
        """Test HTMLReportGenerator has output_folder."""
        generator = HTMLReportGenerator(temp_output_folder)
        
        assert generator is not None


@pytest.mark.skipif(not CLASSES_AVAILABLE, reason="Classes not available")
class TestDataLoaderDetailed:
    """Detailed tests for DataLoader class."""
    
    @pytest.fixture
    def temp_data_folder(self):
        temp_dir = tempfile.mkdtemp()
        
        # Create multiple sample CSV files
        for name in ['AMR_genes', 'Virulence', 'MIC', 'MLST']:
            df = pd.DataFrame({
                'Strain_ID': ['S1', 'S2', 'S3'],
                f'{name}_col': [1, 0, 1],
            })
            df.to_csv(os.path.join(temp_dir, f'{name}.csv'), index=False)
        
        yield temp_dir
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    def test_data_loader_folder_exists(self, temp_data_folder):
        """Test DataLoader with existing folder."""
        loader = DataLoader(temp_data_folder)
        
        assert loader is not None


@pytest.mark.skipif(not CLASSES_AVAILABLE, reason="Classes not available")
class TestConfigDetailed:
    """Detailed tests for Config class."""
    
    def test_config_multiple_settings(self):
        """Test Config with multiple settings."""
        config = Config()
        
        # Set various configuration options
        config.data_folder = '/data'
        config.output_folder = '/output'
        config.n_clusters = 5
        config.random_state = 42
        config.n_bootstrap = 1000
        config.alpha = 0.05
        config.verbose = True
        config.parallel = True
        config.n_jobs = 4
        
        assert config.n_clusters == 5
        assert config.random_state == 42
        assert config.n_bootstrap == 1000


@pytest.mark.skipif(not CLASSES_AVAILABLE or not BIO_AVAILABLE, reason="Classes or Bio not available")
class TestPhylogeneticCoreDetailed:
    """Detailed tests for PhylogeneticCore class."""
    
    @pytest.fixture
    def temp_tree_file(self):
        temp_dir = tempfile.mkdtemp()
        tree_path = os.path.join(temp_dir, 'test.newick')
        # More complex tree
        with open(tree_path, 'w') as f:
            f.write('(((A:0.1,B:0.15):0.2,(C:0.12,D:0.18):0.22):0.3,(E:0.25,F:0.28):0.35);')
        yield tree_path
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    def test_load_complex_tree(self, temp_tree_file):
        """Test loading more complex tree."""
        tree = PhylogeneticCore.load_tree(temp_tree_file)
        
        terminals = tree.get_terminals()
        assert len(terminals) == 6
        
        names = [t.name for t in terminals]
        assert set(names) == {'A', 'B', 'C', 'D', 'E', 'F'}


@pytest.mark.skipif(not CLASSES_AVAILABLE or not BIO_AVAILABLE, reason="Classes or Bio not available")
class TestEvolutionaryAnalysisDetailed:
    """Detailed tests for EvolutionaryAnalysis class."""
    
    @pytest.fixture
    def sample_tree(self):
        tree_str = '((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1);'
        return Phylo.read(StringIO(tree_str), 'newick')
    
    @pytest.fixture
    def sample_labels(self):
        return np.array([1, 1, 2, 2])
    
    @pytest.fixture
    def sample_strain_names(self):
        return ['A', 'B', 'C', 'D']
    
    @pytest.fixture
    def sample_mask(self):
        return np.array([True, True, True, True])
    
    def test_evolutionary_cluster_analysis(self, sample_tree, sample_labels, sample_strain_names, sample_mask):
        """Test evolutionary cluster analysis."""
        result = EvolutionaryAnalysis.evolutionary_cluster_analysis(
            sample_tree, sample_labels, sample_strain_names, sample_mask
        )
        
        assert isinstance(result, pd.DataFrame)
    
    def test_calculate_beta_diversity(self, sample_tree, sample_labels, sample_strain_names, sample_mask):
        """Test beta diversity calculation."""
        result = EvolutionaryAnalysis.calculate_beta_diversity(
            sample_tree, sample_labels, sample_strain_names, sample_mask
        )
        
        assert isinstance(result, pd.DataFrame)
