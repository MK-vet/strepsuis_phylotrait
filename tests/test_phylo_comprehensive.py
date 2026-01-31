#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Comprehensive tests for phylogenetic analysis module.
"""

import os
import tempfile
import shutil
import pytest
import pandas as pd
import numpy as np
import logging

try:
    from Bio import Phylo
    from io import StringIO
    BIO_AVAILABLE = True
except ImportError:
    BIO_AVAILABLE = False

try:
    from strepsuis_phylotrait.phylo_analysis_core import (
        setup_logging,
        print_memory_usage,
        print_section_header,
        print_step,
        create_template_directory,
        PhylogeneticCore,
        ParallelProcessor,
        ClusteringModule,
        TreeAwareClusteringModule,
        EvolutionaryAnalysis,
        DataLoader,
        Visualizer,
        TraitAnalyzer,
        MCAAnalyzer,
        HTMLReportGenerator,
        Config,
    )
    CORE_AVAILABLE = True
except (ImportError, OSError) as e:
    CORE_AVAILABLE = False


@pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
class TestLoggingComprehensive:
    """Comprehensive logging tests."""
    
    @pytest.fixture
    def temp_folder(self):
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    def test_setup_logging_with_subfolder(self, temp_folder):
        """Test logging with nested subfolder."""
        subfolder = os.path.join(temp_folder, 'logs', 'analysis')
        os.makedirs(subfolder, exist_ok=True)
        
        log_file = setup_logging(subfolder)
        
        assert isinstance(log_file, str)
    
    def test_print_functions_sequence(self):
        """Test sequence of print functions."""
        print_section_header("Analysis Start")
        print_step(1, 5, "Loading data")
        mem1 = print_memory_usage()
        print_step(2, 5, "Processing")
        mem2 = print_memory_usage()
        print_step(3, 5, "Clustering")
        print_step(4, 5, "Visualization")
        print_step(5, 5, "Report generation")
        print_section_header("Analysis Complete")
        
        assert mem1 > 0
        assert mem2 > 0


@pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
class TestTemplateComprehensive:
    """Comprehensive template tests."""
    
    @pytest.fixture
    def temp_cwd(self):
        original_cwd = os.getcwd()
        temp_dir = tempfile.mkdtemp()
        os.chdir(temp_dir)
        yield temp_dir
        os.chdir(original_cwd)
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    def test_create_template_multiple_times(self, temp_cwd):
        """Test creating template multiple times."""
        create_template_directory()
        create_template_directory()
        create_template_directory()
        
        assert os.path.exists('templates')
        files = os.listdir('templates')
        assert len(files) > 0


@pytest.mark.skipif(not CORE_AVAILABLE or not BIO_AVAILABLE, reason="Core or Bio not available")
class TestPhylogeneticCoreComprehensive:
    """Comprehensive PhylogeneticCore tests."""
    
    @pytest.fixture
    def simple_tree_file(self):
        temp_dir = tempfile.mkdtemp()
        tree_path = os.path.join(temp_dir, 'simple.newick')
        with open(tree_path, 'w') as f:
            f.write('(A:0.1,B:0.2);')
        yield tree_path
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    @pytest.fixture
    def complex_tree_file(self):
        temp_dir = tempfile.mkdtemp()
        tree_path = os.path.join(temp_dir, 'complex.newick')
        with open(tree_path, 'w') as f:
            f.write('(((A:0.1,B:0.15):0.2,(C:0.12,D:0.18):0.22):0.3,(E:0.25,F:0.28):0.35);')
        yield tree_path
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    def test_load_simple_tree(self, simple_tree_file):
        """Test loading simple tree."""
        tree = PhylogeneticCore.load_tree(simple_tree_file)
        
        terminals = tree.get_terminals()
        assert len(terminals) == 2
    
    def test_load_complex_tree(self, complex_tree_file):
        """Test loading complex tree."""
        tree = PhylogeneticCore.load_tree(complex_tree_file)
        
        terminals = tree.get_terminals()
        assert len(terminals) == 6


@pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
class TestConfigComprehensive:
    """Comprehensive Config tests."""
    
    def test_config_all_attributes(self):
        """Test Config with all common attributes."""
        config = Config()
        
        # Data paths
        config.data_folder = '/data'
        config.output_folder = '/output'
        config.tree_file = '/tree.newick'
        
        # Analysis parameters
        config.n_clusters = 5
        config.min_clusters = 2
        config.max_clusters = 10
        config.random_state = 42
        
        # Bootstrap parameters
        config.n_bootstrap = 1000
        config.alpha = 0.05
        config.ci_level = 0.95
        
        # Execution parameters
        config.verbose = True
        config.parallel = True
        config.n_jobs = 4
        config.timeout = 3600
        
        # Verify all attributes
        assert config.n_clusters == 5
        assert config.random_state == 42
        assert config.n_bootstrap == 1000
        assert config.parallel == True


@pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
class TestVisualizerComprehensive:
    """Comprehensive Visualizer tests."""
    
    @pytest.fixture
    def temp_output_folder(self):
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    def test_visualizer_with_subfolders(self, temp_output_folder):
        """Test Visualizer with subfolders."""
        subfolder = os.path.join(temp_output_folder, 'plots', 'clustering')
        os.makedirs(subfolder, exist_ok=True)
        
        viz = Visualizer(subfolder)
        
        assert viz is not None
        assert viz.output_folder == subfolder


@pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
class TestTraitAnalyzerComprehensive:
    """Comprehensive TraitAnalyzer tests."""
    
    @pytest.fixture
    def large_trait_data(self):
        np.random.seed(42)
        n_samples = 50
        n_traits = 20
        
        data = {}
        for i in range(n_traits):
            data[f'Trait_{i}'] = np.random.binomial(1, np.random.uniform(0.2, 0.8), n_samples)
        
        return pd.DataFrame(data, index=[f'S{i}' for i in range(n_samples)])
    
    def test_trait_analyzer_large_data(self, large_trait_data):
        """Test TraitAnalyzer with large data."""
        analyzer = TraitAnalyzer(large_trait_data)
        
        assert analyzer is not None


@pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
class TestMCAAnalyzerComprehensive:
    """Comprehensive MCAAnalyzer tests."""
    
    @pytest.fixture
    def large_data(self):
        np.random.seed(42)
        n_samples = 50
        n_features = 15
        
        data = {}
        for i in range(n_features):
            data[f'Feature_{i}'] = np.random.binomial(1, np.random.uniform(0.3, 0.7), n_samples)
        
        return pd.DataFrame(data)
    
    def test_mca_analyzer_large_data(self, large_data):
        """Test MCAAnalyzer with large data."""
        mca = MCAAnalyzer(large_data)
        
        assert mca is not None


@pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
class TestHTMLReportGeneratorComprehensive:
    """Comprehensive HTMLReportGenerator tests."""
    
    @pytest.fixture
    def temp_output_folder(self):
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    def test_html_report_generator_with_subfolders(self, temp_output_folder):
        """Test HTMLReportGenerator with subfolders."""
        subfolder = os.path.join(temp_output_folder, 'reports', 'html')
        os.makedirs(subfolder, exist_ok=True)
        
        generator = HTMLReportGenerator(subfolder)
        
        assert generator is not None


@pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
class TestDataLoaderComprehensive:
    """Comprehensive DataLoader tests."""
    
    @pytest.fixture
    def temp_data_folder(self):
        temp_dir = tempfile.mkdtemp()
        
        # Create multiple sample CSV files
        for name, n_cols in [('AMR_genes', 10), ('Virulence', 8), ('MIC', 5), ('MLST', 3)]:
            data = {'Strain_ID': ['S1', 'S2', 'S3', 'S4', 'S5']}
            for i in range(n_cols):
                data[f'{name}_col_{i}'] = np.random.binomial(1, 0.5, 5)
            df = pd.DataFrame(data)
            df.to_csv(os.path.join(temp_dir, f'{name}.csv'), index=False)
        
        yield temp_dir
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    def test_data_loader_multiple_files(self, temp_data_folder):
        """Test DataLoader with multiple files."""
        loader = DataLoader(temp_data_folder)
        
        assert loader is not None


@pytest.mark.skipif(not CORE_AVAILABLE or not BIO_AVAILABLE, reason="Core or Bio not available")
class TestEvolutionaryAnalysisComprehensive:
    """Comprehensive EvolutionaryAnalysis tests."""
    
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
    
    def test_evolutionary_cluster_analysis_basic(self, sample_tree, sample_labels, sample_strain_names, sample_mask):
        """Test basic evolutionary cluster analysis."""
        result = EvolutionaryAnalysis.evolutionary_cluster_analysis(
            sample_tree, sample_labels, sample_strain_names, sample_mask
        )
        
        assert isinstance(result, pd.DataFrame)
    
    def test_calculate_beta_diversity_basic(self, sample_tree, sample_labels, sample_strain_names, sample_mask):
        """Test basic beta diversity calculation."""
        result = EvolutionaryAnalysis.calculate_beta_diversity(
            sample_tree, sample_labels, sample_strain_names, sample_mask
        )
        
        assert isinstance(result, pd.DataFrame)
