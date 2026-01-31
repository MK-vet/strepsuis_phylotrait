#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests that directly call phylo_analysis_core functions to increase coverage.
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
    print(f"Import error: {e}")


@pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
class TestSetupLoggingCoverage:
    """Tests to increase setup_logging coverage."""
    
    @pytest.fixture
    def temp_folder(self):
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    def test_setup_logging_full(self, temp_folder):
        """Test full setup_logging functionality."""
        log_file = setup_logging(temp_folder)
        
        # Log messages at different levels
        logging.info("Info message")
        logging.warning("Warning message")
        logging.debug("Debug message")
        
        assert os.path.exists(log_file)


@pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
class TestPrintFunctionsCoverage:
    """Tests to increase print functions coverage."""
    
    def test_print_memory_usage_multiple(self):
        """Test print_memory_usage multiple times."""
        mem1 = print_memory_usage()
        # Allocate some memory
        _ = [i for i in range(10000)]
        mem2 = print_memory_usage()
        
        assert mem1 > 0
        assert mem2 > 0
    
    def test_print_section_header_various(self):
        """Test print_section_header with various inputs."""
        print_section_header("Test")
        print_section_header("A" * 50)
        print_section_header("Special: !@#$%")
    
    def test_print_step_various(self):
        """Test print_step with various inputs."""
        print_step(0, 10, "Zero step")
        print_step(1, 1, "Single step")
        print_step(100, 100, "Last of many")


@pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
class TestCreateTemplateDirectoryCoverage:
    """Tests to increase create_template_directory coverage."""
    
    @pytest.fixture
    def temp_cwd(self):
        original_cwd = os.getcwd()
        temp_dir = tempfile.mkdtemp()
        os.chdir(temp_dir)
        yield temp_dir
        os.chdir(original_cwd)
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    def test_create_template_directory_twice(self, temp_cwd):
        """Test create_template_directory called twice."""
        create_template_directory()
        create_template_directory()  # Should not raise error
        
        assert os.path.exists('templates')


@pytest.mark.skipif(not CORE_AVAILABLE or not BIO_AVAILABLE, reason="Core or Bio not available")
class TestPhylogeneticCoreCoverage:
    """Tests to increase PhylogeneticCore coverage."""
    
    @pytest.fixture
    def temp_tree_file(self):
        temp_dir = tempfile.mkdtemp()
        tree_path = os.path.join(temp_dir, 'test.newick')
        with open(tree_path, 'w') as f:
            f.write('((A:0.1,B:0.2):0.3,(C:0.15,D:0.25):0.35);')
        yield tree_path
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    def test_load_tree_and_get_terminals(self, temp_tree_file):
        """Test loading tree and getting terminals."""
        tree = PhylogeneticCore.load_tree(temp_tree_file)
        
        terminals = tree.get_terminals()
        assert len(terminals) == 4
        
        # Get terminal names
        names = [t.name for t in terminals]
        assert set(names) == {'A', 'B', 'C', 'D'}


@pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
class TestConfigCoverage:
    """Tests to increase Config coverage."""
    
    def test_config_multiple_attributes(self):
        """Test Config with multiple attributes."""
        config = Config()
        
        config.data_folder = '/path/to/data'
        config.output_folder = '/path/to/output'
        config.n_clusters = 5
        config.random_state = 42
        config.verbose = True
        
        assert config.data_folder == '/path/to/data'
        assert config.n_clusters == 5


@pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
class TestVisualizerCoverage:
    """Tests to increase Visualizer coverage."""
    
    @pytest.fixture
    def temp_output_folder(self):
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    def test_visualizer_init(self, temp_output_folder):
        """Test Visualizer initialization."""
        viz = Visualizer(temp_output_folder)
        
        assert viz is not None
        assert viz.output_folder == temp_output_folder


@pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
class TestTraitAnalyzerCoverage:
    """Tests to increase TraitAnalyzer coverage."""
    
    @pytest.fixture
    def sample_trait_data(self):
        return pd.DataFrame({
            'Trait_A': [1, 0, 1, 0, 1, 0],
            'Trait_B': [0, 1, 0, 1, 0, 1],
            'Trait_C': [1, 1, 0, 0, 1, 1],
        }, index=['S1', 'S2', 'S3', 'S4', 'S5', 'S6'])
    
    def test_trait_analyzer_init(self, sample_trait_data):
        """Test TraitAnalyzer initialization."""
        analyzer = TraitAnalyzer(sample_trait_data)
        
        assert analyzer is not None


@pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
class TestMCAAnalyzerCoverage:
    """Tests to increase MCAAnalyzer coverage."""
    
    @pytest.fixture
    def sample_data(self):
        return pd.DataFrame({
            'A': [1, 0, 1, 0, 1, 0],
            'B': [0, 1, 0, 1, 0, 1],
            'C': [1, 1, 0, 0, 1, 1],
        })
    
    def test_mca_analyzer_init(self, sample_data):
        """Test MCAAnalyzer initialization."""
        mca = MCAAnalyzer(sample_data)
        
        assert mca is not None


@pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
class TestHTMLReportGeneratorCoverage:
    """Tests to increase HTMLReportGenerator coverage."""
    
    @pytest.fixture
    def temp_output_folder(self):
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    def test_html_report_generator_init(self, temp_output_folder):
        """Test HTMLReportGenerator initialization."""
        generator = HTMLReportGenerator(temp_output_folder)
        
        assert generator is not None


@pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
class TestDataLoaderCoverage:
    """Tests to increase DataLoader coverage."""
    
    @pytest.fixture
    def temp_data_folder(self):
        temp_dir = tempfile.mkdtemp()
        
        # Create sample CSV files
        amr_df = pd.DataFrame({
            'Strain_ID': ['S1', 'S2', 'S3', 'S4'],
            'Gene_A': [1, 0, 1, 0],
            'Gene_B': [0, 1, 0, 1],
        })
        amr_df.to_csv(os.path.join(temp_dir, 'AMR_genes.csv'), index=False)
        
        yield temp_dir
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    def test_data_loader_init(self, temp_data_folder):
        """Test DataLoader initialization."""
        loader = DataLoader(temp_data_folder)
        
        assert loader is not None
