#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for core phylogenetic functions.
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
        Config,
    )
    CORE_AVAILABLE = True
except (ImportError, OSError):
    CORE_AVAILABLE = False


class TestSetupLoggingDetailed:
    """Detailed tests for setup_logging function."""
    
    @pytest.fixture
    def temp_folder(self):
        """Create temporary folder."""
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    @pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
    def test_setup_logging_returns_path(self, temp_folder):
        """Test that setup_logging returns log file path."""
        log_file = setup_logging(temp_folder)
        
        assert isinstance(log_file, str)
        assert 'log' in log_file.lower()
    
    @pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
    def test_setup_logging_creates_file(self, temp_folder):
        """Test that setup_logging creates log file."""
        log_file = setup_logging(temp_folder)
        
        # Log something to ensure file is created
        logging.info("Test log message")
        
        assert os.path.exists(log_file)


class TestPrintFunctions:
    """Tests for print utility functions."""
    
    @pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
    def test_print_memory_usage_returns_float(self):
        """Test that print_memory_usage returns float."""
        mem_mb = print_memory_usage()
        
        assert isinstance(mem_mb, float)
        assert mem_mb > 0
    
    @pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
    def test_print_section_header_no_error(self):
        """Test that print_section_header doesn't raise error."""
        # Should not raise exception
        print_section_header("Test Section Header")
        print_section_header("Another Section")
        print_section_header("")
    
    @pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
    def test_print_step_no_error(self):
        """Test that print_step doesn't raise error."""
        # Should not raise exception
        print_step(1, 10, "First Step")
        print_step(5, 10, "Middle Step")
        print_step(10, 10, "Last Step")


class TestCreateTemplateDirectory:
    """Tests for create_template_directory function."""
    
    @pytest.fixture
    def temp_cwd(self):
        """Change to temporary directory."""
        original_cwd = os.getcwd()
        temp_dir = tempfile.mkdtemp()
        os.chdir(temp_dir)
        yield temp_dir
        os.chdir(original_cwd)
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    @pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
    def test_create_template_directory_creates_folder(self, temp_cwd):
        """Test that create_template_directory creates templates folder."""
        create_template_directory()
        
        assert os.path.exists('templates')
    
    @pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
    def test_create_template_directory_creates_file(self, temp_cwd):
        """Test that create_template_directory creates template file."""
        create_template_directory()
        
        # Check if template file exists
        template_files = os.listdir('templates')
        assert len(template_files) > 0


class TestPhylogeneticCoreDetailed:
    """Detailed tests for PhylogeneticCore class."""
    
    @pytest.fixture
    def temp_tree_file(self):
        """Create temporary tree file."""
        temp_dir = tempfile.mkdtemp()
        tree_path = os.path.join(temp_dir, 'test.newick')
        with open(tree_path, 'w') as f:
            f.write('((A:0.1,B:0.2):0.3,(C:0.15,D:0.25):0.35);')
        yield tree_path
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    @pytest.mark.skipif(not CORE_AVAILABLE or not BIO_AVAILABLE, reason="Core or Bio not available")
    def test_load_tree_returns_tree(self, temp_tree_file):
        """Test that load_tree returns a tree object."""
        tree = PhylogeneticCore.load_tree(temp_tree_file)
        
        assert tree is not None
    
    @pytest.mark.skipif(not CORE_AVAILABLE or not BIO_AVAILABLE, reason="Core or Bio not available")
    def test_load_tree_correct_terminals(self, temp_tree_file):
        """Test that loaded tree has correct terminals."""
        tree = PhylogeneticCore.load_tree(temp_tree_file)
        
        terminals = tree.get_terminals()
        terminal_names = [t.name for t in terminals]
        
        assert len(terminals) == 4
        assert 'A' in terminal_names
        assert 'B' in terminal_names
        assert 'C' in terminal_names
        assert 'D' in terminal_names


class TestConfigDetailed:
    """Detailed tests for Config class."""
    
    @pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
    def test_config_initialization(self):
        """Test Config initialization."""
        config = Config()
        
        assert config is not None
    
    @pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
    def test_config_has_attributes(self):
        """Test that Config has expected attributes."""
        config = Config()
        
        # Should have some attributes
        assert hasattr(config, '__dict__')
    
    @pytest.mark.skipif(not CORE_AVAILABLE, reason="Core not available")
    def test_config_can_set_attributes(self):
        """Test that Config can have custom attributes set."""
        config = Config()
        
        config.custom_param = 'test_value'
        config.another_param = 123
        
        assert config.custom_param == 'test_value'
        assert config.another_param == 123
