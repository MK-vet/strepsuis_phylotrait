#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Full tests for config.py.
"""

import os
import pytest
import tempfile

try:
    from strepsuis_phylotrait.config import Config
    CONFIG_AVAILABLE = True
except (ImportError, OSError) as e:
    CONFIG_AVAILABLE = False


@pytest.mark.skipif(not CONFIG_AVAILABLE, reason="Not available")
class TestConfigFull:
    """Full tests for Config class."""
    
    def test_default_config(self):
        """Test default configuration."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config = Config(base_dir=tmpdir)
            assert config is not None
            assert config.base_dir == tmpdir
    
    def test_custom_config(self):
        """Test custom configuration."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config = Config(
                base_dir=tmpdir,
                output_folder='custom_output',
                tree_file='custom_tree.newick',
                n_clusters_range=(3, 10),
                n_ensemble=10,
                dbscan_trials=20,
            )
            assert config.output_folder == 'custom_output'
            assert config.tree_file == 'custom_tree.newick'
            assert config.n_clusters_range == (3, 10)
            assert config.n_ensemble == 10
            assert config.dbscan_trials == 20
    
    def test_config_with_all_files(self):
        """Test configuration with all file paths."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config = Config(
                base_dir=tmpdir,
                output_folder='output',
                tree_file='tree.newick',
                mic_file='MIC.csv',
                amr_genes_file='AMR.csv',
                virulence_genes_file='VIR.csv',
            )
            assert config.mic_file == 'MIC.csv'
            assert config.amr_genes_file == 'AMR.csv'
            assert config.virulence_genes_file == 'VIR.csv'
    
    def test_config_paths(self):
        """Test configuration path methods."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config = Config(
                base_dir=tmpdir,
                output_folder='output',
                tree_file='tree.newick',
            )
            
            if hasattr(config, 'get_tree_path'):
                tree_path = config.get_tree_path()
                assert 'tree.newick' in tree_path
            
            if hasattr(config, 'get_output_path'):
                output_path = config.get_output_path()
                assert 'output' in output_path
    
    def test_config_validation(self):
        """Test configuration validation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config = Config(
                base_dir=tmpdir,
                n_clusters_range=(2, 8),
            )
            
            # Validate range
            assert config.n_clusters_range[0] < config.n_clusters_range[1]
    
    def test_config_to_dict(self):
        """Test configuration to dictionary."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config = Config(
                base_dir=tmpdir,
                output_folder='output',
            )
            
            if hasattr(config, 'to_dict'):
                config_dict = config.to_dict()
                assert 'base_dir' in config_dict
                assert 'output_folder' in config_dict
    
    def test_config_from_dict(self):
        """Test configuration from dictionary."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config_dict = {
                'base_dir': tmpdir,
                'output_folder': 'output',
                'n_clusters_range': (2, 8),
            }
            
            if hasattr(Config, 'from_dict'):
                config = Config.from_dict(config_dict)
                assert config.base_dir == tmpdir
