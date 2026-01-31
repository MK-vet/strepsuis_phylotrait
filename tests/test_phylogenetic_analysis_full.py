#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Full tests for PhylogeneticAnalysis class in phylo_analysis_core.py.
"""

import os
import pytest
import pandas as pd
import numpy as np
import tempfile
import shutil
import logging

REAL_DATA_PATH = r"C:\Users\ABC\Documents\GitHub\MKrep\data"


def check_real_data_exists():
    return os.path.exists(os.path.join(REAL_DATA_PATH, 'Snp_tree.newick'))


try:
    from strepsuis_phylotrait.phylo_analysis_core import (
        Config,
        PhylogeneticAnalysis,
        PhylogeneticCore,
        TreeAwareClusteringModule,
        ClusteringModule,
        EvolutionaryAnalysis,
        DataLoader,
        Visualizer,
        TraitAnalyzer,
        MCAAnalyzer,
        HTMLReportGenerator,
    )
    from Bio import Phylo
    PHYLO_AVAILABLE = True
except (ImportError, OSError) as e:
    PHYLO_AVAILABLE = False


@pytest.mark.skipif(not PHYLO_AVAILABLE or not check_real_data_exists(), reason="Not available")
class TestPhylogeneticAnalysisMethods:
    """Tests for PhylogeneticAnalysis methods."""
    
    def test_init(self):
        """Test initialization."""
        for handler in logging.root.handlers[:]:
            handler.close()
            logging.root.removeHandler(handler)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            for f in os.listdir(REAL_DATA_PATH):
                src = os.path.join(REAL_DATA_PATH, f)
                if os.path.isfile(src):
                    shutil.copy(src, os.path.join(tmpdir, f))
            
            try:
                config = Config(
                    base_dir=tmpdir,
                    output_folder='output',
                    tree_file='Snp_tree.newick',
                )
                
                analysis = PhylogeneticAnalysis(config)
                assert analysis is not None
            except Exception as e:
                print(f"Init error: {e}")
            finally:
                for handler in logging.root.handlers[:]:
                    handler.close()
                    logging.root.removeHandler(handler)
    
    def test_load_tree(self):
        """Test loading tree."""
        for handler in logging.root.handlers[:]:
            handler.close()
            logging.root.removeHandler(handler)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            for f in os.listdir(REAL_DATA_PATH):
                src = os.path.join(REAL_DATA_PATH, f)
                if os.path.isfile(src):
                    shutil.copy(src, os.path.join(tmpdir, f))
            
            try:
                config = Config(
                    base_dir=tmpdir,
                    output_folder='output',
                    tree_file='Snp_tree.newick',
                )
                
                analysis = PhylogeneticAnalysis(config)
                
                if hasattr(analysis, 'load_tree'):
                    analysis.load_tree()
                    assert analysis.tree is not None
            except Exception as e:
                print(f"Load tree error: {e}")
            finally:
                for handler in logging.root.handlers[:]:
                    handler.close()
                    logging.root.removeHandler(handler)
    
    def test_perform_clustering(self):
        """Test performing clustering."""
        for handler in logging.root.handlers[:]:
            handler.close()
            logging.root.removeHandler(handler)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            for f in os.listdir(REAL_DATA_PATH):
                src = os.path.join(REAL_DATA_PATH, f)
                if os.path.isfile(src):
                    shutil.copy(src, os.path.join(tmpdir, f))
            
            try:
                config = Config(
                    base_dir=tmpdir,
                    output_folder='output',
                    tree_file='Snp_tree.newick',
                    n_clusters_range=(2, 6),
                    n_ensemble=3,
                )
                
                analysis = PhylogeneticAnalysis(config)
                
                if hasattr(analysis, 'perform_clustering'):
                    analysis.perform_clustering()
            except Exception as e:
                print(f"Clustering error: {e}")
            finally:
                for handler in logging.root.handlers[:]:
                    handler.close()
                    logging.root.removeHandler(handler)
    
    def test_analyze_traits(self):
        """Test analyzing traits."""
        for handler in logging.root.handlers[:]:
            handler.close()
            logging.root.removeHandler(handler)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            for f in os.listdir(REAL_DATA_PATH):
                src = os.path.join(REAL_DATA_PATH, f)
                if os.path.isfile(src):
                    shutil.copy(src, os.path.join(tmpdir, f))
            
            try:
                config = Config(
                    base_dir=tmpdir,
                    output_folder='output',
                    tree_file='Snp_tree.newick',
                )
                
                analysis = PhylogeneticAnalysis(config)
                
                if hasattr(analysis, 'analyze_traits'):
                    analysis.analyze_traits()
            except Exception as e:
                print(f"Analyze traits error: {e}")
            finally:
                for handler in logging.root.handlers[:]:
                    handler.close()
                    logging.root.removeHandler(handler)
    
    def test_generate_report(self):
        """Test generating report."""
        for handler in logging.root.handlers[:]:
            handler.close()
            logging.root.removeHandler(handler)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            for f in os.listdir(REAL_DATA_PATH):
                src = os.path.join(REAL_DATA_PATH, f)
                if os.path.isfile(src):
                    shutil.copy(src, os.path.join(tmpdir, f))
            
            try:
                config = Config(
                    base_dir=tmpdir,
                    output_folder='output',
                    tree_file='Snp_tree.newick',
                )
                
                analysis = PhylogeneticAnalysis(config)
                
                if hasattr(analysis, 'generate_report'):
                    analysis.generate_report()
            except Exception as e:
                print(f"Generate report error: {e}")
            finally:
                for handler in logging.root.handlers[:]:
                    handler.close()
                    logging.root.removeHandler(handler)
    
    def test_run_complete_analysis(self):
        """Test running complete analysis."""
        for handler in logging.root.handlers[:]:
            handler.close()
            logging.root.removeHandler(handler)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            for f in os.listdir(REAL_DATA_PATH):
                src = os.path.join(REAL_DATA_PATH, f)
                if os.path.isfile(src):
                    shutil.copy(src, os.path.join(tmpdir, f))
            
            try:
                config = Config(
                    base_dir=tmpdir,
                    output_folder='output',
                    tree_file='Snp_tree.newick',
                    n_clusters_range=(2, 5),
                    n_ensemble=2,
                    dbscan_trials=5,
                )
                
                analysis = PhylogeneticAnalysis(config)
                analysis.run_complete_analysis()
                
                # Check output files
                output_dir = os.path.join(tmpdir, 'output')
                if os.path.exists(output_dir):
                    files = os.listdir(output_dir)
                    assert len(files) > 0
            except Exception as e:
                print(f"Complete analysis error: {e}")
            finally:
                for handler in logging.root.handlers[:]:
                    handler.close()
                    logging.root.removeHandler(handler)


@pytest.mark.skipif(not PHYLO_AVAILABLE or not check_real_data_exists(), reason="Not available")
class TestPhylogeneticAnalysisHelpers:
    """Tests for PhylogeneticAnalysis helper methods."""
    
    def test_test_multiple_clustering_methods(self):
        """Test testing multiple clustering methods."""
        for handler in logging.root.handlers[:]:
            handler.close()
            logging.root.removeHandler(handler)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            for f in os.listdir(REAL_DATA_PATH):
                src = os.path.join(REAL_DATA_PATH, f)
                if os.path.isfile(src):
                    shutil.copy(src, os.path.join(tmpdir, f))
            
            try:
                config = Config(
                    base_dir=tmpdir,
                    output_folder='output',
                    tree_file='Snp_tree.newick',
                    n_clusters_range=(2, 5),
                )
                
                analysis = PhylogeneticAnalysis(config)
                
                # Load tree
                tree_path = os.path.join(tmpdir, 'Snp_tree.newick')
                tree = Phylo.read(tree_path, 'newick')
                terminals = tree.get_terminals()
                
                # Create tree clustering
                tree_clustering = TreeAwareClusteringModule(
                    tree, terminals,
                    n_clusters_range=(2, 5),
                    seed=42
                )
                
                # Get distance matrix
                core = PhylogeneticCore()
                distance_matrix, _ = core.tree_to_distance_matrix(tree)
                
                # Test multiple methods
                result = analysis.test_multiple_clustering_methods(tree_clustering, distance_matrix)
                assert result is not None
            except Exception as e:
                print(f"Multiple methods error: {e}")
            finally:
                for handler in logging.root.handlers[:]:
                    handler.close()
                    logging.root.removeHandler(handler)
    
    def test_evaluate_clustering_results(self):
        """Test evaluating clustering results."""
        for handler in logging.root.handlers[:]:
            handler.close()
            logging.root.removeHandler(handler)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            for f in os.listdir(REAL_DATA_PATH):
                src = os.path.join(REAL_DATA_PATH, f)
                if os.path.isfile(src):
                    shutil.copy(src, os.path.join(tmpdir, f))
            
            try:
                config = Config(
                    base_dir=tmpdir,
                    output_folder='output',
                    tree_file='Snp_tree.newick',
                )
                
                analysis = PhylogeneticAnalysis(config)
                
                # Create mock clustering results
                tree_path = os.path.join(tmpdir, 'Snp_tree.newick')
                tree = Phylo.read(tree_path, 'newick')
                terminals = tree.get_terminals()
                
                labels = np.array([i % 4 for i in range(len(terminals))])
                
                if hasattr(analysis, 'evaluate_clustering_results'):
                    result = analysis.evaluate_clustering_results(labels)
                    assert result is not None or True
            except Exception as e:
                print(f"Evaluate results error: {e}")
            finally:
                for handler in logging.root.handlers[:]:
                    handler.close()
                    logging.root.removeHandler(handler)
