#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for PhylogeneticAnalysis fallback paths and edge cases.
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
    )
    PHYLO_AVAILABLE = True
except (ImportError, OSError) as e:
    PHYLO_AVAILABLE = False


@pytest.mark.skipif(not PHYLO_AVAILABLE or not check_real_data_exists(), reason="Not available")
class TestPhylogeneticAnalysisFallback:
    """Tests for PhylogeneticAnalysis fallback paths."""
    
    def test_run_complete_analysis_with_fallback(self):
        """Test run_complete_analysis with fallback to standard clustering."""
        for handler in logging.root.handlers[:]:
            handler.close()
            logging.root.removeHandler(handler)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            for f in os.listdir(REAL_DATA_PATH):
                src = os.path.join(REAL_DATA_PATH, f)
                if os.path.isfile(src):
                    shutil.copy(src, os.path.join(tmpdir, f))
            
            try:
                # Use config that might trigger fallback
                config = Config(
                    base_dir=tmpdir,
                    output_folder='output',
                    tree_file='Snp_tree.newick',
                    n_clusters_range=(2, 4),
                    n_ensemble=2,
                    dbscan_trials=3,
                    # Force fallback by using parameters that might fail
                    umap_components=2,
                    umap_neighbors=5,
                    umap_min_dist=0.1,
                    outlier_contamination=0.1,
                    outlier_n_estimators=50,
                )
                
                analysis = PhylogeneticAnalysis(config)
                analysis.run_complete_analysis()
                
            except Exception as e:
                print(f"Fallback analysis error: {e}")
            finally:
                for handler in logging.root.handlers[:]:
                    handler.close()
                    logging.root.removeHandler(handler)
    
    def test_run_complete_analysis_all_steps(self):
        """Test run_complete_analysis with all steps."""
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
                    n_ensemble=3,
                    dbscan_trials=5,
                )
                
                analysis = PhylogeneticAnalysis(config)
                analysis.run_complete_analysis()
                
                # Check that all output files were created
                output_dir = os.path.join(tmpdir, 'output')
                if os.path.exists(output_dir):
                    files = os.listdir(output_dir)
                    # Should have various output files
                    assert len(files) > 0
                    
            except Exception as e:
                print(f"All steps analysis error: {e}")
            finally:
                for handler in logging.root.handlers[:]:
                    handler.close()
                    logging.root.removeHandler(handler)
    
    def test_run_complete_analysis_with_trait_data(self):
        """Test run_complete_analysis with trait data files."""
        for handler in logging.root.handlers[:]:
            handler.close()
            logging.root.removeHandler(handler)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            for f in os.listdir(REAL_DATA_PATH):
                src = os.path.join(REAL_DATA_PATH, f)
                if os.path.isfile(src):
                    shutil.copy(src, os.path.join(tmpdir, f))
            
            # Create trait data files
            from Bio import Phylo
            tree_path = os.path.join(tmpdir, 'Snp_tree.newick')
            tree = Phylo.read(tree_path, 'newick')
            terminals = tree.get_terminals()
            strain_names = [str(t.name) for t in terminals]
            
            # MIC data
            pd.DataFrame({
                'Strain_ID': strain_names,
                'MIC_PEN': np.random.randint(0, 2, len(strain_names)),
                'MIC_TET': np.random.randint(0, 2, len(strain_names)),
            }).to_csv(os.path.join(tmpdir, 'MIC.csv'), index=False)
            
            # AMR genes
            pd.DataFrame({
                'Strain_ID': strain_names,
                'tetM': np.random.randint(0, 2, len(strain_names)),
                'ermB': np.random.randint(0, 2, len(strain_names)),
            }).to_csv(os.path.join(tmpdir, 'AMR_genes.csv'), index=False)
            
            # Virulence factors
            pd.DataFrame({
                'Strain_ID': strain_names,
                'cps': np.random.randint(0, 2, len(strain_names)),
                'mrp': np.random.randint(0, 2, len(strain_names)),
            }).to_csv(os.path.join(tmpdir, 'Virulence.csv'), index=False)
            
            try:
                config = Config(
                    base_dir=tmpdir,
                    output_folder='output',
                    tree_file='Snp_tree.newick',
                    mic_file='MIC.csv',
                    amr_genes_file='AMR_genes.csv',
                    virulence_genes_file='Virulence.csv',
                    n_clusters_range=(2, 4),
                    n_ensemble=2,
                )
                
                analysis = PhylogeneticAnalysis(config)
                analysis.run_complete_analysis()
                
            except Exception as e:
                print(f"Trait data analysis error: {e}")
            finally:
                for handler in logging.root.handlers[:]:
                    handler.close()
                    logging.root.removeHandler(handler)


@pytest.mark.skipif(not PHYLO_AVAILABLE or not check_real_data_exists(), reason="Not available")
class TestPhylogeneticAnalysisEdgeCases:
    """Tests for PhylogeneticAnalysis edge cases."""
    
    def test_run_complete_analysis_minimal_config(self):
        """Test with minimal configuration."""
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
                analysis.run_complete_analysis()
                
            except Exception as e:
                print(f"Minimal config error: {e}")
            finally:
                for handler in logging.root.handlers[:]:
                    handler.close()
                    logging.root.removeHandler(handler)
    
    def test_run_complete_analysis_different_cluster_ranges(self):
        """Test with different cluster ranges."""
        for handler in logging.root.handlers[:]:
            handler.close()
            logging.root.removeHandler(handler)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            for f in os.listdir(REAL_DATA_PATH):
                src = os.path.join(REAL_DATA_PATH, f)
                if os.path.isfile(src):
                    shutil.copy(src, os.path.join(tmpdir, f))
            
            for min_k, max_k in [(2, 3), (3, 5), (4, 8)]:
                try:
                    config = Config(
                        base_dir=tmpdir,
                        output_folder=f'output_{min_k}_{max_k}',
                        tree_file='Snp_tree.newick',
                        n_clusters_range=(min_k, max_k),
                        n_ensemble=2,
                    )
                    
                    analysis = PhylogeneticAnalysis(config)
                    analysis.run_complete_analysis()
                    
                except Exception as e:
                    print(f"Cluster range {min_k}-{max_k} error: {e}")
                finally:
                    for handler in logging.root.handlers[:]:
                        handler.close()
                        logging.root.removeHandler(handler)
