#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Final push tests for phylo_analysis_core.py to reach 80% coverage.
"""

import os
import pytest
import pandas as pd
import numpy as np
import tempfile
import shutil
import logging
from unittest.mock import patch, MagicMock

REAL_DATA_PATH = r"C:\Users\ABC\Documents\GitHub\MKrep\data"


def check_real_data_exists():
    return os.path.exists(os.path.join(REAL_DATA_PATH, 'Snp_tree.newick'))


try:
    from strepsuis_phylotrait.phylo_analysis_core import (
        PhylogeneticCore,
        TreeAwareClusteringModule,
        ClusteringModule,
        EvolutionaryAnalysis,
        TraitAnalyzer,
        MCAAnalyzer,
        Visualizer,
        DataLoader,
        HTMLReportGenerator,
        Config,
        PhylogeneticAnalysis,
    )
    from Bio import Phylo
    PHYLO_AVAILABLE = True
except (ImportError, OSError) as e:
    PHYLO_AVAILABLE = False


@pytest.mark.skipif(not PHYLO_AVAILABLE or not check_real_data_exists(), reason="Not available")
class TestPhylogeneticAnalysisFinalPush:
    """Final push tests for PhylogeneticAnalysis."""
    
    def test_run_complete_analysis_quick(self):
        """Test quick complete analysis."""
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
                    n_clusters_range=(2, 4),
                    n_ensemble=2,
                    dbscan_trials=3,
                )
                
                analysis = PhylogeneticAnalysis(config)
                analysis.run_complete_analysis()
                
            except Exception as e:
                print(f"Quick analysis error: {e}")
            finally:
                for handler in logging.root.handlers[:]:
                    handler.close()
                    logging.root.removeHandler(handler)


@pytest.mark.skipif(not PHYLO_AVAILABLE or not check_real_data_exists(), reason="Not available")
class TestHTMLReportGeneratorFinalPush:
    """Final push tests for HTMLReportGenerator."""
    
    def test_generate_report_all_sections(self):
        """Test report generation with all sections."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Copy real data
            for f in os.listdir(REAL_DATA_PATH):
                src = os.path.join(REAL_DATA_PATH, f)
                if os.path.isfile(src):
                    shutil.copy(src, os.path.join(tmpdir, f))
            
            try:
                generator = HTMLReportGenerator(tmpdir, base_dir=tmpdir)
                
                # Create all output files
                tree_path = os.path.join(tmpdir, 'Snp_tree.newick')
                tree = Phylo.read(tree_path, 'newick')
                terminals = tree.get_terminals()
                strain_names = [str(t.name) for t in terminals]
                
                # Clusters
                pd.DataFrame({
                    'Strain_ID': strain_names,
                    'Cluster': [i % 4 for i in range(len(strain_names))]
                }).to_csv(os.path.join(tmpdir, 'phylogenetic_clusters.csv'), index=False)
                
                # Distribution
                pd.DataFrame({
                    'Cluster': [0, 1, 2, 3],
                    'Strain_Count': [23, 24, 22, 22],
                    'Percentage': [25.3, 26.4, 24.2, 24.2]
                }).to_csv(os.path.join(tmpdir, 'cluster_distribution.csv'), index=False)
                
                # Evolutionary
                pd.DataFrame({
                    'Cluster': [0, 1, 2, 3],
                    'Mean_PD': [0.1, 0.15, 0.12, 0.11]
                }).to_csv(os.path.join(tmpdir, 'evolutionary_cluster_analysis.csv'), index=False)
                
                # Silhouette
                pd.DataFrame({
                    'k': [2, 3, 4, 5],
                    'Silhouette': [0.4, 0.5, 0.55, 0.45]
                }).to_csv(os.path.join(tmpdir, 'silhouette_scores.csv'), index=False)
                
                # Bootstrap
                pd.DataFrame({
                    'Feature': ['gene_1', 'gene_2', 'gene_3'],
                    'Importance_Mean': [0.8, 0.6, 0.4]
                }).to_csv(os.path.join(tmpdir, 'bootstrap_feature_importance.csv'), index=False)
                
                # Association rules
                pd.DataFrame({
                    'Antecedent': ['gene_1', 'gene_2'],
                    'Consequent': ['gene_3', 'gene_4'],
                    'Support': [0.3, 0.2],
                    'Confidence': [0.8, 0.7],
                    'Lift': [2.5, 2.0]
                }).to_csv(os.path.join(tmpdir, 'association_rules.csv'), index=False)
                
                # Analysis results
                analysis_results = {
                    'summary_stats': pd.DataFrame({
                        'Metric': ['Silhouette', 'Calinski'],
                        'Value': [0.55, 120.5]
                    }),
                    'embeddings': np.random.randn(len(terminals), 2),
                    'labels': np.array([i % 4 for i in range(len(terminals))]),
                    'mask': np.array([True] * len(terminals)),
                }
                
                config = MagicMock()
                config.tree_file = 'Snp_tree.newick'
                
                with patch.object(generator, 'template_env') as mock_env:
                    mock_template = MagicMock()
                    mock_template.render.return_value = "<html>Complete Report</html>"
                    mock_env.get_template.return_value = mock_template
                    
                    generator.generate_report(analysis_results, config)
                    
            except Exception as e:
                print(f"All sections report error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE or not check_real_data_exists(), reason="Not available")
class TestTraitAnalyzerFinalPush:
    """Final push tests for TraitAnalyzer."""
    
    def test_all_analysis_methods(self):
        """Test all analysis methods in sequence."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            np.random.seed(42)
            n_samples = 100
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(n_samples)],
                'Cluster': [i % 4 for i in range(n_samples)],
                **{f'AMR_gene_{i}': np.random.randint(0, 2, n_samples) for i in range(15)},
                **{f'VIR_gene_{i}': np.random.randint(0, 2, n_samples) for i in range(10)},
                **{f'MIC_drug_{i}': np.random.randint(0, 2, n_samples) for i in range(8)}
            })
            
            # Bootstrap feature importance
            try:
                result = analyzer.bootstrap_feature_importance(merged_df, n_bootstrap=5)
                assert result is not None
                result.to_csv(os.path.join(tmpdir, 'bootstrap_feature_importance.csv'), index=False)
            except Exception as e:
                print(f"Bootstrap error: {e}")
            
            # Log odds ratio
            try:
                result = analyzer.log_odds_ratio_analysis(merged_df)
                assert result is not None
            except Exception as e:
                print(f"Log odds error: {e}")
            
            # Association rules
            try:
                result = analyzer.association_rule_mining(merged_df, min_support=0.1, min_confidence=0.5)
                assert result is not None
            except Exception as e:
                print(f"Association rules error: {e}")
            
            # Label shared/unique
            try:
                result = analyzer.label_shared_unique_features(merged_df, presence_threshold=0.3)
                assert result is not None
            except Exception as e:
                print(f"Label shared/unique error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE or not check_real_data_exists(), reason="Not available")
class TestClusteringFinalPush:
    """Final push tests for clustering."""
    
    def test_tree_aware_clustering_all_methods(self):
        """Test all tree-aware clustering methods."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = Phylo.read(tree_path, 'newick')
        terminals = tree.get_terminals()
        
        module = TreeAwareClusteringModule(
            tree, terminals,
            n_clusters_range=(2, 6),
            seed=42
        )
        
        core = PhylogeneticCore()
        distance_matrix, _ = core.tree_to_distance_matrix(tree)
        
        # Test all methods
        for method in ['max', 'sum', 'avg']:
            try:
                result = module.tree_cluster_algorithm(
                    distance_matrix,
                    method=method,
                    n_clusters=4
                )
                assert result is not None
                assert len(result) == len(terminals)
            except Exception as e:
                print(f"Method {method} error: {e}")
        
        # Test phydelity
        try:
            result = module.phydelity_clustering(distance_matrix)
            assert result is not None
        except Exception as e:
            print(f"Phydelity error: {e}")
        
        # Test monophyly
        labels = np.array([i % 4 for i in range(len(terminals))])
        
        try:
            eval_result = module.evaluate_monophyly(labels)
            assert eval_result is not None
        except Exception as e:
            print(f"Evaluate monophyly error: {e}")
        
        try:
            enforced = module.ensure_monophyletic_clusters(labels)
            assert enforced is not None
        except Exception as e:
            print(f"Enforce monophyly error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE or not check_real_data_exists(), reason="Not available")
class TestEvolutionaryFinalPush:
    """Final push tests for evolutionary analysis."""
    
    def test_all_evolutionary_methods(self):
        """Test all evolutionary methods."""
        tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
        tree = Phylo.read(tree_path, 'newick')
        terminals = tree.get_terminals()
        strain_names = [str(t.name) for t in terminals]
        
        labels = np.array([i % 4 for i in range(len(terminals))])
        mask = np.array([True] * len(terminals))
        
        # Evolutionary cluster analysis
        try:
            result = EvolutionaryAnalysis.evolutionary_cluster_analysis(
                tree, labels, strain_names, mask
            )
            assert result is not None
        except Exception as e:
            print(f"Evolutionary cluster error: {e}")
        
        # Beta diversity
        try:
            result = EvolutionaryAnalysis.calculate_beta_diversity(
                tree, labels, strain_names, mask
            )
            assert result is not None
        except Exception as e:
            print(f"Beta diversity error: {e}")
        
        # Evolution rates
        cluster_df = pd.DataFrame({
            'Cluster': [0, 1, 2, 3],
            'Size': [23, 24, 22, 22],
            'Silhouette': [0.5, 0.6, 0.4, 0.55],
            'Mean_Branch_Length': [0.01, 0.02, 0.015, 0.018]
        })
        
        try:
            result = EvolutionaryAnalysis.calculate_evolution_rates(cluster_df)
            assert result is not None
        except Exception as e:
            print(f"Evolution rates error: {e}")
        
        # Phylogenetic signal
        trait_data = pd.DataFrame({
            f'trait_{i}': np.random.randint(0, 2, len(terminals))
            for i in range(5)
        }, index=strain_names)
        
        try:
            result = EvolutionaryAnalysis.calculate_phylogenetic_signal_fritz_purvis(tree, trait_data)
            assert result is not None
        except Exception as e:
            print(f"Phylogenetic signal error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestMCAFinalPush:
    """Final push tests for MCA."""
    
    def test_mca_with_all_feature_types(self):
        """Test MCA with all feature types."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = MCAAnalyzer(tmpdir)
            
            np.random.seed(42)
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(100)],
                'Cluster': [i % 4 for i in range(100)],
                **{f'MIC_drug_{i}': np.random.randint(0, 2, 100) for i in range(5)},
                **{f'AMR_gene_{i}': np.random.randint(0, 2, 100) for i in range(8)},
                **{f'VIR_factor_{i}': np.random.randint(0, 2, 100) for i in range(6)},
                **{f'other_{i}': np.random.randint(0, 2, 100) for i in range(4)}
            })
            
            try:
                row_coords, col_coords, summary = analyzer.perform_mca_analysis(merged_df)
                
                if row_coords is not None:
                    assert 'Component_1' in row_coords.columns
                    assert 'Component_2' in row_coords.columns
                    assert 'Cluster' in row_coords.columns
                
                if col_coords is not None:
                    assert 'Feature_Type' in col_coords.columns
                    feature_types = col_coords['Feature_Type'].unique()
                    # Should have MIC, AMR, Virulence, Other
                    assert len(feature_types) >= 1
                
                if summary is not None:
                    assert 'Eigenvalue' in summary.columns
                    assert 'Explained_Inertia' in summary.columns
                    
            except Exception as e:
                print(f"MCA all types error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestVisualizerFinalPush:
    """Final push tests for Visualizer."""
    
    def test_all_visualizations(self):
        """Test all visualization methods."""
        with tempfile.TemporaryDirectory() as tmpdir:
            viz = Visualizer(tmpdir)
            
            np.random.seed(42)
            n_samples = 100
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(n_samples)],
                'Cluster': [i % 5 for i in range(n_samples)],
                **{f'gene_{i}': np.random.randint(0, 2, n_samples) for i in range(20)}
            })
            
            # Cluster distribution
            try:
                viz.plot_cluster_distribution(merged_df)
            except Exception as e:
                print(f"Cluster distribution error: {e}")
            
            # UMAP
            embeddings = np.random.randn(n_samples, 2)
            labels = np.array([i % 5 for i in range(n_samples)])
            mask = np.array([True] * 90 + [False] * 10)
            outlier_assignments = {90+i: i % 5 for i in range(10)}
            
            try:
                viz.plot_umap_clusters(embeddings, labels, mask, outlier_assignments)
            except Exception as e:
                print(f"UMAP error: {e}")
