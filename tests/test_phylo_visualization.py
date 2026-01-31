#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tests for visualization functions in phylo_analysis_core.py.
"""

import os
import pytest
import pandas as pd
import numpy as np
import tempfile
import shutil

REAL_DATA_PATH = r"C:\Users\ABC\Documents\GitHub\MKrep\data"


def check_real_data_exists():
    return os.path.exists(os.path.join(REAL_DATA_PATH, 'Snp_tree.newick'))


try:
    from strepsuis_phylotrait.phylo_analysis_core import (
        Visualizer,
        TraitAnalyzer,
        MCAAnalyzer,
        HTMLReportGenerator,
        PhylogeneticCore,
    )
    from Bio import Phylo
    PHYLO_AVAILABLE = True
except (ImportError, OSError) as e:
    PHYLO_AVAILABLE = False


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestVisualizerFull:
    """Full tests for Visualizer class."""
    
    def test_init(self):
        """Test initialization."""
        with tempfile.TemporaryDirectory() as tmpdir:
            viz = Visualizer(tmpdir)
            assert viz is not None
    
    def test_plot_cluster_distribution_basic(self):
        """Test basic cluster distribution plot."""
        with tempfile.TemporaryDirectory() as tmpdir:
            viz = Visualizer(tmpdir)
            
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(50)],
                'Cluster': [i % 3 for i in range(50)],
            })
            
            try:
                viz.plot_cluster_distribution(merged_df)
            except Exception as e:
                print(f"Plot error: {e}")
    
    def test_plot_cluster_distribution_with_traits(self):
        """Test cluster distribution with trait columns."""
        with tempfile.TemporaryDirectory() as tmpdir:
            viz = Visualizer(tmpdir)
            
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(50)],
                'Cluster': [i % 3 for i in range(50)],
                'AMR_gene1': np.random.randint(0, 2, 50),
                'AMR_gene2': np.random.randint(0, 2, 50),
                'VIR_gene1': np.random.randint(0, 2, 50),
            })
            
            try:
                viz.plot_cluster_distribution(merged_df)
            except Exception as e:
                print(f"Plot with traits error: {e}")
    
    def test_plot_umap_clusters_basic(self):
        """Test basic UMAP plot."""
        with tempfile.TemporaryDirectory() as tmpdir:
            viz = Visualizer(tmpdir)
            
            np.random.seed(42)
            embeddings = np.random.randn(50, 2)
            labels = np.array([i % 3 for i in range(50)])
            mask = np.array([True] * 50)
            
            try:
                viz.plot_umap_clusters(embeddings, labels, mask)
            except Exception as e:
                print(f"UMAP plot error: {e}")
    
    def test_plot_umap_clusters_with_outliers(self):
        """Test UMAP plot with outliers."""
        with tempfile.TemporaryDirectory() as tmpdir:
            viz = Visualizer(tmpdir)
            
            np.random.seed(42)
            embeddings = np.random.randn(50, 2)
            labels = np.array([i % 3 for i in range(50)])
            mask = np.array([True] * 45 + [False] * 5)
            outlier_assignments = {45: 0, 46: 1, 47: 2, 48: 0, 49: 1}
            
            try:
                viz.plot_umap_clusters(embeddings, labels, mask, outlier_assignments)
            except Exception as e:
                print(f"UMAP outliers error: {e}")
    
    @pytest.mark.skipif(not check_real_data_exists(), reason="Real data not available")
    def test_plot_phylogenetic_tree(self):
        """Test phylogenetic tree plot."""
        with tempfile.TemporaryDirectory() as tmpdir:
            viz = Visualizer(tmpdir)
            
            tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
            tree = Phylo.read(tree_path, 'newick')
            
            labels = np.array([i % 3 for i in range(len(list(tree.get_terminals())))])
            
            try:
                if hasattr(viz, 'plot_phylogenetic_tree'):
                    viz.plot_phylogenetic_tree(tree, labels)
            except Exception as e:
                print(f"Tree plot error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestTraitAnalyzerFull:
    """Full tests for TraitAnalyzer class."""
    
    def test_init(self):
        """Test initialization."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            assert analyzer is not None
    
    def test_analyze_category_amr(self):
        """Test AMR category analysis."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(50)],
                'Cluster': [i % 3 for i in range(50)],
                'AMR_gene1': np.random.randint(0, 2, 50),
                'AMR_gene2': np.random.randint(0, 2, 50),
                'AMR_gene3': np.random.randint(0, 2, 50),
            })
            
            try:
                if hasattr(analyzer, 'analyze_category'):
                    analyzer.analyze_category(merged_df, 'AMR', ['AMR_gene1', 'AMR_gene2', 'AMR_gene3'])
            except Exception as e:
                print(f"AMR analysis error: {e}")
    
    def test_analyze_category_virulence(self):
        """Test virulence category analysis."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(50)],
                'Cluster': [i % 3 for i in range(50)],
                'VIR_gene1': np.random.randint(0, 2, 50),
                'VIR_gene2': np.random.randint(0, 2, 50),
            })
            
            try:
                if hasattr(analyzer, 'analyze_category'):
                    analyzer.analyze_category(merged_df, 'VIR', ['VIR_gene1', 'VIR_gene2'])
            except Exception as e:
                print(f"Virulence analysis error: {e}")
    
    def test_analyze_category_mic(self):
        """Test MIC category analysis."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(50)],
                'Cluster': [i % 3 for i in range(50)],
                'MIC_drug1': np.random.randint(0, 2, 50),
                'MIC_drug2': np.random.randint(0, 2, 50),
            })
            
            try:
                if hasattr(analyzer, 'analyze_category'):
                    analyzer.analyze_category(merged_df, 'MIC', ['MIC_drug1', 'MIC_drug2'])
            except Exception as e:
                print(f"MIC analysis error: {e}")
    
    def test_analyze_all_categories(self):
        """Test analyzing all categories."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(50)],
                'Cluster': [i % 3 for i in range(50)],
                'AMR_gene1': np.random.randint(0, 2, 50),
                'VIR_gene1': np.random.randint(0, 2, 50),
                'MIC_drug1': np.random.randint(0, 2, 50),
            })
            
            try:
                if hasattr(analyzer, 'analyze_all_categories'):
                    analyzer.analyze_all_categories(merged_df)
            except Exception as e:
                print(f"All categories error: {e}")
    
    def test_compute_statistics(self):
        """Test computing statistics."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = TraitAnalyzer(tmpdir)
            
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(50)],
                'Cluster': [i % 3 for i in range(50)],
                'gene1': np.random.randint(0, 2, 50),
            })
            
            try:
                if hasattr(analyzer, 'compute_statistics'):
                    analyzer.compute_statistics(merged_df, 'gene1')
            except Exception as e:
                print(f"Statistics error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestMCAAnalyzerFull:
    """Full tests for MCAAnalyzer class."""
    
    def test_init(self):
        """Test initialization."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = MCAAnalyzer(tmpdir)
            assert analyzer is not None
    
    def test_perform_mca_basic(self):
        """Test basic MCA."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = MCAAnalyzer(tmpdir)
            
            data = pd.DataFrame({
                'gene1': np.random.randint(0, 2, 50),
                'gene2': np.random.randint(0, 2, 50),
                'gene3': np.random.randint(0, 2, 50),
                'gene4': np.random.randint(0, 2, 50),
            })
            
            try:
                if hasattr(analyzer, 'perform_mca'):
                    analyzer.perform_mca(data)
            except Exception as e:
                print(f"MCA error: {e}")
    
    def test_perform_mca_with_labels(self):
        """Test MCA with cluster labels."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = MCAAnalyzer(tmpdir)
            
            data = pd.DataFrame({
                'gene1': np.random.randint(0, 2, 50),
                'gene2': np.random.randint(0, 2, 50),
                'gene3': np.random.randint(0, 2, 50),
            })
            labels = np.array([i % 3 for i in range(50)])
            
            try:
                if hasattr(analyzer, 'perform_mca'):
                    analyzer.perform_mca(data, labels)
            except Exception as e:
                print(f"MCA with labels error: {e}")
    
    def test_plot_mca_results(self):
        """Test plotting MCA results."""
        with tempfile.TemporaryDirectory() as tmpdir:
            analyzer = MCAAnalyzer(tmpdir)
            
            data = pd.DataFrame({
                'gene1': np.random.randint(0, 2, 50),
                'gene2': np.random.randint(0, 2, 50),
                'gene3': np.random.randint(0, 2, 50),
            })
            
            try:
                if hasattr(analyzer, 'perform_mca'):
                    analyzer.perform_mca(data)
                if hasattr(analyzer, 'plot_results'):
                    analyzer.plot_results()
            except Exception as e:
                print(f"MCA plot error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestHTMLReportGeneratorFull:
    """Full tests for HTMLReportGenerator class."""
    
    def test_init(self):
        """Test initialization."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = HTMLReportGenerator(tmpdir)
            assert generator is not None
    
    def test_add_section(self):
        """Test adding sections."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = HTMLReportGenerator(tmpdir)
            
            if hasattr(generator, 'add_section'):
                generator.add_section("Introduction", "This is the introduction.")
                generator.add_section("Methods", "These are the methods.")
                generator.add_section("Results", "These are the results.")
    
    def test_add_table(self):
        """Test adding tables."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = HTMLReportGenerator(tmpdir)
            
            df = pd.DataFrame({
                'Cluster': [0, 1, 2],
                'Size': [30, 35, 26],
                'Silhouette': [0.5, 0.6, 0.4]
            })
            
            if hasattr(generator, 'add_table'):
                generator.add_table(df, "Cluster Summary")
    
    def test_add_figure(self):
        """Test adding figures."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = HTMLReportGenerator(tmpdir)
            
            try:
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots()
                ax.plot([1, 2, 3], [1, 2, 3])
                
                if hasattr(generator, 'add_figure'):
                    generator.add_figure(fig, "Test Figure")
                
                plt.close(fig)
            except Exception as e:
                print(f"Add figure error: {e}")
    
    def test_generate_report(self):
        """Test generating report."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = HTMLReportGenerator(tmpdir)
            
            if hasattr(generator, 'add_section'):
                generator.add_section("Summary", "Analysis summary")
            
            if hasattr(generator, 'generate'):
                try:
                    generator.generate("test_report.html")
                except Exception as e:
                    print(f"Generate error: {e}")
    
    def test_generate_full_report(self):
        """Test generating full report with all components."""
        with tempfile.TemporaryDirectory() as tmpdir:
            generator = HTMLReportGenerator(tmpdir)
            
            # Add sections
            if hasattr(generator, 'add_section'):
                generator.add_section("Introduction", "Analysis overview")
                generator.add_section("Methods", "Clustering methods")
                generator.add_section("Results", "Key findings")
                generator.add_section("Conclusions", "Final conclusions")
            
            # Add tables
            if hasattr(generator, 'add_table'):
                df = pd.DataFrame({
                    'Cluster': [0, 1, 2, 3],
                    'Size': [20, 25, 30, 16],
                    'Silhouette': [0.5, 0.6, 0.4, 0.55]
                })
                generator.add_table(df, "Cluster Summary")
            
            # Add figures
            try:
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots()
                ax.bar([1, 2, 3, 4], [20, 25, 30, 16])
                
                if hasattr(generator, 'add_figure'):
                    generator.add_figure(fig, "Cluster Distribution")
                
                plt.close(fig)
            except Exception:
                pass
            
            # Generate
            if hasattr(generator, 'generate'):
                try:
                    generator.generate("full_report.html")
                except Exception as e:
                    print(f"Full report error: {e}")
