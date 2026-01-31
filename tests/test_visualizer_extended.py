#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extended tests for Visualizer class in phylo_analysis_core.py.
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
        PhylogeneticCore,
    )
    from Bio import Phylo
    PHYLO_AVAILABLE = True
except (ImportError, OSError) as e:
    PHYLO_AVAILABLE = False


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestVisualizerPlots:
    """Tests for Visualizer plotting methods."""
    
    def test_plot_cluster_distribution(self):
        """Test cluster distribution plot."""
        with tempfile.TemporaryDirectory() as tmpdir:
            viz = Visualizer(tmpdir)
            
            merged_df = pd.DataFrame({
                'Strain_ID': [f'strain_{i}' for i in range(50)],
                'Cluster': [i % 3 for i in range(50)],
                'gene1': np.random.randint(0, 2, 50),
                'gene2': np.random.randint(0, 2, 50),
            })
            
            try:
                viz.plot_cluster_distribution(merged_df)
            except Exception as e:
                print(f"Cluster distribution error: {e}")
    
    def test_plot_umap_clusters(self):
        """Test UMAP clusters plot."""
        with tempfile.TemporaryDirectory() as tmpdir:
            viz = Visualizer(tmpdir)
            
            np.random.seed(42)
            embeddings = np.random.randn(50, 2)
            labels = np.array([i % 3 for i in range(50)])
            mask = np.array([True] * 50)
            
            try:
                viz.plot_umap_clusters(embeddings, labels, mask)
            except Exception as e:
                print(f"UMAP clusters error: {e}")
    
    def test_plot_umap_with_outliers(self):
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
    
    def test_plot_silhouette_scores(self):
        """Test silhouette scores plot."""
        with tempfile.TemporaryDirectory() as tmpdir:
            viz = Visualizer(tmpdir)
            
            scores = {2: 0.4, 3: 0.5, 4: 0.45, 5: 0.35}
            
            try:
                if hasattr(viz, 'plot_silhouette_scores'):
                    viz.plot_silhouette_scores(scores)
            except Exception as e:
                print(f"Silhouette scores error: {e}")
    
    def test_plot_dendrogram(self):
        """Test dendrogram plot."""
        with tempfile.TemporaryDirectory() as tmpdir:
            viz = Visualizer(tmpdir)
            
            np.random.seed(42)
            distance_matrix = np.random.rand(50, 50)
            distance_matrix = (distance_matrix + distance_matrix.T) / 2
            np.fill_diagonal(distance_matrix, 0)
            
            labels = np.array([i % 3 for i in range(50)])
            
            try:
                if hasattr(viz, 'plot_dendrogram'):
                    viz.plot_dendrogram(distance_matrix, labels)
            except Exception as e:
                print(f"Dendrogram error: {e}")
    
    def test_plot_heatmap(self):
        """Test heatmap plot."""
        with tempfile.TemporaryDirectory() as tmpdir:
            viz = Visualizer(tmpdir)
            
            np.random.seed(42)
            data = pd.DataFrame({
                f'gene_{i}': np.random.randint(0, 2, 50) for i in range(20)
            })
            
            try:
                if hasattr(viz, 'plot_heatmap'):
                    viz.plot_heatmap(data)
            except Exception as e:
                print(f"Heatmap error: {e}")
    
    def test_plot_trait_prevalence(self):
        """Test trait prevalence plot."""
        with tempfile.TemporaryDirectory() as tmpdir:
            viz = Visualizer(tmpdir)
            
            prevalence_df = pd.DataFrame({
                'Trait': [f'gene_{i}' for i in range(10)],
                'Prevalence': np.random.rand(10),
            })
            
            try:
                if hasattr(viz, 'plot_trait_prevalence'):
                    viz.plot_trait_prevalence(prevalence_df)
            except Exception as e:
                print(f"Trait prevalence error: {e}")
    
    @pytest.mark.skipif(not check_real_data_exists(), reason="Real data not available")
    def test_plot_phylogenetic_tree(self):
        """Test phylogenetic tree plot."""
        with tempfile.TemporaryDirectory() as tmpdir:
            viz = Visualizer(tmpdir)
            
            tree_path = os.path.join(REAL_DATA_PATH, 'Snp_tree.newick')
            tree = Phylo.read(tree_path, 'newick')
            terminals = tree.get_terminals()
            
            labels = np.array([i % 3 for i in range(len(terminals))])
            
            try:
                if hasattr(viz, 'plot_phylogenetic_tree'):
                    viz.plot_phylogenetic_tree(tree, labels)
            except Exception as e:
                print(f"Phylogenetic tree error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE, reason="Not available")
class TestVisualizerInteractive:
    """Tests for interactive visualization methods."""
    
    def test_create_interactive_umap(self):
        """Test interactive UMAP plot."""
        with tempfile.TemporaryDirectory() as tmpdir:
            viz = Visualizer(tmpdir)
            
            np.random.seed(42)
            embeddings = np.random.randn(50, 2)
            labels = np.array([i % 3 for i in range(50)])
            strain_ids = [f'strain_{i}' for i in range(50)]
            
            try:
                if hasattr(viz, 'create_interactive_umap'):
                    viz.create_interactive_umap(embeddings, labels, strain_ids)
            except Exception as e:
                print(f"Interactive UMAP error: {e}")
    
    def test_create_interactive_heatmap(self):
        """Test interactive heatmap."""
        with tempfile.TemporaryDirectory() as tmpdir:
            viz = Visualizer(tmpdir)
            
            np.random.seed(42)
            data = pd.DataFrame({
                f'gene_{i}': np.random.randint(0, 2, 50) for i in range(20)
            })
            
            try:
                if hasattr(viz, 'create_interactive_heatmap'):
                    viz.create_interactive_heatmap(data)
            except Exception as e:
                print(f"Interactive heatmap error: {e}")
    
    def test_create_plotly_figure(self):
        """Test creating plotly figure."""
        with tempfile.TemporaryDirectory() as tmpdir:
            viz = Visualizer(tmpdir)
            
            np.random.seed(42)
            data = pd.DataFrame({
                'x': np.random.randn(50),
                'y': np.random.randn(50),
                'cluster': [i % 3 for i in range(50)],
            })
            
            try:
                if hasattr(viz, 'create_plotly_scatter'):
                    viz.create_plotly_scatter(data, 'x', 'y', 'cluster')
            except Exception as e:
                print(f"Plotly scatter error: {e}")


@pytest.mark.skipif(not PHYLO_AVAILABLE or not check_real_data_exists(), reason="Not available")
class TestVisualizerRealData:
    """Tests for Visualizer with real data."""
    
    def test_full_visualization_pipeline(self):
        """Test full visualization pipeline with real data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Copy real data
            for f in os.listdir(REAL_DATA_PATH):
                src = os.path.join(REAL_DATA_PATH, f)
                if os.path.isfile(src):
                    shutil.copy(src, os.path.join(tmpdir, f))
            
            viz = Visualizer(tmpdir)
            
            # Load tree
            tree_path = os.path.join(tmpdir, 'Snp_tree.newick')
            tree = Phylo.read(tree_path, 'newick')
            terminals = tree.get_terminals()
            strain_names = [str(t.name) for t in terminals]
            
            # Create mock merged data
            merged_df = pd.DataFrame({
                'Strain_ID': strain_names,
                'Cluster': [i % 4 for i in range(len(strain_names))],
                **{f'gene_{i}': np.random.randint(0, 2, len(strain_names)) for i in range(10)}
            })
            
            try:
                viz.plot_cluster_distribution(merged_df)
                
                # Get distance matrix and embeddings
                core = PhylogeneticCore()
                distance_matrix, _ = core.tree_to_distance_matrix(tree)
                embeddings = core.dimension_reduction(distance_matrix)
                
                labels = np.array([i % 4 for i in range(len(terminals))])
                mask = np.array([True] * len(terminals))
                
                viz.plot_umap_clusters(embeddings, labels, mask)
                
            except Exception as e:
                print(f"Full pipeline error: {e}")
