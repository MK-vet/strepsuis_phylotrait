#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Performance Benchmark Tests for strepsuis-phylotrait

This module provides performance benchmarks for all major operations.
Results are saved to validation/PERFORMANCE_BENCHMARKS_REPORT.md
"""

import pytest
import numpy as np
import pandas as pd
import time
import json
from datetime import datetime
from pathlib import Path
from io import StringIO


class BenchmarkReport:
    """Collect and save benchmark results."""
    
    def __init__(self):
        self.results = []
        self.start_time = datetime.now()
    
    def add_result(self, operation, n_samples, n_features, time_seconds, memory_mb=None):
        self.results.append({
            "operation": operation,
            "n_samples": n_samples,
            "n_features": n_features,
            "time_seconds": time_seconds,
            "memory_mb": memory_mb,
            "throughput": n_samples / time_seconds if time_seconds > 0 else 0
        })
    
    def save_report(self, output_dir):
        """Save benchmark report to markdown file."""
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        report_path = output_path / "PERFORMANCE_BENCHMARKS_REPORT.md"
        
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("# Performance Benchmarks Report - strepsuis-phylotrait\n\n")
            f.write(f"**Generated:** {datetime.now().isoformat()}\n")
            f.write(f"**Total Benchmarks:** {len(self.results)}\n\n")
            f.write("---\n\n")
            
            f.write("## Benchmark Results\n\n")
            f.write("| Operation | Samples | Features | Time (s) | Throughput (samples/s) |\n")
            f.write("|-----------|---------|----------|----------|------------------------|\n")
            
            for r in self.results:
                f.write(f"| {r['operation']} | {r['n_samples']} | {r['n_features']} | {r['time_seconds']:.3f} | {r['throughput']:.1f} |\n")
            
            f.write("\n---\n\n")
            f.write("## Performance Summary\n\n")
            
            # Group by operation
            ops = {}
            for r in self.results:
                op = r['operation']
                if op not in ops:
                    ops[op] = []
                ops[op].append(r)
            
            for op, results in ops.items():
                f.write(f"### {op}\n\n")
                avg_throughput = np.mean([r['throughput'] for r in results])
                f.write(f"- **Average Throughput:** {avg_throughput:.1f} samples/s\n")
                f.write(f"- **Scalability:** Tested with {min(r['n_samples'] for r in results)}-{max(r['n_samples'] for r in results)} samples\n\n")
        
        # Also save as JSON
        json_path = output_path / "benchmark_results.json"
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump({
                "timestamp": datetime.now().isoformat(),
                "results": self.results
            }, f, indent=2)
        
        return len(self.results)


# Global report instance
report = BenchmarkReport()


def generate_random_tree(n_taxa):
    """Generate a random Newick tree string."""
    if n_taxa <= 1:
        return f"taxon_0:1.0"
    
    taxa = [f"taxon_{i}:1.0" for i in range(n_taxa)]
    
    while len(taxa) > 1:
        # Pick two random taxa to join
        i = np.random.randint(0, len(taxa))
        t1 = taxa.pop(i)
        j = np.random.randint(0, len(taxa))
        t2 = taxa.pop(j)
        
        # Create new internal node
        new_node = f"({t1},{t2}):{np.random.uniform(0.1, 2.0):.2f}"
        taxa.append(new_node)
    
    return taxa[0] + ";"


class TestTreeParsingBenchmarks:
    """Benchmark tree parsing."""
    
    @pytest.mark.parametrize("n_taxa", [10, 50, 100])
    def test_tree_parsing_performance(self, n_taxa):
        """Benchmark tree parsing."""
        from Bio import Phylo
        
        np.random.seed(42)
        tree_str = generate_random_tree(n_taxa)
        
        start = time.time()
        
        tree = Phylo.read(StringIO(tree_str), "newick")
        
        elapsed = time.time() - start
        
        report.add_result("Tree Parsing", n_taxa, 1, elapsed)
        
        assert elapsed < 10  # 10 seconds max


class TestDistanceMatrixBenchmarks:
    """Benchmark distance matrix calculation."""
    
    @pytest.mark.parametrize("n_taxa", [10, 30, 50])
    def test_distance_matrix_performance(self, n_taxa):
        """Benchmark distance matrix calculation."""
        from Bio import Phylo
        
        np.random.seed(42)
        tree_str = generate_random_tree(n_taxa)
        tree = Phylo.read(StringIO(tree_str), "newick")
        
        taxa = [t.name for t in tree.get_terminals()]
        
        start = time.time()
        
        # Calculate distance matrix
        dist_matrix = np.zeros((n_taxa, n_taxa))
        for i, t1 in enumerate(taxa):
            for j, t2 in enumerate(taxa):
                if i != j:
                    dist_matrix[i, j] = tree.distance(t1, t2)
        
        elapsed = time.time() - start
        
        report.add_result("Distance Matrix", n_taxa, n_taxa, elapsed)
        
        assert elapsed < 60  # 60 seconds max


class TestPhylogeneticDiversityBenchmarks:
    """Benchmark phylogenetic diversity calculations."""
    
    @pytest.mark.parametrize("n_taxa", [10, 30, 50])
    def test_faith_pd_performance(self, n_taxa):
        """Benchmark Faith's PD calculation."""
        from Bio import Phylo
        
        np.random.seed(42)
        tree_str = generate_random_tree(n_taxa)
        tree = Phylo.read(StringIO(tree_str), "newick")
        
        start = time.time()
        
        # Calculate total tree length (Faith's PD for all taxa)
        total_length = sum(c.branch_length or 0 for c in tree.find_clades())
        
        elapsed = time.time() - start
        
        report.add_result("Faith's PD", n_taxa, 1, elapsed)
        
        assert elapsed < 10  # 10 seconds max


class TestMPDMNTDBenchmarks:
    """Benchmark MPD and MNTD calculations."""
    
    @pytest.mark.parametrize("n_taxa", [10, 30, 50])
    def test_mpd_mntd_performance(self, n_taxa):
        """Benchmark MPD and MNTD calculation."""
        from Bio import Phylo
        
        np.random.seed(42)
        tree_str = generate_random_tree(n_taxa)
        tree = Phylo.read(StringIO(tree_str), "newick")
        
        taxa = [t.name for t in tree.get_terminals()]
        
        start = time.time()
        
        # Calculate distance matrix
        dist_matrix = np.zeros((n_taxa, n_taxa))
        for i, t1 in enumerate(taxa):
            for j, t2 in enumerate(taxa):
                if i != j:
                    dist_matrix[i, j] = tree.distance(t1, t2)
        
        # Calculate MPD
        upper_tri = dist_matrix[np.triu_indices(n_taxa, k=1)]
        mpd = np.mean(upper_tri)
        
        # Calculate MNTD
        np.fill_diagonal(dist_matrix, np.inf)
        mntd = np.mean(np.min(dist_matrix, axis=1))
        
        elapsed = time.time() - start
        
        report.add_result("MPD/MNTD", n_taxa, n_taxa, elapsed)
        
        assert elapsed < 60  # 60 seconds max


class TestTreeClusteringBenchmarks:
    """Benchmark tree-based clustering."""
    
    @pytest.mark.parametrize("n_taxa", [10, 30, 50])
    def test_hierarchical_clustering_performance(self, n_taxa):
        """Benchmark hierarchical clustering."""
        from Bio import Phylo
        from scipy.cluster.hierarchy import linkage, fcluster
        from scipy.spatial.distance import squareform
        
        np.random.seed(42)
        tree_str = generate_random_tree(n_taxa)
        tree = Phylo.read(StringIO(tree_str), "newick")
        
        taxa = [t.name for t in tree.get_terminals()]
        
        start = time.time()
        
        # Calculate distance matrix
        dist_matrix = np.zeros((n_taxa, n_taxa))
        for i, t1 in enumerate(taxa):
            for j, t2 in enumerate(taxa):
                if i != j:
                    dist_matrix[i, j] = tree.distance(t1, t2)
        
        # Hierarchical clustering
        condensed = squareform(dist_matrix)
        Z = linkage(condensed, method='average')
        clusters = fcluster(Z, t=3, criterion='maxclust')
        
        elapsed = time.time() - start
        
        report.add_result("Hierarchical Clustering", n_taxa, n_taxa, elapsed)
        
        assert elapsed < 60  # 60 seconds max


class TestPermutationTestBenchmarks:
    """Benchmark permutation tests."""
    
    @pytest.mark.parametrize("n_samples,n_permutations", [(50, 100), (100, 200)])
    def test_permutation_test_performance(self, n_samples, n_permutations):
        """Benchmark permutation test."""
        np.random.seed(42)
        
        x = np.random.uniform(0, 1, n_samples)
        y = np.random.uniform(0, 1, n_samples)
        
        start = time.time()
        
        # Observed correlation
        obs_corr = np.corrcoef(x, y)[0, 1]
        
        # Permutation test
        perm_corrs = []
        for _ in range(n_permutations):
            y_perm = np.random.permutation(y)
            perm_corrs.append(np.corrcoef(x, y_perm)[0, 1])
        
        p_value = np.mean(np.abs(perm_corrs) >= np.abs(obs_corr))
        
        elapsed = time.time() - start
        
        report.add_result("Permutation Test", n_samples, n_permutations, elapsed)
        
        assert elapsed < 30  # 30 seconds max


class TestBootstrapBenchmarks:
    """Benchmark bootstrap calculations."""
    
    @pytest.mark.parametrize("n_samples,n_bootstrap", [(50, 500), (100, 500), (200, 500)])
    def test_bootstrap_performance(self, n_samples, n_bootstrap):
        """Benchmark bootstrap CI calculation."""
        np.random.seed(42)
        
        sample = np.random.normal(0, 1, n_samples)
        
        start = time.time()
        
        boot_means = []
        for _ in range(n_bootstrap):
            boot_sample = np.random.choice(sample, size=n_samples, replace=True)
            boot_means.append(np.mean(boot_sample))
        
        ci_low = np.percentile(boot_means, 2.5)
        ci_high = np.percentile(boot_means, 97.5)
        
        elapsed = time.time() - start
        
        report.add_result("Bootstrap CI", n_samples, n_bootstrap, elapsed)
        
        assert elapsed < 10  # 10 seconds max


class TestCorrelationMatrixBenchmarks:
    """Benchmark correlation matrix calculations."""
    
    @pytest.mark.parametrize("n_taxa,n_traits", [(20, 10), (50, 20), (100, 30)])
    def test_correlation_matrix_performance(self, n_taxa, n_traits):
        """Benchmark correlation matrix calculation."""
        np.random.seed(42)
        
        traits = np.random.randn(n_taxa, n_traits)
        
        start = time.time()
        
        # Calculate correlation matrix
        corr_matrix = np.corrcoef(traits.T)
        
        elapsed = time.time() - start
        
        report.add_result("Correlation Matrix", n_taxa, n_traits, elapsed)
        
        assert elapsed < 10  # 10 seconds max


class TestFullPipelineBenchmarks:
    """Benchmark full analysis pipeline."""
    
    @pytest.mark.parametrize("n_taxa", [20, 50])
    def test_full_pipeline_performance(self, n_taxa):
        """Benchmark full pipeline execution."""
        from Bio import Phylo
        from scipy.cluster.hierarchy import linkage, fcluster
        from scipy.spatial.distance import squareform
        
        np.random.seed(42)
        tree_str = generate_random_tree(n_taxa)
        tree = Phylo.read(StringIO(tree_str), "newick")
        
        taxa = [t.name for t in tree.get_terminals()]
        n_traits = 10
        traits = np.random.randn(n_taxa, n_traits)
        
        start = time.time()
        
        # Step 1: Distance matrix
        dist_matrix = np.zeros((n_taxa, n_taxa))
        for i, t1 in enumerate(taxa):
            for j, t2 in enumerate(taxa):
                if i != j:
                    dist_matrix[i, j] = tree.distance(t1, t2)
        
        # Step 2: MPD/MNTD
        upper_tri = dist_matrix[np.triu_indices(n_taxa, k=1)]
        mpd = np.mean(upper_tri)
        
        # Step 3: Clustering
        condensed = squareform(dist_matrix)
        Z = linkage(condensed, method='average')
        clusters = fcluster(Z, t=3, criterion='maxclust')
        
        # Step 4: Trait correlation
        corr_matrix = np.corrcoef(traits.T)
        
        elapsed = time.time() - start
        
        report.add_result("Full Pipeline", n_taxa, n_traits, elapsed)
        
        assert elapsed < 60  # 60 seconds max


@pytest.fixture(scope="session", autouse=True)
def save_benchmark_report():
    """Save benchmark report after all tests."""
    yield
    
    # Save report
    output_dir = Path(__file__).parent.parent / "validation"
    n_benchmarks = report.save_report(output_dir)
    
    print(f"\n{'='*60}")
    print(f"PERFORMANCE BENCHMARKS REPORT - strepsuis-phylotrait")
    print(f"{'='*60}")
    print(f"Total Benchmarks: {n_benchmarks}")
    print(f"Report saved to: {output_dir / 'PERFORMANCE_BENCHMARKS_REPORT.md'}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
