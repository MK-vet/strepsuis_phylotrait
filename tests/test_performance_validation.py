#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Performance Validation and Benchmarking - strepsuis-phylotrait

This module validates phylogenetic analysis performance characteristics.
Results are saved to validation/PERFORMANCE_VALIDATION_REPORT.md
"""

import pytest
import numpy as np
import pandas as pd
import time
import json
from datetime import datetime
from pathlib import Path
from functools import wraps
import tracemalloc
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist


class PerformanceReport:
    """Collect and save performance validation results."""
    
    def __init__(self):
        self.benchmarks = []
        self.scalability_tests = []
        self.optimizations = []
    
    def add_benchmark(self, name, operation, data_size, time_ms, throughput, memory_mb=None):
        self.benchmarks.append({
            "name": name,
            "operation": operation,
            "data_size": data_size,
            "time_ms": round(time_ms, 2),
            "throughput": throughput,
            "memory_mb": round(memory_mb, 2) if memory_mb else None
        })
    
    def add_scalability(self, operation, sizes, times, memory_usage=None):
        self.scalability_tests.append({
            "operation": operation,
            "sizes": sizes,
            "times_ms": [round(t, 2) for t in times]
        })
    
    def save_report(self, output_dir):
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        report_path = output_path / "PERFORMANCE_VALIDATION_REPORT.md"
        
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("# Performance Validation Report - strepsuis-phylotrait\n\n")
            f.write(f"**Generated:** {datetime.now().isoformat()}\n\n")
            f.write("---\n\n")
            
            f.write("## Benchmark Summary\n\n")
            f.write("| Operation | Data Size | Time (ms) | Throughput |\n")
            f.write("|-----------|-----------|-----------|------------|\n")
            
            for b in self.benchmarks:
                f.write(f"| {b['name']} | {b['data_size']} | {b['time_ms']} | {b['throughput']} |\n")
            
            f.write("\n---\n\n")
            f.write("## Scalability Analysis\n\n")
            
            for s in self.scalability_tests:
                f.write(f"### {s['operation']}\n\n")
                f.write("| Data Size | Time (ms) |\n")
                f.write("|-----------|----------|\n")
                for i, size in enumerate(s['sizes']):
                    f.write(f"| {size} | {s['times_ms'][i]} |\n")
                f.write("\n")
        
        json_path = output_path / "performance_validation_results.json"
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump({
                "timestamp": datetime.now().isoformat(),
                "benchmarks": self.benchmarks,
                "scalability_tests": self.scalability_tests
            }, f, indent=2)
        
        return len(self.benchmarks)


report = PerformanceReport()


def generate_test_data(n_strains, n_features, prevalence=0.3):
    np.random.seed(42)
    data = np.random.binomial(1, prevalence, (n_strains, n_features))
    return data


class TestDistanceMatrixPerformance:
    """Benchmark distance matrix calculations."""
    
    def test_jaccard_distance_performance(self):
        """Benchmark Jaccard distance matrix calculation."""
        sizes = [50, 100, 200, 500]
        times = []
        n_features = 30
        
        for n in sizes:
            data = generate_test_data(n, n_features)
            
            start = time.perf_counter()
            distances = pdist(data, metric='jaccard')
            time_ms = (time.perf_counter() - start) * 1000
            times.append(time_ms)
            
            n_pairs = n * (n - 1) // 2
            report.add_benchmark(
                f"Jaccard Distance (n={n})",
                "jaccard_distance",
                f"{n_pairs} pairs",
                time_ms,
                f"{n_pairs/time_ms*1000:.0f} pairs/s",
                None
            )
        
        report.add_scalability("Jaccard Distance", sizes, times)
        assert len(times) == len(sizes)
    
    def test_hamming_distance_performance(self):
        """Benchmark Hamming distance matrix calculation."""
        sizes = [50, 100, 200, 500]
        times = []
        n_features = 30
        
        for n in sizes:
            data = generate_test_data(n, n_features)
            
            start = time.perf_counter()
            distances = pdist(data, metric='hamming')
            time_ms = (time.perf_counter() - start) * 1000
            times.append(time_ms)
            
            n_pairs = n * (n - 1) // 2
            report.add_benchmark(
                f"Hamming Distance (n={n})",
                "hamming_distance",
                f"{n_pairs} pairs",
                time_ms,
                f"{n_pairs/time_ms*1000:.0f} pairs/s",
                None
            )
        
        report.add_scalability("Hamming Distance", sizes, times)
        assert len(times) == len(sizes)


class TestHierarchicalClusteringPerformance:
    """Benchmark hierarchical clustering performance."""
    
    def test_hierarchical_clustering_performance(self):
        """Benchmark hierarchical clustering."""
        sizes = [50, 100, 200, 500]
        times = []
        n_features = 30
        
        for n in sizes:
            data = generate_test_data(n, n_features)
            distances = pdist(data, metric='jaccard')
            
            start = time.perf_counter()
            Z = linkage(distances, method='average')
            clusters = fcluster(Z, t=0.5, criterion='distance')
            time_ms = (time.perf_counter() - start) * 1000
            times.append(time_ms)
            
            report.add_benchmark(
                f"Hierarchical Clustering (n={n})",
                "hierarchical_clustering",
                f"{n} samples",
                time_ms,
                f"{n/time_ms*1000:.0f} samples/s",
                None
            )
        
        report.add_scalability("Hierarchical Clustering", sizes, times)
        assert len(times) == len(sizes)


class TestCorrelationMatrixPerformance:
    """Benchmark correlation matrix calculations."""
    
    def test_correlation_matrix_performance(self):
        """Benchmark correlation matrix calculation."""
        sizes = [20, 50, 100, 200]
        times = []
        n_strains = 100
        
        for n_features in sizes:
            data = generate_test_data(n_strains, n_features)
            df = pd.DataFrame(data)
            
            start = time.perf_counter()
            corr_matrix = df.corr()
            time_ms = (time.perf_counter() - start) * 1000
            times.append(time_ms)
            
            report.add_benchmark(
                f"Correlation Matrix (f={n_features})",
                "correlation_matrix",
                f"{n_features}x{n_features}",
                time_ms,
                f"{n_features*n_features/time_ms*1000:.0f} elem/s",
                None
            )
        
        report.add_scalability("Correlation Matrix", sizes, times)
        assert len(times) == len(sizes)


class TestDiversityMetricsPerformance:
    """Benchmark diversity metric calculations."""
    
    def test_richness_calculation_performance(self):
        """Benchmark richness calculation."""
        sizes = [100, 500, 1000, 2000]
        times = []
        n_features = 30
        
        for n in sizes:
            data = generate_test_data(n, n_features)
            
            start = time.perf_counter()
            richness = data.sum(axis=1)
            mean_richness = richness.mean()
            std_richness = richness.std()
            time_ms = (time.perf_counter() - start) * 1000
            times.append(time_ms)
            
            report.add_benchmark(
                f"Richness Calculation (n={n})",
                "richness",
                f"{n} samples",
                time_ms,
                f"{n/time_ms*1000:.0f} samples/s",
                None
            )
        
        report.add_scalability("Richness Calculation", sizes, times)
        assert len(times) == len(sizes)
    
    def test_beta_diversity_performance(self):
        """Benchmark beta diversity calculation."""
        sizes = [50, 100, 200]
        times = []
        n_features = 30
        
        for n in sizes:
            data = generate_test_data(n, n_features)
            
            start = time.perf_counter()
            jaccard_distances = []
            for i in range(n):
                for j in range(i+1, n):
                    intersection = np.sum(data[i] & data[j])
                    union = np.sum(data[i] | data[j])
                    if union > 0:
                        jaccard = 1 - intersection / union
                    else:
                        jaccard = 0
                    jaccard_distances.append(jaccard)
            mean_jaccard = np.mean(jaccard_distances)
            time_ms = (time.perf_counter() - start) * 1000
            times.append(time_ms)
            
            n_pairs = n * (n - 1) // 2
            report.add_benchmark(
                f"Beta Diversity (n={n})",
                "beta_diversity",
                f"{n_pairs} pairs",
                time_ms,
                f"{n_pairs/time_ms*1000:.0f} pairs/s",
                None
            )
        
        report.add_scalability("Beta Diversity", sizes, times)
        assert len(times) == len(sizes)


class TestPermutationTestPerformance:
    """Benchmark permutation test performance."""
    
    def test_permutation_test_performance(self):
        """Benchmark permutation test."""
        sizes = [100, 500, 1000]
        times = []
        n_strains = 100
        n_features = 10
        
        for n_permutations in sizes:
            data = generate_test_data(n_strains, n_features)
            trait = data[:, 0]
            
            observed_stat = np.corrcoef(trait, data[:, 1])[0, 1]
            
            start = time.perf_counter()
            null_stats = []
            for _ in range(n_permutations):
                permuted = np.random.permutation(trait)
                null_stat = np.corrcoef(permuted, data[:, 1])[0, 1]
                null_stats.append(null_stat)
            p_value = np.mean(np.abs(null_stats) >= np.abs(observed_stat))
            time_ms = (time.perf_counter() - start) * 1000
            times.append(time_ms)
            
            report.add_benchmark(
                f"Permutation Test (n={n_permutations})",
                "permutation_test",
                f"{n_permutations} permutations",
                time_ms,
                f"{n_permutations/time_ms*1000:.0f} perm/s",
                None
            )
        
        report.add_scalability("Permutation Test", sizes, times)
        assert len(times) == len(sizes)


@pytest.fixture(scope="session", autouse=True)
def save_performance_report():
    yield
    output_dir = Path(__file__).parent.parent / "validation"
    n_benchmarks = report.save_report(output_dir)
    print(f"\n{'='*60}")
    print(f"PERFORMANCE VALIDATION REPORT - strepsuis-phylotrait")
    print(f"Total Benchmarks: {n_benchmarks}")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
