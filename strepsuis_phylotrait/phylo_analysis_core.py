#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Phylogenetic Analysis Core Module
=================================

Phylogenetic clustering and binary trait analysis for bacterial genomics.

Features:
    - Tree-aware clustering respecting phylogenetic structure
    - Phylogenetic diversity metrics (Faith's PD, MPD, MNTD)
    - Phylogenetic signal detection (D-statistic, permutation tests)
    - Binary trait analysis (AMR, virulence factors)
    - Multiple testing correction using FDR
    - Interactive HTML report with DataTables and Plotly
    - Excel report generation

Author: MK-vet
Version: 1.0.0
License: MIT
"""

import base64
import logging
import os
from pathlib import Path
import random
import sys
import time
import warnings
from datetime import datetime
from functools import partial
from io import BytesIO
from itertools import combinations
from multiprocessing import Pool, cpu_count

import jinja2
import matplotlib.pyplot as plt
import numpy as np
import optuna
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import prince
import psutil
import seaborn as sns
try:
    import weasyprint
    WEASYPRINT_AVAILABLE = True
except (ImportError, OSError):
    weasyprint = None
    WEASYPRINT_AVAILABLE = False
from Bio import Phylo

# Excel report generation
try:
    from .excel_report_utils import ExcelReportGenerator, sanitize_sheet_name
except ImportError:
    from strepsuis_phylotrait.excel_report_utils import ExcelReportGenerator, sanitize_sheet_name

# For association rule mining:
from mlxtend.frequent_patterns import apriori, association_rules
from scipy.spatial.distance import cdist
from scipy.stats import chi2_contingency, fisher_exact
from sklearn.cluster import DBSCAN, KMeans
from sklearn.ensemble import IsolationForest, RandomForestClassifier
from sklearn.metrics import silhouette_score
from sklearn.mixture import GaussianMixture
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
from umap import UMAP

# Set global random seeds and control threading for reproducibility
random.seed(42)
np.random.seed(42)
os.environ["PYTHONHASHSEED"] = "0"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

warnings.filterwarnings("ignore", category=UserWarning)
sns.set(style="whitegrid")

LARGE_DATASET_THRESHOLD = 80
LARGE_FEATURE_THRESHOLD = 50
MAX_ITEMSETS_LIMIT = 2000
MAX_BOOTSTRAP_LARGE_DATASET = 30
MIN_ESTIMATORS = 30
BASE_ESTIMATORS = 100
ESTIMATOR_REDUCTION_FACTOR = 10
MAX_LEN_LARGE_DATASET = 2
MAX_LEN_NORMAL = 3
MAX_ACTIVE_TRAITS_LARGE_DATASET = 20
MAX_BOOTSTRAP_LOG_ODDS_LARGE_DATASET = 10


###############################################################################
# LOGGING AND PROGRESS TRACKING UTILITIES
###############################################################################
def setup_logging(output_folder):
    """
    Set up comprehensive logging to file and console.

    Creates log file and configures both file and console handlers with
    timestamp formatting.

    Parameters
    ----------
    output_folder : str
        Directory where log file will be saved

    Returns
    -------
    str
        Path to created log file

    Notes
    -----
    Removes existing handlers to avoid duplicate logging.
    Log file is named 'phylogenetic_analysis.log'.
    """
    os.makedirs(output_folder, exist_ok=True)
    log_file = os.path.join(output_folder, "phylogenetic_analysis.log")

    # Create formatter
    formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )

    # File handler
    file_handler = logging.FileHandler(log_file, mode="w", delay=True)
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)

    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)

    # Configure root logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Remove existing handlers and close them
    for handler in list(logger.handlers):
        try:
            handler.close()
        except Exception:
            pass
    logger.handlers = []

    # Add handlers
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return log_file


def _safe_log(message: str):
    """Log a message, falling back to print if handlers fail."""
    try:
        logging.info(message)
    except (FileNotFoundError, PermissionError, OSError):
        try:
            for handler in list(logging.getLogger().handlers):
                try:
                    handler.close()
                except Exception:
                    pass
            logging.getLogger().handlers = []
        except Exception:
            pass
        print(message)


def print_memory_usage():
    """Log current memory usage of the process."""
    process = psutil.Process(os.getpid())
    mem_mb = process.memory_info().rss / (1024**2)
    _safe_log(f"Memory Usage: {mem_mb:.2f} MB")
    return mem_mb


def print_section_header(title):
    """Print a formatted section header."""
    separator = "=" * 80
    _safe_log("")
    _safe_log(separator)
    _safe_log(f"  {title}")
    _safe_log(separator)
    _safe_log("")


def print_step(step_num, total_steps, description):
    """Print a formatted step indicator."""
    _safe_log(f"\n>>> Step {step_num}/{total_steps}: {description}")
    _safe_log(f"    Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print_memory_usage()


###############################################################################
# 1. CREATE TEMPLATE DIRECTORY & HTML TEMPLATE
###############################################################################
def create_template_directory(base_dir="."):
    """Creates a 'templates' folder if needed and writes a comprehensive HTML template."""
    try:
        resolved_base = os.path.abspath(base_dir)
    except FileNotFoundError:
        resolved_base = None
    if not resolved_base or not os.path.isdir(resolved_base):
        resolved_base = os.path.dirname(os.path.abspath(__file__))
    templates_dir = os.path.join(resolved_base, "templates")
    os.makedirs(templates_dir, exist_ok=True)
    template_content = r"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>{{ title }}</title>
    <!-- Bootstrap CSS -->
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <!-- Plotly JS -->
    <script src="https://cdn.plot.ly/plotly-2.12.1.min.js"></script>
    <!-- jQuery & DataTables -->
    <script src="https://code.jquery.com/jquery-3.6.3.min.js"></script>
    <link rel="stylesheet" href="https://cdn.datatables.net/1.13.4/css/jquery.dataTables.min.css">
    <!-- DataTables Buttons (for export options) -->
    <link rel="stylesheet" href="https://cdn.datatables.net/buttons/2.3.6/css/buttons.dataTables.min.css">
    <script src="https://cdn.datatables.net/1.13.4/js/jquery.dataTables.min.js"></script>
    <script src="https://cdn.datatables.net/buttons/2.3.6/js/dataTables.buttons.min.js"></script>
    <script src="https://cdn.datatables.net/buttons/2.3.6/js/buttons.html5.min.js"></script>
    <style>
        body { padding: 20px; font-family: Arial, sans-serif; }
        .section { margin-bottom: 40px; border-bottom: 1px solid #eee; padding-bottom: 20px; }
        .plot-container { margin-bottom: 20px; }
        .table-container { max-height: 600px; overflow-y: auto; margin-bottom: 20px; }
        table.dataTable { width: 100% !important; }
        h1, h2, h3 { color: #2c3e50; }
        h2 { margin-top: 30px; }
        h3 { margin-top: 20px; }
        .nav-tabs { margin-bottom: 20px; }
        .tab-pane { padding: 15px; }
        .export-btn { margin-left: 5px; }
    </style>
</head>
<body>
<div class="container-fluid">
    <div class="row mb-4">
        <div class="col">
            <h1>{{ title }}</h1>
            <p class="text-muted">Generated on: {{ date }}</p>
        </div>
    </div>
    <!-- 1. Overview -->
    <div class="section">
        <h2>1. Overview</h2>
        <div class="row">
            <div class="col-md-6">
                <h3>Configuration</h3>
                <div class="table-container">
                    <table class="table table-sm table-striped">
                        <tbody>
                        {% for key, value in config.items() %}
                        <tr>
                            <td><strong>{{ key }}</strong></td>
                            <td>{{ value }}</td>
                        </tr>
                        {% endfor %}
                        </tbody>
                    </table>
                </div>
            </div>
            <div class="col-md-6">
                <h3>Summary Statistics</h3>
                <div class="table-container">
                    {{ results.summary_stats|safe }}
                </div>
            </div>
        </div>
    </div>
    <!-- 2. Phylogenetic Tree -->
    <div class="section">
        <h2>2. Phylogenetic Tree</h2>
        <div class="row">
            <div class="col-md-8">
                <h3>Interactive Tree</h3>
                <div class="plot-container" id="plot-tree">
                    {{ results.tree_plot|safe }}
                </div>
            </div>
            <div class="col-md-4">
                <h3>Tree Statistics</h3>
                <div class="table-container">
                    {{ results.tree_stats|safe }}
                </div>
            </div>
        </div>
    </div>
    <!-- 3. Cluster Analysis -->
    <div class="section">
        <h2>3. Cluster Analysis</h2>
        <div class="row">
            <div class="col-md-6">
                <h3>3.1 UMAP Visualization</h3>
                <div class="plot-container" id="plot-umap">
                    {{ results.umap_plot|safe }}
                </div>
            </div>
            <div class="col-md-6">
                <h3>3.2 Cluster Distribution</h3>
                <div class="plot-container" id="plot-cluster-distribution">
                    {{ results.cluster_distribution|safe }}
                </div>
            </div>
        </div>
        <div class="row mt-4">
            <div class="col-md-12">
                <h3>3.3 Cluster Validation</h3>
                <div class="table-container">
                    {{ results.cluster_validation|safe }}
                </div>
            </div>
        </div>
    </div>
    <!-- 4. Evolutionary Analysis -->
    <div class="section">
        <h2>4. Evolutionary Analysis</h2>
        <div class="row">
            <div class="col-md-6">
                <h3>4.1 Evolutionary Metrics</h3>
                <div class="table-container">{{ results.evolutionary_metrics|safe }}</div>
            </div>
            <div class="col-md-6">
                <h3>4.2 Beta Diversity</h3>
                <div class="plot-container" id="plot-beta-diversity">
                    {{ results.beta_diversity|safe }}
                </div>
            </div>
        </div>
        <div class="row mt-4">
            <div class="col-md-6">
                <h3>4.3 Evolution Rates</h3>
                <div class="plot-container" id="plot-evolution-rates">
                    {{ results.evolution_rates|safe }}
                </div>
            </div>
            <div class="col-md-6">
                <h3>4.4 Phylogenetic Signal</h3>
                <div class="table-container">{{ results.phylogenetic_signal|safe }}</div>
            </div>
        </div>
    </div>
    <!-- 5. Trait Analysis -->
    <div class="section">
        <h2>5. Trait Analysis</h2>
        <div class="table-container">
            {{ results.trait_analysis|safe }}
        </div>
    </div>
    <!-- 6. Multiple Correspondence Analysis (MCA) -->
    <div class="section">
        <h2>6. Multiple Correspondence Analysis (MCA)</h2>
        <div class="row">
            <div class="col-md-6">
                <h3>6.1 MCA Row Points (Strains)</h3>
                <div class="plot-container" id="plot-mca-rows">
                    {{ results.mca_row_plot|safe }}
                </div>
            </div>
            <div class="col-md-6">
                <h3>6.2 MCA Column Points (Features)</h3>
                <div class="plot-container" id="plot-mca-cols">
                    {{ results.mca_column_plot|safe }}
                </div>
            </div>
        </div>
        <div class="row mt-4">
            <div class="col-md-12">
                <h3>6.3 MCA Summary</h3>
                <div class="table-container">{{ results.mca_summary|safe }}</div>
            </div>
        </div>
    </div>
    <!-- 7. Log-Odds Ratio Analysis -->
    <div class="section">
        <h2>7. Log-Odds Ratio Analysis</h2>
        <div class="row">
            <div class="col-md-6">
                <h3>7.1 Global Log-Odds</h3>
                <div class="table-container">{{ results.global_log_odds|safe }}</div>
            </div>
            <div class="col-md-6">
                <h3>7.2 Cluster Log-Odds</h3>
                <div class="table-container">{{ results.cluster_log_odds|safe }}</div>
            </div>
        </div>
    </div>
    <!-- 8. Association Rule Mining -->
    <div class="section">
        <h2>8. Association Rule Mining</h2>
        <div class="table-container">{{ results.association_rules|safe }}</div>
    </div>
    <!-- 9. Shared and Unique Features -->
    <div class="section">
        <h2>9. Shared and Unique Features</h2>
        <div class="row">
            <div class="col-md-6">
                <h3>9.1 Shared Features</h3>
                <div class="table-container">{{ results.shared_features|safe }}</div>
            </div>
            <div class="col-md-6">
                <h3>9.2 Unique Features</h3>
                <div class="table-container">{{ results.unique_features|safe }}</div>
            </div>
        </div>
    </div>
    <!-- 10. Pairwise FDR Post-Hoc Analysis -->
    <div class="section">
        <h2>10. Pairwise FDR Post-Hoc Analysis</h2>
        <div class="table-container">{{ results.pairwise_fdr|safe }}</div>
    </div>
    <!-- 11. Bootstrapping Analysis -->
    <div class="section">
        <h2>11. Bootstrapping Analysis</h2>
        <div class="row">
            <div class="col-md-6">
                <h3>11.1 Feature Importance with Bootstrap CI</h3>
                <div class="table-container">{{ results.bootstrap_feature_importance|safe }}</div>
            </div>
            <div class="col-md-6">
                <h3>11.2 Log-Odds with Bootstrap CI</h3>
                <div class="table-container">{{ results.bootstrap_log_odds|safe }}</div>
            </div>
        </div>
    </div>
    <!-- 12. All Raw Data Files -->
    <div class="section">
        <h2>12. All Raw Data Files</h2>
        <div class="table-container">{{ results.download_section|safe }}</div>
    </div>
    <!-- 13. Additional Analysis Results -->
    <div class="section">
        <h2>13. Additional Analysis Results</h2>
        <div class="table-container">{{ results.additional_results|safe }}</div>
    </div>
</div>

<!-- Bootstrap JS -->
<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js"></script>
<script>
    // Initialize DataTables with extended options for all tables with class "dataTable"
    document.addEventListener('DOMContentLoaded', function() {
        $('.dataTable').DataTable({
            paging: true,
            searching: true,
            ordering: true,
            info: true,
            lengthMenu: [10, 25, 50, 100, -1],
            pageLength: 25,
            language: {
                lengthMenu: "Show _MENU_ entries",
                search: "Search:"
            },
            dom: 'Blfrtip',
            buttons: ['copy', 'csv', 'excel']
        });
    });
    // Export table to CSV
    function exportTableToCSV(tableId, filename) {
        var csv = [];
        var rows = document.getElementById(tableId).querySelectorAll('tr');
        for (var i = 0; i < rows.length; i++) {
            var row = [], cols = rows[i].querySelectorAll('td, th');
            for (var j = 0; j < cols.length; j++) {
                var data = cols[j].innerText.replace(/(\r\n|\n|\r)/gm, ' ').replace(/,/g, ';');
                row.push('"' + data + '"');
            }
            csv.push(row.join(','));
        }
        var csv_string = csv.join('\n');
        var a = document.createElement('a');
        a.href = 'data:text/csv;charset=utf-8,' + encodeURIComponent(csv_string);
        a.target = '_blank';
        a.download = filename;
        a.click();
    }
    // Export Plotly chart to PNG
    function exportPlotToPNG(plotId, filename) {
        var gd = document.getElementById(plotId);
        if(gd && gd.data) {
            Plotly.downloadImage(gd, {format: 'png', width: 1200, height: 800, filename: filename});
        } else {
            alert("Plot not found or invalid for exporting.");
        }
    }
</script>
</body>
</html>
"""
    with open(os.path.join(templates_dir, "report_template.html"), "w") as f:
        f.write(template_content)
    print("Template directory and 'report_template.html' created successfully.")


###############################################################################
# 2. HELPER FUNCTION TO FIX PICKLING ISSUE WITH BOOTSTRAP
###############################################################################
def _bootstrap_importance_func(args):
    """
    Top-level function to avoid 'local object' pickling error.
    Receives (seed, X, y).
    """
    seed, X, y = args
    np.random.seed(seed)
    n_samples = len(X)
    indices = np.random.choice(n_samples, size=n_samples, replace=True)
    X_boot, y_boot = X[indices], y[indices]
    rf = RandomForestClassifier(n_estimators=100, random_state=seed)
    rf.fit(X_boot, y_boot)
    return rf.feature_importances_


###############################################################################
# 3. PARALLEL PROCESSOR (OPTIONAL FOR TREE DISTANCES OR OTHER)
###############################################################################
class ParallelProcessor:
    """Handles parallel processing of time-consuming tasks."""

    def __init__(self, n_jobs=None):
        self.n_jobs = n_jobs if n_jobs is not None else max(1, cpu_count() - 1)

    def parallel_map(self, func, iterable):
        """
        Apply function to each item in iterable (sequential implementation).

        Parameters
        ----------
        func : callable
            Function to apply to each item
        iterable : iterable
            Items to process

        Returns
        -------
        list
            Results from applying func to each item
        """
        return [func(item) for item in iterable]

    def map(self, func, iterable):
        return self.parallel_map(func, iterable)

    def parallel_bootstrap(self, data, statistic_func=None, n_bootstrap=1000, **kwargs):
        if statistic_func is None:
            return np.array([])
        return ParallelProcessor.parallel_bootstrap(
            data, statistic_func, n_bootstrap=n_bootstrap, n_jobs=self.n_jobs, **kwargs
        )

    @staticmethod
    def parallel_tree_distance_matrix(tree, terminals, n_jobs=None):
        """
        Compute phylogenetic distance matrix using parallel processing.

        Calculates pairwise distances between all terminal nodes in the tree
        using multiprocessing for improved performance on large trees.

        Parameters
        ----------
        tree : Bio.Phylo.BaseTree.Tree
            Phylogenetic tree object
        terminals : list
            List of terminal nodes (leaves) from the tree
        n_jobs : int, optional
            Number of parallel jobs (default: CPU count - 1)

        Returns
        -------
        np.ndarray
            Symmetric distance matrix (n_terminals, n_terminals)

        Notes
        -----
        Uses multiprocessing.Pool to parallelize distance calculations.
        Significantly faster than sequential computation for trees with >100 taxa.
        """
        if n_jobs is None:
            n_jobs = max(1, cpu_count() - 1)
        n = len(terminals)
        distance_matrix = np.zeros((n, n))
        print(f"Creating a distance matrix for {n} strains using {n_jobs} processes...")

        def process_row(i, terminals, tree):
            row = np.zeros(len(terminals))
            strain1 = terminals[i]
            for j in range(len(terminals)):
                if i == j:
                    continue
                strain2 = terminals[j]
                try:
                    dist_val = tree.distance(str(strain1), str(strain2))
                    row[j] = dist_val
                except Exception as e:
                    print(f"Error computing distance between {strain1} and {strain2}: {e}")
                    row[j] = np.nan
            return i, row

        with Pool(processes=n_jobs) as pool:
            func = partial(process_row, terminals=terminals, tree=tree)
            results = pool.map(func, range(n))

        for i, row in results:
            distance_matrix[i] = row
            for j in range(n):
                if i != j:
                    distance_matrix[j, i] = distance_matrix[i, j]
        return distance_matrix

    @staticmethod
    def parallel_bootstrap(data, statistic_func, n_bootstrap=1000, n_jobs=None, **kwargs):
        """
        Perform parallel bootstrap resampling with custom statistic function.

        Uses multiprocessing to compute bootstrap confidence intervals efficiently.

        Parameters
        ----------
        data : array-like or DataFrame
            Input data to resample
        statistic_func : callable
            Function to compute statistic on each bootstrap sample
        n_bootstrap : int, default=1000
            Number of bootstrap iterations
        n_jobs : int, optional
            Number of parallel jobs (default: CPU count - 1)
        **kwargs : dict
            Additional arguments passed to statistic_func

        Returns
        -------
        np.ndarray
            Array of bootstrap statistics (n_bootstrap,)

        Examples
        --------
        >>> data = np.random.randn(100)
        >>> results = parallel_bootstrap(data, np.mean, n_bootstrap=1000)
        >>> ci_lower, ci_upper = np.percentile(results, [2.5, 97.5])
        """
        if n_jobs is None:
            n_jobs = max(1, cpu_count() - 1)
        n_samples = len(data)

        def bootstrap_sample(seed):
            np.random.seed(seed)
            indices = np.random.choice(n_samples, size=n_samples, replace=True)
            sample_data = data.iloc[indices] if hasattr(data, "iloc") else data[indices]
            return statistic_func(sample_data, **kwargs)

        with Pool(processes=n_jobs) as pool:
            bootstrap_results = pool.map(bootstrap_sample, range(n_bootstrap))
        return np.array(bootstrap_results)

    @staticmethod
    def parallel_feature_importance(X, y, n_bootstrap=1000, n_jobs=None):
        """
        Compute bootstrap feature importance using parallel processing.

        Uses Random Forest classifier to assess feature importance with
        bootstrap confidence intervals.

        Parameters
        ----------
        X : np.ndarray
            Feature matrix (n_samples, n_features)
        y : np.ndarray
            Target labels (n_samples,)
        n_bootstrap : int, default=1000
            Number of bootstrap iterations
        n_jobs : int, optional
            Number of parallel jobs (default: CPU count - 1)

        Returns
        -------
        np.ndarray
            Bootstrap feature importance scores (n_bootstrap, n_features)

        Notes
        -----
        Each bootstrap iteration trains a Random Forest on resampled data.
        Feature importance is measured by Gini importance.
        """
        if n_jobs is None:
            n_jobs = max(1, cpu_count() - 1)
        args_list = [(seed, X, y) for seed in range(n_bootstrap)]
        with Pool(processes=n_jobs) as pool:
            importance_scores = pool.map(_bootstrap_importance_func, args_list)
        return np.array(importance_scores)


###############################################################################
# 4. PHYLOGENETIC CORE
###############################################################################
class PhylogeneticCore:
    """Basic operations on phylogenetic trees."""

    @staticmethod
    def load_tree(tree_path):
        """
        Load phylogenetic tree from Newick format file.

        Parameters
        ----------
        tree_path : str
            Path to Newick format tree file

        Returns
        -------
        Bio.Phylo.BaseTree.Tree
            Parsed phylogenetic tree object

        Examples
        --------
        >>> tree = PhylogeneticCore.load_tree("tree.newick")
        >>> terminals = tree.get_terminals()
        """
        return Phylo.read(tree_path, "newick")

    @staticmethod
    def tree_to_distance_matrix(tree, parallel=False, n_jobs=None):
        """
        Convert phylogenetic tree to pairwise distance matrix.

        Computes patristic distances (sum of branch lengths) between all pairs
        of terminal nodes in the tree.

        Parameters
        ----------
        tree : Bio.Phylo.BaseTree.Tree
            Phylogenetic tree object
        parallel : bool, default=False
            Whether to use parallel processing
        n_jobs : int, optional
            Number of parallel jobs (only used if parallel=True)

        Returns
        -------
        distance_matrix : np.ndarray
            Symmetric distance matrix (n_terminals, n_terminals)
        terminals : list
            List of terminal nodes corresponding to matrix rows/columns

        Notes
        -----
        Parallel processing recommended for trees with >100 taxa.
        Distance values are patristic distances (sum of branch lengths).
        """
        terminals = tree.get_terminals()
        if parallel:
            distance_matrix = ParallelProcessor.parallel_tree_distance_matrix(
                tree, terminals, n_jobs=n_jobs
            )
        else:
            n = len(terminals)
            distance_matrix = np.zeros((n, n))
            for i, strain1 in enumerate(terminals):
                for j, strain2 in enumerate(terminals):
                    if i < j:
                        try:
                            dist_val = tree.distance(str(strain1), str(strain2))
                            distance_matrix[i, j] = dist_val
                            distance_matrix[j, i] = dist_val
                        except Exception as e:
                            print(f"Error computing distance: {e}")
        return distance_matrix, terminals

    @staticmethod
    def dimension_reduction(matrix, n_components=2, n_neighbors=15, min_dist=0.1, random_state=42):
        """
        Reduce dimensionality of distance matrix using UMAP.

        Projects high-dimensional phylogenetic distances to lower-dimensional
        space while preserving local and global structure.

        Parameters
        ----------
        matrix : np.ndarray
            Precomputed distance matrix (n_samples, n_samples)
        n_components : int, default=2
            Number of dimensions in embedding space
        n_neighbors : int, default=15
            Number of neighbors for UMAP graph construction
        min_dist : float, default=0.1
            Minimum distance between points in embedding
        random_state : int, default=42
            Random seed for reproducibility

        Returns
        -------
        np.ndarray
            Low-dimensional embedding (n_samples, n_components)

        Notes
        -----
        Uses UMAP algorithm optimized for precomputed distance matrices.
        Larger n_neighbors preserves more global structure.
        Smaller min_dist creates tighter clusters.
        """
        reducer = UMAP(
            n_components=n_components,
            metric="precomputed",
            random_state=random_state,
            n_neighbors=n_neighbors,
            min_dist=min_dist,
        )
        return reducer.fit_transform(matrix)

    @staticmethod
    def detect_outliers(embeddings, contamination=0.05, n_estimators=200, random_state=42):
        """
        Detect outlier samples using Isolation Forest.

        Identifies anomalous samples that deviate from typical clustering patterns.

        Parameters
        ----------
        embeddings : np.ndarray
            Low-dimensional embedding matrix (n_samples, n_components)
        contamination : float, default=0.05
            Expected proportion of outliers (0.0 to 0.5)
        n_estimators : int, default=200
            Number of isolation trees
        random_state : int, default=42
            Random seed for reproducibility

        Returns
        -------
        clean_embeddings : np.ndarray
            Embeddings with outliers removed
        mask : np.ndarray
            Boolean mask indicating inliers (True) and outliers (False)

        Notes
        -----
        Uses Isolation Forest algorithm for unsupervised outlier detection.
        Contamination parameter should reflect expected outlier proportion.
        """
        iso = IsolationForest(
            contamination=contamination, random_state=random_state, n_estimators=n_estimators
        )
        mask = iso.fit_predict(embeddings) != -1
        print(f"Identified {sum(~mask)} outliers out of {len(mask)} samples")
        return embeddings[mask], mask


###############################################################################
# 5.1 TREE-AWARE CLUSTERING MODULE
###############################################################################
class TreeAwareClusteringModule:
    """
    Perform tree-aware clustering respecting phylogenetic structure.

    Implements TreeCluster algorithm with monophyletic constraints and
    Phydelity-inspired methods for phylogenetically coherent clustering.

    Parameters
    ----------
    tree : Bio.Phylo.BaseTree.Tree
        Phylogenetic tree object
    terminals : list
        List of terminal nodes (leaves)
    n_clusters_range : tuple, default=(2, 20)
        Range of cluster numbers to consider
    seed : int, default=42
        Random seed for reproducibility
    **kwargs : dict
        Additional keyword arguments

    Attributes
    ----------
    tree : Bio.Phylo.BaseTree.Tree
        Input phylogenetic tree
    terminals : list
        Terminal nodes
    n_clusters_range : tuple
        Cluster range
    seed : int
        Random seed

    Notes
    -----
    Tree-aware clustering ensures clusters are monophyletic (single common ancestor).
    Uses patristic distance thresholds for cluster boundary detection.
    """

    def __init__(self, tree, terminals, n_clusters_range=(2, 20), seed=42, **kwargs):
        """
        Initialize tree-aware clustering module.

        Parameters
        ----------
        tree : Bio.Phylo.BaseTree.Tree
            Phylogenetic tree
        terminals : list
            Terminal nodes
        n_clusters_range : tuple, default=(2, 20)
            Range of cluster numbers
        seed : int, default=42
            Random seed
        **kwargs : dict
            Additional arguments
        """
        self.tree = tree
        self.terminals = terminals
        self.n_clusters_range = n_clusters_range
        self.seed = seed
        np.random.seed(seed)

    def tree_cluster_algorithm(self, distance_matrix, method="max", threshold=None):
        print(f"Performing TreeCluster with method: {method}")
        if threshold is None:
            if method == "max":
                threshold = self._auto_threshold_max(distance_matrix=distance_matrix)
            elif method == "sum":
                threshold = self._auto_threshold_sum(distance_matrix=distance_matrix)
            else:
                threshold = self._auto_threshold_avg(distance_matrix=distance_matrix)

        n_terminals = len(self.terminals)
        clusters = [[i] for i in range(n_terminals)]
        cluster_map = list(range(n_terminals))

        internal_nodes = [
            clade for clade in self.tree.find_clades(order="postorder") if not clade.is_terminal()
        ]

        for node in internal_nodes:
            terminals_in_subtree = []
            for leaf in node.get_terminals():
                for i, term in enumerate(self.terminals):
                    if str(term) == str(leaf):
                        terminals_in_subtree.append(i)
                        break

            if method == "max":
                is_valid = self._check_max_constraint(
                    terminals_in_subtree, distance_matrix, threshold
                )
            elif method == "sum":
                is_valid = self._check_sum_constraint(
                    terminals_in_subtree, distance_matrix, threshold
                )
            else:
                is_valid = self._check_avg_constraint(
                    terminals_in_subtree, distance_matrix, threshold
                )

            if is_valid:
                self._merge_clusters(terminals_in_subtree, clusters, cluster_map)

        labels = np.array([cluster_map[i] for i in range(n_terminals)])
        unique_labels = np.unique(labels)
        relabel_map = {old: new for new, old in enumerate(unique_labels)}
        labels = np.array([relabel_map[label] for label in labels])
        return labels

    def _auto_threshold_max(self, conservative_factor=7.0, distance_matrix=None):
        dists = []
        if distance_matrix is None:
            for clade in self.tree.get_nonterminals():
                leaves = clade.get_terminals()
                if len(leaves) >= 2:
                    max_dist = 0
                    for i, leaf1 in enumerate(leaves):
                        for leaf2 in leaves[i + 1 :]:
                            try:
                                dist = self.tree.distance(leaf1, leaf2)
                                max_dist = max(max_dist, dist)
                            except Exception:
                                pass
                    dists.append(max_dist)
        else:
            terminal_index = {str(term): idx for idx, term in enumerate(self.terminals)}
            clade_indices = {}
            for clade in self.tree.get_nonterminals():
                leaves = clade.get_terminals()
                if len(leaves) >= 2:
                    indices = clade_indices.get(id(clade))
                    if indices is None:
                        indices = [
                            terminal_index[str(leaf)]
                            for leaf in leaves
                            if str(leaf) in terminal_index
                        ]
                        clade_indices[id(clade)] = indices
                    if len(indices) >= 2:
                        max_dist = 0
                        has_dist = False
                        for i, t1 in enumerate(indices):
                            for t2 in indices[i + 1 :]:
                                dist = distance_matrix[t1, t2]
                                if np.isnan(dist):
                                    continue
                                has_dist = True
                                max_dist = max(max_dist, dist)
                        if has_dist:
                            dists.append(max_dist)
        if dists:
            q75, q25 = np.percentile(dists, [75, 25])
            iqr = q75 - q25
            threshold = np.median(dists) + conservative_factor * iqr
            return threshold
        return 0.1

    def _auto_threshold_sum(self, conservative_factor=7.0, distance_matrix=None):
        dists = []
        if distance_matrix is None:
            for clade in self.tree.get_nonterminals():
                leaves = clade.get_terminals()
                if len(leaves) >= 2:
                    sum_dist = 0
                    for i, leaf1 in enumerate(leaves):
                        for leaf2 in leaves[i + 1 :]:
                            try:
                                dist = self.tree.distance(leaf1, leaf2)
                                sum_dist += dist
                            except Exception:
                                pass
                    dists.append(sum_dist)
        else:
            terminal_index = {str(term): idx for idx, term in enumerate(self.terminals)}
            clade_indices = {}
            for clade in self.tree.get_nonterminals():
                leaves = clade.get_terminals()
                if len(leaves) >= 2:
                    indices = clade_indices.get(id(clade))
                    if indices is None:
                        indices = [
                            terminal_index[str(leaf)]
                            for leaf in leaves
                            if str(leaf) in terminal_index
                        ]
                        clade_indices[id(clade)] = indices
                    if len(indices) >= 2:
                        sum_dist = 0
                        has_dist = False
                        for i, t1 in enumerate(indices):
                            for t2 in indices[i + 1 :]:
                                dist = distance_matrix[t1, t2]
                                if np.isnan(dist):
                                    continue
                                has_dist = True
                                sum_dist += dist
                        if has_dist:
                            dists.append(sum_dist)
        if dists:
            q75, q25 = np.percentile(dists, [75, 25])
            iqr = q75 - q25
            threshold = np.median(dists) + conservative_factor * iqr
            return threshold
        return 1.0

    def _auto_threshold_avg(self, conservative_factor=7.0, distance_matrix=None):
        dists = []
        if distance_matrix is None:
            for clade in self.tree.get_nonterminals():
                leaves = clade.get_terminals()
                if len(leaves) >= 2:
                    sum_dist = 0
                    count = 0
                    for i, leaf1 in enumerate(leaves):
                        for leaf2 in leaves[i + 1 :]:
                            try:
                                dist = self.tree.distance(leaf1, leaf2)
                                sum_dist += dist
                                count += 1
                            except Exception:
                                pass
                    if count > 0:
                        dists.append(sum_dist / count)
        else:
            terminal_index = {str(term): idx for idx, term in enumerate(self.terminals)}
            clade_indices = {}
            for clade in self.tree.get_nonterminals():
                leaves = clade.get_terminals()
                if len(leaves) >= 2:
                    indices = clade_indices.get(id(clade))
                    if indices is None:
                        indices = [
                            terminal_index[str(leaf)]
                            for leaf in leaves
                            if str(leaf) in terminal_index
                        ]
                        clade_indices[id(clade)] = indices
                    if len(indices) >= 2:
                        sum_dist = 0
                        count = 0
                        for i, t1 in enumerate(indices):
                            for t2 in indices[i + 1 :]:
                                dist = distance_matrix[t1, t2]
                                if np.isnan(dist):
                                    continue
                                sum_dist += dist
                                count += 1
                        if count > 0:
                            dists.append(sum_dist / count)
        if dists:
            q75, q25 = np.percentile(dists, [75, 25])
            iqr = q75 - q25
            threshold = np.median(dists) + conservative_factor * iqr
            return threshold
        return 0.1

    def _check_max_constraint(self, terminals, distance_matrix, threshold):
        max_dist = 0
        for i, t1 in enumerate(terminals):
            for t2 in terminals[i + 1 :]:
                max_dist = max(max_dist, distance_matrix[t1, t2])
        return max_dist <= threshold

    def _check_sum_constraint(self, terminals, distance_matrix, threshold):
        sum_dist = 0
        for i, t1 in enumerate(terminals):
            for t2 in terminals[i + 1 :]:
                sum_dist += distance_matrix[t1, t2]
        return sum_dist <= threshold

    def _check_avg_constraint(self, terminals, distance_matrix, threshold):
        sum_dist = 0
        count = 0
        for i, t1 in enumerate(terminals):
            for t2 in terminals[i + 1 :]:
                sum_dist += distance_matrix[t1, t2]
                count += 1
        if count == 0:
            return True
        return (sum_dist / count) <= threshold

    def _merge_clusters(self, terminals, clusters, cluster_map):
        if not terminals:
            return
        min_cluster = min(cluster_map[t] for t in terminals)
        for t in terminals:
            old_cluster = cluster_map[t]
            if old_cluster != min_cluster:
                for i, cluster_id in enumerate(cluster_map):
                    if cluster_id == old_cluster:
                        cluster_map[i] = min_cluster

    def is_monophyletic(self, cluster_terminals, tree=None):
        if tree is None:
            tree = self.tree
        try:
            mrca = tree.common_ancestor([str(t) for t in cluster_terminals])
            mrca_terminals = [str(t) for t in mrca.get_terminals()]
            return set(mrca_terminals) == set(str(t) for t in cluster_terminals)
        except Exception as e:
            print(f"Error checking monophyly: {e}")
            return False

    def refine_cluster(self, cluster_terminals, tree=None):
        if tree is None:
            tree = self.tree
        if self.is_monophyletic(cluster_terminals, tree):
            return [cluster_terminals]
        try:
            mrca = tree.common_ancestor([str(t) for t in cluster_terminals])
            internal_nodes = [clade for clade in mrca.get_nonterminals() if clade != mrca]
            refined_clusters = []
            for node in internal_nodes:
                node_terminals = [str(t) for t in node.get_terminals()]
                cluster_terminals_set = set(str(t) for t in cluster_terminals)
                if all(t in cluster_terminals_set for t in node_terminals):
                    already_included = False
                    for existing_cluster in refined_clusters:
                        if set(node_terminals).issubset(set(existing_cluster)):
                            already_included = True
                            break
                    if not already_included:
                        refined_clusters.append(node_terminals)
            all_refined = set()
            for cluster in refined_clusters:
                all_refined.update(cluster)
            for terminal in cluster_terminals:
                if str(terminal) not in all_refined:
                    refined_clusters.append([str(terminal)])
            return refined_clusters
        except Exception as e:
            print(f"Error refining cluster: {e}")
            return [cluster_terminals]

    def evaluate_monophyly(self, labels):
        unique_labels = np.unique(labels)
        n_clusters = len(unique_labels)
        monophyletic_count = 0
        non_monophyletic_clusters = []
        for label in unique_labels:
            cluster_indices = np.where(labels == label)[0]
            cluster_terminals = [self.terminals[i] for i in cluster_indices]
            if self.is_monophyletic(cluster_terminals):
                monophyletic_count += 1
            else:
                non_monophyletic_clusters.append(label)
        monophyletic_percentage = (monophyletic_count / n_clusters) * 100 if n_clusters > 0 else 0
        return {
            "total_clusters": n_clusters,
            "monophyletic_clusters": monophyletic_count,
            "monophyletic_percentage": monophyletic_percentage,
            "non_monophyletic_clusters": non_monophyletic_clusters,
        }

    def ensure_monophyletic_clusters(self, labels):
        """Refine clusters to ensure they are monophyletic."""
        unique_labels = np.unique(labels)
        refined_labels = labels.copy()
        for label in unique_labels:
            cluster_indices = np.where(labels == label)[0]
            cluster_terminals = [self.terminals[i] for i in cluster_indices]
            if not self.is_monophyletic(cluster_terminals):
                print(f"Refining non-monophyletic cluster {label}...")
                refined_subclusters = self.refine_cluster(cluster_terminals)
                next_label = max(unique_labels) + 1
                for idx, subcluster in enumerate(refined_subclusters):
                    if subcluster:
                        subcluster_indices = []
                        for terminal in subcluster:
                            for i, term in enumerate(self.terminals):
                                if str(term) == str(terminal):
                                    subcluster_indices.append(i)
                                    break
                        if subcluster_indices:
                            if idx == 0:
                                for i in subcluster_indices:
                                    refined_labels[i] = label
                            else:
                                for i in subcluster_indices:
                                    refined_labels[i] = next_label
                                next_label += 1
        return refined_labels

    def phydelity_clustering(self, distance_matrix):
        print("Performing Phydelity-inspired clustering...")
        n_terminals = len(self.terminals)
        flat_distances = []
        for i in range(n_terminals):
            for j in range(i + 1, n_terminals):
                flat_distances.append(distance_matrix[i, j])
        distances = np.array(flat_distances)
        median_dist = np.median(distances)

        q_distances = []
        for d in distances:
            q_distances.extend([abs(d - d2) for d2 in distances])
        q_distances.sort()
        k = int(len(q_distances) * 0.25)
        qn = q_distances[k] if k < len(q_distances) else 0
        mpl = median_dist + 3.5 * qn
        print(f"Maximal Patristic Distance Limit (MPL): {mpl:.4f}")

        potential_clusters = []
        cluster_labels = np.full(n_terminals, -1)

        for clade in self.tree.get_nonterminals():
            terminal_indices = []
            for leaf in clade.get_terminals():
                for i, term in enumerate(self.terminals):
                    if str(term) == str(leaf):
                        terminal_indices.append(i)
                        break
            if len(terminal_indices) < 2:
                continue

            mean_distance = 0
            count = 0
            for i, t1 in enumerate(terminal_indices):
                for t2 in terminal_indices[i + 1 :]:
                    mean_distance += distance_matrix[t1, t2]
                    count += 1
            if count > 0:
                mean_distance /= count
                if mean_distance <= mpl:
                    potential_clusters.append((clade, terminal_indices, mean_distance))

        potential_clusters.sort(key=lambda x: x[2])
        next_cluster_id = 0

        for _, terminal_indices, _ in potential_clusters:
            if any(cluster_labels[i] != -1 for i in terminal_indices):
                continue
            for i in terminal_indices:
                cluster_labels[i] = next_cluster_id
            next_cluster_id += 1

        for i in range(n_terminals):
            if cluster_labels[i] == -1:
                min_dist = float("inf")
                nearest_cluster = -1
                for j in range(n_terminals):
                    if cluster_labels[j] != -1 and distance_matrix[i, j] < min_dist:
                        min_dist = distance_matrix[i, j]
                        nearest_cluster = cluster_labels[j]
                if min_dist <= mpl and nearest_cluster != -1:
                    cluster_labels[i] = nearest_cluster
                else:
                    cluster_labels[i] = next_cluster_id
                    next_cluster_id += 1

        silhouette = -1
        if len(np.unique(cluster_labels)) > 1:
            try:
                silhouette = silhouette_score(distance_matrix, cluster_labels, metric="precomputed")
            except Exception as e:
                print(f"Error calculating silhouette score: {e}")

        print(
            f"Phydelity identified {len(np.unique(cluster_labels))} clusters with silhouette {silhouette:.2f}"
        )
        return cluster_labels, silhouette


###############################################################################
# 6. CLUSTERING MODULE (STANDARD)
###############################################################################
class ClusteringModule:
    """Performs ensemble clustering using KMeans, DBSCAN, and GaussianMixture."""

    def __init__(self, n_clusters_range=(2, 20), n_ensemble=10, dbscan_trials=50, seed=42):
        self.n_clusters_range = n_clusters_range
        self.n_ensemble = n_ensemble
        self.dbscan_trials = dbscan_trials
        self.seed = seed

    def ensemble_clustering(self, data):
        print("Optimizing DBSCAN parameters (Optuna)...")
        best_labels = None
        best_silhouette = -1
        try:
            best_params = self._optimize_dbscan(data)
            for _ in range(self.n_ensemble):
                for n_clusters in range(self.n_clusters_range[0], self.n_clusters_range[1] + 1):
                    methods = [
                        (KMeans, {"n_clusters": n_clusters, "random_state": self.seed}),
                        (
                            DBSCAN,
                            {
                                "eps": best_params.get("eps", 0.5),
                                "min_samples": best_params.get("min_samples", 5),
                            },
                        ),
                        (GaussianMixture, {"n_components": n_clusters, "random_state": self.seed}),
                    ]
                    for method, params in methods:
                        try:
                            clustering = method(**params)
                            labels = clustering.fit_predict(data)
                            if len(set(labels)) > 1:
                                sil = silhouette_score(data, labels)
                                if sil > best_silhouette:
                                    best_silhouette = sil
                                    # Start cluster numbering from 1
                                    best_labels = labels + 1
                        except Exception as e:
                            print(f"Method {method.__name__} failed: {e}")
                            continue
            print(f"Best silhouette score: {best_silhouette:.2f}")
            return best_labels, best_silhouette
        except Exception as e:
            print(f"Ensemble clustering failed: {e}")
            return None, None

    def _optimize_dbscan(self, data):
        def objective(trial):
            eps = trial.suggest_float("eps", 0.1, 2.0)
            min_samples = trial.suggest_int("min_samples", 3, 8)
            clustering = DBSCAN(eps=eps, min_samples=min_samples)
            labels = clustering.fit_predict(data)
            if len(set(labels)) <= 1:
                return -1.0
            return silhouette_score(data, labels)

        study = optuna.create_study(direction="maximize")
        study.optimize(objective, n_trials=self.dbscan_trials)
        return study.best_params

    def assign_outliers_to_clusters(self, embeddings, mask, labels):
        outlier_indices = np.where(~mask)[0]
        outlier_embeddings = embeddings[~mask]
        unique_labels = np.unique(labels)
        cluster_centers = np.array(
            [embeddings[mask][labels == lbl].mean(axis=0) for lbl in unique_labels]
        )
        outlier_assignments = []
        if len(outlier_embeddings) > 0 and len(cluster_centers) > 0:
            distances = cdist(outlier_embeddings, cluster_centers)
            for idx, distance in zip(outlier_indices, distances):
                closest_cluster = unique_labels[distance.argmin()]
                outlier_assignments.append((idx, closest_cluster))
        return outlier_assignments


###############################################################################
# 7. EVOLUTIONARY ANALYSIS
###############################################################################
class EvolutionaryAnalysis:
    """Handles evolutionary metrics: PD, mean pairwise distance, etc."""

    @staticmethod
    def evolutionary_cluster_analysis(tree, labels, strain_names, mask):
        """
        Compute evolutionary metrics for each phylogenetic cluster.

        Calculates phylogenetic diversity (PD), mean pairwise distance,
        and tree structure metrics for each cluster.

        Parameters
        ----------
        tree : Bio.Phylo.BaseTree.Tree
            Phylogenetic tree
        labels : np.ndarray
            Cluster labels for each strain
        strain_names : list
            List of strain names
        mask : np.ndarray
            Boolean mask for valid strains (inliers)

        Returns
        -------
        pd.DataFrame
            DataFrame with columns:
            - Cluster_ID: Cluster identifier
            - Strains: List of strain names in cluster
            - PD: Faith's Phylogenetic Diversity
            - MeanPairwiseDist: Mean patristic distance
            - InternalNodes: Number of internal nodes
            - SubtreeTerminals: Number of terminal nodes in subtree
        """
        valid_strains = np.array(strain_names)[mask]
        cluster_info = []
        for cluster_id in np.unique(labels):
            cluster_strains = valid_strains[labels == cluster_id]
            if len(cluster_strains) <= 1:
                cluster_info.append(
                    (cluster_id, cluster_strains.tolist(), 0.0, 0.0, 0, len(cluster_strains))
                )
            else:
                cluster_info.append(
                    EvolutionaryAnalysis._analyze_cluster_evolution(
                        tree, cluster_id, cluster_strains
                    )
                )
        df_clusters = pd.DataFrame(
            cluster_info,
            columns=[
                "Cluster_ID",
                "Strains",
                "PD",
                "MeanPairwiseDist",
                "InternalNodes",
                "SubtreeTerminals",
            ],
        )
        return df_clusters

    @staticmethod
    def _analyze_cluster_evolution(tree, cluster_id, cluster_strains):
        try:
            mrca = tree.common_ancestor(cluster_strains)
            pd_value = sum(clade.branch_length or 0 for clade in mrca.find_clades())
            distances = [
                tree.distance(s1, s2)
                for i, s1 in enumerate(cluster_strains)
                for s2 in cluster_strains[i + 1 :]
            ]
            mean_pairwise_dist = np.mean(distances) if distances else 0.0
            internal_nodes = sum(1 for node in mrca.find_clades() if not node.is_terminal())
            return (
                cluster_id,
                cluster_strains.tolist(),
                pd_value,
                mean_pairwise_dist,
                internal_nodes,
                len(mrca.get_terminals()),
            )
        except (ValueError, TypeError, AttributeError):
            return (cluster_id, cluster_strains.tolist(), np.nan, np.nan, np.nan, np.nan)

    @staticmethod
    def calculate_beta_diversity(tree, labels, strain_names, mask):
        """
        Calculate phylogenetic beta diversity between clusters.

        Computes UniFrac-like distances quantifying phylogenetic turnover
        between cluster pairs.

        Parameters
        ----------
        tree : Bio.Phylo.BaseTree.Tree
            Phylogenetic tree
        labels : np.ndarray
            Cluster labels
        strain_names : list
            Strain names
        mask : np.ndarray
            Boolean mask for valid strains

        Returns
        -------
        pd.DataFrame
            Symmetric matrix of beta diversity values between clusters

        Notes
        -----
        Beta diversity measures phylogenetic dissimilarity between communities.
        Values range from 0 (identical) to 1 (completely distinct).
        """
        valid_strains = np.array(strain_names)[mask]
        unique_labels = np.unique(labels)
        n_clusters = len(unique_labels)
        beta_div = np.zeros((n_clusters, n_clusters))
        for i, lbl1 in enumerate(unique_labels):
            s1 = valid_strains[labels == lbl1]
            for j, lbl2 in enumerate(unique_labels[i + 1 :], i + 1):
                s2 = valid_strains[labels == lbl2]
                distances = [tree.distance(a, b) for a in s1 for b in s2]
                beta_div[i, j] = beta_div[j, i] = np.mean(distances) if distances else 0.0
        return pd.DataFrame(beta_div, index=unique_labels, columns=unique_labels)

    @staticmethod
    def calculate_evolution_rates(cluster_df):
        rates = []
        for _, row in cluster_df.iterrows():
            rate = row["PD"] / row["InternalNodes"] if row["InternalNodes"] > 0 else 0.0
            rates.append({"Cluster_ID": row["Cluster_ID"], "EvolutionRate": rate})
        return pd.DataFrame(rates)

    def calculate_phylogenetic_signal_fritz_purvis(self, tree, trait_data):
        """
        Calculate Fritz-Purvis D statistic for binary traits.
        
        D = 0: trait conserved (follows phylogeny)
        D = 1: random distribution
        D < 0: overdispersed (distant taxa similar)
        D > 1: clustered (more than Brownian)
        
        Innovation: Novel application of Fritz-Purvis D to binary AMR/virulence traits.
        
        Args:
            tree: Phylogenetic tree (Bio.Phylo tree object)
            trait_data: DataFrame with binary traits (rows=strains, columns=traits)
        
        Returns:
            DataFrame with D statistics and interpretations
        """
        results = []
        
        for trait in trait_data.columns:
            trait_vector = trait_data[trait].values
            strain_names = trait_data.index.tolist() if hasattr(trait_data, 'index') else list(range(len(trait_vector)))
            
            # Calculate observed sister-clade differences
            scd_obs = self._sister_clade_differences(tree, strain_names, trait_vector)
            
            # Expected under random
            scd_random = self._expected_scd_random(trait_vector)
            
            # Expected under Brownian motion (for binary traits, use parsimony)
            scd_brownian = self._expected_scd_brownian(tree, strain_names, trait_vector)
            
            # D statistic
            if scd_random != scd_brownian:
                D = (scd_obs - scd_brownian) / (scd_random - scd_brownian)
            else:
                D = np.nan
            
            # Interpretation
            if np.isnan(D):
                interpretation = "Cannot calculate (insufficient variation)"
            elif D < 0:
                interpretation = "Overdispersed (distant taxa similar)"
            elif D < 0.5:
                interpretation = "Conserved (follows phylogeny)"
            elif D < 1.5:
                interpretation = "Random distribution"
            else:
                interpretation = "Clustered (more than Brownian)"
            
            results.append({
                'trait': trait,
                'd_statistic': round(D, 4) if not np.isnan(D) else np.nan,
                'interpretation': interpretation,
                'scd_observed': round(scd_obs, 4),
                'scd_random': round(scd_random, 4),
                'scd_brownian': round(scd_brownian, 4),
            })
        
        return pd.DataFrame(results)
    
    def _sister_clade_differences(self, tree, strain_names, trait_vector):
        """Calculate observed sister-clade differences."""
        # Count character changes on tree (parsimony)
        # Simple approach: count transitions along tree edges
        changes = 0
        
        # Map strain names to trait values
        trait_dict = {strain_names[i]: trait_vector[i] for i in range(len(strain_names))}
        
        # Traverse tree and count changes
        for clade in tree.find_clades():
            if clade.is_terminal():
                continue
            
            children = list(clade.clades)
            if len(children) >= 2:
                # Get trait values for children
                child_traits = []
                for child in children:
                    if child.is_terminal():
                        trait_val = trait_dict.get(str(child.name), None)
                        if trait_val is not None:
                            child_traits.append(trait_val)
                    else:
                        # For internal nodes, use majority trait of descendants
                        desc_traits = [trait_dict.get(str(term.name), None) 
                                      for term in child.get_terminals()]
                        desc_traits = [t for t in desc_traits if t is not None]
                        if desc_traits:
                            child_traits.append(int(np.mean(desc_traits) > 0.5))
                
                # Count differences between sister clades
                if len(child_traits) >= 2:
                    for i in range(len(child_traits) - 1):
                        if child_traits[i] != child_traits[i + 1]:
                            changes += 1
        
        return changes
    
    def _expected_scd_random(self, trait_vector):
        """Expected sister-clade differences under random distribution."""
        # Random: each taxon independently assigned trait
        # Expected changes = number of taxa - 1 (maximum possible)
        return len(trait_vector) - 1
    
    def _expected_scd_brownian(self, tree, strain_names, trait_vector):
        """Expected sister-clade differences under Brownian motion (parsimony minimum)."""
        # Minimum changes = number of unique trait states - 1
        unique_states = len(set(trait_vector))
        return max(unique_states - 1, 0)

    def calculate_phylogenetic_signal(cluster_df, output_folder):
        signals = []
        for _, row in cluster_df.iterrows():
            num_strains = len(row["Strains"])
            signal = row["PD"] / num_strains if num_strains > 0 else np.nan
            signals.append({"Cluster_ID": row["Cluster_ID"], "PhylogeneticSignal": signal})
        df_signal = pd.DataFrame(signals)
        df_signal["PhylogeneticSignal"] = df_signal["PhylogeneticSignal"].round(2)
        df_signal.to_csv(os.path.join(output_folder, "phylogenetic_signal.csv"), index=False)
        return df_signal


###############################################################################
# 8. DATA LOADER
###############################################################################
class DataLoader:
    """Loads and merges CSV data from clusters, MIC, AMR, Virulence."""

    def __init__(self, base_dir):
        if hasattr(base_dir, "base_dir"):
            base_dir = getattr(base_dir, "base_dir") or getattr(base_dir, "data_dir", base_dir)
        elif hasattr(base_dir, "data_dir"):
            base_dir = getattr(base_dir, "data_dir")
        self.base_dir = base_dir
        self.data_folder = base_dir

    def load_csv(self, filename):
        """Load a CSV file from the base directory."""
        path = Path(filename)
        if not path.is_absolute():
            path = Path(self.base_dir) / filename
        return pd.read_csv(path)

    def load_and_merge_data(self, clusters_path):
        """
        Load and merge cluster assignments with trait data.

        Combines cluster labels with MIC, AMR genes, and virulence factors
        into a single DataFrame for integrated analysis.

        Parameters
        ----------
        clusters_path : str
            Path to CSV file with cluster assignments (columns: Strain_ID, Cluster)

        Returns
        -------
        pd.DataFrame or None
            Merged DataFrame with columns:
            - Strain_ID: Strain identifier
            - Cluster: Cluster label
            - MIC columns: Antibiotic resistance phenotypes
            - AMR gene columns: Resistance genotypes
            - Virulence columns: Virulence factor presence
            Returns None if loading fails.

        Notes
        -----
        All binary trait columns are filled with 0 for missing values.
        Strain_ID is treated as string to preserve leading zeros.
        """
        try:
            mic_path = os.path.join(self.base_dir, "MIC.csv")
            amr_path = os.path.join(self.base_dir, "AMR_genes.csv")
            vir_path = os.path.join(self.base_dir, "Virulence.csv")

            clusters_df = pd.read_csv(clusters_path, dtype={"Strain_ID": str})
            mic_df = pd.read_csv(mic_path, dtype={"Strain_ID": str})
            amr_genes_df = pd.read_csv(amr_path, dtype={"Strain_ID": str})
            virulence_genes_df = pd.read_csv(vir_path, dtype={"Strain_ID": str})

            merged_df = clusters_df.copy()
            merged_df = pd.merge(merged_df, mic_df, on="Strain_ID", how="left")
            merged_df = pd.merge(merged_df, amr_genes_df, on="Strain_ID", how="left")
            merged_df = pd.merge(merged_df, virulence_genes_df, on="Strain_ID", how="left")

            binary_cols = [c for c in merged_df.columns if c not in ["Strain_ID", "Cluster"]]
            merged_df[binary_cols] = merged_df[binary_cols].fillna(0).astype(int)

            print(f"Data loaded successfully. Total strains: {len(merged_df)}")
            return merged_df
        except Exception as e:
            print(f"Error loading/merging data: {e}")
            return None


###############################################################################
# 9. VISUALIZER
###############################################################################
class Visualizer:
    """Creates PNG plots using Matplotlib/Seaborn."""

    def __init__(self, output_folder):
        self.output_folder = output_folder

    def plot_cluster_distribution(self, merged_df):
        """
        Plot and save cluster size distribution.

        Creates bar plot showing percentage of strains in each cluster
        with strain counts annotated.

        Parameters
        ----------
        merged_df : pd.DataFrame
            Merged data with Cluster and Strain_ID columns

        Returns
        -------
        pd.DataFrame
            Cluster statistics with columns:
            - Cluster: Cluster ID
            - Strain_Count: Number of strains
            - Percentage: Percentage of total strains

        Notes
        -----
        Saves both CSV and 300 DPI PNG to output folder.
        """
        cluster_stats = (
            merged_df.groupby("Cluster")
            .agg({"Strain_ID": ["count", lambda x: len(x) / len(merged_df) * 100]})
            .round(2)
        )
        cluster_stats.columns = ["Strain_Count", "Percentage"]
        cluster_stats = cluster_stats.sort_index()
        cluster_stats.to_csv(os.path.join(self.output_folder, "cluster_distribution.csv"))
        cluster_stats = cluster_stats.reset_index().rename(columns={"Cluster": "Cluster"})

        plt.figure(figsize=(12, 8))
        ax = sns.barplot(
            data=cluster_stats,
            x="Cluster",
            y="Percentage",
            hue="Cluster",
            palette="viridis",
            legend=False,
        )

        for i, p in enumerate(ax.patches):
            ax.annotate(
                f"{int(cluster_stats['Strain_Count'].iloc[i])}",
                (p.get_x() + p.get_width() / 2.0, p.get_height()),
                ha="center",
                va="bottom",
                fontsize=10,
            )

        plt.title("Cluster Distribution (%)", fontsize=16)
        plt.xlabel("Cluster", fontsize=14)
        plt.ylabel("Percentage of Strains", fontsize=14)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_folder, "cluster_distribution.png"), dpi=300)
        plt.close()
        return cluster_stats

    def plot_phylogenetic_tree(
        self, tree, labels, strain_names, mask, output_file="phylogenetic_tree.png"
    ):
        """
        Plot phylogenetic tree with color-coded clusters.

        Generates high-resolution tree visualization with clusters distinguished
        by color for both terminal labels and branches.

        Parameters
        ----------
        tree : Bio.Phylo.BaseTree.Tree
            Phylogenetic tree object
        labels : np.ndarray
            Cluster labels for each strain
        strain_names : list
            List of strain names
        mask : np.ndarray
            Boolean mask for valid strains
        output_file : str, default="phylogenetic_tree.png"
            Output filename

        Notes
        -----
        Saves 300 DPI PNG image to output folder.
        Uses viridis colormap for cluster colors.
        Includes legend showing cluster colors.
        """
        from matplotlib.patches import Patch

        # Ensure cluster numbering starts from 1
        unique_labels = np.unique(labels)
        colors = plt.cm.viridis(np.linspace(0, 1, len(unique_labels)))

        strain_colors = {}
        valid_strains = np.array(strain_names)[mask]
        for i, label in enumerate(unique_labels):
            cluster_strains = valid_strains[labels == label]
            for strain in cluster_strains:
                strain_colors[strain] = colors[i]

        plt.figure(figsize=(15, 10))
        ax = plt.axes()

        # Draw the tree but store the axes
        Phylo.draw(
            tree,
            axes=ax,
            do_show=False,
            label_func=lambda x: x.name if x.name in strain_colors else None,
        )

        # Process text objects to ensure correct coloring
        for text in ax.texts:
            text_content = text.get_text().strip()
            if text_content in strain_colors:
                text.set_color(strain_colors[text_content])
                text.set_fontweight("bold")
                text.set_fontsize(9)

        # Process line objects to color branches
        for line in ax.get_lines():
            line_label = line.get_label()
            if line_label in strain_colors:
                line.set_color(strain_colors[line_label])
                line.set_linewidth(2.5)
                line.set_alpha(0.8)

        # Create legend
        legend_elements = [
            Patch(facecolor=colors[i], label=f"Cluster {lbl}")
            for i, lbl in enumerate(unique_labels)
        ]
        plt.legend(handles=legend_elements, loc="upper right", title="Clusters")
        plt.title("Phylogenetic Tree with Clusters (matplotlib)", fontsize=16)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_folder, output_file), dpi=300)
        plt.close()

    def plot_umap_clusters(self, embeddings, labels, mask, outlier_assignments=None):
        plt.figure(figsize=(12, 10))
        unique_labels = np.unique(labels)
        colors = plt.cm.viridis(np.linspace(0, 1, len(unique_labels)))

        for i, label in enumerate(unique_labels):
            idx = labels == label
            plt.scatter(
                embeddings[mask][idx, 0],
                embeddings[mask][idx, 1],
                c=[colors[i]],
                s=50,
                label=f"Cluster {label}",
                alpha=0.8,
            )

        if outlier_assignments:
            outlier_points = []
            outlier_colors = []
            for idx, cluster in outlier_assignments:
                outlier_points.append(embeddings[idx])
                cluster_idx = np.where(unique_labels == cluster)[0][0]
                outlier_colors.append(colors[cluster_idx])
            if outlier_points:
                outlier_points = np.array(outlier_points)
                plt.scatter(
                    outlier_points[:, 0],
                    outlier_points[:, 1],
                    c=outlier_colors,
                    marker="x",
                    s=100,
                    linewidths=1,
                    edgecolors="black",
                    label="Outliers",
                )

        plt.title("UMAP Clusters (matplotlib)", fontsize=16)
        plt.xlabel("UMAP Dim 1", fontsize=14)
        plt.ylabel("UMAP Dim 2", fontsize=14)
        plt.legend(title="Clusters", title_fontsize=12)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_folder, "umap_clusters.png"), dpi=300)
        plt.close()


###############################################################################
# 10. TRAIT ANALYZER
###############################################################################
class TraitAnalyzer:
    """Analyzes binary traits in the context of clusters."""

    def __init__(self, data=None, output_folder="output"):
        if isinstance(data, str) and output_folder == "output":
            output_folder = data
            data = None
        self.data = data
        self.output_folder = output_folder

    def analyze_all_categories(self, merged_df):
        categories = {
            "MIC": [col for col in merged_df.columns if col.lower().startswith("mic_")],
            "AMR": [col for col in merged_df.columns if col.lower().startswith("amr_")],
            "Virulence": [col for col in merged_df.columns if col.lower().startswith("vir_")],
        }
        combined_html = ""
        for cat_name, cols in categories.items():
            if not cols:
                continue
            print(f"\nCategory Analysis: {cat_name}")
            cat_df = merged_df[["Strain_ID", "Cluster"] + cols].copy()
            self._analyze_trait_frequencies(cat_df, cat_name)
            self._perform_statistical_tests(cat_df, cat_name)
            self._calculate_feature_importance(cat_df, cat_name)
            freq_csv = os.path.join(self.output_folder, f"trait_frequencies_{cat_name.lower()}.csv")
            if os.path.exists(freq_csv):
                df_freq = pd.read_csv(freq_csv)
                combined_html += (
                    f"<h4>{cat_name} Trait Frequencies</h4>"
                    + self._df_to_interactive_table(
                        df_freq,
                        table_id=f"table-{cat_name}-freq",
                        filename=f"{cat_name}_trait_frequencies",
                    )
                )

        with open(os.path.join(self.output_folder, "trait_analysis_summary.html"), "w") as f:
            f.write(combined_html)

        return combined_html

    def _analyze_trait_frequencies(self, df, suffix=""):
        binary_columns = [col for col in df.columns if col not in ["Strain_ID", "Cluster"]]
        trait_freq = (df.groupby("Cluster")[binary_columns].mean() * 100).round(2)
        trait_freq.to_csv(
            os.path.join(self.output_folder, f"trait_frequencies_{suffix.lower()}.csv")
        )
        plt.figure(figsize=(10, 6))
        sns.heatmap(trait_freq, annot=True, fmt=".2f", cmap="viridis")
        plt.title(f"Binary Trait Frequencies (%) - {suffix}")
        plt.ylabel("Cluster")
        plt.xlabel("Trait")
        plt.tight_layout()
        plt.savefig(
            os.path.join(self.output_folder, f"trait_frequencies_heatmap_{suffix.lower()}.png")
        )
        plt.close()

    def _perform_statistical_tests(self, df, suffix="", alpha=0.05):
        binary_columns = [col for col in df.columns if col not in ["Strain_ID", "Cluster"]]
        results = []
        for trait in binary_columns:
            table = pd.crosstab(df["Cluster"], df[trait])
            chi2, pval, *_ = chi2_contingency(table)
            results.append({"Feature": trait, "Chi2": chi2, "p_value": pval})
        results_df = pd.DataFrame(results)
        _, results_df["p_adjusted"], _, _ = multipletests(
            results_df["p_value"], alpha=alpha, method="fdr_bh"
        )
        results_df["significant"] = results_df["p_adjusted"] < alpha
        results_df.to_csv(
            os.path.join(self.output_folder, f"statistical_tests_{suffix.lower()}.csv"), index=False
        )

    def _calculate_feature_importance(self, df, suffix="", n_bootstrap=1000):
        binary_columns = [col for col in df.columns if col not in ["Strain_ID", "Cluster"]]
        X = df[binary_columns].values
        y = df["Cluster"].values
        importance_scores = []

        for _ in tqdm(range(n_bootstrap), desc=f"Bootstrap RF - {suffix}"):
            indices = np.random.choice(len(X), len(X), replace=True)
            X_boot, y_boot = X[indices], y[indices]
            rf = RandomForestClassifier(n_estimators=100, random_state=42)
            rf.fit(X_boot, y_boot)
            importance_scores.append(rf.feature_importances_)

        importance_df = pd.DataFrame(
            {
                "Feature": binary_columns,
                "Importance_Mean": np.mean(importance_scores, axis=0),
                "Importance_Std": np.std(importance_scores, axis=0),
                "Importance_Lower": np.percentile(importance_scores, 2.5, axis=0),
                "Importance_Upper": np.percentile(importance_scores, 97.5, axis=0),
            }
        )
        importance_df = importance_df.sort_values("Importance_Mean", ascending=False)
        importance_df.to_csv(
            os.path.join(self.output_folder, f"feature_importance_{suffix.lower()}.csv"),
            index=False,
        )

        plt.figure(figsize=(10, 8))
        sns.barplot(
            data=importance_df.head(20),
            x="Importance_Mean",
            y="Feature",
            xerr=importance_df.head(20)["Importance_Std"],
        )
        plt.title(f"Top 20 Features (RandomForest) - {suffix}")
        plt.xlabel("Feature Importance (mean ± std)")
        plt.tight_layout()
        plt.savefig(
            os.path.join(self.output_folder, f"feature_importance_{suffix.lower()}.png"), dpi=300
        )
        plt.close()

    def bootstrap_feature_importance(self, df, n_bootstrap=1000):
        binary_columns = [col for col in df.columns if col not in ["Strain_ID", "Cluster"]]
        X = df[binary_columns].values
        y = df["Cluster"].values
        importance_scores = []
        if len(df) > LARGE_DATASET_THRESHOLD:
            n_bootstrap = min(n_bootstrap, MAX_BOOTSTRAP_LARGE_DATASET)
            n_estimators = max(
                MIN_ESTIMATORS,
                BASE_ESTIMATORS - (len(df) - LARGE_DATASET_THRESHOLD) // ESTIMATOR_REDUCTION_FACTOR,
            )
        else:
            n_estimators = BASE_ESTIMATORS

        for _ in tqdm(range(n_bootstrap), desc="Bootstrap Feature Importance"):
            indices = np.random.choice(len(X), len(X), replace=True)
            X_boot, y_boot = X[indices], y[indices]
            rf = RandomForestClassifier(n_estimators=n_estimators, random_state=42)
            rf.fit(X_boot, y_boot)
            importance_scores.append(rf.feature_importances_)

        df_scores = pd.DataFrame(
            {
                "Feature": binary_columns,
                "Importance_Mean": np.mean(importance_scores, axis=0),
                "Importance_Std": np.std(importance_scores, axis=0),
                "Importance_Lower": np.percentile(importance_scores, 2.5, axis=0),
                "Importance_Upper": np.percentile(importance_scores, 97.5, axis=0),
            }
        )
        df_scores = df_scores.sort_values("Importance_Mean", ascending=False)
        df_scores.to_csv(
            os.path.join(self.output_folder, "bootstrap_feature_importance.csv"), index=False
        )
        return df_scores

    def log_odds_ratio_analysis(self, df):
        traits = [col for col in df.columns if col not in ["Strain_ID", "Cluster"]]
        global_data = []
        cluster_data = []
        for trait in traits:
            a = df[trait].sum()
            b = len(df) - a
            global_log_odds = np.log((a + 0.5) / (b + 0.5))
            global_data.append({"Feature": trait, "Log_Odds_Ratio": round(global_log_odds, 2)})

            for cl in sorted(df["Cluster"].unique()):
                subset = df[df["Cluster"] == cl]
                a_c = subset[trait].sum()
                b_c = len(subset) - a_c
                log_odds = np.log((a_c + 0.5) / (b_c + 0.5))
                cluster_data.append(
                    {"Cluster": cl, "Feature": trait, "Log_Odds_Ratio": round(log_odds, 2)}
                )

        df_global = pd.DataFrame(global_data)
        df_cluster = pd.DataFrame(cluster_data)
        df_global.to_csv(os.path.join(self.output_folder, "log_odds_global.csv"), index=False)
        df_cluster.to_csv(os.path.join(self.output_folder, "log_odds_per_cluster.csv"), index=False)
        return df_global, df_cluster

    def association_rule_mining(self, df, min_support=0.05, min_confidence=0.7, max_features=100):
        import traceback

        try:
            binary_cols = [col for col in df.columns if col not in ["Strain_ID", "Cluster"]]
            feature_counts = df[binary_cols].sum().sort_values(ascending=False)
            selected_features = feature_counts.index[:max_features].tolist()
            if len(df) > LARGE_DATASET_THRESHOLD and len(selected_features) > LARGE_FEATURE_THRESHOLD:
                selected_features = selected_features[:LARGE_FEATURE_THRESHOLD]
            print(f"Selected {len(selected_features)}/{len(binary_cols)} most frequent features")
            trait_df = df[selected_features].copy()

            max_len = (
                MAX_LEN_LARGE_DATASET
                if len(trait_df.columns) > LARGE_FEATURE_THRESHOLD
                else MAX_LEN_NORMAL
            )
            if len(trait_df.columns) > LARGE_FEATURE_THRESHOLD:
                adaptive_min_support = max(min_support, 0.1)
                print(f"Adapting min_support to {adaptive_min_support} for large dataset")
            else:
                adaptive_min_support = min_support
            max_itemsets = (
                MAX_ITEMSETS_LIMIT
                if len(df) > LARGE_DATASET_THRESHOLD or len(trait_df.columns) > LARGE_FEATURE_THRESHOLD
                else None
            )

            while len(selected_features) > 10:
                try:
                    print(
                        f"Attempting apriori with {len(selected_features)} features, min_support={adaptive_min_support}"
                    )
                    frequent_itemsets = apriori(
                        trait_df,
                        min_support=adaptive_min_support,
                        use_colnames=True,
                        max_len=max_len,
                    )
                    if max_itemsets and len(frequent_itemsets) > max_itemsets:
                        # Bound itemset volume for large datasets to keep runtime predictable.
                        frequent_itemsets = frequent_itemsets.nlargest(max_itemsets, "support")
                    break
                except MemoryError:
                    selected_features = selected_features[: len(selected_features) // 2]
                    trait_df = df[selected_features].copy()
                    adaptive_min_support += 0.05
                    print(
                        f"MemoryError - reducing to {len(selected_features)} features, min_support={adaptive_min_support}"
                    )

            if "frequent_itemsets" not in locals() or frequent_itemsets.empty:
                print("No frequent itemsets found with the specified support threshold.")
                cols = ["Antecedent", "Consequent", "Support", "Confidence", "Lift"]
                df_rules = pd.DataFrame(columns=cols)
                df_rules.to_csv(
                    os.path.join(self.output_folder, "association_rules.csv"), index=False
                )
                return df_rules

            rules = association_rules(
                frequent_itemsets, metric="confidence", min_threshold=min_confidence
            )
            if rules.empty:
                print("No rules generated that meet the confidence threshold.")
                cols = ["Antecedent", "Consequent", "Support", "Confidence", "Lift"]
                df_rules = pd.DataFrame(columns=cols)
            else:
                rules["Antecedent"] = rules["antecedents"].apply(lambda x: ", ".join(list(x)))
                rules["Consequent"] = rules["consequents"].apply(lambda x: ", ".join(list(x)))
                df_rules = rules[["Antecedent", "Consequent", "support", "confidence", "lift"]]
                df_rules.columns = ["Antecedent", "Consequent", "Support", "Confidence", "Lift"]
                df_rules = df_rules.round(4)

            df_rules.to_csv(os.path.join(self.output_folder, "association_rules.csv"), index=False)
            return df_rules
        except Exception as e:
            print(f"Error in association rule mining: {e}")
            traceback.print_exc()
            cols = ["Antecedent", "Consequent", "Support", "Confidence", "Lift"]
            df_rules = pd.DataFrame(columns=cols)
            df_rules.to_csv(os.path.join(self.output_folder, "association_rules.csv"), index=False)
            return df_rules

    def label_shared_unique_features(self, df, presence_threshold=0.3):
        import traceback

        try:
            if "Cluster" not in df.columns:
                print("Error: 'Cluster' column not found in dataframe.")
                return (
                    pd.DataFrame(
                        columns=[
                            "Feature",
                            "Clusters",
                            "NumClusters",
                            "Count",
                            "Percent_in_Clusters",
                        ]
                    ),
                    pd.DataFrame(columns=["Cluster", "Feature", "Count"]),
                )

            print(f"Analyzing shared/unique features among {len(df['Cluster'].unique())} clusters")
            binary_cols = [col for col in df.columns if col not in ["Strain_ID", "Cluster"]]
            print(f"Found {len(binary_cols)} binary trait columns")

            feature_presence = {}
            clusters = sorted(df["Cluster"].unique())
            for cluster in clusters:
                cluster_df = df[df["Cluster"] == cluster]
                n_samples = len(cluster_df)
                if n_samples == 0:
                    continue
                present_features = []
                for feature in binary_cols:
                    feature_sum = cluster_df[feature].sum()
                    presence = feature_sum / n_samples if n_samples > 0 else 0
                    if feature not in feature_presence:
                        feature_presence[feature] = {}
                    feature_presence[feature][cluster] = presence
                    if presence >= presence_threshold:
                        present_features.append(feature)
                print(f"Cluster {cluster}: {len(present_features)} features above threshold")

            shared_features = []
            unique_features = []

            for feature, presence_map in feature_presence.items():
                present_in = [c for c, v in presence_map.items() if v >= presence_threshold]
                if len(present_in) > 1:
                    total_count = df[feature].sum()
                    percent_in_clusters = total_count / len(df) * 100
                    shared_features.append(
                        {
                            "Feature": feature,
                            "Clusters": ", ".join(map(str, present_in)),
                            "NumClusters": len(present_in),
                            "Count": int(total_count),
                            "Percent_in_Clusters": round(percent_in_clusters, 2),
                        }
                    )
                elif len(present_in) == 1:
                    cluster_val = present_in[0]
                    cluster_df = df[df["Cluster"] == cluster_val]
                    count = int(cluster_df[feature].sum())
                    unique_features.append(
                        {"Cluster": cluster_val, "Feature": feature, "Count": count}
                    )

            df_shared = pd.DataFrame(shared_features)
            df_unique = pd.DataFrame(unique_features)

            if not df_shared.empty:
                df_shared = df_shared.sort_values(
                    by=["NumClusters", "Count"], ascending=[False, False]
                )
                print(f"Found {len(df_shared)} shared features")
            else:
                print("No shared features found!")

            if not df_unique.empty:
                df_unique = df_unique.sort_values(by=["Cluster", "Count"], ascending=[True, False])
                print(f"Found {len(df_unique)} unique features")
            else:
                print("No unique features found!")

            df_shared.to_csv(os.path.join(self.output_folder, "shared_features.csv"), index=False)
            df_unique.to_csv(os.path.join(self.output_folder, "unique_features.csv"), index=False)
            return df_shared, df_unique

        except Exception as e:
            print(f"Error in shared/unique feature analysis: {e}")
            traceback.print_exc()
            return (
                pd.DataFrame(
                    columns=["Feature", "Clusters", "NumClusters", "Count", "Percent_in_Clusters"]
                ),
                pd.DataFrame(columns=["Cluster", "Feature", "Count"]),
            )

    def pairwise_fdr_post_hoc(self, df, alpha=0.05, method="fdr_bh"):
        import traceback

        try:
            binary_cols = [col for col in df.columns if col not in ["Strain_ID", "Cluster"]]
            clusters = sorted(df["Cluster"].unique())
            print(
                f"Performing pairwise post-hoc analysis with {len(clusters)} clusters and {len(binary_cols)} features"
            )
            if len(clusters) < 2:
                print("Not enough clusters for pairwise analysis.")
                df_posthoc = pd.DataFrame(
                    columns=[
                        "Feature",
                        "ClusterA",
                        "ClusterB",
                        "p_value_raw",
                        "p_value_adj",
                        "significant",
                    ]
                )
                df_posthoc.to_csv(
                    os.path.join(self.output_folder, "pairwise_fdr_post_hoc.csv"), index=False
                )
                return df_posthoc

            results = []
            all_pvals = []
            test_indices = []
            active_features = [feat for feat in binary_cols if df[feat].sum() > 0]
            if len(df) > LARGE_DATASET_THRESHOLD and len(active_features) > LARGE_FEATURE_THRESHOLD:
                active_features = active_features[:LARGE_FEATURE_THRESHOLD]
            print(f"Found {len(active_features)} active features for analysis")

            for feature in active_features:
                for i, j in combinations(clusters, 2):
                    cluster_i = df[df["Cluster"] == i]
                    cluster_j = df[df["Cluster"] == j]
                    if len(cluster_i) == 0 or len(cluster_j) == 0:
                        continue
                    data_i = cluster_i[feature]
                    data_j = cluster_j[feature]
                    table = [
                        [int(data_i.sum()), int(data_j.sum())],
                        [len(data_i) - int(data_i.sum()), len(data_j) - int(data_j.sum())],
                    ]
                    if (table[0][0] + table[0][1] == 0) or (table[1][0] + table[1][1] == 0):
                        continue
                    try:
                        _, p_value = fisher_exact(table)
                        results.append(
                            {
                                "Feature": feature,
                                "ClusterA": int(i),
                                "ClusterB": int(j),
                                "p_value_raw": float(p_value),
                            }
                        )
                        all_pvals.append(p_value)
                        test_indices.append(len(results) - 1)
                    except Exception as e:
                        print(f"Error in Fisher test for {feature}, clusters {i}-{j}: {e}")

            print(f"Completed {len(results)} pairwise tests")
            if all_pvals:
                rejected, p_adjusted, _, _ = multipletests(all_pvals, alpha=alpha, method=method)
                for idx, adj_p in zip(test_indices, p_adjusted):
                    results[idx]["p_value_adj"] = float(adj_p)
                    results[idx]["significant"] = bool(adj_p < alpha)

            df_posthoc = pd.DataFrame(results)
            if not df_posthoc.empty:
                df_posthoc = df_posthoc.sort_values(
                    by=["p_value_adj", "Feature", "ClusterA", "ClusterB"]
                )
                df_posthoc = df_posthoc.round(4)
                print(
                    f"Found {sum(df_posthoc['significant'])} significant relationships after FDR correction"
                )
            else:
                print("No significant pairwise relationships found")

            df_posthoc.to_csv(
                os.path.join(self.output_folder, "pairwise_fdr_post_hoc.csv"), index=False
            )
            return df_posthoc

        except Exception as e:
            print(f"Error in pairwise FDR post-hoc analysis: {e}")
            traceback.print_exc()
            df_posthoc = pd.DataFrame(
                columns=[
                    "Feature",
                    "ClusterA",
                    "ClusterB",
                    "p_value_raw",
                    "p_value_adj",
                    "significant",
                ]
            )
            df_posthoc.to_csv(
                os.path.join(self.output_folder, "pairwise_fdr_post_hoc.csv"), index=False
            )
            return df_posthoc

    def bootstrap_log_odds(self, df, n_bootstrap=100):
        import traceback

        try:
            binary_cols = [col for col in df.columns if col not in ["Strain_ID", "Cluster"]]
            active_traits = [col for col in binary_cols if df[col].sum() > 0]
            if len(df) > LARGE_DATASET_THRESHOLD and len(active_traits) > MAX_ACTIVE_TRAITS_LARGE_DATASET:
                active_traits = active_traits[:MAX_ACTIVE_TRAITS_LARGE_DATASET]
            if len(df) > LARGE_DATASET_THRESHOLD:
                n_bootstrap = min(n_bootstrap, MAX_BOOTSTRAP_LOG_ODDS_LARGE_DATASET)
            print(f"Running bootstrap log-odds analysis on {len(active_traits)} active traits")
            results = []
            for trait in active_traits:
                bootstrap_log_odds_values = []
                for i in range(n_bootstrap):
                    bootstrap_indices = np.random.choice(len(df), size=len(df), replace=True)
                    bootstrap_df = df.iloc[bootstrap_indices]
                    a = bootstrap_df[trait].sum()
                    b = len(bootstrap_df) - a
                    log_odds = np.log((a + 0.5) / (b + 0.5))
                    bootstrap_log_odds_values.append(log_odds)
                arr = np.array(bootstrap_log_odds_values)
                mean_log_odds = np.mean(arr)
                std_log_odds = np.std(arr)
                lower_ci = np.percentile(arr, 2.5)
                upper_ci = np.percentile(arr, 97.5)
                results.append(
                    {
                        "Feature": trait,
                        "Bootstrap_Log_Odds_Mean": mean_log_odds,
                        "Bootstrap_Log_Odds_Std": std_log_odds,
                        "Bootstrap_Log_Odds_Lower_CI": lower_ci,
                        "Bootstrap_Log_Odds_Upper_CI": upper_ci,
                        "Significant": (lower_ci * upper_ci) > 0,
                    }
                )

            df_boot = pd.DataFrame(results)
            if not df_boot.empty:
                df_boot["Abs_Log_Odds"] = df_boot["Bootstrap_Log_Odds_Mean"].abs()
                df_boot = df_boot.sort_values("Abs_Log_Odds", ascending=False)
                df_boot.drop(columns=["Abs_Log_Odds"], inplace=True)
                df_boot = df_boot.round(4)
                print(f"Found {len(df_boot)} bootstrap log-odds results")
            else:
                print("No bootstrap log-odds results generated")

            df_boot.to_csv(os.path.join(self.output_folder, "bootstrap_log_odds.csv"), index=False)
            return df_boot

        except Exception as e:
            print(f"Error in bootstrap log-odds analysis: {e}")
            traceback.print_exc()
            df_boot = pd.DataFrame(
                columns=[
                    "Feature",
                    "Bootstrap_Log_Odds_Mean",
                    "Bootstrap_Log_Odds_Std",
                    "Bootstrap_Log_Odds_Lower_CI",
                    "Bootstrap_Log_Odds_Upper_CI",
                    "Significant",
                ]
            )
            df_boot.to_csv(os.path.join(self.output_folder, "bootstrap_log_odds.csv"), index=False)
            return df_boot

    def _df_to_interactive_table(self, df, table_id="table-default", filename="output"):
        if df.empty:
            return "<p>No data available.</p>"
        html = f"""
        <div class="mb-2">
            <button class="btn btn-sm btn-outline-primary" onclick="exportTableToCSV('{table_id}', '{filename}.csv')">Export to CSV</button>
        </div>
        <table id="{table_id}" class="table table-striped table-hover dataTable" border="0">
        <thead><tr>"""
        for col in df.columns:
            html += f"<th>{col}</th>"
        html += "</tr></thead><tbody>"

        for _, row in df.iterrows():
            html += "<tr>"
            for val in row:
                html += f"<td>{val}</td>"
            html += "</tr>"
        html += "</tbody></table>"
        return html


###############################################################################
# 11. MCA ANALYZER
###############################################################################
class MCAAnalyzer:
    """Performs Multiple Correspondence Analysis (MCA)."""

    def __init__(self, data=None, output_folder="output"):
        if isinstance(data, str) and output_folder == "output":
            output_folder = data
            data = None
        self.data = data
        self.output_folder = output_folder

    def perform_mca_analysis(self, merged_df):
        try:
            binary_cols = [c for c in merged_df.columns if c not in ["Strain_ID", "Cluster"]]
            if not binary_cols:
                print("No binary columns for MCA.")
                return None, None, None

            mca_data = merged_df[binary_cols].astype("category")
            clusters = merged_df["Cluster"]
            mca = prince.MCA(n_components=2, random_state=42)
            mca.fit(mca_data)

            row_coords = mca.row_coordinates(mca_data)
            row_coords.columns = ["Component_1", "Component_2"]
            row_coords["Cluster"] = clusters

            col_coords = mca.column_coordinates(mca_data)
            col_coords.columns = ["Component_1", "Component_2"]
            col_coords["Feature_Type"] = "Other"
            for feature in col_coords.index:
                if feature.lower().startswith("mic_"):
                    col_coords.loc[feature, "Feature_Type"] = "MIC"
                elif feature.lower().startswith("amr_"):
                    col_coords.loc[feature, "Feature_Type"] = "AMR"
                elif feature.lower().startswith("vir_"):
                    col_coords.loc[feature, "Feature_Type"] = "Virulence"

            eigenvalues = mca.eigenvalues_
            total_inertia = sum(eigenvalues)
            explained_inertia = [x / total_inertia for x in eigenvalues]

            row_coords.to_csv(
                os.path.join(self.output_folder, "mca_row_coordinates.csv"), index=False
            )
            col_coords.to_csv(
                os.path.join(self.output_folder, "mca_column_coordinates.csv"), index=False
            )

            summary_df = pd.DataFrame(
                {
                    "Component": range(1, len(eigenvalues) + 1),
                    "Eigenvalue": [round(x, 4) for x in eigenvalues],
                    "Explained_Inertia": [round(x, 4) for x in explained_inertia],
                    "Cumulative_Explained_Inertia": [
                        round(sum(explained_inertia[: i + 1]), 4)
                        for i in range(len(explained_inertia))
                    ],
                }
            )
            summary_df.to_csv(os.path.join(self.output_folder, "mca_summary.csv"), index=False)

            mca_guide = pd.DataFrame(
                {
                    "Aspect": [
                        "Column Points",
                        "Row Points",
                        "Proximity (Columns)",
                        "Proximity (Rows)",
                        "Row-Column Relationship",
                        "Component 1",
                        "Component 2",
                    ],
                    "Interpretation": [
                        "Column points represent binary features (MIC, AMR, Virulence genes). Their position shows correlation with other features.",
                        "Row points represent strains. Proximity between points indicates similar trait profiles.",
                        "Features close together tend to co-occur in the same strains.",
                        "Strains close together share similar genetic trait profiles.",
                        "Strains positioned near certain features are more likely to possess those features.",
                        f"Explains {round(explained_inertia[0]*100, 2)}% of variation in the data.",
                        f"Explains {round(explained_inertia[1]*100, 2)}% of variation in the data.",
                    ],
                }
            )
            mca_guide.to_csv(
                os.path.join(self.output_folder, "mca_interpretation_guide.csv"), index=False
            )

            return row_coords, col_coords, summary_df
        except Exception as e:
            print(f"Error in MCA analysis: {e}")
            return None, None, None


###############################################################################
# 12. ADVANCED HTML REPORT GENERATOR
###############################################################################
class HTMLReportGenerator:
    """
    Generates the final HTML report:
    - Rounds numeric columns,
    - Uses DataTables for interactive tables,
    - Integrates Plotly charts,
    - Automatically scans output folder for additional files,
    - Creates a static dendrogram with color-coded clusters.
    """

    def __init__(self, output_folder, base_dir="."):
        self.output_folder = output_folder
        self.base_dir = base_dir
        create_template_directory(base_dir=self.base_dir)
        self.template_env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(searchpath=os.path.join(self.base_dir, "templates")),
            autoescape=jinja2.select_autoescape(["html", "xml"]),
        )

    def _round_numeric_columns(self, df):
        if df is None or df.empty:
            return df
        df_copy = df.copy()
        for col in df_copy.columns:
            if col.lower() != "strain_id":
                if pd.api.types.is_numeric_dtype(df_copy[col]):
                    df_copy[col] = df_copy[col].round(2)
        return df_copy

    def _df_to_interactive_table(self, df, table_id="table-default", filename="output"):
        if df is None or df.empty:
            return "<p>No data available.</p>"
        df_rounded = self._round_numeric_columns(df)
        html = f"""
        <div class="mb-2">
            <button class="btn btn-sm btn-outline-primary" onclick="exportTableToCSV('{table_id}', '{filename}.csv')">Export to CSV</button>
        </div>
        <table id="{table_id}" class="table table-striped table-hover dataTable" border="0">
        <thead><tr>"""
        for col in df_rounded.columns:
            html += f"<th>{col}</th>"
        html += "</tr></thead><tbody>"
        for _, row in df_rounded.iterrows():
            html += "<tr>"
            for val in row:
                html += f"<td>{val}</td>"
            html += "</tr>"
        html += "</tbody></table>"
        return html

    def _create_interactive_tree_plot(self, tree_path):
        """
        Updated to ensure color-coded strain labels and branches in the HTML (static) figure.
        """
        try:
            tree = Phylo.read(tree_path, "newick")
            clusters_file = os.path.join(self.output_folder, "phylogenetic_clusters.csv")
            cluster_map = {}

            if os.path.exists(clusters_file):
                clusters_df = pd.read_csv(clusters_file)
                for _, row in clusters_df.iterrows():
                    strain_id = str(row["Strain_ID"]).strip()
                    cluster_map[strain_id] = int(row["Cluster"])

            leaves = tree.get_terminals()
            leaf_count = len(leaves)
            height = max(10, leaf_count * 0.25)
            width = 14

            plt.figure(figsize=(width, height), dpi=120)
            ax = plt.subplot(1, 1, 1)

            unique_clusters = sorted(set(cluster_map.values())) if cluster_map else []
            color_map = plt.get_cmap("tab10", max(10, len(unique_clusters)))
            cluster_colors = {c: color_map(i % 10) for i, c in enumerate(unique_clusters)}

            # Draw the tree with explicit label function
            Phylo.draw(
                tree,
                axes=ax,
                do_show=False,
                branch_labels=None,
                label_func=lambda x: (
                    x.name if (x.is_terminal() and x.name in cluster_map) else None
                ),
            )

            # Color the text for each leaf according to cluster
            for text in ax.texts:
                node_name = text.get_text().strip()
                if node_name in cluster_map:
                    cluster_val = cluster_map[node_name]
                    color = cluster_colors.get(cluster_val, "black")
                    text.set_color(color)
                    text.set_fontweight("bold")
                    text.set_fontsize(9)

            # Color the branches according to cluster
            for clade in tree.find_clades():
                if clade.is_terminal() and clade.name:
                    node_name = str(clade.name).strip()
                    if node_name in cluster_map:
                        cluster_val = cluster_map[node_name]
                        color = cluster_colors.get(cluster_val, "black")
                        for line in ax.get_lines():
                            if line.get_label() == node_name:
                                line.set_color(color)
                                line.set_linewidth(2.5)
                                line.set_alpha(0.8)

            legend_elements = [
                plt.Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="w",
                    markerfacecolor=cluster_colors.get(c, "black"),
                    markersize=8,
                    label=f"Cluster {c}",
                )
                for c in unique_clusters
            ]
            ax.legend(
                handles=legend_elements,
                loc="upper right",
                title="Cluster Assignments",
                bbox_to_anchor=(1.01, 1),
                fontsize=10,
            )
            plt.title("Phylogenetic Tree with Colored Clusters", fontsize=14)
            plt.tight_layout()

            buffer = BytesIO()
            plt.savefig(buffer, format="png", bbox_inches="tight", dpi=120)
            buffer.seek(0)
            plt.close()

            image_data = base64.b64encode(buffer.getvalue()).decode("utf-8")
            html = f"""
            <div class="text-center">
                <img src="data:image/png;base64,{image_data}" class="img-fluid" alt="Phylogenetic Tree" style="max-width:100%; height:auto;">
            </div>
            <p class="text-center mt-2">
                Phylogenetic tree showing evolutionary relationships among {leaf_count} strains grouped into {len(unique_clusters)} distinct clusters.
                Branches and strain labels are colored according to their cluster assignment.
            </p>
            """
            return html
        except Exception as e:
            import traceback

            trace = traceback.format_exc()
            print(f"Error creating tree plot: {e}")
            return f"<p>Error creating tree visualization: {e}</p><pre>{trace}</pre>"

    def _scan_all_results(self):
        results_dict = {}
        csv_files = [f for f in os.listdir(self.output_folder) if f.endswith(".csv")]
        for csv_file in csv_files:
            try:
                df = pd.read_csv(os.path.join(self.output_folder, csv_file))
                section_name = csv_file.replace(".csv", "").replace(".", "_")
                results_dict[section_name] = {
                    "type": "table",
                    "data": df,
                    "title": section_name.replace("_", " ").title(),
                }
            except Exception as e:
                print(f"Could not process {csv_file}: {e}")

        image_files = [f for f in os.listdir(self.output_folder) if f.endswith((".png", ".jpg"))]
        for img_file in image_files:
            try:
                file_path = os.path.join(self.output_folder, img_file)
                with open(file_path, "rb") as ff:
                    img_data = base64.b64encode(ff.read()).decode("utf-8")
                section_name = img_file.replace(".png", "").replace(".jpg", "").replace(".", "_")
                results_dict[section_name] = {
                    "type": "image",
                    "data": img_data,
                    "title": section_name.replace("_", " ").title(),
                }
            except Exception as e:
                print(f"Could not process {img_file}: {e}")

        return results_dict

    def _generate_download_section(self):
        html = "<h2>All Raw Data Files</h2>"
        html += "<table class='table table-striped'>"
        html += "<thead><tr><th>File</th><th>Size (KB)</th><th>Action</th></tr></thead><tbody>"
        all_files = sorted(os.listdir(self.output_folder))
        for file in all_files:
            file_path = os.path.join(self.output_folder, file)
            if os.path.isfile(file_path):
                size_kb = os.path.getsize(file_path) / 1024
                if size_kb < 5000 and file.lower().endswith((".csv", ".txt")):
                    with open(file_path, "rb") as f:
                        data = base64.b64encode(f.read()).decode("utf-8")
                    download_link = f"<a href='data:text/csv;base64,{data}' download='{file}' class='btn btn-sm btn-primary'>Download</a>"
                else:
                    download_link = (
                        "<span class='btn btn-sm btn-secondary disabled'>No direct download</span>"
                    )
                html += f"<tr><td>{file}</td><td>{size_kb:.2f}</td><td>{download_link}</td></tr>"
        html += "</tbody></table>"
        return html

    def generate_report(self, analysis_results, config, output_file="phylogenetic_report.html"):
        template = self.template_env.get_template("report_template.html")
        results = {}

        tree_file = os.path.join(self.base_dir, config.tree_file)
        if os.path.exists(tree_file):
            results["tree_plot"] = self._create_interactive_tree_plot(tree_file)
        else:
            results["tree_plot"] = "<p>Tree file not found.</p>"

        evo_csv = os.path.join(self.output_folder, "evolutionary_cluster_analysis.csv")
        if os.path.exists(evo_csv):
            df_evo = pd.read_csv(evo_csv)
            results["tree_stats"] = self._df_to_interactive_table(
                df_evo, table_id="table-tree-stats", filename="tree_stats"
            )
        else:
            results["tree_stats"] = "<p>No evolutionary data found.</p>"

        if "summary_stats" in analysis_results and isinstance(
            analysis_results["summary_stats"], pd.DataFrame
        ):
            results["summary_stats"] = self._df_to_interactive_table(
                analysis_results["summary_stats"],
                table_id="table-summary-stats",
                filename="summary_stats",
            )
        else:
            results["summary_stats"] = "<p>No summary stats.</p>"

        if "embeddings" in analysis_results and "labels" in analysis_results:
            mask = analysis_results.get("mask", None)
            if mask is not None:
                full_labels = np.zeros(len(mask), dtype=int)
                full_labels[mask] = analysis_results["labels"]
                full_labels[~mask] = -1
            else:
                full_labels = analysis_results["labels"]

            df_umap = pd.DataFrame(
                {
                    "UMAP1": analysis_results["embeddings"][:, 0],
                    "UMAP2": analysis_results["embeddings"][:, 1],
                    "Cluster": full_labels,
                }
            )
            fig_umap = px.scatter(
                df_umap,
                x="UMAP1",
                y="UMAP2",
                color="Cluster",
                title="UMAP Visualization (Interactive)",
                labels={"UMAP1": "UMAP Dim1", "UMAP2": "UMAP Dim2"},
            )
            results["umap_plot"] = fig_umap.to_html(full_html=False, include_plotlyjs=False)
        else:
            results["umap_plot"] = "<p>No UMAP data integrated yet.</p>"

        dist_csv = os.path.join(self.output_folder, "cluster_distribution.csv")
        if os.path.exists(dist_csv):
            dist_df = pd.read_csv(dist_csv)
            if "Cluster" not in dist_df.columns:
                dist_df.rename(columns={dist_df.columns[0]: "Cluster"}, inplace=True)
            fig_dist = px.bar(
                dist_df,
                x="Cluster",
                y="Percentage",
                title="Cluster Distribution (%) - Interactive",
                labels={"Cluster": "Cluster", "Percentage": "% of Strains"},
                text="Strain_Count",
            )
            fig_dist.update_traces(textposition="outside")
            results["cluster_distribution"] = fig_dist.to_html(
                full_html=False, include_plotlyjs=False
            )
        else:
            results["cluster_distribution"] = "<p>No distribution data integrated yet.</p>"

        if os.path.exists(evo_csv):
            results["evolutionary_metrics"] = self._df_to_interactive_table(
                pd.read_csv(evo_csv), table_id="table-evo-metrics", filename="evolutionary_metrics"
            )
        else:
            results["evolutionary_metrics"] = "<p>No evolutionary metrics integrated yet.</p>"

        beta_csv = os.path.join(self.output_folder, "phylogenetic_beta_diversity.csv")
        if os.path.exists(beta_csv):
            df_beta = pd.read_csv(beta_csv, index_col=0)
            fig_beta = go.Figure(
                data=go.Heatmap(
                    z=df_beta.values,
                    x=df_beta.columns,
                    y=df_beta.index,
                    colorscale="Viridis",
                    colorbar=dict(title="Distance"),
                )
            )
            fig_beta.update_layout(
                title="Beta Diversity (Interactive)",
                xaxis=dict(title="Cluster"),
                yaxis=dict(title="Cluster"),
                height=500,
            )
            results["beta_diversity"] = fig_beta.to_html(full_html=False, include_plotlyjs=False)
        else:
            results["beta_diversity"] = "<p>No beta diversity data.</p>"

        evo_rates_csv = os.path.join(self.output_folder, "evolution_rates.csv")
        if os.path.exists(evo_rates_csv):
            df_rates = pd.read_csv(evo_rates_csv)
            fig_rates = px.bar(
                df_rates,
                x="Cluster_ID",
                y="EvolutionRate",
                title="Evolution Rates by Cluster",
                labels={"Cluster_ID": "Cluster", "EvolutionRate": "Rate"},
            )
            results["evolution_rates"] = fig_rates.to_html(full_html=False, include_plotlyjs=False)
        else:
            results["evolution_rates"] = "<p>No evolution rates data.</p>"

        phylo_signal_csv = os.path.join(self.output_folder, "phylogenetic_signal.csv")
        if os.path.exists(phylo_signal_csv):
            df_signal = pd.read_csv(phylo_signal_csv)
            results["phylogenetic_signal"] = self._df_to_interactive_table(
                df_signal, table_id="table-phylo-signal", filename="phylogenetic_signal"
            )
        else:
            results["phylogenetic_signal"] = "<p>No phylogenetic signal data.</p>"

        if "best_silhouette" in analysis_results:
            df_val = pd.DataFrame(
                {"SilhouetteScore": [round(analysis_results["best_silhouette"], 2)]}
            )
            df_val.to_csv(os.path.join(self.output_folder, "cluster_validation.csv"), index=False)
            results["cluster_validation"] = self._df_to_interactive_table(
                df_val, table_id="table-cluster-validation", filename="cluster_validation"
            )
        else:
            results["cluster_validation"] = "<p>No cluster validation data.</p>"

        trait_summary_path = os.path.join(self.output_folder, "trait_analysis_summary.html")
        if os.path.exists(trait_summary_path):
            with open(trait_summary_path, "r") as f:
                results["trait_analysis"] = f.read()
        else:
            results["trait_analysis"] = "<p>No trait analysis summary available.</p>"

        mca_row_csv = os.path.join(self.output_folder, "mca_row_coordinates.csv")
        mca_col_csv = os.path.join(self.output_folder, "mca_column_coordinates.csv")
        mca_summary_csv = os.path.join(self.output_folder, "mca_summary.csv")

        if os.path.exists(mca_row_csv) and os.path.exists(mca_summary_csv):
            df_mca_row = pd.read_csv(mca_row_csv)
            df_mca_summary = pd.read_csv(mca_summary_csv)
            fig_mca_row = px.scatter(
                df_mca_row,
                x="Component_1",
                y="Component_2",
                color="Cluster",
                title="MCA Row Points (Strains)",
                labels={
                    "Component_1": f'Component 1 ({round(df_mca_summary["Explained_Inertia"][0]*100, 1)}%)',
                    "Component_2": f'Component 2 ({round(df_mca_summary["Explained_Inertia"][1]*100, 1)}%)',
                },
                hover_data=["Cluster"],
            )
            if "Strain_ID" in df_mca_row.columns:
                fig_mca_row.update_traces(
                    hovertemplate="<b>Strain ID:</b> %{customdata[0]}<br><b>Cluster:</b> %{marker.color}<br>"
                    "<b>Component 1:</b> %{x:.3f}<br><b>Component 2:</b> %{y:.3f}",
                    customdata=df_mca_row[["Strain_ID"]],
                )
            fig_mca_row.update_layout(
                height=600,
                width=800,
                legend_title="Cluster",
                hoverlabel=dict(bgcolor="white", font_size=12),
            )
            results["mca_row_plot"] = fig_mca_row.to_html(full_html=False, include_plotlyjs=False)
            results[
                "mca_row_plot"
            ] += """
            <div class="mt-3">
                <h5>How to Read This Plot:</h5>
                <p>Each point represents a bacterial strain. Strains clustering together share similar genetic profiles.
                Hover over points to see the Strain ID and cluster assignment. Colors indicate cluster membership.</p>
            </div>
            """
        else:
            results["mca_row_plot"] = "<p>No MCA row data available.</p>"

        if os.path.exists(mca_col_csv) and os.path.exists(mca_summary_csv):
            df_mca_col = pd.read_csv(mca_col_csv)
            df_mca_col = df_mca_col.reset_index()
            df_mca_summary = pd.read_csv(mca_summary_csv)
            fig_mca_col = px.scatter(
                df_mca_col,
                x="Component_1",
                y="Component_2",
                color="Feature_Type",
                title="MCA Column Points (Features)",
                labels={
                    "Component_1": f'Component 1 ({round(df_mca_summary["Explained_Inertia"][0]*100, 1)}%)',
                    "Component_2": f'Component 2 ({round(df_mca_summary["Explained_Inertia"][1]*100, 1)}%)',
                    "Feature_Type": "Feature Category",
                },
                hover_name="index",
                hover_data=["index", "Feature_Type"],
                text="index",
            )
            fig_mca_col.update_traces(
                textposition="top center",
                textfont=dict(size=8),
                hovertemplate="<b>Feature:</b> %{hovertext}<br><b>Type:</b> %{customdata[1]}<br>"
                "<b>Component 1:</b> %{x:.3f}<br><b>Component 2:</b> %{y:.3f}",
            )
            fig_mca_col.update_layout(
                height=600,
                width=800,
                legend_title="Feature Category",
                hoverlabel=dict(bgcolor="white", font_size=12),
            )
            mca_guide_path = os.path.join(self.output_folder, "mca_interpretation_guide.csv")
            if os.path.exists(mca_guide_path):
                guide_df = pd.read_csv(mca_guide_path)
                guide_html = self._df_to_interactive_table(
                    guide_df, table_id="table-mca-guide", filename="mca_interpretation_guide"
                )
            else:
                guide_html = ""
            results[
                "mca_column_plot"
            ] = f"""
            {fig_mca_col.to_html(full_html=False, include_plotlyjs=False)}
            <div class="mt-3">
                <h5>How to Read This Plot:</h5>
                <p>Each point represents a genetic feature (e.g., MIC, AMR gene, or Virulence factor). 
                Features close together tend to co-occur in the same strains. Hover over points for detailed feature information.</p>
                <h5>Feature Type Explanation:</h5>
                <ul>
                    <li><strong>MIC</strong>: Minimum Inhibitory Concentration</li>
                    <li><strong>AMR</strong>: Antimicrobial Resistance Gene</li>
                    <li><strong>Virulence</strong>: Virulence Factor</li>
                </ul>
                {guide_html}
            </div>
            """
        else:
            results["mca_column_plot"] = "<p>No MCA column data available.</p>"

        if os.path.exists(mca_summary_csv):
            df_mca_sum = pd.read_csv(mca_summary_csv)
            results["mca_summary"] = self._df_to_interactive_table(
                df_mca_sum, table_id="table-mca-summary", filename="mca_summary"
            )
        else:
            results["mca_summary"] = "<p>No MCA summary data.</p>"

        log_odds_global_csv = os.path.join(self.output_folder, "log_odds_global.csv")
        log_odds_cluster_csv = os.path.join(self.output_folder, "log_odds_per_cluster.csv")

        if os.path.exists(log_odds_global_csv):
            df_log_global = pd.read_csv(log_odds_global_csv)
            results["global_log_odds"] = self._df_to_interactive_table(
                df_log_global, table_id="table-log-global", filename="log_odds_global"
            )
        else:
            results["global_log_odds"] = "<p>No global log-odds data.</p>"

        if os.path.exists(log_odds_cluster_csv):
            df_log_cluster = pd.read_csv(log_odds_cluster_csv)
            results["cluster_log_odds"] = self._df_to_interactive_table(
                df_log_cluster, table_id="table-log-cluster", filename="log_odds_cluster"
            )
        else:
            results["cluster_log_odds"] = "<p>No cluster log-odds data.</p>"

        assoc_csv = os.path.join(self.output_folder, "association_rules.csv")
        if os.path.exists(assoc_csv):
            df_assoc = pd.read_csv(assoc_csv)
            results["association_rules"] = self._df_to_interactive_table(
                df_assoc, table_id="table-assoc", filename="association_rules"
            )
        else:
            results["association_rules"] = "<p>No association rules data.</p>"

        shared_csv = os.path.join(self.output_folder, "shared_features.csv")
        unique_csv = os.path.join(self.output_folder, "unique_features.csv")

        if os.path.exists(shared_csv):
            df_shared = pd.read_csv(shared_csv)
            results["shared_features"] = self._df_to_interactive_table(
                df_shared, table_id="table-shared", filename="shared_features"
            )
        else:
            results["shared_features"] = "<p>No shared features data.</p>"

        if os.path.exists(unique_csv):
            df_unique = pd.read_csv(unique_csv)
            results["unique_features"] = self._df_to_interactive_table(
                df_unique, table_id="table-unique", filename="unique_features"
            )
        else:
            results["unique_features"] = "<p>No unique features data.</p>"

        fdr_csv = os.path.join(self.output_folder, "pairwise_fdr_post_hoc.csv")
        if os.path.exists(fdr_csv):
            df_fdr = pd.read_csv(fdr_csv)
            results["pairwise_fdr"] = self._df_to_interactive_table(
                df_fdr, table_id="table-fdr", filename="pairwise_fdr"
            )
        else:
            results["pairwise_fdr"] = "<p>No pairwise FDR data.</p>"

        boot_feat_csv = os.path.join(self.output_folder, "bootstrap_feature_importance.csv")
        boot_log_csv = os.path.join(self.output_folder, "bootstrap_log_odds.csv")

        if os.path.exists(boot_feat_csv):
            df_boot_feat = pd.read_csv(boot_feat_csv)
            results["bootstrap_feature_importance"] = self._df_to_interactive_table(
                df_boot_feat, table_id="table-boot-feat", filename="bootstrap_feature_importance"
            )
        else:
            results["bootstrap_feature_importance"] = "<p>No bootstrap feature importance data.</p>"

        if os.path.exists(boot_log_csv):
            df_boot_log = pd.read_csv(boot_log_csv)
            results["bootstrap_log_odds"] = self._df_to_interactive_table(
                df_boot_log, table_id="table-boot-log", filename="bootstrap_log_odds"
            )
        else:
            results["bootstrap_log_odds"] = "<p>No bootstrap log-odds data.</p>"

        results["download_section"] = self._generate_download_section()

        scanned = self._scan_all_results()
        additional_html = "<h2>Additional Analysis Results</h2>"
        for section_name, info in scanned.items():
            if section_name not in results:
                additional_html += f"<h3>{info['title']}</h3>"
                if info["type"] == "table":
                    additional_html += self._df_to_interactive_table(
                        info["data"], table_id=f"table-{section_name}", filename=section_name
                    )
                elif info["type"] == "image":
                    additional_html += f"<div class='mb-3'><img src='data:image/png;base64,{info['data']}' class='img-fluid' /></div>"

        results["additional_results"] = additional_html

        context = {
            "title": "StrepSuis-PhyloTrait: Integrated Phylogenetic and Binary Trait Analysis of Antimicrobial Resistance and Virulence in Streptococcus suis",
            "date": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
            "config": vars(config),
            "results": results,
        }

        html_rendered = template.render(**context)
        outpath = os.path.join(self.output_folder, output_file)
        with open(outpath, "w", encoding="utf-8") as f:
            f.write(html_rendered)
        print(f"HTML report saved to: {outpath}")
        return outpath

    def generate_excel_report(self, analysis_results, config):
        """
        Generate comprehensive Excel report with all analysis results and PNG charts.

        This function creates a detailed Excel workbook with multiple sheets containing:
        - Metadata and methodology
        - Phylogenetic clustering results
        - Evolutionary metrics (PD, pairwise distances, beta diversity)
        - Binary trait analysis (chi-square, log-odds, feature importance)
        - Association rules and MCA results
        - All CSV files generated during analysis

        Parameters:
            analysis_results (dict): Dictionary containing all analysis results
            config (Config): Configuration object with analysis parameters

        Returns:
            str: Path to generated Excel file
        """
        # Initialize Excel report generator
        excel_gen = ExcelReportGenerator(output_folder=self.output_folder)

        # Save matplotlib figures as PNG
        # Tree plot
        tree_plot_info = analysis_results.get("tree_plot", {})
        if tree_plot_info and "data" in tree_plot_info:
            try:
                import base64
                from io import BytesIO

                from PIL import Image

                img_data = tree_plot_info["data"]
                img = Image.open(BytesIO(base64.b64decode(img_data)))
                filepath = os.path.join(excel_gen.png_folder, "phylogenetic_tree.png")
                img.save(filepath)
                excel_gen.png_files.append(filepath)
                print(f"Saved phylogenetic tree: {filepath}")
            except Exception as e:
                print(f"Could not save tree plot: {e}")

        # UMAP plot
        umap_plot_info = analysis_results.get("umap_plot", {})
        if umap_plot_info and "data" in umap_plot_info:
            try:
                import base64
                from io import BytesIO

                from PIL import Image

                img_data = umap_plot_info["data"]
                img = Image.open(BytesIO(base64.b64decode(img_data)))
                filepath = os.path.join(excel_gen.png_folder, "umap_embedding.png")
                img.save(filepath)
                excel_gen.png_files.append(filepath)
                print(f"Saved UMAP plot: {filepath}")
            except Exception as e:
                print(f"Could not save UMAP plot: {e}")

        # Prepare methodology description
        methodology = {
            "Phylogenetic Clustering": (
                "Tree-aware clustering methods (TreeCluster, Phydelity-inspired) that respect tree structure. "
                "Pairwise patristic distances computed from phylogenetic tree. "
                "Optimal cluster number selected using silhouette analysis."
            ),
            "Outlier Detection": (
                "Isolation Forest algorithm to identify and remove outlier strains. "
                "Contamination parameter controls the expected proportion of outliers."
            ),
            "Ensemble Clustering": (
                "Multiple clustering algorithms tested (KMeans, DBSCAN, Gaussian Mixture). "
                "Best performer selected based on silhouette coefficient. "
                "Optuna used for hyperparameter optimization."
            ),
            "Evolutionary Metrics": (
                "Faith's Phylogenetic Diversity (PD) quantifies evolutionary history captured by clusters. "
                "Pairwise patristic distances measure evolutionary divergence. "
                "Beta diversity metrics assess cluster differentiation."
            ),
            "Binary Trait Analysis": (
                "Chi-square tests for trait-cluster associations with Benjamini-Hochberg FDR correction. "
                "Log-odds ratios with bootstrap confidence intervals for effect size estimation. "
                "Random Forest feature importance for identifying discriminative traits."
            ),
            "Association Rules": (
                "mlxtend library for mining frequent patterns and association rules. "
                "Apriori algorithm identifies trait combinations with high support and confidence."
            ),
            "Multiple Correspondence Analysis": (
                "prince library for dimensionality reduction of categorical trait data. "
                "Produces low-dimensional representation for visualization."
            ),
            "Statistical Corrections": (
                "Benjamini-Hochberg FDR correction for multiple hypothesis testing. "
                "Pairwise comparisons with FDR adjustment for cluster-specific tests."
            ),
        }

        # Prepare sheets data
        sheets_data = {}

        # Scan output folder for CSV files
        csv_files = {}
        for f in os.listdir(self.output_folder):
            if f.lower().endswith(".csv"):
                try:
                    df = pd.read_csv(os.path.join(self.output_folder, f))
                    csv_files[f] = df
                except Exception as e:
                    print(f"Could not load {f}: {e}")

        # Add CSV files as sheets with descriptive names
        for csv_name, df in csv_files.items():
            sheet_name = sanitize_sheet_name(csv_name.replace(".csv", ""))
            description = f"Data from {csv_name}"

            # Add more specific descriptions based on filename
            if "cluster" in csv_name.lower():
                description = "Phylogenetic cluster assignments for all strains"
            elif "outlier" in csv_name.lower():
                description = "Outlier detection results"
            elif "chi" in csv_name.lower() or "chi2" in csv_name.lower():
                description = "Chi-square test results for trait-cluster associations"
            elif "odds" in csv_name.lower() or "logodds" in csv_name.lower():
                description = "Log-odds ratios for trait enrichment in clusters"
            elif "importance" in csv_name.lower() or "feat" in csv_name.lower():
                description = "Random Forest feature importance scores"
            elif "rules" in csv_name.lower() or "assoc" in csv_name.lower():
                description = "Association rules for trait patterns"
            elif "mca" in csv_name.lower():
                description = "Multiple Correspondence Analysis coordinates"
            elif "pd" in csv_name.lower() or "diversity" in csv_name.lower():
                description = "Phylogenetic diversity and evolutionary metrics"
            elif "pairwise" in csv_name.lower() or "distance" in csv_name.lower():
                description = "Pairwise patristic distances between strains"
            elif "trait" in csv_name.lower():
                description = "Binary trait analysis results"

            sheets_data[sheet_name] = (df, description)

        # Prepare metadata
        metadata = {
            "Base_Directory": config.base_dir,
            "Output_Folder": config.output_folder,
            "Tree_File": config.tree_file,
            "MIC_File": config.mic_file,
            "AMR_Genes_File": config.amr_genes_file,
            "Virulence_Genes_File": config.virulence_genes_file,
            "UMAP_Components": config.umap_components,
            "UMAP_Neighbors": config.umap_neighbors,
            "Outlier_Contamination": config.outlier_contamination,
            "Bootstrap_Iterations": config.bootstrap_iterations,
            "FDR_Alpha": config.fdr_alpha,
            "Total_CSV_Files": len(csv_files),
            "Total_PNG_Charts": len(excel_gen.png_files),
        }

        # Generate Excel report
        excel_path = excel_gen.generate_excel_report(
            report_name="Phylogenetic_Clustering_Report",
            sheets_data=sheets_data,
            methodology=methodology,
            **metadata,
        )

        print(f"Excel report saved to: {excel_path}")
        return excel_path

    def convert_html_to_pdf(self, html_path, output_file=None):
        if output_file is None:
            output_file = html_path.replace(".html", ".pdf")
        weasyprint.HTML(filename=html_path).write_pdf(output_file)
        print(f"PDF report saved to: {output_file}")
        return output_file


###############################################################################
# 12. CONFIGURATION
###############################################################################
class Config:
    """Configuration for phylogenetic analysis."""

    def __init__(self, **kwargs):
        self.base_dir = kwargs.get("base_dir", ".")
        self.output_folder = kwargs.get(
            "output_folder", "phylogenetic_clustering_results_advanced17xY"
        )
        self.tree_file = kwargs.get("tree_file", "tree.newick.txt")
        self.mic_file = kwargs.get("mic_file", "MIC.csv")
        self.amr_genes_file = kwargs.get("amr_genes_file", "AMR_genes.csv")
        self.virulence_genes_file = kwargs.get("virulence_genes_file", "Virulence.csv")
        self.umap_components = kwargs.get("umap_components", 2)
        self.umap_neighbors = kwargs.get("umap_neighbors", 15)
        self.umap_min_dist = kwargs.get("umap_min_dist", 0.1)
        self.outlier_contamination = kwargs.get("outlier_contamination", 0.05)
        self.outlier_n_estimators = kwargs.get("outlier_n_estimators", 200)
        self.n_clusters_range = kwargs.get("n_clusters_range", (2, 10))
        self.n_ensemble = kwargs.get("n_ensemble", 5)
        self.dbscan_trials = kwargs.get("dbscan_trials", 20)
        self.bootstrap_iterations = kwargs.get("bootstrap_iterations", 500)
        self.fdr_alpha = kwargs.get("fdr_alpha", 0.05)
        self.parallel_tree = kwargs.get("parallel_tree", False)
        self.parallel_jobs = kwargs.get("parallel_jobs", 1)


###############################################################################
# 13. MAIN PIPELINE
###############################################################################
class PhylogeneticAnalysis:
    """
    Main pipeline orchestrating phylogenetic analysis workflow.

    This class coordinates tree-aware clustering, evolutionary analysis,
    trait analysis, and report generation for bacterial genomics data.

    Parameters
    ----------
    config : Config
        Configuration object containing analysis parameters and file paths.

    Attributes
    ----------
    config : Config
        Analysis configuration
    output_folder : str
        Path to output directory
    core : PhylogeneticCore
        Phylogenetic tree operations handler
    clustering : ClusteringModule
        Clustering algorithm manager
    data_loader : DataLoader
        Data loading and merging handler
    visualizer : Visualizer
        Visualization generation handler
    trait_analyzer : TraitAnalyzer
        Binary trait analysis handler
    mca_analyzer : MCAAnalyzer
        Multiple Correspondence Analysis handler

    Examples
    --------
    >>> config = Config(base_dir=".", tree_file="tree.newick")
    >>> analysis = PhylogeneticAnalysis(config)
    >>> results = analysis.run_complete_analysis()
    """

    def __init__(self, config):
        """
        Initialize phylogenetic analysis pipeline.

        Parameters
        ----------
        config : Config
            Configuration object with analysis parameters.
        """
        self.config = config
        self.base_dir = config.base_dir
        self.output_folder = os.path.join(self.base_dir, config.output_folder)
        os.makedirs(self.output_folder, exist_ok=True)

        # Setup logging
        self.log_file = setup_logging(self.output_folder)
        logging.info("Phylogenetic Analysis Pipeline Initialized")
        logging.info(f"Output folder: {self.output_folder}")
        logging.info(f"Log file: {self.log_file}")

        # Initialize analysis components
        self.core = PhylogeneticCore()
        self.clustering = ClusteringModule(
            n_clusters_range=config.n_clusters_range,
            n_ensemble=config.n_ensemble,
            dbscan_trials=config.dbscan_trials,
            seed=42,
        )
        self.data_loader = DataLoader(self.base_dir)
        self.visualizer = Visualizer(self.output_folder)
        self.trait_analyzer = TraitAnalyzer(self.output_folder)
        self.mca_analyzer = MCAAnalyzer(self.output_folder)

        # Execution metrics
        self.start_time = None
        self.step_times = {}

    def determine_adaptive_cluster_range(self, embeddings, distance_matrix=None):
        """
        Adaptively determines the appropriate cluster range based on dataset characteristics.
        """
        n_samples = embeddings.shape[0]
        # 1. Base min and max on dataset size
        theoretical_max = int(np.sqrt(n_samples / 2))

        # 2. Analyze distance distribution for structure insights
        if distance_matrix is not None:
            flat_distances = distance_matrix[np.triu_indices(n_samples, k=1)]
            q1, median, q3 = np.percentile(flat_distances, [25, 50, 75])
            iqr = q3 - q1
            structure_ratio = iqr / median if median > 0 else 1
            if structure_ratio > 1.5:
                max_k = min(theoretical_max, 20)
            else:
                max_k = min(theoretical_max, 10)
        else:
            max_k = min(theoretical_max, 15)

        # 3. Density-based estimation using nearest neighbors
        from sklearn.neighbors import NearestNeighbors

        nn = NearestNeighbors(n_neighbors=min(20, n_samples - 1))
        nn.fit(embeddings)
        distances, _ = nn.kneighbors(embeddings)
        avg_nn_dist = np.mean(distances[:, 1:], axis=1)
        density_peaks = 0
        for i in range(len(avg_nn_dist)):
            neighbors = nn.kneighbors(
                embeddings[i].reshape(1, -1), n_neighbors=min(10, n_samples), return_distance=False
            )[0]
            if np.all(avg_nn_dist[i] < avg_nn_dist[neighbors[1:]]):
                density_peaks += 1

        min_k = max(2, min(density_peaks, 5))
        max_k = max(min_k + 2, max_k)
        print(
            f"Adaptively determined cluster range: {min_k}-{max_k} (based on {n_samples} samples, {density_peaks} density peaks)"
        )
        return (min_k, max_k)

    def determine_optimal_clusters(self, embeddings, cluster_range=(2, 10)):
        """
        Determine the optimal number of clusters using KMeans and silhouette analysis.
        """
        best_score = -1
        best_n = None
        best_labels = None
        for n in range(cluster_range[0], cluster_range[1] + 1):
            kmeans = KMeans(n_clusters=n, random_state=42)
            labels = kmeans.fit_predict(embeddings)
            score = silhouette_score(embeddings, labels)
            print(f"KMeans with {n} clusters: Silhouette score = {score:.4f}")
            if score > best_score:
                best_score = score
                best_n = n
                best_labels = labels
        return best_n, best_labels, best_score

    def test_multiple_clustering_methods(self, tree_clustering, distance_matrix):
        """
        Test different clustering methods (tree-aware and traditional) and return the best labels.
        Includes an improved approach for automatically determining optimal clusters.
        """
        best_labels = None
        best_silhouette = -1
        best_method = None
        best_params = {}

        # Test tree-aware methods first
        for method in ["max", "sum", "avg"]:
            print(f"\nTesting TreeCluster with method: {method}")
            for threshold_factor in [0.8, 1.0, 1.2, 1.5]:
                if method == "max":
                    threshold = (
                        tree_clustering._auto_threshold_max(
                            conservative_factor=5.0, distance_matrix=distance_matrix
                        )
                        * threshold_factor
                    )
                elif method == "sum":
                    threshold = (
                        tree_clustering._auto_threshold_sum(
                            conservative_factor=5.0, distance_matrix=distance_matrix
                        )
                        * threshold_factor
                    )
                else:
                    threshold = (
                        tree_clustering._auto_threshold_avg(
                            conservative_factor=5.0, distance_matrix=distance_matrix
                        )
                        * threshold_factor
                    )

                labels = tree_clustering.tree_cluster_algorithm(
                    distance_matrix, method=method, threshold=threshold
                )
                if len(np.unique(labels)) > 1:
                    try:
                        sil = silhouette_score(distance_matrix, labels, metric="precomputed")
                        print(
                            f"Method: {method}, Threshold: {threshold:.4f}, Silhouette: {sil:.4f}, Clusters: {len(np.unique(labels))}"
                        )
                        if sil > best_silhouette:
                            best_silhouette = sil
                            best_labels = labels
                            best_method = f"TreeCluster-{method}"
                            best_params = {
                                "threshold": threshold,
                                "threshold_factor": threshold_factor,
                            }
                    except Exception as e:
                        print(f"Error: {e}")

        # Test traditional clustering methods on UMAP embeddings
        print("\nTesting traditional clustering methods")
        embeddings = self.core.dimension_reduction(
            distance_matrix,
            n_components=self.config.umap_components,
            n_neighbors=self.config.umap_neighbors,
            min_dist=self.config.umap_min_dist,
            random_state=42,
        )
        for n_clusters in range(2, 8):
            try:
                kmeans = KMeans(n_clusters=n_clusters, random_state=42)
                labels = kmeans.fit_predict(embeddings)
                sil = silhouette_score(distance_matrix, labels, metric="precomputed")
                print(f"KMeans with {n_clusters} clusters: Silhouette: {sil:.4f}")
                if sil > best_silhouette:
                    best_silhouette = sil
                    best_labels = labels
                    best_method = "KMeans"
            except Exception as e:
                print(f"Error in KMeans: {e}")

        try:
            dbscan = DBSCAN(eps=0.5, min_samples=5)
            labels = dbscan.fit_predict(embeddings)
            if len(set(labels)) > 1:
                sil = silhouette_score(distance_matrix, labels, metric="precomputed")
                print(f"DBSCAN: Silhouette: {sil:.4f}")
                if sil > best_silhouette:
                    best_silhouette = sil
                    best_labels = labels
                    best_method = "DBSCAN"
            else:
                print("DBSCAN produced only one cluster.")
        except Exception as e:
            print(f"Error in DBSCAN: {e}")

        for n_clusters in range(2, 8):
            try:
                gm = GaussianMixture(n_components=n_clusters, random_state=42)
                labels = gm.fit_predict(embeddings)
                if len(set(labels)) > 1:
                    sil = silhouette_score(distance_matrix, labels, metric="precomputed")
                    print(f"GaussianMixture with {n_clusters} clusters: Silhouette: {sil:.4f}")
                    if sil > best_silhouette:
                        best_silhouette = sil
                        best_labels = labels
                        best_method = "GaussianMixture"
            except Exception as e:
                print(f"Error in GaussianMixture: {e}")

        print("\nDetermining optimal number of clusters using silhouette analysis with KMeans")
        optimal_n, optimal_labels, optimal_score = self.determine_optimal_clusters(
            embeddings, cluster_range=self.adaptive_cluster_range  # Using the adaptive range
        )
        print(f"Optimal number of clusters: {optimal_n} with silhouette score: {optimal_score:.4f}")
        if optimal_score > best_silhouette:
            best_silhouette = optimal_score
            best_labels = optimal_labels
            best_method = f"KMeans-optimal-{optimal_n}"

        print(f"\nBest method: {best_method} with silhouette {best_silhouette:.4f}")
        return best_labels, best_silhouette, best_method, best_params

    def run_complete_analysis(self):
        """
        Execute complete phylogenetic analysis pipeline.

        Runs all analysis steps including tree-aware clustering, evolutionary
        analysis, trait analysis, and MCA. Generates comprehensive results
        and visualizations.

        Returns
        -------
        dict or None
            Dictionary containing analysis results with keys:
            - 'embeddings': UMAP embeddings
            - 'labels': Cluster labels
            - 'mask': Inlier mask
            - 'best_silhouette': Best silhouette score
            - 'best_method': Best clustering method name
            - 'best_params': Parameters of best method
            Returns None if analysis fails.

        Notes
        -----
        Pipeline steps:
        1. Phylogenetic clustering using tree-aware methods
        2. Evolutionary analysis (diversity metrics, beta diversity)
        3. Binary trait analysis (log-odds, association rules)
        4. Multiple Correspondence Analysis (MCA)

        All results are saved to output folder specified in config.
        Progress is logged to console and log file.
        """
        self.start_time = time.time()

        print_section_header("PHYLOGENETIC ANALYSIS PIPELINE - STARTING")
        logging.info(f"Analysis Start Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        logging.info("Configuration:")
        logging.info(f"  - Tree file: {self.config.tree_file}")
        logging.info(f"  - MIC file: {self.config.mic_file}")
        logging.info(f"  - AMR genes file: {self.config.amr_genes_file}")
        logging.info(f"  - Virulence genes file: {self.config.virulence_genes_file}")
        logging.info(f"  - UMAP components: {self.config.umap_components}")
        logging.info(f"  - UMAP neighbors: {self.config.umap_neighbors}")
        logging.info(f"  - Outlier contamination: {self.config.outlier_contamination}")
        logging.info(f"  - Bootstrap iterations: {self.config.bootstrap_iterations}")
        logging.info(f"  - FDR alpha: {self.config.fdr_alpha}")
        print_memory_usage()

        analysis_results = {}
        total_steps = 4

        try:
            # Step 1: Phylogenetic Clustering using tree-aware methods
            step_start = time.time()
            print_step(1, total_steps, "Phylogenetic Clustering")

            tree_path = os.path.join(self.base_dir, self.config.tree_file)
            logging.info(f"Loading phylogenetic tree from: {tree_path}")
            tree = self.core.load_tree(tree_path)

            logging.info("Computing distance matrix from tree...")
            distance_matrix, terminals = self.core.tree_to_distance_matrix(
                tree, parallel=False, n_jobs=1
            )
            strain_names = [str(t).strip() for t in terminals]
            logging.info(f"Loaded {len(strain_names)} strains from tree")

            logging.info("Performing UMAP dimensionality reduction...")
            embeddings = self.core.dimension_reduction(
                distance_matrix,
                n_components=self.config.umap_components,
                n_neighbors=self.config.umap_neighbors,
                min_dist=self.config.umap_min_dist,
                random_state=42,
            )
            logging.info(f"UMAP embeddings shape: {embeddings.shape}")

            # Determine adaptive cluster range
            logging.info("Determining adaptive cluster range...")
            self.adaptive_cluster_range = self.determine_adaptive_cluster_range(
                embeddings, distance_matrix
            )
            logging.info(f"Adaptive cluster range: {self.adaptive_cluster_range}")

            # Create tree-aware clustering module
            logging.info("Initializing tree-aware clustering...")
            tree_clustering = TreeAwareClusteringModule(
                tree, terminals, n_clusters_range=self.adaptive_cluster_range, seed=42
            )

            logging.info("Testing multiple clustering methods...")
            best_labels, best_silhouette, best_method, best_params = (
                self.test_multiple_clustering_methods(tree_clustering, distance_matrix)
            )

            if best_labels is not None:
                logging.info(f"Best clustering method: {best_method}")
                logging.info(f"Best silhouette score: {best_silhouette:.4f}")

                pre_monophyly = tree_clustering.evaluate_monophyly(best_labels)
                logging.info(
                    f"Before enforcing monophyly: {pre_monophyly['monophyletic_percentage']:.2f}% clusters are monophyletic"
                )

                logging.info("Enforcing monophyletic clusters...")
                best_labels = tree_clustering.ensure_monophyletic_clusters(best_labels)
                post_monophyly = tree_clustering.evaluate_monophyly(best_labels)
                logging.info(
                    f"After enforcing monophyly: {post_monophyly['monophyletic_percentage']:.2f}% clusters are monophyletic"
                )

                # Offset final labels by +1 so clusters start at 1
                best_labels = best_labels + 1
                logging.info(
                    f"Final cluster labels range: {best_labels.min()} to {best_labels.max()}"
                )

            else:
                logging.warning(
                    "Tree-aware clustering failed. Falling back to standard ensemble clustering."
                )
                embeddings = self.core.dimension_reduction(
                    distance_matrix,
                    n_components=self.config.umap_components,
                    n_neighbors=self.config.umap_neighbors,
                    min_dist=self.config.umap_min_dist,
                    random_state=42,
                )
                clean_embeddings, mask = self.core.detect_outliers(
                    embeddings,
                    contamination=self.config.outlier_contamination,
                    n_estimators=self.config.outlier_n_estimators,
                    random_state=42,
                )
                labels, best_silhouette = self.clustering.ensemble_clustering(clean_embeddings)
                if labels is None:
                    logging.error("Clustering failed completely. Stopping analysis.")
                    return None
                outlier_assignments = self.clustering.assign_outliers_to_clusters(
                    embeddings, mask, labels
                )
                full_labels = np.zeros(len(strain_names), dtype=int)
                full_labels[mask] = labels
                for idx, cluster in outlier_assignments:
                    full_labels[idx] = cluster
                # Also offset final labels by +1
                full_labels = full_labels + 1
                best_labels = full_labels
                mask = np.ones(len(strain_names), dtype=bool)

            mask = np.ones(len(strain_names), dtype=bool)
            logging.info("Saving clustering results...")
            self._save_clustering_results(strain_names, mask, best_labels, [])

            # Plot clusters on UMAP and phylogenetic tree
            logging.info("Generating visualizations...")
            self.visualizer.plot_umap_clusters(embeddings, best_labels, mask, [])
            self.visualizer.plot_phylogenetic_tree(tree, best_labels, strain_names, mask)
            logging.info("Visualizations saved")

            analysis_results["embeddings"] = embeddings
            analysis_results["labels"] = best_labels
            analysis_results["mask"] = mask
            analysis_results["best_silhouette"] = best_silhouette
            analysis_results["best_method"] = best_method
            analysis_results["best_params"] = best_params

            self.step_times["Step 1: Phylogenetic Clustering"] = time.time() - step_start
            logging.info(
                f"Step 1 completed in {self.step_times['Step 1: Phylogenetic Clustering']:.2f} seconds"
            )

            # Step 2: Evolutionary Analysis
            step_start = time.time()
            print_step(2, total_steps, "Evolutionary Analysis")

            logging.info("Computing evolutionary cluster metrics...")
            cluster_df = EvolutionaryAnalysis.evolutionary_cluster_analysis(
                tree, best_labels, strain_names, mask
            )
            cluster_df.to_csv(
                os.path.join(self.output_folder, "evolutionary_cluster_analysis.csv"), index=False
            )
            logging.info(f"Evolutionary cluster analysis saved ({len(cluster_df)} clusters)")

            logging.info("Calculating phylogenetic beta diversity...")
            beta_div = EvolutionaryAnalysis.calculate_beta_diversity(
                tree, best_labels, strain_names, mask
            )
            beta_div.to_csv(os.path.join(self.output_folder, "phylogenetic_beta_diversity.csv"))
            logging.info("Beta diversity matrix saved")

            logging.info("Calculating evolution rates...")
            rates_df = EvolutionaryAnalysis.calculate_evolution_rates(cluster_df)
            rates_df.to_csv(os.path.join(self.output_folder, "evolution_rates.csv"), index=False)
            logging.info(f"Evolution rates saved ({len(rates_df)} clusters)")

            logging.info("Calculating phylogenetic signal...")
            EvolutionaryAnalysis.calculate_phylogenetic_signal(
                cluster_df, self.output_folder
            )
            logging.info("Phylogenetic signal analysis complete")

            self.step_times["Step 2: Evolutionary Analysis"] = time.time() - step_start
            logging.info(
                f"Step 2 completed in {self.step_times['Step 2: Evolutionary Analysis']:.2f} seconds"
            )

            # Step 3: Trait Analysis
            step_start = time.time()
            print_step(3, total_steps, "Binary Trait Analysis")

            clusters_csv = os.path.join(self.output_folder, "phylogenetic_clusters.csv")
            logging.info(f"Loading and merging trait data from: {clusters_csv}")
            merged_df = self.data_loader.load_and_merge_data(clusters_csv)

            if merged_df is not None:
                logging.info(f"Merged data shape: {merged_df.shape}")
                logging.info(f"Unique clusters: {merged_df['Cluster'].nunique()}")

                logging.info("Analyzing trait categories...")
                trait_summary_html = self.trait_analyzer.analyze_all_categories(merged_df)

                logging.info("Plotting cluster distribution...")
                analysis_results["summary_stats"] = self.visualizer.plot_cluster_distribution(
                    merged_df
                )

                logging.info("Performing log-odds ratio analysis...")
                log_odds_global, log_odds_cluster = self.trait_analyzer.log_odds_ratio_analysis(
                    merged_df
                )
                # Save global log-odds results to CSV
                log_odds_global_path = os.path.join(self.output_folder, "log_odds_global.csv")
                if hasattr(log_odds_global, "to_csv"):
                    log_odds_global.to_csv(log_odds_global_path, index=False)
                    logging.info(f"Global log-odds results saved to: {log_odds_global_path}")

                logging.info("Mining association rules...")
                self.trait_analyzer.association_rule_mining(merged_df)

                logging.info("Labeling shared and unique features...")
                self.trait_analyzer.label_shared_unique_features(merged_df, presence_threshold=0.3)

                logging.info("Performing pairwise FDR post-hoc tests...")
                self.trait_analyzer.pairwise_fdr_post_hoc(merged_df)

                logging.info("Computing bootstrap feature importance...")
                self.trait_analyzer.bootstrap_feature_importance(merged_df, n_bootstrap=100)

                logging.info("Computing bootstrap log-odds confidence intervals...")
                self.trait_analyzer.bootstrap_log_odds(merged_df, n_bootstrap=100)

                trait_summary_path = os.path.join(self.output_folder, "trait_analysis_summary.html")
                with open(trait_summary_path, "w") as f:
                    f.write(trait_summary_html)
                logging.info(f"Trait analysis summary saved to: {trait_summary_path}")
            else:
                logging.error("Merged data is None. Skipping trait analysis.")
                return None

            self.step_times["Step 3: Binary Trait Analysis"] = time.time() - step_start
            logging.info(
                f"Step 3 completed in {self.step_times['Step 3: Binary Trait Analysis']:.2f} seconds"
            )

            # Step 4: MCA Analysis
            step_start = time.time()
            print_step(4, total_steps, "Multiple Correspondence Analysis (MCA)")

            if merged_df is not None:
                logging.info("Performing MCA analysis...")
                self.mca_analyzer.perform_mca_analysis(merged_df)
                logging.info("MCA analysis complete")

            self.step_times["Step 4: MCA Analysis"] = time.time() - step_start
            logging.info(
                f"Step 4 completed in {self.step_times['Step 4: MCA Analysis']:.2f} seconds"
            )

            # Print execution summary
            self._print_execution_summary()

            print_section_header("ANALYSIS COMPLETED SUCCESSFULLY")
            return analysis_results

        except Exception as e:
            logging.error(f"Error in analysis pipeline: {e}")
            import traceback

            logging.error(traceback.format_exc())
            return None

    def _print_execution_summary(self):
        """Print comprehensive execution summary."""
        total_time = time.time() - self.start_time

        print_section_header("EXECUTION SUMMARY")
        logging.info(
            f"Total Execution Time: {total_time:.2f} seconds ({total_time/60:.2f} minutes)"
        )
        _safe_log("")
        logging.info("Step-by-Step Timing:")
        for step_name, step_time in self.step_times.items():
            percentage = (step_time / total_time) * 100
            logging.info(f"  {step_name}: {step_time:.2f}s ({percentage:.1f}%)")

        _safe_log("")
        logging.info("Output Files Generated:")
        csv_count = 0
        png_count = 0
        html_count = 0

        for file in os.listdir(self.output_folder):
            if file.endswith(".csv"):
                csv_count += 1
            elif file.endswith(".png"):
                png_count += 1
            elif file.endswith(".html"):
                html_count += 1

        logging.info(f"  - CSV files: {csv_count}")
        logging.info(f"  - PNG images: {png_count}")
        logging.info(f"  - HTML files: {html_count}")
        _safe_log("")
        logging.info(f"Output folder: {self.output_folder}")
        logging.info(f"Log file: {self.log_file}")
        print_memory_usage()

    def _save_clustering_results(self, strain_names, mask, labels, outlier_assignments):
        clean_strain_names = np.array(strain_names)[mask]
        df_inliers = pd.DataFrame({"Strain_ID": clean_strain_names, "Cluster": labels})

        outlier_data = []
        for idx, cluster in outlier_assignments:
            outlier_data.append({"Strain_ID": strain_names[idx], "Cluster": cluster})
        df_outliers = pd.DataFrame(outlier_data)

        df_all = pd.concat([df_inliers, df_outliers], ignore_index=True)
        df_all["Cluster"] = df_all["Cluster"].astype(int)
        df_all.to_csv(os.path.join(self.output_folder, "phylogenetic_clusters.csv"), index=False)


###############################################################################
# 14. MAIN EXECUTION BLOCK
###############################################################################
def main(
    base_dir: str = ".",
    output_folder: str = "output",
    tree_file: str = "Snp_tree.newick",
):
    """
    Main entry point for phylogenetic analysis.
    
    Args:
        base_dir: Base directory containing data files
        output_folder: Output directory for results
        tree_file: Name of the Newick tree file
    """
    print(
        f"""
{'='*80}
PHYLOGENETIC ANALYSIS SCRIPT
Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}
Author: MK-vet (with tree-aware clustering improvements)
{'='*80}
"""
    )

    # Configuration
    config = Config(
        base_dir=base_dir,
        output_folder=output_folder,
        tree_file=tree_file,
        mic_file="MIC.csv",
        amr_genes_file="AMR_genes.csv",
        virulence_genes_file="Virulence.csv",
        umap_components=2,
        umap_neighbors=15,
        umap_min_dist=0.1,
        outlier_contamination=0.05,
        outlier_n_estimators=200,
        n_clusters_range=(2, 10),
        n_ensemble=5,
        dbscan_trials=20,
        bootstrap_iterations=500,
        fdr_alpha=0.05,
        parallel_tree=False,
        parallel_jobs=1,
    )

    # Run analysis pipeline
    analysis = PhylogeneticAnalysis(config)
    final_results = analysis.run_complete_analysis()

    if final_results is not None:
        # Generate comprehensive reports
        print_section_header("GENERATING COMPREHENSIVE REPORTS")
        report_generator = HTMLReportGenerator(config.output_folder, base_dir=config.base_dir)

        # Generate HTML report
        logging.info("Generating interactive HTML report...")
        html_report_path = report_generator.generate_report(
            final_results, config, "phylogenetic_report.html"
        )
        logging.info(f"HTML report saved to: {html_report_path}")

        # Generate Excel report
        logging.info("Generating comprehensive Excel report...")
        excel_report_path = report_generator.generate_excel_report(final_results, config)
        logging.info(f"Excel report saved to: {excel_report_path}")

        # Print final summary
        print_section_header("ANALYSIS COMPLETE - SUMMARY")
        logging.info("All reports generated successfully!")
        _safe_log("")
        logging.info(f"Output Directory: {config.output_folder}")
        logging.info(f"  - HTML Report: {os.path.basename(html_report_path)}")
        logging.info(f"  - Excel Report: {os.path.basename(excel_report_path)}")
        logging.info("  - Log File: phylogenetic_analysis.log")
        _safe_log("")
        logging.info("Please review the reports for detailed analysis results.")
    else:
        logging.error("Analysis was interrupted due to errors.")
        logging.error("Please check the log file for details.")

    _safe_log("")
    logging.info(f"{'='*80}")
    logging.info(
        f"End of Phylogenetic Analysis: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}"
    )
    logging.info(f"{'='*80}")
    logging.shutdown()


if __name__ == "__main__":
    main()
