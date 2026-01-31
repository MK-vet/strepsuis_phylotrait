"""
Large Data Integration for PhyloTrait Module

This module integrates server-side large data pipeline components for
large phylogenetic tree visualization and analysis (5k+ taxa).

Features:
- NetworkX server-side tree layout computation
- Pre-rendered phylogenetic trees (no raw tree data to browser)
- DuckDB-based trait queries
- Support for trees with 10k+ taxa

Author: MK-vet
License: MIT
"""

import numpy as np
import pandas as pd
import networkx as nx
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, Any
import logging

# Import shared large data pipeline components
try:
    from shared.duckdb_handler import DuckDBHandler
    from shared.networkx_server_layout import NetworkXServerLayout
    from shared.datashader_plots import DatashaderPlots
    from shared.holoviews_heatmaps import HoloViewsHeatmaps
    LARGE_DATA_AVAILABLE = True
except ImportError:
    LARGE_DATA_AVAILABLE = False

logger = logging.getLogger(__name__)


class LargeDataPhyloTrait:
    """
    Large-scale phylogenetic tree processing and visualization.

    Handles large phylogenetic trees (5k+ taxa) using server-side processing.
    """

    def __init__(self, output_dir: Optional[Path] = None):
        """
        Initialize large data integration.

        Parameters
        ----------
        output_dir : Path, optional
            Directory for output files
        """
        if not LARGE_DATA_AVAILABLE:
            raise ImportError(
                "Large data pipeline not available. "
                "Install with: pip install duckdb matplotlib networkx datashader"
            )

        self.output_dir = output_dir or Path("phylotrait_large_data_output")
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.db_handler = DuckDBHandler()
        self.layout_renderer = NetworkXServerLayout(figsize=(18, 16), dpi=150)
        self.plotter = DatashaderPlots(width=1400, height=1000)
        self.heatmap = HoloViewsHeatmaps(tile_size=512)

        logger.info(f"Large data PhyloTrait initialized: {self.output_dir}")

    def load_trait_data(
        self,
        data_path: Union[str, Path],
        table_name: str = "trait_data"
    ) -> Dict[str, Any]:
        """
        Load large trait dataset into DuckDB.

        Parameters
        ----------
        data_path : str or Path
            Path to trait data CSV or Parquet
        table_name : str
            Table name in DuckDB

        Returns
        -------
        dict
            Dataset metadata
        """
        data_path = Path(data_path)

        if data_path.suffix == '.csv':
            metadata = self.db_handler.load_csv(data_path, table_name)
        elif data_path.suffix == '.parquet':
            metadata = self.db_handler.load_parquet(data_path, table_name)
        else:
            raise ValueError(f"Unsupported file format: {data_path.suffix}")

        logger.info(f"Loaded {metadata['row_count']:,} taxa with traits")
        return metadata

    def query_taxa_by_trait(
        self,
        table_name: str = "trait_data",
        trait_name: Optional[str] = None,
        trait_value: Optional[Any] = None,
        clade: Optional[str] = None,
        page: int = 1,
        page_size: int = 50
    ) -> Dict[str, Any]:
        """
        Query taxa by trait values.

        Parameters
        ----------
        table_name : str
            Table name
        trait_name : str, optional
            Trait column name
        trait_value : optional
            Trait value to filter
        clade : str, optional
            Clade/lineage filter
        page : int
            Page number
        page_size : int
            Rows per page

        Returns
        -------
        dict
            Paginated query results
        """
        where_clauses = []

        if trait_name and trait_value is not None:
            if isinstance(trait_value, str):
                where_clauses.append(f"{trait_name} = '{trait_value}'")
            else:
                where_clauses.append(f"{trait_name} = {trait_value}")

        if clade:
            where_clauses.append(f"clade = '{clade}'")

        where = " AND ".join(where_clauses) if where_clauses else None

        return self.db_handler.query_table(
            table_name,
            where=where,
            page=page,
            page_size=page_size
        )

    def build_tree_from_newick(
        self,
        newick_path: Union[str, Path]
    ) -> nx.DiGraph:
        """
        Build NetworkX tree from Newick file.

        Parameters
        ----------
        newick_path : str or Path
            Path to Newick format tree file

        Returns
        -------
        nx.DiGraph
            Directed graph representing tree
        """
        try:
            from Bio import Phylo
            from io import StringIO
        except ImportError:
            raise ImportError("Biopython required. Install with: pip install biopython")

        newick_path = Path(newick_path)

        # Read Newick tree
        tree = Phylo.read(newick_path, "newick")

        # Convert to NetworkX
        G = nx.DiGraph()

        def add_edges(clade, parent=None):
            node_name = clade.name or f"internal_{id(clade)}"

            if parent:
                G.add_edge(parent, node_name, length=clade.branch_length or 0)

            for child in clade.clades:
                add_edges(child, node_name)

        add_edges(tree.root)

        logger.info(f"Built tree: {G.number_of_nodes()} nodes")
        return G

    def render_large_tree(
        self,
        G: nx.DiGraph,
        output_path: Optional[Path] = None,
        node_labels: Optional[Dict] = None
    ) -> Path:
        """
        Render large phylogenetic tree.

        Parameters
        ----------
        G : nx.DiGraph
            Tree as directed graph
        output_path : Path, optional
            Output image path
        node_labels : dict, optional
            Node labels to display

        Returns
        -------
        Path
            Path to saved image
        """
        if output_path is None:
            output_path = self.output_dir / "phylogenetic_tree.png"

        logger.info(f"Rendering tree with {G.number_of_nodes()} taxa")

        return self.layout_renderer.layout_hierarchical(
            G,
            output_path=output_path,
            return_base64=False
        )

    def visualize_trait_on_tree(
        self,
        df: pd.DataFrame,
        x_col: str,
        y_col: str,
        trait_col: str,
        output_path: Optional[Path] = None
    ) -> Path:
        """
        Visualize trait distribution on phylogenetic tree coordinates.

        Parameters
        ----------
        df : pd.DataFrame
            Taxa data with tree coordinates and traits
        x_col : str
            X-axis column (tree coordinate)
        y_col : str
            Y-axis column (tree coordinate)
        trait_col : str
            Trait column for coloring
        output_path : Path, optional
            Output image path

        Returns
        -------
        Path
            Path to saved image
        """
        if output_path is None:
            output_path = self.output_dir / "trait_on_tree.png"

        logger.info(f"Visualizing {trait_col} on tree ({len(df):,} taxa)")

        return self.plotter.scatter_plot(
            df, x_col, y_col,
            color_by=trait_col,
            aggregation='count',
            output_path=output_path
        )

    def create_trait_correlation_heatmap(
        self,
        df: pd.DataFrame,
        trait_columns: Optional[List[str]] = None,
        output_dir: Optional[Path] = None
    ) -> Dict[str, Any]:
        """
        Create correlation heatmap for traits (tiled for many traits).

        Parameters
        ----------
        df : pd.DataFrame
            Trait data
        trait_columns : list of str, optional
            Trait columns to include
        output_dir : Path, optional
            Output directory for tiles

        Returns
        -------
        dict
            Heatmap metadata
        """
        if output_dir is None:
            output_dir = self.output_dir / "trait_correlation_tiles"

        if trait_columns:
            df_traits = df[trait_columns]
        else:
            df_traits = df.select_dtypes(include=[np.number])

        logger.info(f"Creating trait correlation heatmap: {df_traits.shape[1]} traits")

        from shared.holoviews_heatmaps import create_correlation_heatmap
        return create_correlation_heatmap(
            df_traits,
            method='spearman',
            output_dir=output_dir
        )

    def export_clade_statistics(
        self,
        table_name: str = "trait_data",
        clade_col: str = "clade",
        output_path: Optional[Path] = None
    ) -> Path:
        """
        Export clade-level statistics.

        Parameters
        ----------
        table_name : str
            Table name
        clade_col : str
            Clade column name
        output_path : Path, optional
            Output CSV path

        Returns
        -------
        Path
            Path to saved CSV
        """
        if output_path is None:
            output_path = self.output_dir / "clade_statistics.csv"

        query = f"""
        SELECT
            {clade_col},
            COUNT(*) as num_taxa,
            COUNT(DISTINCT serotype) as num_serotypes,
            AVG(CAST(pathogenicity_score AS FLOAT)) as avg_pathogenicity
        FROM {table_name}
        GROUP BY {clade_col}
        ORDER BY num_taxa DESC
        """

        return self.db_handler.export_filtered(query, output_path, format='csv')

    def calculate_phylogenetic_signal(
        self,
        G: nx.DiGraph,
        trait_data: pd.DataFrame,
        trait_col: str
    ) -> Dict[str, float]:
        """
        Calculate phylogenetic signal for a trait.

        Parameters
        ----------
        G : nx.DiGraph
            Phylogenetic tree
        trait_data : pd.DataFrame
            Trait data with taxa as index
        trait_col : str
            Trait column name

        Returns
        -------
        dict
            Phylogenetic signal statistics
        """
        logger.info(f"Calculating phylogenetic signal for {trait_col}")

        # Get tips (leaf nodes)
        tips = [n for n in G.nodes() if G.out_degree(n) == 0]

        # Filter trait data to tips
        tip_traits = trait_data.loc[trait_data.index.isin(tips), trait_col]

        # Simple phylogenetic signal: variance within vs between clades
        # This is a simplified version - proper implementation would use Pagel's lambda, etc.

        stats = {
            'num_tips': len(tips),
            'trait_mean': float(tip_traits.mean()),
            'trait_std': float(tip_traits.std()),
            'trait_range': float(tip_traits.max() - tip_traits.min())
        }

        return stats

    def close(self):
        """Close database connection."""
        self.db_handler.close()

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()


def process_large_phylotrait_dataset(
    trait_data_path: Union[str, Path],
    output_dir: Union[str, Path],
    tree_path: Optional[Union[str, Path]] = None
) -> Dict[str, Any]:
    """
    Convenience function to process large phylogenetic trait dataset.

    Parameters
    ----------
    trait_data_path : str or Path
        Path to trait data
    output_dir : str or Path
        Output directory
    tree_path : str or Path, optional
        Path to Newick tree file

    Returns
    -------
    dict
        Processing results and output paths
    """
    output_dir = Path(output_dir)

    with LargeDataPhyloTrait(output_dir) as processor:
        # Load trait data
        metadata = processor.load_trait_data(trait_data_path)

        results = {
            'metadata': metadata,
            'outputs': {}
        }

        # Export clade statistics
        clade_stats = processor.export_clade_statistics(
            output_path=output_dir / "clade_statistics.csv"
        )
        results['outputs']['clade_statistics'] = str(clade_stats)

        # Build and render tree if provided
        if tree_path:
            G = processor.build_tree_from_newick(tree_path)

            tree_viz = processor.render_large_tree(
                G,
                output_path=output_dir / "phylogenetic_tree.png"
            )
            results['outputs']['tree_visualization'] = str(tree_viz)
            results['num_taxa'] = G.number_of_nodes()

        return results


if __name__ == "__main__":
    # Example usage
    logging.basicConfig(level=logging.INFO)

    print("Large Data Integration - PhyloTrait Module")
    print("=" * 60)

    if not LARGE_DATA_AVAILABLE:
        print("ERROR: Large data pipeline not available")
        print("Install with: pip install duckdb matplotlib networkx datashader holoviews")
        exit(1)

    # Create sample large trait dataset
    print("\nCreating sample dataset (5,000 taxa)...")
    np.random.seed(42)

    n_taxa = 5000
    n_clades = 20

    data = {
        'taxon_id': [f"TAXON_{i:05d}" for i in range(n_taxa)],
        'clade': [f"Clade_{i % n_clades}" for i in range(n_taxa)],
        'serotype': np.random.choice(['1', '2', '3', '4', '5', '7', '9'], n_taxa),
        'pathogenicity_score': np.random.uniform(0, 10, n_taxa),
        'host_specificity': np.random.choice(['human', 'swine', 'generalist'], n_taxa),
        'tree_x': np.random.randn(n_taxa),  # Simulated tree coordinates
        'tree_y': np.random.randn(n_taxa),
        'trait1': np.random.randn(n_taxa),
        'trait2': np.random.randn(n_taxa),
        'trait3': np.random.randn(n_taxa)
    }

    df = pd.DataFrame(data)
    test_file = Path("test_phylotrait_large.csv")
    df.to_csv(test_file, index=False)
    print(f"Created: {test_file}")

    # Process dataset
    print("\nProcessing large phylogenetic trait dataset...")
    results = process_large_phylotrait_dataset(
        test_file,
        "test_phylotrait_output"
    )

    print(f"\n✓ Processing complete")
    print(f"Loaded: {results['metadata']['row_count']:,} taxa")
    print(f"Outputs: {list(results['outputs'].keys())}")
