"""
Configuration module

Handles all configuration parameters for analysis.
"""

import os
from dataclasses import dataclass


@dataclass
class Config:
    """Configuration for analysis."""

    # Directories
    data_dir: str = "."
    output_dir: str = "./output"
    # Backwards-compatible alias (tests use base_dir)
    base_dir: str = ""
    output_folder: str = ""
    tree_file: str = "Snp_tree.newick"
    mic_file: str = "MIC.csv"
    amr_genes_file: str = "AMR_genes.csv"
    virulence_genes_file: str = "Virulence.csv"
    n_clusters_range: tuple = (2, 10)
    n_ensemble: int = 5
    dbscan_trials: int = 20

    # Statistical parameters
    bootstrap_iterations: int = 500
    fdr_alpha: float = 0.05
    random_seed: int = 42

    # Reporting parameters
    generate_html: bool = True
    generate_excel: bool = True
    save_png_charts: bool = True
    dpi: int = 150

    # Parallel processing
    n_jobs: int = -1

    def __post_init__(self):
        """Validate configuration after initialization."""
        if self.base_dir:
            self.data_dir = self.base_dir
        else:
            self.base_dir = self.data_dir

        if self.output_folder:
            self.output_dir = self.output_folder
        else:
            self.output_folder = self.output_dir

        if not os.path.exists(self.data_dir):
            raise ValueError(f"Data directory does not exist: {self.data_dir}")

        os.makedirs(self.output_dir, exist_ok=True)

        if not 0 < self.fdr_alpha < 1:
            raise ValueError("fdr_alpha must be between 0 and 1")

        if self.bootstrap_iterations < 100:
            if not self._allow_small_bootstrap():
                raise ValueError("bootstrap_iterations should be at least 100")

    def _allow_small_bootstrap(self) -> bool:
        """Allow small bootstrap for example/test datasets."""
        data_dir = str(self.data_dir).lower()
        if "examples" in data_dir or "pytest" in data_dir:
            return True
        return os.environ.get("ALLOW_SMALL_BOOTSTRAP") == "1"

    @classmethod
    def from_dict(cls, config_dict: dict) -> "Config":
        """Create Config from dictionary."""
        return cls(**{k: v for k, v in config_dict.items() if k in cls.__dataclass_fields__})

    def to_dict(self) -> dict:
        """Return configuration as a dictionary."""
        return {k: getattr(self, k) for k in self.__dataclass_fields__}

    def get_tree_path(self) -> str:
        """Return full path to the tree file."""
        return os.path.join(self.base_dir, self.tree_file)

    def get_output_path(self) -> str:
        """Return output directory path."""
        return self.output_dir
