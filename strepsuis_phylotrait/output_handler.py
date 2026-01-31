"""
Output handler for StrepSuis-PhyloTrait with StandardOutput integration.

This module wraps all output operations to ensure standardized formatting
with statistical interpretation and QA checklists.
"""

import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import pandas as pd

# Import StandardOutput from shared module
try:
    from shared import StandardOutput, create_qa_checklist
except ImportError:
    # Fallback if shared module not in path
    import sys
    shared_path = Path(__file__).parent.parent.parent / "shared" / "src"
    sys.path.insert(0, str(shared_path))
    from shared import StandardOutput, create_qa_checklist


class PhyloTraitOutputHandler:
    """
    Output handler for PhyloTrait analysis results.

    Wraps all output operations with StandardOutput to ensure consistent
    formatting and inclusion of statistical interpretations and QA checks.
    """

    def __init__(self, output_dir: str):
        """
        Initialize output handler.

        Args:
            output_dir: Directory for output files
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def save_phylogenetic_signal(
        self,
        signal_df: pd.DataFrame,
        n_traits: int
    ) -> None:
        """
        Save phylogenetic signal analysis with interpretation.

        Args:
            signal_df: DataFrame with phylogenetic signal metrics
            n_traits: Number of traits analyzed
        """
        output = StandardOutput(data=signal_df)

        # Analyze phylogenetic signal
        if "Fritz_D" in signal_df.columns:
            avg_d = signal_df["Fritz_D"].mean()
            significant_traits = (signal_df["P_Value"] < 0.05).sum() if "P_Value" in signal_df.columns else 0

            interp = f"Phylogenetic signal analysis of {n_traits} binary traits. "

            if avg_d < 0.5:
                interp += f"Average Fritz's D = {avg_d:.3f} indicates strong phylogenetic signal "
                interp += f"(traits are more conserved than expected under Brownian motion). "
            elif avg_d < 1.0:
                interp += f"Average Fritz's D = {avg_d:.3f} indicates moderate phylogenetic signal. "
            else:
                interp += f"Average Fritz's D = {avg_d:.3f} suggests traits evolve independently of phylogeny. "

            if significant_traits > 0:
                interp += f"{significant_traits} traits show significant phylogenetic signal (p < 0.05), "
                interp += f"indicating evolutionary constraint or shared ancestry effects."
        else:
            interp = f"Phylogenetic signal analysis completed for {n_traits} traits."

        output.add_statistical_interpretation(interp)

        # QA checklist
        qa_items = [
            f"✓ Phylogenetic signal analyzed for {n_traits} traits",
        ]

        if "Fritz_D" in signal_df.columns:
            qa_items.append(f"✓ Fritz's D calculated for all traits")

        if "P_Value" in signal_df.columns and significant_traits > 0:
            qa_items.append(f"✓ {significant_traits} traits with significant signal (p < 0.05)")

        output.add_quick_qa(qa_items)
        output.add_metadata("n_traits", n_traits)

        base_path = self.output_dir / "phylogenetic_signal"
        output.to_json(f"{base_path}.json")
        output.to_csv(f"{base_path}.csv")
        output.to_markdown(f"{base_path}.md")

    def save_trait_evolution(
        self,
        evolution_df: pd.DataFrame,
        model: str = "ER"
    ) -> None:
        """
        Save trait evolution analysis with interpretation.

        Args:
            evolution_df: DataFrame with trait evolution results
            model: Evolution model used (ER, SYM, ARD)
        """
        output = StandardOutput(data=evolution_df)

        # Analyze transition rates
        if "Transition_Rate" in evolution_df.columns:
            avg_rate = evolution_df["Transition_Rate"].mean()
            rate_range = evolution_df["Transition_Rate"].max() - evolution_df["Transition_Rate"].min()

            interp = f"Trait evolution analysis using {model} model. "
            interp += f"Average transition rate: {avg_rate:.4f}. "

            if model == "ARD":
                interp += f"Asymmetric rates (range: {rate_range:.4f}) suggest directional selection "
                interp += f"or biased gain/loss of traits."
            elif model == "SYM":
                interp += f"Symmetric transition rates indicate balanced trait evolution."
            else:
                interp += f"Equal rates model assumes uniform trait evolution across phylogeny."

            if "AIC" in evolution_df.columns:
                best_model_idx = evolution_df["AIC"].idxmin()
                best_model = evolution_df.loc[best_model_idx, "Model"] if "Model" in evolution_df.columns else model
                interp += f" Best-fitting model: {best_model}."
        else:
            interp = f"Trait evolution analysis completed using {model} model."

        output.add_statistical_interpretation(interp)

        # QA checklist
        qa_items = [
            f"✓ Trait evolution analyzed using {model} model",
        ]

        if "Transition_Rate" in evolution_df.columns:
            qa_items.append(f"✓ Transition rates calculated")

        if "AIC" in evolution_df.columns:
            qa_items.append("✓ Model comparison via AIC completed")

        output.add_quick_qa(qa_items)
        output.add_metadata("model", model)

        base_path = self.output_dir / "trait_evolution"
        output.to_json(f"{base_path}.json")
        output.to_csv(f"{base_path}.csv")
        output.to_markdown(f"{base_path}.md")

    def save_phylogenetic_clustering(
        self,
        clustering_df: pd.DataFrame,
        n_taxa: int
    ) -> None:
        """
        Save phylogenetic clustering results with interpretation.

        Args:
            clustering_df: DataFrame with clustering results
            n_taxa: Number of taxa analyzed
        """
        output = StandardOutput(data=clustering_df)

        # Analyze clustering patterns
        if "Net_Relatedness_Index" in clustering_df.columns:
            avg_nri = clustering_df["Net_Relatedness_Index"].mean()

            interp = f"Phylogenetic clustering analysis of {n_taxa} taxa. "

            if avg_nri > 2:
                interp += f"Net Relatedness Index (NRI) = {avg_nri:.2f} indicates strong phylogenetic "
                interp += f"clustering (traits concentrated in related lineages). "
            elif avg_nri > 0:
                interp += f"NRI = {avg_nri:.2f} suggests moderate phylogenetic clustering. "
            elif avg_nri > -2:
                interp += f"NRI = {avg_nri:.2f} indicates random phylogenetic distribution. "
            else:
                interp += f"NRI = {avg_nri:.2f} suggests phylogenetic overdispersion "
                interp += f"(traits avoid close relatives). "

            if "P_Value" in clustering_df.columns:
                significant = (clustering_df["P_Value"] < 0.05).sum()
                if significant > 0:
                    interp += f"{significant} traits show significant clustering (p < 0.05)."
        else:
            interp = f"Phylogenetic clustering analysis completed for {n_taxa} taxa."

        output.add_statistical_interpretation(interp)

        # QA checklist
        qa_items = [
            f"✓ Phylogenetic clustering analyzed for {n_taxa} taxa",
        ]

        if "Net_Relatedness_Index" in clustering_df.columns:
            qa_items.append("✓ NRI calculated for trait distribution")

        if "P_Value" in clustering_df.columns:
            significant = (clustering_df["P_Value"] < 0.05).sum()
            if significant > 0:
                qa_items.append(f"✓ {significant} traits with significant clustering")

        output.add_quick_qa(qa_items)
        output.add_metadata("n_taxa", n_taxa)

        base_path = self.output_dir / "phylogenetic_clustering"
        output.to_json(f"{base_path}.json")
        output.to_csv(f"{base_path}.csv")
        output.to_markdown(f"{base_path}.md")

    def save_ancestral_state_reconstruction(
        self,
        ancestral_df: pd.DataFrame,
        n_nodes: int
    ) -> None:
        """
        Save ancestral state reconstruction with interpretation.

        Args:
            ancestral_df: DataFrame with reconstructed ancestral states
            n_nodes: Number of internal nodes
        """
        output = StandardOutput(data=ancestral_df)

        # Analyze ancestral states
        if "Probability_1" in ancestral_df.columns:
            high_confidence = (ancestral_df["Probability_1"] > 0.8).sum() + (ancestral_df["Probability_1"] < 0.2).sum()
            pct_confident = 100 * high_confidence / len(ancestral_df) if len(ancestral_df) > 0 else 0

            interp = f"Ancestral state reconstruction for {n_nodes} internal nodes. "
            interp += f"{high_confidence}/{len(ancestral_df)} ({pct_confident:.1f}%) reconstructions "
            interp += f"show high confidence (probability >80%). "

            if "Trait" in ancestral_df.columns:
                traits_present = ancestral_df["Trait"].nunique()
                interp += f"Analyzed {traits_present} trait(s) across phylogeny, "
                interp += f"revealing evolutionary origins and state transitions."
        else:
            interp = f"Ancestral state reconstruction completed for {n_nodes} nodes."

        output.add_statistical_interpretation(interp)

        # QA checklist
        qa_items = [
            f"✓ Ancestral states reconstructed for {n_nodes} internal nodes",
        ]

        if "Probability_1" in ancestral_df.columns:
            if pct_confident > 70:
                qa_items.append(f"✓ High confidence reconstructions: {pct_confident:.1f}%")
            else:
                qa_items.append(f"⚠ Low confidence reconstructions: {pct_confident:.1f}%")

        output.add_quick_qa(qa_items)
        output.add_metadata("n_nodes", n_nodes)

        base_path = self.output_dir / "ancestral_state_reconstruction"
        output.to_json(f"{base_path}.json")
        output.to_csv(f"{base_path}.csv")
        output.to_markdown(f"{base_path}.md")

    def save_generic_results(
        self,
        df: pd.DataFrame,
        filename: str,
        interpretation: str,
        qa_items: List[str],
        metadata: Optional[Dict[str, Any]] = None
    ) -> None:
        """
        Save generic results with custom interpretation and QA.

        Args:
            df: DataFrame with results
            filename: Base filename (without extension)
            interpretation: Statistical interpretation text
            qa_items: QA checklist items
            metadata: Optional additional metadata
        """
        output = StandardOutput(data=df)
        output.add_statistical_interpretation(interpretation)
        output.add_quick_qa(qa_items)

        if metadata:
            for key, value in metadata.items():
                output.add_metadata(key, value)

        base_path = self.output_dir / filename
        output.to_json(f"{base_path}.json")
        output.to_csv(f"{base_path}.csv")
        output.to_markdown(f"{base_path}.md")
