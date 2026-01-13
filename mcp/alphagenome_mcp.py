"""
Model Context Protocol (MCP) for alphagenome

AlphaGenome is a machine learning model that predicts the functional effects of genetic variants across
diverse genomic modalities including gene expression, chromatin accessibility, histone modifications,
and 3D genome organization. This MCP server provides unified access to tools for variant scoring,
genomic predictions, visualization, and comprehensive analysis workflows.

This MCP Server contains tools extracted from the following tutorial files:
1. batch_variant_scoring
    - alphagenome_score_batch_variants: Score multiple genetic variants across genomic modalities
2. example_analysis_workflow
    - alphagenome_visualize_variant_positions: Visualizes genomic positions of variants near TAL1 gene
    - alphagenome_predict_variant_effects: Predicts functional impact of a variant on gene expression, accessibility, and histone marks
    - alphagenome_compare_variant_scores: Compares predicted TAL1 expression changes between cancer-associated and background variants
3. variant_scoring_ui
    - alphagenome_score_variant: Score a single genetic variant across multiple genomic features
    - alphagenome_visualize_variant_effects: Visualize predicted variant effects across genomic modalities
4. visualization_modality_tour
    - alphagenome_visualize_gene_expression: Visualize RNA-seq and CAGE predictions for gene expression
    - alphagenome_visualize_variant_effect_expression: Visualize variant effects on gene expression
    - alphagenome_visualize_polyadenylation_sites: Visualize custom polyadenylation site annotations
    - alphagenome_visualize_chromatin_accessibility: Visualize DNASE and ATAC chromatin accessibility
    - alphagenome_visualize_splicing: Visualize splicing predictions with junction arcs
    - alphagenome_visualize_chip_histone: Visualize histone modification markers
    - alphagenome_visualize_chip_tf: Visualize transcription factor binding predictions
    - alphagenome_visualize_contact_maps: Visualize DNA-DNA contact map predictions
"""

from fastmcp import FastMCP

from tools.batch_variant_scoring import batch_variant_scoring_mcp
from tools.example_analysis_workflow import example_analysis_workflow_mcp
from tools.variant_scoring_ui import variant_scoring_ui_mcp
from tools.visualization_modality_tour import visualization_modality_tour_mcp

mcp = FastMCP(name="alphagenome")
mcp.mount(batch_variant_scoring_mcp)
mcp.mount(example_analysis_workflow_mcp)
mcp.mount(variant_scoring_ui_mcp)
mcp.mount(visualization_modality_tour_mcp)

if __name__ == "__main__":
    mcp.run()
