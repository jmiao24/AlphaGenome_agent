# [Paper2Agent](https://github.com/jmiao24/Paper2Agent): AlphaGenome Demo

A demonstration of turning the [AlphaGenome paper](https://deepmind.google.com/science/alphagenome/) into an interactive AI agent. This project transforms Google DeepMind's unified DNA sequence-to-function model into a conversational agent that can analyze genetic variants, predict functional effects, and generate visualizations through natural language.

## Folder Structure

```
AlphaGenome_demo/
├── manuscript/
│   └── manuscript.md          # Full AlphaGenome paper manuscript
├── mcp/
│   ├── alphagenome_mcp.py     # MCP server entry point
│   ├── requirements.txt       # Python dependencies
│   └── tools/
│       ├── batch_variant_scoring.py         # Batch variant scoring tools
│       ├── example_analysis_workflow.py     # Analysis workflow examples
│       ├── variant_scoring_ui.py            # Single variant scoring UI
│       └── visualization_modality_tour.py   # Visualization tools
└── supp/
    └── supp_table.xlsx        # Supplementary tables
```

## Quick Start

### 1. Clone the Repository

```bash
git clone https://github.com/jmiao24/AlphaGenome_demo.git
cd AlphaGenome_demo
```

### 2. Install Gemini CLI

Install the [Google Gemini CLI](https://github.com/google-gemini/gemini-cli):

```bash
brew install gemini-cli
```

### 3. Obtain API Key

AlphaGenome requires an API key to access the model. Obtain your API key at:

https://deepmind.google.com/science/alphagenome/

### 4. Install FastMCP

```bash
pip install fastmcp
```

### 5. Install MCP Server

Install the AlphaGenome MCP server using fastmcp:

```bash
fastmcp install gemini-cli ./mcp/alphagenome_mcp.py --with-requirements ./mcp/requirements.txt
```

### 6. Start the Agent

Start Gemini CLI in the repository folder:

```bash
gemini
```

You will now have access to the AlphaGenome agent with all available tools.

## Example Query

```
Analyze heart gene expression data with AlphaGenome MCP to identify the causal gene
for the variant chr11:116837649:T>G, associated with Hypoalphalipoproteinemia.
```

## Available Agent Tools

The agent provides the following capabilities through natural language:

### Batch Variant Scoring
- `alphagenome_score_batch_variants`: Score multiple genetic variants across genomic modalities

### Analysis Workflow
- `alphagenome_visualize_variant_positions`: Visualize genomic positions of variants near genes
- `alphagenome_predict_variant_effects`: Predict functional impact on gene expression, accessibility, and histone marks
- `alphagenome_compare_variant_scores`: Compare predicted expression changes between variant groups

### Single Variant Scoring
- `alphagenome_score_variant`: Score a single genetic variant across multiple genomic features
- `alphagenome_visualize_variant_effects`: Visualize predicted variant effects across modalities

### Visualization Tools
- `alphagenome_visualize_gene_expression`: Visualize RNA-seq and CAGE predictions
- `alphagenome_visualize_variant_effect_expression`: Visualize variant effects on gene expression
- `alphagenome_visualize_polyadenylation_sites`: Visualize polyadenylation site annotations
- `alphagenome_visualize_chromatin_accessibility`: Visualize DNASE and ATAC chromatin accessibility
- `alphagenome_visualize_splicing`: Visualize splicing predictions with junction arcs
- `alphagenome_visualize_chip_histone`: Visualize histone modification markers
- `alphagenome_visualize_chip_tf`: Visualize transcription factor binding predictions
- `alphagenome_visualize_contact_maps`: Visualize DNA-DNA contact map predictions

## About AlphaGenome

AlphaGenome takes 1 megabase of DNA sequence as input and predicts thousands of functional genomic tracks up to single base pair resolution across diverse modalities:

- Gene expression (RNA-seq, CAGE, PRO-cap)
- Splicing (splice sites, splice site usage, splice junctions)
- Chromatin accessibility (DNase-seq, ATAC-seq)
- Histone modifications (ChIP-seq)
- Transcription factor binding (TF ChIP-seq)
- Chromatin contact maps (Hi-C/micro-C)

For more details, see the manuscript in `manuscript/manuscript.md`.
