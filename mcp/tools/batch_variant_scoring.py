"""
Batch variant scoring using AlphaGenome API to predict functional impacts.

This MCP Server provides 1 tool:
1. alphagenome_score_batch_variants: Score multiple genetic variants across genomic modalities

All tools extracted from `google-deepmind/alphagenome/colabs/batch_variant_scoring.ipynb`.
"""

# Standard imports
from typing import Annotated, Literal, Any
import pandas as pd
from pathlib import Path
import os
from fastmcp import FastMCP
from datetime import datetime

# AlphaGenome imports
from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers

# Project structure
PROJECT_ROOT = Path(__file__).parent.parent.parent.resolve()
DEFAULT_INPUT_DIR = PROJECT_ROOT / "tmp" / "inputs"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "tmp" / "outputs"

INPUT_DIR = Path(os.environ.get("BATCH_VARIANT_SCORING_INPUT_DIR", DEFAULT_INPUT_DIR))
OUTPUT_DIR = Path(os.environ.get("BATCH_VARIANT_SCORING_OUTPUT_DIR", DEFAULT_OUTPUT_DIR))

# Ensure directories exist
INPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Timestamp for unique outputs
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

# MCP server instance
batch_variant_scoring_mcp = FastMCP(name="batch_variant_scoring")


@batch_variant_scoring_mcp.tool()
def alphagenome_score_batch_variants(
    vcf_path: Annotated[str | None,
        "Path to VCF (Variant Call Format) file in tab-separated format. "
        "Must include columns: 'variant_id' (unique identifier), 'CHROM' (chromosome as 'chr1', 'chr2', etc.), "
        "'POS' (base pair position, 1-based, hg38 for human or mm10 for mouse), "
        "'REF' (reference nucleotide sequence), 'ALT' (alternate nucleotide sequence). "
        "Example row: chr3_58394738_A_T_b38\tchr3\t58394738\tA\tT"] = None,
    
    api_key: Annotated[str | None,
        "AlphaGenome API key for accessing the model. "
        "Required for all scoring operations. "
        "Obtain from Google AI Studio or AlphaGenome documentation."] = None,
    
    organism: Annotated[Literal["human", "mouse"],
        "Organism for variant scoring. "
        "'human': Uses hg38 reference genome (Homo sapiens). "
        "'mouse': Uses mm10 reference genome (Mus musculus). "
        "Affects available scorers and genomic context."] = "human",
    
    sequence_length: Annotated[Literal["16KB", "100KB", "500KB", "1MB"],
        "Length of genomic sequence context around each variant. "
        "Larger contexts capture long-range regulatory elements but require more computation. "
        "'16KB': Fast, local effects only. "
        "'1MB': Comprehensive, includes distal regulatory regions (tutorial default)."] = "1MB",
    
    score_rna_seq: Annotated[bool,
        "Score impact on RNA-seq gene expression across tissues and cell types."] = True,
    
    score_cage: Annotated[bool,
        "Score impact on CAGE (Cap Analysis Gene Expression) transcription start sites."] = True,
    
    score_procap: Annotated[bool,
        "Score impact on PRO-cap nascent RNA transcription (human only)."] = True,
    
    score_atac: Annotated[bool,
        "Score impact on ATAC-seq chromatin accessibility."] = True,
    
    score_dnase: Annotated[bool,
        "Score impact on DNase-seq chromatin accessibility."] = True,
    
    score_chip_histone: Annotated[bool,
        "Score impact on ChIP-seq histone modifications (H3K27ac, H3K27me3, H3K36me3, H3K4me1, H3K4me3, H3K9me3)."] = True,
    
    score_chip_tf: Annotated[bool,
        "Score impact on ChIP-seq transcription factor binding."] = True,
    
    score_polyadenylation: Annotated[bool,
        "Score impact on RNA polyadenylation sites."] = True,
    
    score_splice_sites: Annotated[bool,
        "Score impact on splice donor and acceptor sites."] = True,
    
    score_splice_site_usage: Annotated[bool,
        "Score impact on splice site usage frequencies."] = True,
    
    score_splice_junctions: Annotated[bool,
        "Score impact on exon-exon junction usage."] = True,
    
    out_prefix: Annotated[str | None,
        "Prefix for output file names. "
        "If None, uses 'batch_variant_scores_YYYYMMDD_HHMMSS'. "
        "Output CSV contains columns: variant_id, scored_interval, gene_id, gene_name, output_type, "
        "variant_scorer, biosample_type, raw_score, quantile_score, and metadata."] = None,
) -> dict:
    """
    Score multiple genetic variants across genomic modalities using AlphaGenome API.
    Input is VCF file with variant information and output is comprehensive scoring table with raw and quantile scores across selected modalities.
    """
    # Input validation
    if vcf_path is None:
        raise ValueError("Path to VCF file must be provided")
    
    if api_key is None:
        raise ValueError("AlphaGenome API key must be provided")
    
    # File existence validation
    vcf_file = Path(vcf_path)
    if not vcf_file.exists():
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")
    
    # Load VCF file
    vcf = pd.read_csv(vcf_path, sep='\t')
    
    # Validate required columns
    required_columns = ['variant_id', 'CHROM', 'POS', 'REF', 'ALT']
    for column in required_columns:
        if column not in vcf.columns:
            raise ValueError(f'VCF file is missing required column: {column}.')
    
    # Load AlphaGenome model
    dna_model = dna_client.create(api_key)
    
    # Parse organism specification
    organism_map = {
        'human': dna_client.Organism.HOMO_SAPIENS,
        'mouse': dna_client.Organism.MUS_MUSCULUS,
    }
    organism_enum = organism_map[organism]
    
    # Parse sequence length
    sequence_length_value = dna_client.SUPPORTED_SEQUENCE_LENGTHS[
        f'SEQUENCE_LENGTH_{sequence_length}'
    ]
    
    # Parse scorer specification
    scorer_selections = {
        'rna_seq': score_rna_seq,
        'cage': score_cage,
        'procap': score_procap,
        'atac': score_atac,
        'dnase': score_dnase,
        'chip_histone': score_chip_histone,
        'chip_tf': score_chip_tf,
        'polyadenylation': score_polyadenylation,
        'splice_sites': score_splice_sites,
        'splice_site_usage': score_splice_site_usage,
        'splice_junctions': score_splice_junctions,
    }
    
    all_scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    selected_scorers = [
        all_scorers[key]
        for key in all_scorers
        if scorer_selections.get(key.lower(), False)
    ]
    
    # Remove any scorers or output types that are not supported for the chosen organism
    unsupported_scorers = [
        scorer
        for scorer in selected_scorers
        if (
            organism_enum.value
            not in variant_scorers.SUPPORTED_ORGANISMS[scorer.base_variant_scorer]
        )
        | (
            (scorer.requested_output == dna_client.OutputType.PROCAP)
            & (organism_enum == dna_client.Organism.MUS_MUSCULUS)
        )
    ]
    if len(unsupported_scorers) > 0:
        print(
            f'Excluding {unsupported_scorers} scorers as they are not supported for'
            f' {organism_enum}.'
        )
        for unsupported_scorer in unsupported_scorers:
            selected_scorers.remove(unsupported_scorer)
    
    # Score variants in the VCF file
    results = []
    
    for i, vcf_row in vcf.iterrows():
        variant = genome.Variant(
            chromosome=str(vcf_row.CHROM),
            position=int(vcf_row.POS),
            reference_bases=vcf_row.REF,
            alternate_bases=vcf_row.ALT,
            name=vcf_row.variant_id,
        )
        interval = variant.reference_interval.resize(sequence_length_value)
        
        variant_scores = dna_model.score_variant(
            interval=interval,
            variant=variant,
            variant_scorers=selected_scorers,
            organism=organism_enum,
        )
        results.append(variant_scores)
    
    # Convert to tidy dataframe
    df_scores = variant_scorers.tidy_scores(results)
    
    # Determine output file prefix
    if out_prefix is None:
        out_prefix = f"batch_variant_scores_{timestamp}"
    
    # Save results
    output_file = OUTPUT_DIR / f"{out_prefix}.csv"
    df_scores.to_csv(output_file, index=False)
    
    # Return standardized format
    return {
        "message": f"Successfully scored {len(vcf)} variants across {len(selected_scorers)} modalities, generating {len(df_scores)} scores",
        "reference": "https://github.com/google-deepmind/alphagenome/blob/main/colabs/batch_variant_scoring.ipynb",
        "artifacts": [
            {
                "description": "Variant scores across genomic modalities",
                "path": str(output_file.resolve())
            }
        ]
    }
