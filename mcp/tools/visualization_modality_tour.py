"""
Comprehensive visualization tutorial for AlphaGenome model predictions across different genomic modalities.

This MCP Server provides 8 tools for visualizing AlphaGenome predictions:
1. alphagenome_visualize_gene_expression: Visualize RNA-seq and CAGE predictions for gene expression
2. alphagenome_visualize_variant_effect_expression: Visualize variant effects on gene expression
3. alphagenome_visualize_polyadenylation_sites: Visualize custom polyadenylation site annotations
4. alphagenome_visualize_chromatin_accessibility: Visualize DNASE and ATAC chromatin accessibility
5. alphagenome_visualize_splicing: Visualize splicing predictions with junction arcs
6. alphagenome_visualize_chip_histone: Visualize histone modification markers
7. alphagenome_visualize_chip_tf: Visualize transcription factor binding predictions
8. alphagenome_visualize_contact_maps: Visualize DNA-DNA contact map predictions

All tools extracted from google-deepmind/alphagenome visualization_modality_tour tutorial.
"""

# Standard imports
from typing import Annotated, Literal, Any
import pandas as pd
import numpy as np
from pathlib import Path
import os
from fastmcp import FastMCP
from datetime import datetime
import matplotlib.pyplot as plt

# AlphaGenome imports
from alphagenome.data import gene_annotation, genome, track_data, transcript
from alphagenome.models import dna_client
from alphagenome.visualization import plot_components

# Project structure
PROJECT_ROOT = Path(__file__).parent.parent.parent.resolve()
DEFAULT_INPUT_DIR = PROJECT_ROOT / "tmp" / "inputs"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "tmp" / "outputs"

INPUT_DIR = Path(os.environ.get("VISUALIZATION_MODALITY_TOUR_INPUT_DIR", DEFAULT_INPUT_DIR))
OUTPUT_DIR = Path(os.environ.get("VISUALIZATION_MODALITY_TOUR_OUTPUT_DIR", DEFAULT_OUTPUT_DIR))

# Ensure directories exist
INPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Configure matplotlib for high-resolution figures
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

# Timestamp for unique outputs
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

# MCP server instance
visualization_modality_tour_mcp = FastMCP(name="visualization_modality_tour")


@visualization_modality_tour_mcp.tool()
def alphagenome_visualize_gene_expression(
    api_key: Annotated[str,
        "AlphaGenome API key for making predictions. "
        "Required for accessing the DNA model and making interval predictions."],
    chromosome: Annotated[str,
        "Chromosome identifier (e.g., 'chr22'). "
        "Must be a valid chromosome name in the reference genome."] = "chr22",
    start_position: Annotated[int,
        "Start position of the genomic interval (0-based). "
        "Example: 36150498 for a region on chr22."] = 36150498,
    end_position: Annotated[int,
        "End position of the genomic interval (0-based). "
        "Example: 36252898. The interval width must be one of the supported sequence lengths."] = 36252898,
    ontology_terms: Annotated[list[str],
        "List of tissue/cell-type ontology IDs to predict expression for. "
        "Examples: ['UBERON:0001159', 'UBERON:0001155'] for sigmoid and transverse colon. "
        "Use UBERON or other standard ontology identifiers."] = ["UBERON:0001159", "UBERON:0001155"],
    out_prefix: Annotated[str | None,
        "Prefix for output file names. If None, uses 'gene_expression_{timestamp}'."] = None,
) -> dict:
    """
    Visualize predicted gene expression (RNA-seq and CAGE) for a genomic interval across tissues.
    Input is genomic coordinates and tissue ontology terms, output is visualization figure showing RNA expression predictions.
    """
    # Set output prefix
    if out_prefix is None:
        out_prefix = f"gene_expression_{timestamp}"
    
    # Create DNA model client
    dna_model = dna_client.create(api_key)
    
    # Load gene annotations from GENCODE
    gtf = pd.read_feather(
        'https://storage.googleapis.com/alphagenome/reference/gencode/'
        'hg38/gencode.v46.annotation.gtf.gz.feather'
    )
    
    # Filter to protein-coding genes and highly supported transcripts
    gtf_transcript = gene_annotation.filter_transcript_support_level(
        gene_annotation.filter_protein_coding(gtf), ['1']
    )
    
    # Extractor for longest transcript per gene
    gtf_longest_transcript = gene_annotation.filter_to_longest_transcript(gtf_transcript)
    longest_transcript_extractor = transcript.TranscriptExtractor(gtf_longest_transcript)
    
    # Define interval and resize to supported sequence length
    interval = genome.Interval(chromosome, start_position, end_position).resize(
        dna_client.SEQUENCE_LENGTH_1MB
    )
    
    # Make predictions
    output = dna_model.predict_interval(
        interval=interval,
        requested_outputs={
            dna_client.OutputType.RNA_SEQ,
            dna_client.OutputType.CAGE,
        },
        ontology_terms=ontology_terms,
    )
    
    # Extract longest transcripts for this interval
    longest_transcripts = longest_transcript_extractor.extract(interval)
    
    # Build plot
    plot = plot_components.plot(
        [
            plot_components.TranscriptAnnotation(longest_transcripts),
            plot_components.Tracks(
                tdata=output.rna_seq,
                ylabel_template='RNA_SEQ: {biosample_name} ({strand})\n{name}',
            ),
            plot_components.Tracks(
                tdata=output.cage,
                ylabel_template='CAGE: {biosample_name} ({strand})\n{name}',
            ),
        ],
        interval=interval,
        title='Predicted RNA Expression (RNA_SEQ, CAGE) for colon tissue',
    )
    
    # Save figure
    figure_path = OUTPUT_DIR / f"{out_prefix}_gene_expression.png"
    plt.savefig(figure_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return {
        "message": "Gene expression visualization completed successfully",
        "reference": "https://github.com/google-deepmind/alphagenome/blob/main/colabs/visualization_modality_tour.ipynb",
        "artifacts": [
            {
                "description": "Gene expression visualization figure",
                "path": str(figure_path.resolve())
            }
        ]
    }


@visualization_modality_tour_mcp.tool()
def alphagenome_visualize_variant_effect_expression(
    api_key: Annotated[str,
        "AlphaGenome API key for making predictions. "
        "Required for accessing the DNA model and making variant predictions."],
    variant_string: Annotated[str,
        "Variant in format 'chr:pos:ref>alt' (e.g., 'chr22:36201698:A>C'). "
        "Specifies the variant to visualize effect for."],
    chromosome: Annotated[str,
        "Chromosome identifier for the prediction interval (e.g., 'chr22')."] = "chr22",
    start_position: Annotated[int,
        "Start position of the genomic interval (0-based)."] = 36150498,
    end_position: Annotated[int,
        "End position of the genomic interval (0-based)."] = 36252898,
    gene_symbol: Annotated[str,
        "Gene symbol to zoom into for visualization (e.g., 'APOL4'). "
        "The plot will focus on the region around this gene."] = "APOL4",
    zoom_offset: Annotated[int,
        "Number of base pairs to add on either side of the gene body for visualization. "
        "Typical values: 500-2000."] = 1000,
    ontology_terms: Annotated[list[str],
        "List of tissue/cell-type ontology IDs. "
        "Examples: ['UBERON:0001159', 'UBERON:0001155'] for colon tissues."] = ["UBERON:0001159", "UBERON:0001155"],
    out_prefix: Annotated[str | None,
        "Prefix for output file names. If None, uses 'variant_effect_{timestamp}'."] = None,
) -> dict:
    """
    Visualize the effect of a genetic variant on predicted gene expression (RNA-seq and CAGE).
    Input is variant specification and genomic interval, output is REF vs ALT comparison visualization.
    """
    # Set output prefix
    if out_prefix is None:
        out_prefix = f"variant_effect_{timestamp}"
    
    # Create DNA model client
    dna_model = dna_client.create(api_key)
    
    # Load gene annotations
    gtf = pd.read_feather(
        'https://storage.googleapis.com/alphagenome/reference/gencode/'
        'hg38/gencode.v46.annotation.gtf.gz.feather'
    )
    
    gtf_transcript = gene_annotation.filter_transcript_support_level(
        gene_annotation.filter_protein_coding(gtf), ['1']
    )
    
    gtf_longest_transcript = gene_annotation.filter_to_longest_transcript(gtf_transcript)
    longest_transcript_extractor = transcript.TranscriptExtractor(gtf_longest_transcript)
    
    # Define interval and variant
    interval = genome.Interval(chromosome, start_position, end_position).resize(
        dna_client.SEQUENCE_LENGTH_1MB
    )
    variant = genome.Variant.from_str(variant_string)
    
    # Make predictions for REF and ALT alleles
    output = dna_model.predict_variant(
        interval=interval,
        variant=variant,
        requested_outputs={
            dna_client.OutputType.RNA_SEQ,
            dna_client.OutputType.CAGE,
        },
        ontology_terms=ontology_terms,
    )
    
    # Get longest transcripts
    longest_transcripts = longest_transcript_extractor.extract(interval)
    
    # Zoom in on the gene region
    gene_interval = gene_annotation.get_gene_interval(gtf, gene_symbol=gene_symbol)
    gene_interval.resize_inplace(gene_interval.width + zoom_offset)
    
    # Define colors for REF and ALT
    ref_alt_colors = {'REF': 'dimgrey', 'ALT': 'red'}
    
    # Build plot
    plot = plot_components.plot(
        [
            plot_components.TranscriptAnnotation(longest_transcripts),
            # RNA-seq tracks
            plot_components.OverlaidTracks(
                tdata={
                    'REF': output.reference.rna_seq.filter_to_nonpositive_strand(),
                    'ALT': output.alternate.rna_seq.filter_to_nonpositive_strand(),
                },
                colors=ref_alt_colors,
                ylabel_template='{biosample_name} ({strand})\n{name}',
            ),
            # CAGE track
            plot_components.OverlaidTracks(
                tdata={
                    'REF': output.reference.cage.filter_to_nonpositive_strand(),
                    'ALT': output.alternate.cage.filter_to_nonpositive_strand(),
                },
                colors=ref_alt_colors,
                ylabel_template='{biosample_name} ({strand})\n{name}',
            ),
        ],
        annotations=[plot_components.VariantAnnotation([variant])],
        interval=gene_interval,
        title='Effect of variant on predicted RNA Expression in colon tissue',
    )
    
    # Save figure
    figure_path = OUTPUT_DIR / f"{out_prefix}_variant_effect.png"
    plt.savefig(figure_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return {
        "message": "Variant effect visualization completed successfully",
        "reference": "https://github.com/google-deepmind/alphagenome/blob/main/colabs/visualization_modality_tour.ipynb",
        "artifacts": [
            {
                "description": "Variant effect visualization figure",
                "path": str(figure_path.resolve())
            }
        ]
    }


@visualization_modality_tour_mcp.tool()
def alphagenome_visualize_polyadenylation_sites(
    api_key: Annotated[str,
        "AlphaGenome API key for making predictions."],
    pa_site_intervals: Annotated[list[tuple[str, int, int, str]],
        "List of polyadenylation (pA) site intervals as tuples (chromosome, start, end, strand). "
        "Example: [('chr22', 36189128, 36189129, '-'), ('chr22', 36190089, 36190090, '-')]. "
        "Each tuple defines a 1bp interval marking a pA site location."] = [
            ("chr22", 36189128, 36189129, "-"),
            ("chr22", 36190089, 36190090, "-"),
            ("chr22", 36190144, 36190145, "-")
        ],
    pa_labels: Annotated[list[str],
        "Labels for each pA site annotation. "
        "Example: ['pA_3', 'pA_2', 'pA_1']. "
        "Must have same length as pa_site_intervals."] = ["pA_3", "pA_2", "pA_1"],
    plot_offset: Annotated[int,
        "Number of base pairs to add on either side of pA site region for plotting. "
        "Typical values: 100-500."] = 200,
    chromosome: Annotated[str,
        "Chromosome identifier for the prediction interval."] = "chr22",
    start_position: Annotated[int,
        "Start position of the prediction interval (0-based)."] = 36150498,
    end_position: Annotated[int,
        "End position of the prediction interval (0-based)."] = 36252898,
    ontology_terms: Annotated[list[str],
        "List of tissue ontology IDs for predictions. "
        "Example: ['UBERON:0001159', 'UBERON:0002048'] for colon and lung."] = ["UBERON:0001159", "UBERON:0002048"],
    out_prefix: Annotated[str | None,
        "Prefix for output file names. If None, uses 'pa_sites_{timestamp}'."] = None,
) -> dict:
    """
    Visualize custom polyadenylation site annotations with RNA-seq predictions.
    Input is pA site coordinates and tissue terms, output is annotated expression visualization.
    """
    # Validate inputs
    if len(pa_site_intervals) != len(pa_labels):
        raise ValueError("pa_site_intervals and pa_labels must have the same length")
    
    # Set output prefix
    if out_prefix is None:
        out_prefix = f"pa_sites_{timestamp}"
    
    # Create DNA model client
    dna_model = dna_client.create(api_key)
    
    # Load gene annotations
    gtf = pd.read_feather(
        'https://storage.googleapis.com/alphagenome/reference/gencode/'
        'hg38/gencode.v46.annotation.gtf.gz.feather'
    )
    
    gtf_transcript = gene_annotation.filter_transcript_support_level(
        gene_annotation.filter_protein_coding(gtf), ['1']
    )
    
    gtf_longest_transcript = gene_annotation.filter_to_longest_transcript(gtf_transcript)
    longest_transcript_extractor = transcript.TranscriptExtractor(gtf_longest_transcript)
    
    # Define prediction interval
    interval = genome.Interval(chromosome, start_position, end_position).resize(
        dna_client.SEQUENCE_LENGTH_1MB
    )
    
    # Make predictions
    output = dna_model.predict_interval(
        interval=interval,
        requested_outputs={
            dna_client.OutputType.RNA_SEQ,
        },
        ontology_terms=ontology_terms,
    )
    
    # Get longest transcripts
    longest_transcripts = longest_transcript_extractor.extract(interval)
    
    # Create pA site intervals from input tuples
    apol4_pAs = [genome.Interval(chrom, start, end, strand) 
                 for chrom, start, end, strand in pa_site_intervals]
    
    # Define plotting interval based on first and last pA site
    first_pa_start = min(start for _, start, _, _ in pa_site_intervals)
    last_pa_end = max(end for _, _, end, _ in pa_site_intervals)
    pA_interval = genome.Interval(
        chromosome,
        first_pa_start - plot_offset,
        last_pa_end + plot_offset,
        '-'
    )
    
    # Build plot
    plot = plot_components.plot(
        [
            plot_components.TranscriptAnnotation(longest_transcripts),
            plot_components.Tracks(
                tdata=output.rna_seq.filter_to_negative_strand(),
                ylabel_template='RNA_SEQ: {biosample_name} ({strand})\n{name}',
                shared_y_scale=True,
            )
        ],
        annotations=[
            plot_components.IntervalAnnotation(
                apol4_pAs,
                alpha=1,
                labels=pa_labels,
                label_angle=90
            )
        ],
        interval=pA_interval,
        title='APOL4 polyadenylation sites annotation',
    )
    
    # Save figure
    figure_path = OUTPUT_DIR / f"{out_prefix}_pa_sites.png"
    plt.savefig(figure_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return {
        "message": "Polyadenylation sites visualization completed successfully",
        "reference": "https://github.com/google-deepmind/alphagenome/blob/main/colabs/visualization_modality_tour.ipynb",
        "artifacts": [
            {
                "description": "pA sites annotation figure",
                "path": str(figure_path.resolve())
            }
        ]
    }


@visualization_modality_tour_mcp.tool()
def alphagenome_visualize_chromatin_accessibility(
    api_key: Annotated[str,
        "AlphaGenome API key for making predictions."],
    variant_string: Annotated[str,
        "Variant in format 'chr:pos:ref>alt' for annotation (e.g., 'chr22:36201698:A>C')."],
    chromosome: Annotated[str,
        "Chromosome identifier."] = "chr22",
    start_position: Annotated[int,
        "Start position of genomic interval (0-based)."] = 36150498,
    end_position: Annotated[int,
        "End position of genomic interval (0-based)."] = 36252898,
    window_size: Annotated[int,
        "Size in base pairs of the window around the variant to plot. "
        "Typical values: 5000-10000."] = 8000,
    promoter_intervals: Annotated[list[tuple[str, int, int, str]] | None,
        "Optional list of promoter intervals as tuples (chr, start, end, name). "
        "Example: [('chr22', 36201799, 36202681, 'Ensembl_promoter:ENSR00001367790')]. "
        "If None, no promoter annotations are shown."] = None,
    ontology_terms: Annotated[list[str],
        "List of intestinal tissue ontology IDs. "
        "Example: ['UBERON:0000317', 'UBERON:0001155'] for various intestinal tissues."] = [
            "UBERON:0000317",
            "UBERON:0001155",
            "UBERON:0001157",
            "UBERON:0001159",
            "UBERON:0004992",
            "UBERON:0008971"
        ],
    out_prefix: Annotated[str | None,
        "Prefix for output file names. If None, uses 'chromatin_accessibility_{timestamp}'."] = None,
) -> dict:
    """
    Visualize predicted chromatin accessibility (DNASE and ATAC) across intestinal tissues.
    Input is genomic interval and tissue terms, output is accessibility prediction visualization.
    """
    # Set output prefix
    if out_prefix is None:
        out_prefix = f"chromatin_accessibility_{timestamp}"
    
    # Create DNA model client
    dna_model = dna_client.create(api_key)
    
    # Load gene annotations
    gtf = pd.read_feather(
        'https://storage.googleapis.com/alphagenome/reference/gencode/'
        'hg38/gencode.v46.annotation.gtf.gz.feather'
    )
    
    gtf_transcript = gene_annotation.filter_transcript_support_level(
        gene_annotation.filter_protein_coding(gtf), ['1']
    )
    
    gtf_longest_transcript = gene_annotation.filter_to_longest_transcript(gtf_transcript)
    longest_transcript_extractor = transcript.TranscriptExtractor(gtf_longest_transcript)
    
    # Define interval
    interval = genome.Interval(chromosome, start_position, end_position).resize(
        dna_client.SEQUENCE_LENGTH_1MB
    )
    
    # Make predictions
    output = dna_model.predict_interval(
        interval,
        requested_outputs={
            dna_client.OutputType.DNASE,
            dna_client.OutputType.ATAC,
        },
        ontology_terms=ontology_terms,
    )
    
    # Get longest transcripts
    longest_transcripts = longest_transcript_extractor.extract(interval)
    
    # Parse variant
    variant = genome.Variant.from_str(variant_string)
    
    # Prepare annotations
    annotations = [plot_components.VariantAnnotation([variant])]
    
    if promoter_intervals is not None:
        promoter_objs = [
            genome.Interval(chrom, start, end, name=name)
            for chrom, start, end, name in promoter_intervals
        ]
        annotations.append(plot_components.IntervalAnnotation(promoter_objs))
    
    # Build plot
    plot = plot_components.plot(
        [
            plot_components.TranscriptAnnotation(longest_transcripts),
            plot_components.Tracks(
                tdata=output.dnase,
                ylabel_template='DNASE: {biosample_name} ({strand})\n{name}',
            ),
            plot_components.Tracks(
                tdata=output.atac,
                ylabel_template='ATAC: {biosample_name} ({strand})\n{name}',
            ),
        ],
        interval=variant.reference_interval.resize(window_size),
        annotations=annotations,
        title='Predicted chromatin accessibility (DNASE, ATAC) for colon tissue',
    )
    
    # Save figure
    figure_path = OUTPUT_DIR / f"{out_prefix}_chromatin_accessibility.png"
    plt.savefig(figure_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return {
        "message": "Chromatin accessibility visualization completed successfully",
        "reference": "https://github.com/google-deepmind/alphagenome/blob/main/colabs/visualization_modality_tour.ipynb",
        "artifacts": [
            {
                "description": "Chromatin accessibility figure",
                "path": str(figure_path.resolve())
            }
        ]
    }


@visualization_modality_tour_mcp.tool()
def alphagenome_visualize_splicing(
    api_key: Annotated[str,
        "AlphaGenome API key for making predictions."],
    variant_string: Annotated[str,
        "Variant in format 'chr:pos:ref>alt' (e.g., 'chr22:36201698:A>C')."],
    chromosome: Annotated[str,
        "Chromosome identifier."] = "chr22",
    start_position: Annotated[int,
        "Start position of genomic interval (0-based)."] = 36150498,
    end_position: Annotated[int,
        "End position of genomic interval (0-based)."] = 36252898,
    gene_symbol: Annotated[str,
        "Gene symbol to zoom into for visualization (e.g., 'APOL4')."] = "APOL4",
    zoom_offset: Annotated[int,
        "Base pairs to add on either side of gene for visualization."] = 1000,
    tissue_filter: Annotated[str,
        "Tissue name to filter splice junctions for sashimi plot (e.g., 'Colon_Transverse')."] = "Colon_Transverse",
    ontology_terms: Annotated[list[str],
        "List of tissue ontology IDs. "
        "Example: ['UBERON:0001157', 'UBERON:0001159'] for intestinal tissues."] = [
            "UBERON:0001157",
            "UBERON:0001159"
        ],
    out_prefix: Annotated[str | None,
        "Prefix for output file names. If None, uses 'splicing_{timestamp}'."] = None,
) -> dict:
    """
    Visualize predicted splicing effects with sashimi plots showing REF vs ALT splice junctions.
    Input is variant and genomic interval, output is comprehensive splicing visualization.
    """
    # Set output prefix
    if out_prefix is None:
        out_prefix = f"splicing_{timestamp}"
    
    # Create DNA model client
    dna_model = dna_client.create(api_key)
    
    # Load gene annotations
    gtf = pd.read_feather(
        'https://storage.googleapis.com/alphagenome/reference/gencode/'
        'hg38/gencode.v46.annotation.gtf.gz.feather'
    )
    
    gtf_transcript = gene_annotation.filter_transcript_support_level(
        gene_annotation.filter_protein_coding(gtf), ['1']
    )
    
    # Create extractors for both all transcripts and longest
    transcript_extractor = transcript.TranscriptExtractor(gtf_transcript)
    
    # Define interval
    interval = genome.Interval(chromosome, start_position, end_position).resize(
        dna_client.SEQUENCE_LENGTH_1MB
    )
    
    # Parse variant
    variant = genome.Variant.from_str(variant_string)
    
    # Make predictions for REF and ALT alleles
    output = dna_model.predict_variant(
        interval=interval,
        variant=variant,
        requested_outputs={
            dna_client.OutputType.RNA_SEQ,
            dna_client.OutputType.SPLICE_SITES,
            dna_client.OutputType.SPLICE_SITE_USAGE,
            dna_client.OutputType.SPLICE_JUNCTIONS,
        },
        ontology_terms=ontology_terms,
    )
    
    # Get all transcripts (not just longest)
    transcripts = transcript_extractor.extract(interval)
    
    # Zoom in on the gene region
    gene_interval = gene_annotation.get_gene_interval(gtf, gene_symbol=gene_symbol)
    gene_interval.resize_inplace(gene_interval.width + zoom_offset)
    
    # Extract REF and ALT outputs
    ref_output = output.reference
    alt_output = output.alternate
    
    # Define colors
    ref_alt_colors = {'REF': 'dimgrey', 'ALT': 'red'}
    
    # Build plot
    plot = plot_components.plot(
        [
            plot_components.TranscriptAnnotation(transcripts),
            plot_components.Sashimi(
                ref_output.splice_junctions
                .filter_to_strand('-')
                .filter_by_tissue(tissue_filter),
                ylabel_template='Reference {biosample_name} ({strand})\n{name}',
            ),
            plot_components.Sashimi(
                alt_output.splice_junctions
                .filter_to_strand('-')
                .filter_by_tissue(tissue_filter),
                ylabel_template='Alternate {biosample_name} ({strand})\n{name}',
            ),
            plot_components.OverlaidTracks(
                tdata={
                    'REF': ref_output.rna_seq.filter_to_nonpositive_strand(),
                    'ALT': alt_output.rna_seq.filter_to_nonpositive_strand(),
                },
                colors=ref_alt_colors,
                ylabel_template='RNA_SEQ: {biosample_name} ({strand})\n{name}',
            ),
            plot_components.OverlaidTracks(
                tdata={
                    'REF': ref_output.splice_sites.filter_to_nonpositive_strand(),
                    'ALT': alt_output.splice_sites.filter_to_nonpositive_strand(),
                },
                colors=ref_alt_colors,
                ylabel_template='SPLICE SITES: {name} ({strand})',
            ),
            plot_components.OverlaidTracks(
                tdata={
                    'REF': (
                        ref_output.splice_site_usage.filter_to_nonpositive_strand()
                    ),
                    'ALT': (
                        alt_output.splice_site_usage.filter_to_nonpositive_strand()
                    ),
                },
                colors=ref_alt_colors,
                ylabel_template=(
                    'SPLICE SITE USAGE: {biosample_name} ({strand})\n{name}'
                ),
            ),
        ],
        interval=gene_interval,
        annotations=[plot_components.VariantAnnotation([variant])],
        title='Predicted REF vs. ALT effects of variant in colon tissue',
    )
    
    # Save figure
    figure_path = OUTPUT_DIR / f"{out_prefix}_splicing.png"
    plt.savefig(figure_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return {
        "message": "Splicing visualization completed successfully",
        "reference": "https://github.com/google-deepmind/alphagenome/blob/main/colabs/visualization_modality_tour.ipynb",
        "artifacts": [
            {
                "description": "Splicing effects figure",
                "path": str(figure_path.resolve())
            }
        ]
    }


@visualization_modality_tour_mcp.tool()
def alphagenome_visualize_chip_histone(
    api_key: Annotated[str,
        "AlphaGenome API key for making predictions."],
    chromosome: Annotated[str,
        "Chromosome identifier."] = "chr22",
    start_position: Annotated[int,
        "Start position of genomic interval (0-based)."] = 36150498,
    end_position: Annotated[int,
        "End position of genomic interval (0-based)."] = 36252898,
    show_tss: Annotated[bool,
        "Whether to show transcription start site annotations as blue intervals."] = True,
    ontology_terms: Annotated[list[str],
        "List of tissue ontology IDs for colon tissues. "
        "Example: ['UBERON:0000317', 'UBERON:0001155'] for various colon regions."] = [
            "UBERON:0000317",
            "UBERON:0001155",
            "UBERON:0001157",
            "UBERON:0001159"
        ],
    out_prefix: Annotated[str | None,
        "Prefix for output file names. If None, uses 'chip_histone_{timestamp}'."] = None,
) -> dict:
    """
    Visualize predicted histone modification markers (H3K4me3, H3K27ac, etc.) across colon tissues.
    Input is genomic interval and tissue terms, output is color-coded histone modification visualization.
    """
    # Set output prefix
    if out_prefix is None:
        out_prefix = f"chip_histone_{timestamp}"
    
    # Create DNA model client
    dna_model = dna_client.create(api_key)
    
    # Load gene annotations
    gtf = pd.read_feather(
        'https://storage.googleapis.com/alphagenome/reference/gencode/'
        'hg38/gencode.v46.annotation.gtf.gz.feather'
    )
    
    gtf_transcript = gene_annotation.filter_transcript_support_level(
        gene_annotation.filter_protein_coding(gtf), ['1']
    )
    
    gtf_longest_transcript = gene_annotation.filter_to_longest_transcript(gtf_transcript)
    longest_transcript_extractor = transcript.TranscriptExtractor(gtf_longest_transcript)
    
    # Define interval
    interval = genome.Interval(chromosome, start_position, end_position).resize(
        dna_client.SEQUENCE_LENGTH_1MB
    )
    
    # Make predictions
    output = dna_model.predict_interval(
        interval=interval,
        requested_outputs={dna_client.OutputType.CHIP_HISTONE},
        ontology_terms=ontology_terms,
    )
    
    # Get longest transcripts
    longest_transcripts = longest_transcript_extractor.extract(interval)
    
    # Reorder tracks by histone mark
    reordered_chip_histone = output.chip_histone.select_tracks_by_index(
        output.chip_histone.metadata.sort_values('histone_mark').index
    )
    
    # Define histone mark colors
    histone_to_color = {
        'H3K27AC': '#e41a1c',
        'H3K36ME3': '#ff7f00',
        'H3K4ME1': '#377eb8',
        'H3K4ME3': '#984ea3',
        'H3K9AC': '#4daf4a',
        'H3K27ME3': '#ffc0cb',
    }
    
    track_colors = (
        reordered_chip_histone.metadata['histone_mark']
        .map(lambda x: histone_to_color.get(x.upper(), '#000000'))
        .values
    )
    
    # Prepare annotations
    annotations = []
    if show_tss:
        gtf_tss = gene_annotation.extract_tss(gtf_longest_transcript)
        tss_as_intervals = [
            genome.Interval(
                chromosome=row.Chromosome,
                start=row.Start,
                end=row.End + 1000,  # Add extra 1Kb so the TSSs are visible
                name=row.gene_name,
            )
            for _, row in gtf_tss.iterrows()
        ]
        annotations.append(
            plot_components.IntervalAnnotation(
                tss_as_intervals, alpha=0.5, colors='blue'
            )
        )
    
    # Build plot
    plot = plot_components.plot(
        [
            plot_components.TranscriptAnnotation(longest_transcripts),
            plot_components.Tracks(
                tdata=reordered_chip_histone,
                ylabel_template=(
                    'CHIP HISTONE: {biosample_name} ({strand})\n{histone_mark}'
                ),
                filled=True,
                track_colors=track_colors,
            ),
        ],
        interval=interval,
        annotations=annotations,
        despine_keep_bottom=True,
        title='Predicted histone modification markers in colon tissue',
    )
    
    # Save figure
    figure_path = OUTPUT_DIR / f"{out_prefix}_chip_histone.png"
    plt.savefig(figure_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return {
        "message": "ChIP-Histone visualization completed successfully",
        "reference": "https://github.com/google-deepmind/alphagenome/blob/main/colabs/visualization_modality_tour.ipynb",
        "artifacts": [
            {
                "description": "Histone modification markers figure",
                "path": str(figure_path.resolve())
            }
        ]
    }


@visualization_modality_tour_mcp.tool()
def alphagenome_visualize_chip_tf(
    api_key: Annotated[str,
        "AlphaGenome API key for making predictions."],
    chromosome: Annotated[str,
        "Chromosome identifier."] = "chr22",
    start_position: Annotated[int,
        "Start position of genomic interval (0-based)."] = 36150498,
    end_position: Annotated[int,
        "End position of genomic interval (0-based)."] = 36252898,
    gene_symbol: Annotated[str | None,
        "Gene symbol to focus visualization on (e.g., 'APOL4'). "
        "If provided, will create a second zoomed-in plot around this gene."] = None,
    zoom_offset: Annotated[int,
        "Base pairs to add on either side of gene for zoomed visualization."] = 1000,
    filter_threshold: Annotated[float,
        "Minimum maximum prediction value for filtering tracks. "
        "Tracks with max value below this are excluded. Typical: 5000-10000."] = 8000,
    top_n_tracks: Annotated[int,
        "Number of top tracks to show in gene-focused plot (ranked by max prediction)."] = 10,
    show_tss: Annotated[bool,
        "Whether to show transcription start site annotations."] = True,
    ontology_terms: Annotated[list[str],
        "List of tissue/cell-line ontology IDs. "
        "Example: ['UBERON:0001159', 'EFO:0002067'] for sigmoid colon and K562."] = [
            "UBERON:0001159",
            "UBERON:0001157",
            "EFO:0002067",
            "EFO:0001187"
        ],
    out_prefix: Annotated[str | None,
        "Prefix for output file names. If None, uses 'chip_tf_{timestamp}'."] = None,
) -> dict:
    """
    Visualize predicted transcription factor binding across cell types and tissues.
    Input is genomic interval and tissue terms, output includes genome-wide and gene-focused TF binding plots.
    """
    # Set output prefix
    if out_prefix is None:
        out_prefix = f"chip_tf_{timestamp}"
    
    # Create DNA model client
    dna_model = dna_client.create(api_key)
    
    # Load gene annotations
    gtf = pd.read_feather(
        'https://storage.googleapis.com/alphagenome/reference/gencode/'
        'hg38/gencode.v46.annotation.gtf.gz.feather'
    )
    
    gtf_transcript = gene_annotation.filter_transcript_support_level(
        gene_annotation.filter_protein_coding(gtf), ['1']
    )
    
    gtf_longest_transcript = gene_annotation.filter_to_longest_transcript(gtf_transcript)
    longest_transcript_extractor = transcript.TranscriptExtractor(gtf_longest_transcript)
    transcript_extractor = transcript.TranscriptExtractor(gtf_transcript)
    
    # Define interval
    interval = genome.Interval(chromosome, start_position, end_position).resize(
        dna_client.SEQUENCE_LENGTH_1MB
    )
    
    # Make predictions
    output = dna_model.predict_interval(
        interval=interval,
        requested_outputs={dna_client.OutputType.CHIP_TF},
        ontology_terms=ontology_terms,
    )
    
    # Get longest transcripts
    longest_transcripts = longest_transcript_extractor.extract(interval)
    
    # Filter tracks by threshold
    output_filtered = output.chip_tf.filter_tracks(
        output.chip_tf.values.max(axis=0) > filter_threshold
    )
    
    # Prepare TSS annotations if requested
    annotations = []
    if show_tss:
        gtf_tss = gene_annotation.extract_tss(gtf_longest_transcript)
        tss_as_intervals = [
            genome.Interval(
                chromosome=row.Chromosome,
                start=row.Start,
                end=row.End + 1000,
                name=row.gene_name,
            )
            for _, row in gtf_tss.iterrows()
        ]
        annotations.append(
            plot_components.IntervalAnnotation(
                tss_as_intervals, alpha=0.3, colors='blue'
            )
        )
    
    # Build main plot
    plot_components.plot(
        components=[
            plot_components.TranscriptAnnotation(longest_transcripts),
            plot_components.Tracks(
                tdata=output_filtered,
                ylabel_template=(
                    'CHIP TF: {biosample_name} ({strand})\n{transcription_factor}'
                ),
                filled=True,
            ),
        ],
        interval=interval,
        title='Predicted TF-binding in K562 and HepG2 cell-lines.',
        despine_keep_bottom=True,
        annotations=annotations,
    )
    
    # Save main figure
    figure_path_main = OUTPUT_DIR / f"{out_prefix}_chip_tf_genome_wide.png"
    plt.savefig(figure_path_main, dpi=300, bbox_inches='tight')
    plt.close()
    
    artifacts = [
        {
            "description": "Genome-wide TF binding figure",
            "path": str(figure_path_main.resolve())
        }
    ]
    
    # Create gene-focused plot if gene symbol provided
    if gene_symbol is not None:
        # Get all transcripts
        transcripts = transcript_extractor.extract(interval)
        
        # Zoom in on the gene region
        gene_interval = gene_annotation.get_gene_interval(gtf, gene_symbol=gene_symbol)
        gene_interval.resize_inplace(gene_interval.width + zoom_offset)
        
        # Compute max predicted values per track in the gene interval
        max_predictions = output.chip_tf.slice_by_interval(
            gene_interval, match_resolution=True
        ).values.max(axis=0)
        
        # Filter to top N tracks
        output_gene_filtered = output.chip_tf.filter_tracks(
            (max_predictions >= np.sort(max_predictions)[-top_n_tracks])
        )
        
        # Build gene-focused plot
        plot_components.plot(
            [
                plot_components.TranscriptAnnotation(transcripts),
                plot_components.Tracks(
                    tdata=output_gene_filtered,
                    ylabel_template=(
                        'CHIP TF: {biosample_name} ({strand})\n{transcription_factor}'
                    ),
                    filled=True,
                ),
            ],
            interval=gene_interval,
            annotations=[],
            despine_keep_bottom=True,
            title=f'Predicted TF-binding around {gene_symbol}',
        )
        
        # Save gene-focused figure
        figure_path_gene = OUTPUT_DIR / f"{out_prefix}_chip_tf_gene_focused.png"
        plt.savefig(figure_path_gene, dpi=300, bbox_inches='tight')
        plt.close()
        
        artifacts.append({
            "description": f"Gene-focused TF binding figure ({gene_symbol})",
            "path": str(figure_path_gene.resolve())
        })
    
    return {
        "message": "ChIP-TF visualization completed successfully",
        "reference": "https://github.com/google-deepmind/alphagenome/blob/main/colabs/visualization_modality_tour.ipynb",
        "artifacts": artifacts
    }


@visualization_modality_tour_mcp.tool()
def alphagenome_visualize_contact_maps(
    api_key: Annotated[str,
        "AlphaGenome API key for making predictions."],
    chromosome: Annotated[str,
        "Chromosome identifier."] = "chr22",
    start_position: Annotated[int,
        "Start position of genomic interval (0-based)."] = 36150498,
    end_position: Annotated[int,
        "End position of genomic interval (0-based)."] = 36252898,
    colormap: Annotated[str,
        "Matplotlib colormap name for contact map visualization. "
        "Common choices: 'autumn_r', 'Reds', 'YlOrRd'."] = "autumn_r",
    vmax: Annotated[float,
        "Maximum value for colormap scaling. Typical: 0.5-2.0."] = 1.0,
    ontology_terms: Annotated[list[str],
        "List of cell-line ontology IDs. "
        "Example: ['EFO:0002824'] for HCT116 colon carcinoma cell line."] = ["EFO:0002824"],
    out_prefix: Annotated[str | None,
        "Prefix for output file names. If None, uses 'contact_maps_{timestamp}'."] = None,
) -> dict:
    """
    Visualize predicted DNA-DNA contact maps showing chromatin interaction patterns.
    Input is genomic interval and cell-line terms, output is contact map heatmap visualization.
    """
    # Set output prefix
    if out_prefix is None:
        out_prefix = f"contact_maps_{timestamp}"
    
    # Create DNA model client
    dna_model = dna_client.create(api_key)
    
    # Load gene annotations
    gtf = pd.read_feather(
        'https://storage.googleapis.com/alphagenome/reference/gencode/'
        'hg38/gencode.v46.annotation.gtf.gz.feather'
    )
    
    gtf_transcript = gene_annotation.filter_transcript_support_level(
        gene_annotation.filter_protein_coding(gtf), ['1']
    )
    
    gtf_longest_transcript = gene_annotation.filter_to_longest_transcript(gtf_transcript)
    longest_transcript_extractor = transcript.TranscriptExtractor(gtf_longest_transcript)
    
    # Define interval
    interval = genome.Interval(chromosome, start_position, end_position).resize(
        dna_client.SEQUENCE_LENGTH_1MB
    )
    
    # Make predictions
    output = dna_model.predict_interval(
        interval=interval,
        requested_outputs={dna_client.OutputType.CONTACT_MAPS},
        ontology_terms=ontology_terms,
    )
    
    # Get longest transcripts
    longest_transcripts = longest_transcript_extractor.extract(interval)
    
    # Build plot
    plot = plot_components.plot(
        [
            plot_components.TranscriptAnnotation(longest_transcripts),
            plot_components.ContactMaps(
                tdata=output.contact_maps,
                ylabel_template='{biosample_name}\n{name}',
                cmap=colormap,
                vmax=vmax,
            ),
        ],
        interval=interval,
        title='Predicted contact maps',
    )
    
    # Save figure
    figure_path = OUTPUT_DIR / f"{out_prefix}_contact_maps.png"
    plt.savefig(figure_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return {
        "message": "Contact maps visualization completed successfully",
        "reference": "https://github.com/google-deepmind/alphagenome/blob/main/colabs/visualization_modality_tour.ipynb",
        "artifacts": [
            {
                "description": "DNA-DNA contact maps figure",
                "path": str(figure_path.resolve())
            }
        ]
    }
