"""
Variant scoring and visualization tools for AlphaGenome.

This tutorial demonstrates how to score and visualize the effects of genetic variants
on different genomic modalities including gene expression and chromatin accessibility.

This MCP Server provides 2 tools:
1. alphagenome_score_variant: Score a single genetic variant across multiple genomic features
2. alphagenome_visualize_variant_effects: Visualize predicted variant effects across genomic modalities

All tools extracted from `google-deepmind/alphagenome/blob/main/colabs/variant_scoring_ui.ipynb`.
"""

# Standard imports
from typing import Annotated, Literal, Any
import pandas as pd
from pathlib import Path
import os
from fastmcp import FastMCP
from datetime import datetime
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

# AlphaGenome imports
from alphagenome.data import gene_annotation, genome, transcript, track_data
from alphagenome.models import dna_client, variant_scorers
from alphagenome.visualization import plot_components

# Project structure
PROJECT_ROOT = Path(__file__).parent.parent.parent.resolve()
DEFAULT_INPUT_DIR = PROJECT_ROOT / "tmp" / "inputs"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "tmp" / "outputs"

INPUT_DIR = Path(os.environ.get("VARIANT_SCORING_UI_INPUT_DIR", DEFAULT_INPUT_DIR))
OUTPUT_DIR = Path(os.environ.get("VARIANT_SCORING_UI_OUTPUT_DIR", DEFAULT_OUTPUT_DIR))

# Ensure directories exist
INPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Timestamp for unique outputs
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

# Reference genome annotations
HG38_GTF_FEATHER = (
    'https://storage.googleapis.com/alphagenome/reference/gencode/'
    'hg38/gencode.v46.annotation.gtf.gz.feather'
)
MM10_GTF_FEATHER = (
    'https://storage.googleapis.com/alphagenome/reference/gencode/'
    'mm10/gencode.vM23.annotation.gtf.gz.feather'
)

# Global caches for efficiency
_prediction_cache = {}
_transcript_extractor_cache = {}

# MCP server instance
variant_scoring_ui_mcp = FastMCP(name="variant_scoring_ui")

# Matplotlib configuration for high-resolution figures
plt.rcParams["figure.dpi"] = 300
plt.rcParams["savefig.dpi"] = 300


@variant_scoring_ui_mcp.tool()
def alphagenome_score_variant(
    variant_chromosome: Annotated[str,
        "Chromosome identifier where the variant is located. "
        "Format: 'chr1', 'chr2', ..., 'chrX', 'chrY'. "
        "Example: 'chr22' for chromosome 22."],
    
    variant_position: Annotated[int,
        "Genomic position of the variant (1-based coordinate). "
        "This is the position on the chromosome where the variant occurs. "
        "Example: 36201698 for position 36,201,698 bp on the chromosome."],
    
    variant_reference_bases: Annotated[str,
        "Reference allele base(s) at this position in the reference genome. "
        "Single nucleotide: 'A', 'C', 'G', or 'T'. "
        "Multi-nucleotide variants: e.g., 'AT', 'GCC'. "
        "Example: 'A' for adenine reference allele."],
    
    variant_alternate_bases: Annotated[str,
        "Alternate allele base(s) representing the variant. "
        "Single nucleotide: 'A', 'C', 'G', or 'T'. "
        "Multi-nucleotide variants: e.g., 'AT', 'GCC'. "
        "Example: 'C' for cytosine alternate allele."],
    
    api_key: Annotated[str,
        "AlphaGenome API key for authentication. "
        "Required to access the AlphaGenome model prediction service. "
        "Obtain from Google Cloud Console with AlphaGenome API enabled."],
    
    organism: Annotated[Literal["human", "mouse"],
        "Organism for which to score the variant. "
        "'human': Homo sapiens (uses hg38 reference genome). "
        "'mouse': Mus musculus (uses mm10 reference genome)."] = "human",
    
    sequence_length: Annotated[Literal["2KB", "16KB", "100KB", "500KB", "1MB"],
        "Length of genomic sequence context around the variant to analyze. "
        "Larger contexts capture long-range regulatory effects but take longer to compute. "
        "'2KB': 2,048 bp - local effects only. "
        "'16KB': 16,384 bp - nearby regulatory elements. "
        "'100KB': 102,400 bp - moderate-range interactions. "
        "'500KB': 512,000 bp - long-range regulatory domains. "
        "'1MB': 1,048,576 bp - comprehensive genomic context (tutorial default)."] = "1MB",
    
    out_prefix: Annotated[str | None,
        "Prefix for output file paths. "
        "If None, uses 'variant_scores_{timestamp}' format. "
        "Output file: '{out_prefix}_scores.csv' containing variant effect scores."] = None,
) -> dict:
    """
    Score a single genetic variant across thousands of genomic features using AlphaGenome.
    Input is variant details (chromosome, position, ref/alt alleles) and output is a CSV file with effect scores across 38,357 genomic features including gene expression, chromatin accessibility, and histone modifications.
    """
    # Input validation
    if not api_key:
        raise ValueError("AlphaGenome API key must be provided")
    
    # Map organism to internal representation
    organism_map = {
        'human': dna_client.Organism.HOMO_SAPIENS,
        'mouse': dna_client.Organism.MUS_MUSCULUS,
    }
    organism_enum = organism_map[organism]
    
    # Initialize DNA model
    dna_model = dna_client.create(api_key)
    
    # Create variant object
    variant = genome.Variant(
        chromosome=variant_chromosome,
        position=variant_position,
        reference_bases=variant_reference_bases,
        alternate_bases=variant_alternate_bases,
    )
    
    # Convert sequence length to internal format
    sequence_length_enum = dna_client.SUPPORTED_SEQUENCE_LENGTHS[
        f'SEQUENCE_LENGTH_{sequence_length}'
    ]
    
    # Create genomic interval centered on variant
    interval = variant.reference_interval.resize(sequence_length_enum)
    
    # Score variant across all recommended scorers
    variant_scores = dna_model.score_variant(
        interval=interval,
        variant=variant,
        variant_scorers=list(variant_scorers.RECOMMENDED_VARIANT_SCORERS.values()),
    )
    
    # Convert to tidy DataFrame format
    df_scores = variant_scorers.tidy_scores(variant_scores)
    
    # Prepare output path
    if out_prefix is None:
        out_prefix = f"variant_scores_{timestamp}"
    
    output_file = OUTPUT_DIR / f"{out_prefix}_scores.csv"
    
    # Save results
    # Remove columns that are not useful for users
    columns = [
        c for c in df_scores.columns if c not in ['variant_id', 'scored_interval']
    ]
    df_scores[columns].to_csv(output_file, index=False)
    
    return {
        "message": f"Scored variant {variant} across {len(df_scores)} genomic features",
        "reference": "https://github.com/google-deepmind/alphagenome/blob/main/colabs/variant_scoring_ui.ipynb",
        "artifacts": [
            {
                "description": "Variant effect scores across genomic features",
                "path": str(output_file.resolve())
            }
        ]
    }


@variant_scoring_ui_mcp.tool()
def alphagenome_visualize_variant_effects(
    variant_chromosome: Annotated[str,
        "Chromosome identifier where the variant is located. "
        "Format: 'chr1', 'chr2', ..., 'chrX', 'chrY'. "
        "Example: 'chr22' for chromosome 22."],
    
    variant_position: Annotated[int,
        "Genomic position of the variant (1-based coordinate). "
        "This is the position on the chromosome where the variant occurs. "
        "Example: 36201698 for position 36,201,698 bp on the chromosome."],
    
    variant_reference_bases: Annotated[str,
        "Reference allele base(s) at this position in the reference genome. "
        "Single nucleotide: 'A', 'C', 'G', or 'T'. "
        "Example: 'A' for adenine reference allele."],
    
    variant_alternate_bases: Annotated[str,
        "Alternate allele base(s) representing the variant. "
        "Single nucleotide: 'A', 'C', 'G', or 'T'. "
        "Example: 'C' for cytosine alternate allele."],
    
    api_key: Annotated[str,
        "AlphaGenome API key for authentication. "
        "Required to access the AlphaGenome model prediction service. "
        "Obtain from Google Cloud Console with AlphaGenome API enabled."],
    
    organism: Annotated[Literal["human", "mouse"],
        "Organism for which to visualize variant effects. "
        "'human': Homo sapiens (uses hg38 reference genome). "
        "'mouse': Mus musculus (uses mm10 reference genome)."] = "human",
    
    sequence_length: Annotated[Literal["2KB", "16KB", "100KB", "500KB", "1MB"],
        "Length of genomic sequence context around the variant to analyze. "
        "Larger contexts capture long-range regulatory effects. "
        "'1MB': 1,048,576 bp - comprehensive genomic context (tutorial default)."] = "1MB",
    
    ontology_terms: Annotated[list[str],
        "List of cell/tissue ontology identifiers to visualize variant effects in specific contexts. "
        "Format: EFO (Experimental Factor Ontology) or CL (Cell Ontology) identifiers. "
        "Examples: ['EFO:0001187'] for brain, ['EFO:0002067'] for heart, ['EFO:0002784'] for liver. "
        "Tutorial uses: ['EFO:0001187', 'EFO:0002067', 'EFO:0002784'] for brain, heart, and liver."] = ['EFO:0001187', 'EFO:0002067', 'EFO:0002784'],
    
    plot_gene_annotation: Annotated[bool,
        "Whether to include gene and transcript annotations in the visualization. "
        "Shows exon-intron structure and gene positions."] = True,
    
    plot_longest_transcript_only: Annotated[bool,
        "Whether to show only the longest transcript per gene (True) or all transcripts (False). "
        "True simplifies visualization, False shows full transcript diversity."] = True,
    
    plot_rna_seq: Annotated[bool,
        "Whether to plot RNA-seq expression predictions for reference and alternate alleles."] = True,
    
    plot_cage: Annotated[bool,
        "Whether to plot CAGE (Cap Analysis of Gene Expression) predictions showing transcription start sites."] = True,
    
    plot_atac: Annotated[bool,
        "Whether to plot ATAC-seq (chromatin accessibility) predictions."] = False,
    
    plot_dnase: Annotated[bool,
        "Whether to plot DNase-seq (chromatin accessibility) predictions."] = False,
    
    plot_chip_histone: Annotated[bool,
        "Whether to plot ChIP-seq histone modification predictions."] = False,
    
    plot_chip_tf: Annotated[bool,
        "Whether to plot ChIP-seq transcription factor binding predictions."] = False,
    
    plot_splice_sites: Annotated[bool,
        "Whether to plot splice site usage predictions."] = True,
    
    plot_splice_site_usage: Annotated[bool,
        "Whether to plot detailed splice site usage patterns."] = False,
    
    plot_contact_maps: Annotated[bool,
        "Whether to plot chromatin contact maps showing 3D genome organization."] = False,
    
    plot_splice_junctions: Annotated[bool,
        "Whether to plot splice junction predictions (sashimi plots)."] = False,
    
    filter_to_positive_strand: Annotated[bool,
        "Filter predictions to show only positive DNA strand features."] = False,
    
    filter_to_negative_strand: Annotated[bool,
        "Filter predictions to show only negative DNA strand features. "
        "Tutorial uses True to focus on negative strand features."] = True,
    
    transcription_factors: Annotated[list[str] | None,
        "List of specific transcription factors to visualize when plot_chip_tf=True. "
        "Must match 'transcription_factor' column values in ChIP-seq metadata. "
        "Example: ['IKZF1', 'CTCF']. If None, shows all available TFs."] = None,
    
    ref_color: Annotated[str,
        "Color for reference allele tracks in visualizations. "
        "Any matplotlib color name or hex code."] = 'dimgrey',
    
    alt_color: Annotated[str,
        "Color for alternate allele tracks in visualizations. "
        "Any matplotlib color name or hex code."] = 'red',
    
    plot_interval_width: Annotated[int,
        "Width (in base pairs) of the genomic region to visualize. "
        "Must be less than sequence_length. Range: 2048-1048576 bp. "
        "Tutorial uses 43,008 bp to focus on local variant effects."] = 43008,
    
    plot_interval_shift: Annotated[int,
        "Shift (in base pairs) of the plot interval from the variant position. "
        "Negative values shift left, positive shift right. Range: -524288 to 524288 bp. "
        "Tutorial uses 0 to center plot on variant."] = 0,
    
    out_prefix: Annotated[str | None,
        "Prefix for output file paths. "
        "If None, uses 'variant_effects_{timestamp}' format. "
        "Output: '{out_prefix}_visualization.png' containing the multi-panel plot."] = None,
) -> dict:
    """
    Visualize predicted effects of a genetic variant across multiple genomic modalities.
    Input is variant details and visualization options, output is a multi-panel figure showing variant effects on gene expression, chromatin accessibility, splice sites, and other genomic features for specified cell/tissue types.
    """
    # Input validation
    if not api_key:
        raise ValueError("AlphaGenome API key must be provided")
    
    if filter_to_positive_strand and filter_to_negative_strand:
        raise ValueError(
            'Cannot specify both filter_to_positive_strand and '
            'filter_to_negative_strand.'
        )
    
    # Map organism to internal representation
    organism_map = {
        'human': dna_client.Organism.HOMO_SAPIENS,
        'mouse': dna_client.Organism.MUS_MUSCULUS,
    }
    organism_enum = organism_map[organism]
    
    # Initialize DNA model
    dna_model = dna_client.create(api_key)
    
    # Create variant object
    variant = genome.Variant(
        chromosome=variant_chromosome,
        position=variant_position,
        reference_bases=variant_reference_bases,
        alternate_bases=variant_alternate_bases,
    )
    
    # Convert sequence length to internal format
    sequence_length_enum = dna_client.SUPPORTED_SEQUENCE_LENGTHS[
        f'SEQUENCE_LENGTH_{sequence_length}'
    ]
    
    # Create genomic interval centered on variant
    interval = variant.reference_interval.resize(sequence_length_enum)
    
    # Load gene annotation (with caching)
    if organism_enum in _transcript_extractor_cache:
        transcript_extractor, longest_transcript_extractor = (
            _transcript_extractor_cache[organism_enum]
        )
    else:
        match organism_enum:
            case dna_client.Organism.HOMO_SAPIENS:
                gtf_path = HG38_GTF_FEATHER
            case dna_client.Organism.MUS_MUSCULUS:
                gtf_path = MM10_GTF_FEATHER
            case _:
                raise ValueError(f'Unsupported organism: {organism_enum}')
        
        gtf = pd.read_feather(gtf_path)
        
        # Filter to protein-coding genes and highly supported transcripts
        gtf_transcript = gene_annotation.filter_transcript_support_level(
            gene_annotation.filter_protein_coding(gtf), ['1']
        )
        
        # Extractor for identifying transcripts in a region
        transcript_extractor = transcript.TranscriptExtractor(gtf_transcript)
        
        # Extractor for longest transcript per gene
        gtf_longest_transcript = gene_annotation.filter_to_longest_transcript(
            gtf_transcript
        )
        longest_transcript_extractor = transcript.TranscriptExtractor(
            gtf_longest_transcript
        )
        _transcript_extractor_cache[organism_enum] = (
            transcript_extractor,
            longest_transcript_extractor,
        )
    
    # Cached prediction function
    def _predict_variant_cached(
        interval, variant, organism, requested_outputs, ontology_terms
    ):
        """Cache wrapper of dna_model.predict_variant."""
        cache_key = (
            str(interval),
            str(variant),
            str(organism),
            tuple(requested_outputs),
            tuple(ontology_terms),
        )
        
        if cache_key in _prediction_cache:
            return _prediction_cache[cache_key]
        
        result = dna_model.predict_variant(
            interval=interval,
            variant=variant,
            organism=organism,
            requested_outputs=requested_outputs,
            ontology_terms=ontology_terms,
        )
        _prediction_cache[cache_key] = result
        return result
    
    # Predict variant effects
    output = _predict_variant_cached(
        interval=interval,
        variant=variant,
        organism=organism_enum,
        requested_outputs=[*dna_client.OutputType],
        ontology_terms=ontology_terms,
    )
    
    # Filter to DNA strand if requested
    ref, alt = output.reference, output.alternate
    
    if filter_to_positive_strand:
        ref = ref.filter_to_strand(strand='+')
        alt = alt.filter_to_strand(strand='+')
    elif filter_to_negative_strand:
        ref = ref.filter_to_strand(strand='-')
        alt = alt.filter_to_strand(strand='-')
    
    # Function to filter to TFs if requested
    def _maybe_filter_to_tfs(
        tdata: track_data.TrackData,
    ) -> track_data.TrackData | None:
        if not plot_chip_tf:
            return tdata
        
        if plot_chip_tf and (transcription_factors is not None):
            if not isinstance(transcription_factors, list):
                raise ValueError('TFs must be input as a list of strings.')
            if not isinstance(transcription_factors[0], str):
                raise ValueError('TFs must be input as a list of strings.')
        
        if transcription_factors is None:
            return tdata
        
        tf_rows = tdata.metadata.index[
            tdata.metadata['transcription_factor'].isin(transcription_factors)
        ]
        if not tf_rows.any():
            print(
                f'No tracks found for specified transcription factors and'
                f' ontology_terms.'
            )
            return None
        missing_tfs = set(transcription_factors).difference(
            set(tdata.metadata['transcription_factor'])
        )
        if missing_tfs:
            print(
                f'Could not find tracks in chip_tf outputs corresponding to the'
                f' following requested TFs'
            )
        return tdata.select_tracks_by_index(tf_rows)
    
    # Build plot components
    components = []
    ref_alt_colors = {'REF': ref_color, 'ALT': alt_color}
    
    # Gene and transcript annotation
    if plot_gene_annotation:
        if plot_longest_transcript_only:
            transcripts = longest_transcript_extractor.extract(interval)
        else:
            transcripts = transcript_extractor.extract(interval)
        components.append(plot_components.TranscriptAnnotation(transcripts))
    
    # Individual output type plots
    plot_map = {
        'plot_atac': (ref.atac, alt.atac, 'ATAC'),
        'plot_cage': (ref.cage, alt.cage, 'CAGE'),
        'plot_chip_histone': (ref.chip_histone, alt.chip_histone, 'CHIP_HISTONE'),
        'plot_chip_tf': (
            _maybe_filter_to_tfs(ref.chip_tf),
            _maybe_filter_to_tfs(alt.chip_tf),
            'CHIP_TF',
        ),
        'plot_contact_maps': (ref.contact_maps, alt.contact_maps, 'CONTACT_MAPS'),
        'plot_dnase': (ref.dnase, alt.dnase, 'DNASE'),
        'plot_rna_seq': (ref.rna_seq, alt.rna_seq, 'RNA_SEQ'),
        'plot_splice_junctions': (
            ref.splice_junctions,
            alt.splice_junctions,
            'SPLICE_JUNCTIONS',
        ),
        'plot_splice_sites': (ref.splice_sites, alt.splice_sites, 'SPLICE_SITES'),
        'plot_splice_site_usage': (
            ref.splice_site_usage,
            alt.splice_site_usage,
            'SPLICE_SITE_USAGE',
        ),
    }
    
    for key, (ref_data, alt_data, output_type) in plot_map.items():
        if eval(key) and ref_data is not None and ref_data.values.shape[-1] == 0:
            print(
                f'Requested plot for output {output_type} but no tracks exist in'
                ' output. This is likely because this output does not exist for your'
                ' ontologies or requested DNA strand.'
            )
        if eval(key) and ref_data and alt_data:
            match output_type:
                case 'CHIP_HISTONE':
                    ylabel_template = (
                        f'{output_type}: {{biosample_name}} ({{strand}})\n{{histone_mark}}'
                    )
                case 'CHIP_TF':
                    ylabel_template = (
                        f'{output_type}: {{biosample_name}}'
                        ' ({strand})\n{transcription_factor}'
                    )
                case 'CONTACT_MAPS':
                    ylabel_template = f'{output_type}: {{biosample_name}} ({{strand}})'
                case 'SPLICE_SITES':
                    ylabel_template = f'{output_type}: {{name}} ({{strand}})'
                case _:
                    ylabel_template = (
                        f'{output_type}: {{biosample_name}} ({{strand}})\n{{name}}'
                    )
            
            if output_type == 'CONTACT_MAPS':
                component = plot_components.ContactMapsDiff(
                    tdata=alt_data - ref_data,
                    ylabel_template=ylabel_template,
                )
                components.append(component)
            elif output_type == 'SPLICE_JUNCTIONS':
                ref_plot = plot_components.Sashimi(
                    ref_data,
                    ylabel_template='REF: ' + ylabel_template,
                )
                alt_plot = plot_components.Sashimi(
                    alt_data,
                    ylabel_template='ALT: ' + ylabel_template,
                )
                components.extend([ref_plot, alt_plot])
            else:
                component = plot_components.OverlaidTracks(
                    tdata={'REF': ref_data, 'ALT': alt_data},
                    colors=ref_alt_colors,
                    ylabel_template=ylabel_template,
                )
                components.append(component)
    
    # Validate plot interval width
    if plot_interval_width > interval.width:
        raise ValueError(
            f'plot_interval_width ({plot_interval_width}) must be less than '
            f'interval.width ({interval.width}).'
        )
    
    # Create the plot
    plot = plot_components.plot(
        components=components,
        interval=interval.shift(plot_interval_shift).resize(plot_interval_width),
        annotations=[
            plot_components.VariantAnnotation([variant]),
        ],
    )
    
    # Save figure
    if out_prefix is None:
        out_prefix = f"variant_effects_{timestamp}"
    
    output_file = OUTPUT_DIR / f"{out_prefix}_visualization.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    return {
        "message": f"Visualized effects of variant {variant} across {len(components)} genomic modalities",
        "reference": "https://github.com/google-deepmind/alphagenome/blob/main/colabs/variant_scoring_ui.ipynb",
        "artifacts": [
            {
                "description": "Multi-panel variant effects visualization",
                "path": str(output_file.resolve())
            }
        ]
    }
