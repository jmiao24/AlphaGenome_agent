"""
Example analysis workflow for TAL1 locus variant analysis.

This MCP Server provides 3 tools:
1. alphagenome_visualize_variant_positions: Visualizes genomic positions of variants near TAL1 gene
2. alphagenome_predict_variant_effects: Predicts functional impact of a variant on gene expression, accessibility, and histone marks
3. alphagenome_compare_variant_scores: Compares predicted TAL1 expression changes between cancer-associated and background variants

All tools extracted from `google-deepmind/alphagenome/blob/main/colabs/example_analysis_workflow.ipynb`.
"""

# Standard imports
from typing import Annotated, Literal, Any
import pandas as pd
import numpy as np
from pathlib import Path
import os
from fastmcp import FastMCP
from datetime import datetime
import io
import itertools
import matplotlib.pyplot as plt

# AlphaGenome imports
from alphagenome.data import gene_annotation
from alphagenome.data import genome
from alphagenome.data import transcript as transcript_utils
from alphagenome.models import dna_client
from alphagenome.models import variant_scorers
from alphagenome.visualization import plot_components
import plotnine as gg

# Project structure
PROJECT_ROOT = Path(__file__).parent.parent.parent.resolve()
DEFAULT_INPUT_DIR = PROJECT_ROOT / "tmp" / "inputs"
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "tmp" / "outputs"

INPUT_DIR = Path(os.environ.get("EXAMPLE_ANALYSIS_WORKFLOW_INPUT_DIR", DEFAULT_INPUT_DIR))
OUTPUT_DIR = Path(os.environ.get("EXAMPLE_ANALYSIS_WORKFLOW_OUTPUT_DIR", DEFAULT_OUTPUT_DIR))

# Ensure directories exist
INPUT_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Configure matplotlib for high-resolution figures
plt.rcParams["figure.dpi"] = 300
plt.rcParams["savefig.dpi"] = 300

# Timestamp for unique outputs
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

# MCP server instance
example_analysis_workflow_mcp = FastMCP(name="example_analysis_workflow")


# Utility functions from tutorial
def oncogenic_tal1_variants() -> pd.DataFrame:
    """Returns a dataframe of oncogenic T-ALL variants that affect TAL1."""
    variant_data = """
ID\tCHROM\tPOS\tREF\tALT\toutput\tStudy ID\tStudy Variant ID
Jurkat\tchr1\t47239296\tC\tCCGTTTCCTAACC\t1\tMansour_2014
MOLT-3\tchr1\t47239296\tC\tACC\t1\tMansour_2014
Patient_1\tchr1\t47239296\tC\tAACG\t1\tMansour_2014
Patient_2\tchr1\t47239291\tCTAACC\tTTTACCGTCTGTTAACGGC\t1\tMansour_2014
Patient_3-5\tchr1\t47239296\tC\tACG\t1\tMansour_2014
Patient_6\tchr1\t47239296\tC\tACC\t1\tMansour_2014
Patient_7\tchr1\t47239295\tAC\tTCAAACTGGTAACC\t1\tMansour_2014
Patient_8\tchr1\t47239296\tC\tAACC\t1\tMansour_2014
new 3' enhancer 1\tchr1\t47212072\tT\tTGGGTAAACCGTCTGTTCAGCG\t1\tSmith_2023\tUPNT802
new 3' enhancer 2\tchr1\t47212074\tG\tGAACGTT\t1\tSmith_2023\tUPNT613
intergenic SNV 1\tchr1\t47230639\tC\tT\t1\tLiu_2020\tSJALL043861_D1
intergenic SNV 2\tchr1\t47230639\tC\tT\t1\tLiu_2020\tSJALL018373_D1
SJALL040467_D1\tchr1\t47239296\tC\tAACC\t1\tLiu_2020\tSJALL040467_D1
PATBGC\tchr1\t47239296\tC\tAACC\t1\tLiu_2017\tPATBGC
PATBTX\tchr1\t47239296\tC\tACGGATATAACC\t1\tLiu_2017\tPATBTX
PARJAY\tchr1\t47239296\tC\tACGGAATTTCTAACC\t1\tLiu_2017\tPARJAY
PARSJG\tchr1\t47239296\tC\tAACC\t1\tLiu_2017\tPARSJG
PASYAJ\tchr1\t47239296\tC\tAACC\t1\tLiu_2017\tPASYAJ
PATRAB\tchr1\t47239293\tTTA\tCTAACGG\t1\tLiu_2017\tPATRAB
PAUBXP\tchr1\t47239296\tC\tACC\t1\tLiu_2017\tPAUBXP
PATENL\tchr1\t47239296\tC\tAACC\t1\tLiu_2017\tPATENL
PARNXJ\tchr1\t47239296\tC\tACG\t1\tLiu_2017\tPARNXJ
PASXSI\tchr1\t47239296\tC\tAACC\t1\tLiu_2017\tPASXSI
PASNEH\tchr1\t47239296\tC\tACC\t1\tLiu_2017\tPASNEH
PAUAFN\tchr1\t47239296\tC\tAACC\t1\tLiu_2017\tPAUAFN
PARASZ\tchr1\t47239296\tC\tACC\t1\tLiu_2017\tPARASZ
PARWNW\tchr1\t47239296\tC\tACC\t1\tLiu_2017\tPARWNW
PASFKA\tchr1\t47239293\tTTA\tACCGTTAATCAA\t1\tLiu_2017\tPASFKA
PATEIT\tchr1\t47239296\tC\tAC\t1\tLiu_2017\tPATEIT
PASMHF\tchr1\t47239296\tC\tAC\t1\tLiu_2017\tPASMHF
PARJNX\tchr1\t47239296\tC\tAC\t1\tLiu_2017\tPARJNX
PASYWF\tchr1\t47239296\tC\tAC\t1\tLiu_2017\tPASYWF
"""
    return pd.read_table(io.StringIO(variant_data), sep='\t')


def vcf_row_to_variant(vcf_row: pd.Series) -> genome.Variant:
    """Parse a row of a vcf df into a genome.Variant."""
    variant = genome.Variant(
        chromosome=str(vcf_row.CHROM),
        position=int(vcf_row.POS),
        reference_bases=vcf_row.REF,
        alternate_bases=vcf_row.ALT,
        name=vcf_row.ID,
    )
    return variant


def generate_background_variants(
    variant: genome.Variant, max_number: int = 100
) -> pd.DataFrame:
    """Generates a dataframe of background variants for a given variant.
    
    This is done by creating new sequences of the same length as the alternate
    allele.
    
    This allows us to test if the specific sequence of the oncogenic variant has a
    greater effect than a random sequence of the same length at the same location.
    
    Args:
        variant: The variant to generate ism variants for.
        max_number: The maximum number of ism variants to generate.
        
    Returns:
        A dataframe of variants.
    """
    nucleotides = np.array(list('ACGT'), dtype='<U1')
    
    def generate_unique_strings(n, max_number, random_seed=42):
        """Generates unique random strings of length n."""
        rng = np.random.default_rng(random_seed)
        
        if 4**n < max_number:
            raise ValueError(
                'Cannot generate that many unique strings for the given length.'
            )
        
        generated_strings = set()
        while len(generated_strings) < max_number:
            indices = rng.integers(0, 4, size=n)
            new_string = ''.join(nucleotides[indices])
            if new_string != variant.alternate_bases:
                generated_strings.add(new_string)
        return list(generated_strings)
    
    permutations = []
    if 4 ** len(variant.alternate_bases) < max_number:
        # Get all
        for p in itertools.product(
            nucleotides, repeat=len(variant.alternate_bases)
        ):
            permutations.append(''.join(p))
    else:
        # Sample some
        permutations = generate_unique_strings(
            len(variant.alternate_bases), max_number
        )
    ism_candidates = pd.DataFrame({
        'ID': ['mut_' + str(variant.position) + '_' + x for x in permutations],
        'CHROM': variant.chromosome,
        'POS': variant.position,
        'REF': variant.reference_bases,
        'ALT': permutations,
        'output': 0.0,
        'original_variant': variant.name,
    })
    return ism_candidates


def oncogenic_and_background_variants(
    input_sequence_length: int, number_of_background_variants: int = 20
) -> pd.DataFrame:
    """Generates a dataframe of all variants for this evaluation."""
    oncogenic_variants = oncogenic_tal1_variants()
    
    variants = []
    for vcf_row in oncogenic_variants.itertuples():
        variants.append(
            genome.Variant(
                chromosome=str(vcf_row.CHROM),
                position=int(vcf_row.POS),
                reference_bases=vcf_row.REF,
                alternate_bases=vcf_row.ALT,
                name=vcf_row.ID,
            )
        )
    
    background_variants = pd.concat([
        generate_background_variants(variant, number_of_background_variants)
        for variant in variants
    ])
    all_variants = pd.concat([oncogenic_variants, background_variants])
    return inference_df(all_variants, input_sequence_length=input_sequence_length)


def inference_df(
    qtl_df: pd.DataFrame,
    input_sequence_length: int,
) -> pd.DataFrame:
    """Returns a pd.DataFrame with variants and intervals ready for inference."""
    df = []
    for _, row in qtl_df.iterrows():
        variant = vcf_row_to_variant(row)
        
        interval = genome.Interval(
            chromosome=row['CHROM'], start=row['POS'], end=row['POS']
        ).resize(input_sequence_length)
        
        df.append({
            'interval': interval,
            'variant': variant,
            'output': row['output'],
            'variant_id': row['ID'],
            'POS': row['POS'],
            'REF': row['REF'],
            'ALT': row['ALT'],
            'CHROM': row['CHROM'],
        })
    return pd.DataFrame(df)


def coarse_grained_mute_groups(eval_df):
    """Assigns variants to coarse-grained groups for plotting."""
    grp = []
    for row in eval_df.itertuples():
        if row.POS >= 47239290:  # MUTE site.
            if row.ALT_len > 4:
                grp.append('MUTE' + '_other')
            else:
                grp.append('MUTE' + '_' + str(row.ALT_len))
        else:
            grp.append(str(row.POS) + '_' + str(row.ALT_len))
    
    grp = pd.Series(grp)
    return pd.Categorical(grp, categories=sorted(grp.unique()), ordered=True)


@example_analysis_workflow_mcp.tool()
def alphagenome_visualize_variant_positions(
    chromosome: Annotated[str,
        "Chromosome identifier (e.g., 'chr1'). "
        "Must match chromosome format used in GENCODE annotation (typically 'chrN' format)."] = "chr1",
    
    interval_start: Annotated[int,
        "Start position of genomic interval to visualize (0-based coordinate). "
        "For TAL1 locus example: 47209255. Should encompass the gene and surrounding regulatory regions."] = 47209255,
    
    interval_end: Annotated[int,
        "End position of genomic interval to visualize (0-based coordinate). "
        "For TAL1 locus example: 47242023. Should encompass the gene and surrounding regulatory regions."] = 47242023,
    
    strand: Annotated[Literal["+", "-"],
        "Strand orientation of the gene. '+' for positive strand, '-' for negative strand. "
        "TAL1 gene is on negative strand ('-')."] = "-",
    
    variant_positions: Annotated[list[int] | None,
        "List of variant positions to mark on the plot (0-based coordinates). "
        "If None, uses oncogenic TAL1 variants from tutorial. "
        "Example: [47212072, 47212074, 47230639, 47239291, 47239293, 47239295, 47239296]"] = None,
    
    position_labels: Annotated[list[str] | None,
        "Custom labels for each variant position. Must match length of variant_positions. "
        "Use empty string '' to skip label for a position. "
        "If None, generates default labels from positions."] = None,
    
    plot_title: Annotated[str,
        "Title for the visualization plot. "
        "Should describe the genomic region and variants being visualized."] = "Positions of variants near TAL1",
    
    api_key: Annotated[str | None,
        "AlphaGenome API key for accessing the model. "
        "If None, reads from ALPHAGENOME_API_KEY environment variable. "
        "Get your key from Google AI Studio."] = None,
    
    out_prefix: Annotated[str | None,
        "Prefix for output file names. If None, uses 'variant_positions_{timestamp}'."] = None,
) -> dict:
    """
    Visualize genomic positions of variants in a specified interval with gene annotations.
    Input is genomic interval coordinates and variant positions, output is annotated genomic map figure.
    """
    # Input validation
    if interval_start >= interval_end:
        raise ValueError(f"interval_start ({interval_start}) must be less than interval_end ({interval_end})")
    
    # Set API key
    if api_key is None:
        api_key = os.environ.get('ALPHAGENOME_API_KEY')
        if api_key is None:
            raise ValueError("API key must be provided via api_key parameter or ALPHAGENOME_API_KEY environment variable")
    
    # Create DNA model client
    dna_model = dna_client.create(api_key)
    
    # Load gene annotations from GENCODE
    gtf = pd.read_feather(
        'https://storage.googleapis.com/alphagenome/reference/gencode/'
        'hg38/gencode.v46.annotation.gtf.gz.feather'
    )
    
    # Filter to MANE select transcripts (one curated transcript per gene)
    gtf_transcript = gene_annotation.filter_protein_coding(gtf)
    gtf_transcript = gene_annotation.filter_to_mane_select_transcript(gtf_transcript)
    transcript_extractor = transcript_utils.TranscriptExtractor(gtf_transcript)
    
    # Define genomic interval
    interval = genome.Interval(
        chromosome=chromosome, 
        start=interval_start, 
        end=interval_end, 
        strand=strand
    )
    
    # Get variant positions
    if variant_positions is None:
        # Use oncogenic TAL1 variants from tutorial
        oncogenic_variants = oncogenic_tal1_variants()
        variant_positions_array = oncogenic_variants['POS'].unique()
        variant_positions = variant_positions_array.tolist()
        variant_positions.sort()
        
        # Use tutorial labels
        position_labels = [
            '47212072, 47212074',
            '',
            '47230639',
            '47239291 - 47239296',
            '',
            '',
            '',
        ]
    else:
        if position_labels is None:
            # Generate default labels from positions
            position_labels = [str(pos) for pos in variant_positions]
        elif len(position_labels) != len(variant_positions):
            raise ValueError(
                f"position_labels length ({len(position_labels)}) must match "
                f"variant_positions length ({len(variant_positions)})"
            )
    
    # Build plot
    fig = plot_components.plot(
        [
            plot_components.TranscriptAnnotation(
                transcript_extractor.extract(interval)
            ),
        ],
        annotations=[
            plot_components.VariantAnnotation(
                [
                    genome.Variant(
                        chromosome=chromosome,
                        position=x,
                        reference_bases='N',
                        alternate_bases='N',
                    )
                    for x in variant_positions
                ],
                labels=position_labels,
                use_default_labels=False,
            )
        ],
        interval=interval,
        title=plot_title,
    )
    
    # Save figure
    if out_prefix is None:
        out_prefix = f"variant_positions_{timestamp}"
    
    output_file = OUTPUT_DIR / f"{out_prefix}.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    return {
        "message": f"Visualized {len(variant_positions)} variant positions in {chromosome}:{interval_start}-{interval_end}",
        "reference": "https://github.com/google-deepmind/alphagenome/blob/main/colabs/example_analysis_workflow.ipynb",
        "artifacts": [
            {
                "description": "Variant position visualization",
                "path": str(output_file.resolve())
            }
        ]
    }


@example_analysis_workflow_mcp.tool()
def alphagenome_predict_variant_effects(
    position: Annotated[int,
        "Genomic position of the variant (0-based coordinate). "
        "Example: 47239296 for the Jurkat variant in TAL1 locus."],
    
    reference_bases: Annotated[str,
        "Reference allele sequence at this position. "
        "Example: 'C' for single nucleotide, 'CTAACC' for deletion. "
        "Must match reference genome at specified position."],
    
    alternate_bases: Annotated[str,
        "Alternate allele sequence (variant). "
        "Example: 'CCGTTTCCTAACC' for insertion, 'T' for SNV. "
        "Can be longer (insertion) or shorter (deletion) than reference."],
    
    chromosome: Annotated[str,
        "Chromosome identifier where variant is located (e.g., 'chr1'). "
        "Must match reference genome format."] = "chr1",
    
    variant_name: Annotated[str,
        "Identifier or name for this variant. "
        "Used in output labels and files. Example: 'Jurkat', 'Patient_1', 'rs123456'."] = "variant",
    
    interval_start: Annotated[int,
        "Start position of genomic interval to analyze (0-based coordinate). "
        "Should encompass the gene of interest and surrounding regulatory regions. "
        "For TAL1 locus: 47209255."] = 47209255,
    
    interval_end: Annotated[int,
        "End position of genomic interval to analyze (0-based coordinate). "
        "For TAL1 locus: 47242023."] = 47242023,
    
    strand: Annotated[Literal["+", "-"],
        "Strand orientation of the gene of interest. "
        "'+' for positive strand, '-' for negative strand (TAL1 is '-')."] = "-",
    
    ontology_terms: Annotated[list[str],
        "Cell Line Ontology (CLO) or Cell Ontology (CL) terms for biological context. "
        "Example: ['CL:0001059'] for CD34-positive common myeloid progenitor cells. "
        "Multiple terms can be provided for different cell types."] = ["CL:0001059"],
    
    input_sequence_length: Annotated[int,
        "Length of sequence for model input (in base pairs). "
        "Must be power of 2. Tutorial uses 2^20 = 1,048,576 bp. "
        "Larger values provide more genomic context but slower inference."] = 1048576,
    
    plot_title: Annotated[str | None,
        "Title for the effect visualization plot. "
        "If None, generates descriptive title with variant information."] = None,
    
    api_key: Annotated[str | None,
        "AlphaGenome API key for accessing the model. "
        "If None, reads from ALPHAGENOME_API_KEY environment variable."] = None,
    
    out_prefix: Annotated[str | None,
        "Prefix for output file names. If None, uses 'variant_effects_{timestamp}'."] = None,
) -> dict:
    """
    Predict functional impact of a specific variant on gene expression, DNA accessibility, and histone modifications.
    Input is variant coordinates and biological context, output is comprehensive multi-track visualization figure.
    """
    # Input validation
    if interval_start >= interval_end:
        raise ValueError(f"interval_start must be less than interval_end")
    
    if not (position >= interval_start and position <= interval_end):
        raise ValueError(
            f"Variant position ({position}) must be within interval "
            f"({interval_start}-{interval_end})"
        )
    
    # Validate input_sequence_length is power of 2
    if input_sequence_length & (input_sequence_length - 1) != 0:
        raise ValueError(f"input_sequence_length must be a power of 2, got {input_sequence_length}")
    
    # Set API key
    if api_key is None:
        api_key = os.environ.get('ALPHAGENOME_API_KEY')
        if api_key is None:
            raise ValueError("API key must be provided via api_key parameter or ALPHAGENOME_API_KEY environment variable")
    
    # Create DNA model client
    dna_model = dna_client.create(api_key)
    
    # Load gene annotations
    gtf = pd.read_feather(
        'https://storage.googleapis.com/alphagenome/reference/gencode/'
        'hg38/gencode.v46.annotation.gtf.gz.feather'
    )
    gtf_transcript = gene_annotation.filter_protein_coding(gtf)
    gtf_transcript = gene_annotation.filter_to_mane_select_transcript(gtf_transcript)
    transcript_extractor = transcript_utils.TranscriptExtractor(gtf_transcript)
    
    # Define genomic interval
    interval = genome.Interval(
        chromosome=chromosome,
        start=interval_start,
        end=interval_end,
        strand=strand
    )
    
    # Create variant object
    variant = genome.Variant(
        chromosome=chromosome,
        position=position,
        reference_bases=reference_bases,
        alternate_bases=alternate_bases,
        name=variant_name,
    )
    
    # Make predictions for REF and ALT alleles
    output = dna_model.predict_variant(
        interval=interval.resize(input_sequence_length),
        variant=variant,
        requested_outputs={
            dna_client.OutputType.RNA_SEQ,
            dna_client.OutputType.CHIP_HISTONE,
            dna_client.OutputType.DNASE,
        },
        ontology_terms=ontology_terms,
    )
    
    # Build visualization plot
    transcripts = transcript_extractor.extract(interval)
    
    if plot_title is None:
        plot_title = (
            f'Effect of variant on predicted RNA Expression, DNAse, and ChIP-Histone\n'
            f'variant={variant}'
        )
    
    fig = plot_components.plot(
        [
            plot_components.TranscriptAnnotation(transcripts),
            # RNA-seq tracks
            plot_components.Tracks(
                tdata=output.alternate.rna_seq.filter_to_nonpositive_strand()
                - output.reference.rna_seq.filter_to_nonpositive_strand(),
                ylabel_template='{biosample_name} ({strand})\n{name}',
                filled=True,
            ),
            # DNAse tracks
            plot_components.Tracks(
                tdata=output.alternate.dnase.filter_to_nonpositive_strand()
                - output.reference.dnase.filter_to_nonpositive_strand(),
                ylabel_template='{biosample_name} ({strand})\n{name}',
                filled=True,
            ),
            # ChIP-Histone tracks
            plot_components.Tracks(
                tdata=output.alternate.chip_histone.filter_to_nonpositive_strand()
                - output.reference.chip_histone.filter_to_nonpositive_strand(),
                ylabel_template='{biosample_name} ({strand})\n{name}',
                filled=True,
            ),
        ],
        annotations=[plot_components.VariantAnnotation([variant])],
        interval=interval,
        title=plot_title,
    )
    
    # Save figure
    if out_prefix is None:
        out_prefix = f"variant_effects_{timestamp}"
    
    output_file = OUTPUT_DIR / f"{out_prefix}.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    return {
        "message": f"Predicted effects of {variant_name} at {chromosome}:{position}",
        "reference": "https://github.com/google-deepmind/alphagenome/blob/main/colabs/example_analysis_workflow.ipynb",
        "artifacts": [
            {
                "description": "Variant effect visualization",
                "path": str(output_file.resolve())
            }
        ]
    }


@example_analysis_workflow_mcp.tool()
def alphagenome_compare_variant_scores(
    target_gene: Annotated[str,
        "Gene symbol to analyze expression changes. "
        "Example: 'TAL1' for T-cell acute lymphoblastic leukemia analysis. "
        "Must be present in GENCODE annotation."] = "TAL1",
    
    ontology_curie: Annotated[str,
        "Cell Ontology CURIE (Compact URI) for cell type of interest. "
        "Example: 'CL:0001059' for CD34-positive common myeloid progenitor cells. "
        "Format: 'CL:NNNNNNN' or 'CLO:NNNNNNN'."] = "CL:0001059",
    
    input_sequence_length: Annotated[int,
        "Length of sequence for model input (in base pairs). "
        "Must be power of 2. Tutorial uses 2^20 = 1,048,576 bp."] = 1048576,
    
    number_of_background_variants: Annotated[int,
        "Number of random background variants to generate per oncogenic variant. "
        "Background variants are random sequences of same length at same position. "
        "Typical range: 3-100. Tutorial uses 3 for demo (fast), 20 for production."] = 3,
    
    max_workers: Annotated[int,
        "Number of parallel workers for variant scoring. "
        "Higher values speed up computation but use more API quota. "
        "Typical range: 1-4."] = 2,
    
    api_key: Annotated[str | None,
        "AlphaGenome API key for accessing the model. "
        "If None, reads from ALPHAGENOME_API_KEY environment variable."] = None,
    
    out_prefix: Annotated[str | None,
        "Prefix for output file names. If None, uses 'variant_comparison_{timestamp}'."] = None,
) -> dict:
    """
    Compare predicted gene expression changes between cancer-associated variants and shuffled background variants.
    Input is target gene and analysis parameters, output is distribution plots and comparison statistics CSV.
    """
    # Input validation
    if input_sequence_length & (input_sequence_length - 1) != 0:
        raise ValueError(f"input_sequence_length must be a power of 2, got {input_sequence_length}")
    
    if number_of_background_variants < 1:
        raise ValueError(f"number_of_background_variants must be >= 1, got {number_of_background_variants}")
    
    if max_workers < 1:
        raise ValueError(f"max_workers must be >= 1, got {max_workers}")
    
    # Set API key
    if api_key is None:
        api_key = os.environ.get('ALPHAGENOME_API_KEY')
        if api_key is None:
            raise ValueError("API key must be provided via api_key parameter or ALPHAGENOME_API_KEY environment variable")
    
    # Create DNA model client
    dna_model = dna_client.create(api_key)
    
    # Prepare variant groups (oncogenic + background)
    eval_df = oncogenic_and_background_variants(
        input_sequence_length=input_sequence_length, 
        number_of_background_variants=number_of_background_variants
    )
    
    # Add variant annotations
    eval_df['ALT_len'] = eval_df['ALT'].str.len()
    eval_df['variant_group'] = (
        eval_df['POS'].astype(str) + '_' + eval_df['ALT_len'].astype(str)
    )
    eval_df['output'] = eval_df['output'].fillna(0) != 0
    eval_df['coarse_grained_variant_group'] = coarse_grained_mute_groups(eval_df)
    
    # Score all variants
    scores = dna_model.score_variants(
        intervals=eval_df['interval'].to_list(),
        variants=eval_df['variant'].to_list(),
        variant_scorers=[variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ']],
        max_workers=max_workers,
    )
    
    # Extract target gene scores for specified cell type
    gene_index = scores[0][0].obs.query(f'gene_name == "{target_gene}"').index[0]
    cell_type_index = scores[0][0].var.query(f'ontology_curie == "{ontology_curie}"').index[0]
    
    def get_gene_score(score_data):
        """Extracts the gene expression score from model output."""
        return score_data[gene_index, cell_type_index].X[0, 0]
    
    eval_df[f'{target_gene.lower()}_diff'] = [
        get_gene_score(x[0]) for x in scores
    ]
    
    # Prepare plot data
    plot_df = eval_df.loc[eval_df.REF != eval_df.ALT].copy()
    plot_df['variant'] = plot_df['variant'].astype(str)
    plot_df = plot_df.loc[
        :,
        [
            'variant',
            'output',
            f'{target_gene.lower()}_diff',
            'coarse_grained_variant_group',
        ],
    ].drop_duplicates()
    
    # Define plot titles
    facet_title_by_group = {
        '47212072_22': 'chr1:47212072\n21 bp ins.',
        '47212074_7': 'chr1:47212072\n21 bp ins.',
        '47230639_1': 'chr1:47230639\nSNV',
        'MUTE_2': 'chr1:47239296\n1 bp ins.',
        'MUTE_3': 'chr1:47239296\n2 bp ins.',
        'MUTE_4': 'chr1:47239296\n3 bp ins.',
        'MUTE_other': 'chr1:47239296\n7-18 bp ins.',
    }
    
    # Generate plots for each variant group
    plt_dict = {}
    for group in plot_df.coarse_grained_variant_group.unique():
        subplot_df = pd.concat(
            [plot_df.assign(plot_group='density'), plot_df.assign(plot_group='rain')]
        )
        subplot_df = subplot_df[subplot_df.coarse_grained_variant_group == group]
        subplot_df = subplot_df[
            ~((subplot_df.plot_group == 'density') & (subplot_df.output))
        ]
        
        col_width = np.ptp(subplot_df[f'{target_gene.lower()}_diff']) / 200
        subplot_df['col_width'] = subplot_df['output'].map(
            {True: 1.5 * col_width, False: 1.25 * col_width}
        )
        
        plt_ = (
            gg.ggplot(subplot_df)
            + gg.aes(x=f'{target_gene.lower()}_diff')
            + gg.geom_col(
                gg.aes(
                    y=1,
                    width='col_width',
                    fill='output',
                    x=f'{target_gene.lower()}_diff',
                    alpha='output',
                ),
                data=subplot_df[subplot_df['plot_group'] == 'rain'],
            )
            + gg.geom_density(
                gg.aes(
                    x=f'{target_gene.lower()}_diff',
                    fill='output',
                ),
                data=subplot_df[subplot_df['plot_group'] == 'density'],
                color='white',
            )
            + gg.facet_wrap('~output + plot_group', nrow=1, scales='free_x')
            + gg.scale_alpha_manual({True: 1, False: 0.3})
            + gg.scale_fill_manual({True: '#FAA41A', False: 'gray'})
            + gg.labs(title=facet_title_by_group.get(group, group))
            + gg.theme_minimal()
            + gg.geom_vline(xintercept=0, linetype='dotted')
            + gg.theme(
                figure_size=(1.2, 3),
                legend_position='none',
                axis_text_x=gg.element_blank(),
                panel_grid_major_x=gg.element_blank(),
                panel_grid_minor_x=gg.element_blank(),
                strip_text=gg.element_blank(),
                axis_title_y=gg.element_blank(),
                axis_title_x=gg.element_blank(),
                plot_title=gg.element_text(size=9),
            )
            + gg.scale_y_reverse()
            + gg.coord_flip()
        )
        
        plt_dict[group] = plt_
    
    # Save plots and data
    if out_prefix is None:
        out_prefix = f"variant_comparison_{timestamp}"
    
    artifacts = []
    
    # Save each plot
    for group, plot_obj in plt_dict.items():
        plot_file = OUTPUT_DIR / f"{out_prefix}_{group}.png"
        plot_obj.save(plot_file, dpi=300, verbose=False)
        artifacts.append({
            "description": f"Comparison plot for {group}",
            "path": str(plot_file.resolve())
        })
    
    # Save comparison data
    data_file = OUTPUT_DIR / f"{out_prefix}_scores.csv"
    plot_df.to_csv(data_file, index=False)
    artifacts.append({
        "description": "Variant scores comparison data",
        "path": str(data_file.resolve())
    })
    
    return {
        "message": f"Compared {len(plot_df)} variants across {len(plt_dict)} variant groups for {target_gene} expression",
        "reference": "https://github.com/google-deepmind/alphagenome/blob/main/colabs/example_analysis_workflow.ipynb",
        "artifacts": artifacts
    }
