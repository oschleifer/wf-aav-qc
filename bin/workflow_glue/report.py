"""Create workflow report."""
import json

from dominate.tags import p
import ezcharts as ezc
from ezcharts.plots import Plot
from ezcharts.components import fastcat
from ezcharts.components.ezchart import EZChart
from ezcharts.components.reports import labs
from ezcharts.layout.snippets import Grid, Tabs
from ezcharts.layout.snippets.table import DataTable
import numpy as np
import pandas as pd

from .util import get_named_logger, wf_parser  # noqa: ABS101


def plot_trucations(report, truncations_file):
    """Make a report section with start and end position histograms.

    The truncations_file contains start and end positions of alignments that are fully
    contained within the ITR-ITR regions.
    """
    df = pd.read_csv(
        truncations_file, sep='\t',
        dtype={
            'Read start': str,
            'Read end': np.uint32,
            'sample_id': str
        }
    )

    with report.add_section("Truncations", "Truncations"):
        p(
            "This plot illustrates the frequency of start and end positions of "
            "alignments that map completely within the transgene plasmid ITR-ITR "
            "region, helping to identify potential truncation hotspots."
        )
        tabs = Tabs()
        with tabs.add_dropdown_menu():

            for sample, df_sample in df.groupby('sample_id'):
                with tabs.add_dropdown_tab(sample):
                    df_sample.drop(columns=['sample_id'], inplace=True)
                    plt = ezc.histplot(data=df_sample, binwidth=5)
                    plt.xAxis = dict(name='Transgene cassette genome position')
                    plt.yAxis = dict(name='Number of alignments')
                    plt.legend = dict(orient='horizontal', top=30)
                    EZChart(plt, theme='epi2melabs')


def plot_itr_coverage(report, coverage_file):
    """Make report section with ITR-ITR coverage of transgene cassette region."""
    df = pd.read_csv(
        coverage_file,
        sep=r"\s+",
        dtype={
            'ref': str,
            'pos': np.uint32,
            'depth': np.uint32,
            'strand': str,
            'sample_id': str
        })

    with report.add_section("ITR-ITR coverage", "Coverage"):
        p(
            "For each transgene reference, sequencing depth is calculated "
            "for both forward and reverse mapping alignments."

        )
        tabs = Tabs()
        with tabs.add_dropdown_menu():

            for sample, df_sample in df.groupby('sample_id'):
                with tabs.add_dropdown_tab(sample):
                    with Grid(columns=1):
                        for ref, df_ref in df_sample.groupby('ref'):
                            plt = ezc.lineplot(
                                data=df_ref, x='pos', y='depth', hue='strand')
                            plt.title = dict(text=ref)
                            plt.legend = dict(
                                orient='horizontal', top=30, icon='rect')
                            for s in plt.series:
                                s.showSymbol = False
                            EZChart(plt, theme='epi2melabs', height='300px')


def plot_contamination(report, class_counts):
    """Make report section with contamination plots.

    Two plots: (1) mapped/unmapped; (2) mapped reads per reference
    """
    df_class_counts = pd.read_csv(
        class_counts,
        sep='\t',
        dtype={
            'Reference': str,
            'Number of alignments': np.uint32,
            'Percentage of alignments': np.float32,
            'sample_id': str
        }
    )

    with report.add_section("Contamination", "Contamination"):
        p(
            "These two plots show mapping summaries that can highlight "
            "potential contamination issues."
        )
        p(
            "The first plot shows the percentage of reads that either map to any "
            "combined reference sequence or are unmapped."
        )
        p(
            "The second plot breaks down the the alignment numbers into the "
            "specific references (host, helper plasmid, Rep-Cap plasmid, and transgene "
            "plasmid)."
        )

        with Grid(columns=2):
            # Prepare data for combined plots
            df_reads = df_class_counts[df_class_counts.Reference.isin(['Mapped', 'Unmapped'])].copy()
            df_reads['Percentage of Reads'] = df_reads['Percentage of alignments']

            df_alns = df_class_counts[~df_class_counts.Reference.isin(['Mapped', 'Unmapped'])].copy()

            # Sort the DataFrame by 'Reference' to have 'Unmapped' first and then 'Mapped'
            df_reads['Reference'] = pd.Categorical(df_reads['Reference'], categories=['Unmapped', 'Mapped'], ordered=True)
            df_reads = df_reads.sort_values(by='Reference')

            # Combine Reference and sample_id for unique x-axis labels
            df_reads['Reference_Sample'] = df_reads['Reference'] + '_' + df_reads['sample_id']
            df_alns['Reference_Sample'] = df_alns['Reference'] + '_' + df_alns['sample_id']


            # Create the first plot: Reads mapped/unmapped
            plt_reads = ezc.barplot(
                df_reads,
                x='Reference_Sample', y='Percentage of Reads',
                color='sample_id', group='Reference'
            )
            plt_reads.title = dict(text='Reads mapped/unmapped')
            plt_reads.xAxis.axisLabel = dict(rotate=45)
            EZChart(plt_reads, theme='epi2melabs', height='400px')

            # Create the second plot: Alignment counts per target
            plt_alns = ezc.barplot(
                df_alns,
                x='Reference_Sample', y='Percentage of alignments',
                color='sample_id', group='Reference'
            )
            plt_alns.title = dict(text='Alignment counts per target')
            plt_alns.xAxis.axisLabel = dict(rotate=45)
            EZChart(plt_alns, theme='epi2melabs', height='400px')

def plot_read_summary(report, stats):
    """Make report section barplots detailing the read quality, read length, and base yield."""
    df_stats = pd.read_csv(
        stats,
        sep='\t',
        dtype={
            'read_id': str,
            'filename': str,
            'sample_name': str,
            'read_length':np.uint32,
            'mean_quality': np.float32,
            'channel': np.uint32,
            'read_number': np.uint32
        })
    with report.add_section("Read Summary", "Read Summary"):
        p(
            "Read quality, read length, base yield"
        )
        
        with Grid(columns=3):
            # Histogram of read quality
            hist_quality, edges_quality = np.histogram(df_stats['mean_quality'], bins=50)
            df_quality = pd.DataFrame({
                'Quality Score': edges_quality[:-1],
                'Number of Reads': hist_quality
            })
            plt_quality = ezc.barplot(
                df_quality,
                x='Quality Score', y='Number of Reads'
            )
            plt_quality.title = dict(text='Read Quality')
            EZChart(plt_quality, theme='epi2melabs', height='400px')

            # Histogram of read lengths
            hist_length, edges_length = np.histogram(df_stats['read_length'], bins=50)
            df_length = pd.DataFrame({
                'Read Length': edges_length[:-1],
                'Number of Reads': hist_length
            })
            plt_length = ezc.barplot(
                df_length,
                x='Read Length', y='Number of Reads'
            )
            plt_length.title = dict(text='Read Length')
            EZChart(plt_length, theme='epi2melabs', height='400px')

            # Line graph of base yield
            hist_yield, edges_yield = np.histogram(df_stats['read_length'], bins=50)
            df_yield = pd.DataFrame({
                'Read Length': edges_yield[:-1],
                'Cumulative Bases': hist_yield
            })
            df_stats['cumulative_bases'] = df_stats['read_length'].cumsum()
            plt_yield = ezc.lineplot(
                df_yield,
                x='Read Length', y='Cumulative Bases'
            )
            plt_yield.title = dict(text='Base yield above read length')
            EZChart(plt_yield, theme='epi2melabs', height='400px')

def plot_aav_structures(report, structures_file):
    """Make report section barplots detailing the AAV structures found."""
    df = pd.read_csv(
        structures_file,
        sep='\t',
        dtype={
            'Assigned_genome_type': str,
            'count': np.uint32,
            'percentage': np.float32,
            'sample_id': str

        })

    with report.add_section("AAV Structures", "Structures"):
        p(
            "The numbers of different of the various AAV transgene genome types  "
            "identified in the sample(s) are summarised here."
        )

        p(
            "A detailed report containing more granular genome type assignments "
            "per read can be found at:"
            " `output/<sample_id>/<sample_id>_aav_per_read_info.tsv` "
        )
        tabs = Tabs()
        with tabs.add_dropdown_menu():

            for sample, df_sample in df.groupby('sample_id'):
                with tabs.add_dropdown_tab(sample):
                    df_sample = df_sample.sort_values('percentage', ascending=False)
                    # Plot of main genome type counts
                    plt = ezc.barplot(
                        df_sample,
                        x='Assigned_genome_type',
                        y='percentage')
                    plt.title = dict(text='Genome types')
                    plt.xAxis.axisLabel = dict(rotate=45)
                    EZChart(plt, theme='epi2melabs')

                    # Table with counts and percentages
                    # (in lieu of being able to annotate bar plots in ezchrts)
                    DataTable.from_pandas(df_sample, use_index=False)


def main(args):
    """Run the entry point."""
    logger = get_named_logger("Report")
    report = labs.LabsReport(
        "AAV QC workflow report", "wf-aav-qc",
        args.params, args.versions)

    with open(args.metadata) as metadata:
        sample_details = sorted([
            {
                'sample': d['alias'],
                'type': d['type'],
                'barcode': d['barcode']
            } for d in json.load(metadata)
        ], key=lambda d: d["sample"])

    if args.stats:
        plot_read_summary(report, args.stats[0])
        # with report.add_section("Read summary", "Read summary"):
        #     # TODO fix this. Do we need o concat stats?
        #     fastcat.SeqSummary(args.stats[0])

    plot_contamination(
        report,
        args.contam_class_counts)
    plot_trucations(report, args.truncations)
    plot_itr_coverage(report, args.itr_coverage)
    plot_aav_structures(report, args.aav_structures)

    with report.add_section("Metadata", "Metadata"):
        tabs = Tabs()
        with tabs.add_dropdown_menu():
            for d in sample_details:
                with tabs.add_dropdown_tab(d["sample"]):
                    df = pd.DataFrame.from_dict(
                        d, orient="index", columns=["Value"])
                    df.index.name = "Key"
                    DataTable.from_pandas(df)

    report.write(args.report)
    logger.info(f"Report written to {args.report}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument(
        "report", help="Report output file")
    parser.add_argument(
        "--stats", nargs='*', help="Fastcat per-read stats file(s).")
    parser.add_argument(
        "--truncations", help="TSV with start and end columns for.")
    parser.add_argument(
        "--itr_coverage", help="TSV with alignment Pos and EndPos columns.")
    parser.add_argument(
        "--contam_class_counts", help="TSV of reference mapping counts.")
    parser.add_argument(
        "--aav_structures", help="TSV of reads with AAV structure assignment.")
    parser.add_argument(
        "--metadata", default='metadata.json',
        help="sample metadata")
    parser.add_argument(
        "--versions", required=True,
        help="directory containing CSVs containing name,version.")
    parser.add_argument(
        "--params", default=None, required=True,
        help="A JSON file containing the workflow parameter key/values")
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")
    return parser
