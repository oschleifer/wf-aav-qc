"""Create workflow report."""
import json

from dominate.tags import p
import ezcharts as ezc
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

def plot_contamination(): return
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

            # Create the first plot: Reads mapped/unmapped
            plt_reads = ezc.barplot(
                df_reads, width=1,
                x='Reference', y='Percentage of Reads',
                hue='sample_id', group='Reference'
            )
            plt_reads.title = dict(text='Reads mapped/unmapped')
            plt_reads.xAxis.axisLabel = dict(rotate=45)
            EZChart(plt_reads, theme='epi2melabs', height='400px')

            # Create the second plot: Alignment counts per target
            plt_alns = ezc.barplot(
                df_alns, width=1,
                x='Reference', y='Percentage of alignments',
                hue='sample_id', group='Reference'
            )
            plt_alns.title = dict(text='Alignment counts per target')
            plt_alns.xAxis.axisLabel = dict(rotate=45)
            EZChart(plt_alns, theme='epi2melabs', height='400px')

def plot_quality(df_stats):
    """Helper function to plot the quality scores."""
    df_quality = df_stats.groupby(['sample_name', 'binned_quality']).size().reset_index(name='read_count')
    df_quality['binned_quality'] = df_quality['binned_quality'].astype(str)
    df_quality['binned_quality'] = df_quality['binned_quality'].apply(lambda x: float(x.split(",")[0][1:]))
    
    combined_qstats = pd.DataFrame()
    for sample_name in df_quality['sample_name'].unique():
        barcode = df_quality[df_quality['sample_name'] == sample_name]
        combined_qstats = pd.concat([combined_qstats, barcode], ignore_index=True)
    
    combined_qstats = combined_qstats.rename(columns={
        'sample_name': 'Barcode',
        'binned_quality': 'Quality score',
        'read_count': 'Number of reads'
    })
    
    plt_quality = ezc.lineplot(
        data=combined_qstats, hue='Barcode',
        x='Quality score', y='Number of reads'
    )
    for series in plt_quality.series:
        series.lineStyle = {'opacity': 0.8}
        series.showSymbol = False
    plt_quality.title = dict(
        text="Read quality",
        subtext=(
            f"Mean: {round(df_stats['mean_quality'].mean(), 1)} "
            f"Median: {round(df_stats['mean_quality'].median())} "
        )
    )
    plt_quality.legend = dict(orient='horizontal', right='right', top=40, icon='rect')
    EZChart(plt_quality, theme='epi2melabs', height='400px')

def plot_length(df_stats):
    """Helper function to plot the read lengths."""
    df_length = df_stats.groupby(['sample_name', 'binned_length']).size().reset_index(name='read_l_count')
    df_length['binned_length'] = df_length['binned_length'].astype(str)
    df_length['binned_length'] = df_length['binned_length'].apply(lambda x: float(x.split(",")[0][1:]))
    
    combined_length = pd.DataFrame()
    for sample_name in df_length['sample_name'].unique():
        barcode = df_length[df_length['sample_name'] == sample_name]
        combined_length = pd.concat([combined_length, barcode], ignore_index=True)
    
    combined_length = combined_length.rename(columns={
        'sample_name': 'Barcode',
        'binned_length': 'Read length / kb',
        'read_l_count': 'Number of reads'
    })
    
    plt_length = ezc.lineplot(
        data=combined_length, hue='Barcode',
        x='Read length / kb', y='Number of reads'
    )
    for series in plt_length.series:
        series.lineStyle = {'opacity': 0.8}
        series.showSymbol = False
    plt_length.title = dict(
        text="Read length",
        subtext=(
            f"Mean: {round(df_stats['read_length'].mean() * 1000)} "
            f"Median: {round(df_stats['read_length'].median() * 1000)} "
            f"Min: {round(df_stats['read_length'].min() * 1000)} "
            f"Max: {round(df_stats['read_length'].max() * 1000)} "
        )
    )
    plt_length.xAxis.min = 0
    plt_length.xAxis.max = max(df_stats['read_length'])
    plt_length.legend = dict(orient='vertical', right='right', top=40, icon='rect')
    EZChart(plt_length, theme='epi2melabs', height='400px')

def plot_yield(df_stats):
    """Helper function to plot the base yield above read length."""
    combined_yield = pd.DataFrame()
    for sample_name, df_sample in df_stats.groupby('sample_name'):
        length = np.concatenate(([0], np.sort(df_sample["read_length"])), dtype="float")
        cumsum = np.cumsum(length[::-1])[::-1]
        df_yield = pd.DataFrame({
            'Read length / kb': length,
            'Yield above length / Gbases': cumsum / 1e9,
            'Barcode': sample_name
        })
        
        combined_yield = pd.concat([combined_yield, df_yield], ignore_index=True)
    
    plt_yield = ezc.lineplot(
        data=combined_yield, hue='Barcode', 
        x='Read length / kb', y='Yield above length / Gbases'
    )
    for series in plt_yield.series:
        series.lineStyle = {'opacity': 0.8}
        series.showSymbol = False
    plt_yield.title = dict(
        text="Base yield above read length"
    )
    plt_yield.legend = dict(orient='horizontal', right='right', top=30, icon='rect')
    EZChart(plt_yield, theme='epi2melabs', height='400px')

def plot_read_summary(report, stats):
    """Make report section barplots detailing the read quality, read length, and base yield."""
    # Load and prepare the data for plotting.
    df_stats = pd.read_csv(
        stats,
        sep='\t',
        dtype={
            'read_id': str,
            'filename': str,
            'sample_name': str,
            'read_length': np.uint32,
            'mean_quality': np.float32,
            'channel': np.uint32,
            'read_number': np.uint32
        }
    )
    
    # bin data
    qbins = np.arange(0, df_stats['mean_quality'].max() + 1, 0.2)
    df_stats['read_length'] = df_stats['read_length'] / 1000
    lbins = np.arange(0, df_stats['read_length'].max() + 1, 0.1)
    df_stats['binned_quality'] = pd.cut(df_stats['mean_quality'], qbins)
    df_stats['binned_length'] = pd.cut(df_stats['read_length'], lbins)
    
    with report.add_section("Read summary", "Read summary"):
        with Grid(columns=3):
            plot_quality(df_stats)
            plot_length(df_stats)
            plot_yield(df_stats)


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
        with report.add_section("Read summary", "Read summary"):
            # TODO fix this. Do we need o concat stats?
            fastcat.SeqSummary(args.stats[0])

    plot_contamination(
        report,
        args.contam_class_counts
        )
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
