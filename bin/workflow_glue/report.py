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
from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource, FactorRange
from bokeh.plotting import figure, show
from bokeh.embed import file_html
from bokeh.resources import CDN
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

        # Prepare data for combined plots
        df_reads = df_class_counts[df_class_counts.Reference.isin(['Mapped', 'Unmapped'])].copy()
        df_reads['Percentage of Reads'] = df_reads['Percentage of alignments']

        df_alns = df_class_counts[~df_class_counts.Reference.isin(['Mapped', 'Unmapped'])].copy()

        # Create the first plot: Reads mapped/unmapped
        source_reads = ColumnDataSource(data=dict(
            x=[(sample, reference) for sample in df_reads['sample_id'].unique() for reference in ['Mapped', 'Unmapped']],
            counts=[df_reads[(df_reads['sample_id'] == sample) & (df_reads['Reference'] == reference)]['Percentage of Reads'].values[0] 
                    if not df_reads[(df_reads['sample_id'] == sample) & (df_reads['Reference'] == reference)].empty else 0
                    for sample in df_reads['sample_id'].unique() for reference in ['Mapped', 'Unmapped']]
        ))

        p1 = figure(x_range=FactorRange(*source_reads.data['x']), height=400, title="Reads mapped/unmapped",
                    toolbar_location=None, tools="")

        p1.vbar(x='x', top='counts', width=0.9, source=source_reads)

        p1.xgrid.grid_line_color = None
        p1.y_range.start = 0

        # Create the second plot: Alignment counts per target
        unique_references = df_alns['Reference'].unique()
        source_alns = ColumnDataSource(data=dict(
            x=[(sample, reference) for sample in df_alns['sample_id'].unique() for reference in unique_references],
            counts=[df_alns[(df_alns['sample_id'] == sample) & (df_alns['Reference'] == reference)]['Percentage of alignments'].values[0] 
                    if not df_alns[(df_alns['sample_id'] == sample) & (df_alns['Reference'] == reference)].empty else 0
                    for sample in df_alns['sample_id'].unique() for reference in unique_references]
        ))

        p2 = figure(x_range=FactorRange(*source_alns.data['x']), height=400, title="Alignment counts per target",
                    toolbar_location=None, tools="")

        p2.vbar(x='x', top='counts', width=0.9, source=source_alns)

        p2.xgrid.grid_line_color = None
        p2.y_range.start = 0

        # Combine the plots into a grid
        grid = gridplot([[p1, p2]], sizing_mode='stretch_both')

        # Output the plots to an HTML file
        html = file_html(grid, CDN, "Contamination Plots")
        with open("contamination_plots.html", "w") as f:
            f.write(html)

def plot_read_summary(report, stats):
    """Make report section barplots detailing the read quality, read length, and base yield."""
    df = pd.read_csv(
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
    with report.add_section("Read Summary", "Reads"):
        p(
            "Read quality, read length, base yield"
        )
        # with Grid(columns=3):
        #     plt = ezc.barplot(
        #         df[['read_length']]
        #     )
        #     plt.title = dict(text='Read lengths')
        #     EZChart(plt,  theme='epi2melabs', height='400px')
            # plt = Plot()
            # for sample, df_sample in df.groupby('sample_name'):
            #     df_read_lengths = df_sample.sort_values('read_length', ascending=True)
            #     plt.dataset(
            #         df_read_lengths[['read_length']]
            #     )
            # plt.series = [dict(type='bar')]
            # plt.title = dict(text='Read lengths')
            # EZChart(plt, theme='epi2melabs', height='400px')

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
