#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: modas_plots.py
Author: CrazyHsu @ crazyhsu9527@gmail.com
Last modified: 2021-09-25 23:45:23
'''

'''
Miscellaneous plot functions to visualize the results of MODAS, which can be invoked by various subcommand. 
Details are: "plot_manhattan", "plot_qq", "plot_genotype_box", "plot_forest", "plot_net"
"plot_scatter", "plot_heatmap", "plot_nucleotide_dst", "plot_scatter_ps", "plot_html"

'''

import argparse, sys
from functions import *

def main(args):
    if args.command == "plot_phylo":
        phylo_plot(args.plink_bed, args.out_prefix, group_file=args.group_file, group_sep=args.group_sep,
                   selected_lines=args.sub_lines, drop_line=args.drop_lines)
    if args.command == "plot_scatterps":
        scatterps_plot(args.input, group_file=args.group_file, ps_type=args.ps_type, out_plot_prefix=args.out_prefix,
                        input_sep=args.input_sep, group_sep=args.group_sep)
    if args.command == "plot_nd":
        nucleotide_dst_plot1(args.plink_bed, args.group_file, args.gff_annotation, gene_list=args.select_genes,
                            chrom=args.chrom, chr_start=int(args.chrom_start)-50, chr_end=int(args.chrom_end)+50,
                            diversity_method=args.method, plot_by_gene=args.plot_by_gene, keep_temp=args.keep_temp,
                            out_dir=args.out_dir, out_prefix=args.out_prefix,
                            left_offset=args.left_offset, right_offset=args.right_offset,
                            window_size=args.window_size, window_step=args.window_step, smooth=args.smooth,
                            plot_left_expand=args.plot_left_expand, plot_right_expand=args.plot_right_expand)

    if args.command == "plot_manhattan":
        thresholdi = [args.thresholdD, args.thresholdL, args.thresholdU]
        manhattan_plot(args.input, args.out_dir, thresholdi=thresholdi, gwas_sep=args.input_sep,
                       data_from=args.data_from, threads=int(args.threads), dpi=int(args.dpi))
    if args.command == "plot_qq":
        qq_plot(args.input, args.out_dir, gwas_sep=args.input_sep, data_from=args.data_from,
                threads=int(args.threads), dpi=int(args.dpi), threshold=float(args.threshold))

    if args.command == "plot_forest":
        forest_plot(args.input, args.out_prefix, order_by_name=args.order_by_name,
                    order_by_effect=args.order_by_effect, order_by_group=args.order_by_group, group=str(args.group),
                    sep=args.input_sep, height=float(args.height), width=float(args.width))
    if args.command == "plot_scattermr":
        scattermr_plot(args.input, args.out_prefix, group_file=args.group_file, sep=args.input_sep,
                       height=float(args.height), width=float(args.width))

    if args.command == "plot_phebox":
        box_plot(args.input, args.snp, args.group, args.selected_genes, args.selected_snps,
                 args.selected_strains, str(args.scale), args.plot_type, args.out_prefix)

    if args.command == "plot_net":
        net_plot(args.input, args.out_prefix, args.input_sep, args.pvalue, module_size=args.module_size)
    if args.command == "plot_go":
        go_plot(args.input, args.gene2go, args.bg_sep, plot_type=args.plot_type, sigp=args.sigp,
                adjustp=args.adjustp, out_prefix=args.out_prefix)

    if args.command == "plot_heatmap":
        heatmap(args.input, args.out_prefix, select_rows=str(args.row_select), select_cols=str(args.col_select),
                cluster_rows=args.cluster_rows, cluster_cols=args.cluster_cols,
                order_col_by=args.order_col_by, order_row_by=args.order_row_by,
                anno_row=args.anno_row, anno_col=args.anno_col, scale_by=args.scale_by,
                show_colnames=args.show_colnames, show_rownames=args.show_rownames,
                height=args.height, width=args.width)
    if args.command == "plot_html":
        pass

if __name__ == '__main__':
    USAGE = ' Miscellaneous plot functions to visualize the results of MODAS. '

    parser = argparse.ArgumentParser(usage='%(prog)s command [options]', description=USAGE)
    subparsers = parser.add_subparsers(title='command', metavar='', dest='command', prog=parser.prog)

    # phylogenetic tree plot
    plot_phylo = subparsers.add_parser('plot_phylo', help="Construct phylogenetic tree",
                                       usage='%(prog)s [options]')
    plot_phylo.add_argument('-bed', dest="plink_bed", help="The plink bed file")
    plot_phylo.add_argument('-o', dest="out_prefix", default="MODASv", help="The output file prefix")
    plot_phylo.add_argument('-g', dest="group_file", default=None,
                            help="The group information used to color the sub-population. Default: None")
    plot_phylo.add_argument('-gs', dest="group_sep", default=",",
                            help="The separator of group file. Default: ','")
    plot_phylo.add_argument('-sub', dest="sub_lines", default=None,
                            help="The selected strains to analysis, with family ID and starin ID per line. Default: None")
    plot_phylo.add_argument('-drop', dest="drop_lines", default=None,
                            help="The drop lines which you don't want to draw in the phylogenetic tree. Default: None")

    # scatter plot for population structure
    plot_scatterps = subparsers.add_parser('plot_scatterps', help='Draw PCA or tSNE plots for population structure',
                                           usage='%(prog)s [options]')
    plot_scatterps.add_argument('-i', dest='input', required=True, help="The input population structure file.")
    plot_scatterps.add_argument('-g', dest="group_file", required=True,
                                help="The group information used to color the sub-population. Default: None")
    plot_scatterps.add_argument('-o', dest="out_prefix", default="MODASv", help="The output file prefix")
    plot_scatterps.add_argument('-ps_type', dest="ps_type", default="pca",
                                help="The population type you want to draw: pca or tsne")
    plot_scatterps.add_argument('-input_sep', dest="input_sep", default=",",
                                help="The separator of input file. Default: ','")
    plot_scatterps.add_argument('-group_sep', dest="group_sep", default=",",
                                help="The separator of input file. Default: ','")

    # manhattan plot
    plot_mht = subparsers.add_parser('plot_manhattan', help='Draw manhattan plots for GWAS results',
                                     usage='%(prog)s [options]')
    plot_mht.add_argument('-i', dest="input", help="The GWAS results")
    plot_mht.add_argument("-od", dest="out_dir", default="MODASv",
                          help="The output directory, and manhattan plot of each "
                               "association file will be suffixed with '_manhattan.jpg'")
    plot_mht.add_argument("-threshold_low", dest="thresholdL", type=float, default=0.000001,
                          help="The low p-value threshold for plotting")
    plot_mht.add_argument("-threshold_upper", dest="thresholdU", type=float, default=0.00001,
                          help="The high p-value threshold for plotting")
    plot_mht.add_argument("-threshold_default", dest="thresholdD", type=float, default=None,
                          help="The data-orient p-value threshold for plotting")
    plot_mht.add_argument("-data_from", dest="data_from", default="file",
                          help="The type of GWAS output stored, 'file' or 'list' or 'directory'. Default: file")
    plot_mht.add_argument("-input_sep", dest="input_sep", default="\t",
                          help="The separator in GWAS file. Default: tab")
    plot_mht.add_argument("-t", dest="threads", default=10,
                          help="The number of threads used to parallel draw the plot. Default: 10")
    plot_mht.add_argument("-dpi", dest="dpi", type=int, default=300,
                          help="The DPI of the plot. Default: 300")

    # QQ plot
    plot_qq = subparsers.add_parser('plot_qq', help='Draw QQ plots for GWAS results',
                                     usage='%(prog)s [options]')
    plot_qq.add_argument('-i', dest="input", help="The GWAS results")
    plot_qq.add_argument("-od", dest="out_dir", default="MODASv",
                         help="The output directory, and qq plot of each "
                              "association file will be suffixed with '_qq.jpg'")
    plot_qq.add_argument("-data_from", dest="data_from", default="file",
                         help="The type of GWAS output stored, 'file' or 'list' or 'directory'. Default: file")
    plot_qq.add_argument("-input_sep", dest="input_sep", default="\t",
                         help="The separator in GWAS file. Default: tab")
    plot_qq.add_argument("-t", dest="threads", type=int, default=10,
                         help="The number of threads used to parallel draw the plot. Default: 10")
    plot_qq.add_argument("-dpi", dest="dpi", type=int, default=300,
                         help="The DPI of the plot. Default: 300")
    plot_qq.add_argument("-threshold", dest="threshold", type=float, default=0.000001,
                         help="The signal threshold for Q-Q plot. Default: 1e-6")

    # phenotype boxplot
    plot_ptb = subparsers.add_parser('plot_phebox', help='Draw box-line plots for genotypes in different strains',
                                     usage='%(prog)s [options]')
    plot_ptb.add_argument('-i', dest="input", help="The expression file")
    plot_ptb.add_argument('-snp', dest="snp", help="The SNP file")
    plot_ptb.add_argument('-group', dest="group", help="The group file")
    plot_ptb.add_argument('-s_genes', dest="selected_genes", default=None,
                          help="The file listing selected genes or a string with comma separated")
    plot_ptb.add_argument('-s_snps', dest="selected_snps", default=None,
                          help="The file listing selected SNPs or a string with comma separated")
    plot_ptb.add_argument('-s_strains', dest="selected_strains", default=None,
                          help="The file listing selected strains or a string with comma separated")
    plot_ptb.add_argument('-scale', dest="scale", default=False, action="store_true", help="Scale the phenotype values")
    plot_ptb.add_argument('-plot_type', dest="plot_type", choices=["single_gene", "multi_genes", "single_pop"],
                          default="single_gene", help="The type of graph to be plotted. "
                                                      "Options are: single_gene, multi_gene, single_pop. "
                                                      "Default: single_gene")
    plot_ptb.add_argument("-o", dest="out_prefix", default="MODASv", help="The output file prefix")

    # forest plot
    plot_forest = subparsers.add_parser('plot_forest', help='Draw forest plots for interested genes under different '
                                                            'traits or conditions using MR results',
                                        usage='%(prog)s [options]')
    plot_forest.add_argument('-i', dest="input", help="The MR results file")
    plot_forest.add_argument("-o", dest="out_prefix", default="MODASv", help="The output file prefix")
    plot_forest.add_argument("-order_by_name", dest="order_by_name", action="store_true", default=False,
                             help="Order the 'exposure' by name. Default: False")
    plot_forest.add_argument("-order_by_effect", dest="order_by_effect", action="store_true", default=False,
                             help="Order the 'exposure' by effect. Default: False")
    plot_forest.add_argument("-order_by_group", dest="order_by_group", action="store_true", default=False,
                             help="Order the 'exposure' by group. Default: False")
    plot_forest.add_argument("-group", dest="group", required="-order_by_group" in sys.argv, default=None,
                             help="The group file used to order the output. Required when '-order_by_group' specified")
    plot_forest.add_argument("-input_sep", dest="input_sep", default=",",
                             help="The separator in MR result. Default: comma")
    plot_forest.add_argument("-height", dest="height", default=6, type=float,
                             help="The height of output plot. Default: 6")
    plot_forest.add_argument("-width", dest="width", default=4, type=float,
                             help="The width of output plot. Default: 4")

    # scatter plot for MR results
    plot_scattermr = subparsers.add_parser('plot_scattermr', help='Draw scatter plots for interested genes under '
                                                                  'different traits or conditions using MR results',
                                           usage='%(prog)s [options]')
    plot_scattermr.add_argument('-i', dest="input", help="The MR results file")
    plot_scattermr.add_argument("-o", dest="out_prefix", default="MODASv", help="The output file prefix")
    plot_scattermr.add_argument("-group_file", dest="group_file", default=None,
                                help="The group file used to color the pTrait")
    plot_scattermr.add_argument("-input_sep", dest="input_sep", default=",",
                                help="The separator in MR result. Default: comma")
    plot_scattermr.add_argument("-height", dest="height", default=5, type=float,
                                help="The height of output plot. Default: ")
    plot_scattermr.add_argument("-width", dest="width", default=9, type=float,
                                help="The width of output plot. Default: ")

    # MR-based network construction and visualization
    plot_net = subparsers.add_parser('plot_net', help='MR-based network construction and visualization',
                                           usage='%(prog)s [options]')
    plot_net.add_argument('-i', dest="input", help="The MR results file")
    plot_net.add_argument("-o", dest="out_prefix", default="MODASv", help="The output file prefix")
    plot_net.add_argument("-input_sep", dest="input_sep", default=",",
                          help="The separator in MR result. Default: comma")
    plot_net.add_argument("-sig_p", dest="pvalue", default=0.001, type=float,
                          help="pvalue cutoff for MR-based network analysis. Default: 1e-3")
    plot_net.add_argument("-ms", dest="module_size", default=5, type=int,
                          help="The minimal size of genes in a module. Default: 5")

    # Enrich plot with various types
    plot_go = subparsers.add_parser('plot_go', help='Draw enrichment plot with various plot type',
                                    usage='%(prog)s [options]')
    plot_go.add_argument('-i', dest="input", required=True, help="The gene list used to be enriched")
    plot_go.add_argument("-o", dest="out_prefix", default="MODASv", help="The output file prefix")
    plot_go.add_argument('-bg', dest="gene2go", required=True, help="The background file containing gene to GO pairs")
    plot_go.add_argument("-bg_sep", dest="bg_sep", default="\t",
                         help="The separator of pairs in gene2go file. Default: tab")
    plot_go.add_argument('-plot_type', dest="plot_type", choices=["barplot", "dotplot", "emap"], default="dotplot",
                         help="The visualization type of plots you want to implement. Default: dotplot"
                              "Choices are: barplot, dotplot, emap")
    plot_go.add_argument('-sigp', dest="sigp", type=float, default=0.05,
                         help="The significant p-value used to filter the results. Default: 0.05")
    plot_go.add_argument('-adjustp', dest="adjustp", action="store_true", default=False,
                         help="Adjust p-value. Default: False")

    # heatmap plot
    plot_heatmap = subparsers.add_parser('plot_heatmap', help='Draw heat map for various genes under different conditions',
                                     usage='%(prog)s [options]')
    plot_heatmap.add_argument('-i', dest="input", help="The matrix used to draw heatmap")
    plot_heatmap.add_argument("-o", dest="out_prefix", default="MODASv", help="The output file prefix")
    plot_heatmap.add_argument("-row_select", dest="row_select", default=None, help="The selected rows")
    plot_heatmap.add_argument("-col_select", dest="col_select", default=None, help="The selected columns")
    plot_heatmap.add_argument("-anno_row", dest="anno_row", default=None, help="Annotate the rows with a file")
    plot_heatmap.add_argument("-anno_col", dest="anno_col", default=None, help="Annotate the columns with a file")
    plot_heatmap.add_argument("-cluster_rows", dest="cluster_rows", action="store_true", default=False,
                              help="Cluster rows")
    plot_heatmap.add_argument("-cluster_cols", dest="cluster_cols", action="store_true", default=False,
                              help="Cluster columns")
    plot_heatmap.add_argument("-order_row_by", dest="order_row_by", default=None,
                              help="Order rows by a file containing a list of names")
    plot_heatmap.add_argument("-order_col_by", dest="order_col_by", default=None,
                              help="Order columns by a file containing a list of names")
    plot_heatmap.add_argument("-scale_by", dest="scale_by", default="none", choices=["none", "row", "column"],
                              help="Scale the values in either the 'row' direction or the 'column' direction, or 'none'. "
                                   "Default: 'none'")
    plot_heatmap.add_argument("-show_rownames", dest="show_rownames", action="store_true", default=False,
                              help="Show row names")
    plot_heatmap.add_argument("-show_colnames", dest="show_colnames", action="store_true", default=False,
                              help="Show column names")
    plot_heatmap.add_argument("-height", dest="height", type=float, default=5,
                              help="The height of output plot. Default: 5")
    plot_heatmap.add_argument("-width", dest="width", type=float, default=5,
                              help="The width of output plot. Default: 5")

    # nucleotide diversity
    plot_nd = subparsers.add_parser('plot_nd', help='Draw nucleotide density plots for a genome region in a population',
                                     usage='%(prog)s [options]')
    plot_nd.add_argument('-bed', dest="plink_bed", help="The plink bed file")
    plot_nd.add_argument("-o", dest="out_prefix", default="MODASv", help="The prefix of output files. Default: MODASv")
    plot_nd.add_argument("-od", dest="out_dir", default="MODASv", help="The output directory. Default: MODASv")
    plot_nd.add_argument("-group", dest="group_file", help="The group file")
    plot_nd.add_argument("-gff_anno", dest="gff_annotation", help="The annotation with format in 'gff'")
    plot_nd.add_argument("-select_genes", dest="select_genes", default=None,
                         help="The file containing a list of genes with one per row or a regular string containing "
                              "genes with comma separated.")
    plot_nd.add_argument("-chrom", dest="chrom", default=None, type=str, help="Specify the chromosome id")
    plot_nd.add_argument("-chrom_start", dest="chrom_start", default=0, type=int, help="Specify the start of interval")
    plot_nd.add_argument("-chrom_end", dest="chrom_end", default=0, type=int, help="Specify the end of interval")
    plot_nd.add_argument("-method", dest="method", default="pi", help="The method used to calculate the diversity")
    plot_nd.add_argument("-plot_by_gene", dest="plot_by_gene", action="store_true", default=False,
                         help="The draw a single plot for every gene. Default: False")
    plot_nd.add_argument("-keep_temp", dest="keep_temp", action="store_true", default=False,
                         help="Keep the temporary files or not. Default: False")
    plot_nd.add_argument("-left_offset", dest="left_offset", type=int, default=0,
                         help="The left offset to expand the gene region. Default: 0")
    plot_nd.add_argument("-right_offset", dest="right_offset", type=int, default=0,
                         help="The left offset to expand the gene region. Default: 0")
    plot_nd.add_argument("-plot_left_expand", dest="plot_left_expand", type=int, default=500,
                         help="The left expanding to plot. Default: 500")
    plot_nd.add_argument("-plot_right_expand", dest="plot_right_expand", type=int, default=500,
                         help="The right expanding to plot. Default: 500")
    plot_nd.add_argument("-window_size", dest="window_size", type=int, default=1000,
                         help="The window size to calculate PI, fst or Tajima's D. Default: 1000")
    plot_nd.add_argument("-window_step", dest="window_step", type=int, default=50,
                         help="The window step to calculate PI, fst or Tajima's D. Default: 50")
    plot_nd.add_argument("-smooth", dest="smooth", action="store_true", default=False,
                         help="Smooth the line plot. Default: False")
    # plink_bed, group_file, gff_annotation, gene_list = "", chrom = None, chr_start = 0, chr_end = 0,
    # diversity_method = "pi", split_by_gene = True, keep_temp = True, out_dir = "nd_output"

    # # html report
    # plot_html = subparsers.add_parser('plot_html', help='Aggregate the results in a html report',
    #                                  usage='%(prog)s [options]')
    # plot_html.add_argument('-c', dest="default_cfg", help="The config file used for init setting.")

    if len(sys.argv) <= 2:
        sys.argv.append('-h')
    args = parser.parse_args()
    main(args)
