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
    if args.command == "plot_manhattan":
        thresholdi = getManhattanThresholdInterval(args.thresholdD, args.thresholdL, args.thresholdU)
        manhattan_plot(args.input, args.out_prefix, thresholdi=thresholdi)
    if args.command == "plot_qq":
        qq_plot(args.input, args.out_prefix)

    if args.command == "plot_forest":
        forest_plot(args.input, args.out_prefix, order=args.order, sep=args.sep)
    if args.command == "plot_scatter":
        scatter_plot(args.input, args.out_prefix, group_file=args.group_file, sep=args.sep)

    if args.command == "plot_genotype_box":
        box_plot(args.input, args.snp_file, args.group_file, args.selected_genes, args.selected_snps,
                 args.selected_strains, args.scale, args.plot_type, args.out_prefix)
    if args.command == "plot_scatter_ps":
        scatter_ps_plot(args.input, group_file=args.group_file, ps_type=args.ps_type, out_plot_prefix=args.out_prefix, sep=args.sep)
    if args.command == "plot_nucleotide_dst":
        nucleotide_dst_plot()

    if args.command == "plot_net":
        pass
    if args.command == "plot_go":
        pass

    if args.command == "plot_heatmap":
        heatmap(args.input, args.out_prefix, cluster_rows=args.cluster_rows, cluster_cols=args.cluster_cols,
                order_col_by=args.order_col_by, order_row_by=args.order_row_by,
                anno_row=args.anno_row, anno_col=args.anno_col, scale=args.scale,
                show_colnames=args.show_colnames, show_rownames=args.show_rownames)
    if args.command == "plot_html":
        pass

if __name__ == '__main__':
    USAGE = ' Miscellaneous plot functions to visualize the results of MODAS. '

    parser = argparse.ArgumentParser(usage='%(prog)s command [options]', description=USAGE)
    subparsers = parser.add_subparsers(title='command', metavar='', dest='command', prog=parser.prog)
    plot_mht = subparsers.add_parser('plot_manhattan', help='Draw manhattan plots for GWAS results',
                                     usage='%(prog)s [options]')
    plot_mht.add_argument('-i', dest="input", help="The GWAS results")
    plot_mht.add_argument("-o", dest="out_prefix", help="The output file prefix")
    plot_mht.add_argument("-threshold_low", dest="thresholdL", help="The low p-value threshold for plotting")
    plot_mht.add_argument("-threshold_upper", dest="thresholdU", help="The high p-value threshold for plotting")
    plot_mht.add_argument("-threshold_default", dest="thresholdD", help="The data-orient p-value threshold for plotting")

    plot_qq = subparsers.add_parser('plot_qq', help='Draw QQ plots for GWAS results',
                                     usage='%(prog)s [options]')
    plot_qq.add_argument('-i', dest="input", help="The GWAS results")
    plot_qq.add_argument("-o", dest="out_prefix", help="The output file prefix")

    plot_gtb = subparsers.add_parser('plot_genotype_box', help='Draw box-line plots for genotypes in different strains',
                                     usage='%(prog)s [options]')
    plot_gtb.add_argument('-c', dest="default_cfg", help="The config file used for init setting.")

    plot_forest = subparsers.add_parser('plot_forest', help='Draw forest plots for interested genes under different traits or conditions using MR results',
                                     usage='%(prog)s [options]')
    plot_forest.add_argument('-c', dest="default_cfg", help="The config file used for init setting.")

    plot_scatter = subparsers.add_parser('plot_scatter', help='Draw scatter plots for interested genes under different traits or conditions using MR results',
                                     usage='%(prog)s [options]')
    plot_scatter.add_argument('-c', dest="default_cfg", help="The config file used for init setting.")

    plot_heatmap = subparsers.add_parser('plot_heatmap', help='Draw heat map for various genes under different conditions',
                                     usage='%(prog)s [options]')
    plot_heatmap.add_argument('-c', dest="default_cfg", help="The config file used for init setting.")

    plot_nuc_dst = subparsers.add_parser('plot_nucleotide_dst', help='Draw nucleotide density plots for a genome region in a population',
                                     usage='%(prog)s [options]')
    plot_nuc_dst.add_argument('-c', dest="default_cfg", help="The config file used for init setting.")

    plot_scatter_ps = subparsers.add_parser('plot_scatter_ps', help='Draw PCA or tSNE plots for population structure',
                                     usage='%(prog)s [options]')
    plot_scatter_ps.add_argument('-c', dest="default_cfg", help="The config file used for init setting.")

    plot_html = subparsers.add_parser('plot_html', help='Aggregate the results in a html report',
                                     usage='%(prog)s [options]')
    plot_html.add_argument('-c', dest="default_cfg", help="The config file used for init setting.")

    if len(sys.argv) <= 2:
        sys.argv.append('-h')
    args = parser.parse_args()
    main(args)
