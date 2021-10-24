#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: functions.py
Author: CrazyHsu @ crazyhsu9527@gmail.com
Last modified: 2021-09-28 13:27:48
'''

import os, re

import pandas as pd
import rpy2.robjects as robjects

############# Functions #############
def validateFile(myFile):
    if not os.path.exists(myFile):
        raise Exception("File '%s' not found! Please input again!" % myFile)

    if not os.path.isfile(myFile):
        raise Exception("File '%s' is not a file! Please input again!" % myFile)

    return True

def validateDir(myDir):
    if not os.path.exists(myDir):
        raise Exception("Dir '%s' not found! Please input again!" % myDir)

    if not os.path.isdir(myDir):
        raise Exception("Dir '%s' is not a directory! Please input again!" % myDir)

    return True

def getManhattanThresholdInterval(thresholdD, thresholdL, thresholdU):
    if not thresholdL:
        thresholdL = 0.000001
    else:
        thresholdL = float(thresholdL)
    if not thresholdU:
        thresholdU = 0.00001
    else:
        thresholdU = float(thresholdU)
    if not thresholdD:
        thresholdD = None
    else:
        thresholdD = float(thresholdD)
    return [thresholdD, thresholdL, thresholdU]

#####################################

############# GWAS results plots #############
'''
### The format of GWAS file (separate can be automated detected):
rs chr ps p_wald
1_1922301 1 1922301 9.121183e-03
1_1928050 1 1928050 1.795902e-03
1_2521954 1 2521954 7.200593e-03
1_2522874 1 2522874 6.791745e-03
'''

def manhattan_plot(gwas_res, out_plot_prefix, thresholdi=None, sep="\t"):
    thresholdi = ",".join(map(str, thresholdi))
    from plotWithR import plotManhattanPlotR
    robjects.r(plotManhattanPlotR)
    robjects.r.plotManhattan(gwas_res, thresholdi, out_plot_prefix, sep=sep)

def qq_plot(gwas_res, out_plot_prefix, sep="\t"):
    from plotWithR import plotQQplotR
    robjects.r(plotQQplotR)
    robjects.r.plotQQplot(gwas_res, out_plot_prefix, sep=sep)

##############################################

############# MR results plots #############
''' the format of MR results
snp,mTrait,pTrait,effect,TMR,pvalue
'''

def forest_plot(mr_res, out_plot_prefix, order=False, sep=","):
    from plotWithR import plotForestPlotR
    robjects.r(plotForestPlotR)
    robjects.r.plotForestPlot(mr_res, order, out_plot_prefix, sep=sep)

def scatter_plot(mr_res, out_plot_prefix, group_file=None, sep=","):
    from plotWithR import plotScatterPlotR
    if group_file != None:
        if not validateFile(group_file):
            group_file = None
    robjects.r(plotScatterPlotR)
    robjects.r.plotScatterPlot(mr_res, str(group_file), out_plot_prefix, sep=sep)

############################################

############# Genotype results plots #############
''' 
### The format of population structure:
RIL,PC1,PC2
GEMS58,38.88,-59.69
CML415,-30.07,1.74
SI273,52.52,11.30
CIMBL135,-38.45,9.53

### The format of group file:
RIL,subpop
GEMS58,NSS
CML415,TST
SI273,NSS
CIMBL135,TST
'''
def scatter_ps_plot(ps_file, group_file, ps_type, out_plot_prefix, sep=","):
    if ps_type.lower() == "pca":
        from plotWithR import plotScatterPsR
        robjects.r(plotScatterPsR)
        robjects.r.plotScatterPs(ps_file, group_file, out_plot_prefix, sep=sep)
    elif ps_type.lower == "tsne":
        pass
    else:
        raise Exception("Please input the correct poplation structure type you want to plot")
##################################################

############# Heatmap #############
def heatmap(input, out_prefix, cluster_rows=True, cluster_cols=True, order_col_by=None, order_row_by=None,
            anno_row=None, anno_col=None, scale="none", show_colnames=False, show_rownames=False):
    from plotWithR import plotHeatmapR
    [cluster_rows, cluster_cols, show_colnames, show_rownames] = ["T" if x else "F" for x in [cluster_rows, cluster_cols, show_colnames, show_rownames]]
    if order_col_by:
        if not validateFile(order_col_by):
            order_col_by = "None"
    if order_row_by:
        if not validateFile(order_row_by):
            order_row_by = "None"
    if anno_row:
        if not validateFile(anno_row):
            anno_row = "None"
    if anno_col:
        if not validateFile(anno_col):
            anno_col = "None"

    robjects.r(plotHeatmapR)
    robjects.r.plotHeatmap(input, out_prefix, cluster_rows, cluster_cols, order_col_by, order_row_by,
                           anno_row, anno_col, scale, show_rownames, show_colnames)

###################################

############# Box plots #############
'''
### The format of input file
RIL	gene1	gene2   ...
CAU1	0.735651576	0.617622876
CAU2	0.345792276	0.043156773
CAU3	0.161468224	0.467046847
CAU4	0.417096623	0.485693651
...

### The format of snp file
RIL	chr4.s_10000	chr4.s_10001    ...
CAU1	A	G
CAU2	A	G
CAU3	A	G
CAU4	A	G
...

### The format of group file
RIL	subpop 
CAU1	temp
CAU2	trop
CAU3	temp
CAU4	trop
...
'''
def box_plot(input_mat_file, snp_file, group_file, selected_genes, selected_snps, selected_strains, scale, plot_type, out_prefix):
    from plotWithR import plotBoxR
    robjects.r(plotBoxR)
    robjects.r.plotBox(input_mat_file, snp_file, group_file, selected_genes, selected_snps, selected_strains, scale, plot_type, out_prefix)
###################################

############# Nucleotide density Plots #############
def getLongestTrans(gpeList):
    pass

def extract_snp(plink_bed, snp_outfile, chrom_range, plink_path=None):
    if not plink_path:
        plink_path = 'plink1.9'
    command = plink_path + ' -bfile ' + plink_bed + ' --out ' + snp_outfile + \
              ' --chr ' + str(chrom_range[0]) + ' --from-bp ' + str(chrom_range[1]) + \
              ' --to-bp ' + str(chrom_range[2]) + ' --recode vcf-iid'
    os.system(command)
    os.system('rm ' + snp_outfile + '.log')
    os.system('rm ' + snp_outfile + '.nosex')
    return snp_outfile + ".vcf"

def merge_nd_file(nd_path, method="pi"):
    if method == 'pi':
        first_columns = ['BIN_START', 'N_VARIANTS', 'PI']
        merge_columns = ['BIN_START', 'PI']
    elif method == 'tajimad':
        first_columns = ['BIN_START', 'N_SNPS', 'TajimaD']
        merge_columns = ['BIN_START', 'TajimaD']
    elif method == 'fst':
        first_columns = ['BIN_START', 'N_VARIANTS', 'MEAN_FST']
        merge_columns = ['BIN_START', 'MEAN_FST']
    else:
        raise TypeError("method = ['pi', 'tajimad', 'fst']")

    nd_data = pd.DataFrame()
    for i, path in enumerate(nd_path):
        data = pd.read_csv(path, header=0, index_col=None, sep='\t')
        if i == 0:
            nd_data = data.loc[:, first_columns]
        else:
            nd_data = nd_data.merge(data.loc[:, merge_columns], on='BIN_START')
    return nd_data

def calc_n_diversity(snp_file, method="pi", vcftools_path=None, sample_file=None, out_path=None, nd_range=None,
                     fst_sample=None, window_size=1000, window_step=50):
    if not vcftools_path:
        vcftools_path = 'vcftools'
    if not out_path:
        dir_path = '/'.join(snp_file.split('/')[:-1])
        file_name = '.'.join(snp_file.split('/')[-1].split('.')[:-1])
        out_path = dir_path + '/' + file_name

    command = vcftools_path + ' --vcf ' + snp_file + ' --out ' + out_path
    if sample_file:
        command = command + ' --keep ' + sample_file
    if nd_range:
        command = command + ' --chr ' + str(nd_range[0]) + ' --from-bp ' + str(nd_range[1]) + ' --to-bp ' + str(nd_range[2])
    if method == 'pi':
        command = command + ' --window-pi ' + str(window_size) + ' --window-pi-step ' + str(window_step)
        suffix = '.windowed.pi'
    elif method == 'tajimad':
        command = command + ' --TajimaD ' + str(window_size)
        suffix = '.Tajima.D'
    elif method == 'fst':
        command = command + ' --weir-fst-pop ' + fst_sample[0] + ' --weir-fst-pop ' + fst_sample[1] \
                  + ' --fst-window-size ' + str(window_size) + ' --fst-window-step ' + str(window_step)
        suffix = '.windowed.weir.fst'
    else:
        raise TypeError("method = ['pi', 'tajimad', 'fst']")
    os.system(command)
    return out_path + suffix

def nucleotide_dst_plot(plink_bed, group_file, gff_annotation, gene_list=None, chrom=None, chr_start=0, chr_end=0, density_method="pi", split_by_gene=True, keep_temp=True, out_dir="nd_output"):
    from commonClasses import GenePredExtLine, GenePredObj
    from dna_features_viewer import GraphicFeature, GraphicRecord

    exon_color = "#F19800"
    gpe_annotation = ""
    if gff_annotation.endswith("gpe"):
        gpe_annotation = gff_annotation
    elif gff_annotation.endswith("gtf"):
        gpe_annotation = ".".join([os.path.splitext(gff_annotation)[0], "gpe"])
        command = "gtfToGenePred {} {} -genePredExt".format(gff_annotation, gpe_annotation)
        os.system(command)
    elif gff_annotation.endswith("gff") or gff_annotation.endswith("gff3"):
        gpe_annotation = ".".join([os.path.splitext(gff_annotation)[0], "gpe"])
        command = "gffread -T -o tmp.gtf {} ".format(gff_annotation)
        os.system(command)
        command = "gtfToGenePred {} {} -genePredExt".format(gff_annotation, gpe_annotation)
        os.system(command)
        os.system('rm tmp.gtf')
    else:
        raise Exception("Please input the annotation with suffix including: gpe, gff3, gff, gtf.")

    group_df = pd.read_csv(group_file, sep=",")
    group_df = group_df.iloc[:, [0, 1]]
    group_df.columns = ["RIL", "Subpop"]
    group_df = group_df.groupby("Subpop")

    gpeObjs = GenePredObj(gpe_annotation)
    final_genes = []
    if gene_list:
        if os.path.isfile(gene_list):
            with open(gene_list) as f:
                for line in f.readlines():
                    final_genes.append(line.strip("\n").split("\t")[0])
        else:
            final_genes = gene_list
    else:
        if chrom and chr_start and chr_end and (split_by_gene == True):
            for strand in gpeObjs.genePredDict[chrom]:
                for transInfo in gpeObjs.genePredDict[chrom][strand]:
                    if transInfo[0] > chr_end: break
                    if transInfo[1] < chr_start: continue
                    if transInfo[0] >= chr_start and transInfo[1] <= chr_end:
                        final_genes.append(transInfo[2].geneName)

    temp_files = []
    if final_genes:
        for gene in final_genes:
            gpeList = gpeObjs.geneName2gpeObj(gene)
            gene_chrom = gpeList[0].chrom
            gene_start = min([x.txStart for x in gpeList])
            gene_end = max([x.txEnd for x in gpeList])
            gene_strand = gpeList[0].strand
            chrom_range = [gene_chrom, gene_start, gene_end]
            snp_outfile = "tmp_snp_outfile.gene_{}".format(gene)
            snp_outfile_vcf = extract_snp(plink_bed, snp_outfile, chrom_range)
            temp_files.append(snp_outfile_vcf)

            all_density_files = []
            for key, item in group_df:
                pop_strains_file = "{}.pop_strains".format(key)
                group_df.get_group(key).to_csv(pop_strains_file, columns=["Subpop"])
                temp_files.append(pop_strains_file)

                density_prefix = "{}.{}_density".format(key, gene)
                density_file = calc_n_diversity(snp_outfile_vcf, method=density_method, sample_file=pop_strains_file, out_path=density_prefix, nd_range=ex_range)
                all_density_files.append(density_file)
            density_data = merge_nd_file(all_density_files, method=density_method)



    elif not final_genes and chrom and chr_start and chr_end and (split_by_gene == False):
        snp_outfile = "tmp_snp_outfile.{}_{}_{}".format(chrom, chr_start, chr_end)
        chrom_range = [chrom, chr_start, chr_end]
        snp_outfile_vcf = extract_snp(plink_bed, snp_outfile, chrom_range)



