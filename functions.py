#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: functions.py
Author: CrazyHsu @ crazyhsu9527@gmail.com
Last modified: 2021-09-28 13:27:48
'''

import os, re, glob

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import rpy2.robjects as robjects
import warnings
warnings.filterwarnings("ignore")

from scipy.sparse import csr_matrix
from scipy.stats import expon
from scipy.interpolate import interp1d
from multiprocessing import Pool

############# Functions #############
def getOneColumnToList(targetFile, sep="\t", col=0):
    df = pd.read_csv(targetFile, sep=sep)
    return df.iloc[:, col].to_list()

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

def resolveDir(dirName, chdir=True):
    if not os.path.exists(dirName):
        os.makedirs(dirName)
    if chdir:
        os.chdir(dirName)

def parallel_plot(myFunc, values, threads=10):
    # for i in values:
    #     print(i)
    # try:
    #     p = Pool(processes=threads)
    #     results = [p.apply_async(myFunc, param) for param in values]
    #     p.close()
    #     p.join()
    # except KeyboardInterrupt as e:
    #     p.terminate()
    #     p.join()
    #     raise e
    # [print(*param, len(param)) for param in values]
    p = Pool(processes=threads)
    results = [p.apply_async(myFunc, param) for param in values]
    p.close()
    p.join()

    return [res.get() for res in results]

    # for i in values:
    #     myFunc(*i)

#####################################

############# Phylogenetic Tree ##############
def ped2fasta(ped_file, out_file):
    out = open(out_file, 'w')
    with open(ped_file) as f:
        for line in f.readlines():
            tmp = line.strip().split(' ')
            name = tmp[1]
            seq = ''.join(tmp[6:])
            seq = re.sub("0", "N", seq)
            out.write('>'+name+"\n"+seq+"\n")
    out.close()

def phylo_plot(plink_bed, output_prefix, group_file=None, group_sep=",", selected_lines=None, drop_line=None):
    command = 'plink --bfile {} --recode --out {} --silent'.format(plink_bed, output_prefix)
    if selected_lines:
        command = command + ' --keep {}'.format(selected_lines)
    os.system(command)
    ped_file = output_prefix + '.ped'
    fasta_file = output_prefix + '.fasta'
    ped2fasta(ped_file, fasta_file)
    nwk_file = '{}.tree.nwk'.format(output_prefix)
    command = 'FastTree -nt -gtr -quiet {} > {}'.format(fasta_file, nwk_file)
    # os.system(command)
    from plotWithR import drawPhyloPlotR
    robjects.r(drawPhyloPlotR)
    robjects.r.drawPhyloPlot(nwk_file, output_prefix, drop_line, group_file, group_sep)

##############################################


############# GWAS results plots #############
'''
### The format of GWAS file (separate can be automated detected):
rs p_wald
1_1922301 9.121183e-03
1_1928050 1.795902e-03
1_2521954 7.200593e-03
1_2522874 6.791745e-03
'''

def manhattan_plot(gwas_res, out_dir, thresholdi=None, gwas_sep="\t", data_from="file", threads=10, dpi=300):
    thresholdi = ",".join(map(str, thresholdi))
    file_path = os.path.dirname(os.path.abspath(__file__))
    from plotWithR import drawManhattanPlotR
    from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
    import logging
    rpy2_logger.setLevel(logging.ERROR)
    cmplot = os.path.join(file_path, "utils", "CMplot.R")
    robjects.r('source("' + cmplot + '")')
    robjects.r(drawManhattanPlotR)

    resolveDir(out_dir, chdir=False)
    if data_from == "file":
        file_name = os.path.basename(gwas_res)
        out_plot_prefix = os.path.join(out_dir, file_name + "_manhattan")
        robjects.r.drawManhattanPlot(gwas_res, thresholdi, out_plot_prefix, gwas_sep, "", dpi)
    elif data_from == "list":
        value_list = []
        with open(gwas_res) as f:
            for i in f.readlines():
                file_name = os.path.basename(i)
                out_plot_prefix = os.path.join(out_dir, file_name + "_manhattan")
                value_list.append((i.strip(), thresholdi, out_plot_prefix, gwas_sep, "", dpi))
        parallel_plot(robjects.r.drawManhattanPlot, value_list, threads=threads)
    elif data_from == "directory":
        value_list = []
        for i in glob.glob(os.path.join(gwas_res + '/*.assoc.txt')):
            file_name = os.path.basename(i)
            out_plot_prefix = os.path.join(out_dir, file_name + "_manhattan")
            value_list.append((i.strip(), thresholdi, out_plot_prefix, gwas_sep, "", dpi))
        parallel_plot(robjects.r.drawManhattanPlot, value_list, threads=threads)
    else:
        return "You should specify the correct type of data from, 'file' or 'list' or 'directory'"

def qq_plot(gwas_res, out_dir, gwas_sep="\t", data_from="file", threads=10, dpi=300, threshold=0.0001):
    file_path = os.path.dirname(os.path.abspath(__file__))
    from plotWithR import drawQQplotR
    from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
    import logging
    rpy2_logger.setLevel(logging.ERROR)
    cmplot = os.path.join(file_path, "utils", "CMplot.R")
    robjects.r('source("' + cmplot + '")')
    robjects.r(drawQQplotR)

    if data_from == "file":
        file_name = os.path.basename(gwas_res)
        out_plot_prefix = os.path.join(out_dir, file_name + "_qq")
        robjects.r.drawQQplot(gwas_res, out_plot_prefix, gwas_sep, "", dpi, threshold)
    elif data_from == "list":
        value_list = []
        with open(gwas_res) as f:
            for i in f.readlines():
                file_name = os.path.basename(i)
                out_plot_prefix = os.path.join(out_dir, file_name + "_qq")
                value_list.append((i.strip(), out_plot_prefix, gwas_sep, "", dpi, threshold))
                # robjects.r.drawQQplot(gwas_res, out_plot_prefix, gwas_sep)
        parallel_plot(robjects.r.drawQQplot, value_list, threads=threads)
    elif data_from == "directory":
        value_list = []
        for i in glob.glob(os.path.join(gwas_res + '/*.assoc.txt')):
            file_name = os.path.basename(i)
            out_plot_prefix = os.path.join(out_dir, file_name + "_qq")
            value_list.append((i.strip(), out_plot_prefix, gwas_sep, "", dpi, threshold))
            # robjects.r.drawQQplot(gwas_res, out_plot_prefix, gwas_sep)
        # print(robjects.r.drawQQplot, value_list, threads)
        parallel_plot(robjects.r.drawQQplot, value_list, threads=threads)
    else:
        return "You should specify the correct type of data from, 'file' or 'list' or 'directory'"

##############################################

############# MR results plots #############
''' the format of MR results
snp,mTrait,pTrait,effect,TMR,pvalue
'''

def forest_plot(mr_res, out_plot_prefix, order_by_name=False, order_by_effect=False,
                order_by_group=False, group=None, sep=",", height=6, width=4):
    from plotWithR import drawForestPlotR
    robjects.r(drawForestPlotR)
    robjects.r.drawForestPlot(mr_res, out_plot_prefix, str(order_by_name), str(order_by_effect), str(order_by_group),
                              group, sep, height, width)

def scattermr_plot(mr_res, out_plot_prefix, group_file=None, sep=",", height=5, width=9):
    from plotWithR import drawScatterMrPlotR
    if group_file != None:
        if not validateFile(group_file):
            group_file = None
    robjects.r(drawScatterMrPlotR)
    robjects.r.drawScatterMrPlot(mr_res, str(group_file), out_plot_prefix, sep, height, width)

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
def scatterps_plot(ps_file, group_file, ps_type, out_plot_prefix, input_sep=",", group_sep=","):
    if ps_type.lower() not in ["pca", "tsne"]:
        print("Please input invalid structure type such as 'pca' or 'tsne'. ")
        return
    from plotWithR import drawScatterPsR
    robjects.r(drawScatterPsR)
    robjects.r.drawScatterPs(ps_file, group_file, out_plot_prefix, ps_type.lower(), input_sep=input_sep, group_sep=group_sep)

    # if ps_type.lower() == "pca":
    #     from plotWithR import drawScatterPsR
    #     robjects.r(drawScatterPsR)
    #     robjects.r.drawScatterPs(ps_file, group_file, out_plot_prefix, input_sep=input_sep, group_sep=group_sep)
    # elif ps_type.lower == "tsne":
    #
    #     pass
    # else:
    #     raise Exception("Please input the correct poplation structure type you want to plot")
##################################################

############# MR network #############
def mr2net(mr_res, mr_sep, out_file):
    id_set = set()
    pair_dict = {}
    out = open(out_file, "w")
    with open(mr_res) as f:
        records = f.readlines()
        print(mr_sep.join(['mTrait', 'pTrait', 'pvalue']), file=out)
        for line in records[1:]:
            infoList = line.strip().split(mr_sep)
            key = "{}_{}".format(infoList[1], infoList[2])
            pair_dict[key] = float(infoList[5])
            id_set.add(infoList[1])
            id_set.add(infoList[2])
    import itertools
    for id_pair in itertools.product(list(id_set), list(id_set)):
        key = "{}_{}".format(id_pair[0], id_pair[1])
        if key in pair_dict:
            print(mr_sep.join(map(str, [id_pair[0], id_pair[1], pair_dict[key]])), file=out)
        else:
            print(mr_sep.join(map(str, [id_pair[0], id_pair[1], "1"])), file=out)
    out.close()

def get_weight(pvalue_list, pvalue):
    pvalue_list = pvalue_list[['mTrait', 'pTrait', 'pvalue']]
    pvalue_matrix = pd.pivot_table(pvalue_list, index='mTrait', columns='pTrait', values='pvalue', fill_value=1)
    id_list = list(set(pvalue_list['mTrait']) | set(pvalue_list['pTrait']))
    pvalue_matrix = pvalue_matrix.reindex(id_list)[id_list].fillna(1)
    pvalue_matrix = -np.log10(pvalue_matrix) + np.log10(pvalue)
    pvalue_matrix[pvalue_matrix < 0] = 0
    weight = 1 - expon.pdf(pvalue_matrix)
    weight = (weight + weight.T) / 2
    weight_csr = csr_matrix(np.triu(weight, k=1))
    weight_pair = list()
    gene_id = pvalue_matrix.columns
    for row in range(weight.shape[0]):
        for col, w in zip(weight_csr.indices[weight_csr.indptr[row]:weight_csr.indptr[row + 1]], weight_csr.data[weight_csr.indptr[row]:weight_csr.indptr[row + 1]]):
            weight_pair.append([gene_id[row], gene_id[col], w])
    weight_pair = pd.DataFrame(weight_pair)
    return weight_pair

def module_identify(edge_weight_fn, module_size):
    prefix = edge_weight_fn.replace('.edge_list', '')
    file_path = os.path.dirname(os.path.abspath(__file__))
    cl1_path = os.path.join(file_path, "utils", "cluster_one-1.0.jar")
    run_command = 'java -jar ' + cl1_path + ' -s ' + str(module_size) + ' -f edge_list -F csv ' + edge_weight_fn + ' >' + prefix + '.cluster_one.result.csv 2>/dev/null'
    os.system(run_command)
    cluster_one_res = pd.read_csv(prefix + '.cluster_one.result.csv')
    cluster_one_res = cluster_one_res.loc[cluster_one_res['P-value'] <= 0.05, :]
    cluster_one_res = cluster_one_res[['Size', 'Members']].sort_values(by='Size', ascending=False)
    cluster_one_res.loc[:, 'module'] = np.arange(1, cluster_one_res.shape[0] + 1)
    cluster_one_res.columns = ['gene_num', 'genes', 'module']
    cluster_one_res = cluster_one_res[['module', 'gene_num', 'genes']]
    return cluster_one_res

def hub_identify(edge_weight, cluster_one_res):
    robjects.r('''
        hub_ide <- function(m_ew){
            g <- igraph::graph.data.frame(m_ew,directed = FALSE)
            hub <- igraph::hub_score(g)
            hub_df <- as.data.frame(hub$vector,stringsAsFactors=F)
            names(hub_df) <- 'hub_score'
            hub_df$gene_id <- rownames(hub_df)
            rownames(hub_df) <- NULL
            return(hub_df)
        }
    ''')
    hub_ide = robjects.r('hub_ide')
    hub_list = list()
    hub_res = pd.DataFrame()
    for index, row in cluster_one_res.iterrows():
        gene_list = row['genes'].split(' ')
        m_ew = edge_weight.loc[(edge_weight['row'].isin(gene_list)) & (edge_weight['col'].isin(gene_list)), :]
        hub = hub_ide(m_ew)
        hub_res = pd.concat([hub_res, hub])
        hub_list.append(' '.join(hub.loc[hub.hub_score >= 0.8, :].sort_values(by='hub_score', ascending=False).apply(lambda x: x['gene_id']+'('+str(x['hub_score'])+')', axis=1).values))
    cluster_one_res.loc[:, 'hub_gene'] = hub_list
    return cluster_one_res, hub_res

def net_plot(input, out_prefix, input_sep, pvalue, module_size=5):
    from plotWithR import drawNetworkPlotR
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()
    mr2net_tmp = "{}.mr2net.tmp".format(out_prefix)
    mr2net(input, input_sep, mr2net_tmp)
    mr_res = pd.read_csv(mr2net_tmp, sep=input_sep)
    ew = get_weight(mr_res, pvalue)
    ew.columns = ['row', 'col', 'weight']
    ew.to_csv(out_prefix + '.edge_list', sep='\t', index=False, header=False)
    cluster_one_res = module_identify(out_prefix + '.edge_list', module_size)
    if cluster_one_res is None:
        raise Exception('cluster_one-1.0.jar file is not in /pwd/MRBIGR/utils, network analysis can not be run.')
    cluster_one_res, hub_res = hub_identify(ew, cluster_one_res)
    cluster_one_res.to_csv(out_prefix + '.module.csv', index=False)

    robjects.r(drawNetworkPlotR)
    for index, row in cluster_one_res.iterrows():
        gene_list = row['genes'].split(' ')
        m_ew = ew.loc[(ew['row'].isin(gene_list)) & (ew['col'].isin(gene_list)), :]
        m_hub = hub_res.loc[hub_res['gene_id'].isin(gene_list), :]
        robjects.r.drawNetworkPlot(m_ew, m_hub, out_prefix + '_module' + str(row['module']) + '.pdf')
    # module_network_plot(ew, cluster_one_res, hub_res, args.o)

###################################

############# GO plots ############
def go_plot(input_list, gene2go_file, bg_sep, plot_type="dotplot", sigp=0.05, adjustp=False, out_prefix="MODASv"):
    from plotWithR import drawEnrichPlotR
    from rpy2.rinterface_lib.callbacks import logger as rpy2_logger
    import logging
    rpy2_logger.setLevel(logging.ERROR)
    robjects.r(drawEnrichPlotR)
    robjects.r.drawEnrichPlot(input_list, gene2go_file, bg_sep, plot_type, sigp, str(adjustp), out_prefix)

###################################

############# Heatmap #############
def heatmap(input, out_prefix, select_rows=None, select_cols=None, cluster_rows=True, cluster_cols=True,
            order_col_by=None, order_row_by=None, anno_row=None, anno_col=None, scale_by="none",
            show_colnames=False, show_rownames=False, height=5, width=5):
    from plotWithR import drawHeatmapR
    [cluster_rows, cluster_cols, show_colnames, show_rownames] = ["T" if x else "F" for x in [cluster_rows, cluster_cols, show_colnames, show_rownames]]
    if order_col_by:
        if not os.path.isfile(order_col_by):
            order_col_by = "None"
    else:
        order_col_by = str(order_col_by)
    if order_row_by:
        if not os.path.isfile(order_row_by):
            order_row_by = "None"
    else:
        order_row_by = str(order_row_by)
    if anno_row:
        if not os.path.isfile(anno_row):
            anno_row = "None"
    else:
        anno_row = str(anno_row)
    if anno_col:
        if not os.path.isfile(anno_col):
            anno_col = "None"
    else:
        anno_col = str(anno_col)
    if select_rows:
        if not os.path.isfile(select_rows):
            select_rows = "None"
    else:
        select_rows = str(select_rows)
    if select_cols:
        if not os.path.isfile(select_cols):
            select_cols = "None"
    else:
        select_cols = str(select_cols)

    robjects.r(drawHeatmapR)
    robjects.r.drawHeatmap(input, out_prefix, select_rows, select_cols, cluster_rows, cluster_cols,
                           order_col_by, order_row_by, anno_row, anno_col, scale_by, show_rownames, show_colnames,
                           height, width)

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
def box_plot(input_mat_file, snp_file, group_file, selected_genes, selected_snps, selected_strains,
             scale, plot_type, out_prefix):
    if selected_genes != None:
        if os.path.isfile(selected_genes):
            selected_genes = ",".join(map(str, getOneColumnToList(selected_genes)))
        else:
            selected_genes = selected_genes
    else:
        selected_genes = "None"
    if selected_snps != None:
        if os.path.isfile(selected_snps):
            selected_snps = ",".join(map(str, getOneColumnToList(selected_snps)))
        else:
            selected_snps = selected_snps
    else:
        selected_snps = "None"
    if selected_strains != None:
        if os.path.isfile(selected_strains):
            selected_strains = ",".join(map(str, getOneColumnToList(selected_strains)))
        else:
            selected_strains = selected_strains
    else:
        selected_strains = "None"

    from plotWithR import drawBoxplotR
    robjects.r(drawBoxplotR)
    robjects.r.drawBoxplot(input_mat_file, snp_file, group_file, selected_genes, selected_snps,
                           selected_strains, scale, plot_type, out_prefix)
###################################

############# Nucleotide diversity Plots #############
def getLongestTrans(gpeList):
    tmp_exon_length = 0
    longestTrans = None
    for gpe in gpeList:
        if gpe.exon_length() > tmp_exon_length:
            tmp_exon_length = gpe.exon_length()
            longestTrans = gpe
    return longestTrans

def extract_snp(plink_bed, snp_outfile, selected_strains=None, chrom_range=None, plink_path=None):
    if not plink_path:
        plink_path = 'plink'
    if chrom_range:
        command = plink_path + ' -bfile ' + plink_bed + ' --out ' + snp_outfile + \
                  ' --chr ' + str(chrom_range[0]) + ' --from-bp ' + str(chrom_range[1]) + \
                  ' --to-bp ' + str(chrom_range[2]) + ' --recode vcf-iid --silent'
    else:
        command = plink_path + ' -bfile ' + plink_bed + ' --out ' + snp_outfile + ' --recode vcf-iid'
    if selected_strains:
        command = command + ' --keep {}'.format(selected_strains)
    os.system(command)
    # os.system('rm ' + snp_outfile + '.log')
    # os.system('rm ' + snp_outfile + '.nosex')
    return snp_outfile + ".vcf"

def merge_nd_file(nd_path, method="pi", group_name=None, window_size=50):
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
            nd_data = nd_data.rename(columns={'PI': group_name[i]})
        else:
            tmp_data = data.loc[:, merge_columns]
            tmp_data = tmp_data.rename(columns={'PI': group_name[i]})
            nd_data = nd_data.merge(tmp_data, on='BIN_START', how="outer")

    pos_list = pd.DataFrame({"BIN_START": range(min(nd_data["BIN_START"]), max(nd_data["BIN_START"]), window_size)})
    nd_data = pos_list.merge(nd_data, on='BIN_START', how="outer")
    nd_data = nd_data.fillna(0)
    return nd_data

def calc_n_diversity(snp_file, method="pi", vcftools_path=None, selected_lines=None, out_path=None, nd_range=None,
                     fst_sample=None, window_size=50, window_step=50):
    if not vcftools_path:
        vcftools_path = 'vcftools'
    if not out_path:
        dir_path = '/'.join(snp_file.split('/')[:-1])
        file_name = '.'.join(snp_file.split('/')[-1].split('.')[:-1])
        out_path = dir_path + '/' + file_name

    command = vcftools_path + ' --vcf ' + snp_file + ' --out ' + out_path
    if selected_lines:
        command = command + ' --keep ' + selected_lines
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
    os.system(command + ' 2>/dev/null')
    return out_path + suffix

def line_smooth(diversity_data):
    new_df = pd.DataFrame()
    new_index = np.arange(min(diversity_data.index), max(diversity_data.index))
    for name in diversity_data.columns:
        new_field = interp1d(diversity_data.index, diversity_data[name], kind="cubic")
        new_df[name] = new_field(new_index)
    new_df.index = new_index
    new_df = new_df.clip(lower=0)
    return new_df

def plot_exon(ax, longestTrans, plot_length, strand, max_y_value):
    import itertools
    exon_list = longestTrans.exons
    exon_height = max_y_value / 6.0
    left_pos = min(itertools.chain(*exon_list))
    right_pos = max(itertools.chain(*exon_list))
    middle_pos = (left_pos+right_pos)/2
    for exon in exon_list:
        if strand == "+":
            exon_start = exon[0]
            exon_length = exon[1] - exon[0]
            fill_color = "#CC79A7"
            edge_color = "#774761"
        else:
            exon_start = exon[1]
            exon_length = exon[0] - exon[1]
            fill_color = "#0072B2"
            edge_color = "#003c5d"
        head_length = min(abs(exon_length) / 2, plot_length * 0.01)
        ax.arrow(exon_start, -exon_height, exon_length, 0, width=exon_height/1.5, head_width=exon_height/1.5,
                 head_length=head_length, shape="full", length_includes_head=True, facecolor=fill_color,
                 edgecolor=edge_color)
    annotation_height = -exon_height/2
    ax.annotate(longestTrans.transName, xy=(middle_pos, 0), xytext=(middle_pos, annotation_height),
                fontsize=12, ha="center")
    return exon_height

def plot_edge(ax, longestTrans, exon_height, strand):
    intron_list = longestTrans.introns
    for intron in intron_list:
        ax.plot(intron, [-exon_height, -exon_height], linestyle="-", color="grey")

def plot_diversity(diversity_data, out_prefix, suffix="pdf", smooth=False, plot_left_lim=0, plot_right_lim=0,
                   longestTrans=None, x_label="", y_label=""):
    if smooth:
        diversity_data = line_smooth(diversity_data)
    fig, ax = plt.subplots(figsize=(12, 9))
    diversity_data.plot(ax=ax, kind="line")
    max_y_value = max(np.max(diversity_data))
    if longestTrans:
        exon_height = plot_exon(ax, longestTrans, plot_right_lim-plot_left_lim, "-", max_y_value)
        plot_edge(ax, longestTrans, exon_height, "-")

    out_file = "{}.{}".format(out_prefix, suffix)
    plt.ticklabel_format(useOffset=False, style='plain')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tick_params(axis="both", length=5)

    y_tick_labels = ax.get_yticks().tolist()
    pasitive_label = [x for i, x in enumerate(y_tick_labels) if x >= 0]
    ax.set_yticks(pasitive_label)

    plt.xlim([plot_left_lim, plot_right_lim])
    plt.xlabel(x_label, fontsize=16)
    plt.ylabel(y_label, fontsize=16)
    plt.savefig(out_file)
    # plt.figure(figsize=(10, 10))
    # plt.plot(diversity_data["loci"].values, diversity_data['Mixed'].values)
    # plt.ticklabel_format(useOffset=False, style='plain')
    # plt.xticks(fontsize=16)
    # plt.yticks(fontsize=16)
    # plt.tick_params(axis="both", length=5)
    # out_file = "{}.{}".format(out_prefix, suffix)
    # plt.savefig(out_file)


    # ax = diversity_data.plot(kind="bar")
    # fig = ax.get_figure()
    # fig.set_size_inches(10, 10)
    # # ax.ticklabel_format(useOffset=False, style='plain')
    # # ax.xticks(fontsize=16)
    # # ax.yticks(fontsize=16)
    # # ax.tick_params(axis="both", length=5)
    # out_file = "{}.{}".format(out_prefix, suffix)
    # fig.savefig(out_file)

# def nucleotide_dst_plot(plink_bed, group_file, gff_annotation, gene_list=None, chrom=None, chr_start=0, chr_end=0,
#                         diversity_method="pi", plot_by_gene=True, keep_temp=True, out_dir="nd_output"):
#     from commonClasses import GenePredExtLine, GenePredObj
#     # from dna_features_viewer import GraphicFeature, GraphicRecord
#
#     exon_color = "#F19800"
#     gpe_annotation = ""
#     if gff_annotation.endswith("gpe"):
#         gpe_annotation = gff_annotation
#     elif gff_annotation.endswith("gtf"):
#         gpe_annotation = ".".join([os.path.splitext(gff_annotation)[0], "gpe"])
#         command = "gtfToGenePred {} {} -genePredExt".format(gff_annotation, gpe_annotation)
#         os.system(command)
#     elif gff_annotation.endswith("gff") or gff_annotation.endswith("gff3"):
#         gpe_annotation = ".".join([os.path.splitext(gff_annotation)[0], "gpe"])
#         command = "gffread -T -o tmp.gtf {} ".format(gff_annotation)
#         os.system(command)
#         command = "gtfToGenePred {} {} -genePredExt".format(gff_annotation, gpe_annotation)
#         os.system(command)
#         os.system('rm tmp.gtf')
#     else:
#         raise Exception("Please input the annotation with suffix including: gpe, gff3, gff, gtf.")
#
#     group_df = pd.read_csv(group_file, sep=",")
#     group_df = group_df.iloc[:, [0, 1]]
#     group_df.columns = ["RIL", "Subpop"]
#     group_df = group_df.groupby("Subpop")
#
#     gpeObjs = GenePredObj(gpe_annotation)
#     final_genes = []
#     if gene_list:
#         if os.path.isfile(gene_list):
#             with open(gene_list) as f:
#                 for line in f.readlines():
#                     final_genes.append(line.strip("\n").split("\t")[0])
#         else:
#             final_genes = gene_list
#     else:
#         if chrom and chr_start and chr_end and plot_by_gene:
#             for strand in gpeObjs.genePredDict[chrom]:
#                 for transInfo in gpeObjs.genePredDict[chrom][strand]:
#                     if transInfo[0] > chr_end: break
#                     if transInfo[1] < chr_start: continue
#                     if transInfo[0] >= chr_start and transInfo[1] <= chr_end:
#                         final_genes.append(transInfo[2].geneName)
#
#     temp_files = []
#     if final_genes:
#         for gene in final_genes:
#             gpeList = gpeObjs.geneName2gpeObj(gene)
#             gene_chrom = gpeList[0].chrom
#             gene_start = min([x.txStart for x in gpeList])
#             gene_end = max([x.txEnd for x in gpeList])
#             gene_strand = gpeList[0].strand
#             chrom_range = [gene_chrom, gene_start, gene_end]
#             snp_outfile = "tmp_snp_outfile.gene_{}".format(gene)
#             snp_outfile_vcf = extract_snp(plink_bed, snp_outfile, chrom_range=chrom_range)
#             temp_files.append(snp_outfile_vcf)
#
#             all_density_files = []
#             for key, item in group_df:
#                 pop_strains_file = "{}.pop_strains".format(key)
#                 group_df.get_group(key).to_csv(pop_strains_file, columns=["Subpop"])
#                 temp_files.append(pop_strains_file)
#
#                 density_prefix = "{}.{}_density".format(key, gene)
#                 density_file = calc_n_diversity(snp_outfile_vcf, method=diversity_method, selected_lines=pop_strains_file, out_path=density_prefix, nd_range=ex_range)
#                 all_density_files.append(density_file)
#             density_data = merge_nd_file(all_density_files, method=diversity_method)
#
#
#
#     elif not final_genes and chrom and chr_start and chr_end and not plot_by_gene:
#         snp_outfile = "tmp_snp_outfile.{}_{}_{}".format(chrom, chr_start, chr_end)
#         chrom_range = [chrom, chr_start, chr_end]
#         snp_outfile_vcf = extract_snp(plink_bed, snp_outfile, chrom_range=chrom_range)

def nucleotide_dst_plot1(plink_bed, group_file, gff_annotation, gene_list="", chrom=None, chr_start=0, chr_end=0,
                         diversity_method="pi", plot_by_gene=True, keep_temp=True, out_dir="MODASv", out_prefix="MODASv",
                         left_offset=0, right_offset=0, window_size=1000, window_step=50, smooth=False,
                         plot_left_expand=500, plot_right_expand=500):
    from commonClasses import GenePredExtLine, GenePredObj
    # from dna_features_viewer import GraphicFeature, GraphicRecord
    resolveDir(out_dir, chdir=False)

    # exon_color = "#F19800"
    # gpe_annotation = ""
    if gff_annotation.endswith("gpe"):
        gpe_annotation = gff_annotation
    elif gff_annotation.endswith("gtf"):
        gpe_annotation = ".".join([os.path.splitext(os.path.basename(gff_annotation))[0], "gpe"])
        gpe_annotation = os.path.join(out_dir, gpe_annotation)
        command = "gtfToGenePred {} {} -genePredExt".format(gff_annotation, gpe_annotation)
        # os.system(command)
    elif gff_annotation.endswith("gff") or gff_annotation.endswith("gff3"):
        gpe_annotation = ".".join([os.path.splitext(os.path.basename(gff_annotation))[0], "gpe"])
        gpe_annotation = os.path.join(out_dir, gpe_annotation)
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
            gene_list = list(map(str.strip, gene_list.strip().split(",")))
            if gene_list:
                final_genes = gene_list
    else:
        if chrom and chr_start and chr_end and not plot_by_gene:
            for strand in gpeObjs.genePredDict[chrom]:
                for transInfo in gpeObjs.genePredDict[chrom][strand]:
                    if transInfo[0] > chr_end: break
                    if transInfo[1] < chr_start: continue
                    if transInfo[0] >= chr_start and transInfo[1] <= chr_end:
                        final_genes.append(transInfo[2].geneName)
    temp_files = []
    if final_genes:
        for gene in final_genes:
            gpeList = gpeObjs.geneName2gpeObj[gene]
            longestTrans = getLongestTrans(gpeList)

            gene_chrom = gpeList[0].chrom
            gene_start = min([x.txStart for x in gpeList]) - left_offset
            gene_end = max([x.txEnd for x in gpeList]) + right_offset
            gene_strand = gpeList[0].strand
            chrom_range = [gene_chrom, gene_start, gene_end]
            snp_outfile = os.path.join(out_dir, "tmp_snp_outfile.gene_{}".format(gene))
            snp_outfile_vcf = extract_snp(plink_bed, snp_outfile, chrom_range=chrom_range)
            temp_files.append(snp_outfile_vcf)

            all_diversity_files = []
            group_name = []
            for key, item in group_df:
                pop_strains_file = os.path.join(out_dir, "{}.pop_strains".format(key))
                group_df.get_group(key).to_csv(pop_strains_file, columns=["RIL"], index=False, header=False)
                temp_files.append(pop_strains_file)

                diversity_prefix = os.path.join(out_dir, "{}.{}_diversity".format(key, gene))
                diversity_file = calc_n_diversity(snp_outfile_vcf, method=diversity_method,
                                                  selected_lines=pop_strains_file, out_path=diversity_prefix,
                                                  nd_range=chrom_range, window_size=window_size, window_step=window_step)
                all_diversity_files.append(diversity_file)
                group_name.append(key)
            diversity_data = merge_nd_file(all_diversity_files, method=diversity_method, group_name=group_name,
                                           window_size=window_size)
            diversity_data = diversity_data.rename(columns={"BIN_START": "loci", "N_VARIANTS": "n_variants"})
            diversity_data = diversity_data.set_index("loci")
            diversity_data = diversity_data.sort_index()
            plot_left_lim = chr_start - plot_left_expand
            plot_right_lim = chr_end + plot_right_expand
            x_label = "Chromosome {}".format(chrom.lower().split("chr")[1] if "chr" in chrom.lower() else chrom)
            y_label = "Nucleotide diversity"
            plot_diversity(diversity_data.iloc[:, 1:], os.path.join(out_dir, out_prefix + "_diversity"), smooth=smooth,
                           plot_left_lim=plot_left_lim, plot_right_lim=plot_right_lim, longestTrans=longestTrans,
                           x_label=x_label, y_label=y_label)


    elif not final_genes and chrom and chr_start and chr_end:
        snp_outfile = "tmp_snp_outfile.{}_{}_{}".format(chrom, chr_start, chr_end)
        snp_outfile = os.path.join(out_dir, snp_outfile)
        chrom_range = [chrom, chr_start, chr_end]
        snp_outfile_vcf = extract_snp(plink_bed, snp_outfile, chrom_range=chrom_range)
        temp_files.append(snp_outfile_vcf)
        all_diversity_files = []
        group_name = []
        for key, item in group_df:
            pop_strains_file = os.path.join(out_dir, "{}.pop_strains".format(key))
            group_df.get_group(key).to_csv(pop_strains_file, columns=["RIL"], index=False, header=False)
            temp_files.append(pop_strains_file)

            diversity_prefix = os.path.join(out_dir, "{}.{}_{}_{}_diversity".format(key, chrom, chr_start, chr_end))
            diversity_file = calc_n_diversity(snp_outfile_vcf, method=diversity_method,
                                              selected_lines=pop_strains_file, out_path=diversity_prefix,
                                              nd_range=chrom_range, window_size=window_size, window_step=window_step)
            all_diversity_files.append(diversity_file)
            group_name.append(key)
        diversity_data = merge_nd_file(all_diversity_files, method=diversity_method, group_name=group_name,
                                       window_size=window_size)
        diversity_data = diversity_data.rename(columns={"BIN_START": "loci", "N_VARIANTS": "n_variants"})
        diversity_data = diversity_data.set_index("loci")
        diversity_data = diversity_data.sort_index()
        plot_left_lim = chr_start - plot_left_expand
        plot_right_lim = chr_end + plot_right_expand
        plot_diversity(diversity_data.iloc[:, 1:], os.path.join(out_dir, out_prefix + "_diversity"), smooth=smooth,
                       plot_left_lim=plot_left_lim, plot_right_lim=plot_right_lim)

