#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: plotPheBox.py
Author: CrazyHsu @ crazyhsu9627@gmail.com
Created on: 2021-12-09 22:01:41
Last modified: 2021-12-09 22:01:41
'''
import argparse, sys, os
import pandas as pd
import rpy2.robjects as robjects

def construct_snp_genotype(bed_file, snps, group_file=None, group_sep=","):
    ped_file = "temp.ped"
    map_file = "temp.map"
    if snps == "None":
        command = "plink --bfile {} --thin 0.00001 --recode --out temp --tab --silent --allow-no-sex"
    else:
        command = "plink --bfile {} --snps {} --recode --out temp --tab --silent --allow-no-sex".format(bed_file, snps)
    os.system(command)
    ped_df = pd.read_csv(ped_file, sep="\t", header=None, dtype=str, index_col=0)
    ped_df = ped_df.applymap(lambda x: x.replace(" ", "") if isinstance(x, str) else x)
    map_df = pd.read_csv(map_file, sep="\t", header=None, dtype=str)

    if snps == "None":
        ped_df = ped_df.iloc[:, [5]]
        ped_df.columns = list(map_df.iloc[:, 1][0])
    else:
        ped_df = ped_df.iloc[:, 5:]
        ped_df.columns = map_df.iloc[:, 1]

    if group_file and os.path.exists(group_file):
        group = pd.read_csv(group_file, sep=group_sep, header=None, dtype=str, index_col=0)
        group.columns = ["group"]
        merged_df = pd.merge(group, ped_df, left_index=True, right_index=True, how="inner")
        ped_df = merged_df
    ped_df.to_csv("tmp_merged.txt", index=True, header=True)
    temp_files = ["temp.ped", "temp.map", "temp.log", "temp.nosex", "tmp_merged.txt"]
    return "tmp_merged.txt", temp_files

drawBoxplotR = '''
    suppressMessages(library(ggplot2))
    suppressMessages(library(data.table))
    suppressMessages(library(ggpubr))
    suppressMessages(library(rstatix))
    suppressMessages(library(R.utils))
    # options(warn=-1)
    removeHeterozygous <- function(df, sep="", col=1){
        g2 <- as.matrix(df[,col,drop=F])
        n <- nrow(g2)
        m <- ncol(g2)
        a1<-matrix(sapply(strsplit(g2,sep),"[",1),nrow = n,ncol=m,byrow = F)
        a2<-matrix(sapply(strsplit(g2,sep),"[",2),nrow = n,ncol=m,byrow = F)
        colnames(a1) <- colnames(g2)
        rownames(a1) <- rownames(g2)
        H<-matrix(as.numeric(!a1==a2),nrow = n,ncol = m,byrow = F)
        H <- as.data.frame(H)
        # colnames(H) <- colnames(df)
        # rownames(H) <- rownames(df)
        return(df[which(H[,1]==0), ,drop=F])
    }

    drawGenoPheBoxplot <- function(input_mat_file, snp_file, group_file, selected_genes, selected_snps, selected_strains, 
        to_scale, rmHeter, plot_type, out_prefix){
        input_mat <- fread(input_mat_file, data.table=getOption("datatable.fread.datatable", FALSE), header=TRUE)
        snp_mat <- fread(snp_file, data.table=getOption("datatable.fread.datatable", FALSE), header=TRUE)
        group_mat <- fread(group_file, data.table=getOption("datatable.fread.datatable", FALSE), header=TRUE)

        rownames(input_mat) <- input_mat[,1]
        input_mat[,1] <- NULL
        rownames(snp_mat) <- snp_mat[,1]
        snp_mat[,1] <- NULL
        rownames(group_mat) <- group_mat[,1]
        group_mat[,1] <- NULL

        if(selected_genes!="None"){
            selected_genes <- trimws(unlist(strsplit(selected_genes, ",")))
        }else{
            selected_genes <- c(colnames(input_mat)[1])
        }
        if(selected_snps!="None"){
            selected_snps <- trimws(unlist(strsplit(selected_snps, ",")))
        }else{
            selected_snps <- c(colnames(snp_mat)[1])
        }
        if(selected_strains!="None"){
            selected_strains <- trimws(unlist(strsplit(selected_strains, ",")))
        }else{
            selected_strains <- rownames(input_mat)
        }

        sub_input_mat <- input_mat[selected_strains, selected_genes, drop=F]
        sub_rownames <- rownames(sub_input_mat)
        sub_input_mat <- data.frame(lapply(sub_input_mat, function(x) {gsub("-", "0", x)}))
        rownames(sub_input_mat) <- sub_rownames
        sub_snp_mat <- snp_mat[selected_strains, selected_snps, drop=F]
        sub_group_mat <- group_mat[selected_strains, , drop=F]
        colnames(sub_group_mat) <- c("Subpop")

        for(snp in selected_snps){
            genotype <- sub_snp_mat[, snp, drop=F]
            if(rmHeter=="True"){
                genotype <- removeHeterozygous(genotype, col=1)
            }
            colnames(genotype) <- c("Genotype")
            if(plot_type=="single_gene"){
                for(gene in selected_genes){
                    single_gene_mat <- sub_input_mat[, gene, drop=F]
                    merged_exp <- merge(x=genotype, y=single_gene_mat, by=0)
                    merged_exp <- na.omit(merged_exp)
                    colnames(merged_exp) <- c("RIL", "Genotype", "Expression")

                    merged_exp$Genotype <- as.factor(merged_exp$Genotype)
                    merged_exp$Expression <- as.numeric(as.vector(merged_exp$Expression))
                    merged_exp <- na.omit(merged_exp)
                    if(to_scale=="True"){
                        merged_exp$Expression <- (merged_exp$Expression-mean(merged_exp$Expression))/sd(merged_exp$Expression)
                    }
                    p <- ggboxplot(data=merged_exp, x="Genotype", y="Expression", fill="Genotype", width=0.5, bxp.errorbar=T) +
                        theme_bw() +
                        theme(axis.text.x = element_text(colour = "black", size = 12, vjust =1 ),
                            axis.text.y = element_text(colour = "black", size = 12, hjust =1 ),
                            axis.title = element_text(margin=margin(10, 5, 0, 0), color = "black", size = 14),
                            axis.title.y = element_text(angle=90),
                            panel.border = element_rect(color = 'black', fill=NA, size = 1),
                            panel.grid.major = element_blank(), panel.grid.minor = element_blank()
                        ) + xlab('')
                    print("1")

                    my_comparisons <- unlist(lapply(2, combn, x = as.vector(unique(merged_exp$Genotype)), simplify = FALSE), recursive = F)
                    print("2")
                    p <- p + stat_compare_means(method="t.test", comparisons = my_comparisons, label="p", tip.length=0.02) +
                        scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
                    ggsave(paste(out_prefix, gene, 'boxplot.pdf', sep='_'), plot=p, width = 3, height = 3, useDingbats=FALSE)
                    print("3")
                }
            }else if(plot_type=="multi_genes") {
                gene_n <- length(selected_genes)
                if(gene_n==1){
                    width <- 3
                }else if(gene_n==2){
                    width <- 3.5
                }else{
                    width <- 3.5 + 0.5 * (gene_n-2)
                }
                if(width >= 8){
                    width <- 8
                }

                merged_exp <- merge(x=genotype, y=sub_input_mat, by=0)
                merged_exp <- na.omit(merged_exp)
                colnames(merged_exp)[1:2] <- c("RIL", "Genotype")
                if(to_scale=="True"){
                    tmp_exp <- merged_exp[,3:dim(merged_exp)[2]]
                    tmp_exp <- lapply(tmp_exp, function(x) as.numeric(as.character(x)))
                    merged_exp[,3:dim(merged_exp)[2]] <- sapply(tmp_exp, function(tmp_exp) (tmp_exp-mean(tmp_exp))/sd(tmp_exp))
                }

                setDT(merged_exp)
                melt_merged_exp <- melt(merged_exp, id.vars = c("RIL", "Genotype"), variable.name = "Gene", value.name = "Expression")                
                melt_merged_exp$Gene <- as.factor(melt_merged_exp$Gene)
                melt_merged_exp$Genotype <- as.factor(melt_merged_exp$Genotype)
                melt_merged_exp$Expression <- as.numeric(as.vector(melt_merged_exp$Expression))
                p <- ggboxplot(data=melt_merged_exp, x="Gene", y="Expression", fill="Genotype", width=0.5, bxp.errorbar=T) + 
                    theme_bw() +
                    theme(axis.text.x = element_text(colour = "black", size = 8, vjust =1, hjust = 1, angle=45),
                        axis.text.y = element_text(colour = "black", size = 12, hjust =1 ),
                        axis.title = element_text(margin=margin(10, 5, 0, 0), color = "black", size = 14),
                        axis.title.y = element_text(angle=90),
                        panel.border = element_rect(color = 'black', fill=NA, size = 1),
                        panel.grid.major = element_blank(), panel.grid.minor = element_blank()
                    ) + xlab('')
                stat.test <- melt_merged_exp %>% group_by(Gene) %>% t_test(Expression ~ Genotype)
                stat.test <- stat.test %>% add_xy_position(x = "Gene", dodge = 0.8)
                p <- p + stat_pvalue_manual(stat.test, label = "p", tip.length = 0.02) + 
                    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
                ggsave(paste(out_prefix, 'multiGenes', snp, 'boxplot.pdf', sep='_'), plot=p, width = width, height = 3, useDingbats=FALSE)
            }else if(plot_type=="single_pop") {
                for(gene in selected_genes){
                    single_gene_mat <- sub_input_mat[, gene, drop=F]
                    colnames(single_gene_mat) <- c("Expression")
                    sub_group_mat_tmp <- sub_group_mat
                    genotype_tmp <- genotype
                    single_gene_mat_tmp <- single_gene_mat

                    sub_group_mat_tmp["RIL"] <- rownames(sub_group_mat_tmp)
                    genotype_tmp["RIL"] <- rownames(genotype_tmp)
                    single_gene_mat_tmp["RIL"] <- rownames(single_gene_mat_tmp)
                    df_lst <- list(sub_group_mat_tmp, genotype_tmp, single_gene_mat_tmp)

                    merged_exp <- Reduce(function(x,y) merge(x,y,by="RIL"), df_lst)
                    merged_exp <- na.omit(merged_exp)
                    if(to_scale=="True"){
                        merged_exp$Expression <- (merged_exp$Expression-mean(merged_exp$Expression))/sd(merged_exp$Expression)
                    }
                    merged_exp = removeHeterozygous(merged_exp, col="Genotype")
                    merged_exp$Subpop <- as.factor(merged_exp$Subpop)
                    merged_exp$Genotype <- as.factor(merged_exp$Genotype)
                    merged_exp$Expression <- as.numeric(as.vector(merged_exp$Expression))
                    p <- ggboxplot(data=merged_exp, x="Subpop", y="Expression", fill="Genotype", width=0.5, 
                        bxp.errorbar=T, outlier.shape=20)+
                        theme_bw() +
                        theme(axis.text.x = element_text(colour = "black", size = 12, vjust =1 ),
                            axis.text.y = element_text(colour = "black", size = 12, hjust =1 ),
                            axis.title = element_text(margin=margin(10, 5, 0, 0), color = "black", size = 14),
                            axis.title.y = element_text(angle=90),
                            panel.border = element_rect(color = 'black', fill=NA, size = 1),
                            panel.grid.major = element_blank(), panel.grid.minor = element_blank()
                        ) + xlab('') + coord_cartesian(ylim = c(0, NA)) 
                    stat.test <- merged_exp %>% group_by(Subpop) %>% t_test(Expression ~ Genotype)
                    stat.test <- stat.test %>% add_xy_position(x = "Subpop", dodge = 0.8)
                    p <- p + stat_pvalue_manual(stat.test, label = "p", tip.length = 0.02, size=3) + 
                        scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
                    ggsave(paste(out_prefix, gene, snp, 'pop.boxplot.pdf', sep='_'), plot=p, width = 3, height = 3, useDingbats=FALSE)
                }
            }
        }
    }
    
    drawPlainPheBoxplot <- function(input_mat_file, group_file, selected_phes, selected_strains, scale, plot_type, out_prefix, height=5, width=8){
        input_mat <- fread(input_mat_file, data.table=getOption("datatable.fread.datatable", FALSE), header=TRUE)
        group_mat <- fread(group_file, data.table=getOption("datatable.fread.datatable", FALSE), header=TRUE)
        
        rownames(input_mat) <- input_mat[,1]
        input_mat[,1] <- NULL
        rownames(group_mat) <- group_mat[,1]
        group_mat[,1] <- NULL
        
        if(selected_phes!="None"){
            selected_phes <- trimws(unlist(strsplit(selected_phes, ",")))
        }else{
            selected_phes <- colnames(input_mat)
        }
        
        if(selected_strains!="None"){
            selected_strains <- trimws(unlist(strsplit(selected_strains, ",")))
        }else{
            selected_strains <- rownames(input_mat)
        }
        
        sub_input_mat <- input_mat[selected_strains, selected_phes, drop=F]
        sub_group_mat <- group_mat[selected_strains, , drop=F]
        
        colnames(group_mat) <- c("group")
        merged_mat <- merge(input_mat, group_mat, by=0)
        colnames(merged_mat)[1] <- c("RIL")
        
        merged_mat_melt <- melt(merged_mat, id.vars = c("RIL", "group"), variable.name = "phe", value.name = "value")
        
        p <- ggboxplot(data=merged_mat_melt, x="phe", y="value", fill="group", width=0.7, bxp.errorbar=T) + 
            theme_bw() + 
            theme(axis.text.x = element_text(colour = "black", size = 12, vjust =1, angle=45, hjust=1),
                axis.text.y = element_text(colour = "black", size = 12, hjust =1 ),
                axis.title = element_text(margin=margin(10, 5, 0, 0), color = "black", size = 14),
                axis.title.y = element_text(angle=90),
                panel.border = element_rect(color = 'black', fill=NA, size = 1),
                panel.grid.major = element_blank(), panel.grid.minor = element_blank()
            ) + xlab('')
        merged_mat_melt$value <- as.numeric(merged_mat_melt$value)
        stat.test <- merged_mat_melt %>% group_by(phe) %>% t_test(value ~ group)
        stat.test <- stat.test %>% add_xy_position(x = "phe", dodge = 0.8)
        p <- p + stat_pvalue_manual(stat.test, label = "p", tip.length = 0.02) + 
            scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
        ggsave(file=paste0(out_prefix, ".group_boxplot.pdf"), plot=p, width = width, height = height, useDingbats=FALSE)
    }
'''

def getOneColumnToList(targetFile, sep="\t", col=0):
    df = pd.read_csv(targetFile, sep=sep)
    return df.iloc[:, col].to_list()

def removeFiles(myDir=None, fileList=None):
    for f in fileList:
        if myDir:
            os.remove(os.path.join(myDir, f.strip("\n")))
        else:
            os.remove(f.strip("\n"))


def genoPhebox(input_mat_file, bed_file, group_file, selected_genes, selected_snps, selected_strains,
               scale, rmHeter, plot_type, out_prefix):
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

    snp_file, temp_files = construct_snp_genotype(bed_file, selected_snps)
    # from plotWithR import drawBoxplotR
    robjects.r(drawBoxplotR)
    robjects.r.drawGenoPheBoxplot(input_mat_file, snp_file, group_file, selected_genes, selected_snps,
                                  selected_strains, scale, rmHeter, plot_type, out_prefix)
    removeFiles(fileList=temp_files)

def plainPheBox(input_mat_file, group_file, selected_phes, selected_strains, scale, plot_type, out_prefix,
                height, width):
    if selected_phes:
        if os.path.isfile(selected_phes):
            selected_phes = ",".join(map(str, getOneColumnToList(selected_phes)))
        else:
            selected_phes = selected_phes
    else:
        selected_phes = "None"
    if selected_strains:
        if os.path.isfile(selected_strains):
            selected_strains = ",".join(map(str, getOneColumnToList(selected_strains)))
        else:
            selected_strains = selected_strains
    else:
        selected_strains = "None"

    robjects.r(drawBoxplotR)
    robjects.r.drawPlainPheBoxplot(input_mat_file, group_file, selected_phes, selected_strains, scale,
                                   plot_type, out_prefix, height=height, width=width)

def main(args):
    if args.command == "plot_genoPhebox":
        genoPhebox(args.input, args.plink_bed, args.group, args.selected_genes, args.selected_snps,
                   args.selected_strains, str(args.scale), str(args.rmHeter), args.plot_type, args.out_prefix)
    if args.command == "plot_plainPhebox":
        plainPheBox(args.input, args.group, args.selected_phes, args.selected_strains, str(args.scale),
                    args.plot_type, args.out_prefix, args.height, args.width)

if __name__ == '__main__':
    USAGE = ' Miscellaneous plot functions to visualize the results of MODAS. '

    parser = argparse.ArgumentParser(usage='%(prog)s command [options]', description=USAGE)
    subparsers = parser.add_subparsers(title='command', metavar='', dest='command', prog=parser.prog)

    plot_gpb = subparsers.add_parser('plot_genoPhebox', help='Draw genotype-dependent boxline plots for the phenotypes in different strains',
                                     usage='%(prog)s [options]')
    plot_gpb.add_argument('-i', dest="input", help="The expression file")
    plot_gpb.add_argument('-bed', dest="plink_bed", help="The plink bed file")
    plot_gpb.add_argument('-group', dest="group", help="The group file")
    plot_gpb.add_argument('-s_genes', dest="selected_genes", default=None,
                          help="The file listing selected genes or a string with comma separated")
    plot_gpb.add_argument('-s_snps', dest="selected_snps", default=None,
                          help="The file listing selected SNPs or a string with comma separated")
    plot_gpb.add_argument('-s_strains', dest="selected_strains", default=None,
                          help="The file listing selected strains or a string with comma separated")
    plot_gpb.add_argument('-scale', dest="scale", default=False, action="store_true", help="Scale the phenotype values")
    plot_gpb.add_argument('-rmHeter', dest="rmHeter", default=False, action="store_true",
                          help="Remove heterozygous genotypes. Default: False")
    plot_gpb.add_argument('-plot_type', dest="plot_type", choices=["single_gene", "multi_genes", "single_pop"],
                          default="single_gene", help="The type of graph to be plotted. "
                                                      "Options are: single_gene, multi_gene, single_pop. "
                                                      "Default: single_gene")
    plot_gpb.add_argument("-o", dest="out_prefix", default="MODASv", help="The output file prefix")

    ###################
    plot_ppb = subparsers.add_parser('plot_plainPhebox', help='Draw plain boxline plots for the phenotypes in different strains',
                                     usage='%(prog)s [options]')
    plot_ppb.add_argument('-i', dest="input", help="The expression file")
    plot_ppb.add_argument('-group', dest="group", help="The group file")
    plot_ppb.add_argument('-s_phes', dest="selected_phes", default=None,
                          help="The file listing selected genes or a string with comma separated")
    plot_ppb.add_argument('-s_strains', dest="selected_strains", default=None,
                          help="The file listing selected strains or a string with comma separated")
    plot_ppb.add_argument('-scale', dest="scale", default=False, action="store_true", help="Scale the phenotype values")
    plot_ppb.add_argument('-plot_type', dest="plot_type", choices=["single", "multi"], default="multi",
                          help="The type of graph to be plotted. Options are: single, multi. Default: multi")
    plot_ppb.add_argument("-o", dest="out_prefix", default="MODASv", help="The output file prefix")
    plot_ppb.add_argument('-width', dest="width", default=8, type=float,
                          help="The width of output plot")
    plot_ppb.add_argument('-height', dest="height", default=5, type=float,
                          help="The width of output plot")

    if len(sys.argv) <= 2:
        sys.argv.append('-h')
    args = parser.parse_args()
    main(args)
