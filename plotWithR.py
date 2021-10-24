#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: plotWithR.py
Author: CrazyHsu @ crazyhsu9527@gmail.com
Last modified: 2021-09-25 23:47:20
'''

plotQQplotR = '''
    library(data.table)
    library(CMplot)
    plotQQplot <- function(gwas_res, out_prefix, sep) {
        sink('/dev/null')
        data <- fread(gwas_res, sep=sep, data.table=getOption("datatable.fread.datatable", FALSE))
        CMplot(data, plot_type='q', col='grey30', conf_int_col='gray', signal_col='red', multracks=False, LOG10=True,
               file='pdf', dpi=300, file_prefix=out_prefix)
        sink()
    }
'''

plotManhattanPlotR = '''
    library(data.table)
    library(CMplot)
    plotManhattan <- function(gwas_res, thresholdi, out_prefix, sep) {
        sink('/dev/null')
        data <- fread(gwas_res, sep=sep, data.table=getOption("datatable.fread.datatable", FALSE))
        thresholdi <- unlist(strsplit(thresholdi,","))
        if(thresholdi[0] == "None"){
            thresholdi[0] <- 1.0/dim(data)[1]
        }
        thresholdi <- as.numeric(thresholdi)
        lim <- -min(d[, 3:dim(d)[2]]) + 2
        CMplot(data, plot_type='m', col=c("grey30", "grey60"), ylim=c(2, lim),
               threshold=thresholdi,
               cex=c(0.5, 0.5, 0.5), signal_cex=c(0.5, 0.5, 0.5),
               threshold_col=c('red', 'green', 'blue'), chr_den_col=NULL,
               amplify=True, signal_pch=c(19, 19, 19), dpi=300,
               signal_col=c('red', 'green', 'blue'), multracks=False, LOG10=True, file='pdf', file_prefix=out_prefix)
        base.sink()
    }
'''

plotForestPlotR = '''
    suppressMessages(library(ggplot2))
    suppressMessages(library(scales))
    suppressMessages(library(data.table))
    options(warn=-1)
    plotForest <- function(mr_res, order, out_prefix, sep){
        data <- fread(mr_res, sep=sep, data.table=getOption("datatable.fread.datatable", FALSE))
        colnames(data) <- gsub("/", ".", colnames(data))\
        colnames(data) <- gsub("^(\\d+)", "X\\1", colnames(data))
        for(mTrait in unique(data[,"mTrait"])){
            sub_data <- data[data$mTrait==mTrait,]
            for(snp in unique(sub_data[,"snp"])){
                d <- sub_data[sub_data$snp==snp,]
                if(order=="True"){
                    d <- d[order(d$pTrait),]
                }
                d$pTrait <- factor(d$pTrait, levels = d$pTrait)
                d$std <- sqrt((d$effect) ^ 2 / d[, 5])
                d$min <- d$effect - d$std
                d$max <- d$effect + d$std
                # effect
                ggplot(data=d, aes(x=effect, y=pTrait, color=pTrait)) +
                    geom_vline(xintercept=0, linetype="dashed", color="gray", size=1) +
                    geom_point(size=1.5, shape=15) +
                    geom_errorbarh(aes(xmin=min, xmax=max), height=.1) +
                    theme_bw() +
                    theme(legend.title = element_blank(),
                          legend.position = 'none',
                          plot.title = element_text(hjust=0.5),
                          panel.grid = element_blank(),
                          axis.line = element_line(colour='black', size=0.1),
                          axis.text.x = element_text(color='black'),
                          axis.text.y = element_text(color=hue_pal()(nrow(d))))+
                          xlab('') + ylab('') +
                          ggtitle(paste(mTrait, snp, sep='_'))
                ggsave(paste(out_prefix, mTrait, snp, 'forestplot.pdf', sep='_'), device='pdf', height=6, width=4)
            }
        }
    }
'''

plotScatterPlotR = '''
    suppressMessages(library(ggplot2))
    suppressMessages(library(ggrepel))
    suppressMessages(library(data.table))
    options(warn=-1)
    plotScatterPlot <- function(mr_res, group, out_prefix, sep){
        data <- fread(mr_res, sep=sep, data.table=getOption("datatable.fread.datatable", FALSE))
        colnames(data) <- gsub("/", ".", colnames(data))\
        colnames(data) <- gsub("^(\\d+)", "X\\1", colnames(data))
        for(mTrait in unique(data[,"mTrait"])){
            sub_data <- data[data$mTrait==mTrait,]
            for(snp in unique(sub_data[,"snp"])){
                d <- sub_data[sub_data$snp==snp,]
                if(group=="None"){
                    res <- d
                    res$group <- 'group'
                }else{
                    groupD <- read_csv(group)
                    if(nrow(groupD)==0){
                        res <- d
                        res$group <- 'group'
                    }else{
                        res <- merge(mr_res, d, by.x ='id', by.y='pTrait')
                        res$group <- factor(res$group)
                        res <- res[order(res$group),]
                    }
                }
                res$log10p <- -log10(res$pvalue)
                res[res$effect<0,'log10p'] <- -res[res$effect<0,'log10p']
                res$pos <- seq(1,nrow(res))
                scatter_plot <- ggplot(data=res,aes(x=pos,y=log10p,color=group))+
                      geom_hline(yintercept=log10(0.05/nrow(res)), linetype="dashed", color = "red", size=1)+
                      geom_hline(yintercept=-log10(0.05/nrow(res)), linetype="dashed", color = "red", size=1)+
                      geom_point(size=1.5,shape=19)+
                      theme_bw()+
                      theme(legend.title = element_blank(),
                            plot.title = element_text(hjust = 0.5),
                            panel.grid = element_blank(),
                            axis.line = element_line(colour = 'black',size=0.1),
                            axis.text = element_text(color = 'black'),
                            axis.text.x = element_blank(),
                            axis.ticks.length.x = unit(0,'cm'))+
                      xlab('')+ylab('')
                if(nrow(group)!=0){
                    scatter_plot <- scatter_plot + geom_label_repel(aes(label=id))
                }
                ggsave(paste(out_prefix, mTrait, snp, 'scatterplot.pdf', sep='_'), plot=scatter_plot, device = 'pdf',height = 5,width = 9)
            }
        }
    }
'''

plotScatterPsR = '''
    suppressMessages(library(ggplot2))
    plotScatterPs <- function(ps_file, group_file, out_prefix, sep) {
        data <- read.csv(ps_file)
        data <- data[c(1,2,3)]
        group <- read.csv(group_file)
        if(nrow(group)!=0){
                data$group <- group[match(data[, 1], group[, 1]), 2]
        }else{
            data$group <- 'one'
        }
        names(data) <- c('sample','PC1','PC2','group')
        scatter_plot <- ggplot(data = data,aes(x=PC1,y=PC2,color=group))+
            geom_point(shape=19, size=1.5)+
            theme_bw() +
            theme(axis.text.x = element_text(colour = "black", size = 12, vjust =1),
                  axis.text.y = element_text(colour = "black", size = 12, hjust =1),
                  axis.title = element_text(margin=margin(10, 5, 0, 0), color = "black", size = 14),
                  axis.title.y = element_text(angle=90),
                  panel.border = element_rect(color = 'black', fill=NA, size = 1)
                 )+
            xlab('')+ylab('')
        if(nrow(group)==0){
            scatter_plot <- scatter_plot + theme(legend.position='none')
        }
        ggsave(paste(out_prefix, 'scatter_ps.pdf', sep='_'), plot=scatter_plot, width = 4.2, height = 3)
    }
'''

plotHeatmapR = '''
    suppressMessages(library(pheatmap))
    plotHeatmap <- function(input, out_prefix, cluster_rows, cluster_cols, order_col_by, order_row_by, anno_row, 
    anno_col, scale, show_rownames, show_colnames) {
        data <- read.csv(input, row.names=1)
        
        if(order_col_by!="None"){
            order_col_by <- read.csv(order_col_by)
            order_col_by <- as.vector(order_col_by[,1])
            data <- data[, order_col_by]
        }
        if(order_row_by!="None"){
            order_row_by <- read.csv(order_row_by)
            order_row_by <- as.vector(order_row_by[,1])
            data <- data[order_row_by, ]
        }
        
        if(anno_col!="None"){
            annotation_col <- read.csv(anno_col)
            colnames(annotation_col) <- annotation_col[,1]
            annotation_col <- annotation_col[,2:dim(annotation_col)[2]]
        }else{
            anno_col <- NA
        }
        if(anno_row!="None"){
            annotation_row <- read.csv(anno_row)
            colnames(annotation_row) <- annotation_row[,1]
            annotation_row <- annotation_row[,2:dim(annotation_row)[2]]
        }else{
            anno_row <- NA
        }
        
        out_file <- paste0(out_prefix, ".pdf")
        pheatmap(data, filename=out_file, scale=scale, cluster_rows=cluster_rows, cluster_cols=cluster_cols,
                 show_rownames=show_rownames, show_colnames=show_colnames, annotation_row=anno_row, 
                 annotation_col=anno_col)
    }
'''

plotBoxR = '''
    suppressMessages(library(ggplot2))
    suppressMessages(library(data.table))
    options(warn=-1)
    plotBox <- function(input_mat_file, snp_file, group_file, selected_genes, selected_snps, selected_strains, scale, plot_type, out_prefix){
        input_mat <- fread(input_mat_file, data.table=getOption("datatable.fread.datatable", FALSE))
        snp_mat <- fread(snp_file, data.table=getOption("datatable.fread.datatable", FALSE))
        group_mat <- read.csv(group_file)
        if(selected_genes!="None"){
            selected_genes <- unlist(strsplit(selected_genes))
        }else{
            selected_genes <- c(colnames(input_mat)[2])
        }
        if(){
            selected_snps <- unlist(strsplit(selected_snps))
        }else{
            selected_snps <- c(colnames(snp_mat)[2])
        }
        if(selected_strains!="None"){
            selected_strains <- unlist(strsplit(selected_strains))
        }else{
            selected_strains <- rownames(input_mat)
        }
         
        sub_input_mat <- input_mat[select_strains, selected_genes]
        sub_snp_mat <- snp_mat[select_strains, select_snps]
        sub_group_mat <- group_mat[select_strains, ]
        conames(sub_group_mat) <- c("RIL", "Subpop")
        
        for(snp in selected_snps){
            genotype <- sub_snp_mat[, snp]
            colnames(genotype) <- c("RIL", "Genotype")
            if(plot_type=="single_gene"){
                for(gene in selected_genes){
                    single_gene_mat <- sub_input_mat[, gene]
                    merged_exp <- merge(x=genotype, y=single_gene_mat, by="RIL")
                    merged_exp <- na.omit(merged_exp)
                    colnames(merged_exp) <- c("RIL", "Genotype", "Expression")
                    merged_exp$Expression <- (merged_exp$Expression-mean(merged_exp$Expression))/sd(merged_exp$Expression)
                    ggplot(data=merged_exp, aes(RIL, Expression, fill=Genotype))+
                        stat_boxplot(geom='errorbar', width=0.2)+
                        geom_boxplot()+
                        theme_bw() +
                        theme(axis.text.x = element_text(colour = "black", size = 12, vjust =1 ),
                            axis.text.y = element_text(colour = "black", size = 12, hjust =1 ),
                            axis.title = element_text(margin=margin(10, 5, 0, 0), color = "black", size = 14),
                            axis.title.y = element_text(angle=90),
                            panel.border = element_rect(color = 'black', fill=NA, size = 1),
                            legend.position = 'none'
                        ) + xlab('')
                    ggsave(paste(out_prefix, gene, 'boxplot.pdf', sep='_'), width = width, height = 3)
                }
            }else if(plot_type=="multi_gene") {
                gene_n <- length(selected_genes)
                if(gene_n==1){
                    width <- 2
                }else if(gene_n==2){
                    width <- 2.5
                }else{
                    width <- 2.5 + 0.5 * (gene_n-2)
                }
                if(width >= 8){
                    width <- 8
                }
            
                merged_exp <- merge(x=genotype, y=sub_input_mat, by="RIL")
                merged_exp <- na.omit(merged_exp)
                tmp_exp <- merged_exp[,2:dim(merged_exp)[2]]
                merged_exp[,2:dim(merged_exp)[2]] <- sapply(tmp_exp, function(tmp_exp) (tmp_exp-mean(tmp_exp))/sd(tmp_exp))
                
                melt_merged_exp <- melt(merged_exp, id.vars = c("RIL", "Genotype"), variable.name = "gene", value.name = "expression")
                colnames(merged_exp) <- c("RIL", "Genotype", "Expression")
                ggplot(data=merged_exp, aes(RIL, Expression, fill=Genotype))+
                    stat_boxplot(geom='errorbar', width=0.2)+
                    geom_boxplot()+
                    theme_bw() +
                    theme(axis.text.x = element_text(colour = "black", size = 12, vjust =1 ),
                        axis.text.y = element_text(colour = "black", size = 12, hjust =1 ),
                        axis.title = element_text(margin=margin(10, 5, 0, 0), color = "black", size = 14),
                        axis.title.y = element_text(angle=90),
                        panel.border = element_rect(color = 'black', fill=NA, size = 1),
                        legend.position = 'none'
                    ) + xlab('')
                ggsave(paste(out_prefix, 'multiGenes', 'boxplot.pdf', sep='_'), width = width, height = 3)
            }else if(plot_type=="single_pop") {
                for(gene in selected_genes){
                    single_gene_mat <- sub_input_mat[, gene]
                    colnames(single_gene_mat) <- c("RIL", "Expression")
                    df_lst <- list(sub_group_mat, genotype, single_gene_mat)
                    merged_exp <- Reduce(function(x,y) merge(x,y,all="RIL“), df_lst)
                    merged_exp <- na.omit(merged_exp)
                    merged_exp$Expression <- (merged_exp$Expression-mean(merged_exp$Expression))/sd(merged_exp$Expression)
                    
                    ggplot(data=merged_exp, aes(Subpop, Expression, fill=Genotype))+
                        stat_boxplot(geom='errorbar', width=0.2)+
                        geom_boxplot()+
                        theme_bw() +
                        theme(axis.text.x = element_text(colour = "black", size = 12, vjust =1 ),
                            axis.text.y = element_text(colour = "black", size = 12, hjust =1 ),
                            axis.title = element_text(margin=margin(10, 5, 0, 0), color = "black", size = 14),
                            axis.title.y = element_text(angle=90),
                            panel.border = element_rect(color = 'black', fill=NA, size = 1),
                            legend.position = 'none'
                        ) + xlab('')
                    ggsave(paste(out_prefix, gene, snp, 'pop.boxplot.pdf', sep='_'), width = width, height = 3)
                }
            }
        }
    }
'''