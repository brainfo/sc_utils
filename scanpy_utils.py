## these are function tools that supplement scanpy ##
########### or useful when you use scanpy ###########
######### edited on 2022-02-09 by Ruby Jiang ########

import magic
import matplotlib.pyplot as plt
import pandas as pd
import scprep
import numpy as np
import scanpy as sc
import anndata
import re

## qc
def qc(data, name, mtprefix):
    """\
        Parameters
        -------
        data
            anndata object
        name
            name to identify your object
        mtprefix
            mitochondria gene prefix
    """
    data.var['mt'] = data.var_names.str.startswith(mtprefix)  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(data, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    sc.pl.violin(data, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True, groupby='status', save='figures/QC_{}.pdf'.format(name))
    fig, (ax1, ax2) = plt.subplots(1,2, constrained_layout=True)
    ax1_dict = sc.pl.scatter(data, x='total_counts', y='pct_counts_mt',ax=ax1, show=False, color='status')
    ax2_dict = sc.pl.scatter(data, x='total_counts', y='n_genes_by_counts', ax=ax2, show=False, color='status')
    plt.savefig('figures/{}_mt_tc.pdf'.format(name),bbox_inches='tight')

def filter_adata(adata, min_genes, percent_mt, percent, filter_mt=True):
    sc.pp.filter_cells(adata, min_genes=min_genes)
    adata_filter = adata[adata.obs.pct_counts_mt < percent_mt, :] 
    if percent>1:
        sc.pp.filter_genes(adata_filter, min_cells=percent)
    else:
        sc.pp.filter_genes(adata_filter, min_cells=percent*adata_filter.n_obs)
    if filter_mt:
        r = re.compile("^MT-")
        mt_adata = list(filter(r.match, adata_filter.var.index.values))
        adata_filter = adata_filter[:,~adata_filter.var.index.isin(mt_adata)]
    return adata_filter
    

def norm_hvg(adata, name):
    """\
        Parameters
        -------
        adata
            anndata object
    """ 
    ## Total-count normalize (library-size correct) the data matrix 𝐗 to 10,000 reads per cell
    sc.pp.normalize_total(adata, target_sum=1e4)
    ## Logarithmize the data:
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    print("Highly variable genes: %d"%sum(adata.var.highly_variable))
    sc.pl.highly_variable_genes(adata, save='{}_hvg.pdf'.format(name))

def cell_cycle_analysis(cell_cycle_genes,adata,name):
    """\
        Parameters
        -------
        cell_cycle_genes
            predefined cell cyle datasets in this util
        adata
            anndata object
        name
            sample name
    """   
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    scdata_cc_genes = adata[:, cell_cycle_genes]
    sc.tl.pca(scdata_cc_genes)
    sc.pl.pca_scatter(scdata_cc_genes, color='phase',save='{}_cell_cycle.pdf'.format(name))

def tsne_and_umap(adata, name):
    # all_data = GSE173193.concatenate(emab_filter, batch_key='trimester', batch_categories=['third','first'], join='inner', index_unique=None)
    ## here need a log file (json format)
    ## here need plot
    #sc.pl.highly_variable_genes(adata)
    sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')
    sc.pp.neighbors(adata, n_pcs =50)
    sc.tl.umap(adata)
    sc.tl.tsne(adata, n_pcs = 50)

    fig, axs = plt.subplots(1, 2, figsize=(8,4),constrained_layout=True)
    #sc.pl.tsne(corr_data, color="batchs", title="MNN tsne", ax=axs[0,0], show=False)
    sc.pl.tsne(adata, color="Cluster", title="tsne", ax=axs[0], show=False)
    #sc.pl.umap(corr_data, color="batchs", title="MNN umap", ax=axs[1,0], show=False)
    sc.pl.umap(adata, color="Cluster", title="umap", ax=axs[1], show=False)
    plt.savefig('figures/{}_projection.pdf'.format(name),bbox_inches='tight')