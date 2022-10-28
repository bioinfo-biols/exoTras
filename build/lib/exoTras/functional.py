## source 
import gseapy as gp
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
from scipy.spatial import cKDTree
import copy
import sys
import os
import pickle

def source_biogenesis(adata_cell, OBScelltype='celltype', Xraw = True, normalW=True):
    if Xraw:
        X_input = adata_cell.raw
    else:
        X_input = adata_cell
        
    if normalW:
        X_norm = sc.pp.scale(X_input.X, zero_center=True, max_value=None, copy=False)
    else:
        X_norm = X_input.X
    
    gsea_pval = []
    num_clusters = len(adata_cell.obs[OBScelltype].cat.categories)

    gmt_path = os.path.dirname(os.path.split(os.path.realpath(__file__))[0]) + '/exosomes.gmt'

    for i in range(num_clusters):
        i = adata_cell.obs[OBScelltype].cat.categories[i]
        gene_rank = pd.DataFrame({'exp': np.array(X_norm[adata_cell.obs[OBScelltype] == str(i), :].mean(axis=0))}, index = X_input.var_names)

        res = gp.prerank(rnk=gene_rank, gene_sets=gmt_path)
        terms = res.res2d.index
        gsea_pval.append([i, res.results[terms[0]]['nes'], res.results[terms[0]]['pval']])        

    gsea_pval_dat = pd.DataFrame(gsea_pval, index=adata_cell.obs[OBScelltype].cat.categories, columns = ['num', 'enrich', 'p'])
    gsea_pval_dat['log1p'] = -np.log10(gsea_pval_dat['p']+1e-4) * np.sign(gsea_pval_dat['enrich'])

    return(gsea_pval_dat)


def near_neighbor(adata_combined, OBSsample='batch', OBSexo='exo', OBScelltype='celltype', OBSMpca='X_pca', cellN=10):
    ## run sc.pp.pca(adata_combined) before
    near_neighbor = []
    for sample in adata_combined.obs[OBSsample].unique().astype(str):
        tse_ref = copy.copy(adata_combined[(adata_combined.obs[OBSexo] == '0') & (adata_combined.obs[OBSsample] == sample),])
        cell_tree = cKDTree(tse_ref.obsm[OBSMpca], leafsize=100)
        
        tmp_umap = adata_combined[(adata_combined.obs[OBSexo] == '1') & (adata_combined.obs[OBSsample] == sample),].obsm[OBSMpca]
        for i in range(tmp_umap.shape[0]):
            TheResult = cell_tree.query(tmp_umap[i,], k=10)

            near_neighbor.append([sample, i] + list(tse_ref.obs[OBScelltype].iloc[TheResult[1]]))#tmp_umap.obs['clusters'][i]] +

    near_neighbor_dat = pd.DataFrame(near_neighbor)
    return(near_neighbor_dat)

def preprocess_source(adata_exo, adata_cell, OBScelltype='celltype', OBSexo='exo'):
    ## cell type
    adata_cell_raw = copy.copy(adata_cell.raw.to_adata())

    adata_exo.obs[OBScelltype] = 'Exosomes'
    adata_exo.obs[OBScelltype] = pd.Series(adata_exo.obs[OBScelltype], dtype="category")

    adata_combined = adata_cell_raw.concatenate(adata_exo, batch_key = OBSexo)

    adata_combined.obs[OBScelltype] = pd.Categorical(adata_combined.obs[OBScelltype], \
        categories = np.append(adata_cell_raw.obs[OBScelltype].cat.categories.values, 'Exosomes'), ordered = False)
    
    adata_combined.raw = adata_combined
    sc.pp.normalize_total(adata_combined, target_sum=1e4)
    sc.pp.log1p(adata_combined)
    sc.pp.highly_variable_genes(adata_combined, min_mean=0.0125, max_mean=3, min_disp=0.5)
    # sc.pl.highly_variable_genes(Normal_combined)
    adata_combined = adata_combined[:, adata_combined.var.highly_variable]#highly_variable

    sc.pp.pca(adata_combined)
    sc.pp.neighbors(adata_combined)
    sc.tl.umap(adata_combined)

    return(adata_combined)

def source_tracker(adata_exo, adata_cell, OBSsample='batch', OBScelltype='celltype', OBSexo='exo', OBSMpca='X_pca', cellN=10, Xraw = True, normalW=True):

    adata_combined = preprocess_source(adata_exo, adata_cell, OBScelltype='celltype', OBSexo='exo')
    gsea_pval_dat = source_biogenesis(adata_cell, OBScelltype='celltype', Xraw = True, normalW=True)
    near_neighbor_dat = near_neighbor(adata_combined, OBSsample='batch', OBSexo='exo', OBScelltype='celltype', OBSMpca='X_pca', cellN=10)
    
    near_neighbor_dat['times'] = ''
    near_neighbor_dat['type'] = ''
    for i in range(near_neighbor_dat.shape[0]):
        ## iteration for all exosomes
        tmp_times = near_neighbor_dat.iloc[i,2:12].value_counts(sort = True)
        near_neighbor_dat.loc[i, 'times'] = tmp_times[0]
        tmp_keys = near_neighbor_dat.iloc[i,2:12].value_counts(sort = True).keys()
        tmp_times_fil = (tmp_times + gsea_pval_dat.loc[tmp_keys, 'log1p']*4).sort_values(ascending = False)
        near_neighbor_dat.loc[i, 'type'] = tmp_times_fil.keys()[0]

    near_neighbor_dat.index = adata_exo.obs.index
    celltype_e_number = pd.DataFrame(near_neighbor_dat.type.value_counts())

    adata_exo.obsm['source'] = near_neighbor_dat[[0, 1, 'type']]
    adata_exo.obsm['source'].columns = ['sample', 'i', 'type']

    return([celltype_e_number, adata_exo, adata_combined])


def ESAI_celltype(adata_exo, adata_cell, OBSsample='batch', OBScelltype='celltype'):
    ## calculate ESAI in cell type
    exo_type_count_total = pd.DataFrame(adata_exo.obsm['source'].groupby(['sample','type'])['sample'].count())
    cell_type_count = pd.DataFrame(adata_cell.obs.groupby([OBSsample, OBScelltype])[OBSsample].count())

    exo_activity = []
    for i, j in cell_type_count.index:
        if (i,j) in exo_type_count_total.index:
            a = exo_type_count_total.loc[(i, j), 'sample']
            b = cell_type_count.loc[(i, j), OBSsample]
            if b > 0:
                exo_activity.append([i, j, a, b, a / b])
            else:
                print([i, j, a, b, a / b])
        else:
            a = 0
            b = cell_type_count.loc[(i, j), OBSsample]
            exo_activity.append([i, j, a, b, 0])   
              
    exo_activity_dat = pd.DataFrame(exo_activity)
    exo_activity_dat_pivot = pd.pivot_table(exo_activity_dat, index=[0], columns= 1, values=4)

    return(exo_activity_dat_pivot)

