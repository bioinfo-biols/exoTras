import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import copy
import sys
import os
import pickle
import statsmodels.api as sm
from multiprocessing import Pool, cpu_count
# import re
# import random
# random.seed( 60 )

def get_sample(sample_log):
    names_list = []
    with open(sample_log, 'r') as  f:
        for line in f.readlines():
            names_list.append(line.strip())
    
    return(names_list)


def read_adata(adata_path, get_only=False):
    if (sys.argv[1]).endswith('.h5'):
        adata = sc.read_10x_h5(adata_path)
    elif (sys.argv[1]).endswith('.h5ad'):
        adata = sc.read_h5ad(adata_path)
    else:
        adata = sc.read_10x_mtx(adata_path + '/outs/raw_feature_bc_matrix/')

    return(adata)

def read_project(output_path, names_list):

    names_list_1 = copy.copy(names_list)
    # adata_list = []
    for i in (names_list):
        if os.path.exists(output_path + '/tmp_out/' + i):
            exec('g{} = sc.read_h5ad("{}/{}/{}/raw_{}.h5ad")'.format(i, output_path+ 'tmp_out', i, i))
            exec('g{}.var_names_make_unique()'.format(i))
            # exec('adata_list.append(g{})'.format(i))
        else:
            print(i)
            names_list_1.remove(i)

    exec('adata_com = g{}.concatenate(g{}, batch_categories = names_list_1)'.format(names_list_1[0], ", g".join(names_list_1[1:])))

    adata_com.write(output_path+'/all.h5ad')

    return(adata_com)

def get_genes(output_path, names_list):
    genes_out = []
    for i in (names_list):
        path_gene = output_path + '/tmp_out/' + i + '/itera_gene.txt'
        if os.path.exists(path_gene):
            with open (path_gene, 'rb') as fp:
                inter_gene = pickle.load(fp)
                genes_out.append(inter_gene[len(inter_gene)-1])
    
    return(genes_out)

def count_genes(genes_out):
    cal_gene_out = []
    for i in genes_out:
        cal_gene_out.extend(list(i))

    pd.Series(cal_gene_out).value_counts().head(15).index.to_series().to_csv('./top15.csv')
    genes = pd.read_csv('./top15.csv')

    return(genes)


def filter_adata(adata):
    adata.var['mt'] = adata.var_names.str.startswith(('mt-', 'MT-'))
    adata.var['ribo'] = adata.var_names.str.startswith(('Rps', 'Rpl', 'RPS', 'RPL'))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], inplace=True, percent_top='', log1p =False)
    adata = adata[(adata.obs.pct_counts_mt < 15) & (adata.obs.pct_counts_ribo < 30), :]
    sc.pp.filter_cells(adata, min_genes=6)
    sc.pp.filter_genes(adata, min_cells=5)

    return(adata)

def get_exo_list(species, score_t='1e-15'):
    ## get initial gene set
    ##exo gene symbols
    if species == 'Homo':
        exo_list = [x.strip() for x in open(os.path.dirname(os.path.split(os.path.realpath(__file__))[0]) + '/id_homo_rna_exo')]
    elif species == 'Mus':
        exo_list = [x.strip() for x in open(os.path.dirname(os.path.split(os.path.realpath(__file__))[0]) + '/id_mus_rna_exo')]
        score_t = '1e-10'
    
    exo_list_count = pd.value_counts(exo_list, sort=False)
    exo_list = exo_list_count[exo_list_count > 0].index
    
    return([exo_list, score_t])


import env
def process_inter(i):
    k = (env.Inter_adata[i, env.Same].X > env.Thershold).sum()
    M = max(env.Inter_adata.obs['n_genes'])
    n = len(env.Same)
    N = int(env.Inter_adata[i, :].obs['n_genes'])
    return([scipy.stats.hypergeom.sf(k, M, n, N), N, k/N, n/M])

def corr_genes(i):
    #i, inter_adata
    if (len(np.unique(env.Inter_adata.X.A[:, i])) > 1):
        tmp = scipy.stats.spearmanr(env.Inter_adata.X.A[:, i], -np.log(env.Inter_adata.obs['score'] + 10**-50))
        return(tmp.correlation)
    else:
        return(-1)

def multi_enrich(inter_adata, same, iteration_list, thershold, threads):
    len((iteration_list)), len(same), max(inter_adata.obs['n_genes'])

    out0 = []
    item_list = range(inter_adata.n_obs)
    pool = Pool(threads, env.initializer, (inter_adata, same, thershold))
    out0 = pool.map(process_inter, item_list)
    pool.close()
    pool.join()

    ##
    out = pd.DataFrame(out0)
    out['score'] = sm.stats.fdrcorrection(out[0])[1]
    # iteration_list_old = iteration_list
    
    inter_adata.obs['score'] = out['score'].tolist()

    return(inter_adata)


def multi_cor(inter_adata, threads):
    result_tmp = []

    pool = Pool(threads, env.initializer_simple, (inter_adata))
    item_list = range(inter_adata.X.shape[1])
    
    result_tmp = pool.map(corr_genes, item_list)
    pool.close()
    pool.join()

    return(result_tmp)

def representative_gene(result_tmp, number_g = 30, alpha = 0.10):
    max_number = []
    max_index = []
    t = copy.deepcopy(result_tmp)
    for _ in range(number_g):
        number = max(t)
        index = t.index(number)
        t[index] = -1
        if (number > alpha):
            max_index.append(index)
            max_number.append(number)

    return([max_index, max_number])

def get_iteration(inter_adata, max_index):
    # p_thre_filter = min(max_number)
    #tt_list  = inter_adata.var.index[[(i >= 0) & (i >= p_thre_filter) for i in result_tmp]]
    tt_list = inter_adata.var.index[max_index]
    iteration_list = list(tt_list)# & exo_list_s  set()

    return(iteration_list)

def iteration(inter_adata, iteration_list, thershold, threads=20, number_g = 30, alpha = 0.10):
    same = [i for i in (iteration_list) if i in inter_adata.var_names]
    inter_adata = multi_enrich(inter_adata, same, iteration_list, thershold, threads)
    result_tmp = multi_cor(inter_adata, threads)

    
    max_index, max_number = representative_gene(result_tmp, number_g = 30, alpha = 0.10)
    iteration_list = get_iteration(inter_adata, max_index)

    return(inter_adata, iteration_list)
    
def process(i):
    k = (env.Adata[i, env.Same].X > env.Thershold).sum()
    M = max(env.Adata.obs['n_genes'])
    n = len(env.Same)
    N = int(env.Adata[i, :].obs['n_genes'])
    return([scipy.stats.hypergeom.sf(k, M, n, N), N, k/N, n/M])


def final_menrich(adata, same, thershold, threads):

    out0 = []
    item_list = range(adata.n_obs)
    pool = Pool(threads, env.initializer_adata, (adata, same, thershold))
    out0 = pool.map(process, item_list)
    pool.close()
    pool.join()

    ##
    out = pd.DataFrame(out0)
    out['score'] = sm.stats.fdrcorrection(out[0])[1]
    #iteration_list_old = iteration_list

    adata.obs['score'] = out['score'].tolist()
    
    return(adata)



## cell free droplets simulation
def zinb_genes(i, gene_exp_a, num):
    cells = num #number of cells
    n = 1
    p = 0.8# 1 - droplet rate
    
    exp = gene_exp_a.iloc[i]#exp_i# gene expression mean
    success = 1
    
    out_zinb = np.random.negative_binomial(success, success/(success+exp), size=cells) * \
        np.random.binomial(n, p, size=cells)
    return(pd.Series(data = out_zinb, name = gene_exp_a.index[i]))
    
def run_weigths(gene_exp_a, num=1000):
    
    weight = []
    def log_result(result):
        # This is called whenever foo_pool(i) returns a result.
        # result_list is modified only by the main process, not the pool workers.
        weight.append(result)

    pool = Pool(cpu_count()-4)

    for i in range(len(gene_exp_a)):
        pool.apply_async(zinb_genes, args = (i, gene_exp_a, num), callback = log_result)
        
    pool.close()
    pool.join()
    
#     weight = (gene_exp_a.apply(zinb_genes, axis = 1))
    
    weight = pd.DataFrame(weight)#.T.fillna(0)
    
    return(weight)


def process_random(gene_exp, count):
    num = len(count)
    out = run_weigths(gene_exp, num)
    
    return(round((out / out.sum(0)) * count))






