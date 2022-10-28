from .utils import *
import os
import sys

def exoTras_command():
    ''' Example of taking inputs for exoTras'''
    args = sys.argv[1:]
    if len(args) < 1:
        print("usage: exosomes_recognizer, source_tracker, ESAI_celltype, cellfree_simulator")


def cellfree_simulator(out_path, gene_exp_exo, gene_exp_cell, expect_UMI = [40, 70, 100, 130], exosomes_fraction = [0.005, 0.01, 0.05, 0.10], exosomes=500):
    all_num = [exosomes/f for f in exosomes_fraction]
    for i in expect_UMI:
        for j in all_num:
            exo_count = np.random.poisson(lam=i, size=exosomes)
            all_count = np.random.poisson(lam=i, size=j)

            exo_data = process_random(gene_exp_exo, exo_count)
            exo_data.columns = ['exo'+str(i) for i in exo_data.columns.values]
            cell_data = process_random(gene_exp_cell, all_count)
            exo_cell_data = pd.concat([cell_data, exo_data], axis=1).fillna(0)

            new_index = [ i.split(' /// ')[-1].split('-')[0] for i in exo_cell_data.index.values]
            exo_cell_data.index = new_index
            exo_cell_data.to_csv(out_path + '/simulation_UMI' + str(i) + '_fraction' + str(int(exosomes/j)) + '.csv')
    
    return('Simulation done!')


def exosomes_recognizer(sample_file, out_path, input_path=None, species='Homo', predefine_threads=-2, get_only=False, score_t = None):
    if predefine_threads == -2:
        threads = cpu_count()-2
    else:
        threads = predefine_threads
    exo_list, score_t = get_exo_list(species, score_t)
    exo_list_s = list(set(exo_list))

    sample_log = get_sample(sample_file)
    len_sample = len(sample_log)

    for sample in sample_log:
        ##read adata
        if input_path!= None:
            each_adata = str(input_path) + '/' + str(sample)
        else:
            each_adata = sample

        adata = read_adata(each_adata, get_only=False)
        adata = filter_adata(adata)

        iteration_list = exo_list_s
        inter_gene = []
        pth = 0.005
        flag = 1
        inter_out_p = pd.DataFrame(list(range(0, adata.obs.shape[0])))
        inter_adata = copy.copy(adata[(adata.obs['total_counts'] < 200) & (adata.obs['total_counts'] > 20),:])
        inter_adata = inter_adata.copy()
        sc.pp.normalize_total(inter_adata, target_sum=1e2)
        # inter_adata.X = (inter_adata.X / inter_adata.X.sum(1).reshape(inter_adata.shape[0], 1))*100
        ##iterations 
        for itera in range(10):
            if(len(iteration_list) == 0):
                print('no genes')
                flag = 0
                break

            thershold = 0
            iteration_list_old = iteration_list
            inter_adata, iteration_list = iteration(inter_adata, iteration_list, thershold, threads=20, number_g = 30, alpha = 0.10)
            print(itera, len(iteration_list))
                
            if(len((iteration_list_old) & (iteration_list)) == max(len(iteration_list_old), len(iteration_list))):
                break
            else:
                inter_gene.append(iteration_list_old)
                inter_out_p = pd.concat((inter_out_p, inter_adata.obs['score']), axis=1, ignore_index=True)

        if flag == 1:
            genes = iteration_list[0:15]
            same = [i for i in (genes) if i in adata.var_names]
            thershold = 0
            adata = final_menrich(adata, same, thershold, threads)

        elif (len(iteration_list_old) > 0) & (len(iteration_list_old) < len(exo_list_s)):
            genes = iteration_list_old[0:15]
            same = [i for i in (genes) if i in adata.var_names]
            thershold = 0
            adata = final_menrich(adata, same, thershold, threads)

        else:
            print('Exosomes undetectable  in ' + str(sample))
            continue
        
        ## write files
        if len_sample > 1:
            adata.write(str(out_path) + '/tmp_out/' + sample + '/raw_' + sample + '.h5ad')
            with open(str(out_path) +'./tmp_out/' + sample + '/itera_gene.txt', 'wb') as fp:
                pickle.dump(inter_gene, fp)
        else:
            adata.obs['exoTras'] = ['Exosomes' if i else 'False' for i in (adata.obs['total_counts'] < 200) & (adata.obs['padjust'] < 1e-15)]
            adata.write(str(out_path) + '/raw_' + sample + '.h5ad')
            with open(str(out_path)+"/itera_gene.txt", "w") as f:
                for s in genes:
                    f.write(str(s) +"\n")

    ## Aggregation
    if len_sample > 1:
        adata_com = read_project(out_path, sample_log)
        genes_out = get_genes(out_path, sample_log)
        genes = count_genes(genes_out)
        iteration_list = set(genes['0'])#exo_list_s
        inter_gene = []
        inter_out_p = pd.DataFrame(list(range(0, adata_com.obs.shape[0])))
        pth = 0.005
        same = [i for i in (iteration_list) if i in adata_com.var_names]
        thershold = 0
        adata_com = final_menrich(adata_com, same, thershold, threads)

        adata_com.obs['exoTras'] = ['Exosomes' if i else 'False' for i in (adata_com.obs['total_counts'] < 200) & (adata_com.obs['padjust'] < score_t)]

        ## write files in this project
        adata_com.write(str(out_path) + './raw_' + 'exoTras' + '.h5ad')
        adata_exo = copy.copy(adata_com[(adata_com.obs['total_counts'] < 200) & (adata_com.obs['padjust'] < score_t), :])
        adata_exo = adata_exo.copy()
        adata_exo.write(str(out_path) + './exosomes_exoTras.h5ad')

