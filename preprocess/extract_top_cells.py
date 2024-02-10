import os
import argparse
import scanpy as sc

argparser = argparse.ArgumentParser()
argparser.add_argument(
    'count_matrix_dir', help='Path to the count matrix directory')
argparser.add_argument('output_dir', help='Path to the output directory')

args = argparser.parse_args()
mtx_dir = args.count_matrix_dir
output_dir = args.output_dir
adata = sc.read_10x_mtx(mtx_dir, cache=True)

forseti_sample_dir = os.getcwd()
output_dir = os.path.join(
    forseti_sample_dir, 'STARsolo_out')
adata = sc.read_10x_mtx(os.path.join(
    forseti_sample_dir, 'STARsolo_out', 'filtered_feature_bc_matrix'), cache=True)

sc.pp.filter_cells(adata, min_genes=3000)
adata.n_obs
with open(os.path.join(output_dir, 'top200cells.txt'), 'w') as f:
    for cellbarcode in adata.obs_names.tolist():
        cellbarcode = cellbarcode.split('-')[0]
        f.write('CB:Z:'+cellbarcode + '\n')
# original cell number: 1222
# min_gene =1000: from 1222 cells -> 1101 cells
# min_gene = 2000: from 1222 cells -> 540 cells
# min_gene = 3000: from 1222 cells -> 214 cells
