import argparse
import makeData as md
import simulator as sl

# 使用argparse 通过命令行设置参数
def parse_args():
    parse = argparse.ArgumentParser(description='Simulation of combining cell data')
    parse.add_argument('-nc', '--num_cells', help='Enter the number of cells to be generated')
    parse.add_argument('-ng', '--num_genes', help='Enter the number of genes to be generated')
    parse.add_argument('-ct', '--cell_types', help='Enter the number of cell types to be generated')
    parse.add_argument('-nmg', '--num_degs', help='Enter the number of DEGs to be generated')
    parse.add_argument('-gec', '--gene_expression_count', help='The average expression of non-zero expression genes')
    parse.add_argument('-dp', '--scSimulate_data_h5ad_path', default='simulate_data.h5ad', help='The data saved in '
                                                                                                'the simulation is in '
                                                                                                'h5ad formate')
    parse.add_argument('-cs', '--chip_size_mm', type=int, help='size of chip, unit is mm')
    parse.add_argument('-minr', '--min_radius_mm', type=int, help='minimum radius of cells, unit is mm')
    parse.add_argument('-maxr', '--max_radius_mm', type=int, help='maximum radius of cells, unit is mm')
    args = parse.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    num_cells = args.num_cells
    num_genes = args.num_genes
    cell_types = args.cell_types
    num_degs = args.num_degs
    gene_expression_count = args.gene_expression_count
    scSimulate_data_h5ad_path = args.scSimulate_data_h5ad_path
    chip_size_mm = args.chip_size_mm
    min_radius_mm = args.min_radius_mm
    max_radius_mm = args.max_radius_mm

    # simulate data
    madata = md.combine_cell_gene_with_markers(num_cells, num_genes, cell_types, num_degs, gene_expression_count)

    # generate h5ad file
    md.h5ad_file(madata, scSimulate_data_h5ad_path)

    # simulate
    sl.simulate(madata, cell_types, chip_size_mm, min_radius_mm, max_radius_mm, num_cells)

# -nc 100 -ng 1000 -ct 5 -nmg 10 -gec 125 -dp 'simulate_data.h5ad' -cs 1 -minr 0.005 -maxr 0.01