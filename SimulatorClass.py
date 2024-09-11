import makeData as md
import simulator as sl
import bin_operation as bo


class Simulator:

    def __init__(self):
        print("调用构造方法")

    def make_data(self,
                  num_cells,
                  num_genes,
                  num_types,
                  num_marker_genes,
                  gene_expression_count):
        adata = md.combine_cell_gene_with_markers(num_cells,
                  num_genes,
                  num_types,
                  num_marker_genes,
                  gene_expression_count)
        return adata

    def generate_h5ad_file(self, adata, h5ad_path):
        md.h5ad_file(adata, h5ad_path)

    def simulate(self,
                 adata,
                 cell_types,
                 chip_size,
                 min_radius,
                 max_radius,
                 num_cells):
        sl.simulate(adata, cell_types, chip_size, min_radius, max_radius, num_cells)


if __name__ == '__main__':
    cell_types = 5
    chip_size = 1.0
    min_radius = 0.005
    max_radius = 0.01
    num_cells = 100
    num_genes = 1000
    num_marker_gene = 10
    gene_expression_count = 125
    h5ad_path = 'simulate_data.h5ad'
    simulator = Simulator()
    adata = simulator.make_data(num_cells, num_genes, cell_types, num_marker_gene, gene_expression_count)
    simulator.generate_h5ad_file(adata, h5ad_path)
    simulator.simulate(adata, cell_types, chip_size, min_radius, max_radius, num_cells)

