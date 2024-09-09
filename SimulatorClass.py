import makeData as md
import simulator as sl
import bin_operation as bo


class Simulator:

    def __init__(self):
        print("调用构造方法")

    def make_data(self,
                  gene_expression,
                  cells,
                  num_marker_gene,
                  gene_expression_count):
        adata = md.combine_cell_gene_with_markers(gene_expression,
                  cells,
                  num_marker_gene,
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
    simulator = Simulator()
    adata = simulator.make_data(100, 1000, 5, 10, 125)
    simulator.generate_h5ad_file(adata, 'simulate_data.h5ad')
    simulator.simulate(adata, 5, 1, 0.005, 0.01, 100)

