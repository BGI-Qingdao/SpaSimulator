import numpy as np
import anndata as ad


def generate_gene_expression(num_genes,
                  num_cells,
                  gene_expression_count):
    """
    生成服从零膨胀负二项分布的基因表达量
    :param num_genes: 基因数量
    :param num_cells: 细胞数量
    :return:
    gene_expression: 基因表达矩阵，维度为num_genes * num_cells
    """
    # 非零表达基因数量为总基因数量的10%-30%
    bio_p_min = 0.09
    bio_p_max = 0.29

    # 为每个基因生成一个独立的表达概率
    p_value = np.random.uniform(bio_p_min, bio_p_max, size=num_genes)

    # 初始化基因表达矩阵，默认为0（对于未表达的基因）
    gene_expression = np.zeros((num_genes, num_cells), dtype=int)

    # 对于每个基因，生成一个长度为num_cell的二项分布数组，然后填充到gene_expression中
    for i in range(num_genes):
        successes = np.random.binomial(n=1, p=p_value[i], size=num_cells)
        gene_expression[i, :] = successes

    # 找到所有非零表达的元素的扁平化索引
    flat_indices = np.flatnonzero(gene_expression)

    # 对非零表达的基因，生成符合正态分布的表达量，并限制在1到250之间
    num_nonzero = len(flat_indices)
    random_values = np.random.normal(loc=gene_expression_count / 2, scale=gene_expression_count / 5, size=num_nonzero)
    random_values = np.clip(random_values, 1, gene_expression_count)

    # 使用扁平化索引来更新gene_expression中的值
    gene_expression.ravel()[flat_indices] = random_values

    return gene_expression


def generate_cell_type(num_types):
    """
    生成细胞类型
    :param num_types: 期望生成细胞类型的数量
    :return:
    cell_types: 所有细胞类型的名称
    """
    cell_types = []
    # 生成细胞类型
    for i in range(num_types):
        cell_types.append(f"cell_type_{i + 1}")
    return cell_types


def generate_cell(num_cells, num_types):
    """
    生成细胞ID以及每个细胞的类型
    :param num_cells:
    :param num_types:
    :return:
    """
    cells = []
    cell_types = generate_cell_type(num_types)
    for i in range(num_cells):
        # 随机选择一个细胞类型
        cell_type = np.random.choice(cell_types)
        # 生成细胞ID
        cell_id = f"cell_{i + 1}"
        # 将细胞ID和类型作为一个元组添加到列表中
        cells.append((cell_id, cell_type))
    return cells


def select_marker_genes(gene_expression,
                        cells,
                        num_degs,
                        gene_expression_count):
    """
    选取marker gene
    :param gene_expression: 基因表达量矩阵，不包含marker gene
    :param cells: 细胞，包含细胞ID和细胞类型
    :param num_degs: 选取DEGs的数量
    :param gene_expression_amount: 基因表达量的平均值
    :return:
    gene_expression: 包含marker gene的基因表达集
    """
    # cells 是一个包含 (cell_id, cell_type) 元组的列表
    cell_types = np.array([cell[1] for cell in cells])
    cell_indices_by_type = {cell_type: np.where(cell_types == cell_type)[0] for cell_type in np.unique(cell_types)}

    selected_genes = set()
    marker_genes = {}
    for cell_type, indices in cell_indices_by_type.items():
        gene_mean_expression = np.mean(gene_expression[:, indices], axis=1)
        available_genes = np.arange(gene_mean_expression.size)
        available_genes = available_genes[~np.isin(available_genes, selected_genes)]
        num_to_select = min(num_degs, len(available_genes))
        if num_to_select == 0:
            continue
        top_genes_indices = np.argsort(gene_mean_expression[available_genes])[::-1][:num_to_select]
        top_genes = available_genes[top_genes_indices]
        selected_genes.update(top_genes)
        marker_genes[cell_type] = top_genes
        for gene in top_genes:
            gene_expression[gene, indices] = np.random.randint(gene_expression_count * 2,
                                                               gene_expression_count * 3 + 1, size=len(indices))
    return gene_expression


def combine_cell_gene_with_markers(num_cells,
                                   num_genes,
                                   num_types,
                                   num_degs,
                                   gene_expression_count):
    """

    :param num_cells: 细胞数量
    :param num_genes: 基因数量
    :param num_types: 细胞类型数量
    :param num_degs: DEGs数量
    :param gene_expression_count: 基因表达量均值
    :return:
    """
    cells = generate_cell(num_cells, num_types)
    gene_expression = generate_gene_expression(num_genes, num_cells,gene_expression_count)
    gene_expression = select_marker_genes(gene_expression, cells, num_degs, gene_expression_count)

    adata = ad.AnnData(X=gene_expression.T)
    adata.obs_names = [cell[0] for cell in cells]
    adata.obs['cell_type'] = [cell[1] for cell in cells]
    adata.var_names = ['Gene_' + str(i + 1) for i in range(num_genes)]

    return adata


def h5ad_file(adata, h5ad_path):
    """生成h5ad文件"""
    adata.write_h5ad(h5ad_path)
    return h5ad_path

num_cells = 100
num_genes = 1000
cell_types = 5
num_degs = 10
gene_expression_count = 125

adata = combine_cell_gene_with_markers(num_cells, num_genes, cell_types, num_degs, 125)
print(adata)