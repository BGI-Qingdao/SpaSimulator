import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import makeData as md

def is_overlapping(new_cell,
                   placed_cells,
                   new_radius):
    """检查新细胞是否与已放置的细胞重叠"""
    for placed_cell in placed_cells:
        distance = np.sqrt((new_cell[0] - placed_cell[0]) ** 2 + (new_cell[1] - placed_cell[1]) ** 2)
        if distance < (new_radius + placed_cell[2]):
            return True
    return False

def save_to_gem(cell_positions, filename):
    """将细胞位置保存为GEM格式文件"""
    with open(filename, 'w') as f:
        f.write("CellID\tX\tY\n")  # 写入表头
        for i, position in enumerate(cell_positions):
            if not np.isnan(position[0]) and not np.isnan(position[1]):  # 确保位置有效
                f.write(f"Cell_{i+1}\t{position[0]}\t{position[1]}\n")

def simulate(adata,
             cell_types,
             chip_size,
             min_radius,
             max_radius,
             num_cells):
    """
    模拟细胞在组织芯片上空间分布
    :param adata: 生成的基因表达矩阵
    :param cell_types: 细胞类型数量
    :param chip_size: 芯片尺寸
    :param min_radius: 细胞最小半径，细胞默认为圆形
    :param max_radius: 细胞最大半径，细胞默认为圆形
    :param num_cells: 细胞数量
    :return:
    """
    # 创建一个图形和轴
    fig, ax = plt.subplots()

    # 为每种细胞类型定义颜色
    unique_cell_types = sorted(set(adata.obs['cell_type']))
    colors = plt.cm.hsv(np.linspace(0, 1, len(unique_cell_types)))

    # 存储所有已放置的细胞位置
    placed_cells = []
    cell_positions = []
    cell_radii = []

    # 模拟每种细胞类型的分布
    for i in range(cell_types):
        # 定义每个细胞类型的区域
        lower_bound = i * chip_size / cell_types
        upper_bound = (i + 1) * chip_size / cell_types
        cell_coords = []
        cell_radii_temp = []

        while len(cell_coords) < num_cells // cell_types:
            # 生成细胞坐标，确保在指定的区域内
            new_cell = np.random.multivariate_normal([(lower_bound + upper_bound) / 2, 0.5], [[0.01, 0], [0, 0.01]], 1).flatten()
            new_cell[0] = np.clip(new_cell[0], lower_bound, upper_bound)  # 确保x坐标在指定区域内

            # 随机选择细胞半径
            cell_radius = np.random.uniform(min_radius, max_radius)

            # 检查是否重叠
            if not is_overlapping(new_cell, placed_cells, cell_radius):
                cell_coords.append(new_cell)
                cell_radii_temp.append(cell_radius)
                placed_cells.append((new_cell[0], new_cell[1], cell_radius))

        # 将当前细胞类型的细胞位置和半径添加到总列表中
        cell_positions.extend(cell_coords)
        cell_radii.extend(cell_radii_temp)
        # 绘制细胞
        for coord, cell_radius in zip(cell_coords, cell_radii_temp):
            circle = patches.Circle((coord[0], coord[1]), cell_radius, color=colors[i])
            ax.add_patch(circle)

    # 确保 cell_positions 和 cell_radii 与 adata.obs 的长度相匹配
    num_obs = len(adata.obs)
    cell_positions_padded = [p if p is not None else [np.nan, np.nan] for p in cell_positions]
    cell_positions_array = np.array(cell_positions_padded + ([[np.nan, np.nan]] * (num_obs - len(cell_positions_padded)))[:num_obs])

    cell_radii_padded = [r if r is not None else np.nan for r in cell_radii]
    cell_radii_array = np.array(cell_radii_padded + ([np.nan] * (num_obs - len(cell_radii_padded)))[:num_obs])

    # 添加图例
    legend_handles = [patches.Patch(color=color, label=cell_type) for color, cell_type in
                      zip(colors, unique_cell_types)]
    ax.legend(handles=legend_handles, loc='best')

    # 将细胞位置和半径信息添加到adata中
    adata.obsm['cell_positions'] = cell_positions_array
    adata.obsm['cell_radii'] = cell_radii_array.reshape(num_obs, 1)

    # 设置图形属性
    ax.set_xlim(0, chip_size)
    ax.set_ylim(0, chip_size)
    ax.set_aspect('equal')  # 确保x和y轴的比例相同
    ax.set_title('Simulated Tissue on Chip with Non-overlapping Cells')
    ax.set_xlabel('mm')
    ax.set_ylabel('mm')
    # 保存图形
    plt.savefig('result.png')
    # 显示图形
    plt.show()
    # 保存细胞位置到GEM文件
    save_to_gem(cell_positions, 'cell_positions.gem')
    return adata

cell_types = 5
chip_size = 1.0
min_radius = 0.005
max_radius = 0.01
num_cells = 100
num_genes = 1000
num_marker_gene = 10
gene_expression_count = 125

adata = md.combine_cell_gene_with_markers(num_cells, num_genes, cell_types, num_marker_gene, gene_expression_count)
adata_simulate = simulate(adata, cell_types, chip_size, min_radius, max_radius, num_cells)