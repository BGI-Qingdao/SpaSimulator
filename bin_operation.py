import pandas as pd
import numpy as np
import anndata as ad


# 设定参数

# spatial_transcriptomics_gem = 'simulated_spatial_transcriptomics.gem'
# bin_size = 20 * 500
# 分箱操作
def bin_operation(spatial_transcriptomics_gem, bin_size):

    # 加载CSV文件
    gem_df = pd.read_csv(spatial_transcriptomics_gem, sep='\t')

    # 设定bin的大小（每个bin包含20个spot点，每个spot点间隔500纳米）
    bin_size_x = bin_size  # 10000纳米或10微米
    bin_size_y = bin_size  # 10000纳米或10微米

    # 计算每个spot的bin坐标
    gem_df['bin_x'] = (gem_df['x'] // bin_size_x).astype(int)
    gem_df['bin_y'] = (gem_df['y'] // bin_size_y).astype(int)

    # 合并bin_x和bin_y以创建一个唯一的bin ID
    gem_df['bin_id'] = gem_df['bin_x'].astype(str) + '_' + gem_df['bin_y'].astype(str)

    # 根据bin_id进行分组，并对每组bin_id进行求和
    binned_data = gem_df.groupby('bin_id')['Count'].sum().reset_index()

    # 计算bin的中心坐标（假设每个bin是正方形，且从0开始编号）
    bin_ids = binned_data['bin_id']
    bin_x_coords = bin_ids.apply(lambda x: int(x.split('_')[0]) * bin_size_x + bin_size_x / 2)
    bin_y_coords = bin_ids.apply(lambda x: int(x.split('_')[1]) * bin_size_y + bin_size_y / 2)

    # 将中心坐标添加到binned_data中
    binned_data['spatial_x'] = bin_x_coords
    binned_data['spatial_y'] = bin_y_coords

    # 计算每个bin中的基因数量
    spot_counts_per_bin = gem_df.groupby('bin_id').size().reset_index(name='gene_count')

    # print(spot_counts_per_bin)

    # 创建Anndata对象
    # Anndata的X应该是二维的，但因为只有一个特征（总表达量），需要扩展维度
    X = np.expand_dims(binned_data['Count'].values, axis=1)
    obs_dict = binned_data[['spatial_x', 'spatial_y']].to_dict(orient='list')
    var_dict = {'gene_ids': ['total_expression']}  # 假设只有一个变量，即总表达量

    adata = ad.AnnData(X=X,
                       obs=pd.DataFrame(obs_dict),
                       var=pd.DataFrame(var_dict))
    adata.uns['bin_gene_counts'] = spot_counts_per_bin.set_index('bin_id').to_dict(orient='index')  # 添加到.uns作为未观测到的集合数据

    return adata


def gengerate_bin_h5ad_file(adata, h5ad_file_path):
    # 保存为.h5ad文件
    adata.write(h5ad_file_path)

