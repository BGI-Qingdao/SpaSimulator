# SpaSimulator

SpaSimulator is a universal simulator for simulating spatially resolved transcriptomics (SRT) data.

## Overview

------

SpaSimulator is a Python library developed by BGI-Qingdao for simulating spatial transcriptomics sequencing data. This software tool is capable of generating data with specified spatial distribution, gene count, cell count, and cell types to simulate and explore the process of transcriptome data generation and analysis under different experimental conditions.

## Features

1. Simulate Spatial Transcriptomics Data: SpaSimulator can generate transcriptomics data with a specified spatial distribution. Users can set the gene count, cell count, and cell types to simulate the gene expression patterns in different biological samples.

2. Experimental Design Optimization: Using SpaSimulator, researchers can optimize experimental design and workflow, reducing resource and cost wastage, and improving experimental efficiency. By simulating transcriptomics data under different conditions, users can assess the effectiveness of experimental design and optimize it.

3. Data Processing and Analysis Validation: SpaSimulator can generate data with known molecular distributions and expression patterns for validating and evaluating the accuracy and effectiveness of data processing and analysis methods. Researchers can generate data with known ground truth using the simulator to evaluate the performance of different analysis algorithms and methods.

4. Algorithm and Method Development: SpaSimulator provides an experimental and evaluation platform for algorithm and method development. Researchers can generate transcriptomics data under different conditions using the simulator to test the performance and effectiveness of new algorithms and methods and make improvements and optimizations.

5. Education and Training: SpaSimulator can be used for educational and training purposes to help students and researchers understand and learn the principles and applications of spatial transcriptomics sequencing technologies such as Stereo-seq. Through the simulator, users can simulate and explore the process of transcriptomics data generation and analysis under different experimental conditions, enhancing their understanding and application capabilities of these technologies.

## Installation

To install and use SpaSimulator, follow these steps:

------

Clone the SpaSimulator repository from GitHub using the following command:

```
git clone https://github.com/BGI-Qingdao/SpaSimulator.git
```

Create new project on Pycharm --> VCS --> Get from Version Control --> Input url
<img src="./resource/github.png" alt="image-20240909171407833" style="zoom:50%;" />

### Dependencies:

```
pandas==2.2.2
numpy==1.26.4
anndata
argparse==1.3.0
matplotlib==3.8.4
python>=3.12.4
```

## Usage

------

This tool is used for spatially resolved gene expression simulation, differential gene simulation, and cell distribution simulation. The main functions include:

- Generate gene expression data
- Combine gene and cell information
- Make h5ad file in anndata format
- Simulate cell distribution through different cell types
- Make pictures to study

### Example workflow:

Windows:

Change parameters in SimulatorClass.py

```
if __name__ == '__main__':
    num_cells = 100
    num_genes = 1000
    num_types = 5
    num_marker_gene = 10
    gene_expression_count = 125
    h5ad_path = 'simulate_data.h5ad'
    chip_size_mm = 1
    min_cell_radius = 0.005
    max_cell_radius = 0.01
    simulator = Simulator()
    adata = simulator.make_data(num_cells, num_genes, num_types, num_marker_gene, gene_expression_count)
    simulator.generate_h5ad_file(adata, h5ad_path)
    simulator.simulate(adata, num_types, chip_size_mm, min_cell_radius, max_cell_radius, num_cells)
```

Linux:

Understand parameters in main.py

```
if __name__ == '__main__':
    args = parse_args()
    num_cells = args.num_cells
    num_genes = args.num_genes
    cell_types = args.cell_types
    num_marker_genes = args.num_marker_genes
    gene_expression_count = args.gene_expression_count
    scSimulate_data_h5ad_path = args.scSimulate_data_h5ad_path
    chip_size_mm = args.chip_size_mm
    min_radius_mm = args.min_radius_mm
    max_radius_mm = args.max_radius_mm

    # simulate data
    madata = md.combine_cell_gene_with_markers(num_cells, num_genes, cell_types, num_marker_genes, gene_expression_count)

    # generate h5ad file
    md.h5ad_file(madata, scSimulate_data_h5ad_path)

    # simulate
    sl.simulate(madata, cell_types, chip_size_mm, min_radius_mm, max_radius_mm, num_cells)
```

Input parameters

```
python main.py -nc 100 -ng 1000 -ct 5 -nmg 10 -gec 125 -dp 'simulate_data.h5ad' -cs 1 -minr 0.005 -maxr 0.01
```

### Parameters explanation:

```
('-nc', '--num_cells', meaning='Enter the number of cells to be generated)
('-ng', '--num_genes', meaning='Enter the number of genes to be generated')
('-ct', '--cell_types', meaning='Enter the number of cell types to be generated)
('-nmg', '--num_marker_genes', meaning='Enter the number of marker genes to be generated')
('-gec', '--gene_expression_count', meaning='The average expression of non-zero expression genes')
('-dp', '--scSimulate_data_h5ad_path', default='simulate_data.h5ad', meaning='The data saved in the simulation is in h5ad formate')
('-cs', '--chip_size_mm', type=int, meaning='size of chip, unit is mm')
('-minr', '--min_radius_mm', type=int, meaning='minimum radius of cells, unit is mm')
('-maxr', '--max_radius_mm', type=int, meaning='maximum radius of cells, unit is mm')
```

All results will be save in a h5ad file, default file name is `simulate_data.h5ad`

## Visualization

------

SpaSimulator offers a picture to show cell distribution.

<img src="./resource/result.png" alt="image-20240909171407833" style="zoom:50%;" />

## Contribution
Contributions to SpaSimulator are welcome. If you have any suggestions, bug fixes, or feature requests, please raise an issue in the GitHub repository.

## License
SpaSimulator is released under the GNU General Public License v3.0 license. For more details, refer to the "LICENSE" file in the repository.

## Contact
For further questions or support regarding SpaSimulator, please contact c-mayuelong@genomics.cn and liyao1@genomics.cn.

