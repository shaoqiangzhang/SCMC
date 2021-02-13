# SCMC
Single Cell Massive Clustering

OS: Linux or MacOS

GPU: Nvidia GPU with CUDA toolkit (https://developer.nvidia.com/cuda-toolkit)

Languages: C++, Perl

C++ libraries: CUDA, openMP


##  compile all c++ files and install MCL 
```
g++ read_mtx_file.cpp -o read_mtx_file
g++  cal_sd.cpp -o cal_sd
g++ data_preprocess.cpp -o data_preprocess
nvcc cal_sim_knn.cu -o cal_sim_knn
g++ -fopenmp all_sim_knn_openmp.cpp -o all_sim_knn_openmp
```

## Format of the input expression file 
```
Gene*Cell	Cell1	Cell2	Cell3	Cell4	Cell5	Cell6
ENSMUSG00000000001	44	67	14	43	55	43
ENSMUSG00000000031	0	0	5	7	0	0
... ... ... ...
```
## If the file is ".mtx" format, please transfer mtx file into the above format

For example, give a file named "matrix.mtx" as follows,
```
%%MatrixMarket matrix coordinate integer general
%
33694 4340 5727695
33665 1 3
33663 1 7
33662 1 11
33661 1 1
33660 1 5
33659 1 17
... ...
```
transfer "matrix.mtx" into a file with name "matrix.txt"
```
./read_mtx_file matrix.mtx > matrix.txt
```

## step 1: normalize cells and calculate standard_deviation and select variable genes


### run the follow command to calculate standard deviations of all genes
```
./cal_sd <input_expression_file>  >data.sd
```
### calculate the standard_deviation distribution of all genes
```
perl distribution_calculation.pl data.sd 1
```
## Step 2: data preprocess (Normalize cells & select top variable genes)

### run the follow command to select genes with a standard deviation cutoff
```
./data_preprocess <input_expression_file> <standard_deviation_cutoff> >new_expression_file
```
## step 3: calculate cell-cell similarity using multiple CPU threads + GPU (CUDA)

### run the following command to calculate similarity among all cells (with similarity cut-off and KNN)
```
./all_sim_knn_openmp <new_expression_file> <cell_number> <K_of_KNN> <similarity_cutoff> <thread_number>
```

## step 4: MCL clustering

### run the MCL program
```
mcl <similarity_file> --abc -o output_cluster_file
```
