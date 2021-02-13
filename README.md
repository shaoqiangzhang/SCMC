# SCMC
Single Cell Massive Clustering

OS: Linux or MacOS
GPU: Nvidia GPU
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
