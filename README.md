# SCMC
Single Cell Massive Clustering

## step 1: normalize cells and calculate standard_deviation and select variable genes

### compile the cpp file
```
g++ -o cal_sd cal_sd.cpp
```
### run the follow command to calculate standard deviations of all genes
```
./cal_sd <input_expression_file>  >data.sd
```
### calculate the standard_deviation distribution of all genes
```
perl distribution_calculation.pl data.sd 1
```
## Step 2: data preprocess (Normalize cells & select top variable genes)

### compile the cpp file
```
g++ data_preprocess.cpp -o data_preprocess
```

### run the follow command to select genes with a standard deviation cutoff
```
./data_preprocess <input_expression_file> <standard_deviation_cutoff> >new_expression_file
```
## step 3: calculate cell-cell similarity

### compile the c++ cuda file 
```
nvcc -o cal_sim_cuda cal_sim_cuda.cu
```
### run the following perl script to calculate similarity among all cells 
```
perl all_sim_cal.pl <new_expression_file> <cell_number>
```
## step 4: select similarity cut-off

```
perl select_sim_cutoff.pl <input_similarity_file> <cutoff>
```
## step 5: MCL clustering

### compile the cpp file
```
g++ MCLustMP.cpp -fopenmp -o mclustmp
```
### run the program
```
./mclustmp input_file <thread_number>  <inflation_number>   > output_file
```
#### note:
1. inflation number 1~2, default=1.5

2. thread number default=1

