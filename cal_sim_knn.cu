//calculate cosine similarity from a cell to all of the other cells parallelized by CUDA
#include<iostream>
#include<iomanip>
#include<vector>
#include<fstream>
#include<sstream>
#include<map>
#include<cmath>
#include<algorithm>
#include<numeric>
#include <cstdlib>
#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include<time.h>
#include <cuda_runtime.h>

using namespace std;
typedef vector<vector<float> > Matrix;
//read a scRNA-seq expression table into a matrix
typedef pair<int, float> intfloatPAIR;

struct CmpByValue {
	bool operator()(const intfloatPAIR& lhs, const intfloatPAIR& rhs) {
		return lhs.second > rhs.second;
	}
};



__global__ void cal_dis(float *train_data, float *test_data, int D, int N, float *dis,int pitch)
{//D is the dimension(#genes) of vector(cell), N is #cell
	//int tid = blockIdx.x;
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	if(tid<N){
		float traintemp = 0;
		float testtemp =0;
		float dottemp =0;
		float trainsum = 0;
		float testsum = 0;
		float dotsum = 0;
		
		for(int i=0;i<D;i++){
			traintemp = *((float*)((char*)train_data + tid * pitch) + i);
			testtemp =test_data[i];
			dottemp = *((float*)((char*)train_data + tid * pitch) + i) * test_data[i];
			dotsum += dottemp;
			trainsum += traintemp*traintemp;
			testsum += testtemp*testtemp;
		}
		dis[tid] = dotsum/(sqrt(testsum)*sqrt(trainsum));
	}
}

void print(float *data,int start_i, int n, float cutoff)
{
   for(int i=start_i+1;i<n;i++)
    {
		if(data[i]>=cutoff){
			cout<<start_i<<"\t"<<i<<"\t"<<data[i]<<"\n";
		}
    }
}

void knn(float *data,int start_i, int n, int k)
{//k-nearest neighbors of start_i
    vector<intfloatPAIR> i2fVector;
	i2fVector.clear();
	for(int i=0;i<n;i++)
    {
        if(i!=start_i){
			intfloatPAIR ifp=make_pair(i,data[i]);
			if(i2fVector.size()<k){
				i2fVector.push_back(ifp);
			}else{
				sort(i2fVector.begin(), i2fVector.end(), CmpByValue());//sort second in descending order
				if(i2fVector[i2fVector.size()-1].second<data[i]){//if the last one less than data[i],replace it by data[i]
					i2fVector.pop_back(); //delete the min(first) one
					i2fVector.push_back(ifp);
				}
			}
		}
    }
	for (int i=0;i<i2fVector.size();i++)
	{//print the KNNs of start_i
		cout<<start_i<<"\t"<<i2fVector[i].first<<"\t"<<i2fVector[i].second<<"\n";
	}
}





int main(int argc, const char** argv){
	int the_cell=0;
	int K=10;
	float cutoff=0.8;
	if(argc!=5){
		cout<<"\nCalculate cosine simialrity scores between a cell and each of the others\n******\nUSAGE:\n******\n";
		cout<<argv[0]<<" <expression_file> <this cell>  <K>  <cutoff>  > OutputFile\n";
		cout<<"\n Note: <this cell> is order of a cell, e.g. 0,1,2,... \n";
		cout<<"K is an integral number (K-nearest neighbors). e.g. K=5\n";
		cout<<"cutoff: is the normalized cell-cell simialrity score between 0 and 1. e.g. cutoff=0.95\n"; 
		cout<<endl;
		exit(1);
	}
	ifstream express_file(argv[1]);

	string scell(argv[2]);
	istringstream iscell(scell);
	iscell>>the_cell;// current cell label
	
	string sk(argv[3]);
	istringstream isk(sk);
	isk>>K;// K of KNN
	
	string scut(argv[4]);
	istringstream iscut(scut);
	iscut>>cutoff;

	string firstline; getline(express_file,firstline);//read the first line
	Matrix ExpressMatrix; //define a matrix to store the expression values
	string GeneName;
	vector<string> allGeneNames;
	for(string s;getline(express_file,s);){//read the input file into a matrix
		vector<float> thisline;
		istringstream sin(s);
		sin>>GeneName;
		allGeneNames.push_back(GeneName);
		for(float a; sin>>a; ){
			thisline.push_back(a);
		}
		ExpressMatrix.push_back(thisline);
		thisline.clear();
	}
	int cell_num=ExpressMatrix[0].size();
	int gene_num=ExpressMatrix.size();
	//int gene_num=allGeneNames.size();
	//cout<<"#cells="<<cell_num<<"; #genes="<<gene_num<<endl;
	
	float *h_train_data, *h_test_data ;
	
	h_train_data=(float *)malloc(gene_num*cell_num*sizeof(float));
	h_test_data=(float *)malloc(gene_num*sizeof(float));
	
	for(int i=0;i<cell_num;i++){//transpose of the Expression Matrix and normalize
		/*float thiscellsum=0;
		for(int j=0;j<gene_num;j++){
			thiscellsum+=ExpressMatrix[j][i];
		}*/
		for(int j=0;j<gene_num;j++){
			h_train_data[i*gene_num+j]=ExpressMatrix[j][i];
			//h_train_data[i][j]=log(((ExpressMatrix[j][i])/thiscellsum)*10000+1);//normalize
		}
	}
	int N=cell_num;
	float distance[N];
	float *d_train_data , *d_test_data , *d_dis;

	size_t pitch_d;
	size_t pitch_h = gene_num * sizeof(int) ;


	//allocate memory on GPU 
	cudaMallocPitch( &d_train_data , &pitch_d , gene_num * sizeof(float) , N ); 
	cudaMalloc( (void**)&d_test_data ,  gene_num*sizeof(float) );
	cudaMalloc( (void**)&d_dis , N*sizeof(float) );
		

	for(int j=0;j<gene_num;j++){
		h_test_data[j]=h_train_data[the_cell*gene_num+j];//
	}
	
		
	//copy training and testing data from host to device
	cudaMemcpy2D( d_train_data , pitch_d , h_train_data , pitch_h , gene_num * sizeof(float) , N , cudaMemcpyHostToDevice );
	cudaMemcpy( d_test_data,  h_test_data ,  gene_num*sizeof(float), cudaMemcpyHostToDevice);

	int threadPerBlock =64;
	int blocksPerGrid=(N+ threadPerBlock-1)/threadPerBlock;


	//calculate the distance
	cal_dis<<<blocksPerGrid, threadPerBlock>>>( d_train_data,d_test_data,gene_num,N,d_dis,pitch_d );

	//copy distance data from device to host
	cudaMemcpy( distance , d_dis  , N*sizeof(float) , cudaMemcpyDeviceToHost);

	//cout<<"distance:"<<endl;;
	print(distance , the_cell, cell_num, cutoff);
	
	knn(distance, the_cell,cell_num, K);
	
	cudaFree(d_train_data);
	cudaFree(d_test_data);
	cudaFree(d_dis);
	
	return 0;
}
