//calculate cosine similarity from a cell to all of the other cells parallelized by CUDA
#include<iostream>
#include<iomanip>
#include<vector>
#include<fstream>
#include<sstream>
#include <cstdlib>
#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include<time.h>

using namespace std;
typedef vector<vector<float> > Matrix;
//read a scRNA-seq expression table into a matrix

__global__ void cal_dis(float *train_data, float *test_data, int startlabel, int D, int N, float *dis,int pitch)
{//D is the dimension(#genes) of vector(cell), N is #cell
	int tid = blockIdx.x;
	if(tid<startlabel){
		dis[tid]=0;
	}else if(tid<N){
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

void print(float *data,int start_i, int n)
{
    for(int i=start_i;i<n;i++)
    {
        cout<<start_i<<"\t"<<i<<"\t"<<data[i]<<"\n";
    }
    //cout<<endl;
}


int main(int argc, const char** argv){
	int the_cell=0;
	if(argc!=3){
		cout<<"\nCalculate cosine simialrity scores between a cell and each of the others\n******\nUSAGE:\n******\n";
		cout<<argv[0]<<" <expression_file> <this cell>  > OutputFile\n";
		cout<<endl;
		exit(1);
	}
	ifstream express_file(argv[1]);

	string scell(argv[2]);
	istringstream iscell(scell);
	iscell>>the_cell;// current cell label

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

	//calculate the distance
	cal_dis<<<N,1>>>( d_train_data,d_test_data,the_cell,gene_num,N,d_dis,pitch_d );

	//copy distance data from device to host
	cudaMemcpy( distance , d_dis  , N*sizeof(float) , cudaMemcpyDeviceToHost);

	//cout<<"distance:"<<endl;;
	print(distance , the_cell, cell_num);

	cudaFree(d_train_data);
	cudaFree(d_test_data);
	cudaFree(d_dis);
	
	return 0;
}
