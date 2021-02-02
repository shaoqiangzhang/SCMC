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

int main(int argc, const char** argv){

	if(argc!=2){
		cout<<"\ncell normalization and calculate standard deviation of each gene\n******\nUSAGE:\n******\n";
		cout<<argv[0]<<" <expression_file>  > OutputFile\n";
		cout<<endl;
		exit(1);
	}
	ifstream express_file(argv[1]);

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
	//cout<<cell_num<<"\t"<<gene_num<<endl;
	Matrix input_data;
	float *genedev;
	for(int i=0;i<cell_num;i++){//transpose of the Expression Matrix and normalize
		float thiscellsum=0;
		for(int j=0;j<gene_num;j++){
			thiscellsum+=ExpressMatrix[j][i];
		}
		vector<float> thiscellvector;
		for(int j=0;j<gene_num;j++){
			if(thiscellsum>0){
				thiscellvector.push_back(log(((ExpressMatrix[j][i])/thiscellsum)*10000+1));
				//input_data[i][j]=log(((ExpressMatrix[j][i])/thiscellsum)*10000+1);//normalize
			}else{
				thiscellvector.push_back(0);
				//input_data[i][j]=0;
			}
		}
		input_data.push_back(thiscellvector);
		thiscellvector.clear();
	}
	genedev=(float *)malloc(gene_num*sizeof(float));
	for(int j=0;j<gene_num;j++){
		float sumvalue=0, meanvalue=0, dev=0;
		for(int i=0;i<cell_num;i++){
			sumvalue=sumvalue+input_data[i][j];
		}
		meanvalue=sumvalue/cell_num;
		for(int i=0;i<cell_num;i++){
			dev=dev+(input_data[i][j]-meanvalue)*(input_data[i][j]-meanvalue);
		}
		genedev[j]=sqrt(dev);
		cout<<genedev[j]<<endl;
	}

}