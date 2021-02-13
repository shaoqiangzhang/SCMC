#include<iostream>
#include<iomanip>
#include<vector>
#include<fstream>
#include<sstream>
#include <cstdlib>
#include<stdio.h>
#include<stdlib.h>
#include <malloc.h>

using namespace std;
typedef vector<vector<int> > Matrix;

int main(int argc, const char** argv){
	int row,col;
	Matrix m;
	if(argc!=2){
		cout<<"\nTransfer mtx file into row*column matrix\n******\nUSAGE:\n******\n";
		cout<<argv[0]<<" <mtx_file>  > OutputFile\n";
		cout<<endl;
		exit(1);
	}
	ifstream express_file(argv[1]);

	string firstline; getline(express_file,firstline);//read the first line
	string secondline; getline(express_file,secondline);//read the second line
	string thirdline; getline(express_file,thirdline);//read the thirdline line
	istringstream sin(thirdline); sin>>row; sin>>col;

	for(int i=0;i<row+1;i++){//initialize matrix
		vector<int> this_row;
		for(int j=0;j<col+1;j++){
			this_row.push_back(0);
		}
		m.push_back(this_row);
		this_row.clear();
	}
	
	for(string s;getline(express_file,s);){//read the input file into a matrix
		istringstream sline(s);
		int i,j,exprvalue;
		sline>>i; sline>>j; sline>>exprvalue;
		m[i][j]=exprvalue;
	}
	
	cout<<"gene*cell";
	for(int j=1;j<col+1;j++){//print the j-th cell
		cout<<"\t"<<j;
	}
	cout<<endl;
	
	for(int i=1; i<row+1; i++)
    {
		cout<<i; //the i-th gene
        for(int j=1; j<col+1; j++)
        {
            cout <<"\t"<< m[i][j];
        }
        cout<<endl;
    }
	
}