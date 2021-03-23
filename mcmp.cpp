#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include <time.h>
using namespace std;

int main(int argc, char **argv)
{
	if(argc<2)
	{
		return 0;
	}
	int N_Thread;
	if(argc<3)
	{
		N_Thread=1;
	}
	else
	{
		string tmp=argv[2];
		N_Thread=atoi(tmp.c_str());
	}
	double N_Inflation;
	if(argc<4){
		N_Inflation=1.5;
	}
	else{
		string tmp2=argv[3];
		N_Inflation=atof(tmp2.c_str());
	}
	const char *File_Table=argv[1];

	int N_Depth=3;
	int N_Power=2;

	//READ TABLE
	ifstream TABLE (File_Table);

	map <pair<string,string>,double> Gene_Map;
	map <string,int> Gene_Count;

	string line;
	while(getline(TABLE,line))
	{
		string Gene1;
		string Gene2;
		string tmp_Score;
		double Score;
		stringstream linestream(line);

		getline(linestream,Gene1,'\t');
		getline(linestream,Gene2,'\t');
		getline(linestream,tmp_Score,'\t');
		Score=atof(tmp_Score.c_str());
		Gene_Count[Gene1]++;
		Gene_Count[Gene2]++;

		Gene_Map[make_pair(Gene1,Gene2)]=Score;
		Gene_Map[make_pair(Gene1,Gene1)]=1;
		Gene_Map[make_pair(Gene2,Gene1)]=Score;
		Gene_Map[make_pair(Gene2,Gene2)]=1;

	}
	//Get All Genes
	const int Total_Gene_Count=Gene_Count.size();
	string All_Genes[Total_Gene_Count];
	map <string,int> Gene2Index;

	int tmp_Gene_Count=0;
	map<string,int>::iterator i;

	for(i=Gene_Count.begin();i!=Gene_Count.end();i++)
	{
		All_Genes[tmp_Gene_Count]=i->first;
		Gene2Index[i->first]=tmp_Gene_Count;
		tmp_Gene_Count++;
	}
	//Make Matrix
	vector<vector<double> > GeneMatrix(Total_Gene_Count,vector<double>(Total_Gene_Count,0));
	//double GeneMatrix[Total_Gene_Count][Total_Gene_Count];
	#pragma omp parallel num_threads(N_Thread)
	for(int i=0;i<Gene_Count.size();i++)
	{
		for(int j=0;j<Gene_Count.size();j++)
		{
			if(Gene_Map.count(make_pair(All_Genes[i],All_Genes[j])))
			{
		
				GeneMatrix[i][j]=Gene_Map[make_pair(All_Genes[i],All_Genes[j])];
		
			}//exitst
			else//not exists
			{
				GeneMatrix[i][j]=0;
		
			}//not exists
		}
		
	}
	
	//Normalization
	//double Matrix_Normal[Total_Gene_Count][Total_Gene_Count];
	vector<vector<double> > Matrix_Normal(Total_Gene_Count,vector<double>(Total_Gene_Count,0));
	#pragma omp for  schedule(dynamic)
	for(int i=0;i<Total_Gene_Count;i++)
	{
		double Sum=0;
		for(int j=0;j<Total_Gene_Count;j++)
		{
			Sum+=GeneMatrix[j][i];
		}
		//#pragma omp for  schedule(dynamic)
		for(int j=0;j<Total_Gene_Count;j++)
		{
			if(Sum==0)
			{
				Matrix_Normal[j][i]=0;
			}
			else
			{
				Matrix_Normal[j][i]=GeneMatrix[j][i]/Sum;
			}
		}
	}
	
	//Iteration
	//cout<<"HI"<<Total_Gene_Count<<endl;
	int Itr=0;
	
	vector<vector<double> > Matrix_Pow(Total_Gene_Count,vector<double>(Total_Gene_Count,0));
	vector<vector<double> > Matrix_Inf(Total_Gene_Count,vector<double>(Total_Gene_Count,0));
	vector<vector<double> > Matrix_New(Total_Gene_Count,vector<double>(Total_Gene_Count,0));

	while(Itr<50)
	//while(Itr<1000)
	{
	#pragma omp for  schedule(dynamic)
		for(int i=0;i<Total_Gene_Count;i++)
		{
			//#pragma omp for  schedule(dynamic)
			for(int j=0;j<Total_Gene_Count;j++)
			{
				double Sum=0;
				for(int k=0;k<Total_Gene_Count;k++)
				{
					Sum+=Matrix_Normal[i][k]*Matrix_Normal[k][j];
				}
				if(Sum<0.001)
				//if(Sum<0.0001)
				{
					Matrix_Pow[i][j]=0;
				}
				else
				{
					Matrix_Pow[i][j]=Sum;
				}
			}
		}
		//Inflate
		#pragma omp for  schedule(dynamic)
		for(int i=0;i<Total_Gene_Count;i++)
		{
			//#pragma omp for  schedule(dynamic)
			for(int j=0;j<Total_Gene_Count;j++)
			{
				Matrix_Inf[i][j]=pow(Matrix_Pow[i][j],N_Inflation);
			}
		}
		//Nomalization again
		#pragma omp for  schedule(dynamic)
		for(int i=0;i<Total_Gene_Count;i++)
		{
			double Sum=0;
			for(int j=0;j<Total_Gene_Count;j++)
			{
				Sum+=Matrix_Inf[j][i];
			}
			//#pragma omp for  schedule(dynamic)
			for(int j=0;j<Total_Gene_Count;j++)
			{
				if(Sum==0)
				{
					Matrix_New[j][i]=0;
				}
				else
				{
					Matrix_New[j][i]=Matrix_Inf[j][i]/Sum;
				}
			}
		}
		//Check Chaos
		double Count_diff=0;
		double Cut_diff=0.0001;
	#pragma omp for  schedule(dynamic)
		for(int i=0;i<Total_Gene_Count;i++)
		{
			for(int j=0;j<Total_Gene_Count;j++)
			{
				Count_diff+=abs(Matrix_New[i][j]-Matrix_Normal[i][j]);
			}
		}
		// Assign New 2 Normal
		#pragma omp for  schedule(dynamic)
		for(int i=0;i<Total_Gene_Count;i++)
		{
			//#pragma omp for  schedule(dynamic)
			for(int j=0;j<Total_Gene_Count;j++)
			{
				Matrix_Normal[i][j]=Matrix_New[i][j];
			}
		}
		if(Count_diff<Cut_diff)//No Chaos
		{
			break;
		}
		Itr++;
	}
	//print Results
	int Check_Used[Total_Gene_Count];
	//Check init
	for(int i=0;i<Total_Gene_Count;i++)
	{
		Check_Used[i]=0;
	}
	for(int i=0;i<Total_Gene_Count;i++)
	{
		int IsCluster=0;
		for(int j=0;j<Total_Gene_Count;j++)
		{
			if(Matrix_Normal[i][j]>0)
			{
				if(!Check_Used[j])
				{
					IsCluster=1;
					cout<<All_Genes[j];
					if(j<Total_Gene_Count-1)
					{
						cout<<"\t";
					}
				}
				Check_Used[j]=1;
			}
		}
		if(IsCluster)
		{
			cout<<endl;
		}
	}
	return 0;
}
