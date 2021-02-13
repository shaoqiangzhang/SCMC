// compute cell-cell similarity using CPU's openMP and GPU's CUDA
#include<fstream>
#include<iostream>
#include<sstream>
#include<ctime>
#include<cstdlib>
#include<omp.h>
#include <stdio.h>
#include<cstring>

using namespace std;

int main(int argc, char** argv){
	string file;
	int cellnum=0;
	int K=10;
	double cutoff=0.8;
	int threadnum=1;
	
	if(argc!=6 ){
		cout<<"\n******\nUSAGE:\n******\n";
		cout<<argv[0]<<" <Input_File> <cell_number> <K_of_KNN> <similarity_cutoff> <number_of_threads> \n";
		cout<<endl;
		exit(1);
	}
	string sfilename(argv[1]);
	istringstream isfilename(sfilename);
	isfilename>>file;
	
	string scell(argv[2]);
	istringstream iscell(scell);
	iscell>>cellnum;//number of cells
	
	string sk(argv[3]);
	istringstream isk(sk);
	isk>>K;// K of KNN
	
	string scut(argv[4]);
	istringstream iscut(scut);
	iscut>>cutoff;
	
	string sthread(argv[5]);
	istringstream isthread(sthread);
	isthread>>threadnum;
	
	time_t tstart,tend;
	time(&tstart);//start to keep running time
	
	omp_set_num_threads(threadnum); 
	
	#pragma omp parallel for 
	for(int i=0;i<threadnum;i++)
	{
		for (int j=0;j<cellnum;j++){
			if(j%threadnum== i){
				ostringstream oss;
				oss<<"./cal_sim_knn "<<file<<" "<<j<<" "<<K<<" "<<cutoff<<" >>"<<file<<".knn"<<K<<".cut"<<cutoff<<".sim"<<i;
				cout<<oss.str()<<"\n";
				string commandline;
				commandline=oss.str();
				system(commandline.c_str());
			}
		}
	}

	ostringstream cat_outfiles;
	cat_outfiles<<"cat "<<file<<".knn"<<K<<".cut"<<cutoff<<".sim* >"<<file<<".K"<<K<<".C"<<cutoff<<".sim";
	string cline;
	cline=cat_outfiles.str();
	system(cline.c_str());
	
	ostringstream rm_outfiles;
	rm_outfiles<<"rm -f "<<file<<".knn"<<K<<".cut"<<cutoff<<".sim*";
	string cdline;
	cdline=rm_outfiles.str();
	system(cdline.c_str());
	
	time(&tend);
	double dif = difftime (tend,tstart);
	cout<<"\n\nTotal running time: "<<dif<<" seconds\n"<<endl;
	
	
}