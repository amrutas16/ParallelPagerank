#include<iostream>
#include<fstream>
#include<math.h>
#include<string>
using namespace std;

#define modified 0
#define orig 1
#define debug 0

void aitken(float *,float*, float*);
int n = 800;
//int n = 4;
int main()
{
	int i, j, k = 0, count = 0;
	string line;
	int node1, node2;
	char nodeString[15];
	char *nodes;
	char delim[3] = "\t";
	ifstream inputFile;	
	float **matrix = new float *[n];
	int *outgoingLinks = new int[n];
	int *allocation = new int[n];
	size_t sz;
	
	/// For Power method
	float *xk = new float[n];
	float *xkPlus1 = new float[n];
	float *xkMinus1 = new float[n];
	float *temp;
	float delta = 1.0;
	float epsilon = 0.0001;
	float d;
	float norm1 = 0, norm2 = 0;
	float alpha = 0.99;

	for(i=0;i<n;i++){
		
		xk[i] = 1.0/n;
		xkPlus1[i] = 0;
		outgoingLinks[i] = 0;
		allocation[i] = 0;
		//matrix[i] = new float[n];
	}

	if( n == 800)
		inputFile.open("dataset800.txt");
	if(n == 4)
		inputFile.open("dataset4.txt");
	
	while(getline(inputFile, line))
	{
		strcpy(nodeString, line.c_str());
		nodes = strtok(nodeString, delim);
		node1 = stoi(nodes, &sz);
		nodes = strtok(NULL, delim);
		node2 = stoi(nodes, &sz);

		++(outgoingLinks[node1 - 1]);

		if(allocation[node2 - 1] == 0){
			matrix[node2 - 1] = new float[n];
			allocation[node2 - 1] = 1;
		}
		matrix[node2 - 1][node1 - 1] = 1;
	}


	for(i=0;i<n;i++)
	{
		
		if(allocation[i] == 0){
			allocation[i] = 1;
			matrix[i] = new float[n];
		}
		
		for(j=0;j<n;j++)
		{
			///Make rest of the elements 0
			if(matrix[i][j] != 1)
				matrix[i][j] = 0;

			///Dangling nodes
			if(outgoingLinks[j] == 0)
				matrix[i][j] = 1/n;
			else
				matrix[i][j] = matrix[i][j] / outgoingLinks[j];

#if orig
			matrix[i][j] = (alpha * matrix[i][j]) + ((1 - alpha) / n);   ///introduce damping factor
#endif
		}
	}
	
#if debug
	for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
		{
			cout<<matrix[i][j]<<" ";
		}
		cout<<"\n";
	}
#endif
	
	///Power method
	while(delta > epsilon)
	{
		//cout<<"\nDelta = "<<delta<<"\n";
		delta = 0;
		norm1 = 0;
		norm2 = 0;
		++k;

		for(i=0;i<n;i++)
		{
			xkPlus1[i] = 0;
			for(j=0;j<n;j++)
			{
				xkPlus1[i] = xkPlus1[i] + (xk[j] * matrix[i][j]);
			}
#if modified
			xkPlus1[i] = xkPlus1[i] * 0.85;  ///Reqd for modified power method
#endif
			
			norm1 += xk[i];
			norm2 += xkPlus1[i];
		}
		d= norm1 - norm2;
		for(i=0;i<n;i++)
		{
#if modified
			xkPlus1[i] = xkPlus1[i] + (d/n);			///Reqd for modified power method
#endif
			delta = delta + fabs(xkPlus1[i] - xk[i]);
		}
		cout<<"\ndelta = "<<delta;
		
#if debug
		cout<<"\nVector = ";
		for(i=0;i<n;i++)
			cout<<xkPlus1[i]<<" ";
		cout<<"\n";
#endif	
		
		if(k % 3 == 0)
			aitken(xkMinus1, xk, xkPlus1);
		
		///for next iteration
		temp = xkMinus1;
		xkMinus1 = xk;
		xk = xkPlus1;
		xkPlus1 = temp;

	}

#if debug
	cout<<"Pagerank vector = ";
	for(i=0;i<n;i++)
		cout<<xk[i]<<"\n";
#endif
	cout<<"\nNumber of iterations = "<<k;
}

void aitken(float *x0, float *x1, float *x2)
{
	float *temp = x2;
	float gi, hi, fi;
	int i;
	for(i=0;i<n;i++)
	{
		gi = (x1[i] - x0[i]) * (x1[i] - x0[i]);
		hi = x2[i] - 2*x1[i] + x0[i];
		x2[i] = x2[i] - (gi/hi);
	}
}
