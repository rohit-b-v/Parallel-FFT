#include <stdio.h>
#include <mpi.h> 
#include <complex.h> 
#include <stdlib.h>
#include <math.h>	

#define PI 3.14159265
#define n 16384 

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void FillArray(double table[][3])
{
    for (int i = 0; i < n; i++)
    {
        table[i][0]=i;
        table[i][1]=fRand(0,2);
        table[i][2]=fRand(0,2);
    }
}

int main()
{
	int my_rank,size;
	MPI_Init(NULL,NULL); 
	MPI_Comm_size(MPI_COMM_WORLD,&size);   
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);   
	double start,finish;
	
	FILE *outfile;
	if(my_rank == 0) 
	{
		outfile = fopen("MPI_output_parallel.txt", "w"); 
	}
	
		if(my_rank == 0) 
		{	
			start = MPI_Wtime();
		}
	     
		double complex evenpart[(n / size / 2)];
		double complex oddpart[(n / size / 2)]; 
		double complex evenpartmaster[ (n / size / 2) * size];
		double complex oddpartmaster[ (n / size / 2) * size]; 
		double storeKsumreal[n]; 
		double storeKsumimag[n]; 
		
		double subtable[(n / size)][3];  
		
        double table[n][3];
        FillArray(table);
		int sendandrecvct = (n / size) * 3; 
		MPI_Scatter(table,sendandrecvct,MPI_DOUBLE,subtable,sendandrecvct,MPI_DOUBLE,0,MPI_COMM_WORLD); 
		for (int k = 0; k < n / 2; k++)  
		{
					
			double sumrealeven = 0.0; 
			double sumimageven = 0.0; 
			double sumrealodd = 0.0; 
			double sumimagodd = 0.0; 
			
			for(int i = 0; i < (n/size)/2; i++) 
			{
				double factoreven , factorodd = 0.0;
				int shiftevenonnonzeroP = my_rank * subtable[2*i][0]; 
				int shiftoddonnonzeroP = my_rank * subtable[2*i + 1][0]; 
				
				
				double realeven = subtable[2*i][1]; 
				double complex imaginaryeven = subtable[2*i][2]; 
				double complex componeeven = (realeven + imaginaryeven * I); 
				if(my_rank == 0) 
				{
					factoreven = ((2*PI)*((2*i)*k))/n; 
				}
				else 
				{
					factoreven = ((2*PI)*((shiftevenonnonzeroP)*k))/n; 
				}
				double complex comptwoeven = (cos(factoreven) - (sin(factoreven)*I)); 
				
				evenpart[i] = (componeeven * comptwoeven); 
				
				
				double realodd = subtable[2*i + 1][1]; 
				double complex imaginaryodd = subtable[2*i + 1][2]; 
				double complex componeodd = (realodd + imaginaryodd * I); 
				if (my_rank == 0)
				{
					factorodd = ((2*PI)*((2*i+1)*k))/n;
				}
				else 
				{
					factorodd = ((2*PI)*((shiftoddonnonzeroP)*k))/n;
				}
							
				double complex comptwoodd = (cos(factorodd) - (sin(factorodd)*I));
				oddpart[i] = (componeodd * comptwoodd); 
				
			}
			
			MPI_Gather(evenpart,(n / size / 2),MPI_DOUBLE_COMPLEX,evenpartmaster,(n / size / 2), MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);
			MPI_Gather(oddpart,(n / size / 2),MPI_DOUBLE_COMPLEX,oddpartmaster,(n / size / 2), MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD);

			if(my_rank == 0)
			{
				for(int i = 0; i < (n / size / 2) * size; i++) 
				{
					sumrealeven += creal(evenpartmaster[i]); 
					sumimageven += cimag(evenpartmaster[i]); 
					sumrealodd += creal(oddpartmaster[i]); 
					sumimagodd += cimag(oddpartmaster[i]); 
				}
				storeKsumreal[k] = sumrealeven + sumrealodd; 
				storeKsumimag[k]  = sumimageven + sumimagodd; 
				storeKsumreal[k + n/2] = sumrealeven - sumrealodd; 
				storeKsumimag[k + n/2] = sumimageven - sumimagodd; 
                if(k <= n) 
				{
					if(k == 0)
					{
						fprintf(outfile," \n\n TOTAL PROCESSED SAMPLES : %d\n",n);
					}
					fprintf(outfile,"================================\n");
					fprintf(outfile,"XR[%d]: %.4f XI[%d]: %.4f \n",k,storeKsumreal[k],k,storeKsumimag[k]);
					fprintf(outfile,"================================\n");
				}
			}
		}
		if(my_rank == 0)
		{
			finish=MPI_Wtime(); 
			double timeElapsed = finish-start;  
			fprintf(outfile,"Time Elaspsed %f Seconds\n",timeElapsed);
		}
	
	
	MPI_Barrier(MPI_COMM_WORLD); 
	MPI_Finalize(); 
	return 0;
}