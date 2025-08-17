#include <stdio.h>
#include <complex.h> 
#include <math.h>	
#include <stdlib.h>
#include <time.h> 
#include <stdlib.h>
#include <string.h>



#define PI 3.14159265
#define n 8192 

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
	double avgtime = 0;
	FILE *outfile;
	outfile = fopen("serial_output.txt", "w"); 

		double start,finish; 
		clock_t t;
   		t = clock(); 
							  
	    double table[n][3];
		double complex evenpart[n/2]; 
		double complex oddpart[n/2]; 
		double storeKsumreal[n]; 
		double storeKsumimag[n]; 
		int k, i ,j; 
		FillArray(table);

		for (k = 0; k < n / 2; k++ ) 
		{	
			
			double sumrealeven = 0.0; 
			double sumimageven = 0.0;
			double sumrealodd = 0.0; 
			double sumimagodd = 0.0; 
			
			for (i = 0; i <= (n/2 - 1); i++) 
			{
				/* -------- EVEN PART -------- */
				double realeven = table[2*i][1]; 
				double complex imaginaryeven = table[2*i][2]; 
				double complex componeeven = (realeven + imaginaryeven * I); 

				double factoreven = ((2*PI)*((2*i)*k))/n; 
				
				double complex comptwoeven = (cos(factoreven) - (sin(factoreven)*I)); 

				evenpart[i] = (componeeven * comptwoeven); 
				
				/* -------- ODD PART -------- */
				double realodd = table[2*i + 1][1]; 
				double complex imaginaryodd = table[2*i + 1][2]; 
				double complex componeodd = (realodd + imaginaryodd * I); 
				
				double factorodd = ((2*PI)*((2*i+1)*k))/n;
															
				double complex comptwoodd = (cos(factorodd) - (sin(factorodd)*I));

				oddpart[i] = (componeodd * comptwoodd);
				
			}
			
			for(i = 0; i < n/2; i++) 
			{
				sumrealeven += creal(evenpart[i]); 
				sumimageven += cimag(evenpart[i]);
				
				sumrealodd += creal(oddpart[i]);
				sumimagodd += cimag(oddpart[i]); 
			}
			
			storeKsumreal[k] = sumrealeven + sumrealodd; 
			storeKsumimag[k] = sumimageven + sumimagodd; 
			
			storeKsumreal[k + n/2] = sumrealeven - sumrealodd; 
			storeKsumimag[k + n/2] = sumimageven - sumimagodd; 
			if(k <= n) 
			{
				if(k == 0)
				{
					fprintf(outfile," \n\nTOTAL PROCESSED SAMPLES : %d\n",n);
				}
				fprintf(outfile,"================================\n");
				fprintf(outfile,"XR[%d]: %.4f XI[%d]: %.4f \n",k,storeKsumreal[k],k,storeKsumimag[k]);
				fprintf(outfile,"================================\n");
			}
		}
		t = clock() - t; 
		double time_taken = ((double)t)/CLOCKS_PER_SEC;
		avgtime = avgtime + time_taken; 
		fprintf(outfile,"Time Elaspsed : %f Seconds\n",time_taken);

	fclose(outfile); 
	return 0;
}