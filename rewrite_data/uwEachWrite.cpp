#include "mex.h"
#include "matrix.h"
#include "stdio.h"

const int NAMELEN = 100;

//////////////////////////////////////////////////////////////////////////////
// uwEachWrite(filename, Ang_plot, RP_dBm)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    (void) nlhs; /* unused parameter */
    // Input arrays
    const mxArray *DATA     = prhs[0];
    char          *filename = (char *)mxGetData(DATA);
    mwSize        nameLen   = mxGetNumberOfElements(DATA);
    // output file name
    char fname[NAMELEN];
    mwSize i, j;
    for (i=0, j=0; i<nameLen; i=i+1, j=j+2)
    {
        fname[i] = filename[j];
    }
    fname[i]='\0';
    
    // Ang_plot
    DATA                = prhs[1];
    double *Ang_plot    = (double *)mxGetData(DATA);
    int N               = mxGetNumberOfElements(DATA);  // 176 or 144
    
    
    // RP_dBm
    DATA    = prhs[2];
    double *RP_dBm      = (double *)mxGetData(DATA);
    int M               = mxGetM(DATA);     // 2048 or 4096
    int MN              = mxGetNumberOfElements(DATA);
    if (mxGetN(DATA) != N)
    {
        printf("read each error, not equal\n");
        exit(1);
    }
    
    // output file
    FILE *fp;
    fp = fopen(fname, "wb");
    fwrite(&N,  sizeof(int), 1, fp);
    fwrite(&M,  sizeof(int), 1, fp);
    fwrite(Ang_plot, sizeof(double), N, fp);
    fwrite(RP_dBm, sizeof(double), MN, fp);   
    fclose(fp);
}
