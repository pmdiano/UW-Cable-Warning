#include "mex.h"
#include "matrix.h"
#include "stdio.h"

const int NAMELEN = 100;

//////////////////////////////////////////////////////////////////////////////
// uwBaseWrite(filename, Ang_endpt, MagMax, MagMin, fft_dR)
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
    
    // Ang_endpt
    DATA                = prhs[1];
    double *data        = (double *)mxGetData(DATA);
    double Ang_endpt    = data[0];
    
    // MagMax
    double MagMax, MagMin, fft_dR;
    DATA    = prhs[2];
    data    = (double *)mxGetData(DATA);
    MagMax  = data[0];
    
    // MagMin
    DATA    = prhs[3];
    data    = (double *)mxGetData(DATA);
    MagMin  = data[0];
    
    // fft_dR
    DATA    = prhs[4];
    data    = (double *)mxGetData(DATA);
    fft_dR  = data[0];
    
    // output file
    FILE *fp;
    fp = fopen(fname, "wb");
    fwrite(&Ang_endpt, sizeof(double), 1, fp);
    fwrite(&MagMax   , sizeof(double), 1, fp);
    fwrite(&MagMin   , sizeof(double), 1, fp);
    fwrite(&fft_dR   , sizeof(double), 1, fp);
    fclose(fp);
}
