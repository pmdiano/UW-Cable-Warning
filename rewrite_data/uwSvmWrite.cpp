#include "mex.h"
#include "matrix.h"
#include "stdio.h"

const int NAMELEN = 100;

//////////////////////////////////////////////////////////////////////////////
// uwSvmWrite(filename, sv, alpha, bias, shift, scale, b)
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
    
    // sv and size
    DATA                = prhs[1];
    double *sv          = (double*)mxGetData(DATA);
    int M               = mxGetM(DATA);
    int N               = mxGetN(DATA);
    
    // alpha
    DATA                = prhs[2];
    double *alpha       = (double*)mxGetData(DATA);
    
    //bias
    DATA                = prhs[3];
    double *bias        = (double*)mxGetData(DATA);
    
    //shift
    DATA                = prhs[4];
    double *shift       = (double*)mxGetData(DATA);
    
    //scale
    DATA                = prhs[5];
    double *scale       = (double*)mxGetData(DATA);
    
    // b
    DATA                = prhs[6];
    double *b           = (double*)mxGetData(DATA);

    
    // output file
    printf("sizeof(double)=%d", sizeof(double));
    FILE *fp;
    fp = fopen(fname, "wb");
    
    fwrite(&M, sizeof(int), 1, fp);
    fwrite(&N, sizeof(int), 1, fp);
    fwrite(sv, sizeof(double), M*N, fp);
    fwrite(alpha, sizeof(double), M, fp);
    fwrite(bias, sizeof(double), 1, fp);
    fwrite(shift, sizeof(double), N, fp);
    fwrite(scale, sizeof(double), N, fp);
    fwrite(b, sizeof(double), 2*N+1, fp);
  
    fclose(fp);
}
