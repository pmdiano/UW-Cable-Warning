#include "inc.h"
using namespace std;

double PTHETA[] = { -90, -89, -88, -87, -86, -85, -84, -83, -82, -81,
                    -80, -79, -78, -77, -76, -75, -74, -73, -72, -71,
                    -70, -69, -68, -67, -66, -65, -64, -63, -62, -61,
                    -60, -59, -58, -57, -56, -55, -54, -53, -52, -51,
                    -50,
                    50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
                    60, 61, 62, 63, 64, 65, 66, 67, 68, 69,
                    70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                    80, 81, 82, 83, 84, 85, 86, 87, 88, 89  };
int PTHETA_LEN = 81;
int PTHETA_END = 50;


double THETA[] = {  -90, -89, -88, -87, -86, -85, -84, -83, -82, -81,
                    -80, -79, -78, -77, -76, -75, -74, -73, -72, -71,
                    -70, -69, -68, -67, -66, -65, -64, -63, -62, -61,
                    -60,
                    60, 61, 62, 63, 64, 65, 66, 67, 68, 69,
                    70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
                    80, 81, 82, 83, 84, 85, 86, 87, 88, 89  };
int THETA_LEN = 61;
double TH_GATE = -72;
double REGIONFILL_T = -40;
int SLICE_HEIGHT = 128;
int SLICE_WIDTH = 64;
int LOWBOUND_PCT = 68;
double TH_PCT  = 97.0;
double TH_PCT2 = 99.7;
int HGATE_FACTOR = 12;
int size_lo = SIZE_LO_BASE;
int size_hi = SIZE_HI_BASE;
int tower_removed = 0;


extern double rhoScore, rhoScorePrev;
extern int rhoNum, rhoNumPrev;
extern cableLine *rhoCableLine;
extern cableLine *rhoCablePrev;
extern int TOTAL_START;
extern int FRAME_I;
extern FILE* pRhoFile;

extern bool g_saveImage;

int X0 = 0;
extern FILE *pTheta;
extern FILE *pFirstRho;

extern char threshName[];
extern char transName[];
extern char houghName[];
extern char thetaDetectName[];
extern char rhoDetectName[];

extern void showRtrackerImage(double *RP_dBm, int *BW, double *Ang_plot, int M, int N, int x0,
    double MagMax, double MagMin, rtracker *rtrackers_, int totalNum, char *nameToSave, bool showCable);

extern int rhoTrackerToCable(rtracker *rtrackers_, int totalNum, cableLine **cls);

// theta particle stuff
extern bool theta_first_frame;
extern tparticle *theta_particles;
extern double theta_previous;

// rho tracker stuff
extern rtracker *g_rtrackers_;

typedef ParticleFilter<ThetaParticle, double, ThetaObservation> ThetaParticleFilter;
unique_ptr<ThetaParticleFilter> tpf;

//////////////////////////////////////////////////////////////////////////////
// slice-thresholding related function
//////////////////////////////////////////////////////////////////////////////
// slice thresholding
void sliceThreshold(double *Bscope, // input Bscope image
                    int M,              // number of rows
                    int N,              // number of columns
                    int *BWBscope,      // output black-n-white image, shall be pre-allocated outside
                    int sliceHeight,    // slice height
                    double percentage,  // percentage of threshold
                    double threshGate,  // threshold gate
                    int y1,             // only process from y1 to y2
                    int y2,
                    int sliceWidth)
{
    if (sliceWidth)
    {
        blockThreshold(Bscope, M, N, BWBscope, sliceHeight, sliceWidth, percentage, threshGate, y1, y2);
        return;
    }

    int numOfSlice, h1, h2, i, ii, j;
    double numOfSliceD, T;

    // get number of slices
    numOfSliceD = (double)M / (double)sliceHeight;
    numOfSlice = (int)(numOfSliceD);
    if (numOfSlice < numOfSliceD)
    {
        numOfSlice++;
    }

    // process each slice
    for (i=0; i<numOfSlice; i++)
    {
        // get the start and end of this slice
        h1 = i*sliceHeight;
        h2 = (i+1)*sliceHeight;
        if (h2>M) h2=M;
        if (y1>=h2 || y2 <= h1) continue;

        // get the threshold in each slice
        T = slicePercent(Bscope, M, N, h1, h2, percentage);
        if (T<threshGate) T=threshGate;

        // compare and set value
        for (ii=h1; ii<h2; ii++)
        {
            for (j=0; j<N; j++)
            {
                BWBscope[ii*N+j] = (Bscope[ii*N+j]>=T);
            }
        }
    }
}

void blockThreshold(double *Bscope,
                    int M,
                    int N,
                    int *BWBscope,
                    int sliceHeight,
                    int sliceWidth,
                    double percentage,
                    double threshGate,
                    int y1,
                    int y2)
{
    int nVer, nHor;
    double nVerD, nHorD, T, T2;
    int i, j, ii, jj, h1, h2, w1, w2;
    double p2 = LOWBOUND_PCT;

    nVerD = (double)M / (double)sliceHeight;
    nVer  = (int)nVerD;
    if ((double)nVer < nVerD) nVer++;

    nHorD = (double)N / (double)sliceWidth;
    nHor  = (int)nHorD;
    if ((double)nHor < nHorD) nHor++;

    for (i=0; i<nVer; i++)
    {
        for (j=0; j<nHor; j++)
        {
            h1 = i*sliceHeight;
            h2 = (i+1)*sliceHeight;
            w1 = j*sliceWidth;
            w2 = (j+1)*sliceWidth;
            if (h2>M) h2 = M;
            if (w2>N) w2 = N;
            if (y1>=h2 || y2<=h1) continue;

            T = blockPercent(Bscope, M, N, h1, h2, w1, w2, percentage, &T2, &p2);
            if (T<threshGate) T=threshGate;
            if (LOWBOUND_PCT<100 && T2<threshGate) T2 = threshGate;

            for (ii=h1; ii<h2; ii++)
            {
                for (jj=w1; jj<w2; jj++)
                {
                    if (LOWBOUND_PCT>=100)
                        BWBscope[ii*N+jj] = (Bscope[ii*N+jj]>=T);
                    else
                        BWBscope[ii*N+jj] = (int)(Bscope[ii*N+jj]-T2+0.5)*(Bscope[ii*N+jj]>=T);
                }
            }
        }
    }
}


double slicePercent(double *Bscope, // input Bscope image
                    int M,              // number of rows
                    int N,              // number of column
                    int h1,             // start row
                    int h2,             // end row
                    double percentage)
{
    double T;
    double *sliceData;
    double *pt;
    int i, j;
    int N2, M2;

    M2 = (int)(double(h2-h1)/2+0.5);
    N2 = (int)(double(N)/2+0.5);
    sliceData = (double*)calloc(sizeof(double), M2*N2);

    pt = sliceData;
    for (i=h1; i<h2; i+=2)
    {
        for (j=0; j<N; j+=2)
        {
            *pt = Bscope[i*N+j];
            pt++;
        }
    }

    // sort the array
    qsort(sliceData, M2*N2, sizeof(double), dcomp);
    i = (int)((double)M2*(double)N2*percentage/100.0);

    // get percentile
    T = sliceData[i];
    free(sliceData);
    return T;
}

int dcomp(const void *a, const void *b)
{
    if (*(double*)a < *(double*)b)
        return -1;
    else if (*(double*)a > *(double*)b)
        return 1;
    else
        return 0;
}

double blockPercent(double *Bscope,
                    int M,  int N,
                    int h1, int h2,
                    int w1, int w2,
                    double percentage,
                    double *T2,
                    double *p2)
{
    double T;
    double *blockData;
    double *pt;
    int i, j;
    int M2, N2;

    M2 = (int)(double(h2-h1)/2+0.5);
    N2 = (int)(double(w2-w1)/2+0.5);
    blockData = (double*)calloc(sizeof(double), M2*N2);

    pt = blockData;
    for (i=h1; i<h2; i+=2)
    {
        for (j=w1; j<w2; j+=2)
        {
            *pt = Bscope[i*N+j];
            pt++;
        }
    }

    qsort(blockData, M2*N2, sizeof(double), dcomp);
    i = (int)((double)M2*(double)N2*percentage/100.0);

    T = blockData[i];

    if (p2 && T2 && *p2 < 100)
    {
        i = (int)((double)M2*(double)N2*(*p2)/100.0);
        *T2 = blockData[i];
    }

    free(blockData);
    return T;
}

//////////////////////////////////////////////////////////////////////////////
// coordinate transformation
//////////////////////////////////////////////////////////////////////////////
void coorTransform(int *BWBscope,       // input b-n-w Bscope image
                   int M,               // # of rows of BWBscope
                   int N,               // # of columns of BWBscope
                   int *BWcart,         // output b-n-w cartesian image, shall be allocated outside
                   int *pNC,            // # of columns of BWcart
                   int *px0,            // middle point in BWcart
                   double Ang_endpt,    // parameter, angle endpoint, in degrees
                   double *Ang_plot,    // 1 by N vector angle value, in degrees
                   int interFac,        // interpolation factor
                   int y1,              // starting and ending row
                   int y2)
{
    double x0D;
    double t2, t;
    int i, x0, NC, m, n;

    // cartesian image size
    x0D = M*sin(Ang_endpt*UW_PI/180);
    x0 = (int)(x0D)+1;
    NC = 2*x0+1;

    // loop through input Bscope
    for (m=y1; m<y2-1; m++)
    {
        for (n=0; n<N; n++)
        {
            if (BWBscope[m*N+n])
            {
                int val = BWBscope[m*N+n];
                coorTransAssign(BWcart, M, NC, m, Ang_plot[n]*UW_PI/180, x0, LOWBOUND_PCT<100 ? val: 0);

                // angle interpolation
                if (interFac >= 1)
                {
                    // interpolate before
                    if (n>0)
                    {
                        t2 = (Ang_plot[n-1]+Ang_plot[n])/2; // midpoint before
                        for (i=0; i<interFac; i++)
                        {
                            t = ((interFac-i)*t2+i*Ang_plot[n])/(double)(interFac);
                            coorTransAssign(BWcart, M, NC, m, t*UW_PI/180, x0, LOWBOUND_PCT<100 ? val: 0);
                        }
                    }

                    // interpolate after
                    if (n<N-1)
                    {
                        t2 = (Ang_plot[n]+Ang_plot[n+1])/2; // midpoint after
                        for (i=1; i<=interFac; i++)
                        {
                            t = ((interFac-i)*Ang_plot[n]+i*t2)/(double)(interFac);
                            coorTransAssign(BWcart, M, NC, m, t*UW_PI/180, x0, LOWBOUND_PCT<100 ? val: 0);
                        }
                    }
                }
            }
        }
    }
}

void coorTransAssign(int *BWcart,       // b-n-w cartesian image
                     int MC,            // # of rows of BWcart
                     int NC,            // # of columns of BWcart
                     double rho,        // rho in Bscope
                     double theta,      // theta in Bscope
                     int x0,
                     int val)
{
    double xD, yD;
    int x, y;

    xD = (rho+1)*sin(theta)+x0;
    yD = (rho+1)*cos(theta);

    x = (int)(xD-0.5);
    y = (int)(yD-0.5);

    if (y>=0 && y<MC)
    {
        if (x>=0 && x<NC)
        {
            if (val)
                BWcart[y*NC+x] = max(val, BWcart[y*NC+x]);
            else
                BWcart[y*NC+x] = 1;
        }
    }
}



//////////////////////////////////////////////////////////////////////////////
// Hough Transform related funtions
//////////////////////////////////////////////////////////////////////////////
// everything shall be allocated outside this function
void houghTransform(int *BWcart,        // input b-n-w cartesian image
                    int MC,             // # of rows of BWcart
                    int NC,             // # of columns of BWcart
                    int y1,             // starting and ending processing region
                    int y2,
                    double *rho,        // rho
                    int rhoLen,
                    double *theta,      // theta
                    int thetaLen,
                    int *H)         // output image, sizeof (rhoLen rows, thetaLen columns)
{
    int i, thetaIdx, rhoIdx;
    int x, y;
    double *cost = NULL, *sint = NULL;
    double firstRho, slope, oneRho;

    // cos/sin lookup table
    cost = (double *)calloc(thetaLen, sizeof(double));
    sint = (double *)calloc(thetaLen, sizeof(double));
    for (i=0; i<thetaLen; i++)
    {
        cost[i] = cos(theta[i]*UW_PI/180);
        sint[i] = sin(theta[i]*UW_PI/180);
    }

    // rho conversion from value to index
    firstRho = rho[0];
    slope = ((double)rhoLen-1) / (rho[rhoLen-1]-firstRho);


    // Hough Transform
    for (y=y1; y<y2; y++)
    {
        for (x=0; x<NC; x++)
        {
            if (BWcart[y*NC+x])
            {
                for (thetaIdx=0; thetaIdx<thetaLen; thetaIdx++)
                {
                    oneRho = (double)x*cost[thetaIdx] + (double)(y-y1)*sint[thetaIdx];
                    rhoIdx = (int)(slope*(oneRho-firstRho)+0.5);
                    H[rhoIdx*thetaLen+thetaIdx] += BWcart[y*NC+x];
                }
            }
        }
    }

    if (cost) free(cost);
    if (sint) free(sint);
}

// identify peaks after HT
int houghPeaks( int *H,             // Hough data
                int rhoLen,             // # of rows in H
                int thetaLen,           // # of columns in H
                int numPeaks,           // # of peaks to find
                double hGate,           // gate value for a peak
                int nHoodRowHalf,       // neighborhood size half row
                int nHoodColHalf,       // neighborhood size half column
                int *peaks)
{
    int num, done, maxRow, maxCol;
    double maxVal;
    int m, n, nn, m1, m2;
    int M, N;
    M = rhoLen;
    N = thetaLen;

    done = 0;
    num = 0;
    while (!done)
    {
        getMaxVal(H, rhoLen, thetaLen, &maxVal, &maxRow, &maxCol);
        if (maxVal < hGate)
        {
            done = 1;
        }
        else
        {
            if (num==0 && maxVal>2*hGate)
            {
                hGate = maxVal/2;
            }

            peaks[2*num] = maxRow;
            peaks[2*num+1] = maxCol;
            num++;

            if (num == numPeaks)
            {
                done = 1;
                break;
            }
            
            // suppress neighborhood
            for (n=maxCol-nHoodColHalf; n<=maxCol+nHoodColHalf; n++)
            {
                if (n<0)
                    nn = n+N;
                else if (n>=N)
                    nn = n-N;
                else
                    nn = n;

                m1 = maxRow-nHoodRowHalf;
                if (m1<0) m1=0;
                m2 = maxRow+nHoodRowHalf;
                if (m2>M-1) m2=M-1;

                for (m=m1; m<=m2; m++)
                {
                    H[m*thetaLen+nn] = 0;
                }
            }
        }
    }

    return num;
}

void getMaxVal( int *H,
                int M,
                int N,
                double *maxVal,
                int *maxRow,
                int *maxCol)
{
    int i, j;
    *maxVal = -1;
    *maxRow = -1;
    *maxCol = -1;

    /*
    for (i=0; i<M; i++)
    {
        for (j=0; j<N; j++)
        {
            if (*maxVal < H[i*N+j] || *maxVal == H[i*N+j] && *maxCol > j || *maxVal == H[i*N+j] && *maxCol == j && *maxRow > i)
            {
                *maxVal = H[i*N+j];
                *maxRow = i;
                *maxCol = j;
            }
        }
    }
    */
    
    int val;
    
    for (i=0; i<M; i++)
    {
        for (j=0; j<N; j++)
        {
            val = H[i*N+j];
            if (*maxVal < val)
            {
                *maxVal = val;
                *maxRow = i;
                *maxCol = j;
            }
        }
    }
}

double getMaxTheta( int *BWcart,            // input b-n-w cartesian image
                    int MC,             // size of BWcart,
                    int NC,
                    double *theta,          // HT parameter
                    int thetaLen,
                    double *RP_dBm,
                    int *BW,
                    double *Ang_plot,
                    int M0, int N0, int x0,
                    svmData *svm,
                    double *rhos, int *status,
                    double MagMax, double MagMin)
{
    int q, rhoLen, M, N, i, j, maxIdx, *H, *H2;
    double d, *rho, maxSum, *thetaSum;

    M = MC;
    N = NC;
    d = sqrt((double)(M-1)*(M-1)+(double)(N-1)*(N-1));
    q = (int)(d)+1;
    rhoLen = 2*q+1;
    rho = (double*)calloc(rhoLen, sizeof(double));
    H = (int*)calloc(rhoLen*thetaLen, sizeof(int));
    thetaSum = (double*)calloc(thetaLen, sizeof(double));

    memset(H, 0, sizeof(int)*rhoLen*thetaLen);
    memset(thetaSum, 0, sizeof(int)*thetaLen);
    for (i=-q; i<=q; i++)
    {
        rho[i+q] = (double)i;
    }

    //
    houghTransform(BWcart, MC, NC, 0, MC, rho, rhoLen, theta, thetaLen, H);

    H2 = (int*)calloc(rhoLen*thetaLen, sizeof(int));
    memcpy(H2, H, rhoLen*thetaLen*sizeof(int));
    

    // get theta sum
    maxSum = -1;
    maxIdx = -1;
    int localMax = -1;
    int localMaxSum = 0;
    int localIdx = -1;
    for (i=0; i<thetaLen; i++)
    {
        // entropy
        thetaSum[i] = 0;        
        for (j=0; j<rhoLen; j++)
        {
            if (H[j*thetaLen+i])
                thetaSum[i] += H[j*thetaLen+i] * log((double)H[j*thetaLen+i]);
        }

        // local max 
        localMaxSum = 0;
        for (int k = 0; k<3; k++)
        {
            localMax = -1;
            localIdx = -1;
            for (j=0; j<rhoLen; j++)
            {
                if (H[j*thetaLen+i]>localMax) 
                {
                    localMax = H[j*thetaLen+i];
                    localIdx = j;
                }
            }
            localMaxSum += localMax;
            for (j=localIdx-50; j<=localIdx+50; j++)
            {
                if (j>=0 && j<rhoLen)
                    H[j*thetaLen+i] = 0;
            }
        }

        // adjust
        thetaSum[i] *= (double)localMaxSum;
        if (maxSum < thetaSum[i])
        {
            maxSum = thetaSum[i];
            maxIdx = i;
        }
    }
    double theta_to_return = theta[maxIdx];
    if (theta_to_return<0) theta_to_return += 180;

    if (theta_first_frame)
    {       
        theta_previous = theta[maxIdx];
        if (theta_previous<0) theta_previous+=180;
        // theta_particles = theta_initialize_particles(THETA_NUM_PARTICLES, theta_previous, H2, rhoLen, thetaLen, maxIdx);

        tpf = unique_ptr<ThetaParticleFilter>(new ThetaParticleFilter(THETA_NUM_PARTICLES, theta_previous));
    }
    else
    {
        // theta_transition_particles(theta_particles, THETA_NUM_PARTICLES, THETA_VARIANCE);
        // theta_compute_likelihood(theta_particles, THETA_NUM_PARTICLES, H2, rhoLen, thetaLen, &similarity_1);
        // theta_normalize_weights(theta_particles, THETA_NUM_PARTICLES);
        // tparticle *theta_particles_n = tresample(theta_particles, THETA_NUM_PARTICLES);
        
        // theta_free_particles(theta_particles, THETA_NUM_PARTICLES);
        // theta_particles = theta_particles_n;
        // theta_previous = theta_particles[0].theta;
        // theta_previous = (int)(theta_previous + 0.5);

        theta_previous = tpf->track(ThetaObservation(H2, rhoLen, thetaLen));
        theta_to_return = theta_previous;
    }
    
    int nCable = 0;
    double rhoCable[8];
    int isCable[8];

    nCable = rhoTrack(theta_to_return, H2, rho, rhoLen, thetaLen, RP_dBm, BW, BWcart, Ang_plot, MagMax, MagMin, M0, N0, NC, x0, svm, rhoCable, isCable, 8, M>4000?30:15);
    int margin = 8;
    if (M>2000) margin = 15;
    if (M>4000) margin = 30;
    nCable = simpleThetaDetect(theta_to_return, H2, rho, rhoLen, thetaLen, RP_dBm, BW, BWcart, Ang_plot, M0, N0, NC, x0, svm, rhoCable, isCable, 8, margin);
    memcpy(rhos, rhoCable, 8*sizeof(double));
    memcpy(status, isCable, 8*sizeof(int));

    if (g_saveImage)
    {
        showImageIntRho(BWcart, M, NC, M > 4000 ? .25 : .5, M > 4000 ? .25 : .5, 0, 1, COLORMAP_JET, thetaDetectName, theta_to_return, rhoCable, isCable, 8);

        showImageInt(H2, rhoLen, thetaLen, rhoLen > 7000 ? 0.125 : 0.25, 8, 0, 0, COLORMAP_JET, houghName, theta[maxIdx], true);
    }

    if (rho) free(rho);
    if (H)   free(H);
    if (H2)  free(H2);
    if (thetaSum) free(thetaSum);
    
    return theta_to_return;
}

// lines is allocated outside
int detectLines( int *BWcart,           // input b-n-w cartesian image
                 int MC,                // size of BWcart,
                 int NC,
                 int ystart,            // the region to detect
                 int yend,
                 double *theta,         // HT parameter
                 int thetaLen,
                 int sliceHeight,       // slice-processing height
                 int sliceOverlap,      // overlapping between slices
                 int numPeaks,          // # of peaks in each slice
                 double *lines,
                 double ratio)
{
    int sliceStep;
    int sliceNum;
    int offset, n, y1, y2, i;
    int M, N;   // size of one slice
    int q;
    double d;
    int *peaks;
    int totalNum, num;
    int nHoodRowHalf, nHoodColHalf;
    double nHood1, nHood2, hGate;
    int maxR, maxC;

    int rhoLen;
    double *rho;
    int *H;

    getMaxVal(BWcart, MC, NC, &hGate, &maxR, &maxC);
    hGate *= ratio;
    hGate *= HGATE_FACTOR;

    sliceStep = sliceHeight - sliceOverlap;
    sliceNum = (int)((double)(MC-sliceStep)/(double)sliceStep + 1.0);

    M = sliceHeight;
    N = NC;
    d = sqrt((double)(M-1)*(M-1)+(double)(N-1)*(N-1));
    q = (int)(d)+1;
    rhoLen = 2*q+1;
    rho = (double*)calloc(rhoLen, sizeof(double));
    H = (int*)calloc(rhoLen*thetaLen, sizeof(int));

    totalNum = 0;
    peaks = (int*)calloc(numPeaks*2, sizeof(int));
    for (n=0; n<sliceNum; n++)
    {
        offset = n*sliceStep;
        y1 = offset;
        y2 = offset + sliceHeight;
        if (y2>MC) y2=MC;

        // not in the region yet
        if (ystart>=y2 || yend<=y1) continue;

        // HT and related stuff, rho, theta, data allocation
        M = y2-y1;
        N = NC;
        d = sqrt((double)(M-1)*(M-1)+(double)(N-1)*(N-1));
        q = (int)(d)+1;
        rhoLen = 2*q+1;
        memset(H, 0, sizeof(int)*rhoLen*thetaLen);
        for (i=-q; i<=q; i++)
        {
            rho[i+q] = (double)i;
        }

        houghTransform(BWcart, MC, NC, y1, y2, rho, rhoLen, theta, thetaLen, H);

        // peaks identification, hGate and nHoodSize, return line
        nHood1 = (double)rhoLen/100.;
        nHood2 = (double)thetaLen/30.;
        nHoodRowHalf = (int)(nHood1);
        nHoodColHalf = max((int)(nHood2), 1);
        //hGate = (double)nHoodRowHalf*10;

        num = houghPeaks(H, rhoLen, thetaLen, numPeaks, hGate, nHoodRowHalf, nHoodColHalf, peaks);

        // from peaks location to rho and theta value
        for (i=0; i<num; i++)
        {
            // rho
            lines[2*totalNum] = rho[ peaks[2*i] ];
            lines[2*totalNum] += (double)(offset)*sin(theta[ peaks[2*i+1] ]*UW_PI/180);

            // theta
            lines[2*totalNum+1] = theta[ peaks[2*i+1] ];            
            totalNum++;
        }
    }

    if (peaks) free(peaks);
    if (rho) free(rho);
    if (H) free(H);

    return totalNum;
}


// retrieve line data
// everything is NOT allocated within
int getLineData(double *Bscope,
                double *Ang_plot,
                int M,                  // # of rows in Bscope
                int N,                  // # of columns; length of Ang_plot
                int *BWBscope,
                int x0,
                double rho,
                double theta,           // theta in degree format
                double *data,           // retrieved data
                int *range,             // range index
                int *angle,             // angle index
                bool makeLong)
{
    int num;    // number of points on that line
    int T1, T2; // first and last non-zero point theta index
    int Ridx, Tidx;
    int i;
    double R, x0cost;
    double *sint;

    if (theta>=90) theta-= 180;

    x0cost = (double)x0 * cos(theta*UW_PI/180);
    sint = (double*)calloc(N, sizeof(double));
    for (i=0; i<N; i++)
    {
        sint[i] = sin((Ang_plot[i]+theta)*UW_PI/180);
    }

    if (BWBscope)
    {
        // first non-zero point
        for (T1=1; T1<N-1; T1++)
        {
            R = (rho-x0cost)/sint[T1];
            Ridx = (int)(R+0.5);
            if (Ridx<1 || Ridx>M) continue;
            if (BWBscope[(Ridx-1)*N+T1]) break;
        }
        
        // last non-zero point
        for (T2=N-2; T2>=1; T2--)
        {
            R = (rho-x0cost)/sint[T2];
            Ridx = (int)(R+0.5);
            if (Ridx<1 || Ridx>M) continue;
            if (BWBscope[(Ridx-1)*N+T2]) break;
        }
        
        // swap if necessary
        if (T1>T2)
        {
            Tidx=T1;
            T1=T2;
            T2=Tidx;
        }

        // make long
        if (makeLong)
        {
            if (T1 > CABLE_TIP_END)
                T1 = CABLE_TIP_END;
            if (N > 300)
            {
                if (T2 < 350-CABLE_TIP_END-1)
                    T2 = 350-CABLE_TIP_END-1;
            }
            else
            {
                if (T2 < 176-CABLE_TIP_END-1)
                    T2 = 176-CABLE_TIP_END-1;
            }
        }
    }
    else
    {
        T1 = CABLE_TIP_END;
        if (N > 300)
            T2 = 350-CABLE_TIP_END-1;
        else
            T2 = 176-CABLE_TIP_END-1;
    }

    // retrieve data
    num = 0;
    for (Tidx=T1; Tidx<=T2; Tidx++)
    {
        R = (rho-x0cost)/sint[Tidx];
        Ridx = (int)(R-0.5);
        if (Ridx<0 || Ridx>=M) continue;

        range[num] = Ridx;
        angle[num] = Tidx;
        data[num]  = Bscope[Ridx*N+Tidx];
        num++;
    }

    if (sint) free(sint);
    return num;
}


//////////////////////////////////////////////////////////////////////////////
// feature-computation related function
//////////////////////////////////////////////////////////////////////////////

// everything is allocated outside
void powerSpectrum(double *input,       // input data
                   int size,            // data size (which is also DFT size)
                   double *spec)
{
    fftw_complex *data, *fftRes;
    fftw_plan planFwd;
    int i;

    data    = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*size);
    fftRes  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*size);

    planFwd = fftw_plan_dft_1d(size, data, fftRes, FFTW_FORWARD, FFTW_ESTIMATE);

    // populate input data
    for (i=0; i<size; i++)
    {
        data[i][0] = input[i];
        data[i][1] = 0.0;
    }

    fftw_execute(planFwd);

    // get spectrum
    for (i=0; i<size; i++)
    {
        spec[i] = fftRes[i][0]*fftRes[i][0] + fftRes[i][1]*fftRes[i][1];
    }

    // free memory
    fftw_destroy_plan(planFwd);

    fftw_free(data);
    fftw_free(fftRes);
}


void autoCorrelation(double *input,
                     int size,
                     double *autoCo)
{
    fftw_complex *data, *fftRes, *ifftRes;
    fftw_plan planFwd, planBwd;
    int i;
    double a1;

    data    = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*size);
    fftRes  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*size);
    ifftRes = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*size);
    
    planFwd = fftw_plan_dft_1d(size, data, fftRes,  FFTW_FORWARD,  FFTW_ESTIMATE);
    planBwd = fftw_plan_dft_1d(size, data, ifftRes, FFTW_BACKWARD, FFTW_ESTIMATE);

    // input data
    for (i=0; i<size; i++)
    {
        data[i][0] = input[i];
        data[i][1] = 0.0;
    }

    fftw_execute(planFwd);

    // compute abs square
    for (i=0; i<size; i++)
    {
        data[i][0] = fftRes[i][0]*fftRes[i][0] + fftRes[i][1]*fftRes[i][1];
        data[i][1] = 0.0;
    }

    fftw_execute(planBwd);

    // auto correlation
    for (i=0; i<size; i++)
    {
        autoCo[i] = ifftRes[i][0];
    }

    a1 = autoCo[0];
    if (a1!=0.0)
    {
        for (i=0; i<size; i++)
        {
            autoCo[i] = autoCo[i] / a1;
        }
    }
    
    fftw_destroy_plan(planFwd);
    fftw_destroy_plan(planBwd);

    fftw_free(data);
    fftw_free(fftRes);
    fftw_free(ifftRes);
}


void computeFeature(double *data,       // input data
                    int dataSize,       // data size
                    int specSize,       // DFT size for computing power spectrum
                    int autoSize,       // DFT size for computing autocorrelation
                    double *feature,    // output feature
                    int featSize)
{
    if (dataSize == 0)
    {
        // Why would dataSize be 0?
        for (int i = 0; i < featSize; i++)
        {
            feature[i] = 0.0;
        }
        return;
    }

    int i;

    double *specData, *autoData, *sortData;
    double *spec, *autoCo;
    double iD;

    specData= (double*)calloc(specSize, sizeof(double));
    autoData= (double*)calloc(autoSize, sizeof(double));
    spec    = (double*)calloc(specSize, sizeof(double));
    autoCo  = (double*)calloc(autoSize, sizeof(double));

    memcpy(specData, data, min(dataSize, specSize)*sizeof(double));
    memcpy(autoData, data, min(dataSize, autoSize)*sizeof(double));

    // feature 1&2: power spectrum distribution
    powerSpectrum(specData, specSize, spec);
    if (spec[0]!=0)
    {
        feature[0] = sumArray(spec, 1, 5)/spec[0]*100.0;
        feature[1] = sumArray(spec, 5, specSize)/spec[0]*100.0;
    }

    // feature 3~6: statistics
    sortData = (double*)calloc(dataSize, sizeof(double));
    memcpy(sortData, data, dataSize*sizeof(double));

    qsort(sortData, dataSize, sizeof(double), dcomp);

    feature[2] = sumArray(data, 0, dataSize)/(double)dataSize;
    feature[3] = sortData[dataSize-1];
    iD = (double)dataSize*0.95;
    i = (int)iD;
    if ((double)i < iD) i++;
    feature[4] = sortData[i-1];
    iD = (double)dataSize*0.68;
    i = (int)iD;
    if ((double)i < iD) i++;
    feature[5] = sortData[i-1];

    // feature 7~9: threshold with 95%
    thresholdFeature(data, dataSize, feature[4], &feature[6], &feature[7], &feature[8]);
    
    // feature 10~12: threshold with 68%
    thresholdFeature(data, dataSize, feature[5], &feature[9], &feature[10], &feature[11]);

    // feature 13~14: autocorrelation
    for (i=0; i<min(dataSize,autoSize); i++)
    {
        autoData[i] += 72.0;
    }
    autoCorrelation(autoData, autoSize, autoCo);
    feature[13] = sumArray(autoCo, 0, autoSize)/(double)autoSize;   // mean
    qsort(autoCo, autoSize, sizeof(double), dcomp);
    
    // median
    if (autoSize % 2)       // exact median
        feature[12] = autoCo[(int)(autoSize/2)];
    else
        feature[12] = (autoCo[autoSize/2] + autoCo[autoSize/2-1])/2.0;


    // clean up
    if (specData) free(specData);
    if (autoData) free(autoData);
    if (spec) free(spec);
    if (autoCo) free (autoCo);
    if (sortData) free(sortData);
}


// sum array from data[n1] to data[n2-1]
double sumArray(double *data,
                int n1,
                int n2)
{
    double s = 0.0;
    int i;
    for (i=n1; i<n2; i++)
    {
        s += data[i];
    }

    return s;
}


// threshold feature group
void thresholdFeature(double *data,     // data and size
                      int size,
                      double T,         // threshold
                      double *f1,       // three features
                      double *f2,
                      double *f3)
{
    int numSeg, locLen, i, distLen;
    int *loc;
    double *dist, *center;

    loc = (int*)calloc(size+1, sizeof(int));
    locLen = compareThreshold(data, size, T, &numSeg, loc);

    for (i=0; i<locLen-1; i++)
    {
        loc[i]=loc[i+1];
    }
    locLen--;

    if (locLen%2)
    {
        loc[locLen]=size;
        locLen++;
    }

    center = (double*)calloc(locLen/2, sizeof(double));
    for (i=0; i<locLen/2; i++)
    {
        center[i] = loc[2*i+1]+loc[2*i];
        center[i] = center[i]/2.0;
    }

    if (locLen/2 <= 1)
    {
        dist = (double*)calloc(1, sizeof(double));
        dist[0] = 0;
        distLen = 1;
    }
    else
    {
        dist = (double*)calloc(locLen/2-1, sizeof(double));
        distLen = locLen/2-1;
        for (i=0; i<distLen; i++)
        {
            dist[i] = center[i+1]-center[i];
        }
    }


    *f1 = (double)numSeg/(double)size*100.0;
    *f2 = sumArray(dist, 0, distLen) / (double)distLen;
    *f3 = uwStd(dist, distLen);

    if (loc) free(loc);
    if (center) free(center);
    if (dist) free(dist);
}

int compareThreshold(double *data,
                     int dataLen,
                     double T,
                     int *num,
                     int *loc)
{
    bool newSeg = 0;
    int locLen = 0;
    int i;
    loc[locLen++] = 0;
    num[0] = 0;

    for (i=0; i<dataLen; i++)
    {
        if (data[i]>=T)
        {
            if (!newSeg)
            {
                newSeg=1;
                num[0]++;
                loc[locLen++]=i+1;
            }
        }
        else
        {
            if (newSeg)
            {
                newSeg=0;
                loc[locLen++]=i+1;
            }
        }
    }

    return locLen;
}

double uwStd(double *data, int dataLen)
{
    int i;
    double m, std;
    std = 0;
    if (dataLen <=1 )
    {
        return 0;
    }
    else
    {
        m = sumArray(data, 0, dataLen) / (double)dataLen;
        for (i=0; i<dataLen; i++)
        {
            std += (data[i]-m)*(data[i]-m);
        }
        std /= (double)(dataLen-1);
        std = sqrt(std);
        return std;
    }
}

void printArray(double *data, int size)
{
    int i;
    for (i=0; i<size; i++)
    {
        printf("%f\n", data[i]);
    }
    printf("\n\n\n");
}

void svmRead(const char* filename, svmData *svm)
{
    FILE *fp;
    int i,j, M, N;
    double *buf = NULL;

    fp = fopen(filename, "rb");
    if (!fp)
    {
        printf("SVM read error");
        exit(1);
    }

    fread(&(svm->M), sizeof(int), 1, fp);
    fread(&(svm->N), sizeof(int), 1, fp);
    M = svm->M;
    N = svm->N;

    // allocate space
    svm->sv = (double**)calloc(sizeof(double*), M);
    for (i=0; i<M; i++)
        (svm->sv)[i] = (double*)calloc(sizeof(double), N);
    svm->alpha = (double*)calloc(sizeof(double), M);
    svm->shift = (double*)calloc(sizeof(double), N);
    svm->scale = (double*)calloc(sizeof(double), N);
    svm->b =(double*)calloc(sizeof(double), 2*N+1);

    // read sv
    buf = (double*)calloc(sizeof(double), M*N);
    fread(buf, sizeof(double), M*N, fp);
    for (i=0; i<M; i++)
    {
        for (j=0; j<N; j++)
        {
            svm->sv[i][j] = buf[j*M+i];
        }
    }

    // alpha, bias, shift, scale, and b
    fread(svm->alpha, sizeof(double), M, fp);
    fread(&(svm->bias), sizeof(double), 1, fp);
    fread(svm->shift, sizeof(double), N, fp);
    fread(svm->scale, sizeof(double), N, fp);
    fread(svm->b, sizeof(double), 2*N+1, fp);

    // get sum of squares of sv
    svm->svss = (double*)calloc(sizeof(double), M);
    for (i=0; i<M; i++)
    {
        svm->svss[i] = 0.0;
        for (j=0; j<N; j++)
        {
            svm->svss[i] += (svm->sv[i][j])*(svm->sv[i][j]);
        }
    }

    if (buf) free(buf);
}

void svmFree(svmData *svm)
{
    int i;
    if (svm->sv)
    {
        for (i=0; i<svm->M; i++)
        {
            free((svm->sv)[i]);
        }
    }

    if (svm->alpha) free(svm->alpha);
    if (svm->shift) free(svm->shift);
    if (svm->scale) free(svm->scale);
    if (svm->svss)  free(svm->svss);
    if (svm->b) free(svm->b);
}

double svmCalc(svmData *svm, double *f)
{
    double sigma, a, val;
    double *s;
    int M = svm->M;
    int N = svm->N;
    sigma = 1.0;
    
    // change the feature
    int i, j;
    for (i=0; i<svm->N; i++)
    {
        f[i] += svm->shift[i];
        f[i] *= svm->scale[i];
    }

    // 
    s = (double*)calloc(sizeof(double), svm->M);
    for (i=0; i<M; i++)
    {
        s[i] = 0;
        for (j=0; j<N; j++)
        {
            s[i] += svm->sv[i][j]*f[j];
        }
        s[i] *= -2.0;
    }
    a = 0;
    for (j=0; j<N; j++)
    {
        a += f[j]*f[j];
    }
    for (i=0; i<M; i++)
    {
        s[i] += svm->svss[i];
        s[i] += a;
        s[i] = exp(-1/(2*sigma*sigma)*s[i]);
    }

    //
    val = 0;
    for (i=0; i<M; i++)
    {
        val += svm->alpha[i] * s[i];
    }

    val += svm->bias;

    if (s) free(s);

    return val;
}

double linearCalc(svmData *svm, double *f)
{
    double val;
    int i;
    val = svm->b[0];

    for (i=0; i<svm->N; i++)
    {
        val += f[i]*svm->b[i+1];
    }

    for (i=0; i<svm->N; i++)
    {
        val += f[i]*f[i]*svm->b[i+1+svm->N];
    }

    return val;
}

int sumArrayInt(int *data, int n1, int n2)
{
    int sum = 0;
    int i;
    for (i=n1; i<n2; i++)
    {
        sum += data[i];
    }

    return sum;
}

// *cLs shall be NULL pointer when calling this function
int singleFrameDetect(double *RP_dBm,
                      int M,            // # of rows of RP_dBm
                      int N,            // # of columns of RP_dBm
                      double *Ang_plot, // angle values
                      double Ang_endpt, // angle value end
                      double MagMax, double MagMin,
                      svmData *svm,
                      int y1,
                      int y2,
                      cableLine **cLs,
                      int *px0,
                      cableLine *tCables,
                      FILE *pFile)
{
    int nCable = 0;
    int nLine, dataNum;
    int i;
    int interFactor, sliceLineNum, sliceSize, sliceOverlap;
    double svmVal, linearVal;
    cableLine *tmp = NULL;
    double *lines = NULL;
    double ratio;

    int NC, x0;
    int *BW = NULL;
    int *BWcart = NULL;

    double lineData[512];
    int lineRange[512];
    int lineAngle[512];
    double f[FEAT_SIZE];

    clock_t t1, t2, t3, t4, t5;
    double tm1, tm2, tm3, tm4;

    if (pFile)
        t1=clock();

    tmp = (cableLine*)calloc(128, sizeof(cableLine));   
    BW = (int*)calloc(M*N, sizeof(int));
    if (M>4000) HGATE_FACTOR = 24;
    sliceThreshold(RP_dBm, M, N, BW, SLICE_HEIGHT, TH_PCT, TH_GATE, y1, y2, SLICE_WIDTH);

    double hFa;
    if (M>4000)
        hFa = 0.25;
    else if (M>2000)
        hFa = 0.5;
    else
        hFa = 1;


    if (g_saveImage)
    {
    #ifdef _DEBUG
        showImageBW(BW, M, N, hFa, N>300 ? 1 : 2, Ang_plot[0]<0);
        showImageInt(BW, M, N, hFa, N>300 ? 1 : 2, Ang_plot[0]<0, 0);
    #endif
        showImageInt(BW, M, N, hFa, N>300 ? 1 : 2, Ang_plot[0]<0, 0, COLORMAP_JET, threshName);
    }

    if (pFile)
        t2 = clock();


    // 2. coordinate transform
    interFactor = (int)((double)M/512.0+0.5);
    // cartesian image size
    double x0D = M*sin(Ang_endpt*UW_PI/180);
    x0 = (int)(x0D)+1;
    NC = 2*x0+1;
    BWcart = (int*)calloc(M*NC, sizeof(int));

    coorTransform(BW, M, N, BWcart, &NC, &x0, Ang_endpt, Ang_plot,
        interFactor, y1, y2);       // check.

    if (pFile)
        t3 = clock();


    // 3. line detection
    sliceSize = (int)((double)M/24.0);
    sliceOverlap = (int)((double)sliceSize/4.0);
    sliceLineNum = 3*(int)(sliceSize/(double)(M/16.0)+0.5);
    int sliceStep = sliceSize - sliceOverlap;
    int sliceNum = (int)((double)(M-sliceStep)/(double)sliceStep + 1.0);
    lines = (double*)malloc(sizeof(double)*2*sliceNum*sliceLineNum);

    ratio = 1.0;

    nLine = detectLines(BWcart, M, NC, y1, y2, THETA, THETA_LEN,
        sliceSize, sliceOverlap, sliceLineNum, lines, ratio);

    double rhos[8];
    int status[8];
    double maxTheta = getMaxTheta(BWcart, M, NC, PTHETA, PTHETA_LEN, RP_dBm, BW, Ang_plot, M, N, x0, svm, rhos, status, MagMax, MagMin);
    // right here get rho score and show 
    rhoScore = doubleFrameScore(rhoCableLine, rhoNum, rhoCablePrev, rhoNumPrev, rhoScorePrev, FRAME_I==TOTAL_START);
    fprintf(pRhoFile, "%f\n", rhoScore);

    int num_trackers2 = get_rtracker_num(g_rtrackers_);

    // This seems to be the final result image
    showRtrackerImage(RP_dBm, BW, Ang_plot, M, N, x0, MagMax, MagMin, g_rtrackers_, num_trackers2, rhoDetectName, rhoScore>=.5);

    fprintf(pTheta, "%4f\n", maxTheta);
    X0 = x0;

    if (g_saveImage)
    {
        //showImageInt(BWcart, M, NC, hFa, hFa, 0, 1, COLORMAP_JET, transName, maxTheta); // show direction
        showImageInt(BWcart, M, NC, hFa, hFa, 0, 1, COLORMAP_JET, transName, -100);
    }

    nLine = removeRepLines(lines, nLine);

    if (pFile)
        t4 = clock();

    // theta cables
    for (i=0; i<8; i++)
    {
        dataNum = getLineData(RP_dBm, Ang_plot, M, N,
            BW, x0, rhos[i], maxTheta,
            lineData, lineRange, lineAngle);
        tCables[i].isCable = status[i]>0;
        tCables[i].rho = rhos[i];
        tCables[i].theta = maxTheta;
        tCables[i].nPix = dataNum;
        tCables[i].score = tCables[i].isCable ? 1 : 0;
        if (tCables[i].data) free(tCables[i].data);
        if (tCables[i].range) free(tCables[i].range);
        if (tCables[i].angle) free(tCables[i].angle);
        tCables[i].data = (double*)calloc(dataNum, sizeof(double));
        tCables[i].range = (int*)calloc(dataNum, sizeof(int));
        tCables[i].angle = (int*)calloc(dataNum, sizeof(int));
        memcpy(tCables[i].data, lineData, dataNum*sizeof(double));
        memcpy(tCables[i].range, lineRange, dataNum*sizeof(int));
        memcpy(tCables[i].angle, lineAngle, dataNum*sizeof(int));
    }
    // 4. get parameters and classification results for each line
    for (i=0; i<nLine; i++)
    {
        dataNum = getLineData(RP_dBm, Ang_plot, M, N,
            BW, x0, lines[2*i], lines[2*i+1],
            lineData, lineRange, lineAngle);
        computeFeature(lineData, dataNum, SPEC_SIZE, AUTO_SIZE, f, FEAT_SIZE);
        // svm
        linearVal = linearCalc(svm, f);
        svmVal = svmCalc(svm, f);

        if (dataNum>=15 && svmVal<0 && linearVal>=0.5)
        {
            tmp[nCable].isCable = true;
            tmp[nCable].score = linearVal;
            tmp[nCable].rho = lines[2*i];
            tmp[nCable].theta = lines[2*i+1];
            tmp[nCable].nPix = dataNum;
            tmp[nCable].data = (double*)calloc(dataNum, sizeof(double));
            tmp[nCable].range = (int*)calloc(dataNum, sizeof(int));
            tmp[nCable].angle = (int*)calloc(dataNum, sizeof(int));
            memcpy(tmp[nCable].data, lineData, dataNum*sizeof(double));
            memcpy(tmp[nCable].range, lineRange, dataNum*sizeof(int));
            memcpy(tmp[nCable].angle, lineAngle, dataNum*sizeof(int));
            nCable++;
        }
    }
    
    // copy back
    if (!nCable)
        *cLs = NULL;
    else
        *cLs = (cableLine*)calloc(nCable, sizeof(cableLine));

    for (i=0; i<nCable; i++)
    {
        (*cLs)[i] = tmp[i];
    }

    // clean up
    if (tmp) free(tmp);
    if (BW) free(BW);
    if (BWcart) free(BWcart);
    if (lines) free(lines);

    if (pFile)
    {
        t5=clock();
        tm1 = double(t2-t1)/(double)CLOCKS_PER_SEC;
        tm2 = double(t3-t2)/(double)CLOCKS_PER_SEC;
        tm3 = double(t4-t3)/(double)CLOCKS_PER_SEC;
        tm4 = double(t5-t4)/(double)CLOCKS_PER_SEC;
        fprintf(pFile, "%f\t%f\t%f\t%f\t", tm1, tm2, tm3, tm4);
    }

    *px0 = x0;
    return nCable;
}

void freeLineArray(cableLine *&cLs, int num)
{
    if (cLs)
    {   
        for (int i=0; i<num; i++)
        {
            if (cLs[i].data) free(cLs[i].data);
            if (cLs[i].angle) free(cLs[i].angle);
            if (cLs[i].range) free(cLs[i].range);
            cLs[i].data = NULL;
            cLs[i].angle = NULL;
            cLs[i].range = NULL;
        }

        free(cLs);
        cLs = NULL;
    }
}

// calculate the frame score from detection results of two frames
double doubleFrameScore(cableLine *current, int currentNum,
                        cableLine *prev, int prevNum, double prevScore,
                        bool firstFrame)    // to indicate whether this is the first frame
{
    static double AN = 0.0;
    static double A = 0.5;
    static double BOOST = 1.5;
    static double DECLINE = 1.0/BOOST;
    static double THETA_BOOST[5] = {1.0, 1.2, 1.3, 1.4, 1.5};
    static double CORRE_BOOST[5] = {1.2, 1.3, 1.4, 1.5, 1.5};
    double s;
    int nParallel;
    int nCorrelated;

    // intital score based on # of cables
    s = AN;
    if (numCable(current, currentNum))
    {
        s = A;
        // parallel check
        nParallel = checkParallel(current, currentNum);
        if (nParallel >= 5)
            s = s*BOOST;
        else
            s = s*THETA_BOOST[nParallel];
    }

    // factor previous score
    if (!firstFrame)
    {
        s = (s+prevScore)/2.0;
    }

    // target correlation
    if (!firstFrame)
    {
        nCorrelated = temporalCorrelate(current, currentNum, prev, prevNum);
        if (nCorrelated == 0)
            s = s*DECLINE;
        else if (nCorrelated >= 5)
            s = s*BOOST;
        else
            s = s*CORRE_BOOST[nCorrelated-1];
    }

    // clip
    if (s>1) s=1.0;

    return s;
}

int numCable(cableLine *cl, int num)
{
    int n(0), i;
    for (i=0; i<num; i++)
    {
        if (cl[i].isCable && cl[i].score>=0.5)
        {
            n++;
        }
    }

    return n;
}

int checkParallel(cableLine *cl, int num)
{
    int i, j;
    int *nArray = (int*)calloc(num, sizeof(int));
    for (i=0; i<num; i++)
    {
        if (cl[i].isCable && cl[i].score>=.5)
        for (j=0; j<num; j++)
        {
            if (!cl[j].isCable) continue;
            if (cl[j].score<.5) continue;
            if (abs(cl[j].theta-cl[i].theta)<=.5)
                nArray[i]++;
        }
    }

    j = 0;
    for (i=0; i<num; i++)
    {
        if (j<nArray[i])
            j = nArray[i];
    }

    free(nArray);
    return (j-1);
}

int temporalCorrelate(cableLine *current, int currentNum,
                     cableLine *prev, int prevNum)
{
    static double RSIGMA = 150.0;
    static double TSIGMA = 3.0;
    static double C = exp(1.0);
    int n(0), i, j;
    double *corre = NULL;

    if (!prev || !prevNum)
    {
        return 0;
    }

    corre = (double*)calloc(prevNum, sizeof(double));

    for (i=0; i<currentNum; i++)
    {
        if (current[i].isCable && current[i].score>=.5)
        {
            memset(corre, 0, sizeof(double)*prevNum);
            for (j=0; j<prevNum; j++)
            {
                corre[j] = lineCorrelate(current[i], prev[j], RSIGMA, TSIGMA, C);
            }

            for (j=1; j<prevNum; j++)
            {
                if (corre[j]>corre[0])
                    corre[0]=corre[j];
            }
            if (corre[0]>=1) n++;
            current[i].score = current[i].score*corre[0];
            current[i].isCable = current[i].isCable && (current[i].score>=.5);
        }
    }

    free(corre);
    return n;
}

double lineCorrelate(cableLine l1, cableLine l2, 
                     double Rsigma, double Tsigma, double C)
{
    double cr, ca, score;
    cr = exp(-(abs(l1.rho)-abs(l2.rho))*(abs(l1.rho)-abs(l2.rho))/(2*Rsigma*Rsigma));
    ca = exp(-(abs(l1.theta)-abs(l2.theta))*(abs(l1.theta)-abs(l2.theta))/(2*Tsigma*Tsigma));
    score = cr*ca*C;

    return score;
}

int removeRepLines(double *cLs, int num)
{
    int *idx = (int*)calloc(num, sizeof(int));
    double *cL2 = (double*)calloc(num, sizeof(double)*2);
    int i, j;
    int newNum;

    for (i=0; i<num; i++)
    {
        idx[i] = 1;
        cL2[2*i] = cLs[2*i];
        cL2[2*i+1] = cLs[2*i+1];
    }

    for (i=num-1; i>=0; i--)
    {
        for (j=i-1; j>=0; j--)
        {
            if (abs(cLs[2*j+1]-cLs[2*i+1])<=1.0) // same line, theta
            {
                if (abs(cLs[2*j]-cLs[2*i])<=5.0)                // rho
                {
                    idx[i] = 0;
                }
            }
        }
    }

    newNum=0;
    for (i=0; i<num; i++)
    {
        if (idx[i])
        {
            cLs[2*newNum] = cL2[2*i];
            cLs[2*newNum+1] = cL2[2*i+1];
            newNum++;
        }
    }

    if (idx) free(idx);
    if (cL2) free(cL2);

    return newNum;
}

void showImageInt(int *data, int height, int width, double hFactor, double wFactor, int flipLR, int scale, unsigned char COLORMAP[64][3], const char *nm, double theta, bool bHough)
{
    int nCh = 3;
    IplImage *img = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, nCh);
    uchar *imgData = (uchar*)img->imageData;
    int wStep = img->widthStep;
    int i, j;
    int idx;
    double scaleD;
    int *dataC = data;

    if (!scale)
    {
        int maxData = -100;
        int *data2 = data;
        for (i=0; i<height; i++)
        {
            for (j=0; j<width; j++)
            {
                if (*data2>maxData)
                    maxData = *data2;
                data2++;
            }
        }

        scaleD = 63.0/double(maxData);
    }

    for (i=0; i<height; i++)
    {
        for (j=0; j<width; j++)
        {
            if (scale)
                idx = *data*scale;
            else
                idx = int(*data*scaleD);
            data++;

            if (idx<0) idx=0;
            if (idx>63) idx=63;

            imgData[j*nCh+0] = COLORMAP[idx][2];
            imgData[j*nCh+1] = COLORMAP[idx][1];
            imgData[j*nCh+2] = COLORMAP[idx][0];
        }
        imgData += wStep;
    }

    // plot lines with theta
    if (!bHough && theta > -99)
    {
        if (theta<0) theta = 180+theta;
        imgData = (uchar*)img->imageData;
        double yD;
        int x, y;
        for (double rho = 500; rho <= 3000; rho += 500)
        {
            for (x = 0; x < width; x++)
            {
                yD = (rho - (double)x * cos(theta*UW_PI/180.)) / sin(theta*UW_PI/180.);
                //if (theta<0) yD = -yD;
                y = (int)(yD+0.5);
                for (i=-3; i<=3; i++)
                {
                    if (y+i>=0 && y+i<height)
                    {
                        imgData[(y+i)*wStep+x*nCh+0] = 0;
                        imgData[(y+i)*wStep+x*nCh+1] = 0;
                        imgData[(y+i)*wStep+x*nCh+2] = 255;
                    }
                }
            }
        }
    }

    if (bHough && theta>-99)    // show the maximum entropy line & maximum position line
    {
        int x, y;
        int tIdx = 0;
        if (theta<0) tIdx = (int)(theta+0.5+90);
        else tIdx = (int)(theta+0.5-(90-PTHETA_LEN));
        imgData = (uchar*)img->imageData;
        for (y=0; y<500; y++)
        {
            for (x=tIdx; x<=tIdx; x++)
            {               
                imgData[y*wStep+x*nCh+0] = 0;
                imgData[y*wStep+x*nCh+1] = 255;
                imgData[y*wStep+x*nCh+2] = 255;
            }               
        }

        double maxVal;
        getMaxVal(dataC, height, width, &maxVal, &y, &x);
        for (y=0; y<500; y++)
        {           
            {               
                imgData[y*wStep+x*nCh+0] = 0;
                imgData[y*wStep+x*nCh+1] = 0;
                imgData[y*wStep+x*nCh+2] = 255;
            }               
        }
    }


    // shift the Hough Transform data
    if (bHough)
    {
        int xsl = (width-1)/2;
        int l;
        IplImage *imgC = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, nCh);
        cvCopy(img, imgC);
        imgData = (uchar*)img->imageData;
        uchar *imgDataC = (uchar*)imgC->imageData;
        for (j=0; j<height; j++)
        {
            for (i=0; i<=xsl; i++)
            {
                for (l=0; l<3; l++)
                    imgData[(i+xsl)*nCh+l] = imgDataC[i*nCh+l];
            }
            for (i=xsl+1; i<=2*xsl; i++)
            {
                for (l=0; l<3; l++)
                    imgData[(i-xsl-1)*nCh+l] = imgDataC[i*nCh+l];
            }
            imgData += wStep;
            imgDataC += wStep;
        }
        cvReleaseImage(&imgC);
    }

    IplImage *img2 = cvCreateImage(cvSize(int(width*wFactor), int(height*hFactor)), IPL_DEPTH_8U, nCh);
    cvResize(img,img2);

    if (!bHough) cvFlip(img2);
    if (flipLR) cvFlip(img2, NULL, 1);

    if (!nm)
    {
        cvNamedWindow("test window");
        cvShowImage("test window", img2);
        cvWaitKey();
        cvReleaseImage(&img);
        cvReleaseImage(&img2);
        cvDestroyWindow("test window");
    }
    else
    {
        cvSaveImage(nm, img2);
        cvReleaseImage(&img);
        cvReleaseImage(&img2);
    }
}


void showImageIntRho(int *data, int height, int width, double hFactor, double wFactor, int flipLR, int scale, unsigned char COLORMAP[64][3], const char *nm, double theta, double rhos[], int isCable[], int rhoNum)
{
    int nCh = 3;
    IplImage *img = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, nCh);
    uchar *imgData = (uchar*)img->imageData;
    int wStep = img->widthStep;
    int i, j;
    int idx;
    double scaleD;

    if (!scale)
    {
        int maxData = -100;
        int *data2 = data;
        for (i=0; i<height; i++)
        {
            for (j=0; j<width; j++)
            {
                if (*data2>maxData)
                    maxData = *data2;
                data2++;
            }
        }

        scaleD = 63.0/double(maxData);
    }

    for (i=0; i<height; i++)
    {
        for (j=0; j<width; j++)
        {
            if (scale)
                idx = *data*scale;
            else
                idx = int(*data*scaleD);
            data++;

            if (idx<0) idx=0;
            if (idx>63) idx=63;

            imgData[j*nCh+0] = COLORMAP[idx][2];
            imgData[j*nCh+1] = COLORMAP[idx][1];
            imgData[j*nCh+2] = COLORMAP[idx][0];
        }
        imgData += wStep;
    }

    // plot lines with theta
    if (theta > -99)
    {
        if (theta>=90) theta = theta-180;
        imgData = (uchar*)img->imageData;
        double yD;
        int x, y;
        for (int k=0; k<rhoNum; k++)
        {
            double rho = rhos[k];
            for (x = 0; x < width; x++)
            {
                yD = (rho - (double)x * cos(theta*UW_PI/180.)) / sin(theta*UW_PI/180.);
                //if (theta<0) yD = -yD;
                y = (int)(yD+0.5);
                for (i=-3; i<=3; i++)
                {
                    if (y+i>=0 && y+i<height)
                    {
                        /*
                        if (isCable[k]==2)  // dormant
                        {
                            imgData[(y+i)*wStep+x*nCh+0] = 0;
                            imgData[(y+i)*wStep+x*nCh+1] = 255;
                            imgData[(y+i)*wStep+x*nCh+2] = 255;
                        }
                        */
                        if (isCable[k]<2)   // candidate + healthy
                        {
                            imgData[(y+i)*wStep+x*nCh+0] = 0;
                            imgData[(y+i)*wStep+x*nCh+1] = isCable[k]?0:255;
                            imgData[(y+i)*wStep+x*nCh+2] = isCable[k]?255:0;
                        }                       
                    }
                }
            }
        }
    }

    IplImage *img2 = cvCreateImage(cvSize(int(width*wFactor), int(height*hFactor)), IPL_DEPTH_8U, nCh);
    cvResize(img,img2);

    if (flipLR) cvFlip(img2, NULL, 1);

    if (!nm)
    {
        cvNamedWindow("test window");
        cvShowImage("test window", img2);
        cvWaitKey();
        cvReleaseImage(&img);
        cvReleaseImage(&img2);
        cvDestroyWindow("test window");
    }
    else
    {
        cvSaveImage(nm, img2);
        cvReleaseImage(&img);
        cvReleaseImage(&img2);
    }
}



void voteTheta(cableLine *cables, int cableNum)
{
    int *countTheta = NULL;
    countTheta = (int*)calloc(THETA_LEN, sizeof(int));
    int i;
    int idx;
    int maxIdx;
    int maxCount, maxCount2;
    double maxTheta, maxTheta2;
    bool mulMax, secondMax;

    // vote
    for (i=0; i<cableNum; i++)
    {
        if (cables[i].isCable && cables[i].score>.5)
        {
            // original position: +6
            idx = thetaToIdx(cables[i].theta);
            countTheta[idx] += 6;

            // +/-1: +2
            idx = thetaToIdx(cables[i].theta+1);
            countTheta[idx] += 2;
            idx = thetaToIdx(cables[i].theta-1);
            countTheta[idx] += 2;

            // +/-2: +1
            idx = thetaToIdx(cables[i].theta+2);
            countTheta[idx] += 1;
            idx = thetaToIdx(cables[i].theta-2);
            countTheta[idx] += 1;
        }
    }

    // find out the maximum; test if multiple max
    mulMax = false;
    maxIdx = 0;
    maxCount = -1;
    for (i=0; i<THETA_LEN; i++)
    {
        if (countTheta[i] > maxCount)
        {
            maxCount = countTheta[i];
            maxIdx   = i;
            mulMax = false;
        }
        else if(countTheta[i] == maxCount)
        {
            mulMax = true;
        }
    }

    // if multiple max, do nothing
    if (!mulMax)
    {
        maxTheta = THETA[maxIdx];
        
        // get the second maximum count theta
        countTheta[maxIdx] = 0;
        countTheta[(maxIdx-1+THETA_LEN)%THETA_LEN] = 0;
        countTheta[(maxIdx+1)%THETA_LEN] = 0;

        maxIdx = 0;
        maxCount2 = -1;
        for (i=0; i<THETA_LEN; i++)
        {
            if (countTheta[i] > maxCount2)
            {
                maxCount2 = countTheta[i];
                maxIdx = i;
            }
        }

        // if second maximum is large enough, keep it
        maxTheta2 = 0;
        if (maxCount2 >= int(double(maxCount)*0.5))
        {
            secondMax = true;
            maxTheta2 = THETA[maxIdx];
        }
        else
            secondMax = false;

        // screen out thetas not in maxTheta or maxTheta2
        for (i=0; i<cableNum; i++)
        {
            if (cables[i].isCable && cables[i].score>.5)
            {
                if (abs(cables[i].theta-maxTheta)>1)
                    if (!secondMax || abs(cables[i].theta-maxTheta2)>1)
                        cables[i].isCable = false;
            }
        }
    }

    free(countTheta);
}

int thetaToIdx(double theta)
{
    int thetaI = (int)theta;
    int idx;

    if (thetaI<0)
        idx = thetaI + 90;
    else
        idx = thetaI - 29;

    if (idx<0 || idx>= THETA_LEN)
        idx = (idx+THETA_LEN) % THETA_LEN;

    return idx;
}

double getAvg(double *RP_dBm, int M, int N,
    int ys, int xs,
    int r, int ss)
{
    double sum = 0;
    int num = 0;
    int x, y;

    for (y = ys-r; y <= ys+r; y += ss)
    {
        for (x = xs-r; x <= xs+r; x += ss)
        {
            if (y>=0 && y<M && x>=0 && x<N)
            {
                sum += RP_dBm[y*N+x];
                num++;
            }
        }
    }

    return sum/num;
}

int rhoTrack(double theta, int *H, double *rho, int rhoLen, int thetaLen, double *RP_dBm, int *BW, int *BWcart, double *Ang_plot, double MagMax, double MagMin, int M, int N, int NC, int x0, svmData *svm, double *rhoCable, int *isCable, int numCable, int margin)
{
    int thetaI, idx, i, *h;
    double theta_show = theta;
    double thresh;

    // get theta index in H
    if (theta<0) theta+=180;
    thetaI = (int)(theta+0.5);
    if (thetaI<PTHETA_END) thetaI = PTHETA_END;
    if (thetaI>180-PTHETA_END) thetaI = 180-PTHETA_END;
    if (thetaI>=90) idx = thetaI-90;
    else idx = thetaI-(90-PTHETA_LEN);

    // copy H(theta) into a separate array
    h = (int*)calloc(rhoLen, sizeof(int));
    for (i=0; i<rhoLen; i++)
    {
        h[i] = H[i*thetaLen+idx];
    }
    thresh = get_rho_thresh(h, rhoLen);

    // process tracker list here
    process_rtrackers(g_rtrackers_, theta, h, rho, rhoLen, thetaLen, RP_dBm, BW, BWcart, Ang_plot, M, N, NC, x0, svm, margin);

    // add new trackers
    int num_trackers = add_rtrackers(&g_rtrackers_, h, rhoLen, theta_show, rho, RP_dBm, BW, Ang_plot, M, N,x0, svm, margin, thresh);
    int num_trackers2 = get_rtracker_num(g_rtrackers_);

    printf("# trackers = %d\n", num_trackers2);
    
    // show tracker image, output rho
    double *rho_show = (double*)calloc(num_trackers2, sizeof(double));
    int *is_show = (int*)calloc(num_trackers2, sizeof(int));
    get_rtracker_rho(g_rtrackers_, rho_show, is_show);
    fprintf(pFirstRho, "%f\n", g_rtrackers_->next_ ? g_rtrackers_->next_->rho : 0);

    if (g_saveImage)
    {
#ifdef _DEBUG
        showImageIntRho(BWcart, M, NC, M > 4000 ? .25 : .5, M > 4000 ? .25 : .5, 0, 1, COLORMAP_JET, NULL, theta_show, rho_show, is_show, num_trackers2);
#endif
    }

    //showRtrackerImage(RP_dBm, BW, Ang_plot, M, N, x0, MagMax, MagMin, g_rtrackers_, num_trackers2, rhoDetectName, rhoScore>=.5);

    rhoNum = rhoTrackerToCable(g_rtrackers_, num_trackers2, &rhoCableLine);

    free(h);
    return num_trackers;
}


int simpleThetaDetect(double theta, int *H, double *rho, int rhoLen, int thetaLen, double *RP_dBm, int *BW, int *BWcart, double *Ang_plot, int M, int N, int NC, int x0, svmData *svm, double *rhoCable, int *isCable, int numCable, int margin)
{
    int thetaI, idx, i, j, *h, iMax, dataNum;
    int numToReturn = 0;
    int lineRange[512];
    int lineAngle[512];
    double vMax, linearVal, svmVal;
    double lineData[512];
    double f[FEAT_SIZE];

    // get theta index in H
    if (theta<0) theta+=180;
    thetaI = (int)(theta+0.5);
    if (thetaI<PTHETA_END) thetaI = PTHETA_END;
    if (thetaI>180-PTHETA_END) thetaI = 180-PTHETA_END;
    if (thetaI>=90) idx = thetaI-90;
    else idx = thetaI-(90-PTHETA_LEN);

    // copy H(theta) into a separate array
    h = (int*)calloc(rhoLen, sizeof(int));
    for (i=0; i<rhoLen; i++)
    {
            h[i] = H[i*thetaLen+idx];
    }
    
    // 
    // for each one
    // find local maximum
    // suppress +/- margin
    // get data from rho/theta, etc
    // classify by SVM
    // if yes, store, +1
    for (i=0; i<numCable; i++)
    {
        // find local maximum
        vMax = -1;
        iMax = -1;
        for (j=0; j<rhoLen; j++)
        {
            if (vMax<h[j])
            {
                vMax = h[j];
                iMax = j;
            }
        }

        // suppress neighborhood
        for (j=iMax-margin; j<=iMax+margin; j++)
        {
            if (j>=0 && j<rhoLen)
                h[j] = 0;
        }

        // retrieve data
        double thetaN = theta;
        if (thetaN>=90) thetaN-=180;
        dataNum = getLineData(RP_dBm, Ang_plot, M, N,
            BW, x0, rho[iMax], thetaN,
            lineData, lineRange, lineAngle);
        computeFeature(lineData, dataNum, SPEC_SIZE, AUTO_SIZE, f, FEAT_SIZE);
        // svm
        linearVal = linearCalc(svm, f);
        svmVal = svmCalc(svm, f);

        rhoCable[i] = rho[iMax];
        isCable[i] = 0;
        if (dataNum>=15 && svmVal<0 && linearVal>=0.5)
        {
            isCable[i] = 1;
            numToReturn++;
        }
    }

    free(h);
    return numToReturn;
}
