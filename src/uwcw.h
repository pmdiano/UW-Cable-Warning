#ifndef __UWCW_H__
#define __UWCW_H__

// common parameters
#define UW_PI           3.14159265
#define CABLE_TIP_END 15    /** number of pixels to exclude on each tip */

#define SIZE_LO_BASE 128
#define SIZE_HI_BASE 512
#define XSIZE_LO 14
#define XSIZE_HI 32
#define YSIZE_LO 18
#define YSIZE_HI 48

extern int size_lo;
extern int size_hi;
extern int tower_removed;

// algorithm parameters
extern double TH_PCT;                   // thresholding percentage
extern double TH_PCT2;
extern double REGIONFILL_T;
extern int SLICE_HEIGHT;
extern int SLICE_WIDTH;                 // in slice-wise thresholding
//#define TH_GATE       -72             // noise floor in thresholding
extern double TH_GATE;
extern double THETA[];
extern int THETA_LEN;
extern int LOWBOUND_PCT;
extern int HGATE_FACTOR;
extern unsigned char COLORMAP_JET[64][3];
#define SPEC_SIZE       2048            // power spectrum DFT size
#define AUTO_SIZE       128             // autocorrelation DFT size
#define FEAT_SIZE       14              // feature size

// svm data structure
typedef struct _svmData
{
    int M;          // # of support vectors
    int N;          // feature dimension
    double **sv;
    double *alpha;
    double bias;
    double *shift;
    double *scale;
    double *svss;   // sum of squares of sv
    double *b;      // numel: 2N+1
} svmData;

// cable line data structure
typedef struct _cableLine
{
    bool isCable;   // svm result
    double score;   // linear combination score
    double rho;     // rho and theta of BWcart
    double theta;
    int nPix;       // number of pixels on this line
    double *data;   // data on this line
    int *range;
    int *angle;
} cableLine;

//////////////////////////////////////////////////////////////////////////////
// slice-thresholding related function
//////////////////////////////////////////////////////////////////////////////

// slice thresholding : main function
void sliceThreshold(double *Bscope, // input Bscope image
                    int M,              // number of rows
                    int N,              // number of columns
                    int *BWBscope,      // output black-n-white image, shall be pre-allocated outside
                    int sliceHeight,    // slice height
                    double percentage,  // percentage of threshold
                    double threshGate,  // threshold gate
                    int y1,             // only process from y1 to y2
                    int y2,
                    int sliceWidth = 0);// slice width

//blockThreahold(Bscope, M, N, BWBscope, sliceHeight, sliceWidth, percentage, threshGate, y1, y2);
void blockThreshold(double *Bscope,
                    int M,
                    int N,
                    int *BWBscope,
                    int sliceHeight,
                    int sliceWidth,
                    double percentage,
                    double threshGate,
                    int y1,
                    int y2);

// get the certain percentage value in a slice
double slicePercent(double *Bscope, // input Bscope image
                    int M,              // number of rows
                    int N,              // number of columns
                    int h1,             // start row
                    int h2,             // end row
                    double percentage);
//T = blockPercent(Bscope, M, N, h1, h2, w1, w2, percentage);
double blockPercent(double *Bscope,
                    int M,  int N,
                    int h1, int h2,
                    int w1, int w2,
                    double percentage,
                    double *T2 = NULL,
                    double *p2 = NULL);

// quick sorting
int partition(double *data, int low, int high);
void QSort(double *data, int low, int high);


//////////////////////////////////////////////////////////////////////////////
// coordinate transformation
//////////////////////////////////////////////////////////////////////////////
void coorTransform(int *BWBscope,       // input b-n-w Bscope image
                   int M,               // # of rows of BWBscope
                   int N,               // # of columns of BWBscope
                   int *BWcart,         // output b-n-w cartesian image, shall be allocated inside
                   int *pNC,            // # of columns of BWcart
                   int *px0,            // middle point in BWcart
                   double Ang_endpt,    // parameter, angle endpoint, in degrees
                   double *Ang_plot,    // 1 by N vector angle value, in degrees
                   int interFac,        // interpolation factor
                   int y1,              // starting and ending row
                   int y2);

// assign a target pixel to 1 in coordinate transform
void coorTransAssign(int *BWcart,       // b-n-w cartesian image
                     int MC,            // # of rows of BWcart
                     int NC,            // # of columns of BWcart
                     double rho,        // rho in Bscope
                     double theta,      // theta in Bscope
                     int x0,            // middle point
                     int val = 0);          


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
                    int *H);            // output image, sizeof (rhoLen rows, thetaLen columns)

// identify peaks after HT
int houghPeaks( int *H,             // Hough data
                int rhoLen,             // # of rows in H
                int thetaLen,           // # of columns in H
                int numPeaks,           // # of peaks to find
                double hGate,           // gate value for a peak
                int nHoodRowHalf,       // neighborhood size half row
                int nHoodColHalf,       // neighborhood size half column
                int *peaks);            // output, length numPeaks*2

// a function used in houghPeaks()
void getMaxVal( int *H,
                int M,
                int N,
                double *maxVal,
                int *maxRow,
                int *maxCol);


// line detection -- connects houghTransform() and houghPeaks()
// return number of lines
// lines is allocated within
int detectLines( int *BWcart,           // input b-n-w cartesian image
                 int MC,                // size of BWcart,
                 int NC,                
                 int ystart,            // region to process
                 int yend,
                 double *theta,         // HT parameter
                 int thetaLen,
                 int sliceHeight,       // slice-processing height
                 int sliceOverlap,      // overlapping between slices
                 int numPeaks,          // # of peaks in each slice
                 double *lines,         // rho and theta for lines
                 double ratio = 1.0);   // multiplication ratio

//PTHETA_LEN, RP_dBm, BW, Ang_plot, M, N, x0, svm
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
                    double MagMax, double MagMin);


// retrieve line data
// everything is NOT allocated within
int getLineData(double *Bscope,
                double *Ang_plot,
                int M,                  // # of rows in Bscope
                int N,                  // # of columns; length of Ang_plot
                int *BWBscope,
                int x0,
                double rho,
                double theta,           // theta in degrees format
                double *data,           // retrieved data
                int *range,             // range index
                int *angle,             // angle index
                bool makeLong = false);



//////////////////////////////////////////////////////////////////////////////
// feature-computation related function
//////////////////////////////////////////////////////////////////////////////

// everything is allocated outside
void powerSpectrum(double *input,       // input data
                   int size,            // data size (which is also DFT size)
                   double *spec);       // output power spectrum


void autoCorrelation(double *input,
                     int size,
                     double *autoCo);

void computeFeature(double *data,       // input data
                    int dataSize,       // data size
                    int specSize,       // DFT size for computing power spectrum
                    int autoSize,       // DFT size for computing autocorrelation
                    double *feature,    // output feature
                    int featSize);      // feature size (14)

// sum array from data[n1] to data[n2-1]
double sumArray(double *data,
                int n1,
                int n2);

// sum array
int sumArrayInt(int *data, int n1, int n2);

// threshold feature group
void thresholdFeature(double *data,     // data and size
                      int size,
                      double T,         // threshold
                      double *f1,       // three features
                      double *f2,
                      double *f3);

int compareThreshold(double *data,
                     int dataLen,
                     double T,
                     int *num,
                     int *loc);

double uwStd(double *data, int dataLen);

void printArray(double *data, int size);


//////////////////////////////////////////////////////////////////////////////
// SVM related function
//////////////////////////////////////////////////////////////////////////////
void svmRead(const char* filename, svmData *svm);
void svmFree(svmData *svm);

// calculate the svm value; notice that feature will change value
double svmCalc(svmData *svm, double *f);
double linearCalc(svmData *svm, double *f);


void freeLines(cableLine *cL);
void freeLineArray(cableLine *cLs, int num);


// cLs will be located within
// returns number of cables found
int singleFrameDetect(
    double *RP_dBm,
    int M,         // # of rows of RP_dBm
    int N,         // # of columns of RP_dBm
    double *Ang_plot,  // angle values
    double Ang_endpt,  // angle value end
    double MagMax, double MagMin,
    svmData *svm,  // svm
    int y1,        // starting and ending point
    int y2,        // predict processing region
    cableLine **cLs,   // detection results
    int *px0,
    cableLine *tCables,
    FILE *pF = NULL // theta tracking cable results, 8 of them
    );

// calculate the frame score from detection results of two frames
double doubleFrameScore(
    cableLine *current,
    int currentNum,
    cableLine *prev,
    int prevNum,
    double prevScore,
    bool firstFrame
    );

// get number of cables
int numCable(cableLine *cl, int num);

// number of parallel cables in one frame
int checkParallel(cableLine *cl, int num);

double lineCorrelate(cableLine l1, cableLine l2, 
                     double Rsigma, double Tsigma, double C);
int temporalCorrelate(cableLine *current, int currentNum,
                     cableLine *prev, int prevNum);

int removeRepLines(double *cLs, int num);

void showImageBW(int *data, int height, int width, double hFactor, double wFactor, int flipLR = 0);
void showImageInt(int *data, int height, int width, double hFactor, double wFactor, int flipLR = 0, int scale = 1, unsigned char COLORMAP[64][3] = COLORMAP_JET, const char *nm = NULL, double theta = -100, bool bHough = false);
void showImageIntRho(int *data, int height, int width, double hFactor, double wFactor, int flipLR, int scale, unsigned char COLORMAP[64][3], const char *nm, double theta, double rhos[], int isCable[], int rhoNum);

void voteTheta(cableLine *cables, int cableNum);

int thetaToIdx(double theta);

void rejectThreshOutlier(int *BW, int *BW2, int M, int N);

int regionFill(double *RP_dBm, int *BW, // input and output; 2 means filled; return region size
    int y1, int y2, int M, int N,       // image size and processing range
    int ys, int xs,                     // seed position
    double T,                           // continuity constraint threshold
    int *ysize, int *xsize,             // return y dimension size and x dimension size
    int Nnb = 8);                       // connectivity, default 8-connected

double getAvg(double *RP_dBm, int M, int N,
    int ys, int xs,
    int r, int ss = 1);

//removeTower(RP_dBm, BW, M, N, M/16, 15, size_lo, size_hi);
void removeTower(double *RP_dBm, int *BW,
    int M, int N, int sliceHeight, double T,
    int size_lo, int size_hi);

int dcomp(const void *a, const void *b);


/// tracking-coupled detection algorithms
int simpleThetaDetect(double theta, int *H, double *rho, int rhoLen, int thetaLen, double *RP_dBm, int *BW, int *BWcart, double *Ang_plot, int M, int N, int NC, int x0, svmData *svm, double *rhoCable, int *isCable, int numCable = 8, int margin = 15);

int rhoTrack(double theta, int *H, double *rho, int rhoLen, int thetaLen, double *RP_dBm, int *BW, int *BWcart, double *Ang_plot, double MagMax, double MagMin, int M, int N, int NC, int x0, svmData *svm, double *rhoCable, int *isCable, int numCable = 8, int margin = 15);

#endif
