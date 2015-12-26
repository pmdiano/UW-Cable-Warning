#include "inc.h"
#include "uwcw.h"
#include "particles.h"
#include "util.h"
using namespace std;

#define SHOW_IMAGE
#define THETA_VOTE
//#define PICK_MAX
//#define THRESH_MAX
//#define REGION_FILL_TEST


unsigned char COLORMAP_JET[64][3]=
{
    {0,     0,   143},
    {0,     0,   159},
    {0,     0,   175},
    {0,     0,   191},
    {0,     0,   207},
    {0,     0,   223},
    {0,     0,   239},
    {0,     0,   255},
    {0,    16,   255},
    {0,    32,   255},
    {0,    48,   255},
    {0,    64,   255},
    {0,    80,   255},
    {0,    96,   255},
    {0,   112,   255},
    {0,   128,   255},
    {0,   143,   255},
    {0,   159,   255},
    {0,   175,   255},
    {0,   191,   255},
    {0,   207,   255},
    {0,   223,   255},
    {0,   239,   255},
    {0,   255,   255},
    {16,   255,   239},
    {32,   255,   223},
    {48,   255,   207},
    {64,   255,   191},
    {80,   255,   175},
    {96,   255,   159},
    {112,   255,   143},
    {128,   255,   128},
    {143,   255,   112},
    {159,   255,    96},
    {175,   255,    80},
    {191,   255,    64},
    {207,   255,    48},
    {223,   255,    32},
    {239,   255,    16},
    {255,   255,     0},
    {255,   239,     0},
    {255,   223,     0},
    {255,   207,     0},
    {255,   191,     0},
    {255,   175,     0},
    {255,   159,     0},
    {255,   143,     0},
    {255,   128,     0},
    {255,   112,     0},
    {255,    96,     0},
    {255,    80,     0},
    {255,    64,     0},
    {255,    48,     0},
    {255,    32,     0},
    {255,    16,     0},
    {255,     0,     0},
    {239,     0,     0},
    {223,     0,     0},
    {207,     0,     0},
    {191,     0,     0},
    {175,     0,     0},
    {159,     0,     0},
    {143,     0,     0},
    {128,     0,     0}
};

int PARTICLE_START = 10001;
int TOTAL_START = 10001;
int FRAME_I = 0;

char threshName[256];
char transName[256];
char thetaName[256];
char houghName[256];
char thetaDetectName[256];
char rhoDetectName[256];
char firstRhoName[256];
FILE *pTheta = NULL;
FILE *pFirstRho = NULL;

extern bool theta_first_frame;
extern tparticle *theta_particles;
extern rtracker *g_rtrackers_;


// rhoTrack result
double rhoScore(0), rhoScorePrev(0);
int rhoNum(0), rhoNumPrev(0);
cableLine *rhoCableLine = NULL;
cableLine *rhoCablePrev = NULL;
FILE *pRhoFile = NULL;


IplImage *pickMax(double *RP_dBm, int M, int N, int sliceHeight, int r);
IplImage *regionFillMax(double *RP_dBm, int M, int N, int sliceHeight, double T);
IplImage *threshMax(double *RP_dBm, int M,int N, double offset);

// load base file
void baseRead(const char* filename, double *Ang_endpt, double *MagMax,
              double *MagMin, double *fft_dR)
{
    FILE *fp;
    fp = fopen(filename, "rb");
    if (fp==NULL)
    {
        printf("File error: %s\n", filename);
        exit(1);
    }

    fread(Ang_endpt, sizeof(double), 1, fp);
    fread(MagMax, sizeof(double), 1, fp);
    fread(MagMin, sizeof(double), 1, fp);
    fread(fft_dR, sizeof(double), 1, fp);
    fclose(fp);
}

// load each frame
// Ang_plot and RP_dBm are allocated inside
void eachRead(const char* filename, int *M, int *N,
              double **Ang_plot, double **RP_dBm)
{
    FILE *fp;
    int i, j;
    int N1, M1;
    double *buf = NULL;
    
    fp = fopen(filename, "rb");
    if (fp==NULL)
    {
        printf("File error: %s\n", filename);
        exit(1);
    }

    fread(&N1, sizeof(int), 1, fp);
    fread(&M1, sizeof(int), 1, fp);

    *Ang_plot = (double *)malloc(sizeof(double) * (N1));
    fread(*Ang_plot, sizeof(double), N1, fp);
    
    // read main data chunk
    buf = (double*)malloc(sizeof(double)*(M1)*(N1));
    fread(buf, sizeof(double), (M1)*(N1), fp);

    *RP_dBm = (double*)malloc(sizeof(double)*M1*N1);
    for (i=0; i<M1; i++)
    {
        for (j=0; j<N1; j++)
        {
            (*RP_dBm)[i*N1+j] = buf[j*M1+i];
        }
    }

    *M = M1;
    *N = N1;

    fclose(fp);
    free(buf);
}

// TODO: CImg is not used
/*
// display an int image
void dispIntImg(int **Img, int M, int N, const char *s)
{
    CImg<double> dispImg(N, M, 1, 1, 0);

    for (int i=0; i<M; i++)
    {
        for (int j=0; j<N; j++)
        {
            dispImg(j,i) = double(Img[i][j]);
        }
    }

    CImgDisplay disp(dispImg, s);
    while (!disp.is_closed() && !disp.is_keyARROWDOWN())
        disp.wait();
}

// display a double image
void dispDoubleImg(double **Img, int M, int N, const char *s)
{
    CImg<double> dispImg(N, M, 1, 1, 0);

    for (int i=0; i<M; i++)
    {
        for (int j=0; j<N; j++)
        {
            dispImg(j,i) = Img[i][j];
        }
    }

    CImgDisplay disp(dispImg, s, 1);
    while (!disp.is_keyARROWDOWN() && !disp.is_closed())
        disp.wait();
}

// show lines overlay
void showHoughLines(int **Img, int M, int N, double *lines, int num)
{
    const unsigned char color[] = {255, 255, 0};


    CImg<unsigned char> dispImg(N, M, 1, 3, 0);
    for (int i=0; i<M; i++)
    {
        for (int j=0; j<N; j++)
        {
            if (Img[i][j])
            {
                dispImg(j, i, 0, 1)=255;
                dispImg(j, i, 0, 2)=255;
                dispImg(j, i, 0, 0)=255;
            }
        }
    }

    // draw lines
    for (int i=0; i<num; i++)
    {
        int x1=0;
        int y1=(int)((lines[2*i]-x1*cos(lines[2*i+1]*UW_PI/180))/sin(lines[2*i+1]*UW_PI/180)+0.5);
        int x2=N-1;
        int y2=(int)((lines[2*i]-x2*cos(lines[2*i+1]*UW_PI/180))/sin(lines[2*i+1]*UW_PI/180)+0.5);
        dispImg.draw_line(x1, y1, x2, y2, color);
    }


    CImgDisplay disp(dispImg);
    while (!disp.is_closed() && !disp.is_keyARROWDOWN())
        disp.wait();
}
*/
// END TODO

/*
void showOverlayCurves(int **BWBscope, int M, int N, double *lines, int num,
                       double **Bscope, double *ang, int x0, svmData *svm)
{
    const unsigned char yellow[] = {255, 255, 0};
    const unsigned char red[] = {255, 0, 0};
    const unsigned char green[] = {0, 255, 0};
    CImg<unsigned char>dispImg(N, M, 1, 3, 0);
    for (int y=0; y<M; y++)
    {
        for (int x=0; x<N; x++)
        {
            
            if (BWBscope[y][x])
            {
                for (int k=0; k<3; k++)
                {
                    dispImg(x,y,0,k)=255;
                }
            }
            
        }
    }

    // draw curves
    double data[512];
    int range[512];
    int angle[512];
    int dataNum;
    double f[FEAT_SIZE];
    double val;
    double lval;
    for (int i=0; i<num; i++)
    {
        dataNum = getLineData(Bscope, ang, M, N,BWBscope, x0,
            lines[2*i], lines[2*i+1], data, range, angle);

        computeFeature(data, dataNum, SPEC_SIZE, AUTO_SIZE, f, FEAT_SIZE);
        lval = linearCalc(svm, f);
        val = svmCalc(svm, f);

        if (val<0 && lval>0.5)
        {
            for (int j=0; j<dataNum; j++)
                dispImg.draw_point(angle[j], range[j], red);
        }
        else
        {
            for (int j=0; j<dataNum; j++)
                dispImg.draw_point(angle[j], range[j], green);
        }

    }

    CImgDisplay disp(dispImg);
    while (!disp.is_closed() && !disp.is_keyARROWDOWN())
        disp.wait();

}

*/
/*
void showOverlays(int **BWBscope, int M, int N, cableLine *cable, int num,
                       double **Bscope, double *ang, svmData *svm)
{
    const unsigned char yellow[] = {255, 255, 0};
    const unsigned char red[] = {255, 0, 0};
    const unsigned char green[] = {0, 255, 0};
    CImg<unsigned char>dispImg(N, M, 1, 3, 0);
    for (int y=0; y<M; y++)
    {
        for (int x=0; x<N; x++)
        {
            
            if (BWBscope[y][x])
            {
                for (int k=0; k<3; k++)
                {
                    dispImg(x,y,0,k)=255;
                }
            }
            
        }
    }

    // draw curves
    int *range;
    int *angle;
    int dataNum;
    for (int i=0; i<num; i++)
    {
        dataNum = cable[i].nPix;
        //dataNum = getLineData(Bscope, ang, M, N,BWBscope, x0,
        //  lines[2*i], lines[2*i+1], data, range, angle);

        //computeFeature(data, dataNum, SPEC_SIZE, AUTO_SIZE, f, FEAT_SIZE);
        //lval = linearCalc(svm, f);
        //val = svmCalc(svm, f);

        angle = cable[i].angle;
        range = cable[i].range;
        if (cable[i].isCable && cable[i].score>0.5)
        {
            for (int j=0; j<dataNum; j++)
                dispImg.draw_point(angle[j], range[j], red);
        }
        else
        {
            for (int j=0; j<dataNum; j++)
                dispImg.draw_point(angle[j], range[j], green);
        }

    }

    CImgDisplay disp(dispImg);
    while (!disp.is_closed() && !disp.is_keyARROWDOWN())
        disp.wait();
}
*/


/*void tester(void)
{
    char basename[] = "D:\\UW_Sun\\Bbase.dat";
    char filename[] = "D:\\UW_Sun\\Datasets\\B_20090901T165047d1_B\\B_20090901T165047d10002.dat";
    char svmname[]  = "D:\\UW_Sun\\svm.dat";

    int M,N;
    int i,num, sum, x0, NC, sliceSize, sliceOverlap, sliceLineNum;
    double Ang_endpt, MagMax, MagMin, fft_dR, sumD;
    double *Ang_plot = NULL;
    double **RP_dBm = NULL;
    int **BW = NULL;
    int **BWcart = NULL;
    cableLine *cable = NULL;
    double *lines = NULL;
    svmData svm;

    baseRead(basename, &Ang_endpt, &MagMax, &MagMin, &fft_dR);
    eachRead(filename, &M, &N, &Ang_plot, &RP_dBm);
    svmRead(svmname, &svm);

    BW = (int**)calloc(M, sizeof(int*));
    for (i=0; i<M; i++)
    {
        BW[i] = (int*)calloc(N, sizeof(int));
    }

    // 1. thresholding
    sliceThreshold(RP_dBm, M, N,BW, SLICE_HEIGHT, TH_PCT, TH_GATE, 0, M, allocateNew);
    sum = sumImage(BW, M, N);


    // 2. coordinate transform
    coorTransform(BW, M, N, &BWcart, &NC, &x0, Ang_endpt, Ang_plot, 4, 0, M, allocateNew);
    sum = sumImage(BWcart, M, N);

    // 3. HT and line detection
    sliceSize = (int)((double)M/16.0);
    sliceOverlap = (int)((double)sliceSize/4.0);
    sliceLineNum = 3*(int)(sliceSize/(double)(M/16.0)+0.5);
    num = detectLines(BWcart, M, NC, 0, M, THETA, THETA_LEN, sliceSize, sliceOverlap, sliceLineNum, &lines, allocateNew);

    num = removeRepLines(lines, num);

    for (i=0; i<num; i++)
    {
        printf("%f %f\n", lines[2*i], lines[2*i+1]);
    }

    double lineData[512];
    int lineRange[512];
    int lineAngle[512];
    double f[FEAT_SIZE];
    double linearVal;
    double svmVal;
    int dataNum;

    for (i=0; i<num; i++)
    {
        dataNum = getLineData(RP_dBm, Ang_plot, M, N,
            BW, x0, lines[2*i], lines[2*i+1],
            lineData, lineRange, lineAngle);
        if (dataNum < 15)
        {
            linearVal = 0;
            svmVal = 10;
        }
        computeFeature(lineData, dataNum, SPEC_SIZE, AUTO_SIZE, f, FEAT_SIZE);
        // svm
        linearVal = linearCalc(&svm, f);
        svmVal = svmCalc(&svm, f);


    }


    for (i=0; i<M; i++)
    {
        if (BW[i]) free(BW[i]);
        if (BWcart[i]) free(BWcart[i]);
    }
    if (BW) free(BW);
    if (BWcart) free(BWcart);
    if (lines) free(lines);
}
*/

/*
void main1(void)
{
    char basename[] = "D:\\UW_Sun\\Bbase.dat";
    char filename[] = "D:\\UW_Sun\\Bframe90.dat";
    char svmName[]  = "D:\\UW_Sun\\svm.dat";
    
    int M, N;
    int i, num;
    double Ang_endpt, MagMax, MagMin, fft_dR;
    double *Ang_plot=NULL;
    double **RP_dBm=NULL;
    int **BW = NULL;
    int **BWcart = NULL;
    cableLine *cable = NULL;
    svmData svm;

    baseRead(basename, &Ang_endpt, &MagMax, &MagMin, &fft_dR);
    eachRead(filename, &M, &N, &Ang_plot, &RP_dBm);
    svmRead(svmName, &svm);


    BW = (int**)malloc(sizeof(int*)*M);
    for (i=0; i<M; i++)
    {
        BW[i] = (int*)malloc(sizeof(int)*N);
    }

    sliceThreshold(RP_dBm, M, N, BW, M/16, 98, -72, 0, M, allocateNew);

    num = singleFrameDetect(RP_dBm, M, N, Ang_plot, Ang_endpt, &svm, 0, M, &cable, allocateNew);
    showOverlays(BW, M, N, cable, num, RP_dBm, Ang_plot, &svm);
    freeLineArray(cable, num);


    return;
}
*/

void paintCvImage(double *RP_dBm, IplImage *img, double magMin, double magMax, unsigned char COLORMAP[64][3])
{
    int width = img->width;
    int height = img->height;
    int widthStep = img->widthStep;
    double dataD;
    int dataI;
    uchar *imgData = (uchar*)img->imageData;
    int nChannels = img->nChannels;

    for (int i=0; i<height; i++)
    {
        for (int j=0; j<width; j++)
        {
            // dataD = *RP_dBm++;
            // dataD -= magMin;
            // dataD /= (magMax-magMin);
            // dataD *= 63;
            // dataI = (int)dataD;
            // if (dataI<0) dataI = 0;
            // if (dataI>63) dataI = 63;
            dataD = (*RP_dBm++ - magMin) / (magMax - magMin) * 63.;
            dataI = clip(dataD, 0, 63);

            imgData[j*nChannels+0] = COLORMAP[dataI][2];
            imgData[j*nChannels+1] = COLORMAP[dataI][1];
            imgData[j*nChannels+2] = COLORMAP[dataI][0];
        }
        imgData += widthStep;
    }
}

void paintCables(IplImage *img, cableLine *cables, int nCable, int flag, int flagWhite = 0, int hFactor = 2)
{
    int width = img->width;
    int height = img->height;
    int widthStep = img->widthStep;
    int nChannels = img->nChannels;
    uchar *data = (uchar*)img->imageData;
    int i, j, x, y;
    uchar R, G, B, R1, G1, B1;
    R = 255;
    G = 0; 
    B = 0;
    if (flagWhite)
    {
        G = 255;
        B = 255;
    }
    R1 = 0;
    G1 = 255;
    B1 = 0;
    if (flagWhite)
    {
        G1 = 0;
    }

    int edge = hFactor;

    #define INRANGE(X, Y) ((Y)>=0 && (Y)<height && (X)>=0 && (X)<width)

    for (i=0; i<nCable; i++)
    {
        if (cables[i].isCable && cables[i].score>.5)    // paint cable as red
        {
        }
        else if (flag)          // paint non-cable as green
        {
            for (j=0; j<cables[i].nPix; j++)
            {
                x = cables[i].angle[j];
                y = cables[i].range[j];
                for (int k=-edge; k<=edge; k++)
                {
                    if (INRANGE(x, y+k)) (data + (y+k)*widthStep)[x*nChannels + 0] = B1;
                    if (INRANGE(x, y+k)) (data + (y+k)*widthStep)[x*nChannels + 1] = G1;
                    if (INRANGE(x, y+k)) (data + (y+k)*widthStep)[x*nChannels + 2] = R1;                
                }
            }
        }
    }

    for (i=0; i<nCable; i++)
    {
        if (cables[i].isCable)  // paint cable as red
        {
            for (j=0; j<cables[i].nPix; j++)
            {
                x = cables[i].angle[j];
                y = cables[i].range[j];
                for (int k=-edge; k<=edge; k++)
                {
                    if (INRANGE(x, y+k)) (data + (y+k)*widthStep)[x*nChannels + 0] = B;
                    if (INRANGE(x, y+k)) (data + (y+k)*widthStep)[x*nChannels + 1] = G;
                    if (INRANGE(x, y+k)) (data + (y+k)*widthStep)[x*nChannels + 2] = R;             
                }
            }
        }
        else if (cables[i].score>1000)
        {
            for (j=0; j<cables[i].nPix; j++)
            {
                x = cables[i].angle[j];
                y = cables[i].range[j];
                for (int k=-edge; k<=edge; k++)
                {
                    if (INRANGE(x, y+k)) (data + (y+k)*widthStep)[x*nChannels + 0] = 0;
                    if (INRANGE(x, y+k)) (data + (y+k)*widthStep)[x*nChannels + 1] = 255;
                    if (INRANGE(x, y+k)) (data + (y+k)*widthStep)[x*nChannels + 2] = 255;               
                }
            }
        }
    }
}

void printCable(cableLine *cables, int num)
{
    for (int i=0; i<num; i++)
    {
        if (cables[i].isCable && cables[i].score > .5)
            printf("rho = %6f, theta = %3f, score = %4f\n", cables[i].rho, cables[i].theta, cables[i].score);
    }
}


void testOneDataset(
    const char  *baseName,
    const char  *eachNameBase,
    const char  *svmName,
    const int   tn,
    const char  *resultName
    )
{
    char eachName[256];
    char timeName[256];
    char oriName[256];
    char resName[256];
    int M, N, x0;
    int wFactor, hFactor;
    int j, num(0), prevNum(0);
    double Ang_endpt, MagMax, MagMin, fft_dR, score(0), prevScore(0);
    double *Ang_plot=NULL;
    double *RP_dBm=NULL;
    cableLine *cable = NULL;
    cableLine *prevCable = NULL;
    cableLine *ptcCables = NULL;
    svmData svm;
    FILE *pFile, *pTime;
    clock_t time1, time2;
    double timeb;

    char rhoResultName[256];

    M=0, N=0;
    baseRead(baseName, &Ang_endpt, &MagMax, &MagMin, &fft_dR);
    svmRead(svmName, &svm);

    pFile = fopen(resultName, "w");
    sprintf(rhoResultName, "%s_rho.txt", resultName);
    pRhoFile = fopen(rhoResultName, "w");

    sprintf(timeName, "%s%s", resultName, "-time.txt");
    pTime = fopen(timeName, "w");

#ifdef SHOW_IMAGE
    //cvNamedWindow("UWCW:original");
    //cvNamedWindow("UWCW:result");
    //cvNamedWindow("UWCW:particle");
    //cvNamedWindow("UWCW:particle_result");
#endif
    particle *ptcs = NULL, *ptcsn = NULL;
    
    int CABLE_IDX = 0;  // 4: 10082, 0; 3: 10120, 2; 0: 10001, 7; 6: 10031, 8
    cableLine *newCables = NULL;
    int maxidx = 0;
    particle *pResult = NULL;
    sprintf(thetaName, "%s_maxTheta.txt", eachNameBase);
    pTheta = fopen(thetaName, "w");
    sprintf(firstRhoName, "%s_firstRho.txt", eachNameBase);
    pFirstRho = fopen(firstRhoName, "w");

    // theta cables
    cableLine tCables[8];
    for (int ii=0; ii<8; ii++)
    {
        tCables[ii].angle = NULL;
        tCables[ii].data = NULL;
        tCables[ii].range = NULL;
    }

    for (FRAME_I=TOTAL_START; FRAME_I<=10000+tn; FRAME_I++)
    {
        if (FRAME_I==TOTAL_START) theta_first_frame = true;
        else theta_first_frame = false;

        sprintf(eachName, "%s%d.dat", eachNameBase, FRAME_I);
        printf("reading %s\n", eachName);
        eachRead(eachName, &M, &N, &Ang_plot, &RP_dBm);

#ifdef SHOW_IMAGE
        sprintf(oriName, "%s%s%d.bmp", eachNameBase, "ori_", FRAME_I);
        sprintf(resName, "%s%s%d.bmp", eachNameBase, "res_", FRAME_I);
        sprintf(threshName, "%s%s%d.bmp", eachNameBase, "thresh_", FRAME_I);
        sprintf(transName, "%s%s%d.bmp", eachNameBase, "trans_", FRAME_I);
        sprintf(houghName, "%s%s%d.bmp", eachNameBase, "hough_", FRAME_I);
        sprintf(thetaDetectName, "%s%s%d.bmp", eachNameBase, "thetaDetect_", FRAME_I);
        sprintf(rhoDetectName, "%s%s%d.bmp", eachNameBase, "rhoTrack_", FRAME_I);
#endif

        if (M>4000)
            hFactor = 4;
        else if (M>2000)
            hFactor = 2;
        else
            hFactor = 1;
        if (N>300)
            wFactor = 1;
        else if (N>150)
            wFactor = 2;
        else
            wFactor = 4;

#ifdef SHOW_IMAGE
        // show original data
        IplImage *oriImg = cvCreateImage(cvSize(N, M), IPL_DEPTH_8U, 3);
        paintCvImage(RP_dBm, oriImg, MagMin, MagMax, COLORMAP_JET);
        IplImage *dstImg = cvCreateImage(cvSize(N*wFactor,M/hFactor), IPL_DEPTH_8U, 3);
        cvResize(oriImg, dstImg);
        if (Ang_plot[0]<0)
            cvFlip(dstImg, NULL, 1);
        cvFlip(dstImg);
        //cvShowImage("UWCW:original", dstImg);
        cvSaveImage(oriName, dstImg);
        //continue;
#endif

        printf("processing..");

        time1 = clock();
        
        TH_GATE = MagMin;
        SLICE_HEIGHT = M / 16;
        SLICE_WIDTH = 0;
        SLICE_HEIGHT = M / 32;
        SLICE_WIDTH = 64/wFactor;
        num = singleFrameDetect(RP_dBm, M, N, Ang_plot, Ang_endpt, MagMax, MagMin, &svm, 0, M, &cable, &x0, tCables, pTime);

#ifdef THRESH_MAX
        IplImage *maxImg = threshMax(RP_dBm, M, N, 35);
        IplImage *dstMax = cvCreateImage(cvSize(N*wFactor,M/hFactor), IPL_DEPTH_8U, 3);
        cvResize(maxImg, dstMax);
        if (Ang_plot[0]<0)
            cvFlip(dstMax, NULL, 1);
        cvFlip(dstMax);
        cvNamedWindow("UWCW:max");
        cvShowImage("UWCW:max", dstMax);
        cvSaveImage("max.png", dstMax);
        cvWaitKey();
        cvReleaseImage(&maxImg);
        cvReleaseImage(&dstMax);
        cvDestroyWindow("UWCW:max");
#endif


#ifdef REGION_FILL_TEST
        IplImage *maxImg = regionFillMax(RP_dBm, M, N, M/16, 15);
        IplImage *dstMax = cvCreateImage(cvSize(N*wFactor,M/hFactor), IPL_DEPTH_8U, 1);
        cvResize(maxImg, dstMax);
        if (Ang_plot[0]<0)
            cvFlip(dstMax, NULL, 1);
        cvFlip(dstMax);
        cvNamedWindow("UWCW:regionFill");
        cvShowImage("UWCW:regionFill", dstMax);
        cvSaveImage("regionFill.png", dstMax);
        cvWaitKey();
        cvReleaseImage(&maxImg);
        cvReleaseImage(&dstMax);
        cvDestroyWindow("UWCW:regionFill");
#endif


#ifdef PICK_MAX
        IplImage *maxImg = pickMax(RP_dBm, M, N, M/24, 10);
        IplImage *dstMax = cvCreateImage(cvSize(N*wFactor,M/hFactor), IPL_DEPTH_8U, 1);
        cvResize(maxImg, dstMax);
        if (Ang_plot[0]<0)
            cvFlip(dstMax, NULL, 1);
        cvFlip(dstMax);
        cvNamedWindow("UWCW:max");
        cvShowImage("UWCW:max", dstMax);
        cvSaveImage("max.png", dstMax);
        cvWaitKey();
        cvReleaseImage(&maxImg);
        cvReleaseImage(&dstMax);
        cvDestroyWindow("UWCW:max");
#endif


        time2 = clock();
        timeb = double(time2-time1) / (double)CLOCKS_PER_SEC;
        /*
        printf("detected %d cables\n", num);
        printCable(cable, num);
        */

#ifdef THETA_VOTE
        voteTheta(cable, num);
#endif


        score = doubleFrameScore(cable, num, prevCable, prevNum, prevScore, FRAME_I==TOTAL_START);
        fprintf(pFile, "%f\n", score);
        fprintf(pTime, "%f\n", timeb);

        

#ifdef SHOW_IMAGE
        // show cable lines overlayyed
        IplImage *oriImg1 = cvCreateImage(cvSize(N, M), IPL_DEPTH_8U, 3);
        paintCvImage(RP_dBm, oriImg1, MagMin, MagMax, COLORMAP_JET);
        if (score>=.5) paintCables(oriImg1, cable, num, 0, 0, hFactor);
        IplImage *dstImg1 = cvCreateImage(cvSize(N*wFactor,M/hFactor), IPL_DEPTH_8U, 3);
        cvResize(oriImg1, dstImg1);
        if (Ang_plot[0]<0)
            cvFlip(dstImg1, NULL, 1);
        cvFlip(dstImg1);
        //cvShowImage(   "UWCW:result", dstImg1);
        cvSaveImage(resName, dstImg1);

        // theta tracking
        IplImage *tImg1 = cvCreateImage(cvSize(N, M), IPL_DEPTH_8U, 3);
        paintCvImage(RP_dBm, tImg1, MagMin, MagMax, COLORMAP_JET);
        paintCables(tImg1, tCables, 8, 0, 0, hFactor);
        IplImage *tImg2 = cvCreateImage(cvSize(N*wFactor, M/hFactor), IPL_DEPTH_8U, 3);
        cvResize(tImg1, tImg2);
        if (Ang_plot[0]<0)
            cvFlip(tImg2, NULL, 1);
        cvFlip(tImg2);
        cvSaveImage(thetaDetectName, tImg2);
#endif

        // particle filtering
#ifdef SHOW_IMAGE
        IplImage *oriImg2, *dstImg2;
        IplImage *oriImg3, *dstImg3;
#endif
        ptcCables = (cableLine*)calloc(NUM_PARTICLES, sizeof(cableLine));
        
        //pResult = NULL;
        if (FRAME_I==PARTICLE_START)
        {
            fprintf(stdout, "\nPARTICLE FILTER initialization\n");
            ptcs = init_distribution(cable[CABLE_IDX].rho, cable[CABLE_IDX].theta, NUM_PARTICLES);
            for (j=0; j<NUM_PARTICLES; j++)
            {
                getParticleWeight(ptcs+j, RP_dBm, Ang_plot, M, N, x0, &svm, ptcCables+j);
            }
            normalizeWeights(ptcs, NUM_PARTICLES);
            pResult = (particle*)calloc(1, sizeof(particle));
            pResult->rho = ptcs[0].rho;
            pResult->theta = ptcs[0].theta;
            pResult->N = ptcs[0].N;
            pResult->data = (double*)calloc(pResult->N, sizeof(double));
            memcpy(pResult->data, ptcs[0].data, sizeof(double)*pResult->N);
        }
        else if(FRAME_I>PARTICLE_START)
        {
            fprintf(stdout, "\nPARTICLE FILTER tracking\n");
            double maxWeight = 0;
            double w;
            int maxl = -1;
            for (j=0; j<NUM_PARTICLES; j++)
            {
                ptcs[j] = transition(ptcs[j]);
                w = getParticleWeight(&(ptcs[j]), RP_dBm, Ang_plot, M, N, x0, &svm, &(ptcCables[j]), pResult);
                if (w>maxWeight)
                {
                    maxWeight = w;
                    maxl = j;
                }
            }

            //justifyWeights(ptcs, NUM_PARTICLES);

            if (pResult)
            {
                freeParticles(pResult, 1);
                pResult = NULL;
            }
            
            int maxJ = -1;
            maxidx = -1;
            for (j=0; j<NUM_PARTICLES; j++)
            {
                if (j != maxl)
                {
                    ptcCables[j].isCable = false;
                }
                else
                {
                    //ptcCables[j].isCable = true;
                    ptcCables[j].score = 1001;
                    newCables = sampleNewCable(ptcCables+j, NUM_NEW_CABLES, RHO_STD, RP_dBm, Ang_plot, M, N, x0, &svm); 
                    if (ptcCables[j].isCable != true)
                    {
                        maxJ = j;                       
                    }
                    else
                    {
                        maxidx = j;
                        pResult = (particle*)calloc(1, sizeof(particle));
                        pResult->rho = ptcCables[j].rho;
                        pResult->theta = ptcCables[j].theta;
                        pResult->N = ptcCables[j].nPix;
                        pResult->data = (double*)calloc(pResult->N, sizeof(double));
                        memcpy(pResult->data, ptcCables[j].data, pResult->N*sizeof(double));
                    }
                }
            }
            
            printCable(ptcCables, NUM_PARTICLES);
            
            
            
            // resample
            normalizeWeights(ptcs, NUM_PARTICLES);
            ptcsn = resample(ptcs, NUM_PARTICLES);
            freeParticles(ptcs, NUM_PARTICLES);
            free(ptcs);
            ptcs = ptcsn;
            
            
            if (maxJ > -1)
            {
                if (updateParticleTracker(ptcs, NUM_PARTICLES, newCables, 2*NUM_NEW_CABLES))
                {
                    for (j=0; j<NUM_PARTICLES; j++)
                    {
                        getParticleWeight(ptcs+j, RP_dBm, Ang_plot, M, N, x0, &svm, ptcCables+j);
                    }
                    normalizeWeights(ptcs, NUM_PARTICLES);
                    printf("\ntracker UPDATED SUCCESS\n");
                }
                else
                {
                    printf("tracker UPDATE FAILURE\n");
                }
            }
            
            
            // only let the first be red
            
        }

#ifdef SHOW_IMAGE
        // show particles
        oriImg2 = cvCreateImage(cvSize(N, M), IPL_DEPTH_8U, 3);
        paintCvImage(RP_dBm, oriImg2, MagMin, MagMax, COLORMAP_JET);
        if (newCables)
            paintCables(oriImg2, newCables, 2*NUM_NEW_CABLES, 1, 0, hFactor);
        paintCables(oriImg2, ptcCables, NUM_PARTICLES, 1, 1, hFactor);
        dstImg2 = cvCreateImage(cvSize(N*wFactor, M/hFactor), IPL_DEPTH_8U, 3);
        cvResize(oriImg2, dstImg2);
        if (Ang_plot[0]<0)
            cvFlip(dstImg2,NULL,1);
        cvFlip(dstImg2);
        //cvShowImage("UWCW:particle", dstImg2);

        oriImg3 = cvCreateImage(cvSize(N, M), IPL_DEPTH_8U, 3);
        paintCvImage(RP_dBm, oriImg3, MagMin, MagMax, COLORMAP_JET);
        if (newCables)
            paintCables(oriImg3, newCables, 2*NUM_NEW_CABLES, 0, 0, hFactor);
        paintCables(oriImg3, ptcCables, NUM_PARTICLES, 0, 1, hFactor);
        dstImg3 = cvCreateImage(cvSize(N*wFactor, M/hFactor), IPL_DEPTH_8U, 3);
        cvResize(oriImg3, dstImg3);
        if (Ang_plot[0]<0)
            cvFlip(dstImg3,NULL,1);
        cvFlip(dstImg3);
        //cvShowImage("UWCW:particle_result", dstImg3);
#endif

        
        freeLineArray(prevCable, prevNum);
        freeLineArray(ptcCables, NUM_PARTICLES);
        freeLineArray(newCables, 2*NUM_NEW_CABLES);
        prevCable = cable;
        cable = NULL;
        prevNum = num;
        prevScore = score;

        // rho track free
        freeLineArray(rhoCablePrev, rhoNumPrev);
        rhoCablePrev = rhoCableLine;
        rhoCableLine = NULL;
        rhoNumPrev = rhoNum;
        rhoScorePrev = rhoScore;

        free(Ang_plot);
        free(RP_dBm);

#ifdef SHOW_IMAGE
        // opencv data
        cvWaitKey(0);
        cvReleaseImage(&oriImg);
        cvReleaseImage(&dstImg);
        cvReleaseImage(&oriImg1);
        cvReleaseImage(&dstImg1);
        cvReleaseImage(&oriImg2);
        cvReleaseImage(&dstImg2);
        cvReleaseImage(&oriImg3);
        cvReleaseImage(&dstImg3);
        cvReleaseImage(&tImg1);
        cvReleaseImage(&tImg2);
#endif
    }

    for (int ii=0; ii<8; ii++)
    {
        if (tCables[ii].angle) free(tCables[ii].angle);
        if (tCables[ii].data) free(tCables[ii].data);
        if (tCables[ii].range) free(tCables[ii].range);
    }

    svmFree(&svm);
    fclose(pFile);
    fclose(pTime);
    fclose(pRhoFile);
    // clean tracker list
    rtracker_free_list(g_rtrackers_);
    g_rtrackers_ = NULL;
    
    freeLineArray(prevCable, prevNum);

    freeParticles(ptcs, NUM_PARTICLES);
    freeParticles(ptcsn, NUM_PARTICLES);
    freeParticles(pResult, 1);
    
#ifdef SHOW_IMAGE
    //cvDestroyWindow("UWCW:original");
    //cvDestroyWindow("UWCW:result");
    //cvDestroyWindow("UWCW:particles");
    //cvDestroyWindow("UWCW:particles_result");
#endif
    fclose(pTheta);
    fclose(pFirstRho);
    if (theta_particles) theta_free_particles(theta_particles, THETA_NUM_PARTICLES);
    theta_particles = NULL;
}
void testAll(int in);


int main(int argc, const char* argv[])
{
    //tester();
    //main1();
    int i;
    i = atoi(argv[1]);
    if (argc > 2)
        TOTAL_START = atoi(argv[2]);
    if (argc > 3)
        PARTICLE_START = atoi(argv[3]);
    //srand(time(NULL));
    testAll(i);

    return 0;
}

void testAll(int in)
{
    const char* baseName[] = {
        "F:\\UW_Sun\\Datasets\\B_20090831T134940d1_A\\Prm_B_20090831T134940d1.dat",
        "F:\\UW_Sun\\Datasets\\B_20090831T135234d1_C\\Prm_B_20090831T135234d1.dat",
        "F:\\UW_Sun\\Datasets\\B_20090901T164646d1_D\\Prm_B_20090901T164646d1.dat",
        "F:\\UW_Sun\\Datasets\\B_20090901T165047d1_B\\Prm_B_20090901T165047d1.dat",
        "F:\\UW_Sun\\Datasets\\B_20090902T101724d1_G\\Prm_B_20090902T101724d1.dat",
        "F:\\UW_Sun\\Datasets\\B_20090902T102617d1_E\\Prm_B_20090902T102617d1.dat",
        "F:\\UW_Sun\\Datasets\\B_20090902T103246d1_F\\Prm_B_20090902T103246d1.dat",
        "F:\\UW_Sun\\Datasets\\B_20090902T103810d1_H\\Prm_B_20090902T103810d1.dat",
        "F:\\UW_Sun\\Datasets\\SBc20110902102717c1_I\\Prm_SBc20110902102717c1.dat",
        "F:\\UW_Sun\\Datasets\\SBc20110902104628c1_J\\Prm_SBc20110902104628c1.dat",
        "F:\\UW_Sun\\Datasets\\CFa20110929162847c1_K\\Prm_CFa20110929162847c1.dat",
        "F:\\UW_Sun\\Datasets\\CFb20110930111304c1_L\\Prm_CFb20110930111304c1.dat",
        "F:\\UW_Sun\\Datasets\\CFb20110930114159c1_M\\Prm_CFb20110930114159c1.dat",
        "F:\\UW_Sun\\Datasets\\CFd20111004141759c1_N\\Prm_CFd20111004141759c1.dat"};

    const char* eachNameBase[] = {
        "F:\\UW_Sun\\Datasets\\B_20090831T134940d1_A\\B_20090831T134940d",
        "F:\\UW_Sun\\Datasets\\B_20090831T135234d1_C\\B_20090831T135234d",
        "F:\\UW_Sun\\Datasets\\B_20090901T164646d1_D\\B_20090901T164646d",
        "F:\\UW_Sun\\Datasets\\B_20090901T165047d1_B\\B_20090901T165047d",
        "F:\\UW_Sun\\Datasets\\B_20090902T101724d1_G\\B_20090902T101724d",
        "F:\\UW_Sun\\Datasets\\B_20090902T102617d1_E\\B_20090902T102617d",
        "F:\\UW_Sun\\Datasets\\B_20090902T103246d1_F\\B_20090902T103246d",
        "F:\\UW_Sun\\Datasets\\B_20090902T103810d1_H\\B_20090902T103810d",
        "F:\\UW_Sun\\Datasets\\SBc20110902102717c1_I\\SBc20110902102717c",
        "F:\\UW_Sun\\Datasets\\SBc20110902104628c1_J\\SBc20110902104628c",
        "F:\\UW_Sun\\Datasets\\CFa20110929162847c1_K\\CFa20110929162847c",
        "F:\\UW_Sun\\Datasets\\CFb20110930111304c1_L\\CFb20110930111304c",
        "F:\\UW_Sun\\Datasets\\CFb20110930114159c1_M\\CFb20110930114159c",
        "F:\\UW_Sun\\Datasets\\CFd20111004141759c1_N\\CFd20111004141759c"};


    const char* resultName[] = {
        "F:\\UW_Sun\\1Aresult.txt",
        "F:\\UW_Sun\\2Cresult.txt",
        "F:\\UW_Sun\\3Dresult.txt",
        "F:\\UW_Sun\\4Bresult.txt",
        "F:\\UW_Sun\\5Gresult.txt",
        "F:\\UW_Sun\\6Eresult.txt",
        "F:\\UW_Sun\\7Fresult.txt",
        "F:\\UW_Sun\\8Hresult.txt",
        "F:\\UW_Sun\\9Iresult.txt",
        "F:\\UW_Sun\\9Jresult.txt",
        "F:\\UW_Sun\\9Kresult.txt",
        "F:\\UW_Sun\\9Lresult.txt",
        "F:\\UW_Sun\\9Mresult.txt",
        "F:\\UW_Sun\\9Nresult.txt"};

    int tn[] = {49, 97, 149, 147, 148, 119, 118, 147, 1193, 1194, 297, 596, 1393, 596};
        
    const char svmName[]  = "F:\\UW_Sun\\svm.dat";

    int i;
    for (i=in; i<in+1; i++)
    {
        testOneDataset(baseName[i], eachNameBase[i], svmName, tn[i], resultName[i]);
    }
}

IplImage *pickMax(double *RP_dBm, int M, int N, int sliceHeight, int r)
{
    int i, j;
    int maxI, maxJ;
    double maxV;
    int m1, m2;

    IplImage *img = cvCreateImage(cvSize(N, M), IPL_DEPTH_8U, 1);
    uchar *imgData = (uchar*)img->imageData;
    int stepSize = img->widthStep;
    
    m1 = 0;
    m2 = m1+sliceHeight;
    m2 = min(m2, M);

    while (m1<M)
    {
        maxV = -100;
        for (i=m1; i<m2; i++)
        {
            for (j=0; j<N; j++)
            {
                imgData[i*stepSize+j] = 0;
                if (RP_dBm[i*N+j] > maxV)
                {
                    maxV = RP_dBm[i*N+j];
                    maxI = i;
                    maxJ = j;
                }
            }
        }

        for (i=max(0, maxI-r); i<=min(M-1, maxI+r); i++)
        {
            for (j=max(0, maxJ-r); j<=min(N-1, maxJ+r); j++)
            {
                imgData[i*stepSize+j] = 255;
            }
        }

        m1 = m2;
        m2 = m1 + sliceHeight;
        m2 = min(m2, M);
    }

    return img;
}

IplImage *threshMax(double *RP_dBm, int M,int N, double offset)
{
    int i, j;
    double T;
    T = -1000;
    double av = 0;

    IplImage *img = cvCreateImage(cvSize(N, M), IPL_DEPTH_8U, 3);
    uchar *imgData = (uchar*)img->imageData;
    int stepSize = img->widthStep;
    int maxI, maxJ;
    
    double *data = RP_dBm;
    for (i=0; i<M; i++)
    {
        for (j=0; j<N; j++)
        {
            if (T<*data)
            {
                T = *data;
                maxI = i;
                maxJ = j;
            }
            av += *data;
            data++;
        }
    }

    av /= (M*N);
    printf("\nmax = %f, avg = %f\n", T, av);

    T = (3.*T+3.*av)/6;

    data = RP_dBm;
    for (i=0; i<M; i++)
    {
        for (j=0; j<N; j++)
        {
            if (*data >= T)
            {
                imgData[i*stepSize+j*3+0] = 255;
                imgData[i*stepSize+j*3+1] = 255;
                imgData[i*stepSize+j*3+2] = 255;
            }
            else
            {
                imgData[i*stepSize+j*3+0] = 0;
                imgData[i*stepSize+j*3+1] = 0;
                imgData[i*stepSize+j*3+2] = 0;
            }
            data++;
        }
    }

    cvCircle(img, cvPoint(maxJ, maxI), 10, cvScalar(0, 0, 255));

    return img;
}

IplImage *regionFillMax(double *RP_dBm, int M, int N, int sliceHeight, double T)
{
    int i, j;
    int maxI, maxJ;
    double maxV;
    int m1, m2;
    int *BW;
    int rsize = 0;

    IplImage *img = cvCreateImage(cvSize(N, M), IPL_DEPTH_8U, 1);
    uchar *imgData = (uchar*)img->imageData;
    int stepSize = img->widthStep;
    int xsize, ysize;
    
    m1 = 0;
    m2 = m1+sliceHeight;
    m2 = min(m2, M);

    BW = (int*)calloc(M*N, sizeof(int));

    if (M>4000 && N > 300)
    {
        size_lo = SIZE_LO_BASE * 4;
        size_hi = SIZE_HI_BASE * 4;
    }
    else if (M>4000 || N>300)
    {
        size_lo = SIZE_LO_BASE * 2;
        size_hi = SIZE_HI_BASE * 2;
    }

    while (m1<M)
    {
        maxV = -1000;
        for (i=m1; i<m2; i++)
        {
            for (j=0; j<N; j++)
            {
                imgData[i*stepSize+j] = 0;
                if (RP_dBm[i*N+j] > maxV)
                {
                    maxV = RP_dBm[i*N+j];
                    maxI = i;
                    maxJ = j;
                }
            }
        }

        if (maxV >= REGIONFILL_T)
        {
            rsize = regionFill(RP_dBm, BW, m1, m2, M, N, maxI, maxJ, T, &ysize, &xsize);
            if (rsize<size_lo || rsize>size_hi) memset(BW+m1*N, 0, (m2-m1)*N*sizeof(int));

            printf("region size = %d, size_lo = %d, size_hi = %d\n", rsize, size_lo, size_hi);
            printf("M=%d, N=%d, ysize=%d, xsize=%d\n\n", M, N, ysize, xsize);
        }


        m1 = m2;
        m2 = m1 + sliceHeight;
        m2 = min(m2, M);

        //cvCircle(img, cvPoint(maxJ, maxI), 10, cvScalar(255), 2);
    }

    for (i=0; i<M; i++)
    {
        for (j=0; j<N; j++)
        {
            if (imgData[j]<128 && BW[i*N+j]>1)
                imgData[j] = 255;
        }
        imgData += stepSize;
    }

    return img;

}

int rhoTrackerToCable(rtracker *rtrackers_, int totalNum, cableLine **cls)
{
    int num(0), i(0), k(0);
    rtracker *current(rtrackers_);

    // get the number of cables
    while (current && !current->rho) current = current->next_;
    for (i=0; i<totalNum; i++)
    {
        if (current->rstatus == dormant)
            ;
        else
            num++;
        
        current = current->next_;
    }

    // transfer data, only status and score are important
    *cls = (cableLine*)calloc(num, sizeof(cableLine));
    current = rtrackers_;
    while (current && !current->rho) current = current->next_;
    for (i=0; i<totalNum; i++)
    {
        if (current->rstatus == dormant)
            ;
        else
        {
            (*cls)[k].isCable = true;
            (*cls)[k].score = current->rpts_[0].linearVal;
            (*cls)[k].rho = current->rho;
            (*cls)[k].theta = current->theta;
            k++;
        }

        current = current->next_;
    }

    return num;
}

void showRtrackerImage(double *RP_dBm, int *BW, double *Ang_plot, int M, int N, int x0,
    double MagMax, double MagMin, rtracker *rtrackers_, int totalNum, char *nameToSave, bool showCable)
{
    int dataNum;
    double lineData[512];
    int lineRange[512];
    int lineAngle[512];
    int wFactor;
    int hFactor;

    if (M>4000)
        hFactor = 4;
    else if (M>2000)
        hFactor = 2;
    else
        hFactor = 1;
    if (N>300)
        wFactor = 1;
    else
        wFactor = 2;

    cableLine *cls = (cableLine*)calloc(totalNum, sizeof(cableLine));

    while(rtrackers_ && !rtrackers_->rho) rtrackers_ = rtrackers_->next_;
    for (int i=0; i<totalNum; i++)
    {
        if (rtrackers_->rstatus == dormant)
        {
            cls[i].isCable = false;
            cls[i].score = 0;
        }
        else
        {
            cls[i].isCable = true;
            cls[i].score = 1;
        }

        dataNum = getLineData(RP_dBm, Ang_plot, M, N, BW, x0, rtrackers_->rho, rtrackers_->theta, lineData, lineRange, lineAngle, true);
        cls[i].rho = rtrackers_->rho;
        cls[i].theta = rtrackers_->theta;
        cls[i].nPix = dataNum;

        cls[i].range = (int*)calloc(dataNum, sizeof(int));
        cls[i].angle = (int*)calloc(dataNum, sizeof(int));
        memcpy(cls[i].range, lineRange, dataNum*sizeof(int));
        memcpy(cls[i].angle, lineAngle, dataNum*sizeof(int));

        rtrackers_ = rtrackers_->next_;
    }

    
    IplImage *tImg1 = cvCreateImage(cvSize(N, M), IPL_DEPTH_8U, 3);
    paintCvImage(RP_dBm, tImg1, MagMin, MagMax, COLORMAP_JET);
    if (showCable) paintCables(tImg1, cls, totalNum, 0, 0, hFactor);
    IplImage *tImg2 = cvCreateImage(cvSize(N*wFactor, M/hFactor), IPL_DEPTH_8U, 3);
    cvResize(tImg1, tImg2);
    if (Ang_plot[0]<0)
        cvFlip(tImg2, NULL, 1);
    cvFlip(tImg2);
    cvSaveImage(nameToSave, tImg2);
}
