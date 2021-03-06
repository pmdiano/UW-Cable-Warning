#include "inc.h"
#include "uwcw.h"
#include "particles.h"
#include "util.h"
using namespace std;

#define THETA_VOTE
//#define PICK_MAX
//#define THRESH_MAX
//#define REGION_FILL_TEST

bool g_saveImage = false;

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
    const char  *resultName,
    const char  *imageNameBase
    )
{
    char eachName[256];
    char timeName[256];
    char oriName[256];
    char resName[256];
    char particleName[256];
    char particleResultName[256];
    int M, N, x0;
    int wFactor, hFactor;
    int num(0), prevNum(0);
    double Ang_endpt, MagMax, MagMin, fft_dR, score(0), prevScore(0);
    double *Ang_plot=NULL;
    double *RP_dBm=NULL;
    cableLine *cable = NULL;
    cableLine *prevCable = NULL;
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

        sprintf(oriName, "%s%s%d.bmp", imageNameBase, "_ori_", FRAME_I);
        sprintf(resName, "%s%s%d.bmp", imageNameBase, "_res_", FRAME_I);
        sprintf(threshName, "%s%s%d.bmp", imageNameBase, "_thresh_", FRAME_I);
        sprintf(transName, "%s%s%d.bmp", imageNameBase, "_trans_", FRAME_I);
        sprintf(houghName, "%s%s%d.bmp", imageNameBase, "_hough_", FRAME_I);
        sprintf(thetaDetectName, "%s%s%d.bmp", imageNameBase, "_thetaDetect_", FRAME_I);
        sprintf(rhoDetectName, "%s%s%d.bmp", imageNameBase, "_rhoTrack_", FRAME_I); // This is the final result
        sprintf(particleName, "%s%s%d.bmp", imageNameBase, "_particle_", FRAME_I);
        sprintf(particleResultName, "%s%s%d.bmp", imageNameBase, "_particleResult_", FRAME_I);

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

        if (g_saveImage)
        {
            // save original data
            IplImage *oriImg = cvCreateImage(cvSize(N, M), IPL_DEPTH_8U, 3);
            paintCvImage(RP_dBm, oriImg, MagMin, MagMax, COLORMAP_JET);
            IplImage *dstImg = cvCreateImage(cvSize(N*wFactor, M / hFactor), IPL_DEPTH_8U, 3);
            cvResize(oriImg, dstImg);
            if (Ang_plot[0] < 0)
                cvFlip(dstImg, NULL, 1);
            cvFlip(dstImg);
            cvSaveImage(oriName, dstImg);

            cvReleaseImage(&oriImg);
            cvReleaseImage(&dstImg);
        }

        printf("processing..");

        time1 = clock();
        
        TH_GATE = MagMin;
        SLICE_HEIGHT = M / 16;
        SLICE_WIDTH = 0;
        SLICE_HEIGHT = M / 32;
        SLICE_WIDTH = 64/wFactor;
        num = singleFrameDetect(RP_dBm, M, N, Ang_plot, Ang_endpt, MagMax, MagMin, &svm, 0, M, &cable, &x0, tCables, pTime);

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

        if (g_saveImage)
        {
            // show cable lines overlayyed
            IplImage *oriImg1 = cvCreateImage(cvSize(N, M), IPL_DEPTH_8U, 3);
            paintCvImage(RP_dBm, oriImg1, MagMin, MagMax, COLORMAP_JET);
            if (score >= .5) paintCables(oriImg1, cable, num, 0, 0, hFactor);
            IplImage *dstImg1 = cvCreateImage(cvSize(N*wFactor, M / hFactor), IPL_DEPTH_8U, 3);
            cvResize(oriImg1, dstImg1);
            if (Ang_plot[0] < 0)
                cvFlip(dstImg1, NULL, 1);
            cvFlip(dstImg1);
            cvSaveImage(resName, dstImg1);

            // theta tracking
            IplImage *tImg1 = cvCreateImage(cvSize(N, M), IPL_DEPTH_8U, 3);
            paintCvImage(RP_dBm, tImg1, MagMin, MagMax, COLORMAP_JET);
            paintCables(tImg1, tCables, 8, 0, 0, hFactor);
            IplImage *tImg2 = cvCreateImage(cvSize(N*wFactor, M / hFactor), IPL_DEPTH_8U, 3);
            cvResize(tImg1, tImg2);
            if (Ang_plot[0] < 0)
                cvFlip(tImg2, NULL, 1);
            cvFlip(tImg2);
            cvSaveImage(thetaDetectName, tImg2);

            cvReleaseImage(&oriImg1);
            cvReleaseImage(&dstImg1);
            cvReleaseImage(&tImg1);
            cvReleaseImage(&tImg2);
        }

        freeLineArray(prevCable, prevNum);
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
    freeLineArray(rhoCablePrev, rhoNumPrev);

    fclose(pTheta);
    fclose(pFirstRho);
    if (theta_particles) theta_free_particles(theta_particles, THETA_NUM_PARTICLES);
    theta_particles = NULL;
}

void testAll(int start, int num);

static const char szHelpInfo[] =
    "uwcw dataset_start dataset_num [total_start] [particle_start] [-s]\n"
    "\n"
    "  dataset_start:\n"
    "    start index of datasets to be tested (from 0)\n"
    "\n"
    "  dataset_num:\n"
    "    total number of datasets to be tested\n"
    "\n"
    "  total_start:\n"
    "    start frame for each dataset to be processed (default: 10001)\n"
    "\n"
    "  particle_start:\n"
    "    start frame for each dataset to use particle filter (default: 10001)\n"
    "\n"
    "  -s\n"
    "    save intermediate result images\n";

int main(int argc, const char* argv[])
{
    int start, num;
    bool parsed = false;

    if (argc >= 3)
    {
        start = atoi(argv[1]);
        num = atoi(argv[2]);
        parsed = true;

        int i = 3;
        if (argc >= 5)
        {
            TOTAL_START = atoi(argv[i++]);
            PARTICLE_START = atoi(argv[i++]);
        }

        if (argc > i && !strcmp(argv[i], "-s"))
        {
            g_saveImage = true;
        }
    }

    if (parsed)
    {
        srand(time(NULL));
        testAll(start, num);
    }
    else
    {
        printf("%s", szHelpInfo);
    }

    return 0;
}

void testAll(int start, int num)
{
    const char* baseName[] = {
        "../Datasets/B_20090831T134940d1_A/Prm_B_20090831T134940d1.dat",
        "../Datasets/B_20090831T135234d1_C/Prm_B_20090831T135234d1.dat",
        "../Datasets/B_20090901T164646d1_D/Prm_B_20090901T164646d1.dat",
        "../Datasets/B_20090901T165047d1_B/Prm_B_20090901T165047d1.dat",
        "../Datasets/B_20090902T101724d1_G/Prm_B_20090902T101724d1.dat",
        "../Datasets/B_20090902T102617d1_E/Prm_B_20090902T102617d1.dat",
        "../Datasets/B_20090902T103246d1_F/Prm_B_20090902T103246d1.dat",
        "../Datasets/B_20090902T103810d1_H/Prm_B_20090902T103810d1.dat",
        "../Datasets/SBc20110902102717c1_I/Prm_SBc20110902102717c1.dat",
        "../Datasets/SBc20110902104628c1_J/Prm_SBc20110902104628c1.dat",
        "../Datasets/CFa20110929162847c1_K/Prm_CFa20110929162847c1.dat",
        "../Datasets/CFb20110930111304c1_L/Prm_CFb20110930111304c1.dat",
        "../Datasets/CFb20110930114159c1_M/Prm_CFb20110930114159c1.dat",
        "../Datasets/CFd20111004141759c1_N/Prm_CFd20111004141759c1.dat"
    };

    const char* eachNameBase[] = {
        "../Datasets/B_20090831T134940d1_A/B_20090831T134940d",
        "../Datasets/B_20090831T135234d1_C/B_20090831T135234d",
        "../Datasets/B_20090901T164646d1_D/B_20090901T164646d",
        "../Datasets/B_20090901T165047d1_B/B_20090901T165047d",
        "../Datasets/B_20090902T101724d1_G/B_20090902T101724d",
        "../Datasets/B_20090902T102617d1_E/B_20090902T102617d",
        "../Datasets/B_20090902T103246d1_F/B_20090902T103246d",
        "../Datasets/B_20090902T103810d1_H/B_20090902T103810d",
        "../Datasets/SBc20110902102717c1_I/SBc20110902102717c",
        "../Datasets/SBc20110902104628c1_J/SBc20110902104628c",
        "../Datasets/CFa20110929162847c1_K/CFa20110929162847c",
        "../Datasets/CFb20110930111304c1_L/CFb20110930111304c",
        "../Datasets/CFb20110930114159c1_M/CFb20110930114159c",
        "../Datasets/CFd20111004141759c1_N/CFd20111004141759c"
    };

    const char* resultNameBase[] = {
        "../test/1A",
        "../test/2C",
        "../test/3D",
        "../test/4B",
        "../test/5G",
        "../test/6E",
        "../test/7F",
        "../test/8H",
        "../test/9I",
        "../test/9J",
        "../test/9K",
        "../test/9L",
        "../test/9M",
        "../test/9N"
    };

    int tn[] = {49, 97, 149, 147, 148, 119, 118, 147, 1193, 1194, 297, 596, 1393, 596};

    const char svmName[]  = "./svm.dat";
    char resultName[256];

    for (int i=start; i<start+num; i++)
    {
        sprintf(resultName, "%sresult.txt", resultNameBase[i]);
        testOneDataset(baseName[i], eachNameBase[i], svmName, tn[i], resultName, resultNameBase[i]);
    }
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

    freeLineArray(cls, totalNum);
}
