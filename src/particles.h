/** @file
    Definitions related to tracking cable line with particle filtering

    Reference: http://blogs.oregonstate.edu/hess/code/particles/
*/

#ifndef __PARTICLES_H__
#define __PARTICLES_H__

#define NUM_PARTICLES   200 /** number of particles */
#define NUM_NEW_CABLES 50
#define RHO_STD         10.0
#define THETA_STD       1.0

/******************************************* Structures ***********************************************/
/**
    A particle is an instantiation of the state variables, in our case, a sample (candidate) cable.
    A collection of particles is essentially a discretization of the posterior probability of the system.
*/
typedef struct _particle
{
    double rho;         /** current rho coordinate */
    double theta;       /** current theta coordinate */
    double rhop;        /** previous rho coordinate */
    double thetap;      /** previous theta coordinate */
    double w;           /** weight */
    double *data;       /** data on the line */
    int N;              /** number of points */
} particle;


//**************************************************************************************************//
//Theta particle
#define THETA_NUM_PARTICLES 80
#define THETA_VARIANCE 0.5
typedef struct _tparticle
{
    double theta;
    double thetap;
    double thetapp;
    double w;
    int *h;
    int n;
} tparticle;

//theta_initialize_particles(maxIdx, theta_particles);
tparticle* theta_initialize_particles(int num, double theta, int* H, int rhoLen, int thetaLen, int idx);

void theta_free_particles(tparticle* theta_particles, int num);

double theta_average(tparticle* theta_particles, int num);

void theta_transition_particles(tparticle* theta_particles, int num, double theta_std);

void theta_compute_likelihood(tparticle* theta_particles, int num, int* H, int rhoLen, int thetaLen,
    double(*similarity)(int *h1, int *h2, int n));

double entropy_max3(int *h, int rhoLen, int ss);

double similarity_1(int *h1, int *h2, int n);
double similarity_cos(int *h1, int *h2, int n);

void theta_normalize_weights(tparticle *theta_particles, int num);

int tparticleCmp(const void *p1, const void *p2);

tparticle* tresample(tparticle* ptcs, int n);


//**************************************************************************************************//
// rho tracker
#define RHO_NUM_PARTICLES 30
#define SVM_BASE 0.3
#define ASSOCIATION_THRESHOLD 0.4
#define DORMANCY_THRESHOLD 1
#define HOUGH_SIDE 3

enum tracker_status{
    head,       // head for the linked list
    candidate,
    healthy,
    dormant,
    dead
};

typedef struct _rparticle
{
    double rho;
    double rhop;
    double rhopp;
    double theta;
    double thetap;
    double thetapp;
    double w;
    double hVal;
    double svmVal;
    double linearVal;
} rparticle;

typedef struct _rtracker
{
    rparticle *rpts_;// particles in this tracker
    int rpts_num;   // num of particles in this tracker

    tracker_status rstatus;
    int life_time;
    int status_time;

    // tracked rho of current time
    double rho;     
    double rhop;
    double rhopp;
    double theta;
    double thetap;
    double thetapp;
    
    // previous tracked data
    double *data_original_;
    int *data_threshold_;
    int data_num;

    // previous tracked Hough data
    int *hough_;
    int hough_num;  // number of hough points saved; be an odd number

    // to make the list of rtrackers
    struct _rtracker *next_;
} rtracker;

// rparticles
rparticle* rho_initialize_particles(int num, double rho, double theta, double hVal, double svmVal, double linearVal);
void rho_transition_particles(rparticle *rho_particles, int num, double rho_std, double theta_new, int ss, rtracker *r_tracker = NULL);
double rho_compute_likelihood(rparticle *rho_particles, int num, double rhop,
    double theta, int *H, double *rho, int rhoLen, int thetaLen,
    double *RP_dBm, int *BW, int *BWcart, double *Ang_plot,
    int M, int N, int NC, int x0, svmData *svm, rtracker *r_track);
void rho_normalize_weights(rparticle *rpts, int num);
int rparticleCmp(const void *p1, const void *p2);
rparticle* rresample(rparticle* ptcs, int n);


int rho_to_idx(double rho, double *rhos, int rhoLen);

// rtracker functions
rtracker *rtracker_init_zero();
rtracker *rtracker_init(double rho, double theta,
    double *lineData, int *lineRange, int *lineAngle, int dataNum,
    int *BW, int M, int N, double hVal, double svmVal, double linearVal,
    int *hough_data, int hnum);
void rtracker_free(rtracker *rtracker_);
int rtracker_free_list(rtracker *head);

int get_rtracker_num(rtracker *rtracker_);
int get_rtracker_rho(rtracker *rtracker_, double *rho, int *is_cable);

void process_rtrackers(rtracker *rtracker_,
    double theta, int *h, double *rho, int rhoLen, int thetaLen,
    double *RP_dBm, int *BW, int *BWcart, double *Ang_plot,
    int M, int N, int NC, int x0, svmData *svm, int margin);

int add_rtrackers(rtracker **rtrackers__, int *h, int rhoLen, double theta, double *rho, double *RP_dBm, int *BW, double *Ang_plot, int M, int N, int x0, svmData *svm, int margin, double thresh);
double get_rho_thresh(int *h, int rhoLen);



/*************************************** Function Prototypes ****************************************/

/**
    Creates an initial distribution of particles
    All the particles are identical to the input rho and theta

    @param rho is the input cable parameter rho
    @param theta is the input cable parameter theta
    @param n is the number of particles to generate; will allocate memory in this function
*/
particle* init_distribution(double rho, double theta, int n);

/**
    Normalizes particle weights so they sum to 1

    @param particles is an array of particles whose weights are to be normalized
    @param n the number of particles
*/
void normalizeWeights(particle* particles, int n);

/**
    Get the weight of a particle by getting data, computing feature, running SVM
    weight is the score of SVM. The higher weight, the more likely to be a cable
    Also store the cable data into cable data structure

    @param ptc the particle whose weight to get
    @param Bscope the Bscope data in a one-dimensional array, row-wise
    @param Ang_plot the angle data in a one-dimensional array
    @param M the height of the Bscope image
    @param N the width of the Bscope image
    @param x0 a parameter for coordinate transformation
    @param svm SVM training result
    @param cable to return the cable data
*/
double getParticleWeight(particle *ptc, double *Bscope, double *Ang_plot,
    int M, int N, int x0, svmData *svm,
    cableLine *cable, particle *pResult = NULL);            /** also store the cable line data */

/**
    generate two random gauss numbers

    @param std1 the standard deviation for y1
    @param std2 the standard deviation for y2
    @param y1 the first output gaussian number
    @param y2 the second output gaussian number
*/
void gauss(double std1, double std2, double *y1, double *y2);

/**
    Sample a transition model for a given particle

    @param p the particle being transitioned
*/
particle transition(particle p);

/**
    resample the particles

    @param ptcs old particles
    @n the number of particles

    @return Returns a new set of particles
*/
particle* resample(particle* ptcs, int n);

/**
    Compare two particles based on weight. For use in qsort.
*/
int particleCmp(const void *p1, const void *p2);

/**
    free a particle array
*/
void freeParticles(particle* ptcs, int n);

/**
    compute normalized cross correlation between two data sequence
*/
double dataCorrelation(double *data1, double *data2, int n, int dataXor = 0);


/**
    another coefficient
*/
double bhattacharyya(double *data1, double *data2, int n);
double simple_c(double *data1, double *data2, int n);

/**
    compute a similarity metric based on Battacharyya similarity coefficient
*/
double similarityCoefficient(double *data1, double *data2, int n);

/**
    sample new cables with the same angle

    @anchor: anchor reference cable line
    @num: number at one side
    @rhoStd: std of rho
*/
cableLine* sampleNewCable(cableLine *anchor, int num, double rhoStd,
    double *Bscope, double *Ang_plot,
    int M, int N, int x0, svmData *svm);

/** 
    update particle tracker
    pick the first cable from new samples to be the seed
    updateParticleTracker(ptcs, NUM_PARTICLES, newCables, 2*NUM_NEW_CABLES);
*/
int updateParticleTracker(particle* ptcs, int nParticle, cableLine* cables, int nCable);

/**
    make all weights to be non-negative
    substract minimum from all weights
*/
void justifyWeights(particle* ptcs, int nParticle);

#endif
