#include "inc.h"
#include "uwcw.h"
#include "particles.h"
using namespace std;

bool theta_first_frame = false;
tparticle *theta_particles = NULL;
double theta_previous = -100;
extern int PTHETA_END;
extern int PTHETA_LEN;
int RHO_VARIANCE = 2;



//*********************************************************//
// rho tracker parameters //
rtracker *g_rtrackers_ = NULL;



particle* init_distribution(double rho, double theta, int n)
{
    particle* ptcs;
    int i;

    ptcs = (particle*)calloc(n, sizeof(particle));
    for (i=0; i<n; i++)
    {
        ptcs[i].rho = rho;
        ptcs[i].rhop = rho;
        ptcs[i].theta = theta;
        ptcs[i].thetap = theta;
        ptcs[i].w = 1.0;
        ptcs[i].data = NULL;
        ptcs[i].N = 0;
    }

    return ptcs;
}

double getParticleWeight(particle *ptc, double *Bscope, double *Ang_plot,
    int M, int N, int x0, svmData *svm,
    cableLine *cable, particle *pResult)
{
    double rho = ptc->rho;
    double theta = ptc->theta;

    int lineRange[512];
    int lineAngle[512];
    double lineData[512];
    double f[FEAT_SIZE];
    double linearVal, svmVal;
    
    int dataNum = getLineData(Bscope, Ang_plot, M, N, NULL, x0,
        rho, theta, lineData, lineRange, lineAngle);
    computeFeature(lineData, dataNum, SPEC_SIZE, AUTO_SIZE, f, FEAT_SIZE);
    linearVal = linearCalc(svm, f);
    svmVal = 0 - svmCalc(svm, f);

    // get weight
    ptc->w = svmVal*linearVal;

    cable->isCable = false;
    if (svmVal>0 && linearVal>=.5)
    {
        cable->isCable = true;
    }

    if (ptc->w < 0 || svmVal < 0 || linearVal < 0.5)
        ptc->w = 0.0001;


    cable->score = linearVal;
    cable->theta = theta;
    cable->rho = rho;
    cable->nPix = dataNum;
    if (cable->data) free(cable->data);
    if (cable->range) free(cable->range);
    if (cable->angle) free(cable->angle);
    cable->data = (double*)calloc(dataNum, sizeof(double));
    cable->range = (int*)calloc(dataNum, sizeof(int));
    cable->angle = (int*)calloc(dataNum, sizeof(int));
    memcpy(cable->data, lineData, dataNum*sizeof(double));
    memcpy(cable->range, lineRange, dataNum*sizeof(int));
    memcpy(cable->angle, lineAngle, dataNum*sizeof(int));


    // update data for future weighting; get real weight
    double w1 = 0;
    int dataXor = 1;
    int smoothRhoTheta = 1;
    if (NULL == ptc->data)
    {
        ptc->N = dataNum;
        ptc->data = (double*)calloc(dataNum, sizeof(double));
        memcpy(ptc->data, lineData, dataNum*sizeof(double));
        if (pResult)
        {
            w1 = dataCorrelation(ptc->data, pResult->data, min(ptc->N, pResult->N), dataXor);
            if (w1 <= 0) w1 = 0.0001;
            ptc->w *= w1;
        }
    }
    else if (ptc->N == dataNum)
    {
        memcpy(ptc->data, lineData, dataNum*sizeof(double));
        if (pResult)
        {
            w1 = dataCorrelation(ptc->data, pResult->data, min(ptc->N, pResult->N), dataXor);
        }
        else
        {
            w1 = dataCorrelation(ptc->data, lineData, dataNum, dataXor);
        }
        if (w1 < 0) w1 = 0.0001;
        ptc->w *= w1;       
    }
    else
    {
        free(ptc->data);
        ptc->data = (double*)calloc(dataNum, sizeof(double));
        memcpy(ptc->data, lineData, dataNum*sizeof(double));
        ptc->N = dataNum;
        // different data size
        if (pResult)
        {
            w1 = dataCorrelation(ptc->data, pResult->data, min(ptc->N, pResult->N), dataXor);
        }
        else
        {
            w1 = dataCorrelation(ptc->data, lineData, min(dataNum, ptc->N), dataXor);
        }
        if (w1 < 0) w1 = 0.0001;
        ptc->w *= w1;
    }


    if (smoothRhoTheta)
    {
        double dtheta;
        double drho;
        if (pResult)
        {
            dtheta  = ptc->theta - pResult->theta;
            drho    = ptc->rho - pResult->rho;
        }
        else
        {
            dtheta  = ptc->theta - ptc->thetap;
            drho    = ptc->rho - ptc->rhop;
        }
        ptc->w *= exp(-dtheta*dtheta/(8*THETA_STD*THETA_STD));
        ptc->w *= exp(-drho*drho/(8*RHO_STD*RHO_STD));
    }

    if (abs(ptc->theta)<60 || abs(ptc->theta)>120)
        ptc->w = 0.0;
    
    return ptc->w;
}

void normalizeWeights(particle* particles, int n)
{
    double sum = 0;
    int i;

    for (i=0; i<n; i++)
    {
        sum += particles[i].w;
    }

    for (i=0; i<n; i++)
    {
        particles[i].w /= sum;
    }
}


void gauss(double std1, double std2, double *y1, double *y2)
{
    double x1, x2, w;

    do
    {
        x1 = 2.0 * rand() / (double)(RAND_MAX) - 1.0;
        x2 = 2.0 * rand() / (double)(RAND_MAX) - 1.0;
        w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = sqrt((-2.0 * log(w)) / w);
    *y1 = x1 * w * std1;
    *y2 = x2 * w * std2;
}

particle transition(particle p)
{
    particle pn;
    double theta, rho, y1, y2;

    gauss(RHO_STD, THETA_STD, &y1, &y2);
        //rho = p.rho + y1;
        //rho = (p.rho + p.rhop)/2 + y1;
    rho = p.rho + (p.rho-p.rhop) + y1;
        //theta = p.theta + y2;
        //theta = (p.theta + p.thetap)/2 + y2;
    theta = p.theta + (p.theta-p.thetap) + y2;

    pn.rho = rho;
    pn.theta = theta;
    pn.rhop = p.rho;
    pn.thetap = p.theta;
    pn.w = 0;
    pn.N = p.N;
    pn.data = p.data;

    return pn;
}

particle* resample(particle* ptcs, int n)
{
    particle *ptcsn;
    int i, j, k=0, np;

    qsort(ptcs, n, sizeof(particle), &particleCmp);
    ptcsn = (particle*)calloc(n, sizeof(particle));
    for (i=0; i<n; i++)
    {
        np = (int)(ptcs[i].w * n + 0.5);
        for (j=0; j<np; j++)
        {
            ptcsn[k] = ptcs[i];
            ptcsn[k].data = (double*)calloc(ptcs[i].N, sizeof(double));
            memcpy(ptcsn[k].data, ptcs[i].data, ptcs[i].N*sizeof(double));
            k++;
            if (k==n)
                return ptcsn;
        }
    }

    while (k<n)
    {
        ptcsn[k] = ptcs[0];
        ptcsn[k].data = (double*)calloc(ptcs[0].N, sizeof(double));
        memcpy(ptcsn[k].data, ptcs[0].data, ptcs[0].N*sizeof(double));
        k++;
    }
    return ptcsn;
}

int particleCmp(const void *p1, const void *p2)
{
    particle *_p1 = (particle*)p1;
    particle *_p2 = (particle*)p2;

    if (_p1->w > _p2->w)
        return -1;
    if (_p1->w < _p2->w)
        return 1;
    return 0;
}

void freeParticles(particle* ptcs, int n)
{
    int i;
    if (ptcs)
    {
        for (i=0; i<n; i++)
        {
            if (ptcs[i].data)
            {
                free(ptcs[i].data);
                ptcs[i].data = NULL;
                ptcs[i].N = 0;
            }
        }
    }
}

double dataCorrelation(double *data1, double *data2, int n, int dataXor)
{
    double m1, m2, s1, s2, c;
    int c1, c2;
    int i;

    m1 = m2 = 0.0;
    for (i=0; i<n; i++)
    {
        m1 += data1[i];
        m2 += data2[i];
    }
    m1 /= n;
    m2 /= n;

    if (dataXor)
    {
        c = 0;
        for (i=0; i<n; i++)
        {
            c1 = (data1[i] >= m1);
            c2 = (data2[i] >= m2);
            c += (1 - (c1 ^ c2));
        }
        c /= (double)n;
    }
    else
    {
        s1 = s2 = 0.0;
        for (i=0; i<n; i++)
        {
            s1 += (data1[i]-m1) * (data1[i]-m1);
            s2 += (data2[i]-m2) * (data2[i]-m2);
        }
        
        c = 0;
        for (i=0; i<n; i++)
        {
            c += (data1[i]-m1) * (data2[i]-m2);
        }
        c /= sqrt(s1*s2);
    }

    return c;
}



double similarityCoefficient(double *data1, double *data2, int n)
{
    double sum1, sum2, sum;
    sum1 = sum2 = sum = 0;
    int i;
    for (i=0; i<n; i++)
    {
        sum1 += data1[i];
        sum2 += data2[i];
    }
    for (i=0; i<n; i++)
    {
        sum += sqrt(abs(data1[i]*data2[i]/sum1/sum2));
    }

    return exp(-sum);
}

cableLine* sampleNewCable(cableLine *anchor, int num, double rhoStd,
    double *Bscope, double *Ang_plot,
    int M, int N, int x0, svmData *svm)
{
    cableLine *cables = (cableLine*)calloc(num*2, sizeof(cableLine));
    double anchorRho1 = anchor->rho;
    double anchorRho2 = anchor->rho;
    double rho1, rho2;
    int i;
    int dataNum;
    int lineRange[512];
    int lineAngle[512];
    double lineData[512];
    double f[FEAT_SIZE];
    double linearVal, svmVal;
    double theta = anchor->theta;

    for (i=0; i<num; i++)
    {
        gauss(rhoStd, rhoStd, &rho1, &rho2);
        rho1 = anchorRho1 + abs(rho1);

        // downward
        dataNum = getLineData(Bscope, Ang_plot, M, N, NULL, x0,
            rho1, theta, lineData, lineRange, lineAngle);
        computeFeature(lineData, dataNum, SPEC_SIZE, AUTO_SIZE, f, FEAT_SIZE);
        linearVal = linearCalc(svm, f);
        svmVal = 0 - svmCalc(svm, f);
        cables[i].isCable = (svmVal>0 && linearVal>=.5);
        cables[i].theta = theta;
        cables[i].rho = rho1;
        cables[i].nPix = dataNum;
        cables[i].data = (double*)calloc(dataNum, sizeof(double));
        cables[i].range = (int*)calloc(dataNum, sizeof(int));
        cables[i].angle = (int*)calloc(dataNum, sizeof(int));
        memcpy(cables[i].data, lineData, dataNum*sizeof(double));
        memcpy(cables[i].range, lineRange, dataNum*sizeof(int));
        memcpy(cables[i].angle, lineAngle, dataNum*sizeof(int));

        anchorRho1 = rho1;
    }

    anchorRho1 = anchorRho2;
    for (i=0; i<num; i++)
    {
        gauss(rhoStd, rhoStd, &rho1, &rho2);
        rho1 = anchorRho1 - abs(rho1);

        // upward
        dataNum = getLineData(Bscope, Ang_plot, M, N, NULL, x0,
            rho1, theta, lineData, lineRange, lineAngle);
        computeFeature(lineData, dataNum, SPEC_SIZE, AUTO_SIZE, f, FEAT_SIZE);
        linearVal = linearCalc(svm, f);
        svmVal = 0 - svmCalc(svm, f);
        cables[i+num].isCable = (svmVal>0 && linearVal>=.5);
        cables[i+num].theta = theta;
        cables[i+num].rho = rho1;
        cables[i+num].nPix = dataNum;
        cables[i+num].data = (double*)calloc(dataNum, sizeof(double));
        cables[i+num].range = (int*)calloc(dataNum, sizeof(int));
        cables[i+num].angle = (int*)calloc(dataNum, sizeof(int));
        memcpy(cables[i+num].data, lineData, dataNum*sizeof(double));
        memcpy(cables[i+num].range, lineRange, dataNum*sizeof(int));
        memcpy(cables[i+num].angle, lineAngle, dataNum*sizeof(int));

        anchorRho1 = rho1;
    }

    return cables;
}


int updateParticleTracker(particle* ptcs, int nParticle, cableLine* cables, int nCable)
{
    int *idx = (int*)calloc(nCable, sizeof(int));
    int i, k, r;

    k = 0;
    for (i=0; i<nCable; i++)
    {
        if (cables[i].isCable)
        {
            idx[k++] = i;
        }
    }

    r = k;
    if (k)
    {
        k = (int)(k/2);
        for (i=0; i<nParticle; i++)
        {
            ptcs[i].rho = cables[ idx[k] ].rho;
            ptcs[i].rhop = cables[ idx[k] ].rho;
            ptcs[i].theta = cables[idx[k]].theta;
            ptcs[i].thetap = cables[idx[k]].theta;
            ptcs[i].w = 1.0;

            if (ptcs[i].data)
            {
                free(ptcs[i].data);
                ptcs[i].data = (double*)calloc(cables[idx[k]].nPix, sizeof(double));
                memcpy(ptcs[i].data, cables[idx[k]].data, sizeof(double)*cables[idx[k]].nPix);
                ptcs[i].N = cables[idx[k]].nPix;
            }
        }
    }

    free(idx);
    return r;
}

void justifyWeights(particle* ptcs, int nParticle)
{
    double minW = 1000;
    int i;
    
    for (i=0; i<nParticle; i++)
    {
        if (ptcs[i].w < minW)
            minW = ptcs[i].w;
    }
    
    if (minW < 0)
    {
        for (i=0; i<nParticle; i++)
        {
            if (ptcs[i].w < 0)
                ptcs[i].w = 0;
            //ptcs[i].w -= minW;
        }
    }
}




// theta stuff
tparticle* theta_initialize_particles(int num, double theta, int* H, int rhoLen, int thetaLen, int idx)
{
    tparticle* theta_particles = (tparticle*)calloc(num, sizeof(tparticle));
    
    if (theta<0) theta = theta+180;
    for (int i=0; i<num; i++)
    {
        theta_particles[i].theta = theta;
        theta_particles[i].thetap = 0;
        theta_particles[i].thetapp = 0;
        theta_particles[i].w = 1;
        theta_particles[i].n = rhoLen;
        theta_particles[i].h = (int*)calloc(rhoLen, sizeof(int));
        for (int j=0; j<rhoLen; j++)
        {
            theta_particles[i].h[j] = H[j*thetaLen+idx];
        }
    }

    return theta_particles;
}

void theta_free_particles(tparticle* theta_particles, int num)
{
    for (int i=0; i<num; i++)
    {
        if (theta_particles[i].h) free(theta_particles[i].h);
        theta_particles[i].h = NULL;
    }

    free(theta_particles);
}

double theta_average(tparticle* theta_particles, int num)
{
    double theta_return = 0;
    int i;
    for (i=0; i<num; i++)
    {
        theta_return += theta_particles[i].w * theta_particles[i].theta;
    }

    return theta_return;
}

void theta_transition_particles(tparticle* theta_particles, int num, double theta_std)
{
    double tn, t, tp, tpp, dum1, dum2;
    int i;

    for (i=0; i<num; i++)
    {
        t   = theta_particles[i].theta;
        tp  = theta_particles[i].thetap;
        tpp = theta_particles[i].thetapp;
        if (tp>0 && tpp>0)
            tn = 3*t - 3*tp + tpp;
        else if (tp>0 && tpp<=0)
            tn = 2*t - tp;
        else
            tn = t;
        gauss(theta_std, theta_std, &dum1, &dum2);
        tn += dum1;

        if (tn<PTHETA_END) tn = PTHETA_END;
        if (tn>180-PTHETA_END) tn = 180-PTHETA_END;

        theta_particles[i].thetapp = tp;
        theta_particles[i].thetap  = t;
        theta_particles[i].theta   = tn;
    }
}

void theta_compute_likelihood(tparticle* theta_particles, int num, int* H, int rhoLen, int thetaLen,
    double(*similarity)(int *h1, int *h2, int n))
{
    int i, j, *h, idx;
    int thetaI;
    double w;

    h = (int*)calloc(rhoLen, sizeof(int));
    for (i=0; i<num; i++)
    {
        thetaI = (int)(theta_particles[i].theta + 0.5);
        if (thetaI<PTHETA_END) thetaI = PTHETA_END;
        if (thetaI>180-PTHETA_END) thetaI = 180-PTHETA_END;
        if (thetaI>=90) idx = thetaI-90;
        else idx = thetaI-(90-PTHETA_LEN);

        for (j=0; j<rhoLen; j++)
        {
            h[j] = H[j*thetaLen+idx];
        }

        w = (*similarity)(theta_particles[i].h, h, rhoLen);
        w *= entropy_max3(h, rhoLen, 1);
        w *= exp(-2*(theta_particles[i].theta-theta_previous)*(theta_particles[i].theta-theta_previous)/(100*0.5*0.5));
        theta_particles[i].w = w;

        for (j=0; j<rhoLen; j++)
        {
            theta_particles[i].h[j] = H[j*thetaLen+idx];
        }
    }

    if (h) free(h);
}

double similarity_1(int *h1, int *h2, int n)
{
    return 1;
}

double similarity_cos(int *h1, int *h2, int n)
{
    double num = 0;
    double dem1 = 0;
    double dem2 = 0;
    int i;

    for (i=0; i<n; i++)
    {
        num  += h1[i]*h2[i];
        dem1 += h1[i]*h1[i];
        dem2 += h2[i]*h2[i];
    }

    return num/sqrt(dem1*dem2);
}

double entropy_max3(int *h, int rhoLen, int ss)
{
    double w = 0;
    int j, k;
    int localMaxSum, localMax, localIdx;
    
    for (j=0; j<rhoLen; j++)
    {
        if (h[j*ss])
            w += h[j*ss] * log((double)h[j*ss]);
    }

    localMaxSum = 0;
    for (k=0; k<3; k++)
    {
        localMax = -1;
        localIdx = -1;
        for (j=0; j<rhoLen; j++)
        {
            if (h[j*ss]>localMax)
            {
                localMax = h[j*ss];
                localIdx = j;
            }
        }
        localMaxSum += localMax;
        for (j=localIdx-50; j<=localIdx+50; j++)
        {
            if (j>=0&& j<rhoLen)
                h[j*ss] = 0;
        }
    }

    w *= (double)localMaxSum;
    return w;
}

void theta_normalize_weights(tparticle *theta_particles, int num)
{
    double wsum = 0;
    int i;

    for (i=0; i<num; i++)
    {
        wsum += theta_particles[i].w;
    }

    for (i=0; i<num; i++)
    {
        theta_particles[i].w /= wsum;
    }
}

int tparticleCmp(const void *p1, const void *p2)
{
    tparticle *_p1 = (tparticle*)p1;
    tparticle *_p2 = (tparticle*)p2;

    if (_p1->w > _p2->w)
        return -1;
    if (_p1->w < _p2->w)
        return 1;
    return 0;
}


tparticle* tresample(tparticle* ptcs, int n)
{
    tparticle *ptcsn;
    int i, j, k=0, np;

    qsort(ptcs, n, sizeof(tparticle), &tparticleCmp);
    ptcsn = (tparticle*)calloc(n, sizeof(tparticle));
    for (i=0; i<n; i++)
    {
        np = (int)(ptcs[i].w * n + 0.5);
        for (j=0; j<np; j++)
        {
            ptcsn[k] = ptcs[i];
            ptcsn[k].h = (int*)calloc(ptcs[i].n, sizeof(int));
            memcpy(ptcsn[k].h, ptcs[i].h, ptcs[i].n*sizeof(int));
            k++;
            if (k==n)
                return ptcsn;
        }
    }

    while (k<n)
    {
        ptcsn[k] = ptcs[0];
        ptcsn[k].h = (int*)calloc(ptcs[0].n, sizeof(int));
        memcpy(ptcsn[k].h, ptcs[0].h, ptcs[0].n*sizeof(int));
        k++;
    }
    return ptcsn;
}




int add_rtrackers(rtracker **rtrackers__, int *h, int rhoLen, double theta, double *rho, double *RP_dBm, int *BW, double *Ang_plot, int M, int N, int x0, svmData *svm, int margin, double thresh)
{
    rtracker *current;
    double vMax, thetaN, linearVal, svmVal;
    int iMax, i, j, dataNum, num_added;
    int lineRange[512];
    int lineAngle[512];
    double lineData[512];
    double f[FEAT_SIZE];

    if (!(*rtrackers__))
    {
        *rtrackers__ = rtracker_init_zero();
    }
    
    // get to the place
    current = *rtrackers__;
    while (current->next_) current = current->next_;
    num_added=0;

    vMax = -1;
    iMax = -1;
    for (i=0; i<rhoLen; i++)
    {
        if (vMax < h[i])
        {
            vMax = h[i];
            iMax = i;
        }
    }


    thetaN = theta;
    if (thetaN>=90) thetaN -= 180;
    // get rho candidates:1. local maxima; 2. over threshold; 3. pass svm and linear
    while (vMax>=thresh)
    {       
        dataNum = getLineData(RP_dBm, Ang_plot, M, N,
            BW, x0, rho[iMax], thetaN,
            lineData, lineRange, lineAngle);
        computeFeature(lineData, dataNum, SPEC_SIZE, AUTO_SIZE, f, FEAT_SIZE);
        linearVal = linearCalc(svm, f);
        svmVal = svmCalc(svm, f);
        if (dataNum>=15 && (svmVal<0 && linearVal>=0.5))
        {
            num_added++;
            dataNum = getLineData(RP_dBm, Ang_plot, M, N,
                NULL, x0, rho[iMax], thetaN,
                lineData, lineRange, lineAngle);
            current->next_ = rtracker_init(rho[iMax], theta,
                lineData, lineRange, lineAngle, dataNum,
                BW, M, N, vMax, svmVal, linearVal, h+iMax-HOUGH_SIDE, 2*HOUGH_SIDE+1);
            current = current->next_;
        }

        // get next one
        for (j=iMax-margin; j<=iMax+margin; j++)
        {
            if (j>=0 && j<rhoLen)
                h[j] = 0;
        }
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
    }

    return num_added;
}


rtracker *rtracker_init_zero()
{
    rtracker *rtracker_ = (rtracker*)calloc(1, sizeof(rtracker));
    rtracker_->data_num = 0;
    rtracker_->data_original_ = NULL;
    rtracker_->data_threshold_ = NULL;
    rtracker_->life_time = 0;
    rtracker_->rho = 0;
    rtracker_->rhop = 0;
    rtracker_->rhopp = 0;
    rtracker_->rpts_num = 0;
    rtracker_->rpts_ = NULL;
    rtracker_->rstatus = head;
    rtracker_->status_time = 0;
    rtracker_->next_ = NULL;


    return rtracker_;
}


rparticle* rho_initialize_particles(int num, double rho, double theta, double hVal, double svmVal, double linearVal)
{
    rparticle *rpts;
    int i;

    rpts = (rparticle*)calloc(num, sizeof(rparticle));
    for (i=0; i<num; i++)
    {
        rpts[i].rho     = rho;//>0 ? rho : -rho;
        rpts[i].rhop    = 0;
        rpts[i].rhopp   = 0;
        rpts[i].theta   = theta;
        rpts[i].thetap  = 0;
        rpts[i].thetapp = 0;
        rpts[i].w       = 1;

        rpts[i].hVal = hVal;
        rpts[i].svmVal = svmVal;
        rpts[i].linearVal = linearVal;
    }

    return rpts;
}


rtracker *rtracker_init(double rho, double theta,
    double *lineData, int *lineRange, int *lineAngle, int dataNum,
    int *BW, int M, int N, double hVal, double svmVal, double linearVal,
    int *hough_data, int hnum)
{
    int i;
    rtracker *rtracker_ = (rtracker*)calloc(1, sizeof(rtracker));
    
    // initialize particles code
    rtracker_->rpts_ = rho_initialize_particles(RHO_NUM_PARTICLES, rho, theta, hVal, svmVal, linearVal);
    rtracker_->rpts_num = RHO_NUM_PARTICLES;
    
    // other stuff
    rtracker_->rstatus = candidate;
    rtracker_->life_time = 1;
    rtracker_->status_time = 1;

    rtracker_->rho = rho;//>0 ? rho : -rho;
    rtracker_->rhop = 0;
    rtracker_->rhopp = 0;
    rtracker_->theta = theta;
    rtracker_->thetap = 0;
    rtracker_->thetapp = 0;

    rtracker_->data_num = dataNum;
    rtracker_->data_original_ = (double*)calloc(dataNum, sizeof(double));
    rtracker_->data_threshold_ = (int*)calloc(dataNum, sizeof(int));
    memcpy(rtracker_->data_original_, lineData, dataNum*sizeof(double));
    for (i=0; i<dataNum; i++)
    {
        rtracker_->data_threshold_[i] = BW[ lineRange[i]*N+lineAngle[i] ];
    }

    // hough data storage
    rtracker_->hough_ = (int*)calloc(hnum, sizeof(int));
    memcpy(rtracker_->hough_, hough_data, hnum*sizeof(int));
    rtracker_->hough_num = hnum;

    rtracker_->next_ = NULL;

    return rtracker_;
}

void rtracker_free(rtracker *rtracker_)
{
    if (rtracker_->data_original_) free(rtracker_->data_original_);
    if (rtracker_->data_threshold_) free(rtracker_->data_threshold_);
    if (rtracker_->rpts_) free(rtracker_->rpts_);
    if (rtracker_->hough_) free(rtracker_->hough_);
    
    rtracker_->data_original_ = NULL;
    rtracker_->data_threshold_ = NULL;
    rtracker_->rpts_ = NULL;
    rtracker_->hough_ = NULL;
    rtracker_->rpts_num = 0;
    rtracker_->hough_num = 0;
    rtracker_->rho = 0;
    rtracker_->rhop = 0;
    rtracker_->rhopp = 0;
    rtracker_->theta = 0;
    rtracker_->thetap = 0;
    rtracker_->thetapp = 0;
}

int rtracker_free_list(rtracker *head)
{
    int num = 0;
    rtracker *current, *tmp;
    current = head;
    while (current)
    {
        rtracker_free(current);
        tmp = current->next_;
        free(current);
        current = tmp;
        num++;
    }
    return num;
}

int get_rtracker_num(rtracker *rtracker_)
{
    int num = 0;
    while (rtracker_)
    {
        if (abs(rtracker_->rho)>0.01)
            num++;
        rtracker_ = rtracker_->next_;
    }
    return num;
}

int get_rtracker_rho(rtracker *rtracker_, double *rho, int *is_cable)
{
    int num = 0;
    while (rtracker_)
    {
        if (abs(rtracker_->rho)>0.01)
        {           
            *rho = rtracker_->rho;
            if (rtracker_->rstatus == healthy)
                *is_cable = 1;
            else if (rtracker_->rstatus == dormant)
                *is_cable = 2;
            else
                *is_cable = 0;
            num++;
            rho++;
            is_cable++;
        }
        rtracker_ = rtracker_->next_;
    }
    return num;
}

void rho_transition_particles(rparticle *rho_particles, int num, double rho_std, double theta_new, int ss, rtracker *r_track)
{
    double rn, r, rp, rpp, dum1, dum2;
    int i;
    bool turn_over = false;

    // condition of turning over 90 degree line
    if ((r_track->theta >= 90 && theta_new < 90)
        || (r_track->theta < 90 && theta_new >= 90))
        turn_over = true;

    if (r_track)
    {
        r   = r_track->rho;
        rp  = r_track->rhop;
        rpp = r_track->rhopp;
        /*
        if (abs(rp)>.1 && abs(rpp)>.1)
            //rn = 3*r - 3*rp + rpp;
            rn = 2*r - rp;
        else if(abs(rp)>.1 && abs(rpp) <.1)
            rn = 2*r - rp;
        else
            rn = r;
            */
    
        if (turn_over)
            rn = -r;
        else
            rn = r;
        
        for (i=0; i<num; i++)
        {
            rho_particles[i].rhopp = rp;
            rho_particles[i].rhop  = r;
            rho_particles[i].rho   = rn + (i-num/2)*ss;

            rho_particles[i].thetapp= rho_particles[i].thetap;
            rho_particles[i].thetap = rho_particles[i].theta;
            rho_particles[i].theta  = theta_new;
        }
    }
    else
    {
        for (i=0; i<num; i++)
        {
            r   = rho_particles[i].rho;
            rp  = rho_particles[i].rhop;
            rpp = rho_particles[i].rhopp;
            /*
            if (abs(rp)>.1 && abs(rpp)>.1)
                //rn = 3*r - 3*rp + rpp;
                rn = 2*r - rp;
            else if(abs(rp)>.1 && abs(rpp) <.1)
                rn = 2*r - rp;
            else
                rn = r;
                */

            if (turn_over)
                rn = -r;
            else
                rn = r;

            gauss(rho_std, rho_std, &dum1, &dum2);
            rn += dum1;

            rho_particles[i].rhopp  = rp;
            rho_particles[i].rhop   = r;
            rho_particles[i].rho    = (int)(rn+0.5);

            rho_particles[i].thetapp= rho_particles[i].thetap;
            rho_particles[i].thetap = rho_particles[i].theta;
            rho_particles[i].theta  = theta_new;
        }
    }
}

int rho_to_idx(double rho, double *rhos, int rhoLen)
{
    double idxD;
    int idx;
    idxD = (rho-rhos[0])/(rhos[1]-rhos[0]);
    idx = (int)(idxD+0.5);
    return idx;
}

double rho_compute_likelihood(rparticle *rho_particles, int num, double rhop,
    double theta, int *h, double *rho, int rhoLen, int thetaLen,
    double *RP_dBm, int *BW, int *BWcart, double *Ang_plot,
    int M, int N, int NC, int x0, svmData *svm, rtracker *r_track)
{
    int i, ridx, j;
    int dataNum;
    double w;
    int lineRange[512];
    int lineAngle[512];
    double linearVal, svmVal;
    double lineData[512];
    double f[FEAT_SIZE];
    double base = 0;

    double data1[512];
    double data2[512];

    double asso_s = -1;
    double rhod;

    for (j=0; j<2*HOUGH_SIDE+1; j++)
    {
        base += r_track->hough_[j] * r_track->hough_[j];
    }

    for (i=0; i<num; i++)
    {
        ridx = rho_to_idx(rho_particles[i].rho, rho, rhoLen);
        w = 1;
        
        // current SVM
        dataNum = getLineData(RP_dBm, Ang_plot, M, N,
            BW, x0, rho_particles[i].rho, rho_particles[i].theta<90 ? rho_particles[i].theta : rho_particles[i].theta-180,
            lineData, lineRange, lineAngle);
        computeFeature(lineData, dataNum, SPEC_SIZE, AUTO_SIZE, f, FEAT_SIZE);
        linearVal = linearCalc(svm, f);
        svmVal = svmCalc(svm, f);
        if (dataNum<15 || linearVal<0) linearVal = 0;
        if (svmVal > SVM_BASE) svmVal = SVM_BASE;

        if (HOUGH_SIDE)
        {
            w = 0;
            for (j=-HOUGH_SIDE; j<=HOUGH_SIDE; j++)
            {
                if (ridx+j>=0 && ridx+j<rhoLen)
                {
                    w += r_track->hough_[j+HOUGH_SIDE] * h[ridx+j];
                }
            }

            w /= base;
            rhod = abs(rho_particles[i].rho) - abs(rhop);
            w *= exp(-2*rhod*rhod/(400*RHO_VARIANCE*RHO_VARIANCE));
            w *= linearVal * (-svmVal + SVM_BASE);
        }
        else
        {
            w = 1;
            w *= linearVal;
            w *= (-svmVal + SVM_BASE);

            // similarity in rho
            rhod = abs(rho_particles[i].rho) - abs(rhop);
            w *= exp(-2*rhod*rhod/(400*RHO_VARIANCE*RHO_VARIANCE));
            
            // original data correlation
            dataNum = getLineData(RP_dBm, Ang_plot, M, N,
                NULL, x0, rho_particles[i].rho, rho_particles[i].theta<90 ? rho_particles[i].theta : rho_particles[i].theta-180,
                lineData, lineRange, lineAngle);
            // TODO: why this is not used?
            //double dc = dataCorrelation(lineData, r_track->data_original_, dataNum);

            memcpy(data1, lineData, dataNum*sizeof(double));
            memcpy(data2, r_track->data_original_, dataNum*sizeof(double));
            double ds = simple_c(data1, data2, dataNum);
            double bc = bhattacharyya(data1, data2, dataNum);
            w *= (ds*bc);
        }

        if (w > asso_s)
            asso_s = w;

        // update
        rho_particles[i].w = w;
        rho_particles[i].hVal = h[ridx];
        //rho_particles[i].svmVal = svmVal;
        //rho_particles[i].linearVal = linearVal;
    }

    return asso_s;
}


void rho_normalize_weights(rparticle *rpts, int num)
{
    double sum = 0;
    int i;
    for (i=0; i<num; i++)
        sum += rpts[i].w;
    for (i=0; i<num; i++)
        rpts[i].w /= sum;
}


int rparticleCmp(const void *p1, const void *p2)
{
    rparticle *_p1 = (rparticle*)p1;
    rparticle *_p2 = (rparticle*)p2;

    if (_p1->w > _p2->w)
        return -1;
    if (_p1->w < _p2->w)
        return 1;
    return 0;
}


rparticle* rresample(rparticle* ptcs, int n)
{
    rparticle *ptcsn;
    int i, j, k=0, np;

    qsort(ptcs, n, sizeof(rparticle), &rparticleCmp);
    ptcsn = (rparticle*)calloc(n, sizeof(rparticle));
    for (i=0; i<n; i++)
    {
        np = (int)(ptcs[i].w * n + 0.5);
        for (j=0; j<np; j++)
        {
            ptcsn[k] = ptcs[i];
            k++;
            if (k==n)
                return ptcsn;
        }
    }

    while (k<n)
    {
        ptcsn[k] = ptcs[0];
        k++;
    }
    return ptcsn;
}



void process_rtrackers(rtracker *rtracker_,
    double theta, int *h, double *rho, int rhoLen, int thetaLen,
    double *RP_dBm, int *BW, int *BWcart, double *Ang_plot,
    int M, int N, int NC, int x0, svmData *svm, int margin)
{
    rtracker *current = rtracker_;
    rtracker *previous;
    rparticle *rptsn;
    double asso_s;  // associatio score
    int dataNum, i, j, ridx;
    double lineData[512];
    int lineRange[512];
    int lineAngle[512];
    bool del = false;

    int ss = 1;
    ss = rhoLen>6000 ? 2 : 1;
    if (rhoLen>6000) RHO_VARIANCE = 4;

    while (current)
    {
        del = false;
        if (current->rstatus == candidate || current->rstatus == healthy)
        {
            // run the particle filtering algorithm
            rho_transition_particles(current->rpts_, current->rpts_num, RHO_VARIANCE, theta, ss, current);
            asso_s = rho_compute_likelihood(current->rpts_, current->rpts_num, current->rho,
                theta, h, rho, rhoLen, thetaLen,
                RP_dBm, BW, BWcart, Ang_plot,
                M, N, NC, x0, svm, current);
            rho_normalize_weights(current->rpts_, current->rpts_num);
            rptsn = rresample(current->rpts_, current->rpts_num);

            free(current->rpts_);
            current->rpts_ = rptsn;

            // update tracker
            current->rhopp = current->rhop;
            current->rhop  = current->rho;
            current->rho   = current->rpts_[0].rho;
            
            current->thetapp = current->thetap;
            current->thetap  = current->theta;
            current->theta   = theta;

            // update hough datas
            int hidx = rho_to_idx(current->rho, rho, rhoLen);
            memcpy(current->hough_, h+hidx-HOUGH_SIDE, (2*HOUGH_SIDE+1)*sizeof(int));
            
            dataNum = getLineData(RP_dBm, Ang_plot, M, N,
                0, x0, current->rho, current->theta<90 ? current->theta : current->theta-180,
                lineData, lineRange, lineAngle);
            //computeFeature(lineData, dataNum, SPEC_SIZE, AUTO_SIZE, f, FEAT_SIZE);
            current->data_num = dataNum;

            free(current->data_original_);
            free(current->data_threshold_);
            current->data_original_ = (double*)calloc(dataNum, sizeof(double));
            current->data_threshold_ = (int*)calloc(dataNum, sizeof(double));
            memcpy(current->data_original_, lineData, sizeof(double)*dataNum);
            for (i=0; i<dataNum; i++)
            {
                current->data_threshold_[i] = BW[ lineRange[i]*N+lineAngle[i] ];
            }
            
            // output results
            if (asso_s >= ASSOCIATION_THRESHOLD)
            {
                current->rstatus = healthy;
                current->life_time++;
                current->status_time++;
                // suppress h
                ridx = rho_to_idx(current->rho, rho, rhoLen);
                for (j=ridx-margin; j<=ridx+margin; j++)
                {
                    if (j>=0 && j<rhoLen)
                        h[j] = 0;
                }
            }
            else
            {
                
                if (current->rstatus == healthy)
                {
                    current->rstatus = dormant;
                    current->life_time++;
                    current->status_time = 1;
                }
                else
                {
                    // candidate with no association: kill it
                    previous->next_ = current->next_;
                    rtracker_free(current);
                    del = true;
                }
            }
        }
        else if (current->rstatus == dormant)
        {
            if (current->status_time >= DORMANCY_THRESHOLD)
            {
                previous->next_ = current->next_;
                rtracker_free(current);
                del = true;
            }
        }

        if (!del) previous = current;
        current = current->next_;
    }
}

double get_rho_thresh(int *h, int rhoLen)
{
    double vMax, thresh;
    int iMax, i, iDum;

    vMax = -1;
    iMax = -1;
    thresh = 0;
    for (i=0; i<rhoLen; i++)
    {
        thresh += h[i];
        if (h[i]>vMax)
        {
            vMax = h[i];
            iMax = i;
        }
    }

    int *h2 = (int*)calloc(rhoLen, sizeof(int));
    memcpy(h2, h, rhoLen*sizeof(int));

    for (int k=0; k<0; k++)
    {
        for (i=iMax-15; i<=iMax+15; i++)
        {
            if (i>=0 && i<rhoLen)
                h2[i] = 0;
        }
        
        getMaxVal(h2, rhoLen, 1, &vMax, &iMax, &iDum);
    }
    
    thresh /= rhoLen;
    thresh += vMax;
    thresh /= 3;

    return thresh;
}

double bhattacharyya(double *data1, double *data2, int n)
{
    double m1, m2, s1, s2, c, max1, max2;
    int i;

    // minimum
    m1 = data1[0];
    m2 = data2[0];
    for (i=0; i<n; i++)
    {
        if (data1[i]<m1)
            m1 = data1[i];
        if (data2[i]<m2)
            m2 = data2[i];
    }

    // offset
    s1 = 0;
    s2 = 0;
    max1 = -1;
    max2 = -1;
    for (i=0; i<n; i++)
    {
        data1[i] -= m1;
        data2[i] -= m2;
        if (data1[i]>max1)
            max1 = data1[i];
        if (data2[i]>max2)
            max2 = data2[i];
        //s1 += data1[i];
        //s2 += data2[i];
    }

    // normalize
    for (i=0; i<n; i++)
    {
        data1[i] /= max1;
        data2[i] /= max2;
    }

    // coefficient
    c = 0;
    s1 = 0;
    s2 = 0;
    for (i=0; i<n; i++)
    {
        c += sqrt(data1[i] * data2[i]);
        //s1 += data1[i] * data1[i];
        //s2 += data2[i] * data2[i];
    }
    return c;
}

double simple_c(double *data1, double *data2, int n)
{
    double m1, m2, M1, M2, a1, a2, dm, dM, da;
    int i;

    m1 = data1[0];
    m2 = data2[0];
    M1 = data1[0];
    M2 = data2[0];
    a1 = 0;
    a2 = 0;
    for (i=0; i<n; i++)
    {
        a1 += data1[i];
        a2 += data2[i];

        if (data1[i]<m1)
            m1 = data1[i];
        if (data1[i]>M1)
            M1 = data1[i];

        if (data2[i]<m2)
            m2 = data2[i];
        if (data2[i]>M2)
            M2 = data2[i];
    }

    a1 /= n;
    a2 /= n;

    dm = exp(-(m1-m2)*(m1-m2)/100);
    dM = exp(-(M1-M2)*(M1-M2)/100);
    da = exp(-(a1-a2)*(a1-a2)/100);

    return dm*dM*da;
}
