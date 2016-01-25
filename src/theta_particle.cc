/**
 * Theta particle class implementation
 */
#include "inc.h"
using namespace std;

extern int PTHETA_END;
extern int PTHETA_LEN;
extern void gauss(double std1, double std2, double *y1, double *y2);

int ThetaParticle::theta_to_index(double theta) {
  int thetaI = (int)(theta + 0.5), index = 0;

  if (thetaI<PTHETA_END) thetaI = PTHETA_END;
  if (thetaI>180-PTHETA_END) thetaI = 180-PTHETA_END;

  if (thetaI>=90) index = thetaI-90;
  else index = thetaI-(90-PTHETA_LEN);

  return index;
}

void ThetaParticle::transition(void) {
  double target_new, dummy1, dummy2;
  bool prev = fabs(m_prev - unitialized_theta) <= epsilon;
  bool prevprev = fabs(m_prevprev - unitialized_theta) <= epsilon;

  if (prevprev)
    target_new = 3 * m_target - 3 * m_prev + m_prevprev;
  else if (prev)
    target_new = 2 * m_target - m_prev;
  else
    target_new = m_target;

  gauss(theta_std, theta_std, &dummy1, &dummy2);
  target_new += dummy1;

  if (target_new < PTHETA_END) target_new = PTHETA_END;
  if (target_new > 180 - PTHETA_END) target_new = 180 - PTHETA_END;

  m_prevprev = m_prev;
  m_prev = m_target;
  m_target = target_new;
}

double ThetaParticle::measure(void) {
  vector<int> h(m_observation.rhoLen);
  int theta_index = theta_to_index(m_target);

  for (int i = 0; i < m_observation.rhoLen; ++i) {
    h[i] = m_observation.hough[i * m_observation.thetaLen + theta_index];
  }

  double weight = entropy_max3(h);
  double delta = m_target - m_prev;
  weight *= exp(-2*delta*delta/(100*0.5*0.5));

  return weight;
}

double ThetaParticle::entropy_max3(vector<int> &h)
{
  double w = 0;
  int j, k;
  int localMaxSum, localMax, localIdx;
  
  for (j=0; j<m_observation.rhoLen; j++) {
    if (h[j]) {
      w += h[j] * log((double)h[j]);
    }
  }

  localMaxSum = 0;
  for (k=0; k<3; k++) {
    localMax = -1;
    localIdx = -1;

    for (j=0; j<m_observation.rhoLen; j++) {
      if (h[j]>localMax) {
        localMax = h[j];
        localIdx = j;
      }
    }

    localMaxSum += localMax;

    for (j=localIdx-50; j<=localIdx+50; j++) {
      if (j>=0&& j<m_observation.rhoLen) {
        h[j] = 0;
      }
    }
  }

  w *= (double)localMaxSum;
  return w;
}
