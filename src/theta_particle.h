/**
 * Theta particle class definition
 */
#ifndef __THETA_PARTICLE_H__
#define __THETA_PARTICLE_H__

struct ThetaObservation {
  int *hough;
  int rhoLen;
  int thetaLen;
  ThetaObservation() : hough(nullptr), rhoLen(0), thetaLen(0) {}
  ThetaObservation(int *h, int rlen, int tlen) : hough(h), rhoLen(rlen), thetaLen(tlen) {}
};

class ThetaParticle : public Particle<double, ThetaObservation> {
  constexpr static const double unitialized_theta = -1.0;
  constexpr static const double theta_std = 1.0;
  constexpr static const double epsilon = 0.0001;

  double m_prev;
  double m_prevprev;
  ThetaObservation m_observation;

  int theta_to_index(double theta);
  double entropy_max3(std::vector<int> &h);

public:
  ThetaParticle(double weight, double target) :
    Particle<double, ThetaObservation>(weight, target),
    m_prev(unitialized_theta),
    m_prevprev(unitialized_theta)
  {}

  virtual ~ThetaParticle() {}

  virtual void observe(ThetaObservation o) override {
    m_observation = o;
  }

  virtual void transition(void) override;
  virtual double measure(void) override;
};

#endif
