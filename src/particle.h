/**
 * Particle base class definition
 */
#ifndef __PARTICLE_H__
#define __PARTICLE_H__

template<typename Target, typename Observation>
class Particle {
protected:
  double m_weight;
  Target m_target;

public:
  Particle(double weight, Target target) :
    m_weight(weight),
    m_target(target)
  {
  }
  virtual ~Particle() {}

  double getWeight(void) { return m_weight; }
  void setWeight(double w) { m_weight = w; }
  Target getTarget(void) { return m_target; }

  virtual void transition(void) = 0;
  virtual void observe(Observation o) = 0;
  virtual double measure(void) = 0;
};

#endif
