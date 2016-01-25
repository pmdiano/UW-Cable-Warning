/**
 * Particle Filter class definition and implementation
 */
#ifndef __PFILTER_H__
#define __PFILTER_H__

template<typename P, typename Target, typename Observation>
class ParticleFilter {
public:
  ParticleFilter(uint32_t num_particles, Target initial_target) :
    m_particles(num_particles, nullptr)
  {
    double weight = 1.0 / (double)num_particles;
    for (int i = 0; i < m_particles.size(); ++i) {
      m_particles[i] = std::make_shared<P>(P(weight, initial_target));
    }
  }

  virtual ~ParticleFilter() {}

  Target track(Observation o)
  {
    for (auto &particle : m_particles) {
      particle->transition();
      particle->observe(o);
      particle->setWeight(particle->measure());
    }

    normalize_weights();
    resample_particles();

    return m_particles[0]->getTarget();
  }

private:
  void normalize_weights(void) {
    double wsum = 0;

    for (auto &particle : m_particles) {
      wsum += particle->getWeight();
    }

    for (auto &particle : m_particles) {
      particle->setWeight(particle->getWeight() / wsum);
    }
  }

  void resample_particles(void) {
    sort(m_particles.begin(), m_particles.end(),
      [](const std::shared_ptr<Particle<Target, Observation>> &lhs,
         const std::shared_ptr<Particle<Target, Observation>> &rhs) {
        return lhs->getWeight() > rhs->getWeight();
      });

    auto replica = m_particles;
    int k = 0, n = m_particles.size();

    for (int i = 0; i < n && k < n; ++i) {
      int np = (int)(replica[i]->getWeight() * n + 0.5);
      while (k < n && np--) {
        m_particles[k++] = replica[i];
      }
    }

    while (k < n) {
      m_particles[k++] = replica[0];
    }
  }

  std::vector<std::shared_ptr<Particle<Target, Observation>>> m_particles;
};

#endif
