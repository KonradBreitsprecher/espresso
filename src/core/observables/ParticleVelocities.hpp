#ifndef OBSERVABLES_PARTICLEVELOCITIES_HPP
#define OBSERVABLES_PARTICLEVELOCITIES_HPP

#include "PidObservable.hpp"

#include <vector>

namespace Observables {

class ParticleVelocities : public PidObservable {
public:
  virtual std::vector<double> operator()(PartCfg &partCfg) const override {
    std::vector<double> res(n_values());
    for (int i = 0; i < ids().size(); i++) {
      res[3 * i + 0] = partCfg[ids()[i]].m.v[0] / time_step;
      res[3 * i + 1] = partCfg[ids()[i]].m.v[1] / time_step;
      res[3 * i + 2] = partCfg[ids()[i]].m.v[2] / time_step;
    }
    return res;
  };
  virtual int n_values() const override { return 3 * ids().size(); }
};

} // Namespace Observables
#endif
