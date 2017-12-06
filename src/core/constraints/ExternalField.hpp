#ifndef CONSTRAINTS_EXTERNALFIELD_HPP
#define CONSTRAINTS_EXTERNALFIELD_HPP

#include "Constraint.hpp"
#include "particle_data.hpp"
#include "LatticeMultiArray.hpp"

namespace Constraints {

namespace Coupling {
 class Direct {
   Vector3d operator()(Particle const&, Vector3d const& field) const {
    return field;
   }
 };

 class Viscous {
   Vector3d operator()(Particle const& p, Vector3d const& field) const {
    return (field - Vector3d{p.m.v});
   }
 };

 class Charge {
   Vector3d operator()(Particle const& p, Vector3d const& field) const {
    return p.p.q * field;
   }
 };
}

//template<typename Coupling>
class ExternalField : public Constraint {
  public:
	//ExternalField(typename LatticeMA<Vector3d,3>::tArrayRef data, std::array<int, 3> halo_size, Coupling = Coupling{})
	ExternalField(typename LatticeMA<Vector3d,3>::tArrayRef data, std::array<int, 3> halo_size)
		: m_lattice(data, halo_size)
	{
		//m_lattice.calculate_gradient();
	}
 
    void set_weights(std::unordered_map<int, double> weights) {
        m_weights = weights;
    }

  virtual void add_energy(Particle *p, double *folded_pos,
      Observable_stat &energy) const override;
 
  virtual void add_force(Particle *p, double *folded_pos) override;

  private:
  	LatticeMA<Vector3d, 3> m_lattice;
    std::unordered_map<int, double> m_weights;
    //Coupling m_coupling;
};

} /* namespace Constraints */

#endif
                                                                                                                                                                                
                                                                                                                                                                                


