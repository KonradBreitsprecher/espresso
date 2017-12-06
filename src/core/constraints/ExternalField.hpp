#ifndef CONSTRAINTS_EXTERNALFIELD_HPP
#define CONSTRAINTS_EXTERNALFIELD_HPP

#include "Constraint.hpp"
#include "particle_data.hpp"
#include "LatticeMultiArray.hpp"

namespace Constraints {

template <typename T, size_t N>
class ExternalField ,: public Constraint {
  public:
	ExternalField(typename Lattice<T,N>::tArrayRef data, std::array<int, N> halo_size)
		: m_lattice(data, halo_size)
	{
		m_lattice.calculate_gradient();
	}
 
 
  void set_field() {
  }

  virtual void add_energy(Particle *p, double *folded_pos,
      Observable_stat &energy) const override;
 
  virtual void add_force(Particle *p, double *folded_pos) override;

  private:
  	Lattice<T,N> m_lattice;

};

} /* namespace Constraints */

#endif
                                                                                                                                                                                
                                                                                                                                                                                


