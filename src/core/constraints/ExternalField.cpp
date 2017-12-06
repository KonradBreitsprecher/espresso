#include "ExternalField.hpp"
#include "energy_inline.hpp"

namespace Constraints {

void ExternalField::add_force(Particle *p, double *folded_pos) {

	boost::array<T,N> grad;
	if (m_lattice.get_gradient_data_from_position(pos, grad))
	{
		for (int i=0; i<3; ++i) {
			p->f.f[i] += grad[i];
		}
	}
	else
		std::cout << "particle not on node" << std::endl;

}

void ExternalField::add_energy(Particle *p, double *folded_pos, Observable_stat &energy) const {
	T pot;
	if (m_lattice.get_lattice_data_from_position(pos, pot))
    	energy.non_bonded[0] += pot;
	else
		std::cout << "particle not on node" << std::endl;
	
}

}
