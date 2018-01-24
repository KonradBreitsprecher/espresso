#include "ExternalField.hpp"
#include "energy_inline.hpp"

namespace Constraints {

void ExternalField::add_force(Particle *p, double *folded_pos) {

	Vector3d force;
	if (m_lattice.get_lattice_data_from_position(folded_pos, force))
	{
		for (int i=0; i<3; ++i) {
			p->f.f[i] += force[i] * m_weights[p->p.identity];
		}
	}
	else
		std::cout << "node " << this_node << " node grid " << node_grid[0] << " " << node_grid[1] << " " << node_grid[2] <<  "  ERROR: cannot map particle " << p->p.identity << " at " << folded_pos[0] << " " <<  folded_pos[1] << " " <<   folded_pos[2] << " on local data" << std::endl;

}

void ExternalField::add_energy(Particle *p, double *folded_pos, Observable_stat &energy) const {
    /*
	T pot;
	if (m_lattice.get_lattice_data_from_position(pos, pot))
    	energy.non_bonded[0] += pot;
	else
		std::cout << "particle not on node" << std::endl;
    */
	
}

}
