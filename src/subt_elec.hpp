/*
  Copyright (C) 2010,2012 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
#ifndef SUBT_ELEC_H
#define SUBT_ELEC_H
/** \file subt_elec.hpp
 *  Routines to subtract the electrostatic Energy and/or the electrostatic force 
 *  for a particle pair.
 *  \ref forces.cpp
*/
#include "utils.hpp"
#include "interaction_data.hpp"
#include "p3m.hpp"
#include "grid.hpp"

#ifdef ELECTROSTATICS

/// set the parameters for the subtract electrostatic potential
int subt_elec_set_params(int bond_type);

/** Computes the negative of the electrostatic pair forces 
    and adds this force to the particle forces (see \ref tclcommand_inter). 
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param iaparams  Parameters of interaction
    @param dx        change in position
    @param force     force on particles
    @return true if bond is broken
*/

inline int calc_subt_elec_pair_force(Particle *p1, Particle *p2, double d[3], double force[3])
{
  double dist2 = sqrlen(d);
  double dist = sqrt(dist2);
  int j;
//  double fac1,fac2, adist, erfc_part_ri;
  double chgfac = p1->p.q*p2->p.q;
/*
  dist2 = pow(p1->r.p[0]-p2->r.p[0],2)+pow(p1->r.p[1]-p2->r.p[1],2)+pow(p1->r.p[2]-p2->r.p[2],2);
  dist = sqrt(dist2);
  if (dist>box_l[0]*0.5)
  	dist -= box_l[0]*0.5;
  dist2 = dist*dist;
*/
  for(j=0;j<3;j++)
	force[j] = -coulomb.prefactor * chgfac * d[j] / dist / dist2;

//  fprintf(stderr,"%e %e\n",d[2],force[2]) ;

/*  if(dist < p3m.params.r_cut) {
    if (dist > 0.0){		//Vincent
      adist = p3m.params.alpha * dist;
#if USE_ERFC_APPROXIMATION
      erfc_part_ri = AS_erfc_part(adist) / dist;
      fac1 = coulomb.prefactor * chgfac  * exp(-adist*adist);
      fac2 = fac1 * (erfc_part_ri + 2.0*p3m.params.alpha*wupii) / dist2;
#else
      erfc_part_ri = erfc(adist) / dist;
      fac1 = coulomb.prefactor * chgfac;
      fac2 = fac1 * (erfc_part_ri + 2.0*p3m.params.alpha*wupii*exp(-adist*adist)) / dist2;
#endif
      for(j=0;j<3;j++)
	force[j] = -fac2 * d[j];
      ESR_TRACE(fprintf(stderr,"%d: RSE: Pair dist=%.3f: force (%.3e,%.3e,%.3e)\n",this_node,
			dist,fac2*d[0],fac2*d[1],fac2*d[2]));
    }
  }
*/
  return 0;
}

inline int subt_elec_pair_energy(Particle *p1, Particle *p2, double d[3], double *_energy)
{
  double dist = normr(d);
  double chgfac = p1->p.q*p2->p.q;

  *_energy=0;

  if(dist < p3m.params.r_cut) {
    *_energy = -coulomb.prefactor * chgfac / dist;
  }

  return 0;
}

#endif

#endif
