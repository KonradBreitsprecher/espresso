/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
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
/** \file global_lattice.hpp 
 *
 * Lattice class definition
 * Contains the lattice layout and pointers to the data fields.
 * For parallelization purposes, it is assumed that a halo region
 * surrounds the local lattice sites.
 */

#ifndef _GLOBAL_LATTICE_HPP
#define _GLOBAL_LATTICE_HPP

#include "grid.hpp"
#include "particle_data.hpp"

#define INTERPOLATION_LINEAR 1
#define index_t long

class GlobalLattice {
public:
    int global_grid[3];
    unsigned int dim;
    double agrid[3];/** lattice constant */

    int halo_grid[3] ;/** number of lattice sites in each direction including halo */
    int halo_size;/** halo size in all directions */

    double offset[3];/** global offset */

    unsigned int interpolation_type;
    char flags;
    size_t element_size;/** Size of each element in size units (=bytes) */
    size_t lattice_dim;/** Dimension of the field, assuming entries are arrays */

    index_t grid_volume;/** total number (volume) of local lattice sites (excluding halo) */
    index_t halo_grid_volume;/** total number (volume) of lattice sites (including halo) */
    index_t halo_grid_surface;/** number of lattice sites in the halo region */
    index_t halo_offset;/** offset for number of halo sites stored in front of the local lattice sites */

    void *_data;/** pointer to the actual lattice data. This can be a contiguous field of arbitrary data. */

    /** particle representation of this lattice. This is needed to
     *  specify interactions between particles and the lattice.
     *  Actually used are only the identity and the type. */
    Particle part_rep;

    /* Constructor */
    GlobalLattice() {}

    /** Initialize lattice.
     *
     * This function initializes the variables describing the lattice
     * layout. Important: The lattice data is <em>not</em> allocated here!
     *
     * \param lattice pointer to the lattice
     * \param agrid   lattice spacing
     */
    int init(double* agrid, double* offset, int halo_size, size_t dim);

    /** lattice memory allocation.
     * \param lattice pointer to the lattice
     */
    void allocate_memory();
    void allocate_memory(size_t element_size);

    void interpolate(double* pos, double* value);

    void interpolate_linear(double* pos, double* value);

    void interpolate_gradient(double* pos, double* value);

    void interpolate_linear_gradient(double* pos, double* value);

    int set_data_for_global_index_with_periodic_image(int* pos, void* data);

    void get_data_for_halo_index(index_t* ind, void** data);

    void get_data_for_linear_index(index_t ind, void** data);

    void set_data_for_halo_grid_index(index_t* ind, void* data);

};

#endif /* GLOBAL_LATTICE_HPP */
