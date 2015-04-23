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
/** \file lattice.cpp 
 *
 * Lattice class definition
 *
 */

#include "global_lattice.hpp"

int GlobalLattice::init(double *agrid, double* offset, int halo_size, size_t dim) {
    this->dim=dim;

    /* determine the number of local lattice nodes */
    for (int d=0; d<3; d++) {
        this->agrid[d] = agrid[d];
        this->global_grid[d] = (int)dround(box_l[d]/agrid[d]);
        this->offset[d]=offset[d];
    }
    
	this->element_size = this->dim*sizeof(double);

    LATTICE_TRACE(fprintf(stderr,"%d: global_grid (%d,%d,%d)  res (%.2f,%.2f,%.2f)  offset(%.2f,%.2f,%.2f)  halosize %d  dim %d\n",this_node, this->global_grid[0],this->global_grid[1],this->global_grid[2],this->agrid[0],this->agrid[1],this->agrid[2], this->offset[0], this->offset[1], this->offset[2], halo_size, this->dim));

    this->halo_size = halo_size;
    /* determine the number of total nodes including halo */
    this->halo_grid[0] = this->global_grid[0] + 2*halo_size ;
    this->halo_grid[1] = this->global_grid[1] + 2*halo_size ;
    this->halo_grid[2] = this->global_grid[2] + 2*halo_size ;

    this->grid_volume = this->global_grid[0]*this->global_grid[1]*this->global_grid[2] ;
    this->halo_grid_volume = this->halo_grid[0]*this->halo_grid[1]*this->halo_grid[2] ;
    this->halo_grid_surface = this->halo_grid_volume - this->grid_volume ;

    this->interpolation_type = INTERPOLATION_LINEAR;
	
    allocate_memory();
    return ES_OK;

}

void GlobalLattice::allocate_memory() {

    this->_data = malloc(this->element_size*this->halo_grid_volume);
	fprintf(stderr, "Lattice Malloc %f Mb\n",this->element_size*this->halo_grid_volume*1e-6);
    memset(this->_data, (unsigned int)(-1), this->element_size*this->halo_grid_volume);
    //memset(this->_data, 0, this->element_size*this->halo_grid_volume);
}

void GlobalLattice::interpolate(double* pos, double* value) {
    if (this->interpolation_type == INTERPOLATION_LINEAR) {
        interpolate_linear(pos, value);
    } else {
        ostringstream msg;
        msg <<"Unknown interpolation type";
        runtimeError(msg);
    }
}

void GlobalLattice::interpolate_linear(double* pos, double* value) {
    int left_halo_index[3];
    double d[3];
    if (this->halo_size <= 0) {
        ostringstream msg;
        msg <<"Error in interpolate_linear: halo size is 0";
        runtimeError(msg);
        return;
    }
    for (int dim = 0; dim<3; dim++) {
        left_halo_index[dim]=(int) floor(pos[dim]/this->agrid[dim]) + this->halo_size;
        d[dim]=(pos[dim]/this->agrid[dim] - floor(pos[dim]/this->agrid[dim]));
        if (left_halo_index[dim] < 0 || left_halo_index[dim] >= this->halo_grid[dim]) {
			fprintf(stderr,"Error in interpolate_linear: Particle out of range left_halo_index[%d] %d  halo_grid %d\nParticle at %f %f %f",dim, left_halo_index[dim],this->halo_grid[dim],pos[0],pos[1],pos[2]);
			fflush(stderr);
            ostringstream msg;
            msg <<"Error in interpolate_linear: Particle out of range: left_halo_index[" << dim << "] = " << left_halo_index[dim] << "\nParticle at " << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
            runtimeError(msg);
            return;
        }
    }
    double w[8];
    index_t index[8];
    w[0] = (1-d[0])*(1-d[1])*(1-d[2]);
    index[0]=get_linear_index(   left_halo_index[0], left_halo_index[1], left_halo_index[2], this->halo_grid);
    w[1] = ( +d[0])*(1-d[1])*(1-d[2]);
    index[1]=get_linear_index(   left_halo_index[0]+1, left_halo_index[1], left_halo_index[2], this->halo_grid);
    w[2] = (1-d[0])*( +d[1])*(1-d[2]);
    index[2]=get_linear_index(   left_halo_index[0], left_halo_index[1]+1, left_halo_index[2], this->halo_grid);
    w[3] = ( +d[0])*( +d[1])*(1-d[2]);
    index[3]=get_linear_index(   left_halo_index[0]+1, left_halo_index[1]+1, left_halo_index[2], this->halo_grid);

    w[4] = (1-d[0])*(1-d[1])*( +d[2]);
    index[4]=get_linear_index(   left_halo_index[0], left_halo_index[1], left_halo_index[2]+1, this->halo_grid);
    w[5] = ( +d[0])*(1-d[1])*( +d[2]);
    index[5]=get_linear_index(   left_halo_index[0]+1, left_halo_index[1], left_halo_index[2]+1, this->halo_grid);
    w[6] = (1-d[0])*( +d[1])*( +d[2]);
    index[6]=get_linear_index(   left_halo_index[0], left_halo_index[1]+1, left_halo_index[2]+1, this->halo_grid);
    w[7] = ( +d[0])*( +d[1])*( +d[2]);
    index[7]=get_linear_index(   left_halo_index[0]+1, left_halo_index[1]+1, left_halo_index[2]+1, this->halo_grid);

    for (unsigned int i = 0; i<this->dim; i++) {
        value[i] = 0;
    }

    double* local_value;
    for (unsigned int i=0; i<8; i++) {
        get_data_for_linear_index(index[i], (void**) &local_value);
        for (unsigned int j = 0; j<this->dim; j++) {
            value[j]+=w[i]*local_value[j];
        }
    }
}

void GlobalLattice::interpolate_gradient(double* pos, double* value) {
    if (this->interpolation_type == INTERPOLATION_LINEAR) {
        interpolate_linear_gradient(pos, value);
    } else {
        ostringstream msg;
        msg <<"Unknown interpolation type";
        runtimeError(msg);
    }
}

void GlobalLattice::interpolate_linear_gradient(double* pos, double* value) {
    int left_halo_index[3];
    double d[3];
    if (this->halo_size <= 0) {
        ostringstream msg;
        msg << "Error in interpolate_linear: halo size is 0";
        runtimeError(msg);
        return;
    }
    for (int dim = 0; dim<3; dim++) {
        left_halo_index[dim]=(int) floor(pos[dim]/this->agrid[dim]) + this->halo_size;
        d[dim]=pos[dim]/this->agrid[dim] - floor(pos[dim]/this->agrid[dim]);
        if (left_halo_index[dim] < 0 || left_halo_index[dim] >= this->halo_grid[dim]) {
			fprintf(stderr,"Error in interpolate_linear_gradient: Particle out of range agrid %f left_halo_index[%d] %d  halo_grid %d\nParticle at %f %f %f",this->agrid[dim], dim, left_halo_index[dim], this->halo_grid[dim],pos[0],pos[1],pos[2]);
			fflush(stderr);
            ostringstream msg;
            msg <<"Error in interpolate_linear_gradient: Particle out of range: left_halo_index[" << dim << "] = " << left_halo_index[dim] << "\nParticle at " << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
            runtimeError(msg);
            return;
        }
    }

    index_t index;
    double* local_value;

    for (unsigned int i = 0; i<3*this->dim; i++) {
        value[i] = 0;
    }

    index=get_linear_index(   left_halo_index[0], left_halo_index[1], left_halo_index[2], this->halo_grid);
    for (unsigned int i = 0; i<this->dim; i++) {
        get_data_for_linear_index(index, (void**) &local_value);
        value[3*i  ]+= (  -1  )*(1-d[1])*(1-d[2]) * local_value[i] / this->agrid[0];
        value[3*i+1]+= (1-d[0])*( -1   )*(1-d[2]) * local_value[i] / this->agrid[1];
        value[3*i+2]+= (1-d[0])*(1-d[1])*(  -1  ) * local_value[i] / this->agrid[2];
    }
    index=get_linear_index(   left_halo_index[0]+1, left_halo_index[1], left_halo_index[2], this->halo_grid);
    for (unsigned int i = 0; i<this->dim; i++) {
        get_data_for_linear_index(index, (void**) &local_value);
        value[3*i  ]+= (  +1  )*(1-d[1])*(1-d[2]) * local_value[i] / this->agrid[0];
        value[3*i+1]+= ( +d[0])*( -1   )*(1-d[2]) * local_value[i] / this->agrid[1];
        value[3*i+2]+= ( +d[0])*(1-d[1])*(  -1  ) * local_value[i] / this->agrid[2];
    }
    index=get_linear_index(   left_halo_index[0], left_halo_index[1]+1, left_halo_index[2], this->halo_grid);
    for (unsigned int i = 0; i<this->dim; i++) {
        get_data_for_linear_index(index, (void**) &local_value);
        value[3*i  ]+= (  -1  )*( +d[1])*(1-d[2]) * local_value[i] / this->agrid[0];
        value[3*i+1]+= (1-d[0])*( +1   )*(1-d[2]) * local_value[i] / this->agrid[1];
        value[3*i+2]+= (1-d[0])*( +d[1])*(  -1  ) * local_value[i] / this->agrid[2];
    }
    index=get_linear_index(   left_halo_index[0]+1, left_halo_index[1]+1, left_halo_index[2], this->halo_grid);
    for (unsigned int i = 0; i<this->dim; i++) {
        get_data_for_linear_index(index, (void**) &local_value);
        value[3*i  ]+= (  +1  )*( +d[1])*(1-d[2]) * local_value[i] / this->agrid[0];
        value[3*i+1]+= ( +d[0])*( +1   )*(1-d[2]) * local_value[i] / this->agrid[1];
        value[3*i+2]+= ( +d[0])*( +d[1])*(  -1  ) * local_value[i] / this->agrid[2];
    }
    index=get_linear_index(   left_halo_index[0]  , left_halo_index[1]  , left_halo_index[2] + 1, this->halo_grid);
    for (unsigned int i = 0; i<this->dim; i++) {
        get_data_for_linear_index(index, (void**) &local_value);
        value[3*i  ]+= (  -1  )*(1-d[1])*( +d[2]) * local_value[i] / this->agrid[0];
        value[3*i+1]+= (1-d[0])*( -1   )*( +d[2]) * local_value[i] / this->agrid[1];
        value[3*i+2]+= (1-d[0])*(1-d[1])*(  +1  ) * local_value[i] / this->agrid[2];
    }
    index=get_linear_index(   left_halo_index[0]+1, left_halo_index[1], left_halo_index[2]+1, this->halo_grid);
    for (unsigned int i = 0; i<this->dim; i++) {
        get_data_for_linear_index(index, (void**) &local_value);
        value[3*i  ]+= (  +1  )*(1-d[1])*( +d[2]) * local_value[i] / this->agrid[0];
        value[3*i+1]+= ( +d[0])*( -1   )*( +d[2]) * local_value[i] / this->agrid[1];
        value[3*i+2]+= ( +d[0])*(1-d[1])*(  +1  ) * local_value[i] / this->agrid[2];
    }
    index=get_linear_index(   left_halo_index[0], left_halo_index[1]+1, left_halo_index[2]+1, this->halo_grid);
    for (unsigned int i = 0; i<this->dim; i++) {
        get_data_for_linear_index(index, (void**) &local_value);
        value[3*i  ]+= (  -1  )*( +d[1])*( +d[2]) * local_value[i] / this->agrid[0];
        value[3*i+1]+= (1-d[0])*( +1   )*( +d[2]) * local_value[i] / this->agrid[1];
        value[3*i+2]+= (1-d[0])*( +d[1])*(  +1  ) * local_value[i] / this->agrid[2];
    }
    index=get_linear_index(   left_halo_index[0]+1, left_halo_index[1]+1, left_halo_index[2]+1, this->halo_grid);
    for (unsigned int i = 0; i<this->dim; i++) {
        get_data_for_linear_index(index, (void**) &local_value);
        value[3*i  ]+= (  +1  )*( +d[1])*( +d[2]) * local_value[i] / this->agrid[0];
        value[3*i+1]+= ( +d[0])*( +1   )*( +d[2]) * local_value[i] / this->agrid[1];
        value[3*i+2]+= ( +d[0])*( +d[1])*(  +1  ) * local_value[i] / this->agrid[2];
    }

}

int GlobalLattice::set_data_for_global_index_with_periodic_image(int* ind, void* data) {

    index_t halo_index[3];

    for (int i = 0; i<3; i++) {
		if (ind[i]<0 || ind[i]>=this->global_grid[i]) {
			fprintf(stderr, "ind[0]=%d  global_grid[0]=%d  agrid[0]=%f  offset[0]=%f\n",ind[0], this->global_grid[0], this->agrid[0], this->offset[0]);
			fprintf(stderr, "ind[1]=%d  global_grid[1]=%d  agrid[1]=%f  offset[1]=%f\n",ind[1], this->global_grid[1], this->agrid[1], this->offset[1]);
			fprintf(stderr, "ind[2]=%d  global_grid[2]=%d  agrid[2]=%f  offset[2]=%f\n",ind[2], this->global_grid[2], this->agrid[2], this->offset[2]);
		}
    }
	
	int out = 0;
	int cnt = 0;
    for (int i=-this->halo_size; i<=this->halo_size; i++) {
        for (int j=-this->halo_size; j<=this->halo_size; j++) {
            for (int k=-this->halo_size; k<=this->halo_size; k++) {
                halo_index[0]=ind[0]+i*this->global_grid[0]+this->halo_size;
                halo_index[1]=ind[1]+j*this->global_grid[1]+this->halo_size;
                halo_index[2]=ind[2]+k*this->global_grid[2]+this->halo_size;
				for (int d=0; d<3; d++) {
					if (halo_index[d] < 0 || halo_index[d] >= this->halo_grid[d])
					{
						out=1;
						break;
					}
				}
				if (out==0) {
					set_data_for_halo_grid_index(halo_index, data);
					cnt++;
				}
				out = 0;

            }
        }
    }
	return cnt;
}


void GlobalLattice::get_data_for_halo_index(index_t* ind, void** data) {
    (*data) = ((char*)this->_data) + get_linear_index(ind[0], ind[1], ind[2], this->halo_grid)*this->element_size;
}

void GlobalLattice::get_data_for_linear_index(index_t ind, void** data) {
	if (ind < 0 || ind >= this->halo_grid_volume) {
		fprintf(stderr,"get_data linear_index out of range: %d\n",ind);
	}
    (*data) = ((char*)this->_data) + ind*this->element_size;
}

void GlobalLattice::set_data_for_halo_grid_index(index_t* ind, void* data) {
    memcpy(((char*)this->_data) + get_linear_index(ind[0], ind[1], ind[2],  this->halo_grid)*this->element_size, data, this->element_size);

}

