/*
  Copyright (C) 2014 The ESPResSo project

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
#include "external_potential.hpp"
#include "global_lattice.hpp"
#include "communication.hpp"
#include "integrate.hpp"

ExternalPotential* external_potentials;
int n_external_potentials;

void external_potential_pre_init() {
  external_potentials = NULL;
  n_external_potentials = 0;
}


int generate_external_potential(ExternalPotential** e) {
  external_potentials = (ExternalPotential*) Utils::realloc(external_potentials,
		  (n_external_potentials+1) * sizeof(ExternalPotential));
  *e = &external_potentials[n_external_potentials];
  n_external_potentials++;
  (*e)->energy = 0;


//  e = &external_potentials[n_external_potentials-1];
  return ES_OK;
}     

int external_potential_tabulated_init(int number, char* filename, int n_particle_types, double* scale) {
  ExternalPotentialTabulated* e = &external_potentials[number].tabulated;

  if (strlen(filename)>MAX_FILENAME_SIZE)
    return ES_ERROR;
  strcpy((char*)&(e->filename), filename);
  external_potentials[number].type=EXTERNAL_POTENTIAL_TYPE_TABULATED;
  external_potentials[number].scale = (double*) Utils::malloc(n_particle_types*sizeof(double));
  external_potentials[number].n_particle_types = n_particle_types;
  for (int i = 0; i < n_particle_types; i++) {
    external_potentials[number].scale[i]=scale[i];
  }
  mpi_external_potential_broadcast(number);
  mpi_external_potential_tabulated_read_potential_file(number);
  return ES_OK;
}

int lattice_read_file(GlobalLattice* lattice, char* filename);

int external_potential_tabulated_read_potential_file(int number) {
  return lattice_read_file(&(external_potentials[number].tabulated.potential),
		  external_potentials[number].tabulated.filename);
}

int lattice_read_file(GlobalLattice* lattice, char* filename) {
 // ExternalPotentialTabulated *e = &(external_potentials[number].e.tabulated);
  FILE* infile = fopen(filename, "r");
  
  if (!infile)  {
      ostringstream msg;
      msg <<"Could not open file "<< filename << "\n";
      runtimeError(msg);
    return ES_ERROR;
  }
  char first_line[100];
  char* token;
  int bins[3];
  double res[3];
  double size[3];
  double offset[3]={0,0,0};
  int dim=0;
  if (fgets(first_line, 100, infile) == NULL) {
      fprintf(stderr, "Nothing read from file\n");
      return ES_ERROR;
  }

  token = strtok(first_line, " \t");
  if (!token) { fprintf(stderr, "Error reading dimensionality\n"); return ES_ERROR; }
  dim = atoi(token);
  if (dim<=0)  { fprintf(stderr, "Error reading dimensionality\n"); return ES_ERROR; }
  
  token = strtok(NULL, " \t");
  if (!token) { fprintf(stderr, "Could not read box_l[0]\n"); return ES_ERROR; }
  size[0] = atof(token);

  token = strtok(NULL, " \t");
  if (!token) { fprintf(stderr, "Could not read box_l[1]\n"); return ES_ERROR; }
  size[1] = atof(token);
  
  token = strtok(NULL, " \t");
  if (!token) { fprintf(stderr, "Could not read box_l[2]\n"); return ES_ERROR;}
  size[2] = atof(token);

  token = strtok(NULL, " \t");
  if (!token) { fprintf(stderr, "Could not read bin[0]\n"); return ES_ERROR;}
  bins[0] = atoi(token);
  
  token = strtok(NULL, " \t");
  if (!token) { fprintf(stderr, "Could not read bin[1]\n"); return ES_ERROR;}
  bins[1] = atoi(token);
  
  token = strtok(NULL, " \t");
  if (!token) { fprintf(stderr, "Could not read bin[2]\n"); return ES_ERROR;}
  bins[2] = atoi(token);

  token = strtok(NULL, " \t");
  if (token) {
    offset[0]=atof(token);
    token = strtok(NULL, " \t");
    if (!token) { fprintf(stderr, "Could not read offset[1]\n"); return ES_ERROR;}
    offset[1] = atof(token);
    token = strtok(NULL, " \t");
    if (!token) { fprintf(stderr, "Could not read offset[2]\n"); return ES_ERROR;}
    offset[2] = atof(token);
  }
  lattice->offset[0]=offset[0];
  lattice->offset[1]=offset[1];
  lattice->offset[2]=offset[2];


  if (size[0] > 0 && abs(size[0] - box_l[0]) > ROUND_ERROR_PREC) {
      ostringstream msg;
      msg <<"Box size in x is wrong "<< size[0] << " vs " << box_l[0] <<"\n";
      runtimeError(msg);
    return ES_ERROR;
  }
  if (size[1] > 0 && abs(size[1] - box_l[1]) > ROUND_ERROR_PREC) {
    ostringstream msg;
    msg <<"Box size in y is wrong "<< size[1] << " vs " << box_l[1] <<"\n";
    runtimeError(msg);
    return ES_ERROR;
  }
  if (size[2] > 0 && abs(size[2] - box_l[2]) > ROUND_ERROR_PREC) {
    ostringstream msg;
    msg <<"Box size in z is wrong "<< size[2] << " vs " << box_l[2] <<"\n";
    runtimeError(msg);
    return ES_ERROR;
  }

  res[0] = box_l[0]/bins[0];
  res[1] = box_l[1]/bins[1];
  res[2] = box_l[2]/bins[2];

  // Now we count how many entries we have:
  int halosize=1;

  lattice->init(res, offset, halosize, dim);
  lattice->interpolation_type = INTERPOLATION_LINEAR;

  char* line = (char*) malloc((3+dim)*ES_DOUBLE_SPACE);
  int ind[3] = {0,0,0};
  double f[3] = {0,0,0};
  int i;
  int cnt = 0;
  while (fgets(line, 200, infile)) {
    token = strtok(line, "\t");
    for (i=0; i<dim;i++) {
      if (!token) { fprintf(stderr, "Could not read f[%d]  index %d %d %d  token %s  got %s\n", i, ind[0], ind[1], ind[2], token, line); return ES_ERROR; }
      f[i] = atof(token);
      token = strtok(NULL, "\t");
    }
	cnt += lattice->set_data_for_global_index_with_periodic_image(ind, f);
	ind[2]++;
	if (ind[2]==bins[2])
	{
		ind[1]++;
		ind[2]=0;
		if (ind[1]==bins[1])
		{
			ind[0]++;
			ind[1]=0;
		}
	}
  }
  free(line);
  fclose(infile);
  fprintf(stderr,"Set %d halo cells, volume %d\n",cnt,(bins[0]+2*halosize)*(bins[1]+2*halosize)*(bins[2]+2*halosize));
  
  if (check_runtime_errors()!=0)
    return ES_ERROR;
  return ES_OK;
}



void add_external_potential_tabulated_forces(ExternalPotential* e, Particle* p) {
  //fprintf(stderr,"pos: %f %f %f  type %d  scale %f n_particle_types %d\n", p->r.p[0],p->r.p[1],p->r.p[2],p->p.type,e->scale[p->p.type], e->n_particle_types);
  if (p->p.type >= e->n_particle_types || e->scale[p->p.type] == 0 ) {
    return;
  }
  double field[3];
  double ppos[3];
  int    img[3];
  memmove(ppos, p->r.p, 3*sizeof(double));
  memmove(img, p->l.i, 3*sizeof(int));
  fold_position(ppos, img);
  e->tabulated.potential.interpolate_gradient(ppos, field);
  p->f.f[0]-=e->scale[p->p.type]*field[0];
  p->f.f[1]-=e->scale[p->p.type]*field[1];
  p->f.f[2]-=e->scale[p->p.type]*field[2];
  //fprintf(stderr,"%d %f force: %f %f %f\n", p->p.type, e->scale[p->p.type], e->scale[p->p.type]*field[0], e->scale[p->p.type]*field[1], e->scale[p->p.type]*field[2]);
  //fprintf(stderr,"pos: %f %f %f\n", p->r.p[0],p->r.p[1],p->r.p[2]);
  //fprintf(stderr,"fpos: %f %f %f\n", ppos[0],ppos[1],ppos[2]);
}

void add_external_potential_forces(Particle* p) {
  for (int i = 0; i < n_external_potentials; i++) {
    if (external_potentials[i].type==EXTERNAL_POTENTIAL_TYPE_TABULATED) {
      add_external_potential_tabulated_forces(&external_potentials[i], p);
    } else {
        ostringstream msg;
        msg <<"unknown external potential type";
        runtimeError(msg);
      return;
    }
  }
}


void add_external_potential_tabulated_energy(ExternalPotential* e, Particle* p) {
  if (p->p.type >= e->n_particle_types || e->scale[p->p.type] == 0) {
    return;
  }
  double potential;
  double ppos[3];
  int img[3];
  memmove(ppos, p->r.p, 3*sizeof(double));
  memmove(img, p->l.i, 3*sizeof(int));
  fold_position(ppos, img);
 
  e->tabulated.potential.interpolate(ppos, &potential);
  e->energy += e->scale[p->p.type] * potential;
}

void add_external_potential_energy(Particle* p) {
  for (int i=0; i<n_external_potentials; i++) {
    if (external_potentials[i].type==EXTERNAL_POTENTIAL_TYPE_TABULATED) {
      add_external_potential_tabulated_energy(&external_potentials[i], p);
    } else {
        ostringstream msg;
        msg <<"unknown external potential type";
        runtimeError(msg);
      return;
    }
  }
}

void external_potential_init_energies() {
  for (int i = 0; i<n_external_potentials; i++) {
    external_potentials[i].energy=0;
  }
}

