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

/** \file iccp3m.cpp
    Detailed Information about the method is included in the corresponding header file \ref iccp3m.hpp.

 */
#include "capacitor.hpp"
#include "utils.hpp"
#include "parser.hpp"
#include <string>
#include <sstream>
#include <iostream>
#include <iterator> 

#ifdef ELECTROSTATICS
/** Parses the generate_potential_from_mesh command.
 */
int tclcommand_gen_pot_from_mesh(ClientData data, Tcl_Interp *interp, int argc, char **argv) 
{
  //char buffer[TCL_DOUBLE_SPACE];
  std::string ext_pot_path; 
  int num_electrodes=0;
  std::vector<std::string> geofiles;
  std::vector<double> potentials;
  std::vector<int> bins;
  
  if(argc < 3) { 
         Tcl_AppendResult(interp, "Usage of generate_potential_from_mesh: NUM_ELECTRODES [(1-N)GEOMETRY-FILES] [(1-N)POTENTIALS] [BINS_X BINS_Y BINS_Z] PATH_TO_EXT_POT", (char *)NULL); 
         return (TCL_ERROR); 
   }
  else {
	  if (ARG_IS_I(1, num_electrodes)) {
		  argc-=2;
		  argv+=2;
	  } else {
		  Tcl_AppendResult(interp, "generate_potential_from_mesh: First argument has to be the number of electrodes", (char *)NULL); 
		  return (TCL_ERROR);
	  }
	  if (argc != 4) {
		  Tcl_AppendResult(interp, "Expecting [(1-N)GEOMETRY-FILES] [(1-N)POTENTIALS] [BIN_X BIN_Y BIN_Z] PATH_TO_EXT_POT",  (char *)NULL); 
		  return (TCL_ERROR);
	  } else {
		  std::istringstream ss(argv[0]);
		  std::istream_iterator<std::string> begin(ss), end;
		  std::vector<std::string> token(begin, end); 
		  for (int i = 0; i < num_electrodes; i++) {
			  if (file_exists(token[i].c_str())) {
				  geofiles.push_back(token[i]);
				  Tcl_AppendResult(interp, "generate_potential_from_mesh: added file ", token[i].c_str(), (char *)NULL); 
			  } else {
				  Tcl_AppendResult(interp, "generate_potential_from_mesh: Could not find file ", token[i].c_str(), (char *)NULL); 
				  return (TCL_ERROR);
			  }
		  }
		  DoubleList potList;
		  init_doublelist(&potList);
		  if (ARG_IS_DOUBLELIST(1, potList)) {
			  for (int i = 0; i < num_electrodes; i++) {
				  potentials.push_back(potList.e[i]);
				  fprintf(stderr,"pot %f\n", potentials[i]);
			  }
		  } else {
			  Tcl_AppendResult(interp, "Expecting potential list\n" , (char *)NULL);
			  return TCL_ERROR;
		  }
		  
		  IntList binList;
		  init_intlist(&binList);
		  if (ARG_IS_INTLIST(2, binList)) {
			  for (int i = 0; i < 3; i++) {
				  bins.push_back(binList.e[i]);
				  fprintf(stderr,"bin %d\n", binList.e[i]);
			  }
		  } else {
			  Tcl_AppendResult(interp, "Expecting bin list\n" , (char *)NULL);
			  return TCL_ERROR;
		  }
		  
		  ext_pot_path.assign(argv[3]);

	  }
  }   
  setup_capacitor(geofiles, potentials, bins, ext_pot_path);

  return TCL_OK;
}

#endif
