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

#ifdef ELECTROSTATICS
/** Parses the generate_potential_from_mesh command.
 */
int tclcommand_gen_pot_from_mesh(ClientData data, Tcl_Interp *interp, int argc, char **argv) 
{
  //char buffer[TCL_DOUBLE_SPACE];
  int num_electrodes=0;
  std::vector<std::string> geofiles;
  std::vector<double> potentials;
  
  if(argc < 3) { 
         Tcl_AppendResult(interp, "Usage of capacitor: NUM_ELECTRODES (1-N)GEOMETRY-FILES (1-N)POTENTIALS", (char *)NULL); 
         return (TCL_ERROR); 
   }
  else {
	  if (ARG_IS_I(1, num_electrodes)) {
		  argc-=2;
		  argv+=2;
	  } else {
		  Tcl_AppendResult(interp, "capacitor: First argument has to be the number of electrodes", (char *)NULL); 
		  return (TCL_ERROR);
	  }
	  fprintf(stderr,"num_electrodes %d argc %d\n", num_electrodes, argc);
	  if (argc != num_electrodes * 2) {
		  Tcl_AppendResult(interp, "Wrong number of arguments in (1-N)GEOMETRY-FILES (1-N)POTENTIALS", (char *)NULL); 
		  return (TCL_ERROR);
	  } else {
		  for (int i = 0; i < num_electrodes; i++) {
			  if (file_exists(argv[i])) {
				  geofiles.push_back(argv[i]);
			  } else {
				  Tcl_AppendResult(interp, "capacitor: Could not find file ", argv[i], (char *)NULL); 
				  return (TCL_ERROR);
			  }
		  }
		  for (int i = num_electrodes; i < 2*num_electrodes; i++) {
			  double pot;
			  if (ARG_IS_D(i,pot)) {
				  potentials.push_back(pot);
			  } else {
				  Tcl_AppendResult(interp, "capacitor: Expecting potential, got ", argv[i], (char *)NULL); 
				  return (TCL_ERROR);
			  }
		  }
	  }
  }   
  setup_capacitor(geofiles, potentials);

  return TCL_OK;
}

#endif
