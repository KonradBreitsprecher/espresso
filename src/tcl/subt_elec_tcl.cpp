/*
  Copyright (C) 2010,2011,2012 The ESPResSo project
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
/** \file subt_elec_tcl.cpp
 *
 *  Implementation of \ref subt_elec_tcl.hpp
 */
#include "subt_elec_tcl.hpp"
#include "subt_elec.hpp"

#ifdef ELECTROSTATICS
    
int tclcommand_inter_parse_subt_elec(Tcl_Interp *interp, int bond_type,
				   int argc, char **argv)
{
  //fprintf(stderr,"%i\n",argc);
  if (argc != 1) {
    Tcl_AppendResult(interp, "subt_elec accepts no further parameters", (char *) NULL);
    return TCL_ERROR;
  }

  CHECK_VALUE(subt_elec_set_params(bond_type), "bond type must be nonnegative");
}

int tclprint_to_result_subt_elecIA(Tcl_Interp *interp)
{
  Tcl_AppendResult(interp, "subt_elec", (char *) NULL);
  return TCL_OK;
}

#endif

