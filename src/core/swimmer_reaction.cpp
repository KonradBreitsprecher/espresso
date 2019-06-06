/*
  Copyright (C) 2010-2018 The ESPResSo project
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
/** \file
 *
 */

#include "swimmer_reaction.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "event.hpp"
#include "forces.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "random.hpp"

#include <utils/math/sgn.hpp>

#include <algorithm>
#include <random>
#include <vector>

reaction_struct reaction;

