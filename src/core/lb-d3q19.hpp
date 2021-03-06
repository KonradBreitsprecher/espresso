/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
/** \file lb-d3q19.hpp
 * Header file for the lattice Boltzmann D3Q19 model.
 *
 * This header file contains the definition of the D3Q19 model.
 */

#ifndef D3Q19_H
#define D3Q19_H

#ifdef LB
#include "lb.hpp"

/** Velocity sub-lattice of the D3Q19 model */
static double d3q19_lattice[19][3] = { {  0.,  0.,  0. },
                                       {  1.,  0.,  0. },
          			       { -1.,  0.,  0. },
          		               {  0.,  1.,  0. }, 
          			       {  0., -1.,  0. },
          			       {  0.,  0.,  1. }, 
          			       {  0.,  0., -1. },
          			       {  1.,  1.,  0. }, 
          			       { -1., -1.,  0. },
          			       {  1., -1.,  0. },
          			       { -1.,  1.,  0. },
          			       {  1.,  0.,  1. },
          			       { -1.,  0., -1. },
          			       {  1.,  0., -1. },
          			       { -1.,  0.,  1. },
          			       {  0.,  1.,  1. },
          			       {  0., -1., -1. },
          			       {  0.,  1., -1. },
          			       {  0., -1.,  1. } } ;

/** Coefficients for pseudo-equilibrium distribution of the D3Q19 model */
static double d3q19_coefficients[19][4] = { { 1./3.,  1.,     3./2., -1./2.  },
					    { 1./18., 1./6.,  1./4., -1./12. },
					    { 1./18., 1./6.,  1./4., -1./12. },
					    { 1./18., 1./6.,  1./4., -1./12. },
					    { 1./18., 1./6.,  1./4., -1./12. },
					    { 1./18., 1./6.,  1./4., -1./12. },
					    { 1./18., 1./6.,  1./4., -1./12. },
					    { 1./36., 1./12., 1./8., -1./24. },
					    { 1./36., 1./12., 1./8., -1./24. },
					    { 1./36., 1./12., 1./8., -1./24. },
					    { 1./36., 1./12., 1./8., -1./24. },
					    { 1./36., 1./12., 1./8., -1./24. },
					    { 1./36., 1./12., 1./8., -1./24. },
					    { 1./36., 1./12., 1./8., -1./24. },
					    { 1./36., 1./12., 1./8., -1./24. },
					    { 1./36., 1./12., 1./8., -1./24. },
					    { 1./36., 1./12., 1./8., -1./24. },
					    { 1./36., 1./12., 1./8., -1./24. },
					    { 1./36., 1./12., 1./8., -1./24. } };

/** Coefficients in the functional for the equilibrium distribution */
static double d3q19_w[19] = {1. / 3.,  1. / 18., 1. / 18., 1. / 18., 1. / 18.,
                             1. / 18., 1. / 18., 1. / 36., 1. / 36., 1. / 36.,
                             1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36.,
                             1. / 36., 1. / 36., 1. / 36., 1. / 36.};

/** Basis of the mode space as described in [Duenweg, Schiller, Ladd] */
static double d3q19_modebase[20][19] = {
  {  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0 },
  {  0.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0 },
  {  0.0,  0.0,  0.0,  1.0, -1.0,  0.0,  0.0,  1.0, -1.0, -1.0,  1.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0 },
  {  0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0, -1.0,  1.0,  1.0, -1.0, -1.0,  1.0 },
  { -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0 },
  {  0.0,  1.0,  1.0, -1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0,  1.0,  1.0, -1.0, -1.0, -1.0, -1.0 },
  { -0.0,  1.0,  1.0,  1.0,  1.0, -2.0, -2.0,  2.0,  2.0,  2.0,  2.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0 },
  {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0, -1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 },
  {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0, -1.0, -1.0,  0.0,  0.0,  0.0,  0.0 },
  {  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0, -1.0, -1.0 },
  {  0.0, -2.0,  2.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0 },
  {  0.0,  0.0,  0.0, -2.0,  2.0,  0.0,  0.0,  1.0, -1.0, -1.0,  1.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0 },
  {  0.0,  0.0,  0.0,  0.0,  0.0, -2.0,  2.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0, -1.0,  1.0,  1.0, -1.0, -1.0,  1.0 },
  {  0.0, -0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0, -1.0,  1.0, -1.0,  1.0,  0.0,  0.0,  0.0,  0.0 },
  {  0.0,  0.0,  0.0, -0.0,  0.0,  0.0,  0.0,  1.0, -1.0, -1.0,  1.0,  0.0,  0.0,  0.0,  0.0, -1.0,  1.0, -1.0,  1.0 },
  {  0.0,  0.0,  0.0,  0.0,  0.0, -0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0, -1.0,  1.0, -1.0,  1.0,  1.0, -1.0 },
  {  1.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0 },
  {  0.0, -1.0, -1.0,  1.0,  1.0, -0.0, -0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0,  1.0,  1.0, -1.0, -1.0, -1.0, -1.0 },
  {  0.0, -1.0, -1.0, -1.0, -1.0,  2.0,  2.0,  2.0,  2.0,  2.0,  2.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0 },
  /* the following values are the (weighted) lengths of the vectors */
  { 1.0, 1./3., 1./3., 1./3., 2./3., 4./9., 4./3., 1./9., 1./9., 1./9., 2./3., 2./3., 2./3., 2./9., 2./9., 2./9., 2.0, 4./9., 4./3. }
};

//LB_Model d3q19_model = { 19, d3q19_lattice, d3q19_coefficients, d3q19_w, nullptr, 1./3. };

#endif /* LB */

#endif /* D3Q19_H */

/*@}*/
