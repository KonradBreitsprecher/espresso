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

#ifndef TRIANGLEMESH_H
#define TRIANGLEMESH_H

#include <fstream>
#include <string>
#include <iterator>
#include <vector>
#include <iostream>
#include "utils.hpp"

typedef struct
{
    public:
        double normal[3];
        double vertices[3][3];
        double edges[3][3];
        double transformationMatrix[4][4];
        double transformedVertices[3][2];
        double edgeStart2D[9][2];
        double helperGradients[9][2];
        double area;
		double center[3];
        double clockwise;
} triangle;

class triangleMesh
{
    public:
		triangleMesh();
        triangleMesh(std::string pathToMeshfile);
		~triangleMesh();
		double sqrDistToMesh(double P[3]);
		bool isInside(double P[3]);
        triangle* _triangles;
        long _numFaces;
    	//static int triangleMesh_to_distanceMapLatticeFile(Constraint_mesh *c, char* filename);
    protected:
    private:
        long getNumFaces(std::string pathToMeshfile);
        void transformPoint(int i, double P[3], double pt[3]);
        void transformPoint2D(int i, double P[3], double pt2D[2]);
        double edgeEquation(int i, int j, double p[2]);
		void vecToFeatureOfFaceIndex(double P[3], int faceIndex, int feature, double distSqr, double res[3]);
		double sqrDistToEdge(double A[3],double a[3], double P[3]);
		double sqrDistToFeatureOfFaceIndex(double P[3], int faceIndex, int feature);	
		double sqrDistToTriangle(int i, double P[3], int *minDistFeature);
		void minVectorToMesh(double P[3], double res[3], int *minDistFaceIndex);
};

#endif // TRIANGLEMESH_H
