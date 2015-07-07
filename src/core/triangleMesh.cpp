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

#include "triangleMesh.hpp"
#define COPYARRAY(dest, src) memcpy((dest), (src), sizeof((src)))

const std::string featureStr[7] = {"A","B","C","E1","E2","E3","F"};
double stringToDouble(std::string str)
{
	std::stringstream ss;
	ss << str;
	double result;
	ss >> result;
	return result;
}

triangleMesh::triangleMesh(std::string pathToMeshfile)
{
	_numFaces = getNumFaces(pathToMeshfile);
	fprintf(stderr, "FOUND %d FACES",_numFaces);
	_triangles = new triangle[_numFaces];
	std::ifstream inFile(pathToMeshfile.c_str());
	std::string line;
	double tmpVertices[3][3] = { {0}, {0}, {0} };
	int vertexCnt = 0;
	int faceCnt = 0;
	while (std::getline(inFile, line))
	{
		std::istringstream ss(line);
		std::istream_iterator < std::string > begin(ss), end;
		std::vector < std::string > arrayTokens(begin, end);
		if (arrayTokens[0] == "facet")
		{
			_triangles[faceCnt].normal[0] = stringToDouble(arrayTokens[2]);
			_triangles[faceCnt].normal[1] = stringToDouble(arrayTokens[3]);
			_triangles[faceCnt].normal[2] = stringToDouble(arrayTokens[4]);
		}
		else if (arrayTokens[0] == "vertex")
		{
			tmpVertices[vertexCnt][0] = stringToDouble(arrayTokens[1]);
			tmpVertices[vertexCnt][1] = stringToDouble(arrayTokens[2]);
			tmpVertices[vertexCnt][2] = stringToDouble(arrayTokens[3]);
			vertexCnt++;
			if (vertexCnt == 3)
			{
				// Vertices
				COPYARRAY(_triangles[faceCnt].vertices, tmpVertices);
				// Edges
				vecsub(tmpVertices[1], tmpVertices[0],_triangles[faceCnt].edges[0]);
				vecsub(tmpVertices[2], tmpVertices[1],_triangles[faceCnt].edges[1]);
				vecsub(tmpVertices[0], tmpVertices[2],_triangles[faceCnt].edges[2]);
				// TransformationMatrix
//				get_n_triangle(_triangles[faceCnt].vertices[2],_triangles[faceCnt].vertices[1],_triangles[faceCnt].vertices[0],_triangles[faceCnt].normal);
//				vecscale(_triangles[faceCnt].normal,1.0/normr(_triangles[faceCnt].normal));
				double nx = _triangles[faceCnt].normal[0];
				double ny = _triangles[faceCnt].normal[1];
				double nz = _triangles[faceCnt].normal[2];

				double s = sin(acos(-nx));
				double transformationMatrix[4][4]= 
					{{-nx,-ny*s,-nz*s,nx*_triangles[faceCnt].vertices[0][0]+ny*s*_triangles[faceCnt].vertices[0][1]+nz*s*_triangles[faceCnt].vertices[0][2]},
					{ny*s,-nx+nz*nz*(1+nx),-nz*ny*(1+nx),-ny*s*_triangles[faceCnt].vertices[0][0]-(-nx+nz*nz*(1+nx))*_triangles[faceCnt].vertices[0][1]-(-nz*ny*(1+nx))*_triangles[faceCnt].vertices[0][2]},
					{nz*s,-ny*nz*(1+nx),-nx+ny*ny*(1+nx),-nz*s*_triangles[faceCnt].vertices[0][0]-(-ny*nz*(1+nx))*_triangles[faceCnt].vertices[0][1]-(-nx+ny*ny*(1+nx))*_triangles[faceCnt].vertices[0][2]},
					{0,0,0,1}};
				COPYARRAY(_triangles[faceCnt].transformationMatrix, transformationMatrix);
				// Transform Vertices to 2D
				transformPoint2D(faceCnt, _triangles[faceCnt].vertices[0], _triangles[faceCnt].transformedVertices[0]);
				transformPoint2D(faceCnt, _triangles[faceCnt].vertices[1], _triangles[faceCnt].transformedVertices[1]);
				transformPoint2D(faceCnt, _triangles[faceCnt].vertices[2], _triangles[faceCnt].transformedVertices[2]);
				// Starting point of 2D edges
				COPYARRAY(_triangles[faceCnt].edgeStart2D[0],_triangles[faceCnt].transformedVertices[0]);
				COPYARRAY(_triangles[faceCnt].edgeStart2D[1],_triangles[faceCnt].transformedVertices[1]);
				COPYARRAY(_triangles[faceCnt].edgeStart2D[2],_triangles[faceCnt].transformedVertices[2]);
				COPYARRAY(_triangles[faceCnt].edgeStart2D[3],_triangles[faceCnt].transformedVertices[0]);
				COPYARRAY(_triangles[faceCnt].edgeStart2D[4],_triangles[faceCnt].transformedVertices[1]);
				COPYARRAY(_triangles[faceCnt].edgeStart2D[5],_triangles[faceCnt].transformedVertices[1]);
				COPYARRAY(_triangles[faceCnt].edgeStart2D[6],_triangles[faceCnt].transformedVertices[2]);
				COPYARRAY(_triangles[faceCnt].edgeStart2D[7],_triangles[faceCnt].transformedVertices[2]);
				COPYARRAY(_triangles[faceCnt].edgeStart2D[8],_triangles[faceCnt].transformedVertices[0]);
				//Direction vector (gradients) of 2D edges
				vecsub(_triangles[faceCnt].transformedVertices[1], _triangles[faceCnt].transformedVertices[0], _triangles[faceCnt].helperGradients[0]);
				vecsub(_triangles[faceCnt].transformedVertices[2], _triangles[faceCnt].transformedVertices[1], _triangles[faceCnt].helperGradients[1]);
				vecsub(_triangles[faceCnt].transformedVertices[0], _triangles[faceCnt].transformedVertices[2], _triangles[faceCnt].helperGradients[2]);
				_triangles[faceCnt].helperGradients[3][0] =  _triangles[faceCnt].helperGradients[0][1];
				_triangles[faceCnt].helperGradients[3][1] = -_triangles[faceCnt].helperGradients[0][0];
				_triangles[faceCnt].helperGradients[4][0] =  _triangles[faceCnt].helperGradients[0][1];
				_triangles[faceCnt].helperGradients[4][1] = -_triangles[faceCnt].helperGradients[0][0];
				_triangles[faceCnt].helperGradients[5][0] =  _triangles[faceCnt].helperGradients[1][1];
				_triangles[faceCnt].helperGradients[5][1] = -_triangles[faceCnt].helperGradients[1][0];
				_triangles[faceCnt].helperGradients[6][0] =  _triangles[faceCnt].helperGradients[1][1];
				_triangles[faceCnt].helperGradients[6][1] = -_triangles[faceCnt].helperGradients[1][0];
				_triangles[faceCnt].helperGradients[7][0] =  _triangles[faceCnt].helperGradients[2][1];
				_triangles[faceCnt].helperGradients[7][1] = -_triangles[faceCnt].helperGradients[2][0];
				_triangles[faceCnt].helperGradients[8][0] =  _triangles[faceCnt].helperGradients[2][1];
				_triangles[faceCnt].helperGradients[8][1] = -_triangles[faceCnt].helperGradients[2][0];
				vertexCnt = 0;
				
				double edgeSum = (_triangles[faceCnt].transformedVertices[1][0]-_triangles[faceCnt].transformedVertices[0][0]) *
						         (_triangles[faceCnt].transformedVertices[1][1]+_triangles[faceCnt].transformedVertices[0][1]) +
			                     (_triangles[faceCnt].transformedVertices[2][0]-_triangles[faceCnt].transformedVertices[1][0]) *  
                                 (_triangles[faceCnt].transformedVertices[2][1]+_triangles[faceCnt].transformedVertices[1][1]) +
								 (_triangles[faceCnt].transformedVertices[0][0]-_triangles[faceCnt].transformedVertices[2][0]) *  
                                 (_triangles[faceCnt].transformedVertices[0][1]+_triangles[faceCnt].transformedVertices[2][1]); 

				if (edgeSum > 0) {
					//fprintf(stderr,"CLOCKWISE\n");
					_triangles[faceCnt].clockwise = -1;
				}
				else {
					//fprintf(stderr,"COUNTER CLOCKWISE\n");
					_triangles[faceCnt].clockwise = 1;
				}

				_triangles[faceCnt].center[0] = (_triangles[faceCnt].vertices[0][0]+_triangles[faceCnt].vertices[1][0]+_triangles[faceCnt].vertices[2][0])/3.0;
				_triangles[faceCnt].center[1] = (_triangles[faceCnt].vertices[0][1]+_triangles[faceCnt].vertices[1][1]+_triangles[faceCnt].vertices[2][1])/3.0;
				_triangles[faceCnt].center[2] = (_triangles[faceCnt].vertices[0][2]+_triangles[faceCnt].vertices[1][2]+_triangles[faceCnt].vertices[2][2])/3.0;

				_triangles[faceCnt].area = area_triangle(_triangles[faceCnt].vertices[0],_triangles[faceCnt].vertices[1],_triangles[faceCnt].vertices[2]);

				faceCnt++;
			}
		}
	} inFile.close();
	fprintf(stderr, "Read in mesh file complete\n");
} 

triangleMesh::~triangleMesh()
{
	fprintf(stderr, "Deleting triangleMesh\n");
	delete[]_triangles;
} 

long triangleMesh::getNumFaces(std::string pathToMeshfile)
{
	std::ifstream inFile(pathToMeshfile.c_str());
	std::string line;
	long faceCnt = 0;
	while (getline(inFile, line))
		if (line.find("facet normal") == 0)
			faceCnt++;
	return faceCnt;
}

void triangleMesh::transformPoint(int i, double P[3], double pt[3])
{
	pt[0] = _triangles[i].transformationMatrix[0][0] * P[0] +
		_triangles[i].transformationMatrix[0][1] * P[1] +
		_triangles[i].transformationMatrix[0][2] * P[2] +
		_triangles[i].transformationMatrix[0][3] * 1;
	pt[1] =	_triangles[i].transformationMatrix[1][0] * P[0] +
		_triangles[i].transformationMatrix[1][1] * P[1] +
		_triangles[i].transformationMatrix[1][2] * P[2] +
		_triangles[i].transformationMatrix[1][3] * 1;
	pt[2] =	_triangles[i].transformationMatrix[2][0] * P[0] +
		_triangles[i].transformationMatrix[2][1] * P[1] +
		_triangles[i].transformationMatrix[2][2] * P[2] +
		_triangles[i].transformationMatrix[2][3] * 1;
} 

void triangleMesh::transformPoint2D(int i, double P[3], double pt2D[2])
{
	double pt[3];
	transformPoint(i, P, pt);
	pt2D[0] = pt[2];
	pt2D[1] = pt[1];
}

double triangleMesh::edgeEquation(int i, int j, double p[2])
{
	return _triangles[i].clockwise * ((p[0] - _triangles[i].edgeStart2D[j][0]) * _triangles[i].helperGradients[j][1] -
							          (p[1] - _triangles[i].edgeStart2D[j][1]) * _triangles[i].helperGradients[j][0]);
}

double triangleMesh::sqrDistToEdge(double A[3],double a[3], double P[3]) 
{ 
	double b[3]; 
	vecsub(P, A, b); 
	return (pow(a[1]*b[2]-a[2]*b[1],2)+pow(a[2]*b[0]-a[0]*b[2],2)+pow(a[0]*b[1]-a[1]*b[0],2))/(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
} 


double triangleMesh::sqrDistToTriangle(int i, double P[3], int *minDistFeature)
{
	//Use 2D representation of triangle and test point to determine feature (A,B,C,E1,E2,E3,F), use 3D coords to calculate distance vector
	double pt[3];
	transformPoint(i, P, pt);
	double p[2];
	p[0] = pt[2];
	p[1] = pt[1];
	/*
	fprintf(stderr,"plot '-' w p ps 5 pt 2,'-' w p ps 5 pt 2,'-' w p ps 5 pt 2,'-' w p ps 5 pt 3\n%.2f %.2f\ne\n%.2f %.2f\ne\n%.2f %.2f\ne\n%.2f %.2f\ne\n",_triangles[i].edgeStart2D[0][0] ,_triangles[i].edgeStart2D[0][1],_triangles[i].edgeStart2D[1][0],_triangles[i].edgeStart2D[1][1],_triangles[i].edgeStart2D[2][0],_triangles[i].edgeStart2D[2][1], p[0],p[1]);
    for (int j = 0; j<9; j++) {	
		fprintf(stderr,"EDGE EQ %d:  %2.2f\n",j,edgeEquation(i, j, p));
	}
    for (int j = 0; j<6; j++) {	
		fprintf(stderr,"F%s %2.2f\n",featureStr[j].c_str(), sqrt(sqrDistToFeatureOfFaceIndex(P,i,j)));
	}
	fprintf(stderr,"F6 %2.2f\n",abs(pt[0]));
	*/
	if (edgeEquation(i, 0, p) <= 0 && edgeEquation(i, 1, p) <= 0 && edgeEquation(i, 2, p) <= 0)
	{
		//fprintf(stderr,"FEATURE F\n");
		*minDistFeature = 6;
		return pt[0]*pt[0];
	}
	else if (edgeEquation(i, 3, p) >= 0 && edgeEquation(i, 8, p) <= 0)
	{
		*minDistFeature = 0;
	}
	else if (edgeEquation(i, 5, p) >= 0 && edgeEquation(i, 4, p) <= 0)
	{
		*minDistFeature = 1;
	}
	else if (edgeEquation(i, 7, p) >= 0 && edgeEquation(i, 6, p) <= 0)
	{
		*minDistFeature = 2;
	}
	else if (edgeEquation(i, 0, p) >= 0 && edgeEquation(i, 3, p) <= 0 && edgeEquation(i, 4, p) >= 0)
	{
		*minDistFeature = 3;
	}
	else if (edgeEquation(i, 1, p) >= 0 && edgeEquation(i, 5, p) <= 0 && edgeEquation(i, 6, p) >= 0)
	{
		*minDistFeature = 4;
	}
	else if (edgeEquation(i, 2, p) >= 0 && edgeEquation(i, 7, p) <= 0 && edgeEquation(i, 8, p) >= 0)
	{
		*minDistFeature = 5;
	}

	//fprintf(stderr,"FEATURE %s\n",featureStr[*minDistFeature].c_str());
	return sqrDistToFeatureOfFaceIndex(P,i,*minDistFeature);
}

double triangleMesh::sqrDistToFeatureOfFaceIndex(double P[3], int faceIndex, int feature)
{
	if (feature < 3) /* Feature is triangle point */
	{
		double res[3];
		vecsub(P, _triangles[faceIndex].vertices[feature], res);
		return sqrlen(res);
	}
	else //if (feature < 6) /* Feature is Edge */
	{
		return sqrDistToEdge(_triangles[faceIndex].vertices[feature-3], _triangles[faceIndex].edges[feature-3], P);
	}
	//else /* Feature is Face */
	//{
	//}
}


void triangleMesh::minVectorToMesh(double P[3], double res[3])
{
	double distSqr=0;
	double minDistSqr = 1e300;
	int minDistFeature=0,minDistFaceIndex=0, feature=0;
	res[0] = res[1] = res[2] = minDistSqr;

	for (int i = 0; i < _numFaces; i++)
	{
		distSqr = sqrDistToTriangle(i, P, &feature);
		if (minDistSqr > distSqr)
		{
			minDistSqr = distSqr;
			minDistFaceIndex = i;
			minDistFeature = feature;
		}
	}
	vecToFeatureOfFaceIndex(P, minDistFaceIndex, minDistFeature, minDistSqr, res);
}

void triangleMesh::vecToFeatureOfFaceIndex(double P[3], int faceIndex, int feature, double distSqr, double res[3])
{

	if (feature < 3) /* Feature is triangle point */
	{
		vecsub(P, _triangles[faceIndex].vertices[feature], res);
	}
	else if (feature < 6) /* Feature is Edge */
	{
		vecdist_line_point(_triangles[faceIndex].vertices[feature-3], _triangles[faceIndex].edges[feature-3], P, res);
	}
	else /* Feature is Face */
	{
		res[0] = _triangles[faceIndex].normal[0]; 
		res[1] = _triangles[faceIndex].normal[1]; 
		res[2] = _triangles[faceIndex].normal[2]; 
		vecscale(res, sqrt(distSqr));
	}
}


double triangleMesh::sqrDistToMesh(double P[3])
{
	int feature,minDistFeature;
	double minDist = 1e10;
	for (int i = 0; i < _numFaces; i++)
	{
		double distSqr = sqrDistToTriangle(i, P, &feature);
		if (minDist > distSqr) {
			minDistFeature = feature;
			minDist = distSqr;
		}
	}

	return minDist;
}

/*
int triangleMesh::triangleMeshs_to_potentialFile(vector<string> filenames)
{
	//TODO: TEST IF POTENTIAL OUTPUT FILE ALREADY EXISTS


	FILE *testMeshFile = fopen(c->filename, "r");
	if (!testMeshFile)
	{
		ostringstream msg;
		msg << "Could not open file " << c->filename << "\n";
		runtimeError(msg);
		return ES_ERROR;
	}
	double P[3], reso[3], minVec[3];
	int i, j, k;
	triangleMesh *t = new triangleMesh(c->filename);
	for (i = 0; i < 3; i++)
	{
		reso[i] = local_box_l[i] / c->bins[i];
	}
	double estimatedTime = t->_numFaces*c->bins[0]*c->bins[1]*c->bins[1]/467440000*0.52;
	fprintf(stderr, "Read in mesh file with %d faces complete. Estimated time for volume data generation: %f min\n",t->_numFaces,estimatedTime);
	fprintf(stderr, "Writing file: %s\n", filename);
	// Header 
	FILE *outfile = fopen(filename, "w");
	fprintf(outfile, "3 %f %f %f %f %f %f 0 0 0\n", local_box_l[0], local_box_l[1], local_box_l[2], reso[0], reso[1], reso[2]);
	for (i = 0; i < c->bins[0]; i++)
	{
		fprintf(stderr, "%d%% complete\n",(int)(100.0*i/c->bins[0]));
		for (j = 0; j < c->bins[1]; j++)
		{
			for (k = 0; k < c->bins[2]; k++)
			{
				P[0] = i * local_box_l[0] / c->bins[0] + c->offset[0];
				P[1] = j * local_box_l[1] / c->bins[1] + c->offset[1];
				P[2] = k * local_box_l[2] / c->bins[2] + c->offset[2];
				t->minVectorToMesh(P, minVec);
				P[0] -= c->offset[0];
				P[1] -= c->offset[1];
				P[2] -= c->offset[2];
//				fprintf(stderr, "\nWrite distance map %f %f %f %f %f %f\n", P[0], P[1], P[2], minVec[0], minVec[1], minVec[2]);
				// Data 
				fprintf(outfile, "%f %f %f %f %f %f %f\n", P[0], P[1], P[2], minVec[0], minVec[1], minVec[2],normr(minVec));
			}
		}
	}
	fclose(outfile);
	delete t;
	fprintf(stderr, "Write distance map complete\n");
	return ES_OK;
}
*/
