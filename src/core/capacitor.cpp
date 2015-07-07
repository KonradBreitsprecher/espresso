#include "capacitor.hpp"

//System size
/*
double* _boxMin;
double* _boxMax;
*/
/*
*/

//double ***_T, ***_Tnew;
//std::vector< std::vector< std::vector<double> > > _T;
//std::vector< std::vector< std::vector<double> > > _Tnew;

inline void print_dVec(double* v)
{
    std::cout << v[0] << " " << v[1] << " " << v[2] << std::endl;
}
inline void print_iVec(int* v)
{
    std::cout << v[0] << " " << v[1] << " " << v[2] << std::endl;
}

inline std::string intToString(int i)
{
    std::string s;
    std::stringstream out;
    out << i;
    return out.str();
}


capacitor::capacitor(std::vector<std::string> geofiles, std::vector<double> potentials, std::vector<int> bins, double surface_prec,int num_iter, double convergence, double eps_0)
{
	_box[0] = box_l[0];
	_box[1] = box_l[1];
	_box[2] = box_l[2];
	
    _bins[0] = bins[0];
    _bins[1] = bins[1];
    _bins[2] = bins[2];

	_surface_prec = surface_prec;
	_num_iter = num_iter;
	_convergence = convergence;
	_eps_0 = eps_0;

	_pref[0] = 1.0/6.0;
	_pref[1] = 1.0/6.0;
	_pref[2] = 1.0/6.0;
    
	/*
    double _reso[] = {_box[0] / _bins[0],_box[1] / _bins[1],_box[2] / _bins[2]};
    double _resoSqr[] = {_reso[0]*_reso[0],
                         _reso[1]*_reso[1],
                         _reso[2]*_reso[2]};

    double twoResoDiag = 1.5 / (_resoSqr[0]+_resoSqr[1]+_resoSqr[2]);

    _pref = new double[3] { _resoSqr[1]* _resoSqr[2] * twoResoDiag,
                            _resoSqr[0]* _resoSqr[2] * twoResoDiag,
                            _resoSqr[0]* _resoSqr[1] * twoResoDiag};
    */
	
	int num_electrodes = geofiles.size();
	_electrodes.reserve(num_electrodes);
	for (int i = 0; i < num_electrodes; i++) {
		electrode *e = new electrode(geofiles[i],potentials[i]);
		_electrodes.push_back(*e);
		fprintf(stderr, "Added electrode %d\n",i);
	}
}

void capacitor::worldToGrid(double worldPoint[3],int gridPoint[3])
{
     gridPoint[0] = (int)(dround(((worldPoint[0]) / _box[0]) * _bins[0])) % _bins[0];
     gridPoint[1] = (int)(dround(((worldPoint[1]) / _box[1]) * _bins[1])) % _bins[1];
     gridPoint[2] = (int)(dround(((worldPoint[2]) / _box[2]) * _bins[2])) % _bins[2];
}

void capacitor::gridToWorld(double worldPoint[3], int gridPoint[3])
{
    worldPoint[0] = gridPoint[0] * _box[0] / _bins[0];
    worldPoint[1] = gridPoint[1] * _box[1] / _bins[1];
    worldPoint[2] = gridPoint[2] * _box[2] / _bins[2];
}

int capacitor::gridToFlatArrayIndex(int* gridPoint)
{
    return gridPoint[2] + gridPoint[1] * _bins[2] + gridPoint[0] * _bins[2] * _bins[1];
}

int capacitor::translatedGrid(int* G, int dim, int ds)
{
    return (G[dim] + _bins[dim] + ds) % _bins[dim];
}

double capacitor::getNeighbourSum(double* data, int* G)
{
    //std::cout << data[gridToFlatArrayIndex(translateGrid(G,2,-1))] << std::endl;
    int GxP[3] = {translatedGrid(G,0,1),G[1],G[2]};
    int GxN[3] = {translatedGrid(G,0,-1),G[1],G[2]};
    int GyP[3] = {G[0],translatedGrid(G,1,1),G[2]};
    int GyN[3] = {G[0],translatedGrid(G,1,-1),G[2]};
    int GzP[3] = {G[0],G[1],translatedGrid(G,2,1)};
    int GzN[3] = {G[0],G[1],translatedGrid(G,2,-1)};

    return _pref[0] * (data[gridToFlatArrayIndex(GxP)] + data[gridToFlatArrayIndex(GxN)] +
					   data[gridToFlatArrayIndex(GyP)] + data[gridToFlatArrayIndex(GyN)] +
					   data[gridToFlatArrayIndex(GzP)] + data[gridToFlatArrayIndex(GzN)]);

                    //double* a = new double[3] {1,1,1};
                    //double b[3] = {1,1,1};
}

double capacitor::interpolatePotOnWorldPoint(double* data, double wP[3])
{
	//The lower grid point without periodic BCs
	int lGg[3] = { (int)(wP[0] / _box[0] * _bins[0]),
		     	   (int)(wP[1] / _box[1] * _bins[1]),
				   (int)(wP[2] / _box[2] * _bins[2])};

	//The lower grid point with periodic BCS
	int lG[3] = { lGg[0] % _bins[0],
			      lGg[1] % _bins[1],
				  lGg[2] % _bins[2]};
    
	//The higher grid points with periodic BCs
	int Gx[3] = {translatedGrid(lG,0,1),lG[1],lG[2]};
    int Gy[3] = {lG[0],translatedGrid(lG,1,1),lG[2]};
    int Gz[3] = {lG[0],lG[1],translatedGrid(lG,2,1)};
    int Gxy[3] = {translatedGrid(lG,0,1),translatedGrid(lG,1,1),lG[2]};
    int Gxz[3] = {translatedGrid(lG,0,1),lG[1],translatedGrid(lG,2,1)};
    int Gyz[3] = {lG[0],translatedGrid(lG,1,1),translatedGrid(lG,2,1)};
    int Gxyz[3] = {translatedGrid(lG,0,1),translatedGrid(lG,1,1),translatedGrid(lG,2,1)};
	
	//Potentials at the grid points
	double pots[8] = {data[gridToFlatArrayIndex(lG)],
			   		  data[gridToFlatArrayIndex(Gx)],
			   		  data[gridToFlatArrayIndex(Gy)],
			   		  data[gridToFlatArrayIndex(Gz)],
			   		  data[gridToFlatArrayIndex(Gxy)],
			   		  data[gridToFlatArrayIndex(Gxz)],
			   		  data[gridToFlatArrayIndex(Gyz)],
			   		  data[gridToFlatArrayIndex(Gxyz)]};
	
	//The grid points without periodic BCs
	int G[8][3] = {{lGg[0],  lGg[1],  lGg[2]  },
				   {lGg[0]+1,lGg[1],  lGg[2]  },
				   {lGg[0],  lGg[1]+1,lGg[2]  },
				   {lGg[0],  lGg[1],  lGg[2]+1},
				   {lGg[0]+1,lGg[1]+1,lGg[2]  },
				   {lGg[0]+1,lGg[1],  lGg[2]+1},
				   {lGg[0],  lGg[1]+1,lGg[2]+1},
				   {lGg[0]+1,lGg[1]+1,lGg[2]+1}};

	//The world points without periodic BCs
	double W[8][3];
	//The weights according to the volume of the associated rectangle
	double we[8];
	//The resulting potential
	double res = 0;
	//std::cout << "Weights: " << std::endl;
	for (int i = 0; i < 8; i++)
	{
		gridToWorld(W[i],G[i]);
		we[i] = abs((W[i][0] - wP[0])*(W[i][1] - wP[1])*(W[i][2] - wP[2])) / _binVolume;
		//std::cout << i << " " << we[i] << std::endl;
		res += pots[7-i]*we[i];
	}

	return res;

}

void capacitor::create_potential_file(std::string ext_pot_path)
{

	double *_T, *_Tnew;
	bool *_TisBoundary;
    int num_gridpoints = _bins[0]*_bins[1]*_bins[2];
	double spacing[3] = {_box[0]/_bins[0],_box[1]/_bins[1],_box[2]/_bins[2]};
	_binVolume = spacing[0]*spacing[1]*spacing[2];
	double minSpacing = sqrt(2) *  min(spacing[0],min(spacing[1],spacing[2]));

    _T =    new double[num_gridpoints]();
    _Tnew = new double[num_gridpoints]();
    _TisBoundary = new bool[num_gridpoints]();

	//Calc surface points and volume distance grid
	std::cout << "Get surface, write distance volume data for each electrode" << std::endl;
    std::ofstream surfaceGridFile;
    surfaceGridFile.open((ext_pot_path + std::string("_surfaceGrid")).c_str());
    for (int i = 0; i < _electrodes.size(); i++)
    {
        std::cout << "Electrode #" << i << std::endl;

        std::ofstream distVolumeGridFile;
        //fname = "./distVolumeGrid" + intToString(i) + ".txt";
        distVolumeGridFile.open((ext_pot_path + std::string("_distMap_el" + intToString(i))).c_str());

        int cnt = 0;
        for (int x = 0; x < _bins[0]; x++)
        {
            std::cout << x*100.0/_bins[0] << "%" << std::endl;
            for (int y = 0; y < _bins[1]; y++)
            {
                for (int z = 0; z < _bins[2]; z++)
                {
                    int G[3] = {x,y,z};
                    double P[3] = {0,0,0};
                    gridToWorld(P, G);

                    double dist = sqrt(_electrodes[i].sqrDistToMesh(P));
                    distVolumeGridFile << P[0] << " " << P[1] << " " << P[2] << " " << dist << "\n";

                    //if (_electrodes[i].sqrDistToMesh(P) <= 0.25)
                    if (dist <= _surface_prec)
                    {
                        _T[cnt] = _electrodes[i].pot;
                        _Tnew[cnt] = _electrodes[i].pot;
						//fprintf(stderr,"SURFACE %d %f\n",cnt,  _electrodes[i].pot);
                        _TisBoundary[cnt] = true;
                        surfaceGridFile << P[0] << " " << P[1] << " " << P[2] << " " << _electrodes[i].pot << "\n";
                    }/* else {
                        _TisBoundary[cnt] = false;
                        _T[cnt] = 0;
                        _Tnew[cnt] = 0;
					}*/
                    cnt++;
                }
            }
        }
        distVolumeGridFile.close();
    }
    surfaceGridFile.close();

	//Calc volume potential
    std::cout << std::endl << "7-Point Stencil Relaxation with boundary values" << std::endl;
	double dMax;
    for (int i = 0; i < _num_iter; i++)
    {
        dMax = 0;
        //std::cout << _T[worldToFlatArrayIndex(new double[3] { 2,2,10})] << std::endl;
        int cnt = 0;
        for (int x = 0; x < _bins[0]; x++)
        {
            for (int y = 0; y < _bins[1]; y++)
            {
                for (int z = 0; z < _bins[2]; z++)
                {
                    int G[3] = {x,y,z};
                    //double* P = gridToWorld(G);
                    //int j = gridToFlatArrayIndex(G);

                    if (!_TisBoundary[cnt])
                    {
                        _Tnew[cnt] = getNeighbourSum(_T, G);
						//if (_Tnew[cnt]<0)
						//	fprintf(stderr,"NEIGHBOUR SUM %f\n",  _Tnew[cnt]);
                    }

                    cnt++;
                }
            }
        }
        cnt = 0;
        for (int x = 0; x < _bins[0]; x++)
        {
            for (int y = 0; y < _bins[1]; y++)
            {
                for (int z = 0; z < _bins[2]; z++)
                {
                    int G[3] =  {x,y,z};
                    //double* P = gridToWorld(G);
                    //int j = gridToFlatArrayIndex(G);

                    if (!_TisBoundary[cnt])
                    {

                        double d =fabs(_T[cnt]-_Tnew[cnt]);
                        if (dMax < d)
                            dMax = d;

                        _T[cnt] = getNeighbourSum(_Tnew, G);
                    }
                    cnt++;
                }
            }
        }
        if (i % 100 == 0)
        {
            std::cout << dMax << std::endl;
            std::cout << "Iteration " << i << std::endl;
        }


        if (dMax < _convergence)
        {
            std::cout << "Convergence after " << i << " Iterations" << std::endl;
            break;
        }

    }
	
	if (dMax > _convergence)
	{
		std::cout << "Reached maximum number of iterations" << std::endl;
	}

	//Save converged volume potential
    std::cout << std::endl << "Save potential mesh to file" << std::endl;
    std::ofstream potVolumeGridFile;
    std::ofstream potVolumeGridFileWCoords;
    potVolumeGridFile.open(ext_pot_path.c_str());
    potVolumeGridFileWCoords.open((ext_pot_path + std::string("_coords")).c_str());
    //potVolumeGridFile << std::setprecision(16) << "1 " << _box[0] << " " << _box[1] << " " << _box[2] << " " << spacing[0] << " " << spacing[1] << " " << spacing[2] << " 0 0 0\n";
    potVolumeGridFile << std::setprecision(16) << "1 " << _box[0] << " " << _box[1] << " " << _box[2] << " " << _bins[0] << " " << _bins[1] << " " << _bins[2] << " 0 0 0\n";
    int cnt = 0;
    for (int x = 0; x < _bins[0]; x++)
    {
        for (int y = 0; y < _bins[1]; y++)
        {
            for (int z = 0; z < _bins[2]; z++)
            {
				//Potential
                int G[3] =  {x,y,z};
                double P[3] = {0,0,0};
                gridToWorld(P, G);

                potVolumeGridFileWCoords << P[0] << " " << P[1] << " " << P[2] << " " << _T[cnt] << "\n";
                potVolumeGridFile << _T[cnt] << "\n";

                cnt++;
            }
        }
    }
    potVolumeGridFile.close();
    potVolumeGridFileWCoords.close();

	//Calculate vacuum surface charge
    std::cout << std::endl << "Save Surface charge to file" << std::endl;
    std::ofstream surfaceChargesFile,normalPointsFile;
    surfaceChargesFile.open((ext_pot_path + std::string("_scharge")).c_str());
    normalPointsFile.open((ext_pot_path + std::string("_normalPoints")).c_str());
	double pref = _eps_0;
	double indChargeSum = 0;
    for (int i = 0; i < _electrodes.size(); i++)
	{
		for (int j = 0; j < _electrodes[i]._numFaces; j++)	
		{
			double nP[3] = {_electrodes[i]._triangles[j].center[0] + _electrodes[i]._triangles[j].normal[0]*(2.0*_surface_prec+minSpacing),
			                _electrodes[i]._triangles[j].center[1] + _electrodes[i]._triangles[j].normal[1]*(2.0*_surface_prec+minSpacing),
			                _electrodes[i]._triangles[j].center[2] + _electrodes[i]._triangles[j].normal[2]*(2.0*_surface_prec+minSpacing)};

			normalPointsFile << nP[0] << " " << nP[1] << " " << nP[2] << std::endl;

            int G[3] =  {0,0,0};
            worldToGrid(nP, G);
			int k = gridToFlatArrayIndex(G);

			double interpolatedPot = interpolatePotOnWorldPoint(_T,nP);

            if (abs(_T[k]-interpolatedPot) > 1)
			{
				std::cout << "Warning: Large difference between potential at next grid point and interpolated potential" << std::endl;
				std::cout << "Pot at next gridpoint: " << _T[k] << std::endl;
				std::cout << "Interpolated Pot: " << interpolatedPot << std::endl;
			}

			if (_TisBoundary[k])
				std::cout << "Warning: Next grid point of surface normal point is boundary" << std::endl;

			double indCharge = (_electrodes[i].pot - interpolatedPot) * _electrodes[i]._triangles[j].area * pref;
			indChargeSum+= indCharge;
            surfaceChargesFile << _electrodes[i]._triangles[j].center[0] << " " << _electrodes[i]._triangles[j].center[1] << " " << _electrodes[i]._triangles[j].center[2] << " " << indCharge << "\n";
		}
	
	}
	surfaceChargesFile.close();
	normalPointsFile.close();
	std::cout << "Total induced charge: " << indChargeSum << std::endl;

	delete[] _T;
	delete[] _Tnew;
	delete[] _TisBoundary; 
}


int setup_capacitor(std::vector<std::string> geofiles, std::vector<double> potentials, std::vector<int> bins, double surface_prec, int num_iter, double convergence, double eps_0, std::string ext_pot_path) 
{
	capacitor cap(geofiles, potentials, bins, surface_prec, num_iter, convergence, eps_0);
	cap.create_potential_file(ext_pot_path);

	return 0;
}
