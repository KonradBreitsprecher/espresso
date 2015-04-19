#ifndef CAPACITOR_H
#define CAPACITOR_H

#include "triangleMesh.hpp"
#include "grid.hpp"
#include <iomanip>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <math.h>

class electrode : public triangleMesh
{
    public:
        electrode();
        electrode(std::string pathToMeshfile, double potential) : triangleMesh(pathToMeshfile) { pot = potential;};
        double pot;
    protected:
    private:
};

class capacitor
{
    public:
		capacitor();
        capacitor(std::vector<std::string> geofiles, std::vector<double> potentials, std::vector<int> bins);
		void create_potential_file(std::string ext_pot_path); 
    protected:
    private:
		void worldToGrid(double worldPoint[3],int gridPoint[3]);
		void gridToWorld(double worldPoint[3], int gridPoint[3]);
		int gridToFlatArrayIndex(int gridPoint[3]);
		int translatedGrid(int* G, int dim, int ds);
		double getNeighbourSum(double* data, int* G);
		
		std::vector<electrode> _electrodes;
		double _box[3], _offset[3], _pref[3];
		int _bins[3];
};

int setup_capacitor(std::vector<std::string> geofiles, std::vector<double> potentials, std::vector<int> bins, std::string ext_pot_path);

#endif // CAPACITOR_H
