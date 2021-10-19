#pragma once
#include "space_grid.h"
#include <armadillo>

using namespace std;
using namespace arma;

enum class ChargedParticleName{
	Proton,
	Electron
};

struct PlasmaParams{
	vector<double> density;
	vector<double> temperature;
	vector<vec3> velocity;
	PlasmaParams(const vector<double>& n, const vector<double>& T, const vector<vec3>& u) :
		density(n), temperature(T), velocity(u){}
	PlasmaParams(const double n, const double T, const vec3 u) :
		density(vector<double>{n}), temperature(vector<double>{T}), velocity(vector<vec3>{u}){}
};

class Plasma{
public:
	Plasma(const SpaceGrid& sg, const PlasmaParams& i_params, const PlasmaParams& e_params);
	Plasma(const PlasmaParams& i_params, const PlasmaParams& e_params);
	Plasma(const SpaceGrid& sg, PlasmaParams& params);
	Plasma(PlasmaParams& params);
	cube GetDistr(const ChargedParticleName cpn, const size_t space_index, const vec& vel_1D);
	vector<vec3> ComputeCurrentDensity();
private:
	SpaceGrid space_grid;
	PlasmaParams ions;
	PlasmaParams electrons;
};
