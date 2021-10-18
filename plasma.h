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
	Plasma(const SpaceGrid& sg, const PlasmaParams& i_params, const PlasmaParams& e_params) : space_grid(sg),
		ions(i_params), electrons(e_params){}
	Plasma(const PlasmaParams& i_params, const PlasmaParams& e_params) : space_grid(),
		ions(i_params), electrons(e_params){}
	Plasma(const SpaceGrid& sg, PlasmaParams& params) : space_grid(sg), ions(params), electrons(params){}
	Plasma(PlasmaParams& params) : space_grid(), ions(params), electrons(params){}
	cube GetDistr(const ChargedParticleName cpn, const size_t space_index, const vector<double>& vel_1D){
		cube distr(vel_1D.size(), vel_1D.size(), vel_1D.size(), fill::zeros);
		double sqr_termal_vel = cpn == ChargedParticleName::Proton ?
				2 * ions.temperature[space_index] / (datum::m_p * datum::c_0 * datum::c_0 / datum::eV) :
				2 * electrons.temperature[space_index] / (datum::m_e * datum::c_0 * datum::c_0 / datum::eV);
		double factor = cpn == ChargedParticleName::Electron ?
				ions.density[space_index] / pow(sqrt(datum::pi * sqr_termal_vel), 3) :
				electrons.density[space_index] / pow(sqrt(datum::pi * sqr_termal_vel), 3);
		for(size_t k = 0; k < vel_1D.size(); ++k){
			for(size_t l = 0; l < vel_1D.size(); ++l){
				for(size_t m = 0; m < vel_1D.size(); ++m){
					distr(m, l, k) = exp(- (vel_1D[k] * vel_1D[k] + vel_1D[l] * vel_1D[l] + vel_1D[m] * vel_1D[m]) / sqr_termal_vel);
				}
			}
		}
		return distr * factor;
	}
private:
	SpaceGrid space_grid;
	PlasmaParams ions;
	PlasmaParams electrons;
};
