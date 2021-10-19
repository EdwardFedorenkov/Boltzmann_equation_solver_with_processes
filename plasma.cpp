#include "plasma.h"

Plasma::Plasma(const SpaceGrid& sg, const PlasmaParams& i_params, const PlasmaParams& e_params) : space_grid(sg),
		ions(i_params), electrons(e_params){}

Plasma::Plasma(const PlasmaParams& i_params, const PlasmaParams& e_params) : space_grid(),
		ions(i_params), electrons(e_params){}

Plasma::Plasma(const SpaceGrid& sg, PlasmaParams& params) : space_grid(sg), ions(params), electrons(params){}

Plasma::Plasma(PlasmaParams& params) : space_grid(), ions(params), electrons(params){}

cube Plasma::GetDistr(const ChargedParticleName cpn, const size_t space_index, const vec& vel_1D){
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

vector<vec3> Plasma::ComputeCurrentDensity(){
	vector<vec3> result(space_grid.GetSize(), vec3(fill::zeros));
	for(size_t i = 0; i < space_grid.GetSize(); ++i){
		vec3 curr(fill::zeros);
		for(size_t j = 0; j < 3; ++j){
			curr(j) = (ions.density[i] * ions.velocity[i](j) - electrons.density[i] * electrons.velocity[i](j));
		}
		curr *= datum::ec * 3e9;
		result[i] = curr;
	}
	return result;
}
