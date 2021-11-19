#include "medium.h"


Plasma::Plasma(const ParticlesType pt, const vector<double>& T, const vector<double>& n) : Tp(T), np(n) {
	if(pt == ParticlesType::H_2Plus)
		ion_mass = 2 * datum::m_p * datum::c_0 * datum::c_0 / datum::eV;
	else if(pt == ParticlesType::Proton)
		ion_mass = datum::m_p * datum::c_0 * datum::c_0 / datum::eV;
	else if(pt == ParticlesType::test)
		ion_mass = datum::c_0 * datum::c_0 / datum::eV;
}

cube Plasma::MakeMaxwellDistr(const size_t sg_idx, const vec& vel_1D) const{
	size_t v_size = vel_1D.size();
	cube distr(v_size,v_size,v_size,fill::zeros);
	double sqr_termal_vel = 2 * Tp[sg_idx] / ion_mass;
	double factor = np[sg_idx] / pow(sqrt(datum::pi * sqr_termal_vel),3);
	for(size_t k = 0; k < v_size; ++k){
		for(size_t l = 0; l < v_size; ++l){
			for(size_t m = 0; m < v_size; ++m){
				distr(m,l,k) =  factor * exp(- (Sqr(vel_1D(m)) + Sqr(vel_1D(k)) + Sqr(vel_1D(l))) / sqr_termal_vel);
			}
		}
	}
	return distr;
}

double Plasma::GetTemperature(const size_t idx) const{
	return Tp[idx];
}

double Plasma::GetDensity(const size_t idx) const{
	return np[idx];
}

size_t Plasma::GetSpaceSize() const{
	return Tp.size();
}
