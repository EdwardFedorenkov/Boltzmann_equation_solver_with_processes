#include "medium.h"

Plasma::Plasma(const double mass, const vector<double>& T, const vector<double>& n) : Tp(T), np(n), ion_mass(mass) {}

Plasma::Plasma(const double mass, const double T, const double n) : Tp(vector<double>{T}), np(vector<double>{n}), ion_mass(mass){}

cube Plasma::MakeMaxwellDistr(const size_t sg_idx, const vec& vel_1D) const{
	size_t v_size = vel_1D.size();
	cube distr(v_size,v_size,v_size,fill::zeros);
	double sqr_termal_vel = 2 * Tp[sg_idx] / ion_mass * datum::c_0 * datum::c_0 * 1e4;
	double factor = np[sg_idx] / pow(sqrt(datum::pi * sqr_termal_vel),3);
	for(size_t k = 0; k < v_size; ++k){
		for(size_t l = 0; l < v_size; ++l){
			for(size_t m = 0; m < v_size; ++m){
				distr(m,l,k) =  factor * exp(- (vel_1D(m)*vel_1D(m) + vel_1D(k)*vel_1D(k) + vel_1D(l)*vel_1D(l)) / sqr_termal_vel);
			}
		}
	}
	return distr;
}

double Plasma::ComputeDencity(const size_t sg_idx, const vec& vel_1D) const{
	return accu(MakeMaxwellDistr(sg_idx, vel_1D)) * pow(abs(vel_1D(1) - vel_1D(0)), 3);
}

double Plasma::GetTemperature(const size_t idx) const{
	return Tp[idx];
}

double Plasma::GetDensity(const size_t idx) const{
	return np[idx];
}

double Plasma::GetTermalVel(const size_t idx) const{
	return sqrt(2 * Tp[idx] / ion_mass) * datum::c_0 * 100;
}

size_t Plasma::GetSpaceSize() const{
	return Tp.size();
}

double Plasma::GetIonMass() const{
	return ion_mass;
}
