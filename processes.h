#pragma once

#include <armadillo>
#include <iostream>
#include <sstream>
#include <iterator>
#include <map>
#include <algorithm>
#include "distribution_func.h"
#include "medium.h"
#include "collisions.h"
#include "elastic_collisions.h"

using namespace arma;
using namespace std;

enum class ProcessType{
	HH_elastic,
	He_elastic,
	Hp_elastic,
	HFastIons_elastic,
	HH_Excitation,
	He_Excitation,
	Hp_Excitation,
	He_Ionization,
	Charge_exchange,
	Hard_spheres_collision
};

enum class DataType{
	Differential_cross_section,
	Rate_coefficient,
	Cross_section_elastic,
	Cross_section_momentum_transport,
	Diffusion_coefficient,
	Hard_spheres_cross_section
};

enum class ParamsType{
	Energy,
	Temperature,
	Angle,
	Density
};

vec BuildLogVec(const double x, const size_t s);

double FittingFunction(const vector<vector<double>>& coeffs, double angle, ProcessType pt);

template<typename func>
double LinearDataApprox(const double E, double angle, const vector<double>& energies,
		const vector<vector<vector<double>>>& E_ixd_to_coeffs, func& f, ProcessType pt){
	auto upper_bound_energy = upper_bound(energies.begin(), energies.end(), E);
	if(upper_bound_energy == energies.end()){
		return f(E_ixd_to_coeffs[energies.size()-1], angle, pt);
	}
	size_t upper_index = upper_bound_energy - energies.begin();
	cout << upper_index << endl;
	if(upper_index == 0){
		return f(E_ixd_to_coeffs[0], angle, pt);
	}
	double upper_cross_section = f(E_ixd_to_coeffs[upper_index], angle, pt);
	double lower_cross_section = f(E_ixd_to_coeffs[upper_index - 1], angle, pt);
	cout << upper_cross_section << ' ' << lower_cross_section << endl;
	return lower_cross_section
			+ (E - energies[upper_index-1]) / (energies[upper_index] - energies[upper_index-1])
			* (upper_cross_section - lower_cross_section);
}

pair<vector<double>, vector<vector<vector<double>>>> ReadElasticCrossSec(const string& path);

class PlasmaGasProcess{
public:
	PlasmaGasProcess(const ProcessType proc, const DataType data_type) : pt(proc), dt(data_type){}

	virtual ~PlasmaGasProcess(){}

	virtual vector<cube> ComputePGRightHandSide(const Plasma& p, const DistributionFunction& df) const = 0;

private:
	ProcessType pt;
	DataType dt;
};

class GasGasProcess{
public:
	GasGasProcess(const ProcessType proc, const DataType data_type) : pt(proc), dt(data_type){}

	virtual ~GasGasProcess(){}

	virtual vector<cube> ComputeGGRightHandSide(const DistributionFunction& df) const = 0;

private:
	ProcessType pt;
	DataType dt;
};

// -----------------------------------------
// ----------Plasma-Gas processes-----------
// -----------------------------------------


class Charge_exchange : public PlasmaGasProcess{
public:
	Charge_exchange(const string& path, const Plasma& p, const DistributionFunction& df) :
		PlasmaGasProcess(ProcessType::Charge_exchange, DataType::Rate_coefficient) {
		fstream file(path, ios_base::in);
		if(file.is_open()){
			string line;
			for(size_t j = 0; j < 4; ++j){
				getline(file, line);
				istringstream buffer(line);
				switch(j){
				case 0:
					T_lims.first = *istream_iterator<double>(buffer);
					break;
				case 1:
					T_lims.second = *istream_iterator<double>(buffer);
					break;
				case 2:
					E_lims.first = *istream_iterator<double>(buffer);
					break;
				case 3:
					E_lims.second = *istream_iterator<double>(buffer);
					break;
				}
			}

			while(getline(file, line)){
				istringstream buffer(line);
				vector<double> param(9, 0.0);
				param[0] = *istream_iterator<double>(buffer);
				for(size_t i = 1; i < 9; ++i){
					getline(file, line);
					istringstream tmp_buffer(line);
					param[i] = *istream_iterator<double>(tmp_buffer);
				}
				fitting_params = join_horiz(fitting_params, vec(param));
			}
		}else
			throw logic_error("No file : " + path);
		size_t v_size = df.GetVelGrid().GetSize();
		vec vel_1D(df.GetVelGrid().Get1DGrid());

		runoff_params = vector<cube>(df.GetSpaceGrid().GetSize(), cube(v_size,v_size,v_size, fill::zeros));
		double gas_mass = df.GetParticleMass() * datum::eV * 1e3 / (datum::c_0 * datum::c_0);
		for(size_t i = 0; i < df.GetSpaceGrid().GetSize(); ++i){
			for(size_t k = 0; k < v_size; ++k){
				double sqr_vz = Sqr(vel_1D(k));
				for(size_t l = 0; l < v_size; ++l){
					double sqr_vy = Sqr(vel_1D(l));
					for(size_t m = 0; m < v_size; ++m){
						double sqr_vx = Sqr(vel_1D(m));
						double energy = gas_mass * (sqr_vz + sqr_vy + sqr_vx) * 0.5;
						runoff_params[i](m,l,k) = ComputeRateCoeff(p.GetTemperature(i), energy);
					}
				}
			}
			runoff_params[i] *= p.ComputeDencity(i, vel_1D);
		}
	}

	vector<cube> ComputePGRightHandSide(const Plasma& p, const DistributionFunction& df) const override{
		size_t v_size = df.GetVelGrid().GetSize();
		size_t x_size = df.GetSpaceGrid().GetSize();
		vector<cube> rhs(x_size, cube(v_size,v_size,v_size, fill::zeros));
		vector<double> source_params = ComputeSourceParam(p, df);
		for(size_t i = 0; i < x_size; ++i){
			rhs[i] = source_params[i] * p.MakeMaxwellDistr(i, df.GetVelGrid().Get1DGrid())
					- runoff_params[i] % df.GetDistrSlice(i);
		}
		return rhs;
	}

	void SaveCXRateCoeff(const double E, const size_t Ne) const{
		vec rate_coeff(Ne, fill::zeros);
		vec vecT(Ne, fill::zeros);
		for(size_t i = 0; i < Ne; ++i){
			vecT(i) = T_lims.first + (T_lims.second - T_lims.first)/(Ne-1)*i;
			rate_coeff(i) = ComputeRateCoeff(vecT(i), E);
		}
		rate_coeff.save("CXRateCoeff.bin", raw_binary);
		vecT.save("T_range.bin", raw_binary);
	}

private:
	double ComputeRateCoeff(const double T, const double E) const{
		double result = 0.0;
		if((T_lims.second - T)*(T - T_lims.first) >= 0 and (E_lims.second - E)*(E - E_lims.first) >= 0){
			result = exp(as_scalar(trans(BuildLogVec(T, fitting_params.n_rows))*fitting_params*BuildLogVec(E, fitting_params.n_cols)));
		}else if(E < E_lims.first){
			result = exp(as_scalar(trans(BuildLogVec(T, fitting_params.n_rows))*fitting_params*
					BuildLogVec(E_lims.first, fitting_params.n_cols)));
		}else
			throw out_of_range("T or E is out of limits : T = " + to_string(T) + " E = " + to_string(E));
		return result;
	}

	vector<double> ComputeSourceParam(const Plasma& p, const DistributionFunction& df) const{
		double phase_volude = pow(df.GetVelGrid().GetGridStep(),3);
		size_t v_size = df.GetVelGrid().GetSize();
		vec vel_1D(df.GetVelGrid().Get1DGrid());

		vector<double> source_params(df.GetSpaceGrid().GetSize(), 0.0);
		double gas_mass = df.GetParticleMass();
		for(size_t i = 0; i < df.GetSpaceGrid().GetSize(); ++i){
			for(size_t k = 0; k < v_size; ++k){
				double sqr_vz = Sqr(vel_1D(k));
				for(size_t l = 0; l < v_size; ++l){
					double sqr_vy = Sqr(vel_1D(l));
					for(size_t m = 0; m < v_size; ++m){
						double sqr_vx = Sqr(vel_1D(m));
						double energy = gas_mass * (sqr_vz + sqr_vy + sqr_vx) / (datum::c_0 * datum::c_0 * 1e4) * 0.25;
						source_params[i] += df.GetDistrSlice(i)(m,l,k) * ComputeRateCoeff(p.GetTemperature(i), energy);
					}
				}
			}
			source_params[i] *= phase_volude;
		}
		return source_params;
	}

	mat fitting_params;
	pair<double, double> T_lims;
	pair<double, double> E_lims;
	vector<cube> runoff_params;
};

class He_ionization : public PlasmaGasProcess{
public:
	He_ionization(const string& path, const Plasma& p) : PlasmaGasProcess(ProcessType::He_Ionization, DataType::Rate_coefficient){
		fstream file(path, ios_base::in);
		if(file.is_open()){
			string line;
			while(getline(file, line)){
				istringstream buffer(line);
				fitting_params.push_back(*istream_iterator<double>(buffer));
			}
		}else
			throw logic_error("No file :" + path);
		runoff_params = vector<double>(p.GetSpaceSize(), 0.0);
		for(size_t i = 0; i < p.GetSpaceSize(); ++i){
			runoff_params[i] = ComputeRateCoeff(p.GetTemperature(i)) * p.GetDensity(i);
		}
	}

	vector<cube> ComputePGRightHandSide(const Plasma& p, const DistributionFunction& df) const override{
		size_t v_size = df.GetVelGrid().GetSize();
		vector<cube> rhs(p.GetSpaceSize(), cube(v_size,v_size,v_size, fill::zeros ));
		for(size_t i = 0; i < p.GetSpaceSize(); ++i){
			rhs[i] = - runoff_params[i] * df.GetDistrSlice(i);
		}
		return rhs;
	}

	void SaveIonizRateCoeff(const size_t N) const{
		double T_min = fitting_params[0];
		double T_max = fitting_params[1];
		vec vecT(N, fill::zeros);
		vec rate_coeff(N, fill::zeros);
		for(size_t i = 0; i < N; ++i){
			vecT(i) = T_min + (T_max - T_min)/(N-1)*i;
			rate_coeff(i) = ComputeRateCoeff(vecT(i));
		}
		vecT.save("T_range.bin", raw_binary);
		rate_coeff.save("IonizRateCoeff.bin", raw_binary);
	}

private:
	double ComputeRateCoeff(const double T) const{
		double T_min = fitting_params[0];
		double T_max = fitting_params[1];
		vec logT_vec;
		if( (T_max - T)*(T - T_min) >= 0 ){
			logT_vec = BuildLogVec(T, fitting_params.size() - 2);
		}else
			throw out_of_range("Electron temperature in He_ionization is out of range");
		vec params(fitting_params.size() - 2, fill::zeros);
		for(size_t i = 0; i < params.n_elem; ++i){
			params(i) = fitting_params[i+2];
		}
		return exp(dot(params, logT_vec));
	}

	vector<double> fitting_params;
	vector<double> runoff_params;
};

class He_elastic : public PlasmaGasProcess{
public:
	He_elastic(const string& path, const double target_mass_, const Plasma& p) : PlasmaGasProcess(ProcessType::He_elastic, DataType::Diffusion_coefficient),
	target_mass(target_mass_) {
		fstream file(path, ios_base::in);
		if(file.is_open()){
			string line;
			while(getline(file, line)){
				istringstream buffer(line);
				vector<double> param((istream_iterator<double>(buffer)),
						istream_iterator<double>());
				energies.push_back(param[0]);
				mt_cross_sections.push_back(param[1] / (datum::a_0 * datum::a_0 * 1e4));
			}
		}else
			throw logic_error(path + " : file not found");
		diffus_coeff = vector<double>(p.GetSpaceSize(), 0.0);
		for(size_t i = 0; i < p.GetSpaceSize(); ++i){
			double T = p.GetTemperature(i);
			diffus_coeff[i] = 8.0 / (3 * sqrt(2 * datum::pi)) * datum::c_0 * 100 * p.GetDensity(i)
					* sqrt(datum::m_e * datum::c_0 * datum::c_0 / datum::eV) / target_mass * pow(sqrt(T), 3)
					* (datum::a_0 * datum::a_0 * 1e4) * ComputeIntOverEnergy(T);
		}
	}

	vector<cube> ComputePGRightHandSide(const Plasma& p, const DistributionFunction& df) const override{
		size_t v_size = df.GetVelGrid().GetSize();
		vec vel_1D(df.GetVelGrid().Get1DGrid());
		double vel_step = df.GetVelGrid().GetGridStep();
		vector<cube> rhs(p.GetSpaceSize(), cube(v_size,v_size,v_size, fill::zeros));
		for(size_t i = 0; i < p.GetSpaceSize(); ++i){
			double sqr_termal_vel = p.GetTemperature(i) / target_mass;
			for(size_t k = 0; k < v_size; ++k){
				double v_z = vel_1D(k);
				for(size_t l = 0; l < v_size; ++l){
					double v_y = vel_1D(l);
					for(size_t m = 0; m < v_size; ++m){
						double v_x = vel_1D(m);
						double f_x_plus = (m == v_size - 1) ? 0 : df.GetDistrSlice(i)(m+1,l,k);
						double f_x_minus = (m == 0) ? 0 : df.GetDistrSlice(i)(m-1,l,k);
						double f_y_plus = (l == v_size - 1) ? 0 : df.GetDistrSlice(i)(m,l+1,k);
						double f_y_minus = (l == 0) ? 0 : df.GetDistrSlice(i)(m,l-1,k);
						double f_z_plus = (k == v_size - 1) ? 0 : df.GetDistrSlice(i)(m,l,k+1);
						double f_z_minus = (k == 0) ? 0 : df.GetDistrSlice(i)(m,l,k-1);
						rhs[i](m,l,k) = v_x / (2 * vel_step) * (f_x_plus - f_x_minus)
								+ v_y / (2 * vel_step) * (f_y_plus - f_y_minus)
								+ v_z / (2 * vel_step) * (f_z_plus - f_z_minus)
								+ sqr_termal_vel * Sqr(datum::c_0 * 100 / vel_step) * (f_x_plus + f_y_plus + f_z_plus + f_x_minus + f_y_minus + f_z_minus);
					}
				}
			}
			rhs[i] += (3 - 6 * sqr_termal_vel * Sqr( datum::c_0 * 100 / vel_step)) * df.GetDistrSlice(i);
			rhs[i] *= diffus_coeff[i] / p.GetTemperature(i);
		}
		return rhs;
	}

	double GetDiffusCoeff(const size_t sg_idx){
		return diffus_coeff[sg_idx];
	}

private:
	double ComputeIntOverEnergy(const double T) const{
		double result = 0;
		for(size_t i = 0; i < energies.size() - 1; ++i){
			double energy_step = energies[i+1] - energies[i];
			double first_param = mt_cross_sections[i]*(1 + energies[i] / energy_step) - energies[i] / energy_step * mt_cross_sections[i+1];
			double second_param = T * (mt_cross_sections[i+1] - mt_cross_sections[i]) / energy_step;
			result += IntegralValue(energies[i+1] / T, first_param, second_param)
					- IntegralValue(energies[i] / T, first_param, second_param);
		}
		return result;
	}

	double IntegralValue(const double x, const double param1, const double param2) const{
		return exp(-x) * (-param1 * (x*x + 2*x + 2) - param2 * (x*x*x + 3*x*x + 6*x + 6));
	}

	double target_mass;
	vector<double> energies;
	vector<double> mt_cross_sections;
	vector<double> diffus_coeff;
};

class HFastIons_elastic : public PlasmaGasProcess{
public:
	HFastIons_elastic(const string& path, const Plasma& p,  const VelocityGrid& v) :
		PlasmaGasProcess(ProcessType::HFastIons_elastic, DataType::Differential_cross_section),
		source_coeff(vector<mat>(p.GetSpaceSize(), mat(pow(v.GetSize(),3),pow(v.GetSize(),3),fill::zeros) )),
		runoff_coeff(vector<vec>(p.GetSpaceSize(), vec(pow(v.GetSize(),3),fill::zeros))){
		auto data = ReadElasticCrossSec(path);
		energies = move(data.first);
		energy_idx_to_coeff = move(data.second);
		size_t v_size = v.GetSize();
		size_t N_plasma_vel = 50;
		vector<double> Vel_1D(v.Get1DGrid());
		double phase_volume = pow(v.GetGridStep(), 3);
		for(size_t i = 0; i < p.GetSpaceSize(); i++){
			for(size_t k1 = 0; k1 < v_size; ++k1){
				for(size_t l1 = 0; l1 < v_size; ++l1){
					for(size_t m1 = 0; m1 < v_size; ++m1){
						for(size_t k2 = 0; k2 < v_size; ++k2){
							for(size_t l2 = 0; l2 < v_size; ++l2){
								for(size_t m2 = 0; m2 < v_size; ++m2){
									if(m1 + l1*v_size + k1*v_size*v_size < m2 + l2*v_size + k2*v_size*v_size){
										vec vel_1 = {Vel_1D[m1], Vel_1D[l1], Vel_1D[k1]};
										vec vel_2 = {Vel_1D[m2], Vel_1D[l2], Vel_1D[k2]};
										double RateCoeff = IntegralOverPlasmaVel(N_plasma_vel, i, p, norm(vel_2 - vel_1), false) * phase_volume
												/ (4 * datum::pi * norm(vel_2 - vel_1) * norm(vel_2 - vel_1));
										source_coeff[i](
													m1 + l1*v_size + k1*v_size*v_size,
													m2 + l2*v_size + k2*v_size*v_size
												)
												= RateCoeff;
										if(isnan(RateCoeff)){
											cout << m1 << ' ' << l1 << ' ' << k1 << ' ' << m2 << ' ' << l2 << ' ' << k2 << endl;
											cout << norm(vel_2 - vel_1) << endl;
											IntegralOverPlasmaVel(N_plasma_vel, i, p, norm(vel_2 - vel_1), true);
										}
									}
								}
							}
						}
					}
				}
			}
			source_coeff[i] = symmatu(source_coeff[i]);
			runoff_coeff[i] = sum(source_coeff[i],1);
		}
	}

	vector<cube> ComputePGRightHandSide(const Plasma& p, const DistributionFunction& df) const override{
		size_t v_size = df.GetVelGrid().GetSize();
		vector<cube> rhs(p.GetSpaceSize(), cube(v_size,v_size,v_size, fill::zeros ));
		for(size_t i = 0; i < p.GetSpaceSize(); ++i){
			vec tmp_rhs = source_coeff[i] * vectorise(df.GetDistrSlice(i)) - runoff_coeff[i] % vectorise(df.GetDistrSlice(i));
			for(size_t k = 0; k < v_size; ++k){
				for(size_t l = 0; l < v_size; ++l){
					for(size_t m = 0; m < v_size; ++m){
						rhs[i](m,l,k) = tmp_rhs(m + l*v_size + k*v_size*v_size);
					}
				}
			}
		}
		return rhs;
	}

private:
	double DiffCrossByGasVel(const double E, const double angle) const{
		return LinearDataApprox(E, angle, energies, energy_idx_to_coeff, FittingFunction, ProcessType::HH_elastic);
	}

	double IntegralOverPlasmaVel(const size_t N, const size_t space_idx, const Plasma& p, const double gas_delta_vel,
			bool flag_to_print) const{
		double result = 0.0;
		double x_min = 2;
		double x_max = 8;
		vec x = linspace(x_min, x_max, N);
		double x_step = x(2) - x(1);
		vec plasma_energys = linspace(x_min*p.GetTemperature(space_idx), p.GetTemperature(space_idx)*x_max, N);
		vec plasma_vel = linspace(sqrt(x_min)*p.GetTermalVel(space_idx), sqrt(x_max) * p.GetTermalVel(space_idx), N);
		if(flag_to_print)
			cout << plasma_vel << '\n';
		for(size_t i = 0; i < N-1; ++i){
			double angle = acos(1 - 2 * gas_delta_vel * gas_delta_vel / (plasma_vel(i) * plasma_vel(i)));
			double angle_next = acos(1 - 2 * gas_delta_vel * gas_delta_vel / (plasma_vel(i+1) * plasma_vel(i+1)));
			result += 0.5 * (DiffCrossByGasVel(plasma_energys(i), angle) * x(i) * exp(-x(i)) +
					   DiffCrossByGasVel(plasma_energys(i+1), angle_next) * x(i+1) * exp(-x(i+1)));
			if(flag_to_print){
				cout << 1 - 2 * gas_delta_vel * gas_delta_vel / (plasma_vel(i) * plasma_vel(i)) << '\n';
			}
		}
		return result * 2 * p.GetDensity(space_idx) * p.GetTermalVel(space_idx) * x_step / sqrt(datum::pi);
	}

	vector<vector<vector<double>>> energy_idx_to_coeff;
	vector<double> energies;
	vector<mat> source_coeff;
	vector<vec> runoff_coeff;
};

// ------------------------------------
// ---------- Gas processes -----------
// ------------------------------------

class Hard_spheres_collision : public GasGasProcess{
public:
	Hard_spheres_collision(const double dcs, const VelocityGrid& v) : GasGasProcess(ProcessType::Hard_spheres_collision,
			DataType::Hard_spheres_cross_section), diff_cross_section(dcs) {
		size_t v_size = v.GetSize();
		size_t N_angle = 500;
		double phase_volume = pow(v.GetGridStep(),3);
		vec vel_1D(v.Get1DGrid());
		mat sourse(v_size * v_size * v_size, v_size * v_size * v_size, fill::zeros);
		mat runoff(sourse);
		double factor = phase_volume * 4 * datum::pi / N_angle * diff_cross_section;
		for(size_t k = 0; k < v_size; ++k){
			for(size_t l = 0; l < v_size; ++l){
				for(size_t m = 0; m < v_size; ++m){
					vec3 second_particle_velocity = {vel_1D(m), vel_1D(l), vel_1D(k)};
					double reletive_velocity = sqrt(Sqr(second_particle_velocity(0))
							+ Sqr(second_particle_velocity(1))
							+ Sqr(second_particle_velocity(2)));
					vector<vec3> sphere_points;
					if(array<size_t,3>({m, l, k}) != array<size_t,3>({0, 0, 0})){
						sphere_points = ScatteringSphere(N_angle, second_particle_velocity, {0,0,0});
					}
					for(auto& point : sphere_points){
						pair<vec3, vec3> post_collision_velocities = PostCollisionVelocities(second_particle_velocity, point, reletive_velocity, 0.0);
						auto Node_1 = FindNearestNode(vel_1D, post_collision_velocities.first);
						auto Node_2 = FindNearestNode(vel_1D, post_collision_velocities.second);
						double elem = factor * reletive_velocity;
						if (Node_1.second && Node_2.second){
							sourse(v_size*v_size*Node_1.first[2] + v_size*Node_1.first[1] + Node_1.first[0],
									v_size*v_size*Node_2.first[2] + v_size*Node_2.first[1] + Node_2.first[0]) += elem;
							runoff(v_size*v_size*k + v_size*l + m, (v_size * v_size * v_size - 1) / 2) += elem;
						}
					}
				}
			}
		}
		collisions_mat = move(sourse) - move(runoff);
	}

	vector<cube> ComputeGGRightHandSide(const DistributionFunction& df) const override{
		size_t x_size = df.GetSpaceGrid().GetSize();
		size_t v_size = df.GetVelGrid().GetSize();
		vector<cube> result(x_size, cube(v_size,v_size,v_size,fill::zeros));
		for(size_t i = 0; i < x_size; ++i){
			result[i] = ComputeCollisionsIntegral(df.GetDistrSlice(i), collisions_mat, df.GetVelGrid(), true);
		}
		return result;
	}

private:
	double diff_cross_section;
	mat collisions_mat;
};

class HH_elastic : public GasGasProcess{
public:
	HH_elastic(const string& path, const VelocityGrid& v) : GasGasProcess(ProcessType::HH_elastic, DataType::Differential_cross_section){
		auto p = ReadElasticCrossSec(path);
		energies = move(p.first);
		energy_idx_to_coeff = move(p.second);
		size_t v_size = v.GetSize();
		size_t N_angle = 500;
		double phase_volume = pow(v.GetGridStep(),3);
		vec vel_1D(v.Get1DGrid());
		mat sourse(v_size * v_size * v_size, v_size * v_size * v_size, fill::zeros);
		mat runoff(sourse);
		double factor = phase_volume * 4 * datum::pi / N_angle;
		for(size_t k = 0; k < v_size; ++k){
			for(size_t l = 0; l < v_size; ++l){
				for(size_t m = 0; m < v_size; ++m){
					vec3 second_particle_velocity = {vel_1D(m), vel_1D(l), vel_1D(k)};
					double reletive_velocity = sqrt(Sqr(second_particle_velocity(0))
							+ Sqr(second_particle_velocity(1))
							+ Sqr(second_particle_velocity(2)));
					double CM_energy = datum::m_p * Sqr(reletive_velocity) * 0.25 / datum::eV;
					vector<vec3> sphere_points;
					if(array<size_t,3>({m, l, k}) != array<size_t,3>({0, 0, 0})){
						sphere_points = ScatteringSphere(N_angle, second_particle_velocity, {0,0,0});
					}
					for(auto& point : sphere_points){
						pair<vec3, vec3> post_collision_velocities = PostCollisionVelocities(second_particle_velocity, point, reletive_velocity, 0.0);
						auto Node_1 = FindNearestNode(vel_1D, post_collision_velocities.first);
						auto Node_2 = FindNearestNode(vel_1D, post_collision_velocities.second);
						double elem = factor * reletive_velocity * ComputeDiffCross(CM_energy, acos(norm_dot(point, -second_particle_velocity)));
						if (Node_1.second && Node_2.second){
							sourse(v_size*v_size*Node_1.first[2] + v_size*Node_1.first[1] + Node_1.first[0],
									v_size*v_size*Node_2.first[2] + v_size*Node_2.first[1] + Node_2.first[0]) += elem;
							runoff(v_size*v_size*k + v_size*l + m, (v_size * v_size * v_size - 1) / 2) += elem;
						}
					}
				}
			}
		}
		collisions_mat = move(sourse) - move(runoff);
	}

	vector<cube> ComputeGGRightHandSide(const DistributionFunction& df) const override{
		size_t x_size = df.GetSpaceGrid().GetSize();
		size_t v_size = df.GetVelGrid().GetSize();
		vector<cube> result(x_size, cube(v_size,v_size,v_size,fill::zeros));
		for(size_t i = 0; i < x_size; ++i){
			result[i] = ComputeCollisionsIntegral(df.GetDistrSlice(i), collisions_mat, df.GetVelGrid(), true);
		}
		return result;
	}

	void SaveHHDiffCross(const double E, const size_t N_angle){
		vec angles(N_angle, fill::zeros);
		vec diff_cross(N_angle, fill::zeros);
		for(size_t i = 0; i < N_angle; ++i){
			angles(i) = i * datum::pi / (N_angle - 1);
			diff_cross(i) = ComputeDiffCross(E, angles(i));
 		}
		angles.save("angles_range.bin", raw_binary);
		diff_cross.save("HH_diff_cross.bin", raw_binary);
	}

private:
	double ComputeDiffCross(const double E, const double angle) const{
		return LinearDataApprox(E, angle, energies, energy_idx_to_coeff, FittingFunction, ProcessType::HH_elastic);
	}

	vector<vector<vector<double>>> energy_idx_to_coeff;
	vector<double> energies;
	mat collisions_mat;
};
