#pragma once

#include <armadillo>
#include <iostream>
#include <sstream>
#include <iterator>
#include <map>
#include <algorithm>
#include "distribution_func.h"
#include "medium.h"

using namespace arma;
using namespace std;

enum class ProcessType{
	HH_elastic,
	He_elastic,
	Hp_elastic,
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

class PlasmaGasProcess{
public:
	PlasmaGasProcess(const ProcessType proc, const DataType data_type) : pt(proc), dt(data_type){}

	virtual ~PlasmaGasProcess();

	virtual vector<cube> ComputeRightHandSide(const Plasma& p, const DistributionFunction& df) const = 0;

private:
	ProcessType pt;
	DataType dt;
};

class GasGasProcess{
public:
	GasGasProcess(const ProcessType proc, const DataType data_type) : pt(proc), dt(data_type){}

	virtual ~GasGasProcess();

	virtual cube ComputeRightHandSide(const DistributionFunction& df) const = 0;

private:
	ProcessType pt;
	DataType dt;
};

class Hard_spheres_collision : public ElementaryProcess{
public:
	Hard_spheres_collision(const double dcs) : ElementaryProcess(ProcessType::Hard_spheres_collision,
			DataType::Hard_spheres_cross_section), diff_cross_section(dcs) {}

	double ComputeDataTypeValue() const override{
		return diff_cross_section;
	}

private:
	double diff_cross_section;
};

class Charge_exchange : public PlasmaGasProcess{
public:
	Charge_exchange(const string& path) : PlasmaGasProcess(ProcessType::Charge_exchange, DataType::Rate_coefficient) {
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
				param.push_back(*istream_iterator<double>(buffer));
				for(size_t i = 0; i < 8; ++i){
					getline(file, line);
					istringstream tmp_buffer(line);
					param.push_back(*istream_iterator<double>(tmp_buffer));
				}
				fitting_params = join_horiz(fitting_params, vec(param));
			}
		}
	}

	void SetRunoffParam(const Plasma& p, const DistributionFunction& df){
		size_t v_size = df.GetVelGrid().GetSize();
		vec vel_1D(df.GetVelGrid().Get1DGrid());

		runoff_params(df.GetSpaceGrid().GetSize(), cube(v_size,v_size,v_size, fill::zeros));
		double gas_mass = df.GetParticle().mass;
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
			runoff_params[i] *= p.GetDensity(i);
		}
	}

	vector<cube> ComputeRightHandSide(const Plasma& p, const DistributionFunction& df) const override{
		size_t v_size = df.GetVelGrid().GetSize();
		size_t x_size = df.GetSpaceGrid().GetSize();
		vector<cube> rhs(x_size, cube(v_size,v_size,v_size, fill::zeros));
		vector<double> source_params = ComputeSourceParam(p, df);
		for(size_t i = 0; i < x_size; ++i){
			rhs[i] = source_params[i] * p.MakeMaxwellDistr(i, df.GetVelGrid().Get1DGrid())
					- runoff_params[i] * df.GetDistrSlice(i);
		}
		return rhs;
	}

private:
	double ComputeRateCoeff(const double T, const double E) const{
		double result = 0.0;
		if((T_lims.second - T)*(T - T_lims.first) >= 0 and (E_lims.second - E)*(E - E_lims.first) >= 0){
			result = exp(as_scalar(trans(BuildLogVec(T, fitting_params.n_cols))*fitting_params*BuildLogVec(E, fitting_params.n_rows)));
		}else
			throw out_of_range("T or E is out of limits");
		return result;
	}

	vector<double> ComputeSourceParam(const Plasma& p, const DistributionFunction& df) const{
		double phase_volude = pow(df.GetVelGrid().GetGridStep(),3);
		size_t v_size = df.GetVelGrid().GetSize();
		vec vel_1D(df.GetVelGrid().Get1DGrid());

		vector<double> source_params(df.GetSpaceGrid().GetSize(), 0.0);
		double gas_mass = df.GetParticle().mass;
		for(size_t i = 0; i < df.GetSpaceGrid().GetSize(); ++i){
			for(size_t k = 0; k < v_size; ++k){
				double sqr_vz = Sqr(vel_1D(k));
				for(size_t l = 0; l < v_size; ++l){
					double sqr_vy = Sqr(vel_1D(l));
					for(size_t m = 0; m < v_size; ++m){
						double sqr_vx = Sqr(vel_1D(m));
						double energy = gas_mass * (sqr_vz + sqr_vy + sqr_vx) * 0.5;
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
	He_ionization(const string& path) : PlasmaGasProcess(ProcessType::He_Ionization, DataType::Rate_coefficient){
		fstream file(path, ios_base::in);
		if(file.is_open()){
			string line;
			while(getline(file, line)){
				istringstream buffer(line);
				fitting_params.push_back(*istream_iterator<double>(buffer));
			}
		}
	}

	void SetRunoffParams(const Plasma& p){
		runoff_params(p.GetSpaceSize(), 0.0);
		for(size_t i = 0; i < p.GetSpaceSize(); ++i){
			runoff_params[i] = ComputeRateCoeff(p.GetTemperature(i)) * p.GetDensity(i);
		}
	}

	vector<cube> ComputeRightHandSide(const Plasma& p, const DistributionFunction& df) const override{
		size_t v_size = df.GetVelGrid().GetSize();
		vector<cube> rhs(p.GetSpaceSize(), cube(v_size,v_size,v_size, fill::zeros ));
		for(size_t i = 0; i < p.GetSpaceSize(); ++i){
			rhs[i] = - runoff_params[i] * df.GetDistrSlice(i);
		}
		return rhs;
	}

private:
	double ComputeRateCoeff(const double T) const{
		double T_min = fitting_params[0];
		double T_max = fitting_params[1];
		vector<double> logT_vec(fitting_params.size() - 2, 0.0);
		if( (T_max - T)*(T - T_min) >= 0 ){
			logT_vec[0] = log(T);
			for(size_t i = 1; i < logT_vec.size(); ++i){
				logT_vec[i] = logT_vec[i-1] * logT_vec[0];
			}
		}else
			throw out_of_range("Electron temperature in He_ionization is out of range");
		return exp(inner_product(fitting_params.begin() + 2, fitting_params.end(), logT_vec.begin(), 0.0));
	}

	vector<double> fitting_params;
	vector<double> runoff_params;
};

class He_elastic : public PlasmaGasProcess{
public:
	He_elastic(const string& path, const double target_mass_) : PlasmaGasProcess(ProcessType::He_elastic, DataType::Diffusion_coefficient),
	target_mass(target_mass_) {
		fstream file(path, ios_base::in);
		if(file.is_open()){
			string line;
			while(getline(file, line)){
				istringstream buffer(line);
				vector<double> param((istream_iterator<double>(buffer)),
						istream_iterator<double>());
				energies.push_back(param[0]);
				mt_cross_sections.push_back(param[1]);
			}
		}else
			throw logic_error(path + " : file not found");
	}

	void SetDiffusCoeff(const Plasma& p){
		diffus_coeff(p.GetSpaceSize(), 0.0);
		for(size_t i = 0; i < p.GetSpaceSize(); ++i){
			double T = p.GetTemperature(i);
			diffus_coeff[i] = 8.0 / (3 * sqrt(2 * datum::pi)) * datum::c_0 * 100 * p.GetDensity(i) * Sqr(T)
					* T / target_mass / sqrt(T / (datum::m_e * datum::c_0 * datum::c_0 / datum::eV) )
					* ComputeIntOverEnergy(T);
		}
	}

	vector<cube> ComputeRightHandSide(const Plasma& p, const DistributionFunction& df) const override{
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
						rhs[i](m,l,k) = v_x / (sqr_termal_vel * 2 * vel_step) * (f_x_plus - f_x_minus)
								+ v_y / (sqr_termal_vel * 2 * vel_step) * (f_y_plus - f_y_minus)
								+ v_z / (sqr_termal_vel * 2 * vel_step) * (f_z_plus - f_z_minus)
								+ Sqr(datum::c_0 * 100 / vel_step) * (f_x_plus + f_y_plus + f_z_plus + f_x_minus + f_y_minus + f_z_minus);
					}
				}
			}
			rhs[i] += (3 / sqr_termal_vel - 2 * Sqr( datum::c_0 * 100 / vel_step)) * df.GetDistrSlice(i);
			rhs[i] *= diffus_coeff[i] / target_mass;
		}
		return rhs;
	}

private:
	double ComputeIntOverEnergy(const double T) const{
		double result = 0;
		for(size_t i = 0; i < energies.size() - 1; ++i){
			double energy_step = energies[i+1] - energies[i];
			double first_param = mt_cross_sections[i]*(1 + energies[i] / energy_step) - energies[i] / energy_step * mt_cross_sections[i+1];
			first_param /= T;
			double second_param = (mt_cross_sections[i+1] - mt_cross_sections[i]) / energy_step;
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

class HHplus_elastic : public PlasmaGasProcess{
public:
	HHplus_elastic(const string& path) : PlasmaGasProcess(ProcessType::Hp_elastic, DataType::Differential_cross_section){
		fstream file(path, ios_base::in);
		if(file.is_open()){
			string line;
			while(getline(file, line)){
				istringstream buffer(line);
				energies.push_back(*istream_iterator<double>(buffer));
				vector<vector<double>> current_params;
				for(size_t i = 0; i < 3; ++i){
					getline(file, line);
					istringstream tmp_buffer(line);
					vector<double> param((istream_iterator<double>(tmp_buffer)),
							istream_iterator<double>());
					current_params.push_back(param);
				}
				energy_idx_to_coeff[energies.size()-1] = current_params;
			}
		}else
			throw logic_error(path + " : file not found");
	}

	void SetCollisions_mat(const vec& vel_1D);

	double ComputeDataTypeValue() const override{
		auto upper_bound_energy = upper_bound(energies.begin(), energies.end(), E);
		if(upper_bound_energy == energies.end()){
			return FittingFunction(energies.size()-1, angle);
		}
		size_t upper_index = upper_bound_energy - energies.begin();
		if(upper_index == 0){
			return FittingFunction(0, angle);
		}
		double upper_cross_section = FittingFunction(upper_index, angle);
		double lower_cross_section = FittingFunction(upper_index - 1, angle);
		return lower_cross_section
				+ (E - energies[upper_index-1]) / (energies[upper_index] - energies[upper_index-1])
				* (upper_cross_section - lower_cross_section);
	}

private:
	double FittingFunction(const size_t E_index, const double angle) const{
		vector<vector<double>> coeffs = energy_idx_to_coeff.at(E_index);
		double cos_angle = cos(angle);
		double sin_angle = sin(angle);
		double ln_angle = log(angle);
		vector<double> first_log_vec(coeffs[0].size(), 1.0);
		vector<double> second_log_vec(coeffs[1].size(), 1.0);
		for(size_t i = 0; i < max(coeffs[0].size(), coeffs[1].size()); ++i){
			if(i < coeffs[0].size() and i != 0)
				first_log_vec[i] = ln_angle * first_log_vec[i-1];
			if(i < coeffs[1].size()){
				if(i != 0)
					second_log_vec[i] = ln_angle * second_log_vec[i-1];
				else
					second_log_vec[i] = ln_angle;
			}
		}
		return (coeffs[2][0] + coeffs[2][1]*(1 - cos_angle) + coeffs[2][1]*sin_angle*sin_angle)
				/ (2*datum::pi*sin_angle)
				* exp(inner_product(coeffs[0].begin(), coeffs[0].end(), first_log_vec.begin(), 0.0)
						/ inner_product(coeffs[1].begin(), coeffs[1].end(), second_log_vec.begin(), 1.0))
				* datum::a_0 * datum::a_0 * 1e4;
	}

	map<size_t, vector<vector<double>>> energy_idx_to_coeff;
	vector<double> energies;
	mat collisions_mat;
	double E = 0.0;
	double angle = 0.0;
};
