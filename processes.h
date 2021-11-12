#pragma once

#include <armadillo>
#include <iostream>
#include <sstream>
#include <iterator>
#include <map>
#include <algorithm>

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
	Charge_exchange
};

enum class DataType{
	Differential_cross_section,
	Rate_coefficient,
	Cross_section_elastic,
	Cross_section_momentum_transport,
	Diffusion_coefficient
};

enum class ParamsType{
	Energy,
	Temperature,
	Angle,
	Density
};

class ElementaryProcess{
public:
	ElementaryProcess(const ProcessType proc, const DataType data_type) : pt(proc), dt(data_type){}

	virtual ~ElementaryProcess();

	virtual double ComputeDataTypeValue(const map<ParamsType, double>& params) const = 0;

private:
	const ProcessType pt;
	const DataType dt;
};

class Charge_exchange : public ElementaryProcess{
public:
	Charge_exchange(const string& path) : ElementaryProcess(ProcessType::Charge_exchange, DataType::Rate_coefficient){

	}

	double ComputeDataTypeValue(const map<ParamsType, double>& params) const override{

	}

private:
	vector<vector<double>> fitting_params;
};

class He_ionization : public ElementaryProcess{
public:
	He_ionization(const string& path) : ElementaryProcess(ProcessType::He_Ionization, DataType::Rate_coefficient){
		fstream file(path, ios_base::in);
		if(file.is_open()){
			string line;
			while(getline(file, line)){
				istringstream buffer(line);
				fitting_params.push_back(*istream_iterator<double>(buffer));
			}
		}
	}

	double ComputeDataTypeValue(const map<ParamsType, double>& params) const override{
		double temperature = params.at(ParamsType::Temperature);
		double T_min = fitting_params[0];
		double T_max = fitting_params[1];
		vector<double> logT_vec(fitting_params.size() - 2, 0.0);
		if(T_max - temperature >= datum::eps and temperature - T_min >= datum::eps){
			logT_vec[0] = log(temperature);
			for(size_t i = 1; i < logT_vec.size(); ++i){
				logT_vec[i] = logT_vec[i-1] * logT_vec[0];
			}
		}else
			throw out_of_range("Electron temperature in He_ionization is out of range");
		return exp(inner_product(fitting_params.begin() + 2, fitting_params.end(), logT_vec.begin(), 0.0));
	}

private:
	vector<double> fitting_params;
};

class He_elastic : public ElementaryProcess{
public:
	He_elastic(const string& path, const double target_mass_) : ElementaryProcess(ProcessType::He_elastic, DataType::Diffusion_coefficient),
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

	double ComputeDataTypeValue(const map<ParamsType, double>& params) const override{
		double temperature = params.at(ParamsType::Temperature);
		double density = params.at(ParamsType::Density);
		return 8.0 / (3 * sqrt(2 * datum::pi)) * datum::c_0 * 100 * density * temperature * temperature
				* temperature / target_mass / sqrt(temperature / (datum::m_e * datum::c_0 * datum::c_0 / datum::eV) )
				* ComputeIntOverEnergy(temperature);
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
};

class HH_elastic : public ElementaryProcess{
public:
	HH_elastic(const ProcessType proc, const string& path) : ElementaryProcess(proc, DataType::Differential_cross_section){
		fstream file(path, ios_base::in);
		if(file.is_open()){
			string line;
			while(getline(file, line)){
				istringstream buffer(line);
				energies.push_back(*istream_iterator<double>(buffer));
				vector<vector<double>> current_params;
				for(size_t i = 0; i < 3; ++i){
					getline(file, line);
					buffer = line;
					vector<double> param((istream_iterator<double>(buffer)),
							istream_iterator<double>());
					current_params.push_back(param);
				}
				energy_idx_to_coeff[energies.size()-1] = current_params;
			}
		}else
			throw logic_error(path + " : file not found");
	}

	double ComputeDataTypeValue(const map<ParamsType, double>& params) const override{
		double energy = params.at(ParamsType::Energy);
		double angle = params.at(ParamsType::Angle);
		auto upper_bound_energy = upper_bound(energies.begin(), energies.end(), energy);
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
				+ (energy - energies[upper_index-1]) / (energy[upper_index] - energies[upper_index-1])
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
};
