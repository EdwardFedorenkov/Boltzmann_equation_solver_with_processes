#pragma once

#include <armadillo>
#include <tuple>
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

class He_elastic : public ElementaryProcess{
public:
	He_elastic(const string& path) : ElementaryProcess(ProcessType::He_elastic, DataType::Diffusion_coefficient) {
		fstream file(path, ios_base::in);
		if(file.is_open()){
			string line;
			while(getline(file, line)){
				istringstream buffer(line);
				vector<double> param((istream_iterator<double>(buffer)),
						istream_iterator<double>());
				energies.push_back(param[0]);
				cross_sections.push_back(param[1]);
			}
		}else
			throw logic_error(path + " : file not found");
	}

	double ComputeDataTypeValue(const map<ParamsType, double>& params) const override{
		double temperature = params.at(ParamsType::Temperature);
		double density = params.at(ParamsType::Density);

	}

private:
	vector<double> energies;
	vector<double> cross_sections;
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
