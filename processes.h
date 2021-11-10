#pragma once

#include <armadillo>
#include <tuple>
#include <fstream>

using namespace arma;
using namespace std;

enum class ProcessType{
	HH_Collisions,
	He_Collisions,
	Hp_Collisions,
	HH_Elastic_and_Excitation,
	He_Excitation,
	Hp_Excitation,
	He_Ionization,
	Charge_exchange
};

enum class DataType{
	Differential_cross_section,
	Rate_coefficient,
	Cross_section_elastic,
	Cross_section_momentum_transport
};

enum class ParamsType{
	Energy,
	Temperature,
	Angle,
	None
};

void fill_params_pair(pair<ParamsType, ParamsType>& ppt, const ParamsType pt, const size_t i){
	if(i == 0)
		ppt.first = pt;
	else
		ppt.second = pt;
}

pair<ParamsType, ParamsType> ReadTextParams(const string& file_name){
	ifstream file_params_names (file_name + ".txt");
	string param_1 = "";
	string param_2 = "";
	if (file_params_names.is_open()){
		getline(file_params_names, param_1);
		getline(file_params_names, param_2);
	    file_params_names.close();
	}
	pair<ParamsType, ParamsType> result = make_pair(ParamsType::None, ParamsType::None);
	for(size_t i = 0; i < 2; ++i){
		switch(i == 0 ? param_1 : param_2){
		case "Energy":
			fill_params_pair(result, ParamsType::Energy, i);
			break;
		case "Temperature":
			fill_params_pair(result, ParamsType::Temperature, i);
			break;
		case "Angle":
			fill_params_pair(result, ParamsType::Angle, i);
			break;
		case "":
			break;
		default:
			throw logic_error("unknown parameter name");
			break;
		}
	}
	return result;
}


class ElementaryProcess{
public:
	ElementaryProcess(const DataType dt, const ProcessType pt, const string& file_name) : data_type(dt), proc(pt),
	param_names(ReadTextParams(file_name)){
		vec params;
		params.load(file_name + ".bin", raw_binary);
	}

	double DiffCrossSectionFitting(const double energy, const double angle) const{
		/*
		 * This part handle HH_Collisions and Hp_Collisions
		 */

	}

	double Rate_coefficient(const double energy, const double temperature) const{

	}


private:
	DataType data_type;
	ProcessType proc;
	pair<ParamsType, ParamsType> param_names;
};
