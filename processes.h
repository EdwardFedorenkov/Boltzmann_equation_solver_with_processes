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
	Hp_Ionization,
	Charge_exchange
};

enum class DataType{
	Differential_cross_section,
	Rate_coefficient,
	Cross_section
};

pair<string, string> ReadTextParams(const string& file_name){
	ifstream file_params_names (file_name + ".txt");
	string param_1 = "";
	string param_2 = "";
	if (file_params_names.is_open()){
		getline(file_params_names, param_1);
		getline(file_params_names, param_2);
	    file_params_names.close();
	}
	return make_pair(param_1, param_2);
}

template<typename T>
class ElementaryProcess{
public:
	ElementaryProcess(const DataType dt, const ProcessType pt, const string& file_name) : data_type(dt), proc(pt),
	param_names(ReadTextParams(file_name)){
		vec params;
		params.load(file_name + ".bin", raw_binary);

	}

private:
	DataType data_type;
	ProcessType proc;
	pair<string, string> param_names;
	T data;
};

