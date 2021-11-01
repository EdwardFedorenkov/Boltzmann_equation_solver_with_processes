#pragma once

#include <armadillo>
#include <map>

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

template<typename T>
class ElementaryProcess{
public:
	ElementaryProcess(const DataType dd, const ProcessType pt, const string& file_name);

private:
	DataType data_type;
	ProcessType proc;
	T data;
	map<size_t, string> param_to_name;
};

