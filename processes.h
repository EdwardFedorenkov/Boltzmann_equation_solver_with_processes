#pragma once

#include <armadillo>

using namespace arma;
using namespace std;

enum class ProcessType{
	HH_Collisions,
	He_Collisions,
	Hp_Collisions,
	HH_Excitation,
	He_Excitation,
	Hp_Excitation,
	He_Ionization,
	Hp_Ionization,
	Charge_exchange
};

enum class DataDescription{
	Differential_cross_section,
	Rate_coefficient,
	Cross_section
};

class ElementaryProcess{
public:
	ElementaryProcess(const DataDescription dd, const ProcessType pt){

	}
private:
	DataDescription data_type;
	ProcessType proc;
	mat data;
};

