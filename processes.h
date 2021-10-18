#pragma once

#include <armadillo>

using namespace arma;
using namespace std;

enum class ProcessType{
	HH_Collisions,
	He_Collisions,
	Hp_Collisions,
	H2H_Collisions,
	H2H2_Collisions,
	Recombination,
	H_Ionization,
	H2_Ionization,
	Dissociation,
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

