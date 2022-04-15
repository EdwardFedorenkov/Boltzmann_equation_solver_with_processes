#pragma once

#include <armadillo>
#include <vector>

using namespace arma;
using namespace std;

class Plasma{
public:
	Plasma(const double mass, const vector<double>& T, const vector<double>& n);

	Plasma(const double mass, const double T, const double n);

	cube MakeMaxwellDistr(const size_t sg_idx, const vec& vel_1D) const;

	double ComputeDencity(const size_t sg_idx, const vec& vel_1D) const;

	size_t GetSpaceSize() const;

	double GetTemperature(const size_t idx) const;

	double GetDensity(const size_t idx) const;

	double GetTermalVel(const size_t idx) const;

	double GetIonMass() const;

private:
	vector<double> Tp;
	vector<double> np;
	double ion_mass;
};
