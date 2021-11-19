#pragma once

#include "distribution_func.h"
#include "space_grid.h"
#include "particles.h"
#include <armadillo>
#include <memory>

using namespace arma;

class Plasma{
public:
	Plasma(const ParticlesType pt, const vector<double>& T, const vector<double>& n);

	cube MakeMaxwellDistr(const size_t sg_idx, const vec& vel_1D) const;

	size_t GetSpaceSize() const;

	double GetTemperature(const size_t idx) const;

	double GetDensity(const size_t idx) const;

private:
	vector<double> Tp;
	vector<double> np;
	double ion_mass;
};

class Gas{
public:
	Gas(const DistributionFunction& distr, const vector<shared_ptr<PlasmaGasProcess>>& pg,
			const vector<shared_ptr<GasGasProcess>>& gg) : df(distr), pg_processes(pg), gg_processes(gg){}

	double ComputeTimeStep() const;

	cube ComputeRightHandSide() const;

	void TimeEvolution_SmartTimeStep() const;


private:
	DistributionFunction df;
	vector<shared_ptr<PlasmaGasProcess>> pg_processes;
	vector<shared_ptr<GasGasProcess>> gg_processes;
};
