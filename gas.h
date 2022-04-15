 #pragma once

#include "distribution_func.h"
#include "processes.h"
#include <armadillo>
#include <memory>

class Gas{
public:
	Gas(const DistributionFunction& distr, const vector<shared_ptr<PlasmaGasProcess>>& pg,
			const vector<shared_ptr<GasGasProcess>>& gg);

	double ComputeTimeStep(const vector<cube>& rhs, const double accuracy) const;

	vector<cube> ComputeRightHandSides(const Plasma& p) const;

	pair<double, double> TimeEvolution_SmartTimeStep(const Plasma& p, const double time_accuracy);

	void TimeEvolution_ConstTimeStep(const Plasma& p, const double time_step);

	void SaveDistr(const size_t space_idx, const size_t time_idx);

private:
	DistributionFunction df;
	vector<shared_ptr<PlasmaGasProcess>> pg_processes;
	vector<shared_ptr<GasGasProcess>> gg_processes;
};
