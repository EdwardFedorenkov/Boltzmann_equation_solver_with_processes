#include "gas.h"

Gas::Gas(const DistributionFunction& distr, const vector<shared_ptr<PlasmaGasProcess>>& pg,
		const vector<shared_ptr<GasGasProcess>>& gg) : df(distr), pg_processes(pg), gg_processes(gg){}

double Gas::ComputeTimeStep(const vector<cube>& rhs, const double accuracy) const{
	size_t x_size = df.GetSpaceGrid().GetSize();
	vec time_step(df.GetSpaceGrid().GetSize(), fill::zeros);
	for(size_t i = 0; i < x_size; ++i){
		cube time_steps = accuracy * df.GetDistrSlice(i) / abs(rhs[i]);
		time_steps.transform([accuracy](double val)
				{ return ((val > datum::eps) and isfinite(val)) ? val : 1; });
		time_step(i) = time_steps.min();
	}
	return min(time_step);
}

vector<cube> Gas::ComputeRightHandSides(const Plasma& p) const{
	size_t x_size = df.GetSpaceGrid().GetSize();
	size_t v_size = df.GetVelGrid().GetSize();
	vector<cube> result(df.GetSpaceGrid().GetSize(), cube(v_size,v_size,v_size,fill::zeros));
	for(const auto& pg_proc : pg_processes){
		vector<cube> tmp = pg_proc->ComputePGRightHandSide(p, df);
		for(size_t i = 0; i < x_size; ++i){
			result[i] += tmp[i];
		}
	}
	for(const auto& gg_proc : gg_processes){
		vector<cube> tmp = gg_proc->ComputeGGRightHandSide(df);
		for(size_t i = 0; i < x_size; ++i){
			result[i] += tmp[i];
		}
	}
	return result;
}

pair<double, double> Gas::TimeEvolution_SmartTimeStep(const Plasma& p, const double time_accuracy){
	vector<cube> rhs = ComputeRightHandSides(p);
	double proc_time_step = ComputeTimeStep(rhs, time_accuracy);
	for(cube& c : rhs){
		c *= proc_time_step;
	}
	double transport_time_step = df.ComputeTransportTimeStep();
	df.ChangeDFbyTransport();
	df.ChangeDFbyProcess(rhs);
	return make_pair(transport_time_step, proc_time_step);
}
