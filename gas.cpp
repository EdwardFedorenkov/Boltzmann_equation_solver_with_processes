#include "gas.h"

Gas::Gas(const DistributionFunction& distr, const vector<shared_ptr<PlasmaGasProcess>>& pg,
		const vector<shared_ptr<GasGasProcess>>& gg) : df(distr), pg_processes(pg), gg_processes(gg){}

double Gas::ComputeTimeStep(const vector<cube>& rhs, const double accuracy) const{
	size_t x_size = df.GetSpaceGrid().GetSize();
	vec time_step(x_size, fill::zeros);
	for(size_t i = 0; i < x_size; ++i){
		double mean_df = mean(mean(mat(mean(df.GetDistrSlice(i),2))));
		cube time_steps = accuracy * df.GetDistrSlice(i) / abs(rhs[i]);
		time_steps.transform([accuracy](double val)
				{ return ((val > datum::eps) and isfinite(val)) ? val : 1; });
		double time_step_condidate = time_steps.min();
		double delta_f = time_step_condidate * mean(mean(mat(mean(abs(rhs[i]),2))));
		if(delta_f < 1e-4*mean_df){
			time_steps = accuracy * df.GetDistrSlice(i) / (abs(rhs[i]) % (rhs[i] < 0));
			time_steps.transform([accuracy](double val)
					{ return ((val > datum::eps) and isfinite(val)) ? val : 1; });
		}
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
	double transport_time_step = df.ComputeTransportTimeStep();
	df.ChangeDFbyTransport();
	df.ChangeDFbyProcess(rhs, proc_time_step);
	return make_pair(transport_time_step, proc_time_step);
}

void Gas::TimeEvolution_ConstTimeStep(const Plasma& p, const double time_step){
	df.ChangeDFbyProcess(ComputeRightHandSides(p), time_step);
}

void Gas::SaveDistr(const size_t space_idx, const size_t time_idx){
	df.Save(space_idx, "DF_z_" + to_string(space_idx) + "_t_" + to_string(time_idx));
}
