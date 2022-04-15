#include "distribution_func.h"
#include "elastic_collisions.h"
#include <vector>
#include <numeric>
#include <algorithm>


DistributionFunction::DistributionFunction(const DistributionType distribution_type,
		const double m, const SpaceGrid& x, const VelocityGrid& v,
		const vector<double>& density, const vector<double>& temperature, const vector<vec3>& mean_vel) :
		mass(m), space_grid(x), velocity_grid(v), distribution_function(field<cube>(x.GetSize())) {
	size_t v_size = v.GetSize();
	size_t x_size = distribution_function.n_elem;
	double termal_velocity;
	double termal_vel_factor = sqrt(2.0 / mass)*datum::c_0*100;
	vector<double> vel_1D = v.Get1DGrid();

	for(size_t i = 0; i < x_size; ++i){
		cube slice(v_size, v_size, v_size);
		termal_velocity = sqrt(temperature[i]) * termal_vel_factor;
		double distr_factor = density[i] / pow(sqrt(datum::pi) * termal_velocity, 3);
		for(size_t k = 0; k < v_size; ++k){
			for(size_t l = 0; l < v_size; ++l){
				for(size_t m = 0; m < v_size; ++m){
					if(distribution_type == DistributionType::Maxwell){
						slice(m,l,k) = distr_factor *
								exp(- ( Sqr(vel_1D[k] - mean_vel[i](2))
									  + Sqr(vel_1D[l] - mean_vel[i](1))
									  + Sqr(vel_1D[m] - mean_vel[i](0)) ) / Sqr(termal_velocity) );
					}else if(distribution_type == DistributionType::TestDistribution_1){
						slice(m,l,k) = distr_factor *
								exp(- ( Sqr(vel_1D[k] - mean_vel[i](2))
									  + Sqr(vel_1D[l] - mean_vel[i](1))
									  + Sqr(vel_1D[m] - mean_vel[i](0)) ) / Sqr(termal_velocity) )
									  * Sqr(vel_1D[l] / termal_velocity);
					}
				}
			}
		}
		distribution_function(i) = move(slice);
	}
}

DistributionFunction::DistributionFunction(const DistributionType distribution_type,
		const double m, const SpaceGrid& x, const VelocityGrid& v,
		vector<double>& density, vector<double>& temperature) :
			mass(m), space_grid(x), velocity_grid(v), distribution_function(field<cube>(x.GetSize())){
	size_t v_size = v.GetSize();
	size_t x_size = distribution_function.n_elem;
	double termal_velocity;
	double termal_vel_factor = sqrt(2.0 / mass)*datum::c_0*100;
	vector<double> vel_1D = v.Get1DGrid();

	for(size_t i = 0; i < x_size; ++i){
		cube slice(v_size, v_size, v_size);
		termal_velocity = sqrt(temperature[i]) * termal_vel_factor;
		double distr_factor = density[i] / pow(sqrt(datum::pi) * termal_velocity, 3);
		for(size_t k = 0; k < v_size; ++k){
			for(size_t l = 0; l < v_size; ++l){
				for(size_t m = 0; m < v_size; ++m){
					if(distribution_type == DistributionType::Maxwell){
						slice(m,l,k) = distr_factor *
								exp(- ( Sqr(vel_1D[k])
									  + Sqr(vel_1D[l])
									  + Sqr(vel_1D[m]) ) / Sqr(termal_velocity) );
					}else if(distribution_type == DistributionType::TestDistribution_1){
						slice(m,l,k) = distr_factor *
								exp(- ( Sqr(vel_1D[k])
									  + Sqr(vel_1D[l])
									  + Sqr(vel_1D[m]) ) / Sqr(termal_velocity) )
									  * Sqr(vel_1D[l] / termal_velocity);
					}
				}
			}
		}
		distribution_function(i) = move(slice);
	}
}

DistributionFunction::DistributionFunction(const field<cube>& df, const double m, const SpaceGrid& x,
		const VelocityGrid& v) : mass(m), space_grid(x), velocity_grid(v), distribution_function(df){}

DistributionFunction::DistributionFunction(const DistributionType distribution_type, const double m, const VelocityGrid& v,
		const double density, const double temperature) :
		mass(m), space_grid(), velocity_grid(v){
	if(distribution_type == DistributionType::Maxwell){
		distribution_function = {Maxwell(density, temperature)};
	}else if(distribution_type == DistributionType::TestDistribution_1){
		distribution_function = {TestDistribution_1(density, temperature)};
	}
}

DistributionFunction::DistributionFunction(const cube& df, const double m, const VelocityGrid& v) : mass(m), space_grid(), velocity_grid(v){
	distribution_function = {df};
};

VelocityGrid DistributionFunction::GetVelGrid() const{
	return velocity_grid;
}

SpaceGrid DistributionFunction::GetSpaceGrid() const{
	return space_grid;
}

vector<double> DistributionFunction::ComputeDensity() const{
	vector<double> density;
	density.reserve(distribution_function.n_elem);
	double phase_volume = pow(velocity_grid.GetGridStep(),3);
	for(const cube& item : distribution_function){
		density.push_back(accu(item) * phase_volume);
	}
	return density;
}

double DistributionFunction::ComputeFallingDensity(bool Is_left_wall) const{
	size_t v_size = velocity_grid.GetSize();
	size_t size_of_integr_part = (v_size - 1) / 2;
	double phase_volume = pow(velocity_grid.GetGridStep(),3);
	double falling_density = 0.0;
	if(Is_left_wall){
		falling_density = accu(distribution_function(0)(span(0,size_of_integr_part - 1), span::all, span::all)) * phase_volume;
	}else{
		falling_density = accu(distribution_function(distribution_function.n_elem - 1)(span(size_of_integr_part + 1, v_size - 1), span::all, span::all)) * phase_volume;
	}
	return falling_density;
}

vector<double> DistributionFunction::ComputeTemperature(const vector<vec3>& mean_vel) const{
	size_t x_size = space_grid.GetSize();
	size_t v_size = velocity_grid.GetSize();
	vector<double> vel_1D = velocity_grid.Get1DGrid();
	vector<double> density = ComputeDensity();
	vector<double> temperature(x_size, 0);
	double phase_volume = pow(velocity_grid.GetGridStep(),3);
	for(size_t i = 0; i < x_size; ++i){
		if(density[i] > datum::eps){
			double factor = mass * phase_volume / (3 * density[i] * Sqr(datum::c_0) * 1e4);
			for(size_t k = 0; k < v_size; ++k){
				for(size_t l = 0; l < v_size; ++l){
					for(size_t m = 0; m < v_size; ++m){
						temperature[i] += ( Sqr(vel_1D[k] - mean_vel[i](2)) + Sqr(vel_1D[l] - mean_vel[i](1)) + Sqr(vel_1D[m] - mean_vel[i](0)) ) *
								distribution_function(i)(m,l,k);
					}
				}
			}
			temperature[i] *= factor;
		}
	}
	return temperature;
}

vector<double> DistributionFunction::ComputeTemperatureStatic() const{
	return ComputeTemperature(vector<vec3>(space_grid.GetSize(), {0,0,0}));
}

vector<vec3> DistributionFunction::ComputeMeanVelocity() const{
	size_t x_size = space_grid.GetSize();
	size_t v_size = velocity_grid.GetSize();
	vec vel_1D(velocity_grid.Get1DGrid());
	vector<double> density = ComputeDensity();
	vector<vec3> mean_vel(x_size, vec({0,0,0}));
	double phase_volume = pow(velocity_grid.GetGridStep(),3);

	for(size_t i = 0; i < x_size; ++i){
		if(density[i] > datum::eps){
			for(size_t k = 0; k < v_size; ++k){
				for(size_t l = 0; l < v_size; ++l){
					for(size_t m = 0; m < v_size; ++m){
						mean_vel[i](0) += distribution_function(i)(m,l,k) * vel_1D(m);
						mean_vel[i](1) += distribution_function(i)(m,l,k) * vel_1D(l);
						mean_vel[i](2) += distribution_function(i)(m,l,k) * vel_1D(k);
					}
				}
			}
			mean_vel[i] *= phase_volume / density[i];
		}
	}
	return mean_vel;
}

field<cube> DistributionFunction::GetFullDistribution() const{
	return distribution_function;
}

cube DistributionFunction::GetDistrSlice(const size_t index) const{
	return distribution_function(index);
}

double DistributionFunction::GetParticleMass() const{
	return mass;
}

double DistributionFunction::ComputeTransportTimeStep() const{
	return 0.5 * space_grid.GetGridStep() / velocity_grid.GetMax();
}

void DistributionFunction::Save(const size_t space_position, const string& file_name) const{
	distribution_function(space_position).save(file_name, raw_binary);
}

void DistributionFunction::ChangeDFbyTransport(){
	if(space_grid.GetSize() != 1){
		for(size_t i = 0; i < space_grid.GetSize(); ++i){
			distribution_function(i) += ComputeFlax(i) * ComputeTransportTimeStep() / space_grid.GetGridStep();
		}
	}
}

void DistributionFunction::ChangeDFbyProcess(const vector<cube>& rhs, const double time_step){
	for(size_t i = 0; i < space_grid.GetSize(); ++i){
		distribution_function(i) += rhs[i] * time_step;
	}
}

cube DistributionFunction::Maxwell(double density, double temperature) const{
	size_t v_size = velocity_grid.GetSize();
	vector<double> vel_1D = velocity_grid.Get1DGrid();
	cube maxwell(v_size, v_size, v_size);
	double sqr_termal_vel = 2 * temperature / mass * datum::c_0 * datum::c_0 * 1e4;
	double factor = density / pow(sqrt(datum::pi * sqr_termal_vel), 3);
	for(size_t k = 0; k < v_size; ++k){
		for(size_t l = 0; l < v_size; ++l){
			for(size_t m = 0; m < v_size; ++m){
				maxwell(m,l,k) = factor * exp(- ( Sqr(vel_1D[k]) + Sqr(vel_1D[l]) + Sqr(vel_1D[m]) ) / sqr_termal_vel);
			}
		}
	}
	return maxwell;
}

cube DistributionFunction::TestDistribution_1(double density, double temperature) const{
	size_t v_size = velocity_grid.GetSize();
	vector<double> vel_1D = velocity_grid.Get1DGrid();
	cube distr(v_size, v_size, v_size);
	double sqr_termal_vel = 2 * temperature / mass * datum::c_0 * datum::c_0 * 1e4;
	double factor = density / pow(sqrt(datum::pi * sqr_termal_vel), 3);
	for(size_t k = 0; k < v_size; ++k){
		for(size_t l = 0; l < v_size; ++l){
			for(size_t m = 0; m < v_size; ++m){
				distr(m,l,k) = factor * exp( -( Sqr(vel_1D[k]) + Sqr(vel_1D[l]) + Sqr(vel_1D[m]) ) / sqr_termal_vel ) *
						Sqr(vel_1D[l]) / sqr_termal_vel;
			}
		}
	}
	return distr;
}

cube DistributionFunction::MaxwellReflectedHalf(double density, double temperature, bool Is_left_wall) const{
	size_t v_size = velocity_grid.GetSize();
	size_t mid_index = (v_size - 1) / 2;
	vector<double> vel_1D = velocity_grid.Get1DGrid();
	cube maxwell(v_size, v_size, v_size, fill::zeros);
	double sqr_termal_vel = 2 * temperature / mass * datum::c_0 * datum::c_0 * 1e4;
	double factor = density / pow(sqrt(datum::pi * sqr_termal_vel), 3);
	for(size_t k = 0; k < v_size; ++k){
		for(size_t l = 0; l < v_size; ++l){
			if(Is_left_wall){
				for(size_t m = mid_index + 1; m < v_size; ++m){
					maxwell(m,l,k) = factor * exp(- ( Sqr(vel_1D[k]) + Sqr(vel_1D[l]) + Sqr(vel_1D[m]) ) / sqr_termal_vel);
				}
			}else{
				for(size_t m = 0; m < mid_index; ++m){
					maxwell(m,l,k) = factor * exp(- ( Sqr(vel_1D[k]) + Sqr(vel_1D[l]) + Sqr(vel_1D[m]) ) / sqr_termal_vel);
				}
			}
		}
	}
	return maxwell;
}

cube DistributionFunction::ConstTemperatureWallReflection(bool Is_left_wall) const{
	size_t v_size = velocity_grid.GetSize();
	cube wall_reflection(v_size, v_size, v_size);
	if(Is_left_wall){
		for(size_t k = 0; k < v_size; ++k){
			wall_reflection = MaxwellReflectedHalf(ComputeFallingDensity(Is_left_wall), space_grid.GetWalls().walls_T.first, Is_left_wall);
		}
	}else{
		for(size_t k = 0; k < v_size; ++k){
			wall_reflection = MaxwellReflectedHalf(ComputeFallingDensity(Is_left_wall), space_grid.GetWalls().walls_T.second, Is_left_wall);
		}
	}
	return wall_reflection;
}

cube DistributionFunction::WallPerfectReflection(bool Is_left_wall) const{
	size_t v_size = velocity_grid.GetSize();
	size_t mid_index = (v_size - 1) / 2;
	cube wall_reflection(v_size, v_size, v_size);
	if(Is_left_wall){
		for(size_t k = 0; k < v_size; ++k){
			mat tmp = flipud(distribution_function(0).slice(k));
			wall_reflection.slice(k).rows(mid_index + 1, v_size - 1) =
				tmp.rows(mid_index + 1, v_size - 1);
		}
	}else{
		for(size_t k = 0; k < v_size; ++k){
			mat tmp = flipud(distribution_function(distribution_function.n_elem - 1).slice(k));
			wall_reflection.slice(k).rows(0, mid_index - 1) =
					tmp.rows(0, mid_index - 1);
		}
	}
	return wall_reflection;
}

cube DistributionFunction::FluxbyThreePoints(const cube& df_left, const cube& df_mid, const cube& df_right) const{
	size_t v_size = velocity_grid.GetSize();
	size_t mid_index = (v_size - 1) / 2;
	vec vel_1D(velocity_grid.Get1DGrid());
	cube flux(v_size,v_size,v_size);
	for(size_t k = 0; k < v_size; ++k){
		for(size_t l = 0; l < v_size; ++l){
			for(size_t m = 0; m < v_size; ++m){
				if(m < mid_index){
					flux(m,l,k) = -vel_1D(m) * ( df_right(m,l,k) - df_mid(m,l,k));
				}else if(m > mid_index){
					flux(m,l,k) = vel_1D(m) * ( df_left(m,l,k) - df_mid(m,l,k));
				}
			}
		}
	}
	return flux;
}

cube DistributionFunction::ComputeFlax(size_t slice_index) const{
	cube flux;
	if(slice_index == 0ul){
		cube df_left;
		if(space_grid.GetWalls().walls_BC.first == BC_Type::ConstantTemperatureWall){
			df_left = ConstTemperatureWallReflection(true);
		}else if(space_grid.GetWalls().walls_BC.first == BC_Type::PerfectReflection){
			df_left = WallPerfectReflection(true);
		}
		flux = FluxbyThreePoints(df_left, distribution_function(0), distribution_function(1));
	}else if(slice_index == space_grid.GetSize() - 1){
		cube df_right;
		if(space_grid.GetWalls().walls_BC.second == BC_Type::ConstantTemperatureWall){
			df_right = ConstTemperatureWallReflection(false);
		}else if(space_grid.GetWalls().walls_BC.second == BC_Type::PerfectReflection){
			df_right = WallPerfectReflection(false);
		}
		flux = FluxbyThreePoints(distribution_function(slice_index - 1), distribution_function(slice_index), df_right);
	}else{
		flux = FluxbyThreePoints(distribution_function(slice_index - 1),
				distribution_function(slice_index), distribution_function(slice_index + 1));
	}
	return flux;
}
