//============================================================================
// Name        : Boltzmann_simple_transport.cpp
// Author      : Fedorenkov Eduard
// Version     : 1.0
// Copyright   : BINP Lab 9-0
// Description : Kinetic code for Gas-pasma interaction.
//============================================================================

#include "profile.h"
#include "test_runner.h"
#include "velocity_grid.h"
#include "space_grid.h"
#include "distribution_func.h"
#include "particles.h"
#include "elastic_collisions.h"

#include <iostream>
#include <armadillo>
#include <vector>
using namespace std;
using namespace arma;

void TestVelocityGrid(){
	VelocityGrid v(11,5);
	ASSERT_EQUAL(v.Get1DGrid(), vector<double>({-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5}));
	ASSERT_EQUAL(v.GetGridStep(), 1);
	ASSERT_EQUAL(v.GetSize(), 11u);
	ASSERT_EQUAL(v.GetMax(), 5);
	VelocityGrid v2(vector<double>({-3.33, -2.22, -1.11, 0, 1.11, 2.22, 3.33}));
	ASSERT_EQUAL(v2.Get1DGrid(), vector<double>({-3.33, -2.22, -1.11, 0, 1.11, 2.22, 3.33}));
	DOUBLES_ASSERT_EQUAL_ACCURATE_TO(datum::eps, v2.GetGridStep(), 1.11);
	ASSERT_EQUAL(v2.GetSize(), 7u);
	ASSERT_EQUAL(v2.GetMax(), 3.33);
}

void TestSpaceGrid(){
	SpaceGrid s;
	ASSERT_EQUAL(s.Get1DGrid(), vector<double>({0}));
	SpaceGrid x(10, 10, make_pair(BC_Type::PerfectReflection, BC_Type::ConstantTemperatureWall), make_pair(0.0, 300 * datum::k_evk));
	ASSERT_EQUAL(x.Get1DGrid(), vector<double>({0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5}));
	ASSERT_EQUAL(x.GetDistance(), 9.5);
	ASSERT_EQUAL(x.GetGridStep(), 1);
	ASSERT_EQUAL(x.GetSize(), 10u);
	ASSERT(x.GetWalls().walls_BC.first == BC_Type::PerfectReflection);
	ASSERT(x.GetWalls().walls_BC.second == BC_Type::ConstantTemperatureWall);
	ASSERT_EQUAL(x.GetWalls().walls_T.first, 0.0);
	ASSERT_EQUAL(x.GetWalls().walls_T.second, 300 * datum::k_evk);
}

void TestParticles(){
	Particle e(ParticlesType::Electron, CollisionModel::test);
	DOUBLES_ASSERT_EQUAL_ACCURATE_TO(1e-5,e.mass, 510999);
	ASSERT(e.particle_type == ParticlesType::Electron);
	ASSERT(e.collision_model == CollisionModel::test);
	Particle h(ParticlesType::H, CollisionModel::hard_spheres);
	Particle h_2(ParticlesType::H_2, CollisionModel::hard_spheres);
	DOUBLES_ASSERT_EQUAL_ACCURATE_TO(datum::eps, 2 * h.mass, h_2.mass);
	DOUBLES_ASSERT_EQUAL_ACCURATE_TO(datum::eps, h.hard_speres_cross_section, h_2.hard_speres_cross_section);
	DOUBLES_ASSERT_EQUAL_ACCURATE_TO(datum::eps, h.hard_speres_cross_section, datum::a_0 * datum::a_0 * 1e4 / 4);
}

void TestDistrFuncSimpleOperations(){
	Particle h(ParticlesType::H, CollisionModel::hard_spheres);
	SpaceGrid x(3, 3, make_pair(BC_Type::PerfectReflection, BC_Type::ConstantTemperatureWall), make_pair(0.0, 1));
	VelocityGrid v(21,10);
	double T = 100.0 / 18 * (h.mass / Sqr(datum::c_0*100));
	vector<double> densitys = {1, 2, 3};
	vector<double> temperatures = {T, T, T};
	vector<vec3> mean_vel = {{0, 0, 0}, {0, 1, 0}, {0, 0, 1}};
	DistributionFunction df(DistributionType::Maxwell, h, x, v, densitys, temperatures, mean_vel);
	ASSERT(df.GetParticle().particle_type == ParticlesType::H);
	ASSERT_EQUAL(df.GetSpaceGrid().Get1DGrid(), vector<double>({0.5, 1.5, 2.5}));
	ASSERT_EQUAL(df.GetSpaceGrid().GetGridStep(), 1);
	ASSERT_EQUAL(df.GetSpaceGrid().GetSize(), 3u);
	ASSERT(df.GetSpaceGrid().GetWalls().walls_BC.first == BC_Type::PerfectReflection);
	ASSERT_EQUAL(df.GetVelGrid().Get1DGrid(), vector<double>({-10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0,
		1, 2, 3, 4, 5, 6, 7, 8, 9, 10}));
	ASSERT_EQUAL(df.GetVelGrid().GetGridStep(), 1);
	ASSERT_EQUAL(df.GetVelGrid().GetMax(), 10);
	ASSERT_EQUAL(df.GetVelGrid().GetSize(), 21u);
	DOUBLES_ASSERT_EQUAL_ACCURATE_TO(datum::eps * 5, df.GetDistrSlice(0)(10,10,10),
			1 / ( pow(sqrt(datum::pi * 2 * T / h.mass) * datum::c_0 *100,3) ) );

	double exp_factor1 = exp( - h.mass / ( 2 * T * Sqr(datum::c_0*100) ) );
	DOUBLES_ASSERT_EQUAL_ACCURATE_TO(datum::eps * 5, df.GetDistrSlice(2)(10,10,10),
			3 / ( pow(sqrt(datum::pi * 2 * T / h.mass) * datum::c_0 *100,3) ) * exp_factor1 );

	double exp_factor2 = exp( - (100 * 2 + 121) * h.mass / ( 2 * T * Sqr(datum::c_0*100) ) );
	DOUBLES_ASSERT_EQUAL_ACCURATE_TO(datum::eps * 100, df.GetDistrSlice(1)(0,0,0),
			2 / ( pow(sqrt(datum::pi * 2 * T / h.mass) * datum::c_0 *100,3) ) * exp_factor2 );

	CONTAINERS_ASSERT_EQUAL_ACCURATE_TO(1e-4, vec(df.ComputeDensity()), vec3({1, 2, 3}));
	CONTAINERS_ASSERT_EQUAL_ACCURATE_TO(1e-3, vec(df.ComputeTemperature(mean_vel)), vec3({T, T, T}));
	//CONTAINERS_ASSERT_EQUAL_ACCURATE_TO(1e-3, vec(df.ComputeMeanVelocity()[0]), vec3({0, 0, 0}));
	//CONTAINERS_ASSERT_EQUAL_ACCURATE_TO(1e-3, vec(df.ComputeMeanVelocity()[1]), vec3({0, 1, 0}));
	//CONTAINERS_ASSERT_EQUAL_ACCURATE_TO(1e-3, vec(df.ComputeMeanVelocity()[2]), vec3({0, 0, 1}));
	DOUBLES_ASSERT_EQUAL_ACCURATE_TO(1e-1, df.GetOneCollisionTime(),
			1 / ( 3 * sqrt(2 * T / h.mass) * datum::c_0 * 100 * h.hard_speres_cross_section * 4 * datum::pi));
}

void TestDistrFuncMoments(){
	cube f(3,3,3);
	VelocityGrid v(3,1);
	Particle p;
	f.slice(0) = {{2.0, 3, 2}, {3, 1, 3}, {2, 3, 2}};
	f.slice(1) = {{4, 4, 4}, {4, 10, 4}, {4, 4, 4}};
	f.slice(2) = f.slice(0);
	DistributionFunction df(f, p, v);
	DOUBLES_ASSERT_EQUAL_ACCURATE_TO(datum::eps, df.ComputeDensity()[0], 84);
	CONTAINERS_ASSERT_EQUAL_ACCURATE_TO(datum::eps, df.ComputeMeanVelocity()[0], vec3({0, 0, 0}));
	DOUBLES_ASSERT_EQUAL_ACCURATE_TO(datum::eps, df.ComputeTemperatureStatic()[0], 146.0 / 252);
	DOUBLES_ASSERT_EQUAL_ACCURATE_TO(datum::eps * 5, df.GetOneCollisionTime(),
			1 / ( 84 * sqrt(2 * 146.0 / 252) * p.hard_speres_cross_section * 4 * datum::pi));
	DOUBLES_ASSERT_EQUAL_ACCURATE_TO(datum::eps, df.ComputeFallingDensity(true), 26);
}

void TestDistrFuncDistrTypes(){
	VelocityGrid v(11,5);
	Particle p;
	double T = 25.0 / 18 * (p.mass / Sqr(datum::c_0*100));
	DistributionFunction df(DistributionType::Maxwell, p, v, 1, T);
	CONTAINERS_ASSERT_EQUAL_ACCURATE_TO(datum::eps * 100, df.Maxwell(1, T), df.GetDistrSlice(0));
	DistributionFunction df_1(DistributionType::TestDistribution_1, p, v, 1, T);
	CONTAINERS_ASSERT_EQUAL_ACCURATE_TO(datum::eps * 100, df_1.TestDistribution_1(1, T), df_1.GetDistrSlice(0));
	CONTAINERS_ASSERT_EQUAL_ACCURATE_TO(datum::eps * 100, df_1.TestDistribution_1(1, T).slice(2), flipud(df_1.GetDistrSlice(0).slice(2)));
	CONTAINERS_ASSERT_EQUAL_ACCURATE_TO(datum::eps * 100, df_1.TestDistribution_1(1, T).slice(2), fliplr(df_1.GetDistrSlice(0).slice(2)));
}

void TestCollisionsIntegral(){
	cube f(3,3,3,fill::ones);
	cube f_shift = DFShiftProcedure(f, {1,1,1});  // have to compare with MatLab variant
	CONTAINERS_ASSERT_EQUAL_ACCURATE_TO(datum::eps, f, f_shift);
	mat coll_mat(27,27,fill::randu);
	cube st = ComputeCollisionsIntegral(f, coll_mat, VelocityGrid(3,1), false);
	DOUBLES_ASSERT_EQUAL_ACCURATE_TO(datum::eps, st(1,1,1), accu(coll_mat));
	DOUBLES_ASSERT_EQUAL_ACCURATE_TO(datum::eps, st(1,1,2), accu(coll_mat(span(0,17), span(0,17))));
}

void TestSimpleTransport(){
	Particle h(ParticlesType::H, CollisionModel::hard_spheres);
	SpaceGrid x(3, 3, make_pair(BC_Type::PerfectReflection, BC_Type::ConstantTemperatureWall), make_pair(0.0, 1));
	VelocityGrid v(21,10);
	double T = 100.0 / 18 * (h.mass / Sqr(datum::c_0*100));
	vector<double> densitys = {1, 2, 3};
	vector<double> temperatures = {T, T, T};
	DistributionFunction df(DistributionType::Maxwell, h, x, v, densitys, temperatures);
	cube maxwell = df.Maxwell(1, T);
	CONTAINERS_ASSERT_EQUAL_ACCURATE_TO(datum::eps * 100, maxwell.slice(7), df.GetDistrSlice(0).slice(7));
	DOUBLES_ASSERT_EQUAL_ACCURATE_TO(1e-5, df.ComputeFallingDensity(true), accu(maxwell(span(0,9), span::all, span::all)) );
	DOUBLES_ASSERT_EQUAL_ACCURATE_TO(1e-5, df.ComputeFallingDensity(false), 3*accu(maxwell(span(0,9), span::all, span::all)) );
	CONTAINERS_ASSERT_EQUAL_ACCURATE_TO(datum::eps, df.WallPerfectReflection(true).slice(5).rows(11, 20),
			df.GetDistrSlice(0).slice(5).rows(11, 20));
	CONTAINERS_ASSERT_EQUAL_ACCURATE_TO(datum::eps, df.WallPerfectReflection(false).slice(5).rows(0, 9),
			df.GetDistrSlice(2).slice(5).rows(0, 9));
}

void RecordPhysParams_LocalCollTest(double temperature, double dencity, double Vmax,
		double time_step, double init_distr_max_val){
	vec params(5, fill::zeros);
	params(0) = temperature;
	params(1) = dencity;
	params(2) = Vmax;
	params(3) = time_step;
	params(4) = init_distr_max_val;
	params.save("params.bin", raw_binary);
}

void LinearDataSaving(mat& data){
	data.save("linear_data.bin", raw_binary);
}

string GetTimeID(const size_t t){
	string result = "";
	vector<size_t> nums;
	size_t i = t;
	while (i != 0){
		nums.push_back(i % 10);
		i /= 10;
	}
	reverse(nums.begin(), nums.end());
	for(size_t j = 0; j < nums.size(); ++j){
		if(j != nums.size() - 1){
			for(size_t k = 0; k < nums[j]; ++k){
				result += "9";
			}
		}else{
			result += to_string(nums[j]);
		}
	}
	return result;
}

int main() {
	// Testing procedure
	TestRunner tr;
	RUN_TEST(tr, TestVelocityGrid);
	RUN_TEST(tr, TestSpaceGrid);
	RUN_TEST(tr, TestParticles);
	RUN_TEST(tr, TestDistrFuncSimpleOperations);
	RUN_TEST(tr, TestDistrFuncMoments);
	RUN_TEST(tr, TestDistrFuncDistrTypes);
	RUN_TEST(tr, TestCollisionsIntegral);
	RUN_TEST(tr, TestSimpleTransport);

	// Creating particles included in model
	Particle H(ParticlesType::H, CollisionModel::hard_spheres);

	// Space and velocity grids parameters
	//size_t x_size = 15;
	//double distance = 100.0;
	size_t v_size = 11;
	//SpaceGrid x(x_size, distance, make_pair(BC_Type::PerfectReflection, BC_Type::PerfectReflection), make_pair(0.0,0.0));

	// Physical problem parameters
	double T = 0.1;   // [eV]
	double n = 5e14;  // [cm^-3]
	//double T_right_wall = 0.1;  // [eV]
	//vec init_temperature = vec(x.Get1DGrid()) * (T_right_wall - T_left_wall) / distance + T_left_wall;         // [eV]
	//vec init_density = vec(x_size, fill::ones) * 5e14;                                                         // [cm^-3]
	//vector<double> init_T(make_move_iterator(init_temperature.begin()), make_move_iterator(init_temperature.end()));
	//vector<double> init_n(make_move_iterator(init_density.begin()), make_move_iterator(init_density.end()));
	double V_max = 2.7 * sqrt(2 * T / H.mass) * datum::c_0 * 100;

	// Velocity grid building
	VelocityGrid v(v_size, V_max);

	// Test Distribution:
	cube df(v_size, v_size, v_size, fill::zeros);
	size_t mid_point = (v_size - 1) / 2;
	df(mid_point, mid_point, mid_point) = n / pow(v.GetGridStep(), 3);

	// Distribution function building
	//DistributionFunction f_H(DistributionType::TestDistribution_1, H, v, n, T);
	DistributionFunction f_H(df, H, v);
	f_H.SaveMatrixes(set<size_t>({0, 1}), (v_size - 1) / 2, 0, "DF0.bin");

	// Time step determination
	double one_collision_time = f_H.GetOneCollisionTime();
	double time_step = 0.05 * one_collision_time;
	size_t time_steps = 10;

	RecordPhysParams_LocalCollTest(T, n, V_max, time_step, max(vectorise(f_H.GetDistrSlice(0))) );

	// Collisions Matrix building
	mat elastic_matrix = BuildCollisionMatrix(v, H);

	for(size_t t = 0; t < time_steps; ++t){
		f_H.RungeKutta2_ElasticCollisons(elastic_matrix, time_step, 0, true);
		/*f_H.RungeKutta2_ElasticCollisons(elastic_matrix, time_step, 0, true);
		dalta_temperature(t) = abs(T0 - f_H.ComputeTemperatureStatic()[0]);
		dalta_density(t) = abs(n0 - f_H.ComputeDensity()[0]);
		mean_velocity_x(t) = f_H.ComputeMeanVelocity()[0](0);
		mean_velocity_y(t) = f_H.ComputeMeanVelocity()[0](1);
		mean_velocity_z(t) = f_H.ComputeMeanVelocity()[0](2);
		*/
		//if( !(t % 10) and t != 0){
			f_H.SaveMatrixes(set<size_t>({0, 1}), (v_size - 1) / 2, 0, "DF" + GetTimeID(t+1) + ".bin");
		//}
	}
	return 0;
}
