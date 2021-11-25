//============================================================================
// Name        : Boltzmann_simple_transport.cpp
// Author      : Fedorenkov Eduard
// Version     : 1.0
// Copyright   : BINP Lab 9-0
// Description : Kinetic code for Gas-pasma interaction.
//============================================================================

#include "profile.h"
#include "test_runner.h"
#include "data_saver.h"
#include "distribution_func.h"
#include "elastic_collisions.h"
#include "collisions.h"
#include "processes.h"
#include "medium.h"
#include "gas.h"

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

void SaveProcessData(const Plasma& p, const DistributionFunction& df){
	Charge_exchange cx("ChargeExchange.txt", p, df);
	He_ionization ioniz("ElectronIonization.txt", p);
}

int main() {
	// Testing procedure
	TestRunner tr;
	RUN_TEST(tr, TestVelocityGrid);
	RUN_TEST(tr, TestSpaceGrid);

	// Gas params
	double Tg = 0.1;
	double ng = 1e14;
	double mg = datum::m_p * datum::c_0 * datum::c_0 / datum::eV;

	// Plasma params
	double Tp = 100;
	double np = 1e14;
	double mi = mg;
	Plasma p(mi, Tp, np);

	// Velocity grid params
	size_t v_size = 11;
	double Vmax = 3 * sqrt(2 * Tg / mg);
	VelocityGrid v(v_size, Vmax);

	DistributionFunction f_H(DistributionType::TestDistribution_1, mg, v, ng, Tg);
	vector<shared_ptr<PlasmaGasProcess>> pg_procs = {
			make_shared<He_ionization>("He_ionization.txt", p)
	};
	vector<shared_ptr<GasGasProcess>> gg_procs = {
			make_shared<Hard_spheres_collision>(datum::a_0 * datum::a_0 * 250, v)
	};
	Gas g(f_H, pg_procs, gg_procs);
	auto time_steps = g.TimeEvolution_SmartTimeStep(p, 0.1);
	double dt = time_steps.first + time_steps.second;

	return 0;
}
