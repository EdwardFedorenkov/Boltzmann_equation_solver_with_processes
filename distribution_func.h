#pragma once

#include "space_grid.h"
#include "velocity_grid.h"
#include "particles.h"
#include <array>
#include <set>
#include <armadillo>

using namespace std;
using namespace arma;


enum class DistributionType{
	Maxwell,              // f ~ exp(-(v - u)^2); v, u - type of vector<double>(3)
	TestDistribution_1,   // f ~ exp(-(v - u1)^2)*u2^2; v, u1, u2 - type of vector<double>(3)
};

class DistributionFunction{
public:
	DistributionFunction(const DistributionType distribution_type,
			const double m, const SpaceGrid& x,
			const VelocityGrid& v,
			const vector<double>& density, const vector<double>& temperature,
			const vector<vec3>& mean_vel);

	DistributionFunction(const DistributionType distribution_type,
			const double m, const SpaceGrid& x,
			const VelocityGrid& v,
			vector<double>& density, vector<double>& temperature);

	DistributionFunction(const field<cube>& df, const double m, const SpaceGrid& x,
			const VelocityGrid& v);

	DistributionFunction(const cube& df, const double m, const VelocityGrid& v);

	DistributionFunction(const DistributionType distribution_type, const double m, const VelocityGrid& v,
			const double density, const double temperature);

	VelocityGrid GetVelGrid() const;

	SpaceGrid GetSpaceGrid() const;

	field<cube> GetFullDistribution() const;

	cube GetDistrSlice(const size_t index) const;

	double GetParticleMass() const;

	vector<double> ComputeDensity() const;

	double ComputeFallingDensity(bool Is_left_wall) const;

	vector<double> ComputeTemperature(const vector<vec3>& first_moment) const;

	vector<double> ComputeTemperatureStatic() const;

	vector<vec3> ComputeMeanVelocity() const;

	double ComputeTransportTimeStep() const;

	void Save(const size_t space_position, const string& file_name) const;

	void ChangeDFbyTransport();

	void ChangeDFbyProcess(const vector<cube>& rhs, const double time_step);

	// ----Can be in private section----

	cube Maxwell(double density, double temperature) const;

	cube TestDistribution_1(double density, double temperature) const;

	cube MaxwellReflectedHalf(double density, double temperature, bool Is_left_wall) const;

	cube ConstTemperatureWallReflection(bool Is_left_wall) const;

	cube WallPerfectReflection(bool Is_left_wall) const;

	cube FluxbyThreePoints(const cube& df_left, const cube& df_mid, const cube& df_right) const;

	cube ComputeFlax(size_t slice_index) const;

private:
	const double mass;             // particle data (mass)
	SpaceGrid space_grid;                // vector size of size_s
	VelocityGrid velocity_grid;          // vector size of size_v
	field<cube> distribution_function;   /* 3D DF preserved in 3D representation
	                                      *size = [size_s X size_v X size_v X size_v].
	                                      */
};
