#pragma once

#include "velocity_grid.h"
#include "particles.h"
#include <armadillo>

using namespace std;
using namespace arma;

vec3 PolarAxisRotation(const vec3& old_axis, const vec3& new_axis, vec3& point);

vector<vec3> ScatteringSphere(size_t N, const vec3& s_pole, const vec3& n_pole);

pair<vec3, vec3> PostCollisionVelocities(const vec3& initial_vel,
		const vec3& direction, const double reletive_velocity, const double energy_param);

pair<array<size_t, 3>, bool> FindNearestNode(const vec& vel_grid, const vec3& point);

bool IndexCheckRange(const array<int, 3>& idx, size_t size);

cube DFShiftProcedure(const cube& df, const uvec3& position);
