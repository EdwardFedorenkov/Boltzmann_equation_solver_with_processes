#pragma once

#include "velocity_grid.h"
#include "particles.h"
#include "collisions.h"
#include <armadillo>

using namespace std;
using namespace arma;

cube ComputeCollisionsIntegral(const cube& df, const mat& coll_mat, const VelocityGrid& v, bool Do_corrections);

mat55 TreatedConservationLawsMatrix(const cube& df, const VelocityGrid& v);

vec5 ComputeCollIntegralMoments(const cube& st, const VelocityGrid& v);

void TreatCollIntegral(cube& st, const cube& df, const VelocityGrid& v);
