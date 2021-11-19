#include "collisions.h"

vec3 PolarAxisRotation(const vec3& old_axis, const vec3& new_axis, vec3& point){
	vec3 ax = normalise(cross(old_axis, new_axis));
	if(norm(ax) >= datum::eps){
		double cos_theta = norm_dot(new_axis, old_axis);
		double theta = acos(cos_theta);
		mat33 M;
		M(0,0) = cos_theta + (1.0 - cos_theta)*Sqr(ax[0]);
		M(0,1) = (1.0 - cos_theta)*ax[0]*ax[1] - sin(theta)*ax[2];
		M(0,2) = (1.0 - cos_theta)*ax[0]*ax[2] + sin(theta)*ax[1];
		M(1,0) = (1.0 - cos_theta)*ax[0]*ax[1] + sin(theta)*ax[2];
		M(1,1) = cos_theta + (1.0 - cos_theta)*Sqr(ax[1]);
		M(1,2) = (1.0 - cos_theta)*ax[1]*ax[2] - sin(theta)*ax[0];
		M(2,0) = (1.0 - cos_theta)*ax[0]*ax[2] - sin(theta)*ax[1];
		M(2,1) = (1.0 - cos_theta)*ax[1]*ax[2] + sin(theta)*ax[0];
		M(2,2) = cos_theta + (1.0 - cos_theta)*Sqr(ax[2]);
		point = M * point;
	}
	return point;
}

vector<vec3> ScatteringSphere(size_t N, const vec3& s_pole, const vec3& n_pole){
	vector<vec3> sphere_points;
	sphere_points.reserve(N);
	double gr = (1.0 + sqrt(5)) / 2.0;
	double ga = 2.0*datum::pi*(1.0 - 1.0 / gr);

	vec3 e_z = {0,0,1};
	vec3 spiral_axis = s_pole - n_pole;
	for(size_t i = 0; i < N; ++i){
		vec3 point;
		point(0) = sin(acos(1.0 - i * (2.0 / (N-1))))*cos(i*ga);
		point(1) = sin(acos(1.0 - i * (2.0 / (N-1))))*sin(i*ga);
		point(2) = 1.0 - i * (2.0 / (N-1));
		PolarAxisRotation(e_z, spiral_axis, point);
		sphere_points.push_back(point);
	}
	return sphere_points;
}

pair<vec3, vec3> PostCollisionVelocities(const vec3& initial_vel,
		const vec3& direction, const double reletive_velocity, const double energy_param){
	vec3 velocity_out_1, velocity_out_2;
	velocity_out_1 = (direction*sqrt(Sqr(reletive_velocity) - 4*Sqr(datum::c_0*100)*energy_param) + initial_vel)*0.5;
	velocity_out_2 = (-direction*sqrt(Sqr(reletive_velocity) - 4*Sqr(datum::c_0*100)*energy_param) + initial_vel)*0.5;
	return make_pair(velocity_out_1,velocity_out_2);
}

pair<array<size_t, 3>, bool> FindNearestNode(const vec& vel_grid, const vec3& point){
	double grid_step = vel_grid(1) - vel_grid(0);
	array<size_t, 3> result = {0, 0, 0};
	for(size_t i = 0; i < 3; ++i){
		result[i] = round((point(i) - vel_grid(0)) / grid_step);
		if(result[i] < -0 || result[i] >= vel_grid.size()){
			return make_pair(result, false);
		}
	}
	return make_pair(result, true);
}

bool IndexCheckRange(const array<int, 3>& idx, size_t size){
	bool result;
	for(const int item : idx){
		if(item >= 0 && static_cast<size_t>(item)< size){
			result = true;
		}else{
			return false;
		}
	}
	return result;
}

cube DFShiftProcedure(const cube& df, const uvec3& position){
	size_t v_size = df.n_slices;
	size_t mid_point = (v_size - 1) / 2;
	cube result(size(df));
	result.fill(0);
	array<int, 3> slice_index = {0,0,0};
	for(size_t k = 0; k < v_size; ++k){
		for(size_t l = 0; l < v_size; ++l){
			for(size_t m = 0; m < v_size; ++m){
				slice_index[0] = m + position(0) - mid_point;
				slice_index[1] = l + position(1) - mid_point;
				slice_index[2] = k + position(2) - mid_point;
				if(IndexCheckRange(slice_index, v_size)){
					result(m,l,k) = df(slice_index[0], slice_index[1], slice_index[2]);
				}
			}
		}
	}
	return result;
}
