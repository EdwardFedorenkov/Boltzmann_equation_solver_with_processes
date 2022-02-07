#include "processes.h"

vec BuildLogVec(const double x, const size_t s){
	vec result(s, fill::ones);
	double logx = log(x);
	for(size_t i = 1; i < s; ++i){
		result[i] = result[i-1] * logx;
	}
	return result;
}

pair<vector<double>, vector<vector<vector<double>>>> ReadElasticCrossSec(const string& path){
	fstream file(path, ios_base::in);
	vector<vector<vector<double>>> result;
	vector<double> energies;
	if(file.is_open()){
		string line;
		while(getline(file, line)){
			istringstream buffer(line);
			energies.push_back(*istream_iterator<double>(buffer));
			vector<vector<double>> current_params;
			for(size_t i = 0; i < 3; ++i){
				getline(file, line);
				istringstream tmp_buffer(line);
				vector<double> param((istream_iterator<double>(tmp_buffer)),
						istream_iterator<double>());
				current_params.push_back(param);
			}
			result.push_back(current_params);
		}
	}else
		throw logic_error(path + " : file not found");
	return make_pair(energies, result);
}

double FittingFunction(const vector<vector<double>>& coeffs, double angle, ProcessType pt){
	angle = (angle < 1e-5) ? 1e-5 : angle;
	angle = (angle > datum::pi - 1e-4) ? datum::pi - 1e-4 : angle;
	double cos_angle = cos(angle);
	double sin_angle = sin(angle);
	double ln_angle = log(angle);
	vector<double> first_log_vec(coeffs[0].size(), 1.0);
	vector<double> second_log_vec(coeffs[1].size(), 1.0);
	for(size_t i = 0; i < max(coeffs[0].size(), coeffs[1].size()); ++i){
		if(i < coeffs[0].size() and i != 0)
			first_log_vec[i] = ln_angle * first_log_vec[i-1];
		if(i < coeffs[1].size()){
			second_log_vec[i] = i != 0 ? ln_angle * second_log_vec[i-1] : ln_angle;
		}
	}
	if(pt == ProcessType::HFastIons_elastic){
		return (coeffs[2][0] + coeffs[2][1]*(1 - cos_angle) + coeffs[2][1]*sin_angle*sin_angle)
				* exp(inner_product(coeffs[0].begin(), coeffs[0].end(), first_log_vec.begin(), 0.0)
				/ inner_product(coeffs[1].begin(), coeffs[1].end(), second_log_vec.begin(), 1.0));
	}else if(pt == ProcessType::HH_elastic){
		return (coeffs[2][0] + coeffs[2][1]*(1 - cos_angle) + coeffs[2][1]*sin_angle*sin_angle)
				/ (2*datum::pi*sin_angle)
				* exp(inner_product(coeffs[0].begin(), coeffs[0].end(), first_log_vec.begin(), 0.0)
						/ inner_product(coeffs[1].begin(), coeffs[1].end(), second_log_vec.begin(), 1.0))
				* datum::a_0 * datum::a_0 * 1e4;
	}
	return 0.0;
}
