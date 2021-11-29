#include "processes.h"

vec BuildLogVec(const double x, const size_t s){
	vec result(s, fill::ones);
	double logx = log(x);
	for(size_t i = 1; i < s; ++i){
		result[i] = result[i-1] * logx;
	}
	return result;
}
