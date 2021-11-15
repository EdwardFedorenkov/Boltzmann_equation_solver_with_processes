#include "processes.h"

vec BuildLogVec(const double x, const size_t s){
	vec result(s, fill::ones);
	double logx = log(x);
	for(auto& item : result){
		item *= logx;
	}
	return result;
}
