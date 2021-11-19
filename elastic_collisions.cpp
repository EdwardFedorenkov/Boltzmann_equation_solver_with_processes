#include "elastic_collisions.h"

cube ComputeCollisionsIntegral(const cube& df, const mat& coll_mat, const VelocityGrid& v, bool Do_corrections){
	size_t v_size = df.n_slices;
	cube collisions_integral(size(df));
	for(size_t k = 0; k < v_size; ++k){
		for(size_t l = 0; l < v_size; ++l){
			for(size_t m = 0; m < v_size; ++m){
				cube df_shift = DFShiftProcedure(df, {m, l, k});
				collisions_integral(m,l,k) = as_scalar(vectorise(df_shift).t()*coll_mat*vectorise(df_shift));
			}
		}
	}
	if(Do_corrections){
		TreatCollIntegral(collisions_integral, df, v);
	}
	return collisions_integral;
}

mat55 TreatedConservationLawsMatrix(const cube& df, const VelocityGrid& v){
	mat::fixed<5,5> m;
	vec vel_1D = v.Get1DGrid();
	vel_1D /= v.GetMax();
	size_t v_size = v.GetSize();

	m(0,0) = accu(df);
	for(size_t k = 0; k < v_size; ++k){
		for(size_t l = 0; l < v_size; ++l){
			for(size_t i = 0; i < v_size; ++i){
				m(0,1) = m(1,0) += df(i,l,k) * vel_1D(i);
				m(0,2) = m(2,0) += df(i,l,k) * vel_1D(l);
				m(0,3) = m(3,0) += df(i,l,k) * vel_1D(k);
				m(1,1) += df(i,l,k) * Sqr(vel_1D(i));
				m(2,2) += df(i,l,k) * Sqr(vel_1D(l));
				m(3,3) += df(i,l,k) * Sqr(vel_1D(k));

				m(2,1) = m(1,2) += df(i,l,k) * vel_1D(i) * vel_1D(l);
				m(3,1) = m(1,3) += df(i,l,k) * vel_1D(i) * vel_1D(k);
				m(3,2) = m(2,3) += df(i,l,k) * vel_1D(l) * vel_1D(k);

				m(0,4) = m(4,0) += ( Sqr(vel_1D(k)) + Sqr(vel_1D(l)) + Sqr(vel_1D(i)) ) * df(i,l,k);
				m(1,4) = m(4,1) += ( Sqr(vel_1D(k)) + Sqr(vel_1D(l)) + Sqr(vel_1D(i)) ) * df(i,l,k) * vel_1D(i);
				m(2,4) = m(4,2) += ( Sqr(vel_1D(k)) + Sqr(vel_1D(l)) + Sqr(vel_1D(i)) ) * df(i,l,k) * vel_1D(l);
				m(3,4) = m(4,3) += ( Sqr(vel_1D(k)) + Sqr(vel_1D(l)) + Sqr(vel_1D(i)) ) * df(i,l,k) * vel_1D(k);
				m(4,4) += Sqr( Sqr(vel_1D(k)) + Sqr(vel_1D(l)) + Sqr(vel_1D(i)) ) * df(i,l,k);
			}
		}
	}
	//m.for_each([](double& x){ if(abs(x) < datum::eps){x = 0;} });
	return m;
}

vec5 ComputeCollIntegralMoments(const cube& st, const VelocityGrid& v){
	vec5 st_moments;
	size_t v_size = v.GetSize();
	vec vel_1D = v.Get1DGrid();
	vel_1D /= v.GetMax();

	st_moments(0) = accu(st);
	for(size_t k = 0; k < v_size; ++k){
		for(size_t l = 0; l < v_size; ++l){
			for(size_t m = 0; m < v_size; ++m){
				st_moments(1) += st(m,l,k) * vel_1D(m);
				st_moments(2) += st(m,l,k) * vel_1D(l);
				st_moments(3) += st(m,l,k) * vel_1D(k);
				st_moments(4) += ( Sqr(vel_1D(k)) + Sqr(vel_1D(l)) + Sqr(vel_1D(m)) ) * st(m,l,k);
			}
		}
	}
	//st_moments.for_each([](double& x){if(abs(x) < datum::eps){ x = 0;}});
	return -st_moments;
}

void TreatCollIntegral(cube& st, const cube& df, const VelocityGrid& v){
	vec coef = solve(TreatedConservationLawsMatrix(df,v),ComputeCollIntegralMoments(st,v));
	size_t v_size = v.GetSize();
	vec vel_1D = v.Get1DGrid();
	vel_1D /= v.GetMax();
	for(size_t k = 0; k < v_size; ++k){
		for(size_t l = 0; l < v_size; ++l){
			for(size_t m = 0; m < v_size; ++m){
				st(m,l,k) += ( coef(0) + coef(1) * vel_1D(m) + coef(2) * vel_1D(l) + coef(3) * vel_1D(k)
						+ coef(4) * ( Sqr(vel_1D(k)) + Sqr(vel_1D(l)) + Sqr(vel_1D(m)) ) ) * df(m,l,k);
			}
		}
	}
}
