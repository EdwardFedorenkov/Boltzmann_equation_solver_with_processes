#pragma once

#include <armadillo>
#include <memory>
#include "processes.h"

using namespace std;
using namespace arma;

template<typename T>
T Sqr(T x){
	return x*x;
}

enum class ParticlesType{
	H,
	H_2,
	Electron,
	Proton,
	H_2Plus,
	test
};

enum class CollisionModel{
	hard_spheres,
	elastic_from_database,
	excitation_from_database,
	test
};

struct Particle{
	Particle() : p_type(ParticlesType::test),
	vector<shared_ptr<ElementaryProcess>>{make_shared<Hard_spheres_collision>(1.0)},
	mass(Sqr(datum::c_0*100)) {}

	Particle(ParticlesType p, vector<shared_ptr<ElementaryProcess>> pl) : p_type(p), processes_list(pl){
		switch(p){
		case ParticlesType::Electron:
			mass = datum::m_e * datum::c_0 * datum::c_0 / datum::eV;
			break;
		case ParticlesType::H or ParticlesType::Proton:
			mass = datum::m_p * datum::c_0 * datum::c_0 / datum::eV;
			break;
		case ParticlesType::H_2 or ParticlesType::H_2Plus:
			mass = 2 * datum::m_p * datum::c_0 * datum::c_0 / datum::eV;
			break;
		case ParticlesType::test:
			mass = Sqr(datum::c_0*100);
			break;
		}
	}

	double GetMass() const{
		return mass;
	}

	ParticlesType GetParticleType() const{
		return p_type;
	}



	ParticlesType p_type;
	vector<shared_ptr<ElementaryProcess>> processes_list;
	const double mass;
};
