#pragma once

using namespace std;

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
