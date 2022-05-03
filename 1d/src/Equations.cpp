#include <cmath>

#include "Equations.h"
#include "Array.h"
#include "ArrayUtils.h"

Euler::Euler(double gamma) {
	mGamma = gamma;
}

double Euler::Gamma() const {
	return mGamma;
}

double Euler::Pressure(DArray &state){
	double rhou_sqrd = state[momx]*state[momx];
	double P = (mGamma - 1)*( state[ener] - 0.5*(rhou_sqrd/state[dens]));

	return P;
}

double Euler::MaxVel(DArray &state){
	double p = Pressure(state);
	double sound = std::sqrt( Gamma()*p / state[dens]); 
	double vel = std::sqrt( state[momx]*state[momx]/state[dens]/state[dens] );

	return vel + sound;
}

DArray Euler::Cons2Prim(DArray &consState){
	DArray primState(nstates);
	primState[0] = consState[dens];
	primState[1] = consState[momx]/consState[dens];
	primState[2] = Pressure(consState);

	return primState;
}

DArray Euler::Prim2Cons(DArray &primState){
	DArray consState(nstates);

	consState[dens] = primState[0];
	consState[momx] = primState[0]*primState[1];
	consState[ener] = primState[2]/( Gamma() - 1 ) + 0.5*primState[0]*primState[1]*primState[1];

	return consState;
}

DArray Euler::Flux(DArray &state){
	double rhou_sqrd = state[momx]*state[momx];
	double p = Pressure(state);
	DArray F(nstates);

	F[dens] = state[momx];
	F[momx] = rhou_sqrd/state[dens] + p;
	F[ener] = (state[ener] + p)*state[momx]/state[dens];

	return F;

}

bool Euler::IsStatePhysical(DArray &state){
	if (state[dens] < 0){
		std::cout << "Density is negative." << std::endl;
		return false;
	}

	if (state[ener] < 0){
		std::cout << "Energy is negative. You need some more crystals." << std::endl;
		return false;
	}

	double p = Pressure(state);

	if (p <= 0){
		std::cout << "Pressure non-positive." << std::endl;
		return false;
	}

	return true;
}

int Euler::NS() const{
	return nstates;
}
