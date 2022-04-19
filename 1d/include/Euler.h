#pragma once

#include "Array.h"
#include "ArrayUtils.h"

#define dens 0
#define momx 1
#define ener 2

#define velx 1
#define pres 2

class Euler{
private:
	double mGamma;
	int nstates = 3;

public:
	Euler(double gamma);
	double Gamma() const;
	double Pressure(DArray &state);
	double MaxVel(DArray &state);
	DArray Cons2Prim(DArray &consState);
	DArray Prim2Cons(DArray &primState);
	DArray Flux(DArray &state);
	bool IsStatePhysical(DArray &state);
};