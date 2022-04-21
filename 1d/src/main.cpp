#include <iostream>
#include <vector>

#include "Array.h"
#include "ArrayUtils.h"
#include "DGOperators.h"
#include "Mesh.h"
#include "Interpolation.h"
#include "Problems.h"

int main(int argc, char *argv[]){

	const int dim = 1;
	const int nstates = 3;
	const int p = 4;
	const int nn = p+1;
	const int nq  = 3*nn/2;
	const int nels = 10;
	double xmin = 0.0;
	double xmax = 1.0;

	Mesh1d M(nels, nn, nq, xmin, xmax);
	Operators dg(nn,nq);
	
	DArray uold(nn, nstates, nels);
	uold = 0.0;
	uold(0,0,0) = 1.0;
	uold(0,1,0) = 2.0;
	uold(0,2,0) = 3.0;
	print(uold);

	EulerProblems ep(p,nq,nels,"sod");
	ep.WhichProblem();

	return 0;
}
