#include <iostream>
#include <vector>

#include "Array.h"
#include "ArrayUtils.h"
#include "DGOperators.h"
#include "Mesh.h"
#include "Interpolation.h"
#include "Problems.h"
#include "Equations.h"

int main(int argc, char *argv[]){

	const int p = 4;
	const int nn = p+1;
	const int nq  = 3*nn/2;
	const int nels = 10;

	EulerProblems ep(p,nq,nels,"sod");
	ep.WhichProblem();
	Euler euler(1.4);
	DG dg(ep.mesh, ep.params, &euler);

	dg.AssembleElement(ep.Uinit);

	return 0;
}
