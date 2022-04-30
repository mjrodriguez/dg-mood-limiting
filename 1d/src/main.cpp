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
	const int nq  = nn;
	const int nels = 2;

	EulerProblems ep(p,nq,nels,"sod");
	ep.WhichProblem();
	Euler euler(1.4);
	DG dg(ep.mesh, ep.params, &euler);

	// print(ep.Uinit);

	DArray uc(ep.Uinit.size());
	uc = 1.0;

	dg.AssembleElement(uc);



	return 0;
}
