#include <cassert>

#include "Problems.h"
#include "Interpolation.h"
#include "Equations.h"

EulerProblems::EulerProblems(int order, int nquads, int nels, std::string problem){
	mProblem = problem;
	mOrder = order;
	mNnodes = mOrder + 1;
	mNquads = nquads;
	mNels = nels;	
	Uinit.alloc(mNnodes,mNeqs,mNels);
	mSetProblem();
}

EulerProblems::~EulerProblems(){
	delete mesh;
	delete params;
}

void EulerProblems::mSetProblem(){
	if (mProblem == "sod"){

		double courantNumber = 0.8;
		double maxTime = 1.3;

		bool limit = true;
		std::string ti = "rk4";
		std::string rs = "llf";
		std::string limtype = "piN";
		std::vector<double> domain{0.0, 1.0};
		std::vector<double> leftBC{1.0, 0.0, 1.0};
		std::vector<double> rightBC{0.125, 0.0, 0.1};

		// Initializing parameters, mesh, and initial condition for Sod shock tube
		params = new Parameters(mNeqs, mOrder, mNquads, mNels, domain, courantNumber, maxTime, ti, rs, limit,limtype, leftBC, rightBC);
		mesh = new Mesh1d(mNels, mNnodes, mNquads, domain[0], domain[1]);
		mComputeSodIC();

	}
	else if (mProblem == "lax"){
		double courantNumber = 0.8;
		double maxTime = 0.2;

		bool limit = true;
		std::string ti = "rk4";
		std::string rs = "llf";
		std::string limtype = "piN";
		std::vector<double> domain{-5.0, 5.0};
		std::vector<double> leftBC{0.445, 0.698, 3.528};
		std::vector<double> rightBC{0.5, 0.0, 0.571};

		// Initializing parameters, mesh, and initial condition for Sod shock tube
		params = new Parameters(mNeqs, mOrder, mNquads, mNels, domain, courantNumber, maxTime, ti, rs, limit,limtype, leftBC, rightBC);
		mesh = new Mesh1d(mNels, mNnodes, mNquads, domain[0], domain[1]);
		mComputeLaxIC();
	}
	else if (mProblem == "shu-osher"){
		double courantNumber = 0.8;
		double maxTime = 1.8;

		bool limit = true;
		std::string ti = "rk4";
		std::string rs = "llf";
		std::string limtype = "piN";
		std::vector<double> domain{-5.0, 5.0};
		std::vector<double> leftBC{3.857143, 2.629369, 10.333333};
		std::vector<double> rightBC{1+0.2*std::sin(5*domain[1]), 0.0, 1.0};

		// Initializing parameters, mesh, and initial condition for Sod shock tube
		params = new Parameters(mNeqs, mOrder, mNquads, mNels, domain, courantNumber, maxTime, ti, rs, limit,limtype, leftBC, rightBC);
		mesh = new Mesh1d(mNels, mNnodes, mNquads, domain[0], domain[1]);
		mComputeShuOsherIC();
	}
	else{
		// std::cout << "Problem not defined. Please define a new problem." << std::endl;
		throw std::runtime_error("The problem is not defined. Please define a new problem.");
	}
}

void EulerProblems::mComputeSodIC(){
	double shockLoc = 0.5;
	Euler euler(1.4);

	// Converting left bc from prim to conservative
	DArray pleftBC(params->leftBC.data(),mNeqs);
	DArray cleftBC = euler.Prim2Cons(pleftBC);
	// Converting right bc from prim to conservative
	DArray prightBC(params->rightBC.data(),mNeqs);
	DArray crightBC = euler.Prim2Cons(prightBC);

	for (int iel = 0; iel < mNels; ++iel){
		DArray X = mesh->X();
		double xc = ( X(0,iel) + X(mNnodes-1,iel) )/2;
		// std::cout << "center of element " << iel << " :  " << xc << std::endl;
		if (xc < shockLoc){
			for (int ieq = 0; ieq < params->leftBC.size(); ieq++){
				for (int i = 0; i < mNnodes; i++){
					Uinit(i,ieq,iel) = cleftBC[ieq];
				}
			}
		}
		else if (xc >= shockLoc){
			for (int ieq = 0; ieq < mNeqs; ieq++){
				for (int i = 0; i < mNnodes; i++){
					Uinit(i,ieq,iel) = crightBC[ieq];
				}
			}
		}
	}

}

void EulerProblems::mComputeLaxIC(){
	double shockLoc = 0.0;
	Euler euler(1.4);

	// Converting left bc from prim to conservative
	DArray pleftBC(params->leftBC.data(),mNeqs);
	DArray cleftBC = euler.Prim2Cons(pleftBC);
	// Converting right bc from prim to conservative
	DArray prightBC(params->rightBC.data(),mNeqs);
	DArray crightBC = euler.Prim2Cons(prightBC);

	for (int iel = 0; iel < mNels; ++iel){
		DArray X = mesh->X();
		double xc = ( X(0,iel) + X(mNnodes-1,iel) )/2;
		// std::cout << "center of element " << iel << " :  " << xc << std::endl;
		if (xc < shockLoc){
			for (int ieq = 0; ieq < params->leftBC.size(); ieq++){
				for (int i = 0; i < mNnodes; i++){
					Uinit(i,ieq,iel) = cleftBC[ieq];
				}
			}
		}
		else if (xc >= shockLoc){
			for (int ieq = 0; ieq < mNeqs; ieq++){
				for (int i = 0; i < mNnodes; i++){
					Uinit(i,ieq,iel) = crightBC[ieq];
				}
			}
		}
	}

}

void EulerProblems::mComputeShuOsherIC(){
	double shockLoc = -4.0;
	Euler euler(1.4);

	// Converting left bc from prim to conservative
	DArray pleftBC(params->leftBC.data(),mNeqs);
	DArray cleftBC = euler.Prim2Cons(pleftBC);
	// Converting right bc from prim to conservative
	DArray prightBC(params->rightBC.data(),mNeqs);
	DArray crightBC = euler.Prim2Cons(prightBC);

	for (int iel = 0; iel < mNels; ++iel){
		DArray X = mesh->X();
		double xc = ( X(0,iel) + X(mNnodes-1,iel) )/2;
		// std::cout << "center of element " << iel << " :  " << xc << std::endl;
		if (xc < shockLoc){
			for (int ieq = 0; ieq < params->leftBC.size(); ieq++){
				for (int i = 0; i < mNnodes; i++){
					Uinit(i,ieq,iel) = cleftBC[ieq];
				}
			}
		}
		else if (xc >= shockLoc){
			for (int ieq = 0; ieq < mNeqs; ieq++){
				for (int i = 0; i < mNnodes; i++){
					if (ieq == 0){
						double xn = X(i,iel);
						Uinit(i,ieq,iel) = 1.0 + 0.2*std::sin(5*xn);
					} else {
						Uinit(i,ieq,iel) = crightBC[ieq];
					}
				}
			}
		}
	}

}

void EulerProblems::WhichProblem(){
	std::cout << "Equations: Euler" << std::endl; 
	std::cout << "Problem: " + mProblem << std::endl;
	std::cout << "Order of Polynomial: " << mOrder << std::endl;
	std::cout << "Number of Quadrature Points: " << mNquads << std::endl;
	std::cout << "Number of elements: " << mNels << std::endl;
}