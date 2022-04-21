//TODO(mjrodriguez): Write implementation for EulerProblems
#include <cassert>

#include "Problems.h"

EulerProblems::EulerProblems(int order, int nquads, int nels, std::string problem){
	mProblem = problem;
	mOrder = order;
	mNnodes = mOrder + 1;
	mNquads = nquads;
	mNels = nels;	

	mSetProblem();
}

EulerProblems::~EulerProblems(){
	delete mesh;
	delete params;
}

void EulerProblems::mSetProblem(){
	if (mProblem == "sod"){

		double courantNumber = 0.8;
		double maxTime = 2.0;

		bool limit = true;
		std::string ti = "rk4";
		std::string rs = "llf";
		std::string limtype = "piN";
		std::vector<double> domain{0.0, 1.0};
		std::vector<double> leftBC{1.0, 0.0, 1.0};
		std::vector<double> rightBC{0.125, 0.0, 0.1};

		params = new Parameters(mNeqs, mOrder, mNquads, mNels, domain, courantNumber, maxTime, ti, rs, limit,limtype, leftBC, rightBC);
		mesh = new Mesh1d(mNels, mNnodes, mNquads, domain[0], domain[1]);
	}
	else{
		// std::cout << "Problem not defined. Please define a new problem." << std::endl;
		throw std::runtime_error("Problem not defined. Please define a new problem.");
	}
}

void EulerProblems::WhichProblem(){
	std::cout << "Equations: Euler" << std::endl; 
	std::cout << "Problem: " + mProblem << std::endl;
	std::cout << "Order of Polynomial: " << mOrder << std::endl;
	std::cout << "Number of elements: " << mNels << std::endl;
}