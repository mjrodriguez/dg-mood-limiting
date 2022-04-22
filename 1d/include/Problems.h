#pragma once

#include <set>
#include <vector>
#include <string>
#include <cassert>

#include "Mesh.h"
#include "Array.h"

struct Parameters{
	int neqs, order, nnodes, nquads, nels;
	double courantNumber, maxTime;
	bool limitBool;
	std::vector<double> domain, leftBC, rightBC;
	std::string timeIntegration, riemannSolver, limiterType;
	
	Parameters(int _neqs, int _order, int _nquads, int _nels, std::vector<double> _domain, double _courantNumber, double _maxTime, std::string _timeIntegration, std::string _riemannSolver, bool _limitBool, std::string _limiterType, std::vector<double> _leftBC, std::vector<double> _rightBC){
		
		assert(_domain.size() == 2);
		assert(_domain[0] < _domain[1]);
		assert(_leftBC.size() == 3);
		assert(_rightBC.size() == 3);

		neqs = _neqs;
		order = _order;
		nnodes = _order+1;
		nquads = _nquads;
		nels = _nels;
		domain.push_back(_domain[0]);
		domain.push_back(_domain[1]);
		courantNumber = _courantNumber;
		maxTime = _maxTime;
		timeIntegration = _timeIntegration;
		riemannSolver = _riemannSolver;
		limitBool =_limitBool;
		limiterType = _limiterType;

		for (int i = 0; i < _leftBC.size(); ++i){
			leftBC.push_back(_leftBC[i]);
			rightBC.push_back(_rightBC[i]);
		}


	}

	
};

class EulerProblems{
private:
	int mNeqs = 3;
	int mOrder, mNnodes, mNquads, mNels;
	double mGamma = 1.4;
	std::string mProblem;
	void mSetProblem();
	void mComputeSodIC();
	void mComputeLaxIC();
	void mComputeShuOsherIC();
public:
	Mesh1d *mesh;
	Parameters *params;
	DArray Uinit;
	EulerProblems(int order, int nquads, int nels, std::string problem);
	~EulerProblems();
	void WhichProblem();

};