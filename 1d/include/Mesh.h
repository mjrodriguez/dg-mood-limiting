#pragma once

#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include "Array.h"
#include "ArrayUtils.h"
#include "Interpolation.h"

class Mesh1d{
private:
	int mDim = 1;
	int mNumFaces = 2;
	int mNumNodes;
	int mNumQuads;
	int mNumEls;
	double mDx;
	double mXmin, mXmax;

	double mJ;
	double mDetJ;
	double mInvJ;

	Nodes* xnodes;
	Nodes* xquads;

public:
	Mesh1d(int nels, int nn, int nq, double xmin, double xmax);
	~Mesh1d();

	int nn() const;
	int nq() const;
	int nels() const;
	double dx() const;

	double xmax() const;
	double xmin() const;

	double X(int iel, int inode) const;

	double J() const;
	double detJ() const;
	double invJ() const;

};