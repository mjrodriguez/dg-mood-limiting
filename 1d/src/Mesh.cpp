#include <stdio.h>
#include <math.h>
#include <cassert>

#include "Mesh.h"

Mesh1d::Mesh1d(int nels, int nn, int nq, double xmin, double xmax){
	assert(nels > 0);
	assert(nn > 1);
	assert(xmin < xmax);

	mNumNodes = nn;
	mNumQuads = nq;
	mNumEls = nels;
	mXmin = xmin;
	mXmax = xmax;

	mDx = (mXmax - mXmin)/(double(mNumEls));
	mJ = mDx;
	mDetJ = std::abs(mJ);
	mInvJ = 1/mJ;

	xnodes = new Nodes(nn);
	xquads = new Nodes(nq);

	mX.alloc(mNumNodes,mNumEls);
	for (int iel = 0; iel < mNumEls; ++iel){
		for (int i = 0; i < mNumNodes; ++i){
			mX(i,iel) = (iel*mDx + mXmin) + xnodes->GetNodes()[i]*mDx;
		}
	}
}

Mesh1d::~Mesh1d(){
	delete xnodes;
	delete xquads;
}

int Mesh1d::nn() const{
	return mNumNodes;
}

int Mesh1d::nq() const{
	return mNumQuads;
}
int Mesh1d::nels() const{
	return mNumEls;
}

double Mesh1d::dx() const{
	return mDx;
}

double Mesh1d::xmax() const {
	return mXmax;
}
double Mesh1d::xmin() const{
	return mXmin;
}

DArray Mesh1d::X() const{
	DArray x(mNumNodes,mNumEls);
	copy(x,mX);
	return x;
}

double Mesh1d::X(int iel, int inode) const{
	return iel*mDx + xnodes->GetNodes()[inode]*mDx;
}

double Mesh1d::J() const{
	return mJ;
}

double Mesh1d::detJ() const {
	return mDetJ;
}

double Mesh1d::invJ() const{{
	return mInvJ;
}}