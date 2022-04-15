#include <stdio.h>
#include <math.h>
#include <cassert>

#include "Mesh.h"

Mesh::Mesh(int nels, int nn, int nq, double xmin, double xmax){
	assert(nels > 0);
	assert(nn > 1);
	assert(xmin < xmax);

	mNumNodes = nn;
	mNumQuads = nq;
	mNumEls = nels;
	mXmin = xmin;
	mXmax = xmax;

	mDx = (mXmax - mXmin)/(double(mNumEls));
	xnodes = new Nodes(nn);
	xquads = new Nodes(nq);

}

Mesh::~Mesh(){
	delete xnodes;
	delete xquads;
}

int Mesh::nn() const{
	return mNumNodes;
}

int Mesh::nq() const{
	return mNumQuads;
}
int Mesh::nels() const{
	return mNumEls;
}
