#include <stdio.h>
#include <math.h>
#include <cassert>
#include "Interpolation.h"
#include "Array.h"

#define PI 3.14159265

/* ===============================================================
Nodes Members
=============================================================== */


// Default Constructor for Nodes
Nodes::Nodes(int numOfNodes){
    assert(numOfNodes > 0);
    mNumOfNodes = numOfNodes;
    mNodes.alloc(numOfNodes);
    mWeights.alloc(numOfNodes);
    lglnodes();
}

// Default Destructor for Nodes
Nodes::~Nodes(){
    // delete mNodes;
    // delete mWeights;
}

// GetSize() returns the number of Nodes
int Nodes::GetSize() const{
    return mNumOfNodes;
}

// Returns the Nodes
DArray Nodes::GetNodes() const{
    DArray x(mNumOfNodes);
    copy(x,mNodes);
    return x;
}

// Returns the weights for the nodes
DArray Nodes::GetWeights() const{
    DArray x(mNumOfNodes);
    copy(x,mWeights);
    return x;
}

void Nodes::lglnodes(){
    int N = mNumOfNodes;
    int order = N-1;
    double error;
    DArray x(N), xold(N), w(N);
    DArray P(N,N);
    for (int i = 0; i < N; i++){
        x[i]    = cos(PI*i/double(N-1));
        // std::cout << x[i] << std::endl;
        xold[i] = 2.0;
    }

    for (int i=0; i<N; i++){
        error = 1.0;
        while (error > 1e-10){
            xold[i] = x[i];
            P(i,0) = 1.0;
            P(i,1) = x[i];

            for (int k=2; k<N; k++){
                P(i,k) = ( (2.0*double(k) - 1.0)*x[i]*P(i,k-1) - double(k-1)*P(i,k-2) ) / double(k);
            }
            x[i] = xold[i] - ( x[i]*P(i,order) - P(i,order-1) )/ ( double(N)*P(i,order) );
            error = fabs(x[i]-xold[i]);
        }
        w[i] = 2.0/( double(order*N)*P(i,order)*P(i,order) );
    }

    for (int i = 0; i < N; i++){
        mNodes[i]   = 0.5*( x[order-i]+1.0 );
        mWeights[i] = 0.5*w[order-i];
    }
};



/* ===============================================================

Interpolation Members

=============================================================== */

Interpolation::Interpolation(int numOfNodes, int numOfQuads){
    assert(numOfNodes > 0 && numOfQuads >= numOfNodes);
    mNumOfNodes = numOfNodes;
    mNumOfQuads = numOfQuads;

    xnodes = new Nodes(mNumOfNodes);
    xquads = new Nodes(mNumOfQuads);
    G.alloc(mNumOfQuads, mNumOfNodes);
    D.alloc(mNumOfQuads, mNumOfNodes);
    // W = new Matrix(mNumOfQuads, mNumOfQuads);
    W.alloc(mNumOfQuads,mNumOfQuads);
    W = 0.0; D = 0.0; G = 0.0;
    for (int i = 0; i<mNumOfQuads; i++){
        W(i,i) = (*xquads).GetWeights()[i];
    }
    lagint();

}

Interpolation::~Interpolation(){
    delete xnodes;
    delete xquads;
    // delete G;
    // delete D;
    // delete W;
}

int Interpolation::GetSizeNodes() const{
    return mNumOfNodes;
}

int Interpolation::GetSizeQuads() const{
    return mNumOfQuads;
}

DArray Interpolation::GetNodes() const{

    return (*xnodes).GetNodes();
}

DArray Interpolation::GetQuads() const{
    return (*xquads).GetNodes();
}

DArray Interpolation::GetG() const{
    DArray x(mNumOfQuads, mNumOfNodes);
    copy(x,G);
    return x;
}

DArray Interpolation::GetD() const{
    DArray x(mNumOfQuads, mNumOfNodes);
    copy(x,D);
    return x;
}

DArray Interpolation::GetW() const{
    DArray x(mNumOfQuads, mNumOfQuads);
    copy(x,W);
    return x;
}


void Interpolation::lagint() {
    int order = mNumOfNodes-1;
    double ell = 1;
    DArray DD(order+1,order+1);
    DArray w(order+1);
    w = 1.0; DD = 0.0;

    for (int j = 0 ; j < order+1; j++){
        for (int k = 0; k < order+1; k++){
            if (j != k){
                w[j] *= (*xnodes).GetNodes()[j]-(*xnodes).GetNodes()[k];
            }
        }
    }
    for (int j = 0; j < order+1; j++){
        w[j] = 1.0/w[j];
    }

    for (int i = 0; i < order+1; i++){
        for (int j = 0; j < order+1; j++){
            if (i != j) DD(i,j) = w[j]/w[i] / ((*xnodes).GetNodes()[i]-(*xnodes).GetNodes()[j]);
        }
        DD(i,i) = 0.0;
        for (int j = 0; j < order+1; j++){
            if (i != j) DD(i,i) -= DD(i,j);
        }
    }

    for (int i=0; i<mNumOfQuads; i++){
        ell = 1.0;
        for (int j=0; j<mNumOfNodes; j++){
            ell *= GetQuads()[i] - GetNodes()[j];
        }
        for (int j=0; j<mNumOfNodes; j++){
            if (GetQuads()[i] == GetNodes()[j]) {
                G(i,j) = 1.0;
            }
            else {
                G(i,j) = ell*w[j]/( GetQuads()[i] - GetNodes()[j] );
            }
        }
    }

    matmat(D, G, DD);
}
