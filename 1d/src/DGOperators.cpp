#include <cassert>
#include "Array.h"
#include "ArrayUtils.h"
#include "Interpolation.h"
#include "DGOperators.h"

Operators::Operators(int numNodes, int numQuads){
    assert(numNodes > 0 && numQuads >= numNodes);
    mNumOfNodes = numNodes;
    mNumOfQuads = numQuads;

    I = new Interpolation(mNumOfNodes, mNumOfQuads);
    mass.alloc(mNumOfNodes, mNumOfNodes);
    massChol.alloc(mNumOfNodes, mNumOfNodes);
    volumeTerm.alloc(mNumOfNodes,mNumOfQuads);
    invMassVol.alloc(mNumOfNodes,mNumOfQuads);
    ComputeMass();
    ComputeVol();
    ComputeInvMassVol();

}

Operators::~Operators(){
    delete I;
}

void Operators::ComputeMass(){
    DArray work(mNumOfQuads,mNumOfNodes); work = 0.0;
    matmat(work,I->GetW(),I->GetG());
    matmat(mass,I->GetG(), work, true, false, 1.0, 0.0);
    // print(mass,"Mass");
}

void Operators::ComputeVol(){
    matmat(volumeTerm, I->GetD(),I->GetW(), true,false);
    // print(volumeTerm,"D^T*W");
}

void Operators::ComputeInvMassVol(){

    copy(massChol,mass);
    copy(invMassVol,volumeTerm);

    chol(massChol);
    cholsolve(invMassVol,massChol);
    // print(invMassVol,"invM*D^T*W");
}

int Operators::Nn()const{
    return mNumOfNodes;
}

int Operators::Nq() const{
    return mNumOfQuads;
}

DArray Operators::GetMass() const{
    DArray x(mNumOfNodes,mNumOfNodes);
    copy(x,mass);
    return x;
}

DArray Operators::GetCholMass() const{
    DArray x(massChol.size());
    copy(x,massChol);
    return x;
}

DArray Operators::GetVolTerm() const{
    DArray x(mNumOfNodes,mNumOfQuads);
    copy(x,volumeTerm);
    return x;
}

DArray Operators::GetInvMassVol() const{
    DArray x(mNumOfNodes,mNumOfQuads);
    copy(x,invMassVol);
    return x;
}

DArray Operators::GetG() const{
    return I->GetG();
}

DG::DG(Mesh1d* m, Parameters* p, Euler* eqn){
    mesh = m;
    params = p;
    equations = eqn;
    
    op = new Operators(params->nnodes, params->nquads);
    rs = new RiemannSolver(equations);
}

DG::~DG(){
    delete op;
    delete rs;
}