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

void DG::AssembleElement(DArray &u){
    DArray Fhat(2, params->neqs);
    DArray q(params->nnodes, params->neqs);
    DArray U(params->nquads, params->neqs, params->nels);
    DArray F(params->nquads, params->neqs, params->nels); 
    for (int iel = 0; iel < params->nels; ++iel){
        // Compute u at the quadrature points
        for (int ieq = 0; ieq < equations->NS(); ++ieq){
            // U = G*u
            DArray Uq(&U(0,ieq,iel),params->nquads,1);
            DArray un(&u(0,ieq,iel),params->nnodes,1);
            matvec(Uq,op->I->G,un);
        }

        // Compute Flux at quadrature points
        for (int i = 0; i < params->nquads; ++i){
            // F = F(U) --> U is computed in previous loop
            DArray s(params->neqs);
            DArray Fq(params->neqs);
            s[0] = U(i,0,iel);
            s[1] = U(i,1,iel);
            s[2] = U(i,2,iel);

            Fq = equations->Flux(s);
            F(i,0,iel) = Fq(0)*mesh->J()*mesh->invJ();
            F(i,1,iel) = Fq(1)*mesh->J()*mesh->invJ();
            F(i,2,iel) = Fq(2)*mesh->J()*mesh->invJ();
        }

        for (int ieq = 0; ieq < params->neqs; ++ieq){
            // Comptue Volumeterm
            DArray qslice(&q(0,ieq), params->nnodes,1);
            DArray Fslice(&F(0,ieq,iel), params->nquads,1);

            matvec(qslice, op->GetVolTerm(),Fslice,false,-1.0,0.0);
            
        }

    }
    
//     std::cout << mesh->J() << std::endl;
// std::cout << mesh->invJ() << std::endl;
    // print(u,"u at nodes");
    // print(U,"u at quads");
}

DG::~DG(){
    delete op;
    delete rs;
}