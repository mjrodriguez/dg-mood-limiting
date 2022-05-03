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
}

DArray DG::AssembleElement(DArray &u){
    DArray Fhat(2, params->neqs);
    DArray q(params->nnodes, params->neqs, params->nels);
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

        // print(DArray (&F(0,0,iel),params->nquads, params->neqs),"F");

        for (int ieq = 0; ieq < params->neqs; ++ieq){
            // Comptue Volumeterm
            DArray qslice(&q(0,ieq, iel), params->nnodes,1);
            DArray Fslice(&F(0,ieq, iel), params->nquads,1);

            matvec(qslice, op->GetVolTerm(),Fslice,false,-1.0,0.0);

            // Compute Numerical Fluxes
            if (iel == 0){
                DArray ul(params->leftBC.data(),params->neqs); // ul will in conservative vars
                DArray ur(params->neqs);
                DArray Fs(params->neqs);

                // Computing interface flux at left face of first element 
                ur[0] = u(0,0,iel);
                ur[1] = u(0,1,iel);
                ur[2] = u(0,2,iel);

                RiemannSolver(ul,ur,Fs);
                Fhat(0,ieq) = Fs[ieq];

                // Computinig interface flux at right face of first element
                ul[0] = u(params->nnodes-1,0,iel);
                ul[1] = u(params->nnodes-1,1,iel);
                ul[2] = u(params->nnodes-1,2,iel);

                ur[0] = u(0,0,iel+1);
                ur[1] = u(0,1,iel+1);
                ur[2] = u(0,2,iel+1);

                RiemannSolver(ul,ur,Fs);
                Fhat(1,ieq) = Fs[ieq];


            }
            else if (iel == params->nels-1){
                DArray ul(params->neqs);
                DArray ur(params->neqs);
                DArray Fs(params->neqs);

                // Computing interface flux at left face of last element
                ul[0] = u(params->nnodes-1,0,iel-1);
                ul[1] = u(params->nnodes-1,1,iel-1);
                ul[2] = u(params->nnodes-1,2,iel-1);

                ur[0] = u(0,0,iel);
                ur[1] = u(0,1,iel);
                ur[2] = u(0,2,iel);
                
                RiemannSolver(ul,ur,Fs);
                Fhat(0,ieq) = Fs[ieq];

                // Computinig interface flux at right face
                ul[0] = u(params->nnodes-1,0,iel);
                ul[1] = u(params->nnodes-1,1,iel);
                ul[2] = u(params->nnodes-1,2,iel);

                ur[0] = params->rightBC[0];
                ur[1] = params->rightBC[1];
                ur[2] = params->rightBC[2];

                RiemannSolver(ul,ur,Fs);
                Fhat(1,ieq) = Fs[ieq];

                print(Fhat, "fhat");
            }
            else {
                DArray ul(params->neqs);
                DArray ur(params->neqs);
                DArray Fs(params->neqs);

                // Computing interface flux at left face
                ul[0] = u(params->nnodes-1,0,iel-1);
                ul[1] = u(params->nnodes-1,1,iel-1);
                ul[2] = u(params->nnodes-1,2,iel-1);

                ur[0] = u(0,0,iel);
                ur[1] = u(0,1,iel);
                ur[2] = u(0,2,iel);
                
                RiemannSolver(ul,ur,Fs);
                Fhat(0,ieq) = Fs[ieq];

                // Computinig interface flux at right face
                ul[0] = u(params->nnodes-1,0,iel);
                ul[1] = u(params->nnodes-1,1,iel);
                ul[2] = u(params->nnodes-1,2,iel);

                ur[0] = u(0,0,iel+1);
                ur[1] = u(0,1,iel+1);
                ur[2] = u(0,2,iel+1);

                RiemannSolver(ul,ur,Fs);
                Fhat(1,ieq) = Fs[ieq];

            }

            qslice[0] -= Fhat(0,ieq);
            qslice[params->nnodes-1] += Fhat(1,ieq);

            cholsolve(qslice,op->GetCholMass());

            // print(qslice);
        }

    }
    
    scale(q,-1.0/mesh->detJ());
    // print(q, "q: "); 
//     std::cout << mesh->J() << std::endl;
// std::cout << mesh->invJ() << std::endl;
    // print(u,"u at nodes");
    // print(U,"u at quads");

    return q;
}

DG::~DG(){
    delete op;
}


void DG::RiemannSolver(DArray &stateL, DArray &stateR, DArray &Fhat){

    if (params->riemannSolver == "llf"){
        LLF(stateL, stateR, Fhat);
    }
}

void DG::LLF(DArray &stateL, DArray & stateR, DArray &Fhat){
	assert(equations->IsStatePhysical(stateL));
	assert(equations->IsStatePhysical(stateR));

	double maxL = equations->MaxVel(stateL);
	double maxR = equations->MaxVel(stateR);
	double Cmax = std::max(maxL, maxR); 

	DArray FL(params->neqs), FR(params->neqs);
    FL = equations->Flux(stateL);
    scale(FL,mesh->J()*mesh->invJ());

	FR = equations->Flux(stateR);
    scale(FR,mesh->J()*mesh->invJ());

	for (int ieq = 0; ieq < equations->NS(); ++ieq){
		Fhat[ieq] = 0.5*(FL[ieq] + FR[ieq]) - 0.5*Cmax*(stateR[ieq] - stateL[ieq]);
	}


}