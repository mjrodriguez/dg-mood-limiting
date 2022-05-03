#pragma once

#include "Array.h"
#include "Interpolation.h"
#include "Problems.h"
#include "Mesh.h"
#include "Equations.h"

class Operators{
private:
    int mNumOfNodes;
    int mNumOfQuads;
    DArray volumeTerm;
    DArray invMassVol;

    void ComputeMass();
    void ComputeVol();
    void ComputeInvMassVol();
public:
    Operators(int numNodes, int numQuads);
    ~Operators();

    Interpolation* I;
    DArray mass, massChol;

    int Nn() const;
    int Nq() const;
    DArray GetMass() const;
    DArray GetCholMass() const;
    DArray GetVolTerm() const;
    DArray GetInvMassVol() const;
    DArray GetG() const;


};

class DG{
private:
    int dim = 1;
    int neqs;
    int order;
    int nels;
    int nquads;
public:
    Operators* op;
    Mesh1d* mesh;
    Parameters* params;
    Euler* equations;

    DG(Mesh1d *m, Parameters *p, Euler *eqn);
    ~DG();

    DArray AssembleElement(DArray &u);

    void RiemannSolver(DArray &stateL, DArray &stateR, DArray &Fhat);
    void LLF(DArray &stateL, DArray &stateR, DArray &Fhat); 

};
