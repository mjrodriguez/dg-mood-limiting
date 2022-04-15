#ifndef DGOPERATORS_H
#define DGOPERATORS_H

#include "Array.h"
#include "Interpolation.h"

class Operators{
private:
    int mNumOfNodes;
    int mNumOfQuads;
    Interpolation* I;
    DArray mass, massChol;
    DArray volumeTerm;
    DArray invMassVol;

    void ComputeMass();
    void ComputeVol();
    void ComputeInvMassVol();
public:
    Operators(int numNodes, int numQuads);
    ~Operators();

    int Nn() const;
    int Nq() const;
    DArray GetMass() const;
    DArray GetCholMass() const;
    DArray GetVolTerm() const;
    DArray GetInvMassVol() const;
    DArray GetG() const;


};


#endif
