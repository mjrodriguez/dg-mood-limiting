#pragma once

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "Array.h"
#include "ArrayUtils.h"

/* ===============================================================
Nodes Class
=============================================================== */

class Nodes{
private:
    int mNumOfNodes;
    DArray mNodes;
    DArray mWeights;

public:
    Nodes(int numOfNodes);
    ~Nodes();
    int GetSize() const;
    DArray GetNodes() const;
    DArray GetWeights() const;

    void lglnodes();
};




/* ===============================================================
Interpolation Class
=============================================================== */

class Interpolation{
private:
    int mNumOfNodes;
    int mNumOfQuads;
    Nodes* xnodes;
    Nodes* xquads;
    DArray G;
    DArray D;
    DArray W;
public:
    Interpolation(int numOfNodes, int numOfQuads);
    ~Interpolation();

    int GetSizeQuads() const;
    int GetSizeNodes() const;
    DArray GetNodes() const;
    DArray GetQuads() const;
    DArray GetG() const;
    DArray GetD() const;
    DArray GetW() const;

    void lagint();
};

