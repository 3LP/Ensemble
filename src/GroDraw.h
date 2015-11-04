//
//  GroDraw.h
//  ofMolecule
//
//  Created by Jackson Chief Elk on 4/30/15.
//
//
#ifndef __ofMolecule__GroDraw__
#define __ofMolecule__GroDraw__
#include <stdio.h>
#include <stdlib.h>
#include "ofMain.h"
#include <string>
#include <vector>
#define MAXBONDLENGTH 0.189 //Max Bond Length of 200picometers, Gro file is in Nanometers
#define MAXAA 1500

class GroDraw {
    
private:
    int i,j;
    
public:
    void setup();
    void update();
    void draw();
    void COM();

public:
    // Variables and Structures
  
    // Gromacs File Parameters
    int Natoms;
    float O1X[MAXAA],O1Y[MAXAA],O1Z[MAXAA]; /*Atom Cartesian Coordinates (INPUT)*/
    char Atom_type[2], ResID[8];
    int AtomID;
    // Structure File
    ofFile myfile;
    //Atomic Vector Arrays and Variables Needed For Transformation
    float VEC1[MAXAA],VEC2[MAXAA],VEC3[MAXAA],Distance[MAXAA][MAXAA]; // COM Vectors
    float ofVecX[MAXAA],ofVecY[MAXAA],ofVecZ[MAXAA]; //OpenFrameWorks Vectors
    int adjgraph[MAXAA][MAXAA];
    float x1,y1,z1; /*Molecular Centroid*/
    float length_increase = 200.0; // Controlls molecular size in OpenFrameworks
    float cumsum; /*Cummulative Sum*/
    // OpenFrameWorks Initial Conditions
    ofVec3f ofmol_com; // Center of Mass of Molecule in OpenF World


    
};

#endif /* defined(__ofMolecule__GroDraw__) */
