//
//  GroDraw.cpp
//  ofMolecule
//
//  Created by Jackson Chief Elk on 4/30/15.
//
//

#include "GroDraw.h"

void GroDraw::setup() {
 // Open FrameWorks Parameters
 //
    ofmol_com.x = 0.0;
    ofmol_com.y = 220.0;
    ofmol_com.z = -100.0;
  //
 // Read in Atomic Coordinates
    Natoms = 53;
    string line;
    for(int j = 0;j< Natoms;j++) {
        O1X[j] = O1Y[j] = O1Z[j] = ofVecX[i] = ofVecY[i] = ofVecZ[i] = 0.0;
    }
  //
  //
    
    myfile.open(ofToDataPath("peptide.gro"), ofFile::ReadWrite, false);
    int j = 0;
    getline(myfile,line);
    // getline(myfile,line);
    while ( j < Natoms )
        {
            getline(myfile,line);
            myfile >> ResID >> Atom_type >> AtomID >> O1X[j] >> O1Y[j] >> O1Z[j];
            // cout << O1X[j] << "\n";
            j++;
        }

        myfile.close();
    //
    
    //Calculate Center of Mass
    //
    x1=y1=z1=0.0;
    for(i=0;i<Natoms;i++) {
        x1+=O1X[i];
        y1+=O1Y[i];
        z1+=O1Z[i];
    }
    //
    // Centers of Mass and distance
    //
    x1/=Natoms;
    y1/=Natoms;
    z1/=Natoms;
    for(i=0;i<Natoms;i++){
        VEC1[i] = O1X[i]-x1;
        VEC2[i] = O1Y[i]-y1;
        VEC3[i] = O1Z[i]-z1;
        ofVecX[i] = length_increase*VEC1[i]+ofmol_com.x;
        ofVecY[i] = length_increase*VEC2[i]+ofmol_com.y;
        ofVecZ[i] = length_increase*VEC3[i]+ofmol_com.z;
        // cout << ofVecX[i] << "\t" << ofVecY[i] << "\t" <<ofVecZ[i] << "\n";
    }
     
    // First Generate Adjacency Graph
    //Initialize
    for(i=0;i<MAXAA;i++){
        for(j=0;j<MAXAA;j++){
            adjgraph[i][j]=0;
        }
    }
    
    //
    for(i=0;i<Natoms;i++){
        
        for(j=0;j<Natoms;j++){
            Distance[i][j] = 0.0;
            // First Calculate Pair-wise Distances
            if (i!=j) Distance[i][j] = sqrt((O1X[i]-O1X[j])*(O1X[i]-O1X[j])+(O1Y[i]-O1Y[j])*(O1Y[i]-O1Y[j])+(O1Z[i]-O1Z[j])*(O1Z[i]-O1Z[j]));
            cout << "Distances:\t" << Distance[i][j] << "\n";
        }
        

        
    }
    
    //
    // Setup Everything For First Update
    //
    // Grompp Steps
    std::string struc_file = ofToDataPath("full.gro");
    std::string input_struc = ofToDataPath("peptide.gro");
    std::string input_top = ofToDataPath("peptide.top");
    std::string input_mdp = ofToDataPath("md-imp.mdp");
    std::string tpr = ofToDataPath("full.tpr");
    std::string edr = ofToDataPath("ener.edr");
    std::string mdout = ofToDataPath("mdout.mdp");
    std::string mdlog = ofToDataPath("md.log");
    std::string cpt = ofToDataPath("state.cpt");
    std::string grompp_path = "/usr/local/bin/grompp";
    char c1[500];
    std::string syscommand;
    syscommand = "export GMX_MAXBACKUP=-1";
    strcpy(c1, syscommand.c_str());
    system(c1);
    syscommand = grompp_path + " -f " + input_mdp + " -c " + input_struc + " -p " +input_top + " -o " + tpr + " -po " + mdout + " -maxwarn 1";
    strcpy(c1, syscommand.c_str());
    system(c1);
    //MDrun steps
    std::string trr = ofToDataPath("full.trr");
    std::string mdrun_path = "/usr/local/bin/mdrun";
    syscommand = mdrun_path + " -s " + tpr + " " +  "-o " + trr + " -nt 1 " + " -g " + mdlog + " -e " + edr + " -c " + struc_file + " -pd";
    strcpy(c1, syscommand.c_str());
    system(c1);
    
}
    

void GroDraw::update() {
    
    // Log Stuff Required - Memory Management
    std::string md_struc = ofToDataPath("full.gro");
    std::string min_struc = ofToDataPath("min.gro");
    std::string mdout = ofToDataPath("mdout.mdp");
    std::string mdlog = ofToDataPath("md.log");
    std::string edr = ofToDataPath("ener.edr");
    std::string cpt = ofToDataPath("state.cpt");
    // Grompp Steps
    std::string input_struc = ofToDataPath("full.gro");
    std::string input_top = ofToDataPath("peptide.top");
    std::string md_mdp = ofToDataPath("md-imp.mdp");
    std::string min_mdp = ofToDataPath("min-imp.mdp");
    std::string output_tpr = ofToDataPath("full.tpr");
    std::string min_tpr = ofToDataPath("min.tpr");
    std::string grompp_path = "/usr/local/bin/grompp";
    std::string mdrun_path = "/usr/local/bin/mdrun";
    std::string mintrr = ofToDataPath("min.trr");
    char c1[500];
    std::string syscommand;
    syscommand = "export GMX_MAXBACKUP=-1";
    strcpy(c1, syscommand.c_str());
    system(c1);
    // Energy Minimization
    syscommand = grompp_path + " -f " + min_mdp + " -c " + md_struc + " -p " + input_top + " -o " + min_tpr + " -po " + mdout + " -maxwarn 1";
    strcpy(c1, syscommand.c_str());
    system(c1);
    syscommand = mdrun_path + " -s " + min_tpr + " -o " + mintrr + " -g " + mdlog + " -nt 1 " + " -e " + edr + " -c " + min_struc  + " -pd";
    strcpy(c1, syscommand.c_str());
    system(c1);
    // Memory Cleanup
    syscommand = "rm " +  mintrr;
    strcpy(c1, syscommand.c_str());
    system(c1);
    syscommand = "rm " +  mdout;
    strcpy(c1, syscommand.c_str());
    system(c1);
    syscommand = "rm " +  mdlog;
    strcpy(c1, syscommand.c_str());
    system(c1);
    syscommand = "rm " + edr;
    strcpy(c1, syscommand.c_str());
    system(c1);
    syscommand = "rm " + min_tpr;
    strcpy(c1, syscommand.c_str());
    system(c1);
    // Test
    std::string pound1 = ofToDataPath("#full.gro.1#");
    syscommand = "rm " + pound1;
    strcpy(c1, syscommand.c_str());
    system(c1);
    std::string pound2 = ofToDataPath("#full.gro.2#");
    syscommand = "rm " + pound2;
    strcpy(c1, syscommand.c_str());
    system(c1);

    
    
    
    
    std::string pound3 = ofToDataPath("#full.gro.3#");
    syscommand = "rm " + pound1;
    strcpy(c1, syscommand.c_str());
    system(c1);
    std::string pound4 = ofToDataPath("#full.gro.4#");
    syscommand = "rm " + pound2;
    strcpy(c1, syscommand.c_str());
    system(c1);
    syscommand = "rm " + md_struc;
    strcpy(c1, syscommand.c_str());
    system(c1);

    

    
    // MD
    syscommand = grompp_path + " -f " + md_mdp + " -c " + min_struc + " -p " +input_top + " -o " + output_tpr + " -po " + mdout + " -maxwarn 1";
    strcpy(c1, syscommand.c_str());
    system(c1);
    std::string tpr = ofToDataPath("full.tpr");
    std::string trr = ofToDataPath("full.trr");
    syscommand = mdrun_path + " -s " + tpr +  " -o " + trr + " -g " + mdlog + " -nt 1 " + " -e " + edr + " -c " + md_struc  + " -pd";
    strcpy(c1, syscommand.c_str());
    system(c1);
    // Coordinate Grab
    myfile.open(ofToDataPath("full.gro"), ofFile::ReadWrite, false);
    int j = 0;
    string line;
    getline(myfile,line);
    // getline(myfile,line);
    while ( j < Natoms )
    {
        getline(myfile,line);
        myfile >> ResID >> Atom_type >> AtomID >> O1X[j] >> O1Y[j] >> O1Z[j];
        // cout << O1X[j] << "\n";
        j++;
    }
    
    myfile.close();
    //
    
    //Calculate Center of Mass
    //
    x1=y1=z1=0.0;
    for(i=0;i<Natoms;i++) {
        x1+=O1X[i];
        y1+=O1Y[i];
        z1+=O1Z[i];
    }
    //
    // Centers of Mass and distance
    //
    x1/=Natoms;
    y1/=Natoms;
    z1/=Natoms;
    for(i=0;i<Natoms;i++){
        VEC1[i] = O1X[i]-x1;
        VEC2[i] = O1Y[i]-y1;
        VEC3[i] = O1Z[i]-z1;
        ofVecX[i] = length_increase*VEC1[i]+ofmol_com.x;
        ofVecY[i] = length_increase*VEC2[i]+ofmol_com.y;
        ofVecZ[i] = length_increase*VEC3[i]+ofmol_com.z;
        // cout << ofVecX[i] << "\t" << ofVecY[i] << "\t" <<ofVecZ[i] << "\n";
    }
    
    // Memory Cleanup
    syscommand = "rm " +  trr;
    strcpy(c1, syscommand.c_str());
    system(c1);
    syscommand = "rm " +  mdout;
    strcpy(c1, syscommand.c_str());
    system(c1);
    syscommand = "rm " +  mdlog;
    strcpy(c1, syscommand.c_str());
    system(c1);
    syscommand = "rm " + edr;
    strcpy(c1, syscommand.c_str());
    system(c1);
    syscommand = "rm " + tpr;
    strcpy(c1, syscommand.c_str());
    system(c1);
    syscommand = "rm " + cpt;
    strcpy(c1, syscommand.c_str());
    system(c1);
    


    
}

void GroDraw::draw() {

    // Draw Sphere in Physics World
    ofSetColor(255,255,0);
    for(i=0;i<Natoms;i++){
        ofDrawSphere(ofVecX[i],ofVecY[i],ofVecZ[i],5);
 
    }
    //
    // Connect atoms with lines
    //
   
    ofSetLineWidth(3);
    for(i = 0;i<Natoms; i++){
        
        for(j=0;j<Natoms;j++){
            
            if (i!=j && Distance[i][j] <= 0.16){
                // If there is bond draw a line between atoms
                ofLine(ofVec3f(ofVecX[i],ofVecY[i],ofVecZ[i]),ofVec3f(ofVecX[j],ofVecY[j],ofVecZ[j]));
                
            }
            
        }
        
    }
    

    

    
}

