#ifndef A_ESDF_CORE_H
#define A_ESDF_CORE_H

#include "geo_sur.h"
#include "geo_curv.h"
#include "geo_vol.h"
namespace n_rf {


class ESDF_Core{

public:
    enum{CURV,SURF,VOL};
    int curstate;
private:
    int sparsecoff;
private:
    string modelname;
    string prepath;
    string ext;
private:
    Surface Surf_core;
    Curve Curv_core;
    Volume Vol_core;


private:
    vector< vector<double>* > display_vertices;
    vector< vector<double>* >display_normal;
    vector< vector<uint>* >display_edges;


    vector< vector<uint>* >display_field_dot;
    vector< vector<unsigned char>* >display_vcolor;
    vector< vector<uint>* >display_faces;


public:

    ESDF_Core(){curstate = SURF;}

    bool clearup();
    //bool ReadFile(string filename);

    int ReadArgs(int argc,char** argv);

    void BuildDisplay(infoSurfDisp info);


//public:
//    vector<vector<double>*>* getDisplayVertices(){return &display_vertices;}
//    vector<vector<double>*>* getDisplayVerticesNormal(){return &display_normal;}
//    vector<vector<uint>*>* getDisplayEdges(){return &display_edges;}
//    vector<vector<uint>*>* getDisplayFaces(){return &display_faces;}
//    vector<vector<uchar>*>* getDisplayColor(){return &display_vcolor;}

public:
    vector<double>* getDisplayVertices();
    vector<double>* getDisplayVerticesNormal();
    vector<uint>* getDisplayEdges();
    vector<uint>* getDisplayFaces();
    vector<unsigned char>* getDisplayColor();

};


}//n_rf













#endif // A_ESDF_CORE_H
