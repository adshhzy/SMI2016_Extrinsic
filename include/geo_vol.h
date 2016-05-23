#ifndef GEO_VOL_H
#define GEO_VOL_H




#include"InfoStruct.h"
#include<vector>
#include"my_mesh.h"
#include "geo_sur.h"
#include<eigen3/Eigen/Dense>
using namespace std;
using namespace Utility;



namespace n_rf {


class Volume{



public:

    int n_vertices, n_edges,n_tetr;
    vector<double>vertices;
    vector<uint>edges;
    vector<uint>tetrahedron2vertices;


private:
    vector<double> vertices_normal;



private:

    vector< uint >vertices2edges;
    vector< uint >vertices2vertices;
    vector< unsigned long >vertices2edges_accumulate;
    vector<uchar>vertices2edgesPN;
    vector< uint >vertices2vertices_inverse;

private:

    vector< uint >vertices2tetr;
    vector< unsigned long >vertices2tetr_accumulate;
    vector< uint >vertices2tetr_inverse;



public:
    Volume();
    void clearup();
    void load();
    void setparameters();


public:

    void ToyGeneration();
    //    void DistanceFieldGeneration();
    //    void ToyGeneration2();
    bool ReadInterface(string filename);
    bool SaveInterface(string filename);


public:
    void BuildEdgesRelations();

    void Rescale_uniform(double lscale);


public:
    double EdgesDifferetialEnergy();

public:

    inline double* v_begin(uint v_ind){ return &(vertices[v_ind * 3]); }
    inline double* v_end(uint v_ind){ return &(vertices[(v_ind+1) * 3]); }
    inline uint* vt_begin(uint v_ind){ return &(vertices2tetr[vertices2tetr_accumulate[v_ind]]); }
    inline uint* vt_end(uint v_ind){ return &(vertices2tetr[vertices2tetr_accumulate[v_ind+1]]); }
    inline uint vt_num(uint v_ind){ return vertices2tetr_accumulate[v_ind+1]-vertices2tetr_accumulate[v_ind]; }
    inline uint* vtinv_begin(uint v_ind){ return &(vertices2tetr_inverse[vertices2tetr_accumulate[v_ind]]); }
    inline uint* vtinv_end(uint v_ind){ return &(vertices2tetr_inverse[vertices2tetr_accumulate[v_ind+1]]); }


    inline uint* tv_begin(uint t_ind){ return &(tetrahedron2vertices[t_ind * 4]); }
    inline uint* tv_end(uint t_ind){ return &(tetrahedron2vertices[t_ind * 4 + 4]); }

    inline uint* ev_begin(uint e_ind){ return &(edges[e_ind * 2]); }
    inline uint* ev_end(uint e_ind){ return &(edges[e_ind * 2 + 2]); }
    inline uint* ve_begin(uint v_ind){ return &(vertices2edges[vertices2edges_accumulate[v_ind]]); }
    inline uint* ve_end(uint v_ind){ return &(vertices2edges[vertices2edges_accumulate[v_ind+1]]); }
    inline uint ve_num(uint v_ind){ return vertices2edges_accumulate[v_ind+1]-vertices2edges_accumulate[v_ind]; }
    inline uchar* vepn_begin(uint v_ind){ return &(vertices2edgesPN[vertices2edges_accumulate[v_ind]]); }
    inline uchar* vepn_end(uint v_ind){ return &(vertices2edgesPN[vertices2edges_accumulate[v_ind+1]]); }


    inline uint* vv_begin(uint v_ind){ return &(vertices2vertices[vertices2edges_accumulate[v_ind]]); }
    inline uint* vv_end(uint v_ind){ return &(vertices2vertices[vertices2edges_accumulate[v_ind+1]]); }
    inline uint* vvinv_begin(uint v_ind){ return &(vertices2vertices_inverse[vertices2edges_accumulate[v_ind]]); }
    inline uint* vvinv_end(uint v_ind){ return &(vertices2vertices_inverse[vertices2edges_accumulate[v_ind+1]]); }
    inline uint vv_num(uint v_ind){ return vertices2edges_accumulate[v_ind+1]-vertices2edges_accumulate[v_ind]; }

    inline double* vnor_begin(uint v_ind){ return &(vertices_normal[v_ind * 3]); }
    inline double* vnor_end(uint v_ind){ return &(vertices_normal[v_ind * 3 + 3]); }


    /********************************************************************************************/
    /********************************************************************************************/

private:
    static bool isload;
    static Mesh sphere;
    static Mesh cylinder;
    static Mesh cone;

private:
    bool isbuilddisp;

private:
    vector<double> display_vertices;
    vector<double> display_normal;
    vector<uint> display_edges;
    vector<unsigned char> display_vcolor;
    vector<uint> display_faces;

    vector<bool> sparseShowField;
    int sparsecoff;


public:
    void sparseSampling(int a);
    void BuildDisplay(infoVolDisp info);
    vector<double>* getDisplayVertices(){return &display_vertices;}
    vector<double>* getDisplayVerticesNormal(){return &display_normal;}
    vector<uint>* getDisplayEdges(){return &display_edges;}
    vector<unsigned char>* getDisplayColor(){return &display_vcolor;}
    vector<uint>* getDisplayFaces(){return &display_faces;}

    bool isBuildDisp(){return isbuilddisp;}

    /********************************************************************************************/
    /********************************************************************************************/






private:
    int width,height,depth;
    double xsa_step,ysa_step,zsa_step;
    double xs,ys,zs;
    double thres_upper,thres_lower;
    vector<vector<double> >volume;
    double (*vfun)(double x,double y,double z,double *grad);
private:
    Surface isosurface;
public:
    void IsoToyGeneration(int w, int h, int d, double wrange, double hrange, double drange);


private:
    vector<double> vertices_field;
    vector<double> vertices_rfield;

private:

    vector<double> vertices_ofield;
    vector<double> vertices_rofield;
public:
    int Initialize(string inputfile, string outputfile = string(), bool isEigenInit = true, bool isGaussIter = true);
    int SpecialCases(string outputfile = string(), bool isEigenInit = true, bool isGaussIter = true);
public:
    double RunGaussSeidelIteration();
    void InitializeFieldByLeastEigenV(bool isunified = true);
    void InitializeField(bool isProj = false, bool isrand = true);
    void RandomInitialize(int Indmethod);
    void BuildOthogonalField(bool isR = false);
public:

    inline double* vvec_begin(int v_ind){ return &(vertices_field[v_ind * 3]); }
    inline double* vvec_end(int v_ind){ return &(vertices_field[(v_ind+1) * 3]); }


    inline double* vovec_begin(int v_ind){ return &(vertices_ofield[v_ind * 3]); }
    inline double* vovec_end(int v_ind){ return &(vertices_ofield[(v_ind+1) * 3]); }


    inline double* vrvec_begin(int v_ind){ return &(vertices_rfield[v_ind * 3]); }
    inline double* vrvec_end(int v_ind){ return &(vertices_rfield[(v_ind+1) * 3]); }

    inline double* vrovec_begin(int v_ind){ return &(vertices_rofield[v_ind * 3]); }
    inline double* vrovec_end(int v_ind){ return &(vertices_rofield[(v_ind+1) * 3]); }



};






}


















#endif // GEO_VOL_H
