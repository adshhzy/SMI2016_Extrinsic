#ifndef GEO_SUR_H
#define GEO_SUR_H


#include"InfoStruct.h"
#include<vector>
#include"my_mesh.h"
#include<eigen3/Eigen/Dense>
using namespace std;
using namespace Utility;

namespace n_rf {


class Surface: public Mesh{
public:
    bool reinitflag;

private:
    string modelname;
    string prepath;
    string ext;
private:
    vector<double> faces_center;
    vector<double> edge_len;
    vector<double> edge_vec;
    vector<double> edge_cot;
    double ave_edge_len;

    vector<int> face_Material1;
    vector<int> face_Material2;
private:

    vector<double>vertices_area;
    vector<double>faces_area;
    vector<double>faces_inverse_area;



private:

    vector<double>id_pick_vertices;
    vector<uint>id_pick_faces;
    vector<unsigned char>id_pick_color;


public:
    void clearup();

private:
    bool saveObjx(string filename);
    bool readSufFile(string filename);

    //bool readEobjfile(string filename);
    //bool readObjxfile(string filename);

    void BuildFacesCenter();
    void ComputeEdgeLength();
    void ComputeArea();


    void meshScalarfunctionSmoothing(vector<double>& scalarF,int maxiter);
    void testSlicer(int thres);
private:
    void BuildPickID();
public:
    vector<double>* getPickIDVertices(){return &id_pick_vertices;}
    vector<uint>* getPickIDFaces(){return &id_pick_faces;}
    vector<unsigned char>* getPickIDColor(){return &id_pick_color;}
public:
    int PickFaceViaColor(unsigned char* pcolor);



public:


    int Initialize(string inputfile, string outputfile = string(), bool isEigenInit = true, bool isGaussIter = true);
    bool ReadFile(string filename);
    bool SaveInterface(string filename);

public:
    Surface();


public:

    inline double* fcent_begin(int f_ind){ return &(faces_center[f_ind * 3]); }
    inline double* fcent_end(int f_ind){ return &(faces_center[(f_ind+1) * 3]); }

    inline double elen_begin(int e_ind){ return (edge_len[e_ind]); }
    inline double *evec_begin(int e_ind){return &(edge_vec[e_ind*3]);}

    inline double ecot_begin(int e_ind){ return (edge_cot[e_ind]); }
    inline double* ecotp_begin(int e_ind){ return &(edge_cot[e_ind]); }


    inline double v_area(int v_ind){return vertices_area[v_ind];}
    inline double f_area(int v_ind){return faces_area[v_ind];}


public:

    void ShowOptInfo(int isacc);

public:

    double computeRotationAlongAxis(const double angle, const double *norAxis, const double *vec, double *vecout);

/********************************************************************************************/
/********************************************************************************************/


private:

    bool isbuilddisp;
    bool isresetColor;
    bool isFaceRenderMode;
    bool ismarkC;
    static bool isload;
    static Mesh sphere;
private:
    double colordegree;
    int sparsecoff;
    vector<uchar> weighted_color;

    vector<uchar> weighted_fcolor;
    vector<uchar> sparseShowField;
    vector<double>invert_faces_colordegree;
    vector<double>invert_vertices_colordegree;
private:
    vector<double> display_vertices;
    vector<double> display_normal;
    vector<uint> display_edges;

    vector<uint>display_field_dot;
    vector<unsigned char> display_vcolor;
    vector<uint> display_faces;
public:
    void load();

public:
    vector<double>* getDisplayVertices(){return &display_vertices;}
    vector<double>* getDisplayVerticesNormal(){return &display_normal;}
    vector<uint>* getDisplayEdges(){return &display_edges;}
    vector<uint>* getDisplayFaces(){return &display_faces;}
    vector<unsigned char>* getDisplayColor(){return &display_vcolor;}


    bool isBuildDisp(){return isbuilddisp;}

    void WeightColor(const double thres,const vector<double>&weights,vector<uchar>&out_color);
    void sparseSampling(int a);
    void BuildDisplay(double colordegree,bool isfield, bool isnormal, bool issurf, bool iswire, bool issingularity, bool ismark, int length, int width, int upnormal);

    void BuildDisplay(infoSurfDisp info);
    void GetPerFacesColorDegree(vector<double>&facedegree);

/********************************************************************************************/
/********************************************************************************************/



private:

    vector<double> vertices_field;
    vector<double> vertices_rfield;

    vector<double> vertices_ofield;
    vector<double> vertices_rofield;

private:
    vector<bool> singularityf;
    vector<uint> singbufferf;

    vector<bool> singularityv;
    vector<uint> singbufferv;
public:
    void InitializeField(bool isProj,bool isrand);
    void BuildOthogonalField();
    int ComputeVerticesFieldSingularity();
    double RunGaussSeidelIteration();
    void InitializeFieldByLeastEigenV(bool isunified = true);
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









}// n_rf




#endif // GEO_SUR_H
