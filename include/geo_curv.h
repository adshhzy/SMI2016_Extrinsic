#ifndef GEORF_H
#define GEORF_H

#include<vector>
#include"my_mesh.h"
#include"InfoStruct.h"
using namespace std;
using namespace Utility;
namespace n_rf {


class Curve{
public:

    int n_vertices, n_edges;
    vector<double>vertices;
    vector<uint>edges;
    vector<double>tangent;


private:

    vector< uint >vertices2edges;
    vector< uint >vertices2vertices;
    vector< unsigned long >vertices2edges_accumulate;
    vector<uchar>vertices2edgesPN;
    vector< uint >vertices2vertices_inverse;


private:
    bool isloop;



private:

    bool isSettangent;
    bool isBuildFrenet;
    double scale;
    double avedis;

    vector<double>linespeed;
    vector<double>curvature;
    vector<double>torsion;
    vector<double>Frenetnormal;
    vector<double>Frenetbinormal;
    vector<double>normCurvature;



public:
    vector<double>edge_len;
    vector<double>vertices_len;
    vector<double>inverse_edge_len;



public:
    void reset();
    void setparameters();

    int Initialize(string inputfile, string outputfile = string(), bool isEigenInit = true, bool isGaussIter = true);
    void BuildToyData(int ind,int numV,double start,double end);
    int SpecialCases(int sInd,string inputfile, string outputfile = string(),bool isEigenInit = true, bool isGaussIter = true);



    bool ReadCurve(string filename);
    bool ReadCurve(ifstream &fin, int p0, int p1, int p2);
    bool SaveCurve(string filename);
    bool ImportCurve(vector<double>& vertices_in, vector<uint>& edges_in);

    void RescaleUniform();
    bool EstimateFirstDerivativeO5(vector<double>& invect, vector<double>& derivative, int interval,bool isloop);

    void BuildEdges(bool isseq = true);
    bool SortLoopEdges();
    void BuildTangent();

    void EstimateFrenetFrame(bool isReSetTangent = true);
    void EliminateDuplication(double thres = 1e-5);
    bool isBuildDisplay();
    bool Isloop(){return isloop;}
    void ChangeLoop(bool isloop){this-> isloop = isloop;}

    void CurveScalarFunctionSmoothing(vector<double>&scalarF,int maxiter);


/********************************************************************************************/
/********************************************************************************************/
private:

    bool isbuild;


    vector<double>display_vertices;
    vector<double>display_vnormals;
    vector<uint>display_edges;
    vector<uint>display_faces;
    vector<unsigned char>v_color;



    static bool isload;
    static Mesh sphere;
    static Mesh cylinder;
    static Mesh cone;
    int n_ver_f1;
    int n_col_f1;
    int n_face_f1;
    int n_vnor_f1;

public:
    void BuildDisplay(bool isrescale = true, bool isuseball=true);
    void BuildDisplay(infoSet info);
    void BuildDisplay(bool isrescale, bool isvmf, bool isrmf, bool isnormal, bool isbinormal,bool istangent,
                      bool isrvector, bool isovector, bool isVector, bool isSurface,
                      double veclen, double thickness,double boxlen,double boxthickness,
                      bool isSurfColor, bool isColorFirst
                      );

private:
    bool load();
public:
    vector<double>* getDisplayVertices(){return &display_vertices;}
    vector<double>* getDisplayVerticesNormal(){return &display_vnormals;}
    vector<uint>* getDisplayEdges(){return &display_edges;}
    vector<uint>* getDisplayFaces(){return &display_faces;}
    vector<unsigned char>* getDisplayColor(){return &v_color;}

/********************************************************************************************/
/********************************************************************************************/


public:
    Curve();

    void vectorize(int first, int second,double *e){
        auto vdim = 3;
        auto p_v1 = &(vertices[first*vdim]);
        auto p_v2 = &(vertices[second*vdim]);
        minusVec(p_v1,p_v2,e,vdim);
    }
    void vectorize(double *v1, int second,double *e){
        auto vdim = 3;
        auto p_v2 = &(vertices[second*vdim]);
        minusVec(v1,p_v2,e,vdim);
    }
    double VerticesDistance(int first, int second)const{
        double dist = 0.0;
        auto p_v1 = &(vertices[first*3]);
        auto p_v2 = &(vertices[second*3]);
        for (int i = 0; i < 3; ++i)
            dist += pow((p_v1[i] - p_v2[i] ), 2);
        return sqrt(dist);
    }


    inline double* v_begin(int v_ind){ return &(vertices[v_ind * 3]); }
    inline double* v_end(int v_ind){ return &(vertices[(v_ind+1) * 3]); }
    inline uint* ev_begin(int e_ind){ return &(edges[e_ind * 2]); }
    inline uint* ev_end(int e_ind){ return &(edges[e_ind * 2 + 2]); }
    inline double* t_begin(int v_ind){ return &(tangent[v_ind * 3]); }
    inline double* t_end(int v_ind){ return &(tangent[(v_ind+1) * 3]); }



    inline double get_curvature(int v_ind){ return curvature[v_ind]; }
    inline double get_linespeed(int v_ind){ return linespeed[v_ind]; }
    inline double* nor_begin(int v_ind){ return &Frenetnormal[v_ind*3]; }
    inline double* nor_end(int v_ind){ return &Frenetnormal[v_ind*3+3]; }
    inline double* binor_begin(int v_ind){ return &Frenetbinormal[v_ind*3]; }
    inline double* binor_end(int v_ind){ return &Frenetbinormal[v_ind*3+3]; }


    inline double elen_begin(int e_ind){ return (edge_len[e_ind]); }
    inline double einvlen_begin(int e_ind){ return (inverse_edge_len[e_ind]); }

public:
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







public:

    vector<double> vertices_field;
    vector<double> vertices_rfield;

    vector<double> vertices_ofield;
    vector<double> vertices_rofield;

    vector<double> vertices_f_field;
private:
    double focusing[3];
    bool isCameraMode;
public:
    void InitializeField(bool isProj,bool isrand);
    void BuildOthogonalField();
    double RunGaussSeidelIteration();
    void InitializeFieldByLeastEigenV(bool isunified = true);
    void InverseField();

public:
    inline double* vvec_begin(int v_ind){ return &(vertices_field[v_ind * 3]); }
    inline double* vvec_end(int v_ind){ return &(vertices_field[(v_ind+1) * 3]); }

    inline double* vovec_begin(int v_ind){ return &(vertices_ofield[v_ind * 3]); }
    inline double* vovec_end(int v_ind){ return &(vertices_ofield[(v_ind+1) * 3]); }

    inline double* vrvec_begin(int v_ind){ return &(vertices_rfield[v_ind * 3]); }
    inline double* vrvec_end(int v_ind){ return &(vertices_rfield[(v_ind+1) * 3]); }

    inline double* vrovec_begin(int v_ind){ return &(vertices_rofield[v_ind * 3]); }
    inline double* vrovec_end(int v_ind){ return &(vertices_rofield[(v_ind+1) * 3]); }

    inline double* vfvec_begin(int v_ind){ return &(vertices_f_field[v_ind * 3]); }
    inline double* vfvec_end(int v_ind){ return &(vertices_f_field[(v_ind+1) * 3]); }
};

























}//n_rf

#endif // GEORF_H
