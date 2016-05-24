#include "geo_curv.h"
#include "readers.h"
#include <eigen3/Eigen/Geometry>



namespace n_rf {


int Curve::Initialize(string inputfile,string outputfile, bool isEigenInit, bool isGaussIter){


    //ReadCurve(string("/Users/Research/Geometry/matlabspace/drawcurveT.cur"));
    //ReadCurve(inputfile);
    if(!readCurfFile(inputfile,vertices,edges,vertices_field,vertices_f_field)){
        cout<<"invalid input: "<< inputfile <<endl;
        return -1;
    }
    setparameters();
    SortLoopEdges();
    //RescaleUniform();
    if(tangent.size()==0)EstimateFrenetFrame();
    BuildEdges(false);
    InitializeField(true,false);

    if(isEigenInit)InitializeFieldByLeastEigenV();
    if(isGaussIter)RunGaussSeidelIteration();


    BuildDisplay();


    BuildDisplay(false,true, false,
                 false,false, false,
                 true,false,
                 true,true,
                 1,1,
                 1.4,1.2,
                 false,false);



    //cout<<"outputfile: "<<outputfile<<endl;
    if(outputfile.size()!=0){
        SaveCurve(outputfile);
        outputDisplay(outputfile);
    }

    return 0;
}

int Curve::SpecialCases(int sInd,string inputfile, string outputfile,bool isEigenInit, bool isGaussIter){
    if(sInd==0){
        isCameraMode =true;
        focusing[0] = 0.94;
        focusing[1] = -1.;
        focusing[2] = 0.0;
        BuildToyData(1,100,-1.5,8.3);
        InitializeField(true,false);
        //if(isEigenInit)InitializeFieldByLeastEigenV();
        if(isGaussIter)RunGaussSeidelIteration();

        BuildDisplay(false,true, false,
                     false,false, false,
                     true,false,
                     true,false,
                     0.6,0.6,
                     1,1,
                     false,false);


    }else if(sInd>0 && sInd<3){

        if(sInd==1)BuildToyData(sInd,100,0,6.30);
        if(sInd==2)BuildToyData(sInd,150,0,6.30);
        InitializeField(true,true);
        if(isEigenInit)InitializeFieldByLeastEigenV();
        if(isGaussIter)RunGaussSeidelIteration();

        if(sInd==2)InverseField();
        BuildDisplay(false,true, false,
                     false,false, false,
                     true,false,
                     true,true,
                     1,1,
                     1.4,1.2,
                     false,false);

    }







    if(outputfile.size()!=0){
        SaveCurve(outputfile);
        outputDisplay(outputfile);
    }

    return 0;
}

bool Curve::SaveCurve(string filename){
    if(n_vertices==0 || n_edges ==0) return false;


    //cout<<"writeCurfFile(filename,vertices,edges,vertices_field,tangent)"<<endl;
    return writeCurfFile(filename,vertices,edges,vertices_field,vertices_f_field);
}

bool Curve::outputDisplay(string filename){

    return writePLYFile(filename,display_vertices,display_faces,display_vnormals,v_color);

}
void Curve::BuildDisplay(infoSet info){

    BuildDisplay(info.isrescale,info.isvmf, info.isrmf,
                 info.isnnormal,info.isbinormal, info.istangent,
                 info.isrvector,info.isovector,
                 info.isVector,info.isSurface,
                 info.veclen,info.vecthickness,
                 info.boxlen,info.boxthickness,
                 info.isSurfColor,info.isColorFirst);


}

double GaussSeidelIteration(const vector<double>&vertices_normal,
                            const vector<double>&edge_weight,
                            const vector<uint>&edges2vertices,const vector<uint>&vertices2edges,
                            const vector<uint>&vertices2vertices, const vector< unsigned long >&vertices2edges_accumulate,
                            vector<double>&vertices_field,
                            const int maxiter = 5000
        );

double Curve::RunGaussSeidelIteration(){


    auto edge_weight = edge_len;
    for(auto &a:edge_weight)a=1/a;


    double en_pre = GaussSeidelIteration(vertices_f_field,edge_weight,edges,vertices2edges,vertices2vertices,vertices2edges_accumulate,vertices_field,50000);


    BuildOthogonalField();
    return en_pre;

}



void EigenInitialization(    const vector<double>&vertices_rfield,const vector<double>&vertices_rofield,
                             const vector<double>&edge_weight,const vector<double>&vertices_Vweight,const vector<double>&edge_Vweight,
                             const vector<uint>&edges2vertices,const vector<uint>&vertices2edges,
                             const vector<uint>&vertices2vertices, const vector< unsigned long >&vertices2edges_accumulate,
                             vector<double>&outx,
                             const int maxIter = 30
                             );

void Curve::InitializeFieldByLeastEigenV(bool isunified){

    if(n_vertices==0 || n_edges ==0 )return;

    //vector<double>edge_weight(n_edges,1.);
    vector<double>edge_Vweight(n_edges,0.);
    vector<double>vertices_Vweight(n_vertices,0.);
    auto edge_weight = edge_len;
    for(auto &a:edge_weight)a=1/a;

    for(int i=0;i<n_vertices;++i){
        for(auto p_ve = ve_begin(i);p_ve != ve_end(i);++p_ve)
            vertices_Vweight[i] += edge_len[*p_ve];
    }
//    for(int i=0;i<n_edges;++i){
//        for(auto p_ef = ef_begin(i);p_ef != ef_end(i);++p_ef)
//            edge_Vweight[i] += edge_len[*p_ef];
//    }
    vector<double>outx;
    EigenInitialization(    vertices_rfield,vertices_rofield,
                                 edge_weight,vertices_Vweight,edge_Vweight,
                                 edges,vertices2edges,
                                 vertices2vertices, vertices2edges_accumulate,
                                 outx
                             );



    vector<double>cosx(n_vertices);
    vector<double>sinx(n_vertices);


    //for(auto a:cosx)cout<<a<<' ';cout<<endl;
    for(int i=0;i<n_vertices;++i){
        //cosx[i] = outx[2*i];sinx[i] = outx[2*i+1];
        cosx[i] = outx[i];sinx[i] = outx[i+n_vertices];
    }

    for(int i=0;i<n_vertices;++i){
        auto p_vro = vrovec_begin(i);
        auto p_vr = vrvec_begin(i);
        auto p_vvec = vvec_begin(i);

        for(int j=0;j<3;++j)p_vvec[j] = p_vr[j]*cosx[i]+p_vro[j]*sinx[i];
        //if(isunified)normalize(p_vvec);
        if(isunified){
            projectVectorNor(p_vvec,t_begin(i),p_vvec);
            if(p_vvec[0]!=p_vvec[0])
            {
                cout<<"Nan Vertices field"<<endl;
                cout<<i<<' ';
                for(int j=0;j<3;++j)cout<< p_vr[j]*cosx[i]+p_vro[j]*sinx[i]<<' ';
                for(int j=0;j<3;++j)cout<< t_begin(i)[j]<<' ';
                for(int j=0;j<3;++j)cout<< p_vr[j]<<' ';
                cout<<endl;

                copyVec(p_vr,p_vvec);


            }
        }

    }

    BuildOthogonalField();
    //vectorfieldcheck();






}

}//n_rf
