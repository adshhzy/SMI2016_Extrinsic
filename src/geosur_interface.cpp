
#include "geo_sur.h"
#include "readers.h"
namespace n_rf {


int Surface::Initialize(string inputfile, string outputfile, bool isEigenInit, bool isGaussIter){
    clearup();
    bool rere = false;

    rere = ReadFile(inputfile);

    if(!rere)return -1;

    setparameters();
    ReScale_uniform(1.0);
    isbuilddisp = false;
    cout<<"read finish!"<<endl;


    BuildNeighborTable();
    //ReOrientFaces();
    ComputeEdgefromFace();

    ComputeEdgeLength();
    ComputeDefectAngle();
    ComputeFaceNormal();

    BuildFacesCenter();
    ComputeArea();


    BuildPickID();


    InitializeField(false,true);
    cout<<"Mesh Initialized!"<<endl;

    if(isEigenInit)InitializeFieldByLeastEigenV();
    if(isGaussIter)RunGaussSeidelIteration();
    ComputeVerticesFieldSingularity();

    sparseSampling(3);

    reinitflag = true;
    isbuilddisp = false;

    if(outputfile.size()!=0)SaveInterface(outputfile);
    return 0;




}






bool Surface::ReadFile(string filename){


    string filepath;
    string modname;
    string extname;

    SplitFileName(filename,filepath,modname,extname);
    cout<<modname<<' '<<extname<<' '<<filepath<<endl;

    bool issuc = false;

    if(extname==".suf"){issuc = readSufFile(filename);}
    else issuc = readfile(filename);

    if(issuc){
        modelname = modname;
        prepath = filepath;
        ext=extname;
    }

    return issuc;
}




bool Surface::SaveInterface(string filename){

//    string prepath,modelname,ext;
//    SplitFileName(filename,prepath,modelname,ext);

//    if(ext==".obj")saveObjFile(filename);
//    else if(ext==".off")writeOffFile(filename,vertices,faces2vertices);
//    else if(ext==".surf")writeSurfFile(filename,vertices,faces2vertices,vertices_field);

//    return true;
    writeSurfFile(filename,vertices,faces2vertices,vertices_field);
    return true;

}

void Surface::testSlicer(int thres){

    cout<<"thres: "<<thres<<endl;

}





}//n_rf
