
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

    if(n_vertices<30000)sparseSampling(1);
    else if(n_vertices<60000)sparseSampling(2);
    else sparseSampling(3);

    reinitflag = true;
    isbuilddisp = false;

    BuildDisplay(infoSurfDisp(true,false,true,false,
                                          true,false,
                                          30,17,130));
    if(outputfile.size()!=0){
        SaveInterface(outputfile);
        outputDisplay(outputfile);
    }
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


    return writeSurfFile(filename,vertices,faces2vertices,vertices_field);


}

bool Surface::outputDisplay(string filename){

    return writePLYFile(filename,display_vertices,display_faces,display_normal,display_vcolor);

}
void Surface::testSlicer(int thres){

    cout<<"thres: "<<thres<<endl;

}





}//n_rf
