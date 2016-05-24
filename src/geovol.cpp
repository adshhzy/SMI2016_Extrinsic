#include "geo_vol.h"
#include<set>






namespace n_rf {
bool Volume::isload = false;
Mesh Volume::sphere;
Mesh Volume::cylinder;
Mesh Volume::cone;

double randf();
Volume::Volume(){

    isbuilddisp = false;
    sparsecoff = 2;
    load();

}

void Volume::setparameters(){

    n_vertices=vertices.size()/3;
    n_tetr = tetrahedron2vertices.size()/4;
}


void Volume::clearup(){


    n_vertices=n_edges=n_tetr=0;
    vertices.clear();
    edges.clear();
    tetrahedron2vertices.clear();


    vertices_normal.clear();

    vertices2edges.clear();
    vertices2vertices.clear();
    vertices2edges_accumulate.clear();
    vertices2edgesPN.clear();
    vertices2vertices_inverse.clear();


    vertices2tetr.clear();
    vertices2tetr_accumulate.clear();
    vertices2tetr_inverse.clear();



    edges2tetr.clear();
    edges2tetr_accumulate.clear();
    edges2tetr_inverse.clear();




    display_vertices.clear();
    display_normal.clear();
    display_edges.clear();
    display_vcolor.clear();
    display_faces.clear();




    vertices_field.clear();
    vertices_rfield.clear();



    vertices_ofield.clear();
    vertices_rofield.clear();
    sparseShowField.clear();

    tetrsvolume.clear();
    edge_cot_weight.clear();
    edge_length.clear();
    edge_M_weight.clear();
    vertices_M_weight.clear();


}


void Volume::BuildEdgesRelations(){
    if(n_vertices==0)return;
    if(n_tetr==0)return;


    vertices2tetr_accumulate.resize(n_vertices+1,0);
    vector<uint>vertices2tetrnum(n_vertices,0);

    for(uint i=0;i<n_tetr;++i){
        auto p_tv = tv_begin(i);
        for(int j=0;j<4;++j){
            vertices2tetrnum[p_tv[j]]++;
        }
    }

    vertices2tetr_accumulate[0]=0;
    for(uint i=0;i<n_vertices;++i){
        vertices2tetr_accumulate[i+1] = vertices2tetr_accumulate[i]+vertices2tetrnum[i];
    }

    vector<uint>curind(n_vertices,0);
    vertices2tetr.resize(vertices2tetr_accumulate[n_vertices]);
    vertices2tetr_inverse.resize(vertices2tetr_accumulate[n_vertices]);
    for(uint i=0;i<n_tetr;++i){
        auto p_tv = tv_begin(i);
        for(int j=0;j<4;++j){
            auto ind = p_tv[j];
            vt_begin(ind)[curind[ind]] = i;
            vtinv_begin(ind)[curind[ind]] = j;
            ++curind[ind];
        }
    }


    vector<bool>ispick(n_vertices,false);
    vector<uint>curbuffer;
    vertices2vertices.clear();
    vertices2vertices.reserve(n_vertices*5);
    vertices2edges_accumulate.resize(n_vertices+1);
    vertices2edges_accumulate[0]=0;
    for(uint i=0;i<n_vertices;++i){
        curbuffer.clear();
        ispick[i]=true;
        auto vtnum = vt_num(i);
        auto p_vt=vt_begin(i);
        for(int j=0;j<vtnum;++j){
            auto tind = p_vt[j];
            auto p_tv = tv_begin(tind);
            for(int k=0;k<4;++k){
                if(!ispick[p_tv[k]]){
                    curbuffer.push_back(p_tv[k]);
                    ispick[p_tv[k]] = true;
                }
            }
        }

        vertices2edges_accumulate[i+1] = vertices2edges_accumulate[i]+curbuffer.size();
        for(auto a:curbuffer)vertices2vertices.push_back(a);
        for(auto a:curbuffer)ispick[a]=false;
        ispick[i]=false;

    }

    vertices2vertices_inverse.resize(vertices2edges_accumulate[n_vertices]);

    for(int i=0;i<n_vertices;++i){
        auto p_vvinv = vvinv_begin(i);
        auto p_vvend = vv_end(i);
        for(auto p_vv = vv_begin(i);p_vv!=p_vvend;++p_vv,++p_vvinv){
            auto p_vvi = vv_begin(*p_vv);
            auto num = vv_num(*p_vv);
            for(int j=0;j<num;++j)if(p_vvi[j]==i){
                (*p_vvinv)=j;break;
            }
        }
    }

    vertices2edges.resize(vertices2edges_accumulate[n_vertices]);
    vertices2edgesPN.resize(vertices2edges_accumulate[n_vertices]);
    for(auto &a:vertices2edges)a=UINTFLAG;
    edges.clear();
    edges.reserve(n_vertices*2);
    uint eind = 0;
    for(int i=0;i<n_vertices;++i){
        auto p_vvinv = vvinv_begin(i);
        auto p_vvend = vv_end(i);
        auto p_ve = ve_begin(i);
        auto p_vepn = vepn_begin(i);
        for(auto p_vv = vv_begin(i);p_vv!=p_vvend;++p_vv,++p_vvinv,++p_ve,++p_vepn){
            if((*p_ve)!=UINTFLAG){continue;}
            auto p_vei = ve_begin(*p_vv);
            auto p_vepni = vepn_begin(*p_vv);
            p_vei[*p_vvinv] = (*p_ve) = eind;
            ++eind;
            edges.push_back(i);edges.push_back(*p_vv);
            (*p_vepn) = true;(*p_vepni) = false;
        }
    }
    n_edges = eind;


    edges2tetr.clear();
    edges2tetr_accumulate.clear();
    edges2tetr_inverse.clear();

    edges2tetr_accumulate.resize(n_edges+1,0);
    vector<bool>tetrpick(n_tetr,false);

    //cout<<"adsadsadasdasda"<<endl;
    for(int i=0;i<n_edges;++i){
        unsigned long acc = 0;
        auto p_ev =ev_begin(i);
        auto vtende0 = vt_end(p_ev[0]);
        auto vtende1 = vt_end(p_ev[1]);
        for(auto p_vt = vt_begin(p_ev[0]);p_vt!=vtende0;++p_vt)tetrpick[*p_vt] = true;

        for(auto p_vt = vt_begin(p_ev[1]);p_vt!=vtende1;++p_vt){
            if(tetrpick[*p_vt]){
                edges2tetr.push_back(*p_vt);++acc;
            }
        }
        for(auto p_vt = vt_begin(p_ev[0]);p_vt!=vtende0;++p_vt)tetrpick[*p_vt] = false;
        edges2tetr_accumulate[i+1] = acc;
        //cout<<acc<<' ';
    }
    for(int i=1;i<=n_edges;++i)edges2tetr_accumulate[i] = edges2tetr_accumulate[i] + edges2tetr_accumulate[i-1];



    //for(int i=1;i<=n_edges;++i)cout<<et_num(i)<<' ';cout<<endl;
    //int aadsa = 0;
    //for(int i=0;i<n_edges;++i)if(et_num(i)==6)aadsa++;
    //cout<<"aadsa "<<aadsa<<endl;
    cout<<"number of vertices: "<<n_vertices<<endl;
    cout<<"number of tetr: "<<n_tetr<<endl;
    cout<<"number of edges: "<<n_edges<<endl;



}


int Volume::Initialize(string inputfile, string outputfile, bool isEigenInit, bool isGaussIter){


    //    infoVolDisp(bool isField,bool isNormal,bool isSurface,bool isWire,bool isSingularity,
    //                bool isIso,int length,int upnormal,
    //                bool isSlice,int xyz,int slicenum,double thickness,bool isLineMode);
    //IsoToyGeneration(30,30,30,1.4,1.4,0.4);
    if(!ReadInterface(inputfile)){

        cout<<"invalid input: "<< inputfile <<endl;
        exit(-1);
    }
    BuildEdgesRelations();
    ComputeEigenWeight();
    InitializeField();
    if(isEigenInit)InitializeFieldByLeastEigenV();
    if(isGaussIter)RunGaussSeidelIteration();
    sparseSampling(2);
    //BuildDisplay(infoVolDisp(true,false,false,false,false,false,30,0,false,0,0,120,false));
    BuildDisplay(infoVolDisp(30,70,30));
    if(outputfile.size()!=0){
        SaveInterface(outputfile);
        outputDisplay(outputfile);
    }
    return 0;

}
int Volume::SpecialCases(string outputfile, bool isEigenInit, bool isGaussIter){
    IsoToyGeneration(30,30,30,1.4,1.4,0.4);
    //IsoToyGeneration(3,3,3,1.,1.,1.);
    BuildEdgesRelations();
    ComputeEigenWeight();

    InitializeField();
    if(isEigenInit)InitializeFieldByLeastEigenV();
    if(isGaussIter)RunGaussSeidelIteration();
    sparseSampling(2);
    //BuildDisplay(infoVolDisp(false,false,false,false,false,false,30,0,false,0,0,120,false));
    BuildDisplay(infoVolDisp(30,70,30));
    if(outputfile.size()!=0){
        SaveInterface(outputfile);
        outputDisplay(outputfile);
    }

    return 0;
}

void Volume::ToyGeneration(){

//    const int numofLayer = 8;

//    const double startscale = 0.6;


//    Mesh sphereRefine;

//    double scalestep = (1.-startscale)/numofLayer;
//    n_vertices = sphereRefine.n_vertices*numofLayer;
//    vertices.resize(n_vertices*3);
//    auto p_v = vertices.data();

//    for(int i=0;i<numofLayer;++i){
//        sphereRefine.ReScale_uniform(startscale+scalestep*(i+1));
//        memcpy((void*)(p_v+sphereRefine.n_vertices*i*3),(void*)(sphereRefine.vertices.data()),sphereRefine.n_vertices*3*sizeof(double));

//    }


//    n_tetr = (numofLayer-1)*(sphereRefine.n_faces*3);
//    tetrahedron2vertices.resize(n_tetr*4);

//    auto p_tv = tetrahedron2vertices.data();
//    int layerjump = sphereRefine.n_vertices;
//    for(int j=0;j<sphereRefine.n_faces;++j){
//        auto ind = j*3*4;
//        auto p_fv = sphereRefine.fv_begin(j);
//        p_tv[ind] = p_fv[0];p_tv[ind+1] = p_fv[1];p_tv[ind+2] = p_fv[2];p_tv[ind+3] = p_fv[0]+layerjump;
//        p_tv[ind+4] = p_fv[0]+layerjump;p_tv[ind+5] = p_fv[1]+layerjump;p_tv[ind+6] = p_fv[2]+layerjump;p_tv[ind+7] = p_fv[1];
//        p_tv[ind+8] = p_fv[1];p_tv[ind+9] = p_fv[2];p_tv[ind+10] = p_fv[0]+layerjump;p_tv[ind+11] = p_fv[2]+layerjump;

//    }



//    for(int i=1;i<numofLayer-1;++i){
//        auto p_tvn = p_tv+(sphereRefine.n_faces*3)*4*i;
//        uint jumpjump = layerjump*i;
//        for(int j=0;j<(sphereRefine.n_faces*3)*4;++j) p_tvn[j] = p_tv[j]+jumpjump;

//    }

//    //cout<<"000: "<<vertices.size()<<endl;

//    uchar colllor[4] = {250,0,0,200};
//    display_vcolor.resize(n_vertices*4);

//    //cout<<"000: "<<vertices.size()<<endl;

//    auto p_color = display_vcolor.data();
//    cout<<"p_color: "<<p_color<<endl;
//    for(int i=0;i<numofLayer;++i){
//        auto p_colorc = p_color+sphereRefine.n_vertices*4*i;
//        colllor[3] = 255-250/numofLayer*i;
//        colllor[0] = 255-250/numofLayer*i;
//        colllor[1] = 250/numofLayer*i;
//        for(int j=0;j<sphereRefine.n_vertices;++j){
//            auto mmm=j*4;
//            for(int k=0;k<4;++k)p_colorc[mmm+k] = colllor[k];

//        }

//    }

//    //cout<<"000: "<<vertices.size()<<endl;
//    vertices_normal.resize(n_vertices*3);

//    sphereRefine.ComputeFaceNormal(true);
//    for(int i=0;i<numofLayer;++i){
//        memcpy((void*)(vertices_normal.data()+sphereRefine.n_vertices*i*3),(void*)(sphereRefine.vertices_normal.data()),sphereRefine.n_vertices*3*sizeof(double));

//    }


//    sphereRefine.ComputeEdgefromFace(true);
//    //edges = sphereRefine.edges;

//    edges.resize(sphereRefine.n_edges*2*numofLayer);

//    auto p_eee = sphereRefine.edges.data();
//    for(int i=0;i<numofLayer;++i){
//        auto p_e = edges.data()+sphereRefine.n_edges*2*i;
//        auto jumpjump = sphereRefine.n_vertices*i;
//        for(int j=0;j<sphereRefine.n_edges*2;++j)p_e[j] = p_eee[j]+jumpjump;
//    }

//    cout<<"n_tetr: "<<n_tetr<<" n_vertices: "<<n_vertices<<endl;


}

void Volume::InitializeField(bool isProj,bool isrand){

    if(n_vertices==0)return;
    if(n_tetr==0)return;

    double iiv[3] = {1/sqrt(3),1/sqrt(3),0};

    double riiv[3];
    vertices_rfield.resize(n_vertices*3);
    for(uint i = 0;i<n_vertices;++i){
        auto p_vr = vrvec_begin(i);
        if(isrand)for(int i =0;i<3;++i)iiv[i]= randf();
        projectVectorNor(iiv,vnor_begin(i),p_vr);
        while(p_vr[0]!=p_vr[0] || p_vr[1]!=p_vr[1] || p_vr[2]!=p_vr[2]){
            for(int j =0;j<3;++j)riiv[j]= randf();
            projectVectorNor(iiv,vnor_begin(i),p_vr);
        }

    }



    if(isProj || vertices_field.size()==0)
    {
        vertices_field = vertices_rfield;

    }else{
        for(uint i = 0;i<n_vertices;++i){
            projectVectorNor(vvec_begin(i),vnor_begin(i),vvec_begin(i));
        }
    }

    BuildOthogonalField(true);

}

void Volume::BuildOthogonalField(bool isR){
    if(isR && vertices_rfield.size()!=0){
        vertices_rofield.resize(vertices_rfield.size());
        for(int i =0;i<n_vertices;++i)cross(vnor_begin(i),vrvec_begin(i),vrovec_begin(i));
    }

    if(vertices_field.size()!=0){
        vertices_ofield.resize(vertices_field.size());
        for(int i =0;i<n_vertices;++i)cross(vnor_begin(i),vvec_begin(i),vovec_begin(i));
    }

}
void Volume::RandomInitialize(int Indmethod){

    if(n_vertices==0)return;
    if(n_tetr==0)return;

    srand(time(NULL));
    double iiv[3];

    auto randomConstProjVertices = [this,&iiv](){
        for(int i =0;i<3;++i)iiv[i]=randf();
        normalize(iiv);
        for(uint i = 0;i<n_vertices;++i){
            projectVectorNor(iiv,vnor_begin(i),vvec_begin(i));
        }
    };
    auto randomProjVertices = [this,&iiv](){
        for(uint i = 0;i<n_vertices;++i){
            for(int j =0;j<3;++j)iiv[j]=randf();
            projectVectorNor(iiv,vnor_begin(i),vvec_begin(i));
        }
    };

    vertices_field.resize(n_vertices*3);
    switch(Indmethod){
    case 0:randomProjVertices();break;
    case 1:randomConstProjVertices();break;

    default:break;
    }


    BuildOthogonalField();


}


double Volume::EdgesDifferetialEnergy(){
    double en = 0;
    for(uint i=0;i<n_edges;++i){
        auto p_ev = ev_begin(i);
        en+=pow(_VerticesDistance(vvec_begin(p_ev[0]),vvec_begin(p_ev[1])),2);
    }
    return en;
}


void Volume::Rescale_uniform(double lscale){

    double xmin = 999e10, ymin=999e10,zmin=999e10,xmax=-999e10,ymax=-999e10,zmax=-999e10;
    double centers[3] = {0,0,0};
    for (int i = 0; i<n_vertices; i++){
        auto point = v_begin(i);
        xmin = min(xmin,point[0]);xmax = max(xmax,point[0]);
        ymin = min(ymin,point[1]);ymax = max(ymax,point[1]);
        zmin = min(zmin,point[2]);zmax = max(zmax,point[2]);
        centers[0]+=point[0];centers[1]+=point[1];centers[2]+=point[2];
        //if(centers[0]!=centers[0]){cout<<centers[0]<<"  "<<point[0]<<" "<<i/dim<<endl;}
    }
    centers[0]/=n_vertices;centers[1]/=n_vertices;centers[2]/=n_vertices;
    double largestdis = max(max(xmax-xmin,ymax-ymin),zmax-zmin)/2 / lscale;
    //cout<<centers[0]<<"  "<<centers[1]<<"  "<<centers[2]<<"   "<<largestdis<<endl;

    largestdis = 1/largestdis;
    for (int j = 0; j<n_vertices*3; j+=3)
        for(int i = 0; i<3; i++)
            vertices[j+i]=(vertices[j+i]-centers[i])*largestdis;


}


double GaussSeidelIteration(const vector<double>&vertices_normal,
                            const vector<double>&edge_weight,
                            const vector<uint>&edges2vertices,const vector<uint>&vertices2edges,
                            const vector<uint>&vertices2vertices, const vector< unsigned long >&vertices2edges_accumulate,
                            vector<double>&vertices_field,
                            const int maxiter = 5000
        );

double Volume::RunGaussSeidelIteration(){

    double en_pre = GaussSeidelIteration(vertices_normal,edge_cot_weight,edges,vertices2edges,vertices2vertices,vertices2edges_accumulate,vertices_field);

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

void Volume::InitializeFieldByLeastEigenV(bool isunified){
    if(n_vertices==0 || n_edges ==0 || n_tetr==0)return;

    vector<double>edge_weight(n_edges,1.);
    vector<double>edge_Vweight(n_edges,0.);
    vector<double>vertices_Vweight(n_vertices,0.);
    //    for(int i=0;i<n_edges;++i){
    //        edge_weight[i] = EdgeCottan(i);
    //    }
    for(int i=0;i<n_vertices;++i){
        for(auto p_vt = vt_begin(i);p_vt != vt_end(i);++p_vt)
            vertices_Vweight[i] += 1.;
    }
    //    for(int i=0;i<n_edges;++i){
    //        for(auto p_et = et_begin(i);p_ef != ef_end(i);++p_ef)
    //            edge_Vweight[i] += faces_area[*p_ef];
    //    }
    vector<double>outx;
    EigenInitialization(    vertices_rfield,vertices_rofield,
                            edge_cot_weight,vertices_M_weight,edge_M_weight,
                            edges,vertices2edges,
                            vertices2vertices, vertices2edges_accumulate,
                            outx,20
                            );



    vector<double>cosx(n_vertices);
    vector<double>sinx(n_vertices);


    for(int i=0;i<n_vertices;++i){

        cosx[i] = outx[i];sinx[i] = outx[i+n_vertices];
    }

    for(int i=0;i<n_vertices;++i){
        auto p_vro = vrovec_begin(i);
        auto p_vr = vrvec_begin(i);
        auto p_vvec = vvec_begin(i);

        for(int j=0;j<3;++j)p_vvec[j] = p_vr[j]*cosx[i]+p_vro[j]*sinx[i];
        //if(isunified)normalize(p_vvec);
        if(isunified){
            projectVectorNor(p_vvec,vnor_begin(i),p_vvec);
            if(p_vvec[0]!=p_vvec[0])
            {
                cout<<"Nan Vertices field"<<endl;
                cout<<i<<' ';
                for(int j=0;j<3;++j)cout<< p_vr[j]*cosx[i]+p_vro[j]*sinx[i]<<' ';
                for(int j=0;j<3;++j)cout<< vnor_begin(i)[j]<<' ';
                for(int j=0;j<3;++j)cout<< p_vr[j]<<' ';
                cout<<endl;

                copyVec(p_vr,p_vvec);


            }
        }

    }

    BuildOthogonalField();



}


void Volume::ComputeEigenWeight(){
    if(n_vertices==0 || n_tetr ==0)return;

    tetrsvolume.resize(n_tetr);
    for(int i=0;i<n_tetr;++i){
        auto p_tv = tv_begin(i);
        tetrsvolume[i] = _TetrahedronVolume(v_begin(p_tv[0]),v_begin(p_tv[1]),v_begin(p_tv[2]),v_begin(p_tv[3]));
    }
    edge_length.resize(n_edges,0);
    double avelen = 0;
    for(int i=0;i<n_edges;++i){
        auto p_ev = ev_begin(i);
        edge_length[i] = _VerticesDistance(v_begin(p_ev[0]),v_begin(p_ev[1]));
    }
    for(auto a:edge_length)avelen+=a;
    avelen/=n_edges;
    for(auto &a:edge_length)a/=avelen;
    edge_cot_weight.resize(n_edges,0);
    vector<bool>vpickbuffer(n_vertices,false);

    for(int i=0;i<n_edges;++i){
        auto p_ev = ev_begin(i);
        vpickbuffer[p_ev[0]] = vpickbuffer[p_ev[1]] = true;
        auto p_et = et_begin(i);
        auto etnum = et_num(i);
        double ecw = 0;
        uint rightvv[2];
        for(int j=0;j<etnum;++j){
            int currrv = 0;
            auto p_tv = tv_begin(p_et[j]);
            for(int k=0;k<4;++k){
                if(!vpickbuffer[p_tv[k]]){
                    rightvv[currrv] = p_tv[k];
                    ++currrv;
                }
            }
            ecw+=_cotBetweenTwoPlanes(v_begin(p_ev[0]),v_begin(p_ev[1]),v_begin(rightvv[0]),v_begin(rightvv[1]));
            //ecw+=_angleBetweenTwoPlanes(v_begin(p_ev[0]),v_begin(p_ev[1]),v_begin(rightvv[0]),v_begin(rightvv[1]));
        }
        edge_cot_weight[i] = ecw*edge_length[i]/2.;

        vpickbuffer[p_ev[0]] = vpickbuffer[p_ev[1]] = false;

    }
    //for(auto a:edge_cot_weight)cout<<a<<' ';cout<<endl;

    edge_M_weight.resize(n_edges,0);

    for(int i=0;i<n_edges;++i){
        auto p_et = et_begin(i);
        auto etnum = et_num(i);
        double emw = 0;
        for(int j=0;j<etnum;++j){
            emw += tetrsvolume[p_et[j]];
        }

        edge_M_weight[i] = emw;

    }
    vertices_M_weight.resize(n_vertices,0);

    for(int i=0;i<n_vertices;++i){
        auto p_vt = vt_begin(i);
        auto vtnum = vt_num(i);
        double vmw = 0;
        for(int j=0;j<vtnum;++j){
            vmw += tetrsvolume[p_vt[j]];
        }

        vertices_M_weight[i] = vmw;

    }
}


}//n_rf





