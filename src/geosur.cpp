#include "geo_sur.h"
#include "readers.h"
#include<iostream>
#include<fstream>
#include <eigen3/Eigen/Geometry>


namespace n_rf {
bool Surface::isload = false;
Mesh Surface::sphere;

double randf() {
    return ((double)(rand()/(double)RAND_MAX)-0.5)*2.;
}
Surface::Surface():Mesh(){
    isbuilddisp = false;

    isresetColor = false;
    reinitflag = true;

    isFaceRenderMode = false;
    ismarkC = false;

    load();

}

void Surface::load(){
    if(!isload){
        sphere.createToy(3);
        sphere.ComputeFaceNormal(true);
        sphere.BuildNeighborTable();
        isload = true;
    }
}


void Surface::BuildFacesCenter(){
    if(n_vertices==0)return;
    if(n_faces==0)return;

    faces_center.resize(n_faces*3);
    for(uint i = 0;i<n_faces;++i){
        auto p_fv = fv_begin(i);
        _TriangleMidpoint(v_begin(p_fv[0]),v_begin(p_fv[1]),v_begin(p_fv[2]),fcent_begin(i));
    }
}
void Surface::ComputeArea(){
    if(n_vertices==0)return;
    if(n_faces==0)return;
    vertices_area.resize(n_vertices);
    faces_area.resize(n_faces);
    faces_inverse_area.resize(n_faces);
    const double scale_up = 1.;
    for(uint i = 0;i<n_faces;++i){
        auto p_fv = fv_begin(i);
        faces_area[i]=TriangleArea(p_fv[0],p_fv[1],p_fv[2]);
    }
    for(int i=0;i<n_vertices;++i){
        double b = 0;
        for(auto p_vf = vf_begin(i);p_vf != vf_end(i);++p_vf)b+=faces_area[*p_vf];
        //cout<<b<<' ';cout<<endl;
    }
    double scaleup = 1e-5;
    for(uint i = 0;i<n_faces;++i){
        faces_inverse_area[i]=scaleup/faces_area[i];
    }


    for(uint i = 0;i<n_vertices;++i){
        double a=0;
        for(auto p_vf = vf_begin(i);p_vf != vf_end(i);++p_vf){
            a+=faces_area[*p_vf];
        }
        vertices_area[i]=a*scale_up/3;
    }


}
void Surface::ComputeEdgeLength(){
    if(n_vertices==0)return;
    if(n_faces==0)return;
    if(n_edges==0)return;

    edge_len.resize(n_edges);
    edge_vec.resize(n_edges*3);
    for(uint i = 0;i<n_edges;++i){
        auto p_ev = ev_begin(i);
        vectorize(p_ev[0],p_ev[1],evec_begin(i));
        edge_len[i]=normVec(evec_begin(i));

    }

    edge_cot.resize(n_edges);
    for(uint i = 0;i<n_edges;++i){
        edge_cot[i] = EdgeCottan(i);
    }

    double ave_edge_lentmp = 0.;
    ave_edge_lentmp = 0.;
    for(auto a:edge_len)ave_edge_lentmp+=a;
    ave_edge_lentmp/=n_edges;
    ave_edge_lentmp*=2.;

    ave_edge_len = 0.;
    double numofE = 0;
    for(auto a:edge_len)if(a<ave_edge_lentmp){ave_edge_len+=a;numofE++;}
    ave_edge_len/=numofE;
    //cout<<"numofE: "<<(numofE==n_edges)<<' '<<ave_edge_len<<' '<<ave_edge_lentmp<<endl;



}



//bool Surface::readObjxfile(string filename){
//    ifstream reader(filename.data(), ifstream::in);
//    if (!reader.good()) {
//        cout << "Can not open the Objx file " << filename << endl;
//        return false;
//    }
//    clearup();
//    auto readVertices = [this](stringstream &ss){
//        double dvalue;
//        for(int i=0;i<3;++i){ss>>dvalue;vertices.push_back(dvalue);}
//    };
//    auto readFace = [this](stringstream &ss){
//        int ivalue;
//        for(int i=0;i<3;++i){
//            ss>>ivalue;faces2vertices.push_back(ivalue-1);
//        }
//    };
//    auto readVerticesField = [this](stringstream &ss){
//        double dvalue;
//        for(int i=0;i<3;++i){ss>>dvalue;vertices_field.push_back(dvalue);}
//    };


//    string oneline;

//    cout<<"reading: "<<filename<<endl;

//    while( getline( reader, oneline ))
//    {
//        stringstream ss( oneline );
//        string token;

//        ss >> token;

//        if( token == "v"  ) { readVertices( ss ); continue; } // vertex
//        if( token == "vt" ) {  continue; } // texture coordinate
//        if( token == "vn" ) {  continue; } // vertex normal
//        if( token == "vf" ) { readVerticesField( ss ); continue; } // tangent vector
//        if( token == "f"  ) { readFace( ss ); continue; } // face
//        if( token[0] == '#' ) continue; // comment
//        if( token == "o" ) continue; // object name
//        if( token == "g" ) continue; // group name
//        if( token == "s" ) continue; // smoothing group
//        if( token == "mtllib" ) continue; // material library
//        if( token == "usemtl" ) continue; // material
//        if( token == "k" ) continue; // field degree
//        if( token == "fs" ) continue; // field singularity
//        if( token == "" ) continue; // empty string
//        if( token == "vp" ) continue;// principal field, ignore here

//        cerr << "Error: does not appear to be a valid Wavefront OBJ file!" << endl;
//        cerr << "(Offending line: " << oneline << ")" << endl;
//        return false;
//    }
//    //cout<<"nfvec "<<nfvec<<endl;


//    reader.close();

//    cout<<vertices_field.size()<<' '<<vertices.size()<<' '<<faces2vertices.size()<<endl;
//    setparameters();
//    cout<<n_vertices<<' '<<n_faces<<endl;

//    auto p_fv = fv_begin(0);
//    cout<<p_fv[0]<<' '<<p_fv[1]<<' '<<p_fv[2]<<endl;

//    isfacesfield = false;
//    return true;

//}

void Surface::clearup(){

    reset();



    faces_center.clear();
    edge_len.clear();
    edge_vec.clear();
    edge_cot.clear();

    vertices_area.clear();
    faces_area.clear();
    faces_inverse_area.clear();
    face_Material1.clear();
    face_Material2.clear();

    weighted_color.clear();
    weighted_fcolor.clear();

    display_vertices.clear();
    display_normal.clear();
    display_edges.clear();

    display_vcolor.clear();
    display_faces.clear();

    invert_faces_colordegree.clear();
    invert_vertices_colordegree.clear();


    vertices_field.clear();
    vertices_rfield.clear();

    vertices_ofield.clear();
    vertices_rofield.clear();

    singularityf.clear();
    singbufferf.clear();

    singularityv.clear();
    singbufferv.clear();




    isbuilddisp = false;
    cout<<"clearup!!!!!!!!"<<endl;


}



double Surface::computeRotationAlongAxis(const double angle, const double *norAxis, const double *vec, double *vecout){

    Eigen::Vector3d axis,invec;
    for(int i=0;i<3;++i)axis(i)= norAxis[i];
    for(int i=0;i<3;++i)invec(i)= vec[i];
    normalize(axis.data());

    Eigen::AngleAxisd tb(angle,axis);
    Eigen::Vector3d outvec = tb*invec;
    for(int i=0;i<3;++i)vecout[i]= outvec(i);


}



void Surface::BuildPickID(){
    if(n_vertices==0)return;
    if(n_faces==0)return;

    id_pick_faces.resize(  n_faces*3);
    for(uint i=0;i<id_pick_faces.size();++i){
        id_pick_faces[i] = i;
    }


    id_pick_vertices.resize(n_faces*3*3);
    for(uint i=0;i<n_faces;++i){
        auto p_fv = fv_begin(i);
        auto ind = i*3*3;

        for(uint j=0;j<3;++j){
            auto p_v = v_begin(p_fv[j]);
            for(int k=0;k<3;++k)id_pick_vertices[ind+j*3+k] = p_v[k];
        }

    }


    id_pick_color.resize(n_faces*3*3*4);
    unsigned char cc[3];

    for(uint i=0;i<n_faces;++i){

        cc[0] = ( i & 0xFF );
        cc[1] = ( i & 0xFF00 ) >> 8;
        cc[2] = ( i & 0xFF0000 ) >> 16;
        cc[3] = ( i & 0xFF000000 ) >> 24;
        auto ind = i*3*4;

        for(uint j=0;j<3;++j){
            id_pick_color[ind+j*4]=cc[0];
            id_pick_color[ind+j*4+1]=cc[1];
            id_pick_color[ind+j*4+2]=cc[2];
            id_pick_color[ind+j*4+3]=cc[3];
        }


    }

    cout<<"BuildPickID!!!!"<<endl;

}





int Surface::PickFaceViaColor(unsigned char* pcolor){


    int pickIndex = pcolor[0] | (pcolor[1]<<8) | (pcolor[2]<<16) | (pcolor[4]<<24);



    cout<<"PickFaceViaColor: "<<pickIndex<<endl;
    if(pickIndex>=n_faces)pickIndex=-1;
    return pickIndex;

}


bool Surface::readSufFile(string filename){

    ifstream reader(filename.data(), ifstream::in);
    if (!reader.good()) {
        cout << "Can not open the Suf file " << filename << endl;
        return false;
    }
    reset();

    cout<<"Reading Suf File"<<endl;



    reader>>n_vertices;
    reader>>n_faces;

    vertices.resize(n_vertices*3);
    face_Material1.clear();
    face_Material2.clear();
    faces2vertices.clear();

    for(int i =0;i<vertices.size();i++){
        reader>>vertices[i];
    }

    int ivalue;
    for(int i =0;i<n_faces;i++){
        for(int j =0;j<3;++j){reader>>ivalue;faces2vertices.push_back(ivalue);}
        reader>>ivalue;face_Material1.push_back(ivalue);
        reader>>ivalue;face_Material2.push_back(ivalue);
    }

    vector<bool>mmmmm(32,false);
    for(int i =0;i<n_faces;i++){
        mmmmm[face_Material1[i]] = true;
    }
    int numofM = 0;
    for(int i =0;i<32;i++){
        if(!mmmmm[i]){numofM = i;break;}
    }

    vector< vector<uint> >mm_faces(numofM);
    for(int i =0;i<n_faces;i++){
        int m1 = face_Material1[i];
        int m2 = face_Material2[i];
        auto p_fv = fv_begin(i);
        for(int j=0;j<3;++j){mm_faces[m1].push_back(p_fv[j]);}
        for(int j=0;j<3;++j){mm_faces[m2].push_back(p_fv[j]);}


    }

    faces2vertices = mm_faces[0];
    n_faces = faces2vertices.size()/3;
    cout<<n_vertices<<n_faces<<endl;




}




void Surface::BuildOthogonalField(){
    if(n_vertices==0)return;
    if(n_faces==0)return;


    if(vertices_field.size()!=0){
        vertices_ofield.resize(vertices_field.size());
        for(int i =0;i<n_vertices;++i)cross(vnor_begin(i),vvec_begin(i),vovec_begin(i));
    }

    if(vertices_rfield.size()!=0){
        vertices_rofield.resize(vertices_rfield.size());
        for(int i =0;i<n_vertices;++i)cross(vnor_begin(i),vrvec_begin(i),vrovec_begin(i));
    }
}
void Surface::InitializeField(bool isProj,bool isrand){

    if(n_vertices==0)return;
    if(n_faces==0)return;

    double iiv[3] = {1/sqrt(3),1/sqrt(3),0};

    vertices_rfield.resize(n_vertices*3);
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


    BuildOthogonalField();
    //ComputeVerticesFieldSingularity();
    //ComputeFacesFieldSingularity();
    //rvectorVMFConstrain = vertices_rfield;

    //cout<<"Face Energy: "<<DualEdgesDifferetialEnergy()<<endl;
    //cout<<"Vertices Energy: "<<EdgesDifferetialEnergy()<<endl;
    //    faces_rfield = faces_field;
    //    faces_rofield = faces_ofield;


}

int Surface::ComputeVerticesFieldSingularity(){

    if(n_vertices==0)return -1;
    if(n_faces==0)return -1;
    if(vertices_field.size()==0)return -1;


    singularityf.resize(n_faces,false);
    singbufferf.clear();



    vector<double>faces_error(n_faces,0);

    double vec1[9],err[3];

    for(uint i =0;i<n_faces;++i){
        auto p_fv = fv_begin(i);
        //cout<<p_fv[0]<<' '<<p_fv[1]<<' '<<p_fv[2]<<endl;
        for(int j=0;j<3;++j)projectVectorNor(vvec_begin(p_fv[j]),fnor_begin(i),vec1+j*3);
        err[0] = angleNor(vec1+0,vec1+3);
        err[1] = angleNor(vec1+3,vec1+6);
        err[2] = angleNor(vec1+6,vec1+0);
        for(int j=0;j<3;++j)faces_error[i]+=err[j];
    }



    for(uint i =0;i<n_faces;++i){
        if(faces_error[i]>1.90*PI){
            //singularityf[i]=true;
            singbufferf.push_back(i);
        }
    }

    //return;
    //for(auto a:singbuffer)cout<<a<<' ';cout<<endl;

    //for(auto a:singbuffer)cout<<ve_num(a)<<' ';cout<<endl;
    //cout<<"Faces singularity: "<<singbufferf.size()<<endl;
    vector<bool> duplicateF(singbufferf.size(),false);
    for(int i =0;i<singbufferf.size();++i){
        if(duplicateF[i])continue;
        auto curF = singbufferf[i];
        for(auto p_ff = ff_begin(curF);p_ff!=ff_end(curF);++p_ff){
            auto f = *p_ff;
            for(int j =0;j<singbufferf.size();++j)if(f==singbufferf[j])duplicateF[j]=true;
        }
    }

    auto temp = singbufferf;
    singbufferf.clear();
    for(int j =0;j<temp.size();++j)if(!duplicateF[j])singbufferf.push_back(temp[j]);
    for(auto a:singbufferf)singularityf[a]=true;
    cout<<"Faces singularity: "<<singbufferf.size()<<endl;


    return singbufferf.size();

}

void Surface::meshScalarfunctionSmoothing(vector<double>& scalarF,int maxiter){

    vector<double>buffer1(scalarF.size());
    vector<double>buffer2(scalarF.size());
    vector<double>*pre = &buffer1,*cur = &buffer2;
    *pre = scalarF;

    double kaipa = 0.1,lamda = 0.5;
    auto gama = 1 / (kaipa - 1 / lamda);

    auto oneiter = [this](double coef,vector<double>*pre,vector<double>*cur){
        double ocoef = 1- coef;
        for(int i =0;i<n_vertices;++i){
            double aaa = 0;
            for(auto p_vv = vv_begin(i);p_vv!=vv_end(i);++p_vv){
                aaa+=pre->at(*p_vv);
            }
            aaa/=(vv_num(i));
            cur->at(i) = ocoef * pre->at(i) + coef * aaa;
        }
    };

    for(int iter=0;iter<maxiter;++iter){

        oneiter(lamda,pre,cur);
        swap(pre,cur);
        oneiter(gama,pre,cur);
        swap(pre,cur);

    }
    scalarF=(*pre);





}



double GaussSeidelIteration(const vector<double>&vertices_normal,
                            const vector<double>&edge_weight,
                            const vector<uint>&edges2vertices,const vector<uint>&vertices2edges,
                            const vector<uint>&vertices2vertices, const vector< unsigned long >&vertices2edges_accumulate,
                            vector<double>&vertices_field,
                            const int maxiter = 5000
                            );

double Surface::RunGaussSeidelIteration(){



    double en_pre = GaussSeidelIteration(vertices_normal,edge_cot,edges,vertices2edges,vertices2vertices,vertices2edges_accumulate,vertices_field);


    return en_pre;

}



void EigenInitialization(    const vector<double>&vertices_rfield,const vector<double>&vertices_rofield,
                             const vector<double>&edge_weight,const vector<double>&vertices_Vweight,const vector<double>&edge_Vweight,
                             const vector<uint>&edges2vertices,const vector<uint>&vertices2edges,
                             const vector<uint>&vertices2vertices, const vector< unsigned long >&vertices2edges_accumulate,
                             vector<double>&outx
                             );

void Surface::InitializeFieldByLeastEigenV(bool isunified){
    if(n_vertices==0 || n_edges ==0 || n_faces==0)return;

    //vector<double>edge_weight(n_edges,1.);
    vector<double>edge_Vweight(n_edges,0.);
    vector<double>vertices_Vweight(n_vertices,0.);
    //    for(int i=0;i<n_edges;++i){
    //        edge_weight[i] = EdgeCottan(i);
    //    }
    for(int i=0;i<n_vertices;++i){
        for(auto p_vf = vf_begin(i);p_vf != vf_end(i);++p_vf)
            vertices_Vweight[i] += faces_area[*p_vf];
    }
    for(int i=0;i<n_edges;++i){
        for(auto p_ef = ef_begin(i);p_ef != ef_end(i);++p_ef)
            edge_Vweight[i] += faces_area[*p_ef];
    }
    vector<double>outx;
    EigenInitialization(    vertices_rfield,vertices_rofield,
                            edge_cot,vertices_Vweight,edge_Vweight,
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
    //vectorfieldcheck();

    int a = ComputeVerticesFieldSingularity();
    //cout<<"Edge Different energy: "<<EdgesDifferetialEnergy()/2.<<endl;
    cout<<"Singularity: "<<a<<endl;


}


}//n_rf
