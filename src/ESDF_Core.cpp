#include "ESDF_Core.h"
#include<sys/time.h>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/IterativeLinearSolvers>
//#include <eigen3/Eigen/CholmodSupport>
//#include <eigen3/Eigen/SparseCholesky>
//#include <eigen3/Eigen/SparseLU>
//#include <eigen3/Eigen/SparseQR>
#include<unistd.h>
namespace n_rf {


int ESDF_Core::ReadArgs(int argc,char** argv){

    bool isGaussIter = false;
    bool isEigenInit = false;

    int isSpecialCase = -1;

    string inputfile;
    string outfile = "output";

    int c;
    optind=1;
    while ((c = getopt(argc, argv, "egs:i:o:")) != -1) {
        switch (c) {
        case 'g':
            isGaussIter = true;
            break;
        case 'e':
            isEigenInit = true;
            break;
        case 'i':
            inputfile = optarg;
            break;
        case 's':
            isSpecialCase = atoi(optarg);
            break;
        case 'o':
            outfile = optarg;
            break;
        case '?':
            cout << "Bad argument setting!" << endl;
            break;
        }
    }



    if(isSpecialCase>=0){
        cout << "Special Case Activated! " << endl;
        cout << "Output File: "<<outfile<<endl;
        cout << "Eigenvector initialization Activation " << isEigenInit << endl;
        cout << "Gauss-Seidel iterations Activation " << isGaussIter << endl;

        if(isSpecialCase<=2){
            Curv_core.SpecialCases(isSpecialCase,inputfile,outfile,isEigenInit,isGaussIter);
            curstate = CURV;return 0;
        }else if(isSpecialCase==3){
            Vol_core.SpecialCases(outfile,isEigenInit,isGaussIter);
            curstate = VOL;
            return 0;
        }
        cout<<"Accepted special case: 0,1,2,3, pleas check."<<endl;
        return -1;

    }else{
        if(inputfile.size()==0){cout << "Invalid Input! " << endl; }
        else{cout << "Input File: "<<inputfile<<endl;}

        cout << "Output File: "<<outfile<<endl;
        cout << "Eigenvector initialization Activation " << isEigenInit << endl;
        cout << "Gauss-Seidel iterations Activation " << isGaussIter << endl;



        SplitFileName(inputfile,prepath,modelname,ext);

        if(ext==".obj" || ext==".off"){
            Surf_core.Initialize(inputfile,outfile,isEigenInit,isGaussIter);
            curstate = SURF;return 0;
        }else if(ext==".curf"){
            Curv_core.Initialize(inputfile,outfile,isEigenInit,isGaussIter);
            curstate = CURV;return 0;
        }else if(ext==".curf"){
            Curv_core.Initialize(inputfile,outfile,isEigenInit,isGaussIter);
            curstate = CURV;return 0;
        }else if(ext == ".volf"){
            Vol_core.Initialize(inputfile,outfile,isEigenInit,isGaussIter);
            curstate = VOL;return 0;
        }



        cout<<"The format of the input file is not support!"<<endl;
    }
    return 0;


}


bool ESDF_Core::clearup(){

    Curv_core.reset();
    Surf_core.clearup();
    Vol_core.clearup();
    return true;

}


void ESDF_Core::BuildDisplay(infoSurfDisp info){

    if(curstate == SURF)Surf_core.BuildDisplay(info);
    else if(curstate == CURV)return Curv_core.BuildDisplay(infoSet(info.length,info.width,
                                                                   info.length,info.upnormal));
    else if(curstate == VOL)return Vol_core.BuildDisplay(infoVolDisp(info.length,info.upnormal,info.width));

}

vector<double>* ESDF_Core::getDisplayVertices(){
    if(curstate == SURF)return Surf_core.getDisplayVertices();
    else if(curstate == CURV)return Curv_core.getDisplayVertices();
    else if(curstate == VOL)return Vol_core.getDisplayVertices();

    return NULL;
}
vector<double>* ESDF_Core::getDisplayVerticesNormal(){
    if(curstate == SURF)return Surf_core.getDisplayVerticesNormal();
    else if(curstate == CURV)return Curv_core.getDisplayVerticesNormal();
    else if(curstate == VOL)return Vol_core.getDisplayVerticesNormal();

    return NULL;
}
vector<uint>* ESDF_Core::getDisplayEdges(){
    if(curstate == SURF)return Surf_core.getDisplayEdges();
    else if(curstate == CURV)return Curv_core.getDisplayEdges();
    else if(curstate == VOL)return Vol_core.getDisplayEdges();

    return NULL;
}
vector<uint>* ESDF_Core::getDisplayFaces(){
    if(curstate == SURF)return Surf_core.getDisplayFaces();
    else if(curstate == CURV)return Curv_core.getDisplayFaces();
    else if(curstate == VOL)return Vol_core.getDisplayFaces();

    return NULL;
}
vector<unsigned char>* ESDF_Core::getDisplayColor(){
    if(curstate == SURF)return Surf_core.getDisplayColor();
    else if(curstate == CURV)return Curv_core.getDisplayColor();
    else if(curstate == VOL)return Vol_core.getDisplayColor();

    return NULL;
}

/****************************************************************************/
/****************************************************************************/

using namespace Eigen;


class mytimer{
public:
    struct timeval tpstart,tpend;
    float timeuse;
    void start(){gettimeofday(&tpstart,NULL);}
    void end(){
        gettimeofday(&tpend,NULL);
               timeuse=(1000000*(tpend.tv_sec-tpstart.tv_sec) + tpend.tv_usec-tpstart.tv_usec)/1000.0;
               cout<<"time: "<<timeuse<<" ms"<<endl;
    }
    mytimer(){}
};

void positiveDefiniteLeastEigenVector(SparseMatrix<double>& A,SparseMatrix<double>& M,vector<double>&outx,const int maxIter = 30){

    int msize = M.cols();
    VectorXd x;
    x.setRandom(msize,1);

    ConjugateGradient<SparseMatrix<double> > solver;
    //CholmodSupernodalLLT<SparseMatrix<double> > solver;
    //SimplicialLDLT<SparseMatrix<double> > solver;
    //BiCGSTAB<SparseMatrix<double>, IncompleteLUT<double>> solver;
    //SparseLU<SparseMatrix<double> > solver;
    //SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > solver;

    solver.compute(A);
    cout<<"MaxIter: "<<maxIter<<endl;
    for(int i=0;i<maxIter;++i){
        cout<<"Iter: "<<i<<endl;
        x = solver.solve(M*x);
        x /= (M*x).dot(x);
        //x/=x.norm();
    }

    outx.resize(x.rows());
    for(int i=0;i<x.rows();++i)outx[i] = x(i);
    cout<<"energy: "<<(A*x).dot(x)<<endl;

}


void EigenInitialization(    const vector<double>&vertices_rfield,const vector<double>&vertices_rofield,
                             const vector<double>&edge_weight,const vector<double>&vertices_Vweight,const vector<double>&edge_Vweight,
                             const vector<uint>&edges2vertices,const vector<uint>&vertices2edges,
                             const vector<uint>&vertices2vertices, const vector< unsigned long >&vertices2edges_accumulate,
                             vector<double>&outx,
                             const int maxIter = 30
                             ){

    cout<<"EigenInitialization Begin"<<endl;
    auto p_vrd = vertices_rfield.data(); auto p_vrod = vertices_rofield.data();
    auto vrvec_begin = [p_vrd](int vind){return p_vrd+vind*3;};
    auto vrovec_begin = [p_vrod](int vind){return p_vrod+vind*3;};

    auto p_vvd = vertices2vertices.data(); auto p_ved = vertices2edges.data();
    auto p_veacd = vertices2edges_accumulate.data();
    //auto p_vsd = vertices2simplex.data();
    //auto p_vsacd = vertices2simplex_accumulate.data();
    auto vv_begin = [p_vvd,p_veacd](uint vind){return p_vvd + p_veacd[vind];};
    auto vv_end = [p_vvd,p_veacd](uint vind){return p_vvd + p_veacd[vind+1];};
    auto ve_begin = [p_ved,p_veacd](uint vind){return p_ved + p_veacd[vind];};
    auto ve_end = [p_ved,p_veacd](uint vind){return p_ved + p_veacd[vind+1];};

    //    auto vs_begin = [p_vsd,p_vsacd](uint vind){return p_vsd + p_vsacd[vind];};
    //    auto vs_end = [p_vsd,p_vsacd](uint vind){return p_vsd + p_vsacd[vind+1];};

    auto p_evd = edges2vertices.data();
    auto ev_begin = [p_evd](uint eind){return p_evd+2*eind;};


    int n_vertices = vertices2edges_accumulate.size()-1;
    int n_edges = edge_weight.size();

    SparseMatrix<double> H(3*n_edges,2*n_vertices);
    SparseMatrix<double> A(2*n_vertices,2*n_vertices);
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;




    tripletList.reserve(n_edges*8+n_vertices*4);
    int offset = n_vertices;
    for(int i=0;i<n_vertices;++i){
        double Aii=0,Aij = 0;
        auto p_ve = ve_begin(i);
        for(auto p_vv = vv_begin(i);p_vv != vv_end(i);++p_vv,++p_ve){
            Aii+=edge_weight[*p_ve];

            Aij = dot(vrvec_begin(i),vrvec_begin(*p_vv))*edge_weight[*p_ve];
            tripletList.push_back(T(i,*p_vv,-Aij));


            Aij = dot(vrovec_begin(i),vrvec_begin(*p_vv))*edge_weight[*p_ve];
            tripletList.push_back(T(i+offset,*p_vv,-Aij));


            Aij = dot(vrvec_begin(i),vrovec_begin(*p_vv))*edge_weight[*p_ve];
            tripletList.push_back(T(i,(*p_vv)+offset,-Aij));


            Aij = dot(vrovec_begin(i),vrovec_begin(*p_vv))*edge_weight[*p_ve];
            tripletList.push_back(T(i+offset,(*p_vv)+offset,-Aij));


        }
        tripletList.push_back(T(i,i,Aii));
        tripletList.push_back(T(i+offset,i+offset,Aii));



    }
    A.setFromTriplets(tripletList.begin(), tripletList.end());


    SparseMatrix<double> M(2*n_vertices,2*n_vertices);
    tripletList.clear();
    tripletList.reserve(2*n_vertices);
    if(0){
        for(int i=0;i<n_vertices;++i){
            auto varea = 1;
            //auto varea = v_area(i);
            tripletList.push_back(T(i*2,i*2,varea));
            tripletList.push_back(T(i*2+1,i*2+1,varea));

        }
    }else{
        for(int i=0;i<n_vertices;++i){
            double Mij = 0;
            auto p_ve = ve_begin(i);
            for(auto p_vv = vv_begin(i);p_vv != vv_end(i);++p_vv,++p_ve){

                Mij = dot(vrvec_begin(i),vrvec_begin(*p_vv))*edge_Vweight[*p_ve];
                tripletList.push_back(T(i,*p_vv,Mij));

                Mij = dot(vrovec_begin(i),vrvec_begin(*p_vv))*edge_Vweight[*p_ve];
                tripletList.push_back(T(i+offset,*p_vv,Mij));


                Mij = dot(vrvec_begin(i),vrovec_begin(*p_vv))*edge_Vweight[*p_ve];
                tripletList.push_back(T(i,(*p_vv)+offset,Mij));


                Mij = dot(vrovec_begin(i),vrovec_begin(*p_vv))*edge_Vweight[*p_ve];
                tripletList.push_back(T(i+offset,(*p_vv)+offset,Mij));


            }
            tripletList.push_back(T(i,i,vertices_Vweight[i]*2));
            tripletList.push_back(T(i+offset,i+offset,vertices_Vweight[i]*2));

        }
    }

    M.setFromTriplets(tripletList.begin(), tripletList.end());


    mytimer mmtt;
    mmtt.start();

    positiveDefiniteLeastEigenVector(A,M,outx,maxIter);

    mmtt.end();



    cout<<"EigenInitialization end"<<endl;
    //void abscccc(int in){cout<<in<<endl;}

}


/****************************************************************************/
/****************************************************************************/

double GaussSeidelIteration(const vector<double>&vertices_normal,
                            const vector<double>&edge_weight,
                            const vector<uint>&edges2vertices,const vector<uint>&vertices2edges,
                            const vector<uint>&vertices2vertices, const vector< unsigned long >&vertices2edges_accumulate,
                            vector<double>&vertices_field,
                            const int maxiter = 5000
        ){



    //for(auto a:edge_weight)cout<<a<<' ';cout<<endl;
    auto p_vvd = vertices2vertices.data(); auto p_ved = vertices2edges.data();
    auto p_veacd = vertices2edges_accumulate.data();

    auto vv_begin = [p_vvd,p_veacd](uint vind){return p_vvd + p_veacd[vind];};
    auto vv_end = [p_vvd,p_veacd](uint vind){return p_vvd + p_veacd[vind+1];};
    auto ve_begin = [p_ved,p_veacd](uint vind){return p_ved + p_veacd[vind];};
    auto ve_end = [p_ved,p_veacd](uint vind){return p_ved + p_veacd[vind+1];};
    auto vv_num = [p_veacd](uint vind){return p_veacd[vind+1]- p_veacd[vind];};

    auto p_evd = edges2vertices.data();
    auto ev_begin = [p_evd](uint eind){return p_evd+2*eind;};

    auto p_vnd = vertices_normal.data();
    auto vnor_begin = [p_vnd](uint vind){return p_vnd+3*vind;};

    int n_vertices = vertices2edges_accumulate.size()-1;
    int n_edges = edge_weight.size();

    auto EdgesDifferetialEnergy=[&edge_weight,ev_begin,n_edges](double *p_field){
        double en = 0;

        for(uint i=0;i<n_edges;++i){
            auto p_ev = ev_begin(i);
            en+=edge_weight[i]*vecSquareDist(p_field+(p_ev[0])*3,p_field+3*(p_ev[1]));
        }
        return en;
    };

    //int maxiter = 10000;
    int iter = 0;


    double alpha = 0.5;

    vector<double> container1  = vertices_field;
    vector<double> container2  = vertices_field;
    double *p_pre = container1.data(),*p_cur = container2.data();
    double en_pre =  EdgesDifferetialEnergy(p_pre),en;

    vector<double>sumweight(n_vertices,0.);
    vector<double>veweight(vertices2edges.size(),1);
    auto p_veweightbegin = veweight.data();
    auto vew_begin = [p_veweightbegin,p_veacd](uint vind){
        return p_veweightbegin+p_veacd[vind];
    };
    for(int i=0;i<n_vertices;++i){
        //auto p_vv = vv_begin(i);
        auto p_ve = ve_begin(i);
        auto vvnum = vv_num(i);
        auto p_vew = vew_begin(i);
        for(int j=0;j<vvnum;++j){
            p_vew[j] = edge_weight[p_ve[j]];
            sumweight[i]+=p_vew[j];
        }
        for(int j=0;j<vvnum;++j){
            p_vew[j] /=  sumweight[i];

        }

    }


    //for(auto a:veweight)cout<<a<<' ';cout<<endl;
    cout<<"GaussSeidelIteration Begin: "<<en_pre<<endl;

    mytimer mmtt;
    mmtt.start();
    while(true){
        memset(p_cur,0,sizeof(double)*n_vertices*3);
        for(int i=0;i<n_vertices;++i){
            auto p_vv = vv_begin(i);
            auto vvnum = vv_num(i);

            auto p_curv = p_cur+i*3;
            auto p_vew = vew_begin(i);


            for(int j=0;j<vvnum;++j){
                weightedAddVec(p_vew[j],p_pre+p_vv[j]*3,p_curv);
            }
            projectVectorNor(p_curv,vnor_begin(i),p_curv);

            weightedAddVec(alpha,1.-alpha,p_pre+i*3,p_curv,p_curv);
            projectVectorNor(p_curv,vnor_begin(i),p_curv);

        }
        en = EdgesDifferetialEnergy(p_cur);
        //cout<<en<<' '<<en_pre<<endl;
        if(fabs(en-en_pre)/en_pre < 1e-6 || iter>maxiter)break;

        en_pre = en;
        swap(p_pre,p_cur);
        ++iter;
        //cout<<"iter: "<<iter<<' '<<en_pre<<endl;

    }

    memcpy(vertices_field.data(),p_pre,n_vertices*3*sizeof(double));

    mmtt.end();
    cout<<"GaussSeidelIteration Finished!"<<endl;
    cout<<"iter: "<<iter<<' '<<en_pre<<endl;

    return en_pre;

}




}//n_rf
