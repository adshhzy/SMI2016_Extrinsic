#include "geo_curv.h"

#include<iostream>
#include<fstream>


namespace n_rf {

bool Curve::isload = false;
Mesh Curve::sphere;
Mesh Curve::cylinder;
Mesh Curve::cone;
Curve::Curve():n_vertices(0),n_edges(0),isSettangent(false),isBuildFrenet(false),isbuild(false)
{

    load();
    isCameraMode = false;

}
void Curve::reset(){
    vertices.clear(); edges.clear(); tangent.clear();

    display_edges.clear();display_vertices.clear();v_color.clear();display_faces.clear();

    n_vertices = n_edges = 0;


    vertices_field.clear();
    vertices_rfield.clear();

    vertices_ofield.clear();
    vertices_rofield.clear();

    vertices_f_field.clear();

    isbuild = false;
    isSettangent = false;
    isBuildFrenet = false;

}
void Curve::setparameters(){
    n_vertices = vertices.size()/3;
    n_edges = edges.size()/2;
    if(n_vertices == n_edges)isloop=true;
    else isloop = false;
}


bool Curve::ImportCurve(vector<double>& vertices_in, vector<uint> &edges_in){

    reset();

    vertices = vertices_in;
    edges = edges_in;

    setparameters();

    BuildEdges(false);

    return true;

}


bool Curve::ReadCurve(string filename){
    cout<<"Reading: "<<filename<<endl;
    ifstream fin(filename.data());
    if(fin.fail())cout<<"error"<<endl;
    reset();
    fin>>n_vertices;fin>>n_edges;
    if(n_vertices==n_edges)isloop = true;
    else if(n_vertices==n_edges+1)isloop = false;
    else {cout<<"Not a curve!"<<endl;return false;}

    cout<<"read: "<<n_vertices<<' '<<n_edges<<endl;
    vertices.resize(n_vertices*3);
    for(int i =0;i<n_vertices*3;++i)fin>>vertices[i];

    edges.resize(n_edges*2);
    bool hasedges = true;
    for(int i =0;i<n_edges*2;++i){
        fin>>edges[i];
        if(fin.eof()){hasedges = false;break;}
    }
    fin.close();
    if(hasedges)SortLoopEdges();
    BuildEdges();
    tangent.clear();

    setparameters();

    //EstimateFrenetFrame();
    //BuildEdges(false);

    //BuildDisplay(false);

    return true;
}
bool Curve::ReadCurve(ifstream &fin, int p0, int p1, int p2){
    n_vertices = p0;
    if(p1==0)isloop = false;
    else if(p1==1)isloop = true;
    else {cout<<"Not a curve!"<<endl;return false;}

    vertices.resize(n_vertices*3);
    for(int i =0;i<n_vertices*3;++i)fin>>vertices[i];


    edges.clear();
    for(int i =0;i<n_vertices-1;i++){
        edges.push_back(i);
        edges.push_back(i+1);
    }
    if(isloop){edges.push_back(n_vertices-1);edges.push_back(0);}


    tangent.clear();
    //BuildDisplay();
    setparameters();
    EstimateFrenetFrame();
    BuildEdges(false);



    return true;


}


void Curve::BuildEdges(bool isseq){

    if(isseq){
        edges.clear();
        for(int i =0;i<n_vertices-1;i++){
            edges.push_back(i);
            edges.push_back(i+1);
        }
        if(isloop){edges.push_back(n_vertices-1);edges.push_back(0);n_edges = n_vertices;}
        else n_edges = n_vertices-1;
    }

    edge_len.resize(n_edges);
    inverse_edge_len.resize(n_edges);
    double evec[3];
    for(uint i = 0;i<n_edges;++i){
        auto p_ev = ev_begin(i);
        vectorize(p_ev[0],p_ev[1],evec);
        edge_len[i]=normVec(evec);
        inverse_edge_len[i] = 1./edge_len[i];

    }

    vertices2edges.clear();
    vertices2edges_accumulate.clear();
    vertices2edges_accumulate.resize(n_vertices,0);

    for (uint i = 0; i < n_edges; i++){
        auto p_ev = ev_begin(i);
        for (uint j = 0; j < 2; ++j){
            ++vertices2edges_accumulate[p_ev[j]];

        }
    }


    //for(a:vertices2edges_accumulate)cout<<a<<' ';cout<<endl;
    for (uint i = 1; i < n_vertices; i++){
        vertices2edges_accumulate[i]+=vertices2edges_accumulate[i-1];
    }
    vertices2edges_accumulate.insert(vertices2edges_accumulate.begin(),0);
    vertices2edges.resize(vertices2edges_accumulate[n_vertices]);
    vector<uchar>cur_num(n_vertices,0);

    for (uint i = 0; i < n_edges; i++){
        auto p_ev = ev_begin(i);
        for (uint j = 0; j < 2; ++j){
            auto v = p_ev[j];
            vertices2edges[vertices2edges_accumulate[v]+cur_num[v]] = i;
            cur_num[v]++;
        }
    }

    vertices2edgesPN.resize(vertices2edges.size(),false);
    for(uint i=0;i<n_vertices;++i){
        auto p_vepn = vepn_begin(i);
        for(auto p_ve = ve_begin(i);p_ve!=ve_end(i);++p_ve,++p_vepn){
            if((*ev_begin(*p_ve))==i)(*p_vepn)=true;
            else (*p_vepn)=false;

        }
    }


    vertices2vertices.resize(vertices2edges.size());
    for(int i=0;i<n_vertices;++i){
        auto p_vepn = vepn_begin(i);
        auto p_vv = vv_begin(i);
        for(auto p_ve = ve_begin(i);p_ve!=ve_end(i);++p_ve,++p_vepn,++p_vv)
            if(*p_vepn)(*p_vv) = ev_begin(*p_ve)[1];
            else (*p_vv) = ev_begin(*p_ve)[0];
    }

    vertices2vertices_inverse.resize(vertices2edges.size());
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

    vertices_len.resize(n_vertices);
    for(int i=0;i<n_vertices;++i){
        auto p_ve = ve_begin(i);
        auto venum = ve_num(i);
        if(venum==2)vertices_len[i] = (elen_begin(p_ve[0])+elen_begin(p_ve[1]))/2.;
        else if(venum==1)vertices_len[i] = elen_begin(p_ve[0]);
    }



}

void Curve::EliminateDuplication(double thres){
    setparameters();
    vector<double>newvertices;
    vector<bool>iskeep(n_vertices,true);


    for(int i=0;i<n_vertices-1;++i){
        if(VerticesDistance(i,i+1)<thres)iskeep[i+1]=false;
    }
    if(isloop)if(VerticesDistance(0,n_vertices-1)<thres)iskeep[n_vertices-1]=false;
    int count = 0;
    for(int i=0;i<n_vertices-1;++i){
        if(iskeep[i]){
            count++;
            auto p_o = v_begin(i);
            newvertices.push_back(p_o[0]);
            newvertices.push_back(p_o[1]);
            newvertices.push_back(p_o[2]);
        }
    }
    vertices=newvertices;
    cout<<"Duplication "<<n_vertices<<' '<<count<<endl;
    setparameters();
    BuildEdges();



}


bool Curve::isBuildDisplay(){return isbuild;}





bool Curve::SortLoopEdges(){

    setparameters();
    auto old_edges = edges;
    auto old_vertices = vertices;
    auto old_vf_field = vertices_f_field;
    bool isRerangeffield = (vertices_f_field.size()==vertices.size());
    vector< vector<int> >neighborhood(n_vertices);
    for (int i = 0; i < edges.size(); i += 2){
        neighborhood[edges[i]].push_back(edges[i + 1]);
        neighborhood[edges[i + 1]].push_back(edges[i]);
    }
    if(isloop)for (auto &v_v : neighborhood)if (v_v.size() != 2){ cout << "Not a closed polygon!" << endl; return false; }
    //for (auto v_v : neighborhood)if (v_v[0] == v_v[1]){ cout << "Not a closed polygon!" << endl; }
    //for (auto v_v : neighborhood){ cout << v_v[0]<<"  "<< v_v[1] << endl; }
    vector<int>orderlist(n_vertices, -1);

    //map<pair<uint, uint>, uint>edges_index;
    vector< vector<int> >edges_index2(n_vertices, vector<int>(n_vertices));

    int seed_index = 0, pre_seed_index = -1, order = 0;
    if(!isloop){
        for(int i=0;i<neighborhood.size();++i)if (neighborhood[i].size() == 1){seed_index=i;break;}
    }
    for (int i = 0; i < n_vertices; ++i){
        orderlist[seed_index] = order;
        auto& v_v = neighborhood[seed_index];
        for (auto next_seed : v_v)if (next_seed != pre_seed_index){
            //edges_index[make_pair(seed_index, next_seed)] = order;
            //edges_index[make_pair(next_seed, seed_index)] = order;
            edges_index2[next_seed][seed_index] = order;
            edges_index2[seed_index][next_seed] = order;
            pre_seed_index = seed_index;
            seed_index = next_seed;
            break;
        }
        ++order;
    }
    //cout << order << endl;
    //for (auto new_order : orderlist)cout << new_order << endl;
    for (auto new_order : orderlist)if (new_order == -1){ cout << "Not a closed polygon!" << endl; return false; }
    for (int i = 0; i < n_vertices; i++){
        for (int j = 0; j < 3; ++j){
            vertices[3 * orderlist[i] + j] = old_vertices[3 * i + j];
            if(isRerangeffield)vertices_f_field[3 * orderlist[i] + j] = old_vf_field[3 * i + j];
        }
    }
    for (int i = 0; i < edges.size(); i += 2){
        //auto e_key = make_pair(old_edges[i], old_edges[i + 1]);
        //auto new_index = 2 * edges_index[e_key];
        auto new_index = 2 * edges_index2[old_edges[i]][old_edges[i + 1]];
        auto p1 = orderlist[old_edges[i]];
        auto p2 = orderlist[old_edges[i + 1]];
        if (p1 < p2){ edges[new_index] = p1; edges[new_index + 1] = p2; }
        else { edges[new_index] = p2; edges[new_index + 1] = p1; }
        if (edges[new_index] == 0 && edges[new_index + 1] == n_vertices - 1)swap(edges[new_index], edges[new_index + 1]);
    }

    BuildEdges(false);

    return true;

}

bool Curve::EstimateFirstDerivativeO5(vector<double>& invect, vector<double>& derivative, int interval,bool isloop){
    if(invect.size()==0 || invect.size()%interval!=0)return false;

    derivative.resize(invect.size());

    auto singleinterval = [this](vector<double>&invec, vector<double>&outvec,bool isloop){
        int vecN = invec.size();
        outvec.resize(vecN);
        if(isloop){
            outvec[0] = 1/12.0*invec[vecN-2] - 8/12.0*invec[vecN-1] + 8/12.0*invec[1] - 1/12.0*invec[2];
            outvec[1] = 1/12.0*invec[vecN-1] - 8/12.0*invec[0] + 8/12.0*invec[2] - 1/12.0*invec[3];
            outvec[vecN-1] = 1/12.0*invec[vecN-3] - 8/12.0*invec[vecN-2] + 8/12.0*invec[0] - 1/12.0*invec[1];
            outvec[vecN-2] = 1/12.0*invec[vecN-4] - 8/12.0*invec[vecN-3] + 8/12.0*invec[vecN-1] - 1/12.0*invec[0];
        }else{
            outvec[0] = -25/12.0*invec[0] +48/12.0*invec[1]-36/12.0*invec[2]+16/12.0*invec[3]-3/12.0*invec[4];
            outvec[1] = (-3/12.0*invec[0] -10/12.0*invec[1]+18/12.0*invec[2]-6/12.0*invec[3]+1/12.0*invec[4]);
            outvec[vecN-1] = -(-25/12.0*invec[vecN-1] +48/12.0*invec[vecN-2]-36/12.0*invec[vecN-3]+16/12.0*invec[vecN-4]-3/12.0*invec[vecN-5]);
            outvec[vecN-2] = -(-3/12.0*invec[vecN-1] -10/12.0*invec[vecN-2]+18/12.0*invec[vecN-3]-6/12.0*invec[vecN-4]+1/12.0*invec[vecN-5]);
        }

        for(int i=2;i<vecN-2;i++){
            outvec[i] = 1/12.0*invec[i-2] - 8/12.0*invec[i-1] + 8/12.0*invec[i+1] - 1/12.0*invec[i+2];
        }

    };
    int vecNN = invect.size()/interval;
    vector<double>invec(vecNN);
    vector<double>outvec(vecNN);
    for(int i=0;i<interval;++i){
        for(int j=0;j<vecNN;++j)invec[j] = invect[i+j*interval];
        singleinterval(invec,outvec,isloop);
        for(int j=0;j<vecNN;++j)derivative[i+j*interval] = outvec[j];
    }

    return true;
}
void Curve::BuildTangent(){

    if(n_vertices==0||n_edges==0)return;
    tangent.resize(n_vertices*3);

    EstimateFirstDerivativeO5(vertices,tangent,3,isloop);
    linespeed.resize(n_vertices);
    for(int i=0;i<n_vertices;++i)linespeed[i] = len(t_begin(i));
    for(int i=0;i<n_vertices;++i)normalize(t_begin(i));
    for(int i=0;i<tangent.size();++i)if(tangent[i]!=tangent[i])exit(121);
    isSettangent = true;
}
void Curve::EstimateFrenetFrame(bool isReSetTangent){

    if(n_vertices==0||n_edges==0)return;

    //BuildTangent();

    if(isReSetTangent || tangent.size()!=n_vertices*3){
        BuildTangent();
    }


    Frenetnormal.resize(n_vertices*3);
    EstimateFirstDerivativeO5(tangent,Frenetnormal,3,isloop);

    curvature.resize(n_vertices);
    for(int i=0;i<n_vertices;++i)curvature[i] = len(nor_begin(i))/linespeed[i];

    normCurvature.resize(n_vertices);
    double cmaxC = 3.5/(*max_element(curvature.begin(),curvature.end()));
    for(int i=0;i<n_vertices;++i)normCurvature[i] =min(1., curvature[i] *cmaxC);
    CurveScalarFunctionSmoothing(normCurvature,20);

    //for(int i=0;i<n_vertices;++i)cout<<curvature[i]<<' ';
    //cout<<endl;

    for(int i=0;i<n_vertices;++i)normalize(nor_begin(i));
    for(int i=0;i<tangent.size();++i)if(Frenetnormal[i]!=Frenetnormal[i])exit(121);


    Frenetbinormal.resize(n_vertices*3);
    for(int i=0;i<n_vertices;++i)cross(t_begin(i),nor_begin(i),binor_begin(i));
    for(int i=0;i<n_vertices;++i)normalize(binor_begin(i));


    vector<double>dbinormal(n_vertices*3);
    torsion.resize(n_vertices);
    EstimateFirstDerivativeO5(Frenetbinormal,dbinormal,3,isloop);

    for(int i=0;i<n_vertices;++i)torsion[i] = -dot(nor_begin(i),&(dbinormal[i*3]));
    for(int i=0;i<n_vertices;++i)torsion[i] /=linespeed[i];

    if(1)for(int i=0;i<n_vertices;++i)inversevec(binor_begin(i),binor_begin(i));
    //for(int i=0;i<n_vertices;++i)cout<<torsion[i]<<' ';
    //cout<<endl;





    isBuildFrenet = true;
    isSettangent = true;
    //cout<<"EstimateFrenetFrame"<<endl;
    //for(auto j:Frenetbinormal)cout<<j<<' ';cout<<endl;


}
void Curve::RescaleUniform(){
    double xmin = 999, ymin=999,zmin=999,xmax=-999,ymax=-999,zmax=-999;
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
    //double largestdis = max(max(xmax-xmin,ymax-ymin),zmax-zmin)/2;
    vector<double>dis(3);
    dis[0] = (xmax-xmin)/2;dis[1] = (ymax-ymin)/2;dis[2] = (zmax-zmin)/2;
    //cout<<centers[0]<<"  "<<centers[1]<<"  "<<centers[2]<<"   "<<largestdis<<endl;

    for (int j = 0; j<n_vertices*3; j+=3)
        for(int i = 0; i<3; i++)
            vertices[j+i]=(vertices[j+i]-centers[i])/dis[i];

}


void Curve::CurveScalarFunctionSmoothing(vector<double>&scalarF,int maxiter){


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

void Biarc(double u, double *p_v, double *p_t, double *p_b){

    if(u<PI){
        p_v[0] = sin(u);
        p_v[1] = cos(u);
        p_v[2] = 0;
    }else{
        p_v[0] = sin(u);
        p_v[1] = -1;
        p_v[2] = cos(u)+1;
    }

}
void MonkeyCurve3(double u, double *p_v, double *p_t, double *p_b){
    p_v[0] = cos(u);
    p_v[1] = sin(u);
    p_v[2] = pow(cos(u),4) - 4 *pow(cos(u),2)*pow(sin(u),2);

    p_t[0] = -sin(u);
    p_t[1] = cos(u);
    p_t[2] = -4.0*sin(u)*pow(cos(u),3) + 4*(2*pow(sin(u),3)*cos(u) - 2*pow(cos(u),3)*sin(u));
}
void Curve::BuildToyData(int ind, int numV, double start, double end){

    reset();

    auto ptr = Biarc;
    switch (ind) {
    case 2: ptr = MonkeyCurve3;isloop = true;break;
    case 1: ptr = Biarc;isloop = false;break;
    default:
        break;
    }

    double L = end-start;
    if(L<0){L=-L;start = end;}
    int numofpoint = numV;
    double step;
    cout<<L<<' '<<numofpoint<<endl;
    if(isloop){
        if(L>3.1415926535898*2){
            step = 3.1415926535898 * 2.0/numofpoint;
            start = 0.0;
        }
        else {isloop = false;step = L/(numofpoint-1);}
    }else{
        step = L/(numofpoint-1);
    }
    vertices.resize(numofpoint*3);
    tangent.resize(numofpoint*3,1);
    Frenetbinormal.resize(numofpoint*3,1);
    for(int i =0;i<numofpoint;i++){
        auto p_v = v_begin(i);
        auto p_t = t_begin(i);
        auto p_b = binor_begin(i);
        ptr(start+i*step,p_v,p_t,p_b);
        normalize(p_t);
        normalize(p_b);
    }
    //edges.clear();
    for(int i =0;i<numofpoint-1;i++){
        edges.push_back(i);
        edges.push_back(i+1);
    }
    if(isloop){edges.push_back(numofpoint-1);edges.push_back(0);}

    isbuild = false;
    isSettangent = false;
    setparameters();
    BuildEdges(false);
    EstimateFrenetFrame();


}

/************************************************************************/
/************************************************************************/
void Curve::BuildOthogonalField(){
    if(n_vertices==0)return;
    if(n_edges==0)return;


    if(vertices_field.size()!=0){
        vertices_ofield.resize(vertices_field.size());
        for(int i =0;i<n_vertices;++i)cross(vfvec_begin(i),vvec_begin(i),vovec_begin(i));
    }

    if(vertices_rfield.size()!=0){
        vertices_rofield.resize(vertices_rfield.size());
        for(int i =0;i<n_vertices;++i)cross(vfvec_begin(i),vrvec_begin(i),vrovec_begin(i));
    }
}
double randf();
void Curve::InitializeField(bool isProj,bool isrand){

    if(n_vertices==0)return;
    if(n_edges==0)return;

    if(vertices_f_field.size()!=n_vertices*3){
        if(isCameraMode){
            cout<<"CameraMode: calculate LookAt Vector."<<endl;
            cout<<"focusing: "<<focusing[0]<<' '<<focusing[1]<<' '<<focusing[2]<<endl;
            vertices_f_field.resize(n_vertices*3);
            for(int i=0;i<n_vertices;++i){
                minusVec(focusing,v_begin(i),vfvec_begin(i));
                normalize(vfvec_begin(i));
            }
        }else {
            cout<<"Curve mode: using tangent field as default f."<<endl;
            vertices_f_field = tangent;
        }

    }
    double iiv[3] = {1/sqrt(3),1/sqrt(3),0};

    vertices_rfield.resize(n_vertices*3);
    double riiv[3];
    vertices_rfield.resize(n_vertices*3);

    for(uint i = 0;i<n_vertices;++i){
        auto p_vr = vrvec_begin(i);
        if(isrand)for(int i =0;i<3;++i)iiv[i]= randf();
        projectVectorNor(iiv,vfvec_begin(i),p_vr);
        while(p_vr[0]!=p_vr[0] || p_vr[1]!=p_vr[1] || p_vr[2]!=p_vr[2]){
            for(int j =0;j<3;++j)riiv[j]= randf();
            projectVectorNor(iiv,vfvec_begin(i),p_vr);
        }

    }




    if(isProj || vertices_field.size()==0)
    {
        vertices_field = vertices_rfield;

    }else{
        for(uint i = 0;i<n_vertices;++i){
            projectVectorNor(vvec_begin(i),vfvec_begin(i),vvec_begin(i));
        }
    }



    BuildOthogonalField();

    cout<<"curve field initialize!"<<endl;



}

void Curve::InverseField(){

    for(auto &a:vertices_field)a=-a;
    for(auto &a:vertices_ofield)a=-a;

}

}//n_rf

