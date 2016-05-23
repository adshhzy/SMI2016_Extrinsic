#include "geo_vol.h"
#include "readers.h"
namespace n_rf {



double vfun1(double x,double y,double z,double *grad){

    //x^2+2y^2+3z^2

    grad[0]=2*x;
    grad[1]=4*y;
    grad[2]=6*z;


    normalize(grad);

    return x*x+2*y*y+3*z*z;


}

double vfun2(double x,double y,double z,double *grad){

    //x^2+2y^2+3z^2

    grad[0]=2*x;
    grad[1]=2*y;
    grad[2]=2*z;


    normalize(grad);

    return x*x+y*y+z*z;


}

double vfun3(double x,double y,double z,double *grad){

    //x^2+2y^2+3z^2

    grad[0]=2*x;
    grad[1]=4*y;
    grad[2]=3*z*z;


    normalize(grad);

    return x*x+2*y*y+z*z*z;


}

double isotorus(double x,double y,double z,double *grad){

    //x^2+2y^2+3z^2

    double circleradiu = 0.5;


    double nnorm = sqrt(x*x+y*y);

    double xx = circleradiu*x/nnorm;double yy = circleradiu*y/nnorm;



    grad[0]=x - xx;
    grad[1]=y - yy;
    grad[2]=z;
    double rev = len(grad);


    normalize(grad);

    return rev;


}



void Volume::IsoToyGeneration(int w, int h, int d,double wrange,double hrange,double drange){

    width = w;height = h;depth = d;

    //double wrange = 1.4,hrange = 1.4,drange = 0.4;
    xs = wrange/2.;ys = hrange/2.;zs = drange/2.;

    xsa_step = wrange/(width);
    ysa_step = hrange/(height);
    zsa_step = drange/(depth);

    vfun = isotorus;

    volume.resize(depth,vector<double >(width*height) );
    vertices.resize(width*height*depth*3);
    vertices_normal.resize(width*height*depth*3);
    double *p_v = vertices.data();
    double *p_vn = vertices_normal.data();
    for(int i=0;i<depth;++i){
        double *p_vol = volume[i].data();
        double z = -zs+i*zsa_step;
        for(int j=0;j<height;++j){
            double y = -ys+j*ysa_step;
            double x = -xs;
            for(int k=0;k<width;++k){
                (*p_vol) = vfun(x,y,z,p_vn);
                p_v[0]=x;p_v[1]=y;p_v[2]=z;
                p_v+=3;++p_vol;x+=xsa_step; p_vn+=3;
            }
        }
    }

    n_vertices = width*height*depth;
    display_vcolor.resize(n_vertices*4,155);

    int subdivct = 6,subdivcv = subdivct*4;
    n_tetr = (width-1)*(height-1)*(depth-1)*subdivct;
    tetrahedron2vertices.resize(n_tetr*4);


    int wh = width*height;
    //    int bias[20];

    //    bias[0] = 0;bias[1] = 0+1;bias[2] = 0+1+width;bias[3] = 0+1+wh;

    //    bias[4] = 0;bias[5] = 0+wh;bias[6] = 0+wh+1;bias[7] = 0+wh+width;

    //    bias[8] = 0;bias[9] = 0+width+1;bias[10] = 0+width;bias[11] = 0+wh+width;

    //    bias[12] = 0+wh+width+1;bias[13] = 0+width+1;bias[14] = 0+wh+1;bias[15] =  0+wh+width;

    //    bias[16] = 0;bias[17] = 0+width+1;bias[18] = 0+1+wh;bias[19] = 0+width +wh;

    //1248,1258,2348,2378,2568,2678
    int cubeInd[8];

    cubeInd[0] = 0;cubeInd[1] = 1;cubeInd[2] = width+1;cubeInd[3] = width;
    cubeInd[4] = wh;cubeInd[5] = wh+1;cubeInd[6] = wh+width+1;cubeInd[7] = wh+width;

    int bias[24] = {1,2,4,8,1,2,5,8,2,3,4,8,2,3,7,8,2,5,6,8,2,6,7,8};
    for(int i=0;i<24;++i)bias[i] = cubeInd[bias[i]-1];





    auto p_tv = tetrahedron2vertices.data();
    for(int k=0;k<depth-1;++k){
        int indde = k*(width-1)*(height-1)*subdivcv;
        int vindk = k*wh;
        for(int i = 0;i<height-1;++i){
            int indhe = indde+i*(width-1)*subdivcv;
            int vindh = vindk+i*width;
            for(int j=0;j<width-1;++j){
                int ind = indhe+j*subdivcv;
                int vind = vindh+j;
                for(int l=0;l<subdivcv;l++)p_tv[ind+l] = vind+bias[l];
            }
        }
    }


    vector<double>mmmax(depth);
    for(int i=0;i<depth;++i)mmmax[i] = *(max_element(volume[i].begin(),volume[i].end()));
    thres_upper = *(max_element(mmmax.begin(),mmmax.end()));
    for(int i=0;i<depth;++i)mmmax[i] = *(min_element(volume[i].begin(),volume[i].end()));
    thres_lower = *(min_element(mmmax.begin(),mmmax.end()));

    cout<<"thres_lower: "<<thres_lower<<' '<<"thres_upper: "<<thres_upper<<endl;

    //ContourVolume(100);

}

bool Volume::ReadInterface(string filename){
    bool a = readVolfFile(filename,vertices,tetrahedron2vertices,vertices_normal,vertices_field);
    if(a)setparameters();

    return a;
}
bool Volume::SaveInterface(string filename){


    if(n_vertices==0||n_tetr==0)return false;
    return writeVolfFile(filename,vertices,tetrahedron2vertices,vertices_normal,vertices_field);

}
bool Volume::outputDisplay(string filename){

    return writePLYFile(filename,display_vertices,display_faces,display_normal,display_vcolor);

}
}//n_rf
