#include "geo_vol.h"

#include<list>



namespace n_rf {

void ComputeAngleAxisMatrix(Eigen::AngleAxisd &m_aa,double *oriVec,double *newVec);
void Volume::load(){
    if(!isload){
        sphere.createToy(3);
        sphere.ReScale_uniform(0.005);
        sphere.BuildNeighborTable();
        sphere.ComputeFaceNormal(true);

        cylinder.createToy(2);
        cylinder.BuildNeighborTable();
        cylinder.ComputeFaceNormal(true);
        cylinder.MoveBottom(false,false,true);

        cone.createToy(1);
        cone.ReScale(2.5,2.5,0.5);

        cone.BuildNeighborTable();
        cone.ComputeFaceNormal(true);




        isload = true;
    }
}
void Volume::BuildDisplay(infoVolDisp info){

    if(n_vertices==0)return;
    //if(n_tetr==0)return;

    //cout<<n_tetr<<endl;
    vector<uchar> isShowVertices(n_vertices,1);

    bool ispoint = info.isSurface;
    bool isfield = info.isField;
    bool isVMF = info.isSingularity;
    double widthscale = info.upnormal/100.*0.03;
    double lengthscale = info.length/200.;
    double thicknessscale = info.thickness/10000.;
    unsigned char red[4] = {255,0,0,255};
    unsigned char green[4] = {0,255,0,255};
    unsigned char dark[4] = {0,0,0,255};
    unsigned char blue[4] = {0,0,255,255};
    unsigned char pink[4] = {255,204,204,255};

    display_faces.clear();
    display_vertices.clear();
    display_normal.clear();
    display_edges.clear();
    display_vcolor.clear();


    if(!info.isIso){
        display_faces.clear();
        display_vertices.clear();
        display_normal.clear();
        display_edges.clear();
        display_vcolor.clear();
        double *p_spv = NULL, *p_v = NULL,*p_vn = NULL;
        uint *p_sf = NULL;
        auto n_fs = sphere.n_faces*3;
        auto n_vs = sphere.n_vertices*3;

        if(info.isSlice){

            //isShowVertices.resize(n_vertices);
            int snum = int(double(info.slicenum)/100.*double(depth));
            memset(isShowVertices.data(),0,n_vertices*sizeof(uchar));
            int layernum = width*height;
            memset(isShowVertices.data()+snum*layernum,1,layernum*sizeof(uchar));

        }


        if(ispoint){
            for(int iv = 0;iv<n_vertices;++iv){
                if(!isShowVertices[iv])continue;
                if(sparseShowField.size()!=0)if(!sparseShowField[iv])continue;
                int offset = display_vertices.size()/3;
                p_v = v_begin(iv);
                for(int j = 0;j<sphere.n_vertices;++j){
                    p_spv = sphere.v_begin(j);
                    for(int k=0;k<3;++k)display_vertices.push_back(p_v[k]+p_spv[k]);
                }
                p_sf = sphere.fv_begin(0);
                for(int j = 0;j<n_fs;++j){
                    display_faces.push_back(offset+p_sf[j]);
                }
                p_vn = sphere.vnor_begin(0);
                for(int j = 0;j<n_vs;++j){
                    display_normal.push_back(p_vn[j]);
                }
            }

            display_vcolor.resize(display_vertices.size()/3*4,105);
        }

        uint rect[6] = {0,1,3,0,3,2};

        //double thicknessscale = 0.01;
        auto constructRect = [this,&widthscale,&lengthscale,&thicknessscale,&isShowVertices,rect](double *p_vvec,double *p_ovvec,unsigned char*cc){

            uint offset;
            auto p_vec = p_vvec;
            auto p_ovec = p_ovvec;
            double vecw[3],vecl[3],vect[3],orth[3],kkk[3];
            //box1
            for(int i=0;i<n_vertices;++i){
                if(!isShowVertices[i])continue;
                if(sparseShowField.size()!=0)if(!sparseShowField[i])continue;
                p_vec=p_vvec+i*3;p_ovec=p_ovvec+i*3;
                offset = display_vertices.size()/3;
                auto p_v = v_begin(i);
                auto p_vn = vnor_begin(i);
                cross(p_vec,p_ovec,orth);
                weightedAddVec(thicknessscale,1,orth,p_v,kkk);
                product(widthscale,p_ovec,vecw);
                product(lengthscale,p_vec,vecl);
                for(int j=0;j<3;++j)display_vertices.push_back(kkk[j]+vecw[j]);
                for(int j=0;j<3;++j)display_vertices.push_back(kkk[j]-vecw[j]);
                for(int j=0;j<3;++j)display_vertices.push_back(kkk[j]+vecw[j]+vecl[j]);
                for(int j=0;j<3;++j)display_vertices.push_back(kkk[j]-vecw[j]+vecl[j]);

                for(int j=0;j<6;++j)display_faces.push_back(offset+rect[j]);
                for(int k=0;k<4;++k)for(int j=0;j<3;++j)display_normal.push_back(orth[j]);

                for(int k=0;k<4;++k)for(int j=0;j<4;++j)display_vcolor.push_back(cc[j]);

            }

            //box2
            if(1)for(int i=0;i<n_vertices;++i){
                if(!isShowVertices[i])continue;
                if(sparseShowField.size()!=0)if(!sparseShowField[i])continue;
                p_vec=p_vvec+i*3;p_ovec=p_ovvec+i*3;
                offset = display_vertices.size()/3;
                auto p_v = v_begin(i);
                auto p_vn = vnor_begin(i);
                cross(p_vec,p_ovec,orth);
                weightedAddVec(-thicknessscale,1,orth,p_v,kkk);
                product(widthscale,p_ovec,vecw);
                product(lengthscale,p_vec,vecl);
                for(int j=0;j<3;++j)display_vertices.push_back(kkk[j]+vecw[j]);
                for(int j=0;j<3;++j)display_vertices.push_back(kkk[j]-vecw[j]);
                for(int j=0;j<3;++j)display_vertices.push_back(kkk[j]+vecw[j]+vecl[j]);
                for(int j=0;j<3;++j)display_vertices.push_back(kkk[j]-vecw[j]+vecl[j]);

                for(int j=0;j<6;++j)display_faces.push_back(offset+rect[j]);
                for(int k=0;k<4;++k)for(int j=0;j<3;++j)display_normal.push_back(-orth[j]);

                for(int k=0;k<4;++k)for(int j=0;j<4;++j)display_vcolor.push_back(cc[j]);

            }
            //box3
            if(1)for(int i=0;i<n_vertices;++i){
                if(!isShowVertices[i])continue;
                if(sparseShowField.size()!=0)if(!sparseShowField[i])continue;
                p_vec=p_vvec+i*3;p_ovec=p_ovvec+i*3;
                offset = display_vertices.size()/3;
                auto p_v = v_begin(i);
                cross(p_vec,p_ovec,orth);

                product(thicknessscale,orth,vect);
                product(widthscale,p_ovec,vecw);
                product(lengthscale,p_vec,vecl);

                add(p_v,vecw,kkk);

                for(int j=0;j<3;++j)display_vertices.push_back(kkk[j]+vect[j]);
                for(int j=0;j<3;++j)display_vertices.push_back(kkk[j]-vect[j]);


                for(int j=0;j<3;++j)display_vertices.push_back(kkk[j]+vect[j]+vecl[j]);
                for(int j=0;j<3;++j)display_vertices.push_back(kkk[j]-vect[j]+vecl[j]);

                for(int j=0;j<6;++j)display_faces.push_back(offset+rect[j]);
                for(int k=0;k<4;++k)for(int j=0;j<3;++j)display_normal.push_back(p_ovec[j]);

                for(int k=0;k<4;++k)for(int j=0;j<4;++j)display_vcolor.push_back(cc[j]);

            }
            //box4
            if(1)for(int i=0;i<n_vertices;++i){
                if(!isShowVertices[i])continue;
                if(sparseShowField.size()!=0)if(!sparseShowField[i])continue;
                p_vec=p_vvec+i*3;p_ovec=p_ovvec+i*3;
                offset = display_vertices.size()/3;
                auto p_v = v_begin(i);
                cross(p_vec,p_ovec,orth);

                product(thicknessscale,orth,vect);
                product(widthscale,p_ovec,vecw);
                product(lengthscale,p_vec,vecl);

                minusVec(p_v,vecw,kkk);

                for(int j=0;j<3;++j)display_vertices.push_back(kkk[j]+vect[j]);
                for(int j=0;j<3;++j)display_vertices.push_back(kkk[j]-vect[j]);


                for(int j=0;j<3;++j)display_vertices.push_back(kkk[j]+vect[j]+vecl[j]);
                for(int j=0;j<3;++j)display_vertices.push_back(kkk[j]-vect[j]+vecl[j]);

                for(int j=0;j<6;++j)display_faces.push_back(offset+rect[j]);
                for(int k=0;k<4;++k)for(int j=0;j<3;++j)display_normal.push_back(-p_ovec[j]);

                for(int k=0;k<4;++k)for(int j=0;j<4;++j)display_vcolor.push_back(cc[j]);

            }

            //box5
            if(1)for(int i=0;i<n_vertices;++i){
                if(!isShowVertices[i])continue;
                if(sparseShowField.size()!=0)if(!sparseShowField[i])continue;
                p_vec=p_vvec+i*3;p_ovec=p_ovvec+i*3;
                offset = display_vertices.size()/3;
                auto p_v = v_begin(i);
                cross(p_vec,p_ovec,orth);

                product(thicknessscale,orth,vect);
                product(widthscale,p_ovec,vecw);
                product(lengthscale,p_vec,vecl);

                add(p_v,vecl,kkk);

                for(int j=0;j<3;++j)display_vertices.push_back(kkk[j]+vect[j]+vecw[j]);
                for(int j=0;j<3;++j)display_vertices.push_back(kkk[j]-vect[j]+vecw[j]);


                for(int j=0;j<3;++j)display_vertices.push_back(kkk[j]+vect[j]-vecw[j]);
                for(int j=0;j<3;++j)display_vertices.push_back(kkk[j]-vect[j]-vecw[j]);

                for(int j=0;j<6;++j)display_faces.push_back(offset+rect[j]);
                for(int k=0;k<4;++k)for(int j=0;j<3;++j)display_normal.push_back(p_vec[j]);

                for(int k=0;k<4;++k)for(int j=0;j<4;++j)display_vcolor.push_back(cc[j]);

            }
            //box6
            if(1)for(int i=0;i<n_vertices;++i){
                if(!isShowVertices[i])continue;
                if(sparseShowField.size()!=0)if(!sparseShowField[i])continue;
                p_vec=p_vvec+i*3;p_ovec=p_ovvec+i*3;
                offset = display_vertices.size()/3;
                auto p_v = v_begin(i);
                cross(p_vec,p_ovec,orth);

                product(thicknessscale,orth,vect);
                product(widthscale,p_ovec,vecw);
                product(lengthscale,p_vec,vecl);

                //minusVec(p_v,vecl,kkk);

                for(int j=0;j<3;++j)display_vertices.push_back(p_v[j]+vect[j]+vecw[j]);
                for(int j=0;j<3;++j)display_vertices.push_back(p_v[j]-vect[j]+vecw[j]);


                for(int j=0;j<3;++j)display_vertices.push_back(p_v[j]+vect[j]-vecw[j]);
                for(int j=0;j<3;++j)display_vertices.push_back(p_v[j]-vect[j]-vecw[j]);

                for(int j=0;j<6;++j)display_faces.push_back(offset+rect[j]);
                for(int k=0;k<4;++k)for(int j=0;j<3;++j)display_normal.push_back(-p_vec[j]);

                for(int k=0;k<4;++k)for(int j=0;j<4;++j)display_vcolor.push_back(cc[j]);

            }



        };

        auto constructRect2 = [this,&widthscale,&lengthscale,&isShowVertices,rect](double *p_vvec,double *p_ovvec,unsigned char*cc){

            uint offset;
            auto p_vec = p_vvec;
            auto p_ovec = p_ovvec;
            double vecw[3],vecl[3];

            double vec[3],endp[3],nz[3] = {0,0,1};
            Eigen::AngleAxisd t;
            vector<Eigen::Vector3d>ori(cylinder.n_vertices);
            vector<Eigen::Vector3d>orivnor(cylinder.n_vertices);
            vector<Eigen::Vector3d>tranformed(cylinder.n_vertices);

            vector<Eigen::Vector3d>ori_cone(cone.n_vertices);
            vector<Eigen::Vector3d>orivnor_cone(cone.n_vertices);
            vector<Eigen::Vector3d>tranformed_cone(cone.n_vertices);

            for(int j = 0;j<cylinder.n_vertices;++j){
                auto p_spv = cylinder.v_begin(j);
                for(int k=0;k<3;++k)ori[j](k) = p_spv[k];
            }
            for(int j = 0;j<cylinder.n_vertices;++j){
                auto p_spv = cylinder.vnor_begin(j);
                for(int k=0;k<3;++k)orivnor[j](k) = p_spv[k];
            }

            for(int j = 0;j<cone.n_vertices;++j){
                auto p_spv = cone.v_begin(j);
                for(int k=0;k<3;++k)ori_cone[j](k) = p_spv[k];
            }
            for(int j = 0;j<cone.n_vertices;++j){
                auto p_spv = cone.vnor_begin(j);
                for(int k=0;k<3;++k)orivnor_cone[j](k) = p_spv[k];
            }


            for(int i=0;i<n_vertices;++i){
                if(!isShowVertices[i])continue;
                if(sparseShowField.size()!=0)if(!sparseShowField[i])continue;

                p_vec=p_vvec+i*3;p_ovec=p_ovvec+i*3;

                auto p_v = v_begin(i);
                product(lengthscale*0.25,p_vec,endp);
                add(p_v,endp,endp);


                //auto p_vn = vnor_begin(i);
                ComputeAngleAxisMatrix(t,nz,p_vec);
                Eigen::Transform<double,3,Eigen::Affine> t_m(t),t_m_cone(t);

                //t_m_cone.scale(Eigen::Vector3d(scalecone*thickness,scalecone*thickness,scalelen*veclen));


                for(int j = 0;j<cylinder.n_vertices;++j){
                    tranformed[j] = t_m*orivnor[j];
                }

                for(int j = 0;j<cylinder.n_vertices;++j){
                    for(int k=0;k<3;++k)display_normal.push_back(tranformed[j](k));
                }

                for(int j = 0;j<cone.n_vertices;++j){
                    tranformed_cone[j] = t_m_cone*orivnor_cone[j];
                }

                for(int j = 0;j<cone.n_vertices;++j){
                    for(int k=0;k<3;++k)display_normal.push_back(tranformed_cone[j](k));
                }

                t_m.scale(Eigen::Vector3d(widthscale*3,widthscale*3,lengthscale*0.5));
                for(int j = 0;j<cylinder.n_vertices;++j){
                    tranformed[j] = t_m*ori[j];
                }



                t_m_cone.scale(Eigen::Vector3d(widthscale*1,widthscale*1,lengthscale*0.3));
                for(int j = 0;j<cone.n_vertices;++j){
                    tranformed_cone[j] = t_m_cone*ori_cone[j];
                }



                auto offset = display_vertices.size()/3;
                for(int j = 0;j<cylinder.n_vertices;++j){
                    for(int k=0;k<3;++k)display_vertices.push_back(tranformed[j](k)+p_v[k]);
                }
                auto p_sf = cylinder.fv_begin(0);
                for(int j = 0;j<cylinder.n_faces*3;++j){
                    display_faces.push_back(offset+p_sf[j]);
                }
                offset = display_vertices.size()/3;
                for(int j = 0;j<cone.n_vertices;++j){
                    for(int k=0;k<3;++k)display_vertices.push_back(tranformed_cone[j](k)+endp[k]);
                }
                p_sf = cone.fv_begin(0);
                for(int j = 0;j<cone.n_faces*3;++j){
                    display_faces.push_back(offset+p_sf[j]);
                }






                for(int k=0;k<cylinder.n_vertices;++k)for(int j=0;j<4;++j)display_vcolor.push_back(cc[j]);
                for(int k=0;k<cone.n_vertices;++k)for(int j=0;j<4;++j)display_vcolor.push_back(cc[j]);
            }

        };

        auto constructRect3 = [this,&widthscale,&lengthscale,&isShowVertices,rect](double *p_vvec,double *p_ovvec,unsigned char*cc){

            uint offset;
            auto p_vec = p_vvec;
            auto p_ovec = p_ovvec;
            double vecw[3],vecl[3];

            double vec[3],endp[3],nz[3] = {0,0,1};
            for(int i=0;i<n_vertices;++i){
                if(!isShowVertices[i])continue;
                if(sparseShowField.size()!=0)if(!sparseShowField[i])continue;
                auto offset = display_vertices.size()/3;
                auto p_vec = p_vvec+i*3;
                auto p_v = v_begin(i);
                product(lengthscale,p_vec,vecl);
                for(int j = 0;j<3;++j)display_vertices.push_back(p_v[j]);
                for(int j = 0;j<3;++j)display_vertices.push_back(p_v[j]+vecl[j]);


                display_edges.push_back(offset);
                display_edges.push_back(offset+1);

                for(int k=0;k<2;++k)for(int j=0;j<4;++j)display_vcolor.push_back(cc[j]);


            }
        };

        auto constructRect4 = [this,&widthscale,&lengthscale,&isShowVertices,rect](double *p_vvec,double *p_ovvec,unsigned char*cc){

            uint offset;
            auto p_vec = p_vvec;
            auto p_ovec = p_ovvec;
            double vecw[3],vecl[3];

            double vec[3],endp[3],nz[3] = {0,0,1};
            Eigen::AngleAxisd t;
            vector<Eigen::Vector3d>ori(cylinder.n_vertices);
            vector<Eigen::Vector3d>orivnor(cylinder.n_vertices);
            vector<Eigen::Vector3d>tranformed(cylinder.n_vertices);


            for(int j = 0;j<cylinder.n_vertices;++j){
                auto p_spv = cylinder.v_begin(j);
                for(int k=0;k<3;++k)ori[j](k) = p_spv[k];
            }
            for(int j = 0;j<cylinder.n_vertices;++j){
                auto p_spv = cylinder.vnor_begin(j);
                for(int k=0;k<3;++k)orivnor[j](k) = p_spv[k];
            }


            for(int i=0;i<n_vertices;++i){
                if(!isShowVertices[i])continue;
                if(sparseShowField.size()!=0)if(!sparseShowField[i])continue;

                p_vec=p_vvec+i*3;p_ovec=p_ovvec+i*3;

                auto p_v = v_begin(i);
                product(lengthscale*0.25,p_vec,endp);
                add(p_v,endp,endp);


                //auto p_vn = vnor_begin(i);
                ComputeAngleAxisMatrix(t,nz,p_vec);
                Eigen::Transform<double,3,Eigen::Affine> t_m(t),t_m_cone(t);

                //t_m_cone.scale(Eigen::Vector3d(scalecone*thickness,scalecone*thickness,scalelen*veclen));


                for(int j = 0;j<cylinder.n_vertices;++j){
                    tranformed[j] = t_m*orivnor[j];
                }

                for(int j = 0;j<cylinder.n_vertices;++j){
                    for(int k=0;k<3;++k)display_normal.push_back(tranformed[j](k));
                }



                t_m.scale(Eigen::Vector3d(widthscale*3,widthscale*3,lengthscale*0.5));
                for(int j = 0;j<cylinder.n_vertices;++j){
                    tranformed[j] = t_m*ori[j];
                }


                auto offset = display_vertices.size()/3;
                for(int j = 0;j<cylinder.n_vertices;++j){
                    for(int k=0;k<3;++k)display_vertices.push_back(tranformed[j](k)+p_v[k]);
                }
                auto p_sf = cylinder.fv_begin(0);
                for(int j = 0;j<cylinder.n_faces*3;++j){
                    display_faces.push_back(offset+p_sf[j]);
                }


                for(int k=0;k<cylinder.n_vertices;++k)for(int j=0;j<4;++j)display_vcolor.push_back(cc[j]);

            }

        };
        //cylinder.ReScale(widthscale,widthscale,lengthscale*0.2);

        if(info.isLineMode){
            if(isVMF && vertices_field.size()!=0)constructRect4(vvec_begin(0),vovec_begin(0),blue);
            if(isfield && vertices_normal.size()!=0)constructRect4(vnor_begin(0),vvec_begin(0),red);
            if(info.isNormal && vertices_normal.size()!=0)constructRect4(vnor_begin(0),vvec_begin(0),red);
        }else{
            if(isVMF && vertices_field.size()!=0)constructRect2(vvec_begin(0),vnor_begin(0),green);
            if(isfield && vertices_field.size()!=0)constructRect(vnor_begin(0),vovec_begin(0),green);
            if(info.isNormal && vertices_normal.size()!=0)constructRect(vnor_begin(0),vvec_begin(0),green);
        }

    }else{


        //isosurface.BuildDisplay(info);
        isosurface.BuildDisplay(infoSurfDisp(info.isField,info.isNormal,info.isSurface,info.isWire,
                                             info.isSingularity,false,info.length,info.length,info.upnormal));


        display_vertices = *(isosurface.getDisplayVertices());
        display_faces = *(isosurface.getDisplayFaces());
        display_normal = *(isosurface.getDisplayVerticesNormal());
        display_vcolor = *(isosurface.getDisplayColor());

        display_edges = *(isosurface.getDisplayEdges());

        //cout<<"ccc"<<endl;

    }
    isbuilddisp = true;
}


void Volume::sparseSampling(int a){
    if(n_vertices==0)return;

    sparsecoff = a;

    list<uint>ll1,ll2;

    auto issample = [this,a,&ll1,&ll2](uint ind){

        //if(vertices_flinge[ind])return false;

        auto p_pre = &ll1;auto p_cur = &ll2;

        uint *pp;
        int num;
        p_pre->clear();p_cur->clear();
        p_pre->push_back(ind);

        for(int i=0;i<a;++i){
            while(!p_pre->empty()){
                ind = p_pre->front();
                p_pre->pop_front();
                num = vv_num(ind);
                pp=vv_begin(ind);
                for(int j=0;j<num;++j){
                    if(sparseShowField[pp[j]])return false;
                    else p_cur->push_back(pp[j]);
                }

            }
            swap(p_pre,p_cur);
        }

        return true;


    };

    sparseShowField.clear();
    sparseShowField.resize(n_vertices,false);
    vector<bool>visited(n_vertices,false);
    uint seed = 100;
    list<uint>ll3,*p_ppre,*p_pcur;
    list<uint>ll4;



    p_ppre = &ll3;p_pcur = &ll4;

    p_ppre->clear();p_pcur->clear();
    p_ppre->push_back(seed);
    visited[seed] = true;
    int iter = 0,knn = 0;
    while(true){
        ++iter;

        while(!p_ppre->empty()){
            auto ind = p_ppre->front();
            p_ppre->pop_front();

            auto vvum = vv_num(ind);
            auto p_vv = vv_begin(ind);
            if(issample(ind)){sparseShowField[ind] = true;++knn;}
            for(int i=0; i<vvum; ++i){
                if(!visited[p_vv[i]]){
                    visited[p_vv[i]] = true;
                    p_pcur->push_back(p_vv[i]);
                }
            }
        }
        //cout<<p_pcur->size()<<' ';cout<<endl;
        if(p_pcur->empty())break;
        swap(p_ppre,p_pcur);


    }

    //cout<<"sparse: "<<iter<<' '<<knn<<endl;
    //for(auto a:visited)cout<<a<<' ';cout<<endl;
}



}//n_rf
