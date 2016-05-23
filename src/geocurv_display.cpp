#include "geo_curv.h"

#include <eigen3/Eigen/Geometry>

namespace n_rf {
void ComputeAngleAxisMatrix(Eigen::AngleAxisd &m_aa,double *oriVec,double *newVec){
    Eigen::Vector3d oriVecE,newVecE,nc;
    for(int j=0;j<3;++j)oriVecE(j) = oriVec[j];
    for(int j=0;j<3;++j)newVecE(j) = newVec[j];
    nc = oriVecE.cross(newVecE);
    if(len(nc.data())<1e-6){nc(0) = 1.0;nc(1) = 0.0;nc(2) = 0.0;}
    nc.normalize();
    double angle = acos(cosine(oriVec,newVec));
    m_aa = Eigen::AngleAxisd(angle,nc);
}
void Curve::BuildDisplay(bool isrescale,bool isuseball){
    setparameters();
    //scale = 3;
    double sphere_scale = 0.20;
    //cout<<"isrescale "<<isrescale<<endl;
    //isrescale = false;
    //bool isuseball = false;


    if(isrescale){
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
        double largestdis = max(max(xmax-xmin,ymax-ymin),zmax-zmin)/2;
        //cout<<centers[0]<<"  "<<centers[1]<<"  "<<centers[2]<<"   "<<largestdis<<endl;

        //for(auto& v:vertices)v/=largestdis;
        for (int j = 0; j<n_vertices*3; j+=3)
            for(int i = 0; i<3; i++)
                vertices[j+i]=(vertices[j+i]-centers[i]);


    }

    avedis = 0.0;
    for(int i =0; i< n_edges; i++){
        auto p_v = ev_begin(i);
        avedis+=VerticesDistance(*p_v,*(p_v+1));
    }
    avedis/=n_edges;


    scale=avedis;
    sphere_scale*=avedis;
    sphere.ReScale_uniform(sphere_scale);




    display_vertices = vertices;
    //display_edges = edges;
    display_faces.clear();
    //rvector = vertices;

    display_vnormals=display_vertices;



    double *p_spv = NULL, *p_v = NULL,*p_vn = NULL;
    uint *p_sf = NULL;
    auto n_fs = sphere.n_faces*3;
    auto n_vs = sphere.n_vertices*3;
    //cout<<"ppppp: "<<display_vertices.size()<<endl;

    if(isuseball)for(int iv = 0;iv<n_vertices;++iv){
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
            display_vnormals.push_back(p_vn[j]);
        }
    }
    double vec[3],nz[3] = {0,0,1};
    Eigen::AngleAxisd t;
    vector<Eigen::Vector3d>ori(cylinder.n_vertices);
    vector<Eigen::Vector3d>tranformed(cylinder.n_vertices);
    for(int j = 0;j<cylinder.n_vertices;++j){
        p_spv = cylinder.v_begin(j);
        for(int k=0;k<3;++k)ori[j](k) = p_spv[k];
    }
    vector<Eigen::Vector3d>orinor(cylinder.n_vertices);
    for(int j = 0;j<cylinder.n_vertices;++j){
        p_spv = cylinder.vnor_begin(j);
        for(int k=0;k<3;++k)orinor[j](k) = p_spv[k];
    }

    if(isuseball)for(int iv = 0;iv<n_edges;++iv){
        int offset = display_vertices.size()/3;
        auto p_ev = ev_begin(iv);
        minusVec(v_begin(p_ev[0]),v_begin(p_ev[1]),vec);
        normalize(vec);
        ComputeAngleAxisMatrix(t,nz,vec);
        Eigen::Transform<double,3,Eigen::Affine> t_m(t);
        _VerticesMidpoint(v_begin(p_ev[0]),v_begin(p_ev[1]),vec);
        //t_m.translate(Eigen::Vector3d(vec[0],vec[1],10*vec[2]));
        double d = _VerticesDistance(v_begin(p_ev[0]),v_begin(p_ev[1]));
        t_m.scale(Eigen::Vector3d(sphere_scale,sphere_scale,d));
        for(int j = 0;j<cylinder.n_vertices;++j){
            tranformed[j] = t_m*ori[j];
        }

        for(int j = 0;j<cylinder.n_vertices;++j){
            for(int k=0;k<3;++k)display_vertices.push_back(tranformed[j](k)+vec[k]);
        }

        for(int j = 0;j<cylinder.n_vertices;++j){
            tranformed[j] = t_m*orinor[j];
        }

        for(int j = 0;j<cylinder.n_vertices;++j){
            for(int k=0;k<3;++k)display_vnormals.push_back(tranformed[j](k));
        }

        p_sf = cylinder.fv_begin(0);
        for(int j = 0;j<cylinder.n_faces*3;++j){
            display_faces.push_back(offset+p_sf[j]);
        }


    }

    n_ver_f1 = display_vertices.size();
    //cout<<"fff: "<<display_faces.size()<<endl;

    n_face_f1 = display_faces.size();

    n_vnor_f1 = display_vnormals.size();





    unsigned char red[4] = {225,0,0,255};
    unsigned char blue[4] = {0,0,225,255};
    unsigned char green[4] = {0,255,0,255};
    unsigned char gray[4] = {225,225,225,255};
    v_color.clear();
    for (int i = 0; i < n_vertices; i++)for(int j =0;j<4;++j)v_color.push_back(red[j]);

    if(isuseball){
        for (int i = 0; i <  sphere.n_vertices; i++)for(int j =0;j<4;++j)v_color.push_back(gray[j]);
        for (int i = sphere.n_vertices; i < n_vertices * sphere.n_vertices; i++)for(int j =0;j<4;++j)v_color.push_back(gray[j]);
        for (int i = 0; i < n_edges * cylinder.n_vertices; i++)for(int j =0;j<4;++j)v_color.push_back(gray[j]);
    }



    n_col_f1 = v_color.size();
    isbuild = true;


}



void Curve::BuildDisplay(bool isrescale, bool isvmf, bool isrmf, bool isnormal, bool isbinormal, bool istangent, bool isrvector, bool isovector,
                         bool isVector, bool isSurface,
                         double veclen, double thickness, double boxlen, double boxthickness, bool isSurfColor, bool isColorFirst)

{

    if(0){
        if(!isbuild)BuildDisplay(isrescale);
    }else BuildDisplay(isrescale,!isSurface);;
    if(!isBuildFrenet)EstimateFrenetFrame();
    if(!isSettangent)BuildTangent();
    display_vertices.resize(n_ver_f1);
    display_faces.resize(n_face_f1);
    display_vnormals.resize(n_vnor_f1);

    //display_edges.resize(2*n_edges);
    v_color.resize(n_col_f1);

    unsigned char red[4] = {255,0,0,255};
    unsigned char green[4] = {0,255,0,255};
    unsigned char dark[4] = {0,0,0,255};
    unsigned char blue[4] = {0,0,255,255};
    unsigned char pink[4] = {255,204,204,255};
    unsigned char gray[4] = {225,225,225,255};
    thickness*=1.5;
    veclen*=0.9;
    boxthickness*=0.15;
    boxlen*=0.5;
    bool isunitcolor = true;

    if(isCameraMode){
        cout<<"focusing: "<<focusing[0]<<' '<<focusing[1]<<' '<<focusing[2]<<endl;
        int offset = display_vertices.size()/3;
        for(int j = 0;j<sphere.n_vertices;++j){
            auto p_spv = sphere.v_begin(j);
            for(int k=0;k<3;++k)display_vertices.push_back(focusing[k]+3.*p_spv[k]);
        }
        auto p_sf = sphere.fv_begin(0);
        for(int j = 0;j<sphere.n_faces*3;++j){
            display_faces.push_back(offset+p_sf[j]);
        }
        auto p_vn = sphere.vnor_begin(0);
        for(int j = 0;j<sphere.n_vertices*3;++j){
            display_vnormals.push_back(p_vn[j]);
        }
        for(int j = 0;j<sphere.n_vertices;++j){
            for(int j =0;j<4;++j)v_color.push_back(green[j]);
        }


    }


    if(isVector){

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
        for(int j = 0;j<cone.n_vertices;++j){
            auto p_spv = cone.v_begin(j);
            for(int k=0;k<3;++k)ori_cone[j](k) = p_spv[k];
        }
        for(int j = 0;j<cylinder.n_vertices;++j){
            auto p_spv = cylinder.vnor_begin(j);
            for(int k=0;k<3;++k)orivnor[j](k) = p_spv[k];
        }
        for(int j = 0;j<cone.n_vertices;++j){
            auto p_spv = cone.vnor_begin(j);
            for(int k=0;k<3;++k)orivnor_cone[j](k) = p_spv[k];
        }

        auto constructFunction = [this,&ori,&tranformed,&ori_cone,&tranformed_cone,&veclen,&thickness,
                &orivnor,&orivnor_cone,isColorFirst,green](double *p_fn,unsigned char*cc){
            int offset = display_vertices.size()/3;
            double scalelen = scale *3, scalecyl = scale * 0.2,scalecone = scale* 0.09;
            Eigen::AngleAxisd t;
            double vec[3],endp[3],nz[3] = {0,0,1};

            for (int i = 0; i < n_vertices; i++){
                offset = display_vertices.size()/3;
                auto p_v = v_begin(i);
                auto p_rv = p_fn + i*3;
                product(scalelen*veclen,p_rv,vec);
                add(p_v,vec,vec);
                _VerticesMidpoint(p_v,vec,vec);
                product(scalelen*veclen,p_rv,endp);
                add(p_v,endp,endp);
                ComputeAngleAxisMatrix(t,nz,p_rv);
                Eigen::Transform<double,3,Eigen::Affine> t_m_cone(t),t_m_cyl(t),t_m(t);
                t_m_cone.scale(Eigen::Vector3d(scalecone*thickness,scalecone*thickness,scalelen*veclen));
                t_m_cyl.scale(Eigen::Vector3d(scalecyl*thickness,scalecyl*thickness,scalelen*veclen));
                for(int j = 0;j<cylinder.n_vertices;++j){
                    tranformed[j] = t_m_cyl*ori[j];
                }

                for(int j = 0;j<cone.n_vertices;++j){
                    tranformed_cone[j] = t_m_cone*ori_cone[j];
                }
                for(int j = 0;j<cylinder.n_vertices;++j){
                    for(int k=0;k<3;++k)display_vertices.push_back(tranformed[j](k)+vec[k]);
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


                for(int j = 0;j<cylinder.n_vertices;++j){
                    tranformed[j] = t_m*orivnor[j];
                }
                for(int j = 0;j<cone.n_vertices;++j){
                    tranformed_cone[j] = t_m*orivnor_cone[j];
                }
                for(int j = 0;j<cylinder.n_vertices;++j){
                    for(int k=0;k<3;++k)display_vnormals.push_back(tranformed[j](k));
                }
                for(int j = 0;j<cone.n_vertices;++j){
                    for(int k=0;k<3;++k)display_vnormals.push_back(tranformed_cone[j](k));
                }

            }
            //cout<<"cone: "<<cone.n_vertices<<" "<<cone.n_faces<<endl;
            if(isColorFirst){
                for (int i = 0; i < cylinder.n_vertices; i++)for(int j =0;j<4;++j)v_color.push_back(green[j]);
                for (int i = 0; i < cone.n_vertices; i++)for(int j =0;j<4;++j)v_color.push_back(green[j]);
                for (int i = cylinder.n_vertices; i < n_vertices * cylinder.n_vertices; i++)for(int j =0;j<4;++j)v_color.push_back(cc[j]);
                for (int i = cone.n_vertices; i < n_vertices * cone.n_vertices; i++)for(int j =0;j<4;++j)v_color.push_back(cc[j]);
            }else{
                for (int i = 0; i < n_vertices * cylinder.n_vertices; i++)for(int j =0;j<4;++j)v_color.push_back(cc[j]);
                for (int i = 0; i < n_vertices * cone.n_vertices; i++)for(int j =0;j<4;++j)v_color.push_back(cc[j]);
            }
        };
        auto constructFunction_single = [this,&ori,&tranformed,&ori_cone,&tranformed_cone,&veclen,&thickness,
                &orivnor,&orivnor_cone](double *p_rv,int ind, unsigned char*cc){
            int offset = display_vertices.size()/3;
            double scale1 = scale *3.2, scale2 = scale * 0.13;
            Eigen::AngleAxisd t;
            double vec[3],endp[3],nz[3] = {0,0,1};

            offset = display_vertices.size()/3;
            auto p_v = v_begin(ind);

            product(scale1*veclen,p_rv,vec);
            add(p_v,vec,vec);
            _VerticesMidpoint(p_v,vec,vec);
            product(scale1*veclen,p_rv,endp);
            add(p_v,endp,endp);
            ComputeAngleAxisMatrix(t,nz,p_rv);
            Eigen::Transform<double,3,Eigen::Affine> t_m(t);
            t_m.scale(Eigen::Vector3d(scale2*thickness,scale2*thickness,scale1*veclen));
            for(int j = 0;j<cylinder.n_vertices;++j){
                tranformed[j] = t_m*ori[j];
            }
            for(int j = 0;j<cone.n_vertices;++j){
                tranformed_cone[j] = t_m*ori_cone[j];
            }
            for(int j = 0;j<cylinder.n_vertices;++j){
                for(int k=0;k<3;++k)display_vertices.push_back(tranformed[j](k)+vec[k]);
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


            for(int j = 0;j<cylinder.n_vertices;++j){
                tranformed[j] = t_m*orivnor[j];
            }
            for(int j = 0;j<cone.n_vertices;++j){
                tranformed_cone[j] = t_m*orivnor_cone[j];
            }
            for(int j = 0;j<cylinder.n_vertices;++j){
                for(int k=0;k<3;++k)display_vnormals.push_back(tranformed[j](k));
            }
            for(int j = 0;j<cone.n_vertices;++j){
                for(int k=0;k<3;++k)display_vnormals.push_back(tranformed_cone[j](k));
            }
            //cout<<"cone: "<<cone.n_vertices<<" "<<cone.n_faces<<endl;
            for (int i = 0; i < cylinder.n_vertices; i++)for(int j =0;j<4;++j)v_color.push_back(cc[j]);
            for (int i = 0; i < cone.n_vertices; i++)for(int j =0;j<4;++j)v_color.push_back(cc[j]);
        };

        if(isrvector){
            //cout<<"isVector"<<endl;
            if(isunitcolor){
                if(isvmf && vertices_field.size()!=0)constructFunction(vvec_begin(0),blue);
                //if(isvmf && vertices_f_field.size()!=0)constructFunction(vfvec_begin(0),blue);
                //if(isrmf && rvectorRMF.size()!=0)constructFunction(rr_begin(0),blue);
                if(isnormal && Frenetnormal.size()!=0)constructFunction(nor_begin(0),blue);
                if(isbinormal && Frenetbinormal.size()!=0)constructFunction(binor_begin(0),blue);
            }else{
                if(isvmf && vertices_field.size()!=0)constructFunction(vvec_begin(0),green);
                //if(isrmf && rvectorRMF.size()!=0)constructFunction(rr_begin(0),red);
                if(isnormal && Frenetnormal.size()!=0)constructFunction(nor_begin(0),pink);
                if(isbinormal && Frenetbinormal.size()!=0)constructFunction(binor_begin(0),blue);
            }
        }
        if(isovector){
            if(isvmf && vertices_ofield.size()!=0)constructFunction(vovec_begin(0),green);
            // if(isrmf && rvectorRMF.size()!=0)constructFunction(or_begin(0),red);
            if(isnormal && Frenetnormal.size()!=0)constructFunction(nor_begin(0),pink);
            if(isbinormal && Frenetbinormal.size()!=0)constructFunction(binor_begin(0),blue);
        }

        if(istangent && tangent.size()!=0)constructFunction(t_begin(0),red);
        //        for(int i =0;i<n_vertices;++i)if(is_constrained(i))
        //        {
        //            //cout<<"constrain: "<<i<<endl;
        //            constructFunction_single(rvc_begin(i),i,dark);
        //        }

    }

    if(isSurface){
        auto constructSurface2 = [this,&veclen,&thickness](double *p_fn,unsigned char*cc){

            double scale1 = scale *3*veclen;
            double vec[3];
            int  offset = display_vertices.size()/3;
            for (int i = 0; i < n_vertices; i++){
                auto p_v = v_begin(i);
                auto p_rv = p_fn + i*3;
                product(scale1,p_rv,vec);
                for(int k=0;k<3;++k)display_vertices.push_back(p_v[k]+vec[k]);
                for(int k=0;k<3;++k)display_vertices.push_back(p_v[k]-vec[k]);
            }

            for (int i = 0; i < n_vertices-1; i++){
                int b = offset+2*i;
                display_faces.push_back(b);
                display_faces.push_back(b+1);
                display_faces.push_back(b+2);
                display_faces.push_back(b+2);
                display_faces.push_back(b+1);
                display_faces.push_back(b+3);
            }
            if(isloop){
                int b = offset+2*(n_vertices-1);
                display_faces.push_back(b);
                display_faces.push_back(b+1);
                display_faces.push_back(offset);
                display_faces.push_back(offset);
                display_faces.push_back(b+1);
                display_faces.push_back(offset+1);
            }


            for (int i = 0; i < n_vertices * 2; i++)for(int j =0;j<4;++j)v_color.push_back(cc[j]);
        };


        auto constructSurface = [this,&veclen,&thickness](double *p_fn,unsigned char*cc){

            double scale1 = scale *3*veclen;
            double vec[3],ovec[3];
            double nor1[3],nor2[3];
            double scale2 = scale*1*thickness;
            int  offset = display_vertices.size()/3;


            for (int i = 0; i < n_vertices; i++){
                auto p_v = v_begin(i);
                auto p_rv = p_fn + i*3;
                cross(t_begin(i),p_rv,ovec);
                product(scale1,p_rv,vec);
                product(scale2,ovec,ovec);
                add(vec,ovec,nor1);minusVec(ovec,vec,nor2);

                for(int k=0;k<3;++k)display_vertices.push_back(p_v[k]+nor1[k]);
                for(int k=0;k<3;++k)display_vertices.push_back(p_v[k]+nor2[k]);
                for(int k=0;k<3;++k)display_vertices.push_back(p_v[k]-nor2[k]);
                for(int k=0;k<3;++k)display_vertices.push_back(p_v[k]-nor1[k]);

                //cout<<"adsadsadasdasasdasd"<<endl;
                normalize(nor1);
                normalize(nor2);


                for(int k=0;k<3;++k)display_vnormals.push_back(nor1[k]);
                for(int k=0;k<3;++k)display_vnormals.push_back(nor2[k]);
                for(int k=0;k<3;++k)display_vnormals.push_back(-nor2[k]);
                for(int k=0;k<3;++k)display_vnormals.push_back(-nor1[k]);
            }

            int  ff[] = {4,1,0,5,4,1,2,3,6,3,6,7,0,2,4,2,4,6,1,3,5,3,5,7};
            //int  ff[] = {0,1,4,4,1,5,2,3,6,6,3,7,0,2,4,4,2,6,1,3,5,5,3,7};
            for (int i = 0; i < n_vertices-1; i++){
                int b = offset+4*i;
                for(int j=0;j<24;++j)display_faces.push_back(b+ff[j]);


            }
            if(isloop){
                int b = offset+2*(n_vertices-1);
                display_faces.push_back(b);
                display_faces.push_back(b+1);
                display_faces.push_back(offset);
                display_faces.push_back(offset);
                display_faces.push_back(b+1);
                display_faces.push_back(offset+1);


            }else{
                int b = offset;
                display_faces.push_back(b);display_faces.push_back(b+1);display_faces.push_back(b+2);
                display_faces.push_back(b+2);display_faces.push_back(b+1);display_faces.push_back(b+3);


                b = offset+4*(n_vertices-1);
                display_faces.push_back(b);display_faces.push_back(b+1);display_faces.push_back(b+2);
                display_faces.push_back(b+2);display_faces.push_back(b+1);display_faces.push_back(b+3);

            }


            for (int i = 0; i < n_vertices * 4; i++)for(int j =0;j<4;++j)v_color.push_back(cc[j]);
        };

        auto constructSurface3 = [this,&boxlen,&boxthickness,isSurfColor](double *p_fn,unsigned char*cc){

            double scale1 = scale *3*boxlen;
            double vec[3],ovec[3];
            double nor1[3],nor2[3];
            double scale2 = scale*1*boxthickness;
            int  offset = display_vertices.size()/3;

            for (int i = 0; i < n_vertices; i++){
                auto p_v = v_begin(i);
                auto p_rv = p_fn + i*3;
                cross(t_begin(i),p_rv,ovec);
                product(scale1,p_rv,vec);
                product(scale2,ovec,ovec);
                add(vec,ovec,nor1);
                minusVec(vec,ovec,nor2);

                for(int j=0;j<2;++j){
                    for(int k=0;k<3;++k)display_vertices.push_back(p_v[k]+nor1[k]);
                    for(int k=0;k<3;++k)display_vertices.push_back(p_v[k]+nor2[k]);
                    for(int k=0;k<3;++k)display_vertices.push_back(p_v[k]-nor2[k]);
                    for(int k=0;k<3;++k)display_vertices.push_back(p_v[k]-nor1[k]);
                }


                for(int k=0;k<3;++k)display_vnormals.push_back(p_rv[k]);
                for(int k=0;k<3;++k)display_vnormals.push_back(p_rv[k]);
                for(int k=0;k<3;++k)display_vnormals.push_back(-p_rv[k]);
                for(int k=0;k<3;++k)display_vnormals.push_back(-p_rv[k]);

                for(int k=0;k<3;++k)display_vnormals.push_back(ovec[k]);
                for(int k=0;k<3;++k)display_vnormals.push_back(-ovec[k]);
                for(int k=0;k<3;++k)display_vnormals.push_back(ovec[k]);
                for(int k=0;k<3;++k)display_vnormals.push_back(-ovec[k]);
            }

            int  ff[] = {0,1,8,8,1,9,2,3,10,10,3,11,
                         4,6,12,12,6,14,5,7,13,13,7,15};
            for (int i = 0; i < n_vertices-1; i++){
                int b = offset+8*i;
                for(int j=0;j<24;++j)display_faces.push_back(b+ff[j]);


            }
            int newvnum = n_vertices*8;
            if(1)
                if(isloop){
                    for(int i=0;i<2;++i){
                        int e = offset+8*(n_vertices-1)+i*2;
                        int b = offset +i*2;
                        display_faces.push_back(b);
                        display_faces.push_back(b+1);
                        display_faces.push_back(e);
                        display_faces.push_back(e);
                        display_faces.push_back(b+1);
                        display_faces.push_back(e+1);
                    }
                    for(int i=0;i<2;++i){
                        int e = offset+8*(n_vertices-1)+i+4;
                        int b = offset + i +4;
                        display_faces.push_back(b);
                        display_faces.push_back(b+2);
                        display_faces.push_back(e);
                        display_faces.push_back(e);
                        display_faces.push_back(b+2);
                        display_faces.push_back(e+2);
                    }


                }else{
                    int b = offset;

                    display_faces.push_back(b);display_faces.push_back(b+1);display_faces.push_back(b+2);
                    display_faces.push_back(b+2);display_faces.push_back(b+1);display_faces.push_back(b+3);


                    b = offset+8*(n_vertices-1);
                    display_faces.push_back(b);display_faces.push_back(b+1);display_faces.push_back(b+2);
                    display_faces.push_back(b+2);display_faces.push_back(b+1);display_faces.push_back(b+3);

                }


            if(!isSurfColor)for (int i = 0; i < newvnum; i++)for(int j =0;j<4;++j)v_color.push_back(cc[j]);
            else{
                uchar colo[4] = {0,100,0,255};
                for(int i=0;i<n_vertices;++i){
                    //int kkkccc = min(255,int(normCurvature[i]*255.));
                    //colo[0] = kkkccc;colo[1] = 255-kkkccc;
                    //ThresholdColoring(160.*(1.-normCurvature[i]),colo);

                    for(int k=0;k<8;++k)for(int j =0;j<4;++j)v_color.push_back(colo[j]);
                }
            }
        };






        auto *pcolor = gray;
        if(isovector){
            if(1){

                if(isvmf && vertices_field.size()!=0)constructSurface3(vvec_begin(0),pcolor);
                //if(isrmf && rvectorRMF.size()!=0)constructSurface3(rr_begin(0),pcolor);
                if(isnormal && Frenetnormal.size()!=0)constructSurface3(nor_begin(0),pcolor);
                if(isbinormal && Frenetbinormal.size()!=0)constructSurface3(binor_begin(0),pcolor);
            }else{
                if(isvmf && vertices_field.size()!=0)constructSurface3(vvec_begin(0),green);
                //if(isrmf && rvectorRMF.size()!=0)constructSurface3(rr_begin(0),red);
                if(isnormal && Frenetnormal.size()!=0)constructSurface3(nor_begin(0),blue);
                if(isbinormal && Frenetbinormal.size()!=0)constructSurface3(binor_begin(0),pink);
            }
        }
        if(isrvector){
            if(1){
                if(isvmf && vertices_ofield.size()!=0)constructSurface3(vovec_begin(0),pcolor);
                //if(isrmf && rvectorRMF.size()!=0)constructSurface3(or_begin(0),pcolor);
                if(isnormal && Frenetnormal.size()!=0)constructSurface3(binor_begin(0),pcolor);
                if(isbinormal && Frenetbinormal.size()!=0)constructSurface3(nor_begin(0),pcolor);
            }else{
                if(isvmf && vertices_ofield.size()!=0)constructSurface3(vovec_begin(0),green);
                //if(isrmf && rvectorRMF.size()!=0)constructSurface3(or_begin(0),red);
                if(isnormal && Frenetnormal.size()!=0)constructSurface3(nor_begin(0),blue);
                if(isbinormal && Frenetbinormal.size()!=0)constructSurface3(binor_begin(0),pink);
            }
        }
    }

    //    display_vnormals = display_vertices;
    //    for(int i=0;i<display_vnormals.size();i+=3){
    //        normalize(display_vnormals.data()+i);
    //    }

    isbuild = true;


}


bool Curve::load(){
    if(!isload){
        //sphere.readfile(string("/Users/Research/Geometry/MM/QTProgram/sphere.off"));
        sphere.createToy(3);
        //cylinder.readfile(string("/Users/Research/Geometry/MM/QTProgram/cylinder.off"));
        cylinder.createToy(2);
        //cone.readfile(string("/Users/Research/Geometry/MM/QTProgram/cone.off"));
        cone.createToy(1);
        cone.ReScale(2.5,2.5,0.5);
        cylinder.ReScale(0.5,0.5,1.0);

        sphere.BuildNeighborTable();
        cylinder.BuildNeighborTable();
        cone.BuildNeighborTable();

        sphere.ComputeFaceNormal(true);
        cylinder.ComputeFaceNormal(true);
        cone.ComputeFaceNormal(true);


        isload = true;
    }

    return true;
    //sphere.readfile(string("/Users/Research/Geometry/RMFpro/cube.off"));
}

}//n_rf
