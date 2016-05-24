#ifndef INFOSTRUCT_H
#define INFOSTRUCT_H



class infoSet{

public:
    bool isrescale;
    bool isvmf;
    bool isrmf;
    bool isnnormal;
    bool isbinormal;
    bool istangent;
    bool isrvector;
    bool isovector;
    bool isVector;
    bool isSurface;
    double veclen;
    double vecthickness;
    double boxlen;
    double boxthickness;

    bool isshowall;
    int pickid;

    bool isSurfColor;
    bool isColorFirst;

    infoSet(){}
    infoSet(bool isrescale,bool isvmf, bool isrmf, bool isnnormal,bool isbinormal, bool istangent,
            bool isrvector,bool isovector,
            bool isVector, bool isSurface,
            double veclen,double vecthickness,
            double boxlen,double boxthickness,
            bool isshowall,int pickid,
            bool isSurfColor, bool isColorFirst):
        isrescale(isrescale),
        isvmf(isvmf),isrmf(isrmf),
        isnnormal(isnnormal),isbinormal(isbinormal),istangent(istangent),
        isrvector(isrvector),isovector(isovector),
        isVector(isVector),isSurface(isSurface),
        veclen(veclen),vecthickness(vecthickness),
        boxlen(boxlen),boxthickness(boxthickness),
        isshowall(isshowall),pickid(pickid),
        isSurfColor(isSurfColor),isColorFirst(isColorFirst)
    {}

    infoSet(double veclen,double vecthickness,
            double boxlen,double boxthickness):
        isrescale(false),
        isvmf(true),isrmf(false),
        isnnormal(false),isbinormal(false),istangent(false),
        isrvector(true),isovector(false),
        isVector(true),isSurface(true),
        veclen(veclen),vecthickness(vecthickness),
        boxlen(boxlen),boxthickness(boxthickness),
        isshowall(true),pickid(false),
        isSurfColor(false),isColorFirst(false)
    {}


    infoSet(int veclen,int vecthickness,
            int boxlen,int boxthickness):
        isrescale(false),
        isvmf(true),isrmf(false),
        isnnormal(false),isbinormal(false),istangent(false),
        isrvector(true),isovector(false),
        isVector(true),isSurface(true),
        veclen(veclen/100.),vecthickness(vecthickness/100.),
        boxlen(boxlen/100.),boxthickness(boxthickness/100.),
        isshowall(true),pickid(false),
        isSurfColor(false),isColorFirst(false)
    {}



};

class infoFrame{
public:
    int init_div;
    int subdiv;
    int iter;
    bool isrefine;
    bool isrotate;
    bool isopt;
    int mode;
    double alpha;
    bool isInverseInit;
    bool isRandomInit;
    bool isReInit;

    infoFrame(){}

    infoFrame(int init_div, int subdiv,int iter,bool isrefine,bool isrotate,bool isopt,int mode,double alpha,
              bool isReInit = false,bool isRandomInit = false,bool isInverseInit = false):
        init_div(init_div),subdiv(subdiv),iter(iter),
        isrefine(isrefine),isrotate(isrotate),isopt(isopt),
        mode(mode),alpha(alpha),isInverseInit(isInverseInit),
        isRandomInit(isRandomInit),isReInit(isReInit)
    {}


};



class infoSurfDisp{
public:
    bool isField;
    bool isNormal;
    bool isSurface;
    bool isWire;
    bool isSingularity;
    bool ismark;
    int length;
    int width;
    int upnormal;

    infoSurfDisp():length(30),upnormal(30){}

    infoSurfDisp(bool isField,bool isNormal,bool isSurface,bool isWire,bool isSingularity,bool ismark,int length,int width,int upnormal):
        isField(isField),isNormal(isNormal),isSurface(isSurface),isWire(isWire),isSingularity(isSingularity),ismark(ismark),length(length),width(width),upnormal(upnormal)
    {}


};

class infoSurfCompute{
public:
    bool isGS;

    infoSurfCompute():isGS(false){}

    infoSurfCompute(bool isGS):
        isGS(isGS)
    {}


};

class infoSurfOutput{
public:
    int n_vertices;
    int n_faces;
    int n_edges;
    double OptTime;
    double FinalEnergy;
    int numOfSingularity;

};

class infoVolDisp{
public:
    bool isField;
    bool isNormal;
    bool isSurface;
    bool isWire;
    bool isSingularity;
    bool isIso;
    int length;
    int upnormal;
    bool isSlice;
    int xyz;
    int slicenum;
    double thickness;
    bool isLineMode;
    //int upnormal;

    infoVolDisp():isIso(false),length(30){}

    infoVolDisp(bool isField,bool isNormal,bool isSurface,bool isWire,bool isSingularity,
                bool isIso,int length,int upnormal,
                bool isSlice,int xyz,int slicenum,double thickness,bool isLineMode):
        isField(isField),isNormal(isNormal),isSurface(isSurface),isWire(isWire),isSingularity(isSingularity),
        isIso(isIso),length(length),upnormal(upnormal),
        isSlice(isSlice),xyz(xyz),slicenum(slicenum),thickness(thickness),isLineMode(isLineMode)
    {}

    infoVolDisp(int length,int upnormal,double thickness,bool isLineMode = false):
        isField(false),isNormal(true),isSurface(false),isWire(false),isSingularity(false),
        isIso(false),length(length),upnormal(upnormal),
        isSlice(false),xyz(0),slicenum(0),thickness(thickness),isLineMode(isLineMode)
    {}


};

#endif // INFOSTRUCT_H
