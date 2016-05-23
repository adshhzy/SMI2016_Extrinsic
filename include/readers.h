#ifndef READERS_H
#define READERS_H


#include<vector>
#include<string>
using namespace std;

bool readOffFile(string filename,vector<double>&vertices,vector<unsigned int>&faces2vertices);


bool writeOffFile(string filename,const vector<double>&vertices,const vector<unsigned int>&faces2vertices);



bool readObjFile(string filename,vector<double>&vertices,vector<unsigned int>&faces2vertices,vector<double>vertices_normal);

bool writeObjFile(string filename,const vector<double>&vertices,const vector<unsigned int>&faces2vertices);


bool readSurfFile(string filename,vector<double>&vertices,vector<unsigned int>&faces2vertices,vector<double>&vertices_field);

bool writeSurfFile(string filename,const vector<double>&vertices,const vector<unsigned int>&faces2vertices,const vector<double>&vertices_field);


bool readCurfFile(string filename, vector<double>&vertices, vector<unsigned int>&edges2vertices, vector<double>&vertices_field, vector<double> &vertices_tangent);


bool writeCurfFile(string filename, const vector<double>&vertices, const vector<unsigned int>&edges2vertices, vector<double>&vertices_field , vector<double> &vertices_tangent);


bool readVolfFile(string filename, vector<double>&vertices, vector<unsigned int>&tets2vertices, vector<double> &vertices_normal, vector<double> &vertices_field);


bool writeVolfFile(string filename, const vector<double>&vertices, const vector<unsigned int>&tets2vertices, vector<double> &vertices_normal, vector<double> &vertices_field);





#endif // READERS_H
