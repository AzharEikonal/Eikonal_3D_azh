#include <iostream>
#include <fstream> 
#include <vector> 
#include <sstream>
#include <set>
#include <bits/stdc++.h>
#include <Eigen/Dense>
// #include <openacc.h>

using namespace Eigen;
using namespace std;

struct Vertex{
    double x,y,z;
};

struct Tetrahedron{
    int v1, v2, v3, v4;
};

// read values as doubles separated by space in a line.

vector<double> splitstring(const string& str){
    std::vector<double> values;
    std::istringstream iss(str);
    double val;
    while (iss >>val) {
        values.push_back(val);
    }
    return values;
}

bool readVTKfile(string filename, vector<Vertex> &vertices, vector<Tetrahedron> &Tetrahedra){
    ifstream file(filename);
    if (!file.is_open()){
        cerr<< "error opening file: "<< filename<<endl;
        return false;
    }
    string line;
    while(getline(file, line)){
        if(line.find("POINTS")!= string::npos){
            int numvertices;
            istringstream iss(line);
            iss >> line>> numvertices;
            vertices.reserve(numvertices);

            for (int i=0; i<numvertices; i++){
                getline(file, line);
                vector<double> coords =splitstring(line);
                if (coords.size()!=3){
                    cerr<< "Not able to read the line"<< i+1<<endl;
                    return false;
                }

                Vertex vertex;
                vertex.x= coords[0];
                vertex.y= coords[1];
                vertex.z = coords[2];
                vertices.push_back(vertex);
            }
        }
        else if (line.find("CELLS") !=string:: npos){
            int numtets;
            istringstream iss(line);
            iss>> line>> numtets;
            Tetrahedra.reserve(numtets);
            for(int i=0; i<numtets; i++){
                getline(file, line);
                vector<double> cellindices= splitstring(line);
                if (cellindices.size()!=5 && cellindices[0] !=4){
                    cerr<< "Not able to read tetrahedra at line: "<< i+1<<endl;
                }

                //we cast indices as integers.
                Tetrahedron Tet;
                Tet.v1= static_cast<int> (cellindices[1]);
                Tet.v2 =static_cast<int> (cellindices[2]);
                Tet.v3= static_cast<int> (cellindices[3]);
                Tet. v4= static_cast<int> (cellindices[4]);
                Tetrahedra.push_back(Tet);
            }
        }
    }
    return true;
}

bool readvertices(string filename, vector<Vertex> &vertices){
    ifstream file(filename);
    if(!file.is_open()){
        cerr<< "error in opening fiel: "<< filename<<endl;
        return false;
    }
    string line;
    while(getline(file, line)){
        vector<double> coords= splitstring(line);
        if (coords.size()!=3){
            cerr<< "Not able to read the line"<<endl;
            return false;
        }
        Vertex vertex;
        vertex.x=coords[0];
        vertex.y= coords[1];
        vertex.z= coords[2];
        vertices.push_back(vertex);
    }
    return true;
}

bool readtets(string filename, vector<Tetrahedron> &Tetrahedra){
    ifstream file(filename);
    if(!file.is_open()){
        cerr<< "error in opening file: "<< filename<<endl;
        return false;
    }
    string line;
    while(getline(file, line)){
        vector<double> cellindices= splitstring(line);
        if (cellindices.size()!=4){
            cerr<< "Not able to read tetrahedra at line: "<<endl;
        }

        //we cast indices as integers.
        Tetrahedron Tet;
        Tet.v1= static_cast<int> (cellindices[0]);
        Tet.v2 =static_cast<int> (cellindices[1]);
        Tet.v3= static_cast<int>(cellindices[2]);
        Tet. v4= static_cast<int> (cellindices[3]);
        Tetrahedra.push_back(Tet);
    }
    return true;
}

//one ring of a vertex in mesh, returns indices of vertices. 
vector<int> onering1(vector<Tetrahedron> &tetrahedra, int ind_vertex){
        vector<int> one_ring;
        // #pragma acc parallel loop
        for(int i=0; i< tetrahedra.size(); i++){
            if (ind_vertex == tetrahedra[i].v1){
                one_ring.push_back(tetrahedra[i].v1);
                one_ring.push_back(tetrahedra[i].v2);
                one_ring.push_back(tetrahedra[i].v3);
                one_ring.push_back(tetrahedra[i].v4);
            }
            else if (ind_vertex == tetrahedra[i].v2){
                one_ring.push_back(tetrahedra[i].v1);
                one_ring.push_back(tetrahedra[i].v2);
                one_ring.push_back(tetrahedra[i].v3);
                one_ring.push_back(tetrahedra[i].v4);
            }
            else if (ind_vertex== tetrahedra[i].v3){
                one_ring.push_back(tetrahedra[i].v1);
                one_ring.push_back(tetrahedra[i].v2);
                one_ring.push_back(tetrahedra[i].v3);
                one_ring.push_back(tetrahedra[i].v4);
            }
            else if (ind_vertex ==tetrahedra[i].v4){
                one_ring.push_back(tetrahedra[i].v1);
                one_ring.push_back(tetrahedra[i].v2);
                one_ring.push_back(tetrahedra[i].v3);
                one_ring.push_back(tetrahedra[i].v4);
            }
            else {
                continue;
            }
        }
        return one_ring;
}
vector<int> onering2(vector<Tetrahedron> &tetrahedra, int index_v){
    set<int> one_ring;
    // #pragma acc parallel loop
    for(Tetrahedron t: tetrahedra){
        if(t.v1== index_v){
            one_ring.insert(t.v2);
            one_ring.insert(t.v3);
            one_ring.insert(t.v4);
        }
        else if(t.v2 == index_v){
            one_ring.insert(t.v1);
            one_ring.insert(t.v3);
            one_ring.insert(t.v4);
        }
        else if(t.v3 == index_v){
            one_ring.insert(t.v1);
            one_ring.insert(t.v2);
            one_ring.insert(t.v4);
        }
        else if(t.v4 == index_v){
            one_ring.insert(t.v1);
            one_ring.insert(t.v2);
            one_ring.insert(t.v3);
        }
    }
    vector<int> one_ring_vec;
    for(int i: one_ring){
        one_ring_vec.push_back(i);
    }
    return one_ring_vec;
}

double sq(double n){
    return n*n;
}

double norm(Vertex v) {
    return sqrt(v.x*v.x+ v.y*v.y + v.z*v.z);
}

// double f1(double x, double y, Vertex v1, Vertex v2, Vertex v3, Vertex v4, double phi1, double phi2, double phi3){
//     Vector3d u1, u2, u3, u4;
//     u1 << v1.x, v1.y, v1.z;
//     u2 <<v2.x, v2.y, v2.z;
//     u3<< v3.x, v3.y, v3.z;
//     u4 <<v4.x, v4.y, v4.z;

//     Vector3d e54;
//     Vector3d e13;
//     e13= u3-u1;
//     e54= u4-x*u1-y*u2-(1-x-y)*u3;
//     double e54te54= e54.transpose()*e54;
//     return (phi3-phi1)*sqrt(e54te54)- e54.transpose()*e13;
// }
// double f2(double x, double y, Vertex v1, Vertex v2, Vertex v3, Vertex v4, double phi1, double phi2, double phi3){
//     Vector3d u1, u2, u3, u4;
//     u1 << v1.x, v1.y, v1.z;
//     u2 <<v2.x, v2.y, v2.z;
//     u3<< v3.x, v3.y, v3.z;
//     u4 <<v4.x, v4.y, v4.z;

//     Vector3d e54= u4-x*u1-y*u2-(1-x-y)*u3;
//     Vector3d e13= u3-u1;
//     Vector3d e23= u3-u2;
//     double phi32= phi3-phi2;
//     double phi13= phi3-phi1;

//     double e54te13=e54.transpose()*e13;
//     double e54te23= e54.transpose()*e23;
//     return phi32*e54te13- phi13*e54te23;
// }

// double f1(double x, double y, 
//     Vertex v1, Vertex v2, Vertex v3, Vertex v4, double phi1, double phi2, double phi3){
//     Vertex e_1_3; 
//     e_1_3.x= v3.x-v1.x;
//     e_1_3.y = v3.y-v1.y;
//     e_1_3.z = v3.z-v1.z;

//     Vertex e_2_3;
//     e_2_3.x= v3.x-v2.x;
//     e_2_3.y = v3.y-v2.y;
//     e_2_3.z = v3.z-v2.z;

//     Vertex e_3_4;
//     e_3_4.x= v4.x-v3.x;
//     e_3_4.y = v4.y-v3.y;
//     e_3_4.z = v4.z-v3.z;
//     double phi13= phi3-phi1;

//     double a= sq(phi13)*(e_1_3.x*e_1_3.x + e_1_3.y*e_1_3.y + e_1_3.z*e_1_3.z)- sq(e_1_3.x*e_1_3.x+ e_1_3.y*e_1_3.y+ e_1_3.z*e_1_3.z);
//     double b = sq(phi13)*(e_2_3.x* e_2_3.x+ e_2_3.y*e_2_3.y+ e_2_3.z*e_2_3.z) -sq(e_2_3.x*e_1_3.x + e_2_3.y*e_1_3.y + e_2_3.z*e_1_3.z);
//     double c= sq(phi13)*(2*(e_1_3.x*e_2_3.x+ e_1_3.y*e_2_3.y+ e_1_3.z*e_2_3.z)) - 2*((e_1_3.x*e_1_3.x+ e_1_3.y*e_1_3.y+ e_1_3.z*e_1_3.z)*(e_2_3.x*e_1_3.x + e_2_3.y*e_1_3.y + e_2_3.z*e_1_3.z));
//     double d= sq(phi13)*(e_1_3.x* e_3_4.x+ e_1_3.y*e_3_4.y+ e_1_3.z*e_3_4.z) - 2*((e_3_4.x*e_1_3.x+ e_3_4.y*e_1_3.y+ e_3_4.z*e_1_3.z)*(e_1_3.x*e_1_3.x+ e_1_3.y*e_1_3.y+ e_1_3.z*e_1_3.z));
//     double e= sq(phi13)*(e_2_3.x* e_3_4.x+ e_2_3.y*e_3_4.y+ e_2_3.z*e_3_4.z) - 2*((e_2_3.x*e_1_3.x+ e_2_3.y*e_1_3.y+ e_2_3.z*e_1_3.z)*(e_3_4.x*e_1_3.x+ e_3_4.y*e_1_3.y+ e_3_4.z*e_1_3.z));
//     double f= sq(phi13)*(e_3_4.x*e_3_4.x+ e_3_4.y*e_3_4.y+ e_3_4.z*e_3_4.z) - sq(e_3_4.x*e_1_3.x+ e_3_4.y*e_1_3.y+ e_3_4.z*e_1_3.z);

//     return a*x*x+ b*y*y+c*x*y+d*x+e*y+f;
// }

// double f2(double x, double y, 
//     Vertex v1, Vertex v2, Vertex v3, Vertex v4, double phi1, double phi2, double phi3){
//     Vertex e_1_3; 
//     e_1_3.x= v3.x-v1.x;
//     e_1_3.y = v3.y-v1.y;
//     e_1_3.z = v3.z-v1.z;

//     Vertex e_2_3;
//     e_2_3.x= v3.x-v2.x;
//     e_2_3.y = v3.y-v2.y;
//     e_2_3.z = v3.z-v2.z;

//     Vertex e_3_4;
//     e_3_4.x= v4.x-v3.x;
//     e_3_4.y = v4.y-v3.y;
//     e_3_4.z = v4.z-v3.z;

//     double d= (phi3-phi2)*(e_1_3.x*e_1_3.x+ e_1_3.y*e_1_3.y+ e_1_3.z*e_1_3.z) - (phi3-phi1)*(e_1_3.x*e_2_3.x+ e_1_3.y*e_2_3.y+ e_1_3.z*e_2_3.z);
//     double e= (phi3-phi2)*(e_2_3.x*e_1_3.x+ e_2_3.y*e_1_3.y+ e_2_3.z*e_1_3.z) - (phi3-phi1)*(e_2_3.x*e_2_3.x+ e_2_3.y*e_2_3.y+ e_2_3.z*e_2_3.z);
//     double f= (phi3-phi2)*(e_3_4.x*e_1_3.x+ e_3_4.y*e_1_3.y+ e_3_4.z*e_1_3.z) - (phi3-phi1)*(e_3_4.x*e_2_3.x+ e_3_4.y*e_2_3.y+ e_3_4.z*e_2_3.z);

//     return d*x+e*y+f;
// }

// double df1_dx(double x, double y, Vertex v1, Vertex v2, Vertex v3, Vertex v4, double phi1, double phi2, double phi3){
//     Vertex e_1_3; 
//     e_1_3.x= v3.x-v1.x;
//     e_1_3.y = v3.y-v1.y;
//     e_1_3.z = v3.z-v1.z;

//     Vertex e_2_3;
//     e_2_3.x= v3.x-v2.x;
//     e_2_3.y = v3.y-v2.y;
//     e_2_3.z = v3.z-v2.z;

//     Vertex e_3_4;
//     e_3_4.x= v4.x-v3.x;
//     e_3_4.y = v4.y-v3.y;
//     e_3_4.z = v4.z-v3.z;
//     double phi13= phi3-phi1;

//     double a= sq(phi13)*(e_1_3.x*e_1_3.x + e_1_3.y*e_1_3.y + e_1_3.z*e_1_3.z)- sq(e_1_3.x*e_1_3.x+ e_1_3.y*e_1_3.y+ e_1_3.z*e_1_3.z);
//     double b = sq(phi13)*(e_2_3.x* e_2_3.x+ e_2_3.y*e_2_3.y+ e_2_3.z*e_2_3.z) -sq(e_2_3.x*e_1_3.x + e_2_3.y*e_1_3.y + e_2_3.z*e_1_3.z);
//     double c= sq(phi13)*(2*(e_1_3.x*e_2_3.x+ e_1_3.y*e_2_3.y+ e_1_3.z*e_2_3.z)) - 2*((e_1_3.x*e_1_3.x+ e_1_3.y*e_1_3.y+ e_1_3.z*e_1_3.z)*(e_2_3.x*e_1_3.x + e_2_3.y*e_1_3.y + e_2_3.z*e_1_3.z));
//     double d= sq(phi13)*(e_1_3.x* e_3_4.x+ e_1_3.y*e_3_4.y+ e_1_3.z*e_3_4.z) - 2*((e_3_4.x*e_1_3.x+ e_3_4.y*e_1_3.y+ e_3_4.z*e_1_3.z)*(e_1_3.x*e_1_3.x+ e_1_3.y*e_1_3.y+ e_1_3.z*e_1_3.z));
//     double e= sq(phi13)*(e_2_3.x* e_3_4.x+ e_2_3.y*e_3_4.y+ e_2_3.z*e_3_4.z) - 2*((e_2_3.x*e_1_3.x+ e_2_3.y*e_1_3.y+ e_2_3.z*e_1_3.z)*(e_3_4.x*e_1_3.x+ e_3_4.y*e_1_3.y+ e_3_4.z*e_1_3.z));
//     double f= sq(phi13)*(e_3_4.x*e_3_4.x+ e_3_4.y*e_3_4.y+ e_3_4.z*e_3_4.z) - sq(e_3_4.x*e_1_3.x+ e_3_4.y*e_1_3.y+ e_3_4.z*e_1_3.z);
//     return 2*a*x + c*y+d;
// }

// double df2_dx(double x, double y, Vertex v1, Vertex v2, Vertex v3, Vertex v4, double phi1, double phi2, double phi3){
//     Vertex e_1_3; 
//     e_1_3.x= v3.x-v1.x;
//     e_1_3.y = v3.y-v1.y;
//     e_1_3.z = v3.z-v1.z;

//     Vertex e_2_3;
//     e_2_3.x= v3.x-v2.x;
//     e_2_3.y = v3.y-v2.y;
//     e_2_3.z = v3.z-v2.z;

//     Vertex e_3_4;
//     e_3_4.x= v4.x-v3.x;
//     e_3_4.y = v4.y-v3.y;
//     e_3_4.z = v4.z-v3.z;

//     double f= (phi3-phi2)*(e_1_3.x*e_1_3.x+ e_1_3.y*e_1_3.y+ e_1_3.z*e_1_3.z) - (phi3-phi1)*(e_1_3.x*e_2_3.x+ e_1_3.y*e_2_3.y+ e_1_3.z*e_2_3.z);

//     return f;
// }

// double df1_dy(double x, double y, Vertex v1, Vertex v2, Vertex v3, Vertex v4, double phi1, double phi2, double phi3){
//       Vertex e_1_3; 
//     e_1_3.x= v3.x-v1.x;
//     e_1_3.y = v3.y-v1.y;
//     e_1_3.z = v3.z-v1.z;

//     Vertex e_2_3;
//     e_2_3.x= v3.x-v2.x;
//     e_2_3.y = v3.y-v2.y;
//     e_2_3.z = v3.z-v2.z;

//     Vertex e_3_4;
//     e_3_4.x= v4.x-v3.x;
//     e_3_4.y = v4.y-v3.y;
//     e_3_4.z = v4.z-v3.z;
//     double phi13= phi3-phi1;

//     double a= sq(phi13)*(e_1_3.x*e_1_3.x + e_1_3.y*e_1_3.y + e_1_3.z*e_1_3.z)- sq(e_1_3.x*e_1_3.x+ e_1_3.y*e_1_3.y+ e_1_3.z*e_1_3.z);
//     double b = sq(phi13)*(e_2_3.x* e_2_3.x+ e_2_3.y*e_2_3.y+ e_2_3.z*e_2_3.z) -sq(e_2_3.x*e_1_3.x + e_2_3.y*e_1_3.y + e_2_3.z*e_1_3.z);
//     double c= sq(phi13)*(2*(e_1_3.x*e_2_3.x+ e_1_3.y*e_2_3.y+ e_1_3.z*e_2_3.z)) - 2*((e_1_3.x*e_1_3.x+ e_1_3.y*e_1_3.y+ e_1_3.z*e_1_3.z)*(e_2_3.x*e_1_3.x + e_2_3.y*e_1_3.y + e_2_3.z*e_1_3.z));
//     double d= sq(phi13)*(e_1_3.x* e_3_4.x+ e_1_3.y*e_3_4.y+ e_1_3.z*e_3_4.z) - 2*((e_3_4.x*e_1_3.x+ e_3_4.y*e_1_3.y+ e_3_4.z*e_1_3.z)*(e_1_3.x*e_1_3.x+ e_1_3.y*e_1_3.y+ e_1_3.z*e_1_3.z));
//     double e= sq(phi13)*(e_2_3.x* e_3_4.x+ e_2_3.y*e_3_4.y+ e_2_3.z*e_3_4.z) - 2*((e_2_3.x*e_1_3.x+ e_2_3.y*e_1_3.y+ e_2_3.z*e_1_3.z)*(e_3_4.x*e_1_3.x+ e_3_4.y*e_1_3.y+ e_3_4.z*e_1_3.z));
//     double f= sq(phi13)*(e_3_4.x*e_3_4.x+ e_3_4.y*e_3_4.y+ e_3_4.z*e_3_4.z) - sq(e_3_4.x*e_1_3.x+ e_3_4.y*e_1_3.y+ e_3_4.z*e_1_3.z);
//     return 2*b*y + c*x+e;
// }

// double df2_dy(double x, double y, Vertex v1, Vertex v2, Vertex v3, Vertex v4, double phi1, double phi2, double phi3){
//     Vertex e_1_3; 
//     e_1_3.x= v3.x-v1.x;
//     e_1_3.y = v3.y-v1.y;
//     e_1_3.z = v3.z-v1.z;

//     Vertex e_2_3;
//     e_2_3.x= v3.x-v2.x;
//     e_2_3.y = v3.y-v2.y;
//     e_2_3.z = v3.z-v2.z;

//     Vertex e_3_4;
//     e_3_4.x= v4.x-v3.x;
//     e_3_4.y = v4.y-v3.y;
//     e_3_4.z = v4.z-v3.z;

//     double f= (phi3-phi2)*(e_2_3.x*e_1_3.x+ e_2_3.y*e_1_3.y+ e_2_3.z*e_1_3.z) - (phi3-phi1)*(e_2_3.x*e_2_3.x+ e_2_3.y*e_2_3.y+ e_2_3.z*e_2_3.z);

//     return f;

// }
double f1(double x, double y, Vertex v1, Vertex v2, Vertex v3, Vertex v4, double phi1, double phi2, double phi3){
    Vertex e13;
    e13.x= v3.x-v1.x;
    e13.y= v3.y-v1.y;
    e13.z= v3.z-v1.z;
    Vertex e23;
    e23.x= v3.x-v2.x;
    e23.y= v3.y-v2.y;
    e23.z= v3.z-v2.z;
    double a= (phi3-phi2)*sq(norm(e13))- (phi3-phi1)*((v3.x-v2.x)*(v3.x-v1.x)+ (v3.y-v2.y)*(v3.y-v1.y)+ (v3.z-v2.z)*(v3.z-v1.z));
    double b= (phi3-phi2)*((v3.x-v2.x)*(v3.x-v1.x)+ (v3.y-v2.y)*(v3.y-v1.y)+ (v3.z-v2.z)*(v3.z-v1.z))- (phi3-phi1)*sq(norm(e23));
    double c= (phi3-phi2)*((v3.x-v1.x)*(v3.x+v4.x)+ (v3.y-v1.y)*(v3.y+v4.y)+ (v3.z-v1.z)*(v3.z+v4.z))- (phi3-phi1)*((v3.x-v2.x)*(v3.x+v4.x)+ (v3.y-v2.y)*(v3.y+v4.y)+ (v3.z-v2.z)*(v3.z+v4.z));
    return a*x+b*y+c;
}
double f2(double x, double y, Vertex v1, Vertex v2, Vertex v3, Vertex v4, double phi1, double phi2, double phi3){
    Vertex e54;
    e54.x= v4.x-x*v1.x-y*v2.x-(1-x-y)*v3.x;
    e54.y= v4.y-x*v1.y-y*v2.y-(1-x-y)*v3.y;
    e54.z= v4.z-x*v1.z-y*v2.z-(1-x-y)*v3.z;
    Vertex e13;
    e13.x= v3.x-v1.x;
    e13.y= v3.y-v1.y;
    e13.z= v3.z-v1.z;
    double a = sq(norm(e13));
    double b= (v3.x-v1.x)*(v3.x-v2.x)+ (v3.y-v1.y)*(v3.y-v2.y)+ (v3.z-v1.z)*(v3.z-v2.z);
    double c= (v3.x- v1.x)*(v3.x+v4.x)+ (v3.y-v1.y)*(v3.y+v4.y)+ (v3.z-v1.z)*(v3.z+v4.z);
    return (phi3-phi1)*norm(e54)-a*x-b*y-c;
}
// double df1_dx(double x, double y, Vertex v1, Vertex v2, Vertex v3, Vertex v4, double phi1, double phi2, double phi3){
//     Vertex e13;
//     e13.x= v3.x-v1.x;
//     e13.y= v3.y-v1.y;
//     e13.z= v3.z-v1.z;
//     double a= (phi3-phi2)*sq(norm(e13))- (phi3-phi1)*((v3.x-v2.x)*(v3.x-v1.x)+ (v3.y-v2.y)*(v3.y-v1.y)+ (v3.z-v2.z)*(v3.z-v1.z));
//     return a; 
// }
// double df1_dy(double x, double y, Vertex v1, Vertex v2, Vertex v3, Vertex v4, double phi1, double phi2, double phi3){
//     Vertex e23;
//     e23.x= v3.x-v2.x;
//     e23.y= v3.y-v2.y;
//     e23.z= v3.z-v2.z;
//     double b= (phi3-phi2)*((v3.x-v2.x)*(v3.x-v1.x)+ (v3.y-v2.y)*(v3.y-v1.y)+ (v3.z-v2.z)*(v3.z-v1.z))- (phi3-phi1)*sq(norm(e23));
//     return b;
// }
// double df2_dx(double x, double y, Vertex v1, Vertex v2, Vertex v3, Vertex v4, double phi1, double phi2, double phi3){
//     Vertex e54;
//     e54.x= v4.x-x*v1.x-y*v2.x-(1-x-y)*v3.x;
//     e54.y= v4.y-x*v1.y-y*v2.y-(1-x-y)*v3.y;
//     e54.z= v4.z-x*v1.z-y*v2.z-(1-x-y)*v3.z;
//     Vertex e13;
//     e13.x= v3.x-v1.x;
//     e13.y= v3.y-v1.y;
//     e13.z= v3.z-v1.z;
//     double a = sq(norm(e13));
//     return (phi3-phi1)*((e54.x*e13.x)+ (e54.y*e13.y)+ (e54.z*e13.z))/(norm(e54))-a;
// }
// double df2_dy(double x, double y, Vertex v1, Vertex v2, Vertex v3, Vertex v4, double phi1, double phi2, double phi3){
//     Vertex e54;
//     e54.x= v4.x-x*v1.x-y*v2.x-(1-x-y)*v3.x;
//     e54.y= v4.y-x*v1.y-y*v2.y-(1-x-y)*v3.y;
//     e54.z= v4.z-x*v1.z-y*v2.z-(1-x-y)*v3.z;
//     Vertex e23;
//     e23.x= v3.x-v2.x;
//     e23.y= v3.y-v2.y;
//     e23.z= v3.z-v2.z;
//     double b= (v3.x-v1.x)*(v3.x-v2.x)+ (v3.y-v1.y)*(v3.y-v2.y)+ (v3.z-v1.z)*(v3.z-v2.z);
//     return (phi3-phi1)*((e54.x*e23.x)+ (e54.y*e23.y)+ (e54.z*e23.z))/(norm(e54)) -b;
// }

//solving for lambda 1 and lambda 2 by Broyden's method method:
vector<double> Broyden_method(Vertex v1, Vertex v2, Vertex v3, Vertex v4, double phi1, double phi2, double phi3, int* flag){
    //initial iteration
    double x=0.1;
    double y=0.5;
    int max_iterations=50;
    double tolerance= 1e-8;

    double inv_J11=1;
    double inv_J12=0;
    double inv_J21=0;
    double inv_J22=1;

    double J11_inv_new;
    double J12_inv_new;
    double J21_inv_new;
    double J22_inv_new;

    double x_new, y_new;
    #pragma acc parallel loop
    for (int i=0; i<max_iterations; i++){
        double F1=f1(x, y, v1, v2, v3, v4, phi1, phi2, phi3);
        double F2= f2(x, y, v1, v2, v3, v4, phi1, phi2, phi3);

        x_new= x-(inv_J11*F1+inv_J12*F2);
        y_new= y-(inv_J21*F1+ inv_J22*F2);

        double delta_x= x_new-x;
        double delta_y= y_new-y;

        double delta_F1= f1(x_new, y_new, v1, v2, v3, v4, phi1, phi2, phi3)- F1;
        double delta_F2= f2(x_new, y_new, v1, v2, v3, v4, phi1, phi2, phi3)- F2;

        if ((abs(delta_x)<tolerance && abs(delta_y)<tolerance)){
            *flag=1;
            return {x_new, y_new};
        }
        //Good Broyden's method: 
        J11_inv_new= inv_J11+((delta_x- inv_J11*delta_F1- inv_J12* delta_F2)*(inv_J11*delta_x + inv_J21*delta_y))/((delta_x*inv_J11 + delta_y*inv_J21)*delta_F1 +(delta_x*inv_J12 +delta_y*inv_J22)*delta_F2);
        J12_inv_new= inv_J12+((delta_x- inv_J11*delta_F1-inv_J12* delta_F2)*(inv_J12*delta_x + inv_J22*delta_y))/((delta_x*inv_J11 + delta_y*inv_J21)*delta_F1 +(delta_x*inv_J12 + delta_y*inv_J22)*delta_F2);
        J21_inv_new= inv_J21+((delta_y- inv_J21*delta_F1- inv_J22* delta_F2)*(inv_J11*delta_x +inv_J21*delta_y))/((delta_x*inv_J11 +delta_y*inv_J21)*delta_F1+(delta_x*inv_J12 +delta_y*inv_J22)*delta_F2);
        J22_inv_new= inv_J22+((delta_y-inv_J21*delta_F1- inv_J22* delta_F2)*(inv_J12*delta_x + inv_J22*delta_y))/((delta_x*inv_J11 + delta_y*inv_J21)*delta_F1 +(delta_x*inv_J12 +delta_y*inv_J22)*delta_F2);
        
        // // Bad broyden's method:
        // double denominator;
        // denominator= pow(delta_F1,2) + pow(delta_F2,2);
        // J11_inv_new= inv_J11 + (delta_F1*(delta_x- inv_J11*delta_F1- inv_J12* delta_F2))/denominator;
        // J12_inv_new= inv_J12 + (delta_F2*(delta_x- inv_J11*delta_F1- inv_J12* delta_F2))/denominator;
        // J21_inv_new= inv_J21 + (delta_F1*(delta_y- inv_J21*delta_F1- inv_J22* delta_F2))/denominator;
        // J22_inv_new= inv_J22 + (delta_F2*(delta_y- inv_J21*delta_F1- inv_J22* delta_F2))/denominator;

        inv_J11= J11_inv_new;
        inv_J12= J12_inv_new;
        inv_J21=J21_inv_new;
        inv_J22= J22_inv_new;

        x= x_new;
        y= y_new; 
    }
    *flag=0;
    return {NAN, NAN};
}

// 2D local solver
double localsolver_2D(Vertex v1, Vertex v2, Vertex v3, double phi1, double phi2){
    double tol = 1e-12;
    double phi= phi2-phi1;
    Vertex e13;
    e13.x= v3.x- v1.x;
    e13.y=v3.y- v1.y;
    e13.z= v3.z- v1.z;

    Vertex e12;
    e12.x = v2.x- v1.x;
    e12.y= v2.y- v1.y;
    e12.z= v2.z- v1.z;

    double a1= e13.x;
    double b1= -e12.x;
    double a2=e13.y;
    double b2= -e12.y;
    double a3= e13.z;
    double b3= -e12.z;

    Vertex e32;
    e32.x= v3.x-v2.x;
    e32.y= v3.y-v2.y;
    e32.z= v3.z-v2.z;

    if (phi1 == INFINITY && phi2 < INFINITY) {
        double phi3 = phi2 + norm(e32);
        return phi3;
    }
    else if (phi1 < INFINITY && phi2 == INFINITY) {
        double phi3 = phi1 + norm(e13);
        return phi3;
    }
    else if (phi1 == INFINITY && phi2 == INFINITY) {
        return INFINITY;
    } 
    else {
        double A= (sq(phi)-sq(b1))*sq(b1) + (sq(phi)-sq(b2))*sq(b2) + (sq(phi)- sq(b3))*sq(b3) - 2*sq(b1)*sq(b2) - 2*sq(b1)*sq(b3) - 2*sq(b2)*sq(b3);
        double B= (sq(phi)- sq(b1))*2*a1*b1 + (sq(phi)- sq(b2))*2*a2*b2 + (sq(phi)-sq(b3))*2*a3*b3 - sq(b3)*2*a2*b2 - sq(b2)*2*a3*b3 - sq(b1)*2*a3*b3 - sq(b3)*2*a1*b1 - sq(b2)*2*a1*b1 - sq(b1)*2*a2*b2;
        double C= (sq(phi)- sq(b1))*sq(a1) + (sq(phi)- sq(b2))*sq(a2) + (sq(phi)- sq(b3))*sq(a3) - 2*a1*a2*b1*b2 - 2*a2*a3*b2*b3 - 2*a1*a3*b1*b3;
        if (abs(A) < tol) {
            if (abs(B) < tol) {
                double phi3 = min(phi1 + norm(e13), phi2 + norm(e32));
                return phi3;
            } 
            else {
                double lambda1 = -C/B;
                Vertex t;
                t.x= v3.x- v1.x -lambda1*(v2.x-v1.x);
                t.y= v3.y- v1.y -lambda1*(v2.y-v1.y);
                t.z= v3.z- v1.z- lambda1*(v2.z-v1.z);
                double phi3 =lambda1 * phi + phi1 + norm(t);
                return phi3;
            }
        }
        if (B*B-4*A*C<0){
            double phi3 = min(phi1 + norm(e13), phi2 + norm(e32));
            return phi3;
        }
        else {
            double lambda1 = (-B+sqrt(B*B-4*A*C))/(2*A);
            double lambda2 = (-B-sqrt(B*B-4*A*C))/(2*A);

            if (lambda1 > 0 && lambda1<1 && lambda2 > 0 && lambda2<1) {
                Vertex t1;
                t1.x= v3.x- v1.x -lambda1*(v2.x-v1.x);
                t1.y= v3.y-v1.y - lambda1*(v2.y-v1.y);
                t1.z= v3.z- v1.z - lambda1*(v2.z-v1.z);

                Vertex t2;
                t2.x= v3.x- v1.x - lambda2*(v2.x-v1.x);
                t2.y= v3.y- v1.y - lambda2*(v2.y-v1.y);
                t2.z= v3.z- v1.z -lambda2*(v2.z-v1.z);

                double phi3 =min(lambda1 * phi + phi1 + norm(t1), lambda2 * phi + phi1 + norm(t2));
                return phi3;
            }
            else if (lambda1 > 0 && lambda1<1) {
                Vertex t; 
                t.x= v3.x- v1.x -lambda1*(v2.x-v1.x);
                t.y= v3.y-v1.y - lambda1*(v2.y-v1.y);
                t.z= v3.z- v1.z - lambda1*(v2.z-v1.z);

                double phi3 = lambda1 * phi + phi1 + norm(t);
                return phi3;
            }
            else if (lambda2 > 0 && lambda2<1) {
                Vertex t; 
                t.x= v3.x- v1.x - lambda2*(v2.x-v1.x);
                t.y= v3.y - v1.y - lambda2*(v2.y-v1.y);
                t.z= v3.z- v1.z -lambda2*(v2.z-v1.z);

                double phi3 =lambda2 * phi +phi1 + norm(t);
                return phi3;
            }
            else {
                double phi3 =min(phi1 + norm(e13), phi2 + norm(e32));
                return phi3;
            }
        }
    }

}
// VectorXd equations(Vector2d& variables,  Vertex v1,  Vertex v2,  Vertex v3,  Vertex v4, double phi1, double phi2, double phi3) {
//     double x = variables[0];
//     double y =variables[1];
//     VectorXd result(2);
//     result[0] =f1(x, y, v1, v2, v3, v4, phi1, phi2, phi3);
//     result[1] = f2(x, y, v1, v2, v3, v4, phi1, phi2, phi3);
//     return result;
// }
// MatrixXd jacobian(Vector2d& variables, Vertex v1, Vertex v2, Vertex v3, Vertex v4, double phi1, double phi2, double phi3) {
//     double x = variables[0];
//     double y = variables[1];

//     MatrixXd jacobian(2, 2);
//     jacobian(0, 0)= df1_dx(x, y, v1, v2, v3, v4, phi1, phi2, phi3);
//     jacobian(0, 1)= df1_dy(x, y, v1, v2, v3, v4, phi1, phi2, phi3);
//     jacobian(1, 0)=df2_dx(x, y, v1, v2, v3, v4, phi1, phi2, phi3);
//     jacobian(1, 1)= df2_dy(x, y, v1, v2, v3, v4, phi1, phi2, phi3);

//     return jacobian;
// }
// vector<double> Newton_method(const Vertex v1, const Vertex v2, const Vertex v3, const Vertex v4, double phi1, double phi2, double phi3, int *flag){
//     Vector2d variables(0.9, 0.8);
//     double tolerance = 1e-4;
//     int max_iterations = 100;
//     *flag=0;

//     for (int i=0; i<max_iterations; i++){
//         VectorXd F= equations(variables, v1, v2, v3, v4, phi1, phi2, phi3);
//         MatrixXd J= jacobian(variables, v1, v2, v3, v4, phi1, phi2, phi3);
//         Vector2d change= J.fullPivLu().solve(-F);
//         variables += change;
//         if (abs(variables[0])<tolerance){ 
//             variables(0)=0.0;
//         }
//         if (abs(variables[1])<tolerance){
//             variables(1)=0.0;
//         }
//         if (change.norm() <tolerance){
//             *flag=1;
//             return {variables[0], variables[1]};
//         }
//     }
//     // return nan
//     *flag=0;
//     return {NAN, NAN};   
// }
// 3Dlocal solver
double localsolver(Vertex v1, Vertex v2, Vertex v3, Vertex v4, double phi1, double phi2, double phi3){
    double phi4;
    if (phi1== INFINITY && phi2== INFINITY && phi3== INFINITY){
        return INFINITY;
    }
    else if(phi1== INFINITY && phi2== INFINITY && phi3 <INFINITY){
        // double phia= localsolver_2D(v2, v3, v4, phi2, phi3);
        // double phib= localsolver_2D(v1, v3, v4, phi1, phi3);
        // phi4= min(phia, phib);
        Vertex t;
        t.x= v4.x- v3.x;
        t.y= v4.y- v3.y;
        t.z= v4.z- v3.z;
        phi4= phi3 + norm(t);  
        //cout<<"read 1"<<endl;
        return phi4;
    }
    else if(phi1 ==INFINITY && phi3== INFINITY && phi2 <INFINITY){
        // double phia= localsolver_2D(v2, v3, v4, phi2, phi3);
        // double phib= localsolver_2D(v1, v2, v4, phi1, phi2);
        // phi4= min(phia, phib);
        Vertex t;
        t.x= v4.x- v2.x;
        t.y= v4.y- v2.y;
        t.z= v4.z- v2.z;
        phi4= phi2 + norm(t);
        //cout<<"read 2"<<endl;
        return phi4;
    }
    else if(phi2== INFINITY && phi3== INFINITY && phi1<INFINITY){
        // double phia= localsolver_2D(v1, v3, v4, phi1, phi3);
        // double phib= localsolver_2D(v1, v2, v4, phi1, phi2);
        // phi4= min(phia, phib);
        Vertex t;
        t.x= v4.x- v1.x;
        t.y= v4.y- v1.y;
        t.z= v4.z- v1.z;
        phi4= phi1 + norm(t);
        //cout<<"read 3"<<endl;
        return phi4;
       
    }
    else if(phi1< INFINITY && phi2< INFINITY && phi3== INFINITY){
        phi4= localsolver_2D(v1, v2, v4, phi1, phi2);
        return phi4;
    }
    else if(phi1 < INFINITY && phi3< INFINITY && phi2==INFINITY){
        phi4= localsolver_2D(v1, v3, v4, phi1, phi3);
        return phi4;
    }
    else if (phi2 < INFINITY && phi3 < INFINITY && phi1==INFINITY){
        phi4= localsolver_2D(v2, v3, v4, phi2, phi3);
        return phi4;
    }
    else{
    int *flag= new int;
    *flag=0;
    vector<double> lambdas = Broyden_method(v1, v2, v3, v4, phi1, phi2, phi3, flag);
    double lambda1= lambdas[0];
    double lambda2= lambdas[1];
    if (((lambda1>0 && lambda1<1) && (lambda2>0 && lambda2<1)) && *flag==1){
        Vertex v5;
        v5.x= lambda1*v1.x+ lambda2*v2.x+ (1-lambda1-lambda2)*v3.x;
        v5.y= lambda1*v1.y+ lambda2*v2.y+ (1-lambda1-lambda2)*v3.y;
        v5.z= lambda1*v1.z+ lambda2*v2.z+ (1-lambda1-lambda2)*v3.z;
        Vertex e_5_4;
        e_5_4.x= v4.x-v5.x;
        e_5_4.y= v4.y-v5.y;
        e_5_4.z= v4.z-v5.z;

        phi4= lambda1*phi1+ lambda2*phi2 + (1-lambda1-lambda2)*phi3+ norm(e_5_4);
        return phi4; 
    }
    else if((lambda1==0 && lambda2==1) && *flag==1){
        Vertex v5;
        v5.x= lambda2*v2.x+ (1-lambda2)*v3.x;
        v5.y= lambda2*v2.y+ (1-lambda2)*v3.y;
        v5.z= lambda2*v2.z+ (1-lambda2)*v3.z;
        Vertex e_5_4;
        e_5_4.x= v4.x-v5.x;
        e_5_4.y= v4.y-v5.y;
        e_5_4.z= v4.z-v5.z;

        phi4= phi2 + norm(e_5_4);
        return phi4;
    }
    else if((lambda1==1 && lambda2==0) && *flag==1){
        Vertex v5;
        v5.x= lambda1*v1.x+ (1-lambda1)*v3.x;
        v5.y= lambda1*v1.y+ (1-lambda1)*v3.y;
        v5.z= lambda1*v1.z+ (1-lambda1)*v3.z;
        Vertex e_5_4;
        e_5_4.x= v4.x-v5.x;
        e_5_4.y= v4.y-v5.y;
        e_5_4.z= v4.z-v5.z;
        phi4= phi1+ norm(e_5_4);
        return phi4;
    }
    else if((lambda1==0 && lambda2==0) && *flag==1){
        Vertex v5;
        v5.x= v3.x;
        v5.y= v3.y;
        v5.z= v3.z;
        Vertex e_5_4;
        e_5_4.x= v4.x-v5.x;
        e_5_4.y= v4.y-v5.y;
        e_5_4.z= v4.z-v5.z;
        phi4 = phi3+ norm(e_5_4);
        return phi4;
    }
    else{
        //cout<<"Exceptional case"<<endl;
        //consider 3 triangles namely tri_124, tri_134, tri_234
        vector<double> values_2d;
        // solve for each triangle and push the value in vector. 
        values_2d.push_back(localsolver_2D(v1, v2, v4, phi1, phi2));
        values_2d.push_back(localsolver_2D(v1,v3, v4, phi1, phi3));
        values_2d.push_back(localsolver_2D(v2, v3, v4, phi2, phi3));
        phi4= *min_element(values_2d.begin(), values_2d.end());
        return phi4;
    }
    }
}

// bool checklie(int index1, int index2, int index3,int index4 , vector<Tetrahedron> &tetrahedra){
//     set<int> s;
//     s={index1, index2, index3, index4};
//     for (Tetrahedron tet: tetrahedra){
//         set<int> t;
//         t={tet.v1, tet.v2, tet.v3, tet.v4};
//         if (s==t){
//             return true;
//         }
//         else continue;
//     }
//     return false;
// }

//update function
double update(int v, vector<double> &values, vector<Vertex> &vertices, vector<Tetrahedron> &tetrahedra){
    vector<int> oner = onering1(tetrahedra, v);

    vector<double> valueslist;
    #pragma acc parallel loop
    for (int i=0; i< oner.size(); i=i+4){
        int indexv1= oner[i];
        int indexv2 = oner[(i+1)%oner.size()];
        int indexv3 = oner[(i+2)%oner.size()];
        int indexv4 = oner[(i+3)%oner.size()];
        //bool check= checklie(indexv1, indexv2, indexv3, indexv4, tetrahedra);
            if (v== indexv1){
                double phi1= values[indexv2];
                double phi2= values[indexv3];
                double phi3= values[indexv4];
                Vertex v1= vertices[indexv2];
                Vertex v2= vertices[indexv3];
                Vertex v3= vertices[indexv4];
                Vertex v4= vertices[v];
                double val;
                val= localsolver(v1, v2, v3, v4, phi1, phi2, phi3);
                valueslist.push_back(val);
            }
            else if(v== indexv2){
                double phi1= values[indexv1];
                double phi2= values[indexv3];
                double phi3= values[indexv4];
                Vertex v1= vertices[indexv1];
                Vertex v2= vertices[indexv3];
                Vertex v3= vertices[indexv4];
                Vertex v4= vertices[v];
                double val;
                val= localsolver(v1, v2, v3, v4, phi1, phi2, phi3);
                valueslist.push_back(val);
            }
            else if(v== indexv3){
                double phi1= values[indexv1];
                double phi2= values[indexv2];
                double phi3= values[indexv4];
                Vertex v1= vertices[indexv1];
                Vertex v2= vertices[indexv2];
                Vertex v3= vertices[indexv4];
                Vertex v4= vertices[v];
                double val;
                val= localsolver(v1, v2, v3, v4, phi1, phi2, phi3);
                valueslist.push_back(val);
            }
            else {
                double phi1= values[indexv1];
                double phi2= values[indexv2];
                double phi3= values[indexv3];
                Vertex v1= vertices[indexv1];
                Vertex v2= vertices[indexv2];
                Vertex v3= vertices[indexv3];
                Vertex v4= vertices[v];
                double val;
                val= localsolver(v1, v2, v3, v4, phi1, phi2, phi3);
                valueslist.push_back(val);
            }
    }
    double phi4= *min_element(valueslist.begin(), valueslist.end());
    return phi4;
}

bool IsvinL(int v, vector<int>L){
    for(int i=0; i<L.size(); i++){
        if (v== L[i]){
            return true;
        }
        else continue;
    }
    return false;
}

int main(){
    vector<Vertex> vertices;
    vector<Tetrahedron> tetrahedra;

    // string filename1 = "vertices2.txt";
    // string filename2 = "tets2.txt";
    // if (readvertices(filename1, vertices)) {
    //     cout<< "Sucessfully loaded"<<endl;
    // } else {
    //     cout <<"Failed to load the vertices file." << endl;
    // }
    // if (readtets(filename2, tetrahedra)) {
    //     cout<< "Sucessfully loaded"<<endl;
    // } else {
    //     cout <<"Failed to load the tetrahedra file." << endl;
    // }

    string file1= "sphere.vtk";
    readVTKfile(file1, vertices, tetrahedra);

    cout<<vertices.size()<<endl;
    cout<<tetrahedra.size()<<endl;
    //set B
    vector<int> B;
    B={2,3,4,5,6,7,8};
    //check which vertices is in B
    // for (int i=0; i<vertices.size(); i++){
    //     if (IsvinL(i,B)){
    //         cout<<i<<" ";
    //     }
    // }
    // Active list
    vector<int> L;

    //Initialise values
    vector<double> values;
    for (int i=0; i<vertices.size(); i++){
        if (IsvinL(i,B)){
            values.push_back(0);
        }
        else values.push_back(INFINITY);
    }

    //Initialise active list
    for(int i=0; i<vertices.size(); i++){
        if (!IsvinL(i,B)){
            if((onering2(tetrahedra, i)).size()!=0){
                for(int j: onering2(tetrahedra, i)){
                    if(IsvinL(j,B)){
                        if(!IsvinL(i,L)){
                            L.push_back(i);
                        }
                    }
                }
            }
        }
    }

    cout<<"size of L is :"<< L.size()<<endl;
    double tolerance= 1e-12; 
    clock_t start, end;
    start =clock();
    //#pragma acc kernels
    int count =0;
    while(!L.empty()){
        int ls_old= L.size();
        for (int i=0; i<vertices.size(); i++){
            //cout<<"started vertex: "<< i<<endl;
            if (!IsvinL(i,B)){
                double p =values[i];
                double q= update(i, values, vertices, tetrahedra);
                if (abs(p-q)<tolerance){
                    for(int j: onering2(tetrahedra, i)){
                        if (!IsvinL(j,L)){
                            p=values[j];
                            q= update(j, values, vertices, tetrahedra);
                            if (p>q){
                                values[j]=q;
                                L.push_back(j);
                            }
                        }
                    }
                    if (IsvinL(i,L)){
                        L.erase(remove(L.begin(), L.end(), i), L.end());
                    }
                }
                else {
                    values[i]= min(p,q);
                }
            }
        }
        int ls_new= L.size();
        // print values
        // for (int i=0; i<values.size(); i++){
        //     cout<< values[i]<<" ";
        // }
        // cout<<endl;
        //print L
        // for(int i=0; i<L.size(); i++){
        //     cout<< L[i]<<" ";
        // }
        // cout<<endl;
        cout<<"size of L is :"<< L.size()<<endl;
        int diff= ls_old-ls_new;
        if(diff==0){
            count++;
        }
        if (count>=5){
            break;
        }
    }
    end= clock();

    double time_taken= double (end-start)/ double(CLOCKS_PER_SEC);
    cout<< "Time taken by the cpu is "<< time_taken<<endl;
    ofstream outputfile1("valuescalc.txt");
    if (outputfile1.is_open()){
        for (int i=0; i<values.size(); i++){
            outputfile1<<values[i]<<endl;
        }
    }
    outputfile1.close();

    vector<double> values1;
    for (int i=0; i<vertices.size(); i++){
        vector<double> err;
        for (int j=0; j<B.size(); j++){
            Vertex diff;
            diff.x= vertices[i].x- vertices[B[j]].x;
            diff.y= vertices[i].y- vertices[B[j]].y;
            diff.z= vertices[i].z- vertices[B[j]].z;
            err.push_back(norm(diff));
        }
        values1.push_back(*min_element(err.begin(), err.end()));
    }
    ofstream outputfile2("valuestrue.txt");
    if (outputfile2.is_open()){
        for (int i=0; i<values1.size(); i++){
            outputfile2<<values1[i]<<endl;
        }
    }
    outputfile2.close();
    vector<double> errors;
    for (int i=0; i<values.size(); i++){
        errors.push_back(abs(values[i]-values1[i]));
    }
    ofstream outputfile3("errors.txt");
    if (outputfile3.is_open()){
        for (int i=0; i<values.size(); i++){
            outputfile3<<errors[i]<<endl;
        }
    }
    outputfile3.close();
    for (int i=0; i<vertices.size(); i++){
        cout<< values[i]<<" ";
    }
    cout<<endl;
    
    cout<<"errors are: "<<endl;
    double sums=0;
    for (int i=0; i<errors.size(); i++){
        cout<< errors[i]<<" ";
        sums =sums + errors[i]*errors[i];
    }
    cout<< "max of error is :"<<endl;
    cout << *max_element(errors.begin(), errors.end())<<endl;
    cout<< "norm of the error is: "<<endl;
    cout<< sqrt(sums)/vertices.size()<<endl;
    cout<< "Time taken by the cpu is "<< time_taken<<endl;
    vector<double> relerror;
    for (int i=0; i<vertices.size(); i++){
        if (values1[i]==0){
            relerror.push_back(0.0);
        }
        else {
            double rat=errors[i]/values1[i];
            relerror.push_back(rat);
        }
    }
    double sum3=0; 
    for (int i=0; i<relerror.size(); i++){
        sum3+=relerror[i]*relerror[i];
    }
    sum3= sqrt(sum3);
    cout<< "relative error is coming out to be: "<<" for the "<<file1<<" "<<sum3/vertices.size()<<endl;
    //print one ring of every vertex
    // for (int i=0; i<vertices.size(); i++){
    //     cout<<"vertex: "<<i<<endl;
    //     for (int j: onering(tetrahedra, i)){
    //         cout<<j<<" ";
    //     }
    //     cout<<endl;
    // }
    // for(int i=0; i<vertices.size(); i++){
    //     vector<int> oner= onering(tetrahedra, i);
    //     for(int j=0; j<oner.size(); j=j+4){
    //         int indexv1= oner[j];
    //         int indexv2= oner[(j+1)%oner.size()];
    //         int indexv3= oner[(j+2)%oner.size()];
    //         int indexv4= oner[(j+3)%oner.size()];
    //         // set<int> s;
    //         // s={indexv1, indexv2, indexv3, indexv4};
    //         bool check= checklie(indexv1, indexv2, indexv3, indexv4, tetrahedra);
    //         // for (Tetrahedron tet: tetrahedra){
    //         //     set<int> t;
    //         //     t={tet.v1, tet.v2, tet.v3, tet.v4};
    //         //     if (s==t){
    //         //         check=true;
    //         //         break;
    //         //     }
    //         //     else continue;
    //         // }
    //         if (check==true){
    //            cout<<"good vertex: "<<i<<endl;
    //            cout<<"good triangle: "<<indexv1<<" "<<indexv2<<" "<<indexv3<<" "<<indexv4<<endl;
    //         }
    //         else continue;
    //     }
    // }

}
