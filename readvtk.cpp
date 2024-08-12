#include <iostream> 
#include <fstream>
#include <sstream>
#include <vector> 
using namespace std;

struct Vertex{
    double x,y,z;
};
struct Tetrahedron{
    int v1, v2, v3, v4;
};
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
int main() {
    string filename= "cone.vtk";
    vector<Vertex> vertices;
    vector<Tetrahedron> tetrahedra;
    readVTKfile(filename, vertices, tetrahedra);
    for (const auto& v: vertices){
        cout<<v.x<<" "<<v.y<<" "<<v.z<<endl;
    }
    for (const auto& t: tetrahedra){
        cout<<t.v1<<" "<<t.v2<<" "<<t.v3<<" "<<t.v4<<endl;
    }
}