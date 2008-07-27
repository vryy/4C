//"pre_exodus_centerline.H"

#ifndef PRE_EXODUS_CENTERLINE_H_
#define PRE_EXODUS_CENTERLINE_H_


#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "pre_exodus_reader.H" //contains class Mesh (necessary for method - element_cosys)
using namespace std;


namespace EXODUS
{

map<int,map<int,vector<vector<double> > > > EleCenterlineInfo(string& cline,EXODUS::Mesh& mymesh);

class Centerline
{
private:
	RCP<map<int,vector<double> > > points_;
public:
	//ctor
	Centerline(string);
	Centerline(const EXODUS::NodeSet& ns, const RCP<map<int,vector<double> > > nodes);
	//dtor
	~Centerline();
	//returns points_
	RCP<map<int,vector<double> > > GetPoints() const {return points_;};
	//displays points_ on console
	void PrintPoints();
	//creates gmsh-file to visualize points
	void PlotCL_Gmsh();
	
};

//calculates distance between two 3-dim vectors
double distance3d(vector<double>,vector<double>);
//calculates difference of two 3-dim vectors
vector<double> substract3d(vector<double>,vector<double>);
//calculates sum of two 3-dim vectors
vector<double> add3d(vector<double>,vector<double>);
//calculates cross product of two 3-dim vectors
vector<double> cross_product3d(vector<double>,vector<double>);
//normalizes a vector
void normalize3d(vector<double>&);
//creates local coordinate systems for each element and returns resulting map
map<int,map<int,vector<vector<double> > > > element_cosys(EXODUS::Centerline&,const EXODUS::Mesh&);
//creates local coordinate systems for each element and creates gmsh-file for visualizing
void PlotCosys(EXODUS::Centerline&,const EXODUS::Mesh&);

}// of namespace EXODUS

#endif /*PRE_EXODUS_CENTERLINE_H_*/
