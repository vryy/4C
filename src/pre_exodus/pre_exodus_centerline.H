//"pre_exodus_centerline.H"

#ifndef PRE_EXODUS_CENTERLINE_H_
#define PRE_EXODUS_CENTERLINE_H_


#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "pre_exodus_reader.H" //contains class Mesh (necessary for method - element_cosys)

// Helper function to sort list of pairs according to second entry of the pair
bool MyDataSortPredicate(pair<int, double> lhs,pair<int, double> rhs);

namespace EXODUS
{

map<int,map<int,vector<vector<double> > > > EleCenterlineInfo(string& cline,EXODUS::Mesh& mymesh, const vector<double>coordcorr);

map<int,double> NdCenterlineThickness(string cline,const std::set<int>& nodes, const map<int,vector<int> >& conn, const EXODUS::Mesh& mesh, const double ratio, const vector<double> coordcorr);

class Centerline
{
private:
  RCP<std::map<int,vector<double> > > points_; /// centerline points
  RCP<std::map<int,double> > diam_; /// diameter of the best fit circle in centerline point

public:
  //ctor
  Centerline(string,vector<double> coordcorr);
  Centerline(const EXODUS::NodeSet& ns, const RCP<map<int,vector<double> > > nodes);
  //dtor
  ~Centerline();
  //returns points_
  RCP<map<int,vector<double> > > GetPoints() const {return points_;};
  //returns diam_
  RCP<map<int,double> > GetDiams() const {return diam_;};
  //displays points_ on console
  void PrintPoints();
  //creates gmsh-file to visualize points
  void PlotCL_Gmsh();

};

//! calculates distance between two 3-dim vectors
double distance3d(vector<double>,vector<double>);
//! calculates difference of two 3-dim vectors
vector<double> substract3d(vector<double>,vector<double>);
//! calculates sum of two 3-dim vectors
vector<double> add3d(vector<double>,vector<double>);
//! calculates cross product of two 3-dim vectors
vector<double> cross_product3d(vector<double>,vector<double>);
//! calculates scalar product of two 3-dim vectors
double scalar_product3d(vector<double> v1,vector<double> v2);
//! normalizes a vector
void normalize3d(vector<double>&);
//! creates local coordinate systems for each element and returns resulting map
map<int,map<int,vector<vector<double> > > > element_cosys(EXODUS::Centerline&,const EXODUS::Mesh&, const vector<int>&);
/// creates local coordinate systems for each element and returns resulting map 
///different calculation method for element normals applied in this function
map<int,map<int,vector<vector<double> > > > element_cosys(EXODUS::Centerline&,const EXODUS::Mesh&, const vector<int>&, std::set<int>& );

/*! creates degenerated local coordinate system in case of given centerpoint for each element and returns resulting map.
 * CoSys contains only three copies of the same vector!  */
map<int,map<int,vector<vector<double> > > > element_degcosys
(
  EXODUS::Centerline&,  ///< centerpoints
  const EXODUS::Mesh&,  ///< mech
  const vector<int>&    ///< element ids
);

//! creates local coordinate systems for each element and creates gmsh-file for visualizing
void PlotCosys(EXODUS::Centerline&,const EXODUS::Mesh&,const vector<int>&);

// Helper function to sort list of pairs according to second entry of the pair
//bool MyDataSortPredicate(const pair<int, double>& lhs, const pair<int, double>& rhs);

}// of namespace EXODUS

#endif /*PRE_EXODUS_CENTERLINE_H_*/

