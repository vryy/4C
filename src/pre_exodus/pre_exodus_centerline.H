/*----------------------------------------------------------------------*/
/*!
\file pre_exodus_centerline.H

\brief pre_exodus centerline definition

\level 1

\maintainer Martin Kronbichler

 */
/*----------------------------------------------------------------------*/

#ifndef PRE_EXODUS_CENTERLINE_H_
#define PRE_EXODUS_CENTERLINE_H_


#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "pre_exodus_reader.H" //contains class Mesh (necessary for method - element_cosys)

// Helper function to sort list of pairs according to second entry of the pair
bool MyDataSortPredicate(std::pair<int, double> lhs,std::pair<int, double> rhs);

namespace EXODUS
{

std::map<int,std::map<int,std::vector<std::vector<double> > > > EleCenterlineInfo(std::string& cline,EXODUS::Mesh& mymesh, const std::vector<double> coordcorr);

std::map<int,double> NdCenterlineThickness(std::string cline,const std::set<int>& nodes, const std::map<int,std::vector<int> >& conn, const EXODUS::Mesh& mesh, const double ratio, const std::vector<double> coordcorr);

class Centerline
{
private:
  Teuchos::RCP<std::map<int,std::vector<double> > > points_; /// centerline points
  Teuchos::RCP<std::map<int,double> > diam_; /// diameter of the best fit circle in centerline point

public:
  //ctor
  Centerline(std::string,std::vector<double> coordcorr);
  Centerline(const EXODUS::NodeSet& ns, const Teuchos::RCP<std::map<int,std::vector<double> > > nodes);
  //dtor
  ~Centerline();
  //returns points_
  Teuchos::RCP<std::map<int,std::vector<double> > > GetPoints() const {return points_;};
  //returns diam_
  Teuchos::RCP<std::map<int,double> > GetDiams() const {return diam_;};
  //displays points_ on console
  void PrintPoints();
  //creates gmsh-file to visualize points
  void PlotCL_Gmsh();

};

//! calculates distance between two 3-dim vectors
double distance3d(std::vector<double>,std::vector<double>);
//! calculates difference of two 3-dim vectors
std::vector<double> substract3d(std::vector<double>,std::vector<double>);
//! calculates sum of two 3-dim vectors
std::vector<double> add3d(std::vector<double>,std::vector<double>);
//! calculates cross product of two 3-dim vectors
std::vector<double> cross_product3d(std::vector<double>,std::vector<double>);
//! calculates scalar product of two 3-dim vectors
double scalar_product3d(std::vector<double> v1,std::vector<double> v2);
//! normalizes a vector
void normalize3d(std::vector<double>&);
//! creates local coordinate systems for each element and returns resulting map
std::map<int,std::map<int,std::vector<std::vector<double> > > > element_cosys(EXODUS::Centerline&,const EXODUS::Mesh&, const std::vector<int>&);
/// creates local coordinate systems for each element and returns resulting map
///different calculation method for element normals applied in this function
std::map<int,std::map<int,std::vector<std::vector<double> > > > element_cosys(EXODUS::Centerline&,const EXODUS::Mesh&, const std::vector<int>&, std::set<int>& );

/*! creates degenerated local coordinate system in case of given centerpoint for each element and returns resulting map.
 * CoSys contains only three copies of the same vector!  */
std::map<int,std::map<int,std::vector<std::vector<double> > > > element_degcosys
(
  EXODUS::Centerline&,  ///< centerpoints
  const EXODUS::Mesh&,  ///< mech
  const std::vector<int>&    ///< element ids
);

//! creates local coordinate systems for each element and creates gmsh-file for visualizing
void PlotCosys(EXODUS::Centerline&,const EXODUS::Mesh&,const std::vector<int>&);

// Helper function to sort list of pairs according to second entry of the pair
//bool MyDataSortPredicate(const pair<int, double>& lhs, const pair<int, double>& rhs);

}// of namespace EXODUS

#endif /*PRE_EXODUS_CENTERLINE_H_*/
