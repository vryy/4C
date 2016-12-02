/*----------------------------------------------------------------------*/
/*!
\file pre_exodus_reader.H

\brief preprocessor reader for exodusII format

\maintainer Martin Kronbichler

\level 1

Here everything related with the exodus format and the accessible data
is handed to a c++ object mesh.
*/
/*----------------------------------------------------------------------*/
#ifndef PRE_EXODUS_READER_H_
#define PRE_EXODUS_READER_H_

#include <string>
#include <iostream>
#include <vector>
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_element.H"

namespace EXODUS
{

//forward declaration
class ElementBlock;
class NodeSet;
class SideSet;

/*!
\brief Mesh will in future store all information necessary to build a mesh

\author maf (frenzel@lnm.mw.tum.de)
*/
class Mesh
{
public:
  //! constructor
  Mesh(std::string exofilename);

  //! extension constructor, adds Elementblock and nodes to a basemesh
  Mesh(const Mesh& basemesh,
       const Teuchos::RCP<std::map<int,std::vector<double> > > extNodes,
       const std::map<int,Teuchos::RCP<ElementBlock> >& extBlocks,
       const std::map<int,NodeSet>& extNodesets,
       const std::map<int,SideSet>& extSidesets,
       const std::string newtitle);

  //! empty constructor
  Mesh();

  //! destructor
  virtual ~Mesh();

  //! Print mesh info
  void Print(std::ostream& os, bool verbose = false) const;

  //! Print Nodes and Coords
  void PrintNodes(std::ostream& os, bool storeid=false) const;

  //! Get numer of nodes in mesh
  int GetNumNodes() const {return nodes_->size();}

  //! Get number of elements in mesh
  int GetNumEle() const {return num_elem_;}

  //! Get number of dimensions
  int GetNumDim() const {return num_dim_;}

  //! Get number of dimensions
  int GetBACIDim() const {return baci_dim_;}

  //! Get exodus file id
  int GetExoId() const {return exoid_;}

  //! Get mesh title
  std::string GetTitle() const;

  //! Get ElementBlock map
  std::map<int,Teuchos::RCP<ElementBlock> > GetElementBlocks() const {return elementBlocks_;}

  //! Get Number of ElementBlocks
  int GetNumElementBlocks() const {return elementBlocks_.size();}

  //! Get one ElementBlock
  Teuchos::RCP<ElementBlock> GetElementBlock(const int id) const;

  //! Get NodeSet map
  std::map<int,NodeSet> GetNodeSets() const {return nodeSets_;}

  //! Get Number of NodeSets
  int GetNumNodeSets() const {return nodeSets_.size();}

  //! Get one NodeSet
  NodeSet GetNodeSet(const int id) const;

  //! Get SideSet map
  std::map<int,SideSet> GetSideSets() const {return sideSets_;}

  //! Get Number of SideSets
  int GetNumSideSets() const {return sideSets_.size();}

  //! Get one SideSet
  SideSet GetSideSet(const int id) const;

  //! Get Side Set Connectivity with Global Nodes
  std::map<int,std::vector<int> > GetSideSetConn(const SideSet sideset) const;

  //! Get Side Set Connectivity with Global Nodes
  std::map<int,std::vector<int> > GetSideSetConn(const SideSet sideset, bool checkoutside) const;

  //! Make sure child ele (SideSet) is outward oriented w.r.t. parent ele
  std::vector<int> OutsideOrientedSide(const std::vector<int> parentele, const std::vector<int> sidemap) const;

  //! Get egde Normal at node
  std::vector<double> Normal(const int head1, const int origin, const int head2) const;

  //! Get normalized Vector between 2 nodes
  std::vector<double> NodeVec(const int tail, const int head) const;

  //! Transform SideSet into ElementBlock
  std::vector<ElementBlock> SideSetToEBlocks(const SideSet& sideset, const std::map<int,std::vector<int> >& sidesetconn) const;

  //! Transform SideSet into NodeSet
  NodeSet SideSetToNodeSet(const SideSet& sideset, const std::map<int,std::vector<int> >& sideconn) const;

  //! Get Set of Nodes in SideSet
  std::set<int> GetSideSetNodes(const EXODUS::SideSet& sideset, const std::map<int,std::vector<int> >& sideconn) const;

  //! Get Node map
  Teuchos::RCP<std::map<int,std::vector<double> > > GetNodes() const {return nodes_;}

  //! Get one nodal coords
  std::vector<double> GetNode(const int NodeID) const;

  //! Set one nodal coords
  void SetNode(const int NodeID, const std::vector<double> coord);

  //! Set number of space dimensions
  void SetNsd (const int nsd);

  //! Close Exodus File
  void CloseExo() const;

  //! Write Mesh into exodus file
  void WriteMesh(const std::string newexofilename) const;

  //! Add Element Block to mesh
  void AddElementBlock(const Teuchos::RCP<EXODUS::ElementBlock> eblock) const;

  //! Erase Element Block from mesh
  void EraseElementBlock(const int id);

  //! Erase SideSet from mesh
  void EraseSideSet(const int id);

  //! Calculate the midpoint of all elements and return map<midpoint,map<eb,ele> >
  std::map<int,std::pair<int,int> > createMidpoints(std::map<int,std::vector<double> >& midpoints, const std::vector<int>& eb_ids) const;

  //! Adjust local element ids referenced in SideSet to global ids
  std::map<int,std::vector<int> > GlobalifySSeleids(const int ssid) const;

  //! Plot Nodes in Gmsh-file
  void PlotNodesGmsh() const;

  //! Plot all ElementBlocks into Gmsh-file
  void PlotElementBlocksGmsh(const std::string fname,const EXODUS::Mesh& mymesh) const;
  void PlotElementBlocksGmsh(const std::string fname,const EXODUS::Mesh& mymesh, const std::vector<int>& ebids) const;

  //! Plot Connectivity into Gmsh-file
  void PlotConnGmsh(const std::string fname,const EXODUS::Mesh& mymesh, const std::map<int,std::vector<int> >& conn) const;

private:
  Teuchos::RCP<std::map<int,std::vector<double> > > nodes_;

  std::map<int,Teuchos::RCP<ElementBlock> > elementBlocks_;

  std::map<int,NodeSet> nodeSets_;

  std::map<int,SideSet> sideSets_;

  //! number of dimensions
  int num_dim_;
  //! number of dimensions for BACI problem (wall and fluid2 elements require 2d, although we have spatial dimensions)
  int baci_dim_;
  //! number of elements
  int num_elem_;
  //! exoid
  int exoid_;
  //! title
  std::string title_;
};


/*!
\brief ElementBlock is a set of Elements of same discretization Type

A Element Block is a tiny class storing element-type, name, etc. of a ElementBlock
It implements its printout.

\author maf (frenzel@lnm.mw.tum.de)
*/
class ElementBlock
{
public:
  enum Shape{
    dis_none,                   ///< unknown dis type
    quad4,                      ///< 4 noded quadrilateral
    quad8,                      ///< 8 noded quadrilateral
    quad9,                      ///< 9 noded quadrilateral
    shell4,
    shell8,
    shell9,
    tri3,                       ///< 3 noded triangle
    tri6,                       ///< 6 noded triangle
    hex8,                       ///< 8 noded hexahedra
    hex20,                      ///< 20 noded hexahedra
    hex27,                      ///< 27 noded hexahedra
    tet4,                       ///< 4 noded tetrahedra
    tet10,                      ///< 10 noded tetrahedra
    wedge6,                     ///< 6 noded wedge
    wedge15,                    ///< 15 noded wedge
    pyramid5,                   ///< 5 noded pyramid
    bar2,                      ///< 2 noded line
    bar3,                      ///< 3 noded line
    point1,                     ///< 1 noded point
    max_distype                 ///<  end marker. must be the last entry
  };

  ElementBlock(ElementBlock::Shape DisType,
               Teuchos::RCP<std::map<int,std::vector<int> > >& eleconn,    // Element connectivity
               std::string name);

  virtual ~ElementBlock();

  ElementBlock::Shape GetShape() const {return distype_;}

  int GetNumEle() const {return eleconn_->size();}

  Teuchos::RCP<std::map<int,std::vector<int> > > GetEleConn() const {return eleconn_;}

  std::vector<int> GetEleNodes(int i) const;

  std::string GetName() const {return name_;}

  int GetEleNode(int ele, int node) const;

  void FillEconnArray(int *connarray) const;

  void Print(std::ostream& os, bool verbose=false) const;

private:
  Shape distype_;

  //! Element Connectivity
  Teuchos::RCP<std::map<int,std::vector<int> > > eleconn_;

  std::string name_;
};

class NodeSet
{
public:
  NodeSet(const std::set<int>& nodeids, const std::string& name, const std::string& propname);

  virtual ~NodeSet();

  std::set<int> GetNodeSet() const {return nodeids_;};

  std::string GetName() const {return name_;};

  std::string GetPropName() const {return propname_;};

  void FillNodelistArray(int *nodelist) const;

  inline int GetNumNodes() const {return nodeids_.size();}

  void Print(std::ostream& os, bool verbose=false) const;

private:
  std::set<int> nodeids_;  // nodids in NodeSet
  std::string name_;       // NodeSet name
  std::string propname_;   // Icem assignes part names as property names
};

class SideSet
{
public:
  SideSet(const std::map<int,std::vector<int> >& sides, const std::string& name);

  virtual ~SideSet();

  inline int GetNumSides() const {return sides_.size();}

  std::string GetName() const {return name_;}

  std::map<int,std::vector<int> > GetSideSet() const {return sides_;}

  void ReplaceSides(std::map<int,std::vector<int> >newsides) {sides_ = newsides; return;};

  std::vector<int> GetFirstSideSet() const {return sides_.begin()->second;}

  void FillSideLists(int* elemlist, int* sidelist) const;
  void FillSideLists(int* elemlist, int* sidelist, const std::map<int,std::vector<int> >& sides) const;

  void Print(std::ostream& os, bool verbose=false) const;

private:
  std::map<int,std::vector<int> > sides_;
  std::string name_;
};
// *********** end of classes

Mesh QuadtoTri(EXODUS::Mesh& basemesh);

inline ElementBlock::Shape StringToShape(const std::string shape)
{
  if      (shape.compare("SPHERE")== 0)   return ElementBlock::point1;
  else if (shape.compare("QUAD4")== 0)    return ElementBlock::quad4;
  else if (shape.compare("QUAD8")== 0)    return ElementBlock::quad8;
  else if (shape.compare("QUAD9")== 0)    return ElementBlock::quad9;
  else if (shape.compare("SHELL4")== 0)   return ElementBlock::shell4;
  else if (shape.compare("SHELL8")== 0)   return ElementBlock::shell8;
  else if (shape.compare("SHELL9")== 0)   return ElementBlock::shell9;
  else if (shape.compare("TRI3")== 0)     return ElementBlock::tri3;
  else if (shape.compare("TRI6")== 0)     return ElementBlock::tri6;
  else if (shape.compare("HEX8")== 0)     return ElementBlock::hex8;
  else if (shape.compare("HEX20")== 0)    return ElementBlock::hex20;
  else if (shape.compare("HEX27")== 0)    return ElementBlock::hex27;
  else if (shape.compare("HEX")== 0)      return ElementBlock::hex8; // really needed????? a.g. 08/08
  else if (shape.compare("TET4")== 0)     return ElementBlock::tet4; // TODO:: gibts das eigentlich?
  else if (shape.compare("TETRA4")== 0)   return ElementBlock::tet4;
  else if (shape.compare("TETRA10")== 0)  return ElementBlock::tet10;
  else if (shape.compare("TETRA")== 0)    return ElementBlock::tet4;  // really needed????? a.g. 08/08
  else if (shape.compare("WEDGE6")== 0)   return ElementBlock::wedge6;
  else if (shape.compare("WEDGE15")== 0)  return ElementBlock::wedge15;
  else if (shape.compare("WEDGE")== 0)    return ElementBlock::wedge6;   // really needed????? a.g. 08/08
  else if (shape.compare("PYRAMID5")== 0) return ElementBlock::pyramid5;
  else if (shape.compare("PYRAMID")== 0)  return ElementBlock::pyramid5; // really needed????? a.g. 08/08
  else if (shape.compare("BAR2")== 0)     return ElementBlock::bar2;
  else if (shape.compare("BAR3")== 0)     return ElementBlock::bar3;
  else{
    std::cout << "Unknown Exodus Element Shape Name: " << shape;
    dserror("Unknown Exodus Element Shape Name!"); return ElementBlock::dis_none;
  }
}

inline std::string ShapeToString(const ElementBlock::Shape shape)
{
  switch (shape){
  case ElementBlock::point1:   return "SPHERE"; break;
  case ElementBlock::quad4:    return "QUAD4";  break;
  case ElementBlock::quad8:    return "QUAD8";  break;
  case ElementBlock::quad9:    return "QUAD9";  break;
  case ElementBlock::shell4:   return "SHELL4";  break;
  case ElementBlock::shell8:   return "SHELL8"; break;
  case ElementBlock::shell9:   return "SHELL9"; break;
  case ElementBlock::tri3:     return "TRI3";   break;
  case ElementBlock::tri6:     return "TRI6";   break;
  case ElementBlock::hex8:     return "HEX8";   break;
  case ElementBlock::hex20:    return "HEX20";  break;
  case ElementBlock::hex27:    return "HEX27";  break;
  case ElementBlock::tet4:     return "TET4";   break;
  case ElementBlock::tet10:    return "TET10";  break;
  case ElementBlock::wedge6:   return "WEDGE6"; break;
  case ElementBlock::wedge15:  return "WEDGE15";break;
  case ElementBlock::pyramid5: return "PYRAMID5";break;
  case ElementBlock::bar2:     return "BAR2";   break;
  case ElementBlock::bar3:     return "BAR3";   break;
  default: dserror("Unknown ElementBlock::Shape");
  }
  return "xxx";
}

inline DRT::Element::DiscretizationType PreShapeToDrt(const ElementBlock::Shape shape)
{
  switch (shape){
  case ElementBlock::point1:   return DRT::Element::point1; break;
  case ElementBlock::quad4:    return DRT::Element::quad4;  break;
  case ElementBlock::quad8:    return DRT::Element::quad8;  break;
  case ElementBlock::quad9:    return DRT::Element::quad9;  break;
  case ElementBlock::shell4:   return DRT::Element::quad4;  break;
  case ElementBlock::shell8:   return DRT::Element::quad8;  break;
  case ElementBlock::shell9:   return DRT::Element::quad9;  break;
  case ElementBlock::tri3:     return DRT::Element::tri3;   break;
  case ElementBlock::tri6:     return DRT::Element::tri6;   break;
  case ElementBlock::hex8:     return DRT::Element::hex8;   break;
  case ElementBlock::hex20:    return DRT::Element::hex20;  break;
  case ElementBlock::hex27:    return DRT::Element::hex27;  break;
  case ElementBlock::tet4:     return DRT::Element::tet4;   break;
  case ElementBlock::tet10:    return DRT::Element::tet10;  break;
  case ElementBlock::wedge6:   return DRT::Element::wedge6; break;
  case ElementBlock::wedge15:  return DRT::Element::wedge15;break;
  case ElementBlock::pyramid5: return DRT::Element::pyramid5;break;
  case ElementBlock::bar2:     return DRT::Element::line2;   break;
  case ElementBlock::bar3:     return DRT::Element::line3;   break;
  default: dserror("Unknown ElementBlock::Shape");
  }
  return DRT::Element::max_distype;

}

void PrintMap(std::ostream& os,const std::map<int,std::vector<int> > mymap);
void PrintMap(std::ostream& os, const std::map<int,std::vector<double> >mymap);
void PrintMap(std::ostream& os,const std::map<int,std::set<int> > mymap);
void PrintMap(std::ostream& os,const std::map<int,std::map<int,int> > mymap);
void PrintMap(std::ostream& os,const std::map<int,std::pair<int,int> > mymap);
void PrintMap(std::ostream& os,const std::map<int,int> mymap);
void PrintMap(std::ostream& os,const std::map<int,double> mymap);
void PrintMap(std::ostream& os,const std::map<double,int> mymap);
void PrintVec(std::ostream& os, const std::vector<int> actvec);
void PrintVec(std::ostream& os, const std::vector<double> actvec);
void PrintSet(std::ostream& os, const std::set<int> actset);


int HexSideNumberExoToBaci(const int exoface);
int PyrSideNumberExoToBaci(const int exoface);



} // of namespace EXODUS

#endif /*PRE_EXODUS_READER_H_*/
