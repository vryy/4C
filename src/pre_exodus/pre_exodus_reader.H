/*----------------------------------------------------------------------*/
/*!
\file pre_exodus_reader.H

\brief preprocessor reader for exodusII format 

<pre>
Maintainer: Moritz & Georg
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/frenzel
            089 - 289-15240
</pre>

Here everything related with the exodus format and the accessible data
is handed to a c++ object mesh.
*/
/*----------------------------------------------------------------------*/
#ifdef D_EXODUS
#ifndef PRE_EXODUS_READER_H_
#define PRE_EXODUS_READER_H_

#include <string>
#include <iostream>
#include <vector>
#include <Teuchos_RefCountPtr.hpp>
#include "../drt_lib/drt_dserror.H"
//include "pre_node.H"

extern "C" /* stuff which is c and is accessed from c++ */
{
#include <exodusII.h>
}

using namespace std;
using namespace Teuchos;

namespace EXODUS
{

// this minifunction maps exodus ids starting from 1 to store ids starting from 0
inline const int eid(const int storeid) {return storeid-1;}


//forward declaration
class ElementBlock;
class NodeSet;
class SideSet;
class PreNode;

/*!
\class Mesh

\brief Mesh will in future store all information necessary to build a mesh


\author maf (frenzel@lnm.mw.tum.de)
*/
class Mesh
{
public:
  //! constructor
  Mesh(string exofilename);
  
  //! destructor
  virtual ~Mesh();
  
  //! Print mesh info
  void Print(ostream& os, bool verbose = false) const;
  
  //! Get numer of nodes in mesh
  int GetNumNodes() const {return num_nodes_;}
  
  //! Get number of elements in mesh
  int GetNumEle() const {return num_elem_;}
  
  //! Get exodus file id
  int GetExoId() const {return exoid_;}
  
  //! Get ElementBlock map
  map<int,ElementBlock> GetElementBlocks() const {return elementBlocks_;}
  
  //! Get Number of ElementBlocks
  int GetNumElementBlocks() const {return elementBlocks_.size();}
  
  //! Get NodeSet map
  map<int,NodeSet> GetNodeSets() const {return nodeSets_;}
  
  //! Get Number of ElementBlocks
  int GetNumNodeSets() const {return nodeSets_.size();}
  
  //! Get SideSet map
  map<int,SideSet> GetSideSets() const {return sideSets_;}
  
  //! Get Number of ElementBlocks
  int GetNumSideSets() const {return sideSets_.size();}
  
  RCP<PreNode> GetNode (const int NodeID) const { const int id = eid(NodeID); return myNodes_[id];};

  //! Close Exodus File
  void CloseExo();
  
  //! Write Mesh into exodus file
  void WriteMesh(string newexofilename);
  
  //! Add Element Block to mesh
  void AddElementBlock(const EXODUS::ElementBlock eblock) const;

private:
  vector<RCP<PreNode> > myNodes_;
  
  map<int,ElementBlock> elementBlocks_;
  
  map<int,NodeSet> nodeSets_;
  
  map<int,SideSet> sideSets_;
  
  //! number of dimensions
  int num_dim_;
  //! number of nodes
  int num_nodes_;
  //! number of elements
  int num_elem_;
  //! exoid
  int exoid_;
  //! title
  char title_[MAX_LINE_LENGTH+1];
};


/*!
\class ElementBlock

\brief ElementBlock is a set of Elements of same discretization Type

A Element Block is a tiny class storing element-type, name, etc. of a ElementBlock
It implements its printout.

\author maf (frenzel@lnm.mw.tum.de)
*/
class ElementBlock
{
public:
  enum Shape{
    dis_none,                   // unknown dis type
    quad4,                      // 4 noded quadrilateral
    quad8,                      // 8 noded quadrilateral
    quad9,                      // 9 noded quadrilateral
    shell4,
    shell8,
    shell9,
    tri3,                       // 3 noded triangle
    tri6,                       // 6 noded triangle
    hex8,                       // 8 noded hexahedra
    hex20,                      // 20 noded hexahedra
    hex27,                      // 27 noded hexahedra
    tet4,                       // 4 noded tetrahedra
    tet10,                      // 10 noded tetrahedra
    wedge6,                     // 6 noded wedge
    wedge15,                    // 15 noded wedge
    pyramid5,                   // 5 noded pyramid
    bar2,                      // 2 noded line
    bar3,                      // 3 noded line
    point1,                     // 1 noded point
    max_distype                 //  end marker. must be the last entry
  };

  ElementBlock(ElementBlock::Shape DisType,
               map<int,vector<int> > &eleconn,    // Element connectivity
               string name);
               
  virtual ~ElementBlock();
  
  ElementBlock::Shape GetShape() const {return distype_;}
  
  int GetNumEle() const {return eleconn_.size();}
  
  map<int,vector<int> > GetEleConn() const {return eleconn_;}
  
  vector<int> GetEleNodes(int i) const;
  
  string GetName() const {return name_;}
  
  int GetEleNode(int ele, int node) const;
  
  void FillEconnArray(int *connarray) const;
  
  void Print(ostream& os, bool verbose=false) const;
  
private:
  Shape distype_;
  
  // Element Connectivity
  map<int,vector<int> > eleconn_;
  
  string name_;
};

class NodeSet
{
public:
  NodeSet(set<int> nodeids, string name, string propname);
  
  virtual ~NodeSet();
  
  set<int> GetNodeSet() const {return nodeids_;};
  
  string GetName() const {return name_;};
  
  string GetPropName() const {return propname_;};
  
  void FillNodelistArray(int *nodelist) const;
  
  inline int GetNumNodes() const {return nodeids_.size();}
  
  void Print(ostream& os, bool verbose=false) const;

private:
  set<int> nodeids_;  // nodids in NodeSet
  string name_;       // NodeSet name
  string propname_;   // Icem assignes part names as property names
};

class SideSet
{
public:
  SideSet(string name);
  virtual ~SideSet();
  void Print(ostream& os) const {os << "SideSet, named: " << name_.c_str() << endl;};
private:
  string name_;
};


inline ElementBlock::Shape StringToShape(const string shape)
{
  if      (shape.compare(0,6,"SPHERE")== 0)   return ElementBlock::point1;
  else if (shape.compare(0,5,"QUAD4")== 0)    return ElementBlock::quad4;
  else if (shape.compare(0,5,"QUAD8")== 0)    return ElementBlock::quad8;
  else if (shape.compare(0,5,"QUAD9")== 0)    return ElementBlock::quad9;
  else if (shape.compare(0,6,"SHELL4")== 0)   return ElementBlock::shell4;
  else if (shape.compare(0,6,"SHELL8")== 0)   return ElementBlock::shell8;
  else if (shape.compare(0,6,"SHELL9")== 0)   return ElementBlock::shell9;
  else if (shape.compare(0,4,"TRI3")== 0)     return ElementBlock::tri3;
  else if (shape.compare(0,4,"TRI6")== 0)     return ElementBlock::tri6;
  else if (shape.compare(0,4,"HEX8")== 0)     return ElementBlock::hex8;
  else if (shape.compare(0,5,"HEX20")== 0)    return ElementBlock::hex20;
  else if (shape.compare(0,5,"HEX27")== 0)    return ElementBlock::hex27;
  else if (shape.compare(0,4,"TET4")== 0)     return ElementBlock::tet4;
  else if (shape.compare(0,5,"TET10")== 0)    return ElementBlock::tet10;
  else if (shape.compare(0,6,"WEDGE6")== 0)   return ElementBlock::wedge6;
  else if (shape.compare(0,7,"WEDGE15")== 0)  return ElementBlock::wedge15;
  else if (shape.compare(0,8,"PYRAMID5")== 0) return ElementBlock::pyramid5;
  else if (shape.compare(0,4,"BAR2")== 0)     return ElementBlock::bar2;
  else if (shape.compare(0,4,"BAR3")== 0)     return ElementBlock::bar3;
  else{
    dserror("Unknown Exodus Element Shape Name"); return ElementBlock::dis_none;
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
  case ElementBlock::tet4:     return "HEX27";  break;
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

} // of namespace EXODUS

#endif /*PRE_EXODUS_READER_H_*/
#endif
