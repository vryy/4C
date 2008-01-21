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
#ifdef EXODUS
#ifndef PRE_EXODUS_READER_H_
#define PRE_EXODUS_READER_H_

#include <string>
#include <iostream>
#include <vector>
#include <Teuchos_RefCountPtr.hpp>
#include "../drt_lib/drt_dserror.H"

extern "C" /* stuff which is c and is accessed from c++ */
{
#include <exodusII.h>
}

using namespace std;
using namespace Teuchos;

//forward declaration
class Entity;
class PreNode;

/*!
\class Mesh

\brief Mesh will in future store all information necessary to build a mesh


\author maf (frenzel@lnm.mw.tum.de)
*/
class Mesh
{
public:
  enum EntityType
  {
    elem_blk,
    node_set,
    side_set,
  };
  //! constructor
  Mesh(string exofilename);
  
  //! destructor
  virtual ~Mesh();
  
  //! get number of entities
  int GetNumberEntities(){ return num_entities_;};
  
  //! get one entitiy
  RCP<Entity> GetEntity(int EntityID){ return myEntities_[EntityID];};
  
  //! Print mesh info
  void Print(ostream& os) const;
  
  //! Get numer of nodes in mesh
  int GetNumNodes() const {return num_nodes_;}
  
  //! Get number of elements in mesh
  int GetNumEle() const {return num_elem_;}
  
  //! Get exodus file id
  int GetExoId() const {return exoid_;}
  
  //! Get number of Entities in Mesh 
  int GetNumEntities() const {return num_entities_;}
  
//  //! Get nodal coordinates
//  int GetNodalCorrd_x() const {return x_;}
//
//  //! Get nodal coordinates
//  int GetNodalCorrd_y() const {return y_;}
//
//  //! Get nodal coordinates
//  int GetNodalCorrd_z() const {return z_;}
  
  RCP<PreNode> GetNode(int NodeID){ return myNodes_[NodeID];};

  void CloseExo();

private:
  
  vector<RCP<Entity> > myEntities_;
  
  vector<RCP<PreNode> > myNodes_;
  
  //! number of dimensions
  int num_dim_;
  //! number of nodes
  int num_nodes_;
  //! number of elements
  int num_elem_;
  //! number of element blocks
  int num_elem_blk_;
  //! number of node sets
  int num_node_sets_;
  //! number of side sets
  int num_side_sets_;
  //! exoid
  int exoid_;
  //! title
  char title_[MAX_LINE_LENGTH+1];
  
  //! number of my entities
  int num_entities_;
  
  //! nodal coordinates
  float x_, y_, z_;
  
  
};

/*!
\class Entity

\brief Entity of mesh such as ElementBlock, NodeSet, SideSet

A mesh entity stores as many properties as possible,
such as names, number of nodes, number of elements, element types, ...
It implements its printout.

\author maf (frenzel@lnm.mw.tum.de)
*/
class Entity
{
public:
  //! Type of entity (ElementBlock, SideSet, or NodeSet)
  enum EntityType
  {
    elem_blk,
    node_set,
    side_set,
  };
    
  //! constructor
  Entity(int exoid, int entityID, int typeID, EntityType entitytype);
  
  //! destructor
  virtual ~Entity();
  
  //! Set Property name, designated by ICEM to available parts
  void SetPropertyName(string prop_name);
  
  //! Get element type
  string GetElementType() const {return elem_type_;};
  
  //! Get num_nodes_
  int GetNumNodes() const {return num_nodes_;};
  
  //! GetNumNodPerElem
  int GetNumNodpElem() const {return num_nod_per_elem_;};
  
  //! GetNumEle
  int GetNumEle() const {return num_el_in_blk_;};
  
  //! GetEntityId()
  int GetEntityId() const {return entityID_;};

  //! Print entity
  void Print(ostream& os) const;
  
private:
  //! entity ID
  int entityID_;
  //! element type
  string elem_type_;
  //! number of elements in block
  int num_el_in_blk_;
  //! number of nodes per element
  int num_nod_per_elem_;
  //! number of attributes
  int num_attr_;
  //! number of nodes in entity
  int num_nodes_;
  //! entity name
  string entity_name_;
  //! entity type name
  string entity_type_;
  //! property name (so far provided from ICEM for NodeSets
  string entity_prop_name_;
  //! element connectivity (only for element blocks)
  vector<int> elem_conn_;

};


class PreElement
{
public:
  enum ShapeType{
    quad4,
    tri3,
  };
  PreElement(int id, ShapeType shape);
  
  PreElement(const PreElement& old);
  
  virtual ~PreElement();
  
  inline int Id() const {return id_; }
  
  inline ShapeType Shape() const {return shape_;}
  
  inline int NumNode() const {return nodeid_.size(); }
  
  inline const int* NodeIds() const {if (nodeid_.size()) return &nodeid_[0];
                                     else                return NULL;     }
  
  PreNode** Nodes() { if (node_.size()) return &node_[0]; else return NULL; }
  
  void Print(ostream& os) const;
  
  void SetNodeIds(const int nnode, const int* nodes);
  
private:
  int id_;
  
  ShapeType shape_;
  
  vector<int> nodeid_;
  
  vector<PreNode*> node_;
  
};

class PreNode
{
public:
  PreNode(int id, const double* coords);
  
  PreNode(const PreNode& old);
  
  virtual ~PreNode();
  
  inline int Id() const {return id_;}
  
  inline const double* X() const {return x_; }
  
  inline int NumElement() const {return element_.size(); }
  
  PreElement** Elements()   { if (NumElement()) return &element_[0]; else return NULL; }
  
  void Print(ostream& os) const;
  
  void AddElementPtr(PreElement* eleptr)
  { const int size = element_.size(); element_.resize(size+1);
    element_[size] = eleptr; return;}

private:
  int id_;
  
  double x_[3];
  
  vector<PreElement*> element_;
  
};



#endif /*PRE_EXODUS_READER_H_*/
#endif
