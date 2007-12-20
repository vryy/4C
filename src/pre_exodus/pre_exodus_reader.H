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

extern "C" /* stuff which is c and is accessed from c++ */
{
#include <exodusII.h>
}

using namespace std;

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

  //! Print mesh info
  void Print(ostream& os) const;
  
private:
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

};



#endif /*PRE_EXODUS_READER_H_*/
#endif
