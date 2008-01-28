/*----------------------------------------------------------------------*/
/*!
\file pre_node.H

\brief preprocessor reader for exodusII format 

<pre>
Maintainer: Moritz
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/frenzel
            089 - 289-15240
</pre>

Here everything related with the preprocess node
*/
/*----------------------------------------------------------------------*/
#ifdef D_EXODUS
#ifndef PRE_ELEMENT_H_
#define PRE_ELEMENT_H_

#include <string>
#include <iostream>
#include <vector>
#include <Teuchos_RefCountPtr.hpp>
#include "../drt_lib/drt_dserror.H"
#include "pre_exodus_reader.H"

using namespace std;
using namespace Teuchos;

namespace EXODUS
{

/*!
\class PreElement

\brief PreElement is an element class for preprocess handling of elements. It is 
adopted to DRT::Elements

\author maf (frenzel@lnm.mw.tum.de)
*/

// forward declaration
class PreNode;

class PreElement
{
public:
//  enum ShapeType{
//    quad4,
//    tri3,
//  };
  PreElement(int id, ElementBlock::Shape distype);
  
  PreElement(const PreElement& old);
  
  virtual ~PreElement();
  
  inline int Id() const {return id_; }
  
  inline ElementBlock::Shape GetShape() const {return distype_;}
  
  inline int NumNode() const {return nodeid_.size(); }
  
  inline const int* NodeIds() const {if (nodeid_.size()) return &nodeid_[0];
                                     else                return NULL;     }
  
  PreNode** Nodes() { if (node_.size()) return &node_[0]; else return NULL; }
  
  void Print(ostream& os) const;
  
  void SetNodeIds(const int nnode, const int* nodes);
  
private:
  int id_;
  
  ElementBlock::Shape distype_;
  
  vector<int> nodeid_;
  
  vector<PreNode*> node_;
  
};
} // of namespace EXODUS

#endif /*PRE_ELEMENT_H_*/
#endif // D_EXODUS

