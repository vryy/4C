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
#ifndef PRE_NODE_H_
#define PRE_NODE_H_

#include <string>
#include <iostream>
#include <vector>
#include <Teuchos_RefCountPtr.hpp>
#include "../drt_lib/drt_dserror.H"
//#include "pre_element.H"

using namespace std;
using namespace Teuchos;

namespace EXODUS
{

// forward declaration
class PreElement;

/*!
\class PreNode

\brief PreNode is a node class for preprocess handling of nodes. It is 
adopted to DRT::Node

\author maf (frenzel@lnm.mw.tum.de)
*/


class PreNode
{
public:
  PreNode(int id, const double* coords);
  
  PreNode(const PreNode& old);
  
  virtual ~PreNode();
  
  inline int ExoId() const {return id_+1;}

  inline int StoreId() const {return id_;}

  inline const double* X() const {return x_; }
  
  inline int NumElement() const {return element_.size(); }
  
  PreElement** Elements()   { if (NumElement()) return &element_[0]; else return NULL; }
  
  void Print(ostream& os, bool storeid=false) const;
  
  void AddElementPtr(PreElement* eleptr)
  { const int size = element_.size(); element_.resize(size+1);
    element_[size] = eleptr; return;}

private:
  int id_;
  
  double x_[3];
  
  vector<PreElement*> element_;
  
};

} // of namespace EXODUS

#endif /*PRE_NODE_H_*/
#endif // D_EXODUS

