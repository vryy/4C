/*----------------------------------------------------------------------*/
/*!
\file pre_exodus_soshextrusion.H

\brief solid-shell body creation by extruding surface 

<pre>
Maintainer: Moritz
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/frenzel
            089 - 289-15240
</pre>

Here everything related with solid-shell body extrusion
*/
/*----------------------------------------------------------------------*/
#ifdef EXODUS
#ifndef PRE_EXODUS_SOSHEXTRUSION_H_
#define PRE_EXODUS_SOSHEXTRUSION_H_

#include <string>
#include <iostream>
#include <vector>
#include <Teuchos_RefCountPtr.hpp>
#include "../drt_lib/drt_dserror.H"
#include "pre_exodus_reader.H"

extern "C" /* stuff which is c and is accessed from c++ */
{
#include <exodusII.h>
}

using namespace std;
using namespace Teuchos;

/*!
\class Soshextrusion

\brief Soshextrusion handles the extrusion of a solid-shell mesh body based
on a surface mesh


\author maf (frenzel@lnm.mw.tum.de)
*/
class Soshextrusion : public Mesh
{
public:
  //! constructor
  Soshextrusion(string exofilename, double thickness, int layers);
  
  //! destructor
  virtual ~Soshextrusion();
};

// forward declaration
class Node;

class ShellElement
{
public:
  enum ShellShape{
    quad4,
    tri3,
  };
  ShellElement(int id, ShellShape shape);
  
  ShellElement(const ShellElement& old);
  
  virtual ~ShellElement();
  
  inline int Id() const {return id_; }
  
  inline ShellShape Shape() const {return shellshape_;}
  
  inline int NumNode() const {return nodeid_.size(); }
  
  inline const int* NodeIds() const {if (nodeid_.size()) return &nodeid_[0];
                                     else                return NULL;     }
  
  Node** Nodes() { if (node_.size()) return &node_[0]; else return NULL; }
  
  void Print(ostream& os) const;
  
private:
  int id_;
  
  ShellShape shellshape_;
  
  vector<int> nodeid_;
  
  vector<Node*> node_;
  
};

class Node
{
  Node(int id, const double* coords);
  
  Node(const Node& old);
  
  virtual ~Node();
  
  inline int Id() const {return id_;}
  
  inline const double* X() const {return x_; }
  
  inline int NumElement() const {return element_.size(); }
  
  ShellElement** Elements()   { if (NumElement()) return &element_[0]; else return NULL; }
  
  void Print(ostream& os) const;
  
  void AddElementPtr(ShellElement* eleptr)
  { const int size = element_.size(); element_.resize(size+1);
    element_[size] = eleptr; return;}

private:
  int id_;
  
  double x_[3];
  
  vector<ShellElement*> element_;
  
};




#endif /*PRE_EXODUS_SOSHEXTRUSION_H_*/
#endif /*EXODUS*/

