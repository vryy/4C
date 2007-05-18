/*!----------------------------------------------------------------------
\file
\brief connect nodes from two nodesets by their distance


<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef NODEMATCHINGOCTREE_H
#define NODEMATCHINGOCTREE_H

#include "drt_utils.H"
#include "drt_discret.H"
#include "drt_parobject.H"
#include "drt_exporter.H"
#include "Teuchos_RefCountPtr.hpp"

extern "C"
{
#include "../headers/standardtypes.h"
}
class OctreeElement;


class NodeMatchingOctree
{
public:
  // @{ \name Constructors and destructors

  /*!
  \brief Standard Constructor

  */

  NodeMatchingOctree(
    RefCountPtr<DRT::Discretization> actdis,
    const vector <int> &             masternodeids,
    int                              maxnodeperleaf,
    double                           tol
    );

  /*!
  \brief Destructor

  */
  virtual ~NodeMatchingOctree();

   //@}

  void CreateGlobalNodeMatching(
    const vector <int> slavenodeids,
    vector <int>       dofsforpbcplane,
    map<int,int> &     sidtomid
    );
  
  bool SearchClosestNodeOnThisProc(
    vector <double> x,
    int           & idofclosestpoint,
    double        & distofclosestpoint
  );
  
protected:
private:
  RefCountPtr<DRT::Discretization> discret_;
  RefCountPtr<OctreeElement>       octreeroot_;
  vector <double>                  masterplanecoords_;
}; // class NodeMatchingOctree





class OctreeElement
{

public:
  // @{ \name Constructors and destructors

  /*!
  \brief Standard Constructor

  */

  OctreeElement(
    RefCountPtr<DRT::Discretization> actdis,
    vector <int> &                   nodeids,
    Epetra_SerialDenseMatrix&        boundingbox,
    int                              layer,
    int                              maxnodeperleaf,
    double                           tol
    );

  /*!
  \brief Destructor

  */
  virtual ~OctreeElement();

  //@}

  bool IsPointInBoundingBox(
    vector <double> &x
    );

  bool IsLeaf();
  
  RefCountPtr<OctreeElement> ReturnChildContainingPoint(
    vector <double> &x
    );
  
  void SearchClosestNodeInLeaf(
    vector <double> & x,
    int             & idofclosestpoint,
    double          & distofclosestpoint
  );

  
  void Print(ostream& os) const;
protected:

private:
  RefCountPtr<DRT::Discretization> discret_;
  Epetra_SerialDenseMatrix         boundingbox_;
  vector <int>                     nodeids_;
  int                              layer_;
  RefCountPtr<OctreeElement>       octreechild1_;
  RefCountPtr<OctreeElement>       octreechild2_;
  
}; // class OctreeElement


#endif  // not NODEMATCHINGOCTREE_H

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
