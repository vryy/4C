/*!----------------------------------------------------------------------
\file beam3contact_octtree.H

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef BEAM3OCTTREE_H
#define BEAM3OCTTREE_H

#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_lib/drt_discret.H"
#include "beam3contact.H"

using namespace std;
using namespace Teuchos;
using namespace CONTACT;

/*!
 \brief Octtree for beam contact search...
 Refer also to the Semesterarbeit of Christian Roth, 2011
*/
class Beam3ContactOctTree
{
  public:
    //!\brief Constructor
    Beam3ContactOctTree(ParameterList& params,
												DRT::Discretization& discret,
                        DRT::Discretization& searchdis,
                        const int& dofoffset);
    
    //!\brief Destructor
    virtual ~Beam3ContactOctTree() {};

    // call octtree search routine
    vector<RCP<Beam3contact> > OctTreeSearch(std::map<int, LINALG::Matrix<3,1> >&  currentpositions);

  private:
    // Bounding Box Types
    enum BboxType
    {
    	none,
      axisaligned,
      cyloriented
    };

    // generator of extended Bounding Boxes (axis aligned as well as cylindrical oriented)
    void CreateBoundingBoxes(std::map<int, LINALG::Matrix<3,1> >&  currentpositions,RCP<Epetra_MultiVector> allBBoxes);
    // create axis aligned bounding boxes
    void CreateAABB(const Epetra_SerialDenseMatrix& coord, RCP<Epetra_MultiVector> allBBoxes, const int& elecolid, const double& PeriodLength);
    // create coylindrical oriented bounding boxes
    void CreateCOBB(const Epetra_SerialDenseMatrix&  coord, RCP<Epetra_MultiVector> allBBoxes, const int& elecolid, const double& PeriodLength);

    // spatial division of 3-dimensional volume
    void locateAll(RCP<Epetra_MultiVector> allBBoxes,
                   std::vector<std::vector<int> >& bboxesinoctants);

    // primitive for locateAll
    void locateBox(std::vector<std::vector<double> > allBBoxesstdvec,
                   LINALG::Matrix<1,6> lim,
                   std::vector<std::vector<double> >* OctreeLimits,
                   std::vector<std::vector<int> >& bboxesinoctants,
                   int& treedepth);

    // build list
    void BoundingBoxIntersection(std::map<int, LINALG::Matrix<3,1> >&  currentpositions,
                          vector<RCP<Beam3contact> >* contactpairs,
                          std::vector<std::vector<int> >* bboxesinoctants,
                          RCP<Epetra_MultiVector> allBBoxes);

    // intersection method applying axis-aligned bounding boxes
    bool IntersectionAABB(std::vector<int>& bboxIDs, RCP<Epetra_MultiVector> allBBoxes);
    // intersection method applying cylindrical oriented bounding boxes
    bool IntersectionCOBB(std::vector<int>& bboxIDs, RCP<Epetra_MultiVector> allBBoxes);
    
    bool periodicBC_;

    // problem discretization
    DRT::Discretization& discret_;
    
    // contact discretization
    DRT::Discretization& searchdis_;
    
    // number of initial nodes
    int basisnodes_;
    
    // dof offset between discret_ and searchdis_
    int dofoffset_;
    
    // max tree depth
    int maxtreedepth_;

    // diameters of all beam elements
    RCP<Epetra_Vector> diameter_;

    // redundant map of bounding boxes
    RCP<Epetra_Map> AABBmap_;

    // Bounding Box type
    Beam3ContactOctTree::BboxType boundingbox_;
};

#endif  // #ifndef BEAM3OCTTREE_H
#endif  // #ifdef CCADISCRET
