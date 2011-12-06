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

    std::vector<int> InWhichOctantLies(const int& thisBBoxID);
    // intersection test of all elements in the octant in which a given bounding box lies
    bool IsecBBoxesOfOctantWith(const int& bboxid, Epetra_SerialDenseMatrix& nodecoords);

  private:
    // Bounding Box Types
    enum BboxType
    {
    	none,
      axisaligned,
      cyloriented
    };

    // generator of extended Bounding Boxes (axis aligned as well as cylindrical oriented)
    void CreateBoundingBoxes(std::map<int, LINALG::Matrix<3,1> >&  currentpositions);
    // get the dimensions of the root octant
    LINALG::Matrix<1,6> GetRootOctant();
    // create axis aligned bounding boxes
    void CreateAABB(const Epetra_SerialDenseMatrix& coord, const int& elecolid);
    // create coylindrical oriented bounding boxes
    void CreateCOBB(const Epetra_SerialDenseMatrix&  coord, const int& elecolid);
    // create axis aligned bounding boxes
    RCP<Epetra_SerialDenseMatrix> CreateAABBNoElement(const Epetra_SerialDenseMatrix& coord);
    // create coylindrical oriented bounding boxes
    RCP<Epetra_SerialDenseMatrix> CreateCOBBNoElement(const Epetra_SerialDenseMatrix&  coord);

    // spatial division of 3-dimensional volume
    void locateAll(LINALG::Matrix<1,6>& rootoctantlim, std::vector<std::vector<int> >& bboxesinoctants);

    // primitive for locateAll
    void locateBox(std::vector<std::vector<double> > allBBoxesstdvec,
                   LINALG::Matrix<1,6> lim,
                   std::vector<std::vector<double> >* OctreeLimits,
                   std::vector<std::vector<int> >& bboxesinoctants,
                   int& treedepth);

    // build list
    void BoundingBoxIntersection(std::map<int, LINALG::Matrix<3,1> >&  currentpositions, std::vector<std::vector<int> >& bboxesinoctants, vector<RCP<Beam3contact> >* contactpairs);

    // intersection method applying axis-aligned bounding boxes when both boxes belong to existing elements
    bool IntersectionAABB(std::vector<int>& bboxIDs);
    // intersection method applying cylindrical oriented bounding boxes when both boxes belong to existing elements
    bool IntersectionCOBB(std::vector<int>& bboxIDs);
    // intersection method applying axis-aligned bounding boxes when only one of the boxes belongs to an existing element
    bool IsecAABBNoElement(const int& bboxinoct, RCP<Epetra_SerialDenseMatrix> bboxlimits);
    // intersection method applying cylindrical oriented bounding boxes when only one of the boxes belongs to an existing element
    bool IsecCOBBNoElement(const int& bboxinoct, RCP<Epetra_SerialDenseMatrix> bboxlimits);
    
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

    // stores the IDs and the coordinates of all bounding boxes
    RCP<Epetra_MultiVector> allbboxes_;

    // vector listing the bounding boxes located in the octants
    RCP<Epetra_MultiVector> bboxesinoctants_;

    // mapping bounding boxes to octants they lie in
    RCP<Epetra_MultiVector> bbox2octant_;

    // vector holding information on how many times a bounding box is shifted due to periodic boundary conditions
    RCP<Epetra_Vector> numshifts_;

    // Bounding Box type
    Beam3ContactOctTree::BboxType boundingbox_;
};

#endif  // #ifndef BEAM3OCTTREE_H
#endif  // #ifdef CCADISCRET
