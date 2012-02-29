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
#include "beam3contact.H"

using namespace std;
using namespace Teuchos;
using namespace CONTACT;

namespace DRT
{
  class Discretization;
}

namespace LINALG
{
  class SparseMatrix;
}

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

    //!\brief call octtree search routine
    vector<RCP<Beam3contact> > OctTreeSearch(std::map<int, LINALG::Matrix<3,1> >&  currentpositions, int step=-1);

    //!\brief checks in which octant a given bounding box lies
    std::vector<int> InWhichOctantLies(const int& thisBBoxID);

    //!\brief intersection test of all elements in the octant in which a given bounding box lies
    bool IntersectBBoxesWith(Epetra_SerialDenseMatrix& nodecoords, Epetra_SerialDenseMatrix& nodeLID);

    //!\brief output of octree discretization, bounding boxes and contact pairs
    void OctreeOutput(std::vector<RCP<Beam3contact> >& cpairs, int step);

  private:
    // Bounding Box Types
    enum BboxType
    {
    	none,
      axisaligned,
      cyloriented,
      spherical
    };

    //!\brief Initialize class vectors for new Octree search
    void InitializeOctreeSearch();
    //!\brief generator of extended Bounding Boxes (axis aligned as well as cylindrical oriented)
    void CreateBoundingBoxes(std::map<int, LINALG::Matrix<3,1> >&  currentpositions);
    //!\brief get the dimensions of the root octant
    LINALG::Matrix<6,1> GetRootBox();
    //!\brief create axis aligned bounding boxes
    void CreateAABB(Epetra_SerialDenseMatrix& coord, const int& elecolid, RCP<Epetra_SerialDenseMatrix> bboxlimits=Teuchos::null);
    //!\brief create coylindrical oriented bounding boxes
    void CreateCOBB(Epetra_SerialDenseMatrix& coord, const int& elecolid, RCP<Epetra_SerialDenseMatrix> bboxlimits=Teuchos::null);
    // create spherical bounding boxes for crosslinker
    void CreateSPBB(Epetra_SerialDenseMatrix& coord, const int& elecolid, RCP<Epetra_SerialDenseMatrix> bboxlimits=Teuchos::null);


    //!\brief spatial division of 3-dimensional volume. Returns false if no bounding boxes were found in any octant
    bool locateAll();

    //!\brief primitive for locateAll
    void locateBox(std::vector<std::vector<double> >& allbboxesstdvec,
                   LINALG::Matrix<6,1>& lim,
                   std::vector<LINALG::Matrix<6,1> >& OctreeLimits,
                   std::vector<std::vector<int> >& bboxesinoctants,
                   int& treedepth);

    //!\brief build list
    void BoundingBoxIntersection(std::map<int, LINALG::Matrix<3,1> >&  currentpositions, vector<RCP<Beam3contact> >* contactpairs);

    //!\brief intersection method applying axis-aligned bounding boxes when both boxes belong to existing elements
    bool IntersectionAABB(const std::vector<int>& bboxIDs, RCP<Epetra_SerialDenseMatrix> bboxlimits=Teuchos::null);
    //!\brief intersection method applying cylindrical oriented bounding boxes when both boxes belong to existing elements
    bool IntersectionCOBB(const std::vector<int>& bboxIDs, RCP<Epetra_SerialDenseMatrix> bboxlimits=Teuchos::null);
    // intersection method applying spherical bounding boxes for crosslinker when both boxes belong to existing elements
    bool IntersectionSPBB(const std::vector<int>& bboxIDs, RCP<Epetra_SerialDenseMatrix> bboxlimits=Teuchos::null);
    
    //! \brief comunicate Vector to all participating processors
    void CommunicateVector(Epetra_Vector& InVec, Epetra_Vector& OutVec, bool zerofy=false, bool doexport=true, bool doimport=true);
    //! \brief comunicate MultiVector to all participating processors
    void CommunicateMultiVector(Epetra_MultiVector& InVec, Epetra_MultiVector& OutVec, bool zerofy=false, bool doexport=true, bool doimport=true);

    //!\brief flag indicating the use of periodic boundary conditions
    bool periodicBC_;

    //!\brief vector holding the edge lengths of the cuboid periodic volume
    Teuchos::RCP<std::vector<double> > periodlength_;

    //!\brief Matrix holding the spatial limits of the root box
    LINALG::Matrix<6,1> rootbox_;

    //!\brief Matrix holding the spatial limits of the root box in reference configuration
    LINALG::Matrix<6,1> initbox_;

    //!\brief problem discretization
    DRT::Discretization& discret_;
    
    //!\brief contact discretization
    DRT::Discretization& searchdis_;
    
    //!\brief number of initial nodes
    int basisnodes_;
    
    //!\brief dof offset between discret_ and searchdis_
    int dofoffset_;

    //!\brief maximum tree depth
    int maxtreedepth_;

    //!\brief minimum number of BBs per octant
    int minbboxesinoctant_;

    //!\brief scalar factor for diameters
    double extrusionfactor_;

    //!\brief scalar factor for radial scaling of the bounding box
    double radialextrusion_;

    //!\brief diameters of all beam elements
    RCP<Epetra_Vector> diameter_;

    //!\brief stores the IDs and the coordinates of all bounding boxes
    RCP<Epetra_MultiVector> allbboxes_;

    //!\brief vector listing the bounding boxes located in the octants
    RCP<Epetra_MultiVector> bboxesinoctants_;

    //!\brief mapping bounding boxes to octants they lie in
    RCP<Epetra_MultiVector> bbox2octant_;

    //!\brief storage vector for octree octant limits
    std::vector<LINALG::Matrix<6,1> > octreelimits_;

    //!\brief vector holding information on how many times a bounding box is shifted due to periodic boundary conditions
    RCP<Epetra_Vector> numshifts_;

    //!\brief vector holding information to which line (design line) a bounding box belongs
    RCP<Epetra_Vector> bbox2line_;

    //!\brief Bounding Box type
    Beam3ContactOctTree::BboxType boundingbox_;
};

#endif  // #ifndef BEAM3OCTTREE_H
#endif  // #ifdef CCADISCRET
