/*!----------------------------------------------------------------------
\file statmech_search.H

<pre>
Maintainer: Kei Müller
            mueller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>

*----------------------------------------------------------------------*/


#include <Teuchos_RCP.hpp>
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"
#include "../linalg/linalg_fixedsizematrix.H"

namespace DRT
{
  class Discretization;
}

namespace LINALG
{
  class SparseMatrix;
}

namespace STATMECH
{
namespace SEARCH
{
// object for octree-type search
class Octree
{
	public:
	//! \brief constructor
  Octree(Teuchos::RCP<std::vector<double> > periodlength,
         Teuchos::RCP<DRT::Discretization> discret,
         Teuchos::RCP<Epetra_Map> bspotrowmap,
         Teuchos::RCP<Epetra_Map> bspotcolmap,
         int maxtreedepth,
         int minbboxesperoct,
         double bindingradius);

  //!\brief Destructor
  virtual ~Octree() {};

  class OctreeNode
  {
    private:
    int OctantID;
    LINALG::Matrix<6,1> limits;
    public:
    //Constructor
    OctreeNode(LINALG::Matrix<6,1>& inlimits); //
    //Destructor
    virtual ~OctreeNode() {};

    std::vector<OctreeNode*>  children;
    int HaveChildren();
    void AddChild(OctreeNode* new_node);
    //BspotOctreeNode* GetChild(int i);
    void SetOctantID(int id);
    int GetOctantID();
    LINALG::Matrix<6,1> LimitsOfNode();
  };

  //! \brief return row map format vector holding all bounding boxes in the octants
  const Teuchos::RCP<Epetra_MultiVector> BBoxesInOctRow() const{return bboxinoctrow_;};
  //! \brief return row map format vector holding all crosslinkers in the octants
  const Teuchos::RCP<Epetra_MultiVector> CrosslinkerInOctants() const{return crosslinkerinoctants_;};

  //! \brief Build Octree with spherical bounding boxes
  void BuildOctree(const Teuchos::RCP<Epetra_MultiVector>& bspotpositions);
  //! \brief Locate the Oktant of the Octree each position belongs to
  void LocatePositions(Teuchos::RCP<Epetra_MultiVector>& positions, Teuchos::RCP<Epetra_Map> transfermap);

  //! \brief Locate bounding boxes within octree
  void LocateBoundingBox(std::vector<std::vector<double> >& allbboxesstdvec,
                         LINALG::Matrix<6,1>&               lim,
                         std::vector<LINALG::Matrix<6,1> >& OctreeLimits,
                         std::vector<std::vector<int> >&    bboxesinoctants,
                         int&                               treedepth,
                         OctreeNode&                        parentNode);

private:
  //!\brief map of binding spots, row format
  Teuchos::RCP<Epetra_Map> bspotrowmap_;
  //!\brief map of binding spots, column format
  Teuchos::RCP<Epetra_Map> bspotcolmap_;
  //!\brief pointer to discretization
  Teuchos::RCP<DRT::Discretization> discret_;
  //! \brief comunicate MultiVector to all participating processors
  void CommunicateMultiVector(Teuchos::RCP<Epetra_MultiVector> InVec,
                              Teuchos::RCP<Epetra_MultiVector> OutVec,
                              bool                             doexport=true,
                              bool                             doimport=true,
                              bool                             zerofy=true,
                              bool                             exportinsert=false);

  //!\brief stores the IDs and the coordinates of all bounding boxes
  Teuchos::RCP<Epetra_MultiVector> allbboxes_;
  //!\brief vector listing the bounding boxes located in the octants
  Teuchos::RCP<Epetra_MultiVector> bboxesinoctants_;
  //!\brief vector listing the bounding boxes located in the octants in rowmapformat (necessary for parallel binding partner detection)
  Teuchos::RCP<Epetra_MultiVector> bboxinoctrow_;
  //!\brief storage vector for octree octant limits
  std::vector<LINALG::Matrix<6,1> > octreelimits_;
  //!\brief length of bounding boxes around the binding spot (cubes)
  double bsr_;
  //!\brief maximum tree depth
  int maxtreedepth_;
  //!\brief minimum number of BBs per octants
  int minbboxesinoctant_;
  //!\brief flag indicating the use of periodic boundary conditions
  bool periodicBC_;
  //!\brief mapping bounding boxes to octants they lie in (not needed for the search algorithm as it is implemented in statmechmanager atm)
  Teuchos::RCP<Epetra_MultiVector> bbox2octant_;
  //!\brief listing all crosslinkers in each oktant
  Teuchos::RCP<Epetra_MultiVector> crosslinkerinoctants_;
  //! \brief vector holding the dimensions of the periodic box
  Teuchos::RCP<std::vector<double> > periodlength_;
  //! \brief rootnode of the tree of RCP pointers that allow to reconstruct thedimensions of the octree. This is crucial in order to determine the Octant that belongs to a position (e.g. a crosslinker position)
  Octree::OctreeNode* RootNode_;
};

// Binning strategy search
class BinSearch
{
  public:
    //! \brief constructor
    BinSearch(Teuchos::RCP<DRT::Discretization> discret,
              const std::vector<double>&        rootboxdim,
              const std::vector<int>&           numbins);

    //!\brief Destructor
    virtual ~BinSearch() {};

    class Bin
    {
      public:
        Bin(const int& binid);
        ~Bin() {};
        const int& GetBinID() const {return binid_;}
        const std::vector<int>& GetBinMembers() const {return binmembers_;}
        void AddMember(const int& memberid);
      private:
        int binid_;
        std::vector<int> binmembers_;
    };

    //!\brief assign given positions to bins
    void AssignPositionsToBins(const Teuchos::RCP<Epetra_MultiVector>& positions);

    const std::vector<int>& GetResolution() const { return resolution_;}

    const LINALG::Matrix<6,1>& GetRootLimits() const { return rootlimits_;}

    const Bin& GetBin(const int& binid) const {return bins_.at(binfullmap_->LID(binid));}

    const int CalculateBinId(const std::vector<int>& indices);

    const std::vector<int> GetSurroundingBins(const std::vector<int> indices);

    bool HaveBin(const int& binid);

  private:

    void SetupRootBin(const std::vector<double>& rootboxdim);

    void BuildBinMaps();

    //! \brief comunicate MultiVector to all participating processors
    void CommunicateMultiVector(Teuchos::RCP<Epetra_MultiVector> InVec,
                                Teuchos::RCP<Epetra_MultiVector> OutVec,
                                bool                             doexport=true,
                                bool                             doimport=true,
                                bool                             zerofy=true,
                                bool                             exportinsert=false);

    Teuchos::RCP<DRT::Discretization> discret_;
    std::vector<int> resolution_;
    Teuchos::RCP<Epetra_Map> binrowmap_;
    Teuchos::RCP<Epetra_Map> binfullmap_;
    LINALG::Matrix<6,1> rootlimits_;

    std::vector<Bin> bins_;
};
}
}
