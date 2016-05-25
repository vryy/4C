/*----------------------------------------------------------------------*/
/*!
\file binning_strategy.H

\brief Binning strategy for neighborhood search

\level 2

\maintainer Georg Hammerl
*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                              ghamm 09/12 |
 *----------------------------------------------------------------------*/
#ifndef BINNING_STRATEGY_H
#define BINNING_STRATEGY_H

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 09/12 |
 *----------------------------------------------------------------------*/
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_inpar/inpar_cavitation.H"

#include "Epetra_MpiComm.h"
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_CrsGraph.h>
#include <vector>

/*----------------------------------------------------------------------*
 | forward declarations                                     ghamm 09/12 |
 *----------------------------------------------------------------------*/
namespace DRT
{
  class Discretization;
  class Node;
}


/*----------------------------------------------------------------------*
 | binning strategy                                         ghamm 11/13 |
 *----------------------------------------------------------------------*/
/// BINSTRATEGY namespace
namespace BINSTRATEGY
{
  // forward declaration
  class Less;

  class BinningStrategy
  {
  public:

    /// constructor
    BinningStrategy(
      const Epetra_Comm& lcomm,  ///< local epetra communicator
      double cutoff_radius,   ///< radius of influence
      LINALG::Matrix<3,2> XAABB   ///< radius of influence
      );

    /// constructor
    BinningStrategy(
      const Epetra_Comm& lcomm  ///< local epetra communicator
      );

    /// repartitioning constructor
    BinningStrategy(std::vector<Teuchos::RCP<DRT::Discretization> > dis,
                    std::vector<Teuchos::RCP<Epetra_Map> >& stdelecolmap,
                    std::vector<Teuchos::RCP<Epetra_Map> >& stdnodecolmap);

    /// destructor
    virtual ~BinningStrategy(){};

    /// elements are assigned to bins (either slave or master side of mortar interface)
    virtual void DistributeElesToBins(
      const DRT::Discretization& mortardis,   ///< mortar interface discretization
      std::map<int, std::set<int> >& binelemap,   ///< map of bins and assigned elements
      bool isslave   ///< decide whether slave or master side is processed
      );

    /// elements from the underlying discretization are assigned to bins
    virtual void DistributeElesToBins(
      Teuchos::RCP<DRT::Discretization> underlyingdis,   ///< discret of which the elements are distributed to bins
      std::map<int, std::set<int> >& rowelesinbin   ///< map of bins and assigned row elements
      );

    /// nodes are assigned to bins
    void DistributeNodesToBins(
      Teuchos::RCP<DRT::Discretization> underlyingdis,
      std::map<int, std::vector<int> >& nodesinbin
      );

    /// elements/nodes from all provided discr. are distributed according to bin distrib.
    void WeightedRepartitioning(
      std::vector<Teuchos::RCP<DRT::Discretization> > dis,
      std::vector<Teuchos::RCP<Epetra_Map> >& stdelecolmap,
      std::vector<Teuchos::RCP<Epetra_Map> >& stdnodecolmap
      );

    /// element map with extended ghosting for master side of mortar problem is returned
    virtual Teuchos::RCP<Epetra_Map> ExtendGhosting(
      DRT::Discretization& mortardis,   ///< mortar interface discretization
      Teuchos::RCP<Epetra_Map> initial_elecolmap,   ///< initial element col map for standard ghosting
      std::map<int, std::set<int> >& slavebinelemap,    /// map of bins and assigned slave elements
      std::map<int, std::set<int> >& masterbinelemap   /// map of bins and assigned master elements
      );

    /// extend ghosting according to scatra elements
    virtual void ExtendGhosting(
      Teuchos::RCP<DRT::Discretization> scatradis,       ///< scatra discretization
      std::map<int, std::set<int> >& escapedpartelemap,  ///< map of escaped particles and corresponding scatra elements
      std::map<int, std::set<int> >& myescapedpartelemap ///< map of escaped particles and corresponding scatra elements including particles from other procs
      );

    /// extend ghosting according to bin distribution
    Teuchos::RCP<Epetra_Map> ExtendGhosting(
      const Epetra_Map* initial_elecolmap,
      std::map<int, std::set<int> >& binelemap,
      std::map<int, std::set<int> >& extendedghosting,
      Teuchos::RCP<Epetra_Map> bincolmap = Teuchos::null
      );

    virtual void ExtendGhosting(
        std::vector<Teuchos::RCP<DRT::Discretization> > dis
      );

    /// collect information of content in bins from other procs via Round Robin loop
    virtual void CollectInformation(
      Teuchos::RCP<Epetra_Map> rowbins,
      std::map<int, std::vector<int> >& nodesinbin,
      std::map<int, std::vector<int> >& mynodesinbin
      );

    /// return flag whether pbcs are applied to the problem
    virtual bool HavePBCs(const int dim = -1) { if(dim == -1) return havepbc_; else return pbconoff_[dim]; };

    /// return delta for pbc in x, y, or z direction
    virtual double PBCDelta(const int dim) { return pbcdeltas_[dim]; };

    /// return dimension of the problem
    virtual INPAR::PARTICLE::ParticleDim ParticleDim() { return particle_dim_; };

    /// convert i,j,k into bin id
    virtual int ConvertijkToGid(
      int* ijk  ///< i,j,k to be converted into a bin id
      );

    /// get all bin ids for given range of ijk
    virtual void GidsInijkRange(
      const int* ijk_bound,   ///< given range of ijk
      std::set<int>& binIds,   ///< all bin ids in the specified range
      bool checkexistence   ///< check can be added whether the gids are on this proc
      );

    /// get all bin ids for given range of ijk
    virtual void GidsInijkRange(
      const int* ijk_bound,   ///< given range of ijk
      std::vector<int>& binIds,   ///< all bin ids in the specified range
      bool checkexistence   ///< check can be added whether the gids are on this proc
      );

    ///convert bin id into i,j,k
    virtual void ConvertGidToijk(
      int gid,  ///< bin id to be converted into i,j,k
      int* ijk  ///< resulting i,j,k
      );

    /// convert position into bin id
    virtual int ConvertPosToGid(
      const std::vector<double>& pos  ///< position of which the corresponding bin id is asked for
      );

    /// convert position into bin id
    virtual int ConvertPosToGid(
      const double* pos   ///< position of which the corresponding bin id is asked for
      );

    ///convert position into i,j,k
    virtual void ConvertPosToijk(
      const double* pos,   ///< position to be converted into i,j,k
      int* ijk   ///< resulting i,j,k
      );

    /// convert position into bin id
    virtual int ConvertPosToGid(
      const LINALG::Matrix<3,1>& pos   ///< position of which the corresponding bin id is asked for
      );

    ///convert position into i,j,k
    virtual void ConvertPosToijk(
      const LINALG::Matrix<3,1>& pos,   ///< position to be converted into i,j,k
      int* ijk  ///< resulting i,j,k
      );

    /// get 26 neighboring bin ids to binId (if existing)
    virtual void GetBinConnectivity(
      const int binId,   ///< bin id whose connectivity is asked for
      std::vector<int>& binIds   ///< all neighboring bins on axes
      );

    /// get nodal coordinates of bin with given id
    virtual void GetBinCorners(
      const int binId,   ///< bin id of which corners are calculated
      std::vector<LINALG::Matrix<3,1> >& bincorners   ///< corner position of given bin
      );

    /// centroid position for given bin id is returned
    virtual LINALG::Matrix<3,1> GetBinCentroid(
      const int binId   ///< bin id of which centroid is calculated
      );

    virtual void RevertExtendedGhosting(std::vector<Teuchos::RCP<DRT::Discretization> > dis,
        std::vector<Teuchos::RCP<Epetra_Map> >& stdelecolmap,
        std::vector<Teuchos::RCP<Epetra_Map> >& stdnodecolmap);

  protected:

    /// bins are created based on XAABB and cutoff radius
    virtual void CreateBins(
      Teuchos::RCP<DRT::Discretization> dis   ///< underlying discret which can be used to compute bounding box
      );

    // create XAABB for discretization
    virtual void CreateXAABB(Teuchos::RCP<DRT::Discretization> dis);

    // create XAABB for discretization and compute cutoff with largest element in discret
    virtual void CreateXAABB(Teuchos::RCP<DRT::Discretization> dis, LINALG::Matrix<3,2>& XAABB, double& cutoff);

    /// build periodic boundary conditions
    virtual void BuildParticlePeriodicBC();

    Teuchos::RCP<DRT::Discretization> particledis_;     ///< particle discretization with bins

    double cutoff_radius_;    ///< maximum influence radius

    LINALG::Matrix<3,2> XAABB_;    ///< extendedAxisAlignedBoundingBox of bins

    double bin_size_[3];    ///< size of each bin in Cartesian coordinates

    double inv_bin_size_[3];   ///< inverse of size of each bin in Cartesian coordinates

    int bin_per_dir_[3];    ///< number of bins per direction

    double inv_bin_per_dir_[3];    ///< inverse of number of bins per direction

    bool havepbc_;   ///< flag whether periodic boundary conditions are specified for particles

    bool pbconoff_[3];   ///< flags for existence of pbcs in x, y, z direction

    double pbcdeltas_[3];   ///< deltas for pbcs in x, y, z direction

    const INPAR::PARTICLE::ParticleDim particle_dim_;   ///< dimension of problem

    const bool sparse_binning_;   ///< flag whether bin scheme is sparse

    //! @name miscellaneous stuff
    //@{

    const int myrank_;    ///< myrank

    //@}

  }; // BinningStrategy

  // class for comparing Teuchos::RCP<DRT::Node> in a std::set
  class Less
  {
  public:
    bool operator()(const Teuchos::RCP<const DRT::Node>& first, const Teuchos::RCP<const DRT::Node>& second) const;
  };

} // namespace BINSTRATEGY


/*----------------------------------------------------------------------*/
#endif  // BINNING_STRATEGY_H
