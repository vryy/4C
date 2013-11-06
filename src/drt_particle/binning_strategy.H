/*----------------------------------------------------------------------*/
/*!
\file binning_strategy.H

\brief Binning strategy for neighborhood search

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                              ghamm 09/12 |
 *----------------------------------------------------------------------*/
#ifndef BINNING_STRATEGY_H
#define BINNING_STRATEGY_H

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 09/12 |
 *----------------------------------------------------------------------*/
#include "../drt_adapter/adapter_algorithmbase.H"
#include "../linalg/linalg_fixedsizematrix.H"

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

    BinningStrategy(
      const Epetra_Comm& lcomm,  ///< local epetra communicator
      double cutoff_radius,   ///< radius of influence
      LINALG::Matrix<3,2> XAABB   ///< radius of influence
      );

    /// elements are assigned to bins (either slave or master side of mortar interface)
    void DistributeElesToBins(
      DRT::Discretization& mortardis,   ///< mortar interface discretization
      std::map<int, std::set<int> >& binelemap,   /// map of bins and assigned elements
      bool isslave   ///< decide whether slave or master side is processed
      );

    /// element map with extended ghosting for master side of mortar problem is returned
    Teuchos::RCP<Epetra_Map> ExtendGhosting(
      DRT::Discretization& mortardis,   ///< mortar interface discretization
      Teuchos::RCP<Epetra_Map> initial_elecolmap,   ///< initial element col map for standard ghosting
      std::map<int, std::set<int> >& slavebinelemap,    /// map of bins and assigned slave elements
      std::map<int, std::set<int> >& masterbinelemap   /// map of bins and assigned master elements
      );

    virtual ~BinningStrategy(){};

  protected:

    // bins are created based on XAABB and cutoff radius
    void CreateBins(Teuchos::RCP<DRT::Discretization> dis);

    /// convert i,j,k into bin id
    virtual int ConvertijkToGid(
      int* ijk  ///< i,j,k to be converted into a bin id
      );

    /// get all bin ids for given range of ijk
    void GidsInijkRange(
      int* ijk_bound,   ///< given range of ijk
      std::set<int>& binIds,   ///< all bin ids in the specified range
      DRT::Discretization& dis,   /// underlying discretization
      bool checkexistence   ///< check can be added whether the gids are on this proc
      );

    ///convert bin id into i,j,k
    virtual void ConvertGidToijk(
      int gid,  ///< bin id to be converted into i,j,k
      int* ijk  ///< resulting i,j,k
      );

    /// convert position into bin id
    int ConvertPosToGid(
      const std::vector<double>& pos  ///< position of which the corresponding bin id is asked for
      );

    /// convert position into bin id
    int ConvertPosToGid(
      const double* pos   ///< position of which the corresponding bin id is asked for
      );

    ///convert position into i,j,k
    void ConvertPosToijk(
      const double* pos,   ///< position to be converted into i,j,k
      int* ijk   ///< resulting i,j,k
      );

    /// convert position into bin id
    int ConvertPosToGid(
      const LINALG::Matrix<3,1> pos   ///< position of which the corresponding bin id is asked for
      );

    ///convert position into i,j,k
    void ConvertPosToijk(
      const LINALG::Matrix<3,1> pos,   ///< position to be converted into i,j,k
      int* ijk  ///< resulting i,j,k
      );

    /// get 26 neighboring bin ids to binId (if existing)
    void GetBinConnectivity(
      int binId,   ///< bin id whose connectivity is asked for
      std::vector<int>& binIds   ///< all neighboring bins on axes
      );  

    /// get nodal coordinates of bin with given id
    void GetBinCorners(
      int binId,   ///< bin id of which corners are calculated
      std::vector<LINALG::Matrix<3,1> >& bincorners   ///< corner position of given bin
      );

    /// centroid position for given bin id is returned
    LINALG::Matrix<3,1> GetBinCentroid(
      int binId   ///< bin id of which centroid is calculated
      );

    //! @name binning strategy related stuff
    //@{

    double cutoff_radius_;    ///< maximum influence radius

    LINALG::Matrix<3,2> XAABB_;    ///< extendedAxisAlignedBoundingBox of bins

    double bin_size_[3];    ///< size of each bin in Cartesian coordinates

    int bin_per_dir_[3];    ///< number of bins per direction

    //@}

    //! @name miscellaneous stuff
    //@{

    int myrank_;    ///< myrank

    //@}

  }; // BinningStrategy

  // class for comparing Teuchos::RCP<DRT::Node> in a std::set
  class Less
  {
  public:
    bool operator()(const Teuchos::RCP<const DRT::Node>& first, const Teuchos::RCP<const DRT::Node>& second) const;
  };

} // namespace PARTICLE


/*----------------------------------------------------------------------*/
#endif  // PARTICLE_ALGORITHM_H
