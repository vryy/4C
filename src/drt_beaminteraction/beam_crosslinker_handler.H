/*----------------------------------------------------------------------*/
/*! \file

\brief Handler to control beam crosslinker simulations

\level 2

\maintainer Jonas Eichinger
*----------------------------------------------------------------------*/


#ifndef BEAM_CROSSLINKER_HANDLER_H
#define BEAM_CROSSLINKER_HANDLER_H

#include "../drt_binstrategy/binning_strategy.H"

// forward declarations
namespace DRT
{
  class Exporter;
  class Element;
  class Discretization;
  class Node;
}  // namespace DRT

namespace IO
{
  class DiscretizationWriter;
}

namespace LINALG
{
  class MapExtractor;
}

namespace BEAMINTERACTION
{
  class BeamCrosslinkerHandler
  {
   public:
    // constructor
    BeamCrosslinkerHandler();

    /// initialize linker handler
    void Init(int myrank, Teuchos::RCP<BINSTRATEGY::BinningStrategy> binstrategy);

    /// setup linker handler
    void Setup();

    // destructor
    virtual ~BeamCrosslinkerHandler(){};

    /// get binning strategy
    virtual inline Teuchos::RCP<BINSTRATEGY::BinningStrategy>& BinStrategy()
    {
      return binstrategy_;
    }

    virtual inline BINSTRATEGY::BinningStrategy const& BinStrategy() const { return *binstrategy_; }

    /// initial distribution of linker ( more general nodes of bindis_ ) to bins
    virtual void DistributeLinkerToBins(Teuchos::RCP<Epetra_Map> const& linkerrowmap);

    /// remove all linker
    virtual void RemoveAllLinker();

    /// get bin colume map
    virtual inline Teuchos::RCP<Epetra_Map>& BinColMap() { return bincolmap_; }

    /// get myrank
    virtual inline int MyRank() { return myrank_; }

    /// linker are checked whether they have moved out of their current bin
    /// and transferred if necessary
    virtual Teuchos::RCP<std::list<int>> TransferLinker(bool const fill_using_ghosting = true);

    /// node is placed into the correct row bin or put into the list of homeless linker
    virtual bool PlaceNodeCorrectly(Teuchos::RCP<DRT::Node> node,  ///< node to be placed
        const double* currpos,                              ///< current position of this node
        std::list<Teuchos::RCP<DRT::Node>>& homelesslinker  ///< list of homeless linker
    );

    /// round robin loop to fill linker into its correct bin on according proc
    virtual void FillLinkerIntoBinsRoundRobin(
        std::list<Teuchos::RCP<DRT::Node>>& homelesslinker  ///< list of homeless linker
    );

    /// get neighbouring bins of linker containing boundary row bins
    virtual void GetNeighbouringBinsOfLinkerContainingBoundaryRowBins(std::set<int>& colbins) const;

   protected:
    /// fill linker into their correct bin on according proc using remote id list
    virtual Teuchos::RCP<std::list<int>> FillLinkerIntoBinsRemoteIdList(
        std::list<Teuchos::RCP<DRT::Node>>& homelesslinker  ///< set of homeless linker
    );

    /// fill linker into their correct bin on according proc using one layer ghosting
    /// note, this is faster than the other two method as there is no communication required
    /// to find new owner ( complete one layer bin ghosting is required though)
    virtual Teuchos::RCP<std::list<int>> FillLinkerIntoBinsUsingGhosting(
        std::list<Teuchos::RCP<DRT::Node>>& homelesslinker  ///< set of homeless linker
    );

    /// receive linker and fill them in correct bin
    virtual void ReceiveLinkerAndFillThemInBins(int const numrec, DRT::Exporter& exporter,
        std::list<Teuchos::RCP<DRT::Node>>& homelesslinker);

   private:
    /// binning strategy
    Teuchos::RCP<BINSTRATEGY::BinningStrategy> binstrategy_;

    /// myrank
    int myrank_;

    /// colmap of bins
    Teuchos::RCP<Epetra_Map> bincolmap_;
  };

}  // namespace BEAMINTERACTION


/*----------------------------------------------------------------------*/
#endif
