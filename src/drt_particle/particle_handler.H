/*----------------------------------------------------------------------*/
/*!
\file particle_handler.H

\brief Handler to control particle simulations

\level 2

\maintainer Georg Hammerl
*----------------------------------------------------------------------*/


#ifndef PARTICLE_HANDLER_H
#define PARTICLE_HANDLER_H

#include "binning_strategy.H"

#include <Epetra_MpiComm.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_CrsGraph.h>
#include <vector>


namespace DRT
{
  class Element;
  class Discretization;
  class Node;
}
namespace IO
{
  class DiscretizationWriter;
}
namespace LINALG
{
  class MapExtractor;
}

/// for particle simulation without using any algorithm related stuff
namespace PARTICLE
{
  class ParticleHandler
  {
  public:

    // constructor
    ParticleHandler(const int myrank);

    // constructor
    ParticleHandler(const Epetra_Comm& comm);

    // destructor
    virtual ~ParticleHandler(){};

    /// get binning strategy
    virtual inline Teuchos::RCP<BINSTRATEGY::BinningStrategy>& BinStrategy() {return binstrategy_;}

    virtual inline BINSTRATEGY::BinningStrategy const& BinStrategy() const {return *binstrategy_;}

    /// initial distribution of particles ( more general nodes of bindis_ ) to bins
    virtual void DistributeParticlesToBins(Teuchos::RCP<Epetra_Map> const& particlerowmap);

    /// get bin colume map
    virtual inline Teuchos::RCP<Epetra_Map>& BinColMap(){return bincolmap_;}

    /// get myrank
    virtual inline int MyRank() { return myrank_;}

    /// particles are checked whether they have moved out of their current bin
    /// and transferred if necessary
    virtual Teuchos::RCP<std::list<int> > TransferParticles(Teuchos::RCP<Epetra_Vector> disnp = Teuchos::null);

    /// node is placed into the correct row bin or put into the list of homeless particles
    virtual bool PlaceNodeCorrectly(
      Teuchos::RCP<DRT::Node> node,  ///< node to be placed
      const double* currpos,  ///< current position of this node
      std::list<Teuchos::RCP<DRT::Node> >& homelessparticles  ///< list of homeless particles
      );

    /// round robin loop to fill particles into its correct bin on according proc
    virtual void FillParticlesIntoBinsRoundRobin(
      std::list<Teuchos::RCP<DRT::Node> >& homelessparticles  ///< list of homeless particles
      );

    /// bins are distributed to the processors
    virtual Teuchos::RCP<Epetra_Map> DistributeBinsToProcs();

    /// setup ghosting of bins and particles
    virtual void SetupGhosting(
      Teuchos::RCP<Epetra_Map> binrowmap  ///< rowmap of bins
      );

    /// get neighbouring bins of particle containing boundary row bins
    virtual void GetNeighbouringBinsOfParticleContainingBoundaryRowBins(
        std::set<int>& colbins
      ) const;

  protected:

    /// create graph for bin distribution
    Teuchos::RCP<const Epetra_CrsGraph> CreateGraph();

    /// fill particles into their correct bin on according proc
    virtual Teuchos::RCP<std::list<int> > FillParticlesIntoBinsRemoteIdList(
      std::list<Teuchos::RCP<DRT::Node> >& homelessparticles   ///< set of homeless particles
      );

    /// bins are distributed to the processors
    virtual Teuchos::RCP<Epetra_Map> DistributeBinsToProcsBasedOnUnderlyingDiscret(
      Teuchos::RCP<DRT::Discretization> underlyingdis,  ///< underlying discretization on which distribution depends on
      std::map<int, std::set<int> >& fluideles  ///< map from bin id to row fluid elements in this bin
      );

  private:

    /// binning strategy
    Teuchos::RCP<BINSTRATEGY::BinningStrategy> binstrategy_;

    /// myrank
    const int myrank_;

    /// colmap of bins
    Teuchos::RCP<Epetra_Map> bincolmap_;

  };

}


/*----------------------------------------------------------------------*/
#endif
