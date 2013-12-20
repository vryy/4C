/*----------------------------------------------------------------------*/
/*!
\file particle_algorithm.H

\brief Algorithm to control particle simulations

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
#ifndef PARTICLE_ALGORITHM_H
#define PARTICLE_ALGORITHM_H

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
namespace ADAPTER
{
  class Structure;
}
namespace IO
{
  class DiscretizationWriter;
}
namespace LINALG
{
  class MapExtractor;
}

/*----------------------------------------------------------------------*
 | particle algorithm                                       ghamm 09/12 |
 *----------------------------------------------------------------------*/
/// PARTICLE: Particle simulation
namespace PARTICLE
{
  // forward declaration
  class Less;

  class Algorithm : public ADAPTER::AlgorithmBase
  {
  public:

    Algorithm(
      const Epetra_Comm& lcomm,  ///< local epetra communicator
      const Teuchos::ParameterList& params   ///< problem parameters
      );

    virtual ~Algorithm(){};

    /// do initialization of problem, reduced version for restart
    virtual void Init(bool restarted);

    /// read restart
    virtual void ReadRestart(
      int restart  ///< time step from which to restart
      );

    /// setup of the system
    virtual void SetupSystem();

    /// outer time loop for particle problem
    virtual void Timeloop();

    /// test results (if necessary)
    virtual void TestResults(
      const Epetra_Comm& comm  ///< local epetra communicator
      );

    /// get wall discretization 
    Teuchos::RCP<DRT::Discretization> WallDiscret() { return particlewalldis_; }

    /// get underlyling structural problem
    Teuchos::RCP<ADAPTER::Structure> Structure() { return structure_; }

    /// convert i,j,k into bin id
    virtual int ConvertijkToGid(
      int* ijk  ///< i,j,k to be converted into a bin id
      );

    /// get all bin ids for given range of ijk
    void GidsInijkRange(
      int* ijk_bound,   ///< given range of ijk
      std::set<int>& binIds,   ///< all bin ids in the specified range
      bool checkexistence   ///< check can be added whether the gids are on this proc
      );

    ///convert bin id into i,j,k
    virtual void ConvertGidToijk(
      int gid,  ///< bin id to be converted into i,j,k
      int* ijk  ///< resulting i,j,k
      );

    // return gravity acceleration
    LINALG::Matrix<3,1> GetGravityAcc() {return gravity_acc_ ; }


  protected:

    /// prepare time step
    virtual void PrepareTimeStep();

    /// solve the current particle time step
    virtual void Integrate();

    /// force on particles is calculated and applied
    virtual void CalculateAndApplyForcesToParticles();

    /// update the current time step
    virtual void Update();

    /// calculate stresses, strains, energies
    virtual void PrepareOutput();

    /// output particle time step
    virtual void Output();

    /// find XAABB and divide into bins using the cutoff radius
    virtual void CreateBins(
      Teuchos::RCP<DRT::Discretization> dis   ///< based on this discretization the bins are created
      );

    /// bins are distributed to the processors
    virtual Teuchos::RCP<Epetra_Map> DistributeBinsToProcs();

    // bins are redistributed according to the current fill level with particles
    virtual void DynamicLoadBalancing();

    // create graph for bin distribution
    Teuchos::RCP<const Epetra_CrsGraph> CreateGraph();

    /// round robin loop to fill particles into its correct bin on according proc
    virtual void FillParticlesIntoBins(
      std::set<Teuchos::RCP<DRT::Node>,Less>& homelessparticles  ///< set of homeless particles
      );

    /// node is placed into the correct row bin or put into the set of homeless particles
    virtual bool PlaceNodeCorrectly(
      Teuchos::RCP<DRT::Node> node,  ///< node to be placed
      const double* currpos,  ///< current position of this node
      std::set<Teuchos::RCP<DRT::Node>,Less>& homelessparticles  ///< set of homeless particles
      );

    /// setup ghosting of bins and particles
    virtual void SetupGhosting(
      Teuchos::RCP<Epetra_Map> binrowmap  ///< rowmap of bins
      );

    /// particles are checked whether they have moved out of their current bin
    /// and transferred if necessary
    virtual void TransferParticles();

    /// particle walls are set up and assigned to bins
    virtual void SetupParticleWalls(
      Teuchos::RCP<DRT::Discretization> basediscret
      );

    /// build connectivity from particle wall elements to bins
    virtual void BuildElementToBinPointers(
      bool wallpointer   ///< bool whether wall element pointer need rebuild
      );

    /// bins are distributed to the processors
    virtual Teuchos::RCP<Epetra_Map> DistributeBinsToProcsBasedOnUnderlyingDiscret(
      Teuchos::RCP<DRT::Discretization> underlyingdis,  ///< underlying discretization on which distribution depends on
      std::map<int, std::set<int> >& fluideles  ///< map from bin id to fluid elements in this bin
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

    /// assign wall elements to bins
    void AssignWallElesToBins();

    //! @name particle related stuff
    //@{

    Teuchos::RCP<DRT::Discretization> particledis_;     ///< particle discretization

    Teuchos::RCP<ADAPTER::Structure> particles_;     ///< particle time integration

    Teuchos::RCP<Epetra_Map> bincolmap_;     ///< colmap of bins

    double cutoff_radius_;    ///< maximum influence radius of particles

    LINALG::Matrix<3,2> XAABB_;    ///< extendedAxisAlignedBoundingBox of all particles

    double bin_size_[3];    ///< size of each bin in Cartesian coordinates

    int bin_per_dir_[3];    ///< number of bins per direction

    LINALG::Matrix<3,1> gravity_acc_;   ///< acceleration due to gravity

    Teuchos::RCP<ADAPTER::Structure> structure_;     ///< structural (wall) time integration

    Teuchos::RCP<DRT::Discretization> particlewalldis_;     ///< fully redundant particle wall discretization

    Teuchos::RCP<LINALG::MapExtractor> wallextractor_;   ///< map extractor for particle wall

    bool moving_walls_;

    //@}

    //! @name miscellaneous stuff
    //@{

    int myrank_;    ///< myrank

    //@}

  }; // Algorithm

  // class for comparing Teuchos::RCP<DRT::Node> in a std::set
  class Less
  {
  public:
    bool operator()(const Teuchos::RCP<const DRT::Node>& first, const Teuchos::RCP<const DRT::Node>& second) const;
  };

} // namespace PARTICLE


/*----------------------------------------------------------------------*/
#endif  // PARTICLE_ALGORITHM_H
