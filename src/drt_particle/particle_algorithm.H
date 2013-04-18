/*----------------------------------------------------------------------*/
/*!
\file particle_algorithm.H

\brief Algorithm to control particle simulations

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-152537
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


  protected:

    /// prepare time step
    virtual void PrepareTimeStep();

    /// solve the current particle time step
    virtual void Integrate();

    /// force on particles is calculated and applied
    virtual void CalculateAndApplyForcesToParticles();

    /// update the current time step
    virtual void Update();

    /// output particle time step
    virtual void Output();

    /// find XAABB and divide into bins using the cutoff radius
    virtual void CreateBins();

    /// bins are distributed to the processors
    virtual Teuchos::RCP<Epetra_Map> DistributeBinsToProcs();

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
      Teuchos::RCP<DRT::Discretization> structdis
      );

    /// convert i,j,k into bin id
    virtual int ConvertijkToGid(
      int* ijk  ///< i,j,k to be converted into a bin id
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

    ///convert bin id into i,j,k
    virtual void ConvertGidToijk(
      int gid,  ///< bin id to be converted into i,j,k
      int* ijk  ///< resulting i,j,k
      );

    /// get all bin ids for given range of ijk
    void GidsInijkRange(
      int* ijk_bound,   ///< given range of ijk
      std::set<int>& binIds   ///< all bin ids in the specified range
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

    Teuchos::RCP<DRT::Discretization> particlewalldis_;     ///< fully redundant particle wall discretization

    //! Info: Currently it is assumed that bins and wall elements are fixed in space (initial content is stored here)
    //! Future: Dynamical registering of wall elements in bins when they are moving
    std::map<int, std::set<int> > particlewalls_;   ///< map that contains all wall element ids of each bin

    //@}

    //! @name output discretization writer
    Teuchos::RCP<IO::DiscretizationWriter> wall_output_;
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
