/*----------------------------------------------------------------------*/
/*!
\file scatra_particle_coupling.H

\brief Algorithm to control cavitation simulations

\level 3

<pre>
\maintainer Magnus Winter
            winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*----------------------------------------------------------------------*/

#ifndef SCATRA_PARTICLE_COUPLING_H
#define SCATRA_PARTICLE_COUPLING_H

#include "particle_algorithm.H"
#include "../drt_inpar/inpar_levelset.H"
#include "particle_timint_centrdiff.H"
#include "../drt_structure/strtimint_centrdiff.H"

#include "Epetra_MpiComm.h"
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>

namespace SCATRA
{
  class ScaTraTimIntImpl;
}
namespace DRT
{
  class Discretization;
  class Element;

  namespace MESHFREE
  {
    class MeshfreeMultiBin;
  }
}


namespace PARTICLE
{
  class ScatraParticleCoupling : public PARTICLE::Algorithm
  {
  public:

    ScatraParticleCoupling(
      Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,   ///< underlying scatra field
      Teuchos::RCP<Teuchos::ParameterList> params      ///< problem parameters
      );

    /// do initialization of problem
    virtual void Init(bool restarted);

    /// initialize particle field with particles
    virtual void InitialSeeding();

    /// solve the current particle time step
    virtual void Integrate();

    /// transfer particles into their correct bins and update solution
    virtual void Update();

    /// correct level-set function using escaped particles
    Teuchos::RCP<Epetra_Vector> CorrectionStep();

    /// compute velocities of particles at intermediate time n+theta
    Teuchos::RCP<Epetra_Vector> GetVelocity(const double theta); /// intermediate point in time

    /// set time and step from outside (i.e., particle algorithm is not main algorithm)
    virtual void SetTimeStepAdditionalParticles(const double time, const int step)
    {
      SetTimeStep(time,step);
      particles_->SetTimeStep(time,step);
      return;
    }

  protected:

    /// setup ghosting of bins, particles & underlying scatra
    virtual void SetupGhosting(
      Teuchos::RCP<Epetra_Map> binrowmap,  ///< rowmap of bins
      std::map<int, std::set<int> >& fluideles   ///< map from bin id to scatra elements in this bin
    );

    /// build connectivity from fluid elements to bins
    virtual void BuildElementToBinPointers();

  private:

    /// adjust particle radii to current interface position
    void AdjustParticleRadii();

    /// reseeding of particles
    void Reseeding();

    /// attract particle to their correct side of interface
    void Attraction(
      const int global_num_particles,
      Teuchos::RCP<Epetra_Vector> inc_dis,
      Teuchos::RCP<Epetra_Vector> lambda,
      Teuchos::RCP<Epetra_Vector> phi_target,
      Teuchos::RCP<Epetra_Vector> phi_particle,
      Teuchos::RCP<Epetra_Vector> norm_gradphi_particle,
      Teuchos::RCP<Epetra_Vector> disn
    );

    /// delete particles
    void DeleteParticles(std::vector<int> part_del); ///< global ids of particles to be deleted

    /// get underlying element as well as particle position
    DRT::Element* GetEleCoordinates(
      const int particle_gid,             ///< global particle id
      LINALG::Matrix<3,1>& elecoord       ///< matrix to be filled with particle coordinates in element space
    );

    /// get underlying element as well as position in element space
    DRT::Element* GetEleCoordinatesFromPosition(
      LINALG::Matrix<3,1>& myposition,          ///< position
      DRT::MESHFREE::MeshfreeMultiBin* currbin, ///< corresponding bin
      LINALG::Matrix<3,1>& elecoord             ///< matrix to be filled with particle coordinates in element space
    );

    /// compute level-set value of particle
    void GetLSValParticle(
      double& phi_particle,               ///< phi value at particle position
      LINALG::Matrix<3,1>& gradphi,       ///< gradient of level-set field at particle position
      DRT::Element* scatraele,            ///< corresponding scatra element
      const LINALG::Matrix<3,1> elecoord, ///< matrix containing particle coordinates in element space
      const Teuchos::RCP<const Epetra_Vector> phinp, ///< vector containing level-set values (col map)
      bool compute_normal = true          ///< additionally compute normal vector
    );


    /// compute level-set value of particle
    void GetVelParticle(
      LINALG::Matrix<3,1>& vel,           ///< velocity of level-set field at particle position
      DRT::Element* scatraele,            ///< corresponding scatra element
      const LINALG::Matrix<3,1> elecoord, ///< matrix containing particle coordinates in element space
      const Teuchos::RCP<Epetra_Vector> lsvel ///< vector containing level-set values (colmap)
    );

    /// delete particles outside of the physical domain
    void DeleteParticlesOutOfPhysicalDomain();

    /// compares the second entry of a pair<int,double>. To be passed to the sorting algo (helper function for Reseeding())
    /// descending order
    static bool MyComparePairs (const std::pair<int,double>& first, const std::pair<int,double>& second)
    {
      if (first.second > second.second)
        return true;
      else
        return false;
    };

    BINSTRATEGY::UTILS::BinContentType bin_transpcontent_;

    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_;     ///< underlying scatra time integrator

    Teuchos::RCP<DRT::Discretization> scatradis_;       ///< underlying scatra discretization

    const Teuchos::RCP<Teuchos::ParameterList> params_;       ///< the parameter list

    INPAR::PARTICLE::Escaped escaped_;                  ///< escaped criterion

    int reseeding_;                                     ///< reseeding frequency

    int num_particles_per_bin_;                         ///< number of particles per bin of each type

    double binlength_max_;                              ///< characteristic bin length: maximal edge length

    double binlength_min_;                              ///< characteristic bin length: minimal edge length

    double b_max_;                                      ///< width of band with particles

    double r_min_;                                      ///< minimal radius of particles

    double r_max_;                                      ///< maximal radius of particles

    bool fast_;                                         ///< use fast check based on nodal values of scatra elements to determine far-away bins

    double delete_more_;                                ///< allows for deleting also escaped particles that moved to far from the interface to be useful

  }; // ScatraParticleCoupling


} // namespace PARTICLE


/*----------------------------------------------------------------------*/
#endif  // SCATRA_PARTICLE_COUPLING_H
