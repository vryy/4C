/*----------------------------------------------------------------------*/
/*!
\file scatra_particle_coupling.H

\brief Algorithm to control cavitation simulations

<pre>
Maintainer: Georg Hammerl / Ursula Rasthofer
            hammerl/rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237/15236
</pre>
*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                              ghamm 11/12 |
 *----------------------------------------------------------------------*/
#ifndef SCATRA_PARTICLE_COUPLING_H
#define SCATRA_PARTICLE_COUPLING_H

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 11/12 |
 *----------------------------------------------------------------------*/
#include "particle_algorithm.H"
#include "../drt_inpar/inpar_levelset.H"
#include "particle_timint_centrdiff.H"
#include "../drt_structure/strtimint_centrdiff.H"

#include "Epetra_MpiComm.h"
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>

/*----------------------------------------------------------------------*
 | forward declarations                                     ghamm 11/12 |
 *----------------------------------------------------------------------*/
namespace SCATRA
{
  class ScaTraTimIntImpl;
}
namespace DRT
{
  class Discretization;
  class Element;
}

/*----------------------------------------------------------------------*
 | scatra-particle algorithm                                ghamm 11/12 |
 *----------------------------------------------------------------------*/
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
    void Init(bool restarted);

    /// initialize particle field with particles
    void InitialSeeding();

    /// prepare time step
    void PrepareTimeStep();

    /// solve the current particle time step
    void Integrate();

    /// transfer particles into their correct bins and update solution
    void Update();

    /// correct level set function with escaped particles
    Teuchos::RCP<Epetra_Vector> CorrectionStep();

    /// output particle time step
    void Output();

  protected:

    /// setup ghosting of bins, particles & underlying scatra
    virtual void SetupGhosting(
      Teuchos::RCP<Epetra_Map> binrowmap,  ///< rowmap of bins
      std::map<int, std::set<int> >& fluideles   ///< map from bin id to scatra elements in this bin
      );

    /// build connectivity from fluid elements to bins
    virtual void BuildElementToBinPointers();

    /// set velocity of particle
    void SetVelocity();

     //@}

  private:

    /// compute level-set value of particle with local id
    double GetPhiParticle(
      const int particle_lid              ///< local particle id
      );

    /// get underlying element as well as particle position
    DRT::Element* GetEleCoordinates(
      const int particle_lid,             ///< local particle id
      LINALG::Matrix<3,1>& elecoord       ///< matrix to be filled with particle coordinates in element space
    );

    /// compute level-set value of particle with local id
    void GetLSValParticle(
      double& phi_particle,               ///< phi value at particle position
      LINALG::Matrix<3,1>& gradphi,       ///< gradient of level-set field at particle position
      DRT::Element* scatraele,            ///< corresponding scatra element
      const LINALG::Matrix<3,1> elecoord  ///< matrix containing particle coordinates in element space
    );

    //! @name scatra related stuff
    //@{

    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_;     ///< underlying scatra time integrator

    Teuchos::RCP<DRT::Discretization> scatradis_;       ///< underlying scatra discretization

    const Teuchos::RCP<Teuchos::ParameterList> params_;       ///< the parameter list

    INPAR::PARTICLE::ParticleDim particle_dim_;         ///< dimension of problem

    //@}

  }; // Algorithm


} // namespace PARTICLE


/*----------------------------------------------------------------------*/
#endif  // SCATRA_PARTICLE_COUPLING_H
