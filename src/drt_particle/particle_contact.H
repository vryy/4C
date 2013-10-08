/*----------------------------------------------------------------------*/
/*!
\file particle_contact.H
\brief Particle collision handling

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef PARTICLE_CONTACT_H
#define PARTICLE_CONTACT_H

/*----------------------------------------------------------------------*/
/* headers */
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_inpar/inpar_cavitation.H"
#include "particle_node.H"

/*----------------------------------------------------------------------*/
/* defines */
#define energyoutput

// forward declarations
namespace DRT
{
	class Element;
	class Node;
}
namespace PARTICLE
{
  class Algorithm;
}
/*----------------------------------------------------------------------*/
/* belongs to particle namespace */
namespace PARTICLE
{

  /*====================================================================*/
  /*!
   * \brief Collision handling with the discrete element method
   *  assuming rigid spherical particles
   *  For details see Bachelor thesis Niklas Fehn
   *
   * \author hammerl
   * \date 09/13
   */
  class ParticleCollisionHandler
  {

  public:

    //! @name Life
    //@{

    //! Constructor
    ParticleCollisionHandler(Teuchos::RCP<DRT::Discretization> discret, Teuchos::RCP<PARTICLE::Algorithm> particlealgorithm);

    //@}

    //! @name Actions
    //@{

    /// set states from time integrator to prepare collision calculation
    void SetState(
      Teuchos::RCP<Epetra_Vector> radius,
      Teuchos::RCP<Epetra_Vector> mass
      );

    /// compute collisions (inter-particle and particle-wall), returns contact energy
    double ComputeCollisions(double dt, Teuchos::RCP<Epetra_Vector>& f_contact, Teuchos::RCP<Epetra_Vector>& m_contact);

    /// get neighboring particles and wall elements
    void FindNeighbors(
      std::set<DRT::Node*>& NeighboringParticles,
      std::set<DRT::Element*>& NeighboringWallEles,
      std::set<int>& binIds
      );

    /// calculate normal contact force for single contact pair
    void CalculateNormalContactForce(
        double g,
        double v_rel_normal,
        double m_eff,
        double& ContactForce,
        int owner_i,
        int owner_j
        );

    /// calculate tangential contact force
    void CalculateTangentialContactForce(
        double ContactForce,
        double *normal,
        double *TangentialContactForce,
        PARTICLE::Collision &currentColl,
        double *v_rel_tangential,
        double m_eff,
        const double dt,
        int owner_i,
        int owner_j
        );

    /// read initial contact parameters and validate them
    void ReadContactParameters(double density);

    /// assemble energies of particles
    double EnergyAssemble(double owner_i, double owner_j);

    /// print maximum penetration of particles
    void PrintMaxPenetration(int step, double time);

    //@}

    //! @name Access methods
    //@{

    /// return minimal radius to be expected
    double GetMinRadius() { return r_min_; }

    // return maximum radius to be expected
    double GetMaxRadius() { return r_max_; }

    //@}

  private:

    /// ID of actual processor in parallel
    int myrank_;

    //! @name helper classes
    //@{

    /// particle discretization
    Teuchos::RCP<DRT::Discretization> discret_;

    /// particle algorithm that is connected to the time integration
    Teuchos::RCP<PARTICLE::Algorithm> particle_algorithm_;

    //@}

    //! @name Collision parameters
    //@{

    /// overall contact strategy
    INPAR::PARTICLE::ContactStrategy contact_strategy_;

    /// normal contact model
    INPAR::PARTICLE::NormalContact normal_contact_;

    /// particle Poisson ratio
    double nue_;
    /// particle Young's modulus
    double young_;

    /// minimum expected radius
    double r_min_;
    /// maximum expected radius
    double r_max_;
    /// max expected particle velocity
    double v_max_;

    /// max allowed relative penetration
    double c_;
    /// critical time step for contact
    double dt_krit_;

    /// coefficient of restitution (inter-particle)
    double e_;
    /// coefficient of restitution (particle-wall)
    double e_wall_;

    /// normal contact stiffness (inter-particle)
    double k_normal_;
    /// tangential contact stiffness (inter-particle)
    double k_tang_;
    /// tangential contact stiffness (particle-wall)
    double k_tang_wall_;

    /// tangential to normal stiffness ratio (inter-particle)
    double kappa_;
    /// tangential to normal stiffness ratio (particle-wall)
    double kappa_wall_;

    /// damping parameter normal contact (inter-particle)
    double d_normal_;
    /// damping parameter tangential contact (inter-particle)
    double d_tang_;
    /// damping parameter normal contact (particle-wall)
    double d_normal_wall_;
    /// damping parameter tangential contact (particle-wall)
    double d_tang_wall_;

    /// frictional coefficient (inter-particle)
    double mu_;
    /// frictional coefficient (particle-wall)
    double mu_wall_;

    // tension cutoff
    bool tension_cutoff_;

    /// energy stored in elastic springs due to contact
    double contact_energy_;
    /// maximum penetration occured
    double g_max_;

    //@}

    //! @name Global properties in col layout at \f$t_{n+1}\f$ for collision detection
    //@{

    // node-based
    Teuchos::RCP<Epetra_Vector> radiusncol_;
    Teuchos::RCP<Epetra_Vector> masscol_;

    // dof-based
    Teuchos::RCP<const Epetra_Vector> disncol_;
    Teuchos::RCP<const Epetra_Vector> velncol_;
    Teuchos::RCP<const Epetra_Vector> ang_velncol_;

    //@}

    /// struct to store wall contact point and additional information about corresponding wall ele
    struct WallContactPoint {
      LINALG::Matrix<3,1> point;
      double penetration;
      std::map<int,LINALG::Matrix<3,1> > nodalCoordinates;
      std::vector<int> lm;
      std::vector<int> lmowner;
    };

  };  // class ParticleCollisionHandler

}  // namespace PARTICLE

/*----------------------------------------------------------------------*/
#endif  // #ifndef PARTICLE_CONTACT_H
