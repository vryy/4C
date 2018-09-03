/*----------------------------------------------------------------------*/
/*!
\file particle_contact.H

\brief Particle collision handling

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*-----------------------------------------------------------------------*/
#ifndef PARTICLE_CONTACT_H
#define PARTICLE_CONTACT_H

/*----------------------------------------------------------------------*/
/* headers */
#include "../drt_geometry/searchtree_nearestobject.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_inpar/inpar_particle.H"
#include "particle_node.H"
#include "../drt_lib/drt_element.H"
#include <Epetra_FEVector.h>
#include <boost/unordered_map.hpp>


// forward declarations
namespace DRT
{
  class Element;
  class Node;
}
namespace PARTICLE
{
  class Algorithm;
  struct WallContactPoint;
  struct ParticleCollData;
}
namespace BINSTRATEGY
{
  class Less;
}

/*----------------------------------------------------------------------*/
/* belongs to particle namespace */
namespace PARTICLE
{
  /*====================================================================*/
  /*!
   * \brief Collision handler base for particle contact
   *
   * \author hammerl
   * \date 09/13
   */
  class ParticleCollisionHandlerBase
  {

  public:

    //! @name Life
    //@{

    //! Constructor
    ParticleCollisionHandlerBase(
      Teuchos::RCP<DRT::Discretization> discret,
      Teuchos::RCP<PARTICLE::Algorithm> particlealgorithm,
      const Teuchos::ParameterList& particledynparams
      );

    //! Destructor
    virtual ~ParticleCollisionHandlerBase(){ ; }

    //@}

    //! @name Actions
    //@{

    /// set states from time integrator to prepare collision calculation
    virtual void Init(
        Teuchos::RCP<Epetra_Vector> disn,
        Teuchos::RCP<Epetra_Vector> veln,
        Teuchos::RCP<Epetra_Vector> angVeln,
        Teuchos::RCP<const Epetra_Vector> radiusn,
        Teuchos::RCP<Epetra_Vector> anglen,
        Teuchos::RCP<Epetra_Vector> mass
        );

    /// compute collisions (inter-particle and particle-wall), returns contact energy
    virtual double EvaluateParticleContact(
      double dt,
      Teuchos::RCP<Epetra_Vector> sysvec1,
      Teuchos::RCP<Epetra_Vector> sysvec2,
      Teuchos::RCP<Epetra_FEVector> f_structure = Teuchos::null,
      int timestep = 0
      ) = 0;

    /// assemble energies of particles
    virtual double EnergyAssemble(double owner_i, double owner_j);

    //@}

    //! @name Access methods
    //@{

    /// return minimal radius to be expected
    virtual double GetMinRadius() { return r_min_; }

    // return maximum radius to be expected
    virtual double GetMaxRadius() { return r_max_; }

    /// return maximum penetration of particles
    virtual double GetMaxPenetration() {return g_max_; }

    /// return maximum penetration of particles
    virtual double GetMaxPenetrationParticle() {return g_max_particle_; }

    /// return maximum penetration of particles
    virtual double GetMaxPenetrationWall() {return g_max_wall_; }

    /// return maximum penetration of particles
    virtual double GetMaxPenetrationParticleRel() {return g_max_particle_rel_; }

    /// return maximum penetration of particles
    virtual double GetMaxPenetrationWallRel() {return g_max_wall_rel_; }

    /// return maximum occurring relative stress
    virtual double GetMaxYieldStressRel() {return yield_max_rel_; }

    /// return maximum penetration of particles
    virtual double GetAllTimeMaxPenetrationParticle() {return g_alltimemax_particle_; }

    /// return maximum penetration of particles
    virtual double GetAllTimeMaxPenetrationWall() {return g_alltimemax_wall_; }

    /// return maximum penetration of particles
    virtual double GetAllTimeMaxPenetrationParticleRel() {return g_alltimemax_particle_rel_; }

    /// return maximum penetration of particles
    virtual double GetAllTimeMaxPenetrationWallRel() {return g_alltimemax_wall_rel_; }

    /// return maximum occurring relative stress
    virtual double GetAllTimeMaxYieldStressRel() {return yield_alltimemax_rel_; }


    //@}

  protected:

    //! @name Collision parameters
    //@{

    /// normal contact model
    const INPAR::PARTICLE::NormalContact normal_contact_;

    /// rolling contact model
    const INPAR::PARTICLE::RollingContact rolling_contact_;

    /// adhesion law governing normal contact of particles
    const INPAR::PARTICLE::NormalAdhesion normal_adhesion_;

    /// particle Poisson ratio
    double nue_;
    /// particle Young's modulus
    double young_;
    /// particle yield strength
    double yield_;

    /// minimum expected radius
    double r_min_;
    /// maximum expected radius
    double r_max_;
    /// max expected particle velocity
    const double v_max_;
    /// max allowed relative particle-particle penetration
    const double c_;
    /// max allowed relative particle-wall penetration
    const double c_wall_;
    /// critical time step for contact
    double dt_krit_;

    /// coefficient of restitution (inter-particle)
    const double e_;
    /// coefficient of restitution (particle-wall)
    const double e_wall_;

    /// normal contact stiffness (inter-particle)
    double k_normal_;
    /// normal contact stiffness (particle-wall)
    double k_normal_wall_;
    /// tangential contact stiffness (inter-particle)
    double k_tang_;
    /// tangential contact stiffness (particle-wall)
    double k_tang_wall_;
    /// rolling stiffness (particle-wall)
    double k_roll_wall_;
    /// rolling stiffness (inter-particle)
    double k_roll_;
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
    /// damping parameter rolling contact (particle-wall)
    double d_roll_wall_;
    /// damping parameter rolling contact (inter-particle)
    double d_roll_;

    /// frictional coefficient (inter-particle)
    const double mu_;
    /// frictional coefficient (particle-wall)
    const double mu_wall_;
    /// rolling frictional coefficient (particle-wall)
    const double mu_roll_wall_;
    /// rolling frictional coefficient (inter-particle)
    const double mu_roll_;

    // tension cutoff
    const bool tension_cutoff_;

    // parameter for regularization of damping contact force
    double damp_reg_fac_;

    /// energy stored in elastic springs due to contact
    double contact_energy_;

    /// maximum penetration occurred
    double g_max_;

    /// maximum penetration between particles
    double g_max_particle_;
    /// maximum penetration between particle and wall
    double g_max_wall_;
    /// maximum penetration between particles
    double g_max_particle_rel_;
    /// maximum penetration between particle and wall
    double g_max_wall_rel_;
    /// maximum percentage of yield strength occurred
    double yield_max_rel_;
    /// maximum penetration between particles
    double g_alltimemax_particle_;
    /// maximum penetration between particle and wall
    double g_alltimemax_wall_;
    /// maximum penetration between particles
    double g_alltimemax_particle_rel_;
    /// maximum penetration between particle and wall
    double g_alltimemax_wall_rel_;
    /// maximum percentage of yield strength occurred
    double yield_alltimemax_rel_;

    /// Wallcontact counter
    LINALG::Matrix<30,1> wallcontact_counter_;

    /// gap between two particles or between particle and wall in adhesion equilibrium
    const double adhesion_eq_gap_;
    /// stiffness for normal adhesion force
    const double adhesion_normal_stiff_;
    /// damping coefficient for normal adhesion force
    const double adhesion_normal_damp_;
    /// depth of Lennard-Jones adhesion potential well
    const double adhesion_normal_eps_;
    /// maximum adhesion force
    const double adhesion_max_force_;
    /// maximum displacement from adhesion equilibrium
    const double adhesion_max_disp_;
    /// surface energy density for the calculation of the pull-out force
    const double adhesion_surface_energy_;
    const double adhesion_surface_energy_factor_;
    double adhesion_surface_energy_wall_;
    double adhesion_surface_energy_wall_factor_;
    /// surface energy density for the calculation of the pull-out force
    const int adhesion_maxwalleleID_;
    const int adhesion_minwalleleID_;
    /// surface energy distribution
    const INPAR::PARTICLE::AdhesionSurfaceEnergyDistribution adhesion_surface_energy_distribution_;
    /// standard deviation of adhesion surface energy
    const double adhesion_surface_energy_distribution_sigma_;
    /// surface energy distribution limited by multiple of standard deviation
    const double adhesion_surface_energy_cutoff_factor_;
    /// contact force where maximum pullout force is achieved
    const double adhesion_max_contact_pressure_;
    double adhesion_max_contact_pressure_wall_;
    bool adhesion_variant_max_contact_force_;
    /// shift van-der-Waals-curve to g = 0
    int adhesion_vdW_curve_shift_;

    //! write system energy every given step
    const int writeenergyevery_;

    //@}

    /// ID of actual processor in parallel
    int myrank_;

    //! @name helper classes
    //@{

    /// particle discretization
    const Teuchos::RCP<const DRT::Discretization> discret_;

    /// particle algorithm that is connected to the time integration
    const Teuchos::RCP<const PARTICLE::Algorithm> particle_algorithm_;

    //@}

    //! @name Global properties in col layout at \f$t_{n+1}\f$ for collision detection
    //@{
    std::vector<ParticleCollData> particleData_;   ///< temporary storage for all particle collision relevant data
  };  // class ParticleCollisionHandlerBase

  /*====================================================================*/
  /*!
   * \brief Collision handling with the discrete element method
   *  assuming rigid spherical particles
   *  For details see Bachelor thesis Niklas Fehn
   *
   * \author hammerl
   * \date 09/13
   */
  class ParticleCollisionHandlerDEM : public PARTICLE::ParticleCollisionHandlerBase
  {

  public:

  //! @name Life
  //@{

  //! Constructor
  ParticleCollisionHandlerDEM(
    Teuchos::RCP<DRT::Discretization> discret,
    Teuchos::RCP<PARTICLE::Algorithm> particlealgorithm,
    const Teuchos::ParameterList& particledynparams
    );

  //! Destructor
  virtual ~ParticleCollisionHandlerDEM(){ ; }

  //@}

  //! @name Actions
  //@{

  /// compute collisions (inter-particle and particle-wall), returns contact energy
  virtual double EvaluateParticleContact(
    const double dt,
    Teuchos::RCP<Epetra_Vector> f_contact,
    Teuchos::RCP<Epetra_Vector> m_contact,
    Teuchos::RCP<Epetra_FEVector> f_structure,
    int timestep);

  /// calculate contact with neighboring particles
  void CalcNeighboringParticlesContact(
      DRT::Node*                           particle_i,              //!< i-th particle
      ParticleCollData&                    data_i,                  //!< collision data associated with i-th particle
      const std::list<DRT::Node*>&         neighboring_particles,   //!< list of particles adjacent to i-th particle
      const bool                           havepbc,                 //!< flag indicating periodic boundaries
      const double                         dt,                      //!< time step size
      const Teuchos::RCP<Epetra_Vector>&   f_contact,               //!< global force vector
      const Teuchos::RCP<Epetra_Vector>&   m_contact,               //!< global moment vector
      int timestep
      );

  /// calculate contact with neighboring walls
  void CalcNeighboringWallsContact(
      DRT::Node*                                        particle_i,          //!< i-th particle
      const ParticleCollData&                           data_i,              //!< collision data associated with i-th particle
      const boost::unordered_map<int,DRT::Element*>&    neighboring_walls,   //!< map of wall elements adjacent to i-th particle
      const double                                      dt,                  //!< time step size
      const Teuchos::RCP<DRT::Discretization>&          walldiscret,         //!< wall discretization
      const Teuchos::RCP<const Epetra_Vector>&          walldisn,            //!< wall displacement
      const Teuchos::RCP<const Epetra_Vector>&          wallveln,            //!< wall velocity
      const Teuchos::RCP<Epetra_Vector>&                f_contact,           //!< global force vector
      const Teuchos::RCP<Epetra_Vector>&                m_contact,           //!< global moment vector
      const Teuchos::RCP<Epetra_FEVector>&              f_structure,          //!< global wall force vector
      int timestep
      );

  /// calculate normal contact force for single contact pair
  void CalculateNormalContactForce(
    const double g,
    const double radius_i,
    const double radius_j,
    const double v_rel_normal,
    const double m_eff,
    PARTICLE::Adhesion &currentAdhe,
    double& totalnormalforce,
    double& normalcontactforce,
    const int owner_i,
    const int owner_j,
    bool iswallele=false
    );

  /// calculate tangential contact force
  virtual void CalculateTangentialContactForce(
    const double normalcontactforce,
    const LINALG::Matrix<3,1>& normal,
    LINALG::Matrix<3,1>& tangentcontactforce,
    PARTICLE::Collision &currentColl,
    const LINALG::Matrix<3,1>& v_rel_tangential,
    const double m_eff,
    const double dt,
    const int owner_i,
    const int owner_j
    );

  /// calculate tangential contact force
  virtual void CalculateRollingContactForce(
    const double normalcontactforce,
    const LINALG::Matrix<3,1>& normal,
    LINALG::Matrix<3,1>& tangentrollingforce,
    PARTICLE::Collision &currentColl,
    const LINALG::Matrix<3,1>& v_roll,
    const double m_eff,
    const double dt,
    const int owner_i,
    const int owner_j
    );

  /// calculate gap at which vdW-curve intersects linear ramp
  virtual void CalculateIntersectionGap(
      double a,
      double b,
      double c,
      double d,
      double& gap_intersect);

  /// calculate virtual quasi-force for rolling resistance
  virtual void CalculateVirtualRollingForce(
    const double normalcontactforce,
    const LINALG::Matrix<3,1>& normal,
    LINALG::Matrix<3,1>& tangentrollingforce,
  const LINALG::Matrix<3,1> angvel_i,
  const LINALG::Matrix<3,1> angvel_j,
    const double radius_i,
    const double radius_j,
    const int owner_i,
    const int owner_j,
    bool iswallele=false
    );

  /// check whether normal particle adhesion needs to be evaluated
  bool ConsiderNormalAdhesion(const double normalgap) const;
  //@}

  protected:

  /// compute contact normal, contact gap, and vectors from particle centers to contact point
  virtual void ComputeContactPointAndNormalAndGap(
      LINALG::Matrix<3,1>&      r_iC,     //!< vector from center of particle i to contact point C
      LINALG::Matrix<3,1>&      r_jC,     //!< vector from center of particle j to contact point C
      LINALG::Matrix<3,1>&      normal,   //!< contact normal
      double&                   gap,      //!< contact gap
      const ParticleCollData&   data_i,   //!< collision data for particle i
      const ParticleCollData&   data_j    //!< collision data for particle j
      ) const;

  /// compute contact normal, contact gap, contact point on wall element, and type of nearest contact object
  virtual void ComputeContactPointAndNormalAndGapAndObjectType(
      LINALG::Matrix<3,1>&                        C,              //!< contact point on wall element
      GEO::ObjectType&                            objecttype,     //!< type of nearest contact object
      LINALG::Matrix<3,1>&                        normal,         //!< contact normal
      double&                                     gap,            //!< contact gap
      DRT::Element*                               wallele,        //!< wall element
      const std::map<int,LINALG::Matrix<3,1> >&   nodecoord,      //!< node coordinates of wall element
      const ParticleCollData&                     data            //!< particle collision data
      ) const;

  /// compute relative velocity in contact point between two particles
  void ComputeRelativeVelocity(
      LINALG::Matrix<3,1>&         v_rel,    //!< relative velocity between two particles
      const LINALG::Matrix<3,1>&   r_iC,     //!< vector from center of particle i to contact point C
      const LINALG::Matrix<3,1>&   r_jC,     //!< vector from center of particle j to contact point C
      const ParticleCollData&      data_i,   //!< collision data for particle i
      const ParticleCollData&      data_j    //!< collision data for particle j
      ) const;

  /// compute relative velocity in contact point between particle and wall element
  void ComputeRelativeVelocity(
      LINALG::Matrix<3,1>&         v_rel,   //!< relative velocity between particle and wall element
      const LINALG::Matrix<3,1>&   r_C,     //!< vector from center of particle to contact point C on wall element
      const LINALG::Matrix<3,1>&   v_C,     //!< velocity of contact point C on wall element
      const ParticleCollData&      data     //!< particle collision data
      ) const;

  /// compute rolling velocity in contact point between two particles
  void ComputeRollingVelocity(
      LINALG::Matrix<3,1>&         v_roll,   //!< rolling velocity between two particles
      const LINALG::Matrix<3,1>&   normal,   //!< contact normal
      const LINALG::Matrix<3,1>&   r_iC,     //!< vector from center of particle i to contact point C
      const LINALG::Matrix<3,1>&   r_jC,     //!< vector from center of particle j to contact point C
      const ParticleCollData&      data_i,   //!< collision data for particle i
      const ParticleCollData&      data_j    //!< collision data for particle j
      ) const;

  /// compute rolling velocity in contact point between particle and wall element
  void ComputeRollingVelocity(
      LINALG::Matrix<3,1>&         v_roll,   //!< rolling velocity between two particles
      const LINALG::Matrix<3,1>&   r_C,     //!< vector from center of particle to contact point C
      const ParticleCollData&      data   //!< collision data for particle
      ) const;

  void ComputeAdhesionSurfaceEnergyDistribution(
      PARTICLE::Adhesion &currentAdhe,
      bool iswallele
  ) const;

  //! gather *.csv files with normal particle-particle and particle-wall contact forces across all processors
  void GatherNormalContactForcesToFile() const;

  //! write normal contact force between two particles or between particle and wall into *.csv file
  void OutputNormalContactForceToFile(
      const double&                normalcontactforce,     //!< normal contact force
      const LINALG::Matrix<3,1>&   forceapplicationpoint   //!< force application point
      ) const;

  //! write particle-particle and particle-wall contact forces into *.csv file every given step
  const int writecontactforcesevery_;
   };  // class ParticleCollisionHandlerDEM


  /*====================================================================*/
  /*!
   * \brief collision handling with the discrete element method, assuming rigid, ellipsoidal particles
   *        (for details see bachelor thesis by Christoph Haslinger)
   * \author fang
   * \date 09/17
   */
  class ParticleCollisionHandlerDEMEllipsoids : public PARTICLE::ParticleCollisionHandlerDEM
  {
    public:

      //! constructor
      ParticleCollisionHandlerDEMEllipsoids(
          Teuchos::RCP<DRT::Discretization>   discret,             //!< particle discretization
          Teuchos::RCP<PARTICLE::Algorithm>   particlealgorithm,   //!< particle algorithm
          const Teuchos::ParameterList&       particledynparams    //!< parameter list
          );

      //! destructor
      virtual ~ParticleCollisionHandlerDEMEllipsoids(){return;};

    private:

      //! check whether two ellipsoids collide
      bool CollisionCheck(
          const LINALG::Matrix<4,4>&   E_i,   //!< first ellipsoid matrix in homogeneous coordinates
          const LINALG::Matrix<4,4>&   E_j    //!< second ellipsoid matrix in homogeneous coordinates
          ) const;

      //! determine whether the generalized eigenproblem A*x = lambda B*x exhibits complex eigenvalues
      bool ComplexEigenValues(
          const LINALG::Matrix<4,4>&   A,   //!< first matrix of generalized eigenproblem
          const LINALG::Matrix<4,4>&   B    //!< second matrix of generalized eigenproblem
          ) const;

      //! compute contact point on particle surface
      void ComputeContactPoint(
          LINALG::Matrix<1,4>&         C_particle,   //!< contact point on particle surface
          const LINALG::Matrix<3,1>&   semiaxes,     //!< particle semi-axes
          const LINALG::Matrix<4,4>&   T,            //!< translation matrix
          const LINALG::Matrix<4,4>&   R,            //!< rotation matrix
          const LINALG::Matrix<1,4>&   C,            //!< contact point inside contact plane
          const LINALG::Matrix<1,4>&   d             //!< direction vector of contact line
          ) const;

      //! compute contact normal, contact gap, and vectors from particle centers to contact point
      void ComputeContactPointAndNormalAndGap(
          LINALG::Matrix<3,1>&      r_iC,     //!< vector from center of particle i to contact point C
          LINALG::Matrix<3,1>&      r_jC,     //!< vector from center of particle j to contact point C
          LINALG::Matrix<3,1>&      normal,   //!< contact normal
          double&                   gap,      //!< contact gap
          const ParticleCollData&   data_i,   //!< collision data for particle i
          const ParticleCollData&   data_j    //!< collision data for particle j
          ) const;

      //! compute contact normal, contact gap, contact point on wall element, and type of nearest contact object
      void ComputeContactPointAndNormalAndGapAndObjectType(
          LINALG::Matrix<3,1>&                        C,              //!< contact point on wall element
          GEO::ObjectType&                            objecttype,     //!< type of nearest contact object
          LINALG::Matrix<3,1>&                        normal,         //!< contact normal
          double&                                     gap,            //!< contact gap
          DRT::Element*                               wallele,        //!< wall element
          const std::map<int,LINALG::Matrix<3,1> >&   nodecoord,      //!< node coordinates of wall element
          const ParticleCollData&                     data            //!< particle collision data
          ) const;

      //! compute eigenvectors of B^-1 * A and determine whether there are complex eigenvalues
      bool ComputeEigenVectors(
          const LINALG::Matrix<4,4>&   A,    //!< first matrix of generalized eigenproblem
          const LINALG::Matrix<4,4>&   B,    //!< second matrix of generalized eigenproblem
          LINALG::Matrix<3,1>&         V0,   //!< first eigenvector of generalized eigenproblem
          LINALG::Matrix<3,1>&         V1,   //!< second eigenvector of generalized eigenproblem
          LINALG::Matrix<3,1>&         V2,   //!< third eigenvector of generalized eigenproblem
          LINALG::Matrix<3,1>&         V3    //!< fourth eigenvector of generalized eigenproblem
          ) const;

      //! compute ellipsoid matrix expressed in homogeneous coordinates
      void ComputeEllipsoidMatrix(
          LINALG::Matrix<4,4>&         E,          //!< ellipsoid matrix
          const LINALG::Matrix<3,1>&   semiaxes,   //!< ellipsoid semi-axes
          const LINALG::Matrix<4,4>&   T,          //!< translation matrix
          const LINALG::Matrix<4,4>&   R           //!< rotation matrix
          ) const;

      //! compute rotation matrix expressed in homogeneous coordinates
      void ComputeRotationMatrix(
          LINALG::Matrix<4,4>&         R,      //!< rotation matrix
          const LINALG::Matrix<3,1>&   angle   //!< rotation angle vector
          ) const;

      //! compute matrix transforming an ellipsoid to a unit sphere centered at the origin of the coordinate system
      void ComputeTransformationMatrix(
          LINALG::Matrix<4,4>&         M,          //!< transformation matrix
          const LINALG::Matrix<3,1>&   semiaxes,   //!< ellipsoid semi-axes
          const LINALG::Matrix<4,4>&   T,          //!< translation matrix
          const LINALG::Matrix<4,4>&   R           //!< rotation matrix
          ) const;

      //! compute matrix transforming a unit sphere centered at the origin of the coordinate system back to an ellipsoid
      void ComputeTransformationMatrixInverse(
          LINALG::Matrix<4,4>&         M_inv,      //!< inverse transformation matrix
          const LINALG::Matrix<3,1>&   semiaxes,   //!< ellipsoid semi-axes
          const LINALG::Matrix<4,4>&   T,          //!< translation matrix
          const LINALG::Matrix<4,4>&   R           //!< rotation matrix
          ) const;

      //! compute translation matrix expressed in homogeneous coordinates
      void ComputeTranslationMatrix(
          LINALG::Matrix<4,4>&         T,             //!< translation matrix
          const LINALG::Matrix<3,1>&   displacement   //!< displacement vector
          ) const;
  }; // class ParticleCollisionHandlerDEMEllipsoids


  // event to store future inter-particle collision
  class Event
  {
  public:

     /// Constructor
    explicit Event(
      INPAR::PARTICLE::coltype coltype,
      double time,
      DRT::Node* part_1,
      DRT::Node* part_2
      ) : coltype(coltype),
      time(time),
      particle_1(part_1),
      particle_2(part_2)
      { ; }

    /// Destructor
    virtual ~Event(){ ; }

    INPAR::PARTICLE::coltype coltype;
    double time;
    DRT::Node* particle_1;
    DRT::Node* particle_2;

    struct Helper
    {
      /// operator to sort two events in increasing time
      bool operator()(Teuchos::RCP<Event> event1, Teuchos::RCP<Event> event2);
    };
  };

  // event to store future particle-wall collision
  class WallEvent : public PARTICLE::Event
  {
  public:

    /// Constructor
    explicit WallEvent(
      INPAR::PARTICLE::coltype coltype,
      double time,
      DRT::Node* part_1,
      DRT::Element* wall,
      LINALG::Matrix<3,1> wallcollpoint_pos,
      LINALG::Matrix<3,1> wallcollpoint_vel
      ) : Event(coltype, time, part_1, NULL),
      wall(wall),
      wallcollpoint_pos(wallcollpoint_pos),
      wallcollpoint_vel(wallcollpoint_vel)
      { ; }

    /// Destructor
    virtual ~WallEvent(){ ; }

    DRT::Element* wall;
    LINALG::Matrix<3,1> wallcollpoint_pos;
    LINALG::Matrix<3,1> wallcollpoint_vel;
  };

  /*====================================================================*/
  /*!
   * \brief Collision handling with the molecular dynamics like type
   *  assuming rigid spherical particles
   *
   * \author hammerl
   * \date 09/13
   */
  class ParticleCollisionHandlerMD : public PARTICLE::ParticleCollisionHandlerBase
  {
  public:

    //! @name Life
    //@{

    //! Constructor
    ParticleCollisionHandlerMD(
      Teuchos::RCP<DRT::Discretization> discret,
      Teuchos::RCP<PARTICLE::Algorithm> particlealgorithm,
      const Teuchos::ParameterList& particledynparams
      );

    //! Destructor
    virtual ~ParticleCollisionHandlerMD(){ ; }

    //@}

    //! @name Actions
    //@{

    /// compute series of collisions (inter-particle and particle-wall) during on time step, returns contact energy
    virtual double EvaluateParticleContact(
      double dt,
      Teuchos::RCP<Epetra_Vector> disn,
      Teuchos::RCP<Epetra_Vector> veln,
      Teuchos::RCP<Epetra_FEVector> f_structure,
      int timestep);

    /// contact for one event is evaluated
    void HandleCollision(
      Teuchos::RCP<Event> next_event,
      const double dt
      );

    /// returns time to inter-particle collision in event
    Teuchos::RCP<Event> ComputeCollisionWithParticle(
      DRT::Node* particle1,
      DRT::Node* particle2,
      const double remaining_dt
      );

    /// returns time to particle-wall collision in event
    Teuchos::RCP<WallEvent> ComputeCollisionWithWall(
      DRT::Node* particle1,
      DRT::Element* wall,
      const double dt
      );

    /// setup of initial event queue at time step begin
    void InitializeEventQueue(
      std::set<Teuchos::RCP<Event>, PARTICLE::Event::Helper>& eventqueue,
      const double dt
      );

    /// event queue is filled with new events based on particles that have just collided
    void SearchForNewCollisions(
      Teuchos::RCP<Event> event,
      std::set<Teuchos::RCP<Event>, Event::Helper>& eventqueue,
      const double dt
      );

    /// computes time to collision between a particle and an element for hard sphere
    /// particles: searches hierarchically element, edges, corners
    bool ComputeCollisionOfParticleWithWall(
        DRT::Element* wallele,
        const Epetra_SerialDenseMatrix& xyze_current,
        const Epetra_SerialDenseMatrix& xyze_final,
        const LINALG::Matrix<3,1>& particle_pos,
        const LINALG::Matrix<3,1>& particle_vel,
        const double radius,
        double& timetocollision,
        LINALG::Matrix<3,1>& wall_pos,
        LINALG::Matrix<3,1>& wall_vel,
        const double remaining_dt,
        const double dt);

    /// computes time to collision between a particle and an element for hard sphere
    /// particles (templated on distype of element) using FAD
    template<DRT::Element::DiscretizationType DISTYPE>
    bool ComputeCollisionOfParticleWithElementT_FAD(
        DRT::Element* wallele,
        const Epetra_SerialDenseMatrix& xyze_n,
        const Epetra_SerialDenseMatrix& xyze_final,
        const LINALG::Matrix<3,1>& position,
        const LINALG::Matrix<3,1>& velocity,
        const double radius,
        double& timetocollision,
        LINALG::Matrix<3,1>& wall_pos,
        LINALG::Matrix<3,1>& wall_vel,
        const double remaining_dt,
        const double dt,
        bool& checkedges);

    /// computes time to collision between a particle and an element for hard sphere particles
    bool ComputeCollisionOfParticleWithElement(
        DRT::Element* wallele,
        const Epetra_SerialDenseMatrix& xyze_current,
        const Epetra_SerialDenseMatrix& xyze_final,
        const LINALG::Matrix<3,1>& position,
        const LINALG::Matrix<3,1>& velocity,
        const double radius,
        double& timetocollision,
        LINALG::Matrix<3,1>& wall_pos,
        LINALG::Matrix<3,1>& wall_vel,
        const double remaining_dt,
        const double dt,
        bool& checkedges);

    /// computes time to collision between a particle and a line for hard sphere
    /// particles (templated on distype of line)
    template<DRT::Element::DiscretizationType DISTYPE>
    bool ComputeCollisionOfParticleWithLineT(
        const Epetra_SerialDenseMatrix& xyze_line_n,
        const Epetra_SerialDenseMatrix& xyze_line_final,
        const LINALG::Matrix<3,1>& particle_pos,
        const LINALG::Matrix<3,1>& particle_vel,
        const double radius,
        double& timetocollision,
        LINALG::Matrix<3,1>& wallcollpoint_pos,
        LINALG::Matrix<3,1>& wallcollpoint_vel,
        bool& checkcorners);

    /// computes time to collision between a particle and a line for hard sphere
    /// particles (templated on distype of line) using FAD
    template<DRT::Element::DiscretizationType DISTYPE>
    bool ComputeCollisionOfParticleWithLineT_FAD(
        const Epetra_SerialDenseMatrix& xyze_line_n,
        const Epetra_SerialDenseMatrix& xyze_line_final,
        const LINALG::Matrix<3,1>& particle_pos,
        const LINALG::Matrix<3,1>& particle_vel,
        const double radius,
        double& timetocollision,
        LINALG::Matrix<3,1>& wallcollpoint_pos,
        LINALG::Matrix<3,1>& wallcollpoint_vel,
        const double remaining_dt,
        const double dt,
        bool& checkcorners);

    /// computes time to collision between a particle and a line for hard sphere particles
    bool ComputeCollisionOfParticleWithLine(
        const DRT::Element::DiscretizationType distype,
        const Epetra_SerialDenseMatrix& xyze_line_n,
        const Epetra_SerialDenseMatrix& xyze_line_final,
        const LINALG::Matrix<3,1>& particle_pos,
        const LINALG::Matrix<3,1>& particle_vel,
        const double radius,
        double& timetocollision,
        LINALG::Matrix<3,1>& wallcollpoint_pos,
        LINALG::Matrix<3,1>& wallcollpoint_vel,
        const double remaining_dt,
        const double dt,
        bool& checkcorners);

    /// computes time to collision between a particle and a corner of an element
    bool ComputeCollisionOfParticleWithCorner(
        const Epetra_SerialDenseMatrix& xyze_corner_n,
        const Epetra_SerialDenseMatrix& xyze_corner_np,
        const LINALG::Matrix<3,1>& particle_pos,
        const LINALG::Matrix<3,1>& particle_vel,
        const double radius,
        double& timetocollision,
        LINALG::Matrix<3,1>& wallcollpoint_pos,
        LINALG::Matrix<3,1>& wallcollpoint_vel,
        const double remaining_dt,
        const double dt);

    /// necessary data for computing inter-particle collisions is returned
    void GetCollisionData(
      const DRT::Node* particle1,
      const DRT::Node* particle2,
      LINALG::Matrix<3,1>& pos_1,
      LINALG::Matrix<3,1>& pos_2,
      LINALG::Matrix<3,1>& vel_1,
      LINALG::Matrix<3,1>& vel_2,
      double& rad_1,
      double& rad_2,
      double& ddt_1,
      double& ddt_2
      );

    /// necessary data for computing particle collision is returned
    void GetCollisionData(
      const DRT::Node* particle,
      LINALG::Matrix<3,1>& pos,
      LINALG::Matrix<3,1>& vel,
      double& rad,
      double& ddt
      );

     //@}

   };  // class ParticleCollisionHandlerMD

  //! @name Container structs for faster access
  //@{

  /// struct to store wall contact point and additional information about corresponding wall ele
  struct WallContactPoint {
   int eleid;
   LINALG::Matrix<3,1> point;
   LINALG::Matrix<3,1> normal;
   double penetration;
   std::map<int,LINALG::Matrix<3,1> > nodalCoordinates;
   std::vector<int> lm;
   std::vector<int> lmowner;
  };

  /// struct to store particle data for faster access
  struct ParticleCollData {
    LINALG::Matrix<3,1> dis;
    LINALG::Matrix<3,1> vel;
    LINALG::Matrix<3,1> ang;
    LINALG::Matrix<3,1> angvel;
    LINALG::Matrix<3,1> semiaxes;
    double ddt;
    double rad;
    double mass;
    double density;
    std::vector<int> lm;
    int owner;
    int id;
  };

  //@}

}  // namespace PARTICLE

/*----------------------------------------------------------------------*/
#endif  // #ifndef PARTICLE_CONTACT_H
