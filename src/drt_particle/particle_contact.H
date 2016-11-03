/*----------------------------------------------------------------------*/
/*!
\file particle_contact.H

\brief Particle collision handling

\level 3

\maintainer Georg Hammerl
*/

/*----------------------------------------------------------------------*/
#ifndef PARTICLE_CONTACT_H
#define PARTICLE_CONTACT_H

/*----------------------------------------------------------------------*/
/* headers */
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_inpar/inpar_particle.H"
#include "particle_node.H"
#include "../drt_lib/drt_element.H"
#include <Epetra_FEVector.h>


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
  struct HeatSource;
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
        Teuchos::RCP<Epetra_Vector> radiusn,
        Teuchos::RCP<Epetra_Vector> mass,
        Teuchos::RCP<Epetra_Vector> densityn = Teuchos::null,
        Teuchos::RCP<Epetra_Vector> specEnthalpyn = Teuchos::null);

    /// compute collisions (inter-particle and particle-wall), returns contact energy
    virtual double EvaluateParticleContact(
      double dt,
      Teuchos::RCP<Epetra_Vector> sysvec1,
      Teuchos::RCP<Epetra_Vector> sysvec2,
      Teuchos::RCP<Epetra_Vector> sysvec3 = Teuchos::null
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

    //@}

  protected:

    //! @name Collision parameters
    //@{

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
    /// dismember radius in case of thermodynamics
    double r_dismember_;
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

    //! write system energy every given step
    int writeenergyevery_;

    //@}

    /// ID of actual processor in parallel
    int myrank_;

    //! @name helper classes
    //@{

    /// particle discretization
    Teuchos::RCP<DRT::Discretization> discret_;

    /// particle algorithm that is connected to the time integration
    Teuchos::RCP<PARTICLE::Algorithm> particle_algorithm_;

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
    Teuchos::RCP<Epetra_Vector> specEnthalpyDotn = Teuchos::null);

  /// calculate contact with neighboring particles
  virtual void CalcNeighboringParticlesContact(
    DRT::Node* particle_i,
    const ParticleCollData& data_i,
    const std::list<DRT::Node*>& neighboring_particles,
    const bool havepbc,
    const double dt,
    const Teuchos::RCP<Epetra_Vector>& f_contact,
    const Teuchos::RCP<Epetra_Vector>& m_contact
    );

  /// calculate contact with neighboring walls
  virtual void CalcNeighboringWallsContact(
    DRT::Node* particle_i,
    const ParticleCollData& data_i,
    const std::set<DRT::Element*>& neighboring_walls,
    const double dt,
    const Teuchos::RCP<DRT::Discretization>& walldiscret,
    const Teuchos::RCP<const Epetra_Vector>& walldisn,
    const Teuchos::RCP<const Epetra_Vector>& wallveln,
    const Teuchos::RCP<Epetra_Vector>& f_contact,
    const Teuchos::RCP<Epetra_Vector>& m_contact,
    const Teuchos::RCP<Epetra_FEVector>& f_structure
    );

  ///calculate contact with neighboring heat sources
  void CalcNeighboringHeatSourcesContact(
    DRT::Node* particle_i,
    const ParticleCollData& data_i,
    const Teuchos::RCP<std::set<Teuchos::RCP<HeatSource>, BINSTRATEGY::Less> > neighboring_heatSources,
    const Teuchos::RCP<Epetra_Vector>& specEnthalpyDotn);

  /// calculate normal contact force for single contact pair
  virtual void CalculateNormalContactForce(
    const double g,
    const double v_rel_normal,
    const double m_eff,
    double& normalcontactforce,
    const int owner_i,
    const int owner_j
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

  //@}

   };  // class ParticleCollisionHandlerDEM


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
      Teuchos::RCP<Epetra_Vector> specEnthalpyn = Teuchos::null);

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
   double penetration;
   std::map<int,LINALG::Matrix<3,1> > nodalCoordinates;
   std::vector<int> lm;
   std::vector<int> lmowner;
  };

  /// struct to store particle data for faster access
  struct ParticleCollData {
    LINALG::Matrix<3,1> dis;
    LINALG::Matrix<3,1> vel;
    LINALG::Matrix<3,1> angvel;
    double ddt;
    double rad;
    double mass;
    double density;
    double specEnthalpy;
    std::vector<int> lm;
    int owner;
  };

  //@}

}  // namespace PARTICLE

/*----------------------------------------------------------------------*/
#endif  // #ifndef PARTICLE_CONTACT_H
