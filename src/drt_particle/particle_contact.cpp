/*----------------------------------------------------------------------*/
/*!
\file particle_contact.cpp

\brief Particle collision handling

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*-----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "particle_contact.H"
#include "particle_algorithm.H"
#include "particle_timint_centrdiff.H"
#include "particle_heatSource.H"
#include "../drt_adapter/ad_str_structure.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_mat/particle_mat.H"
#include "../drt_mat/particle_mat_ellipsoids.H"
#include "../drt_mat/extparticle_mat.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_binstrategy/drt_meshfree_multibin.H"
#include "../drt_geometry/searchtree_geometry_service.H"
#include "../drt_geometry/position_array.H"
#include "../drt_geometry/element_coordtrafo.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/largerotations.H"
#include "../drt_io/io_control.H"
#include "../headers/compiler_definitions.h"

#include "../drt_mat/stvenantkirchhoff.H"

#include <Teuchos_TimeMonitor.hpp>
#include <Sacado.hpp>

typedef Sacado::Fad::DFad<double> FAD;

//#define OUTPUT

/*----------------------------------------------------------------------*
 | constructor for particle contact                        ghamm 09/13  |
 *----------------------------------------------------------------------*/

PARTICLE::ParticleCollisionHandlerBase::ParticleCollisionHandlerBase(
  Teuchos::RCP<DRT::Discretization> discret,
  Teuchos::RCP<PARTICLE::Algorithm> particlealgorithm,
  const Teuchos::ParameterList& particledynparams
  ) :
  normal_contact_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::NormalContact>(particledynparams,"NORMAL_CONTACT_LAW")),
  normal_adhesion_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::NormalAdhesion>(particledynparams,"NORMAL_ADHESION_LAW")),
  nue_(0.),
  young_(0.),
  r_min_(particledynparams.get<double>("MIN_RADIUS")),
  r_max_(particledynparams.get<double>("MAX_RADIUS")),
  r_dismember_(0.),
  v_max_(particledynparams.get<double>("MAX_VELOCITY")),
  c_(particledynparams.get<double>("REL_PENETRATION")),
  c_wall_(particledynparams.get<double>("REL_PENETRATION_WALL")),
  dt_krit_(0.),
  e_(particledynparams.get<double>("COEFF_RESTITUTION")),
  e_wall_(particledynparams.get<double>("COEFF_RESTITUTION_WALL")),
  k_normal_(particledynparams.get<double>("NORMAL_STIFF")),
  k_normal_wall_(particledynparams.get<double>("NORMAL_STIFF_WALL")),
  k_tang_(0.),
  k_tang_wall_(0.),
  kappa_(0.),
  kappa_wall_(0.),
  d_normal_(0.),
  d_tang_(0.),
  d_normal_wall_(0.),
  d_tang_wall_(0.),
  mu_(particledynparams.get<double>("FRICT_COEFF")),
  mu_wall_(particledynparams.get<double>("FRICT_COEFF_WALL")),
  tension_cutoff_(DRT::INPUT::IntegralValue<int>(particledynparams,"TENSION_CUTOFF") == 1),
  contact_energy_(0.),
  g_max_(0.),
  adhesion_eq_gap_(particledynparams.get<double>("ADHESION_EQ_GAP")),
  adhesion_normal_stiff_(particledynparams.get<double>("ADHESION_NORMAL_STIFF")),
  adhesion_normal_damp_(particledynparams.get<double>("ADHESION_NORMAL_DAMP")),
  adhesion_normal_eps_(particledynparams.get<double>("ADHESION_NORMAL_EPS")),
  adhesion_max_force_(particledynparams.get<double>("ADHESION_MAX_FORCE")),
  adhesion_max_disp_(particledynparams.get<double>("ADHESION_MAX_DISP")),
  writeenergyevery_(particledynparams.get<int>("RESEVRYERGY")),
  myrank_(discret->Comm().MyPID()),
  discret_(discret),
  particle_algorithm_(particlealgorithm),
  particleData_(0)
{
  // extract the material
  const MAT::PAR::ParticleMat* particleMat = particle_algorithm_->ParticleMat();
  // currently all particles have identical density and radius
  double density = particleMat->initDensity_;
  nue_ = particleMat->poissonRatio_;
  young_ = particleMat->youngModulus_;

  if(particle_algorithm_->ParticleInteractionType() != INPAR::PARTICLE::None)
  {
    // safety checks for particle-particle contact
    if(r_min_<0.0 or r_max_<0.0 or v_max_<0.0 or young_<0.0)
      dserror("Invalid input parameter (MIN_RADIUS, MAX_RADIUS, MAX_VELOCITY, YOUNG's modulus have to be larger than zero)");
    if(nue_<=-1.0 or nue_>0.5)
      dserror("Invalid input parameter (NUE has to be in the range ]-1.0;0.5])");
    if(not((c_ <= 0. and k_normal_ > 0.) or (c_ > 0. and v_max_ > 0. and k_normal_ <= 0.)))
      dserror("For particle-particle contact, you must correctly specify either the relative penetration (along with the maximum velocity), or the normal stiffness, but neither both nor none of them!");
    if(r_min_>r_max_)
      dserror("inversed radii (MIN_RADIUS > MAX_RADIUS)");
    if (particleMat->initRadius_ < r_min_)
      dserror("INITRADIUS too small (it should be >= MIN_RADIUS)");
    if (particleMat->initRadius_ > r_max_)
      dserror("INITRADIUS too big (it should be <= MAX_RADIUS)");
    if(e_<0.0 and (normal_contact_ == INPAR::PARTICLE::LinSpringDamp or particle_algorithm_->ParticleInteractionType()==INPAR::PARTICLE::Normal_MD))
      dserror("Invalid input parameter COEFF_RESTITUTION for this kind of contact law!");

    // no further information necessary for MD like contact
    if(particle_algorithm_->ParticleInteractionType() == INPAR::PARTICLE::Normal_MD)
      return;

    // compute minimum particle mass
    double mass_min(0.);
    const MAT::PAR::ParticleMatEllipsoids* const particlematellipsoids = dynamic_cast<const MAT::PAR::ParticleMatEllipsoids* const>(particleMat);
    if(particlematellipsoids == NULL)
      mass_min = density * 4.0/3.0 * M_PI * pow( r_min_ ,3.0 );
    else
      mass_min = density * 4.0/3.0 * M_PI * particlematellipsoids->SemiAxes()(0) * particlematellipsoids->SemiAxes()(1) * particlematellipsoids->SemiAxes()(2);

    // data for critical time step computation for DEM like contact
    double k_tkrit = 0.0;

    const double G = young_ / (2.0*(1.0+nue_));

    // kappa - tangential to normal stiffness ratio
    kappa_ = (1.0-nue_)/(1.0-0.5*nue_);

    if(particle_algorithm_->WallDiscret() != Teuchos::null)
    {
      // safety checks for particle-wall contact
      if(not((c_wall_ <= 0. and k_normal_wall_ > 0.) or (c_wall_ > 0. and v_max_ > 0. and k_normal_wall_ <= 0.)))
        dserror("For particle-wall contact, you must correctly specify either the relative penetration (along with the maximum velocity), or the normal stiffness, but neither both nor none of them!");
      if(e_wall_<0.0 and (normal_contact_ == INPAR::PARTICLE::LinSpringDamp or particle_algorithm_->ParticleInteractionType()==INPAR::PARTICLE::Normal_MD))
        dserror("Invalid input parameter COEFF_RESTITUTION_WALL for this kind of contact law!");

      // wall material properties: young's modulus, poisson's ratio and shear modulus
      const double young_wall = particledynparams.get<double>("YOUNG_WALL");
      const double nue_wall = particledynparams.get<double>("NUE_WALL");
      if ( young_wall<0.0 or nue_wall<=-1.0 or nue_wall>0.5 )
        dserror("Invalid input parameter (YOUNG_WALL has to be larger than zero, NUE_WALL has to be in the range ]-1.0;0.5])");
      const double G_wall = 0.5*young_wall/(1.0+nue_wall);

      if(G_wall<0.0)
        dserror("Wall shear modulus has to be larger than zero");

      // kappa - tangential to normal stiffness ratio
      kappa_wall_ = ( (1.0-nue_)/G + (1.0-nue_wall)/G_wall ) / ( (1.0-0.5*nue_)/G + (1.0-0.5*nue_wall)/G_wall );
    }

    //------------stiffness----------------------------
    switch(normal_contact_)
    {
    case INPAR::PARTICLE::LinSpring:
    case INPAR::PARTICLE::LinSpringDamp:
    {
      // calculate normal stiffness from relative penetration and other input parameters if necessary
      if(c_ > 0.)
        k_normal_ = 2./3. * r_max_ * M_PI * density * v_max_ * v_max_ / (c_ * c_);
      if(c_wall_ > 0.)
        k_normal_wall_ = 2./3. * r_max_ * M_PI * density * v_max_ * v_max_ / (c_wall_ * c_wall_);

      // for tangential contact same stiffness is used
      k_tang_ = kappa_ * k_normal_;
      k_tang_wall_ = kappa_wall_ * k_normal_wall_;

      // critical time step size is calculated based on particle-particle (not particle-wall) contact
      k_tkrit = k_normal_;

      break;
    }

    case INPAR::PARTICLE::Hertz:
    case INPAR::PARTICLE::LeeHerrmann:
    case INPAR::PARTICLE::KuwabaraKono:
    case INPAR::PARTICLE::Tsuji:
    {
      if(particle_algorithm_->ParticleInteractionType()==INPAR::PARTICLE::NormalAndTang_DEM)
        dserror("tangential contact only with linear normal model implemented");

      if(c_ > 0.)
      {
        // calculate normal stiffness from relative penetration and other input parameters if necessary
        k_normal_ = 10./3. * M_PI * density * v_max_ * v_max_ * pow(r_max_,0.5) / pow(2.*c_,2.5);

        // critical time step size is calculated based on particle-particle (not particle-wall) contact
        k_tkrit = 2./3. * r_max_ * M_PI * density * v_max_ * v_max_ / (c_ * c_);
      }
      else
        // critical time step size is calculated based on particle-particle (not particle-wall) contact
        k_tkrit = pow(2048./1875. * density * v_max_ * v_max_ * M_PI * pow(r_max_,3) * pow(k_normal_,4),0.2);

      if(c_wall_ > 0.)
        // calculate normal stiffness from relative penetration and other input parameters if necessary
        k_normal_wall_ = 10./3. * M_PI * density * v_max_ * v_max_ * pow(r_max_,0.5) / pow(2.*c_wall_,2.5);

      // for tangential contact the user specified (nonlinear) normal stiffness which has to be transformed into a linear normal
      // stiffness with the same relative penetration which is used as (linear) tangential stiffness afterwards

      break;
    }

    default:
    {
      dserror("Normal contact law does not exist!");
      break;
    }
    }
    //---------------------------------------------------------

    //------------------damping--------------------------------
    d_normal_ = -1.0;
    d_tang_ = -1.0;

    if(normal_contact_ == INPAR::PARTICLE::LinSpringDamp)
    {
      double user_normal_damping = particledynparams.get<double>("NORMAL_DAMP");
      if(user_normal_damping >= 0.0)
      {
        dserror("Invalid input parameter NORMAL_DAMP for this kind of contact law");
      }
      double user_tang_damping = particledynparams.get<double>("TANG_DAMP");
      if(user_tang_damping >= 0.0)
      {
        dserror("Invalid input parameter TANG_DAMP for this kind of contact law");
      }
    }

    if(normal_contact_ == INPAR::PARTICLE::LeeHerrmann || normal_contact_ == INPAR::PARTICLE::KuwabaraKono ||
        normal_contact_ == INPAR::PARTICLE::Tsuji)
    {
      double user_normal_damping = particledynparams.get<double>("NORMAL_DAMP");
      if(user_normal_damping >= 0.0)
      {
        // user has to specify normal damping coefficient
        d_normal_ = user_normal_damping;
      }
      else
      {
        dserror("For this kind of contact law the input parameter NORMAL_DAMP is invalid");
      }
      double user_tang_damping = particledynparams.get<double>("TANG_DAMP");
      if(user_tang_damping >= 0.0)
      {
        // user has to specify tangential damping coefficient
        d_tang_ = user_tang_damping;
      }
      else
      {
        if(particle_algorithm_->ParticleInteractionType()==INPAR::PARTICLE::NormalAndTang_DEM)
          dserror("For this kind of contact law the input parameter TANG_DAMP is invalid");
      }
    }
    //------------------damping (wall)--------------------------------
    if(particle_algorithm_->WallDiscret() != Teuchos::null)
    {
      d_normal_wall_ = -1.0;
      d_tang_wall_ = -1.0;

      if(normal_contact_ == INPAR::PARTICLE::LinSpringDamp)
      {
        double user_normal_damping = particledynparams.get<double>("NORMAL_DAMP_WALL");
        if(user_normal_damping >= 0.0)
        {
          dserror("Invalid input parameter NORMAL_DAMP_WALL for this kind of contact law");
        }
        double user_tang_damping = particledynparams.get<double>("TANG_DAMP_WALL");
        if(user_tang_damping >= 0.0)
        {
          dserror("Invalid input parameter TANG_DAMP_WALL for this kind of contact law");
        }
      }

      if(normal_contact_ == INPAR::PARTICLE::LeeHerrmann || normal_contact_ == INPAR::PARTICLE::KuwabaraKono ||
          normal_contact_ == INPAR::PARTICLE::Tsuji)
      {
        double user_normal_damping = particledynparams.get<double>("NORMAL_DAMP_WALL");
        if(user_normal_damping >= 0.0)
        {
          // user has to specify normal damping coefficient
          d_normal_wall_ = user_normal_damping;
        }
        else
        {
          dserror("For this kind of contact law the input parameter NORMAL_DAMP_WALL is invalid");
        }
        double user_tang_damping = particledynparams.get<double>("TANG_DAMP_WALL");
        if(user_tang_damping >= 0.0)
        {
          // user has to specify tangential damping coefficient
          d_tang_wall_ = user_tang_damping;
        }
        else
        {
          if(particle_algorithm_->ParticleInteractionType()==INPAR::PARTICLE::NormalAndTang_DEM)
            dserror("For this kind of contact law the input parameter TANG_DAMP_WALL is invalid");
        }
      }
    }
    //---------------------------------------------------------------

    double factor = 0.0;
    const double safety = 0.75;

    // initialize factor
    if(particle_algorithm_->ParticleInteractionType()==INPAR::PARTICLE::Normal_DEM || particle_algorithm_->ParticleInteractionType()==INPAR::PARTICLE::Normal_DEM_thermo)
    { factor = 0.34; }
    if(particle_algorithm_->ParticleInteractionType()==INPAR::PARTICLE::NormalAndTang_DEM)
    { factor = 0.22; }

    // calculate critical time step for DEM like contact
    dt_krit_ = safety * factor * sqrt( mass_min / k_tkrit );

    double dt = particledynparams.get<double>("TIMESTEP");
    if(dt>dt_krit_ and myrank_==0)
    {
      std::cout << "\n\nW A R N I N G : time step larger than critical time step!" << std::endl;
      std::cout << "W A R N I N G : calculated critical time step: "<<dt_krit_<<" !\n\n" << std::endl;
    }

    // check frictional coefficient
    if(particle_algorithm_->ParticleInteractionType() == INPAR::PARTICLE::NormalAndTang_DEM)
    {
      if(mu_<=0.0 or (particle_algorithm_->WallDiscret() != Teuchos::null and mu_wall_<=0.0))
       dserror("Friction coefficient invalid");
    }
  }
  const MAT::PAR::ExtParticleMat* extParticleMat = particle_algorithm_->ExtParticleMat();
  if (extParticleMat != NULL)
  {
    r_dismember_ = extParticleMat->dismemberRadius_;
    if (r_dismember_ < 0)
      dserror("Invalid or unset input parameter (DISMEMBER_RADIUS must be larger than zero)");
    if (r_dismember_ < r_min_)
      dserror("DISMEMBER_RADIUS < MIN_RADIUS -> DISMEMBER_RADIUS is too small!");

    if (r_dismember_ > r_max_)
            dserror("DISMEMBER_RADIUS > MAX_RADIUS -> DISMEMBER_RADIUS is too big!");
  }

  // safety checks for particle adhesion
  switch(normal_adhesion_)
  {
    case INPAR::PARTICLE::adhesion_none:
    {
      if(std::abs(adhesion_eq_gap_) > GEO::TOL14
          or adhesion_normal_stiff_ >= 0.
          or adhesion_normal_damp_ >= 0.
          or adhesion_normal_eps_ >= 0.
          or adhesion_max_force_ >= 0.
          or adhesion_max_disp_ >= 0.)
        dserror("If you do not want to consider particle adhesion, please do not specify any associated parameters in the input file!");

      break;
    }

    case INPAR::PARTICLE::adhesion_lennardjones:
    {
      if(adhesion_normal_stiff_ >= 0.
          or adhesion_normal_damp_ >= 0.
          or adhesion_normal_eps_ <= 0.)
        dserror("Invalid particle adhesion parameters specified in input file for Lennard-Jones potential!");

      break;
    }

    case INPAR::PARTICLE::adhesion_linspring:
    {
      if(adhesion_normal_stiff_ <= 0.
          or adhesion_normal_damp_ >= 0.
          or adhesion_normal_eps_ >= 0.)
        dserror("Invalid particle adhesion parameters specified in input file for linear-spring model!");

      break;
    }

    case INPAR::PARTICLE::adhesion_linspringdamp:
    {
      if(adhesion_normal_stiff_ <= 0.
          or adhesion_normal_damp_ <= 0.
          or adhesion_normal_eps_ >= 0.)
        dserror("Invalid particle adhesion parameters specified in input file for linear-spring-damp model!");

      break;
    }

    default:
    {
      dserror("Unknown adhesion law governing normal contact of particles!");
      break;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | set states from time integrator to prepare collisions   catta 09/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerBase::Init(
    Teuchos::RCP<Epetra_Vector> disn,
    Teuchos::RCP<Epetra_Vector> veln,
    Teuchos::RCP<Epetra_Vector> angVeln,
    Teuchos::RCP<const Epetra_Vector> radiusn,
    Teuchos::RCP<Epetra_Vector> orientn,
    Teuchos::RCP<Epetra_Vector> mass,
    Teuchos::RCP<Epetra_Vector> densityn,
    Teuchos::RCP<Epetra_Vector> specEnthalpyn)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleCollisionHandlerBase::ContactInit");
  // export everything in col layout

  // dof based vectors
  Teuchos::RCP<Epetra_Vector> disnCol = LINALG::CreateVector(*discret_->DofColMap(),false);
  Teuchos::RCP<Epetra_Vector> velnCol = LINALG::CreateVector(*discret_->DofColMap(),false);
  Teuchos::RCP<Epetra_Vector> angVelnCol = LINALG::CreateVector(*discret_->DofColMap(),false);
  Teuchos::RCP<Epetra_Vector> orientnCol(Teuchos::null);
  if(orientn != Teuchos::null)
    orientnCol = LINALG::CreateVector(*discret_->DofColMap(),false);

  // setup importer for dof based vectors once in the beginning and reuse it
  Epetra_Import dofimporter(*discret_->DofColMap(), *discret_->DofRowMap());
  int err = 0;
  err += disnCol->Import(*disn, dofimporter, Insert);
  err += velnCol->Import(*veln, dofimporter, Insert);
  err += angVelnCol->Import(*angVeln, dofimporter, Insert);
  if(orientn != Teuchos::null)
    err += orientnCol->Import(*orientn, dofimporter, Insert);
  if (err)
    dserror("Export using importer failed for dof based Epetra_Vector: return value != 0");

  // node based vector
  Teuchos::RCP<Epetra_Vector> massCol = LINALG::CreateVector(*discret_->NodeColMap(),false);

  // setup importer for node based vector once in the beginning and reuse it
  Epetra_Import nodeimporter(*discret_->NodeColMap(), *discret_->NodeRowMap());
  err += massCol->Import(*mass, nodeimporter, Insert);

  // export radius vector
  bool radius_nodebased(true);
  Teuchos::RCP<Epetra_Vector> radiusnCol(Teuchos::null);
  if(radiusn->Map().SameAs(*discret_->NodeRowMap()))
  {
    radiusnCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
    err += radiusnCol->Import(*radiusn, nodeimporter, Insert);
  }
  else if(radiusn->Map().SameAs(*discret_->DofRowMap()))
  {
    radius_nodebased = false;
    radiusnCol = LINALG::CreateVector(*discret_->DofColMap(),false);
    err += radiusnCol->Import(*radiusn, dofimporter, Insert);
  }
  else
    dserror("You have a very strange radius vector...");

  Teuchos::RCP<Epetra_Vector> densityCol, specEnthalpyCol;
  if (densityn != Teuchos::null)
  {
    densityCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
    specEnthalpyCol = LINALG::CreateVector(*discret_->NodeColMap(),false);
    err += densityCol->Import(*densityn, nodeimporter, Insert);
    err += specEnthalpyCol->Import(*specEnthalpyn, nodeimporter, Insert);
  }

  if (err)
    dserror("Export using importer failed for node based Epetra_Vector: return value != 0");

  // fill particleData_
  const int numcolparticles = discret_->NodeColMap()->NumMyElements();
  particleData_.resize(numcolparticles);

  for (int i=0; i<numcolparticles; ++i)
  {
    // particle for which data will be collected
    DRT::Node *particle = discret_->lColNode(i);

    std::vector<int> lm;
    lm.reserve(3);

    // extract global dof ids and fill into lm_i
    discret_->Dof(particle, lm);

    ParticleCollData& data = particleData_[particle->LID()];

    //position, velocity and angular velocity of particle
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*disnCol,data.dis,lm);
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*velnCol,data.vel,lm);
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*angVelnCol,data.angvel,lm);
    if(orientn != Teuchos::null)
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*orientnCol,data.ang,lm);

    const int lid = particle->LID();
    if(radius_nodebased)
      data.rad = (*radiusnCol)[lid];
    else
    {
      data.rad = r_max_;
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*radiusnCol,data.semiaxes,lm);
    }
    data.mass = (*massCol)[lid];
    if (densityn != Teuchos::null)
    {
      data.density = (*densityCol)[lid];
      data.specEnthalpy = (*specEnthalpyCol)[lid];
    }

    // set ddt
    data.ddt = 0;

    // lm vector and owner
    data.lm.swap(lm);
    data.owner = particle->Owner();

    // particle ID
    data.id = particle->Id();

#ifdef DEBUG
  if(particle->NumElement() != 1)
    dserror("More than one element for this particle");
#endif
  }
}

/*----------------------------------------------------------------------*
 | assemble energies of particles                          ghamm 09/13  |
 *----------------------------------------------------------------------*/
double PARTICLE::ParticleCollisionHandlerBase::EnergyAssemble(double owner_i,double owner_j)
{
  // initialize scaling factor to unity
  double value = 1.0;

  // extract ID of current processor
  const int myrank = discret_->Comm().MyPID();

  // adjust scaling factor in parallel case if necessary
  // contact with wall
  if(owner_j<0)
  {
    if(owner_i!=myrank)
     value = 0.0;
  }
  // contact between particle_i and particle_j
  else
  {
    if((owner_i==myrank && owner_j!=myrank) or (owner_i!=myrank && owner_j==myrank))
      value = 0.5;
    else if(owner_i!=myrank and owner_j!=myrank)
      value = 0.0;
  }

  // return scaling factor
  return value;
}


/*----------------------------------------------------------------------*
 | constructor for DEM based particle contact              ghamm 09/13  |
 *----------------------------------------------------------------------*/
PARTICLE::ParticleCollisionHandlerDEM::ParticleCollisionHandlerDEM(
  Teuchos::RCP<DRT::Discretization> discret,
  Teuchos::RCP<PARTICLE::Algorithm> particlealgorithm,
  const Teuchos::ParameterList& particledynparams
  ) :
  PARTICLE::ParticleCollisionHandlerBase(
    discret,
    particlealgorithm,
    particledynparams
    ),
  writecontactforcesevery_(particledynparams.get<int>("CONTACTFORCESEVRY"))
{
  return;
}


/*----------------------------------------------------------------------*
 | compute collisions (inter-particle and particle-wall)   ghamm 09/13  |
 *----------------------------------------------------------------------*/
double PARTICLE::ParticleCollisionHandlerDEM::EvaluateParticleContact(
  const double dt,
  Teuchos::RCP<Epetra_Vector> f_contact,
  Teuchos::RCP<Epetra_Vector> m_contact,
  Teuchos::RCP<Epetra_Vector> specEnthalpyDotn,
  Teuchos::RCP<Epetra_FEVector> f_structure)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleCollisionHandlerDEM::ContactSearchAndCalculation");

  // initialize internal variables
  contact_energy_ = 0.0;
  g_max_ = 0.0;

  // get wall discretization and states for particles
  Teuchos::RCP<DRT::Discretization> walldiscret = particle_algorithm_->WallDiscret();
  Teuchos::RCP<const Epetra_Vector> walldisn(Teuchos::null);
  Teuchos::RCP<const Epetra_Vector> wallveln(Teuchos::null);
  if(walldiscret != Teuchos::null)
  {
    walldisn = walldiscret->GetState("walldisnp");
    wallveln = walldiscret->GetState("wallvelnp");
  }

  const bool havepbc = particle_algorithm_->BinStrategy().HavePBC();

  // store bins, which have already been examined
  std::set<int> examinedbins;

  // list of all particles in the neighborhood of currparticle
  std::list<DRT::Node*> neighboring_particles;

  // loop over all row particles
  const int numrowparticles = discret_->NodeRowMap()->NumMyElements();
  for(int rowPar_i=0; rowPar_i<numrowparticles; ++rowPar_i)
  {
    // extract the particle
    DRT::Node *currparticle = discret_->lRowNode(rowPar_i);

    //find the bin it belongs to
    DRT::Element* currentBin = currparticle->Elements()[0];

    // check I own the bin
    assert(currentBin->Owner() == myrank_);

    const int binId = currentBin->Id();
    // if a bin has already been examined --> continue with next particle
    if( examinedbins.find(binId) != examinedbins.end() )
      continue;
    //else: bin is examined for the first time --> new entry in examinedbins_
    examinedbins.insert(binId);

    // extract the pointer to the particles
    DRT::Node** currentBinParticles = currentBin->Nodes();

    // remove current content but keep memory
    neighboring_particles.clear();

    // list of walls that border on the CurrentBin
    boost::unordered_map<int, DRT::Element*> neighboring_walls;

    // list of heat sources that border on the CurrentBin
    const Teuchos::RCP<boost::unordered_map<int , Teuchos::RCP<HeatSource> > > neighboring_heatSources = particle_algorithm_->ParticleInteractionType() == INPAR::PARTICLE::Normal_DEM_thermo ? Teuchos::rcp(new boost::unordered_map<int , Teuchos::RCP<HeatSource> >) : Teuchos::null;

    particle_algorithm_->GetNeighbouringItems(binId, neighboring_particles, &neighboring_walls, neighboring_heatSources);

    // loop over all particles in CurrentBin
    for(int i=0; i<currentBin->NumNode(); ++i)
    {
      // determine the particle we are analizing
      DRT::Node* particle_i = currentBinParticles[i];

      // extract data
      ParticleCollData& data_i = particleData_[particle_i->LID()];

      // compute contact with neighboring walls
      CalcNeighboringWallsContact(particle_i, data_i, neighboring_walls, dt,
          walldiscret, walldisn, wallveln, f_contact, m_contact, f_structure);

      if (f_contact != Teuchos::null and m_contact != Teuchos::null)
      {
        // compute contact with neighboring particles
        CalcNeighboringParticlesContact(particle_i, data_i, neighboring_particles,
            havepbc, dt, f_contact, m_contact, specEnthalpyDotn);

        if (particle_algorithm_->ParticleInteractionType() == INPAR::PARTICLE::Normal_DEM_thermo)
        {
          CalcNeighboringHeatSourcesContact(data_i, neighboring_heatSources, specEnthalpyDotn);
        }
      }
    }
  }

  // communicate maximum particle-particle or particle-wall penetration
  double g_max(g_max_);
  particle_algorithm_->Comm().MaxAll(&g_max,&g_max_,1);

  // gather *.csv files with normal particle-particle and particle-wall contact forces across all processors if applicable
  if(writecontactforcesevery_ and particle_algorithm_->Step() % writecontactforcesevery_ == 0)
    GatherNormalContactForcesToFile();

  // erase temporary storage for collision data
  particleData_.clear();

  return contact_energy_;
}


/*----------------------------------------------------------------------*
 | calculate contact with neighboring heat sources         katta 10/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerDEM::CalcNeighboringHeatSourcesContact(
  const ParticleCollData& data_i,
  const Teuchos::RCP<boost::unordered_map<int, Teuchos::RCP<HeatSource> > > neighboring_heatSources,
  const Teuchos::RCP<Epetra_Vector>& specEnthalpyDotn)
{
  double specEnthalpyDot_i = 0.0;
  boost::unordered_map<int , Teuchos::RCP<HeatSource> >::const_iterator ii;
  for(ii = neighboring_heatSources->begin(); ii != neighboring_heatSources->end();  ++ii)
  {
    const Teuchos::RCP<HeatSource> hs = ii->second;
    if (hs->minVerZone_[0]<=data_i.dis(0) &&
        hs->minVerZone_[1]<=data_i.dis(1) &&
        hs->minVerZone_[2]<=data_i.dis(2) &&
        hs->maxVerZone_[0]>=data_i.dis(0) &&
        hs->maxVerZone_[1]>=data_i.dis(1) &&
        hs->maxVerZone_[2]>=data_i.dis(2))
    {
      specEnthalpyDot_i += (hs->QDot_)/data_i.density;
    }
  }

  LINALG::Assemble(*specEnthalpyDotn, specEnthalpyDot_i, data_i.id, data_i.owner);
}


/*----------------------------------------------------------------------*
 | calculate contact with neighboring particles             ghamm 04/16 |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerDEM::CalcNeighboringParticlesContact(
    DRT::Node*                           particle_i,              //!< i-th particle
    const ParticleCollData&              data_i,                  //!< collision data associated with i-th particle
    const std::list<DRT::Node*>&         neighboring_particles,   //!< list of particles adjacent to i-th particle
    const bool                           havepbc,                 //!< flag indicating periodic boundaries
    const double                         dt,                      //!< time step size
    const Teuchos::RCP<Epetra_Vector>&   f_contact,               //!< global force vector
    const Teuchos::RCP<Epetra_Vector>&   m_contact,               //!< global moment vector
    const Teuchos::RCP<Epetra_Vector>&   specEnthalpyDotn         //!< global enthalpy vector
    )
{
  std::map<int, PARTICLE::Collision>& history_particle_i = static_cast<PARTICLE::ParticleNode*>(particle_i)->Get_history_particle();
  if(history_particle_i.size() > 20)
    dserror("Particle i is in contact with more than 20 particles! Check whether history is deleted correctly!");

  // check whether there is contact between particle i and all other particles in the neighborhood except those which
  // have a lower or equal ID than particle i (--> ignoring self-contact)
  for(std::list<DRT::Node*>::const_iterator j=neighboring_particles.begin(); j!=neighboring_particles.end(); ++j)
  {
    DRT::Node* neighborparticle = (*j);

    // extract data
    ParticleCollData& data_j = particleData_[neighborparticle->LID()];
    std::map<int, PARTICLE::Collision>& history_particle_j = static_cast<PARTICLE::ParticleNode*>(neighborparticle)->Get_history_particle();
    if(history_particle_j.size() > 20)
      dserror("Particle j is in contact with more than 20 particles! Check whether history is deleted correctly!");

    // evaluate contact only once in case we own particle j. Otherwise compute everything,
    // the assemble method does not write in case the particle is a ghost
    if(data_i.id < data_j.id || data_j.owner != myrank_)
    {
      // normalized mass
      const double m_eff = data_i.mass * data_j.mass / (data_i.mass + data_j.mass);

      // might need to shift position of particle j in the presence of periodic boundary conditions
      if(havepbc)
      {
        static int ijk_i[3], ijk_j[3];
        particle_algorithm_->BinStrategy().ConvertPosToijk(data_i.dis,ijk_i);
        particle_algorithm_->BinStrategy().ConvertPosToijk(data_j.dis,ijk_j);
        for(unsigned idim=0; idim<3; ++idim)
        {
          if(particle_algorithm_->BinStrategy().HavePBC(idim))
          {
            if(ijk_i[idim] - ijk_j[idim] < -1)
              data_j.dis(idim) -= particle_algorithm_->BinStrategy().PBCDelta(idim);
            else if(ijk_i[idim] - ijk_j[idim] > 1)
              data_j.dis(idim) += particle_algorithm_->BinStrategy().PBCDelta(idim);
          }
        }
      }

      // compute contact normal, contact gap, and vectors from particle centers to contact point C
      static LINALG::Matrix<3,1> r_iC, r_jC, normal;
      double g(0.);
      ComputeContactPointAndNormalAndGap(r_iC,r_jC,normal,g,data_i,data_j);

      // in case of penetration or adhesion, contact forces and moments are calculated
      if(g <= 0.0 or ConsiderNormalAdhesion(g))
      {
        // contact forces
        double normalcontactforce = 0.0;
        static LINALG::Matrix<3,1> tangentcontactforce(true);

        if(g<=0. and std::abs(g)>g_max_)
          g_max_ = std::abs(g);

        // relative velocity in contact point between particles i and j
        static LINALG::Matrix<3,1> v_rel;
        ComputeRelativeVelocity(v_rel,r_iC,r_jC,data_i,data_j);

        // part of v_rel in normal direction: v_rel * n
        const double v_rel_normal(v_rel.Dot(normal));

        // calculation of normal contact force
        CalculateNormalContactForce(g, data_i.rad, data_j.rad, v_rel_normal, m_eff, normalcontactforce, data_i.owner, data_j.owner);

        // output normal contact force between the two particles to *.csv file if applicable
        if(writecontactforcesevery_ and particle_algorithm_->Step() % writecontactforcesevery_ == 0 and data_i.owner == myrank_ and data_j.owner <= myrank_)
        {
          LINALG::Matrix<3,1> forceapplicationpoint(true);
          forceapplicationpoint.Update(1.,data_i.dis,data_i.rad+0.5*g,normal);
          OutputNormalContactForceToFile(normalcontactforce,forceapplicationpoint);
        };

        // calculation of tangential contact force
        if(g <= 0. and particle_algorithm_->ParticleInteractionType() == INPAR::PARTICLE::NormalAndTang_DEM)
        {
          // velocity v_rel_tangential
          static LINALG::Matrix<3,1> v_rel_tangential;
          v_rel_tangential.Update(1.,v_rel,-v_rel_normal,normal);

          // if history variables do not exist -> create them
          if(history_particle_i.find(data_j.id) == history_particle_i.end())
          {
            // safety check
            if(history_particle_j.find(data_i.id) != history_particle_j.end())
              dserror("History of particle j contains particle i, but not vice versa! Something must be wrong...");

            // initialize history of particle i
            PARTICLE::Collision col;
            col.stick = true;
            memset(col.g_t, 0., sizeof(col.g_t));

            // insert new entries
            history_particle_i.insert(std::pair<int,PARTICLE::Collision>(data_j.id,col));
            history_particle_j.insert(std::pair<int,PARTICLE::Collision>(data_i.id,col));
          }

          // safety check
          else if(history_particle_j.find(data_i.id) == history_particle_j.end())
            dserror("History of particle i contains particle j, but not vice versa! Something must be wrong...");

          // extract histories
          PARTICLE::Collision& history_particle_i_j = history_particle_i[data_j.id];
          PARTICLE::Collision& history_particle_j_i = history_particle_j[data_i.id];

          CalculateTangentialContactForce(normalcontactforce, normal, tangentcontactforce,
              history_particle_i_j, v_rel_tangential, m_eff, dt, data_i.owner, data_j.owner);

          // obtain history of particle j by copying from history of particle i
          for(unsigned dim=0; dim<3; ++dim)
            history_particle_j_i.g_t[dim] = -history_particle_i_j.g_t[dim];
          history_particle_j_i.stick = history_particle_i_j.stick;
        }

        // calculation of overall contact forces and moments acting on particles i and j
        static LINALG::Matrix<3,1> contactforce_i, contactforce_j, contactmoment_i, contactmoment_j;
        contactforce_i.Update(normalcontactforce,normal,1.,tangentcontactforce);
        contactforce_j.Update(-1.,contactforce_i); // actio = reactio
        contactmoment_i.CrossProduct(r_iC,contactforce_i);
        contactmoment_j.CrossProduct(r_jC,contactforce_j);

        // assembly contact forces and moments for particles i and j
        LINALG::Assemble(*f_contact, contactforce_i, data_i.lm, data_i.owner);
        LINALG::Assemble(*f_contact, contactforce_j, data_j.lm, data_j.owner);
        LINALG::Assemble(*m_contact, contactmoment_i, data_i.lm, data_i.owner);
        LINALG::Assemble(*m_contact, contactmoment_j, data_j.lm, data_j.owner);

        // --- calculate thermodinamic exchange ---//
        if (specEnthalpyDotn != Teuchos::null)
        {
          // extract the material
          const MAT::PAR::ExtParticleMat* extParticleMat = particle_algorithm_->ExtParticleMat();
          // find the interesting quantities
          static LINALG::Matrix<3,1> r_ij;
          r_ij.Update(1.,data_j.dis,-1.,data_i.dis);
          const double norm_r_contact = r_ij.Norm2();
          const double intersectionArea = PARTICLE::Utils::IntersectionAreaPvsP(data_i.rad,data_j.rad,norm_r_contact);
          const double temperature_i = PARTICLE::Utils::SpecEnthalpy2Temperature(data_i.specEnthalpy,extParticleMat);
          const double temperature_j = PARTICLE::Utils::SpecEnthalpy2Temperature(data_j.specEnthalpy,extParticleMat);

          const double enthalpyDotn2i = extParticleMat->thermalConductivity_ * intersectionArea  * (temperature_j - temperature_i)/norm_r_contact;

          double specEnthalpyDotn2i = enthalpyDotn2i/data_i.mass;
          double specEnthalpyDotn2j = -enthalpyDotn2i/data_j.mass; // actio = - reactio

          LINALG::Assemble(*specEnthalpyDotn, specEnthalpyDotn2i, data_i.id, data_i.owner);
          LINALG::Assemble(*specEnthalpyDotn, specEnthalpyDotn2j, data_j.id, data_j.owner);
        }
      }
      else if(particle_algorithm_->ParticleInteractionType()==INPAR::PARTICLE::NormalAndTang_DEM)// g > 0.0 and no adhesion --> no contact
      {
        // erase entries in histories if still existing
        history_particle_i.erase(data_j.id);
        history_particle_j.erase(data_i.id);
      }
    }
  }  // loop over neighboring particles

  return;
}


/*----------------------------------------------------------------------*
 | calculate contact with neighboring walls                ghamm 04/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerDEM::CalcNeighboringWallsContact(
    DRT::Node*                                        particle_i,          //!< i-th particle
    const ParticleCollData&                           data_i,              //!< collision data associated with i-th particle
    const boost::unordered_map<int,DRT::Element*>&    neighboring_walls,   //!< map of wall elements adjacent to i-th particle
    const double                                      dt,                  //!< time step size
    const Teuchos::RCP<DRT::Discretization>&          walldiscret,         //!< wall discretization
    const Teuchos::RCP<const Epetra_Vector>&          walldisn,            //!< wall displacement
    const Teuchos::RCP<const Epetra_Vector>&          wallveln,            //!< wall velocity
    const Teuchos::RCP<Epetra_Vector>&                f_contact,           //!< global force vector
    const Teuchos::RCP<Epetra_Vector>&                m_contact,           //!< global moment vector
    const Teuchos::RCP<Epetra_FEVector>&              f_structure          //!< global wall force vector
    )
{
  const double radius_i = data_i.rad;
  const double mass_i = data_i.mass;
  const std::vector<int>& lm_i = data_i.lm;
  const int owner_i = data_i.owner;

  // evaluate contact with walls first
  std::vector<WallContactPoint> surfaces;
  std::vector<WallContactPoint> lines;
  std::vector<WallContactPoint> nodes;

  std::set<int> unusedIds;

  // check whether there is contact between particle i and neighboring walls
  for(boost::unordered_map<int, DRT::Element*>::const_iterator w=neighboring_walls.begin(); w!=neighboring_walls.end();  ++w)
  {
    DRT::Element* neighboringwallele = w->second;
    const int numnodes = neighboringwallele->NumNode();
    std::vector<int> lm_wall;
    lm_wall.reserve(numnodes * 3);

    std::vector<int> lmowner;
    std::vector<int> lmstride;
    neighboringwallele->LocationVector(*walldiscret,lm_wall,lmowner,lmstride);

    // nodal displacements
    std::vector<double> nodal_disp(numnodes * 3);
    DRT::UTILS::ExtractMyValues(*walldisn,nodal_disp,lm_wall);

    // get current position of nodes: x = X + u
    std::map<int,LINALG::Matrix<3,1> > nodeCoord;
    DRT::Node** wallnodes = neighboringwallele->Nodes();
    for(int counter=0; counter<numnodes; ++counter)
    {
      static LINALG::Matrix<3,1> currpos;
      const double* X = wallnodes[counter]->X();
      currpos(0) = X[0] + nodal_disp[counter*3+0];
      currpos(1) = X[1] + nodal_disp[counter*3+1];
      currpos(2) = X[2] + nodal_disp[counter*3+2];
      nodeCoord[wallnodes[counter]->Id()] = currpos;
    }

    // compute contact normal, contact gap, contact point on wall element, and type of nearest contact object
    static LINALG::Matrix<3,1> nearestPoint, normal;
    GEO::ObjectType objecttype(GEO::NOTYPE_OBJECT);
    double penetration(0.);
    ComputeContactPointAndNormalAndGapAndObjectType(nearestPoint,objecttype,normal,penetration,neighboringwallele,nodeCoord,data_i);

    if(penetration <= 0.0 or ConsiderNormalAdhesion(penetration))
    {
      // get pointer to the current object type of closest point
      std::vector<WallContactPoint> *pointer=0;
      switch(objecttype)
      {
      case GEO::SURFACE_OBJECT:
      {
        pointer = &surfaces;
      }
      break;
      case GEO::LINE_OBJECT:
      {
        pointer = &lines;
      }
      break;
      case GEO::NODE_OBJECT:
      {
        pointer = &nodes;
      }
      break;
      default:
        dserror("unknown object type");
      break;
      }

      // check, whether point has already been detected (e.g. one line element between two surfaces)
      bool insert = true;
      for(size_t i=0; i<(*pointer).size(); ++i)
      {
        static LINALG::Matrix<3,1> distance_vector;
        distance_vector.Update(1.0, nearestPoint, -1.0, (*pointer)[i].point);
        const double distance = distance_vector.Norm2();
        const double adaptedtol = GEO::TOL7 * radius_i;

        if (distance < adaptedtol)
        {
          // point has already been detected --> do not insert
          insert = false;
          unusedIds.insert(neighboringwallele->Id());
          break;
        }
      }

      // insert contact point with current surface in corresponding map (surf, line, node)
      if(insert)
      {
        WallContactPoint currentContact = { neighboringwallele->Id(), nearestPoint, normal, penetration, nodeCoord, lm_wall, lmowner };
        (*pointer).push_back(currentContact);
      }
    }
    // penetration > 0.0 and no adhesion --> no contact
    else
    {
      unusedIds.insert(neighboringwallele->Id());
    }
  }

  // find entries of lines and nodes which are within the penetration volume of the current particle
  // hierarchical: surfaces first
  for(size_t s=0; s<surfaces.size(); ++s)
  {
    if(surfaces[s].penetration < 0.)
    {
      // within this radius no other contact point can lie: radius = sqrt(r_i^2 - (r_i-|g|)^2)
      const double rminusg = radius_i-std::abs(surfaces[s].penetration);
      const double radius_surface = sqrt(radius_i*radius_i - rminusg*rminusg);

      for(size_t l=0; l<lines.size(); ++l)
      {
        static LINALG::Matrix<3,1> distance_vector;
        distance_vector.Update(1.0, surfaces[s].point, -1.0, lines[l].point);
        const double distance = distance_vector.Norm2();
        if(distance <= radius_surface)
          unusedIds.insert(lines[l].eleid);
      }
      for(size_t p=0; p<nodes.size(); ++p)
      {
        static LINALG::Matrix<3,1> distance_vector;
        distance_vector.Update(1.0, surfaces[s].point, -1.0, nodes[p].point);
        const double distance = distance_vector.Norm2();
        if(distance <= radius_surface)
          unusedIds.insert(nodes[p].eleid);
      }
    }
  }
  // find entries of nodes which are within the penetration volume of the current particle
  // hierarchical: lines next
  for(size_t l=0; l<lines.size(); ++l)
  {
    if(lines[l].penetration < 0.)
    {
      // radius = sqrt(r_i^2 - (r_i-|g|)^2)
      const double rminusg = radius_i-std::abs(lines[l].penetration);
      const double radius_line = sqrt(radius_i*radius_i - rminusg*rminusg);

      for(size_t p=0; p<nodes.size(); ++p)
      {
        static LINALG::Matrix<3,1> distance_vector;
        distance_vector.Update(1.0, lines[l].point, -1.0, nodes[p].point);
        const double distance = distance_vector.Norm2();
        if(distance <= radius_line)
          unusedIds.insert(nodes[p].eleid);
      }
    }
  }

  // write entries of lines and nodes to surfaces if contact has to be evaluated
  for(size_t l=0; l<lines.size(); ++l)
  {
    if( !unusedIds.count(lines[l].eleid) )
      surfaces.push_back(lines[l]);
  }
  for(size_t p=0; p<nodes.size(); ++p)
  {
    if( !unusedIds.count(nodes[p].eleid) )
      surfaces.push_back(nodes[p]);
  }

  // evaluate contact between particle_i and entries of surfaces
  std::map<int, PARTICLE::Collision>& history_wall = static_cast<PARTICLE::ParticleNode*>(particle_i)->Get_history_wall();
  if(history_wall.size() > 3)
    std::cout << "ATTENTION: Contact of particle " << data_i.id << " with " << history_wall.size() << " wall elements." << std::endl;

  for(size_t s=0; s<surfaces.size(); ++s)
  {
    // gid of wall element
    const WallContactPoint& wallcontact = surfaces[s];
    const int gid_wall = wallcontact.eleid;

    //-------get velocity of contact point-----------------------
    static LINALG::Matrix<3,1> vel_nearestPoint;
    vel_nearestPoint.PutScalar(0.0);
    static LINALG::Matrix<2,1> elecoord;
    DRT::Element *CurrentEle = walldiscret->gElement(gid_wall);
    const LINALG::SerialDenseMatrix xyze(GEO::getCurrentNodalPositions(CurrentEle, wallcontact.nodalCoordinates));

    // get coordinates of the projection point in parameter space of the element (xi_coordinates)
    GEO::CurrentToSurfaceElementCoordinates(CurrentEle->Shape(), xyze, wallcontact.point, elecoord);

    const int numnodes = CurrentEle->NumNode();
    Epetra_SerialDenseVector funct(numnodes);

    // get shape functions of the element evaluated at the projection point
    DRT::UTILS::shape_function_2D(funct,elecoord(0,0),elecoord(1,0),CurrentEle->Shape());

    std::vector<double> nodal_vel(numnodes * 3);
    DRT::UTILS::ExtractMyValues(*wallveln,nodal_vel,wallcontact.lm);
    for(int node=0; node<numnodes; ++node)
    {
      for(int dim=0; dim<3; ++dim)
      {
        vel_nearestPoint(dim) += funct[node] * nodal_vel[node * 3 + dim];
      }
    }
    //-----------------------------------------------------------

    // relative velocity in contact point between particle and wall element
    static LINALG::Matrix<3,1> r_iC;
    r_iC.Update(1.,wallcontact.point,-1.,data_i.dis);
    static LINALG::Matrix<3,1> v_rel;
    ComputeRelativeVelocity(v_rel,r_iC,vel_nearestPoint,data_i);

    // part of v_rel in normal-direction: v_rel * n
    const double v_rel_normal(v_rel.Dot(wallcontact.normal));

    // penetration
    // g = norm_r_contact - radius_i;
    const double g(wallcontact.penetration);
    if(g<=0. and std::abs(g)>g_max_)
      g_max_ = std::abs(g);
    //-------------------------------------------------------

    // normalized mass
    const double m_eff = mass_i;

    // contact force
    double normalcontactforce = 0.0;
    static LINALG::Matrix<3,1> tangentcontactforce;
    tangentcontactforce.PutScalar(0.0);

    // normal contact force between particle and wall (note: owner_j = -1)
    CalculateNormalContactForce(g, data_i.rad, 0., v_rel_normal, m_eff, normalcontactforce, owner_i, -1);

    // output normal contact force between particle and wall to *.csv file if applicable
    if(writecontactforcesevery_ and particle_algorithm_->Step() % writecontactforcesevery_ == 0 and owner_i == myrank_)
      OutputNormalContactForceToFile(normalcontactforce,wallcontact.point);

    if(particle_algorithm_->ParticleInteractionType() == INPAR::PARTICLE::NormalAndTang_DEM)
    {
      // velocity v_rel_tangential
      static LINALG::Matrix<3,1> v_rel_tangential;
      v_rel_tangential.Update(1., v_rel, -v_rel_normal, wallcontact.normal);

      // if g < 0 and g_lasttimestep > 0 -> create history variables
      if(history_wall.find(gid_wall) == history_wall.end())
      {
         PARTICLE::Collision col;
         // initialize with stick
         col.stick = true;
         // initialize g_t[3]
         for(int dim=0; dim<3; ++dim)
         {
           col.g_t[dim] = 0.0;
         }

         // insert new entry
         history_wall.insert(std::pair<int,PARTICLE::Collision>(gid_wall,col));
      }

      // calculation of tangential contact force
      CalculateTangentialContactForce(normalcontactforce, wallcontact.normal, tangentcontactforce,
                  history_wall[gid_wall], v_rel_tangential, m_eff, dt, owner_i, -1);
    }

    // calculation of overall contact force
    static LINALG::Matrix<3,1> contactforce;
    contactforce.Update(normalcontactforce,wallcontact.normal,1.,tangentcontactforce);

    // calculation of overall contact moment: m_i = (r_i * n) x F_t
    static LINALG::Matrix<3,1> contactmoment;
    contactmoment.CrossProduct(r_iC,contactforce);

    if (f_contact != Teuchos::null and m_contact != Teuchos::null)
    {
      // do assembly of contact forces
      LINALG::Assemble(*f_contact,contactforce,lm_i,owner_i);
      // do assembly of contact moments
      LINALG::Assemble(*m_contact,contactmoment,lm_i,owner_i);
    }

    if (f_structure != Teuchos::null)
    {
      // forces on wall elements
      double nodal_forces[numnodes * 3];
      for(int node=0; node<numnodes; ++node)
      {
        for(int dim=0; dim<3; ++dim)
        {
          nodal_forces[node * 3 + dim] = funct[node] *(- contactforce(dim));
        }
      }

      // assembly of contact forces on walls
      if(owner_i == myrank_)
      {
        const int err = f_structure->SumIntoGlobalValues(numnodes * 3, &(wallcontact.lm)[0], &nodal_forces[0]);
        if (err<0)
          dserror("summing into Epetra_FEVector failed");
      }
    }
  } // end for contact points on surfaces

  if(particle_algorithm_->ParticleInteractionType() == INPAR::PARTICLE::NormalAndTang_DEM)
  {
    //delete those entries in history_wall_ which are no longer in contact with particle_i in current time step
    for(boost::unordered_map<int, DRT::Element*>::const_iterator w=neighboring_walls.begin(); w != neighboring_walls.end(); ++w)
    {
      const int gid_wall = w->first;
      if( unusedIds.find(gid_wall) != unusedIds.end() and history_wall.find(gid_wall) != history_wall.end() )
        history_wall.erase(gid_wall);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | calculate normal contact force for single contact pair  ghamm 09/13  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerDEM::CalculateNormalContactForce(
  const double g,
  const double radius_i,
  const double radius_j,
  const double v_rel_normal,
  const double m_eff,
  double& normalcontactforce,
  const int owner_i,
  const int owner_j
  )
{
  // evaluate normal contact force arising from penalty contact
  if(g <= 0.)
  {
    // damping parameter
    double d = 0.0;

    //--------------------------------------------------------
    // which contact law: for details see Bachelor thesis Niklas Fehn
    // LinSpring = linear spring
    // Hertz = normal force law of Hertz
    // LinSpringDamp = linear spring damper element
    // LeeHerrmann = nonlinear normal force law of Lee and Herrmann
    // KuwabaraKono = nonlinear normal force law of Kuwabara und Kono
    // Tsuji = nonlinear normal force law of Tsuji
    //---------------------------------------------------------

    // damping parameter
    // TODO: for uni-sized particles this can be done once in the beginning -> efficiency
    // TODO: different m_eff for particle-wall and particel-particle needed when uni-size assumption
    if(owner_j < 0) // contact particle-wall
    {
      if(normal_contact_==INPAR::PARTICLE::LinSpringDamp)
      {
        if(e_wall_ != 0.0)
        {
          const double lnewall = log(e_wall_);
          d = 2.0 * std::abs(lnewall) * sqrt(k_normal_wall_ * m_eff / (lnewall*lnewall + M_PI*M_PI));
        }
        else
        {
          d = 2.0 * sqrt(k_normal_wall_ * m_eff);
        }
      }
      else
      {
        d = d_normal_wall_;
      }
    }
    else // contact particle-particle
    {
      if(normal_contact_==INPAR::PARTICLE::LinSpringDamp)
      {
        if(e_ != 0.0)
        {
          const double lne = log(e_);
          d = 2.0 * std::abs(lne) * sqrt(k_normal_ * m_eff / (lne*lne + M_PI*M_PI));
        }
        else
        {
          d = 2.0 * sqrt(k_normal_ * m_eff);
        }
      }
      else
      {
        d = d_normal_;
      }
    }

    // contact force
    const double k = owner_j < 0 ? k_normal_wall_ : k_normal_;
    switch(normal_contact_)
    {
    case INPAR::PARTICLE::LinSpring:
    {
      normalcontactforce = k * g;

      if(writeenergyevery_)
      {
        //monitor E N E R G Y: here: calculate energy of elastic contact
        contact_energy_ += EnergyAssemble(owner_i,owner_j)* 0.5 * k * g * g;
      }
    }
    break;
    case INPAR::PARTICLE::Hertz:
    {
      normalcontactforce = - k * pow(-g,1.5);

      if(writeenergyevery_)
      {
        //monitor E N E R G Y: here: calculate energy of elastic contact
        contact_energy_ += EnergyAssemble(owner_i,owner_j)* 0.4 * k * pow(-g,2.5);
      }
    }
    break;
    case INPAR::PARTICLE::LinSpringDamp:
    {
      normalcontactforce = k * g - d * v_rel_normal;

      // tension-cutoff
      if(tension_cutoff_ && normalcontactforce > 0.0)
      {
        normalcontactforce = 0.0;
      }

      if(writeenergyevery_)
      {
        //monitor E N E R G Y: here: calculate energy of elastic contact
        contact_energy_ += EnergyAssemble(owner_i,owner_j)* 0.5 * k * g * g;
      }
    }
    break;
    case INPAR::PARTICLE::LeeHerrmann:
    {
      // m_eff = m_i * m_j / ( m_i + m_j)

      normalcontactforce = - k * pow(-g,1.5) - m_eff * d * v_rel_normal;

      // tension-cutoff
      if(tension_cutoff_ && normalcontactforce > 0.0)
      {
        normalcontactforce = 0.0;
      }

      if(writeenergyevery_)
      {
        //monitor E N E R G Y: here: calculate energy of elastic contact
        contact_energy_ += EnergyAssemble(owner_i,owner_j)* 0.4 * k * pow(-g,2.5);
      }
    }
    break;
    case INPAR::PARTICLE::KuwabaraKono:
    {
      normalcontactforce = - k * pow(-g,1.5) - d * v_rel_normal * pow(-g,0.5);

      // tension-cutoff
      if(tension_cutoff_ && normalcontactforce > 0.0)
      {
        normalcontactforce = 0.0;
      }

      if(writeenergyevery_)
      {
        //monitor E N E R G Y: here: calculate energy of elastic contact
        contact_energy_ += EnergyAssemble(owner_i,owner_j)* 0.4 * k * pow(-g,2.5);
      }
    }
    break;
    case INPAR::PARTICLE::Tsuji:
    {
      normalcontactforce = - k * pow(-g,1.5) - d * v_rel_normal * pow(-g,0.25);

      // tension-cutoff
      if(tension_cutoff_ && normalcontactforce > 0.0)
      {
        normalcontactforce = 0.0;
      }

      if(writeenergyevery_)
      {
        //monitor E N E R G Y: here: calculate energy of elastic contact
        contact_energy_ += EnergyAssemble(owner_i,owner_j)* 0.4 * k * pow(-g,2.5);
      }
    }
    break;
    default:
      dserror("specified normal contact law does not exist");
      break;
    }
  }

  // evaluate normal contact force arising from particle adhesion
  if(normal_adhesion_ != INPAR::PARTICLE::adhesion_none)
  {
    // initialize variables for normal adhesion force and energy
    double normaladhesionforce(0.), normaladhesionenergy(0.);

    // calculate normal adhesion force and update internal system energy if applicable
    switch(normal_adhesion_)
    {
      case INPAR::PARTICLE::adhesion_lennardjones:
      {
        const double eq_dist = radius_i+radius_j+adhesion_eq_gap_;
        const double fraction = eq_dist/(radius_i+radius_j+g);
        normaladhesionforce = 12.*adhesion_normal_eps_/eq_dist*(pow(fraction,7)-pow(fraction,13));

        if(writeenergyevery_)
          normaladhesionenergy = EnergyAssemble(owner_i,owner_j)*adhesion_normal_eps_*(pow(fraction,12)-2.*pow(fraction,6));

        break;
      }

      case INPAR::PARTICLE::adhesion_linspring:
      case INPAR::PARTICLE::adhesion_linspringdamp:
      {
        normaladhesionforce = adhesion_normal_stiff_*(g-adhesion_eq_gap_);

        if(normal_adhesion_ == INPAR::PARTICLE::adhesion_linspringdamp)
          normaladhesionforce -= 2.*adhesion_normal_damp_*sqrt(adhesion_normal_stiff_*m_eff)*v_rel_normal;

        if(writeenergyevery_)
          normaladhesionenergy = EnergyAssemble(owner_i,owner_j)*0.5*adhesion_normal_stiff_*pow(g-adhesion_eq_gap_,2);

        break;
      }

      default:
      {
        dserror("Unknown adhesion law governing normal contact of particles!");
        break;
      }
    }

    // add normal adhesion force and energy to total normal contact force and energy if upper force threshold is not exceeded
    if(adhesion_max_force_ < 0. or std::abs(normaladhesionforce) <= adhesion_max_force_)
    {
      normalcontactforce += normaladhesionforce;
      contact_energy_ += normaladhesionenergy;
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | calculate tangential contact force for single           ghamm 09/13  |
 | contact pair                                                         |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerDEM::CalculateTangentialContactForce(
  const double normalcontactforce,
  const LINALG::Matrix<3,1>& normal,
  LINALG::Matrix<3,1>& tangentcontactforce,
  PARTICLE::Collision &currentColl,
  const LINALG::Matrix<3,1>& v_rel_tangential,
  const double m_eff,
  const double dt,
  const int owner_i,
  const int owner_j
  )
{
  // damping parameter
  double d = -0.0;
  // stiffness
  double k = -1.0;
  // frictional coefficient
  double mu = -1.0;

  if(owner_j<0) //contact particle-wall
  {
    // friction
    mu = mu_wall_;
    // stiffness
    k = k_tang_wall_;
    // damping
    if(d_tang_wall_ < 0.0)
    {
      if(e_wall_ != 0.0)
      {
        const double lnewall = log(e_wall_);
        d = 2.0 * std::abs(lnewall) * sqrt(k_normal_wall_ * m_eff / (lnewall*lnewall + M_PI*M_PI));
      }
      else
      {
        d = 2.0 * sqrt(k_normal_wall_ * m_eff);
      }
    }
    else
    {
      d = d_tang_wall_;
    }
  }
  else // contact particle-particle
  {
    // friction
    mu = mu_;
    // stiffness
    k = k_tang_;
    // damping
    if(d_tang_ < 0.0)
    {
      if(e_ != 0.0)
      {
        const double lne = log(e_);
        d = 2.0 * std::abs(lne) * sqrt(k_normal_ * m_eff / (lne*lne + M_PI*M_PI));
      }
      else
      {
        d = 2.0 * sqrt(k_normal_ * m_eff);
      }
    }
    else
    {
      d = d_tang_;
    }
  }

  // store length of g_t at time n
  double old_length = 0.0;
  double interime = 0.0;
  for(int n=0; n<3; ++n)
  {
    old_length += currentColl.g_t[n] * currentColl.g_t[n];
    interime += normal(n) * currentColl.g_t[n];
  }
  old_length = sqrt(old_length);

  // projection of g_t onto current normal at time n+1
  double new_length = 0.0;
  for(int n=0; n<3; ++n)
  {
    currentColl.g_t[n] += - interime * normal(n);
    new_length += currentColl.g_t[n] * currentColl.g_t[n];
  }
  new_length = sqrt(new_length);

  // ensure that g_t has the same length as before projection
  // if almost no tangential spring elongation, neglect it
  if(new_length > 1.0E-14)
  {
    const double scale = old_length/new_length;
    for(int n=0; n<3; ++n)
    {
      currentColl.g_t[n] = scale * currentColl.g_t[n];
    }
  }

  // update of elastic tangential displacement if stick is true
  if(currentColl.stick == true)
  {
    for(int n=0; n<3; ++n)
    {
      currentColl.g_t[n] += v_rel_tangential(n) * dt;
    }
  }

  // calculate tangential test force
  for(int n=0; n<3; ++n)
    tangentcontactforce(n) = - k * currentColl.g_t[n] - d * v_rel_tangential(n);

  // norm of tangential contact force
  const double norm_f_t(tangentcontactforce.Norm2());

  // Coulomb friction law

  // tangential contact force for "stick" - case----------------------
  if( norm_f_t <= (mu * std::abs(normalcontactforce)) )
  {
    currentColl.stick = true;
    //tangential contact force already calculated
  }
  else //"slip"-case
  {
    currentColl.stick = false;
    //calculate tangent vector ( unit vector in (test-)tangentcontactforce-direction )
    static LINALG::Matrix<3,1> tangent;
    tangent.Update(1.0/norm_f_t,tangentcontactforce);

    // calculate tangent contact force and tangential displacements
    tangentcontactforce.Update(mu*std::abs(normalcontactforce),tangent);
    const double kinv = 1.0/k;
    for(int n=0; n<3; ++n)
      currentColl.g_t[n] = - kinv * (tangentcontactforce(n) + d * v_rel_tangential(n));
  }
  //---------------------------------------------------------------

  if(writeenergyevery_)
  {
    new_length = 0.0;
    for(int n=0; n<3; ++n)
    {
      new_length += currentColl.g_t[n] * currentColl.g_t[n];
    }
    new_length = sqrt(new_length);

    contact_energy_ += EnergyAssemble(owner_i,owner_j)* 0.5 * k * new_length * new_length;
  }

  return;
}


/*------------------------------------------------------------------------------------------------------*
 | compute contact normal, contact gap, and vectors from particle centers to contact point   fang 10/17 |
 *------------------------------------------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerDEM::ComputeContactPointAndNormalAndGap(
    LINALG::Matrix<3,1>&      r_iC,     //!< vector from center of particle i to contact point C
    LINALG::Matrix<3,1>&      r_jC,     //!< vector from center of particle j to contact point C
    LINALG::Matrix<3,1>&      normal,   //!< contact normal
    double&                   gap,      //!< contact gap
    const ParticleCollData&   data_i,   //!< collision data for particle i
    const ParticleCollData&   data_j    //!< collision data for particle j
    ) const
{
  // compute distance vector between particles i and j
  normal.Update(1.,data_j.dis,-1.,data_i.dis);

  // compute contact gap, i.e., negative contact penetration
  gap = normal.Norm2()-data_i.rad-data_j.rad;

  // compute contact normal vector
  normal.Scale(1./normal.Norm2());

  // compute vectors from particle centers to contact point
  r_iC.Update(data_i.rad+.5*gap,normal);
  r_jC.Update(-data_j.rad-.5*gap,normal);

  return;
}


/*---------------------------------------------------------------------------------------------------------------------*
 | compute contact normal, contact gap, contact point on wall element, and type of nearest contact object   fang 10/17 |
 *---------------------------------------------------------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerDEM::ComputeContactPointAndNormalAndGapAndObjectType(
    LINALG::Matrix<3,1>&                        C,              //!< contact point on wall element
    GEO::ObjectType&                            objecttype,     //!< type of nearest contact object
    LINALG::Matrix<3,1>&                        normal,         //!< contact normal
    double&                                     gap,            //!< contact gap
    DRT::Element*                               wallele,        //!< wall element
    const std::map<int,LINALG::Matrix<3,1> >&   nodecoord,      //!< node coordinates of wall element
    const ParticleCollData&                     data            //!< particle collision data
    ) const
{
  // find contact point on wall element, i.e., point with smallest distance to particle center
  objecttype = GEO::nearest3DObjectOnElement(wallele,nodecoord,data.dis,C);

  // compute vector from particle center to contact point on wall element
  normal.Update(1.,C,-1.,data.dis);

  // compute distance between particle center and contact point on wall element
  const double distance = normal.Norm2();
  if(distance < GEO::TOL14)
    dserror("Particle center and contact point on wall element coincide! Bad initialization?");

  // compute contact normal vector
  normal.Scale(1./distance);

  // compute contact gap
  gap = distance - data.rad;

  return;
}


/*-------------------------------------------------------------------------------*
 | compute relative velocity in contact point between two particles   fang 10/17 |
 *-------------------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerDEM::ComputeRelativeVelocity(
    LINALG::Matrix<3,1>&         v_rel,    //!< relative velocity between two particles
    const LINALG::Matrix<3,1>&   r_iC,     //!< vector from center of particle i to contact point C
    const LINALG::Matrix<3,1>&   r_jC,     //!< vector from center of particle j to contact point C
    const ParticleCollData&      data_i,   //!< collision data for particle i
    const ParticleCollData&      data_j    //!< collision data for particle j
    ) const
{
  // compute relative velocity between two particles:
  // v_rel = v_i - v_j + omega_i x r_iC - omega_j x r_jC

  // compute translational part
  v_rel.Update(1.,data_i.vel,-1.,data_j.vel);

  // add rotational part, taking account of roundoff in parallel settings
  static LINALG::Matrix<3,1> v_rel_rot_i, v_rel_rot_j;
  v_rel_rot_i.CrossProduct(data_i.angvel,r_iC);
  v_rel_rot_j.CrossProduct(data_j.angvel,r_jC);
  if(data_i.id < data_j.id)
  {
    v_rel += v_rel_rot_i;
    v_rel -= v_rel_rot_j;
  }
  else
  {
    v_rel -= v_rel_rot_j;
    v_rel += v_rel_rot_i;
  }

  return;
}


/*-------------------------------------------------------------------------------------------*
 | compute relative velocity in contact point between particle and wall element   fang 10/17 |
 *-------------------------------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerDEM::ComputeRelativeVelocity(
    LINALG::Matrix<3,1>&         v_rel,   //!< relative velocity between particle and wall element
    const LINALG::Matrix<3,1>&   r_C,     //!< vector from center of particle to contact point C on wall element
    const LINALG::Matrix<3,1>&   v_C,     //!< velocity of contact point C on wall element
    const ParticleCollData&      data     //!< particle collision data
    ) const
{
  // compute relative velocity between particle and wall element:
  // v_rel = v_i + omega_i x r_iC - v_C = (v_i - v_C) + omega_i x r_iC

  // compute translational part
  v_rel.Update(1.,data.vel,-1.,v_C);

  // add rotational part
  static LINALG::Matrix<3,1> v_rel_rot;
  v_rel_rot.CrossProduct(data.angvel,r_C);
  v_rel += v_rel_rot;

  return;
}


/*---------------------------------------------------------------------------*
 | check whether normal particle adhesion needs to be evaluated   fang 02/17 |
 *---------------------------------------------------------------------------*/
bool PARTICLE::ParticleCollisionHandlerDEM::ConsiderNormalAdhesion(const double normalgap) const
{
  return normal_adhesion_ != INPAR::PARTICLE::adhesion_none and (adhesion_max_disp_ < 0. or normalgap <= adhesion_eq_gap_ + adhesion_max_disp_);
}


/*----------------------------------------------------------------------------------------------------------------------*
 | gather *.csv files with normal particle-particle and particle-wall contact forces across all processors   fang 02/17 |
 *----------------------------------------------------------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerDEM::GatherNormalContactForcesToFile() const
{
  // gathering is performed by first processor only
  if(myrank_ == 0)
  {
    // open destination file
    std::ostringstream step;
    step << particle_algorithm_->Step();
    const std::string filename(DRT::Problem::Instance()->OutputControlFile()->FileName()+".particle_contact_forces.csv."+step.str());
    std::ofstream file;
    file.open(filename.c_str());

    // loop over all processors
    for(int iproc=0; iproc<discret_->Comm().NumProc(); ++iproc)
    {
      // open source file associated with current processor
      std::ostringstream iprocstr;
      iprocstr << iproc;
      const std::string procfilename(filename+"."+iprocstr.str());
      std::ifstream procfile;
      procfile.open(procfilename.c_str());

      // gather current source file into destination file
      if(procfile)
        file << procfile.rdbuf();

      // close and remove current source file
      procfile.close();
      if(remove(procfilename.c_str()))
        dserror("Couldn't remove file named "+procfilename+"!");
    }

    // close destination file
    file.close();
  }

  return;
}


/*------------------------------------------------------------------------------------------------------------*
 | write normal contact force between two particles or between particle and wall into *.csv file   fang 02/17 |
 *------------------------------------------------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerDEM::OutputNormalContactForceToFile(
    const double&                normalcontactforce,     //!< normal contact force
    const LINALG::Matrix<3,1>&   forceapplicationpoint   //!< force application point
    ) const
{
  // set file name
  std::ostringstream step;
  step << particle_algorithm_->Step();
  std::ostringstream myrank;
  myrank << myrank_;
  const std::string filename(DRT::Problem::Instance()->OutputControlFile()->FileName()+".particle_contact_forces.csv."+step.str()+"."+myrank.str());

  // open file in appropriate mode
  std::ofstream file;
  FILE* temp = fopen(filename.c_str(),"r");
  if(!temp)
  {
    file.open(filename.c_str(),std::fstream::trunc);

    // write header at beginning of file
    if(myrank_ == 0)
      file << "x,y,z,F" << std::endl;
  }
  else
  {
    fclose(temp);
    file.open(filename.c_str(),std::fstream::app);
  }

  // write application point and value of normal contact force to file
  file << std::setprecision(16) << std::fixed << forceapplicationpoint(0) << "," << forceapplicationpoint(1) << "," << forceapplicationpoint(2) << "," << normalcontactforce << std::endl;

  // close file
  file.close();

  return;
}


/*----------------------------------------------------------------------*
 | constructor                                               fang 09/17 |
 *----------------------------------------------------------------------*/
PARTICLE::ParticleCollisionHandlerDEMEllipsoids::ParticleCollisionHandlerDEMEllipsoids(
    Teuchos::RCP<DRT::Discretization>   discret,             //!< particle discretization
    Teuchos::RCP<PARTICLE::Algorithm>   particlealgorithm,   //!< particle algorithm
    const Teuchos::ParameterList&       particledynparams    //!< parameter list
    ) :
    // call base class constructor
    ParticleCollisionHandlerDEM(discret,particlealgorithm,particledynparams)
{
  // safety check
  if(normal_adhesion_ != INPAR::PARTICLE::adhesion_none)
    dserror("Adhesion not yet implemented for ellipsoidal particles!");

  // store minimum and maximum radii
  particle_algorithm_->AdapterParticle()->Radiusn()->MinValue(&r_min_);
  particle_algorithm_->AdapterParticle()->Radiusn()->MaxValue(&r_max_);

  // the following code solely serves to estimate the penalty contact stiffness
  // to be explicitly prescribed in the *.dat file prior to the actual simulation

  // extract material parameters for ellipsoidal particles
  const MAT::PAR::ParticleMatEllipsoids* const particlematellipsoids = dynamic_cast<const MAT::PAR::ParticleMatEllipsoids* const>(particle_algorithm_->ParticleMat());
  if(particlematellipsoids == NULL)
    dserror("Couldn't extract material parameters for ellipsoidal particles!");

  // extract semi-axes of ellipsoidal particles
  const LINALG::Matrix<3,1>& semiaxes = particlematellipsoids->SemiAxes();

  // safety checks
  switch(normal_contact_)
  {
    case INPAR::PARTICLE::LinSpring:
    case INPAR::PARTICLE::LinSpringDamp:
    {
      // stiffness value not explicitly prescribed
      if(particledynparams.get<double>("NORMAL_STIFF") <= 0. or
        (particle_algorithm_->WallDiscret() != Teuchos::null and particledynparams.get<double>("NORMAL_STIFF_WALL") <= 0.))
      {
        // extract maximum angular velocity of ellipsoidal particles
        const double angVel_max = particledynparams.get<double>("MAX_ANGULAR_VELOCITY");

        // compute equivalent radius
        const double r_equiv = pow(semiaxes(0)*semiaxes(1)*semiaxes(2),1./3.);

        // compute square of maximum rotational velocity
        const double v_rot_sq = .2*angVel_max*angVel_max*(r_max_*r_max_+r_equiv*r_equiv/(r_max_*r_max_*r_min_*r_min_));

        // safety check for particle-particle contact
        if(particledynparams.get<double>("NORMAL_STIFF") <= 0.)
        {
          // minimum stiffness value in conservative case
          const double k_normal_conservative = 2./3.*r_equiv*r_equiv*r_equiv*M_PI*particlematellipsoids->initDensity_*(v_max_*v_max_+v_rot_sq)/(c_*c_*r_max_*r_max_);

          // stiffness value must be explicitly prescribed
          dserror("For ellipsoidal particles, the penalty contact stiffness for particle-particle contact "
                  "must be prescribed manually! The value should be at least greater than %lf (no rotation considered) "
                  "or even greater than %lf (worst-case scenario with rotation, very conservative).",k_normal_,k_normal_conservative);
        }

        // safety check for particle-wall contact
        else
        {
          // minimum stiffness value in conservative case
          const double k_normal_wall_conservative = 2./3.*r_equiv*r_equiv*r_equiv*M_PI*particlematellipsoids->initDensity_*(v_max_*v_max_+v_rot_sq)/(c_wall_*c_wall_*r_max_*r_max_);

          // stiffness value must be explicitly prescribed
          dserror("For ellipsoidal particles, the penalty contact stiffness for particle-wall contact "
                  "must be prescribed manually! The value should be at least greater than %lf (no rotation considered) "
                  "or even greater than %lf (worst-case scenario with rotation, very conservative).",k_normal_wall_,k_normal_wall_conservative);
        }
      }

      break;
    }

    case INPAR::PARTICLE::Hertz:
    case INPAR::PARTICLE::LeeHerrmann:
    case INPAR::PARTICLE::KuwabaraKono:
    case INPAR::PARTICLE::Tsuji:
    {
      dserror("Chosen normal contact model not available for ellipsoidal particles!");
      break;
    }

    default:
    {
      dserror("Invalid particle contact interaction type!");
      break;
    }
  }

  return;
}


/*----------------------------------------------------------------------------*
 | check whether two ellipsoids collide                            fang 10/17 |
 *----------------------------------------------------------------------------*/
bool PARTICLE::ParticleCollisionHandlerDEMEllipsoids::CollisionCheck(
    const LINALG::Matrix<4,4>&   E_i,   //!< first ellipsoid matrix in homogeneous coordinates
    const LINALG::Matrix<4,4>&   E_j    //!< second ellipsoid matrix in homogeneous coordinates
    ) const
{
  // We consider two ellipsoids described by the ellipsoid matrices E_i and E_j, i.e., by quadric equations in
  // homogeneous coordinates. The generalized eigenproblem det(lambda E_i + E_j) = 0 leads to the characteristic
  // equation a*lambda^4 + b*lambda^3 + c*lambda^2 + d*lambda + e. The two ellipsoids are separated if and
  // only if the characteristic equation has two distinct positive roots. One can show that this is equivalent
  // to all four roots being real.

  // change sign of matrix E_i to obtain standard generalized eigenproblem E_j*x = lambda E_i_temp*x
  static LINALG::Matrix<4,4> E_i_temp;
  E_i_temp = E_i;
  E_i_temp.Scale(-1.);

  // check whether standard generalized eigenproblem exhibits complex eigenvalues
  static LINALG::Matrix<3,1> dummy;
  return ComputeEigenVectors(E_j,E_i_temp,dummy,dummy,dummy,dummy);

  // alternative way of checking whether standard generalized eigenproblem exhibits complex eigenvalues:
  // only compute eigenvalues, but not the corresponding eigenvectors, employing a different set of LAPACK routines
  // however: method not reliable, sometimes gives small, artificial complex parts of actually purely real eigenvalues
  //      --> method deactivated at the moment...
  // return ComplexEigenValues(E_j,E_i_temp);

  // alternative way of checking whether the characteristic equation has two distinct positive roots:
  // discriminant of characteristic equation must be positive and auxiliary quantity P (see below) must be negative
  // however: method not reliable in case of small discriminant (absolute value < 1.e-7), most probably due to roundoff
  //      --> method deactivated at the moment...
  /* START OF ALTERNATIVE METHOD
  // compute coefficients of characteristic equation
  double a = E_i(0,0)*E_i(1,1)*E_i(2,2)*E_i(3,3) - E_i(0,0)*E_i(1,1)*E_i(2,3)*E_i(3,2) - E_i(0,0)*E_i(1,2)*E_i(2,1)*E_i(3,3)
           + E_i(0,0)*E_i(1,2)*E_i(2,3)*E_i(3,1) + E_i(0,0)*E_i(1,3)*E_i(2,1)*E_i(3,2) - E_i(0,0)*E_i(1,3)*E_i(2,2)*E_i(3,1)
           - E_i(0,1)*E_i(1,0)*E_i(2,2)*E_i(3,3) + E_i(0,1)*E_i(1,0)*E_i(2,3)*E_i(3,2) + E_i(0,1)*E_i(1,2)*E_i(2,0)*E_i(3,3)
           - E_i(0,1)*E_i(1,2)*E_i(2,3)*E_i(3,0) - E_i(0,1)*E_i(1,3)*E_i(2,0)*E_i(3,2) + E_i(0,1)*E_i(1,3)*E_i(2,2)*E_i(3,0)
           + E_i(0,2)*E_i(1,0)*E_i(2,1)*E_i(3,3) - E_i(0,2)*E_i(1,0)*E_i(2,3)*E_i(3,1) - E_i(0,2)*E_i(1,1)*E_i(2,0)*E_i(3,3)
           + E_i(0,2)*E_i(1,1)*E_i(2,3)*E_i(3,0) + E_i(0,2)*E_i(1,3)*E_i(2,0)*E_i(3,1) - E_i(0,2)*E_i(1,3)*E_i(2,1)*E_i(3,0)
           - E_i(0,3)*E_i(1,0)*E_i(2,1)*E_i(3,2) + E_i(0,3)*E_i(1,0)*E_i(2,2)*E_i(3,1) + E_i(0,3)*E_i(1,1)*E_i(2,0)*E_i(3,2)
           - E_i(0,3)*E_i(1,1)*E_i(2,2)*E_i(3,0) - E_i(0,3)*E_i(1,2)*E_i(2,0)*E_i(3,1) + E_i(0,3)*E_i(1,2)*E_i(2,1)*E_i(3,0);
  double b = E_i(0,0)*E_i(1,1)*E_i(2,2)*E_j(3,3) - E_i(0,0)*E_i(1,1)*E_i(2,3)*E_j(3,2) - E_i(0,0)*E_i(1,1)*E_i(3,2)*E_j(2,3)
           + E_i(0,0)*E_i(1,1)*E_i(3,3)*E_j(2,2) - E_i(0,0)*E_i(1,2)*E_i(2,1)*E_j(3,3) + E_i(0,0)*E_i(1,2)*E_i(2,3)*E_j(3,1)
           + E_i(0,0)*E_i(1,2)*E_i(3,1)*E_j(2,3) - E_i(0,0)*E_i(1,2)*E_i(3,3)*E_j(2,1) + E_i(0,0)*E_i(1,3)*E_i(2,1)*E_j(3,2)
           - E_i(0,0)*E_i(1,3)*E_i(2,2)*E_j(3,1) - E_i(0,0)*E_i(1,3)*E_i(3,1)*E_j(2,2) + E_i(0,0)*E_i(1,3)*E_i(3,2)*E_j(2,1)
           + E_i(0,0)*E_i(2,1)*E_i(3,2)*E_j(1,3) - E_i(0,0)*E_i(2,1)*E_i(3,3)*E_j(1,2) - E_i(0,0)*E_i(2,2)*E_i(3,1)*E_j(1,3)
           + E_i(0,0)*E_i(2,2)*E_i(3,3)*E_j(1,1) + E_i(0,0)*E_i(2,3)*E_i(3,1)*E_j(1,2) - E_i(0,0)*E_i(2,3)*E_i(3,2)*E_j(1,1)
           - E_i(0,1)*E_i(1,0)*E_i(2,2)*E_j(3,3) + E_i(0,1)*E_i(1,0)*E_i(2,3)*E_j(3,2) + E_i(0,1)*E_i(1,0)*E_i(3,2)*E_j(2,3)
           - E_i(0,1)*E_i(1,0)*E_i(3,3)*E_j(2,2) + E_i(0,1)*E_i(1,2)*E_i(2,0)*E_j(3,3) - E_i(0,1)*E_i(1,2)*E_i(2,3)*E_j(3,0)
           - E_i(0,1)*E_i(1,2)*E_i(3,0)*E_j(2,3) + E_i(0,1)*E_i(1,2)*E_i(3,3)*E_j(2,0) - E_i(0,1)*E_i(1,3)*E_i(2,0)*E_j(3,2)
           + E_i(0,1)*E_i(1,3)*E_i(2,2)*E_j(3,0) + E_i(0,1)*E_i(1,3)*E_i(3,0)*E_j(2,2) - E_i(0,1)*E_i(1,3)*E_i(3,2)*E_j(2,0)
           - E_i(0,1)*E_i(2,0)*E_i(3,2)*E_j(1,3) + E_i(0,1)*E_i(2,0)*E_i(3,3)*E_j(1,2) + E_i(0,1)*E_i(2,2)*E_i(3,0)*E_j(1,3)
           - E_i(0,1)*E_i(2,2)*E_i(3,3)*E_j(1,0) - E_i(0,1)*E_i(2,3)*E_i(3,0)*E_j(1,2) + E_i(0,1)*E_i(2,3)*E_i(3,2)*E_j(1,0)
           + E_i(0,2)*E_i(1,0)*E_i(2,1)*E_j(3,3) - E_i(0,2)*E_i(1,0)*E_i(2,3)*E_j(3,1) - E_i(0,2)*E_i(1,0)*E_i(3,1)*E_j(2,3)
           + E_i(0,2)*E_i(1,0)*E_i(3,3)*E_j(2,1) - E_i(0,2)*E_i(1,1)*E_i(2,0)*E_j(3,3) + E_i(0,2)*E_i(1,1)*E_i(2,3)*E_j(3,0)
           + E_i(0,2)*E_i(1,1)*E_i(3,0)*E_j(2,3) - E_i(0,2)*E_i(1,1)*E_i(3,3)*E_j(2,0) + E_i(0,2)*E_i(1,3)*E_i(2,0)*E_j(3,1)
           - E_i(0,2)*E_i(1,3)*E_i(2,1)*E_j(3,0) - E_i(0,2)*E_i(1,3)*E_i(3,0)*E_j(2,1) + E_i(0,2)*E_i(1,3)*E_i(3,1)*E_j(2,0)
           + E_i(0,2)*E_i(2,0)*E_i(3,1)*E_j(1,3) - E_i(0,2)*E_i(2,0)*E_i(3,3)*E_j(1,1) - E_i(0,2)*E_i(2,1)*E_i(3,0)*E_j(1,3)
           + E_i(0,2)*E_i(2,1)*E_i(3,3)*E_j(1,0) + E_i(0,2)*E_i(2,3)*E_i(3,0)*E_j(1,1) - E_i(0,2)*E_i(2,3)*E_i(3,1)*E_j(1,0)
           - E_i(0,3)*E_i(1,0)*E_i(2,1)*E_j(3,2) + E_i(0,3)*E_i(1,0)*E_i(2,2)*E_j(3,1) + E_i(0,3)*E_i(1,0)*E_i(3,1)*E_j(2,2)
           - E_i(0,3)*E_i(1,0)*E_i(3,2)*E_j(2,1) + E_i(0,3)*E_i(1,1)*E_i(2,0)*E_j(3,2) - E_i(0,3)*E_i(1,1)*E_i(2,2)*E_j(3,0)
           - E_i(0,3)*E_i(1,1)*E_i(3,0)*E_j(2,2) + E_i(0,3)*E_i(1,1)*E_i(3,2)*E_j(2,0) - E_i(0,3)*E_i(1,2)*E_i(2,0)*E_j(3,1)
           + E_i(0,3)*E_i(1,2)*E_i(2,1)*E_j(3,0) + E_i(0,3)*E_i(1,2)*E_i(3,0)*E_j(2,1) - E_i(0,3)*E_i(1,2)*E_i(3,1)*E_j(2,0)
           - E_i(0,3)*E_i(2,0)*E_i(3,1)*E_j(1,2) + E_i(0,3)*E_i(2,0)*E_i(3,2)*E_j(1,1) + E_i(0,3)*E_i(2,1)*E_i(3,0)*E_j(1,2)
           - E_i(0,3)*E_i(2,1)*E_i(3,2)*E_j(1,0) - E_i(0,3)*E_i(2,2)*E_i(3,0)*E_j(1,1) + E_i(0,3)*E_i(2,2)*E_i(3,1)*E_j(1,0)
           - E_i(1,0)*E_i(2,1)*E_i(3,2)*E_j(0,3) + E_i(1,0)*E_i(2,1)*E_i(3,3)*E_j(0,2) + E_i(1,0)*E_i(2,2)*E_i(3,1)*E_j(0,3)
           - E_i(1,0)*E_i(2,2)*E_i(3,3)*E_j(0,1) - E_i(1,0)*E_i(2,3)*E_i(3,1)*E_j(0,2) + E_i(1,0)*E_i(2,3)*E_i(3,2)*E_j(0,1)
           + E_i(1,1)*E_i(2,0)*E_i(3,2)*E_j(0,3) - E_i(1,1)*E_i(2,0)*E_i(3,3)*E_j(0,2) - E_i(1,1)*E_i(2,2)*E_i(3,0)*E_j(0,3)
           + E_i(1,1)*E_i(2,2)*E_i(3,3)*E_j(0,0) + E_i(1,1)*E_i(2,3)*E_i(3,0)*E_j(0,2) - E_i(1,1)*E_i(2,3)*E_i(3,2)*E_j(0,0)
           - E_i(1,2)*E_i(2,0)*E_i(3,1)*E_j(0,3) + E_i(1,2)*E_i(2,0)*E_i(3,3)*E_j(0,1) + E_i(1,2)*E_i(2,1)*E_i(3,0)*E_j(0,3)
           - E_i(1,2)*E_i(2,1)*E_i(3,3)*E_j(0,0) - E_i(1,2)*E_i(2,3)*E_i(3,0)*E_j(0,1) + E_i(1,2)*E_i(2,3)*E_i(3,1)*E_j(0,0)
           + E_i(1,3)*E_i(2,0)*E_i(3,1)*E_j(0,2) - E_i(1,3)*E_i(2,0)*E_i(3,2)*E_j(0,1) - E_i(1,3)*E_i(2,1)*E_i(3,0)*E_j(0,2)
           + E_i(1,3)*E_i(2,1)*E_i(3,2)*E_j(0,0) + E_i(1,3)*E_i(2,2)*E_i(3,0)*E_j(0,1) - E_i(1,3)*E_i(2,2)*E_i(3,1)*E_j(0,0);
  double c = E_i(0,0)*E_i(1,1)*E_j(2,2)*E_j(3,3) - E_i(0,0)*E_i(1,1)*E_j(2,3)*E_j(3,2) - E_i(0,0)*E_i(1,2)*E_j(2,1)*E_j(3,3)
           + E_i(0,0)*E_i(1,2)*E_j(2,3)*E_j(3,1) + E_i(0,0)*E_i(1,3)*E_j(2,1)*E_j(3,2) - E_i(0,0)*E_i(1,3)*E_j(2,2)*E_j(3,1)
           - E_i(0,0)*E_i(2,1)*E_j(1,2)*E_j(3,3) + E_i(0,0)*E_i(2,1)*E_j(1,3)*E_j(3,2) + E_i(0,0)*E_i(2,2)*E_j(1,1)*E_j(3,3)
           - E_i(0,0)*E_i(2,2)*E_j(1,3)*E_j(3,1) - E_i(0,0)*E_i(2,3)*E_j(1,1)*E_j(3,2) + E_i(0,0)*E_i(2,3)*E_j(1,2)*E_j(3,1)
           + E_i(0,0)*E_i(3,1)*E_j(1,2)*E_j(2,3) - E_i(0,0)*E_i(3,1)*E_j(1,3)*E_j(2,2) - E_i(0,0)*E_i(3,2)*E_j(1,1)*E_j(2,3)
           + E_i(0,0)*E_i(3,2)*E_j(1,3)*E_j(2,1) + E_i(0,0)*E_i(3,3)*E_j(1,1)*E_j(2,2) - E_i(0,0)*E_i(3,3)*E_j(1,2)*E_j(2,1)
           - E_i(0,1)*E_i(1,0)*E_j(2,2)*E_j(3,3) + E_i(0,1)*E_i(1,0)*E_j(2,3)*E_j(3,2) + E_i(0,1)*E_i(1,2)*E_j(2,0)*E_j(3,3)
           - E_i(0,1)*E_i(1,2)*E_j(2,3)*E_j(3,0) - E_i(0,1)*E_i(1,3)*E_j(2,0)*E_j(3,2) + E_i(0,1)*E_i(1,3)*E_j(2,2)*E_j(3,0)
           + E_i(0,1)*E_i(2,0)*E_j(1,2)*E_j(3,3) - E_i(0,1)*E_i(2,0)*E_j(1,3)*E_j(3,2) - E_i(0,1)*E_i(2,2)*E_j(1,0)*E_j(3,3)
           + E_i(0,1)*E_i(2,2)*E_j(1,3)*E_j(3,0) + E_i(0,1)*E_i(2,3)*E_j(1,0)*E_j(3,2) - E_i(0,1)*E_i(2,3)*E_j(1,2)*E_j(3,0)
           - E_i(0,1)*E_i(3,0)*E_j(1,2)*E_j(2,3) + E_i(0,1)*E_i(3,0)*E_j(1,3)*E_j(2,2) + E_i(0,1)*E_i(3,2)*E_j(1,0)*E_j(2,3)
           - E_i(0,1)*E_i(3,2)*E_j(1,3)*E_j(2,0) - E_i(0,1)*E_i(3,3)*E_j(1,0)*E_j(2,2) + E_i(0,1)*E_i(3,3)*E_j(1,2)*E_j(2,0)
           + E_i(0,2)*E_i(1,0)*E_j(2,1)*E_j(3,3) - E_i(0,2)*E_i(1,0)*E_j(2,3)*E_j(3,1) - E_i(0,2)*E_i(1,1)*E_j(2,0)*E_j(3,3)
           + E_i(0,2)*E_i(1,1)*E_j(2,3)*E_j(3,0) + E_i(0,2)*E_i(1,3)*E_j(2,0)*E_j(3,1) - E_i(0,2)*E_i(1,3)*E_j(2,1)*E_j(3,0)
           - E_i(0,2)*E_i(2,0)*E_j(1,1)*E_j(3,3) + E_i(0,2)*E_i(2,0)*E_j(1,3)*E_j(3,1) + E_i(0,2)*E_i(2,1)*E_j(1,0)*E_j(3,3)
           - E_i(0,2)*E_i(2,1)*E_j(1,3)*E_j(3,0) - E_i(0,2)*E_i(2,3)*E_j(1,0)*E_j(3,1) + E_i(0,2)*E_i(2,3)*E_j(1,1)*E_j(3,0)
           + E_i(0,2)*E_i(3,0)*E_j(1,1)*E_j(2,3) - E_i(0,2)*E_i(3,0)*E_j(1,3)*E_j(2,1) - E_i(0,2)*E_i(3,1)*E_j(1,0)*E_j(2,3)
           + E_i(0,2)*E_i(3,1)*E_j(1,3)*E_j(2,0) + E_i(0,2)*E_i(3,3)*E_j(1,0)*E_j(2,1) - E_i(0,2)*E_i(3,3)*E_j(1,1)*E_j(2,0)
           - E_i(0,3)*E_i(1,0)*E_j(2,1)*E_j(3,2) + E_i(0,3)*E_i(1,0)*E_j(2,2)*E_j(3,1) + E_i(0,3)*E_i(1,1)*E_j(2,0)*E_j(3,2)
           - E_i(0,3)*E_i(1,1)*E_j(2,2)*E_j(3,0) - E_i(0,3)*E_i(1,2)*E_j(2,0)*E_j(3,1) + E_i(0,3)*E_i(1,2)*E_j(2,1)*E_j(3,0)
           + E_i(0,3)*E_i(2,0)*E_j(1,1)*E_j(3,2) - E_i(0,3)*E_i(2,0)*E_j(1,2)*E_j(3,1) - E_i(0,3)*E_i(2,1)*E_j(1,0)*E_j(3,2)
           + E_i(0,3)*E_i(2,1)*E_j(1,2)*E_j(3,0) + E_i(0,3)*E_i(2,2)*E_j(1,0)*E_j(3,1) - E_i(0,3)*E_i(2,2)*E_j(1,1)*E_j(3,0)
           - E_i(0,3)*E_i(3,0)*E_j(1,1)*E_j(2,2) + E_i(0,3)*E_i(3,0)*E_j(1,2)*E_j(2,1) + E_i(0,3)*E_i(3,1)*E_j(1,0)*E_j(2,2)
           - E_i(0,3)*E_i(3,1)*E_j(1,2)*E_j(2,0) - E_i(0,3)*E_i(3,2)*E_j(1,0)*E_j(2,1) + E_i(0,3)*E_i(3,2)*E_j(1,1)*E_j(2,0)
           + E_i(1,0)*E_i(2,1)*E_j(0,2)*E_j(3,3) - E_i(1,0)*E_i(2,1)*E_j(0,3)*E_j(3,2) - E_i(1,0)*E_i(2,2)*E_j(0,1)*E_j(3,3)
           + E_i(1,0)*E_i(2,2)*E_j(0,3)*E_j(3,1) + E_i(1,0)*E_i(2,3)*E_j(0,1)*E_j(3,2) - E_i(1,0)*E_i(2,3)*E_j(0,2)*E_j(3,1)
           - E_i(1,0)*E_i(3,1)*E_j(0,2)*E_j(2,3) + E_i(1,0)*E_i(3,1)*E_j(0,3)*E_j(2,2) + E_i(1,0)*E_i(3,2)*E_j(0,1)*E_j(2,3)
           - E_i(1,0)*E_i(3,2)*E_j(0,3)*E_j(2,1) - E_i(1,0)*E_i(3,3)*E_j(0,1)*E_j(2,2) + E_i(1,0)*E_i(3,3)*E_j(0,2)*E_j(2,1)
           - E_i(1,1)*E_i(2,0)*E_j(0,2)*E_j(3,3) + E_i(1,1)*E_i(2,0)*E_j(0,3)*E_j(3,2) + E_i(1,1)*E_i(2,2)*E_j(0,0)*E_j(3,3)
           - E_i(1,1)*E_i(2,2)*E_j(0,3)*E_j(3,0) - E_i(1,1)*E_i(2,3)*E_j(0,0)*E_j(3,2) + E_i(1,1)*E_i(2,3)*E_j(0,2)*E_j(3,0)
           + E_i(1,1)*E_i(3,0)*E_j(0,2)*E_j(2,3) - E_i(1,1)*E_i(3,0)*E_j(0,3)*E_j(2,2) - E_i(1,1)*E_i(3,2)*E_j(0,0)*E_j(2,3)
           + E_i(1,1)*E_i(3,2)*E_j(0,3)*E_j(2,0) + E_i(1,1)*E_i(3,3)*E_j(0,0)*E_j(2,2) - E_i(1,1)*E_i(3,3)*E_j(0,2)*E_j(2,0)
           + E_i(1,2)*E_i(2,0)*E_j(0,1)*E_j(3,3) - E_i(1,2)*E_i(2,0)*E_j(0,3)*E_j(3,1) - E_i(1,2)*E_i(2,1)*E_j(0,0)*E_j(3,3)
           + E_i(1,2)*E_i(2,1)*E_j(0,3)*E_j(3,0) + E_i(1,2)*E_i(2,3)*E_j(0,0)*E_j(3,1) - E_i(1,2)*E_i(2,3)*E_j(0,1)*E_j(3,0)
           - E_i(1,2)*E_i(3,0)*E_j(0,1)*E_j(2,3) + E_i(1,2)*E_i(3,0)*E_j(0,3)*E_j(2,1) + E_i(1,2)*E_i(3,1)*E_j(0,0)*E_j(2,3)
           - E_i(1,2)*E_i(3,1)*E_j(0,3)*E_j(2,0) - E_i(1,2)*E_i(3,3)*E_j(0,0)*E_j(2,1) + E_i(1,2)*E_i(3,3)*E_j(0,1)*E_j(2,0)
           - E_i(1,3)*E_i(2,0)*E_j(0,1)*E_j(3,2) + E_i(1,3)*E_i(2,0)*E_j(0,2)*E_j(3,1) + E_i(1,3)*E_i(2,1)*E_j(0,0)*E_j(3,2)
           - E_i(1,3)*E_i(2,1)*E_j(0,2)*E_j(3,0) - E_i(1,3)*E_i(2,2)*E_j(0,0)*E_j(3,1) + E_i(1,3)*E_i(2,2)*E_j(0,1)*E_j(3,0)
           + E_i(1,3)*E_i(3,0)*E_j(0,1)*E_j(2,2) - E_i(1,3)*E_i(3,0)*E_j(0,2)*E_j(2,1) - E_i(1,3)*E_i(3,1)*E_j(0,0)*E_j(2,2)
           + E_i(1,3)*E_i(3,1)*E_j(0,2)*E_j(2,0) + E_i(1,3)*E_i(3,2)*E_j(0,0)*E_j(2,1) - E_i(1,3)*E_i(3,2)*E_j(0,1)*E_j(2,0)
           + E_i(2,0)*E_i(3,1)*E_j(0,2)*E_j(1,3) - E_i(2,0)*E_i(3,1)*E_j(0,3)*E_j(1,2) - E_i(2,0)*E_i(3,2)*E_j(0,1)*E_j(1,3)
           + E_i(2,0)*E_i(3,2)*E_j(0,3)*E_j(1,1) + E_i(2,0)*E_i(3,3)*E_j(0,1)*E_j(1,2) - E_i(2,0)*E_i(3,3)*E_j(0,2)*E_j(1,1)
           - E_i(2,1)*E_i(3,0)*E_j(0,2)*E_j(1,3) + E_i(2,1)*E_i(3,0)*E_j(0,3)*E_j(1,2) + E_i(2,1)*E_i(3,2)*E_j(0,0)*E_j(1,3)
           - E_i(2,1)*E_i(3,2)*E_j(0,3)*E_j(1,0) - E_i(2,1)*E_i(3,3)*E_j(0,0)*E_j(1,2) + E_i(2,1)*E_i(3,3)*E_j(0,2)*E_j(1,0)
           + E_i(2,2)*E_i(3,0)*E_j(0,1)*E_j(1,3) - E_i(2,2)*E_i(3,0)*E_j(0,3)*E_j(1,1) - E_i(2,2)*E_i(3,1)*E_j(0,0)*E_j(1,3)
           + E_i(2,2)*E_i(3,1)*E_j(0,3)*E_j(1,0) + E_i(2,2)*E_i(3,3)*E_j(0,0)*E_j(1,1) - E_i(2,2)*E_i(3,3)*E_j(0,1)*E_j(1,0)
           - E_i(2,3)*E_i(3,0)*E_j(0,1)*E_j(1,2) + E_i(2,3)*E_i(3,0)*E_j(0,2)*E_j(1,1) + E_i(2,3)*E_i(3,1)*E_j(0,0)*E_j(1,2)
           - E_i(2,3)*E_i(3,1)*E_j(0,2)*E_j(1,0) - E_i(2,3)*E_i(3,2)*E_j(0,0)*E_j(1,1) + E_i(2,3)*E_i(3,2)*E_j(0,1)*E_j(1,0);
  double d = E_i(0,0)*E_j(1,1)*E_j(2,2)*E_j(3,3) - E_i(0,0)*E_j(1,1)*E_j(2,3)*E_j(3,2) - E_i(0,0)*E_j(1,2)*E_j(2,1)*E_j(3,3)
           + E_i(0,0)*E_j(1,2)*E_j(2,3)*E_j(3,1) + E_i(0,0)*E_j(1,3)*E_j(2,1)*E_j(3,2) - E_i(0,0)*E_j(1,3)*E_j(2,2)*E_j(3,1)
           - E_i(0,1)*E_j(1,0)*E_j(2,2)*E_j(3,3) + E_i(0,1)*E_j(1,0)*E_j(2,3)*E_j(3,2) + E_i(0,1)*E_j(1,2)*E_j(2,0)*E_j(3,3)
           - E_i(0,1)*E_j(1,2)*E_j(2,3)*E_j(3,0) - E_i(0,1)*E_j(1,3)*E_j(2,0)*E_j(3,2) + E_i(0,1)*E_j(1,3)*E_j(2,2)*E_j(3,0)
           + E_i(0,2)*E_j(1,0)*E_j(2,1)*E_j(3,3) - E_i(0,2)*E_j(1,0)*E_j(2,3)*E_j(3,1) - E_i(0,2)*E_j(1,1)*E_j(2,0)*E_j(3,3)
           + E_i(0,2)*E_j(1,1)*E_j(2,3)*E_j(3,0) + E_i(0,2)*E_j(1,3)*E_j(2,0)*E_j(3,1) - E_i(0,2)*E_j(1,3)*E_j(2,1)*E_j(3,0)
           - E_i(0,3)*E_j(1,0)*E_j(2,1)*E_j(3,2) + E_i(0,3)*E_j(1,0)*E_j(2,2)*E_j(3,1) + E_i(0,3)*E_j(1,1)*E_j(2,0)*E_j(3,2)
           - E_i(0,3)*E_j(1,1)*E_j(2,2)*E_j(3,0) - E_i(0,3)*E_j(1,2)*E_j(2,0)*E_j(3,1) + E_i(0,3)*E_j(1,2)*E_j(2,1)*E_j(3,0)
           - E_i(1,0)*E_j(0,1)*E_j(2,2)*E_j(3,3) + E_i(1,0)*E_j(0,1)*E_j(2,3)*E_j(3,2) + E_i(1,0)*E_j(0,2)*E_j(2,1)*E_j(3,3)
           - E_i(1,0)*E_j(0,2)*E_j(2,3)*E_j(3,1) - E_i(1,0)*E_j(0,3)*E_j(2,1)*E_j(3,2) + E_i(1,0)*E_j(0,3)*E_j(2,2)*E_j(3,1)
           + E_i(1,1)*E_j(0,0)*E_j(2,2)*E_j(3,3) - E_i(1,1)*E_j(0,0)*E_j(2,3)*E_j(3,2) - E_i(1,1)*E_j(0,2)*E_j(2,0)*E_j(3,3)
           + E_i(1,1)*E_j(0,2)*E_j(2,3)*E_j(3,0) + E_i(1,1)*E_j(0,3)*E_j(2,0)*E_j(3,2) - E_i(1,1)*E_j(0,3)*E_j(2,2)*E_j(3,0)
           - E_i(1,2)*E_j(0,0)*E_j(2,1)*E_j(3,3) + E_i(1,2)*E_j(0,0)*E_j(2,3)*E_j(3,1) + E_i(1,2)*E_j(0,1)*E_j(2,0)*E_j(3,3)
           - E_i(1,2)*E_j(0,1)*E_j(2,3)*E_j(3,0) - E_i(1,2)*E_j(0,3)*E_j(2,0)*E_j(3,1) + E_i(1,2)*E_j(0,3)*E_j(2,1)*E_j(3,0)
           + E_i(1,3)*E_j(0,0)*E_j(2,1)*E_j(3,2) - E_i(1,3)*E_j(0,0)*E_j(2,2)*E_j(3,1) - E_i(1,3)*E_j(0,1)*E_j(2,0)*E_j(3,2)
           + E_i(1,3)*E_j(0,1)*E_j(2,2)*E_j(3,0) + E_i(1,3)*E_j(0,2)*E_j(2,0)*E_j(3,1) - E_i(1,3)*E_j(0,2)*E_j(2,1)*E_j(3,0)
           + E_i(2,0)*E_j(0,1)*E_j(1,2)*E_j(3,3) - E_i(2,0)*E_j(0,1)*E_j(1,3)*E_j(3,2) - E_i(2,0)*E_j(0,2)*E_j(1,1)*E_j(3,3)
           + E_i(2,0)*E_j(0,2)*E_j(1,3)*E_j(3,1) + E_i(2,0)*E_j(0,3)*E_j(1,1)*E_j(3,2) - E_i(2,0)*E_j(0,3)*E_j(1,2)*E_j(3,1)
           - E_i(2,1)*E_j(0,0)*E_j(1,2)*E_j(3,3) + E_i(2,1)*E_j(0,0)*E_j(1,3)*E_j(3,2) + E_i(2,1)*E_j(0,2)*E_j(1,0)*E_j(3,3)
           - E_i(2,1)*E_j(0,2)*E_j(1,3)*E_j(3,0) - E_i(2,1)*E_j(0,3)*E_j(1,0)*E_j(3,2) + E_i(2,1)*E_j(0,3)*E_j(1,2)*E_j(3,0)
           + E_i(2,2)*E_j(0,0)*E_j(1,1)*E_j(3,3) - E_i(2,2)*E_j(0,0)*E_j(1,3)*E_j(3,1) - E_i(2,2)*E_j(0,1)*E_j(1,0)*E_j(3,3)
           + E_i(2,2)*E_j(0,1)*E_j(1,3)*E_j(3,0) + E_i(2,2)*E_j(0,3)*E_j(1,0)*E_j(3,1) - E_i(2,2)*E_j(0,3)*E_j(1,1)*E_j(3,0)
           - E_i(2,3)*E_j(0,0)*E_j(1,1)*E_j(3,2) + E_i(2,3)*E_j(0,0)*E_j(1,2)*E_j(3,1) + E_i(2,3)*E_j(0,1)*E_j(1,0)*E_j(3,2)
           - E_i(2,3)*E_j(0,1)*E_j(1,2)*E_j(3,0) - E_i(2,3)*E_j(0,2)*E_j(1,0)*E_j(3,1) + E_i(2,3)*E_j(0,2)*E_j(1,1)*E_j(3,0)
           - E_i(3,0)*E_j(0,1)*E_j(1,2)*E_j(2,3) + E_i(3,0)*E_j(0,1)*E_j(1,3)*E_j(2,2) + E_i(3,0)*E_j(0,2)*E_j(1,1)*E_j(2,3)
           - E_i(3,0)*E_j(0,2)*E_j(1,3)*E_j(2,1) - E_i(3,0)*E_j(0,3)*E_j(1,1)*E_j(2,2) + E_i(3,0)*E_j(0,3)*E_j(1,2)*E_j(2,1)
           + E_i(3,1)*E_j(0,0)*E_j(1,2)*E_j(2,3) - E_i(3,1)*E_j(0,0)*E_j(1,3)*E_j(2,2) - E_i(3,1)*E_j(0,2)*E_j(1,0)*E_j(2,3)
           + E_i(3,1)*E_j(0,2)*E_j(1,3)*E_j(2,0) + E_i(3,1)*E_j(0,3)*E_j(1,0)*E_j(2,2) - E_i(3,1)*E_j(0,3)*E_j(1,2)*E_j(2,0)
           - E_i(3,2)*E_j(0,0)*E_j(1,1)*E_j(2,3) + E_i(3,2)*E_j(0,0)*E_j(1,3)*E_j(2,1) + E_i(3,2)*E_j(0,1)*E_j(1,0)*E_j(2,3)
           - E_i(3,2)*E_j(0,1)*E_j(1,3)*E_j(2,0) - E_i(3,2)*E_j(0,3)*E_j(1,0)*E_j(2,1) + E_i(3,2)*E_j(0,3)*E_j(1,1)*E_j(2,0)
           + E_i(3,3)*E_j(0,0)*E_j(1,1)*E_j(2,2) - E_i(3,3)*E_j(0,0)*E_j(1,2)*E_j(2,1) - E_i(3,3)*E_j(0,1)*E_j(1,0)*E_j(2,2)
           + E_i(3,3)*E_j(0,1)*E_j(1,2)*E_j(2,0) + E_i(3,3)*E_j(0,2)*E_j(1,0)*E_j(2,1) - E_i(3,3)*E_j(0,2)*E_j(1,1)*E_j(2,0);
  double e = E_j(0,0)*E_j(1,1)*E_j(2,2)*E_j(3,3) - E_j(0,0)*E_j(1,1)*E_j(2,3)*E_j(3,2) - E_j(0,0)*E_j(1,2)*E_j(2,1)*E_j(3,3)
           + E_j(0,0)*E_j(1,2)*E_j(2,3)*E_j(3,1) + E_j(0,0)*E_j(1,3)*E_j(2,1)*E_j(3,2) - E_j(0,0)*E_j(1,3)*E_j(2,2)*E_j(3,1)
           - E_j(0,1)*E_j(1,0)*E_j(2,2)*E_j(3,3) + E_j(0,1)*E_j(1,0)*E_j(2,3)*E_j(3,2) + E_j(0,1)*E_j(1,2)*E_j(2,0)*E_j(3,3)
           - E_j(0,1)*E_j(1,2)*E_j(2,3)*E_j(3,0) - E_j(0,1)*E_j(1,3)*E_j(2,0)*E_j(3,2) + E_j(0,1)*E_j(1,3)*E_j(2,2)*E_j(3,0)
           + E_j(0,2)*E_j(1,0)*E_j(2,1)*E_j(3,3) - E_j(0,2)*E_j(1,0)*E_j(2,3)*E_j(3,1) - E_j(0,2)*E_j(1,1)*E_j(2,0)*E_j(3,3)
           + E_j(0,2)*E_j(1,1)*E_j(2,3)*E_j(3,0) + E_j(0,2)*E_j(1,3)*E_j(2,0)*E_j(3,1) - E_j(0,2)*E_j(1,3)*E_j(2,1)*E_j(3,0)
           - E_j(0,3)*E_j(1,0)*E_j(2,1)*E_j(3,2) + E_j(0,3)*E_j(1,0)*E_j(2,2)*E_j(3,1) + E_j(0,3)*E_j(1,1)*E_j(2,0)*E_j(3,2)
           - E_j(0,3)*E_j(1,1)*E_j(2,2)*E_j(3,0) - E_j(0,3)*E_j(1,2)*E_j(2,0)*E_j(3,1) + E_j(0,3)*E_j(1,2)*E_j(2,1)*E_j(3,0);

  // normalize coefficients w.r.t. first coefficient a for better conditioning
  if(abs(a) > 1.e-16)
  {
    const double inverse = 1./a;
    a = 1.;
    b *= inverse;
    c *= inverse;
    d *= inverse;
    e *= inverse;
  }

  // compute discriminant of characteristic equation
  const double discriminant = 256.*a*a*a*e*e*e - 192.*a*a*b*d*e*e - 128.*a*a*c*c*e*e + 144.*a*a*c*d*d*e
                            -  27.*a*a*d*d*d*d + 144.*a*b*b*c*e*e -   6.*a*b*b*d*d*e -  80.*a*b*c*c*d*e
                            +  18.*a*b*c*d*d*d +  16.*a*c*c*c*c*e -   4.*a*c*c*c*d*d -  27.*b*b*b*b*e*e
                            +  18.*b*b*b*c*d*e -   4.*b*b*b*d*d*d -   4.*b*b*c*c*c*e +      b*b*c*c*d*d;

  // compute auxiliary quantity
  const double P = 8.*a*c - 3.*b*b;

  // collision check via discriminant and auxiliary quantity
  return (discriminant > 0. and P < 0.) ? false : true;
  END OF ALTERNATIVE METHOD */
}


/*----------------------------------------------------------------------------*
 | compute contact point on particle surface                       fang 10/17 |
 *----------------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerDEMEllipsoids::ComputeContactPoint(
    LINALG::Matrix<1,4>&         C_particle,   //!< contact point on particle surface
    const LINALG::Matrix<3,1>&   semiaxes,     //!< particle semi-axes
    const LINALG::Matrix<4,4>&   T,            //!< translation matrix
    const LINALG::Matrix<4,4>&   R,            //!< rotation matrix
    const LINALG::Matrix<1,4>&   C,            //!< contact point inside contact plane
    const LINALG::Matrix<1,4>&   d             //!< direction vector of contact line
    ) const
{
  // procedure:
  // 1.) affinely transform the ellipsoidal particle to a unit sphere centered at the origin
  //     of the coordinate system, employing homogeneous coordinates and the following steps:
  //     a) translate each point X on the surface of the original ellipsoid by the same distance,
  //        such that the center of the ellipsoid coincides with the origin of the coordinate system,
  //        and obtain the translated point Xt = X * T on the surface of the translated ellipsoid
  //     b) rotate the translated ellipsoid, such that its semi-axes are parallel to the coordinate axes,
  //        and thereby transform each point Xt into Xtt = Xt * R = Xt * T * R
  //     c) scale the translated and rotated ellipsoid, such that all of its semi-axes have unit length,
  //        and thereby transform each point Xtt into X' = Xtt * S = Xt * R * S = X * T * R * S = X * M,
  //        where M = T * R * S is the overall transformation matrix
  // 2.) apply the same transformation to the contact line L = C + t * d = C + t * (V3 - V2), where C is
  //     the contact point inside the contact plane, t is the line parameter, d is the direction vector,
  //     and V2 and V3 are the two last eigenvectors of the generalized eigenproblem
  // 3.) compute the two intersection points between the unit sphere and the transformed line L', given by
  //     L' = C' + t * d' = C' + t * (V3' - V2')
  // 4.) transform back the point closer to C' to obtain the sought contact point on the particle surface

  // step 1.)

  // compute overall transformation matrix M
  static LINALG::Matrix<4,4> M;
  ComputeTransformationMatrix(M,semiaxes,T,R);

  // step 2.)

  // transform contact point C
  static LINALG::Matrix<1,4> C_t;
  C_t.Multiply(C,M);
  C_t(3) = 0.;

  // transform direction vector d
  static LINALG::Matrix<1,4> d_t;
  d_t.Multiply(d,M);

  // normalize transformed direction vector d'
  d_t.Scale(1./d_t.Norm2());

  // step 3.)
  // each of the two intersection points lies on both the surface of the unit sphere and the transformed line
  // we thus solve the quadratic equation || C' + t * d' ||^2 = 1 for the line parameter t

  // precompute components
  double Cd = C_t.Dot(d_t);
  double root = sqrt(Cd*Cd-(C_t.Dot(C_t)-1.));

  // determine solutions of quadratic equation
  double t1 = -Cd+root;
  double t2 = -Cd-root;

  // step 4.)

  // compute intersection point closer to C'
  static LINALG::Matrix<1,4> C_particle_t;
  C_particle_t.Clear();
  double t = abs(t1) < abs(t2) ? t1 : t2;
  for(int i=0; i<3; ++i)
    C_particle_t(i) = C_t(i)+t*d_t(i);
  C_particle_t(3) = 1.;

  // transform intersection point back
  LINALG::FixedSizeSerialDenseSolver<4,4> inverter;
  inverter.SetMatrix(M);
  if(inverter.Invert())
    dserror("Inversion of 4x4 matrix failed!");
  C_particle.Multiply(C_particle_t,M);
  C_particle(3) = 0.;

  return;
}


/*---------------------------------------------------------------------------------------------------*
 | compute eigenvectors of B^-1 * A and determine whether there are complex eigenvalues   fang 10/17 |
 *---------------------------------------------------------------------------------------------------*/
bool PARTICLE::ParticleCollisionHandlerDEMEllipsoids::ComputeEigenVectors(
    const LINALG::Matrix<4,4>&   A,    //!< first matrix of generalized eigenproblem
    const LINALG::Matrix<4,4>&   B,    //!< second matrix of generalized eigenproblem
    LINALG::Matrix<3,1>&         V0,   //!< first eigenvector of generalized eigenproblem
    LINALG::Matrix<3,1>&         V1,   //!< second eigenvector of generalized eigenproblem
    LINALG::Matrix<3,1>&         V2,   //!< third eigenvector of generalized eigenproblem
    LINALG::Matrix<3,1>&         V3    //!< fourth eigenvector of generalized eigenproblem
    ) const
{
  // compute matrix B^-1 * A
  static LINALG::Matrix<4,4> temp4x4;
  temp4x4 = B;
  LINALG::FixedSizeSerialDenseSolver<4,4> inverter;
  inverter.SetMatrix(temp4x4);
  if(inverter.Invert())
    dserror("Inversion of 4x4 matrix failed!");
  static LINALG::Matrix<4,4> Binv_A;
  Binv_A.Multiply(temp4x4,A);

  // compute eigenvectors of matrix B^-1 * A
  LINALG::Matrix<1,4> WR(true);        // real eigenvalue parts
  static LINALG::Matrix<1,4> WI;   // imaginary eigenvalue parts
  double work[16];
  int info;
  Epetra_LAPACK lapack;
  lapack.GEEV('N','V',4,Binv_A.A(),4,WR.A(),WI.A(),NULL,4,temp4x4.A(),4,work,16,&info);
  if(info < 0)
    dserror("LAPACK algorithm GEEV: The %d-th argument has an illegal value!",info);

  // normalize eigenvectors w.r.t. respective last homogeneous coordinates
  for(unsigned i=0; i<4 ; ++i)
  {
    // compute inverse of current last homogeneous coordinate, avoiding division by zero
    const double inverse = 1./(abs(temp4x4(3,i)) > 1.e-16 ? temp4x4(3,i) : 1.e-16);

    // normalize current eigenvector
    for(unsigned j=0; j<4 ; ++j)
      temp4x4(j,i) = temp4x4(j,i)*inverse;
  }

  // sort eigenvectors in ascending order of real parts of associated eigenvalues
  static std::vector<int> indices(4,0);
  std::iota(indices.begin(),indices.end(),0);
  std::sort(indices.begin(),indices.end(),[&WR](int a, int b){return WR(a) < WR(b);});
  for(unsigned i=0; i<3 ; ++i)
  {
    V0(i) = temp4x4(i,indices[0]);
    V1(i) = temp4x4(i,indices[1]);
    V2(i) = temp4x4(i,indices[2]);
    V3(i) = temp4x4(i,indices[3]);
  }

  // determine whether there are complex eigenvalues
  return WI.Norm2() > 1.e-12 ? true : false;
}


/*----------------------------------------------------------------------------*
 | compute ellipsoid matrix expressed in homogeneous coordinates   fang 10/17 |
 *----------------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerDEMEllipsoids::ComputeEllipsoidMatrix(
    LINALG::Matrix<4,4>&         E,          //!< ellipsoid matrix
    const LINALG::Matrix<3,1>&   semiaxes,   //!< ellipsoid semi-axes
    const LINALG::Matrix<4,4>&   T,          //!< translation matrix
    const LINALG::Matrix<4,4>&   R           //!< rotation matrix
    ) const
{
  // clear ellipsoid matrix
  E.Clear();

  // assemble plain ellipsoid matrix without considering translation and rotation
  E(0,0) = 1./(semiaxes(0)*semiaxes(0));
  E(1,1) = 1./(semiaxes(1)*semiaxes(1));
  E(2,2) = 1./(semiaxes(2)*semiaxes(2));
  E(3,3) = -1.;

  // compute ellipsoid matrix expressed in homogeneous coordinates:
  // E = T * R * S * R^T * T^T, where S is the plain ellipsoid matrix without considering translation and rotation
  static LINALG::Matrix<4,4> temp4x4;
  temp4x4.MultiplyNT(E,R);
  E.Multiply(R,temp4x4);
  temp4x4.MultiplyNT(E,T);
  E.Multiply(T,temp4x4);

  return;
}


/*------------------------------------------------------------------------------------------------------*
 | compute contact normal, contact gap, and vectors from particle centers to contact point   fang 10/17 |
 *------------------------------------------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerDEMEllipsoids::ComputeContactPointAndNormalAndGap(
    LINALG::Matrix<3,1>&      r_iC,     //!< vector from center of particle i to contact point C
    LINALG::Matrix<3,1>&      r_jC,     //!< vector from center of particle j to contact point C
    LINALG::Matrix<3,1>&      normal,   //!< contact normal
    double&                   gap,      //!< contact gap
    const ParticleCollData&   data_i,   //!< collision data for particle i
    const ParticleCollData&   data_j    //!< collision data for particle j
    ) const
{
  // initialize matrices for collision detection for ellipsoidal particles (involving homogeneous coordinates)
  static LINALG::Matrix<4,4> T_i, T_j;   // translation matrices
  static LINALG::Matrix<4,4> R_i, R_j;   // rotation matrices
  static LINALG::Matrix<4,4> E_i, E_j;   // final ellipsoid matrices

  // assemble translation, rotation, and ellipsoid matrices for ellipsoids i and j
  ComputeTranslationMatrix(T_i,data_i.dis);
  ComputeRotationMatrix(R_i,data_i.ang);
  ComputeEllipsoidMatrix(E_i,data_i.semiaxes,T_i,R_i);
  ComputeTranslationMatrix(T_j,data_j.dis);
  ComputeRotationMatrix(R_j,data_j.ang);
  ComputeEllipsoidMatrix(E_j,data_j.semiaxes,T_j,R_j);

  // perform collision check
  if(CollisionCheck(E_i,E_j))
  {
    // collision case: gradually shrink ellipsoids until they are separated, then compute eigenvectors to determine contact point
    // initialize scaling factor
    double W(1.);

    // begin shrinking procedure
    do
    {
      // step size for scaling factor modification
      const double W_delta(1.e-2);

      // decrease scaling factor to shrink ellipsoids
      W -= W_delta;

      // safety check
      if(W < .9)
        dserror("Ellipsoids shrunk more than 10%%... maybe something is wrong...");

      // update matrices for shrunk ellipsoids
      // note: only the bottom right entry of either E_i or E_j depends on the scaling factor W,
      //       according to [... some terms independent of W ...] - W^2
      const double E_3_3_delta = 2.*W*W_delta-W_delta*W_delta;
      E_i(3,3) += E_3_3_delta;
      E_j(3,3) += E_3_3_delta;
    }

    // perform collision check for shrunk ellipsoids
    while(CollisionCheck(E_i,E_j));

    // solve generalized eigenproblem as soon as shrunk ellipsoids no longer collide:
    //     lambda * E_i * V + E_j * V = 0   or   lambda * E_j * V + E_i * V = 0
    // --> compute eigenvectors V0, V1, V2, and V3 of:
    //     ((-E_i)^-1 * E_j)   or   (E_j^-1 * (-E_i))
    // although the eigenvectors of both matrices are analytically identical, computationally they are not,
    // and therefore we check the global IDs of particles i and j to always solve the same eigenproblem,
    // no matter how many processors are employed
    static LINALG::Matrix<3,1> V0, V1, V2, V3;
    E_i.Scale(-1.);
    if(data_i.id < data_j.id ? ComputeEigenVectors(E_j,E_i,V0,V1,V2,V3) : ComputeEigenVectors(E_i,E_j,V0,V1,V2,V3))
      dserror("Shrunk ellipsoids still collide!");

    // compute contact point inside contact plane
    static LINALG::Matrix<3,1> C;
    C.Update(V2,V3);
    C.Scale(.5);

    // compute vectors from particle centers to contact point
    r_iC.Update(1.,C,-1.,data_i.dis);
    r_jC.Update(1.,C,-1.,data_j.dis);

    // compute normal vector of contact plane
    V0 -= C;
    V1 -= C;
    normal.CrossProduct(V1,V0);
    normal.Scale(1./normal.Norm2());

    // check whether normal vector points from particle i to particle j, and flip if not
    static LINALG::Matrix<3,1> r_ij;
    r_ij.Update(-1.,data_i.dis,1.,data_j.dis);
    if(normal.Dot(r_ij) < 0)
      normal.Scale(-1.);

    // compute contact points on the surfaces of particles i and j:
    // 1.) define the line L = C + t * d = C + t * (V3 - V2), where C is the contact point inside the contact plane,
    //     t is the line parameter, d is the direction vector, and V2 and V3 are the two last eigenvectors of the
    //     generalized eigenproblem
    // 2.) compute contact point on the surface of particle i by intersecting particle i and line L after an affine
    //     transformation
    // 3.) repeat step 2.) for particle j to obtain the sought contact point on the surface of particle j

    // step 1.)
    // express contact point C in homogeneous coordinates
    // the last homogeneous coordinate of any point is unity by definition
    static LINALG::Matrix<1,4> C_hom;
    for(int i=0; i<3; ++i)
      C_hom(i) = C(i);
    C_hom(3) = 1.;

    // express direction vector d = V3 - V2 of line L in homogeneous coordinates
    // note that the last homogeneous coordinate is zero, since d is a vector and not a point
    // think of d as the difference of the two points V3 and V2, each with last homogeneous coordinate 1
    static LINALG::Matrix<1,4> d_hom;
    for(int i=0; i<3; ++i)
      d_hom(i)= V3(i)-V2(i);
    d_hom(3) = 0.;

    // step 2.)
    static LINALG::Matrix<1,4> C_i;
    ComputeContactPoint(C_i,data_i.semiaxes,T_i,R_i,C_hom,d_hom);

    // step 3.)
    static LINALG::Matrix<1,4> C_j;
    ComputeContactPoint(C_j,data_j.semiaxes,T_j,R_j,C_hom,d_hom);

    // compute contact gap: gap = || C_i - C_j ||
    C_i -= C_j;
    gap = -C_i.Norm2();
  }

  else
    // non-collision case: set gap to infinity
    gap = std::numeric_limits<double>::infinity();

  return;
}


/*---------------------------------------------------------------------------------------------------------------------*
 | compute contact normal, contact gap, contact point on wall element, and type of nearest contact object   fang 10/17 |
 *---------------------------------------------------------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerDEMEllipsoids::ComputeContactPointAndNormalAndGapAndObjectType(
    LINALG::Matrix<3,1>&                        C,              //!< contact point on wall element
    GEO::ObjectType&                            objecttype,     //!< type of nearest contact object
    LINALG::Matrix<3,1>&                        normal,         //!< contact normal
    double&                                     gap,            //!< contact gap
    DRT::Element*                               wallele,        //!< wall element
    const std::map<int,LINALG::Matrix<3,1> >&   nodecoord,      //!< node coordinates of wall element
    const ParticleCollData&                     data            //!< particle collision data
    ) const
{
  // for ellipsoidal particles, particle-wall contact is evaluated as follows:
  // 1.) affinely transform a given ellipsoid to a unit sphere centered at the origin of the coordinate system,
  //     employing homogeneous coordinates
  // 2.) apply the same affine transformation to the wall element
  // 3.) project the center of the unit sphere onto the transformed wall element to obtain the contact point
  //     on the transformed wall element
  // 4.) compute the contact point on the surface of the unit sphere, i.e., the intersection point between
  //     the surface of the unit sphere and the line through the center of the unit sphere and the contact point
  //     on the transformed wall element
  // 5.) transform back the contact point on the surface of the unit sphere to obtain the contact point on the
  //     surface of the original ellipsoid
  // 6.) project the contact point on the surface of the original ellipsoid onto the original wall element to
  //     obtain the contact point on the original wall element
  // 7.) calculate the contact normal and the contact gap from the contact points on the original wall element
  //     and the surface of the original ellipsoid

  // step 1.)
  // compute transformation matrix M in homogeneous coordinates
  static LINALG::Matrix<4,4> T, R, M;
  ComputeTranslationMatrix(T,data.dis);
  ComputeRotationMatrix(R,data.ang);
  ComputeTransformationMatrix(M,data.semiaxes,T,R);

  // step 2.)
  // transform the node coordinates of the wall element
  std::map<int,LINALG::Matrix<3,1> > nodecoord_t;
  static LINALG::Matrix<3,1> coord, coord_t;
  static LINALG::Matrix<1,4> coord_hom, coord_hom_t;
  for(int inode=0; inode<wallele->NumNode(); ++inode)
  {
    const LINALG::Matrix<3,1>& currnodecoord = nodecoord.at(wallele->Nodes()[inode]->Id());
    for(int i=0; i<3; ++i)
      coord_hom(i) = currnodecoord(i);
    coord_hom(3) = 1.;
    coord_hom_t.Multiply(coord_hom,M);
    for(int i=0; i<3; ++i)
      coord_t(i) = coord_hom_t(i);
    nodecoord_t[wallele->Nodes()[inode]->Id()] = coord_t;
  }

  // step 3.)
  coord_t.Clear();
  GEO::nearest3DObjectOnElement(wallele,nodecoord_t,coord_t,C);

  // step 4.)
  const double norm1 = C.Norm2();
  C.Scale(1./norm1);

  // step 5.)
  for(int i=0; i<3; ++i)
    coord_hom_t(i) = C(i);
  coord_hom_t(3) = 1.;
  LINALG::FixedSizeSerialDenseSolver<4,4> inverter;
  inverter.SetMatrix(M);
  if(inverter.Invert())
    dserror("Inversion of 4x4 matrix failed!");
  coord_hom.Multiply(coord_hom_t,M);
  for(int i=0; i<3; ++i)
    coord(i) = coord_hom(i);

  // step 6.)
  objecttype = GEO::nearest3DObjectOnElement(wallele,nodecoord,coord,C);
  if(objecttype != GEO::SURFACE_OBJECT)
    dserror("For ellipsoidal particles, particle-wall contact is only implemented for surfaces, but not for lines and corners!");

  // step 7.)
  normal.Update(1.,coord,-1.,C);
  const double norm2 = normal.Norm2();
  normal.Scale(1./norm2);
  gap = norm1 < 1. ? -norm2 : norm2;

  return;
}


/*---------------------------------------------------------------------------*
 | compute rotation matrix expressed in homogeneous coordinates   fang 10/17 |
 *---------------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerDEMEllipsoids::ComputeRotationMatrix(
    LINALG::Matrix<4,4>&         R,      //!< rotation matrix
    const LINALG::Matrix<3,1>&   angle   //!< rotation angle vector
    ) const
{
  // clear rotation matrix
  R.Clear();

  // compute rotation matrix expressed in homogeneous coordinates
  static LINALG::Matrix<3,3> temp3x3;
  LARGEROTATIONS::angletotriad(angle,temp3x3);
  for(int i=0; i<3; ++i)
    for(int j=0; j<3; ++j)
      R(i,j) = temp3x3(i,j);
  R(3,3) = 1.;

  return;
}


/*------------------------------------------------------------------------------------------------------------------------*
 | compute matrix transforming an ellipsoid to a unit sphere centered at the origin of the coordinate system   fang 10/17 |
 *------------------------------------------------------------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerDEMEllipsoids::ComputeTransformationMatrix(
    LINALG::Matrix<4,4>&         M,          //!< transformation matrix
    const LINALG::Matrix<3,1>&   semiaxes,   //!< ellipsoid semi-axes
    const LINALG::Matrix<4,4>&   T,          //!< translation matrix
    const LINALG::Matrix<4,4>&   R           //!< rotation matrix
    ) const
{
  // clear transformation matrix
  M.Clear();

  // assemble scaling matrix to affinely transform a plain ellipsoid to a plain unit sphere
  M(0,0) = 1./(semiaxes(0));
  M(1,1) = 1./(semiaxes(1));
  M(2,2) = 1./(semiaxes(2));
  M(3,3) = 1.;

  // multiply rotation matrix R by scaling matrix
  static LINALG::Matrix<4,4> temp4x4;
  temp4x4.Multiply(R,M);

  // multiply translation matrix T by result R * M
  M.Multiply(T,temp4x4);

  return;
}


/*------------------------------------------------------------------------------*
 | compute translation matrix expressed in homogeneous coordinates   fang 10/17 |
 *------------------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerDEMEllipsoids::ComputeTranslationMatrix(
    LINALG::Matrix<4,4>&         T,             //!< translation matrix
    const LINALG::Matrix<3,1>&   displacement   //!< displacement vector
    ) const
{
  // clear translation matrix
  T.Clear();

  // compute translation matrix expressed in homogeneous coordinates
  for(int i=0; i<4; ++i)
    T(i,i) = 1.;
  for(int dim=0; dim<3; ++dim)
    T(3,dim) = -displacement(dim);

  return;
}


/*-----------------------------------------------------------------------------------------------------------*
 | determine whether the generalized eigenproblem A*x = lambda B*x exhibits complex eigenvalues   fang 10/17 |
 *-----------------------------------------------------------------------------------------------------------*/
bool PARTICLE::ParticleCollisionHandlerDEMEllipsoids::ComplexEigenValues(
    const LINALG::Matrix<4,4>&   A,   //!< first matrix of generalized eigenproblem
    const LINALG::Matrix<4,4>&   B    //!< second matrix of generalized eigenproblem
    ) const
{
  // set number of matrix rows or columns
  int N(4);

  // copy matrices A and B
  Epetra_SerialDenseMatrix A_temp(N,N), B_temp(N,N);
  for(int i=0; i<N; ++i)
    for(int j=0; j<N; ++j)
    {
      A_temp(i,j) = A(i,j);
      B_temp(i,j) = B(i,j);
    }

  //
  // step 1: transform matrix B to upper triangluar matrix via QR factorization
  //

  int jpvt[N];     // order of permutation matrix
  double tau[N];   // factor for calculation of orthogonal matrix Q
  int lwork(3*N+1), info;
  double dummy_lwork[lwork];
  dgeqp3(&N,&N,B_temp.A(),&N,jpvt,tau,dummy_lwork,&lwork,&info);
  if(info < 0)
    dserror("LAPACK algorithm dgeqp3: The %d-th argument has an illegal value!",info);

  // calculate matrix Q from multiplying Householder transformations
  Epetra_SerialDenseMatrix Q(N,N), temp(N,N);
  for(int i=0; i<N; ++i)
    Q(i,i) = 1.;
  for(int i=0; i<N; ++i)
  {
    Epetra_SerialDenseVector v(N);
    v(i) = 1.;
    for(int j=i+1; j<N; ++j)
      v(j) = B_temp(j,i);
    Epetra_SerialDenseMatrix H(N,N);
    H.Multiply('N','T',tau[i],v,v,0.);
    H.Scale(-1.);
    for(int j=0; j<N; ++j)
      H(j,j) += 1.;
    Q.Apply(H,temp);
    Q = temp;
  }

  // calculate permutation matrix P
  Epetra_SerialDenseMatrix P(N,N);
  for(int i=0; i<N; ++i)
    P(jpvt[i]-1,i) = 1.;

  // annul lower diagonal entries of B
  for(int i=0; i<N; ++i)
    for(int j=i+1; j<N; ++j)
      B_temp(j,i) = 0.;

  // transform matrix A according to Q^T * A * P
  temp.Multiply('T','N',1.,Q,A_temp,0.);
  A_temp.Multiply('N','N',1.,temp,P,0.);

  //
  // step 2: transform matrix A to upper Hessenberg matrix and keep matrix B as upper diagonal matrix
  //

  char job = 'P';
  int ilo, ihi;
  double dummy_N[N], dummy_6N[6*N];
  dggbal(&job,&N,A_temp.A(),&N,B_temp.A(),&N,&ilo,&ihi,dummy_N,dummy_N,dummy_6N,&info);
  if(info)
    dserror("LAPACK algorithm dggbal threw error code %d!",info);

  char comp = 'I';
  Epetra_SerialDenseMatrix Z(N,N);
  dgghrd(&comp,&comp,&N,&ilo,&ihi,A_temp.A(),&N,B_temp.A(),&N,Q.A(),&N,Z.A(),&N,&info);
  if(info < 0)
    dserror("LAPACK algorithm dgghrd: The %d-th argument has an illegal value!",info);

  //
  // step 3: transform upper Hessenberg matrix A to upper diagonal matrix and keep matrix B as upper diagonal matrix via QZ transformation
  //

  // vectors containing eigenvalues of generalized eigenproblem
  double alphar[N], alphai[N], beta[N];

  job ='E';
  comp = 'V';
  dhgeqz(&job,&comp,&comp,&N,&ilo,&ihi,A_temp.A(),&N,B_temp.A(),&N,alphar,alphai,beta,Q.A(),&N,Z.A(),&N,dummy_N,&N,&info);
  if(info < 0)
    dserror("LAPACK algorithm dhgeqz: The %d-th argument has an illegal value!",info);
  else if(info > N)
    dserror("LAPACK algorithm dhgeqz: The QZ iteration did not converge! (H,T) is not in Schur form, but the eigenvalues should be correct!");

  bool complex(false);
  for(int i=0; i<N; ++i)
    if(abs(alphai[i]) > 1.e-16)
      complex = true;

  return complex;
}


/*----------------------------------------------------------------------*
 | constructor for MD based particle contact               ghamm 09/13  |
 *----------------------------------------------------------------------*/
PARTICLE::ParticleCollisionHandlerMD::ParticleCollisionHandlerMD(
  Teuchos::RCP<DRT::Discretization> discret,
  Teuchos::RCP<PARTICLE::Algorithm> particlealgorithm,
  const Teuchos::ParameterList& particledynparams
  ) :
  PARTICLE::ParticleCollisionHandlerBase(
    discret,
    particlealgorithm,
    particledynparams
    )
{
  // safety check
  if(particlealgorithm->BinStrategy()->HavePBC())
    dserror("Periodic boundary conditions not yet implemented for molecular dynamics!");

  return;
}


/*----------------------------------------------------------------------*
 | compute collisions (inter-particle and particle-wall)   ghamm 09/13  |
 *----------------------------------------------------------------------*/
double PARTICLE::ParticleCollisionHandlerMD::EvaluateParticleContact(
  double dt,
  Teuchos::RCP<Epetra_Vector> disn,
  Teuchos::RCP<Epetra_Vector> veln,
  Teuchos::RCP<Epetra_Vector> specEnthalpyn,
  Teuchos::RCP<Epetra_FEVector> f_structure)
{

  // ddt is the most advanced collision time (between [0,dt))
  double ddt = 0.0;
  std::set<Teuchos::RCP<Event>, Event::Helper> eventqueue;
  // initializing the event queue with future collisions

  InitializeEventQueue(eventqueue, dt);



  // collisions are still evaluated at the end of the time step
  while ((!eventqueue.empty()) and ((ddt - dt) <= GEO::TOL14 or ddt < dt))
  {
    // find the next event in the eventqueue and deal with it
    Teuchos::RCP<Event> next_event = *eventqueue.begin();

#ifdef OUTPUT
    std::cout << " Dealing with the collision event at time " << next_event->time
        << " between the two collision partners " << next_event->particle_1->Id() << " & ";
    if (next_event->coltype == INPAR::PARTICLE::particle_particle)
    {
      std::cout << "GID of second particle " << next_event->particle_2->Id() << std::endl;
    }
    else if (next_event->coltype == INPAR::PARTICLE::particle_wall)
    {
      std::cout << "GID of wall " << Teuchos::rcp_dynamic_cast<WallEvent>(next_event)->wall->Id() << std::endl;
    }
#endif

    // time variable is set to collision time of the next event
    ddt = next_event->time;

    // collision still happening in this time step?
    if ((ddt - dt) < GEO::TOL14 or ddt < dt)
    {
      // HandleCollision computes new velocities and updates position of colliding particles
      HandleCollision(next_event, dt);

      // Updating the event queue
#ifdef OUTPUT
      std::cout << " ErasingInvalidCollisions " << std::endl;
#endif

      // loop over event queue and erase invalid collisions
      for (std::set<Teuchos::RCP<Event>, Event::Helper>::iterator iter=eventqueue.begin(); iter!=eventqueue.end(); /*no ++iter*/)
      {
#ifdef OUTPUT
        Teuchos::RCP<Event> output1 = *iter;
        std::cout << "GID of first particle " << output1->particle_1->Id() << std::endl;
        if (output1->coltype == INPAR::PARTICLE::particle_particle)
        {
          std::cout << "GID of second particle " << output1->particle_2->Id() << std::endl;
        }
        else if (output1->coltype == INPAR::PARTICLE::particle_wall)
        {
          std::cout << "GID of wall " << Teuchos::rcp_dynamic_cast<WallEvent>(output1)->wall->Id() << std::endl;
        }
#endif
        Teuchos::RCP<Event> invalid_col = *iter;
        if (next_event->coltype == INPAR::PARTICLE::particle_particle)
        {
          if ((invalid_col->coltype == INPAR::PARTICLE::particle_particle) and (invalid_col->particle_1->Id() == next_event->particle_1->Id()
                  or invalid_col->particle_1->Id() == next_event->particle_2->Id()
                  or invalid_col->particle_2->Id() == next_event->particle_1->Id()
                  or invalid_col->particle_2->Id() == next_event->particle_2->Id()))
          {
#ifdef OUTPUT
            std::cout << " Erasing the following elements from the eventqueue: " << std::endl;
            std::cout << " happening at time:  " << invalid_col->time << std::endl;
            std::cout << "GID of first particle " << invalid_col->particle_1->Id() << std::endl;
            std::cout << "GID of second particle " << invalid_col->particle_2->Id() << std::endl;
#endif

            // go to next element and erase invalid one
            std::set<Teuchos::RCP<Event>, Event::Helper>::iterator tmp = iter;
            ++iter;
            eventqueue.erase(tmp);
          }
          else if ((invalid_col->coltype == INPAR::PARTICLE::particle_wall) and (invalid_col->particle_1->Id() == next_event->particle_1->Id()
                  or invalid_col->particle_1->Id() == next_event->particle_2->Id()))
          {
#ifdef OUTPUT
            std::cout << " Erasing the following elements from the eventqueue: " << std::endl;
            std::cout << " happening at time:  " << invalid_col->time << std::endl;
            std::cout << "GID of first particle " << invalid_col->particle_1->Id() << std::endl;
            std::cout << "GID of wall " << Teuchos::rcp_dynamic_cast<WallEvent>(invalid_col)->wall->Id() << std::endl;
#endif

            std::set<Teuchos::RCP<Event>, Event::Helper>::iterator tmp = iter;
            ++iter;
            eventqueue.erase(tmp);
          }
          else
          {
            ++iter;
          }
        }
        else if (next_event->coltype == INPAR::PARTICLE::particle_wall)
        {
          if ((invalid_col->coltype == INPAR::PARTICLE::particle_particle) and (invalid_col->particle_1->Id() == next_event->particle_1->Id()
                  or invalid_col->particle_2->Id() == next_event->particle_1->Id()))

          {
#ifdef OUTPUT
            std::cout << " Erasing the following elements from the eventqueue: " << std::endl;
            std::cout << " happening at time:  " << invalid_col->time << std::endl;
            std::cout << "GID of first particle " << invalid_col->particle_1->Id() << std::endl;
            std::cout << "GID of second particle " << invalid_col->particle_2->Id() << std::endl;
#endif

            std::set<Teuchos::RCP<Event>, Event::Helper>::iterator tmp = iter;
            ++iter;
            eventqueue.erase(tmp);
          }
          else if ((invalid_col->coltype == INPAR::PARTICLE::particle_wall) and (invalid_col->particle_1->Id() == next_event->particle_1->Id()))
          {
#ifdef OUTPUT
            std::cout << " Erasing the following elements from the eventqueue: " << std::endl;
            std::cout << " happening at time:  " << invalid_col->time << std::endl;
            std::cout << "GID of first particle " << invalid_col->particle_1->Id() << std::endl;
            std::cout << "GID of wall " << Teuchos::rcp_dynamic_cast<WallEvent>(invalid_col)->wall->Id() << std::endl;
#endif

            std::set<Teuchos::RCP<Event>, Event::Helper>::iterator tmp = iter;
            ++iter;
            eventqueue.erase(tmp);
          }
          else
          {
            ++iter;
          }
        }
        else
          dserror("you should not show up here");
      }

#ifdef OUTPUT
        std::cout << "The eventqueue contains after deleting invalid events: ";
        int count = 0;
        for (std::set<Teuchos::RCP<Event>, Event::Helper>::iterator iter2 = eventqueue.begin(); iter2 != eventqueue.end(); ++iter2)
        {
          std::cout << "This is entry number " << count << " of the eventqueue:" << std::endl;

          Teuchos::RCP<Event> output2 = *iter2;
          std::cout << " happening at time:  " << output2->time << std::endl;
          std::cout << "GID of first particle " << output2->particle_1->Id() << std::endl;
          if (output2->coltype == INPAR::PARTICLE::particle_particle)
          {
            std::cout << "GID of second particle " << output2->particle_2->Id() << std::endl;
          }
          else if (output2->coltype == INPAR::PARTICLE::particle_wall)
          {
            std::cout << "GID of wall " << Teuchos::rcp_dynamic_cast<WallEvent>(output2)->wall->Id() << std::endl;
          }
          count++;
        }
        std::cout << std::endl;
#endif

      SearchForNewCollisions(next_event, eventqueue, dt);
    }

  }


  // updating of all particles to the end of the current time step
  for (int i=0; i<discret_->NodeColMap()->NumMyElements(); ++i)
    (particleData_[i].dis).Update(dt - particleData_[i].ddt,particleData_[i].vel,1);

  // copy values from col to row layout
  for(int i=0; i<discret_->NumMyColNodes(); ++i)
  {
    if (particleData_[i].owner == myrank_)
    {
      for (int dim=0;dim<3;++dim)
      {
        int lidDof = disn->Map().LID((particleData_[i].lm)[dim]);

        (*disn)[lidDof] = particleData_[i].dis(dim);
        (*veln)[lidDof] = particleData_[i].vel(dim);
      }
    }
  }

#ifdef DEBUG
  // test for inter-particle penetration --> endless loops can be caused by overlapping particles
  for (int i=0; i<discret_->NodeColMap()->NumMyElements(); ++i)
  {
    DRT::Node* currparticle = discret_->lColNode(i);

    static LINALG::Matrix<3,1> currposition;
    static LINALG::Matrix<3,1> currvelocity;
    double currradius;
    double currtime;

    GetCollisionData(currparticle, currposition, currvelocity, currradius, currtime);

    // gather all particles and walls in the vicinity of currparticle
    std::list<DRT::Node*> neighboring_particles;
    boost::unordered_map<int, DRT::Element*> neighboring_walls;
    particle_algorithm_->GetNeighbouringItems(currparticle, neighboring_particles, neighboring_walls);

    // loop over all neighbouring particles and check if the sum of their radii is larger than their distance
    for (std::list<DRT::Node*>::iterator iter=neighboring_particles.begin(); iter!=neighboring_particles.end(); ++iter)
    {
      if (currparticle->Id() >= (*iter)->Id())
      {
        continue;
      }

      static LINALG::Matrix<3,1> iterposition;
      static LINALG::Matrix<3,1> itervelocity;
      double iterradius;
      double itertime;

      GetCollisionData(*iter, iterposition, itervelocity, iterradius, itertime);

      static LINALG::Matrix<3,1> distance;
      distance.Update(1.0, currposition, -1.0, iterposition);
      if (distance.Norm2() - (currradius + iterradius) < -GEO::TOL14)
      {
        dserror("Particles penetrated!");
      }
    }

    // test for particle-wall penetration
    for (boost::unordered_map<int, DRT::Element*>::iterator iter=neighboring_walls.begin(); iter!=neighboring_walls.end(); ++iter)
    {
      Teuchos::RCP<DRT::Discretization> walldiscret = particle_algorithm_->WallDiscret();
      Teuchos::RCP<const Epetra_Vector> walldisnp = walldiscret->GetState("walldisnp");

      DRT::Element* currElem = iter->second;
      DRT::Node** nodes = currElem->Nodes();
      int numnodes = currElem->NumNode();

      std::vector<int> lm;
      lm.reserve(numnodes*3);

      std::vector<int> lmowner;
      std::vector<int> lmstride;
      currElem->LocationVector(*walldiscret, lm, lmowner, lmstride);

      std::map<int, LINALG::Matrix<3,1> > wallpositions;
      // nodal position of wall at the end of the time step
      std::vector<double> nodaldisnp(numnodes*3);
      DRT::UTILS::ExtractMyValues(*walldisnp, nodaldisnp, lm);

      for (int j=0; j<numnodes; ++j)
      {
        const DRT::Node* node = nodes[j];

        static LINALG::Matrix<3,1> nodepos;
        for (int i=0; i<3; ++i)
        {
          nodepos(i) = node->X()[i] + nodaldisnp[3*j + i];
        }
        wallpositions[node->Id()] = nodepos;
      }

      LINALG::Matrix<3,1> minDistCoords;
      GEO::nearest3DObjectOnElement(currElem, wallpositions, currposition, minDistCoords);
      static LINALG::Matrix<3,1> distance;
      distance.Update(1.0, currposition, -1.0, minDistCoords);

      if (distance.Norm2() < (currradius - GEO::TOL14))
      {
        std::cout << "particle " << currparticle->Id() << std::endl;
        std::cout << "wall element " << currElem->Id() << std::endl;
        std::cout << "distance " << distance.Norm2() << std::endl;
        std::cout << "currentradius " << currradius << std::endl;
        std::cout << "particle is penetrating the wall" << std::endl;
        dserror("Particle is penetrating the wall");
      }
    }
  }
#endif

  // reset vector because it needs to be rebuild in the next time step
  particleData_.clear();

  // no internal energy is stored in contact because no overlap is allowed
  return 0.0;
}


/*----------------------------------------------------------------------*
 | contact for one event is evaluated                      ghamm 09/13  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerMD::HandleCollision(
  Teuchos::RCP<Event> next_event,
  const double dt
  )
{
  // STARTING FROM HERE:
  // event.time in event queue is related to beginning of time step

  // check for Particle Particle Collision or Particle Wall Collision
  switch(next_event->coltype)
  {
  case INPAR::PARTICLE::particle_particle:
  {
#ifdef OUTPUT
    std::cout << "Handle collision of particle " << next_event->particle_1->Id()
        << " with partice " << next_event->particle_2->Id() << std::endl;
#endif

    // get particle data
    DRT::Node* particle_1 = next_event->particle_1;
    DRT::Node* particle_2 = next_event->particle_2;

    int lid_1 = particle_1->LID();
    int lid_2 = particle_2->LID();

    static LINALG::Matrix<3,1> pos_1, pos_2, vel_1, vel_2, pos_1_new, pos_2_new, vel_1_new, vel_2_new;
    pos_1 = particleData_[lid_1].dis;
    pos_2 = particleData_[lid_2].dis;
    vel_1 = particleData_[lid_1].vel;
    vel_2 = particleData_[lid_2].vel;

    // compute particle positions and collision normal at collision time
    static LINALG::Matrix<3,1> unitcollnormal;
    pos_1_new.Update(1.,pos_1,next_event->time-particleData_[lid_1].ddt,vel_1);
    pos_2_new.Update(1.,pos_2,next_event->time-particleData_[lid_2].ddt,vel_2);
    unitcollnormal.Update(1.,pos_2_new,-1.,pos_1_new);
    unitcollnormal.Scale(1.0/unitcollnormal.Norm2());

    // compute velocities of particles in normal direction
    const double veln1(vel_1.Dot(unitcollnormal));
    const double veln2(vel_2.Dot(unitcollnormal));

    // check for collision: normal velocity of particle_1 must be greater than of particle_2
    const double deltaveln = veln1 - veln2;
    if (deltaveln > GEO::TOL14)
    {
      // get masses
      const double mass_1 = particleData_[lid_1].mass;
      const double mass_2 = particleData_[lid_2].mass;
      const double invmass = 1.0 / (mass_1 + mass_2);


      // compute new velocities in normal direction
      const double veln1_new = veln1 - (1.0 + e_) * mass_2 * deltaveln * invmass;
      const double veln2_new = veln2 + (1.0 + e_) * mass_1 * deltaveln * invmass;

      // compute new velocities
      vel_1_new.Update(1.,vel_1,veln1_new-veln1,unitcollnormal);
      vel_2_new.Update(1.,vel_2,veln2_new-veln2,unitcollnormal);

      // update ddt
      particleData_[lid_1].ddt = next_event->time;
      particleData_[lid_2].ddt = next_event->time;

      particleData_[lid_1].dis = pos_1_new;
      particleData_[lid_2].dis = pos_2_new;

      particleData_[lid_1].vel = vel_1_new;
      particleData_[lid_2].vel = vel_2_new;

#ifdef OUTPUT
      std::cout << "New position of particle with GID " << particle_1->Id() << "  x: " << pos_1_new(0) << "  y: " << pos_1_new(1) << "  z: " << pos_1_new(2) << std::endl;
      std::cout << "New position of particle with GID " << particle_2->Id() << "  x: " << pos_2_new(0) << "  y: " << pos_2_new(1) << "  z: " << pos_2_new(2) << std::endl;
#endif
    }
  }
  break;
  case INPAR::PARTICLE::particle_wall:
  {
#ifdef OUTPUT
    std::cout << "Handle collision of particle " << next_event->particle_1->Id()
        << " with wall id: " << Teuchos::rcp_static_cast<WallEvent>(next_event)->wall->Id() << std::endl;
#endif

    // collision time
    double colltime = next_event->time;

    static LINALG::Matrix<3,1> initposition;
    static LINALG::Matrix<3,1> initvelocity;
    double radius;
    double particle_time;

    GetCollisionData(next_event->particle_1, initposition, initvelocity, radius, particle_time);

    // advance particle in time to collision time
    static LINALG::Matrix<3,1> newpos;
    newpos.Update(1.0, initposition, colltime - particle_time, initvelocity);

    static LINALG::Matrix<3,1> collnormal;
    collnormal.Update(1.0, newpos, -1.0, Teuchos::rcp_static_cast<WallEvent>(next_event)->wallcollpoint_pos);

    // safety check
    double normallength = collnormal.Norm2();
    if (std::abs(normallength - radius) > GEO::TOL7)
    {
      std::cout << "ran into error :" << std::endl;
      std::cout << "particle pos: " << newpos << std::endl;
      std::cout << "wallcollpoint_pos: " << Teuchos::rcp_static_cast<WallEvent>(next_event)->wallcollpoint_pos << std::endl;
      std::cout << "normallength: " << normallength << std::endl;
      std::cout << "radius: " << radius << std::endl;
      std::cout << "colltime: " << colltime << std::endl;
      std::cout << "particle_time: " << particle_time << std::endl;
      dserror("Particle and wall collision detected but distance does not match radius");
    }

    // normalize colnormal
    collnormal.Scale(1.0 / normallength);
    double veln_particle = collnormal.Dot(initvelocity);
    double veln_wall = Teuchos::rcp_static_cast<WallEvent>(next_event)->wallcollpoint_vel.Dot(collnormal);

    // walls have infinite mass
    static LINALG::Matrix<3,1> newvel;
    newvel.Update(1.0, initvelocity, (1.0 + e_wall_) * (veln_wall - veln_particle), collnormal);

    // write particle data
    int lid = next_event->particle_1->LID();
    particleData_[lid].ddt = next_event->time;
    particleData_[lid].dis = newpos;
    particleData_[lid].vel = newvel;

#ifdef OUTPUT
    std::cout << "New position of particle " << next_event->particle_1->Id() << "  x: " << newpos(0) << "  y: " << newpos(1) << "  z: " << newpos(2) << std::endl;
    std::cout << "New velocity of particle " << next_event->particle_1->Id() << "  x: " << newvel(0) << "  y: " << newvel(1) << "  z: " << newvel(2) << std::endl;
#endif
  }
  break;
  default:
    dserror("collision type does not exist");
  break;
  }

  return;
}


/*----------------------------------------------------------------------*
 | returns time to inter-particle collision in event       ghamm 09/13  |
 *----------------------------------------------------------------------*/
Teuchos::RCP<PARTICLE::Event> PARTICLE::ParticleCollisionHandlerMD::ComputeCollisionWithParticle(
  DRT::Node* particle_1,
  DRT::Node* particle_2,
  const double remaining_dt
  )
{
  Teuchos::RCP<Event> newevent = Teuchos::null;
  // IMPORTANT: First Particle is the colliding particle which carries the current time in the ddt_-vector

  // leave in case of self-contact
  if (particle_1->Id() == particle_2->Id())
    return newevent;

  static LINALG::Matrix<3,1> pos_1, pos_2, vel_1, vel_2;
  double rad_1, rad_2, ddt_1, ddt_2;
  if(not particleData_.empty())
  {
    const ParticleCollData& particle1 = particleData_[particle_1->LID()];
    pos_1.Update(particle1.dis);
    vel_1.Update(particle1.vel);
    rad_1 = particle1.rad;
    ddt_1 = particle1.ddt;
    const ParticleCollData& particle2 = particleData_[particle_2->LID()];
    pos_2.Update(particle2.dis);
    vel_2.Update(particle2.vel);
    rad_2 = particle2.rad;
    ddt_2 = particle2.ddt;
  }
  else
  {
    GetCollisionData(particle_1, particle_2, pos_1, pos_2, vel_1, vel_2, rad_1, rad_2, ddt_1, ddt_2);
  }

  static LINALG::Matrix<3,1> deltax, deltav;

  for (int i=0; i<3; ++i)
  {
    // first particle is the colliding particle thus it carries the actual time in the ddt_-vector
    // additionally the second particle needs to be updated to the time of particle_1
    deltax(i) = pos_1(i) - pos_2(i) - (ddt_1 - ddt_2) * vel_2(i);

    deltav(i) = vel_1(i) - vel_2(i);
  }

  const double sigma = rad_1 + rad_2;

  // finding possible collision times
  const double a = deltav.Dot(deltav);
  const double b = 2.0 * deltav.Dot(deltax);
  const double c = deltax.Dot(deltax) - sigma * sigma;

  const double discriminant = b * b - 4.0 * a * c;

  if (discriminant >= 0.0 and a>0.0)
  {
    const double sqrdiscr = sqrt(discriminant);
    const double inv2a = 1.0 / (2.0 * a);
    const double tc1 = (-b - sqrdiscr) * inv2a;
    const double tc2 = (-b + sqrdiscr) * inv2a;

    // immediate collision of particles expected
    if (std::abs(tc1) <= GEO::TOL14 or std::abs(tc2) <= GEO::TOL14)
    {
      // compute collision normal to detect whether collision has already happened at the end of the last time step
      static LINALG::Matrix<3,1> colnormal;
      colnormal.Update(1.0, pos_2, -1.0, pos_1);

      // velocities of particles in normal direction
      const double vel_col_1 =   vel_1.Dot(colnormal);
      const double vel_col_2 = - vel_2.Dot(colnormal);

      if ((vel_col_1 + vel_col_2) > GEO::TOL14)
      {
        newevent = Teuchos::rcp(new Event(INPAR::PARTICLE::particle_particle, ddt_1+0.0, particle_1, particle_2));
      }
    }
    // tc1 is negative
    else if (tc1 < -GEO::TOL14 and tc2 > GEO::TOL14 and tc2 < 1.1*remaining_dt)
    {
      newevent = Teuchos::rcp(new Event(INPAR::PARTICLE::particle_particle, ddt_1+tc2, particle_1, particle_2));
    }
    // tc2 is negative
    else if (tc1 > GEO::TOL14 and tc2 < -GEO::TOL14 and tc1 < 1.1*remaining_dt)
    {
      newevent = Teuchos::rcp(new Event(INPAR::PARTICLE::particle_particle, ddt_1+tc1, particle_1, particle_2));
    }
    // both positive, smaller one is chosen (tc1 is almost identical to tc2)
    else if (tc1 > GEO::TOL14 and tc2 > GEO::TOL14 and std::min(tc1, tc2) < 1.1*remaining_dt)
    {
      newevent = Teuchos::rcp(new Event(INPAR::PARTICLE::particle_particle, ddt_1+std::min(tc1, tc2), particle_1, particle_2));
    }
  }

  // event carries time to collision plus time of the last event
  return newevent;
}


/*----------------------------------------------------------------------*
 | returns time to particle-wall collision in event        ghamm 09/13  |
 *----------------------------------------------------------------------*/
Teuchos::RCP<PARTICLE::WallEvent> PARTICLE::ParticleCollisionHandlerMD::ComputeCollisionWithWall(
  DRT::Node* particle,
  DRT::Element* wall,
  const double dt
  )
{
  static LINALG::Matrix<3,1> position;
  static LINALG::Matrix<3,1> velocity;
  double radius;
  double particle_time;

  if(not particleData_.empty())
  {
    const ParticleCollData& data1 = particleData_[particle->LID()];
    position.Update(data1.dis);
    velocity.Update(data1.vel);
    radius = data1.rad;
    particle_time = data1.ddt;
  }
  else
  {
    GetCollisionData(particle, position, velocity, radius, particle_time);
  }

  // variables to be filled with collision data
  double timetocollision = 0.0;;
  static LINALG::Matrix<3,1> wallcoll_pos;
  static LINALG::Matrix<3,1> wallcoll_vel;

  // get wall discretization and displacement states
  Teuchos::RCP<DRT::Discretization> walldiscret = particle_algorithm_->WallDiscret();
  Teuchos::RCP<const Epetra_Vector> walldisn = walldiscret->GetState("walldisn");
  Teuchos::RCP<const Epetra_Vector> walldisnp = walldiscret->GetState("walldisnp");

  DRT::Node** nodes = wall->Nodes();
  const int numnodes = wall->NumNode();

  std::vector<int> lm;
  lm.reserve(numnodes * 3);
  std::vector<int> lmowner;
  std::vector<int> lmstride;
  wall->LocationVector(*walldiscret, lm, lmowner, lmstride);

  // initial nodal displacements
  std::vector<double> nodaldisn(numnodes * 3);
  DRT::UTILS::ExtractMyValues(*walldisn, nodaldisn, lm);

  // final nodal displacements
  std::vector<double> nodaldisnp(numnodes * 3);
  DRT::UTILS::ExtractMyValues(*walldisnp, nodaldisnp, lm);

  Epetra_SerialDenseMatrix xyze_n(3,numnodes);
  Epetra_SerialDenseMatrix xyze_np(3,numnodes);

  for (int inode=0; inode<numnodes; ++inode)
  {
    const DRT::Node* node = nodes[inode];
    const double* x_refe = node->X();
    for (int a=0; a<3; ++a)
    {
      xyze_n(a,inode) = x_refe[a] + nodaldisn[3 * inode + a];
      xyze_np(a,inode) = x_refe[a] + nodaldisnp[3 * inode + a];
    }
  }

  // Note: particles can already be advanced to some point in [0,dt)
  // Hence, their position is already updated to that point in time (done in HandleCollision) while wall positions are always
  // interpolated to the respective time using displacements from time n and n+1
  const bool validcollision = ComputeCollisionOfParticleWithWall(wall, xyze_n, xyze_np, position, velocity,
      radius, timetocollision, wallcoll_pos, wallcoll_vel, dt - particle_time, dt);

  // fill particle-wall event in which time to collision is inserted here so that current time needs to be added
  Teuchos::RCP<WallEvent> wallevent = Teuchos::null;
  if(validcollision == true)
  {
    wallevent = Teuchos::rcp(new WallEvent(INPAR::PARTICLE::particle_wall, particle_time+timetocollision, particle, wall, wallcoll_pos, wallcoll_vel));

#ifdef DEBUG
    // safety check
    double colltime = wallevent->time;

    static LINALG::Matrix<3,1> initposition;
    static LINALG::Matrix<3,1> initvelocity;

    double radius;
    double particle_time;

    GetCollisionData(particle, initposition, initvelocity, radius, particle_time);

    // advance particle in time to collision time
    static LINALG::Matrix<3,1> newpos;
    newpos.Update(1.0, initposition, colltime - particle_time, initvelocity);

    static LINALG::Matrix<3,1> collnormal;
    collnormal.Update(1.0, newpos, -1.0, Teuchos::rcp_static_cast<WallEvent>(wallevent)->wallcollpoint_pos);

    double normallength = collnormal.Norm2();
    if (std::abs(normallength - radius) > GEO::TOL7)
    {
      std::cout << "ran into error :" << std::endl;
      std::cout << "particle pos: " << newpos << std::endl;
      std::cout << "wallcollpoint_pos: " << Teuchos::rcp_static_cast<WallEvent>(wallevent)->wallcollpoint_pos << std::endl;
      std::cout << "normallength: " << normallength << std::endl;
      std::cout << "radius: " << radius << std::endl;
      std::cout << "wall id: " << wall->Id() << " and particle lid: " << particle->LID() << " and particle gid: " << particle->Id() << std::endl;
      std::cout << "colltime: " << colltime << std::endl;
      std::cout << "particle_time: " << particle_time << std::endl;
      dserror("Particle and wall collision detected but distance does not match radius");
    }
#endif
  }

  return wallevent;
}


/*----------------------------------------------------------------------*
 | setup of initial event queue at time step begin         ghamm 09/13  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerMD::InitializeEventQueue(
  std::set<Teuchos::RCP<Event>, Event::Helper>& eventqueue,
  const double dt
  )
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleCollisionHandlerMD::InitializingEventqueue");

  // setup of initial event queue
  std::set<int> examinedbins;
  // list of all particles in the neighborhood of currparticle
  std::list<DRT::Node*> neighboring_particles;

  for (int i=0; i<discret_->NodeColMap()->NumMyElements(); ++i)
  {
    // particle for which contact will be detected
    DRT::Node *currparticle = discret_->lColNode(i);

    DRT::Element** currele = currparticle->Elements();
    DRT::Element* currbin = currele[0];
    const int binId = currbin->Id();

    // if a bin has already been examined --> continue with next particle
    if( examinedbins.count(binId) == 1 )
    {
      continue;
    }
    // else: bin is examined for the first time --> new entry in examinedbins
    else
    {
      examinedbins.insert(binId);
    }

    // remove current content but keep memory
    neighboring_particles.clear();
    boost::unordered_map<int, DRT::Element*> neighbouring_walls;
    //std::set<Teuchos::RCP<HeatSource>, BINSTRATEGY::Less> neighboring_heatSources;

    // gather all neighboring particles and wall elements
    particle_algorithm_->GetNeighbouringItems(currparticle, neighboring_particles, neighbouring_walls);

    DRT::Node** particles = currbin->Nodes();
    for(int iparticle=0; iparticle<currbin->NumNode(); ++iparticle)
    {
      DRT::Node* currnode = particles[iparticle];

      // loop over all neighboring particles and check if they collide with currnode
      for (std::list<DRT::Node*>::iterator iter=neighboring_particles.begin(); iter!=neighboring_particles.end(); ++iter)
      {
        // in order to avoid double evaluation of the same contact, only search for contact when id_1 < id_2
        if (currnode->Id() > (*iter)->Id())
        {
          continue; // with next particle in neighborhood
        }

        Teuchos::RCP<PARTICLE::Event> newevent = ComputeCollisionWithParticle(currnode, *iter, dt);

        // insert event into event queue if collision is valid
        if (newevent != Teuchos::null)
        {
#ifdef OUTPUT
          std::cout << "inserting inter-particle collision in the event queue" << std::endl;
#endif
          eventqueue.insert(newevent);
        }
      }

      // loop over all neighbouring wall elements and check if they collide with currparticle
      for (boost::unordered_map<int, DRT::Element*>::iterator iter=neighbouring_walls.begin(); iter!=neighbouring_walls.end(); ++iter)
      {
        Teuchos::RCP<WallEvent> newevent = ComputeCollisionWithWall(currnode, iter->second, dt);

        // insert event into event queue if collision is valid
        if (newevent != Teuchos::null)
        {
#ifdef OUTPUT
          std::cout << "inserting particle-wall collision in the event queue" << std::endl;
#endif
          // here we are at the beginning of the time step so that event.time does not need an update
          eventqueue.insert(newevent);
        }
      }
    }
  }

  // print event queue
#ifdef OUTPUT
  std::cout << "The eventqueue contains after initializing: ";
  int count = 0;
  for (std::set<Teuchos::RCP<Event>, Event::Helper>::iterator iter2 = eventqueue.begin(); iter2 != eventqueue.end(); ++iter2)
  {
    std::cout << "This is entry number " << count << " of the eventqueue:" << std::endl;

    Teuchos::RCP<Event> output2 = *iter2;
    std::cout << " happening at time:  " << output2->time << std::endl;
    std::cout << "GID of first particle " << output2->particle_1->Id() << std::endl;
    if (output2->coltype == INPAR::PARTICLE::particle_particle)
    {
      std::cout << "GID of second particle " << output2->particle_2->Id() << std::endl;
    }
    else if (output2->coltype == INPAR::PARTICLE::particle_wall)
    {
      std::cout << "GID of wall " << Teuchos::rcp_dynamic_cast<WallEvent>(output2)->wall->Id() << std::endl;
    }
    count++;
  }
  std::cout << std::endl;
#endif

  return;
}


/*----------------------------------------------------------------------*
 | event queue is filled with new events based on          ghamm 09/13  |
 | particles that have just collided                                    |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerMD::SearchForNewCollisions(
  Teuchos::RCP<Event> lastevent,
  std::set<Teuchos::RCP<Event>, Event::Helper>& eventqueue,
  const double dt
  )
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::ParticleCollisionHandlerMD::UpdatdingEventqueue");

  std::list<DRT::Node*> neighbouring_particles;
  boost::unordered_map<int, DRT::Element*> neighbouring_walls;
  //std::set<Teuchos::RCP<HeatSource>, BINSTRATEGY::Less> neighboring_heatSources;

  //
  // searching for new collisions for the first particle of the last collision
  //
  particle_algorithm_->GetNeighbouringItems(lastevent->particle_1, neighbouring_particles, neighbouring_walls);

  // particle-particle collision
  for (std::list<DRT::Node*>::iterator iter=neighbouring_particles.begin(); iter!=neighbouring_particles.end(); ++iter)
  {
    // IMPORTANT: first particle must be the particle, that collided in this time step
    Teuchos::RCP<Event> newevent = ComputeCollisionWithParticle(lastevent->particle_1, *iter, dt-lastevent->time);

    if(newevent != Teuchos::null)
    {
      if (lastevent->coltype == INPAR::PARTICLE::particle_particle)
      {
        // do not add event of particles that has already been processed in lastevent
        if (newevent->particle_2->Id() != lastevent->particle_2->Id())
        {
          eventqueue.insert(newevent);
        }
      }
      else
      {
        eventqueue.insert(newevent);
      }
    }
  }

  // particle-wall collision
  for(boost::unordered_map<int, DRT::Element*>::iterator iter=neighbouring_walls.begin(); iter!=neighbouring_walls.end(); ++iter)
  {
    Teuchos::RCP<WallEvent> newevent = ComputeCollisionWithWall(lastevent->particle_1, iter->second, dt);

    if (newevent != Teuchos::null)
      eventqueue.insert(newevent);
  }

  //
  // searching for new collisions for the second particle of the last collision (if existing)
  //
  if (lastevent->coltype == INPAR::PARTICLE::particle_particle)
  {
    // reuse neighborhood in case both particles reside in the same bin
    const int binId_part1 = lastevent->particle_1->Elements()[0]->Id();
    const int binId_part2 = lastevent->particle_2->Elements()[0]->Id();
    if( binId_part1 != binId_part2)
      particle_algorithm_->GetNeighbouringItems(lastevent->particle_2, neighbouring_particles, neighbouring_walls);

    // particle-particle collision
    for(std::list<DRT::Node*>::iterator iter=neighbouring_particles.begin(); iter!=neighbouring_particles.end(); ++iter)
    {
      // IMPORTANT: first particle must be the particle, that collided in this time step
      Teuchos::RCP<Event> newevent = ComputeCollisionWithParticle(lastevent->particle_2, *iter, dt-lastevent->time);

      // do not add event of particles that has already been processed in lastevent
      if (newevent != Teuchos::null and newevent->particle_2->Id() != lastevent->particle_1->Id())
      {
        eventqueue.insert(newevent);
      }
    }

    // particle-wall collision
    for(boost::unordered_map<int, DRT::Element*>::iterator iter=neighbouring_walls.begin(); iter!=neighbouring_walls.end(); ++iter)
    {
      Teuchos::RCP<WallEvent> newevent = ComputeCollisionWithWall(lastevent->particle_2, iter->second, dt);

      if (newevent != Teuchos::null)
      {
        eventqueue.insert(newevent);
      }
    }
  }

#ifdef OUTPUT
  std::cout << "The eventqueue contains after searching for new events: ";
  int count = 0;
  for (std::set<Teuchos::RCP<Event>, Event::Helper>::iterator iter2 = eventqueue.begin(); iter2 != eventqueue.end(); ++iter2)
  {
    std::cout << "This is entry number " << count << " of the eventqueue:" << std::endl;

    Teuchos::RCP<Event> output2 = *iter2;
    std::cout << " happening at time:  " << output2->time << std::endl;
    std::cout << "GID of first particle " << output2->particle_1->Id() << std::endl;
    if (output2->coltype == INPAR::PARTICLE::particle_particle)
    {
      std::cout << "GID of second particle " << output2->particle_2->Id() << std::endl;
    }
    else if (output2->coltype == INPAR::PARTICLE::particle_wall)
    {
      std::cout << "GID of wall " << Teuchos::rcp_dynamic_cast<WallEvent>(output2)->wall->Id() << std::endl;
    }
    count++;
  }
  std::cout << std::endl;
#endif

  return;
}


/*----------------------------------------------------------------------*
 | computes time to collision between a particle and a     ghamm 05/14  |
 | wall: searches hierarchically element, edges, corners                |
 *----------------------------------------------------------------------*/
bool PARTICLE::ParticleCollisionHandlerMD::ComputeCollisionOfParticleWithWall(
    DRT::Element* wallele,
    const Epetra_SerialDenseMatrix& xyze_n,
    const Epetra_SerialDenseMatrix& xyze_np,
    const LINALG::Matrix<3,1>& particle_pos,
    const LINALG::Matrix<3,1>& particle_vel,
    const double radius,
    double& timetocollision,
    LINALG::Matrix<3,1>& wall_pos,
    LINALG::Matrix<3,1>& wall_vel,
    const double remaining_dt,
    const double dt)
{
  bool checkedges = false;
  bool checkcorners = false;

  // check collision with element itself first
  bool validcollision = ComputeCollisionOfParticleWithElement(wallele, xyze_n, xyze_np, particle_pos, particle_vel,
      radius, timetocollision, wall_pos, wall_vel, remaining_dt, dt, checkedges);

  // if necessary, check for collision of particle with edges of element
  if(checkedges == true)
  {
    DRT::Element::DiscretizationType eleshape = wallele->Shape();
    // some variables to store data to find closest collision point with edges of wall element
    double timetocollision_line = GEO::LARGENUMBER;;
    static LINALG::Matrix<3,1> wallcollpoint_pos_line;
    static LINALG::Matrix<3,1> wallcollpoint_vel_line;
    // set time to collision to a large number and search for edge collision points which are closer than this number
    timetocollision = GEO::LARGENUMBER;

    // run over all line elements
    const std::vector<Teuchos::RCP<DRT::Element> > eleLines = wallele->Lines();

    const int numnodes_line = eleLines[0]->NumNode();
    Epetra_SerialDenseMatrix xyze_line_n(3,numnodes_line);
    Epetra_SerialDenseMatrix xyze_line_np(3,numnodes_line);

    bool validlinecollision = false;
    bool checkcorners_iter = false;
    for(int i=0; i<wallele->NumLine(); ++i)
    {
      for(int inode=0; inode<numnodes_line; ++inode)
      {
        switch(eleshape)
        {
        case DRT::Element::tri3:
        case DRT::Element::tri6:
        {
          for (int a=0; a<3; ++a)
          {
            xyze_line_n(a,inode) = xyze_n(a,DRT::UTILS::eleNodeNumbering_tri6_lines[i][inode]);
            xyze_line_np(a,inode) = xyze_np(a,DRT::UTILS::eleNodeNumbering_tri6_lines[i][inode]);
          }
          break;
        }
        case DRT::Element::quad4:
        case DRT::Element::quad8:
        case DRT::Element::quad9:
        {
          for (int a=0; a<3; ++a)
          {
            xyze_line_n(a,inode) = xyze_n(a,DRT::UTILS::eleNodeNumbering_quad9_lines[i][inode]);
            xyze_line_np(a,inode) = xyze_np(a,DRT::UTILS::eleNodeNumbering_quad9_lines[i][inode]);
          }
          break;
        }
        default: dserror("intface type not supported %d", wallele->Shape());
          break;
        }
      }

      // find possible collision point for this line
      validlinecollision = ComputeCollisionOfParticleWithLine(eleLines[i]->Shape(), xyze_line_n, xyze_line_np, particle_pos, particle_vel,
          radius, timetocollision_line, wallcollpoint_pos_line, wallcollpoint_vel_line, remaining_dt, dt, checkcorners_iter);

      // find closest valid edge contact point if more than one exists
      // check whether collision time is reasonable
      if (validlinecollision && timetocollision_line < timetocollision)
      {
        validcollision = true;
        timetocollision = timetocollision_line;
        wall_pos.Update(wallcollpoint_pos_line);
        wall_vel.Update(wallcollpoint_vel_line);
#ifdef OUTPUT
        std::cout << "found the following valid line pos: " <<wall_pos(0) << " " << wall_pos(1) << " " << wall_pos(2) << std::endl;
        std::cout << "found the following valid line time to coll: " << timetocollision << std::endl;
#endif
      }
      if(checkcorners_iter == true)
        checkcorners = true;
    }

    // if something was found, corners do not need to be checked
    if(validcollision == true)
    {
      return validcollision;
    }
  }

  // if necessary, check for collision of particle with corners of element
  if(checkcorners == true)
  {
    // some variables to store data to find closest collision point with edges of wall element
    double timetocollision_corner = GEO::LARGENUMBER;
    static LINALG::Matrix<3,1> wallcollpoint_pos_corner;
    static LINALG::Matrix<3,1> wallcollpoint_vel_corner;

    // set time to collision to a large number and search for edge collision points which are closer than this number
    timetocollision = GEO::LARGENUMBER;

    // loop all corner nodes
    bool validcornercollision = false;
    for(int inode=0; inode<DRT::UTILS::getNumberOfElementCornerNodes(wallele->Shape()); ++inode)
    {
      Epetra_SerialDenseMatrix xyze_corner_n(3,1);
      Epetra_SerialDenseMatrix xyze_corner_np(3,1);

      for (int a=0; a<3; ++a)
      {
        xyze_corner_n(a,0) = xyze_n(a,inode);
        xyze_corner_np(a,0) = xyze_np(a,inode);
      }

      // find possible collision point for this corner
      validcornercollision = ComputeCollisionOfParticleWithCorner(xyze_corner_n, xyze_corner_np, particle_pos, particle_vel,
          radius, timetocollision_corner, wallcollpoint_pos_corner, wallcollpoint_vel_corner, remaining_dt, dt);

      // find closest valid corner contact point if more than one exists
      // check whether collision time is reasonable
      if (validcornercollision && timetocollision_corner < timetocollision)
      {
        validcollision = true;
        timetocollision = timetocollision_corner;
        wall_pos.Update(wallcollpoint_pos_corner);
        wall_vel.Update(wallcollpoint_vel_corner);
      }
    }
  }

  return validcollision;
}

/*----------------------------------------------------------------------*
 | computes time to collision between a particle  and      ghamm 04/14  |
 | an element for hard sphere particles (templated on distype of ele)   |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType DISTYPE>
bool PARTICLE::ParticleCollisionHandlerMD::ComputeCollisionOfParticleWithElementT_FAD(
    DRT::Element* wallele,
    const Epetra_SerialDenseMatrix& xyze_n,
    const Epetra_SerialDenseMatrix& xyze_final,
    const LINALG::Matrix<3,1>& particle_pos,
    const LINALG::Matrix<3,1>& particle_vel,
    const double radius,
    double& timetocollision,
    LINALG::Matrix<3,1>& wallcollpoint_pos,
    LINALG::Matrix<3,1>& wallcollpoint_vel,
    const double remaining_dt,
    const double dt,
    bool& checkedges)
{
  bool validcollision = false;

  const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // solution vector
  static LINALG::TMatrix<FAD,3,1> coll_solution;

  static LINALG::TMatrix<FAD,3,1> wallcollpoint_pos_fad;
  static LINALG::TMatrix<FAD,3,1> wallcollpoint_vel_fad;

  static LINALG::TMatrix<FAD,3,1> particle_pos_fad;
  static LINALG::TMatrix<FAD,3,1> particle_vel_fad;
  for(int i=0; i<3; ++i)
  {
    particle_pos_fad(i) = particle_pos(i);
    particle_vel_fad(i) = particle_vel(i);
  }

  // initial guess for wall collision point
  GEO::startingValueCurrentToElementCoords<DISTYPE>(coll_solution);
  // starting time is zero
  coll_solution(2) = 0.0;

  // setup FAD
  for(int i=0; i<3; ++i)
    coll_solution(i).diff(i,3);

  // unit normal at collision point
  static LINALG::TMatrix<FAD,3,1> unitnormal;

  // velocity of wall element is constant over the time step
  Epetra_SerialDenseMatrix vele(xyze_n);
  vele.Scale(-1.0);
  vele += xyze_final;
  vele.Scale(1.0 / dt);

  static LINALG::TMatrix<FAD,3,numnode> vele_fad;
  static LINALG::TMatrix<FAD,3,numnode> xyze_n_fad;
  static LINALG::TMatrix<FAD,3,numnode> xyze_colltime;
  for(int i=0; i<3; ++i)
  {
    for(int j=0; j<numnode; ++j)
    {
      vele_fad(i,j) = vele(i,j);
      xyze_n_fad(i,j) = xyze_n(i,j);
    }
  }

  // the following equation is solved iteratively w.r.t. ele coords and time to collision (ttc):
  // particle_pos + ttc * particle_vel + UnitNormalAtWallCollPoint(t)*radius = WallCollPoint(t)
  // with t = dt - remaining_dt + ttc

  // iteration for contact search
  const int maxiter = 10;
  int iter = 0;
  while (iter < maxiter)
  {
    ++iter;

    // compute wall positions at collision time: xyze_colltime = xyze_n + colltime*vele
    FAD colltime = dt-remaining_dt+coll_solution(2);
    xyze_colltime.Update(1.0, xyze_n_fad, colltime, vele_fad);

    // position and velocity of wall collision point at collision time
    static LINALG::TMatrix<FAD,numnode,1> funct;
    DRT::UTILS::shape_function<DISTYPE>(coll_solution, funct);
    wallcollpoint_pos_fad.Clear();
    wallcollpoint_vel_fad.Clear();
    for(int i=0; i<numnode; ++i)
      for(int j=0; j<3; ++j)
      {
        wallcollpoint_pos_fad(j) += xyze_colltime(j,i) * funct(i);
        wallcollpoint_vel_fad(j) += vele_fad(j,i) * funct(i);
      }

    // get unit normal of wall element at collision point
    {
      static LINALG::TMatrix<FAD,2,numnode> deriv;
      DRT::UTILS::shape_function_2D_deriv1(deriv,coll_solution(0),coll_solution(1),DISTYPE);

      // compute dXYZ / drs
      static LINALG::TMatrix<FAD,3,2> dxyzdrs;
      dxyzdrs.Clear();
      for (int i=0; i<3; ++i)
        for (int j=0; j<2; ++j)
          for (int k=0; k<numnode; ++k)
            dxyzdrs(i,j) += xyze_colltime(i,k)*deriv(j,k);

      // compute normal at collision point
      unitnormal(0) = dxyzdrs(1,0) * dxyzdrs(2,1) - dxyzdrs(2,0) * dxyzdrs(1,1);
      unitnormal(1) = dxyzdrs(2,0) * dxyzdrs(0,1) - dxyzdrs(0,0) * dxyzdrs(2,1);
      unitnormal(2) = dxyzdrs(0,0) * dxyzdrs(1,1) - dxyzdrs(1,0) * dxyzdrs(0,1);

      FAD norm = unitnormal(0)*unitnormal(0) + unitnormal(1)*unitnormal(1) + unitnormal(2)*unitnormal(2);

      // normalize
      unitnormal.Scale(1.0 / std::pow(norm,0.5));
    }

    // update particle position to current time
    static LINALG::TMatrix<FAD,3,1> particle_pos_colltime;
    particle_pos_colltime.Update(1.0, particle_pos_fad, coll_solution(2), particle_vel_fad);

    // check whether normal is pointing outward (particle is inside) and adapt it if necessary
    static LINALG::TMatrix<FAD,3,1> testvector_fad;
    testvector_fad.Update(1.0, wallcollpoint_pos_fad, -1.0, particle_pos_colltime);

    if (unitnormal.Dot(testvector_fad) < 0.0)
      unitnormal.Scale(-1.0);

    // compute residual
    static LINALG::TMatrix<FAD,3,1> residual;
    double bnorm = 0.0;
    for (int i=0; i<3; ++i)
    {
      residual(i) = particle_pos_colltime(i)  + unitnormal(i) * radius - wallcollpoint_pos_fad(i);
      bnorm += residual(i).val() * residual(i).val();
    }

    if (std::sqrt(bnorm) < GEO::TOL14)
    {
      break;
    }

    // compute dF/dx
    static LINALG::Matrix<3,3> A;

    for (int j=0; j<3; ++j)
    {
      A(0, j) = residual(0).dx(j);
      A(1, j) = residual(1).dx(j);
      A(2, j) = residual(2).dx(j);
    }

    // compute rhs
    static LINALG::Matrix<3,1> b;
    for (int i=0; i<3; ++i)
    {
      b(i) = - residual(i).val();
    }

    // solve linear problem A dx = b
    static LINALG::Matrix<3,1> dx;
    const double det = LINALG::scaledGaussElimination<3>(A, b, dx);

    if (std::abs(det) < GEO::TOL14)
    {
      if(unitnormal.Dot(particle_vel_fad) < GEO::TOL10 and iter>1)
      {
#ifdef OUTPUT
        std::cout << "particle path is parallel to wall --> left iteration" << std::endl;
#endif
        return validcollision;
      }
    }

    // update of coll_solution
    for (int i=0; i<3; ++i)
    {
      coll_solution(i) += dx(i);
    }

    // in case of parallel movement of particle to wall, element coords get extremely large
    if (std::abs(coll_solution(0)) > 1.0e3 or std::abs(coll_solution(1)) > 1.0e3)
    {
#ifdef OUTPUT
      std::cout << "elecoord of wall is extremely large --> left iteration" << std::endl;
#endif
      return validcollision;
    }
  }

  FAD scalar_partvel = unitnormal.Dot(particle_vel_fad);
  FAD scalar_wallvel = unitnormal.Dot(wallcollpoint_vel_fad);

  // check whether collision position is valid and otherwise check for edge collision points
  if (GEO::checkPositionWithinElementParameterSpace(coll_solution, DISTYPE) == true)
  {
    // check whether collision time is reasonable
    if (coll_solution(2) >= -GEO::TOL14 and coll_solution(2) < 1.1 * remaining_dt)
    {
      // decide if collision is still going to happen (--> valid) or has already happened in the last time step (--> invalid)
      if ((scalar_partvel - scalar_wallvel) > GEO::TOL14)
      {
        validcollision = true;
        timetocollision = coll_solution(2).val();
        for(int j=0; j<3; ++j)
        {
          wallcollpoint_pos(j) = wallcollpoint_pos_fad(j).val();
          wallcollpoint_vel(j) = wallcollpoint_vel_fad(j).val();
        }
#ifdef OUTPUT
        std::cout << "valid wall coll: time to collision: " << timetocollision << std::endl;
        std::cout << "wallcollpoint_pos: " << wallcollpoint_pos << std::endl;
#endif
      }
    }
  }
  else
  {
    // pre computations
    static LINALG::Matrix<3,1> particevel_cpy;
    static LINALG::Matrix<3,1> wallvel_cpy;
    for(int j=0; j<3; ++j)
    {
      particevel_cpy(j) = particle_vel_fad(j).val();
      wallvel_cpy(j) = wallcollpoint_vel_fad(j).val();
    }

    // particle path points in element direction
    // slightly flying away is still treated as possible edge contact scenario
    if ((scalar_partvel - scalar_wallvel) > -0.1*(particevel_cpy.Norm1()+wallvel_cpy.Norm1()))
    {
      checkedges = true;
    }
    else
    {
#ifdef OUTPUT
      std::cout << "particle flys away from element -> do not check edges and corners" << std::endl;
#endif
#ifdef DEBUG
      // check edges
      checkedges = true;
#endif
    }
  }

  return validcollision;
}


/*----------------------------------------------------------------------*
 | computes time to collision between a particle and       ghamm 09/13  |
 | an element for hard sphere particles                                 |
 *----------------------------------------------------------------------*/
bool PARTICLE::ParticleCollisionHandlerMD::ComputeCollisionOfParticleWithElement(
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
    const double dt,
    bool& checkedges)
{
  bool validcollision = false;
  switch (wallele->Shape())
  {
  case DRT::Element::quad4:
    validcollision = ComputeCollisionOfParticleWithElementT_FAD<DRT::Element::quad4>(
        wallele, xyze_current, xyze_final, particle_pos, particle_vel, radius,
        timetocollision, wall_pos, wall_vel, remaining_dt, dt, checkedges);
    break;
  case DRT::Element::quad8:
    validcollision = ComputeCollisionOfParticleWithElementT_FAD<DRT::Element::quad8>(
        wallele, xyze_current, xyze_final, particle_pos, particle_vel, radius,
        timetocollision, wall_pos, wall_vel, remaining_dt, dt, checkedges);
    break;
  case DRT::Element::quad9:
    validcollision = ComputeCollisionOfParticleWithElementT_FAD<DRT::Element::quad9>(
        wallele, xyze_current, xyze_final, particle_pos, particle_vel, radius,
        timetocollision, wall_pos, wall_vel, remaining_dt, dt, checkedges);
    break;
  case DRT::Element::tri3:
    validcollision = ComputeCollisionOfParticleWithElementT_FAD<DRT::Element::tri3>(
        wallele, xyze_current, xyze_final, particle_pos, particle_vel, radius,
        timetocollision, wall_pos, wall_vel, remaining_dt, dt, checkedges);
    break;
  case DRT::Element::tri6:
    validcollision = ComputeCollisionOfParticleWithElementT_FAD<DRT::Element::tri6>(
        wallele, xyze_current, xyze_final, particle_pos, particle_vel, radius,
        timetocollision, wall_pos, wall_vel, remaining_dt, dt, checkedges);
    break;
  default:
    dserror("please add your surface element type here");
    break;
  }

  return validcollision;
}


/*----------------------------------------------------------------------*
 | computes time to collision between a particle and       ghamm 04/14  |
 | an element edge for hard sphere particles (templated on distype      |
 | of wall element edge)                                                |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType DISTYPE>
bool PARTICLE::ParticleCollisionHandlerMD::ComputeCollisionOfParticleWithLineT(
    const Epetra_SerialDenseMatrix& xyze_line_n,
    const Epetra_SerialDenseMatrix& xyze_line_np,
    const LINALG::Matrix<3,1>& particle_pos,
    const LINALG::Matrix<3,1>& particle_vel,
    const double radius,
    double& timetocollision,
    LINALG::Matrix<3,1>& wallcollpoint_pos,
    LINALG::Matrix<3,1>& wallcollpoint_vel,
    bool& checkcorners)
{
  // TODO:
  // TODO: this method only works for fixed walls, terms with time derivatives of wall position are missing
  // TODO:
  bool validlinecollision = false;

  const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // coll_solution contains: r and time to collision (r is element coord)
  static LINALG::Matrix<2,1> coll_solution;

  // initial guess for edge collision point
  GEO::startingValueCurrentToElementCoords<DISTYPE>(coll_solution);
  // starting time is zero
  coll_solution(1) = 0.0;

  // connection vector between edge coll point and particle position
  static LINALG::Matrix<3,1> F;
  // compute first derivative of r
  static LINALG::Matrix<3,1> F_deriv1;
  // compute second derivative of r
  static LINALG::Matrix<3,1> F_deriv2;

  double distance = 0.0;

  // the following functional is minimized iteratively (variables are ele coord and time to collision (ttc)):
  // [ { EdgeCollPoint(t) - (particle_pos + ttc * particle_vel) }^2 - radius^2 ]^2 +
  //     0.5 * ( [EdgeCollPoint(t)- (particle_pos + ttc * particle_vel)] dot EdgeCollPoint_deriv(t) )^2
  // with t = dt - remaining_dt + ttc

  // iteration for contact search
  const int maxiter = 20;
  int iter = 0;
  while (iter < maxiter)
  {
    ++iter;

    // update particle position to current time
    static LINALG::Matrix<3,1> particle_pos_t;
    particle_pos_t.Update(1.0, particle_pos, coll_solution(1), particle_vel);

    // determine shapefunction, 1. and 2. derivative at current solution
    static LINALG::Matrix<numnodes,1> funct;
    DRT::UTILS::shape_function_1D(funct, coll_solution(0), DISTYPE);

    static LINALG::Matrix<1,numnodes> deriv1;
    DRT::UTILS::shape_function_1D_deriv1(deriv1, coll_solution(0), DISTYPE);

    static LINALG::Matrix<1,numnodes> deriv2;
    DRT::UTILS::shape_function_1D_deriv2(deriv2, coll_solution(0), DISTYPE);

    // compute linear system

    wallcollpoint_pos.Clear();
    F_deriv1.Clear();
    F_deriv2.Clear();

    for(int i=0; i<3; ++i)
      for(int inode=0; inode<numnodes; ++inode)
      {
        wallcollpoint_pos(i) += xyze_line_np(i,inode) * funct(inode);
        F_deriv1(i)          += xyze_line_np(i,inode) * deriv1(0,inode);
        F_deriv2(i)          += xyze_line_np(i,inode) * deriv2(0,inode);
      }

    // subtract current particle position from the edge collision point
    F.Update(1.0, wallcollpoint_pos, -1.0, particle_pos_t);

    // compute rhs
    static LINALG::Matrix<2,1> b_line;

    double dotprod_pos_pos = 0.0;
    double dotprod_pos_deriv1 = 0.0;
    double dotprod_deriv1_deriv1 = 0.0;
    double dotprod_pos_deriv2 = 0.0;
    double dotprod_pos_v = 0.0;
    double dotprod_v_deriv1 = 0.0;
    for(int i=0; i<3; ++i)
    {
      dotprod_pos_pos    += F(i)*F(i);
      dotprod_pos_deriv1 += F(i)*F_deriv1(i);
      dotprod_deriv1_deriv1 += F_deriv1(i)*F_deriv1(i);
      dotprod_pos_deriv2 += F(i)*F_deriv2(i);
      dotprod_pos_v      += - F(i)*particle_vel(i);
      dotprod_v_deriv1   += - particle_vel(i)*F_deriv1(i);
    }

    const double pos2subtractrad2 = dotprod_pos_pos - radius*radius;

    b_line(0) = pos2subtractrad2*2.0*dotprod_pos_deriv1 +
                  dotprod_pos_deriv1*(dotprod_deriv1_deriv1 + dotprod_pos_deriv2);

    b_line(1) = pos2subtractrad2*2.0*dotprod_pos_v +
                  dotprod_pos_deriv1*dotprod_v_deriv1;

    // distance between particle and edge
    distance = F.Norm2() - radius;

    // rhs is negative residual
    b_line.Scale(-1.0);

    if (b_line.Norm2() < GEO::TOL14)
    {
      break;
    }


    // determine system matrix A_line
    // compute dF/dx
    static LINALG::Matrix<2,2> A_line;

    double dotprod_deriv1_deriv2 = 0.0;
    double dotprod_v_deriv2 = 0.0;
    double dotprod_v_v = 0.0;
    for(int i=0; i<3; ++i)
    {
      dotprod_deriv1_deriv2 += F_deriv1(i)*F_deriv2(i);
      dotprod_v_deriv2   += - particle_vel(i)*F_deriv2(i);
      dotprod_v_v += particle_vel(i)*particle_vel(i);
    }

    // A_line(0,0)
    A_line(0,0) = 2.0*dotprod_pos_deriv1*2.0*dotprod_pos_deriv1 +
                    pos2subtractrad2*2.0*(dotprod_pos_deriv2 + dotprod_deriv1_deriv1) +
                    (dotprod_deriv1_deriv1 + dotprod_pos_deriv2)*(dotprod_deriv1_deriv1 + dotprod_pos_deriv2) +
                    dotprod_pos_deriv1*3.0*dotprod_deriv1_deriv2;

    // A_line(0,1)
    A_line(0,1) = 2.0*dotprod_pos_v*2.0*dotprod_pos_deriv1 +
                    pos2subtractrad2*2.0*dotprod_v_deriv1 +
                    dotprod_v_deriv1*(dotprod_deriv1_deriv1 + dotprod_pos_deriv2) +
                    dotprod_pos_deriv1*dotprod_v_deriv2;

    // A_line(1,0)
    A_line(1,0) = A_line(0,1);

    // A_line(1,1)
    A_line(1,1) = 2.0*dotprod_pos_v*2.0*dotprod_pos_v +
                    pos2subtractrad2*2.0*dotprod_v_v +
                    dotprod_v_deriv1*dotprod_v_deriv1;


    // solve linear problem A dx = b
    static LINALG::Matrix<2,1> dx_line;
    const double det = LINALG::scaledGaussElimination<2>(A_line, b_line, dx_line);

    if (std::abs(det) < GEO::TOL14)
    {
      static LINALG::Matrix<3,1> crossproduct;
      crossproduct.CrossProduct(particle_vel, F_deriv1);
      if(crossproduct.Norm1() < GEO::TOL10 and iter>1)
      {
#ifdef OUTPUT
        std::cout << "particle path is parallel to line --> left iteration" << std::endl;
#endif
        // leave here
        return validlinecollision;
      }
    }

    // update of coll_solution
    coll_solution.Update(1.0, dx_line, 1.0);
  }

  // check whether collision position is valid and otherwise check for corner collision points
  if (GEO::checkPositionWithinElementParameterSpace(coll_solution, DISTYPE) == true and distance/radius < GEO::TOL7)
  {
    // check whether collision time is reasonable
    if (coll_solution(1) >= -GEO::TOL14)
    {
      const double scalar_partvel = F.Dot(particle_vel);
      const double scalar_wallvel = F.Dot(wallcollpoint_vel);

      // decide if collision is still going to happen (--> valid) or has already happened in the last time step (--> invalid)
      if ((scalar_partvel - scalar_wallvel) > GEO::TOL14)
      {
        timetocollision = coll_solution(1);
        validlinecollision = true;
      }
    }
  }
  else
  {
    checkcorners = true;
  }

  return validlinecollision;
}


/*----------------------------------------------------------------------*
 | computes time to collision between a particle and       ghamm 04/14  |
 | an element edge for hard sphere particles (templated on distype      |
 | of wall element edge)                                                |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType DISTYPE>
bool PARTICLE::ParticleCollisionHandlerMD::ComputeCollisionOfParticleWithLineT_FAD(
    const Epetra_SerialDenseMatrix& xyze_line_n,
    const Epetra_SerialDenseMatrix& xyze_line_np,
    const LINALG::Matrix<3,1>& particle_pos,
    const LINALG::Matrix<3,1>& particle_vel,
    const double radius,
    double& timetocollision,
    LINALG::Matrix<3,1>& wallcollpoint_pos,
    LINALG::Matrix<3,1>& wallcollpoint_vel,
    const double remaining_dt,
    const double dt,
    bool& checkcorners)
{
#ifdef DEBUG
  bool heuristic_break = false;
#endif
  bool validlinecollision = false;

  const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // solution vector
  static LINALG::TMatrix<FAD,2,1> coll_solution;

  static LINALG::TMatrix<FAD,3,1> wallcollpoint_pos_fad;
  static LINALG::TMatrix<FAD,3,1> wallcollpoint_vel_fad;
  static LINALG::TMatrix<FAD,3,1> wallcollpoint_deriv_vel_fad;

  static LINALG::TMatrix<FAD,3,1> particle_pos_fad;
  static LINALG::TMatrix<FAD,3,1> particle_vel_fad;
  for(int i=0; i<3; ++i)
  {
    particle_pos_fad(i) = particle_pos(i);
    particle_vel_fad(i) = particle_vel(i);
  }

  // initial guess for wall collision point
  GEO::startingValueCurrentToElementCoords<DISTYPE>(coll_solution);
  // starting time is zero
  coll_solution(1) = 0.0;

  // setup FAD
  for(int i=0; i<2; ++i)
    coll_solution(i).diff(i,2);

  // unit normal at collision point
  static LINALG::TMatrix<FAD,3,1> unitnormal;

  // velocity of wall element is constant over the time step
  static LINALG::TMatrix<FAD,3,numnodes> vele_fad;
  static LINALG::TMatrix<FAD,3,numnodes> xyze_n_fad;
  static LINALG::TMatrix<FAD,3,numnodes> xyze_colltime;
  const double invdt = 1.0 / dt;
  for(int i=0; i<3; ++i)
  {
    for(int j=0; j<numnodes; ++j)
    {
      vele_fad(i,j) = (xyze_line_np(i,j) - xyze_line_n(i,j))*invdt;
      xyze_n_fad(i,j) = xyze_line_n(i,j);
    }
  }

  // connection vector between edge coll point and particle position
  static LINALG::TMatrix<FAD,3,1> F;
  // compute first derivative of r
  static LINALG::TMatrix<FAD,3,1> F_deriv1;
  // compute second derivative of r
  static LINALG::TMatrix<FAD,3,1> F_deriv2;

  double distance = 0.0;
  double ttc_old = GEO::LARGENUMBER;

  // the following functional is minimized iteratively (variables are ele coord and time to collision (ttc)):
  // [ { EdgeCollPoint(t) - (particle_pos + ttc * particle_vel) }^2 - radius^2 ]^2 +
  //     0.5 * ( [EdgeCollPoint(t)- (particle_pos + ttc * particle_vel)] dot EdgeCollPoint_deriv(t) )^2
  // with t = dt - remaining_dt + ttc

  // iteration for contact search
  const int maxiter = 20;
  int iter = 0;
  while (iter < maxiter)
  {
    ++iter;

    // compute wall positions at collision time: xyze_colltime = xyze_n + colltime*vele
    static FAD colltime;
    colltime = dt-remaining_dt+coll_solution(1);
    xyze_colltime.Update(1.0, xyze_n_fad, colltime, vele_fad);

    // update particle position to current time
    static LINALG::TMatrix<FAD,3,1> particle_pos_colltime;
    particle_pos_colltime.Update(1.0, particle_pos_fad, coll_solution(1), particle_vel_fad);

    // determine shape function, 1. and 2. derivative at current solution
    static LINALG::TMatrix<FAD,numnodes,1> funct;
    DRT::UTILS::shape_function_1D(funct, coll_solution(0), DISTYPE);

    static LINALG::TMatrix<FAD,1,numnodes> deriv1;
    DRT::UTILS::shape_function_1D_deriv1(deriv1, coll_solution(0), DISTYPE);

    static LINALG::TMatrix<FAD,1,numnodes> deriv2;
    DRT::UTILS::shape_function_1D_deriv2(deriv2, coll_solution(0), DISTYPE);

    // position and velocity of wall collision point at collision time and derivs (unrolled due to FAD)
    switch(numnodes)
    {
    case 2:
      wallcollpoint_pos_fad(0) = xyze_colltime(0,0) * funct(0) + xyze_colltime(0,1) * funct(1);
      wallcollpoint_pos_fad(1) = xyze_colltime(1,0) * funct(0) + xyze_colltime(1,1) * funct(1);
      wallcollpoint_pos_fad(2) = xyze_colltime(2,0) * funct(0) + xyze_colltime(2,1) * funct(1);

      wallcollpoint_vel_fad(0) = vele_fad(0,0) * funct(0) + vele_fad(0,1) * funct(1);
      wallcollpoint_vel_fad(1) = vele_fad(1,0) * funct(0) + vele_fad(1,1) * funct(1);
      wallcollpoint_vel_fad(2) = vele_fad(2,0) * funct(0) + vele_fad(2,1) * funct(1);

      wallcollpoint_deriv_vel_fad(0) = vele_fad(0,0) * deriv1(0,0) + vele_fad(0,1) * deriv1(0,1);
      wallcollpoint_deriv_vel_fad(1) = vele_fad(1,0) * deriv1(0,0) + vele_fad(1,1) * deriv1(0,1);
      wallcollpoint_deriv_vel_fad(2) = vele_fad(2,0) * deriv1(0,0) + vele_fad(2,1) * deriv1(0,1);

      F_deriv1(0) = xyze_colltime(0,0) * deriv1(0,0) + xyze_colltime(0,1) * deriv1(0,1);
      F_deriv1(1) = xyze_colltime(1,0) * deriv1(0,0) + xyze_colltime(1,1) * deriv1(0,1);
      F_deriv1(2) = xyze_colltime(2,0) * deriv1(0,0) + xyze_colltime(2,1) * deriv1(0,1);

      F_deriv2(0) = xyze_colltime(0,0) * deriv2(0,0) + xyze_colltime(0,1) * deriv2(0,1);
      F_deriv2(1) = xyze_colltime(1,0) * deriv2(0,0) + xyze_colltime(1,1) * deriv2(0,1);
      F_deriv2(2) = xyze_colltime(2,0) * deriv2(0,0) + xyze_colltime(2,1) * deriv2(0,1);
      break;
    case 3:
      wallcollpoint_pos_fad(0) = xyze_colltime(0,0) * funct(0) + xyze_colltime(0,1) * funct(1) + xyze_colltime(0,2) * funct(2);
      wallcollpoint_pos_fad(1) = xyze_colltime(1,0) * funct(0) + xyze_colltime(1,1) * funct(1) + xyze_colltime(1,2) * funct(2);
      wallcollpoint_pos_fad(2) = xyze_colltime(2,0) * funct(0) + xyze_colltime(2,1) * funct(1) + xyze_colltime(2,2) * funct(2);

      wallcollpoint_vel_fad(0) = vele_fad(0,0) * funct(0) + vele_fad(0,1) * funct(1) + vele_fad(0,2) * funct(2);
      wallcollpoint_vel_fad(1) = vele_fad(1,0) * funct(0) + vele_fad(1,1) * funct(1) + vele_fad(1,2) * funct(2);
      wallcollpoint_vel_fad(2) = vele_fad(2,0) * funct(0) + vele_fad(2,1) * funct(1) + vele_fad(2,2) * funct(2);

      wallcollpoint_deriv_vel_fad(0) = vele_fad(0,0) * deriv1(0,0) + vele_fad(0,1) * deriv1(0,1) + vele_fad(0,2) * deriv1(0,2);
      wallcollpoint_deriv_vel_fad(1) = vele_fad(1,0) * deriv1(0,0) + vele_fad(1,1) * deriv1(0,1) + vele_fad(1,2) * deriv1(0,2);
      wallcollpoint_deriv_vel_fad(2) = vele_fad(2,0) * deriv1(0,0) + vele_fad(2,1) * deriv1(0,1) + vele_fad(2,2) * deriv1(0,2);

      F_deriv1(0) = xyze_colltime(0,0) * deriv1(0,0) + xyze_colltime(0,1) * deriv1(0,1) + xyze_colltime(0,2) * deriv1(0,2);
      F_deriv1(1) = xyze_colltime(1,0) * deriv1(0,0) + xyze_colltime(1,1) * deriv1(0,1) + xyze_colltime(1,2) * deriv1(0,2);
      F_deriv1(2) = xyze_colltime(2,0) * deriv1(0,0) + xyze_colltime(2,1) * deriv1(0,1) + xyze_colltime(2,2) * deriv1(0,2);

      F_deriv2(0) = xyze_colltime(0,0) * deriv2(0,0) + xyze_colltime(0,1) * deriv2(0,1) + xyze_colltime(0,2) * deriv2(0,2);
      F_deriv2(1) = xyze_colltime(1,0) * deriv2(0,0) + xyze_colltime(1,1) * deriv2(0,1) + xyze_colltime(1,2) * deriv2(0,2);
      F_deriv2(2) = xyze_colltime(2,0) * deriv2(0,0) + xyze_colltime(2,1) * deriv2(0,1) + xyze_colltime(2,2) * deriv2(0,2);
      break;
    default:
      dserror("not yet implemented");
      break;
    }

    // subtract current particle position from the edge collision point
    F.Update(1.0, wallcollpoint_pos_fad, -1.0, particle_pos_colltime);

    // compute rhs
    static LINALG::TMatrix<FAD,2,1> residual;
    static FAD dotprod_pos_pos, dotprod_pos_deriv1, dotprod_deriv1_deriv1, dotprod_pos_deriv2,
               dotprod_pos_v, dotprod_v_deriv1, dotprod_pos_derivv, pos2subtractrad2;

    dotprod_pos_pos = F(0)*F(0) + F(1)*F(1) + F(2)*F(2);
    dotprod_pos_deriv1 = F(0)*F_deriv1(0) + F(1)*F_deriv1(1) + F(2)*F_deriv1(2);
    dotprod_deriv1_deriv1 = F_deriv1(0)*F_deriv1(0) + F_deriv1(1)*F_deriv1(1) + F_deriv1(2)*F_deriv1(2);
    dotprod_pos_deriv2 = F(0)*F_deriv2(0) + F(1)*F_deriv2(1) + F(2)*F_deriv2(2);
    dotprod_pos_v = F(0)*(wallcollpoint_vel_fad(0) - particle_vel(0)) + F(1)*(wallcollpoint_vel_fad(1) - particle_vel(1)) + F(2)*(wallcollpoint_vel_fad(2) - particle_vel(2));
    dotprod_v_deriv1 = (wallcollpoint_vel_fad(0) - particle_vel(0))*F_deriv1(0) + (wallcollpoint_vel_fad(1) - particle_vel(1))*F_deriv1(1) +(wallcollpoint_vel_fad(2) - particle_vel(2))*F_deriv1(2);
    dotprod_pos_derivv = F(0)*wallcollpoint_deriv_vel_fad(0) + F(1)*wallcollpoint_deriv_vel_fad(1) + F(2)*wallcollpoint_deriv_vel_fad(2);

    pos2subtractrad2 = dotprod_pos_pos - radius*radius;

    residual(0) = pos2subtractrad2*2.0*dotprod_pos_deriv1 +
                  dotprod_pos_deriv1*(dotprod_deriv1_deriv1 + dotprod_pos_deriv2);

    residual(1) = pos2subtractrad2*2.0*dotprod_pos_v +
                  dotprod_pos_deriv1*(dotprod_v_deriv1+dotprod_pos_derivv);

    // distance between particle and edge
    const double distance_old = distance;
    distance = std::pow(dotprod_pos_pos, 0.5).val() - radius;
    const double residualnorm1 = std::abs(residual(0).val()) + std::abs(residual(1).val());

    // break if collision found or residual is small and collision not possible
    if (distance<GEO::TOL12 || (residualnorm1<1.0e-12*radius*radius && distance>GEO::TOL2*radius))
    {
      break;
    }

    // break here if there will be definitely no contact using some heuristic
    if (iter>2 && ((std::abs(distance-distance_old)<0.1*radius && distance>radius) || std::abs(coll_solution(0).val())>1.5) )
    {
#ifdef OUTPUT
    std::cout << "break due to far away of line" << std::endl;
#endif
#ifndef DEBUG
      break;
#else
      heuristic_break = true;
#endif
    }
    if (residualnorm1<1.0e-12*radius*radius && distance>GEO::TOL2*radius)
    {
#ifdef OUTPUT
    std::cout << "break due to reached tolerance and medium far away of line" << std::endl;
#endif
#ifndef DEBUG
      break;
#else
      heuristic_break = true;
#endif
    }
    const double remainingtime = std::max(2.0*remaining_dt,GEO::TOL2*dt);
    if (iter>2 && (coll_solution(1).val()>remainingtime && std::abs(coll_solution(1).val()-ttc_old)<0.2*remainingtime))
    {
#ifdef OUTPUT
    std::cout << "break due to time to collision with line too large" << std::endl;
#endif
#ifndef DEBUG
      break;
#else
      heuristic_break = true;
      if (iter == maxiter)
      {
        std::cout << "INFO: max iterations in line contact searching reached" << std::endl;
      }
#endif
    }

    // store old time to collision
    ttc_old = coll_solution(1).val();

    // compute dF/dx
    static LINALG::Matrix<2,2> A;

    for (int j=0; j<2; ++j)
    {
      A(0, j) = residual(0).dx(j);
      A(1, j) = residual(1).dx(j);
    }

    // compute rhs
    static LINALG::Matrix<2,1> b;
    for (int i=0; i<2; ++i)
    {
      b(i) = - residual(i).val();
    }

    // solve linear problem A dx = b
    static LINALG::Matrix<2,1> dx;
    const double det = LINALG::scaledGaussElimination<2>(A, b, dx);

    if (std::abs(det) < GEO::TOL14)
    {
      static LINALG::Matrix<3,1> copy_particle_vel, copy_F_deriv1;
      for(int i=0; i<3; ++i)
      {
        copy_particle_vel(i) = particle_vel_fad(i).val();
        copy_F_deriv1(i) = F_deriv1(i).val();
      }
      static LINALG::Matrix<3,1> crossproduct;
      crossproduct.CrossProduct(copy_particle_vel, copy_F_deriv1);
      if(crossproduct.Norm1() < GEO::TOL10 and iter>1)
      {
#ifdef OUTPUT
        std::cout << "particle path is parallel to line --> left iteration" << std::endl;
#endif
        // leave here
        return validlinecollision;
      }
    }

    // update of coll_solution
    for (int i=0; i<2; ++i)
    {
      coll_solution(i) += dx(i);
    }

    // in case of parallel movement of particle to wall, element coord get extremely large
    if (std::abs(coll_solution(0)) > 1.0e3)
    {
#ifdef OUTPUT
      std::cout << "elecoord of line is extremely large --> left iteration" << std::endl;
#endif
      return validlinecollision;
    }
  }

  // check whether collision position is valid and otherwise check for corner collision points
  if (GEO::checkPositionWithinElementParameterSpace(coll_solution, DISTYPE) == true and distance < GEO::TOL7*radius)
  {
#ifdef OUTPUT
      std::cout << "elecoord is within parameter space and distance is small" << std::endl;
#endif
    // check whether collision time is reasonable
    if (coll_solution(1) >= -GEO::TOL14 && coll_solution(1) < 1.1 * remaining_dt)
    {
      FAD scalar_partvel = F.Dot(particle_vel_fad);
      FAD scalar_wallvel = F.Dot(wallcollpoint_vel_fad);

      // decide if collision is still going to happen (--> valid) or has already happened in the last time step (--> invalid)
      if ((scalar_partvel - scalar_wallvel) > GEO::TOL14)
      {
        timetocollision = coll_solution(1).val();
        validlinecollision = true;
        for(int j=0; j<3; ++j)
        {
          wallcollpoint_pos(j) = wallcollpoint_pos_fad(j).val();
          wallcollpoint_vel(j) = wallcollpoint_vel_fad(j).val();
        }

#ifdef DEBUG
        {
          // advance particle in time to collision time
          static LINALG::Matrix<3,1> newpos;
          newpos.Update(1.0, particle_pos, timetocollision, particle_vel);

          static LINALG::Matrix<3,1> collnormal;
          collnormal.Update(1.0, newpos, -1.0, wallcollpoint_pos);

          // safety check
          double normallength = collnormal.Norm2();
          if (std::abs(normallength - radius) > GEO::TOL7)
          {
            std::cout << "ran into error in line coll finding :" << std::endl;
            std::cout << "particle pos: " << newpos << std::endl;
            std::cout << "wallcollpoint_pos: " << wallcollpoint_pos << std::endl;
            std::cout << "normallength: " << normallength << std::endl;
            std::cout << "radius: " << radius << std::endl;
            std::cout << "timetocollision: " << timetocollision << std::endl;
            dserror("Particle and wall collision detected but distance does not match radius");
          }
        }
        if(heuristic_break == true)
          dserror("in release version, the loop would have been left due to heuristic criterion");
#endif
      }
    }
  }
  else
  {
    checkcorners = true;
  }

  return validlinecollision;
}


/*----------------------------------------------------------------------*
 | computes time to collision between a particle and       ghamm 04/14  |
 | an element edge for hard sphere particles                            |
 *----------------------------------------------------------------------*/
bool PARTICLE::ParticleCollisionHandlerMD::ComputeCollisionOfParticleWithLine(
    const DRT::Element::DiscretizationType distype,
    const Epetra_SerialDenseMatrix& xyze_line_n,
    const Epetra_SerialDenseMatrix& xyze_line_np,
    const LINALG::Matrix<3,1>& particle_pos,
    const LINALG::Matrix<3,1>& particle_vel,
    const double radius,
    double& timetocollision,
    LINALG::Matrix<3,1>& wallcollpoint_pos,
    LINALG::Matrix<3,1>& wallcollpoint_vel,
    const double remaining_dt,
    const double dt,
    bool& checkcorners)
{
  bool validlinecollision = false;
  switch (distype)
  {
  case DRT::Element::line2:
//    validlinecollision = ComputeCollisionOfParticleWithLineT<DRT::Element::line2>(
//        xyze_line_n, xyze_line_np, particle_pos, particle_vel, radius,
//        timetocollision, wallcollpoint_pos, wallcollpoint_vel, checkcorners);
    validlinecollision = ComputeCollisionOfParticleWithLineT_FAD<DRT::Element::line2>(
        xyze_line_n, xyze_line_np, particle_pos, particle_vel, radius,
        timetocollision, wallcollpoint_pos, wallcollpoint_vel,remaining_dt, dt, checkcorners);
    break;
  case DRT::Element::line3:
//    validlinecollision = ComputeCollisionOfParticleWithLineT<DRT::Element::line3>(
//        xyze_line_n, xyze_line_np, particle_pos, particle_vel, radius,
//        timetocollision, wallcollpoint_pos, wallcollpoint_vel, checkcorners);
    validlinecollision = ComputeCollisionOfParticleWithLineT_FAD<DRT::Element::line3>(
        xyze_line_n, xyze_line_np, particle_pos, particle_vel, radius,
        timetocollision, wallcollpoint_pos, wallcollpoint_vel,remaining_dt, dt, checkcorners);
    break;
  default:
    dserror("please add your line element type here");
    break;
  }

  return validlinecollision;
}


/*----------------------------------------------------------------------*
 | computes time to collision between a particle and       ghamm 05/14  |
 | an element corner for hard sphere particles                          |
 *----------------------------------------------------------------------*/
bool PARTICLE::ParticleCollisionHandlerMD::ComputeCollisionOfParticleWithCorner(
    const Epetra_SerialDenseMatrix& xyze_corner_n,
    const Epetra_SerialDenseMatrix& xyze_corner_np,
    const LINALG::Matrix<3,1>& particle_pos,
    const LINALG::Matrix<3,1>& particle_vel,
    const double radius,
    double& timetocollision,
    LINALG::Matrix<3,1>& wallcollpoint_pos,
    LINALG::Matrix<3,1>& wallcollpoint_vel,
    const double remaining_dt,
    const double dt)
{
  bool validcornercollision = false;

  // velocity of wall element is constant over the time step
  Epetra_SerialDenseMatrix vcorner(xyze_corner_n);
  vcorner.Scale(-1.0);
  vcorner += xyze_corner_np;
  static LINALG::Matrix<3,1> corner_vel;
  for(int i=0; i<3; ++i)
    corner_vel(i) = vcorner(i,0);
  corner_vel.Scale(1.0 / dt);

  static LINALG::Matrix<3,1> deltax, deltav;

  for (int i=0; i<3; ++i)
  {
    deltax(i) = xyze_corner_n(i,0) + (dt - remaining_dt) *  corner_vel(i) - particle_pos(i) ;
    deltav(i) = corner_vel(i) - particle_vel(i);
  }

  const double sigma = radius;

  // finding possible collision times
  const double a = deltav.Dot(deltav);
  const double b = 2.0 * deltav.Dot(deltax);
  const double c = deltax.Dot(deltax) - sigma * sigma;

  const double discriminant = b * b - 4.0 * a * c;

  if (discriminant >= 0.0 and a>0.0)
  {
    const double sqrdiscr = sqrt(discriminant);
    const double inv2a = 1.0 / (2.0 * a);
    const double tc1 = (-b - sqrdiscr) * inv2a;
    const double tc2 = (-b + sqrdiscr) * inv2a;

    // tc1 is negative
    if (tc1 < -GEO::TOL14 and tc2 > -GEO::TOL14)
    {
      timetocollision = tc2;
    }
    // tc2 is negative
    else if (tc1 > -GEO::TOL14 and tc2 < -GEO::TOL14)
    {
      timetocollision = tc1;
    }
    // both positive, smaller one is chosen (tc1 is almost identical to tc2)
    else if (tc1 > -GEO::TOL14 and tc2 > -GEO::TOL14)
    {
      timetocollision = std::min(tc1, tc2);
    }

    // check whether collision time is reasonable
    if (timetocollision >= -GEO::TOL14 && timetocollision < 1.1*remaining_dt)
    {
      const double scalar_partvel = particle_vel.Dot(deltax);
      const double scalar_wallvel = corner_vel.Dot(deltax);

      // decide if collision is still going to happen (--> valid) or has already happened in the last time step (--> invalid)
      if ((scalar_partvel - scalar_wallvel) > GEO::TOL14)
      {
        validcornercollision = true;
        for(int j=0; j<3; ++j)
        {
          wallcollpoint_pos(j) = xyze_corner_n(j,0) + (dt - remaining_dt + timetocollision) *  corner_vel(j);
          wallcollpoint_vel(j) = corner_vel(j);
        }
      }
    }
  }

  return validcornercollision;
}


/*----------------------------------------------------------------------*
 | necessary data for computing particle collisions        ghamm 09/13  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerMD::GetCollisionData(
  const DRT::Node* particle_1,
  const DRT::Node* particle_2,
  LINALG::Matrix<3,1>& pos_1,
  LINALG::Matrix<3,1>& pos_2,
  LINALG::Matrix<3,1>& vel_1,
  LINALG::Matrix<3,1>& vel_2,
  double& rad_1,
  double& rad_2,
  double& ddt_1,
  double& ddt_2
  )
{
  // data of particle 1
  if(particle_1 != NULL)
  {
    GetCollisionData(particle_1, pos_1, vel_1, rad_1, ddt_1);
  }

  // data of particle 2
  if(particle_2 != NULL)
  {
    GetCollisionData(particle_2, pos_2, vel_2, rad_2, ddt_2);
  }

  return;
}


/*----------------------------------------------------------------------*
 | necessary data for computing particle collisions        ghamm 09/13  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerMD::GetCollisionData(
  const DRT::Node* particle,
  LINALG::Matrix<3,1>& dis,
  LINALG::Matrix<3,1>& vel,
  double& rad,
  double& ddt
  )
{
  // data of particle 1
  const int lid = particle->LID();
  dis = particleData_[lid].dis;
  vel = particleData_[lid].vel;
  rad = particleData_[lid].rad;
  ddt = particleData_[lid].ddt;

  return;
}


/*----------------------------------------------------------------------*
 | operator to sort two events in increasing time          ghamm 09/13  |
 *----------------------------------------------------------------------*/
bool PARTICLE::Event::Helper::operator()(Teuchos::RCP<Event> event1, Teuchos::RCP<Event> event2)
{
  if (event1->coltype == INPAR::PARTICLE::particle_particle and event2->coltype == INPAR::PARTICLE::particle_particle)
  {
    // compare two events with particle-particle-collision
    // check if these two events are happening simultaneously
    if (std::abs(event1->time - event2->time) < GEO::TOL14)
    {
      int gid1 = event1->particle_1->Id();
      int gid2 = event1->particle_2->Id();
      int gid3 = event2->particle_1->Id();
      int gid4 = event2->particle_2->Id();

      // check if particles are involved in multiple events
      if (gid1 == gid3 || gid1 == gid4 || gid2 == gid3 || gid2 == gid4)
      {
        std::cout << ("ERROR: THREE PARTICLES COLLIDING AT THE SAME TIME!") << std::endl;
      }
      // in order to return proper sorted events use gids although binary collision is not true
      if (gid1 == gid3)
      {
        return gid2 < gid4;
      }
      else if (gid1 == gid4)
      {
        return gid2 < gid3;
      }
      else if (gid2 == gid3)
      {
        return gid1 < gid4;
      }
      else if (gid2 == gid4)
      {
        return gid1 < gid3;
      }
      else
        return gid1 < gid3;

    }
    else
    {
      // return correct order with increasing time to event
      return event1->time < event2->time;
    }
  }
  else if (event1->coltype == INPAR::PARTICLE::particle_particle and event2->coltype == INPAR::PARTICLE::particle_wall)
  {
    // compare event particle-particle-collision with event particle-wall-collision
    // check if these two events are happening simultaneously
    if (std::abs(event1->time - event2->time) < GEO::TOL14)
    {
      int gid1 = event1->particle_1->Id();
      int gid2 = event1->particle_2->Id();
      int gid3 = event2->particle_1->Id();

      if (gid1 == gid3 || gid2 == gid3)
      {
        std::cout << ("NOTE: TWO PARTICLES AND WALL COLLIDING AT THE SAME TIME!") << std::endl;
      }

      if (gid1 == gid3)
      {
        // wall collision is inserted before inter-particle collision
        return false;
      }
      else
        return gid1 < gid3;
    }
    else
    {
      // return correct order with increasing time to event
      return event1->time < event2->time;
    }
  }
  else if (event1->coltype == INPAR::PARTICLE::particle_wall and event2->coltype == INPAR::PARTICLE::particle_particle)
  {
    // compare event particle-wall-collision with event particle-particle-collision
    // check if these two events are happening simultaneously
    if (std::abs(event1->time - event2->time) < GEO::TOL14)
    {
      int gid1 = event1->particle_1->Id();
      int gid3 = event2->particle_1->Id();
      int gid4 = event2->particle_2->Id();

      if (gid1 == gid3 || gid1 == gid4)
      {
        std::cout << ("NOTE: TWO PARTICLES AND WALL COLLIDING AT THE SAME TIME!") << std::endl;
      }

      if (gid1 == gid3)
      {
        // wall collision is inserted before inter-particle collision
        return true;
      }
      else
        return gid1 < gid3;
    }
    else
    {
      // return correct order with increasing time to event
      return event1->time < event2->time;
    }
  }
  else if (event1->coltype == INPAR::PARTICLE::particle_wall and event2->coltype == INPAR::PARTICLE::particle_wall)
  {
    // compare event particle-wall-collision with event particle-wall-collision
    // check if these two events are happening simultaneously
    if (std::abs(event1->time - event2->time) < GEO::TOL14)
    {
      int gid1 = event1->particle_1->Id();
      int gid3 = event2->particle_1->Id();

      if (gid1 == gid3)
      {
        std::cout << ("INFO: PARTICLE AND TWO WALLS COLLIDING AT THE SAME TIME OR EDGE CONTACT!") << std::endl;
        return Teuchos::rcp_dynamic_cast<WallEvent>(event1,true)->wall->Id() < Teuchos::rcp_dynamic_cast<WallEvent>(event2,true)->wall->Id();
      }
      else
        return gid1 < gid3;
    }
    else
    {
      // return correct order with increasing time to event
      return event1->time < event2->time;
    }
  }

  dserror("you should not show up here");

  // order in set time increase
  return event1->time < event2->time;
}
