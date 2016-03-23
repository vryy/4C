/*----------------------------------------------------------------------*/
/*!
\file particle_contact.cpp

\brief Particle collision handling

\maintainer Georg Hammerl
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "particle_contact.H"
#include "particle_algorithm.H"
#include "../drt_adapter/ad_str_structure.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_mat/particle_mat.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_meshfree_discret/drt_meshfree_multibin.H"
#include "../drt_geometry/searchtree_geometry_service.H"
#include "../drt_geometry/position_array.H"
#include "../drt_geometry/element_coordtrafo.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include "../drt_mat/stvenantkirchhoff.H"

#include <Teuchos_TimeMonitor.hpp>
#include <Epetra_FEVector.h>
#include <Sacado.hpp>

typedef Sacado::Fad::DFad<double> FAD;

//#define OUTPUT


PARTICLE::ParticleCollisionHandlerBase::ParticleCollisionHandlerBase(
  Teuchos::RCP<DRT::Discretization> discret,
  Teuchos::RCP<PARTICLE::Algorithm> particlealgorithm,
  const Teuchos::ParameterList& particledynparams
  ) :
  myrank_(discret->Comm().MyPID()),
  discret_(discret),
  particle_algorithm_(particlealgorithm),
  contact_energy_(0.0),
  g_max_(0.0),
  writeenergyevery_(particledynparams.get<int>("RESEVRYERGY")),
  radiusncol_(Teuchos::null),
  masscol_(Teuchos::null),
  disncol_(Teuchos::null),
  velncol_(Teuchos::null),
  ang_velncol_(Teuchos::null)
{
  // make sure that a particle material is defined in the dat-file
  int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_particlemat);
  if (id==-1)
    dserror("Could not find particle material");

  const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
  const MAT::PAR::ParticleMat* actmat = static_cast<const MAT::PAR::ParticleMat*>(mat);
  // currently all particles have identical density and radius
  double density = actmat->density_;
  nue_ = actmat->poissonratio_;
  young_ = actmat->young_;

  ReadContactParameters(density);

  return;
}


/*----------------------------------------------------------------------*
 | set states from time integrator to prepare collisions   ghamm 09/13  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerBase::SetState(
  Teuchos::RCP<Epetra_Vector> radius,
  Teuchos::RCP<Epetra_Vector> mass)
{
  // node based vectors
  radiusncol_ = LINALG::CreateVector(*discret_->NodeColMap(),false);
  LINALG::Export(*radius,*radiusncol_);
  masscol_ = LINALG::CreateVector(*discret_->NodeColMap(),false);
  LINALG::Export(*mass,*masscol_);

  // miraculous transformation from row to col layout ...
  disncol_ = Teuchos::rcp(new Epetra_Vector(*discret_->GetState("bubblepos")));
  velncol_ = Teuchos::rcp(new Epetra_Vector(*discret_->GetState("bubblevel")));
  ang_velncol_ = discret_->GetState("bubbleangvel");

  return;
}


/*----------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerBase::GetNeighbouringParticlesAndWalls(
    DRT::Node* particle,
    std::set<DRT::Node*>& neighboringparticles,
    std::set<DRT::Element*>& neighboringwalls)
{
  if (particle->NumElement() != 1)
    dserror("More than one element for this particle");

  DRT::Element** CurrentBin = particle->Elements();
  const int binId = CurrentBin[0]->Id();

  int ijk[3];
  particle_algorithm_->ConvertGidToijk(binId,ijk);

  // ijk_range contains: i_min   i_max     j_min     j_max    k_min     k_max
  int ijk_range[] = {ijk[0]-1, ijk[0]+1, ijk[1]-1, ijk[1]+1, ijk[2]-1, ijk[2]+1};
  std::vector<int> binIds;
  binIds.reserve(27);

  // do not check on existence here -> shifted to GetBinContent
  particle_algorithm_->GidsInijkRange(ijk_range,binIds,false);

  GetBinContent(neighboringparticles, neighboringwalls, binIds);

  // delete particle (there is no contact between particle i and particle i)
  neighboringparticles.erase(particle);

  return;
}


/*----------------------------------------------------------------------*
 | get particles and wall elements in given bins           ghamm 09/13  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerBase::GetBinContent(
  std::set<DRT::Node*> &particles,
  std::set<DRT::Element*> &walls,
  std::vector<int> &binIds
  )
{
  // loop over all bins
  for(std::vector<int>::const_iterator bin=binIds.begin(); bin!=binIds.end(); ++bin)
  {
    // extract bins from discretization after checking on existence
    const int lid = discret_->ElementColMap()->LID(*bin);
    if(lid<0)
      continue;
    DRT::MESHFREE::MeshfreeMultiBin* neighboringbin =
        dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(discret_->lColElement(lid));

    // gather wall elements
    DRT::Element** walleles = neighboringbin->AssociatedWallEles();
    int numwalls = neighboringbin->NumAssociatedWallEle();
    for(int iwall=0;iwall<numwalls; ++iwall)
      walls.insert(walleles[iwall]);

    // gather particles
    DRT::Node **nodes = neighboringbin->Nodes();
    int numparticles = neighboringbin->NumNode();
    for(int inode=0; inode<numparticles; ++inode)
      particles.insert(nodes[inode]);
  }

  return;
}


/*----------------------------------------------------------------------*
 | read initial contact parameters and validate them       ghamm 09/13  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerBase::ReadContactParameters(double density)
{
  // extract input parameters
  const Teuchos::ParameterList& particleparams = DRT::Problem::Instance()->ParticleParams();
  contact_strategy_ = DRT::INPUT::IntegralValue<INPAR::PARTICLE::ContactStrategy>(particleparams,"CONTACT_STRATEGY");
  normal_contact_ = DRT::INPUT::IntegralValue<INPAR::PARTICLE::NormalContact>(particleparams,"NORMAL_CONTACT_LAW");

  if(contact_strategy_ != INPAR::PARTICLE::None)
  {
    r_min_ = particleparams.get<double>("MIN_RADIUS");
    r_max_ = particleparams.get<double>("MAX_RADIUS");
    v_max_ = particleparams.get<double>("MAX_VELOCITY");
    c_ = particleparams.get<double>("REL_PENETRATION");
    e_ = particleparams.get<double>("COEFF_RESTITUTION");
    e_wall_ = particleparams.get<double>("COEFF_RESTITUTION_WALL");
    mu_wall_ = particleparams.get<double>("FRICT_COEFF_WALL");
    mu_ = particleparams.get<double>("FRICT_COEFF");
    tension_cutoff_ = (DRT::INPUT::IntegralValue<int>(particleparams,"TENSION_CUTOFF") == 1);

    if(r_min_<0.0 or r_max_<0.0 or v_max_<0.0 or c_<0.0 or young_<0.0)
      dserror("Invalid input parameter (MIN_RADIUS,MAX_RADIUS,MAX_VELOCITY,REL_PENETRATION, YOUNG's modulus have to be larger than zero)");

    if((e_<0.0 or e_wall_<0.0) and (normal_contact_ == INPAR::PARTICLE::LinSpringDamp or contact_strategy_==INPAR::PARTICLE::Normal_MD))
      dserror("Invalid input parameter COEFF_RESTITUTION for this kind of contact law");

    // no further information necessary for MD like contact
    if(contact_strategy_ == INPAR::PARTICLE::Normal_MD)
      return;

    // data for critical time step computation for DEM like contact
    const double mass_min = density * 4.0/3.0 * M_PI * pow( r_min_ ,3.0 );

    double k_tkrit = 0.0;

    // wall material properties are always taken from the first available St. Venant Kirchhoff material
    int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_stvenant);
    if (id==-1)
      dserror("Could not find wall material which is assumed to be the first St. Venant Kirchhoff material");
    Teuchos::RCP<MAT::StVenantKirchhoff> wallmat = Teuchos::rcp_dynamic_cast<MAT::StVenantKirchhoff>(MAT::Material::Factory(id));

    const double G_wall = wallmat->ShearMod();
    const double nue_wall = wallmat->PoissonRatio();
    if(G_wall<0.0)
      dserror("Wall shear modulus hasto be larger than zero");

    const double G = young_ / (2.0*(1.0+nue_));

    // kappa - tangential to normal stiffness ratio
    kappa_ = (1.0-nue_)/(1.0-0.5*nue_);
    kappa_wall_ = ( (1.0-nue_)/G + (1.0-nue_wall)/G_wall ) / ( (1.0-0.5*nue_)/G + (1.0-0.5*nue_wall)/G_wall );

    //------------stiffness----------------------------
    switch(normal_contact_)
    {
    case INPAR::PARTICLE::LinSpring:
    case INPAR::PARTICLE::LinSpringDamp:
    {
      // stiffness calculated from relative penetration and some other input parameters (linear spring)
      k_normal_ = 2.0/3.0 * r_max_ * M_PI * density * v_max_ * v_max_ / (c_ * c_);
      // for tangential contact same stiffness is used
      k_tang_ = kappa_ * k_normal_;
      k_tang_wall_ = kappa_wall_ * k_normal_;
      k_tkrit = k_normal_;

      double user_normal_stiffness = particleparams.get<double>("NORMAL_STIFF");
      if(user_normal_stiffness > 0.0)
      {
        // if user specifies normal stiffness, this stiffness will be used as normal and tangential stiffness for simulation
        k_normal_ = user_normal_stiffness;
        // for tangential contact same stiffness is used
        k_tang_ =kappa_ * k_normal_;
        k_tang_wall_ =kappa_wall_ * k_normal_;
        // stiffness used for calculation of critical time step
        k_tkrit = k_normal_;

        std::cout<<"WARNING: stiffness calculated from relative penetration will be overwritten by input NORMAL_STIFF!!!"<<std::endl;
      }
    }
    break;
    case INPAR::PARTICLE::Hertz:
    case INPAR::PARTICLE::LeeHerrmann:
    case INPAR::PARTICLE::KuwabaraKono:
    case INPAR::PARTICLE::Tsuji:
    {
      if(contact_strategy_==INPAR::PARTICLE::NormalAndTang_DEM)
        dserror("tangential contact only with linear normal model implemented");

      // stiffness calculated from relative penetration and some other input parameters (Hertz)
      k_normal_ = 10.0/3.0 * M_PI * density * v_max_ * v_max_ * pow(r_max_,0.5) / pow(2.0*c_,2.5);
      // stiffness used for calculation of critical time step (linear spring stiffness needed!)
      k_tkrit = 2.0/3.0 * r_max_ * M_PI * density * v_max_ * v_max_ / (c_ * c_);

      double user_normal_stiffness = particleparams.get<double>("NORMAL_STIFF");
      if(user_normal_stiffness > 0.0)
      {
        // if user specifies normal stiffness, this stiffness will be used as normal stiffness for simulation
        k_normal_ = user_normal_stiffness;
        // for tangential contact the user specified (nonlinear) normal stiffness which has to be transformed into a linear normal
        // stiffness with the same relative penetration which is used as (linear) tangential stiffness afterwards
        const double value = 2048.0/1875.0 * density * v_max_ * v_max_ * M_PI * pow(r_max_,3.0) * pow(k_normal_,4.0);
        // stiffness used for calculation of critical time step (linear spring stiffness needed!)
        k_tkrit = pow(value,0.2);

        std::cout<<"WARNING: stiffness calculated from relative penetration will be overwritten by input NORMAL_STIFF!!!"<<std::endl;
      }
    }
    break;
    default:
      dserror("normal contact law does not exist");
    break;
    }
    //---------------------------------------------------------

    //------------------damping--------------------------------
    d_normal_ = -1.0;
    d_tang_ = -1.0;

    if(normal_contact_ == INPAR::PARTICLE::LinSpringDamp)
    {
      double user_normal_damping = particleparams.get<double>("NORMAL_DAMP");
      if(user_normal_damping >= 0.0)
      {
        dserror("Invalid input parameter NORMAL_DAMP for this kind of contact law");
      }
      double user_tang_damping = particleparams.get<double>("TANG_DAMP");
      if(user_tang_damping >= 0.0)
      {
        dserror("Invalid input parameter TANG_DAMP for this kind of contact law");
      }
    }

    if(normal_contact_ == INPAR::PARTICLE::LeeHerrmann || normal_contact_ == INPAR::PARTICLE::KuwabaraKono ||
        normal_contact_ == INPAR::PARTICLE::Tsuji)
    {
      double user_normal_damping = particleparams.get<double>("NORMAL_DAMP");
      if(user_normal_damping >= 0.0)
      {
        // user has to specify normal damping coefficient
        d_normal_ = user_normal_damping;
      }
      else
      {
        dserror("For this kind of contact law the input parameter NORMAL_DAMP is invalid");
      }
      double user_tang_damping = particleparams.get<double>("TANG_DAMP");
      if(user_tang_damping >= 0.0)
      {
        // user has to specify tangential damping coefficient
        d_tang_ = user_tang_damping;
      }
      else
      {
        if(contact_strategy_==INPAR::PARTICLE::NormalAndTang_DEM)
          dserror("For this kind of contact law the input parameter TANG_DAMP is invalid");
      }
    }
    //------------------damping (wall)--------------------------------
    d_normal_wall_ = -1.0;
    d_tang_wall_ = -1.0;

    if(normal_contact_ == INPAR::PARTICLE::LinSpringDamp)
    {
      double user_normal_damping = particleparams.get<double>("NORMAL_DAMP_WALL");
      if(user_normal_damping >= 0.0)
      {
        dserror("Invalid input parameter NORMAL_DAMP_WALL for this kind of contact law");
      }
      double user_tang_damping = particleparams.get<double>("TANG_DAMP_WALL");
      if(user_tang_damping >= 0.0)
      {
        dserror("Invalid input parameter TANG_DAMP_WALL for this kind of contact law");
      }
    }

    if(normal_contact_ == INPAR::PARTICLE::LeeHerrmann || normal_contact_ == INPAR::PARTICLE::KuwabaraKono ||
        normal_contact_ == INPAR::PARTICLE::Tsuji)
    {
      double user_normal_damping = particleparams.get<double>("NORMAL_DAMP_WALL");
      if(user_normal_damping >= 0.0)
      {
        // user has to specify normal damping coefficient
        d_normal_wall_ = user_normal_damping;
      }
      else
      {
        dserror("For this kind of contact law the input parameter NORMAL_DAMP_WALL is invalid");
      }
      double user_tang_damping = particleparams.get<double>("TANG_DAMP_WALL");
      if(user_tang_damping >= 0.0)
      {
        // user has to specify tangential damping coefficient
        d_tang_wall_ = user_tang_damping;
      }
      else
      {
        if(contact_strategy_==INPAR::PARTICLE::NormalAndTang_DEM)
          dserror("For this kind of contact law the input parameter TANG_DAMP_WALL is invalid");
      }

    }
    //---------------------------------------------------------------

    double factor = 0.0;
    double safety = 0.75;

    // initialize factor
    if(contact_strategy_==INPAR::PARTICLE::Normal_DEM)
    { factor = 0.34; }
    if(contact_strategy_==INPAR::PARTICLE::NormalAndTang_DEM)
    { factor = 0.22; }

    // calculate critical time step for DEM like contact
    dt_krit_ = safety * factor * sqrt( mass_min / k_tkrit );

    double dt = particleparams.get<double>("TIMESTEP");
    if(dt>dt_krit_ and myrank_==0)
    {
      std::cout << "\n\nW A R N I N G : time step larger than critical time step!" << std::endl;
      std::cout << "W A R N I N G : calculated critical time step: "<<dt_krit_<<" !\n\n" << std::endl;
    }

    // check frictional coefficient
    if(contact_strategy_ == INPAR::PARTICLE::NormalAndTang_DEM)
    {
      if(mu_<=0.0 or mu_wall_<=0.0)
       dserror("Friction coefficient invalid");
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | assemble energies of particles                          ghamm 09/13  |
 *----------------------------------------------------------------------*/
double PARTICLE::ParticleCollisionHandlerBase::EnergyAssemble(double owner_i,double owner_j)
{
  double value = -1.0;

  const int myrank = discret_->Comm().MyPID();

  //contact  with wall
  if(owner_j<0)
  {
    if(owner_i!=myrank)
     value = 0.0;
  }
  //contact between particle_i and particle_j
  else
  {
    if((owner_i==myrank && owner_j!=myrank) or (owner_i!=myrank && owner_j==myrank))
      value = 0.5;
    else if(owner_i!=myrank and owner_j!=myrank)
      value = 0.0;
  }

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
    )
{
  return;
}


/*----------------------------------------------------------------------*
 | compute collisions (inter-particle and particle-wall)   ghamm 09/13  |
 *----------------------------------------------------------------------*/
double PARTICLE::ParticleCollisionHandlerDEM::EvaluateParticleContact(
  const double dt,
  Teuchos::RCP<Epetra_Vector> f_contact,
  Teuchos::RCP<Epetra_Vector> m_contact
  )
{
  TEUCHOS_FUNC_TIME_MONITOR("TimeforContactSearchAndCalculation");

  contact_energy_ = 0.0;

  // get wall discretization and states for particles
  Teuchos::RCP<DRT::Discretization> walldiscret = particle_algorithm_->WallDiscret();
  Teuchos::RCP<const Epetra_Vector> walldisn = walldiscret->GetState("walldisnp");
  Teuchos::RCP<const Epetra_Vector> wallveln = walldiscret->GetState("wallvelnp");

  // define vector for contact force
  Teuchos::RCP<Epetra_FEVector> f_structure = Teuchos::rcp(new Epetra_FEVector(*discret_->DofRowMap()));


  // gather data for all particles initially
  const int numcolparticles = discret_->NodeColMap()->NumMyElements();
  particledata_.resize(numcolparticles);
  {
    for (int i=0; i<numcolparticles; ++i)
    {
      // particle for which data will be collected
      DRT::Node *particle = discret_->lColNode(i);

      static LINALG::Matrix<3,1> position;
      static LINALG::Matrix<3,1> velocity;
      static LINALG::Matrix<3,1> angvel;

      std::vector<int> lm;
      lm.reserve(3);

      // extract global dof ids and fill into lm_i
      discret_->Dof(particle, lm);

      //position, velocity and angular velocity of particle
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*disncol_,position,lm);
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*velncol_,velocity,lm);
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*ang_velncol_,angvel,lm);

      const int lid = particle->LID();
      // radius and mass of particle
      const double radius = (*radiusncol_)[lid];
      const double mass = (*masscol_)[lid];

      ParticleCollData& data = particledata_[particle->LID()];
      data.pos = position;
      data.vel = velocity;
      data.angvel = angvel;
      data.rad = radius;
      data.mass = mass;
      data.lm = lm;
    }
  }


  // store bins, which have already been examined
  std::set<int> examinedbins;

  // loop over all particles
  for(int i=0; i<numcolparticles; ++i)
  {
    DRT::Node *currparticle = discret_->lColNode(i);

    if(currparticle->NumElement() != 1)
      dserror("More than one element for this particle");

    DRT::Element** CurrentBin = currparticle->Elements();
    int binId = CurrentBin[0]->Id();

    // if a bin has already been examined --> continue with next particle
    if( examinedbins.count(binId) == 1 )
    {
      continue;
    }
    //else: bin is examined for the first time --> new entry in examinedbins_
    else
    {
      examinedbins.insert(binId);
    }

    // list of all particles in the neighborhood of currparticle
    std::set<DRT::Node*> neighboringparticles;

    // list of walls that border on the CurrentBin
    std::set<DRT::Element*> neighboringwalls;

    GetNeighbouringParticlesAndWalls(currparticle, neighboringparticles, neighboringwalls);

    DRT::Node **NodesInCurrentBin = CurrentBin[0]->Nodes();
    const int numparticle = CurrentBin[0]->NumNode();

    // loop over all particles in CurrentBin
    for(int i=0; i<numparticle; ++i)
    {
      DRT::Node *particle_i = NodesInCurrentBin[i];

      // extract data
      static LINALG::Matrix<3,1> position_i, vel_i, angvel_i;
      ParticleCollData& data = particledata_[particle_i->LID()];
      position_i = data.pos;
      vel_i = data.vel;
      angvel_i = data.angvel;
      const double radius_i = data.rad;
      const double mass_i = data.mass;
      std::vector<int> lm_i = data.lm;
      const int owner_i = particle_i->Owner();

      // evaluate contact with walls first
      std::vector<WallContactPoint> surfaces;
      std::vector<WallContactPoint> lines;
      std::vector<WallContactPoint> nodes;

      std::set<int> unusedIds;

      // check whether there is contact between particle i and neighboring walls
      for(std::set<DRT::Element*>::const_iterator w=neighboringwalls.begin(); w!=neighboringwalls.end();  ++w)
      {
        DRT::Element* neighboringwallele = (*w);
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

        LINALG::Matrix<3,1> nearestPoint;

        //-------find point on wall element with smallest distance to particle_i-------------------
        GEO::ObjectType objecttype = GEO::nearest3DObjectOnElement(neighboringwallele,nodeCoord,position_i,nearestPoint);
        //-----------------------------------------------------------------------------------------

        static LINALG::Matrix<3,1> r_i_wall;
        r_i_wall.Update(1.0, nearestPoint, -1.0, position_i);
        const double distance_i_wall = r_i_wall.Norm2();
        const double penetration = distance_i_wall-radius_i;

        if(penetration <= 0.0)
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
            WallContactPoint currentContact = { neighboringwallele->Id(), nearestPoint, penetration, nodeCoord, lm_wall, lmowner };
            (*pointer).push_back(currentContact);
          }
        }
        // penetration > 0.0 --> contact impossible
        else
        {
          unusedIds.insert(neighboringwallele->Id());
        }
      }

      // find entries of lines and nodes which are within the penetration volume of the current particle
      // hierarchical: surfaces first
      for(size_t s=0; s<surfaces.size(); ++s)
      {
        // within this radius no other contact point can lie: radius = sqrt(r_i^2 - (r_i-|g|)^2)
        const double rminusg = radius_i-std::fabs(surfaces[s].penetration);
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
      // find entries of nodes which are within the penetration volume of the current particle
      // hierarchical: lines next
      for(size_t l=0; l<lines.size(); ++l)
      {
        // radius = sqrt(r_i^2 - (r_i-|g|)^2)
        const double rminusg = radius_i-std::fabs(lines[l].penetration);
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
        dserror("Contact with more than 3 wall elements. Check whether history is deleted correctly.");

      for(size_t s=0; s<surfaces.size(); ++s)
      {
        // gid of wall element
        WallContactPoint wallcontact = surfaces[s];
        const int gid_wall = wallcontact.eleid;

        // distance-vector
        static LINALG::Matrix<3,1> r_contact;
        r_contact.Update(1.0, wallcontact.point, -1.0, position_i);

        // distance between centre of mass of two particles
        const double norm_r_contact(r_contact.Norm2());

        // normal vector
        static LINALG::Matrix<3,1> normal;
        normal.Update(1.0/norm_r_contact, r_contact);

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

        // velocity v_rel = v_i - v_wall
        static LINALG::Matrix<3,1> v_rel;
        v_rel.Update(1.,vel_i,-1.,vel_nearestPoint);

        // part of v_rel in normal-direction: v_rel * n
        const double v_rel_normal(v_rel.Dot(normal));

        // penetration
        // g = norm_r_contact - radius_i;
        const double g(wallcontact.penetration);
        if(std::fabs(g)>g_max_)
          g_max_ = std::fabs(g);
        //-------------------------------------------------------

        // normalized mass
        const double m_eff = mass_i;

        // contact force
        double normalcontactforce = 0.0;
        static LINALG::Matrix<3,1> tangentcontactforce;
        tangentcontactforce.PutScalar(0.0);

        // normal contact force between particle and wall (note: owner_j = -1)
        CalculateNormalContactForce(g, v_rel_normal, mass_i, normalcontactforce, owner_i, -1);

        if(contact_strategy_ == INPAR::PARTICLE::NormalAndTang_DEM)
        {
          // velocity v_rel = v_i + omega_i x (r'_i n) - v_wall = (v_i - v_wall) + omega_i x (r'_i n)
          static LINALG::Matrix<3,1> v_rel_rot;
          v_rel_rot.CrossProduct(angvel_i,normal);
          v_rel.Update(radius_i+g,v_rel_rot,1.);

          // velocity v_rel_tangential
          static LINALG::Matrix<3,1> v_rel_tangential;
          v_rel_tangential.Update(1., v_rel, -v_rel_normal, normal);

          // if g < 0 and g_lasttimestep > 0 -> create history variables
          if(!history_wall.count(gid_wall))
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
          CalculateTangentialContactForce(normalcontactforce, normal, tangentcontactforce,
                      history_wall[gid_wall], v_rel_tangential, m_eff, dt, owner_i, -1);
        }

        // calculation of overall contact force
        static LINALG::Matrix<3,1> contactforce;
        contactforce.Update(normalcontactforce,normal,1.,tangentcontactforce);

        // calculation of overall contact moment: m_i = (r_i * n) x F_t
        static LINALG::Matrix<3,1> contactmoment;
        contactmoment.CrossProduct(normal,tangentcontactforce);
        contactmoment.Scale(radius_i+g);

        // assembly of contact forces and moments
        static Epetra_SerialDenseVector val_i(3), m_i(3);
        static std::vector<int> lmowner_i(3);

        for(int n=0; n<3; ++n)
        {
          val_i[n] = contactforce(n);
          m_i[n] = contactmoment(n);
          lmowner_i[n] = owner_i;
        }

        // do assembly of contact moments
        LINALG::Assemble(*m_contact,m_i,lm_i,lmowner_i);

        // do assembly of contact forces
        LINALG::Assemble(*f_contact,val_i,lm_i,lmowner_i);

        // forces on wall elements
        double nodal_forces[numnodes * 3];
        for(int node=0; node<numnodes; ++node)
        {
          for(int dim=0; dim<3; ++dim)
          {
            nodal_forces[node * 3 + dim] = funct[node] *(- val_i[dim]);
          }
        }

        // assembly of contact forces on walls
        if(owner_i == myrank_)
        {
          const int err = f_structure->SumIntoGlobalValues(numnodes * 3, &(wallcontact.lm)[0], &nodal_forces[0]);
          if (err<0)
            dserror("summing into Epetra_FEVector failed");
        }
      } // end for contact points on surfaces

      if(contact_strategy_ == INPAR::PARTICLE::NormalAndTang_DEM)
      {
        //delete those entries in history_wall_ which are no longer in contact with particle_i in current time step
        for(std::set<DRT::Element*>::const_iterator w=neighboringwalls.begin(); w != neighboringwalls.end(); ++w)
        {
          const int gid_wall = (*w)->Id();
          if( unusedIds.count(gid_wall) and history_wall.count(gid_wall) )
            history_wall.erase(gid_wall);
        }
      }

      // check whether there is contact between particle i and all other particles in the neighborhood except those which
      // have a lower or equal ID than particle i (--> ignoring self-contact)
      std::map<int, PARTICLE::Collision>& history_particle = static_cast<PARTICLE::ParticleNode*>(particle_i)->Get_history_particle();
      if(history_particle.size() > 12)
        dserror("Contact with more than 12 particles particles. Check whether history is deleted correctly.");

      for(std::set<DRT::Node*>::const_iterator j=neighboringparticles.begin(); j!=neighboringparticles.end(); ++j)
      {
        DRT::Node* neighborparticle = (*j);
        const int gid_j = neighborparticle->Id();
        // evaluate contact only once!
        if(particle_i->Id() >= gid_j)
          continue;

        // extract data
        static LINALG::Matrix<3,1> position_j, vel_j, angvel_j;
        ParticleCollData& data = particledata_[neighborparticle->LID()];
        position_j = data.pos;
        vel_j = data.vel;
        angvel_j = data.angvel;
        const double radius_j = data.rad;
        const double mass_j = data.mass;
        std::vector<int> lm_j = data.lm;
        const int owner_j = neighborparticle->Owner();


        // normalized mass
        const double m_eff = mass_i * mass_j / (mass_i + mass_j);

        // distance vector and distance between two particles
        static LINALG::Matrix<3,1> r_contact;
        r_contact.Update(1.,position_j,-1,position_i);
        const double norm_r_contact(r_contact.Norm2());

        // penetration
        const double g = norm_r_contact - radius_i - radius_j;
        // in case of penetration contact forces and moments are calculated
        if(g <= 0.0)
        {
          // contact forces
          double normalcontactforce = 0.0;
          static LINALG::Matrix<3,1> tangentcontactforce(true);

          if(std::fabs(g)>g_max_)
            g_max_ = std::fabs(g);

          // velocity v_rel = v_i - v_j
          static LINALG::Matrix<3,1> v_rel;
          v_rel.Update(1.,vel_i,-1.,vel_j);

          // normal vector
          static LINALG::Matrix<3,1> normal;
          normal.Update(1./norm_r_contact,r_contact);

          // part of v_rel in normal- irection: v_rel * n
          const double v_rel_normal(v_rel.Dot(normal));

          // calculation of normal contact force
          CalculateNormalContactForce(g, v_rel_normal, m_eff, normalcontactforce, owner_i, owner_j);

          // calculation of tangential contact force
          if(contact_strategy_ == INPAR::PARTICLE::NormalAndTang_DEM)
          {
            // velocity v_rel = v_i - v_j + omega_i x (r'_i n) + omega_j x (r'_j n)
            static LINALG::Matrix<3,1> v_rel_rot;
            v_rel_rot.CrossProduct(angvel_i, normal);
            v_rel.Update(radius_i+g*0.5, v_rel_rot,1.);
            v_rel_rot.CrossProduct(angvel_j, normal);
            v_rel.Update(radius_j+g*0.5, v_rel_rot,1.);

            // velocity v_rel_tangential
            static LINALG::Matrix<3,1> v_rel_tangential;
            v_rel_tangential.Update(1.,v_rel,-v_rel_normal,normal);

            // if history variables does not exist -> create it
            if(!history_particle.count(gid_j))
            {
              PARTICLE::Collision col;
              // initialize with stick
              col.stick = true;
              //initialize g_t[3]
              for(int dim=0; dim<3; ++dim)
              {
                col.g_t[dim] = 0.0;
              }

              //insert new entry
              history_particle.insert(std::pair<int,PARTICLE::Collision>(gid_j,col));
            }

            CalculateTangentialContactForce(normalcontactforce, normal, tangentcontactforce,
                      history_particle[gid_j], v_rel_tangential, m_eff, dt, owner_i, owner_j);
          }

          // calculation of overall contact force and moment
          static LINALG::Matrix<3,1> contactforce_i, contactforce_j, contactmoment_i, contactmoment_j;
          const double r_i = radius_i + g*0.5;
          const double r_j = radius_j + g*0.5;
          contactforce_i.Update(normalcontactforce,normal,1.,tangentcontactforce);
          contactforce_j.Update(-1.,contactforce_i); // actio = reactio
          contactmoment_i.CrossProduct(normal,tangentcontactforce);
          contactmoment_i.Scale(r_i); // m_i = (r_i * n) x F_t
          contactmoment_j.Update(r_j/r_i,contactmoment_i); // m_j = r_j/r_i * m_i

          static Epetra_SerialDenseVector val_i(3), val_j(3), m_i(3), m_j(3);
          static std::vector<int> lmowner_i(3), lmowner_j(3);

          for(unsigned dim=0; dim<3; ++dim)
          {
            val_i[dim] = contactforce_i(dim);
            val_j[dim] = contactforce_j(dim);
            m_i[dim] = contactmoment_i(dim);
            m_j[dim] = contactmoment_j(dim);
            lmowner_i[dim] = owner_i;
            lmowner_j[dim] = owner_j;
          }

          // assembly of contact moments
          LINALG::Assemble(*m_contact,m_i,lm_i,lmowner_i);
          LINALG::Assemble(*m_contact,m_j,lm_j,lmowner_j);

          // assembly of contact forces
          LINALG::Assemble(*f_contact,val_i,lm_i,lmowner_i);
          LINALG::Assemble(*f_contact,val_j,lm_j,lmowner_j);

        }
        else // g > 0.0 --> no contact
        {
          // erase entry in history if still existing
          if(history_particle.count(gid_j))
          {
           history_particle.erase(gid_j);
          }
        }

      }
    }
  }

  // erase temporary storage for collision data
  particledata_.clear();

  // call global assemble for particle forces on walls
  const int err = f_structure->GlobalAssemble(Add, false);
  if (err<0)
    dserror("global assemble into fluidforces failed");


  radiusncol_ = Teuchos::null;
  masscol_ =Teuchos::null;
  velncol_ = Teuchos::null;
  disncol_ = Teuchos::null;
  ang_velncol_ = Teuchos::null;

  //apply contact-forces
// structure->SetForceInterface(f_structure);

//  Teuchos::RCP<Epetra_Vector> forces_col_vector = LINALG::CreateVector(*structure->Discretization()->DofColMap());
//  LINALG::Export(*f_structure, *forces_col_vector);
//
//  double force=0.0;
//  double total_force=0.0;
//  for(int n=0; n<structure->Discretization()->NumMyRowElements(); n++)
//  {
//    DRT::Element* Ele = structure->Discretization()->lRowElement(n);
//    //only second element
//    if(Ele->Id()==1)
//    {
//      DRT::Node** Nodes=Ele->Nodes();
//      int numnodes = Ele->NumNode();
//      for(int i=0;i<numnodes;++i)
//      {
//        std::vector<int> lm;
//        lm.reserve(3);
//        structure->Discretization()->Dof(Nodes[i], lm);
//        std::vector<double> nodal_force_vector(3);
//        DRT::UTILS::ExtractMyValues(*forces_col_vector,nodal_force_vector,lm);
//        //force in y-direction
//        force += nodal_force_vector[1];
//      }
//    }
//  }

//  structure->Discretization()->Comm().SumAll(&force,&total_force,1);
//
//  if(timen_>=0.25 and timen_ <= 1.0)
//  {
//    //delta_work= force * velocity * delta_t
//    work_ += total_force*0.16*dt;
//  }

  return contact_energy_;
}


/*----------------------------------------------------------------------*
 | calculate normal contact force for single contact pair  ghamm 09/13  |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerDEM::CalculateNormalContactForce(
  const double g,
  const double v_rel_normal,
  const double m_eff,
  double& normalcontactforce,
  const int owner_i,
  const int owner_j
  )
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
      d = 2.0 * std::fabs(log(e_wall_)) * sqrt(k_normal_ * m_eff / (pow(log(e_wall_),2.0) + M_PI*M_PI));
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
      d = 2.0 * std::fabs(log(e_)) * sqrt(k_normal_ * m_eff / (pow(log(e_),2.0) + M_PI*M_PI));
    }
    else
    {
      d = d_normal_;
    }
  }

  // contact force
  switch(normal_contact_)
  {
  case INPAR::PARTICLE::LinSpring:
  {
    normalcontactforce = k_normal_ * g;

    if(writeenergyevery_)
    {
      //monitor E N E R G Y: here: calculate energy of elastic contact
      contact_energy_ += EnergyAssemble(owner_i,owner_j)* 1.0/2.0 * k_normal_ * g * g;
    }
  }
  break;
  case INPAR::PARTICLE::Hertz:
  {
    normalcontactforce = - k_normal_ * pow(-g,1.5);

    if(writeenergyevery_)
    {
      //monitor E N E R G Y: here: calculate energy of elastic contact
      contact_energy_ += EnergyAssemble(owner_i,owner_j)* 0.4 * k_normal_ * pow(-g,2.5);
    }
  }
  break;
  case INPAR::PARTICLE::LinSpringDamp:
  {
    normalcontactforce = k_normal_ * g - d * v_rel_normal;

    // tension-cutoff
    if(tension_cutoff_)
    {
      if(normalcontactforce>0.0)
      {
        normalcontactforce = 0.0;
      }
    }

    if(writeenergyevery_)
    {
      //monitor E N E R G Y: here: calculate energy of elastic contact
      contact_energy_ += EnergyAssemble(owner_i,owner_j)* 1.0/2.0 * k_normal_ * g * g;
    }
  }
  break;
  case INPAR::PARTICLE::LeeHerrmann:
  {
    // m_eff = m_i * m_j / ( m_i + m_j)

    normalcontactforce = - k_normal_ * pow(-g,1.5) - m_eff * d * v_rel_normal;

    // tension-cutoff
    if(tension_cutoff_)
    {
      if(normalcontactforce>0.0)
      {
        normalcontactforce = 0.0;
      }
    }

    if(writeenergyevery_)
    {
      //monitor E N E R G Y: here: calculate energy of elastic contact
      contact_energy_ += EnergyAssemble(owner_i,owner_j)* 0.4 * k_normal_ * pow(-g,2.5);
    }
  }
  break;
  case INPAR::PARTICLE::KuwabaraKono:
  {
    normalcontactforce = - k_normal_ * pow(-g,1.5) - d * v_rel_normal * pow(-g,0.5);

    // tension-cutoff
    if(tension_cutoff_)
    {
      if(normalcontactforce>0.0)
      {
        normalcontactforce = 0.0;
      }
    }

    if(writeenergyevery_)
    {
      //monitor E N E R G Y: here: calculate energy of elastic contact
      contact_energy_ += EnergyAssemble(owner_i,owner_j)* 0.4 * k_normal_ * pow(-g,2.5);
    }
  }
  break;
  case INPAR::PARTICLE::Tsuji:
  {
    normalcontactforce = - k_normal_ * pow(-g,1.5) - d * v_rel_normal * pow(-g,0.25);

    // tension-cutoff
    if(tension_cutoff_)
    {
      if(normalcontactforce>0.0)
      {
        normalcontactforce = 0.0;
      }
    }

    if(writeenergyevery_)
    {
      //monitor E N E R G Y: here: calculate energy of elastic contact
      contact_energy_ += EnergyAssemble(owner_i,owner_j)* 0.4 * k_normal_ * pow(-g,2.5);
    }
  }
  break;
  default:
    dserror("specified normal contact law does not exist");
    break;
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
      d = 2.0 * std::fabs(log(e_wall_)) * sqrt(k_normal_ * m_eff / (pow(log(e_wall_),2.0) + M_PI*M_PI));
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
      d = 2.0 * std::fabs(log(e_)) * sqrt(k_normal_ * m_eff / (pow(log(e_),2.0) + M_PI*M_PI));
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
  if( norm_f_t <= (mu * std::fabs(normalcontactforce)) )
  {
    currentColl.stick = true;
    //tangential contact force already calculated
  }
  else //"slip"-case
  {
    currentColl.stick = false;
    //calculate tangent vector ( unit vector in (test-)tangentcontactforce-direction )
    static LINALG::Matrix<3,1> tangent;
    tangent.Update(1./norm_f_t,tangentcontactforce);

    // calculate tangent contact force and tangential displacements
    tangentcontactforce.Update(mu*std::fabs(normalcontactforce),tangent);
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
    ),
    ddt_(Teuchos::null)
{
  return;
}


/*----------------------------------------------------------------------*
 | compute collisions (inter-particle and particle-wall)   ghamm 09/13  |
 *----------------------------------------------------------------------*/
double PARTICLE::ParticleCollisionHandlerMD::EvaluateParticleContact(
  double dt,
  Teuchos::RCP<Epetra_Vector> disn,
  Teuchos::RCP<Epetra_Vector> veln
  )
{
  // setup empty vector which contains current time between [0,dt) for each particle
  ddt_ = LINALG::CreateVector(*discret_->NodeColMap(), true);

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

#ifdef OUTPUT
        std::cout << "The eventqueue contains: ";
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

      }

      SearchForNewCollisions(next_event, eventqueue, dt);
    }

  }

  // updating of all particles to the end of the current time step
  for (int i=0; i<discret_->NumMyColNodes(); ++i)
  {
    for (int dim=0; dim<3; ++dim)
    {
      (*disncol_)[3*i + dim] += (*velncol_)[3*i + dim] * (dt - (*ddt_)[i]);
    }
  }

  // copy values from col to row layout
  for(int i=0; i<disn->MyLength(); ++i)
  {
    int gid = disn->Map().GID(i);
    int lid = disncol_->Map().LID(gid);
    (*disn)[i] = (*disncol_)[lid];
    (*veln)[i] = (*velncol_)[lid];
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
    std::set<DRT::Node*> neighboringparticles;
    std::set<DRT::Element*> neighboringwalls;
    GetNeighbouringParticlesAndWalls(currparticle, neighboringparticles, neighboringwalls);

    // loop over all neighbouring particles and check if the sum of their radii is larger than their distance
    for (std::set<DRT::Node*>::iterator iter=neighboringparticles.begin(); iter!=neighboringparticles.end(); ++iter)
    {
      if (currparticle->Id() > (*iter)->Id())
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
    for (std::set<DRT::Element*>::iterator iter=neighboringwalls.begin(); iter!=neighboringwalls.end(); ++iter)
    {
      Teuchos::RCP<DRT::Discretization> walldiscret = particle_algorithm_->WallDiscret();
      Teuchos::RCP<const Epetra_Vector> walldisnp = walldiscret->GetState("walldisnp");

      DRT::Node** nodes = (*iter)->Nodes();
      int numnodes = (*iter)->NumNode();

      std::vector<int> lm;
      lm.reserve(numnodes*3);

      std::vector<int> lmowner;
      std::vector<int> lmstride;
      (*iter)->LocationVector(*walldiscret, lm, lmowner, lmstride);

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
      GEO::nearest3DObjectOnElement(*iter, wallpositions, currposition, minDistCoords);
      static LINALG::Matrix<3,1> distance;
      distance.Update(1.0, currposition, -1.0, minDistCoords);

      if (distance.Norm2() < (currradius - GEO::TOL14))
      {
        std::cout << "particle " << currparticle->Id() << std::endl;
        std::cout << "wall element " << (*iter)->Id() << std::endl;
        std::cout << "distance " << distance.Norm2() << std::endl;
        std::cout << "currentradius " << currradius << std::endl;
        std::cout << "particle is penetrating the wall" << std::endl;
        dserror("Particle is penetrating the wall");
      }
    }
  }
#endif

  // reset vector because it needs to be rebuild in the next time step
  ddt_ = Teuchos::null;

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
    // get particle data
    DRT::Node* particle_1 = next_event->particle_1;
    DRT::Node* particle_2 = next_event->particle_2;

    std::vector<int> lm_1, lm_2;
    lm_1.reserve(3);
    lm_2.reserve(3);

    discret_->Dof(particle_1, lm_1);
    discret_->Dof(particle_2, lm_2);

    static LINALG::Matrix<3,1> pos_1, pos_2, vel_1, vel_2, pos_1_new, pos_2_new, vel_1_new, vel_2_new;
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*disncol_, pos_1, lm_1);
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*disncol_, pos_2, lm_2);
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*velncol_, vel_1, lm_1);
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*velncol_, vel_2, lm_2);

    int lid_1 = particle_1->LID();
    int lid_2 = particle_2->LID();

    // compute particle positions and collision normal at collision time
    static LINALG::Matrix<3,1> unitcollnormal;
    pos_1_new.Update(1.,pos_1,next_event->time-(*ddt_)[lid_1],vel_1);
    pos_2_new.Update(1.,pos_2,next_event->time-(*ddt_)[lid_2],vel_2);
    unitcollnormal.Update(1.,pos_2_new,-1.,pos_1_new);
    unitcollnormal.Scale(1./unitcollnormal.Norm2());

    // compute velocities of particles in normal direction
    const double veln1(vel_1.Dot(unitcollnormal));
    const double veln2(vel_2.Dot(unitcollnormal));

    // check for collision: normal velocity of particle_1 must be greater than of particle_2
    const double deltaveln = veln1 - veln2;
    if (deltaveln > GEO::TOL14)
    {
      // get masses
      const double mass_1 = (*masscol_)[lid_1];
      const double mass_2 = (*masscol_)[lid_2];
      const double mass = mass_1 + mass_2;

      // compute new velocities in normal direction
      const double veln1_new = veln1 - (1.0 + e_) * mass_2 * deltaveln / mass;
      const double veln2_new = veln2 + (1.0 + e_) * mass_1 * deltaveln / mass;

      // compute new velocities
      vel_1_new.Update(1.,vel_1,veln1_new-veln1,unitcollnormal);
      vel_2_new.Update(1.,vel_2,veln2_new-veln2,unitcollnormal);

      // update ddt
      (*ddt_)[lid_1] = next_event->time;
      (*ddt_)[lid_2] = next_event->time;

      for (int i=0; i<3; ++i)
      {
        lid_1 = disncol_->Map().LID(lm_1[i]);
        lid_2 = disncol_->Map().LID(lm_2[i]);
        // update particle positions
        (*disncol_)[lid_1] = pos_1_new(i);
        (*disncol_)[lid_2] = pos_2_new(i);
        // update particle velocities
        (*velncol_)[lid_1] = vel_1_new(i);
        (*velncol_)[lid_2] = vel_2_new(i);
      }
#ifdef OUTPUT
      std::cout << "New position of particle with GID " << particle_1->Id() << "  x: " << pos_1_new[0] << "  y: " << pos_1_new[1] << "  z: " << pos_1_new[2] << std::endl;
      std::cout << "New position of particle with GID " << particle_2->Id() << "  x: " << pos_2_new[0] << "  y: " << pos_2_new[1] << "  z: " << pos_2_new[2] << std::endl;
#endif
    }
  }
  break;
  case INPAR::PARTICLE::particle_wall:
  {
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
    if (std::fabs(normallength - radius) > GEO::TOL7)
    {
      std::cout << "ran into error :" << std::endl;
      std::cout << "particle pos: " << newpos << std::endl;
      std::cout << "wallcollpoint_pos: " << Teuchos::rcp_static_cast<WallEvent>(next_event)->wallcollpoint_pos << std::endl;
      std::cout << "normallength: " << normallength << std::endl;
      std::cout << "radius: " << radius << std::endl;
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
    (*ddt_)[lid] = next_event->time;

    std::vector<int> lm;
    lm.reserve(3);
    discret_->Dof(next_event->particle_1, lm);

    for (int i=0; i<3; ++i)
    {
      lid = disncol_->Map().LID(lm[i]);
      (*disncol_)[lid] = newpos(i);
      (*velncol_)[lid] = newvel(i);
    }
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

  if (particle_1->Id() == particle_2->Id())
    dserror("Particle %i cannot collide with itself!", particle_1->Id());

  static LINALG::Matrix<3,1> pos_1, pos_2, vel_1, vel_2;
  double rad_1, rad_2, ddt_1, ddt_2;
  if(particledata_.size() != 0)
  {
    const ParticleCollData& particle1 = particledata_[particle_1->LID()];
    pos_1 = particle1.pos;
    vel_1 = particle1.vel;
    rad_1 = particle1.rad;
    ddt_1 = particle1.ddt;
    const ParticleCollData& particle2 = particledata_[particle_2->LID()];
    pos_2 = particle2.pos;
    vel_2 = particle2.vel;
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
    if (std::fabs(tc1) <= GEO::TOL14 or std::fabs(tc2) <= GEO::TOL14)
    {
      // compute collision normal to detect whether collision has already happened at the end of the last time step
      static LINALG::Matrix<3,1> colnormal;
      colnormal.Update(1.0, pos_2, -1.0, pos_1);

      // velocities of particles in normal direction
      const double vel_col_1 =   vel_1.Dot(colnormal);
      const double vel_col_2 = - vel_2.Dot(colnormal);

      if ((vel_col_1 + vel_col_2) > GEO::TOL14)
      {
        // create event
        newevent = Teuchos::rcp(new Event(INPAR::PARTICLE::particle_particle, 0.0, particle_1, particle_2));
      }
    }
    // tc1 is negative
    else if (tc1 < -GEO::TOL14 and tc2 > GEO::TOL14 and tc2 < 1.1*remaining_dt)
    {
      // create event
      newevent = Teuchos::rcp(new Event(INPAR::PARTICLE::particle_particle, tc2, particle_1, particle_2));
    }
    // tc2 is negative
    else if (tc1 > GEO::TOL14 and tc2 < -GEO::TOL14 and tc1 < 1.1*remaining_dt)
    {
      newevent = Teuchos::rcp(new Event(INPAR::PARTICLE::particle_particle, tc1, particle_1, particle_2));
    }
    // both positive, smaller one is chosen (tc1 is almost identical to tc2)
    else if (tc1 > GEO::TOL14 and tc2 > GEO::TOL14 and std::min(tc1, tc2) < 1.1*remaining_dt)
    {
      newevent = Teuchos::rcp(new Event(INPAR::PARTICLE::particle_particle, std::min(tc1, tc2), particle_1, particle_2));
    }
  }

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

  if(particledata_.size() != 0)
  {
    const ParticleCollData& particle_data = particledata_[particle->LID()];
    position = particle_data.pos;
    velocity = particle_data.vel;
    radius = particle_data.rad;
    particle_time = particle_data.ddt;
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
    wallevent = Teuchos::rcp(new WallEvent(INPAR::PARTICLE::particle_wall, timetocollision, particle, wall, wallcoll_pos, wallcoll_vel));

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
  TEUCHOS_FUNC_TIME_MONITOR("InitializingEventqueue");


  // gather data for all particles initially
  const int numparticles = discret_->NodeColMap()->NumMyElements();
  particledata_.resize(numparticles);
  {
    for (int i=0; i<numparticles; ++i)
    {
      // particle for which data will be collected
      DRT::Node *particle = discret_->lColNode(i);

      static LINALG::Matrix<3,1> position;
      static LINALG::Matrix<3,1> velocity;
      double radius;
      double particle_time;

      GetCollisionData(particle, position, velocity, radius, particle_time);

      ParticleCollData& data = particledata_[particle->LID()];
      data.pos = position;
      data.vel = velocity;
      data.rad = radius;
      data.ddt = particle_time;
    }
  }

  // setup of initial event queue
  std::set<int> examinedbins;
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

    // now all particles in this bin are processed
    std::set<DRT::Node*> neighbouring_particles;
    std::set<DRT::Element*> neighbouring_walls;

    // gather all neighboring particles and wall elements
    GetNeighbouringParticlesAndWalls(currparticle, neighbouring_particles, neighbouring_walls);
    // add currparticle to neighbor list again
    neighbouring_particles.insert(currparticle);

    DRT::Node** particles = currbin->Nodes();
    for(int iparticle=0; iparticle<currbin->NumNode(); ++iparticle)
    {
      DRT::Node* currnode = particles[iparticle];
      // remove current node from neighbor list
      neighbouring_particles.erase(currnode);

      // loop over all neighbouring particles and check if they collide with currnode
      for (std::set<DRT::Node*>::iterator iter=neighbouring_particles.begin(); iter!=neighbouring_particles.end(); ++iter)
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

      // add current node to neighbor list again
      neighbouring_particles.insert(currnode);

      // loop over all neighbouring wall elements and check if they collide with currparticle
      for (std::set<DRT::Element*>::iterator iter=neighbouring_walls.begin(); iter!=neighbouring_walls.end(); ++iter)
      {
        Teuchos::RCP<WallEvent> newevent = ComputeCollisionWithWall(currnode, *iter, dt);

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

  // important to erase initial data here as content is no longer valid during collisions
  particledata_.clear();

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
  TEUCHOS_FUNC_TIME_MONITOR("UpdatdingEventqueue");

  std::set<DRT::Node*> neighbouring_particles;
  std::set<DRT::Element*> neighbouring_walls;

  //
  // searching for new collisions for the first particle of the last collision
  //
  GetNeighbouringParticlesAndWalls(lastevent->particle_1, neighbouring_particles, neighbouring_walls);

  // particle-particle collision
  for (std::set<DRT::Node*>::iterator iter=neighbouring_particles.begin(); iter!=neighbouring_particles.end(); ++iter)
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
          // only time to collision is returned --> time of the last event needs to be added
          newevent->time += lastevent->time;
          eventqueue.insert(newevent);
        }
      }
      else
      {
        // only time to collision is returned --> time of the last event needs to be added
        newevent->time += lastevent->time;
        eventqueue.insert(newevent);
      }
    }
  }

  // particle-wall collision
  for(std::set<DRT::Element*>::iterator iter=neighbouring_walls.begin(); iter!=neighbouring_walls.end(); ++iter)
  {
    Teuchos::RCP<WallEvent> newevent = ComputeCollisionWithWall(lastevent->particle_1, *iter, dt);

    if (newevent != Teuchos::null)
    {
      // only time to collision is returned --> time of the last event needs to be added
      newevent->time += lastevent->time;
      eventqueue.insert(newevent);
    }
  }

  //
  // searching for new collisions for the second particle of the last collision (if existing)
  //
  if (lastevent->coltype == INPAR::PARTICLE::particle_particle)
  {
    GetNeighbouringParticlesAndWalls(lastevent->particle_2, neighbouring_particles, neighbouring_walls);

    // particle-particle collision
    for(std::set<DRT::Node*>::iterator iter=neighbouring_particles.begin(); iter!=neighbouring_particles.end(); ++iter)
    {
      // IMPORTANT: first particle must be the particle, that collided in this time step
      Teuchos::RCP<Event> newevent = ComputeCollisionWithParticle(lastevent->particle_2, *iter, dt-lastevent->time);

      // do not add event of particles that has already been processed in lastevent
      if (newevent != Teuchos::null and newevent->particle_2->Id() != lastevent->particle_1->Id())
      {
        // only time to collision is returned --> time of the last event needs to be added
        newevent->time += lastevent->time;
        eventqueue.insert(newevent);
      }
    }

    // particle-wall collision
    for(std::set<DRT::Element*>::iterator iter=neighbouring_walls.begin(); iter!=neighbouring_walls.end(); ++iter)
    {
      Teuchos::RCP<WallEvent> newevent = ComputeCollisionWithWall(lastevent->particle_2, *iter, dt);

      if (newevent != Teuchos::null)
      {
        // only time to collision is returned --> time of the last event needs to be added
        newevent->time += lastevent->time;
        eventqueue.insert(newevent);
      }
    }
  }

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
    double timetocollision_line;
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
    double timetocollision_corner;
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
    const double det = LINALG::gaussElimination<true,3>(A, b, dx);

    if (std::fabs(det) < GEO::TOL14)
    {
      if(unitnormal.Dot(particle_vel_fad) < GEO::TOL10 and iter>1)
      {
//        std::cout << "particle path is parallel to wall --> left iteration" << std::endl;
        return validcollision;
      }
    }

    // update of coll_solution
    for (int i=0; i<3; ++i)
    {
      coll_solution(i) += dx(i);
    }

    // in case of parallel movement of particle to wall, element coords get extremely large
    if (std::fabs(coll_solution(0)) > 1.0e3 or std::fabs(coll_solution(1)) > 1.0e3)
    {
//      std::cout << "elecoord is extremely large --> left iteration" << std::endl;
      return validcollision;
    }
  }

  // check whether collision position is valid and otherwise check for edge collision points
  if (GEO::checkPositionWithinElementParameterSpace(coll_solution, DISTYPE) == true)
  {
    // check whether collision time is reasonable
    if (coll_solution(2) >= -GEO::TOL14 and coll_solution(2) < 1.1 * remaining_dt)
    {
      FAD scalar_partvel = unitnormal.Dot(particle_vel_fad);
      FAD scalar_wallvel = unitnormal.Dot(wallcollpoint_vel_fad);
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
      }
    }
  }
  else
  {
    checkedges = true;
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
    const double det = LINALG::gaussElimination<true,2>(A_line, b_line, dx_line);

    if (std::fabs(det) < GEO::TOL14)
    {
      static LINALG::Matrix<3,1> crossproduct;
      crossproduct.CrossProduct(particle_vel, F_deriv1);
      if(crossproduct.Norm1() < GEO::TOL10 and iter>1)
      {
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
    const double residualnorm1 = std::fabs(residual(0).val()) + std::fabs(residual(1).val());

#ifdef OUTPUT
    std::cout << "iteration: " << iter << "residualnorm: " << residualnorm1 << " distance: " << distance <<
        " time to collision: " << coll_solution(1).val() << " remaining dt: " << remaining_dt << std::endl;
#endif

    // break if collision found or residual is small and collision not possible
    if (distance<GEO::TOL12 || (residualnorm1<1.0e-12*radius*radius && distance>GEO::TOL2*radius))
    {
      break;
    }

    // break here if there will be definitely no contact using some heuristic
    if (iter>2 && ((std::fabs(distance-distance_old)<0.1*radius && distance>radius) || std::fabs(coll_solution(0).val())>1.5) )
    {
#ifdef OUTPUT
    std::cout << "break due to far away" << std::endl;
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
    std::cout << "break due to reached tolerance and medium far away" << std::endl;
#endif
#ifndef DEBUG
      break;
#else
      heuristic_break = true;
#endif
    }
    const double remainingtime = std::max(2.0*remaining_dt,GEO::TOL2*dt);
#ifdef OUTPUT
  std::cout << "remainingtime: " << remainingtime << " coll_solution(1).val(): " << coll_solution(1).val() << " coll_solution(1).val()-ttc_old: " << (coll_solution(1).val()-ttc_old) << std::endl;
#endif
    if (iter>2 && (coll_solution(1).val()>remainingtime && std::fabs(coll_solution(1).val()-ttc_old)<0.2*remainingtime))
    {
#ifdef OUTPUT
    std::cout << "break due to time to collision too large" << std::endl;
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
    const double det = LINALG::gaussElimination<true,2>(A, b, dx);

    if (std::fabs(det) < GEO::TOL14)
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
    if (std::fabs(coll_solution(0)) > 1.0e3)
    {
//      std::cout << "elecoord is extremely large --> left iteration" << std::endl;
      return validlinecollision;
    }
  }

  // check whether collision position is valid and otherwise check for corner collision points
  if (GEO::checkPositionWithinElementParameterSpace(coll_solution, DISTYPE) == true and distance < GEO::TOL7*radius)
  {
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
  LINALG::Matrix<3,1>& pos,
  LINALG::Matrix<3,1>& vel,
  double& rad,
  double& ddt
  )
{
  // data of particle 1
  std::vector<int> lm;
  lm.reserve(3);
  discret_->Dof(particle, lm);
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*disncol_, pos, lm);
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<3,1> >(*velncol_, vel, lm);
  const int lid = particle->LID();
  rad = (*radiusncol_)[lid];
  ddt = (*ddt_)[lid];

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
    if (std::fabs(event1->time - event2->time) < GEO::TOL14)
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
    if (std::fabs(event1->time - event2->time) < GEO::TOL14)
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
    if (std::fabs(event1->time - event2->time) < GEO::TOL14)
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
    if (std::fabs(event1->time - event2->time) < GEO::TOL14)
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

