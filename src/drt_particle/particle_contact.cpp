/*----------------------------------------------------------------------*/
/*!
\file particle_contact.cpp
\brief Particle collision handling

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
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
#include "../drt_geometry/element_normals.H"

#include "../drt_mat/stvenantkirchhoff.H"

#include <Teuchos_TimeMonitor.hpp>
#include <Epetra_FEVector.h>


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
  int binId = CurrentBin[0]->Id();

  int ijk[3];
  particle_algorithm_->ConvertGidToijk(binId,ijk);

  // ijk_range contains: i_min   i_max     j_min     j_max    k_min     k_max
  int ijk_range[] = {ijk[0]-1, ijk[0]+1, ijk[1]-1, ijk[1]+1, ijk[2]-1, ijk[2]+1};
  std::set<int> binIds;

  particle_algorithm_->GidsInijkRange(ijk_range,binIds,true);

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
  std::set<int> &binIds
  )
{
  // loop over all bins
  for(std::set<int>::const_iterator bin=binIds.begin(); bin!=binIds.end(); ++bin)
  {
    DRT::Element *neighboringbin = discret_->gElement(*bin);

    // gather wall elements
    DRT::Element** walleles = static_cast<DRT::MESHFREE::MeshfreeMultiBin*>(neighboringbin)->AssociatedWallEles();
    int numwalls = static_cast<DRT::MESHFREE::MeshfreeMultiBin*>(neighboringbin)->NumAssociatedWallEle();
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
  //extract input-parameters
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

    if(r_min_<0.0 or r_max_<0.0 or v_max_<0.0 or c_<0.0)
      dserror("Invalid input parameter (MIN_RADIUS,MAX_RADIUS,MAX_VELOCITY,REL_PENETRATION have to be larger than zero)");

    if((e_<0.0 or e_wall_<0.0) && normal_contact_ == INPAR::PARTICLE::LinSpringDamp)
      dserror("Invalid input parameter COEFF_RESTITUTION for this kind of contact law");

    //critical time step
    double mass_min = 0.0;
    mass_min = density * 4.0/3.0 * M_PI * pow( r_min_ ,3.0 );

    double k_tkrit = 0.0;

    // here the first element of the structural problem is asked for its material and assumed that wall element has same properties
    Teuchos::RCP<DRT::Discretization> structdis = particle_algorithm_->Structure()->Discretization();
    int local_structmatid = -1;
    // in order to allow procs without wall element communication is necessary
    if(structdis->NumMyColElements() != 0)
    {
      Teuchos::RCP<MAT::Material> mat = structdis->lColElement(0)->Material();
      MAT::PAR::Parameter* params = mat->Parameter();
      local_structmatid = params->Id();
    }
    int structmatid = -1;
    structdis->Comm().MaxAll(&local_structmatid, &structmatid, 1);

    Teuchos::RCP<MAT::StVenantKirchhoff> structmat = Teuchos::rcp_dynamic_cast<MAT::StVenantKirchhoff>(MAT::Material::Factory(structmatid));
    if(structmat == Teuchos::null)
      dserror("only stvenantkirchhoff material is supported so far :-(");

    double G_wall = structmat->ShearMod();
    double nue_wall = structmat->PoissonRatio();

    double G= young_ / (2*(1+nue_));

    //kappa - tangential to normal stiffness ratio
    kappa_ = (1-nue_)/(1-0.5*nue_);
    kappa_wall_ = ( (1-nue_)/G + (1-nue_wall)/G_wall ) / ( (1-0.5*nue_)/G + (1-0.5*nue_wall)/G_wall );

    //------------stiffness----------------------------
    switch(normal_contact_)
    {
    case INPAR::PARTICLE::LinSpring:
    case INPAR::PARTICLE::LinSpringDamp:
    {
      //stiffness calculated from relative penetration and some other input parameters (linear spring)
      k_normal_ = 2.0/3.0 * r_max_ * M_PI * density * pow(v_max_,2.0) / pow(c_,2.0);
      //for tangential contact same stiffness is used
      k_tang_ = kappa_ * k_normal_;
      k_tang_wall_ = kappa_wall_ * k_normal_;
      k_tkrit = k_normal_;

      double user_normal_stiffness = particleparams.get<double>("NORMAL_STIFF");
      if(user_normal_stiffness > 0.0)
      {
        //if user specifies normal stiffness, this stiffness will be used as normal and tangential stiffness for simulation
        k_normal_ = user_normal_stiffness;
        //for tangential contact same stiffness is used
        k_tang_ =kappa_ * k_normal_;
        k_tang_wall_ =kappa_wall_ * k_normal_;
        //stiffness used for calculation of critical time step
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

      //stiffness calculated from relative penetration and some other input parameters (Hertz)
      k_normal_ = 10.0/3.0 * M_PI * density * pow(v_max_,2.0) * pow(r_max_,0.5) / pow(2*c_,2.5);
      //stiffness used for calculation of critical time step (linear spring stiffness needed!)
      k_tkrit = 2.0/3.0 * r_max_ * M_PI * density * pow(v_max_,2.0) / pow(c_,2.0);

      double user_normal_stiffness = particleparams.get<double>("NORMAL_STIFF");
      if(user_normal_stiffness > 0.0)
      {
        //if user specifies normal stiffness, this stiffness will be used as normal stiffness for simulation
        k_normal_ = user_normal_stiffness;
        //for tangential contact the user specified (nonlinear) normal stiffness which has to be transformed into a linear normal
        //stiffness with the same relative penetration which is used as (linear) tangential stiffness afterwards
        double value = 2048.0/1875.0 * density * pow(v_max_,2.0) * M_PI * pow(r_max_,3.0) * pow(k_normal_,4.0);
        //stiffness used for calculation of critical time step (linear spring stiffness needed!)
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
        //user has to specify normal damping coefficient
        d_normal_ = user_normal_damping;
      }
      else
      {
        dserror("For this kind of contact law the input parameter NORMAL_DAMP is invalid");
      }
      double user_tang_damping = particleparams.get<double>("TANG_DAMP");
      if(user_tang_damping >= 0.0)
      {
        //user has to specify tangential damping coefficient
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
        //user has to specify normal damping coefficient
        d_normal_wall_ = user_normal_damping;
      }
      else
      {
        dserror("For this kind of contact law the input parameter NORMAL_DAMP_WALL is invalid");
      }
      double user_tang_damping = particleparams.get<double>("TANG_DAMP_WALL");
      if(user_tang_damping >= 0.0)
      {
        //user has to specify tangential damping coefficient
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

    //initialize factor
    if(contact_strategy_==INPAR::PARTICLE::Normal_DEM)
    { factor = 0.34; }
    if(contact_strategy_==INPAR::PARTICLE::NormalAndTang_DEM)
    { factor = 0.22; }

    //calculate critical time step
    dt_krit_ = safety * factor * sqrt( mass_min / k_tkrit );

    //check frictional coefficient
    if(contact_strategy_ == INPAR::PARTICLE::NormalAndTang_DEM)
    {
      if(mu_<=0.0 or mu_wall_ <=0.0)
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
  double dt,
  Teuchos::RCP<Epetra_Vector> f_contact,
  Teuchos::RCP<Epetra_Vector> m_contact
  )
{
  TEUCHOS_FUNC_TIME_MONITOR("TimeforContactSearchAndCalculation");

  if(dt>dt_krit_ and myrank_==0)
  {
    std::cout << "W A R N I N G : time step larger than critical time step!" << std::endl;
    std::cout << "W A R N I N G : calculated critical time step: "<<dt_krit_<<" !" << std::endl;
  }

  contact_energy_ = 0.0;

  // get wall discretization and states for particles
  Teuchos::RCP<DRT::Discretization> walldiscret = particle_algorithm_->WallDiscret();
  Teuchos::RCP<const Epetra_Vector> walldisn = walldiscret->GetState("walldisnp");
  Teuchos::RCP<const Epetra_Vector> wallveln = walldiscret->GetState("wallvelnp");

  // define vector for contact force
  Teuchos::RCP<Epetra_FEVector> f_structure = Teuchos::rcp(new Epetra_FEVector(*discret_->DofRowMap()));

  // store bins, which have already been examined
  std::set<int> examinedbins;

  // loop over all particles
  for(int i=0; i<discret_->NodeColMap()->NumMyElements(); ++i)
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
    int numparticle = CurrentBin[0]->NumNode();

    // loop over all particles in CurrentBin
    for(int i=0; i<numparticle; ++i)
    {
      DRT::Node *particle_i = NodesInCurrentBin[i];

      std::vector<int> lm_i;
      lm_i.reserve(3);

      // extract global dof ids and fill into lm_i
      discret_->Dof(particle_i, lm_i);
      int owner_i = particle_i->Owner();

      //position of particle i
      std::vector<double> myposition_i(3);
      DRT::UTILS::ExtractMyValues(*disncol_,myposition_i,lm_i);

      //velocity of particle i
      std::vector<double> myvel_i(3);
      DRT::UTILS::ExtractMyValues(*velncol_,myvel_i,lm_i);

      //angular-velocity of particle i
      std::vector<double> myangvel_i(3);
      DRT::UTILS::ExtractMyValues(*ang_velncol_,myangvel_i,lm_i);

      int lid = discret_->NodeColMap()->LID(particle_i->Id());

      // radius of particle i
      double radius_i = (*radiusncol_)[lid];

      // mass of particle i
      double mass_i = (*masscol_)[lid];

      // evaluate contact with walls first
      std::map<int,WallContactPoint> surfaces;
      std::map<int,WallContactPoint> lines;
      std::map<int,WallContactPoint> nodes;

      std::set<int> unusedIds;

      // check whether there is contact between particle i and neighboring walls
      for(std::set<DRT::Element*>::const_iterator w=neighboringwalls.begin(); w!=neighboringwalls.end();  ++w)
      {
        int numnodes = (*w)->NumNode();
        std::vector<int> lm_wall;
        lm_wall.reserve(numnodes * 3);

        std::vector<int> lmowner;
        std::vector<int> lmstride;
        (*w)->LocationVector(*walldiscret,lm_wall,lmowner,lmstride);

        // nodal displacements
        std::vector<double> nodal_disp(numnodes * 3);
        DRT::UTILS::ExtractMyValues(*walldisn,nodal_disp,lm_wall);

        // get current position of nodes: x = X + u
        std::map<int,LINALG::Matrix<3,1> > nodeCoord;
        DRT::Node** wallnodes = (*w)->Nodes();
        for(int counter=0; counter<numnodes; ++counter)
        {
          LINALG::Matrix<3,1> currpos;
          const double* X = wallnodes[counter]->X();
          currpos(0) = X[0] + nodal_disp[counter*3+0];
          currpos(1) = X[1] + nodal_disp[counter*3+1];
          currpos(2) = X[2] + nodal_disp[counter*3+2];
          nodeCoord[wallnodes[counter]->Id()] = currpos;
        }

        LINALG::Matrix<3,1> nearestPoint;
        LINALG::Matrix<3,1> position_i;

        //transfer entries from myposition_i to position_i
        for(int n=0; n<3; ++n)
        {
          position_i(n)=myposition_i[n];
        }

        //-------find point on wall element with smallest distance to particle_i-------------------
        GEO::ObjectType objecttype = GEO::nearest3DObjectOnElement((*w),nodeCoord,position_i,nearestPoint);
        //-----------------------------------------------------------------------------------------

        LINALG::Matrix<3,1> r_i_wall;
        r_i_wall.Update(1.0, nearestPoint, -1.0, position_i);
        double distance_i_wall = r_i_wall.Norm2();
        double penetration = distance_i_wall-radius_i;

        if(penetration <= 0.0)
        {
          // get pointer to the current object type of closest point
          std::map<int,WallContactPoint> *pointer=0;
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
          for(std::map<int,WallContactPoint>::const_iterator iter = (*pointer).begin(); iter != (*pointer).end();++iter)
          {
            LINALG::Matrix<3,1> distance_vector;
            distance_vector.Update(1.0, nearestPoint, -1.0, (iter->second).point);
            double distance = distance_vector.Norm2();
            double adaptedtol = GEO::TOL7 * radius_i;

            if (distance < adaptedtol)
            {
              // point has already been detected --> do not insert
              insert = false;
              unusedIds.insert((*w)->Id());
              break;
            }
          }

          // insert contact point with current surface in corresponding map (surf, line, node)
          if(insert)
          {
            WallContactPoint currentContact = { nearestPoint, penetration , nodeCoord ,lm_wall , lmowner };
            (*pointer).insert(std::pair<int,WallContactPoint>((*w)->Id(),currentContact));
          }
        }
        // penetration > 0.0 --> contact impossible
        else
        {
          unusedIds.insert((*w)->Id());
        }
      }

      // find entries of lines and nodes which are within the penetration volume of the current particle
      // hierarchical: surfaces first
      for(std::map<int,WallContactPoint>::const_iterator surfiter = surfaces.begin(); surfiter != surfaces.end(); ++surfiter)
      {
        // within this radius no other contact point can lie: radius = sqrt(r_i^2 - (r_i-|g|)^2)
        double radius_surface = sqrt(pow(radius_i,2.0) - pow(radius_i-fabs((surfiter->second).penetration),2.0));

        for(std::map<int,WallContactPoint>::const_iterator lineiter = lines.begin(); lineiter != lines.end();++lineiter)
        {
          LINALG::Matrix<3,1> distance_vector;
          distance_vector.Update(1.0, (surfiter->second).point, -1.0, (lineiter->second).point);
          double distance = distance_vector.Norm2();
          if(distance <= radius_surface)
            unusedIds.insert(lineiter->first);
        }
        for(std::map<int,WallContactPoint>::const_iterator nodeiter=nodes.begin(); nodeiter!= nodes.end(); ++nodeiter)
        {
          LINALG::Matrix<3,1> distance_vector;
          distance_vector.Update(1.0, (surfiter->second).point, -1.0, (nodeiter->second).point);
          double distance = distance_vector.Norm2();
          if(distance <= radius_surface)
            unusedIds.insert(nodeiter->first);
        }
      }
      // find entries of nodes which are within the penetration volume of the current particle
      // hierarchical: lines next
      for(std::map<int,WallContactPoint>::const_iterator lineiter=lines.begin(); lineiter != lines.end(); ++lineiter)
      {
        // radius = sqrt(r_i^2 - (r_i-|g|)^2)
        double radius_line = sqrt(pow(radius_i,2.0) - pow(radius_i-fabs((lineiter->second).penetration),2.0));

        for(std::map<int,WallContactPoint>::const_iterator nodeiter=nodes.begin(); nodeiter!=nodes.end(); ++nodeiter)
        {
          LINALG::Matrix<3,1> distance_vector;
          distance_vector.Update(1.0, (lineiter->second).point, -1.0, (nodeiter->second).point);
          double distance = distance_vector.Norm2();
          if(distance <= radius_line)
            unusedIds.insert(nodeiter->first);
        }
      }

      // write entries of lines and nodes to surfaces if contact has to be evaluated
      for(std::map<int,WallContactPoint>::const_iterator iter = lines.begin(); iter != lines.end() ;++iter)
      {
        if( !unusedIds.count(iter->first) )
          surfaces.insert(std::pair<int,WallContactPoint>(iter->first,iter->second));
      }
      for(std::map<int,WallContactPoint>::const_iterator iter=nodes.begin(); iter!=nodes.end(); ++iter)
      {
        if( !unusedIds.count(iter->first) )
          surfaces.insert(std::pair<int,WallContactPoint>(iter->first,iter->second));
      }

      // evaluate contact between particle_i and entries of surfaces
      std::map<int, PARTICLE::Collision>& history_wall = static_cast<PARTICLE::ParticleNode*>(particle_i)->Get_history_wall();
      if(history_wall.size() > 3)
        dserror("Contact with more than 3 wall elements. Check whether history is deleted correctly.");

      for(std::map<int,WallContactPoint>::const_iterator iter = surfaces.begin(); iter != surfaces.end(); ++iter)
      {
        // gid of wall element
        int gid_wall = iter->first;

        // distance-vector
        double r_contact[3] = {0.0};

        // normal-vector
        double normal[3] = {0.0};

        // distance between centre of mass of two particles
        double norm_r_contact = 0.0;

        // velocity v_rel = v_i - v_j
        double v_rel[3] = {0.0};

        // velocity v_rel_tangential
        double v_rel_tangential[3] = {0.0};

        // part of v_rel in normal-direction
        double v_rel_normal = 0.0;

        // penetration
        double g = 0.0;

        // normalised mass
        double m_eff = mass_i;

        // contact force
        double normalcontactforce = 0.0;
        double tangentcontactforce[3] = {0.0};

        //distance-vector and distance--------------------------
        for(int n=0; n<3; ++n)
        {
          //calculate entries of r_contact
          r_contact[n] = (iter->second).point(n) - myposition_i[n];
          //length of r_contact
          norm_r_contact += r_contact[n]*r_contact[n];
        }
        norm_r_contact = sqrt(norm_r_contact);
        //--------------------------------------------------------

        //normal-vector-------------------------------------------
        for(int n=0; n<3; ++n)
        {
          normal[n] = r_contact[n]/norm_r_contact;
        }
        //-------------------------------------------------------

        // penetration--------------------------------------------
        // g = norm_r_contact - radius_i;
        g = (iter->second).penetration;

        if(fabs(g)>g_max_)
          g_max_ = fabs(g);
        //-------------------------------------------------------

        //-------get velocity of contact point-----------------------
        LINALG::Matrix<3,1> vel_nearestPoint(true);
        LINALG::Matrix<2,1> elecoord(true);
        DRT::Element *CurrentEle = walldiscret->gElement(gid_wall);
        const LINALG::SerialDenseMatrix xyze(GEO::getCurrentNodalPositions(CurrentEle, (iter->second).nodalCoordinates));

        // get coordinates of the projection point in parameter space of the element (xi_coordinates)
        GEO::CurrentToSurfaceElementCoordinates(CurrentEle->Shape(), xyze, (iter->second).point, elecoord);

        int numnodes = CurrentEle->NumNode();
        Epetra_SerialDenseVector funct(numnodes);

        // get shape functions of the element evaluated at the projection point
        DRT::UTILS::shape_function_2D(funct,elecoord(0,0),elecoord(1,0),CurrentEle->Shape());

        std::vector<double> nodal_vel(numnodes * 3);
        DRT::UTILS::ExtractMyValues(*wallveln,nodal_vel,(iter->second).lm);
        for(int node=0; node<numnodes; ++node)
        {
          for(int dim=0; dim<3; ++dim)
          {
            vel_nearestPoint(dim) += funct[node] * nodal_vel[node * 3 + dim];
          }
        }
        //-----------------------------------------------------------

        // velocity v_rel = v_i - v_wall and v_rel in normal-direction: v_rel * n
        for(int n=0; n<3; ++n)
        {
          v_rel[n] = myvel_i[n] - vel_nearestPoint(n);
          v_rel_normal += v_rel[n]*normal[n];
        }

        // normal contact force between particle and wall (note: owner_j = -1)
        CalculateNormalContactForce(g, v_rel_normal, mass_i, normalcontactforce, owner_i, -1);

        if(contact_strategy_ == INPAR::PARTICLE::NormalAndTang_DEM)
        {
          // velocity v_rel = v_i + omega_i x (r'_i n) - v_wall = (v_i - v_wall) + omega_i x (r'_i n)
          // and velocity v_rel_tangential
          for(int n=0; n<3; ++n)
          {
            v_rel[n] += (radius_i+g) * (myangvel_i[(n+1)%3] * normal[(n+2)%3] - myangvel_i[(n+2)%3] * normal[(n+1)%3]);
            v_rel_tangential[n] = v_rel[n] - v_rel_normal * normal[n];
          }

          // if g < 0 and g_lasttimestep > 0 -> create history variables
          if(!history_wall.count(gid_wall))
          {
             PARTICLE::Collision col;
             // initialize with stick
             col.stick = true;
             // initialize g_t[3]
             for(int n=0; n<3; ++n)
             {
               col.g_t[n] = 0.0;
             }

             // insert new entry
             history_wall.insert(std::pair<int,PARTICLE::Collision>(gid_wall,col));
          }

          // calculation of tangential contact force
          CalculateTangentialContactForce(normalcontactforce, normal, tangentcontactforce,
                      history_wall[gid_wall], v_rel_tangential, m_eff, dt,owner_i, -1);
        }

        // assembly of contact forces and moments
        Epetra_SerialDenseVector val_i(3);
        Epetra_SerialDenseVector m_i(3);
        std::vector<int> lmowner_i(3);

        double r_i = radius_i + g;
        int owner_i = particle_i->Owner();

        for(int n=0; n<3; ++n)
        {
          lmowner_i[n] = owner_i;
          // contact forces
          val_i[n] = normalcontactforce * normal[n] + tangentcontactforce[n];
          // moments: m_i = (r_i * n) x F_t
          m_i[n] = r_i * (normal[(n+1)%3] * tangentcontactforce[(n+2)%3] - normal[(n+2)%3] * tangentcontactforce[(n+1)%3]);
        }

        // do assembly of contact moments
        LINALG::Assemble(*m_contact,m_i,lm_i,lmowner_i);

        // do assembly of contact forces
        LINALG::Assemble(*f_contact,val_i,lm_i,lmowner_i);

        // forces on wall elements
        double nodal_forces[numnodes * 3];
        for(int node=0; node<numnodes; ++node)
        {
          for(int n=0; n<3; ++n)
          {
            nodal_forces[node * 3 + n] = funct[node] *(- val_i[n]);
          }
        }

        // assembly of contact forces on walls
        if(owner_i == myrank_)
        {
          int err = f_structure->SumIntoGlobalValues(numnodes * 3, &((iter->second).lm)[0], &nodal_forces[0]);
          if (err<0)
            dserror("summing into Epetra_FEVector failed");
        }
      } // end for contact points on surfaces

      if(contact_strategy_ == INPAR::PARTICLE::NormalAndTang_DEM)
      {
        //delete those entries in history_wall_ which are no longer in contact with particle_i in current time step
        for(std::set<DRT::Element*>::const_iterator w=neighboringwalls.begin(); w != neighboringwalls.end(); ++w)
        {
          int gid_wall = (*w)->Id();
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
        int gid_j = (*j)->Id();
        //evaluate contact only once!
        if(particle_i->Id() >= gid_j)
          continue;

        int owner_j = (*j)->Owner();
        std::vector<int> lm_j;
        lm_j.reserve(3);

        // extract global dof ids
        discret_->Dof((*j), lm_j);

        // position of particle j
        std::vector<double> myposition_j(3);
        DRT::UTILS::ExtractMyValues(*disncol_,myposition_j,lm_j);

        // velocity of particle j
        std::vector<double> myvel_j(3);
        DRT::UTILS::ExtractMyValues(*velncol_,myvel_j,lm_j);

        // angular velocity of particle j
        std::vector<double> myangvel_j(3);
        DRT::UTILS::ExtractMyValues(*ang_velncol_,myangvel_j,lm_j);

        lid = discret_->NodeColMap()->LID(gid_j);

        // radius of particle j
        double radius_j = (*radiusncol_)[lid];

        // mass of particle j
        double mass_j = (*masscol_)[lid];

        // normalized mass
        double m_eff = mass_i * mass_j / (mass_i + mass_j);

        // contact force
        double normalcontactforce = 0.0;
        double tangentcontactforce[3] = {0.0};

        // distance vector and distance between two particles
        double r_contact[3];
        double norm_r_contact = 0.0;
        for(int n=0; n<3; ++n)
        {
          r_contact[n] = myposition_j[n]-myposition_i[n];
          norm_r_contact += r_contact[n]*r_contact[n];
        }
        norm_r_contact = sqrt(norm_r_contact);

        // penetration
        double g = norm_r_contact - radius_i - radius_j;
        // in case of penetration contact forces and moments are calculated
        if(g <= 0.0)
        {
          if(fabs(g)>g_max_)
            g_max_ = fabs(g);

          // normal vector and velocity v_rel = v_i - v_j and part of v_rel in normal-direction: v_rel * n
          // velocity v_rel
          double v_rel[3];
          // normal vector
          double normal[3];
          // part of v_rel in normal-direction
          double v_rel_normal = 0.0;
          for(int n=0; n<3; ++n)
          {
            normal[n] = r_contact[n]/norm_r_contact;
            v_rel[n] = myvel_i[n] - myvel_j[n];
            v_rel_normal += v_rel[n]*normal[n];
          }

          // calculation of normal contact force
          CalculateNormalContactForce(g, v_rel_normal, m_eff, normalcontactforce, owner_i, owner_j);

          // calculation of tangential contact force
          if(contact_strategy_ == INPAR::PARTICLE::NormalAndTang_DEM)
          {
            // velocity v_rel = v_i - v_j + omega_i x (r'_i n) + omega_j x (r'_j n)
            // and velocity v_rel_tangential
            double v_rel_tangential[3];
            for(int n=0; n<3; ++n)
            {
              v_rel[n] += (radius_i+g/2.0) * (myangvel_i[(n+1)%3] * normal[(n+2)%3] - myangvel_i[(n+2)%3] * normal[(n+1)%3]) + (radius_j+g/2.0) * (myangvel_j[(n+1)%3] * normal[(n+2)%3] - myangvel_j[(n+2)%3] * normal[(n+1)%3]);
              v_rel_tangential[n] = v_rel[n] - v_rel_normal * normal[n];
            }

            // if history variables does not exist -> create it
            if(!history_particle.count(gid_j))
            {
              PARTICLE::Collision col;
              // initialize with stick
              col.stick = true;
              //initialize g_t[3]
              for(int n=0; n<3; ++n)
              {
                col.g_t[n] = 0.0;
              }

              //insert new entry
              history_particle.insert(std::pair<int,PARTICLE::Collision>(gid_j,col));
            }

            CalculateTangentialContactForce(normalcontactforce, normal, tangentcontactforce,
                      history_particle[gid_j], v_rel_tangential, m_eff, dt, owner_i, owner_j);
          }

          //----------ASSEMBLY---------------------------------------

          Epetra_SerialDenseVector val_i(3);
          Epetra_SerialDenseVector m_i(3);
          std::vector<int> lmowner_i(3);

          Epetra_SerialDenseVector val_j(3);
          Epetra_SerialDenseVector m_j(3);
          std::vector<int> lmowner_j(3);

          double r_i = radius_i + g/2.0;
          double r_j = radius_j + g/2.0;

          // forces
          for(int n=0; n<3; ++n)
          {
            lmowner_i[n] = owner_i;
            lmowner_j[n] = owner_j;

            val_i[n] = normalcontactforce * normal[n] + tangentcontactforce[n];
            // actio = reactio
            val_j[n] = - val_i[n];
          }

          // moments
          for(int n=0; n<3; ++n)
          {
            // m_i = (r_i * n) x F_t
            m_i[n] = r_i * (normal[(n+1)%3] * tangentcontactforce[(n+2)%3] - normal[(n+2)%3] * tangentcontactforce[(n+1)%3]);

            // m_j = r_j/r_i * m_i
            m_j[n] = r_j/r_i * m_i[n];
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

  // call global assemble for particle forces on walls
  int err = f_structure->GlobalAssemble(Add, false);
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
  double g,
  double v_rel_normal,
  double m_eff,
  double& normalcontactforce,
  int owner_i,
  int owner_j
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
      d = 2 * fabs(log(e_wall_)) * sqrt(k_normal_ * m_eff / (pow(log(e_wall_),2.0)+ pow(M_PI,2.0)));
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
	    d = 2 * fabs(log(e_)) * sqrt(k_normal_ * m_eff / (pow(log(e_),2.0)+ pow(M_PI,2.0)));
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
  double normalcontactforce,
  double *normal,
  double *tangentcontactforce,
  PARTICLE::Collision &currentColl,
  double *v_rel_tangential,
  double m_eff,
  const double dt,
  int owner_i,
  int owner_j
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
      d = 2 * fabs(log(e_wall_)) * sqrt(k_normal_ * m_eff / (pow(log(e_wall_),2.0)+ pow(M_PI,2.0)));
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
	    d = 2 * fabs(log(e_)) * sqrt(k_normal_ * m_eff / (pow(log(e_),2.0)+ pow(M_PI,2.0)));
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
		interime += normal[n] * currentColl.g_t[n];
	}
	old_length = sqrt(old_length);

	// projection of g_t onto current normal at time n+1
  double new_length = 0.0;
	for(int n=0; n<3; ++n)
	{
		currentColl.g_t[n] += - interime * normal[n];
		new_length += currentColl.g_t[n] * currentColl.g_t[n];
	}
	new_length = sqrt(new_length);

	// ensure that g_t has the same length as before projection
	// if almost no tangential spring elongation, neglect it
	if(new_length > 1.0E-14)
	{
		for(int n=0; n<3; ++n)
		{
			currentColl.g_t[n] = old_length/new_length * currentColl.g_t[n];
		}
	}

	// update of elastic tangential displacement if stick is true
  if(currentColl.stick == true)
  {
    for(int n=0; n<3; ++n)
    {
      currentColl.g_t[n] += v_rel_tangential[n] * dt;
    }
  }

	// calculate tangential test force
  // norm of tangential contact force
  double norm_f_t = 0.0;
	for(int n=0; n<3; ++n)
	{
		tangentcontactforce[n] = - k * currentColl.g_t[n] - d * v_rel_tangential[n];
		norm_f_t += tangentcontactforce[n] * tangentcontactforce[n];
	}
	norm_f_t = sqrt(norm_f_t);

	// Coulomb friction law

	// tangential contact force for "stick" - case----------------------
	if( norm_f_t <= (mu * fabs(normalcontactforce)) )
	{
	  currentColl.stick = true;
		//tangential contact force already calculated
	}
	else //"slip"-case
	{
	  currentColl.stick = false;
		//calculate tangent vector ( unit vector in (test-)tangentcontactforce-direction )
		double tangent[3];
		for(int n=0; n<3; ++n)
		{
			tangent[n] = tangentcontactforce[n] / norm_f_t;
		}

		// calculate tangent contact force and tangential displacements
		for(int n=0; n<3; ++n)
		{
			tangentcontactforce[n] = mu * fabs(normalcontactforce) * tangent[n];
			currentColl.g_t[n] = - 1/k * (tangentcontactforce[n] + d * v_rel_tangential[n]);
		}
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

    contact_energy_ += EnergyAssemble(owner_i,owner_j)* 0.5 * k * pow(new_length,2.0);
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

#ifdef DEBUG
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
      std::cout << " ErasingInvalidCollisions " << std::endl;

      // loop over event queue and erase invalid collisions
      for (std::set<Teuchos::RCP<Event>, Event::Helper>::iterator iter=eventqueue.begin(); iter!=eventqueue.end(); /*no ++iter*/)
      {
#ifdef DEBUG
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
#ifdef DEBUG
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
#ifdef DEBUG
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
#ifdef DEBUG
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
#ifdef DEBUG
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

#ifdef DEBUG
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

    LINALG::Matrix<3,1> currposition;
    LINALG::Matrix<3,1> currvelocity;
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

      LINALG::Matrix<3,1> iterposition;
      LINALG::Matrix<3,1> itervelocity;
      double iterradius;
      double itertime;

      GetCollisionData(*iter, iterposition, itervelocity, iterradius, itertime);

      LINALG::Matrix<3,1> distance;
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

        LINALG::Matrix<3,1> nodepos;
        for (int i=0; i<3; ++i)
        {
          nodepos(i) = node->X()[i] + nodaldisnp[3*j + i];
        }
        wallpositions[node->Id()] = nodepos;
      }

      LINALG::Matrix<3,1> dummyvec;
      double distance;
      GEO::getDistanceToSurface(*iter, wallpositions, currposition, dummyvec, distance);

      if (distance < (currradius - GEO::TOL14))
      {
        std::cout << "distance " << distance << std::endl;
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

    std::vector<double> pos_1(3), pos_2(3), vel_1(3), vel_2(3), pos_1_new(3), pos_2_new(3), vel_1_new(3), vel_2_new(3);
    DRT::UTILS::ExtractMyValues(*disncol_, pos_1, lm_1);
    DRT::UTILS::ExtractMyValues(*disncol_, pos_2, lm_2);
    DRT::UTILS::ExtractMyValues(*velncol_, vel_1, lm_1);
    DRT::UTILS::ExtractMyValues(*velncol_, vel_2, lm_2);

    int lid_1 = ddt_->Map().LID(particle_1->Id());
    int lid_2 = ddt_->Map().LID(particle_2->Id());

    // compute particle positions and collision normal at collision time
    LINALG::Matrix<3,1> unitcollnormal;
    for (int i=0; i<3; ++i)
    {
      pos_1_new[i] = pos_1[i] + (next_event->time - (*ddt_)[lid_1]) * vel_1[i];
      pos_2_new[i] = pos_2[i] + (next_event->time - (*ddt_)[lid_2]) * vel_2[i];
      unitcollnormal(i) = pos_2_new[i] - pos_1_new[i];
    }
    unitcollnormal.Scale(1.0/unitcollnormal.Norm2());

    // compute velocities of particles in normal direction
    double veln1 = vel_1[0] * unitcollnormal(0) + vel_1[1] * unitcollnormal(1) + vel_1[2] * unitcollnormal(2);
    double veln2 = vel_2[0] * unitcollnormal(0) + vel_2[1] * unitcollnormal(1) + vel_2[2] * unitcollnormal(2);

    // check for collision: normal velocity of particle_1 must be greater than of particle_2
    double deltaveln = veln1 - veln2;
    if (deltaveln > GEO::TOL14)
    {
      // get masses
      double mass_1 = (*masscol_)[lid_1];
      double mass_2 = (*masscol_)[lid_2];
      double mass = mass_1 + mass_2;

      // compute new velocities in normal direction
      double veln1_new = veln1 - (1.0 + e_) * mass_2 * deltaveln / mass;
      double veln2_new = veln2 + (1.0 + e_) * mass_1 * deltaveln / mass;

      // compute new velocities
      for(int i=0; i<3; ++i)
      {
        vel_1_new[i] = vel_1[i] + (veln1_new - veln1) * unitcollnormal(i);
        vel_2_new[i] = vel_2[i] + (veln2_new - veln2) * unitcollnormal(i);
      }

      // update ddt
      (*ddt_)[lid_1] = next_event->time;
      (*ddt_)[lid_2] = next_event->time;

      for (int i=0; i<3; ++i)
      {
        lid_1 = disncol_->Map().LID(lm_1[i]);
        lid_2 = disncol_->Map().LID(lm_2[i]);
        // update particle positions
        (*disncol_)[lid_1] = pos_1_new[i];
        (*disncol_)[lid_2] = pos_2_new[i];
        // update particle velocities
        (*velncol_)[lid_1] = vel_1_new[i];
        (*velncol_)[lid_2] = vel_2_new[i];
      }
      std::cout << "New position of particle with GID " << particle_1->Id() << "  x: " << pos_1_new[0] << "  y: " << pos_1_new[1] << "  z: " << pos_1_new[2] << std::endl;
      std::cout << "New position of particle with GID " << particle_2->Id() << "  x: " << pos_2_new[0] << "  y: " << pos_2_new[1] << "  z: " << pos_2_new[2] << std::endl;
    }
  }
  break;
  case INPAR::PARTICLE::particle_wall:
  {
    // collision time
    double colltime = next_event->time;

    LINALG::Matrix<3,1> initposition;
    LINALG::Matrix<3,1> initvelocity;
    double radius;
    double particle_time;

    GetCollisionData(next_event->particle_1, initposition, initvelocity, radius, particle_time);

    // advance particle in time to collision time
    LINALG::Matrix<3,1> newpos;
    newpos.Update(1.0, initposition, colltime - particle_time, initvelocity);

    LINALG::Matrix<3,1> collnormal;
    collnormal.Update(1.0, newpos, -1.0, Teuchos::rcp_static_cast<WallEvent>(next_event)->wallcollpoint_pos);

    // safety check
    double normallength = collnormal.Norm2();
    if (std::abs(normallength - radius) > GEO::TOL7)
    {
      dserror("Particle and wall collision detected but distance does not match radius");
    }

    // normalize colnormal
    collnormal.Scale(1.0 / normallength);
    double veln_particle = collnormal.Dot(initvelocity);
    double veln_wall = Teuchos::rcp_static_cast<WallEvent>(next_event)->wallcollpoint_vel.Dot(collnormal);

    // walls have infinite mass
    LINALG::Matrix<3,1> newvel;
    newvel.Update(1.0, initvelocity, (1.0 + e_wall_) * (veln_wall - veln_particle), collnormal);

    // write particle data
    int lid = ddt_->Map().LID(next_event->particle_1->Id());
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
    std::cout << "New position of particle " << next_event->particle_1->Id() << "  x: " << newpos(0) << "  y: " << newpos(1) << "  z: " << newpos(2) << std::endl;

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
  // IMPORTANT: First Particle is the colliding particle which carries the current time in the ddt_-vector

  // create event in which the time still needs to be set correctly
  Teuchos::RCP<Event> newevent = Teuchos::rcp(new Event(INPAR::PARTICLE::particle_particle, -1000.0, particle_1, particle_2));

  if (particle_1->Id() == particle_2->Id())
    dserror("Particle %i cannot collide with itself!", particle_1->Id());

  LINALG::Matrix<3,1> pos_1, pos_2, vel_1, vel_2;
  double rad_1, rad_2, ddt_1, ddt_2;
  GetCollisionData(particle_1, particle_2, pos_1, pos_2, vel_1, vel_2, rad_1, rad_2, ddt_1, ddt_2);

  LINALG::Matrix<3,1> deltax, deltav;

  for (int i=0; i<3; ++i)
  {
    // first particle is the colliding particle thus it carries the actual time in the ddt_-vector
    // additionally the second particle needs to be updated to the time of particle_1
    deltax(i) = pos_1(i) - pos_2(i) - (ddt_1 - ddt_2) * vel_2(i);

    deltav(i) = vel_1(i) - vel_2(i);
  }

  double sigma = rad_1 + rad_2;

  // finding possible collision times
  double a = deltav.Dot(deltav);
  double b = 2.0 * deltav.Dot(deltax);
  double c = deltax.Dot(deltax) - sigma * sigma;

  double discriminant = b * b - 4.0 * a * c;

  if (discriminant >= 0.0)
  {
    double tc1 = (-b - sqrt(discriminant)) / (2 * a);
    double tc2 = (-b + sqrt(discriminant)) / (2 * a);

    // immediate collision of particles expected
    if (abs(tc1) <= GEO::TOL14 or abs(tc2) <= GEO::TOL14)
    {
      // compute collision normal to detect whether collision has already happened at the end of the last time step
      LINALG::Matrix<3,1> colnormal;
      colnormal.Update(1.0, pos_2, -1.0, pos_1);

      // velocities of particles in normal direction
      double vel_col_1 =   vel_1.Dot(colnormal);
      double vel_col_2 = - vel_2.Dot(colnormal);

      if ((vel_col_1 + vel_col_2) > GEO::TOL14)
      {
        newevent->time = 0.0;
      }
    }
    // tc1 is negative
    else if (tc1 < -GEO::TOL14 and tc2 > GEO::TOL14)
    {
      newevent->time = tc2;
    }
    // tc2 is negative
    else if (tc1 > GEO::TOL14 and tc2 < -GEO::TOL14)
    {
      newevent->time = tc1;
    }
    // both positive, smaller one is chosen (tc1 is almost identical to tc2)
    else if (tc1 > GEO::TOL14 and tc2 > GEO::TOL14)
    {
      newevent->time = tc1 < tc2 ? tc1 : tc2;
    }

    if(newevent->time > 1.1*remaining_dt)
      newevent->time = -1000.0;
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
  LINALG::Matrix<3,1> position;
  LINALG::Matrix<3,1> velocity;
  double radius;
  double particle_time;

  GetCollisionData(particle, position, velocity, radius, particle_time);

  // variables to be filled with collision data
  double timetocollision = 0.0;;
  LINALG::Matrix<3,1> wallcoll_pos(true);
  LINALG::Matrix<3,1> wallcoll_vel(true);

  // get wall discretization and displacement states
  Teuchos::RCP<DRT::Discretization> walldiscret = particle_algorithm_->WallDiscret();
  Teuchos::RCP<const Epetra_Vector> walldisn = walldiscret->GetState("walldisn");
  Teuchos::RCP<const Epetra_Vector> walldisnp = walldiscret->GetState("walldisnp");

  DRT::Node** nodes = wall->Nodes();
  int numnodes = wall->NumNode();

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
  ComputeCollisionOfParticleWithElement(wall->Shape(), xyze_n, xyze_np, position, velocity,
      radius, timetocollision, wallcoll_pos, wallcoll_vel, dt - particle_time, dt);

  // fill particle-wall event in which time to collision is inserted here so that current time needs to be added
  Teuchos::RCP<WallEvent> wallevent =
      Teuchos::rcp(new WallEvent(INPAR::PARTICLE::particle_wall, timetocollision, particle, wall, wallcoll_pos, wallcoll_vel));

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

  for (int i=0; i<discret_->NodeColMap()->NumMyElements(); ++i)
  {
    // particle for which contact will be detected
    DRT::Node *currparticle = discret_->lColNode(i);

    std::set<DRT::Node*> neighbouring_particles;
    std::set<DRT::Element*> neighbouring_walls;

    // gather all neighbouring particles and wall elements
    GetNeighbouringParticlesAndWalls(currparticle, neighbouring_particles, neighbouring_walls);

    // loop over all neighbouring particles and check if they collide with currparticle
    for (std::set<DRT::Node*>::iterator iter=neighbouring_particles.begin(); iter!=neighbouring_particles.end(); ++iter)
    {
      // in order to avoid double evaluation of the same contact, only search for contact when id_1 < id_2
      if (currparticle->Id() > (*iter)->Id())
      {
        continue; // with next particle in neighbourhood
      }

      Teuchos::RCP<PARTICLE::Event> newevent = ComputeCollisionWithParticle(currparticle, *iter, dt);

      // insert event into event queue if collision is valid
      if (newevent->time >= 0.0)
      {
        std::cout << "inserting inter-particle collision in the event queue" << std::endl;
        eventqueue.insert(newevent);
      }
    }

    // loop over all neighbouring wall elements and check if they collide with currparticle
    for (std::set<DRT::Element*>::iterator iter=neighbouring_walls.begin(); iter!=neighbouring_walls.end(); ++iter)
    {
      Teuchos::RCP<WallEvent> newevent = ComputeCollisionWithWall(currparticle, *iter, dt);

      // insert event into event queue if collision is valid
      if (newevent->time >= -GEO::TOL14)
      {
        std::cout << "inserting particle-wall collision in the event queue" << std::endl;
        // here we are at the beginning of the time step so that event.time does not need an update
        eventqueue.insert(newevent);
      }
    }
  }

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

    if (lastevent->coltype == INPAR::PARTICLE::particle_particle)
    {
      // do not add event of particles that has already been processed in lastevent
      if (newevent->time >= -GEO::TOL14 and newevent->particle_2->Id() != lastevent->particle_2->Id())
      {
        // only time to collision is returned --> time of the last event needs to be added
        newevent->time += lastevent->time;
        eventqueue.insert(newevent);
      }
    }
    else
    {
      if (newevent->time >= -GEO::TOL14)
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

    if (newevent->time >= -GEO::TOL14)
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
      if (newevent->time >= -GEO::TOL14 and newevent->particle_2->Id() != lastevent->particle_1->Id())
      {
        // only time to collision is returned --> time of the last event needs to be added
        newevent->time += lastevent->time;
        eventqueue.insert(newevent);
      }
    }

    // particle-wall collision
    for(std::set<DRT::Element*>::iterator iter=neighbouring_walls.begin(); iter!=neighbouring_walls.end(); ++iter)
    {
      Teuchos::RCP<WallEvent> newevent = ComputeCollisionWithWall(lastevent->particle_2, *iter, dt-lastevent->time);

      if (newevent->time >= -GEO::TOL14)
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
 | computes time to collision between a particle and and   ghamm 09/13  |
 | an element for hard sphere particles (templated on distype of ele)   |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType DISTYPE>
void PARTICLE::ParticleCollisionHandlerMD::ComputeCollisionOfParticleWithElementT(
    const Epetra_SerialDenseMatrix& xyze_n,
    const Epetra_SerialDenseMatrix& xyze_final,
    const LINALG::Matrix<3,1>& particle_pos,
    const LINALG::Matrix<3,1>& particle_vel,
    const double radius,
    double& timetocollision,
    LINALG::Matrix<3,1>& wallcollpoint_pos,
    LINALG::Matrix<3,1>& wallcollpoint_vel,
    const double remaining_dt,
    const double dt)
{
  const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

  // coll_solution contains: r, s and time to collision (r and s are element coords)
  static LINALG::Matrix<3,1> coll_solution;
  // initial guess for wall collision point
  GEO::startingValueCurrentToElementCoords<DISTYPE>(coll_solution);
  // starting time is zero
  coll_solution(2) = 0.0;

  // unit normal at collision point
  LINALG::Matrix<3,1> unitnormal;

  // velocity of wall element is constant over the time step
  Epetra_SerialDenseMatrix vele_colltime(xyze_n);
  vele_colltime.Scale(-1.0);
  vele_colltime += xyze_final;
  vele_colltime.Scale(1.0 / dt);

  Epetra_SerialDenseMatrix deltaxyze(xyze_n);
  deltaxyze.Scale(-1.0);
  deltaxyze += xyze_final;

  // the following equation is solved iteratively w.r.t. ele coords and time to collision (ttc):
  // particle_pos + ttc * particle_vel + UnitNormalAtWallCollPoint*radius = WallCollPoint at

  // iteration for contact search
  const int maxiter = 10;
  int iter = 0;
  while (iter < maxiter)
  {
    iter++;

    // compute wall position at collision time: xyze_time = xyze_current + (xyze_final - xyze_current)*coltime/timestep
    double colltime = dt-remaining_dt+coll_solution(2);
    Epetra_SerialDenseMatrix xyze_colltime(deltaxyze);
    xyze_colltime.Scale(colltime / dt);
    xyze_colltime += xyze_n;

    // get unit normal of wall element at collision point
    GEO::computeNormalToSurfaceElement(DISTYPE, xyze_colltime, coll_solution, unitnormal);

    // position of wall collision point
    GEO::elementToCurrentCoordinates(DISTYPE, xyze_colltime, coll_solution, wallcollpoint_pos);

    // check whether normal is pointing outward (particle is inside) and adapt it if necessary
    static LINALG::Matrix<3,1> testvector;
    testvector.Update(1.0, wallcollpoint_pos, -1.0, particle_pos);

    if (unitnormal.Dot(testvector) < 0.0)
      unitnormal.Scale(-1.0);

    // velocity of wall collision point
    GEO::elementToCurrentCoordinates(DISTYPE, vele_colltime, coll_solution, wallcollpoint_vel);

    // compute rhs
    static LINALG::Matrix<3,1> b;
    for (int i=0; i<3; ++i)
    {
      b(i) = -( particle_pos(i) + coll_solution(2) * particle_vel(i)  + unitnormal(i) * radius - wallcollpoint_pos(i) );
    }

    if (b.Norm2() < GEO::TOL14)
    {
      break;
    }

    // compute dxyzdrs at collision point and length of normal
    static LINALG::Matrix<2,numnode> deriv;
    DRT::UTILS::shape_function_2D_deriv1(deriv, coll_solution(0), coll_solution(1), DISTYPE);

    static LINALG::Matrix<3,2> dxyzdrs;
    dxyzdrs.Clear();
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 2; j++)
        for (int k = 0; k < numnode; k++)
          dxyzdrs(i, j) += xyze_colltime(i, k) * deriv(j, k);

    static LINALG::Matrix<3,1> nrm;
    nrm(0) = dxyzdrs(1, 0) * dxyzdrs(2, 1) - dxyzdrs(2, 0) * dxyzdrs(1, 1);
    nrm(1) = dxyzdrs(2, 0) * dxyzdrs(0, 1) - dxyzdrs(0, 0) * dxyzdrs(2, 1);
    nrm(2) = dxyzdrs(0, 0) * dxyzdrs(1, 1) - dxyzdrs(1, 0) * dxyzdrs(0, 1);
    double normallength = nrm.Norm2();

    // compute d2xyzdrs at collision point
    static LINALG::Matrix<3,numnode> deriv2;
    DRT::UTILS::shape_function_2D_deriv2(deriv2, coll_solution(0), coll_solution(1), DISTYPE);

    static LINALG::Matrix<3,3> d2xyzdrs;
    d2xyzdrs.Clear();
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        for (int k = 0; k < numnode; k++)
          d2xyzdrs(i, j) += xyze_colltime(i, k) * deriv2(j, k);

    // compute dnormal/drs at collision point
    static LINALG::Matrix<3,2> dnormal;

    // second derivative in r
    dnormal(0, 0) = d2xyzdrs(1, 0) * dxyzdrs(2, 1) + d2xyzdrs(2, 2) * dxyzdrs(1, 0)
                  - d2xyzdrs(2, 0) * dxyzdrs(1, 1) - d2xyzdrs(1, 2) * dxyzdrs(2, 0);
    dnormal(1, 0) = d2xyzdrs(2, 0) * dxyzdrs(2, 0) + d2xyzdrs(0, 2) * dxyzdrs(0, 1)
                  - d2xyzdrs(0, 0) * dxyzdrs(2, 1) - d2xyzdrs(2, 2) * dxyzdrs(0, 0);
    dnormal(2, 0) = d2xyzdrs(0, 0) * dxyzdrs(0, 0) + d2xyzdrs(1, 2) * dxyzdrs(1, 1)
                  - d2xyzdrs(1, 0) * dxyzdrs(0, 1) - d2xyzdrs(0, 2) * dxyzdrs(1, 0);

    // second derivative in s
    dnormal(0, 1) = d2xyzdrs(1, 2) * dxyzdrs(2, 1) + d2xyzdrs(2, 1) * dxyzdrs(1, 0)
                  - d2xyzdrs(2, 2) * dxyzdrs(1, 1) - d2xyzdrs(1, 1) * dxyzdrs(2, 0);
    dnormal(1, 1) = d2xyzdrs(2, 2) * dxyzdrs(2, 0) + d2xyzdrs(0, 1) * dxyzdrs(0, 1)
                  - d2xyzdrs(0, 2) * dxyzdrs(2, 1) - d2xyzdrs(2, 1) * dxyzdrs(0, 0);
    dnormal(2, 1) = d2xyzdrs(0, 2) * dxyzdrs(0, 0) + d2xyzdrs(1, 1) * dxyzdrs(1, 1)
                  - d2xyzdrs(1, 2) * dxyzdrs(0, 1) - d2xyzdrs(0, 1) * dxyzdrs(1, 0);

    // compute dF/dx
    static LINALG::Matrix<3,3> A;

    // first/second column: derivative w. r.t. r/s
    for (int j = 0; j < 2; j++)
    {
      for (int i = 0; i < 3; i++)
      {
        A(i, j) = dnormal(i, j) * (1.0 / normallength) * radius - dxyzdrs(i, j);
      }
    }
    // third column: derivative with respect to time
    for (int i = 0; i < 3; i++)
    {
      A(i, 2) = particle_vel(i) - wallcollpoint_vel(i);
    }

    // solve linear problem A dx = b
    static LINALG::Matrix<3,1> dx;
    double det = LINALG::gaussElimination<true,3>(A, b, dx);

    if (fabs(det) < GEO::TOL14)
    {
//      std::cout << "Determinant near zero meaning, the particle path is parallel to the wall" << std::endl;
      if(unitnormal.Dot(particle_vel) < GEO::TOL10 and iter>1)
      {
//        std::cout << "particle path is parallel to wall --> left iteration" << std::endl;
        coll_solution(2) = -1000.0;
        break;
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
      coll_solution(2) = -1000.0;
      break;
    }
  }

//  std::cout << "needed " << iter << " iterations" << " with solution: " <<  coll_solution(0) << " " << coll_solution(1) << " " <<  coll_solution(2) << std::endl;

  // check if collision is valid
  timetocollision = -1000.0;

  // check whether collision position lies on element
  if (GEO::checkPositionWithinElementParameterSpace(coll_solution, DISTYPE) == false)
  {
    return;
  }

  // check whether collision time is reasonable
  if (coll_solution(2) >= -GEO::TOL14 and coll_solution(2) < 1.1 * remaining_dt)
  {
    double scalar_partvel = unitnormal.Dot(particle_vel);
    double scalar_wallvel = unitnormal.Dot(wallcollpoint_vel);
    // decide if collision is still going to happen (--> valid) or has already happened in the last time step (--> invalid)
    if ((scalar_partvel - scalar_wallvel) > GEO::TOL14)
    {
      timetocollision = coll_solution(2);
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | computes time to collision between a particle and       ghamm 09/13  |
 | an element for hard sphere particles                                 |
 *----------------------------------------------------------------------*/
void PARTICLE::ParticleCollisionHandlerMD::ComputeCollisionOfParticleWithElement(
    const DRT::Element::DiscretizationType distype,
    const Epetra_SerialDenseMatrix& xyze_current,
    const Epetra_SerialDenseMatrix& xyze_final,
    const LINALG::Matrix<3,1>& particle_pos,
    const LINALG::Matrix<3,1>& particle_vel,
    const double radius,
    double& timetocollision,
    LINALG::Matrix<3,1>& wall_pos,
    LINALG::Matrix<3,1>& wall_vel,
    const double remaining_dt,
    const double dt)
{
  switch (distype)
  {
  case DRT::Element::quad4:
    ComputeCollisionOfParticleWithElementT<DRT::Element::quad4>(
        xyze_current, xyze_final, particle_pos, particle_vel, radius, timetocollision, wall_pos, wall_vel, remaining_dt, dt);
    break;
  case DRT::Element::quad8:
    ComputeCollisionOfParticleWithElementT<DRT::Element::quad8>(
        xyze_current, xyze_final, particle_pos, particle_vel, radius, timetocollision, wall_pos, wall_vel, remaining_dt, dt);
    break;
  case DRT::Element::quad9:
    ComputeCollisionOfParticleWithElementT<DRT::Element::quad9>(
        xyze_current, xyze_final, particle_pos, particle_vel, radius, timetocollision, wall_pos, wall_vel, remaining_dt, dt);
    break;
  case DRT::Element::tri3:
    ComputeCollisionOfParticleWithElementT<DRT::Element::tri3>(
        xyze_current, xyze_final, particle_pos, particle_vel, radius, timetocollision, wall_pos, wall_vel, remaining_dt, dt);
    break;
  case DRT::Element::tri6:
    ComputeCollisionOfParticleWithElementT<DRT::Element::tri6>(
        xyze_current, xyze_final, particle_pos, particle_vel, radius, timetocollision, wall_pos, wall_vel, remaining_dt, dt);
    break;
  default:
    dserror("please add your surface element type here");
    break;
  }

  return;
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
  std::vector<double> position(3);
  std::vector<double> velocity(3);
  DRT::UTILS::ExtractMyValues(*disncol_, position, lm);
  DRT::UTILS::ExtractMyValues(*velncol_, velocity, lm);
  int lid = radiusncol_->Map().LID(particle->Id());
  rad = (*radiusncol_)[lid];
  ddt = (*ddt_)[lid];

  for(int d=0; d<3; ++d)
  {
   pos(d) = position[d];
   vel(d) = velocity[d];
  }

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
    if (abs(event1->time - event2->time) < GEO::TOL14)
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
    if (abs(event1->time - event2->time) < GEO::TOL14)
    {
      int gid1 = event1->particle_1->Id();
      int gid2 = event1->particle_2->Id();
      int gid3 = event2->particle_1->Id();

      if (gid1 == gid3 || gid2 == gid3)
      {
        std::cout << ("ERROR: TWO PARTICLES AND WALL COLLIDING AT THE SAME TIME!") << std::endl;
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
    if (abs(event1->time - event2->time) < GEO::TOL14)
    {
      int gid1 = event1->particle_1->Id();
      int gid3 = event2->particle_1->Id();
      int gid4 = event2->particle_2->Id();

      if (gid1 == gid3 || gid1 == gid4)
      {
        std::cout << ("ERROR: TWO PARTICLES AND WALL COLLIDING AT THE SAME TIME!") << std::endl;
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
    if (abs(event1->time - event2->time) < GEO::TOL14)
    {
      int gid1 = event1->particle_1->Id();
      int gid3 = event2->particle_1->Id();

      if (gid1 == gid3)
      {
        std::cout << ("ERROR: PARTICLE AND TWO WALLS COLLIDING AT THE SAME TIME!") << std::endl;
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

