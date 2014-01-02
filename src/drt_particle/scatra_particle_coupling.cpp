/*----------------------------------------------------------------------*/
/*!
\file scatra_particle_coupling.cpp

\brief Algorithm to track particles for level-set problems

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 11/12 |
 *----------------------------------------------------------------------*/

#include "particle_algorithm.H"
#include "scatra_particle_coupling.H"
#include "particle_timint_rk.H"
#include "../drt_adapter/adapter_particle.H"
#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_meshfree_discret/drt_meshfree_multibin.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_parallel.H"

#include "../drt_mortar/mortar_utils.H"
#include "../drt_mortar/mortar_coupling3d_classes.H"

#include "../drt_geometry/searchtree.H"
#include "../drt_geometry/searchtree_geometry_service.H"
#include "../drt_geometry/intersection_math.H"
#include "../drt_geometry/element_coordtrafo.H"
#include "../drt_geometry/position_array.H"
#include "../drt_geometry/element_volume.H"
#include "../drt_cut/cut_position.H"
#include "../linalg/linalg_utils.H"

#include "../drt_io/io_pstream.H"
#include "../headers/definitions.h"

#include <Epetra_FEVector.h>
#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 | algorithm constructor                                rasthofer 09/13 |
 *----------------------------------------------------------------------*/
PARTICLE::ScatraParticleCoupling::ScatraParticleCoupling(
  Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,
  Teuchos::RCP<Teuchos::ParameterList> params
  ) : PARTICLE::Algorithm(scatra->Discretization()->Comm(),*params),
  scatra_(scatra),
  scatradis_(scatra->Discretization()),
  params_(params),
  particle_dim_(DRT::INPUT::IntegralValue<INPAR::PARTICLE::ParticleDim>(params->sublist("PARTICLE"),"DIMENSION"))
{

  Init(false);

  return;
}


/*----------------------------------------------------------------------*
 | initialization of the system                             ghamm 11/12 |
 *----------------------------------------------------------------------*/
void PARTICLE::ScatraParticleCoupling::Init(bool restarted)
{
  // -------------------------------------------------------------------
  //               setup particle discretization
  // -------------------------------------------------------------------

  // FillComplete() necessary for DRT::Geometry .... could be removed perhaps
  particledis_->FillComplete(false,false,false);
  // extract noderowmap because it will be called Reset() after adding elements
  Teuchos::RCP<Epetra_Map> particlerowmap = Teuchos::rcp(new Epetra_Map(*particledis_->NodeRowMap()));
  CreateBins(scatradis_);

  // gather all scatra coleles in each bin for proper extended ghosting
  std::map<int, std::set<int> > scatraelesinbins;
  Teuchos::RCP<Epetra_Map> binrowmap = DistributeBinsToProcsBasedOnUnderlyingDiscret(scatradis_, scatraelesinbins);

  //--------------------------------------------------------------------
  // -> 1) create a set of homeless particles that are not in a bin on this proc
  std::set<Teuchos::RCP<DRT::Node>, BINSTRATEGY::Less> homelessparticles;

  for (int lid = 0; lid < particlerowmap->NumMyElements(); ++lid)
  {
    DRT::Node* node = particledis_->gNode(particlerowmap->GID(lid));
    const double* currpos = node->X();
    PlaceNodeCorrectly(Teuchos::rcp(node,false), currpos, homelessparticles);
  }

  // start round robin loop to fill particles into their correct bins
  FillParticlesIntoBins(homelessparticles);

  // ghost bins, particles and fluid elements according to the bins
  SetupGhosting(binrowmap, scatraelesinbins);

  // some output
  IO::cout << "after ghosting" << IO::endl;
  DRT::UTILS::PrintParallelDistribution(*particledis_);
  DRT::UTILS::PrintParallelDistribution(*scatradis_);

  // -------------------------------------------------------------------
  //               setup time integration
  // -------------------------------------------------------------------

  // the following has only to be done once --> skip in case of restart
  // remark:
  // read restart deletes all initial particles, here the ones of the inital seeding
  // this is ok
  // then Init() is called once more
  // therefore, we have to ensure that the time integrator is only setup once
  if(not restarted)
  {
    // create time integrator
    Teuchos::RCP<ADAPTER::ParticleBaseAlgorithm> particles =
        Teuchos::rcp(new ADAPTER::ParticleBaseAlgorithm((*params_), particledis_));
    particles_ = particles->ParticleField();
    // set particle algorithm into time integration
    particles_->SetParticleAlgorithm(Teuchos::rcp(this,false));
  }

  return;
}


/*----------------------------------------------------------------------*
 | initialize particle field with particles             rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::ScatraParticleCoupling::InitialSeeding()
{

  // -------------------------------------------------------------------
  //               setup initial seeding
  // -------------------------------------------------------------------

  // get number of particles of each sign per bin around interface
  const int num_particles_per_bin = params_->sublist("PARTICLE").get<int>("NUMPARTICLE");

  // get characteristic element length of scatra discretization
  // assumed equal to bin edge length
  // TODO: remove
//  std::cout << bin_size_[0] << "  " << bin_size_[1] << "  " << bin_size_[2] << std::endl;
//  if (std::abs(bin_size_[0]-bin_size_[1]) > 10e-9 or std::abs(bin_size_[0]-bin_size_[2]) > 10e-9)
//    dserror("Cubic bins expected");
  const double binlength = cutoff_radius_; // should be equal to bin_size_[0] for the present settings

  // define band of elements around interface to be filled with particles, i.e., band of size bandwith * h
  // on both sides of the interface
  const double bandwith = params_->sublist("PARTICLE").get<double>("PARTICLEBANDWIDTH");

  // -------------------------------------------------------------------
  //               find bins to be filled
  // -------------------------------------------------------------------

  // vector containing id of bins (equal to elements) which will be filled with particles
  std::vector<int> particlebins;

  // get phinp
  const Teuchos::RCP<const Epetra_Vector> phinp = scatra_->Phinp();

  // loop all row bins on this proc
  for (int ibin = 0; ibin < particledis_->NumMyRowElements(); ibin++)
  {
    // get pointer to current bin
    DRT::Element* actele = particledis_->lRowElement(ibin);
    DRT::MESHFREE::MeshfreeMultiBin* actbin = static_cast<DRT::MESHFREE::MeshfreeMultiBin*>(actele);

    // get pointer to associated scatra elements
    DRT::Element** scatraelesinbin = actbin->AssociatedFluidEles();

    // bool to go to next bin if this bin is to far away from the interface
    // or as soon as this bin is identified to be sufficiently close to set particles
    bool next_bin = false;

    // loop all elements in bin
    for(int iele = 0; iele < actbin->NumAssociatedFluidEle(); iele++)
    {
      DRT::Element* scatraele = scatraelesinbin[iele];
      if (scatraele->Shape() != DRT::Element::hex8)
       dserror("Other element than hex8 not yet supported!");

      // get nodal phi values
      std::vector<int> lm;
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      lm.clear();
      lmowner.clear();
      lmstride.clear();
      scatraele->LocationVector(*scatradis_,lm,lmowner,lmstride);
      std::vector<double> myphinp(lm.size());
      DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);

      // check, if bin is close to the interface
      if (myphinp[0] > (2.0 * bandwith * binlength))
      {
        // bin is located away from the interface
        // abort loop elements and go to the next bin
        break; // loop all elements
      }

      // get element center
      const DRT::Element::DiscretizationType distype = DRT::Element::hex8;
      const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      // get center coordinates of current element
      LINALG::Matrix<3,1> elecenter(true);
      LINALG::Matrix<3,numnode> xyzele(true);
      GEO::fillInitialPositionArray<distype>(scatraele, xyzele);

      for (int inode = 0; inode < numnode; inode++)
        for (int idim = 0; idim < 3; idim++)
          elecenter(idim,0) += xyzele(idim,inode);

      for (int idim = 0; idim < 3; idim++)
        elecenter(idim,0) /= numnode;

      // assuming bin length to be chosen such that only one element has its center inside the bin,
      // we only check the distance of this element to the interface
      // get center of current bin
      LINALG::Matrix<3,1> center(true);
      center = GetBinCentroid(ibin);
      if ((elecenter(0,0) > (center(0,0) - (bin_size_[0]/2.0))) and (elecenter(0,0) < (center(0,0) + (bin_size_[0]/2.0))))
      {
        if ((elecenter(1,0) > (center(1,0) - (bin_size_[1]/2.0))) and (elecenter(1,0) < (center(1,0) + (bin_size_[1]/2.0))))
        {
          if ((elecenter(2,0) > (center(2,0) - (bin_size_[1]/2.0))) and (elecenter(2,0) < (center(2,0) + (bin_size_[2]/2.0))))
          {
            // loop all nodes of this element
            for (int inode = 0; inode < numnode; inode++)
              if (std::abs(myphinp[inode]) < (bandwith * binlength))
              {
                // bin is close to the interface and particles should be set
                // add bin to list
                particlebins.push_back(ibin);
                // element with center in bin is found and we may go to the next bin
                next_bin = true;
                break;
              }
          }
        }
      }

      // note: we only abort the element loop if the element with center in the bin
      //       is sufficiently close to the interface; if this is not the case,
      //       we loop the remaining elements although we expect only one element
      //       with center in the bin; with respect to inhomogeneous meshes, for which
      //       more elements with center in the bin cannot be excluded, this choice represents
      //       a saver strategy
      if (next_bin == true)
        break; // loop all elements

    } // end loop all elements

  }  // end loop all bins

  // Output for testing
  std::cout << "bins with particles " << particlebins.size() << std::endl;
//  for (std::size_t i=0; i<particlebins.size(); i++)
//    std::cout << particlebins[i] << std::endl;
//  dserror("ENDE");

  // -------------------------------------------------------------------
  //               seed particles in selected bins
  // -------------------------------------------------------------------

  // initialize particle id with largest particle id in use + 1 (on each proc)
  int maxparticleid = particledis_->NodeRowMap()->MaxAllGID();
  int currentparticleid = maxparticleid;

  // list of all particles and their position
  std::map<int,std::vector<double> > particle_positions;
  // list of all particles and their sign
  std::map<int,double > particle_sign;

  // loop all bins which should be filled
  for (std::size_t ibin = 0; ibin < particlebins.size(); ibin++)
  {
    // set particle id
    currentparticleid += 1;

    // get center of current bin
    LINALG::Matrix<3,1> center(true);
    center = GetBinCentroid(ibin);
    // get max and min values of bin corners
    std::vector<LINALG::Matrix<3,1> > bincorners(2);
    for (int rr = 0; rr<2; rr++)
    {
      for (int idim=0; idim<3; idim++)
        (bincorners[rr])(idim,0) = center(idim,0) + pow(-1.0,(double)rr+1.0) * 0.5 * bin_size_[idim];
    }
//    double xmin = center(0,0)-(0.5*bin_size_[0]);
//    double xmax = center(0,0)+(0.5*bin_size_[0]);
//
//    double ymin = center(1,0)-(0.5*bin_size_[1]);
//    double ymax = center(1,0)+(0.5*bin_size_[1]);
//
//    double zmin = center(2,0)-(0.5*bin_size_[2]);
//    double zmax = center(2,0)+(0.5*bin_size_[2]);

    // loop number of particles to be seeded
    for (int ipart = 0; ipart < (2 * num_particles_per_bin); ipart++)
    {
      // get position randomly in bin
      std::vector<double> position(3);

      DRT::UTILS::Random* random = DRT::Problem::Instance()->Random();

      for (int idim=0; idim<3; idim++)
      {
        // set range (default: [-1;1])
        random->SetRandRange((bincorners[0])(idim,0),(bincorners[1])(idim,0));
        // get position
        position[idim] = random->Uni();
      }

      // set particle to mid plane in case of quasi-2D simulations
      switch (particle_dim_)
      {
        case INPAR::PARTICLE::particle_2Dx:
        {
          position[0] = 0.0;
          break;
        }
        case INPAR::PARTICLE::particle_2Dy:
        {
          position[1] = 0.0;
          break;
        }
        case INPAR::PARTICLE::particle_2Dz:
        {
          position[2] = 0.0;
          break;
        }
        default:
          break;
      }

      // place particle with id and position
      std::set<Teuchos::RCP<DRT::Node>,BINSTRATEGY::Less> homelessparticles;
      Teuchos::RCP<DRT::Node> newparticle = Teuchos::rcp(new DRT::Node(currentparticleid, &position[0], myrank_));
      PlaceNodeCorrectly(newparticle, &position[0], homelessparticles);
      if (homelessparticles.size() != 0)
        dserror("New particle could not be inserted on this proc!.");

      // add particle to list
      particle_positions.insert(std::pair<int,std::vector<double> >(currentparticleid,position));
      // set sign (+ for the first half, - for the second one)
      if (ipart < num_particles_per_bin)
        particle_sign.insert(std::pair<int,double>(currentparticleid,1.0));
      else
        particle_sign.insert(std::pair<int,double>(currentparticleid,-1.0));
    }
  }

  // rebuild connectivity and assign degrees of freedom (note: IndependentDofSet)
  particledis_->FillComplete(true, false, true);

  // update of state vectors to the new maps
  Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->UpdateStatesAfterParticleTransfer();

  // -------------------------------------------------------------------
  //               initialize state vectors
  // -------------------------------------------------------------------

  // get position vector at time_(n+1)
  Teuchos::RCP<Epetra_Vector> disnp = particles_->WriteAccessDispnp();
  // get sign vector
  Teuchos::RCP<Epetra_Vector> sign = Teuchos::rcp_dynamic_cast<PARTICLE::TimIntRK>(particles_)->WriteAccessSign();

  // get maps
  const Epetra_Map* dofrowmap = particledis_->DofRowMap();
  const Epetra_Map* noderowmap = particledis_->NodeRowMap();

  for (int inode=0; inode<particledis_->NumMyRowNodes(); inode++)
  {
    // get global id of current particle
    const int gid = noderowmap->GID(inode);

    DRT::Node* currparticle = particledis_->gNode(gid);
    // get the first dof gid of a particle and convert it into a LID
    int doflid = dofrowmap->LID(particledis_->Dof(currparticle, 0));

    // insert values into state vectors
    (*sign)[inode] = particle_sign[gid];
    for (int idim=0; idim<3; ++idim)
      (*disnp)[doflid+idim] = (particle_positions[gid])[idim];
  }

  // get position vector at time_(n)
  Teuchos::RCP<Epetra_Vector> disn = particles_->WriteAccessDispn();
  disn->Update(1.0,*disnp,0.0);
  // reset to zero just to be sure
  Teuchos::RCP<Epetra_Vector> radius = particles_->WriteAccessRadius();
  radius->Scale(0.0);

  // -------------------------------------------------------------------
  //               attraction step
  // -------------------------------------------------------------------
  // move particle to the correct side of the interface
  
  // prepare required vectors
  // noderowmap-based vectors
  // predictor
  Teuchos::RCP<Epetra_Vector> lambda = Teuchos::rcp(new Epetra_Vector(*noderowmap,true));
  // intialize with 1.0
  lambda->PutScalar(1.0);
  // target phi_goal
  Teuchos::RCP<Epetra_Vector> phi_target = Teuchos::rcp(new Epetra_Vector(*noderowmap,true));
  // current phi of particle
  Teuchos::RCP<Epetra_Vector> phi_particle = Teuchos::rcp(new Epetra_Vector(*noderowmap,true));
  // dofrowmap-based vectors
  // normal vector of level-set field at particle position
  Teuchos::RCP<Epetra_Vector> norm_gradphi_particle = Teuchos::rcp(new Epetra_Vector(*dofrowmap,true));
    // increment position vector for attraction loop
  Teuchos::RCP<Epetra_Vector> inc_dis = Teuchos::rcp(new Epetra_Vector(*dofrowmap,true));

  // get minimal radius of particles
  const double r_min = params_->sublist("PARTICLE").get<double>("MIN_RADIUS");
  if (r_min < 0.0) dserror("Set minimal particle radius in input file!");

  // loop all particles and intialize phi_target, phi_particle and gradphi_particle
  for (int inode=0; inode<particledis_->NumMyRowNodes(); inode++)
  {
    // get sign of particle
    const double signP = (*sign)[inode];

    // get target phi for this particle
    DRT::UTILS::Random* random = DRT::Problem::Instance()->Random();
    random->SetRandRange(r_min,(bandwith * binlength));
    const double local_phi_target = signP * random->Uni();
    (*phi_target)[inode] = local_phi_target;

    // get phi and gradient of particle at current particle position
    double current_phi_particle = 0.0;
//    current_phi_particle = GetPhiParticle(inode);
    LINALG::Matrix<3,1> normal_particle(true);

    // element coordinates of particle position in scatra element
    LINALG::Matrix<3,1> elecoord(true);
    DRT::Element* scatraele = GetEleCoordinates(inode,elecoord);

    // compute level-set value
    GetLSValParticle(current_phi_particle,normal_particle,scatraele,elecoord);
    // GetLSValParticle return gardient of phi, which we have to normalize here
    const double norm = normal_particle.Norm2();
    normal_particle.Scale(1.0/norm);

    // store values
    (*phi_particle)[inode] = current_phi_particle;

    // get global id of current particle
    const int gid = noderowmap->GID(inode);

    DRT::Node* currparticle = particledis_->gNode(gid);
    // get the first dof gid of a particle and convert it into a LID
    int doflid = dofrowmap->LID(particledis_->Dof(currparticle, 0));
    for (int idim=0; idim<3; ++idim)
      (*norm_gradphi_particle)[doflid+idim] = normal_particle(idim,0);
    
    // initialize increment position vector, i.e., multiply by lambda*(phi_target-current_phi_particle)
    // note lambda is set to 1.0 at the beginning and skiped here
    const double fac = local_phi_target - current_phi_particle;
    normal_particle.Scale(fac);
    for (int idim=0; idim<3; ++idim)
      (*inc_dis)[doflid+idim] = normal_particle(idim,0);
  }

  // number of placed particles should finally be equal to the number of seeded particles
  // important: this number must be global and takes into account the placed particles on
  // all procs
  int num_placed_particles = 0;
  // total (global) number of seeded particles
  const int global_num_particles = particledis_->NumGlobalNodes();
  // counter for number of attraction steps
  int step_counter = 0;
  // define vector for update of map-based vectors
  Teuchos::RCP<Epetra_Vector> old;

  // attraction loop
  while ((num_placed_particles < global_num_particles) and step_counter < 15)
  {
    // increment counter
    step_counter += 1;

    // compute new position
    disnp = particles_->WriteAccessDispnp();
    disn = particles_->WriteAccessDispn();
    sign = Teuchos::rcp_dynamic_cast<PARTICLE::TimIntRK>(particles_)->WriteAccessSign();
    // x_P^new = x_P + lambda*(phi_target-current_phi_particle) * n_p
    disnp->Update(1.0,*inc_dis,1.0);

    // transfer particles to bins
    TransferParticles(false);
    // update of state vectors to the new maps
    Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->UpdateStatesAfterParticleTransfer();
    // likewise update present vectors according to the new distribution of particles
    old = lambda;
    lambda = LINALG::CreateVector(*particledis_->NodeRowMap(),true);
    LINALG::Export(*old, *lambda);

    old = phi_target;
    phi_target = LINALG::CreateVector(*particledis_->NodeRowMap(),true);
    LINALG::Export(*old, *phi_target);

    old = phi_particle;
    phi_particle = LINALG::CreateVector(*particledis_->NodeRowMap(),true);
    LINALG::Export(*old, *phi_particle);

    old = norm_gradphi_particle;
    norm_gradphi_particle = LINALG::CreateVector(*particledis_->DofRowMap(),true);
    LINALG::Export(*old, *norm_gradphi_particle);

    old = inc_dis;
    inc_dis = LINALG::CreateVector(*particledis_->DofRowMap(),true);
    LINALG::Export(*old, *inc_dis);

    // check new position
    // counter for local placed particles
    int local_num_placed_particles = 0;
    
    // loop all particles
    for (int inode=0; inode<particledis_->NumMyRowNodes(); inode++)
    {
      // get phi and gradient of particle at current particle position
      double current_phi_particle = 0.0;
      LINALG::Matrix<3,1> normal_particle(true);

      // element coordinates of particle position in scatra element
      LINALG::Matrix<3,1> elecoord(true);
      DRT::Element* scatraele = GetEleCoordinates(inode,elecoord);

      // get dof lid
      const int gid = noderowmap->GID(inode);
      DRT::Node* currparticle = particledis_->gNode(gid);
      int doflid = dofrowmap->LID(particledis_->Dof(currparticle, 0));

      if (scatraele != NULL) // particle has not left domain
      {
        // compute level-set value
        GetLSValParticle(current_phi_particle,normal_particle,scatraele,elecoord);
        // GetLSValParticle return gradient of phi, which we have to normalize here
        const double norm = normal_particle.Norm2();
        normal_particle.Scale(1.0/norm);
        
        const double signP = (*sign)[inode];
        if ((current_phi_particle >= (signP * r_min))
            and (current_phi_particle <= (signP * bandwith * binlength))
            and ((*lambda)[inode] == 1.0)) // particle in correct band
        {
          // this means particle is placed correctly
          local_num_placed_particles += 1;
          // store actual phi value
          (*phi_particle)[inode] = current_phi_particle;
          // set lambda to zero to avoid further movement of this particle
          (*lambda)[inode] = 0.0;
        }
        else
        {
          if ((*lambda)[inode] == 1.0) // we could not move the particle to the correct band
          {
            // we go back to our initial position
            // we do not overwrite phi_particle and norm_gradphi_particle
            // we set disnp back to disn
            for (int idim=0; idim<3; ++idim)
              (*disnp)[doflid+idim] = (*disn)[doflid+idim];
            // and only go half of the distance in the next step
            (*lambda)[inode] = 0.5;
          }
          else // we accept the new position and start the attraction step from this point with lambda again set to 1.0
          {
            // store new phi value
            (*phi_particle)[inode] = current_phi_particle;
            // set lamda to 1.0
            (*lambda)[inode] = 1.0;
            // store new gradient
            for (int idim=0; idim<3; ++idim)
              (*norm_gradphi_particle)[doflid+idim] = normal_particle(idim,0);
          }
        }
      }
      else
      {
        // if we are outside of the domain, we go back to our initial position
        // we do not overwrite phi_particle and norm_gradphi_particle
        // we set disnp back to disn
        for (int idim=0; idim<3; ++idim)
          (*disnp)[doflid+idim] = (*disn)[doflid+idim];
        // and start another attraction step with lambda halved
        double mylambda = (*lambda)[inode];
        (*lambda)[inode] = 0.5 * mylambda;
      }

      // set new increment vector
      const double fac = (*lambda)[inode] * ((*phi_target)[inode] - (*phi_particle)[inode]);
      for (int idim=0; idim<3; ++idim)
        (*inc_dis)[doflid+idim] = fac * (*norm_gradphi_particle)[doflid+idim];

    } // end loop all particles

    // collect correctly placed particles of this step form all procs
    int global_num_placed_particles = 0;
    scatradis_->Comm().SumAll(&local_num_placed_particles,&global_num_placed_particles,1);
    // update number of placed particles
    num_placed_particles += global_num_placed_particles;

    // possibly position have been reverted and particles changed the bin or even the processor
    // therefore, we have to recall the following steps
    TransferParticles(false);
    // update of state vectors to the new maps
    Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->UpdateStatesAfterParticleTransfer();
    // likewise update present vectors according to the new distribution of particles
    Teuchos::RCP<Epetra_Vector> old;
    old = lambda;
    lambda = LINALG::CreateVector(*particledis_->NodeRowMap(),true);
    LINALG::Export(*old, *lambda);

    old = phi_target;
    phi_target = LINALG::CreateVector(*particledis_->NodeRowMap(),true);
    LINALG::Export(*old, *phi_target);

    old = phi_particle;
    phi_particle = LINALG::CreateVector(*particledis_->NodeRowMap(),true);
    LINALG::Export(*old, *phi_particle);

    old = norm_gradphi_particle;
    norm_gradphi_particle = LINALG::CreateVector(*particledis_->DofRowMap(),true);
    LINALG::Export(*old, *norm_gradphi_particle);

    old = inc_dis;
    inc_dis = LINALG::CreateVector(*particledis_->DofRowMap(),true);
    LINALG::Export(*old, *inc_dis);

    // update state vectors
    disnp = particles_->WriteAccessDispnp();
    disn = particles_->WriteAccessDispn();
    disn->Update(1.0,*disnp,0.0);
  }
  
  // delete particles which could not be attracted correctly -> all particles with lambda != 0.0
  if (num_placed_particles != global_num_particles) // we have particles with lambda != 0.0
  {
    // loop all particles
    for (int inode=0; inode<particledis_->NumMyRowNodes(); inode++)
    {
      // particles with lambda != 0.0 will be deleted
      if ((*lambda)[inode] > (0.0+1.0e-8))
      {
        // get particle
        DRT::Node *currparticle = particledis_->lRowNode(inode);
        // get gid
        const double gid = currparticle->Id();
        if(currparticle->NumElement() != 1)
          dserror("ERROR: A particle is assigned to more than one bin!");
        // get corresponding bin
        DRT::Element** currele = currparticle->Elements();
        DRT::Element* currbin = currele[0];

        // remove particle from bin
        static_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currbin)->DeleteNode(gid);

        // remove particle from discretization
        particledis_->DeleteNode(gid);
      }
    }

    // finish discretization
    particledis_->FillComplete(true,false,false);
    // adapt state vectors to new maps
    Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->UpdateStatesAfterParticleTransfer();
    // adapt also phi_particle which will be used below to set the radius
    old = phi_particle;
    phi_particle = LINALG::CreateVector(*particledis_->NodeRowMap(),true);
    LINALG::Export(*old, *phi_particle);
  }

  // -------------------------------------------------------------------
  //               set radius
  // -------------------------------------------------------------------
  
  // get radius
  radius = particles_->WriteAccessRadius();

  // get maximal radius of particles
  const double r_max = params_->sublist("PARTICLE").get<double>("MAX_RADIUS");
  if (r_max < 0.0) dserror("Set maximal particle radius in input file!");

  for (int inode=0; inode<particledis_->NumMyRowNodes(); inode++)
  {
    // get phi
    const double myphi = (*phi_particle)[inode];
    if (std::abs(myphi) > r_max)
      (*radius)[inode] = r_max;
    else if (std::abs(myphi) < r_min)
      (*radius)[inode] = r_min;
    else
      (*radius)[inode] = std::abs(myphi);
  }

  return;
}


/*----------------------------------------------------------------------*
 | prepare time step                                   rasthofer 09/13  |
 *----------------------------------------------------------------------*/
void PARTICLE::ScatraParticleCoupling::PrepareTimeStep()
{
  // update time and step
  Algorithm::PrepareTimeStep();
  return;
}


/*----------------------------------------------------------------------*
 | solve the current particle time step                 rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::ScatraParticleCoupling::Integrate()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::Algorithm::Integrate");
  particles_->IntegrateStep();

  return;
}


void PARTICLE::ScatraParticleCoupling::SetVelocity()
{
  // TODO: insert velocity in particles: veln_, vel_(0) whatever you need :-)
#if 0
  //Define pointers to vectors veln and disn (which are defined in timint.cpp)
  Teuchos::RCP<Epetra_Vector> velnp = particles_->ExtractVelnp();
  Teuchos::RCP<Epetra_Vector> disnp = particles_->ExtractDispnp();

  const int nsd = 3;

  LINALG::Matrix<nsd,1> particleposition;
  LINALG::Matrix<nsd,1> elecoord(true);
  DRT::Element* targetscatraele = NULL;

  int NumParticles = particledis_->NumMyRowNodes();
  std::cout << "Number of Nodes (=Particles) is " << NumParticles << endl;

  int maxparticleid = particledis_->NodeRowMap()->MaxAllGID();
  std::cout << "maxparticleid is " << maxparticleid << endl;

  for (int k = 0; k < NumParticles; k++)
  {

    DRT::Node* currentparticle = particledis_->lRowNode(k);  //input: LID!

    if(currentparticle->NumElement() != 1)
      dserror("ERROR: A particle is assigned to more than one bin!");
    DRT::Element** currelement = currentparticle->Elements();
  #ifdef DEBUG
      DRT::MESHFREE::MeshfreeMultiBin* test = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currelement[0]);
      if(test == NULL) dserror("dynamic cast from DRT::Element to DRT::MESHFREE::MeshfreeMultiBin failed");
  #endif

    DRT::MESHFREE::MeshfreeMultiBin* currbin = static_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currelement[0]);
    DRT::Element** scatraelesinbin = currbin->AssociatedFluidEles();
    int numscatraelesinbin = currbin->NumAssociatedFluidEle();

    std::set<int>::const_iterator eleiter;

    // search for underlying scatra element with standard search in case nothing was found
    for(int ele=0; ele<numscatraelesinbin; ++ele)
    {
      DRT::Element* scatraele = scatraelesinbin[ele];
      const DRT::Element::DiscretizationType distype = DRT::Element::hex8;
      const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      // get node coordinates of the current element
      static LINALG::Matrix<nsd,numnode> xyze;
      GEO::fillInitialPositionArray<distype>(scatraele, xyze);

      std::vector<int> lm_b = particledis_->Dof(currentparticle);
      int posx = disnp->Map().LID(lm_b[0]);
      for (int d=0; d<3; ++d)
        particleposition(d) = (*disnp)[posx+d];

      //get coordinates of the particle position in parameter space of the element
      bool insideele = GEO::currentToVolumeElementCoordinates(scatraele->Shape(), xyze, particleposition, elecoord, false);

      if(insideele == true)
      {
        targetscatraele = scatraele;

        // shape functions
        LINALG::Matrix<numnode,1> funct;

        //fill vector
        DRT::UTILS::shape_function_3D(funct,elecoord(0),elecoord(1),elecoord(2),targetscatraele->Shape());

        //Find values of Phi on 8 nodes of each element
        // extract local values from the global vectors
        // get element location vector and ownerships
        std::vector<int> lm;
        std::vector<int> lmowner;
        std::vector<int> lmstride;
        targetscatraele->LocationVector(*scatradis_,lm,lmowner,lmstride);

        const Teuchos::RCP<Epetra_MultiVector> vel_field = scatra_->ConVel();  //Address External Velocity Field
        Epetra_SerialDenseMatrix my_velocities(nsd, numnode);

        DRT::UTILS::ExtractMyNodeBasedValues(targetscatraele, my_velocities, vel_field, nsd); //fill my_velocities with x,y,z-velocities on nodes

        //convert SerialDenseMatrix to LINALG::Matrix
        LINALG::Matrix<nsd,numnode> my_velocities_in_matrix;
        for (int i = 0; i < nsd; i ++)
        {
          for (int j = 0; j < numnode; j ++)
          {
            my_velocities_in_matrix(i,j) = my_velocities(i,j);
          }
        }

        LINALG::Matrix<nsd,1> particle_velocity;
        particle_velocity.Multiply(my_velocities_in_matrix,funct); //(3x8) * (8x1) = (3x1)-Matrix (x_vel, y_vel, z_vel)^T

        //fill veln-vector with external velocity field values (on particle positions)
        for(int dim=0; dim<3; ++dim)
        {
          (*velnp)[3*k+dim] = particle_velocity(dim,0);
        }
      } //end: if inside element
    } //end: for all elements
  } //end: for all particles

#endif

  return;
}



Teuchos::RCP<Epetra_Vector> PARTICLE::ScatraParticleCoupling::CorrectionStep()
{
  //Loop over all particles to find those which have escaped
    Teuchos::RCP<Epetra_Vector> phinp_copy = Teuchos::rcp(new Epetra_Vector(*(scatradis_->DofRowMap()),true));
#if 0
  //Define pointers to vectors (which are defined in timint.cpp)
  Teuchos::RCP<Epetra_Vector> disn = particles_->ExtractDispnp();
  Teuchos::RCP<Epetra_Vector> radiusn = Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->ExtractRadiusnp();
  const Teuchos::RCP<Epetra_IntVector> signn = Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->ExtractSignVector();

  const Teuchos::RCP<const Epetra_Vector> phinp = scatra_->Phinp();

  Teuchos::RCP<Epetra_Vector> phinp_copy = Teuchos::rcp(new Epetra_Vector(*(scatradis_->DofRowMap()),true));
  phinp_copy->Update(1.0,*(scatra_->Phinp()),0.0);

  const int nsd = 3;

  LINALG::Matrix<nsd,1> particleposition;
  LINALG::Matrix<nsd,1> elecoord(true);
  DRT::Element* targetscatraele = NULL;

  int NumParticles = particledis_->NumMyRowNodes();
  std::cout << "Number of Nodes (=Particles) is " << NumParticles << endl;

  int maxparticleid = particledis_->NodeRowMap()->MaxAllGID();
  std::cout << "maxparticleid is " << maxparticleid << endl;

  //loop over all particles
  for (int k = 0; k < maxparticleid; k ++)
  {
    //extract values from disn-epetra-vector, fill into standard matrix
    std::vector <double> currparticlepos(nsd);
  //  LINALG::Matrix<nsd,1> currparticlepos(false);
    for(int dim=0; dim < nsd; ++dim)
    {
      currparticlepos[dim] = (*disn)[3*k+dim];
    }
    DRT::Node* currentparticle = particledis_->lRowNode(k);  //input: LID!

    if(currentparticle->NumElement() != 1)
      dserror("ERROR: A particle is assigned to more than one bin!");
    DRT::Element** currelement = currentparticle->Elements();
  #ifdef DEBUG
      DRT::MESHFREE::MeshfreeMultiBin* test = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currelement[0]);
      if(test == NULL) dserror("dynamic cast from DRT::Element to DRT::MESHFREE::MeshfreeMultiBin failed");
  #endif

    DRT::MESHFREE::MeshfreeMultiBin* currbin = static_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currelement[0]);
    DRT::Element** scatraelesinbin = currbin->AssociatedFluidEles();
    int numscatraelesinbin = currbin->NumAssociatedFluidEle();

    std::set<int>::const_iterator eleiter;
    // search for underlying scatra element with standard search in case nothing was found
    for(int ele=0; ele<numscatraelesinbin; ++ele)
    {
      DRT::Element* scatraele = scatraelesinbin[ele];
      const DRT::Element::DiscretizationType distype = DRT::Element::hex8;
      const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      // get node coordinates of the current element
      static LINALG::Matrix<nsd,numnode> currelecoords;
      GEO::fillInitialPositionArray<distype>(scatraele, currelecoords);

      std::vector<int> lm_b = particledis_->Dof(currentparticle);
      int posx = disn->Map().LID(lm_b[0]);
      for (int d=0; d<3; ++d)
        particleposition(d) = (*disn)[posx+d];

      //get coordinates of the particle position in parameter space of the element
      bool insideele = GEO::currentToVolumeElementCoordinates(scatraele->Shape(), currelecoords, particleposition, elecoord, false);

      if(insideele == true)
      {
        targetscatraele = scatraele;

        // shape functions
        LINALG::Matrix<numnode,1> funct;

        //fill vector
        DRT::UTILS::shape_function_3D(funct,elecoord(0),elecoord(1),elecoord(2),targetscatraele->Shape());

        //Find values of Phi on 8 nodes of each element
        // extract local values from the global vectors
        // get element location vector and ownerships
        std::vector<int> lm;
        std::vector<int> lmowner;
        std::vector<int> lmstride;
        targetscatraele->LocationVector(*scatradis_,lm,lmowner,lmstride);

        if (phinp_copy==Teuchos::null)
          dserror("Cannot get state vector 'hist' and/or 'phinp'");
        std::vector<double> myphinp(numnode);
        DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm); //fill phinp with phi values

        LINALG::Matrix<numnode,1> phi_vector_values;
        for (int s = 0; s < numnode; s++)
        {
          phi_vector_values(s,0) = myphinp[s];  //fill values in myphinp in new vector
        }

        //value of phi on current particle position
        LINALG::Matrix<1,1> phi_particle;
        phi_particle.MultiplyTN(funct,phi_vector_values);

        double phi_particle_position = phi_particle(0,0);

        //if a particle has escaped: address 8 nodes
//        if (phi_particle_position*((*signn)[k]) < 0)
//        {
//          double Krit = phi_particle_position*((*signn)[k]);
//          std::cout << "Kriterium ist " << Krit << endl;
//          //for all corners of the underlying element (in which particle is located)
//          for (int corner = 0; corner < numnode; corner++)
//          {
//            double local_phi = (*signn)[k]*((*radiusn)[k] - sqrt(
//                                        (currparticlepos[0]-currelecoords(0,corner))*(currparticlepos[0]-currelecoords(0,corner))
//                                      + (currparticlepos[1]-currelecoords(1,corner))*(currparticlepos[1]-currelecoords(1,corner))
//                                      + (currparticlepos[2]-currelecoords(2,corner))*(currparticlepos[2]-currelecoords(2,corner))));
//            if ((*signn)[k] > 0)
//            {
//              double phi_plus = max(local_phi, myphinp[corner]);
//              //renew phi on node: extract LID of node!
//              const int nodegid = (targetscatraele->Nodes()[corner])->Id();
//              const int lid = phinp_copy->Map().LID(nodegid);
//              if (abs(phi_plus < abs(myphinp[corner])))
//              {
//                (*phinp_copy)[lid] = phi_plus;
//                std::cout << "corrected!!!" << endl;
//              }
//            }
//
//            if ((*signn)[k] < 0)
//            {
//              double phi_minus = min(local_phi, myphinp[corner]);
//              //renew phi on node: extract LID of node!
//              const int nodegid = (targetscatraele->Nodes()[corner])->Id();
//              const int lid = phinp_copy->Map().LID(nodegid);
//              if (abs(phi_minus < abs(myphinp[corner])))
//              {
//                (*phinp_copy)[lid] = phi_minus;
//                std::cout << "corrected!!!" << endl;
//              }
//            }
//          }
//        }

        /// improved correction algorithm: Wang et al., Journal of Computational Physics
        if (phi_particle_position*((*signn)[k]) < 0)
        {
          double Krit = phi_particle_position*((*signn)[k]);
          std::cout << "Kriterium ist " << Krit << endl;
          //for all corners of the underlying element (in which particle is located)
          for (int corner = 0; corner < numnode; corner++)
          {
            double local_phi_minus = (*signn)[k]*((*radiusn)[k] - sqrt(
                                        (currparticlepos[0]-currelecoords(0,corner))*(currparticlepos[0]-currelecoords(0,corner))
                                      + (currparticlepos[1]-currelecoords(1,corner))*(currparticlepos[1]-currelecoords(1,corner))
                                      + (currparticlepos[2]-currelecoords(2,corner))*(currparticlepos[2]-currelecoords(2,corner))));
            double local_phi_plus = (*signn)[k]*((*radiusn)[k] + sqrt(
                                        (currparticlepos[0]-currelecoords(0,corner))*(currparticlepos[0]-currelecoords(0,corner))
                                      + (currparticlepos[1]-currelecoords(1,corner))*(currparticlepos[1]-currelecoords(1,corner))
                                      + (currparticlepos[2]-currelecoords(2,corner))*(currparticlepos[2]-currelecoords(2,corner))));
            if ((*signn)[k] > 0)
            {
              //renew phi on node: extract LID of node!
              const int nodegid = (targetscatraele->Nodes()[corner])->Id();
              const int lid = phinp_copy->Map().LID(nodegid);

              if (myphinp[corner] <= 0)
              {
                if ((abs(myphinp[corner]-(*phinp_copy)[lid]) < 0.001) || (abs(local_phi_minus) < abs(((*phinp_copy)[lid]))))
                  (*phinp_copy)[lid] = local_phi_minus;
              }
              if (myphinp[corner] > 0)
              {
                if ((abs(myphinp[corner]-(*phinp_copy)[lid]) < 0.001) || (abs(local_phi_plus) < abs(((*phinp_copy)[lid]))))
                  (*phinp_copy)[lid] = local_phi_plus;
              }
            }

            if ((*signn)[k] < 0)
            {
              //renew phi on node: extract LID of node!
              const int nodegid = (targetscatraele->Nodes()[corner])->Id();
              const int lid = phinp_copy->Map().LID(nodegid);

              if (myphinp[corner] <= 0)
              {
                if ((abs(myphinp[corner]-(*phinp_copy)[lid]) < 0.001) || (abs(local_phi_plus) < abs(((*phinp_copy)[lid]))))
                  (*phinp_copy)[lid] = local_phi_plus;
              }
              if (myphinp[corner] > 0)
              {
                if ((abs(myphinp[corner]-(*phinp_copy)[lid]) < 0.001) || (abs(local_phi_minus) < abs(((*phinp_copy)[lid]))))
                  (*phinp_copy)[lid] = local_phi_minus;
              }
            }
          }
        }

      } //end: if inside element
    } //end: for all elements
  } //end: for all particles

  //for each node: work out distance from escaped particle to node

  //compare this distance to local phi-value, decide if phi-value on that node has to be updated

  //if an escaped particle is detected in the same bin, we can easily proceed with our algorithm
 // double rsdst = 352.3;
# endif

  return phinp_copy;
}


/*----------------------------------------------------------------------*
 | transfer particles into their correct bins and update solution       |
 |                                                      rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::ScatraParticleCoupling::Update()
{
  // transfer particles into their correct bins
  // note: this function has to be called for each time step
  // TODO: einmal sollte reichen, ggf. nach correction und nochmal nach reseeding
  //Algorithm::TransferParticles();

  // update displacements, velocities, accelerations
  // after this call we will have disn_==dis_, etc
  // update time and step
  //Algorithm::Update();
  particles_->Update();

  return;
}


/*----------------------------------------------------------------------*
 | output particle time step                            rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::ScatraParticleCoupling::Output()
{
  // write data to file
  Algorithm::Output();
  return;
}


/*----------------------------------------------------------------------*
 | compute level-set value of particle with local id    rasthofer 12/13 |
 *----------------------------------------------------------------------*/
double PARTICLE::ScatraParticleCoupling::GetPhiParticle(
  const int particle_lid              ///< local particle id
  )
{
  //------------------------------------------------------------------------
  // element coordinates of particle position in scatra element
  //------------------------------------------------------------------------

  // get particle corresponding to lid
  DRT::Node* currparticle = particledis_->lRowNode(particle_lid);

  Teuchos::RCP<const Epetra_Vector> particlepos = particles_->Dispnp();
  // fill particle position
  LINALG::Matrix<3,1> particleposition;
  std::vector<int> lm_b = particledis_->Dof(currparticle);
  int posx = particlepos->Map().LID(lm_b[0]);
  for (int idim=0; idim<3; ++idim)
    particleposition(idim,0) = (*particlepos)[posx+idim];

  // variables to store information about element in which the particle is located
  DRT::Element* targetscatraele = NULL;
  LINALG::Matrix<3,1> elecoord(true);

  // find out in which scatra element the current particle is located
  if(currparticle->NumElement() != 1)
      dserror("ERROR: A particle is assigned to more than one bin!");

  DRT::Element** currele = currparticle->Elements();
#ifdef DEBUG
  DRT::MESHFREE::MeshfreeMultiBin* test = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currele[0]);
  if(test == NULL) dserror("dynamic cast from DRT::Element to DRT::MESHFREE::MeshfreeMultiBin failed");
#endif
  DRT::MESHFREE::MeshfreeMultiBin* currbin = static_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currele[0]);
  DRT::Element** scatraelesinbin = currbin->AssociatedFluidEles();
  int numscatraelesinbin = currbin->NumAssociatedFluidEle();

  std::set<int>::const_iterator eleiter;

  // search for underlying scatra element with standard search in case nothing was found
  for(int ele=0; ele<numscatraelesinbin; ++ele)
  {
    DRT::Element* scatraele = scatraelesinbin[ele];
    const LINALG::SerialDenseMatrix xyze(GEO::InitialPositionArray(scatraele));

    //get coordinates of the particle position in parameter space of the element
    bool insideele = GEO::currentToVolumeElementCoordinates(scatraele->Shape(), xyze, particleposition, elecoord, false);

    if(insideele == true)
    {
      targetscatraele = scatraele;
      // leave loop over all scatra eles in bin
      break;
    }
  }

  if(targetscatraele == NULL) dserror("Could not found corresponding element!");

  //------------------------------------------------------------------------
  // compute level-set value
  //------------------------------------------------------------------------

  if (targetscatraele->Shape() != DRT::Element::hex8)
    dserror("Other element than hex8 not yet supported!");
  const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;

  // get phinp
  const Teuchos::RCP<const Epetra_Vector> phinp = scatra_->Phinp();

  // get nodal phi values
  std::vector<int> lm;
  std::vector<int> lmowner;
  std::vector<int> lmstride;
  lm.clear();
  lmowner.clear();
  lmstride.clear();
  targetscatraele->LocationVector(*scatradis_,lm,lmowner,lmstride);
  std::vector<double> myphinp(lm.size());
  DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);
  // fill all element arrays
  LINALG::Matrix<numnode,1> ephinp(true);
  for (int i=0;i<numnode;++i)
  {
    // split for each transported scalar, insert into element arrays
    ephinp(i,0) = myphinp[i];
  } // for i

  // shape functions
  LINALG::Matrix<numnode,1> funct(true);
  //fill vectors
  DRT::UTILS::shape_function_3D(funct,elecoord(0),elecoord(1),elecoord(2),targetscatraele->Shape());

  // finally compute phi
  double phi = funct.Dot(ephinp);

  LINALG::Matrix<3,numnode> xyze(true);
  GEO::fillInitialPositionArray<DRT::Element::hex8,3,LINALG::Matrix<3,numnode> >(targetscatraele,xyze);

  LINALG::Matrix<3,numnode> deriv(true);
  DRT::UTILS::shape_function_3D_deriv1(deriv,elecoord(0),elecoord(1),elecoord(2),targetscatraele->Shape());
  // get transposed of the jacobian matrix d x / d \xi
  // xjm(i,j) = deriv(i,k)*xyze(j,k)
  static LINALG::Matrix<3,3> xjm(true);
  xjm.MultiplyNT(deriv,xyze);
  const double det = xjm.Determinant();
  if (det < 0.0)
    dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", targetscatraele->Id(), det);
  // inverse of jacobian
  LINALG::Matrix<3,3> xji(true);
  xji.Invert(xjm);
  // compute global derivates
  LINALG::Matrix<3,numnode> derxy(true);
  // derxy(i,j) = xji(i,k) * deriv(k,j)
  derxy.Multiply(xji,deriv);

  // get gradient of phi
  LINALG::Matrix<3,1> gradphi(true);
  gradphi.Multiply(derxy,ephinp);

  return phi;
}


/*----------------------------------------------------------------------*
 | get underlying element as well as particle position   rasthofer 12/13 |
 *----------------------------------------------------------------------*/
DRT::Element* PARTICLE::ScatraParticleCoupling::GetEleCoordinates(
  const int particle_lid,             ///< local particle id
  LINALG::Matrix<3,1>& elecoord       ///< matrix to be filled with particle coordinates in element space
  )
{
  // get particle corresponding to lid
  DRT::Node* currparticle = particledis_->lRowNode(particle_lid);

  Teuchos::RCP<const Epetra_Vector> particlepos = particles_->Dispnp();
  // fill particle position
  LINALG::Matrix<3,1> particleposition;
  std::vector<int> lm_b = particledis_->Dof(currparticle);
  int posx = particlepos->Map().LID(lm_b[0]);
  for (int idim=0; idim<3; ++idim)
    particleposition(idim,0) = (*particlepos)[posx+idim];

  // variables to store information about element in which the particle is located
  DRT::Element* targetscatraele = NULL;
  elecoord.Clear();

  // find out in which scatra element the current particle is located
  if(currparticle->NumElement() != 1)
      dserror("ERROR: A particle is assigned to more than one bin!");

  DRT::Element** currele = currparticle->Elements();
#ifdef DEBUG
  DRT::MESHFREE::MeshfreeMultiBin* test = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currele[0]);
  if(test == NULL) dserror("dynamic cast from DRT::Element to DRT::MESHFREE::MeshfreeMultiBin failed");
#endif
  DRT::MESHFREE::MeshfreeMultiBin* currbin = static_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currele[0]);
  DRT::Element** scatraelesinbin = currbin->AssociatedFluidEles();
  int numscatraelesinbin = currbin->NumAssociatedFluidEle();

  std::set<int>::const_iterator eleiter;

  // search for underlying scatra element with standard search in case nothing was found
  for(int ele=0; ele<numscatraelesinbin; ++ele)
  {
    DRT::Element* scatraele = scatraelesinbin[ele];
    const LINALG::SerialDenseMatrix xyze(GEO::InitialPositionArray(scatraele));

    //get coordinates of the particle position in parameter space of the element
    bool insideele = GEO::currentToVolumeElementCoordinates(scatraele->Shape(), xyze, particleposition, elecoord, false);

    if(insideele == true)
    {
      targetscatraele = scatraele;
      // leave loop over all scatra eles in bin
      break;
    }
  }

  if(targetscatraele == NULL) dserror("Could not found corresponding element!"); // TODO bei erweitertem gebiet sollte das bei partikel ausserhalb auftreten

  return targetscatraele;
}


/*----------------------------------------------------------------------*
 | compute level-set value of particle with local id    rasthofer 12/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::ScatraParticleCoupling::GetLSValParticle(
  double& phi_particle,               ///< phi value at particle position
  LINALG::Matrix<3,1>& gradphi,       ///< gradient of level-set field at particle position
  DRT::Element* scatraele,            ///< corresponding scatra element
  const LINALG::Matrix<3,1> elecoord  ///< matrix containing particle coordinates in element space
  )
{
  if (scatraele->Shape() != DRT::Element::hex8)
    dserror("Other element than hex8 not yet supported!");
  const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;

  // get phinp
  const Teuchos::RCP<const Epetra_Vector> phinp = scatra_->Phinp();

  // get nodal phi values
  std::vector<int> lm;
  std::vector<int> lmowner;
  std::vector<int> lmstride;
  lm.clear();
  lmowner.clear();
  lmstride.clear();
  scatraele->LocationVector(*scatradis_,lm,lmowner,lmstride);
  std::vector<double> myphinp(lm.size());
  DRT::UTILS::ExtractMyValues(*phinp,myphinp,lm);
  // fill all element arrays
  LINALG::Matrix<numnode,1> ephinp(true);
  for (int i=0;i<numnode;++i)
  {
    // split for each transported scalar, insert into element arrays
    ephinp(i,0) = myphinp[i];
  } // for i

  // shape functions
  LINALG::Matrix<numnode,1> funct(true);
  //fill vectors
  DRT::UTILS::shape_function_3D(funct,elecoord(0),elecoord(1),elecoord(2),scatraele->Shape());

  // finally compute phi
  phi_particle = funct.Dot(ephinp);

  LINALG::Matrix<3,numnode> xyze(true);
  GEO::fillInitialPositionArray<DRT::Element::hex8,3,LINALG::Matrix<3,numnode> >(scatraele,xyze);

  LINALG::Matrix<3,numnode> deriv(true);
  DRT::UTILS::shape_function_3D_deriv1(deriv,elecoord(0),elecoord(1),elecoord(2),scatraele->Shape());
  // get transposed of the jacobian matrix d x / d \xi
  // xjm(i,j) = deriv(i,k)*xyze(j,k)
  static LINALG::Matrix<3,3> xjm(true);
  xjm.MultiplyNT(deriv,xyze);
  const double det = xjm.Determinant();
  if (det < 0.0)
    dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", scatraele->Id(), det);
  // inverse of jacobian
  LINALG::Matrix<3,3> xji(true);
  xji.Invert(xjm);
  // compute global derivates
  LINALG::Matrix<3,numnode> derxy(true);
  // derxy(i,j) = xji(i,k) * deriv(k,j)
  derxy.Multiply(xji,deriv);

  // get gradient of phi
  gradphi.Clear();
  gradphi.Multiply(derxy,ephinp);

  return;
}


/*----------------------------------------------------------------------*
| setup ghosting of bins, particles & underlying fluid      ghamm 11/12 |
 *----------------------------------------------------------------------*/
void PARTICLE::ScatraParticleCoupling::SetupGhosting(Teuchos::RCP<Epetra_Map> binrowmap, std::map<int, std::set<int> >& scatraelesinbins)
{
  //--------------------------------------------------------------------
  // 1st step: no ghosting for bins and particles
  //--------------------------------------------------------------------

  // Note: IndependentDofSet is used; new dofs are numbered from zero, minnodgid is ignored and it does not register in static_dofsets_
  particledis_->FillComplete(true, false, true);

  // set bin row map to be equal to bin col map because we do not have ghosting
  bincolmap_ = binrowmap;

  //--------------------------------------------------------------------
  // 2nd step: extend ghosting of underlying scatra discretization according to bin distribution
  //--------------------------------------------------------------------
  std::map<int, std::set<int> > extendedscatraghosting;
  {
    // do communication to gather all elements for extended ghosting
    const int numproc = scatradis_->Comm().NumProc();

    for (int iproc = 0; iproc < numproc; ++iproc)
    {
      // first: proc i tells all procs how many row bins it has
      int numbin = binrowmap->NumMyElements();
      scatradis_->Comm().Broadcast(&numbin, 1, iproc);
      // second: proc i tells all procs which row bins it has
      std::vector<int> binid(numbin,0);
      if(iproc == myrank_)
      {
        int* binrowmapentries = binrowmap->MyGlobalElements();
        for (int i=0; i<numbin; ++i)
          binid[i] = binrowmapentries[i];
      }
      scatradis_->Comm().Broadcast(&binid[0], numbin, iproc);

      // loop over all own bins and find requested ones
      std::map<int, std::set<int> > sdata;
      std::map<int, std::set<int> > rdata;

      for(int i=0; i<numbin; ++i)
      {
        sdata[binid[i]].insert(scatraelesinbins[binid[i]].begin(),scatraelesinbins[binid[i]].end());
      }

      LINALG::Gather<int>(sdata, rdata, 1, &iproc, scatradis_->Comm());

      // proc i has to store the received data
      if(iproc == myrank_)
      {
        extendedscatraghosting = rdata;
      }
    }

    // reduce map of sets to one set and copy to a vector to create extended scatra colmap
    std::set<int> reduscatraeleset;
    std::map<int, std::set<int> >::iterator iter;
    for(iter=extendedscatraghosting.begin(); iter!= extendedscatraghosting.end(); ++iter)
    {
      reduscatraeleset.insert(iter->second.begin(),iter->second.end());
    }
    std::vector<int> scatracolgids(reduscatraeleset.begin(),reduscatraeleset.end());
    Teuchos::RCP<Epetra_Map> scatracolmap = Teuchos::rcp(new Epetra_Map(-1,(int)scatracolgids.size(),&scatracolgids[0],0,Comm()));

    scatradis_->ExtendedGhosting(*scatracolmap,true,true,true,false);

  }

  //--------------------------------------------------------------------
  // 4th step: assign scatra elements to bins
  //--------------------------------------------------------------------
  {
    for(std::map<int, std::set<int> >::const_iterator biniter=extendedscatraghosting.begin(); biniter!=extendedscatraghosting.end(); ++biniter)
    {
      DRT::MESHFREE::MeshfreeMultiBin* currbin = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(particledis_->gElement(biniter->first));
      for(std::set<int>::const_iterator scatraeleiter=biniter->second.begin(); scatraeleiter!=biniter->second.end(); ++scatraeleiter)
      {
        int scatraeleid = *scatraeleiter;
        currbin->AddAssociatedFluidEle(scatraeleid, scatradis_->gElement(scatraeleid));
//          cout << "in bin with id:" << currbin->Id() << " is fluid ele with id" << fluideleid << "with pointer" << fluiddis_->gElement(fluideleid) << endl;
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*
| build connectivity from fluid elements to bins            ghamm 07/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::ScatraParticleCoupling::BuildElementToBinPointers()
{
  // loop over column bins and fill scatra elements
  const int numcolbin = particledis_->NumMyColElements();
  for (int ibin=0; ibin<numcolbin; ++ibin)
  {
    DRT::Element* actele = particledis_->lColElement(ibin);
    DRT::MESHFREE::MeshfreeMultiBin* actbin = static_cast<DRT::MESHFREE::MeshfreeMultiBin*>(actele);
    const int numfluidele = actbin->NumAssociatedFluidEle();
    const int* fluideleids = actbin->AssociatedFluidEleIds();
    std::vector<DRT::Element*> fluidelements(numfluidele);
    for(int iele=0; iele<numfluidele; ++iele)
    {
      const int fluideleid = fluideleids[iele];
      fluidelements[iele] = scatradis_->gElement(fluideleid);
    }
    actbin->BuildFluidElePointers(&fluidelements[0]);
  }

  return;
}
