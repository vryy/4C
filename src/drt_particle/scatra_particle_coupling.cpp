/*----------------------------------------------------------------------*/
/*! \file

\level 3

\brief Algorithm to track particles for level-set problems

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
*----------------------------------------------------------------------*/

#include "particle_algorithm.H"
#include "scatra_particle_coupling.H"
#include "particle_timint_rk.H"
#include "../drt_adapter/adapter_particle.H"
#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_levelset/levelset_algorithm.H"
#include "../drt_binstrategy/drt_meshfree_multibin.H"
#include "../drt_binstrategy/binning_strategy_utils.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_parallel.H"

#include "../drt_mortar/mortar_utils.H"
#include "../drt_mortar/mortar_coupling3d_classes.H"

#include "../drt_geometry/searchtree.H"
#include "../drt_geometry/searchtree_geometry_service.H"
#include "../drt_geometry/element_coordtrafo.H"
#include "../drt_geometry/position_array.H"
#include "../drt_geometry/element_volume.H"
#include "../drt_cut/cut_position.H"
#include "../linalg/linalg_utils.H"

#include "../drt_io/io_pstream.H"
#include "../headers/definitions.h"

#include <Epetra_FEVector.h>
#include <Teuchos_TimeMonitor.hpp>

#define MULTIPLE_ELE

// TODO: if 0 and ifndef MULTIPLE_ELE (i.e., gray) code is outdated and may be removed
//       I keep it for the time being to have a backup!

/*----------------------------------------------------------------------*
 | algorithm constructor                                rasthofer 09/13 |
 *----------------------------------------------------------------------*/
PARTICLE::ScatraParticleCoupling::ScatraParticleCoupling(
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra, Teuchos::RCP<Teuchos::ParameterList> params)
    : PARTICLE::Algorithm(scatra->Discretization()->Comm(), *params),
      bin_transpcontent_(BINSTRATEGY::UTILS::Scatra),
      scatra_(scatra),
      scatradis_(scatra->Discretization()),
      params_(params),
      escaped_(DRT::INPUT::IntegralValue<INPAR::PARTICLEOLD::Escaped>(
          params->sublist("PARTICLE"), "ESCAPED")),
      reseeding_(params->sublist("PARTICLE").get<int>("RESEEDING")),
      fast_(DRT::INPUT::IntegralValue<bool>(params->sublist("PARTICLE"), "FAST_CHECK")),
      delete_more_(params->sublist("PARTICLE").get<double>("DELETE_CRITICAL_PARTICLES"))
{
  if (rep_strategy_ != INPAR::PARTICLEOLD::repstr_everydt)
    dserror("REPARTITIONSTRATEGY must be set to Everydt");

  Init(false);

  return;
}


/*----------------------------------------------------------------------*
 | initialization of the system                         rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::ScatraParticleCoupling::Init(bool restarted)
{
  // -------------------------------------------------------------------
  //               setup particle discretization
  // -------------------------------------------------------------------

  // this is done here this way to avoid overloading of ReadRestart()
  if (restarted)
  {
    // FillComplete() necessary for DRT::Geometry .... could be removed perhaps
    BinStrategy()->BinDiscret()->FillComplete(false, false, false);

    std::list<Teuchos::RCP<DRT::Node>> homelessparticles;
    Teuchos::RCP<Epetra_Map> particlerowmap =
        Teuchos::rcp(new Epetra_Map(*BinStrategy()->BinDiscret()->NodeRowMap()));
    for (int lid = 0; lid < particlerowmap->NumMyElements(); ++lid)
    {
      DRT::Node* node = BinStrategy()->BinDiscret()->gNode(particlerowmap->GID(lid));
      const double* currpos = node->X();
      PlaceNodeCorrectly(Teuchos::rcp(node, false), currpos, homelessparticles);
    }

    // start round robin loop to fill particles into their correct bins
    FillParticlesIntoBinsRoundRobin(homelessparticles);

    BinStrategy()->BinDiscret()->FillComplete(true, false, true);
    return;
  }

  // FillComplete() necessary for DRT::Geometry .... could be removed perhaps
  BinStrategy()->BinDiscret()->FillComplete(false, false, false);
  // extract noderowmap because it will be called Reset() after adding elements
  Teuchos::RCP<Epetra_Map> particlerowmap =
      Teuchos::rcp(new Epetra_Map(*BinStrategy()->BinDiscret()->NodeRowMap()));
  BinStrategy()->CreateBinsScatra(scatradis_);

  // setup pbcs after bins have been created
  BinStrategy()->BuildPeriodicBC();

  // gather all scatra col eles in each bin for proper extended ghosting
  std::map<int, std::set<int>> scatraelesinbins;
  Teuchos::RCP<Epetra_Map> binrowmap =
      DistributeBinsToProcsBasedOnUnderlyingDiscret(scatradis_, scatraelesinbins);

  //--------------------------------------------------------------------
  // -> 1) create a list of homeless particles that are not in a bin on this proc
  std::list<Teuchos::RCP<DRT::Node>> homelessparticles;

  for (int lid = 0; lid < particlerowmap->NumMyElements(); ++lid)
  {
    DRT::Node* node = BinStrategy()->BinDiscret()->gNode(particlerowmap->GID(lid));
    const double* currpos = node->X();
    PlaceNodeCorrectly(Teuchos::rcp(node, false), currpos, homelessparticles);
  }

  // start round robin loop to fill particles into their correct bins
  FillParticlesIntoBinsRoundRobin(homelessparticles);

  // ghost bins, particles and fluid elements according to the bins
  SetupGhosting(binrowmap, scatraelesinbins);

  // some output
  if (MyRank() == 0) IO::cout << "after ghosting" << IO::endl;
  DRT::UTILS::PrintParallelDistribution(*BinStrategy()->BinDiscret());
  DRT::UTILS::PrintParallelDistribution(*scatradis_);

  // -------------------------------------------------------------------
  //               setup time integration
  // -------------------------------------------------------------------

  // create time integrator
  Teuchos::ParameterList timintpara;
  timintpara.set<double>("TIMESTEP", params_->get<double>("TIMESTEP"));
  timintpara.set<int>("NUMSTEP", params_->get<int>("NUMSTEP"));
  timintpara.set<int>("RESTARTEVRY", params_->get<int>("RESTARTEVRY"));
  // for large problems, it might be useful to only write the restart data
  if (DRT::INPUT::IntegralValue<bool>(params_->sublist("PARTICLE"), "RESTARTDATA"))
    timintpara.set<int>("RESULTSEVRY", params_->get<int>("RESTARTEVRY"));
  else
    timintpara.set<int>("RESULTSEVRY", params_->get<int>("RESULTSEVRY"));
  Teuchos::RCP<ADAPTER::ParticleBaseAlgorithm> particles =
      Teuchos::rcp(new ADAPTER::ParticleBaseAlgorithm(timintpara, BinStrategy()->BinDiscret()));
  particles_ = particles->ParticleField();
  // set particle algorithm into time integration
  particles_->SetParticleAlgorithm(Teuchos::rcp(this, false));

  // initialize particle field
  particles_->Init();

  // -------------------------------------------------------------------
  //               initialize all members
  // -------------------------------------------------------------------

  // get number of particles of each sign per bin around interface
  num_particles_per_bin_ = params_->sublist("PARTICLE").get<int>("NUMPARTICLE");

  // get characteristic element length of scatra discretization
  // assumed equal to bin edge length
  //  std::cout << "number of bins " << BinStrategy()->BinDiscret()->NumMyRowElements() <<
  //  std::endl; std::cout << bin_size_[0] << "  " << bin_size_[1] << "  " << bin_size_[2] <<
  //  std::endl; if (std::abs(bin_size_[0]-bin_size_[1]) > 10e-9 or
  //  std::abs(bin_size_[0]-bin_size_[2]) > 10e-9)
  //    dserror("Cubic bins expected");

  // get bin edge length
  // const double binlength = cutoff_radius_; // should be equal to bin_size_[0] for the present
  // settings: only valid for cubic bins
  switch (BinStrategy()->ParticleDim())
  {
    case INPAR::PARTICLEOLD::particle_3D:
    {
      // get maximal bin edge length
      binlength_max_ = std::max(BinStrategy()->BinSize()[0], BinStrategy()->BinSize()[1]);
      binlength_max_ = std::max(binlength_max_, BinStrategy()->BinSize()[2]);
      // get minimal bin edge length
      binlength_min_ = std::min(BinStrategy()->BinSize()[0], BinStrategy()->BinSize()[1]);
      binlength_min_ = std::min(binlength_min_, BinStrategy()->BinSize()[2]);
      break;
    }
    case INPAR::PARTICLEOLD::particle_2Dx:
    {
      // get maximal bin edge length
      binlength_max_ = std::max(BinStrategy()->BinSize()[1], BinStrategy()->BinSize()[2]);
      // get minimal bin edge length
      binlength_min_ = std::min(BinStrategy()->BinSize()[1], BinStrategy()->BinSize()[2]);
      break;
    }
    case INPAR::PARTICLEOLD::particle_2Dy:
    {
      // get maximal bin edge length
      binlength_max_ = std::max(BinStrategy()->BinSize()[0], BinStrategy()->BinSize()[2]);
      // get minimal bin edge length
      binlength_min_ = std::min(BinStrategy()->BinSize()[0], BinStrategy()->BinSize()[2]);
      break;
    }
    case INPAR::PARTICLEOLD::particle_2Dz:
    {
      // get maximal bin edge length
      binlength_max_ = std::max(BinStrategy()->BinSize()[0], BinStrategy()->BinSize()[1]);
      // get minimal bin edge length
      binlength_min_ = std::min(BinStrategy()->BinSize()[0], BinStrategy()->BinSize()[1]);
      break;
    }
    default:
      break;
  }
  // write info to screen
  if (MyRank() == 0)
    std::cout << "Bin info: b_min " << std::setprecision(12) << binlength_min_ << " b_max "
              << binlength_max_ << std::endl;

  // in the following, the variables are named as given in Enright et al. 2002
  // define band of elements around interface to be filled with particles, i.e., band of size
  // bandwith * h on both sides of the interface
  const double bandwith = params_->sublist("PARTICLE").get<double>("PARTICLEBANDWIDTH");
  b_max_ = bandwith * binlength_max_;
  // get minimal radius of particles
  const double r_min_fac = params_->sublist("PARTICLE").get<double>("MIN_RADIUS");
  r_min_ = r_min_fac * binlength_min_;
  // get maximal radius of particles
  const double r_max_fac = params_->sublist("PARTICLE").get<double>("MAX_RADIUS");
  r_max_ = r_max_fac * binlength_min_;

  return;
}


/*----------------------------------------------------------------------*
 | initialize particle field with particles             rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::ScatraParticleCoupling::InitialSeeding()
{
  // remark: Since the particle Id depends on the number of processors and various random numbers
  //         are used, this function cannot be tested in parallel by usual means.
  //         I thus carefully checked this function as far as possible also in parallel.
  //         The following tests were performed:
  //         - for a plane interface, particles were seeded at the bin center and then attracted to
  //         a fixed
  //           target level-set value, this position was checked exactly by values
  //         - for a circular interface, particles were seeded at the bin center and then attracted
  //         to a fixed
  //           target level-set value, the positions were the visually compared in paraview and
  //           matched exactly
  //         - the correct numbering of the particles by gids was checked via screen output
  //         - no segmentations faults were obtained for up to 11 procs
  //         - the final particle field for zalesaks disk was checked in paraview for up to 5 procs,
  //         particle
  //           distribution in paraview was as expected
  //         - for a plane interface, 1 attraction step is need as expected for all number of procs,
  //         likewise the
  //           circular interface showed 5 required steps to attract all particles, for zalesaks
  //           disk, an almost constant number of particles for all procs remained after 15 steps

  // -------------------------------------------------------------------
  //               setup initial seeding
  // -------------------------------------------------------------------

  if (MyRank() == 0) std::cout << "-------- INITIAL SEEDING --------- " << std::endl;

  // -------------------------------------------------------------------
  //               find bins to be filled
  // -------------------------------------------------------------------

  // vector containing id of bins (equal to elements) which will be filled with particles
  std::vector<int> particlebins;

  // get phinp
  const Teuchos::RCP<const Epetra_Vector> row_phinp = scatra_->Phinp();
  // export phi vector to col map
  // this is necessary here, since the elements adjacent to a row node may have
  // ghosted nodes, which are only contained in the col format
  const Epetra_Map* scatra_dofcolmap = scatradis_->DofColMap();
  Teuchos::RCP<Epetra_Vector> phinp = Teuchos::rcp(new Epetra_Vector(*scatra_dofcolmap));
  LINALG::Export(*row_phinp, *phinp);

  // loop all row bins on this proc
  for (int ibin = 0; ibin < BinStrategy()->BinDiscret()->NumMyRowElements(); ibin++)
  {
    // get pointer to current bin
    DRT::Element* actele = BinStrategy()->BinDiscret()->lRowElement(ibin);
    DRT::MESHFREE::MeshfreeMultiBin* actbin =
        dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(actele);

    // get pointer to associated scatra elements
    DRT::Element** scatraelesinbin = actbin->AssociatedEles(bin_transpcontent_);
    // check for null-pointer in case of holes in the domain
    if (scatraelesinbin == NULL) continue;

#ifndef MULTIPLE_ELE

    // bool to go to next bin if this bin is to far away from the interface
    // or as soon as this bin is identified to be sufficiently close to set particles
    bool next_bin = false;

    // loop all elements in bin
    for (int iele = 0; iele < actbin->NumAssociatedEle(bin_transpcontent_); iele++)
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
      scatraele->LocationVector(*scatradis_, lm, lmowner, lmstride);
      std::vector<double> myphinp(lm.size());
      DRT::UTILS::ExtractMyValues(*phinp, myphinp, lm);

      // check, if bin is close to the interface
      if (std::abs(myphinp[0]) > (10.0 * b_max_))
      {
        // bin is located away from the interface
        // abort loop elements and go to the next bin
        break;  // loop all elements
      }

      // get element center
      const DRT::Element::DiscretizationType distype = DRT::Element::hex8;
      const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      // get center coordinates of current element
      LINALG::Matrix<3, 1> elecenter(true);
      LINALG::Matrix<3, numnode> xyzele(true);
      GEO::fillInitialPositionArray<distype>(scatraele, xyzele);

      for (int inode = 0; inode < numnode; inode++)
        for (int idim = 0; idim < 3; idim++) elecenter(idim, 0) += xyzele(idim, inode);

      const double inv_numnode = 1.0 / ((double)numnode);
      for (int idim = 0; idim < 3; idim++) elecenter(idim, 0) *= inv_numnode;

      // assuming bin length to be chosen such that only one element has its center inside the bin,
      // we only check the distance of this element to the interface
      // get center of current bin
      LINALG::Matrix<3, 1> center(true);
      center = GetBinCentroid(actbin->Id());

      if ((elecenter(0, 0) > (center(0, 0) - (bin_size_[0] / 2.0))) and
          (elecenter(0, 0) < (center(0, 0) + (bin_size_[0] / 2.0))))
      {
        if ((elecenter(1, 0) > (center(1, 0) - (bin_size_[1] / 2.0))) and
            (elecenter(1, 0) < (center(1, 0) + (bin_size_[1] / 2.0))))
        {
          if ((elecenter(2, 0) > (center(2, 0) - (bin_size_[2] / 2.0))) and
              (elecenter(2, 0) < (center(2, 0) + (bin_size_[2] / 2.0))))
          {
            // loop all nodes of this element
            for (int inode = 0; inode < numnode; inode++)
              if (std::abs(myphinp[inode]) < b_max_)
              {
                // bin is close to the interface and particles should be set
                // add bin to list
                particlebins.push_back(actbin->Id());
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
      if (next_bin == true) break;  // loop all elements

    }  // end loop all elements
#else

    // coarse check
    if (fast_)
    {
      // if first node of first associated element is far from the interface the bin is also
      // expected to be far from the interface
      DRT::Element* scatraele = scatraelesinbin[0];
      if (scatraele->Shape() != DRT::Element::hex8)
        dserror("Other element than hex8 not yet supported!");
      // This check should also work for all other element types, but has not yet been tested or
      // applied! Thus, check first! In general, you should first test entire algorithm for other
      // element types than hex8 before using it! But should work!

      // get nodal phi values
      std::vector<int> lm;
      std::vector<int> lmowner;
      std::vector<int> lmstride;
      lm.clear();
      lmowner.clear();
      lmstride.clear();
      scatraele->LocationVector(*scatradis_, lm, lmowner, lmstride);
      std::vector<double> myphinp(lm.size());
      DRT::UTILS::ExtractMyValues(*phinp, myphinp, lm);

      // check, if bin is close to the interface
      if (std::abs(myphinp[0]) > (10.0 * b_max_))
      {
        // bin is located away from the interface
        // go to the next bin
        continue;
      }
    }

    // get bin corners
    std::vector<LINALG::Matrix<3, 1>> actbincorner;
    BinStrategy()->GetBinCorners(actbin->Id(), actbincorner);

    for (std::size_t icor = 0; icor < actbincorner.size(); icor++)
    {
      // get element coordinates of bin corner and associated scatra element
      LINALG::Matrix<3, 1> elecoord(true);
      DRT::Element* scatraele = GetEleCoordinatesFromPosition(actbincorner[icor], actbin, elecoord);

      // compute level-set value
      // get phi and gradient of particle at current particle position
      double phi_corner = 0.0;
      // we do not need the gradient here, so this is just a dummy
      LINALG::Matrix<3, 1> normal_particle(true);
      GetLSValParticle(phi_corner, normal_particle, scatraele, elecoord, phinp, false);

      // is this corner in band
      if (std::abs(phi_corner) < b_max_)
      {
        // bin is close to the interface and particles should be set
        // add bin to list
        particlebins.push_back(actbin->Id());
        break;
      }
    }
#endif

  }  // end loop all bins

  // output for testing
  // std::cout << "myrank  " << MyRank() << "bins with particles " << particlebins.size() <<
  // std::endl;
  //  for (std::size_t i=0; i<particlebins.size(); i++)
  //    std::cout << particlebins[i] << std::endl;

  // -------------------------------------------------------------------
  //               seed particles in selected bins
  // -------------------------------------------------------------------

  // initialize particle id with largest particle id in use + 1 (on each proc)
  int maxparticleid = BinStrategy()->BinDiscret()->NodeRowMap()->MaxAllGID();
  int currentparticleid = maxparticleid;

  // compute local offset for global id numbering
  int myoffset = 0;

  // communicate number of bins of each proc to all procs
  const int numproc = BinStrategy()->BinDiscret()->Comm().NumProc();

  for (int iproc = 0; iproc < numproc; ++iproc)
  {
    int numbins = particlebins.size();
    BinStrategy()->BinDiscret()->Comm().Broadcast(&numbins, 1, iproc);
    if (MyRank() > iproc) myoffset += numbins;
  }
  // std::cout << "myrank  " << MyRank() << "   offset  " << myoffset << std::endl;

  // some output
  int mynumbins = particlebins.size();
  int allnumbins = 0;
  BinStrategy()->BinDiscret()->Comm().SumAll(&mynumbins, &allnumbins, 1);
  if (MyRank() == 0)
    std::cout << "--- total number of expected particles  "
              << (2 * num_particles_per_bin_ * allnumbins) << std::endl;

  // added offset to current particle id
  currentparticleid += (2 * num_particles_per_bin_ * myoffset);

  // list of all particles and their position
  std::map<int, std::vector<double>> particle_positions;
  // list of all particles and their sign
  std::map<int, double> particle_sign;

  // loop all bins which should be filled
  for (std::size_t ibin = 0; ibin < particlebins.size(); ibin++)
  {
    // get center of current bin
    LINALG::Matrix<3, 1> center(true);
    center = BinStrategy()->GetBinCentroid(particlebins[ibin]);

    // get max and min values of bin corners
    std::vector<LINALG::Matrix<3, 1>> bincorners(2);
    for (int rr = 0; rr < 2; rr++)
    {
      for (int idim = 0; idim < 3; idim++)
        (bincorners[rr])(idim, 0) =
            center(idim, 0) + pow(-1.0, (double)rr + 1.0) * 0.5 * BinStrategy()->BinSize()[idim];
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
    for (int ipart = 0; ipart < (2 * num_particles_per_bin_); ipart++)
    {
      // set particle id
      currentparticleid += 1;

      // get position randomly in bin
      std::vector<double> position(3);

      DRT::UTILS::Random* random = DRT::Problem::Instance()->Random();

      for (int idim = 0; idim < 3; idim++)
      {
        // set range (default: [-1;1])
        random->SetRandRange((bincorners[0])(idim, 0), (bincorners[1])(idim, 0));
        // get position
        position[idim] = random->Uni();
      }

      // set particle to mid plane in case of quasi-2D simulations
      switch (BinStrategy()->ParticleDim())
      {
        case INPAR::PARTICLEOLD::particle_2Dx:
        {
          position[0] = 0.0;
          break;
        }
        case INPAR::PARTICLEOLD::particle_2Dy:
        {
          position[1] = 0.0;
          break;
        }
        case INPAR::PARTICLEOLD::particle_2Dz:
        {
          position[2] = 0.0;
          break;
        }
        default:
          break;
      }

      // place particle with id and position
      std::list<Teuchos::RCP<DRT::Node>> homelessparticles;
      Teuchos::RCP<DRT::Node> newparticle =
          Teuchos::rcp(new DRT::Node(currentparticleid, &position[0], MyRank()));
      PlaceNodeCorrectly(newparticle, &position[0], homelessparticles);
      if (homelessparticles.size() != 0)
        dserror("New particle could not be inserted on this proc!.");

      // add particle to list
      particle_positions.insert(std::pair<int, std::vector<double>>(currentparticleid, position));
      // set sign (+ for the first half, - for the second one)
      if (ipart < num_particles_per_bin_)
        particle_sign.insert(std::pair<int, double>(currentparticleid, 1.0));
      else
        particle_sign.insert(std::pair<int, double>(currentparticleid, -1.0));
    }
  }

  // rebuild connectivity and assign degrees of freedom (note: IndependentDofSet)
  BinStrategy()->BinDiscret()->FillComplete(true, false, true);

  // update of state vectors to the new maps
  particles_->UpdateStatesAfterParticleTransfer();

  // -------------------------------------------------------------------
  //               initialize state vectors
  // -------------------------------------------------------------------

  // get position vector at time_(n+1)
  Teuchos::RCP<Epetra_Vector> disnp = particles_->WriteAccessDispnp();
  // get sign vector
  Teuchos::RCP<Epetra_Vector> sign =
      Teuchos::rcp_dynamic_cast<PARTICLE::TimIntRK>(particles_)->WriteAccessSign();

  // get maps
  const Epetra_Map* dofrowmap = BinStrategy()->BinDiscret()->DofRowMap();
  const Epetra_Map* noderowmap = BinStrategy()->BinDiscret()->NodeRowMap();

  for (int inode = 0; inode < BinStrategy()->BinDiscret()->NumMyRowNodes(); inode++)
  {
    // get global id of current particle
    const int gid = noderowmap->GID(inode);

    DRT::Node* currparticle = BinStrategy()->BinDiscret()->gNode(gid);
    // get the first dof gid of a particle and convert it into a LID
    int doflid = dofrowmap->LID(BinStrategy()->BinDiscret()->Dof(currparticle, 0));

    // insert values into state vectors
    (*sign)[inode] = particle_sign[gid];
    for (int idim = 0; idim < 3; ++idim) (*disnp)[doflid + idim] = (particle_positions[gid])[idim];
  }

  // reset to zero just to be sure
  Teuchos::RCP<Epetra_Vector> radius = particles_->WriteAccessRadiusn();
  radius->Scale(0.0);

  // -------------------------------------------------------------------
  //               attraction step
  // -------------------------------------------------------------------
  // move particle to the correct side of the interface

  // prepare required vectors
  // noderowmap-based vectors
  // predictor
  Teuchos::RCP<Epetra_Vector> lambda = Teuchos::rcp(new Epetra_Vector(*noderowmap, true));
  // initialize with 1.0
  lambda->PutScalar(1.0);
  // target phi_goal
  Teuchos::RCP<Epetra_Vector> phi_target = Teuchos::rcp(new Epetra_Vector(*noderowmap, true));
  // current phi of particle
  Teuchos::RCP<Epetra_Vector> phi_particle = Teuchos::rcp(new Epetra_Vector(*noderowmap, true));
  // dofrowmap-based vectors
  // normal vector of level-set field at particle position
  Teuchos::RCP<Epetra_Vector> norm_gradphi_particle =
      Teuchos::rcp(new Epetra_Vector(*dofrowmap, true));
  // increment position vector for attraction loop
  Teuchos::RCP<Epetra_Vector> inc_dis = Teuchos::rcp(new Epetra_Vector(*dofrowmap, true));
  // initial or intermediate position vector
  Teuchos::RCP<Epetra_Vector> disn = Teuchos::rcp(new Epetra_Vector(*disnp));

  // loop all particles and initialize phi_target, phi_particle and gradphi_particle
  for (int inode = 0; inode < BinStrategy()->BinDiscret()->NumMyRowNodes(); inode++)
  {
    // get sign of particle
    const double signP = (*sign)[inode];

    // get target phi for this particle
    DRT::UTILS::Random* random = DRT::Problem::Instance()->Random();
    random->SetRandRange(r_min_, b_max_);
    const double local_phi_target = signP * random->Uni();
    (*phi_target)[inode] = local_phi_target;

    // get phi and gradient of particle at current particle position
    double current_phi_particle = 0.0;
    LINALG::Matrix<3, 1> normal_particle(true);

    // get global id of current particle
    const int gid = noderowmap->GID(inode);

    // element coordinates of particle position in scatra element
    LINALG::Matrix<3, 1> elecoord(true);
    DRT::Element* scatraele = GetEleCoordinates(gid, elecoord);

    // compute level-set value
    GetLSValParticle(current_phi_particle, normal_particle, scatraele, elecoord, phinp);
    // GetLSValParticle returns gradient of phi, which we have to normalize here
    const double norm = normal_particle.Norm2();
    if (norm > 1.0e-9)
      normal_particle.Scale(1.0 / norm);
    else
      normal_particle.PutScalar(0.0);

    // store values
    (*phi_particle)[inode] = current_phi_particle;

    DRT::Node* currparticle = BinStrategy()->BinDiscret()->gNode(gid);
    // get the first dof gid of a particle and convert it into a LID
    int doflid = dofrowmap->LID(BinStrategy()->BinDiscret()->Dof(currparticle, 0));
    for (int idim = 0; idim < 3; ++idim)
      (*norm_gradphi_particle)[doflid + idim] = normal_particle(idim, 0);

    // initialize increment position vector, i.e., multiply by
    // lambda*(phi_target-current_phi_particle) note lambda is set to 1.0 at the beginning and
    // skipped here
    const double fac = local_phi_target - current_phi_particle;
    normal_particle.Scale(fac);
    for (int idim = 0; idim < 3; ++idim) (*inc_dis)[doflid + idim] = normal_particle(idim, 0);
  }

  // total (global) number of seeded particles
  const int global_num_particles = BinStrategy()->BinDiscret()->NumGlobalNodes();

  // some output and checks
  if (MyRank() == 0)
    std::cout << "--- total number of seeded particles  " << global_num_particles << std::endl;
  if (global_num_particles != (2 * num_particles_per_bin_ * allnumbins))
    dserror("Number of seeded particles does not match number of expected particles");

#if 0
  // number of placed particles should finally be equal to the number of seeded particles
  // important: this number must be global and takes into account the placed particles on
  // all procs
  int num_placed_particles = 0;
  // counter for number of attraction steps
  int step_counter = 0;
  // define vector for update of map-based vectors
  Teuchos::RCP<Epetra_Vector> old;

  // attraction loop
  while ((num_placed_particles < global_num_particles) and step_counter < 15)
  {
    // increment counter
    step_counter += 1;
    if (MyRank() == 0)
      std::cout << "--- attraction step  "<< step_counter << std::endl;

    // compute new position
    disnp = particles_->WriteAccessDispnp();
    // x_P^new = x_P + lambda*(phi_target-current_phi_particle) * n_p
    disnp->Update(1.0,*inc_dis,1.0);

    // transfer particles to bins
    TransferParticles(false);
    // update of state vectors to the new maps
    particles_->UpdateStatesAfterParticleTransfer();
    // likewise update present vectors according to the new distribution of particles
    old = lambda;
    lambda = LINALG::CreateVector(*BinStrategy()->BinDiscret()->NodeRowMap(),true);
    LINALG::Export(*old, *lambda);

    old = phi_target;
    phi_target = LINALG::CreateVector(*BinStrategy()->BinDiscret()->NodeRowMap(),true);
    LINALG::Export(*old, *phi_target);

    old = phi_particle;
    phi_particle = LINALG::CreateVector(*BinStrategy()->BinDiscret()->NodeRowMap(),true);
    LINALG::Export(*old, *phi_particle);

    old = norm_gradphi_particle;
    norm_gradphi_particle = LINALG::CreateVector(*BinStrategy()->BinDiscret()->DofRowMap(),true);
    LINALG::Export(*old, *norm_gradphi_particle);

    old = inc_dis;
    inc_dis = LINALG::CreateVector(*BinStrategy()->BinDiscret()->DofRowMap(),true);
    LINALG::Export(*old, *inc_dis);

    old = disn;
    disn = LINALG::CreateVector(*BinStrategy()->BinDiscret()->DofRowMap(),true);
    LINALG::Export(*old, *disn);

    // check new position
    // counter for local placed particles
    int local_num_placed_particles = 0;

    // get new distribution of vectors
    disnp = particles_->WriteAccessDispnp();
    sign = Teuchos::rcp_dynamic_cast<PARTICLE::TimIntRK>(particles_)->WriteAccessSign();

    // loop all particles
    for (int inode=0; inode<BinStrategy()->BinDiscret()->NumMyRowNodes(); inode++)
    {
      // get new maps
      const Epetra_Map* att_dofrowmap = BinStrategy()->BinDiscret()->DofRowMap();
      const Epetra_Map* att_noderowmap = BinStrategy()->BinDiscret()->NodeRowMap();

      // get dof lid
      const int gid = att_noderowmap->GID(inode);
      DRT::Node* currparticle = BinStrategy()->BinDiscret()->gNode(gid);
      int doflid = att_dofrowmap->LID(BinStrategy()->BinDiscret()->Dof(currparticle, 0));

      if ((*lambda)[inode] > 0.0)
      {
        // element coordinates of particle position in scatra element
        LINALG::Matrix<3,1> elecoord(true);
        DRT::Element* scatraele = GetEleCoordinates(gid,elecoord);

        if (scatraele != NULL) // particle has not left domain
        {
          // get phi and gradient of particle at current particle position
          double current_phi_particle = 0.0;
          LINALG::Matrix<3,1> normal_particle(true);
          // compute level-set value
          GetLSValParticle(current_phi_particle,normal_particle,scatraele,elecoord,phinp);
          // GetLSValParticle returns gradient of phi, which we have to normalize here
          const double norm = normal_particle.Norm2();
          if (norm > 1.0e-9)
            normal_particle.Scale(1.0/norm);
          else
            normal_particle.PutScalar(0.0);

          const double signP = (*sign)[inode];
          if ( (( signP > 0.0 and (current_phi_particle >= r_min_)
                              and (current_phi_particle <= b_max_))
              or
                ( signP < 0.0 and (current_phi_particle <= (signP * r_min_))
                              and (current_phi_particle >= (signP * b_max_))) )
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

      } // lambda != 0

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
    if (MyRank() == 0)
      std::cout << "--- placed particles  " << num_placed_particles << std::endl;

    // possibly position have been reverted and particles changed the bin or even the processor
    // therefore, we have to recall the following steps
    TransferParticles(false);
    // update of state vectors to the new maps
    particles_->UpdateStatesAfterParticleTransfer();
    // likewise update present vectors according to the new distribution of particles
    Teuchos::RCP<Epetra_Vector> old;
    old = lambda;
    lambda = LINALG::CreateVector(*BinStrategy()->BinDiscret()->NodeRowMap(),true);
    LINALG::Export(*old, *lambda);

    old = phi_target;
    phi_target = LINALG::CreateVector(*BinStrategy()->BinDiscret()->NodeRowMap(),true);
    LINALG::Export(*old, *phi_target);

    old = phi_particle;
    phi_particle = LINALG::CreateVector(*BinStrategy()->BinDiscret()->NodeRowMap(),true);
    LINALG::Export(*old, *phi_particle);

    old = norm_gradphi_particle;
    norm_gradphi_particle = LINALG::CreateVector(*BinStrategy()->BinDiscret()->DofRowMap(),true);
    LINALG::Export(*old, *norm_gradphi_particle);

    old = inc_dis;
    inc_dis = LINALG::CreateVector(*BinStrategy()->BinDiscret()->DofRowMap(),true);
    LINALG::Export(*old, *inc_dis);

    old = disn;
    disn = LINALG::CreateVector(*BinStrategy()->BinDiscret()->DofRowMap(),true);
    LINALG::Export(*old, *disn);

    // update state vectors
    disnp = particles_->WriteAccessDispnp();
    disn->Update(1.0,*disnp,0.0);
  }

  // delete particles which could not be attracted correctly -> all particles with lambda != 0.0
  if (num_placed_particles != global_num_particles) // we have particles with lambda != 0.0
  {
    if (MyRank() == 0)
      std::cout << "--- delete unplaced particles  "<< global_num_particles-num_placed_particles << std::endl;

    // vector with particles to delete
    std::vector<int> part_del;

    // loop all particles
    for (int inode=0; inode<BinStrategy()->BinDiscret()->NumMyRowNodes(); inode++)
    {
      // particles with lambda != 0.0 will be deleted
      if ((*lambda)[inode] > (0.0+1.0e-8))
      {
        // get particle
        DRT::Node *currparticle = BinStrategy()->BinDiscret()->lRowNode(inode);
        // store gid
        part_del.push_back(currparticle->Id());
      }
    }

    // remove particles
    DeleteParticles(part_del);

    // transfer phi_particle, which will be used below to set the radius, to new map
    old = phi_particle;
    phi_particle = LINALG::CreateVector(*BinStrategy()->BinDiscret()->NodeRowMap(),true);
    LINALG::Export(*old, *phi_particle);
  }

#else
  // do attraction loop
  Attraction(
      global_num_particles, inc_dis, lambda, phi_target, phi_particle, norm_gradphi_particle, disn);
  // caution: this function changes pointers
#endif

    // -------------------------------------------------------------------
    //               set radius
    // -------------------------------------------------------------------

#if 0
  // get radius and sign vector
  radius = particles_->WriteAccessRadius();
  sign = Teuchos::rcp_dynamic_cast<PARTICLE::TimIntRK>(particles_)->WriteAccessSign();

  for (int inode=0; inode<BinStrategy()->BinDiscret()->NumMyRowNodes(); inode++)
  {
    // get phi and sign
    const double myphi = (*phi_particle)[inode];
    const double mysign = (*sign)[inode];

    // set radius
    if ((mysign*myphi) > r_max_)
      (*radius)[inode] = r_max_;
    else if ((mysign*myphi) < r_min_)
      (*radius)[inode] = r_min_;
    else
      (*radius)[inode] = (mysign*myphi);
  }
#else
  AdjustParticleRadii();
#endif

  if (MyRank() == 0) std::cout << "----------------------------- " << std::endl;

  // testing of convergence of rk: do not delete this, it might be helpful
  //  // get position vector
  //  disnp = particles_->WriteAccessDispnp();
  //  // assume at least one particle and set desired position
  //  (*disnp)[0]=-6.0;
  //  (*disnp)[1]=0.0;
  //  disn->Update(1.0,*disnp,0.0);
  //  // adapt state vectors to new maps
  //  TransferParticles(false);
  //  particles_->UpdateStatesAfterParticleTransfer();
  //  // the print current position of first particle in rk scheme
  //  // comment correction and reseeding step in level-set algorithm

  return;
}


/*----------------------------------------------------------------------*
 | solve the current particle time step                 rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::ScatraParticleCoupling::Integrate()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::Algorithm::Integrate");
  particles_->IntegrateStep();

  // in case of holes in the domain, inaccuracies in the particle approach
  // may transport particles into this hole
  // these particles are lost and have to be deleted
  DeleteParticlesOutOfPhysicalDomain();

  return;
}


/*----------------------------------------------------------------------*
 | solve the current particle time step                 rasthofer 02/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> PARTICLE::ScatraParticleCoupling::CorrectionStep()
{
  // -------------------------------------------------------------------
  //                  setup correction by particles
  // -------------------------------------------------------------------

  // get level-set values on nodes
  const Teuchos::RCP<const Epetra_Vector> row_phinp = scatra_->Phinp();
  // export phi vector to col map
  // this is necessary here, since the elements adjacent to a row node may have
  // ghosted nodes, which are only contained in the col format
  const Epetra_Map* scatra_dofcolmap = scatradis_->DofColMap();
  Teuchos::RCP<Epetra_Vector> phinp = Teuchos::rcp(new Epetra_Vector(*scatra_dofcolmap));
  LINALG::Export(*row_phinp, *phinp);

  // get sign vector
  Teuchos::RCP<const Epetra_Vector> sign =
      Teuchos::rcp_dynamic_cast<PARTICLE::TimIntRK>(particles_)->Sign();
  // and radius vector
  Teuchos::RCP<const Epetra_Vector> radius =
      Teuchos::rcp_dynamic_cast<PARTICLE::TimIntRK>(particles_)->Radiusn();

  // -------------------------------------------------------------------
  //                     find escaped particles
  // -------------------------------------------------------------------

  if (MyRank() == 0) std::cout << "-------- PARTICLE CORRECTION --------- " << std::endl;

  // prepare map for elements and corresponding escaped particles (particles on this proc only)
  // int: global id element
  // set<int>: global id of escaped particles for this element
  std::map<int, std::set<int>> escaped_particles_list;

  // loop all particles
  for (int inode = 0; inode < BinStrategy()->BinDiscret()->NumMyRowNodes(); inode++)
  {
    // get phi and gradient of particle at current particle position
    double phi_particle = 0.0;
    // we do not need the gradient here, so this is just a dummy
    LINALG::Matrix<3, 1> normal_particle(true);

    // element coordinates of particle position in scatra element
    LINALG::Matrix<3, 1> elecoord(true);
    DRT::Element* scatraele =
        GetEleCoordinates(BinStrategy()->BinDiscret()->NodeRowMap()->GID(inode), elecoord);

    // compute level-set value (false=no computation of gradient)
    GetLSValParticle(phi_particle, normal_particle, scatraele, elecoord, phinp, false);

    // check, whether particle is escaped or not
    // criterion
    // sign(P) * phi(x_p) < 0.0, -> half
    // or this particle is on wrong side by more than its radius -> full

    // get sign
    const double signP = (*sign)[inode];
    if (escaped_ == INPAR::PARTICLEOLD::half)
    {
      if (signP * phi_particle < 0.0)
        // add escaped particle
        escaped_particles_list[scatraele->Id()].insert(
            (BinStrategy()->BinDiscret()->lRowNode(inode))->Id());
    }
    else if (escaped_ == INPAR::PARTICLEOLD::full)
    {
      // get radius
      const double radiusP = (*radius)[inode];
      if ((signP * phi_particle < 0.0) and
          (std::abs(phi_particle) >
              radiusP))  // phi_particle equal to distance to interface (at least approximately)
        // add escaped particle
        escaped_particles_list[scatraele->Id()].insert(
            (BinStrategy()->BinDiscret()->lRowNode(inode))->Id());
    }
    else
      dserror("Unknown escape criterion!");

  }  // end loop all particles

  //  if (MyRank()==1)
  //  {
  //  for (std::map<int,std::set<int> >::iterator it=escaped_particles_list.begin();
  //  it!=escaped_particles_list.end(); it++)
  //     for (std::set<int>::iterator ite=it->second.begin(); ite!=it->second.end(); ite++)
  //         std::cout << "ele  " << it->first << "  particles  " << *ite << std::endl;
  //  }

  // -------------------------------------------------------------------
  //                     communicate escaped particles
  // -------------------------------------------------------------------

  // complete map of elements and corresponding escaped particles (including ghosted scatra
  // elements)
  std::map<int, std::set<int>> extended_escaped_particles_list;
  // fill map and setup particle ghosting (provides column map)
  BinStrategy()->ExtendGhosting(
      scatradis_, escaped_particles_list, extended_escaped_particles_list);
  // since all row/col maps and dofs have been killed, we have to rebuild state vectors
  particles_->UpdateStatesAfterParticleTransfer();

  //  if (MyRank()==1)
  //  {
  //  for (std::map<int,std::set<int> >::iterator it=extended_escaped_particles_list.begin();
  //  it!=extended_escaped_particles_list.end(); it++)
  //     for (std::set<int>::iterator ite=it->second.begin(); ite!=it->second.end(); ite++)
  //       std::cout << "ele  " << it->first << "  particles  " << *ite << std::endl;
  //  }

  // -------------------------------------------------------------------
  //                         perform correction
  // -------------------------------------------------------------------

  // vector for corrected values
  Teuchos::RCP<Epetra_Vector> corrected_phinp =
      Teuchos::rcp(new Epetra_Vector(*(scatradis_->DofRowMap()), true));

  // get sign, radius and position vector in col layout
  Teuchos::RCP<Epetra_Vector> signcol =
      Teuchos::rcp(new Epetra_Vector(*(BinStrategy()->BinDiscret()->NodeColMap()), true));
  LINALG::Export(*(Teuchos::rcp_dynamic_cast<PARTICLE::TimIntRK>(particles_)->Sign()), *signcol);

  Teuchos::RCP<Epetra_Vector> radiuscol =
      Teuchos::rcp(new Epetra_Vector(*(BinStrategy()->BinDiscret()->NodeColMap()), true));
  LINALG::Export(*(particles_->Radiusn()), *radiuscol);

  Teuchos::RCP<Epetra_Vector> disnpcol =
      Teuchos::rcp(new Epetra_Vector(*(BinStrategy()->BinDiscret()->DofColMap()), true));
  LINALG::Export(*(particles_->Dispnp()), *disnpcol);

#if 0
  // loop all nodes
  for (int inode=0; inode<scatradis_->NumMyRowNodes(); inode++)
  {
    // get node
    DRT::Node* actnode = scatradis_->lRowNode(inode);
    // number of adjacent elements
    const int numele = actnode->NumElement();
    // get list of adjacent elements
    const DRT::Element*const* adjele = actnode->Elements();

    // get coordinates of node
    LINALG::Matrix<3,1> nodecoords(true);
    for (int idim=0; idim<3; idim++)
      nodecoords(idim,0) = (actnode->X())[idim];

    // get dof of node
    std::vector<int> dofs = scatradis_->Dof(0,actnode);
    if (dofs.size() != 1)
      dserror("Only one dof expected for level-set node!");
    const int doflid = scatradis_->DofRowMap()->LID(dofs[0]);
    // initialize nodal values for correction
    double phi_node_plus = (*row_phinp)[doflid];
    double phi_node_minus =(*row_phinp)[doflid];
    double corrected_phi = 0.0;

    // loop all adjacent elements
    for(int iele=0; iele<numele; iele++)
    {
      // check list of elements with escaped particles
      std::map<int,std::set<int> >::iterator iter = extended_escaped_particles_list.find(adjele[iele]->Id());

      // do we have escaped particles for this element
      if (iter != extended_escaped_particles_list.end())
      {
//        std::cout << "node " << actnode->Id() << " is corrected" << std::endl;
        // get vector of escaped particles
        std::set<int> particleset = iter->second;

        // loop all escaped particles
        for (std::set<int>::iterator ipart=particleset.begin(); ipart!=particleset.end(); ipart++)
        {
          // get current particle id
          const int partgid = *ipart;
          // transfer to lid
          const int partlid = BinStrategy()->BinDiscret()->NodeColMap()->LID(partgid);
          // get sign of particle
          const double actsignP = (*signcol)[partlid];
          // get radius of particle
          const double actradiusP =(*radiuscol)[partlid];

          // get position of particle
          LINALG::Matrix<3,1> actpositionP(true);
          DRT::Node* currparticle = BinStrategy()->BinDiscret()->gNode(partgid);
          // get the first dof gid of a particle and convert it into a LID
          const int partdoflid = BinStrategy()->BinDiscret()->DofColMap()->LID(BinStrategy()->BinDiscret()->Dof(currparticle, 0));
          for (int idim=0; idim<3; ++idim)
            actpositionP(idim,0) = (*disnpcol)[partdoflid+idim];

          //std::cout << " particle with id   " << partgid << " has sign " << actsignP << "  radius  " << actradiusP << " and position  " << actpositionP << std::endl;

          // temporary vector for prediction level-set value by particle
          LINALG::Matrix<3,1> tmp(true);
          tmp.Update(1.0,actpositionP,0.0);
          tmp.Update(1.0,nodecoords,-1.0);

          // adapt distance vector in case of quasi-2D simulations
          switch (particle_dim_)
          {
            case INPAR::PARTICLEOLD::particle_2Dx:
            {
              tmp(0,0) = 0.0;
              break;
            }
            case INPAR::PARTICLEOLD::particle_2Dy:
            {
              tmp(1,0) = 0.0;
              break;
            }
            case INPAR::PARTICLEOLD::particle_2Dz:
            {
              tmp(2,0) = 0.0;
              break;
            }
            default:
              break;
          }

          // predict level-set value by particle

          const double norm_tmp = tmp.Norm2();
          double predicted_phi = actsignP * (actradiusP - norm_tmp);

          // update phi_node_plus and phi_node_minus
          if (actsignP > 0.0)
            phi_node_plus = std::max(phi_node_plus,predicted_phi);
          else
            phi_node_minus = std::min(phi_node_minus,predicted_phi);
        } // end loop escaped particles

//        std::cout << "plus  " << phi_node_plus << "  minus  " << phi_node_minus << std::endl;
      }
    } // end loop all adjacent elements

    // compute corrected phi
    if (std::abs(phi_node_plus) <= std::abs(phi_node_minus))
      corrected_phi = phi_node_plus;
    else
      corrected_phi = phi_node_minus;

    // store new value
    int err = corrected_phinp->ReplaceMyValues(1,&corrected_phi,&inode);
    if (err != 0) dserror("Could not store corrected value");

  }// end loop all nodes
#endif

  // initialize corrected phi values with values of level-set solution
  corrected_phinp->Update(1.0, *row_phinp, 0.0);

  // loop all elements with escaped particles
  for (std::map<int, std::set<int>>::iterator iter = extended_escaped_particles_list.begin();
       iter != extended_escaped_particles_list.end(); iter++)
  {
    // get element from map
    DRT::Element* actele = scatradis_->gElement(iter->first);

    // loop all nodes of element
    for (int inode = 0; inode < actele->NumNode(); inode++)
    {
      // get node
      DRT::Node* actnode = (actele->Nodes())[inode];

      // we only correct row nodes
      if (actnode->Owner() == MyRank())
      {
        // get coordinates of node
        LINALG::Matrix<3, 1> nodecoords(true);
        for (int idim = 0; idim < 3; idim++) nodecoords(idim, 0) = (actnode->X())[idim];

        // get dof of node
        std::vector<int> dofs = scatradis_->Dof(0, actnode);
        if (dofs.size() != 1) dserror("Only one dof expected for level-set node!");
        const int doflid = scatradis_->DofRowMap()->LID(dofs[0]);
        // initialize nodal values for correction
        double phi_node_plus = (*corrected_phinp)[doflid];
        double phi_node_minus = (*corrected_phinp)[doflid];
        double corrected_phi = 0.0;

        // get set of escaped particles
        std::set<int> particleset = iter->second;

        // loop all escaped particles
        for (std::set<int>::iterator ipart = particleset.begin(); ipart != particleset.end();
             ipart++)
        {
          // get current particle id
          const int partgid = *ipart;
          // transfer to lid
          const int partlid = BinStrategy()->BinDiscret()->NodeColMap()->LID(partgid);
          // get sign of particle
          const double actsignP = (*signcol)[partlid];
          // get radius of particle
          const double actradiusP = (*radiuscol)[partlid];

          // get position of particle
          LINALG::Matrix<3, 1> actpositionP(true);
          DRT::Node* currparticle = BinStrategy()->BinDiscret()->gNode(partgid);
          // get the first dof gid of a particle and convert it into a LID
          const int partdoflid = BinStrategy()->BinDiscret()->DofColMap()->LID(
              BinStrategy()->BinDiscret()->Dof(currparticle, 0));
          for (int idim = 0; idim < 3; ++idim)
            actpositionP(idim, 0) = (*disnpcol)[partdoflid + idim];

          // temporary vector for prediction level-set value by particle
          LINALG::Matrix<3, 1> tmp(true);
          tmp.Update(1.0, actpositionP, 0.0);
          tmp.Update(1.0, nodecoords, -1.0);

          // adapt distance vector in case of quasi-2D simulations
          switch (BinStrategy()->ParticleDim())
          {
            case INPAR::PARTICLEOLD::particle_2Dx:
            {
              tmp(0, 0) = 0.0;
              break;
            }
            case INPAR::PARTICLEOLD::particle_2Dy:
            {
              tmp(1, 0) = 0.0;
              break;
            }
            case INPAR::PARTICLEOLD::particle_2Dz:
            {
              tmp(2, 0) = 0.0;
              break;
            }
            default:
              break;
          }

          // predict level-set value by particle
          const double norm_tmp = tmp.Norm2();
          double predicted_phi = actsignP * (actradiusP - norm_tmp);

          // update phi_node_plus and phi_node_minus
          if (actsignP > 0.0)
            phi_node_plus = std::max(phi_node_plus, predicted_phi);
          else
            phi_node_minus = std::min(phi_node_minus, predicted_phi);

        }  // end loop escaped particles

        // compute corrected phi
        if (std::abs(phi_node_plus) <= std::abs(phi_node_minus))
          corrected_phi = phi_node_plus;
        else
          corrected_phi = phi_node_minus;

        // store new value
        int err = corrected_phinp->ReplaceMyValues(1, &corrected_phi, &doflid);
        if (err != 0) dserror("Could not store corrected value");
      }
    }  // loop all nodes of element

  }  // end loop all elements with escaped particles


  return corrected_phinp;
}


/*----------------------------------------------------------------------*
 | adjust particle radii to current interface position  rasthofer 02/14 |
 *----------------------------------------------------------------------*/
void PARTICLE::ScatraParticleCoupling::AdjustParticleRadii()
{
  // get level-set values on nodes
  const Teuchos::RCP<const Epetra_Vector> row_phinp = scatra_->Phinp();
  // export phi vector to col map
  // this is necessary here, since the elements adjacent to a row node may have
  // ghosted nodes, which are only contained in the col format
  const Epetra_Map* scatra_dofcolmap = scatradis_->DofColMap();
  Teuchos::RCP<Epetra_Vector> phinp = Teuchos::rcp(new Epetra_Vector(*scatra_dofcolmap));
  LINALG::Export(*row_phinp, *phinp);

  // get radius and sign vector
  Teuchos::RCP<Epetra_Vector> radius = particles_->WriteAccessRadiusn();
  Teuchos::RCP<const Epetra_Vector> sign =
      Teuchos::rcp_dynamic_cast<PARTICLE::TimIntRK>(particles_)->Sign();

  // loop all particles
  for (int inode = 0; inode < BinStrategy()->BinDiscret()->NumMyRowNodes(); inode++)
  {
    // get phi of particle at current particle position
    double myphi = 0.0;
    // we do not need the gradient here, so this is just a dummy
    LINALG::Matrix<3, 1> normal_particle(true);

    // element coordinates of particle position in scatra element
    LINALG::Matrix<3, 1> elecoord(true);
    DRT::Element* scatraele =
        GetEleCoordinates(BinStrategy()->BinDiscret()->NodeRowMap()->GID(inode), elecoord);

    // compute level-set value (false=no computation of gradient)
    GetLSValParticle(myphi, normal_particle, scatraele, elecoord, phinp, false);

    // get sign
    const double mysign = (*sign)[inode];

    if ((mysign * myphi) > r_max_)
      (*radius)[inode] = r_max_;
    else if ((mysign * myphi) < r_min_)
      (*radius)[inode] = r_min_;
    else
      (*radius)[inode] = (mysign * myphi);
  }  // end loop all particles

  return;
}


/*----------------------------------------------------------------------*
 | transfer particles into their correct bins and update solution       |
 |                                                      rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::ScatraParticleCoupling::Update()
{
  // adjust particle radii to current interface position
  // the order AdjustParticleRadii() and then Reseeding() is important,
  // since radii are used in the sorting algorithms
  AdjustParticleRadii();

  // we take step from scatra instead of particles, since particles use structure
  // style, where time and step are updated at the end of the time step (see below),
  // whereas scatra uses fluid style, where time and step are set at the beginning
  const int step = scatra_->Step();
  if ((step % reseeding_) == 0) Reseeding();

  // update displacements, velocities, accelerations
  // after this call we will have disn_==dis_, etc
  // update time and step
  // for Runge-Kutta schemes only time and step are updated
  particles_->Update();

  return;
}


/*----------------------------------------------------------------------*
 | reseeding of particles                               rasthofer 01/14 |
 *----------------------------------------------------------------------*/
void PARTICLE::ScatraParticleCoupling::Reseeding()
{
  // if the level-set field is not reinitialized periodically, a
  // forced reinitialization is required at this place
  // otherwise, the phi-criterions concerning the interface band will fail

  // -------------------------------------------------------------------
  //               prepare reseeding
  // -------------------------------------------------------------------

  if (MyRank() == 0) std::cout << "-------- RESEEDING --------- " << std::endl;

  // -------------------------------------------------------------------
  //               determine useless particles to be deleted
  //                  determine new particles to be seeded
  // -------------------------------------------------------------------

  // get node row map
  const Epetra_Map* noderowmap = BinStrategy()->BinDiscret()->NodeRowMap();
  // get sign and radius vector
  Teuchos::RCP<Epetra_Vector> sign =
      Teuchos::rcp_dynamic_cast<PARTICLE::TimIntRK>(particles_)->WriteAccessSign();
  Teuchos::RCP<Epetra_Vector> radius = particles_->WriteAccessRadiusn();

  // map containing ids of bins (equal to elements) which will get additional particles
  // as well as the required number of positive and negative particles
  // pair<plus,minus> particles
  std::map<int, std::pair<int, int>> particlebins;
  // vector containing ids of particles to be deleted
  std::vector<int> deleteparticles;

  // get phinp
  const Teuchos::RCP<const Epetra_Vector> row_phinp = scatra_->Phinp();
  // export phi vector to col map
  // this is necessary here, since the elements adjacent to a row node may have
  // ghosted nodes, which are only contained in the col format
  const Epetra_Map* scatra_dofcolmap = scatradis_->DofColMap();
  Teuchos::RCP<Epetra_Vector> phinp = Teuchos::rcp(new Epetra_Vector(*scatra_dofcolmap));
  LINALG::Export(*row_phinp, *phinp);

  // loop all row bins on this proc
  for (int ibin = 0; ibin < BinStrategy()->BinDiscret()->NumMyRowElements(); ibin++)
  {
    // get pointer to current bin
    DRT::Element* actele = BinStrategy()->BinDiscret()->lRowElement(ibin);
    DRT::MESHFREE::MeshfreeMultiBin* actbin =
        dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(actele);

    // get pointer to associated scatra elements
    DRT::Element** scatraelesinbin = actbin->AssociatedEles(bin_transpcontent_);
    // check for null-pointer in case of holes in the domain
    if (scatraelesinbin == NULL) continue;

#ifndef MULTIPLE_ELE
    // do a first fast check, if bin is near interface
    {
      DRT::Element* first_scatraele = scatraelesinbin[0];
      // get first node of first scatra element
      const DRT::Node* first_scatranode = (first_scatraele->Nodes())[0];
      const int scatra_dofgid = scatradis_->Dof(0, first_scatranode, 0);
      const int lid = scatra_dofcolmap->LID(scatra_dofgid);
      const double first_phi = (*phinp)[lid];

      // if bin is far away from the interface,
      // all non-escaped particles will be deleted
      if (std::abs(first_phi) > (10.0 * b_max_))
      {
        // loop particles if available
        for (int ipart = 0; ipart < actbin->NumNode(); ipart++)
        {
          // get current particle
          const DRT::Node* currentparticle = (actbin->Nodes())[ipart];
          // get gid and lid
          const int gid = currentparticle->Id();
          const int lid = noderowmap->LID(gid);
          // check for non-escaped particles only
          if (((*sign)[lid] * first_phi) > 0.0) deleteparticles.push_back(gid);
        }

        // everything is done for this bin: go to the next one
        continue;
      }
    }

    // everything that ends up here is near the interface, and we have to have a
    // closer look at these bins

    // get bin center
    LINALG::Matrix<3, 1> center(true);
    center = GetBinCentroid(actbin->Id());

    // variables to store information about the element with center in bin
    DRT::Element* targetscatraele = NULL;
    bool found_ele = false;

    // search for underlying scatra element
    for (int iele = 0; iele < actbin->NumAssociatedEle(bin_transpcontent_); ++iele)
    {
      DRT::Element* scatraele = scatraelesinbin[iele];
      if (scatraele->Shape() != DRT::Element::hex8)
        dserror("Other element than hex8 not yet supported!");

      // get element center
      const DRT::Element::DiscretizationType distype = DRT::Element::hex8;
      const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      // get center coordinates of current element
      LINALG::Matrix<3, 1> elecenter(true);
      LINALG::Matrix<3, numnode> xyzele(true);
      GEO::fillInitialPositionArray<distype>(scatraele, xyzele);

      for (int inode = 0; inode < numnode; inode++)
        for (int idim = 0; idim < 3; idim++) elecenter(idim, 0) += xyzele(idim, inode);

      const double inv_numnode = 1.0 / ((double)numnode);
      for (int idim = 0; idim < 3; idim++) elecenter(idim, 0) *= inv_numnode;

      if ((elecenter(0, 0) > (center(0, 0) - (bin_size_[0] / 2.0))) and
          (elecenter(0, 0) < (center(0, 0) + (bin_size_[0] / 2.0))))
      {
        if ((elecenter(1, 0) > (center(1, 0) - (bin_size_[1] / 2.0))) and
            (elecenter(1, 0) < (center(1, 0) + (bin_size_[1] / 2.0))))
        {
          if ((elecenter(2, 0) > (center(2, 0) - (bin_size_[2] / 2.0))) and
              (elecenter(2, 0) < (center(2, 0) + (bin_size_[2] / 2.0))))
          {
            targetscatraele = scatraele;
            if (found_ele == false)
              found_ele = true;
            else
              dserror("More than one element found for this bin!");
          }
        }
      }
    }

    if (targetscatraele == NULL) dserror("Could not found corresponding element!");
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;

    // get nodal phi values
    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    lm.clear();
    lmowner.clear();
    lmstride.clear();
    targetscatraele->LocationVector(*scatradis_, lm, lmowner, lmstride);
    std::vector<double> myphinp(lm.size());
    DRT::UTILS::ExtractMyValues(*phinp, myphinp, lm);
#else

    // do a first fast check, if bin is near interface
    if (fast_)
    {
      DRT::Element* scatraele = scatraelesinbin[0];
      // get first node of first scatra element
      const DRT::Node* first_scatranode = (scatraele->Nodes())[0];
      const int scatra_dofgid = scatradis_->Dof(0, first_scatranode, 0);
      const int lid = scatra_dofcolmap->LID(scatra_dofgid);
      const double first_phi = (*phinp)[lid];

      // if bin is far away from the interface,
      // all non-escaped particles will be deleted
      if (std::abs(first_phi) > (10.0 * b_max_))
      {
        // loop particles if available
        for (int ipart = 0; ipart < actbin->NumNode(); ipart++)
        {
          // get current particle
          const DRT::Node* currentparticle = (actbin->Nodes())[ipart];
          // get gid and lid
          const int gid = currentparticle->Id();
          const int lid = noderowmap->LID(gid);
          // check for non-escaped particles only
          if (((*sign)[lid] * first_phi) > 0.0) deleteparticles.push_back(gid);
        }

        // everything is done for this bin: go to the next one
        continue;
      }
    }

    // everything that ends up here is near the interface, and we have to have a
    // closer look at these bins

    // get bin corners
    std::vector<LINALG::Matrix<3, 1>> actbincorner;
    BinStrategy()->GetBinCorners(actbin->Id(), actbincorner);

    // vector of corner values
    std::vector<double> myphinp;

    for (std::size_t icor = 0; icor < actbincorner.size(); icor++)
    {
      // get element coordinates of bin corner and associated scatra element
      LINALG::Matrix<3, 1> elecoord(true);
      DRT::Element* scatraele = GetEleCoordinatesFromPosition(actbincorner[icor], actbin, elecoord);
      // this is a simple solution for holes in the domain
      if (scatraele == NULL) continue;

      // compute level-set value
      // get phi and gradient of particle at current particle position
      double phi_corner = 0.0;
      // we do not need the gradient here, so this is just a dummy
      LINALG::Matrix<3, 1> normal_particle(true);
      GetLSValParticle(phi_corner, normal_particle, scatraele, elecoord, phinp, false);

      // store value
      myphinp.push_back(phi_corner);
    }
    // this is simple solution for holes in the domain
    if (myphinp.size() < 8) continue;

    // a bin always has eight corners!
    const int numnode = 8;
#endif

    // determine number of nodes in band and intersection status
    int numnodes_in_band = 0;
    bool intersected = false;
    double phi_ref = 0.0;
    bool set_phi_ref = false;

    for (std::size_t rr = 0; rr < myphinp.size(); rr++)
    {
      // is node in band
      if (std::abs(myphinp[rr]) <= b_max_) numnodes_in_band += 1;

      if (not set_phi_ref)
      {
        // set phi_ref to first non-zero value
        if (std::abs(myphinp[rr]) > 1.0e-9)
        {
          phi_ref = myphinp[rr];
          set_phi_ref = true;
        }
      }
      else
      {
        // change in sign
        if ((myphinp[rr] * phi_ref) < 0.0)  // this criterion excludes all touched cases
          intersected = true;
      }
    }

    // bin is not in band
    if (numnodes_in_band == 0 and intersected == false)
    {
      // loop particles if available
      for (int ipart = 0; ipart < actbin->NumNode(); ipart++)
      {
        // get current particle
        const DRT::Node* currentparticle = (actbin->Nodes())[ipart];
        // get gid and lid
        const int gid = currentparticle->Id();
        const int lid = noderowmap->LID(gid);
        // check for non-escaped particles only
        if (((*sign)[lid] * phi_ref) > 0.0) deleteparticles.push_back(gid);
      }
    }
    else if (numnodes_in_band == 0 and intersected == true)
    {
      dserror("Too large bin!");
    }
    // bin is intersected
    else if (numnodes_in_band == numnode and intersected == true)
    {
#ifndef MULTIPLE_ELE
      // first, get phi at cell center (assumed equal to element center) (Hex8 assumed!!!!!)
      LINALG::Matrix<3, 1> centerelecoord(true);
      // fill all element arrays
      LINALG::Matrix<numnode, 1> ephinp(true);
      for (int i = 0; i < numnode; ++i) ephinp(i, 0) = myphinp[i];

      // shape functions
      LINALG::Matrix<numnode, 1> funct(true);
      // fill vectors
      DRT::UTILS::shape_function_3D(
          funct, centerelecoord(0), centerelecoord(1), centerelecoord(2), targetscatraele->Shape());

      // finally compute phi
      double phi = funct.Dot(ephinp);
#else

      // first, get phi at center
      // get bin center
      LINALG::Matrix<3, 1> center(true);
      center = BinStrategy()->GetBinCentroid(actbin->Id());

      // get element coordinates of bin center and associated scatra element
      LINALG::Matrix<3, 1> elecoord(true);
      DRT::Element* scatraele = GetEleCoordinatesFromPosition(center, actbin, elecoord);

      // compute level-set value
      // get phi and gradient of particle at current particle position
      double phi = 0.0;
      // we do not need the gradient here, so this is just a dummy
      LINALG::Matrix<3, 1> normal_particle(true);
      GetLSValParticle(phi, normal_particle, scatraele, elecoord, phinp, false);

#endif

      // compute portion of minus domain

      // compute largest distance in bin
      double h_ref = 0.0;
      if (BinStrategy()->ParticleDim() == INPAR::PARTICLEOLD::particle_3D)
      {
        // characteristic bin length
        const double h = std::pow(
            BinStrategy()->BinSize()[0] * BinStrategy()->BinSize()[1] * BinStrategy()->BinSize()[2],
            1.0 / 3.0);
        h_ref = std::sqrt(3.0) * h;
      }
      else
      {
        // characteristic bin length
        double h = -1.0;
        if (BinStrategy()->ParticleDim() == INPAR::PARTICLEOLD::particle_2Dx)
          h = std::sqrt(BinStrategy()->BinSize()[1] * BinStrategy()->BinSize()[2]);
        if (BinStrategy()->ParticleDim() == INPAR::PARTICLEOLD::particle_2Dy)
          h = std::sqrt(BinStrategy()->BinSize()[0] * BinStrategy()->BinSize()[2]);
        if (BinStrategy()->ParticleDim() == INPAR::PARTICLEOLD::particle_2Dz)
          h = std::sqrt(BinStrategy()->BinSize()[0] * BinStrategy()->BinSize()[1]);
        h_ref = std::sqrt(2.0) * h;
      }

      double vol_center = 1.0 - (std::abs(phi) / h_ref + 0.5);
      // safety check
      // if (vol_center<0.0) dserror("Volume should be positive!");
      if (vol_center < 0.0)
      {
        // be tolerant
        if (std::abs(vol_center) < 0.75)
          vol_center = 0.0;
        else
        {
          std::cout << "vol " << vol_center << std::endl;
          for (std::size_t rr = 0; rr < myphinp.size(); rr++) std::cout << myphinp[rr] << std::endl;
          dserror("Volume should be positive!");
        }
      }

      // number of desired negative and positive particles
      int num_particles_minus = 0.0;
      int num_particles_plus = 0.0;

      if (phi >= 0.0)
      {
        num_particles_minus = (int)(((double)(2 * num_particles_per_bin_)) * vol_center);
        num_particles_plus = (2 * num_particles_per_bin_) - num_particles_minus;
      }
      else
      {
        num_particles_plus = (int)(((double)(2 * num_particles_per_bin_)) * vol_center);
        num_particles_minus = (2 * num_particles_per_bin_) - num_particles_plus;
      }

      int current_number_particles_minus = 0;
      int current_number_particles_plus = 0;
      // loop particles of intersected bin
      for (int ipart = 0; ipart < actbin->NumNode(); ipart++)
      {
        // get current particle
        const DRT::Node* currentparticle = (actbin->Nodes())[ipart];
        // get gid and lid
        const int gid = currentparticle->Id();
        const int lid = noderowmap->LID(gid);
        // get sign
        const double signP = (*sign)[lid];

        // determine number of positive and negative particles
        if (signP > 0.0)
          current_number_particles_plus += 1;
        else
          current_number_particles_minus += 1;
      }

      // add bin to list if new particles have to be inserted
      if ((current_number_particles_plus < num_particles_plus) or
          (current_number_particles_minus < num_particles_minus))
      {
        std::pair<int, int> seed_particles(0, 0);
        if (current_number_particles_plus < num_particles_plus)
          seed_particles.first = num_particles_plus - current_number_particles_plus;
        if (current_number_particles_minus < num_particles_minus)
          seed_particles.second = num_particles_minus - current_number_particles_minus;

        std::pair<int, std::pair<int, int>> fillbin(actbin->Id(), seed_particles);
        particlebins.insert(fillbin);
      }
      else
      {
        // delete particles if there are too many positive or negative ones
        // list for particles sorted by distance from interface
        // pair< int eleGID, double sortcriterion >
        std::list<std::pair<int, double>> particle_list_plus;
        std::list<std::pair<int, double>> particle_list_minus;

        for (int ipart = 0; ipart < actbin->NumNode(); ipart++)
        {
          // get current particle
          const DRT::Node* currentparticle = (actbin->Nodes())[ipart];
          // get gid and lid
          const int gid = currentparticle->Id();
          const int lid = noderowmap->LID(gid);
          // get sign
          const double signP = (*sign)[lid];
          // and radius
          const double radiusP = (*radius)[lid];

          // get phi and gradient of particle at current particle position
          double current_phi_particle = 0.0;
          // dummy vector
          LINALG::Matrix<3, 1> normal_particle(true);

          // element coordinates of particle position in scatra element
          LINALG::Matrix<3, 1> elecoord(true);
          DRT::Element* scatraele = GetEleCoordinates(gid, elecoord);

          // compute level-set value
          GetLSValParticle(
              current_phi_particle, normal_particle, scatraele, elecoord, phinp, false);

          // insert particle
          std::pair<int, double> thispair;
          thispair.first = gid;
          // compute sort criterion
          thispair.second = signP * current_phi_particle - radiusP;

          if (signP > 0.0)
          {
            // std::cout << "plus particle " << current_phi_particle << " sort " <<
            // signP*current_phi_particle-radiusP << std::endl;
            particle_list_plus.push_back(thispair);
          }
          else
          {
            // std::cout << "minus particle " << current_phi_particle << " sort " <<
            // signP*current_phi_particle-radiusP << std::endl;
            particle_list_minus.push_back(thispair);
          }
        }

        // sort particles in descending order
        // for escaped particles, which are not excluded, sort criterion yields negative value,
        // and they are stored at the bottom
        // this is the STL sorting, which is pretty fast
        particle_list_minus.sort(MyComparePairs);
        particle_list_plus.sort(MyComparePairs);

        // number of particles to be deleted
        const int numdelpartplus = current_number_particles_plus - num_particles_plus;
        std::list<std::pair<int, double>>::iterator iterplus = particle_list_plus.begin();
        // std::cout << " numdelpartplus " << numdelpartplus << " current_number_particles_plus " <<
        // current_number_particles_plus << " num_particles_plus " << num_particles_plus <<
        // std::endl;
        // take first values of list
        for (int ipart = 0; ipart < numdelpartplus; ipart++)
        {
          if (iterplus->second >=
              0.0)  // also delete particles with r_min, as they can be reconstructed
          {
            deleteparticles.push_back(iterplus->first);
            // std::cout << " plus " << iterplus->first << std::endl;
            iterplus++;
          }
          else
            break;  // all the remaining particles are either escaped or are closer to the interface
                    // than r_min
        }

        const int numdelpartminus = current_number_particles_minus - num_particles_minus;
        std::list<std::pair<int, double>>::iterator iterminus = particle_list_minus.begin();
        // take first values of list
        // std::cout << " numdelpartminus " << numdelpartminus << " current_number_particles_minus "
        // << current_number_particles_minus << " num_particles_minus " << num_particles_minus <<
        // std::endl;
        for (int ipart = 0; ipart < numdelpartminus; ipart++)
        {
          if (iterminus->second >=
              0.0)  // also delete particles with r_min, as they can be reconstructed
          {
            deleteparticles.push_back(iterminus->first);
            // std::cout << " minus " << iterminus->first << std::endl;
            iterminus++;
          }
          else
            break;  // all the remaining particles are either escaped or are closer to the interface
                    // than r_min
        }
      }
    }
    else if ((numnodes_in_band == numnode and intersected == false)  // bin is completely in band
             or (numnodes_in_band > 0 and numnodes_in_band < 8 and
                    intersected == false))  // bin is partially in band
    {
      int numparticle = actbin->NumNode();

      // delete particles outside of band (partially-in-band bins only)
      if (numnodes_in_band > 0 and numnodes_in_band < 8 and intersected == false)
      {
        // loop particles to determine useless particles that should be deleted
        // also count number of particles that will be deleted
        int counter = 0;
        for (int ipart = 0; ipart < actbin->NumNode(); ipart++)
        {
          // get current particle
          const DRT::Node* currentparticle = (actbin->Nodes())[ipart];
          // get gid and lid
          const int gid = currentparticle->Id();
          const int lid = noderowmap->LID(gid);

          // get phi and gradient of particle at current particle position
          double current_phi_particle = 0.0;
          // dummy vector
          LINALG::Matrix<3, 1> normal_particle(true);

          // element coordinates of particle position in scatra element
          LINALG::Matrix<3, 1> elecoord(true);
          DRT::Element* scatraele = GetEleCoordinates(gid, elecoord);

          // compute level-set value
          GetLSValParticle(
              current_phi_particle, normal_particle, scatraele, elecoord, phinp, false);

          // check for non-escaped particles out of band only
          if (((*sign)[lid] * current_phi_particle) > b_max_)
          {
            deleteparticles.push_back(gid);
            counter += 1;
          }
        }

        // get remaining number of particles
        numparticle = actbin->NumNode() - counter;
      }

      // particles are missing (it is not distinguished between escaped and non-escaped particles)
      if (numparticle < (2 * num_particles_per_bin_))
      {
        std::pair<int, int> seed_particles(0, 0);
        if (phi_ref > 0.0)
          seed_particles.first = (2 * num_particles_per_bin_) - numparticle;
        else
          seed_particles.second = (2 * num_particles_per_bin_) - numparticle;

        std::pair<int, std::pair<int, int>> fillbin(actbin->Id(), seed_particles);
        particlebins.insert(fillbin);
      }
      // too many particles
      else if (numparticle > (2 * num_particles_per_bin_))
      {
        // list for particles sorted by distance from interface
        // pair< int eleGID, double sortcriterion >
        std::list<std::pair<int, double>> particle_list;

        for (int ipart = 0; ipart < actbin->NumNode(); ipart++)
        {
          // get current particle
          const DRT::Node* currentparticle = (actbin->Nodes())[ipart];
          // get gid and lid
          const int gid = currentparticle->Id();
          const int lid = noderowmap->LID(gid);
          // get sign
          const double signP = (*sign)[lid];
          // and radius
          const double radiusP = (*radius)[lid];

          // get phi and gradient of particle at current particle position
          double current_phi_particle = 0.0;
          // dummy vector
          LINALG::Matrix<3, 1> normal_particle(true);

          // element coordinates of particle position in scatra element
          LINALG::Matrix<3, 1> elecoord(true);
          DRT::Element* scatraele = GetEleCoordinates(gid, elecoord);

          // compute level-set value
          GetLSValParticle(
              current_phi_particle, normal_particle, scatraele, elecoord, phinp, false);

          if (numnodes_in_band > 0 and numnodes_in_band < 8 and
              intersected == false)  //(partially-in-band bins only)
          {
            // is particle already added to list of deleted particles
            if (((*sign)[lid] * current_phi_particle) > b_max_) continue;
          }

          // insert particle
          std::pair<int, double> thispair;
          thispair.first = gid;
          // compute sort criterion
          thispair.second = signP * current_phi_particle - radiusP;
          particle_list.push_back(thispair);
        }

        // sort particles in descending order
        // for escaped particles, which are not excluded, sort criterion yields negative value,
        // and they are stored at the bottom
        // this is the STL sorting, which is pretty fast
        particle_list.sort(MyComparePairs);
        // for (std::list< std::pair< int, double > >::iterator rr=particle_list.begin();
        // rr!=particle_list.end(); rr++)
        //   std::cout << "id  " << rr->first << " val " << rr->second << std::endl;

        // number of particles to be deleted
        const int numdelpart = numparticle - (2 * num_particles_per_bin_);
        std::list<std::pair<int, double>>::iterator iter = particle_list.begin();
        // take first values of list
        for (int ipart = 0; ipart < numdelpart; ipart++)
        {
          if (iter->second >=
              0.0)  // also delete particles with r_min, as they can be reconstructed
          {
            deleteparticles.push_back(iter->first);
            iter++;
          }
          else
            break;  // all the remaining particles are either escaped or are closer to the interface
                    // than r_min
        }
      }
    }
    // else dserror("Unknown situation!");
    else
    {
      std::cout << "## WARNING: special case in reseeding function" << std::endl;
      // remark: this special case has to taken into account, since it may occur in underresolved
      // thin filaments as
      //         they may arise for the single vortex problem; in this case proper reinitialization
      //         can not be ensured, and hence such cases emerge

      for (std::size_t rr = 0; rr < myphinp.size(); rr++) std::cout << myphinp[rr] << std::endl;

      std::cout << "number of nodes in band  " << numnodes_in_band << "  number of nodes  "
                << numnode << "  intersected  " << intersected << std::endl;
    }

  }  // end: loop all bins

  // delete also particles that have escaped by more than (three times) their radius
  // in case of bubbles their rather destroy the interface than being helpful
  // in particular, if spurious velocities cause their escape
  // for problems with merging interfaces
  if (delete_more_ > 1.0)
  {
    for (int inode = 0; inode < BinStrategy()->BinDiscret()->NumMyRowNodes(); inode++)
    {
      // get gid
      const int mygid = BinStrategy()->BinDiscret()->NodeRowMap()->GID(inode);
      // get phi of particle at current particle position
      double myphi = 0.0;
      // we do not need the gradient here, so this is just a dummy
      LINALG::Matrix<3, 1> normal_particle(true);

      // element coordinates of particle position in scatra element
      LINALG::Matrix<3, 1> elecoord(true);
      DRT::Element* scatraele = GetEleCoordinates(mygid, elecoord);

      // compute level-set value (false=no computation of gradient)
      GetLSValParticle(myphi, normal_particle, scatraele, elecoord, phinp, false);

      // get sign
      const double mysign = (*sign)[inode];

      // check for escaped particles only
      if ((mysign * myphi) < 0.0)
      {
        // delete particles that have escaped by more than "delete_more_" times their radius
        if (std::abs(myphi) >
            (delete_more_ * r_min_))  // note: escaped particles have minimal radius
        {
          bool add = true;
          for (std::size_t kk = 0; kk < deleteparticles.size(); kk++)
          {
            if (deleteparticles[kk] == mygid)
            {
              std::cout << "Particle " << mygid << " with phi " << myphi << " and sign " << mysign
                        << "is already in list!" << std::endl;
              add = false;
              // taking the dserror solution would be better
              // dserror("already in list!");
            }
          }

          if (add) deleteparticles.push_back(mygid);
        }
      }
    }
  }


  // -------------------------------------------------------------------
  //           delete useless particles
  // -------------------------------------------------------------------

  // some output and checks
  int mydelpart = deleteparticles.size();
  int alldelpart = 0;
  scatradis_->Comm().SumAll(&mydelpart, &alldelpart, 1);
  if (MyRank() == 0)
  {
    std::cout << "--- current number of particles " << BinStrategy()->BinDiscret()->NumGlobalNodes()
              << std::endl;
    std::cout << "--- number of particles to be deleted " << alldelpart << std::endl;
  }

  // remove particles (note: calls FillComplete)
  DeleteParticles(deleteparticles);

  // some output and checks
  if (MyRank() == 0)
    std::cout << "--- current number of particles " << BinStrategy()->BinDiscret()->NumGlobalNodes()
              << std::endl;

  // -------------------------------------------------------------------
  //             seed new particles
  // -------------------------------------------------------------------

  // initialize particle id with largest particle id in use + 1 (on each proc)
  int maxparticleid = BinStrategy()->BinDiscret()->NodeRowMap()->MaxAllGID();
  int currentparticleid = maxparticleid;

  // compute local offset for global id numbering
  int myoffset = 0;

  // communicate number of expected new particles of each proc to all procs
  const int numproc = BinStrategy()->BinDiscret()->Comm().NumProc();

  for (int iproc = 0; iproc < numproc; ++iproc)
  {
    int numparts = 0;
    std::map<int, std::pair<int, int>>::iterator iter = particlebins.begin();
    for (iter = particlebins.begin(); iter != particlebins.end(); iter++)
      numparts += (iter->second).first + (iter->second).second;

    // std::cout << "myrank  " << MyRank() << "   numparts  " << numparts << std::endl;

    BinStrategy()->BinDiscret()->Comm().Broadcast(&numparts, 1, iproc);
    if (MyRank() > iproc) myoffset += numparts;
  }
  // std::cout << "myrank  " << MyRank() << "   offset  " << myoffset << std::endl;

  // added offset to current particle id
  currentparticleid += myoffset;

  // list of all particles and their position
  std::map<int, std::vector<double>> particle_positions;
  // list of all particles and their sign
  std::map<int, double> particle_sign;

  // loop all bins which should be filled
  for (std::map<int, std::pair<int, int>>::iterator it = particlebins.begin();
       it != particlebins.end(); it++)
  {
    // get center of current bin
    LINALG::Matrix<3, 1> center(true);
    center = BinStrategy()->GetBinCentroid(it->first);

    // get max and min values of bin corners
    std::vector<LINALG::Matrix<3, 1>> bincorners(2);
    for (int rr = 0; rr < 2; rr++)
    {
      for (int idim = 0; idim < 3; idim++)
        (bincorners[rr])(idim, 0) =
            center(idim, 0) + pow(-1.0, (double)rr + 1.0) * 0.5 * BinStrategy()->BinSize()[idim];
    }

    // loop number of particles to be seeded
    const int mypartplus = (it->second).first;
    const int mypartminus = (it->second).second;
    const int mypart = mypartplus + mypartminus;
    for (int ipart = 0; ipart < mypart; ipart++)
    {
      // set particle id
      currentparticleid += 1;

      // get position randomly in bin
      std::vector<double> position(3);

      DRT::UTILS::Random* random = DRT::Problem::Instance()->Random();

      for (int idim = 0; idim < 3; idim++)
      {
        // set range (default: [-1;1])
        random->SetRandRange((bincorners[0])(idim, 0), (bincorners[1])(idim, 0));
        // get position
        position[idim] = random->Uni();
      }

      // set particle to mid plane in case of quasi-2D simulations
      switch (BinStrategy()->ParticleDim())
      {
        case INPAR::PARTICLEOLD::particle_2Dx:
        {
          position[0] = 0.0;
          break;
        }
        case INPAR::PARTICLEOLD::particle_2Dy:
        {
          position[1] = 0.0;
          break;
        }
        case INPAR::PARTICLEOLD::particle_2Dz:
        {
          position[2] = 0.0;
          break;
        }
        default:
          break;
      }

      // place particle with id and position
      std::list<Teuchos::RCP<DRT::Node>> homelessparticles;
      Teuchos::RCP<DRT::Node> newparticle =
          Teuchos::rcp(new DRT::Node(currentparticleid, &position[0], MyRank()));
      PlaceNodeCorrectly(newparticle, &position[0], homelessparticles);
      if (homelessparticles.size() != 0)
        dserror("New particle could not be inserted on this proc!.");

      // add particle to list
      particle_positions.insert(std::pair<int, std::vector<double>>(currentparticleid, position));
      // set sign (+ for the first half, - for the second one)
      if (ipart < mypartplus)
        particle_sign.insert(std::pair<int, double>(currentparticleid, 1.0));
      else
        particle_sign.insert(std::pair<int, double>(currentparticleid, -1.0));
    }
  }

  // rebuild connectivity and assign degrees of freedom (note: IndependentDofSet)
  BinStrategy()->BinDiscret()->FillComplete(true, false, true);

  // update of state vectors to the new maps
  particles_->UpdateStatesAfterParticleTransfer();

  // -------------------------------------------------------------------
  //               initialize state vectors
  // -------------------------------------------------------------------

  // get position vector at time_(n+1)
  Teuchos::RCP<Epetra_Vector> disnp = particles_->WriteAccessDispnp();
  // get sign vector
  sign = Teuchos::rcp_dynamic_cast<PARTICLE::TimIntRK>(particles_)->WriteAccessSign();
  // get radius vector
  radius = particles_->WriteAccessRadiusn();

  // get maps
  const Epetra_Map* dofrowmap = BinStrategy()->BinDiscret()->DofRowMap();
  // update node row map due to FillCompleteCall
  noderowmap = BinStrategy()->BinDiscret()->NodeRowMap();

  // prepare also attraction
  // noderowmap-based vectors
  // predictor
  Teuchos::RCP<Epetra_Vector> lambda = Teuchos::rcp(new Epetra_Vector(*noderowmap, true));
  // initialize with 0.0 to ensure that old particles will not be moved
  lambda->PutScalar(0.0);
  // target phi_goal
  Teuchos::RCP<Epetra_Vector> phi_target = Teuchos::rcp(new Epetra_Vector(*noderowmap, true));
  // current phi of particle
  Teuchos::RCP<Epetra_Vector> phi_particle = Teuchos::rcp(new Epetra_Vector(*noderowmap, true));
  // dofrowmap-based vectors
  // normal vector of level-set field at particle position
  Teuchos::RCP<Epetra_Vector> norm_gradphi_particle =
      Teuchos::rcp(new Epetra_Vector(*dofrowmap, true));
  // increment position vector for attraction loop
  Teuchos::RCP<Epetra_Vector> inc_dis = Teuchos::rcp(new Epetra_Vector(*dofrowmap, true));
  // initialize with 0.0 to ensure that old particles will not be moved
  inc_dis->PutScalar(0.0);

  // check sizes particle_positions and particle_sign
  if (particle_positions.size() != particle_sign.size()) dserror("Same length expected");

  for (std::map<int, std::vector<double>>::iterator it = particle_positions.begin();
       it != particle_positions.end(); it++)
  {
    // get global id of current particle
    const int gid = it->first;

    DRT::Node* currparticle = BinStrategy()->BinDiscret()->gNode(gid);
    // get the first dof gid of a particle and convert it into a LID
    const int doflid = dofrowmap->LID(BinStrategy()->BinDiscret()->Dof(currparticle, 0));
    // get particle lid
    const int partlid = noderowmap->LID(gid);

    // initialize state vectors

    // insert values into state vectors
    (*sign)[partlid] = particle_sign[gid];
    for (int idim = 0; idim < 3; ++idim) (*disnp)[doflid + idim] = (particle_positions[gid])[idim];
    // reset to zero just to be sure
    (*radius)[partlid] = 0.0;

    // prepare attraction

    // set lambda to 1.0 to start attraction
    (*lambda)[partlid] = 1.0;

    // get sign of particle
    const double signP = (*sign)[partlid];

    // get target phi for this particle
    DRT::UTILS::Random* random = DRT::Problem::Instance()->Random();
    random->SetRandRange(r_min_, b_max_);
    const double local_phi_target = signP * random->Uni();
    (*phi_target)[partlid] = local_phi_target;

    // get phi and gradient of particle at current particle position
    double current_phi_particle = 0.0;
    LINALG::Matrix<3, 1> normal_particle(true);

    // element coordinates of particle position in scatra element
    LINALG::Matrix<3, 1> elecoord(true);
    DRT::Element* scatraele = GetEleCoordinates(gid, elecoord);

    // compute level-set value
    GetLSValParticle(current_phi_particle, normal_particle, scatraele, elecoord, phinp);
    // GetLSValParticle returns gradient of phi, which we have to normalize here
    const double norm = normal_particle.Norm2();
    if (norm > 1.0e-9)
      normal_particle.Scale(1.0 / norm);
    else
      normal_particle.PutScalar(0.0);

    // store values
    (*phi_particle)[partlid] = current_phi_particle;

    for (int idim = 0; idim < 3; ++idim)
      (*norm_gradphi_particle)[doflid + idim] = normal_particle(idim, 0);

    // initialize increment position vector, i.e., multiply by
    // lambda*(phi_target-current_phi_particle) note lambda is set to 1.0 at the beginning and
    // skipped here
    const double fac = local_phi_target - current_phi_particle;
    normal_particle.Scale(fac);
    for (int idim = 0; idim < 3; ++idim) (*inc_dis)[doflid + idim] = normal_particle(idim, 0);
  }

  // initial or intermediate position vector
  Teuchos::RCP<Epetra_Vector> disn = Teuchos::rcp(new Epetra_Vector(*disnp));

  // -------------------------------------------------------------------
  //               attraction step
  // -------------------------------------------------------------------
  // move particle to the correct side of the interface

  // total (global) number of reseeded particles
  int my_new_particles = particle_positions.size();
  int global_num_particles = 0;
  scatradis_->Comm().SumAll(&my_new_particles, &global_num_particles, 1);

  // some output and checks
  if (MyRank() == 0)
  {
    std::cout << "--- total number of reseeded particles  " << global_num_particles << std::endl;
    std::cout << "--- total number of particles  " << BinStrategy()->BinDiscret()->NumGlobalNodes()
              << std::endl;
  }

#if 0
  // number of placed particles should finally be equal to the number of seeded particles
  // important: this number must be global and takes into account the placed particles on
  // all procs
  int num_placed_particles = 0;
  // counter for number of attraction steps
  int step_counter = 0;
  // define vector for update of map-based vectors
  Teuchos::RCP<Epetra_Vector> old;

  // attraction loop
  while ((num_placed_particles < global_num_particles) and step_counter < 15)
  {
    // increment counter
    step_counter += 1;
    if (MyRank() == 0)
      std::cout << "--- attraction step  "<< step_counter << std::endl;

    // compute new position
    disnp = particles_->WriteAccessDispnp();
    // x_P^new = x_P + lambda*(phi_target-current_phi_particle) * n_p
    disnp->Update(1.0,*inc_dis,1.0);

    // transfer particles to bins
    TransferParticles(false);
    // update of state vectors to the new maps
    particles_->UpdateStatesAfterParticleTransfer();
    // likewise update present vectors according to the new distribution of particles
    old = lambda;
    lambda = LINALG::CreateVector(*BinStrategy()->BinDiscret()->NodeRowMap(),true);
    LINALG::Export(*old, *lambda);

    old = phi_target;
    phi_target = LINALG::CreateVector(*BinStrategy()->BinDiscret()->NodeRowMap(),true);
    LINALG::Export(*old, *phi_target);

    old = phi_particle;
    phi_particle = LINALG::CreateVector(*BinStrategy()->BinDiscret()->NodeRowMap(),true);
    LINALG::Export(*old, *phi_particle);

    old = norm_gradphi_particle;
    norm_gradphi_particle = LINALG::CreateVector(*BinStrategy()->BinDiscret()->DofRowMap(),true);
    LINALG::Export(*old, *norm_gradphi_particle);

    old = inc_dis;
    inc_dis = LINALG::CreateVector(*BinStrategy()->BinDiscret()->DofRowMap(),true);
    LINALG::Export(*old, *inc_dis);

    old = disn;
    disn = LINALG::CreateVector(*BinStrategy()->BinDiscret()->DofRowMap(),true);
    LINALG::Export(*old, *disn);

    // check new position
    // counter for local placed particles
    int local_num_placed_particles = 0;

    // get new distribution of vectors
    disnp = particles_->WriteAccessDispnp();
    sign = Teuchos::rcp_dynamic_cast<PARTICLE::TimIntRK>(particles_)->WriteAccessSign();

    // loop all particles
    for (int inode=0; inode<BinStrategy()->BinDiscret()->NumMyRowNodes(); inode++)
    {
      // get new maps
      const Epetra_Map* att_dofrowmap = BinStrategy()->BinDiscret()->DofRowMap();
      const Epetra_Map* att_noderowmap = BinStrategy()->BinDiscret()->NodeRowMap();

      // get dof lid
      const int gid = att_noderowmap->GID(inode);
      DRT::Node* currparticle = BinStrategy()->BinDiscret()->gNode(gid);
      int doflid = att_dofrowmap->LID(BinStrategy()->BinDiscret()->Dof(currparticle, 0));

      if ((*lambda)[inode] > 0.0)
      {
        // element coordinates of particle position in scatra element
        LINALG::Matrix<3,1> elecoord(true);
        DRT::Element* scatraele = GetEleCoordinates(gid,elecoord);

        if (scatraele != NULL) // particle has not left domain
        {
          // get phi and gradient of particle at current particle position
          double current_phi_particle = 0.0;
          LINALG::Matrix<3,1> normal_particle(true);
          // compute level-set value
          GetLSValParticle(current_phi_particle,normal_particle,scatraele,elecoord,phinp);
          // GetLSValParticle returns gradient of phi, which we have to normalize here
          const double norm = normal_particle.Norm2();
          if (norm > 1.0e-9)
            normal_particle.Scale(1.0/norm);
          else
            normal_particle.PutScalar(0.0);

          const double signP = (*sign)[inode];
          if ( (( signP > 0.0 and (current_phi_particle >= r_min_)
                              and (current_phi_particle <= b_max_))
              or
                ( signP < 0.0 and (current_phi_particle <= (signP * r_min_))
                              and (current_phi_particle >= (signP * b_max_))) )
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

      } // lambda != 0

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
    if (MyRank() == 0)
      std::cout << "--- placed particles  " << num_placed_particles << std::endl;

    // possibly position have been reverted and particles changed the bin or even the processor
    // therefore, we have to recall the following steps
    TransferParticles(false);
    // update of state vectors to the new maps
    particles_->UpdateStatesAfterParticleTransfer();
    // likewise update present vectors according to the new distribution of particles
    Teuchos::RCP<Epetra_Vector> old;
    old = lambda;
    lambda = LINALG::CreateVector(*BinStrategy()->BinDiscret()->NodeRowMap(),true);
    LINALG::Export(*old, *lambda);

    old = phi_target;
    phi_target = LINALG::CreateVector(*BinStrategy()->BinDiscret()->NodeRowMap(),true);
    LINALG::Export(*old, *phi_target);

    old = phi_particle;
    phi_particle = LINALG::CreateVector(*BinStrategy()->BinDiscret()->NodeRowMap(),true);
    LINALG::Export(*old, *phi_particle);

    old = norm_gradphi_particle;
    norm_gradphi_particle = LINALG::CreateVector(*BinStrategy()->BinDiscret()->DofRowMap(),true);
    LINALG::Export(*old, *norm_gradphi_particle);

    old = inc_dis;
    inc_dis = LINALG::CreateVector(*BinStrategy()->BinDiscret()->DofRowMap(),true);
    LINALG::Export(*old, *inc_dis);

    old = disn;
    disn = LINALG::CreateVector(*BinStrategy()->BinDiscret()->DofRowMap(),true);
    LINALG::Export(*old, *disn);

    // update state vectors
    disnp = particles_->WriteAccessDispnp();
    disn->Update(1.0,*disnp,0.0);
  }

  // delete particles which could not be attracted correctly -> all particles with lambda != 0.0
  if (num_placed_particles != global_num_particles) // we have particles with lambda != 0.0
  {
    if (MyRank() == 0)
      std::cout << "--- delete unplaced particles  "<< global_num_particles-num_placed_particles << std::endl;

    // vector with particles to delete
    std::vector<int> part_del;

    // loop all particles
    for (int inode=0; inode<BinStrategy()->BinDiscret()->NumMyRowNodes(); inode++)
    {
      // particles with lambda != 0.0 will be deleted
      if ((*lambda)[inode] > (0.0+1.0e-8))
      {
        // get particle
        DRT::Node *currparticle = BinStrategy()->BinDiscret()->lRowNode(inode);
        // store gid
        part_del.push_back(currparticle->Id());
      }
    }

    // remove particles
    DeleteParticles(part_del);
  }
#else
  // do attraction loop
  Attraction(
      global_num_particles, inc_dis, lambda, phi_target, phi_particle, norm_gradphi_particle, disn);
  // caution: this function changes pointers
#endif

  // -------------------------------------------------------------------
  //               set radius
  // -------------------------------------------------------------------

  AdjustParticleRadii();

  if (MyRank() == 0) std::cout << "----------------------------- " << std::endl;

  return;
}


/*----------------------------------------------------------------------*
 | output particle time step                            rasthofer 03/14 |
 *----------------------------------------------------------------------*/
void PARTICLE::ScatraParticleCoupling::Attraction(const int global_num_particles,
    Teuchos::RCP<Epetra_Vector> inc_dis, Teuchos::RCP<Epetra_Vector> lambda,
    Teuchos::RCP<Epetra_Vector> phi_target, Teuchos::RCP<Epetra_Vector> phi_particle,
    Teuchos::RCP<Epetra_Vector> norm_gradphi_particle, Teuchos::RCP<Epetra_Vector> disn)
{
  // get vectors from time integration
  Teuchos::RCP<Epetra_Vector> disnp = particles_->WriteAccessDispnp();
  Teuchos::RCP<Epetra_Vector> sign =
      Teuchos::rcp_dynamic_cast<PARTICLE::TimIntRK>(particles_)->WriteAccessSign();

  // get phinp
  const Teuchos::RCP<const Epetra_Vector> row_phinp = scatra_->Phinp();
  // export phi vector to col map
  // this is necessary here, since the elements adjacent to a row node may have
  // ghosted nodes, which are only contained in the col format
  const Epetra_Map* scatra_dofcolmap = scatradis_->DofColMap();
  Teuchos::RCP<Epetra_Vector> phinp = Teuchos::rcp(new Epetra_Vector(*scatra_dofcolmap));
  LINALG::Export(*row_phinp, *phinp);

  // number of placed particles should finally be equal to the number of seeded particles (here
  // global_num_particles) important: this number must be global and takes into account the placed
  // particles on all procs
  int num_placed_particles = 0;
  // counter for number of attraction steps
  int step_counter = 0;
  // define vector for update of map-based vectors
  Teuchos::RCP<Epetra_Vector> old;

  // attraction loop
  while ((num_placed_particles < global_num_particles) and step_counter < 15)
  {
    // increment counter
    step_counter += 1;
    if (MyRank() == 0) std::cout << "--- attraction step  " << step_counter << std::endl;

    // compute new position
    disnp = particles_->WriteAccessDispnp();
    // x_P^new = x_P + lambda*(phi_target-current_phi_particle) * n_p
    disnp->Update(1.0, *inc_dis, 1.0);

    // set current position to particle nodes
    SetParticleNodePos();

    // transfer particles to bins
    TransferParticles(true, false);

    // update present vectors according to the new distribution of particles
    old = lambda;
    lambda = LINALG::CreateVector(*BinStrategy()->BinDiscret()->NodeRowMap(), true);
    LINALG::Export(*old, *lambda);

    old = phi_target;
    phi_target = LINALG::CreateVector(*BinStrategy()->BinDiscret()->NodeRowMap(), true);
    LINALG::Export(*old, *phi_target);

    old = phi_particle;
    phi_particle = LINALG::CreateVector(*BinStrategy()->BinDiscret()->NodeRowMap(), true);
    LINALG::Export(*old, *phi_particle);

    old = norm_gradphi_particle;
    norm_gradphi_particle = LINALG::CreateVector(*BinStrategy()->BinDiscret()->DofRowMap(), true);
    LINALG::Export(*old, *norm_gradphi_particle);

    old = inc_dis;
    inc_dis = LINALG::CreateVector(*BinStrategy()->BinDiscret()->DofRowMap(), true);
    LINALG::Export(*old, *inc_dis);

    old = disn;
    disn = LINALG::CreateVector(*BinStrategy()->BinDiscret()->DofRowMap(), true);
    LINALG::Export(*old, *disn);

    // check new position
    // counter for local placed particles
    int local_num_placed_particles = 0;

    // get new distribution of vectors
    disnp = particles_->WriteAccessDispnp();
    sign = Teuchos::rcp_dynamic_cast<PARTICLE::TimIntRK>(particles_)->WriteAccessSign();

    // loop all particles
    for (int inode = 0; inode < BinStrategy()->BinDiscret()->NumMyRowNodes(); inode++)
    {
      // get new maps
      const Epetra_Map* att_dofrowmap = BinStrategy()->BinDiscret()->DofRowMap();
      const Epetra_Map* att_noderowmap = BinStrategy()->BinDiscret()->NodeRowMap();

      // get dof lid
      const int gid = att_noderowmap->GID(inode);
      DRT::Node* currparticle = BinStrategy()->BinDiscret()->gNode(gid);
      int doflid = att_dofrowmap->LID(BinStrategy()->BinDiscret()->Dof(currparticle, 0));

      if ((*lambda)[inode] > 0.0)
      {
        // element coordinates of particle position in scatra element
        LINALG::Matrix<3, 1> elecoord(true);
        DRT::Element* scatraele = GetEleCoordinates(gid, elecoord);

        if (scatraele != NULL)  // particle has not left domain
        {
          // get phi and gradient of particle at current particle position
          double current_phi_particle = 0.0;
          LINALG::Matrix<3, 1> normal_particle(true);
          // compute level-set value
          GetLSValParticle(current_phi_particle, normal_particle, scatraele, elecoord, phinp);
          // GetLSValParticle returns gradient of phi, which we have to normalize here
          const double norm = normal_particle.Norm2();
          if (norm > 1.0e-9)
            normal_particle.Scale(1.0 / norm);
          else
            normal_particle.PutScalar(0.0);

          const double signP = (*sign)[inode];
          if (((signP > 0.0 and (current_phi_particle >= r_min_) and
                   (current_phi_particle <= b_max_)) or
                  (signP < 0.0 and (current_phi_particle <= (signP * r_min_)) and
                      (current_phi_particle >= (signP * b_max_)))) and
              ((*lambda)[inode] == 1.0))  // particle in correct band
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
            if ((*lambda)[inode] == 1.0)  // we could not move the particle to the correct band
            {
              // we go back to our initial position
              // we do not overwrite phi_particle and norm_gradphi_particle
              // we set disnp back to disn
              for (int idim = 0; idim < 3; ++idim) (*disnp)[doflid + idim] = (*disn)[doflid + idim];
              // and only go half of the distance in the next step
              (*lambda)[inode] = 0.5;
            }
            else  // we accept the new position and start the attraction step from this point with
                  // lambda again set to 1.0
            {
              // store new phi value
              (*phi_particle)[inode] = current_phi_particle;
              // set lamda to 1.0
              (*lambda)[inode] = 1.0;
              // store new gradient
              for (int idim = 0; idim < 3; ++idim)
                (*norm_gradphi_particle)[doflid + idim] = normal_particle(idim, 0);
            }
          }
        }
        else
        {
          std::cout << "Special case during attraction: particle left domain!" << std::endl;
          // if we are outside of the domain, we go back to our initial position
          // we do not overwrite phi_particle and norm_gradphi_particle
          // we set disnp back to disn
          for (int idim = 0; idim < 3; ++idim) (*disnp)[doflid + idim] = (*disn)[doflid + idim];
          // and start another attraction step with lambda halved
          double mylambda = (*lambda)[inode];
          (*lambda)[inode] = 0.5 * mylambda;
        }

      }  // lambda != 0

      // set new increment vector
      const double fac = (*lambda)[inode] * ((*phi_target)[inode] - (*phi_particle)[inode]);
      for (int idim = 0; idim < 3; ++idim)
        (*inc_dis)[doflid + idim] = fac * (*norm_gradphi_particle)[doflid + idim];

    }  // end loop all particles

    // collect correctly placed particles of this step form all procs
    int global_num_placed_particles = 0;
    scatradis_->Comm().SumAll(&local_num_placed_particles, &global_num_placed_particles, 1);
    // update number of placed particles
    num_placed_particles += global_num_placed_particles;
    if (MyRank() == 0) std::cout << "--- placed particles  " << num_placed_particles << std::endl;

    // possibly position have been reverted and particles changed the bin or even the processor
    // therefore, we have to recall the following steps
    SetParticleNodePos();
    TransferParticles(true, false);
    // likewise update present vectors according to the new distribution of particles
    Teuchos::RCP<Epetra_Vector> old;
    old = lambda;
    lambda = LINALG::CreateVector(*BinStrategy()->BinDiscret()->NodeRowMap(), true);
    LINALG::Export(*old, *lambda);

    old = phi_target;
    phi_target = LINALG::CreateVector(*BinStrategy()->BinDiscret()->NodeRowMap(), true);
    LINALG::Export(*old, *phi_target);

    old = phi_particle;
    phi_particle = LINALG::CreateVector(*BinStrategy()->BinDiscret()->NodeRowMap(), true);
    LINALG::Export(*old, *phi_particle);

    old = norm_gradphi_particle;
    norm_gradphi_particle = LINALG::CreateVector(*BinStrategy()->BinDiscret()->DofRowMap(), true);
    LINALG::Export(*old, *norm_gradphi_particle);

    old = inc_dis;
    inc_dis = LINALG::CreateVector(*BinStrategy()->BinDiscret()->DofRowMap(), true);
    LINALG::Export(*old, *inc_dis);

    old = disn;
    disn = LINALG::CreateVector(*BinStrategy()->BinDiscret()->DofRowMap(), true);
    LINALG::Export(*old, *disn);

    // update state vectors
    disnp = particles_->WriteAccessDispnp();
    disn->Update(1.0, *disnp, 0.0);
  }

  // delete particles which could not be attracted correctly -> all particles with lambda != 0.0
  if (num_placed_particles != global_num_particles)  // we have particles with lambda != 0.0
  {
    if (MyRank() == 0)
      std::cout << "--- delete unplaced particles  " << global_num_particles - num_placed_particles
                << std::endl;

    // vector with particles to delete
    std::vector<int> part_del;

    // loop all particles
    for (int inode = 0; inode < BinStrategy()->BinDiscret()->NumMyRowNodes(); inode++)
    {
      // particles with lambda != 0.0 will be deleted
      if ((*lambda)[inode] > (0.0 + 1.0e-8))
      {
        // get particle
        DRT::Node* currparticle = BinStrategy()->BinDiscret()->lRowNode(inode);
        // store gid
        part_del.push_back(currparticle->Id());
      }
    }

    // remove particles
    DeleteParticles(part_del);
  }

  return;
}


/*----------------------------------------------------------------------*
 | delete particles                                     rasthofer 03/14 |
 *----------------------------------------------------------------------*/
void PARTICLE::ScatraParticleCoupling::DeleteParticles(
    std::vector<int> part_del  ///< global ids of particles to be deleted
)
{
  for (std::size_t ipart = 0; ipart < part_del.size(); ipart++)
  {
    // get gid
    const double gid = part_del[ipart];
    // get particle
    DRT::Node* currparticle = BinStrategy()->BinDiscret()->gNode(gid);
    if (currparticle->NumElement() != 1)
      dserror("ERROR: A particle is assigned to more than one bin!");
    // get corresponding bin
    DRT::Element** currele = currparticle->Elements();
    DRT::Element* currbin = currele[0];

    // remove particle from bin
    dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currbin)->DeleteNode(gid);

    // remove particle from discretization
    BinStrategy()->BinDiscret()->DeleteNode(gid);
  }

  // finish discretization
  BinStrategy()->BinDiscret()->FillComplete(true, false, false);
  // adapt state vectors to new maps
  particles_->UpdateStatesAfterParticleTransfer();

  return;
}


#if 0
/*----------------------------------------------------------------------*
 | get underlying element as well as particle position   rasthofer 12/13 |
 *----------------------------------------------------------------------*/
DRT::Element* PARTICLE::ScatraParticleCoupling::GetEleCoordinates(
  const int particle_gid,             ///< global particle id
  LINALG::Matrix<3,1>& elecoord       ///< matrix to be filled with particle coordinates in element space
  )
{
  // get particle corresponding to lid
  DRT::Node* currparticle = BinStrategy()->BinDiscret()->gNode(particle_gid);

  Teuchos::RCP<const Epetra_Vector> particlepos = particles_->Dispnp();
  // fill particle position
  LINALG::Matrix<3,1> particleposition;
  std::vector<int> lm_b = BinStrategy()->BinDiscret()->Dof(currparticle);
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
  DRT::MESHFREE::MeshfreeMultiBin* currbin = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currele[0]);
  DRT::Element** scatraelesinbin = currbin->AssociatedEles(bin_transpcontent_);
  int numscatraelesinbin = currbin->NumAssociatedEle(bin_transpcontent_);

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
  // remark: if an extended (bin) domain is considered to keep particles that leave the domain
  //         during attraction, this error should be provoked

  return targetscatraele;
}
#endif


/*----------------------------------------------------------------------*
 | get underlying element as well as particle position   rasthofer 12/13 |
 *----------------------------------------------------------------------*/
DRT::Element* PARTICLE::ScatraParticleCoupling::GetEleCoordinates(
    const int particle_gid,  ///< global particle id
    LINALG::Matrix<3, 1>&
        elecoord  ///< matrix to be filled with particle coordinates in element space
)
{
  // get particle corresponding to lid
  DRT::Node* currparticle = BinStrategy()->BinDiscret()->gNode(particle_gid);

  Teuchos::RCP<const Epetra_Vector> particlepos = particles_->Dispnp();
  // fill particle position
  LINALG::Matrix<3, 1> particleposition;
  std::vector<int> lm_b = BinStrategy()->BinDiscret()->Dof(currparticle);
  int posx = particlepos->Map().LID(lm_b[0]);
  for (int idim = 0; idim < 3; ++idim) particleposition(idim, 0) = (*particlepos)[posx + idim];

  // find out in which scatra element the current particle is located
  if (currparticle->NumElement() != 1)
    dserror("ERROR: A particle is assigned to more than one bin!");

  DRT::Element** currele = currparticle->Elements();
#ifdef DEBUG
  DRT::MESHFREE::MeshfreeMultiBin* test =
      dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currele[0]);
  if (test == NULL)
    dserror("dynamic cast from DRT::Element to DRT::MESHFREE::MeshfreeMultiBin failed");
#endif
  DRT::MESHFREE::MeshfreeMultiBin* currbin =
      dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currele[0]);

  return GetEleCoordinatesFromPosition(particleposition, currbin, elecoord);
}


/*----------------------------------------------------------------------*
 | get underlying element as well as position in element space          |
 |                                                      rasthofer 06/14 |
 *----------------------------------------------------------------------*/
DRT::Element* PARTICLE::ScatraParticleCoupling::GetEleCoordinatesFromPosition(
    LINALG::Matrix<3, 1>& myposition,          ///< position
    DRT::MESHFREE::MeshfreeMultiBin* currbin,  ///< corresponding bin
    LINALG::Matrix<3, 1>&
        elecoord  ///< matrix to be filled with particle coordinates in element space
)
{
  // variables to store information about element in which the particle is located
  DRT::Element* targetscatraele = NULL;
  elecoord.Clear();

  DRT::Element** scatraelesinbin = currbin->AssociatedEles(bin_transpcontent_);
  int numscatraelesinbin = currbin->NumAssociatedEle(bin_transpcontent_);

  // search for underlying scatra element with standard search in case nothing was found
  for (int ele = 0; ele < numscatraelesinbin; ++ele)
  {
    DRT::Element* scatraele = scatraelesinbin[ele];
    const LINALG::SerialDenseMatrix xyze(GEO::InitialPositionArray(scatraele));

    // get coordinates of the particle position in parameter space of the element
    bool insideele = GEO::currentToVolumeElementCoordinates(
        scatraele->Shape(), xyze, myposition, elecoord, false);

    if (insideele == true)
    {
      targetscatraele = scatraele;
      // leave loop over all scatra eles in bin
      break;
    }
  }

  if (targetscatraele == NULL)
  {
    //    std::cout << XAABB_ << "   " << myposition << std::endl;
    dserror("Could not found corresponding element!");
  }
  // remark: if an extended (bin) domain is considered to keep particles that leave the domain
  //         during attraction, this error should be provoked

  return targetscatraele;
}


// FOR EXTENSION TO OTHER ELEMENT TYPES: consider template version
// void PARTICLE::ScatraParticleCoupling::GetLSValParticle(
//  double& phi_particle,               ///< phi value at particle position
//  LINALG::Matrix<3,1>& gradphi,       ///< gradient of level-set field at particle position
//  DRT::Element* scatraele,            ///< corresponding scatra element
//  const LINALG::Matrix<3,1> elecoord, ///< matrix containing particle coordinates in element space
//  const Teuchos::RCP<const Epetra_Vector> phinp, ///< vector containing level-set values (col map)
//  bool compute_normal                 ///< additionally compute normal vector Optional -> default
//  = true)
//  )
//{
//  switch (scatraele->Shape())
//  {
//    case  DRT::Element::hex8:
//    {
//      GetLSValParticle<DRT::Element::hex8>(phi_particle,gradphi,scatraele,elecoord,phinp,compute_normal);
//      break;
//    }
//    default:
//    {
//      dserror("Add your element type!");
//      break;
//    }
//  }
//
//  return;
//}
// also for GetVelParticle


/*----------------------------------------------------------------------*
 | compute level-set value of particle                  rasthofer 12/13 |
 *----------------------------------------------------------------------*/
void PARTICLE::ScatraParticleCoupling::GetLSValParticle(
    double& phi_particle,           ///< phi value at particle position
    LINALG::Matrix<3, 1>& gradphi,  ///< gradient of level-set field at particle position
    DRT::Element* scatraele,        ///< corresponding scatra element
    const LINALG::Matrix<3, 1>
        elecoord,  ///< matrix containing particle coordinates in element space
    const Teuchos::RCP<const Epetra_Vector>
        phinp,           ///< vector containing level-set values (col map)
    bool compute_normal  ///< additionally compute normal vector Optional -> default = true)
)
{
  if (scatraele->Shape() != DRT::Element::hex8)
    dserror("Other element than hex8 not yet supported!");  // -> for adaption see comment above
  static const int numnode =
      DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;

  // get nodal phi values
  std::vector<int> lm;
  std::vector<int> lmowner;
  std::vector<int> lmstride;
  lm.clear();
  lmowner.clear();
  lmstride.clear();
  scatraele->LocationVector(*scatradis_, lm, lmowner, lmstride);
  static LINALG::Matrix<numnode, 1> ephinp;
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<numnode, 1>>(*phinp, ephinp, lm);

  // shape functions
  static LINALG::Matrix<numnode, 1> funct;
  DRT::UTILS::shape_function_3D(funct, elecoord(0), elecoord(1), elecoord(2), scatraele->Shape());

  // finally compute phi
  phi_particle = funct.Dot(ephinp);

  // compute normal vector if required
  if (compute_normal)
  {
    static LINALG::Matrix<3, numnode> xyze;
    GEO::fillInitialPositionArray<DRT::Element::hex8, 3, LINALG::Matrix<3, numnode>>(
        scatraele, xyze);

    static LINALG::Matrix<3, numnode> deriv;
    DRT::UTILS::shape_function_3D_deriv1(
        deriv, elecoord(0), elecoord(1), elecoord(2), scatraele->Shape());
    // get transposed of the jacobian matrix d x / d \xi
    // xjm(i,j) = deriv(i,k)*xyze(j,k)
    static LINALG::Matrix<3, 3> xjm(true);
    xjm.MultiplyNT(deriv, xyze);
    const double det = xjm.Determinant();
    if (det < 0.0)
      dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", scatraele->Id(), det);
    // inverse of jacobian
    static LINALG::Matrix<3, 3> xji;
    xji.Invert(xjm);
    // compute global derivates
    static LINALG::Matrix<3, numnode> derxy;
    // derxy(i,j) = xji(i,k) * deriv(k,j)
    derxy.Multiply(xji, deriv);

    // get gradient of phi
    gradphi.Clear();
    gradphi.Multiply(derxy, ephinp);
  }

  return;
}


/*----------------------------------------------------------------------*
 | compute velocity of particle                         rasthofer 01/14 |
 *----------------------------------------------------------------------*/
void PARTICLE::ScatraParticleCoupling::GetVelParticle(
    LINALG::Matrix<3, 1>& vel,  ///< velocity of level-set field at particle position
    DRT::Element* scatraele,    ///< corresponding scatra element
    const LINALG::Matrix<3, 1>
        elecoord,  ///< matrix containing particle coordinates in element space
    const Teuchos::RCP<Epetra_Vector> lsvel  ///< vector containing level-set values (colmap)
)
{
  if (scatraele->Shape() != DRT::Element::hex8)
    dserror("Other element than hex8 not yet supported!");
  static const int numnode =
      DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;

  // construct location vector for velocity related dofs on scatra discretization
  DRT::Element::LocationArray la(scatra_->Discretization()->NumDofSets());
  scatraele->LocationVector(*scatra_->Discretization(), la, false);
  const int numveldofpernode = la[scatra_->NdsVel()].lm_.size() / numnode;
  std::vector<int> lmvel(3 * numnode, -1);
  for (int inode = 0; inode < numnode; ++inode)
    for (int idim = 0; idim < 3; ++idim)
      lmvel[inode * 3 + idim] = la[scatra_->NdsVel()].lm_[inode * numveldofpernode + idim];

  // get nodal values of velocity field from secondary dofset on scatra discretization
  // extract local values from global vector
  static LINALG::Matrix<3, numnode> convel;
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<3, numnode>>(*lsvel, convel, lmvel);

  // shape functions
  static LINALG::Matrix<numnode, 1> funct;
  // fill vectors
  DRT::UTILS::shape_function_3D(funct, elecoord(0), elecoord(1), elecoord(2), scatraele->Shape());

  // finally compute vel
  vel.Multiply(convel, funct);

  return;
}


/*----------------------------------------------------------------------*
 | compute velocities of particles at intermediate time n+theta         |
 |                                                      rasthofer 01/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> PARTICLE::ScatraParticleCoupling::GetVelocity(const double theta)
{
  // inialize vector for particle velocities to be filled
  const Epetra_Map* dofrowmap = BinStrategy()->BinDiscret()->DofRowMap();
  const Epetra_Map* noderowmap = BinStrategy()->BinDiscret()->NodeRowMap();
  Teuchos::RCP<Epetra_Vector> vel = LINALG::CreateVector(*dofrowmap, true);

  // get velocity field from scatra
  const Teuchos::RCP<Epetra_Vector> lsvel =
      Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(scatra_)->ConVelTheta(theta);

  // loop all particles
  for (int inode = 0; inode < BinStrategy()->BinDiscret()->NumMyRowNodes(); inode++)
  {
    // get global id of current particle
    const int gid = noderowmap->GID(inode);

    // element coordinates of particle position in scatra element
    LINALG::Matrix<3, 1> elecoord(true);
    DRT::Element* scatraele = GetEleCoordinates(gid, elecoord);
    // particle is in hole of domain, don't move it any more
    // this results in zero velocity, since vel is initialized with zeros
    if (scatraele == NULL) continue;

    // compute velocity
    LINALG::Matrix<3, 1> vel_particle(true);
    GetVelParticle(vel_particle, scatraele, elecoord, lsvel);

    DRT::Node* currparticle = BinStrategy()->BinDiscret()->gNode(gid);
    // get the first dof gid of a particle and convert it into a LID
    int doflid = dofrowmap->LID(BinStrategy()->BinDiscret()->Dof(currparticle, 0));
    for (int idim = 0; idim < 3; ++idim) (*vel)[doflid + idim] = vel_particle(idim, 0);
  }

  return vel;
}


/*----------------------------------------------------------------------*
| setup ghosting of bins, particles & underlying fluid      ghamm 11/12 |
 *----------------------------------------------------------------------*/
void PARTICLE::ScatraParticleCoupling::SetupGhosting(
    Teuchos::RCP<Epetra_Map> binrowmap, std::map<int, std::set<int>>& scatraelesinbins)
{
  //--------------------------------------------------------------------
  // 1st step: no ghosting for bins and particles
  //--------------------------------------------------------------------

  // Note: IndependentDofSet is used; new dofs are numbered from zero, minnodgid is ignored and it
  // does not register in static_dofsets_
  BinStrategy()->BinDiscret()->FillComplete(true, false, true);

  // set bin row map to be equal to bin col map because we do not have ghosting
  BinColMap() = binrowmap;

  //--------------------------------------------------------------------
  // 2nd step: extend ghosting of underlying scatra discretization according to bin distribution
  //--------------------------------------------------------------------
  std::map<int, std::set<int>> extendedscatraghosting;
  {
    // do communication to gather all elements for extended ghosting
    const int numproc = scatradis_->Comm().NumProc();

    for (int iproc = 0; iproc < numproc; ++iproc)
    {
      // first: proc i tells all procs how many row bins it has
      int numbin = binrowmap->NumMyElements();
      scatradis_->Comm().Broadcast(&numbin, 1, iproc);
      // second: proc i tells all procs which row bins it has
      std::vector<int> binid(numbin, 0);
      if (iproc == MyRank())
      {
        int* binrowmapentries = binrowmap->MyGlobalElements();
        for (int i = 0; i < numbin; ++i) binid[i] = binrowmapentries[i];
      }
      scatradis_->Comm().Broadcast(&binid[0], numbin, iproc);

      // loop over all own bins and find requested ones
      std::map<int, std::set<int>> sdata;
      std::map<int, std::set<int>> rdata;

      for (int i = 0; i < numbin; ++i)
      {
        sdata[binid[i]].insert(
            scatraelesinbins[binid[i]].begin(), scatraelesinbins[binid[i]].end());
      }

      LINALG::Gather<int>(sdata, rdata, 1, &iproc, scatradis_->Comm());

      // proc i has to store the received data
      if (iproc == MyRank())
      {
        extendedscatraghosting = rdata;
      }
    }

    // reduce map of sets to one set and copy to a vector to create extended scatra colmap
    std::set<int> reduscatraeleset;
    std::map<int, std::set<int>>::iterator iter;
    for (iter = extendedscatraghosting.begin(); iter != extendedscatraghosting.end(); ++iter)
    {
      reduscatraeleset.insert(iter->second.begin(), iter->second.end());
    }

    // insert existing ghosting of scatra
    const Epetra_Map* elecolmap = scatradis_->ElementColMap();
    for (int lid = 0; lid < elecolmap->NumMyElements(); ++lid)
      reduscatraeleset.insert(elecolmap->GID(lid));

    std::vector<int> scatracolgids(reduscatraeleset.begin(), reduscatraeleset.end());
    Teuchos::RCP<Epetra_Map> scatracolmap =
        Teuchos::rcp(new Epetra_Map(-1, (int)scatracolgids.size(), &scatracolgids[0], 0, Comm()));

    scatradis_->ExtendedGhosting(*scatracolmap, true, true, true, false);
  }

  //--------------------------------------------------------------------
  // 4th step: assign scatra elements to bins
  //--------------------------------------------------------------------
  {
    for (std::map<int, std::set<int>>::const_iterator biniter = extendedscatraghosting.begin();
         biniter != extendedscatraghosting.end(); ++biniter)
    {
      DRT::MESHFREE::MeshfreeMultiBin* currbin = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(
          BinStrategy()->BinDiscret()->gElement(biniter->first));
      for (std::set<int>::const_iterator scatraeleiter = biniter->second.begin();
           scatraeleiter != biniter->second.end(); ++scatraeleiter)
      {
        int scatraeleid = *scatraeleiter;
        currbin->AddAssociatedEle(
            bin_transpcontent_, scatraeleid, scatradis_->gElement(scatraeleid));
        //          cout << "in bin with id:" << currbin->Id() << " is fluid ele with id" <<
        //          fluideleid << "with pointer" << fluiddis_->gElement(fluideleid) << endl;
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
  const int numcolbin = BinStrategy()->BinDiscret()->NumMyColElements();
  for (int ibin = 0; ibin < numcolbin; ++ibin)
  {
    DRT::Element* actele = BinStrategy()->BinDiscret()->lColElement(ibin);
    DRT::MESHFREE::MeshfreeMultiBin* actbin =
        dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(actele);
    const int numfluidele = actbin->NumAssociatedEle(bin_transpcontent_);
    const int* fluideleids = actbin->AssociatedEleIds(bin_transpcontent_);
    std::vector<DRT::Element*> fluidelements(numfluidele);
    for (int iele = 0; iele < numfluidele; ++iele)
    {
      const int fluideleid = fluideleids[iele];
      fluidelements[iele] = scatradis_->gElement(fluideleid);
    }
    actbin->BuildElePointers(bin_transpcontent_, &fluidelements[0]);
  }

  return;
}

/*----------------------------------------------------------------------*
| delete particles out of physical domain               rasthofer 11/14 |
 *----------------------------------------------------------------------*/
void PARTICLE::ScatraParticleCoupling::DeleteParticlesOutOfPhysicalDomain()
{
  // delete particles that are out of the physical domain
  // possible situation: there is a hole somewhere
  // or the obvious extension of a row of bins around the domain is considered

  // vector containing ids of particles to be deleted
  std::vector<int> deleteparticles;

  // loop all particles
  for (int inode = 0; inode < BinStrategy()->BinDiscret()->NumMyRowNodes(); inode++)
  {
    // get global id of current particle
    const int gid = (BinStrategy()->BinDiscret()->NodeRowMap())->GID(inode);

    // element coordinates of particle position in scatra element
    LINALG::Matrix<3, 1> elecoord(true);
    DRT::Element* scatraele = GetEleCoordinates(gid, elecoord);
    // particle is in hole of domain, don't move it any more
    // this results in zero velocity, since vel is initialized with zeros
    if (scatraele == NULL) deleteparticles.push_back(gid);
  }

  if (deleteparticles.size() > 0)
  {
    dserror(
        "There seems to be a particle outside of the physical domain! Check situation and remove "
        "this dserror if everything seems correct!");
    // remove particles (note: calls FillComplete)
    DeleteParticles(deleteparticles);
  }

  return;
}
