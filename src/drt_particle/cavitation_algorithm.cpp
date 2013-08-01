/*----------------------------------------------------------------------*/
/*!
\file cavitation_algorithm.cpp

\brief Algorithm to control cavitation simulations

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 11/12 |
 *----------------------------------------------------------------------*/
#include "cavitation_algorithm.H"
#include "particle_timint_centrdiff.H"
#include "../drt_adapter/ad_fld_base_algorithm.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_fluid_ele/fluid_ele_calc.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/particle_mat.H"
#include "../drt_meshfree_discret/drt_meshfree_multibin.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_parallel.H"

#include "../drt_geometry/searchtree.H"
#include "../drt_geometry/searchtree_geometry_service.H"
#include "../drt_geometry/intersection_math.H"
#include "../drt_geometry/element_coordtrafo.H"
#include "../drt_geometry/position_array.H"
#include "../drt_cut/cut_position.H"
#include "../linalg/linalg_utils.H"

#include "../drt_inpar/inpar_cavitation.H"
#include "../drt_io/io_pstream.H"
#include "../headers/definitions.h"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

/*----------------------------------------------------------------------*
 | Algorithm constructor                                    ghamm 11/12 |
 *----------------------------------------------------------------------*/
CAVITATION::Algorithm::Algorithm(
  const Epetra_Comm& comm,
  const Teuchos::ParameterList& params
  ) : PARTICLE::Algorithm(comm,params),
  coupalgo_(DRT::INPUT::IntegralValue<INPAR::CAVITATION::CouplingStrategyOverFields>(params,"COUPALGO")),
  fluiddis_(Teuchos::null)
{
  // setup fluid time integrator
  fluiddis_ = DRT::Problem::Instance()->GetDis("fluid");
  // ask base algorithm for the fluid time integrator
  Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluid =
      Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(DRT::Problem::Instance()->CavitationParams(),false));
  fluid_ = fluid->FluidFieldrcp();

  return;
}


/*----------------------------------------------------------------------*
 | time loop of the cavitation algorithm                    ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Timeloop()
{
  // time loop
  while (NotFinished())
  {
    // counter and print header; predict solution of both fields
    PrepareTimeStep();

    // particle time step is solved
    Integrate();

    // deal with particle inflow
    ParticleInflow();

    // transfer particles into their correct bins
    TransferParticles();

    // update displacements, velocities, accelerations
    // after this call we will have disn_==dis_, etc
    // update time and step
    Update();

    // write output to screen and files
    Output();

  }  // NotFinished

}


/*----------------------------------------------------------------------*
 | setup of the system                                      ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::SetupSystem()
{
  return;
}


/*----------------------------------------------------------------------*
 | initialization of the system                             ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Init(bool restarted)
{
  // FillComplete() necessary for DRT::Geometry .... could be removed perhaps
  particledis_->FillComplete(false,false,false);
  // extract noderowmap because it will be called Reset() after adding elements
  Teuchos::RCP<Epetra_Map> particlerowmap = Teuchos::rcp(new Epetra_Map(*particledis_->NodeRowMap()));
  CreateBins();

  // gather all fluid coleles in each bin for proper extended ghosting
  std::map<int, std::set<int> > fluideles;
  Teuchos::RCP<Epetra_Map> binrowmap = DistributeBinsToProcs(fluideles);

  // the following has only to be done once --> skip in case of restart
  if(not restarted)
  {
    // read out bubble inflow condition and set bubble inflows in corresponding bins
    // assumption: only row bins are available up to here
    BuildBubbleInflowCondition();
  }

  //--------------------------------------------------------------------
  // -> 1) create a set of homeless particles that are not in a bin on this proc
  std::set<Teuchos::RCP<DRT::Node>,PARTICLE::Less> homelessparticles;

  for (int lid = 0; lid < particlerowmap->NumMyElements(); ++lid)
  {
    DRT::Node* node = particledis_->gNode(particlerowmap->GID(lid));
    const double* currpos = node->X();
    PlaceNodeCorrectly(Teuchos::rcp(node,false), currpos, homelessparticles);
  }

  // start round robin loop to fill particles into their correct bins
  FillParticlesIntoBins(homelessparticles);

  // ghost bins, particles and fluid elements according to the bins
  SetupGhosting(binrowmap, fluideles);

  // the following has only to be done once --> skip in case of restart
  if(not restarted)
  {
    // access structural dynamic params list which will be possibly modified while creating the time integrator
    const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

    // create time integrator based on structural time integration
    Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> particles =
        Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(DRT::Problem::Instance()->CavitationParams(), const_cast<Teuchos::ParameterList&>(sdyn), particledis_));
    particles_ = particles->StructureFieldrcp();

    // determine consistent initial acceleration for the particles
    CalculateAndApplyForcesToParticles();
    Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->DetermineMassDampConsistAccel();
  }

  // some output
  IO::cout << "after ghosting" << IO::endl;
  DRT::UTILS::PrintParallelDistribution(*particledis_);
  DRT::UTILS::PrintParallelDistribution(*fluiddis_);

  return;
}


/*----------------------------------------------------------------------*
 | prepare time step                                       ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::PrepareTimeStep()
{
  PARTICLE::Algorithm::PrepareTimeStep();
  fluid_->PrepareTimeStep();
  return;
}


/*----------------------------------------------------------------------*
 | solve the current particle time step                    ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Integrate()
{
  // apply forces and solve particle time step
  PARTICLE::Algorithm::Integrate();

  {
    Teuchos::RCP<Teuchos::Time> t = Teuchos::TimeMonitor::getNewTimer("CAVITATION::Algorithm::IntegrateFluid");
    Teuchos::TimeMonitor monitor(*t);
    fluid_->Solve();
  }

  return;
}


/*----------------------------------------------------------------------*
 | calculate fluid forces on particle and apply it         ghamm 01/13  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::CalculateAndApplyForcesToParticles()
{
  Teuchos::RCP<Teuchos::Time> t = Teuchos::TimeMonitor::getNewTimer("CAVITATION::Algorithm::CalculateAndApplyForcesToParticles");
  Teuchos::TimeMonitor monitor(*t);
  const int dim = 3;

  fluiddis_->ClearState();
  particledis_->ClearState();

  //TODO 1: Test whether these states are updated enough
  // at the beginning of the time step: veln := velnp(previous step)
  fluiddis_->SetState("veln",fluid_->Veln());
  fluiddis_->SetState("velnm",fluid_->Velnm());

  particledis_->SetState("bubblepos", particles_->ExtractDispn());
  particledis_->SetState("bubblevel", particles_->ExtractVeln());
  particledis_->SetState("bubbleacc", particles_->ExtractAccn());
  // miraculous transformation from row to col layout ...
  Teuchos::RCP<const Epetra_Vector> bubblepos = particledis_->GetState("bubblepos");
  Teuchos::RCP<const Epetra_Vector> bubblevel = particledis_->GetState("bubblevel");
  Teuchos::RCP<const Epetra_Vector> bubbleacc = particledis_->GetState("bubbleacc");

  // bubble radius of layout nodal col map --> SetState not possible
  Teuchos::RCP<Epetra_Vector> bubbleradius = LINALG::CreateVector(*particledis_->NodeColMap(),false);
  LINALG::Export(*Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->ExtractRadiusnp(),*bubbleradius);

  // vectors to be filled with forces
  Teuchos::RCP<Epetra_Vector> bubbleforces = LINALG::CreateVector(*particledis_->DofColMap(),true);
  Teuchos::RCP<Epetra_Vector> fluidforces = LINALG::CreateVector(*fluiddis_->DofRowMap(),true);

  if(fluiddis_->NumMyColElements() <= 0)
    dserror("there is no fluid element to ask for material parameters");
  Teuchos::RCP<MAT::NewtonianFluid> actmat = Teuchos::rcp_dynamic_cast<MAT::NewtonianFluid>(fluiddis_->lColElement(0)->Material());
  if(actmat == Teuchos::null)
    dserror("type cast of fluid material failed");

  double rho_l = actmat->Density();
  double mu_l = actmat->Viscosity();
  double rho_b = Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->ParticleDensity();

  // define element matrices and vectors
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector1;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;

  // all particles (incl ghost) are evaluated
  for(int i=0;i<particledis_->NodeColMap()->NumMyElements();++i)
  {
    DRT::Node* currparticle = particledis_->lColNode(i);
    // fill particle position
    LINALG::Matrix<3,1> particleposition;
    std::vector<int> lm_b = particledis_->Dof(currparticle);
    int posx = bubblepos->Map().LID(lm_b[0]);
    for (int d=0; d<dim; ++d)
      particleposition(d) = (*bubblepos)[posx+d];


    //--------------------------------------------------------------------
    // 1st step: element coordinates of particle position in fluid element
    //--------------------------------------------------------------------

    // variables to store information about element in which the particle is located
    DRT::Element* targetfluidele = NULL;
    LINALG::Matrix<3,1> elecoord(true);

    // find out in which fluid element the current particle is located
    if(currparticle->NumElement() != 1)
      dserror("ERROR: A particle is assigned to more than one bin!");
    DRT::Element** currele = currparticle->Elements();
#ifdef DEBUG
    DRT::MESHFREE::MeshfreeMultiBin* test = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currele[0]);
    if(test == NULL) dserror("dynamic cast from DRT::Element to DRT::MESHFREE::MeshfreeMultiBin failed");
#endif
    DRT::MESHFREE::MeshfreeMultiBin* currbin = static_cast<DRT::MESHFREE::MeshfreeMultiBin*>(currele[0]);
    DRT::Element** fluidelesinbin = currbin->AssociatedFluidEles();
    int numfluidelesinbin = currbin->NumAssociatedFluidEle();

    std::set<int>::const_iterator eleiter;
    for(int ele=0; ele<numfluidelesinbin; ++ele)
    {
      DRT::Element* fluidele = fluidelesinbin[ele];
      DRT::Node** fluidnodes = fluidele->Nodes();
      const int numnode = fluidele->NumNode();

      // fill currentpositions of fluid element
      std::map<int,LINALG::Matrix<3,1> > currentfluidpositions;
      for (int j=0; j < numnode; j++)
      {
        const DRT::Node* node = fluidnodes[j];
        LINALG::Matrix<3,1> currpos;
        for (int d=0; d<dim; ++d)
        {
          currpos(d,0) = node->X()[d];
        }
        currentfluidpositions[node->Id()] = currpos;
      }
      const LINALG::SerialDenseMatrix xyze(GEO::getCurrentNodalPositions(fluidele, currentfluidpositions));

      //get coordinates of the particle position in parameter space of the element
      bool insideele = GEO::currentToVolumeElementCoordinates(fluidele->Shape(), xyze, particleposition, elecoord);

      if(insideele == true)
      {
        targetfluidele = fluidele;
        // leave loop over all fluid eles in bin
        break;
      }
    }


    //--------------------------------------------------------------------
    // 2nd step: forces on this bubble are calculated
    //--------------------------------------------------------------------

    if(targetfluidele == NULL)
    {
      std::cout << "INFO: currparticle with Id: " << currparticle->Id() << " and position: " << particleposition(0) << " "
          << particleposition(1) << " " << particleposition(2) << " " << " does not have an underlying fluid element -> no forces calculated" << std::endl;

      std::vector<double> tmpposition(dim);
      for(int d=0; d<dim; ++d)
        tmpposition[d] = particleposition(d);
      int bubbleBinId = ConvertPosToGid(tmpposition);
      cout << "particle is in binId: " << bubbleBinId << " while currbin->Id() is " << currbin->Id() <<
          " . The following number of fluid eles is in this bin:" << numfluidelesinbin << endl;

      // do not assemble forces for this bubble and continue with next bubble
      continue;
    }

    // get element location vector and ownerships
    std::vector<int> lm_f;
    std::vector<int> lmowner_f;
    std::vector<int> lmstride;
    targetfluidele->LocationVector(*fluiddis_,lm_f,lmowner_f,lmstride);

    // Reshape element matrices and vectors and initialize to zero
    elevector1.Size(dim);
    elevector2.Size(dim);
    elevector3.Size(dim);

    // set action in order to calculate the velocity and material derivative of the velocity
    Teuchos::ParameterList params;
    params.set<int>("action",FLD::calc_mat_deriv_u_and_rot_u);
    params.set<double>("timestep",Dt());
    params.set<LINALG::Matrix<3,1> >("elecoords", elecoord);

    // call the element specific evaluate method (elevec1 = fluid vel u; elevec2 = mat deriv of fluid vel, elevec3 = rot of fluid vel)
    targetfluidele->Evaluate(params,*fluiddis_,lm_f,elematrix1,elematrix2,elevector1,elevector2,elevector3);

    // get bubble velocity and acceleration
    std::vector<double> v_bub(lm_b.size());
    DRT::UTILS::ExtractMyValues(*bubblevel,v_bub,lm_b);

    // get bubble radius
    double r_bub = (*bubbleradius)[ particledis_->NodeColMap()->LID(currparticle->Id()) ];

    // bubble Reynolds number
    LINALG::Matrix<3,1> v_rel(false);
    for (int d=0; d<dim; ++d)
      v_rel(d) = elevector1[d] - v_bub[d];

    double v_relabs = v_rel.Norm2();
    double Re_b = 2.0 * r_bub * v_relabs * rho_l / mu_l;

    bool output = false;
    if(output)
    {
      cout << "v_rel: " << v_rel(0) << "  " << v_rel(1) << "  " << v_rel(2) << "  " << endl;
      cout << "v_relabs: " << v_relabs << endl;
      cout << "bubble Reynolds number: " << Re_b << endl;
    }

    // variable to sum forces for the current bubble under observation
    LINALG::Matrix<3,1> sumforces;
    /*------------------------------------------------------------------*/
    //// 2.1) drag force = 0.5 * c_d * rho_l * Pi * r_b^2 * |u-v| * (u-v) or
    //// Stokes law for very small Re: drag force = 6.0 * Pi * mu_l * r_b * (u-v)
    double coeff1 = 0.0;
    if(Re_b < 0.1)
    {
      coeff1 = 6.0 * M_PI * mu_l * r_bub;
    }
    else
    {
      double c_d = 0.0;
      if(Re_b < 1000.0)
        c_d = 24.0 * (1.0 + 0.15 * pow(Re_b,0.687)) / Re_b;
      else
        c_d = 0.44;

      coeff1 = 0.5 * c_d * rho_l * M_PI * r_bub * r_bub * v_relabs;
    }

    LINALG::Matrix<3,1> dragforce(v_rel);
    dragforce.Scale(coeff1);
    //assemble
    sumforces.Update(1.0, dragforce);
    /*------------------------------------------------------------------*/

    /*------------------------------------------------------------------*/
    //// 2.2) lift force = c_l * rho_l * volume_b * (u-v) x rot_u   with rot_u = nabla x u
    double c_l = 0.5;
    double vol_b = 4.0 / 3.0 * M_PI * r_bub * r_bub* r_bub;
    LINALG::Matrix<3,1> rot_u;
    for (int d=0; d<dim; ++d)
      rot_u(d) = elevector3(d);

    LINALG::Matrix<3,1> liftforce;
    liftforce(0) = v_rel(1) * rot_u(2) - v_rel(2) * rot_u(1);
    liftforce(1) = v_rel(2) * rot_u(0) - v_rel(0) * rot_u(2);
    liftforce(2) = v_rel(0) * rot_u(1) - v_rel(1) * rot_u(0);

    double coeff2 = c_l * rho_l * vol_b;
    liftforce.Scale(coeff2);
    //assemble
    sumforces.Update(1.0, liftforce, 1.0);
    // store forces for coupling to fluid
    LINALG::Matrix<3,1> couplingforce(sumforces);
    /*------------------------------------------------------------------*/

    /*------------------------------------------------------------------*/
    //// 2.3) gravity and buoyancy forces = volume_b * rho_bub * g - volume_b * rho_l * ( g - Du/Dt )
    LINALG::Matrix<3,1> Du_Dt;
    for (int d=0; d<dim; ++d)
      Du_Dt(d) = elevector2[d];

    LINALG::Matrix<3,1> grav_buoy_force;
    grav_buoy_force.Update(rho_b, gravity_acc_);
    grav_buoy_force.Update(-rho_l, gravity_acc_, rho_l, Du_Dt, 1.0);
    grav_buoy_force.Scale(vol_b);
    //assemble
    sumforces.Update(1.0, grav_buoy_force, 1.0);
    /*------------------------------------------------------------------*/

    /*------------------------------------------------------------------*/
    //// 2.4) virtual/added mass = c_VM * rho_l * volume_b * ( Du/Dt - Dv/Dt )
    //// Note: implicit treatment of bubble acceleration in added mass, other forces explicit
    //// final force = \frac{ sum all forces (2.1, 2.2, 2.3) + c_VM * rho_l * volume_b * Du/Dt }{ 1 + c_VM * rho_l / rho_b }
    double c_VM = 0.5;
    double coeff3 = c_VM * rho_l * vol_b;
    double coeff4 = 1.0 + c_VM * rho_l / rho_b;

    LINALG::Matrix<3,1> bubbleforce;
    bubbleforce.Update(1.0/coeff4, sumforces, coeff3/coeff4, Du_Dt);
    /*------------------------------------------------------------------*/


    //--------------------------------------------------------------------
    // 3rd step: assemble bubble forces
    //--------------------------------------------------------------------

    // assemble of bubble forces can be done without further need of LINALG::Export (ghost bubbles also evaluated)
    Epetra_SerialDenseVector forcecurrbubble(3);
    for(int d=0; d<dim; ++d)
      forcecurrbubble[d] = bubbleforce(d);
    // in order to assemble all entries even for col layout, it is necessary to use myrank_ here
    std::vector<int> lmowner_b(lm_b.size(), myrank_);
    LINALG::Assemble(*bubbleforces,forcecurrbubble,lm_b,lmowner_b);

    // coupling forces between fluid and particle only include certain forces
    // calculate added mass force
    LINALG::Matrix<3,1> addedmassforce;
    double m_b = vol_b * rho_b;
    addedmassforce.Update(coeff3, Du_Dt, -coeff3/m_b, bubbleforce);
    //// coupling force = -(dragforce + liftforce + addedmassforce); actio = reactio --> minus sign
    couplingforce.Update(-1.0, addedmassforce, -1.0);

    // assemble of fluid forces can only be done on row elements (all bubbles in this element must be considered: ghost near boundary)
    const int numnode = targetfluidele->NumNode();
    Epetra_SerialDenseVector funct(numnode);
    // get shape functions of the element; evaluated at the bubble position --> distribution
    DRT::UTILS::shape_function_3D(funct,elecoord(0,0),elecoord(1,0),elecoord(2,0),targetfluidele->Shape());
    // prepare assembly for fluid forces (pressure degrees do not have to be filled)
    Epetra_SerialDenseVector val(numnode*(dim+1));
    for(int iter=0; iter<numnode; iter++)
    {
      for(int d=0; d<dim; ++d)
      {
        val[iter*(dim+1) + d] = funct[iter] * couplingforce(d);
      }
    }
    // do assembly of bubble forces on fluid
    LINALG::Assemble(*fluidforces,val,lm_f,lmowner_f);


    //--------------------------------------------------------------------
    // 4th step: output
    //--------------------------------------------------------------------
    if(output)
    {
      // gravity
      LINALG::Matrix<3,1> gravityforce(gravity_acc_);
      gravityforce.Scale(rho_b*vol_b);
      cout << "gravity force       : " << gravityforce << endl;

      // buoyancy
      double coeff5 = - vol_b * rho_l;
      LINALG::Matrix<3,1> buoyancyforce(true);
      buoyancyforce.Update(coeff5, gravity_acc_);
      cout << "buoyancy force      : " << buoyancyforce << endl;

      // effective buoyancy / inertia term
      LINALG::Matrix<3,1> effectbuoyancyforce;
      effectbuoyancyforce.Update(-coeff5, Du_Dt);
      cout << "effective buoy force: " << effectbuoyancyforce << endl;

      // drag, lift and added mass force
      cout << "dragforce force     : " << dragforce << endl;
      cout << "liftforce force     : " << liftforce << endl;
      cout << "added mass force    : " << addedmassforce << endl;

      // sum over all bubble forces
      cout << "sum over all forces : " << bubbleforce << endl;

      // fluid force
      cout << "fluid force         : " << couplingforce << endl;
    }

  } // end iparticle

  //--------------------------------------------------------------------
  // 4th step: apply forces to bubbles and fluid field
  //--------------------------------------------------------------------
  particledis_->SetState("particleforces", bubbleforces);

  if(coupalgo_ != INPAR::CAVITATION::OneWay)
    fluid_->ApplyExternalForces(fluidforces);

  return;
}


/*----------------------------------------------------------------------*
 | particles are inserted into domain                      ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::ParticleInflow()
{
  std::map<int, std::list<Teuchos::RCP<BubbleSource> > >::const_iterator biniter;

  int timeforinflow = 0;
  for(biniter=bubble_source_.begin(); biniter!=bubble_source_.end(); ++biniter)
  {
    // all particles have the same inflow frequency --> it is enough to test one
    // assumption only valid in case of one condition or conditions with identical inflow frequency
    if(biniter->second.size() != 0)
    {
      double inflowtime = 1.0 / biniter->second.front()->inflow_freq_;
      if(Step() % ((int)(inflowtime/Dt())) == 0)
      {
        timeforinflow = 1;
        break;
      }
    }
  }

  int globaltimeforinflow = 0;
  particledis_->Comm().MaxAll(&timeforinflow, &globaltimeforinflow, 1);
  if(globaltimeforinflow == 0)
    return; // no inflow detected


  // initialize bubble id with largest bubble id in use + 1 (on each proc)
  int maxbubbleid = particledis_->NodeRowMap()->MaxAllGID()+1;

  // start filling particles
  int inflowcounter = 0;
  for(biniter=bubble_source_.begin(); biniter!=bubble_source_.end(); ++biniter)
  {
    std::list<Teuchos::RCP<BubbleSource> >::const_iterator particleiter;
    for(particleiter=biniter->second.begin(); particleiter!=biniter->second.end(); ++particleiter)
    {
      std::vector<double> inflow_position = (*particleiter)->inflow_position_;
      std::set<Teuchos::RCP<DRT::Node>,PARTICLE::Less> homelessparticles;
      int newbubbleid = maxbubbleid + (*particleiter)->inflowid_;
      Teuchos::RCP<DRT::Node> newparticle = Teuchos::rcp(new DRT::Node(newbubbleid, &inflow_position[0], myrank_));
      PlaceNodeCorrectly(newparticle, &inflow_position[0], homelessparticles);
      if(homelessparticles.size() != 0)
        dserror("New bubble could not be inserted on this proc! Bubble inflow broken.");
    }
    inflowcounter += (int)biniter->second.size();
  }

  std::cout << "Inflow of " << inflowcounter << " bubbles on proc " << myrank_ << std::endl;

  // rebuild connectivity and assign degrees of freedom (note: IndependentDofSet)
  particledis_->FillComplete(true, false, true);

  // update of state vectors to the new maps
  Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->UpdateStatesAfterParticleTransfer();

  // insert data for new bubbles into state vectors
  const Epetra_Map* dofrowmap = particledis_->DofRowMap();
  const Epetra_Map* noderowmap = particledis_->NodeRowMap();
  Teuchos::RCP<Epetra_Vector> disn = particles_->ExtractDispnp();
  Teuchos::RCP<Epetra_Vector> veln = particles_->ExtractVelnp();
  Teuchos::RCP<Epetra_Vector> radiusn = Teuchos::rcp_dynamic_cast<PARTICLE::TimIntCentrDiff>(particles_)->ExtractRadiusnp();

  for(biniter=bubble_source_.begin(); biniter!=bubble_source_.end(); ++biniter)
  {
    std::list<Teuchos::RCP<BubbleSource> >::const_iterator particleiter;
    for(particleiter=biniter->second.begin(); particleiter!=biniter->second.end(); ++particleiter)
    {
      std::vector<double> inflow_position = (*particleiter)->inflow_position_;
      std::vector<double> inflow_vel = (*particleiter)->inflow_vel_;
      int inflow_vel_curve = (*particleiter)->inflow_vel_curve_;
      double inflow_radius = (*particleiter)->inflow_radius_;
      int newbubbleid = maxbubbleid + (*particleiter)->inflowid_;

      double curvefac = 1.0;
      // curves are numbered starting with 1 in the input file
      if(inflow_vel_curve > 0)
        curvefac = DRT::Problem::Instance()->Curve(inflow_vel_curve-1).f(Time());

      DRT::Node* currparticle = particledis_->gNode(newbubbleid);
      // get the first gid of a particle and convert it into a LID
      int lid = dofrowmap->LID(particledis_->Dof(currparticle, 0));
      for(int dim=0; dim<3; dim++)
      {
        (*disn)[lid+dim] = inflow_position[dim];
        (*veln)[lid+dim] = inflow_vel[dim] * curvefac;
      }
      lid = noderowmap->LID(newbubbleid);
      (*radiusn)[lid] = inflow_radius;
    }
  }

  return;
}
/*----------------------------------------------------------------------*
 | update the current time step                            ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Update()
{
  PARTICLE::Algorithm::Update();
  fluid_->Update();
  return;
}


/*----------------------------------------------------------------------*
| read restart information for given time step              ghamm 03/13 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::ReadRestart(int restart)
{
  PARTICLE::Algorithm::ReadRestart(restart);
  fluid_->ReadRestart(restart);
  return;
}


/*----------------------------------------------------------------------*
| find XAABB and divide into bins                           ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::CreateBins()
{
  // if not yet specified, get XAABB_ from underlying discretization
  if( XAABB_(2,1) > 0.9e12  and  XAABB_(2,1) < 1.1e12 )
  {
    IO::cout << "XAABB is computed based on the underlying fluid discretization" << IO::endl;
    XAABB_ = GEO::getXAABBofDis(*fluiddis_);
    // local bounding box
    double locmin[3] = {XAABB_(0,0), XAABB_(1,0), XAABB_(2,0)};
    double locmax[3] = {XAABB_(0,1), XAABB_(1,1), XAABB_(2,1)};
    // global bounding box
    double globmin[3];
    double globmax[3];
    // do the necessary communication
    Comm().MinAll(&locmin[0], &globmin[0], 3);
    Comm().MaxAll(&locmax[0], &globmax[0], 3);

    for(int dim=0; dim<3; dim++)
    {
      XAABB_(dim,0) = globmin[dim];
      XAABB_(dim,1) = globmax[dim];
    }
  }

  // divide global bounding box into bins
  for (int dim = 0; dim < 3; dim++)
  {
    // std::floor leads to bins that are at least of size cutoff_radius
    bin_per_dir_[dim] = std::max(1, (int)((XAABB_(dim,1)-XAABB_(dim,0))/cutoff_radius_));
    bin_size_[dim] = (XAABB_(dim,1)-XAABB_(dim,0))/bin_per_dir_[dim];
  }

  IO::cout << "Global bounding box size: " << XAABB_;
  IO::cout << "bins per direction: " << "x = " << bin_per_dir_[0] << " y = " << bin_per_dir_[1] << " z = " << bin_per_dir_[2] << IO::endl;

  return;
}


/*----------------------------------------------------------------------*
| bins are distributed to the processors                    ghamm 11/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> CAVITATION::Algorithm::DistributeBinsToProcs(std::map<int, std::set<int> >& fluideles)
{
  //--------------------------------------------------------------------
  // 1st step: exploiting bounding box idea for fluid elements and bins
  //--------------------------------------------------------------------
  {
    Teuchos::RCP<Teuchos::Time> t = Teuchos::TimeMonitor::getNewTimer("CAVITATION::Algorithm::DistributeBinsToProcs_step1");
    Teuchos::TimeMonitor monitor(*t);
    // loop over row fluid elements is enough, extended ghosting always includes standard ghosting
    for (int lid = 0; lid < fluiddis_->NumMyRowElements(); ++lid)
    {
      DRT::Element* fluidele = fluiddis_->lRowElement(lid);
      DRT::Node** fluidnodes = fluidele->Nodes();
      const int numnode = fluidele->NumNode();

      // initialize ijk_range with ijk of first node of fluid element
      int ijk[3];
      {
        const DRT::Node* node = fluidnodes[0];
        const double* coords = node->X();
        ConvertPosToijk(coords, ijk);
      }

      // ijk_range contains: i_min i_max j_min j_max k_min k_max
      int ijk_range[] = {ijk[0], ijk[0], ijk[1], ijk[1], ijk[2], ijk[2]};

      // fill in remaining nodes
      for (int j=1; j<numnode; j++)
      {
        const DRT::Node* node = fluidnodes[j];
        const double* coords = node->X();
        int ijk[3];
        ConvertPosToijk(coords, ijk);

        for(int dim=0; dim<3; dim++)
        {
          if(ijk[dim]<ijk_range[dim*2])
            ijk_range[dim*2]=ijk[dim];
          if(ijk[dim]>ijk_range[dim*2+1])
            ijk_range[dim*2+1]=ijk[dim];
        }
      }

      // get corresponding bin ids in ijk range
      std::set<int> binIds;
      GidsInijkRange(&ijk_range[0], binIds, false);

      // assign fluid element to bins
      for(std::set<int>::const_iterator biniter=binIds.begin(); biniter!=binIds.end(); ++biniter)
        fluideles[*biniter].insert(fluidele->Id());
    }
  }


  //--------------------------------------------------------------------
  // 2nd step: decide which proc will be owner of each bin
  //--------------------------------------------------------------------

  std::vector<int> rowbins;
  {
    Teuchos::RCP<Teuchos::Time> t = Teuchos::TimeMonitor::getNewTimer("CAVITATION::Algorithm::DistributeBinsToProcs_step2");
    Teuchos::TimeMonitor monitor(*t);
    // NOTE: This part of the setup can be the bottleneck because vectors of all bins
    // are needed on each proc (memory issue!!); std::map could perhaps help when gathering
    // num fluid nodes in each bin, then block wise communication after copying data to vector

    int numbins = bin_per_dir_[0]*bin_per_dir_[1]*bin_per_dir_[2];
    std::vector<int> mynumeles_per_bin(numbins,0);

    std::map<int, std::set<int> >::const_iterator iter;
    for(iter=fluideles.begin(); iter!=fluideles.end(); ++ iter)
    {
      mynumeles_per_bin[iter->first] = iter->second.size();
    }

    // find maximum number of eles in each bin over all procs (init with -1)
    std::vector<int> maxnumeles_per_bin(numbins,-1);
    fluiddis_->Comm().MaxAll(&mynumeles_per_bin[0], &maxnumeles_per_bin[0], numbins);

    // it is possible that several procs have the same number of eles in a bin
    // only proc which has maximum number of eles in a bin writes its rank
    std::vector<int> myrank_per_bin(numbins,-1);
    for(int i=0; i<numbins; i++)
    {
      if(mynumeles_per_bin[i] == maxnumeles_per_bin[i])
        myrank_per_bin[i] = myrank_;
    }

    mynumeles_per_bin.clear();
    maxnumeles_per_bin.clear();

    // find maximum myrank for each bin over all procs (init with -1)
    std::vector<int> maxmyrank_per_bin(numbins,-1);
    fluiddis_->Comm().MaxAll(&myrank_per_bin[0], &maxmyrank_per_bin[0], numbins);

    // distribute bins to proc with highest rank
    for(int gid=0; gid<numbins; gid++)
    {
      if(myrank_ == maxmyrank_per_bin[gid])
      {
        Teuchos::RCP<DRT::Element> bin = DRT::UTILS::Factory("MESHFREEMULTIBIN","dummy", gid, myrank_);
        particledis_->AddElement(bin);
        rowbins.push_back(gid);
      }
    }

    myrank_per_bin.clear();
    maxmyrank_per_bin.clear();
  }

  // return binrowmap (without having called FillComplete on particledis_ so far)
  return Teuchos::rcp(new Epetra_Map(-1,(int)rowbins.size(),&rowbins[0],0,Comm()));
}


/*----------------------------------------------------------------------*
| setup ghosting of bins, particles & underlying fluid      ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::SetupGhosting(Teuchos::RCP<Epetra_Map> binrowmap, std::map<int, std::set<int> >& fluideles)
{
  //--------------------------------------------------------------------
  // 1st and 2nd step
  //--------------------------------------------------------------------

  PARTICLE::Algorithm::SetupGhosting(binrowmap);


  //--------------------------------------------------------------------
  // 3st step: extend ghosting of underlying fluid discretization according to bin distribution
  //--------------------------------------------------------------------
  std::map<int, std::set<int> > extendedfluidghosting;
  {
    // do communication to gather all elements for extended ghosting
    const int numproc = fluiddis_->Comm().NumProc();

    for (int iproc = 0; iproc < numproc; ++iproc)
    {
      // first: proc i tells all procs how many col bins it has
      int numbin = bincolmap_->NumMyElements();
      fluiddis_->Comm().Broadcast(&numbin, 1, iproc);
      // second: proc i tells all procs which col bins it has
      std::vector<int> binid(numbin,0);
      if(iproc == myrank_)
      {
        int* bincolmap = bincolmap_->MyGlobalElements();
        for (int i=0; i<numbin; ++i)
          binid[i] = bincolmap[i];
      }
      fluiddis_->Comm().Broadcast(&binid[0], numbin, iproc);

      // loop over all own bins and find requested ones
      std::map<int, std::set<int> > sdata;
      std::map<int, std::set<int> > rdata;

      for(int i=0; i<numbin; i++)
      {
        sdata[binid[i]].insert(fluideles[binid[i]].begin(),fluideles[binid[i]].end());
      }

      LINALG::Gather<int>(sdata, rdata, 1, &iproc, fluiddis_->Comm());

      // proc i has to store the received data
      if(iproc == myrank_)
      {
        extendedfluidghosting = rdata;
      }
    }

    //reduce map of sets to one set and copy to a vector to create fluidcolmap
    std::set<int> redufluideleset;
    std::map<int, std::set<int> >::iterator iter;
    for(iter=extendedfluidghosting.begin(); iter!= extendedfluidghosting.end(); ++iter)
    {
      redufluideleset.insert(iter->second.begin(),iter->second.end());
    }
    std::vector<int> fluidcolgids(redufluideleset.begin(),redufluideleset.end());
    Teuchos::RCP<Epetra_Map> fluidcolmap = Teuchos::rcp(new Epetra_Map(-1,(int)fluidcolgids.size(),&fluidcolgids[0],0,Comm()));

    // create ghosting for fluid eles (each knowing its node ids)
    fluiddis_->ExportColumnElements(*fluidcolmap);

    // create a set of node IDs for each proc (row + ghost)
    std::set<int> nodes;
    for (int lid=0;lid<fluidcolmap->NumMyElements();++lid)
    {
      DRT::Element* ele = fluiddis_->gElement(fluidcolmap->GID(lid));
      const int* nodeids = ele->NodeIds();
      for(int inode=0;inode<ele->NumNode();inode++)
        nodes.insert(nodeids[inode]);
    }

    // copy nodegids to a vector and create nodecolmap
    std::vector<int> colnodes(nodes.begin(),nodes.end());
    Teuchos::RCP<Epetra_Map> nodecolmap = Teuchos::rcp(new Epetra_Map(-1,(int)colnodes.size(),&colnodes[0],0,Comm()));

    // create ghosting for nodes
    fluiddis_->ExportColumnNodes(*nodecolmap);

    // do a final fillcomplete to build connectivity
    fluiddis_->FillComplete(true,true,true);

  }

  //--------------------------------------------------------------------
  // 4th step: assign fluid elements to bins
  //--------------------------------------------------------------------
  {
    for(std::map<int, std::set<int> >::const_iterator biniter=extendedfluidghosting.begin(); biniter!=extendedfluidghosting.end(); ++biniter)
    {
      DRT::MESHFREE::MeshfreeMultiBin* currbin = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(particledis_->gElement(biniter->first));
      for(std::set<int>::const_iterator fluideleiter=biniter->second.begin(); fluideleiter!=biniter->second.end(); ++fluideleiter)
      {
        int fluideleid = *fluideleiter;
        currbin->AddAssociatedFluidEle(fluideleid, fluiddis_->gElement(fluideleid));
//          cout << "in bin with id:" << currbin->Id() << " is fluid ele with id" << fluideleid << "with pointer" << fluiddis_->gElement(fluideleid) << endl;
      }
    }
  }

#ifdef DEBUG
  // check whether each particle has an underlying fluid element
  std::map<int,LINALG::Matrix<3,1> > currentpositions;
  for(int i=0; i<fluiddis_->NumMyColNodes(); i++)
  {
    DRT::Node* node = fluiddis_->lColNode(i);
    LINALG::Matrix<3,1> currpos;

    for (int a=0; a<3; a++)
    {
      currpos(a) = node->X()[a];
    }
    currentpositions.insert(std::pair<int,LINALG::Matrix<3,1> >(node->Id(),currpos));
  }
  // start loop over all particles
  for(int k=0; k<particledis_->NumMyColNodes(); k++)
  {
    DRT::Node* particle = particledis_->lColNode(k);
    const double* pos = particle->X();
    LINALG::Matrix<3,1> projpoint;
    for(int dim=0; dim<3; dim++)
      projpoint(dim) = pos[dim];
    bool foundele = false;
    for(int i=0; i<fluiddis_->NumMyColElements(); i++)
    {
      DRT::Element* fluidele = fluiddis_->lColElement(i);

      LINALG::Matrix<3,1> elecoord(true);
      const LINALG::SerialDenseMatrix xyze(GEO::getCurrentNodalPositions(Teuchos::rcp(fluidele,false), currentpositions));

      //get coordinates of the particle position in parameter space of the element
      foundele = GEO::currentToVolumeElementCoordinates(fluidele->Shape(), xyze, projpoint, elecoord);

      if(foundele == true)
        break;
    }
    if(foundele == false)
      dserror("particle (Id:%d) was found which does not have fluid support", particle->Id());
  }
#endif

  return;
}


/*----------------------------------------------------------------------*
| build connectivity from fluid elements to bins            ghamm 07/13 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::BuildElementToBinPointers()
{
  // first call base class to associate potential particle walls
  PARTICLE::Algorithm::BuildElementToBinPointers();

  // loop over column bins and fill fluid elements
  const int numcolbin = particledis_->NumMyColElements();
  for (int ibin=0; ibin<numcolbin; ++ibin)
  {
    DRT::Element* actele = particledis_->lColElement(ibin);
    DRT::MESHFREE::MeshfreeMultiBin* actbin = dynamic_cast<DRT::MESHFREE::MeshfreeMultiBin*>(actele);
    const int numfluidele = actbin->NumAssociatedFluidEle();
    const int* fluideleids = actbin->AssociatedFluidEleIds();
    std::vector<DRT::Element*> fluidelements(numfluidele);
    for(int iele=0; iele<numfluidele; iele++)
    {
      const int fluideleid = fluideleids[iele];
      fluidelements[iele] = fluiddis_->gElement(fluideleid);
    }
    actbin->BuildFluidElePointers(&fluidelements[0]);
  }

  return;
}


/*----------------------------------------------------------------------*
| single fields are tested                                  ghamm 11/12 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::TestResults(const Epetra_Comm& comm)
{
  DRT::Problem::Instance()->AddFieldTest(fluid_->CreateFieldTest());
  PARTICLE::Algorithm::TestResults(comm);
  return;
}


/*----------------------------------------------------------------------*
 | output particle time step                                ghamm 11/12  |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::Output()
{
  fluid_->Output();
  PARTICLE::Algorithm::Output();
  return;
}


/*----------------------------------------------------------------------*
 | get adjacent bins to corner, where ijk is in 1st octant ghamm 02/13  |
 *----------------------------------------------------------------------*/
std::vector<int> CAVITATION::Algorithm::AdjacentBinstoCorner(int* ijk)
{
  std::vector<int> adjbins;
  adjbins.reserve(8);

  // get all adjacent bins to the current corner, including the bin itself
  for(int i=-1;i<1;i++)
  {
    for(int j=-1;j<1;j++)
    {
      for(int k=-1;k<1;k++)
      {
        int ijk_neighbor[3] = {ijk[0]+i, ijk[1]+j, ijk[2]+k};

        int neighborgid = ConvertijkToGid(&ijk_neighbor[0]);
        if(neighborgid != -1)
        {
          adjbins.push_back(neighborgid);
        }
      } // end for int k
    } // end for int j
  } // end for int i

  return adjbins;
}


/*----------------------------------------------------------------------*
| setup of bubble inflow                                    ghamm 01/13 |
 *----------------------------------------------------------------------*/
void CAVITATION::Algorithm::BuildBubbleInflowCondition()
{
  // build inflow boundary condition
  std::vector<DRT::Condition*> conds;
  particledis_->GetCondition("ParticleInflow", conds);
  // unique bubbleinflow id over all inflow conditions
  int bubbleinflowid = 0;
  for (size_t i=0; i<conds.size(); i++)
  {
    if(i>0)
      dserror("only taken care of one particle inflow condition so far. "
          "Remedy: bubble_source_ needs to be a vector of the current layout");
    /*
     * inflow condition --> bubble sources
     *
     *  example: num_per_dir = {4, 5, 1}
     *
     *       <-> (dist_x = (vertex2_x-vertex1_x)/(num_per_dir_x-1))
     *
     *   x  x  x  x<-------- vertex2
     *
     *   x  x  x  x
     *
     *   x  x  x  x   ^
     *                | (dist_y = (vertex2_y-vertex1_y)/(num_per_dir_y-1) )
     *   x  x  x  x   ^
     *
     *   x  x  x  x
     *   ^
     *   |
     * vertex1
     *
     */

    // extract data from inflow condition
    const std::vector<double>* vertex1 = conds[i]->Get<std::vector<double> >("vertex1");
    const std::vector<double>* vertex2 = conds[i]->Get<std::vector<double> >("vertex2");
    const std::vector<int>* num_per_dir = conds[i]->Get<std::vector<int> >("num_per_dir");
    const std::vector<double>* inflow_vel = conds[i]->Get<std::vector<double> >("inflow_vel");
    int inflow_vel_curve = conds[i]->GetInt("inflow_vel_curve");
    double inflow_freq = conds[i]->GetDouble("inflow_freq");

    // make sure that a particle material is defined in the dat-file
    int id = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_particlemat);
    if (id==-1)
      dserror("Could not find particle material");

    const MAT::PAR::Parameter* mat = DRT::Problem::Instance()->Materials()->ParameterById(id);
    const MAT::PAR::ParticleMat* actmat = static_cast<const MAT::PAR::ParticleMat*>(mat);
    double initial_radius = actmat->initialradius_;

    double inflowtime = 1.0 / inflow_freq;
    if(std::abs(inflowtime/Dt() - (int)(inflowtime/Dt())) > EPS9)
      dserror("1/inflow_freq with inflow_freq = %f cannot be divided by fluid time step %f", inflowtime, Dt());
/* MUST BE ADDED WHEN PARTICLE CONTACT IS CONSIDERED
    double inflow_vel_mag = sqrt((*inflow_vel)[0]*(*inflow_vel)[0] + (*inflow_vel)[1]*(*inflow_vel)[1] + (*inflow_vel)[2]*(*inflow_vel)[2]);
    if(initial_radius/inflow_vel_mag > inflowtime)
      dserror("Overlap for inflowing bubbles expected: initial_radius/inflow_vel_mag = %f s > inflow_freq = %f s", initial_radius/inflow_vel_mag, inflowtime);
*/
    // loop over all bubble inflow positions and fill them into bin when they are on this proc;
    // up to here, only row bins are available
    std::vector<double> source_pos(3);
    for(int z=0; z<(*num_per_dir)[2]; ++z)
    {
      double dist_z = ((*vertex2)[2] - (*vertex1)[2]) / ((((*num_per_dir)[2]-1)!=0) ? ((*num_per_dir)[2]-1) : 1);
      source_pos[2] = (*vertex1)[2] + z * dist_z;
      for(int y=0; y<(*num_per_dir)[1]; ++y)
      {
        double dist_y = ((*vertex2)[1] - (*vertex1)[1]) / ((((*num_per_dir)[1]-1)!=0) ? ((*num_per_dir)[1]-1) : 1);
        source_pos[1] = (*vertex1)[1] + y * dist_y;
        for(int x=0; x<(*num_per_dir)[0]; ++x)
        {
          double dist_x = ((*vertex2)[0] - (*vertex1)[0]) / ((((*num_per_dir)[0]-1)!=0) ? ((*num_per_dir)[0]-1) : 1);
          source_pos[0] = (*vertex1)[0] + x * dist_x;
          // check whether this source position is on this proc
          int binId = ConvertPosToGid(source_pos);
          bool found = particledis_->HaveGlobalElement(binId);
          if(found == true)
          {
            Teuchos::RCP<BubbleSource> bubbleinflow = Teuchos::rcp(new BubbleSource(
                                                                          bubbleinflowid,
                                                                          source_pos,
                                                                          *inflow_vel,
                                                                          inflow_vel_curve,
                                                                          initial_radius,
                                                                          inflow_freq));
            bubble_source_[binId].push_back(bubbleinflow);
#ifdef DEBUG
            if(particledis_->gElement(binId)->Owner() != myrank_)
              dserror("Only row bins should show up here. Either add additional if-case or move ghosting to a later point in time.");
#endif
          }
          bubbleinflowid++;
        }
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | particle source                                         ghamm 02/13  |
 *----------------------------------------------------------------------*/
CAVITATION::BubbleSource::BubbleSource(
  int bubbleinflowid,
  std::vector<double> inflow_position,
  std::vector<double> inflow_vel,
  int inflow_vel_curve,
  double inflow_radius,
  double inflow_freq
  ) :
  inflowid_(bubbleinflowid),
  inflow_position_(inflow_position),
  inflow_vel_(inflow_vel),
  inflow_vel_curve_(inflow_vel_curve),
  inflow_radius_(inflow_radius),
  inflow_freq_(inflow_freq)
{
}
