/*----------------------------------------------------------------------*/
/*!
\file particleMeshFree_rendering.cpp

\brief rendering routine for smoothed particle hydrodynamics

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/
#include "particleMeshFree_rendering.H"
#include "particleMeshFree_rendering_resulttest.H"
#include "particleMeshFree_interaction.H"
#include "particleMeshFree_weightFunction.H"
#include "particle_algorithm.H"

#include "../drt_io/io.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../linalg/linalg_utils.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 | constructor                                             sfuchs 06/17 |
 *----------------------------------------------------------------------*/
PARTICLE::Rendering::Rendering(
    Teuchos::RCP<DRT::Discretization> discret,
    Teuchos::RCP<PARTICLE::Algorithm> particleAlgorithm,
    Teuchos::RCP<PARTICLE::ParticleMeshFreeInteractionHandler> interHandler,
    Teuchos::RCP<PARTICLE::WeightFunction_Base> weightFunctionHandler
    ) :
    discret_(discret),
    particle_algorithm_(particleAlgorithm),
    interHandler_(interHandler),
    weightFunctionHandler_(weightFunctionHandler),
    vel_(Teuchos::null),
    acc_(Teuchos::null),
    density_(Teuchos::null),
    specEnthalpy_(Teuchos::null),
    temperature_(Teuchos::null),
    pressure_(Teuchos::null)
{
  // row node map of rendering discretization
  Teuchos::RCP<Epetra_Map> noderowmap = Teuchos::rcp(new Epetra_Map(*discret_->NodeRowMap()));
  // fully overlapping node map
  Teuchos::RCP<Epetra_Map> rednodecolmap = LINALG::AllreduceEMap(*noderowmap);

  // row ele map of rendering discretization
  Teuchos::RCP<Epetra_Map> elerowmap = Teuchos::rcp(new Epetra_Map(*discret_->ElementRowMap()));
  // fully overlapping ele map
  Teuchos::RCP<Epetra_Map> redelecolmap = LINALG::AllreduceEMap(*elerowmap);

  // do the fully overlapping ghosting of the rendering discretization to have everything redundant
  discret_->ExportColumnNodes(*rednodecolmap);
  discret_->ExportColumnElements(*redelecolmap);

  // final fill complete
  discret_->FillComplete(true, false, false);

#ifdef DEBUG
  // some output to screen
  if (discret_->Comm().MyPID() == 0)
    std::cout << "after fully overlapping ghosting of the rendering discretization" << std::endl;
  DRT::UTILS::PrintParallelDistribution(*discret_);
#endif

  // set up connectivity bins to rendering nodes (overlapping)
  for (int lidNode = 0; lidNode < discret_->NumMyColNodes(); ++lidNode)
  {
    DRT::Node* currentNode = discret_->lColNode(lidNode);
    const int binID = particle_algorithm_->BinStrategy()->ConvertPosToGid(currentNode->X());
    binsToRenderingNodes_[binID].push_back(currentNode);
  }

  // write mesh file once
  Teuchos::RCP<IO::DiscretizationWriter> output = discret_->Writer();
  output->WriteMesh(particle_algorithm_->Step(), particle_algorithm_->Time());

  return;
} // PARTICLE::Rendering::Rendering

/*----------------------------------------------------------------------*
 | update the rendering vectors                            sfuchs 06/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::Rendering::UpdateRenderingVectors(
    Teuchos::RCP<DRT::Discretization> pDiscret,
    Teuchos::RCP<const Epetra_Vector> pDis,
    Teuchos::RCP<const Epetra_Vector> pVel,
    Teuchos::RCP<const Epetra_Vector> pAcc,
    Teuchos::RCP<const Epetra_Vector> pDensity,
    Teuchos::RCP<const Epetra_Vector> pSpecEnthalpy,
    Teuchos::RCP<const Epetra_Vector> pTemperature,
    Teuchos::RCP<const Epetra_Vector> pRadius,
    Teuchos::RCP<const Epetra_Vector> pPressure,
    Teuchos::RCP<const Epetra_Vector> pMass)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::Rendering::UpdateStateVectors");

  // re-create the dof-based rendering vectors
  vel_ = Teuchos::rcp(new Epetra_FEVector(*discret_->DofRowMap(), true));
  acc_ = Teuchos::rcp(new Epetra_FEVector(*discret_->DofRowMap(), true));

  // re-create the node-based rendering vectors
  density_ = Teuchos::rcp(new Epetra_FEVector(*discret_->NodeRowMap(), true));
  specEnthalpy_ = Teuchos::rcp(new Epetra_FEVector(*discret_->NodeRowMap(), true));
  temperature_ = Teuchos::rcp(new Epetra_FEVector(*discret_->NodeRowMap(), true));
  pressure_ = Teuchos::rcp(new Epetra_FEVector(*discret_->NodeRowMap(), true));

  // set that keeps track of the bins that have been already examined
  std::set<int> examinedbins;

  // loop over the particles (to check only the bins that own particles)
  for(int rowPar_i=0; rowPar_i<pDiscret->NodeRowMap()->NumMyElements(); ++rowPar_i)
  {
    // extract the particle
    DRT::Node *currparticle = pDiscret->lRowNode(rowPar_i);

    // find the bin it belongs to
    DRT::Element* currentBin = currparticle->Elements()[0];
    const int binId = currentBin->Id();

    // if a bin has already been examined --> continue with next particle
    if( examinedbins.find(binId) != examinedbins.end() )
      continue;
    //else: bin is examined for the first time --> new entry in examinedbins_
    examinedbins.insert(binId);

    // create and fill the list of neighboring rendering nodes
    std::list<DRT::Node*> nbrRdgNodes = GetNeighboringRenderingNodes(binId);

    // extract the pointer to the particles and loop over all particles in CurrentBin
    DRT::Node** currentBinParticles = currentBin->Nodes();
    for (int i=0; i<currentBin->NumNode(); ++i)
    {
      // determine the particle we are analizing
      DRT::Node* particle_i = currentBinParticles[i];
      // extract the gid
      const int gidNode_i = particle_i->Id();
      // determine the lid of the particle i (row style)
      const int lidRowNode_i = pDiscret->NodeRowMap()->LID(gidNode_i);

      // --- extract states of particle i ---

      const double density_i = (*pDensity)[lidRowNode_i];

      // skip non-interacting boundary particles
      if (density_i <= 0.0)
        continue;

      const double radius_i = (*pRadius)[lidRowNode_i];

      LINALG::Matrix<3,1> dis_i(true);
      LINALG::Matrix<3,1> vel_i(true);
      LINALG::Matrix<3,1> acc_i(true);
      for (int dim=0;dim<3;++dim)
      {
        dis_i(dim) = (*pDis)[3*lidRowNode_i + dim];
        vel_i(dim) = (*pVel)[3*lidRowNode_i + dim];
        acc_i(dim) = (*pAcc)[3*lidRowNode_i + dim];
      }

      const double mass_i = (*pMass)[lidRowNode_i];
      const double massOverDensity_i = mass_i / density_i;

      const double specEnthalpy_i = (*pSpecEnthalpy)[lidRowNode_i];
      const double temperature_i = (*pTemperature)[lidRowNode_i];
      const double pressure_i = (*pPressure)[lidRowNode_i];

      // loop over all neighboring rendering nodes
      std::list<DRT::Node*>::const_iterator nbrRdgNode_k;
      for (nbrRdgNode_k = nbrRdgNodes.begin(); nbrRdgNode_k != nbrRdgNodes.end(); ++nbrRdgNode_k)
      {
        LINALG::Matrix<3,1> disNbrRdgNode_k(true);
        for (int dim=0; dim<3; ++dim)
        {
          disNbrRdgNode_k(dim) = ((*nbrRdgNode_k)->X())[dim];
        }

        //In case particle i and rendering node k are close to two opposite periodic boundaries, the position of rendering node k is shifted correspondingly.
        //It is sufficient to only shift disNbrRdgNode_k since the terms evaluated in the following do only depend on the relative distance between these two
        //positions but not on the absolute positions.
        interHandler_->ShiftPeriodicBoundaryPair(dis_i, disNbrRdgNode_k, radius_i, 0.0);

        LINALG::Matrix<3,1> rRel(true);
        for (int dim=0; dim<3; ++dim)
        {
          rRel(dim) = dis_i(dim) - disNbrRdgNode_k(dim);
        }
        const double rRelNorm2 = rRel.Norm2();

        // skip in case particles are too apart
        if (rRelNorm2 > radius_i)
          continue;

        // rendering currently via Monaghan 2005, equation (2.8)
        // using the same kernel as in the simulation

        // determine weight function
        const double weight_ik = weightFunctionHandler_->W(rRelNorm2, radius_i);
        const double massOverDensityWeight_ik = massOverDensity_i*weight_ik;

        // render states of particle i to rendering node k
        LINALG::Matrix<3,1> vel_ik(vel_i);
        LINALG::Matrix<3,1> acc_ik(acc_i);
        vel_ik.Scale(massOverDensityWeight_ik);
        acc_ik.Scale(massOverDensityWeight_ik);

        double density_ik = density_i * massOverDensityWeight_ik; // = mass_i * weight_ik;
        double specEnthalpy_ik = specEnthalpy_i * massOverDensityWeight_ik;
        double temperature_ik = temperature_i * massOverDensityWeight_ik;
        double pressure_ik = pressure_i * massOverDensityWeight_ik;

        // extract the gid
        const int gidRdgNode_k = (*nbrRdgNode_k)->Id();

        // extract global dof ids and fill into lmRdgNode_k
        std::vector<int> lmRdgNode_k;
        lmRdgNode_k.reserve(3);
        discret_->Dof(*nbrRdgNode_k, lmRdgNode_k);

        // write the rendering vectors
        // use of Epetra_FEVectors to take care of ghosted rendering nodes
        vel_->SumIntoGlobalValues(3, &lmRdgNode_k[0], &vel_ik(0));
        acc_->SumIntoGlobalValues(3, &lmRdgNode_k[0], &acc_ik(0));

        density_->SumIntoGlobalValues(1, &gidRdgNode_k, &density_ik);
        specEnthalpy_->SumIntoGlobalValues(1, &gidRdgNode_k, &specEnthalpy_ik);
        temperature_->SumIntoGlobalValues(1, &gidRdgNode_k, &temperature_ik);
        pressure_->SumIntoGlobalValues(1, &gidRdgNode_k, &pressure_ik);
      }
    }
  }

  // assemble rendering vectors
  int err = 0;
  err += vel_->GlobalAssemble(Add, false);
  err += acc_->GlobalAssemble(Add, true);
  err += density_->GlobalAssemble(Add, true);
  err += specEnthalpy_->GlobalAssemble(Add, true);
  err += temperature_->GlobalAssemble(Add, true);
  err += pressure_->GlobalAssemble(Add, true);
  if (err!=0)
  {
    dserror("global assemble of rendering vectors failed!");
  }

  return;
} // PARTICLE::Rendering::UpdateStateVectors

/*----------------------------------------------------------------------*
 | get neighboring rendering nodes                        sfuchs 06/17 |
 *----------------------------------------------------------------------*/
std::list<DRT::Node*> PARTICLE::Rendering::GetNeighboringRenderingNodes(const int binId)
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::Rendering::GetNeighboringRenderingNodes");

  // create an empty list
  std::list<DRT::Node*> neighboringRenderingNodes;

  // find the neighboring bins and fill binIds
  std::vector<int> binIds;
  binIds.reserve(27);
  particle_algorithm_->BinStrategy()->GetNeighborAndOwnBinIds(binId,binIds);

  // concatenate the various rendering node lists
  for (std::vector<int>::const_iterator ii = binIds.begin(); ii != binIds.end(); ++ii)
  {
    neighboringRenderingNodes.insert(neighboringRenderingNodes.end(), binsToRenderingNodes_[*ii].begin(), binsToRenderingNodes_[*ii].end());
  }

  return neighboringRenderingNodes;
} // PARTICLE::Rendering::GetNeighboringRenderingNodes

/*----------------------------------------------------------------------*
 | rendering output state                                  sfuchs 06/17 |
 *----------------------------------------------------------------------*/
void PARTICLE::Rendering::OutputState()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLE::Rendering::OutputState");

  Teuchos::RCP<IO::DiscretizationWriter> output = discret_->Writer();

  // extract step and time
  const int step = particle_algorithm_->Step();
  const double time = particle_algorithm_->Time();

  output->NewStep(step, time);

  // write dof-based vectors
  output->WriteVector("velocity", vel_, IO::dofvector);
  output->WriteVector("acceleration", acc_, IO::dofvector);

  // write node-based vectors
  output->WriteVector("density", density_, IO::nodevector);
  output->WriteVector("specEnthalpy", specEnthalpy_, IO::nodevector);
  output->WriteVector("temperature", temperature_, IO::nodevector);
  output->WriteVector("pressure", pressure_, IO::nodevector);

  return;
} // PARTICLE::Rendering::OutputState

/*----------------------------------------------------------------------*
 | create result test for rendering                        sfuchs 06/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> PARTICLE::Rendering::CreateFieldTest()
{
  return Teuchos::rcp(new ParticleMeshfreeRenderingResultTest(*this));
} // PARTICLE::Rendering::CreateFieldTest
