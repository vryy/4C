/*----------------------------------------------------------------------*/
/*!
\file particleMeshFree_rendering.cpp

\brief Handler of the rendering-related stuff

\level 1

\maintainer Cattabiani Alessandro
*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
/* headers */
#include "../drt_lib/drt_discret.H"
#include "../drt_io/io.H"

#include "particleMeshFree_weightFunction.H"
#include "particleMeshFree_rendering.H"
#include "particle_algorithm.H"

/*----------------------------------------------------------------------*
 | constructor                                             katta 10/16  |
 *----------------------------------------------------------------------*/
PARTICLE::Rendering::Rendering(Teuchos::RCP<PARTICLE::Algorithm> particleAlgorithm) :
  discret_(DRT::Problem::Instance()->GetDis("rendering")),
  particle_algorithm_(particleAlgorithm),
  weightFunctionHandler_(Teuchos::null),
  trg_writeMesh_(true)
{
  // --- finalize the discret set up --- //

  discret_->FillComplete(false,false,false);

  // --- set up connectivity --- //

  // create and fill bins2renderingNodes
  for (int lidNode = 0; lidNode < discret_->NumMyRowNodes(); ++lidNode)
  {
    DRT::Node* currentNode = discret_->lRowNode(lidNode);
    const int binID = particle_algorithm_->BinStrategy()->ConvertPosToGid(currentNode->X());
    bins2renderingNodes_[binID].push_back(currentNode);
  }

  // --- set up state vectors --- //

  // create the dof-based vectors
  vel_ = Teuchos::rcp(new Epetra_FEVector(*discret_->DofRowMap(), true));
  acc_ = Teuchos::rcp(new Epetra_FEVector(*discret_->DofRowMap(), true));

  // create the node-based vectors
  density_ = Teuchos::rcp(new Epetra_FEVector(*discret_->NodeRowMap(), true));
  specEnthalpy_ = Teuchos::rcp(new Epetra_FEVector(*discret_->NodeRowMap(), true));
  temperature_ = Teuchos::rcp(new Epetra_FEVector(*discret_->NodeRowMap(), true));
  pressure_ = Teuchos::rcp(new Epetra_FEVector(*discret_->NodeRowMap(), true));

  // set the correct WeightFunction

  switch (DRT::INPUT::IntegralValue<INPAR::PARTICLE::WeightFunction>(DRT::Problem::Instance()->ParticleParams(),"WEIGHT_FUNCTION"))
  {
  case INPAR::PARTICLE::CubicBspline :
  {
    weightFunctionHandler_ = Teuchos::rcp(new PARTICLE::WeightFunction_CubicBspline);
    break;
  }
  case INPAR::PARTICLE::SqrtHyperbola :
  {
    weightFunctionHandler_ = Teuchos::rcp(new PARTICLE::WeightFunction_SqrtHyperbola);
    break;
  }
  case INPAR::PARTICLE::HyperbolaNoRsz :
  {
    weightFunctionHandler_ = Teuchos::rcp(new PARTICLE::WeightFunction_HyperbolaNoRsz);
    break;
  }
  }

}

/*----------------------------------------------------------------------*
 | update state vectors                                    katta 10/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::Rendering::UpdateStateVectors(
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
  // set to 0 the renderingStateVectors
  vel_->PutScalar(0.0);
  acc_->PutScalar(0.0);
  density_->PutScalar(0.0);
  specEnthalpy_->PutScalar(0.0);
  temperature_->PutScalar(0.0);
  pressure_->PutScalar(0.0);

  // set that keeps track of the bins that have been already examined
  std::set<int> examinedbins;



  // loop over the particles (to check only the bins that own particles)
  for(int rowPar_i=0; rowPar_i<pDiscret->NodeRowMap()->NumMyElements(); ++rowPar_i)
  {
    // extract the particle
    DRT::Node *currparticle = pDiscret->lRowNode(rowPar_i);

    //find the bin it belongs to
    DRT::Element* currentBin = currparticle->Elements()[0];

    const int binId = currentBin->Id();
    // if a bin has already been examined --> continue with next particle
    if( examinedbins.find(binId) != examinedbins.end() )
      continue;
    //else: bin is examined for the first time --> new entry in examinedbins_
    examinedbins.insert(binId);

    // create and fill the list of neighbouring nodes
    std::list<DRT::Node*> neighboringRenderingNodes = GetNeighbouringRenderingNodes(binId);

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

      // extract particle state vectors
        // useful state vectors that does not appear in the rendering output
        const double massOverDensity_i =  (*pMass)[lidRowNode_i] / (*pDensity)[lidRowNode_i];
        const double radius_i = (*pRadius)[lidRowNode_i];
        //rendering output
        LINALG::Matrix<3,1> vel, acc;
        for (int dim=0;dim<3;++dim)
        {
          vel(dim) = massOverDensity_i * (*(*vel_)(0))[3*lidRowNode_i + dim];
          acc(dim) = massOverDensity_i * (*(*acc_)(0))[3*lidRowNode_i + dim];
        }
        const double density_i = (*pMass)[lidRowNode_i]; // multiplying massOverDensity * density is the same as considering just the mass
        const double specEnthalpy_i = massOverDensity_i * (*pSpecEnthalpy)[lidRowNode_i];
        const double temperature_i = massOverDensity_i * (*pTemperature)[lidRowNode_i];
        const double pressure_i = massOverDensity_i * (*pPressure)[lidRowNode_i];

      for (std::list<DRT::Node*>::const_iterator neighboringRenderingNode_i  = neighboringRenderingNodes.begin();
                                                 neighboringRenderingNode_i != neighboringRenderingNodes.end();
                                                 ++neighboringRenderingNode_i)
      {
        const double* disNeighboringRenderingNode_i = (*neighboringRenderingNode_i)->X();

        LINALG::Matrix<3,1> rRel;
        for (int dim=0; dim<3; ++dim)
        {
          rRel(dim) = (*pDis)[3*lidRowNode_i + dim] - disNeighboringRenderingNode_i[dim];
        }
        const double rRelNorm2 = rRel.Norm2();

        // skip in case particles are too apart
        if (rRelNorm2< radius_i)
        {
          // give me the proper weight
          const double weight = weightFunctionHandler_->Weight(rRelNorm2, radius_i);

          // specify the variables with the particular weight
          LINALG::Matrix<3,1> velWeight_i(vel), accWeight_i(acc);
          velWeight_i.Scale(weight);
          accWeight_i.Scale(weight);
          const double densityWeight_i = density_i * weight; // multiplying massOverDensity * density is the same as considering just the mass
          const double specEnthalpyWeight_i = specEnthalpy_i * weight;
          const double temperatureWeight_i = temperature_i * weight;
          const double pressureWeight_i = pressure_i * weight;

          // determine the renderingNode we are analizing
          DRT::Node* renderingNode_i = *neighboringRenderingNode_i;
          // extract the gid
          const int gidRenderingNode_i = renderingNode_i->Id();

          std::vector<int> lmRenderingNode_i;
          lmRenderingNode_i.reserve(3);
          // extract global dof ids and fill into lmRenderingNode_i
          discret_->Dof(renderingNode_i, lmRenderingNode_i);

          // write the renderingVectors (even if the particles are ghosted, we use FE vectors to take care of that)
          vel_->SumIntoGlobalValues(3,&lmRenderingNode_i[0],&velWeight_i(0));
          acc_->SumIntoGlobalValues(3,&lmRenderingNode_i[0],&accWeight_i(0));

          density_->SumIntoGlobalValues(1,&gidRenderingNode_i,&densityWeight_i);
          specEnthalpy_->SumIntoGlobalValues(1,&gidRenderingNode_i,&specEnthalpyWeight_i);
          temperature_->SumIntoGlobalValues(1,&gidRenderingNode_i,&temperatureWeight_i);
          pressure_->SumIntoGlobalValues(1,&gidRenderingNode_i,&pressureWeight_i);
        }
      }

    }

  }

  // assemble everything
  int err = 0;
  err += vel_->GlobalAssemble(Add, false);
  err += acc_->GlobalAssemble(Add, true);
  err += density_->GlobalAssemble(Add, true);
  err += specEnthalpy_->GlobalAssemble(Add, true);
  err += temperature_->GlobalAssemble(Add, true);
  err += pressure_->GlobalAssemble(Add, true);
  if (err!=0)
  {
    dserror("global assembling in the output failed");
  }

}


/*----------------------------------------------------------------------*
 | get neighbouring rendering nodes                        katta 10/16  |
 *----------------------------------------------------------------------*/
std::list<DRT::Node*> PARTICLE::Rendering::GetNeighbouringRenderingNodes(const int binId)
{
  // create an empty list
  std::list<DRT::Node*> neighboringRenderingParticles;

  // find the neighbouring bins and fill binIds
  std::vector<int> binIds;
  binIds.reserve(27);
  particle_algorithm_->BinStrategy()->GetNeighborAndOwnBinIds(binId,binIds);

  // concatenate the various rendering node lists
  for (std::vector<int>::const_iterator ii = binIds.begin(); ii != binIds.end(); ++ii)
  {
    neighboringRenderingParticles.insert(neighboringRenderingParticles.end(), bins2renderingNodes_[*ii].begin(), bins2renderingNodes_[*ii].end());
  }

  return neighboringRenderingParticles;
}


/*----------------------------------------------------------------------*
 | Rendering output state                                  katta 10/16  |
 *----------------------------------------------------------------------*/
void PARTICLE::Rendering::OutputState()
{
  Teuchos::RCP<IO::DiscretizationWriter> output = discret_->Writer();

  // extract step and time
  const int step = particle_algorithm_->Step();
  const double time = particle_algorithm_->Time();

  if (trg_writeMesh_)
  {
    output->WriteMesh(step, time);
  }

  output->NewStep(step, time);

  // write dof-based vectors
  output->WriteVector("velocity", vel_);
  output->WriteVector("acceleration", acc_);

  // write node-based vectors
  output->WriteVector("density", density_, output->nodevector);
  output->WriteVector("specEnthalpy", specEnthalpy_, output->nodevector);
  output->WriteVector("temperature", temperature_, output->nodevector);
  output->WriteVector("pressure", pressure_, output->nodevector);

  trg_writeMesh_ = false;
}

