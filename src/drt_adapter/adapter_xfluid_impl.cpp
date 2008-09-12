/*----------------------------------------------------------------------*/
/*!
\file adapter_xfluid_impl.cpp

\brief Fluid field adapter

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "adapter_xfluid_impl.H"

#include "../drt_fluid/xfluidresulttest.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/linalg_blocksparsematrix.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>

extern struct _FILES  allfiles;
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::XFluidImpl::XFluidImpl(
        Teuchos::RCP<DRT::Discretization> dis,
        const Teuchos::RCP<DRT::Discretization> soliddis,
        Teuchos::RCP<LINALG::Solver> solver,
        Teuchos::RCP<ParameterList> params,
        Teuchos::RCP<IO::DiscretizationWriter> output,
        bool isale)
  : fluid_(dis, *solver, *params, *output, isale),
    dis_(dis),
    solver_(solver),
    params_(params),
    output_(output)
{
  vector<string> conditions_to_copy;
  conditions_to_copy.push_back("FSICoupling");
  conditions_to_copy.push_back("XFEMCoupling");
  boundarydis_ = DRT::UTILS::CreateDiscretizationFromCondition(soliddis, "FSICoupling", "Boundary", "BELE3", conditions_to_copy);
  dsassert(boundarydis_->NumGlobalNodes() > 0, "empty discretization detected. FSICoupling condition applied?");
  
//  boundaryoutput_ = rcp(new IO::DiscretizationWriter(boundarydis_));
//  boundaryoutput_->WriteMesh(0,0.0);

  // create node and element distribution with elements and nodes ghosted on all processors
  const Epetra_Map noderowmap = *boundarydis_->NodeRowMap();
  std::cout << "noderowmap->UniqueGIDs(): " << noderowmap.UniqueGIDs() << endl;
  //std::cout << noderowmap << endl;

  Teuchos::RCP<Epetra_Map> newnodecolmap = LINALG::AllreduceEMap(noderowmap);
  //std::cout << *newnodecolmap << endl;

  DRT::UTILS::RedistributeWithNewNodalDistribution(*boundarydis_, noderowmap, *newnodecolmap);

  DRT::UTILS::SetupNDimExtractor(*boundarydis_,"FSICoupling",interface_);
  DRT::UTILS::SetupNDimExtractor(*boundarydis_,"FREESURFCoupling",freesurface_);

  // create interface DOF vectors using the solid parallel distribution
  const Epetra_Map* fluidsurface_dofrowmap = boundarydis_->DofRowMap();
  idispnp_    = LINALG::CreateVector(*fluidsurface_dofrowmap,true);
  ivelnp_     = LINALG::CreateVector(*fluidsurface_dofrowmap,true);
  itrueresnp_ = LINALG::CreateVector(*fluidsurface_dofrowmap,true);

  idispn_   = LINALG::CreateVector(*fluidsurface_dofrowmap,true);
  iveln_    = LINALG::CreateVector(*fluidsurface_dofrowmap,true);
  iaccn_    = LINALG::CreateVector(*fluidsurface_dofrowmap,true);

  fluid_.SetFreeSurface(&freesurface_);
  std::cout << "XFluidImpl constructor done" << endl;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluidImpl::InitialGuess()
{
  dserror("not implemented");
  return fluid_.InitialGuess();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluidImpl::RHS()
{
  dserror("not implemented");
  return fluid_.Residual();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluidImpl::Velnp()
{
  dserror("not implemented");
  return fluid_.Velnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluidImpl::Veln()
{
  dserror("not implemented");
  return fluid_.Veln();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluidImpl::Dispnp()
{
  dserror("not implemented");
//  return fluid_.Dispnp();
  return null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::XFluidImpl::DofRowMap()
{
  dserror("not implemented");
  const Epetra_Map* dofrowmap = dis_->DofRowMap();
  return Teuchos::rcp(dofrowmap, false);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::XFluidImpl::SystemMatrix()
{
  dserror("not implemented");
  // if anything (e.g. monolithic FSI) we give fluid coupling and interface DOF combined back
  return fluid_.SystemMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase> ADAPTER::XFluidImpl::BlockSystemMatrix()
{
  dserror("no block matrix here");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase> ADAPTER::XFluidImpl::MeshMoveMatrix()
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::XFluidImpl::Discretization()
{
  return boundarydis_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Teuchos::RCP<Epetra_Vector> ADAPTER::XFluidImpl::StructCondRHS() const
// {
//   return interface_.ExtractCondVector(Velnp());
// }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::PrepareTimeStep()
{
  fluid_.PrepareTimeStep();

  // we add the whole fluid mesh displacement later on?
  //fluid_.Dispnp()->PutScalar(0.);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::Evaluate(Teuchos::RCP<const Epetra_Vector> vel)
{
  if (vel!=Teuchos::null)
  {
    fluid_.Evaluate(vel);
  }
  else
  {
    fluid_.Evaluate(Teuchos::null);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::Update()
{
  fluid_.TimeUpdate();

  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
  const double dt = fsidyn.get<double>("TIMESTEP");

  // compute acceleration at timestep n+1
  Teuchos::RCP<Epetra_Vector> iaccnp = rcp(new Epetra_Vector(iaccn_->Map()));
  Teuchos::RCP<Epetra_Vector> ivelnp = rcp(new Epetra_Vector(iveln_->Map()));
  const double theta = 1.0;
  iaccnp->Update(-(1.0-theta)/(theta),*iaccn_,0.0);
  iaccnp->Update(1.0/(theta*dt),*ivelnp_,-1.0/(theta*dt),*iveln_,1.0);


//  double beta;
//  double gamma;
//  if (Teuchos::getIntegralValue<int>(fsidyn,"SECONDORDER") == 1 and not Step()==1)
//  {
//    gamma = 0.5;
//    beta = gamma/2.0;
//  }
//  else
//  {
//    cout << "first step" << endl;
//    //dserror("nope");
//    gamma = 1.0;
//    beta = gamma/2.0;
//  }
//  
//  // compute acceleration at timestep n+1
//  iaccnp->Update(-(1.0-(2.0*beta))/(2.0*beta),*iaccn_,0.0);
//  iaccnp->Update(-1.0/(beta*dt),*iveln_,1.0);
//  iaccnp->Update(1.0/(beta*dt*dt),*idispnp_,-1.0/(beta*dt*dt),*idispn_,1.0);
//  
//  // compute velocity at timestep n+1
//  ivelnp->Update(1.0,*iveln_,0.0);
//  ivelnp->Update(gamma*dt,*iaccnp,(1-gamma)*dt,*iaccn_,1.0);
  
  // update acceleration n
  iaccn_->Update(1.0,*iaccnp,0.0);

  // update velocity n
  iveln_->Update(1.0,*ivelnp_,0.0);
  
  // update displacement n
  idispn_->Update(1.0,*idispnp_,0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::Output()
{
  fluid_.Output();

//  boundaryoutput_->NewStep(Step(),Time());
//  boundaryoutput_->WriteVector("interface displacement", idisp_);
//  boundaryoutput_->WriteVector("interface velocity", ivel_);
//  boundaryoutput_->WriteVector("interface velocity (n)", iveln_);
//  boundaryoutput_->WriteVector("interface velocity (n-1)", ivelnm_);
//  boundaryoutput_->WriteVector("interface acceleration (n)", iaccn_);
//  boundaryoutput_->WriteVector("interface force", itrueres_);

  // create interface DOF vectors using the fluid parallel distribution
  const Epetra_Map* fluidsurface_dofcolmap = boundarydis_->DofColMap();
  Teuchos::RCP<Epetra_Vector> ivelnpcol   = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  Teuchos::RCP<Epetra_Vector> ivelncol    = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  Teuchos::RCP<Epetra_Vector> iaccncol    = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  Teuchos::RCP<Epetra_Vector> idispnpcol  = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  Teuchos::RCP<Epetra_Vector> itruerescol = LINALG::CreateVector(*fluidsurface_dofcolmap,true);

  // map to fluid parallel distribution
  LINALG::Export(*idispnp_ ,*idispnpcol);
  LINALG::Export(*ivelnp_  ,*ivelnpcol);
  LINALG::Export(*iveln_   ,*ivelncol);
  LINALG::Export(*iaccn_   ,*iaccncol);

  LINALG::Export(*itrueresnp_,*itruerescol);

  PrintInterfaceVectorField(idispnpcol, itruerescol, "_solution_iforce_", "interface force");
  PrintInterfaceVectorField(idispnpcol, ivelnpcol, "_solution_ivel_"  , "interface velocity n+1");
  PrintInterfaceVectorField(idispnpcol, ivelncol , "_solution_iveln_" , "interface velocity n");
  PrintInterfaceVectorField(idispnpcol, iaccncol , "_solution_iaccn_" , "interface acceleration n");

}

void ADAPTER::XFluidImpl::PrintInterfaceVectorField(
    const Teuchos::RCP<Epetra_Vector>   displacementfield,
    const Teuchos::RCP<Epetra_Vector>   vectorfield,
    const std::string filestr,
    const std::string name_in_gmsh
    )
{
  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = (bool)getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT");
  if (gmshdebugout)
  {
    std::stringstream filename;
    std::stringstream filenamedel;
    filename << allfiles.outputfile_kenner << filestr << std::setw(5) << setfill('0') << Step() << ".pos";
    filenamedel << allfiles.outputfile_kenner << filestr << std::setw(5) << setfill('0') << Step()-5 << ".pos";
    std::remove(filenamedel.str().c_str());
    std::cout << "writing " << left << std::setw(50) <<filename.str()<<"...";
    std::ofstream f_system(filename.str().c_str());
    
    {
      stringstream gmshfilecontent;
      gmshfilecontent << "View \" " << name_in_gmsh << " \" {" << endl;
      for (int i=0; i<boundarydis_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = boundarydis_->lColElement(i);
        //      cout << *actele << endl;
        vector<int> lm;
        vector<int> lmowner;
        actele->LocationVector(*boundarydis_, lm, lmowner);
        
        // extract local values from the global vector
        vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*vectorfield, myvelnp, lm);
        
        vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*displacementfield, mydisp, lm);
        
        const int nsd = 3;
        const int numnode = actele->NumNode();
        BlitzMat elementvalues(nsd,numnode);
        BlitzMat elementpositions(nsd,numnode);
        int counter = 0;
        for (int iparam=0; iparam<numnode; ++iparam)
        {
          const DRT::Node* node = actele->Nodes()[iparam];
          //        cout << *node << endl;
          const double* pos = node->X();
          for (int isd = 0; isd < nsd; ++isd)
          {
            elementvalues(isd,iparam) = myvelnp[counter];
            elementpositions(isd,iparam) = pos[isd] + mydisp[counter];
            counter++;
          }
        }
        //      cout << elementpositions << endl;
        //      exit(1);
        
        gmshfilecontent << IO::GMSH::cellWithVectorFieldToString(
            actele->Shape(), elementvalues, elementpositions) << endl;
      }
      gmshfilecontent << "};" << endl;
      f_system << gmshfilecontent.str();
    }
    f_system.close();
    std::cout << " done" << endl;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::NonlinearSolve()
{

  //cout << "XFluidImpl::NonlinearSolve()" << endl;

  // create interface DOF vectors using the fluid parallel distribution
  const Epetra_Map* fluidsurface_dofcolmap = boundarydis_->DofColMap();
  Teuchos::RCP<Epetra_Vector> idispcolnp = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  Teuchos::RCP<Epetra_Vector> idispcoln  = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  Teuchos::RCP<Epetra_Vector> ivelcolnp  = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  Teuchos::RCP<Epetra_Vector> ivelcoln   = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  Teuchos::RCP<Epetra_Vector> iacccoln   = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
//  Teuchos::RCP<Epetra_Vector> itruerescol = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
//
//  Teuchos::RCP<Epetra_Vector> ivelncol = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
//  Teuchos::RCP<Epetra_Vector> iaccncol = LINALG::CreateVector(*fluidsurface_dofcolmap,true);

  // map to fluid parallel distribution
  LINALG::Export(*idispnp_,*idispcolnp);
  LINALG::Export(*idispn_ ,*idispcoln);
  LINALG::Export(*ivelnp_ ,*ivelcolnp);
  LINALG::Export(*iveln_  ,*ivelcoln);
  LINALG::Export(*iaccn_  ,*iacccoln);


//  cout << "SetState" << endl;
  boundarydis_->SetState("idispcolnp",idispcolnp);
  boundarydis_->SetState("idispcoln" ,idispcoln);
  
  boundarydis_->SetState("ivelcolnp" ,ivelcolnp);
  boundarydis_->SetState("ivelcoln"  ,ivelcoln);
  boundarydis_->SetState("iacccoln"  ,iacccoln);
  
  fluid_.NonlinearSolve(boundarydis_);
  
  Teuchos::RCP<const Epetra_Vector> itruerescol = boundarydis_->GetState("iforcenp");
  
  boundarydis_->ClearState();
  

  // map back to solid parallel distribution
  LINALG::Export(*itruerescol,*itrueresnp_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::XFluidImpl::InnerVelocityRowMap()
{
  // build inner velocity map
  // dofs at the interface are excluded
  // we use only velocity dofs and only those without Dirichlet constraint

  Teuchos::RCP<const Epetra_Map> velmap = fluid_.VelocityRowMap(); //???
  Teuchos::RCP<Epetra_Vector> dirichtoggle = fluid_.Dirichlet();   //???
  Teuchos::RCP<const Epetra_Map> fullmap = DofRowMap();            //???

  const int numvelids = velmap->NumMyElements();
  std::vector<int> velids;
  velids.reserve(numvelids);
  for (int i=0; i<numvelids; ++i)
  {
    int gid = velmap->GID(i);
    // NOTE: in xfem, there are no interface dofs in the fluid field
    if ((*dirichtoggle)[fullmap->LID(gid)]==0.)
    {
      velids.push_back(gid);
    }
  }

  innervelmap_ = Teuchos::rcp(new Epetra_Map(-1,velids.size(), &velids[0], 0, velmap->Comm()));

  return innervelmap_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::XFluidImpl::VelocityRowMap()
{
  return fluid_.VelocityRowMap();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::XFluidImpl::PressureRowMap()
{
  return fluid_.PressureRowMap();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::SetMeshMap(Teuchos::RCP<const Epetra_Map> mm)
{
  dserror("makes no sense here");
  //meshmap_.Setup(*dis_->DofRowMap(),mm,LINALG::SplitMap(*dis_->DofRowMap(),*mm));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::XFluidImpl::ResidualScaling() const
{
  return fluid_.ResidualScaling();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::XFluidImpl::TimeScaling() const
{
  return 1./fluid_.Dt();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::ReadRestart(int step)
{
  fluid_.ReadRestart(step);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::XFluidImpl::Time()
{
  return fluid_.Time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int ADAPTER::XFluidImpl::Step()
{
  return fluid_.Step();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::LiftDrag()
{
  fluid_.LiftDrag();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::XFluidImpl::ExtractInterfaceForces()
{
  return interface_.ExtractCondVector(itrueresnp_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::XFluidImpl::ExtractInterfaceForcesRobin()
{
  dserror("no Robin coupling here");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::XFluidImpl::ExtractInterfaceFluidVelocity()
{
  dserror("no Robin coupling here");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::XFluidImpl::ExtractInterfaceVeln()
{
  return interface_.ExtractCondVector(iveln_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
{
  interface_.InsertCondVector(ivel,ivelnp_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::ApplyInterfaceRobinValue(Teuchos::RCP<Epetra_Vector> ivel, Teuchos::RCP<Epetra_Vector> iforce)
{
  dserror("robin unimplemented");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::ApplyMeshDisplacement(Teuchos::RCP<Epetra_Vector> idisp)
{
  interface_.InsertCondVector(idisp,idispnp_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::ApplyMeshVelocity(Teuchos::RCP<Epetra_Vector> gridvel)
{
  dserror("makes no sense here!");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx)
{
  dserror("not implemented!");
  // get interface velocity at t(n)
  const Teuchos::RCP<Epetra_Vector> veln = Interface().ExtractCondVector(Veln());

  // We convert Delta d(n+1,i+1) to Delta u(n+1,i+1) here.
  //
  // Delta d(n+1,i+1) = ( Delta u(n+1,i+1) + u(n) ) * dt
  //
  fcx->Update(-1.,*veln,TimeScaling());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::VelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx)
{
  dserror("not implemented!");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int ADAPTER::XFluidImpl::Itemax() const
{
  return fluid_.Itemax();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::SetItemax(int itemax)
{
  fluid_.SetItemax(itemax);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::XFluidImpl::IntegrateInterfaceShape()
{
  dserror("not implemented!");
  return interface_.ExtractCondVector(fluid_.IntegrateInterfaceShape("FSICoupling"));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::UseBlockMatrix(const LINALG::MultiMapExtractor& domainmaps,
                                         const LINALG::MultiMapExtractor& rangemaps,
                                         bool splitmatrix)
{
  dserror("no, probably not");
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::XFluidImpl::RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel)
{
  dserror("not implemented!");
  const Epetra_Map* dofrowmap = Discretization()->DofRowMap();
  Teuchos::RCP<Epetra_Vector> relax = LINALG::CreateVector(*dofrowmap,true);
  interface_.InsertCondVector(ivel,relax);
  //fluid_.LinearRelaxationSolve(relax);
  return ExtractInterfaceForces();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::XFluidImpl::CreateFieldTest()
{
  return Teuchos::rcp(new FLD::XFluidResultTest(fluid_));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluidImpl::ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres)
{
  dserror("not implemented!");
  return (fluid_.VelPresSplitter()).ExtractOtherVector(velpres);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::SetInitialFlowField(int whichinitialfield,int startfuncno)
{
  dserror("not implemented!");
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::SetTimeLomaFields(RCP<const Epetra_Vector> noddensn,RCP<const Epetra_Vector> noddensnm)
{
   dserror("not implemented!");
   return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::SetIterLomaFields(RCP<const Epetra_Vector> noddensnp)
{
   dserror("not implemented!");
   return;
}

#endif  // #ifdef CCADISCRET
