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

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Epetra_Export.h>

#include "../drt_fluid/xfluidresulttest.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/linalg_blocksparsematrix.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_control.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_geometry/intersection_service.H"
#include "../drt_geometry/element_volume.H"
#include "../drt_fem_general/debug_nan.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::XFluidImpl::XFluidImpl(
        Teuchos::RCP<DRT::Discretization> dis,
        const Teuchos::RCP<DRT::Discretization> soliddis,
        Teuchos::RCP<ParameterList> params)
  : fluid_(dis, *params),
    params_(params)
{

  //if (!soliddis->Filled()) cout << " solid dis is not filled"<< endl;
  // test intersection without XFEM
  //Teuchos::RCP< UTILS::PotentialManager > pot_man =rcp(new UTILS::PotentialManager(dis, soliddis, *soliddis));
  // test intersection without XFEM

  vector<string> conditions_to_copy;
  conditions_to_copy.push_back("FSICoupling");
  conditions_to_copy.push_back("XFEMCoupling");
  boundarydis_ = DRT::UTILS::CreateDiscretizationFromCondition(soliddis, "FSICoupling", "boundary", "BELE3", conditions_to_copy);
  if (boundarydis_->NumGlobalNodes() == 0)
  {
    cout << "Empty boundary discretization detected. No FSI coupling will be performed..." << endl;
  }

  // sanity check (does not work for thin objects)
  vector< DRT::Condition * >      conditions;
  boundarydis_->GetCondition ("XFEMCoupling", conditions);
  const std::size_t numxfemcond = conditions.size();
  boundarydis_->GetCondition ("FSICoupling", conditions);
  const std::size_t numfsicond = conditions.size();
  if (numxfemcond != numfsicond)
    dserror("number of xfem conditions has to match number of fsi conditions");

  // remove internal surface elements that occur for flat 3D hex8 meshes
  RemoveInternalSurfElements(soliddis);
  // TODO: should we worry about internal surface elements within the soliddis?

  // create node and element distribution with elements and nodes ghosted on all processors
  const Epetra_Map noderowmap = *boundarydis_->NodeRowMap();
  const Epetra_Map elemrowmap = *boundarydis_->ElementRowMap();

  // put all boundary nodes and elements onto all processors
  const Epetra_Map nodecolmap = *LINALG::AllreduceEMap(noderowmap);
  const Epetra_Map elemcolmap = *LINALG::AllreduceEMap(elemrowmap);

  // redistribute nodes and elements to column (ghost) map
  boundarydis_->ExportColumnNodes(nodecolmap);
  boundarydis_->ExportColumnElements(elemcolmap);

  // Now we are done. :)
  const int err = boundarydis_->FillComplete();
  if (err) dserror("FillComplete() returned err=%d",err);

  DRT::UTILS::PrintParallelDistribution(*boundarydis_);

  boundaryoutput_ = rcp(new IO::DiscretizationWriter(boundarydis_));
  // mesh output is only needed for post processing in paraview
  // or if element date is read during restart
  // if some procs have no row elements, this WriteMesh() call fails currently
//  if (boundarydis_->NumMyRowElements() > 0)
//    boundaryoutput_->WriteMesh(0,0.0);

  interface_.Setup(*boundarydis_);
  fluid_.SetSurfaceSplitter(&interface_);

  // create interface DOF vectors using the solid parallel distribution
  idispnp_    = LINALG::CreateVector(*boundarydis_->DofRowMap(),true);
  ivelnp_     = LINALG::CreateVector(*boundarydis_->DofRowMap(),true);
  itrueresnp_ = LINALG::CreateVector(*boundarydis_->DofRowMap(),true);

  idispn_   = LINALG::CreateVector(*boundarydis_->DofRowMap(),true);
  iveln_    = LINALG::CreateVector(*boundarydis_->DofRowMap(),true);
  ivelnm_   = LINALG::CreateVector(*boundarydis_->DofRowMap(),true);
  iaccnp_   = LINALG::CreateVector(*boundarydis_->DofRowMap(),true);
  iaccn_    = LINALG::CreateVector(*boundarydis_->DofRowMap(),true);

  idispstepinc_ = LINALG::CreateVector(*boundarydis_->DofRowMap(),true);
  idispiterinc_ = LINALG::CreateVector(*boundarydis_->DofRowMap(),true);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluidImpl::InitialGuess()
{
  return fluid_.InitialGuess();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluidImpl::RHS()
{
  return fluid_.Residual();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluidImpl::TrueResidual()
{
  dserror("stop here");
  return fluid_.TrueResidual();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluidImpl::Velnp()
{
  return fluid_.Velnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluidImpl::Velaf()
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluidImpl::Veln()
{
  return fluid_.Veln();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluidImpl::Accam()
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluidImpl::Hist()
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluidImpl::Dispnp()
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluidImpl::ConvectiveVel()
{
  return Teuchos::null; // no moving mesh present
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::XFluidImpl::DofRowMap()
{
  Teuchos::RCP<const Epetra_Map> dofrowmap = rcp(fluid_.Discretization()->DofRowMap(),false);
  return dofrowmap;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::XFluidImpl::SystemMatrix()
{
  return fluid_.SystemMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<std::string, Teuchos::RCP<LINALG::SparseMatrix> > ADAPTER::XFluidImpl::CouplingMatrices()
{
  return fluid_.CouplingMatrices();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<std::string, Teuchos::RCP<Epetra_Vector> > ADAPTER::XFluidImpl::CouplingVectors()
{
  return fluid_.CouplingVectors();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase> ADAPTER::XFluidImpl::BlockSystemMatrix()
{
  dserror("no block matrix here");
  return fluid_.BlockSystemMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase> ADAPTER::XFluidImpl::ShapeDerivatives()
{
  dserror("nope");
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
Teuchos::RCP<const LINALG::MapExtractor> ADAPTER::XFluidImpl::GetDBCMapExtractor()
{
  return fluid_.DirichMaps();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Teuchos::RCP<Epetra_Vector> ADAPTER::XFluidImpl::StructCondRHS() const
// {
//   return interface_.ExtractCondVector(Velnp());
// }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::TimeLoop()
{
  dserror("Not implemented.");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::PrepareTimeStep()
{
  // update acceleration n
  iaccn_->Update(1.0,*iaccnp_,0.0);

  // update velocity n-1
  ivelnm_->Update(1.0,*iveln_,0.0);

  // update velocity n
  iveln_->Update(1.0,*ivelnp_,0.0);

  // update displacement n
  idispn_->Update(1.0,*idispnp_,0.0);

  idispstepinc_->PutScalar(0.0);
  idispiterinc_->PutScalar(0.0);

  ComputeInterfaceVelocities();
  ComputeInterfaceAccelerations();

  PrepareBoundaryDis();

  fluid_.PrepareTimeStep(boundarydis_);

  velpresstepinc_ = Teuchos::rcp(new Epetra_Vector(*Discretization()->DofRowMap(), true));
  velpresiterinc_ = Teuchos::rcp(new Epetra_Vector(*Discretization()->DofRowMap(), true));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::Evaluate(Teuchos::RCP<const Epetra_Vector> velpresstepinc)
{
  ComputeInterfaceVelocities();
  ComputeInterfaceAccelerations();

  PrepareBoundaryDis();

  if (idispiterinc_ == Teuchos::null)
    dserror("idispiterinc_ == Teuchos::null");
  boundarydis_->SetState("veliface nodal iterinc", idispiterinc_);

  // The field solver expects an iteration increment only. And
  // there are Dirichlet conditions that need to be preserved. So take
  // the sum of increments we get from NOX and apply the latest iteration
  // increment only.
  // Naming:
  //
  // x^n+1_i+1 = x^n+1_i + disiterinc  (sometimes referred to as residual increment), and
  //
  // x^n+1_i+1 = x^n     + disstepinc

  // reset iterative increment
  Teuchos::RCP<Epetra_Vector> velpresiterinc = Teuchos::null;

  if (velpresstepinc==Teuchos::null)
    dserror("velpresstepinc==Teuchos::null");

  DRT::DEBUGGING::NaNChecker(*velpresstepinc);

  if (velpresstepinc_==Teuchos::null)
    dserror("velpresstepinc_==Teuchos::null");
  if (velpresiterinc_==Teuchos::null)
    dserror("velpresiterinc_==Teuchos::null");

  // iteration increments
  if (not velpresstepinc_->Map().SameAs(velpresstepinc->Map()))  // global map size has changed
  {
    cout << "dummy" << endl;
    velpresiterinc_ = Teuchos::rcp(new Epetra_Vector(*fluid_.Discretization()->DofRowMap(), true));
    velpresstepinc_ = Teuchos::rcp(new Epetra_Vector(*velpresstepinc));
  }
  else // update stepinc and compute iterinc
  {
    cout << "real" << endl;
    velpresiterinc_->Update(1.0,*velpresstepinc,-1.0,*velpresstepinc_,0.0);

    // update incremental displacement member to provided step increments
    // shortly: disinc_^<i> := disp^<i+1>
    velpresstepinc_->Update(1.0,*velpresstepinc,0.0);
  }


  // do fluid evaluation provided residual displacements - iteration increment
  fluid_.Evaluate(boundarydis_, velpresiterinc_);

  itrueresnp_ = LINALG::CreateVector(*boundarydis_->DofRowMap(),true);
//  // get surface force
//  Teuchos::RCP<const Epetra_Vector> itruerescol = boundarydis_->GetState("iforcenp");
//
//  // dump all vectors in the boundary discretization
//  boundarydis_->ClearState();
//
//  // map back to solid parallel distribution
//  Teuchos::RCP<Epetra_Export> conimpo = Teuchos::rcp (new Epetra_Export(itruerescol->Map(),itrueresnp_->Map()));
//
//  itrueresnp_->PutScalar(0.0);
//  itrueresnp_->Export(*itruerescol,*conimpo,Add);

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::Update()
{
  fluid_.TimeUpdate();

  ComputeInterfaceVelocities();
  ComputeInterfaceAccelerations();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::ComputeInterfaceVelocities()
{
  // The displacement - velocity relation can be chosen freely and is used as defined
  // in the FSI parameter list (see "SECONDORDER")

  if (fluid_.TimIntScheme() == timeint_stationary)
  {
    ivelnp_->PutScalar(0.0);
    iveln_->PutScalar(0.0);
    ivelnm_->PutScalar(0.0);
    cout << "stationary FSI" << endl;
  }
  else
  {
    double thetaiface = 0.0;
    if (params_->get<bool>("interface second order"))  thetaiface = 0.5;
    else                                               thetaiface = 1.0;

    ivelnp_->Update(1.0/(thetaiface*Dt()),*idispnp_,-1.0/(thetaiface*Dt()),*idispn_,0.0);
    ivelnp_->Update(-(1.0-thetaiface)/thetaiface,*iveln_,1.0);

    cout << "instationary FSI" << endl;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::ComputeInterfaceAccelerations()
{
  // for interface acceleration, the SAME factor has to be used as for the fluid
  // time discretization. This way, we can extrapolate the fluid
  // material velocity and acceleration into the previously unknown
  // fluid domain.

  if (fluid_.TimIntScheme() == timeint_stationary)
  {
    iaccnp_->PutScalar(0.0);
    iaccn_->PutScalar(0.0);
  }
  else
  {
    const double dt = Dt();
    const double theta = fluid_.Theta();

    // compute acceleration at timestep n+1
    iaccnp_->Update(-(1.0-theta)/(theta),*iaccn_,0.0);
    iaccnp_->Update(1.0/(theta*dt),*ivelnp_,-1.0/(theta*dt),*iveln_,1.0);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::StatisticsAndOutput()
{
  fluid_.StatisticsAndOutput();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::Output()
{
  // now the interface output
  boundaryoutput_->NewStep(Step(),Time());
  boundaryoutput_->WriteVector("idispnp", idispnp_);
  boundaryoutput_->WriteVector("idispn", idispn_);
  boundaryoutput_->WriteVector("ivelnp", ivelnp_);
  boundaryoutput_->WriteVector("iveln", iveln_);
  boundaryoutput_->WriteVector("ivelnm", ivelnm_);
  boundaryoutput_->WriteVector("iaccnp", iaccnp_);
  boundaryoutput_->WriteVector("iaccn", iaccn_);
  boundaryoutput_->WriteVector("itrueresnp", itrueresnp_);

  // first fluid output
  fluid_.StatisticsAndOutput();

  // create interface DOF vectors using the fluid parallel distribution
  Teuchos::RCP<const Epetra_Vector> ivelnpcol   = DRT::UTILS::GetColVersionOfRowVector(boundarydis_, ivelnp_);
  Teuchos::RCP<const Epetra_Vector> ivelncol    = DRT::UTILS::GetColVersionOfRowVector(boundarydis_, iveln_);
  Teuchos::RCP<const Epetra_Vector> iaccnpcol   = DRT::UTILS::GetColVersionOfRowVector(boundarydis_, iaccnp_);
  Teuchos::RCP<const Epetra_Vector> iaccncol    = DRT::UTILS::GetColVersionOfRowVector(boundarydis_, iaccn_);
  Teuchos::RCP<const Epetra_Vector> idispnpcol  = DRT::UTILS::GetColVersionOfRowVector(boundarydis_, idispnp_);
  Teuchos::RCP<const Epetra_Vector> itruerescol = DRT::UTILS::GetColVersionOfRowVector(boundarydis_, itrueresnp_);

  // print redundant arrays on proc 0
  if (boundarydis_->Comm().MyPID() == 0)
  {
    PrintInterfaceVectorField(idispnpcol, ivelnpcol  , "sol_iface_velnp", "interface velocity n+1");
    PrintInterfaceVectorField(idispnpcol, ivelncol   , "sol_iface_veln" , "interface velocity n");
    PrintInterfaceVectorField(idispnpcol, iaccnpcol  , "sol_iface_accnp", "interface acceleration n+1");
    PrintInterfaceVectorField(idispnpcol, iaccncol   , "sol_iface_accn" , "interface acceleration n");
    PrintInterfaceVectorField(idispnpcol, idispnpcol , "sol_iface_disp" , "interface displacement n+1");
    PrintInterfaceVectorField(idispnpcol, itruerescol, "sol_iface_force", "interface force");
  }

  if (boundarydis_->Comm().MyPID() == 0 && itruerescol->MyLength() >= 3)
  {
    std::ofstream f;
    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                            + ".outifaceforce.txt";
    if (Step() <= 1)
      f.open(fname.c_str(),std::fstream::trunc);
    else
      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);

    f << Time() << " " << (*itruerescol)[0] << "  " << "\n";

    f.close();
  }

  if (boundarydis_->Comm().MyPID() == 0 && idispnpcol->MyLength() >= 3)
  {
    std::ofstream f;
    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                            + ".outifacedispnp.txt";
    if (Step() <= 1)
      f.open(fname.c_str(),std::fstream::trunc);
    else
      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);

    f << Time() << " " << (*idispnpcol)[0] << "  " << "\n";

    f.close();
  }

  if (boundarydis_->Comm().MyPID() == 0 && ivelnpcol->MyLength() >= 3)
  {
    std::ofstream f;
    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                            + ".outifacevelnp.txt";
    if (Step() <= 1)
      f.open(fname.c_str(),std::fstream::trunc);
    else
      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);

    f << Time() << " " << (*ivelnpcol)[0] << "  " << "\n";

    f.close();
  }

  if (boundarydis_->Comm().MyPID() == 0 && ivelnpcol->MyLength() >= 3)
  {
    std::ofstream f;
    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                            + ".outifaceveln.txt";
    if (Step() <= 1)
      f.open(fname.c_str(),std::fstream::trunc);
    else
      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);

    f << Time() << " " << (*ivelncol)[0] << "  " << "\n";

    f.close();
  }

  if (boundarydis_->Comm().MyPID() == 0 && iaccnpcol->MyLength() >= 3)
  {
    std::ofstream f;
    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                            + ".outifaceaccnp.txt";
    if (Step() <= 1)
      f.open(fname.c_str(),std::fstream::trunc);
    else
      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);

    f << Time() << " " << (*iaccnpcol)[0] << "  " << "\n";

    f.close();
  }

  if (boundarydis_->Comm().MyPID() == 0 && iaccncol->MyLength() >= 3)
  {
    std::ofstream f;
    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                            + ".outifaceaccn.txt";
    if (Step() <= 1)
      f.open(fname.c_str(),std::fstream::trunc);
    else
      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);

    f << Time() << " " << (*iaccncol)[0] << "  " << "\n";

    f.close();
  }
}


void ADAPTER::XFluidImpl::PrintInterfaceVectorField(
    const Teuchos::RCP<const Epetra_Vector>   displacementfield,
    const Teuchos::RCP<const Epetra_Vector>   vectorfield,
    const std::string& filestr,
    const std::string& name_in_gmsh
    ) const
{
  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT")==1;
  const bool screen_out = getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT_SCREEN")==1;
  if (gmshdebugout)
  {
    const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles(filestr, Step(), 5, screen_out, boundarydis_->Comm().MyPID());
    std::ofstream gmshfilecontent(filename.c_str());

    {
      gmshfilecontent << "View \" " << name_in_gmsh << " \" {\n";
      for (int i=0; i<boundarydis_->NumMyColElements(); ++i)
      {
        const DRT::Element* bele = boundarydis_->lColElement(i);
        vector<int> lm;
        vector<int> lmowner;
        bele->LocationVector(*boundarydis_, lm, lmowner);

        // extract local values from the global vector
        vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*vectorfield, myvelnp, lm);

        vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*displacementfield, mydisp, lm);

        const int nsd = 3;
        const int numnode = bele->NumNode();
        LINALG::SerialDenseMatrix elementvalues(nsd,numnode);
        LINALG::SerialDenseMatrix elementpositions(nsd,numnode);
        int counter = 0;
        for (int iparam=0; iparam<numnode; ++iparam)
        {
          const DRT::Node* node = bele->Nodes()[iparam];
          const double* pos = node->X();
          for (int isd = 0; isd < nsd; ++isd)
          {
            elementvalues(isd,iparam) = myvelnp[counter];
            elementpositions(isd,iparam) = pos[isd] + mydisp[counter];
            counter++;
          }
        }

        IO::GMSH::cellWithVectorFieldToStream(
            bele->Shape(), elementvalues, elementpositions, gmshfilecontent);
      }
      gmshfilecontent << "};\n";
    }
    gmshfilecontent.close();
    if (screen_out) std::cout << " done" << endl;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::PrepareBoundaryDis()
{
  // put vectors into boundary discretization (SetState generates col vector automatically)
  boundarydis_->SetState("idispcolnp",idispnp_);
  boundarydis_->SetState("idispcoln" ,idispn_);
  boundarydis_->SetState("ivelcolnp" ,ivelnp_);
  boundarydis_->SetState("ivelcoln"  ,iveln_);
  boundarydis_->SetState("ivelcolnm" ,ivelnm_);
  boundarydis_->SetState("iacccoln"  ,iaccn_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::NonlinearSolve()
{
  cout << "ADAPTER::XFluidImpl::NonlinearSolve()" << endl;

  ComputeInterfaceVelocities();
  ComputeInterfaceAccelerations();

  PrepareBoundaryDis();
  idispiterinc_ = rcp(new Epetra_Vector(idispnp_->Map(),true));
  boundarydis_->SetState("veliface nodal iterinc", idispiterinc_);

  // solve
  fluid_.NonlinearSolve(boundarydis_);

  // get surface force
  Teuchos::RCP<const Epetra_Vector> itruerescol = boundarydis_->GetState("iforcenp");

  // dump all vectors in the boundary discretization
  boundarydis_->ClearState();

  // map back to solid parallel distribution
  Teuchos::RCP<Epetra_Export> conimpo = Teuchos::rcp (new Epetra_Export(itruerescol->Map(),itrueresnp_->Map()));
  itrueresnp_->PutScalar(0.0);
  itrueresnp_->Export(*itruerescol,*conimpo,Add);

  if (TimIntScheme() == timeint_stationary)
  {
    LiftDrag();
  }

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::XFluidImpl::InnerVelocityRowMap()
{
  // build inner velocity map
  // dofs at the interface are excluded
  // we use only velocity dofs and only those without Dirichlet constraint
  const Teuchos::RCP<const LINALG::MapExtractor> dbcmaps = fluid_.DirichMaps();
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
  maps.push_back(interface_.OtherMap());
  maps.push_back(dbcmaps->OtherMap());
  innervelmap_ = LINALG::MultiMapExtractor::MergeMaps(maps);

  // deliver pizza
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
  if (params_->get<bool>("interface second order"))
  {
    return 2./fluid_.Dt();
  }
  else
    return 1./fluid_.Dt();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::ReadRestart(int step)
{
  IO::DiscretizationReader reader(boundarydis_,step);
  reader.ReadDouble("time");
  reader.ReadInt("step");
  // read interface position
  reader.ReadVector(idispnp_,"idispnp");
  reader.ReadVector(idispn_, "idispn");
  reader.ReadVector(ivelnp_, "ivelnp");
  reader.ReadVector(iveln_, "iveln");
  reader.ReadVector(ivelnm_, "ivelnm");
  reader.ReadVector(iaccn_, "iaccn");
  reader.ReadVector(itrueresnp_, "itrueresnp");

  // SetState generates colvectors automatically
  boundarydis_->SetState("idispcolnp",idispnp_);
  boundarydis_->SetState("idispcoln" ,idispn_);

  fluid_.ReadRestart(step,boundarydis_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::XFluidImpl::Time() const
{
  return fluid_.Time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int ADAPTER::XFluidImpl::Step() const
{
  return fluid_.Step();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::XFluidImpl::Dt() const
{
  return fluid_.Dt();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::LiftDrag()
{
  // get forces on all procs
  // create interface DOF vectors using the fluid parallel distribution
  Teuchos::RCP<const Epetra_Vector> iforcecol = DRT::UTILS::GetColVersionOfRowVector(boundarydis_, itrueresnp_);

  if (boundarydis_->Comm().MyPID() == 0)
  {
    // compute force components
    const int nsd = 3;
    const Epetra_Map* dofcolmap = boundarydis_->DofColMap();
    LINALG::Matrix<3,1> c(true);
    for (int inode = 0; inode < boundarydis_->NumMyColNodes(); ++inode)
    {
      const DRT::Node* node = boundarydis_->lColNode(inode);
      const std::vector<int> dof = boundarydis_->Dof(node);
      for (int isd = 0; isd < nsd; ++isd)
      {
        // minus to get correct sign of lift and drag (force acting on the body)
        c(isd) -= (*iforcecol)[dofcolmap->LID(dof[isd])];
      }
    }

    // print to file
    std::ostringstream s;
    std::ostringstream header;

    header << left  << std::setw(10) << "Time"
           << right << std::setw(16) << "F_x"
           << right << std::setw(16) << "F_y"
           << right << std::setw(16) << "F_z";
    s << left  << std::setw(10) << scientific << Time()
      << right << std::setw(16) << scientific << c(0)
      << right << std::setw(16) << scientific << c(1)
      << right << std::setw(16) << scientific << c(2);

    std::ofstream f;
    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
                            + ".liftdrag.txt";
    if (Step() <= 1)
    {
      f.open(fname.c_str(),std::fstream::trunc);
      //f << header.str() << endl;
    }
    else
    {
      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
    }
    f << s.str() << "\n";
    f.close();

//    std::cout << header.str() << endl << s.str() << endl;
  }

  fluid_.LiftDrag();
}

void ADAPTER::XFluidImpl::RemoveInternalSurfElements(
    const Teuchos::RCP<DRT::Discretization> soliddis
    )
{
//if (boundarydis_->Comm().NumProc() == 1)
{
  // for structures consisting of one layer of hex8 elements, internal surface elements appear
  // find internal elements
  std::set<int> internal_surfeles;
  for (int isurf=0; isurf<boundarydis_->NumMyRowElements(); ++isurf)
  {
    const DRT::Element* surfele = boundarydis_->lRowElement(isurf);

    LINALG::SerialDenseMatrix          xyze_surf(3,surfele->NumNode());
    GEO::fillInitialPositionArray(surfele,xyze_surf);

    // center in local coordinates
    const LINALG::Matrix<2,1> localcenterpos(DRT::UTILS::getLocalCenterPosition<2>(surfele->Shape()));
    // center in physical coordinates
    LINALG::Matrix<3,1> physicalcenterpos;
    GEO::elementToCurrentCoordinates(surfele->Shape(), xyze_surf, localcenterpos, physicalcenterpos);

    LINALG::Matrix<3,1> unitnormalvec;
    GEO::computeNormalToSurfaceElement(surfele->Shape(), xyze_surf, localcenterpos, unitnormalvec);

    for (int ivol=0; ivol<soliddis->NumMyColElements(); ++ivol)
    {
      const DRT::Element* solidele = soliddis->lColElement(ivol);

      // check if solid ele is connected with surf
      bool connected = false;
      for (int inode = 0; inode < surfele->NumNode(); ++inode)
        for (int jnode = 0; jnode < solidele->NumNode(); ++jnode)
          if (surfele->NodeIds()[inode] == solidele->NodeIds()[jnode])
          {
            connected = true;
            break;
          }


      if (connected)
      {
        LINALG::SerialDenseMatrix          solidxyze(3,solidele->NumNode());
        GEO::fillInitialPositionArray(solidele,solidxyze);
        const double volume = GEO::ElementVolume(solidele->Shape(), solidxyze);
        const double length = pow(volume,1.0/3.0);

        LINALG::Matrix<3,1> normalvector = unitnormalvec;
        normalvector.Scale(length/10.0);

        LINALG::Matrix<3,1> testpos(true);
        testpos += physicalcenterpos;
        testpos += normalvector;

        const bool inside = GEO::checkPositionWithinElement(solidele, solidxyze, testpos);
        if (inside)
        {
          internal_surfeles.insert(surfele->Id());
          break;
        }
      }
    }
  }

  if (not internal_surfeles.empty())
    cout << RED_LIGHT << "Found " << internal_surfeles.size() << " internal surface elements! Removing them..." << END_COLOR << endl;
  // remove internal elements
  for (std::set<int>::const_iterator ele=internal_surfeles.begin();
       ele!=internal_surfeles.end();
       ++ele)
  {
    boundarydis_->DeleteElement(*ele);
  }
  const int err2 = boundarydis_->FillComplete(false,true,true);
  if (err2) dserror("FillComplete() returned err=%d",err2);
}
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::XFluidImpl::ExtractInterfaceForces()
{
  return interface_.ExtractFSICondVector(itrueresnp_);
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
  // it depends, when this method is called, and when velnp is updated
  // the FSI algorithm expects first an time update and then asks for the old time step velocity
  // meaning that it gets the velocity from the new time step
  // not clear? exactly! thats why the FSI time update should be more clear about it
  // needs discussion with the FSI people
  return interface_.ExtractFSICondVector(ivelnp_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
{
  interface_.InsertFSICondVector(ivel,ivelnp_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::ApplyInterfaceRobinValue(Teuchos::RCP<Epetra_Vector> ivel, Teuchos::RCP<Epetra_Vector> iforce)
{
  dserror("robin unimplemented");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> idisp)
{
  interface_.InsertFSICondVector(idisp,idispnp_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::ApplyMeshDisplacementIncrement(Teuchos::RCP<const Epetra_Vector> idispstepinc)
{
  // The field solver always expects an iteration increment only. And
  // there are Dirichlet conditions that need to be preserved. So take
  // the sum of increments (disstepinc) and compute disiterinc
  // Naming:
  //
  // x^n+1_i+1 = x^n+1_i + disiterinc  (sometimes referred to as residual increment), and
  //
  // x^n+1_i+1 = x^n     + disstepinc

  if (idispstepinc == Teuchos::null) // first step
  {
    dserror("schimpf");
    idispnp_      = Teuchos::rcp(new Epetra_Vector(*idispn_));
    idispiterinc_ = Teuchos::rcp(new Epetra_Vector(idispn_->Map(),true));
    idispstepinc_ = Teuchos::rcp(new Epetra_Vector(idispn_->Map(),true));
  }
  else
  {
    idispnp_->Update(1.0, *idispn_, 1.0, *idispstepinc, 0.0);

    // iteration increments
    idispiterinc_ = Teuchos::rcp(new Epetra_Vector(*idispstepinc));
    idispiterinc_->Update(-1.0,*idispstepinc_,1.0);

    // update step increment
    // shortly: disinc_^<i> := stepincdisp^<i+1>
    idispstepinc_->Update(1.0,*idispstepinc,0.0);
  }
  DRT::DEBUGGING::NaNChecker(*idispn_);
  DRT::DEBUGGING::NaNChecker(*idispnp_);
  DRT::DEBUGGING::NaNChecker(*idispiterinc_);
  DRT::DEBUGGING::NaNChecker(*idispstepinc_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::ApplyMeshVelocity(Teuchos::RCP<const Epetra_Vector> gridvel)
{
  dserror("makes no sense here!");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx)
{
  dserror("not needed");
  // get interface velocity at t(n)
  const Teuchos::RCP<Epetra_Vector> veln = Interface().ExtractFSICondVector(Veln());

  // We convert Delta d(n+1,i+1) to Delta u(n+1,i+1) here.
  //
  // Delta d(n+1,i+1) = ( theta Delta u(n+1,i+1) + u(n) ) * dt
  //
  double timescale = TimeScaling();
  fcx->Update(-timescale*fluid_.Dt(),*veln,timescale);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::VelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx)
{
  dserror("not needed");
  // get interface velocity at t(n)
  const Teuchos::RCP<Epetra_Vector> veln = Interface().ExtractFSICondVector(Veln());

  // We convert Delta u(n+1,i+1) to Delta d(n+1,i+1) here.
  //
  // Delta d(n+1,i+1) = ( theta Delta u(n+1,i+1) + u(n) ) * dt
  //
  double timescale = 1./TimeScaling();
  fcx->Update(fluid_.Dt(),*veln,timescale);
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
  return interface_.ExtractFSICondVector(fluid_.IntegrateInterfaceShape("FSICoupling"));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::XFluidImpl::RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel)
{
  dserror("not implemented!");
  const Epetra_Map* dofrowmap = Discretization()->DofRowMap();
  Teuchos::RCP<Epetra_Vector> relax = LINALG::CreateVector(*dofrowmap,true);
  interface_.InsertFSICondVector(ivel,relax);
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
void ADAPTER::XFluidImpl::Predictor()
{
  dserror("not implemented!");
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::MultiCorrector()
{
  dserror("not implemented!");
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::SetInitialFlowField(int whichinitialfield,int startfuncno)
{
  fluid_.SetInitialFlowField(boundarydis_,whichinitialfield,startfuncno);
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::SetIterLomaFields(RCP<const Epetra_Vector> scalaraf,RCP<const Epetra_Vector> scalaram,RCP<const Epetra_Vector> scalardtam,const double thermpressaf,const double thermpressam,const double thermpressdtam,Teuchos::RCP<DRT::Discretization> scatradis)
{
   dserror("not implemented!");
   return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::SetTimeLomaFields(RCP<const Epetra_Vector> scalarnp,const double thermpressnp,RCP<const Epetra_Vector> scatraresidual,Teuchos::RCP<DRT::Discretization> scatradis)
{
   dserror("not implemented!");
   return;
}

#endif  // #ifdef CCADISCRET
