//*----------------------------------------------------------------------*/
/*!
\file adapter_xfluid2_impl.cpp

\brief XFluid field adapter

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "adapter_xfluid2_impl.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Epetra_Export.h>


#include "../drt_lib/drt_condition_utils.H"
#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_utils.H"
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
ADAPTER::XFluid2Impl::XFluid2Impl(
        Teuchos::RCP<DRT::Discretization> xfluiddis,
        Teuchos::RCP<DRT::Discretization> soliddis,
        Teuchos::RCP<LINALG::Solver> solver,
        Teuchos::RCP<ParameterList> params,
        Teuchos::RCP<IO::DiscretizationWriter> output )
  : fluid_(xfluiddis, soliddis, *solver, *params, output),
    xfluiddis_(xfluiddis),
    soliddis_(soliddis),
    solver_(solver),
    params_(params)//,
//    output_(output)
{
  //Assign after calling the Xfluid-constructor
  boundarydis_ = fluid_.Boundary_Dis();

  // the solid mesh has to match the interface mesh
  // so we have to compute a interface true residual vector itrueresidual_
  interface_.Setup(*boundarydis_);
  fluid_.SetSurfaceSplitter(&interface_);
//
//  itrueresnp_ = LINALG::CreateVector(*boundarydis_->DofRowMap(),true);


}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluid2Impl::Dispnp()
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluid2Impl::ConvectiveVel()
{
  return Teuchos::null; // no moving mesh present
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::XFluid2Impl::DofRowMap()
{
  Teuchos::RCP<const Epetra_Map> dofrowmap = rcp(fluid_.Discretization()->DofRowMap(),false);
  return dofrowmap;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::XFluid2Impl::SystemMatrix()
{
  return fluid_.SystemMatrix();
}


///*----------------------------------------------------------------------*/
///*----------------------------------------------------------------------*/
//std::map<std::string, Teuchos::RCP<LINALG::SparseMatrix> > ADAPTER::XFluid2Impl::CouplingMatrices()
//{
//  dserror("no coupling matrices here");
//
////  return fluid_.CouplingMatrices();
//}


///*----------------------------------------------------------------------*/
///*----------------------------------------------------------------------*/
//std::map<std::string, Teuchos::RCP<Epetra_Vector> > ADAPTER::XFluid2Impl::CouplingVectors()
//{
//  dserror("no coupling vectors here");
//
////  return fluid_.CouplingVectors();
//}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase> ADAPTER::XFluid2Impl::BlockSystemMatrix()
{
  dserror("no block matrix here");

  return fluid_.BlockSystemMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase> ADAPTER::XFluid2Impl::ShapeDerivatives()
{
  dserror("nope");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::XFluid2Impl::Discretization()
{
  // here we have to return the boundarydis for adapter_coupling
  // REMARK: coupling only possible between structure and boundarydis, these discretizations have matching nodes
  return boundarydis_ ;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const LINALG::MapExtractor> ADAPTER::XFluid2Impl::GetDBCMapExtractor()
{
  dserror("in XFluid2Impl::GetDBCMapExtractor");
  return fluid_.DirichMaps();
}


/*----------------------------------------------------------------------*/
// Integrate routine for STATIONARY INTERFACES
/*----------------------------------------------------------------------*/
void ADAPTER::XFluid2Impl::Integrate()
{

    // call standard integrate routine for stationary interfaces
    fluid_.Integrate();

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluid2Impl::TimeLoop()
{
  // integrate routine for Xfluid computations with stationary interfaces
  Integrate();

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluid2Impl::PrepareTimeStep()
{

//  cout << "XFluid2Impl::prepare time step" << endl;
//  dserror("XFluid2Impl::PrepareTimeStep");

//  // update acceleration n
//  fluid_.iaccn_->Update(1.0,*fluid_.iaccnp_,0.0);

  // update velocity n-1
  fluid_.ivelnm_->Update(1.0,*fluid_.iveln_,0.0);

  // update velocity n
  fluid_.iveln_->Update(1.0,*fluid_.ivelnp_,0.0);

  // update displacement n
  fluid_.idispn_->Update(1.0,*fluid_.idispnp_,0.0);
//
//  idispstepinc_->PutScalar(0.0);
//  idispiterinc_->PutScalar(0.0);
//
//  SetSolidDisplacement( Time() );
//
//  ComputeSolidVelocities();


//  ComputeInterfaceAccelerations();
//
//  PrepareBoundaryDis();




  fluid_.PrepareTimeStep();
//
//  velpresstepinc_ = Teuchos::rcp(new Epetra_Vector(*Discretization()->DofRowMap(), true));
//  velpresiterinc_ = Teuchos::rcp(new Epetra_Vector(*Discretization()->DofRowMap(), true));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluid2Impl::Evaluate(Teuchos::RCP<const Epetra_Vector> velpresstepinc)
{
  dserror("XFluid2Impl::Evaluate");
//
//  ComputeInterfaceVelocities();
//  ComputeInterfaceAccelerations();
//
//  PrepareBoundaryDis();
//
//  if (idispiterinc_ == Teuchos::null)
//    dserror("idispiterinc_ == Teuchos::null");
//  fluid_.Boundary_Dis()->SetState("veliface nodal iterinc", idispiterinc_);
//
//  // The field solver expects an iteration increment only. And
//  // there are Dirichlet conditions that need to be preserved. So take
//  // the sum of increments we get from NOX and apply the latest iteration
//  // increment only.
//  // Naming:
//  //
//  // x^n+1_i+1 = x^n+1_i + disiterinc  (sometimes referred to as residual increment), and
//  //
//  // x^n+1_i+1 = x^n     + disstepinc
//
//  // reset iterative increment
//  Teuchos::RCP<Epetra_Vector> velpresiterinc = Teuchos::null;
//
//  if (velpresstepinc==Teuchos::null)
//    dserror("velpresstepinc==Teuchos::null");
//
//  DRT::DEBUGGING::NaNChecker(*velpresstepinc);
//
//  if (velpresstepinc_==Teuchos::null)
//    dserror("velpresstepinc_==Teuchos::null");
//  if (velpresiterinc_==Teuchos::null)
//    dserror("velpresiterinc_==Teuchos::null");
//
//  // iteration increments
//  if (not velpresstepinc_->Map().SameAs(velpresstepinc->Map()))  // global map size has changed
//  {
//    cout << "dummy" << endl;
//    velpresiterinc_ = Teuchos::rcp(new Epetra_Vector(*fluid_.Discretization()->DofRowMap(), true));
//    velpresstepinc_ = Teuchos::rcp(new Epetra_Vector(*velpresstepinc));
//  }
//  else // update stepinc and compute iterinc
//  {
//    cout << "real" << endl;
//    velpresiterinc_->Update(1.0,*velpresstepinc,-1.0,*velpresstepinc_,0.0);
//
//    // update incremental displacement member to provided step increments
//    // shortly: disinc_^<i> := disp^<i+1>
//    velpresstepinc_->Update(1.0,*velpresstepinc,0.0);
//  }
//
//
//  // do fluid evaluation provided residual displacements - iteration increment
////  fluid_.Evaluate(boundarydis_, velpresiterinc_);
//
//  itrueresnp_ = LINALG::CreateVector(*fluid_.Boundary_Dis()->DofRowMap(),true);
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
void ADAPTER::XFluid2Impl::Update()
{

  cout << "ADAPTER::XFluid2Impl::Update" << endl;

  // call fluid time update
  fluid_.TimeUpdate();

}




/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluid2Impl::ComputeInterfaceVelocities()
{

  INPAR::XFEM::MovingBoundary xfluid_mov_bound = DRT::INPUT::get<INPAR::XFEM::MovingBoundary>(params_->sublist("XFEM"), "XFLUID_BOUNDARY");

  // -------------------------------------------------------------------
  //       update interface (solid) displacement and velocity
  // -------------------------------------------------------------------
  if(xfluid_mov_bound == INPAR::XFEM::XFluidMovingBoundary)
  {
    // set interface displacement for moving boundary cases in XFluid, not in Xfsi
    fluid_.SetInterfaceDisplacement( Time() );
  }
  else if(xfluid_mov_bound == INPAR::XFEM::XFluidStationaryBoundary)
  {
    dserror("stationary boundary should be called via Xfluid->Integrate ...");
  }

  // compute the interface velocity for all moving boundary cases
  fluid_.ComputeInterfaceVelocities();

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluid2Impl::ComputeInterfaceAccelerations()
{
  cout << "XFluid2Impl::ComputeInterfaceAcceleration" << " -> Do we need this function???" << endl;


  // for interface acceleration, the SAME factor has to be used as for the fluid
  // time discretization. This way, we can extrapolate the fluid
  // material velocity and acceleration into the previously unknown
  // fluid domain.

//  if (fluid_.TimIntScheme() == INPAR::FLUID::timeint_stationary)
//  {
//    iaccnp_->PutScalar(0.0);
//    iaccn_->PutScalar(0.0);
//  }
//  else
//  {
//    const double dt = Dt();
//    const double theta = fluid_.Theta();
//
//    // compute acceleration at timestep n+1
//    iaccnp_->Update(-(1.0-theta)/(theta),*iaccn_,0.0);
//    iaccnp_->Update(1.0/(theta*dt),*ivelnp_,-1.0/(theta*dt),*iveln_,1.0);
//
//    // note: for BDF2, the acceleration is ignored, so no need for another if here
//  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluid2Impl::StatisticsAndOutput()
{
  fluid_.StatisticsAndOutput();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluid2Impl::Output()
{

  fluid_.StatisticsAndOutput();

//  dserror("XFluid2Impl::Output");

//
//  // now the interface output
//  boundaryoutput_->NewStep(Step(),Time());
//  boundaryoutput_->WriteVector("idispnp", idispnp_);
//  boundaryoutput_->WriteVector("idispn", idispn_);
//  boundaryoutput_->WriteVector("ivelnp", ivelnp_);
//  boundaryoutput_->WriteVector("iveln", iveln_);
//  boundaryoutput_->WriteVector("ivelnm", ivelnm_);
//  boundaryoutput_->WriteVector("iaccnp", iaccnp_);
//  boundaryoutput_->WriteVector("iaccn", iaccn_);
//  boundaryoutput_->WriteVector("itrueresnp", itrueresnp_);
//
//  // first fluid output
//  fluid_.StatisticsAndOutput();
//
//  // create interface DOF vectors using the fluid parallel distribution
//  Teuchos::RCP<const Epetra_Vector> ivelnpcol   = DRT::UTILS::GetColVersionOfRowVector(boundarydis_, ivelnp_);
//  Teuchos::RCP<const Epetra_Vector> ivelncol    = DRT::UTILS::GetColVersionOfRowVector(boundarydis_, iveln_);
//  Teuchos::RCP<const Epetra_Vector> iaccnpcol   = DRT::UTILS::GetColVersionOfRowVector(boundarydis_, iaccnp_);
//  Teuchos::RCP<const Epetra_Vector> iaccncol    = DRT::UTILS::GetColVersionOfRowVector(boundarydis_, iaccn_);
//  Teuchos::RCP<const Epetra_Vector> idispnpcol  = DRT::UTILS::GetColVersionOfRowVector(boundarydis_, idispnp_);
//  Teuchos::RCP<const Epetra_Vector> itruerescol = DRT::UTILS::GetColVersionOfRowVector(boundarydis_, itrueresnp_);
//
//  // print redundant arrays on proc 0
//  if (boundarydis_->Comm().MyPID() == 0)
//  {
//    PrintInterfaceVectorField(idispnpcol, ivelnpcol  , "sol_iface_velnp", "interface velocity n+1");
//    PrintInterfaceVectorField(idispnpcol, ivelncol   , "sol_iface_veln" , "interface velocity n");
//    PrintInterfaceVectorField(idispnpcol, iaccnpcol  , "sol_iface_accnp", "interface acceleration n+1");
//    PrintInterfaceVectorField(idispnpcol, iaccncol   , "sol_iface_accn" , "interface acceleration n");
//    PrintInterfaceVectorField(idispnpcol, idispnpcol , "sol_iface_disp" , "interface displacement n+1");
//    PrintInterfaceVectorField(idispnpcol, itruerescol, "sol_iface_force", "interface force");
//  }
//
//  if (boundarydis_->Comm().MyPID() == 0 && itruerescol->MyLength() >= 3)
//  {
//    std::ofstream f;
//    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
//                            + ".outifaceforce.txt";
//    if (Step() <= 1)
//      f.open(fname.c_str(),std::fstream::trunc);
//    else
//      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
//
//    f << Time() << " " << (*itruerescol)[0] << "  " << "\n";
//
//    f.close();
//  }
//
//  if (boundarydis_->Comm().MyPID() == 0 && idispnpcol->MyLength() >= 3)
//  {
//    std::ofstream f;
//    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
//                            + ".outifacedispnp.txt";
//    if (Step() <= 1)
//      f.open(fname.c_str(),std::fstream::trunc);
//    else
//      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
//
//    f << Time() << " " << (*idispnpcol)[0] << "  " << "\n";
//
//    f.close();
//  }
//
//  if (boundarydis_->Comm().MyPID() == 0 && ivelnpcol->MyLength() >= 3)
//  {
//    std::ofstream f;
//    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
//                            + ".outifacevelnp.txt";
//    if (Step() <= 1)
//      f.open(fname.c_str(),std::fstream::trunc);
//    else
//      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
//
//    f << Time() << " " << (*ivelnpcol)[0] << "  " << "\n";
//
//    f.close();
//  }
//
//  if (boundarydis_->Comm().MyPID() == 0 && ivelnpcol->MyLength() >= 3)
//  {
//    std::ofstream f;
//    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
//                            + ".outifaceveln.txt";
//    if (Step() <= 1)
//      f.open(fname.c_str(),std::fstream::trunc);
//    else
//      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
//
//    f << Time() << " " << (*ivelncol)[0] << "  " << "\n";
//
//    f.close();
//  }
//
//  if (boundarydis_->Comm().MyPID() == 0 && iaccnpcol->MyLength() >= 3)
//  {
//    std::ofstream f;
//    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
//                            + ".outifaceaccnp.txt";
//    if (Step() <= 1)
//      f.open(fname.c_str(),std::fstream::trunc);
//    else
//      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
//
//    f << Time() << " " << (*iaccnpcol)[0] << "  " << "\n";
//
//    f.close();
//  }
//
//  if (boundarydis_->Comm().MyPID() == 0 && iaccncol->MyLength() >= 3)
//  {
//    std::ofstream f;
//    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
//                            + ".outifaceaccn.txt";
//    if (Step() <= 1)
//      f.open(fname.c_str(),std::fstream::trunc);
//    else
//      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
//
//    f << Time() << " " << (*iaccncol)[0] << "  " << "\n";
//
//    f.close();
//  }
}


void ADAPTER::XFluid2Impl::PrintInterfaceVectorField(
    const Teuchos::RCP<const Epetra_Vector>   displacementfield,
    const Teuchos::RCP<const Epetra_Vector>   vectorfield,
    const std::string& filestr,
    const std::string& name_in_gmsh
    ) const
{
  dserror("XFluid2Impl::PrintInterfaceVectorField");


  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = DRT::INPUT::IntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT")==1;
  const bool screen_out = DRT::INPUT::IntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT_SCREEN")==1;
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
        vector<int> lmstride;
        bele->LocationVector(*boundarydis_, lm, lmowner,lmstride);

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
void ADAPTER::XFluid2Impl::PrepareBoundaryDis()
{
  dserror("XFluid2Impl::PrepareBoundaryDis");

//  // put vectors into boundary discretization (SetState generates col vector automatically)
//  boundarydis_->SetState("idispcolnp",idispnp_);
//  boundarydis_->SetState("idispcoln" ,idispn_);
//  boundarydis_->SetState("ivelcolnp" ,ivelnp_);
//  boundarydis_->SetState("ivelcoln"  ,iveln_);
//  boundarydis_->SetState("ivelcolnm" ,ivelnm_);
//  boundarydis_->SetState("iacccoln"  ,iaccn_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluid2Impl::NonlinearSolve()
{
  // compute or set interface velocities just for XFluidMovingBoundary
  // REMARK: for XFSI this is done by ApplyMeshDisplacement and ApplyInterfaceVelocities
  INPAR::XFEM::MovingBoundary xfluid_mov_bound = DRT::INPUT::get<INPAR::XFEM::MovingBoundary>(params_->sublist("XFEM"), "XFLUID_BOUNDARY");
  if( xfluid_mov_bound == INPAR::XFEM::XFluidMovingBoundary)
  {
    ComputeInterfaceVelocities();
    ComputeInterfaceAccelerations();
  }


  // prepare nonlinear solve (includes cut, compute values for newly created nodes, set old part of righthandside, Neumann
  fluid_.PrepareNonlinearSolve();

  // solve
  fluid_.NonlinearSolve();


}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::XFluid2Impl::InnerVelocityRowMap()
{
  dserror("ADAPTER::XFluid2Impl::InnerVelocityRowMap");

//  // build inner velocity map
//  // dofs at the interface are excluded
//  // we use only velocity dofs and only those without Dirichlet constraint
//  const Teuchos::RCP<const LINALG::MapExtractor> dbcmaps = fluid_.DirichMaps();
//  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
//  maps.push_back(interface_.OtherMap());
//  maps.push_back(dbcmaps->OtherMap());
//  innervelmap_ = LINALG::MultiMapExtractor::MergeMaps(maps);
//
//  // deliver pizza
//  return innervelmap_;

  return null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::XFluid2Impl::VelocityRowMap()
{
  return fluid_.VelocityRowMap();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::XFluid2Impl::PressureRowMap()
{
  return fluid_.PressureRowMap();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::XFluid2Impl::SetMeshMap(Teuchos::RCP<const Epetra_Map> mm)
{
  dserror("makes no sense here");
  //meshmap_.Setup(*dis_->DofRowMap(),mm,LINALG::SplitMap(*dis_->DofRowMap(),*mm));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::XFluid2Impl::ResidualScaling() const
{
  return fluid_.ResidualScaling();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::XFluid2Impl::TimeScaling() const
{
  dserror("XFluid2Impl::TimeScaling");


//  if (params_->get<bool>("interface second order"))
//  {
//    return 2./fluid_.Dt();
//  }
//  else
//    return 1./fluid_.Dt();
  return 1.0;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluid2Impl::ReadRestart(int step)
{
  dserror("XFluid2Impl::ReadRestart");

//  IO::DiscretizationReader reader(boundarydis_,step);
//  reader.ReadDouble("time");
//  reader.ReadInt("step");
//  // read interface position
//  reader.ReadVector(idispnp_,"idispnp");
//  reader.ReadVector(idispn_, "idispn");
//  reader.ReadVector(ivelnp_, "ivelnp");
//  reader.ReadVector(iveln_, "iveln");
//  reader.ReadVector(ivelnm_, "ivelnm");
//  reader.ReadVector(iaccn_, "iaccn");
//  reader.ReadVector(itrueresnp_, "itrueresnp");
//
//  // SetState generates colvectors automatically
//  boundarydis_->SetState("idispcolnp",idispnp_);
//  boundarydis_->SetState("idispcoln" ,idispn_);
//
//  fluid_.ReadRestart(step,boundarydis_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::XFluid2Impl::Time() const
{
  return fluid_.Time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int ADAPTER::XFluid2Impl::Step() const
{
  return fluid_.Step();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::XFluid2Impl::Dt() const
{
  return fluid_.Dt();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluid2Impl::LiftDrag()
{
  fluid_.LiftDrag();
}

void ADAPTER::XFluid2Impl::RemoveInternalSurfElements(
    const Teuchos::RCP<DRT::Discretization> soliddis
    )
{
  dserror("XFluid2Impl::RemoveInternalSurfElements");
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
Teuchos::RCP<Epetra_Vector> ADAPTER::XFluid2Impl::ExtractInterfaceForces()
{

  cout << "ExtractInterfaceForces (itrueresnp)" << endl;

  // the trueresidual vector has to match the solid dis
  cout << "itrueresidual_" << *fluid_.itrueresidual_ << endl;

  return interface_.ExtractFSICondVector(fluid_.itrueresidual_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::XFluid2Impl::ExtractInterfaceForcesRobin()
{
  dserror("no Robin coupling here");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::XFluid2Impl::ExtractInterfaceFluidVelocity()
{
  dserror("no Robin coupling here");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::XFluid2Impl::ExtractInterfaceVeln()
{
  cout << "call ExtractInterfaceVeln() "<< endl;

  // it depends, when this method is called, and when velnp is updated
  // the FSI algorithm expects first an time update and then asks for the old time step velocity
  // meaning that it gets the velocity from the new time step
  // not clear? exactly! thats why the FSI time update should be more clear about it
  // needs discussion with the FSI people
  return interface_.ExtractFSICondVector(fluid_.ivelnp_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluid2Impl::ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
{
  cout << "ApplyInterfaceVelocities" << endl;
  interface_.InsertFSICondVector(ivel,fluid_.ivelnp_);

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluid2Impl::ApplyInterfaceRobinValue(Teuchos::RCP<Epetra_Vector> ivel, Teuchos::RCP<Epetra_Vector> iforce)
{
  dserror("XFluid2Impl::ApplyInterfaceRobinValue");
  dserror("robin unimplemented");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluid2Impl::ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> idisp)
{
  cout << "ApplyMeshDisplacement" << endl;
  interface_.InsertFSICondVector(idisp,fluid_.idispnp_);

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluid2Impl::ApplyMeshDisplacementIncrement(Teuchos::RCP<const Epetra_Vector> idispstepinc)
{
  dserror("XFluid2Impl::ApplyMeshDisplacementIncrement");
//
//  // The field solver always expects an iteration increment only. And
//  // there are Dirichlet conditions that need to be preserved. So take
//  // the sum of increments (disstepinc) and compute disiterinc
//  // Naming:
//  //
//  // x^n+1_i+1 = x^n+1_i + disiterinc  (sometimes referred to as residual increment), and
//  //
//  // x^n+1_i+1 = x^n     + disstepinc
//
//  if (idispstepinc == Teuchos::null) // first step
//  {
//    dserror("schimpf");
//    idispnp_      = Teuchos::rcp(new Epetra_Vector(*idispn_));
//    idispiterinc_ = Teuchos::rcp(new Epetra_Vector(idispn_->Map(),true));
//    idispstepinc_ = Teuchos::rcp(new Epetra_Vector(idispn_->Map(),true));
//  }
//  else
//  {
//    idispnp_->Update(1.0, *idispn_, 1.0, *idispstepinc, 0.0);
//
//    // iteration increments
//    idispiterinc_ = Teuchos::rcp(new Epetra_Vector(*idispstepinc));
//    idispiterinc_->Update(-1.0,*idispstepinc_,1.0);
//
//    // update step increment
//    // shortly: disinc_^<i> := stepincdisp^<i+1>
//    idispstepinc_->Update(1.0,*idispstepinc,0.0);
//  }
//  DRT::DEBUGGING::NaNChecker(*idispn_);
//  DRT::DEBUGGING::NaNChecker(*idispnp_);
//  DRT::DEBUGGING::NaNChecker(*idispiterinc_);
//  DRT::DEBUGGING::NaNChecker(*idispstepinc_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluid2Impl::ApplyMeshVelocity(Teuchos::RCP<const Epetra_Vector> gridvel)
{
  dserror("makes no sense here!");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluid2Impl::DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx)
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
void ADAPTER::XFluid2Impl::VelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx)
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
int ADAPTER::XFluid2Impl::Itemax() const
{
  return fluid_.Itemax();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluid2Impl::SetItemax(int itemax)
{
  fluid_.SetItemax(itemax);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::XFluid2Impl::IntegrateInterfaceShape()
{
  //return interface_.ExtractFSICondVector(fluid_.IntegrateInterfaceShape("FSICoupling"));
  dserror("IntegrateInterfaceShape not implemented!");
  return null;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::XFluid2Impl::RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel)
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
Teuchos::RCP<DRT::ResultTest> ADAPTER::XFluid2Impl::CreateFieldTest()
{
  return Teuchos::rcp(new FLD::XFluidResultTest2(&fluid_));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluid2Impl::ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres)
{
  dserror("not implemented!");
  return (fluid_.VelPresSplitter())->ExtractOtherVector(velpres);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluid2Impl::Predictor()
{
  dserror("not implemented!");
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluid2Impl::MultiCorrector()
{
  dserror("not implemented!");
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluid2Impl::SetInitialFlowField(const INPAR::FLUID::InitialField initfield,const int startfuncno)
{
  fluid_.SetInitialFlowField(initfield,startfuncno);
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluid2Impl::SetIterLomaFields(RCP<const Epetra_Vector> scalaraf,RCP<const Epetra_Vector> scalaram,RCP<const Epetra_Vector> scalardtam,RCP<const Epetra_Vector> fsscalaraf,const double thermpressaf,const double thermpressam,const double thermpressdtaf,const double thermpressdtam,Teuchos::RCP<DRT::Discretization> scatradis)
{
   dserror("not implemented!");
   return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluid2Impl::SetTimeLomaFields(RCP<const Epetra_Vector> scalarnp,const double thermpressnp,RCP<const Epetra_Vector> scatraresidual,Teuchos::RCP<DRT::Discretization> scatradis, const int whichscalar)
{
   dserror("not implemented!");
   return;
}

#endif  // #ifdef CCADISCRET
