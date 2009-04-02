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
#include "../drt_io/io_control.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_geometry/position_array.H"
#include "../drt_geometry/element_volume.H"
#include "../drt_geometry/intersection_service.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Epetra_Export.h>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::XFluidImpl::XFluidImpl(
        Teuchos::RCP<DRT::Discretization> dis,
        const Teuchos::RCP<DRT::Discretization> soliddis,
        Teuchos::RCP<ParameterList> params)
  : fluid_(dis, *params),
    dis_(dis),
    params_(params)
{
  vector<string> conditions_to_copy;
  conditions_to_copy.push_back("FSICoupling");
  conditions_to_copy.push_back("XFEMCoupling");
  boundarydis_ = DRT::UTILS::CreateDiscretizationFromCondition(soliddis, "FSICoupling", "boundary", "BELE3", conditions_to_copy);
  dsassert(boundarydis_->NumGlobalNodes() > 0, "empty discretization detected. FSICoupling condition applied?");
  
  // sanity check
  vector< DRT::Condition * >      conditions;
  boundarydis_->GetCondition ("XFEMCoupling", conditions);
  const unsigned numxfemcond = conditions.size();
  boundarydis_->GetCondition ("FSICoupling", conditions);
  const unsigned numfsicond = conditions.size();
  if (numxfemcond != numfsicond)
    dserror("number of xfem conditions has to match number of fsi conditions");
  
  // remove internal surface elements that occur for flat 3D hex8 meshes
  RemoveInternalSurfElements(soliddis);

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
  
  boundaryoutput_ = rcp(new IO::DiscretizationWriter(boundarydis_));
  boundaryoutput_->WriteMesh(0,0.0);

  DRT::UTILS::SetupNDimExtractor(*boundarydis_,"FSICoupling",interface_);
  DRT::UTILS::SetupNDimExtractor(*boundarydis_,"FREESURFCoupling",freesurface_);

  // create interface DOF vectors using the solid parallel distribution
  idispnp_    = LINALG::CreateVector(*boundarydis_->DofRowMap(),true);
  ivelnp_     = LINALG::CreateVector(*boundarydis_->DofRowMap(),true);
  itrueresnp_ = LINALG::CreateVector(*boundarydis_->DofRowMap(),true);

  idispn_   = LINALG::CreateVector(*boundarydis_->DofRowMap(),true);
  iveln_    = LINALG::CreateVector(*boundarydis_->DofRowMap(),true);
  ivelnm_   = LINALG::CreateVector(*boundarydis_->DofRowMap(),true);
  iaccn_    = LINALG::CreateVector(*boundarydis_->DofRowMap(),true);

  fluid_.SetFreeSurface(&freesurface_);
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
  return fluid_.SystemMatrix();
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
Teuchos::RCP<LINALG::BlockSparseMatrixBase> ADAPTER::XFluidImpl::MeshMoveMatrix()
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

  // update velocity n-1
  ivelnm_->Update(1.0,*iveln_,0.0);
  
  // update velocity n
  iveln_->Update(1.0,*ivelnp_,0.0);
  
  // update displacement n
  idispn_->Update(1.0,*idispnp_,0.0);
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
  // first fluid output
  fluid_.StatisticsAndOutput();
  
  // now the interface output
  boundaryoutput_->NewStep(Step(),Time());
  boundaryoutput_->WriteVector("idispnp", idispnp_);
  boundaryoutput_->WriteVector("idispn", idispn_);
  boundaryoutput_->WriteVector("ivelnp", ivelnp_);
  boundaryoutput_->WriteVector("iveln", iveln_);
  boundaryoutput_->WriteVector("ivelnm", ivelnm_);
  boundaryoutput_->WriteVector("iaccn", iaccn_);
//  boundaryoutput_->WriteVector("interface force", itrueres_);

  // now interface gmsh output
  
  // create interface DOF vectors using the fluid parallel distribution
  Teuchos::RCP<Epetra_Vector> ivelnpcol   = LINALG::CreateVector(*boundarydis_->DofColMap(),true);
//  Teuchos::RCP<Epetra_Vector> ivelncol    = LINALG::CreateVector(*boundarydis_->DofColMap(),true);
//  Teuchos::RCP<Epetra_Vector> iaccncol    = LINALG::CreateVector(*boundarydis_->DofColMap(),true);
  Teuchos::RCP<Epetra_Vector> idispnpcol  = LINALG::CreateVector(*boundarydis_->DofColMap(),true);
  Teuchos::RCP<Epetra_Vector> itruerescol = LINALG::CreateVector(*boundarydis_->DofColMap(),true);
  
  // map to fluid parallel distribution
  LINALG::Export(*idispnp_ ,*idispnpcol);
  LINALG::Export(*ivelnp_  ,*ivelnpcol);
//  LINALG::Export(*iveln_   ,*ivelncol);
//  LINALG::Export(*iaccn_   ,*iaccncol);
  LINALG::Export(*itrueresnp_,*itruerescol);
  
  // print redundant arrays on proc 0
  if (boundarydis_->Comm().MyPID() == 0)
  {
    PrintInterfaceVectorField(idispnpcol, itruerescol, ".solution_iforce_", "interface force");
    PrintInterfaceVectorField(idispnpcol, ivelnpcol, ".solution_ivel_"  , "interface velocity n+1");
    PrintInterfaceVectorField(idispnpcol, idispnpcol, ".solution_idisp_"  , "interface displacement n+1");
//    PrintInterfaceVectorField(idispnpcol, ivelncol , ".solution_iveln_" , "interface velocity n");
//    PrintInterfaceVectorField(idispnpcol, iaccncol , ".solution_iaccn_" , "interface acceleration n");
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

}

void ADAPTER::XFluidImpl::PrintInterfaceVectorField(
    const Teuchos::RCP<Epetra_Vector>   displacementfield,
    const Teuchos::RCP<Epetra_Vector>   vectorfield,
    const std::string filestr,
    const std::string name_in_gmsh
    ) const
{
  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = (bool)getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT");
  const bool screen_out = false;
  if (gmshdebugout)
  {
    std::stringstream filename;
    std::stringstream filenamedel;
    const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileName();
    filename << filebase << filestr << std::setw(5) << setfill('0') << Step() << ".pos";
    filenamedel << filebase << filestr << std::setw(5) << setfill('0') << Step()-5 << ".pos";
    std::remove(filenamedel.str().c_str());
    if (screen_out) std::cout << "writing " << left << std::setw(50) <<filename.str()<<"...";
    std::ofstream f_system(filename.str().c_str());
    
    {
      stringstream gmshfilecontent;
      gmshfilecontent << "View \" " << name_in_gmsh << " \" {\n";
      for (int i=0; i<boundarydis_->NumMyColElements(); ++i)
      {
        const DRT::Element* actele = boundarydis_->lColElement(i);
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
        LINALG::SerialDenseMatrix elementvalues(nsd,numnode);
        LINALG::SerialDenseMatrix elementpositions(nsd,numnode);
        int counter = 0;
        for (int iparam=0; iparam<numnode; ++iparam)
        {
          const DRT::Node* node = actele->Nodes()[iparam];
          const double* pos = node->X();
          for (int isd = 0; isd < nsd; ++isd)
          {
            elementvalues(isd,iparam) = myvelnp[counter];
            elementpositions(isd,iparam) = pos[isd] + mydisp[counter];
            counter++;
          }
        }
        
        gmshfilecontent << IO::GMSH::cellWithVectorFieldToString(
            actele->Shape(), elementvalues, elementpositions) << "\n";
      }
      gmshfilecontent << "};\n";
      f_system << gmshfilecontent.str();
    }
    f_system.close();
    if (screen_out) std::cout << " done" << endl;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::NonlinearSolve()
{

  // create interface DOF vectors using the fluid parallel distribution
  Teuchos::RCP<Epetra_Vector> idispcolnp = LINALG::CreateVector(*boundarydis_->DofColMap(),true);
  Teuchos::RCP<Epetra_Vector> idispcoln  = LINALG::CreateVector(*boundarydis_->DofColMap(),true);
  Teuchos::RCP<Epetra_Vector> ivelcolnp  = LINALG::CreateVector(*boundarydis_->DofColMap(),true);
  Teuchos::RCP<Epetra_Vector> ivelcoln   = LINALG::CreateVector(*boundarydis_->DofColMap(),true);
  Teuchos::RCP<Epetra_Vector> ivelcolnm  = LINALG::CreateVector(*boundarydis_->DofColMap(),true);
  Teuchos::RCP<Epetra_Vector> iacccoln   = LINALG::CreateVector(*boundarydis_->DofColMap(),true);
//  Teuchos::RCP<Epetra_Vector> itruerescol = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
//
//  Teuchos::RCP<Epetra_Vector> ivelncol = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
//  Teuchos::RCP<Epetra_Vector> iaccncol = LINALG::CreateVector(*fluidsurface_dofcolmap,true);

  // map to fluid parallel distribution
  LINALG::Export(*idispnp_,*idispcolnp);
  LINALG::Export(*idispn_ ,*idispcoln);
  LINALG::Export(*ivelnp_ ,*ivelcolnp);
  LINALG::Export(*iveln_  ,*ivelcoln);
  LINALG::Export(*ivelnm_ ,*ivelcolnm);
  LINALG::Export(*iaccn_  ,*iacccoln);


  boundarydis_->SetState("idispcolnp",idispcolnp);
  boundarydis_->SetState("idispcoln" ,idispcoln);
  
  boundarydis_->SetState("ivelcolnp" ,ivelcolnp);
  boundarydis_->SetState("ivelcoln"  ,ivelcoln);
  boundarydis_->SetState("ivelcolnm" ,ivelcolnm);
  boundarydis_->SetState("iacccoln"  ,iacccoln);
  
  fluid_.NonlinearSolve(boundarydis_);
  
  Teuchos::RCP<const Epetra_Vector> itruerescol = boundarydis_->GetState("iforcenp");
  
  boundarydis_->ClearState();

  // map back to solid parallel distribution
  Teuchos::RCP<Epetra_Export> conimpo = Teuchos::rcp (new Epetra_Export(itruerescol->Map(),itrueresnp_->Map()));
  itrueresnp_->PutScalar(0.0);
  itrueresnp_->Export(*itruerescol,*conimpo,Add); 
  //LINALG::Export(*itruerescol,*itrueresnp_);
  
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
  // create interface DOF vectors using the fluid parallel distribution
  Teuchos::RCP<Epetra_Vector> idispcolnp = LINALG::CreateVector(*boundarydis_->DofColMap(),true);
  Teuchos::RCP<Epetra_Vector> idispcoln  = LINALG::CreateVector(*boundarydis_->DofColMap(),true);

  // map to fluid parallel distribution
  LINALG::Export(*idispnp_,*idispcolnp);
  LINALG::Export(*idispn_ ,*idispcoln);

  boundarydis_->SetState("idispcolnp",idispcolnp);
  boundarydis_->SetState("idispcoln" ,idispcoln);
  
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
void ADAPTER::XFluidImpl::LiftDrag()
{
  // get forces on all procs
  // create interface DOF vectors using the fluid parallel distribution
  Teuchos::RCP<Epetra_Vector> iforcecol = LINALG::CreateVector(*boundarydis_->DofColMap(),true);

  // map to fluid parallel distribution
  LINALG::Export(*itrueresnp_,*iforcecol);
  
  
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
    std::stringstream s;
    std::stringstream header;

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

    std::cout << header.str() << endl << s.str() << endl;
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
  set<int> internal_surfeles;
  for (int isurf=0; isurf<boundarydis_->NumMyRowElements(); ++isurf)
  {
    const DRT::Element* surfele = boundarydis_->lRowElement(isurf);
    
    LINALG::SerialDenseMatrix          xyze_surf(3,surfele->NumNode());
    GEO::fillInitialPositionArray(surfele,xyze_surf);
   
    // center in local coordinates
    const LINALG::Matrix<2,1> localcenterpos(DRT::UTILS::getLocalCenterPosition<2>(surfele->Shape()));
    // center in physical coordinates
    static LINALG::Matrix<3,1> physicalcenterpos;
    GEO::elementToCurrentCoordinates(surfele->Shape(), xyze_surf, localcenterpos, physicalcenterpos);
    
    LINALG::Matrix<3,1> unitnormalvec;
    GEO::computeNormalToSurfaceElement(surfele, xyze_surf, localcenterpos, unitnormalvec);
    
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
        const double length = pow(GEO::ElementVolume(*solidele),1.0/3.0);
        
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
  for (set<int>::const_iterator ele=internal_surfeles.begin();
       ele!=internal_surfeles.end();
       ++ele)
  {
    boundarydis_->DeleteElement(*ele);
  }
  const int err2 = boundarydis_->FillComplete();
  if (err2) dserror("FillComplete() returned err=%d",err2);
}
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
  // Delta d(n+1,i+1) = ( theta Delta u(n+1,i+1) + u(n) ) * dt
  //
  double timescale = TimeScaling();
  fcx->Update(-timescale*fluid_.Dt(),*veln,timescale);
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
  dserror("not implemented!");
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::SetTimeLomaFields(RCP<const Epetra_Vector> densnp,RCP<const Epetra_Vector> densn,RCP<const Epetra_Vector> densnm,RCP<const Epetra_Vector> scatraresidual,const double eosfac)
{
   dserror("not implemented!");
   return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::SetIterLomaFields(RCP<const Epetra_Vector> densnp,RCP<const Epetra_Vector> densdtnp,const double eosfac)
{
   dserror("not implemented!");
   return;
}

#endif  // #ifdef CCADISCRET
