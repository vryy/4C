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

#include "../drt_lib/drt_condition_utils.H"
#include "../drt_io/io_gmsh.H"

extern "C" /* stuff which is c and is accessed from c++ */
{
#include "../headers/standardtypes.h"
}
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::XFluidImpl::XFluidImpl(
        RCP<DRT::Discretization> dis,
        const RCP<DRT::Discretization> soliddis,
        RCP<LINALG::Solver> solver,
        RCP<ParameterList> params,
        RCP<IO::DiscretizationWriter> output,
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

  // create node and element distribution with elements and nodes ghosted on all processors
  const Epetra_Map noderowmap = *boundarydis_->NodeRowMap();
  std::cout << "noderowmap->UniqueGIDs(): " << noderowmap.UniqueGIDs() << endl;
  //std::cout << noderowmap << endl;

  Teuchos::RCP<Epetra_Map> newnodecolmap = LINALG::AllreduceEMap(noderowmap);
  //std::cout << *newnodecolmap << endl;

  DRT::UTILS::RedistributeWithNewNodalDistribution(*boundarydis_, noderowmap, *newnodecolmap);



  UTILS::SetupNDimExtractor(*boundarydis_,"FSICoupling",interface_);
  UTILS::SetupNDimExtractor(*boundarydis_,"FREESURFCoupling",freesurface_);

  // create interface DOF vectors using the solid parallel distribution
  const Epetra_Map* fluidsurface_dofrowmap = boundarydis_->DofRowMap();
  ivel_     = LINALG::CreateVector(*fluidsurface_dofrowmap,true);
  idisp_    = LINALG::CreateVector(*fluidsurface_dofrowmap,true);
  itrueres_ = LINALG::CreateVector(*fluidsurface_dofrowmap,true);

  // create interface DOF vectors using the fluid parallel distribution
  const Epetra_Map* fluidsurface_dofcolmap = boundarydis_->DofColMap();
  ivelcol_     = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  idispcol_    = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  itruerescol_ = LINALG::CreateVector(*fluidsurface_dofcolmap,true);

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
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::Output()
{
  fluid_.Output();
  
  static int step_counter = 1;
  
  
  std::stringstream filename;
  filename << allfiles.outputfile_kenner << "_solution_interfaceforce_" << std::setw(5) << setfill('0') << step_counter << ".pos";
  std::cout << "writing '"<<filename.str()<<"'...";
  std::ofstream f_system(filename.str().c_str());

  {
    stringstream gmshfilecontent;
    gmshfilecontent << "View \" " << "Interface Force \" {" << endl;
    for (int i=0; i<boundarydis_->NumMyColElements(); ++i)
    {
      const DRT::Element* actele = boundarydis_->lColElement(i);
//      cout << *actele << endl;
      vector<int> lm;
      vector<int> lmowner;
      actele->LocationVector(*boundarydis_, lm, lmowner);

      // extract local values from the global vector
      vector<double> myvelnp(lm.size());
      DRT::UTILS::ExtractMyValues(*itruerescol_, myvelnp, lm);
      
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*idispcol_, mydisp, lm);

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
  
  step_counter++;
  
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::XFluidImpl::NonlinearSolve()
{

  cout << "XFluidImpl::NonlinearSolve()" << endl;
  fluid_.NonlinearSolve(boundarydis_,idispcol_,ivelcol_,itruerescol_);
  LINALG::Export(*itruerescol_,*itrueres_);
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
  return interface_.ExtractCondVector(itrueres_);
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
void ADAPTER::XFluidImpl::ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
{
//  cout << "applying interface velocity" << endl;

  interface_.InsertCondVector(ivel,ivel_);
  LINALG::Export(*ivel_,*ivelcol_);
//
//  for (int entry = 0; entry < ivelcol_->MyLength(); ++entry)
//  {
//    if (entry%3 == 0)
//    {
//      //(*ivelcol)[entry] = (-1.5*std::sin(2.0*time_* PI) * PI);
//      (*ivelcol_)[entry] = -5.0;
//    }
//    else
//    {
//      (*ivelcol_)[entry] = 0.0;
//    }
//
//  }

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
  interface_.InsertCondVector(idisp,idisp_);
  LINALG::Export(*idisp_,*idispcol_);


//  //idispcol_->PutScalar( 0.0); // ganz schlecht
//  //idispcol_->PutScalar( 0.049); // naja
//  //idispcol_->PutScalar(-0.149); // works best
//  for (int entry = 0; entry < idispcol_->MyLength(); ++entry)
//  {
//    if (entry%3 == 0)
//    {
//      //(*ivelcol)[entry] = (-1.5*std::sin(2.0*time_* PI) * PI);
//      //(*idispcol_)[entry] = -0.149;
//      (*idispcol_)[entry] = 0.049;
//    }
//    else
//    {
//      (*idispcol_)[entry] = 0.0;
//    }
//
//  }

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
  fluid_.LinearRelaxationSolve(relax);
  return ExtractInterfaceForces();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::XFluidImpl::CreateFieldTest()
{
  return Teuchos::rcp(new XFluidResultTest(fluid_));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::XFluidImpl::ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres)
{
    dserror("not implemented!");
   return (fluid_.VelPresSplitter()).ExtractOtherVector(velpres);
}


#endif  // #ifdef CCADISCRET
