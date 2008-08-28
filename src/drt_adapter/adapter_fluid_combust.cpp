/*------------------------------------------------------------------------------------------------*/
/*!
\file adapter_fluid_combust.cpp

\brief Fluid field adapter

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*/
/*------------------------------------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "adapter_fluid_combust.H"

#include "../drt_lib/drt_condition_utils.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_lib/drt_globalproblem.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>

extern "C" /* stuff which is c and is accessed from c++ */
{
#include "../headers/standardtypes.h"
}
extern struct _FILES  allfiles;
extern struct _GENPROB     genprob;
/*------------------------------------------------------------------------------------------------*
 | constructor                                                                        henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
ADAPTER::FluidCombust::FluidCombust(Teuchos::RCP<DRT::Discretization> dis,
                                    Teuchos::RCP<LINALG::Solver> solver,
                                    Teuchos::RCP<ParameterList> params,
                                    Teuchos::RCP<IO::DiscretizationWriter> output)
  : fluid_(dis, *solver, *params, *output),
    dis_(dis),
    solver_(solver),
    params_(params),
    output_(output)
{
  vector<string> conditions_to_copy;
  conditions_to_copy.push_back("FSICoupling");
  conditions_to_copy.push_back("XFEMCoupling");
//  boundarydis_ = DRT::UTILS::CreateDiscretizationFromCondition(soliddis, "FSICoupling", "Boundary", "BELE3", conditions_to_copy);
  boundarydis_ = DRT::UTILS::CreateDiscretizationFromCondition(null, "FSICoupling", "Boundary", "BELE3", conditions_to_copy);
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
//  DRT::UTILS::SetupNDimExtractor(*boundarydis_,"FREESURFCoupling",freesurface_);

  // create interface DOF vectors using the solid parallel distribution
  const Epetra_Map* fluidsurface_dofrowmap = boundarydis_->DofRowMap();
//  idispnp_    = LINALG::CreateVector(*fluidsurface_dofrowmap,true);
  ivelnp_     = LINALG::CreateVector(*fluidsurface_dofrowmap,true);
  itrueresnp_ = LINALG::CreateVector(*fluidsurface_dofrowmap,true);

//  idispn_   = LINALG::CreateVector(*fluidsurface_dofrowmap,true);
  iveln_    = LINALG::CreateVector(*fluidsurface_dofrowmap,true);
  ivelnm_   = LINALG::CreateVector(*fluidsurface_dofrowmap,true);
  iaccn_    = LINALG::CreateVector(*fluidsurface_dofrowmap,true);

//  fluid_.SetFreeSurface(&freesurface_);
  std::cout << "FluidCombust constructor done" << endl;
}
/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ADAPTER::FluidCombust::Discretization()
{
  return boundarydis_;
}
/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void ADAPTER::FluidCombust::PrepareTimeStep()
{
  fluid_.PrepareTimeStep();
}
/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void ADAPTER::FluidCombust::Evaluate(Teuchos::RCP<const Epetra_Vector> vel)
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
/*------------------------------------------------------------------------------------------------*
 | henke 08/08  |
 *------------------------------------------------------------------------------------------------*/
void ADAPTER::FluidCombust::Update()
{
  fluid_.TimeUpdate();

//  henke 08/08 const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
  const double dt = params_->get<double>("TIMESTEP");

  // compute acceleration at timestep n
  Teuchos::RCP<Epetra_Vector> iaccnp = rcp(new Epetra_Vector(iaccn_->Map()));
//  Teuchos::RCP<Epetra_Vector> ivelnp = rcp(new Epetra_Vector(iveln_->Map()));
  const double theta = 1.0;
  iaccnp->Update(-(1.0-theta)/(theta),*iaccn_,0.0);
  iaccnp->Update(1.0/(theta*dt),*ivelnp_,-1.0/(theta*dt),*iveln_,1.0);

//  const double beta = 1.0/4.0;
//  const double gamma = 1.0/2.0;
//  
//  iaccnp->Update(-(1.0-(2.0*beta))/(2.0*beta),*iaccn_,0.0);
//  iaccnp->Update(-1.0/(beta*dt),*iveln_,1.0);
//  iaccnp->Update(1.0/(beta*dt*dt),*idispnp_,-1.0/(beta*dt*dt),*idispn_,1.0);
//  
//  ivelnp->Update(1.0,*iveln_,0.0);
//  ivelnp->Update(gamma*dt,*iaccnp,(1-gamma)*dt,*iaccn_,1.0);
  
  // update acceleration n
  iaccn_->Update(1.0,*iaccnp,0.0);

  // update velocity n-1
  ivelnm_->Update(1.0,*iveln_,0.0);

  // update velocity n
  iveln_->Update(1.0,*ivelnp_,0.0);
//  iveln_->Update(1.0,*ivelnp,0.0);
  
  // update displacement n
//  idispn_->Update(1.0,*idispnp_,0.0);
}
/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void ADAPTER::FluidCombust::Output()
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
  Teuchos::RCP<Epetra_Vector> ivelcol     = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  Teuchos::RCP<Epetra_Vector> ivelncol    = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  Teuchos::RCP<Epetra_Vector> ivelnmcol   = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  Teuchos::RCP<Epetra_Vector> iaccncol    = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  Teuchos::RCP<Epetra_Vector> idispcol    = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
  Teuchos::RCP<Epetra_Vector> itruerescol = LINALG::CreateVector(*fluidsurface_dofcolmap,true);

  // map to fluid parallel distribution
//  LINALG::Export(*idispnp_ ,*idispcol);
  LINALG::Export(*ivelnp_  ,*ivelcol);
  LINALG::Export(*iveln_   ,*ivelncol);
  LINALG::Export(*ivelnm_  ,*ivelnmcol);
  LINALG::Export(*iaccn_   ,*iaccncol);

  LINALG::Export(*itrueresnp_,*itruerescol);

  PrintInterfaceVectorField(idispcol, itruerescol, "_solution_iforce_", "interface force");
  PrintInterfaceVectorField(idispcol, ivelcol  , "_solution_ivel_"  , "interface velocity n+1");
  PrintInterfaceVectorField(idispcol, ivelncol , "_solution_iveln_" , "interface velocity n");
  PrintInterfaceVectorField(idispcol, ivelnmcol, "_solution_ivelnm_", "interface velocity n-1");
  PrintInterfaceVectorField(idispcol, iaccncol , "_solution_iaccn_" , "interface acceleration n");

}
/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void ADAPTER::FluidCombust::NonlinearSolve()
{

  cout << "FluidCombust::NonlinearSolve()" << endl;

  // create interface DOF vectors using the fluid parallel distribution
//  const Epetra_Map* fluidsurface_dofcolmap = boundarydis_->DofColMap();
//  Teuchos::RCP<Epetra_Vector> ivelcol     = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
//  Teuchos::RCP<Epetra_Vector> idispcol    = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
//  Teuchos::RCP<Epetra_Vector> itruerescol = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
//
//  Teuchos::RCP<Epetra_Vector> ivelncol = LINALG::CreateVector(*fluidsurface_dofcolmap,true);
//  Teuchos::RCP<Epetra_Vector> iaccncol = LINALG::CreateVector(*fluidsurface_dofcolmap,true);

  // map to fluid parallel distribution
//  LINALG::Export(*ivel_,*ivelcol);
//  LINALG::Export(*idisp_,*idispcol);
//
//  LINALG::Export(*iveln_,*ivelncol);
//  LINALG::Export(*iaccn_,*iaccncol);

  cout << "SetState" << endl;
//  boundarydis_->SetState("idispcolnp",idispnp_);
  boundarydis_->SetState("ivelcolnp",ivelnp_);
  
//  boundarydis_->SetState("idispcoln",idispn_);
  boundarydis_->SetState("ivelcoln",iveln_);
  boundarydis_->SetState("iacccoln",iaccn_);
  
  fluid_.NonlinearSolve(boundarydis_);
  
  Teuchos::RCP<const Epetra_Vector> itruerescol = boundarydis_->GetState("iforcenp");
  
  boundarydis_->ClearState();
  

  // map back to solid parallel distribution
  LINALG::Export(*itruerescol,*itrueresnp_);
}
/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidCombust::InnerVelocityRowMap()
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
/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidCombust::VelocityRowMap()
{
  return fluid_.VelocityRowMap();
}
/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidCombust::PressureRowMap()
{
  return fluid_.PressureRowMap();
}
/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
double ADAPTER::FluidCombust::ResidualScaling() const
{
  return fluid_.ResidualScaling();
}
/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
double ADAPTER::FluidCombust::TimeScaling() const
{
  return 1./fluid_.Dt();
}
/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void ADAPTER::FluidCombust::ReadRestart(int step)
{
  fluid_.ReadRestart(step);
}
/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
double ADAPTER::FluidCombust::Time()
{
  return fluid_.Time();
}
/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
int ADAPTER::FluidCombust::Step()
{
  return fluid_.Step();
}
/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
int ADAPTER::FluidCombust::Itemax() const
{
  return fluid_.Itemax();
}
/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void ADAPTER::FluidCombust::SetItemax(int itemax)
{
  fluid_.SetItemax(itemax);
}
/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void ADAPTER::FluidCombust::LiftDrag()
{
  fluid_.LiftDrag();
}
/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidCombust::ExtractInterfaceForces()
{
  return interface_.ExtractCondVector(itrueresnp_);
}
/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidCombust::ExtractInterfaceVeln()
{
  return interface_.ExtractCondVector(iveln_);
}
/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void ADAPTER::FluidCombust::ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
{
  interface_.InsertCondVector(ivel,ivelnp_);
}
/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::FluidCombust::CreateFieldTest()
{
  return Teuchos::rcp(new FLD::CombustFluidResultTest(fluid_));
}
/*------------------------------------------------------------------------------------------------*
 | Es sieht so aus als hätte Axel diese Funktion mittlerweile wieder eingeführt       henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
//Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidCombust::ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres)
//{
//  dserror("not implemented!");
//  return (fluid_.VelPresSplitter()).ExtractOtherVector(velpres);
//}
/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void ADAPTER::FluidCombust::PrintInterfaceVectorField(
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

#endif  // #ifdef CCADISCRET
