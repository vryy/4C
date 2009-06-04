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
#include "../drt_lib/linalg_blocksparsematrix.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>

/*------------------------------------------------------------------------------------------------*
 | constructor                                                                        henke 08/08 |
 |
 | It initializes its associated time integration scheme fluid_.
 *------------------------------------------------------------------------------------------------*/
ADAPTER::FluidCombust::FluidCombust(Teuchos::RCP<DRT::Discretization> dis,
                                    Teuchos::RCP<LINALG::Solver> solver,
                                    Teuchos::RCP<ParameterList> params,
                                    Teuchos::RCP<IO::DiscretizationWriter> output)
  : fluid_(dis, *solver, *params, *output), // calls COMBUST::CombustFluidImplicitTimeInt()
    fluiddis_(dis),
    solver_(solver),
    params_(params),
    output_(output)
{
  std::cout << "Proc " << fluiddis_->Comm().MyPID() << "ADAPTER::FluidCombust constructor done \n" << endl;
}

void ADAPTER::FluidCombust::SetInitialFlowField(int whichinitialfield, int startfuncno)
{
  return fluid_.SetInitialFlowField(whichinitialfield, startfuncno);
}

Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidCombust::TrueResidual()
{
  return fluid_.TrueResidual();
}

Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidCombust::Veln()
{
  return fluid_.Veln();
}

Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidCombust::Velnp()
{
  return fluid_.Velnp();
}

/*------------------------------------------------------------------------------------------------*
 | return null pointer instead of subgrid velocity/viscosity vector                   henke 05/09 |
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidCombust::SgVelVisc()
{
  /* This function is merely required to be able to call SetVelocityField(), in the constructor of
   * ScaTraFluidCouplingAlgorithm. There, a vector sgvelvisc_ is transferred to the initial scalar
   * transport field by default. The combustion fluid does not have this vector and thus this
   * function returns a null pointer
   */
  return Teuchos::null;
}
/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void ADAPTER::FluidCombust::PrepareTimeStep()
{
  fluid_.PrepareTimeStep();
}

/*------------------------------------------------------------------------------------------------*
 | Wozu ist diese Abfrage nötig? henke 10/08 |
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
 | In der Zeitintegration gibt es keine Update Funktion, das passiert hier!          henke 10/08  |
 *------------------------------------------------------------------------------------------------*/
void ADAPTER::FluidCombust::Update()
{
  fluid_.TimeUpdate();

// henke 08/08 Ich denke die Parameterliste müsste mit params_ eigentlich schon alles beeinhalten!
// -> wozu ist das hier dann gut?
//  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
  const double dt = params_->get<double>("TIMESTEP");

  // compute acceleration at timestep n+1
  Teuchos::RCP<Epetra_Vector> iaccnp = rcp(new Epetra_Vector(iaccn_->Map()));
  Teuchos::RCP<Epetra_Vector> ivelnp = rcp(new Epetra_Vector(iveln_->Map()));
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
  iveln_->Update(1.0,*ivelnp_,0.0);

  // update velocity n
  iveln_->Update(1.0,*ivelnp_,0.0);

}

/*------------------------------------------------------------------------------------------------*
 | Was davon brauche ich? henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void ADAPTER::FluidCombust::Output()
{
  fluid_.Output();
}

/*------------------------------------------------------------------------------------------------*
 | Baue diese Funktion eventuell für dich um!                                         henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void ADAPTER::FluidCombust::PrintInterfaceVectorField(
    const Teuchos::RCP<Epetra_Vector>   displacementfield,
    const Teuchos::RCP<Epetra_Vector>   vectorfield,
    const std::string filestr,
    const std::string name_in_gmsh
    )
{
/*  const Teuchos::ParameterList& xfemparams = DRT::Problem::Instance()->XFEMGeneralParams();
  const bool gmshdebugout = (bool)getIntegralValue<int>(xfemparams,"GMSH_DEBUG_OUT");
  if (gmshdebugout)
  {
    std::stringstream filename;
    std::stringstream filenamedel;
    const std::string filebase = DRT::Problem::Instance()->OutputControlFile()->FileName();
    filename << filebase << filestr << std::setw(5) << setfill('0') << Step() << ".pos";
    filenamedel << filebase << filestr << std::setw(5) << setfill('0') << Step()-5 << ".pos";
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
        LINALG::SerialDenseMatrix elementvalues(nsd,numnode);
        LINALG::SerialDenseMatrix elementpositions(nsd,numnode);
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
*/
}

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void ADAPTER::FluidCombust::NonlinearSolve()
{

  std::cout << "Proc " << fluiddis_->Comm().MyPID()<< "FluidCombust::NonlinearSolve()" << endl;

  fluid_.NonlinearSolve();
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
double ADAPTER::FluidCombust::Time() const
{
  return fluid_.Time();
}

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
int ADAPTER::FluidCombust::Step() const
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
  dserror("Thou shalt not call this function!");
  //fluid_.SetItemax(itemax);
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
  dserror("Don't use this function, I don't know what it does!");
  return interface_.ExtractCondVector(itrueresnp_);
}

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidCombust::ExtractInterfaceVeln()
{
  dserror("Don't use this function, I don't know what it does!");
  return interface_.ExtractCondVector(iveln_);
}

/*------------------------------------------------------------------------------------------------*
 | henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
void ADAPTER::FluidCombust::ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
{
  dserror("Don't use this function, I don't know what it does!");
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
 | Zum Exportieren der Konvektionsgeschw. Warum die Bezeichnung velpres? Druck?       henke 08/08 |
 *------------------------------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidCombust::ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres)
{
  return (fluid_.VelPresSplitter()).ExtractOtherVector(velpres);
}

void ADAPTER::FluidCombust::ImportDiscretization(Teuchos::RCP<DRT::Discretization> importdis)
{
  /* This function is not needed, because the ImportInterface is used instead. Besides, it does not work
   * anymore, since the InterfaceHandle receives a FlameFront and the discretizations. henke 01/09 */
  dserror("Thou shalt not call this function!");

  // import level set discretization from combustion algorithm
  gfuncdis_ = importdis;

  // construct a combustion interface handle (couldn't we also construct a flame front here!) and
  // thus automatically process the flame front
//  interfacehandle = rcp(new COMBUST::InterfaceHandleCombust(fluiddis_,gfuncdis_));

  // pass geometrical information aboout flame front to fluid time integration scheme
//  fluid_.IncorporateInterface(interfacehandle);
}

void ADAPTER::FluidCombust::ImportInterface(const Teuchos::RCP<COMBUST::InterfaceHandleCombust>& interfacehandle)
{
  // pass geometrical information aboout flame front to fluid time integration scheme
  fluid_.IncorporateInterface(interfacehandle);
}

#endif  // #ifdef CCADISCRET
