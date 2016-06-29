/*----------------------------------------------------------------------*/
/*!
\file fsi_xfem_algorithm.cpp

\brief Basis of monolithic XFSI algorithm that performs a coupling between the
       structural field equation and XFEM fluid field equations

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/


#include "../drt_adapter/ad_fld_base_algorithm.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_fluid_xfluid/xfluid.H"
#include "../drt_adapter/ad_ale.H"
#include "../drt_adapter/ad_ale_fpsi.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_discret.H"

#include "fsi_xfem_algorithm.H"



//! Note: The order of calling the two BaseAlgorithm-constructors is
//! important here! In here control file entries are written. And these entries
//! define the order in which the filters handle the Discretizations, which in
//! turn defines the dof number ordering of the Discretizations.
/*----------------------------------------------------------------------*
 | constructor (public)                                    schott 08/14 |
 *----------------------------------------------------------------------*/
FSI::AlgorithmXFEM::AlgorithmXFEM(const Epetra_Comm& comm,
                                  const Teuchos::ParameterList& timeparams)
  : AlgorithmBase(comm, timeparams),
    structp_block_(0),
    fluid_block_(1),
    ale_i_block_(2)
{
  // access structural dynamic params list which will be possibly modified while creating the time integrator
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();
  const Teuchos::ParameterList& xfdyn = DRT::Problem::Instance()->XFluidDynamicParams();
  bool ale = DRT::INPUT::IntegralValue<bool>((xfdyn.sublist("GENERAL")),"ALE_XFluid");

  //--------------------------------------------
  // ask base algorithm for the structural time integrator
  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");
  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure = Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(timeparams, const_cast<Teuchos::ParameterList&>(sdyn), structdis));
  structure_ = Teuchos::rcp_dynamic_cast<ADAPTER::FSIStructureWrapper>(structure->StructureField());

  if(structure_ == Teuchos::null)
    dserror("cast from ADAPTER::Structure to ADAPTER::FSIStructureWrapper failed");


  //--------------------------------------------
  // ask base algorithm for the fluid time integrator
  //do not init in ale case!!! (will be done in MonolithicAFSI_XFEM::Setup System())
  Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluid = Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(timeparams,fdyn,"fluid",ale,!ale));
  fluid_ = Teuchos::rcp_dynamic_cast<FLD::XFluid>(fluid->FluidField());
  if (fluid_ == Teuchos::null)
    dserror("Cast of Fluid to XFluid failed! - Everything fine in SetupFluid()?");

  if (ale)
  {
    // ask base algorithm for the ale time integrator
    Teuchos::RCP<DRT::Discretization> aledis = DRT::Problem::Instance()->GetDis("ale");
     Teuchos::RCP<ADAPTER::AleBaseAlgorithm> ale = Teuchos::rcp(new ADAPTER::AleBaseAlgorithm(timeparams, aledis));
     ale_ =  Teuchos::rcp_dynamic_cast<ADAPTER::AleFpsiWrapper>(ale->AleField());
     if(ale_ == Teuchos::null)
        dserror("cast from ADAPTER::Ale to ADAPTER::AleFpsiWrapper failed");

    // build a proxy of the fluid discretization for the structure field
    Teuchos::RCP<DRT::DofSet> aledofset = ale_->WriteAccessDiscretization()->GetDofSetProxy();
    if (fluid_->Discretization()->AddDofSet(aledofset) != 1)
      dserror("Fluid Discretization does not have two Dofsets (Fluid/Ale)!");
  }
  else
    ale_ = Teuchos::null;


  return;
}



/*----------------------------------------------------------------------*
 | destructor (public)                                     schott 08/14 |
 *----------------------------------------------------------------------*/
FSI::AlgorithmXFEM::~AlgorithmXFEM()
{
}



/*----------------------------------------------------------------------*
 | update (protected)                                      schott 08/14 |
 *----------------------------------------------------------------------*/
void FSI::AlgorithmXFEM::Update()
{
  dserror("currently unused");

  StructureField()->Update();
  FluidField()->Update();
  if (HaveAle()) AleField()->Update();

  return;
}



/*----------------------------------------------------------------------*
 | calculate stresses, strains, energies                   schott 08/14 |
 *----------------------------------------------------------------------*/
void FSI::AlgorithmXFEM::PrepareOutput()
{
  StructureField()->PrepareOutput();
}
