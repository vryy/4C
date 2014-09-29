/*----------------------------------------------------------------------*/
/*!
\file xfluid_levelset_coupling_algorithm.cpp

\brief Basis of xfluid-levelset coupling.

<pre>
Maintainer: Magnus Winter
            winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089/28915245
</pre>
*/
/*----------------------------------------------------------------------*/

#include "../drt_scatra/scatra_timint_ost.H"

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fluid_xfluid/xfluid.H"

#include "xfluid_levelset_coupling_algorithm.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
XFLUIDLEVELSET::Algorithm::Algorithm(
    const Epetra_Comm& comm,
    const Teuchos::ParameterList& prbdyn,
    const Teuchos::ParameterList& solverparams
    )
:  ScaTraFluidCouplingAlgorithm(comm,prbdyn,false,"scatra",solverparams)
{
  //Note: The ScaTra base algorithm is initialized without a velocity field. This works in this setting, as the level set
  //      field is not propagated. For future use the velocity field should be initialized in this algorithm. However, the
  //      XFluid class will contain enriched nodes. This will have to be dealt with. Inspiration can be found in the combust
  //      algorithm, where this is done.

  //Give Scatra Values to fluid. The cut is performed and the XFluidState-class is initialized.
  //Note: Before this call, the FluidField(), does not have acces to Velnp(). For propagation of the ScatraField(),
  //      a way to provide the ScaTra field with the fluid velocity is necessary for time dependent calculations.
  SetScaTraValuesInFluid();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
XFLUIDLEVELSET::Algorithm::~Algorithm()
{
  return;
}


/*----------------------------------------------------------------------*/
// Stores the phi values on the node map. This is necessary to call the
// level set cut correctly.
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> XFLUIDLEVELSET::Algorithm::StorePhiVectors(Teuchos::RCP<DRT::Discretization> fluiddis, Teuchos::RCP<DRT::Discretization> scatradis, Teuchos::RCP<const Epetra_Vector> phinp)
{
  /* In the processing of the flame front the fluid discretization is cut by the level set
   * function (G-function). This is the reason why fluid elements need to be able to access the
   * corresponding scalar G-function values for all their nodes. Therefore, the ScaTra dof-based
   * vector phinp has to be rearranged in a parallel environment to represent a Fluid node-based
   * vector. This involves two steps:
   * 1. ScaTra DofRowMap -> Fluid NodeRowMap
   * 2. Fluid NodeRowMap -> Fluid NodeColMap
   *
   * henke 02/09
   */

  // reset vectors
  // phin_  = Teuchos::null;
  // phinp_ = Teuchos::null;

  /// Levelset values on node map.
//  Teuchos::RCP<Epetra_Vector> phinpnode=Teuchos::null;

  //phinpnode_=Teuchos::null;


  //------------------------------------------------------------------------------------------------
  // Rearranging phi vectors from ScaTra DofRowMap to fluiddis NodeRowMap
  //------------------------------------------------------------------------------------------------
//  const Teuchos::RCP<Epetra_Vector> phinrow  = Teuchos::rcp(new Epetra_Vector(*fluiddis_->NodeRowMap()));
  const Teuchos::RCP<Epetra_Vector> phinprow = Teuchos::rcp(new Epetra_Vector(*fluiddis->NodeRowMap()));

  //----------------------------------------------------------------------------------------------
  // congruent (matching) discretizations (Fluid == ScaTra)
  //----------------------------------------------------------------------------------------------
  if(true)
  {
    /* The premise of geometrically identical discretizations leads to a significant simplification
     * of the rearrangment process. If nodes of both discretizations are distributed in the same way
     * (which they are) and the dofs belonging to nodes lie on the same processors, then the
     * following holds:
     * G-function NodeMap   identical to   Fluid NodeMap
     * G-function NodeMap     similar to   G-function DofMap   and therefore
     * G-function DofMap      similar to   Fluid NodeMap
     *
     * This means that vector phinp living on the G-function DofRowMap can be directly copied to the
     * Fluid NodeRowMap without any rearrangement for congruent discretizations.
     *
     * Here is what we used to do:
     *
     *  // get the G-function values corresponding to fluid nodes
     *  *phinmrow = *phinm;
     *  *phinrow  = *phin;
     *  *phinprow = *phinp;
     *  // ausgeschrieben, was *phinprow = *phinp in einer Zeile macht:
     *  for(int lnodeid=0;lnodeid<fluiddis_->NumMyRowNodes();lnodeid++)
     *  {
     *    // get the G-function value corresponding to this fluid node
     *    double value = (*phinp)[lnodeid];
     *    phinprow->ReplaceMyValue(lnodeid,value);
     *    }
     *
     * henke 02/09
     *
     * This is not possible anymore if we use periodic boundary conditions. The number of nodes remains
     * the same, but dofs are removed, since master and slave nodes share a dof. Things are not so
     * easy anymore and we have to pick the corresponding G-function dof by looping over all fluid
     * nodes.
     *
     * henke 04/11
     */

    //DRT::UTILS::PrintParallelDistribution(*fluiddis_);
    //DRT::UTILS::PrintParallelDistribution(*gfuncdis_);

    // loop all nodes on the processor
    for(int lfluidnodeid=0;lfluidnodeid<fluiddis->NumMyRowNodes();lfluidnodeid++)
    {
      // get the processor's scatra node
      // remark: we rely on identical parallel distributions of fluid and G-function discretizations, here
      DRT::Node* scatranode = scatradis->lRowNode(lfluidnodeid);
#ifdef DEBUG
      if(scatradis_->NumDof(gfuncnode)!=1) dserror("Levelset node should have 1 dof");
#endif
      // find out the local dof id of the dof at the scatra node
      // this assumes the levelset is set at the dofid = 0! Possible to extend this to more scalars.
      const int scatradofgid = scatradis->Dof(scatranode,0);

//      const int lscatradofidn = phin->Map().LID(gfuncdofgid);
      const int lscatradofidnp = phinp->Map().LID(scatradofgid);

//      if (lscatradofidn < 0) dserror("local dof id not found in map for given global dof id");
      if (lscatradofidnp < 0) dserror("local dof id not found in map for given global dof id");

      // now copy the values
//      const double valuen = (*phin)[lscatradofidn];
//      const int errn = phinrow->ReplaceMyValue(lfluidnodeid,0,valuen);
//      if (errn) dserror("error while inserting value into phinrow");

      const double valuenp = (*phinp)[lscatradofidnp];
      const int errnp = phinprow->ReplaceMyValue(lfluidnodeid,0,valuenp);
      if (errnp) dserror("error while inserting value into phinprow");
    }

  }
  //----------------------------------------------------------------------------------------------
  // no congruent (matching) discretizations (Fluid != G-function)
  //----------------------------------------------------------------------------------------------
  else
  {
    /* In general this operation is very complex. It involves parallel communication and
     * interpolation of the G-function at fluid nodes, if non-matching discretizations for fluid
     * and G-function are used.
     *
     * henke 02/09
     */
    dserror("non-congruent discretizations cannot be handled yet!");
  }

  //------------------------------------------------------------------------------------------------
  // Export vector phinp from fluiddis NodeRowMap to fluiddis NodeColMap for parallel
  // accessibility.
  // remark: SetState() can not be used here, because it is designed for dof-based vectors only.
  //------------------------------------------------------------------------------------------------
//  const Teuchos::RCP<Epetra_Vector> phincol = Teuchos::rcp(new Epetra_Vector(*fluiddis_->NodeColMap()));
//  LINALG::Export(*phinrow,*phincol);
  const Teuchos::RCP<Epetra_Vector> phinpcol = Teuchos::rcp(new Epetra_Vector(*fluiddis->NodeColMap()));
  LINALG::Export(*phinprow,*phinpcol);

  // store vector on fluiddis NodeColMap holding G-function values in member variable
//  phin_  = phincol;
//  phinpnode = phinpcol;
  return phinpcol;
}

/*------------------------------------------------------------------------------------------------*
 | public: algorithm for a stationary TPF problem                                    winter 07/14 |
 *------------------------------------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::SolveStationaryProblem()
{
  if (Comm().MyPID()==0)
  {
    printf("------Stationary-Xfluid-LevelSet-XFEM------  time step ----------------------------------------\n");
  }

  // check time integration schemes of single fields
  // remark: this was already done in ScaTraFluidCouplingAlgorithm() before
  if (FluidField().TimIntScheme() != INPAR::FLUID::timeint_stationary)
    dserror("Fluid time integration scheme is not stationary");
  if (ScaTraField()->MethodName() != INPAR::SCATRA::timeint_stationary)
    dserror("Scatra time integration scheme is not stationary");

  // write Scatra output
  Output();

  // run the simulation, calls the xfluid-"integrate()" routine
  FluidField().Integrate();

//  // solve level set equation
  if (Comm().MyPID()==0)
    std::cout << "/!\\ warning === Level-set field not solved for Fluid_XFEM_LevelSet problems" << std::endl;


  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::SetScaTraValuesInFluid()
{
  // set level set in fluid field

  Teuchos::RCP<Epetra_Vector> phinpnode = StorePhiVectors(FluidField().Discretization(), ScaTraField()->Discretization(), ScaTraField()->Phinp());

  switch(FluidField().TimIntScheme())
  {
  case INPAR::FLUID::timeint_stationary:
  {
    //The state class is initialized in the Fluid().Init() call in FluidBaseAlgorithm().
    //This State is set without any information about the level set. Thus Set StateLS is called,
    //  to set the state correctly for further calculations.
    Teuchos::rcp_dynamic_cast<FLD::XFluid>(FluidFieldrcp())->SetStateLS(*phinpnode);
  }
  break;
  default:
    dserror("Time integration scheme not supported");
    break;
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::Output()
{

  //Xfluid handles output for itself currently. Only ScaTra field output is done here.
  ScaTraField()->Output();

  return;
}


/*----------------------------------------------------------------------*
 | perform result test                                     winter 06/14 |
 *----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::TestResults()
{
  //Currently no testing for ScaTra values, as nothing is changed from input.
  std::cout << "------------TestResults:--Fluid values are compared.\n";
  std::cout << "--------------------------ScaTra values are not compared as of now.\n";

  // perform result tests if required
  DRT::Problem::Instance()->AddFieldTest(FluidField().CreateFieldTest());
  DRT::Problem::Instance()->TestAll(Comm());


//  if (ScaTraField()->MethodName() != INPAR::SCATRA::timeint_gen_alpha)
//  {
//    // DRT::Problem::Instance()->TestAll() is called in level-set field after adding particles
//    Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())->TestResults();
//  }
//  else
//  {
//    DRT::Problem::Instance()->AddFieldTest(CreateScaTraFieldTest());
//    DRT::Problem::Instance()->TestAll(Comm());
//  }

  return;
}
