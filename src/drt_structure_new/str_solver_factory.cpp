/*-----------------------------------------------------------*/
/*!
\file str_solver_factory.cpp

\brief str_solver_factory

\maintainer Michael Hiermeier

\date Sep 9, 2015

\level 3

*/
/*-----------------------------------------------------------*/


#include "str_solver_factory.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"

#include "../drt_io/io_control.H"

#include "../drt_inpar/inpar_structure.H"
#include "../drt_inpar/inpar_contact.H"

#include <Teuchos_ParameterList.hpp>


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::SOLVER::Factory::Factory()
{
  // empty
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::SOLVER::Factory::LinSolMap>
STR::SOLVER::Factory::BuildLinSolvers(
    const std::set<enum INPAR::STR::ModelType>& modeltypes,
    const Teuchos::ParameterList& sdyn,
    DRT::Discretization& actdis) const
{
  // create a new standard map
  Teuchos::RCP<LinSolMap> linsolvers = Teuchos::rcp(new LinSolMap());

  std::set<enum INPAR::STR::ModelType>::const_iterator mt_iter;
  // loop over all model types
  for (mt_iter = modeltypes.begin(); mt_iter != modeltypes.end(); ++mt_iter)
  {
    switch(*mt_iter)
    {
    case INPAR::STR::model_structure:
    case INPAR::STR::model_springdashpot:
    {
      /* Check if the structural linear solver was already added and skip
       * if true. */
      LinSolMap::iterator iter =
          linsolvers->find(INPAR::STR::model_structure);
      if (iter==linsolvers->end())
        (*linsolvers)[INPAR::STR::model_structure] =
            BuildStructureLinSolver(sdyn,actdis);
      break;
    }
    /* ToDo Check if this makes sense for simulations where both, meshtying and
     *      contact, are present. If we need two linsolvers, please adjust the
     *      implementation (maps for pre-conditioning, etc.). */
    case INPAR::STR::model_contact:
    case INPAR::STR::model_meshtying:
      (*linsolvers)[*mt_iter] =
          BuildMeshtyingContactLinSolver(sdyn,actdis);
      break;
    case INPAR::STR::model_lag_pen_constraint:
      (*linsolvers)[*mt_iter] =
          BuildLagPenConstraintLinSolver(sdyn,actdis);
      break;
    case INPAR::STR::model_cardiovascular0d:
      (*linsolvers)[*mt_iter] = BuildWindkesselLinSolver(sdyn,actdis);
      break;
    default:
      dserror("No idea which solver to use for the given model type %s",
          ModelTypeString(*mt_iter).c_str());
      break;
    }
  }

  return linsolvers;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::Solver> STR::SOLVER::Factory::BuildStructureLinSolver(
        const Teuchos::ParameterList& sdyn,
        DRT::Discretization& actdis) const
{
  Teuchos::RCP<LINALG::Solver> linsolver = Teuchos::null;

  // get the linear solver number used for structural problems
  const int linsolvernumber = sdyn.get<int>("LINEAR_SOLVER");
  // check if the structural solver has a valid solver number
  if (linsolvernumber == (-1))
    dserror("no linear solver defined for structural field. "
        "Please set LINEAR_SOLVER in STRUCTURAL DYNAMIC to a valid number!");

  linsolver = Teuchos::rcp(
      new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
      actdis.Comm(),
      DRT::Problem::Instance()->ErrorFile()->Handle()));

  actdis.ComputeNullSpaceIfNecessary(linsolver->Params());

  if((linsolver->Params().isSublist("Aztec Parameters") or
      linsolver->Params().isSublist("Belos Parameters")) and
      linsolver->Params().isSublist("ML Parameters") and                  // TODO what about MueLu?
      DRT::INPUT::IntegralValue<INPAR::STR::STC_Scale>(sdyn,"STC_SCALING")!=INPAR::STR::stc_none)
  {
    Teuchos::ParameterList& mllist = linsolver->Params().sublist("ML Parameters");
    Teuchos::RCP<std::vector<double> > ns = mllist.get<Teuchos::RCP<std::vector<double> > >("nullspace");

    const int size=actdis.DofRowMap()->NumMyElements();

    // extract the six nullspace vectors corresponding to the modes
    // trans x, trans y, trans z, rot x, rot y, rot z
    // Note: We assume 3d here!
    Teuchos::RCP<Epetra_Vector> nsv1=
        Teuchos::rcp(new Epetra_Vector(View,*(actdis.DofRowMap()),&((*ns)[0])));
    Teuchos::RCP<Epetra_Vector> nsv2=
        Teuchos::rcp(new Epetra_Vector(View,*(actdis.DofRowMap()),&((*ns)[size])));
    Teuchos::RCP<Epetra_Vector> nsv3=
        Teuchos::rcp(new Epetra_Vector(View,*(actdis.DofRowMap()),&((*ns)[2*size])));
    Teuchos::RCP<Epetra_Vector> nsv4=
        Teuchos::rcp(new Epetra_Vector(View,*(actdis.DofRowMap()),&((*ns)[3*size])));
    Teuchos::RCP<Epetra_Vector> nsv5=
        Teuchos::rcp(new Epetra_Vector(View,*(actdis.DofRowMap()),&((*ns)[4*size])));
    Teuchos::RCP<Epetra_Vector> nsv6=
        Teuchos::rcp(new Epetra_Vector(View,*(actdis.DofRowMap()),&((*ns)[5*size])));


    //prepare matrix for scaled thickness business of thin shell structures
    Teuchos::RCP<LINALG::SparseMatrix> stcinv=
      Teuchos::rcp(new LINALG::SparseMatrix(*actdis.DofRowMap(), 81, true, true));

    stcinv->Zero();
    // create the parameters for the discretization
    Teuchos::ParameterList p;
    // action for elements
    const std::string action = "calc_stc_matrix_inverse";
    p.set("action", action);
    p.set<int>("stc_scaling",DRT::INPUT::IntegralValue<INPAR::STR::STC_Scale>(sdyn,"STC_SCALING"));
    p.set("stc_layer",1);

    actdis.Evaluate(p, stcinv, Teuchos::null,  Teuchos::null, Teuchos::null, Teuchos::null);

    stcinv->Complete();

    for (int lay = 2; lay <= sdyn.get<int>("STC_LAYER"); ++lay)
    {
      Teuchos::ParameterList pe;

      p.set("stc_layer", lay);

      Teuchos::RCP<LINALG::SparseMatrix> tmpstcmat=
        Teuchos::rcp(new LINALG::SparseMatrix(*(actdis.DofRowMap()),81,true,true));
      tmpstcmat->Zero();

      actdis.Evaluate(p, tmpstcmat, Teuchos::null,  Teuchos::null, Teuchos::null, Teuchos::null);
      tmpstcmat->Complete();

      stcinv = MLMultiply(*stcinv,*tmpstcmat,false,false,true);
    }

    Teuchos::RCP<Epetra_Vector> temp =
        LINALG::CreateVector(*(actdis.DofRowMap()),false);

    stcinv->Multiply(false,*nsv1,*temp);
    nsv1->Update(1.0,*temp,0.0);
    stcinv->Multiply(false,*nsv2,*temp);
    nsv2->Update(1.0,*temp,0.0);
    stcinv->Multiply(false,*nsv3,*temp);
    nsv3->Update(1.0,*temp,0.0);
    stcinv->Multiply(false,*nsv4,*temp);
    nsv4->Update(1.0,*temp,0.0);
    stcinv->Multiply(false,*nsv5,*temp);
    nsv5->Update(1.0,*temp,0.0);
    stcinv->Multiply(false,*nsv6,*temp);
    nsv6->Update(1.0,*temp,0.0);
  }

  return linsolver;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::Solver> STR::SOLVER::Factory::BuildMeshtyingContactLinSolver(
        const Teuchos::ParameterList& sdyn,
        DRT::Discretization& actdis) const
{
  Teuchos::RCP<LINALG::Solver> linsolver = Teuchos::null;

  // get mortar information
  std::vector<DRT::Condition*> mtcond(0);
  std::vector<DRT::Condition*> ccond(0);
  actdis.GetCondition("Mortar", mtcond);
  actdis.GetCondition("Contact",ccond);
  bool onlymeshtying       = false;
  bool onlycontact         = false;
  bool meshtyingandcontact = false;
  if(mtcond.size()!=0 and ccond.size()!=0)
    meshtyingandcontact = true;
  if(mtcond.size()!=0 and ccond.size()==0)
    onlymeshtying = true;
  if(mtcond.size()==0 and ccond.size()!=0)
    onlycontact = true;

  const Teuchos::ParameterList& mcparams     = DRT::Problem::Instance()->ContactDynamicParams();
  switch(DRT::INPUT::IntegralValue<int>(mcparams,"SYSTEM"))
  {
    case INPAR::CONTACT::system_saddlepoint:
    {
      // meshtying/contact for structure
      // get the solver number used for meshtying/contact problems
      const int linsolvernumber = mcparams.get<int>("LINEAR_SOLVER");
      // check if the meshtying/contact solver has a valid solver number
      if (linsolvernumber == (-1))
        dserror("no linear solver defined for meshtying/contact problem. Please set LINEAR_SOLVER in CONTACT DYNAMIC to a valid number!");

      // plausibility check

      // solver can be either UMFPACK (direct solver) or an Aztec_MSR/Belos (iterative solver)
      INPAR::SOLVER::SolverType sol  = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(DRT::Problem::Instance()->SolverParams(linsolvernumber),"SOLVER");
      INPAR::SOLVER::AzPrecType prec = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(DRT::Problem::Instance()->SolverParams(linsolvernumber),"AZPREC");
      if (sol != INPAR::SOLVER::umfpack && sol != INPAR::SOLVER::superlu) {
        // if an iterative solver is chosen we need a block preconditioner like CheapSIMPLE
        if (prec != INPAR::SOLVER::azprec_CheapSIMPLE &&
            prec != INPAR::SOLVER::azprec_TekoSIMPLE  &&
            prec != INPAR::SOLVER::azprec_MueLuAMG_contactSP)  // TODO adapt error message
          dserror("You have chosen an iterative linear solver. For mortar/Contact in saddlepoint formulation you have to choose a block preconditioner such as SIMPLE. Choose CheapSIMPLE, TekoSIMPLE (if Teko is available) or MueLu_contactSP (if MueLu is available) in the SOLVER %i block in your dat file.",linsolvernumber);
      }

      // build meshtying/contact solver
      linsolver =
      Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
                             actdis.Comm(),
                             DRT::Problem::Instance()->ErrorFile()->Handle()));

      actdis.ComputeNullSpaceIfNecessary(linsolver->Params());

      // feed the solver object with additional information
      if     (onlycontact or meshtyingandcontact) linsolver->Params().set<bool>("CONTACT",true);
      else if(onlymeshtying)                      linsolver->Params().set<bool>("MESHTYING",true);
      else dserror("this cannot be: no saddlepoint problem for beamcontact or pure structure problem.");

      INPAR::CONTACT::SolvingStrategy soltype = DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(mcparams,"STRATEGY");
      if (soltype==INPAR::CONTACT::solution_lagmult)
      {
        // get the solver number used for structural problems
        const int linsolvernumber = sdyn.get<int>("LINEAR_SOLVER");
        // check if the structural solver has a valid solver number
        if (linsolvernumber == (-1))
          dserror("no linear solver defined for structural field. Please set LINEAR_SOLVER in STRUCTURAL DYNAMIC to a valid number!");

        // provide null space information
        if (prec == INPAR::SOLVER::azprec_CheapSIMPLE ||
                  prec == INPAR::SOLVER::azprec_TekoSIMPLE) {
          actdis.ComputeNullSpaceIfNecessary(linsolver->Params().sublist("CheapSIMPLE Parameters").sublist("Inverse1")); // Inverse2 is created within blockpreconditioners.cpp
        } else if (prec == INPAR::SOLVER::azprec_MueLuAMG_contactSP) { /* do nothing here */    }
      }
    }
    break;
    default:
    {
      // meshtying/contact for structure
       // get the solver number used for meshtying/contact problems
       const int linsolvernumber = mcparams.get<int>("LINEAR_SOLVER");
       // check if the meshtying/contact solver has a valid solver number
       if (linsolvernumber == (-1))
         dserror("no linear solver defined for meshtying/contact problem. Please set LINEAR_SOLVER in CONTACT DYNAMIC to a valid number!");

      // build meshtying solver
      linsolver =
      Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
                             actdis.Comm(),
                             DRT::Problem::Instance()->ErrorFile()->Handle()));
      actdis.ComputeNullSpaceIfNecessary(linsolver->Params());
    }
    break;
  }

  return linsolver;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::Solver> STR::SOLVER::Factory::BuildLagPenConstraintLinSolver(
        const Teuchos::ParameterList& sdyn,
        DRT::Discretization& actdis) const
{
  Teuchos::RCP<LINALG::Solver> linsolver = Teuchos::null;

  dserror("Not yet implemented!");

  return linsolver;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::Solver> STR::SOLVER::Factory::BuildWindkesselLinSolver(
        const Teuchos::ParameterList& sdyn,
        DRT::Discretization& actdis) const
{
  Teuchos::RCP<LINALG::Solver> linsolver = Teuchos::null;

  dserror("Not yet implemented!");

  return linsolver;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<std::map<enum INPAR::STR::ModelType, Teuchos::RCP<LINALG::Solver> > >
    STR::SOLVER::BuildLinSolvers(
    const std::set<enum INPAR::STR::ModelType>& modeltypes,
    const Teuchos::ParameterList& sdyn,
    DRT::Discretization& actdis)
{
  Factory factory;
  return factory.BuildLinSolvers(modeltypes,sdyn,actdis);
}
