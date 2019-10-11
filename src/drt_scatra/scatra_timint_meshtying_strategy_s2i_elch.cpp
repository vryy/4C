/*----------------------------------------------------------------------*/
/*! \file

\brief Scatra-scatra interface coupling strategy for electrochemistry problems

\level 2

\maintainer Christoph Schmidt

*----------------------------------------------------------------------*/
#include "scatra_timint_meshtying_strategy_s2i_elch.H"

#include "../drt_adapter/adapter_coupling.H"

#include "../drt_lib/drt_discret.H"

#include "../drt_mat/electrode.H"
#include "../drt_mat/soret.H"

#include "../drt_mortar/mortar_element.H"

#include "../drt_scatra_ele/scatra_ele_boundary_calc_elch_electrode.H"
#include "../drt_scatra_ele/scatra_ele_boundary_calc_elch_electrode_sti_thermo.H"
#include "../drt_scatra_ele/scatra_ele_boundary_calc_sti_electrode.H"
#include "../drt_scatra_ele/scatra_ele_calc_utils.H"
#include "../drt_scatra_ele/scatra_ele_parameter_elch.H"
#include "../drt_scatra_ele/scatra_ele_parameter_timint.H"
#include "../drt_scatra_ele/scatra_ele_parameter_boundary.H"

#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_solver.H"

/*----------------------------------------------------------------------*
 | constructor                                               fang 12/14 |
 *----------------------------------------------------------------------*/
SCATRA::MeshtyingStrategyS2IElch::MeshtyingStrategyS2IElch(
    SCATRA::ScaTraTimIntElch* elchtimint,  //!< elch time integrator
    const Teuchos::ParameterList&
        parameters  //!< input parameters for scatra-scatra interface coupling
    )
    : MeshtyingStrategyS2I(elchtimint, parameters),
      etagrowthmin_(0.),
      intlayergrowth_startstep_(-1),
      intlayergrowth_timestep_active_(false)
{
  return;
}  // SCATRA::MeshtyingStrategyS2IElch::MeshtyingStrategyS2IElch


/*---------------------------------------------------------------------------*
 | compute time step size                                         fang 02/18 |
 *---------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2IElch::ComputeTimeStepSize(double& dt)
{
  // consider adaptive time stepping for scatra-scatra interface layer growth if necessary
  if (intlayergrowth_timestep_ > 0.)
  {
    // add state vectors to discretization
    scatratimint_->Discretization()->ClearState();
    scatratimint_->AddTimeIntegrationSpecificVectors();

    // create parameter list
    Teuchos::ParameterList condparams;

    // action for elements
    condparams.set<int>("action", SCATRA::bd_calc_elch_minmax_overpotential);

    // initialize results
    condparams.set<double>("etagrowthmin", std::numeric_limits<double>::infinity());
    condparams.set<double>("etagrowthmax", -std::numeric_limits<double>::infinity());

    // evaluate minimum and maximum interfacial overpotential associated with scatra-scatra
    // interface layer growth
    scatratimint_->Discretization()->EvaluateCondition(condparams, "S2ICouplingGrowth");
    scatratimint_->Discretization()->ClearState();

    // communicate minimum interfacial overpotential associated with scatra-scatra interface layer
    // growth
    double etagrowthmin(0.);
    scatratimint_->Discretization()->Comm().MinAll(
        &condparams.get<double>("etagrowthmin"), &etagrowthmin, 1);

    // adaptive time stepping for scatra-scatra interface layer growth is currently inactive
    if (not intlayergrowth_timestep_active_)
    {
      // check whether adaptive time stepping for scatra-scatra interface layer growth needs to be
      // activated this is the case if the minimum interfacial overpotential is currently positive,
      // but would turn negative after adding twice the change in the minimum interfacial
      // overpotential during the previous time step, i.e., eta - 2*(eta_old - eta) < 0, so that
      // lithium plating could take place after the current time step
      if (etagrowthmin > 0. and etagrowthmin - 2. * (etagrowthmin_ - etagrowthmin) < 0.)
        // activate adaptive time stepping for scatra-scatra interface layer growth
        intlayergrowth_timestep_active_ = true;
    }

    // adaptive time stepping for scatra-scatra interface layer growth is currently active
    else
    {
      // communicate maximum interfacial overpotential associated with scatra-scatra interface layer
      // growth
      double etagrowthmax(0.);
      scatratimint_->Discretization()->Comm().MaxAll(
          &condparams.get<double>("etagrowthmax"), &etagrowthmax, 1);

      // check whether maximum interfacial overpotential has become negative
      if (etagrowthmax < 0. and intlayergrowth_startstep_ < 0)
        // store current time step as indicator for completed onset of scatra-scatra interface layer
        // growth
        intlayergrowth_startstep_ = scatratimint_->Step();

      // check whether adaptive time stepping for scatra-scatra interface layer growth needs to be
      // deactivated this is the case if ten time steps have passed since the completed onset of
      // scatra-scatra interface layer growth or if the minimum interfacial overpotential is
      // positive and increasing
      if (scatratimint_->Step() == intlayergrowth_startstep_ + 10 or
          (etagrowthmin > 0. and etagrowthmin > etagrowthmin_))
      {
        // deactivate adaptive time stepping for scatra-scatra interface layer growth
        intlayergrowth_timestep_active_ = false;

        // reset time step tracker
        intlayergrowth_startstep_ = -1;
      }
    }

    // update minimum interfacial overpotential associated with scatra-scatra interface layer growth
    etagrowthmin_ = etagrowthmin;

    // reduce time step size if necessary
    if (dt > intlayergrowth_timestep_ and intlayergrowth_timestep_active_)
      dt = intlayergrowth_timestep_;
  }

  return;
}  // SCATRA::MeshtyingStrategyS2IElch::ComputeTimeStepSize


/*--------------------------------------------------------------------------------------*
 | evaluate scatra-scatra interface coupling conditions (electrochemistry)   fang 10/14 |
 *--------------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2IElch::EvaluateMeshtying()
{
  // safety check
  if (DRT::INPUT::IntegralValue<int>(*(ElchTimInt()->ElchParameterList()), "BLOCKPRECOND"))
    dserror("Block preconditioning doesn't work for scatra-scatra interface coupling yet!");

  // call base class routine
  SCATRA::MeshtyingStrategyS2I::EvaluateMeshtying();

  return;
}  // SCATRA::MeshtyingStrategyS2IElch::EvaluateMeshtying


/*----------------------------------------------------------------------------*
 | build maps associated with blocks of global system matrix       fang 06/15 |
 *----------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2IElch::BuildBlockMaps(
    const std::vector<Teuchos::RCP<DRT::Condition>>&
        partitioningconditions,                             //!< domain partitioning conditions
    std::vector<Teuchos::RCP<const Epetra_Map>>& blockmaps  //!< empty vector for maps to be built
    ) const
{
  if (matrixtype_ == INPAR::S2I::matrix_block_condition_dof)
  {
    // safety check
    if (DRT::INPUT::IntegralValue<int>(
            ElchTimInt()->ElchParameterList()->sublist("DIFFCOND"), "CURRENT_SOLUTION_VAR"))
      dserror(
          "For chosen type of global block system matrix, current must not constitute solution "
          "variable!");

    // extract number of domain partitioning conditions
    const unsigned ncond = partitioningconditions.size();

    // prepare vector for maps to be built
    blockmaps.resize(ncond * 2, Teuchos::null);

    // loop over all domain partitioning conditions
    for (unsigned icond = 0; icond < ncond; ++icond)
    {
      // initialize sets for dof IDs associated with current partitioning condition
      std::vector<std::set<int>> dofids(2);

      // extract nodes associated with current domain partitioning condition
      const std::vector<int>* nodegids = partitioningconditions[icond]->Nodes();

      // loop over all nodes associated with current domain partitioning condition
      for (unsigned inode = 0; inode < nodegids->size(); ++inode)
      {
        // extract global ID of current node
        const int nodegid = (*nodegids)[inode];

        // consider current node only if node is owned by current processor
        // need to make sure that node is stored on current processor, otherwise cannot resolve
        // "->Owner()"
        if (scatratimint_->Discretization()->HaveGlobalNode(nodegid) and
            scatratimint_->Discretization()->gNode(nodegid)->Owner() ==
                scatratimint_->Discretization()->Comm().MyPID())
        {
          // extract dof IDs associated with current node
          const std::vector<int> nodedofs = scatratimint_->Discretization()->Dof(
              0, scatratimint_->Discretization()->gNode(nodegid));

          // add concentration dof IDs to associated set
          std::copy(nodedofs.begin(), --nodedofs.end(), std::inserter(dofids[0], dofids[0].end()));

          // add electric potential dof ID to associated set
          dofids[1].insert(nodedofs.back());
        }
      }

      // transform sets for dof IDs into vectors and then into Epetra maps
      for (unsigned iset = 0; iset < 2; ++iset)
      {
        int nummyelements(0);
        int* myglobalelements(NULL);
        std::vector<int> dofidvec;
        if (dofids[iset].size() > 0)
        {
          dofidvec.reserve(dofids[iset].size());
          dofidvec.assign(dofids[iset].begin(), dofids[iset].end());
          nummyelements = dofidvec.size();
          myglobalelements = &(dofidvec[0]);
        }
        blockmaps[2 * icond + iset] =
            Teuchos::rcp(new Epetra_Map(-1, nummyelements, myglobalelements,
                scatratimint_->DofRowMap()->IndexBase(), scatratimint_->DofRowMap()->Comm()));
      }
    }
  }

  // call base class routine for other types of global system matrix
  else
    SCATRA::MeshtyingStrategyS2I::BuildBlockMaps(partitioningconditions, blockmaps);

  return;
}  // SCATRA::MeshtyingStrategyS2I::BuildBlockMaps


/*-------------------------------------------------------------------------------*
 | build null spaces associated with blocks of global system matrix   fang 07/15 |
 *-------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2IElch::BuildBlockNullSpaces() const
{
  // call base class routine
  SCATRA::MeshtyingStrategyS2I::BuildBlockNullSpaces();

  if (matrixtype_ == INPAR::S2I::matrix_block_condition_dof)
  {
    // loop over blocks of global system matrix
    for (int iblock = 0; iblock < blockmaps_->NumMaps(); ++iblock)
    {
      // store number of current block as string, starting from 1
      std::ostringstream iblockstr;
      iblockstr << iblock + 1;

      // access parameter sublist associated with smoother for current matrix block
      Teuchos::ParameterList& mueluparams = scatratimint_->Solver()
                                                ->Params()
                                                .sublist("Inverse" + iblockstr.str())
                                                .sublist("MueLu Parameters");

      // extract already reduced null space associated with current matrix block
      std::vector<double>& nullspace =
          *mueluparams.get<Teuchos::RCP<std::vector<double>>>("nullspace");

      // Each matrix block is associated with either concentration dofs or electric potential dofs
      // only. However, since the original full null space was computed for all degrees of freedom
      // on the discretization, the reduced null spaces still have the full dimension, i.e., the
      // full number of null space vectors equaling the total number of primary variables. Hence, we
      // need to decrease the dimension of each null space by one and remove the corresponding zero
      // null space vector from the null space.
      if (iblock % 2 == 0)
        // null space associated with concentration dofs
        // remove zero null space vector associated with electric potential dofs by truncating null
        // space
        nullspace.resize(blockmaps_->Map(iblock)->NumMyElements());

      else
        // null space associated with electric potential dofs
        // remove zero null space vector(s) associated with concentration dofs and retain only the
        // last null space vector associated with electric potential dofs
        nullspace.erase(
            nullspace.begin(), nullspace.end() - blockmaps_->Map(iblock)->NumMyElements());

      // decrease null space dimension and number of partial differential equations by one
      --mueluparams.get<int>("null space: dimension");
      --mueluparams.get<int>("PDE equations");
    }
  }

  return;
}  // SCATRA::MeshtyingStrategyS2IElch::BuildBlockNullSpaces


/*------------------------------------------------------------------------*
 | instantiate strategy for Newton-Raphson convergence check   fang 02/16 |
 *------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2IElch::InitConvCheckStrategy()
{
  if (couplingtype_ == INPAR::S2I::coupling_mortar_saddlepoint_petrov or
      couplingtype_ == INPAR::S2I::coupling_mortar_saddlepoint_bubnov)
    convcheckstrategy_ = Teuchos::rcp(new SCATRA::ConvCheckStrategyS2ILMElch(
        scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));
  else
    convcheckstrategy_ = Teuchos::rcp(new SCATRA::ConvCheckStrategyStdElch(
        scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));

  return;
}  // SCATRA::MeshtyingStrategyS2IElch::InitConvCheckStrategy


/*------------------------------------------------------------------------------------------*
 | update solution after convergence of the nonlinear Newton-Raphson iteration   fang 01/17 |
 *------------------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2IElch::Update() const
{
  // update scatra-scatra interface layer thicknesses in case of semi-implicit solution approach
  if (intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_semi_implicit)
  {
    // extract boundary conditions for scatra-scatra interface layer growth
    std::vector<DRT::Condition*> conditions;
    scatratimint_->Discretization()->GetCondition("S2ICouplingGrowth", conditions);

    // loop over all conditions
    for (unsigned icond = 0; icond < conditions.size(); ++icond)
    {
      // extract current condition
      const DRT::Condition* const condition = conditions[icond];

      // extract kinetic model from current condition
      switch (condition->GetInt("kinetic model"))
      {
        case INPAR::S2I::growth_kinetics_butlervolmer:
        {
          // extract parameters from current condition
          const double kr = condition->GetDouble("k_r");
          const double alphaa = condition->GetDouble("alpha_a");
          const double alphac = condition->GetDouble("alpha_c");
          const double frt = ElchTimInt()->FRT();
          const double conductivity_inverse = 1. / condition->GetDouble("conductivity");
          const double faraday =
              DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->Faraday();

          // pre-compute integration factor
          const double integrationfac(condition->GetDouble("molar mass") * scatratimint_->Dt() /
                                      (condition->GetDouble("density") * faraday));

          // extract nodal cloud from current condition
          const std::vector<int>* nodegids = condition->Nodes();

          // loop over all nodes
          for (unsigned inode = 0; inode < nodegids->size(); ++inode)
          {
            // extract global ID of current node
            const int nodegid((*nodegids)[inode]);

            // process only nodes stored by current processor
            if (scatratimint_->Discretization()->HaveGlobalNode(nodegid))
            {
              // extract current node
              const DRT::Node* const node = scatratimint_->Discretization()->gNode(nodegid);

              // process only nodes owned by current processor
              if (node->Owner() == scatratimint_->Discretization()->Comm().MyPID())
              {
                // extract local ID of first scalar transport degree of freedom associated with
                // current node
                const int doflid_scatra = scatratimint_->Discretization()->DofRowMap()->LID(
                    scatratimint_->Discretization()->Dof(
                        0, node, 0));  // Do not remove the first zero, i.e., the first function
                                       // argument, otherwise an error is thrown in debug mode!
                if (doflid_scatra < 0)
                  dserror("Couldn't extract local ID of scalar transport degree of freedom!");

                // extract local ID of scatra-scatra interface layer thickness variable associated
                // with current node
                const int doflid_growth = scatratimint_->Discretization()->DofRowMap(2)->LID(
                    scatratimint_->Discretization()->Dof(2, node, 0));
                if (doflid_growth < 0)
                  dserror("Couldn't extract local ID of scatra-scatra interface layer thickness!");

                // extract slave-side electric potential associated with current node
                const double slavepot = (*scatratimint_->Phiafnp())[doflid_scatra + 1];

                // extract master-side lithium concentration associated with current node
                const double masterphi = (*imasterphinp_)[doflid_scatra];

                // extract master-side electric potential associated with current node
                const double masterpot = (*imasterphinp_)[doflid_scatra + 1];

                // compute interface layer resistance associated with current node
                const double resistance = (*growthn_)[doflid_growth] * conductivity_inverse;

                // check existence of interface layer and set Heaviside value accordingly
                const unsigned heaviside(resistance > 0. ? 1 : 0);

                // compute exchange current density
                const double i0 = kr * faraday * pow(masterphi, alphaa);

                // compute initial guess of Butler-Volmer current density associated with lithium
                // plating, neglecting overpotential due to resistance of plated lithium
                double eta = slavepot - masterpot;
                double i = i0 * (heaviside * exp(alphaa * frt * eta) - exp(-alphac * frt * eta));

                // initialize Newton-Raphson iteration counter
                unsigned iternum(0);

                // apply Newton-Raphson method to compute Butler-Volmer current density associated
                // with lithium plating, involving overpotential due to resistance of plated lithium
                while (true)
                {
                  // increment counter
                  ++iternum;

                  // compute current Newton-Raphson residual
                  eta = slavepot - masterpot -
                        resistance *
                            i;  // open-circuit potential is zero for lithium plating reaction
                  const double expterm1 = heaviside * exp(alphaa * frt * eta);
                  const double expterm2 = exp(-alphac * frt * eta);
                  const double residual = i0 * (expterm1 - expterm2) - i;

                  // convergence check
                  if (std::abs(residual) < intlayergrowth_convtol_)
                    break;
                  else if (iternum == intlayergrowth_itemax_)
                    dserror(
                        "Local Newton-Raphson iteration for scatra-scatra interface layer growth "
                        "did not converge!");

                  // compute linearization of current Newton-Raphson residual w.r.t. Butler-Volmer
                  // current density associated with lithium plating
                  const double linearization =
                      -i0 * resistance * frt * (alphaa * expterm1 + alphac * expterm2) - 1.;

                  // update Butler-Volmer current density
                  i -= residual / linearization;
                }

                // enforce plating condition, i.e., consider initial lithium plating only in case of
                // negative overpotential
                if (!heaviside and eta >= 0.) i = 0.;

                // update lithium plating variable
                (*growthn_)[doflid_growth] -= i * integrationfac;
              }  // nodes owned by current processor
            }    // nodes stored by current processor
          }      // loop over all nodes

          break;
        }

        default:
        {
          dserror("Kinetic model for scatra-scatra interface layer growth is not yet implemented!");
          break;
        }
      }  // kinetic models
    }    // loop over all conditions
  }      // semi-implicit evaluation of scatra-scatra interface layer growth

  else
    // call base class routine
    MeshtyingStrategyS2I::Update();

  return;
}


/*----------------------------------------------------------------------*
 | singleton access method                                   fang 01/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
SCATRA::MortarCellCalcElch<distypeS, distypeM>*
SCATRA::MortarCellCalcElch<distypeS, distypeM>::Instance(
    const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
    const INPAR::S2I::InterfaceSides&
        lmside,  //!< flag for interface side underlying Lagrange multiplier definition
    const int& numdofpernode_slave,      //!< number of slave-side degrees of freedom per node
    const int& numdofpernode_master,     //!< number of master-side degrees of freedom per node
    const std::string& disname,          //!< name of mortar discretization
    const MortarCellCalcElch* delete_me  //!< pointer to instance to be deleted
)
{
  // static map assigning mortar discretization names to class instances
  static std::map<std::string, MortarCellCalcElch<distypeS, distypeM>*> instances;

  // create new instance or return existing one
  if (!delete_me)
  {
    // create new instance if not yet available
    if (instances.find(disname) == instances.end())
      instances[disname] = new MortarCellCalcElch<distypeS, distypeM>(
          couplingtype, lmside, numdofpernode_slave, numdofpernode_master);
  }

  // delete existing instance
  else
  {
    // loop over all existing instances
    for (typename std::map<std::string, MortarCellCalcElch<distypeS, distypeM>*>::iterator i =
             instances.begin();
         i != instances.end(); ++i)
    {
      // check whether current instance should be deleted
      if (i->second == delete_me)
      {
        // delete current instance
        delete i->second;

        // remove deleted instance from map
        instances.erase(i);

        // return null pointer
        return NULL;
      }
    }

    // catch internal error
    dserror("Instance to be deleted couldn't be found in static map!");
  }

  // return existing instance
  return instances[disname];
}


/*----------------------------------------------------------------------*
 | singleton destruction                                     fang 01/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalcElch<distypeS, distypeM>::Done()
{
  // delete singleton
  Instance(INPAR::S2I::coupling_undefined, INPAR::S2I::side_undefined, 0, 0, "", this);

  return;
}


/*----------------------------------------------------------------------*
 | protected constructor for singletons                      fang 01/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
SCATRA::MortarCellCalcElch<distypeS, distypeM>::MortarCellCalcElch(
    const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
    const INPAR::S2I::InterfaceSides&
        lmside,  //!< flag for interface side underlying Lagrange multiplier definition
    const int& numdofpernode_slave,  //!< number of slave-side degrees of freedom per node
    const int& numdofpernode_master  //!< number of master-side degrees of freedom per node
    )
    : my::MortarCellCalc(couplingtype, lmside, numdofpernode_slave, numdofpernode_master)
{
  return;
}


/*---------------------------------------------------------------------------*
 | evaluate and assemble interface linearizations and residuals   fang 01/16 |
 *---------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalcElch<distypeS, distypeM>::EvaluateCondition(
    const DRT::Discretization& idiscret,     //!< interface discretization
    MORTAR::IntCell& cell,                   //!< mortar integration cell
    MORTAR::MortarElement& slaveelement,     //!< slave-side mortar element
    MORTAR::MortarElement& masterelement,    //!< master-side mortar element
    DRT::Element::LocationArray& la_slave,   //!< slave-side location array
    DRT::Element::LocationArray& la_master,  //!< master-side location array
    const Teuchos::ParameterList& params,    //!< parameter list
    Epetra_SerialDenseMatrix&
        k_ss,  //!< linearizations of slave-side residuals w.r.t. slave-side dofs
    Epetra_SerialDenseMatrix&
        k_sm,  //!< linearizations of slave-side residuals w.r.t. master-side dofs
    Epetra_SerialDenseMatrix&
        k_ms,  //!< linearizations of master-side residuals w.r.t. slave-side dofs
    Epetra_SerialDenseMatrix&
        k_mm,  //!< linearizations of master-side residuals w.r.t. master-side dofs
    Epetra_SerialDenseVector& r_s,  //!< slave-side residual vector
    Epetra_SerialDenseVector& r_m   //!< master-side residual vector
)
{
  // safety checks
  if (my::numdofpernode_slave_ != 2 or my::numdofpernode_master_ != 2)
    dserror("Invalid number of degrees of freedom per node!");
  if (DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->EquPot() !=
      INPAR::ELCH::equpot_divi)
    dserror("Invalid closing equation for electric potential!");

  // extract condition from parameter list
  DRT::Condition* condition = params.get<DRT::Condition*>("condition");
  if (condition == NULL) dserror("Cannot access scatra-scatra interface coupling condition!");

  // access material of slave element
  Teuchos::RCP<const MAT::Electrode> matelectrode =
      Teuchos::rcp_dynamic_cast<const MAT::Electrode>(slaveelement.Material());
  if (matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // extract nodal state variables associated with slave and master elements
  this->ExtractNodeValues(idiscret, la_slave, la_master);

  // determine quadrature rule
  const DRT::UTILS::IntPointsAndWeights<2> intpoints(DRT::UTILS::intrule_tri_7point);

  const int kineticmodel = my::scatraparamsboundary_->KineticModel();
  const int numelectrons = my::scatraparamsboundary_->NumElectrons();
  const std::vector<int>* stoichiometries = my::scatraparamsboundary_->Stoichiometries();
  const double kr = my::scatraparamsboundary_->Kr();
  const double alphaa = my::scatraparamsboundary_->AlphaA();
  const double alphac = my::scatraparamsboundary_->AlphaC();
  const double resistance = my::scatraparamsboundary_->Resistance();
  const double itemaxmodifiedBV = my::scatraparamsboundary_->ItemaxmodifiedBV();
  const double convtolmodifiedBV = my::scatraparamsboundary_->ConvtolmodifiedBV();

  // loop over all integration points
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    // evaluate shape functions and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndDomIntFacAtIntPoint(
        slaveelement, masterelement, cell, intpoints, iquad);

    // overall integration factors
    const double timefacfac =
        DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->TimeFac() * fac;
    const double timefacrhsfac =
        DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->TimeFacRhs() * fac;
    if (timefacfac < 0. or timefacrhsfac < 0.) dserror("Integration factor is negative!");

    DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<
        distypeS>::template EvaluateS2ICouplingAtIntegrationPoint<distypeM>(matelectrode,
        my::ephinp_slave_, my::ephinp_master_, my::funct_slave_, my::funct_master_,
        my::test_lm_slave_, my::test_lm_master_, kineticmodel, numelectrons, stoichiometries, kr,
        alphaa, alphac, resistance, itemaxmodifiedBV, convtolmodifiedBV, timefacfac, timefacrhsfac,
        GetFRT(), k_ss, k_sm, k_ms, k_mm, r_s, r_m);
  }

  return;
}


/*--------------------------------------------------------------------------------------------------------*
 | evaluate and assemble interface linearizations and residuals for node-to-segment coupling   fang
 08/16 |
 *--------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalcElch<distypeS, distypeM>::EvaluateConditionNTS(
    DRT::Condition& condition,            //!< scatra-scatra interface coupling condition
    const MORTAR::MortarNode& slavenode,  //!< slave-side node
    const double& lumpedarea,  //!< lumped interface area fraction associated with slave-side node
    MORTAR::MortarElement& slaveelement,   //!< slave-side mortar element
    MORTAR::MortarElement& masterelement,  //!< master-side mortar element
    const std::vector<LINALG::Matrix<my::nen_slave_, 1>>&
        ephinp_slave,  //!< state variables at slave-side nodes
    const std::vector<LINALG::Matrix<my::nen_master_, 1>>&
        ephinp_master,  //!< state variables at master-side nodes
    Epetra_SerialDenseMatrix&
        k_ss,  //!< linearizations of slave-side residuals w.r.t. slave-side dofs
    Epetra_SerialDenseMatrix&
        k_sm,  //!< linearizations of slave-side residuals w.r.t. master-side dofs
    Epetra_SerialDenseMatrix&
        k_ms,  //!< linearizations of master-side residuals w.r.t. slave-side dofs
    Epetra_SerialDenseMatrix&
        k_mm,  //!< linearizations of master-side residuals w.r.t. master-side dofs
    Epetra_SerialDenseVector& r_s,  //!< slave-side residual vector
    Epetra_SerialDenseVector& r_m   //!< master-side residual vector
)
{
  // safety checks
  if (my::numdofpernode_slave_ != 2 or my::numdofpernode_master_ != 2)
    dserror("Invalid number of degrees of freedom per node!");
  if (DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->EquPot() !=
      INPAR::ELCH::equpot_divi)
    dserror("Invalid closing equation for electric potential!");

  // access material of slave element
  Teuchos::RCP<const MAT::Electrode> matelectrode = Teuchos::rcp_dynamic_cast<const MAT::Electrode>(
      Teuchos::rcp_dynamic_cast<DRT::FaceElement>(condition.Geometry()[slaveelement.Id()])
          ->ParentElement()
          ->Material());
  if (matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // evaluate shape functions at position of slave-side node
  my::EvalShapeFuncAtSlaveNode(slavenode, slaveelement, masterelement);

  const int kineticmodel = my::scatraparamsboundary_->KineticModel();
  const int numelectrons = my::scatraparamsboundary_->NumElectrons();
  const std::vector<int>* stoichiometries = my::scatraparamsboundary_->Stoichiometries();
  const double kr = my::scatraparamsboundary_->Kr();
  const double alphaa = my::scatraparamsboundary_->AlphaA();
  const double alphac = my::scatraparamsboundary_->AlphaC();
  const double resistance = my::scatraparamsboundary_->Resistance();
  const double itemaxmodifiedBV = my::scatraparamsboundary_->ItemaxmodifiedBV();
  const double convtolmodifiedBV = my::scatraparamsboundary_->ConvtolmodifiedBV();

  // overall integration factors
  const double timefacfac =
      DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->TimeFac() * lumpedarea;
  const double timefacrhsfac =
      DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->TimeFacRhs() * lumpedarea;
  if (timefacfac < 0. or timefacrhsfac < 0.) dserror("Integration factor is negative!");

  DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<
      distypeS>::template EvaluateS2ICouplingAtIntegrationPoint<distypeM>(matelectrode,
      ephinp_slave, ephinp_master, my::funct_slave_, my::funct_master_, my::funct_slave_,
      my::funct_master_, kineticmodel, numelectrons, stoichiometries, kr, alphaa, alphac,
      resistance, itemaxmodifiedBV, convtolmodifiedBV, timefacfac, timefacrhsfac,
      DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->FRT(), k_ss, k_sm, k_ms, k_mm, r_s,
      r_m);

  return;
}


/*----------------------------------------------------------------------*
 | evaluate factor F/RT                                      fang 01/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
double SCATRA::MortarCellCalcElch<distypeS, distypeM>::GetFRT() const
{
  // fetch factor F/RT from electrochemistry parameter list
  return DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->FRT();
};


/*----------------------------------------------------------------------*
 | singleton access method                                   fang 01/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
SCATRA::MortarCellCalcElchSTIThermo<distypeS, distypeM>*
SCATRA::MortarCellCalcElchSTIThermo<distypeS, distypeM>::Instance(
    const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
    const INPAR::S2I::InterfaceSides&
        lmside,  //!< flag for interface side underlying Lagrange multiplier definition
    const int& numdofpernode_slave,   //!< number of slave-side degrees of freedom per node
    const int& numdofpernode_master,  //!< number of master-side degrees of freedom per node
    const std::string& disname,       //!< name of mortar discretization
    const MortarCellCalcElchSTIThermo* delete_me  //!< pointer to instance to be deleted
)
{
  // static map assigning mortar discretization names to class instances
  static std::map<std::string, MortarCellCalcElchSTIThermo<distypeS, distypeM>*> instances;

  // create new instance or return existing one
  if (!delete_me)
  {
    // create new instance if not yet available
    if (instances.find(disname) == instances.end())
      instances[disname] = new MortarCellCalcElchSTIThermo<distypeS, distypeM>(
          couplingtype, lmside, numdofpernode_slave, numdofpernode_master);
  }

  // delete existing instance
  else
  {
    // loop over all existing instances
    for (typename std::map<std::string, MortarCellCalcElchSTIThermo<distypeS, distypeM>*>::iterator
             i = instances.begin();
         i != instances.end(); ++i)
    {
      // check whether current instance should be deleted
      if (i->second == delete_me)
      {
        // delete current instance
        delete i->second;

        // remove deleted instance from map
        instances.erase(i);

        // return null pointer
        return NULL;
      }
    }

    // catch internal error
    dserror("Instance to be deleted couldn't be found in static map!");
  }

  // return existing instance
  return instances[disname];
}

/*----------------------------------------------------------------------*
 | singleton destruction                                     fang 01/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalcElchSTIThermo<distypeS, distypeM>::Done()
{
  // delete singleton
  Instance(INPAR::S2I::coupling_undefined, INPAR::S2I::side_undefined, 0, 0, "", this);

  return;
};


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 01/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
SCATRA::MortarCellCalcElchSTIThermo<distypeS, distypeM>::MortarCellCalcElchSTIThermo(
    const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
    const INPAR::S2I::InterfaceSides&
        lmside,  //!< flag for interface side underlying Lagrange multiplier definition
    const int& numdofpernode_slave,  //!< number of slave-side degrees of freedom per node
    const int& numdofpernode_master  //!< number of master-side degrees of freedom per node
    )
    :  // call base class constructor
      myelch::MortarCellCalcElch(couplingtype, lmside, numdofpernode_slave, numdofpernode_master),

      // initialize member variable
      etempnp_slave_(true)
{
  return;
}


/*--------------------------------------------------------------------------------------------------------------------*
 | evaluate single mortar integration cell of particular slave-side and master-side discretization
 types   fang 01/17 |
 *--------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalcElchSTIThermo<distypeS, distypeM>::Evaluate(
    const DRT::Discretization& idiscret,     //!< interface discretization
    MORTAR::IntCell& cell,                   //!< mortar integration cell
    MORTAR::MortarElement& slaveelement,     //!< slave-side mortar element
    MORTAR::MortarElement& masterelement,    //!< master-side mortar element
    DRT::Element::LocationArray& la_slave,   //!< slave-side location array
    DRT::Element::LocationArray& la_master,  //!< master-side location array
    const Teuchos::ParameterList& params,    //!< parameter list
    Epetra_SerialDenseMatrix& cellmatrix1,   //!< cell matrix 1
    Epetra_SerialDenseMatrix& cellmatrix2,   //!< cell matrix 2
    Epetra_SerialDenseMatrix& cellmatrix3,   //!< cell matrix 3
    Epetra_SerialDenseMatrix& cellmatrix4,   //!< cell matrix 4
    Epetra_SerialDenseVector& cellvector1,   //!< cell vector 1
    Epetra_SerialDenseVector& cellvector2    //!< cell vector 2
)
{
  // extract and evaluate action
  switch (DRT::INPUT::get<INPAR::S2I::EvaluationActions>(params, "action"))
  {
    // evaluate and assemble off-diagonal interface linearizations
    case INPAR::S2I::evaluate_condition_od:
    {
      EvaluateConditionOD(idiscret, cell, slaveelement, masterelement, la_slave, la_master, params,
          cellmatrix1, cellmatrix3);

      break;
    }

    // call base class routine
    default:
    {
      my::Evaluate(idiscret, cell, slaveelement, masterelement, la_slave, la_master, params,
          cellmatrix1, cellmatrix2, cellmatrix3, cellmatrix4, cellvector1, cellvector2);

      break;
    }
  }

  return;
}


/*---------------------------------------------------------------------------*
 | evaluate and assemble off-diagonal interface linearizations    fang 01/17 |
 *---------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalcElchSTIThermo<distypeS, distypeM>::EvaluateConditionOD(
    const DRT::Discretization& idiscret,     //!< interface discretization
    MORTAR::IntCell& cell,                   //!< mortar integration cell
    MORTAR::MortarElement& slaveelement,     //!< slave-side mortar element
    MORTAR::MortarElement& masterelement,    //!< master-side mortar element
    DRT::Element::LocationArray& la_slave,   //!< slave-side location array
    DRT::Element::LocationArray& la_master,  //!< master-side location array
    const Teuchos::ParameterList& params,    //!< parameter list
    Epetra_SerialDenseMatrix&
        k_ss,  //!< linearizations of slave-side residuals w.r.t. slave-side dofs
    Epetra_SerialDenseMatrix&
        k_ms  //!< linearizations of master-side residuals w.r.t. slave-side dofs
)
{
  // safety checks
  if (my::numdofpernode_slave_ != 2 or my::numdofpernode_master_ != 2)
    dserror("Invalid number of degrees of freedom per node!");
  if (DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->EquPot() !=
      INPAR::ELCH::equpot_divi)
    dserror("Invalid closing equation for electric potential!");

  // extract condition from parameter list
  DRT::Condition* s2icondition = params.get<DRT::Condition*>("condition");
  if (s2icondition == NULL) dserror("Cannot access scatra-scatra interface coupling condition!");

  // access material of slave element
  Teuchos::RCP<const MAT::Electrode> matelectrode = Teuchos::rcp_dynamic_cast<const MAT::Electrode>(
      Teuchos::rcp_dynamic_cast<DRT::FaceElement>(s2icondition->Geometry()[slaveelement.Id()])
          ->ParentElement()
          ->Material());
  if (matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // extract nodal state variables associated with slave and master elements
  ExtractNodeValues(idiscret, la_slave, la_master);

  // determine quadrature rule
  const DRT::UTILS::IntPointsAndWeights<2> intpoints(DRT::UTILS::intrule_tri_7point);

  const int kineticmodel = my::scatraparamsboundary_->KineticModel();
  const int numelectrons = my::scatraparamsboundary_->NumElectrons();
  const std::vector<int>* stoichiometries = my::scatraparamsboundary_->Stoichiometries();
  const double kr = my::scatraparamsboundary_->Kr();
  const double alphaa = my::scatraparamsboundary_->AlphaA();
  const double alphac = my::scatraparamsboundary_->AlphaC();

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac =
        my::EvalShapeFuncAndDomIntFacAtIntPoint(slaveelement, masterelement, cell, intpoints, gpid);

    // evaluate overall integration factor
    const double timefacfac =
        DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->TimeFac() * fac;
    if (timefacfac < 0.) dserror("Integration factor is negative!");

    DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<
        distypeS>::template EvaluateS2ICouplingODAtIntegrationPoint<distypeM>(matelectrode,
        my::ephinp_slave_, etempnp_slave_, my::ephinp_master_, my::funct_slave_, my::funct_master_,
        my::test_lm_slave_, my::test_lm_master_, kineticmodel, numelectrons, stoichiometries, kr,
        alphaa, alphac, timefacfac, k_ss, k_ms);
  }  // loop over integration points

  return;
}


/*------------------------------------------------------------------------------------*
 | extract nodal state variables associated with mortar integration cell   fang 01/17 |
 *------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalcElchSTIThermo<distypeS, distypeM>::ExtractNodeValues(
    const DRT::Discretization& idiscret,    //!< interface discretization
    DRT::Element::LocationArray& la_slave,  //!< slave-side location array
    DRT::Element::LocationArray& la_master  //!< master-side location array
)
{
  // call base class routine
  my::ExtractNodeValues(idiscret, la_slave, la_master);

  // extract nodal temperature variables associated with mortar integration cell
  my::ExtractNodeValues(etempnp_slave_, idiscret, la_slave, "thermo", 1);

  return;
}


/*----------------------------------------------------------------------*
 | evaluate factor F/RT                                      fang 01/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
double SCATRA::MortarCellCalcElchSTIThermo<distypeS, distypeM>::GetFRT() const
{
  // evaluate local temperature value
  const double temperature = my::funct_slave_.Dot(etempnp_slave_);

  // safety check
  if (temperature <= 0.) dserror("Temperature is non-positive!");

  const double faraday = DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->Faraday();
  const double gasconstant =
      DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->GasConstant();

  // evaluate factor F/RT
  return faraday / (gasconstant * temperature);
}


/*----------------------------------------------------------------------*
 | singleton access method                                   fang 01/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
SCATRA::MortarCellCalcSTIElch<distypeS, distypeM>*
SCATRA::MortarCellCalcSTIElch<distypeS, distypeM>::Instance(
    const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
    const INPAR::S2I::InterfaceSides&
        lmside,  //!< flag for interface side underlying Lagrange multiplier definition
    const int& numdofpernode_slave,         //!< number of slave-side degrees of freedom per node
    const int& numdofpernode_master,        //!< number of master-side degrees of freedom per node
    const std::string& disname,             //!< name of mortar discretization
    const MortarCellCalcSTIElch* delete_me  //!< pointer to instance to be deleted
)
{
  // static map assigning mortar discretization names to class instances
  static std::map<std::string, MortarCellCalcSTIElch<distypeS, distypeM>*> instances;

  // create new instance or return existing one
  if (!delete_me)
  {
    // create new instance if not yet available
    if (instances.find(disname) == instances.end())
      instances[disname] = new MortarCellCalcSTIElch<distypeS, distypeM>(
          couplingtype, lmside, numdofpernode_slave, numdofpernode_master);
  }

  // delete existing instance
  else
  {
    // loop over all existing instances
    for (typename std::map<std::string, MortarCellCalcSTIElch<distypeS, distypeM>*>::iterator i =
             instances.begin();
         i != instances.end(); ++i)
    {
      // check whether current instance should be deleted
      if (i->second == delete_me)
      {
        // delete current instance
        delete i->second;

        // remove deleted instance from map
        instances.erase(i);

        // return null pointer
        return NULL;
      }
    }

    // catch internal error
    dserror("Instance to be deleted couldn't be found in static map!");
  }

  // return existing instance
  return instances[disname];
}


/*----------------------------------------------------------------------*
 | singleton destruction                                     fang 01/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalcSTIElch<distypeS, distypeM>::Done()
{
  // delete singleton
  Instance(INPAR::S2I::coupling_undefined, INPAR::S2I::side_undefined, 0, 0, "", this);

  return;
};


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 01/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
SCATRA::MortarCellCalcSTIElch<distypeS, distypeM>::MortarCellCalcSTIElch(
    const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
    const INPAR::S2I::InterfaceSides&
        lmside,  //!< flag for interface side underlying Lagrange multiplier definition
    const int& numdofpernode_slave,  //!< number of slave-side degrees of freedom per node
    const int& numdofpernode_master  //!< number of master-side degrees of freedom per node
    )
    :  // call base class constructor
      my::MortarCellCalc(couplingtype, lmside, numdofpernode_slave, numdofpernode_master),

      // initialize member variables
      eelchnp_slave_(2, LINALG::Matrix<my::nen_slave_, 1>(true)),
      eelchnp_master_(2, LINALG::Matrix<my::nen_master_, 1>(true))
{
  return;
}


/*--------------------------------------------------------------------------------------------------------------------*
 | evaluate single mortar integration cell of particular slave-side and master-side discretization
 types   fang 01/17 |
 *--------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalcSTIElch<distypeS, distypeM>::Evaluate(
    const DRT::Discretization& idiscret,     //!< interface discretization
    MORTAR::IntCell& cell,                   //!< mortar integration cell
    MORTAR::MortarElement& slaveelement,     //!< slave-side mortar element
    MORTAR::MortarElement& masterelement,    //!< master-side mortar element
    DRT::Element::LocationArray& la_slave,   //!< slave-side location array
    DRT::Element::LocationArray& la_master,  //!< master-side location array
    const Teuchos::ParameterList& params,    //!< parameter list
    Epetra_SerialDenseMatrix& cellmatrix1,   //!< cell matrix 1
    Epetra_SerialDenseMatrix& cellmatrix2,   //!< cell matrix 2
    Epetra_SerialDenseMatrix& cellmatrix3,   //!< cell matrix 3
    Epetra_SerialDenseMatrix& cellmatrix4,   //!< cell matrix 4
    Epetra_SerialDenseVector& cellvector1,   //!< cell vector 1
    Epetra_SerialDenseVector& cellvector2    //!< cell vector 2
)
{
  // extract and evaluate action
  switch (DRT::INPUT::get<INPAR::S2I::EvaluationActions>(params, "action"))
  {
    // evaluate and assemble interface linearizations and residuals
    case INPAR::S2I::evaluate_condition:
    {
      EvaluateCondition(idiscret, cell, slaveelement, masterelement, la_slave, la_master, params,
          cellmatrix1, cellvector1);

      break;
    }

    // evaluate and assemble off-diagonal interface linearizations
    case INPAR::S2I::evaluate_condition_od:
    {
      EvaluateConditionOD(idiscret, cell, slaveelement, masterelement, la_slave, la_master, params,
          cellmatrix1, cellmatrix2);

      break;
    }

    // call base class routine
    default:
    {
      my::Evaluate(idiscret, cell, slaveelement, masterelement, la_slave, la_master, params,
          cellmatrix1, cellmatrix2, cellmatrix3, cellmatrix4, cellvector1, cellvector2);

      break;
    }
  }

  return;
}


/*---------------------------------------------------------------------------*
 | evaluate and assemble interface linearizations and residuals   fang 01/17 |
 *---------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalcSTIElch<distypeS, distypeM>::EvaluateCondition(
    const DRT::Discretization& idiscret,     //!< interface discretization
    MORTAR::IntCell& cell,                   //!< mortar integration cell
    MORTAR::MortarElement& slaveelement,     //!< slave-side mortar element
    MORTAR::MortarElement& masterelement,    //!< master-side mortar element
    DRT::Element::LocationArray& la_slave,   //!< slave-side location array
    DRT::Element::LocationArray& la_master,  //!< master-side location array
    const Teuchos::ParameterList& params,    //!< parameter list
    Epetra_SerialDenseMatrix&
        k_ss,  //!< linearizations of slave-side residuals w.r.t. slave-side dofs
    Epetra_SerialDenseVector& r_s  //!< slave-side residual vector
)
{
  // safety check
  if (my::numdofpernode_slave_ != 1 or my::numdofpernode_master_ != 1)
    dserror("Invalid number of degrees of freedom per node!");

  // extract condition from parameter list
  DRT::Condition* s2icondition = params.get<DRT::Condition*>("condition");
  if (s2icondition == NULL) dserror("Cannot access scatra-scatra interface coupling condition!");

  // access primary and secondary materials of slave element
  const Teuchos::RCP<const MAT::Soret> matsoret = Teuchos::rcp_dynamic_cast<const MAT::Soret>(
      Teuchos::rcp_dynamic_cast<DRT::FaceElement>(s2icondition->Geometry()[slaveelement.Id()])
          ->ParentElement()
          ->Material());
  const Teuchos::RCP<const MAT::Electrode> matelectrode =
      Teuchos::rcp_dynamic_cast<const MAT::Electrode>(
          Teuchos::rcp_dynamic_cast<DRT::FaceElement>(s2icondition->Geometry()[slaveelement.Id()])
              ->ParentElement()
              ->Material(1));
  if (matsoret == Teuchos::null or matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // extract nodal state variables associated with slave and master elements
  ExtractNodeValues(idiscret, la_slave, la_master);

  // determine quadrature rule
  const DRT::UTILS::IntPointsAndWeights<2> intpoints(DRT::UTILS::intrule_tri_7point);

  const int kineticmodel = my::scatraparamsboundary_->KineticModel();
  const double kr = my::scatraparamsboundary_->Kr();
  const double alphaa = my::scatraparamsboundary_->AlphaA();
  const double alphac = my::scatraparamsboundary_->AlphaC();
  const double peltier = my::scatraparamsboundary_->Peltier();

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac =
        my::EvalShapeFuncAndDomIntFacAtIntPoint(slaveelement, masterelement, cell, intpoints, gpid);

    // evaluate overall integration factors
    const double timefacfac =
        DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("thermo")->TimeFac() * fac;
    const double timefacrhsfac =
        DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("thermo")->TimeFacRhs() * fac;
    if (timefacfac < 0. or timefacrhsfac < 0.) dserror("Integration factor is negative!");

    DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<
        distypeS>::template EvaluateS2ICouplingAtIntegrationPoint<distypeM>(matelectrode,
        my::ephinp_slave_[0], eelchnp_slave_, eelchnp_master_, my::funct_slave_, my::funct_master_,
        kineticmodel, kr, alphaa, alphac, peltier, timefacfac, timefacrhsfac, k_ss, r_s);
  }  // loop over integration points

  return;
}


/*---------------------------------------------------------------------------*
 | evaluate and assemble off-diagonal interface linearizations    fang 01/17 |
 *---------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalcSTIElch<distypeS, distypeM>::EvaluateConditionOD(
    const DRT::Discretization& idiscret,     //!< interface discretization
    MORTAR::IntCell& cell,                   //!< mortar integration cell
    MORTAR::MortarElement& slaveelement,     //!< slave-side mortar element
    MORTAR::MortarElement& masterelement,    //!< master-side mortar element
    DRT::Element::LocationArray& la_slave,   //!< slave-side location array
    DRT::Element::LocationArray& la_master,  //!< master-side location array
    const Teuchos::ParameterList& params,    //!< parameter list
    Epetra_SerialDenseMatrix&
        k_ss,  //!< linearizations of slave-side residuals w.r.t. slave-side dofs
    Epetra_SerialDenseMatrix&
        k_sm  //!< linearizations of slave-side residuals w.r.t. master-side dofs
)
{
  // safety check
  if (my::numdofpernode_slave_ != 1 or my::numdofpernode_master_ != 1)
    dserror("Invalid number of degrees of freedom per node!");

  // extract condition from parameter list
  DRT::Condition* s2icondition = params.get<DRT::Condition*>("condition");
  if (s2icondition == NULL) dserror("Cannot access scatra-scatra interface coupling condition!");

  // access primary and secondary materials of parent element
  Teuchos::RCP<const MAT::Soret> matsoret = Teuchos::rcp_dynamic_cast<const MAT::Soret>(
      Teuchos::rcp_dynamic_cast<DRT::FaceElement>(s2icondition->Geometry()[slaveelement.Id()])
          ->ParentElement()
          ->Material());
  Teuchos::RCP<const MAT::Electrode> matelectrode = Teuchos::rcp_dynamic_cast<const MAT::Electrode>(
      Teuchos::rcp_dynamic_cast<DRT::FaceElement>(s2icondition->Geometry()[slaveelement.Id()])
          ->ParentElement()
          ->Material(1));
  if (matsoret == Teuchos::null or matelectrode == Teuchos::null)
    dserror("Invalid electrode or soret material for scatra-scatra interface coupling!");

  // extract nodal state variables associated with slave and master elements
  ExtractNodeValues(idiscret, la_slave, la_master);

  // determine quadrature rule
  const DRT::UTILS::IntPointsAndWeights<2> intpoints(DRT::UTILS::intrule_tri_7point);

  const int kineticmodel = my::scatraparamsboundary_->KineticModel();
  const double kr = my::scatraparamsboundary_->Kr();
  const double alphaa = my::scatraparamsboundary_->AlphaA();
  const double alphac = my::scatraparamsboundary_->AlphaC();
  const double peltier = my::scatraparamsboundary_->Peltier();

  // loop over all integration points
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    // evaluate shape functions and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndDomIntFacAtIntPoint(
        slaveelement, masterelement, cell, intpoints, iquad);

    // overall integration factors
    const double timefacfac =
        DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("thermo")->TimeFac() * fac;
    if (timefacfac < 0.) dserror("Integration factor is negative!");

    DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<
        distypeS>::template EvaluateS2ICouplingODAtIntegrationPoint<distypeM>(matelectrode,
        my::ephinp_slave_[0], eelchnp_slave_, eelchnp_master_, my::funct_slave_, my::funct_master_,
        kineticmodel, kr, alphaa, alphac, peltier, timefacfac, k_ss, k_sm);
  }  // loop over integration points

  return;
}


/*------------------------------------------------------------------------------------*
 | extract nodal state variables associated with mortar integration cell   fang 01/17 |
 *------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distypeS, DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalcSTIElch<distypeS, distypeM>::ExtractNodeValues(
    const DRT::Discretization& idiscret,    //!< interface discretization
    DRT::Element::LocationArray& la_slave,  //!< slave-side location array
    DRT::Element::LocationArray& la_master  //!< master-side location array
)
{
  // extract nodal temperature variables associated with slave element
  my::ExtractNodeValues(my::ephinp_slave_[0], idiscret, la_slave);

  // extract nodal electrochemistry variables associated with mortar integration cell
  my::ExtractNodeValues(
      eelchnp_slave_, eelchnp_master_, idiscret, la_slave, la_master, "scatra", 1);

  return;
}


// forward declarations
template class SCATRA::MortarCellCalcElch<DRT::Element::tri3, DRT::Element::tri3>;
template class SCATRA::MortarCellCalcElch<DRT::Element::tri3, DRT::Element::quad4>;
template class SCATRA::MortarCellCalcElch<DRT::Element::quad4, DRT::Element::tri3>;
template class SCATRA::MortarCellCalcElch<DRT::Element::quad4, DRT::Element::quad4>;
template class SCATRA::MortarCellCalcElchSTIThermo<DRT::Element::tri3, DRT::Element::tri3>;
template class SCATRA::MortarCellCalcElchSTIThermo<DRT::Element::tri3, DRT::Element::quad4>;
template class SCATRA::MortarCellCalcElchSTIThermo<DRT::Element::quad4, DRT::Element::tri3>;
template class SCATRA::MortarCellCalcElchSTIThermo<DRT::Element::quad4, DRT::Element::quad4>;
template class SCATRA::MortarCellCalcSTIElch<DRT::Element::tri3, DRT::Element::tri3>;
template class SCATRA::MortarCellCalcSTIElch<DRT::Element::tri3, DRT::Element::quad4>;
template class SCATRA::MortarCellCalcSTIElch<DRT::Element::quad4, DRT::Element::tri3>;
template class SCATRA::MortarCellCalcSTIElch<DRT::Element::quad4, DRT::Element::quad4>;
