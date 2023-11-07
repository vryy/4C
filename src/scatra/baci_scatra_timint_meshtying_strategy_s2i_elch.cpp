/*----------------------------------------------------------------------*/
/*! \file

\brief Scatra-scatra interface coupling strategy for electrochemistry problems

\level 2


*----------------------------------------------------------------------*/
#include "baci_scatra_timint_meshtying_strategy_s2i_elch.H"

#include "baci_coupling_adapter.H"
#include "baci_lib_discret.H"
#include "baci_lib_globalproblem.H"
#include "baci_lib_utils_gid_vector.H"
#include "baci_lib_utils_parameter_list.H"
#include "baci_linalg_mapextractor.H"
#include "baci_linalg_sparseoperator.H"
#include "baci_linear_solver_method_linalg.H"
#include "baci_mat_electrode.H"
#include "baci_mat_soret.H"
#include "baci_mortar_element.H"
#include "baci_scatra_ele_boundary_calc_elch_electrode.H"
#include "baci_scatra_ele_boundary_calc_elch_electrode_sti_thermo.H"
#include "baci_scatra_ele_boundary_calc_elch_electrode_utils.H"
#include "baci_scatra_ele_boundary_calc_sti_electrode.H"
#include "baci_scatra_ele_parameter_boundary.H"
#include "baci_scatra_ele_parameter_elch.H"
#include "baci_scatra_ele_parameter_timint.H"
#include "baci_utils_singleton_owner.H"

/*----------------------------------------------------------------------*
 | constructor                                               fang 12/14 |
 *----------------------------------------------------------------------*/
SCATRA::MeshtyingStrategyS2IElch::MeshtyingStrategyS2IElch(
    SCATRA::ScaTraTimIntElch* elchtimint, const Teuchos::ParameterList& parameters)
    : MeshtyingStrategyS2I(elchtimint, parameters),
      etagrowthmin_(0.),
      intlayergrowth_startstep_(-1),
      intlayergrowth_timestep_active_(false)
{
}  // SCATRA::MeshtyingStrategyS2IElch::MeshtyingStrategyS2IElch


/*---------------------------------------------------------------------------*
 | compute time step size                                         fang 02/18 |
 *---------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2IElch::ComputeTimeStepSize(double& dt)
{
  // consider adaptive time stepping for scatra-scatra interface layer growth if necessary
  if (intlayergrowth_timestep_ > 0.0)
  {
    // add state vectors to discretization
    scatratimint_->Discretization()->ClearState();
    scatratimint_->AddTimeIntegrationSpecificVectors();

    // create parameter list
    Teuchos::ParameterList condparams;

    // action for elements
    DRT::UTILS::AddEnumClassToParameterList<SCATRA::BoundaryAction>(
        "action", SCATRA::BoundaryAction::calc_elch_minmax_overpotential, condparams);

    // initialize results
    condparams.set<double>("etagrowthmin", std::numeric_limits<double>::infinity());
    condparams.set<double>("etagrowthmax", -std::numeric_limits<double>::infinity());

    // extract boundary conditions for scatra-scatra interface layer growth
    std::vector<DRT::Condition*> conditions;
    scatratimint_->Discretization()->GetCondition("S2ICouplingGrowth", conditions);

    // collect condition specific data and store to scatra boundary parameter class
    SetConditionSpecificScaTraParameters(*conditions[0]);
    // evaluate minimum and maximum interfacial overpotential associated with scatra-scatra
    // interface layer growth
    scatratimint_->Discretization()->EvaluateCondition(condparams, Teuchos::null, Teuchos::null,
        Teuchos::null, Teuchos::null, Teuchos::null, "S2ICouplingGrowth");
    scatratimint_->Discretization()->ClearState();

    // communicate minimum interfacial overpotential associated with scatra-scatra interface layer
    // growth
    double etagrowthmin(0.0);
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
      if (etagrowthmin > 0.0 and etagrowthmin - 2.0 * (etagrowthmin_ - etagrowthmin) < 0.0)
        // activate adaptive time stepping for scatra-scatra interface layer growth
        intlayergrowth_timestep_active_ = true;
    }

    // adaptive time stepping for scatra-scatra interface layer growth is currently active
    else
    {
      // communicate maximum interfacial overpotential associated with scatra-scatra interface layer
      // growth
      double etagrowthmax(0.0);
      scatratimint_->Discretization()->Comm().MaxAll(
          &condparams.get<double>("etagrowthmax"), &etagrowthmax, 1);

      // check whether maximum interfacial overpotential has become negative
      if (etagrowthmax < 0.0 and intlayergrowth_startstep_ < 0)
      {
        // store current time step as indicator for completed onset of scatra-scatra interface layer
        // growth
        intlayergrowth_startstep_ = scatratimint_->Step();
      }

      // check whether adaptive time stepping for scatra-scatra interface layer growth needs to be
      // deactivated this is the case if ten time steps have passed since the completed onset of
      // scatra-scatra interface layer growth or if the minimum interfacial overpotential is
      // positive and increasing
      if (scatratimint_->Step() == intlayergrowth_startstep_ + 10 or
          (etagrowthmin > 0.0 and etagrowthmin > etagrowthmin_))
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
}  // SCATRA::MeshtyingStrategyS2IElch::EvaluateMeshtying

/*-----------------------------------------------------------------------*
 *-----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2IElch::EvaluatePointCoupling()
{
  // extract multi-scale coupling conditions
  // loop over conditions
  for (const auto& slave_condition : KineticsConditionsMeshtyingSlaveSide())
  {
    auto* cond_slave = slave_condition.second;

    // only evaluate point coupling conditions
    if (cond_slave->GType() != DRT::Condition::GeometryType::Point) continue;

    auto* cond_master = MasterConditions()[slave_condition.first];

    // extract nodal cloud
    const std::vector<int>* const nodeids_slave = cond_slave->Nodes();
    const std::vector<int>* const nodeids_master = cond_master->Nodes();

    if (nodeids_slave->size() != 1 or nodeids_master->size() != 1)
      dserror("only one node per condition allowed");

    const int nodeid_slave = (*nodeids_slave)[0];
    const int nodeid_master = (*nodeids_master)[0];

    auto dis = scatratimint_->Discretization();

    auto* slave_node = dis->gNode(nodeid_slave);
    auto* master_node = dis->gNode(nodeid_master);

    // extract degrees of freedom from node
    const std::vector<int> slave_dofs = dis->Dof(0, slave_node);
    const std::vector<int> master_dofs = dis->Dof(0, master_node);

    const int ed_conc_gid = slave_dofs[0];
    const int ed_pot_gid = slave_dofs[1];
    const int el_conc_gid = master_dofs[0];
    const int el_pot_gid = master_dofs[1];

    auto dof_row_map = scatratimint_->DofRowMap();
    const int ed_conc_lid = dof_row_map->LID(ed_conc_gid);
    const int ed_pot_lid = dof_row_map->LID(ed_pot_gid);
    const int el_conc_lid = dof_row_map->LID(el_conc_gid);
    const int el_pot_lid = dof_row_map->LID(el_pot_gid);

    // extract electrode-side and electrolyte-side values at coupling point
    auto phinp = scatratimint_->Phinp();
    const double ed_conc = (*phinp)[ed_conc_lid];
    const double ed_pot = (*phinp)[ed_pot_lid];
    const double el_conc = (*phinp)[el_conc_lid];
    const double el_pot = (*phinp)[el_pot_lid];

    // compute matrix and vector contributions according to kinetic model for current point coupling
    // condition
    const int kinetic_model = cond_slave->GetInt("kinetic model");
    switch (kinetic_model)
    {
      case INPAR::S2I::kinetics_butlervolmer:
      case INPAR::S2I::kinetics_butlervolmerreduced:
      {
        // access material of electrode
        auto matelectrode =
            Teuchos::rcp_dynamic_cast<const MAT::Electrode>(slave_node->Elements()[0]->Material());
        if (matelectrode == Teuchos::null)
          dserror("Invalid electrode material for multi-scale coupling!");

        // access input parameters associated with current condition
        const int nume = cond_slave->GetInt("e-");
        if (nume != 1)
        {
          dserror(
              "Invalid number of electrons involved in charge transfer at "
              "electrode-electrolyte interface!");
        }
        const std::vector<int>* stoichiometries =
            cond_slave->Get<std::vector<int>>("stoichiometries");
        if (stoichiometries == nullptr)
        {
          dserror(
              "Cannot access vector of stoichiometric coefficients for multi-scale "
              "coupling!");
        }
        if (stoichiometries->size() != 1)
          dserror("Number of stoichiometric coefficients does not match number of scalars!");
        if ((*stoichiometries)[0] != -1) dserror("Invalid stoichiometric coefficient!");
        const double faraday =
            DRT::Problem::Instance(0)->ELCHControlParams().get<double>("FARADAY_CONSTANT");
        const double gasconstant =
            DRT::Problem::Instance(0)->ELCHControlParams().get<double>("GAS_CONSTANT");
        const double frt =
            faraday / (gasconstant * (DRT::Problem::Instance(0)->ELCHControlParams().get<double>(
                                         "TEMPERATURE")));
        const double alphaa = cond_slave->GetDouble("alpha_a");
        const double alphac = cond_slave->GetDouble("alpha_c");
        const double kr = cond_slave->GetDouble("k_r");
        if (kr < 0.0) dserror("Charge transfer constant k_r is negative!");

        // extract saturation value of intercalated lithium concentration from electrode material
        const double cmax = matelectrode->CMax();
        if (cmax < 1.0e-12)
          dserror("Saturation value c_max of intercalated lithium concentration is too small!");

        // compute domain integration factor
        constexpr double four_pi = 4.0 * M_PI;
        const double fac = DRT::INPUT::IntegralValue<bool>(
                               *scatratimint_->ScatraParameterList(), "SPHERICALCOORDS")
                               ? *slave_node->X().data() * *slave_node->X().data() * four_pi
                               : 1.0;
        const double timefacfac =
            DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance(dis->Name())->TimeFac() * fac;
        const double timefacrhsfac =
            DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance(dis->Name())->TimeFacRhs() * fac;
        if (timefacfac < 0.0 or timefacrhsfac < 0.0) dserror("Integration factor is negative!");
        // no deformation available
        const double dummy_detF(1.0);

        // equilibrium electric potential difference and its derivative w.r.t. concentration
        // at electrode surface
        const double epd =
            matelectrode->ComputeOpenCircuitPotential(ed_conc, faraday, frt, dummy_detF);
        const double epdderiv = matelectrode->ComputeDOpenCircuitPotentialDConcentration(
            ed_conc, faraday, frt, dummy_detF);

        // overpotential
        const double eta = ed_pot - el_pot - epd;

        // Butler-Volmer exchange mass flux density
        const double j0 =
            cond_slave->GetInt("kinetic model") == INPAR::S2I::kinetics_butlervolmerreduced
                ? kr
                : kr * std::pow(el_conc, alphaa) * std::pow(cmax - ed_conc, alphaa) *
                      std::pow(ed_conc, alphac);

        // exponential Butler-Volmer terms
        const double expterm1 = std::exp(alphaa * frt * eta);
        const double expterm2 = std::exp(-alphac * frt * eta);
        const double expterm = expterm1 - expterm2;

        // safety check
        if (std::abs(expterm) > 1.0e5)
        {
          dserror(
              "Overflow of exponential term in Butler-Volmer formulation detected! Value: "
              "%lf",
              expterm);
        }

        // core residual term associated with Butler-Volmer mass flux density
        const double j = j0 * expterm;

        // initialize a dummy resistance as the method below requires a resistance which is not
        // relevant in this case
        const double dummyresistance(0.0);
        // define flux linearization terms
        double dj_ded_conc(0.0), dj_del_conc(0.0), dj_ded_pot(0.0), dj_del_pot(0.0);
        // calculate flux linearizations
        DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeUtils::
            CalculateButlerVolmerElchLinearizations(kinetic_model, j0, frt, epdderiv, alphaa,
                alphac, dummyresistance, expterm1, expterm2, kr, faraday, el_conc, ed_conc, cmax,
                eta, dj_ded_conc, dj_del_conc, dj_ded_pot, dj_del_pot);

        // assemble concentration residuals
        auto residual = scatratimint_->Residual();
        (*residual)[ed_conc_lid] -= timefacrhsfac * j;
        (*residual)[el_conc_lid] -= timefacrhsfac * j * -1.0;

        // assemble potential residuals
        (*residual)[ed_pot_lid] -= timefacrhsfac * nume * j;
        (*residual)[el_pot_lid] -= timefacrhsfac * nume * j * -1.0;

        // assemble concentration linearizations
        auto sys_mat = scatratimint_->SystemMatrixOperator();
        sys_mat->Assemble(timefacfac * dj_ded_conc, ed_conc_gid, ed_conc_gid);
        sys_mat->Assemble(timefacfac * dj_del_conc, ed_conc_gid, el_conc_gid);
        sys_mat->Assemble(timefacfac * dj_ded_pot, ed_conc_gid, ed_pot_gid);
        sys_mat->Assemble(timefacfac * dj_del_pot, ed_conc_gid, el_pot_gid);

        sys_mat->Assemble(timefacfac * dj_ded_conc * -1.0, el_conc_gid, ed_conc_gid);
        sys_mat->Assemble(timefacfac * dj_del_conc * -1.0, el_conc_gid, el_conc_gid);
        sys_mat->Assemble(timefacfac * dj_ded_pot * -1.0, el_conc_gid, ed_pot_gid);
        sys_mat->Assemble(timefacfac * dj_del_pot * -1.0, el_conc_gid, el_pot_gid);

        // assemble potential linearizations
        sys_mat->Assemble(timefacfac * nume * dj_ded_conc, ed_pot_gid, ed_conc_gid);
        sys_mat->Assemble(timefacfac * nume * dj_del_conc, ed_pot_gid, el_conc_gid);
        sys_mat->Assemble(timefacfac * nume * dj_ded_pot, ed_pot_gid, ed_pot_gid);
        sys_mat->Assemble(timefacfac * nume * dj_del_pot, ed_pot_gid, el_pot_gid);

        sys_mat->Assemble(timefacfac * nume * dj_ded_conc * -1.0, el_pot_gid, ed_conc_gid);
        sys_mat->Assemble(timefacfac * nume * dj_del_conc * -1.0, el_pot_gid, el_conc_gid);
        sys_mat->Assemble(timefacfac * nume * dj_ded_pot * -1.0, el_pot_gid, ed_pot_gid);
        sys_mat->Assemble(timefacfac * nume * dj_del_pot * -1.0, el_pot_gid, el_pot_gid);

        break;
      }
      case INPAR::S2I::kinetics_nointerfaceflux:
        break;

      default:
      {
        dserror("Kinetic model for s2i coupling not yet implemented!");
        break;
      }
    }
  }
}

/*------------------------------------------------------------------------*
 | instantiate strategy for Newton-Raphson convergence check   fang 02/16 |
 *------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2IElch::InitConvCheckStrategy()
{
  if (couplingtype_ == INPAR::S2I::coupling_mortar_saddlepoint_petrov or
      couplingtype_ == INPAR::S2I::coupling_mortar_saddlepoint_bubnov)
  {
    convcheckstrategy_ = Teuchos::rcp(new SCATRA::ConvCheckStrategyS2ILMElch(
        scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));
  }
  else if (ElchTimInt()->MacroScale())
  {
    convcheckstrategy_ = Teuchos::rcp(new SCATRA::ConvCheckStrategyStdMacroScaleElch(
        scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));
  }
  else
    convcheckstrategy_ = Teuchos::rcp(new SCATRA::ConvCheckStrategyStdElch(
        scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));
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
    for (const auto& condition : conditions)
    {
      // extract current condition
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
          for (int nodegid : *nodegids)
          {
            // extract global ID of current node
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
                const double masterphi = (*imasterphi_on_slave_side_np_)[doflid_scatra];

                // extract master-side electric potential associated with current node
                const double masterpot = (*imasterphi_on_slave_side_np_)[doflid_scatra + 1];

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
                  {
                    dserror(
                        "Local Newton-Raphson iteration for scatra-scatra interface layer growth "
                        "did not converge!");
                  }

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
}


/*----------------------------------------------------------------------*
 | singleton access method                                   fang 01/16 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
SCATRA::MortarCellCalcElch<distypeS, distypeM>*
SCATRA::MortarCellCalcElch<distypeS, distypeM>::Instance(
    const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
    const INPAR::S2I::InterfaceSides&
        lmside,  //!< flag for interface side underlying Lagrange multiplier definition
    const int& numdofpernode_slave,   //!< number of slave-side degrees of freedom per node
    const int& numdofpernode_master,  //!< number of master-side degrees of freedom per node
    const std::string& disname        //!< name of mortar discretization
)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::string>(
      [](const INPAR::S2I::CouplingType& couplingtype, const INPAR::S2I::InterfaceSides& lmside,
          const int& numdofpernode_slave, const int& numdofpernode_master)
      {
        return std::unique_ptr<MortarCellCalcElch<distypeS, distypeM>>(
            new MortarCellCalcElch<distypeS, distypeM>(
                couplingtype, lmside, numdofpernode_slave, numdofpernode_master));
      });

  return singleton_map[disname].Instance(CORE::UTILS::SingletonAction::create, couplingtype, lmside,
      numdofpernode_slave, numdofpernode_master);
}


/*----------------------------------------------------------------------*
 | protected constructor for singletons                      fang 01/16 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
SCATRA::MortarCellCalcElch<distypeS, distypeM>::MortarCellCalcElch(
    const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
    const INPAR::S2I::InterfaceSides&
        lmside,  //!< flag for interface side underlying Lagrange multiplier definition
    const int& numdofpernode_slave,  //!< number of slave-side degrees of freedom per node
    const int& numdofpernode_master  //!< number of master-side degrees of freedom per node
    )
    : my::MortarCellCalc(couplingtype, lmside, numdofpernode_slave, numdofpernode_master)
{
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
void SCATRA::MortarCellCalcElch<distypeS, distypeM>::EvaluateCondition(
    const DRT::Discretization& idiscret, MORTAR::IntCell& cell, MORTAR::MortarElement& slaveelement,
    MORTAR::MortarElement& masterelement, DRT::Element::LocationArray& la_slave,
    DRT::Element::LocationArray& la_master, const Teuchos::ParameterList& params,
    CORE::LINALG::SerialDenseMatrix& k_ss, CORE::LINALG::SerialDenseMatrix& k_sm,
    CORE::LINALG::SerialDenseMatrix& k_ms, CORE::LINALG::SerialDenseMatrix& k_mm,
    CORE::LINALG::SerialDenseVector& r_s, CORE::LINALG::SerialDenseVector& r_m)
{
  // safety checks
  if (my::numdofpernode_slave_ != 2 or my::numdofpernode_master_ != 2)
    dserror("Invalid number of degrees of freedom per node!");
  if (DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->EquPot() !=
      INPAR::ELCH::equpot_divi)
    dserror("Invalid closing equation for electric potential!");

  // extract condition from parameter list
  DRT::Condition* condition = params.get<DRT::Condition*>("condition");
  if (condition == nullptr) dserror("Cannot access scatra-scatra interface coupling condition!");

  // access material of slave element
  Teuchos::RCP<const MAT::Electrode> matelectrode =
      Teuchos::rcp_dynamic_cast<const MAT::Electrode>(slaveelement.Material());
  if (matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // extract nodal state variables associated with slave and master elements
  this->ExtractNodeValues(idiscret, la_slave, la_master);

  // determine quadrature rule
  const CORE::DRT::UTILS::IntPointsAndWeights<2> intpoints(
      CORE::DRT::UTILS::GaussRule2D::tri_7point);

  // dummy matrix of nodal temperature values
  CORE::LINALG::Matrix<nen_slave_, 1> dummy_slave_temp(true);
  CORE::LINALG::Matrix<nen_master_, 1> dummy_master_temp(true);
  // always in contact
  const double pseudo_contact_fac = 1.0;

  // loop over all integration points
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    // evaluate shape functions and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndDomIntFacAtIntPoint(
        slaveelement, masterelement, cell, intpoints, iquad);
    // no deformation available
    const double dummy_detF(1.0);

    // overall integration factors
    const double timefacfac =
        DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->TimeFac() * fac;
    const double timefacrhsfac =
        DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->TimeFacRhs() * fac;
    if (timefacfac < 0.0 or timefacrhsfac < 0.0) dserror("Integration factor is negative!");

    DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<
        distypeS>::template EvaluateS2ICouplingAtIntegrationPoint<distypeM>(matelectrode,
        my::ephinp_slave_, my::ephinp_master_, dummy_slave_temp, dummy_master_temp,
        pseudo_contact_fac, my::funct_slave_, my::funct_master_, my::test_lm_slave_,
        my::test_lm_master_, my::scatraparamsboundary_, timefacfac, timefacrhsfac, dummy_detF,
        GetFRT(), my::numdofpernode_slave_, k_ss, k_sm, k_ms, k_mm, r_s, r_m);
  }
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
void SCATRA::MortarCellCalcElch<distypeS, distypeM>::EvaluateConditionNTS(DRT::Condition& condition,
    const MORTAR::MortarNode& slavenode, const double& lumpedarea,
    MORTAR::MortarElement& slaveelement, MORTAR::MortarElement& masterelement,
    const std::vector<CORE::LINALG::Matrix<nen_slave_, 1>>& ephinp_slave,
    const std::vector<CORE::LINALG::Matrix<nen_master_, 1>>& ephinp_master,
    CORE::LINALG::SerialDenseMatrix& k_ss, CORE::LINALG::SerialDenseMatrix& k_sm,
    CORE::LINALG::SerialDenseMatrix& k_ms, CORE::LINALG::SerialDenseMatrix& k_mm,
    CORE::LINALG::SerialDenseVector& r_s, CORE::LINALG::SerialDenseVector& r_m)
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

  // dummy matrix of nodal temperature values
  CORE::LINALG::Matrix<nen_slave_, 1> dummy_slave_temp(true);
  CORE::LINALG::Matrix<nen_master_, 1> dummy_master_temp(true);
  // always in contact
  const double pseudo_contact_fac = 1.0;

  // overall integration factors
  const double timefacfac =
      DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->TimeFac() * lumpedarea;
  const double timefacrhsfac =
      DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->TimeFacRhs() * lumpedarea;
  if (timefacfac < 0. or timefacrhsfac < 0.) dserror("Integration factor is negative!");

  // no deformation available
  const double dummy_detF(1.0);

  DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<
      distypeS>::template EvaluateS2ICouplingAtIntegrationPoint<distypeM>(matelectrode,
      ephinp_slave, ephinp_master, dummy_slave_temp, dummy_master_temp, pseudo_contact_fac,
      my::funct_slave_, my::funct_master_, my::funct_slave_, my::funct_master_,
      my::scatraparamsboundary_, timefacfac, timefacrhsfac, dummy_detF,
      DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->FRT(), my::numdofpernode_slave_,
      k_ss, k_sm, k_ms, k_mm, r_s, r_m);
}


/*----------------------------------------------------------------------*
 | evaluate factor F/RT                                      fang 01/17 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
double SCATRA::MortarCellCalcElch<distypeS, distypeM>::GetFRT() const
{
  // fetch factor F/RT from electrochemistry parameter list
  return DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->FRT();
}


/*----------------------------------------------------------------------*
 | singleton access method                                   fang 01/17 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
SCATRA::MortarCellCalcElchSTIThermo<distypeS, distypeM>*
SCATRA::MortarCellCalcElchSTIThermo<distypeS, distypeM>::Instance(
    const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
    const INPAR::S2I::InterfaceSides&
        lmside,  //!< flag for interface side underlying Lagrange multiplier definition
    const int& numdofpernode_slave,   //!< number of slave-side degrees of freedom per node
    const int& numdofpernode_master,  //!< number of master-side degrees of freedom per node
    const std::string& disname        //!< name of mortar discretization
)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::string>(
      [](const INPAR::S2I::CouplingType& couplingtype, const INPAR::S2I::InterfaceSides& lmside,
          const int& numdofpernode_slave, const int& numdofpernode_master)
      {
        return std::unique_ptr<MortarCellCalcElchSTIThermo<distypeS, distypeM>>(
            new MortarCellCalcElchSTIThermo<distypeS, distypeM>(
                couplingtype, lmside, numdofpernode_slave, numdofpernode_master));
      });

  return singleton_map[disname].Instance(CORE::UTILS::SingletonAction::create, couplingtype, lmside,
      numdofpernode_slave, numdofpernode_master);
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 01/17 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
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
}


/*--------------------------------------------------------------------------------------------------------------------*
 | evaluate single mortar integration cell of particular slave-side and master-side discretization
 types   fang 01/17 |
 *--------------------------------------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
void SCATRA::MortarCellCalcElchSTIThermo<distypeS, distypeM>::Evaluate(
    const DRT::Discretization& idiscret,           //!< interface discretization
    MORTAR::IntCell& cell,                         //!< mortar integration cell
    MORTAR::MortarElement& slaveelement,           //!< slave-side mortar element
    MORTAR::MortarElement& masterelement,          //!< master-side mortar element
    DRT::Element::LocationArray& la_slave,         //!< slave-side location array
    DRT::Element::LocationArray& la_master,        //!< master-side location array
    const Teuchos::ParameterList& params,          //!< parameter list
    CORE::LINALG::SerialDenseMatrix& cellmatrix1,  //!< cell matrix 1
    CORE::LINALG::SerialDenseMatrix& cellmatrix2,  //!< cell matrix 2
    CORE::LINALG::SerialDenseMatrix& cellmatrix3,  //!< cell matrix 3
    CORE::LINALG::SerialDenseMatrix& cellmatrix4,  //!< cell matrix 4
    CORE::LINALG::SerialDenseVector& cellvector1,  //!< cell vector 1
    CORE::LINALG::SerialDenseVector& cellvector2   //!< cell vector 2
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
}


/*---------------------------------------------------------------------------*
 | evaluate and assemble off-diagonal interface linearizations    fang 01/17 |
 *---------------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
void SCATRA::MortarCellCalcElchSTIThermo<distypeS, distypeM>::EvaluateConditionOD(
    const DRT::Discretization& idiscret,     //!< interface discretization
    MORTAR::IntCell& cell,                   //!< mortar integration cell
    MORTAR::MortarElement& slaveelement,     //!< slave-side mortar element
    MORTAR::MortarElement& masterelement,    //!< master-side mortar element
    DRT::Element::LocationArray& la_slave,   //!< slave-side location array
    DRT::Element::LocationArray& la_master,  //!< master-side location array
    const Teuchos::ParameterList& params,    //!< parameter list
    CORE::LINALG::SerialDenseMatrix&
        k_ss,  //!< linearizations of slave-side residuals w.r.t. slave-side dofs
    CORE::LINALG::SerialDenseMatrix&
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
  if (s2icondition == nullptr) dserror("Cannot access scatra-scatra interface coupling condition!");

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
  const CORE::DRT::UTILS::IntPointsAndWeights<2> intpoints(
      CORE::DRT::UTILS::GaussRule2D::tri_7point);

  // dummy matrix of nodal master temperature values and shape derivatives
  CORE::LINALG::Matrix<nen_master_, 1> dummy_master_temp(true);
  CORE::LINALG::Matrix<nsd_slave_ + 1, nen_slave_> dummy_shapederivatives(true);
  // always in contact
  const double pseudo_contact_fac = 1.0;

  // loop over integration points
  for (int gpid = 0; gpid < intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac =
        my::EvalShapeFuncAndDomIntFacAtIntPoint(slaveelement, masterelement, cell, intpoints, gpid);

    // evaluate overall integration factor
    const double timefac = DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->TimeFac();
    const double timefacfac = timefac * fac;
    if (timefacfac < 0.) dserror("Integration factor is negative!");

    const double timefacwgt = timefac * intpoints.IP().qwgt[gpid];

    // no deformation available
    const double dummy_detF(1.0);

    DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<
        distypeS>::template EvaluateS2ICouplingODAtIntegrationPoint<distypeM>(matelectrode,
        my::ephinp_slave_, etempnp_slave_, dummy_master_temp, my::ephinp_master_,
        pseudo_contact_fac, my::funct_slave_, my::funct_master_, my::test_lm_slave_,
        my::test_lm_master_, dummy_shapederivatives, dummy_shapederivatives,
        my::scatraparamsboundary_, SCATRA::DifferentiationType::temp, timefacfac, timefacwgt,
        dummy_detF, my::numdofpernode_slave_, k_ss, k_ms);
  }  // loop over integration points
}


/*------------------------------------------------------------------------------------*
 | extract nodal state variables associated with mortar integration cell   fang 01/17 |
 *------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
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
}


/*----------------------------------------------------------------------*
 | evaluate factor F/RT                                      fang 01/17 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
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
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
SCATRA::MortarCellCalcSTIElch<distypeS, distypeM>*
SCATRA::MortarCellCalcSTIElch<distypeS, distypeM>::Instance(
    const INPAR::S2I::CouplingType& couplingtype,  //!< flag for meshtying method
    const INPAR::S2I::InterfaceSides&
        lmside,  //!< flag for interface side underlying Lagrange multiplier definition
    const int& numdofpernode_slave,   //!< number of slave-side degrees of freedom per node
    const int& numdofpernode_master,  //!< number of master-side degrees of freedom per node
    const std::string& disname        //!< name of mortar discretization
)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::string>(
      [](const INPAR::S2I::CouplingType& couplingtype, const INPAR::S2I::InterfaceSides& lmside,
          const int& numdofpernode_slave, const int& numdofpernode_master)
      {
        return std::unique_ptr<MortarCellCalcSTIElch<distypeS, distypeM>>(
            new MortarCellCalcSTIElch<distypeS, distypeM>(
                couplingtype, lmside, numdofpernode_slave, numdofpernode_master));
      });

  return singleton_map[disname].Instance(CORE::UTILS::SingletonAction::create, couplingtype, lmside,
      numdofpernode_slave, numdofpernode_master);
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 01/17 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
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
      eelchnp_slave_(2, CORE::LINALG::Matrix<nen_slave_, 1>(true)),
      eelchnp_master_(2, CORE::LINALG::Matrix<nen_master_, 1>(true))
{
}


/*--------------------------------------------------------------------------------------------------------------------*
 | evaluate single mortar integration cell of particular slave-side and master-side discretization
 types   fang 01/17 |
 *--------------------------------------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
void SCATRA::MortarCellCalcSTIElch<distypeS, distypeM>::Evaluate(
    const DRT::Discretization& idiscret,           //!< interface discretization
    MORTAR::IntCell& cell,                         //!< mortar integration cell
    MORTAR::MortarElement& slaveelement,           //!< slave-side mortar element
    MORTAR::MortarElement& masterelement,          //!< master-side mortar element
    DRT::Element::LocationArray& la_slave,         //!< slave-side location array
    DRT::Element::LocationArray& la_master,        //!< master-side location array
    const Teuchos::ParameterList& params,          //!< parameter list
    CORE::LINALG::SerialDenseMatrix& cellmatrix1,  //!< cell matrix 1
    CORE::LINALG::SerialDenseMatrix& cellmatrix2,  //!< cell matrix 2
    CORE::LINALG::SerialDenseMatrix& cellmatrix3,  //!< cell matrix 3
    CORE::LINALG::SerialDenseMatrix& cellmatrix4,  //!< cell matrix 4
    CORE::LINALG::SerialDenseVector& cellvector1,  //!< cell vector 1
    CORE::LINALG::SerialDenseVector& cellvector2   //!< cell vector 2
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
}


/*---------------------------------------------------------------------------*
 | evaluate and assemble interface linearizations and residuals   fang 01/17 |
 *---------------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
void SCATRA::MortarCellCalcSTIElch<distypeS, distypeM>::EvaluateCondition(
    const DRT::Discretization& idiscret,     //!< interface discretization
    MORTAR::IntCell& cell,                   //!< mortar integration cell
    MORTAR::MortarElement& slaveelement,     //!< slave-side mortar element
    MORTAR::MortarElement& masterelement,    //!< master-side mortar element
    DRT::Element::LocationArray& la_slave,   //!< slave-side location array
    DRT::Element::LocationArray& la_master,  //!< master-side location array
    const Teuchos::ParameterList& params,    //!< parameter list
    CORE::LINALG::SerialDenseMatrix&
        k_ss,  //!< linearizations of slave-side residuals w.r.t. slave-side dofs
    CORE::LINALG::SerialDenseVector& r_s  //!< slave-side residual vector
)
{
  // safety check
  if (my::numdofpernode_slave_ != 1 or my::numdofpernode_master_ != 1)
    dserror("Invalid number of degrees of freedom per node!");

  // extract condition from parameter list
  DRT::Condition* s2icondition = params.get<DRT::Condition*>("condition");
  if (s2icondition == nullptr) dserror("Cannot access scatra-scatra interface coupling condition!");

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
  const CORE::DRT::UTILS::IntPointsAndWeights<2> intpoints(
      CORE::DRT::UTILS::GaussRule2D::tri_7point);

  // dummy matrix for derivative of slave fluxes w.r.t. master side temperatures
  CORE::LINALG::SerialDenseMatrix dummy_ksm;
  // always in contact
  const double pseudo_contact_fac = 1.0;

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

    // no deformation available
    const double dummy_detF(1.0);

    DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<
        distypeS>::template EvaluateS2ICouplingAtIntegrationPoint<distypeM>(matelectrode,
        my::ephinp_slave_[0], my::ephinp_master_[0], eelchnp_slave_, eelchnp_master_,
        pseudo_contact_fac, my::funct_slave_, my::funct_master_, my::scatraparamsboundary_,
        timefacfac, timefacrhsfac, dummy_detF, k_ss, dummy_ksm, r_s);
  }  // loop over integration points
}


/*---------------------------------------------------------------------------*
 | evaluate and assemble off-diagonal interface linearizations    fang 01/17 |
 *---------------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
void SCATRA::MortarCellCalcSTIElch<distypeS, distypeM>::EvaluateConditionOD(
    const DRT::Discretization& idiscret,     //!< interface discretization
    MORTAR::IntCell& cell,                   //!< mortar integration cell
    MORTAR::MortarElement& slaveelement,     //!< slave-side mortar element
    MORTAR::MortarElement& masterelement,    //!< master-side mortar element
    DRT::Element::LocationArray& la_slave,   //!< slave-side location array
    DRT::Element::LocationArray& la_master,  //!< master-side location array
    const Teuchos::ParameterList& params,    //!< parameter list
    CORE::LINALG::SerialDenseMatrix&
        k_ss,  //!< linearizations of slave-side residuals w.r.t. slave-side dofs
    CORE::LINALG::SerialDenseMatrix&
        k_sm  //!< linearizations of slave-side residuals w.r.t. master-side dofs
)
{
  // safety check
  if (my::numdofpernode_slave_ != 1 or my::numdofpernode_master_ != 1)
    dserror("Invalid number of degrees of freedom per node!");

  // extract condition from parameter list
  DRT::Condition* s2icondition = params.get<DRT::Condition*>("condition");
  if (s2icondition == nullptr) dserror("Cannot access scatra-scatra interface coupling condition!");

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
  const CORE::DRT::UTILS::IntPointsAndWeights<2> intpoints(
      CORE::DRT::UTILS::GaussRule2D::tri_7point);

  // dummy matrix for shape derivatives
  CORE::LINALG::Matrix<3, nen_slave_> dummy_shape_deriv;
  // always in contact
  const double pseudo_contact_fac = 1.0;

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

    // no deformation available
    const double dummy_detF(1.0);

    DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<
        distypeS>::template EvaluateS2ICouplingODAtIntegrationPoint<distypeM>(matelectrode,
        my::ephinp_slave_[0], my::ephinp_master_[0], eelchnp_slave_, eelchnp_master_,
        pseudo_contact_fac, my::funct_slave_, my::funct_master_, my::scatraparamsboundary_,
        timefacfac, fac, dummy_detF, SCATRA::DifferentiationType::elch, dummy_shape_deriv,
        dummy_shape_deriv, k_ss, k_sm);
  }  // loop over integration points
}


/*------------------------------------------------------------------------------------*
 | extract nodal state variables associated with mortar integration cell   fang 01/17 |
 *------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
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
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SCATRA::MeshtyingStrategyS2IElchSCL::MeshtyingStrategyS2IElchSCL(
    SCATRA::ScaTraTimIntElch* elchtimint, const Teuchos::ParameterList& parameters)
    : MeshtyingStrategyS2IElch(elchtimint, parameters)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2IElchSCL::SetupMeshtying()
{
  // extract scatra-scatra coupling conditions from discretization
  std::vector<DRT::Condition*> s2imeshtying_conditions(0, nullptr);
  scatratimint_->Discretization()->GetCondition("S2IMeshtying", s2imeshtying_conditions);

  std::set<int> islavenodegidset;
  std::set<int> imasternodegidset;

  for (const auto& s2imeshtying_condition : s2imeshtying_conditions)
  {
    if (s2imeshtying_condition->GetInt("S2IKineticsID") != -1)
      dserror("No kinetics condition is allowed for the coupled space-charge layer problem.");

    switch (s2imeshtying_condition->GetInt("interface side"))
    {
      case INPAR::S2I::side_slave:
      {
        DRT::UTILS::AddOwnedNodeGIDFromList(
            *scatratimint_->Discretization(), *s2imeshtying_condition->Nodes(), islavenodegidset);
        break;
      }
      case INPAR::S2I::side_master:
      {
        DRT::UTILS::AddOwnedNodeGIDFromList(
            *scatratimint_->Discretization(), *s2imeshtying_condition->Nodes(), imasternodegidset);
        break;
      }
      default:
      {
        dserror("interface side must bee slave or master");
        break;
      }
    }
  }

  std::vector<int> islavenodegidvec(islavenodegidset.begin(), islavenodegidset.end());
  std::vector<int> imasternodegidvec(imasternodegidset.begin(), imasternodegidset.end());

  icoup_ = Teuchos::rcp(new CORE::ADAPTER::Coupling());
  icoup_->SetupCoupling(*(scatratimint_->Discretization()), *(scatratimint_->Discretization()),
      imasternodegidvec, islavenodegidvec, 2, true, 1.0e-8);
}

/*------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2IElchSCL::Solve(const Teuchos::RCP<CORE::LINALG::Solver>& solver,
    const Teuchos::RCP<CORE::LINALG::SparseOperator>& systemmatrix,
    const Teuchos::RCP<Epetra_Vector>& increment, const Teuchos::RCP<Epetra_Vector>& residual,
    const Teuchos::RCP<Epetra_Vector>& phinp, const int& iteration,
    const Teuchos::RCP<CORE::LINALG::KrylovProjector>& projector) const
{
  solver->Solve(
      systemmatrix->EpetraOperator(), increment, residual, true, iteration == 1, projector);
}



// forward declarations
template class SCATRA::MortarCellCalcElch<CORE::FE::CellType::tri3, CORE::FE::CellType::tri3>;
template class SCATRA::MortarCellCalcElch<CORE::FE::CellType::tri3, CORE::FE::CellType::quad4>;
template class SCATRA::MortarCellCalcElch<CORE::FE::CellType::quad4, CORE::FE::CellType::tri3>;
template class SCATRA::MortarCellCalcElch<CORE::FE::CellType::quad4, CORE::FE::CellType::quad4>;
template class SCATRA::MortarCellCalcElchSTIThermo<CORE::FE::CellType::tri3,
    CORE::FE::CellType::tri3>;
template class SCATRA::MortarCellCalcElchSTIThermo<CORE::FE::CellType::tri3,
    CORE::FE::CellType::quad4>;
template class SCATRA::MortarCellCalcElchSTIThermo<CORE::FE::CellType::quad4,
    CORE::FE::CellType::tri3>;
template class SCATRA::MortarCellCalcElchSTIThermo<CORE::FE::CellType::quad4,
    CORE::FE::CellType::quad4>;
template class SCATRA::MortarCellCalcSTIElch<CORE::FE::CellType::tri3, CORE::FE::CellType::tri3>;
template class SCATRA::MortarCellCalcSTIElch<CORE::FE::CellType::tri3, CORE::FE::CellType::quad4>;
template class SCATRA::MortarCellCalcSTIElch<CORE::FE::CellType::quad4, CORE::FE::CellType::tri3>;
template class SCATRA::MortarCellCalcSTIElch<CORE::FE::CellType::quad4, CORE::FE::CellType::quad4>;
