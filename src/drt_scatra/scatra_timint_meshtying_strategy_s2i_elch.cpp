/*!----------------------------------------------------------------------
\file scatra_timint_meshtying_strategy_s2i_elch.cpp

\brief Scatra-scatra interface coupling strategy for electrochemistry problems

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/
#include "scatra_timint_meshtying_strategy_s2i_elch.H"

#include "../drt_lib/drt_discret.H"

#include "../drt_mat/electrode.H"

#include "../drt_mortar/mortar_element.H"

#include "../drt_scatra_ele/scatra_ele_boundary_calc_elch_electrode.H"
#include "../drt_scatra_ele/scatra_ele_calc_utils.H"
#include "../drt_scatra_ele/scatra_ele_parameter_elch.H"
#include "../drt_scatra_ele/scatra_ele_parameter_timint.H"

#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_solver.H"

/*----------------------------------------------------------------------*
 | constructor                                               fang 12/14 |
 *----------------------------------------------------------------------*/
SCATRA::MeshtyingStrategyS2IElch::MeshtyingStrategyS2IElch(
    SCATRA::ScaTraTimIntElch*       elchtimint,   //!< elch time integrator
    const Teuchos::ParameterList&   parameters    //!< input parameters for scatra-scatra interface coupling
    ) :
MeshtyingStrategyS2I(elchtimint,parameters)
{
  return;
} // SCATRA::MeshtyingStrategyS2IElch::MeshtyingStrategyS2IElch


/*--------------------------------------------------------------------------------------*
 | evaluate scatra-scatra interface coupling conditions (electrochemistry)   fang 10/14 |
 *--------------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2IElch::EvaluateMeshtying()
{
  // safety check
  if(DRT::INPUT::IntegralValue<int>(*(ElchTimInt()->ElchParameterList()),"BLOCKPRECOND"))
    dserror("Block preconditioning doesn't work for scatra-scatra interface coupling yet!");

  // call base class routine
  SCATRA::MeshtyingStrategyS2I::EvaluateMeshtying();

  return;
} // SCATRA::MeshtyingStrategyS2IElch::EvaluateMeshtying


/*----------------------------------------------------------------------------*
 | build maps associated with blocks of global system matrix       fang 06/15 |
 *----------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2IElch::BuildBlockMaps(
    const std::vector<Teuchos::RCP<DRT::Condition> >&   partitioningconditions,   //!< domain partitioning conditions
    std::vector<Teuchos::RCP<const Epetra_Map> >&       blockmaps                 //!< empty vector for maps to be built
    ) const
{
  if(matrixtype_ == INPAR::S2I::matrix_block_condition_dof)
  {
    // safety check
    if(DRT::INPUT::IntegralValue<int>(ElchTimInt()->ElchParameterList()->sublist("DIFFCOND"),"CURRENT_SOLUTION_VAR"))
      dserror("For chosen type of global block system matrix, current must not constitute solution variable!");

    // extract number of domain partitioning conditions
    const unsigned ncond = partitioningconditions.size();

    // prepare vector for maps to be built
    blockmaps.resize(ncond*2,Teuchos::null);

    // loop over all domain partitioning conditions
    for(unsigned icond=0; icond<ncond; ++icond)
    {
      // initialize sets for dof IDs associated with current partitioning condition
      std::vector<std::set<int> > dofids(2);

      // extract nodes associated with current domain partitioning condition
      const std::vector<int>* nodegids = partitioningconditions[icond]->Nodes();

      // loop over all nodes associated with current domain partitioning condition
      for (unsigned inode=0; inode<nodegids->size(); ++inode)
      {
        // extract global ID of current node
        const int nodegid = (*nodegids)[inode];

        // consider current node only if node is owned by current processor
        // need to make sure that node is stored on current processor, otherwise cannot resolve "->Owner()"
        if(scatratimint_->Discretization()->HaveGlobalNode(nodegid) and scatratimint_->Discretization()->gNode(nodegid)->Owner() == scatratimint_->Discretization()->Comm().MyPID())
        {
          // extract dof IDs associated with current node
          const std::vector<int> nodedofs = scatratimint_->Discretization()->Dof(scatratimint_->Discretization()->gNode(nodegid));

          // add concentration dof IDs to associated set
          std::copy(nodedofs.begin(),--nodedofs.end(),std::inserter(dofids[0],dofids[0].end()));

          // add electric potential dof ID to associated set
          dofids[1].insert(nodedofs.back());
        }
      }

      // transform sets for dof IDs into vectors and then into Epetra maps
      for(unsigned iset=0; iset<2; ++iset)
      {
        int nummyelements(0);
        int* myglobalelements(NULL);
        std::vector<int> dofidvec;
        if(dofids[iset].size() > 0)
        {
          dofidvec.reserve(dofids[iset].size());
          dofidvec.assign(dofids[iset].begin(),dofids[iset].end());
          nummyelements = dofidvec.size();
          myglobalelements = &(dofidvec[0]);
        }
        blockmaps[2*icond+iset] = Teuchos::rcp(new Epetra_Map(-1,nummyelements,myglobalelements,scatratimint_->DofRowMap()->IndexBase(),scatratimint_->DofRowMap()->Comm()));
      }
    }
  }

  // call base class routine for other types of global system matrix
  else
    SCATRA::MeshtyingStrategyS2I::BuildBlockMaps(partitioningconditions,blockmaps);

  return;
} // SCATRA::MeshtyingStrategyS2I::BuildBlockMaps


/*-------------------------------------------------------------------------------*
 | build null spaces associated with blocks of global system matrix   fang 07/15 |
 *-------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2IElch::BuildBlockNullSpaces() const
{
  // call base class routine
  SCATRA::MeshtyingStrategyS2I::BuildBlockNullSpaces();

  if(matrixtype_ == INPAR::S2I::matrix_block_condition_dof)
  {
    // loop over blocks of global system matrix
    for(int iblock=0; iblock<blockmaps_->NumMaps(); ++iblock)
    {
      // store number of current block as string, starting from 1
      std::ostringstream iblockstr;
      iblockstr << iblock+1;

      // access parameter sublist associated with smoother for current matrix block
      Teuchos::ParameterList& mueluparams = scatratimint_->Solver()->Params().sublist("Inverse"+iblockstr.str()).sublist("MueLu Parameters");

      // extract already reduced null space associated with current matrix block
      std::vector<double>& nullspace = *mueluparams.get<Teuchos::RCP<std::vector<double> > >("nullspace");

      // Each matrix block is associated with either concentration dofs or electric potential dofs only. However, since the original
      // full null space was computed for all degrees of freedom on the discretization, the reduced null spaces still have the full
      // dimension, i.e., the full number of null space vectors equaling the total number of primary variables. Hence, we need to
      // decrease the dimension of each null space by one and remove the corresponding zero null space vector from the null space.
      if(iblock%2 == 0)
        // null space associated with concentration dofs
        // remove zero null space vector associated with electric potential dofs by truncating null space
        nullspace.resize(blockmaps_->Map(iblock)->NumMyElements());

      else
        // null space associated with electric potential dofs
        // remove zero null space vector(s) associated with concentration dofs and retain only the last null space vector associated with electric potential dofs
        nullspace.erase(nullspace.begin(),nullspace.end()-blockmaps_->Map(iblock)->NumMyElements());

      // decrease null space dimension and number of partial differential equations by one
      --mueluparams.get<int>("null space: dimension");
      --mueluparams.get<int>("PDE equations");
    }
  }

  return;
} // SCATRA::MeshtyingStrategyS2IElch::BuildBlockNullSpaces


/*------------------------------------------------------------------------*
 | instantiate strategy for Newton-Raphson convergence check   fang 02/16 |
 *------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2IElch::InitConvCheckStrategy()
{
  if(couplingtype_ == INPAR::S2I::coupling_mortar_saddlepoint_petrov or couplingtype_ == INPAR::S2I::coupling_mortar_saddlepoint_bubnov)
    convcheckstrategy_ = Teuchos::rcp(new SCATRA::ConvCheckStrategyS2ILMElch(scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));
  else
    convcheckstrategy_ = Teuchos::rcp(new SCATRA::ConvCheckStrategyStdElch(scatratimint_->ScatraParameterList()->sublist("NONLINEAR")));

  return;
} // SCATRA::MeshtyingStrategyS2IElch::InitConvCheckStrategy


/*------------------------------------------------------------------------------------------*
 | update solution after convergence of the nonlinear Newton-Raphson iteration   fang 01/17 |
 *------------------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2IElch::Update() const
{
  // update scatra-scatra interface layer thicknesses in case of semi-implicit solution approach
  if(intlayergrowth_evaluation_ == INPAR::S2I::growth_evaluation_semi_implicit)
  {
    // extract boundary conditions for scatra-scatra interface layer growth
    std::vector<DRT::Condition*> conditions;
    scatratimint_->Discretization()->GetCondition("S2ICouplingGrowth",conditions);

    // loop over all conditions
    for(unsigned icond=0; icond<conditions.size(); ++icond)
    {
      // extract current condition
      const DRT::Condition* const condition = conditions[icond];

      // extract kinetic model from current condition
      switch(condition->GetInt("kinetic model"))
      {
        case INPAR::S2I::growth_kinetics_butlervolmer:
        {
          // extract parameters from current condition
          const double kr = condition->GetDouble("k_r");
          const double alphaa = condition->GetDouble("alpha_a");
          const double alphac = condition->GetDouble("alpha_c");
          const double frt = ElchTimInt()->FRT();
          const double conductivity_inverse = 1./condition->GetDouble("conductivity");

          // pre-compute integration factor
          const double integrationfac(condition->GetDouble("molar mass")*scatratimint_->Dt()/(condition->GetDouble("density")*INPAR::ELCH::faraday_const));

          // extract nodal cloud from current condition
          const std::vector<int>* nodegids = condition->Nodes();

          // loop over all nodes
          for(unsigned inode=0; inode<nodegids->size(); ++inode)
          {
            // extract global ID of current node
            const int nodegid((*nodegids)[inode]);

            // extract current node
            const DRT::Node* const node = scatratimint_->Discretization()->gNode(nodegid);

            // process only nodes owned by current processor
            if(scatratimint_->Discretization()->HaveGlobalNode(nodegid) and node->Owner() == scatratimint_->Discretization()->Comm().MyPID())
            {
              // extract local ID of first scalar transport degree of freedom associated with current node
              const int doflid_scatra = scatratimint_->Discretization()->DofRowMap()->LID(scatratimint_->Discretization()->Dof(node,0));
              if(doflid_scatra < 0)
                dserror("Couldn't extract local ID of scalar transport degree of freedom!");

              // extract local ID of scatra-scatra interface layer thickness variable associated with current node
              const int doflid_growth = scatratimint_->Discretization()->DofRowMap(2)->LID(scatratimint_->Discretization()->Dof(2,node,0));
              if(doflid_growth < 0)
                dserror("Couldn't extract local ID of scatra-scatra interface layer thickness!");

              // extract slave-side electric potential associated with current node
              const double slavepot = (*scatratimint_->Phiafnp())[doflid_scatra+1];

              // extract master-side lithium concentration associated with current node
              const double masterphi = (*imasterphinp_)[doflid_scatra];

              // extract master-side electric potential associated with current node
              const double masterpot = (*imasterphinp_)[doflid_scatra+1];

              // compute interface layer resistance associated with current node
              const double resistance = (*growthn_)[doflid_growth]*conductivity_inverse;

              // check existence of interface layer and set Heaviside value accordingly
              const unsigned heaviside(resistance > 0. ? 1 : 0);

              // compute exchange current density
              const double i0 = kr*INPAR::ELCH::faraday_const*pow(masterphi,alphaa);

              // compute initial guess of Butler-Volmer current density associated with lithium plating, neglecting overpotential due to resistance of plated lithium
              double eta = slavepot-masterpot;
              double i = i0*(heaviside*exp(alphaa*frt*eta)-exp(-alphac*frt*eta));

              // initialize Newton-Raphson iteration counter
              unsigned iternum(0);

              // apply Newton-Raphson method to compute Butler-Volmer current density associated with lithium plating, involving overpotential due to resistance of plated lithium
              while(true)
              {
                // increment counter
                ++iternum;

                // compute current Newton-Raphson residual
                eta = slavepot-masterpot-resistance*i;   // open-circuit potential is zero for lithium plating reaction
                const double expterm1 = heaviside*exp(alphaa*frt*eta);
                const double expterm2 = exp(-alphac*frt*eta);
                const double residual = i0*(expterm1-expterm2)-i;

                // convergence check
                if(std::abs(residual) < intlayergrowth_convtol_)
                  break;
                else if(iternum == intlayergrowth_itemax_)
                  dserror("Local Newton-Raphson iteration for scatra-scatra interface layer growth did not converge!");

                // compute linearization of current Newton-Raphson residual w.r.t. Butler-Volmer current density associated with lithium plating
                const double linearization = -i0*resistance*frt*(alphaa*expterm1+alphac*expterm2)-1.;

                // update Butler-Volmer current density
                i -= residual/linearization;
              }

              // enforce plating condition, i.e., consider initial lithium plating only in case of negative overpotential
              if(!heaviside and eta >= 0.)
                i = 0.;

              // update lithium plating variable
              (*growthn_)[doflid_growth] -= i*integrationfac;
            }
          } // loop over all nodes

          break;
        }

        default:
        {
          dserror("Kinetic model for scatra-scatra interface layer growth is not yet implemented!");
          break;
        }
      } // kinetic models
    } // loop over all conditions
  } // semi-implicit evaluation of scatra-scatra interface layer growth

  else
    // call base class routine
    MeshtyingStrategyS2I::Update();

  return;
}


/*----------------------------------------------------------------------*
 | singleton access method                                   fang 01/16 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,DRT::Element::DiscretizationType distypeM>
SCATRA::MortarCellCalcElch<distypeS,distypeM>* SCATRA::MortarCellCalcElch<distypeS,distypeM>::Instance(
    const INPAR::S2I::CouplingType&     couplingtype,           //!< flag for meshtying method
    const INPAR::S2I::InterfaceSides&   lmside,                 //!< flag for interface side underlying Lagrange multiplier definition
    const int&                          numdofpernode_slave,    //!< number of slave-side degrees of freedom per node
    const int&                          numdofpernode_master,   //!< number of master-side degrees of freedom per node
    bool                                create                  //!< creation flag
    )
{
  static MortarCellCalcElch<distypeS,distypeM>* instance;

  if(create)
  {
    if(instance == NULL)
      instance = new MortarCellCalcElch<distypeS,distypeM>(couplingtype,lmside,numdofpernode_slave,numdofpernode_master);
  }

  else if(instance != NULL)
  {
    delete instance;
    instance = NULL;
  }

  return instance;
}


/*----------------------------------------------------------------------*
 | singleton destruction                                     fang 01/16 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalcElch<distypeS,distypeM>::Done()
{
  // delete singleton
  Instance(INPAR::S2I::coupling_undefined,INPAR::S2I::side_undefined,0,0,false);

  return;
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 01/16 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,DRT::Element::DiscretizationType distypeM>
SCATRA::MortarCellCalcElch<distypeS,distypeM>::MortarCellCalcElch(
    const INPAR::S2I::CouplingType&     couplingtype,          //!< flag for meshtying method
    const INPAR::S2I::InterfaceSides&   lmside,                //!< flag for interface side underlying Lagrange multiplier definition
    const int&                          numdofpernode_slave,   //!< number of slave-side degrees of freedom per node
    const int&                          numdofpernode_master   //!< number of master-side degrees of freedom per node
    ) :
    my::MortarCellCalc(couplingtype,lmside,numdofpernode_slave,numdofpernode_master)
{
  return;
}


/*---------------------------------------------------------------------------*
 | evaluate and assemble interface linearizations and residuals   fang 01/16 |
 *---------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalcElch<distypeS,distypeM>::EvaluateCondition(
    DRT::Condition&                                          condition,       //!< scatra-scatra interface coupling condition
    MORTAR::IntCell&                                         cell,            //!< mortar integration cell
    MORTAR::MortarElement&                                   slaveelement,    //!< slave-side mortar element
    MORTAR::MortarElement&                                   masterelement,   //!< master-side mortar element
    const std::vector<LINALG::Matrix<my::nen_slave_,1> >&    ephinp_slave,    //!< state variables at slave-side nodes
    const std::vector<LINALG::Matrix<my::nen_master_,1> >&   ephinp_master,   //!< state variables at master-side nodes
    Epetra_SerialDenseMatrix&                                k_ss,            //!< linearizations of slave-side residuals w.r.t. slave-side dofs
    Epetra_SerialDenseMatrix&                                k_sm,            //!< linearizations of slave-side residuals w.r.t. master-side dofs
    Epetra_SerialDenseMatrix&                                k_ms,            //!< linearizations of master-side residuals w.r.t. slave-side dofs
    Epetra_SerialDenseMatrix&                                k_mm,            //!< linearizations of master-side residuals w.r.t. master-side dofs
    Epetra_SerialDenseVector&                                r_s,             //!< slave-side residual vector
    Epetra_SerialDenseVector&                                r_m              //!< master-side residual vector
    )
{
  // safety checks
  if(my::numdofpernode_slave_ != 2 or my::numdofpernode_master_ != 2)
    dserror("Invalid number of degrees of freedom per node!");
  if(DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->EquPot() != INPAR::ELCH::equpot_divi)
    dserror("Invalid closing equation for electric potential!");

  // access material of slave element
  Teuchos::RCP<const MAT::Electrode> matelectrode = Teuchos::rcp_dynamic_cast<const MAT::Electrode>(Teuchos::rcp_dynamic_cast<DRT::FaceElement>(condition.Geometry()[slaveelement.Id()])->ParentElement()->Material());
  if(matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // determine quadrature rule
  const DRT::UTILS::IntPointsAndWeights<2> intpoints(DRT::UTILS::intrule_tri_7point);

  // loop over all integration points
  for(int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    // evaluate shape functions and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndDomIntFacAtIntPoint(slaveelement,masterelement,cell,intpoints,iquad);

    // overall integration factors
    const double timefacfac = DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->TimeFac()*fac;
    const double timefacrhsfac = DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->TimeFacRhs()*fac;
    if(timefacfac < 0. or timefacrhsfac < 0.)
      dserror("Integration factor is negative!");

    DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distypeS>::template EvaluateS2ICouplingAtIntegrationPoint<distypeM>(
        condition,
        matelectrode,
        ephinp_slave,
        ephinp_master,
        my::funct_slave_,
        my::funct_master_,
        my::test_lm_slave_,
        my::test_lm_master_,
        timefacfac,
        timefacrhsfac,
        DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->FRT(),
        k_ss,
        k_sm,
        k_ms,
        k_mm,
        r_s,
        r_m
        );
  }

  return;
}


/*--------------------------------------------------------------------------------------------------------*
 | evaluate and assemble interface linearizations and residuals for node-to-segment coupling   fang 08/16 |
 *--------------------------------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalcElch<distypeS,distypeM>::EvaluateConditionNTS(
    DRT::Condition&                                          condition,       //!< scatra-scatra interface coupling condition
    const MORTAR::MortarNode&                                slavenode,       //!< slave-side node
    const double&                                            lumpedarea,      //!< lumped interface area fraction associated with slave-side node
    MORTAR::MortarElement&                                   slaveelement,    //!< slave-side mortar element
    MORTAR::MortarElement&                                   masterelement,   //!< master-side mortar element
    const std::vector<LINALG::Matrix<my::nen_slave_,1> >&    ephinp_slave,    //!< state variables at slave-side nodes
    const std::vector<LINALG::Matrix<my::nen_master_,1> >&   ephinp_master,   //!< state variables at master-side nodes
    Epetra_SerialDenseMatrix&                                k_ss,            //!< linearizations of slave-side residuals w.r.t. slave-side dofs
    Epetra_SerialDenseMatrix&                                k_sm,            //!< linearizations of slave-side residuals w.r.t. master-side dofs
    Epetra_SerialDenseMatrix&                                k_ms,            //!< linearizations of master-side residuals w.r.t. slave-side dofs
    Epetra_SerialDenseMatrix&                                k_mm,            //!< linearizations of master-side residuals w.r.t. master-side dofs
    Epetra_SerialDenseVector&                                r_s,             //!< slave-side residual vector
    Epetra_SerialDenseVector&                                r_m              //!< master-side residual vector
    )
{
  // safety checks
  if(my::numdofpernode_slave_ != 2 or my::numdofpernode_master_ != 2)
    dserror("Invalid number of degrees of freedom per node!");
  if(DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->EquPot() != INPAR::ELCH::equpot_divi)
    dserror("Invalid closing equation for electric potential!");

  // access material of slave element
  Teuchos::RCP<const MAT::Electrode> matelectrode = Teuchos::rcp_dynamic_cast<const MAT::Electrode>(Teuchos::rcp_dynamic_cast<DRT::FaceElement>(condition.Geometry()[slaveelement.Id()])->ParentElement()->Material());
  if(matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // evaluate shape functions at position of slave-side node
  my::EvalShapeFuncAtSlaveNode(slavenode,slaveelement,masterelement);

  // overall integration factors
  const double timefacfac = DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->TimeFac()*lumpedarea;
  const double timefacrhsfac = DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->TimeFacRhs()*lumpedarea;
  if(timefacfac < 0. or timefacrhsfac < 0.)
    dserror("Integration factor is negative!");

  DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distypeS>::template EvaluateS2ICouplingAtIntegrationPoint<distypeM>(
      condition,
      matelectrode,
      ephinp_slave,
      ephinp_master,
      my::funct_slave_,
      my::funct_master_,
      my::funct_slave_,
      my::funct_master_,
      timefacfac,
      timefacrhsfac,
      DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->FRT(),
      k_ss,
      k_sm,
      k_ms,
      k_mm,
      r_s,
      r_m
      );

  return;
}


// forward declarations
template class SCATRA::MortarCellCalcElch<DRT::Element::tri3,DRT::Element::tri3>;
template class SCATRA::MortarCellCalcElch<DRT::Element::tri3,DRT::Element::quad4>;
template class SCATRA::MortarCellCalcElch<DRT::Element::quad4,DRT::Element::tri3>;
template class SCATRA::MortarCellCalcElch<DRT::Element::quad4,DRT::Element::quad4>;
