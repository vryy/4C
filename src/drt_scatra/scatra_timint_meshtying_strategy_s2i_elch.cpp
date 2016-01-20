/*!----------------------------------------------------------------------
\file scatra_timint_meshtying_strategy_s2i_elch.cpp

\brief Scatra-scatra interface coupling strategy for electrochemistry problems

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/

#include "scatra_timint_meshtying_strategy_s2i_elch.H"

#include "../drt_lib/drt_discret.H"

#include "../drt_mat/electrode.H"

#include "../drt_mortar/mortar_element.H"

#include "../drt_scatra_ele/scatra_ele_calc_utils.H"
#include "../drt_scatra_ele/scatra_ele_parameter_elch.H"
#include "../drt_scatra_ele/scatra_ele_parameter_timint.H"

#include "../linalg/linalg_solver.H"

/*----------------------------------------------------------------------*
 | constructor                                               fang 12/14 |
 *----------------------------------------------------------------------*/
SCATRA::MeshtyingStrategyS2IElch::MeshtyingStrategyS2IElch(
    SCATRA::ScaTraTimIntElch*       elchtimint,   //! elch time integrator
    const Teuchos::ParameterList&   parameters    //! input parameters for scatra-scatra interface coupling
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

  SCATRA::MeshtyingStrategyS2I::EvaluateMeshtying();

  return;
} // SCATRA::MeshtyingStrategyS2IElch::EvaluateMeshtying


/*----------------------------------------------------------------------------*
 | build maps associated with blocks of global system matrix       fang 06/15 |
 *----------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2IElch::BuildBlockMaps(
    const std::vector<Teuchos::RCP<DRT::Condition> >&   partitioningconditions,   //! domain partitioning conditions
    std::vector<Teuchos::RCP<const Epetra_Map> >&       blockmaps                 //! empty vector for maps to be built
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


/*----------------------------------------------------------------------------------------*
 | evaluate single mortar integration cell                                     fang 01/16 |
 *----------------------------------------------------------------------------------------*/
void SCATRA::MeshtyingStrategyS2IElch::MortarCellCalc(
    MORTAR::MortarElement&         slaveelement,    //!< slave-side mortar element
    MORTAR::MortarElement&         masterelement,   //!< master-side mortar element
    MORTAR::IntCell&               cell,            //!< mortar integration cell
    DRT::Condition&                condition,       //!< scatra-scatra interface coupling condition
    const Epetra_Vector&           iphinp,          //!< interface state vector
    const DRT::Discretization&     idiscret,        //!< interface discretization
    DRT::Element::LocationArray&   la_slave,        //!< slave-side location array
    DRT::Element::LocationArray&   la_master        //!< master-side location array
    ) const
{
  switch(slaveelement.Shape())
  {
    case DRT::Element::tri3:
    {
      MortarCellCalcElch<DRT::Element::tri3>(
          slaveelement,
          masterelement,
          cell,
          condition,
          iphinp,
          idiscret,
          la_slave,
          la_master
          );
      break;
    }

    default:
    {
      dserror("Invalid slave-side discretization type!");
      break;
    }
  }

  return;
}


/*---------------------------------------------------------------------------------------------------*
 | evaluate single mortar integration cell of particular slave-side discretization type   fang 01/16 |
 *---------------------------------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS>
void SCATRA::MeshtyingStrategyS2IElch::MortarCellCalcElch(
    MORTAR::MortarElement&         slaveelement,    //!< slave-side mortar element
    MORTAR::MortarElement&         masterelement,   //!< master-side mortar element
    MORTAR::IntCell&               cell,            //!< mortar integration cell
    DRT::Condition&                condition,       //!< scatra-scatra interface coupling condition
    const Epetra_Vector&           iphinp,          //!< interface state vector
    const DRT::Discretization&     idiscret,        //!< interface discretization
    DRT::Element::LocationArray&   la_slave,        //!< slave-side location array
    DRT::Element::LocationArray&   la_master        //!< master-side location array
    ) const
{
  switch(masterelement.Shape())
  {
    case DRT::Element::tri3:
    {
      SCATRA::MortarCellCalcElch<distypeS,DRT::Element::tri3>::Instance()->Evaluate(
          *islavematrix_,
          *imastermatrix_,
          *islaveresidual_,
          *imasterresidual_,
          iphinp,
          idiscret,
          condition,
          slaveelement,
          masterelement,
          cell,
          la_slave,
          la_master
          );

      break;
    }

    default:
    {
      dserror("Invalid master-side discretization type!");
      break;
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | singleton access method                                   fang 01/16 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,DRT::Element::DiscretizationType distypeM>
SCATRA::MortarCellCalcElch<distypeS,distypeM>* SCATRA::MortarCellCalcElch<distypeS,distypeM>::Instance(bool create)
{
  static MortarCellCalcElch<distypeS,distypeM>* instance;

  if(create)
  {
    if(instance == NULL)
      instance = new MortarCellCalcElch<distypeS,distypeM>();
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
  Instance(false);

  return;
}


/*---------------------------------------------------------------------------*
 | evaluate and assemble interface linearizations and residuals   fang 01/16 |
 *---------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,DRT::Element::DiscretizationType distypeM>
void SCATRA::MortarCellCalcElch<distypeS,distypeM>::CalcMatAndRhs(
    LINALG::SparseMatrix&                                                                              islavematrix,      //!< linearizations of slave-side residuals
    LINALG::SparseMatrix&                                                                              imastermatrix,     //!< linearizations of master-side residuals
    Epetra_Vector&                                                                                     islaveresidual,    //!< slave-side residual vector
    Epetra_FEVector&                                                                                   imasterresidual,   //!< master-side residual vector
    const DRT::Discretization&                                                                         idiscret,          //!< interface discretization
    DRT::Condition&                                                                                    condition,         //!< scatra-scatra interface coupling condition
    MORTAR::MortarElement&                                                                             slaveelement,      //!< slave-side mortar element
    MORTAR::MortarElement&                                                                             masterelement,     //!< master-side mortar element
    MORTAR::IntCell&                                                                                   cell,              //!< mortar integration cell
    std::vector<LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distypeS>::numNodePerElement,1> >&   ephinp_slave,      //!< state variables at slave-side nodes
    std::vector<LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distypeM>::numNodePerElement,1> >&   ephinp_master,     //!< state variables at master-side nodes
    DRT::Element::LocationArray&                                                                       la_slave,          //!< slave-side location array
    DRT::Element::LocationArray&                                                                       la_master          //!< master-side location array
    ) const
{
  // safety check
  if(slaveelement.NumDofPerNode(*slaveelement.Nodes()[0]) != 2 or masterelement.NumDofPerNode(*masterelement.Nodes()[0]) != 2)
    dserror("Invalid number of degrees of freedom per node!");

  // number of element nodes
  static const int nen_slave = DRT::UTILS::DisTypeToNumNodePerEle<distypeS>::numNodePerElement;
  static const int nen_master = DRT::UTILS::DisTypeToNumNodePerEle<distypeM>::numNodePerElement;

  // access input parameters associated with current condition
  const int kineticmodel = condition.GetInt("kinetic model");
  if(kineticmodel != INPAR::S2I::kinetics_butlervolmer)
    dserror("Invalid kinetic model for scatra-scatra interface coupling!");
  const int nume = condition.GetInt("e-");
  if(not (nume > 0))
    dserror("Charge transfer at electrode-electrolyte interface must involve a positive number of electrons!");
  const std::vector<int>* stoichiometries = condition.GetMutable<std::vector<int> >("stoichiometries");
  if(stoichiometries == NULL)
    dserror("Cannot access vector of stoichiometric coefficients for scatra-scatra interface coupling!");
  if(stoichiometries->size() != 1)
    dserror("Invalid number of stoichiometric coefficients!");
  if((*stoichiometries)[0] != -1)
    dserror("Invalid stoichiometric coefficient!");
  int reactivespecies(abs((*stoichiometries)[0]));
  if(reactivespecies > 1)
    dserror("Charge transfer at electrode-electrolyte interface must not involve more than one reactive species!");
  const double faraday = INPAR::ELCH::faraday_const;
  const double alphaa = condition.GetDouble("alpha_a");
  const double alphac = condition.GetDouble("alpha_c");
  const double kr = condition.GetDouble("k_r");
  if(kr < 0.)
    dserror("Charge transfer constant k_r is negative!");
  const double fns = -1./faraday/nume*(*stoichiometries)[0];

  // access material of slave element
  Teuchos::RCP<const MAT::Electrode> matelectrode = Teuchos::rcp_dynamic_cast<const MAT::Electrode>(Teuchos::rcp_dynamic_cast<DRT::FaceElement>(condition.Geometry()[slaveelement.Id()])->ParentElement()->Material());

  // safety check
  if(matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // extract saturation value of intercalated lithium concentration from electrode material
  const double cmax = matelectrode->CMax();
  if(cmax < 1.e-12)
    dserror("Saturation value c_max of intercalated lithium concentration is too small!");

  // determine quadrature rule
  const DRT::UTILS::IntPointsAndWeights<2> intpoints(DRT::UTILS::intrule_tri_7point);

  // loop over all integration points
  for(int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    // evaluate shape functions and domain integration factor at current integration point
    static LINALG::Matrix<nen_slave,1> funct_slave;
    static LINALG::Matrix<nen_master,1> funct_master;
    const double fac = MortarCellCalc<distypeS,distypeM>::EvalShapeFuncAndDomIntFacAtIntPoint(funct_slave,funct_master,slaveelement,masterelement,cell,intpoints,iquad);

    // overall integration factors
    const double timefacfac = DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->TimeFac()*fac;
    const double timefacrhsfac = DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->TimeFacRhs()*fac;
    if(timefacfac < 0. or timefacrhsfac < 0.)
      dserror("Integration factor is negative!");

    // evaluate state variables at current integration point
    const double eslavephiint = funct_slave.Dot(ephinp_slave[0]);
    const double eslavepotint = funct_slave.Dot(ephinp_slave[1]);
    const double emasterphiint = funct_master.Dot(ephinp_master[0]);
    const double emasterpotint = funct_master.Dot(ephinp_master[1]);

    // extract factor F/RT
    const double frt = DRT::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->FRT();

    // equilibrium electric potential difference and its derivative w.r.t. concentration at electrode surface
    const double epd = matelectrode->ComputeOpenCircuitPotential(eslavephiint,faraday,frt);
    const double epdderiv = matelectrode->ComputeFirstDerivOpenCircuitPotential(eslavephiint,faraday,frt);

    // electrode-electrolyte overpotential at integration point
    const double eta = eslavepotint-emasterpotint-epd;

    const double i0 = kr*faraday*pow(emasterphiint,alphaa)*pow(cmax-eslavephiint,alphaa)*pow(eslavephiint,alphac);
    const double expterm1 = exp(alphaa*frt*eta);
    const double expterm2 = exp(-alphac*frt*eta);
    const double expterm = expterm1-expterm2;

    // safety check
    if(abs(expterm)>1.e5)
      dserror("Overflow of exponential term in Butler-Volmer formulation detected! Value: %lf",expterm);

    // core residual
    const double i = fns*i0*expterm*timefacrhsfac;

    // core linearizations
    const double di_dc_slave = fns*timefacfac*(kr*faraday*pow(emasterphiint,alphaa)*pow(cmax-eslavephiint,alphaa-1.)*pow(eslavephiint,alphac-1.)*(-alphaa*eslavephiint+alphac*(cmax-eslavephiint))*expterm+i0*(-alphaa*frt*epdderiv*expterm1-alphac*frt*epdderiv*expterm2));
    const double di_dc_master = fns*timefacfac*kr*faraday*alphaa*pow(emasterphiint,alphaa-1.)*pow(cmax-eslavephiint,alphaa)*pow(eslavephiint,alphac)*expterm;
    const double di_dpot_slave = fns*timefacfac*i0*(alphaa*frt*expterm1+alphac*frt*expterm2);
    const double di_dpot_master = -di_dpot_slave;

    for (int vi=0; vi<nen_slave; ++vi)
    {
      if(la_slave[0].lmowner_[vi*2] == idiscret.Comm().MyPID())
      {
        const int row_conc_slave = la_slave[0].lm_[vi*2];
        const int row_pot_slave = row_conc_slave+1;

        for (int ui=0; ui<nen_slave; ++ui)
        {
          const int col_conc_slave = la_slave[0].lm_[ui*2];
          const int col_pot_slave = col_conc_slave+1;

          islavematrix.Assemble(funct_slave(vi)*di_dc_slave*funct_slave(ui),row_conc_slave,col_conc_slave);
          islavematrix.Assemble(funct_slave(vi)*nume*di_dc_slave*funct_slave(ui),row_pot_slave,col_conc_slave);
          islavematrix.Assemble(funct_slave(vi)*di_dpot_slave*funct_slave(ui),row_conc_slave,col_pot_slave);
          islavematrix.Assemble(funct_slave(vi)*nume*di_dpot_slave*funct_slave(ui),row_pot_slave,col_pot_slave);
        }

        for(int ui=0; ui<nen_master; ++ui)
        {
          const int col_conc_master = la_master[0].lm_[ui*2];
          const int col_pot_master = col_conc_master+1;

          islavematrix.Assemble(funct_slave(vi)*di_dc_master*funct_master(ui),row_conc_slave,col_conc_master);
          islavematrix.Assemble(funct_slave(vi)*nume*di_dc_master*funct_master(ui),row_pot_slave,col_conc_master);
          islavematrix.Assemble(funct_slave(vi)*di_dpot_master*funct_master(ui),row_conc_slave,col_pot_master);
          islavematrix.Assemble(funct_slave(vi)*nume*di_dpot_master*funct_master(ui),row_pot_slave,col_pot_master);
        }

        if(islaveresidual.SumIntoGlobalValue(row_conc_slave,0,-funct_slave(vi)*i))
          dserror("Assembly into slave-side residual vector not successful!");
        if(islaveresidual.SumIntoGlobalValue(row_pot_slave,0,-funct_slave(vi)*nume*i))
          dserror("Assembly into slave-side residual vector not successful!");
      }
    }

    for(int vi=0; vi<nen_master; ++vi)
    {
      if(slaveelement.Owner() == idiscret.Comm().MyPID())
      {
        const int row_conc_master = la_master[0].lm_[vi*2];
        const int row_pot_master = row_conc_master+1;

        for (int ui=0; ui<nen_slave; ++ui)
        {
          const int col_conc_slave = la_slave[0].lm_[ui*2];
          const int col_pot_slave = col_conc_slave+1;

          imastermatrix.FEAssemble(-funct_master(vi)*di_dc_slave*funct_slave(ui),row_conc_master,col_conc_slave);
          imastermatrix.FEAssemble(-funct_master(vi)*nume*di_dc_slave*funct_slave(ui),row_pot_master,col_conc_slave);
          imastermatrix.FEAssemble(-funct_master(vi)*di_dpot_slave*funct_slave(ui),row_conc_master,col_pot_slave);
          imastermatrix.FEAssemble(-funct_master(vi)*nume*di_dpot_slave*funct_slave(ui),row_pot_master,col_pot_slave);
        }

        for(int ui=0; ui<nen_master; ++ui)
        {
          const int col_conc_master = la_master[0].lm_[ui*2];
          const int col_pot_master = col_conc_master+1;

          imastermatrix.FEAssemble(-funct_master(vi)*di_dc_master*funct_master(ui),row_conc_master,col_conc_master);
          imastermatrix.FEAssemble(-funct_master(vi)*nume*di_dc_master*funct_master(ui),row_pot_master,col_conc_master);
          imastermatrix.FEAssemble(-funct_master(vi)*di_dpot_master*funct_master(ui),row_conc_master,col_pot_master);
          imastermatrix.FEAssemble(-funct_master(vi)*nume*di_dpot_master*funct_master(ui),row_pot_master,col_pot_master);
        }

        const double residual_conc_master = funct_master(vi)*i;
        if(imasterresidual.SumIntoGlobalValues(1,&row_conc_master,&residual_conc_master))
          dserror("Assembly into master-side residual vector not successful!");
        const double residual_pot_master = nume*residual_conc_master;
        if(imasterresidual.SumIntoGlobalValues(1,&row_pot_master,&residual_pot_master))
          dserror("Assembly into master-side residual vector not successful!");
      }
    }
  }

  return;
}
