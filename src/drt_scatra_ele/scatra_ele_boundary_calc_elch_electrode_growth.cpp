/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_elch_electrode_growth.cpp

\brief evaluation of ScaTra boundary elements for isothermal electrodes exhibiting surface layer growth, e.g., lithium plating

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
 */
/*----------------------------------------------------------------------*/
#include "scatra_ele_boundary_calc_elch_electrode_growth.H"

#include "scatra_ele_parameter_elch.H"
#include "scatra_ele_parameter_std.H"
#include "scatra_ele_parameter_timint.H"

#include "../drt_inpar/inpar_s2i.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_mat/electrode.H"

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 01/17 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype>* DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype>::Instance(
    const int                                         numdofpernode,
    const int                                         numscal,
    const std::string&                                disname,
    const ScaTraEleBoundaryCalcElchElectrodeGrowth*   delete_me
    )
{
  static std::map<std::string,ScaTraEleBoundaryCalcElchElectrodeGrowth<distype>* > instances;

  if(delete_me == NULL)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleBoundaryCalcElchElectrodeGrowth<distype>(numdofpernode,numscal,disname);
  }

  else
  {
    for(typename std::map<std::string,ScaTraEleBoundaryCalcElchElectrodeGrowth<distype>* >::iterator i=instances.begin(); i!=instances.end(); ++i)
      if(i->second == delete_me)
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
  }

  return instances[disname];
}


/*----------------------------------------------------------------------*
 | singleton destruction                                     fang 01/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype>::Done()
{
  // delete singleton
  Instance(0,0,"",this);

  return;
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 01/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype>::ScaTraEleBoundaryCalcElchElectrodeGrowth(
    const int            numdofpernode,
    const int            numscal,
    const std::string&   disname
    )
  : // constructor of base class
    myelectrode::ScaTraEleBoundaryCalcElchElectrode(numdofpernode,numscal,disname),

    // initialize member variable
    egrowth_(true)
{
  return;
}


/*--------------------------------------------------------------------------------------------------------------------------*
 | evaluate minimum and maximum interfacial overpotential associated with scatra-scatra interface layer growth   fang 02/18 |
 *--------------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype>::EvaluateMinMaxOverpotential(
    const DRT::FaceElement*        ele,              //!< current boundary element
    Teuchos::ParameterList&        params,           //!< parameter list
    DRT::Discretization&           discretization,   //!< discretization
    DRT::Element::LocationArray&   la                //!< location array
    )
{
  // access material of parent element
  Teuchos::RCP<const MAT::Electrode> matelectrode = Teuchos::rcp_dynamic_cast<const MAT::Electrode>(ele->ParentElement()->Material());
  if (matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  ExtractNodeValues(discretization,la);
  std::vector<LINALG::Matrix<my::nen_,1> > emasterphinp(my::numdofpernode_,LINALG::Matrix<my::nen_,1>(true));
  my::ExtractNodeValues(emasterphinp,discretization,la,"imasterphinp");

  // get scatra-scatra interface coupling condition
  Teuchos::RCP<DRT::Condition> s2icondition = params.get<Teuchos::RCP<DRT::Condition> >("condition");
  if(s2icondition == Teuchos::null)
    dserror("Cannot access scatra-scatra interface coupling condition!");
  if(s2icondition->Type() != DRT::Condition::S2ICouplingGrowth)
    dserror("Received illegal condition type!");

  // access input parameters associated with condition
  const int kineticmodel = s2icondition->GetInt("kinetic model");
  if(kineticmodel != INPAR::S2I::growth_kinetics_butlervolmer)
    dserror("Received illegal kinetic model for scatra-scatra interface coupling involving interface layer growth!");
  const double faraday = myelch::elchparams_->Faraday();
  const double alphaa = s2icondition->GetDouble("alpha_a");
  const double alphac = s2icondition->GetDouble("alpha_c");
  const double kr = s2icondition->GetDouble("k_r");
  if (kr < 0.)
    dserror("Charge transfer constant k_r is negative!");
  const double conductivity_inverse = 1./s2icondition->GetDouble("conductivity");

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid=0; gpid<intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions at current integration point
    my::EvalShapeFuncAndIntFac(intpoints,gpid);

    // evaluate factor F/RT
    const double frt = myelectrode::GetFRT();

    // evaluate dof values at current integration point on present and opposite side of scatra-scatra interface
    const double eslavepotint = my::funct_.Dot(my::ephinp_[1]);
    const double eslavegrowthint = my::funct_.Dot(egrowth_);
    const double emasterphiint = my::funct_.Dot(emasterphinp[0]);
    const double emasterpotint = my::funct_.Dot(emasterphinp[1]);

    // evaluate scatra-scatra interface layer resistance at current integration point
    const double eslaveresistanceint = eslavegrowthint*conductivity_inverse;

    // compute exchange current density
    double i0 = kr*faraday*pow(emasterphiint,alphaa);

    // compute Butler-Volmer current density via Newton-Raphson iteration
    const double i = GetButlerVolmerCurrentDensity(i0,alphaa,alphac,frt,eslavepotint,emasterpotint,0.,eslaveresistanceint,eslavegrowthint,s2icondition);

    // calculate electrode-electrolyte overpotential at integration point
    const double eta = eslavepotint-emasterpotint-i*eslaveresistanceint;

    // check for minimality and update result if applicable
    double& etagrowthmin = params.get<double>("etagrowthmin");
    if(eta < etagrowthmin)
      etagrowthmin = eta;

    // check for maximality and update result if applicable
    double& etagrowthmax = params.get<double>("etagrowthmax");
    if(eta > etagrowthmax)
      etagrowthmax = eta;
  }

  return;
}


/*-------------------------------------------------------------------------------------------*
 | evaluate scatra-scatra interface coupling condition (electrochemistry)         fang 01/17 |
 *-------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype>::EvaluateS2ICoupling(
    const DRT::FaceElement*        ele,              ///< current boundary element
    Teuchos::ParameterList&        params,           ///< parameter list
    DRT::Discretization&           discretization,   ///< discretization
    DRT::Element::LocationArray&   la,               ///< location array
    Epetra_SerialDenseMatrix&      eslavematrix,     ///< element matrix for slave side
    Epetra_SerialDenseMatrix&      emastermatrix,    ///< element matrix for master side
    Epetra_SerialDenseVector&      eslaveresidual    ///< element residual for slave side
    )
{
  // safety checks
  if(my::numscal_ != 1)
    dserror("Invalid number of transported scalars!");
  if(my::numdofpernode_ != 2)
    dserror("Invalid number of degrees of freedom per node!");
  if(myelch::elchparams_->EquPot() != INPAR::ELCH::equpot_divi)
    dserror("Invalid closing equation for electric potential!");

  // access material of parent element
  Teuchos::RCP<const MAT::Electrode> matelectrode = Teuchos::rcp_dynamic_cast<const MAT::Electrode>(ele->ParentElement()->Material());
  if (matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  ExtractNodeValues(discretization,la);
  std::vector<LINALG::Matrix<my::nen_,1> > emasterphinp(my::numdofpernode_,LINALG::Matrix<my::nen_,1>(true));
  my::ExtractNodeValues(emasterphinp,discretization,la,"imasterphinp");

  // get scatra-scatra interface coupling condition
  Teuchos::RCP<DRT::Condition> s2icondition = params.get<Teuchos::RCP<DRT::Condition> >("condition");
  if(s2icondition == Teuchos::null)
    dserror("Cannot access scatra-scatra interface coupling condition!");

  // extract condition type
  const DRT::Condition::ConditionType& s2iconditiontype = s2icondition->Type();
  if(s2iconditiontype != DRT::Condition::S2ICoupling and s2iconditiontype != DRT::Condition::S2ICouplingGrowth)
    dserror("Received illegal condition type!");

  // access input parameters associated with condition
  const int kineticmodel = s2icondition->GetInt("kinetic model");
  if((s2iconditiontype == DRT::Condition::S2ICoupling and kineticmodel != INPAR::S2I::kinetics_butlervolmer) or
     (s2iconditiontype == DRT::Condition::S2ICouplingGrowth and kineticmodel != INPAR::S2I::growth_kinetics_butlervolmer))
    dserror("Received illegal kinetic model for scatra-scatra interface coupling involving interface layer growth!");
  const int nume = s2icondition->GetInt("e-");
  if(nume != 1)
    dserror("Invalid number of electrons involved in charge transfer at electrode-electrolyte interface!");
  const std::vector<int>* stoichiometries = s2icondition->GetMutable<std::vector<int> >("stoichiometries");
  if(stoichiometries == NULL)
    dserror("Cannot access vector of stoichiometric coefficients for scatra-scatra interface coupling!");
  if(stoichiometries->size() != 1)
    dserror("Number of stoichiometric coefficients does not match number of scalars!");
  if((*stoichiometries)[0] != -1)
    dserror("Invalid stoichiometric coefficient!");
  const double faraday = myelch::elchparams_->Faraday();
  const double invF = 1./faraday;
  const double alphaa = s2icondition->GetDouble("alpha_a");
  const double alphac = s2icondition->GetDouble("alpha_c");
  const double kr = s2icondition->GetDouble("k_r");
  if (kr < 0.)
    dserror("Charge transfer constant k_r is negative!");
  const double conductivity_inverse = 1./s2icondition->GetDouble("conductivity");

  // extract saturation value of intercalated lithium concentration from electrode material
  const double cmax = matelectrode->CMax();
  if(cmax < 1.e-12)
    dserror("Saturation value c_max of intercalated lithium concentration is too small!");

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid=0; gpid<intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndIntFac(intpoints,gpid);

    // evaluate overall integration factors
    const double timefacfac = my::scatraparamstimint_->TimeFac()*fac;
    const double timefacrhsfac = my::scatraparamstimint_->TimeFacRhs()*fac;
    if (timefacfac < 0. or timefacrhsfac < 0.)
      dserror("Integration factor is negative!");

    // evaluate factor F/RT
    const double frt = myelectrode::GetFRT();

    // evaluate dof values at current integration point on present and opposite side of scatra-scatra interface
    const double eslavephiint = my::funct_.Dot(my::ephinp_[0]);
    const double eslavepotint = my::funct_.Dot(my::ephinp_[1]);
    const double eslavegrowthint = my::funct_.Dot(egrowth_);
    const double emasterphiint = my::funct_.Dot(emasterphinp[0]);
    const double emasterpotint = my::funct_.Dot(emasterphinp[1]);

    // evaluate scatra-scatra interface layer resistance at current integration point
    const double eslaveresistanceint = eslavegrowthint*conductivity_inverse;

    // equilibrium electric potential difference and its derivative w.r.t. concentration at electrode surface
    const double epd = s2iconditiontype == DRT::Condition::S2ICoupling ? matelectrode->ComputeOpenCircuitPotential(eslavephiint,faraday,frt) : 0.;
    const double epdderiv = matelectrode->ComputeFirstDerivOpenCircuitPotential(eslavephiint,faraday,frt);

    // compute exchange current density
    double i0 = kr*faraday*pow(emasterphiint,alphaa);
    if(s2iconditiontype == DRT::Condition::S2ICoupling and not std::isinf(epd))
      i0 *= pow(cmax-eslavephiint,alphaa)*pow(eslavephiint,alphac);

    // compute Butler-Volmer current density via Newton-Raphson iteration
    const double i = GetButlerVolmerCurrentDensity(i0,alphaa,alphac,frt,eslavepotint,emasterpotint,epd,eslaveresistanceint,eslavegrowthint,s2icondition);

    // continue with evaluation of linearizations and residual contributions only in case of non-zero Butler-Volmer current density
    // to avoid unnecessary effort and to consistently enforce the lithium plating condition
    if(std::abs(i) > 1.e-16)
    {
      // calculate electrode-electrolyte overpotential at integration point and regularization factor
      const double eta = eslavepotint-emasterpotint-epd-i*eslaveresistanceint;
      const double regfac = GetRegularizationFactor(eslavegrowthint,s2icondition,eta);

      // exponential Butler-Volmer terms
      const double expterm1 = exp(alphaa*frt*eta);
      const double expterm2 = exp(-alphac*frt*eta);
      const double expterm = regfac*(expterm1-expterm2);

      // safety check
      if(abs(expterm)>1.e5)
        dserror("Overflow of exponential term in Butler-Volmer formulation detected! Value: %lf",expterm);

      // compute linearizations of Butler-Volmer current density via implicit differentiation, where F = i0*expterm - i = 0
      const double dF_dc_slave = s2iconditiontype == DRT::Condition::S2ICoupling ? kr*faraday*pow(emasterphiint,alphaa)*pow(cmax-eslavephiint,alphaa-1.)*pow(eslavephiint,alphac-1.)*(-alphaa*eslavephiint+alphac*(cmax-eslavephiint))*expterm - i0*(alphaa*frt*epdderiv*expterm1 + alphac*frt*epdderiv*expterm2) : 0.;
      const double dF_dc_master = i0*alphaa/emasterphiint*expterm;
      const double dF_dpot_slave = i0*frt*regfac*(alphaa*expterm1+alphac*expterm2);
      const double dF_dpot_master = -dF_dpot_slave;
      const double dF_di_inverse = -1./(i0*frt*eslaveresistanceint*regfac*(alphaa*expterm1+alphac*expterm2)+1.);
      const double di_dc_slave = -dF_dc_slave*dF_di_inverse;
      const double di_dc_master = -dF_dc_master*dF_di_inverse;
      const double di_dpot_slave = -dF_dpot_slave*dF_di_inverse;
      const double di_dpot_master = -dF_dpot_master*dF_di_inverse;

      // compute linearizations and residual contributions associated with equations for lithium transport
      for(int irow=0; irow<my::nen_; ++irow)
      {
        const int row_conc = irow*2;
        const double funct_irow_invF_timefacfac = my::funct_(irow)*invF*timefacfac;

        for(int icol=0; icol<my::nen_; ++icol)
        {
          const int col_conc = icol*2;
          const int col_pot = col_conc+1;

          eslavematrix(row_conc,col_conc) += funct_irow_invF_timefacfac*di_dc_slave*my::funct_(icol);
          eslavematrix(row_conc,col_pot) += funct_irow_invF_timefacfac*di_dpot_slave*my::funct_(icol);
          emastermatrix(row_conc,col_conc) += funct_irow_invF_timefacfac*di_dc_master*my::funct_(icol);
          emastermatrix(row_conc,col_pot) += funct_irow_invF_timefacfac*di_dpot_master*my::funct_(icol);
        }

        eslaveresidual[row_conc] -= my::funct_(irow)*invF*timefacrhsfac*i;
      }
    } // if(std::abs(i) > 1.e-16)
  } // loop over integration points

  // compute linearizations and residual contributions associated with closing equations for electric potential
  for (int irow=0; irow<my::nen_; ++irow)
  {
    const int row_conc = irow*2;
    const int row_pot = row_conc+1;

    for (int icol=0; icol<my::nen_; ++icol)
    {
      const int col_conc = icol*2;
      const int col_pot = col_conc+1;

      eslavematrix(row_pot,col_conc) += nume*eslavematrix(row_conc,col_conc);
      eslavematrix(row_pot,col_pot) += nume*eslavematrix(row_conc,col_pot);
      emastermatrix(row_pot,col_conc) += nume*emastermatrix(row_conc,col_conc);
      emastermatrix(row_pot,col_pot) += nume*emastermatrix(row_conc,col_pot);
    }

    eslaveresidual[row_pot] += nume*eslaveresidual[row_conc];
  }

  return;
} // DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype>::EvaluateS2ICoupling


/*----------------------------------------------------------------------*
 | evaluate action                                           fang 01/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype>::EvaluateAction(
    DRT::FaceElement*              ele,              //!< boundary element
    Teuchos::ParameterList&        params,           //!< parameter list
    DRT::Discretization&           discretization,   //!< discretization
    SCATRA::BoundaryAction         action,           //!< action
    DRT::Element::LocationArray&   la,               //!< location array
    Epetra_SerialDenseMatrix&      elemat1_epetra,   //!< element matrix 1
    Epetra_SerialDenseMatrix&      elemat2_epetra,   //!< element matrix 2
    Epetra_SerialDenseVector&      elevec1_epetra,   //!< element right-hand side vector 1
    Epetra_SerialDenseVector&      elevec2_epetra,   //!< element right-hand side vector 2
    Epetra_SerialDenseVector&      elevec3_epetra    //!< element right-hand side vector 3
    )
{
  // determine and evaluate action
  switch(action)
  {
    case SCATRA::bd_calc_s2icoupling_growthgrowth:
    {
      EvaluateS2ICouplingGrowthGrowth(
          ele,
          params,
          discretization,
          la,
          elemat1_epetra,
          elevec1_epetra
          );
      break;
    }

    case SCATRA::bd_calc_s2icoupling_growthscatra:
    {
      EvaluateS2ICouplingGrowthScatra(
          ele,
          params,
          discretization,
          la,
          elemat1_epetra,
          elemat2_epetra
          );
      break;
    }

    case SCATRA::bd_calc_s2icoupling_scatragrowth:
    {
      EvaluateS2ICouplingScatraGrowth(
          ele,
          params,
          discretization,
          la,
          elemat1_epetra
          );
      break;
    }

    case SCATRA::bd_calc_elch_minmax_overpotential:
    {
      EvaluateMinMaxOverpotential(
          ele,
          params,
          discretization,
          la
          );
      break;
    }

    default:
    {
      myelch::EvaluateAction(
          ele,
          params,
          discretization,
          action,
          la,
          elemat1_epetra,
          elemat2_epetra,
          elevec1_epetra,
          elevec2_epetra,
          elevec3_epetra
          );
     break;
    }
  } // switch action

  return 0;
}


/*-------------------------------------------------------------------------------------------------------------------------------*
 | evaluate global scatra-growth matrix block for scatra-scatra interface coupling involving interface layer growth   fang 01/17 |
 *-------------------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype>::EvaluateS2ICouplingScatraGrowth(
    const DRT::FaceElement*        ele,              ///< current boundary element
    Teuchos::ParameterList&        params,           ///< parameter list
    DRT::Discretization&           discretization,   ///< discretization
    DRT::Element::LocationArray&   la,               ///< location array
    Epetra_SerialDenseMatrix&      eslavematrix      ///< element matrix for slave side
    )
{
  // access material of parent element
  Teuchos::RCP<const MAT::Electrode> matelectrode = Teuchos::rcp_dynamic_cast<const MAT::Electrode>(ele->ParentElement()->Material());
  if (matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  ExtractNodeValues(discretization,la);
  std::vector<LINALG::Matrix<my::nen_,1> > emasterphinp(my::numdofpernode_,LINALG::Matrix<my::nen_,1>(true));
  my::ExtractNodeValues(emasterphinp,discretization,la,"imasterphinp");

  // get scatra-scatra interface coupling condition
  Teuchos::RCP<DRT::Condition> s2icondition = params.get<Teuchos::RCP<DRT::Condition> >("condition");
  if(s2icondition == Teuchos::null)
    dserror("Cannot access scatra-scatra interface coupling condition!");

  // extract condition type
  const DRT::Condition::ConditionType& s2iconditiontype = s2icondition->Type();
  if(s2iconditiontype != DRT::Condition::S2ICoupling and s2iconditiontype != DRT::Condition::S2ICouplingGrowth)
    dserror("Received illegal condition type!");

  // access input parameters associated with condition
  const int kineticmodel = s2icondition->GetInt("kinetic model");
  if((s2iconditiontype == DRT::Condition::S2ICoupling and kineticmodel != INPAR::S2I::kinetics_butlervolmer) or
     (s2iconditiontype == DRT::Condition::S2ICouplingGrowth and kineticmodel != INPAR::S2I::growth_kinetics_butlervolmer))
    dserror("Received illegal kinetic model for scatra-scatra interface coupling involving interface layer growth!");
  const int nume = s2icondition->GetInt("e-");
  if(nume != 1)
    dserror("Invalid number of electrons involved in charge transfer at electrode-electrolyte interface!");
  const std::vector<int>* stoichiometries = s2icondition->GetMutable<std::vector<int> >("stoichiometries");
  if(stoichiometries == NULL)
    dserror("Cannot access vector of stoichiometric coefficients for scatra-scatra interface coupling!");
  if(stoichiometries->size() != 1)
    dserror("Number of stoichiometric coefficients does not match number of scalars!");
  if((*stoichiometries)[0] != -1)
    dserror("Invalid stoichiometric coefficient!");
  const double faraday = myelch::elchparams_->Faraday();
  const double invF = 1./faraday;
  const double alphaa = s2icondition->GetDouble("alpha_a");
  const double alphac = s2icondition->GetDouble("alpha_c");
  const double kr = s2icondition->GetDouble("k_r");
  if (kr < 0.)
    dserror("Charge transfer constant k_r is negative!");
  const double conductivity_inverse = 1./s2icondition->GetDouble("conductivity");

  // extract saturation value of intercalated lithium concentration from electrode material
  const double cmax = matelectrode->CMax();
  if(cmax < 1.e-12)
    dserror("Saturation value c_max of intercalated lithium concentration is too small!");

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for(int gpid=0; gpid<intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndIntFac(intpoints,gpid);

    // evaluate overall integration factors
    const double timefacfac = my::scatraparamstimint_->TimeFac()*fac;
    const double timefacrhsfac = my::scatraparamstimint_->TimeFacRhs()*fac;
    if (timefacfac < 0. or timefacrhsfac < 0.)
      dserror("Integration factor is negative!");

    // evaluate factor F/RT
    const double frt = myelectrode::GetFRT();

    // evaluate dof values at current integration point on present and opposite side of scatra-scatra interface
    const double eslavephiint = my::funct_.Dot(my::ephinp_[0]);
    const double eslavepotint = my::funct_.Dot(my::ephinp_[1]);
    const double eslavegrowthint = my::funct_.Dot(egrowth_);
    const double emasterphiint = my::funct_.Dot(emasterphinp[0]);
    const double emasterpotint = my::funct_.Dot(emasterphinp[1]);

    // evaluate scatra-scatra interface layer resistance at current integration point
    const double eslaveresistanceint = eslavegrowthint*conductivity_inverse;

    // equilibrium electric potential difference at electrode surface
    const double epd = s2iconditiontype == DRT::Condition::S2ICoupling ? matelectrode->ComputeOpenCircuitPotential(eslavephiint,faraday,frt) : 0.;

    // compute exchange current density
    double i0 = kr*faraday*pow(emasterphiint,alphaa);
    if(s2iconditiontype == DRT::Condition::S2ICoupling and not std::isinf(epd))
      i0 *= pow(cmax-eslavephiint,alphaa)*pow(eslavephiint,alphac);

    // compute Butler-Volmer current density via Newton-Raphson iteration
    const double i = GetButlerVolmerCurrentDensity(i0,alphaa,alphac,frt,eslavepotint,emasterpotint,epd,eslaveresistanceint,eslavegrowthint,s2icondition);

    // continue with evaluation of linearizations only in case of non-zero Butler-Volmer current density
    // to avoid unnecessary effort and to consistently enforce the lithium plating condition
    if(std::abs(i) > 1.e-16)
    {
      // calculate electrode-electrolyte overpotential at integration point, regularization factor and derivative of regularization factor
      const double eta = eslavepotint-emasterpotint-epd-i*eslaveresistanceint;
      const double regfac = GetRegularizationFactor(eslavegrowthint,s2icondition,eta);
      const double regfacderiv = GetRegularizationFactorDerivative(eslavegrowthint,s2icondition,eta);

      // exponential Butler-Volmer terms
      const double expterm1 = exp(alphaa*frt*eta);
      const double expterm2 = exp(-alphac*frt*eta);
      const double expterm = regfac*(expterm1-expterm2);

      // safety check
      if(abs(expterm)>1.e5)
        dserror("Overflow of exponential term in Butler-Volmer formulation detected! Value: %lf",expterm);

      // compute linearization of Butler-Volmer current density w.r.t. scatra-scatra interface layer thickness
      // via implicit differentiation, where F = i0*expterm - i = 0
      const double dF_dgrowth = -i0*i*frt*regfac*conductivity_inverse*(alphaa*expterm1+alphac*expterm2)+regfacderiv*i0*(expterm1+expterm2);
      const double dF_di_inverse = -1./(i0*frt*eslaveresistanceint*regfac*(alphaa*expterm1+alphac*expterm2)+1.);
      const double di_dgrowth = -dF_dgrowth*dF_di_inverse;

      // compute linearizations associated with equations for lithium transport
      for(int irow=0; irow<my::nen_; ++irow)
      {
        const int row_conc = irow*2;
        const double funct_irow_invF_timefacfac = my::funct_(irow)*invF*timefacfac;

        for(int icol=0; icol<my::nen_; ++icol)
          eslavematrix(row_conc,icol) += funct_irow_invF_timefacfac*di_dgrowth*my::funct_(icol);
      }
    } // if(std::abs(i) > 1.e-16)
  } // loop over integration points

  // compute linearizations associated with closing equations for electric potential
  for(int irow=0; irow<my::nen_; ++irow)
  {
    const int row_conc = irow*2;
    const int row_pot = row_conc+1;

    for(int icol=0; icol<my::nen_; ++icol)
      eslavematrix(row_pot,icol) += nume*eslavematrix(row_conc,icol);
  }

  return;
} // DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype>::EvaluateS2ICouplingScatraGrowth


/*-------------------------------------------------------------------------------------------------------------------------------*
 | evaluate global growth-scatra matrix block for scatra-scatra interface coupling involving interface layer growth   fang 01/17 |
 *-------------------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype>::EvaluateS2ICouplingGrowthScatra(
    const DRT::FaceElement*        ele,              ///< current boundary element
    Teuchos::ParameterList&        params,           ///< parameter list
    DRT::Discretization&           discretization,   ///< discretization
    DRT::Element::LocationArray&   la,               ///< location array
    Epetra_SerialDenseMatrix&      eslavematrix,     ///< element matrix for slave side
    Epetra_SerialDenseMatrix&      emastermatrix     ///< element matrix for master side
    )
{
  // access material of parent element
  Teuchos::RCP<const MAT::Electrode> matelectrode = Teuchos::rcp_dynamic_cast<const MAT::Electrode>(ele->ParentElement()->Material());
  if (matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  ExtractNodeValues(discretization,la);
  std::vector<LINALG::Matrix<my::nen_,1> > emasterphinp(my::numdofpernode_,LINALG::Matrix<my::nen_,1>(true));
  my::ExtractNodeValues(emasterphinp,discretization,la,"imasterphinp");

  // get scatra-scatra interface coupling condition
  Teuchos::RCP<DRT::Condition> s2icondition = params.get<Teuchos::RCP<DRT::Condition> >("condition");
  if(s2icondition == Teuchos::null)
    dserror("Cannot access scatra-scatra interface coupling condition!");
  if(s2icondition->Type() != DRT::Condition::S2ICouplingGrowth)
    dserror("Received illegal condition type!");

  // access input parameters associated with condition
  const int kineticmodel = s2icondition->GetInt("kinetic model");
  if(kineticmodel != INPAR::S2I::growth_kinetics_butlervolmer)
    dserror("Received illegal kinetic model for scatra-scatra interface coupling involving interface layer growth!");
  const int nume = s2icondition->GetInt("e-");
  if(nume != 1)
    dserror("Invalid number of electrons involved in charge transfer at electrode-electrolyte interface!");
  const std::vector<int>* stoichiometries = s2icondition->GetMutable<std::vector<int> >("stoichiometries");
  if(stoichiometries == NULL)
    dserror("Cannot access vector of stoichiometric coefficients for scatra-scatra interface coupling!");
  if(stoichiometries->size() != 1)
    dserror("Number of stoichiometric coefficients does not match number of scalars!");
  if((*stoichiometries)[0] != -1)
    dserror("Invalid stoichiometric coefficient!");
  const double faraday = myelch::elchparams_->Faraday();
  const double alphaa = s2icondition->GetDouble("alpha_a");
  const double alphac = s2icondition->GetDouble("alpha_c");
  const double kr = s2icondition->GetDouble("k_r");
  if (kr < 0.)
    dserror("Charge transfer constant k_r is negative!");
  const double conductivity_inverse = 1./s2icondition->GetDouble("conductivity");
  const double factor = s2icondition->GetDouble("molar mass")/(s2icondition->GetDouble("density")*faraday);

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int gpid=0; gpid<intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndIntFac(intpoints,gpid);

    // evaluate overall integration factors
    const double timefacfac = my::scatraparamstimint_->TimeFac()*fac;
    const double timefacrhsfac = my::scatraparamstimint_->TimeFacRhs()*fac;
    if (timefacfac < 0. or timefacrhsfac < 0.)
      dserror("Integration factor is negative!");

    // evaluate factor F/RT
    const double frt = myelectrode::GetFRT();

    // evaluate dof values at current integration point on present and opposite side of scatra-scatra interface
    const double eslavepotint = my::funct_.Dot(my::ephinp_[1]);
    const double eslavegrowthint = my::funct_.Dot(egrowth_);
    const double emasterphiint = my::funct_.Dot(emasterphinp[0]);
    const double emasterpotint = my::funct_.Dot(emasterphinp[1]);

    // evaluate scatra-scatra interface layer resistance at current integration point
    const double eslaveresistanceint = eslavegrowthint*conductivity_inverse;

    // compute exchange current density
    const double i0 = kr*faraday*pow(emasterphiint,alphaa);

    // compute Butler-Volmer current density via Newton-Raphson iteration
    const double i = GetButlerVolmerCurrentDensity(i0,alphaa,alphac,frt,eslavepotint,emasterpotint,0.,eslaveresistanceint,eslavegrowthint,s2icondition);

    // continue with evaluation of linearizations only in case of non-zero Butler-Volmer current density
    // to avoid unnecessary effort and to consistently enforce the lithium plating condition
    if(std::abs(i) > 1.e-16)
    {
      // calculate electrode-electrolyte overpotential at integration point and regularization factor
      const double eta = eslavepotint-emasterpotint-i*eslaveresistanceint;
      const double regfac = GetRegularizationFactor(eslavegrowthint,s2icondition,eta);

      // exponential Butler-Volmer terms
      const double expterm1 = exp(alphaa*frt*eta);
      const double expterm2 = exp(-alphac*frt*eta);
      const double expterm = regfac*(expterm1-expterm2);

      // safety check
      if(abs(expterm)>1.e5)
        dserror("Overflow of exponential term in Butler-Volmer formulation detected! Value: %lf",expterm);

      // compute linearizations of Butler-Volmer current density via implicit differentiation, where F = i0*expterm - i = 0
      const double dF_dc_master = i0*alphaa/emasterphiint*expterm;
      const double dF_dpot_slave = i0*frt*regfac*(alphaa*expterm1+alphac*expterm2);
      const double dF_dpot_master = -dF_dpot_slave;
      const double dF_di_inverse = -1./(i0*frt*eslaveresistanceint*regfac*(alphaa*expterm1+alphac*expterm2)+1.);
      const double di_dc_master = -dF_dc_master*dF_di_inverse;
      const double di_dpot_slave = -dF_dpot_slave*dF_di_inverse;
      const double di_dpot_master = -dF_dpot_master*dF_di_inverse;

      // compute linearizations associated with equation for scatra-scatra interface layer growth
      for(int irow=0; irow<my::nen_; ++irow)
      {
        const double funct_irow_factor_timefacfac = my::funct_(irow)*factor*timefacfac;

        for(int icol=0; icol<my::nen_; ++icol)
        {
          const int col_conc = icol*2;
          const int col_pot = col_conc+1;

          eslavematrix(irow,col_pot) += funct_irow_factor_timefacfac*di_dpot_slave*my::funct_(icol);
          emastermatrix(irow,col_conc) += funct_irow_factor_timefacfac*di_dc_master*my::funct_(icol);
          emastermatrix(irow,col_pot) += funct_irow_factor_timefacfac*di_dpot_master*my::funct_(icol);
        }
      }
    } // if(std::abs(i) > 1.e-16)
  } // loop over integration points

  return;
}


/*-------------------------------------------------------------------------------------------------------------------------------*
 | evaluate global growth-growth matrix block for scatra-scatra interface coupling involving interface layer growth   fang 01/17 |
 *-------------------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype>::EvaluateS2ICouplingGrowthGrowth(
    const DRT::FaceElement*        ele,              ///< current boundary element
    Teuchos::ParameterList&        params,           ///< parameter list
    DRT::Discretization&           discretization,   ///< discretization
    DRT::Element::LocationArray&   la,               ///< location array
    Epetra_SerialDenseMatrix&      eslavematrix,     ///< element matrix for slave side
    Epetra_SerialDenseVector&      eslaveresidual    ///< element residual for slave side
    )
{
  // access material of parent element
  Teuchos::RCP<const MAT::Electrode> matelectrode = Teuchos::rcp_dynamic_cast<const MAT::Electrode>(ele->ParentElement()->Material());
  if (matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  ExtractNodeValues(discretization,la);
  std::vector<LINALG::Matrix<my::nen_,1> > emasterphinp(my::numdofpernode_,LINALG::Matrix<my::nen_,1>(true));
  LINALG::Matrix<my::nen_,1> eslavegrowthhist(true);
  my::ExtractNodeValues(emasterphinp,discretization,la,"imasterphinp");
  my::ExtractNodeValues(eslavegrowthhist,discretization,la,"growthhist",2);

  // get scatra-scatra interface coupling condition
  Teuchos::RCP<DRT::Condition> s2icondition = params.get<Teuchos::RCP<DRT::Condition> >("condition");
  if(s2icondition == Teuchos::null)
    dserror("Cannot access scatra-scatra interface coupling condition!");
  if(s2icondition->Type() != DRT::Condition::S2ICouplingGrowth)
    dserror("Received illegal condition type!");

  // access input parameters associated with condition
  const int kineticmodel = s2icondition->GetInt("kinetic model");
  if(kineticmodel != INPAR::S2I::growth_kinetics_butlervolmer)
    dserror("Received illegal kinetic model for scatra-scatra interface coupling involving interface layer growth!");
  const int nume = s2icondition->GetInt("e-");
  if(nume != 1)
    dserror("Invalid number of electrons involved in charge transfer at electrode-electrolyte interface!");
  const std::vector<int>* stoichiometries = s2icondition->GetMutable<std::vector<int> >("stoichiometries");
  if(stoichiometries == NULL)
    dserror("Cannot access vector of stoichiometric coefficients for scatra-scatra interface coupling!");
  if(stoichiometries->size() != 1)
    dserror("Number of stoichiometric coefficients does not match number of scalars!");
  if((*stoichiometries)[0] != -1)
    dserror("Invalid stoichiometric coefficient!");
  const double faraday = myelch::elchparams_->Faraday();
  const double alphaa = s2icondition->GetDouble("alpha_a");
  const double alphac = s2icondition->GetDouble("alpha_c");
  const double kr = s2icondition->GetDouble("k_r");
  if (kr < 0.)
    dserror("Charge transfer constant k_r is negative!");
  const double conductivity_inverse = 1./s2icondition->GetDouble("conductivity");
  const double factor = s2icondition->GetDouble("molar mass")/(s2icondition->GetDouble("density")*faraday);

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for(int gpid=0; gpid<intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndIntFac(intpoints,gpid);

    // evaluate mass matrix
    for(int irow=0; irow<my::nen_; ++irow)
       for (int icol=0; icol<my::nen_; ++icol)
         eslavematrix(irow,icol) += my::funct_(irow)*my::funct_(icol)*fac;

    // evaluate overall integration factors
    const double timefacfac = my::scatraparamstimint_->TimeFac()*fac;
    const double timefacrhsfac = my::scatraparamstimint_->TimeFacRhs()*fac;
    if (timefacfac < 0. or timefacrhsfac < 0.)
      dserror("Integration factor is negative!");

    // evaluate factor F/RT
    const double frt = myelectrode::GetFRT();

    // evaluate dof values at current integration point on present and opposite side of scatra-scatra interface
    const double eslavepotint = my::funct_.Dot(my::ephinp_[1]);
    const double eslavegrowthint = my::funct_.Dot(egrowth_);
    const double eslavegrowthhistint = my::funct_.Dot(eslavegrowthhist);
    const double emasterphiint = my::funct_.Dot(emasterphinp[0]);
    const double emasterpotint = my::funct_.Dot(emasterphinp[1]);

    // evaluate scatra-scatra interface layer resistance at current integration point
    const double eslaveresistanceint = eslavegrowthint*conductivity_inverse;

    // compute exchange current density
    const double i0 = kr*faraday*pow(emasterphiint,alphaa);

    // compute Butler-Volmer current density via Newton-Raphson iteration
    const double i = GetButlerVolmerCurrentDensity(i0,alphaa,alphac,frt,eslavepotint,emasterpotint,0.,eslaveresistanceint,eslavegrowthint,s2icondition);

    // continue with evaluation of linearizations and residual contributions only in case of non-zero Butler-Volmer current density
    // to avoid unnecessary effort and to consistently enforce the lithium plating condition
    if(std::abs(i) > 1.e-16)
    {
      // calculate electrode-electrolyte overpotential at integration point, regularization factor and derivative of regularization factor
      const double eta = eslavepotint-emasterpotint-i*eslaveresistanceint;
      const double regfac = GetRegularizationFactor(eslavegrowthint,s2icondition,eta);
      const double regfacderiv = GetRegularizationFactorDerivative(eslavegrowthint,s2icondition,eta);

      // exponential Butler-Volmer terms
      const double expterm1 = exp(alphaa*frt*eta);
      const double expterm2 = exp(-alphac*frt*eta);
      const double expterm = regfac*(expterm1-expterm2);

      // safety check
      if (abs(expterm)>1.e5)
        dserror("Overflow of exponential term in Butler-Volmer formulation detected! Value: %lf",expterm);

      // compute linearization of Butler-Volmer current density w.r.t. scatra-scatra interface layer thickness
      // via implicit differentiation, where F = i0*expterm - i = 0
      const double dF_dgrowth = -i0*i*frt*regfac*conductivity_inverse*(alphaa*expterm1+alphac*expterm2)+regfacderiv*i0*(expterm1+expterm2);
      const double dF_di_inverse = -1./(i0*frt*eslaveresistanceint*regfac*(alphaa*expterm1+alphac*expterm2)+1.);
      const double di_dgrowth = -dF_dgrowth*dF_di_inverse;

      // compute linearizations and residual contributions associated with equation for scatra-scatra interface layer growth
      for(int irow=0; irow<my::nen_; ++irow)
      {
        const double funct_irow_factor_timefacfac = my::funct_(irow)*factor*timefacfac;

        for(int icol=0; icol<my::nen_; ++icol)
          eslavematrix(irow,icol) += funct_irow_factor_timefacfac*di_dgrowth*my::funct_(icol);

        eslaveresidual[irow] -= my::funct_(irow)*(eslavegrowthint-eslavegrowthhistint)*fac;
        eslaveresidual[irow] -= my::funct_(irow)*factor*i*timefacrhsfac;
      }
    } // if(std::abs(i) > 1.e-16)
  } // loop over integration points

  return;
}


/*-----------------------------------------------------------------------------*
 | extract nodal state variables associated with boundary element   fang 01/17 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype>::ExtractNodeValues(
    const DRT::Discretization&     discretization,  //!< discretization
    DRT::Element::LocationArray&   la               //!< location array
    )
{
  // call base class routine
  my::ExtractNodeValues(discretization,la);

  // extract nodal growth variables associated with boundary element
  my::ExtractNodeValues(egrowth_,discretization,la,"growth",2);

  return;
}


/*---------------------------------------------------------------------------------*
 | compute Butler-Volmer current density via Newton-Raphson iteration   fang 01/17 |
 *---------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype>::GetButlerVolmerCurrentDensity(
    const double&                         i0,           //!< exchange current density
    const double&                         alphaa,       //!< anodic transfer coefficient
    const double&                         alphac,       //!< cathodic transfer coefficient
    const double&                         frt,          //!< factor F/RT
    const double&                         pot_ed,       //!< electrode-side electric potential
    const double&                         pot_el,       //!< electrolyte-side electric potential
    const double&                         epd,          //!< half-cell open-circuit potential
    const double&                         resistance,   //!< scatra-scatra interface layer resistance
    const double&                         thickness,    //!< scatra-scatra interface layer thickness
    const Teuchos::RCP<DRT::Condition>&   condition     //!< scatra-scatra interface coupling condition
    ) const
{
  // initialize Butler-Volmer current density
  double i(0.);

  // compute Butler-Volmer current density in case of physically reasonable half-cell open-circuit potential
  if(not std::isinf(epd))
  {
    // compute initial guess of Butler-Volmer current density, neglecting overpotential due to scatra-scatra interface layer resistance
    double eta = pot_ed-pot_el-epd;
    double regfac = GetRegularizationFactor(thickness,condition,eta);
    i = i0*regfac*(exp(alphaa*frt*eta)-exp(-alphac*frt*eta));

    // initialize Newton-Raphson iteration counter
    unsigned iternum(0);

    // apply Newton-Raphson method to compute Butler-Volmer current density, involving overpotential due to scatra-scatra interface layer resistance
    while(true)
    {
      // increment counter
      ++iternum;

      // compute current Newton-Raphson residual
      eta = pot_ed-pot_el-epd-resistance*i;
      regfac = GetRegularizationFactor(thickness,condition,eta);
      const double expterm1 = exp(alphaa*frt*eta);
      const double expterm2 = exp(-alphac*frt*eta);
      const double residual = i0*regfac*(expterm1-expterm2)-i;

      // convergence check
      if(std::abs(residual) < my::scatraparams_->IntLayerGrowthConvTol())
        break;
      else if(iternum == my::scatraparams_->IntLayerGrowthIteMax())
        dserror("Local Newton-Raphson iteration for Butler-Volmer current density did not converge!");

      // compute linearization of current Newton-Raphson residual w.r.t. Butler-Volmer current density
      const double linearization = -i0*resistance*frt*regfac*(alphaa*expterm1+alphac*expterm2)-1.;

      // update Butler-Volmer current density
      i -= residual/linearization;
    }

    // enforce plating condition, i.e., consider initial lithium plating only in case of negative overpotential
    if(std::abs(regfac) < 1.e-16 and eta >= 0.)
      i = 0.;
  }

  return i;
}


/*----------------------------------------------------------------------*
 | compute regularization factor for lithium stripping       fang 01/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype>::GetRegularizationFactor(
    const double&                         thickness,   //!< scatra-scatra interface layer thickness
    const Teuchos::RCP<DRT::Condition>&   condition,   //!< scatra-scatra interface coupling condition
    const double                          eta          //!< electrode-electrolyte overpotential at integration point
    ) const
{
  // initialize regularization factor
  double regfac(1.);

  // actually compute regularization factor if lithium stripping is relevant
  if(condition->Type() == DRT::Condition::S2ICouplingGrowth and eta > 0.)
  {
    // extract regularization type from scatra-scatra interface coupling condition
    const std::string* regtype = condition->Get<std::string>("regularization type");
    if(!regtype)
      dserror("Can't extract regularization type from condition for scatra-scatra interface layer growth!");

    // extract regularization parameter from scatra-scatra interface coupling condition
    const double regpar = condition->GetDouble("regularization parameter");
    if(regpar < 0.)
      dserror("Regularization parameter for lithium stripping must not be negative!");

    // polynomial regularization, cf. Hein, Latz, Electrochimica Acta 201 (2016) 354-365
    if(*regtype == "polynomial")
    {
      // use regularization parameter if specified, otherwise take default value according to reference
      const double thickness0 = regpar > 0. ? regpar : 4.8e-7;

      // compute regularization factor
      regfac = thickness <= 0. ? 0. : pow(thickness,4)/(pow(thickness,4)+pow(thickness0,4));
    }

    // trigonometrical regularization involving (co)sine half-wave
    else if(*regtype == "trigonometrical")
    {
      // use regularization parameter if specified, otherwise take lithium atom diameter as default value
      const double thickness_regend = regpar > 0. ? regpar : 2.9e-7;

      // compute regularization factor
      if(thickness <= 0.)
        regfac = 0.;
      else if(thickness < thickness_regend)
        regfac = 0.5*cos(thickness/thickness_regend*M_PI-M_PI)+0.5;
    }

    // non-regularized Heaviside function
    else if(*regtype == "none")
    {
      if(thickness <= 0.)
        regfac = 0.;
    }

    // safety check
    else
      dserror("Invalid type of regularization for lithium stripping!");
  }

  return regfac;
}


/*-------------------------------------------------------------------------------------------------------------------*
 | compute derivative of regularization factor for lithium stripping w.r.t. thickness of plated lithium   fang 01/17 |
 *-------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<distype>::GetRegularizationFactorDerivative(
    const double&                         thickness,   //!< scatra-scatra interface layer thickness
    const Teuchos::RCP<DRT::Condition>&   condition,   //!< scatra-scatra interface coupling condition
    const double                          eta          //!< electrode-electrolyte overpotential at integration point
    ) const
{
  // initialize derivative of regularization factor
  double regfacderiv(0.);

  // actually compute derivative of regularization factor if lithium stripping is relevant
  if (condition->Type() == DRT::Condition::S2ICouplingGrowth and thickness > 0. and eta > 0.)
  {
    // extract regularization type from scatra-scatra interface coupling condition
    const std::string* regtype = condition->Get<std::string>("regularization type");
    if(!regtype)
      dserror("Can't extract regularization type from condition for scatra-scatra interface layer growth!");

    // extract regularization parameter from scatra-scatra interface coupling condition
    const double regpar = condition->GetDouble("regularization parameter");
    if(regpar < 0.)
      dserror("Regularization parameter for lithium stripping must not be negative!");

    // polynomial regularization, cf. Hein, Latz, Electrochimica Acta 201 (2016) 354-365
    if (*regtype == "polynomial")
    {
      // use regularization parameter if specified, otherwise take default value according to reference
      const double thickness0 = regpar > 0. ? regpar : 4.8e-7;

      // compute derivative of regularization factor
      const double thickness0_pow4 = pow(thickness0,4);
      regfacderiv = 4.*pow(thickness,3)*thickness0_pow4/pow(thickness0_pow4+pow(thickness,4),2);
    }

    // trigonometrical regularization involving (co)sine half-wave
    else if(*regtype == "trigonometrical")
    {
      // use regularization parameter if specified, otherwise take lithium atom diameter as default value
      const double thickness_regend = regpar > 0. ? regpar : 2.9e-7;

      // compute derivative of regularization factor
      const double thickness_regend_inverse = 1./thickness_regend;
      if(thickness < thickness_regend)
        regfacderiv = 0.5*sin(thickness*thickness_regend_inverse*M_PI)*M_PI*thickness_regend_inverse;
    }

    // non-regularized Heaviside function
    else if(*regtype == "none")
    {
      // do nothing and retain derivative as initialized, since non-regularized Heaviside function cannot be properly differentiated
    }

    // safety check
    else
      dserror("Invalid type of regularization for lithium stripping!");
  }

  return regfacderiv;
}

// template classes
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<DRT::Element::quad4>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<DRT::Element::line3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<DRT::Element::nurbs3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeGrowth<DRT::Element::nurbs9>;
