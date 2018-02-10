/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_elch_electrode.cpp

\brief evaluation of ScaTra boundary elements for isothermal electrodes

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
 */
/*----------------------------------------------------------------------*/
#include "scatra_ele_boundary_calc_elch_electrode.H"
#include "scatra_ele_parameter_elch.H"
#include "scatra_ele_parameter_timint.H"

#include "../drt_inpar/inpar_s2i.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_mat/electrode.H"

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 02/15 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>* DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::Instance(
    const int numdofpernode,
    const int numscal,
    const std::string& disname,
    const ScaTraEleBoundaryCalcElchElectrode* delete_me
    )
{
  static std::map<std::string,ScaTraEleBoundaryCalcElchElectrode<distype>* >  instances;

  if(delete_me == NULL)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleBoundaryCalcElchElectrode<distype>(numdofpernode,numscal,disname);
  }

  else
  {
    for( typename std::map<std::string,ScaTraEleBoundaryCalcElchElectrode<distype>* >::iterator i=instances.begin(); i!=instances.end(); ++i )
      if ( i->second == delete_me )
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
  }

  return instances[disname];
}


/*----------------------------------------------------------------------*
 | singleton destruction                                     fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::Done()
{
  // delete singleton
  Instance(0,0,"",this);

  return;
}


/*----------------------------------------------------------------------*
 | protected constructor for singletons                      fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::ScaTraEleBoundaryCalcElchElectrode(
    const int numdofpernode,
    const int numscal,
    const std::string& disname)
  : // constructor of base class
    myelch::ScaTraEleBoundaryCalcElch(numdofpernode,numscal,disname)
{
  return;
}


/*-------------------------------------------------------------------------------------*
 | evaluate scatra-scatra interface coupling condition (electrochemistry)   fang 04/15 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::EvaluateS2ICoupling(
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
  if(matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  this->ExtractNodeValues(discretization,la);
  std::vector<LINALG::Matrix<my::nen_,1> > emasterphinp(my::numdofpernode_,LINALG::Matrix<my::nen_,1>(true));
  my::ExtractNodeValues(emasterphinp,discretization,la,"imasterphinp");

  // get current scatra-scatra interface coupling condition
  Teuchos::RCP<DRT::Condition> s2icondition = params.get<Teuchos::RCP<DRT::Condition> >("condition");
  if(s2icondition == Teuchos::null)
    dserror("Cannot access scatra-scatra interface coupling condition!");

  // dummy element matrix and vector
  Epetra_SerialDenseMatrix dummymatrix;
  Epetra_SerialDenseVector dummyvector;

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

    EvaluateS2ICouplingAtIntegrationPoint<distype>(
        *s2icondition,
        matelectrode,
        my::ephinp_,
        emasterphinp,
        my::funct_,
        my::funct_,
        my::funct_,
        my::funct_,
        timefacfac,
        timefacrhsfac,
        GetFRT(),
        eslavematrix,
        emastermatrix,
        dummymatrix,
        dummymatrix,
        eslaveresidual,
        dummyvector
        );
  } // loop over integration points

  return;
} // DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::EvaluateS2ICoupling


/*---------------------------------------------------------------------------------------*
 | evaluate scatra-scatra interface coupling condition at integration point   fang 05/16 |
 *---------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
template <DRT::Element::DiscretizationType distype_master>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::EvaluateS2ICouplingAtIntegrationPoint(
    DRT::Condition&                                                                                                s2icondition,    //!< scatra-scatra interface coupling condition
    const Teuchos::RCP<const MAT::Electrode>&                                                                      matelectrode,    //!< electrode material
    const std::vector<LINALG::Matrix<my::nen_,1> >&                                                                eslavephinp,     //!< state variables at slave-side nodes
    const std::vector<LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement,1> >&   emasterphinp,    //!< state variables at master-side nodes
    const LINALG::Matrix<my::nen_,1>&                                                                              funct_slave,     //!< slave-side shape function values
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement,1>&                 funct_master,    //!< master-side shape function values
    const LINALG::Matrix<my::nen_,1>&                                                                              test_slave,      //!< slave-side test function values
    const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement,1>&                 test_master,     //!< master-side test function values
    const double                                                                                                   timefacfac,      //!< time-integration factor times domain-integration factor
    const double                                                                                                   timefacrhsfac,   //!< time-integration factor for right-hand side times domain-integration factor
    const double                                                                                                   frt,             //!< factor F/(RT)
    Epetra_SerialDenseMatrix&                                                                                      k_ss,            //!< linearizations of slave-side residuals w.r.t. slave-side dofs
    Epetra_SerialDenseMatrix&                                                                                      k_sm,            //!< linearizations of slave-side residuals w.r.t. master-side dofs
    Epetra_SerialDenseMatrix&                                                                                      k_ms,            //!< linearizations of master-side residuals w.r.t. slave-side dofs
    Epetra_SerialDenseMatrix&                                                                                      k_mm,            //!< linearizations of master-side residuals w.r.t. master-side dofs
    Epetra_SerialDenseVector&                                                                                      r_s,             //!< slave-side residual vector
    Epetra_SerialDenseVector&                                                                                      r_m              //!< master-side residual vector
    )
{
  // number of nodes of master-side mortar element
  const int nen_master = DRT::UTILS::DisTypeToNumNodePerEle<distype_master>::numNodePerElement;

  // evaluate dof values at current integration point on present and opposite side of scatra-scatra interface
  const double eslavephiint = funct_slave.Dot(eslavephinp[0]);
  const double eslavepotint = funct_slave.Dot(eslavephinp[1]);
  const double emasterphiint = funct_master.Dot(emasterphinp[0]);
  const double emasterpotint = funct_master.Dot(emasterphinp[1]);

  // compute matrix and vector contributions according to kinetic model for current scatra-scatra interface coupling condition
  const int kinmodel = s2icondition.GetInt("kinetic model");
  switch(kinmodel)
  {
    // Butler-Volmer kinetics
    case INPAR::S2I::kinetics_butlervolmer:
    case INPAR::S2I::kinetics_butlervolmerpeltier:
    case INPAR::S2I::kinetics_butlervolmerreduced:
    {
      // access input parameters associated with current condition
      const int nume = s2icondition.GetInt("e-");
      if(nume != 1)
        dserror("Invalid number of electrons involved in charge transfer at electrode-electrolyte interface!");
      const std::vector<int>* stoichiometries = s2icondition.GetMutable<std::vector<int> >("stoichiometries");
      if(stoichiometries == NULL)
        dserror("Cannot access vector of stoichiometric coefficients for scatra-scatra interface coupling!");
      if(stoichiometries->size() != 1)
        dserror("Number of stoichiometric coefficients does not match number of scalars!");
      if((*stoichiometries)[0] != -1)
        dserror("Invalid stoichiometric coefficient!");
      const double faraday = INPAR::ELCH::faraday_const;
      const double alphaa = s2icondition.GetDouble("alpha_a");
      const double alphac = s2icondition.GetDouble("alpha_c");
      const double kr = s2icondition.GetDouble("k_r");
      if(kr < 0.)
        dserror("Charge transfer constant k_r is negative!");

      // extract saturation value of intercalated lithium concentration from electrode material
      const double cmax = matelectrode->CMax();
      if(cmax < 1.e-12)
        dserror("Saturation value c_max of intercalated lithium concentration is too small!");

      // equilibrium electric potential difference at electrode surface
      const double epd = matelectrode->ComputeOpenCircuitPotential(eslavephiint,faraday,frt);

      // skip further computation in case equilibrium electric potential difference is outside physically meaningful range
      if(not std::isinf(epd))
      {
        // derivative of equilibrium electric potential difference w.r.t. concentration at electrode surface
        const double epdderiv = matelectrode->ComputeFirstDerivOpenCircuitPotential(eslavephiint,faraday,frt);

        // electrode-electrolyte overpotential at integration point
        const double eta = eslavepotint-emasterpotint-epd;

        // Butler-Volmer exchange mass flux density
        const double j0(kinmodel == INPAR::S2I::kinetics_butlervolmerreduced ? kr : kr*pow(emasterphiint,alphaa)*pow(cmax-eslavephiint,alphaa)*pow(eslavephiint,alphac));

        // exponential Butler-Volmer terms
        const double expterm1 = exp(alphaa*frt*eta);
        const double expterm2 = exp(-alphac*frt*eta);
        const double expterm = expterm1-expterm2;

        // safety check
        if(abs(expterm)>1.e5)
          dserror("Overflow of exponential term in Butler-Volmer formulation detected! Value: %lf",expterm);

        // core residual term associated with Butler-Volmer mass flux density
        const double j = j0*expterm*timefacrhsfac;

        // core linearizations associated with Butler-Volmer mass flux density
        const double dj_dc_slave(kinmodel == INPAR::S2I::kinetics_butlervolmerreduced ? timefacfac*j0*frt*epdderiv*(-alphaa*expterm1-alphac*expterm2) : timefacfac*(kr*pow(emasterphiint,alphaa)*pow(cmax-eslavephiint,alphaa-1.)*pow(eslavephiint,alphac-1.)*(-alphaa*eslavephiint+alphac*(cmax-eslavephiint))*expterm+j0*frt*epdderiv*(-alphaa*expterm1-alphac*expterm2)));
        const double dj_dc_master(kinmodel == INPAR::S2I::kinetics_butlervolmerreduced ? 0.0 : timefacfac*j0*alphaa/emasterphiint*expterm);
        const double dj_dpot_slave = timefacfac*j0*(alphaa*frt*expterm1+alphac*frt*expterm2);
        const double dj_dpot_master = -dj_dpot_slave;

        if(k_ss.M() and k_sm.M() and r_s.Length())
        {
          for(int vi=0; vi<my::nen_; ++vi)
          {
            const int row_conc = vi*2;
            const int row_pot = row_conc+1;

            for(int ui=0; ui<my::nen_; ++ui)
            {
              const int col_conc = ui*2;
              const int col_pot = col_conc+1;

              k_ss(row_conc,col_conc) += test_slave(vi)*dj_dc_slave*funct_slave(ui);
              k_ss(row_conc,col_pot) += test_slave(vi)*dj_dpot_slave*funct_slave(ui);
              k_ss(row_pot,col_conc) += nume*test_slave(vi)*dj_dc_slave*funct_slave(ui);
              k_ss(row_pot,col_pot) += nume*test_slave(vi)*dj_dpot_slave*funct_slave(ui);
            }

            for(int ui=0; ui<nen_master; ++ui)
            {
              const int col_conc = ui*2;
              const int col_pot = col_conc+1;

              k_sm(row_conc,col_conc) += test_slave(vi)*dj_dc_master*funct_master(ui);
              k_sm(row_conc,col_pot) += test_slave(vi)*dj_dpot_master*funct_master(ui);
              k_sm(row_pot,col_conc) += nume*test_slave(vi)*dj_dc_master*funct_master(ui);
              k_sm(row_pot,col_pot) += nume*test_slave(vi)*dj_dpot_master*funct_master(ui);
            }

            r_s[row_conc] -= test_slave(vi)*j;
            r_s[row_pot] -= nume*test_slave(vi)*j;
          }
        }
        else if(k_ss.M() or k_sm.M() or r_s.Length())
          dserror("Must provide both slave-side matrices and slave-side vector or none of them!");


        if(k_ms.M() and k_mm.M() and r_m.Length())
        {
          for(int vi=0; vi<nen_master; ++vi)
          {
            const int row_conc = vi*2;
            const int row_pot = row_conc+1;

            for(int ui=0; ui<my::nen_; ++ui)
            {
              const int col_conc = ui*2;
              const int col_pot = col_conc+1;

              k_ms(row_conc,col_conc) -= test_master(vi)*dj_dc_slave*funct_slave(ui);
              k_ms(row_conc,col_pot) -= test_master(vi)*dj_dpot_slave*funct_slave(ui);
              k_ms(row_pot,col_conc) -= nume*test_master(vi)*dj_dc_slave*funct_slave(ui);
              k_ms(row_pot,col_pot) -= nume*test_master(vi)*dj_dpot_slave*funct_slave(ui);
            }

            for(int ui=0; ui<nen_master; ++ui)
            {
              const int col_conc = ui*2;
              const int col_pot = col_conc+1;

              k_mm(row_conc,col_conc) -= test_master(vi)*dj_dc_master*funct_master(ui);
              k_mm(row_conc,col_pot) -= test_master(vi)*dj_dpot_master*funct_master(ui);
              k_mm(row_pot,col_conc) -= nume*test_master(vi)*dj_dc_master*funct_master(ui);
              k_mm(row_pot,col_pot) -= nume*test_master(vi)*dj_dpot_master*funct_master(ui);
            }

            r_m[row_conc] += test_master(vi)*j;
            r_m[row_pot] += nume*test_master(vi)*j;
          }
        }
        else if(k_ms.M() or k_mm.M() or r_m.Length())
          dserror("Must provide both master-side matrices and master-side vector or none of them!");
      }

      break;
    }

    default:
    {
      dserror("Kinetic model for scatra-scatra interface coupling is not yet implemented!");
      break;
    }
  } // switch(kineticmodel)

  return;
}


/*---------------------------------------------------------------------------------------------------------------------------*
 | evaluate off-diagonal system matrix contributions associated with scatra-scatra interface coupling condition   fang 11/17 |
 *---------------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::EvaluateS2ICouplingOD(
    const DRT::FaceElement*        ele,              ///< current boundary element
    Teuchos::ParameterList&        params,           ///< parameter list
    DRT::Discretization&           discretization,   ///< discretization
    DRT::Element::LocationArray&   la,               ///< location array
    Epetra_SerialDenseMatrix&      eslavematrix      ///< element matrix for slave side
    )
{
  // access material of parent element
  Teuchos::RCP<const MAT::Electrode> matelectrode = Teuchos::rcp_dynamic_cast<const MAT::Electrode>(ele->ParentElement()->Material());
  if(matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  this->ExtractNodeValues(discretization,la);
  std::vector<LINALG::Matrix<my::nen_,1> > emasterphinp(my::numdofpernode_,LINALG::Matrix<my::nen_,1>(true));
  my::ExtractNodeValues(emasterphinp,discretization,la,"imasterphinp");

  // get current scatra-scatra interface coupling condition
  Teuchos::RCP<DRT::Condition> s2icondition = params.get<Teuchos::RCP<DRT::Condition> >("condition");
  if(s2icondition == Teuchos::null)
    dserror("Cannot access scatra-scatra interface coupling condition!");

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for(int gpid=0; gpid<intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions at current integration point
    my::EvalShapeFuncAndIntFac(intpoints,gpid);

    // evaluate shape derivatives
    static LINALG::Matrix<my::nsd_+1,my::nen_> shapederivatives;
    my::EvalShapeDerivatives(shapederivatives);

    // evaluate overall integration factor
    const double timefacwgt = my::scatraparamstimint_->TimeFac()*intpoints.IP().qwgt[gpid];
    if(timefacwgt < 0.)
      dserror("Integration factor is negative!");

    // evaluate dof values at current integration point on present and opposite side of scatra-scatra interface
    const double eslavephiint = my::funct_.Dot(my::ephinp_[0]);
    const double eslavepotint = my::funct_.Dot(my::ephinp_[1]);
    const double emasterphiint = my::funct_.Dot(emasterphinp[0]);
    const double emasterpotint = my::funct_.Dot(emasterphinp[1]);

    // compute matrix and vector contributions according to kinetic model for current scatra-scatra interface coupling condition
    const int kinmodel = s2icondition->GetInt("kinetic model");
    switch(kinmodel)
    {
      // Butler-Volmer kinetics
      case INPAR::S2I::kinetics_butlervolmer:
      case INPAR::S2I::kinetics_butlervolmerreduced:
      {
        // access input parameters associated with current condition
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
        const double faraday = INPAR::ELCH::faraday_const;
        const double alphaa = s2icondition->GetDouble("alpha_a");
        const double alphac = s2icondition->GetDouble("alpha_c");
        const double kr = s2icondition->GetDouble("k_r");
        if(kr < 0.)
          dserror("Charge transfer constant k_r is negative!");

        // extract saturation value of intercalated lithium concentration from electrode material
        const double cmax = matelectrode->CMax();
        if(cmax < 1.e-12)
          dserror("Saturation value c_max of intercalated lithium concentration is too small!");

        // compute factor F/(RT)
        const double frt = GetFRT();

        // equilibrium electric potential difference at electrode surface
        const double epd = matelectrode->ComputeOpenCircuitPotential(eslavephiint,faraday,frt);

        // skip further computation in case equilibrium electric potential difference is outside physically meaningful range
        if(not std::isinf(epd))
        {
          // electrode-electrolyte overpotential at integration point
          const double eta = eslavepotint-emasterpotint-epd;

          // Butler-Volmer exchange mass flux density
          const double j0(kinmodel == INPAR::S2I::kinetics_butlervolmerreduced ? kr : kr*pow(emasterphiint,alphaa)*pow(cmax-eslavephiint,alphaa)*pow(eslavephiint,alphac));

          // exponential Butler-Volmer terms
          const double expterm1 = exp(alphaa*frt*eta);
          const double expterm2 = exp(-alphac*frt*eta);
          const double expterm = expterm1-expterm2;

          // safety check
          if(abs(expterm)>1.e5)
            dserror("Overflow of exponential term in Butler-Volmer formulation detected! Value: %lf",expterm);

          // core linearization associated with Butler-Volmer mass flux density
          const double dj_dd_slave = timefacwgt*j0*expterm;

          // loop over matrix columns
          for(int ui=0; ui<my::nen_; ++ui)
          {
            const int fui = ui*3;

            // loop over matrix rows
            for(int vi=0; vi<my::nen_; ++vi)
            {
              const int row_conc = vi*2;
              const int row_pot = row_conc+1;
              const double vi_dj_dd_slave = my::funct_(vi)*dj_dd_slave;

              // loop over spatial dimensions
              for(unsigned dim=0; dim<3; ++dim)
              {
                // compute linearizations w.r.t. slave-side structural displacements
                eslavematrix(row_conc,fui+dim) += vi_dj_dd_slave*shapederivatives(dim,ui);
                eslavematrix(row_pot,fui+dim) += nume*vi_dj_dd_slave*shapederivatives(dim,ui);
              }
            }
          }
        }

        break;
      }

      default:
      {
        dserror("Kinetic model for scatra-scatra interface coupling is not yet implemented!");
        break;
      }
    } // switch(kineticmodel)
  } // loop over integration points

  return;
} // DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::EvaluateS2ICouplingOD


/*-------------------------------------------------------------------------------------*
 | extract valence of species k from element material                       fang 02/15 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::GetValence(
    const Teuchos::RCP<const MAT::Material>&   material,   // element material
    const int                                  k           // species number
    ) const
{
  // valence cannot be computed for electrode material
  dserror("Valence cannot be computed for electrode material!");

  return 0.;
}


/*----------------------------------------------------------------------*
 | evaluate factor F/RT                                      fang 08/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::GetFRT() const
{
  // fetch factor F/RT from electrochemistry parameter list in isothermal case
  return myelch::elchparams_->FRT();
};


// template classes
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::quad4>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::tri3>;
// explicit instantiation of template methods
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::tri3>::EvaluateS2ICouplingAtIntegrationPoint<DRT::Element::tri3>(DRT::Condition&,const Teuchos::RCP<const MAT::Electrode>&,const std::vector<LINALG::Matrix<my::nen_,1> >&,const std::vector<LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement,1> >&,const LINALG::Matrix<my::nen_,1>&,const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement,1>&,const LINALG::Matrix<my::nen_,1>&,const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement,1>&,const double,const double,const double,Epetra_SerialDenseMatrix&,Epetra_SerialDenseMatrix&,Epetra_SerialDenseMatrix&,Epetra_SerialDenseMatrix&,Epetra_SerialDenseVector&,Epetra_SerialDenseVector&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::tri3>::EvaluateS2ICouplingAtIntegrationPoint<DRT::Element::quad4>(DRT::Condition&,const Teuchos::RCP<const MAT::Electrode>&,const std::vector<LINALG::Matrix<my::nen_,1> >&,const std::vector<LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement,1> >&,const LINALG::Matrix<my::nen_,1>&,const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement,1>&,const LINALG::Matrix<my::nen_,1>&,const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement,1>&,const double,const double,const double,Epetra_SerialDenseMatrix&,Epetra_SerialDenseMatrix&,Epetra_SerialDenseMatrix&,Epetra_SerialDenseMatrix&,Epetra_SerialDenseVector&,Epetra_SerialDenseVector&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::quad4>::EvaluateS2ICouplingAtIntegrationPoint<DRT::Element::tri3>(DRT::Condition&,const Teuchos::RCP<const MAT::Electrode>&,const std::vector<LINALG::Matrix<my::nen_,1> >&,const std::vector<LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement,1> >&,const LINALG::Matrix<my::nen_,1>&,const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement,1>&,const LINALG::Matrix<my::nen_,1>&,const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tri3>::numNodePerElement,1>&,const double,const double,const double,Epetra_SerialDenseMatrix&,Epetra_SerialDenseMatrix&,Epetra_SerialDenseMatrix&,Epetra_SerialDenseMatrix&,Epetra_SerialDenseVector&,Epetra_SerialDenseVector&);
template void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::quad4>::EvaluateS2ICouplingAtIntegrationPoint<DRT::Element::quad4>(DRT::Condition&,const Teuchos::RCP<const MAT::Electrode>&,const std::vector<LINALG::Matrix<my::nen_,1> >&,const std::vector<LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement,1> >&,const LINALG::Matrix<my::nen_,1>&,const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement,1>&,const LINALG::Matrix<my::nen_,1>&,const LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::quad4>::numNodePerElement,1>&,const double,const double,const double,Epetra_SerialDenseMatrix&,Epetra_SerialDenseMatrix&,Epetra_SerialDenseMatrix&,Epetra_SerialDenseMatrix&,Epetra_SerialDenseVector&,Epetra_SerialDenseVector&);
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::line3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::nurbs3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::nurbs9>;
