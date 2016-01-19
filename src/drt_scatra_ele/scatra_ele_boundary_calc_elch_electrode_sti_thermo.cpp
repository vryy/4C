/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_elch_electrode_sti_thermo.cpp

\brief evaluation of ScaTra boundary elements for thermodynamic electrodes

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
 */
/*----------------------------------------------------------------------*/
#include "scatra_ele_boundary_calc_elch_electrode_sti_thermo.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_mat/electrode.H"
#include "../drt_mat/soret.H"

#include "../drt_inpar/inpar_s2i.H"

#include "../drt_scatra_ele/scatra_ele_parameter_elch.H"
#include "../drt_scatra_ele/scatra_ele_parameter_timint.H"

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 08/15 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>* DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>::Instance(
    const int numdofpernode,
    const int numscal,
    const std::string& disname,
    const ScaTraEleBoundaryCalcElchElectrodeSTIThermo* delete_me
    )
{
  static std::map<std::string,ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>* > instances;

  if(delete_me == NULL)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>(numdofpernode,numscal,disname);
  }

  else
  {
    for(typename std::map<std::string,ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>* >::iterator i=instances.begin(); i!=instances.end(); ++i)
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
 | singleton destruction                                     fang 08/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>::Done()
{
  // delete singleton
  Instance(0,0,"",this);

  return;
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 08/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>::ScaTraEleBoundaryCalcElchElectrodeSTIThermo(
    const int numdofpernode,
    const int numscal,
    const std::string& disname
    ) :
    // constructor of base class
    myelectrode::ScaTraEleBoundaryCalcElchElectrode(numdofpernode,numscal,disname)
{
  return;
}


/*---------------------------------------------------------------------------------------------------------------------------*
 | evaluate off-diagonal system matrix contributions associated with scatra-scatra interface coupling condition   fang 08/15 |
 *---------------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>::EvaluateS2ICouplingOD(
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

  // get global and interface state vectors
  Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  Teuchos::RCP<const Epetra_Vector> imasterphinp = discretization.GetState("imasterphinp");
  if(phinp == Teuchos::null or imasterphinp == Teuchos::null)
    dserror("Cannot get state vector \"phinp\" or \"imasterphinp\"!");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  std::vector<LINALG::Matrix<my::nen_,1> > eslavephinp(my::numdofpernode_,LINALG::Matrix<my::nen_,1>(true));
  std::vector<LINALG::Matrix<my::nen_,1> > emasterphinp(my::numdofpernode_,LINALG::Matrix<my::nen_,1>(true));
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_,1> >(*phinp,eslavephinp,la[0].lm_);
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_,1> >(*imasterphinp,emasterphinp,la[0].lm_);

  // get current scatra-scatra interface coupling condition
  Teuchos::RCP<DRT::Condition> s2icondition = params.get<Teuchos::RCP<DRT::Condition> >("condition");
  if(s2icondition == Teuchos::null)
    dserror("Cannot access scatra-scatra interface coupling condition!");

  // access input parameters associated with current condition
  const int kineticmodel = s2icondition->GetInt("kinetic model");
  const int nume = s2icondition->GetInt("e-");
  if(not (nume > 0))
    dserror("Charge transfer at electrode-electrolyte interface must involve a positive number of electrons!");
  const std::vector<int>* stoichiometries = s2icondition->GetMutable<std::vector<int> >("stoichiometries");
  if(stoichiometries == NULL)
    dserror("Cannot access vector of stoichiometric coefficients for scatra-scatra interface coupling!");
  if(stoichiometries->size() != (unsigned) my::numscal_)
    dserror("Number of stoichiometric coefficients does not match number of scalars!");
  int reactivespecies = 0;
  for(int k=0; k<my::numscal_; ++k)
    reactivespecies += abs((*stoichiometries)[k]);
  if(reactivespecies > 1)
    dserror("Charge transfer at electrode-electrolyte interface must not involve more than one reactive species!");
  const double faraday = INPAR::ELCH::faraday_const;
  const double alphaa = s2icondition->GetDouble("alpha_a");
  const double alphac = s2icondition->GetDouble("alpha_c");
  const double kr = s2icondition->GetDouble("k_r");
  if(kr < 0.)
    dserror("Charge transfer constant k_r is negative!");
  const double cmax = matelectrode->CMax();
  if(cmax < 1.e-12)
    dserror("Saturation value c_max of intercalated Lithium concentration is too small!");

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over scalars
  for(int k=0; k<my::numscal_; ++k)
  {
    if(!(*stoichiometries)[k])
      continue;

    const double fns = -1./faraday/nume*(*stoichiometries)[k];

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
      const double frt = GetFRT(discretization,la);

      // evaluate factor F/RT²
      const double frtt = GetFRTT(discretization,la);

      // evaluate dof values at current integration point on present and opposite side of scatra-scatra interface
      const double eslavephiint = my::funct_.Dot(eslavephinp[k]);
      const double eslavepotint = my::funct_.Dot(eslavephinp[my::numscal_]);
      const double emasterphiint = my::funct_.Dot(emasterphinp[k]);
      const double emasterpotint = my::funct_.Dot(emasterphinp[my::numscal_]);

      // equilibrium electric potential difference at electrode surface
      const double epd = matelectrode->ComputeOpenCircuitPotential(eslavephiint,faraday,frt);

      // electrode-electrolyte overpotential at integration point
      const double eta = eslavepotint-emasterpotint-epd;

      // compute derivatives of scatra-scatra interface coupling residuals w.r.t. thermo dofs according to kinetic model for current scatra-scatra interface coupling condition
      switch(kineticmodel)
      {
        // Butler-Volmer-Peltier kinetics
        case INPAR::S2I::kinetics_butlervolmerpeltier:
        {
          const double i0 = kr*faraday*pow(emasterphiint,alphaa)*pow(cmax-eslavephiint,alphaa)*pow(eslavephiint,alphac);
          const double expterm1 = exp(alphaa*frt*eta);
          const double expterm2 = exp(-alphac*frt*eta);
          const double expterm = expterm1-expterm2;

          // safety check
          if(abs(expterm)>1.e5)
            dserror("Overflow of exponential term in Butler-Volmer formulation detected! Value: %lf",expterm);

          for (int vi=0; vi<my::nen_; ++vi)
            for (int ui=0; ui<my::nen_; ++ui)
              eslavematrix(vi*my::numdofpernode_+k,ui) -= my::funct_(vi)*fns*timefacfac*i0*frtt*eta*(alphaa*expterm1+alphac*expterm2)*my::funct_(ui);

          break;
        }

        default:
        {
          dserror("Kinetic model for scatra-scatra interface coupling is not yet implemented!");
          break;
        }
      } // switch(kineticmodel)
    } // loop over integration points

    // compute matrix contributions arising from closing equation for electric potential
    switch(myelch::elchparams_->EquPot())
    {
      case INPAR::ELCH::equpot_divi:
      {
        for (int vi=0; vi<my::nen_; ++vi)
          for (int ui=0; ui<my::nen_; ++ui)
            eslavematrix(vi*my::numdofpernode_+my::numscal_,ui) += nume*eslavematrix(vi*my::numdofpernode_+k,ui);

        break;
      }

      default:
      {
        dserror("Closing equation for electric potential not implemented for thermodynamic model!");
        break;
      }
    }  // switch(myelch::elchparams_->EquPot())
  } // loop over scalars

  return;
} // DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>::EvaluateS2ICouplingOD


/*----------------------------------------------------------------------*
 | evaluate action                                           fang 08/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>::EvaluateAction(
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
    case SCATRA::bd_calc_s2icoupling_od:
    {
      EvaluateS2ICouplingOD(ele,params,discretization,la,elemat1_epetra);
      break;
    }

    default:
    {
      myelectrode::EvaluateAction(
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


/*----------------------------------------------------------------------*
 | evaluate factor F/RT                                      fang 08/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>::GetFRT(
    DRT::Discretization&           discretization,   ///< discretization
    DRT::Element::LocationArray&   la                ///< location array
    ) const
{
  // evaluate factor F/RT
  return INPAR::ELCH::faraday_const/(INPAR::ELCH::gas_const*GetTemperature(discretization,la));
}


/*----------------------------------------------------------------------*
 | evaluate factor F/RT²                                     fang 08/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>::GetFRTT(
    DRT::Discretization&           discretization,   ///< discretization
    DRT::Element::LocationArray&   la                ///< location array
    ) const
{
  // evaluate local temperature value
  const double temperature = GetTemperature(discretization,la);

  // evaluate factor F/RT²
  return INPAR::ELCH::faraday_const/(INPAR::ELCH::gas_const*temperature*temperature);
}


/*----------------------------------------------------------------------*
 | evaluate local temperature value                          fang 08/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<distype>::GetTemperature(
    DRT::Discretization&           discretization,   ///< discretization
    DRT::Element::LocationArray&   la                ///< location array
    ) const
{
  // extract global state vector with thermo dofs from discretization
  Teuchos::RCP<const Epetra_Vector> tempnp = discretization.GetState(2,"thermo");
  if(tempnp == Teuchos::null)
    dserror("Cannot read state vector \"thermo\" from discretization!");

  // extract local nodal temperature values from global state vector
  LINALG::Matrix<my::nen_,1> etempnp(true);
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_,1> >(*tempnp,etempnp,la[2].lm_);

  // evaluate local temperature value
  const double temperature = my::funct_.Dot(etempnp);

  // safety check
  if(temperature <= 0.)
    dserror("Temperature is non-positive!");

  return temperature;
}


// template classes
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::quad4>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::line3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::nurbs3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrodeSTIThermo<DRT::Element::nurbs9>;
