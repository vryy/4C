/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_sti_electrode.cpp

\brief evaluation of ScaTra boundary elements for heat transport within electrodes

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
 */
/*----------------------------------------------------------------------*/
#include "scatra_ele_boundary_calc_sti_electrode.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_mat/electrode.H"
#include "../drt_mat/soret.H"

#include "../drt_inpar/inpar_s2i.H"

#include "../drt_scatra_ele/scatra_ele_parameter_elch.H"
#include "../drt_scatra_ele/scatra_ele_parameter_timint.H"

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 11/15 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype>* DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype>::Instance(
    const int numdofpernode,
    const int numscal,
    const std::string& disname,
    const ScaTraEleBoundaryCalcSTIElectrode* delete_me
    )
{
  static std::map<std::string,ScaTraEleBoundaryCalcSTIElectrode<distype>* > instances;

  if(delete_me == NULL)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleBoundaryCalcSTIElectrode<distype>(numdofpernode,numscal,disname);
  }

  else
  {
    for(typename std::map<std::string,ScaTraEleBoundaryCalcSTIElectrode<distype>* >::iterator i=instances.begin(); i!=instances.end(); ++i)
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
 | singleton destruction                                     fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype>::Done()
{
  // delete singleton
  Instance(0,0,"",this);

  return;
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype>::ScaTraEleBoundaryCalcSTIElectrode(
    const int numdofpernode,
    const int numscal,
    const std::string& disname
    ) :
    // constructor of base class
    my::ScaTraEleBoundaryCalc(numdofpernode,numscal,disname)
{
  return;
}


/*----------------------------------------------------------------------------------------------------------------------------*
 | evaluate main-diagonal system matrix contributions associated with scatra-scatra interface coupling condition   fang 08/15 |
 *----------------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype>::EvaluateS2ICoupling(
    const DRT::FaceElement*        ele,              ///< current boundary element
    Teuchos::ParameterList&        params,           ///< parameter list
    DRT::Discretization&           discretization,   ///< discretization
    DRT::Element::LocationArray&   la,               ///< location array
    Epetra_SerialDenseMatrix&      eslavematrix,     ///< element matrix for slave side
    Epetra_SerialDenseVector&      eslaveresidual    ///< element residual for slave side
    )
{
  // safety check
  if(my::numscal_ != 1 or my::numdofpernode_ != 1)
    dserror("Invalid number of transported scalars or degrees of freedom per node!");

  // access primary and secondary materials of parent element
  Teuchos::RCP<const MAT::Soret> matsoret = Teuchos::rcp_dynamic_cast<const MAT::Soret>(ele->ParentElement()->Material());
  Teuchos::RCP<const MAT::Electrode> matelectrode = Teuchos::rcp_dynamic_cast<const MAT::Electrode>(ele->ParentElement()->Material(1));
  if(matsoret == Teuchos::null or matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // get global and interface state vectors
  Teuchos::RCP<const Epetra_Vector> thermo = discretization.GetState("phinp");
  Teuchos::RCP<const Epetra_Vector> scatra = discretization.GetState(2,"scatra");
  Teuchos::RCP<const Epetra_Vector> imasterscatra = discretization.GetState(2,"imasterscatra");
  if(thermo == Teuchos::null or scatra == Teuchos::null or imasterscatra == Teuchos::null)
    dserror("Cannot get state vector \"phinp\" or \"scatra\" or \"imasterscatra\"!");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  LINALG::Matrix<my::nen_,1> ethermo(true);
  std::vector<LINALG::Matrix<my::nen_,1> > eslavescatra(2,LINALG::Matrix<my::nen_,1>(true));
  std::vector<LINALG::Matrix<my::nen_,1> > emasterscatra(2,LINALG::Matrix<my::nen_,1>(true));
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_,1> >(*thermo,ethermo,la[0].lm_);
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_,1> >(*scatra,eslavescatra,la[2].lm_);
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_,1> >(*imasterscatra,emasterscatra,la[2].lm_);

  // get current scatra-scatra interface coupling condition
  Teuchos::RCP<DRT::Condition> s2icondition = params.get<Teuchos::RCP<DRT::Condition> >("condition");
  if(s2icondition == Teuchos::null)
    dserror("Cannot access scatra-scatra interface coupling condition!");

  // access input parameters associated with current condition
  const double faraday = INPAR::ELCH::faraday_const;
  const double alphaa = s2icondition->GetDouble("alpha_a");
  const double alphac = s2icondition->GetDouble("alpha_c");
  const double kr = s2icondition->GetDouble("k_r");
  if(kr < 0.)
    dserror("Charge transfer constant k_r is negative!");
  const double cmax = matelectrode->CMax();
  if(cmax < 1.e-12)
    dserror("Saturation value c_max of intercalated Lithium concentration is too small!");
  const double peltier = s2icondition->GetDouble("peltier");

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

    // evaluate dof values at current integration point on present and opposite side of scatra-scatra interface
    const double etempint = my::funct_.Dot(ethermo);
    const double eslavephiint = my::funct_.Dot(eslavescatra[0]);
    const double eslavepotint = my::funct_.Dot(eslavescatra[1]);
    const double emasterphiint = my::funct_.Dot(emasterscatra[0]);
    const double emasterpotint = my::funct_.Dot(emasterscatra[1]);

    // evaluate factor F/RT
    const double frt = faraday/(INPAR::ELCH::gas_const*etempint);

    // evaluate factor F/RT²
    const double frtt = frt/etempint;

    // equilibrium electric potential difference at electrode surface
    const double epd = matelectrode->ComputeOpenCircuitPotential(eslavephiint,faraday,frt);

    // electrode-electrolyte overpotential at integration point
    const double eta = eslavepotint-emasterpotint-epd;

    // compute matrix and vector contributions according to kinetic model for current scatra-scatra interface coupling condition
    switch(s2icondition->GetInt("kinetic model"))
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

        for(int vi=0; vi<my::nen_; ++vi)
        {
          for(int ui=0; ui<my::nen_; ++ui)
            eslavematrix(vi,ui) += timefacfac*i0*frtt*eta*(alphaa*expterm1+alphac*expterm2)*(eta+peltier)*my::funct_(vi)*my::funct_(ui);
          eslaveresidual[vi] += timefacrhsfac*my::funct_(vi)*i0*expterm*(eta+peltier);
        }

        break;
      }

      default:
      {
        dserror("Kinetic model for scatra-scatra interface coupling is not yet implemented!");
        break;
      }
    } // switch(s2icondition->GetInt("kinetic model"))
  } // loop over integration points

  return;
}


/*---------------------------------------------------------------------------------------------------------------------------*
 | evaluate off-diagonal system matrix contributions associated with scatra-scatra interface coupling condition   fang 08/15 |
 *---------------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype>::EvaluateS2ICouplingOD(
    const DRT::FaceElement*        ele,              ///< current boundary element
    Teuchos::ParameterList&        params,           ///< parameter list
    DRT::Discretization&           discretization,   ///< discretization
    DRT::Element::LocationArray&   la,               ///< location array
    Epetra_SerialDenseMatrix&      eslavematrix,     ///< element matrix for slave side
    Epetra_SerialDenseMatrix&      emastermatrix     ///< element matrix for master side
    )
{
  // safety check
  if(my::numscal_ != 1 or my::numdofpernode_ != 1)
    dserror("Invalid number of transported scalars or degrees of freedom per node!");

  // access primary and secondary materials of parent element
  Teuchos::RCP<const MAT::Soret> matsoret = Teuchos::rcp_dynamic_cast<const MAT::Soret>(ele->ParentElement()->Material());
  Teuchos::RCP<const MAT::Electrode> matelectrode = Teuchos::rcp_dynamic_cast<const MAT::Electrode>(ele->ParentElement()->Material(1));
  if(matsoret == Teuchos::null or matelectrode == Teuchos::null)
   dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // get global and interface state vectors
  Teuchos::RCP<const Epetra_Vector> thermo = discretization.GetState("phinp");
  Teuchos::RCP<const Epetra_Vector> scatra = discretization.GetState(2,"scatra");
  Teuchos::RCP<const Epetra_Vector> imasterscatra = discretization.GetState(2,"imasterscatra");
  if(thermo == Teuchos::null or scatra == Teuchos::null or imasterscatra == Teuchos::null)
    dserror("Cannot get state vector \"phinp\" or \"scatra\" or \"imasterscatra\"!");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  LINALG::Matrix<my::nen_,1> ethermo(true);
  std::vector<LINALG::Matrix<my::nen_,1> > eslavescatra(2,LINALG::Matrix<my::nen_,1>(true));
  std::vector<LINALG::Matrix<my::nen_,1> > emasterscatra(2,LINALG::Matrix<my::nen_,1>(true));
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_,1> >(*thermo,ethermo,la[0].lm_);
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_,1> >(*scatra,eslavescatra,la[2].lm_);
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_,1> >(*imasterscatra,emasterscatra,la[2].lm_);

  // get current scatra-scatra interface coupling condition
  Teuchos::RCP<DRT::Condition> s2icondition = params.get<Teuchos::RCP<DRT::Condition> >("condition");
  if(s2icondition == Teuchos::null)
   dserror("Cannot access scatra-scatra interface coupling condition!");

  // access input parameters associated with current condition
  const double faraday = INPAR::ELCH::faraday_const;
  const double alphaa = s2icondition->GetDouble("alpha_a");
  const double alphac = s2icondition->GetDouble("alpha_c");
  const double kr = s2icondition->GetDouble("k_r");
  if(kr < 0.)
   dserror("Charge transfer constant k_r is negative!");
  const double cmax = matelectrode->CMax();
  if(cmax < 1.e-12)
   dserror("Saturation value c_max of intercalated Lithium concentration is too small!");
  const double peltier = s2icondition->GetDouble("peltier");

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

    // evaluate dof values at current integration point on present and opposite side of scatra-scatra interface
    const double etempint = my::funct_.Dot(ethermo);
    const double eslavephiint = my::funct_.Dot(eslavescatra[0]);
    const double eslavepotint = my::funct_.Dot(eslavescatra[1]);
    const double emasterphiint = my::funct_.Dot(emasterscatra[0]);
    const double emasterpotint = my::funct_.Dot(emasterscatra[1]);

    // evaluate factor F/RT
    const double frt = faraday/(INPAR::ELCH::gas_const*etempint);

    // equilibrium electric potential difference and its derivative w.r.t. concentration at electrode surface
    const double epd = matelectrode->ComputeOpenCircuitPotential(eslavephiint,faraday,frt);
    const double epdderiv = matelectrode->ComputeFirstDerivOpenCircuitPotential(eslavephiint,faraday,frt);

    // electrode-electrolyte overpotential at integration point
    const double eta = eslavepotint-emasterpotint-epd;

    // compute derivatives of scatra-scatra interface coupling residuals w.r.t. concentration and electric potential according to kinetic model for current thermo-thermo interface coupling condition
    switch(s2icondition->GetInt("kinetic model"))
    {
      // Butler-Volmer-Peltier kinetics
      case INPAR::S2I::kinetics_butlervolmerpeltier:
      {
        const double i0 = kr*faraday*pow(emasterphiint,alphaa)*pow(cmax-eslavephiint,alphaa)*pow(eslavephiint,alphac);
        const double di0_dcs = kr*faraday*pow(emasterphiint,alphaa)*(alphac*pow(eslavephiint,alphac-1.)*pow(cmax-eslavephiint,alphaa)-alphaa*pow(eslavephiint,alphac)*pow(cmax-eslavephiint,alphaa-1.));
        const double di0_dce = kr*faraday*alphaa*pow(emasterphiint,alphaa-1.)*pow(eslavephiint,alphac)*pow(cmax-eslavephiint,alphaa);

        const double expterm1 = exp(alphaa*frt*eta);
        const double dexpterm1_dcs = -epdderiv*alphaa*frt*expterm1;

        const double expterm2 = exp(-alphac*frt*eta);
        const double dexpterm2_dcs = epdderiv*alphac*frt*expterm2;

        const double expterm = expterm1-expterm2;
        if(abs(expterm)>1.e5)
         dserror("Overflow of exponential term in Butler-Volmer formulation detected! Value: %lf",expterm);
        const double dexpterm_dcs = dexpterm1_dcs-dexpterm2_dcs;

        const double di_dcs = di0_dcs*expterm+i0*dexpterm_dcs;
        const double di_dce = di0_dce*expterm;
        const double di_dphis = i0*(alphaa*frt*expterm1+alphac*frt*expterm2);
        const double di_dphie = -1.*di_dphis;

        for (int vi=0; vi<my::nen_; ++vi)
          for (int ui=0; ui<my::nen_; ++ui)
          {
            // recurring indices
            const int colconc(ui*2);
            const int colpot(colconc+1);

            // linearizations w.r.t. concentration on slave side
            eslavematrix(vi,colconc) += -timefacfac*(di_dcs*(eta+peltier)-i0*expterm*epdderiv)*my::funct_(vi)*my::funct_(ui);

            // linearizations w.r.t. electric potential on slave side
            eslavematrix(vi,colpot) += -timefacfac*(di_dphis*(eta+peltier)+i0*expterm)*my::funct_(vi)*my::funct_(ui);

            // linearizations w.r.t. concentration on master side
            emastermatrix(vi,colconc) += -timefacfac*di_dce*(eta+peltier)*my::funct_(vi)*my::funct_(ui);

            // linearizations w.r.t. electric potential on master side
            emastermatrix(vi,colpot) += -timefacfac*(di_dphie*(eta+peltier)-i0*expterm)*my::funct_(vi)*my::funct_(ui);
          }

        break;
      }

      default:
      {
        dserror("Kinetic model for scatra-scatra interface coupling is not yet implemented!");
        break;
      }
    } // switch(s2icondition->GetInt("kinetic model"))
  } // loop over integration points

  return;
}


/*----------------------------------------------------------------------*
 | evaluate action                                           fang 11/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<distype>::EvaluateAction(
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
    case SCATRA::bd_calc_s2icoupling:
    {
      EvaluateS2ICoupling(ele,params,discretization,la,elemat1_epetra,elevec1_epetra);
      break;
    }

    case SCATRA::bd_calc_s2icoupling_od:
    {
      EvaluateS2ICouplingOD(ele,params,discretization,la,elemat1_epetra,elemat2_epetra);
      break;
    }

    default:
    {
      my::EvaluateAction(
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


// template classes
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::quad4>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::line3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::nurbs3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcSTIElectrode<DRT::Element::nurbs9>;
