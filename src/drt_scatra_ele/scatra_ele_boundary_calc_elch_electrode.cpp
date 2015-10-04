/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_elch_electrode.cpp

\brief evaluation of ScaTra boundary elements for electrodes

<pre>
Maintainer: Rui Fang
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
    const DRT::FaceElement*     ele,              ///< current boundary element
    Teuchos::ParameterList&     params,           ///< parameter list
    DRT::Discretization&        discretization,   ///< discretization
    std::vector<int>&           lm,               ///< location vector
    Epetra_SerialDenseMatrix&   eslavematrix,     ///< element matrix for slave side
    Epetra_SerialDenseMatrix&   emastermatrix,    ///< element matrix for master side
    Epetra_SerialDenseVector&   eslaveresidual    ///< element residual for slave side
    )
{
  // access material of parent element
  Teuchos::RCP<const MAT::Electrode> matelectrode = Teuchos::rcp_dynamic_cast<const MAT::Electrode>(ele->ParentElement()->Material());

  // safety check
  if(matelectrode == Teuchos::null)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // get global and interface state vectors
  Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  Teuchos::RCP<const Epetra_Vector> imasterphinp = discretization.GetState("imasterphinp");
  if (phinp == Teuchos::null or imasterphinp == Teuchos::null)
    dserror("Cannot get state vector \"phinp\" or \"imasterphinp\"!");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  std::vector<LINALG::Matrix<my::nen_,1> > eslavephinp(my::numdofpernode_,LINALG::Matrix<my::nen_,1>(true));
  std::vector<LINALG::Matrix<my::nen_,1> > emasterphinp(my::numdofpernode_,LINALG::Matrix<my::nen_,1>(true));
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_,1> >(*phinp,eslavephinp,lm);
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_,1> >(*imasterphinp,emasterphinp,lm);

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
  const double frt = myelch::elchparams_->FRT();
  if(frt <= 0.)
    dserror("Factor F/RT is negative!");
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

      // evaluate dof values at current integration point on present and opposite side of scatra-scatra interface
      const double eslavephiint = my::funct_.Dot(eslavephinp[k]);
      const double eslavepotint = my::funct_.Dot(eslavephinp[my::numscal_]);
      const double emasterphiint = my::funct_.Dot(emasterphinp[k]);
      const double emasterpotint = my::funct_.Dot(emasterphinp[my::numscal_]);

      // equilibrium electric potential difference and its derivative w.r.t. concentration at electrode surface
      const double epd = matelectrode->ComputeOpenCircuitPotential(eslavephiint,faraday,frt);
      const double epdderiv = matelectrode->ComputeFirstDerivOpenCircuitPotential(eslavephiint,faraday,frt);

      // electrode-electrolyte overpotential at integration point
      const double eta = eslavepotint-emasterpotint-epd;

      // compute matrix and vector contributions according to kinetic model for current scatra-scatra interface coupling condition
      switch(kineticmodel)
      {
        // Butler-Volmer kinetics
        case INPAR::S2I::kinetics_butlervolmer:
        {
          const double i0 = kr*faraday*pow(emasterphiint,alphaa)*pow(cmax-eslavephiint,alphaa)*pow(eslavephiint,alphac);
          const double expterm1 = exp(alphaa*frt*eta);
          const double expterm2 = exp(-alphac*frt*eta);
          const double expterm = expterm1-expterm2;

          // safety check
          if(abs(expterm)>1.e5)
            dserror("Overflow of exponential term in Butler-Volmer formulation detected! Value: %lf",expterm);

          for (int vi=0; vi<my::nen_; ++vi)
          {
            const double funct_vi_fns_timefacfac = my::funct_(vi)*fns*timefacfac;
            const double funct_vi_fns_timefacfac_i0 = funct_vi_fns_timefacfac*i0;
            const int fvi = vi*my::numdofpernode_+k;

            for (int ui=0; ui<my::nen_; ++ui)
            {
              eslavematrix(fvi,ui*my::numdofpernode_+k) += funct_vi_fns_timefacfac*(kr*faraday*pow(emasterphiint,alphaa)*pow(cmax-eslavephiint,alphaa-1.)*pow(eslavephiint,alphac-1.)*(-alphaa*eslavephiint+alphac*(cmax-eslavephiint))*expterm+i0*(-alphaa*frt*epdderiv*expterm1-alphac*frt*epdderiv*expterm2))*my::funct_(ui);
              eslavematrix(fvi,ui*my::numdofpernode_+my::numscal_) += funct_vi_fns_timefacfac_i0*(alphaa*frt*expterm1+alphac*frt*expterm2)*my::funct_(ui);
              emastermatrix(fvi,ui*my::numdofpernode_+k) += funct_vi_fns_timefacfac*kr*faraday*alphaa*pow(emasterphiint,alphaa-1.)*pow(cmax-eslavephiint,alphaa)*pow(eslavephiint,alphac)*expterm*my::funct_(ui);
              emastermatrix(fvi,ui*my::numdofpernode_+my::numscal_) += funct_vi_fns_timefacfac_i0*(-alphaa*frt*expterm1-alphac*frt*expterm2)*my::funct_(ui);
            }

            eslaveresidual[fvi] -= my::funct_(vi)*fns*i0*expterm*timefacrhsfac;
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

    // compute matrix and rhs contributions arising from closing equation for electric potential
    switch(myelch::elchparams_->EquPot())
    {
    case INPAR::ELCH::equpot_enc:
    {
      // do nothing, since no boundary integral present
      break;
    }

    case INPAR::ELCH::equpot_enc_pde:
    case INPAR::ELCH::equpot_enc_pde_elim:
    case INPAR::ELCH::equpot_divi:
    {
      for (int vi=0; vi<my::nen_; ++vi)
      {
        for (int ui=0; ui<my::nen_; ++ui)
        {
          eslavematrix(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+k) += nume*eslavematrix(vi*my::numdofpernode_+k,ui*my::numdofpernode_+k);
          eslavematrix(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+my::numscal_) += nume*eslavematrix(vi*my::numdofpernode_+k,ui*my::numdofpernode_+my::numscal_);
          emastermatrix(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+k) += nume*emastermatrix(vi*my::numdofpernode_+k,ui*my::numdofpernode_+k);
          emastermatrix(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+my::numscal_) += nume*emastermatrix(vi*my::numdofpernode_+k,ui*my::numdofpernode_+my::numscal_);
        }

        eslaveresidual[vi*my::numdofpernode_+my::numscal_] += nume*eslaveresidual[vi*my::numdofpernode_+k];
      }

      break;
    }

    case INPAR::ELCH::equpot_laplace:
    {
      dserror("Laplace equation combined with scatra-scatra interface coupling not implemented!");
      break;
    }

    case INPAR::ELCH::equpot_poisson:
    {
      dserror("Poisson equation combined with scatra-scatra interface coupling not implemented!");
      break;
    }

    default:
    {
      dserror("Unknown closing equation for electric potential!");
      break;
    }
    } // switch(equpot_)
  } // loop over scalars

  return;
} // DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::EvaluateS2ICoupling


/*-------------------------------------------------------------------------------------*
 | extract valence of species k from element material                       fang 02/15 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
const double DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::GetValence(
    const Teuchos::RCP<const MAT::Material>&   material,   // element material
    const int                                  k           // species number
    ) const
{
  // valence cannot be computed for electrode material
  dserror("Valence cannot be computed for electrode material!");

  return 0.;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::quad4>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::line3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::nurbs3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<DRT::Element::nurbs9>;
