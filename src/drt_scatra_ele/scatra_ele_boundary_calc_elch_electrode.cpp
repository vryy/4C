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

#include "../drt_mat/material.H"

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 02/15 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>* DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::Instance(
    const int numdofpernode,
    const int numscal,
    const std::string& disname,
    bool create
    )
{
  static std::map<std::string,ScaTraEleBoundaryCalcElchElectrode<distype>* >  instances;

  if(create)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleBoundaryCalcElchElectrode<distype>(numdofpernode,numscal,disname);
  }

  else if(instances.find(disname) != instances.end())
  {
    for( typename std::map<std::string,ScaTraEleBoundaryCalcElchElectrode<distype>* >::iterator i=instances.begin(); i!=instances.end(); ++i )
     {
      delete i->second;
      i->second = NULL;
     }

    instances.clear();
    return NULL;
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
  Instance(0,0,"",false);

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
  // check material of parent element
  if(ele->ParentElement()->Material()->MaterialType() != INPAR::MAT::m_electrode)
    dserror("Invalid electrode material for scatra-scatra interface coupling!");

  // get global and interface state vectors
  Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  Teuchos::RCP<const Epetra_Vector> imasterphinp = discretization.GetState("imasterphinp");
  if (phinp == Teuchos::null or imasterphinp == Teuchos::null)
    dserror("Cannot get state vector \"phinp\" or \"imasterphinp\"!");

  // extract local nodal values on present and opposite side of scatra-scatra interface
  std::vector<double> eslavephinpvec(lm.size());
  DRT::UTILS::ExtractMyValues(*phinp,eslavephinpvec,lm);
  std::vector<LINALG::Matrix<my::nen_,1> > eslavephinp(my::numscal_);
  LINALG::Matrix<my::nen_,1> eslavepotnp(true);
  std::vector<double> emasterphinpvec(lm.size());
  DRT::UTILS::ExtractMyValues(*imasterphinp,emasterphinpvec,lm);
  std::vector<LINALG::Matrix<my::nen_,1> > emasterphinp(my::numscal_);
  LINALG::Matrix<my::nen_,1> emasterpotnp(true);
  for(int inode=0; inode<my::nen_; ++inode)
  {
    for(int k=0; k<my::numscal_; ++k)
    {
      eslavephinp[k](inode,0) = eslavephinpvec[inode*my::numdofpernode_+k];
      emasterphinp[k](inode,0) = emasterphinpvec[inode*my::numdofpernode_+k];
    }

    eslavepotnp(inode,0) = eslavephinpvec[inode*my::numdofpernode_+my::numscal_];
    emasterpotnp(inode,0) = emasterphinpvec[inode*my::numdofpernode_+my::numscal_];
  }

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
  const double cmax = s2icondition->GetDouble("c_max");
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
      const double eslavepotint = my::funct_.Dot(eslavepotnp);
      const double emasterphiint = my::funct_.Dot(emasterphinp[k]);
      const double emasterpotint = my::funct_.Dot(emasterpotnp);

      // equilibrium electric potential difference and its derivative w.r.t. concentration at electrode surface
      double epd(0.);
      double epdderiv(0.);
      EquilibriumPotentialDifference(s2icondition,eslavephiint,epd,epdderiv);

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


/*-------------------------------------------------------------------------------------------*
 | equilibrium electric potential difference at electrode-electrolyte interface   fang 01/15 |
 *-------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::EquilibriumPotentialDifference(
    const Teuchos::RCP<DRT::Condition>&   condition,       //! boundary condition
    const double&                         eslavephiint,   //! concentration of intercalated Lithium at electrode surface at current Gauss point
    double&                               epd,             //! equilibrium electric potential difference to be computed at current Gauss point
    double&                               epdderiv         //! derivative of equilibrium electric potential difference to be computed at current Gauss point
    ) const
{
  // access model for equilibrium electric potential difference
  const int epdmodel = condition->GetInt("epd model");

  // compute corresponding equilibrium electric potential difference
  switch(epdmodel)
  {
    // Redlich-Kister expansion
    case INPAR::S2I::epd_redlichkister:
    {
      // extract relevant parameters from boundary condition
      const double DeltaG = condition->GetDouble("DeltaG");
      const int numcoeff = condition->GetInt("numcoeff");
      const std::vector<double>* coefficients = condition->GetMutable<std::vector<double> >("coefficients");
      if((int) coefficients->size() != numcoeff)
        dserror("Length of Redlich-Kister coefficient vector doesn't match prescribed number of coefficients!");
      const double faraday = INPAR::ELCH::faraday_const;
      const double frt = myelch::elchparams_->FRT();
      const double cmax = condition->GetDouble("c_max");

      // intercalation fraction at electrode surface
      double X = eslavephiint/cmax;

      // need to avoid intercalation fraction of exactly 0.5 due to singularity in Redlich-Kister expansion
      if(X == 0.5)
        X = 0.499999;

      // equilibrium electric potential difference according to Redlich-Kister expansion
      epd = DeltaG + faraday/frt*log((1.-X)/X);
      for(int i=0; i<numcoeff; ++i)
        epd += (*coefficients)[i]*(pow(2.*X-1.,i+1)-2.*i*X*(1.-X)*pow(2.*X-1.,i-1));
      epd /= faraday;

      // derivative of equilibrium electric potential difference w.r.t. concentration at electrode surface
      epdderiv = faraday/(2.*frt*X*(X-1.));
      for(int i=0; i<numcoeff; ++i)
        epdderiv += (*coefficients)[i]*((2.*i+1.)*pow(2.*X-1.,i)+2.*X*i*(X-1.)*(i-1.)*pow(2.*X-1.,i-2));
      epdderiv *= 2./(faraday*cmax);

      break;
    }

    default:
    {
      dserror("Unknown model for equilibrium electric potential difference!");
      break;
    }
  }

  return;
} // DRT::ELEMENTS::ScaTraEleBoundaryCalcElchElectrode<distype>::EquilibriumPotentialDifference


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
