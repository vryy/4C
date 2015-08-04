/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_elch_diffcond.cpp

\brief evaluation of ScaTra boundary elements for diffusion-conduction formulation

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
 */
/*----------------------------------------------------------------------*/
#include "../drt_inpar/inpar_s2i.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_mat/elchmat.H"
#include "../drt_mat/elchphase.H"
#include "../drt_mat/ion.H"
#include "../drt_mat/material.H"
#include "../drt_mat/newman.H"

#include "scatra_ele.H"
#include "scatra_ele_calc_elch_diffcond.H" // for diffusion manager
#include "scatra_ele_boundary_calc_elch_diffcond.H"

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 02/15 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype>* DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype>::Instance(
    const int numdofpernode,
    const int numscal,
    const std::string& disname,
    bool create
    )
{
  static std::map<std::string,ScaTraEleBoundaryCalcElchDiffCond<distype>* >  instances;

  if(create)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleBoundaryCalcElchDiffCond<distype>(numdofpernode,numscal,disname);
  }

  else if(instances.find(disname) != instances.end())
  {
    for( typename std::map<std::string,ScaTraEleBoundaryCalcElchDiffCond<distype>* >::iterator i=instances.begin(); i!=instances.end(); ++i )
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
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype>::Done()
{
  // delete singleton
  Instance(0,0,"",false);

  return;
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype>::ScaTraEleBoundaryCalcElchDiffCond(const int numdofpernode,const int numscal,const std::string& disname)
  : // constructor of base class
    myelectrode::ScaTraEleBoundaryCalcElchElectrode(numdofpernode,numscal,disname),
    // initialization of diffusion manager
    dmedc_(Teuchos::rcp(new ScaTraEleDiffManagerElchDiffCond(my::numscal_)))
{
  // replace standard electrochemistry parameter class by electrochemistry parameter class for diffusion-conduction formulation
  my::scatraparams_ = DRT::ELEMENTS::ScaTraEleParameterElchDiffCond::Instance(disname);

  return;
}


/*----------------------------------------------------------------------*
 | evaluate Neumann boundary condition                       fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype>::EvaluateNeumann(
    DRT::FaceElement*                   ele,
    Teuchos::ParameterList&             params,
    DRT::Discretization&                discretization,
    DRT::Condition&                     condition,
    DRT::Element::LocationArray&        la,
    Epetra_SerialDenseVector&           elevec1,
    const double                        scalar
    )
{
  // get material of parent element
  Teuchos::RCP<MAT::Material> mat = ele->ParentElement()->Material();

  if(mat->MaterialType() == INPAR::MAT::m_elchmat)
  {
    const MAT::ElchMat* actmat = static_cast<const MAT::ElchMat*>(mat.get());

    for (int iphase=0; iphase < actmat->NumPhase();++iphase)
    {
      const int phaseid = actmat->PhaseID(iphase);
      Teuchos::RCP<const MAT::Material> singlemat = actmat->PhaseById(phaseid);

      if(singlemat->MaterialType() == INPAR::MAT::m_elchphase)
      {
        const MAT::ElchPhase* actsinglemat = static_cast<const MAT::ElchPhase*>(singlemat.get());

        dmedc_->SetPhasePoro(actsinglemat->Epsilon(),iphase);
      }
      else
        dserror("Invalid material!");
    }
  }
  else
    dserror("Invalid material!");

  // call base class routine
  my::EvaluateNeumann(ele,params,discretization,condition,la,elevec1,dmedc_->GetPhasePoro(0));

  // add boundary flux contributions to potential equation
  switch(ElchParams()->EquPot())
  {
  case INPAR::ELCH::equpot_divi:
  {
    for(int k=0; k<my::numscal_; ++k)
    {
      // get valence
      const double valence_k = GetValence(mat,k);

      for(int vi=0; vi<my::nen_; ++vi)
        elevec1[vi*my::numdofpernode_+my::numscal_] += valence_k*elevec1[vi*my::numdofpernode_+k];;
    } // loop over scalars

    break;
  }

  default:
  {
    dserror("Closing equation for electric potential not recognized!");
    break;
  }
  }

  return 0;
}


/*----------------------------------------------------------------------*
 | evaluate electrode kinetics boundary condition            fang 02/15 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype>::CalcElchBoundaryKinetics(
    DRT::FaceElement*                 ele,
    Teuchos::ParameterList&           params,
    DRT::Discretization&              discretization,
    std::vector<int>&                 lm,
    Epetra_SerialDenseMatrix&         elemat1_epetra,
    Epetra_SerialDenseVector&         elevec1_epetra,
    const double                      scalar
    )
{
  Teuchos::RCP<MAT::Material> material = ele->ParentElement()->Material();

  if(material->MaterialType() == INPAR::MAT::m_elchmat)
  {
    const MAT::ElchMat* elchmat = static_cast<const MAT::ElchMat*>(material.get());

    for(int iphase=0; iphase<elchmat->NumPhase(); ++iphase)
    {
      Teuchos::RCP<const MAT::Material> phase = elchmat->PhaseById(elchmat->PhaseID(iphase));

      if(phase->MaterialType() == INPAR::MAT::m_elchphase)
        dmedc_->SetPhasePoro((static_cast<const MAT::ElchPhase*>(phase.get()))->Epsilon(),iphase);

      else
        dserror("Invalid material!");
    }
  }

  else
    dserror("Invalid material!");

  // call base class routine
  myelch::CalcElchBoundaryKinetics(ele,params,discretization,lm,elemat1_epetra,elevec1_epetra,dmedc_->GetPhasePoro(0));

  return;
}


/*----------------------------------------------------------------------*
 | evaluate electrode kinetics boundary condition            fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype>::EvaluateElchBoundaryKinetics(
    const DRT::Element*                 ele,        ///< the actual boundary element
    Epetra_SerialDenseMatrix&           emat,       ///< element-matrix
    Epetra_SerialDenseVector&           erhs,       ///< element-rhs
    const std::vector<double>&          ephinp,     ///< actual conc. and pot. values
    const std::vector<double>&          ehist,      ///< element history vector
    double                              timefac,    ///< time factor
    Teuchos::RCP<const MAT::Material>   material,   ///< the material
    Teuchos::RCP<DRT::Condition>        cond,       ///< the condition
    const int                           nume,       ///< number of transferred electrons
    const std::vector<int>              stoich,     ///< stoichiometry of the reaction
    const int                           kinetics,   ///< desired electrode kinetics model
    const double                        pot0,       ///< actual electrode potential on metal side
    const double                        frt,        ///< factor F/RT
    const double                        scalar      ///< scaling factor for element matrix and residual contributions
)
{
  // call base class routine
  myelch::EvaluateElchBoundaryKinetics(ele,emat,erhs,ephinp,ehist,timefac,material,cond,nume,stoich,kinetics,pot0,frt,scalar);

  // compute matrix and residual contributions arising from closing equation for electric potential
  switch(ElchParams()->EquPot())
  {
  case INPAR::ELCH::equpot_enc:
  {
    // do nothing, since no boundary integral present
    break;
  }

  case INPAR::ELCH::equpot_divi:
  {
    for(int k=0; k<my::numscal_; ++k)
    {
      for (int vi=0; vi<my::nen_; ++vi)
      {
        for (int ui=0; ui<my::nen_; ++ui)
        {
          emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+k) += nume*emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+k);
          emat(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+my::numscal_) += nume*emat(vi*my::numdofpernode_+k,ui*my::numdofpernode_+my::numscal_);
        }

        erhs[vi*my::numdofpernode_+my::numscal_] += nume*erhs[vi*my::numdofpernode_+k];
      }
    }

    break;
  }

  default:
  {
    dserror("Unknown closing equation for electric potential!");
    break;
  }
  } // switch(ElchParams()->EquPot())

  return;
} // DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype>::EvaluateElchBoundaryKinetics


/*-------------------------------------------------------------------------------------*
 | evaluate scatra-scatra interface coupling condition (electrochemistry)   fang 12/14 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype>::EvaluateS2ICoupling(
    const DRT::FaceElement*     ele,              ///< current boundary element
    Teuchos::ParameterList&     params,           ///< parameter list
    DRT::Discretization&        discretization,   ///< discretization
    std::vector<int>&           lm,               ///< location vector
    Epetra_SerialDenseMatrix&   eslavematrix,     ///< element matrix for slave side
    Epetra_SerialDenseMatrix&   emastermatrix,    ///< element matrix for master side
    Epetra_SerialDenseVector&   eslaveresidual    ///< element residual for slave side
    )
{
  // get material of parent element
  Teuchos::RCP<MAT::Material> material = ele->ParentElement()->Material();

  // safety check
  if(material->MaterialType() != INPAR::MAT::m_elchmat)
    dserror("Invalid electrolyte material for scatra-scatra interface coupling!");

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
  const double frt = myelch::ElchParams()->FRT();
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
    const double valence_k = GetValence(material,k);
    if(valence_k != nume)
      dserror("Number of transferred electrons must equal charge number of reacting species!");

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
      EquilibriumPotentialDifference(s2icondition,emasterphiint,epd,epdderiv);

      // electrode-electrolyte overpotential at integration point
      const double eta = emasterpotint-eslavepotint-epd;

      // compute matrix and vector contributions according to kinetic model for current scatra-scatra interface coupling condition
      switch(kineticmodel)
      {
        // Butler-Volmer kinetics
        case INPAR::S2I::kinetics_butlervolmer:
        {
          const double i0 = kr*faraday*pow(eslavephiint,alphaa)*pow(cmax-emasterphiint,alphaa)*pow(emasterphiint,alphac);
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
              eslavematrix(fvi,ui*my::numdofpernode_+k) -= funct_vi_fns_timefacfac*kr*faraday*alphaa*pow(eslavephiint,alphaa-1.)*pow(cmax-emasterphiint,alphaa)*pow(emasterphiint,alphac)*expterm*my::funct_(ui);
              eslavematrix(fvi,ui*my::numdofpernode_+my::numscal_) -= funct_vi_fns_timefacfac_i0*(-alphaa*frt*expterm1-alphac*frt*expterm2)*my::funct_(ui);
              emastermatrix(fvi,ui*my::numdofpernode_+k) -= funct_vi_fns_timefacfac*(kr*faraday*pow(eslavephiint,alphaa)*pow(cmax-emasterphiint,alphaa-1.)*pow(emasterphiint,alphac-1.)*(-alphaa*emasterphiint+alphac*(cmax-emasterphiint))*expterm+i0*(-alphaa*frt*epdderiv*expterm1-alphac*frt*epdderiv*expterm2))*my::funct_(ui);
              emastermatrix(fvi,ui*my::numdofpernode_+my::numscal_) -= funct_vi_fns_timefacfac_i0*(alphaa*frt*expterm1+alphac*frt*expterm2)*my::funct_(ui);
            }

            eslaveresidual[fvi] += my::funct_(vi)*fns*i0*expterm*timefacrhsfac;
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
    switch(ElchParams()->EquPot())
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
} // DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype>::EvaluateS2ICoupling


/*-------------------------------------------------------------------------------------*
 | extract valence of species k from element material                       fang 12/14 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
const double DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype>::GetValence(
    const Teuchos::RCP<const MAT::Material>&   material,   // element material
    const int                                  k           // species number
    ) const
{
  double valence(0.);

  if(material->MaterialType() == INPAR::MAT::m_elchmat)
  {
    const Teuchos::RCP<const MAT::ElchMat> elchmat = Teuchos::rcp_dynamic_cast<const MAT::ElchMat>(material);

    // safety check
    if(elchmat->NumPhase() != 1)
      dserror("Only one material phase is allowed at the moment!");

    // loop over phases
    for(int iphase=0; iphase<elchmat->NumPhase(); ++iphase)
    {
      const Teuchos::RCP<const MAT::ElchPhase> phase = Teuchos::rcp_dynamic_cast<const MAT::ElchPhase>(elchmat->PhaseById(elchmat->PhaseID(iphase)));

      // loop over species within phase
      for(int imat=0; imat<phase->NumMat(); ++imat)
      {
        const Teuchos::RCP<const MAT::Material> species = phase->MatById(phase->MatID(imat));

        if(species->MaterialType() == INPAR::MAT::m_newman)
        {
          valence = Teuchos::rcp_static_cast<const MAT::Newman>(species)->Valence();
          if(abs(valence) < 1.e-14)
            dserror("Received zero valence!");
        }
        else if(species->MaterialType() == INPAR::MAT::m_ion)
        {
          valence = Teuchos::rcp_static_cast<const MAT::Ion>(species)->Valence();
          if(abs(valence) < 1.e-14)
            dserror ("Received zero valence!");
        }
        else
          dserror("Unknown material species!");
      }
    }
  }

  else
    dserror("Unknown material!");

  return valence;
}


/*-------------------------------------------------------------------------------------------*
 | equilibrium electric potential difference at electrode-electrolyte interface   fang 01/15 |
 *-------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype>::EquilibriumPotentialDifference(
    const Teuchos::RCP<DRT::Condition>&   condition,       //! boundary condition
    const double&                         emasterphiint,   //! concentration of intercalated Lithium at electrode surface at current Gauss point
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
      const double frt = myelch::ElchParams()->FRT();
      const double cmax = condition->GetDouble("c_max");

      // intercalation fraction at electrode surface
      double X = emasterphiint/cmax;

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
}


// template classes
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<DRT::Element::quad4>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<DRT::Element::line3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<DRT::Element::nurbs3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<DRT::Element::nurbs9>;
