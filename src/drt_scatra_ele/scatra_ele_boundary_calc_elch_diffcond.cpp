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
#include "scatra_ele_boundary_calc_elch_diffcond.H"
#include "scatra_ele_calc_elch_diffcond.H" // for diffusion manager

#include "../drt_mat/elchmat.H"
#include "../drt_mat/elchphase.H"
#include "../drt_mat/ion.H"
#include "../drt_mat/newman.H"

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
  switch(myelch::elchparams_->EquPot())
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
  switch(myelch::elchparams_->EquPot())
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
  } // switch(myelch::elchparams_->EquPot())

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
  // this function should never be called
  dserror("Each scatra-scatra interface for electrochemistry problems with conforming interface discretization "
          "must have an electrode on the slave side and the electrolyte on the master side, not the other way around!");

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
