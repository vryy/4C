/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_elch_diffcond.cpp

\brief evaluation of ScaTra boundary elements for diffusion-conduction formulation

\level 2

<pre>
\maintainer Rui Fang
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
    const ScaTraEleBoundaryCalcElchDiffCond *delete_me
    )
{
  static std::map<std::string,ScaTraEleBoundaryCalcElchDiffCond<distype>* >  instances;

  if(delete_me == NULL)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleBoundaryCalcElchDiffCond<distype>(numdofpernode,numscal,disname);
  }

  else
  {
    for( typename std::map<std::string,ScaTraEleBoundaryCalcElchDiffCond<distype>* >::iterator i=instances.begin(); i!=instances.end(); ++i )
      if ( i->second == delete_me )
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
    dserror("Could not locate the desired instance. Internal error.");
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
  Instance(0,0,"",this);

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
 | evaluate action                                           fang 08/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype>::EvaluateAction(
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
  // extract location vector associated with primary dofset
  std::vector<int>& lm = la[0].lm_;

  // determine and evaluate action
  switch(action)
  {
  case SCATRA::bd_calc_elch_boundary_kinetics:
  {
    // access material of parent element
    Teuchos::RCP<MAT::Material> material = ele->ParentElement()->Material();

    // extract porosity from material and store in diffusion manager
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

    // process electrode kinetics boundary condition
    myelch::CalcElchBoundaryKinetics(
        ele,
        params,
        discretization,
        lm,
        elemat1_epetra,
        elevec1_epetra,
        dmedc_->GetPhasePoro(0)
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
  if(myelch::elchparams_->BoundaryFluxCoupling())
  {
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
  }

  return 0;
}


/*----------------------------------------------------------------------*
 | evaluate an electrode kinetics boundary condition         fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype>::EvaluateElchBoundaryKinetics(
    const DRT::Element*                               ele,        ///< current element
    Epetra_SerialDenseMatrix&                         emat,       ///< element matrix
    Epetra_SerialDenseVector&                         erhs,       ///< element right-hand side vector
    const std::vector<LINALG::Matrix<my::nen_,1> >&   ephinp,     ///< nodal values of concentration and electric potential
    const std::vector<LINALG::Matrix<my::nen_,1> >&   ehist,      ///< nodal history vector
    double                                            timefac,    ///< time factor
    Teuchos::RCP<const MAT::Material>                 material,   ///< material
    Teuchos::RCP<DRT::Condition>                      cond,       ///< electrode kinetics boundary condition
    const int                                         nume,       ///< number of transferred electrons
    const std::vector<int>                            stoich,     ///< stoichiometry of the reaction
    const int                                         kinetics,   ///< desired electrode kinetics model
    const double                                      pot0,       ///< electrode potential on metal side
    const double                                      frt,        ///< factor F/RT
    const double                                      scalar      ///< scaling factor for element matrix and right-hand side contributions
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
    const DRT::FaceElement*        ele,              ///< current boundary element
    Teuchos::ParameterList&        params,           ///< parameter list
    DRT::Discretization&           discretization,   ///< discretization
    DRT::Element::LocationArray&   la,               ///< location array
    Epetra_SerialDenseMatrix&      eslavematrix,     ///< element matrix for slave side
    Epetra_SerialDenseMatrix&      emastermatrix,    ///< element matrix for master side
    Epetra_SerialDenseVector&      eslaveresidual    ///< element residual for slave side
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
double DRT::ELEMENTS::ScaTraEleBoundaryCalcElchDiffCond<distype>::GetValence(
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
