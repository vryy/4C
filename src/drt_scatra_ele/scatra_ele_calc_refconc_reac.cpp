/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_refconc_reac.cpp

\brief main file containing routines for calculation of scatra element formulated in reference concentrations
       and with advanced reaction terms

\level 3

<pre>
  \maintainer Moritz Thon
              thon@mhpc.mw.tum.de
              http://www.mhpc.mw.tum.de
              089 - 289-10364
</pre>
 *----------------------------------------------------------------------*/


#include "scatra_ele_calc_refconc_reac.H"
#include "../drt_mat/matlist_reactions.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcRefConcReac<distype>::ScaTraEleCalcRefConcReac(const int numdofpernode,const int numscal,const std::string& disname)
: DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode,numscal,disname),
  DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::ScaTraEleCalcAdvReac(numdofpernode,numscal,disname),
  J_(1.0),
  C_inv_(true),
  dJdX_(true)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcRefConcReac<distype> * DRT::ELEMENTS::ScaTraEleCalcRefConcReac<distype>::Instance(
  const int numdofpernode,
  const int numscal,
  const std::string& disname,
  const ScaTraEleCalcRefConcReac* delete_me )
{
  static std::map<std::string,ScaTraEleCalcRefConcReac<distype>* >  instances;

  if(delete_me == NULL)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleCalcRefConcReac<distype>(numdofpernode,numscal,disname);
  }

  else
  {
    for( typename std::map<std::string,ScaTraEleCalcRefConcReac<distype>* >::iterator i=instances.begin(); i!=instances.end(); ++i )
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
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcRefConcReac<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( 0, 0, "", this );
}

//!
/*----------------------------------------------------------------------------*
 |  Set reac. body force, reaction coefficient and derivatives     thon 02/16 |
 *---------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcRefConcReac<distype>::SetAdvancedReactionTerms(
      const int                                 k,          //!< index of current scalar
      const Teuchos::RCP<MAT::MatListReactions> matreaclist //!< index of current scalar
  )
{
  const Teuchos::RCP<ScaTraEleReaManagerAdvReac> remanager = advreac::ReaManager();

  remanager->AddToReaBodyForce( matreaclist->CalcReaBodyForceTerm(k,my::scatravarmanager_->Phinp(),1.0/J_)*J_ ,k );

  matreaclist->CalcReaBodyForceDerivMatrix(k,remanager->GetReaBodyForceDerivVector(k),my::scatravarmanager_->Phinp());

}

/*------------------------------------------------------------------------------------------*
 |  calculation of convective element matrix: add conservative contributions     thon 02/16 |
 *------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcRefConcReac<distype>::CalcMatConvAddCons(
  Epetra_SerialDenseMatrix&     emat,
  const int                     k,
  const double                  timefacfac,
  const double                  vdiv,
  const double                  densnp
  )
{
  dserror("If you want to calculate the reference concentrations the CONVFORM must be 'convective'!");
}

/*------------------------------------------------------------------------------*
 | set internal variables                                           thon 02/16  |
 *------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcRefConcReac<distype>::SetInternalVariablesForMatAndRHS()
{
  //do the usual and...
  advreac::SetInternalVariablesForMatAndRHS();

  /////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////
  //spatial node coordinates
  LINALG::Matrix<my::nsd_,my::nen_> xyze(my::xyze_);
  xyze+=my::edispnp_;

  //! transposed jacobian "dx/ds"
  LINALG::Matrix<my::nsd_,my::nsd_> dxds(true);
  dxds.MultiplyNT(my::deriv_,xyze);

  // deformation gradtient dx/dX = dx/ds * ds/dX = dx/ds * (dX/ds)^(-1)
  LINALG::Matrix<my::nsd_,my::nsd_> F(true);
  F.MultiplyTT(dxds,my::xij_);

  // inverse of jacobian "dx/dX"
  LINALG::Matrix<my::nsd_,my::nsd_> F_inv(true);
  J_ = F_inv.Invert(F);

  //calculate inverse of cauchy-green stress tensor
  C_inv_.MultiplyNT(F_inv,F_inv);

  ////////////////////////////////////////////////////////////////////////////////////////////////
  //calculate derivative dJ/dX by finite differences
  ////////////////////////////////////////////////////////////////////////////////////////////////
  const double epsilon= 1.0e-8;

  for (int i=0;i<3;i++)
  {
    LINALG::Matrix<my::nsd_,my::nen_> xyze_epsilon(my::xyze_);
    for (int j=0; j<my::nen_; ++j)
      xyze_epsilon(i,j)=xyze_epsilon(i,j)+epsilon;

    LINALG::Matrix<my::nsd_,my::nsd_> xjm_epsilon(true);
    xjm_epsilon.MultiplyNT(my::deriv_,xyze_epsilon);

    LINALG::Matrix<my::nsd_,my::nsd_> xij_epsilon(true);
    xij_epsilon.Invert(xjm_epsilon);

    // dx/dX = dx/ds * ds/dX = dx/ds * (dX/ds)^(-1)
    LINALG::Matrix<my::nsd_,my::nsd_> F_epsilon(true);
    F_epsilon.MultiplyTT(dxds,xij_epsilon);

    // inverse of transposed jacobian "ds/dX"
    const double J_epsilon= F_epsilon.Determinant();
    const double dJdX_i = (J_epsilon-J_)/epsilon;

    dJdX_(i,0)= dJdX_i;
  }

  return;
}

/*------------------------------------------------------------------- *
 |  calculation of diffusive element matrix                thon 02/16 |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcRefConcReac<distype>::CalcMatDiff(
  Epetra_SerialDenseMatrix&     emat,
  const int                     k,
  const double                  timefacfac
  )
{
  LINALG::Matrix<my::nsd_,my::nsd_> Diff_tens(C_inv_);
  Diff_tens.Scale(my::diffmanager_->GetIsotropicDiff(k));

  for (int vi=0; vi<my::nen_; ++vi)
  {
    const int fvi = vi*my::numdofpernode_+k;

    for (int ui=0; ui<my::nen_; ++ui)
    {
      const int fui = ui*my::numdofpernode_+k;

      double laplawf = 0.0;
//      GetLaplacianWeakForm(laplawf,Diff_tens,vi,ui);
      for (int j = 0; j<my::nsd_; j++)
      {
        for (int i = 0; i<my::nsd_; i++)
        {
          laplawf += my::derxy_(j, vi)*Diff_tens(j,i)*my::derxy_(i, ui);
        }
      }

      emat(fvi,fui) += timefacfac*laplawf;
    }
  }



  LINALG::Matrix<my::nsd_,my::nsd_> Diff_tens2(C_inv_);
  Diff_tens2.Scale(my::diffmanager_->GetIsotropicDiff(k)/J_);

  for (int vi=0; vi<my::nen_; ++vi)
  {
    const int fvi = vi*my::numdofpernode_+k;

    double laplawf2 = 0.0;
    for (int j = 0; j<my::nsd_; j++)
    {
      for (int i = 0; i<my::nsd_; i++)
      {
        laplawf2 += my::derxy_(j,vi)*Diff_tens2(j,i)*dJdX_(i);
      }
    }

    for (int ui=0; ui<my::nen_; ++ui)
    {
      const int fui = ui*my::numdofpernode_+k;

      emat(fvi,fui) -= timefacfac*laplawf2*my::funct_(ui);
    }
  }

  return;
}

/*-------------------------------------------------------------------- *
 |  standard Galerkin diffusive term on right hand side     ehrl 11/13 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcRefConcReac<distype>::CalcRHSDiff(
  Epetra_SerialDenseVector&     erhs,
  const int                     k,
  const double                  rhsfac
  )
{
  /////////////////////////////////////////////////////////////////////
  // \D* \grad c_0 \times \grad \phi ...
  /////////////////////////////////////////////////////////////////////
  LINALG::Matrix<my::nsd_,my::nsd_> Diff_tens(C_inv_);
  Diff_tens.Scale(my::diffmanager_->GetIsotropicDiff(k));

  const LINALG::Matrix<my::nsd_,1>&  gradphi = my::scatravarmanager_->GradPhi(k);

  for (int vi=0; vi<my::nen_; ++vi)
  {
    const int fvi = vi*my::numdofpernode_+k;

    double laplawf(0.0);
//    GetLaplacianWeakFormRHS(laplawf,Diff_tens,gradphi,vi);
    for (int j = 0; j<my::nsd_; j++)
    {
      for (int i = 0; i<my::nsd_; i++)
      {
        laplawf += my::derxy_(j,vi)*Diff_tens(j,i)*gradphi(i);
      }
    }

    erhs[fvi] -= rhsfac*laplawf;
  }

  /////////////////////////////////////////////////////////////////////
  // ... + \D* c_0/J * \grad J \times \grad \phi
  /////////////////////////////////////////////////////////////////////
  LINALG::Matrix<my::nsd_,my::nsd_> Diff_tens2(C_inv_);
  Diff_tens2.Scale(my::diffmanager_->GetIsotropicDiff(k)/J_*my::scatravarmanager_->Phinp(k));

  for (int vi=0; vi<my::nen_; ++vi)
  {
    const int fvi = vi*my::numdofpernode_+k;

    double laplawf2(0.0);
//    GetLaplacianWeakFormRHS(laplawf2,Diff_tens2,dJdX_,vi);
    for (int j = 0; j<my::nsd_; j++)
    {
      for (int i = 0; i<my::nsd_; i++)
      {
        laplawf2 += my::derxy_(j,vi)*Diff_tens2(j,i)*dJdX_(i);
      }
    }

    erhs[fvi] += rhsfac*laplawf2;
  }

  return;
}




// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<DRT::Element::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<DRT::Element::quad9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<DRT::Element::nurbs9>;
//template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<DRT::Element::nurbs27>;
