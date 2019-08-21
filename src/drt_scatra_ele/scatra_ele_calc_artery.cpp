/*----------------------------------------------------------------------*/
/*! \file
 \brief evaluation class containing routines for calculation of scalar transport
        within 1D-arteries (blood vessels)
        only pressure-based formulation supports this

   \level 3

   \maintainer  Johannes Kremheller
 *----------------------------------------------------------------------*/

#include "scatra_ele_calc_artery.H"
#include "scatra_ele_parameter_std.H"
#include "scatra_ele_parameter_timint.H"
#include "../drt_mat/scatra_mat.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcArtery<distype, probdim>::ScaTraEleCalcArtery(
    const int numdofpernode, const int numscal, const std::string& disname)
    : my::ScaTraEleCalc(numdofpernode, numscal, disname)
{
  // safety check
  // Note: if higher order 1D elements should be used, the approach with adding the length from
  // below
  //        will not work anymore
  if (my::nen_ != 2)
    dserror(
        "Only line2 elements supported so far, you have %d nodes, if called with 2D or 3D element, "
        "think again",
        my::nen_);
  // replace internal variable manager by internal variable manager for artery
  my::scatravarmanager_ =
      Teuchos::rcp(new ScaTraEleInternalVariableManagerArtery<my::nsd_, my::nen_>(my::numscal_));
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcArtery<distype, probdim>*
DRT::ELEMENTS::ScaTraEleCalcArtery<distype, probdim>::Instance(const int numdofpernode,
    const int numscal, const std::string& disname, const ScaTraEleCalcArtery* delete_me)
{
  static std::map<std::string, ScaTraEleCalcArtery<distype, probdim>*> instances;

  if (delete_me == NULL)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] =
          new ScaTraEleCalcArtery<distype, probdim>(numdofpernode, numscal, disname);
  }

  else
  {
    for (typename std::map<std::string, ScaTraEleCalcArtery<distype, probdim>*>::iterator i =
             instances.begin();
         i != instances.end(); ++i)
      if (i->second == delete_me)
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
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcArtery<distype, probdim>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(0, 0, "", this);
}

/*----------------------------------------------------------------------*
 | setup element evaluation                            kremheller 03/18 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleCalcArtery<distype, probdim>::SetupCalc(
    DRT::Element* ele, DRT::Discretization& discretization)
{
  // base class
  my::SetupCalc(ele, discretization);

  // set the artery material in the variable manager
  VarManager()->SetArteryMaterial(ele);

  return 0;
}

/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)              kremheller 04/18 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcArtery<distype, probdim>::Materials(
    const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
    const int k,                                       //!< id of current scalar
    double& densn,                                     //!< density at t_(n)
    double& densnp,                                    //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,                                    //!< density at t_(n+alpha_M)
    double& visc,                                      //!< fluid viscosity
    const int iquad                                    //!< id of current gauss point

)
{
  switch (material->MaterialType())
  {
    case INPAR::MAT::m_scatra:
    {
      const Teuchos::RCP<const MAT::ScatraMat>& actmat =
          Teuchos::rcp_dynamic_cast<const MAT::ScatraMat>(material);

      densn = 1.0;
      densam = 1.0;
      densnp = 1.0;

      my::diffmanager_->SetIsotropicDiff(actmat->Diffusivity(), k);

      break;
    }
    default:
    {
      dserror("Material type %i is not supported!", material->MaterialType());
      break;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | set internal variables                              kremheller 03/18 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcArtery<distype, probdim>::SetInternalVariablesForMatAndRHS()
{
  VarManager()->SetInternalVariablesArtery(my::funct_, my::derxy_, my::deriv_, my::xjm_,
      my::ephinp_, my::ephin_, my::ehist_, earterypressurenp_);

  return;
}

/*----------------------------------------------------------------------*
 | extract element based or nodal values               kremheller 03/18 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcArtery<distype, probdim>::ExtractElementAndNodeValues(
    DRT::Element* ele, Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la)
{
  //---------------------------------------------------------------------------------------------
  //                                 SCATRA
  //---------------------------------------------------------------------------------------------

  // extract local values from the global vectors
  Teuchos::RCP<const Epetra_Vector> hist = discretization.GetState("hist");
  Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  if (hist == Teuchos::null || phinp == Teuchos::null)
    dserror("Cannot get state vector 'hist' and/or 'phinp'");

  // values of scatra field are always in first dofset
  const std::vector<int>& lm = la[0].lm_;
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_, 1>>(*hist, my::ehist_, lm);
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_, 1>>(*phinp, my::ephinp_, lm);

  if (my::scatraparatimint_->IsGenAlpha() and not my::scatraparatimint_->IsIncremental())
  {
    // extract additional local values from global vector
    Teuchos::RCP<const Epetra_Vector> phin = discretization.GetState("phin");
    if (phin == Teuchos::null) dserror("Cannot get state vector 'phin'");
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_, 1>>(*phin, my::ephin_, lm);
  }

  //---------------------------------------------------------------------------------------------
  //                                 CURRENT LENGTH
  //---------------------------------------------------------------------------------------------
  // extract element and node values of the artery
  if (discretization.HasState(1, "curr_seg_lengths"))
  {
    Teuchos::RCP<const Epetra_Vector> curr_seg_lengths =
        discretization.GetState(1, "curr_seg_lengths");
    std::vector<double> seglengths(la[1].lm_.size());

    DRT::UTILS::ExtractMyValues(*curr_seg_lengths, seglengths, la[1].lm_);

    const double curr_ele_length = std::accumulate(seglengths.begin(), seglengths.end(), 0.0);

    static LINALG::Matrix<probdim, 1> arteryrefpos0;
    for (unsigned int d = 0; d < probdim; ++d) arteryrefpos0(d) = my::xyze_(d, 0);
    static LINALG::Matrix<probdim, 1> arteryrefpos1;
    for (unsigned int d = 0; d < probdim; ++d) arteryrefpos1(d) = my::xyze_(d, 1);

    static LINALG::Matrix<probdim, 1> dist0;
    dist0.Update(-1.0, arteryrefpos0, 1.0, arteryrefpos1, 0.0);
    const double arteryreflength = dist0.Norm2();

    // this is a hack
    // will not work for anything else but line2 elements
    // change in length is simply added to displacement of second node
    for (unsigned int d = 0; d < probdim; ++d)
    {
      my::edispnp_(d, 0) = 0.0;
      my::edispnp_(d, 1) = (curr_ele_length / arteryreflength - 1.0) * dist0(d);
    }

    my::UpdateNodeCoordinates();
  }

  int ndsscatra_artery = 1;
  if (discretization.NumDofSets() == 3) ndsscatra_artery = 2;

  //---------------------------------------------------------------------------------------------
  //                                 ARTERY
  //---------------------------------------------------------------------------------------------
  // extract element and node values of the artery
  if (discretization.HasState(ndsscatra_artery, "one_d_artery_pressure"))
  {
    Teuchos::RCP<const Epetra_Vector> arterypn =
        discretization.GetState(ndsscatra_artery, "one_d_artery_pressure");
    // values of scatra field are always in first dofset
    const std::vector<int>& lm_artery = la[ndsscatra_artery].lm_;
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_, 1>>(
        *arterypn, earterypressurenp_, lm_artery);
  }
  else
    dserror("Something went wrong here, scatra-dis does not have artery primary variable");

  // ---------------------------------------------------------------------
  // call routine for calculation of body force in element nodes
  // (time n+alpha_F for generalized-alpha scheme, at time n+1 otherwise)
  // ---------------------------------------------------------------------
  my::BodyForce(ele);
  //--------------------------------------------------------------------------------
  // further node-based source terms not given via Neumann volume condition
  // i.e., set special body force for homogeneous isotropic turbulence
  //--------------------------------------------------------------------------------
  my::OtherNodeBasedSourceTerms(lm, discretization, params);
}

/*-----------------------------------------------------------------------------*
 |  calculation of convective element matrix                                    |
 |  in convective form (OD fluid)                              kremheller 05/18 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcArtery<distype, probdim>::CalcMatConvODFluid(
    Epetra_SerialDenseMatrix& emat, const int k, const int ndofpernodefluid,
    const double timefacfac, const double densnp, const LINALG::Matrix<my::nsd_, 1>& gradphi)
{
  const double prefac = timefacfac * VarManager()->Diam() * VarManager()->Diam() / 32.0 /
                        VarManager()->Visc() * (-1.0);
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    const int fvi = vi * my::numdofpernode_ + k;
    const double v = prefac * my::funct_(vi);


    for (unsigned ui = 0; ui < my::nen_; ++ui)
    {
      // get correct factor
      double laplawf(0.0);
      for (unsigned j = 0; j < my::nsd_; j++) laplawf += my::derxy_(j, ui) * gradphi(j);
      const int fui = ui;
      emat(fvi, fui) += v * laplawf;
    }
  }
  return;
}

// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcArtery<DRT::Element::line2, 1>;
template class DRT::ELEMENTS::ScaTraEleCalcArtery<DRT::Element::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcArtery<DRT::Element::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcArtery<DRT::Element::line3, 1>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::line3,2>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::line3,3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcArtery<DRT::Element::tri3, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcArtery<DRT::Element::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcArtery<DRT::Element::tri6, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcArtery<DRT::Element::quad4, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcArtery<DRT::Element::quad4, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcArtery<DRT::Element::quad9, 2>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::quad9,3>;
template class DRT::ELEMENTS::ScaTraEleCalcArtery<DRT::Element::nurbs9, 2>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::nurbs9,3>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcArtery<DRT::Element::hex8, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcArtery<DRT::Element::hex27, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcArtery<DRT::Element::tet4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcArtery<DRT::Element::tet10, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcArtery<DRT::Element::pyramid5, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcStd<DRT::Element::nurbs27>;s
