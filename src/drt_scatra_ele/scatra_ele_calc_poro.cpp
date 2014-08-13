/*----------------------------------------------------------------------*/
/*!
 \file scatra_ele_calc_poro.cpp

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/


#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"

#include "../drt_mat/structporo.H"
#include "../drt_mat/scatra_mat.H"

//#include "scatra_ele_parameter.H"
//
//#include "../drt_nurbs_discret/drt_nurbs_utils.H"
//#include "../drt_geometry/position_array.H"
//
//#include "scatra_ele_action.H"
//#include "scatra_ele.H"

#include "scatra_ele_calc_poro.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcPoro<distype> * DRT::ELEMENTS::ScaTraEleCalcPoro<distype>::Instance(
  const int numdofpernode,
  const int numscal,
  bool create )
{
  static ScaTraEleCalcPoro<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new ScaTraEleCalcPoro<distype>(numdofpernode,numscal);
    }
  }
  else
  {
    if ( instance!=NULL )
      delete instance;
    instance = NULL;
  }
  return instance;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcPoro<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( 0, 0, false );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcPoro<distype>::ScaTraEleCalcPoro(const int numdofpernode,const int numscal)
  : DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode,numscal)
{

}

///*----------------------------------------------------------------------*
// * Action type: Evaluate
// *----------------------------------------------------------------------*/
//template <DRT::Element::DiscretizationType distype>
//int DRT::ELEMENTS::ScaTraEleCalcPoro<distype>::Evaluate(
//  DRT::ELEMENTS::Transport*  ele,
//  Teuchos::ParameterList&    params,
//  DRT::Discretization&       discretization,
//  const std::vector<int>&    lm,
//  Epetra_SerialDenseMatrix&  elemat1_epetra,
//  Epetra_SerialDenseMatrix&  elemat2_epetra,
//  Epetra_SerialDenseVector&  elevec1_epetra,
//  Epetra_SerialDenseVector&  elevec2_epetra,
//  Epetra_SerialDenseVector&  elevec3_epetra
//  )
//{
//  // check for the action parameter
//  const SCATRA::Action action = DRT::INPUT::get<SCATRA::Action>(params,"action");
//  switch(action)
//  {
//    case SCATRA::calc_scatra_mono_odblock_mesh:
//    {
//      return EvaluateODMesh(
//          ele,
//          params,
//          discretization,
//          lm,
//          elemat1_epetra,
//          elemat2_epetra,
//          elevec1_epetra,
//          elevec2_epetra,
//          elevec3_epetra
//          );
//      break;
//    }
//    case SCATRA::calc_scatra_mono_odblock_fluid:
//    {
//      return EvaluateODFluid(
//          ele,
//          params,
//          discretization,
//          lm,
//          elemat1_epetra,
//          elemat2_epetra,
//          elevec1_epetra,
//          elevec2_epetra,
//          elevec3_epetra
//          );
//      break;
//    }
//    default:
//    {
//      return my::Evaluate(
//          ele,
//          params,
//          discretization,
//          lm,
//          elemat1_epetra,
//          elemat2_epetra,
//          elevec1_epetra,
//          elevec2_epetra,
//          elevec3_epetra
//          );
//      break;
//    }
//  }
//
//  //you should no turn up here -> return error code
//  return -1;
//}

/*----------------------------------------------------------------------*
 |  Material ScaTra                                           |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcPoro<distype>::MatScaTra(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
  Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
  double&                                 visc,     //!< fluid viscosity
  const int                               iquad     //!< id of current gauss point
  )
{
  if(iquad==-1)
    dserror("no gauss point given for evaluation of scatra material. Check your input file.");

  const double porosity = GetPorosityAtGP(iquad);

  const Teuchos::RCP<const MAT::ScatraMat>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::ScatraMat>(material);

  // set diffusivity (scaled with porosity)
  SetDiffusivity(actmat,k,diffmanager,porosity);

  // set reaction coefficient
  SetReaCoefficient(actmat,k,reamanager,porosity);

  // set densities (scaled with porosity)
  SetDensities(porosity,densn,densnp,densam);

  return;
} // ScaTraEleCalcPoro<distype>::MatScaTra

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleCalcPoro<distype>::GetPorosityAtGP(const int  iquad)
{
  //access structure discretization
  Teuchos::RCP<DRT::Discretization> structdis = Teuchos::null;
  structdis = DRT::Problem::Instance()->GetDis("structure");
  //get corresponding fluid element (it has the same global ID as the scatra element)
  DRT::Element* structele = structdis->gElement(my::eid_);
  if (structele == NULL)
    dserror("Structure element %i not on local processor", my::eid_);

  const Teuchos::RCP<const MAT::StructPoro>& structmat
            = Teuchos::rcp_dynamic_cast<const MAT::StructPoro>(structele->Material());
  if(structmat == Teuchos::null)
    dserror("invalid structure material for poroelasticity");

  return structmat->GetPorosityAtGP(iquad);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
inline void DRT::ELEMENTS::ScaTraEleCalcPoro<distype>::SetDiffusivity(
    const Teuchos::RCP<const MAT::ScatraMat>& material,
    const int                                 k,
    Teuchos::RCP<ScaTraEleDiffManager>        diffmanager,
    const double                              porosity)
{
  diffmanager->SetIsotropicDiff(material->Diffusivity()*porosity,k);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
inline void DRT::ELEMENTS::ScaTraEleCalcPoro<distype>::SetReaCoefficient(
    const Teuchos::RCP<const MAT::ScatraMat>& material,
    const int                                k,
    Teuchos::RCP<ScaTraEleReaManager>        reamanager,
    const double                             porosity)
{
  //set reaction coefficient (no scaling with porosity)
  reamanager->SetReaCoeff(material->ReaCoeff(),k);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
inline void DRT::ELEMENTS::ScaTraEleCalcPoro<distype>::SetDensities(
    double  porosity,
    double& densn,
    double& densnp,
    double& densam
    )
{
  densn = porosity;
  densnp = porosity;
  densam = porosity;

  return;
}


// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::line3>;

// 2D elements
//template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::quad9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::wedge6>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::nurbs9>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoro<DRT::Element::nurbs27>;

