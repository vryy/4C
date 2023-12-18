/*----------------------------------------------------------------------*/
/*! \file

\brief Internal implementation of PressureBased artery element


\level 3
*/
/*----------------------------------------------------------------------*/

#include "baci_art_net_artery_ele_calc_pres_based.H"

#include "baci_art_net_art_junction.H"
#include "baci_art_net_art_terminal_bc.H"
#include "baci_art_net_artery_ele_calc.H"
#include "baci_discretization_fem_general_utils_fem_shapefunctions.H"
#include "baci_lib_globalproblem.H"
#include "baci_lib_utils.H"
#include "baci_mat_cnst_1d_art.H"
#include "baci_utils_function.H"
#include "baci_utils_singleton_owner.H"

#include <fstream>
#include <iomanip>

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::ArteryEleCalcPresBased<distype>::ArteryEleCalcPresBased(
    const int numdofpernode, const std::string& disname)
    : DRT::ELEMENTS::ArteryEleCalc<distype>(numdofpernode, disname)
{
}

/*----------------------------------------------------------------------*
 | singleton access method                                   vuong 08/16 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::ArteryEleCalcPresBased<distype>*
DRT::ELEMENTS::ArteryEleCalcPresBased<distype>::Instance(
    const int numdofpernode, const std::string& disname)
{
  using Key = std::pair<std::string, int>;
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<Key>(
      [](const int numdofpernode, const std::string& disname)
      {
        return std::unique_ptr<ArteryEleCalcPresBased<distype>>(
            new ArteryEleCalcPresBased<distype>(numdofpernode, disname));
      });

  std::pair<std::string, int> key(disname, numdofpernode);

  return singleton_map[key].Instance(CORE::UTILS::SingletonAction::create, numdofpernode, disname);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
int DRT::ELEMENTS::ArteryEleCalcPresBased<distype>::Evaluate(Artery* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra, Teuchos::RCP<MAT::Material> mat)
{
  // the number of nodes
  const int numnode = my::iel_;

  // construct views
  CORE::LINALG::Matrix<numnode, numnode> elemat1(elemat1_epetra.values(), true);
  CORE::LINALG::Matrix<numnode, 1> elevec1(elevec1_epetra.values(), true);

  // ---------------------------------------------------------------------
  // call routine for calculating element matrix and right hand side
  // ---------------------------------------------------------------------
  Sysmat(ele, discretization, la, elemat1, elevec1, mat);

  return 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
int DRT::ELEMENTS::ArteryEleCalcPresBased<distype>::EvaluateService(Artery* ele,
    const ARTERY::Action action, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra, Teuchos::RCP<MAT::Material> mat)
{
  switch (action)
  {
    case ARTERY::calc_flow_pressurebased:
      EvaluateFlow(ele, discretization, la, elevec1_epetra, mat);
      break;
    default:
      dserror("Unkown type of action %d for Artery (PressureBased formulation)", action);
  }

  return 0;
}

template <CORE::FE::CellType distype>
int DRT::ELEMENTS::ArteryEleCalcPresBased<distype>::ScatraEvaluate(Artery* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra, Teuchos::RCP<MAT::Material> mat)
{
  dserror(
      "not implemented by pressure-based formulation, should be done by cloned "
      "ScaTra-Discretization");

  return 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ArteryEleCalcPresBased<distype>::Sysmat(Artery* ele,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    CORE::LINALG::Matrix<my::iel_, my::iel_>& sysmat, CORE::LINALG::Matrix<my::iel_, 1>& rhs,
    Teuchos::RCP<const MAT::Material> material)
{
  // clear
  rhs.Clear();
  sysmat.Clear();

  // set element data
  const int numnode = my::iel_;

  // get pressure
  Teuchos::RCP<const Epetra_Vector> pressnp = discretization.GetState(0, "pressurenp");
  if (pressnp == Teuchos::null) dserror("could not get pressure inside artery element");

  // extract local values of pressure field from global state vector
  CORE::LINALG::Matrix<my::iel_, 1> mypress(true);
  DRT::UTILS::ExtractMyValues<CORE::LINALG::Matrix<my::iel_, 1>>(*pressnp, mypress, la[0].lm_);

  // calculate the element length
  const double L = CalculateEleLength(ele, discretization, la);

  // check here, if we really have an artery !!
  if (material->MaterialType() != INPAR::MAT::m_cnst_art) dserror("Wrong material type for artery");

  // cast the material to artery material material
  const MAT::Cnst_1d_art* actmat = static_cast<const MAT::Cnst_1d_art*>(material.get());

  if (actmat->IsCollapsed()) return;

  // Read in diameter
  const double diam = actmat->Diam();
  // Read in blood viscosity
  const double visc = actmat->Viscosity();

  const double hag_pois = M_PI * pow(diam, 4) / 128.0 / visc;
  // gaussian points
  const CORE::FE::IntegrationPoints1D intpoints(ele->GaussRule());

  // get Jacobian matrix and determinant
  // actually compute its transpose....
  /*
                                        _____________________________________
      ds     L      dxi    2           /         2            2            2
      --- = ---   ; --- = ---   ; L = / ( x - x )  + ( y - y )  + ( z - z )
      dxi    2      ds     L         v     1   2        2    2       1   2

  */
  my::xji_ = 2.0 / L;

  const double prefac = hag_pois * my::xji_(0, 0);

  // integration loop
  for (int iquad = 0; iquad < intpoints.nquad; ++iquad)
  {
    // coordinates of the current integration point
    const double xi = intpoints.qxg[iquad][0];
    const double wgt = intpoints.qwgt[iquad];

    const double fac = prefac * wgt;

    // shape functions and their derivatives
    CORE::FE::shape_function_1D_deriv1(my::deriv_, xi, distype);

    for (int inode = 0; inode < numnode; inode++)
      for (int jnode = 0; jnode < numnode; jnode++)
        sysmat(inode, jnode) += my::deriv_(0, inode) * fac * my::deriv_(0, jnode);

    // note: incremental form since rhs-coupling with poromultielastscatra-framework might be
    //       nonlinear
    CORE::LINALG::Matrix<1, 1> pressgrad;
    pressgrad.Multiply(my::deriv_, mypress);
    for (int inode = 0; inode < numnode; inode++)
      rhs(inode) -= my::deriv_(0, inode) * fac * pressgrad(0, 0);
  }


  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ArteryEleCalcPresBased<distype>::EvaluateFlow(Artery* ele,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    CORE::LINALG::SerialDenseVector& flowVec, Teuchos::RCP<const MAT::Material> material)
{
  // get pressure
  Teuchos::RCP<const Epetra_Vector> pressnp = discretization.GetState(0, "pressurenp");
  if (pressnp == Teuchos::null) dserror("could not get pressure inside artery element");

  // extract local values of pressure field from global state vector
  CORE::LINALG::Matrix<my::iel_, 1> mypress(true);
  DRT::UTILS::ExtractMyValues<CORE::LINALG::Matrix<my::iel_, 1>>(*pressnp, mypress, la[0].lm_);

  // calculate the element length
  const double L = CalculateEleLength(ele, discretization, la);

  // check here, if we really have an artery !!
  if (material->MaterialType() != INPAR::MAT::m_cnst_art) dserror("Wrong material type for artery");

  // cast the material to artery material material
  const MAT::Cnst_1d_art* actmat = static_cast<const MAT::Cnst_1d_art*>(material.get());

  // Read in diameter
  const double diam = actmat->Diam();
  // Read in blood viscosity
  const double visc = actmat->Viscosity();

  const double hag_pois = M_PI * pow(diam, 4) / 128.0 / visc;

  // TODO: this works only for line 2 elements
  flowVec(0) = -hag_pois * (mypress(1) - mypress(0)) / L;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
double DRT::ELEMENTS::ArteryEleCalcPresBased<distype>::CalculateEleLength(
    Artery* ele, DRT::Discretization& discretization, DRT::Element::LocationArray& la)
{
  double length;
  // get current element length
  if (discretization.NumDofSets() > 1 && discretization.HasState(1, "curr_seg_lengths"))
  {
    Teuchos::RCP<const Epetra_Vector> curr_seg_lengths =
        discretization.GetState(1, "curr_seg_lengths");
    std::vector<double> seglengths(la[1].lm_.size());

    DRT::UTILS::ExtractMyValues(*curr_seg_lengths, seglengths, la[1].lm_);

    length = std::accumulate(seglengths.begin(), seglengths.end(), 0.0);
  }
  else
    length = my::CalculateEleLength(ele);

  return length;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes

// 1D elements
template class DRT::ELEMENTS::ArteryEleCalcPresBased<CORE::FE::CellType::line2>;

BACI_NAMESPACE_CLOSE
