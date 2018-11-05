/*----------------------------------------------------------------------*/
/*!
\file artery_ele_calc_pres_based.cpp

\brief Internal implementation of PressureBased artery element

\maintainer Johannes Kremheller
            kremheller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/

\level 3
*/
/*----------------------------------------------------------------------*/

#include "artery_ele_calc_pres_based.H"


#include "artery_ele_calc.H"

#include "../drt_mat/cnst_1d_art.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_globalproblem.H"
#include "art_junction.H"
#include "art_terminal_bc.H"
#include <fstream>
#include <iomanip>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ArteryEleCalcPresBased<distype>::ArteryEleCalcPresBased(
    const int numdofpernode, const std::string& disname)
    : DRT::ELEMENTS::ArteryEleCalc<distype>(numdofpernode, disname)
{
}

/*----------------------------------------------------------------------*
 | singleton access method                                   vuong 08/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ArteryEleCalcPresBased<distype>*
DRT::ELEMENTS::ArteryEleCalcPresBased<distype>::Instance(
    const int numdofpernode, const std::string& disname, const ArteryEleCalcPresBased* delete_me)
{
  static std::map<std::pair<std::string, int>, ArteryEleCalcPresBased<distype>*> instances;

  std::pair<std::string, int> key(disname, numdofpernode);

  if (delete_me == NULL)
  {
    if (instances.find(key) == instances.end())
      instances[key] = new ArteryEleCalcPresBased<distype>(numdofpernode, disname);
  }

  else
  {
    // since we keep several instances around in the general case, we need to
    // find which of the instances to delete with this call. This is done by
    // letting the object to be deleted hand over the 'this' pointer, which is
    // located in the map and deleted
    for (typename std::map<std::pair<std::string, int>, ArteryEleCalcPresBased<distype>*>::iterator
             i = instances.begin();
         i != instances.end(); ++i)
      if (i->second == delete_me)
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
    dserror("Could not locate the desired instance. Internal error.");
  }

  return instances[key];
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ArteryEleCalcPresBased<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(0, "", this);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ArteryEleCalcPresBased<distype>::Evaluate(Artery* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& elemat1_epetra,
    Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
    Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
    Teuchos::RCP<MAT::Material> mat)
{
  // the number of nodes
  const int numnode = my::iel_;

  // construct views
  LINALG::Matrix<numnode, numnode> elemat1(elemat1_epetra.A(), true);
  LINALG::Matrix<numnode, 1> elevec1(elevec1_epetra.A(), true);

  // ---------------------------------------------------------------------
  // call routine for calculating element matrix and right hand side
  // ---------------------------------------------------------------------
  Sysmat(ele, discretization, la, elemat1, elevec1, mat);

  return 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ArteryEleCalcPresBased<distype>::EvaluateService(Artery* ele,
    const ARTERY::Action action, Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra, Teuchos::RCP<MAT::Material> mat)
{
  switch (action)
  {
    default:
      dserror("Unkown type of action %d for Artery (PressureBased formulation)", action);
  }

  return 0;
}

template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ArteryEleCalcPresBased<distype>::ScatraEvaluate(Artery* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization, std::vector<int>& lm,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra, Teuchos::RCP<MAT::Material> mat)
{
  dserror(
      "not implemented by pressure-based formulation, should be done by cloned "
      "ScaTra-Discretization");

  return 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ArteryEleCalcPresBased<distype>::Sysmat(Artery* ele,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    LINALG::Matrix<my::iel_, my::iel_>& sysmat, LINALG::Matrix<my::iel_, 1>& rhs,
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
  LINALG::Matrix<my::iel_, 1> mypress(true);
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::iel_, 1>>(*pressnp, mypress, la[0].lm_);

  double L;
  // get current element length
  if (discretization.NumDofSets() > 1 && discretization.HasState(1, "curr_seg_lengths"))
  {
    Teuchos::RCP<const Epetra_Vector> curr_seg_lengths =
        discretization.GetState(1, "curr_seg_lengths");
    std::vector<double> seglengths(la[1].lm_.size());

    DRT::UTILS::ExtractMyValues(*curr_seg_lengths, seglengths, la[1].lm_);

    L = std::accumulate(seglengths.begin(), seglengths.end(), 0.0);
  }
  else
    L = my::CalculateEleLength(ele);

  // check here, if we really have an artery !!
  if (material->MaterialType() != INPAR::MAT::m_cnst_art) dserror("Wrong material type for artery");

  // cast the material to artery material material
  const MAT::Cnst_1d_art* actmat = static_cast<const MAT::Cnst_1d_art*>(material.get());

  // Read in blood viscosity
  const double diam = actmat->Diam();
  // Read in blood viscosity
  const double visc = actmat->Viscosity();

  const double hag_pois = M_PI * pow(diam, 4) / 128.0 / visc;
  // gaussian points
  const DRT::UTILS::IntegrationPoints1D intpoints(ele->GaussRule());

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
    DRT::UTILS::shape_function_1D_deriv1(my::deriv_, xi, distype);

    for (int inode = 0; inode < numnode; inode++)
      for (int jnode = 0; jnode < numnode; jnode++)
        sysmat(inode, jnode) += my::deriv_(0, inode) * fac * my::deriv_(0, jnode);

    LINALG::Matrix<1, 1> pressgrad;
    pressgrad.Multiply(my::deriv_, mypress);
    for (int inode = 0; inode < numnode; inode++)
      rhs(inode) -= my::deriv_(0, inode) * fac * pressgrad(0, 0);
  }


  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes

// 1D elements
template class DRT::ELEMENTS::ArteryEleCalcPresBased<DRT::Element::line2>;
