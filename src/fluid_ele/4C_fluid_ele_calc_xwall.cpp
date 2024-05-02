/*----------------------------------------------------------------------*/
/*! \file

\brief main file containing routines for calculation of fluid element with xfem wall modeling

\level 2


*/
/*----------------------------------------------------------------------*/

#include "4C_fluid_ele_calc_xwall.hpp"

#include "4C_discretization_fem_general_extract_values.hpp"
#include "4C_discretization_fem_general_utils_gder2.hpp"
#include "4C_fluid_ele.hpp"
#include "4C_fluid_ele_action.hpp"
#include "4C_fluid_ele_parameter_std.hpp"
#include "4C_fluid_ele_parameter_timint.hpp"
#include "4C_fluid_rotsym_periodicbc.hpp"
#include "4C_lib_condition_utils.hpp"
#include "4C_linalg_utils_densematrix_eigen.hpp"
#include "4C_mat_newtonianfluid.hpp"

FOUR_C_NAMESPACE_OPEN


/*-----------------------------------------------------------------------------*
 | Constructor                                                      bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::FluidEleCalcXWall()
    : DRT::ELEMENTS::FluidEleCalc<distype, enrtype>::FluidEleCalc(),
      ewdist_(true),
      etauw_(true),
      einctauw_(true),
      eramp_(true),
      epsi_(true),
      epsinew_(true),
      epsiold_(true),
      eincwdist_(true),
      visc_(0.0),
      viscinv_(0.0),
      xyze_(true),
      functenr_(true),
      funct_(true),
      derxyenr_(true),
      derxy_(true),
      derxyenr2_(true),
      derxy2_(true),
      deriv_(true),
      deriv2_(true),
      k_(0.41),
      b_(5.17),
      expmkmb_(exp(-k_ * b_)),
      mk_(-1.0)
{
}

template <CORE::FE::CellType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>*
DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::Instance(CORE::UTILS::SingletonAction action)
{
  static CORE::UTILS::SingletonOwner<DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>>
      singleton_owner(
          []()
          {
            return std::unique_ptr<DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>>(
                new DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>());
          });

  return singleton_owner.Instance(action);
}

/*-----------------------------------------------------------------------------*
 | Entry supporting methods of the element                          bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::EvaluateService(DRT::ELEMENTS::Fluid* ele,
    Teuchos::ParameterList& params, Teuchos::RCP<CORE::MAT::Material>& mat,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  calcoldandnewpsi_ = false;
  const FLD::Action act = CORE::UTILS::GetAsEnum<FLD::Action>(params, "action");
  if (act == FLD::xwall_l2_projection) calcoldandnewpsi_ = true;
  GetEleProperties(ele, discretization, lm, params, mat);

  // non-enriched case, solve problem as usual
  if (enrtype == DRT::ELEMENTS::Fluid::xwall)  // this element has the same no of dofs on each node
  {
    std::vector<int> assembletoggle;
    int nodecount = 0;
    for (std::vector<int>::const_iterator i = lm.begin(); i != lm.end(); ++i)
    {
      ++nodecount;
      if (nodecount == 8)  // change here if I use a different number of dofs
      {
        assembletoggle.push_back(0);
        nodecount = 0;
      }
      else
      {
        assembletoggle.push_back(1);
      }
    }

    // the last dof must have been an unused pressure dof
    if (nodecount != 0)
      FOUR_C_THROW("something is wrong in this element with the number of virtual nodes vs dofs");

    int err = EvaluateServiceXWall(
        ele, params, mat, discretization, lm, elemat1, elemat2, elevec1, elevec2, elevec3);

    // for some EvaluateService actions, elevec1 is not necessary
    if (elevec1.length() != 0 && act != FLD::tauw_via_gradient && act != FLD::calc_div_u &&
        act != FLD::calc_dt_via_cfl && act != FLD::xwall_calc_mk &&
        act != FLD::calc_mass_flow_periodic_hill)
    {
      int row1 = 0;
      // assembly back into the old vector
      for (std::vector<int>::const_iterator i = assembletoggle.begin(); i != assembletoggle.end();
           ++i)
      {
        if (*i == 0) elevec1[row1] = 0.0;
        ++row1;
      }
    }

    return err;
  }
  else
    FOUR_C_THROW("not xwall element");


  return 1;
}


/*----------------------------------------------------------------------*
 * Evaluate supporting methods of the element
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::EvaluateServiceXWall(
    DRT::ELEMENTS::Fluid* ele, Teuchos::ParameterList& params,
    Teuchos::RCP<CORE::MAT::Material>& mat, DRT::Discretization& discretization,
    std::vector<int>& lm, CORE::LINALG::SerialDenseMatrix& elemat1,
    CORE::LINALG::SerialDenseMatrix& elemat2, CORE::LINALG::SerialDenseVector& elevec1,
    CORE::LINALG::SerialDenseVector& elevec2, CORE::LINALG::SerialDenseVector& elevec3)
{
  // get the action required
  const FLD::Action act = CORE::UTILS::GetAsEnum<FLD::Action>(params, "action");

  switch (act)
  {
    case FLD::xwall_l2_projection:
    {
      // this action only considers enriched elements
      // I only need funct_ and maybe the first derivative
      my::is_higher_order_ele_ = false;
      return XWallProjection(ele, params, discretization, lm, mat, elemat1,
          elemat2);  // here we can go further in the element and implement the matrix and rhs terms
    }
    break;
    case FLD::tauw_via_gradient:
    {
      // this action only considers enriched elements
      // I only need funct_ and maybe the first derivative
      my::is_higher_order_ele_ = false;
      return TauWViaGradient(ele, params, discretization, lm, mat, elevec1,
          elevec2);  // here we can go further in the element and implement the matrix and rhs terms
    }
    break;
    case FLD::xwall_calc_mk:
    {
      // this action only considers enriched elements
      // It is essential that the second derivatives exist!
      my::is_higher_order_ele_ = true;
      return CalcMK(ele, params, discretization, lm, mat, elevec1,
          elevec2);  // here we can go further in the element and implement the matrix and rhs terms
    }
    break;
    default:
      return my::EvaluateService(
          ele, params, mat, discretization, lm, elemat1, elemat2, elevec1, elevec2, elevec3);
      break;
  }  // end of switch(act)

  return 0;
}

/*-----------------------------------------------------------------------------*
 | Action type: Evaluate                                            bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::Evaluate(DRT::ELEMENTS::Fluid* ele,
    DRT::Discretization& discretization, const std::vector<int>& lm, Teuchos::ParameterList& params,
    Teuchos::RCP<CORE::MAT::Material>& mat, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra, bool offdiag)
{
  calcoldandnewpsi_ = false;

  GetEleProperties(ele, discretization, lm, params, mat);

  if (enrtype == DRT::ELEMENTS::Fluid::xwall)  // this element has the same no of dofs on each node
  {
    std::vector<int> assembletoggle;
    int nodecount = 0;
    for (std::vector<int>::const_iterator i = lm.begin(); i != lm.end(); ++i)
    {
      ++nodecount;
      if (nodecount == 8)
      {
        nodecount = 0;
        assembletoggle.push_back(0);
      }
      else
      {
        assembletoggle.push_back(1);
      }
    }

    // the last dof must have been an unused pressure dof
    if (nodecount != 0)
      FOUR_C_THROW("something is wrong in this element with the number of virtual nodes vs dofs");

    int err = my::Evaluate(ele, discretization, lm, params, mat, elemat1_epetra, elemat2_epetra,
        elevec1_epetra, elevec2_epetra, elevec3_epetra, my::intpoints_);

    int row1 = 0;
    int col1 = 0;
    // assembly back into the old matrix
    for (std::vector<int>::const_iterator i = assembletoggle.begin(); i != assembletoggle.end();
         ++i)
    {
      row1 = 0;
      for (std::vector<int>::const_iterator j = assembletoggle.begin(); j != assembletoggle.end();
           ++j)
      {
        if (*i == 0 && *j == 0 && row1 == col1)
          elemat1_epetra[col1][row1] = 1.0;
        else if (*i == 0 or *j == 0)
          elemat1_epetra[col1][row1] = 0.0;
        ++row1;
      }
      ++col1;
    }

    row1 = 0;
    // assembly back into the old matrix
    for (std::vector<int>::const_iterator i = assembletoggle.begin(); i != assembletoggle.end();
         ++i)
    {
      if (*i == 0) elevec1_epetra[row1] = 0.0;
      ++row1;
    }

    return err;
  }
  else
    FOUR_C_THROW(
        "this should not have happended: some nodes have too many dofs in the LM vector, because "
        "they are dof-blending nodes and the wrong LocationVector() function is called");


  return 1;
}

/*-----------------------------------------------------------------------------*
 | Get properties for this element                                  bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::GetEleProperties(DRT::ELEMENTS::Fluid* ele,
    DRT::Discretization& discretization, const std::vector<int>& lm, Teuchos::ParameterList& params,
    Teuchos::RCP<CORE::MAT::Material>& mat)
{
  is_blending_ele_ = false;
  visc_ = 0.0;

  if (!(enrtype == DRT::ELEMENTS::Fluid::xwall))
    FOUR_C_THROW("This class is exclusively for the xwall enrichment type up to now");

  // rotate the vector field in the case of rotationally symmetric boundary conditions
  if (my::rotsymmpbc_->HasRotSymmPBC()) FOUR_C_THROW("rotsymm pbc don't work with xwall");

  // get xwall toggle
  {
    const Teuchos::RCP<Epetra_Vector> xwalltoggle =
        params.get<Teuchos::RCP<Epetra_Vector>>("xwalltoggle");

    std::vector<double> mylocal(ele->NumNode());
    CORE::FE::ExtractMyNodeBasedValues(ele, mylocal, *xwalltoggle);

    for (unsigned inode = 0; inode < (unsigned)enren_; ++inode)  // number of nodes
    {
      etoggle_(inode) = mylocal.at(inode);
    }
  }

  int xwallnodes = 0;
  // init nodal ramp values
  for (int pq = 0; pq < enren_; ++pq)
  {
    if (etoggle_(pq) > 0.5)
    {
      ++xwallnodes;
      eramp_(pq) = 1.0;
    }
    else
      eramp_(pq) = 0.0;
  }

  if (xwallnodes > 0 && xwallnodes < enren_) is_blending_ele_ = true;

  // get wall distance
  {
    const Teuchos::RCP<Epetra_Vector> walldist =
        params.get<Teuchos::RCP<Epetra_Vector>>("walldist");
    //      std::cout << *walldist << std::endl;
    std::vector<double> mylocal(ele->NumNode());
    CORE::FE::ExtractMyNodeBasedValues(ele, mylocal, *walldist);

    for (unsigned inode = 0; inode < (unsigned)enren_; ++inode)  // number of nodes
    {
      ewdist_(inode) = mylocal.at(inode);
    }
  }

  // get tauw
  {
    const Teuchos::RCP<Epetra_Vector> tauw = params.get<Teuchos::RCP<Epetra_Vector>>("tauw");

    std::vector<double> mylocal(ele->NumNode());
    CORE::FE::ExtractMyNodeBasedValues(ele, mylocal, *tauw);

    for (unsigned inode = 0; inode < (unsigned)enren_; ++inode)  // number of nodes
    {
      etauw_(inode) = mylocal.at(inode);
    }
  }

  // get increment of tauw
  {
    const Teuchos::RCP<Epetra_Vector> inctauw = params.get<Teuchos::RCP<Epetra_Vector>>("inctauw");

    std::vector<double> mylocal(ele->NumNode());
    CORE::FE::ExtractMyNodeBasedValues(ele, mylocal, *inctauw);

    for (unsigned inode = 0; inode < (unsigned)enren_; ++inode)  // number of nodes
    {
      einctauw_(inode) = mylocal.at(inode);
    }
  }

  // get viscosity and density
  {
    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
    if (!actmat) FOUR_C_THROW("not a newtonian fluid");
    // get constant dynamic viscosity
    dens_ = actmat->Density();
    densinv_ = 1.0 / dens_;
    visc_ = actmat->Viscosity() * densinv_;
    viscinv_ = 1.0 / visc_;
  }

  // calculate nodal shape values of psi
  for (int inode = 0; inode < enren_; inode++)
  {
    double utaunode = sqrt(etauw_(inode) * densinv_);
    double psinode = SpaldingsLaw(ewdist_(inode), utaunode);

    epsi_(inode) = psinode;
  }

  if (calcoldandnewpsi_ == true)
  {
    // get old wall distance in case of ale
    if (ele->IsAle())
    {
      const Teuchos::RCP<Epetra_Vector> incwdist =
          params.get<Teuchos::RCP<Epetra_Vector>>("incwalldist");

      std::vector<double> mylocal(ele->NumNode());
      CORE::FE::ExtractMyNodeBasedValues(ele, mylocal, *incwdist);

      for (unsigned inode = 0; inode < (unsigned)enren_; ++inode)  // number of nodes
      {
        eincwdist_(inode) = mylocal.at(inode);
      }
    }

    for (int inode = 0; inode < enren_; inode++)
    {
      epsinew_(inode) = epsi_(inode);
      double utaunode = sqrt((etauw_(inode) - einctauw_(inode)) * densinv_);
      double psinode = SpaldingsLaw(ewdist_(inode) - eincwdist_(inode), utaunode);

      epsiold_(inode) = psinode;
    }
  }


  // get element mk for stabilization
  const Teuchos::RCP<Epetra_Vector> mkvec = params.get<Teuchos::RCP<Epetra_Vector>>("mk");
  mk_ = (*mkvec)[mkvec->Map().LID(ele->Id())];

  numgpnorm_ = params.get<int>("gpnorm");
  numgpnormow_ = params.get<int>("gpnormow");
  numgpplane_ = params.get<int>("gppar");
  // get node coordinates and number of elements per node
  CORE::GEO::fillInitialPositionArray<distype, nsd_, CORE::LINALG::Matrix<nsd_, nen_>>(
      ele, my::xyze_);
  CORE::LINALG::Matrix<nsd_, nen_> edispnp(true);
  if (ele->IsAle()) GetGridDispALE(discretization, lm, edispnp);
  PrepareGaussRule();

  return;
}

/*-----------------------------------------------------------------------------*
 | Go wall shear stress increment backwards                         bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::XWallTauWIncBack()
{
  for (unsigned inode = 0; inode < (unsigned)enren_; ++inode)  // number of nodes
  {
    etauw_(inode) -= einctauw_(inode);
    epsi_(inode) = epsiold_(inode);
    ewdist_(inode) -= eincwdist_(inode);
  }

  return;
}

/*-----------------------------------------------------------------------------*
 | Go wall shear stress increment forward                           bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::XWallTauWIncForward()
{
  for (unsigned inode = 0; inode < (unsigned)enren_; ++inode)  // number of nodes
  {
    etauw_(inode) += einctauw_(inode);
    epsi_(inode) = epsinew_(inode);
    ewdist_(inode) += eincwdist_(inode);
  }

  return;
}

/*-----------------------------------------------------------------------------*
 | Calculate shape functions at integration point                   bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::EvalShapeFuncAndDerivsAtIntPoint(
    const double* gpcoord,  // actual integration point (coords)
    double gpweight         // actual integration point (weight)
)
{
  EvalStdShapeFuncAndDerivsAtIntPoint(gpcoord, gpweight);

  EvalEnrichment();
  return;
}

/*-----------------------------------------------------------------------------*
 | Calculate shape functions at integration point                   bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::EvalStdShapeFuncAndDerivsAtIntPoint(
    const double* gpcoord,  // actual integration point (coords)
    double gpweight         // actual integration point (weight)
)
{
  funct_.Clear();
  derxy_.Clear();
  derxy2_.Clear();

  // put the geometry in the local copy
  for (int inode = 0; inode < enren_; ++inode)
  {
    for (int sdm = 0; sdm < nsd_; ++sdm) xyze_(sdm, inode) = my::xyze_(sdm, inode);
  }

  // coordinates of the current integration point
  for (int idim = 0; idim < nsd_; idim++)
  {
    my::xsi_(idim) = gpcoord[idim];
  }

  // shape functions and their first derivatives
  CORE::FE::shape_function<distype>(my::xsi_, funct_);
  CORE::FE::shape_function_deriv1<distype>(my::xsi_, deriv_);
  if (my::is_higher_order_ele_ && distype == CORE::FE::CellType::hex8)
  {
    // get the second derivatives of standard element at current GP
    CORE::FE::shape_function_deriv2<distype>(my::xsi_, deriv2_);
  }
  else
    deriv2_.Clear();

  // get Jacobian matrix and determinant
  // actually compute its transpose....
  /*
    +-            -+ T      +-            -+
    | dx   dx   dx |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dr   dr   dr |
    |              |        |              |
    | dy   dy   dy |        | dx   dy   dz |
    | --   --   -- |   =    | --   --   -- |
    | dr   ds   dt |        | ds   ds   ds |
    |              |        |              |
    | dz   dz   dz |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dt   dt   dt |
    +-            -+        +-            -+
  */
  my::xjm_.MultiplyNT(deriv_, xyze_);

  my::det_ = my::xji_.Invert(my::xjm_);

  if (my::det_ < 1E-16)
    FOUR_C_THROW(
        "GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", my::eid_, my::det_);

  // compute integration factor
  my::fac_ = gpweight * my::det_;

  // compute global first derivates
  derxy_.Multiply(my::xji_, deriv_);

  //--------------------------------------------------------------
  //             compute global second derivatives
  //--------------------------------------------------------------
  if (my::is_higher_order_ele_ && distype == CORE::FE::CellType::hex8)
  {
    CORE::FE::gder2<distype, enren_>(my::xjm_, derxy_, deriv2_, xyze_, derxy2_);
  }
  else
    derxy2_.Clear();

  return;
}

/*-----------------------------------------------------------------------------*
 | Calculate enrichment shape functions                             bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::EvalEnrichment()
{
  // first clear everything
  functenr_.Clear();
  derxyenr_.Clear();
  derxyenr2_.Clear();

  CORE::LINALG::Matrix<nsd_, 1> derpsigp(true);
  CORE::LINALG::Matrix<numderiv2_, 1> der2psigp(true);

  double psigp = EnrichmentShapeDer(derpsigp, der2psigp);

  // shape function
  for (int inode = 0; inode < enren_; ++inode)
    functenr_(inode) = funct_(inode) * (psigp - epsi_(inode));

  // first derivative
  for (int inode = 0; inode < enren_; ++inode)
  {
    for (int sdm = 0; sdm < nsd_; ++sdm)
    {
      derxyenr_(sdm, inode) =
          derxy_(sdm, inode) * (psigp - epsi_(inode)) + funct_(inode) * derpsigp(sdm);
    }
  }

  if (my::is_higher_order_ele_)
  {
    for (int inode = 0; inode < enren_; ++inode)
    {
      for (int sdm = 0; sdm < numderiv2_; ++sdm)
      {
        const std::array<int, 6> i = {0, 1, 2, 0, 0, 1};
        const std::array<int, 6> j = {0, 1, 2, 1, 2, 2};
        derxyenr2_(sdm, inode) += derxy2_(sdm, inode) * (psigp - epsi_(inode)) +
                                  derxy_(i[sdm], inode) * derpsigp(j[sdm]) +
                                  derxy_(j[sdm], inode) * derpsigp(i[sdm]) +
                                  funct_(inode) * der2psigp(sdm);
      }
    }
  }

  // treat blending elements with ramp functions
  if (is_blending_ele_)
  {
    CORE::LINALG::Matrix<nsd_, 1> derramp(true);
    derramp.Multiply(derxy_, eramp_);
    CORE::LINALG::Matrix<numderiv2_, 1> der2ramp(true);
    der2ramp.Multiply(derxy2_, eramp_);
    double ramp = eramp_.Dot(funct_);

    for (int inode = 0; inode < enren_; ++inode)
    {
      if (my::is_higher_order_ele_)
      {
        for (int sdm = 0; sdm < numderiv2_; ++sdm)
        {
          const std::array<int, 6> i = {0, 1, 2, 0, 0, 1};
          const std::array<int, 6> j = {0, 1, 2, 1, 2, 2};

          derxyenr2_(sdm, inode) *= ramp;
          derxyenr2_(sdm, inode) += derxyenr_(i[sdm], inode) * derramp(j[sdm]) +
                                    derxyenr_(j[sdm], inode) * derramp(i[sdm]);
          derxyenr2_(sdm, inode) += functenr_(inode) * der2ramp(sdm);
        }
      }
      for (int sdm = 0; sdm < nsd_; ++sdm)
      {
        derxyenr_(sdm, inode) = derxyenr_(sdm, inode) * ramp + functenr_(inode) * derramp(sdm);
      }
      // treat derivative first, because we use functenr_ without ramp_
      functenr_(inode) *= ramp;
    }
  }

  // put everything in standard shape function
  for (int inode = 0; inode < enren_; ++inode)
  {
    const int inodetwo = inode * 2;
    my::funct_(inodetwo) = funct_(inode);
    my::funct_(inodetwo + 1) = functenr_(inode);
    for (int sdm = 0; sdm < nsd_; ++sdm)
    {
      my::derxy_(sdm, inodetwo) = derxy_(sdm, inode);
      my::derxy_(sdm, inodetwo + 1) = derxyenr_(sdm, inode);
    }

    if (my::is_higher_order_ele_)
    {
      for (int sdm = 0; sdm < numderiv2_; ++sdm)
      {
        my::derxy2_(sdm, inodetwo) = derxy2_(sdm, inode);
        my::derxy2_(sdm, inodetwo + 1) = derxyenr2_(sdm, inode);
      }
    }
  }

  return;
}

/*-----------------------------------------------------------------------------*
 | Calculate enrichment shape functions and derivatives                        |
 | including transformation from y+ to (x,y,z)                      bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
double DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::EnrichmentShapeDer(
    CORE::LINALG::Matrix<nsd_, 1>& derpsigp, CORE::LINALG::Matrix<numderiv2_, 1>& der2psigp)
{
  // calculate transformation ---------------------------------------
  double wdist = ewdist_.Dot(funct_);
  double tauw = etauw_.Dot(funct_);
  CORE::LINALG::Matrix<nsd_, 1> derwdist(true);
  derwdist.Multiply(derxy_, ewdist_);
  CORE::LINALG::Matrix<nsd_, 1> dertauw(true);
  dertauw.Multiply(derxy_, etauw_);
  CORE::LINALG::Matrix<numderiv2_, 1> der2wdist(true);
  if (my::is_higher_order_ele_) der2wdist.Multiply(derxy2_, ewdist_);
  CORE::LINALG::Matrix<numderiv2_, 1> der2tauw(true);
  if (my::is_higher_order_ele_) der2tauw.Multiply(derxy2_, etauw_);
  CORE::LINALG::Matrix<nsd_, 1> dertrans(true);
  CORE::LINALG::Matrix<numderiv2_, 1> der2trans_1(true);
  CORE::LINALG::Matrix<numderiv2_, 1> der2trans_2(true);

  if (tauw < 1.0e-10) FOUR_C_THROW("tauw is almost zero");
  if (dens_ < 1.0e-10) FOUR_C_THROW("density is almost zero");

  const double utau = sqrt(tauw * densinv_);
  const double fac = 1.0 / (2.0 * sqrt(dens_ * tauw));
  const double wdistfac = wdist * fac;

  for (int sdm = 0; sdm < nsd_; ++sdm)
    dertrans(sdm) = (utau * derwdist(sdm) + wdistfac * dertauw(sdm)) * viscinv_;


  // second derivative, first part: to be multiplied with der2psigpsc
  // second derivative, second part: to be multiplied with derpsigpsc
  if (my::is_higher_order_ele_)
  {
    const double wdistfactauwtwoinv = wdistfac / (tauw * 2.0);

    for (int sdm = 0; sdm < numderiv2_; ++sdm)
    {
      const std::array<int, 6> i = {0, 1, 2, 0, 0, 1};
      const std::array<int, 6> j = {0, 1, 2, 1, 2, 2};

      der2trans_1(sdm) = dertrans(i[sdm]) * dertrans(j[sdm]);

      der2trans_2(sdm) = (derwdist(j[sdm]) * fac * dertauw(i[sdm]) + wdistfac * der2tauw(sdm) -
                             wdistfactauwtwoinv * dertauw(i[sdm]) * dertauw(j[sdm]) +
                             dertauw(j[sdm]) * fac * derwdist(i[sdm]) + utau * der2wdist(sdm)) *
                         viscinv_;
    }
  }
  // calculate transformation done ----------------------------------

  // get enrichment function and scalar derivatives
  const double psigp = SpaldingsLaw(wdist, utau);
  const double derpsigpsc = DerSpaldingsLaw(wdist, utau, psigp);
  const double der2psigpsc = Der2SpaldingsLaw(wdist, utau, psigp, derpsigpsc);

  // calculate final derivatives
  for (int sdm = 0; sdm < nsd_; ++sdm)
  {
    derpsigp(sdm) = derpsigpsc * dertrans(sdm);
  }
  if (my::is_higher_order_ele_)
    for (int sdm = 0; sdm < numderiv2_; ++sdm)
    {
      der2psigp(sdm) = der2psigpsc * der2trans_1(sdm);
      der2psigp(sdm) += derpsigpsc * der2trans_2(sdm);
    }

  return psigp;
}

/*-----------------------------------------------------------------------------*
 | Enrichment function (modification of Spalding's law)             bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
double DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::SpaldingsLaw(double dist, double utau)
{
  // watch out, this is not exactly Spalding's law but psi=u_+*k, which saves quite some
  // multiplications

  double yplus = dist * utau * viscinv_;
  double psi = 0.0;

  const double km1 = 1.0 / k_;

  if (yplus >
      11.0)  // this is approximately where the intersection of log law and linear region lies
    psi = log(yplus) + b_ * k_;
  else
    psi = yplus * k_;

  double inc = 10.0;
  double fn = 10.0;
  int count = 0;
  while (abs(inc) > 1.0E-14 && abs(fn) > 1.0E-14 && 1000 > count++)
  {
    double psiquad = psi * psi;
    double exppsi = exp(psi);
    fn = -yplus + psi * km1 +
         expmkmb_ *
             (exppsi - 1.0 - psi - psiquad * 0.5 - psiquad * psi / 6.0 - psiquad * psiquad / 24.0);
    double dfn = km1 + expmkmb_ * (exppsi - 1.0 - psi - psiquad * 0.5 - psiquad * psi / 6.0);

    inc = fn / dfn;

    psi -= inc;
  }

  return psi;

  // Reichardt's law 1951
  // return (1.0/k_*log(1.0+0.4*yplus)+7.8*(1.0-exp(-yplus/11.0)-(yplus/11.0)*exp(-yplus/3.0)))*k_;
}

/*-----------------------------------------------------------------------------*
 | Derivative of enrichment function w.r.t. y+                         bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
double DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::DerSpaldingsLaw(
    double dist, double utau, double psi)
{
  // derivative with respect to y+!
  // spaldings law according to paper (derivative)
  return 1.0 /
         (1.0 / k_ + expmkmb_ * (exp(psi) - 1.0 - psi - psi * psi * 0.5 - psi * psi * psi / 6.0));

  // Reichardt's law
  //  double yplus=dist*utau*viscinv_;
  //  return
  //  (0.4/(k_*(1.0+0.4*yplus))+7.8*(1.0/11.0*exp(-yplus/11.0)-1.0/11.0*exp(-yplus/3.0)+yplus/33.0*exp(-yplus/3.0)))*k_;
}

/*-----------------------------------------------------------------------------*
 | Second derivative of enrichment function w.r.t. y+               bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
double DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::Der2SpaldingsLaw(
    double dist, double utau, double psi, double derpsi)
{
  // derivative with respect to y+!
  // spaldings law according to paper (2nd derivative)
  return -expmkmb_ * (exp(psi) - 1 - psi - psi * psi * 0.5) * derpsi * derpsi * derpsi;

  // Reichardt's law
  //  double yplus=dist*utau*viscinv_;
  //  return
  //  (-0.4*0.4/(k_*(1.0+0.4*yplus)*(1.0+0.4*yplus))+7.8*(-1.0/121.0*exp(-yplus/11.0)+(2.0/33.0-yplus/99.0)*exp(-yplus/3.0)))*k_;
}

/*-----------------------------------------------------------------------------*
 | Calculate matrix for l2 projection                               bk 07/2014 |
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::TauWViaGradient(DRT::ELEMENTS::Fluid* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::RCP<CORE::MAT::Material>& mat, CORE::LINALG::SerialDenseVector& elevec1,
    CORE::LINALG::SerialDenseVector& elevec2)
{
  //----------------------------------------------------------------------------
  //   Extract velocity/pressure from global vectors
  //----------------------------------------------------------------------------

  CORE::LINALG::Matrix<nsd_, nen_> evel(true);
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &evel, nullptr, "vel");

  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  CORE::GEO::fillInitialPositionArray<distype, nsd_, CORE::LINALG::Matrix<nsd_, nen_>>(
      ele, my::xyze_);

  if (ele->IsAle())
  {
    CORE::LINALG::Matrix<nsd_, nen_> edispnp(true);
    CORE::LINALG::Matrix<nsd_, nen_> egridv(true);
    my::GetGridDispVelALE(discretization, lm, edispnp, egridv);
    evel -= egridv;
  }


  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------


  // number of nodes, the derivative should be calculated at the nodes
  for (int inode = 0; inode < enren_; ++inode)
  {
    // calculate only for the wall nodes
    if (ewdist_(inode) < 1e-4)
    {
      CORE::LINALG::Matrix<3, 1> test = CORE::FE::GetNodeCoordinates(inode, distype);
      const std::array<double, 3> gp = {test(0, 0), test(1, 0), test(2, 0)};
      const double* gpc = gp.data();
      // evaluate shape functions and derivatives at integration point
      EvalShapeFuncAndDerivsAtIntPoint(gpc, 1.0);

      // calculate wall-normal vector
      CORE::LINALG::Matrix<nsd_, 1> normwall(true);
      normwall.Multiply(derxy_, ewdist_);

      // at certain corner elements, it can happen, that the normal vector calculated at the
      // boundary is zero. so we move a bit inside the element and calculate the properties there
      // instead
      if (normwall.Norm2() < 1.0e-10)
      {
        test.Scale(0.95);
        const std::array<double, 3> gp = {test(0, 0), test(1, 0), test(2, 0)};
        const double* gpc = gp.data();
        // evaluate shape functions and derivatives at integration point
        EvalShapeFuncAndDerivsAtIntPoint(gpc, 1.0);

        normwall.Multiply(derxy_, ewdist_);
        if (normwall.Norm2() < 1.0e-10)
          FOUR_C_THROW("normal vector has length zero, even in the second try");
      }

      // unit vector
      normwall.Scale(1.0 / normwall.Norm2());

      CORE::LINALG::Matrix<nsd_, nsd_> velderxy(true);
      velderxy.MultiplyNT(evel, my::derxy_);

      // remove normal part

      //      normwall.Scale(1.0/normwall.Norm2());
      CORE::LINALG::Matrix<nsd_, nsd_> velderxywoun(true);
      for (int idim = 0; idim < nsd_; idim++)
        for (int jdim = 0; jdim < nsd_; jdim++)
          velderxywoun(idim, jdim) = velderxy(idim, jdim) * (1.0 - abs(normwall(idim)));

      // now transform to derivative w.r.t. n
      CORE::LINALG::Matrix<nsd_, 1> veldern(true);
      for (int idim = 0; idim < nsd_; idim++)
        for (int jdim = 0; jdim < nsd_; jdim++)
          veldern(idim) += velderxywoun(idim, jdim) * normwall(jdim);

      // calculate wall shear stress with dynamic! viscosity
      elevec1(inode) = visc_ * dens_ * veldern.Norm2();
      // the following vector counts, how often this node is assembled
      //(e.g. from neighboring elements)
      elevec2(inode) = 1.0;
    }
  }  // end of integration loop

  return 0;
}

/*-----------------------------------------------------------------------------*
 | Calculate stabilization parameter mk                             bk 07/2014 |
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
double DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::GetMK()
{
  if (mk_ < 0.0)
    return CalcMK();
  else
    return mk_;
  FOUR_C_THROW("mk could not be determined for xwall");
  return 0.0;
}


/*-----------------------------------------------------------------------------*
 | Calculate stabilization parameter mk                             bk 07/2014 |
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
double DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::CalcMK()
{
  if (my::is_higher_order_ele_ == false)
    FOUR_C_THROW("It is essential that the second derivatives exist!");

  CORE::FE::GaussIntegration intpoints(CORE::FE::CellType::hex8, 1);
  if (distype == CORE::FE::CellType::hex8)
  {
    CORE::FE::GaussIntegration intpointstmp(cgp_);
    intpoints = intpointstmp;
  }
  else if (distype == CORE::FE::CellType::tet4)
  {
    CORE::FE::GaussIntegration intpointsplane(CORE::FE::CellType::tet4, 2 * numgpnorm_ - 1);
    intpoints = intpointsplane;
  }

  CORE::LINALG::SerialDenseMatrix elemat_epetra1;
  CORE::LINALG::SerialDenseMatrix elemat_epetra2;
  elemat_epetra1.shape(nen_, nen_);
  elemat_epetra2.shape(nen_, nen_);
  CORE::LINALG::Matrix<nen_, nen_> Amat(elemat_epetra1.values(), true);
  CORE::LINALG::Matrix<nen_, nen_> Bmat(elemat_epetra2.values(), true);

  double vol = 0.0;
  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------
  for (CORE::FE::GaussIntegration::iterator iquad = intpoints.begin(); iquad != intpoints.end();
       ++iquad)
  {
    // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(), iquad.Weight());

    const unsigned Velx = 0;
    const unsigned Vely = 1;
    const unsigned Velz = 2;

    for (int vi = 0; vi < nen_; ++vi)
    {
      for (int ui = 0; ui < nen_; ++ui)
      {
        //(x,x)
        // mixed derivatives are neglected such that the standard values are recovered for ideal
        // higher order cube elements, i.e. C_l={12,60} for polynomial orders p={2,3}. See Harari &
        // Hughes 1991 mixed derivatives practically don't have any influence
        Amat(vi, ui) += my::fac_ * (my::derxy2_(Velx, vi) * my::derxy2_(Velx, ui)) +
                        my::fac_ * (my::derxy2_(Vely, vi) * my::derxy2_(Vely, ui)) +
                        my::fac_ * (my::derxy2_(Velz, vi) * my::derxy2_(Velz, ui));
      }
    }

    /*
  //    /                \
  //   |                  |
  //   |Nabla(v),Nabla(u) |
  //   |                  |
  //    \                / volume
     */

    for (int vi = 0; vi < nen_; ++vi)
    {
      for (int ui = 0; ui < nen_; ++ui)
      {
        //(x,x)
        Bmat(vi, ui) += my::fac_ * (my::derxy_(Velx, vi) * my::derxy_(Velx, ui)) +
                        my::fac_ * (my::derxy_(Vely, vi) * my::derxy_(Vely, ui)) +
                        my::fac_ * (my::derxy_(Velz, vi) * my::derxy_(Velz, ui));
      }
    }

    vol += my::fac_;
  }  // gauss loop

  const double maxeigenvalue = CORE::LINALG::GeneralizedEigen(elemat_epetra1, elemat_epetra2);

  double h_u = 0.0;
  if (my::fldpara_->WhichTau() == INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall ||
      my::fldpara_->WhichTau() == INPAR::FLUID::tau_codina ||
      my::fldpara_->WhichTau() == INPAR::FLUID::tau_codina_convscaled)
  {
    if (!(my::fldpara_->CharEleLengthU() == INPAR::FLUID::volume_equivalent_diameter_u))
      FOUR_C_THROW("only volume equivalent diameter defined up to now");

    // volume equivalent diameter
    h_u = std::pow((6. * vol / M_PI), (1.0 / 3.0)) / sqrt(3.0);
  }
  else
    FOUR_C_THROW(
        "Element length not defined for dynamic determination of mk for your stabilization "
        "parameter");

  if (abs(maxeigenvalue) < 1.e-9)
  {
    std::cout << "Warning: maxeigenvalue zero:  " << maxeigenvalue << std::endl;
    return 0.33333333333;
  }
  else if (1.0 / (maxeigenvalue * h_u * h_u) > 0.33)
  {
    std::cout << "Warning: mk larger than 0.33:  " << maxeigenvalue * h_u * h_u << std::endl;

    return 0.33333333333;
  }

  // safety factor
  const double sfac = 1.0;
  //  std::cout << sfac/(maxeigenvalue*h_u*h_u) << std::endl;
  return sfac / (maxeigenvalue * h_u * h_u);
}

/*-----------------------------------------------------------------------------*
 | Calculate stabilization parameter mk                             bk 07/2014 |
 | (call for action type)                                                      |
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::CalcMK(DRT::ELEMENTS::Fluid* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::RCP<CORE::MAT::Material>& mat, CORE::LINALG::SerialDenseVector& elevec1,
    CORE::LINALG::SerialDenseVector& elevec2)
{
  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  CORE::GEO::fillInitialPositionArray<distype, nsd_, CORE::LINALG::Matrix<nsd_, nen_>>(
      ele, my::xyze_);

  CORE::LINALG::Matrix<nsd_, nen_> edispnp(true);
  if (ele->IsAle()) GetGridDispALE(discretization, lm, edispnp);

  elevec1[0] = CalcMK();
  return 0;
}

/*-----------------------------------------------------------------------------*
 | Calculate matrix for l2 projection                               bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::XWallProjection(DRT::ELEMENTS::Fluid* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::RCP<CORE::MAT::Material>& mat, CORE::LINALG::SerialDenseMatrix& elemat1,
    CORE::LINALG::SerialDenseMatrix& elemat2)
{
  const int numdof = 3;

  //----------------------------------------------------------------------------
  //   Extract velocity/pressure from global vectors
  //----------------------------------------------------------------------------

  CORE::LINALG::Matrix<nsd_, nen_> eveln(true);
  my::ExtractValuesFromGlobalVector(discretization, lm, *my::rotsymmpbc_, &eveln, nullptr, "veln");

  CORE::LINALG::Matrix<nsd_, nen_> eaccn(true);
  bool switchonaccn = discretization.HasState("accn");
  if (switchonaccn)
    my::ExtractValuesFromGlobalVector(
        discretization, lm, *my::rotsymmpbc_, &eaccn, nullptr, "accn");

  CORE::LINALG::Matrix<nsd_, nen_> evelnp(true);
  bool switchonvelnp = discretization.HasState("velnp");
  if (switchonvelnp)
    my::ExtractValuesFromGlobalVector(
        discretization, lm, *my::rotsymmpbc_, &evelnp, nullptr, "velnp");

  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  CORE::GEO::fillInitialPositionArray<distype, nsd_, CORE::LINALG::Matrix<nsd_, nen_>>(
      ele, my::xyze_);

  CORE::LINALG::Matrix<nsd_, nen_> edispnp(true);
  if (ele->IsAle()) GetGridDispALE(discretization, lm, edispnp);

  //  CORE::FE::GaussIntegration intpoints(CORE::FE::CellType::line6);

  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------
  for (CORE::FE::GaussIntegration::iterator iquad = my::intpoints_.begin();
       iquad != my::intpoints_.end(); ++iquad)
  {
    // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(), iquad.Weight());

    CORE::LINALG::Matrix<nen_, 1> newfunct(my::funct_);

    XWallTauWIncBack();

    // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(), iquad.Weight());

    CORE::LINALG::Matrix<nen_, 1> oldfunct(my::funct_);

    XWallTauWIncForward();

    //----------------------------------------------------------------------------
    //                         MASS MATRIX
    //----------------------------------------------------------------------------
    std::array<int, nsd_> idim_nsd_p_idim;
    CORE::LINALG::Matrix<nsd_ * nsd_, enren_> lin_resM_Du(true);
    CORE::LINALG::Matrix<enren_ * nsd_, enren_ * nsd_> estif_u(true);

    for (int idim = 0; idim < nsd_; ++idim)
    {
      idim_nsd_p_idim[idim] = idim * nsd_ + idim;
    }

    for (int ui = 1; ui < nen_; ui += 2)
    {
      const double v = my::fac_ * newfunct(ui);
      int nodeui = 0;
      if (ui % 2 == 1)
        nodeui = (ui - 1) / 2;
      else
        FOUR_C_THROW("something wrong with indices. also correct in all following terms");

      for (int idim = 0; idim < nsd_; ++idim)
      {
        lin_resM_Du(idim_nsd_p_idim[idim], nodeui) += v;
      }
    }

    for (int ui = 1; ui < nen_; ui += 2)
    {
      const int nodeui = (ui - 1) / 2;
      const int nsd_nodeui = nsd_ * nodeui;
      for (int vi = 1; vi < nen_; vi += 2)
      {
        const int nsd_nodevi = nsd_ * (vi - 1) / 2;

        for (int idim = 0; idim < nsd_; ++idim)
        {
          estif_u(nsd_nodevi + idim, nsd_nodeui + idim) +=
              newfunct(vi) * lin_resM_Du(idim_nsd_p_idim[idim], nodeui);
        }  // end for (idim)
      }    // vi
    }      // ui

    // add velocity-velocity part to matrix
    for (int ui = 0; ui < enren_; ++ui)
    {
      const int numdof_ui = numdof * ui;
      const int nsd_ui = nsd_ * ui;

      for (int jdim = 0; jdim < nsd_; ++jdim)
      {
        const int numdof_ui_jdim = numdof_ui + jdim;
        const int nsd_ui_jdim = nsd_ui + jdim;

        for (int vi = 0; vi < enren_; ++vi)
        {
          const int numdof_vi = numdof * vi;
          const int nsd_vi = nsd_ * vi;

          for (int idim = 0; idim < nsd_; ++idim)
          {
            elemat1(numdof_vi + idim, numdof_ui_jdim) += estif_u(nsd_vi + idim, nsd_ui_jdim);
          }
        }
      }
    }

    //----------------------------------------------------------------------------
    //                         RHS
    //----------------------------------------------------------------------------
    estif_u.Clear();
    lin_resM_Du.Clear();

    for (int ui = 1; ui < nen_; ui += 2)
    {
      const double v = my::fac_ * (oldfunct(ui) - newfunct(ui));
      const int nodeui = (ui - 1) / 2;

      for (int idim = 0; idim < nsd_; ++idim)
      {
        lin_resM_Du(idim_nsd_p_idim[idim], nodeui) += v;
      }
    }

    for (int ui = 1; ui < nen_; ui += 2)
    {
      const int nodeui = (ui - 1) / 2;
      const int nsd_nodeui = nsd_ * nodeui;
      for (int vi = 1; vi < nen_; vi += 2)
      {
        const int nsd_nodevi = nsd_ * (vi - 1) / 2;
        for (int idim = 0; idim < nsd_; ++idim)
        {
          estif_u(nsd_nodevi + idim, nsd_nodeui + idim) +=
              newfunct(vi) * lin_resM_Du(idim_nsd_p_idim[idim], nodeui);
        }  // end for (idim)
      }    // vi
    }      // ui

    // veln
    // add velocity-velocity part to rhs
    for (int ui = 0; ui < enren_; ++ui)
    {
      const int nsd_ui = nsd_ * ui;
      const int uix = 2 * ui + 1;

      for (int jdim = 0; jdim < nsd_; ++jdim)
      {
        const int nsd_ui_jdim = nsd_ui + jdim;

        for (int vi = 0; vi < enren_; ++vi)
        {
          const int numdof_vi = numdof * vi;
          const int nsd_vi = nsd_ * vi;

          for (int idim = 0; idim < nsd_; ++idim)
          {
            elemat2(numdof_vi + idim, 0) += estif_u(nsd_vi + idim, nsd_ui_jdim) * eveln(jdim, uix);
          }
        }
      }
    }

    // accn
    if (switchonaccn)
    {
      // add velocity-velocity part to rhs
      for (int ui = 0; ui < enren_; ++ui)
      {
        const int nsd_ui = nsd_ * ui;
        const int uix = 2 * ui + 1;

        for (int jdim = 0; jdim < nsd_; ++jdim)
        {
          const int nsd_ui_jdim = nsd_ui + jdim;

          for (int vi = 0; vi < enren_; ++vi)
          {
            const int numdof_vi = numdof * vi;
            const int nsd_vi = nsd_ * vi;

            for (int idim = 0; idim < nsd_; ++idim)
            {
              elemat2(numdof_vi + idim, 1) +=
                  estif_u(nsd_vi + idim, nsd_ui_jdim) * eaccn(jdim, uix);
            }
          }
        }
      }
    }

    if (switchonvelnp)
    {
      // velnp
      // add velocity-velocity part to rhs
      for (int ui = 0; ui < enren_; ++ui)
      {
        const int nsd_ui = nsd_ * ui;
        const int uix = 2 * ui + 1;

        for (int jdim = 0; jdim < nsd_; ++jdim)
        {
          const int nsd_ui_jdim = nsd_ui + jdim;

          for (int vi = 0; vi < enren_; ++vi)
          {
            const int numdof_vi = numdof * vi;
            const int nsd_vi = nsd_ * vi;

            for (int idim = 0; idim < nsd_; ++idim)
            {
              elemat2(numdof_vi + idim, 2) +=
                  estif_u(nsd_vi + idim, nsd_ui_jdim) * evelnp(jdim, uix);
            }
          }
        }
      }
    }
  }


  return 0;
}

/*---------------------------------------------------------------------------*
 | get ALE grid displacements only for element                      bk 02/15 |
 *---------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::GetGridDispALE(
    DRT::Discretization& discretization, const std::vector<int>& lm,
    CORE::LINALG::Matrix<nsd_, nen_>& edispnp)
{
  my::ExtractValuesFromGlobalVector(
      discretization, lm, *my::rotsymmpbc_, &edispnp, nullptr, "dispnp");

  // add displacement when fluid nodes move in the ALE case
  // xyze_ does only know 8 nodes
  // edispnp also knows the virtual ones but doens't do anything with them
  for (unsigned inode = 0; inode < (unsigned)enren_; ++inode)  // number of nodes
    for (int sdm = 0; sdm < nsd_; ++sdm) my::xyze_(sdm, inode) += edispnp(sdm, inode * 2);
}

template <CORE::FE::CellType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::LinMeshMotion_3D(
    CORE::LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& emesh,
    const CORE::LINALG::Matrix<nsd_, nen_>& evelaf, const double& press, const double& timefac,
    const double& timefacfac)
{
  // xGderiv_ = sum(gridx(k,i) * deriv_(j,k), k);
  // xGderiv_ == xjm_
  FOUR_C_THROW("wrong");

  return;
}


/*-----------------------------------------------------------------------------*
 | Prepare custom (direction-dependent) Gauss rule                  bk 06/2014 |
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::PrepareGaussRule()
{
  // which is the wall-normal element direction?
  // calculate jacobian at element center
  my::is_higher_order_ele_ = false;
  my::EvalShapeFuncAndDerivsAtEleCenter();
  my::is_higher_order_ele_ = true;

  // the derivative of the wall distance with respect to the local coordinates
  // shows how the local axes are oriented with respect to the wall-normal vector
  CORE::LINALG::Matrix<nsd_, 1> normwallrst(true);
  normwallrst.Multiply(deriv_, ewdist_);
  double normwallrstnorm2 = normwallrst.Norm2();
  const double dot1 = abs(normwallrst(0) / normwallrstnorm2);
  const double dot2 = abs(normwallrst(1) / normwallrstnorm2);
  const double dot3 = abs(normwallrst(2) / normwallrstnorm2);

  double minyp = 1e9;
  for (int inode = 0; inode < enren_; inode++)
  {
    double utaunode = sqrt(etauw_(inode) * densinv_);
    double yp = ewdist_(inode) * viscinv_ * utaunode;
    if (yp < minyp) minyp = yp;
  }

  if (minyp > 15.0) numgpnorm_ = numgpnormow_;

  if (distype == CORE::FE::CellType::tet4)
  {
    CORE::FE::GaussIntegration intpointsplane(CORE::FE::CellType::tet4, 2 * numgpnorm_ - 1);
    my::intpoints_ = intpointsplane;
  }
  else  // hex8
  {
    cgp_ = Teuchos::rcp(new CORE::FE::CollectedGaussPoints(numgpnorm_ * numgpplane_ * numgpplane_));
    // get the quad9 gaussrule for the in plane integration
    CORE::FE::GaussIntegration intpointsplane(CORE::FE::CellType::quad8, 2 * numgpplane_ - 1);
    // get the quad9 gaussrule for the in normal integration
    CORE::FE::GaussIntegration intpointsnormal(CORE::FE::CellType::line3, 2 * numgpnorm_ - 1);

    // 0.9 corresponds to an angle of 25.8 deg
    if (dot1 < 0.90 && dot2 < 0.90 && dot3 < 0.90)
    {  // element, where the wall normal direction does not point in one specific element direction,
       // e.g. in corners
      cgp_->IncreaseReserved(
          (numgpnorm_ * numgpnorm_ * numgpnorm_) - (numgpnorm_ * numgpplane_ * numgpplane_));
      CORE::FE::GaussIntegration intpointsplane(CORE::FE::CellType::quad8, 2 * numgpnorm_ - 1);
      // start loop over integration points in layer
      for (CORE::FE::GaussIntegration::iterator iquadplane = intpointsplane.begin();
           iquadplane != intpointsplane.end(); ++iquadplane)
      {
        // start loop over integration points in layer
        for (CORE::FE::GaussIntegration::iterator iquadnorm = intpointsnormal.begin();
             iquadnorm != intpointsnormal.end(); ++iquadnorm)
        {
          cgp_->Append(iquadnorm.Point()[0], iquadplane.Point()[0], iquadplane.Point()[1],
              iquadplane.Weight() * iquadnorm.Weight());
        }
      }
    }
    else if (dot1 > dot2 && dot1 > dot3)
    {
      // start loop over integration points in layer
      for (CORE::FE::GaussIntegration::iterator iquadplane = intpointsplane.begin();
           iquadplane != intpointsplane.end(); ++iquadplane)
      {
        // start loop over integration points in layer
        for (CORE::FE::GaussIntegration::iterator iquadnorm = intpointsnormal.begin();
             iquadnorm != intpointsnormal.end(); ++iquadnorm)
        {
          cgp_->Append(iquadnorm.Point()[0], iquadplane.Point()[0], iquadplane.Point()[1],
              iquadplane.Weight() * iquadnorm.Weight());
        }
      }
    }
    else if (dot2 > dot3)
    {
      // start loop over integration points in layer
      for (CORE::FE::GaussIntegration::iterator iquadplane = intpointsplane.begin();
           iquadplane != intpointsplane.end(); ++iquadplane)
      {
        // start loop over integration points in layer
        for (CORE::FE::GaussIntegration::iterator iquadnorm = intpointsnormal.begin();
             iquadnorm != intpointsnormal.end(); ++iquadnorm)
        {
          cgp_->Append(iquadplane.Point()[0], iquadnorm.Point()[0], iquadplane.Point()[1],
              iquadplane.Weight() * iquadnorm.Weight());
        }
      }
    }
    else
    {
      // start loop over integration points in layer
      for (CORE::FE::GaussIntegration::iterator iquadplane = intpointsplane.begin();
           iquadplane != intpointsplane.end(); ++iquadplane)
      {
        // start loop over integration points in layer
        for (CORE::FE::GaussIntegration::iterator iquadnorm = intpointsnormal.begin();
             iquadnorm != intpointsnormal.end(); ++iquadnorm)
        {
          cgp_->Append(iquadplane.Point()[0], iquadplane.Point()[1], iquadnorm.Point()[0],
              iquadplane.Weight() * iquadnorm.Weight());
        }
      }
    }
    CORE::FE::GaussIntegration grule(cgp_);
    my::intpoints_ = grule;
  }

  return;
}

template <CORE::FE::CellType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalcXWall<distype, enrtype>::Sysmat(
    const CORE::LINALG::Matrix<nsd_, nen_>& ebofoaf,
    const CORE::LINALG::Matrix<nsd_, nen_>& eprescpgaf,
    const CORE::LINALG::Matrix<nsd_, nen_>& ebofon,
    const CORE::LINALG::Matrix<nsd_, nen_>& eprescpgn,
    const CORE::LINALG::Matrix<nsd_, nen_>& evelaf, const CORE::LINALG::Matrix<nsd_, nen_>& evelam,
    const CORE::LINALG::Matrix<nsd_, nen_>& eveln, const CORE::LINALG::Matrix<nsd_, nen_>& evelnp,
    const CORE::LINALG::Matrix<nsd_, nen_>& fsevelaf, const CORE::LINALG::Matrix<nen_, 1>& fsescaaf,
    const CORE::LINALG::Matrix<nsd_, nen_>& evel_hat,
    const CORE::LINALG::Matrix<nsd_ * nsd_, nen_>& ereynoldsstress_hat,
    const CORE::LINALG::Matrix<nen_, 1>& epreaf, const CORE::LINALG::Matrix<nen_, 1>& epream,
    const CORE::LINALG::Matrix<nen_, 1>& epren, const CORE::LINALG::Matrix<nen_, 1>& eprenp,
    const CORE::LINALG::Matrix<nsd_, nen_>& eaccam, const CORE::LINALG::Matrix<nen_, 1>& escaaf,
    const CORE::LINALG::Matrix<nen_, 1>& escaam, const CORE::LINALG::Matrix<nen_, 1>& escadtam,
    const CORE::LINALG::Matrix<nsd_, nen_>& eveldtam, const CORE::LINALG::Matrix<nen_, 1>& epredtam,
    const CORE::LINALG::Matrix<nen_, 1>& escabofoaf, const CORE::LINALG::Matrix<nen_, 1>& escabofon,
    const CORE::LINALG::Matrix<nsd_, nen_>& emhist, const CORE::LINALG::Matrix<nsd_, nen_>& edispnp,
    const CORE::LINALG::Matrix<nsd_, nen_>& egridv,
    CORE::LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& estif,
    CORE::LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& emesh,
    CORE::LINALG::Matrix<(nsd_ + 1) * nen_, 1>& eforce, const CORE::LINALG::Matrix<nen_, 1>& eporo,
    const CORE::LINALG::Matrix<nsd_, 2 * nen_>& egradphi,
    const CORE::LINALG::Matrix<nen_, 2 * 1>& ecurvature, const double thermpressaf,
    const double thermpressam, const double thermpressdtaf, const double thermpressdtam,
    Teuchos::RCP<const CORE::MAT::Material> material, double& Cs_delta_sq, double& Ci_delta_sq,
    double& Cv, bool isale, double* saccn, double* sveln, double* svelnp,
    const CORE::FE::GaussIntegration& intpoints)
{
  my::Sysmat(ebofoaf, eprescpgaf, ebofon, eprescpgn, evelaf, evelam, eveln, evelnp, fsevelaf,
      fsescaaf, evel_hat, ereynoldsstress_hat, epreaf, epream, epren, eprenp, eaccam, escaaf,
      escaam, escadtam, eveldtam, epredtam, escabofoaf, escabofon, emhist, edispnp, egridv, estif,
      emesh,  // -> emesh
      eforce, eporo, egradphi, ecurvature, thermpressaf, thermpressam, thermpressdtaf,
      thermpressdtam, material, Cs_delta_sq, Ci_delta_sq, Cv, isale, saccn, sveln, svelnp,
      intpoints);
  return;
}

template class DRT::ELEMENTS::FluidEleCalcXWall<CORE::FE::CellType::tet4,
    DRT::ELEMENTS::Fluid::xwall>;
template class DRT::ELEMENTS::FluidEleCalcXWall<CORE::FE::CellType::hex8,
    DRT::ELEMENTS::Fluid::xwall>;

FOUR_C_NAMESPACE_CLOSE
