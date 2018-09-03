/*--------------------------------------------------------------------------*/
/*!
\file lubrication_ele_calc.cpp

\brief main file containing routines for calculation of lubrication element

\level 3

\maintainer Alexander Seitz

*/
/*--------------------------------------------------------------------------*/

#include "../drt_lubrication_ele/lubrication_ele_calc.H"

#include "../drt_geometry/position_array.H"

#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_lubrication_ele/lubrication_ele_parameter.H"

#include "../drt_inpar/inpar_lubrication.H"

#include "../drt_lubrication_ele/lubrication_ele_calc_utils.H"
#include "../drt_lubrication_ele/lubrication_ele_action.H"

#include "../drt_lib/standardtypes_cpp.H"

#include "../drt_mat/lubrication_mat.H"


/*----------------------------------------------------------------------*
 * Constructor
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::LubricationEleCalc(const std::string& disname)
    : lubricationpara_(
          DRT::ELEMENTS::LubricationEleParameter::Instance(disname)),  // standard parameter list
      eprenp_(true),                                                   // initialized to zero
      xsi_(true),                                                      // initialized to zero
      xyze_(true),                                                     // initialized to zero
      funct_(true),                                                    // initialized to zero
      deriv_(true),                                                    // initialized to zero
      derxy_(true),                                                    // initialized to zero
      xjm_(true),                                                      // initialized to zero
      xij_(true),                                                      // initialized to zero
      eheinp_(true),
      eheidotnp_(true),
      edispnp_(true),
      viscmanager_(
          Teuchos::rcp(new LubricationEleViscManager())),  // viscosity manager for viscosity
      lubricationvarmanager_(Teuchos::rcp(
          new LubricationEleInternalVariableManager<nsd_, nen_>())),  // internal variable manager
      eid_(0),
      ele_(NULL),
      Dt_(0.0)
{
  dsassert(
      nsd_ >= nsd_ele_, "problem dimension has to be equal or larger than the element dimension!");

  return;
}

/*----------------------------------------------------------------------*
 | singleton access method                                  wirtz 10/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::LubricationEleCalc<distype, probdim>*
DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::Instance(
    const std::string& disname, const LubricationEleCalc* delete_me)
{
  static std::map<std::string, LubricationEleCalc<distype, probdim>*> instances;

  if (delete_me == NULL)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] = new LubricationEleCalc<distype, probdim>(disname);
  }
  else
  {
    // since we keep several instances around in the general case, we need to
    // find which of the instances to delete with this call. This is done by
    // letting the object to be deleted hand over the 'this' pointer, which is
    // located in the map and deleted
    for (typename std::map<std::string, LubricationEleCalc<distype, probdim>*>::iterator i =
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
 | singleton destruction                                    wirtz 10/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::Done()
{
  // delete instance
  Instance("", this);

  return;
}

/*----------------------------------------------------------------------*
 * Action type: Evaluate                                    wirtz 10/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
int DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::Evaluate(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& elemat1_epetra,
    Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
    Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra)
{
  //--------------------------------------------------------------------------------
  // preparations for element
  //--------------------------------------------------------------------------------

  if (SetupCalc(ele, discretization) == -1) return 0;

  //--------------------------------------------------------------------------------
  // extract element based or nodal values
  //--------------------------------------------------------------------------------

  ExtractElementAndNodeValues(ele, params, discretization, la);

  //--------------------------------------------------------------------------------
  // calculate element coefficient matrix and rhs
  //--------------------------------------------------------------------------------

  Sysmat(ele, elemat1_epetra, elevec1_epetra);

  return 0;
}

/*----------------------------------------------------------------------*
 * Action type: Evaluate element matrix associated to the linearization
 * of the residual wrt. the discrete film height (only used for monolithic
 * EHL problems)
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
int DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::EvaluateEHLMon(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& elemat1_epetra,
    Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
    Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra)
{
  //--------------------------------------------------------------------------------
  // preparations for element
  //--------------------------------------------------------------------------------

  if (SetupCalc(ele, discretization) == -1) return 0;

  // set time step as a class member
  double dt = params.get<double>("delta time");
  Dt_ = dt;

  //--------------------------------------------------------------------------------
  // extract element based or nodal values
  //--------------------------------------------------------------------------------

  ExtractElementAndNodeValues(ele, params, discretization, la);

  //--------------------------------------------------------------------------------
  // calculate element off-diagonal-matrix for height linearization in monolithic EHL
  //--------------------------------------------------------------------------------

  MatrixforEHLMon(ele, elemat1_epetra, elemat2_epetra);

  return 0;
}

/*----------------------------------------------------------------------*
 | setup element evaluation                                 wirtz 10/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
int DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::SetupCalc(
    DRT::Element* ele, DRT::Discretization& discretization)
{
  // get element coordinates
  ReadElementCoordinates(ele);

  // set element id
  eid_ = ele->Id();
  // set element
  ele_ = ele;

  return 0;
}

/*----------------------------------------------------------------------*
 | read element coordinates                                 wirtz 10/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::ReadElementCoordinates(
    const DRT::Element* ele)
{
  // Directly copy the coordinates since in 3D the transformation is just the identity
  GEO::fillInitialPositionArray<distype, nsd_, LINALG::Matrix<nsd_, nen_>>(ele, xyze_);

  return;
}  // LubricationEleCalc::ReadElementCoordinates

/*----------------------------------------------------------------------*
 | extract element based or nodal values
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::ExtractElementAndNodeValues(
    DRT::Element* ele, Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la)
{
  // 1. Extract the tangential velocity

  // nodeset the velocity is defined on
  const int ndsvel = 1;

  // get the global vector
  Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState(ndsvel, "av_tang_vel");
  if (vel.is_null()) dserror("got NULL pointer for \"av_tang_vel\"");

  const int numveldofpernode = la[ndsvel].lm_.size() / nen_;

  // construct location vector for velocity related dofs
  std::vector<int> lmvel(nsd_ * nen_, -1);
  for (int inode = 0; inode < nen_; ++inode)
    for (int idim = 0; idim < nsd_; ++idim)
      lmvel[inode * nsd_ + idim] = la[ndsvel].lm_[inode * numveldofpernode + idim];

  DRT::UTILS::ExtractMyValues<LINALG::Matrix<nsd_, nen_>>(*vel, eAvTangVel_, lmvel);

  // 2. In case of ale, extract the displacements of the element nodes and update the nodal
  // coordinates

  // Only required, in case of ale:
  if (params.get<bool>("isale"))
  {
    // get number of dofset associated with displacement related dofs
    const int ndsdisp = 1;  // needs further implementation: params.get<int>("ndsdisp");
    Teuchos::RCP<const Epetra_Vector> dispnp = discretization.GetState(ndsdisp, "dispnp");
    if (dispnp == Teuchos::null) dserror("Cannot get state vector 'dispnp'");

    // determine number of displacement related dofs per node
    const int numdispdofpernode = la[ndsdisp].lm_.size() / nen_;

    // construct location vector for displacement related dofs
    std::vector<int> lmdisp(nsd_ * nen_, -1);
    for (int inode = 0; inode < nen_; ++inode)
      for (int idim = 0; idim < nsd_; ++idim)
        lmdisp[inode * nsd_ + idim] = la[ndsdisp].lm_[inode * numdispdofpernode + idim];

    // extract local values of displacement field from global state vector
    DRT::UTILS::ExtractMyValues<LINALG::Matrix<nsd_, nen_>>(*dispnp, edispnp_, lmdisp);

    // add nodal displacements to point coordinates
    xyze_ += edispnp_;
  }

  // 3. Extract the film height at the element nodes

  // get number of dofset associated with height dofs
  const int ndsheight = 1;  // needs further implementation: params.get<int>("ndsheight");

  // get the global vector containing the heights
  Teuchos::RCP<const Epetra_Vector> height = discretization.GetState(ndsheight, "height");
  if (height == Teuchos::null) dserror("Cannot get state vector height");

  // determine number of height related dofs per node
  const int numheightdofpernode = la[ndsheight].lm_.size() / nen_;

  // construct location vector for height related dofs
  std::vector<int> lmheight(nsd_ * nen_, -1);
  for (int inode = 0; inode < nen_; ++inode)
    for (int idim = 0; idim < nsd_; ++idim)
      lmheight[inode * nsd_ + idim] = la[ndsheight].lm_[inode * numheightdofpernode + idim];

  // extract local height from global state vector
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<nsd_, nen_>>(*height, eheinp_, la[ndsheight].lm_);

  // 3.1. Extract the film height time derivative at the element node
  const int ndsheightdot = 1;

  // get the global vector containing the heightdots
  Teuchos::RCP<const Epetra_Vector> heightdot = discretization.GetState(ndsheightdot, "heightdot");
  if (heightdot == Teuchos::null) dserror("Cannot get state vector heightdot");

  // determine number of heightdot related dofs per node
  const int numheightdotdofpernode = la[ndsheightdot].lm_.size() / nen_;

  // construct location vector for heightdot related dofs
  std::vector<int> lmheightdot(nsd_ * nen_, -1);
  for (int inode = 0; inode < nen_; ++inode)
    for (int idim = 0; idim < nsd_; ++idim)
      lmheightdot[inode * nsd_ + idim] =
          la[ndsheightdot].lm_[inode * numheightdotdofpernode + idim];

  // extract local height from global state vector
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<nsd_, nen_>>(
      *heightdot, eheidotnp_, la[ndsheightdot].lm_);

  // 4. Extract the pressure field at the element nodes

  // get the global vector containing the pressure
  Teuchos::RCP<const Epetra_Vector> prenp = discretization.GetState("prenp");
  if (prenp == Teuchos::null) dserror("Cannot get state vector 'prenp'");

  // values of pressure field are always in first dofset
  const std::vector<int>& lm = la[0].lm_;

  // extract the local values at the element nodes
  DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_, 1>>(*prenp, eprenp_, lm);

  return;
}

/*----------------------------------------------------------------------*
|  calculate system matrix and rhs (public)                 wirtz 10/15 |
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::Sysmat(
    DRT::Element* ele,               ///< the element whose matrix is calculated
    Epetra_SerialDenseMatrix& emat,  ///< element matrix to calculate
    Epetra_SerialDenseVector& erhs   ///< element rhs to calculate
)
{
  //----------------------------------------------------------------------
  // get material parameters (evaluation at element center)
  //----------------------------------------------------------------------
  // density at t_(n)
  double densn(1.0);
  // density at t_(n+1) or t_(n+alpha_F)
  double densnp(1.0);
  // density at t_(n+alpha_M)
  double densam(1.0);

  // fluid viscosity
  double visc(0.0);

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      LUBRICATION::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    // set gauss point variables needed for evaluation of mat and rhs
    SetInternalVariablesForMatAndRHS();

    // calculate height (i.e. the distance of the contacting bodies) at Integration point
    double heightint(0.0);
    CalcHeightAtIntPoint(heightint);

    // calculate heightDot (i.e. the distance of the contacting bodies) at Integration point
    double heightdotint(0.0);
    CalcHeightDotAtIntPoint(heightdotint);

    // calculate average surface velocity of the contacting bodies at Integration point
    LINALG::Matrix<nsd_, 1> avrvel(true);  // average surface velocity, initialized to zero
    CalcAvrVelAtIntPoint(avrvel);

    //----------------------------------------------------------------------
    // get material parameters (evaluation at integration point)
    //----------------------------------------------------------------------
    GetMaterialParams(ele, densn, densnp, densam, visc, iquad);

    // integration factors
    const double timefacfac = fac;  // works only for stationary problems!

    double rhsfac = fac;  // works only for stationary problems!

    // 1) element matrix

    // 1.1) calculation of Poiseuille contribution of element matrix

    CalcMatPsl(emat, timefacfac, visc, heightint);

    // 2) rhs matrix

    // 2.1) calculation of Poiseuille contribution of rhs matrix

    CalcRhsPsl(erhs, rhsfac, visc, heightint);

    // 2.2) calculation of Wedge contribution of rhs matrix

    CalcRhsWdg(erhs, rhsfac, heightint, avrvel);

    // 2.3) calculation of squeeze contribution to RHS matrix
    CalcRhsSqz(erhs, rhsfac, heightdotint);

  }  // end loop Gauss points

  return;
}

template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::MatrixforEHLMon(
    DRT::Element* ele,  ///< the element whose matrix is calculated
    Epetra_SerialDenseMatrix& ematheight, Epetra_SerialDenseMatrix& ematvel)
{
  //----------------------------------------------------------------------
  // get material parameters (evaluation at element center)
  //----------------------------------------------------------------------
  // density at t_(n)
  double densn(1.0);
  // density at t_(n+1) or t_(n+alpha_F)
  double densnp(1.0);
  // density at t_(n+alpha_M)
  double densam(1.0);

  // fluid viscosity
  double visc(0.0);

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      LUBRICATION::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    // set gauss point variables needed for evaluation of mat and rhs
    SetInternalVariablesForMatAndRHS();

    // calculate height (i.e. the distance of the contacting bodies) at Integration point
    double heightint(0.0);
    CalcHeightAtIntPoint(heightint);

    // calculate average surface velocity of the contacting bodies at Integration point
    LINALG::Matrix<nsd_, 1> avrvel(true);  // average surface velocity, initialized to zero
    CalcAvrVelAtIntPoint(avrvel);

    //----------------------------------------------------------------------
    // get material parameters (evaluation at integration point)
    //----------------------------------------------------------------------
    GetMaterialParams(ele, densn, densnp, densam, visc, iquad);

    const LINALG::Matrix<nsd_, 1>& gradpre = lubricationvarmanager_->GradPre();

    // Linearization of Poiseuille term wrt the film height
    for (int vi = 0; vi < nen_; vi++)
    {
      double laplawf(0.0);
      GetLaplacianWeakFormRHS(laplawf, gradpre, vi);
      for (int ui = 0; ui < nen_; ui++)
      {
        double val = fac * (1 / (12 * visc)) * 3 * heightint * heightint * laplawf * funct_(ui);
        ematheight(vi, (ui * nsd_)) -= val;
      }
    }  // end loop for linearization of Poiseuille term wrt the film height

    // Linearization of Couette term wrt the film height
    for (int vi = 0; vi < nen_; vi++)
    {
      double val(0.0);

      for (int idim = 0; idim < nsd_; idim++)
      {
        val += derxy_(idim, vi) * avrvel(idim);
      }

      for (int ui = 0; ui < nen_; ui++)
      {
        ematheight(vi, (ui * nsd_)) += fac * val * funct_(ui);
      }
    }  // end loop for linearization of Couette term wrt the film height

    // Linearization of Squeeze term wrt the film height
    for (int vi = 0; vi < nen_; vi++)
    {
      for (int ui = 0; ui < nen_; ui++)
      {
        ematheight(vi, (ui * nsd_)) -= fac * (1.0 / Dt_) * funct_(ui) * funct_(vi);
      }
    }  // end loop for linearization of Squeeze term wrt the film height

    // Linearization of Couette term wrt the velocities
    for (int vi = 0; vi < nen_; vi++)
    {
      for (int ui = 0; ui < nen_; ui++)
      {
        for (int idim = 0; idim < nsd_; idim++)
        {
          ematvel(vi, (ui * nsd_ + idim)) -= fac * heightint * derxy_(idim, vi) * funct_(ui);
        }
      }
    }  // end loop for linearization of Couette term wrt the velocities

  }  // end gauß point loop

  return;
}


/*------------------------------------------------------------------------------*
 | set internal variables                                           wirtz 10/15 |
 *------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::SetInternalVariablesForMatAndRHS()
{
  lubricationvarmanager_->SetInternalVariables(funct_, derxy_, eprenp_);
  return;
}

/*----------------------------------------------------------------------*
 |  get the material constants  (private)                     wirtz 10/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::GetMaterialParams(
    const DRT::Element* ele,  //!< the element we are dealing with
    double& densn,            //!< density at t_(n)
    double& densnp,           //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,           //!< density at t_(n+alpha_M)
    double& visc,             //!< fluid viscosity
    const int iquad           //!< id of current gauss point
)
{
  // get the material
  Teuchos::RCP<MAT::Material> material = ele->Material();

  Materials(material, densn, densnp, densam, visc, iquad);

  return;
}  // LubricationEleCalc::GetMaterialParams

/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                   wirtz 10/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::Materials(
    const Teuchos::RCP<MAT::Material> material,  //!< pointer to current material
    double& densn,                               //!< density at t_(n)
    double& densnp,                              //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,                              //!< density at t_(n+alpha_M)
    double& visc,                                //!< fluid viscosity
    const int iquad                              //!< id of current gauss point

)
{
  switch (material->MaterialType())
  {
    case INPAR::MAT::m_lubrication:
      MatLubrication(material, densn, densnp, densam, visc, iquad);
      break;
    default:
      dserror("Material type %i is not supported", material->MaterialType());
      break;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Material Lubrication                                    wirtz 10/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::MatLubrication(
    const Teuchos::RCP<MAT::Material> material,  //!< pointer to current material
    double& densn,                               //!< density at t_(n)
    double& densnp,                              //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,                              //!< density at t_(n+alpha_M)
    double& visc,                                //!< fluid viscosity
    const int iquad                              //!< id of current gauss point (default = -1)
)
{
  const Teuchos::RCP<MAT::LubricationMat>& actmat =
      Teuchos::rcp_dynamic_cast<MAT::LubricationMat>(material);
  // get constant viscosity

  double pressure = 0.0;
  //  const double pres = my::eprenp_.Dot(my::funct_);
  //  const double pre = lubricationvarmanager_->Prenp();
  //  const double p = eprenp_.Dot(funct_);

  visc = actmat->ComputeViscosity(pressure);

  viscmanager_->SetIsotropicVisc(visc);
  return;
}  // LubricationEleCalc<distype>::MatLubrication

/*------------------------------------------------------------------- *
 |  calculation of Poiseuille element matrix
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::CalcMatPsl(Epetra_SerialDenseMatrix& emat,
    const double timefacfac, const double viscosity, const double height)
{
  // Poiseuille term
  const double fac_psl = timefacfac * (1 / (12 * viscosity)) * height * height * height;
  for (int vi = 0; vi < nen_; ++vi)
  {
    for (int ui = 0; ui < nen_; ++ui)
    {
      double laplawf(0.0);
      GetLaplacianWeakForm(laplawf, ui, vi);
      emat(vi, ui) -= fac_psl * laplawf;
    }
  }
  return;
}

/*------------------------------------------------------------------- *
 |  calculation of Poiseuille rhs matrix
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::CalcRhsPsl(Epetra_SerialDenseVector& erhs,
    const double rhsfac, const double viscosity, const double height)
{
  // Poiseuille rhs term
  const double fac_rhs_psl = rhsfac * (1 / (12 * viscosity)) * height * height * height;

  const LINALG::Matrix<nsd_, 1>& gradpre = lubricationvarmanager_->GradPre();

  for (int vi = 0; vi < nen_; ++vi)
  {
    double laplawf(0.0);
    GetLaplacianWeakFormRHS(laplawf, gradpre, vi);
    erhs[vi] += fac_rhs_psl * laplawf;
  }
  return;
}

/*------------------------------------------------------------------- *
 |  calculation of Wedge rhs matrix
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::CalcRhsWdg(Epetra_SerialDenseVector& erhs,
    const double rhsfac, const double height, const LINALG::Matrix<nsd_, 1> velocity)
{
  // Wedge rhs term
  const double fac_rhs_wdg = rhsfac * height;

  for (int vi = 0; vi < nen_; ++vi)
  {
    double val(0.0);

    for (int i = 0; i < nsd_; ++i)
    {
      val += derxy_(i, vi) * velocity(i);
    }
    erhs[vi] -= fac_rhs_wdg * val;
  }
  return;
}

/*------------------------------------------------------------------- *
 |  calculation of Squeeze rhs matrix                         Faraji  |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::CalcRhsSqz(
    Epetra_SerialDenseVector& erhs, const double rhsfac, const double heightdot)
{
  // Squeeze rhs term
  const double fac_rhs_sqz = rhsfac * heightdot;

  for (int vi = 0; vi < nen_; ++vi)
  {
    erhs[vi] += fac_rhs_sqz * funct_(vi);
  }
  return;
}

/*----------------------------------------------------------------------*
 | evaluate shape functions and derivatives at int. point   wirtz 10/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
double DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::EvalShapeFuncAndDerivsAtIntPoint(
    const DRT::UTILS::IntPointsAndWeights<nsd_ele_>& intpoints,  ///< integration points
    const int iquad                                              ///< id of current Gauss point
)
{
  // coordinates of the current integration point
  const double* gpcoord = (intpoints.IP().qxg)[iquad];
  for (int idim = 0; idim < nsd_ele_; idim++) xsi_(idim) = gpcoord[idim];

  const double det = EvalShapeFuncAndDerivsInParameterSpace();

  if (det < 1E-16)
    dserror("GLOBAL ELEMENT NO. %d \nZERO OR NEGATIVE JACOBIAN DETERMINANT: %lf", eid_, det);

  // compute global spatial derivatives
  derxy_.Multiply(xij_, deriv_);

  // set integration factor: fac = Gauss weight * det(J)
  const double fac = intpoints.IP().qwgt[iquad] * det;

  // return integration factor for current GP: fac = Gauss weight * det(J)
  return fac;

}  // LubricationImpl::EvalShapeFuncAndDerivsAtIntPoint

/*----------------------------------------------------------------------*
 | evaluate shape functions and derivatives in parameter space wirtz 10/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
double DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::EvalShapeFuncAndDerivsInParameterSpace()
{
  double det = 0.0;

  if (nsd_ == nsd_ele_)  // standard case
  {
    // shape functions and their first derivatives
    DRT::UTILS::shape_function<distype>(xsi_, funct_);
    DRT::UTILS::shape_function_deriv1<distype>(xsi_, deriv_);


    // compute Jacobian matrix and determinant
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

    xjm_.MultiplyNT(deriv_, xyze_);
    det = xij_.Invert(xjm_);
  }
  else  // element dimension is smaller than problem dimension -> mannifold
  {
    static LINALG::Matrix<nsd_ele_, nen_> deriv_red;

    // shape functions and their first derivatives
    DRT::UTILS::shape_function<distype>(xsi_, funct_);
    DRT::UTILS::shape_function_deriv1<distype>(xsi_, deriv_red);

    //! metric tensor at integration point
    static LINALG::Matrix<nsd_ele_, nsd_ele_> metrictensor;
    static LINALG::Matrix<nsd_, 1> normalvec;

    // the metric tensor and the area of an infinitesimal surface/line element
    // optional: get unit normal at integration point as well
    DRT::UTILS::ComputeMetricTensorForBoundaryEle<distype, nsd_>(
        xyze_, deriv_red, metrictensor, det, &normalvec);

    if (det < 1E-16)
      dserror("GLOBAL ELEMENT NO. %d \nZERO OR NEGATIVE JACOBIAN DETERMINANT: %lf", eid_, det);

    // transform the derivatives and Jacobians to the higher dimensional coordinates(problem
    // dimension)
    static LINALG::Matrix<nsd_ele_, nsd_> xjm_red;
    xjm_red.MultiplyNT(deriv_red, xyze_);

    for (int i = 0; i < nsd_; i++)
    {
      for (int j = 0; j < nsd_ele_; j++) xjm_(j, i) = xjm_red(j, i);
      xjm_(nsd_ele_, i) = normalvec(i, 0);
    }

    for (int i = 0; i < nen_; i++)
    {
      for (int j = 0; j < nsd_ele_; j++) deriv_(j, i) = deriv_red(j, i);
      deriv_(nsd_ele_, i) = 0.0;
    }

    // special case: 1D element embedded in 3D problem
    if (nsd_ele_ == 1 and nsd_ == 3)
    {
      // compute second unit normal
      const double normalvec2_0 = xjm_red(0, 1) * normalvec(2, 0) - normalvec(1, 0) * xjm_red(0, 2);
      const double normalvec2_1 = xjm_red(0, 2) * normalvec(0, 0) - normalvec(2, 0) * xjm_red(0, 0);
      const double normalvec2_2 = xjm_red(0, 0) * normalvec(1, 0) - normalvec(0, 0) * xjm_red(0, 1);

      // norm
      const double norm2 = std::sqrt(
          normalvec2_0 * normalvec2_0 + normalvec2_1 * normalvec2_1 + normalvec2_2 * normalvec2_2);

      xjm_(2, 0) = normalvec2_0 / norm2;
      xjm_(2, 1) = normalvec2_1 / norm2;
      xjm_(2, 2) = normalvec2_2 / norm2;

      for (int i = 0; i < nen_; i++) deriv_(2, i) = 0.0;
    }

    xij_.Invert(xjm_);
  }

  return det;
}

/*-----------------------------------------------------------------------*
  |  get the lubrication height interpolated at the Int Point
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::CalcHeightAtIntPoint(
    double& heightint  //!< lubrication height at Int point
)
{
  // interpolate the height at the integration point
  for (int j = 0; j < nen_; j++)
  {
    heightint +=
        funct_(j) * eheinp_(0, j);  // Note that the same value is stored for all space dimensions
  }

  return;
}  // ReynoldsEleCalc::CalcHeightAtIntPoint

/*-----------------------------------------------------------------------*
  |  get the lubrication heightDot interpolated at the Int Point    Faraji|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::CalcHeightDotAtIntPoint(
    double& heightdotint  //!< lubrication heightDot at Int point
)
{
  // interpolate the heightDot at the integration point
  for (int j = 0; j < nen_; j++)
  {
    heightdotint +=
        funct_(j) *
        eheidotnp_(0, j);  // Note that the same value is stored for all space dimensions
  }

  return;
}

/*-----------------------------------------------------------------------*
  |  get the average velocity of the contacting bodies interpolated at   |
  |  the Int Point
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::CalcAvrVelAtIntPoint(
    LINALG::Matrix<nsd_, 1>& avrvel  //!< average surface velocity at Int point
)
{
  // interpolate the velocities at the integration point
  avrvel.Multiply(1., eAvTangVel_, funct_, 0.);

  return;
}


/*----------------------------------------------------------------------*
 | evaluate service routine                                 wirtz 10/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
int DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::EvaluateService(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& elemat1_epetra,
    Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
    Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra)
{
  // setup
  if (SetupCalc(ele, discretization) == -1) return 0;

  // check for the action parameter
  const LUBRICATION::Action action = DRT::INPUT::get<LUBRICATION::Action>(params, "action");

  // evaluate action
  EvaluateAction(ele, params, discretization, action, la, elemat1_epetra, elemat2_epetra,
      elevec1_epetra, elevec2_epetra, elevec3_epetra);

  return 0;
}

/*----------------------------------------------------------------------*
 | evaluate action                                          wirtz 10/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
int DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::EvaluateAction(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    const LUBRICATION::Action& action, DRT::Element::LocationArray& la,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra)
{
  //(for now) only first dof set considered
  const std::vector<int>& lm = la[0].lm_;
  // determine and evaluate action
  switch (action)
  {
    case LUBRICATION::calc_error:
    {
      // check if length suffices
      if (elevec1_epetra.Length() < 1) dserror("Result vector too short");

      // need current solution
      Teuchos::RCP<const Epetra_Vector> prenp = discretization.GetState("prenp");
      if (prenp == Teuchos::null) dserror("Cannot get state vector 'prenp'");
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_, 1>>(*prenp, eprenp_, lm);

      CalErrorComparedToAnalytSolution(ele, params, elevec1_epetra);

      break;
    }

    case LUBRICATION::calc_mean_pressures:
    {
      // get flag for inverting
      bool inverting = params.get<bool>("inverting");

      // need current pressure vector
      // -> extract local values from the global vectors
      Teuchos::RCP<const Epetra_Vector> prenp = discretization.GetState("prenp");
      if (prenp == Teuchos::null) dserror("Cannot get state vector 'prenp'");
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<nen_, 1>>(*prenp, eprenp_, lm);

      // calculate pressures and domain integral
      CalculatePressures(ele, elevec1_epetra, inverting);

      break;
    }

    default:
    {
      dserror("Not acting on this action. Forgot implementation?");
      break;
    }
  }  // switch(action)

  return 0;
}

/*----------------------------------------------------------------------*
  |  calculate error compared to analytical solution        wirtz 10/15 |
  *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::CalErrorComparedToAnalytSolution(
    const DRT::Element* ele, Teuchos::ParameterList& params, Epetra_SerialDenseVector& errors)
{
  if (DRT::INPUT::get<LUBRICATION::Action>(params, "action") != LUBRICATION::calc_error)
    dserror("How did you get here?");

  // -------------- prepare common things first ! -----------------------
  // set constants for analytical solution
  const double t = lubricationpara_->Time();

  // integration points and weights
  // more GP than usual due to (possible) cos/exp fcts in analytical solutions
  const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      LUBRICATION::DisTypeToGaussRuleForExactSol<distype>::rule);

  const INPAR::LUBRICATION::CalcError errortype =
      DRT::INPUT::get<INPAR::LUBRICATION::CalcError>(params, "calcerrorflag");
  switch (errortype)
  {
    case INPAR::LUBRICATION::calcerror_byfunction:
    {
      const int errorfunctno = params.get<int>("error function number");

      // analytical solution
      double pre_exact(0.0);
      double deltapre(0.0);
      //! spatial gradient of current pressure value
      LINALG::Matrix<nsd_, 1> gradpre(true);
      LINALG::Matrix<nsd_, 1> gradpre_exact(true);
      LINALG::Matrix<nsd_, 1> deltagradpre(true);

      // start loop over integration points
      for (int iquad = 0; iquad < intpoints.IP().nquad; iquad++)
      {
        const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

        // get coordinates at integration point
        // gp reference coordinates
        LINALG::Matrix<nsd_, 1> xyzint(true);
        xyzint.Multiply(xyze_, funct_);

        // function evaluation requires a 3D position vector!!
        double position[3] = {0.0, 0.0, 0.0};

        for (int dim = 0; dim < nsd_; ++dim) position[dim] = xyzint(dim);

        // pressure at integration point at time step n+1
        const double prenp = funct_.Dot(eprenp_);
        // spatial gradient of current pressure value
        gradpre.Multiply(derxy_, eprenp_);

        pre_exact = DRT::Problem::Instance()->Funct(errorfunctno - 1).Evaluate(0, position, t);

        std::vector<double> gradpre_exact_vec = DRT::Problem::Instance()
                                                    ->Funct(errorfunctno - 1)
                                                    .EvaluateSpatialDerivative(0, position, t);

        if (gradpre_exact_vec.size())
        {
          if (nsd_ == nsd_ele_)
            for (int dim = 0; dim < nsd_; ++dim) gradpre_exact(dim) = gradpre_exact_vec[dim];
          else
          {
            std::cout << "Warning: Gradient of analytical solution cannot be evaluated correctly "
                         "for lubrication on curved surfaces!"
                      << std::endl;
            gradpre_exact.Clear();
          }
        }
        else
        {
          std::cout << "Warning: Gradient of analytical solution was not evaluated!" << std::endl;
          gradpre_exact.Clear();
        }

        // error at gauss point
        deltapre = prenp - pre_exact;
        deltagradpre.Update(1.0, gradpre, -1.0, gradpre_exact);

        // 0: delta pressure for L2-error norm
        // 1: delta pressure for H1-error norm
        // 2: analytical pressure for L2 norm
        // 3: analytical pressure for H1 norm

        // the error for the L2 and H1 norms are evaluated at the Gauss point

        // integrate delta pressure for L2-error norm
        errors(0) += deltapre * deltapre * fac;
        // integrate delta pressure for H1-error norm
        errors(1) += deltapre * deltapre * fac;
        // integrate analytical pressure for L2 norm
        errors(2) += pre_exact * pre_exact * fac;
        // integrate analytical pressure for H1 norm
        errors(3) += pre_exact * pre_exact * fac;

        // integrate delta pressure derivative for H1-error norm
        errors(1) += deltagradpre.Dot(deltagradpre) * fac;
        // integrate analytical pressure derivative for H1 norm
        errors(3) += gradpre_exact.Dot(gradpre_exact) * fac;
      }  // loop over integration points
    }
    break;
    default:
      dserror("Unknown analytical solution!");
      break;
  }  // switch(errortype)

  return;
}  // DRT::ELEMENTS::LubricationEleCalc<distype,probdim>::CalErrorComparedToAnalytSolution

/*---------------------------------------------------------------------*
|  calculate pressure(s) and domain integral               wirtz 10/15 |
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::LubricationEleCalc<distype, probdim>::CalculatePressures(
    const DRT::Element* ele, Epetra_SerialDenseVector& pressures, const bool inverting)
{
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      LUBRICATION::DisTypeToOptGaussRule<distype>::rule);

  // integration loop
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    // calculate integrals of (inverted) pressure(s) and domain
    if (inverting)
    {
      for (int i = 0; i < nen_; i++)
      {
        const double fac_funct_i = fac * funct_(i);
        if (std::abs(eprenp_(i, 0)) > EPS14)
          pressures[0] += fac_funct_i / eprenp_(i, 0);
        else
          dserror("Division by zero");
        // for domain volume
        pressures[1] += fac_funct_i;
      }
    }
    else
    {
      for (int i = 0; i < nen_; i++)
      {
        const double fac_funct_i = fac * funct_(i);
        pressures[0] += fac_funct_i * eprenp_(i, 0);
        // for domain volume
        pressures[1] += fac_funct_i;
      }
    }
  }  // loop over integration points

  return;
}  // LubricationEleCalc::CalculatePressures

// template classes

// 1D elements
template class DRT::ELEMENTS::LubricationEleCalc<DRT::Element::line2, 1>;
template class DRT::ELEMENTS::LubricationEleCalc<DRT::Element::line2, 2>;
template class DRT::ELEMENTS::LubricationEleCalc<DRT::Element::line2, 3>;
template class DRT::ELEMENTS::LubricationEleCalc<DRT::Element::line3, 1>;
template class DRT::ELEMENTS::LubricationEleCalc<DRT::Element::line3, 2>;
template class DRT::ELEMENTS::LubricationEleCalc<DRT::Element::line3, 3>;

// 2D elements
template class DRT::ELEMENTS::LubricationEleCalc<DRT::Element::tri3, 2>;
template class DRT::ELEMENTS::LubricationEleCalc<DRT::Element::tri3, 3>;
template class DRT::ELEMENTS::LubricationEleCalc<DRT::Element::tri6, 2>;
template class DRT::ELEMENTS::LubricationEleCalc<DRT::Element::tri6, 3>;
template class DRT::ELEMENTS::LubricationEleCalc<DRT::Element::quad4, 2>;
template class DRT::ELEMENTS::LubricationEleCalc<DRT::Element::quad4, 3>;
template class DRT::ELEMENTS::LubricationEleCalc<DRT::Element::quad8, 2>;
template class DRT::ELEMENTS::LubricationEleCalc<DRT::Element::quad8, 3>;
template class DRT::ELEMENTS::LubricationEleCalc<DRT::Element::quad9, 2>;
template class DRT::ELEMENTS::LubricationEleCalc<DRT::Element::quad9, 3>;
