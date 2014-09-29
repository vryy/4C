/*--------------------------------------------------------------------------*/
/*!
 \file acou_ele_calc.cpp
 \brief

 <pre>
 Maintainer: Svenja Schoeder
 schoeder@lnm.mw.tum.de
 http://www.lnm.mw.tum.de
 089 - 289-15271
 </pre>
 */
/*--------------------------------------------------------------------------*/

#include "acou_ele_calc.H"
#include "acou_ele_action.H"
#include "acou_utils.H"

#include "../drt_fem_general/drt_utils_boundary_integration.H"

#include "../drt_geometry/position_array.H"
#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_elementtype.H"

#include "../drt_mat/acoustic.H"

#include <Epetra_SerialDenseSolver.h>
#include <Teuchos_TimeMonitor.hpp>

namespace
{
  void zeroMatrix(Epetra_SerialDenseMatrix &mat)
  {
    std::memset(mat.A(), 0, sizeof(double) * mat.M() * mat.N());
  }

  void reshapeMatrixIfNecessary(Epetra_SerialDenseMatrix &matrix, const int nrows,const int ncols)
  {
    if (nrows != matrix.M() || ncols != matrix.N())
      matrix.Shape(nrows, ncols);
  }
}

/*----------------------------------------------------------------------*
 * Constructor
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::AcouEleCalc<distype>::AcouEleCalc()
{}

/*----------------------------------------------------------------------*
 * Action type: Evaluate
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcouEleCalc<distype>::Evaluate(DRT::ELEMENTS::Acou* ele,
    DRT::Discretization & discretization, const std::vector<int> & lm,
    Teuchos::ParameterList& params, Teuchos::RCP<MAT::Material> & mat,
    Epetra_SerialDenseMatrix& elemat1_epetra,
    Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra,
    Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra,
    const DRT::UTILS::GaussIntegration &, bool offdiag)
{
  return this->Evaluate(ele, discretization, lm, params, mat, elemat1_epetra,
      elemat2_epetra, elevec1_epetra, elevec2_epetra, elevec3_epetra, offdiag);
}

/*----------------------------------------------------------------------*
 * Evaluate
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcouEleCalc<distype>::Evaluate(DRT::ELEMENTS::Acou* ele,
    DRT::Discretization & discretization, const std::vector<int> & lm,
    Teuchos::ParameterList& params, Teuchos::RCP<MAT::Material> & mat,
    Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix&,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector&, bool offdiag)
{
  // check if this is an hdg element and init completepoly
  if (const DRT::ELEMENTS::Acou * hdgele = dynamic_cast<const DRT::ELEMENTS::Acou*>(ele))
    usescompletepoly_ = hdgele->UsesCompletePolynomialSpace();
  else
    dserror("cannot cast element to acou element");
  const ACOU::Action action = DRT::INPUT::get<ACOU::Action>(params, "action");
  InitializeShapes(ele);

  bool updateonly = false;
  shapes_->Evaluate(*ele);

  switch (action)
  {
  case ACOU::project_field:
  {
    if (mat->MaterialType() != INPAR::MAT::m_acousticmat)
      dserror("for physical type 'lossless' please supply MAT_Acoustic");
    ElementInit(ele, params);
    localSolver_->ProjectField(ele, params, mat, discretization, lm, elevec1,elevec2, dyna_);
    break;
  }
  case ACOU::project_optical_field:
  {
    if (mat->MaterialType() != INPAR::MAT::m_acousticmat)
      dserror("for physical type 'lossless' please supply MAT_Acoustic");
    localSolver_->ProjectOpticalField(ele, params, mat, discretization, lm,elevec1, elevec2, dyna_);
    break;
  }
  case ACOU::ele_init:
  {
    ElementInit(ele, params);
    break;
  }
  case ACOU::fill_restart_vecs:
  {
    bool padapty = params.get<bool>("padaptivity");
    ReadGlobalVectors(ele, discretization, lm, padapty);
    FillRestartVectors(ele, discretization);
    break;
  }
  case ACOU::ele_init_from_restart:
  {
    ElementInit(ele, params);
    ElementInitFromRestart(ele, discretization);
    break;
  }
  case ACOU::interpolate_hdg_to_node:
  {
    const bool padapty = params.get<bool>("padaptivity");
    ReadGlobalVectors(ele, discretization, lm, padapty);
    NodeBasedValues(ele, discretization, lm, elevec1, padapty);
    break;
  }
  case ACOU::interpolate_psi_to_node:
  {
    const bool padapty = params.get<bool>("padaptivity");
    double dt = params.get<double>("dt");
    ReadGlobalVectors(ele, discretization, lm, padapty);
    NodeBasedPsi(mat, ele, discretization, lm, elevec1, dt);
    break;
  }
  case ACOU::calc_acou_error:
  {
    const bool padapty = params.get<bool>("padaptivity");
    ReadGlobalVectors(ele, discretization, lm, padapty);
    ComputeError(ele, params, mat, discretization, lm, elevec1);
    break;
  }
  case ACOU::calc_abc:
  {
    int face = params.get<int>("face");
    int sumindex = 0;
    for (int i = 0; i < face; ++i)
    {
      DRT::UTILS::PolynomialSpaceParams params(DRT::UTILS::DisTypeToFaceShapeType<distype>::shape,ele->Faces()[i]->Degree(), usescompletepoly_);
      int nfdofs = DRT::UTILS::PolynomialSpaceCache<nsd_ - 1>::Instance().Create(params)->Size();
      sumindex += nfdofs;
    }

    if (!params.isParameter("nodeindices"))
      localSolver_->ComputeAbsorbingBC(ele, params, mat, face, elemat1, elevec1,sumindex);
    else
      dserror("why would you set an absorbing LINE in THREE dimensions?");

    break;
  }
  case ACOU::calc_pressuremon:
  {
    int face = params.get<int>("face");
    int sumindex = 0;
    for (int i = 0; i < face; ++i)
    {
      DRT::UTILS::PolynomialSpaceParams params(DRT::UTILS::DisTypeToFaceShapeType<distype>::shape,ele->Faces()[i]->Degree(), usescompletepoly_);
      int nfdofs = DRT::UTILS::PolynomialSpaceCache<nsd_ - 1>::Instance().Create(params)->Size();
      sumindex += nfdofs;
    }

    if (!params.isParameter("nodeindices"))
      localSolver_->ComputeSourcePressureMonitor(ele, params, mat, face, elemat1, elevec1, sumindex);
    else
      localSolver_->ComputeSourcePressureMonitorLine3D(ele, params, mat, face, elemat1, elevec1, sumindex);

    break;
  }
  case ACOU::calc_pmon_nodevals:
  {
    int face = params.get<int>("face");
    const bool padapty = params.get<bool>("padaptivity");
    int sumindex = 0;
    for (int i = 0; i < face; ++i)
    {
      DRT::UTILS::PolynomialSpaceParams params(DRT::UTILS::DisTypeToFaceShapeType<distype>::shape,ele->Faces()[i]->Degree(), usescompletepoly_);
      int nfdofs = DRT::UTILS::PolynomialSpaceCache<nsd_ - 1>::Instance().Create(params)->Size();
      sumindex += nfdofs;
    }

    ReadGlobalVectors(ele, discretization, lm, padapty);
    ComputePMonNodeVals(ele, params, mat, face, elemat1, elevec1, sumindex);

    break;
  }
  case ACOU::calc_systemmat_and_residual:
  {
    const bool resonly = params.get<bool>("resonly");
    const bool adjoint = params.get<bool>("adjoint");
    const bool padapty = params.get<bool>("padaptivity");
    double dt = params.get<double>("dt");
    dyna_ = params.get<INPAR::ACOU::DynamicType>("dynamic type");

    ReadGlobalVectors(ele, discretization, lm, padapty);
    zeroMatrix(elevec1);
    localSolver_->ComputeMatrices(mat, *ele, dt, dyna_, adjoint);

    if (!resonly)
      localSolver_->CondenseLocalPart(elemat1, dyna_);

    VectorHandling(ele, params, dt);
    localSolver_->ComputeResidual(elevec1, interiorVelnp_, interiorPressnp_, traceVal_, dyna_);

    break;
  }
  case ACOU::update_secondary_solution:
    updateonly = true; // no break here!!!
  case ACOU::update_secondary_solution_and_calc_residual:
  {
    bool adjoint = params.get<bool>("adjoint");
    bool errormaps = params.get<bool>("errormaps");
    const bool padapty = params.get<bool>("padaptivity");
    double dt = params.get<double>("dt");
    dyna_ = params.get<INPAR::ACOU::DynamicType>("dynamic type");

    ReadGlobalVectors(ele, discretization, lm, padapty);

    zeroMatrix(elevec1);
    localSolver_->ComputeMatrices(mat, *ele, dt, dyna_, adjoint);

    // this happens for DIRK time integration and is disadvantageous, since history variables are shaped according to elevec size in UpdateInteriorVariablesAndComputeResidual
    if (unsigned(elevec1.M()) != lm.size())
      elevec1.Size(lm.size());

    UpdateInteriorVariablesAndComputeResidual(discretization, params, *ele, mat, elevec1, dt, errormaps, updateonly);
    break;
  }
  default:
    dserror("unknown action supplied");
    break;
  } // switch(action)
  return 0;
}

template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::InitializeShapes(
    const DRT::ELEMENTS::Acou* ele)
{
  // TODO: check distype

  if (shapes_ == Teuchos::null)
    shapes_ = Teuchos::rcp(
        new DRT::UTILS::ShapeValues<distype>(ele->Degree(), usescompletepoly_,2 * ele->Degree()));
  else if (shapes_->degree_ != unsigned(ele->Degree()) || shapes_->usescompletepoly_ != usescompletepoly_)
    shapes_ = Teuchos::rcp(new DRT::UTILS::ShapeValues<distype>(ele->Degree(), usescompletepoly_,2 * ele->Degree()));

  if (localSolver_ == Teuchos::null)
    localSolver_ = Teuchos::rcp(new LocalSolver(*shapes_));
  else if (localSolver_->ndofs_ != shapes_->ndofs_)
    localSolver_ = Teuchos::rcp(new LocalSolver(*shapes_));

  localSolver_->FaceSpecificConstruction(ele, usescompletepoly_);
}

/*----------------------------------------------------------------------*
 * DIRKVectorHandling
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::VectorHandling(
    DRT::ELEMENTS::Acou * ele, Teuchos::ParameterList& params, double& dt)
{
  if (dyna_ == INPAR::ACOU::acou_dirk23 || dyna_ == INPAR::ACOU::acou_dirk33 ||
      dyna_ == INPAR::ACOU::acou_dirk34 || dyna_ == INPAR::ACOU::acou_dirk54)
  {
    // when we are here, the time integrator is dirk, and we have to do the update things
    double dirk_a[6][6];
    double dirk_b[6];
    int dirk_q;
    ACOU::FillDIRKValues(dyna_, dirk_a, dirk_b, dirk_q);
    int stage = params.get<int>("stage");

    if (stage == 0)
    {
      ele->eleinteriorPressnp_ = ele->eleinteriorPressn_;
      ele->eleinteriorVelnp_ = ele->eleinteriorVeln_;
    }

    Epetra_SerialDenseVector tempVecp1(shapes_->ndofs_);
    Epetra_SerialDenseVector tempVecp2(shapes_->ndofs_);
    Epetra_SerialDenseVector tempVecv1(shapes_->ndofs_ * nsd_);
    Epetra_SerialDenseVector tempVecv2(shapes_->ndofs_ * nsd_);

    ele->elesp_[stage] = ele->eleinteriorPressn_;
    ele->elesp_[stage].Scale(1.0 / dt); // dt includes a_ii
    ele->elesv_[stage] = ele->eleinteriorVeln_;
    ele->elesv_[stage].Scale(1.0 / dt); // dt includes a_ii

    for (int i = 0; i < stage; ++i)
    {
      tempVecp1 = ele->elesp_[i];
      tempVecp1.Scale(-1.0);
      tempVecp2 = ele->eleyp_[i];
      tempVecp2.Scale(1.0 / dt);
      tempVecp1 += tempVecp2;
      tempVecp1.Scale(dirk_a[stage][i] / dirk_a[stage][stage]);
      ele->elesp_[stage] += tempVecp1;
      tempVecv1 = ele->elesv_[i];
      tempVecv1.Scale(-1.0);
      tempVecv2 = ele->eleyv_[i];
      tempVecv2.Scale(1.0 / dt);
      tempVecv1 += tempVecv2;
      tempVecv1.Scale(dirk_a[stage][i] / dirk_a[stage][stage]);
      ele->elesv_[stage] += tempVecv1;
    }

    // these vectors s are used for the calculation of the residual -> write them to used local solver variable
    interiorPressnp_ = ele->elesp_[stage];
    interiorVelnp_ = ele->elesv_[stage];
    interiorPressnp_.Scale(dt);
    interiorVelnp_.Scale(dt);
  }
  else if (dyna_ == INPAR::ACOU::acou_bdf2)
  {
    for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_; ++i)
      interiorVelnp_[i] = 4.0 / 3.0 * interiorVelnp_[i] - 1.0 / 3.0 * interiorVelnm_[i];
    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      interiorPressnp_[i] = 4.0 / 3.0 * interiorPressnp_[i] - 1.0 / 3.0 * interiorPressnm_[i];
  }
  else if (dyna_ == INPAR::ACOU::acou_bdf3)
  {
    for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_; ++i)
      interiorVelnp_[i] = 18.0 / 11.0 * interiorVelnp_[i] - 9.0 / 11.0 * interiorVelnm_[i] + 2.0 / 11.0 * interiorVelnmm_[i];
    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      interiorPressnp_[i] = 18.0 / 11.0 * interiorPressnp_[i] - 9.0 / 11.0 * interiorPressnm_[i] + 2.0 / 11.0 * interiorPressnmm_[i];
  }
  else if (dyna_ == INPAR::ACOU::acou_bdf4)
  {
    for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_; ++i)
      interiorVelnp_[i] = 48.0 / 25.0 * interiorVelnp_[i] - 36.0 / 25.0 * interiorVelnm_[i] + 16.0 / 25.0 * interiorVelnmm_[i] - 3.0 / 25.0 * interiorVelnmmm_[i];
    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      interiorPressnp_[i] = 48.0 / 25.0 * interiorPressnp_[i] - 36.0 / 25.0 * interiorPressnm_[i] + 16.0 / 25.0 * interiorPressnmm_[i] - 3.0 / 25.0 * interiorPressnmmm_[i];
  }
  return;
}

/*----------------------------------------------------------------------*
 * ReadGlobalVectors
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::ReadGlobalVectors(DRT::Element * ele,
    DRT::Discretization & discretization, const std::vector<int> & lm,
    const bool padaptivity)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::ReadGlobalVectors");
  DRT::ELEMENTS::Acou * acouele = dynamic_cast<DRT::ELEMENTS::Acou*>(ele);

  // read vectors from element storage
  reshapeMatrixIfNecessary(interiorVelnp_, acouele->eleinteriorVelnp_.M(), 1);
  reshapeMatrixIfNecessary(interiorPressnp_, acouele->eleinteriorPressnp_.M(),1);
  reshapeMatrixIfNecessary(interiorVeln_, acouele->eleinteriorVeln_.M(), 1);
  reshapeMatrixIfNecessary(interiorPressn_, acouele->eleinteriorPressn_.M(), 1);
  reshapeMatrixIfNecessary(interiorVelnm_, acouele->eleinteriorVelnm_.M(), 1);
  reshapeMatrixIfNecessary(interiorPressnm_, acouele->eleinteriorPressnm_.M(),1);
  reshapeMatrixIfNecessary(interiorVelnmm_, acouele->eleinteriorVelnmm_.M(), 1);
  reshapeMatrixIfNecessary(interiorPressnmm_, acouele->eleinteriorPressnmm_.M(),1);
  reshapeMatrixIfNecessary(interiorVelnmmm_, acouele->eleinteriorVelnmmm_.M(),1);
  reshapeMatrixIfNecessary(interiorPressnmmm_,acouele->eleinteriorPressnmmm_.M(), 1);

  interiorVelnp_ = acouele->eleinteriorVelnp_;
  interiorPressnp_ = acouele->eleinteriorPressnp_;
  interiorVeln_ = acouele->eleinteriorVeln_;
  interiorPressn_ = acouele->eleinteriorPressn_;
  interiorVelnm_ = acouele->eleinteriorVelnm_;
  interiorPressnm_ = acouele->eleinteriorPressnm_;
  interiorVelnmm_ = acouele->eleinteriorVelnmm_;
  interiorPressnmm_ = acouele->eleinteriorPressnmm_;
  interiorVelnmmm_ = acouele->eleinteriorVelnmmm_;
  interiorPressnmmm_ = acouele->eleinteriorPressnmmm_;

  // read vectors from time integrator
  if (discretization.HasState("trace")) // in case of "update interior variables"
  {
    traceVal_.resize(lm.size());
    Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState("trace");
    DRT::UTILS::ExtractMyValues(*matrix_state, traceVal_, lm);
  }
  if (discretization.HasState("trace_m")) // in case of trapezoidal time integration
  {
    traceValm_.resize(lm.size());
    Teuchos::RCP<const Epetra_Vector> matrix_state_m = discretization.GetState("trace_m");
    DRT::UTILS::ExtractMyValues(*matrix_state_m, traceValm_, lm);
  }

  return;
} // ReadGlobalVectors

/*----------------------------------------------------------------------*
 * FillRestartVectors
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::FillRestartVectors(DRT::Element * ele,
    DRT::Discretization & discretization)
{
  // sort this back to the interior values vector
  std::vector<double> interiorValnp(shapes_->ndofs_ * (nsd_ + 1));
  for (unsigned int i = 0; i < interiorValnp.size(); ++i)
  {
    if ((i + 1) % (nsd_ + 1) == 0)
      interiorValnp[i] = interiorPressnp_((i + 1) / (nsd_ + 1) - 1);
    else
    {
      int xyz = i % (nsd_ + 1); // 0 for x, 1 for y and 2 for z (for 3D)
      interiorValnp[i] = interiorVelnp_(xyz * shapes_->ndofs_ + i / (nsd_ + 1));
    }
  }

  // tell this change in the interior variables the discretization
  std::vector<int> localDofs = discretization.Dof(1, ele);
  const Epetra_Map* intdofcolmap = discretization.DofColMap(1);
  {
    Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1,"intvelnp");
    Epetra_Vector& secondary = const_cast<Epetra_Vector&>(*matrix_state);
    for (unsigned int i = 0; i < localDofs.size(); ++i)
    {
      const int lid = intdofcolmap->LID(localDofs[i]);
      secondary[lid] = interiorValnp[i];
    }
  }

  {
    Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1,"intveln");
    Epetra_Vector& secondary = const_cast<Epetra_Vector&>(*matrix_state);
    for (unsigned int i = 0; i < interiorValnp.size(); ++i)
    {
      if ((i + 1) % (nsd_ + 1) == 0)
        interiorValnp[i] = interiorPressn_((i + 1) / (nsd_ + 1) - 1);
      else
      {
        int xyz = i % (nsd_ + 1); // 0 for x, 1 for y and 2 for z (for 3D)
        interiorValnp[i] = interiorVeln_(xyz * shapes_->ndofs_ + i / (nsd_ + 1));
      }
    }
    for (unsigned int i = 0; i < localDofs.size(); ++i)
    {
      const int lid = intdofcolmap->LID(localDofs[i]);
      secondary[lid] = interiorValnp[i];
    }
  }

  // only do this if desired!
  if (interiorPressnm_.M() > 0 && discretization.HasState(1, "intvelnm"))
  {
    Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1,"intvelnm");
    Epetra_Vector& secondary = const_cast<Epetra_Vector&>(*matrix_state);
    for (unsigned int i = 0; i < interiorValnp.size(); ++i)
    {
      if ((i + 1) % (nsd_ + 1) == 0)
        interiorValnp[i] = interiorPressnm_((i + 1) / (nsd_ + 1) - 1);
      else
      {
        int xyz = i % (nsd_ + 1); // 0 for x, 1 for y and 2 for z (for 3D)
        interiorValnp[i] = interiorVelnm_(xyz * shapes_->ndofs_ + i / (nsd_ + 1));
      }
    }
    for (unsigned int i = 0; i < localDofs.size(); ++i)
    {
      const int lid = intdofcolmap->LID(localDofs[i]);
      secondary[lid] = interiorValnp[i];
    }
  }

  // only do this if desired!
  if (interiorPressnmm_.M() > 0 && discretization.HasState(1, "intvelnmm"))
  {
    Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1,"intvelnmm");
    Epetra_Vector& secondary = const_cast<Epetra_Vector&>(*matrix_state);
    for (unsigned int i = 0; i < interiorValnp.size(); ++i)
    {
      if ((i + 1) % (nsd_ + 1) == 0)
        interiorValnp[i] = interiorPressnmm_((i + 1) / (nsd_ + 1) - 1);
      else
      {
        int xyz = i % (nsd_ + 1); // 0 for x, 1 for y and 2 for z (for 3D)
        interiorValnp[i] = interiorVelnmm_(xyz * shapes_->ndofs_ + i / (nsd_ + 1));
      }
    }
    for (unsigned int i = 0; i < localDofs.size(); ++i)
    {
      const int lid = intdofcolmap->LID(localDofs[i]);
      secondary[lid] = interiorValnp[i];
    }
  }

  // only do this if desired!
  if (interiorPressnmmm_.M() > 0 && discretization.HasState(1, "intvelnmmm"))
  {
    Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1,"intvelnmmm");
    Epetra_Vector& secondary = const_cast<Epetra_Vector&>(*matrix_state);
    for (unsigned int i = 0; i < interiorValnp.size(); ++i)
    {
      if ((i + 1) % (nsd_ + 1) == 0)
        interiorValnp[i] = interiorPressnmmm_((i + 1) / (nsd_ + 1) - 1);
      else
      {
        int xyz = i % (nsd_ + 1); // 0 for x, 1 for y and 2 for z (for 3D)
        interiorValnp[i] = interiorVelnmmm_(xyz * shapes_->ndofs_ + i / (nsd_ + 1));
      }
    }
    for (unsigned int i = 0; i < localDofs.size(); ++i)
    {
      const int lid = intdofcolmap->LID(localDofs[i]);
      secondary[lid] = interiorValnp[i];
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 * ElementInitFromRestart
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::ElementInitFromRestart(
    DRT::Element * ele, DRT::Discretization & discretization)
{
  DRT::ELEMENTS::Acou * acouele = dynamic_cast<DRT::ELEMENTS::Acou*>(ele);
  std::vector<double> interiorValnp(shapes_->ndofs_ * (nsd_ + 1));

  Teuchos::RCP<const Epetra_Vector> intvel = discretization.GetState(1,"intvelnp");
  std::vector<int> localDofs1 = discretization.Dof(1, ele);
  DRT::UTILS::ExtractMyValues(*intvel, interiorValnp, localDofs1);

  // now write this in corresponding interiorVelnp_ and interiorPressnp_
  for (unsigned int i = 0; i < interiorValnp.size(); ++i)
  {
    if ((i + 1) % (nsd_ + 1) == 0)
      acouele->eleinteriorPressnp_((i + 1) / (nsd_ + 1) - 1) = interiorValnp[i];
    else
    {
      int xyz = i % (nsd_ + 1);
      acouele->eleinteriorVelnp_(xyz * shapes_->ndofs_ + i / (nsd_ + 1)) = interiorValnp[i];
    }
  }

  intvel = discretization.GetState(1, "intveln");
  DRT::UTILS::ExtractMyValues(*intvel, interiorValnp, localDofs1);
  for (unsigned int i = 0; i < interiorValnp.size(); ++i)
  {
    if ((i + 1) % (nsd_ + 1) == 0)
      acouele->eleinteriorPressn_((i + 1) / (nsd_ + 1) - 1) = interiorValnp[i];
    else
    {
      int xyz = i % (nsd_ + 1);
      acouele->eleinteriorVeln_(xyz * shapes_->ndofs_ + i / (nsd_ + 1)) = interiorValnp[i];
    }
  }
  if (discretization.HasState(1, "intvelnm") && acouele->eleinteriorPressnm_.M() > 0)
  {
    intvel = discretization.GetState(1, "intvelnm");
    DRT::UTILS::ExtractMyValues(*intvel, interiorValnp, localDofs1);
    for (unsigned int i = 0; i < interiorValnp.size(); ++i)
    {
      if ((i + 1) % (nsd_ + 1) == 0)
        acouele->eleinteriorPressnm_((i + 1) / (nsd_ + 1) - 1) = interiorValnp[i];
      else
      {
        int xyz = i % (nsd_ + 1);
        acouele->eleinteriorVelnm_(xyz * shapes_->ndofs_ + i / (nsd_ + 1)) = interiorValnp[i];
      }
    }
  }

  if (discretization.HasState(1, "intvelnmm") && acouele->eleinteriorPressnmm_.M() > 0)
  {
    intvel = discretization.GetState(1, "intvelnmm");
    DRT::UTILS::ExtractMyValues(*intvel, interiorValnp, localDofs1);
    for (unsigned int i = 0; i < interiorValnp.size(); ++i)
    {
      if ((i + 1) % (nsd_ + 1) == 0)
        acouele->eleinteriorPressnmm_((i + 1) / (nsd_ + 1) - 1) = interiorValnp[i];
      else
      {
        int xyz = i % (nsd_ + 1);
        acouele->eleinteriorVelnmm_(xyz * shapes_->ndofs_ + i / (nsd_ + 1)) = interiorValnp[i];
      }
    }
  }

  if (discretization.HasState(1, "intvelnmmm") && acouele->eleinteriorPressnmmm_.M() > 0)
  {
    intvel = discretization.GetState(1, "intvelnmmm");
    DRT::UTILS::ExtractMyValues(*intvel, interiorValnp, localDofs1);
    for (unsigned int i = 0; i < interiorValnp.size(); ++i)
    {
      if ((i + 1) % (nsd_ + 1) == 0)
        acouele->eleinteriorPressnmmm_((i + 1) / (nsd_ + 1) - 1) = interiorValnp[i];
      else
      {
        int xyz = i % (nsd_ + 1);
        acouele->eleinteriorVelnmmm_(xyz * shapes_->ndofs_ + i / (nsd_ + 1)) = interiorValnp[i];
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 * ComputeError
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::ComputeError(DRT::ELEMENTS::Acou* ele,
    Teuchos::ParameterList& params, Teuchos::RCP<MAT::Material>& mat,
    DRT::Discretization& discretization, const std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec)
{
  shapes_->Evaluate(*ele);

  double time = params.get<double>("time");

  if (params.get<INPAR::ACOU::CalcError>("error calculation") != INPAR::ACOU::calcerror_1d)
    dserror("no analytical solution available");

  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double c = actmat->SpeedofSound();

  // get function
  const int *start_func = params.getPtr<int>("funct");

  double err_p = 0.0,norm_p = 0.0;
  for(unsigned  int q=0; q<shapes_->nqpoints_; ++q)
  {
    double numerical = 0.0;
    for (unsigned int i=0; i<shapes_->ndofs_; ++i)
    numerical += shapes_->shfunct(i,q) * interiorPressnp_(i);

    double exact = 0.0;
    double xyz[nsd_];
    for (unsigned int d=0; d<nsd_; ++d)
    xyz[d] = shapes_->xyzreal(d,q);

    // careful: we assume that the x-axis is the orientation of your 1D problem
    double xyzp[nsd_];
    double xyzm[nsd_];
    xyzp[0]=xyz[0]+c*time;
    xyzm[0]=xyz[0]-c*time;
    for (unsigned int d=1; d<nsd_; ++d)
    {
      xyzp[d] = xyz[d];
      xyzm[d] = xyz[d];
    }

    // calculation of the analytical solution
    exact = 0.5 * DRT::Problem::Instance()->Funct(*start_func-1).Evaluate(0,xyzp,0.0,NULL)
          + 0.5 * DRT::Problem::Instance()->Funct(*start_func-1).Evaluate(0,xyzm,0.0,NULL);

    err_p += ( exact - numerical ) * ( exact - numerical ) * shapes_->jfac(q);
    norm_p += exact * exact * shapes_->jfac(q);
  }

  elevec[1] += err_p;
  elevec[3] += norm_p;

  return;
}

/*----------------------------------------------------------------------*
 * Element init
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::ElementInit(DRT::ELEMENTS::Acou* ele,
    Teuchos::ParameterList& params)
{
  shapes_->Evaluate(*ele);

  // each element has to store the interior vectors by itseld, p-adaptivity or not
  // so, shape it, as you need it
  dyna_ = params.get<INPAR::ACOU::DynamicType>("dynamic type");
  ele->eleinteriorVelnp_.Shape(shapes_->ndofs_ * nsd_, 1);
  ele->eleinteriorPressnp_.Shape(shapes_->ndofs_, 1);
  ele->eleinteriorVeln_.Shape(shapes_->ndofs_ * nsd_, 1);
  ele->eleinteriorPressn_.Shape(shapes_->ndofs_, 1);
  switch (dyna_)
  {
  case INPAR::ACOU::acou_bdf4:
    ele->eleinteriorVelnmmm_.Shape(shapes_->ndofs_ * nsd_, 1);
    ele->eleinteriorPressnmmm_.Shape(shapes_->ndofs_, 1); // no break here!
  case INPAR::ACOU::acou_bdf3:
    ele->eleinteriorVelnmm_.Shape(shapes_->ndofs_ * nsd_, 1); // no break here!
    ele->eleinteriorPressnmm_.Shape(shapes_->ndofs_, 1);
  case INPAR::ACOU::acou_bdf2:
    ele->eleinteriorVelnm_.Shape(shapes_->ndofs_ * nsd_, 1);
    ele->eleinteriorPressnm_.Shape(shapes_->ndofs_, 1);
    break; // here you go!
  case INPAR::ACOU::acou_dirk23:
    ele->elesp_.resize(2);
    ele->eleyp_.resize(2);
    ele->elefp_.resize(2);
    ele->elesv_.resize(2);
    ele->eleyv_.resize(2);
    ele->elefv_.resize(2);
    break;
  case INPAR::ACOU::acou_dirk33:
  case INPAR::ACOU::acou_dirk34:
    ele->elesv_.resize(3);
    ele->eleyv_.resize(3);
    ele->elefv_.resize(3);
    ele->elesp_.resize(3);
    ele->eleyp_.resize(3);
    ele->elefp_.resize(3);
    break;
  case INPAR::ACOU::acou_dirk54:
    ele->elesv_.resize(5);
    ele->eleyv_.resize(5);
    ele->elefv_.resize(5);
    ele->elesp_.resize(5);
    ele->eleyp_.resize(5);
    ele->elefp_.resize(5);
    break;
  default:
    break; // do nothing for impl and trap
  }
  for (unsigned int i = 0; i < ele->elesp_.size(); ++i) // zero size if not dirk
  {
    ele->elesp_[i].Shape(shapes_->ndofs_, 1);
    ele->eleyp_[i].Shape(shapes_->ndofs_, 1);
    ele->elefp_[i].Shape(shapes_->ndofs_, 1);
    ele->elesv_[i].Shape(shapes_->ndofs_ * nsd_, 1);
    ele->eleyv_[i].Shape(shapes_->ndofs_ * nsd_, 1);
    ele->elefv_[i].Shape(shapes_->ndofs_ * nsd_, 1);
  }

  if (params.get<bool>("padaptivity"))
    ele->elenodeTrace_.Shape(nen_, 1);

  return;
}

/*----------------------------------------------------------------------*
 * ProjectField
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::ProjectField(
    DRT::ELEMENTS::Acou* ele, Teuchos::ParameterList& params,
    Teuchos::RCP<MAT::Material>& mat, DRT::Discretization& discretization,
    const std::vector<int>& lm, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, INPAR::ACOU::DynamicType dyna)
{
  shapes_.Evaluate(*ele);

  // reshape elevec2 as matrix
  dsassert(elevec2.M() == 0 ||
      unsigned(elevec2.M()) == (nsd_+1)*shapes_.ndofs_, "Wrong size in project vector 2");

  // get function
  const int *start_func = params.getPtr<int>("funct");

  // internal variables
  if (elevec2.M() > 0)
  {
    Epetra_SerialDenseMatrix localMat(shapes_.ndofs_, 1);

    for (unsigned int q = 0; q < shapes_.nqpoints_; ++q)
    {
      const double fac = shapes_.jfac(q);
      const double sqrtfac = std::sqrt(fac);
      double xyz[nsd_];
      for (unsigned int d = 0; d < nsd_; ++d)
        xyz[d] = shapes_.xyzreal(d, q); // coordinates of quadrature point in real coordinates
      double p;
      dsassert(start_func != NULL,"funct not set for initial value");
      EvaluateAll(*start_func, xyz, p, 1.0); // u and p at quadrature point

      // now fill the components in the one-sided mass matrix and the right hand side
      for (unsigned int i = 0; i < shapes_.ndofs_; ++i)
      {
        massPart(i, q) = shapes_.shfunct(i, q) * sqrtfac;
        localMat(i, 0) += shapes_.shfunct(i, q) * p * fac;
      }
    }

    Mmat.Multiply('N', 'T', 1., massPart, massPart, 0.);
    {
      Epetra_SerialDenseSolver inverseMass;
      inverseMass.SetMatrix(Mmat);
      inverseMass.SetVectors(localMat, localMat);
      inverseMass.Solve();
    }

    for (unsigned int r = 0; r < shapes_.ndofs_; ++r)
      ele->eleinteriorPressnp_(r) += localMat(r, 0); // pressure

  }
  switch (dyna)
  {
  case INPAR::ACOU::acou_bdf4:
    ele->eleinteriorPressnmmm_ = ele->eleinteriorPressnp_; // no break here
  case INPAR::ACOU::acou_bdf3:
    ele->eleinteriorPressnmm_ = ele->eleinteriorPressnp_; // no break here
  case INPAR::ACOU::acou_bdf2:
    ele->eleinteriorPressnm_ = ele->eleinteriorPressnp_; // no break here                                          // here you go
  default:
    ele->eleinteriorPressn_ = ele->eleinteriorPressnp_;
    break;
  }

  // in case this paramter "faceconsider" is set, we are applying Dirichlet values
  // and have to evaluate the trace field for the given face!
  if (params.isParameter("faceconsider"))
  {
    Teuchos::Array<int> *functno = params.getPtr<Teuchos::Array<int> >("funct");
    const unsigned int *faceConsider = params.getPtr<unsigned int>("faceconsider");
    //Teuchos::Array<int> *onoff = params.getPtr<Teuchos::Array<int> >("onoff");
    //double *time = params.getPtr<double>("time");

    DRT::UTILS::ShapeValuesFaceParams svfparams(
        ele->Faces()[*faceConsider]->Degree(),
        shapes_.usescompletepoly_,
        2 * ele->Faces()[*faceConsider]->Degree());
    shapesface_ = DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
    shapesface_->EvaluateFace(*ele, *faceConsider, shapes_);

    Epetra_SerialDenseMatrix mass(shapesface_->nfdofs_, shapesface_->nfdofs_);
    Epetra_SerialDenseVector trVec(shapesface_->nfdofs_);

    for (unsigned int q = 0; q < shapesface_->nqpoints_; ++q)
    {
      const double fac = shapesface_->jfac(q);
      double xyz[nsd_];
      for (unsigned int d = 0; d < nsd_; ++d)
        xyz[d] = shapesface_->xyzreal(d, q);
      double p;

      EvaluateAll((*functno)[0], xyz, p, 1.0);

      // now fill the components in the mass matrix and the right hand side
      for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
      {
        // mass matrix
        for (unsigned int j = 0; j < shapesface_->nfdofs_; ++j)
          mass(i, j) += shapesface_->shfunct(i, q) * shapesface_->shfunct(j, q) * fac;
        trVec(i) += shapesface_->shfunct(i, q) * p * fac;
      }
    }

    Epetra_SerialDenseSolver inverseMass;
    inverseMass.SetMatrix(mass);
    inverseMass.SetVectors(trVec, trVec);
    inverseMass.Solve();

    for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
      elevec1(i) = trVec(i);
  }

  return 0;
}

/*----------------------------------------------------------------------*
 * ProjectOpticalField
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::ProjectOpticalField(
    DRT::ELEMENTS::Acou* ele, Teuchos::ParameterList& params,
    Teuchos::RCP<MAT::Material>& mat, DRT::Discretization& discretization,
    const std::vector<int>& lm, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, INPAR::ACOU::DynamicType dyna)
{

  shapes_.Evaluate(*ele);

  bool meshconform = params.get<bool>("mesh conform");

  // this is the easy case: the corresponding optical element has exactly the same nodes
  // and hence we can easily get the scatra dofs of the nodes and substitute the pressure
  // values
  int numlightnode = nen_;

  double lightxyz[numlightnode][nsd_];
  double values[numlightnode];
  double absorptioncoeff = -1.0;

  if (meshconform)
  {
    // get the absorption coefficient
    absorptioncoeff = params.get<double>("absorption");

    Teuchos::RCP<std::vector<double> > nodevals = params.get<Teuchos::RCP<std::vector<double> > >("nodevals");

    if (int((*nodevals).size() / (nsd_ + 1)) != numlightnode)
      dserror("node number not matching");

    // meanwhile: check if we are at least near the element, otherwise we can already go home
    for (int i = 0; i < numlightnode; ++i)
    {
      for (unsigned int j = 0; j < nsd_; ++j)
        lightxyz[i][j] = (*nodevals)[i * (nsd_ + 1) + j];
      values[i] = (*nodevals)[i * (nsd_ + 1) + nsd_];
    }

  }
  else
  {
    // get node based values!
    Teuchos::RCP<const Epetra_Vector> matrix_state = params.get<Teuchos::RCP<Epetra_Vector> >("pressurenode");
    std::vector<double> nodevals;
    DRT::UTILS::ExtractMyNodeBasedValues(ele, nodevals, *matrix_state);
    if (nodevals.size() != nen_)
      dserror("node number not matching: %d vs. %d", nodevals.size(), nen_);

    for (int i = 0; i < numlightnode; ++i)
    {
      for (unsigned int j = 0; j < nsd_; ++j)
        lightxyz[i][j] = shapes_.xyze(j, i);
      values[i] = nodevals[i];
    }
  }

  // now we have locations "lightxyz" and corresponding phis in "values"
  // what we want to do now is to evaluate the field in the locations of the
  // Gauss points in this element, to assign an initial distribution
  // so we do exactly the same as in ProjectInitalField BUT call an other
  // evaluate function which will interpolate from given lightxyz and values

  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double rho = actmat->Density();
  double pulse = params.get<double>("pulse");

  // internal variables
  if (elevec2.M() > 0)
  {
    Epetra_SerialDenseMatrix localMat(shapes_.ndofs_, 1);

    for (unsigned int q = 0; q < shapes_.nqpoints_; ++q)
    {
      const double fac = shapes_.jfac(q);
      const double sqrtfac = std::sqrt(fac);
      double xyz[nsd_];
      for (unsigned int d = 0; d < nsd_; ++d)
        xyz[d] = shapes_.xyzreal(d, q); // coordinates of quadrature point in real coordinates
      double p = 0.0;
      EvaluateLight(lightxyz, values, numlightnode, xyz, p, rho / pulse, absorptioncoeff); // p at quadrature point

      for (unsigned int i = 0; i < shapes_.ndofs_; ++i)
      {
        massPart(i, q) = shapes_.shfunct(i, q) * sqrtfac;
        localMat(i, 0) += shapes_.shfunct(i, q) * p * fac;
      }
    }

    Mmat.Multiply('N', 'T', 1., massPart, massPart, 0.);
    {
      Epetra_SerialDenseSolver inverseMass;
      inverseMass.SetMatrix(Mmat);
      inverseMass.SetVectors(localMat, localMat);
      inverseMass.Solve();
    }

    for (unsigned int r = 0; r < shapes_.ndofs_; ++r)
    {
      elevec2[r * (nsd_ + 1) + nsd_] += localMat(r, 0); // pressure
      ele->eleinteriorPressnp_(r) += localMat(r, 0); // pressure
    }
  }
  switch (dyna)
  {
  case INPAR::ACOU::acou_bdf4:
    ele->eleinteriorPressnmmm_ = ele->eleinteriorPressnp_; // no break here
  case INPAR::ACOU::acou_bdf3:
    ele->eleinteriorPressnmm_ = ele->eleinteriorPressnp_; // no break here
  case INPAR::ACOU::acou_bdf2:
    ele->eleinteriorPressnm_ = ele->eleinteriorPressnp_; // no break here                                          // here you go
  default:
    ele->eleinteriorPressn_ = ele->eleinteriorPressnp_;
    break;
  }

  return 0;
}

/*----------------------------------------------------------------------*
 * EvaluateAll
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::EvaluateAll(
    const int start_func, const double(&xyz)[nsd_], double &p, double rho) const
{
  p = DRT::Problem::Instance()->Funct(start_func - 1).Evaluate(0, xyz, 0.0,
      NULL);
  return;
}

/*----------------------------------------------------------------------*
 * EvaluateLight
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::EvaluateLight(
    double lightxyz[][nsd_], double values[], int numnode,
    const double(&xyz)[nsd_], double &p, double rho,
    double absorptioncoeff) const
{

  // interpolation from nodes
  if (distype == DRT::Element::quad4)
  {
    if (numnode != 4)
      dserror("wrong number of nodes given");

      //*******************************************************************
      LINALG::Matrix<4, 4> coeff(true);
      for(int i=0; i<4; ++i)
      {
        coeff(i,0)=1.0;
        coeff(i,1)=lightxyz[i][0];
        coeff(i,2)=lightxyz[i][1];
        coeff(i,3)=lightxyz[i][0]*lightxyz[i][1];
      }
      LINALG::Matrix<4,4> coeff_N(true);
      coeff_N(0,0)=1.0;
      coeff_N(1,1)=1.0;
      coeff_N(2,2)=1.0;
      coeff_N(3,3)=1.0;
      {
        LINALG::FixedSizeSerialDenseSolver<4,4,4> inverseCoeff;
        inverseCoeff.SetMatrix(coeff);
        inverseCoeff.SetVectors(coeff_N,coeff_N);
        int err = inverseCoeff.Solve();
        if(err != 0) dserror("Inversion of matrix in light evaluation failed with errorcode %d",err);
      }

      p = ( coeff_N(0,0) + coeff_N(1,0) * xyz[0] + coeff_N(2,0) * xyz[1] + coeff_N(3,0) * xyz[0] * xyz[1] ) * values[0]
      + ( coeff_N(0,1) + coeff_N(1,1) * xyz[0] + coeff_N(2,1) * xyz[1] + coeff_N(3,1) * xyz[0] * xyz[1] ) * values[1]
      + ( coeff_N(0,2) + coeff_N(1,2) * xyz[0] + coeff_N(2,2) * xyz[1] + coeff_N(3,2) * xyz[0] * xyz[1] ) * values[2]
      + ( coeff_N(0,3) + coeff_N(1,3) * xyz[0] + coeff_N(2,3) * xyz[1] + coeff_N(3,3) * xyz[0] * xyz[1] ) * values[3];
      p *= -absorptioncoeff;

    }
    else if (distype == DRT::Element::hex8)
    {
      if(numnode != 8) dserror("wrong number of nodes given");

      //*******************************************************************
      LINALG::Matrix<8,8> coeff(true);
      for(int i=0; i<8; ++i)
      {
        coeff(i,0)=1.0;
        coeff(i,1)=lightxyz[i][0];
        coeff(i,2)=lightxyz[i][1];
        coeff(i,3)=lightxyz[i][2];
        coeff(i,4)=lightxyz[i][0]*lightxyz[i][1];
        coeff(i,5)=lightxyz[i][0]*lightxyz[i][2];
        coeff(i,6)=lightxyz[i][1]*lightxyz[i][2];
        coeff(i,7)=lightxyz[i][0]*lightxyz[i][1]*lightxyz[i][2];
      }
      LINALG::Matrix<8,8> coeff_N(true);
      for(int i=0; i<8; ++i)
      coeff_N(i,i)=1.0;

      {
        LINALG::FixedSizeSerialDenseSolver<8,8,8> inverseCoeff;
        inverseCoeff.SetMatrix(coeff);
        inverseCoeff.SetVectors(coeff_N,coeff_N);
        int err = inverseCoeff.Solve();
        if(err != 0) dserror("Inversion of matrix in light evaluation failed with errorcode %d",err);
      }

      p = ( coeff_N(0,0) + coeff_N(1,0) * xyz[0] + coeff_N(2,0) * xyz[1] + coeff_N(3,0) * xyz[2] + coeff_N(4,0) * xyz[0] * xyz[1] + coeff_N(5,0) * xyz[0] * xyz[2] + coeff_N(6,0) * xyz[1] * xyz[2] + coeff_N(7,0) * xyz[0] * xyz[1] * xyz[2] ) * values[0]
      + ( coeff_N(0,1) + coeff_N(1,1) * xyz[0] + coeff_N(2,1) * xyz[1] + coeff_N(3,1) * xyz[2] + coeff_N(4,1) * xyz[0] * xyz[1] + coeff_N(5,1) * xyz[0] * xyz[2] + coeff_N(6,1) * xyz[1] * xyz[2] + coeff_N(7,1) * xyz[0] * xyz[1] * xyz[2] ) * values[1]
      + ( coeff_N(0,2) + coeff_N(1,2) * xyz[0] + coeff_N(2,2) * xyz[1] + coeff_N(3,2) * xyz[2] + coeff_N(4,2) * xyz[0] * xyz[1] + coeff_N(5,2) * xyz[0] * xyz[2] + coeff_N(6,2) * xyz[1] * xyz[2] + coeff_N(7,2) * xyz[0] * xyz[1] * xyz[2] ) * values[2]
      + ( coeff_N(0,3) + coeff_N(1,3) * xyz[0] + coeff_N(2,3) * xyz[1] + coeff_N(3,3) * xyz[2] + coeff_N(4,3) * xyz[0] * xyz[1] + coeff_N(5,3) * xyz[0] * xyz[2] + coeff_N(6,3) * xyz[1] * xyz[2] + coeff_N(7,3) * xyz[0] * xyz[1] * xyz[2] ) * values[3]
      + ( coeff_N(0,4) + coeff_N(1,4) * xyz[0] + coeff_N(2,4) * xyz[1] + coeff_N(3,4) * xyz[2] + coeff_N(4,4) * xyz[0] * xyz[1] + coeff_N(5,4) * xyz[0] * xyz[2] + coeff_N(6,4) * xyz[1] * xyz[2] + coeff_N(7,4) * xyz[0] * xyz[1] * xyz[2] ) * values[4]
      + ( coeff_N(0,5) + coeff_N(1,5) * xyz[0] + coeff_N(2,5) * xyz[1] + coeff_N(3,5) * xyz[2] + coeff_N(4,5) * xyz[0] * xyz[1] + coeff_N(5,5) * xyz[0] * xyz[2] + coeff_N(6,5) * xyz[1] * xyz[2] + coeff_N(7,5) * xyz[0] * xyz[1] * xyz[2] ) * values[5]
      + ( coeff_N(0,6) + coeff_N(1,6) * xyz[0] + coeff_N(2,6) * xyz[1] + coeff_N(3,6) * xyz[2] + coeff_N(4,6) * xyz[0] * xyz[1] + coeff_N(5,6) * xyz[0] * xyz[2] + coeff_N(6,6) * xyz[1] * xyz[2] + coeff_N(7,6) * xyz[0] * xyz[1] * xyz[2] ) * values[6]
      + ( coeff_N(0,7) + coeff_N(1,7) * xyz[0] + coeff_N(2,7) * xyz[1] + coeff_N(3,7) * xyz[2] + coeff_N(4,7) * xyz[0] * xyz[1] + coeff_N(5,7) * xyz[0] * xyz[2] + coeff_N(6,7) * xyz[1] * xyz[2] + coeff_N(7,7) * xyz[0] * xyz[1] * xyz[2] ) * values[7];
      p *= -absorptioncoeff;
    }
    else
    dserror("not yet implemented");

    return;
  }

      /*----------------------------------------------------------------------*
       * NodeBasedValues
       *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::NodeBasedValues(
    DRT::ELEMENTS::Acou* ele, DRT::Discretization& discretization,
    const std::vector<int>& lm, Epetra_SerialDenseVector& elevec1,
    bool padaptivity)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::NodeBasedValues");

  dsassert(elevec1.M() == (int)nen_*(nsd_+2)+1, "Vector does not have correct size");
  elevec1.Scale(0.0);
  Epetra_SerialDenseMatrix locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);
  Epetra_SerialDenseVector values(shapes_->ndofs_);

  double norm_p = 0.0, area = 0.0;
  for (unsigned int q = 0; q < shapes_->nqpoints_; ++q)
  {
    double numerical_p = 0.0;
    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      numerical_p += shapes_->shfunct(i, q) * interiorPressnp_(i);
    norm_p += numerical_p * shapes_->jfac(q);
    area += shapes_->jfac(q);
  }
  elevec1((nsd_ + 2) * nen_) = norm_p / area;

  for (unsigned int i = 0; i < nen_; ++i)
  {
    // evaluate shape polynomials in node
    for (unsigned int idim = 0; idim < nsd_; idim++)
      shapes_->xsi(idim) = locations(idim, i);
    shapes_->polySpace_->Evaluate(shapes_->xsi, values);

    // compute values for velocity and pressure by summing over all basis functions
    double sum = 0;
    for (unsigned int k = 0; k < shapes_->ndofs_; ++k)
    {
      sum += values(k) * interiorPressnp_(k);
    }
    elevec1(nsd_ * nen_ + i) = sum;

    for (unsigned int d = 0; d < nsd_; ++d)
    {
      sum = 0.0;
      for (unsigned int k = 0; k < shapes_->ndofs_; ++k)
        sum += values(k) * interiorVelnp_(d * shapes_->ndofs_ + k);
      elevec1(d * nen_ + i) = sum;
    }
  }

  if (!padaptivity)
  {
    Epetra_SerialDenseVector fvalues(1);
    Epetra_SerialDenseVector touchcount(nen_);
    int sumindex = 0;
    for (unsigned int face = 0; face < nfaces_; ++face)
    {
      DRT::UTILS::ShapeValuesFaceParams svfparams(
          ele->Faces()[face]->Degree(),
          usescompletepoly_,
          2 * ele->Faces()[face]->Degree());
      localSolver_->shapesface_ = DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
      localSolver_->shapesface_->EvaluateFace(*ele, face, *shapes_);

      fvalues.Resize(localSolver_->shapesface_->nfdofs_);

      // const int * fnodeIds = ele->Faces()[face]->NodeIds();
      for (int i = 0;i < DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace; ++i)
      {
        // evaluate shape polynomials in node
        for (unsigned int idim = 0; idim < nsd_ - 1; idim++)
          localSolver_->shapesface_->xsi(idim) = locations(idim, i);

        localSolver_->shapesface_->polySpace_->Evaluate(localSolver_->shapesface_->xsi, fvalues); // TODO: fix face orientation here

        // compute values for velocity and pressure by summing over all basis functions
        double sum = 0.0;
        for (unsigned int k = 0; k < localSolver_->shapesface_->nfdofs_; ++k)
          sum += fvalues(k) * traceVal_[sumindex + k];

        elevec1((nsd_ + 1) * nen_ + localSolver_->shapesface_->faceNodeOrder[face][i]) += sum;
        touchcount(localSolver_->shapesface_->faceNodeOrder[face][i])++;
      }
      sumindex += localSolver_->shapesface_->nfdofs_;
    }

    for (unsigned int i = 0; i < nen_; ++i)
      elevec1((nsd_ + 1) * nen_ + i) /= touchcount(i);
  }
  else
    for (unsigned int i = 0; i < nen_; ++i)
      elevec1((nsd_ + 1) * nen_ + i) = ele->elenodeTrace_(i);

  return;
} // NodeBasedValues

/*----------------------------------------------------------------------*
 * NodeBasedPsi
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::NodeBasedPsi(
    const Teuchos::RCP<MAT::Material> &mat, DRT::ELEMENTS::Acou* ele,
    DRT::Discretization& discretization, const std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, double dt)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::NodeBasedValues");

  dsassert(elevec1.M() == int(nen_), "Vector does not have correct size");
  elevec1.Scale(0.0);
  Epetra_SerialDenseMatrix locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);
  Epetra_SerialDenseVector values(shapes_->ndofs_);

  // calculate mass matrix
  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  // double rho = actmat->Density();
  double c = actmat->SpeedofSound();

  for (unsigned int i = 0; i < nen_; ++i)
  {
    // evaluate shape polynomials in node
    for (unsigned int idim = 0; idim < nsd_; idim++)
      shapes_->xsi(idim) = locations(idim, i);
    shapes_->polySpace_->Evaluate(shapes_->xsi, values);

    // compute values for velocity and pressure by summing over all basis functions
    double sum = 0;
    for (unsigned int k = 0; k < shapes_->ndofs_; ++k)
      sum += values(k) * interiorPressnp_(k);

    elevec1(i) = sum / c / c / dt;
  }

  return;
}

/*----------------------------------------------------------------------*
 * Instance
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::AcouEleCalc<distype> *
DRT::ELEMENTS::AcouEleCalc<distype>::Instance(bool create)
{
  static AcouEleCalc<distype> * instance;
  if (create)
  {
    if (instance == NULL)
      instance = new AcouEleCalc<distype>();
  }
  else
  {
    if (instance != NULL)
      delete instance;
    instance = NULL;
  }
  return instance;
}

/*----------------------------------------------------------------------*
 * Done
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(false);
}

/*----------------------------------------------------------------------*
 * Constructor LocalSolver
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::LocalSolver(
    DRT::UTILS::ShapeValues<distype> &shapeValues) :
    ndofs_(shapeValues.ndofs_), shapes_(shapeValues)
{
  // shape all matrices
  reshapeMatrixIfNecessary(Amat, nsd_ * ndofs_, nsd_ * ndofs_);
  reshapeMatrixIfNecessary(invAmat, nsd_ * ndofs_, nsd_ * ndofs_);
  reshapeMatrixIfNecessary(Bmat, nsd_ * ndofs_, ndofs_);
  reshapeMatrixIfNecessary(Mmat, ndofs_, ndofs_);
  reshapeMatrixIfNecessary(Dmat, ndofs_, ndofs_);
  reshapeMatrixIfNecessary(Hmat, ndofs_, nsd_ * ndofs_);
  reshapeMatrixIfNecessary(massPart, ndofs_, shapeValues.nqpoints_);
}

/*----------------------------------------------------------------------*
 * FaceSpecificConstruction
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::FaceSpecificConstruction(
    const DRT::ELEMENTS::Acou* ele, bool completepoly)
{
  int onfdofs = 0;
  int nfdofs = 0;
  for (unsigned int i = 0; i < nfaces_; ++i)
  {
    DRT::UTILS::PolynomialSpaceParams params(DRT::UTILS::DisTypeToFaceShapeType<distype>::shape,ele->Faces()[i]->Degree(), completepoly);
    nfdofs = DRT::UTILS::PolynomialSpaceCache<nsd_ - 1>::Instance().Create(params)->Size();
    onfdofs += nfdofs;
  }

  reshapeMatrixIfNecessary(Cmat, nsd_ * ndofs_, onfdofs);
  reshapeMatrixIfNecessary(Emat, ndofs_, onfdofs);
  reshapeMatrixIfNecessary(Gmat, onfdofs, onfdofs);
  reshapeMatrixIfNecessary(Imat, onfdofs, nsd_ * ndofs_);
  reshapeMatrixIfNecessary(Jmat, onfdofs, ndofs_);
}

/*----------------------------------------------------------------------*
 * UpdateInteriorVariablesAndComputeResidual
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::UpdateInteriorVariablesAndComputeResidual(
    DRT::Discretization & discretization, Teuchos::ParameterList& params,
    DRT::ELEMENTS::Acou & ele, const Teuchos::RCP<MAT::Material> &mat,
    Epetra_SerialDenseVector & elevec, double dt, bool errormaps,
    bool updateonly)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::UpdateInteriorVariablesAndComputeResidual");

  int onfdofs = elevec.M();
  double theta = 1.0;
  if (dyna_ == INPAR::ACOU::acou_trapezoidal)
    theta = 0.66;
  int stage = -1;
  int dirk_q = -1;

  Epetra_SerialDenseVector tempVelnp;
  Epetra_SerialDenseVector tempPressnp;
  if (dyna_ == INPAR::ACOU::acou_bdf2)
  {
    tempVelnp.Shape(shapes_->ndofs_ * nsd_, 1);
    tempVelnp = interiorVelnp_;
    tempPressnp.Shape(shapes_->ndofs_, 1);
    tempPressnp = interiorPressnp_;
    for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_; ++i)
      interiorVeln_[i] = interiorVelnp_[i] * 4.0 / 3.0 - interiorVelnm_[i] / 3.0;
    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      interiorPressn_[i] = interiorPressnp_[i] * 4.0 / 3.0 - interiorPressnm_[i] / 3.0;
    dt *= 2.0 / 3.0;
  }
  else if (dyna_ == INPAR::ACOU::acou_bdf3)
  {
    tempVelnp.Shape(shapes_->ndofs_ * nsd_, 1);
    tempVelnp = interiorVelnp_;
    tempPressnp.Shape(shapes_->ndofs_, 1);
    tempPressnp = interiorPressnp_;
    for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_; ++i)
      interiorVeln_[i] = interiorVelnp_[i] * 18.0 / 11.0 - interiorVelnm_[i] * 9.0 / 11.0 + interiorVelnmm_[i] * 2.0 / 11.0;
    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      interiorPressn_[i] = interiorPressnp_[i] * 18.0 / 11.0 - interiorPressnm_[i] * 9.0 / 11.0 + interiorPressnmm_[i] * 2.0 / 11.0;
    dt *= 6.0 / 11.0;
  }
  else if (dyna_ == INPAR::ACOU::acou_bdf4)
  {
    tempVelnp.Shape(shapes_->ndofs_ * nsd_, 1);
    tempVelnp = interiorVelnp_;
    tempPressnp.Shape(shapes_->ndofs_, 1);
    tempPressnp = interiorPressnp_;
    for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_; ++i)
      interiorVeln_[i] = interiorVelnp_[i] * 48.0 / 25.0 - interiorVelnm_[i] * 36.0 / 25.0 + interiorVelnmm_[i] * 16.0 / 25.0 - interiorVelnmmm_[i] * 3.0 / 25.0;
    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      interiorPressn_[i] = interiorPressnp_[i] * 48.0 / 25.0 - interiorPressnm_[i] * 36.0 / 25.0 + interiorPressnmm_[i] * 16.0 / 25.0 - interiorPressnmmm_[i] * 3.0 / 25.0;
    dt *= 12.0 / 25.0;
  }
  else if (dyna_ == INPAR::ACOU::acou_dirk23 || dyna_ == INPAR::ACOU::acou_dirk33 ||
           dyna_ == INPAR::ACOU::acou_dirk34 || dyna_ == INPAR::ACOU::acou_dirk54)
  {
    // do the dirk
    stage = params.get<int>("stage");
    interiorVeln_ = ele.elesv_[stage];
    interiorPressn_ = ele.elesp_[stage];
    interiorVeln_.Scale(dt);
    interiorPressn_.Scale(dt);
  }
  else
  {
    interiorVeln_ = interiorVelnp_;
    interiorPressn_ = interiorPressnp_;
  }

  Epetra_SerialDenseVector traceVal_SDV(onfdofs);
  for (int i = 0; i < onfdofs; ++i)
    traceVal_SDV(i) = traceVal_[i];

  Epetra_SerialDenseVector traceVal_SDV_m(onfdofs);
  if (dyna_ == INPAR::ACOU::acou_trapezoidal)
  {
    for (int i = 0; i < onfdofs; ++i)
      traceVal_SDV_m(i) = traceValm_[i];
  }

  // *****************************************************
  // update interior variables first
  // *****************************************************

  Epetra_SerialDenseVector tempVec1(shapes_->ndofs_ * nsd_);
  tempVec1.Multiply('N', 'N', 1.0, localSolver_->Amat, interiorVeln_, 0.0);
  if (theta != 1.0)
  {
    tempVec1.Multiply('N', 'N', -(1.0 - theta), localSolver_->Bmat,interiorPressn_, 1.0);
    tempVec1.Multiply('N', 'N', -(1.0 - theta), localSolver_->Cmat,traceVal_SDV_m, 1.0);
  }
  tempVec1.Multiply('N', 'N', -theta, localSolver_->Cmat, traceVal_SDV, 1.0);

  Epetra_SerialDenseVector tempVec2(shapes_->ndofs_);

  tempVec2.Multiply('N', 'N', 1.0, localSolver_->Mmat, interiorPressn_, 0.0);
  if (theta != 1.0)
  {
    tempVec2.Multiply('N', 'N', -(1.0 - theta), localSolver_->Hmat,interiorVeln_, 1.0);
    tempVec2.Multiply('N', 'N', -(1.0 - theta), localSolver_->Dmat,interiorPressn_, 1.0);
    tempVec2.Multiply('N', 'N', -(1.0 - theta), localSolver_->Emat,traceVal_SDV_m, 1.0);
  }
  tempVec2.Multiply('N', 'N', -(theta), localSolver_->Emat, traceVal_SDV, 1.0);

  // now, we have to do the Schur complement thing, just as in "CondenseLocalPart" but without C and E
  Epetra_SerialDenseMatrix tempMat1(shapes_->ndofs_, nsd_ * shapes_->ndofs_);
  tempMat1.Multiply('N', 'N', theta, localSolver_->Hmat, localSolver_->invAmat,0.0);
  tempVec2.Multiply('N', 'N', -1.0, tempMat1, tempVec1, 1.0);

  Epetra_SerialDenseMatrix tempMat2(shapes_->ndofs_, shapes_->ndofs_);
  tempMat2 = localSolver_->Dmat;
  tempMat2.Scale(theta);
  tempMat2.Multiply('N', 'N', -theta, tempMat1, localSolver_->Bmat, 1.0);
  tempMat2 += localSolver_->Mmat;

  {
    Epetra_SerialDenseSolver inverseDmat;
    inverseDmat.SetMatrix(tempMat2);
    int err = inverseDmat.Invert();
    if (err != 0)
      dserror("Inversion of temporary matrix for Schur complement failed with errorcode %d",err);
  }

  interiorPressnp_.Multiply('N', 'N', 1.0, tempMat2, tempVec2, 0.0);
  tempVec1.Multiply('N', 'N', -theta, localSolver_->Bmat, interiorPressnp_,1.0);
  interiorVelnp_.Multiply('N', 'N', 1.0, localSolver_->invAmat, tempVec1, 0.0);

  // tell this change in the interior variables the discretization
  bool padaptivity = params.get<bool>("padaptivity");

  if (dyna_ == INPAR::ACOU::acou_dirk23 || dyna_ == INPAR::ACOU::acou_dirk33 ||
      dyna_ == INPAR::ACOU::acou_dirk34 || dyna_ == INPAR::ACOU::acou_dirk54)
  {
    double dirk_a[6][6];
    double dirk_b[6];
    ACOU::FillDIRKValues(dyna_, dirk_a, dirk_b, dirk_q);

    ele.eleyp_[stage] = interiorPressnp_;
    ele.eleyv_[stage] = interiorVelnp_;

    ele.elefp_[stage] = ele.eleinteriorPressn_;
    ele.elefp_[stage].Scale(-1.0);
    ele.elefp_[stage] += ele.eleyp_[stage];
    ele.elefp_[stage].Scale(1.0 / dt); // dt includes a_ii
    ele.elefv_[stage] = ele.eleinteriorVeln_;
    ele.elefv_[stage].Scale(-1.0);
    ele.elefv_[stage] += ele.eleyv_[stage];
    ele.elefv_[stage].Scale(1.0 / dt);

    Epetra_SerialDenseVector tempVecp1(shapes_->ndofs_);
    Epetra_SerialDenseVector tempVecv1(shapes_->ndofs_ * nsd_);
    for (int i = 0; i < stage; ++i)
    {
      tempVecp1 = ele.elefp_[i];
      tempVecp1.Scale(-dirk_a[stage][i] / dirk_a[stage][stage]);
      ele.elefp_[stage] += tempVecp1;
      tempVecv1 = ele.elefv_[i];
      tempVecv1.Scale(-dirk_a[stage][i] / dirk_a[stage][stage]);
      ele.elefv_[stage] += tempVecv1;
    }
    tempVecp1 = ele.elefp_[stage];
    tempVecp1.Scale(dt * dirk_b[stage] / dirk_a[stage][stage]);
    ele.eleinteriorPressnp_ += tempVecp1;
    tempVecv1 = ele.elefv_[stage];
    tempVecv1.Scale(dt * dirk_b[stage] / dirk_a[stage][stage]);
    ele.eleinteriorVelnp_ += tempVecv1;

    if (stage == dirk_q - 1)
    {
      ele.eleinteriorPressn_ = ele.eleinteriorPressnp_;
      ele.eleinteriorVeln_ = ele.eleinteriorVelnp_;
    }
  } // if(dyna_ == INPAR::ACOU::acou_dirk?? )
  else if (!padaptivity)
  {
    ele.eleinteriorPressnp_ = interiorPressnp_;
    ele.eleinteriorVelnp_ = interiorVelnp_;
  }

  // *****************************************************
  // local postprocessing to calculate error maps
  // *****************************************************

  if (errormaps)
  {
    const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
    double rho = actmat->Density();

    // first step: calculate the gradient we need (bold p in paper by Nguyen)
    Epetra_SerialDenseVector temp(shapes_->ndofs_ * nsd_);
    temp.Multiply('N', 'N', 1.0, localSolver_->Bmat, interiorPressnp_, 0.0);
    temp.Multiply('N', 'N', 1.0, localSolver_->Cmat, traceVal_SDV, 1.0);
    Epetra_SerialDenseVector p(shapes_->ndofs_ * nsd_);
    p.Multiply('N', 'N', rho / dt, localSolver_->invAmat, temp, 0.0);

    // second step: postprocess the pressure field: therefore we need some additional matrices!
    unsigned int ndofspost = 1;
    for (unsigned int i = 0; i < nsd_; ++i)
      ndofspost *= (ele.Degree() + 2);

    Epetra_SerialDenseMatrix H(ndofspost, ndofspost);
    Epetra_SerialDenseVector R(ndofspost);
    double err_p = CalculateError(ele, H, R, p);

    Teuchos::RCP<std::vector<double> > values = params.get<Teuchos::RCP<std::vector<double> > >("elevals");

    if ((stage == -1 || stage == dirk_q - 1) && padaptivity) // time integrators or last stage of DIRK
    {
      double padaptivitytol = params.get<double>("padaptivitytol");
      double delta_k = 0.0;
      if (err_p != 0.0)
        delta_k = std::ceil(std::log10(err_p / padaptivitytol));

      int newdeg = (ele.Degree() + delta_k);
      if (newdeg < 1)
        newdeg = 1;
      if (newdeg > 6)
        newdeg = 6;
      (*values)[ele.Id()] = newdeg;

      // projection of the internal field (cannot yet project trace field, since the
      // face degree depends on all neighbors! (also not needed for DIRK)
      ProjectPadapField(ele, newdeg);
      updateonly = true;
      return;
    }
    else
      (*values)[ele.Id()] = err_p;
  }

  if (updateonly)
    return;

  // *****************************************************
  // compute residual second (reuse intermediate matrices)
  // *****************************************************
  if (dyna_ == INPAR::ACOU::acou_dirk23 || dyna_ == INPAR::ACOU::acou_dirk33 ||
      dyna_ == INPAR::ACOU::acou_dirk34 || dyna_ == INPAR::ACOU::acou_dirk54)
  {
    double dirk_a[6][6];
    double dirk_b[6];
    ACOU::FillDIRKValues(dyna_, dirk_a, dirk_b, dirk_q);

    stage++;
    if (stage == dirk_q)
    {
      stage = 0;
      ele.eleinteriorPressnp_ = ele.eleinteriorPressn_;
      ele.eleinteriorVelnp_ = ele.eleinteriorVeln_;
    }

    Epetra_SerialDenseVector tempVecp1(shapes_->ndofs_);
    Epetra_SerialDenseVector tempVecp2(shapes_->ndofs_);
    Epetra_SerialDenseVector tempVecv1(shapes_->ndofs_ * nsd_);
    Epetra_SerialDenseVector tempVecv2(shapes_->ndofs_ * nsd_);

    ele.elesp_[stage] = ele.eleinteriorPressn_;
    ele.elesp_[stage].Scale(1.0 / dt); // dt includes a_ii
    ele.elesv_[stage] = ele.eleinteriorVeln_;
    ele.elesv_[stage].Scale(1.0 / dt); // dt includes a_ii

    for (int i = 0; i < stage; ++i)
    {
      tempVecp1 = ele.elesp_[i];
      tempVecp1.Scale(-1.0);
      tempVecp2 = ele.eleyp_[i];
      tempVecp2.Scale(1.0 / dt);
      tempVecp1 += tempVecp2;
      tempVecp1.Scale(dirk_a[stage][i] / dirk_a[stage][stage]);
      ele.elesp_[stage] += tempVecp1;
      tempVecv1 = ele.elesv_[i];
      tempVecv1.Scale(-1.0);
      tempVecv2 = ele.eleyv_[i];
      tempVecv2.Scale(1.0 / dt);
      tempVecv1 += tempVecv2;
      tempVecv1.Scale(dirk_a[stage][i] / dirk_a[stage][stage]);
      ele.elesv_[stage] += tempVecv1;
    }

    // these vectors s are used for the calculation of the residual -> write them to used local solver variable
    interiorPressnp_ = ele.elesp_[stage];
    interiorVelnp_ = ele.elesv_[stage];
    interiorPressnp_.Scale(dt);
    interiorVelnp_.Scale(dt);
  }
  else if (dyna_ == INPAR::ACOU::acou_bdf2)
  {
    for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_; ++i)
      interiorVelnp_[i] = 4.0 / 3.0 * interiorVelnp_[i] - 1.0 / 3.0 * tempVelnp[i];
    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      interiorPressnp_[i] = 4.0 / 3.0 * interiorPressnp_[i] - 1.0 / 3.0 * tempPressnp[i];
    ele.eleinteriorPressnm_ = tempPressnp;
    ele.eleinteriorVelnm_ = tempVelnp;
  }
  else if (dyna_ == INPAR::ACOU::acou_bdf3)
  {
    for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_; ++i)
      interiorVelnp_[i] = 18.0 / 11.0 * interiorVelnp_[i] - 9.0 / 11.0 * tempVelnp[i] + 2.0 / 11.0 * interiorVelnm_[i];
    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      interiorPressnp_[i] = 18.0 / 11.0 * interiorPressnp_[i] - 9.0 / 11.0 * tempPressnp[i] + 2.0 / 11.0 * interiorPressnm_[i];
    ele.eleinteriorPressnmm_ = ele.eleinteriorPressnm_;
    ele.eleinteriorPressnm_ = tempPressnp;
    ele.eleinteriorVelnmm_ = ele.eleinteriorVelnm_;
    ele.eleinteriorVelnm_ = tempVelnp;
  }
  else if (dyna_ == INPAR::ACOU::acou_bdf4)
  {
    for (unsigned int i = 0; i < shapes_->ndofs_ * nsd_; ++i)
      interiorVelnp_[i] = 48.0 / 25.0 * interiorVelnp_[i] - 36.0 / 25.0 * tempVelnp[i] + 16.0 / 25.0 * interiorVelnm_[i] - 3.0 / 25.0 * interiorVelnmm_[i];
    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      interiorPressnp_[i] = 48.0 / 25.0 * interiorPressnp_[i] - 36.0 / 25.0 * tempPressnp[i] + 16.0 / 25.0 * interiorPressnm_[i] - 3.0 / 25.0 * interiorPressnmm_[i];
    ele.eleinteriorPressnmmm_ = ele.eleinteriorPressnmm_;
    ele.eleinteriorPressnmm_ = ele.eleinteriorPressnm_;
    ele.eleinteriorPressnm_ = tempPressnp;
    ele.eleinteriorVelnmmm_ = ele.eleinteriorVelnmm_;
    ele.eleinteriorVelnmm_ = ele.eleinteriorVelnm_;
    ele.eleinteriorVelnm_ = tempVelnp;
  }

  tempVec1.Resize(shapes_->ndofs_);
  tempVec1.Multiply('N', 'N', 1.0, localSolver_->Mmat, interiorPressnp_, 0.0);
  if (theta != 1.0)
  {
    tempVec1.Multiply('N', 'N', -(1.0 - theta), localSolver_->Hmat, interiorVelnp_, 1.0);
    tempVec1.Multiply('N', 'N', -(1.0 - theta), localSolver_->Dmat, interiorPressnp_, 1.0);
    tempVec1.Multiply('N', 'N', -(1.0 - theta), localSolver_->Emat, traceVal_SDV, 1.0);
  }

  Epetra_SerialDenseVector f(shapes_->ndofs_ * nsd_);
  f.Multiply('N', 'N', 1.0, localSolver_->Amat, interiorVelnp_, 0.0);
  if (theta != 1.0)
  {
    f.Multiply('N', 'N', -(1.0 - theta), localSolver_->Bmat, interiorPressnp_,1.0);
    f.Multiply('N', 'N', -(1.0 - theta), localSolver_->Cmat, traceVal_SDV, 1.0);
  }
  tempVec1.Multiply('N', 'N', -1.0, tempMat1, f, 1.0);

  tempVec2.Resize(shapes_->ndofs_);
  tempVec2.Multiply('N', 'N', 1.0, tempMat2, tempVec1, 0.0);

  elevec.Multiply('N', 'N', -theta, localSolver_->Jmat, tempVec2, 0.0);

  tempVec1.Resize(shapes_->ndofs_ * nsd_);
  tempVec1.Multiply('N', 'N', -theta, localSolver_->Bmat, tempVec2, 0.0);
  tempVec1 += f;
  tempVec2.Resize(shapes_->ndofs_ * nsd_);
  tempVec2.Multiply('N', 'N', 1.0, localSolver_->invAmat, tempVec1, 0.0);

  elevec.Multiply('N', 'N', -theta, localSolver_->Imat, tempVec2, 1.0);

  if (theta != 1.0)
  {
    elevec.Multiply('N', 'N', -(1.0 - theta), localSolver_->Imat, interiorVelnp_, 1.0);
    elevec.Multiply('N', 'N', -(1.0 - theta), localSolver_->Jmat, interiorPressnp_, 1.0);
    elevec.Multiply('N', 'N', -(1.0 - theta), localSolver_->Gmat, traceVal_SDV, 1.0);
  }
  return;
} // UpdateInteriorVariablesAndComputeResidual

/*----------------------------------------------------------------------*
 * ProjectPadapField
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::ProjectPadapField(
    DRT::ELEMENTS::Acou & ele, int newdeg)
{

  // create shape values for projection
  Teuchos::RCP<DRT::UTILS::ShapeValues<distype> > proj_shapes = Teuchos::rcp(
      new DRT::UTILS::ShapeValues<distype>(newdeg, usescompletepoly_,2 * newdeg));
  proj_shapes->Evaluate(ele);

  // calculate mass matrix with new degree
  Epetra_SerialDenseMatrix massPart(proj_shapes->ndofs_,proj_shapes->nqpoints_);
  Epetra_SerialDenseMatrix massMat(proj_shapes->ndofs_, proj_shapes->ndofs_);
  for (unsigned int q = 0; q < proj_shapes->nqpoints_; ++q)
  {
    const double sqrtfac = std::sqrt(proj_shapes->jfac(q));
    for (unsigned int i = 0; i < proj_shapes->ndofs_; ++i)
      massPart(i, q) = proj_shapes->shfunct(i, q) * sqrtfac;
  }
  for (unsigned int i = 0; i < proj_shapes->ndofs_; ++i)
    for (unsigned int j = 0; j < proj_shapes->ndofs_; ++j)
    {
      double sum = 0.0;
      for (unsigned int k = 0; k < proj_shapes->nqpoints_; ++k)
        sum += massPart(i, k) * massPart(j, k);
      massMat(i, j) = sum;
    }
  // invert the mass matrix
  {
    Epetra_SerialDenseSolver inversemat;
    inversemat.SetMatrix(massMat);
    inversemat.Invert();
  }

  // calculate rectangular mass matrix with mixed degree
  LINALG::Matrix<nsd_, 1> xsi;
  Epetra_SerialDenseVector values(shapes_->ndofs_);
  Epetra_SerialDenseMatrix rectMat(proj_shapes->ndofs_, shapes_->ndofs_);

  for (unsigned int i = 0; i < proj_shapes->ndofs_; ++i)
    for (unsigned int k = 0; k < proj_shapes->nqpoints_; ++k)
    {
      // calculate values of original shape functions at new quadrature point
      const double* gpcoord = proj_shapes->quadrature_->Point(k);
      for (unsigned int idim = 0; idim < nsd_; idim++)
        xsi(idim) = gpcoord[idim];

      shapes_->polySpace_->Evaluate(xsi, values);
      for (unsigned int j = 0; j < shapes_->ndofs_; ++j)
      {
        rectMat(i, j) += proj_shapes->jfac(k) * proj_shapes->shfunct(i, k) * values(j);
      }
    }
  Epetra_SerialDenseMatrix rectMat2(proj_shapes->ndofs_, shapes_->ndofs_); // TODO this can be done nicer
  rectMat2.Multiply('N', 'N', 1.0, massMat, rectMat, 0.0);

  // build a nsd_ rectMat
  Epetra_SerialDenseMatrix rectMat3(proj_shapes->ndofs_ * nsd_,
      shapes_->ndofs_ * nsd_);
  for (unsigned d = 0; d < nsd_; ++d)
    for (unsigned int j = 0; j < proj_shapes->ndofs_; ++j)
      for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
        rectMat3(j + proj_shapes->ndofs_ * d, i + shapes_->ndofs_ * d) =
            rectMat2(j, i);

  // in case of bdf, we have to project the history as well
  switch (dyna_)
  {
  case INPAR::ACOU::acou_bdf4:
    reshapeMatrixIfNecessary(ele.eleinteriorPressnmmm_, proj_shapes->ndofs_, 1);
    reshapeMatrixIfNecessary(ele.eleinteriorVelnmmm_, proj_shapes->ndofs_ * nsd_, 1);
    ele.eleinteriorPressnmmm_.Multiply('N', 'N', 1.0, rectMat2, ele.eleinteriorPressnmm_, 0.0);
    ele.eleinteriorVelnmmm_.Multiply('N', 'N', 1.0, rectMat3, ele.eleinteriorVelnmm_, 0.0); // no break here
  case INPAR::ACOU::acou_bdf3:
    reshapeMatrixIfNecessary(ele.eleinteriorPressnmm_, proj_shapes->ndofs_, 1);
    reshapeMatrixIfNecessary(ele.eleinteriorVelnmm_, proj_shapes->ndofs_ * nsd_,1);
    ele.eleinteriorPressnmm_.Multiply('N', 'N', 1.0, rectMat2, ele.eleinteriorPressn_, 0.0);
    ele.eleinteriorVelnmm_.Multiply('N', 'N', 1.0, rectMat3, ele.eleinteriorVeln_, 0.0); // no break here
  case INPAR::ACOU::acou_bdf2:
    reshapeMatrixIfNecessary(ele.eleinteriorPressnm_, proj_shapes->ndofs_, 1);
    reshapeMatrixIfNecessary(ele.eleinteriorVelnm_, proj_shapes->ndofs_ * nsd_, 1);
    ele.eleinteriorPressnm_.Multiply('N', 'N', 1.0, rectMat2, ele.eleinteriorPressnp_, 0.0);
    ele.eleinteriorVelnm_.Multiply('N', 'N', 1.0, rectMat3, ele.eleinteriorVelnp_, 0.0); // no break here
  default:
    reshapeMatrixIfNecessary(ele.eleinteriorPressn_, proj_shapes->ndofs_, 1);
    reshapeMatrixIfNecessary(ele.eleinteriorVeln_, proj_shapes->ndofs_ * nsd_, 1);
    ele.eleinteriorPressn_.Multiply('N', 'N', 1.0, rectMat2, ele.eleinteriorPressnp_, 0.0);
    ele.eleinteriorVeln_.Multiply('N', 'N', 1.0, rectMat3, ele.eleinteriorVelnp_, 0.0);
    break; // there you go
  }

  // save the projected pressure to the element!
  reshapeMatrixIfNecessary(ele.eleinteriorPressnp_, proj_shapes->ndofs_, 1);
  ele.eleinteriorPressnp_.Multiply('N', 'N', 1.0, rectMat2, interiorPressnp_, 0.0);

  // project the velocity field as well!! (reuse matrices)
  reshapeMatrixIfNecessary(ele.eleinteriorVelnp_, proj_shapes->ndofs_ * nsd_, 1);
  ele.eleinteriorVelnp_.Multiply('N', 'N', 1.0, rectMat3, interiorVelnp_, 0.0);

  // calculation of node based trace values, since this is no longer possible afterwards!!
  {
    //ele.elenodeTrace_.Shape(nen_,1);
    reshapeMatrixIfNecessary(ele.elenodeTrace_, nen_, 1);
    Epetra_SerialDenseVector fvalues(1);
    Epetra_SerialDenseVector touchcount(nen_);
    Epetra_SerialDenseMatrix locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);

    int sumindex = 0;
    for (unsigned int face = 0; face < nfaces_; ++face)
    {
      DRT::UTILS::ShapeValuesFaceParams svfparams(
          ele.Faces()[face]->Degree(),
          usescompletepoly_,
          2 * ele.Faces()[face]->Degree());
      localSolver_->shapesface_ =DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
      localSolver_->shapesface_->EvaluateFace(ele, face, *shapes_);
      fvalues.Resize(localSolver_->shapesface_->nfdofs_);

      for (int i = 0; i < DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace; ++i)
      {
        // evaluate shape polynomials in node
        for (unsigned int idim = 0; idim < nsd_ - 1; idim++)
          localSolver_->shapesface_->xsi(idim) = locations(idim, i);
        localSolver_->shapesface_->polySpace_->Evaluate(localSolver_->shapesface_->xsi, fvalues); // TODO: fix face orientation here

        // compute values for velocity and pressure by summing over all basis functions
        double sum = 0.0;
        for (unsigned int k = 0; k < localSolver_->shapesface_->nfdofs_; ++k)
          sum += fvalues(k) * traceVal_[sumindex + k];

        ele.elenodeTrace_(localSolver_->shapesface_->faceNodeOrder[face][i]) += sum;
        touchcount(localSolver_->shapesface_->faceNodeOrder[face][i])++;
      }
      sumindex += localSolver_->shapesface_->nfdofs_;
    }

    for (unsigned int i = 0; i < nen_; ++i)
      ele.elenodeTrace_(i) /= touchcount(i);
  }

  return;
}

/*----------------------------------------------------------------------*
 * CalculateError
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::AcouEleCalc<distype>::CalculateError(
    DRT::ELEMENTS::Acou & ele, Epetra_SerialDenseMatrix & h,
    Epetra_SerialDenseVector & rhs, Epetra_SerialDenseVector & p)
{
  DRT::UTILS::PolynomialSpace<nsd_> postpoly(distype, ele.Degree() + 1, ele.UsesCompletePolynomialSpace());
  LINALG::Matrix<nsd_, 1> xsi;
  int ndofspost = 1;
  for (unsigned int i = 0; i < nsd_; ++i)
    ndofspost *= (ele.Degree() + 2);

  Epetra_SerialDenseMatrix derivs(nsd_, ndofspost);
  Epetra_SerialDenseVector myvalues(shapes_->ndofs_);

  LINALG::Matrix<nsd_, nen_> deriv;
  LINALG::Matrix<nsd_, nsd_> xjm, xji;

  Epetra_SerialDenseVector values(ndofspost);

  Teuchos::RCP<DRT::UTILS::GaussPoints> postquad;
  postquad = DRT::UTILS::GaussPointCache::Instance().Create(distype,(ele.Degree() + 1) * 2);

  LINALG::Matrix<nsd_, nen_> xyze;
  GEO::fillInitialPositionArray<distype, nsd_, LINALG::Matrix<nsd_, nen_> >(&ele, xyze);

  for (int q = 0; q < postquad->NumPoints(); ++q)
  {
    const double* gpcoord = postquad->Point(q);
    for (unsigned int idim = 0; idim < nsd_; idim++)
      xsi(idim) = gpcoord[idim];

    postpoly.Evaluate(xsi, values);
    postpoly.Evaluate_deriv1(xsi, derivs);

    DRT::UTILS::shape_function_deriv1<distype>(xsi, deriv);
    xjm.MultiplyNT(deriv, xyze);
    const double jfac = xji.Invert(xjm) * postquad->Weight(q);

    // transform shape functions derivatives
    for (int i = 0; i < ndofspost; ++i)
    {
      double res[nsd_];
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        res[d] = xji(d, 0) * derivs(0, i);
        for (unsigned int e = 1; e < nsd_; ++e)
          res[d] += xji(d, e) * derivs(e, i);
      }
      for (unsigned int d = 0; d < nsd_; ++d)
        derivs(d, i) = res[d];
    }

    shapes_->polySpace_->Evaluate(xsi, myvalues);

    for (int j = 0; j < ndofspost; ++j)
      h(0, j) += values(j) * jfac;
    for (unsigned int j = 0; j < shapes_->ndofs_; ++j)
      rhs(0) += myvalues(j) * jfac * interiorPressnp_[j];

    for (int i = 1; i < ndofspost; ++i)
    {
      for (int j = 0; j < ndofspost; ++j)
      {
        double t = 0;
        for (unsigned int d = 0; d < nsd_; ++d)
          t += derivs(d, i) * derivs(d, j);
        h(i, j) += t * jfac;
      }
    }
    double ugrad[nsd_];
    for (unsigned int d = 0; d < nsd_; ++d)
      ugrad[d] = 0.0;
    for (unsigned int d = 0; d < nsd_; ++d)
      for (unsigned int j = 0; j < shapes_->ndofs_; ++j)
        ugrad[d] += myvalues(j) * p(j + d * shapes_->ndofs_);
    for (int i = 1; i < ndofspost; ++i)
      for (unsigned int d = 0; d < nsd_; ++d)
        rhs(i) += ugrad[d] * derivs(d, i) * jfac;

  } // for (int q=0; q<postquad->NumPoints(); ++q)

  {
    Epetra_SerialDenseSolver inverseH;
    inverseH.SetMatrix(h);
    inverseH.SetVectors(rhs, rhs);
    inverseH.Solve();
  }

  double err_p = 0.0;
  double numerical_post = 0.0;
  double numerical = 0.0;
  double area = 0.0;
  for (unsigned int q = 0; q < shapes_->nqpoints_; ++q)
  {
    numerical_post = 0.0;
    numerical = 0.0;
    const double* gpcoord = shapes_->quadrature_->Point(q);
    for (unsigned int idim = 0; idim < nsd_; idim++)
      xsi(idim) = gpcoord[idim];

    postpoly.Evaluate(xsi, values);
    for (int i = 0; i < ndofspost; ++i)
      numerical_post += values(i) * rhs(i);
    for (unsigned int i = 0; i < shapes_->ndofs_; ++i)
      numerical += shapes_->shfunct(i, q) * interiorPressnp_(i);
    area += shapes_->jfac(q);
    err_p += (numerical_post - numerical) * (numerical_post - numerical)
        * shapes_->jfac(q);
  } // for (int q=0; q<postquad->NumPoints(); ++q)

  err_p /= area;

  return err_p;
} // FillMatrices

/*----------------------------------------------------------------------*
 * ComputeABCNodeVals
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::ComputePMonNodeVals(
    DRT::ELEMENTS::Acou* ele, Teuchos::ParameterList& params,
    Teuchos::RCP<MAT::Material> & mat, int face,
    Epetra_SerialDenseMatrix &elemat, Epetra_SerialDenseVector &elevec,
    int indexstart)
{

  localSolver_->shapesface_->EvaluateFace(*ele, face, *shapes_);

  Teuchos::RCP<std::vector<int> > indices;
  if (params.isParameter("nodeindices"))
  {
    indices = params.get<Teuchos::RCP<std::vector<int> > >("nodeindices");
    if (int(indices->size()) != elevec.M())
      dserror("Vector does not have correct size for line");
  }
  else if (elevec.M() != DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace)
    dserror("Vector does not have correct size");

  Epetra_SerialDenseMatrix locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);

  Epetra_SerialDenseVector fvalues(localSolver_->shapesface_->nfdofs_);
  LINALG::Matrix<nsd_-1, 1> xsiFl;

  for (int i = 0; i < DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace; ++i)
  {
    // evaluate shape polynomials in node
    for (unsigned int idim = 0; idim < nsd_ - 1; idim++)
      xsiFl(idim) = locations(idim, i);
    localSolver_->shapesface_->polySpace_->Evaluate(xsiFl, fvalues); // TODO: fix face orientation here

    // compute values for velocity and pressure by summing over all basis functions
    double sum = 0.0;
    for (unsigned int k = 0; k < localSolver_->shapesface_->nfdofs_; ++k)
      sum += fvalues(k) * traceVal_[indexstart + k];

    if (params.isParameter("nodeindices"))
    {
      // is node i part of the line?
      for (unsigned int j = 0; j < indices->size(); ++j)
        if ((*indices)[j] == i)
        {
          elevec(j) = sum;
        }
    }
    else
      elevec(i) = sum;
  }

  return;
}

/*----------------------------------------------------------------------*
 * ComputeSourcePressureMonitor
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::ComputeSourcePressureMonitor(
    DRT::ELEMENTS::Acou* ele, Teuchos::ParameterList& params,
    Teuchos::RCP<MAT::Material> & mat, int face,
    Epetra_SerialDenseMatrix &elemat, Epetra_SerialDenseVector &elevec,
    int indexstart)
{

  shapesface_->EvaluateFace(*ele, face,shapes_);

  // get the values for the source term!
  Teuchos::RCP<Epetra_MultiVector> adjointrhs = params.get<Teuchos::RCP<Epetra_MultiVector> >("adjointrhs");
  int step = params.get<int>("step");
  int stepmax = adjointrhs->NumVectors();
  if (step < stepmax)
  {
    bool smoothing = true;
    if (!smoothing)
    {
      // adjointrhs is a nodebased vector. Which nodes do belong to this element and which values are associated to those?
      const int * fnodeIds = ele->Faces()[face]->NodeIds();
      int numfnode = ele->Faces()[face]->NumNode();
      double values[numfnode];

      for (int i = 0; i < numfnode; ++i)
      {
        int localnodeid = adjointrhs->Map().LID(fnodeIds[i]);
        if (localnodeid >= 0)
          values[i] = adjointrhs->operator ()(stepmax - step - 1)->operator [](localnodeid); // in inverse order -> we're integrating backwards in time
        else
          dserror("could not get value from node %d of element %d on proc %d ", fnodeIds[i], ele->Id(), ele->Owner());

      } // for(int i=0; i<numfnode; ++i)

      Epetra_SerialDenseMatrix locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);

      Epetra_SerialDenseVector fvalues(shapesface_->nfdofs_);
      LINALG::Matrix<nsd_ - 1, 1> xsitemp;
      for (int i = 0;i < DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace; ++i)
      {
        // evaluate shape polynomials in node
        for (unsigned int idim = 0; idim < nsd_ - 1; idim++)
          xsitemp(idim) = locations(idim, i); // TODO: fix face orientation here
        shapesface_->polySpace_->Evaluate(xsitemp, fvalues);

        // compute values for velocity and pressure by summing over all basis functions
        for (unsigned int k = 0; k < shapesface_->nfdofs_; ++k)
          elevec(indexstart + k) += fvalues(k) * values[i];
      }
    }
    else
    {
      // adjointrhs is a nodebased vector. Which nodes do belong to this element and which values are associated to those?
      const int * fnodeIds = ele->Faces()[face]->NodeIds();
      int numfnode = ele->Faces()[face]->NumNode();
      double fnodexyz[numfnode][nsd_];
      double values[numfnode];

      for (int i = 0; i < numfnode; ++i)
      {
        int localnodeid = adjointrhs->Map().LID(fnodeIds[i]);
        if (localnodeid >= 0)
          values[i] = adjointrhs->operator ()(stepmax - step - 1)->operator [](localnodeid); // in inverse order -> we're integrating backwards in time
        else
          dserror("could not get value from node %d of element %d on proc %d ",fnodeIds[i], ele->Id(), ele->Owner());
        for (unsigned int d = 0; d < nsd_; ++d)
          fnodexyz[i][d] = shapesface_->xyze(d, i);

      } // for(int i=0; i<numfnode; ++i)

      // now get the nodevalues to the dof values just as in ProjectOpticalField function
      Epetra_SerialDenseMatrix mass(shapesface_->nfdofs_, shapesface_->nfdofs_);
      Epetra_SerialDenseVector trVec(shapesface_->nfdofs_);

      for (unsigned int q = 0; q < shapesface_->nqpoints_; ++q)
      {
        const double fac = shapesface_->jfac(q);
        double xyz[nsd_];
        for (unsigned int d = 0; d < nsd_; ++d)
          xyz[d] = shapesface_->xyzreal(d, q);
        double val = 0.0;

        EvaluateFaceAdjoint(fnodexyz, values, numfnode, xyz, val);

        for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
        {
          // mass matrix
          for (unsigned int j = 0; j < shapesface_->nfdofs_; ++j)
            mass(i, j) += shapesface_->shfunct(i, q) * shapesface_->shfunct(j, q) * fac;
          trVec(i) += shapesface_->shfunct(i, q) * val * fac * double(numfnode) / double(shapesface_->nfdofs_);
        }
      } // for(unsigned int q=0; q<nfdofs_; ++q)

      {
        Epetra_SerialDenseSolver inverseMass;
        inverseMass.SetMatrix(mass);
        inverseMass.SetVectors(trVec, trVec);
        int err = inverseMass.Solve();
        if (err != 0)
          dserror("Inversion of matrix in source computation failed with errorcode %d",err);
      }

      for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
        elevec(indexstart + i) = trVec(i);
    }
    elevec.Scale(-1.0);

  } // if(step<stepmax)

  return;
}

/*----------------------------------------------------------------------*
 * ComputeAbsorbingBC
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::ComputeAbsorbingBC(
    DRT::ELEMENTS::Acou* ele, Teuchos::ParameterList& params,
    Teuchos::RCP<MAT::Material> & mat, int face,
    Epetra_SerialDenseMatrix &elemat, Epetra_SerialDenseVector &elevec,
    int indexstart)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::ComputeAbsorbingBC");

  shapesface_->EvaluateFace(*ele, face,shapes_);

  bool resonly = params.get<bool>("resonly");

  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double c = actmat->SpeedofSound();

  if (!resonly)
  {
    // loop over number of shape functions
    for (unsigned int p = 0; p < shapesface_->nfdofs_; ++p)
    {
      // loop over number of shape functions
      for (unsigned int q = 0; q <= p; ++q)
      {
        double tempG = 0.0;
        for (unsigned int i = 0; i < shapesface_->nqpoints_; ++i)
          tempG += shapesface_->jfac(i) * shapesface_->shfunct(p, i) * shapesface_->shfunct(q, i);

        elemat(indexstart + p, indexstart + q) = elemat(indexstart + q,indexstart + p) -= tempG / c;

      } // for (unsigned int q=0; q<nfdofs_; ++q)
    } // for (unsigned int p=0; p<nfdofs_; ++p)
  } // if(!resonly)

  return;
}

/*----------------------------------------------------------------------*
 * ComputeAbsorbingBC
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::ComputeSourcePressureMonitorLine3D(
    DRT::ELEMENTS::Acou* ele, Teuchos::ParameterList& params,
    Teuchos::RCP<MAT::Material> & mat, int face,
    Epetra_SerialDenseMatrix &elemat, Epetra_SerialDenseVector &elevec,
    int indexstart)
{

  Teuchos::RCP<Epetra_MultiVector> adjointrhs = params.get<Teuchos::RCP<Epetra_MultiVector> >("adjointrhs");
  int step = params.get<int>("step");
  int stepmax = adjointrhs->NumVectors();

  elevec.Scale(0.0);

  if (step < stepmax)
  {
    Teuchos::RCP<std::vector<int> > indices = params.get<Teuchos::RCP<std::vector<int> > >("nodeindices");
    int numlinenode = indices->size();
    double linevalues[numlinenode];

    for (int i = 0; i < numlinenode; ++i)
    {
      int nodeid = ele->Faces()[face]->NodeIds()[(*indices)[i]];
      int lnodeid = adjointrhs->Map().LID(nodeid);
      if (lnodeid >= 0)
        linevalues[i] = adjointrhs->operator ()(stepmax - step - 1)->operator [](lnodeid);
      else
        dserror("could not get value from node %d of element %d on proc %d ",nodeid, ele->Id(), ele->Owner());
    }
    // transform line values to face values, then we can do the standard procedure
    int numfnode = ele->Faces()[face]->NumNode();
    double fnodexyz[numfnode][nsd_];
    double values[numfnode];

    for (int i = 0; i < numfnode; ++i)
    {
      values[i] = 0.0;
      for (unsigned int d = 0; d < nsd_; ++d)
        fnodexyz[i][d] = shapesface_->xyze(d, i);
    }
    for (int i = 0; i < numlinenode; ++i)
      values[(*indices)[i]] = linevalues[i];

    // now get the nodevalues to the dof values just as in ProjectOpticalField function
    Epetra_SerialDenseMatrix mass(shapesface_->nfdofs_, shapesface_->nfdofs_);
    Epetra_SerialDenseVector trVec(shapesface_->nfdofs_);

    for (unsigned int q = 0; q < shapesface_->nqpoints_; ++q)
    {
      const double fac = shapesface_->jfac(q);
      double xyz[nsd_];
      for (unsigned int d = 0; d < nsd_; ++d)
        xyz[d] = shapesface_->xyzreal(d, q);
      double val = 0.0;

      EvaluateFaceAdjoint(fnodexyz, values, numfnode, xyz, val);

      for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
      {
        // mass matrix
        for (unsigned int j = 0; j < shapesface_->nfdofs_; ++j)
          mass(i, j) += shapesface_->shfunct(i, q) * shapesface_->shfunct(j, q) * fac;
        trVec(i, 0) += shapesface_->shfunct(i, q) * val * fac * double(numfnode) / double(shapesface_->nfdofs_);
      }
    } // for(unsigned int q=0; q<nfdofs_; ++q)
    {
      Epetra_SerialDenseSolver inverseMass;
      inverseMass.SetMatrix(mass);
      inverseMass.SetVectors(trVec, trVec);
      int err = inverseMass.Solve();
      if (err != 0)
        dserror("Inversion of matrix in source line 3D failed with errorcode %d",err);
    }
    for (unsigned int i = 0; i < shapesface_->nfdofs_; ++i)
      elevec(indexstart + i) = trVec(i);
    elevec.Scale(-1.0);
  }

  return;
}

/*----------------------------------------------------------------------*
 * EvaluateFaceAdjoint
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::EvaluateFaceAdjoint(
    double fnodexyz[][nsd_], double values[], int numfnode,
    const double(&xyz)[nsd_], double &val) const
{
  const DRT::Element::DiscretizationType facedis =
      DRT::UTILS::DisTypeToFaceShapeType<distype>::shape;
  if (facedis == DRT::Element::line2)
  {
    if (numfnode != 2)
      dserror("number of nodes per face should be 2 for face discretization = line2");
    val = values[1] * sqrt( (fnodexyz[0][0] - xyz[0]) * (fnodexyz[0][0] - xyz[0])
                          + (fnodexyz[0][1] - xyz[1]) * (fnodexyz[0][1] - xyz[1]))
        + values[0] * sqrt( (fnodexyz[1][0] - xyz[0]) * (fnodexyz[1][0] - xyz[0])
                          + (fnodexyz[1][1] - xyz[1]) * (fnodexyz[1][1] - xyz[1]));
    double dist = sqrt( (fnodexyz[0][0] - fnodexyz[1][0]) * (fnodexyz[0][0] - fnodexyz[1][0])
                      + (fnodexyz[0][1] - fnodexyz[1][1]) * (fnodexyz[0][1] - fnodexyz[1][1]));
    val /= dist;
  }
  else if (facedis == DRT::Element::quad4)
  {
    if (numfnode != 4)
      dserror("number of nodes per face should be 4 for face discretization = quad4");
      LINALG::Matrix<4, 4> mat(true);
      for(int i=0; i<4; ++i)
      {
        mat(i,0) = 1.0;
        mat(i,1) = fnodexyz[i][0];
        mat(i,2) = fnodexyz[i][1];
        mat(i,3) = fnodexyz[i][2];
      }
      LINALG::FixedSizeSerialDenseSolver<4,4> inversemat;
      inversemat.SetMatrix(mat);
      int err = inversemat.Invert();
      if(err != 0) dserror("Inversion of 4x4 matrix failed with errorcode %d",err);

      val = 0.0;
      for(int i=0; i<4; ++i)
      {
        LINALG::Matrix<4,1> lvec(true);
        LINALG::Matrix<4,1> rvec(true);
        rvec.PutScalar(0.0);
        rvec(i) = 1.0;
        lvec.Multiply(mat,rvec);
        val += ( lvec(0) + lvec(1) * xyz[0] + lvec(2) * xyz[1] + lvec(3) * xyz[2] ) * values[i];
      }
    }
    else
    dserror("not yet implemented");

  return;
}

/*----------------------------------------------------------------------*
 * ComputeInteriorMatrices
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::ComputeInteriorMatrices(
    const Teuchos::RCP<MAT::Material> &mat, double dt,
    INPAR::ACOU::DynamicType dyna)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::ComputeInteriorMatrices");

  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double rho = actmat->Density();
  double c = actmat->SpeedofSound();

  Epetra_SerialDenseMatrix gradPart(ndofs_ * nsd_, shapes_.nqpoints_);

  // loop quadrature points
  for (unsigned int q = 0; q < shapes_.nqpoints_; ++q)
  {
    const double sqrtfac = std::sqrt(shapes_.jfac(q));
    // loop shape functions
    for (unsigned int i = 0; i < ndofs_; ++i)
    {
      const double valf = shapes_.shfunct(i, q) * sqrtfac;
      massPart(i, q) = valf;
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        const double vald = shapes_.shderxy(i * nsd_ + d, q) * sqrtfac;
        gradPart(d * ndofs_ + i, q) = vald;
      }
    }
  }

  // multiply matrices to perform summation over quadrature points
  Mmat.Multiply('N', 'T', 1.0 / c / c / dt, massPart, massPart, 0.0);
  Dmat.Multiply('N', 'T', rho / dt, massPart, massPart, 0.0); // only temporary for smaller inversion
  Bmat.Multiply('N', 'T', -1.0, gradPart, massPart, 0.0);
  for (unsigned int j = 0; j < ndofs_; ++j)
    for (unsigned int i = 0; i < ndofs_; ++i)
    {
      for (unsigned int d = 0; d < nsd_; ++d)
      {
        Amat(d * ndofs_ + i, d * ndofs_ + j) = rho * c * c * Mmat(i, j);
        Hmat(j, d * ndofs_ + i) = -rho * Bmat(d * ndofs_ + i, j);
      }
    }

  // we will need the inverse of A
  {
    Epetra_SerialDenseSolver inverseAmat;
    inverseAmat.SetMatrix(Dmat);
    int err = inverseAmat.Invert();
    if (err != 0)
      dserror("Inversion for Amat failed with errorcode %d", err);
    for (unsigned int j = 0; j < ndofs_; ++j)
      for (unsigned int i = 0; i < ndofs_; ++i)
        for (unsigned int d = 0; d < nsd_; ++d)
          invAmat(d * ndofs_ + i, d * ndofs_ + j) = Dmat(i, j);
  }
  zeroMatrix(Dmat); // delete this, since we used Dmat only for the calculation of invAmat
  return;
}

/*----------------------------------------------------------------------*
 * ComputeResidual
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::ComputeResidual(
    Epetra_SerialDenseVector & elevec, Epetra_SerialDenseVector & interiorVeln,
    Epetra_SerialDenseVector & interiorPressn, std::vector<double> traceVal,
    INPAR::ACOU::DynamicType dyna)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::ComputeResidual");

  /*
   for implicit Euler
   -1
   +---------+    +---------+
   |         |    |    n-1  |
   n              | A    B  |    | A V     |
   R  = - [ I  J ] |         |    |    n-1  |
   | H   D+M |    | M P     |
   +---------+    +---------+

   for trapezoidal rule
   -1
   +---------+    +----------------------------------------------------------+
   |         |    |    n-1                     n-1          n-1              |
   n              | A    B  |    | A V     - (1 - theta) ( B P   + C Lambda    )            |                    n-1     n-1           n-1
   R  = - [ I  J ] |         |    |    n-1                     n-1       n-1           n-1   | - (1 - theta) ( I V   + J P    + G Lambda    )
   | H   D+M |    | M P     - (1 - theta) ( H V   + DmM P    + E Lambda    ) |
   +---------+    +----------------------------------------------------------+
   */

  int onfdofs = elevec.M();
  Epetra_SerialDenseVector traceVal_SDV(onfdofs);
  if (dyna == INPAR::ACOU::acou_trapezoidal)
    for (int i = 0; i < onfdofs; ++i)
      traceVal_SDV(i) = traceVal[i];

  double theta = 1.0;
  if (dyna == INPAR::ACOU::acou_trapezoidal)
    theta = 0.66;

  Epetra_SerialDenseVector tempVec1(ndofs_);
  tempVec1.Multiply('N', 'N', 1.0, Mmat, interiorPressn, 0.0);

  if (theta != 1.0)
  {
    tempVec1.Multiply('N', 'N', -(1.0 - theta), Hmat, interiorVeln, 1.0);
    tempVec1.Multiply('N', 'N', -(1.0 - theta), Dmat, interiorPressn, 1.0);
    tempVec1.Multiply('N', 'N', -(1.0 - theta), Emat, traceVal_SDV, 1.0);
  }

  Epetra_SerialDenseVector f(ndofs_ * nsd_);
  f.Multiply('N', 'N', 1.0, Amat, interiorVeln, 0.0);
  f.Multiply('N', 'N', -(1.0 - theta), Bmat, interiorPressn, 1.0);
  f.Multiply('N', 'N', -(1.0 - theta), Cmat, traceVal_SDV, 1.0);

  Epetra_SerialDenseMatrix tempMat1(ndofs_, ndofs_ * nsd_);
  tempMat1.Multiply('N', 'N', theta, Hmat, invAmat, 0.0);
  tempVec1.Multiply('N', 'N', -1.0, tempMat1, f, 1.0); // right part of w

  Epetra_SerialDenseMatrix tempMat2(ndofs_, ndofs_);
  tempMat2 = Dmat;
  tempMat2.Scale(theta);
  tempMat2 += Mmat;

  tempMat2.Multiply('N', 'N', -theta, tempMat1, Bmat, 1.0); // = D + M - H A^{-1} B
  {
    Epetra_SerialDenseSolver inverseinW;
    inverseinW.SetMatrix(tempMat2);
    int err = inverseinW.Invert();
    if (err != 0)
      dserror("Inversion of temporary matrix for Schur complement failed with errorcode %d",err);
  }
  // tempMat2 = ( D + M - H A^{-1} B )^{-1}

  Epetra_SerialDenseVector tempVec2(ndofs_);
  tempVec2.Multiply('N', 'N', 1.0, tempMat2, tempVec1, 0.0); // = w = ( D + M - H A^{-1} B )^{-1} ( ... )

  elevec.Multiply('N', 'N', -theta, Jmat, tempVec2, 0.0);

  tempVec1.Resize(ndofs_ * nsd_);
  tempVec1.Multiply('N', 'N', -theta, Bmat, tempVec2, 0.0);
  tempVec1 += f;
  tempVec2.Shape(ndofs_ * nsd_, 1);
  tempVec2.Multiply('N', 'N', 1.0, invAmat, tempVec1, 0.0);

  elevec.Multiply('N', 'N', -theta, Imat, tempVec2, 1.0);
  if (theta != 1.0)
  {
    elevec.Multiply('N', 'N', -(1.0 - theta), Imat, interiorVeln, 1.0);
    elevec.Multiply('N', 'N', -(1.0 - theta), Jmat, interiorPressn, 1.0);
    elevec.Multiply('N', 'N', -(1.0 - theta), Gmat, traceVal_SDV, 1.0);
  }
  return;
} // ComputeResidual

/*----------------------------------------------------------------------*
 * ComputeFaceMatrices
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::ComputeFaceMatrices(
    const int face, const Teuchos::RCP<MAT::Material> & mat, double dt,
    int indexstart)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::ComputeFaceMatrices");

  // Compute the matrices C, E and G
  // Here, we don't consider material properties! Don't forget this during condensation

  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double c = actmat->SpeedofSound();
  double tau = 1.0 / (c) ;/// dt; // / 1.0e6;

  // loop over number of shape functions
  for (unsigned int p = 0; p < shapesface_->nfdofs_; ++p)
  {
    // loop over number of shape functions
    for (unsigned int q = 0; q < ndofs_; ++q)
    {
      // C and E

      // numerical integration: sum over quadrature points
      double tempE = 0.0;
      for (unsigned int i = 0; i < shapesface_->nqpoints_; ++i)
      {
        double temp = shapesface_->jfac(i) * shapesface_->shfunct(p, i)
            * shapesface_->shfunctI(q, i);
        tempE += temp;
        for (unsigned int j = 0; j < nsd_; ++j)
        {
          double temp_d = temp * shapesface_->normals(j, i);
          Cmat(j * ndofs_ + q, indexstart + p) += temp_d;
          Imat(indexstart + p, j * ndofs_ + q) += temp_d;
        }
      }
      Emat(q, indexstart + p) = tempE;
      Jmat(indexstart + p, q) = tempE;

    } // for (unsigned int q=0; q<ndofs_; ++q)
  } // for (unsigned int p=0; p<nfdofs_; ++p)

  // G

  // loop over number of shape functions
  for (unsigned int p = 0; p < shapesface_->nfdofs_; ++p)
  {
    // loop over number of shape functions
    for (unsigned int q = 0; q <= p; ++q)
    {
      double tempG = 0.0;
      for (unsigned int i = 0; i < shapesface_->nqpoints_; ++i)
      {
        tempG += shapesface_->jfac(i) * shapesface_->shfunct(p, i) * shapesface_->shfunct(q, i);
      }
      Gmat(indexstart + p, indexstart + q) = Gmat(indexstart + q,indexstart + p) -= tau * tempG;

    } // for (unsigned int q=0; q<nfdofs_; ++q)
  } // for (unsigned int p=0; p<nfdofs_; ++p)

  // one term is still missing in D!!
  for (unsigned int p = 0; p < ndofs_; ++p)
  {
    for (unsigned int q = 0; q <= p; ++q)
    {
      double tempD = 0.0;
      for (unsigned int i = 0; i < shapesface_->nqpoints_; ++i)
      {
        tempD += shapesface_->jfac(i) * shapesface_->shfunctI(p, i) * shapesface_->shfunctI(q, i);
      }
      Dmat(p, q) = Dmat(q, p) += tau * tempD;
    }
  }

  return;
} // ComputeFaceMatrices

/*----------------------------------------------------------------------*
 * CondenseLocalPart
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::CondenseLocalPart(
    Epetra_SerialDenseMatrix &eleMat, INPAR::ACOU::DynamicType dyna)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::CondenseLocalPart");

  /*
   THE MATRIX
   -1
   +--------+    +-----+
   |        |    |     |
   | A   B  |    |  C  |
   K = - [ I  J ] |        |    |     | + G  = - I x - J y + G
   | H  D+M |    |  E  |
   +--------+    +-----+

   -1     -1           -1
   y = [ D - H  A   B  ]   [ E - H  A   C  ]

   -1
   x = A   [  C - B  y ]

   */

  int onfdofs = eleMat.M();
  double theta = 1.0;
  if (dyna == INPAR::ACOU::acou_trapezoidal)
    theta = 0.66;

  Epetra_SerialDenseMatrix tempMat1(ndofs_, ndofs_ * nsd_);
  tempMat1.Multiply('N', 'N', theta, Hmat, invAmat, 0.0); // =  H A^{-1}

  Epetra_SerialDenseMatrix tempMat2(ndofs_, ndofs_);

  // for nonlinear solves, this is not D+M but D+M(P^i)+dMdPP^i
  tempMat2 = Dmat;
  tempMat2.Scale(theta);
  tempMat2 += Mmat;

  tempMat2.Multiply('N', 'N', -theta, tempMat1, Bmat, 1.0); // = D - H A^{-1} B
  Epetra_SerialDenseMatrix tempMat3(ndofs_, onfdofs);
  tempMat3 = Emat;
  tempMat3.Multiply('N', 'N', -theta, tempMat1, Cmat, theta); // = E - H A^{-1} C

  {
    Epetra_SerialDenseSolver inverseinW;
    inverseinW.SetMatrix(tempMat2);
    int err = inverseinW.Invert();
    if (err != 0)
      dserror("Inversion of temporary matrix for Schur complement failed with errorcode %d", err);
  }
  // tempMat2 = (  D - H A^{-1} B )^{-1}

  eleMat = Gmat; // = G
  eleMat.Scale(theta);
  tempMat1.Shape(ndofs_, onfdofs);
  tempMat1.Multiply('N', 'N', 1.0, tempMat2, tempMat3, 0.0); // = y
  eleMat.Multiply('N', 'N', -theta, Jmat, tempMat1, 1.0); // = - J y + G

  tempMat2.Shape(ndofs_ * nsd_, onfdofs);
  tempMat2 = Cmat;
  tempMat2.Multiply('N', 'N', -theta, Bmat, tempMat1, theta); // = - C^T + B^T W

  tempMat3.Shape(ndofs_ * nsd_, onfdofs);
  tempMat3.Multiply('N', 'N', 1.0, invAmat, tempMat2, 0.0); // = x = A^{-1} ( C - B y )

  eleMat.Multiply('N', 'N', -theta, Imat, tempMat3, 1.0); // = K = G - I x - J y

  return;
} // CondenseLocalPart

/*----------------------------------------------------------------------*
 * Compute internal and face matrices
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::ComputeMatrices(
    const Teuchos::RCP<MAT::Material> &mat, DRT::ELEMENTS::Acou & ele,
    double dt, INPAR::ACOU::DynamicType dyna, bool adjoint)
{
  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double rho = actmat->Density();
  double c = actmat->SpeedofSound();

  // init face matrices
  zeroMatrix(Cmat);
  zeroMatrix(Emat);
  zeroMatrix(Gmat);
  zeroMatrix(Imat);
  zeroMatrix(Jmat);

  switch (dyna)
  {
  case INPAR::ACOU::acou_bdf4:
    dt *= 12.0 / 25.0;
    break;
  case INPAR::ACOU::acou_bdf3:
    dt *= 6.0 / 11.0;
    break;
  case INPAR::ACOU::acou_bdf2:
    dt *= 2.0 / 3.0;
    break;
  default:
    break;
  }
  double tau = 1.0 / (c)  ;/// dt;

  ComputeInteriorMatrices(mat, dt, dyna);

  int sumindex = 0;
  for (unsigned int face = 0; face < nfaces_; ++face)
  {
    DRT::UTILS::ShapeValuesFaceParams svfparams(
        ele.Faces()[face]->Degree(),
        shapes_.usescompletepoly_, 2 * ele.Faces()[face]->Degree());
    shapesface_ = DRT::UTILS::ShapeValuesFaceCache<distype>::Instance().Create(svfparams);
    shapesface_->EvaluateFace(ele, face, shapes_);

    ComputeFaceMatrices(face, mat, dt, sumindex);
    sumindex += shapesface_->nfdofs_;
  }

  // the matrices are slightly different for the adjoint problem
  if (adjoint)
  {
    Bmat.Scale(-1.0 * rho);
    Cmat.Scale(rho);
    Emat.Scale(tau);
    Hmat.Scale(-1.0 / rho);
    Jmat.Scale(-tau);
  }
  else
  {
    Imat.Scale(rho);
    Emat.Scale(-tau);
    Jmat.Scale(tau);
  }

  return;
}

// template classes
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::hex8>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::hex20>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::hex27>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::tet4>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::tet10>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::wedge6>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::quad4>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::quad8>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::quad9>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::tri3>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::tri6>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::nurbs9>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::nurbs27>;
