/*----------------------------------------------------------------------------*/
/*! \file
\brief three dimensional total Lagrange truss element (can be connected to beam3 elements and
adapts assembly automatically according to the thereby changed number of nodal degrees of freedom)

\level 3


*/
/*---------------------------------------------------------------------------*/

#include "truss3.H"

#include "../drt_beam3/beam3eb.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"

#include "../drt_mat/stvenantkirchhoff.H"

#include "../drt_structure_new/str_elements_paramsinterface.H"

/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public) cyron 08/08|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Truss3::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  SetParamsInterfacePtr(params);

  ELEMENTS::ActionType act = ELEMENTS::none;

  if (IsParamsInterface())
  {
    act = ParamsInterface().GetActionType();
  }
  else  // Todo remove as soon as old structural time integration is gone
  {
    // get the action required
    std::string action = params.get<std::string>("action", "calc_none");
    if (action == "calc_none")
      dserror("No action supplied");
    else if (action == "calc_struct_linstiff")
      act = ELEMENTS::struct_calc_linstiff;
    else if (action == "calc_struct_nlnstiff")
      act = ELEMENTS::struct_calc_nlnstiff;
    else if (action == "calc_struct_internalforce")
      act = ELEMENTS::struct_calc_internalforce;
    else if (action == "calc_struct_linstiffmass")
      act = ELEMENTS::struct_calc_linstiffmass;
    else if (action == "calc_struct_nlnstiffmass")
      act = ELEMENTS::struct_calc_nlnstiffmass;
    else if (action == "calc_struct_nlnstifflmass")
      act = ELEMENTS::struct_calc_nlnstifflmass;
    else if (action == "calc_struct_stress")
      act = ELEMENTS::struct_calc_stress;
    else if (action == "calc_struct_update_istep")
      act = ELEMENTS::struct_calc_update_istep;
    else if (action == "calc_struct_reset_istep")
      act = ELEMENTS::struct_calc_reset_istep;
    else if (action == "calc_struct_ptcstiff")
      act = ELEMENTS::struct_calc_ptcstiff;
    else
    {
      std::cout << action << std::endl;
      dserror("Unknown type of action for Truss3");
    }
  }

  switch (act)
  {
    case ELEMENTS::struct_calc_ptcstiff:
    {
      EvaluatePTC<2, 3, 3>(params, elemat1);

      break;
    }
    /*in case that only linear stiffness matrix is required b3_nlstiffmass is called with zero
       displacement and residual values*/
    case ELEMENTS::struct_calc_linstiff:
    {
      // only nonlinear case implemented!
      dserror("linear stiffness matrix called, but not implemented");

      break;
    }
    // calculate internal energy
    case ELEMENTS::struct_calc_energy:
    {
      // need current global displacement and get them from discretization
      // making use of the local-to-global map lm one can extract current displacement and residual
      // values for each degree of freedom
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

      t3_energy(params, mydisp, &elevec1);

      break;
    }
    // nonlinear stiffness and mass matrix are calculated even if only nonlinear stiffness matrix is
    // required
    case ELEMENTS::struct_calc_nlnstiffmass:
    case ELEMENTS::struct_calc_nlnstifflmass:
    case ELEMENTS::struct_calc_nlnstiff:
    case ELEMENTS::struct_calc_internalforce:
    {
      // need current global displacement and residual forces and get them from discretization
      // making use of the local-to-global map lm one can extract current displacement and residual
      // values for each degree of freedom
      //
      // get element displacements
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      // get residual displacements
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (res == Teuchos::null) dserror("Cannot get state vectors 'residual displacement'");
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res, myres, lm);

      /*first displacement vector is modified for proper element evaluation in case of periodic
       *boundary conditions; in case that no periodic boundary conditions are to be applied the
       *following code line may be ignored or deleted*/
      if (params.get<Teuchos::RCP<Epetra_MultiVector>>("RandomNumbers", Teuchos::null) !=
          Teuchos::null)
        NodeShift<2, 3>(params, mydisp);


      // only if random numbers for Brownian dynamics are passed to element, get element velocities
      std::vector<double> myvel(lm.size());
      if (params.get<Teuchos::RCP<Epetra_MultiVector>>("RandomNumbers", Teuchos::null) !=
          Teuchos::null)
      {
        Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState("velocity");
        DRT::UTILS::ExtractMyValues(*vel, myvel, lm);
      }

      // for engineering strains instead of total lagrange use t3_nlnstiffmass2
      if (act == ELEMENTS::struct_calc_nlnstiffmass)
        t3_nlnstiffmass(params, myvel, mydisp, &elemat1, &elemat2, &elevec1);
      else if (act == ELEMENTS::struct_calc_nlnstifflmass)
      {
        t3_nlnstiffmass(params, myvel, mydisp, &elemat1, &elemat2, &elevec1);
        t3_lumpmass(&elemat2);
      }
      else if (act == ELEMENTS::struct_calc_nlnstiff)
        t3_nlnstiffmass(params, myvel, mydisp, &elemat1, nullptr, &elevec1);
      else if (act == ELEMENTS::struct_calc_internalforce)
        t3_nlnstiffmass(params, myvel, mydisp, nullptr, nullptr, &elevec1);

      break;
    }
    case ELEMENTS::struct_calc_update_istep:
    {
      break;
    }
    case ELEMENTS::struct_postprocess_stress:
    {
      // no stress calculation for postprocess. Does not really make sense!
      dserror("No stress output for Truss3!");

      break;
    }
    case ELEMENTS::struct_calc_recover:
    case ELEMENTS::struct_calc_predict:
    case ELEMENTS::struct_calc_reset_istep:
    case ELEMENTS::struct_calc_stress:
    case ELEMENTS::none:
    case ELEMENTS::struct_calc_internalinertiaforce:
    case ELEMENTS::struct_calc_linstiffmass:
    case ELEMENTS::struct_calc_nlnstiff_gemm:
    case ELEMENTS::struct_calc_thickness:
    case ELEMENTS::struct_calc_eleload:
    case ELEMENTS::struct_calc_fsiload:
    case ELEMENTS::struct_calc_store_istep:
    case ELEMENTS::struct_calc_recover_istep:
    case ELEMENTS::struct_calc_reset_all:
    case ELEMENTS::struct_calc_errornorms:
    case ELEMENTS::struct_postprocess_thickness:
    case ELEMENTS::struct_init_gauss_point_data_output:
    case ELEMENTS::struct_gauss_point_data_output:
    case ELEMENTS::struct_update_prestress:
    case ELEMENTS::analyse_jacobian_determinant:
    case ELEMENTS::inversedesign_update:
    case ELEMENTS::inversedesign_switch:
    case ELEMENTS::multi_readrestart:
    case ELEMENTS::multi_init_eas:
    case ELEMENTS::multi_set_eas:
    case ELEMENTS::multi_calc_dens:
    case ELEMENTS::shell_calc_stc_matrix:
    case ELEMENTS::shell_calc_stc_matrix_inverse:
    case ELEMENTS::struct_calc_stifftemp:
    case ELEMENTS::struct_calc_global_gpstresses_map:
    case ELEMENTS::struct_interpolate_velocity_to_point:
    case ELEMENTS::struct_calc_mass_volume:
    case ELEMENTS::struct_calc_brownianforce:
    case ELEMENTS::struct_calc_brownianstiff:
    case ELEMENTS::struct_poro_calc_fluidcoupling:
    case ELEMENTS::struct_poro_calc_scatracoupling:
    case ELEMENTS::struct_poro_calc_prescoupling:
    case ELEMENTS::struct_calc_addjacPTC:
    case ELEMENTS::struct_create_backup:
    case ELEMENTS::struct_recover_from_backup:
    {
      break;
    }
    default:
    {
      std::cout << "\ncalled element with action type " << ActionType2String(act);
      dserror("Unknown type of action for Truss3");
      break;
    }
  }
  return 0;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public) cyron 03/08|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Truss3::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  dserror("This method needs to be modified for bio-polymer networks!");

  return 0;
}


/*-----------------------------------------------------------------------------------------------------------*
 | Evaluate PTC damping (public) cyron 04/10|
 *----------------------------------------------------------------------------------------------------------*/
template <int nnode, int ndim, int dof>
void DRT::ELEMENTS::Truss3::EvaluatePTC(
    Teuchos::ParameterList& params, Epetra_SerialDenseMatrix& elemat1)
{
  dserror(
      "Truss3::EvaluatePTC is deprecated; if needed adapt parameter handling according to "
      "parameter interface pointer first!");
}

/*--------------------------------------------------------------------------------------*
 | calculation of elastic energy                                             cyron 12/10|
 *--------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::t3_energy(
    Teuchos::ParameterList& params, std::vector<double>& disp, Epetra_SerialDenseVector* intenergy)
{
  intenergy->Resize(1);

  // get the material law
  Teuchos::RCP<const MAT::Material> currmat = Material();
  double ym = 0.0;

  // assignment of material parameters; only St.Venant material is accepted for this truss
  switch (currmat->MaterialType())
  {
    case INPAR::MAT::m_stvenant:  // only linear elastic material supported
    {
      const auto* actmat = static_cast<const MAT::StVenantKirchhoff*>(currmat.get());
      ym = actmat->Youngs();

      break;
    }
    default:
    {
      dserror("unknown or improper type of material law");
      break;
    }
  }

  // current node position (first entries 0 .. 2 for first node, 3 ..5 for second node)
  LINALG::Matrix<6, 1> xcurr;

  // auxiliary vector for both internal force and stiffness matrix: N^T_(,xi)*N_(,xi)*xcurr
  LINALG::Matrix<3, 1> aux;

  // strain
  double epsilon;

  // current nodal position (first
  for (int j = 0; j < 3; ++j)
  {
    xcurr(j) = Nodes()[0]->X()[j] + disp[j];          // first node
    xcurr(j + 3) = Nodes()[1]->X()[j] + disp[3 + j];  // second node
  }

  // computing auxiliary vector aux = 4.0*N^T_{,xi} * N_{,xi} * xcurr
  aux(0) = (xcurr(0) - xcurr(3));
  aux(1) = (xcurr(1) - xcurr(4));
  aux(2) = (xcurr(2) - xcurr(5));

  double lcurr = std::sqrt(std::pow(aux(0), 2) + std::pow(aux(1), 2) + std::pow(aux(2), 2));


  switch (kintype_)
  {
    case tr3_totlag:
    {
      // calculating Green-Lagrange strain epsilon
      epsilon = 0.5 * (std::pow(lcurr / lrefe_, 2) - 1.0);

      // W_int = 1/2*E*A*lrefe*\epsilon^2
      (*intenergy)(0) = 0.5 * (ym * crosssec_ * lrefe_ * std::pow(epsilon, 2));
    }
    break;
    case tr3_engstrain:
    {
      // calculating strain epsilon from node position by scalar product:
      epsilon = (lcurr - lrefe_) / lrefe_;

      // W_int = 1/2*E*A*lrefe*\epsilon^2
      (*intenergy)(0) = 0.5 * (ym * crosssec_ * lrefe_ * std::pow(epsilon, 2));
    }
    break;
    default:
      dserror("Unknown type kintype_ for Truss3");
      break;
  }
}

/*--------------------------------------------------------------------------------------*
 | switch between kintypes                                                      tk 11/08|
 *--------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::t3_nlnstiffmass(Teuchos::ParameterList& params,
    std::vector<double>& vel, std::vector<double>& disp, Epetra_SerialDenseMatrix* stiffmatrix,
    Epetra_SerialDenseMatrix* massmatrix, Epetra_SerialDenseVector* force)
{
  /*
   * It is observed that for a mixed problems, such is the case for biopolymer network simulations
   * (), the method "Evaluate" hands in the larger matrices and vectors of size of element described
   * in the input file. For example, if the computational volume contains both Beam and Truss
   * elements. The Evaluate hand into the method a 12x12 matrix. However, for truss element we need
   * only 6x6. Therefore, an appropriate mapping needs to be established to ensure proper assemblies
   * for corresponding DOFs. The algorithm implemented here is valid only for linear elements i.e.
   * element containing two nodes.
   */
  // 6x6 Stiffness Matrix of the Truss
  Epetra_SerialDenseMatrix DummyStiffMatrix;
  DummyStiffMatrix.Shape(6, 6);
  DummyStiffMatrix.Scale(0);
  // 6x6 force vector of the Truss
  Epetra_SerialDenseVector DummyForce;
  DummyForce.Size(6);
  DummyForce.Scale(0);
  // 1x6 velocity vector
  LINALG::Matrix<1, 6> DummyVel;
  DummyVel.Clear();
  // 1x6 displacement vector
  LINALG::Matrix<1, 6> DummyDisp;
  DummyDisp.Clear();
  // Map velocity global level into element level
  if (vel.size() > 12)
    dserror("Vector is larger than 12. Please use different mapping strategy!");
  else if (vel.size() == 6)
  {
    for (int i = 0; i < 6; i++) DummyVel(i) += vel[i];
  }
  else if (vel.size() == 12)
  {
    for (int i = 0; i < 3; i++)
    {
      DummyVel(i) += vel[i];
      DummyVel(i + 3) += vel[i + 6];
    }
  }
  // Map displacement global level into element level
  if (disp.size() > 12)
    dserror("Vector is larger than 12. Please use different mapping strategy!");
  else if (disp.size() == 6)
  {
    for (int i = 0; i < 6; i++) DummyDisp(i) += disp[i];
  }
  else if (disp.size() == 12)
  {
    for (int i = 0; i < 3; i++)
    {
      DummyDisp(i) += disp[i];
      DummyDisp(i + 3) += disp[i + 6];
    }
  }
  switch (kintype_)
  {
    case tr3_totlag:
      t3_nlnstiffmass_totlag(DummyDisp, DummyStiffMatrix, massmatrix, DummyForce);
      break;
    case tr3_engstrain:
      t3_nlnstiffmass_engstr(DummyDisp, DummyStiffMatrix, massmatrix, DummyForce);
      break;
    default:
      dserror("Unknown type kintype_ for Truss3");
      break;
  }

  /*the following function call applies statistical forces and damping matrix according to the
   * fluctuation dissipation theorem; it is dedicated to the application of truss3 elements in the
   * frame of statistical mechanics problems; for these problems a special vector has to be passed
   * to the element packed in the params parameter list; in case that the control routine calling
   * the element does not attach this special vector to params the following method is just doing
   * nothing, which means that for any ordinary problem of structural mechanics it may be ignored*/
  CalcBrownian<2, 3, 3, 3>(params, DummyVel, DummyDisp, DummyStiffMatrix, DummyForce);

  // Map element level into global 12 by 12 element
  if (force->Length() > 12)
    dserror("Vector is larger than 12. Please use different mapping strategy!");
  else if (force->Length() == 6)
  {
    for (int i = 0; i < 6; i++) (*force)(i) += DummyForce(i);
  }
  else if (force->Length() == 12)
  {
    for (int i = 0; i < 3; i++)
    {
      (*force)(i) += DummyForce(i);
      (*force)(i + 6) += DummyForce(i + 3);
    }
  }

  if (stiffmatrix != nullptr)
  {
    // Map element level into global 12 by 12 element
    if (stiffmatrix->RowDim() > 12)
      dserror("Matrix is larger than 12. Please use different mapping strategy!");
    else if (stiffmatrix->RowDim() == 6)
    {
      for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++) (*stiffmatrix)(i, j) += DummyStiffMatrix(i, j);
    }
    else if (stiffmatrix->RowDim() == 12)
    {
      for (int i = 0; i < 3; i++)
      {
        for (int j = 0; j < 3; j++)
        {
          (*stiffmatrix)(i, j) += DummyStiffMatrix(i, j);
          (*stiffmatrix)(i, j + 6) += DummyStiffMatrix(i, j + 3);
          (*stiffmatrix)(i + 6, j + 6) += DummyStiffMatrix(i + 3, j + 3);
          (*stiffmatrix)(i + 6, j) += DummyStiffMatrix(i + 3, j);
        }
      }
    }
  }
}

/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private) cyron 08/08|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::t3_nlnstiffmass_totlag(LINALG::Matrix<1, 6>& DummyDisp,
    Epetra_SerialDenseMatrix& DummyStiffMatrix, Epetra_SerialDenseMatrix* massmatrix,
    Epetra_SerialDenseVector& DummyForce)
{
  // current node position (first entries 0 .. 2 for first node, 3 ..5 for second node)
  LINALG::Matrix<6, 1> xcurr;

  // current nodal position
  for (int j = 0; j < 3; ++j)
  {
    xcurr(j) = Nodes()[0]->X()[j] + DummyDisp(j);          // first node
    xcurr(j + 3) = Nodes()[1]->X()[j] + DummyDisp(3 + j);  // second node
  }

  // calculate force vector and stiffness matrix
  CalcInternalForceStiffTotLag(xcurr, DummyForce, DummyStiffMatrix);

  // get the material law
  Teuchos::RCP<const MAT::Material> currmat = Material();
  double density = 0;

  // assignment of material parameters; only St.Venant material is accepted for this truss
  switch (currmat->MaterialType())
  {
    case INPAR::MAT::m_stvenant:  // only linear elastic material supported
    {
      const auto* actmat = static_cast<const MAT::StVenantKirchhoff*>(currmat.get());
      density = actmat->Density();
    }
    break;
    default:
      dserror("unknown or improper type of material law");
      break;
  }

  // calculating mass matrix
  if (massmatrix != nullptr)
  {
    for (int i = 0; i < 3; ++i)
    {
      (*massmatrix)(i, i) = density * lrefe_ * crosssec_ / 3;
      (*massmatrix)(i + 3, i + 3) = density * lrefe_ * crosssec_ / 3;
      (*massmatrix)(i, i + 3) = density * lrefe_ * crosssec_ / 6;
      (*massmatrix)(i + 3, i) = density * lrefe_ * crosssec_ / 6;
    }
  }
}

/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private) tk 10/08| | engineering strain measure, large
 displacements and rotations                                                |
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::t3_nlnstiffmass_engstr(const LINALG::Matrix<1, 6>& DummyDisp,
    Epetra_SerialDenseMatrix& DummyStiffMatrix, Epetra_SerialDenseMatrix* massmatrix,
    Epetra_SerialDenseVector& DummyForce)
{
  // current node position (first entries 0 .. 2 for first node, 3 ..5 for second node)
  LINALG::Matrix<6, 1> xcurr;

  // Green-Lagrange strain
  double epsilon;

  // auxiliary vector for both internal force and stiffness matrix: N^T_(,xi)*N_(,xi)*xcurr
  LINALG::Matrix<6, 1> aux;

  // current nodal position (first
  for (int j = 0; j < 3; ++j)
  {
    xcurr(j) = Nodes()[0]->X()[j] + DummyDisp(j);          // first node
    xcurr(j + 3) = Nodes()[1]->X()[j] + DummyDisp(3 + j);  // second node
  }

  // computing auxiliary vector aux = 4.0*N^T_{,xi} * N_{,xi} * xcurr
  aux(0) = (xcurr(0) - xcurr(3));
  aux(1) = (xcurr(1) - xcurr(4));
  aux(2) = (xcurr(2) - xcurr(5));
  aux(3) = (xcurr(3) - xcurr(0));
  aux(4) = (xcurr(4) - xcurr(1));
  aux(5) = (xcurr(5) - xcurr(2));

  double lcurr = std::sqrt(std::pow(aux(0), 2) + std::pow(aux(1), 2) + std::pow(aux(2), 2));

  // calculating strain epsilon from node position by scalar product:
  epsilon = (lcurr - lrefe_) / lrefe_;

  // get the material law
  Teuchos::RCP<const MAT::Material> currmat = Material();
  double ym = 0;
  double density = 0;

  // assignment of material parameters; only St.Venant material is accepted for this truss
  switch (currmat->MaterialType())
  {
    case INPAR::MAT::m_stvenant:  // only linear elastic material supported
    {
      const auto* actmat = static_cast<const MAT::StVenantKirchhoff*>(currmat.get());
      ym = actmat->Youngs();
      density = actmat->Density();
    }
    break;
    default:
      dserror("unknown or improper type of material law");
      break;
  }

  // resulting force scaled by current length
  double forcescalar = (ym * crosssec_ * epsilon) / lcurr;

  // computing global internal forces
  for (int i = 0; i < 6; ++i) DummyForce(i) = forcescalar * aux(i);

  // computing linear stiffness matrix
  for (int i = 0; i < 3; ++i)
  {
    // stiffness entries for first node
    DummyStiffMatrix(i, i) = forcescalar;
    DummyStiffMatrix(i, 3 + i) = -forcescalar;
    // stiffness entries for second node
    DummyStiffMatrix(i + 3, i + 3) = forcescalar;
    DummyStiffMatrix(i + 3, i) = -forcescalar;
  }

  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j)
      DummyStiffMatrix(i, j) += (ym * crosssec_ / std::pow(lcurr, 3)) * aux(i) * aux(j);

  // calculating mass matrix.
  if (massmatrix != nullptr)
  {
    for (int i = 0; i < 3; ++i)
    {
      (*massmatrix)(i, i) = density * lrefe_ * crosssec_ / 3;
      (*massmatrix)(i + 3, i + 3) = density * lrefe_ * crosssec_ / 3;
      (*massmatrix)(i, i + 3) = density * lrefe_ * crosssec_ / 6;
      (*massmatrix)(i + 3, i) = density * lrefe_ * crosssec_ / 6;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::CalcInternalForceStiffTotLag(
    const LINALG::Matrix<6, 1>& nodal_positions_totlag, Epetra_SerialDenseVector& forcevec,
    Epetra_SerialDenseMatrix& stiffmat)
{
  /* current nodal displacement (first entries 0 .. 2 for first node, 3 ..5 for second node)
   * compared to reference configuration; note: in general this is not equal to the values in disp
   * since the latter one refers to a nodal displacement compared to a reference configuration
   * before the first time step whereas the following variable refers to the displacement with
   * respect to a reference configuration which may have been set up at any point of time during the
   * simulation (usually this
   * is only important if an element changes its reference position during simulation)*/

  // current displacement = current position - reference position
  LINALG::Matrix<6, 1> ucurr(nodal_positions_totlag);
  ucurr -= X_;

  // auxiliary vector for both internal force and stiffness matrix: N^T_(,xi)*N_(,xi)*xcurr
  LINALG::Matrix<6, 1> aux;

  // computing auxiliary vector aux = N^T_{,xi} * N_{,xi} * nodal_positions_totlag
  aux(0) = nodal_positions_totlag(0) - nodal_positions_totlag(3);
  aux(1) = nodal_positions_totlag(1) - nodal_positions_totlag(4);
  aux(2) = nodal_positions_totlag(2) - nodal_positions_totlag(5);
  aux(3) = nodal_positions_totlag(3) - nodal_positions_totlag(0);
  aux(4) = nodal_positions_totlag(4) - nodal_positions_totlag(1);
  aux(5) = nodal_positions_totlag(5) - nodal_positions_totlag(2);

  // current length
  lcurr_ = std::sqrt(aux(0) * aux(0) + aux(1) * aux(1) + aux(2) * aux(2));

  // calculating strain epsilon from node position by scalar product:
  // epsilon = (xrefe + 0.5*ucurr)^T * N_{,s}^T * N_{,s} * d
  // Green-Lagrange strain ( 1D truss: epsilon = 0.5 (l² - L²)/L² =( 2 L Δl + Δl² )/L² = Δl/L +
  // 0.5(Δl/L)² )
  double epsilon = 0;
  epsilon += (X_(0) + 0.5 * ucurr(0)) * (ucurr(0) - ucurr(3));
  epsilon += (X_(3) + 0.5 * ucurr(3)) * (ucurr(3) - ucurr(0));
  epsilon += (X_(1) + 0.5 * ucurr(1)) * (ucurr(1) - ucurr(4));
  epsilon += (X_(4) + 0.5 * ucurr(4)) * (ucurr(4) - ucurr(1));
  epsilon += (X_(2) + 0.5 * ucurr(2)) * (ucurr(2) - ucurr(5));
  epsilon += (X_(5) + 0.5 * ucurr(5)) * (ucurr(5) - ucurr(2));
  epsilon /= lrefe_ * lrefe_;

  // get the material law
  Teuchos::RCP<const MAT::Material> currmat = Material();
  double ym = 0;

  // assignment of material parameters; only St.Venant material is accepted for this truss
  switch (currmat->MaterialType())
  {
    case INPAR::MAT::m_stvenant:  // only linear elastic material supported
    {
      const auto* actmat = static_cast<const MAT::StVenantKirchhoff*>(currmat.get());
      ym = actmat->Youngs();
      break;
    }

    default:
    {
      dserror("unknown or improper type of material law");
      break;
    }
  }

  eint_ = 0.5 * ym * crosssec_ * lrefe_ * epsilon * epsilon;

  double lrefeinv = 1.0 / lrefe_;
  // computing global internal forces
  for (int i = 0; i < 6; ++i) forcevec(i) = ym * crosssec_ * epsilon * lrefeinv * aux(i);

  // computing linear stiffness matrix
  for (int i = 0; i < 3; ++i)
  {
    // stiffness entries for first node
    stiffmat(i, i) = ym * crosssec_ * epsilon * lrefeinv;
    stiffmat(i, 3 + i) = -1.0 * ym * crosssec_ * epsilon * lrefeinv;
    // stiffness entries for second node
    stiffmat(i + 3, i + 3) = ym * crosssec_ * epsilon * lrefeinv;
    stiffmat(i + 3, i) = -1.0 * ym * crosssec_ * epsilon * lrefeinv;
  }

  double lrefepow3inv = lrefeinv * lrefeinv * lrefeinv;
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j) stiffmat(i, j) += ym * crosssec_ * lrefepow3inv * aux(i) * aux(j);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Truss3::t3_lumpmass(Epetra_SerialDenseMatrix* emass)
{
  // lump mass matrix
  if (emass != nullptr)
  {
    // we assume #elemat2 is a square matrix
    for (int c = 0; c < (*emass).N(); ++c)  // parse columns
    {
      double d = 0.0;
      for (int r = 0; r < (*emass).M(); ++r)  // parse rows
      {
        d += (*emass)(r, c);  // accumulate row entries
        (*emass)(r, c) = 0.0;
      }
      (*emass)(c, c) = d;  // apply sum of row entries on diagonal
    }
  }
}

/*-----------------------------------------------------------------------------------------------------------*
 |computes the number of different random numbers required in each time step for generation of
 stochastic    | |forces; (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Truss3::HowManyRandomNumbersINeed()
{
  /*at each Gauss point one needs as many random numbers as randomly excited degrees of freedom,
   *i.e. three random numbers for the translational degrees of freedom*/
  return (3 * 2);

}  // DRT::ELEMENTS::Beam3::HowManyRandomNumbersINeed

/*-----------------------------------------------------------------------------------------------------------*
 | Assemble stochastic and viscous forces and respective stiffness according to fluctuation
 dissipation      | | theorem (public) cyron 03/10|
 *----------------------------------------------------------------------------------------------------------*/
template <int nnode, int ndim, int dof, int randompergauss>
inline void DRT::ELEMENTS::Truss3::CalcBrownian(Teuchos::ParameterList& params,
    const LINALG::Matrix<1, 6>& DummyVel, const LINALG::Matrix<1, 6>& DummyDisp,
    Epetra_SerialDenseMatrix& DummyStiffMatrix, Epetra_SerialDenseVector& DummyForce)
{
  // if no random numbers for generation of stochastic forces are passed to the element no Brownian
  // dynamics calculations are conducted
  if (params.get<Teuchos::RCP<Epetra_MultiVector>>("RandomNumbers", Teuchos::null) == Teuchos::null)
    return;

  dserror(
      "Truss3::CalcBrownian is deprecated; if needed adapt parameter handling according to "
      "parameter interface pointer first! Furthermore introduce own action types "
      "struct_calc_brownianforce and struct_calc_brownianstiff");
}

/*-----------------------------------------------------------------------------------------------------------*
 | shifts nodes so that proper evaluation is possible even in case of periodic boundary conditions;
 if two   | | nodes within one element are separated by a periodic boundary, one of them is shifted
 such that the final | | distance in R^3 is the same as the initial distance in the periodic space;
 the shift affects computation  | | on element level within that very iteration step, only (no
 change in global variables performed)          |                                 | | (public) cyron
 10/09|
 *----------------------------------------------------------------------------------------------------------*/
template <int nnode, int ndim>
inline void DRT::ELEMENTS::Truss3::NodeShift(
    Teuchos::ParameterList& params, std::vector<double>& disp)
{
  dserror(
      "Truss3::NodeShift is deprecated; if needed adapt parameter handling according to parameter "
      "interface pointer and new sections in input file (statmech section is no longer existent ) "
      "first!");
}
