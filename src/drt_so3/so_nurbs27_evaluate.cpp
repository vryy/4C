/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation methods of the quadratic NURBS 27 element

\level 2

\maintainer Christoph Meier
*----------------------------------------------------------------------*/
#include "so_nurbs27.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_mat/so3_material.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_elements_paramsinterface.H"


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::NURBS::So_nurbs27::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra)
{
  LINALG::Matrix<81, 81> elemat1(elemat1_epetra.A(), true);
  LINALG::Matrix<81, 81> elemat2(elemat2_epetra.A(), true);
  LINALG::Matrix<81, 1> elevec1(elevec1_epetra.A(), true);
  LINALG::Matrix<81, 1> elevec2(elevec2_epetra.A(), true);

  // start with "none"
  DRT::ELEMENTS::NURBS::So_nurbs27::ActionType act = So_nurbs27::none;

  // get the required action
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    dserror("No action supplied");
  else if (action == "calc_struct_linstiff")
    act = So_nurbs27::calc_struct_linstiff;
  else if (action == "calc_struct_nlnstiff")
    act = So_nurbs27::calc_struct_nlnstiff;
  else if (action == "calc_struct_internalforce")
    act = So_nurbs27::calc_struct_internalforce;
  else if (action == "calc_struct_linstiffmass")
    act = So_nurbs27::calc_struct_linstiffmass;
  else if (action == "calc_struct_nlnstiffmass")
    act = So_nurbs27::calc_struct_nlnstiffmass;
  else if (action == "calc_struct_eleload")
    act = So_nurbs27::calc_struct_eleload;
  else if (action == "calc_struct_fsiload")
    act = So_nurbs27::calc_struct_fsiload;
  else if (action == "calc_struct_update_istep")
    act = So_nurbs27::calc_struct_update_istep;
  else if (action == "calc_stc_matrix")
    act = So_nurbs27::calc_stc_matrix;
  else if (action == "calc_stc_matrix_inverse")
    act = So_nurbs27::calc_stc_matrix_inverse;
  else if (action == "calc_struct_reset_istep")
    act = So_nurbs27::calc_struct_reset_istep;
  else if (action == "calc_struct_energy")
    act = So_nurbs27::calc_struct_energy;
  else if (action == "calc_struct_nlnstifflmass")
    act = So_nurbs27::calc_struct_nlnstifflmass;
  else if (action == "calc_struct_recover")
    return 0;
  else if (action == "calc_struct_predict")
    return 0;
  else
    dserror("Unknown type of action '%s' for So_nurbs27", action.c_str());
  // what should the element do
  switch (act)
  {
    // linear stiffness
    case calc_struct_linstiff:
    {
      // need current displacement and residual forces
      std::vector<double> mydisp(lm.size());
      for (unsigned i = 0; i < mydisp.size(); ++i) mydisp[i] = 0.0;
      std::vector<double> myres(lm.size());
      for (unsigned i = 0; i < myres.size(); ++i) myres[i] = 0.0;
      sonurbs27_nlnstiffmass(lm, discretization, mydisp, myres, &elemat1, NULL, &elevec1, params);
    }
    break;

    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (disp == Teuchos::null || res == Teuchos::null)
        dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res, myres, lm);
      LINALG::Matrix<81, 81>* matptr = NULL;
      if (elemat1.IsInitialized()) matptr = &elemat1;

      sonurbs27_nlnstiffmass(lm, discretization, mydisp, myres, matptr, NULL, &elevec1, params);
    }
    break;

    // internal force vector only
    case calc_struct_internalforce:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (disp == Teuchos::null || res == Teuchos::null)
        dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res, myres, lm);
      // create a dummy element matrix to apply linearised EAS-stuff onto
      LINALG::Matrix<81, 81> myemat(true);
      sonurbs27_nlnstiffmass(lm, discretization, mydisp, myres, &myemat, NULL, &elevec1, params);
    }
    break;

    // nonlinear stiffness, internal force vector, and consistent mass matrix
    case calc_struct_nlnstiffmass:
    case calc_struct_nlnstifflmass:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (disp == Teuchos::null || res == Teuchos::null)
        dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res, myres, lm);

      sonurbs27_nlnstiffmass(
          lm, discretization, mydisp, myres, &elemat1, &elemat2, &elevec1, params);

      if (act == calc_struct_nlnstifflmass) lumpmass(&elemat2);
    }
    break;

    case calc_struct_eleload:
      dserror("this method is not supposed to evaluate a load, use EvaluateNeumann(...)");
      break;

    case calc_struct_fsiload:
      dserror("Case not yet implemented");
      break;

    case calc_struct_update_istep:
    {
      // Update of history for materials
      SolidMaterial()->Update();
    }
    break;

    case calc_struct_reset_istep:
    {
      // Reset of history (if needed)
      SolidMaterial()->ResetStep();
    }
    break;

    case calc_stc_matrix_inverse:
    {
      const INPAR::STR::STC_Scale stc_scaling =
          DRT::INPUT::get<INPAR::STR::STC_Scale>(params, "stc_scaling");
      if (stc_scaling == INPAR::STR::stc_none)
        dserror("To scale or not to scale, that's the query!");
      else
      {
        CalcSTCMatrix(elemat1, stc_scaling, params.get<int>("stc_layer"), lm, discretization, true);
      }
    }
    break;

    case calc_stc_matrix:
    {
      const INPAR::STR::STC_Scale stc_scaling =
          DRT::INPUT::get<INPAR::STR::STC_Scale>(params, "stc_scaling");
      if (stc_scaling == INPAR::STR::stc_none)
        dserror("To scale or not to scale, that's the query!");
      else
      {
        CalcSTCMatrix(
            elemat1, stc_scaling, params.get<int>("stc_layer"), lm, discretization, false);
      }
    }
    break;

    case calc_struct_energy:
    {
      if (elevec1_epetra.Length() < 1) dserror("The given result vector is too short.");

      // need current displacement
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

      elevec1_epetra(0) = CalcIntEnergy(discretization, mydisp, params);
      break;
    }

    default:
      dserror("Unknown type of action for So_nurbs27");
  }
  return 0;
}  // DRT::ELEMENTS::So_nurbs27::Evaluate


/*----------------------------------------------------------------------*
 | calc. scaled thickness matrix for thin shell-like structs   (public) |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NURBS::So_nurbs27::CalcSTCMatrix(LINALG::Matrix<81, 81>& elemat1,
    const INPAR::STR::STC_Scale stc_scaling, const int stc_layer, std::vector<int>& lm,
    DRT::Discretization& discretization, bool do_inverse)
{
  // --------------------------------------------------
  // Initialisation of nurbs specific stuff
  std::vector<Epetra_SerialDenseVector> myknots(3);

  // for isogeometric elements:
  //     o get knots
  //     o get weights
  DRT::NURBS::NurbsDiscretization* nurbsdis =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

  if (nurbsdis == NULL)
  {
    dserror("So_nurbs27 appeared in non-nurbs discretisation\n");
  }

  bool zero_ele = (*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots, Id());

  // there is nothing to be done for zero sized elements in knotspan
  if (zero_ele)
  {
    return;
  }

  LINALG::Matrix<27, 1> weights;
  DRT::Node** nodes = Nodes();
  for (int inode = 0; inode < 27; inode++)
  {
    DRT::NURBS::ControlPoint* cp = dynamic_cast<DRT::NURBS::ControlPoint*>(nodes[inode]);

    weights(inode) = cp->W();
  }


  // --------------------------------------------------
  // determine the lengths in r-, s- and t-direction

  // compute coordinates  of corners 0,2,6,18

  LINALG::Matrix<27, 1> funct;

  LINALG::Matrix<3, 1> x0;
  LINALG::Matrix<3, 1> x2;
  LINALG::Matrix<3, 1> x6;
  LINALG::Matrix<3, 1> x18;

  {
    LINALG::Matrix<3, 1> gpa;
    gpa(0) = -1.0;
    gpa(1) = -1.0;
    gpa(2) = -1.0;

    DRT::NURBS::UTILS::nurbs_get_3D_funct(funct, gpa, myknots, weights, DRT::Element::nurbs27);

    for (int isd = 0; isd < 3; ++isd)
    {
      double val = 0;
      for (int inode = 0; inode < 27; ++inode)
      {
        val += (((nodes[inode])->X())[isd]) * funct(inode);
      }
      x0(isd) = val;
    }
  }

  {
    LINALG::Matrix<3, 1> gpa;
    gpa(0) = 1.0;
    gpa(1) = -1.0;
    gpa(2) = -1.0;

    DRT::NURBS::UTILS::nurbs_get_3D_funct(funct, gpa, myknots, weights, DRT::Element::nurbs27);

    for (int isd = 0; isd < 3; ++isd)
    {
      double val = 0;
      for (int inode = 0; inode < 27; ++inode)
      {
        val += (((nodes[inode])->X())[isd]) * funct(inode);
      }
      x2(isd) = val;
    }
  }
  {
    LINALG::Matrix<3, 1> gpa;
    gpa(0) = 1.0;
    gpa(1) = 1.0;
    gpa(2) = -1.0;

    DRT::NURBS::UTILS::nurbs_get_3D_funct(funct, gpa, myknots, weights, DRT::Element::nurbs27);

    for (int isd = 0; isd < 3; ++isd)
    {
      double val = 0;
      for (int inode = 0; inode < 27; ++inode)
      {
        val += (((nodes[inode])->X())[isd]) * funct(inode);
      }
      x6(isd) = val;
    }
  }
  {
    LINALG::Matrix<3, 1> gpa;
    gpa(0) = -1.0;
    gpa(1) = -1.0;
    gpa(2) = 1.0;

    DRT::NURBS::UTILS::nurbs_get_3D_funct(funct, gpa, myknots, weights, DRT::Element::nurbs27);

    for (int isd = 0; isd < 3; ++isd)
    {
      double val = 0;
      for (int inode = 0; inode < 27; ++inode)
      {
        val += (((nodes[inode])->X())[isd]) * funct(inode);
      }
      x18(isd) = val;
    }
  }

  LINALG::Matrix<3, 1> deltaX;

  deltaX.Update(1.0, x2, -1.0, x0);
  const double length_r = deltaX.Norm2();
  deltaX.Update(1.0, x6, -1.0, x0);
  const double length_s = deltaX.Norm2();
  deltaX.Update(1.0, x18, -1.0, x0);
  const double length_t = deltaX.Norm2();

  double ratio = 1.0;

  std::vector<int> topnodeids;
  std::vector<int> midnodeids;
  std::vector<int> botnodeids;

  if (length_t <= length_r && length_t <= length_s)
  {
    for (int i = 0; i < 9; ++i) botnodeids.push_back(i);
    for (int i = 9; i < 18; ++i) midnodeids.push_back(i);
    for (int i = 18; i < 27; ++i) topnodeids.push_back(i);

    ratio = (length_r + length_s) / (2.0 * length_t);
  }
  else if (length_s <= length_r && length_s <= length_t)
  {
    botnodeids.push_back(0);
    botnodeids.push_back(1);
    botnodeids.push_back(2);
    botnodeids.push_back(9);
    botnodeids.push_back(10);
    botnodeids.push_back(11);
    botnodeids.push_back(18);
    botnodeids.push_back(19);
    botnodeids.push_back(20);

    midnodeids.push_back(3);
    midnodeids.push_back(4);
    midnodeids.push_back(5);
    midnodeids.push_back(12);
    midnodeids.push_back(13);
    midnodeids.push_back(14);
    midnodeids.push_back(21);
    midnodeids.push_back(22);
    midnodeids.push_back(23);

    topnodeids.push_back(6);
    topnodeids.push_back(7);
    topnodeids.push_back(8);
    topnodeids.push_back(15);
    topnodeids.push_back(16);
    topnodeids.push_back(17);
    topnodeids.push_back(24);
    topnodeids.push_back(25);
    topnodeids.push_back(26);

    ratio = (length_r + length_t) / (2.0 * length_s);
  }
  else if (length_r <= length_s && length_r <= length_t)
  {
    for (int i = 0; i < 27; i += 3) botnodeids.push_back(i);

    for (int i = 1; i < 27; i += 3) midnodeids.push_back(i);

    for (int i = 2; i < 27; i += 3) topnodeids.push_back(i);

    ratio = (length_t + length_s) / (2.0 * length_r);
  }


  double C = 1.0;
  if (stc_scaling == INPAR::STR::stc_currsym)
  {
    C = ratio;
  }
  else
  {
    C = ratio * ratio;
  }


  double fac1 = 0.0;
  double fac2 = 0.0;

  if (do_inverse)
  {
    fac1 = (1.0 - C);
    fac2 = C;
  }
  else
  {
    fac1 = (C - 1.0) / (C);
    fac2 = 1.0 / C;
  }

  LINALG::Matrix<27, 1> adjele(true);

  for (int i = 0; i < 27; i++)
  {
    adjele(i, 0) = nodes[i]->NumElement();
  }
  /*
    // loop row midnode
    for(int i=0; i<9; i++)
      {
        int dvi=3*midnodeids[i];
        int dui=3*topnodeids[i];
        int dwi=3*botnodeids[i];

        for(int j=0; j<3; j++)
        {
          elemat1(dvi+j,dvi+j)+=fac2/adjele(midnodeids[i],0);
          elemat1(dvi+j,dui+j)+=fac1/adjele(midnodeids[i],0);
          elemat1(dvi+j,dwi+j)+=fac1/adjele(midnodeids[i],0);
        }
      }

    // loop row botnode
    for(int i=0; i<9; i++)
      {
        int dvi=3*botnodeids[i];

        for(int j=0; j<3; j++)
          {
            elemat1(dvi+j,dvi+j)+=1.0/adjele(botnodeids[i],0);
          }
      }

    // loop row topnode
    for(int i=0; i<9; i++)
      {
        int dvi=3*topnodeids[i];

        for(int j=0; j<3; j++)
          {
            elemat1(dvi+j,dvi+j)+=1.0/adjele(topnodeids[i],0);
          }
      }

  */

  // loop row midnode
  for (int i = 0; i < 9; i++)
  {
    int dvi = 3 * midnodeids[i];

    for (int j = 0; j < 3; j++) elemat1(dvi + j, dvi + j) += 1.0 / adjele(midnodeids[i], 0);
  }

  // loop row botnode
  for (int i = 0; i < 9; i++)
  {
    int dvi = 3 * botnodeids[i];
    int dui = 3 * midnodeids[i];

    for (int j = 0; j < 3; j++)
    {
      elemat1(dvi + j, dvi + j) += fac2 * 1.0 / adjele(botnodeids[i], 0);
      elemat1(dvi + j, dui + j) += fac1 * 1.0 / adjele(botnodeids[i], 0);
    }
  }

  // loop row topnode
  for (int i = 0; i < 9; i++)
  {
    int dvi = 3 * topnodeids[i];
    int dui = 3 * midnodeids[i];

    for (int j = 0; j < 3; j++)
    {
      elemat1(dvi + j, dvi + j) += fac2 * 1.0 / adjele(topnodeids[i], 0);
      elemat1(dvi + j, dui + j) += fac1 * 1.0 / adjele(topnodeids[i], 0);
    }
  }

  return;
}  // CalcSTCMatrix



/*----------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)              |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::NURBS::So_nurbs27::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  SetParamsInterfacePtr(params);
  // get values and switches from the condition
  const std::vector<int>* onoff = condition.Get<std::vector<int>>("onoff");
  const std::vector<double>* val = condition.Get<std::vector<double>>("val");

  /*
   **    TIME CURVE BUSINESS
   */
  // find out whether we will use a time curve
  double time = -1.0;
  if (IsParamsInterface())
    time = ParamsInterface().GetTotalTime();
  else
    time = params.get("total time", -1.0);

  // ensure that at least as many curves/functs as dofs are available
  if (int(onoff->size()) < NUMDIM_SONURBS27)
    dserror("Fewer functions or curves defined than the element has dofs.");

  for (int checkdof = NUMDIM_SONURBS27; checkdof < int(onoff->size()); ++checkdof)
  {
    if ((*onoff)[checkdof] != 0)
      dserror("Number of Dimensions in Neumann_Evalutaion is 3. Further DoFs are not considered.");
  }

  // (SPATIAL) FUNCTION BUSINESS
  const std::vector<int>* funct = condition.Get<std::vector<int>>("funct");
  LINALG::Matrix<NUMDIM_SONURBS27, 1> xrefegp(false);
  bool havefunct = false;
  if (funct)
    for (int dim = 0; dim < NUMDIM_SONURBS27; dim++)
      if ((*funct)[dim] > 0) havefunct = havefunct or true;

  // --------------------------------------------------
  // Initialisation of nurbs specific stuff
  std::vector<Epetra_SerialDenseVector> myknots(3);

  // for isogeometric elements:
  //     o get knots
  //     o get weights
  DRT::NURBS::NurbsDiscretization* nurbsdis =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

  if (nurbsdis == NULL) dserror("So_nurbs27 appeared in non-nurbs discretisation\n");

  // there is nothing to be done for zero sized elements in knotspan
  if ((*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots, Id())) return (0);

  LINALG::Matrix<27, 1> weights;
  DRT::Node** nodes = Nodes();
  for (int inode = 0; inode < 27; inode++)
    weights(inode) = dynamic_cast<DRT::NURBS::ControlPoint*>(nodes[inode])->W();

  /*------------------------------------------------------------------*/
  /*                   update element geometry                        */
  /*------------------------------------------------------------------*/

  // material coord. of element
  LINALG::Matrix<27, 3> xrefe;
  for (int i = 0; i < 27; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];
  }
  /* ================================================= Loop over Gauss Points */
  const int numgp = 27;
  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_hex_27point;
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
  LINALG::Matrix<3, 1> gpa;

  LINALG::Matrix<27, 1> shape;
  LINALG::Matrix<3, 27> deriv;

  for (int gp = 0; gp < numgp; ++gp)
  {
    gpa(0) = intpoints.qxg[gp][0];
    gpa(1) = intpoints.qxg[gp][1];
    gpa(2) = intpoints.qxg[gp][2];

    DRT::NURBS::UTILS::nurbs_get_3D_funct_deriv(
        shape, deriv, gpa, myknots, weights, DRT::Element::nurbs27);

    // compute the Jacobian matrix
    LINALG::Matrix<NUMDIM_SONURBS27, NUMDIM_SONURBS27> jac;
    jac.Multiply(deriv, xrefe);

    // compute determinant of Jacobian
    const double detJ = jac.Determinant();
    if (detJ == 0.0)
      dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0)
      dserror("NEGATIVE JACOBIAN DETERMINANT");

    // material/reference co-ordinates of Gauss point
    if (havefunct)
    {
      for (int dim = 0; dim < NUMDIM_SONURBS27; dim++)
      {
        xrefegp(dim) = 0.0;
        for (int nodid = 0; nodid < NUMNOD_SONURBS27; ++nodid)
          xrefegp(dim) += shape(nodid) * xrefe(nodid, dim);
      }
    }

    // integration factor
    const double fac = intpoints.qwgt[gp] * detJ;
    // distribute/add over element load vector
    for (int dim = 0; dim < NUMDIM_SONURBS27; dim++)
    {
      if ((*onoff)[dim])
      {
        // function evaluation
        const int functnum = (funct) ? (*funct)[dim] : -1;
        const double functfac =
            (functnum > 0)
                ? DRT::Problem::Instance()->Funct(functnum - 1).Evaluate(dim, xrefegp.A(), time)
                : 1.0;
        const double dim_fac = (*val)[dim] * fac * functfac;
        for (int nodid = 0; nodid < NUMNOD_SONURBS27; ++nodid)
        {
          elevec1[nodid * NUMDIM_SONURBS27 + dim] += shape(nodid) * dim_fac;
        }
      }
    }

  } /* end of Loop over GP */

  return 0;
}  // DRT::ELEMENTS::So_nurbs27::EvaluateNeumann


/*----------------------------------------------------------------------*
 |  init the element jacobian mapping (protected)                       |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NURBS::So_nurbs27::InitJacobianMapping(DRT::Discretization& dis)
{
  // --------------------------------------------------
  // Initialisation of nurbs specific stuff
  std::vector<Epetra_SerialDenseVector> myknots(3);

  // for isogeometric elements:
  //     o get knots
  //     o get weights
  DRT::NURBS::NurbsDiscretization* nurbsdis =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(dis));

  if (nurbsdis == NULL)
  {
    dserror("So_nurbs27 appeared in non-nurbs discretisation\n");
  }

  bool zero_ele = (*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots, Id());

  // there is nothing to be done for zero sized elements in knotspan
  if (zero_ele)
  {
    return;
  }

  LINALG::Matrix<27, 1> weights;
  DRT::Node** nodes = Nodes();
  for (int inode = 0; inode < 27; inode++)
  {
    DRT::NURBS::ControlPoint* cp = dynamic_cast<DRT::NURBS::ControlPoint*>(nodes[inode]);

    weights(inode) = cp->W();
  }

  const static std::vector<LINALG::Matrix<3, 27>> derivs = sonurbs27_derivs(myknots, weights);
  LINALG::Matrix<27, 3> xrefe;
  for (int i = 0; i < 27; ++i)
  {
    xrefe(i, 0) = Nodes()[i]->X()[0];
    xrefe(i, 1) = Nodes()[i]->X()[1];
    xrefe(i, 2) = Nodes()[i]->X()[2];
  }

  const int numgp = 27;

  invJ_.resize(numgp);
  detJ_.resize(numgp);
  for (int gp = 0; gp < numgp; ++gp)
  {
    invJ_[gp].Multiply(derivs[gp], xrefe);
    detJ_[gp] = invJ_[gp].Invert();
    if (detJ_[gp] == 0.0)
      dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ_[gp] < 0.0)
      dserror("NEGATIVE JACOBIAN DETERMINANT %12.5e IN ELEMENT ID %d, gauss point %d", detJ_[gp],
          Id(), gp);
  }
  return;
}  // DRT::ELEMENTS::So_nurbs27::InitJacobianMapping()

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                                      |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NURBS::So_nurbs27::sonurbs27_nlnstiffmass(
    std::vector<int>& lm,                 // location matrix
    DRT::Discretization& discretization,  // discretisation to extract knot vector
    std::vector<double>& disp,            // current displacements
    std::vector<double>& residual,        // current residual displ
    LINALG::Matrix<81, 81>* stiffmatrix,  // element stiffness matrix
    LINALG::Matrix<81, 81>* massmatrix,   // element mass matrix
    LINALG::Matrix<81, 1>* force,         // element internal force vector
    Teuchos::ParameterList& params)       // strain output option
{
  // --------------------------------------------------
  // Initialisation of nurbs specific stuff
  std::vector<Epetra_SerialDenseVector> myknots(3);

  // for isogeometric elements:
  //     o get knots
  //     o get weights
  DRT::NURBS::NurbsDiscretization* nurbsdis =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

  if (nurbsdis == NULL)
  {
    dserror("So_nurbs27 appeared in non-nurbs discretisation\n");
  }

  bool zero_ele = (*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots, Id());

  // there is nothing to be done for zero sized elements in knotspan
  if (zero_ele)
  {
    return;
  }

  LINALG::Matrix<27, 1> weights;
  DRT::Node** nodes = Nodes();
  for (int inode = 0; inode < 27; inode++)
  {
    DRT::NURBS::ControlPoint* cp = dynamic_cast<DRT::NURBS::ControlPoint*>(nodes[inode]);

    weights(inode) = cp->W();
  }

  // update element geometry
  LINALG::Matrix<27, 3> xrefe;  // material coord. of element
  LINALG::Matrix<27, 3> xcurr;  // current  coord. of element
  for (int i = 0; i < 27; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];

    xcurr(i, 0) = xrefe(i, 0) + disp[i * 3];
    xcurr(i, 1) = xrefe(i, 1) + disp[i * 3 + 1];
    xcurr(i, 2) = xrefe(i, 2) + disp[i * 3 + 2];
  }

  /*------------------------------------------------------------------*/
  /*                    Loop over Gauss Points                        */
  /*------------------------------------------------------------------*/
  const int numgp = 27;

  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_hex_27point;
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);

  invJ_.resize(numgp);
  detJ_.resize(numgp);

  LINALG::Matrix<27, 1> funct;
  LINALG::Matrix<3, 27> deriv;

  LINALG::Matrix<3, 27> N_XYZ;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  LINALG::Matrix<3, 3> defgrd(true);
  for (int gp = 0; gp < numgp; ++gp)
  {
    LINALG::Matrix<3, 1> gpa;
    gpa(0) = intpoints.qxg[gp][0];
    gpa(1) = intpoints.qxg[gp][1];
    gpa(2) = intpoints.qxg[gp][2];

    DRT::NURBS::UTILS::nurbs_get_3D_funct_deriv(
        funct, deriv, gpa, myknots, weights, DRT::Element::nurbs27);

    /* get the inverse of the Jacobian matrix which looks like:
    **            [ x_,r  y_,r  z_,r ]^-1
    **     J^-1 = [ x_,s  y_,s  z_,s ]
    **            [ x_,t  y_,t  z_,t ]
    */
    LINALG::Matrix<3, 3> invJac(true);

    invJac.Multiply(deriv, xrefe);
    double detJ = invJac.Invert();

    if (detJ == 0.0)
      dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0)
      dserror("NEGATIVE JACOBIAN DETERMINANT %12.5e IN ELEMENT ID %d, gauss point %d", detJ_[gp],
          Id(), gp);

    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.Multiply(invJac, deriv);

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
    defgrd.MultiplyTT(xcurr, N_XYZ);

    // Right Cauchy-Green tensor = F^T * F
    LINALG::Matrix<3, 3> cauchygreen;
    cauchygreen.MultiplyTN(defgrd, defgrd);

    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    Epetra_SerialDenseVector glstrain_epetra(6);
    LINALG::Matrix<6, 1> glstrain(glstrain_epetra.A(), true);
    glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
    glstrain(3) = cauchygreen(0, 1);
    glstrain(4) = cauchygreen(1, 2);
    glstrain(5) = cauchygreen(2, 0);

    /* non-linear B-operator (may so be called, meaning
    ** of B-operator is not so sharp in the non-linear realm) *
    ** B = F . Bl *
    **
    **      [ ... | F_11*N_{,1}^k  F_21*N_{,1}^k  F_31*N_{,1}^k | ... ]
    **      [ ... | F_12*N_{,2}^k  F_22*N_{,2}^k  F_32*N_{,2}^k | ... ]
    **      [ ... | F_13*N_{,3}^k  F_23*N_{,3}^k  F_33*N_{,3}^k | ... ]
    ** B =  [ ~~~   ~~~~~~~~~~~~~  ~~~~~~~~~~~~~  ~~~~~~~~~~~~~   ~~~ ]
    **      [       F_11*N_{,2}^k+F_12*N_{,1}^k                       ]
    **      [ ... |          F_21*N_{,2}^k+F_22*N_{,1}^k        | ... ]
    **      [                       F_31*N_{,2}^k+F_32*N_{,1}^k       ]
    **      [                                                         ]
    **      [       F_12*N_{,3}^k+F_13*N_{,2}^k                       ]
    **      [ ... |          F_22*N_{,3}^k+F_23*N_{,2}^k        | ... ]
    **      [                       F_32*N_{,3}^k+F_33*N_{,2}^k       ]
    **      [                                                         ]
    **      [       F_13*N_{,1}^k+F_11*N_{,3}^k                       ]
    **      [ ... |          F_23*N_{,1}^k+F_21*N_{,3}^k        | ... ]
    **      [                       F_33*N_{,1}^k+F_31*N_{,3}^k       ]
    */
    LINALG::Matrix<6, 81> bop;
    for (int i = 0; i < 27; ++i)
    {
      bop(0, 3 * i) = defgrd(0, 0) * N_XYZ(0, i);
      bop(0, 3 * i + 1) = defgrd(1, 0) * N_XYZ(0, i);
      bop(0, 3 * i + 2) = defgrd(2, 0) * N_XYZ(0, i);
      bop(1, 3 * i) = defgrd(0, 1) * N_XYZ(1, i);
      bop(1, 3 * i + 1) = defgrd(1, 1) * N_XYZ(1, i);
      bop(1, 3 * i + 2) = defgrd(2, 1) * N_XYZ(1, i);
      bop(2, 3 * i) = defgrd(0, 2) * N_XYZ(2, i);
      bop(2, 3 * i + 1) = defgrd(1, 2) * N_XYZ(2, i);
      bop(2, 3 * i + 2) = defgrd(2, 2) * N_XYZ(2, i);
      /* ~~~ */
      bop(3, 3 * i) = defgrd(0, 0) * N_XYZ(1, i) + defgrd(0, 1) * N_XYZ(0, i);
      bop(3, 3 * i + 1) = defgrd(1, 0) * N_XYZ(1, i) + defgrd(1, 1) * N_XYZ(0, i);
      bop(3, 3 * i + 2) = defgrd(2, 0) * N_XYZ(1, i) + defgrd(2, 1) * N_XYZ(0, i);
      bop(4, 3 * i) = defgrd(0, 1) * N_XYZ(2, i) + defgrd(0, 2) * N_XYZ(1, i);
      bop(4, 3 * i + 1) = defgrd(1, 1) * N_XYZ(2, i) + defgrd(1, 2) * N_XYZ(1, i);
      bop(4, 3 * i + 2) = defgrd(2, 1) * N_XYZ(2, i) + defgrd(2, 2) * N_XYZ(1, i);
      bop(5, 3 * i) = defgrd(0, 2) * N_XYZ(0, i) + defgrd(0, 0) * N_XYZ(2, i);
      bop(5, 3 * i + 1) = defgrd(1, 2) * N_XYZ(0, i) + defgrd(1, 0) * N_XYZ(2, i);
      bop(5, 3 * i + 2) = defgrd(2, 2) * N_XYZ(0, i) + defgrd(2, 0) * N_XYZ(2, i);
    }

    // call material law
    LINALG::Matrix<6, 6> cmat(true);
    LINALG::Matrix<6, 1> stress(true);
    params.set<int>("gp", gp);
    SolidMaterial()->Evaluate(&defgrd, &glstrain, params, &stress, &cmat, Id());
    // end of call material law

    double detJ_w = detJ * intpoints.qwgt[gp];
    // update internal force vector
    if (force != NULL)
    {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      force->MultiplyTN(detJ_w, bop, stress, 1.0);
    }

    // update stiffness matrix
    if (stiffmatrix != NULL)
    {
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      LINALG::Matrix<6, 81> cb;
      cb.Multiply(cmat, bop);
      stiffmatrix->MultiplyTN(detJ_w, bop, cb, 1.0);

      // integrate `geometric' stiffness matrix and add to keu *****************
      LINALG::Matrix<6, 1> sfac(stress);  // auxiliary integrated stress
      sfac.Scale(detJ_w);                 // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
      std::vector<double> SmB_L(3);       // intermediate Sm.B_L
      // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
      for (int inod = 0; inod < 27; ++inod)
      {
        SmB_L[0] = sfac(0) * N_XYZ(0, inod) + sfac(3) * N_XYZ(1, inod) + sfac(5) * N_XYZ(2, inod);
        SmB_L[1] = sfac(3) * N_XYZ(0, inod) + sfac(1) * N_XYZ(1, inod) + sfac(4) * N_XYZ(2, inod);
        SmB_L[2] = sfac(5) * N_XYZ(0, inod) + sfac(4) * N_XYZ(1, inod) + sfac(2) * N_XYZ(2, inod);
        for (int jnod = 0; jnod < 27; ++jnod)
        {
          double bopstrbop = 0.0;  // intermediate value
          for (int idim = 0; idim < 3; ++idim)
          {
            bopstrbop += N_XYZ(idim, jnod) * SmB_L[idim];
          }

          (*stiffmatrix)(3 * inod, 3 * jnod) += bopstrbop;
          (*stiffmatrix)(3 * inod + 1, 3 * jnod + 1) += bopstrbop;
          (*stiffmatrix)(3 * inod + 2, 3 * jnod + 2) += bopstrbop;
        }
      }  // end of integrate `geometric' stiffness
    }    // if (stiffmatrix)

    if (massmatrix != NULL)  // evaluate mass matrix
    {
      double density = Material()->Density(gp);
      // integrate consistent mass matrix
      const double factor = detJ_w * density;
      double ifactor, massfactor;
      for (int inod = 0; inod < 27; ++inod)
      {
        ifactor = funct(inod) * factor;
        for (int jnod = 0; jnod < 27; ++jnod)
        {
          massfactor = funct(jnod) * ifactor;  // intermediate factor
          (*massmatrix)(3 * inod, 3 * jnod) += massfactor;
          (*massmatrix)(3 * inod + 1, 3 * jnod + 1) += massfactor;
          (*massmatrix)(3 * inod + 2, 3 * jnod + 2) += massfactor;
        }
      }
    }  // end of mass matrix

  } /* end of Loop over GP */

  return;
}  // DRT::ELEMENTS::So_nurbs27::sonurbs27_nlnstiffmass

/*----------------------------------------------------------------------*
 |  Evaluate nurbs27 Shape fcts at all 27 Gauss Points                     |
 *----------------------------------------------------------------------*/
const std::vector<LINALG::Matrix<27, 1>> DRT::ELEMENTS::NURBS::So_nurbs27::sonurbs27_shapefcts(
    const std::vector<Epetra_SerialDenseVector>& myknots, const LINALG::Matrix<27, 1>& weights)
{
  const int numgp = 27;

  std::vector<LINALG::Matrix<27, 1>> shapefcts(numgp);
  // (r,s,t) gp-locations of fully integrated quadratic Nurbs 27
  // fill up nodal f at each gp
  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_hex_27point;
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
  for (int igp = 0; igp < intpoints.nquad; ++igp)
  {
    LINALG::Matrix<3, 1> gp;
    gp(0) = intpoints.qxg[igp][0];
    gp(1) = intpoints.qxg[igp][1];
    gp(2) = intpoints.qxg[igp][2];

    DRT::NURBS::UTILS::nurbs_get_3D_funct(
        shapefcts[igp], gp, myknots, weights, DRT::Element::nurbs27);
  }
  return shapefcts;
}


/*----------------------------------------------------------------------*
 |  Evaluate nurbs27 Shape fct derivs at all 27 Gauss Points              |
 *----------------------------------------------------------------------*/
const std::vector<LINALG::Matrix<3, 27>> DRT::ELEMENTS::NURBS::So_nurbs27::sonurbs27_derivs(
    const std::vector<Epetra_SerialDenseVector>& myknots, const LINALG::Matrix<27, 1>& weights)
{
  const int numgp = 27;

  std::vector<LINALG::Matrix<3, 27>> derivs(numgp);
  // (r,s,t) gp-locations of fully integrated quadratic Nurbs 27
  // fill up df w.r.t. rst directions (NUMDIM) at each gp
  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_hex_27point;
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
  for (int igp = 0; igp < intpoints.nquad; ++igp)
  {
    LINALG::Matrix<3, 1> gp;
    gp(0) = intpoints.qxg[igp][0];
    gp(1) = intpoints.qxg[igp][1];
    gp(2) = intpoints.qxg[igp][2];

    LINALG::Matrix<27, 1> dummyfct;

    DRT::NURBS::UTILS::nurbs_get_3D_funct_deriv(
        dummyfct, derivs[igp], gp, myknots, weights, DRT::Element::nurbs27);
  }
  return derivs;
}

/*----------------------------------------------------------------------*
 |  Evaluate nurbs27 Weights at all 27 Gauss Points                     |
 *----------------------------------------------------------------------*/
const std::vector<double> DRT::ELEMENTS::NURBS::So_nurbs27::sonurbs27_gpweights()
{
  const int numgp = 27;

  std::vector<double> gpweights(numgp);
  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_hex_27point;
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
  for (int i = 0; i < numgp; ++i)
  {
    gpweights[i] = intpoints.qwgt[i];
  }
  return gpweights;
}


/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::NURBS::So_nurbs27Type::Initialize(DRT::Discretization& dis)
{
  for (int i = 0; i < dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::NURBS::So_nurbs27* actele =
        dynamic_cast<DRT::ELEMENTS::NURBS::So_nurbs27*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_nurbs27* failed");
    actele->InitJacobianMapping(dis);
  }
  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate internal energy of the element (private)                  |
 *----------------------------------------------------------------------*/
double DRT::ELEMENTS::NURBS::So_nurbs27::CalcIntEnergy(
    DRT::Discretization& discretization,  // discretisation to extract knot vector
    std::vector<double>& disp,            // current displacements
    Teuchos::ParameterList& params)       // strain output option
{
  double energy = 0.;

  // --------------------------------------------------
  // Initialisation of nurbs specific stuff
  std::vector<Epetra_SerialDenseVector> myknots(3);

  // for isogeometric elements:
  //     o get knots
  //     o get weights
  DRT::NURBS::NurbsDiscretization* nurbsdis =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));
  if (nurbsdis == NULL) dserror("So_nurbs27 appeared in non-nurbs discretisation\n");

  bool zero_ele = (*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots, Id());

  // there is nothing to be done for zero sized elements in knotspan
  if (zero_ele) return 0.;

  LINALG::Matrix<27, 1> weights;
  DRT::Node** nodes = Nodes();
  for (int inode = 0; inode < 27; inode++)
  {
    DRT::NURBS::ControlPoint* cp = dynamic_cast<DRT::NURBS::ControlPoint*>(nodes[inode]);

    weights(inode) = cp->W();
  }

  // update element geometry
  LINALG::Matrix<27, 3> xrefe;  // material coord. of element
  LINALG::Matrix<27, 3> xcurr;  // current  coord. of element
  for (int i = 0; i < 27; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];

    xcurr(i, 0) = xrefe(i, 0) + disp[i * 3];
    xcurr(i, 1) = xrefe(i, 1) + disp[i * 3 + 1];
    xcurr(i, 2) = xrefe(i, 2) + disp[i * 3 + 2];
  }
  /*------------------------------------------------------------------*/
  /*                    Loop over Gauss Points                        */
  /*------------------------------------------------------------------*/
  const int numgp = 27;

  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_hex_27point;
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);

  invJ_.resize(numgp);
  detJ_.resize(numgp);

  LINALG::Matrix<27, 1> funct;
  LINALG::Matrix<3, 27> deriv;

  LINALG::Matrix<3, 27> N_XYZ;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  LINALG::Matrix<3, 3> defgrd(true);
  for (int gp = 0; gp < numgp; ++gp)
  {
    LINALG::Matrix<3, 1> gpa;
    gpa(0) = intpoints.qxg[gp][0];
    gpa(1) = intpoints.qxg[gp][1];
    gpa(2) = intpoints.qxg[gp][2];

    DRT::NURBS::UTILS::nurbs_get_3D_funct_deriv(
        funct, deriv, gpa, myknots, weights, DRT::Element::nurbs27);

    /* get the inverse of the Jacobian matrix which looks like:
    **            [ x_,r  y_,r  z_,r ]^-1
    **     J^-1 = [ x_,s  y_,s  z_,s ]
    **            [ x_,t  y_,t  z_,t ]
    */
    LINALG::Matrix<3, 3> invJac(true);

    invJac.Multiply(deriv, xrefe);
    double detJ = invJac.Invert();

    if (detJ == 0.0)
      dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0)
      dserror("NEGATIVE JACOBIAN DETERMINANT %12.5e IN ELEMENT ID %d, gauss point %d", detJ_[gp],
          Id(), gp);

    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.Multiply(invJac, deriv);

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
    defgrd.MultiplyTT(xcurr, N_XYZ);

    // Right Cauchy-Green tensor = F^T * F
    LINALG::Matrix<3, 3> cauchygreen;
    cauchygreen.MultiplyTN(defgrd, defgrd);

    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    Epetra_SerialDenseVector glstrain_epetra(6);
    LINALG::Matrix<6, 1> glstrain(glstrain_epetra.A(), true);
    glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
    glstrain(3) = cauchygreen(0, 1);
    glstrain(4) = cauchygreen(1, 2);
    glstrain(5) = cauchygreen(2, 0);

    double psi = 0.0;
    SolidMaterial()->StrainEnergy(glstrain, psi, Id());

    double detJ_w = detJ * intpoints.qwgt[gp];
    energy += detJ_w * psi;
  }

  return energy;
}

/*----------------------------------------------------------------------*
 |  lump mass matrix (private)                               bborn 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NURBS::So_nurbs27::lumpmass(
    LINALG::Matrix<NUMDOF_SONURBS27, NUMDOF_SONURBS27>* emass)
{
  // lump mass matrix
  if (emass != NULL)
  {
    // we assume #elemat2 is a square matrix
    for (unsigned int c = 0; c < (*emass).N(); ++c)  // parse columns
    {
      double d = 0.0;
      for (unsigned int r = 0; r < (*emass).M(); ++r)  // parse rows
      {
        d += (*emass)(r, c);  // accumulate row entries
        (*emass)(r, c) = 0.0;
      }
      (*emass)(c, c) = d;  // apply sum of row entries on diagonal
    }
  }
}
