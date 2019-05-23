/*!----------------------------------------------------------------------

\brief Solid-scatra elements evaluate

\level 2

   \maintainer Christoph Schmidt

*----------------------------------------------------------------------*/

#include "so3_scatra.H"

#include "../drt_lib/drt_utils.H"
#include "../drt_mat/so3_material.H"

/*----------------------------------------------------------------------*
 |  preevaluate the element (public)                                       |
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Scatra<so3_ele, distype>::PreEvaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la)
{
  if (la.Size() > 1)
  {
    // ask for the number of dofs of second dofset (scatra)
    const int numscal = discretization.NumDof(1, Nodes()[0]);

    if (la[1].Size() != numnod_ * numscal)
      dserror("So3_Scatra: PreEvaluate: Location vector length for concentrations does not match!");

    if (discretization.HasState(1, "scalarfield"))  // if concentrations were set
    {
      if (not(distype == DRT::Element::hex8 or distype == DRT::Element::hex27 or
              distype == DRT::Element::tet4 or distype == DRT::Element::tet10))
        dserror(
            "The Solidscatra elements are only tested for the Hex8, Hex27 and Tet4 case. The "
            "following should work, but keep your eyes open (especially with the order of the Gauß "
            "points");

      /* =========================================================================*/
      // start concentration business
      /* =========================================================================*/
      Teuchos::RCP<std::vector<std::vector<double>>> gpconc = Teuchos::rcp(
          new std::vector<std::vector<double>>(numgpt_, std::vector<double>(numscal, 0.0)));

      // check if you can get the scalar state
      Teuchos::RCP<const Epetra_Vector> concnp = discretization.GetState(1, "scalarfield");

      if (concnp == Teuchos::null)
        dserror("calc_struct_nlnstiff: Cannot get state vector 'temperature' ");

      // extract local values of the global vectors
      Teuchos::RCP<std::vector<double>> myconc =
          Teuchos::rcp(new std::vector<double>(la[1].lm_.size(), 0.0));

      DRT::UTILS::ExtractMyValues(*concnp, *myconc, la[1].lm_);

      // element vector for k-th scalar
      std::vector<LINALG::Matrix<numnod_, 1>> econc(numscal);
      for (int k = 0; k < numscal; ++k)
      {
        for (int i = 0; i < numnod_; ++i)
        {
          (econc.at(k))(i, 0) = myconc->at(numscal * i + k);
        }
      }

      /* =========================================================================*/
      /* ================================================= Loop over Gauss Points */
      /* =========================================================================*/
      // volume of current element in reference configuration
      double volume_ref = 0.0;
      // mass in current element in reference configuration
      std::vector<double> mass_ref(numscal, 0.0);

      for (int igp = 0; igp < numgpt_; ++igp)
      {
        // detJrefpar_wgp = det(dX/dr) * w_gp to calculate volume in reference configuration
        const double detJrefpar_wgp = detJ_[igp] * intpoints_.qwgt[igp];

        volume_ref += detJrefpar_wgp;

        // concentrations at current gauß point
        std::vector<double> conc_gp_k(numscal, 0.0);

        // shape functions evaluated at current gauß point
        LINALG::Matrix<numnod_, 1> shapefunct_gp(true);
        DRT::UTILS::shape_function<distype>(xsi_[igp], shapefunct_gp);

        for (int k = 0; k < numscal; ++k)
        {
          // identical shapefunctions for displacements and temperatures
          conc_gp_k.at(k) = shapefunct_gp.Dot(econc.at(k));

          mass_ref.at(k) += conc_gp_k.at(k) * detJrefpar_wgp;
        }

        gpconc->at(igp) = conc_gp_k;
      }

      params.set<Teuchos::RCP<std::vector<std::vector<double>>>>("gp_conc", gpconc);

      // compute average concentrations
      for (int k = 0; k < numscal; ++k)
      {
        // now mass_ref is the element averaged concentration
        mass_ref.at(k) /= volume_ref;
      }
      Teuchos::RCP<std::vector<std::vector<double>>> avgconc =
          Teuchos::rcp(new std::vector<std::vector<double>>(numgpt_, mass_ref));

      params.set<Teuchos::RCP<std::vector<std::vector<double>>>>("avg_conc", avgconc);

    }  // if (discretization.HasState(1,"scalarfield"))

    // If you need a pointer to the scatra material, use these lines:
    // we assume that the second material of the structure is the scatra element material
    // Teuchos::RCP<MAT::Material> scatramat = so3_ele::Material(1);
    // params.set< Teuchos::RCP<MAT::Material> >("scatramat",scatramat);
  }

  // TODO: (thon) actually we do not want this here, since it has nothing to do with scatra specific
  // stuff. But for now we let it be...
  std::vector<double> center = DRT::UTILS::ElementCenterRefeCoords(this);
  Teuchos::RCP<std::vector<double>> xrefe = Teuchos::rcp(new std::vector<double>(center));
  params.set<Teuchos::RCP<std::vector<double>>>("position", xrefe);

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                                       |
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::So3_Scatra<so3_ele, distype>::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra)
{
  // start with ActionType "none"
  typename So3_Scatra::ActionType act = So3_Scatra::none;

  // get the required action
  std::string action = params.get<std::string>("action", "none");

  // get the required action and safety check
  if (action == "none")
    dserror("No action supplied");
  else if (action == "calc_struct_stiffscalar")
    act = So3_Scatra::calc_struct_stiffscalar;

  // at the moment all cases need the PreEvaluate routine, since we always need the concentration
  // value at the gp
  PreEvaluate(params, discretization, la);

  // what action shall be performed
  switch (act)
  {
    // coupling terms K_dS of stiffness matrix K^{SSI} for monolithic SSI
    case So3_Scatra::calc_struct_stiffscalar:
    {
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState(0, "displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");

      // get my displacement vector
      std::vector<double> mydisp((la[0].lm_).size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, la[0].lm_);

      // calculate the stiffness matrix
      nln_kdS_ssi(la, mydisp, elemat1_epetra, params);

      break;
    }

    default:
    {
      // call the base class routine
      so3_ele::Evaluate(params, discretization, la[0].lm_, elemat1_epetra, elemat2_epetra,
          elevec1_epetra, elevec2_epetra, elevec3_epetra);
      break;
    }  // default
  }    // switch(act)

  return 0;
}  // Evaluate


/*----------------------------------------------------------------------*
 | evaluate only the mechanical-scatra stiffness term     schmidt 10/17 |
 | for monolithic SSI, contribution to k_dS (private)                   |
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Scatra<so3_ele, distype>::nln_kdS_ssi(DRT::Element::LocationArray& la,
    std::vector<double>& disp,                  // current displacement
    Epetra_SerialDenseMatrix& stiffmatrix_kdS,  // (numdim_*numnod_ ; numnod_)
    Teuchos::ParameterList& params)
{
  // calculate current and material coordinates of element
  LINALG::Matrix<numnod_, numdim_> xrefe(true);  // X, material coord. of element
  LINALG::Matrix<numnod_, numdim_> xcurr(true);  // x, current  coord. of element
  DRT::Node** nodes = Nodes();
  for (int i = 0; i < numnod_; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i, 0) = x[0];
    xrefe(i, 1) = x[1];
    xrefe(i, 2) = x[2];

    xcurr(i, 0) = xrefe(i, 0) + disp[i * numdofpernode_ + 0];
    xcurr(i, 1) = xrefe(i, 1) + disp[i * numdofpernode_ + 1];
    xcurr(i, 2) = xrefe(i, 2) + disp[i * numdofpernode_ + 2];
  }

  // shape functions and their first derivatives
  LINALG::Matrix<numnod_, 1> shapefunct(true);
  LINALG::Matrix<numdim_, numnod_> deriv(true);
  // compute derivatives N_XYZ at gp w.r.t. material coordinates
  LINALG::Matrix<numdim_, numnod_> N_XYZ(true);
  // compute deformation gradient w.r.t. to material configuration
  LINALG::Matrix<numdim_, numdim_> defgrad(true);

  // get numscatradofspernode from parameter list
  const int numscatradofspernode = params.get<int>("numscatradofspernode", -1);
  if (numscatradofspernode == -1)
    dserror("Could not read 'numscatradofspernode' from parameter list!");

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  for (int gp = 0; gp < numgpt_; ++gp)
  {
    // get shape functions and their derivatives
    DRT::UTILS::shape_function<distype>(xsi_[gp], shapefunct);
    DRT::UTILS::shape_function_deriv1<distype>(xsi_[gp], deriv);

    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 . N_rst
    N_XYZ.Multiply(invJ_[gp], deriv);

    // (material) deformation gradient
    // F = d xcurr / d xrefe = xcurr^T . N_XYZ^T
    defgrad.MultiplyTT(xcurr, N_XYZ);

    // right Cauchy-Green tensor = F^T . F
    LINALG::Matrix<3, 3> cauchygreen;
    cauchygreen.MultiplyTN(defgrad, defgrad);

    // calculate vector of right Cauchy-Green tensor
    LINALG::Matrix<numstr_, 1> cauchygreenvec;
    cauchygreenvec(0) = cauchygreen(0, 0);
    cauchygreenvec(1) = cauchygreen(1, 1);
    cauchygreenvec(2) = cauchygreen(2, 2);
    cauchygreenvec(3) = 2 * cauchygreen(0, 1);
    cauchygreenvec(4) = 2 * cauchygreen(1, 2);
    cauchygreenvec(5) = 2 * cauchygreen(2, 0);

    // Green Lagrange strain
    LINALG::Matrix<numstr_, 1> glstrain;
    // Green-Lagrange strain matrix E = 0.5 * (Cauchygreen - Identity)
    glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
    glstrain(3) = cauchygreen(0, 1);
    glstrain(4) = cauchygreen(1, 2);
    glstrain(5) = cauchygreen(2, 0);

    // calculate nonlinear B-operator
    LINALG::Matrix<numstr_, numdofperelement_> bop(true);
    CalculateBop(&bop, &defgrad, &N_XYZ);

    /*==== call material law ======================================================*/
    // init derivative of second Piola-Kirchhoff stresses w.r.t. concentrations dSdc
    LINALG::Matrix<numstr_, 1> dSdc(true);
    // set current gauss point
    params.set<int>("gp", gp);

    // get dSdc, hand in NULL as 'cmat' to evaluate the off-diagonal block
    Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_static_cast<MAT::So3Material>(Material());
    so3mat->Evaluate(&defgrad, &glstrain, params, &dSdc, NULL, Id());

    /*==== end of call material law ===============================================*/

    // k_dS = B^T . dS/dc * detJ * N * w(gp)
    const double detJ_w = detJ_[gp] * intpoints_.qwgt[gp];
    LINALG::Matrix<numdofperelement_, 1> BdSdc(true);
    BdSdc.MultiplyTN(detJ_w, bop, dSdc);

    // loop over rows
    for (unsigned rowi = 0; rowi < numdofperelement_; ++rowi)
    {
      const double BdSdc_rowi = BdSdc(rowi, 0);
      // loop over columns
      for (unsigned coli = 0; coli < numnod_; ++coli)
      {
        stiffmatrix_kdS(rowi, coli * numscatradofspernode) += BdSdc_rowi * shapefunct(coli, 0);
      }
    }
  }  // gauss point loop

  return;
}  // nln_kdS_ssi


/*----------------------------------------------------------------------*
 | calculate the nonlinear B-operator (private)           schmidt 10/17 |
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Scatra<so3_ele, distype>::CalculateBop(
    LINALG::Matrix<numstr_, numdofperelement_>* bop,  //!< (o): nonlinear B-operator
    const LINALG::Matrix<numdim_, numdim_>* defgrad,  //!< (i): deformation gradient
    const LINALG::Matrix<numdim_, numnod_>* N_XYZ)
    const  //!< (i): (material) derivative of shape functions
{
  // calc bop matrix if provided
  if (bop != NULL)
  {
    /* non-linear B-operator (may so be called, meaning of B-operator is not so
    **  sharp in the non-linear realm) *
    **   B = F^{i,T} . B_L *
    ** with linear B-operator B_L =  N_XYZ (6x24) = (3x8)
    **
    **   B    =   F^T  . N_XYZ
    ** (6x24)    (3x3)   (3x8)
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
    for (int i = 0; i < numnod_; ++i)
    {
      (*bop)(0, numdofpernode_ * i + 0) = (*defgrad)(0, 0) * (*N_XYZ)(0, i);
      (*bop)(0, numdofpernode_ * i + 1) = (*defgrad)(1, 0) * (*N_XYZ)(0, i);
      (*bop)(0, numdofpernode_ * i + 2) = (*defgrad)(2, 0) * (*N_XYZ)(0, i);
      (*bop)(1, numdofpernode_ * i + 0) = (*defgrad)(0, 1) * (*N_XYZ)(1, i);
      (*bop)(1, numdofpernode_ * i + 1) = (*defgrad)(1, 1) * (*N_XYZ)(1, i);
      (*bop)(1, numdofpernode_ * i + 2) = (*defgrad)(2, 1) * (*N_XYZ)(1, i);
      (*bop)(2, numdofpernode_ * i + 0) = (*defgrad)(0, 2) * (*N_XYZ)(2, i);
      (*bop)(2, numdofpernode_ * i + 1) = (*defgrad)(1, 2) * (*N_XYZ)(2, i);
      (*bop)(2, numdofpernode_ * i + 2) = (*defgrad)(2, 2) * (*N_XYZ)(2, i);
      /* ~~~ */
      (*bop)(3, numdofpernode_ * i + 0) =
          (*defgrad)(0, 0) * (*N_XYZ)(1, i) + (*defgrad)(0, 1) * (*N_XYZ)(0, i);
      (*bop)(3, numdofpernode_ * i + 1) =
          (*defgrad)(1, 0) * (*N_XYZ)(1, i) + (*defgrad)(1, 1) * (*N_XYZ)(0, i);
      (*bop)(3, numdofpernode_ * i + 2) =
          (*defgrad)(2, 0) * (*N_XYZ)(1, i) + (*defgrad)(2, 1) * (*N_XYZ)(0, i);
      (*bop)(4, numdofpernode_ * i + 0) =
          (*defgrad)(0, 1) * (*N_XYZ)(2, i) + (*defgrad)(0, 2) * (*N_XYZ)(1, i);
      (*bop)(4, numdofpernode_ * i + 1) =
          (*defgrad)(1, 1) * (*N_XYZ)(2, i) + (*defgrad)(1, 2) * (*N_XYZ)(1, i);
      (*bop)(4, numdofpernode_ * i + 2) =
          (*defgrad)(2, 1) * (*N_XYZ)(2, i) + (*defgrad)(2, 2) * (*N_XYZ)(1, i);
      (*bop)(5, numdofpernode_ * i + 0) =
          (*defgrad)(0, 2) * (*N_XYZ)(0, i) + (*defgrad)(0, 0) * (*N_XYZ)(2, i);
      (*bop)(5, numdofpernode_ * i + 1) =
          (*defgrad)(1, 2) * (*N_XYZ)(0, i) + (*defgrad)(1, 0) * (*N_XYZ)(2, i);
      (*bop)(5, numdofpernode_ * i + 2) =
          (*defgrad)(2, 2) * (*N_XYZ)(0, i) + (*defgrad)(2, 0) * (*N_XYZ)(2, i);
    }
  }

  return;
}  // CalculateBop


/*----------------------------------------------------------------------*
 | initialize element (private)                            schmidt 10/17|
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Scatra<so3_ele, distype>::InitElement()
{
  // resize gauss point coordinates, inverse of the jacobian and determinant of the jacobian
  xsi_.resize(numgpt_);
  invJ_.resize(numgpt_);
  detJ_.resize(numgpt_);

  // calculate coordinates in reference (material) configuration
  LINALG::Matrix<numnod_, numdim_> xrefe;
  for (int i = 0; i < numnod_; ++i)
  {
    Node** nodes = Nodes();
    if (!nodes) dserror("Nodes() returned null pointer");
    xrefe(i, 0) = Nodes()[i]->X()[0];
    xrefe(i, 1) = Nodes()[i]->X()[1];
    xrefe(i, 2) = Nodes()[i]->X()[2];
  }

  // calculate gauss point coordinates, the inverse jacobian and the determinant of the jacobian
  for (int gp = 0; gp < numgpt_; ++gp)
  {
    // gauss point coordinates
    const double* gpcoord = intpoints_.Point(gp);
    for (int idim = 0; idim < numdim_; idim++) xsi_[gp](idim) = gpcoord[idim];

    // get derivative of shape functions w.r.t. parameter coordinates, needed for calculation of the
    // inverse of the jacobian
    LINALG::Matrix<numdim_, numnod_> deriv;
    DRT::UTILS::shape_function_deriv1<distype>(xsi_[gp], deriv);

    // get the inverse of the Jacobian matrix which looks like:
    /*
                 [ X_,r  Y_,r  Z_,r ]^-1
          J^-1 = [ X_,s  Y_,s  Z_,s ]
                 [ X_,t  Y_,t  Z_,t ]
     */

    invJ_[gp].Multiply(deriv, xrefe);
    // here Jacobian is inverted and det(J) is calculated
    detJ_[gp] = invJ_[gp].Invert();

    // make sure determinant of jacobian is positive
    if (detJ_[gp] <= 0.0) dserror("Element Jacobian mapping %10.5e <= 0.0", detJ_[gp]);
  }

  return;
}


template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8, DRT::Element::hex8>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex27, DRT::Element::hex27>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8fbar, DRT::Element::hex8>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet4, DRT::Element::tet4>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet10, DRT::Element::tet10>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_weg6, DRT::Element::wedge6>;
