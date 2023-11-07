/*----------------------------------------------------------------------*/
/*! \file
\brief NStet5 element

\level 3

*----------------------------------------------------------------------*/

#include "baci_lib_discret.H"
#include "baci_lib_utils.H"
#include "baci_linalg_utils_densematrix_inverse.H"
#include "baci_mat_aaaneohooke.H"
#include "baci_mat_elasthyper.H"
#include "baci_mat_micromaterial.H"
#include "baci_mat_stvenantkirchhoff.H"
#include "baci_so3_nstet5.H"
#include "baci_so3_prestress.H"
#include "baci_utils_exceptions.H"

#include <Teuchos_SerialDenseSolver.hpp>


/*----------------------------------------------------------------------*
 |  init the element jacobian mapping (protected)              gee 03/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet5::InitElement()
{
  CORE::LINALG::Matrix<4, 3> xrefe;
  CORE::LINALG::Matrix<4, 3 + 1> J;
  {
    // compute element volume and center node coordinate
    DRT::Node** nodes = Nodes();  // outer nodes only
    for (double& i : midX_) i = 0.0;
    for (int i = 0; i < 4; ++i)
    {
      const auto& x = nodes[i]->X();
      J(i, 0) = 1.0;
      J(i, 1) = xrefe(i, 0) = x[0];
      J(i, 2) = xrefe(i, 1) = x[1];
      J(i, 3) = xrefe(i, 2) = x[2];
      midX_[0] += x[0];
      midX_[1] += x[1];
      midX_[2] += x[2];
    }
    for (double& i : midX_) i /= 4;

    V_ = J.Determinant() / 6.0;
    if (V_ == 0.0)
      dserror("Element volume is zero");
    else if (V_ < 0.0)
      dserror("Element volume is negative");
  }

  //------------------------------------------------------------------subtets
  /* structure of nxyz_:
  **             [   dN_1     dN_1     dN_1   ]
  **             [  ------   ------   ------  ]
  **             [    dX       dY       dZ    ]
  **    nxyz_ =  [     |        |        |    ]
  **             [                            ]
  **             [   dN_4     dN_4     dN_4   ]
  **             [  -------  -------  ------- ]
  **             [    dX       dY       dZ    ]
  */
  const double gploc = 0.25;
  CORE::LINALG::Matrix<4, 1> funct;
  ShapeFunction(funct, gploc, gploc, gploc, gploc);
  CORE::LINALG::Matrix<4, 4> deriv(true);
  ShapeFunctionDerivatives(deriv);
  CORE::LINALG::Matrix<3, 4> tmp;
  CORE::LINALG::Matrix<4, 3> Iaug;  // initialize to zero
  CORE::LINALG::Matrix<4, 3> partials;
  CORE::LINALG::FixedSizeSerialDenseSolver<4, 4, 3> solver;
  // loop i subelements
  for (int i = 0; i < 4; ++i)
  {
    // master tet has node numbering [0 1 2 3]

    // subtets have node numberings  [0 1 2 4]
    //                               [1 3 2 4]
    //                               [0 3 1 4]
    //                               [0 2 3 4]
    for (int k = 0; k < 3; ++k)
    {
      xrefe(0, k) = Nodes()[SubLM(i)[0]]->X()[k];
      xrefe(1, k) = Nodes()[SubLM(i)[1]]->X()[k];
      xrefe(2, k) = Nodes()[SubLM(i)[2]]->X()[k];
      xrefe(3, k) = MidX()[k];
    }

    // volume of subelements
    for (int j = 0; j < 4; ++j)
    {
      J(j, 0) = 1.0;
      J(j, 1) = xrefe(j, 0);
      J(j, 2) = xrefe(j, 1);
      J(j, 3) = xrefe(j, 2);
    }
    subV_[i] = J.Determinant() / 6.0;
    if (subV_[i] == 0.0)
      dserror("NSTET5 %d Subelement %d volume is zero %10.6e", Id(), i, subV_[i]);
    else if (subV_[i] < 0.0)
      dserror("NSTET5 %d Subelement %d volume is negative %10.6e", Id(), i, subV_[i]);

    // spatial derivatives of shape functions
    tmp.MultiplyTN(xrefe, deriv);
    for (int j = 0; j < 4; j++) J(0, j) = 1;
    for (int row = 0; row < 3; row++)
      for (int col = 0; col < 4; col++) J(row + 1, col) = tmp(row, col);
    Iaug = 0.0;
    Iaug(1, 0) = 1;
    Iaug(2, 1) = 1;
    Iaug(3, 2) = 1;
    partials = 0.0;
    solver.SetMatrix(J);
    solver.SetVectors(partials, Iaug);
    solver.FactorWithEquilibration(true);
    int err = solver.Factor();
    int err2 = solver.Solve();
    if (err || err2) dserror("Inversion of Jacobian failed");
    subnxyz_[i].Multiply(deriv, partials);
  }  // for (int i=0; i<4; ++i)


  return;
}  // DRT::ELEMENTS::NStet5::InitElement



/*----------------------------------------------------------------------*
 |  evaluate the element (public)                              gee 03/12|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::NStet5::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra)
{
  CORE::LINALG::Matrix<15, 15> elemat1(elemat1_epetra.values(), true);
  CORE::LINALG::Matrix<15, 15> elemat2(elemat2_epetra.values(), true);
  CORE::LINALG::Matrix<15, 1> elevec1(elevec1_epetra.values(), true);

  // start with "none"
  DRT::ELEMENTS::NStet5::ActionType act = NStet5::none;

  // get the required action
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    dserror("No action supplied");
  else if (action == "calc_struct_linstiff")
    act = NStet5::calc_struct_linstiff;
  else if (action == "calc_struct_nlnstiff")
    act = NStet5::calc_struct_nlnstiff;
  else if (action == "calc_struct_internalforce")
    act = NStet5::calc_struct_internalforce;
  else if (action == "calc_struct_linstiffmass")
    act = NStet5::calc_struct_linstiffmass;
  else if (action == "calc_struct_nlnstiffmass")
    act = NStet5::calc_struct_nlnstiffmass;
  else if (action == "calc_struct_nlnstifflmass")
    act = NStet5::calc_struct_nlnstifflmass;
  else if (action == "calc_struct_stress")
    act = NStet5::calc_struct_stress;
  else if (action == "calc_struct_eleload")
    act = NStet5::calc_struct_eleload;
  else if (action == "calc_struct_fsiload")
    act = NStet5::calc_struct_fsiload;
  else if (action == "calc_struct_update_istep")
    act = NStet5::calc_struct_update_istep;
  else if (action == "calc_struct_reset_istep")
    act = NStet5::calc_struct_reset_istep;
  else if (action == "multi_calc_dens")
    act = NStet5::multi_calc_dens;
  else if (action == "multi_readrestart")
    act = NStet5::multi_readrestart;
  else if (action == "calc_struct_recover")
    return 0;
  else if (action == "calc_struct_predict")
    return 0;
  else
    dserror("Unknown type of action for NStet5");

  // what should the element do
  switch (act)
  {
    //==================================================================================
    // nonlinear stiffness, internal force vector, and consistent mass matrix
    case calc_struct_nlnstiffmass:
    case calc_struct_nlnstifflmass:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      nstet5nlnstiffmass(lm, mydisp, &elemat1, &elemat2, &elevec1, nullptr, nullptr,
          INPAR::STR::stress_none, INPAR::STR::strain_none);
      if (act == calc_struct_nlnstifflmass) nstet5lumpmass(&elemat2);
    }
    break;

    //==================================================================================
    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      CORE::LINALG::Matrix<15, 15>* elemat1ptr = nullptr;
      if (elemat1.IsInitialized()) elemat1ptr = &elemat1;
      nstet5nlnstiffmass(lm, mydisp, elemat1ptr, nullptr, &elevec1, nullptr, nullptr,
          INPAR::STR::stress_none, INPAR::STR::strain_none);
    }
    break;

    //==================================================================================
    // internal force vector only
    case calc_struct_internalforce:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      // CORE::LINALG::Matrix<15,15> myemat(true);
      // nstet5nlnstiffmass(lm,mydisp,&myemat,nullptr,&elevec1,
      //                  nullptr,nullptr,INPAR::STR::stress_none,INPAR::STR::strain_none);
      nstet5nlnstiffmass(lm, mydisp, nullptr, nullptr, &elevec1, nullptr, nullptr,
          INPAR::STR::stress_none, INPAR::STR::strain_none);
    }
    break;

    //==================================================================================
    // evaluate stresses and strains at gauss point
    case calc_struct_stress:
    {
      // nothing to do for ghost elements
      if (discretization.Comm().MyPID() == Owner())
      {
        //------------------------------- compute element stress from stabilization
        Teuchos::RCP<std::vector<char>> stressdata =
            params.get<Teuchos::RCP<std::vector<char>>>("stress", Teuchos::null);
        Teuchos::RCP<std::vector<char>> straindata =
            params.get<Teuchos::RCP<std::vector<char>>>("strain", Teuchos::null);
        if (stressdata == Teuchos::null) dserror("Cannot get stress 'data'");
        if (straindata == Teuchos::null) dserror("Cannot get strain 'data'");
        auto iostress =
            DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress", INPAR::STR::stress_none);
        auto iostrain =
            DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain", INPAR::STR::strain_none);
        Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
        if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
        std::vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
        CORE::LINALG::Matrix<1, 6> stress(true);
        CORE::LINALG::Matrix<1, 6> strain(true);
        CORE::LINALG::Matrix<1, 6> elestress(true);
        CORE::LINALG::Matrix<1, 6> elestrain(true);
        nstet5nlnstiffmass(
            lm, mydisp, nullptr, nullptr, nullptr, &elestress, &elestrain, iostress, iostrain);

        //--------------------------------- interpolate nodal stress from every node
        Teuchos::RCP<Epetra_MultiVector> nodestress = ElementType().nstress_;
        Teuchos::RCP<Epetra_MultiVector> nodestrain = ElementType().nstrain_;
        const int numnode = NumNode();
        for (int i = 0; i < numnode; ++i)
        {
          const int gid = Nodes()[i]->Id();
          const int lid = nodestress->Map().LID(gid);
          if (lid == -1) dserror("Cannot find matching nodal stresses/strains");
          for (int j = 0; j < 6; ++j)
          {
            stress(0, j) += (*(*nodestress)(j))[lid];
            strain(0, j) += (*(*nodestrain)(j))[lid];
          }
        }

        for (int j = 0; j < 6; ++j)
        {
          stress(0, j) /= numnode;
          strain(0, j) /= numnode;
        }

        //-------------------------------------------------------- add element stress

        for (int j = 0; j < 6; ++j)
        {
          stress(0, j) += elestress(0, j);
          strain(0, j) += elestrain(0, j);
        }

        //----------------------------------------------- add final stress to storage
        {
          DRT::PackBuffer data;
          AddtoPack(data, stress);
          data.StartPacking();
          AddtoPack(data, stress);
          std::copy(data().begin(), data().end(), std::back_inserter(*stressdata));
        }
        {
          DRT::PackBuffer data;
          AddtoPack(data, strain);
          data.StartPacking();
          AddtoPack(data, strain);
          std::copy(data().begin(), data().end(), std::back_inserter(*straindata));
        }
      }  // if (discretization.Comm().MyPID()==Owner())
    }
    break;

    //==================================================================================
    case calc_struct_eleload:
      dserror("this class is not supposed to evaluate a load, use EvaluateNeumann(...)");
      break;

    //==================================================================================
    case calc_struct_fsiload:
      dserror("Case not yet implemented");
      break;

    //==================================================================================
    case calc_struct_update_istep:
    {
      Teuchos::RCP<MAT::Material> mat = Material();
      if (mat->MaterialType() == INPAR::MAT::m_struct_multiscale)
      {
        auto* micro = dynamic_cast<MAT::MicroMaterial*>(mat.get());
        micro->Update();
      }
    }
    break;

    //==================================================================================
    case calc_struct_reset_istep:
    {
      ;  // there is nothing to do here at the moment
    }
    break;

    //==================================================================================
    // linear stiffness and consistent mass matrix
    case calc_struct_linstiffmass:
      dserror("Case 'calc_struct_linstiffmass' not implemented");
      break;

    //==================================================================================
    // linear stiffness
    case calc_struct_linstiff:
      dserror("action calc_struct_linstiff currently not supported");
      break;

    //==================================================================================
    case multi_calc_dens:
    {
      nstet5_homog(params);
    }
    break;

    //==================================================================================
    case multi_readrestart:
    {
      nstet5_read_restart_multi();
    }
    break;

    //==================================================================================
    default:
      dserror("Unknown type of action for NStet5");
  }
  return 0;
}


/*----------------------------------------------------------------------*
 |  evaluate the element (private)                             gee 03/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet5::nstet5nlnstiffmass(std::vector<int>& lm,  // location matrix
    std::vector<double>& disp,                                        // current displacements
    CORE::LINALG::Matrix<15, 15>* stiffmatrix,                        // element stiffness matrix
    CORE::LINALG::Matrix<15, 15>* massmatrix,                         // element mass matrix
    CORE::LINALG::Matrix<15, 1>* force,                               // stress output options
    CORE::LINALG::Matrix<1, 6>* elestress,                            // stress output
    CORE::LINALG::Matrix<1, 6>* elestrain,                            // strain output
    const INPAR::STR::StressType iostress,                            // type of stress
    const INPAR::STR::StrainType iostrain)                            // type of strain
{
  if (elestrain) (*elestrain) = 0.0;
  if (elestress) (*elestress) = 0.0;


  for (int sub = 0; sub < 4; ++sub)  // loop subelements
  {
    // subelement deformation gradient previously computed in PreEvaluate
    CORE::LINALG::Matrix<3, 3>& F = SubF(sub);

    //--------------------------- Right Cauchy-Green tensor C = = F^T * F
    CORE::LINALG::Matrix<3, 3> cauchygreen;
    cauchygreen.MultiplyTN(F, F);

    // --Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    CORE::LINALG::Matrix<6, 1> glstrain;
    glstrain(0) = 0.5 * (cauchygreen(0, 0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1, 1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2, 2) - 1.0);
    glstrain(3) = cauchygreen(0, 1);
    glstrain(4) = cauchygreen(1, 2);
    glstrain(5) = cauchygreen(2, 0);

    //------------------------------------ B=operator ( same as in hex8 !)
    /*
    ** B = F : N,xyz
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

    // 6x15 n_stresses * number degrees of freedom per element
    CORE::LINALG::Matrix<6, 12> bop;
    const CORE::LINALG::Matrix<4, 3>& nxyz = SubNxyz(sub);
    for (int i = 0; i < 4; i++)
    {
      bop(0, 3 * i + 0) = F(0, 0) * nxyz(i, 0);
      bop(0, 3 * i + 1) = F(1, 0) * nxyz(i, 0);
      bop(0, 3 * i + 2) = F(2, 0) * nxyz(i, 0);
      bop(1, 3 * i + 0) = F(0, 1) * nxyz(i, 1);
      bop(1, 3 * i + 1) = F(1, 1) * nxyz(i, 1);
      bop(1, 3 * i + 2) = F(2, 1) * nxyz(i, 1);
      bop(2, 3 * i + 0) = F(0, 2) * nxyz(i, 2);
      bop(2, 3 * i + 1) = F(1, 2) * nxyz(i, 2);
      bop(2, 3 * i + 2) = F(2, 2) * nxyz(i, 2);
      /* ~~~ */
      bop(3, 3 * i + 0) = F(0, 0) * nxyz(i, 1) + F(0, 1) * nxyz(i, 0);
      bop(3, 3 * i + 1) = F(1, 0) * nxyz(i, 1) + F(1, 1) * nxyz(i, 0);
      bop(3, 3 * i + 2) = F(2, 0) * nxyz(i, 1) + F(2, 1) * nxyz(i, 0);
      bop(4, 3 * i + 0) = F(0, 1) * nxyz(i, 2) + F(0, 2) * nxyz(i, 1);
      bop(4, 3 * i + 1) = F(1, 1) * nxyz(i, 2) + F(1, 2) * nxyz(i, 1);
      bop(4, 3 * i + 2) = F(2, 1) * nxyz(i, 2) + F(2, 2) * nxyz(i, 1);
      bop(5, 3 * i + 0) = F(0, 2) * nxyz(i, 0) + F(0, 0) * nxyz(i, 2);
      bop(5, 3 * i + 1) = F(1, 2) * nxyz(i, 0) + F(1, 0) * nxyz(i, 2);
      bop(5, 3 * i + 2) = F(2, 2) * nxyz(i, 0) + F(2, 0) * nxyz(i, 2);
    }

    //------------------------------------------------- call material law
    CORE::LINALG::Matrix<6, 6> cmat(true);
    CORE::LINALG::Matrix<6, 1> stress(true);
    double density = -999.99;
#ifndef PUSO_NSTET5
    {
      SelectMaterial(stress, cmat, density, glstrain, F, 0);

      // define stuff we need to do the split
      CORE::LINALG::Matrix<6, 6> cmatdev;
      CORE::LINALG::Matrix<6, 1> stressdev;

      // do just the deviatoric components
      NStet5Type::DevStressTangent(stressdev, cmatdev, cmat, stress, cauchygreen);
      stress = stressdev;
      cmat = cmatdev;

      stress.Scale(ALPHA_NSTET5);
      cmat.Scale(ALPHA_NSTET5);
    }
#else
    {
      SelectMaterial(stress, cmat, density, glstrain, F, 0);
      stress.Scale(ALPHA_NSTET5);
      cmat.Scale(ALPHA_NSTET5);
      glstrainbar = glstrain;
    }
#endif

    //---------------------------------------------- output of stress and strain
    {
      CORE::LINALG::Matrix<6, 1> glstrainbar(false);
      if (iostrain != INPAR::STR::strain_none)
      {
        // do deviatoric F, C, E
        const double J = F.Determinant();
        CORE::LINALG::Matrix<3, 3> Cbar(cauchygreen);
        Cbar.Scale(pow(J, -2. / 3.));
        glstrainbar(0) = 0.5 * (Cbar(0, 0) - 1.0);
        glstrainbar(1) = 0.5 * (Cbar(1, 1) - 1.0);
        glstrainbar(2) = 0.5 * (Cbar(2, 2) - 1.0);
        glstrainbar(3) = Cbar(0, 1);
        glstrainbar(4) = Cbar(1, 2);
        glstrainbar(5) = Cbar(2, 0);
      }
      //-----------------------------------------------------------------strain
      switch (iostrain)
      {
        case INPAR::STR::strain_gl:
        {
          if (elestrain == nullptr) dserror("no strain data available");
          for (int i = 0; i < 3; ++i)
            (*elestrain)(0, i) += (SubV(sub) / Vol() * ALPHA_NSTET5 * glstrainbar(i));
          for (int i = 0; i < 3; ++i)
            (*elestrain)(0, i) += (SubV(sub) / Vol() * ALPHA_NSTET5 * 0.5 * glstrainbar(i));
        }
        break;
        case INPAR::STR::strain_ea:
        {
          if (elestrain == nullptr) dserror("no strain data available");
          CORE::LINALG::Matrix<3, 3> gl;
          gl(0, 0) = glstrainbar(0);  // divide off-diagonals by 2
          gl(0, 1) = 0.5 * glstrainbar(3);
          gl(0, 2) = 0.5 * glstrainbar(5);
          gl(1, 0) = gl(0, 1);
          gl(1, 1) = glstrainbar(1);
          gl(1, 2) = 0.5 * glstrainbar(4);
          gl(2, 0) = gl(0, 2);
          gl(2, 1) = gl(1, 2);
          gl(2, 2) = glstrainbar(2);

          CORE::LINALG::Matrix<3, 3> Fbar(true);
          Fbar.SetCopy(F.A());
#ifndef PUSO_NSTET5
          Fbar.Scale(pow(F.Determinant(), -1. / 3.));
#endif
          CORE::LINALG::Matrix<3, 3> invdefgrd;
          invdefgrd.Invert(Fbar);

          CORE::LINALG::Matrix<3, 3> temp;
          CORE::LINALG::Matrix<3, 3> euler_almansi;
          temp.Multiply(gl, invdefgrd);
          euler_almansi.MultiplyTN(invdefgrd, temp);


          (*elestrain)(0, 0) += (SubV(sub) / Vol() * ALPHA_NSTET5 * euler_almansi(0, 0));
          (*elestrain)(0, 1) += (SubV(sub) / Vol() * ALPHA_NSTET5 * euler_almansi(1, 1));
          (*elestrain)(0, 2) += (SubV(sub) / Vol() * ALPHA_NSTET5 * euler_almansi(2, 2));
          (*elestrain)(0, 3) += (SubV(sub) / Vol() * ALPHA_NSTET5 * euler_almansi(0, 1));
          (*elestrain)(0, 4) += (SubV(sub) / Vol() * ALPHA_NSTET5 * euler_almansi(1, 2));
          (*elestrain)(0, 5) += (SubV(sub) / Vol() * ALPHA_NSTET5 * euler_almansi(0, 2));
        }
        break;
        case INPAR::STR::strain_none:
          break;
        default:
          dserror("requested strain option not available");
      }
      //-----------------------------------------------------------------stress
      switch (iostress)
      {
        case INPAR::STR::stress_2pk:
        {
          if (elestress == nullptr) dserror("no stress data available");

          for (int i = 0; i < 6; ++i)
            (*elestress)(0, i) += (SubV(sub) / Vol() * stress(i));  // ALPHA_NSTET already in stress
        }
        break;
        case INPAR::STR::stress_cauchy:
        {
          if (elestress == nullptr) dserror("no stress data available");

          CORE::LINALG::Matrix<3, 3> pkstress;
          pkstress(0, 0) = stress(0);  // ALPHA_NSTET already in stress
          pkstress(0, 1) = stress(3);
          pkstress(0, 2) = stress(5);
          pkstress(1, 0) = pkstress(0, 1);
          pkstress(1, 1) = stress(1);
          pkstress(1, 2) = stress(4);
          pkstress(2, 0) = pkstress(0, 2);
          pkstress(2, 1) = pkstress(1, 2);
          pkstress(2, 2) = stress(2);

          CORE::LINALG::Matrix<3, 3> temp;
          CORE::LINALG::Matrix<3, 3> cauchystress;

          CORE::LINALG::Matrix<3, 3> Fbar(true);
          Fbar.SetCopy(F.A());
#ifndef PUSO_NSTET5
          Fbar.Scale(pow(F.Determinant(), -1. / 3.));
#endif
          temp.Multiply(1.0 / Fbar.Determinant(), Fbar, pkstress);
          cauchystress.MultiplyNT(temp, Fbar);

          (*elestress)(0, 0) += (SubV(sub) / Vol() * cauchystress(0, 0));
          (*elestress)(0, 1) += (SubV(sub) / Vol() * cauchystress(1, 1));
          (*elestress)(0, 2) += (SubV(sub) / Vol() * cauchystress(2, 2));
          (*elestress)(0, 3) += (SubV(sub) / Vol() * cauchystress(0, 1));
          (*elestress)(0, 4) += (SubV(sub) / Vol() * cauchystress(1, 2));
          (*elestress)(0, 5) += (SubV(sub) / Vol() * cauchystress(0, 2));
        }
        break;
        case INPAR::STR::stress_none:
          break;
        default:
          dserror("requested stress type not available");
      }
    }

    //----------------------------------------------- internal force and tangent
    // update internal force vector
    if (force)
    {
      CORE::LINALG::Matrix<12, 1> subforce(true);
      // integrate internal force vector f = f + (B^T . sigma) * V
      subforce.MultiplyTN(SubV(sub), bop, stress, 0.0);

      for (int i = 0; i < 4; ++i)  // loop 4 nodes of subelement
      {
        (*force)(SubLM(sub)[i] * 3 + 0) += subforce(i * 3 + 0);
        (*force)(SubLM(sub)[i] * 3 + 1) += subforce(i * 3 + 1);
        (*force)(SubLM(sub)[i] * 3 + 2) += subforce(i * 3 + 2);
      }
    }  // if (force)

    // update stiffness matrix
    if (stiffmatrix)
    {
      const double V = SubV(sub);
      CORE::LINALG::Matrix<12, 12> substiffmatrix(true);
      // integrate elastic stiffness matrix
      // keu = keu + (B^T . C . B) * V
      CORE::LINALG::Matrix<6, 12> cb;
      cb.Multiply(cmat, bop);
      substiffmatrix.MultiplyTN(V, bop, cb, 0.0);

      // integrate `geometric' stiffness matrix and add to keu
      double sBL[3];
      for (int i = 0; i < 4; ++i)
      {
        sBL[0] = V * (stress(0) * nxyz(i, 0) + stress(3) * nxyz(i, 1) + stress(5) * nxyz(i, 2));
        sBL[1] = V * (stress(3) * nxyz(i, 0) + stress(1) * nxyz(i, 1) + stress(4) * nxyz(i, 2));
        sBL[2] = V * (stress(5) * nxyz(i, 0) + stress(4) * nxyz(i, 1) + stress(2) * nxyz(i, 2));
        for (int j = 0; j < 4; ++j)
        {
          double BsB = 0.0;
          for (int dim = 0; dim < 3; ++dim) BsB += nxyz(j, dim) * sBL[dim];
          substiffmatrix(3 * i + 0, 3 * j + 0) += BsB;
          substiffmatrix(3 * i + 1, 3 * j + 1) += BsB;
          substiffmatrix(3 * i + 2, 3 * j + 2) += BsB;
        }
      }

      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
        {
          (*stiffmatrix)(SubLM(sub)[i] * 3 + 0, SubLM(sub)[j] * 3 + 0) +=
              substiffmatrix(i * 3 + 0, j * 3 + 0);
          (*stiffmatrix)(SubLM(sub)[i] * 3 + 0, SubLM(sub)[j] * 3 + 1) +=
              substiffmatrix(i * 3 + 0, j * 3 + 1);
          (*stiffmatrix)(SubLM(sub)[i] * 3 + 0, SubLM(sub)[j] * 3 + 2) +=
              substiffmatrix(i * 3 + 0, j * 3 + 2);

          (*stiffmatrix)(SubLM(sub)[i] * 3 + 1, SubLM(sub)[j] * 3 + 0) +=
              substiffmatrix(i * 3 + 1, j * 3 + 0);
          (*stiffmatrix)(SubLM(sub)[i] * 3 + 1, SubLM(sub)[j] * 3 + 1) +=
              substiffmatrix(i * 3 + 1, j * 3 + 1);
          (*stiffmatrix)(SubLM(sub)[i] * 3 + 1, SubLM(sub)[j] * 3 + 2) +=
              substiffmatrix(i * 3 + 1, j * 3 + 2);

          (*stiffmatrix)(SubLM(sub)[i] * 3 + 2, SubLM(sub)[j] * 3 + 0) +=
              substiffmatrix(i * 3 + 2, j * 3 + 0);
          (*stiffmatrix)(SubLM(sub)[i] * 3 + 2, SubLM(sub)[j] * 3 + 1) +=
              substiffmatrix(i * 3 + 2, j * 3 + 1);
          (*stiffmatrix)(SubLM(sub)[i] * 3 + 2, SubLM(sub)[j] * 3 + 2) +=
              substiffmatrix(i * 3 + 2, j * 3 + 2);
        }
    }  // if (stiffmatrix)

    if (massmatrix)
    {
      CORE::LINALG::Matrix<12, 12> submassmatrix;
      submassmatrix = 0.0;

      // for mass matrix use a 4 gauss points integration:
      // ( 1 gauss point is not enough!)
      const double alpha = (5.0 + 3.0 * sqrt(5.0)) / 20.0;
      const double beta = (5.0 - sqrt(5.0)) / 20.0;
      const double weight = 0.25;
      const double V = SubV(sub);
      double xsi[4][4];
      xsi[0][0] = alpha;
      xsi[0][1] = beta;
      xsi[0][2] = beta;
      xsi[0][3] = beta;
      xsi[1][0] = beta;
      xsi[1][1] = alpha;
      xsi[1][2] = beta;
      xsi[1][3] = beta;
      xsi[2][0] = beta;
      xsi[2][1] = beta;
      xsi[2][2] = alpha;
      xsi[2][3] = beta;
      xsi[3][0] = beta;
      xsi[3][1] = beta;
      xsi[3][2] = beta;
      xsi[3][3] = alpha;
      for (auto& gp : xsi)
      {
        CORE::LINALG::Matrix<4, 1> funct;
        ShapeFunction(funct, gp[0], gp[1], gp[2], gp[3]);
        const double f = density * V * weight;
        for (int i = 0; i < 4; ++i)
          for (int j = 0; j < 4; ++j)
          {
            const double fac = funct(i) * funct(j) * f;
            submassmatrix(3 * i + 0, 3 * j + 0) += fac;
            submassmatrix(3 * i + 1, 3 * j + 1) += fac;
            submassmatrix(3 * i + 2, 3 * j + 2) += fac;
          }
      }  // for (int gp=0; gp<4; ++gp)
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
        {
          (*massmatrix)(SubLM(sub)[i] * 3 + 0, SubLM(sub)[j] * 3 + 0) +=
              submassmatrix(i * 3 + 0, j * 3 + 0);
          (*massmatrix)(SubLM(sub)[i] * 3 + 0, SubLM(sub)[j] * 3 + 1) +=
              submassmatrix(i * 3 + 0, j * 3 + 1);
          (*massmatrix)(SubLM(sub)[i] * 3 + 0, SubLM(sub)[j] * 3 + 2) +=
              submassmatrix(i * 3 + 0, j * 3 + 2);

          (*massmatrix)(SubLM(sub)[i] * 3 + 1, SubLM(sub)[j] * 3 + 0) +=
              submassmatrix(i * 3 + 1, j * 3 + 0);
          (*massmatrix)(SubLM(sub)[i] * 3 + 1, SubLM(sub)[j] * 3 + 1) +=
              submassmatrix(i * 3 + 1, j * 3 + 1);
          (*massmatrix)(SubLM(sub)[i] * 3 + 1, SubLM(sub)[j] * 3 + 2) +=
              submassmatrix(i * 3 + 1, j * 3 + 2);

          (*massmatrix)(SubLM(sub)[i] * 3 + 2, SubLM(sub)[j] * 3 + 0) +=
              submassmatrix(i * 3 + 2, j * 3 + 0);
          (*massmatrix)(SubLM(sub)[i] * 3 + 2, SubLM(sub)[j] * 3 + 1) +=
              submassmatrix(i * 3 + 2, j * 3 + 1);
          (*massmatrix)(SubLM(sub)[i] * 3 + 2, SubLM(sub)[j] * 3 + 2) +=
              submassmatrix(i * 3 + 2, j * 3 + 2);
        }
    }  // if (massmatrix)
  }    // for (int sub=0; sub<4; ++sub) loop subelements


  return;
}  // DRT::ELEMENTS::NStet5::nstet5nlnstiffmass


/*----------------------------------------------------------------------*
 |  lump mass matrix                                           gee 03/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet5::nstet5lumpmass(CORE::LINALG::Matrix<15, 15>* emass)
{
  // lump mass matrix
  if (emass != nullptr)
  {
    // we assume #elemat2 is a square matrix
    for (unsigned c = 0; c < (*emass).numCols(); ++c)  // parse columns
    {
      double d = 0.0;
      for (unsigned r = 0; r < (*emass).numRows(); ++r)  // parse rows
      {
        d += (*emass)(r, c);  // accumulate row entries
        (*emass)(r, c) = 0.0;
      }
      (*emass)(c, c) = d;  // apply sum of row entries on diagonal
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | material laws for NStet5 (protected)                        gee 03/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet5::SelectMaterial(CORE::LINALG::Matrix<6, 1>& stress,
    CORE::LINALG::Matrix<6, 6>& cmat, double& density, CORE::LINALG::Matrix<6, 1>& glstrain,
    CORE::LINALG::Matrix<3, 3>& defgrd, int gp)
{
  CORE::LINALG::SerialDenseVector stress_e(Teuchos::View, stress.A(), stress.numRows());
  CORE::LINALG::SerialDenseMatrix cmat_e(
      Teuchos::View, cmat.A(), cmat.numRows(), cmat.numRows(), cmat.numCols());
  const CORE::LINALG::SerialDenseVector glstrain_e(Teuchos::View, glstrain.A(), glstrain.numRows());

  Teuchos::RCP<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    case INPAR::MAT::m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      auto* stvk = dynamic_cast<MAT::StVenantKirchhoff*>(mat.get());
      stvk->Evaluate(&glstrain_e, &cmat_e, &stress_e);
      density = stvk->Density();
    }
    break;
    case INPAR::MAT::m_aaaneohooke: /*-- special case of generalised NeoHookean material see
                                       Raghavan, Vorp */
    {
      auto* aaa = dynamic_cast<MAT::AAAneohooke*>(mat.get());
      Teuchos::ParameterList params;
      aaa->Evaluate(&defgrd, &glstrain, params, &stress, &cmat, gp, Id());
      density = aaa->Density();
    }
    break;
    case INPAR::MAT::m_elasthyper: /*----------- general hyperelastic matrial */
    {
      auto* hyper = dynamic_cast<MAT::ElastHyper*>(mat.get());
      Teuchos::ParameterList params;
      hyper->Evaluate(&defgrd, &glstrain, params, &stress, &cmat, gp, Id());
      density = hyper->Density();
      return;
      break;
    }
    default:
      dserror("Illegal type %d of material for element NStet5 tet4", mat->MaterialType());
      break;
  }

  /*--------------------------------------------------------------------*/
  return;
}  // DRT::ELEMENTS::NStet5::SelectMaterial



/*----------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)     gee 03/12|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::NStet5::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseMatrix* elemat1)
{
  dserror("DRT::ELEMENTS::NStet5::EvaluateNeumann not implemented");
  return -1;
}  // DRT::ELEMENTS::NStet5::EvaluateNeumann
