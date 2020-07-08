/*----------------------------------------------------------------------------*/
/*! \file

\brief three dimensional nonlinear corotational Reissner beam element

\level 2

*/
/*----------------------------------------------------------------------------*/


#include "beam3.H"
#include "beam3_base.H"
#include "beam3r.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_utils_sparse_algebra_math.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_fem_general/largerotations.H"
#include "../drt_structure_new/str_elements_paramsinterface.H"


/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public) cyron 01/08|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  SetParamsInterfacePtr(params);

  if (IsParamsInterface()) SetBrownianDynParamsInterfacePtr();

  // start with "none"
  ELEMENTS::ActionType act = ELEMENTS::none;

  if (IsParamsInterface())
  {
    act = ParamsInterface().GetActionType();
  }
  else
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
      act = ELEMENTS::struct_calc_nlnstifflmass;  // with lumped mass matrix
    else if (action == "calc_struct_stress")
      act = ELEMENTS::struct_calc_stress;
    else if (action == "calc_struct_eleload")
      act = ELEMENTS::struct_calc_eleload;
    else if (action == "calc_struct_fsiload")
      act = ELEMENTS::struct_calc_fsiload;
    else if (action == "calc_struct_update_istep")
      act = ELEMENTS::struct_calc_update_istep;
    else if (action == "calc_struct_reset_istep")
      act = ELEMENTS::struct_calc_reset_istep;
    else if (action == "calc_struct_ptcstiff")
      act = ELEMENTS::struct_calc_ptcstiff;
    else if (action == "calc_struct_energy")
      act = ELEMENTS::struct_calc_energy;
    else
      dserror("Unknown type of action for Beam3");
  }

  switch (act)
  {
    case ELEMENTS::struct_calc_ptcstiff:
    {
      const int nnode = NumNode();

      switch (nnode)
      {
        case 2:
          EvaluatePTC<2>(params, elemat1);
          break;
        case 3:
          EvaluatePTC<3>(params, elemat1);
          break;
        case 4:
          EvaluatePTC<4>(params, elemat1);
          break;
        case 5:
          EvaluatePTC<5>(params, elemat1);
          break;
        case 6:
          EvaluatePTC<6>(params, elemat1);
          break;
        default:
          dserror("Only Line2, Line3, Line4, Line5 and Line6 Elements implemented.");
      }
    }
    break;
    /*in case that only linear stiffness matrix is required b3_nlstiffmass is called with zero
     dispalcement and residual values*/
    case ELEMENTS::struct_calc_linstiff:
    {
      // only nonlinear case implemented!
      dserror("linear stiffness matrix called, but not implemented");
    }
    break;
    case ELEMENTS::struct_calc_energy:
    {
      if (elevec1 != Teuchos::null)  // old structural time integration
      {
        if (elevec1.M() != 1)
          dserror(
              "energy vector of invalid size %i, expected row dimension 1 (total elastic energy of "
              "element)!",
              elevec1.M());

        // need current global displacement and and get them from discretization
        // making use of the local-to-global map lm one can extract current displacemnet and
        // residual values for each degree of freedom
        Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
        if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
        std::vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

        const int nnode = NumNode();

        switch (nnode)
        {
          case 2:
          {
            b3_energy<2>(params, mydisp, &elevec1);
            break;
          }
          case 3:
          {
            b3_energy<3>(params, mydisp, &elevec1);
            break;
          }
          case 4:
          {
            b3_energy<4>(params, mydisp, &elevec1);
            break;
          }
          case 5:
          {
            b3_energy<5>(params, mydisp, &elevec1);
            break;
          }
          case 6:
          {
            b3_energy<6>(params, mydisp, &elevec1);
            break;
          }
          default:
            dserror("Only Line2, Line3, Line4, Line5 and Line6 Elements implemented.");
        }
      }
      else if (IsParamsInterface())  // new structural time integration
      {
        dserror("not implemented yet");
      }

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
      // making use of the local-to-global map lm one can extract current displacemnet and residual
      // values for each degree of freedom
      //
      // get element displcements
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

      // get residual displacements
      Teuchos::RCP<const Epetra_Vector> res = discretization.GetState("residual displacement");
      if (res == Teuchos::null) dserror("Cannot get state vectors 'residual displacement'");
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res, myres, lm);

      const int nnode = NumNode();

      if (act == ELEMENTS::struct_calc_nlnstiffmass)
      {
        switch (nnode)
        {
          case 2:
          {
            b3_nlnstiffmass<2>(params, mydisp, &elemat1, &elemat2, &elevec1);
            break;
          }
          case 3:
          {
            b3_nlnstiffmass<3>(params, mydisp, &elemat1, &elemat2, &elevec1);
            break;
          }
          case 4:
          {
            b3_nlnstiffmass<4>(params, mydisp, &elemat1, &elemat2, &elevec1);
            break;
          }
          case 5:
          {
            b3_nlnstiffmass<5>(params, mydisp, &elemat1, &elemat2, &elevec1);
            break;
          }
          case 6:
          {
            b3_nlnstiffmass<6>(params, mydisp, &elemat1, &elemat2, &elevec1);
            break;
          }
          default:
            dserror("Only Line2, Line3, Line4, Line5 and Line6 Elements implemented.");
        }
      }
      else if (act == ELEMENTS::struct_calc_nlnstifflmass)
      {
        switch (nnode)
        {
          case 2:
          {
            b3_nlnstiffmass<2>(params, mydisp, &elemat1, &elemat2, &elevec1);
            lumpmass<2>(&elemat2);
            break;
          }
          case 3:
          {
            b3_nlnstiffmass<3>(params, mydisp, &elemat1, &elemat2, &elevec1);
            lumpmass<3>(&elemat2);
            break;
          }
          case 4:
          {
            b3_nlnstiffmass<4>(params, mydisp, &elemat1, &elemat2, &elevec1);
            lumpmass<4>(&elemat2);
            break;
          }
          case 5:
          {
            b3_nlnstiffmass<5>(params, mydisp, &elemat1, &elemat2, &elevec1);
            lumpmass<5>(&elemat2);
            break;
          }
          case 6:
          {
            b3_nlnstiffmass<6>(params, mydisp, &elemat1, &elemat2, &elevec1);
            lumpmass<6>(&elemat2);
            break;
          }
          default:
            dserror("Only Line2, Line3, Line4, Line5 and Line6 Elements implemented.");
        }
      }
      else if (act == ELEMENTS::struct_calc_nlnstiff)
      {
        switch (nnode)
        {
          case 2:
          {
            b3_nlnstiffmass<2>(params, mydisp, &elemat1, NULL, &elevec1);
            break;
          }
          case 3:
          {
            b3_nlnstiffmass<3>(params, mydisp, &elemat1, NULL, &elevec1);
            break;
          }
          case 4:
          {
            b3_nlnstiffmass<4>(params, mydisp, &elemat1, NULL, &elevec1);
            break;
          }
          case 5:
          {
            b3_nlnstiffmass<5>(params, mydisp, &elemat1, NULL, &elevec1);
            break;
          }
          case 6:
          {
            b3_nlnstiffmass<6>(params, mydisp, &elemat1, NULL, &elevec1);
            break;
          }
          default:
            dserror("Only Line2, Line3, Line4, Line5 and Line6 Elements implemented.");
        }
      }

      else if (act == ELEMENTS::struct_calc_internalforce)
      {
        switch (nnode)
        {
          case 2:
          {
            b3_nlnstiffmass<2>(params, mydisp, NULL, NULL, &elevec1);
            break;
          }
          case 3:
          {
            b3_nlnstiffmass<3>(params, mydisp, NULL, NULL, &elevec1);
            break;
          }
          case 4:
          {
            b3_nlnstiffmass<4>(params, mydisp, NULL, NULL, &elevec1);
            break;
          }
          case 5:
          {
            b3_nlnstiffmass<5>(params, mydisp, NULL, NULL, &elevec1);
            break;
          }
          case 6:
          {
            b3_nlnstiffmass<6>(params, mydisp, NULL, NULL, &elevec1);
            break;
          }
          default:
            dserror("Only Line2, Line3, Line4, Line5 and Line6 Elements implemented.");
        }
      }

      /*at the end of an iteration step the geometric configuration has to be updated: the starting
       * point for the next iteration step is the configuration at the end of the current step */
      Qold_ = Qnew_;
      QoldNode_ = QnewNode_;
      curvold_ = curvnew_;
      thetaold_ = thetanew_;
      ThetaOldNode_ = ThetaNewNode_;
      thetaprimeold_ = thetaprimenew_;

      /*
          //the following code block can be used to check quickly whether the nonlinear stiffness
         matrix is calculated
          //correctly or not by means of a numerically approximated stiffness matrix
          //The code block will work for all higher order elements.
          if(Id() == 3) //limiting the following tests to certain element numbers
          {
            //variable to store numerically approximated stiffness matrix
            Epetra_SerialDenseMatrix stiff_approx;
            stiff_approx.Shape(6*nnode,6*nnode);


            //relative error of numerically approximated stiffness matrix
            Epetra_SerialDenseMatrix stiff_relerr;
            stiff_relerr.Shape(6*nnode,6*nnode);

            //characteristic length for numerical approximation of stiffness
            double h_rel = 1e-8;

            //flag indicating whether approximation leads to significant relative error
            int outputflag = 0;

            //calculating strains in new configuration
            for(int i=0; i<6; i++) //for all dof
            {
              for(int k=0; k<nnode; k++)//for all nodes
              {

                Epetra_SerialDenseVector force_aux;
                force_aux.Size(6*nnode);

                //create new displacement and velocity vectors in order to store artificially
         modified displacements vector<double> vel_aux(myvel); vector<double> disp_aux(mydisp);

                //modifying displacement artificially (for numerical derivative of internal forces):
                disp_aux[6*k + i] += h_rel;
                vel_aux[6*k + i] += h_rel / params.get<double>("delta time",0.01);

            //b3_nlnstiffmass is a templated function. therefore we need to point out the number of
         nodes in advance switch(nnode)
                  {
                      case 2:
                      {
                        b3_nlnstiffmass<2>(params,vel_aux,disp_aux,NULL,NULL,&force_aux);
                        break;
                      }
                      case 3:
                      {
                        b3_nlnstiffmass<3>(params,vel_aux,disp_aux,NULL,NULL,&force_aux);
                        break;
                      }
                      case 4:
                      {
                        b3_nlnstiffmass<4>(params,vel_aux,disp_aux,NULL,NULL,&force_aux);
                        break;
                      }
                      case 5:
                      {
                        b3_nlnstiffmass<5>(params,vel_aux,disp_aux,NULL,NULL,&force_aux);
                        break;
                      }
                      case 6:
                      {
                        b3_nlnstiffmass<6>(params,vel_aux,disp_aux,NULL,NULL,&force_aux);
                        break;
                      }
                      default:
                        dserror("Only Line2, Line3, Line4, Line5 and Line6 Elements implemented.");
                  }

                //computing derivative d(fint)/du numerically by finite difference
                for(int u = 0 ; u < 6*nnode ; u++ )
                  stiff_approx(u,k*6+i)= ( pow(force_aux[u],2) - pow(elevec1(u),2) )/ (h_rel *
         (force_aux[u] + elevec1(u) ) );

              } //for(int k=0; k<nnode; k++)//for all nodes

            } //for(int i=0; i<3; i++) //for all dof


            for(int line=0; line<6*nnode; line++)
            {
              for(int col=0; col<6*nnode; col++)
              {
                stiff_relerr(line,col)= fabs( ( pow(elemat1(line,col),2) -
         pow(stiff_approx(line,col),2) )/ ( (elemat1(line,col) + stiff_approx(line,col)) *
         elemat1(line,col) ));

                //suppressing small entries whose effect is only confusing and NaN entires (which
         arise due to zero entries) if ( fabs( stiff_relerr(line,col) ) < h_rel*500 || isnan(
         stiff_relerr(line,col)) || elemat1(line,col) == 0) //isnan = is not a number
                  stiff_relerr(line,col) = 0;

                //if ( stiff_relerr(line,col) > 0)
                  outputflag = 1;
              } //for(int col=0; col<3*nnode; col++)

            } //for(int line=0; line<3*nnode; line++)

            if(outputflag ==1)
            {
              std::cout<<"\n\n acutally calculated stiffness matrix"<< elemat1;
              std::cout<<"\n\n approximated stiffness matrix"<< stiff_approx;
              std::cout<<"\n\n rel error stiffness matrix"<< stiff_relerr;
            }

          } //end of section in which numerical approximation for stiffness matrix is computed
      */
    }
    break;
    case ELEMENTS::struct_calc_update_istep:
    {
      /*the action calc_struct_update_istep is called in the very end of a time step when the new
       * dynamic equilibrium has finally been found; this is the point where the variable
       * representing the geometric status of the beam have to be updated; the geometric status is
       * represented by means of the triads Tnew_, the curvatures curvnew_ and the angular values
       * thetaanew_ and thetaprimenew_*/
      Qconv_ = Qnew_;
      curvconv_ = curvnew_;
      thetaconv_ = thetanew_;
      thetaprimeconv_ = thetaprimenew_;
    }
    break;
    case ELEMENTS::struct_calc_reset_istep:
    {
      /*the action calc_struct_reset_istep is called by the adaptive time step controller; carries
       * out one test step whose purpose is only figuring out a suitabel timestep; thus this step
       * may be a very bad one in order to iterated towards the new dynamic equilibrium and the
       * thereby gained new geometric configuration should not be applied as starting point for any
       * further iteration step; as a consequence the thereby generated change of the geometric
       * configuration should be canceled and the configuration should be reset to the value at the
       * beginning of the time step*/
      Qold_ = Qconv_;
      curvold_ = curvconv_;
      thetaold_ = thetaconv_;
      thetaprimeold_ = thetaprimeconv_;
      Qnew_ = Qconv_;
      curvnew_ = curvconv_;
      thetanew_ = thetaconv_;
      thetaprimenew_ = thetaprimeconv_;
    }
    break;
    case ELEMENTS::struct_calc_stress:
      dserror("No stress output implemented for beam3 elements");
      break;
    case ELEMENTS::struct_calc_brownianforce:
    case ELEMENTS::struct_calc_brownianstiff:
    {
      // number of nodes
      const int nnode = NumNode();

      if (nnode != 2)
        dserror("Only 2 noded elements elements implemented so far for statmech problems.");

      // get element displacements
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);

      // get element velocity
      Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState("velocity");
      if (vel == Teuchos::null) dserror("Cannot get state vectors 'velocity'");
      std::vector<double> myvel(lm.size());
      DRT::UTILS::ExtractMyValues(*vel, myvel, lm);

      if (act == ELEMENTS::struct_calc_brownianforce)
        CalcBrownianForcesAndStiff<2, 3, 6, 4>(params, myvel, mydisp, NULL,
            &elevec1);  // todo template 3,6,4 are hard coded, is that always correct?
      else if (act == ELEMENTS::struct_calc_brownianstiff)
        CalcBrownianForcesAndStiff<2, 3, 6, 4>(params, myvel, mydisp, &elemat1, &elevec1);
      else
        dserror("You shouldn't be here.");

      break;
    }
    case ELEMENTS::struct_calc_recover:
    {
      // do nothing here
      break;
    }
    case ELEMENTS::struct_calc_predict:
    {
      // do nothing here
      break;
    }

    default:
      dserror("Unknown type of action for Beam3 %d", act);
      break;
  }
  return 0;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  Integrate a Line Neumann boundary condition (public)                                       cyron
 03/08|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Beam3::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  SetParamsInterfacePtr(params);

  // find out whether we will use a time curve
  double time = -1.0;
  if (this->IsParamsInterface())
    time = this->ParamsInterfacePtr()->GetTotalTime();
  else
    time = params.get("total time", -1.0);

  // no. of nodes on this element; the following line is only valid for elements with constant
  // number of degrees of freedom per node
  const int numdf = 6;
  const DiscretizationType distype = this->Shape();

  // gaussian points
  const DRT::UTILS::IntegrationPoints1D intpoints(MyGaussRule(NumNode(), gaussunderintegration));

  // declaration of variable in order to store shape function
  Epetra_SerialDenseVector funct(NumNode());

  // get values and switches from the condition

  // onoff is related to the first 6 flags of a line Neumann condition in the input file;
  // value 1 for flag i says that condition is active for i-th degree of freedom
  const std::vector<int>* onoff = condition.Get<std::vector<int>>("onoff");
  // val is related to the 6 "val" fields after the onoff flags of the Neumann condition
  // in the input file; val gives the values of the force as a multiple of the prescribed load curve
  const std::vector<double>* val = condition.Get<std::vector<double>>("val");
  // funct is related to the 6 "funct" fields after the val field of the Neumann condition
  // in the input file; funct gives the number of the function defined in the section FUNCT
  const std::vector<int>* functions = condition.Get<std::vector<int>>("funct");

#ifndef BEAM3DISCRETELINENEUMANN
  // integration loops
  for (int numgp = 0; numgp < intpoints.nquad; ++numgp)
  {
    // integration points in parameter space and weights
    const double xi = intpoints.qxg[numgp][0];
    const double wgt = intpoints.qwgt[numgp];
    // evaluation of shape funcitons at Gauss points
    DRT::UTILS::shape_function_1D(funct, xi, distype);

    // position vector at the gauss point at reference configuration needed for function evaluation
    std::vector<double> X_ref(3, 0.0);
    // calculate coordinates of corresponding Guass point in reference configuration
    for (int node = 0; node < NumNode(); node++)
    {
      for (int dof = 0; dof < 3; dof++)
      {
        X_ref[dof] += funct[node] * Nodes()[node]->X()[dof];
      }
    }

    double fac = 0;
    fac = wgt * jacobi_[numgp];

    // load vector ar
    double ar[numdf];

    // loop the dofs of a node
    for (int dof = 0; dof < numdf; ++dof) ar[dof] = fac * (*onoff)[dof] * (*val)[dof];
    double functionfac = 1.0;
    int functnum = -1;

    // sum up load components
    for (int dof = 0; dof < 6; ++dof)
    {
      if (functions)
        functnum = (*functions)[dof];
      else
        functnum = -1;

      if (functnum > 0)
      {
        // evaluate function at the position of the current node       --> dof here correct?
        functionfac = DRT::Problem::Instance()->Funct(functnum - 1).Evaluate(dof, &X_ref[0], time);
      }
      else
        functionfac = 1.0;

      for (int node = 0; node < NumNode(); ++node)
      {
        elevec1[node * numdf + dof] += funct[node] * ar[dof] * functionfac;
      }
    }
  }  // for (int numgp=0; numgp<intpoints.nquad; ++numgp)
#else
  {
    // integration points in parameter space and weights
    const double xi = BEAM3DISCRETELINENEUMANN;
    // evaluation of shape funcitons at Gauss points
    DRT::UTILS::shape_function_1D(funct, xi, distype);

    // load vector ar
    double ar[numdf];

    // loop the dofs of a node
    for (int dof = 0; dof < numdf; ++dof) ar[dof] = (*onoff)[dof] * (*val)[dof] * curvefac[dof];

    // sum up load components
    for (int dof = 0; dof < 6; ++dof)
    {
      for (int node = 0; node < NumNode(); ++node)
      {
        elevec1[node * numdf + dof] += funct[node] * ar[dof];
      }
    }
  }
#endif
  return 0;
}

/*-----------------------------------------------------------------------------------------------------------*
 | auxiliary functions for dealing with large rotations and nonlinear stiffness cyron 04/08|
 *----------------------------------------------------------------------------------------------------------*/
// computing basis of stiffness matrix of Crisfield, Vol. 2, equation (17.105)
template <int nnode>
inline void DRT::ELEMENTS::Beam3::computestiffbasis(const LINALG::Matrix<3, 3>& Tnew,
    const LINALG::Matrix<3, 3>& Cm, const LINALG::Matrix<3, 3>& Cb, const LINALG::Matrix<3, 3>& S,
    LINALG::Matrix<6 * nnode, 6 * nnode>& stiffmatrix, const LINALG::Matrix<1, nnode>& funct,
    const LINALG::Matrix<1, nnode>& deriv)
{
  // calculating the first matrix of (17.105) directly involves multiplications of large matrices
  // (e.g. with the 6*nnode x 6 - matrix X) application of the definitions in (17.100) allows
  // blockwise evaluation with multiplication and addition of 3x3-matrices only in the follwoing all
  // the blocks are treated separately; their name is related directly to their content: e.g. TCmTt
  // is the product of the 3 matrices T * C_m * T^t (with T and C_m according to (17.74) and (17.76)
  // for the blockwise calculation on which the following steps are based on the relation S^t = -S
  // for spin matrices was applied

  LINALG::Matrix<3, 3> TCmTt;
  LINALG::Matrix<3, 3> TCbTt;
  LINALG::Matrix<3, 3> StTCmTt;
  LINALG::Matrix<3, 3> StTCmTtS;
  LINALG::Matrix<3, 3> TCmTtS;

  // calculating TCmTt and TCbTt
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      TCmTt(i, j) = 0.0;
      TCbTt(i, j) = 0.0;
      for (int k = 0; k < 3; ++k)
      {
        TCmTt(i, j) += Tnew(i, k) * Cm(k, k) * Tnew(j, k);
        TCbTt(i, j) += Tnew(i, k) * Cb(k, k) * Tnew(j, k);
      }
    }
  }

  // calculating StTCmTt
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      StTCmTt(i, j) = 0.0;
      for (int k = 0; k < 3; ++k) StTCmTt(i, j) += S(k, i) * TCmTt(k, j);
    }
  }

  // calculating StTCmTtS
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      StTCmTtS(i, j) = 0.0;
      for (int k = 0; k < 3; ++k) StTCmTtS(i, j) += StTCmTt(i, k) * S(k, j);
    }
  }

  // calculating TCmTtS
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      TCmTtS(i, j) = 0.0;
      for (int k = 0; k < 3; ++k) TCmTtS(i, j) += TCmTt(i, k) * S(k, j);
    }
  }

  // calculating basis of stiffness matrix by means of above blocks
  // the calculations are based on the work by jbueckle
  for (int n = 0; n < nnode; ++n)
  {
    for (int m = 0; m < nnode; ++m)
    {
      for (int i = 0; i < 3; ++i)
      {
        for (int j = 0; j < 3; ++j)
        {
          stiffmatrix(6 * n + i, 6 * m + j) = deriv(n) * deriv(m) * TCmTt(i, j);
          stiffmatrix(6 * n + 3 + i, 6 * m + 3 + j) =
              funct(n) * funct(m) * StTCmTtS(i, j) + deriv(n) * deriv(m) * TCbTt(i, j);
          stiffmatrix(6 * n + 3 + i, 6 * m + j) = funct(n) * deriv(m) * StTCmTt(i, j);
          stiffmatrix(6 * n + i, 6 * m + 3 + j) = deriv(n) * funct(m) * TCmTtS(i, j);
        }
      }
    }
  }

  return;
}  // DRT::ELEMENTS::Beam3::computestiffbasis

/*---------------------------------------------------------------------------*
 |this function performs an update of the rotation (in quaterion form) at the|
 |numgp-th Gauss point by the incremental rotation deltatheta, by means of a |
 |quaternion product and then computes the respective new triad Tnew at the  |
 | Gauss point                                            (public) cyron02/09|
 *---------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Beam3::updatetriad(
    const LINALG::Matrix<3, 1>& deltatheta, LINALG::Matrix<3, 3>& Tnew, const int numgp)
{
  // computing quaternion equivalent to rotation by deltatheta
  LINALG::Matrix<4, 1> Qrot;
  LARGEROTATIONS::angletoquaternion(deltatheta, Qrot);

  // computing quaternion Qnew_ for new configuration of Qold_ for old configuration by means of a
  // quaternion product
  LARGEROTATIONS::quaternionproduct(Qold_[numgp], Qrot, Qnew_[numgp]);

  // normalizing quaternion in order to make sure that it keeps unit absolute values through time
  // stepping
  double abs = Qnew_[numgp].Norm2();
  for (int i = 0; i < 4; i++) Qnew_[numgp](i) = Qnew_[numgp](i) / abs;

  LARGEROTATIONS::quaterniontotriad(Qnew_[numgp], Tnew);

  return;
}  // DRT::ELEMENTS::Beam3::updatetriad

/*----------------------------------------------------------------------------------------*
 |this function performs an update of the rotation (in quaterion form) at the node by     |
 |the incremental rotation deltatheta, by means of a quaternion product and then computes |
 |the respective new triad Tnew at node                          (public) mukherjee 10/14 |
 *----------------------------------------------------------------------------------------*/

inline void DRT::ELEMENTS::Beam3::UpdateNodalTriad(
    std::vector<double>& disp, std::vector<LINALG::Matrix<3, 1>>& Tnew)
{
  LINALG::Matrix<3, 1> thetanew_node;

  for (int node = 0; node < 2; node++)
  {
    /*compute interpolated angle displacemnt at specific Gauss point;
     *angle displacement is taken from discretization*/
    for (int dof = 0; dof < 3; dof++) thetanew_node(dof) += disp[6 * node + 3 + dof];

    ThetaNewNode_[node] = thetanew_node;
    // computing quaternion equivalent to rotation by deltatheta
    LINALG::Matrix<4, 1> Qcurr;
    LARGEROTATIONS::angletoquaternion(ThetaNewNode_[node], Qcurr);


    LINALG::Matrix<3, 3> Triad(true);
    LARGEROTATIONS::quaterniontotriad(Qcurr, Triad);
    for (int i = 0; i < 3; i++) Tnew[node](i) = Triad(0, i);

  }  // for (int node=0; node<nnode; ++node)

  return;
}  // DRT::ELEMENTS::Beam3::UpdateNodalTriad

/*-------------------------------------------------------------------------------------------*
 | Calculate change in angle from reference configuration                     mukherjee 01/15|
 *-------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3::CalcDeltaTheta(std::vector<double>& disp)
{
  // current tangential vector
  std::vector<LINALG::Matrix<3, 1>> tcurrNode;
  std::vector<LINALG::Matrix<3, 1>> trefNode;
  trefNode.resize(2);
  tcurrNode.resize(2);
  // Update nodal triad
  UpdateNodalTriad(disp, tcurrNode);

  for (int node = 0; node < 2; node++)
  {
    LINALG::Matrix<3, 3> Triad(true);
    LARGEROTATIONS::quaterniontotriad(Qref_[node], Triad);
    for (int i = 0; i < 3; i++) trefNode[node](i) = Triad(0, i);
  }

  // Calculate current angle
  double dotprod = 0.0;
  LINALG::Matrix<1, 3> crossprod(true);
  double CosTheta = 0.0;
  double SinTheta = 0.0;

  double norm_t1_curr = tcurrNode[0].Norm2();
  double norm_t2_curr = tcurrNode[1].Norm2();
  for (int j = 0; j < 3; ++j) dotprod += tcurrNode[0](j) * tcurrNode[1](j);

  CosTheta = dotprod / (norm_t1_curr * norm_t2_curr);

  // cross product
  crossprod(0) = tcurrNode[0](1) * tcurrNode[1](2) - tcurrNode[0](2) * tcurrNode[1](1);
  crossprod(1) = tcurrNode[0](2) * tcurrNode[1](0) - tcurrNode[0](0) * tcurrNode[1](2);
  crossprod(2) = tcurrNode[0](0) * tcurrNode[1](1) - tcurrNode[0](1) * tcurrNode[1](0);

  double norm = crossprod.Norm2();
  SinTheta = norm / (norm_t1_curr * norm_t2_curr);

  double ThetaBoundary1 = M_PI / 4;
  double ThetaBoundary2 = 3 * M_PI / 4;
  double thetacurr = 0;
  if (SinTheta >= 0)
  {
    if (CosTheta >= cos(ThetaBoundary1))
      thetacurr = asin(SinTheta);
    else if (CosTheta <= cos(ThetaBoundary2))
      thetacurr = M_PI - asin(SinTheta);
    else
      thetacurr = acos(CosTheta);
  }
  else
    dserror("Angle more than 180 degrees!");


  // Calculate reference angle
  double norm_t1_ref = trefNode[0].Norm2();
  double norm_t2_ref = trefNode[1].Norm2();
  for (int j = 0; j < 3; ++j) dotprod += trefNode[0](j) * trefNode[1](j);

  CosTheta = dotprod / (norm_t1_ref * norm_t2_ref);

  // cross product
  crossprod(0) = trefNode[0](1) * trefNode[1](2) - trefNode[0](2) * trefNode[1](1);
  crossprod(1) = trefNode[0](2) * trefNode[1](0) - trefNode[0](0) * trefNode[1](2);
  crossprod(2) = trefNode[0](0) * trefNode[1](1) - trefNode[0](1) * trefNode[1](0);

  norm = crossprod.Norm2();
  SinTheta = norm / (norm_t1_ref * norm_t2_ref);

  double thetaref = 0;
  if (SinTheta >= 0)
  {
    if (CosTheta >= cos(ThetaBoundary1))
      thetaref = asin(SinTheta);
    else if (CosTheta <= cos(ThetaBoundary2))
      thetaref = M_PI - asin(SinTheta);
    else
      thetaref = acos(CosTheta);
  }
  else
    dserror("Angle more than 180 degrees!");

  deltatheta_ = abs(thetacurr - thetaref);

  return;
}

/*-------------------------------------------------------------------------------------------------------*
 |updating local curvature according to Crisfield, Vol. 2, pages 209 - 210; an exact update of | |
 the curvature is computed by means of equation (16.148) instead of an approximated one as given by
 | | equs. (17.72) and (17.73); note that this function requires as input parameters the angle delta
 theta | | from (16.146), which is responsible for the rotation of the triad at the Gauss point as
 well as its   | | derivative with respect to the curve parameter s, i.e. d(delta theta)/ds. This
 derivative is denoted  | | as deltathetaprime (public)cyron02/09|
 *-------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Beam3::updatecurvature(const LINALG::Matrix<3, 3>& Tnew,
    LINALG::Matrix<3, 1>& deltatheta, LINALG::Matrix<3, 1>& deltathetaprime, const int numgp)
{
  // declaration of local variables:
  LINALG::Matrix<3, 1> omega;       // omega according to (16.146)
  LINALG::Matrix<3, 1> omegaprime;  // omega' according to (16.147)
  LINALG::Matrix<3, 1> chignl;      // chi_gnl according to (16.148)

  // absolute value of rotation vector delta theta
  double abs_theta = deltatheta.Norm2();

  // as in (16.147) division by delta theta arises we have to assume that this angle is unequal to
  // zero
  if (abs_theta > 0)
  {
    // compute omega according to (16.146)
    omega = deltatheta;
    omega.Scale(2 * tan(0.5 * abs_theta) / abs_theta);

    // in (16.147) the brackets contain a 3x3 matrix which is denoted as Aux here and computed now:
    LINALG::Matrix<3, 3> Aux;
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        Aux(i, j) =
            -(1 - abs_theta / sin(abs_theta)) * deltatheta(i) * deltatheta(j) / pow(abs_theta, 2);
        if (i == j) Aux(i, j) += 1;
      }
    }

    omegaprime.Multiply(Aux, deltathetaprime);
    /*we apply the prefactor of (16.147); here we need an angle -PI < theta < PI; note that theta
     * may be assumed to lie within this domain as otherwise the element would rotate at a specific
     * Gauss point by more that 180° within one single iteration step which would disrupt
     * convergence anyway; note that one could ensure -PI < theta < PI easily by applying the three
     * code lines
     *
     *
     *  LINALG::Matrix<4,1> q;
     *  LARGEROTATIONS::angletoquaternion(deltatheta,q);
     *  LARGEROTATIONS::quaterniontoangle(q,deltatheta);
     *
     * in the very beginning of this method. However, for the above reason we assume that theta lies
     * always in the proper region and thus save the related compuational cost and just throw an
     * error if this prerequesite is unexpectedly not satisfied */
    if (abs_theta > M_PI)
      dserror("delta theta exceeds region for which equation (16.147) is valid");
    omegaprime.Scale(2 * tan(abs_theta / 2) / abs_theta);

    // compute chignl from omega and omega' according to (16.148)
    chignl(0) = 0.5 * (omega(1) * omegaprime(2) - omega(2) * omegaprime(1));
    chignl(1) = 0.5 * (omega(2) * omegaprime(0) - omega(0) * omegaprime(2));
    chignl(2) = 0.5 * (omega(0) * omegaprime(1) - omega(1) * omegaprime(0));
    chignl += omegaprime;
    chignl.Scale(1 / (1 + pow(tan(abs_theta / 2), 2)));
  }
  // with delta theta == 0 we get omega == 0 and thus according to (16.147) omega' = d(delta
  // theta)/ds
  else
    chignl = deltathetaprime;

  curvnew_[numgp].MultiplyTN(Tnew, chignl);
  curvnew_[numgp] += curvold_[numgp];

  return;
}  // DRT::ELEMENTS::Beam3::updatecurvature



/*------------------------------------------------------------------------------------------------------*
 |updating local curvature according approximately to Crisfield, Vol. 2, eqs. (17.72) and (17.73);
 note:| |this update is suitable for beams with linear shape functions, only (public)cyron02/09|
 *------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Beam3::approxupdatecurvature(const LINALG::Matrix<3, 3>& Tnew,
    LINALG::Matrix<3, 1> deltatheta, LINALG::Matrix<3, 1> deltathetaprime, const int numgp)
{
  // old triad
  LINALG::Matrix<3, 3> Told;
  LARGEROTATIONS::quaterniontotriad(Qold_[numgp], Told);

  // compute spin matrix from eq. (17.73)
  LINALG::Matrix<3, 3> spin;
  LARGEROTATIONS::computespin(spin, deltatheta);

  // turning spin matrix to left right hand side matrix of eq. (17.73)
  for (int i = 0; i < 3; i++) spin(i, i) += 1;

  // complete right hand side matrix of eq. (17.73)
  // mid point triad
  LINALG::Matrix<3, 3> Tmid;
  Tmid.Multiply(spin, Told);

  // eq. (17.72)
  curvnew_[numgp].MultiplyTN(Tmid, deltathetaprime);

  curvnew_[numgp] += curvold_[numgp];

  return;
}  // DRT::ELEMENTS::Beam3::approxupdatecurvature



/*-----------------------------------------------------------------------------------------------------*
 |computes stiffness matrix Ksigma1 according to Crisfield, Vol. 2, equation (17.106)
 (public)cyron02/09|
 *-----------------------------------------------------------------------------------------------------*/
template <int nnode>
inline void DRT::ELEMENTS::Beam3::computeKsig1(LINALG::Matrix<6 * nnode, 6 * nnode>& Ksig1,
    const LINALG::Matrix<3, 1>& stressn, const LINALG::Matrix<3, 1>& stressm,
    const LINALG::Matrix<1, nnode>& funct, const LINALG::Matrix<1, nnode>& deriv)
{
  LINALG::Matrix<3, 3> Sn;
  LINALG::Matrix<3, 3> Sm;
  // first we caluclate S(n) and S(m)
  LARGEROTATIONS::computespin(Sn, stressn);
  LARGEROTATIONS::computespin(Sm, stressm);

  // then we insert them blockwise in Ksig1
  for (int n = 0; n < nnode; ++n)
  {
    for (int m = 0; m < nnode; ++m)
    {
      for (int i = 0; i < 3; ++i)
      {
        for (int j = 0; j < 3; ++j)
        {
          Ksig1(6 * n + i, 6 * m + j) = 0;
          Ksig1(6 * n + 3 + i, 6 * m + 3 + j) = -deriv(n) * funct(m) * Sm(i, j);
          Ksig1(6 * n + 3 + i, 6 * m + j) = 0;
          Ksig1(6 * n + i, 6 * m + 3 + j) = -deriv(n) * funct(m) * Sn(i, j);
        }
      }
    }
  }
  return;
}  // DRT::ELEMENTS::Beam3::computeKsig1



/*----------------------------------------------------------------------------------------------------------------------*
 |computes stiffness matrix Ksigma2 according to Crisfield, Vol. 2, equation (17.107a) and (17.107b)
 (public)  cyron02/09|
 *----------------------------------------------------------------------------------------------------------------------*/
template <int nnode>
inline void DRT::ELEMENTS::Beam3::computeKsig2(LINALG::Matrix<6 * nnode, 6 * nnode>& Ksig2,
    const LINALG::Matrix<3, 1>& stressn, const LINALG::Matrix<3, 3>& S,
    const LINALG::Matrix<1, nnode>& funct, const LINALG::Matrix<1, nnode>& deriv)
{
  LINALG::Matrix<3, 3> Sn;
  LINALG::Matrix<3, 3> Y;
  // first we compute S(n) and Y according to (17.107b)
  LARGEROTATIONS::computespin(Sn, stressn);
  Y.Multiply(S, Sn);

  // then we insert them into Ksig2 by means of a blockwise algorithm
  for (int n = 0; n < nnode; ++n)
  {
    for (int m = 0; m < nnode; ++m)
    {
      for (int i = 0; i < 3; ++i)
      {
        for (int j = 0; j < 3; ++j)
        {
          Ksig2(6 * n + i, 6 * m + j) = 0;
          Ksig2(6 * n + 3 + i, 6 * m + 3 + j) = funct(n) * funct(m) * Y(i, j);
          Ksig2(6 * n + 3 + i, 6 * m + j) = funct(n) * deriv(m) * Sn(i, j);
          Ksig2(6 * n + i, 6 * m + 3 + j) = 0;
        }
      }
    }
  }
  return;
}  // DRT::ELEMENTS::Beam3::computeKsig2

/*------------------------------------------------------------------------------------------------------------*
 | calculation of elastic energy (private) cyron 12/10|
 *-----------------------------------------------------------------------------------------------------------*/
template <int nnode>
void DRT::ELEMENTS::Beam3::b3_energy(
    Teuchos::ParameterList& params, std::vector<double>& disp, Epetra_SerialDenseVector* intenergy)
{
  // initialize energies (only one kind of energy computed here
  intenergy->Scale(0.0);
  bool calcenergy = false;
  if (params.isParameter("energyoftype") == false)
    calcenergy = true;
  else if (params.get<std::string>("energyoftype") == "beam3")
    calcenergy = true;

  if (calcenergy)
  {
    // constitutive matrices from Crisfield, Vol. 2, equation (17.76)
    LINALG::Matrix<3, 3> Cm;
    LINALG::Matrix<3, 3> Cb;

    // normal/shear strain and bending strain(curvature)
    LINALG::Matrix<3, 1> epsilonn;
    LINALG::Matrix<3, 1> epsilonm;

    // derivative of x with respect to xi
    LINALG::Matrix<3, 1> dxdxi_gp;

    // triad at GP, Crisfiel Vol. 2, equation (17.73)
    LINALG::Matrix<3, 3> Tnew;

    /*first displacement vector is modified for proper element evaluation in case of periodic
     *boundary conditions; in case that no periodic boundary conditions are to be applied the
     *following code line may be ignored or deleted*/

    // Get integrationpoints for underintegration
    DRT::UTILS::IntegrationPoints1D gausspoints(MyGaussRule(nnode, gaussunderintegration));

    // Get DiscretizationType
    const DRT::Element::DiscretizationType distype = Shape();

    // Matrices for h and h,xi
    LINALG::Matrix<1, nnode> deriv;

    // Loop through all GP and calculate their contribution to the forcevector and stiffnessmatrix
    for (int numgp = 0; numgp < gausspoints.nquad; numgp++)
    {
      // Get location and weight of GP in parameter space
      const double xi = gausspoints.qxg[numgp][0];
      const double wgt = gausspoints.qwgt[numgp];

      // Get h,xi
      DRT::UTILS::shape_function_1D_deriv1(deriv, xi, distype);

      // set up current dxdxi, theta, and thetaprime at the GP
      dxdxi_gp.Clear();
      for (int dof = 0; dof < 3; ++dof)  // j
        for (int node = 0; node < nnode; ++node)
          dxdxi_gp(dof) += (Nodes()[node]->X()[dof] + disp[6 * node + dof]) * deriv(node);


      // compute current triad at numgp-th Gauss point
      LARGEROTATIONS::quaterniontotriad(Qnew_[numgp], Tnew);

      // setting constitutive parameters , Crisfield, Vol. 2, equation (17.76)
      GetConstitutiveMatrices(Cm, Cb);

      // computing current axial and shear strain epsilon, Crisfield, Vol. 2, equation (17.97)
      epsilonn.Clear();
      epsilonn.MultiplyTN(Tnew, dxdxi_gp);
      epsilonn.Scale(1 / jacobi_[numgp]);
      epsilonn(0) -= 1.0;

      epsilonm.Clear();
      epsilonm = curvnew_[numgp];

      // adding elastic energy from epsilonn at this Gauss point
      if (intenergy->M() == 1)
      {
        for (int i = 0; i < 3; i++)
        {
          (*intenergy)(0) += 0.5 * epsilonn(i) * epsilonn(i) * Cm(i, i) * wgt * jacobi_[numgp];
          (*intenergy)(0) += 0.5 * epsilonm(i) * epsilonm(i) * Cb(i, i) * wgt * jacobi_[numgp];
        }
      }
      else if (intenergy->M() == 6)
      {
        for (int i = 0; i < 3; i++)
        {
          (*intenergy)(i) += 0.5 * epsilonn(i) * epsilonn(i) * Cm(i, i) * wgt * jacobi_[numgp];
          (*intenergy)(i + 3) += 0.5 * epsilonm(i) * epsilonm(i) * Cb(i, i) * wgt * jacobi_[numgp];
        }
      }
      else
        dserror("energy vector of invalid size %i!", intenergy->M());
    }
  }
  return;
}  // DRT::ELEMENTS::Beam3::b3_energy


/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private) cyron 01/08|
 *-----------------------------------------------------------------------------------------------------------*/
template <int nnode>
void DRT::ELEMENTS::Beam3::b3_nlnstiffmass(Teuchos::ParameterList& params,
    std::vector<double>& disp, Epetra_SerialDenseMatrix* stiffmatrix,
    Epetra_SerialDenseMatrix* massmatrix, Epetra_SerialDenseVector* force)
{
  // Calculate change in angle from reference configuration
  //  CalcDeltaTheta(disp);

  // current node position (first entries 0 .. 2 for first node, 3 ..5 for second node)
  LINALG::Matrix<6, 1> xcurr;

  // auxiliary vector for both internal force and stiffness matrix: N^T_(,xi)*N_(,xi)*xcurr
  LINALG::Matrix<3, 1> aux;

  // current nodal position (first
  for (int j = 0; j < 3; ++j)
  {
    xcurr(j) = Nodes()[0]->X()[j] + disp[j];          // first node
    xcurr(j + 3) = Nodes()[1]->X()[j] + disp[6 + j];  // second node
  }

  // computing auxiliary vector aux = 4.0*N^T_{,xi} * N_{,xi} * xcurr
  aux(0) = (xcurr(0) - xcurr(3));
  aux(1) = (xcurr(1) - xcurr(4));
  aux(2) = (xcurr(2) - xcurr(5));

  // current length
  lcurr_ = sqrt(pow(aux(0), 2) + pow(aux(1), 2) + pow(aux(2), 2));

  // constitutive laws from Crisfield, Vol. 2, equation (17.76)
  LINALG::Matrix<3, 3> Cm;
  LINALG::Matrix<3, 3> Cb;

  // normal/shear strain and bending strain(curvature)
  LINALG::Matrix<3, 1> epsilonn;
  LINALG::Matrix<3, 1> epsilonm;

  // stress values n and m, Crisfield, Vol. 2, equation (17.78)
  LINALG::Matrix<3, 1> stressn;
  LINALG::Matrix<3, 1> stressm;

  // derivative of x with respect to xi
  LINALG::Matrix<3, 1> dxdxi_gp;

  // theta interpolated at the GP. Note: theta is not the current angle theta but only the
  // difference to the reference configuration
  LINALG::Matrix<3, 1> thetanew_gp;


  // derivative of theta with respect to xi
  LINALG::Matrix<3, 1> thetaprimenew_gp;

  // stiffness base of stiffness matrix, Crisfield Vol. 2, equation (17.105)
  LINALG::Matrix<6 * nnode, 6 * nnode> Kstiff_gp;

  // nonlinear parts of stiffness matrix, Crisfield Vol. 2, equation (17.83) and (17.87)
  LINALG::Matrix<6 * nnode, 6 * nnode> Ksig1_gp;
  LINALG::Matrix<6 * nnode, 6 * nnode> Ksig2_gp;


  // triad at GP, Crisfiel Vol. 2, equation (17.73)
  LINALG::Matrix<3, 3> Tnew;

  // assignment of material parameters; only Reissner beam material based on hyperelastic
  // stored energy function is accepted for this beam (analog to St. Venant-Kirchhoff for 3D
  // continua)
  GetConstitutiveMatrices(Cm, Cb);

  //"new" variables have to be adopted to current discplacement

  // Get integrationpoints for underintegration
  DRT::UTILS::IntegrationPoints1D gausspoints(MyGaussRule(nnode, gaussunderintegration));

  // Get DiscretizationType
  const DRT::Element::DiscretizationType distype = Shape();

  // Matrices for h and h,xi
  LINALG::Matrix<1, nnode> funct;
  LINALG::Matrix<1, nnode> deriv;



  // Loop through all GP and calculate their contribution to the forcevector and stiffnessmatrix
  for (int numgp = 0; numgp < gausspoints.nquad; numgp++)
  {
    // Get location and weight of GP in parameter space
    const double xi = gausspoints.qxg[numgp][0];
    const double wgt = gausspoints.qwgt[numgp];

    // Get h and h,xi
    DRT::UTILS::shape_function_1D(funct, xi, distype);
    DRT::UTILS::shape_function_1D_deriv1(deriv, xi, distype);

    dxdxi_gp.Clear();
    thetanew_gp.Clear();
    thetaprimenew_gp.Clear();

    // set up current dxdxi, theta, and thetaprime at the GP
    for (int dof = 0; dof < 3; ++dof)  // j
    {
      for (int node = 0; node < nnode; ++node)
      {
        dxdxi_gp(dof) += (Nodes()[node]->X()[dof] + disp[6 * node + dof]) * deriv(node);

        /*compute interpolated angle displacemnt at specific Gauss point; angle displacement is
         * taken from discretization and interpolated with the according shapefunctions*/
        thetanew_gp(dof) += disp[6 * node + 3 + dof] * funct(node);

        /*compute interpolated angle displacemnt at specific Gauss point;
         *angle displacement is taken from discretization*/
        //        thetanew_node(dof)        += disp[6*node+3+dof];

        /*compute derivative of interpolated angle displacement with respect to curve parameter at
         * specific Gauss point; angle displacement is taken from discretization and interpolated
         * with the according shapefunctions*/
        thetaprimenew_gp(dof) += disp[6 * node + 3 + dof] * deriv(node) / jacobi_[numgp];

      }  // for (int node=0; node<nnode; ++node)
    }    // for (int dof=0; dof<3; ++dof)

    // update theta and thetaprime
    thetanew_[numgp] = thetanew_gp;
    thetaprimenew_[numgp] = thetaprimenew_gp;



    /*to perform a curvature update at a specific Gauss point we need to know the angle deltatheta
     * by which the triad at that Gauss point is rotated and furthermore the derivative of this
     * rotation angle along the curve. The latter one is denoted as d(delta theta)/ds in Crisfield,
     * Vol. 2, page 209, and can be computed according to the comments in the bottom of this page*/

    // compute delta theta for (16.146) at the Gauss point only according to remark in the bottom of
    // page 209
    LINALG::Matrix<3, 1> deltatheta_gp = thetanew_[numgp];
    deltatheta_gp -= thetaold_[numgp];

    // compute d(delta theta)/ds for (16.147) at the Gauss point only according to remark in the
    // bottom of page 209
    LINALG::Matrix<3, 1> deltathetaprime_gp = thetaprimenew_[numgp];
    deltathetaprime_gp -= thetaprimeold_[numgp];

    /*update triad at Gauss point as shown in Crisfield, Vol. 2, equation (17.65) Note: instead of
     *the matrix multiplication of (17.65) we use a quaternion product*/
    updatetriad(deltatheta_gp, Tnew, numgp);

    // updating local curvature according
    // updatecurvature(Tnew,deltatheta_gp,deltathetaprime_gp,numgp);
    approxupdatecurvature(Tnew, deltatheta_gp, deltathetaprime_gp, numgp);

    epsilonn.Clear();

    // computing current axial and shear strain epsilon, Crisfield, Vol. 2, equation (17.97)
    epsilonn.MultiplyTN(Tnew, dxdxi_gp);
    epsilonn.Scale(1 / jacobi_[numgp]);
    epsilonn(0) -= 1.0;

    epsilonm.Clear();

    // computing spin matrix S(dxdxi_gp) according to Crisfield, Vol. 2, equation (17.100)
    LINALG::Matrix<3, 3> S_gp;

    S_gp.Clear();

    LARGEROTATIONS::computespin(S_gp, dxdxi_gp);

    // stress values n and m, Crisfield, Vol. 2, equation (17.76) and (17.78)
    epsilonn(0) *= Cm(0, 0);
    epsilonn(1) *= Cm(1, 1);
    epsilonn(2) *= Cm(2, 2);

    stressn.Clear();

    stressn.Multiply(Tnew, epsilonn);

    // turning bending strain epsilonm into bending stress stressm
    epsilonm = curvnew_[numgp];

    epsilonm(0) *= Cb(0, 0);
    epsilonm(1) *= Cb(1, 1);
    epsilonm(2) *= Cb(2, 2);

    stressm.Clear();

    stressm.Multiply(Tnew, epsilonm);

    // computing global internal forces, Crisfield Vol. 2, equation (17.102)
    // note: X = [h,xi(1)*I 0; h(1)S h,xi(1)I;h,xi(2)*I 0; h(2)S h,xi(2)I......] with S =
    // S(dxdx_gp);
    if (force != NULL)
    {
      for (int node = 0; node < nnode; ++node)
      {
        for (int j = 0; j < 3; ++j)
        {
          (*force)(6 * node + j) += wgt * deriv(node) * stressn(j);
          (*force)(6 * node + 3 + j) += wgt * deriv(node) * stressm(j);

          for (int i = 0; i < 3; ++i)
            (*force)(6 * node + 3 + j) += wgt * funct(node) * S_gp(i, j) * stressn(i);
        }
      }
    }  // if (force != NULL)


    // computing nonlinear stiffness matrix
    if (stiffmatrix != NULL)
    {
      Kstiff_gp.Clear();
      Ksig1_gp.Clear();
      Ksig2_gp.Clear();


      // setting up basis of stiffness matrix according to Crisfield, Vol. 2, equation (17.105)
      computestiffbasis<nnode>(Tnew, Cm, Cb, S_gp, Kstiff_gp, funct, deriv);

      Kstiff_gp.Scale(wgt / jacobi_[numgp]);

      // adding nonlinear (stress dependent) parts to tangent stiffness matrix, Crisfield, Vol. 2
      // equs. (17.106), (17.107) and (17.107b)
      computeKsig1<nnode>(Ksig1_gp, stressn, stressm, funct, deriv);
      computeKsig2<nnode>(Ksig2_gp, stressn, S_gp, funct, deriv);

      Ksig1_gp.Scale(wgt);
      Ksig2_gp.Scale(wgt);

      // shifting values from fixed size matrix to epetra matrix *stiffmatrix
      for (int i = 0; i < 6 * nnode; i++)
      {
        for (int j = 0; j < 6 * nnode; j++)
        {
          (*stiffmatrix)(i, j) += Kstiff_gp(i, j);
          (*stiffmatrix)(i, j) += Ksig1_gp(i, j);
          (*stiffmatrix)(i, j) += Ksig2_gp(i, j);
        }
      }
    }  // if (stiffmatrix != NULL)
  }    // for(int numgp=0; numgp < gausspoints.nquad; numgp++)


  // calculating mass matrix (local version = global version)
  // We use a consistent Timoshenko mass matrix here
  if (massmatrix != NULL)
  {
    // tensor of mass moments of inertia for translational and rotational motion
    double mass_inertia_translational = 0.0;
    LINALG::Matrix<3, 3> Jp(true);

    GetTranslationalAndRotationalMassInertiaTensor(mass_inertia_translational, Jp);

    // Get the applied integrationpoints for complete integration
    DRT::UTILS::IntegrationPoints1D gausspointsmass(MyGaussRule(nnode, gaussexactintegration));

    // Matrix to store Shape functions as defined throughout the FE lecture
    LINALG::Matrix<6, 6 * nnode> N;

    // Matrix to store mass matrix at each GP in
    LINALG::Matrix<6 * nnode, 6 * nnode> massmatrixgp;

    for (int numgp = 0; numgp < gausspointsmass.nquad; numgp++)
    {
      // Get location and weight of GP in parameter space
      const double xi = gausspointsmass.qxg[numgp][0];
      const double wgt = gausspointsmass.qwgt[numgp];

      // Get h
      DRT::UTILS::shape_function_1D(funct, xi, distype);

      // Set N to zeros
      N.Clear();

      // Fill up N as defined in the FE lecture
      for (int node = 0; node < nnode; node++)
        for (int dof = 0; dof < 6; dof++) N(dof, 6 * node + dof) = funct(node);

      // m= integraloverxi [N_t*N]
      massmatrixgp.MultiplyTN(1.0, N, N);

      std::cout << "*** WARNING ***: Incorrect mass matrix implementation. This beam element is "
                   "only applicable to static problems so far!"
                << std::endl;
      // According to the paper of Jelenic and Crisfield "Geometrically exact 3D beam theory:
      // implementation of a strain-invariant finite element for statics and dynamics, 1999, page
      // 146, a time integration scheme that delivers angular velocities and angular accelerations
      // as needed for the inertia terms of geometrically exact beams has to be based on
      // multiplicative rotation angle increments between two successive time steps. Since BACI does
      // all displacement updates in an additive manner, the global vector of rotational
      // displacements has no physical meaning and, consequently the global velocity and
      // acceleration vectors resulting from the BACI time integration schemes have no physical
      // meaning, too. Therefore, a mass matrix in combination with this global acceleration vector
      // is meaningless from a physical point of view. (Christoph Meier, 04.14)

      for (int i = 0; i < 6 * nnode; i++)
      {
        for (int j = 0; j < nnode; j++)
        {
          // These lines multiply all entries associated with translational motion with the
          // corresponging pre-factor (density*crosssec)
          massmatrixgp(i, 6 * j + 0) = mass_inertia_translational * massmatrixgp(i, 6 * j + 0);
          massmatrixgp(i, 6 * j + 1) = mass_inertia_translational * massmatrixgp(i, 6 * j + 1);
          massmatrixgp(i, 6 * j + 2) = mass_inertia_translational * massmatrixgp(i, 6 * j + 2);

          // These lines multiply all entries associated with
          // the rotations. The Irr_ comes from calculations considering
          // the moment created by a rotation theta and refers to the assumed direction
          // Note: This is an approximation for the rotational values around t2 and t3. For exact
          // one would have to differentiate again.

          /* update:
           * I had to adapt this to the more general material parameter definition from input
           * file. The rotational mass inertia factor which was used before is now extracted as
           * the axial component of the mass moment of inertia tensor */
          massmatrixgp(i, 6 * j + 3) = Jp(0, 0) * massmatrixgp(i, 6 * j + 3);
          massmatrixgp(i, 6 * j + 4) = Jp(0, 0) * massmatrixgp(i, 6 * j + 4);
          massmatrixgp(i, 6 * j + 5) = Jp(0, 0) * massmatrixgp(i, 6 * j + 5);
        }
      }

      // Sum up the massmatrices calculated at each gp using the lengthfactor and the weight
      for (int i = 0; i < 6 * nnode; i++)
        for (int j = 0; j < 6 * nnode; j++)
          (*massmatrix)(i, j) += jacobimass_[numgp] * wgt * massmatrixgp(i, j);


    }  // for(int numgp=0; numgp < gausspointsmass.nquad; numgp++)

  }  // if (massmatrix != NULL)

  // in statistical mechanics simulations, a deletion influenced by the values of the internal force
  // vector might occur
  if (params.get<std::string>("internalforces", "no") == "yes" && force != NULL)
  {
    eps_ = epsilonn(0);
    f_ = *force;
    Ngp_ = stressn;  // correct?
  }

  LINALG::Matrix<6, 1> axialforce;
  axialforce.Clear();
  LINALG::Matrix<6, 1> moment;
  moment.Clear();
  for (int i = 0; i < 3; i++)
  {
    axialforce(i) = (*force)(i);
    axialforce(i + 3) = (*force)(i + 6);
    moment(i) = (*force)(i + 3);
    moment(i + 3) = (*force)(i + 9);
  }

  NormForce = axialforce.Norm2();
  NormMoment = moment.Norm2();
  if (NormMoment != 0)
    RatioNormForceMoment = NormForce / NormMoment;
  else
    RatioNormForceMoment = 0;

  //  std::cout<<"Beam3 Element ID="<<this->Id()<<std::endl;
  //  std::cout<<"Norm Force="<<NormForce<<std::endl;
  //  std::cout<<"Norm Moment="<<NormMoment<<std::endl;
  //  std::cout<<"Norm force/Norm Moment="<<axialforce.Norm2()/moment.Norm2()<<std::endl;

  //  std::cout<<"beam3 force="<<*force<<std::endl;

  return;

}  // DRT::ELEMENTS::Beam3::b3_nlnstiffmass

/*------------------------------------------------------------------------------------------------------------*
 | lump mass matrix             (private)                                                   cyron
 01/08|
 *------------------------------------------------------------------------------------------------------------*/
template <int nnode>
void DRT::ELEMENTS::Beam3::lumpmass(Epetra_SerialDenseMatrix* emass)
{
  // lump mass matrix
  if (emass != NULL)
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
 | Evaluate PTC damping (public) cyron 10/08|
 *----------------------------------------------------------------------------------------------------------*/
template <int nnode>
void DRT::ELEMENTS::Beam3::EvaluatePTC(
    Teuchos::ParameterList& params, Epetra_SerialDenseMatrix& elemat1)
{
  // Get the applied integrationpoints for underintegration
  DRT::UTILS::IntegrationPoints1D gausspointsptc(MyGaussRule(nnode, gaussunderintegration));
  // Get discretization typ
  const DRT::Element::DiscretizationType distype = Shape();
  // matrix to store Ansatzfunktionen
  LINALG::Matrix<1, nnode> funct;

  for (int gp = 0; gp < gausspointsptc.nquad; gp++)
  {
    // Get location and weight of GP in parameter space
    const double xi = gausspointsptc.qxg[gp][0];
    const double wgt = gausspointsptc.qwgt[gp];

    DRT::UTILS::shape_function_1D(funct, xi, distype);

    // computing angle increment from current position in comparison with last converged position
    // for damping
    LINALG::Matrix<4, 1> deltaQ;
    LARGEROTATIONS::quaternionproduct(
        LARGEROTATIONS::inversequaternion(Qconv_[gp]), Qnew_[gp], deltaQ);
    LINALG::Matrix<3, 1> deltatheta;
    LARGEROTATIONS::quaterniontoangle(deltaQ, deltatheta);

    // isotropic artificial stiffness
    LINALG::Matrix<3, 3> artstiff;
    artstiff = LARGEROTATIONS::Tmatrix(deltatheta);

    // scale artificial damping with crotptc parameter for PTC method
    artstiff.Scale(params.get<double>("crotptc", 0.0));

    for (int i = 0; i < nnode; i++)
      for (int j = 0; j < nnode; j++)
        for (int k = 0; k < 3; k++)
          for (int l = 0; l < 3; l++)
            elemat1(i * 6 + 3 + k, j * 6 + 3 + l) +=
                artstiff(k, l) * funct(i) * funct(j) * wgt * jacobi_[gp];

    // PTC for translational degrees of freedom
    for (int i = 0; i < nnode; i++)
      for (int j = 0; j < nnode; j++)
        for (int k = 0; k < 3; k++)
          elemat1(i * 6 + k, j * 6 + k) +=
              params.get<double>("ctransptc", 0.0) * funct(i) * funct(j) * wgt * jacobi_[gp];
  }

  return;
}  // DRT::ELEMENTS::Beam3::EvaluatePTC

/*-----------------------------------------------------------------------------------------------------------*
 | computes damping coefficients per lengthand stores them in a matrix in the following order:
 damping of    | | translation parallel to filament axis, damping of translation orthogonal to
 filament axis, damping of     | | rotation around filament axis (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Beam3::MyDampingConstants(LINALG::Matrix<3, 1>& gamma)
{
  // translational damping coefficients according to Howard, p. 107, table 6.2;
  gamma(0) = 2 * PI * BrownianDynParamsInterface().GetViscosity();
  gamma(1) = 4 * PI * BrownianDynParamsInterface().GetViscosity();

  /*damping coefficient of rigid straight rod spinning around its own axis according to Howard, p.
   *107, table 6.2;
   *as this coefficient is very small for thin rods it is increased artificially by a factor for
   *numerical convencience*/
  double rsquare = std::pow(GetCircularCrossSectionRadiusForInteractions(), 2.0);
  double artificial =
      4000;  // 1000;  //1000 not bad for standard Actin3D_10.dat files; for 40 elements also 1
             // seems to work really well; for large networks 4000 seems good (artificial
             // contribution then still just ~0.1 % of nodal moments)
  gamma(2) = 4 * PI * BrownianDynParamsInterface().GetViscosity() * rsquare * artificial;

  /* in the following section damping coefficients are replaced by those suggested in Ortega2003,
   * which allows for a comparison of the finite element simulation with the results of that
   * article; note that we assume that the element length is equivalent to the particle length in
   * the following when computing the length to diameter ratio p*/
  /*
  double lrefe=0;
  for (int gp=0; gp<nnode-1; gp++)
    lrefe += gausspointsdamping.qwgt[gp]*jacobi_[gp];

  double p=lrefe/(pow(crosssec_*4.0/PI,0.5));
  double Ct=0.312+0.565/p-0.100/pow(p,2.0);
  double Cr=-0.662+0.917/p-0.05/pow(p,2.0);
  gamma(0) = 2.0*PI*BrownianDynParamsInterface().GetViscosity()/(log(p) + 2*Ct - Cr);
  gamma(1) = 4.0*PI*BrownianDynParamsInterface().GetViscosity()/(log(p)+Cr);
  gamma(3) = 4.0*PI*BrownianDynParamsInterface().GetViscosity()*rsquare*artificial*(0.96 + 0.64992/p
  - 0.17568/p^2);
  */
}

/*-----------------------------------------------------------------------------------------------------------*
 |computes the number of different random numbers required in each time step for generation of
 stochastic    | |forces; (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3::HowManyRandomNumbersINeed() const
{
  /*at each Gauss point one needs as many random numbers as randomly excited degrees of freedom,
   *i.e. three
   *random numbers for the translational degrees of freedom and one random number for the rotation
   *around the element axis*/
  return (4 * NumNode());
}

/*-----------------------------------------------------------------------------------------------------------*
 |computes velocity of background fluid and gradient of that velocity at a certain evaluation point
 in       | |the physical space                                                         (public)
 cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template <int ndim>  // number of dimensions of embedding space
void DRT::ELEMENTS::Beam3::MyBackgroundVelocity(Teuchos::ParameterList& params,  //!< parameter list
    const LINALG::Matrix<ndim, 1>& evaluationpoint,  //!< point at which background velocity and its
                                                     //!< gradient has to be computed
    LINALG::Matrix<ndim, 1>& velbackground,          //!< velocity of background fluid
    LINALG::Matrix<ndim, ndim>& velbackgroundgrad)   //!< gradient of velocity of background fluid
{
  /*note: this function is not yet a general one, but always assumes a shear flow, where the
   * velocity of the background fluid is always directed in direction
   * params.get<int>("DBCDISPDIR",0) and orthogonal to z-axis. In 3D the velocity increases linearly
   * in z and equals zero for z = 0.
   * In 2D the velocity increases linearly in y and equals zero for y = 0. */

  // velocity at upper boundary of domain
  //  double uppervel = 0.0;

  // default values for background velocity and its gradient
  velbackground.PutScalar(0);
  velbackgroundgrad.PutScalar(0);

  // fixme: this needs to go somewhere else, outside of element
  //  double time = -1.0;
  //
  //
  //  double shearamplitude = BrownianDynParamsInterface().GetShearAmplitude();
  //  int curvenumber = BrownianDynParamsInterface().GetCurveNumber() - 1;
  //  int dbcdispdir = BrownianDynParamsInterface().GetDbcDispDir() - 1;
  //
  //  Teuchos::RCP<std::vector<double> > periodlength =
  //  BrownianDynParamsInterface().GetPeriodLength(); INPAR::STATMECH::DBCType dbctype =
  //  BrownianDynParamsInterface().GetDbcType(); bool shearflow = false;
  //  if(dbctype==INPAR::STATMECH::dbctype_shearfixed ||
  //     dbctype==INPAR::STATMECH::dbctype_shearfixeddel ||
  //     dbctype==INPAR::STATMECH::dbctype_sheartrans ||
  //     dbctype==INPAR::STATMECH::dbctype_affineshear||
  //     dbctype==INPAR::STATMECH::dbctype_affinesheardel)
  //    shearflow = true;
  //
  //  if(periodlength->at(0) > 0.0)
  //    if(shearflow  && curvenumber >=  0 && dbcdispdir >= 0 )
  //    {
  //      uppervel = shearamplitude *
  //      (DRT::Problem::Instance()->Funct(curvenumber).EvaluateTimeDerivative(time,1))[1];
  //
  //      //compute background velocity
  //      velbackground(dbcdispdir) = (evaluationpoint(ndim-1) / periodlength->at(ndim-1)) *
  //      uppervel;
  //
  //      //compute gradient of background velocity
  //      velbackgroundgrad(dbcdispdir,ndim-1) = uppervel / periodlength->at(ndim-1);
  //    }
}
/*-----------------------------------------------------------------------------------------------------------*
 | computes rotational damping forces and stiffness (public) cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template <int nnode>  // number of nodes
inline void DRT::ELEMENTS::Beam3::MyRotationalDamping(
    Teuchos::ParameterList& params,         //!< parameter list
    const std::vector<double>& vel,         //!< element velocity vector
    const std::vector<double>& disp,        //!< element disp vector
    Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
    Epetra_SerialDenseVector* force)        //!< element internal force vector
{
  // get time step size
  double dt = 1000;
  if (IsParamsInterface())
    dt = ParamsInterface().GetDeltaTime();
  else
    dt = params.get<double>("delta time", 1000);

  // integration points for underintegration
  DRT::UTILS::IntegrationPoints1D gausspoints(MyGaussRule(nnode, gaussunderintegration));

  // damping coefficients for translational and rotatinal degrees of freedom
  LINALG::Matrix<3, 1> gamma(true);
  MyDampingConstants(gamma);


  // matrix to store basis functions evaluated at a certain Gauss point
  LINALG::Matrix<1, nnode> funct;

  for (int gp = 0; gp < gausspoints.nquad; gp++)  // loop through Gauss points
  {
    // get evaluated basis functions at current Gauss point
    DRT::UTILS::shape_function_1D(funct, gausspoints.qxg[gp][0], Shape());

    // rotation between last converged position and current position expressend as a quaternion
    LINALG::Matrix<4, 1> deltaQ;
    LARGEROTATIONS::quaternionproduct(
        LARGEROTATIONS::inversequaternion(Qconv_[gp]), Qnew_[gp], deltaQ);

    // rotation between last converged position and current position expressed as a three element
    // rotation vector
    LINALG::Matrix<3, 1> deltatheta;
    LARGEROTATIONS::quaterniontoangle(deltaQ, deltatheta);

    // angular velocity at this Gauss point
    LINALG::Matrix<3, 1> omega(true);
    omega += deltatheta;
    omega.Scale(1 / dt);

    // compute matrix T*W*T^t
    LINALG::Matrix<3, 3> Tnew;
    LINALG::Matrix<3, 3> TWTt;
    LARGEROTATIONS::quaterniontotriad(Qnew_[gp], Tnew);
    for (int k = 0; k < 3; k++)
      for (int j = 0; j < 3; j++) TWTt(k, j) = Tnew(k, 0) * Tnew(j, 0);

    // compute vector T*W*T^t*\omega
    LINALG::Matrix<3, 1> TWTtomega;
    TWTtomega.Multiply(TWTt, omega);

    // compute matrix T*W*T^t*H^(-1)
    LINALG::Matrix<3, 3> TWTtHinv;
    TWTtHinv.Multiply(TWTt, LARGEROTATIONS::Tmatrix(deltatheta));

    // compute spin matrix S(\omega)
    LINALG::Matrix<3, 3> Sofomega;
    LARGEROTATIONS::computespin(Sofomega, omega);

    // compute matrix T*W*T^t*S(\omega)
    LINALG::Matrix<3, 3> TWTtSofomega;
    TWTtSofomega.Multiply(TWTt, Sofomega);

    // compute spin matrix S(T*W*T^t*\omega)
    LINALG::Matrix<3, 3> SofTWTtomega;
    LARGEROTATIONS::computespin(SofTWTtomega, TWTtomega);

    // loop over all line nodes
    for (int i = 0; i < nnode; i++)
      // loop over three dimensions in line direction
      for (int k = 0; k < 3; k++)
      {
        if (force != NULL)
          (*force)(i * 6 + 3 + k) +=
              gamma(2) * TWTtomega(k) * funct(i) * gausspoints.qwgt[gp] * jacobi_[gp];

        if (stiffmatrix != NULL)
          // loop over all column nodes
          for (int j = 0; j < nnode; j++)
            // loop over three dimensions in column direction
            for (int l = 0; l < 3; l++)
              (*stiffmatrix)(i * 6 + 3 + k, j * 6 + 3 + l) +=
                  gamma(2) * (TWTtHinv(k, l) / dt + TWTtSofomega(k, l) - SofTWTtomega(k, l)) *
                  funct(i) * funct(j) * gausspoints.qwgt[gp] * jacobi_[gp];
      }
  }


  return;
}  // DRT::ELEMENTS::Beam3::MyRotationalDamping(.)

/*-----------------------------------------------------------------------------------------------------------*
 | computes translational damping forces and stiffness (public) cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template <int nnode, int ndim,
    int dof>  // number of nodes, number of dimensions of embedding space, number of degrees of
              // freedom per node
inline void
DRT::ELEMENTS::Beam3::MyTranslationalDamping(Teuchos::ParameterList& params,  //!< parameter list
    const std::vector<double>& vel,         //!< element velocity vector
    const std::vector<double>& disp,        //!< element disp vector
    Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
    Epetra_SerialDenseVector* force)        //!< element internal force vector
{
  // get time step size
  double dt = 1000;
  if (IsParamsInterface())
    dt = ParamsInterface().GetDeltaTime();
  else
    dt = params.get<double>("delta time", 1000);

  // velocity and gradient of background velocity field
  LINALG::Matrix<ndim, 1> velbackground;
  LINALG::Matrix<ndim, ndim> velbackgroundgrad;

  // evaluation point in physical space corresponding to a certain Gauss point in parameter space
  LINALG::Matrix<ndim, 1> evaluationpoint;

  // damping coefficients for translational and rotatinal degrees of freedom
  LINALG::Matrix<3, 1> gamma(true);
  MyDampingConstants(gamma);

  // get vector jacobi with Jacobi determinants at each integration point (gets by default those
  // values required for consistent damping matrix)
  std::vector<double> jacobi(jacobimass_);

  // determine type of numerical integration performed (lumped damping matrix via lobatto
  // integration!)
  IntegrationType integrationtype = gaussexactintegration;

  // get Gauss points and weights for evaluation of damping matrix
  DRT::UTILS::IntegrationPoints1D gausspoints(MyGaussRule(nnode, integrationtype));

  // matrix to store basis functions and their derivatives evaluated at a certain Gauss point
  LINALG::Matrix<1, nnode> funct;
  LINALG::Matrix<1, nnode> deriv;

  for (int gp = 0; gp < gausspoints.nquad; gp++)
  {
    // evaluate basis functions and their derivatives at current Gauss point
    DRT::UTILS::shape_function_1D(funct, gausspoints.qxg[gp][0], Shape());
    DRT::UTILS::shape_function_1D_deriv1(deriv, gausspoints.qxg[gp][0], Shape());

    // compute point in phyiscal space corresponding to Gauss point
    evaluationpoint.PutScalar(0);
    // loop over all line nodes
    for (int i = 0; i < nnode; i++)
      // loop over all dimensions
      for (int j = 0; j < ndim; j++)
        evaluationpoint(j) += funct(i) * (Nodes()[i]->X()[j] + disp[dof * i + j]);

    // compute velocity and gradient of background flow field at evaluationpoint
    MyBackgroundVelocity<ndim>(params, evaluationpoint, velbackground, velbackgroundgrad);


    // compute tangent vector t_{\par} at current Gauss point
    LINALG::Matrix<ndim, 1> tpar(true);
    for (int i = 0; i < nnode; i++)
      for (int k = 0; k < ndim; k++)
        tpar(k) += deriv(i) * (Nodes()[i]->X()[k] + disp[dof * i + k]) / jacobi[gp];

    // compute velocity vector at this Gauss point
    LINALG::Matrix<ndim, 1> velgp(true);
    for (int i = 0; i < nnode; i++)
      for (int l = 0; l < ndim; l++) velgp(l) += funct(i) * vel[dof * i + l];

    // compute matrix product (t_{\par} \otimes t_{\par}) \cdot velbackgroundgrad
    LINALG::Matrix<ndim, ndim> tpartparvelbackgroundgrad(true);
    for (int i = 0; i < ndim; i++)
      for (int j = 0; j < ndim; j++)
        for (int k = 0; k < ndim; k++)
          tpartparvelbackgroundgrad(i, j) += tpar(i) * tpar(k) * velbackgroundgrad(k, j);

    // loop over all line nodes
    for (int i = 0; i < nnode; i++)
      // loop over lines of matrix t_{\par} \otimes t_{\par}
      for (int k = 0; k < ndim; k++)
        // loop over columns of matrix t_{\par} \otimes t_{\par}
        for (int l = 0; l < ndim; l++)
        {
          if (force != NULL)
            (*force)(i * dof + k) +=
                funct(i) * jacobi[gp] * gausspoints.qwgt[gp] *
                ((k == l) * gamma(1) + (gamma(0) - gamma(1)) * tpar(k) * tpar(l)) *
                (velgp(l) - velbackground(l));

          if (stiffmatrix != NULL)
            // loop over all column nodes
            for (int j = 0; j < nnode; j++)
            {
              (*stiffmatrix)(i * dof + k, j * dof + l) +=
                  gausspoints.qwgt[gp] * funct(i) * funct(j) * jacobi[gp] *
                  ((k == l) * gamma(1) + (gamma(0) - gamma(1)) * tpar(k) * tpar(l)) / dt;
              (*stiffmatrix)(i * dof + k, j * dof + l) -=
                  gausspoints.qwgt[gp] * funct(i) * funct(j) * jacobi[gp] *
                  (velbackgroundgrad(k, l) * gamma(1) +
                      (gamma(0) - gamma(1)) * tpartparvelbackgroundgrad(k, l));
              (*stiffmatrix)(i * dof + k, j * dof + k) += gausspoints.qwgt[gp] * funct(i) *
                                                          deriv(j) * (gamma(0) - gamma(1)) *
                                                          tpar(l) * (velgp(l) - velbackground(l));
              (*stiffmatrix)(i * dof + k, j * dof + l) += gausspoints.qwgt[gp] * funct(i) *
                                                          deriv(j) * (gamma(0) - gamma(1)) *
                                                          tpar(k) * (velgp(l) - velbackground(l));
            }
        }
  }

  return;
}  // DRT::ELEMENTS::Beam3::MyTranslationalDamping(.)

/*-----------------------------------------------------------------------------------------------------------*
 | computes stochastic forces and resulting stiffness (public) cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template <int nnode, int ndim, int dof,
    int randompergauss>  // number of nodes, number of dimensions of embedding space, number of
                         // degrees of freedom per node, number of random numbers required per Gauss
                         // point
inline void
DRT::ELEMENTS::Beam3::MyStochasticForces(Teuchos::ParameterList& params,  //!< parameter list
    const std::vector<double>& vel,         //!< element velocity vector
    const std::vector<double>& disp,        //!< element disp vector
    Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
    Epetra_SerialDenseVector* force)        //!< element internal force vector
{
  // damping coefficients for three translational and one rotatinal degree of freedom
  LINALG::Matrix<3, 1> gamma(true);
  MyDampingConstants(gamma);


  // get vector jacobi with Jacobi determinants at each integration point (gets by default those
  // values required for consistent damping matrix)
  std::vector<double> jacobi(jacobimass_);

  // determine type of numerical integration performed (lumped damping matrix via lobatto
  // integration!)
  IntegrationType integrationtype = gaussexactintegration;

  // get Gauss points and weights for evaluation of damping matrix
  DRT::UTILS::IntegrationPoints1D gausspoints(MyGaussRule(nnode, integrationtype));

  // matrix to store basis functions and their derivatives evaluated at a certain Gauss point
  LINALG::Matrix<1, nnode> funct;
  LINALG::Matrix<1, nnode> deriv;


  /*get pointer at Epetra multivector in parameter list linking to random numbers for stochastic
   * forces with zero mean and standard deviation (2*kT / dt)^0.5; note carefully: a space between
   * the two subsequal ">" signs is mandatory for the C++ parser in order to avoid confusion with
   * ">>" for streams*/
  Teuchos::RCP<Epetra_MultiVector> randomforces = BrownianDynParamsInterface().GetRandomForces();



  for (int gp = 0; gp < gausspoints.nquad; gp++)
  {
    // evaluate basis functions and their derivatives at current Gauss point
    DRT::UTILS::shape_function_1D(funct, gausspoints.qxg[gp][0], Shape());
    DRT::UTILS::shape_function_1D_deriv1(deriv, gausspoints.qxg[gp][0], Shape());

    // compute tangent vector t_{\par} at current Gauss point
    LINALG::Matrix<ndim, 1> tpar(true);
    for (int i = 0; i < nnode; i++)
      for (int k = 0; k < ndim; k++)
        tpar(k) += deriv(i) * (Nodes()[i]->X()[k] + disp[dof * i + k]) / jacobi[gp];


    // loop over all line nodes
    for (int i = 0; i < nnode; i++)
      // loop dimensions with respect to lines
      for (int k = 0; k < ndim; k++)
        // loop dimensions with respect to columns
        for (int l = 0; l < ndim; l++)
        {
          if (force != NULL)
            (*force)(i * dof + k) -= funct(i) *
                                     (sqrt(gamma(1)) * (k == l) +
                                         (sqrt(gamma(0)) - sqrt(gamma(1))) * tpar(k) * tpar(l)) *
                                     (*randomforces)[gp * randompergauss + l][LID()] *
                                     sqrt(jacobi[gp] * gausspoints.qwgt[gp]);

          if (stiffmatrix != NULL)
            // loop over all column nodes
            for (int j = 0; j < nnode; j++)
            {
              (*stiffmatrix)(i * dof + k, j * dof + k) -=
                  funct(i) * deriv(j) * tpar(l) * (*randomforces)[gp * randompergauss + l][LID()] *
                  sqrt(gausspoints.qwgt[gp] / jacobi[gp]) * (sqrt(gamma(0)) - sqrt(gamma(1)));
              (*stiffmatrix)(i * dof + k, j * dof + l) -=
                  funct(i) * deriv(j) * tpar(k) * (*randomforces)[gp * randompergauss + l][LID()] *
                  sqrt(gausspoints.qwgt[gp] / jacobi[gp]) * (sqrt(gamma(0)) - sqrt(gamma(1)));
            }
        }
  }



  return;
}  // DRT::ELEMENTS::Beam3::MyStochasticForces(.)

/*-----------------------------------------------------------------------------------------------------------*
 | computes stochastic moments and (if required) resulting stiffness (public) cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template <int nnode,
    int randompergauss>  // number of nodes, number of random numbers required per Gauss point,
                         // number of random numbers required per Gauss point
inline void
DRT::ELEMENTS::Beam3::MyStochasticMoments(Teuchos::ParameterList& params,  //!< parameter list
    const std::vector<double>& vel,         //!< element velocity vector
    const std::vector<double>& disp,        //!< element disp vector
    Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
    Epetra_SerialDenseVector* force)        //!< element internal force vector
{
  // damping coefficients for three translational and one rotatinal degree of freedom
  LINALG::Matrix<3, 1> gamma(true);
  MyDampingConstants(gamma);

  // determine type of numerical integration performed (note: underintegration applied as for
  // related points triads already known from elasticity)
  IntegrationType integrationtype = gaussunderintegration;

  // get Gauss points and weights for evaluation of damping matrix
  DRT::UTILS::IntegrationPoints1D gausspoints(MyGaussRule(nnode, integrationtype));

  // matrix to store basis functions and their derivatives evaluated at a certain Gauss point
  LINALG::Matrix<1, nnode> funct;

  /*get pointer at Epetra multivector in parameter list linking to random numbers for stochastic
   * forces with zero mean and standard deviation (2*kT / dt)^0.5; note carefully: a space between
   * the two subsequal ">" signs is mandatory for the C++ parser in order to avoid confusion with
   * ">>" for streams*/
  Teuchos::RCP<Epetra_MultiVector> randomforces = BrownianDynParamsInterface().GetRandomForces();

  for (int gp = 0; gp < gausspoints.nquad; gp++)
  {
    // evaluate basis functions and their derivatives at current Gauss point
    DRT::UTILS::shape_function_1D(funct, gausspoints.qxg[gp][0], Shape());

    // get current triad at this Gauss point:
    LINALG::Matrix<3, 3> Tnew;
    LARGEROTATIONS::quaterniontotriad(Qnew_[gp], Tnew);

    // get first column out of Tnew
    LINALG::Matrix<3, 1> t1;
    for (int i = 0; i < 3; i++) t1(i) = Tnew(i, 0);

    // compute spin matrix from first column of Tnew times random number
    LINALG::Matrix<3, 3> S;
    LARGEROTATIONS::computespin(S, t1);
    S.Scale((*randomforces)[gp * randompergauss + 3][LID()]);


    // loop over all line nodes
    for (int i = 0; i < nnode; i++)
      // loop over lines of matrix t_{\par} \otimes t_{\par}
      for (int k = 0; k < 3; k++)
      {
        if (force != NULL)
          (*force)(i * 6 + 3 + k) -= funct(i) * t1(k) *
                                     (*randomforces)[gp * randompergauss + 3][LID()] *
                                     sqrt(jacobi_[gp] * gausspoints.qwgt[gp] * gamma(2));

        if (stiffmatrix != NULL)
          // loop over all column nodes
          for (int j = 0; j < nnode; j++)
            // loop over three dimensions with respect to columns
            for (int l = 0; l < 3; l++)
              (*stiffmatrix)(i * 6 + 3 + k, j * 6 + 3 + l) +=
                  funct(i) * funct(j) * S(k, l) *
                  sqrt(jacobi_[gp] * gausspoints.qwgt[gp] * gamma(2));
      }
  }
  return;
}  // DRT::ELEMENTS::Beam3::MyStochasticMoments(.)

/*-----------------------------------------------------------------------------------------------------------*
 | Assemble stochastic and viscous forces and respective stiffness according to fluctuation
 dissipation      | | theorem (public) cyron 10/09|
 *----------------------------------------------------------------------------------------------------------*/
template <int nnode, int ndim, int dof,
    int randompergauss>  // number of nodes, number of dimensions of embedding space, number of
                         // degrees of freedom per node, number of random numbers required per Gauss
                         // point
inline void
DRT::ELEMENTS::Beam3::CalcBrownianForcesAndStiff(Teuchos::ParameterList& params,
    const std::vector<double>& vel,         //!< element velocity vector
    const std::vector<double>& disp,        //!< element displacement vector
    Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
    Epetra_SerialDenseVector* force)        //!< element internal force vector
{
  // add stiffness and forces due to translational damping effects
  MyTranslationalDamping<nnode, ndim, dof>(params, vel, disp, stiffmatrix, force);

  // add stiffness and forces (i.e. moments) due to rotational damping effects
  MyRotationalDamping<nnode>(params, vel, disp, stiffmatrix, force);

  // add stochastic forces and (if required) resulting stiffness
  MyStochasticForces<nnode, ndim, dof, randompergauss>(params, vel, disp, stiffmatrix, force);

  // add stochastic moments and resulting stiffness
  MyStochasticMoments<nnode, randompergauss>(params, vel, disp, stiffmatrix, force);

  return;

}  // CalcBrownian()
