/*----------------------------------------------------------------------*/
/*! \file
\brief A 2D constraint element with no physics attached
\level 2

\maintainer Matthias Mayr

*----------------------------------------------------------------------*/

#include <Epetra_SerialDenseSolver.h>

#include "constraint_element2.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../linalg/linalg_utils_sparse_algebra_math.H"
#include "../linalg/linalg_serialdensevector.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::ConstraintElement2::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  ActionType act = none;

  // get the required action and distinguish between 2d and 3d MPC's
  std::string action = params.get<std::string>("action", "none");
  if (action == "none")
    return 0;
  else if (action == "calc_MPC_stiff")
  {
    Teuchos::RCP<DRT::Condition> condition = params.get<Teuchos::RCP<DRT::Condition>>("condition");
    const std::string* type = condition->Get<std::string>("control value");

    if (*type == "dist")
      act = calc_MPC_dist_stiff;
    else if (*type == "angle")
      act = calc_MPC_angle_stiff;
    else
      dserror(
          "No constraint type in 2d MPC specified. Value to control should by either be 'dist' or "
          "'angle'!");
  }
  else
    dserror("Unknown type of action for ConstraintElement2");

  switch (act)
  {
    case none:
    {
      return (0);
    }
    break;
    case calc_MPC_dist_stiff:
    {
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vector 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      const int numnode = 3;
      const int numdim = 2;
      LINALG::Matrix<numnode, numdim> xscurr;  // material coord. of element
      SpatialConfiguration(xscurr, mydisp);
      LINALG::Matrix<numdim, 1> elementnormal;
      ComputeNormal(xscurr, elementnormal);
      double normaldistance = ComputeNormalDist(xscurr, elementnormal);
      ComputeFirstDerivDist(xscurr, elevec1, elementnormal);
      ComputeSecondDerivDist(xscurr, elemat1, elementnormal);
      // update corresponding column in "constraint" matrix
      elevec2 = elevec1;
      elevec3[0] = normaldistance;
    }
    break;
    case calc_MPC_angle_stiff:
    {
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp == Teuchos::null) dserror("Cannot get state vector 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp, mydisp, lm);
      const int numnode = 3;
      const int numdim = 2;
      LINALG::Matrix<numnode, numdim> xscurr;  // material coord. of element
      SpatialConfiguration(xscurr, mydisp);

      double angle = ComputeAngle(xscurr);

      ComputeFirstDerivAngle(xscurr, elevec1);
      ComputeSecondDerivAngle(xscurr, elemat1);

      // update corresponding column in "constraint" matrix
      elevec2 = elevec1;
      elevec3[0] = angle;
    }
    break;
    default:
      dserror("Unimplemented type of action");
  }
  return 0;


}  // end of DRT::ELEMENTS::ConstraintElement2::Evaluate

/*----------------------------------------------------------------------*
 * Evaluate Neumann (->dserror) */
int DRT::ELEMENTS::ConstraintElement2::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  dserror("You called Evaluate Neumann of constraint element.");
  return 0;
}


/*----------------------------------------------------------------------*
 * compute 2d normal */
void DRT::ELEMENTS::ConstraintElement2::ComputeNormal(
    const LINALG::Matrix<3, 2>& xc, LINALG::Matrix<2, 1>& elenorm)
{
  elenorm(0, 0) = xc(0, 1) - xc(1, 1);
  elenorm(1, 0) = -xc(0, 0) + xc(1, 0);
  return;
}


/*----------------------------------------------------------------------*
 * normal distance between third point and line */
double DRT::ELEMENTS::ConstraintElement2::ComputeNormalDist(
    const LINALG::Matrix<3, 2>& xc, const LINALG::Matrix<2, 1>& normal)
{
  return (normal(0, 0) * (-xc(0, 0) + xc(2, 0)) - normal(1, 0) * (xc(0, 1) - xc(2, 1))) /
         normal.Norm2();
}

/*----------------------------------------------------------------------*
 * compute angle at second point */
double DRT::ELEMENTS::ConstraintElement2::ComputeAngle(const LINALG::Matrix<3, 2>& xc)
{
  return (acos((xc(0, 1) * (xc(1, 0) - xc(2, 0)) + xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1) +
                   xc(0, 0) * (-xc(1, 1) + xc(2, 1))) /
               sqrt((pow(xc(0, 0) - xc(1, 0), 2) + pow(xc(0, 1) - xc(1, 1), 2)) *
                    (pow(xc(1, 0) - xc(2, 0), 2) + pow(xc(1, 1) - xc(2, 1), 2)))) +
          acos(0.0));
}


/*----------------------------------------------------------------------*
 * second derivatives */
void DRT::ELEMENTS::ConstraintElement2::ComputeFirstDerivDist(const LINALG::Matrix<3, 2>& xc,
    Epetra_SerialDenseVector& elevector, const LINALG::Matrix<2, 1>& normal)
{
  double normcube = pow(normal.Norm2(), 3);

  elevector[0] = (normal(0, 0) * (-pow(xc(1, 0), 2) + xc(0, 0) * (xc(1, 0) - xc(2, 0)) +
                                     xc(1, 0) * xc(2, 0) + normal(0, 0) * (xc(1, 1) - xc(2, 1)))) /
                 normcube;

  elevector[1] = (normal(1, 0) * (-pow(xc(1, 0), 2) + xc(0, 0) * (xc(1, 0) - xc(2, 0)) +
                                     xc(1, 0) * xc(2, 0) + normal(0, 0) * (xc(1, 1) - xc(2, 1)))) /
                 normcube;

  elevector[2] = -((normal(0, 0) * (pow(xc(0, 0), 2) + pow(xc(0, 1), 2) + xc(1, 0) * xc(2, 0) -
                                       xc(0, 0) * (xc(1, 0) + xc(2, 0)) + xc(1, 1) * xc(2, 1) -
                                       xc(0, 1) * (xc(1, 1) + xc(2, 1)))) /
                   normcube);

  elevector[3] = -((normal(1, 0) * (pow(xc(0, 0), 2) + pow(xc(0, 1), 2) + xc(1, 0) * xc(2, 0) -
                                       xc(0, 0) * (xc(1, 0) + xc(2, 0)) + xc(1, 1) * xc(2, 1) -
                                       xc(0, 1) * (xc(1, 1) + xc(2, 1)))) /
                   normcube);

  elevector[4] = normal(0, 0) / normal.Norm2();

  elevector[5] = normal(1, 0) / normal.Norm2();
  elevector.Scale(-1.0);
  return;
}

/*----------------------------------------------------------------------*
 * first derivatives */
void DRT::ELEMENTS::ConstraintElement2::ComputeFirstDerivAngle(
    const LINALG::Matrix<3, 2>& xc, Epetra_SerialDenseVector& elevector)
{
  LINALG::SerialDenseVector vec1(2);
  vec1[1] = xc(0, 0) - xc(1, 0);
  vec1[0] = -(xc(0, 1) - xc(1, 1));

  LINALG::SerialDenseVector vec2(2);
  vec2[0] = -xc(1, 0) + xc(2, 0);
  vec2[1] = -xc(1, 1) + xc(2, 1);

  const double vec1normsquare = pow(vec1.Norm2(), 2);
  const double vec2normsquare = pow(vec2.Norm2(), 2);

  elevector[0] = -((vec2[1] / sqrt(vec1normsquare * vec2normsquare) -
                       (vec2normsquare * vec1[1] *
                           (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                               xc(1, 0) * xc(2, 1))) /
                           pow(vec1normsquare * vec2normsquare, 1.5)) /
                   sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                    xc(1, 0) * xc(2, 1),
                                2) /
                                (vec1normsquare * vec2normsquare)));
  elevector[1] = ((-(vec2[0] / sqrt(vec1normsquare * vec2normsquare)) +
                      (vec2normsquare * vec1[0] *
                          (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                              xc(1, 0) * xc(2, 1))) /
                          pow(vec1normsquare * vec2normsquare, 1.5)) /
                  sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                   xc(1, 0) * xc(2, 1),
                               2) /
                               (vec1normsquare * vec2normsquare)));
  elevector[2] = (((xc(0, 1) - xc(2, 1)) / sqrt(vec1normsquare * vec2normsquare) -
                      ((-2 * vec2normsquare * vec1[1] - 2 * vec1normsquare * vec2[0]) *
                          (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                              xc(1, 0) * xc(2, 1))) /
                          (2. * pow(vec1normsquare * vec2normsquare, 1.5))) /
                  sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                   xc(1, 0) * xc(2, 1),
                               2) /
                               (vec1normsquare * vec2normsquare)));
  elevector[3] = (((-xc(0, 0) + xc(2, 0)) / sqrt(vec1normsquare * vec2normsquare) -
                      ((2 * vec2normsquare * vec1[0] - 2 * vec1normsquare * vec2[1]) *
                          (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                              xc(1, 0) * xc(2, 1))) /
                          (2. * pow(vec1normsquare * vec2normsquare, 1.5))) /
                  sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                   xc(1, 0) * xc(2, 1),
                               2) /
                               (vec1normsquare * vec2normsquare)));
  elevector[4] = (((-xc(0, 1) + xc(1, 1)) / sqrt(vec1normsquare * vec2normsquare) -
                      (vec1normsquare * vec2[0] *
                          (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                              xc(1, 0) * xc(2, 1))) /
                          pow(vec1normsquare * vec2normsquare, 1.5)) /
                  sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                   xc(1, 0) * xc(2, 1),
                               2) /
                               (vec1normsquare * vec2normsquare)));
  elevector[5] = ((vec1[1] / sqrt(vec1normsquare * vec2normsquare) -
                      (vec1normsquare * vec2(1) *
                          (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                              xc(1, 0) * xc(2, 1))) /
                          pow(vec1normsquare * vec2normsquare, 1.5)) /
                  sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                   xc(1, 0) * xc(2, 1),
                               2) /
                               (vec1normsquare * vec2normsquare)));
}

/*----------------------------------------------------------------------*
 * second derivatives */
void DRT::ELEMENTS::ConstraintElement2::ComputeSecondDerivDist(const LINALG::Matrix<3, 2>& xc,
    Epetra_SerialDenseMatrix& elematrix, const LINALG::Matrix<2, 1>& normal)
{
  double normsquare = pow(normal.Norm2(), 2);
  double normcube = pow(normal.Norm2(), 3);
  double normpowfour = pow(normal.Norm2(), 4);
  double normpowfive = pow(normal.Norm2(), 5);

  elematrix(0, 0) =
      (normal(0, 0) *
          (-2 * pow(xc(1, 0), 3) - 2 * xc(1, 0) * pow(xc(1, 1), 2) -
              2 * pow(xc(0, 0), 2) * (xc(1, 0) - xc(2, 0)) +
              pow(xc(0, 1), 2) * (xc(1, 0) - xc(2, 0)) + 2 * pow(xc(1, 0), 2) * xc(2, 0) -
              pow(xc(1, 1), 2) * xc(2, 0) +
              xc(0, 1) * (2 * xc(1, 1) * xc(2, 0) + xc(1, 0) * (xc(1, 1) - 3 * xc(2, 1))) +
              3 * xc(1, 0) * xc(1, 1) * xc(2, 1) +
              xc(0, 0) * (4 * pow(xc(1, 0), 2) - 4 * xc(1, 0) * xc(2, 0) +
                             3 * normal(0, 0) * (-xc(1, 1) + xc(2, 1))))) /
      normpowfive;

  elematrix(0, 1) = (3 * pow(normal(0, 0), 2) * normal(1, 0) * (xc(0, 0) - xc(2, 0)) +
                        normsquare * normal(1, 0) * (-xc(1, 0) + xc(2, 0)) +
                        normal(0, 0) * (3 * pow(normal(1, 0), 2) * (xc(0, 1) - xc(2, 1)) +
                                           normsquare * (-xc(1, 1) + xc(2, 1)))) /
                    normpowfive;

  elematrix(0, 2) =
      (normal(0, 0) *
          (pow(xc(0, 0), 3) + pow(xc(1, 0), 3) + xc(1, 0) * pow(xc(1, 1), 2) -
              2 * pow(xc(1, 0), 2) * xc(2, 0) + pow(xc(1, 1), 2) * xc(2, 0) +
              pow(xc(0, 1), 2) * (-2 * xc(1, 0) + xc(2, 0)) -
              pow(xc(0, 0), 2) * (xc(1, 0) + 2 * xc(2, 0)) - 3 * xc(1, 0) * xc(1, 1) * xc(2, 1) +
              xc(0, 0) * (pow(xc(0, 1), 2) - pow(xc(1, 0), 2) - 2 * pow(xc(1, 1), 2) +
                             4 * xc(1, 0) * xc(2, 0) + xc(0, 1) * (xc(1, 1) - 3 * xc(2, 1)) +
                             3 * xc(1, 1) * xc(2, 1)) +
              xc(0, 1) * (-2 * xc(1, 1) * xc(2, 0) + xc(1, 0) * (xc(1, 1) + 3 * xc(2, 1))))) /
      normpowfive;

  elematrix(0, 3) =
      (normpowfour +
          normsquare *
              (normal(1, 0) * (xc(0, 0) - xc(2, 0)) + normal(0, 0) * (xc(1, 1) - xc(2, 1))) +
          3 * normal(0, 0) * normal(1, 0) *
              (normal(0, 0) * (-xc(0, 0) + xc(2, 0)) + normal(1, 0) * (-xc(0, 1) + xc(2, 1)))) /
      normpowfive;

  elematrix(0, 4) = (normal(0, 0) * normal(1, 0)) / normcube;

  elematrix(0, 5) = -(pow(normal(0, 0), 2) / normcube);

  elematrix(1, 0) = (3 * pow(normal(0, 0), 2) * normal(1, 0) * (xc(0, 0) - xc(2, 0)) +
                        normsquare * normal(1, 0) * (-xc(1, 0) + xc(2, 0)) +
                        normal(0, 0) * (3 * pow(normal(1, 0), 2) * (xc(0, 1) - xc(2, 1)) +
                                           normsquare * (-xc(1, 1) + xc(2, 1)))) /
                    normpowfive;

  elematrix(1, 1) =
      (normal(1, 0) *
          (-2 * pow(xc(1, 0), 2) * xc(1, 1) - 2 * pow(xc(1, 1), 3) +
              3 * xc(1, 0) * xc(1, 1) * xc(2, 0) +
              xc(0, 1) * (3 * pow(xc(1, 0), 2) - 3 * xc(1, 0) * xc(2, 0) +
                             4 * xc(1, 1) * (xc(1, 1) - xc(2, 1))) +
              pow(xc(0, 0), 2) * (xc(1, 1) - xc(2, 1)) -
              2 * pow(xc(0, 1), 2) * (xc(1, 1) - xc(2, 1)) - pow(xc(1, 0), 2) * xc(2, 1) +
              2 * pow(xc(1, 1), 2) * xc(2, 1) +
              xc(0, 0) * (-3 * xc(0, 1) * (xc(1, 0) - xc(2, 0)) - 3 * xc(1, 1) * xc(2, 0) +
                             xc(1, 0) * (xc(1, 1) + 2 * xc(2, 1))))) /
      normpowfive;

  elematrix(1, 2) =
      (-normpowfour +
          normsquare *
              (normal(1, 0) * (xc(1, 0) - xc(2, 0)) + normal(0, 0) * (xc(0, 1) - xc(2, 1))) +
          3 * normal(0, 0) * normal(1, 0) *
              (normal(0, 0) * (-xc(0, 0) + xc(2, 0)) + normal(1, 0) * (-xc(0, 1) + xc(2, 1)))) /
      normpowfive;

  elematrix(1, 3) =
      (normal(1, 0) *
          (pow(xc(0, 1), 3) + pow(xc(1, 0), 2) * xc(1, 1) + pow(xc(1, 1), 3) -
              3 * xc(1, 0) * xc(1, 1) * xc(2, 0) -
              xc(0, 1) * (2 * pow(xc(1, 0), 2) - 3 * xc(1, 0) * xc(2, 0) +
                             xc(1, 1) * (xc(1, 1) - 4 * xc(2, 1))) +
              xc(0, 0) * (xc(0, 1) * (xc(1, 0) - 3 * xc(2, 0)) + 3 * xc(1, 1) * xc(2, 0) +
                             xc(1, 0) * (xc(1, 1) - 2 * xc(2, 1))) +
              pow(xc(1, 0), 2) * xc(2, 1) - 2 * pow(xc(1, 1), 2) * xc(2, 1) +
              pow(xc(0, 0), 2) * (xc(0, 1) - 2 * xc(1, 1) + xc(2, 1)) -
              pow(xc(0, 1), 2) * (xc(1, 1) + 2 * xc(2, 1)))) /
      normpowfive;

  elematrix(1, 4) = pow(normal(1, 0), 2) / normcube;

  elematrix(1, 5) = -((normal(0, 0) * normal(1, 0)) / normcube);

  elematrix(2, 0) =
      (normal(0, 0) *
          (pow(xc(0, 0), 3) + pow(xc(1, 0), 3) + xc(1, 0) * pow(xc(1, 1), 2) -
              2 * pow(xc(1, 0), 2) * xc(2, 0) + pow(xc(1, 1), 2) * xc(2, 0) +
              pow(xc(0, 1), 2) * (-2 * xc(1, 0) + xc(2, 0)) -
              pow(xc(0, 0), 2) * (xc(1, 0) + 2 * xc(2, 0)) - 3 * xc(1, 0) * xc(1, 1) * xc(2, 1) +
              xc(0, 0) * (pow(xc(0, 1), 2) - pow(xc(1, 0), 2) - 2 * pow(xc(1, 1), 2) +
                             4 * xc(1, 0) * xc(2, 0) + xc(0, 1) * (xc(1, 1) - 3 * xc(2, 1)) +
                             3 * xc(1, 1) * xc(2, 1)) +
              xc(0, 1) * (-2 * xc(1, 1) * xc(2, 0) + xc(1, 0) * (xc(1, 1) + 3 * xc(2, 1))))) /
      normpowfive;

  elematrix(2, 1) =
      (-normpowfour +
          normsquare *
              (normal(1, 0) * (xc(1, 0) - xc(2, 0)) + normal(0, 0) * (xc(0, 1) - xc(2, 1))) +
          3 * normal(0, 0) * normal(1, 0) *
              (normal(0, 0) * (-xc(0, 0) + xc(2, 0)) + normal(1, 0) * (-xc(0, 1) + xc(2, 1)))) /
      normpowfive;

  elematrix(2, 2) =
      -((normal(0, 0) *
            (2 * pow(xc(0, 0), 3) - 2 * pow(xc(1, 0), 2) * xc(2, 0) + pow(xc(1, 1), 2) * xc(2, 0) +
                pow(xc(0, 1), 2) * (-3 * xc(1, 0) + xc(2, 0)) -
                2 * pow(xc(0, 0), 2) * (2 * xc(1, 0) + xc(2, 0)) -
                3 * xc(1, 0) * xc(1, 1) * xc(2, 1) +
                xc(0, 1) * (-2 * xc(1, 1) * xc(2, 0) + 3 * xc(1, 0) * (xc(1, 1) + xc(2, 1))) +
                xc(0, 0) * (2 * pow(xc(0, 1), 2) + 2 * pow(xc(1, 0), 2) - pow(xc(1, 1), 2) +
                               4 * xc(1, 0) * xc(2, 0) + 3 * xc(1, 1) * xc(2, 1) -
                               xc(0, 1) * (xc(1, 1) + 3 * xc(2, 1))))) /
          normpowfive);

  elematrix(2, 3) =
      (3 * normal(0, 0) * normal(1, 0) *
              (normal(0, 0) * (xc(0, 0) - xc(2, 0)) + normal(1, 0) * (xc(0, 1) - xc(2, 1))) +
          normsquare *
              (normal(1, 0) * (-xc(0, 0) + xc(2, 0)) + normal(0, 0) * (-xc(0, 1) + xc(2, 1)))) /
      normpowfive;

  elematrix(2, 4) = -((normal(0, 0) * normal(1, 0)) / normcube);

  elematrix(2, 5) = pow(normal(0, 0), 2) / normcube;

  elematrix(3, 0) =
      (normpowfour +
          normsquare *
              (normal(1, 0) * (xc(0, 0) - xc(2, 0)) + normal(0, 0) * (xc(1, 1) - xc(2, 1))) +
          3 * normal(0, 0) * normal(1, 0) *
              (normal(0, 0) * (-xc(0, 0) + xc(2, 0)) + normal(1, 0) * (-xc(0, 1) + xc(2, 1)))) /
      normpowfive;

  elematrix(3, 1) =
      (normal(1, 0) *
          (pow(xc(0, 1), 3) + pow(xc(1, 0), 2) * xc(1, 1) + pow(xc(1, 1), 3) -
              3 * xc(1, 0) * xc(1, 1) * xc(2, 0) -
              xc(0, 1) * (2 * pow(xc(1, 0), 2) - 3 * xc(1, 0) * xc(2, 0) +
                             xc(1, 1) * (xc(1, 1) - 4 * xc(2, 1))) +
              xc(0, 0) * (xc(0, 1) * (xc(1, 0) - 3 * xc(2, 0)) + 3 * xc(1, 1) * xc(2, 0) +
                             xc(1, 0) * (xc(1, 1) - 2 * xc(2, 1))) +
              pow(xc(1, 0), 2) * xc(2, 1) - 2 * pow(xc(1, 1), 2) * xc(2, 1) +
              pow(xc(0, 0), 2) * (xc(0, 1) - 2 * xc(1, 1) + xc(2, 1)) -
              pow(xc(0, 1), 2) * (xc(1, 1) + 2 * xc(2, 1)))) /
      normpowfive;

  elematrix(3, 2) =
      (3 * normal(0, 0) * normal(1, 0) *
              (normal(0, 0) * (xc(0, 0) - xc(2, 0)) + normal(1, 0) * (xc(0, 1) - xc(2, 1))) +
          normsquare *
              (normal(1, 0) * (-xc(0, 0) + xc(2, 0)) + normal(0, 0) * (-xc(0, 1) + xc(2, 1)))) /
      normpowfive;

  elematrix(3, 3) =
      -((normal(1, 0) *
            (2 * pow(xc(0, 1), 3) - 3 * xc(1, 0) * xc(1, 1) * xc(2, 0) +
                pow(xc(1, 0), 2) * xc(2, 1) - 2 * pow(xc(1, 1), 2) * xc(2, 1) +
                pow(xc(0, 0), 2) * (2 * xc(0, 1) - 3 * xc(1, 1) + xc(2, 1)) -
                2 * pow(xc(0, 1), 2) * (2 * xc(1, 1) + xc(2, 1)) -
                xc(0, 0) * (-3 * xc(1, 1) * xc(2, 0) + xc(0, 1) * (xc(1, 0) + 3 * xc(2, 0)) +
                               xc(1, 0) * (-3 * xc(1, 1) + 2 * xc(2, 1))) +
                xc(0, 1) * (-pow(xc(1, 0), 2) + 3 * xc(1, 0) * xc(2, 0) +
                               2 * xc(1, 1) * (xc(1, 1) + 2 * xc(2, 1))))) /
          normpowfive);

  elematrix(3, 4) = -(pow(normal(1, 0), 2) / normcube);

  elematrix(3, 5) = (normal(0, 0) * (-xc(0, 0) + xc(1, 0))) / normcube;

  elematrix(4, 0) = (normal(0, 0) * normal(1, 0)) / normcube;

  elematrix(4, 1) = pow(normal(1, 0), 2) / normcube;

  elematrix(4, 2) = -((normal(0, 0) * normal(1, 0)) / normcube);

  elematrix(4, 3) = -(pow(normal(1, 0), 2) / normcube);

  elematrix(4, 4) = 0;

  elematrix(4, 5) = 0;

  elematrix(5, 0) = -(pow(normal(0, 0), 2) / normcube);

  elematrix(5, 1) = -((normal(0, 0) * normal(1, 0)) / normcube);

  elematrix(5, 2) = pow(normal(0, 0), 2) / normcube;

  elematrix(5, 3) = (normal(0, 0) * (-xc(0, 0) + xc(1, 0))) / normcube;

  elematrix(5, 4) = 0;

  elematrix(5, 5) = 0;
  return;

  elematrix.Scale(-1.0);
}

/*----------------------------------------------------------------------*
 * second derivatives */
void DRT::ELEMENTS::ConstraintElement2::ComputeSecondDerivAngle(
    const LINALG::Matrix<3, 2>& xc, Epetra_SerialDenseMatrix& elematrix)
{
  LINALG::SerialDenseVector vec1(2);
  vec1[1] = xc(0, 0) - xc(1, 0);
  vec1[0] = -(xc(0, 1) - xc(1, 1));

  LINALG::SerialDenseVector vec2(2);
  vec2[0] = -xc(1, 0) + xc(2, 0);
  vec2[1] = -xc(1, 1) + xc(2, 1);

  const double vec1sq = pow(vec1.Norm2(), 2);
  const double vec2sq = pow(vec2.Norm2(), 2);

  elematrix(0, 0) =
      -(((-2 * vec2sq * vec1[1] * vec2[1]) / pow(vec1sq * vec2sq, 1.5) -
            (vec2sq * (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 1.5) +
            (3 * pow(vec2sq, 2) * pow(vec1[1], 2) *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 2.5)) /
          sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                           xc(1, 0) * xc(2, 1),
                       2) /
                       (vec1sq * vec2sq))) +
      ((vec2[1] / sqrt(vec1sq * vec2sq) - (vec2sq * vec1[1] *
                                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                                  xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1))) /
                                              pow(vec1sq * vec2sq, 1.5)) *
          ((-2 * vec2[1] *
               (vec2(1) * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                   xc(1, 0) * xc(2, 1))) /
                  (vec1sq * vec2sq) +
              (2 * vec1[1] *
                  pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1),
                      2)) /
                  (pow(vec1sq, 2) * vec2sq))) /
          (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                xc(1, 0) * xc(2, 1),
                            2) /
                            (vec1sq * vec2sq),
                    1.5));
  elematrix(0, 1) =
      -(((vec2sq * vec1[1] * vec2[0]) / pow(vec1sq * vec2sq, 1.5) +
            (vec2sq * vec1[0] * vec2[1]) / pow(vec1sq * vec2sq, 1.5) -
            (3 * pow(vec2sq, 2) * vec1[0] * vec1[1] *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 2.5)) /
          sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                           xc(1, 0) * xc(2, 1),
                       2) /
                       (vec1sq * vec2sq))) +
      ((vec2[1] / sqrt(vec1sq * vec2sq) - (vec2sq * vec1[1] *
                                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                                  xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1))) /
                                              pow(vec1sq * vec2sq, 1.5)) *
          ((2 * vec2[0] *
               (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                   xc(1, 0) * xc(2, 1))) /
                  (vec1sq * vec2sq) -
              (2 * vec1[0] *
                  pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1),
                      2)) /
                  (pow(vec1sq, 2) * vec2sq))) /
          (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                xc(1, 0) * xc(2, 1),
                            2) /
                            (vec1sq * vec2sq),
                    1.5));
  elematrix(0, 2) =
      -((-((-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0]) * vec2[1]) /
                (2. * pow(vec1sq * vec2sq, 1.5)) -
            (vec2sq * vec1[1] * (xc(0, 1) - xc(2, 1))) / pow(vec1sq * vec2sq, 1.5) +
            (vec2sq * (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 1.5) +
            (2 * vec1[1] * vec2[0] *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 1.5) +
            (3 * vec2sq * vec1[1] * (-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0]) *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                (2. * pow(vec1sq * vec2sq, 2.5))) /
          sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                           xc(1, 0) * xc(2, 1),
                       2) /
                       (vec1sq * vec2sq))) +
      ((vec2[1] / sqrt(vec1sq * vec2sq) - (vec2sq * vec1[1] *
                                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                                  xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1))) /
                                              pow(vec1sq * vec2sq, 1.5)) *
          ((-2 * (xc(0, 1) - xc(2, 1)) *
               (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                   xc(1, 0) * xc(2, 1))) /
                  (vec1sq * vec2sq) -
              (2 * vec1[1] *
                  pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1),
                      2)) /
                  (pow(vec1sq, 2) * vec2sq) -
              (2 * vec2[0] *
                  pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1),
                      2)) /
                  (vec1sq * pow(vec2sq, 2)))) /
          (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                xc(1, 0) * xc(2, 1),
                            2) /
                            (vec1sq * vec2sq),
                    1.5));
  elematrix(0, 3) =
      -((-(1 / sqrt(vec1sq * vec2sq)) -
            (vec2[1] * (2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1])) /
                (2. * pow(vec1sq * vec2sq, 1.5)) -
            (vec2sq * vec1[1] * (-xc(0, 0) + xc(2, 0))) / pow(vec1sq * vec2sq, 1.5) +
            (2 * vec1[1] * vec2(1) *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 1.5) +
            (3 * vec2sq * vec1[1] * (2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1]) *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                (2. * pow(vec1sq * vec2sq, 2.5))) /
          sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                           xc(1, 0) * xc(2, 1),
                       2) /
                       (vec1sq * vec2sq))) +
      ((vec2[1] / sqrt(vec1sq * vec2sq) - (vec2sq * vec1[1] *
                                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                                  xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1))) /
                                              pow(vec1sq * vec2sq, 1.5)) *
          ((-2 * (-xc(0, 0) + xc(2, 0)) *
               (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                   xc(1, 0) * xc(2, 1))) /
                  (vec1sq * vec2sq) +
              (2 * vec1[0] *
                  pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1),
                      2)) /
                  (pow(vec1sq, 2) * vec2sq) -
              (2 * vec2[1] *
                  pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1),
                      2)) /
                  (vec1sq * pow(vec2sq, 2)))) /
          (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                xc(1, 0) * xc(2, 1),
                            2) /
                            (vec1sq * vec2sq),
                    1.5));
  elematrix(0, 4) =
      -((-((vec1sq * vec2[0] * vec2[1]) / pow(vec1sq * vec2sq, 1.5)) -
            (vec2sq * vec1[1] * (-xc(0, 1) + xc(1, 1))) / pow(vec1sq * vec2sq, 1.5) +
            (3 * vec1sq * vec2sq * vec1[1] * vec2[0] *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 2.5) -
            (2 * vec1[1] * vec2[0] *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 1.5)) /
          sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                           xc(1, 0) * xc(2, 1),
                       2) /
                       (vec1sq * vec2sq))) +
      ((vec2[1] / sqrt(vec1sq * vec2sq) - (vec2sq * vec1[1] *
                                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                                  xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1))) /
                                              pow(vec1sq * vec2sq, 1.5)) *
          ((-2 * (-xc(0, 1) + xc(1, 1)) *
               (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                   xc(1, 0) * xc(2, 1))) /
                  (vec1sq * vec2sq) +
              (2 * vec2[0] *
                  pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1),
                      2)) /
                  (vec1sq * pow(vec2sq, 2)))) /
          (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                xc(1, 0) * xc(2, 1),
                            2) /
                            (vec1sq * vec2sq),
                    1.5));
  elematrix(0, 5) =
      -((1 / sqrt(vec1sq * vec2sq) - (vec2sq * pow(vec1[1], 2)) / pow(vec1sq * vec2sq, 1.5) -
            (vec1sq * pow(vec2[1], 2)) / pow(vec1sq * vec2sq, 1.5) +
            (3 * vec1sq * vec2sq * vec1[1] * vec2(1) *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 2.5) -
            (2 * vec1[1] * vec2(1) *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 1.5)) /
          sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                           xc(1, 0) * xc(2, 1),
                       2) /
                       (vec1sq * vec2sq))) +
      ((vec2[1] / sqrt(vec1sq * vec2sq) - (vec2sq * vec1[1] *
                                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                                  xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1))) /
                                              pow(vec1sq * vec2sq, 1.5)) *
          ((-2 * vec1[1] *
               (vec2(1) * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                   xc(1, 0) * xc(2, 1))) /
                  (vec1sq * vec2sq) +
              (2 * vec2[1] *
                  pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1),
                      2)) /
                  (vec1sq * pow(vec2sq, 2)))) /
          (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                xc(1, 0) * xc(2, 1),
                            2) /
                            (vec1sq * vec2sq),
                    1.5));
  elematrix(1, 0) =
      -(((vec2sq * vec1[1] * vec2[0]) / pow(vec1sq * vec2sq, 1.5) +
            (vec2sq * vec1[0] * vec2[1]) / pow(vec1sq * vec2sq, 1.5) -
            (3 * pow(vec2sq, 2) * vec1[0] * vec1[1] *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 2.5)) /
          sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                           xc(1, 0) * xc(2, 1),
                       2) /
                       (vec1sq * vec2sq))) +
      ((-(vec2[0] / sqrt(vec1sq * vec2sq)) + (vec2sq * vec1[0] *
                                                 (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                                     xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1))) /
                                                 pow(vec1sq * vec2sq, 1.5)) *
          ((-2 * vec2[1] *
               (vec2(1) * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                   xc(1, 0) * xc(2, 1))) /
                  (vec1sq * vec2sq) +
              (2 * vec1[1] *
                  pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1),
                      2)) /
                  (pow(vec1sq, 2) * vec2sq))) /
          (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                xc(1, 0) * xc(2, 1),
                            2) /
                            (vec1sq * vec2sq),
                    1.5));
  elematrix(1, 1) =
      -(((-2 * vec2sq * vec1[0] * vec2[0]) / pow(vec1sq * vec2sq, 1.5) -
            (vec2sq * (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 1.5) +
            (3 * pow(vec2sq, 2) * pow(vec1[0], 2) *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 2.5)) /
          sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                           xc(1, 0) * xc(2, 1),
                       2) /
                       (vec1sq * vec2sq))) +
      ((-(vec2[0] / sqrt(vec1sq * vec2sq)) + (vec2sq * vec1[0] *
                                                 (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                                     xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1))) /
                                                 pow(vec1sq * vec2sq, 1.5)) *
          ((2 * vec2[0] *
               (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                   xc(1, 0) * xc(2, 1))) /
                  (vec1sq * vec2sq) -
              (2 * vec1[0] *
                  pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1),
                      2)) /
                  (pow(vec1sq, 2) * vec2sq))) /
          (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                xc(1, 0) * xc(2, 1),
                            2) /
                            (vec1sq * vec2sq),
                    1.5));
  elematrix(1, 2) =
      -((1 / sqrt(vec1sq * vec2sq) +
            (vec2[0] * (-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0])) /
                (2. * pow(vec1sq * vec2sq, 1.5)) +
            (vec2sq * vec1[0] * (xc(0, 1) - xc(2, 1))) / pow(vec1sq * vec2sq, 1.5) -
            (2 * vec1[0] * vec2[0] *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 1.5) -
            (3 * vec2sq * vec1[0] * (-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0]) *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                (2. * pow(vec1sq * vec2sq, 2.5))) /
          sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                           xc(1, 0) * xc(2, 1),
                       2) /
                       (vec1sq * vec2sq))) +
      ((-(vec2[0] / sqrt(vec1sq * vec2sq)) + (vec2sq * vec1[0] *
                                                 (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                                     xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1))) /
                                                 pow(vec1sq * vec2sq, 1.5)) *
          ((-2 * (xc(0, 1) - xc(2, 1)) *
               (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                   xc(1, 0) * xc(2, 1))) /
                  (vec1sq * vec2sq) -
              (2 * vec1[1] *
                  pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1),
                      2)) /
                  (pow(vec1sq, 2) * vec2sq) -
              (2 * vec2[0] *
                  pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1),
                      2)) /
                  (vec1sq * pow(vec2sq, 2)))) /
          (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                xc(1, 0) * xc(2, 1),
                            2) /
                            (vec1sq * vec2sq),
                    1.5));
  elematrix(1, 3) =
      -(((vec2[0] * (2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1])) /
                (2. * pow(vec1sq * vec2sq, 1.5)) +
            (vec2sq * vec1[0] * (-xc(0, 0) + xc(2, 0))) / pow(vec1sq * vec2sq, 1.5) +
            (vec2sq * (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 1.5) -
            (2 * vec1[0] * vec2(1) *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 1.5) -
            (3 * vec2sq * vec1[0] * (2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1]) *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                (2. * pow(vec1sq * vec2sq, 2.5))) /
          sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                           xc(1, 0) * xc(2, 1),
                       2) /
                       (vec1sq * vec2sq))) +
      ((-(vec2[0] / sqrt(vec1sq * vec2sq)) + (vec2sq * vec1[0] *
                                                 (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                                     xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1))) /
                                                 pow(vec1sq * vec2sq, 1.5)) *
          ((-2 * (-xc(0, 0) + xc(2, 0)) *
               (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                   xc(1, 0) * xc(2, 1))) /
                  (vec1sq * vec2sq) +
              (2 * vec1[0] *
                  pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1),
                      2)) /
                  (pow(vec1sq, 2) * vec2sq) -
              (2 * vec2[1] *
                  pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1),
                      2)) /
                  (vec1sq * pow(vec2sq, 2)))) /
          (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                xc(1, 0) * xc(2, 1),
                            2) /
                            (vec1sq * vec2sq),
                    1.5));
  elematrix(1, 4) =
      -((-(1 / sqrt(vec1sq * vec2sq)) + (vec1sq * pow(vec2[0], 2)) / pow(vec1sq * vec2sq, 1.5) +
            (vec2sq * vec1[0] * (-xc(0, 1) + xc(1, 1))) / pow(vec1sq * vec2sq, 1.5) -
            (3 * vec1sq * vec2sq * vec1[0] * vec2[0] *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 2.5) +
            (2 * vec1[0] * vec2[0] *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 1.5)) /
          sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                           xc(1, 0) * xc(2, 1),
                       2) /
                       (vec1sq * vec2sq))) +
      ((-(vec2[0] / sqrt(vec1sq * vec2sq)) + (vec2sq * vec1[0] *
                                                 (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                                     xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1))) /
                                                 pow(vec1sq * vec2sq, 1.5)) *
          ((-2 * (-xc(0, 1) + xc(1, 1)) *
               (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                   xc(1, 0) * xc(2, 1))) /
                  (vec1sq * vec2sq) +
              (2 * vec2[0] *
                  pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1),
                      2)) /
                  (vec1sq * pow(vec2sq, 2)))) /
          (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                xc(1, 0) * xc(2, 1),
                            2) /
                            (vec1sq * vec2sq),
                    1.5));
  elematrix(1, 5) =
      -(((vec2sq * vec1[0] * vec1[1]) / pow(vec1sq * vec2sq, 1.5) +
            (vec1sq * vec2[0] * vec2[1]) / pow(vec1sq * vec2sq, 1.5) -
            (3 * vec1sq * vec2sq * vec1[0] * vec2(1) *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 2.5) +
            (2 * vec1[0] * vec2(1) *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 1.5)) /
          sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                           xc(1, 0) * xc(2, 1),
                       2) /
                       (vec1sq * vec2sq))) +
      ((-(vec2[0] / sqrt(vec1sq * vec2sq)) + (vec2sq * vec1[0] *
                                                 (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                                     xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1))) /
                                                 pow(vec1sq * vec2sq, 1.5)) *
          ((-2 * vec1[1] *
               (vec2(1) * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                   xc(1, 0) * xc(2, 1))) /
                  (vec1sq * vec2sq) +
              (2 * vec2[1] *
                  pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1),
                      2)) /
                  (vec1sq * pow(vec2sq, 2)))) /
          (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                xc(1, 0) * xc(2, 1),
                            2) /
                            (vec1sq * vec2sq),
                    1.5));
  elematrix(2, 0) = -((-((-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0]) * vec2[1]) /
                              (2. * pow(vec1sq * vec2sq, 1.5)) -
                          (vec2sq * vec1[1] * (xc(0, 1) - xc(2, 1))) / pow(vec1sq * vec2sq, 1.5) +
                          (3 * vec2sq * vec1[1] * (-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0]) *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              (2. * pow(vec1sq * vec2sq, 2.5)) -
                          ((-2 * vec2sq - 4 * vec1[1] * vec2[0]) *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              (2. * pow(vec1sq * vec2sq, 1.5))) /
                        sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                         xc(1, 0) * xc(2, 1),
                                     2) /
                                     (vec1sq * vec2sq))) +
                    (((xc(0, 1) - xc(2, 1)) / sqrt(vec1sq * vec2sq) -
                         ((-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0]) *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                             (2. * pow(vec1sq * vec2sq, 1.5))) *
                        ((-2 * vec2[1] *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                                (vec1sq * vec2sq) +
                            (2 * vec1[1] *
                                pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                        xc(1, 0) * xc(2, 1),
                                    2)) /
                                (pow(vec1sq, 2) * vec2sq))) /
                        (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                              xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1),
                                          2) /
                                          (vec1sq * vec2sq),
                                  1.5));
  elematrix(2, 1) = -((1 / sqrt(vec1sq * vec2sq) +
                          (vec2[0] * (-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0])) /
                              (2. * pow(vec1sq * vec2sq, 1.5)) +
                          (vec2sq * vec1[0] * (xc(0, 1) - xc(2, 1))) / pow(vec1sq * vec2sq, 1.5) -
                          (2 * vec1[0] * vec2[0] *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              pow(vec1sq * vec2sq, 1.5) -
                          (3 * vec2sq * vec1[0] * (-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0]) *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              (2. * pow(vec1sq * vec2sq, 2.5))) /
                        sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                         xc(1, 0) * xc(2, 1),
                                     2) /
                                     (vec1sq * vec2sq))) +
                    (((xc(0, 1) - xc(2, 1)) / sqrt(vec1sq * vec2sq) -
                         ((-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0]) *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                             (2. * pow(vec1sq * vec2sq, 1.5))) *
                        ((2 * vec2[0] *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                                (vec1sq * vec2sq) -
                            (2 * vec1[0] *
                                pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                        xc(1, 0) * xc(2, 1),
                                    2)) /
                                (pow(vec1sq, 2) * vec2sq))) /
                        (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                              xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1),
                                          2) /
                                          (vec1sq * vec2sq),
                                  1.5));
  elematrix(2, 2) = -((-(((-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0]) * (xc(0, 1) - xc(2, 1))) /
                           pow(vec1sq * vec2sq, 1.5)) +
                          (3 * pow(-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0], 2) *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              (4. * pow(vec1sq * vec2sq, 2.5)) -
                          ((2 * vec1sq + 2 * vec2sq + 8 * vec1[1] * vec2[0]) *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              (2. * pow(vec1sq * vec2sq, 1.5))) /
                        sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                         xc(1, 0) * xc(2, 1),
                                     2) /
                                     (vec1sq * vec2sq))) +
                    (((xc(0, 1) - xc(2, 1)) / sqrt(vec1sq * vec2sq) -
                         ((-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0]) *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                             (2. * pow(vec1sq * vec2sq, 1.5))) *
                        ((-2 * (xc(0, 1) - xc(2, 1)) *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                                (vec1sq * vec2sq) -
                            (2 * vec1[1] *
                                pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                        xc(1, 0) * xc(2, 1),
                                    2)) /
                                (pow(vec1sq, 2) * vec2sq) -
                            (2 * vec2[0] *
                                pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                        xc(1, 0) * xc(2, 1),
                                    2)) /
                                (vec1sq * pow(vec2sq, 2)))) /
                        (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                              xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1),
                                          2) /
                                          (vec1sq * vec2sq),
                                  1.5));
  elematrix(2, 3) = -((-((-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0]) * (-xc(0, 0) + xc(2, 0))) /
                              (2. * pow(vec1sq * vec2sq, 1.5)) -
                          ((2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1]) * (xc(0, 1) - xc(2, 1))) /
                              (2. * pow(vec1sq * vec2sq, 1.5)) +
                          (3 * (-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0]) *
                              (2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1]) *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              (4. * pow(vec1sq * vec2sq, 2.5)) -
                          ((-4 * vec1[0] * vec2[0] + 4 * vec1[1] * vec2[1]) *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              (2. * pow(vec1sq * vec2sq, 1.5))) /
                        sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                         xc(1, 0) * xc(2, 1),
                                     2) /
                                     (vec1sq * vec2sq))) +
                    (((xc(0, 1) - xc(2, 1)) / sqrt(vec1sq * vec2sq) -
                         ((-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0]) *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                             (2. * pow(vec1sq * vec2sq, 1.5))) *
                        ((-2 * (-xc(0, 0) + xc(2, 0)) *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                                (vec1sq * vec2sq) +
                            (2 * vec1[0] *
                                pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                        xc(1, 0) * xc(2, 1),
                                    2)) /
                                (pow(vec1sq, 2) * vec2sq) -
                            (2 * vec2[1] *
                                pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                        xc(1, 0) * xc(2, 1),
                                    2)) /
                                (vec1sq * pow(vec2sq, 2)))) /
                        (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                              xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1),
                                          2) /
                                          (vec1sq * vec2sq),
                                  1.5));
  elematrix(2, 4) = -((-((-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0]) * (-xc(0, 1) + xc(1, 1))) /
                              (2. * pow(vec1sq * vec2sq, 1.5)) -
                          (vec1sq * vec2[0] * (xc(0, 1) - xc(2, 1))) / pow(vec1sq * vec2sq, 1.5) +
                          (3 * vec1sq * vec2[0] * (-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0]) *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              (2. * pow(vec1sq * vec2sq, 2.5)) -
                          ((-2 * vec1sq - 4 * vec1[1] * vec2[0]) *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              (2. * pow(vec1sq * vec2sq, 1.5))) /
                        sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                         xc(1, 0) * xc(2, 1),
                                     2) /
                                     (vec1sq * vec2sq))) +
                    (((xc(0, 1) - xc(2, 1)) / sqrt(vec1sq * vec2sq) -
                         ((-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0]) *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                             (2. * pow(vec1sq * vec2sq, 1.5))) *
                        ((-2 * (-xc(0, 1) + xc(1, 1)) *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                                (vec1sq * vec2sq) +
                            (2 * vec2[0] *
                                pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                        xc(1, 0) * xc(2, 1),
                                    2)) /
                                (vec1sq * pow(vec2sq, 2)))) /
                        (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                              xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1),
                                          2) /
                                          (vec1sq * vec2sq),
                                  1.5));
  elematrix(2, 5) = -((-(1 / sqrt(vec1sq * vec2sq)) -
                          (vec1[1] * (-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0])) /
                              (2. * pow(vec1sq * vec2sq, 1.5)) -
                          (vec1sq * vec2[1] * (xc(0, 1) - xc(2, 1))) / pow(vec1sq * vec2sq, 1.5) +
                          (2 * vec1[1] * vec2(1) *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              pow(vec1sq * vec2sq, 1.5) +
                          (3 * vec1sq * (-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0]) * vec2(1) *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              (2. * pow(vec1sq * vec2sq, 2.5))) /
                        sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                         xc(1, 0) * xc(2, 1),
                                     2) /
                                     (vec1sq * vec2sq))) +
                    (((xc(0, 1) - xc(2, 1)) / sqrt(vec1sq * vec2sq) -
                         ((-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0]) *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                             (2. * pow(vec1sq * vec2sq, 1.5))) *
                        ((-2 * vec1[1] *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                                (vec1sq * vec2sq) +
                            (2 * vec2[1] *
                                pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                        xc(1, 0) * xc(2, 1),
                                    2)) /
                                (vec1sq * pow(vec2sq, 2)))) /
                        (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                              xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1),
                                          2) /
                                          (vec1sq * vec2sq),
                                  1.5));
  elematrix(3, 0) = -((-(1 / sqrt(vec1sq * vec2sq)) -
                          (vec2[1] * (2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1])) /
                              (2. * pow(vec1sq * vec2sq, 1.5)) -
                          (vec2sq * vec1[1] * (-xc(0, 0) + xc(2, 0))) / pow(vec1sq * vec2sq, 1.5) +
                          (2 * vec1[1] * vec2(1) *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              pow(vec1sq * vec2sq, 1.5) +
                          (3 * vec2sq * vec1[1] * (2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1]) *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              (2. * pow(vec1sq * vec2sq, 2.5))) /
                        sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                         xc(1, 0) * xc(2, 1),
                                     2) /
                                     (vec1sq * vec2sq))) +
                    (((-xc(0, 0) + xc(2, 0)) / sqrt(vec1sq * vec2sq) -
                         ((2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1]) *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                             (2. * pow(vec1sq * vec2sq, 1.5))) *
                        ((-2 * vec2[1] *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                                (vec1sq * vec2sq) +
                            (2 * vec1[1] *
                                pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                        xc(1, 0) * xc(2, 1),
                                    2)) /
                                (pow(vec1sq, 2) * vec2sq))) /
                        (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                              xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1),
                                          2) /
                                          (vec1sq * vec2sq),
                                  1.5));
  elematrix(3, 1) = -(((vec2[0] * (2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1])) /
                              (2. * pow(vec1sq * vec2sq, 1.5)) +
                          (vec2sq * vec1[0] * (-xc(0, 0) + xc(2, 0))) / pow(vec1sq * vec2sq, 1.5) -
                          (3 * vec2sq * vec1[0] * (2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1]) *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              (2. * pow(vec1sq * vec2sq, 2.5)) -
                          ((-2 * vec2sq + 4 * vec1[0] * vec2[1]) *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              (2. * pow(vec1sq * vec2sq, 1.5))) /
                        sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                         xc(1, 0) * xc(2, 1),
                                     2) /
                                     (vec1sq * vec2sq))) +
                    (((-xc(0, 0) + xc(2, 0)) / sqrt(vec1sq * vec2sq) -
                         ((2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1]) *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                             (2. * pow(vec1sq * vec2sq, 1.5))) *
                        ((2 * vec2[0] *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                                (vec1sq * vec2sq) -
                            (2 * vec1[0] *
                                pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                        xc(1, 0) * xc(2, 1),
                                    2)) /
                                (pow(vec1sq, 2) * vec2sq))) /
                        (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                              xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1),
                                          2) /
                                          (vec1sq * vec2sq),
                                  1.5));
  elematrix(3, 2) = -((-((-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0]) * (-xc(0, 0) + xc(2, 0))) /
                              (2. * pow(vec1sq * vec2sq, 1.5)) -
                          ((2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1]) * (xc(0, 1) - xc(2, 1))) /
                              (2. * pow(vec1sq * vec2sq, 1.5)) +
                          (3 * (-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0]) *
                              (2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1]) *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              (4. * pow(vec1sq * vec2sq, 2.5)) -
                          ((-4 * vec1[0] * vec2[0] + 4 * vec1[1] * vec2[1]) *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              (2. * pow(vec1sq * vec2sq, 1.5))) /
                        sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                         xc(1, 0) * xc(2, 1),
                                     2) /
                                     (vec1sq * vec2sq))) +
                    (((-xc(0, 0) + xc(2, 0)) / sqrt(vec1sq * vec2sq) -
                         ((2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1]) *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                             (2. * pow(vec1sq * vec2sq, 1.5))) *
                        ((-2 * (xc(0, 1) - xc(2, 1)) *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                                (vec1sq * vec2sq) -
                            (2 * vec1[1] *
                                pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                        xc(1, 0) * xc(2, 1),
                                    2)) /
                                (pow(vec1sq, 2) * vec2sq) -
                            (2 * vec2[0] *
                                pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                        xc(1, 0) * xc(2, 1),
                                    2)) /
                                (vec1sq * pow(vec2sq, 2)))) /
                        (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                              xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1),
                                          2) /
                                          (vec1sq * vec2sq),
                                  1.5));
  elematrix(3, 3) = -((-(((2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1]) * (-xc(0, 0) + xc(2, 0))) /
                           pow(vec1sq * vec2sq, 1.5)) +
                          (3 * pow(2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1], 2) *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              (4. * pow(vec1sq * vec2sq, 2.5)) -
                          ((2 * vec1sq + 2 * vec2sq - 8 * vec1[0] * vec2[1]) *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              (2. * pow(vec1sq * vec2sq, 1.5))) /
                        sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                         xc(1, 0) * xc(2, 1),
                                     2) /
                                     (vec1sq * vec2sq))) +
                    (((-xc(0, 0) + xc(2, 0)) / sqrt(vec1sq * vec2sq) -
                         ((2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1]) *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                             (2. * pow(vec1sq * vec2sq, 1.5))) *
                        ((-2 * (-xc(0, 0) + xc(2, 0)) *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                                (vec1sq * vec2sq) +
                            (2 * vec1[0] *
                                pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                        xc(1, 0) * xc(2, 1),
                                    2)) /
                                (pow(vec1sq, 2) * vec2sq) -
                            (2 * vec2[1] *
                                pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                        xc(1, 0) * xc(2, 1),
                                    2)) /
                                (vec1sq * pow(vec2sq, 2)))) /
                        (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                              xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1),
                                          2) /
                                          (vec1sq * vec2sq),
                                  1.5));
  elematrix(3, 4) = -((1 / sqrt(vec1sq * vec2sq) -
                          ((2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1]) * (-xc(0, 1) + xc(1, 1))) /
                              (2. * pow(vec1sq * vec2sq, 1.5)) -
                          (vec1sq * vec2[0] * (-xc(0, 0) + xc(2, 0))) / pow(vec1sq * vec2sq, 1.5) -
                          (2 * vec1[0] * vec2[0] *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              pow(vec1sq * vec2sq, 1.5) +
                          (3 * vec1sq * vec2[0] * (2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1]) *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              (2. * pow(vec1sq * vec2sq, 2.5))) /
                        sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                         xc(1, 0) * xc(2, 1),
                                     2) /
                                     (vec1sq * vec2sq))) +
                    (((-xc(0, 0) + xc(2, 0)) / sqrt(vec1sq * vec2sq) -
                         ((2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1]) *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                             (2. * pow(vec1sq * vec2sq, 1.5))) *
                        ((-2 * (-xc(0, 1) + xc(1, 1)) *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                                (vec1sq * vec2sq) +
                            (2 * vec2[0] *
                                pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                        xc(1, 0) * xc(2, 1),
                                    2)) /
                                (vec1sq * pow(vec2sq, 2)))) /
                        (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                              xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1),
                                          2) /
                                          (vec1sq * vec2sq),
                                  1.5));
  elematrix(3, 5) = -((-(vec1[1] * (2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1])) /
                              (2. * pow(vec1sq * vec2sq, 1.5)) -
                          (vec1sq * vec2[1] * (-xc(0, 0) + xc(2, 0))) / pow(vec1sq * vec2sq, 1.5) +
                          (3 * vec1sq * vec2[1] * (2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1]) *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              (2. * pow(vec1sq * vec2sq, 2.5)) -
                          ((-2 * vec1sq + 4 * vec1[0] * vec2[1]) *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              (2. * pow(vec1sq * vec2sq, 1.5))) /
                        sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                         xc(1, 0) * xc(2, 1),
                                     2) /
                                     (vec1sq * vec2sq))) +
                    (((-xc(0, 0) + xc(2, 0)) / sqrt(vec1sq * vec2sq) -
                         ((2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1]) *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                             (2. * pow(vec1sq * vec2sq, 1.5))) *
                        ((-2 * vec1[1] *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                                (vec1sq * vec2sq) +
                            (2 * vec2[1] *
                                pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                        xc(1, 0) * xc(2, 1),
                                    2)) /
                                (vec1sq * pow(vec2sq, 2)))) /
                        (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                              xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1),
                                          2) /
                                          (vec1sq * vec2sq),
                                  1.5));
  elematrix(4, 0) = -((-((vec1sq * vec2[0] * vec2[1]) / pow(vec1sq * vec2sq, 1.5)) -
                          (vec2sq * vec1[1] * (-xc(0, 1) + xc(1, 1))) / pow(vec1sq * vec2sq, 1.5) +
                          (3 * vec1sq * vec2sq * vec1[1] * vec2[0] *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              pow(vec1sq * vec2sq, 2.5) -
                          (2 * vec1[1] * vec2[0] *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              pow(vec1sq * vec2sq, 1.5)) /
                        sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                         xc(1, 0) * xc(2, 1),
                                     2) /
                                     (vec1sq * vec2sq))) +
                    (((-xc(0, 1) + xc(1, 1)) / sqrt(vec1sq * vec2sq) -
                         (vec1sq * vec2[0] *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                             pow(vec1sq * vec2sq, 1.5)) *
                        ((-2 * vec2[1] *
                             (vec2(1) * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                                (vec1sq * vec2sq) +
                            (2 * vec1[1] *
                                pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                        xc(1, 0) * xc(2, 1),
                                    2)) /
                                (pow(vec1sq, 2) * vec2sq))) /
                        (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                              xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1),
                                          2) /
                                          (vec1sq * vec2sq),
                                  1.5));
  elematrix(4, 1) =
      -((-(1 / sqrt(vec1sq * vec2sq)) + (vec1sq * pow(vec2[0], 2)) / pow(vec1sq * vec2sq, 1.5) +
            (vec2sq * vec1[0] * (-xc(0, 1) + xc(1, 1))) / pow(vec1sq * vec2sq, 1.5) -
            (3 * vec1sq * vec2sq * vec1[0] * vec2[0] *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 2.5) +
            (2 * vec1[0] * vec2[0] *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 1.5)) /
          sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                           xc(1, 0) * xc(2, 1),
                       2) /
                       (vec1sq * vec2sq))) +
      (((-xc(0, 1) + xc(1, 1)) / sqrt(vec1sq * vec2sq) -
           (vec1sq * vec2[0] *
               (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                   xc(1, 0) * xc(2, 1))) /
               pow(vec1sq * vec2sq, 1.5)) *
          ((2 * vec2[0] *
               (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                   xc(1, 0) * xc(2, 1))) /
                  (vec1sq * vec2sq) -
              (2 * vec1[0] *
                  pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1),
                      2)) /
                  (pow(vec1sq, 2) * vec2sq))) /
          (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                xc(1, 0) * xc(2, 1),
                            2) /
                            (vec1sq * vec2sq),
                    1.5));
  elematrix(4, 2) = -((-((-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0]) * (-xc(0, 1) + xc(1, 1))) /
                              (2. * pow(vec1sq * vec2sq, 1.5)) -
                          (vec1sq * vec2[0] * (xc(0, 1) - xc(2, 1))) / pow(vec1sq * vec2sq, 1.5) +
                          (vec1sq * (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                        xc(1, 0) * xc(2, 1))) /
                              pow(vec1sq * vec2sq, 1.5) +
                          (2 * vec1[1] * vec2[0] *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              pow(vec1sq * vec2sq, 1.5) +
                          (3 * vec1sq * vec2[0] * (-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0]) *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              (2. * pow(vec1sq * vec2sq, 2.5))) /
                        sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                         xc(1, 0) * xc(2, 1),
                                     2) /
                                     (vec1sq * vec2sq))) +
                    (((-xc(0, 1) + xc(1, 1)) / sqrt(vec1sq * vec2sq) -
                         (vec1sq * vec2[0] *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                             pow(vec1sq * vec2sq, 1.5)) *
                        ((-2 * (xc(0, 1) - xc(2, 1)) *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                                (vec1sq * vec2sq) -
                            (2 * vec1[1] *
                                pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                        xc(1, 0) * xc(2, 1),
                                    2)) /
                                (pow(vec1sq, 2) * vec2sq) -
                            (2 * vec2[0] *
                                pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                        xc(1, 0) * xc(2, 1),
                                    2)) /
                                (vec1sq * pow(vec2sq, 2)))) /
                        (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                              xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1),
                                          2) /
                                          (vec1sq * vec2sq),
                                  1.5));
  elematrix(4, 3) = -((1 / sqrt(vec1sq * vec2sq) -
                          ((2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1]) * (-xc(0, 1) + xc(1, 1))) /
                              (2. * pow(vec1sq * vec2sq, 1.5)) -
                          (vec1sq * vec2[0] * (-xc(0, 0) + xc(2, 0))) / pow(vec1sq * vec2sq, 1.5) -
                          (2 * vec1[0] * vec2[0] *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              pow(vec1sq * vec2sq, 1.5) +
                          (3 * vec1sq * vec2[0] * (2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1]) *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              (2. * pow(vec1sq * vec2sq, 2.5))) /
                        sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                         xc(1, 0) * xc(2, 1),
                                     2) /
                                     (vec1sq * vec2sq))) +
                    (((-xc(0, 1) + xc(1, 1)) / sqrt(vec1sq * vec2sq) -
                         (vec1sq * vec2[0] *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                             pow(vec1sq * vec2sq, 1.5)) *
                        ((-2 * (-xc(0, 0) + xc(2, 0)) *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                                (vec1sq * vec2sq) +
                            (2 * vec1[0] *
                                pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                        xc(1, 0) * xc(2, 1),
                                    2)) /
                                (pow(vec1sq, 2) * vec2sq) -
                            (2 * vec2[1] *
                                pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                        xc(1, 0) * xc(2, 1),
                                    2)) /
                                (vec1sq * pow(vec2sq, 2)))) /
                        (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                              xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1),
                                          2) /
                                          (vec1sq * vec2sq),
                                  1.5));
  elematrix(4, 4) =
      -(((-2 * vec1sq * vec2[0] * (-xc(0, 1) + xc(1, 1))) / pow(vec1sq * vec2sq, 1.5) -
            (vec1sq * (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 1.5) +
            (3 * pow(vec1sq, 2) * pow(vec2[0], 2) *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 2.5)) /
          sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                           xc(1, 0) * xc(2, 1),
                       2) /
                       (vec1sq * vec2sq))) +
      (((-xc(0, 1) + xc(1, 1)) / sqrt(vec1sq * vec2sq) -
           (vec1sq * vec2[0] *
               (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                   xc(1, 0) * xc(2, 1))) /
               pow(vec1sq * vec2sq, 1.5)) *
          ((-2 * (-xc(0, 1) + xc(1, 1)) *
               (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                   xc(1, 0) * xc(2, 1))) /
                  (vec1sq * vec2sq) +
              (2 * vec2[0] *
                  pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1),
                      2)) /
                  (vec1sq * pow(vec2sq, 2)))) /
          (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                xc(1, 0) * xc(2, 1),
                            2) /
                            (vec1sq * vec2sq),
                    1.5));
  elematrix(4, 5) = -((-((vec1sq * vec1[1] * vec2[0]) / pow(vec1sq * vec2sq, 1.5)) -
                          (vec1sq * vec2[1] * (-xc(0, 1) + xc(1, 1))) / pow(vec1sq * vec2sq, 1.5) +
                          (3 * pow(vec1sq, 2) * vec2[0] * vec2(1) *
                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                  xc(1, 0) * xc(2, 1))) /
                              pow(vec1sq * vec2sq, 2.5)) /
                        sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                         xc(1, 0) * xc(2, 1),
                                     2) /
                                     (vec1sq * vec2sq))) +
                    (((-xc(0, 1) + xc(1, 1)) / sqrt(vec1sq * vec2sq) -
                         (vec1sq * vec2[0] *
                             (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                             pow(vec1sq * vec2sq, 1.5)) *
                        ((-2 * vec1[1] *
                             (vec2(1) * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                 xc(1, 0) * xc(2, 1))) /
                                (vec1sq * vec2sq) +
                            (2 * vec2[1] *
                                pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                        xc(1, 0) * xc(2, 1),
                                    2)) /
                                (vec1sq * pow(vec2sq, 2)))) /
                        (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                              xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1),
                                          2) /
                                          (vec1sq * vec2sq),
                                  1.5));
  elematrix(5, 0) =
      -((1 / sqrt(vec1sq * vec2sq) - (vec2sq * pow(vec1[1], 2)) / pow(vec1sq * vec2sq, 1.5) -
            (vec1sq * pow(vec2[1], 2)) / pow(vec1sq * vec2sq, 1.5) +
            (3 * vec1sq * vec2sq * vec1[1] * vec2(1) *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 2.5) -
            (2 * vec1[1] * vec2(1) *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 1.5)) /
          sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                           xc(1, 0) * xc(2, 1),
                       2) /
                       (vec1sq * vec2sq))) +
      ((vec1[1] / sqrt(vec1sq * vec2sq) - (vec1sq * vec2(1) *
                                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                                  xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1))) /
                                              pow(vec1sq * vec2sq, 1.5)) *
          ((-2 * vec2[1] *
               (vec2(1) * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                   xc(1, 0) * xc(2, 1))) /
                  (vec1sq * vec2sq) +
              (2 * vec1[1] *
                  pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1),
                      2)) /
                  (pow(vec1sq, 2) * vec2sq))) /
          (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                xc(1, 0) * xc(2, 1),
                            2) /
                            (vec1sq * vec2sq),
                    1.5));
  elematrix(5, 1) =
      -(((vec2sq * vec1[0] * vec1[1]) / pow(vec1sq * vec2sq, 1.5) +
            (vec1sq * vec2[0] * vec2[1]) / pow(vec1sq * vec2sq, 1.5) -
            (3 * vec1sq * vec2sq * vec1[0] * vec2(1) *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 2.5) +
            (2 * vec1[0] * vec2(1) *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 1.5)) /
          sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                           xc(1, 0) * xc(2, 1),
                       2) /
                       (vec1sq * vec2sq))) +
      ((vec1[1] / sqrt(vec1sq * vec2sq) - (vec1sq * vec2(1) *
                                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                                  xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1))) /
                                              pow(vec1sq * vec2sq, 1.5)) *
          ((2 * vec2[0] *
               (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                   xc(1, 0) * xc(2, 1))) /
                  (vec1sq * vec2sq) -
              (2 * vec1[0] *
                  pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1),
                      2)) /
                  (pow(vec1sq, 2) * vec2sq))) /
          (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                xc(1, 0) * xc(2, 1),
                            2) /
                            (vec1sq * vec2sq),
                    1.5));
  elematrix(5, 2) =
      -((-(1 / sqrt(vec1sq * vec2sq)) -
            (vec1[1] * (-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0])) /
                (2. * pow(vec1sq * vec2sq, 1.5)) -
            (vec1sq * vec2[1] * (xc(0, 1) - xc(2, 1))) / pow(vec1sq * vec2sq, 1.5) +
            (2 * vec1[1] * vec2(1) *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 1.5) +
            (3 * vec1sq * (-2 * vec2sq * vec1[1] - 2 * vec1sq * vec2[0]) * vec2(1) *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                (2. * pow(vec1sq * vec2sq, 2.5))) /
          sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                           xc(1, 0) * xc(2, 1),
                       2) /
                       (vec1sq * vec2sq))) +
      ((vec1[1] / sqrt(vec1sq * vec2sq) - (vec1sq * vec2(1) *
                                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                                  xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1))) /
                                              pow(vec1sq * vec2sq, 1.5)) *
          ((-2 * (xc(0, 1) - xc(2, 1)) *
               (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                   xc(1, 0) * xc(2, 1))) /
                  (vec1sq * vec2sq) -
              (2 * vec1[1] *
                  pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1),
                      2)) /
                  (pow(vec1sq, 2) * vec2sq) -
              (2 * vec2[0] *
                  pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1),
                      2)) /
                  (vec1sq * pow(vec2sq, 2)))) /
          (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                xc(1, 0) * xc(2, 1),
                            2) /
                            (vec1sq * vec2sq),
                    1.5));
  elematrix(5, 3) =
      -((-(vec1[1] * (2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1])) /
                (2. * pow(vec1sq * vec2sq, 1.5)) -
            (vec1sq * vec2[1] * (-xc(0, 0) + xc(2, 0))) / pow(vec1sq * vec2sq, 1.5) +
            (vec1sq * (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 1.5) -
            (2 * vec1[0] * vec2(1) *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 1.5) +
            (3 * vec1sq * vec2[1] * (2 * vec2sq * vec1[0] - 2 * vec1sq * vec2[1]) *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                (2. * pow(vec1sq * vec2sq, 2.5))) /
          sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                           xc(1, 0) * xc(2, 1),
                       2) /
                       (vec1sq * vec2sq))) +
      ((vec1[1] / sqrt(vec1sq * vec2sq) - (vec1sq * vec2(1) *
                                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                                  xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1))) /
                                              pow(vec1sq * vec2sq, 1.5)) *
          ((-2 * (-xc(0, 0) + xc(2, 0)) *
               (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                   xc(1, 0) * xc(2, 1))) /
                  (vec1sq * vec2sq) +
              (2 * vec1[0] *
                  pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1),
                      2)) /
                  (pow(vec1sq, 2) * vec2sq) -
              (2 * vec2[1] *
                  pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1),
                      2)) /
                  (vec1sq * pow(vec2sq, 2)))) /
          (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                xc(1, 0) * xc(2, 1),
                            2) /
                            (vec1sq * vec2sq),
                    1.5));
  elematrix(5, 4) =
      -((-((vec1sq * vec1[1] * vec2[0]) / pow(vec1sq * vec2sq, 1.5)) -
            (vec1sq * vec2[1] * (-xc(0, 1) + xc(1, 1))) / pow(vec1sq * vec2sq, 1.5) +
            (3 * pow(vec1sq, 2) * vec2[0] * vec2(1) *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 2.5)) /
          sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                           xc(1, 0) * xc(2, 1),
                       2) /
                       (vec1sq * vec2sq))) +
      ((vec1[1] / sqrt(vec1sq * vec2sq) - (vec1sq * vec2(1) *
                                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                                  xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1))) /
                                              pow(vec1sq * vec2sq, 1.5)) *
          ((-2 * (-xc(0, 1) + xc(1, 1)) *
               (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                   xc(1, 0) * xc(2, 1))) /
                  (vec1sq * vec2sq) +
              (2 * vec2[0] *
                  pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1),
                      2)) /
                  (vec1sq * pow(vec2sq, 2)))) /
          (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                xc(1, 0) * xc(2, 1),
                            2) /
                            (vec1sq * vec2sq),
                    1.5));
  elematrix(5, 5) =
      -(((-2 * vec1sq * vec1[1] * vec2[1]) / pow(vec1sq * vec2sq, 1.5) -
            (vec1sq * (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 1.5) +
            (3 * pow(vec1sq, 2) * pow(vec2[1], 2) *
                (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                    xc(1, 0) * xc(2, 1))) /
                pow(vec1sq * vec2sq, 2.5)) /
          sqrt(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                           xc(1, 0) * xc(2, 1),
                       2) /
                       (vec1sq * vec2sq))) +
      ((vec1[1] / sqrt(vec1sq * vec2sq) - (vec1sq * vec2(1) *
                                              (vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) +
                                                  xc(1, 1) * xc(2, 0) - xc(1, 0) * xc(2, 1))) /
                                              pow(vec1sq * vec2sq, 1.5)) *
          ((-2 * vec1[1] *
               (vec2(1) * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                   xc(1, 0) * xc(2, 1))) /
                  (vec1sq * vec2sq) +
              (2 * vec2[1] *
                  pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                          xc(1, 0) * xc(2, 1),
                      2)) /
                  (vec1sq * pow(vec2sq, 2)))) /
          (2. * pow(1 - pow(vec2[1] * xc(0, 0) - vec2[0] * xc(0, 1) + xc(1, 1) * xc(2, 0) -
                                xc(1, 0) * xc(2, 1),
                            2) /
                            (vec1sq * vec2sq),
                    1.5));

  elematrix.Scale(-1.0);
}
