/*!----------------------------------------------------------------------
\file constraint_element_evaluate.cpp
\brief

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "constraint_element.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "Epetra_SerialDenseSolver.h"

using namespace DRT::UTILS;



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::ConstraintElement::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  ActionType act = none;

  // get the required action and distinguish between 2d and 3d MPC's
  string action = params.get<string>("action","none");
  if (action == "none") return 0;
  else if (action=="calc_MPC3D_stiff")
  {
    act=calc_MPC3D_stiff;
  }
  else if (action=="calc_MPC3D_state")
  {
  act=calc_MPC3D_state;
  }
  else if (action=="calc_MPC2D_stiff")
  {

    RCP<DRT::Condition> condition = params.get<RefCountPtr<DRT::Condition> >("condition");
    const string* type = condition->Get<string>("control value");

    if (*type == "dist") act = calc_MPC2D_dist_stiff;
    else if (*type == "angle") act = calc_MPC2D_angle_stiff;
    else dserror("No constraint type in 2d MPC specified. Value to control should by either be 'dist' or 'angle'!");

  }
  else
    dserror("Unknown type of action for ConstraintElement");

  switch (act)
  {
    case none:
    {
      return(0);
    }
    break;
    case calc_MPC3D_state:
    {
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get state vector 'displacement'");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      const int numnod = NumNode(); 

      if (numnod == 4)
      {
        const int numdim = 3;
        const int numnode = 4;
        LINALG::Matrix<numnode,numdim> xscurr;  // material coord. of element
        SpatialConfiguration(xscurr,mydisp);
        LINALG::Matrix<numdim,1> elementnormal;

        ComputeNormal3D(xscurr,elementnormal);
        if(abs(elementnormal.Norm2())<1E-6)
        {
          dserror("Bad plane, points almost on a line!");
        }
 
        elevec3[0] =ComputeNormalDist3D(xscurr,elementnormal);
      }
      else if (numnod == 2)
      {
        RCP<DRT::Condition> condition = params.get<RefCountPtr<DRT::Condition> >("condition");
        const vector<double>*  direct = condition->Get<vector<double> > ("direction");
               
        // compute difference
        elevec3[0] = ComputeWeightedDistance(mydisp,*direct);
      }
    }
    break;
    case calc_MPC3D_stiff:
    {
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get state vector 'displacement'");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      const int numnod = NumNode();
      
      if (numnod==4)
      {
        const int numdim = 3;
        const int numnode = 4;

        LINALG::Matrix<numnode,numdim> xscurr;  // material coord. of element
        SpatialConfiguration(xscurr,mydisp);
  
        LINALG::Matrix<numdim,1> elementnormal;
        ComputeNormal3D(xscurr,elementnormal);
        if(abs(elementnormal.Norm2())<1E-6)
        {
          dserror("Bad plane, points almost on a line!");
        }
        double normaldistance =ComputeNormalDist3D(xscurr,elementnormal);
  
        ComputeFirstDeriv3D(xscurr,elevec1,elementnormal);
        ComputeSecondDeriv3D(xscurr,elemat1,elementnormal);
  
        //update corresponding column in "constraint" matrix
        elevec2=elevec1;
        elevec3[0]=normaldistance;
      }
      else if (numnod == 2)
      {
        RCP<DRT::Condition> condition = params.get<RefCountPtr<DRT::Condition> >("condition");
        const vector<double>*  direct = condition->Get<vector<double> > ("direction");  
        
        //Compute weighted difference between masternode and other node and it's derivative
        ComputeFirstDerivWeightedDistance(elevec1,*direct);
        elevec3[0] = ComputeWeightedDistance(mydisp,*direct);
        elevec2=elevec1;
      }
      
    }
    break;
    case calc_MPC2D_dist_stiff:
    {
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get state vector 'displacement'");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      const int numnode = 3;
      const int numdim = 2;
      LINALG::Matrix<numnode,numdim> xscurr;  // material coord. of element
      SpatialConfiguration(xscurr,mydisp);
      LINALG::Matrix<numdim,1> elementnormal;
      ComputeNormal2D(xscurr,elementnormal);
      double normaldistance =ComputeNormalDist2D(xscurr,elementnormal);
      ComputeFirstDerivDist2D(xscurr,elevec1,elementnormal);
      ComputeSecondDerivDist2D(xscurr,elemat1,elementnormal);
      //update corresponding column in "constraint" matrix
      elevec2=elevec1;
      elevec3[0]=normaldistance;
    }
    break;
    case calc_MPC2D_angle_stiff:
    {
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get state vector 'displacement'");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      const int numnode = 3;
      const int numdim = 2;
      LINALG::Matrix<numnode,numdim> xscurr;  // material coord. of element
      SpatialConfiguration(xscurr,mydisp);

      double angle=ComputeAngle2D(xscurr);

      ComputeFirstDerivAngle2D(xscurr,elevec1);
      ComputeSecondDerivAngle2D(xscurr,elemat1);

      //update corresponding column in "constraint" matrix
      elevec2=elevec1;
      elevec3[0]=angle;

    }
    break;
    default:
      dserror("Unimplemented type of action");
  }
  return 0;


} // end of DRT::ELEMENTS::ConstraintElement::Evaluate

/*----------------------------------------------------------------------*
 * Evaluate Neumann (->dserror) */
int DRT::ELEMENTS::ConstraintElement::EvaluateNeumann
(
  ParameterList& params,
  DRT::Discretization&      discretization,
  DRT::Condition&           condition,
  vector<int>&              lm,
  Epetra_SerialDenseVector& elevec1
)
{
  dserror("You called Evaluate Neumann of constraint element.");
  return 0;
}

/*----------------------------------------------------------------------*
 * compute 3d normal */
void DRT::ELEMENTS::ConstraintElement::ComputeNormal3D
(
  const LINALG::Matrix<4,3>& xc,
  LINALG::Matrix<3,1>& elenorm
)
{
  elenorm(0,0)=-(xc(0,2)*xc(1,1)) + xc(0,1)*xc(1,2) + xc(0,2)*xc(2,1) -
    xc(1,2)*xc(2,1) - xc(0,1)*xc(2,2) + xc(1,1)*xc(2,2);
  elenorm(1,0)=xc(0,2)*xc(1,0) - xc(0,0)*xc(1,2) - xc(0,2)*xc(2,0) +
    xc(1,2)*xc(2,0) + xc(0,0)*xc(2,2) - xc(1,0)*xc(2,2);
  elenorm(2,0)=-(xc(0,1)*xc(1,0)) + xc(0,0)*xc(1,1) + xc(0,1)*xc(2,0) -
    xc(1,1)*xc(2,0) - xc(0,0)*xc(2,1) + xc(1,0)*xc(2,1);
  return ;
}

/*----------------------------------------------------------------------*
 * compute 2d normal */
void DRT::ELEMENTS::ConstraintElement::ComputeNormal2D
(
  const LINALG::Matrix<3,2>& xc,
  LINALG::Matrix<2,1>& elenorm
)
{
  elenorm(0,0)=xc(0,1) - xc(1,1);
  elenorm(1,0)=-xc(0,0) + xc(1,0);
  return ;
}

/*----------------------------------------------------------------------*
 * normal distance between fourth point and plane */
double DRT::ELEMENTS::ConstraintElement::ComputeNormalDist3D
(
  const LINALG::Matrix<4,3>& xc,
  const LINALG::Matrix<3,1>& normal
)
{
  return (-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*
      (xc(0,2) - xc(3,2)))/(-normal.Norm2());
}

/*----------------------------------------------------------------------*
 * normal distance between third point and line */
double DRT::ELEMENTS::ConstraintElement::ComputeNormalDist2D
(
  const LINALG::Matrix<3,2>& xc,
  const LINALG::Matrix<2,1>& normal
)
{
  return (normal(0,0)*(-xc(0,0) + xc(2,0)) - normal(1,0)*(xc(0,1) - xc(2,1)))/normal.Norm2();
}

/*----------------------------------------------------------------------*
 * compute angle at second point */
double DRT::ELEMENTS::ConstraintElement::ComputeAngle2D
(
  const LINALG::Matrix<3,2>& xc
)
{
  return (acos((xc(0,1)*(xc(1,0) - xc(2,0)) + xc(1,1)*xc(2,0) - xc(1,0)*xc(2,1) + xc(0,0)*(-xc(1,1) + xc(2,1)))/
       sqrt((pow(xc(0,0) - xc(1,0),2) + pow(xc(0,1) - xc(1,1),2))*(pow(xc(1,0) - xc(2,0),2) + pow(xc(1,1) - xc(2,1),2))))+acos(0.0));

}

/*----------------------------------------------------------------------*
 * first derivatives */
void DRT::ELEMENTS::ConstraintElement::ComputeFirstDeriv3D
(
  const LINALG::Matrix<4,3>& xc,
  Epetra_SerialDenseVector& elevector,
  const LINALG::Matrix<3,1>& normal
)
{

  double normsquare=pow(normal.Norm2(),2);
  double normcube=pow(normal.Norm2(),3);

  elevector[0]=
    (-((-2*normal(2,0)*(-xc(1,1) + xc(2,1)) + 2*normal(1,0)*(-xc(1,2) +
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2)))) + 2*normsquare*(xc(1,2)*(xc(2,1) - xc(3,1)) +
    xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) + xc(1,1)*(-xc(2,2) +
    xc(3,2))))/(2.*normcube);

  elevector[1]=
    (-((-2*(normal(2,0))*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(-xc(1,2) +
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2)))) + 2*normsquare*(-(xc(2,2)*xc(3,0)) +
    xc(1,2)*(-xc(2,0) + xc(3,0)) + xc(1,0)*(xc(2,2) - xc(3,2)) +
    xc(2,0)*xc(3,2)))/(2.*normcube);

  elevector[2]=
    (2*normsquare*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) - (2*normal(1,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(xc(1,1) - xc(2,1)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) -
    xc(3,2))))/(2.*normcube);

  elevector[3]=
    (-((-2*normal(2,0)*(xc(0,1) - xc(2,1)) + 2*normal(1,0)*(xc(0,2) -
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2)))) + 2*normsquare*((xc(0,2) - xc(2,2))*(-xc(0,1) +
    xc(3,1)) + (xc(0,1) - xc(2,1))*(xc(0,2) - xc(3,2))))/(2.*normcube);

  elevector[4]=
    (-((-2*normal(2,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(xc(0,2) -
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2)))) + 2*normsquare*((xc(0,2) - xc(2,2))*(xc(0,0) -
    xc(3,0)) + (-xc(0,0) + xc(2,0))*(xc(0,2) - xc(3,2))))/(2.*normcube);

  elevector[5]=
    (2*normsquare*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) - (2*normal(1,0)*(-xc(0,0) + xc(2,0)) -
    2*normal(0,0)*(-xc(0,1) + xc(2,1)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) -
    xc(3,2))))/(2.*normcube);

  elevector[6]=
    (-((-2*normal(2,0)*(-xc(0,1) + xc(1,1)) + 2*normal(1,0)*(-xc(0,2) +
    xc(1,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2)))) + 2*normsquare*((xc(0,2) - xc(1,2))*(xc(0,1) -
    xc(3,1)) + (-xc(0,1) + xc(1,1))*(xc(0,2) - xc(3,2))))/(2.*normcube);

  elevector[7]=
    (-((-2*normal(2,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(-xc(0,2) +
    xc(1,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2)))) + 2*normsquare*((-xc(0,2) + xc(1,2))*(xc(0,0) -
    xc(3,0)) + (xc(0,0) - xc(1,0))*(xc(0,2) - xc(3,2))))/(2.*normcube);

  elevector[8]=
    (2*normsquare*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(-xc(0,1) + xc(3,1))) - (2*normal(1,0)*(xc(0,0) - xc(1,0)) -
    2*normal(0,0)*(xc(0,1) - xc(1,1)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) -
    xc(3,2))))/(2.*normcube);

  elevector[9]=
    (-(xc(1,2)*xc(2,1)) + xc(0,2)*(-xc(1,1) + xc(2,1)) + xc(0,1)*(xc(1,2) - xc(2,2))
    + xc(1,1)*xc(2,2))/normal.Norm2();

  elevector[10]=
    normal(1,0)/normal.Norm2();

  elevector[11]=
    (-(xc(1,1)*xc(2,0)) + xc(0,1)*(-xc(1,0) + xc(2,0)) + xc(0,0)*(xc(1,1) - xc(2,1))
    + xc(1,0)*xc(2,1))/normal.Norm2();

  return;
}

/*----------------------------------------------------------------------*
 * second derivatives */
void DRT::ELEMENTS::ConstraintElement::ComputeFirstDerivDist2D
(
  const LINALG::Matrix<3,2>& xc,
  Epetra_SerialDenseVector& elevector, 
  const LINALG::Matrix<2,1>& normal
)
{
  double normcube=pow(normal.Norm2(),3);

  elevector[0]=(normal(0,0)*(-pow(xc(1,0),2) + xc(0,0)*(xc(1,0) - xc(2,0)) +
      xc(1,0)*xc(2,0) + normal(0,0)*(xc(1,1) - xc(2,1))))/normcube;

  elevector[1]=
    (normal(1,0)*(-pow(xc(1,0),2) + xc(0,0)*(xc(1,0) - xc(2,0)) + xc(1,0)*xc(2,0) +
    normal(0,0)*(xc(1,1) - xc(2,1))))/normcube;

  elevector[2]=
    -((normal(0,0)*(pow(xc(0,0),2) + pow(xc(0,1),2) + xc(1,0)*xc(2,0) -
    xc(0,0)*(xc(1,0) + xc(2,0)) + xc(1,1)*xc(2,1) - xc(0,1)*(xc(1,1) +
    xc(2,1))))/normcube);

  elevector[3]=
    -((normal(1,0)*(pow(xc(0,0),2) + pow(xc(0,1),2) + xc(1,0)*xc(2,0) -
    xc(0,0)*(xc(1,0) + xc(2,0)) + xc(1,1)*xc(2,1) - xc(0,1)*(xc(1,1) +
    xc(2,1))))/normcube);

  elevector[4]=normal(0,0)/normal.Norm2();

  elevector[5]=normal(1,0)/normal.Norm2();
  elevector.Scale(-1.0);
  return;
}

/*----------------------------------------------------------------------*
 * first derivatives */
void DRT::ELEMENTS::ConstraintElement::ComputeFirstDerivAngle2D
(
  const LINALG::Matrix<3,2>& xc,
  Epetra_SerialDenseVector& elevector
)
{
  LINALG::SerialDenseVector vec1(2);
  vec1[1]=xc(0,0) - xc(1,0);
  vec1[0]=-(xc(0,1) - xc(1,1));

  LINALG::SerialDenseVector vec2(2);
  vec2[0]=-xc(1,0) + xc(2,0);
  vec2[1]=-xc(1,1) + xc(2,1);

  const double vec1normsquare=pow(vec1.Norm2(),2);
  const double vec2normsquare=pow(vec2.Norm2(),2);

  elevector[0]
  =
  -((vec2[1]/sqrt(vec1normsquare*vec2normsquare) -
  (vec2normsquare*vec1[1]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1normsquare*vec2normsquare,1.5))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1normsquare*vec2normsquare)))
  ;
  elevector[1]
  =
  ((-(vec2[0]/sqrt(vec1normsquare*vec2normsquare)) +
  (vec2normsquare*vec1[0]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1normsquare*vec2normsquare,1.5))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1normsquare*vec2normsquare)))
  ;
  elevector[2]
  =
  (((xc(0,1) - xc(2,1))/sqrt(vec1normsquare*vec2normsquare) -
  ((-2*vec2normsquare*vec1[1] - 2*vec1normsquare*vec2[0])*(vec2[1]*xc(0,0) -
  vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1normsquare*vec2normsquare,1.5)))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1normsquare*vec2normsquare)))
  ;
  elevector[3]
  =
  (((-xc(0,0) + xc(2,0))/sqrt(vec1normsquare*vec2normsquare) -
  ((2*vec2normsquare*vec1[0] - 2*vec1normsquare*vec2[1])*(vec2[1]*xc(0,0) -
  vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1normsquare*vec2normsquare,1.5)))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1normsquare*vec2normsquare)))
  ;
  elevector[4]
  =
  (((-xc(0,1) + xc(1,1))/sqrt(vec1normsquare*vec2normsquare) -
  (vec1normsquare*vec2[0]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1normsquare*vec2normsquare,1.5))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1normsquare*vec2normsquare)))
  ;
  elevector[5]
  =
  ((vec1[1]/sqrt(vec1normsquare*vec2normsquare) -
  (vec1normsquare*vec2(1)*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1normsquare*vec2normsquare,1.5))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1normsquare*vec2normsquare)))
  ;
}

/*----------------------------------------------------------------------*
 * second derivatives */
void DRT::ELEMENTS::ConstraintElement::ComputeSecondDeriv3D
(
  const LINALG::Matrix<4,3>& xc,
  Epetra_SerialDenseMatrix& elematrix,
  const LINALG::Matrix<3,1>& normal
)
{

  double normsquare=pow(normal.Norm2(),2);
  double normcube=pow(normal.Norm2(),3);
  double normpowfour=pow(normal.Norm2(),4);
  double normpowfive=pow(normal.Norm2(),5);

  elematrix(0,0)=
    (-4*normsquare*(pow(xc(1,1) - xc(2,1),2) + pow(xc(1,2) -
    xc(2,2),2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) + 3*pow(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) +
    2*normal(1,0)*(-xc(1,2) + xc(2,2)),2)*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) + 2*normal(1,0)*(-xc(1,2) +
    xc(2,2)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) +
    xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(0,1)=
    (-4*normsquare*(xc(1,0) - xc(2,0))*(-xc(1,1) + xc(2,1))*(-(normal(0,0)*(xc(0,0) -
    xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) +
    3*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(-xc(1,2) +
    xc(2,2)))*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) + 2*normal(1,0)*(-xc(1,2) +
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(-xc(1,1) + xc(2,1))
    + 2*normal(1,0)*(-xc(1,2) + xc(2,2)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) +
    xc(3,0)) + xc(1,0)*(xc(2,2) - xc(3,2)) + xc(2,0)*xc(3,2)) -
    2*normsquare*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(-xc(1,2) +
    xc(2,2)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) +
    xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(0,2)=
    (-2*normsquare*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) + 2*normal(1,0)*(-xc(1,2) +
    xc(2,2)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) - 4*normsquare*(xc(1,0) - xc(2,0))*(-xc(1,2) +
    xc(2,2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) + 3*(2*normal(1,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(xc(1,1) - xc(2,1)))*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) +
    2*normal(1,0)*(-xc(1,2) + xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal(1,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(xc(1,1) -
    xc(2,1)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) +
    xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(0,3)=
    (3*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) + 2*normal(1,0)*(xc(0,2) -
    xc(2,2)))*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) + 2*normal(1,0)*(-xc(1,2) +
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*(xc(0,1) - xc(2,1))*(-xc(1,1) +
    xc(2,1)) + 2*(xc(0,2) - xc(2,2))*(-xc(1,2) + xc(2,2)))*(-(normal(0,0)*(xc(0,0) -
    xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) + 2*normal(1,0)*(-xc(1,2) +
    xc(2,2)))*((xc(0,2) - xc(2,2))*(-xc(0,1) + xc(3,1)) + (xc(0,1) -
    xc(2,1))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) +
    2*normal(1,0)*(xc(0,2) - xc(2,2)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1)
    - xc(2,1)*xc(3,2) + xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(0,4)=
    (-4*normsquare*(-2*xc(1,1)*xc(2,0) + xc(0,1)*(-xc(1,0) + xc(2,0)) +
    2*xc(0,0)*(xc(1,1) - xc(2,1)) + xc(1,0)*xc(2,1) +
    xc(2,0)*xc(2,1))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) +
    xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) + 3*(-2*normal(2,0)*(-xc(0,0) + xc(2,0))
    - 2*normal(0,0)*(xc(0,2) - xc(2,2)))*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) +
    2*normal(1,0)*(-xc(1,2) + xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) + 2*normal(1,0)*(-xc(1,2) +
    xc(2,2)))*((xc(0,2) - xc(2,2))*(xc(0,0) - xc(3,0)) + (-xc(0,0) +
    xc(2,0))*(xc(0,2) - xc(3,2))) + 4*normpowfour*(-xc(2,2) + xc(3,2)) -
    2*normsquare*(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(xc(0,2) -
    xc(2,2)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) +
    xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(0,5)=
    (-2*normsquare*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) + 2*normal(1,0)*(-xc(1,2) +
    xc(2,2)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) + 4*normpowfour*(xc(2,1) - xc(3,1)) -
    4*normsquare*(-2*xc(1,2)*xc(2,0) + xc(0,2)*(-xc(1,0) + xc(2,0)) +
    2*xc(0,0)*(xc(1,2) - xc(2,2)) + xc(1,0)*xc(2,2) +
    xc(2,0)*xc(2,2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) +
    xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) + 3*(2*normal(1,0)*(-xc(0,0) + xc(2,0))
    - 2*normal(0,0)*(-xc(0,1) + xc(2,1)))*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) +
    2*normal(1,0)*(-xc(1,2) + xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(-xc(0,1) +
    xc(2,1)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) +
    xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(0,6)=
    (-2*normsquare*(2*(xc(0,1) - xc(1,1))*(xc(1,1) - xc(2,1)) + 2*(xc(0,2) -
    xc(1,2))*(xc(1,2) - xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) +
    3*(-2*normal(2,0)*(-xc(0,1) + xc(1,1)) + 2*normal(1,0)*(-xc(0,2) +
    xc(1,2)))*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) + 2*normal(1,0)*(-xc(1,2) +
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(-xc(1,1) + xc(2,1))
    + 2*normal(1,0)*(-xc(1,2) + xc(2,2)))*((xc(0,2) - xc(1,2))*(xc(0,1) - xc(3,1)) +
    (-xc(0,1) + xc(1,1))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(-xc(0,1)
    + xc(1,1)) + 2*normal(1,0)*(-xc(0,2) + xc(1,2)))*(xc(1,2)*(xc(2,1) - xc(3,1)) +
    xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) + xc(1,1)*(-xc(2,2) +
    xc(3,2))))/(4.*normpowfive);

  elematrix(0,7)=
    (-4*normsquare*(xc(1,0)*xc(1,1) + xc(0,1)*(xc(1,0) - xc(2,0)) + xc(1,1)*xc(2,0)
    - 2*xc(0,0)*(xc(1,1) - xc(2,1)) - 2*xc(1,0)*xc(2,1))*(-(normal(0,0)*(xc(0,0) -
    xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) +
    3*(-2*normal(2,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(-xc(0,2) +
    xc(1,2)))*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) + 2*normal(1,0)*(-xc(1,2) +
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(-xc(1,1) + xc(2,1))
    + 2*normal(1,0)*(-xc(1,2) + xc(2,2)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) +
    (xc(0,0) - xc(1,0))*(xc(0,2) - xc(3,2))) + 4*normpowfour*(xc(1,2) -
    xc(3,2)) - 2*normsquare*(-2*normal(2,0)*(xc(0,0) - xc(1,0)) -
    2*normal(0,0)*(-xc(0,2) + xc(1,2)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1)
    - xc(2,1)*xc(3,2) + xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(0,8)=
    (4*normpowfour*(-xc(1,1) + xc(3,1)) -
    2*normsquare*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) + 2*normal(1,0)*(-xc(1,2) +
    xc(2,2)))*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(-xc(0,1) + xc(3,1))) - 4*normsquare*(xc(1,0)*xc(1,2) +
    xc(0,2)*(xc(1,0) - xc(2,0)) + xc(1,2)*xc(2,0) - 2*xc(0,0)*(xc(1,2) - xc(2,2)) -
    2*xc(1,0)*xc(2,2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) +
    xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) + 3*(2*normal(1,0)*(xc(0,0) - xc(1,0)) -
    2*normal(0,0)*(xc(0,1) - xc(1,1)))*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) +
    2*normal(1,0)*(-xc(1,2) + xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal(1,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(xc(0,1) -
    xc(1,1)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) +
    xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(0,9)=
    -((-(xc(1,2)*xc(2,1)) + xc(0,2)*(-xc(1,1) + xc(2,1)) + xc(0,1)*(xc(1,2) -
    xc(2,2)) + xc(1,1)*xc(2,2))*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) +
    2*normal(1,0)*(-xc(1,2) + xc(2,2))))/(2.*normcube);

  elematrix(0,10)=
    (normal(0,0)*(xc(0,2)*xc(1,1)*xc(1,2) - xc(1,0)*xc(1,1)*xc(2,0) +
    xc(1,1)*pow(xc(2,0),2) + xc(0,0)*(xc(1,0) - xc(2,0))*(xc(1,1) - xc(2,1)) +
    pow(xc(1,0),2)*xc(2,1) - xc(0,2)*xc(1,2)*xc(2,1) + pow(xc(1,2),2)*xc(2,1) -
    xc(1,0)*xc(2,0)*xc(2,1) - xc(0,2)*xc(1,1)*xc(2,2) - xc(1,1)*xc(1,2)*xc(2,2) +
    xc(0,2)*xc(2,1)*xc(2,2) - xc(1,2)*xc(2,1)*xc(2,2) + xc(1,1)*pow(xc(2,2),2) -
    xc(0,1)*(pow(xc(1,0),2) + pow(xc(1,2),2) - 2*xc(1,0)*xc(2,0) +
    pow(xc(2,0),2) - 2*xc(1,2)*xc(2,2) + pow(xc(2,2),2))))/normcube;

  elematrix(0,11)=
    -((normal(0,0)*(-(xc(0,1)*xc(1,1)*xc(1,2)) + xc(1,0)*xc(1,2)*xc(2,0) -
    xc(1,2)*pow(xc(2,0),2) + xc(0,1)*xc(1,2)*xc(2,1) + xc(1,1)*xc(1,2)*xc(2,1) -
    xc(1,2)*pow(xc(2,1),2) + xc(0,2)*(pow(xc(1,0),2) + pow(xc(1,1),2) -
    2*xc(1,0)*xc(2,0) + pow(xc(2,0),2) - 2*xc(1,1)*xc(2,1) + pow(xc(2,1),2)) -
    xc(0,0)*(xc(1,0) - xc(2,0))*(xc(1,2) - xc(2,2)) - pow(xc(1,0),2)*xc(2,2) +
    xc(0,1)*xc(1,1)*xc(2,2) - pow(xc(1,1),2)*xc(2,2) + xc(1,0)*xc(2,0)*xc(2,2) -
    xc(0,1)*xc(2,1)*xc(2,2) + xc(1,1)*xc(2,1)*xc(2,2)))/normcube);

  elematrix(1,0)=
    (-4*normsquare*(xc(1,0) - xc(2,0))*(-xc(1,1) + xc(2,1))*(-(normal(0,0)*(xc(0,0) -
    xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) +
    3*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(-xc(1,2) +
    xc(2,2)))*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) + 2*normal(1,0)*(-xc(1,2) +
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(-xc(1,1) + xc(2,1))
    + 2*normal(1,0)*(-xc(1,2) + xc(2,2)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) +
    xc(3,0)) + xc(1,0)*(xc(2,2) - xc(3,2)) + xc(2,0)*xc(3,2)) -
    2*normsquare*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(-xc(1,2) +
    xc(2,2)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) +
    xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(1,1)=
    (-4*normsquare*(pow(xc(1,0) - xc(2,0),2) + pow(xc(1,2) -
    xc(2,2),2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) + 3*pow(-2*normal(2,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(-xc(1,2) + xc(2,2)),2)*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(-xc(1,2) +
    xc(2,2)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) + xc(3,0)) + xc(1,0)*(xc(2,2)
    - xc(3,2)) + xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(1,2)=
    (-2*normsquare*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(-xc(1,2) +
    xc(2,2)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) - 4*normsquare*(xc(1,1) - xc(2,1))*(-xc(1,2) +
    xc(2,2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) + 3*(2*normal(1,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(xc(1,1) - xc(2,1)))*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(-xc(1,2) + xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal(1,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(xc(1,1) -
    xc(2,1)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) + xc(3,0)) + xc(1,0)*(xc(2,2)
    - xc(3,2)) + xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(1,3)=
    (-4*normsquare*(2*xc(0,1)*(xc(1,0) - xc(2,0)) + xc(1,1)*xc(2,0) -
    2*xc(1,0)*xc(2,1) + xc(2,0)*xc(2,1) + xc(0,0)*(-xc(1,1) +
    xc(2,1)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) + 3*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) +
    2*normal(1,0)*(xc(0,2) - xc(2,2)))*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(-xc(1,2) + xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(-xc(1,2) +
    xc(2,2)))*((xc(0,2) - xc(2,2))*(-xc(0,1) + xc(3,1)) + (xc(0,1) -
    xc(2,1))*(xc(0,2) - xc(3,2))) + 4*normpowfour*(xc(2,2) - xc(3,2)) -
    2*normsquare*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) + 2*normal(1,0)*(xc(0,2) -
    xc(2,2)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) + xc(3,0)) + xc(1,0)*(xc(2,2)
    - xc(3,2)) + xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(1,4)=
    (3*(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(xc(0,2) -
    xc(2,2)))*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(-xc(1,2) +
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*(xc(1,0) - xc(2,0))*(-xc(0,0) +
    xc(2,0)) + 2*(xc(0,2) - xc(2,2))*(-xc(1,2) + xc(2,2)))*(-(normal(0,0)*(xc(0,0) -
    xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(-xc(1,2) +
    xc(2,2)))*((xc(0,2) - xc(2,2))*(xc(0,0) - xc(3,0)) + (-xc(0,0) +
    xc(2,0))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(-xc(0,0) + xc(2,0))
    - 2*normal(0,0)*(xc(0,2) - xc(2,2)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) +
    xc(3,0)) + xc(1,0)*(xc(2,2) - xc(3,2)) + xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(1,5)=
    (4*normpowfour*(-xc(2,0) + xc(3,0)) -
    2*normsquare*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(-xc(1,2) +
    xc(2,2)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) - 4*normsquare*(-2*xc(1,2)*xc(2,1) +
    xc(0,2)*(-xc(1,1) + xc(2,1)) + 2*xc(0,1)*(xc(1,2) - xc(2,2)) + xc(1,1)*xc(2,2) +
    xc(2,1)*xc(2,2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) +
    xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) + 3*(2*normal(1,0)*(-xc(0,0) + xc(2,0))
    - 2*normal(0,0)*(-xc(0,1) + xc(2,1)))*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(-xc(1,2) + xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(-xc(0,1) +
    xc(2,1)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) + xc(3,0)) + xc(1,0)*(xc(2,2)
    - xc(3,2)) + xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(1,6)=
    (-4*normsquare*(xc(1,0)*xc(1,1) - 2*xc(0,1)*(xc(1,0) - xc(2,0)) -
    2*xc(1,1)*xc(2,0) + xc(0,0)*(xc(1,1) - xc(2,1)) +
    xc(1,0)*xc(2,1))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) +
    xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) + 3*(-2*normal(2,0)*(-xc(0,1) + xc(1,1))
    + 2*normal(1,0)*(-xc(0,2) + xc(1,2)))*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(-xc(1,2) + xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(-xc(1,2) +
    xc(2,2)))*((xc(0,2) - xc(1,2))*(xc(0,1) - xc(3,1)) + (-xc(0,1) +
    xc(1,1))*(xc(0,2) - xc(3,2))) + 4*normpowfour*(-xc(1,2) + xc(3,2)) -
    2*normsquare*(-2*normal(2,0)*(-xc(0,1) + xc(1,1)) + 2*normal(1,0)*(-xc(0,2) +
    xc(1,2)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) + xc(3,0)) + xc(1,0)*(xc(2,2)
    - xc(3,2)) + xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(1,7)=
    (-2*normsquare*(2*(xc(0,0) - xc(1,0))*(xc(1,0) - xc(2,0)) + 2*(xc(0,2) -
    xc(1,2))*(xc(1,2) - xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) +
    3*(-2*normal(2,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(-xc(0,2) +
    xc(1,2)))*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(-xc(1,2) +
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(xc(1,0) - xc(2,0))
    - 2*normal(0,0)*(-xc(1,2) + xc(2,2)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) +
    (xc(0,0) - xc(1,0))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(xc(0,0) -
    xc(1,0)) - 2*normal(0,0)*(-xc(0,2) + xc(1,2)))*(-(xc(2,2)*xc(3,0)) +
    xc(1,2)*(-xc(2,0) + xc(3,0)) + xc(1,0)*(xc(2,2) - xc(3,2)) +
    xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(1,8)=
    (4*normpowfour*(xc(1,0) - xc(3,0)) - 2*normsquare*(-2*normal(2,0)*(xc(1,0)
    - xc(2,0)) - 2*normal(0,0)*(-xc(1,2) + xc(2,2)))*((xc(0,1) - xc(1,1))*(xc(0,0) -
    xc(3,0)) + (xc(0,0) - xc(1,0))*(-xc(0,1) + xc(3,1))) -
    4*normsquare*(xc(1,1)*xc(1,2) + xc(0,2)*(xc(1,1) - xc(2,1)) + xc(1,2)*xc(2,1) -
    2*xc(0,1)*(xc(1,2) - xc(2,2)) - 2*xc(1,1)*xc(2,2))*(-(normal(0,0)*(xc(0,0) -
    xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) +
    3*(2*normal(1,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(xc(0,1) -
    xc(1,1)))*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(-xc(1,2) +
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal(1,0)*(xc(0,0) - xc(1,0)) -
    2*normal(0,0)*(xc(0,1) - xc(1,1)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) +
    xc(3,0)) + xc(1,0)*(xc(2,2) - xc(3,2)) +xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(1,9)=
    (normal(1,0)*(xc(0,2)*xc(1,0)*xc(1,2) + pow(xc(1,1),2)*xc(2,0) -
    xc(0,2)*xc(1,2)*xc(2,0) + pow(xc(1,2),2)*xc(2,0) + xc(0,1)*(xc(1,0) -
    xc(2,0))*(xc(1,1) - xc(2,1)) - xc(1,0)*xc(1,1)*xc(2,1) - xc(1,1)*xc(2,0)*xc(2,1)
    + xc(1,0)*pow(xc(2,1),2) - xc(0,2)*xc(1,0)*xc(2,2) - xc(1,0)*xc(1,2)*xc(2,2) +
    xc(0,2)*xc(2,0)*xc(2,2) - xc(1,2)*xc(2,0)*xc(2,2) + xc(1,0)*pow(xc(2,2),2) -
    xc(0,0)*(pow(xc(1,1),2) + pow(xc(1,2),2) - 2*xc(1,1)*xc(2,1) +
    pow(xc(2,1),2) - 2*xc(1,2)*xc(2,2) + pow(xc(2,2),2))))/normcube;

  elematrix(1,10)=
    -(normal(1,0)*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(-xc(1,2) +
    xc(2,2))))/(2.*normcube);

  elematrix(1,11)=
    -((normal(1,0)*(-(xc(0,1)*xc(1,1)*xc(1,2)) + xc(1,0)*xc(1,2)*xc(2,0) -
    xc(1,2)*pow(xc(2,0),2) + xc(0,1)*xc(1,2)*xc(2,1) + xc(1,1)*xc(1,2)*xc(2,1) -
    xc(1,2)*pow(xc(2,1),2) + xc(0,2)*(pow(xc(1,0),2) + pow(xc(1,1),2) -
    2*xc(1,0)*xc(2,0) + pow(xc(2,0),2) - 2*xc(1,1)*xc(2,1) + pow(xc(2,1),2)) -
    xc(0,0)*(xc(1,0) - xc(2,0))*(xc(1,2) - xc(2,2)) - pow(xc(1,0),2)*xc(2,2) +
    xc(0,1)*xc(1,1)*xc(2,2) - pow(xc(1,1),2)*xc(2,2) + xc(1,0)*xc(2,0)*xc(2,2) -
    xc(0,1)*xc(2,1)*xc(2,2) + xc(1,1)*xc(2,1)*xc(2,2)))/normcube);

  elematrix(2,0)=
    (-2*normsquare*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) + 2*normal(1,0)*(-xc(1,2) +
    xc(2,2)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) - 4*normsquare*(xc(1,0) - xc(2,0))*(-xc(1,2) +
    xc(2,2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) + 3*(2*normal(1,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(xc(1,1) - xc(2,1)))*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) +
    2*normal(1,0)*(-xc(1,2) + xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal(1,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(xc(1,1) -
    xc(2,1)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) +
    xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(2,1)=
    (-2*normsquare*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(-xc(1,2) +
    xc(2,2)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) - 4*normsquare*(xc(1,1) - xc(2,1))*(-xc(1,2) +
    xc(2,2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) + 3*(2*normal(1,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(xc(1,1) - xc(2,1)))*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(-xc(1,2) + xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal(1,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(xc(1,1) -
    xc(2,1)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) + xc(3,0)) + xc(1,0)*(xc(2,2)
    - xc(3,2)) + xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(2,2)=
    (-4*normsquare*(2*normal(1,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(xc(1,1) -
    xc(2,1)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) + 3*pow(2*normal(1,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(xc(1,1) - xc(2,1)),2)*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    4*normsquare*(pow(xc(1,0) - xc(2,0),2) + pow(xc(1,1) -
    xc(2,1),2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(2,3)=
    (4*normpowfour*(-xc(2,1) + xc(3,1)) -
    2*normsquare*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) + 2*normal(1,0)*(xc(0,2) -
    xc(2,2)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) + 3*(2*normal(1,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(xc(1,1) - xc(2,1)))*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) +
    2*normal(1,0)*(xc(0,2) - xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    4*normsquare*(2*xc(0,2)*(xc(1,0) - xc(2,0)) + xc(1,2)*xc(2,0) -
    2*xc(1,0)*xc(2,2) + xc(2,0)*xc(2,2) + xc(0,0)*(-xc(1,2) +
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal(1,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(xc(1,1) - xc(2,1)))*((xc(0,2) - xc(2,2))*(-xc(0,1) + xc(3,1)) +
    (xc(0,1) - xc(2,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(2,4)=
    (4*normpowfour*(xc(2,0) - xc(3,0)) -
    2*normsquare*(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(xc(0,2) -
    xc(2,2)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) + 3*(2*normal(1,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(xc(1,1) - xc(2,1)))*(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) -
    2*normal(0,0)*(xc(0,2) - xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    4*normsquare*(2*xc(0,2)*(xc(1,1) - xc(2,1)) + xc(1,2)*xc(2,1) -
    2*xc(1,1)*xc(2,2) + xc(2,1)*xc(2,2) + xc(0,1)*(-xc(1,2) +
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal(1,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(xc(1,1) - xc(2,1)))*((xc(0,2) - xc(2,2))*(xc(0,0) - xc(3,0)) +
    (-xc(0,0) + xc(2,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(2,5)=
    (-2*normsquare*(2*normal(1,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(xc(1,1) -
    xc(2,1)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) - 2*normsquare*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) -
    2*normal(0,0)*(-xc(0,1) + xc(2,1)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0)
    - xc(2,0)*xc(3,1) + xc(1,0)*(-xc(2,1) + xc(3,1))) + 3*(2*normal(1,0)*(xc(1,0) -
    xc(2,0)) - 2*normal(0,0)*(xc(1,1) - xc(2,1)))*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) -
    2*normal(0,0)*(-xc(0,1) + xc(2,1)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*(xc(1,0) - xc(2,0))*(-xc(0,0) + xc(2,0)) + 2*(xc(1,1) -
    xc(2,1))*(-xc(0,1) + xc(2,1)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) -
    xc(3,2))))/(4.*normpowfive);

  elematrix(2,6)=
    (4*normpowfour*(xc(1,1) - xc(3,1)) -
    2*normsquare*(-2*normal(2,0)*(-xc(0,1) + xc(1,1)) + 2*normal(1,0)*(-xc(0,2) +
    xc(1,2)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) + 3*(-2*normal(2,0)*(-xc(0,1) + xc(1,1)) +
    2*normal(1,0)*(-xc(0,2) + xc(1,2)))*(2*normal(1,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(xc(1,1) - xc(2,1)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    4*normsquare*(xc(1,0)*xc(1,2) - 2*xc(0,2)*(xc(1,0) - xc(2,0)) -
    2*xc(1,2)*xc(2,0) + xc(0,0)*(xc(1,2) - xc(2,2)) +
    xc(1,0)*xc(2,2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) +
    xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal(1,0)*(xc(1,0) -
    xc(2,0)) - 2*normal(0,0)*(xc(1,1) - xc(2,1)))*((xc(0,2) - xc(1,2))*(xc(0,1) -
    xc(3,1)) + (-xc(0,1) + xc(1,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(2,7)=
    (4*normpowfour*(-xc(1,0) + xc(3,0)) -
    2*normsquare*(-2*normal(2,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(-xc(0,2) +
    xc(1,2)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) + 3*(-2*normal(2,0)*(xc(0,0) - xc(1,0)) -
    2*normal(0,0)*(-xc(0,2) + xc(1,2)))*(2*normal(1,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(xc(1,1) - xc(2,1)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    4*normsquare*(xc(1,1)*xc(1,2) - 2*xc(0,2)*(xc(1,1) - xc(2,1)) -
    2*xc(1,2)*xc(2,1) + xc(0,1)*(xc(1,2) - xc(2,2)) +
    xc(1,1)*xc(2,2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) +
    xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal(1,0)*(xc(1,0) -
    xc(2,0)) - 2*normal(0,0)*(xc(1,1) - xc(2,1)))*((-xc(0,2) + xc(1,2))*(xc(0,0) -
    xc(3,0)) + (xc(0,0) - xc(1,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(2,8)=
    (-2*normsquare*(2*normal(1,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(xc(1,1) -
    xc(2,1)))*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(-xc(0,1) + xc(3,1))) - 2*normsquare*(2*normal(1,0)*(xc(0,0) - xc(1,0)) -
    2*normal(0,0)*(xc(0,1) - xc(1,1)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0)
    - xc(2,0)*xc(3,1) + xc(1,0)*(-xc(2,1) + xc(3,1))) + 3*(2*normal(1,0)*(xc(0,0) -
    xc(1,0)) - 2*normal(0,0)*(xc(0,1) - xc(1,1)))*(2*normal(1,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(xc(1,1) - xc(2,1)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*(xc(0,0) - xc(1,0))*(xc(1,0) - xc(2,0)) + 2*(xc(0,1) -
    xc(1,1))*(xc(1,1) - xc(2,1)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) -
    xc(3,2))))/(4.*normpowfive);

  elematrix(2,9)=
    (normal(2,0)*(xc(0,2)*xc(1,0)*xc(1,2) + pow(xc(1,1),2)*xc(2,0) -
    xc(0,2)*xc(1,2)*xc(2,0) + pow(xc(1,2),2)*xc(2,0) + xc(0,1)*(xc(1,0) -
    xc(2,0))*(xc(1,1) - xc(2,1)) - xc(1,0)*xc(1,1)*xc(2,1) - xc(1,1)*xc(2,0)*xc(2,1)
    + xc(1,0)*pow(xc(2,1),2) - xc(0,2)*xc(1,0)*xc(2,2) - xc(1,0)*xc(1,2)*xc(2,2) +
    xc(0,2)*xc(2,0)*xc(2,2) - xc(1,2)*xc(2,0)*xc(2,2) + xc(1,0)*pow(xc(2,2),2) -
    xc(0,0)*(pow(xc(1,1),2) + pow(xc(1,2),2) - 2*xc(1,1)*xc(2,1) +
    pow(xc(2,1),2) - 2*xc(1,2)*xc(2,2) + pow(xc(2,2),2))))/normcube;

  elematrix(2,10)=
    -((normal(2,0)*(-(xc(0,2)*xc(1,1)*xc(1,2)) + xc(1,0)*xc(1,1)*xc(2,0) -
    xc(1,1)*pow(xc(2,0),2) - xc(0,0)*(xc(1,0) - xc(2,0))*(xc(1,1) - xc(2,1)) -
    pow(xc(1,0),2)*xc(2,1) + xc(0,2)*xc(1,2)*xc(2,1) - pow(xc(1,2),2)*xc(2,1) +
    xc(1,0)*xc(2,0)*xc(2,1) + xc(0,2)*xc(1,1)*xc(2,2) + xc(1,1)*xc(1,2)*xc(2,2) -
    xc(0,2)*xc(2,1)*xc(2,2) + xc(1,2)*xc(2,1)*xc(2,2) - xc(1,1)*pow(xc(2,2),2) +
    xc(0,1)*(pow(xc(1,0),2) + pow(xc(1,2),2) - 2*xc(1,0)*xc(2,0) +
    pow(xc(2,0),2) - 2*xc(1,2)*xc(2,2) +pow(xc(2,2),2))))/normcube);

  elematrix(2,11)=
    -((2*normal(1,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(xc(1,1) -
    xc(2,1)))*(-(xc(1,1)*xc(2,0)) + xc(0,1)*(-xc(1,0) + xc(2,0)) + xc(0,0)*(xc(1,1)
    - xc(2,1)) + xc(1,0)*xc(2,1)))/(2.*normcube);

  elematrix(3,0)=
    (3*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) + 2*normal(1,0)*(xc(0,2) -
    xc(2,2)))*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) + 2*normal(1,0)*(-xc(1,2) +
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*(xc(0,1) - xc(2,1))*(-xc(1,1) +
    xc(2,1)) + 2*(xc(0,2) - xc(2,2))*(-xc(1,2) + xc(2,2)))*(-(normal(0,0)*(xc(0,0) -
    xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) + 2*normal(1,0)*(-xc(1,2) +
    xc(2,2)))*((xc(0,2) - xc(2,2))*(-xc(0,1) + xc(3,1)) + (xc(0,1) -
    xc(2,1))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) +
    2*normal(1,0)*(xc(0,2) - xc(2,2)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1)
    - xc(2,1)*xc(3,2) + xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(3,1)=
    (-4*normsquare*(2*xc(0,1)*(xc(1,0) - xc(2,0)) + xc(1,1)*xc(2,0) -
    2*xc(1,0)*xc(2,1) + xc(2,0)*xc(2,1) + xc(0,0)*(-xc(1,1) +
    xc(2,1)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) + 3*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) +
    2*normal(1,0)*(xc(0,2) - xc(2,2)))*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(-xc(1,2) + xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(-xc(1,2) +
    xc(2,2)))*((xc(0,2) - xc(2,2))*(-xc(0,1) + xc(3,1)) + (xc(0,1) -
    xc(2,1))*(xc(0,2) - xc(3,2))) + 4*normpowfour*(xc(2,2) - xc(3,2)) -
    2*normsquare*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) + 2*normal(1,0)*(xc(0,2) -
    xc(2,2)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) + xc(3,0)) + xc(1,0)*(xc(2,2)
    - xc(3,2)) + xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(3,2)=
    (4*normpowfour*(-xc(2,1) + xc(3,1)) -
    2*normsquare*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) + 2*normal(1,0)*(xc(0,2) -
    xc(2,2)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) + 3*(2*normal(1,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(xc(1,1) - xc(2,1)))*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) +
    2*normal(1,0)*(xc(0,2) - xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    4*normsquare*(2*xc(0,2)*(xc(1,0) - xc(2,0)) + xc(1,2)*xc(2,0) -
    2*xc(1,0)*xc(2,2) + xc(2,0)*xc(2,2) + xc(0,0)*(-xc(1,2) +
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal(1,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(xc(1,1) - xc(2,1)))*((xc(0,2) - xc(2,2))*(-xc(0,1) + xc(3,1)) +
    (xc(0,1) - xc(2,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(3,3)=
    (3*pow(-2*normal(2,0)*(xc(0,1) - xc(2,1)) + 2*normal(1,0)*(xc(0,2) -
    xc(2,2)),2)*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 4*normsquare*(pow(xc(0,1) - xc(2,1),2) +
    pow(xc(0,2) - xc(2,2),2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) + 2*normal(1,0)*(xc(0,2) -
    xc(2,2)))*((xc(0,2) - xc(2,2))*(-xc(0,1) + xc(3,1)) + (xc(0,1) -
    xc(2,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(3,4)=
    (-4*normsquare*(-xc(0,0) + xc(2,0))*(xc(0,1) - xc(2,1))*(-(normal(0,0)*(xc(0,0) -
    xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) +
    3*(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(xc(0,2) -
    xc(2,2)))*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) + 2*normal(1,0)*(xc(0,2) -
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(xc(0,1) - xc(2,1))
    + 2*normal(1,0)*(xc(0,2) - xc(2,2)))*((xc(0,2) - xc(2,2))*(xc(0,0) - xc(3,0)) +
    (-xc(0,0) + xc(2,0))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(-xc(0,0)
    + xc(2,0)) - 2*normal(0,0)*(xc(0,2) - xc(2,2)))*((xc(0,2) - xc(2,2))*(-xc(0,1) +
    xc(3,1)) + (xc(0,1) - xc(2,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(3,5)=
    (-2*normsquare*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) + 2*normal(1,0)*(xc(0,2) -
    xc(2,2)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) + 3*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) -
    2*normal(0,0)*(-xc(0,1) + xc(2,1)))*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) +
    2*normal(1,0)*(xc(0,2) - xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-xc(0,0) + xc(2,0))*(xc(0,2) - xc(2,2))*(-(normal(0,0)*(xc(0,0) -
    xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(-xc(0,1) +
    xc(2,1)))*((xc(0,2) - xc(2,2))*(-xc(0,1) + xc(3,1)) + (xc(0,1) -
    xc(2,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(3,6)=
    (3*(-2*normal(2,0)*(-xc(0,1) + xc(1,1)) + 2*normal(1,0)*(-xc(0,2) +
    xc(1,2)))*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) + 2*normal(1,0)*(xc(0,2) -
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*(-xc(0,1) + xc(1,1))*(xc(0,1) -
    xc(2,1)) + 2*(-xc(0,2) + xc(1,2))*(xc(0,2) - xc(2,2)))*(-(normal(0,0)*(xc(0,0) -
    xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) + 2*normal(1,0)*(xc(0,2) -
    xc(2,2)))*((xc(0,2) - xc(1,2))*(xc(0,1) - xc(3,1)) + (-xc(0,1) +
    xc(1,1))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(-xc(0,1) + xc(1,1))
    + 2*normal(1,0)*(-xc(0,2) + xc(1,2)))*((xc(0,2) - xc(2,2))*(-xc(0,1) + xc(3,1)) +
    (xc(0,1) - xc(2,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(3,7)=
    (-4*normsquare*(-(xc(1,1)*xc(2,0)) + xc(0,1)*(-2*xc(1,0) + xc(2,0)) +
    xc(0,0)*(xc(0,1) + xc(1,1) - 2*xc(2,1)) +
    2*xc(1,0)*xc(2,1))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) +
    xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) + 3*(-2*normal(2,0)*(xc(0,0) - xc(1,0))
    - 2*normal(0,0)*(-xc(0,2) + xc(1,2)))*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) +
    2*normal(1,0)*(xc(0,2) - xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) + 2*normal(1,0)*(xc(0,2) -
    xc(2,2)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(xc(0,0) - xc(1,0)) -
    2*normal(0,0)*(-xc(0,2) + xc(1,2)))*((xc(0,2) - xc(2,2))*(-xc(0,1) + xc(3,1)) +
    (xc(0,1) - xc(2,1))*(xc(0,2) - xc(3,2))) + 4*normpowfour*(-xc(0,2) +
    xc(3,2)))/(4.*normpowfive);

  elematrix(3,8)=
    (4*normpowfour*(xc(0,1) - xc(3,1)) - 2*normsquare*(-2*normal(2,0)*(xc(0,1)
    - xc(2,1)) + 2*normal(1,0)*(xc(0,2) - xc(2,2)))*((xc(0,1) - xc(1,1))*(xc(0,0) -
    xc(3,0)) + (xc(0,0) - xc(1,0))*(-xc(0,1) + xc(3,1))) + 3*(2*normal(1,0)*(xc(0,0) -
    xc(1,0)) - 2*normal(0,0)*(xc(0,1) - xc(1,1)))*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) +
    2*normal(1,0)*(xc(0,2) - xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-(xc(1,2)*xc(2,0)) + xc(0,2)*(-2*xc(1,0) + xc(2,0)) +
    xc(0,0)*(xc(0,2) + xc(1,2) - 2*xc(2,2)) +
    2*xc(1,0)*xc(2,2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) +
    xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal(1,0)*(xc(0,0) -
    xc(1,0)) - 2*normal(0,0)*(xc(0,1) - xc(1,1)))*((xc(0,2) - xc(2,2))*(-xc(0,1) +
    xc(3,1)) + (xc(0,1) - xc(2,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(3,9)=
    -((-2*normal(2,0)*(xc(0,1) - xc(2,1)) + 2*normal(1,0)*(xc(0,2) -
    xc(2,2)))*(-(xc(1,2)*xc(2,1)) + xc(0,2)*(-xc(1,1) + xc(2,1)) + xc(0,1)*(xc(1,2)
    - xc(2,2)) + xc(1,1)*xc(2,2)))/(2.*normcube);

  elematrix(3,10)=
    -((normal(0,0)*(xc(0,1)*xc(1,0)*xc(2,0) - xc(0,1)*pow(xc(2,0),2) +
    xc(1,1)*pow(xc(2,0),2) + pow(xc(0,0),2)*(xc(1,1) - xc(2,1)) +
    pow(xc(0,2),2)*(xc(1,1) - xc(2,1)) - xc(1,0)*xc(2,0)*xc(2,1) +
    xc(0,0)*(-2*xc(1,1)*xc(2,0) + xc(0,1)*(-xc(1,0) + xc(2,0)) + (xc(1,0) +
    xc(2,0))*xc(2,1)) + xc(0,1)*xc(1,2)*xc(2,2) - xc(1,2)*xc(2,1)*xc(2,2) -
    xc(0,1)*pow(xc(2,2),2) + xc(1,1)*pow(xc(2,2),2) + xc(0,2)*(xc(1,2)*xc(2,1) +
    (-2*xc(1,1) + xc(2,1))*xc(2,2) + xc(0,1)*(-xc(1,2) +
    xc(2,2)))))/normcube);

  elematrix(3,11)=
    -((normal(0,0)*(xc(0,2)*xc(1,0)*xc(2,0) - xc(0,2)*pow(xc(2,0),2) +
    xc(1,2)*pow(xc(2,0),2) + xc(0,2)*xc(1,1)*xc(2,1) - xc(0,2)*pow(xc(2,1),2) +
    xc(1,2)*pow(xc(2,1),2) + pow(xc(0,0),2)*(xc(1,2) - xc(2,2)) +
    pow(xc(0,1),2)*(xc(1,2) - xc(2,2)) - xc(1,0)*xc(2,0)*xc(2,2) -
    xc(1,1)*xc(2,1)*xc(2,2) + xc(0,0)*(-2*xc(1,2)*xc(2,0) + xc(0,2)*(-xc(1,0) +
    xc(2,0)) + (xc(1,0) + xc(2,0))*xc(2,2)) + xc(0,1)*(-2*xc(1,2)*xc(2,1) +
    xc(0,2)*(-xc(1,1) + xc(2,1)) + (xc(1,1) +
    xc(2,1))*xc(2,2))))/normcube);

  elematrix(4,0)=
    (-4*normsquare*(-2*xc(1,1)*xc(2,0) + xc(0,1)*(-xc(1,0) + xc(2,0)) +
    2*xc(0,0)*(xc(1,1) - xc(2,1)) + xc(1,0)*xc(2,1) +
    xc(2,0)*xc(2,1))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) +
    xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) + 3*(-2*normal(2,0)*(-xc(0,0) + xc(2,0))
    - 2*normal(0,0)*(xc(0,2) - xc(2,2)))*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) +
    2*normal(1,0)*(-xc(1,2) + xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) + 2*normal(1,0)*(-xc(1,2) +
    xc(2,2)))*((xc(0,2) - xc(2,2))*(xc(0,0) - xc(3,0)) + (-xc(0,0) +
    xc(2,0))*(xc(0,2) - xc(3,2))) + 4*normpowfour*(-xc(2,2) + xc(3,2)) -
    2*normsquare*(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(xc(0,2) -
    xc(2,2)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) +
    xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(4,1)=
    (3*(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(xc(0,2) -
    xc(2,2)))*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(-xc(1,2) +
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*(xc(1,0) - xc(2,0))*(-xc(0,0) +
    xc(2,0)) + 2*(xc(0,2) - xc(2,2))*(-xc(1,2) + xc(2,2)))*(-(normal(0,0)*(xc(0,0) -
    xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(-xc(1,2) +
    xc(2,2)))*((xc(0,2) - xc(2,2))*(xc(0,0) - xc(3,0)) + (-xc(0,0) +
    xc(2,0))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(-xc(0,0) + xc(2,0))
    - 2*normal(0,0)*(xc(0,2) - xc(2,2)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) +
    xc(3,0)) + xc(1,0)*(xc(2,2) - xc(3,2)) +xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(4,2)=
    (4*normpowfour*(xc(2,0) - xc(3,0)) -
    2*normsquare*(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(xc(0,2) -
    xc(2,2)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) + 3*(2*normal(1,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(xc(1,1) - xc(2,1)))*(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) -
    2*normal(0,0)*(xc(0,2) - xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    4*normsquare*(2*xc(0,2)*(xc(1,1) - xc(2,1)) + xc(1,2)*xc(2,1) -
    2*xc(1,1)*xc(2,2) + xc(2,1)*xc(2,2) + xc(0,1)*(-xc(1,2) +
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal(1,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(xc(1,1) - xc(2,1)))*((xc(0,2) - xc(2,2))*(xc(0,0) - xc(3,0)) +
    (-xc(0,0) + xc(2,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(4,3)=
    (-4*normsquare*(-xc(0,0) + xc(2,0))*(xc(0,1) - xc(2,1))*(-(normal(0,0)*(xc(0,0) -
    xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) +
    3*(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(xc(0,2) -
    xc(2,2)))*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) + 2*normal(1,0)*(xc(0,2) -
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(xc(0,1) - xc(2,1))
    + 2*normal(1,0)*(xc(0,2) - xc(2,2)))*((xc(0,2) - xc(2,2))*(xc(0,0) - xc(3,0)) +
    (-xc(0,0) + xc(2,0))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(-xc(0,0)
    + xc(2,0)) - 2*normal(0,0)*(xc(0,2) - xc(2,2)))*((xc(0,2) - xc(2,2))*(-xc(0,1) +
    xc(3,1)) + (xc(0,1) - xc(2,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(4,4)=
    (3*pow(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(xc(0,2) -
    xc(2,2)),2)*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 4*normsquare*(pow(xc(0,0) - xc(2,0),2) +
    pow(xc(0,2) - xc(2,2),2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(xc(0,2) -
    xc(2,2)))*((xc(0,2) - xc(2,2))*(xc(0,0) - xc(3,0)) + (-xc(0,0) +
    xc(2,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(4,5)=
    (-2*normsquare*(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(xc(0,2) -
    xc(2,2)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) + 3*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) -
    2*normal(0,0)*(-xc(0,1) + xc(2,1)))*(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) -
    2*normal(0,0)*(xc(0,2) - xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-xc(0,1) + xc(2,1))*(xc(0,2) - xc(2,2))*(-(normal(0,0)*(xc(0,0) -
    xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(-xc(0,1) +
    xc(2,1)))*((xc(0,2) - xc(2,2))*(xc(0,0) - xc(3,0)) + (-xc(0,0) +
    xc(2,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(4,6)=
    (-4*normsquare*(xc(0,1)*(xc(1,0) - 2*xc(2,0)) + 2*xc(1,1)*xc(2,0) -
    xc(1,0)*xc(2,1) + xc(0,0)*(xc(0,1) - 2*xc(1,1) + xc(2,1)))*(-(normal(0,0)*(xc(0,0)
    - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) +
    3*(-2*normal(2,0)*(-xc(0,1) + xc(1,1)) + 2*normal(1,0)*(-xc(0,2) +
    xc(1,2)))*(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(xc(0,2) -
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(-xc(0,0) + xc(2,0))
    - 2*normal(0,0)*(xc(0,2) - xc(2,2)))*((xc(0,2) - xc(1,2))*(xc(0,1) - xc(3,1)) +
    (-xc(0,1) + xc(1,1))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(-xc(0,1)
    + xc(1,1)) + 2*normal(1,0)*(-xc(0,2) + xc(1,2)))*((xc(0,2) - xc(2,2))*(xc(0,0) -
    xc(3,0)) + (-xc(0,0) + xc(2,0))*(xc(0,2) - xc(3,2))) +
    4*normpowfour*(xc(0,2) - xc(3,2)))/(4.*normpowfive);

  elematrix(4,7)=
    (3*(-2*normal(2,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(-xc(0,2) +
    xc(1,2)))*(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(xc(0,2) -
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*(xc(0,0) - xc(1,0))*(-xc(0,0) +
    xc(2,0)) + 2*(-xc(0,2) + xc(1,2))*(xc(0,2) - xc(2,2)))*(-(normal(0,0)*(xc(0,0) -
    xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(xc(0,2) -
    xc(2,2)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(xc(0,0) - xc(1,0)) -
    2*normal(0,0)*(-xc(0,2) + xc(1,2)))*((xc(0,2) - xc(2,2))*(xc(0,0) - xc(3,0)) +
    (-xc(0,0) + xc(2,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(4,8)=
    (4*normpowfour*(-xc(0,0) + xc(3,0)) -
    2*normsquare*(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(xc(0,2) -
    xc(2,2)))*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(-xc(0,1) + xc(3,1))) + 3*(2*normal(1,0)*(xc(0,0) - xc(1,0)) -
    2*normal(0,0)*(xc(0,1) - xc(1,1)))*(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) -
    2*normal(0,0)*(xc(0,2) - xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-(xc(1,2)*xc(2,1)) + xc(0,2)*(-2*xc(1,1) + xc(2,1)) +
    xc(0,1)*(xc(0,2) + xc(1,2) - 2*xc(2,2)) +
    2*xc(1,1)*xc(2,2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) +
    xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal(1,0)*(xc(0,0) -
    xc(1,0)) - 2*normal(0,0)*(xc(0,1) - xc(1,1)))*((xc(0,2) - xc(2,2))*(xc(0,0) -
    xc(3,0)) + (-xc(0,0) + xc(2,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(4,9)=
    -((normal(1,0)*(pow(xc(0,1),2)*(xc(1,0) - xc(2,0)) + pow(xc(0,2),2)*(xc(1,0) -
    xc(2,0)) + xc(0,0)*xc(1,1)*xc(2,1) - xc(1,1)*xc(2,0)*xc(2,1) -
    xc(0,0)*pow(xc(2,1),2) + xc(1,0)*pow(xc(2,1),2) + xc(0,1)*(xc(1,1)*xc(2,0) +
    (-2*xc(1,0) + xc(2,0))*xc(2,1) + xc(0,0)*(-xc(1,1) + xc(2,1))) +
    xc(0,0)*xc(1,2)*xc(2,2) - xc(1,2)*xc(2,0)*xc(2,2) - xc(0,0)*pow(xc(2,2),2) +
    xc(1,0)*pow(xc(2,2),2) + xc(0,2)*(xc(1,2)*xc(2,0) + (-2*xc(1,0) +
    xc(2,0))*xc(2,2) + xc(0,0)*(-xc(1,2) + xc(2,2)))))/normcube);

  elematrix(4,10)=
    -(normal(1,0)*(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(xc(0,2) -
    xc(2,2))))/(2.*normcube);

  elematrix(4,11)=
    -((normal(1,0)*(xc(0,2)*xc(1,0)*xc(2,0) - xc(0,2)*pow(xc(2,0),2) +
    xc(1,2)*pow(xc(2,0),2) + xc(0,2)*xc(1,1)*xc(2,1) - xc(0,2)*pow(xc(2,1),2) +
    xc(1,2)*pow(xc(2,1),2) + pow(xc(0,0),2)*(xc(1,2) - xc(2,2)) +
    pow(xc(0,1),2)*(xc(1,2) - xc(2,2)) - xc(1,0)*xc(2,0)*xc(2,2) -
    xc(1,1)*xc(2,1)*xc(2,2) + xc(0,0)*(-2*xc(1,2)*xc(2,0) + xc(0,2)*(-xc(1,0) +
    xc(2,0)) + (xc(1,0) + xc(2,0))*xc(2,2)) + xc(0,1)*(-2*xc(1,2)*xc(2,1) +
    xc(0,2)*(-xc(1,1) + xc(2,1)) + (xc(1,1) +xc(2,1))*xc(2,2))))/normcube);

  elematrix(5,0)=
    (-2*normsquare*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) + 2*normal(1,0)*(-xc(1,2) +
    xc(2,2)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) + 4*normpowfour*(xc(2,1) - xc(3,1)) -
    4*normsquare*(-2*xc(1,2)*xc(2,0) + xc(0,2)*(-xc(1,0) + xc(2,0)) +
    2*xc(0,0)*(xc(1,2) - xc(2,2)) + xc(1,0)*xc(2,2) +
    xc(2,0)*xc(2,2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) +
    xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) + 3*(2*normal(1,0)*(-xc(0,0) + xc(2,0))
    - 2*normal(0,0)*(-xc(0,1) + xc(2,1)))*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) +
    2*normal(1,0)*(-xc(1,2) + xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(-xc(0,1) +
    xc(2,1)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) +
    xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(5,1)=
    (4*normpowfour*(-xc(2,0) + xc(3,0)) -
    2*normsquare*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(-xc(1,2) +
    xc(2,2)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) - 4*normsquare*(-2*xc(1,2)*xc(2,1) +
    xc(0,2)*(-xc(1,1) + xc(2,1)) + 2*xc(0,1)*(xc(1,2) - xc(2,2)) + xc(1,1)*xc(2,2) +
    xc(2,1)*xc(2,2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) +
    xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) + 3*(2*normal(1,0)*(-xc(0,0) + xc(2,0))
    - 2*normal(0,0)*(-xc(0,1) + xc(2,1)))*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(-xc(1,2) + xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(-xc(0,1) +
    xc(2,1)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) + xc(3,0)) + xc(1,0)*(xc(2,2)
    - xc(3,2)) + xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(5,2)=
    (-2*normsquare*(2*normal(1,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(xc(1,1) -
    xc(2,1)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) - 2*normsquare*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) -
    2*normal(0,0)*(-xc(0,1) + xc(2,1)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0)
    - xc(2,0)*xc(3,1) + xc(1,0)*(-xc(2,1) + xc(3,1))) + 3*(2*normal(1,0)*(xc(1,0) -
    xc(2,0)) - 2*normal(0,0)*(xc(1,1) - xc(2,1)))*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) -
    2*normal(0,0)*(-xc(0,1) + xc(2,1)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*(xc(1,0) - xc(2,0))*(-xc(0,0) + xc(2,0)) + 2*(xc(1,1) -
    xc(2,1))*(-xc(0,1) + xc(2,1)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) -
    xc(3,2))))/(4.*normpowfive);

  elematrix(5,3)=
    (-2*normsquare*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) + 2*normal(1,0)*(xc(0,2) -
    xc(2,2)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) + 3*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) -
    2*normal(0,0)*(-xc(0,1) + xc(2,1)))*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) +
    2*normal(1,0)*(xc(0,2) - xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-xc(0,0) + xc(2,0))*(xc(0,2) - xc(2,2))*(-(normal(0,0)*(xc(0,0) -
    xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(-xc(0,1) +
    xc(2,1)))*((xc(0,2) - xc(2,2))*(-xc(0,1) + xc(3,1)) + (xc(0,1) -
    xc(2,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(5,4)=
    (-2*normsquare*(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(xc(0,2) -
    xc(2,2)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) + 3*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) -
    2*normal(0,0)*(-xc(0,1) + xc(2,1)))*(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) -
    2*normal(0,0)*(xc(0,2) - xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-xc(0,1) + xc(2,1))*(xc(0,2) - xc(2,2))*(-(normal(0,0)*(xc(0,0) -
    xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(-xc(0,1) +
    xc(2,1)))*((xc(0,2) - xc(2,2))*(xc(0,0) - xc(3,0)) + (-xc(0,0) +
    xc(2,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(5,5)=
    (-4*normsquare*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(-xc(0,1) +
    xc(2,1)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) - 4*normsquare*(pow(xc(0,0) - xc(2,0),2) +
    pow(xc(0,1) - xc(2,1),2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) +
    3*pow(2*normal(1,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(-xc(0,1) +
    xc(2,1)),2)*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(5,6)=
    (-2*normsquare*(-2*normal(2,0)*(-xc(0,1) + xc(1,1)) + 2*normal(1,0)*(-xc(0,2) +
    xc(1,2)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) + 4*normpowfour*(-xc(0,1) + xc(3,1)) +
    3*(-2*normal(2,0)*(-xc(0,1) + xc(1,1)) + 2*normal(1,0)*(-xc(0,2) +
    xc(1,2)))*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(-xc(0,1) +
    xc(2,1)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 4*normsquare*(xc(0,2)*(xc(1,0) - 2*xc(2,0)) +
    2*xc(1,2)*xc(2,0) - xc(1,0)*xc(2,2) + xc(0,0)*(xc(0,2) - 2*xc(1,2) +
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal(1,0)*(-xc(0,0) + xc(2,0))
    - 2*normal(0,0)*(-xc(0,1) + xc(2,1)))*((xc(0,2) - xc(1,2))*(xc(0,1) - xc(3,1)) +
    (-xc(0,1) + xc(1,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(5,7)=
    (4*normpowfour*(xc(0,0) - xc(3,0)) - 2*normsquare*(-2*normal(2,0)*(xc(0,0)
    - xc(1,0)) - 2*normal(0,0)*(-xc(0,2) + xc(1,2)))*((-xc(0,1) + xc(2,1))*(xc(0,0) -
    xc(3,0)) + (xc(0,0) - xc(2,0))*(xc(0,1) - xc(3,1))) + 3*(-2*normal(2,0)*(xc(0,0) -
    xc(1,0)) - 2*normal(0,0)*(-xc(0,2) + xc(1,2)))*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) -
    2*normal(0,0)*(-xc(0,1) + xc(2,1)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    4*normsquare*(xc(0,2)*(xc(1,1) - 2*xc(2,1)) + 2*xc(1,2)*xc(2,1) -
    xc(1,1)*xc(2,2) + xc(0,1)*(xc(0,2) - 2*xc(1,2) + xc(2,2)))*(-(normal(0,0)*(xc(0,0)
    - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(-xc(0,1) +
    xc(2,1)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(5,8)=
    (-2*normsquare*(2*normal(1,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(xc(0,1) -
    xc(1,1)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) - 2*normsquare*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) -
    2*normal(0,0)*(-xc(0,1) + xc(2,1)))*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) +
    (xc(0,0) - xc(1,0))*(-xc(0,1) + xc(3,1))) + 3*(2*normal(1,0)*(xc(0,0) - xc(1,0)) -
    2*normal(0,0)*(xc(0,1) - xc(1,1)))*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) -
    2*normal(0,0)*(-xc(0,1) + xc(2,1)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*(xc(0,0) - xc(1,0))*(-xc(0,0) + xc(2,0)) + 2*(xc(0,1) -
    xc(1,1))*(-xc(0,1) + xc(2,1)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) -
    xc(3,2))))/(4.*normpowfive);

  elematrix(5,9)=
    -((normal(2,0)*(pow(xc(0,1),2)*(xc(1,0) - xc(2,0)) + pow(xc(0,2),2)*(xc(1,0) -
    xc(2,0)) + xc(0,0)*xc(1,1)*xc(2,1) - xc(1,1)*xc(2,0)*xc(2,1) -
    xc(0,0)*pow(xc(2,1),2) + xc(1,0)*pow(xc(2,1),2) + xc(0,1)*(xc(1,1)*xc(2,0) +
    (-2*xc(1,0) + xc(2,0))*xc(2,1) + xc(0,0)*(-xc(1,1) + xc(2,1))) +
    xc(0,0)*xc(1,2)*xc(2,2) - xc(1,2)*xc(2,0)*xc(2,2) - xc(0,0)*pow(xc(2,2),2) +
    xc(1,0)*pow(xc(2,2),2) + xc(0,2)*(xc(1,2)*xc(2,0) + (-2*xc(1,0) +
    xc(2,0))*xc(2,2) + xc(0,0)*(-xc(1,2) + xc(2,2)))))/normcube);

  elematrix(5,10)=
    -((normal(2,0)*(xc(0,1)*xc(1,0)*xc(2,0) - xc(0,1)*pow(xc(2,0),2) +
    xc(1,1)*pow(xc(2,0),2) + pow(xc(0,0),2)*(xc(1,1) - xc(2,1)) +
    pow(xc(0,2),2)*(xc(1,1) - xc(2,1)) - xc(1,0)*xc(2,0)*xc(2,1) +
    xc(0,0)*(-2*xc(1,1)*xc(2,0) + xc(0,1)*(-xc(1,0) + xc(2,0)) + (xc(1,0) +
    xc(2,0))*xc(2,1)) + xc(0,1)*xc(1,2)*xc(2,2) - xc(1,2)*xc(2,1)*xc(2,2) -
    xc(0,1)*pow(xc(2,2),2) + xc(1,1)*pow(xc(2,2),2) + xc(0,2)*(xc(1,2)*xc(2,1) +
    (-2*xc(1,1) + xc(2,1))*xc(2,2) + xc(0,1)*(-xc(1,2) +
    xc(2,2)))))/normcube);

  elematrix(5,11)=
    -((-(xc(1,1)*xc(2,0)) + xc(0,1)*(-xc(1,0) + xc(2,0)) + xc(0,0)*(xc(1,1) -
    xc(2,1)) + xc(1,0)*xc(2,1))*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) -
    2*normal(0,0)*(-xc(0,1) + xc(2,1))))/(2.*normcube);

  elematrix(6,0)=
    (-2*normsquare*(2*(xc(0,1) - xc(1,1))*(xc(1,1) - xc(2,1)) + 2*(xc(0,2) -
    xc(1,2))*(xc(1,2) - xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) +
    3*(-2*normal(2,0)*(-xc(0,1) + xc(1,1)) + 2*normal(1,0)*(-xc(0,2) +
    xc(1,2)))*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) + 2*normal(1,0)*(-xc(1,2) +
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(-xc(1,1) + xc(2,1))
    + 2*normal(1,0)*(-xc(1,2) + xc(2,2)))*((xc(0,2) - xc(1,2))*(xc(0,1) - xc(3,1)) +
    (-xc(0,1) + xc(1,1))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(-xc(0,1)
    + xc(1,1)) + 2*normal(1,0)*(-xc(0,2) + xc(1,2)))*(xc(1,2)*(xc(2,1) - xc(3,1)) +
    xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) + xc(1,1)*(-xc(2,2) +
    xc(3,2))))/(4.*normpowfive);

  elematrix(6,1)=
    (-4*normsquare*(xc(1,0)*xc(1,1) - 2*xc(0,1)*(xc(1,0) - xc(2,0)) -
    2*xc(1,1)*xc(2,0) + xc(0,0)*(xc(1,1) - xc(2,1)) +
    xc(1,0)*xc(2,1))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) +
    xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) + 3*(-2*normal(2,0)*(-xc(0,1) + xc(1,1))
    + 2*normal(1,0)*(-xc(0,2) + xc(1,2)))*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(-xc(1,2) + xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(-xc(1,2) +
    xc(2,2)))*((xc(0,2) - xc(1,2))*(xc(0,1) - xc(3,1)) + (-xc(0,1) +
    xc(1,1))*(xc(0,2) - xc(3,2))) + 4*normpowfour*(-xc(1,2) + xc(3,2)) -
    2*normsquare*(-2*normal(2,0)*(-xc(0,1) + xc(1,1)) + 2*normal(1,0)*(-xc(0,2) +
    xc(1,2)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) + xc(3,0)) + xc(1,0)*(xc(2,2)
    - xc(3,2)) + xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(6,2)=
    (4*normpowfour*(xc(1,1) - xc(3,1)) -
    2*normsquare*(-2*normal(2,0)*(-xc(0,1) + xc(1,1)) + 2*normal(1,0)*(-xc(0,2) +
    xc(1,2)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) + 3*(-2*normal(2,0)*(-xc(0,1) + xc(1,1)) +
    2*normal(1,0)*(-xc(0,2) + xc(1,2)))*(2*normal(1,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(xc(1,1) - xc(2,1)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    4*normsquare*(xc(1,0)*xc(1,2) - 2*xc(0,2)*(xc(1,0) - xc(2,0)) -
    2*xc(1,2)*xc(2,0) + xc(0,0)*(xc(1,2) - xc(2,2)) +
    xc(1,0)*xc(2,2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) +
    xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal(1,0)*(xc(1,0) -
    xc(2,0)) - 2*normal(0,0)*(xc(1,1) - xc(2,1)))*((xc(0,2) - xc(1,2))*(xc(0,1) -
    xc(3,1)) + (-xc(0,1) + xc(1,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(6,3)=
    (3*(-2*normal(2,0)*(-xc(0,1) + xc(1,1)) + 2*normal(1,0)*(-xc(0,2) +
    xc(1,2)))*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) + 2*normal(1,0)*(xc(0,2) -
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*(-xc(0,1) + xc(1,1))*(xc(0,1) -
    xc(2,1)) + 2*(-xc(0,2) + xc(1,2))*(xc(0,2) - xc(2,2)))*(-(normal(0,0)*(xc(0,0) -
    xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) + 2*normal(1,0)*(xc(0,2) -
    xc(2,2)))*((xc(0,2) - xc(1,2))*(xc(0,1) - xc(3,1)) + (-xc(0,1) +
    xc(1,1))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(-xc(0,1) + xc(1,1))
    + 2*normal(1,0)*(-xc(0,2) + xc(1,2)))*((xc(0,2) - xc(2,2))*(-xc(0,1) + xc(3,1)) +
    (xc(0,1) - xc(2,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(6,4)=
    (-4*normsquare*(xc(0,1)*(xc(1,0) - 2*xc(2,0)) + 2*xc(1,1)*xc(2,0) -
    xc(1,0)*xc(2,1) + xc(0,0)*(xc(0,1) - 2*xc(1,1) + xc(2,1)))*(-(normal(0,0)*(xc(0,0)
    - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) +
    3*(-2*normal(2,0)*(-xc(0,1) + xc(1,1)) + 2*normal(1,0)*(-xc(0,2) +
    xc(1,2)))*(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(xc(0,2) -
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(-xc(0,0) + xc(2,0))
    - 2*normal(0,0)*(xc(0,2) - xc(2,2)))*((xc(0,2) - xc(1,2))*(xc(0,1) - xc(3,1)) +
    (-xc(0,1) + xc(1,1))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(-xc(0,1)
    + xc(1,1)) + 2*normal(1,0)*(-xc(0,2) + xc(1,2)))*((xc(0,2) - xc(2,2))*(xc(0,0) -
    xc(3,0)) + (-xc(0,0) + xc(2,0))*(xc(0,2) - xc(3,2))) +
    4*normpowfour*(xc(0,2) - xc(3,2)))/(4.*normpowfive);

  elematrix(6,5)=
    (-2*normsquare*(-2*normal(2,0)*(-xc(0,1) + xc(1,1)) + 2*normal(1,0)*(-xc(0,2) +
    xc(1,2)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) + 4*normpowfour*(-xc(0,1) + xc(3,1)) +
    3*(-2*normal(2,0)*(-xc(0,1) + xc(1,1)) + 2*normal(1,0)*(-xc(0,2) +
    xc(1,2)))*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(-xc(0,1) +
    xc(2,1)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 4*normsquare*(xc(0,2)*(xc(1,0) - 2*xc(2,0)) +
    2*xc(1,2)*xc(2,0) - xc(1,0)*xc(2,2) + xc(0,0)*(xc(0,2) - 2*xc(1,2) +
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal(1,0)*(-xc(0,0) + xc(2,0))
    - 2*normal(0,0)*(-xc(0,1) + xc(2,1)))*((xc(0,2) - xc(1,2))*(xc(0,1) - xc(3,1)) +
    (-xc(0,1) + xc(1,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(6,6)=
    (-4*normsquare*(pow(xc(0,1) - xc(1,1),2) + pow(xc(0,2) -
    xc(1,2),2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) + 3*pow(-2*normal(2,0)*(-xc(0,1) + xc(1,1)) +
    2*normal(1,0)*(-xc(0,2) + xc(1,2)),2)*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-2*normal(2,0)*(-xc(0,1) + xc(1,1)) + 2*normal(1,0)*(-xc(0,2) +
    xc(1,2)))*((xc(0,2) - xc(1,2))*(xc(0,1) - xc(3,1)) + (-xc(0,1) +
    xc(1,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(6,7)=
    (-4*normsquare*(xc(0,0) - xc(1,0))*(-xc(0,1) + xc(1,1))*(-(normal(0,0)*(xc(0,0) -
    xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) +
    3*(-2*normal(2,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(-xc(0,2) +
    xc(1,2)))*(-2*normal(2,0)*(-xc(0,1) + xc(1,1)) + 2*normal(1,0)*(-xc(0,2) +
    xc(1,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(-xc(0,1) + xc(1,1))
    + 2*normal(1,0)*(-xc(0,2) + xc(1,2)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) +
    (xc(0,0) - xc(1,0))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(xc(0,0) -
    xc(1,0)) - 2*normal(0,0)*(-xc(0,2) + xc(1,2)))*((xc(0,2) - xc(1,2))*(xc(0,1) -
    xc(3,1)) + (-xc(0,1) + xc(1,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(6,8)=
    (-2*normsquare*(-2*normal(2,0)*(-xc(0,1) + xc(1,1)) + 2*normal(1,0)*(-xc(0,2) +
    xc(1,2)))*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(-xc(0,1) + xc(3,1))) - 4*normsquare*(xc(0,0) - xc(1,0))*(-xc(0,2) +
    xc(1,2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) + 3*(2*normal(1,0)*(xc(0,0) - xc(1,0)) -
    2*normal(0,0)*(xc(0,1) - xc(1,1)))*(-2*normal(2,0)*(-xc(0,1) + xc(1,1)) +
    2*normal(1,0)*(-xc(0,2) + xc(1,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal(1,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(xc(0,1) -
    xc(1,1)))*((xc(0,2) - xc(1,2))*(xc(0,1) - xc(3,1)) + (-xc(0,1) +
    xc(1,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(6,9)=
    -((-2*normal(2,0)*(-xc(0,1) + xc(1,1)) + 2*normal(1,0)*(-xc(0,2) +
    xc(1,2)))*(-(xc(1,2)*xc(2,1)) + xc(0,2)*(-xc(1,1) + xc(2,1)) + xc(0,1)*(xc(1,2)
    - xc(2,2)) + xc(1,1)*xc(2,2)))/(2.*normcube);

  elematrix(6,10)=
    (normal(0,0)*(pow(xc(0,2),2)*xc(1,1) - xc(0,2)*xc(1,1)*xc(1,2) +
    xc(1,0)*xc(1,1)*xc(2,0) - xc(0,0)*(xc(0,1)*(xc(1,0) - xc(2,0)) + xc(1,1)*xc(2,0)
    + xc(1,0)*(xc(1,1) - 2*xc(2,1))) + pow(xc(0,0),2)*(xc(1,1) - xc(2,1)) -
    pow(xc(0,2),2)*xc(2,1) - pow(xc(1,0),2)*xc(2,1) + 2*xc(0,2)*xc(1,2)*xc(2,1)
    - pow(xc(1,2),2)*xc(2,1) + xc(0,1)*(pow(xc(1,0),2) - xc(1,0)*xc(2,0) -
    (xc(0,2) - xc(1,2))*(xc(1,2) - xc(2,2))) - xc(0,2)*xc(1,1)*xc(2,2) +
    xc(1,1)*xc(1,2)*xc(2,2)))/normcube;

  elematrix(6,11)=
    ((-(xc(1,2)*xc(2,1)) + xc(0,2)*(-xc(1,1) + xc(2,1)) + xc(0,1)*(xc(1,2) -
    xc(2,2)) + xc(1,1)*xc(2,2))*(pow(xc(0,1),2)*xc(1,2) - xc(0,1)*xc(1,1)*xc(1,2)
    + xc(1,0)*xc(1,2)*xc(2,0) + xc(0,2)*(pow(xc(1,0),2) - xc(1,0)*xc(2,0) -
    (xc(0,1) - xc(1,1))*(xc(1,1) - xc(2,1))) - xc(0,1)*xc(1,2)*xc(2,1) +
    xc(1,1)*xc(1,2)*xc(2,1) - xc(0,0)*(xc(0,2)*(xc(1,0) - xc(2,0)) + xc(1,2)*xc(2,0)
    + xc(1,0)*(xc(1,2) - 2*xc(2,2))) + pow(xc(0,0),2)*(xc(1,2) - xc(2,2)) -
    pow(xc(0,1),2)*xc(2,2) - pow(xc(1,0),2)*xc(2,2) + 2*xc(0,1)*xc(1,1)*xc(2,2)
    - pow(xc(1,1),2)*xc(2,2)))/normcube;

  elematrix(7,0)=
    (-4*normsquare*(xc(1,0)*xc(1,1) + xc(0,1)*(xc(1,0) - xc(2,0)) + xc(1,1)*xc(2,0)
    - 2*xc(0,0)*(xc(1,1) - xc(2,1)) - 2*xc(1,0)*xc(2,1))*(-(normal(0,0)*(xc(0,0) -
    xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) +
    3*(-2*normal(2,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(-xc(0,2) +
    xc(1,2)))*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) + 2*normal(1,0)*(-xc(1,2) +
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(-xc(1,1) + xc(2,1))
    + 2*normal(1,0)*(-xc(1,2) + xc(2,2)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) +
    (xc(0,0) - xc(1,0))*(xc(0,2) - xc(3,2))) + 4*normpowfour*(xc(1,2) -
    xc(3,2)) - 2*normsquare*(-2*normal(2,0)*(xc(0,0) - xc(1,0)) -
    2*normal(0,0)*(-xc(0,2) + xc(1,2)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1)
    - xc(2,1)*xc(3,2) + xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(7,1)=
    (-2*normsquare*(2*(xc(0,0) - xc(1,0))*(xc(1,0) - xc(2,0)) + 2*(xc(0,2) -
    xc(1,2))*(xc(1,2) - xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) +
    3*(-2*normal(2,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(-xc(0,2) +
    xc(1,2)))*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(-xc(1,2) +
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(xc(1,0) - xc(2,0))
    - 2*normal(0,0)*(-xc(1,2) + xc(2,2)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) +
    (xc(0,0) - xc(1,0))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(xc(0,0) -
    xc(1,0)) - 2*normal(0,0)*(-xc(0,2) + xc(1,2)))*(-(xc(2,2)*xc(3,0)) +
    xc(1,2)*(-xc(2,0) + xc(3,0)) + xc(1,0)*(xc(2,2) - xc(3,2)) +
    xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(7,2)=
    (4*normpowfour*(-xc(1,0) + xc(3,0)) -
    2*normsquare*(-2*normal(2,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(-xc(0,2) +
    xc(1,2)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) + 3*(-2*normal(2,0)*(xc(0,0) - xc(1,0)) -
    2*normal(0,0)*(-xc(0,2) + xc(1,2)))*(2*normal(1,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(xc(1,1) - xc(2,1)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    4*normsquare*(xc(1,1)*xc(1,2) - 2*xc(0,2)*(xc(1,1) - xc(2,1)) -
    2*xc(1,2)*xc(2,1) + xc(0,1)*(xc(1,2) - xc(2,2)) +
    xc(1,1)*xc(2,2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) +
    xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal(1,0)*(xc(1,0) -
    xc(2,0)) - 2*normal(0,0)*(xc(1,1) - xc(2,1)))*((-xc(0,2) + xc(1,2))*(xc(0,0) -
    xc(3,0)) + (xc(0,0) - xc(1,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(7,3)=
    (-4*normsquare*(-(xc(1,1)*xc(2,0)) + xc(0,1)*(-2*xc(1,0) + xc(2,0)) +
    xc(0,0)*(xc(0,1) + xc(1,1) - 2*xc(2,1)) +
    2*xc(1,0)*xc(2,1))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) +
    xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) + 3*(-2*normal(2,0)*(xc(0,0) - xc(1,0))
    - 2*normal(0,0)*(-xc(0,2) + xc(1,2)))*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) +
    2*normal(1,0)*(xc(0,2) - xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) + 2*normal(1,0)*(xc(0,2) -
    xc(2,2)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(xc(0,0) - xc(1,0)) -
    2*normal(0,0)*(-xc(0,2) + xc(1,2)))*((xc(0,2) - xc(2,2))*(-xc(0,1) + xc(3,1)) +
    (xc(0,1) - xc(2,1))*(xc(0,2) - xc(3,2))) + 4*normpowfour*(-xc(0,2) +
    xc(3,2)))/(4.*normpowfive);

  elematrix(7,4)=
    (3*(-2*normal(2,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(-xc(0,2) +
    xc(1,2)))*(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(xc(0,2) -
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*(xc(0,0) - xc(1,0))*(-xc(0,0) +
    xc(2,0)) + 2*(-xc(0,2) + xc(1,2))*(xc(0,2) - xc(2,2)))*(-(normal(0,0)*(xc(0,0) -
    xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(xc(0,2) -
    xc(2,2)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(xc(0,0) - xc(1,0)) -
    2*normal(0,0)*(-xc(0,2) + xc(1,2)))*((xc(0,2) - xc(2,2))*(xc(0,0) - xc(3,0)) +
    (-xc(0,0) + xc(2,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(7,5)=
    (4*normpowfour*(xc(0,0) - xc(3,0)) - 2*normsquare*(-2*normal(2,0)*(xc(0,0)
    - xc(1,0)) - 2*normal(0,0)*(-xc(0,2) + xc(1,2)))*((-xc(0,1) + xc(2,1))*(xc(0,0) -
    xc(3,0)) + (xc(0,0) - xc(2,0))*(xc(0,1) - xc(3,1))) + 3*(-2*normal(2,0)*(xc(0,0) -
    xc(1,0)) - 2*normal(0,0)*(-xc(0,2) + xc(1,2)))*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) -
    2*normal(0,0)*(-xc(0,1) + xc(2,1)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    4*normsquare*(xc(0,2)*(xc(1,1) - 2*xc(2,1)) + 2*xc(1,2)*xc(2,1) -
    xc(1,1)*xc(2,2) + xc(0,1)*(xc(0,2) - 2*xc(1,2) + xc(2,2)))*(-(normal(0,0)*(xc(0,0)
    - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(-xc(0,1) +
    xc(2,1)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(7,6)=
    (-4*normsquare*(xc(0,0) - xc(1,0))*(-xc(0,1) + xc(1,1))*(-(normal(0,0)*(xc(0,0) -
    xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) +
    3*(-2*normal(2,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(-xc(0,2) +
    xc(1,2)))*(-2*normal(2,0)*(-xc(0,1) + xc(1,1)) + 2*normal(1,0)*(-xc(0,2) +
    xc(1,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(-xc(0,1) + xc(1,1))
    + 2*normal(1,0)*(-xc(0,2) + xc(1,2)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) +
    (xc(0,0) - xc(1,0))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal(2,0)*(xc(0,0) -
    xc(1,0)) - 2*normal(0,0)*(-xc(0,2) + xc(1,2)))*((xc(0,2) - xc(1,2))*(xc(0,1) -
    xc(3,1)) + (-xc(0,1) + xc(1,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(7,7)=
    (-4*normsquare*(pow(xc(0,0) - xc(1,0),2) + pow(xc(0,2) -
    xc(1,2),2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) + 3*pow(-2*normal(2,0)*(xc(0,0) - xc(1,0)) -
    2*normal(0,0)*(-xc(0,2) + xc(1,2)),2)*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-2*normal(2,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(-xc(0,2) +
    xc(1,2)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(7,8)=
    (-2*normsquare*(-2*normal(2,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(-xc(0,2) +
    xc(1,2)))*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(-xc(0,1) + xc(3,1))) - 4*normsquare*(xc(0,1) - xc(1,1))*(-xc(0,2) +
    xc(1,2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) + 3*(2*normal(1,0)*(xc(0,0) - xc(1,0)) -
    2*normal(0,0)*(xc(0,1) - xc(1,1)))*(-2*normal(2,0)*(xc(0,0) - xc(1,0)) -
    2*normal(0,0)*(-xc(0,2) + xc(1,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal(1,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(xc(0,1) -
    xc(1,1)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(7,9)=
    (normal(1,0)*(xc(0,0)*pow(xc(1,1),2) + xc(0,0)*pow(xc(1,2),2) +
    pow(xc(0,1),2)*(xc(1,0) - xc(2,0)) + pow(xc(0,2),2)*(xc(1,0) - xc(2,0)) -
    pow(xc(1,1),2)*xc(2,0) - pow(xc(1,2),2)*xc(2,0) - xc(0,0)*xc(1,1)*xc(2,1) +
    xc(1,0)*xc(1,1)*xc(2,1) - xc(0,1)*(-2*xc(1,1)*xc(2,0) + xc(0,0)*(xc(1,1) -
    xc(2,1)) + xc(1,0)*(xc(1,1) + xc(2,1))) - xc(0,0)*xc(1,2)*xc(2,2) +
    xc(1,0)*xc(1,2)*xc(2,2) - xc(0,2)*(-2*xc(1,2)*xc(2,0) + xc(0,0)*(xc(1,2) -
    xc(2,2)) + xc(1,0)*(xc(1,2) + xc(2,2)))))/normcube;

  elematrix(7,10)=
    -(normal(1,0)*(-2*normal(2,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(-xc(0,2) +
    xc(1,2))))/(2.*normcube);

  elematrix(7,11)=
    (normal(1,0)*(pow(xc(0,1),2)*xc(1,2) - xc(0,1)*xc(1,1)*xc(1,2) +
    xc(1,0)*xc(1,2)*xc(2,0) + xc(0,2)*(pow(xc(1,0),2) - xc(1,0)*xc(2,0) - (xc(0,1)
    - xc(1,1))*(xc(1,1) - xc(2,1))) - xc(0,1)*xc(1,2)*xc(2,1) +
    xc(1,1)*xc(1,2)*xc(2,1) - xc(0,0)*(xc(0,2)*(xc(1,0) - xc(2,0)) + xc(1,2)*xc(2,0)
    + xc(1,0)*(xc(1,2) - 2*xc(2,2))) + pow(xc(0,0),2)*(xc(1,2) - xc(2,2)) -
    pow(xc(0,1),2)*xc(2,2) - pow(xc(1,0),2)*xc(2,2) + 2*xc(0,1)*xc(1,1)*xc(2,2)
    - pow(xc(1,1),2)*xc(2,2)))/normcube;

  elematrix(8,0)=
    (4*normpowfour*(-xc(1,1) + xc(3,1)) -
    2*normsquare*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) + 2*normal(1,0)*(-xc(1,2) +
    xc(2,2)))*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(-xc(0,1) + xc(3,1))) - 4*normsquare*(xc(1,0)*xc(1,2) +
    xc(0,2)*(xc(1,0) - xc(2,0)) + xc(1,2)*xc(2,0) - 2*xc(0,0)*(xc(1,2) - xc(2,2)) -
    2*xc(1,0)*xc(2,2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) +
    xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) + 3*(2*normal(1,0)*(xc(0,0) - xc(1,0)) -
    2*normal(0,0)*(xc(0,1) - xc(1,1)))*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) +
    2*normal(1,0)*(-xc(1,2) + xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal(1,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(xc(0,1) -
    xc(1,1)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) +
    xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(8,1)=
    (4*normpowfour*(xc(1,0) - xc(3,0)) - 2*normsquare*(-2*normal(2,0)*(xc(1,0)
    - xc(2,0)) - 2*normal(0,0)*(-xc(1,2) + xc(2,2)))*((xc(0,1) - xc(1,1))*(xc(0,0) -
    xc(3,0)) + (xc(0,0) - xc(1,0))*(-xc(0,1) + xc(3,1))) -
    4*normsquare*(xc(1,1)*xc(1,2) + xc(0,2)*(xc(1,1) - xc(2,1)) + xc(1,2)*xc(2,1) -
    2*xc(0,1)*(xc(1,2) - xc(2,2)) - 2*xc(1,1)*xc(2,2))*(-(normal(0,0)*(xc(0,0) -
    xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) +
    3*(2*normal(1,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(xc(0,1) -
    xc(1,1)))*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(-xc(1,2) +
    xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal(1,0)*(xc(0,0) - xc(1,0)) -
    2*normal(0,0)*(xc(0,1) - xc(1,1)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) +
    xc(3,0)) + xc(1,0)*(xc(2,2) - xc(3,2)) +
    xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(8,2)=
    (-2*normsquare*(2*normal(1,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(xc(1,1) -
    xc(2,1)))*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(-xc(0,1) + xc(3,1))) - 2*normsquare*(2*normal(1,0)*(xc(0,0) - xc(1,0)) -
    2*normal(0,0)*(xc(0,1) - xc(1,1)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0)
    - xc(2,0)*xc(3,1) + xc(1,0)*(-xc(2,1) + xc(3,1))) + 3*(2*normal(1,0)*(xc(0,0) -
    xc(1,0)) - 2*normal(0,0)*(xc(0,1) - xc(1,1)))*(2*normal(1,0)*(xc(1,0) - xc(2,0)) -
    2*normal(0,0)*(xc(1,1) - xc(2,1)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*(xc(0,0) - xc(1,0))*(xc(1,0) - xc(2,0)) + 2*(xc(0,1) -
    xc(1,1))*(xc(1,1) - xc(2,1)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) -
    xc(3,2))))/(4.*normpowfive);

  elematrix(8,3)=
    (4*normpowfour*(xc(0,1) - xc(3,1)) - 2*normsquare*(-2*normal(2,0)*(xc(0,1)
    - xc(2,1)) + 2*normal(1,0)*(xc(0,2) - xc(2,2)))*((xc(0,1) - xc(1,1))*(xc(0,0) -
    xc(3,0)) + (xc(0,0) - xc(1,0))*(-xc(0,1) + xc(3,1))) + 3*(2*normal(1,0)*(xc(0,0) -
    xc(1,0)) - 2*normal(0,0)*(xc(0,1) - xc(1,1)))*(-2*normal(2,0)*(xc(0,1) - xc(2,1)) +
    2*normal(1,0)*(xc(0,2) - xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-(xc(1,2)*xc(2,0)) + xc(0,2)*(-2*xc(1,0) + xc(2,0)) +
    xc(0,0)*(xc(0,2) + xc(1,2) - 2*xc(2,2)) +
    2*xc(1,0)*xc(2,2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) +
    xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal(1,0)*(xc(0,0) -
    xc(1,0)) - 2*normal(0,0)*(xc(0,1) - xc(1,1)))*((xc(0,2) - xc(2,2))*(-xc(0,1) +
    xc(3,1)) + (xc(0,1) - xc(2,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(8,4)=
    (4*normpowfour*(-xc(0,0) + xc(3,0)) -
    2*normsquare*(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(xc(0,2) -
    xc(2,2)))*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(-xc(0,1) + xc(3,1))) + 3*(2*normal(1,0)*(xc(0,0) - xc(1,0)) -
    2*normal(0,0)*(xc(0,1) - xc(1,1)))*(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) -
    2*normal(0,0)*(xc(0,2) - xc(2,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-(xc(1,2)*xc(2,1)) + xc(0,2)*(-2*xc(1,1) + xc(2,1)) +
    xc(0,1)*(xc(0,2) + xc(1,2) - 2*xc(2,2)) +
    2*xc(1,1)*xc(2,2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) +
    xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal(1,0)*(xc(0,0) -
    xc(1,0)) - 2*normal(0,0)*(xc(0,1) - xc(1,1)))*((xc(0,2) - xc(2,2))*(xc(0,0) -
    xc(3,0)) + (-xc(0,0) + xc(2,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(8,5)=
    (-2*normsquare*(2*normal(1,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(xc(0,1) -
    xc(1,1)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) - 2*normsquare*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) -
    2*normal(0,0)*(-xc(0,1) + xc(2,1)))*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) +
    (xc(0,0) - xc(1,0))*(-xc(0,1) + xc(3,1))) + 3*(2*normal(1,0)*(xc(0,0) - xc(1,0)) -
    2*normal(0,0)*(xc(0,1) - xc(1,1)))*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) -
    2*normal(0,0)*(-xc(0,1) + xc(2,1)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*(xc(0,0) - xc(1,0))*(-xc(0,0) + xc(2,0)) + 2*(xc(0,1) -
    xc(1,1))*(-xc(0,1) + xc(2,1)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) -
    xc(3,2))))/(4.*normpowfive);

  elematrix(8,6)=
    (-2*normsquare*(-2*normal(2,0)*(-xc(0,1) + xc(1,1)) + 2*normal(1,0)*(-xc(0,2) +
    xc(1,2)))*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(-xc(0,1) + xc(3,1))) - 4*normsquare*(xc(0,0) - xc(1,0))*(-xc(0,2) +
    xc(1,2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) + 3*(2*normal(1,0)*(xc(0,0) - xc(1,0)) -
    2*normal(0,0)*(xc(0,1) - xc(1,1)))*(-2*normal(2,0)*(-xc(0,1) + xc(1,1)) +
    2*normal(1,0)*(-xc(0,2) + xc(1,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal(1,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(xc(0,1) -
    xc(1,1)))*((xc(0,2) - xc(1,2))*(xc(0,1) - xc(3,1)) + (-xc(0,1) +
    xc(1,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(8,7)=
    (-2*normsquare*(-2*normal(2,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(-xc(0,2) +
    xc(1,2)))*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(-xc(0,1) + xc(3,1))) - 4*normsquare*(xc(0,1) - xc(1,1))*(-xc(0,2) +
    xc(1,2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))) + 3*(2*normal(1,0)*(xc(0,0) - xc(1,0)) -
    2*normal(0,0)*(xc(0,1) - xc(1,1)))*(-2*normal(2,0)*(xc(0,0) - xc(1,0)) -
    2*normal(0,0)*(-xc(0,2) + xc(1,2)))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal(1,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(xc(0,1) -
    xc(1,1)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(8,8)=
    (-4*normsquare*(2*normal(1,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(xc(0,1) -
    xc(1,1)))*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(-xc(0,1) + xc(3,1))) + 3*pow(2*normal(1,0)*(xc(0,0) - xc(1,0)) -
    2*normal(0,0)*(xc(0,1) - xc(1,1)),2)*(-(normal(0,0)*(xc(0,0) - xc(3,0))) +
    normal(1,0)*(-xc(0,1) + xc(3,1)) - normal(2,0)*(xc(0,2) - xc(3,2))) -
    4*normsquare*(pow(xc(0,0) - xc(1,0),2) + pow(xc(0,1) -
    xc(1,1),2))*(-(normal(0,0)*(xc(0,0) - xc(3,0))) + normal(1,0)*(-xc(0,1) + xc(3,1)) -
    normal(2,0)*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(8,9)=
    -((normal(2,0)*(-(xc(0,0)*pow(xc(1,1),2)) - xc(0,0)*pow(xc(1,2),2) +
    pow(xc(1,1),2)*xc(2,0) + pow(xc(1,2),2)*xc(2,0) + pow(xc(0,1),2)*(-xc(1,0)
    + xc(2,0)) + pow(xc(0,2),2)*(-xc(1,0) + xc(2,0)) + xc(0,0)*xc(1,1)*xc(2,1) -
    xc(1,0)*xc(1,1)*xc(2,1) + xc(0,1)*(-2*xc(1,1)*xc(2,0) + xc(0,0)*(xc(1,1) -
    xc(2,1)) + xc(1,0)*(xc(1,1) + xc(2,1))) + xc(0,0)*xc(1,2)*xc(2,2) -
    xc(1,0)*xc(1,2)*xc(2,2) + xc(0,2)*(-2*xc(1,2)*xc(2,0) + xc(0,0)*(xc(1,2) -
    xc(2,2)) + xc(1,0)*(xc(1,2) + xc(2,2)))))/normcube);

  elematrix(8,10)=
    (normal(2,0)*(pow(xc(0,2),2)*xc(1,1) - xc(0,2)*xc(1,1)*xc(1,2) +
    xc(1,0)*xc(1,1)*xc(2,0) - xc(0,0)*(xc(0,1)*(xc(1,0) - xc(2,0)) + xc(1,1)*xc(2,0)
    + xc(1,0)*(xc(1,1) - 2*xc(2,1))) + pow(xc(0,0),2)*(xc(1,1) - xc(2,1)) -
    pow(xc(0,2),2)*xc(2,1) - pow(xc(1,0),2)*xc(2,1) + 2*xc(0,2)*xc(1,2)*xc(2,1)
    - pow(xc(1,2),2)*xc(2,1) + xc(0,1)*(pow(xc(1,0),2) - xc(1,0)*xc(2,0) -
    (xc(0,2) - xc(1,2))*(xc(1,2) - xc(2,2))) - xc(0,2)*xc(1,1)*xc(2,2) +
    xc(1,1)*xc(1,2)*xc(2,2)))/normcube;

  elematrix(8,11)=
    -((2*normal(1,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(xc(0,1) -
    xc(1,1)))*(-(xc(1,1)*xc(2,0)) + xc(0,1)*(-xc(1,0) + xc(2,0)) + xc(0,0)*(xc(1,1)
    - xc(2,1)) + xc(1,0)*xc(2,1)))/(2.*normcube);

  elematrix(9,0)=
    -((-(xc(1,2)*xc(2,1)) + xc(0,2)*(-xc(1,1) + xc(2,1)) + xc(0,1)*(xc(1,2) -
    xc(2,2)) + xc(1,1)*xc(2,2))*(-2*normal(2,0)*(-xc(1,1) + xc(2,1)) +
    2*normal(1,0)*(-xc(1,2) + xc(2,2))))/(2.*normcube);

  elematrix(9,1)=
    (normal(1,0)*(xc(0,2)*xc(1,0)*xc(1,2) + pow(xc(1,1),2)*xc(2,0) -
    xc(0,2)*xc(1,2)*xc(2,0) + pow(xc(1,2),2)*xc(2,0) + xc(0,1)*(xc(1,0) -
    xc(2,0))*(xc(1,1) - xc(2,1)) - xc(1,0)*xc(1,1)*xc(2,1) - xc(1,1)*xc(2,0)*xc(2,1)
    + xc(1,0)*pow(xc(2,1),2) - xc(0,2)*xc(1,0)*xc(2,2) - xc(1,0)*xc(1,2)*xc(2,2) +
    xc(0,2)*xc(2,0)*xc(2,2) - xc(1,2)*xc(2,0)*xc(2,2) + xc(1,0)*pow(xc(2,2),2) -
    xc(0,0)*(pow(xc(1,1),2) + pow(xc(1,2),2) - 2*xc(1,1)*xc(2,1) +
    pow(xc(2,1),2) - 2*xc(1,2)*xc(2,2) + pow(xc(2,2),2))))/normcube;

  elematrix(9,2)=
    (normal(2,0)*(xc(0,2)*xc(1,0)*xc(1,2) + pow(xc(1,1),2)*xc(2,0) -
    xc(0,2)*xc(1,2)*xc(2,0) + pow(xc(1,2),2)*xc(2,0) + xc(0,1)*(xc(1,0) -
    xc(2,0))*(xc(1,1) - xc(2,1)) - xc(1,0)*xc(1,1)*xc(2,1) - xc(1,1)*xc(2,0)*xc(2,1)
    + xc(1,0)*pow(xc(2,1),2) - xc(0,2)*xc(1,0)*xc(2,2) - xc(1,0)*xc(1,2)*xc(2,2) +
    xc(0,2)*xc(2,0)*xc(2,2) - xc(1,2)*xc(2,0)*xc(2,2) + xc(1,0)*pow(xc(2,2),2) -
    xc(0,0)*(pow(xc(1,1),2) + pow(xc(1,2),2) - 2*xc(1,1)*xc(2,1) +
    pow(xc(2,1),2) - 2*xc(1,2)*xc(2,2) + pow(xc(2,2),2))))/normcube;

  elematrix(9,3)=
    -((-2*normal(2,0)*(xc(0,1) - xc(2,1)) + 2*normal(1,0)*(xc(0,2) -
    xc(2,2)))*(-(xc(1,2)*xc(2,1)) + xc(0,2)*(-xc(1,1) + xc(2,1)) + xc(0,1)*(xc(1,2)
    - xc(2,2)) + xc(1,1)*xc(2,2)))/(2.*normcube);

  elematrix(9,4)=
    -((normal(1,0)*(pow(xc(0,1),2)*(xc(1,0) - xc(2,0)) + pow(xc(0,2),2)*(xc(1,0) -
    xc(2,0)) + xc(0,0)*xc(1,1)*xc(2,1) - xc(1,1)*xc(2,0)*xc(2,1) -
    xc(0,0)*pow(xc(2,1),2) + xc(1,0)*pow(xc(2,1),2) + xc(0,1)*(xc(1,1)*xc(2,0) +
    (-2*xc(1,0) + xc(2,0))*xc(2,1) + xc(0,0)*(-xc(1,1) + xc(2,1))) +
    xc(0,0)*xc(1,2)*xc(2,2) - xc(1,2)*xc(2,0)*xc(2,2) - xc(0,0)*pow(xc(2,2),2) +
    xc(1,0)*pow(xc(2,2),2) + xc(0,2)*(xc(1,2)*xc(2,0) + (-2*xc(1,0) +
    xc(2,0))*xc(2,2) + xc(0,0)*(-xc(1,2) + xc(2,2)))))/normcube);

  elematrix(9,5)=
    -((normal(2,0)*(pow(xc(0,1),2)*(xc(1,0) - xc(2,0)) + pow(xc(0,2),2)*(xc(1,0) -
    xc(2,0)) + xc(0,0)*xc(1,1)*xc(2,1) - xc(1,1)*xc(2,0)*xc(2,1) -
    xc(0,0)*pow(xc(2,1),2) + xc(1,0)*pow(xc(2,1),2) + xc(0,1)*(xc(1,1)*xc(2,0) +
    (-2*xc(1,0) + xc(2,0))*xc(2,1) + xc(0,0)*(-xc(1,1) + xc(2,1))) +
    xc(0,0)*xc(1,2)*xc(2,2) - xc(1,2)*xc(2,0)*xc(2,2) - xc(0,0)*pow(xc(2,2),2) +
    xc(1,0)*pow(xc(2,2),2) + xc(0,2)*(xc(1,2)*xc(2,0) + (-2*xc(1,0) +
    xc(2,0))*xc(2,2) + xc(0,0)*(-xc(1,2) + xc(2,2)))))/normcube);

  elematrix(9,6)=
    -((-2*normal(2,0)*(-xc(0,1) + xc(1,1)) + 2*normal(1,0)*(-xc(0,2) +
    xc(1,2)))*(-(xc(1,2)*xc(2,1)) + xc(0,2)*(-xc(1,1) + xc(2,1)) + xc(0,1)*(xc(1,2)
    - xc(2,2)) + xc(1,1)*xc(2,2)))/(2.*normcube);

  elematrix(9,7)=
    (normal(1,0)*(xc(0,0)*pow(xc(1,1),2) + xc(0,0)*pow(xc(1,2),2) +
    pow(xc(0,1),2)*(xc(1,0) - xc(2,0)) + pow(xc(0,2),2)*(xc(1,0) - xc(2,0)) -
    pow(xc(1,1),2)*xc(2,0) - pow(xc(1,2),2)*xc(2,0) - xc(0,0)*xc(1,1)*xc(2,1) +
    xc(1,0)*xc(1,1)*xc(2,1) - xc(0,1)*(-2*xc(1,1)*xc(2,0) + xc(0,0)*(xc(1,1) -
    xc(2,1)) + xc(1,0)*(xc(1,1) + xc(2,1))) - xc(0,0)*xc(1,2)*xc(2,2) +
    xc(1,0)*xc(1,2)*xc(2,2) - xc(0,2)*(-2*xc(1,2)*xc(2,0) + xc(0,0)*(xc(1,2) -
    xc(2,2)) + xc(1,0)*(xc(1,2) + xc(2,2)))))/normcube;

  elematrix(9,8)=
    -((normal(2,0)*(-(xc(0,0)*pow(xc(1,1),2)) - xc(0,0)*pow(xc(1,2),2) +
    pow(xc(1,1),2)*xc(2,0) + pow(xc(1,2),2)*xc(2,0) + pow(xc(0,1),2)*(-xc(1,0)
    + xc(2,0)) + pow(xc(0,2),2)*(-xc(1,0) + xc(2,0)) + xc(0,0)*xc(1,1)*xc(2,1) -
    xc(1,0)*xc(1,1)*xc(2,1) + xc(0,1)*(-2*xc(1,1)*xc(2,0) + xc(0,0)*(xc(1,1) -
    xc(2,1)) + xc(1,0)*(xc(1,1) + xc(2,1))) + xc(0,0)*xc(1,2)*xc(2,2) -
    xc(1,0)*xc(1,2)*xc(2,2) + xc(0,2)*(-2*xc(1,2)*xc(2,0) + xc(0,0)*(xc(1,2) -
    xc(2,2)) + xc(1,0)*(xc(1,2) + xc(2,2)))))/normcube);

  elematrix(9,9)=0;

  elematrix(9,10)=0;

  elematrix(9,11)=0;

  elematrix(10,0)=
    (normal(0,0)*(xc(0,2)*xc(1,1)*xc(1,2) - xc(1,0)*xc(1,1)*xc(2,0) +
    xc(1,1)*pow(xc(2,0),2) + xc(0,0)*(xc(1,0) - xc(2,0))*(xc(1,1) - xc(2,1)) +
    pow(xc(1,0),2)*xc(2,1) - xc(0,2)*xc(1,2)*xc(2,1) + pow(xc(1,2),2)*xc(2,1) -
    xc(1,0)*xc(2,0)*xc(2,1) - xc(0,2)*xc(1,1)*xc(2,2) - xc(1,1)*xc(1,2)*xc(2,2) +
    xc(0,2)*xc(2,1)*xc(2,2) - xc(1,2)*xc(2,1)*xc(2,2) + xc(1,1)*pow(xc(2,2),2) -
    xc(0,1)*(pow(xc(1,0),2) + pow(xc(1,2),2) - 2*xc(1,0)*xc(2,0) +
    pow(xc(2,0),2) - 2*xc(1,2)*xc(2,2) + pow(xc(2,2),2))))/normcube;

  elematrix(10,1)=
    -(normal(1,0)*(-2*normal(2,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(-xc(1,2) +
    xc(2,2))))/(2.*normcube);

  elematrix(10,2)=
    -((normal(2,0)*(-(xc(0,2)*xc(1,1)*xc(1,2)) + xc(1,0)*xc(1,1)*xc(2,0) -
    xc(1,1)*pow(xc(2,0),2) - xc(0,0)*(xc(1,0) - xc(2,0))*(xc(1,1) - xc(2,1)) -
    pow(xc(1,0),2)*xc(2,1) + xc(0,2)*xc(1,2)*xc(2,1) - pow(xc(1,2),2)*xc(2,1) +
    xc(1,0)*xc(2,0)*xc(2,1) + xc(0,2)*xc(1,1)*xc(2,2) + xc(1,1)*xc(1,2)*xc(2,2) -
    xc(0,2)*xc(2,1)*xc(2,2) + xc(1,2)*xc(2,1)*xc(2,2) - xc(1,1)*pow(xc(2,2),2) +
    xc(0,1)*(pow(xc(1,0),2) + pow(xc(1,2),2) - 2*xc(1,0)*xc(2,0) +
    pow(xc(2,0),2) - 2*xc(1,2)*xc(2,2) +pow(xc(2,2),2))))/normcube);

  elematrix(10,3)=
    -((normal(0,0)*(xc(0,1)*xc(1,0)*xc(2,0) - xc(0,1)*pow(xc(2,0),2) +
    xc(1,1)*pow(xc(2,0),2) + pow(xc(0,0),2)*(xc(1,1) - xc(2,1)) +
    pow(xc(0,2),2)*(xc(1,1) - xc(2,1)) - xc(1,0)*xc(2,0)*xc(2,1) +
    xc(0,0)*(-2*xc(1,1)*xc(2,0) + xc(0,1)*(-xc(1,0) + xc(2,0)) + (xc(1,0) +
    xc(2,0))*xc(2,1)) + xc(0,1)*xc(1,2)*xc(2,2) - xc(1,2)*xc(2,1)*xc(2,2) -
    xc(0,1)*pow(xc(2,2),2) + xc(1,1)*pow(xc(2,2),2) + xc(0,2)*(xc(1,2)*xc(2,1) +
    (-2*xc(1,1) + xc(2,1))*xc(2,2) + xc(0,1)*(-xc(1,2) +xc(2,2)))))/normcube);

  elematrix(10,4)=
    -(normal(1,0)*(-2*normal(2,0)*(-xc(0,0) + xc(2,0)) - 2*normal(0,0)*(xc(0,2) -
    xc(2,2))))/(2.*normcube);

  elematrix(10,5)=
    -((normal(2,0)*(xc(0,1)*xc(1,0)*xc(2,0) - xc(0,1)*pow(xc(2,0),2) +
    xc(1,1)*pow(xc(2,0),2) + pow(xc(0,0),2)*(xc(1,1) - xc(2,1)) +
    pow(xc(0,2),2)*(xc(1,1) - xc(2,1)) - xc(1,0)*xc(2,0)*xc(2,1) +
    xc(0,0)*(-2*xc(1,1)*xc(2,0) + xc(0,1)*(-xc(1,0) + xc(2,0)) + (xc(1,0) +
    xc(2,0))*xc(2,1)) + xc(0,1)*xc(1,2)*xc(2,2) - xc(1,2)*xc(2,1)*xc(2,2) -
    xc(0,1)*pow(xc(2,2),2) + xc(1,1)*pow(xc(2,2),2) + xc(0,2)*(xc(1,2)*xc(2,1) +
    (-2*xc(1,1) + xc(2,1))*xc(2,2) + xc(0,1)*(-xc(1,2) +xc(2,2)))))/normcube);

  elematrix(10,6)=
    (normal(0,0)*(pow(xc(0,2),2)*xc(1,1) - xc(0,2)*xc(1,1)*xc(1,2) +
    xc(1,0)*xc(1,1)*xc(2,0) - xc(0,0)*(xc(0,1)*(xc(1,0) - xc(2,0)) + xc(1,1)*xc(2,0)
    + xc(1,0)*(xc(1,1) - 2*xc(2,1))) + pow(xc(0,0),2)*(xc(1,1) - xc(2,1)) -
    pow(xc(0,2),2)*xc(2,1) - pow(xc(1,0),2)*xc(2,1) + 2*xc(0,2)*xc(1,2)*xc(2,1)
    - pow(xc(1,2),2)*xc(2,1) + xc(0,1)*(pow(xc(1,0),2) - xc(1,0)*xc(2,0) -
    (xc(0,2) - xc(1,2))*(xc(1,2) - xc(2,2))) - xc(0,2)*xc(1,1)*xc(2,2) +
    xc(1,1)*xc(1,2)*xc(2,2)))/normcube;

  elematrix(10,7)=
    -(normal(1,0)*(-2*normal(2,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(-xc(0,2) +
    xc(1,2))))/(2.*normcube);

  elematrix(10,8)=
    (normal(2,0)*(pow(xc(0,2),2)*xc(1,1) - xc(0,2)*xc(1,1)*xc(1,2) +
    xc(1,0)*xc(1,1)*xc(2,0) - xc(0,0)*(xc(0,1)*(xc(1,0) - xc(2,0)) + xc(1,1)*xc(2,0)
    + xc(1,0)*(xc(1,1) - 2*xc(2,1))) + pow(xc(0,0),2)*(xc(1,1) - xc(2,1)) -
    pow(xc(0,2),2)*xc(2,1) - pow(xc(1,0),2)*xc(2,1) + 2*xc(0,2)*xc(1,2)*xc(2,1)
    - pow(xc(1,2),2)*xc(2,1) + xc(0,1)*(pow(xc(1,0),2) - xc(1,0)*xc(2,0) -
    (xc(0,2) - xc(1,2))*(xc(1,2) - xc(2,2))) - xc(0,2)*xc(1,1)*xc(2,2) +
    xc(1,1)*xc(1,2)*xc(2,2)))/normcube;

  elematrix(10,9)=0;

  elematrix(10,10)=0;

  elematrix(10,11)=0;

  elematrix(11,0)=
    -((normal(0,0)*(-(xc(0,1)*xc(1,1)*xc(1,2)) + xc(1,0)*xc(1,2)*xc(2,0) -
    xc(1,2)*pow(xc(2,0),2) + xc(0,1)*xc(1,2)*xc(2,1) + xc(1,1)*xc(1,2)*xc(2,1) -
    xc(1,2)*pow(xc(2,1),2) + xc(0,2)*(pow(xc(1,0),2) + pow(xc(1,1),2) -
    2*xc(1,0)*xc(2,0) + pow(xc(2,0),2) - 2*xc(1,1)*xc(2,1) + pow(xc(2,1),2)) -
    xc(0,0)*(xc(1,0) - xc(2,0))*(xc(1,2) - xc(2,2)) - pow(xc(1,0),2)*xc(2,2) +
    xc(0,1)*xc(1,1)*xc(2,2) - pow(xc(1,1),2)*xc(2,2) + xc(1,0)*xc(2,0)*xc(2,2) -
    xc(0,1)*xc(2,1)*xc(2,2) + xc(1,1)*xc(2,1)*xc(2,2)))/normcube);

  elematrix(11,1)=
    -((normal(1,0)*(-(xc(0,1)*xc(1,1)*xc(1,2)) + xc(1,0)*xc(1,2)*xc(2,0) -
    xc(1,2)*pow(xc(2,0),2) + xc(0,1)*xc(1,2)*xc(2,1) + xc(1,1)*xc(1,2)*xc(2,1) -
    xc(1,2)*pow(xc(2,1),2) + xc(0,2)*(pow(xc(1,0),2) + pow(xc(1,1),2) -
    2*xc(1,0)*xc(2,0) + pow(xc(2,0),2) - 2*xc(1,1)*xc(2,1) + pow(xc(2,1),2)) -
    xc(0,0)*(xc(1,0) - xc(2,0))*(xc(1,2) - xc(2,2)) - pow(xc(1,0),2)*xc(2,2) +
    xc(0,1)*xc(1,1)*xc(2,2) - pow(xc(1,1),2)*xc(2,2) + xc(1,0)*xc(2,0)*xc(2,2) -
    xc(0,1)*xc(2,1)*xc(2,2) + xc(1,1)*xc(2,1)*xc(2,2)))/normcube);

  elematrix(11,2)=
    -((2*normal(1,0)*(xc(1,0) - xc(2,0)) - 2*normal(0,0)*(xc(1,1) -
    xc(2,1)))*(-(xc(1,1)*xc(2,0)) + xc(0,1)*(-xc(1,0) + xc(2,0)) + xc(0,0)*(xc(1,1)
    - xc(2,1)) + xc(1,0)*xc(2,1)))/(2.*normcube);

  elematrix(11,3)=
    -((normal(0,0)*(xc(0,2)*xc(1,0)*xc(2,0) - xc(0,2)*pow(xc(2,0),2) +
    xc(1,2)*pow(xc(2,0),2) + xc(0,2)*xc(1,1)*xc(2,1) - xc(0,2)*pow(xc(2,1),2) +
    xc(1,2)*pow(xc(2,1),2) + pow(xc(0,0),2)*(xc(1,2) - xc(2,2)) +
    pow(xc(0,1),2)*(xc(1,2) - xc(2,2)) - xc(1,0)*xc(2,0)*xc(2,2) -
    xc(1,1)*xc(2,1)*xc(2,2) + xc(0,0)*(-2*xc(1,2)*xc(2,0) + xc(0,2)*(-xc(1,0) +
    xc(2,0)) + (xc(1,0) + xc(2,0))*xc(2,2)) + xc(0,1)*(-2*xc(1,2)*xc(2,1) +
    xc(0,2)*(-xc(1,1) + xc(2,1)) + (xc(1,1) +xc(2,1))*xc(2,2))))/normcube);

  elematrix(11,4)=
    -((normal(1,0)*(xc(0,2)*xc(1,0)*xc(2,0) - xc(0,2)*pow(xc(2,0),2) +
    xc(1,2)*pow(xc(2,0),2) + xc(0,2)*xc(1,1)*xc(2,1) - xc(0,2)*pow(xc(2,1),2) +
    xc(1,2)*pow(xc(2,1),2) + pow(xc(0,0),2)*(xc(1,2) - xc(2,2)) +
    pow(xc(0,1),2)*(xc(1,2) - xc(2,2)) - xc(1,0)*xc(2,0)*xc(2,2) -
    xc(1,1)*xc(2,1)*xc(2,2) + xc(0,0)*(-2*xc(1,2)*xc(2,0) + xc(0,2)*(-xc(1,0) +
    xc(2,0)) + (xc(1,0) + xc(2,0))*xc(2,2)) + xc(0,1)*(-2*xc(1,2)*xc(2,1) +
    xc(0,2)*(-xc(1,1) + xc(2,1)) + (xc(1,1) +xc(2,1))*xc(2,2))))/normcube);

  elematrix(11,5)=
    -((-(xc(1,1)*xc(2,0)) + xc(0,1)*(-xc(1,0) + xc(2,0)) + xc(0,0)*(xc(1,1) -
    xc(2,1)) + xc(1,0)*xc(2,1))*(2*normal(1,0)*(-xc(0,0) + xc(2,0)) -
    2*normal(0,0)*(-xc(0,1) + xc(2,1))))/(2.*normcube);

  elematrix(11,6)=
    ((-(xc(1,2)*xc(2,1)) + xc(0,2)*(-xc(1,1) + xc(2,1)) + xc(0,1)*(xc(1,2) -
    xc(2,2)) + xc(1,1)*xc(2,2))*(pow(xc(0,1),2)*xc(1,2) - xc(0,1)*xc(1,1)*xc(1,2)
    + xc(1,0)*xc(1,2)*xc(2,0) + xc(0,2)*(pow(xc(1,0),2) - xc(1,0)*xc(2,0) -
    (xc(0,1) - xc(1,1))*(xc(1,1) - xc(2,1))) - xc(0,1)*xc(1,2)*xc(2,1) +
    xc(1,1)*xc(1,2)*xc(2,1) - xc(0,0)*(xc(0,2)*(xc(1,0) - xc(2,0)) + xc(1,2)*xc(2,0)
    + xc(1,0)*(xc(1,2) - 2*xc(2,2))) + pow(xc(0,0),2)*(xc(1,2) - xc(2,2)) -
    pow(xc(0,1),2)*xc(2,2) - pow(xc(1,0),2)*xc(2,2) + 2*xc(0,1)*xc(1,1)*xc(2,2)
    - pow(xc(1,1),2)*xc(2,2)))/normcube;

  elematrix(11,7)=
    (normal(1,0)*(pow(xc(0,1),2)*xc(1,2) - xc(0,1)*xc(1,1)*xc(1,2) +
    xc(1,0)*xc(1,2)*xc(2,0) + xc(0,2)*(pow(xc(1,0),2) - xc(1,0)*xc(2,0) - (xc(0,1)
    - xc(1,1))*(xc(1,1) - xc(2,1))) - xc(0,1)*xc(1,2)*xc(2,1) +
    xc(1,1)*xc(1,2)*xc(2,1) - xc(0,0)*(xc(0,2)*(xc(1,0) - xc(2,0)) + xc(1,2)*xc(2,0)
    + xc(1,0)*(xc(1,2) - 2*xc(2,2))) + pow(xc(0,0),2)*(xc(1,2) - xc(2,2)) -
    pow(xc(0,1),2)*xc(2,2) - pow(xc(1,0),2)*xc(2,2) + 2*xc(0,1)*xc(1,1)*xc(2,2)
    - pow(xc(1,1),2)*xc(2,2)))/normcube;

  elematrix(11,8)=
    -((2*normal(1,0)*(xc(0,0) - xc(1,0)) - 2*normal(0,0)*(xc(0,1) -
    xc(1,1)))*(-(xc(1,1)*xc(2,0)) + xc(0,1)*(-xc(1,0) + xc(2,0)) + xc(0,0)*(xc(1,1)
    - xc(2,1)) + xc(1,0)*xc(2,1)))/(2.*normcube);

  elematrix(11,9)=0;

  elematrix(11,10)=0;

  elematrix(11,11)=0;
  return;
}

/*----------------------------------------------------------------------*
 * second derivatives */
void DRT::ELEMENTS::ConstraintElement::ComputeSecondDerivDist2D
(
  const LINALG::Matrix<3,2>& xc,
  Epetra_SerialDenseMatrix& elematrix,
  const LINALG::Matrix<2,1>& normal
)
{

  double normsquare=pow(normal.Norm2(),2);
  double normcube=pow(normal.Norm2(),3);
  double normpowfour=pow(normal.Norm2(),4);
  double normpowfive=pow(normal.Norm2(),5);

  elematrix(0,0)=
  (normal(0,0)*(-2*pow(xc(1,0),3) - 2*xc(1,0)*pow(xc(1,1),2) -
  2*pow(xc(0,0),2)*(xc(1,0) - xc(2,0)) + pow(xc(0,1),2)*(xc(1,0) - xc(2,0)) +
  2*pow(xc(1,0),2)*xc(2,0) - pow(xc(1,1),2)*xc(2,0) +
  xc(0,1)*(2*xc(1,1)*xc(2,0) + xc(1,0)*(xc(1,1) - 3*xc(2,1))) +
  3*xc(1,0)*xc(1,1)*xc(2,1) + xc(0,0)*(4*pow(xc(1,0),2) - 4*xc(1,0)*xc(2,0) +
  3*normal(0,0)*(-xc(1,1) + xc(2,1)))))/normpowfive;

  elematrix(0,1)=
  (3*pow(normal(0,0),2)*normal(1,0)*(xc(0,0) - xc(2,0)) +
  normsquare*normal(1,0)*(-xc(1,0) + xc(2,0)) +
  normal(0,0)*(3*pow(normal(1,0),2)*(xc(0,1) - xc(2,1)) + normsquare*(-xc(1,1) +
  xc(2,1))))/normpowfive;

  elematrix(0,2)=
  (normal(0,0)*(pow(xc(0,0),3) + pow(xc(1,0),3) + xc(1,0)*pow(xc(1,1),2) -
  2*pow(xc(1,0),2)*xc(2,0) + pow(xc(1,1),2)*xc(2,0) +
  pow(xc(0,1),2)*(-2*xc(1,0) + xc(2,0)) - pow(xc(0,0),2)*(xc(1,0) + 2*xc(2,0))
  - 3*xc(1,0)*xc(1,1)*xc(2,1) + xc(0,0)*(pow(xc(0,1),2) - pow(xc(1,0),2) -
  2*pow(xc(1,1),2) + 4*xc(1,0)*xc(2,0) + xc(0,1)*(xc(1,1) - 3*xc(2,1)) +
  3*xc(1,1)*xc(2,1)) + xc(0,1)*(-2*xc(1,1)*xc(2,0) + xc(1,0)*(xc(1,1) +
  3*xc(2,1)))))/normpowfive;

  elematrix(0,3)=
  (normpowfour + normsquare*(normal(1,0)*(xc(0,0) - xc(2,0)) +
  normal(0,0)*(xc(1,1) - xc(2,1))) + 3*normal(0,0)*normal(1,0)*(normal(0,0)*(-xc(0,0) +
  xc(2,0)) + normal(1,0)*(-xc(0,1) + xc(2,1))))/normpowfive;

  elematrix(0,4)=(normal(0,0)*normal(1,0))/normcube;

  elematrix(0,5)=-(pow(normal(0,0),2)/normcube);

  elematrix(1,0)=
  (3*pow(normal(0,0),2)*normal(1,0)*(xc(0,0) - xc(2,0)) +
  normsquare*normal(1,0)*(-xc(1,0) + xc(2,0)) +
  normal(0,0)*(3*pow(normal(1,0),2)*(xc(0,1) - xc(2,1)) + normsquare*(-xc(1,1) +
  xc(2,1))))/normpowfive;

  elematrix(1,1)=
  (normal(1,0)*(-2*pow(xc(1,0),2)*xc(1,1) - 2*pow(xc(1,1),3) +
  3*xc(1,0)*xc(1,1)*xc(2,0) + xc(0,1)*(3*pow(xc(1,0),2) - 3*xc(1,0)*xc(2,0) +
  4*xc(1,1)*(xc(1,1) - xc(2,1))) + pow(xc(0,0),2)*(xc(1,1) - xc(2,1)) -
  2*pow(xc(0,1),2)*(xc(1,1) - xc(2,1)) - pow(xc(1,0),2)*xc(2,1) +
  2*pow(xc(1,1),2)*xc(2,1) + xc(0,0)*(-3*xc(0,1)*(xc(1,0) - xc(2,0)) -
  3*xc(1,1)*xc(2,0) + xc(1,0)*(xc(1,1) + 2*xc(2,1)))))/normpowfive;

  elematrix(1,2)=
  (-normpowfour + normsquare*(normal(1,0)*(xc(1,0) - xc(2,0)) +
  normal(0,0)*(xc(0,1) - xc(2,1))) + 3*normal(0,0)*normal(1,0)*(normal(0,0)*(-xc(0,0) +
  xc(2,0)) + normal(1,0)*(-xc(0,1) + xc(2,1))))/normpowfive;

  elematrix(1,3)=
  (normal(1,0)*(pow(xc(0,1),3) + pow(xc(1,0),2)*xc(1,1) + pow(xc(1,1),3) -
  3*xc(1,0)*xc(1,1)*xc(2,0) - xc(0,1)*(2*pow(xc(1,0),2) - 3*xc(1,0)*xc(2,0) +
  xc(1,1)*(xc(1,1) - 4*xc(2,1))) + xc(0,0)*(xc(0,1)*(xc(1,0) - 3*xc(2,0)) +
  3*xc(1,1)*xc(2,0) + xc(1,0)*(xc(1,1) - 2*xc(2,1))) + pow(xc(1,0),2)*xc(2,1) -
  2*pow(xc(1,1),2)*xc(2,1) + pow(xc(0,0),2)*(xc(0,1) - 2*xc(1,1) + xc(2,1)) -
  pow(xc(0,1),2)*(xc(1,1) + 2*xc(2,1))))/normpowfive;

  elematrix(1,4)=pow(normal(1,0),2)/normcube;

  elematrix(1,5)=-((normal(0,0)*normal(1,0))/normcube);

  elematrix(2,0)=
  (normal(0,0)*(pow(xc(0,0),3) + pow(xc(1,0),3) + xc(1,0)*pow(xc(1,1),2) -
  2*pow(xc(1,0),2)*xc(2,0) + pow(xc(1,1),2)*xc(2,0) +
  pow(xc(0,1),2)*(-2*xc(1,0) + xc(2,0)) - pow(xc(0,0),2)*(xc(1,0) + 2*xc(2,0))
  - 3*xc(1,0)*xc(1,1)*xc(2,1) + xc(0,0)*(pow(xc(0,1),2) - pow(xc(1,0),2) -
  2*pow(xc(1,1),2) + 4*xc(1,0)*xc(2,0) + xc(0,1)*(xc(1,1) - 3*xc(2,1)) +
  3*xc(1,1)*xc(2,1)) + xc(0,1)*(-2*xc(1,1)*xc(2,0) + xc(1,0)*(xc(1,1) +
  3*xc(2,1)))))/normpowfive;

  elematrix(2,1)=
  (-normpowfour + normsquare*(normal(1,0)*(xc(1,0) - xc(2,0)) +
  normal(0,0)*(xc(0,1) - xc(2,1))) + 3*normal(0,0)*normal(1,0)*(normal(0,0)*(-xc(0,0) +
  xc(2,0)) + normal(1,0)*(-xc(0,1) + xc(2,1))))/normpowfive;

  elematrix(2,2)=
  -((normal(0,0)*(2*pow(xc(0,0),3) - 2*pow(xc(1,0),2)*xc(2,0) +
  pow(xc(1,1),2)*xc(2,0) + pow(xc(0,1),2)*(-3*xc(1,0) + xc(2,0)) -
  2*pow(xc(0,0),2)*(2*xc(1,0) + xc(2,0)) - 3*xc(1,0)*xc(1,1)*xc(2,1) +
  xc(0,1)*(-2*xc(1,1)*xc(2,0) + 3*xc(1,0)*(xc(1,1) + xc(2,1))) +
  xc(0,0)*(2*pow(xc(0,1),2) + 2*pow(xc(1,0),2) - pow(xc(1,1),2) +
  4*xc(1,0)*xc(2,0) + 3*xc(1,1)*xc(2,1) - xc(0,1)*(xc(1,1) +
  3*xc(2,1)))))/normpowfive);

  elematrix(2,3)=
  (3*normal(0,0)*normal(1,0)*(normal(0,0)*(xc(0,0) - xc(2,0)) + normal(1,0)*(xc(0,1) -
  xc(2,1))) + normsquare*(normal(1,0)*(-xc(0,0) + xc(2,0)) + normal(0,0)*(-xc(0,1) +
  xc(2,1))))/normpowfive;

  elematrix(2,4)=-((normal(0,0)*normal(1,0))/normcube);

  elematrix(2,5)=pow(normal(0,0),2)/normcube;

  elematrix(3,0)=
  (normpowfour + normsquare*(normal(1,0)*(xc(0,0) - xc(2,0)) +
  normal(0,0)*(xc(1,1) - xc(2,1))) + 3*normal(0,0)*normal(1,0)*(normal(0,0)*(-xc(0,0) +
  xc(2,0)) + normal(1,0)*(-xc(0,1) + xc(2,1))))/normpowfive;

  elematrix(3,1)=
  (normal(1,0)*(pow(xc(0,1),3) + pow(xc(1,0),2)*xc(1,1) + pow(xc(1,1),3) -
  3*xc(1,0)*xc(1,1)*xc(2,0) - xc(0,1)*(2*pow(xc(1,0),2) - 3*xc(1,0)*xc(2,0) +
  xc(1,1)*(xc(1,1) - 4*xc(2,1))) + xc(0,0)*(xc(0,1)*(xc(1,0) - 3*xc(2,0)) +
  3*xc(1,1)*xc(2,0) + xc(1,0)*(xc(1,1) - 2*xc(2,1))) + pow(xc(1,0),2)*xc(2,1) -
  2*pow(xc(1,1),2)*xc(2,1) + pow(xc(0,0),2)*(xc(0,1) - 2*xc(1,1) + xc(2,1)) -
  pow(xc(0,1),2)*(xc(1,1) + 2*xc(2,1))))/normpowfive;

  elematrix(3,2)=
  (3*normal(0,0)*normal(1,0)*(normal(0,0)*(xc(0,0) - xc(2,0)) + normal(1,0)*(xc(0,1) -
  xc(2,1))) + normsquare*(normal(1,0)*(-xc(0,0) + xc(2,0)) + normal(0,0)*(-xc(0,1) +
  xc(2,1))))/normpowfive;

  elematrix(3,3)=
  -((normal(1,0)*(2*pow(xc(0,1),3) - 3*xc(1,0)*xc(1,1)*xc(2,0) +
  pow(xc(1,0),2)*xc(2,1) - 2*pow(xc(1,1),2)*xc(2,1) +
  pow(xc(0,0),2)*(2*xc(0,1) - 3*xc(1,1) + xc(2,1)) -
  2*pow(xc(0,1),2)*(2*xc(1,1) + xc(2,1)) - xc(0,0)*(-3*xc(1,1)*xc(2,0) +
  xc(0,1)*(xc(1,0) + 3*xc(2,0)) + xc(1,0)*(-3*xc(1,1) + 2*xc(2,1))) +
  xc(0,1)*(-pow(xc(1,0),2) + 3*xc(1,0)*xc(2,0) + 2*xc(1,1)*(xc(1,1) +
  2*xc(2,1)))))/normpowfive);

  elematrix(3,4)=-(pow(normal(1,0),2)/normcube);

  elematrix(3,5)=(normal(0,0)*(-xc(0,0) + xc(1,0)))/normcube;

  elematrix(4,0)=(normal(0,0)*normal(1,0))/normcube;

  elematrix(4,1)=pow(normal(1,0),2)/normcube;

  elematrix(4,2)=-((normal(0,0)*normal(1,0))/normcube);

  elematrix(4,3)=-(pow(normal(1,0),2)/normcube);

  elematrix(4,4)=0;

  elematrix(4,5)=0;

  elematrix(5,0)=-(pow(normal(0,0),2)/normcube);

  elematrix(5,1)=-((normal(0,0)*normal(1,0))/normcube);

  elematrix(5,2)=pow(normal(0,0),2)/normcube;

  elematrix(5,3)=(normal(0,0)*(-xc(0,0) + xc(1,0)))/normcube;

  elematrix(5,4)=0;

  elematrix(5,5)=0;
  return;

  elematrix.Scale(-1.0);

}

/*----------------------------------------------------------------------*
 * second derivatives */
void DRT::ELEMENTS::ConstraintElement::ComputeSecondDerivAngle2D
(
  const LINALG::Matrix<3,2>& xc,
  Epetra_SerialDenseMatrix& elematrix
)
{
  LINALG::SerialDenseVector vec1(2);
  vec1[1]=xc(0,0) - xc(1,0);
  vec1[0]=-(xc(0,1) - xc(1,1));

  LINALG::SerialDenseVector vec2(2);
  vec2[0]=-xc(1,0) + xc(2,0);
  vec2[1]=-xc(1,1) + xc(2,1);

  const double vec1sq=pow(vec1.Norm2(),2);
  const double vec2sq=pow(vec2.Norm2(),2);

  elematrix(0,0)
  =
  -(((-2*vec2sq*vec1[1]*vec2[1])/pow(vec1sq*vec2sq,1.5)
  - (vec2sq*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5) +
  (3*pow(vec2sq,2)*pow(vec1[1],2)*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) +
  xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,2.5))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) +
  ((vec2[1]/sqrt(vec1sq*vec2sq) -
  (vec2sq*vec1[1]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))*((-2*vec2[1]*(vec2(1
  )*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) +
  (2*vec1[1]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(pow(vec1sq,2)*vec2sq)))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(0,1)
  =
  -(((vec2sq*vec1[1]*vec2[0])/pow(vec1sq*vec2sq,1.5) +
  (vec2sq*vec1[0]*vec2[1])/pow(vec1sq*vec2sq,1.5) -
  (3*pow(vec2sq,2)*vec1[0]*vec1[1]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) +
  xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,2.5))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) +
  ((vec2[1]/sqrt(vec1sq*vec2sq) -
  (vec2sq*vec1[1]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))*((2*vec2[0]*(vec2[1]
  *xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) -
  (2*vec1[0]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(pow(vec1sq,2)*vec2sq)))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(0,2)
  =
  -((-((-2*vec2sq*vec1[1] -
  2*vec1sq*vec2[0])*vec2[1])/(2.*pow(vec1sq*vec2sq,1.5))
  - (vec2sq*vec1[1]*(xc(0,1) -
  xc(2,1)))/pow(vec1sq*vec2sq,1.5) +
  (vec2sq*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5) +
  (2*vec1[1]*vec2[0]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5) +
  (3*vec2sq*vec1[1]*(-2*vec2sq*vec1[1] -
  2*vec1sq*vec2[0])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,2.5)))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) +
  ((vec2[1]/sqrt(vec1sq*vec2sq) -
  (vec2sq*vec1[1]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))*((-2*(xc(0,1) -
  xc(2,1))*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) -
  (2*vec1[1]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(pow(vec1sq,2)*vec2sq) -
  (2*vec2[0]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(vec1sq*pow(vec2sq,2))))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(0,3)
  =
  -((-(1/sqrt(vec1sq*vec2sq)) - (vec2[1]*(2*vec2sq*vec1[0]
  - 2*vec1sq*vec2[1]))/(2.*pow(vec1sq*vec2sq,1.5)) -
  (vec2sq*vec1[1]*(-xc(0,0) +
  xc(2,0)))/pow(vec1sq*vec2sq,1.5) +
  (2*vec1[1]*vec2(1)*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5) +
  (3*vec2sq*vec1[1]*(2*vec2sq*vec1[0] -
  2*vec1sq*vec2[1])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,2.5)))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) +
  ((vec2[1]/sqrt(vec1sq*vec2sq) -
  (vec2sq*vec1[1]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))*((-2*(-xc(0,0) +
  xc(2,0))*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) +
  (2*vec1[0]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(pow(vec1sq,2)*vec2sq) -
  (2*vec2[1]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(vec1sq*pow(vec2sq,2))))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(0,4)
  =
  -((-((vec1sq*vec2[0]*vec2[1])/pow(vec1sq*vec2sq,1.5))
  - (vec2sq*vec1[1]*(-xc(0,1) +
  xc(1,1)))/pow(vec1sq*vec2sq,1.5) +
  (3*vec1sq*vec2sq*vec1[1]*vec2[0]*(vec2[1]*xc(0,0) -
  vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,2.5) -
  (2*vec1[1]*vec2[0]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) +
  ((vec2[1]/sqrt(vec1sq*vec2sq) -
  (vec2sq*vec1[1]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))*((-2*(-xc(0,1) +
  xc(1,1))*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) +
  (2*vec2[0]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(vec1sq*pow(vec2sq,2))))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(0,5)
  =
  -((1/sqrt(vec1sq*vec2sq) -
  (vec2sq*pow(vec1[1],2))/pow(vec1sq*vec2sq,1.5) -
  (vec1sq*pow(vec2[1],2))/pow(vec1sq*vec2sq,1.5) +
  (3*vec1sq*vec2sq*vec1[1]*vec2(1)*(vec2[1]*xc(0,0) -
  vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,2.5) -
  (2*vec1[1]*vec2(1)*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) +
  ((vec2[1]/sqrt(vec1sq*vec2sq) -
  (vec2sq*vec1[1]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))*((-2*vec1[1]*(vec2(1
  )*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) +
  (2*vec2[1]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(vec1sq*pow(vec2sq,2))))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(1,0)
  =
  -(((vec2sq*vec1[1]*vec2[0])/pow(vec1sq*vec2sq,1.5) +
  (vec2sq*vec1[0]*vec2[1])/pow(vec1sq*vec2sq,1.5) -
  (3*pow(vec2sq,2)*vec1[0]*vec1[1]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) +
  xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,2.5))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) +
  ((-(vec2[0]/sqrt(vec1sq*vec2sq)) +
  (vec2sq*vec1[0]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))*((-2*vec2[1]*(vec2(1
  )*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) +
  (2*vec1[1]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(pow(vec1sq,2)*vec2sq)))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(1,1)
  =
  -(((-2*vec2sq*vec1[0]*vec2[0])/pow(vec1sq*vec2sq,1.5)
  - (vec2sq*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5) +
  (3*pow(vec2sq,2)*pow(vec1[0],2)*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) +
  xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,2.5))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) +
  ((-(vec2[0]/sqrt(vec1sq*vec2sq)) +
  (vec2sq*vec1[0]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))*((2*vec2[0]*(vec2[1]
  *xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) -
  (2*vec1[0]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(pow(vec1sq,2)*vec2sq)))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(1,2)
  =
  -((1/sqrt(vec1sq*vec2sq) + (vec2[0]*(-2*vec2sq*vec1[1] -
  2*vec1sq*vec2[0]))/(2.*pow(vec1sq*vec2sq,1.5)) +
  (vec2sq*vec1[0]*(xc(0,1) -
  xc(2,1)))/pow(vec1sq*vec2sq,1.5) -
  (2*vec1[0]*vec2[0]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5) -
  (3*vec2sq*vec1[0]*(-2*vec2sq*vec1[1] -
  2*vec1sq*vec2[0])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,2.5)))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) +
  ((-(vec2[0]/sqrt(vec1sq*vec2sq)) +
  (vec2sq*vec1[0]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))*((-2*(xc(0,1) -
  xc(2,1))*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) -
  (2*vec1[1]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(pow(vec1sq,2)*vec2sq) -
  (2*vec2[0]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(vec1sq*pow(vec2sq,2))))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(1,3)
  =
  -(((vec2[0]*(2*vec2sq*vec1[0] -
  2*vec1sq*vec2[1]))/(2.*pow(vec1sq*vec2sq,1.5)) +
  (vec2sq*vec1[0]*(-xc(0,0) +
  xc(2,0)))/pow(vec1sq*vec2sq,1.5) +
  (vec2sq*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5) -
  (2*vec1[0]*vec2(1)*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5) -
  (3*vec2sq*vec1[0]*(2*vec2sq*vec1[0] -
  2*vec1sq*vec2[1])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,2.5)))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) +
  ((-(vec2[0]/sqrt(vec1sq*vec2sq)) +
  (vec2sq*vec1[0]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))*((-2*(-xc(0,0) +
  xc(2,0))*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) +
  (2*vec1[0]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(pow(vec1sq,2)*vec2sq) -
  (2*vec2[1]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(vec1sq*pow(vec2sq,2))))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(1,4)
  =
  -((-(1/sqrt(vec1sq*vec2sq)) +
  (vec1sq*pow(vec2[0],2))/pow(vec1sq*vec2sq,1.5) +
  (vec2sq*vec1[0]*(-xc(0,1) +
  xc(1,1)))/pow(vec1sq*vec2sq,1.5) -
  (3*vec1sq*vec2sq*vec1[0]*vec2[0]*(vec2[1]*xc(0,0) -
  vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,2.5) +
  (2*vec1[0]*vec2[0]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) +
  ((-(vec2[0]/sqrt(vec1sq*vec2sq)) +
  (vec2sq*vec1[0]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))*((-2*(-xc(0,1) +
  xc(1,1))*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) +
  (2*vec2[0]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(vec1sq*pow(vec2sq,2))))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(1,5)
  =
  -(((vec2sq*vec1[0]*vec1[1])/pow(vec1sq*vec2sq,1.5) +
  (vec1sq*vec2[0]*vec2[1])/pow(vec1sq*vec2sq,1.5) -
  (3*vec1sq*vec2sq*vec1[0]*vec2(1)*(vec2[1]*xc(0,0) -
  vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,2.5) +
  (2*vec1[0]*vec2(1)*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) +
  ((-(vec2[0]/sqrt(vec1sq*vec2sq)) +
  (vec2sq*vec1[0]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))*((-2*vec1[1]*(vec2(1
  )*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) +
  (2*vec2[1]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(vec1sq*pow(vec2sq,2))))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(2,0)
  =
  -((-((-2*vec2sq*vec1[1] -
  2*vec1sq*vec2[0])*vec2[1])/(2.*pow(vec1sq*vec2sq,1.5))
  - (vec2sq*vec1[1]*(xc(0,1) -
  xc(2,1)))/pow(vec1sq*vec2sq,1.5) +
  (3*vec2sq*vec1[1]*(-2*vec2sq*vec1[1] -
  2*vec1sq*vec2[0])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,2.5)) -
  ((-2*vec2sq - 4*vec1[1]*vec2[0])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) +
  xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,1.5)))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) + (((xc(0,1) -
  xc(2,1))/sqrt(vec1sq*vec2sq) - ((-2*vec2sq*vec1[1] -
  2*vec1sq*vec2[0])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,1.5)))*((-2*vec2[1]*(
  vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) +
  (2*vec1[1]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(pow(vec1sq,2)*vec2sq)))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(2,1)
  =
  -((1/sqrt(vec1sq*vec2sq) + (vec2[0]*(-2*vec2sq*vec1[1] -
  2*vec1sq*vec2[0]))/(2.*pow(vec1sq*vec2sq,1.5)) +
  (vec2sq*vec1[0]*(xc(0,1) -
  xc(2,1)))/pow(vec1sq*vec2sq,1.5) -
  (2*vec1[0]*vec2[0]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5) -
  (3*vec2sq*vec1[0]*(-2*vec2sq*vec1[1] -
  2*vec1sq*vec2[0])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,2.5)))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) + (((xc(0,1) -
  xc(2,1))/sqrt(vec1sq*vec2sq) - ((-2*vec2sq*vec1[1] -
  2*vec1sq*vec2[0])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,1.5)))*((2*vec2[0]*(
  vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) -
  (2*vec1[0]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(pow(vec1sq,2)*vec2sq)))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(2,2)
  =
  -((-(((-2*vec2sq*vec1[1] - 2*vec1sq*vec2[0])*(xc(0,1) -
  xc(2,1)))/pow(vec1sq*vec2sq,1.5)) +
  (3*pow(-2*vec2sq*vec1[1] -
  2*vec1sq*vec2[0],2)*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0)
  - xc(1,0)*xc(2,1)))/(4.*pow(vec1sq*vec2sq,2.5)) -
  ((2*vec1sq + 2*vec2sq + 8*vec1[1]*vec2[0])*(vec2[1]*xc(0,0) -
  vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,1.5)))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) + (((xc(0,1) -
  xc(2,1))/sqrt(vec1sq*vec2sq) - ((-2*vec2sq*vec1[1] -
  2*vec1sq*vec2[0])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,1.5)))*((-2*(xc(0,1) -
  xc(2,1))*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) -
  (2*vec1[1]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(pow(vec1sq,2)*vec2sq) -
  (2*vec2[0]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(vec1sq*pow(vec2sq,2))))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(2,3)
  =
  -((-((-2*vec2sq*vec1[1] - 2*vec1sq*vec2[0])*(-xc(0,0) +
  xc(2,0)))/(2.*pow(vec1sq*vec2sq,1.5)) -
  ((2*vec2sq*vec1[0] - 2*vec1sq*vec2[1])*(xc(0,1) -
  xc(2,1)))/(2.*pow(vec1sq*vec2sq,1.5)) +
  (3*(-2*vec2sq*vec1[1] -
  2*vec1sq*vec2[0])*(2*vec2sq*vec1[0] -
  2*vec1sq*vec2[1])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(4.*pow(vec1sq*vec2sq,2.5)) -
  ((-4*vec1[0]*vec2[0] + 4*vec1[1]*vec2[1])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) +
  xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,1.5)))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) + (((xc(0,1) -
  xc(2,1))/sqrt(vec1sq*vec2sq) - ((-2*vec2sq*vec1[1] -
  2*vec1sq*vec2[0])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,1.5)))*((-2*(-xc(0,0)
  + xc(2,0))*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) +
  (2*vec1[0]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(pow(vec1sq,2)*vec2sq) -
  (2*vec2[1]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(vec1sq*pow(vec2sq,2))))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(2,4)
  =
  -((-((-2*vec2sq*vec1[1] - 2*vec1sq*vec2[0])*(-xc(0,1) +
  xc(1,1)))/(2.*pow(vec1sq*vec2sq,1.5)) -
  (vec1sq*vec2[0]*(xc(0,1) -
  xc(2,1)))/pow(vec1sq*vec2sq,1.5) +
  (3*vec1sq*vec2[0]*(-2*vec2sq*vec1[1] -
  2*vec1sq*vec2[0])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,2.5)) -
  ((-2*vec1sq - 4*vec1[1]*vec2[0])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) +
  xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,1.5)))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) + (((xc(0,1) -
  xc(2,1))/sqrt(vec1sq*vec2sq) - ((-2*vec2sq*vec1[1] -
  2*vec1sq*vec2[0])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,1.5)))*((-2*(-xc(0,1)
  + xc(1,1))*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) +
  (2*vec2[0]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(vec1sq*pow(vec2sq,2))))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(2,5)
  =
  -((-(1/sqrt(vec1sq*vec2sq)) -
  (vec1[1]*(-2*vec2sq*vec1[1] -
  2*vec1sq*vec2[0]))/(2.*pow(vec1sq*vec2sq,1.5)) -
  (vec1sq*vec2[1]*(xc(0,1) -
  xc(2,1)))/pow(vec1sq*vec2sq,1.5) +
  (2*vec1[1]*vec2(1)*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5) +
  (3*vec1sq*(-2*vec2sq*vec1[1] -
  2*vec1sq*vec2[0])*vec2(1)*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) +
  xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,2.5)))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) + (((xc(0,1) -
  xc(2,1))/sqrt(vec1sq*vec2sq) - ((-2*vec2sq*vec1[1] -
  2*vec1sq*vec2[0])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,1.5)))*((-2*vec1[1]*(
  vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) +
  (2*vec2[1]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(vec1sq*pow(vec2sq,2))))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(3,0)
  =
  -((-(1/sqrt(vec1sq*vec2sq)) - (vec2[1]*(2*vec2sq*vec1[0]
  - 2*vec1sq*vec2[1]))/(2.*pow(vec1sq*vec2sq,1.5)) -
  (vec2sq*vec1[1]*(-xc(0,0) +
  xc(2,0)))/pow(vec1sq*vec2sq,1.5) +
  (2*vec1[1]*vec2(1)*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5) +
  (3*vec2sq*vec1[1]*(2*vec2sq*vec1[0] -
  2*vec1sq*vec2[1])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,2.5)))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) + (((-xc(0,0) +
  xc(2,0))/sqrt(vec1sq*vec2sq) - ((2*vec2sq*vec1[0] -
  2*vec1sq*vec2[1])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,1.5)))*((-2*vec2[1]*(
  vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) +
  (2*vec1[1]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(pow(vec1sq,2)*vec2sq)))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(3,1)
  =
  -(((vec2[0]*(2*vec2sq*vec1[0] -
  2*vec1sq*vec2[1]))/(2.*pow(vec1sq*vec2sq,1.5)) +
  (vec2sq*vec1[0]*(-xc(0,0) +
  xc(2,0)))/pow(vec1sq*vec2sq,1.5) -
  (3*vec2sq*vec1[0]*(2*vec2sq*vec1[0] -
  2*vec1sq*vec2[1])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,2.5)) -
  ((-2*vec2sq + 4*vec1[0]*vec2[1])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) +
  xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,1.5)))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) + (((-xc(0,0) +
  xc(2,0))/sqrt(vec1sq*vec2sq) - ((2*vec2sq*vec1[0] -
  2*vec1sq*vec2[1])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,1.5)))*((2*vec2[0]*(
  vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) -
  (2*vec1[0]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(pow(vec1sq,2)*vec2sq)))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(3,2)
  =
  -((-((-2*vec2sq*vec1[1] - 2*vec1sq*vec2[0])*(-xc(0,0) +
  xc(2,0)))/(2.*pow(vec1sq*vec2sq,1.5)) -
  ((2*vec2sq*vec1[0] - 2*vec1sq*vec2[1])*(xc(0,1) -
  xc(2,1)))/(2.*pow(vec1sq*vec2sq,1.5)) +
  (3*(-2*vec2sq*vec1[1] -
  2*vec1sq*vec2[0])*(2*vec2sq*vec1[0] -
  2*vec1sq*vec2[1])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(4.*pow(vec1sq*vec2sq,2.5)) -
  ((-4*vec1[0]*vec2[0] + 4*vec1[1]*vec2[1])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) +
  xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,1.5)))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) + (((-xc(0,0) +
  xc(2,0))/sqrt(vec1sq*vec2sq) - ((2*vec2sq*vec1[0] -
  2*vec1sq*vec2[1])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,1.5)))*((-2*(xc(0,1) -
  xc(2,1))*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) -
  (2*vec1[1]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(pow(vec1sq,2)*vec2sq) -
  (2*vec2[0]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(vec1sq*pow(vec2sq,2))))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(3,3)
  =
  -((-(((2*vec2sq*vec1[0] - 2*vec1sq*vec2[1])*(-xc(0,0) +
  xc(2,0)))/pow(vec1sq*vec2sq,1.5)) +
  (3*pow(2*vec2sq*vec1[0] - 2*vec1sq*vec2[1],2)*(vec2[1]*xc(0,0)
  - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(4.*pow(vec1sq*vec2sq,2.5)) -
  ((2*vec1sq + 2*vec2sq - 8*vec1[0]*vec2[1])*(vec2[1]*xc(0,0) -
  vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,1.5)))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) + (((-xc(0,0) +
  xc(2,0))/sqrt(vec1sq*vec2sq) - ((2*vec2sq*vec1[0] -
  2*vec1sq*vec2[1])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,1.5)))*((-2*(-xc(0,0)
  + xc(2,0))*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) +
  (2*vec1[0]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(pow(vec1sq,2)*vec2sq) -
  (2*vec2[1]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(vec1sq*pow(vec2sq,2))))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(3,4)
  =
  -((1/sqrt(vec1sq*vec2sq) - ((2*vec2sq*vec1[0] -
  2*vec1sq*vec2[1])*(-xc(0,1) +
  xc(1,1)))/(2.*pow(vec1sq*vec2sq,1.5)) -
  (vec1sq*vec2[0]*(-xc(0,0) +
  xc(2,0)))/pow(vec1sq*vec2sq,1.5) -
  (2*vec1[0]*vec2[0]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5) +
  (3*vec1sq*vec2[0]*(2*vec2sq*vec1[0] -
  2*vec1sq*vec2[1])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,2.5)))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) + (((-xc(0,0) +
  xc(2,0))/sqrt(vec1sq*vec2sq) - ((2*vec2sq*vec1[0] -
  2*vec1sq*vec2[1])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,1.5)))*((-2*(-xc(0,1)
  + xc(1,1))*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) +
  (2*vec2[0]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(vec1sq*pow(vec2sq,2))))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(3,5)
  =
  -((-(vec1[1]*(2*vec2sq*vec1[0] -
  2*vec1sq*vec2[1]))/(2.*pow(vec1sq*vec2sq,1.5)) -
  (vec1sq*vec2[1]*(-xc(0,0) +
  xc(2,0)))/pow(vec1sq*vec2sq,1.5) +
  (3*vec1sq*vec2[1]*(2*vec2sq*vec1[0] -
  2*vec1sq*vec2[1])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,2.5)) -
  ((-2*vec1sq + 4*vec1[0]*vec2[1])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) +
  xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,1.5)))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) + (((-xc(0,0) +
  xc(2,0))/sqrt(vec1sq*vec2sq) - ((2*vec2sq*vec1[0] -
  2*vec1sq*vec2[1])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,1.5)))*((-2*vec1[1]*(
  vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) +
  (2*vec2[1]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(vec1sq*pow(vec2sq,2))))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(4,0)
  =
  -((-((vec1sq*vec2[0]*vec2[1])/pow(vec1sq*vec2sq,1.5))
  - (vec2sq*vec1[1]*(-xc(0,1) +
  xc(1,1)))/pow(vec1sq*vec2sq,1.5) +
  (3*vec1sq*vec2sq*vec1[1]*vec2[0]*(vec2[1]*xc(0,0) -
  vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,2.5) -
  (2*vec1[1]*vec2[0]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) + (((-xc(0,1) +
  xc(1,1))/sqrt(vec1sq*vec2sq) -
  (vec1sq*vec2[0]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))*((-2*vec2[1]*(vec2(1
  )*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) +
  (2*vec1[1]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(pow(vec1sq,2)*vec2sq)))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(4,1)
  =
  -((-(1/sqrt(vec1sq*vec2sq)) +
  (vec1sq*pow(vec2[0],2))/pow(vec1sq*vec2sq,1.5) +
  (vec2sq*vec1[0]*(-xc(0,1) +
  xc(1,1)))/pow(vec1sq*vec2sq,1.5) -
  (3*vec1sq*vec2sq*vec1[0]*vec2[0]*(vec2[1]*xc(0,0) -
  vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,2.5) +
  (2*vec1[0]*vec2[0]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) + (((-xc(0,1) +
  xc(1,1))/sqrt(vec1sq*vec2sq) -
  (vec1sq*vec2[0]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))*((2*vec2[0]*(vec2[1]
  *xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) -
  (2*vec1[0]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(pow(vec1sq,2)*vec2sq)))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(4,2)
  =
  -((-((-2*vec2sq*vec1[1] - 2*vec1sq*vec2[0])*(-xc(0,1) +
  xc(1,1)))/(2.*pow(vec1sq*vec2sq,1.5)) -
  (vec1sq*vec2[0]*(xc(0,1) -
  xc(2,1)))/pow(vec1sq*vec2sq,1.5) +
  (vec1sq*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5) +
  (2*vec1[1]*vec2[0]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5) +
  (3*vec1sq*vec2[0]*(-2*vec2sq*vec1[1] -
  2*vec1sq*vec2[0])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,2.5)))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) + (((-xc(0,1) +
  xc(1,1))/sqrt(vec1sq*vec2sq) -
  (vec1sq*vec2[0]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))*((-2*(xc(0,1) -
  xc(2,1))*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) -
  (2*vec1[1]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(pow(vec1sq,2)*vec2sq) -
  (2*vec2[0]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(vec1sq*pow(vec2sq,2))))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(4,3)
  =
  -((1/sqrt(vec1sq*vec2sq) - ((2*vec2sq*vec1[0] -
  2*vec1sq*vec2[1])*(-xc(0,1) +
  xc(1,1)))/(2.*pow(vec1sq*vec2sq,1.5)) -
  (vec1sq*vec2[0]*(-xc(0,0) +
  xc(2,0)))/pow(vec1sq*vec2sq,1.5) -
  (2*vec1[0]*vec2[0]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5) +
  (3*vec1sq*vec2[0]*(2*vec2sq*vec1[0] -
  2*vec1sq*vec2[1])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,2.5)))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) + (((-xc(0,1) +
  xc(1,1))/sqrt(vec1sq*vec2sq) -
  (vec1sq*vec2[0]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))*((-2*(-xc(0,0) +
  xc(2,0))*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) +
  (2*vec1[0]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(pow(vec1sq,2)*vec2sq) -
  (2*vec2[1]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(vec1sq*pow(vec2sq,2))))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(4,4)
  =
  -(((-2*vec1sq*vec2[0]*(-xc(0,1) +
  xc(1,1)))/pow(vec1sq*vec2sq,1.5) -
  (vec1sq*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5) +
  (3*pow(vec1sq,2)*pow(vec2[0],2)*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) +
  xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,2.5))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) + (((-xc(0,1) +
  xc(1,1))/sqrt(vec1sq*vec2sq) -
  (vec1sq*vec2[0]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))*((-2*(-xc(0,1) +
  xc(1,1))*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) +
  (2*vec2[0]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(vec1sq*pow(vec2sq,2))))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(4,5)
  =
  -((-((vec1sq*vec1[1]*vec2[0])/pow(vec1sq*vec2sq,1.5))
  - (vec1sq*vec2[1]*(-xc(0,1) +
  xc(1,1)))/pow(vec1sq*vec2sq,1.5) +
  (3*pow(vec1sq,2)*vec2[0]*vec2(1)*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) +
  xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,2.5))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) + (((-xc(0,1) +
  xc(1,1))/sqrt(vec1sq*vec2sq) -
  (vec1sq*vec2[0]*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))*((-2*vec1[1]*(vec2(1
  )*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) +
  (2*vec2[1]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(vec1sq*pow(vec2sq,2))))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(5,0)
  =
  -((1/sqrt(vec1sq*vec2sq) -
  (vec2sq*pow(vec1[1],2))/pow(vec1sq*vec2sq,1.5) -
  (vec1sq*pow(vec2[1],2))/pow(vec1sq*vec2sq,1.5) +
  (3*vec1sq*vec2sq*vec1[1]*vec2(1)*(vec2[1]*xc(0,0) -
  vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,2.5) -
  (2*vec1[1]*vec2(1)*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) +
  ((vec1[1]/sqrt(vec1sq*vec2sq) -
  (vec1sq*vec2(1)*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))*((-2*vec2[1]*(vec2(1
  )*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) +
  (2*vec1[1]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(pow(vec1sq,2)*vec2sq)))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(5,1)
  =
  -(((vec2sq*vec1[0]*vec1[1])/pow(vec1sq*vec2sq,1.5) +
  (vec1sq*vec2[0]*vec2[1])/pow(vec1sq*vec2sq,1.5) -
  (3*vec1sq*vec2sq*vec1[0]*vec2(1)*(vec2[1]*xc(0,0) -
  vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,2.5) +
  (2*vec1[0]*vec2(1)*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) +
  ((vec1[1]/sqrt(vec1sq*vec2sq) -
  (vec1sq*vec2(1)*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))*((2*vec2[0]*(vec2[1]
  *xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) -
  (2*vec1[0]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(pow(vec1sq,2)*vec2sq)))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(5,2)
  =
  -((-(1/sqrt(vec1sq*vec2sq)) -
  (vec1[1]*(-2*vec2sq*vec1[1] -
  2*vec1sq*vec2[0]))/(2.*pow(vec1sq*vec2sq,1.5)) -
  (vec1sq*vec2[1]*(xc(0,1) -
  xc(2,1)))/pow(vec1sq*vec2sq,1.5) +
  (2*vec1[1]*vec2(1)*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5) +
  (3*vec1sq*(-2*vec2sq*vec1[1] -
  2*vec1sq*vec2[0])*vec2(1)*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) +
  xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,2.5)))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) +
  ((vec1[1]/sqrt(vec1sq*vec2sq) -
  (vec1sq*vec2(1)*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))*((-2*(xc(0,1) -
  xc(2,1))*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) -
  (2*vec1[1]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(pow(vec1sq,2)*vec2sq) -
  (2*vec2[0]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(vec1sq*pow(vec2sq,2))))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(5,3)
  =
  -((-(vec1[1]*(2*vec2sq*vec1[0] -
  2*vec1sq*vec2[1]))/(2.*pow(vec1sq*vec2sq,1.5)) -
  (vec1sq*vec2[1]*(-xc(0,0) +
  xc(2,0)))/pow(vec1sq*vec2sq,1.5) +
  (vec1sq*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5) -
  (2*vec1[0]*vec2(1)*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5) +
  (3*vec1sq*vec2[1]*(2*vec2sq*vec1[0] -
  2*vec1sq*vec2[1])*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(2.*pow(vec1sq*vec2sq,2.5)))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) +
  ((vec1[1]/sqrt(vec1sq*vec2sq) -
  (vec1sq*vec2(1)*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))*((-2*(-xc(0,0) +
  xc(2,0))*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) +
  (2*vec1[0]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(pow(vec1sq,2)*vec2sq) -
  (2*vec2[1]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(vec1sq*pow(vec2sq,2))))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(5,4)
  =
  -((-((vec1sq*vec1[1]*vec2[0])/pow(vec1sq*vec2sq,1.5))
  - (vec1sq*vec2[1]*(-xc(0,1) +
  xc(1,1)))/pow(vec1sq*vec2sq,1.5) +
  (3*pow(vec1sq,2)*vec2[0]*vec2(1)*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) +
  xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,2.5))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) +
  ((vec1[1]/sqrt(vec1sq*vec2sq) -
  (vec1sq*vec2(1)*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))*((-2*(-xc(0,1) +
  xc(1,1))*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) +
  (2*vec2[0]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(vec1sq*pow(vec2sq,2))))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;
  elematrix(5,5)
  =
  -(((-2*vec1sq*vec1[1]*vec2[1])/pow(vec1sq*vec2sq,1.5)
  - (vec1sq*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5) +
  (3*pow(vec1sq,2)*pow(vec2[1],2)*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) +
  xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,2.5))/sqrt(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq))) +
  ((vec1[1]/sqrt(vec1sq*vec2sq) -
  (vec1sq*vec2(1)*(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/pow(vec1sq*vec2sq,1.5))*((-2*vec1[1]*(vec2(1
  )*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1)))/(vec1sq*vec2sq) +
  (2*vec2[1]*pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2))/(vec1sq*pow(vec2sq,2))))/(2.*pow(1 -
  pow(vec2[1]*xc(0,0) - vec2[0]*xc(0,1) + xc(1,1)*xc(2,0) -
  xc(1,0)*xc(2,1),2)/(vec1sq*vec2sq),1.5))
  ;

  elematrix.Scale(-1.0);
}

double DRT::ELEMENTS::ConstraintElement::ComputeWeightedDistance
(
  const vector<double> disp,
  const vector<double> direct
)
{

  // norm of direct
  double norm = sqrt(pow(direct.at(0),2)+pow(direct.at(1),2)+pow(direct.at(2),2));
  double result=0.0;
  
  for(int i = 0;i<3;i++)
  {
    result += (disp.at(i)-disp.at(i+3))*direct.at(i);
  }
  result/=norm;
  return result;
}

void DRT::ELEMENTS::ConstraintElement::ComputeFirstDerivWeightedDistance
 (
   Epetra_SerialDenseVector& elevector,
   const vector<double> direct  
 )
{
  // norm of direct
  double norm = sqrt(pow(direct.at(0),2)+pow(direct.at(1),2)+pow(direct.at(2),2));
  
  for(int i = 0;i<3;i++)
  {
    elevector(i)=-direct.at(i)/norm;
    elevector(3+i)=direct.at(i)/norm;
  }
  
  return;
}

//=======================================================================


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::ConstraintElementRegister::Initialize(DRT::Discretization& dis)
{
  return 0;
}

#endif  // #ifdef CCADISCRET
