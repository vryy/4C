/*!----------------------------------------------------------------------
\file constraint_element3_evaluate.cpp
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
#include "constraint_element3.H"
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
int DRT::ELEMENTS::ConstraintElement3::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  ActionType act = none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") return 0;
  else if (action=="calc_MPC_stiff")       act = calc_MPC_stiffness;
  else
    dserror("Unknown type of action for ConstraintElement3");

  switch (act)
  {
    case none:
    {
    }
    break;
    case calc_MPC_stiffness:
    {
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get state vector 'displacement'");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      const int numnod=NumNode();//numnod = 4;
      const int numdim=3;
      const int numdofpernode=3;
      LINALG::SerialDenseMatrix xsrefe(numnod,numdim);  // material coord. of element
      LINALG::SerialDenseMatrix xscurr(numnod,numdim);  // material coord. of element
      for (int i=0; i<numnod; ++i)
      {
        for (int j = 0; j < numdofpernode; ++j)
        {
          xsrefe(i,j) = Nodes()[i]->X()[j];
          xscurr(i,j) = xsrefe(i,j) + mydisp[i*numdofpernode+j];
        }
      }

      LINALG::SerialDenseVector elementnormal(numdim);
      ComputeElementNormal(xscurr,elementnormal);
      if(abs(elementnormal.Norm2())<1E-6)
      {
        dserror("Bad plane, points almost on a line!");
      }
      double normaldistance =ComputeNormalDist(xscurr,elementnormal);
      ComputeFirstDeriv(xscurr,elevec1,elementnormal);
      ComputeSecondDeriv(xscurr,elemat1,elementnormal);
      const int ID =params.get("ConditionID",-1);
      RCP<Epetra_Vector> lambdav=rcp(new Epetra_Vector(*(params.get<RCP<Epetra_Vector> >("LagrMultVector"))));

      if (ID<0)
      {
        dserror("Condition ID for volume constraint missing!");
      }

      const int minID =params.get("MinID",0);
      //update corresponding column in "constraint" matrix
      elevec2=elevec1;
      elevec1.Scale(1*(*lambdav)[ID-minID]);
      elemat1.Scale(1*(*lambdav)[ID-minID]);
      //call submethod for volume evaluation
      if(discretization.Comm().MyPID()==Owner())
      {
        // write normal distance to parameter list
        char ndistname[30];
        sprintf(ndistname,"computed normal distance %d",ID);
        //update volume in parameter list
        params.set(ndistname, normaldistance);
      }
    }

    break;
    default:
      dserror("Unimplemented type of action");
  }
  return 0;


} // end of DRT::ELEMENTS::ConstraintElement3::Evaluate




int DRT::ELEMENTS::ConstraintElement3::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{

  return 0;
}

void DRT::ELEMENTS::ConstraintElement3::ComputeElementNormal(const LINALG::SerialDenseMatrix& xc,
    LINALG::SerialDenseVector& elenorm)
{
  elenorm[0]=-(xc(0,2)*xc(1,1)) + xc(0,1)*xc(1,2) + xc(0,2)*xc(2,1) -
    xc(1,2)*xc(2,1) - xc(0,1)*xc(2,2) + xc(1,1)*xc(2,2);
  elenorm[1]=xc(0,2)*xc(1,0) - xc(0,0)*xc(1,2) - xc(0,2)*xc(2,0) +
    xc(1,2)*xc(2,0) + xc(0,0)*xc(2,2) - xc(1,0)*xc(2,2);
  elenorm[2]=-(xc(0,1)*xc(1,0)) + xc(0,0)*xc(1,1) + xc(0,1)*xc(2,0) -
    xc(1,1)*xc(2,0) - xc(0,0)*xc(2,1) + xc(1,0)*xc(2,1);
  return ;
}

double DRT::ELEMENTS::ConstraintElement3::ComputeNormalDist(const LINALG::SerialDenseMatrix& xc,
    const LINALG::SerialDenseVector& normal)
{
  return (-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*
      (xc(0,2) - xc(3,2)))/normal.Norm2();
}

void DRT::ELEMENTS::ConstraintElement3::ComputeFirstDeriv(const LINALG::SerialDenseMatrix& xc,
    Epetra_SerialDenseVector& elevector,
                                                          const LINALG::SerialDenseVector& normal)
{ 
  double normsquare=pow(normal.Norm2(),2);
  double normcube=pow(normal.Norm2(),3);
  
  elevector[0]=
    (-((-2*normal[3]*(-xc(1,1) + xc(2,1)) + 2*normal[2]*(-xc(1,2) +
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2)))) + 2*normsquare*(xc(1,2)*(xc(2,1) - xc(3,1)) +
    xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) + xc(1,1)*(-xc(2,2) +
    xc(3,2))))/(2.*normcube);

  elevector[1]=
    (-((-2*(normal[3])*(xc(1,0) - xc(2,0)) - 2*normal[1]*(-xc(1,2) +
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2)))) + 2*normsquare*(-(xc(2,2)*xc(3,0)) +
    xc(1,2)*(-xc(2,0) + xc(3,0)) + xc(1,0)*(xc(2,2) - xc(3,2)) +
    xc(2,0)*xc(3,2)))/(2.*normcube);

  elevector[2]=
    (2*normsquare*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) - (2*normal[2]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(xc(1,1) - xc(2,1)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) -
    xc(3,2))))/(2.*normcube);

  elevector[3]=
    (-((-2*normal[3]*(xc(0,1) - xc(2,1)) + 2*normal[2]*(xc(0,2) -
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2)))) + 2*normsquare*((xc(0,2) - xc(2,2))*(-xc(0,1) +
    xc(3,1)) + (xc(0,1) - xc(2,1))*(xc(0,2) - xc(3,2))))/(2.*normcube);

  elevector[4]=
    (-((-2*normal[3]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(xc(0,2) -
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2)))) + 2*normsquare*((xc(0,2) - xc(2,2))*(xc(0,0) -
    xc(3,0)) + (-xc(0,0) + xc(2,0))*(xc(0,2) - xc(3,2))))/(2.*normcube);

  elevector[5]=
    (2*normsquare*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) - (2*normal[2]*(-xc(0,0) + xc(2,0)) -
    2*normal[1]*(-xc(0,1) + xc(2,1)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) -
    xc(3,2))))/(2.*normcube);

  elevector[6]=
    (-((-2*normal[3]*(-xc(0,1) + xc(1,1)) + 2*normal[2]*(-xc(0,2) +
    xc(1,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2)))) + 2*normsquare*((xc(0,2) - xc(1,2))*(xc(0,1) -
    xc(3,1)) + (-xc(0,1) + xc(1,1))*(xc(0,2) - xc(3,2))))/(2.*normcube);

  elevector[7]=
    (-((-2*normal[3]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(-xc(0,2) +
    xc(1,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2)))) + 2*normsquare*((-xc(0,2) + xc(1,2))*(xc(0,0) -
    xc(3,0)) + (xc(0,0) - xc(1,0))*(xc(0,2) - xc(3,2))))/(2.*normcube);

  elevector[8]=
    (2*normsquare*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(-xc(0,1) + xc(3,1))) - (2*normal[2]*(xc(0,0) - xc(1,0)) -
    2*normal[1]*(xc(0,1) - xc(1,1)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) -
    xc(3,2))))/(2.*normcube);

  elevector[9]=
    (-(xc(1,2)*xc(2,1)) + xc(0,2)*(-xc(1,1) + xc(2,1)) + xc(0,1)*(xc(1,2) - xc(2,2))
    + xc(1,1)*xc(2,2))/normal.Norm2();

  elevector[10]=
    normal[2]/normal.Norm2();

  elevector[11]=
    (-(xc(1,1)*xc(2,0)) + xc(0,1)*(-xc(1,0) + xc(2,0)) + xc(0,0)*(xc(1,1) - xc(2,1))
    + xc(1,0)*xc(2,1))/normal.Norm2();
  
  return;
}

void DRT::ELEMENTS::ConstraintElement3::ComputeSecondDeriv(const LINALG::SerialDenseMatrix& xc,
    Epetra_SerialDenseMatrix& elematrix,
    const LINALG::SerialDenseVector& normal)
{
  
  double normsquare=pow(normal.Norm2(),2);
  double normcube=pow(normal.Norm2(),3);
  double normpowfour=pow(normal.Norm2(),4);
  double normpowfive=pow(normal.Norm2(),5);
  
  elematrix(0,0)=
    (-4*normsquare*(pow(xc(1,1) - xc(2,1),2) + pow(xc(1,2) -
    xc(2,2),2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) + 3*pow(-2*normal[3]*(-xc(1,1) + xc(2,1)) +
    2*normal[2]*(-xc(1,2) + xc(2,2)),2)*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-2*normal[3]*(-xc(1,1) + xc(2,1)) + 2*normal[2]*(-xc(1,2) +
    xc(2,2)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) +
    xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(0,1)=
    (-4*normsquare*(xc(1,0) - xc(2,0))*(-xc(1,1) + xc(2,1))*(-(normal[1]*(xc(0,0) -
    xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) +
    3*(-2*normal[3]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(-xc(1,2) +
    xc(2,2)))*(-2*normal[3]*(-xc(1,1) + xc(2,1)) + 2*normal[2]*(-xc(1,2) +
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(-xc(1,1) + xc(2,1))
    + 2*normal[2]*(-xc(1,2) + xc(2,2)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) +
    xc(3,0)) + xc(1,0)*(xc(2,2) - xc(3,2)) + xc(2,0)*xc(3,2)) -
    2*normsquare*(-2*normal[3]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(-xc(1,2) +
    xc(2,2)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) +
    xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(0,2)=
    (-2*normsquare*(-2*normal[3]*(-xc(1,1) + xc(2,1)) + 2*normal[2]*(-xc(1,2) +
    xc(2,2)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) - 4*normsquare*(xc(1,0) - xc(2,0))*(-xc(1,2) +
    xc(2,2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) + 3*(2*normal[2]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(xc(1,1) - xc(2,1)))*(-2*normal[3]*(-xc(1,1) + xc(2,1)) +
    2*normal[2]*(-xc(1,2) + xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal[2]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(xc(1,1) -
    xc(2,1)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) +
    xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(0,3)=
    (3*(-2*normal[3]*(xc(0,1) - xc(2,1)) + 2*normal[2]*(xc(0,2) -
    xc(2,2)))*(-2*normal[3]*(-xc(1,1) + xc(2,1)) + 2*normal[2]*(-xc(1,2) +
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*(xc(0,1) - xc(2,1))*(-xc(1,1) +
    xc(2,1)) + 2*(xc(0,2) - xc(2,2))*(-xc(1,2) + xc(2,2)))*(-(normal[1]*(xc(0,0) -
    xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal[3]*(-xc(1,1) + xc(2,1)) + 2*normal[2]*(-xc(1,2) +
    xc(2,2)))*((xc(0,2) - xc(2,2))*(-xc(0,1) + xc(3,1)) + (xc(0,1) -
    xc(2,1))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(xc(0,1) - xc(2,1)) +
    2*normal[2]*(xc(0,2) - xc(2,2)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1)
    - xc(2,1)*xc(3,2) + xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(0,4)=
    (-4*normsquare*(-2*xc(1,1)*xc(2,0) + xc(0,1)*(-xc(1,0) + xc(2,0)) +
    2*xc(0,0)*(xc(1,1) - xc(2,1)) + xc(1,0)*xc(2,1) +
    xc(2,0)*xc(2,1))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) +
    xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) + 3*(-2*normal[3]*(-xc(0,0) + xc(2,0))
    - 2*normal[1]*(xc(0,2) - xc(2,2)))*(-2*normal[3]*(-xc(1,1) + xc(2,1)) +
    2*normal[2]*(-xc(1,2) + xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal[3]*(-xc(1,1) + xc(2,1)) + 2*normal[2]*(-xc(1,2) +
    xc(2,2)))*((xc(0,2) - xc(2,2))*(xc(0,0) - xc(3,0)) + (-xc(0,0) +
    xc(2,0))*(xc(0,2) - xc(3,2))) + 4*normpowfour*(-xc(2,2) + xc(3,2)) -
    2*normsquare*(-2*normal[3]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(xc(0,2) -
    xc(2,2)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) +
    xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(0,5)=
    (-2*normsquare*(-2*normal[3]*(-xc(1,1) + xc(2,1)) + 2*normal[2]*(-xc(1,2) +
    xc(2,2)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) + 4*normpowfour*(xc(2,1) - xc(3,1)) -
    4*normsquare*(-2*xc(1,2)*xc(2,0) + xc(0,2)*(-xc(1,0) + xc(2,0)) +
    2*xc(0,0)*(xc(1,2) - xc(2,2)) + xc(1,0)*xc(2,2) +
    xc(2,0)*xc(2,2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) +
    xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) + 3*(2*normal[2]*(-xc(0,0) + xc(2,0))
    - 2*normal[1]*(-xc(0,1) + xc(2,1)))*(-2*normal[3]*(-xc(1,1) + xc(2,1)) +
    2*normal[2]*(-xc(1,2) + xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal[2]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(-xc(0,1) +
    xc(2,1)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) +
    xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(0,6)=
    (-2*normsquare*(2*(xc(0,1) - xc(1,1))*(xc(1,1) - xc(2,1)) + 2*(xc(0,2) -
    xc(1,2))*(xc(1,2) - xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) +
    3*(-2*normal[3]*(-xc(0,1) + xc(1,1)) + 2*normal[2]*(-xc(0,2) +
    xc(1,2)))*(-2*normal[3]*(-xc(1,1) + xc(2,1)) + 2*normal[2]*(-xc(1,2) +
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(-xc(1,1) + xc(2,1))
    + 2*normal[2]*(-xc(1,2) + xc(2,2)))*((xc(0,2) - xc(1,2))*(xc(0,1) - xc(3,1)) +
    (-xc(0,1) + xc(1,1))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(-xc(0,1)
    + xc(1,1)) + 2*normal[2]*(-xc(0,2) + xc(1,2)))*(xc(1,2)*(xc(2,1) - xc(3,1)) +
    xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) + xc(1,1)*(-xc(2,2) +
    xc(3,2))))/(4.*normpowfive);

  elematrix(0,7)=
    (-4*normsquare*(xc(1,0)*xc(1,1) + xc(0,1)*(xc(1,0) - xc(2,0)) + xc(1,1)*xc(2,0)
    - 2*xc(0,0)*(xc(1,1) - xc(2,1)) - 2*xc(1,0)*xc(2,1))*(-(normal[1]*(xc(0,0) -
    xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) +
    3*(-2*normal[3]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(-xc(0,2) +
    xc(1,2)))*(-2*normal[3]*(-xc(1,1) + xc(2,1)) + 2*normal[2]*(-xc(1,2) +
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(-xc(1,1) + xc(2,1))
    + 2*normal[2]*(-xc(1,2) + xc(2,2)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) +
    (xc(0,0) - xc(1,0))*(xc(0,2) - xc(3,2))) + 4*normpowfour*(xc(1,2) -
    xc(3,2)) - 2*normsquare*(-2*normal[3]*(xc(0,0) - xc(1,0)) -
    2*normal[1]*(-xc(0,2) + xc(1,2)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1)
    - xc(2,1)*xc(3,2) + xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(0,8)=
    (4*normpowfour*(-xc(1,1) + xc(3,1)) -
    2*normsquare*(-2*normal[3]*(-xc(1,1) + xc(2,1)) + 2*normal[2]*(-xc(1,2) +
    xc(2,2)))*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(-xc(0,1) + xc(3,1))) - 4*normsquare*(xc(1,0)*xc(1,2) +
    xc(0,2)*(xc(1,0) - xc(2,0)) + xc(1,2)*xc(2,0) - 2*xc(0,0)*(xc(1,2) - xc(2,2)) -
    2*xc(1,0)*xc(2,2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) +
    xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) + 3*(2*normal[2]*(xc(0,0) - xc(1,0)) -
    2*normal[1]*(xc(0,1) - xc(1,1)))*(-2*normal[3]*(-xc(1,1) + xc(2,1)) +
    2*normal[2]*(-xc(1,2) + xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal[2]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(xc(0,1) -
    xc(1,1)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) +
    xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(0,9)=
    -((-(xc(1,2)*xc(2,1)) + xc(0,2)*(-xc(1,1) + xc(2,1)) + xc(0,1)*(xc(1,2) -
    xc(2,2)) + xc(1,1)*xc(2,2))*(-2*normal[3]*(-xc(1,1) + xc(2,1)) +
    2*normal[2]*(-xc(1,2) + xc(2,2))))/(2.*normcube);

  elematrix(0,10)=
    (normal[1]*(xc(0,2)*xc(1,1)*xc(1,2) - xc(1,0)*xc(1,1)*xc(2,0) +
    xc(1,1)*pow(xc(2,0),2) + xc(0,0)*(xc(1,0) - xc(2,0))*(xc(1,1) - xc(2,1)) +
    pow(xc(1,0),2)*xc(2,1) - xc(0,2)*xc(1,2)*xc(2,1) + pow(xc(1,2),2)*xc(2,1) -
    xc(1,0)*xc(2,0)*xc(2,1) - xc(0,2)*xc(1,1)*xc(2,2) - xc(1,1)*xc(1,2)*xc(2,2) +
    xc(0,2)*xc(2,1)*xc(2,2) - xc(1,2)*xc(2,1)*xc(2,2) + xc(1,1)*pow(xc(2,2),2) -
    xc(0,1)*(pow(xc(1,0),2) + pow(xc(1,2),2) - 2*xc(1,0)*xc(2,0) +
    pow(xc(2,0),2) - 2*xc(1,2)*xc(2,2) + pow(xc(2,2),2))))/normcube;

  elematrix(0,11)=
    -((normal[1]*(-(xc(0,1)*xc(1,1)*xc(1,2)) + xc(1,0)*xc(1,2)*xc(2,0) -
    xc(1,2)*pow(xc(2,0),2) + xc(0,1)*xc(1,2)*xc(2,1) + xc(1,1)*xc(1,2)*xc(2,1) -
    xc(1,2)*pow(xc(2,1),2) + xc(0,2)*(pow(xc(1,0),2) + pow(xc(1,1),2) -
    2*xc(1,0)*xc(2,0) + pow(xc(2,0),2) - 2*xc(1,1)*xc(2,1) + pow(xc(2,1),2)) -
    xc(0,0)*(xc(1,0) - xc(2,0))*(xc(1,2) - xc(2,2)) - pow(xc(1,0),2)*xc(2,2) +
    xc(0,1)*xc(1,1)*xc(2,2) - pow(xc(1,1),2)*xc(2,2) + xc(1,0)*xc(2,0)*xc(2,2) -
    xc(0,1)*xc(2,1)*xc(2,2) + xc(1,1)*xc(2,1)*xc(2,2)))/normcube);

  elematrix(1,0)=
    (-4*normsquare*(xc(1,0) - xc(2,0))*(-xc(1,1) + xc(2,1))*(-(normal[1]*(xc(0,0) -
    xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) +
    3*(-2*normal[3]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(-xc(1,2) +
    xc(2,2)))*(-2*normal[3]*(-xc(1,1) + xc(2,1)) + 2*normal[2]*(-xc(1,2) +
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(-xc(1,1) + xc(2,1))
    + 2*normal[2]*(-xc(1,2) + xc(2,2)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) +
    xc(3,0)) + xc(1,0)*(xc(2,2) - xc(3,2)) + xc(2,0)*xc(3,2)) -
    2*normsquare*(-2*normal[3]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(-xc(1,2) +
    xc(2,2)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) +
    xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(1,1)=
    (-4*normsquare*(pow(xc(1,0) - xc(2,0),2) + pow(xc(1,2) -
    xc(2,2),2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) + 3*pow(-2*normal[3]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(-xc(1,2) + xc(2,2)),2)*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-2*normal[3]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(-xc(1,2) +
    xc(2,2)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) + xc(3,0)) + xc(1,0)*(xc(2,2)
    - xc(3,2)) + xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(1,2)=
    (-2*normsquare*(-2*normal[3]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(-xc(1,2) +
    xc(2,2)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) - 4*normsquare*(xc(1,1) - xc(2,1))*(-xc(1,2) +
    xc(2,2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) + 3*(2*normal[2]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(xc(1,1) - xc(2,1)))*(-2*normal[3]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(-xc(1,2) + xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal[2]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(xc(1,1) -
    xc(2,1)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) + xc(3,0)) + xc(1,0)*(xc(2,2)
    - xc(3,2)) + xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(1,3)=
    (-4*normsquare*(2*xc(0,1)*(xc(1,0) - xc(2,0)) + xc(1,1)*xc(2,0) -
    2*xc(1,0)*xc(2,1) + xc(2,0)*xc(2,1) + xc(0,0)*(-xc(1,1) +
    xc(2,1)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) + 3*(-2*normal[3]*(xc(0,1) - xc(2,1)) +
    2*normal[2]*(xc(0,2) - xc(2,2)))*(-2*normal[3]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(-xc(1,2) + xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal[3]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(-xc(1,2) +
    xc(2,2)))*((xc(0,2) - xc(2,2))*(-xc(0,1) + xc(3,1)) + (xc(0,1) -
    xc(2,1))*(xc(0,2) - xc(3,2))) + 4*normpowfour*(xc(2,2) - xc(3,2)) -
    2*normsquare*(-2*normal[3]*(xc(0,1) - xc(2,1)) + 2*normal[2]*(xc(0,2) -
    xc(2,2)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) + xc(3,0)) + xc(1,0)*(xc(2,2)
    - xc(3,2)) + xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(1,4)=
    (3*(-2*normal[3]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(xc(0,2) -
    xc(2,2)))*(-2*normal[3]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(-xc(1,2) +
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*(xc(1,0) - xc(2,0))*(-xc(0,0) +
    xc(2,0)) + 2*(xc(0,2) - xc(2,2))*(-xc(1,2) + xc(2,2)))*(-(normal[1]*(xc(0,0) -
    xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal[3]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(-xc(1,2) +
    xc(2,2)))*((xc(0,2) - xc(2,2))*(xc(0,0) - xc(3,0)) + (-xc(0,0) +
    xc(2,0))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(-xc(0,0) + xc(2,0))
    - 2*normal[1]*(xc(0,2) - xc(2,2)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) +
    xc(3,0)) + xc(1,0)*(xc(2,2) - xc(3,2)) + xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(1,5)=
    (4*normpowfour*(-xc(2,0) + xc(3,0)) -
    2*normsquare*(-2*normal[3]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(-xc(1,2) +
    xc(2,2)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) - 4*normsquare*(-2*xc(1,2)*xc(2,1) +
    xc(0,2)*(-xc(1,1) + xc(2,1)) + 2*xc(0,1)*(xc(1,2) - xc(2,2)) + xc(1,1)*xc(2,2) +
    xc(2,1)*xc(2,2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) +
    xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) + 3*(2*normal[2]*(-xc(0,0) + xc(2,0))
    - 2*normal[1]*(-xc(0,1) + xc(2,1)))*(-2*normal[3]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(-xc(1,2) + xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal[2]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(-xc(0,1) +
    xc(2,1)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) + xc(3,0)) + xc(1,0)*(xc(2,2)
    - xc(3,2)) + xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(1,6)=
    (-4*normsquare*(xc(1,0)*xc(1,1) - 2*xc(0,1)*(xc(1,0) - xc(2,0)) -
    2*xc(1,1)*xc(2,0) + xc(0,0)*(xc(1,1) - xc(2,1)) +
    xc(1,0)*xc(2,1))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) +
    xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) + 3*(-2*normal[3]*(-xc(0,1) + xc(1,1))
    + 2*normal[2]*(-xc(0,2) + xc(1,2)))*(-2*normal[3]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(-xc(1,2) + xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal[3]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(-xc(1,2) +
    xc(2,2)))*((xc(0,2) - xc(1,2))*(xc(0,1) - xc(3,1)) + (-xc(0,1) +
    xc(1,1))*(xc(0,2) - xc(3,2))) + 4*normpowfour*(-xc(1,2) + xc(3,2)) -
    2*normsquare*(-2*normal[3]*(-xc(0,1) + xc(1,1)) + 2*normal[2]*(-xc(0,2) +
    xc(1,2)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) + xc(3,0)) + xc(1,0)*(xc(2,2)
    - xc(3,2)) + xc(2,0)*xc(3,2)))/(4.*normpowfive);
  
  elematrix(1,7)=
    (-2*normsquare*(2*(xc(0,0) - xc(1,0))*(xc(1,0) - xc(2,0)) + 2*(xc(0,2) -
    xc(1,2))*(xc(1,2) - xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) +
    3*(-2*normal[3]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(-xc(0,2) +
    xc(1,2)))*(-2*normal[3]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(-xc(1,2) +
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(xc(1,0) - xc(2,0))
    - 2*normal[1]*(-xc(1,2) + xc(2,2)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) +
    (xc(0,0) - xc(1,0))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(xc(0,0) -
    xc(1,0)) - 2*normal[1]*(-xc(0,2) + xc(1,2)))*(-(xc(2,2)*xc(3,0)) +
    xc(1,2)*(-xc(2,0) + xc(3,0)) + xc(1,0)*(xc(2,2) - xc(3,2)) +
    xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(1,8)=
    (4*normpowfour*(xc(1,0) - xc(3,0)) - 2*normsquare*(-2*normal[3]*(xc(1,0)
    - xc(2,0)) - 2*normal[1]*(-xc(1,2) + xc(2,2)))*((xc(0,1) - xc(1,1))*(xc(0,0) -
    xc(3,0)) + (xc(0,0) - xc(1,0))*(-xc(0,1) + xc(3,1))) -
    4*normsquare*(xc(1,1)*xc(1,2) + xc(0,2)*(xc(1,1) - xc(2,1)) + xc(1,2)*xc(2,1) -
    2*xc(0,1)*(xc(1,2) - xc(2,2)) - 2*xc(1,1)*xc(2,2))*(-(normal[1]*(xc(0,0) -
    xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) +
    3*(2*normal[2]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(xc(0,1) -
    xc(1,1)))*(-2*normal[3]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(-xc(1,2) +
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal[2]*(xc(0,0) - xc(1,0)) -
    2*normal[1]*(xc(0,1) - xc(1,1)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) +
    xc(3,0)) + xc(1,0)*(xc(2,2) - xc(3,2)) +xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(1,9)=
    (normal[2]*(xc(0,2)*xc(1,0)*xc(1,2) + pow(xc(1,1),2)*xc(2,0) -
    xc(0,2)*xc(1,2)*xc(2,0) + pow(xc(1,2),2)*xc(2,0) + xc(0,1)*(xc(1,0) -
    xc(2,0))*(xc(1,1) - xc(2,1)) - xc(1,0)*xc(1,1)*xc(2,1) - xc(1,1)*xc(2,0)*xc(2,1)
    + xc(1,0)*pow(xc(2,1),2) - xc(0,2)*xc(1,0)*xc(2,2) - xc(1,0)*xc(1,2)*xc(2,2) +
    xc(0,2)*xc(2,0)*xc(2,2) - xc(1,2)*xc(2,0)*xc(2,2) + xc(1,0)*pow(xc(2,2),2) -
    xc(0,0)*(pow(xc(1,1),2) + pow(xc(1,2),2) - 2*xc(1,1)*xc(2,1) +
    pow(xc(2,1),2) - 2*xc(1,2)*xc(2,2) + pow(xc(2,2),2))))/normcube;

  elematrix(1,10)=
    -(normal[2]*(-2*normal[3]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(-xc(1,2) +
    xc(2,2))))/(2.*normcube);

  elematrix(1,11)=
    -((normal[2]*(-(xc(0,1)*xc(1,1)*xc(1,2)) + xc(1,0)*xc(1,2)*xc(2,0) -
    xc(1,2)*pow(xc(2,0),2) + xc(0,1)*xc(1,2)*xc(2,1) + xc(1,1)*xc(1,2)*xc(2,1) -
    xc(1,2)*pow(xc(2,1),2) + xc(0,2)*(pow(xc(1,0),2) + pow(xc(1,1),2) -
    2*xc(1,0)*xc(2,0) + pow(xc(2,0),2) - 2*xc(1,1)*xc(2,1) + pow(xc(2,1),2)) -
    xc(0,0)*(xc(1,0) - xc(2,0))*(xc(1,2) - xc(2,2)) - pow(xc(1,0),2)*xc(2,2) +
    xc(0,1)*xc(1,1)*xc(2,2) - pow(xc(1,1),2)*xc(2,2) + xc(1,0)*xc(2,0)*xc(2,2) -
    xc(0,1)*xc(2,1)*xc(2,2) + xc(1,1)*xc(2,1)*xc(2,2)))/normcube);

  elematrix(2,0)=
    (-2*normsquare*(-2*normal[3]*(-xc(1,1) + xc(2,1)) + 2*normal[2]*(-xc(1,2) +
    xc(2,2)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) - 4*normsquare*(xc(1,0) - xc(2,0))*(-xc(1,2) +
    xc(2,2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) + 3*(2*normal[2]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(xc(1,1) - xc(2,1)))*(-2*normal[3]*(-xc(1,1) + xc(2,1)) +
    2*normal[2]*(-xc(1,2) + xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal[2]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(xc(1,1) -
    xc(2,1)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) +
    xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(2,1)=
    (-2*normsquare*(-2*normal[3]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(-xc(1,2) +
    xc(2,2)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) - 4*normsquare*(xc(1,1) - xc(2,1))*(-xc(1,2) +
    xc(2,2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) + 3*(2*normal[2]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(xc(1,1) - xc(2,1)))*(-2*normal[3]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(-xc(1,2) + xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal[2]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(xc(1,1) -
    xc(2,1)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) + xc(3,0)) + xc(1,0)*(xc(2,2)
    - xc(3,2)) + xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(2,2)=
    (-4*normsquare*(2*normal[2]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(xc(1,1) -
    xc(2,1)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) + 3*pow(2*normal[2]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(xc(1,1) - xc(2,1)),2)*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    4*normsquare*(pow(xc(1,0) - xc(2,0),2) + pow(xc(1,1) -
    xc(2,1),2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(2,3)=
    (4*normpowfour*(-xc(2,1) + xc(3,1)) -
    2*normsquare*(-2*normal[3]*(xc(0,1) - xc(2,1)) + 2*normal[2]*(xc(0,2) -
    xc(2,2)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) + 3*(2*normal[2]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(xc(1,1) - xc(2,1)))*(-2*normal[3]*(xc(0,1) - xc(2,1)) +
    2*normal[2]*(xc(0,2) - xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    4*normsquare*(2*xc(0,2)*(xc(1,0) - xc(2,0)) + xc(1,2)*xc(2,0) -
    2*xc(1,0)*xc(2,2) + xc(2,0)*xc(2,2) + xc(0,0)*(-xc(1,2) +
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal[2]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(xc(1,1) - xc(2,1)))*((xc(0,2) - xc(2,2))*(-xc(0,1) + xc(3,1)) +
    (xc(0,1) - xc(2,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(2,4)=
    (4*normpowfour*(xc(2,0) - xc(3,0)) -
    2*normsquare*(-2*normal[3]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(xc(0,2) -
    xc(2,2)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) + 3*(2*normal[2]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(xc(1,1) - xc(2,1)))*(-2*normal[3]*(-xc(0,0) + xc(2,0)) -
    2*normal[1]*(xc(0,2) - xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    4*normsquare*(2*xc(0,2)*(xc(1,1) - xc(2,1)) + xc(1,2)*xc(2,1) -
    2*xc(1,1)*xc(2,2) + xc(2,1)*xc(2,2) + xc(0,1)*(-xc(1,2) +
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal[2]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(xc(1,1) - xc(2,1)))*((xc(0,2) - xc(2,2))*(xc(0,0) - xc(3,0)) +
    (-xc(0,0) + xc(2,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(2,5)=
    (-2*normsquare*(2*normal[2]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(xc(1,1) -
    xc(2,1)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) - 2*normsquare*(2*normal[2]*(-xc(0,0) + xc(2,0)) -
    2*normal[1]*(-xc(0,1) + xc(2,1)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0)
    - xc(2,0)*xc(3,1) + xc(1,0)*(-xc(2,1) + xc(3,1))) + 3*(2*normal[2]*(xc(1,0) -
    xc(2,0)) - 2*normal[1]*(xc(1,1) - xc(2,1)))*(2*normal[2]*(-xc(0,0) + xc(2,0)) -
    2*normal[1]*(-xc(0,1) + xc(2,1)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*(xc(1,0) - xc(2,0))*(-xc(0,0) + xc(2,0)) + 2*(xc(1,1) -
    xc(2,1))*(-xc(0,1) + xc(2,1)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) -
    xc(3,2))))/(4.*normpowfive);

  elematrix(2,6)=
    (4*normpowfour*(xc(1,1) - xc(3,1)) -
    2*normsquare*(-2*normal[3]*(-xc(0,1) + xc(1,1)) + 2*normal[2]*(-xc(0,2) +
    xc(1,2)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) + 3*(-2*normal[3]*(-xc(0,1) + xc(1,1)) +
    2*normal[2]*(-xc(0,2) + xc(1,2)))*(2*normal[2]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(xc(1,1) - xc(2,1)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    4*normsquare*(xc(1,0)*xc(1,2) - 2*xc(0,2)*(xc(1,0) - xc(2,0)) -
    2*xc(1,2)*xc(2,0) + xc(0,0)*(xc(1,2) - xc(2,2)) +
    xc(1,0)*xc(2,2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) +
    xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal[2]*(xc(1,0) -
    xc(2,0)) - 2*normal[1]*(xc(1,1) - xc(2,1)))*((xc(0,2) - xc(1,2))*(xc(0,1) -
    xc(3,1)) + (-xc(0,1) + xc(1,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(2,7)=
    (4*normpowfour*(-xc(1,0) + xc(3,0)) -
    2*normsquare*(-2*normal[3]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(-xc(0,2) +
    xc(1,2)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) + 3*(-2*normal[3]*(xc(0,0) - xc(1,0)) -
    2*normal[1]*(-xc(0,2) + xc(1,2)))*(2*normal[2]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(xc(1,1) - xc(2,1)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    4*normsquare*(xc(1,1)*xc(1,2) - 2*xc(0,2)*(xc(1,1) - xc(2,1)) -
    2*xc(1,2)*xc(2,1) + xc(0,1)*(xc(1,2) - xc(2,2)) +
    xc(1,1)*xc(2,2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) +
    xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal[2]*(xc(1,0) -
    xc(2,0)) - 2*normal[1]*(xc(1,1) - xc(2,1)))*((-xc(0,2) + xc(1,2))*(xc(0,0) -
    xc(3,0)) + (xc(0,0) - xc(1,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(2,8)=
    (-2*normsquare*(2*normal[2]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(xc(1,1) -
    xc(2,1)))*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(-xc(0,1) + xc(3,1))) - 2*normsquare*(2*normal[2]*(xc(0,0) - xc(1,0)) -
    2*normal[1]*(xc(0,1) - xc(1,1)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0)
    - xc(2,0)*xc(3,1) + xc(1,0)*(-xc(2,1) + xc(3,1))) + 3*(2*normal[2]*(xc(0,0) -
    xc(1,0)) - 2*normal[1]*(xc(0,1) - xc(1,1)))*(2*normal[2]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(xc(1,1) - xc(2,1)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*(xc(0,0) - xc(1,0))*(xc(1,0) - xc(2,0)) + 2*(xc(0,1) -
    xc(1,1))*(xc(1,1) - xc(2,1)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) -
    xc(3,2))))/(4.*normpowfive);

  elematrix(2,9)=
    (normal[3]*(xc(0,2)*xc(1,0)*xc(1,2) + pow(xc(1,1),2)*xc(2,0) -
    xc(0,2)*xc(1,2)*xc(2,0) + pow(xc(1,2),2)*xc(2,0) + xc(0,1)*(xc(1,0) -
    xc(2,0))*(xc(1,1) - xc(2,1)) - xc(1,0)*xc(1,1)*xc(2,1) - xc(1,1)*xc(2,0)*xc(2,1)
    + xc(1,0)*pow(xc(2,1),2) - xc(0,2)*xc(1,0)*xc(2,2) - xc(1,0)*xc(1,2)*xc(2,2) +
    xc(0,2)*xc(2,0)*xc(2,2) - xc(1,2)*xc(2,0)*xc(2,2) + xc(1,0)*pow(xc(2,2),2) -
    xc(0,0)*(pow(xc(1,1),2) + pow(xc(1,2),2) - 2*xc(1,1)*xc(2,1) +
    pow(xc(2,1),2) - 2*xc(1,2)*xc(2,2) + pow(xc(2,2),2))))/normcube;

  elematrix(2,10)=
    -((normal[3]*(-(xc(0,2)*xc(1,1)*xc(1,2)) + xc(1,0)*xc(1,1)*xc(2,0) -
    xc(1,1)*pow(xc(2,0),2) - xc(0,0)*(xc(1,0) - xc(2,0))*(xc(1,1) - xc(2,1)) -
    pow(xc(1,0),2)*xc(2,1) + xc(0,2)*xc(1,2)*xc(2,1) - pow(xc(1,2),2)*xc(2,1) +
    xc(1,0)*xc(2,0)*xc(2,1) + xc(0,2)*xc(1,1)*xc(2,2) + xc(1,1)*xc(1,2)*xc(2,2) -
    xc(0,2)*xc(2,1)*xc(2,2) + xc(1,2)*xc(2,1)*xc(2,2) - xc(1,1)*pow(xc(2,2),2) +
    xc(0,1)*(pow(xc(1,0),2) + pow(xc(1,2),2) - 2*xc(1,0)*xc(2,0) +
    pow(xc(2,0),2) - 2*xc(1,2)*xc(2,2) +pow(xc(2,2),2))))/normcube);

  elematrix(2,11)=
    -((2*normal[2]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(xc(1,1) -
    xc(2,1)))*(-(xc(1,1)*xc(2,0)) + xc(0,1)*(-xc(1,0) + xc(2,0)) + xc(0,0)*(xc(1,1)
    - xc(2,1)) + xc(1,0)*xc(2,1)))/(2.*normcube);

  elematrix(3,0)=
    (3*(-2*normal[3]*(xc(0,1) - xc(2,1)) + 2*normal[2]*(xc(0,2) -
    xc(2,2)))*(-2*normal[3]*(-xc(1,1) + xc(2,1)) + 2*normal[2]*(-xc(1,2) +
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*(xc(0,1) - xc(2,1))*(-xc(1,1) +
    xc(2,1)) + 2*(xc(0,2) - xc(2,2))*(-xc(1,2) + xc(2,2)))*(-(normal[1]*(xc(0,0) -
    xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal[3]*(-xc(1,1) + xc(2,1)) + 2*normal[2]*(-xc(1,2) +
    xc(2,2)))*((xc(0,2) - xc(2,2))*(-xc(0,1) + xc(3,1)) + (xc(0,1) -
    xc(2,1))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(xc(0,1) - xc(2,1)) +
    2*normal[2]*(xc(0,2) - xc(2,2)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1)
    - xc(2,1)*xc(3,2) + xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(3,1)=
    (-4*normsquare*(2*xc(0,1)*(xc(1,0) - xc(2,0)) + xc(1,1)*xc(2,0) -
    2*xc(1,0)*xc(2,1) + xc(2,0)*xc(2,1) + xc(0,0)*(-xc(1,1) +
    xc(2,1)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) + 3*(-2*normal[3]*(xc(0,1) - xc(2,1)) +
    2*normal[2]*(xc(0,2) - xc(2,2)))*(-2*normal[3]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(-xc(1,2) + xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal[3]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(-xc(1,2) +
    xc(2,2)))*((xc(0,2) - xc(2,2))*(-xc(0,1) + xc(3,1)) + (xc(0,1) -
    xc(2,1))*(xc(0,2) - xc(3,2))) + 4*normpowfour*(xc(2,2) - xc(3,2)) -
    2*normsquare*(-2*normal[3]*(xc(0,1) - xc(2,1)) + 2*normal[2]*(xc(0,2) -
    xc(2,2)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) + xc(3,0)) + xc(1,0)*(xc(2,2)
    - xc(3,2)) + xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(3,2)=
    (4*normpowfour*(-xc(2,1) + xc(3,1)) -
    2*normsquare*(-2*normal[3]*(xc(0,1) - xc(2,1)) + 2*normal[2]*(xc(0,2) -
    xc(2,2)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) + 3*(2*normal[2]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(xc(1,1) - xc(2,1)))*(-2*normal[3]*(xc(0,1) - xc(2,1)) +
    2*normal[2]*(xc(0,2) - xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    4*normsquare*(2*xc(0,2)*(xc(1,0) - xc(2,0)) + xc(1,2)*xc(2,0) -
    2*xc(1,0)*xc(2,2) + xc(2,0)*xc(2,2) + xc(0,0)*(-xc(1,2) +
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal[2]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(xc(1,1) - xc(2,1)))*((xc(0,2) - xc(2,2))*(-xc(0,1) + xc(3,1)) +
    (xc(0,1) - xc(2,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(3,3)=
    (3*pow(-2*normal[3]*(xc(0,1) - xc(2,1)) + 2*normal[2]*(xc(0,2) -
    xc(2,2)),2)*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 4*normsquare*(pow(xc(0,1) - xc(2,1),2) +
    pow(xc(0,2) - xc(2,2),2))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-2*normal[3]*(xc(0,1) - xc(2,1)) + 2*normal[2]*(xc(0,2) -
    xc(2,2)))*((xc(0,2) - xc(2,2))*(-xc(0,1) + xc(3,1)) + (xc(0,1) -
    xc(2,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(3,4)=
    (-4*normsquare*(-xc(0,0) + xc(2,0))*(xc(0,1) - xc(2,1))*(-(normal[1]*(xc(0,0) -
    xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) +
    3*(-2*normal[3]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(xc(0,2) -
    xc(2,2)))*(-2*normal[3]*(xc(0,1) - xc(2,1)) + 2*normal[2]*(xc(0,2) -
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(xc(0,1) - xc(2,1))
    + 2*normal[2]*(xc(0,2) - xc(2,2)))*((xc(0,2) - xc(2,2))*(xc(0,0) - xc(3,0)) +
    (-xc(0,0) + xc(2,0))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(-xc(0,0)
    + xc(2,0)) - 2*normal[1]*(xc(0,2) - xc(2,2)))*((xc(0,2) - xc(2,2))*(-xc(0,1) +
    xc(3,1)) + (xc(0,1) - xc(2,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(3,5)=
    (-2*normsquare*(-2*normal[3]*(xc(0,1) - xc(2,1)) + 2*normal[2]*(xc(0,2) -
    xc(2,2)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) + 3*(2*normal[2]*(-xc(0,0) + xc(2,0)) -
    2*normal[1]*(-xc(0,1) + xc(2,1)))*(-2*normal[3]*(xc(0,1) - xc(2,1)) +
    2*normal[2]*(xc(0,2) - xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-xc(0,0) + xc(2,0))*(xc(0,2) - xc(2,2))*(-(normal[1]*(xc(0,0) -
    xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal[2]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(-xc(0,1) +
    xc(2,1)))*((xc(0,2) - xc(2,2))*(-xc(0,1) + xc(3,1)) + (xc(0,1) -
    xc(2,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(3,6)=
    (3*(-2*normal[3]*(-xc(0,1) + xc(1,1)) + 2*normal[2]*(-xc(0,2) +
    xc(1,2)))*(-2*normal[3]*(xc(0,1) - xc(2,1)) + 2*normal[2]*(xc(0,2) -
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*(-xc(0,1) + xc(1,1))*(xc(0,1) -
    xc(2,1)) + 2*(-xc(0,2) + xc(1,2))*(xc(0,2) - xc(2,2)))*(-(normal[1]*(xc(0,0) -
    xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal[3]*(xc(0,1) - xc(2,1)) + 2*normal[2]*(xc(0,2) -
    xc(2,2)))*((xc(0,2) - xc(1,2))*(xc(0,1) - xc(3,1)) + (-xc(0,1) +
    xc(1,1))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(-xc(0,1) + xc(1,1))
    + 2*normal[2]*(-xc(0,2) + xc(1,2)))*((xc(0,2) - xc(2,2))*(-xc(0,1) + xc(3,1)) +
    (xc(0,1) - xc(2,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(3,7)=
    (-4*normsquare*(-(xc(1,1)*xc(2,0)) + xc(0,1)*(-2*xc(1,0) + xc(2,0)) +
    xc(0,0)*(xc(0,1) + xc(1,1) - 2*xc(2,1)) +
    2*xc(1,0)*xc(2,1))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) +
    xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) + 3*(-2*normal[3]*(xc(0,0) - xc(1,0))
    - 2*normal[1]*(-xc(0,2) + xc(1,2)))*(-2*normal[3]*(xc(0,1) - xc(2,1)) +
    2*normal[2]*(xc(0,2) - xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal[3]*(xc(0,1) - xc(2,1)) + 2*normal[2]*(xc(0,2) -
    xc(2,2)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(xc(0,0) - xc(1,0)) -
    2*normal[1]*(-xc(0,2) + xc(1,2)))*((xc(0,2) - xc(2,2))*(-xc(0,1) + xc(3,1)) +
    (xc(0,1) - xc(2,1))*(xc(0,2) - xc(3,2))) + 4*normpowfour*(-xc(0,2) +
    xc(3,2)))/(4.*normpowfive);

  elematrix(3,8)=
    (4*normpowfour*(xc(0,1) - xc(3,1)) - 2*normsquare*(-2*normal[3]*(xc(0,1)
    - xc(2,1)) + 2*normal[2]*(xc(0,2) - xc(2,2)))*((xc(0,1) - xc(1,1))*(xc(0,0) -
    xc(3,0)) + (xc(0,0) - xc(1,0))*(-xc(0,1) + xc(3,1))) + 3*(2*normal[2]*(xc(0,0) -
    xc(1,0)) - 2*normal[1]*(xc(0,1) - xc(1,1)))*(-2*normal[3]*(xc(0,1) - xc(2,1)) +
    2*normal[2]*(xc(0,2) - xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-(xc(1,2)*xc(2,0)) + xc(0,2)*(-2*xc(1,0) + xc(2,0)) +
    xc(0,0)*(xc(0,2) + xc(1,2) - 2*xc(2,2)) +
    2*xc(1,0)*xc(2,2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) +
    xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal[2]*(xc(0,0) -
    xc(1,0)) - 2*normal[1]*(xc(0,1) - xc(1,1)))*((xc(0,2) - xc(2,2))*(-xc(0,1) +
    xc(3,1)) + (xc(0,1) - xc(2,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(3,9)=
    -((-2*normal[3]*(xc(0,1) - xc(2,1)) + 2*normal[2]*(xc(0,2) -
    xc(2,2)))*(-(xc(1,2)*xc(2,1)) + xc(0,2)*(-xc(1,1) + xc(2,1)) + xc(0,1)*(xc(1,2)
    - xc(2,2)) + xc(1,1)*xc(2,2)))/(2.*normcube);

  elematrix(3,10)=
    -((normal[1]*(xc(0,1)*xc(1,0)*xc(2,0) - xc(0,1)*pow(xc(2,0),2) +
    xc(1,1)*pow(xc(2,0),2) + pow(xc(0,0),2)*(xc(1,1) - xc(2,1)) +
    pow(xc(0,2),2)*(xc(1,1) - xc(2,1)) - xc(1,0)*xc(2,0)*xc(2,1) +
    xc(0,0)*(-2*xc(1,1)*xc(2,0) + xc(0,1)*(-xc(1,0) + xc(2,0)) + (xc(1,0) +
    xc(2,0))*xc(2,1)) + xc(0,1)*xc(1,2)*xc(2,2) - xc(1,2)*xc(2,1)*xc(2,2) -
    xc(0,1)*pow(xc(2,2),2) + xc(1,1)*pow(xc(2,2),2) + xc(0,2)*(xc(1,2)*xc(2,1) +
    (-2*xc(1,1) + xc(2,1))*xc(2,2) + xc(0,1)*(-xc(1,2) +
    xc(2,2)))))/normcube);

  elematrix(3,11)=
    -((normal[1]*(xc(0,2)*xc(1,0)*xc(2,0) - xc(0,2)*pow(xc(2,0),2) +
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
    xc(2,0)*xc(2,1))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) +
    xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) + 3*(-2*normal[3]*(-xc(0,0) + xc(2,0))
    - 2*normal[1]*(xc(0,2) - xc(2,2)))*(-2*normal[3]*(-xc(1,1) + xc(2,1)) +
    2*normal[2]*(-xc(1,2) + xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal[3]*(-xc(1,1) + xc(2,1)) + 2*normal[2]*(-xc(1,2) +
    xc(2,2)))*((xc(0,2) - xc(2,2))*(xc(0,0) - xc(3,0)) + (-xc(0,0) +
    xc(2,0))*(xc(0,2) - xc(3,2))) + 4*normpowfour*(-xc(2,2) + xc(3,2)) -
    2*normsquare*(-2*normal[3]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(xc(0,2) -
    xc(2,2)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) +
    xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(4,1)=
    (3*(-2*normal[3]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(xc(0,2) -
    xc(2,2)))*(-2*normal[3]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(-xc(1,2) +
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*(xc(1,0) - xc(2,0))*(-xc(0,0) +
    xc(2,0)) + 2*(xc(0,2) - xc(2,2))*(-xc(1,2) + xc(2,2)))*(-(normal[1]*(xc(0,0) -
    xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal[3]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(-xc(1,2) +
    xc(2,2)))*((xc(0,2) - xc(2,2))*(xc(0,0) - xc(3,0)) + (-xc(0,0) +
    xc(2,0))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(-xc(0,0) + xc(2,0))
    - 2*normal[1]*(xc(0,2) - xc(2,2)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) +
    xc(3,0)) + xc(1,0)*(xc(2,2) - xc(3,2)) +xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(4,2)=
    (4*normpowfour*(xc(2,0) - xc(3,0)) -
    2*normsquare*(-2*normal[3]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(xc(0,2) -
    xc(2,2)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) + 3*(2*normal[2]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(xc(1,1) - xc(2,1)))*(-2*normal[3]*(-xc(0,0) + xc(2,0)) -
    2*normal[1]*(xc(0,2) - xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    4*normsquare*(2*xc(0,2)*(xc(1,1) - xc(2,1)) + xc(1,2)*xc(2,1) -
    2*xc(1,1)*xc(2,2) + xc(2,1)*xc(2,2) + xc(0,1)*(-xc(1,2) +
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal[2]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(xc(1,1) - xc(2,1)))*((xc(0,2) - xc(2,2))*(xc(0,0) - xc(3,0)) +
    (-xc(0,0) + xc(2,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(4,3)=
    (-4*normsquare*(-xc(0,0) + xc(2,0))*(xc(0,1) - xc(2,1))*(-(normal[1]*(xc(0,0) -
    xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) +
    3*(-2*normal[3]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(xc(0,2) -
    xc(2,2)))*(-2*normal[3]*(xc(0,1) - xc(2,1)) + 2*normal[2]*(xc(0,2) -
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(xc(0,1) - xc(2,1))
    + 2*normal[2]*(xc(0,2) - xc(2,2)))*((xc(0,2) - xc(2,2))*(xc(0,0) - xc(3,0)) +
    (-xc(0,0) + xc(2,0))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(-xc(0,0)
    + xc(2,0)) - 2*normal[1]*(xc(0,2) - xc(2,2)))*((xc(0,2) - xc(2,2))*(-xc(0,1) +
    xc(3,1)) + (xc(0,1) - xc(2,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(4,4)=
    (3*pow(-2*normal[3]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(xc(0,2) -
    xc(2,2)),2)*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 4*normsquare*(pow(xc(0,0) - xc(2,0),2) +
    pow(xc(0,2) - xc(2,2),2))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-2*normal[3]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(xc(0,2) -
    xc(2,2)))*((xc(0,2) - xc(2,2))*(xc(0,0) - xc(3,0)) + (-xc(0,0) +
    xc(2,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(4,5)=
    (-2*normsquare*(-2*normal[3]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(xc(0,2) -
    xc(2,2)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) + 3*(2*normal[2]*(-xc(0,0) + xc(2,0)) -
    2*normal[1]*(-xc(0,1) + xc(2,1)))*(-2*normal[3]*(-xc(0,0) + xc(2,0)) -
    2*normal[1]*(xc(0,2) - xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-xc(0,1) + xc(2,1))*(xc(0,2) - xc(2,2))*(-(normal[1]*(xc(0,0) -
    xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal[2]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(-xc(0,1) +
    xc(2,1)))*((xc(0,2) - xc(2,2))*(xc(0,0) - xc(3,0)) + (-xc(0,0) +
    xc(2,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(4,6)=
    (-4*normsquare*(xc(0,1)*(xc(1,0) - 2*xc(2,0)) + 2*xc(1,1)*xc(2,0) -
    xc(1,0)*xc(2,1) + xc(0,0)*(xc(0,1) - 2*xc(1,1) + xc(2,1)))*(-(normal[1]*(xc(0,0)
    - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) +
    3*(-2*normal[3]*(-xc(0,1) + xc(1,1)) + 2*normal[2]*(-xc(0,2) +
    xc(1,2)))*(-2*normal[3]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(xc(0,2) -
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(-xc(0,0) + xc(2,0))
    - 2*normal[1]*(xc(0,2) - xc(2,2)))*((xc(0,2) - xc(1,2))*(xc(0,1) - xc(3,1)) +
    (-xc(0,1) + xc(1,1))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(-xc(0,1)
    + xc(1,1)) + 2*normal[2]*(-xc(0,2) + xc(1,2)))*((xc(0,2) - xc(2,2))*(xc(0,0) -
    xc(3,0)) + (-xc(0,0) + xc(2,0))*(xc(0,2) - xc(3,2))) +
    4*normpowfour*(xc(0,2) - xc(3,2)))/(4.*normpowfive);

  elematrix(4,7)=
    (3*(-2*normal[3]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(-xc(0,2) +
    xc(1,2)))*(-2*normal[3]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(xc(0,2) -
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*(xc(0,0) - xc(1,0))*(-xc(0,0) +
    xc(2,0)) + 2*(-xc(0,2) + xc(1,2))*(xc(0,2) - xc(2,2)))*(-(normal[1]*(xc(0,0) -
    xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal[3]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(xc(0,2) -
    xc(2,2)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(xc(0,0) - xc(1,0)) -
    2*normal[1]*(-xc(0,2) + xc(1,2)))*((xc(0,2) - xc(2,2))*(xc(0,0) - xc(3,0)) +
    (-xc(0,0) + xc(2,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(4,8)=
    (4*normpowfour*(-xc(0,0) + xc(3,0)) -
    2*normsquare*(-2*normal[3]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(xc(0,2) -
    xc(2,2)))*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(-xc(0,1) + xc(3,1))) + 3*(2*normal[2]*(xc(0,0) - xc(1,0)) -
    2*normal[1]*(xc(0,1) - xc(1,1)))*(-2*normal[3]*(-xc(0,0) + xc(2,0)) -
    2*normal[1]*(xc(0,2) - xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-(xc(1,2)*xc(2,1)) + xc(0,2)*(-2*xc(1,1) + xc(2,1)) +
    xc(0,1)*(xc(0,2) + xc(1,2) - 2*xc(2,2)) +
    2*xc(1,1)*xc(2,2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) +
    xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal[2]*(xc(0,0) -
    xc(1,0)) - 2*normal[1]*(xc(0,1) - xc(1,1)))*((xc(0,2) - xc(2,2))*(xc(0,0) -
    xc(3,0)) + (-xc(0,0) + xc(2,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(4,9)=
    -((normal[2]*(pow(xc(0,1),2)*(xc(1,0) - xc(2,0)) + pow(xc(0,2),2)*(xc(1,0) -
    xc(2,0)) + xc(0,0)*xc(1,1)*xc(2,1) - xc(1,1)*xc(2,0)*xc(2,1) -
    xc(0,0)*pow(xc(2,1),2) + xc(1,0)*pow(xc(2,1),2) + xc(0,1)*(xc(1,1)*xc(2,0) +
    (-2*xc(1,0) + xc(2,0))*xc(2,1) + xc(0,0)*(-xc(1,1) + xc(2,1))) +
    xc(0,0)*xc(1,2)*xc(2,2) - xc(1,2)*xc(2,0)*xc(2,2) - xc(0,0)*pow(xc(2,2),2) +
    xc(1,0)*pow(xc(2,2),2) + xc(0,2)*(xc(1,2)*xc(2,0) + (-2*xc(1,0) +
    xc(2,0))*xc(2,2) + xc(0,0)*(-xc(1,2) + xc(2,2)))))/normcube);

  elematrix(4,10)=
    -(normal[2]*(-2*normal[3]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(xc(0,2) -
    xc(2,2))))/(2.*normcube);

  elematrix(4,11)=
    -((normal[2]*(xc(0,2)*xc(1,0)*xc(2,0) - xc(0,2)*pow(xc(2,0),2) +
    xc(1,2)*pow(xc(2,0),2) + xc(0,2)*xc(1,1)*xc(2,1) - xc(0,2)*pow(xc(2,1),2) +
    xc(1,2)*pow(xc(2,1),2) + pow(xc(0,0),2)*(xc(1,2) - xc(2,2)) +
    pow(xc(0,1),2)*(xc(1,2) - xc(2,2)) - xc(1,0)*xc(2,0)*xc(2,2) -
    xc(1,1)*xc(2,1)*xc(2,2) + xc(0,0)*(-2*xc(1,2)*xc(2,0) + xc(0,2)*(-xc(1,0) +
    xc(2,0)) + (xc(1,0) + xc(2,0))*xc(2,2)) + xc(0,1)*(-2*xc(1,2)*xc(2,1) +
    xc(0,2)*(-xc(1,1) + xc(2,1)) + (xc(1,1) +xc(2,1))*xc(2,2))))/normcube);

  elematrix(5,0)=
    (-2*normsquare*(-2*normal[3]*(-xc(1,1) + xc(2,1)) + 2*normal[2]*(-xc(1,2) +
    xc(2,2)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) + 4*normpowfour*(xc(2,1) - xc(3,1)) -
    4*normsquare*(-2*xc(1,2)*xc(2,0) + xc(0,2)*(-xc(1,0) + xc(2,0)) +
    2*xc(0,0)*(xc(1,2) - xc(2,2)) + xc(1,0)*xc(2,2) +
    xc(2,0)*xc(2,2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) +
    xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) + 3*(2*normal[2]*(-xc(0,0) + xc(2,0))
    - 2*normal[1]*(-xc(0,1) + xc(2,1)))*(-2*normal[3]*(-xc(1,1) + xc(2,1)) +
    2*normal[2]*(-xc(1,2) + xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal[2]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(-xc(0,1) +
    xc(2,1)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) +
    xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(5,1)=
    (4*normpowfour*(-xc(2,0) + xc(3,0)) -
    2*normsquare*(-2*normal[3]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(-xc(1,2) +
    xc(2,2)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) - 4*normsquare*(-2*xc(1,2)*xc(2,1) +
    xc(0,2)*(-xc(1,1) + xc(2,1)) + 2*xc(0,1)*(xc(1,2) - xc(2,2)) + xc(1,1)*xc(2,2) +
    xc(2,1)*xc(2,2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) +
    xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) + 3*(2*normal[2]*(-xc(0,0) + xc(2,0))
    - 2*normal[1]*(-xc(0,1) + xc(2,1)))*(-2*normal[3]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(-xc(1,2) + xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal[2]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(-xc(0,1) +
    xc(2,1)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) + xc(3,0)) + xc(1,0)*(xc(2,2)
    - xc(3,2)) + xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(5,2)=
    (-2*normsquare*(2*normal[2]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(xc(1,1) -
    xc(2,1)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) - 2*normsquare*(2*normal[2]*(-xc(0,0) + xc(2,0)) -
    2*normal[1]*(-xc(0,1) + xc(2,1)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0)
    - xc(2,0)*xc(3,1) + xc(1,0)*(-xc(2,1) + xc(3,1))) + 3*(2*normal[2]*(xc(1,0) -
    xc(2,0)) - 2*normal[1]*(xc(1,1) - xc(2,1)))*(2*normal[2]*(-xc(0,0) + xc(2,0)) -
    2*normal[1]*(-xc(0,1) + xc(2,1)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*(xc(1,0) - xc(2,0))*(-xc(0,0) + xc(2,0)) + 2*(xc(1,1) -
    xc(2,1))*(-xc(0,1) + xc(2,1)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) -
    xc(3,2))))/(4.*normpowfive);

  elematrix(5,3)=
    (-2*normsquare*(-2*normal[3]*(xc(0,1) - xc(2,1)) + 2*normal[2]*(xc(0,2) -
    xc(2,2)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) + 3*(2*normal[2]*(-xc(0,0) + xc(2,0)) -
    2*normal[1]*(-xc(0,1) + xc(2,1)))*(-2*normal[3]*(xc(0,1) - xc(2,1)) +
    2*normal[2]*(xc(0,2) - xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-xc(0,0) + xc(2,0))*(xc(0,2) - xc(2,2))*(-(normal[1]*(xc(0,0) -
    xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal[2]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(-xc(0,1) +
    xc(2,1)))*((xc(0,2) - xc(2,2))*(-xc(0,1) + xc(3,1)) + (xc(0,1) -
    xc(2,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(5,4)=
    (-2*normsquare*(-2*normal[3]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(xc(0,2) -
    xc(2,2)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) + 3*(2*normal[2]*(-xc(0,0) + xc(2,0)) -
    2*normal[1]*(-xc(0,1) + xc(2,1)))*(-2*normal[3]*(-xc(0,0) + xc(2,0)) -
    2*normal[1]*(xc(0,2) - xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-xc(0,1) + xc(2,1))*(xc(0,2) - xc(2,2))*(-(normal[1]*(xc(0,0) -
    xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal[2]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(-xc(0,1) +
    xc(2,1)))*((xc(0,2) - xc(2,2))*(xc(0,0) - xc(3,0)) + (-xc(0,0) +
    xc(2,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(5,5)=
    (-4*normsquare*(2*normal[2]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(-xc(0,1) +
    xc(2,1)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) - 4*normsquare*(pow(xc(0,0) - xc(2,0),2) +
    pow(xc(0,1) - xc(2,1),2))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) +
    3*pow(2*normal[2]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(-xc(0,1) +
    xc(2,1)),2)*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(5,6)=
    (-2*normsquare*(-2*normal[3]*(-xc(0,1) + xc(1,1)) + 2*normal[2]*(-xc(0,2) +
    xc(1,2)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) + 4*normpowfour*(-xc(0,1) + xc(3,1)) +
    3*(-2*normal[3]*(-xc(0,1) + xc(1,1)) + 2*normal[2]*(-xc(0,2) +
    xc(1,2)))*(2*normal[2]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(-xc(0,1) +
    xc(2,1)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 4*normsquare*(xc(0,2)*(xc(1,0) - 2*xc(2,0)) +
    2*xc(1,2)*xc(2,0) - xc(1,0)*xc(2,2) + xc(0,0)*(xc(0,2) - 2*xc(1,2) +
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal[2]*(-xc(0,0) + xc(2,0))
    - 2*normal[1]*(-xc(0,1) + xc(2,1)))*((xc(0,2) - xc(1,2))*(xc(0,1) - xc(3,1)) +
    (-xc(0,1) + xc(1,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(5,7)=
    (4*normpowfour*(xc(0,0) - xc(3,0)) - 2*normsquare*(-2*normal[3]*(xc(0,0)
    - xc(1,0)) - 2*normal[1]*(-xc(0,2) + xc(1,2)))*((-xc(0,1) + xc(2,1))*(xc(0,0) -
    xc(3,0)) + (xc(0,0) - xc(2,0))*(xc(0,1) - xc(3,1))) + 3*(-2*normal[3]*(xc(0,0) -
    xc(1,0)) - 2*normal[1]*(-xc(0,2) + xc(1,2)))*(2*normal[2]*(-xc(0,0) + xc(2,0)) -
    2*normal[1]*(-xc(0,1) + xc(2,1)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    4*normsquare*(xc(0,2)*(xc(1,1) - 2*xc(2,1)) + 2*xc(1,2)*xc(2,1) -
    xc(1,1)*xc(2,2) + xc(0,1)*(xc(0,2) - 2*xc(1,2) + xc(2,2)))*(-(normal[1]*(xc(0,0)
    - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal[2]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(-xc(0,1) +
    xc(2,1)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(5,8)=
    (-2*normsquare*(2*normal[2]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(xc(0,1) -
    xc(1,1)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) - 2*normsquare*(2*normal[2]*(-xc(0,0) + xc(2,0)) -
    2*normal[1]*(-xc(0,1) + xc(2,1)))*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) +
    (xc(0,0) - xc(1,0))*(-xc(0,1) + xc(3,1))) + 3*(2*normal[2]*(xc(0,0) - xc(1,0)) -
    2*normal[1]*(xc(0,1) - xc(1,1)))*(2*normal[2]*(-xc(0,0) + xc(2,0)) -
    2*normal[1]*(-xc(0,1) + xc(2,1)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*(xc(0,0) - xc(1,0))*(-xc(0,0) + xc(2,0)) + 2*(xc(0,1) -
    xc(1,1))*(-xc(0,1) + xc(2,1)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) -
    xc(3,2))))/(4.*normpowfive);

  elematrix(5,9)=
    -((normal[3]*(pow(xc(0,1),2)*(xc(1,0) - xc(2,0)) + pow(xc(0,2),2)*(xc(1,0) -
    xc(2,0)) + xc(0,0)*xc(1,1)*xc(2,1) - xc(1,1)*xc(2,0)*xc(2,1) -
    xc(0,0)*pow(xc(2,1),2) + xc(1,0)*pow(xc(2,1),2) + xc(0,1)*(xc(1,1)*xc(2,0) +
    (-2*xc(1,0) + xc(2,0))*xc(2,1) + xc(0,0)*(-xc(1,1) + xc(2,1))) +
    xc(0,0)*xc(1,2)*xc(2,2) - xc(1,2)*xc(2,0)*xc(2,2) - xc(0,0)*pow(xc(2,2),2) +
    xc(1,0)*pow(xc(2,2),2) + xc(0,2)*(xc(1,2)*xc(2,0) + (-2*xc(1,0) +
    xc(2,0))*xc(2,2) + xc(0,0)*(-xc(1,2) + xc(2,2)))))/normcube);

  elematrix(5,10)=
    -((normal[3]*(xc(0,1)*xc(1,0)*xc(2,0) - xc(0,1)*pow(xc(2,0),2) +
    xc(1,1)*pow(xc(2,0),2) + pow(xc(0,0),2)*(xc(1,1) - xc(2,1)) +
    pow(xc(0,2),2)*(xc(1,1) - xc(2,1)) - xc(1,0)*xc(2,0)*xc(2,1) +
    xc(0,0)*(-2*xc(1,1)*xc(2,0) + xc(0,1)*(-xc(1,0) + xc(2,0)) + (xc(1,0) +
    xc(2,0))*xc(2,1)) + xc(0,1)*xc(1,2)*xc(2,2) - xc(1,2)*xc(2,1)*xc(2,2) -
    xc(0,1)*pow(xc(2,2),2) + xc(1,1)*pow(xc(2,2),2) + xc(0,2)*(xc(1,2)*xc(2,1) +
    (-2*xc(1,1) + xc(2,1))*xc(2,2) + xc(0,1)*(-xc(1,2) +
    xc(2,2)))))/normcube);

  elematrix(5,11)=
    -((-(xc(1,1)*xc(2,0)) + xc(0,1)*(-xc(1,0) + xc(2,0)) + xc(0,0)*(xc(1,1) -
    xc(2,1)) + xc(1,0)*xc(2,1))*(2*normal[2]*(-xc(0,0) + xc(2,0)) -
    2*normal[1]*(-xc(0,1) + xc(2,1))))/(2.*normcube);

  elematrix(6,0)=
    (-2*normsquare*(2*(xc(0,1) - xc(1,1))*(xc(1,1) - xc(2,1)) + 2*(xc(0,2) -
    xc(1,2))*(xc(1,2) - xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) +
    3*(-2*normal[3]*(-xc(0,1) + xc(1,1)) + 2*normal[2]*(-xc(0,2) +
    xc(1,2)))*(-2*normal[3]*(-xc(1,1) + xc(2,1)) + 2*normal[2]*(-xc(1,2) +
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(-xc(1,1) + xc(2,1))
    + 2*normal[2]*(-xc(1,2) + xc(2,2)))*((xc(0,2) - xc(1,2))*(xc(0,1) - xc(3,1)) +
    (-xc(0,1) + xc(1,1))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(-xc(0,1)
    + xc(1,1)) + 2*normal[2]*(-xc(0,2) + xc(1,2)))*(xc(1,2)*(xc(2,1) - xc(3,1)) +
    xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) + xc(1,1)*(-xc(2,2) +
    xc(3,2))))/(4.*normpowfive);

  elematrix(6,1)=
    (-4*normsquare*(xc(1,0)*xc(1,1) - 2*xc(0,1)*(xc(1,0) - xc(2,0)) -
    2*xc(1,1)*xc(2,0) + xc(0,0)*(xc(1,1) - xc(2,1)) +
    xc(1,0)*xc(2,1))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) +
    xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) + 3*(-2*normal[3]*(-xc(0,1) + xc(1,1))
    + 2*normal[2]*(-xc(0,2) + xc(1,2)))*(-2*normal[3]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(-xc(1,2) + xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal[3]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(-xc(1,2) +
    xc(2,2)))*((xc(0,2) - xc(1,2))*(xc(0,1) - xc(3,1)) + (-xc(0,1) +
    xc(1,1))*(xc(0,2) - xc(3,2))) + 4*normpowfour*(-xc(1,2) + xc(3,2)) -
    2*normsquare*(-2*normal[3]*(-xc(0,1) + xc(1,1)) + 2*normal[2]*(-xc(0,2) +
    xc(1,2)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) + xc(3,0)) + xc(1,0)*(xc(2,2)
    - xc(3,2)) + xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(6,2)=
    (4*normpowfour*(xc(1,1) - xc(3,1)) -
    2*normsquare*(-2*normal[3]*(-xc(0,1) + xc(1,1)) + 2*normal[2]*(-xc(0,2) +
    xc(1,2)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) + 3*(-2*normal[3]*(-xc(0,1) + xc(1,1)) +
    2*normal[2]*(-xc(0,2) + xc(1,2)))*(2*normal[2]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(xc(1,1) - xc(2,1)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    4*normsquare*(xc(1,0)*xc(1,2) - 2*xc(0,2)*(xc(1,0) - xc(2,0)) -
    2*xc(1,2)*xc(2,0) + xc(0,0)*(xc(1,2) - xc(2,2)) +
    xc(1,0)*xc(2,2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) +
    xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal[2]*(xc(1,0) -
    xc(2,0)) - 2*normal[1]*(xc(1,1) - xc(2,1)))*((xc(0,2) - xc(1,2))*(xc(0,1) -
    xc(3,1)) + (-xc(0,1) + xc(1,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(6,3)=
    (3*(-2*normal[3]*(-xc(0,1) + xc(1,1)) + 2*normal[2]*(-xc(0,2) +
    xc(1,2)))*(-2*normal[3]*(xc(0,1) - xc(2,1)) + 2*normal[2]*(xc(0,2) -
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*(-xc(0,1) + xc(1,1))*(xc(0,1) -
    xc(2,1)) + 2*(-xc(0,2) + xc(1,2))*(xc(0,2) - xc(2,2)))*(-(normal[1]*(xc(0,0) -
    xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal[3]*(xc(0,1) - xc(2,1)) + 2*normal[2]*(xc(0,2) -
    xc(2,2)))*((xc(0,2) - xc(1,2))*(xc(0,1) - xc(3,1)) + (-xc(0,1) +
    xc(1,1))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(-xc(0,1) + xc(1,1))
    + 2*normal[2]*(-xc(0,2) + xc(1,2)))*((xc(0,2) - xc(2,2))*(-xc(0,1) + xc(3,1)) +
    (xc(0,1) - xc(2,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(6,4)=
    (-4*normsquare*(xc(0,1)*(xc(1,0) - 2*xc(2,0)) + 2*xc(1,1)*xc(2,0) -
    xc(1,0)*xc(2,1) + xc(0,0)*(xc(0,1) - 2*xc(1,1) + xc(2,1)))*(-(normal[1]*(xc(0,0)
    - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) +
    3*(-2*normal[3]*(-xc(0,1) + xc(1,1)) + 2*normal[2]*(-xc(0,2) +
    xc(1,2)))*(-2*normal[3]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(xc(0,2) -
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(-xc(0,0) + xc(2,0))
    - 2*normal[1]*(xc(0,2) - xc(2,2)))*((xc(0,2) - xc(1,2))*(xc(0,1) - xc(3,1)) +
    (-xc(0,1) + xc(1,1))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(-xc(0,1)
    + xc(1,1)) + 2*normal[2]*(-xc(0,2) + xc(1,2)))*((xc(0,2) - xc(2,2))*(xc(0,0) -
    xc(3,0)) + (-xc(0,0) + xc(2,0))*(xc(0,2) - xc(3,2))) +
    4*normpowfour*(xc(0,2) - xc(3,2)))/(4.*normpowfive);

  elematrix(6,5)=
    (-2*normsquare*(-2*normal[3]*(-xc(0,1) + xc(1,1)) + 2*normal[2]*(-xc(0,2) +
    xc(1,2)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) + 4*normpowfour*(-xc(0,1) + xc(3,1)) +
    3*(-2*normal[3]*(-xc(0,1) + xc(1,1)) + 2*normal[2]*(-xc(0,2) +
    xc(1,2)))*(2*normal[2]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(-xc(0,1) +
    xc(2,1)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 4*normsquare*(xc(0,2)*(xc(1,0) - 2*xc(2,0)) +
    2*xc(1,2)*xc(2,0) - xc(1,0)*xc(2,2) + xc(0,0)*(xc(0,2) - 2*xc(1,2) +
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal[2]*(-xc(0,0) + xc(2,0))
    - 2*normal[1]*(-xc(0,1) + xc(2,1)))*((xc(0,2) - xc(1,2))*(xc(0,1) - xc(3,1)) +
    (-xc(0,1) + xc(1,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(6,6)=
    (-4*normsquare*(pow(xc(0,1) - xc(1,1),2) + pow(xc(0,2) -
    xc(1,2),2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) + 3*pow(-2*normal[3]*(-xc(0,1) + xc(1,1)) +
    2*normal[2]*(-xc(0,2) + xc(1,2)),2)*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-2*normal[3]*(-xc(0,1) + xc(1,1)) + 2*normal[2]*(-xc(0,2) +
    xc(1,2)))*((xc(0,2) - xc(1,2))*(xc(0,1) - xc(3,1)) + (-xc(0,1) +
    xc(1,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(6,7)=
    (-4*normsquare*(xc(0,0) - xc(1,0))*(-xc(0,1) + xc(1,1))*(-(normal[1]*(xc(0,0) -
    xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) +
    3*(-2*normal[3]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(-xc(0,2) +
    xc(1,2)))*(-2*normal[3]*(-xc(0,1) + xc(1,1)) + 2*normal[2]*(-xc(0,2) +
    xc(1,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(-xc(0,1) + xc(1,1))
    + 2*normal[2]*(-xc(0,2) + xc(1,2)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) +
    (xc(0,0) - xc(1,0))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(xc(0,0) -
    xc(1,0)) - 2*normal[1]*(-xc(0,2) + xc(1,2)))*((xc(0,2) - xc(1,2))*(xc(0,1) -
    xc(3,1)) + (-xc(0,1) + xc(1,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(6,8)=
    (-2*normsquare*(-2*normal[3]*(-xc(0,1) + xc(1,1)) + 2*normal[2]*(-xc(0,2) +
    xc(1,2)))*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(-xc(0,1) + xc(3,1))) - 4*normsquare*(xc(0,0) - xc(1,0))*(-xc(0,2) +
    xc(1,2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) + 3*(2*normal[2]*(xc(0,0) - xc(1,0)) -
    2*normal[1]*(xc(0,1) - xc(1,1)))*(-2*normal[3]*(-xc(0,1) + xc(1,1)) +
    2*normal[2]*(-xc(0,2) + xc(1,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal[2]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(xc(0,1) -
    xc(1,1)))*((xc(0,2) - xc(1,2))*(xc(0,1) - xc(3,1)) + (-xc(0,1) +
    xc(1,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(6,9)=
    -((-2*normal[3]*(-xc(0,1) + xc(1,1)) + 2*normal[2]*(-xc(0,2) +
    xc(1,2)))*(-(xc(1,2)*xc(2,1)) + xc(0,2)*(-xc(1,1) + xc(2,1)) + xc(0,1)*(xc(1,2)
    - xc(2,2)) + xc(1,1)*xc(2,2)))/(2.*normcube);

  elematrix(6,10)=
    (normal[1]*(pow(xc(0,2),2)*xc(1,1) - xc(0,2)*xc(1,1)*xc(1,2) +
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
    - 2*xc(0,0)*(xc(1,1) - xc(2,1)) - 2*xc(1,0)*xc(2,1))*(-(normal[1]*(xc(0,0) -
    xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) +
    3*(-2*normal[3]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(-xc(0,2) +
    xc(1,2)))*(-2*normal[3]*(-xc(1,1) + xc(2,1)) + 2*normal[2]*(-xc(1,2) +
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(-xc(1,1) + xc(2,1))
    + 2*normal[2]*(-xc(1,2) + xc(2,2)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) +
    (xc(0,0) - xc(1,0))*(xc(0,2) - xc(3,2))) + 4*normpowfour*(xc(1,2) -
    xc(3,2)) - 2*normsquare*(-2*normal[3]*(xc(0,0) - xc(1,0)) -
    2*normal[1]*(-xc(0,2) + xc(1,2)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1)
    - xc(2,1)*xc(3,2) + xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(7,1)=
    (-2*normsquare*(2*(xc(0,0) - xc(1,0))*(xc(1,0) - xc(2,0)) + 2*(xc(0,2) -
    xc(1,2))*(xc(1,2) - xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) +
    3*(-2*normal[3]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(-xc(0,2) +
    xc(1,2)))*(-2*normal[3]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(-xc(1,2) +
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(xc(1,0) - xc(2,0))
    - 2*normal[1]*(-xc(1,2) + xc(2,2)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) +
    (xc(0,0) - xc(1,0))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(xc(0,0) -
    xc(1,0)) - 2*normal[1]*(-xc(0,2) + xc(1,2)))*(-(xc(2,2)*xc(3,0)) +
    xc(1,2)*(-xc(2,0) + xc(3,0)) + xc(1,0)*(xc(2,2) - xc(3,2)) +
    xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(7,2)=
    (4*normpowfour*(-xc(1,0) + xc(3,0)) -
    2*normsquare*(-2*normal[3]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(-xc(0,2) +
    xc(1,2)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0) - xc(2,0)*xc(3,1) +
    xc(1,0)*(-xc(2,1) + xc(3,1))) + 3*(-2*normal[3]*(xc(0,0) - xc(1,0)) -
    2*normal[1]*(-xc(0,2) + xc(1,2)))*(2*normal[2]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(xc(1,1) - xc(2,1)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    4*normsquare*(xc(1,1)*xc(1,2) - 2*xc(0,2)*(xc(1,1) - xc(2,1)) -
    2*xc(1,2)*xc(2,1) + xc(0,1)*(xc(1,2) - xc(2,2)) +
    xc(1,1)*xc(2,2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) +
    xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal[2]*(xc(1,0) -
    xc(2,0)) - 2*normal[1]*(xc(1,1) - xc(2,1)))*((-xc(0,2) + xc(1,2))*(xc(0,0) -
    xc(3,0)) + (xc(0,0) - xc(1,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(7,3)=
    (-4*normsquare*(-(xc(1,1)*xc(2,0)) + xc(0,1)*(-2*xc(1,0) + xc(2,0)) +
    xc(0,0)*(xc(0,1) + xc(1,1) - 2*xc(2,1)) +
    2*xc(1,0)*xc(2,1))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) +
    xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) + 3*(-2*normal[3]*(xc(0,0) - xc(1,0))
    - 2*normal[1]*(-xc(0,2) + xc(1,2)))*(-2*normal[3]*(xc(0,1) - xc(2,1)) +
    2*normal[2]*(xc(0,2) - xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal[3]*(xc(0,1) - xc(2,1)) + 2*normal[2]*(xc(0,2) -
    xc(2,2)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(xc(0,0) - xc(1,0)) -
    2*normal[1]*(-xc(0,2) + xc(1,2)))*((xc(0,2) - xc(2,2))*(-xc(0,1) + xc(3,1)) +
    (xc(0,1) - xc(2,1))*(xc(0,2) - xc(3,2))) + 4*normpowfour*(-xc(0,2) +
    xc(3,2)))/(4.*normpowfive);

  elematrix(7,4)=
    (3*(-2*normal[3]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(-xc(0,2) +
    xc(1,2)))*(-2*normal[3]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(xc(0,2) -
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*(xc(0,0) - xc(1,0))*(-xc(0,0) +
    xc(2,0)) + 2*(-xc(0,2) + xc(1,2))*(xc(0,2) - xc(2,2)))*(-(normal[1]*(xc(0,0) -
    xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(-2*normal[3]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(xc(0,2) -
    xc(2,2)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(xc(0,0) - xc(1,0)) -
    2*normal[1]*(-xc(0,2) + xc(1,2)))*((xc(0,2) - xc(2,2))*(xc(0,0) - xc(3,0)) +
    (-xc(0,0) + xc(2,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(7,5)=
    (4*normpowfour*(xc(0,0) - xc(3,0)) - 2*normsquare*(-2*normal[3]*(xc(0,0)
    - xc(1,0)) - 2*normal[1]*(-xc(0,2) + xc(1,2)))*((-xc(0,1) + xc(2,1))*(xc(0,0) -
    xc(3,0)) + (xc(0,0) - xc(2,0))*(xc(0,1) - xc(3,1))) + 3*(-2*normal[3]*(xc(0,0) -
    xc(1,0)) - 2*normal[1]*(-xc(0,2) + xc(1,2)))*(2*normal[2]*(-xc(0,0) + xc(2,0)) -
    2*normal[1]*(-xc(0,1) + xc(2,1)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    4*normsquare*(xc(0,2)*(xc(1,1) - 2*xc(2,1)) + 2*xc(1,2)*xc(2,1) -
    xc(1,1)*xc(2,2) + xc(0,1)*(xc(0,2) - 2*xc(1,2) + xc(2,2)))*(-(normal[1]*(xc(0,0)
    - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal[2]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(-xc(0,1) +
    xc(2,1)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(7,6)=
    (-4*normsquare*(xc(0,0) - xc(1,0))*(-xc(0,1) + xc(1,1))*(-(normal[1]*(xc(0,0) -
    xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) +
    3*(-2*normal[3]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(-xc(0,2) +
    xc(1,2)))*(-2*normal[3]*(-xc(0,1) + xc(1,1)) + 2*normal[2]*(-xc(0,2) +
    xc(1,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(-xc(0,1) + xc(1,1))
    + 2*normal[2]*(-xc(0,2) + xc(1,2)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) +
    (xc(0,0) - xc(1,0))*(xc(0,2) - xc(3,2))) - 2*normsquare*(-2*normal[3]*(xc(0,0) -
    xc(1,0)) - 2*normal[1]*(-xc(0,2) + xc(1,2)))*((xc(0,2) - xc(1,2))*(xc(0,1) -
    xc(3,1)) + (-xc(0,1) + xc(1,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(7,7)=
    (-4*normsquare*(pow(xc(0,0) - xc(1,0),2) + pow(xc(0,2) -
    xc(1,2),2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) + 3*pow(-2*normal[3]*(xc(0,0) - xc(1,0)) -
    2*normal[1]*(-xc(0,2) + xc(1,2)),2)*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-2*normal[3]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(-xc(0,2) +
    xc(1,2)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(7,8)=
    (-2*normsquare*(-2*normal[3]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(-xc(0,2) +
    xc(1,2)))*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(-xc(0,1) + xc(3,1))) - 4*normsquare*(xc(0,1) - xc(1,1))*(-xc(0,2) +
    xc(1,2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) + 3*(2*normal[2]*(xc(0,0) - xc(1,0)) -
    2*normal[1]*(xc(0,1) - xc(1,1)))*(-2*normal[3]*(xc(0,0) - xc(1,0)) -
    2*normal[1]*(-xc(0,2) + xc(1,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal[2]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(xc(0,1) -
    xc(1,1)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(7,9)=
    (normal[2]*(xc(0,0)*pow(xc(1,1),2) + xc(0,0)*pow(xc(1,2),2) +
    pow(xc(0,1),2)*(xc(1,0) - xc(2,0)) + pow(xc(0,2),2)*(xc(1,0) - xc(2,0)) -
    pow(xc(1,1),2)*xc(2,0) - pow(xc(1,2),2)*xc(2,0) - xc(0,0)*xc(1,1)*xc(2,1) +
    xc(1,0)*xc(1,1)*xc(2,1) - xc(0,1)*(-2*xc(1,1)*xc(2,0) + xc(0,0)*(xc(1,1) -
    xc(2,1)) + xc(1,0)*(xc(1,1) + xc(2,1))) - xc(0,0)*xc(1,2)*xc(2,2) +
    xc(1,0)*xc(1,2)*xc(2,2) - xc(0,2)*(-2*xc(1,2)*xc(2,0) + xc(0,0)*(xc(1,2) -
    xc(2,2)) + xc(1,0)*(xc(1,2) + xc(2,2)))))/normcube;

  elematrix(7,10)=
    -(normal[2]*(-2*normal[3]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(-xc(0,2) +
    xc(1,2))))/(2.*normcube);

  elematrix(7,11)=
    (normal[2]*(pow(xc(0,1),2)*xc(1,2) - xc(0,1)*xc(1,1)*xc(1,2) +
    xc(1,0)*xc(1,2)*xc(2,0) + xc(0,2)*(pow(xc(1,0),2) - xc(1,0)*xc(2,0) - (xc(0,1)
    - xc(1,1))*(xc(1,1) - xc(2,1))) - xc(0,1)*xc(1,2)*xc(2,1) +
    xc(1,1)*xc(1,2)*xc(2,1) - xc(0,0)*(xc(0,2)*(xc(1,0) - xc(2,0)) + xc(1,2)*xc(2,0)
    + xc(1,0)*(xc(1,2) - 2*xc(2,2))) + pow(xc(0,0),2)*(xc(1,2) - xc(2,2)) -
    pow(xc(0,1),2)*xc(2,2) - pow(xc(1,0),2)*xc(2,2) + 2*xc(0,1)*xc(1,1)*xc(2,2)
    - pow(xc(1,1),2)*xc(2,2)))/normcube;

  elematrix(8,0)=
    (4*normpowfour*(-xc(1,1) + xc(3,1)) -
    2*normsquare*(-2*normal[3]*(-xc(1,1) + xc(2,1)) + 2*normal[2]*(-xc(1,2) +
    xc(2,2)))*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(-xc(0,1) + xc(3,1))) - 4*normsquare*(xc(1,0)*xc(1,2) +
    xc(0,2)*(xc(1,0) - xc(2,0)) + xc(1,2)*xc(2,0) - 2*xc(0,0)*(xc(1,2) - xc(2,2)) -
    2*xc(1,0)*xc(2,2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) +
    xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) + 3*(2*normal[2]*(xc(0,0) - xc(1,0)) -
    2*normal[1]*(xc(0,1) - xc(1,1)))*(-2*normal[3]*(-xc(1,1) + xc(2,1)) +
    2*normal[2]*(-xc(1,2) + xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal[2]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(xc(0,1) -
    xc(1,1)))*(xc(1,2)*(xc(2,1) - xc(3,1)) + xc(2,2)*xc(3,1) - xc(2,1)*xc(3,2) +
    xc(1,1)*(-xc(2,2) + xc(3,2))))/(4.*normpowfive);

  elematrix(8,1)=
    (4*normpowfour*(xc(1,0) - xc(3,0)) - 2*normsquare*(-2*normal[3]*(xc(1,0)
    - xc(2,0)) - 2*normal[1]*(-xc(1,2) + xc(2,2)))*((xc(0,1) - xc(1,1))*(xc(0,0) -
    xc(3,0)) + (xc(0,0) - xc(1,0))*(-xc(0,1) + xc(3,1))) -
    4*normsquare*(xc(1,1)*xc(1,2) + xc(0,2)*(xc(1,1) - xc(2,1)) + xc(1,2)*xc(2,1) -
    2*xc(0,1)*(xc(1,2) - xc(2,2)) - 2*xc(1,1)*xc(2,2))*(-(normal[1]*(xc(0,0) -
    xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) +
    3*(2*normal[2]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(xc(0,1) -
    xc(1,1)))*(-2*normal[3]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(-xc(1,2) +
    xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal[2]*(xc(0,0) - xc(1,0)) -
    2*normal[1]*(xc(0,1) - xc(1,1)))*(-(xc(2,2)*xc(3,0)) + xc(1,2)*(-xc(2,0) +
    xc(3,0)) + xc(1,0)*(xc(2,2) - xc(3,2)) +
    xc(2,0)*xc(3,2)))/(4.*normpowfive);

  elematrix(8,2)=
    (-2*normsquare*(2*normal[2]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(xc(1,1) -
    xc(2,1)))*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(-xc(0,1) + xc(3,1))) - 2*normsquare*(2*normal[2]*(xc(0,0) - xc(1,0)) -
    2*normal[1]*(xc(0,1) - xc(1,1)))*(xc(1,1)*(xc(2,0) - xc(3,0)) + xc(2,1)*xc(3,0)
    - xc(2,0)*xc(3,1) + xc(1,0)*(-xc(2,1) + xc(3,1))) + 3*(2*normal[2]*(xc(0,0) -
    xc(1,0)) - 2*normal[1]*(xc(0,1) - xc(1,1)))*(2*normal[2]*(xc(1,0) - xc(2,0)) -
    2*normal[1]*(xc(1,1) - xc(2,1)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*(xc(0,0) - xc(1,0))*(xc(1,0) - xc(2,0)) + 2*(xc(0,1) -
    xc(1,1))*(xc(1,1) - xc(2,1)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) -
    xc(3,2))))/(4.*normpowfive);

  elematrix(8,3)=
    (4*normpowfour*(xc(0,1) - xc(3,1)) - 2*normsquare*(-2*normal[3]*(xc(0,1)
    - xc(2,1)) + 2*normal[2]*(xc(0,2) - xc(2,2)))*((xc(0,1) - xc(1,1))*(xc(0,0) -
    xc(3,0)) + (xc(0,0) - xc(1,0))*(-xc(0,1) + xc(3,1))) + 3*(2*normal[2]*(xc(0,0) -
    xc(1,0)) - 2*normal[1]*(xc(0,1) - xc(1,1)))*(-2*normal[3]*(xc(0,1) - xc(2,1)) +
    2*normal[2]*(xc(0,2) - xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-(xc(1,2)*xc(2,0)) + xc(0,2)*(-2*xc(1,0) + xc(2,0)) +
    xc(0,0)*(xc(0,2) + xc(1,2) - 2*xc(2,2)) +
    2*xc(1,0)*xc(2,2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) +
    xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal[2]*(xc(0,0) -
    xc(1,0)) - 2*normal[1]*(xc(0,1) - xc(1,1)))*((xc(0,2) - xc(2,2))*(-xc(0,1) +
    xc(3,1)) + (xc(0,1) - xc(2,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(8,4)=
    (4*normpowfour*(-xc(0,0) + xc(3,0)) -
    2*normsquare*(-2*normal[3]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(xc(0,2) -
    xc(2,2)))*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(-xc(0,1) + xc(3,1))) + 3*(2*normal[2]*(xc(0,0) - xc(1,0)) -
    2*normal[1]*(xc(0,1) - xc(1,1)))*(-2*normal[3]*(-xc(0,0) + xc(2,0)) -
    2*normal[1]*(xc(0,2) - xc(2,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    4*normsquare*(-(xc(1,2)*xc(2,1)) + xc(0,2)*(-2*xc(1,1) + xc(2,1)) +
    xc(0,1)*(xc(0,2) + xc(1,2) - 2*xc(2,2)) +
    2*xc(1,1)*xc(2,2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) +
    xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) - 2*normsquare*(2*normal[2]*(xc(0,0) -
    xc(1,0)) - 2*normal[1]*(xc(0,1) - xc(1,1)))*((xc(0,2) - xc(2,2))*(xc(0,0) -
    xc(3,0)) + (-xc(0,0) + xc(2,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(8,5)=
    (-2*normsquare*(2*normal[2]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(xc(0,1) -
    xc(1,1)))*((-xc(0,1) + xc(2,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(2,0))*(xc(0,1) - xc(3,1))) - 2*normsquare*(2*normal[2]*(-xc(0,0) + xc(2,0)) -
    2*normal[1]*(-xc(0,1) + xc(2,1)))*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) +
    (xc(0,0) - xc(1,0))*(-xc(0,1) + xc(3,1))) + 3*(2*normal[2]*(xc(0,0) - xc(1,0)) -
    2*normal[1]*(xc(0,1) - xc(1,1)))*(2*normal[2]*(-xc(0,0) + xc(2,0)) -
    2*normal[1]*(-xc(0,1) + xc(2,1)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*(xc(0,0) - xc(1,0))*(-xc(0,0) + xc(2,0)) + 2*(xc(0,1) -
    xc(1,1))*(-xc(0,1) + xc(2,1)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) -
    xc(3,2))))/(4.*normpowfive);

  elematrix(8,6)=
    (-2*normsquare*(-2*normal[3]*(-xc(0,1) + xc(1,1)) + 2*normal[2]*(-xc(0,2) +
    xc(1,2)))*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(-xc(0,1) + xc(3,1))) - 4*normsquare*(xc(0,0) - xc(1,0))*(-xc(0,2) +
    xc(1,2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) + 3*(2*normal[2]*(xc(0,0) - xc(1,0)) -
    2*normal[1]*(xc(0,1) - xc(1,1)))*(-2*normal[3]*(-xc(0,1) + xc(1,1)) +
    2*normal[2]*(-xc(0,2) + xc(1,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal[2]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(xc(0,1) -
    xc(1,1)))*((xc(0,2) - xc(1,2))*(xc(0,1) - xc(3,1)) + (-xc(0,1) +
    xc(1,1))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(8,7)=
    (-2*normsquare*(-2*normal[3]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(-xc(0,2) +
    xc(1,2)))*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(-xc(0,1) + xc(3,1))) - 4*normsquare*(xc(0,1) - xc(1,1))*(-xc(0,2) +
    xc(1,2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))) + 3*(2*normal[2]*(xc(0,0) - xc(1,0)) -
    2*normal[1]*(xc(0,1) - xc(1,1)))*(-2*normal[3]*(xc(0,0) - xc(1,0)) -
    2*normal[1]*(-xc(0,2) + xc(1,2)))*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    2*normsquare*(2*normal[2]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(xc(0,1) -
    xc(1,1)))*((-xc(0,2) + xc(1,2))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(8,8)=
    (-4*normsquare*(2*normal[2]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(xc(0,1) -
    xc(1,1)))*((xc(0,1) - xc(1,1))*(xc(0,0) - xc(3,0)) + (xc(0,0) -
    xc(1,0))*(-xc(0,1) + xc(3,1))) + 3*pow(2*normal[2]*(xc(0,0) - xc(1,0)) -
    2*normal[1]*(xc(0,1) - xc(1,1)),2)*(-(normal[1]*(xc(0,0) - xc(3,0))) +
    normal[2]*(-xc(0,1) + xc(3,1)) - normal[3]*(xc(0,2) - xc(3,2))) -
    4*normsquare*(pow(xc(0,0) - xc(1,0),2) + pow(xc(0,1) -
    xc(1,1),2))*(-(normal[1]*(xc(0,0) - xc(3,0))) + normal[2]*(-xc(0,1) + xc(3,1)) -
    normal[3]*(xc(0,2) - xc(3,2))))/(4.*normpowfive);

  elematrix(8,9)=
    -((normal[3]*(-(xc(0,0)*pow(xc(1,1),2)) - xc(0,0)*pow(xc(1,2),2) +
    pow(xc(1,1),2)*xc(2,0) + pow(xc(1,2),2)*xc(2,0) + pow(xc(0,1),2)*(-xc(1,0)
    + xc(2,0)) + pow(xc(0,2),2)*(-xc(1,0) + xc(2,0)) + xc(0,0)*xc(1,1)*xc(2,1) -
    xc(1,0)*xc(1,1)*xc(2,1) + xc(0,1)*(-2*xc(1,1)*xc(2,0) + xc(0,0)*(xc(1,1) -
    xc(2,1)) + xc(1,0)*(xc(1,1) + xc(2,1))) + xc(0,0)*xc(1,2)*xc(2,2) -
    xc(1,0)*xc(1,2)*xc(2,2) + xc(0,2)*(-2*xc(1,2)*xc(2,0) + xc(0,0)*(xc(1,2) -
    xc(2,2)) + xc(1,0)*(xc(1,2) + xc(2,2)))))/normcube);

  elematrix(8,10)=
    (normal[3]*(pow(xc(0,2),2)*xc(1,1) - xc(0,2)*xc(1,1)*xc(1,2) +
    xc(1,0)*xc(1,1)*xc(2,0) - xc(0,0)*(xc(0,1)*(xc(1,0) - xc(2,0)) + xc(1,1)*xc(2,0)
    + xc(1,0)*(xc(1,1) - 2*xc(2,1))) + pow(xc(0,0),2)*(xc(1,1) - xc(2,1)) -
    pow(xc(0,2),2)*xc(2,1) - pow(xc(1,0),2)*xc(2,1) + 2*xc(0,2)*xc(1,2)*xc(2,1)
    - pow(xc(1,2),2)*xc(2,1) + xc(0,1)*(pow(xc(1,0),2) - xc(1,0)*xc(2,0) -
    (xc(0,2) - xc(1,2))*(xc(1,2) - xc(2,2))) - xc(0,2)*xc(1,1)*xc(2,2) +
    xc(1,1)*xc(1,2)*xc(2,2)))/normcube;

  elematrix(8,11)=
    -((2*normal[2]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(xc(0,1) -
    xc(1,1)))*(-(xc(1,1)*xc(2,0)) + xc(0,1)*(-xc(1,0) + xc(2,0)) + xc(0,0)*(xc(1,1)
    - xc(2,1)) + xc(1,0)*xc(2,1)))/(2.*normcube);

  elematrix(9,0)=
    -((-(xc(1,2)*xc(2,1)) + xc(0,2)*(-xc(1,1) + xc(2,1)) + xc(0,1)*(xc(1,2) -
    xc(2,2)) + xc(1,1)*xc(2,2))*(-2*normal[3]*(-xc(1,1) + xc(2,1)) +
    2*normal[2]*(-xc(1,2) + xc(2,2))))/(2.*normcube);

  elematrix(9,1)=
    (normal[2]*(xc(0,2)*xc(1,0)*xc(1,2) + pow(xc(1,1),2)*xc(2,0) -
    xc(0,2)*xc(1,2)*xc(2,0) + pow(xc(1,2),2)*xc(2,0) + xc(0,1)*(xc(1,0) -
    xc(2,0))*(xc(1,1) - xc(2,1)) - xc(1,0)*xc(1,1)*xc(2,1) - xc(1,1)*xc(2,0)*xc(2,1)
    + xc(1,0)*pow(xc(2,1),2) - xc(0,2)*xc(1,0)*xc(2,2) - xc(1,0)*xc(1,2)*xc(2,2) +
    xc(0,2)*xc(2,0)*xc(2,2) - xc(1,2)*xc(2,0)*xc(2,2) + xc(1,0)*pow(xc(2,2),2) -
    xc(0,0)*(pow(xc(1,1),2) + pow(xc(1,2),2) - 2*xc(1,1)*xc(2,1) +
    pow(xc(2,1),2) - 2*xc(1,2)*xc(2,2) + pow(xc(2,2),2))))/normcube;

  elematrix(9,2)=
    (normal[3]*(xc(0,2)*xc(1,0)*xc(1,2) + pow(xc(1,1),2)*xc(2,0) -
    xc(0,2)*xc(1,2)*xc(2,0) + pow(xc(1,2),2)*xc(2,0) + xc(0,1)*(xc(1,0) -
    xc(2,0))*(xc(1,1) - xc(2,1)) - xc(1,0)*xc(1,1)*xc(2,1) - xc(1,1)*xc(2,0)*xc(2,1)
    + xc(1,0)*pow(xc(2,1),2) - xc(0,2)*xc(1,0)*xc(2,2) - xc(1,0)*xc(1,2)*xc(2,2) +
    xc(0,2)*xc(2,0)*xc(2,2) - xc(1,2)*xc(2,0)*xc(2,2) + xc(1,0)*pow(xc(2,2),2) -
    xc(0,0)*(pow(xc(1,1),2) + pow(xc(1,2),2) - 2*xc(1,1)*xc(2,1) +
    pow(xc(2,1),2) - 2*xc(1,2)*xc(2,2) + pow(xc(2,2),2))))/normcube;

  elematrix(9,3)=
    -((-2*normal[3]*(xc(0,1) - xc(2,1)) + 2*normal[2]*(xc(0,2) -
    xc(2,2)))*(-(xc(1,2)*xc(2,1)) + xc(0,2)*(-xc(1,1) + xc(2,1)) + xc(0,1)*(xc(1,2)
    - xc(2,2)) + xc(1,1)*xc(2,2)))/(2.*normcube);

  elematrix(9,4)=
    -((normal[2]*(pow(xc(0,1),2)*(xc(1,0) - xc(2,0)) + pow(xc(0,2),2)*(xc(1,0) -
    xc(2,0)) + xc(0,0)*xc(1,1)*xc(2,1) - xc(1,1)*xc(2,0)*xc(2,1) -
    xc(0,0)*pow(xc(2,1),2) + xc(1,0)*pow(xc(2,1),2) + xc(0,1)*(xc(1,1)*xc(2,0) +
    (-2*xc(1,0) + xc(2,0))*xc(2,1) + xc(0,0)*(-xc(1,1) + xc(2,1))) +
    xc(0,0)*xc(1,2)*xc(2,2) - xc(1,2)*xc(2,0)*xc(2,2) - xc(0,0)*pow(xc(2,2),2) +
    xc(1,0)*pow(xc(2,2),2) + xc(0,2)*(xc(1,2)*xc(2,0) + (-2*xc(1,0) +
    xc(2,0))*xc(2,2) + xc(0,0)*(-xc(1,2) + xc(2,2)))))/normcube);

  elematrix(9,5)=
    -((normal[3]*(pow(xc(0,1),2)*(xc(1,0) - xc(2,0)) + pow(xc(0,2),2)*(xc(1,0) -
    xc(2,0)) + xc(0,0)*xc(1,1)*xc(2,1) - xc(1,1)*xc(2,0)*xc(2,1) -
    xc(0,0)*pow(xc(2,1),2) + xc(1,0)*pow(xc(2,1),2) + xc(0,1)*(xc(1,1)*xc(2,0) +
    (-2*xc(1,0) + xc(2,0))*xc(2,1) + xc(0,0)*(-xc(1,1) + xc(2,1))) +
    xc(0,0)*xc(1,2)*xc(2,2) - xc(1,2)*xc(2,0)*xc(2,2) - xc(0,0)*pow(xc(2,2),2) +
    xc(1,0)*pow(xc(2,2),2) + xc(0,2)*(xc(1,2)*xc(2,0) + (-2*xc(1,0) +
    xc(2,0))*xc(2,2) + xc(0,0)*(-xc(1,2) + xc(2,2)))))/normcube);

  elematrix(9,6)=
    -((-2*normal[3]*(-xc(0,1) + xc(1,1)) + 2*normal[2]*(-xc(0,2) +
    xc(1,2)))*(-(xc(1,2)*xc(2,1)) + xc(0,2)*(-xc(1,1) + xc(2,1)) + xc(0,1)*(xc(1,2)
    - xc(2,2)) + xc(1,1)*xc(2,2)))/(2.*normcube);

  elematrix(9,7)=
    (normal[2]*(xc(0,0)*pow(xc(1,1),2) + xc(0,0)*pow(xc(1,2),2) +
    pow(xc(0,1),2)*(xc(1,0) - xc(2,0)) + pow(xc(0,2),2)*(xc(1,0) - xc(2,0)) -
    pow(xc(1,1),2)*xc(2,0) - pow(xc(1,2),2)*xc(2,0) - xc(0,0)*xc(1,1)*xc(2,1) +
    xc(1,0)*xc(1,1)*xc(2,1) - xc(0,1)*(-2*xc(1,1)*xc(2,0) + xc(0,0)*(xc(1,1) -
    xc(2,1)) + xc(1,0)*(xc(1,1) + xc(2,1))) - xc(0,0)*xc(1,2)*xc(2,2) +
    xc(1,0)*xc(1,2)*xc(2,2) - xc(0,2)*(-2*xc(1,2)*xc(2,0) + xc(0,0)*(xc(1,2) -
    xc(2,2)) + xc(1,0)*(xc(1,2) + xc(2,2)))))/normcube;

  elematrix(9,8)=
    -((normal[3]*(-(xc(0,0)*pow(xc(1,1),2)) - xc(0,0)*pow(xc(1,2),2) +
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
    (normal[1]*(xc(0,2)*xc(1,1)*xc(1,2) - xc(1,0)*xc(1,1)*xc(2,0) +
    xc(1,1)*pow(xc(2,0),2) + xc(0,0)*(xc(1,0) - xc(2,0))*(xc(1,1) - xc(2,1)) +
    pow(xc(1,0),2)*xc(2,1) - xc(0,2)*xc(1,2)*xc(2,1) + pow(xc(1,2),2)*xc(2,1) -
    xc(1,0)*xc(2,0)*xc(2,1) - xc(0,2)*xc(1,1)*xc(2,2) - xc(1,1)*xc(1,2)*xc(2,2) +
    xc(0,2)*xc(2,1)*xc(2,2) - xc(1,2)*xc(2,1)*xc(2,2) + xc(1,1)*pow(xc(2,2),2) -
    xc(0,1)*(pow(xc(1,0),2) + pow(xc(1,2),2) - 2*xc(1,0)*xc(2,0) +
    pow(xc(2,0),2) - 2*xc(1,2)*xc(2,2) + pow(xc(2,2),2))))/normcube;

  elematrix(10,1)=
    -(normal[2]*(-2*normal[3]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(-xc(1,2) +
    xc(2,2))))/(2.*normcube);

  elematrix(10,2)=
    -((normal[3]*(-(xc(0,2)*xc(1,1)*xc(1,2)) + xc(1,0)*xc(1,1)*xc(2,0) -
    xc(1,1)*pow(xc(2,0),2) - xc(0,0)*(xc(1,0) - xc(2,0))*(xc(1,1) - xc(2,1)) -
    pow(xc(1,0),2)*xc(2,1) + xc(0,2)*xc(1,2)*xc(2,1) - pow(xc(1,2),2)*xc(2,1) +
    xc(1,0)*xc(2,0)*xc(2,1) + xc(0,2)*xc(1,1)*xc(2,2) + xc(1,1)*xc(1,2)*xc(2,2) -
    xc(0,2)*xc(2,1)*xc(2,2) + xc(1,2)*xc(2,1)*xc(2,2) - xc(1,1)*pow(xc(2,2),2) +
    xc(0,1)*(pow(xc(1,0),2) + pow(xc(1,2),2) - 2*xc(1,0)*xc(2,0) +
    pow(xc(2,0),2) - 2*xc(1,2)*xc(2,2) +pow(xc(2,2),2))))/normcube);

  elematrix(10,3)=
    -((normal[1]*(xc(0,1)*xc(1,0)*xc(2,0) - xc(0,1)*pow(xc(2,0),2) +
    xc(1,1)*pow(xc(2,0),2) + pow(xc(0,0),2)*(xc(1,1) - xc(2,1)) +
    pow(xc(0,2),2)*(xc(1,1) - xc(2,1)) - xc(1,0)*xc(2,0)*xc(2,1) +
    xc(0,0)*(-2*xc(1,1)*xc(2,0) + xc(0,1)*(-xc(1,0) + xc(2,0)) + (xc(1,0) +
    xc(2,0))*xc(2,1)) + xc(0,1)*xc(1,2)*xc(2,2) - xc(1,2)*xc(2,1)*xc(2,2) -
    xc(0,1)*pow(xc(2,2),2) + xc(1,1)*pow(xc(2,2),2) + xc(0,2)*(xc(1,2)*xc(2,1) +
    (-2*xc(1,1) + xc(2,1))*xc(2,2) + xc(0,1)*(-xc(1,2) +xc(2,2)))))/normcube);

  elematrix(10,4)=
    -(normal[2]*(-2*normal[3]*(-xc(0,0) + xc(2,0)) - 2*normal[1]*(xc(0,2) -
    xc(2,2))))/(2.*normcube);

  elematrix(10,5)=
    -((normal[3]*(xc(0,1)*xc(1,0)*xc(2,0) - xc(0,1)*pow(xc(2,0),2) +
    xc(1,1)*pow(xc(2,0),2) + pow(xc(0,0),2)*(xc(1,1) - xc(2,1)) +
    pow(xc(0,2),2)*(xc(1,1) - xc(2,1)) - xc(1,0)*xc(2,0)*xc(2,1) +
    xc(0,0)*(-2*xc(1,1)*xc(2,0) + xc(0,1)*(-xc(1,0) + xc(2,0)) + (xc(1,0) +
    xc(2,0))*xc(2,1)) + xc(0,1)*xc(1,2)*xc(2,2) - xc(1,2)*xc(2,1)*xc(2,2) -
    xc(0,1)*pow(xc(2,2),2) + xc(1,1)*pow(xc(2,2),2) + xc(0,2)*(xc(1,2)*xc(2,1) +
    (-2*xc(1,1) + xc(2,1))*xc(2,2) + xc(0,1)*(-xc(1,2) +xc(2,2)))))/normcube);

  elematrix(10,6)=
    (normal[1]*(pow(xc(0,2),2)*xc(1,1) - xc(0,2)*xc(1,1)*xc(1,2) +
    xc(1,0)*xc(1,1)*xc(2,0) - xc(0,0)*(xc(0,1)*(xc(1,0) - xc(2,0)) + xc(1,1)*xc(2,0)
    + xc(1,0)*(xc(1,1) - 2*xc(2,1))) + pow(xc(0,0),2)*(xc(1,1) - xc(2,1)) -
    pow(xc(0,2),2)*xc(2,1) - pow(xc(1,0),2)*xc(2,1) + 2*xc(0,2)*xc(1,2)*xc(2,1)
    - pow(xc(1,2),2)*xc(2,1) + xc(0,1)*(pow(xc(1,0),2) - xc(1,0)*xc(2,0) -
    (xc(0,2) - xc(1,2))*(xc(1,2) - xc(2,2))) - xc(0,2)*xc(1,1)*xc(2,2) +
    xc(1,1)*xc(1,2)*xc(2,2)))/normcube;

  elematrix(10,7)=
    -(normal[2]*(-2*normal[3]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(-xc(0,2) +
    xc(1,2))))/(2.*normcube);

  elematrix(10,8)=
    (normal[3]*(pow(xc(0,2),2)*xc(1,1) - xc(0,2)*xc(1,1)*xc(1,2) +
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
    -((normal[1]*(-(xc(0,1)*xc(1,1)*xc(1,2)) + xc(1,0)*xc(1,2)*xc(2,0) -
    xc(1,2)*pow(xc(2,0),2) + xc(0,1)*xc(1,2)*xc(2,1) + xc(1,1)*xc(1,2)*xc(2,1) -
    xc(1,2)*pow(xc(2,1),2) + xc(0,2)*(pow(xc(1,0),2) + pow(xc(1,1),2) -
    2*xc(1,0)*xc(2,0) + pow(xc(2,0),2) - 2*xc(1,1)*xc(2,1) + pow(xc(2,1),2)) -
    xc(0,0)*(xc(1,0) - xc(2,0))*(xc(1,2) - xc(2,2)) - pow(xc(1,0),2)*xc(2,2) +
    xc(0,1)*xc(1,1)*xc(2,2) - pow(xc(1,1),2)*xc(2,2) + xc(1,0)*xc(2,0)*xc(2,2) -
    xc(0,1)*xc(2,1)*xc(2,2) + xc(1,1)*xc(2,1)*xc(2,2)))/normcube);

  elematrix(11,1)=
    -((normal[2]*(-(xc(0,1)*xc(1,1)*xc(1,2)) + xc(1,0)*xc(1,2)*xc(2,0) -
    xc(1,2)*pow(xc(2,0),2) + xc(0,1)*xc(1,2)*xc(2,1) + xc(1,1)*xc(1,2)*xc(2,1) -
    xc(1,2)*pow(xc(2,1),2) + xc(0,2)*(pow(xc(1,0),2) + pow(xc(1,1),2) -
    2*xc(1,0)*xc(2,0) + pow(xc(2,0),2) - 2*xc(1,1)*xc(2,1) + pow(xc(2,1),2)) -
    xc(0,0)*(xc(1,0) - xc(2,0))*(xc(1,2) - xc(2,2)) - pow(xc(1,0),2)*xc(2,2) +
    xc(0,1)*xc(1,1)*xc(2,2) - pow(xc(1,1),2)*xc(2,2) + xc(1,0)*xc(2,0)*xc(2,2) -
    xc(0,1)*xc(2,1)*xc(2,2) + xc(1,1)*xc(2,1)*xc(2,2)))/normcube);

  elematrix(11,2)=
    -((2*normal[2]*(xc(1,0) - xc(2,0)) - 2*normal[1]*(xc(1,1) -
    xc(2,1)))*(-(xc(1,1)*xc(2,0)) + xc(0,1)*(-xc(1,0) + xc(2,0)) + xc(0,0)*(xc(1,1)
    - xc(2,1)) + xc(1,0)*xc(2,1)))/(2.*normcube);

  elematrix(11,3)=
    -((normal[1]*(xc(0,2)*xc(1,0)*xc(2,0) - xc(0,2)*pow(xc(2,0),2) +
    xc(1,2)*pow(xc(2,0),2) + xc(0,2)*xc(1,1)*xc(2,1) - xc(0,2)*pow(xc(2,1),2) +
    xc(1,2)*pow(xc(2,1),2) + pow(xc(0,0),2)*(xc(1,2) - xc(2,2)) +
    pow(xc(0,1),2)*(xc(1,2) - xc(2,2)) - xc(1,0)*xc(2,0)*xc(2,2) -
    xc(1,1)*xc(2,1)*xc(2,2) + xc(0,0)*(-2*xc(1,2)*xc(2,0) + xc(0,2)*(-xc(1,0) +
    xc(2,0)) + (xc(1,0) + xc(2,0))*xc(2,2)) + xc(0,1)*(-2*xc(1,2)*xc(2,1) +
    xc(0,2)*(-xc(1,1) + xc(2,1)) + (xc(1,1) +xc(2,1))*xc(2,2))))/normcube);

  elematrix(11,4)=
    -((normal[2]*(xc(0,2)*xc(1,0)*xc(2,0) - xc(0,2)*pow(xc(2,0),2) +
    xc(1,2)*pow(xc(2,0),2) + xc(0,2)*xc(1,1)*xc(2,1) - xc(0,2)*pow(xc(2,1),2) +
    xc(1,2)*pow(xc(2,1),2) + pow(xc(0,0),2)*(xc(1,2) - xc(2,2)) +
    pow(xc(0,1),2)*(xc(1,2) - xc(2,2)) - xc(1,0)*xc(2,0)*xc(2,2) -
    xc(1,1)*xc(2,1)*xc(2,2) + xc(0,0)*(-2*xc(1,2)*xc(2,0) + xc(0,2)*(-xc(1,0) +
    xc(2,0)) + (xc(1,0) + xc(2,0))*xc(2,2)) + xc(0,1)*(-2*xc(1,2)*xc(2,1) +
    xc(0,2)*(-xc(1,1) + xc(2,1)) + (xc(1,1) +xc(2,1))*xc(2,2))))/normcube);

  elematrix(11,5)=
    -((-(xc(1,1)*xc(2,0)) + xc(0,1)*(-xc(1,0) + xc(2,0)) + xc(0,0)*(xc(1,1) -
    xc(2,1)) + xc(1,0)*xc(2,1))*(2*normal[2]*(-xc(0,0) + xc(2,0)) -
    2*normal[1]*(-xc(0,1) + xc(2,1))))/(2.*normcube);

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
    (normal[2]*(pow(xc(0,1),2)*xc(1,2) - xc(0,1)*xc(1,1)*xc(1,2) +
    xc(1,0)*xc(1,2)*xc(2,0) + xc(0,2)*(pow(xc(1,0),2) - xc(1,0)*xc(2,0) - (xc(0,1)
    - xc(1,1))*(xc(1,1) - xc(2,1))) - xc(0,1)*xc(1,2)*xc(2,1) +
    xc(1,1)*xc(1,2)*xc(2,1) - xc(0,0)*(xc(0,2)*(xc(1,0) - xc(2,0)) + xc(1,2)*xc(2,0)
    + xc(1,0)*(xc(1,2) - 2*xc(2,2))) + pow(xc(0,0),2)*(xc(1,2) - xc(2,2)) -
    pow(xc(0,1),2)*xc(2,2) - pow(xc(1,0),2)*xc(2,2) + 2*xc(0,1)*xc(1,1)*xc(2,2)
    - pow(xc(1,1),2)*xc(2,2)))/normcube;

  elematrix(11,8)=
    -((2*normal[2]*(xc(0,0) - xc(1,0)) - 2*normal[1]*(xc(0,1) -
    xc(1,1)))*(-(xc(1,1)*xc(2,0)) + xc(0,1)*(-xc(1,0) + xc(2,0)) + xc(0,0)*(xc(1,1)
    - xc(2,1)) + xc(1,0)*xc(2,1)))/(2.*normcube);

  elematrix(11,9)=0;

  elematrix(11,10)=0;

  elematrix(11,11)=0;
  return;
}


//=======================================================================


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::ConstraintElement3Register::Initialize(DRT::Discretization& dis)
{
  return 0;
}

#endif  // #ifdef CCADISCRET
