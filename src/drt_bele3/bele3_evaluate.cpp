/*!----------------------------------------------------------------------
\file bele3_evaluate.cpp

\maintainer Martin Pfaller

\brief dummy 3D boundary element without any physics
*----------------------------------------------------------------------*/

#include "bele3.H"
#include "../drt_so3/so_surface.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_timecurve.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_mat/newtonianfluid.H"



/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            gammi 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Bele3::Evaluate(Teuchos::ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    std::vector<int>&         lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  // get the required action
  std::string action = params.get<std::string>("action","none");

  if (action=="calc_struct_constrvol")
  {
    //create communicator
    const Epetra_Comm& Comm = discretization.Comm();

    //We are not interested in volume of ghosted elements
    if(Comm.MyPID()==Owner())
    {
      // element geometry update
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==Teuchos::null) dserror("Cannot get state vector 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      const int numdim = 3;
      LINALG::SerialDenseMatrix xscurr(NumNode(),numdim);  // material coord. of element
      SpatialConfiguration(xscurr,mydisp);
      //call submethod for volume evaluation and store rseult in third systemvector
      double volumeele = ComputeConstrVols(xscurr,NumNode());
      elevec3[0]= volumeele;
    }
  }
  return 0;
}


/*----------------------------------------------------------------------*
 |  do nothing (public)                                      a.ger 07/07|
 |                                                                      |
 |  The function is just a dummy.                                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Bele3::EvaluateNeumann(Teuchos::ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           std::vector<int>&         lm,
                                           Epetra_SerialDenseVector& elevec1,
                                           Epetra_SerialDenseMatrix* elemat1)
{
  return 0;
}


/*----------------------------------------------------------------------*
 * Compute Volume enclosed by surface.                          tk 10/07*
 * ---------------------------------------------------------------------*/
double DRT::ELEMENTS::Bele3::ComputeConstrVols
(
    const LINALG::SerialDenseMatrix& xc,
    const int numnode
)
{
  double V = 0.0;

  //Volume is calculated by evaluating the integral
  // 1/3*int_A(x dydz + y dxdz + z dxdy)

  // we compute the three volumes separately
  for (int indc = 0; indc < 3; indc++)
  {
    //split current configuration between "ab" and "c"
    // where a!=b!=c and a,b,c are in {x,y,z}
    LINALG::SerialDenseMatrix ab= xc;
    LINALG::SerialDenseVector c (numnode);
    for (int i = 0; i < numnode; i++)
    {
      ab(i,indc) = 0.0; // project by z_i = 0.0
      c(i) = xc(i,indc); // extract z coordinate
    }
    // index of variables a and b
    int inda = (indc+1)%3;
    int indb = (indc+2)%3;

    // get gaussrule
    const DRT::UTILS::IntegrationPoints2D  intpoints(getOptimalGaussrule(Shape()));
    int ngp = intpoints.nquad;

    // allocate vector for shape functions and matrix for derivatives
    LINALG::SerialDenseVector  funct(numnode);
    LINALG::SerialDenseMatrix  deriv(2,numnode);

    /*----------------------------------------------------------------------*
     |               start loop over integration points                     |
     *----------------------------------------------------------------------*/
    for (int gpid = 0; gpid < ngp; ++gpid)
    {
      const double e0 = intpoints.qxg[gpid][0];
      const double e1 = intpoints.qxg[gpid][1];

      // get shape functions and derivatives of shape functions in the plane of the element
      DRT::UTILS::shape_function_2D(funct,e0,e1,Shape());
      DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,Shape());

      double detA;
      // compute "metric tensor" deriv*ab, which is a 2x3 matrix with zero indc'th column
      LINALG::SerialDenseMatrix metrictensor(2,3);
      metrictensor.Multiply('N','N',1.0,deriv,ab,0.0);
      //LINALG::SerialDenseMatrix metrictensor(2,2);
      //metrictensor.Multiply('N','T',1.0,dxyzdrs,dxyzdrs,0.0);
      detA =  metrictensor(0,inda)*metrictensor(1,indb)-metrictensor(0,indb)*metrictensor(1,inda);
      const double dotprodc = funct.Dot(c);
      // add weighted volume at gausspoint
      V -= dotprodc*detA*intpoints.qwgt[gpid];

    }
  }
  return V/3.0;
}
