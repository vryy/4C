/*!-----------------------------------------------------------------------------------------------------------
 \file beam3_evaluate.cpp
 \brief

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

 *-----------------------------------------------------------------------------------------------------------*/
#ifdef D_BEAM3
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "beam3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../drt_lib/linalg_fixedsizematrix.H"
#include "../drt_inpar/inpar_statmech.H"



/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                                 cyron 01/08|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3::Evaluate(ParameterList& params,
    DRT::Discretization& discretization,
    vector<int>& lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::Beam3::ActionType act = Beam3::calc_none;
  // get the action required
  string action = params.get<string>("action","calc_none");
  if (action == "calc_none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff") act = Beam3::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff") act = Beam3::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = Beam3::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass") act = Beam3::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass") act = Beam3::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass") act = Beam3::calc_struct_nlnstifflmass; //with lumped mass matrix
  else if (action=="calc_struct_stress") act = Beam3::calc_struct_stress;
  else if (action=="calc_struct_eleload") act = Beam3::calc_struct_eleload;
  else if (action=="calc_struct_fsiload") act = Beam3::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep") act = Beam3::calc_struct_update_istep;
  else if (action=="calc_struct_update_imrlike") act = Beam3::calc_struct_update_imrlike;
  else if (action=="calc_struct_reset_istep") act = Beam3::calc_struct_reset_istep;
  else if (action=="calc_struct_ptcstiff")        act = Beam3::calc_struct_ptcstiff;
  else dserror("Unknown type of action for Beam3");

  switch(act)
  {
    case Beam3::calc_struct_ptcstiff:
    {
      const int nnode = NumNode();
      
    	switch(nnode)
   	  {
   	  		case 2:EvaluatePTC<2>(params, elemat1); break;
   	  		case 3:EvaluatePTC<3>(params, elemat1); break;
   	  		case 4:EvaluatePTC<4>(params, elemat1); break;  		
   	  		case 5:EvaluatePTC<5>(params, elemat1); break;  		
   	  		default:dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
   	  }
    }
    break;
    /*in case that only linear stiffness matrix is required b3_nlstiffmass is called with zero dispalcement and
     residual values*/
    case Beam3::calc_struct_linstiff:
    {
      //only nonlinear case implemented!
      dserror("linear stiffness matrix called, but not implemented");

    }
    break;

    //nonlinear stiffness and mass matrix are calculated even if only nonlinear stiffness matrix is required
    case Beam3::calc_struct_nlnstiffmass:
    case Beam3::calc_struct_nlnstifflmass:
    case Beam3::calc_struct_nlnstiff:
    case Beam3::calc_struct_internalforce:
    {
      // need current global displacement and residual forces and get them from discretization
      // making use of the local-to-global map lm one can extract current displacemnet and residual values for each degree of freedom
      //
      // get element displcements
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get state vectors 'displacement'");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

      // get residual displacements
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (res==null) dserror("Cannot get state vectors 'residual displacement'");
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);

      // get element velocities
      RefCountPtr<const Epetra_Vector> vel  = discretization.GetState("velocity");
      if (vel==null) dserror("Cannot get state vectors 'velocity'");
      vector<double> myvel(lm.size());
      DRT::UTILS::ExtractMyValues(*vel,myvel,lm);

      const int nnode = NumNode();

      if (act == Beam3::calc_struct_nlnstiffmass)
      {
    	  switch(nnode)
    	  {
  	  		case 2:
  	  		{
  	  			b3_nlnstiffmass<2>(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
  	  			break;
  	  		}
  	  		case 3:
  	  		{
  	  			b3_nlnstiffmass<3>(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
  	  			break;
  	  		}
  	  		case 4:
  	  		{
  	  			b3_nlnstiffmass<4>(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
  	  			break;
  	  		}
  	  		case 5:
  	  		{
  	  			b3_nlnstiffmass<5>(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
  	  			break;
  	  		}
  	  		default:
  	  			dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
    	  }
      }
      else if (act == Beam3::calc_struct_nlnstifflmass)
      {
    	  switch(nnode)
    	  {
  	  		case 2:
  	  		{
  	  			b3_nlnstiffmass<2>(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
  	  			lumpmass<2>(&elemat2);
  	  			break;
  	  		}
  	  		case 3:
  	  		{
  	  			b3_nlnstiffmass<3>(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
  	  			lumpmass<3>(&elemat2);
  	  			break;
  	  		}
  	  		case 4:
  	  		{
  	  			b3_nlnstiffmass<4>(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
  	  			lumpmass<4>(&elemat2);
  	  			break;
  	  		}
  	  		case 5:
  	  		{
  	  			b3_nlnstiffmass<5>(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
  	  			lumpmass<5>(&elemat2);
  	  			break;
  	  		}
  	  		default:
  	  			dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
    	  }
      }
      else if (act == Beam3::calc_struct_nlnstiff)
      {
    	  switch(nnode)
    	  {
  	  		case 2:
  	  		{
  	  			b3_nlnstiffmass<2>(params,myvel,mydisp,&elemat1,NULL,&elevec1);
  	  			break;
  	  		}
  	  		case 3:
  	  		{
  	  			b3_nlnstiffmass<3>(params,myvel,mydisp,&elemat1,NULL,&elevec1);
  	  			break;
  	  		}
  	  		case 4:
  	  		{
  	  			b3_nlnstiffmass<4>(params,myvel,mydisp,&elemat1,NULL,&elevec1);
  	  			break;
  	  		}
  	  		case 5:
  	  		{
  	  			b3_nlnstiffmass<5>(params,myvel,mydisp,&elemat1,NULL,&elevec1);
  	  			break;
  	  		}
  	  		default:
  	  			dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
    	  }
      }

      else if (act == Beam3::calc_struct_internalforce)
      {
    	  switch(nnode)
    	  {
  	  		case 2:
  	  		{
  	  			b3_nlnstiffmass<2>(params,myvel,mydisp,NULL,NULL,&elevec1);
  	  			break;
  	  		}
  	  		case 3:
  	  		{
  	  			b3_nlnstiffmass<3>(params,myvel,mydisp,NULL,NULL,&elevec1);
  	  			break;
  	  		}
  	  		case 4:
  	  		{
  	  			b3_nlnstiffmass<4>(params,myvel,mydisp,NULL,NULL,&elevec1);
  	  			break;
  	  		}
  	  		case 5:
  	  		{
  	  			b3_nlnstiffmass<5>(params,myvel,mydisp,NULL,NULL,&elevec1);
  	  			break;
  	  		}
  	  		default:
  	  			dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
    	  }
      }

  /*at the end of an iteration step the geometric configuration has to be updated: the starting point for the
   * next iteration step is the configuration at the end of the current step */
  Qold_ = Qnew_;
  curvold_ = curvnew_;
  thetaold_ = thetanew_;
  thetaprimeold_ = thetaprimenew_;

/*
    //the following code block can be used to check quickly whether the nonlinear stiffness matrix is calculated
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

      		//create new displacement and velocity vectors in order to store artificially modified displacements
      		vector<double> vel_aux(6*nnode);
      		vector<double> disp_aux(6*nnode);

      			DRT::UTILS::ExtractMyValues(*disp,disp_aux,lm);
      			DRT::UTILS::ExtractMyValues(*vel,vel_aux,lm);

      		//modifying displacement artificially (for numerical derivative of internal forces):
      		disp_aux[6*k + i] += h_rel;
      		vel_aux[6*k + i] += h_rel / params.get<double>("delta time",0.01);
      		
		  //b3_nlnstiffmass is a templated function. therefore we need to point out the number of nodes in advance
        	  switch(nnode)
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
        	  		default:
        	  			dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
        	  }

        	//computing derivative d(fint)/du numerically by finite difference
      		for(int u = 0 ; u < 6*nnode ; u++ )
      			stiff_approx(u,k*6+i)= ( pow(force_aux[u],2) - pow(elevec1(u),2) )/ (h_rel * (force_aux[u] + elevec1(u) ) );

      	} //for(int k=0; k<nnode; k++)//for all nodes

      } //for(int i=0; i<3; i++) //for all dof
      

      for(int line=0; line<6*nnode; line++)
      {
      	for(int col=0; col<6*nnode; col++)
      	{
      		stiff_relerr(line,col)= fabs( ( pow(elemat1(line,col),2) - pow(stiff_approx(line,col),2) )/ ( (elemat1(line,col) + stiff_approx(line,col)) * elemat1(line,col) ));

      		//suppressing small entries whose effect is only confusing and NaN entires (which arise due to zero entries)
      		if ( fabs( stiff_relerr(line,col) ) < h_rel*500 || isnan( stiff_relerr(line,col)) || elemat1(line,col) == 0) //isnan = is not a number
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
    case calc_struct_update_istep:
    case calc_struct_update_imrlike:
    {
      /*the action calc_struct_update_istep is called in the very end of a time step when the new dynamic
       * equilibrium has finally been found; this is the point where the variable representing the geometric
       * status of the beam have to be updated; the geometric status is represented by means of the triads Tnew_,
       * the curvatures curvnew_ and the angular values thetaanew_ and thetaprimenew_*/
      Qconv_ = Qnew_;
      curvconv_ = curvnew_;
      thetaconv_ = thetanew_;
      thetaprimeconv_ = thetaprimenew_;
    }
    break;
    case calc_struct_reset_istep:
    {
      /*the action calc_struct_reset_istep is called by the adaptive time step controller; carries out one test
       * step whose purpose is only figuring out a suitabel timestep; thus this step may be a very bad one in order
       * to iterated towards the new dynamic equilibrium and the thereby gained new geometric configuration should
       * not be applied as starting point for any further iteration step; as a consequence the thereby generated change
       * of the geometric configuration should be canceled and the configuration should be reset to the value at the
       * beginning of the time step*/
      Qold_ = Qconv_;
      curvold_ = curvconv_;
      thetaold_ = thetaconv_;
      thetaprimeold_ = thetaprimeconv_;
    }
    break;
    case calc_struct_stress:
      dserror("No stress output implemented for beam3 elements");
    default:
      dserror("Unknown type of action for Beam3 %d", act);
  }
  return 0;

}

/*-----------------------------------------------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)                                       cyron 03/08|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Beam3::EvaluateNeumann(ParameterList& params,
                                        DRT::Discretization& discretization,
                                        DRT::Condition& condition,
                                        vector<int>& lm,
                                        Epetra_SerialDenseVector& elevec1,
                                        Epetra_SerialDenseMatrix* elemat1)
{  
  // get element displacements
  RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
  if (disp==null) dserror("Cannot get state vector 'displacement'");
  vector<double> mydisp(lm.size());
  DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
  // get element velocities (UNCOMMENT IF NEEDED)
  /*
  RefCountPtr<const Epetra_Vector> vel  = discretization.GetState("velocity");
  if (vel==null) dserror("Cannot get state vectors 'velocity'");
  vector<double> myvel(lm.size());
  DRT::UTILS::ExtractMyValues(*vel,myvel,lm);
  */

  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const vector<int>* curve = condition.Get<vector<int> >("curve");
  int curvenum = -1;
  // number of the load curve related with a specific line Neumann condition called
  if (curve) curvenum = (*curve)[0];
  // amplitude of load curve at current time called
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)//notation for this function similar to Crisfield, Volume 1;
  curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

  // no. of nodes on this element; the following line is only valid for elements with constant number of
  // degrees of freedom per node
  const int numdf = 6;
  const DiscretizationType distype = this->Shape();

  // gaussian points
  const DRT::UTILS::IntegrationPoints1D intpoints(MyGaussRule(NumNode(),gaussunderintegration));

  //declaration of variable in order to store shape function
  Epetra_SerialDenseVector funct(NumNode());

  // get values and switches from the condition

  // onoff is related to the first 6 flags of a line Neumann condition in the input file;
  // value 1 for flag i says that condition is active for i-th degree of freedom
  const vector<int>* onoff = condition.Get<vector<int> >("onoff");
  // val is related to the 6 "val" fields after the onoff flags of the Neumann condition
  // in the input file; val gives the values of the force as a multiple of the prescribed load curve
  const vector<double>* val = condition.Get<vector<double> >("val");

  //integration loops
  for (int numgp=0; numgp<intpoints.nquad; ++numgp)
  {
    //integration points in parameter space and weights
    const double xi = intpoints.qxg[numgp][0];
    const double wgt = intpoints.qwgt[numgp];

    //evaluation of shape funcitons at Gauss points
    DRT::UTILS::shape_function_1D(funct,xi,distype);

    double fac=0;
    fac = wgt * alpha_[numgp];

    // load vector ar
    double ar[numdf];

    // loop the dofs of a node
    for (int dof=0; dof<numdf; ++dof)
      ar[dof] = fac * (*onoff)[dof]*(*val)[dof]*curvefac;


    //sum up load components
    for (int node=0; node<NumNode(); ++node)
      for (int dof=0; dof<numdf; ++dof)
        elevec1[node*numdf+dof] += funct[node] *ar[dof];

  } // for (int numgp=0; numgp<intpoints.nquad; ++numgp)

  return 0;
}

/*-----------------------------------------------------------------------------------------------------------*
 | auxiliary functions for dealing with large rotations and nonlinear stiffness                    cyron 04/08|
 *----------------------------------------------------------------------------------------------------------*/
//computing basis of stiffness matrix of Crisfield, Vol. 2, equation (17.105)
template<int nnode>
inline void DRT::ELEMENTS::Beam3::computestiffbasis(const LINALG::Matrix<3,3>& Tnew, const LINALG::Matrix<3,1>& Cm, const LINALG::Matrix<3,1>& Cb, const LINALG::Matrix<3,3>& S, LINALG::Matrix<6*nnode,6*nnode>& stiffmatrix, const LINALG::Matrix<1,nnode>& funct, const LINALG::Matrix<1,nnode>& deriv)
{
  //calculating the first matrix of (17.105) directly involves multiplications of large matrices (e.g. with the 6*nnode x 6 - matrix X)
  //application of the definitions in (17.100) allows blockwise evaluation with multiplication and addition of 3x3-matrices only
  //in the follwoing all the blocks are treated separately; their name is related directly to their content:
  //e.g. TCmTt is the product of the 3 matrices T * C_m * T^t (with T and C_m according to (17.74) and (17.76)
  //for the blockwise calculation on which the following steps are based on the relation S^t = -S for spin matrices was applied

  LINALG::Matrix<3,3> TCmTt;
  LINALG::Matrix<3,3> TCbTt;
  LINALG::Matrix<3,3> StTCmTt;
  LINALG::Matrix<3,3> StTCmTtS;
  LINALG::Matrix<3,3> TCmTtS;

  //calculating TCmTt and TCbTt
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      TCmTt(i,j) = 0.0;
      TCbTt(i,j) = 0.0;
      for (int k = 0; k < 3; ++k)
      {
        TCmTt(i,j) += Tnew(i,k)*Cm(k)*Tnew(j,k);
        TCbTt(i,j) += Tnew(i,k)*Cb(k)*Tnew(j,k);
      }
    }
  }

  //calculating StTCmTt
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      StTCmTt(i,j) = 0.0;
      for (int k = 0; k < 3; ++k)
        StTCmTt(i,j) += S(k,i)*TCmTt(k,j);
    }
  }

  //calculating StTCmTtS
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      StTCmTtS(i,j) = 0.0;
      for (int k = 0; k < 3; ++k)
        StTCmTtS(i,j) += StTCmTt(i,k)*S(k,j);
    }
  }

  //calculating TCmTtS
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      TCmTtS(i,j) = 0.0;
      for (int k = 0; k < 3; ++k)
        TCmTtS(i,j) += TCmTt(i,k)*S(k,j);
    }
  }

  //calculating basis of stiffness matrix by means of above blocks
  //the calculations are based on the work by jbueckle
  for (int n = 0; n < nnode; ++n)
  {
    for (int m = 0; m < nnode; ++m)
    {
    	for (int i = 0; i < 3; ++i)
    	{
    		for (int j = 0; j < 3; ++j)
    		{
    		      stiffmatrix(6 * n + i  ,6 * m + j)   =  deriv(n) * deriv(m) * TCmTt(i,j);
    		      stiffmatrix(6 * n + 3 + i  ,6 * m + 3 + j)   =  funct(n) * funct(m) * StTCmTtS(i,j) + deriv(n) * deriv(m) * TCbTt(i,j);
    		      stiffmatrix(6 * n + 3 + i  ,6 * m + j)   =  funct(n) * deriv(m) * StTCmTt(i,j);
    		      stiffmatrix(6 * n + i  ,6 * m + 3 + j)   =  deriv(n) * funct(m) * TCmTtS(i,j);
    		}
    	}
    }
  }
 
  return;
} // DRT::ELEMENTS::Beam3::computestiffbasis




/*---------------------------------------------------------------------------*
 |computes from a quaternion q rodrigues parameters omega (public)cyron02/09|
 *---------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Beam3::quaterniontorodrigues(const LINALG::Matrix<4,1>& q, LINALG::Matrix<3,1>& omega)
{
  /*the Rodrigues parameters are defined only for angles whose absolute valued is smaller than PI, i.e. for which
   * the fourth component of the quaternion is unequal zero; if this is not satisfied for the quaternion passed into
   * this method an error is thrown*/
  if(q(3) == 0)
    dserror("cannot compute Rodrigues parameters for angles with absolute valued PI !!!");

  //in any case except for the one dealt with above the angle can be computed from a quaternion via Crisfield, Vol. 2, eq. (16.79)
  for(int i = 0; i<3; i++)
    omega(i) = q(i)*2/q(3);


  return;
} //DRT::ELEMENTS::Beam3::quaterniontorodrigues




/*----------------------------------------------------------------------*
 |computes from a quaternion q the related angle theta (public)cyron10/08|
 *----------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Beam3::quaterniontoangle(const LINALG::Matrix<4,1>& q, LINALG::Matrix<3,1>& theta)
{
  /*the following function computes from a quaternion q an angle theta within [-PI; PI]; such an interval is
   * imperative for the use of the resulting angle together with formulae like Crisfield, Vol. 2, equation (16.90);
   * note that these formulae comprise not only trigonometric functions, but rather the angle theta directly. Hence
   * they are not 2*PI-invariant !!! */

  //first we consider the case that the absolute value of the rotation angle equals zero
  if(q(0) == 0 && q(1) == 0 && q(2) == 0 )
  {
    for(int i = 0; i<3; i++)
      theta(i) = 0;

    return;
  }

  //second we consider the case that the abolute value of the rotation angle equals PI
  if(q(3) == 0)
  {
    //note that with q(3) == 0 the first three elements of q represent the unit direction vector of the angle
    //according to Crisfield, Vol. 2, equation (16.67)
    for(int i = 0; i<3; i++)
      theta(i) = q(i) * PI;

    return;
  }

  //in any case except for the one dealt with above the angle can be computed from a quaternion via Crisfield, Vol. 2, eq. (16.79)
  LINALG::Matrix<3,1> omega;
  for(int i = 0; i<3; i++)
    omega(i) = q(i)*2/q(3);

  double tanhalf = omega.Norm2() / 2;

  double thetaabs = atan(tanhalf)*2;

  for(int i = 0; i<3; i++)
      theta(i) = thetaabs* omega(i) / omega.Norm2();

  return;
} //DRT::ELEMENTS::Beam3::quaterniontoangle()




/*---------------------------------------------------------------------------*
 |computes a spin matrix out of a rotation vector 		   (public)cyron02/09|
 *---------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Beam3::computespin(LINALG::Matrix<3,3>& spin, LINALG::Matrix<3,1> rotationangle, const double& spinscale)
{
  //function based on Crisfield Vol. 2, Section 16 (16.8)
  rotationangle.Scale(spinscale);
  spin(0,0) = 0;
  spin(0,1) = -rotationangle(2);
  spin(0,2) = rotationangle(1);
  spin(1,0) = rotationangle(2);
  spin(1,1) = 0;
  spin(1,2) = -rotationangle(0);
  spin(2,0) = -rotationangle(1);
  spin(2,1) = rotationangle(0);
  spin(2,2) = 0;

  return;
} // DRT::ELEMENTS::Beam3::computespin




/*----------------------------------------------------------------------*
 |computes a rotation matrix R from a quaternion q						|
 |cf. Crisfield, Vol. 2, equation (16.70) 			  (public)cyron10/08|
 *----------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Beam3::quaterniontotriad(const LINALG::Matrix<4,1>& q, LINALG::Matrix<3,3>& R)
{
  //separate storage of vector part of q
  LINALG::Matrix<3,1> qvec;
  for(int i = 0; i<3; i++)
    qvec(i) = q(i);

  //setting R to third summand of equation (16.70)
  computespin(R, qvec, 2*q(3));

  //adding second summand of equation (16.70)
  for(int i = 0; i<3; i++)
    for(int j = 0; j<3; j++)
      R(i,j) += 2*q(i)*q(j);

  //adding diagonal entries according to first summand of equation (16.70)
  R(0,0) = 1 - 2*(q(1)*q(1) + q(2)*q(2));
  R(1,1) = 1 - 2*(q(0)*q(0) + q(2)*q(2));
  R(2,2) = 1 - 2*(q(0)*q(0) + q(1)*q(1));

  return;
} // DRT::ELEMENTS::Beam3::quaterniontotriad




/*---------------------------------------------------------------------------*
 |computes a quaternion from an angle vector 		  	   (public)cyron02/09|
 *---------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3::angletoquaternion(const LINALG::Matrix<3,1>& theta, LINALG::Matrix<4,1>& q)
{
  //absolute value of rotation angle theta
  double abs_theta = theta.Norm2();

  //computing quaterion for rotation by angle theta, Crisfield, Vol. 2, equation (16.67)
  if (abs_theta > 0)
  {
    q(0) = theta(0) * sin(abs_theta / 2) / abs_theta;
    q(1) = theta(1) * sin(abs_theta / 2) / abs_theta;
    q(2) = theta(2) * sin(abs_theta / 2) / abs_theta;
    q(3) = cos(abs_theta / 2);
  }
  else
  {
    q.PutScalar(0.0);
    q(3) = 1;
  }

  return;
}// DRT::ELEMENTS::Beam3::angletoquaternion




/*---------------------------------------------------------------------------*
 |computes a quaternion q from a rotation matrix R; all operations are      |
 |performed according to Crisfield, Vol. 2, section 16.10 and the there      |
 |described Spurrier's algorithm		  	   			   (public)cyron02/09|
 *---------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3::triadtoquaternion(const LINALG::Matrix<3,3>& R, LINALG::Matrix<4,1>& q)
{
  double trace = R(0,0) + R(1,1) + R(2,2);
  if(trace>R(0,0)  && trace>R(1,1) && trace>R(2,2))
  {
    q(3) = 0.5 * pow(1 + trace, 0.5);
    q(0) = (R(2,1) - R(1,2)) / (4*q(3));
    q(1) = (R(0,2) - R(2,0)) / (4*q(3));
    q(2) = (R(1,0) - R(0,1)) / (4*q(3));
  }
  else
  {
    for(int i = 0 ; i<3 ; i++)
    {
      int j = (i+1)% 3;
      int k = (i+2)% 3;

      if(R(i,i) >= R(j,j) && R(i,i) >= R(k,k))
      {
        //equation (16.78a)
        q(i) = pow(0.5*R(i,i) + 0.25*(1 - trace) , 0.5);

        //equation (16.78b)
        q(3) = 0.25*(R(k,j) - R(j,k)) / q(i);

        //equation (16.78c)
        q(j) = 0.25*(R(j,i) + R(i,j)) / q(i);
        q(k) = 0.25*(R(k,i) + R(i,k)) / q(i);
       }
     }
   }
  return;
}// DRT::ELEMENTS::Beam3::TriadToQuaternion




/*---------------------------------------------------------------------------*
 |matrix H^(-1) which turns non-additive spin variables into additive ones   |
 |according to Crisfield, Vol. 2, equation (16.93)	  	   (public)cyron02/09|
 *---------------------------------------------------------------------------*/
LINALG::Matrix<3,3> DRT::ELEMENTS::Beam3::Hinv(LINALG::Matrix<3,1> theta)
{
  LINALG::Matrix<3,3> result;
  double theta_abs = pow(theta(0)*theta(0) + theta(1)*theta(1) + theta(2)*theta(2) ,0.5);

  //in case of theta_abs == 0 the following computation has problems with singularities
  if(theta_abs > 0)
  {
    computespin(result, theta, -0.5);

    for(int i = 0; i<3; i++)
      result(i,i) += theta_abs/( 2*tan(theta_abs/2) );

    for(int i = 0; i<3; i++)
      for(int j=0; j<3; j++)
        result(i,j) += theta(i) * theta(j) * (1 - theta_abs/(2*tan(theta_abs/2)) )/pow(theta_abs,2);
  }
  //in case of theta_abs == 0 H(theta) is the identity matrix and hence also Hinv
  else
  {
    result.PutScalar(0.0);
    for(int j=0; j<3; j++)
      result(j,j) = 1;
  }

  return result;
}// DRT::ELEMENTS::Beam3::Hinv




/*---------------------------------------------------------------------------*
 |this function performs an update of the central triad as in principle 	 |
 |given in Crisfield, Vol. 2, equation (17.65), but by means of a			 |
 |quaternion product and then calculation of the equivalent rotation matrix   |
 | according to eq. (16.70)								   (public)cyron02/09|
 *---------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Beam3::updatetriad(const LINALG::Matrix<3,1>& deltatheta, LINALG::Matrix<3,3>& Tnew, const int numgp)
{
  //computing quaternion equivalent to rotation by deltatheta
  LINALG::Matrix<4,1> Qrot;
  angletoquaternion(deltatheta,Qrot);

  //computing quaternion Qnew_ for new configuration of Qold_ for old configuration by means of a quaternion product
  quaternionproduct(Qold_[numgp],Qrot,Qnew_[numgp]);

  //normalizing quaternion in order to make sure that it keeps unit absolute values through time stepping
  double abs = Qnew_[numgp].Norm2();
  for(int i = 0; i<4; i++)
    Qnew_[numgp](i) = Qnew_[numgp](i) / abs;

  quaterniontotriad(Qnew_[numgp],Tnew);

  return;
} //DRT::ELEMENTS::Beam3::updatetriad

/*-----------------------------------------------------------------------------------*
 |computes inverse quaternion q^{-1} for input quaternion q 	   (public)cyron02/09|
 *-----------------------------------------------------------------------------------*/
inline LINALG::Matrix<4,1> DRT::ELEMENTS::Beam3::inversequaternion(const LINALG::Matrix<4,1>& q)
{
  //square norm ||q||^2 of quaternion q
  double qnormsq = q.Norm2() * q.Norm2();

  //declaration of variable for inverse quaternion
  LINALG::Matrix<4,1> qinv;

  //inverse quaternion q^(-1) = [-q0, -q1, -q2, q3] / ||q||^2;
  for(int i = 0; i<3; i++)
    qinv(i) = -q(i) / qnormsq;

  qinv(3) = q(3) / qnormsq;

  return qinv;

} //DRT::ELEMENTS::Beam3::inversequaternion

/*---------------------------------------------------------------------------------------------------*
 |quaternion product q12 = q2*q1, Crisfield, Vol. 2, equation (16.71)		  	   (public)cyron02/09|
 *---------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Beam3::quaternionproduct(const LINALG::Matrix<4,1>& q1,const LINALG::Matrix<4,1>& q2,LINALG::Matrix<4,1>& q12)
{
  q12(0) = q2(3)*q1(0) + q1(3)*q2(0) + q2(1)*q1(2) - q1(1)*q2(2);
  q12(1) = q2(3)*q1(1) + q1(3)*q2(1) + q2(2)*q1(0) - q1(2)*q2(0);
  q12(2) = q2(3)*q1(2) + q1(3)*q2(2) + q2(0)*q1(1) - q1(0)*q2(1);
  q12(3) = q2(3)*q1(3) - q2(2)*q1(2) - q2(1)*q1(1) - q2(0)*q1(0);
} //DRT::ELEMENTS::Beam3::quaternionproduct


/*-------------------------------------------------------------------------------------------------------*
 |updating local curvature according to Crisfield, Vol. 2, pages 209 - 210; an exact update of  	     |
 | the curvature is computed by means of equation (16.148) instead of an approximated one as given by 	 |
 | equs. (17.72) and (17.73); note that this function requires as input parameters the angle delta theta |
 | from (16.146), which is responsible for the rotation of the triad at the Gauss point as well as its   |
 | derivative with respect to the curve parameter s, i.e. d(delta theta)/ds. This derivative is denoted  |
 | as deltathetaprime																   (public)cyron02/09|
 *-------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Beam3::updatecurvature(const LINALG::Matrix<3,3>& Tnew, LINALG::Matrix<3,1>& deltatheta,LINALG::Matrix<3,1>& deltathetaprime, const int numgp)
{
  //declaration of local variables:
  LINALG::Matrix<3,1> omega; //omega according to (16.146)
  LINALG::Matrix<3,1> omegaprime; //omega' according to (16.147)
  LINALG::Matrix<3,1> chignl; //chi_gnl according to (16.148)

  //absolute value of rotation vector delta theta
  double abs_theta = deltatheta.Norm2();

  //as in (16.147) division by delta theta arises we have to assume that this angle is unequal to zero
  if (abs_theta > 0)
  {
    //compute omega according to (16.146)
    omega = deltatheta;
    omega.Scale(2*tan(0.5*abs_theta) / abs_theta);

    //in (16.147) the brackets contain a 3x3 matrix which is denoted as Aux here and computed now:
    LINALG::Matrix<3,3> Aux;
    for(int i = 0; i<3; i++)
    {
      for(int j = 0; j<3; j++)
      {
        Aux(i,j) = -(1 - abs_theta / sin(abs_theta) ) * deltatheta(i)*deltatheta(j) / pow(abs_theta,2);
        if(i==j)
          Aux(i,j) += 1;
      }
    }

    omegaprime.Multiply(Aux,deltathetaprime);
    /*we apply the prefactor of (16.147); here we need an angle -PI < theta < PI; note that theta may be assumed to
     * lie within this domain as otherwise the element would rotate at a specific Gauss point by more that 180° within
     * one single iteration step which would disrupt convergence anyway; note that one could ensure -PI < theta < PI
     * easily by applying the three code lines
     *
     *
     *  LINALG::Matrix<4,1> q;
     *   angletoquaternion(deltatheta,q);
     *  quaterniontoangle(q,deltatheta);
     *
     * in the very beginning of this method. However, for the above reason we assume that theta lies always in the proper
     * region and thus save the related compuational cost and just throw an error if this prerequesite is unexpectedly not
     * satisfied */
    if(abs_theta > PI)
      dserror("delta theta exceeds region for which equation (16.147) is valid");
    omegaprime.Scale(2*tan(abs_theta / 2) / abs_theta);

    //compute chignl from omega and omega' according to (16.148)
    chignl(0) = 0.5*(omega(1)*omegaprime(2) - omega(2)*omegaprime(1)) ;
    chignl(1) = 0.5*(omega(2)*omegaprime(0) - omega(0)*omegaprime(2)) ;
    chignl(2) = 0.5*(omega(0)*omegaprime(1) - omega(1)*omegaprime(0)) ;
    chignl += omegaprime;
    chignl.Scale( 1/(1 + pow(tan(abs_theta/2),2) ));

  }
  //with delta theta == 0 we get omega == 0 and thus according to (16.147) omega' = d(delta theta)/ds
  else
    chignl = deltathetaprime;

  curvnew_[numgp].MultiplyTN(Tnew,chignl);
  curvnew_[numgp] += curvold_[numgp];

  return;
} //DRT::ELEMENTS::Beam3::updatecurvature




/*------------------------------------------------------------------------------------------------------*
 |updating local curvature according approximately to Crisfield, Vol. 2, eqs. (17.72) and (17.73); note:|
 |this update is suitable for beams with linear shape functions, only 				  (public)cyron02/09|
 *------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Beam3::approxupdatecurvature(const LINALG::Matrix<3,3>& Tnew, LINALG::Matrix<3,1> deltatheta,LINALG::Matrix<3,1> deltathetaprime, const int numgp)
{
  //old triad
  LINALG::Matrix<3,3> Told;
  quaterniontotriad(Qold_[numgp],Told);

  //compute spin matrix from eq. (17.73)
  LINALG::Matrix<3,3> spin;
  computespin(spin, deltatheta, 1.0);

  //turning spin matrix to left right hand side matrix of eq. (17.73)
  for(int i = 0; i<3; i++)
    spin(i,i) += 1;

  //complete right hand side matrix of eq. (17.73)
  //mid point triad
  LINALG::Matrix<3,3> Tmid;
  Tmid.Multiply(spin,Told);

  //eq. (17.72)
  curvnew_[numgp].MultiplyTN(Tmid,deltathetaprime);

  curvnew_[numgp] += curvold_[numgp];

  return;
} //DRT::ELEMENTS::Beam3::approxupdatecurvature




/*-----------------------------------------------------------------------------------------------------*
 |computes stiffness matrix Ksigma1 according to Crisfield, Vol. 2, equation (17.106) (public)cyron02/09|
 *-----------------------------------------------------------------------------------------------------*/
template<int nnode>
inline void DRT::ELEMENTS::Beam3::computeKsig1(LINALG::Matrix<6*nnode,6*nnode>& Ksig1, const LINALG::Matrix<3,1>& stressn, const LINALG::Matrix<3,1>& stressm, const LINALG::Matrix<1,nnode>& funct, const LINALG::Matrix<1,nnode>& deriv)
{

  LINALG::Matrix<3,3> Sn;
  LINALG::Matrix<3,3> Sm;
  //first we caluclate S(n) and S(m)
  computespin(Sn,stressn,1.0);
  computespin(Sm,stressm,1.0);

  //then we insert them blockwise in Ksig1
  for (int n = 0; n < nnode; ++n)
  {
    for (int m = 0; m < nnode; ++m)
    {
    	for (int i = 0; i < 3; ++i)
    	{
    		for (int j = 0; j < 3; ++j)
    		{
		      Ksig1(6 * n + i  ,6 * m + j)   =  0;
		      Ksig1(6 * n + 3 + i  ,6 * m + 3 + j)   =  -deriv(n) * funct(m) * Sm(i,j);
		      Ksig1(6 * n + 3 + i  ,6 * m + j)   =  0;
		      Ksig1(6 * n + i  ,6 * m + 3 + j)   =  -deriv(n) * funct(m) * Sn(i,j);
    		}
    	}
    }
  }
  return;
} //DRT::ELEMENTS::Beam3::computeKsig1





/*----------------------------------------------------------------------------------------------------------------------*
 |computes stiffness matrix Ksigma2 according to Crisfield, Vol. 2, equation (17.107a) and (17.107b) (public)  cyron02/09|
 *----------------------------------------------------------------------------------------------------------------------*/
template<int nnode>
inline void DRT::ELEMENTS::Beam3::computeKsig2(LINALG::Matrix<6*nnode,6*nnode>& Ksig2, const LINALG::Matrix<3,1>& stressn, const LINALG::Matrix<3,3>& S, const LINALG::Matrix<1,nnode>& funct, const LINALG::Matrix<1,nnode>& deriv)
{

  LINALG::Matrix<3,3> Sn;
  LINALG::Matrix<3,3> Y;
  //first we compute S(n) and Y according to (17.107b)
  computespin(Sn,stressn, 1.0);
  Y.Multiply(S,Sn);

  //then we insert them into Ksig2 by means of a blockwise algorithm
  for (int n = 0; n < nnode; ++n)
  {
    for (int m = 0; m < nnode; ++m)
    {
    	for (int i = 0; i < 3; ++i)
    	{
    		for (int j = 0; j < 3; ++j)
    		{
		      Ksig2(6 * n + i  ,6 * m + j)   =  0;
		      Ksig2(6 * n + 3 + i  ,6 * m + 3 + j)   =  funct(n)*funct(m)*Y(i,j);
		      Ksig2(6 * n + 3 + i  ,6 * m + j)   =  funct(n)*deriv(m)*Sn(i,j);
		      Ksig2(6 * n + i  ,6 * m + 3 + j)   =  0;
    		}
    	}
    }
  }
  return;
} // DRT::ELEMENTS::Beam3::computeKsig2


/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   cyron 01/08|
 *-----------------------------------------------------------------------------------------------------------*/
template<int nnode>
void DRT::ELEMENTS::Beam3::b3_nlnstiffmass( ParameterList& params,
                                            vector<double>&           vel,
                                            vector<double>&           disp,
                                            Epetra_SerialDenseMatrix* stiffmatrix,
                                            Epetra_SerialDenseMatrix* massmatrix,
                                            Epetra_SerialDenseVector* force)
{
  //constitutive laws from Crisfield, Vol. 2, equation (17.76)
  LINALG::Matrix<3,1> Cm;
  LINALG::Matrix<3,1> Cb;

  //normal/shear strain and bending strain(curvature)
  LINALG::Matrix<3,1> epsilonn;
  LINALG::Matrix<3,1> epsilonm;

  //stress values n and m, Crisfield, Vol. 2, equation (17.78)
  LINALG::Matrix<3,1> stressn;
  LINALG::Matrix<3,1> stressm;

  //derivative of x with respect to xi
  LINALG::Matrix<3,1> dxdxi_gp;

  //theta interpolated at the GP. Note: theta is not the current angle theta but only the difference to the reference configuration
  LINALG::Matrix<3,1> thetanew_gp;

  //derivative of theta with respect to xi
  LINALG::Matrix<3,1> thetaprimenew_gp;

  //stiffness base of stiffness matrix, Crisfield Vol. 2, equation (17.105)
  LINALG::Matrix<6*nnode,6*nnode> Kstiff_gp;

  //nonlinear parts of stiffness matrix, Crisfield Vol. 2, equation (17.83) and (17.87)
  LINALG::Matrix<6*nnode,6*nnode> Ksig1_gp;
  LINALG::Matrix<6*nnode,6*nnode> Ksig2_gp;


  //triad at GP, Crisfiel Vol. 2, equation (17.73)
  LINALG::Matrix<3,3> Tnew;


  //first of all we get the material law
  Teuchos::RCP<const MAT::Material> currmat = Material();
  double ym = 0;
  double sm = 0;
  double density = 0;

  //assignment of material parameters; only St.Venant material is accepted for this beam
  switch(currmat->MaterialType())
  {
    case INPAR::MAT::m_stvenant:// only linear elastic material supported
    {
      const MAT::StVenantKirchhoff* actmat = static_cast<const MAT::StVenantKirchhoff*>(currmat.get());
      ym = actmat->Youngs();
      sm = ym / (2*(1 + actmat->PoissonRatio()));
      density = actmat->Density();
    }
    break;
    default:
    dserror("unknown or improper type of material law");
  }

  //"new" variables have to be adopted to displacement passed in from BACI driver

  //Get integrationpoints for underintegration
  DRT::UTILS::IntegrationPoints1D gausspoints(MyGaussRule(nnode,gaussunderintegration));

  //Get DiscretizationType
  const DRT::Element::DiscretizationType distype = Shape();

  //Matrices for h and h,xi
  LINALG::Matrix<1,nnode> funct;
  LINALG::Matrix<1,nnode> deriv;

  //Loop through all GP and calculate their contribution to the forcevector and stiffnessmatrix
  for(int numgp=0; numgp < gausspoints.nquad; numgp++)
  {
  	//Get location and weight of GP in parameter space
  	const double xi = gausspoints.qxg[numgp][0];
  	const double wgt = gausspoints.qwgt[numgp];

  	//Get h and h,xi
  	DRT::UTILS::shape_function_1D(funct,xi,distype);
  	DRT::UTILS::shape_function_1D_deriv1(deriv,xi,distype);

  	dxdxi_gp.Clear();
  	thetanew_gp.Clear();
  	thetaprimenew_gp.Clear();

	  //set up current dxdxi, theta, and thetaprime at the GP
		for (int dof=0; dof<3; ++dof)//j
		{
			for (int node=0; node<nnode; ++node)
			{
				 dxdxi_gp(dof) += (Nodes()[node]->X()[dof]+disp[6*node+dof])*deriv(node); 

				/*compute interpolated angle displacemnt at specific Gauss point; angle displacement is
				 * taken from discretization and interpolated with the according shapefunctions*/
				thetanew_gp(dof)          += disp[6*node+3+dof]*funct(node);

				/*compute derivative of interpolated angle displacement with respect to curve parameter at specific Gauss point; angle displacement is
				 * taken from discretization and interpolated with the according shapefunctions*/
				thetaprimenew_gp(dof)     += disp[6*node+3+dof]*deriv(node)/alpha_[numgp];

			}//for (int node=0; node<nnode; ++node)
		}//for (int dof=0; dof<3; ++dof)

		//update theta and thetaprime
		thetanew_[numgp]= thetanew_gp;
		thetaprimenew_[numgp]= thetaprimenew_gp;

		/*to perform a curvature update at a specific Gauss point we need to know the angle deltatheta by which the triad at that
		 * Gauss point is rotated and furthermore the derivative of this rotation angle along the curve. The latter one is
		 * denoted as d(delta theta)/ds in Crisfield, Vol. 2, page 209, and can be computed according to the comments in the bottom
		 * of this page*/

		//compute delta theta for (16.146) at the Gauss point only according to remark in the bottom of page 209
		LINALG::Matrix<3,1> deltatheta_gp = thetanew_[numgp];
		deltatheta_gp -= thetaold_[numgp];

		//compute d(delta theta)/ds for (16.147) at the Gauss point only according to remark in the bottom of page 209
		LINALG::Matrix<3,1> deltathetaprime_gp = thetaprimenew_[numgp];
		deltathetaprime_gp -= thetaprimeold_[numgp];

		/*update triad at Gauss point as shown in Crisfield, Vol. 2, equation (17.65) Note: instead of the matrix multiplication of (17.65) we use
		 *a quaternion product*/
		updatetriad(deltatheta_gp,Tnew,numgp);

		//updating local curvature according to Crisfield, Vol. 2, pages 209 - 210
		updatecurvature(Tnew,deltatheta_gp,deltathetaprime_gp,numgp);

		epsilonn.Clear();

		//computing current axial and shear strain epsilon, Crisfield, Vol. 2, equation (17.97)
		epsilonn.MultiplyTN(Tnew,dxdxi_gp);
		epsilonn.Scale(1/alpha_[numgp]);
		epsilonn(0) -=  1.0;

		epsilonm.Clear();

		//computing spin matrix S(dxdxi_gp) according to Crisfield, Vol. 2, equation (17.100)
		LINALG::Matrix<3,3> S_gp;

		S_gp.Clear();

		computespin(S_gp,dxdxi_gp,1.0);

		//stress values n and m, Crisfield, Vol. 2, equation (17.76) and (17.78)
		epsilonn(0) *= ym*crosssec_;
		epsilonn(1) *= sm*crosssecshear_;
		epsilonn(2) *= sm*crosssecshear_;

		stressn.Clear();

		stressn.Multiply(Tnew,epsilonn);

		//turning bending strain epsilonm into bending stress stressm
		epsilonm = curvnew_[numgp];
		epsilonm(0) *= sm*Irr_;
		epsilonm(1) *= ym*Iyy_;
		epsilonm(2) *= ym*Izz_;

		stressm.Clear();

		stressm.Multiply(Tnew,epsilonm);

		//computing global internal forces, Crisfield Vol. 2, equation (17.102)
		//note: X = [h,xi(1)*I 0; h(1)S h,xi(1)I;h,xi(2)*I 0; h(2)S h,xi(2)I......] with S = S(dxdx_gp);
		if (force != NULL)
		{
			for (int node=0; node<nnode; ++node)
			{
				for (int j=0; j<3; ++j)
				{
					(*force)(6*node+j) += wgt * deriv(node) * stressn(j);
					(*force)(6*node+3+j) += wgt * deriv(node) * stressm(j);

					for (int i=0; i<3; ++i)
						(*force)(6*node+3+j) += wgt * funct(node) * S_gp(i,j) * stressn(i);
				}
			}
		}//if (force != NULL)


		//computing nonlinear stiffness matrix
		if (stiffmatrix != NULL)
		{

			Kstiff_gp.Clear();
			Ksig1_gp.Clear();
			Ksig2_gp.Clear();

			//setting constitutive parameters , Crisfield, Vol. 2, equation (17.76)
	    Cm(0) = ym*crosssec_;
	    Cm(1) = sm*crosssecshear_;
	    Cm(2) = sm*crosssecshear_;
	    Cb(0) = sm*Irr_;
	    Cb(1) = ym*Iyy_;
	    Cb(2) = ym*Izz_;

			//setting up basis of stiffness matrix according to Crisfield, Vol. 2, equation (17.105)
	    computestiffbasis<nnode>(Tnew,Cm,Cb,S_gp,Kstiff_gp,funct,deriv);

	    Kstiff_gp.Scale(wgt/alpha_[numgp]);

	    //adding nonlinear (stress dependent) parts to tangent stiffness matrix, Crisfield, Vol. 2 equs. (17.106), (17.107) and (17.107b)
	    computeKsig1<nnode>(Ksig1_gp,stressn,stressm,funct,deriv);
	    computeKsig2<nnode>(Ksig2_gp,stressn,S_gp,funct,deriv);	

	    Ksig1_gp.Scale(wgt);
	    Ksig2_gp.Scale(wgt);

	    //shifting values from fixed size matrix to epetra matrix *stiffmatrix
     	for(int i = 0; i < 6*nnode; i++)
     	{
     		for(int j = 0; j < 6*nnode; j++)
     		{
     			(*stiffmatrix)(i,j) += Kstiff_gp(i,j);
     			(*stiffmatrix)(i,j) += Ksig1_gp(i,j);
     			(*stiffmatrix)(i,j) += Ksig2_gp(i,j);
     		}
     	}
	  }//if (stiffmatrix != NULL)
  }//for(int numgp=0; numgp < gausspoints.nquad; numgp++)


  //calculating mass matrix (local version = global version)
  //We use a consistent Timoshenko mass matrix here
  if (massmatrix != NULL)
  {

	  //Get the applied integrationpoints for complete integration
	  DRT::UTILS::IntegrationPoints1D gausspointsmass(MyGaussRule(nnode,gaussexactintegration));

	  //Matrix to store Shape functions as defined throughout the FE lecture
	  LINALG::Matrix<6,6*nnode> N;

	  //Matrix to store mass matrix at each GP in
	  LINALG::Matrix<6*nnode,6*nnode> massmatrixgp;

	  for(int numgp=0; numgp < gausspointsmass.nquad; numgp++)
	  {
		  //Get location and weight of GP in parameter space
		  const double xi = gausspointsmass.qxg[numgp][0];
		  const double wgt = gausspointsmass.qwgt[numgp];

		  //Get h
		  DRT::UTILS::shape_function_1D(funct,xi,distype);

		  //Set N to zeros
		  N.Clear();

		  //Fill up N as defined in the FE lecture
		  for (int node=0; node<nnode; node++)
			  for (int dof=0; dof<6; dof++)
				  N(dof,6*node+dof)=funct(node);

		  //m= density*crosssec* integraloverxi [N_t*N]
     	massmatrixgp.MultiplyTN(density*crosssec_,N,N);

   	  for (int i=0; i<6*nnode; i++)
   	  {
   		  for (int j=0; j<nnode; j++)
   		  {
   			  massmatrixgp(i,6*j+3)= Irr_/crosssec_ * massmatrixgp(i,6*j+3);
   			  massmatrixgp(i,6*j+4)= Irr_/crosssec_ * massmatrixgp(i,6*j+4);
   			  massmatrixgp(i,6*j+5)= Irr_/crosssec_ * massmatrixgp(i,6*j+5);
   			  //This function multiplies all entries associated with
   			  //the rotations. The Irr_ comes from calculations considering
   			  //the moment created by a rotation theta and refers to the assumed direction
   			  //Note: This is an approximation for the rotational values around t2 and t3. For exact
   			  //one would have to differntiate again.
   		  }
   	  }

   	  //Sum up the massmatrices calculated at each gp using the lengthfactor and the weight
   	  for (int i=0; i<6*nnode; i++)
   		  for (int j=0; j<6*nnode; j++)
   			  (*massmatrix)(i,j)+= alphamass_[numgp]*wgt*massmatrixgp(i,j);


	  }//for(int numgp=0; numgp < gausspointsmass.nquad; numgp++)

  }//if (massmatrix != NULL)
  

  /*the following function call applied statistical forces and damping matrix according to the fluctuation dissipation theorem;
   * it is dedicated to the application of beam2 elements in the frame of statistical mechanics problems; for these problems a
   * special vector has to be passed to the element packed in the params parameter list; in case that the control routine calling
   * the element does not attach this special vector to params the following method is just doing nothing, which means that for
   * any ordinary problem of structural mechanics it may be ignored*/
   CalcBrownian<nnode,3,6,4>(params,vel,disp,stiffmatrix,force);


  return;

} // DRT::ELEMENTS::Beam3::b3_nlnstiffmass

/*------------------------------------------------------------------------------------------------------------*
 | lump mass matrix					   (private)                                                   cyron 01/08|
 *------------------------------------------------------------------------------------------------------------*/
template<int nnode>
void DRT::ELEMENTS::Beam3::lumpmass(Epetra_SerialDenseMatrix* emass)
{
  // lump mass matrix
  if (emass != NULL)
  {
	  // we assume #elemat2 is a square matrix
	  for (int c=0; c<(*emass).N(); ++c) // parse columns
	  {
		  double d = 0.0;
		  for (int r=0; r<(*emass).M(); ++r) // parse rows
		  {
			  d += (*emass)(r,c); // accumulate row entries
			  (*emass)(r,c) = 0.0;
		  }

		  (*emass)(c,c) = d; // apply sum of row entries on diagonal
	  }
  }
}

/*-----------------------------------------------------------------------------------------------------------*
 | Evaluate PTC damping (public)                                                                  cyron 10/08|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode>
void DRT::ELEMENTS::Beam3::EvaluatePTC(ParameterList& params,
                                      Epetra_SerialDenseMatrix& elemat1)
{
  /*comment on the calculation of a proper PTC parameter: it is not a priori clear which phyiscal parameters
   * should be incalculated when computing a proper PTC parameter; the reason for instability without PTC is
   * seemingly the vast difference in stiffness of different eigenmodes (namely those related to bending and
   * stretching, respectively). To understand the influence of different parameters to the numerical stability
   * we study a simple model problem: we consider a clamped beam of length L with horizontal and vertical tip
   * loads F_h and F_v, respectively. These tip loads go along horizontal and vertical displacements u and w of
   * the tip point, respectively. Let the beam be undeformed at the beginning of a time step and let the
   * displacements u and w arise at the end of a time step of lenght dt. Then simple calculus shows that with a
   * damping constat gamma ~ eta*L there holds true F_h = EIw/L^3 + gamma*w/dt and F_v = EAu/L + gamma*u/dt.
   * Stability is assumed to be preserved if the ratio between beding u and w is close to one for F_h = F_v.
   * Thus we expect stability if either EI/L^3 ~ EA/L or EI/L^3, EA/L << gamma/dt. In the first case the
   * elastic resistence against bending and streching is comparable, in the second case the problem is dominated
   * by viscous instead of elastic forces. In practice time step size is oriented to either to the bending or
   * stretching time constants tau of the system with dt ~ const1*tau and typically const1 ~ 1e-3. The bending time
   * constants are given by tau_EI = const2*eta*L^4 / EI and the stretching time constants by tau_EA = const3*eta*L^4 / EA,
   * with constant expressions const2, const3 < 1. If dt is chosen according to tau_EI we get gamma /dt ~ const1*const2*EI / L^3,
   * which is always much smaller than the stiffness expression EI/L^3 related to bending and if dt is chosen
   * according to tau_EA the same rationale applies. Therefore EI/L^3, EA/L << gamma/dt never arises in practice
   * and stability depends on the requirement EI/L^3 ~ EA/L. If this requirement is naturally violated an
   * artificial PTC damping has to be employed, which increases the damping stiffness that far that the ratio
   * EI/L^3 ~ EA/L can no longer destabilize the system.
   *
   * The crucial question is obviously how the PTC damping parameter scales with different simulation parameters.
   * In the following we discuss the effect of variations of different parameters:
   *
   * Young's modulus E: As both bending and axial stiffness scale linearly with the Young's modulus E  one may assume
   * that the PTC parameter may be calculated independently on this parameter; this was indeed found in practice:
   * varying E over 3 orders of magnitude upwards and downwards (and scaling time step size by the same factor
   * as all eigenfrequencies depend linearly on Young's modulus) did not affect the PTC parameter required for
   * stabilization.For too small values of E instability was found due to too large curvature in the beam elements,
   * however, this is expected as the beam formulation is valid for moderate curvature only and small values
   * of E naturally admit increasing curvature.
   *
   * Viscosity eta: a similar rationale as for Young's modulus E applies. All the system time constants depend
   * linearly on eta. On the other hand the critical ratio between bending and axial stiffness does not depend
   * on eta. Thus scaling eta and time step size dt by the same factor does not change the PTC factor required
   * for stabilization.
   *
   * Numerical tests revealed that refining the discretization by factor const and at the same time the time
   * step size by a factor const^2 (because the critical axial eigenfrequencies scale with L^2 for element
   * length L) did not change the required PTC parameter. One and the same parameter could be used for a wide
   * range of element lengths up to a scale where the element length became comparable with the persistnece
   * length l_p. For L >= l_p / 2 the simulation became unstable, however, this is supposed to happen not due
   * to an improper PTC parameter, but rather due to the large deformations arsing then, which violated the
   * small strain assumption of this Reissner element. Thus the PTC parameter depends rather on physical parameters
   * than on the choice of the discretization.
   *
   * The above parameter discussion reveals how to adapt the PTC factor in case of changes of the environment of
   * a structure with fixed cross section A, moment of inertia I and length L. However, how to choose the PTC
   * factor and time step size dt for a first discretization and parameter set up has not been discussed so far.
   * Indeed the latter step can be done heuristically once for
   *
   * Cross section A, moment of inertia I: from the above discussed physics one might assume a dependence of the
   * PTC parameter on the ratio of bending and strechting stiffness, i.e. on EI / EA. Such a dependence might
   * considerably exacerbate the application of the PTC algorithm. However, by means of numerical experiments a
   * different rule to deterime the PTC parameter was found: Beyond some ratio EI / EA simulations were found to
   * be unstable without PTC damping. However, a constant PTC damping factor was capable of stabilizing the system
   * over a large range of ratios EI / EA, if the time step size was adopted accordingly. The time step size
   * has to be determined both with respect to bending and stretching time constants. When scaling I by a factor
   * const_I and A by a factor const_A, one first has to decide which of both types of time constants may become
   * critical by the parameter change. Subsequently one has to scale the time step size either by 1/const_A if
   * the stretching time constants are the critical ones or by 1/const_I otherwise.
   *
   *
   * Length L: reduing
   */



  double basisdamp   = (20e-2)*PI*3; //(20e-2)*PI for A = 1.9e-8, (20e-2)*PI*3 for A = 1.9e-6
  double anisofactor = 10;


  
   
  
  //Get the applied integrationpoints for underintegration
  DRT::UTILS::IntegrationPoints1D gausspointsptc(MyGaussRule(nnode,gaussunderintegration));
  //Get discretization typ
  const DRT::Element::DiscretizationType distype = Shape();  
  //matrix to store Ansatzfunktionen
  LINALG::Matrix<1,nnode> funct;
  
  for (int gp=0; gp<gausspointsptc.nquad; gp++)
  {
    
    //Get location and weight of GP in parameter space  
    const double xi = gausspointsptc.qxg[gp][0];
    const double wgt = gausspointsptc.qwgt[gp];
    
    DRT::UTILS::shape_function_1D(funct,xi,distype);
    
    //computing angle increment from current position in comparison with last converged position for damping
    LINALG::Matrix<4,1> deltaQ;
    LINALG::Matrix<3,1> deltatheta;
    quaternionproduct(inversequaternion(Qconv_[gp]),Qnew_[gp],deltaQ);
    quaterniontoangle(deltaQ,deltatheta);
  
    //computing special matrix for anisotropic damping
    LINALG::Matrix<3,3> Tconv;
    LINALG::Matrix<3,3> Theta;
    quaterniontotriad(Qconv_[gp],Tconv);    
    for(int k=0; k<3; k++)
      for(int j = 0; j<3; j++)
        Theta(k,j) = Tconv(k,0)*Tconv(j,0);

    //inverse exponential map
    LINALG::Matrix<3,3> Hinverse;
    Hinverse = Hinv(deltatheta);

    //isotropic artificial stiffness
    LINALG::Matrix<3,3> artstiff;
    artstiff = Hinverse;
    artstiff.Scale(basisdamp);

    //anisotropic artificial stiffness
    LINALG::Matrix<3,3> auxstiff;
    auxstiff.Multiply(Theta,Hinverse);
    auxstiff.Scale(anisofactor*basisdamp);
    artstiff += auxstiff;

    //scale artificial damping with dti parameter for PTC method
    artstiff.Scale( params.get<double>("dti",0.0) );

    for(int i=0; i<nnode; i++)
      for (int j=0; j<nnode; j++)
        for(int k=0; k<3; k++)
          for (int l=0; l<3; l++)
            elemat1(i*6+3+k,j*6+3+l) += artstiff(k,l)*funct(i)*funct(j)*wgt*alpha_[gp];
  }
  
    
  return;
} //DRT::ELEMENTS::Beam3::EvaluatePTC

/*-----------------------------------------------------------------------------------------------------------*
 | computes damping coefficients per lengthand stores them in a matrix in the following order: damping of    |
 | translation parallel to filament axis, damping of translation orthogonal to filament axis, damping of     |
 | rotation around filament axis                                             (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Beam3::MyDampingConstants(ParameterList& params,LINALG::Matrix<3,1>& gamma, const INPAR::STATMECH::FrictionModel& frictionmodel)
{  
  //translational damping coefficients according to Howard, p. 107, table 6.2;
  gamma(0) = 2*PI*params.get<double>("ETA",0.0);
  gamma(1) = 4*PI*params.get<double>("ETA",0.0);
  
  /*damping coefficient of rigid straight rod spinning around its own axis according to Howard, p. 107, table 6.2;
   *as this coefficient is very small for thin rods it is increased artificially by a factor for numerical convencience*/
  double rsquare = pow((4*Iyy_/PI),0.5);
  double artificial = 60*16*2; 
  gamma(2) = 4*PI*params.get<double>("ETA",0.0)*rsquare*artificial;
  
  //in case of an isotropic friction model the same damping coefficients are applied parallel to the polymer axis as perpendicular to it
  if(frictionmodel == INPAR::STATMECH::frictionmodel_isotropicconsistent || frictionmodel == INPAR::STATMECH::frictionmodel_isotropiclumped)
    gamma(0) = gamma(1);


   /* in the following section damping coefficients are replaced by those suggested in Ortega2003, which allows for a 
    * comparison of the finite element simulation with the results of that article; note that we assume that the element
    * length is equivalent to the particle length in the following when computing the length to diameter ratio p*/
   /*
   double lrefe=0;
   for (int gp=0; gp<nnode-1; gp++)
     lrefe += gausspointsdamping.qwgt[gp]*alpha_[gp];
   
   double p=lrefe/(pow(crosssec_*4.0/PI,0.5));
   double Ct=0.312+0.565/p-0.100/pow(p,2.0);
   double Cr=-0.662+0.917/p-0.05/pow(p,2.0);
   gamma(0) = 2.0*PI*params.get<double>("ETA",0.0)/(log(p) + 2*Ct - Cr);
   gamma(1) = 4.0*PI*params.get<double>("ETA",0.0)/(log(p)+Cr);
   gamma(3) = 4.0*PI*params.get<double>("ETA",0.0)*rsquare*artificial*(0.96 + 0.64992/p - 0.17568/p^2);
   */    
}

/*-----------------------------------------------------------------------------------------------------------*
 |computes the number of different random numbers required in each time step for generation of stochastic    |
 |forces;                                                                    (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3::HowManyRandomNumbersINeed()
{
  /*at each Gauss point one needs as many random numbers as randomly excited degrees of freedom, i.e. three
   *random numbers for the translational degrees of freedom and one random number for the rotation around the element axis*/
  return (4*NumNode());

}

/*-----------------------------------------------------------------------------------------------------------*
 |computes velocity of background fluid and gradient of that velocity at a certain evaluation point in       |
 |the physical space                                                         (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int ndim> //number of dimensions of embedding space
void DRT::ELEMENTS::Beam3::MyBackgroundVelocity(ParameterList& params,  //!<parameter list
                                                const LINALG::Matrix<ndim,1>& evaluationpoint,  //!<point at which background velocity and its gradient has to be computed
                                                LINALG::Matrix<ndim,1>& velbackground,  //!< velocity of background fluid
                                                LINALG::Matrix<ndim,ndim>& velbackgroundgrad) //!<gradient of velocity of background fluid
{
  
  /*note: this function is not yet a general one, but always assumes a shear flow, where the velocity of the
   * background fluid is always directed in x-direction. In 3D the velocity increases linearly in z and equals zero for z = 0.
   * In 2D the velocity increases linearly in y and equals zero for y = 0. */
  
  velbackground.PutScalar(0);
  velbackground(0) = evaluationpoint(ndim-1) * params.get<double>("CURRENTSHEAR",0.0);
  
  velbackgroundgrad.PutScalar(0);
  velbackgroundgrad(0,ndim-1) = params.get<double>("CURRENTSHEAR",0.0);

}
/*-----------------------------------------------------------------------------------------------------------*
 | computes rotational damping forces and stiffness (public)                                    cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode> //number of nodes
inline void DRT::ELEMENTS::Beam3::MyRotationalDamping(ParameterList& params,  //!<parameter list
                                              const vector<double>&     vel,  //!< element velocity vector
                                              const vector<double>&     disp, //!<element disp vector
                                              Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseVector* force)//!< element internal force vector
{
  //get time step size
  double dt = params.get<double>("delta time",0.0);
  
  //integration points for underintegration
  DRT::UTILS::IntegrationPoints1D gausspoints(MyGaussRule(nnode,gaussunderintegration));
  
  //get friction model according to which forces and damping are applied
  INPAR::STATMECH::FrictionModel frictionmodel = Teuchos::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL");

  //damping coefficients for translational and rotatinal degrees of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma,frictionmodel);
  

  //matrix to store basis functions evaluated at a certain Gauss point
  LINALG::Matrix<1,nnode> funct;
  
  for (int gp=0; gp<gausspoints.nquad; gp++)//loop through Gauss points
  {    
    //get evaluated basis functions at current Gauss point
    DRT::UTILS::shape_function_1D(funct,gausspoints.qxg[gp][0],Shape());
    
    //rotation between last converged position and current position expressend as a quaternion
    LINALG::Matrix<4,1>  deltaQ;
    quaternionproduct(inversequaternion(Qconv_[gp]),Qnew_[gp],deltaQ);
    
    //rotation between last converged position and current position expressed as a three element rotation vector
    LINALG::Matrix<3,1> deltatheta;
    quaterniontoangle(deltaQ,deltatheta);
    
    //angular velocity at this Gauss point
    LINALG::Matrix<3,1> omega(true);
    omega += deltatheta;
    omega.Scale(1/dt);
    
    //compute matrix T*W*T^t
    LINALG::Matrix<3,3> Tnew;
    LINALG::Matrix<3,3> TWTt;
    quaterniontotriad(Qnew_[gp],Tnew);    
    for(int k=0; k<3; k++)
      for(int j = 0; j<3; j++)
        TWTt(k,j) = Tnew(k,0)*Tnew(j,0);
    
    //compute vector T*W*T^t*\omega 
    LINALG::Matrix<3,1> TWTtomega;
    TWTtomega.Multiply(TWTt,omega);
    
    //compute matrix T*W*T^t*H^(-1)
    LINALG::Matrix<3,3> TWTtHinv;
    TWTtHinv.Multiply(TWTt,Hinv(deltatheta));
    
    //compute spin matrix S(\omega)
    LINALG::Matrix<3,3> Sofomega;
    computespin(Sofomega,omega,1);
    
    //compute matrix T*W*T^t*S(\omega)
    LINALG::Matrix<3,3> TWTtSofomega;
    TWTtSofomega.Multiply(TWTt,Sofomega);
    
    //compute spin matrix S(T*W*T^t*\omega) 
    LINALG::Matrix<3,3> SofTWTtomega;
    computespin(SofTWTtomega,TWTtomega,1);
         
    //loop over all line nodes
    for(int i=0; i<nnode; i++)
      //loop over three dimensions in line direction
      for(int k=0; k<3; k++)
      {
        if(force != NULL)
          (*force)(i*6+3+k) += gamma(2)*TWTtomega(k)*funct(i)*gausspoints.qwgt[gp]*alpha_[gp]; 
        
        if(stiffmatrix != NULL)
          //loop over all column nodes
          for (int j=0; j<nnode; j++)
            //loop over three dimensions in column direction
            for (int l=0; l<3; l++)
              (*stiffmatrix)(i*6+3+k,j*6+3+l) += gamma(2)*( TWTtHinv(k,l) / dt + TWTtSofomega(k,l) - SofTWTtomega(k,l) )*funct(i)*funct(j)*gausspoints.qwgt[gp]*alpha_[gp]; 
      }     
  }
  
  
  
  /*
  
  //neuer Code so umgewandelt, dass explizite Rotationsdämpfung
  
  for (int gp=0; gp<gausspoints.nquad; gp++)//loop through Gauss points
  {    
    //get evaluated basis functions at current Gauss point
    DRT::UTILS::shape_function_1D(funct,gausspoints.qxg[gp][0],Shape());
    
    //rotation between last converged position and current position expressend as a quaternion
    LINALG::Matrix<4,1>  deltaQ;
    quaternionproduct(inversequaternion(Qconv_[gp]),Qnew_[gp],deltaQ);
    
    //rotation between last converged position and current position expressed as a three element rotation vector
    LINALG::Matrix<3,1> deltatheta;
    quaterniontoangle(deltaQ,deltatheta);
    
    //angular velocity at this Gauss point
    LINALG::Matrix<3,1> omega(true);
    omega += deltatheta;
    omega.Scale(1/dt);
    
    //compute matrix T*W*T^t
    LINALG::Matrix<3,3> Tconv;
    LINALG::Matrix<3,3> TWTt;
    quaterniontotriad(Qconv_[gp],Tconv);    
    for(int k=0; k<3; k++)
      for(int j = 0; j<3; j++)
        TWTt(k,j) = Tconv(k,0)*Tconv(j,0);
    
    //compute vector T*W*T^t*\omega 
    LINALG::Matrix<3,1> TWTtomega;
    TWTtomega.Multiply(TWTt,omega);
    
    //compute matrix T*W*T^t*H^(-1)
    LINALG::Matrix<3,3> TWTtHinv;
    TWTtHinv.Multiply(TWTt,Hinv(deltatheta));
      
    //loop over all line nodes
    for(int i=0; i<nnode; i++)
      //loop over three dimensions in line direction
      for(int k=0; k<3; k++)
      {
        if(force != NULL)
          (*force)(i*6+3+k) += gamma(2)*TWTtomega(k)*funct(i)*gausspoints.qwgt[gp]*alpha_[gp]; 
        
        if(stiffmatrix != NULL)
          //loop over all column nodes
          for (int j=0; j<nnode; j++)
            //loop over three dimensions in column direction
            for (int l=0; l<3; l++)
              (*stiffmatrix)(i*6+3+k,j*6+3+l) += gamma(2)*( TWTtHinv(k,l) / dt  )*funct(i)*funct(j)*gausspoints.qwgt[gp]*alpha_[gp]; 
      }     
  }
  
  
  
  */
  
  
  
  
  
  
  
  
 /*

     
    double rsquare = pow((4*Iyy_/PI),0.5);
    //gamma_a artificially increased by factor artificial
    double artificial = 60*16*2;//knyrim Rotationsdämpfung alt:60*16
    double gammaa = 4*PI*params.get<double>("ETA",0.0)*(rsquare)*artificial;
    
    
    //aux variables for rot damp calculation,vector contains matrix for each Gausspoint
    vector<LINALG::Matrix<4,1> > deltaQ;
    deltaQ.resize(nnode-1);
    vector<LINALG::Matrix<3,1> > deltatheta;
    deltatheta.resize(nnode-1);  
    vector<LINALG::Matrix<3,1> > omega;
    omega.resize(nnode-1);
    vector<LINALG::Matrix<3,3> > Tconv;
    Tconv.resize(nnode-1);
    vector<LINALG::Matrix<3,3> >Theta;
    Theta.resize(nnode-1);
    vector<LINALG::Matrix<3,3> > Hinverse;
    Hinverse.resize(nnode-1);
    vector<LINALG::Matrix<3,3> > artstiff;
    artstiff.resize(nnode-1);
    vector<LINALG::Matrix<3,1> > artforce;
    artforce.resize(nnode-1);
    
    const DRT::Element::DiscretizationType distype = Shape(); //Get discretization typ 
   
    for (int gp=0; gp<nnode-1; gp++)//loop through Gauss points
    {
      
      //Get location and weight of GP in parameter space  
      const double xi = gausspoints.qxg[gp][0];
      const double wgt = gausspoints.qwgt[gp];
      
      DRT::UTILS::shape_function_1D(funct,xi,distype);//get evaluated ansatzfunktionen at gausspoints
      
      //computing angle increment from current position in comparison with last converged position for damping
      quaternionproduct(inversequaternion(Qconv_[gp]),Qnew_[gp],deltaQ[gp]);
      quaterniontoangle(deltaQ[gp],deltatheta[gp]);
      
      //angular velocity
      for(int j=0; j<3; j++)
        omega[gp](j,0)=deltatheta[gp](j,0);
      omega[gp].Scale(1/dt);
      
      //computing special matrix for anisotropic damping
      quaterniontotriad(Qconv_[gp],Tconv[gp]);    
      for(int k=0; k<3; k++)
        for(int j = 0; j<3; j++)
          Theta[gp](k,j) = Tconv[gp](k,0)*Tconv[gp](j,0);
      
      //inverse exponential map
      Hinverse[gp]=Hinv(deltatheta[gp]);
    
      //stiffness due to rotational damping
      artstiff[gp].Multiply(Theta[gp],Hinverse[gp]);
      artstiff[gp].Scale(gammaa/dt); 
    
      //forces due to rotational damping  
      artforce[gp].Multiply(Theta[gp],omega[gp]);
      artforce[gp].Scale(gammaa);
    
      
      for(int i=0; i<nnode; i++)//loop twice over nodes (integration N_xi*N_xi)
      {
        //add rot damp forces
        for (int k=0; k<3; k++)
          if(force != NULL)
            (*force)(i*6+3+k) +=artforce[gp](k)*funct(i)*wgt*alpha_[gp];    
        
        for (int j=0; j<nnode; j++)
        {
          for(int k=0; k<3; k++)//loop over rotation dof blocks
            for (int l=0; l<3; l++)
              if(stiffmatrix != NULL)
                (*stiffmatrix)(i*6+3+k,j*6+3+l) += artstiff[gp](k,l)*funct(i)*funct(j)*wgt*alpha_[gp]; 
       
        }
      }   
    }
  
  
 */ 
  
  
  return;
}//DRT::ELEMENTS::Beam3::MyRotationalDamping(.)

/*-----------------------------------------------------------------------------------------------------------*
 | computes translational damping forces and stiffness (public)                                 cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim, int dof> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node
inline void DRT::ELEMENTS::Beam3::MyTranslationalDamping(ParameterList& params,  //!<parameter list
                                                  const vector<double>&     vel,  //!< element velocity vector
                                                  const vector<double>&     disp, //!<element disp vector
                                                  Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                                                  Epetra_SerialDenseVector* force)//!< element internal force vector
{  
  //get time step size
  double dt = params.get<double>("delta time",0.0);
  
  //velocity and gradient of background velocity field
  LINALG::Matrix<ndim,1> velbackground;
  LINALG::Matrix<ndim,ndim> velbackgroundgrad;
  
  //evaluation point in physical space corresponding to a certain Gauss point in parameter space
  LINALG::Matrix<ndim,1> evaluationpoint;
  
  //get friction model according to which forces and damping are applied
  INPAR::STATMECH::FrictionModel frictionmodel = Teuchos::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL");

  //damping coefficients for translational and rotatinal degrees of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma,frictionmodel);
  
  //get vector alpha with Jacobi determinants at each integration point (gets by default those values required for consistent damping matrix)
  vector<double> alpha(alphamass_);
  
  //determine type of numerical integration performed (lumped damping matrix via lobatto integration!)
  IntegrationType integrationtype = gaussexactintegration;
  if(frictionmodel == INPAR::STATMECH::frictionmodel_isotropiclumped)
  {
    integrationtype = lobattointegration;
    alpha = alphanode_;
  }
  
  //get Gauss points and weights for evaluation of damping matrix
  DRT::UTILS::IntegrationPoints1D gausspoints(MyGaussRule(nnode,integrationtype));
  
  //matrix to store basis functions and their derivatives evaluated at a certain Gauss point
  LINALG::Matrix<1,nnode> funct;
  LINALG::Matrix<1,nnode> deriv;

  for(int gp=0; gp < gausspoints.nquad; gp++)
  {    
    //evaluate basis functions and their derivatives at current Gauss point
    DRT::UTILS::shape_function_1D(funct,gausspoints.qxg[gp][0],Shape());
    DRT::UTILS::shape_function_1D_deriv1(deriv,gausspoints.qxg[gp][0],Shape());
     
    //compute point in phyiscal space corresponding to Gauss point
    evaluationpoint.PutScalar(0);
    //loop over all line nodes
    for(int i=0; i<nnode; i++)
      //loop over all dimensions
      for(int j=0; j<ndim; j++)
        evaluationpoint(j) += funct(i)*(Nodes()[i]->X()[j]+disp[dof*i+j]);
    
    //compute velocity and gradient of background flow field at evaluationpoint
    MyBackgroundVelocity<ndim>(params,evaluationpoint,velbackground,velbackgroundgrad);

 
    //compute tangent vector t_{\par} at current Gauss point
    LINALG::Matrix<ndim,1> tpar(true);
    for(int i=0; i<nnode; i++)
      for(int k=0; k<ndim; k++)
        tpar(k) += deriv(i)*(Nodes()[i]->X()[k]+disp[dof*i+k]) / alpha[gp];
    
    //compute velocity vector at this Gauss point
    LINALG::Matrix<ndim,1> velgp(true);
    for(int i=0; i<nnode; i++)
      for(int l=0; l<ndim; l++)
        velgp(l) += funct(i)*vel[dof*i+l]; 
    
    //compute matrix product (t_{\par} \otimes t_{\par}) \cdot velbackgroundgrad
    LINALG::Matrix<ndim,ndim> tpartparvelbackgroundgrad(true);
    for(int i=0; i<ndim; i++)
      for(int j=0; j<ndim; j++)
        for(int k=0; k<ndim; k++)
          tpartparvelbackgroundgrad(i,j) += tpar(i)*tpar(k)*velbackgroundgrad(k,j);
        
    //loop over all line nodes
    for(int i=0; i<nnode; i++)            
      //loop over lines of matrix t_{\par} \otimes t_{\par}
      for(int k=0; k<ndim; k++)
        //loop over columns of matrix t_{\par} \otimes t_{\par}
        for(int l=0; l<ndim; l++)           
        {               
          if(force != NULL)
            (*force)(i*dof+k)+= funct(i)*alpha[gp]*gausspoints.qwgt[gp]*( (k==l)*gamma(1) + (gamma(0) - gamma(1))*tpar(k)*tpar(l) ) *(velgp(l)- velbackground(l));
          
          if(stiffmatrix != NULL)
            //loop over all column nodes
            for (int j=0; j<nnode; j++) 
            {
              (*stiffmatrix)(i*dof+k,j*dof+l) += gausspoints.qwgt[gp]*funct(i)*funct(j)*alpha[gp]*(                 (k==l)*gamma(1) + (gamma(0) - gamma(1))*tpar(k)*tpar(l) ) / dt;
              (*stiffmatrix)(i*dof+k,j*dof+l) -= gausspoints.qwgt[gp]*funct(i)*funct(j)*alpha[gp]*( velbackgroundgrad(k,l)*gamma(1) + (gamma(0) - gamma(1))*tpartparvelbackgroundgrad(k,l) ) ;             
              (*stiffmatrix)(i*dof+k,j*dof+k) += gausspoints.qwgt[gp]*funct(i)*deriv(j)*                                                   (gamma(0) - gamma(1))*tpar(l)*(velgp(l) - velbackground(l));
              (*stiffmatrix)(i*dof+k,j*dof+l) += gausspoints.qwgt[gp]*funct(i)*deriv(j)*                                                   (gamma(0) - gamma(1))*tpar(k)*(velgp(l) - velbackground(l));
            }    
        }   
  }
 
  return;
}//DRT::ELEMENTS::Beam3::MyTranslationalDamping(.)

/*-----------------------------------------------------------------------------------------------------------*
 | computes stochastic forces and resulting stiffness (public)                                  cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
inline void DRT::ELEMENTS::Beam3::MyStochasticForces(ParameterList& params,  //!<parameter list
                                              const vector<double>&     vel,  //!< element velocity vector
                                              const vector<double>&     disp, //!<element disp vector
                                              Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseVector* force)//!< element internal force vector
{
  //get friction model according to which forces and damping are applied
  INPAR::STATMECH::FrictionModel frictionmodel = Teuchos::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL");
  
  //damping coefficients for three translational and one rotatinal degree of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma,frictionmodel);
  

  //get vector alpha with Jacobi determinants at each integration point (gets by default those values required for consistent damping matrix)
  vector<double> alpha(alphamass_);
  
  //determine type of numerical integration performed (lumped damping matrix via lobatto integration!)
  IntegrationType integrationtype = gaussexactintegration;
  if(frictionmodel == INPAR::STATMECH::frictionmodel_isotropiclumped)
  {
    integrationtype = lobattointegration;
    alpha = alphanode_;
  }
  
  //get Gauss points and weights for evaluation of damping matrix
  DRT::UTILS::IntegrationPoints1D gausspoints(MyGaussRule(nnode,integrationtype));
  
  //matrix to store basis functions and their derivatives evaluated at a certain Gauss point
  LINALG::Matrix<1,nnode> funct;
  LINALG::Matrix<1,nnode> deriv;
  
  
  /*get pointer at Epetra multivector in parameter list linking to random numbers for stochastic forces with zero mean
   * and standard deviation (2*kT / dt)^0.5; note carefully: a space between the two subsequal ">" signs is mandatory
   * for the C++ parser in order to avoid confusion with ">>" for streams*/
   RCP<Epetra_MultiVector> randomnumbers = params.get<  RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null);
   


  for(int gp=0; gp < gausspoints.nquad; gp++)
  {
    //evaluate basis functions and their derivatives at current Gauss point
    DRT::UTILS::shape_function_1D(funct,gausspoints.qxg[gp][0],Shape());
    DRT::UTILS::shape_function_1D_deriv1(deriv,gausspoints.qxg[gp][0],Shape());
    
    //compute tangent vector t_{\par} at current Gauss point
    LINALG::Matrix<ndim,1> tpar(true);
    for(int i=0; i<nnode; i++)
      for(int k=0; k<ndim; k++)
        tpar(k) += deriv(i)*(Nodes()[i]->X()[k]+disp[dof*i+k]) / alpha[gp];
     
    
    //loop over all line nodes
    for(int i=0; i<nnode; i++)             
      //loop dimensions with respect to lines
      for(int k=0; k<ndim; k++)
        //loop dimensions with respect to columns
        for(int l=0; l<ndim; l++)           
        {
          if(force != NULL)
            (*force)(i*dof+k) -= funct(i)*(sqrt(gamma(1))*(k==l) + (sqrt(gamma(0)) - sqrt(gamma(1)))*tpar(k)*tpar(l))*(*randomnumbers)[gp*randompergauss+l][LID()]*sqrt(alpha[gp]*gausspoints.qwgt[gp]);          

          if(stiffmatrix != NULL)
            //loop over all column nodes
            for (int j=0; j<nnode; j++) 
            {            
              (*stiffmatrix)(i*dof+k,j*dof+k) -= funct(i)*deriv(j)*tpar(l)*(*randomnumbers)[gp*randompergauss+l][LID()]*sqrt(gausspoints.qwgt[gp]/ alpha[gp])*(sqrt(gamma(0)) - sqrt(gamma(1)));   
              (*stiffmatrix)(i*dof+k,j*dof+l) -= funct(i)*deriv(j)*tpar(k)*(*randomnumbers)[gp*randompergauss+l][LID()]*sqrt(gausspoints.qwgt[gp]/ alpha[gp])*(sqrt(gamma(0)) - sqrt(gamma(1)));  
            }
        }  
  }
  
  
  
  return;
}//DRT::ELEMENTS::Beam3::MyStochasticForces(.)

/*-----------------------------------------------------------------------------------------------------------*
 | computes stochastic moments and (if required) resulting stiffness (public)                   cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int randompergauss> //number of nodes, number of random numbers required per Gauss point, number of random numbers required per Gauss point
inline void DRT::ELEMENTS::Beam3::MyStochasticMoments(ParameterList& params,  //!<parameter list
                                              const vector<double>&     vel,  //!< element velocity vector
                                              const vector<double>&     disp, //!<element disp vector
                                              Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseVector* force)//!< element internal force vector
{

  //get friction model according to which forces and damping are applied
  INPAR::STATMECH::FrictionModel frictionmodel = Teuchos::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL");
  
  //damping coefficients for three translational and one rotatinal degree of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma,frictionmodel);
  
  //determine type of numerical integration performed (note: underintegration applied as for related points triads already known from elasticity)
  IntegrationType integrationtype = gaussunderintegration;

  //get Gauss points and weights for evaluation of damping matrix
  DRT::UTILS::IntegrationPoints1D gausspoints(MyGaussRule(nnode,integrationtype));
  
  //matrix to store basis functions and their derivatives evaluated at a certain Gauss point
  LINALG::Matrix<1,nnode> funct;
   
  /*get pointer at Epetra multivector in parameter list linking to random numbers for stochastic forces with zero mean
   * and standard deviation (2*kT / dt)^0.5; note carefully: a space between the two subsequal ">" signs is mandatory
   * for the C++ parser in order to avoid confusion with ">>" for streams*/
   RCP<Epetra_MultiVector> randomnumbers = params.get<  RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null);

  for(int gp=0; gp < gausspoints.nquad; gp++)
  {
    //evaluate basis functions and their derivatives at current Gauss point
    DRT::UTILS::shape_function_1D(funct,gausspoints.qxg[gp][0],Shape());
    
    //get current triad at this Gauss point:
    LINALG::Matrix<3,3> Tnew;
    quaterniontotriad(Qnew_[gp],Tnew); 
    
    //get first column out of Tnew 
    LINALG::Matrix<3,1> t1;
    for(int i=0; i<3; i++)
      t1(i) = Tnew(i,0);
    
    //compute spin matrix from first column of Tnew times random number
    LINALG::Matrix<3,3> S;
    computespin(S,t1,(*randomnumbers)[gp*randompergauss+3][LID()]);
    
    
    //loop over all line nodes
    for(int i=0; i<nnode; i++)             
      //loop over lines of matrix t_{\par} \otimes t_{\par}
      for(int k=0; k<3; k++)
      {
        if(force != NULL)
          (*force)(i*6+3+k) -= funct(i)*t1(k)*(*randomnumbers)[gp*randompergauss+3][LID()]*sqrt(alpha_[gp]*gausspoints.qwgt[gp]*gamma(2));
      
        if(stiffmatrix != NULL)
          //loop over all column nodes
          for (int j=0; j<nnode; j++) 
            //loop over three dimensions with respect to columns
            for(int l=0; l<3; l++)           
              (*stiffmatrix)(i*6+3+k,j*6+3+l) += funct(i)*funct(j)*S(k,l)*sqrt(alpha_[gp]*gausspoints.qwgt[gp]*gamma(2)); 
               
    }
  }
  return;
}//DRT::ELEMENTS::Beam3::MyStochasticMoments(.)

/*-----------------------------------------------------------------------------------------------------------*
 | Assemble stochastic and viscous forces and respective stiffness according to fluctuation dissipation      |
 | theorem                                                                               (public) cyron 10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
inline void DRT::ELEMENTS::Beam3::CalcBrownian(ParameterList& params,
                                              const vector<double>&           vel,  //!< element velocity vector
                                              const vector<double>&           disp, //!< element displacement vector
                                              Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseVector* force) //!< element internal force vector
{   
  //if no random numbers for generation of stochastic forces are passed to the element no Brownian dynamics calculations are conducted
  if( params.get<  RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null) == Teuchos::null)
    return;
  
  //add stiffness and forces due to translational damping effects
  MyTranslationalDamping<nnode,3,6>(params,vel,disp,stiffmatrix,force); 

  //add stiffness and forces (i.e. moments) due to rotational damping effects
  MyRotationalDamping<nnode>(params,vel,disp,stiffmatrix,force);

  //add stochastic forces and (if required) resulting stiffness
  MyStochasticForces<nnode,3,6,4>(params,vel,disp,stiffmatrix,force);
  
  //add stochastic moments and resulting stiffness
  //MyStochasticMoments<nnode,4>(params,vel,disp,stiffmatrix,force);


return;

}//DRT::ELEMENTS::Beam3::CalcBrownian(.)



#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_BEAM3
