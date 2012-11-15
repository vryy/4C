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



#include "beam3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_fem_general/largerotations.H"
#include "../drt_inpar/inpar_statmech.H"

/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                                 cyron 01/08|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3::Evaluate(Teuchos::ParameterList& params,
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
  else if (action=="calc_struct_energy")        act = Beam3::calc_struct_energy;
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
    case Beam3::calc_struct_energy:
    {
      // need current global displacement and and get them from discretization
      // making use of the local-to-global map lm one can extract current displacemnet and residual values for each degree of freedom
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get state vectors 'displacement'");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

      const int nnode = NumNode();

      switch(nnode)
      {
        case 2:
        {
          b3_energy<2>(params,mydisp,&elevec1);
          break;
        }
        case 3:
        {
          b3_energy<3>(params,mydisp,&elevec1);
          break;
        }
        case 4:
        {
          b3_energy<4>(params,mydisp,&elevec1);
          break;
        }
        case 5:
        {
          b3_energy<5>(params,mydisp,&elevec1);
          break;
        }
        default:
          dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
      }
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
      RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get state vectors 'displacement'");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

      // get residual displacements
      RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (res==null) dserror("Cannot get state vectors 'residual displacement'");
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);

      //only if random numbers for Brownian dynamics are passed to element, get element velocities
      vector<double> myvel(lm.size());
      if( params.get<  RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null) != Teuchos::null)
      {
        RCP<const Epetra_Vector> vel  = discretization.GetState("velocity");
        DRT::UTILS::ExtractMyValues(*vel,myvel,lm);
      }

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
      		vector<double> vel_aux(myvel);
      		vector<double> disp_aux(mydisp);

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
      Qnew_ = Qconv_;
      curvnew_ = curvconv_;
      thetanew_ = thetaconv_;
      thetaprimenew_ = thetaprimeconv_;
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

int DRT::ELEMENTS::Beam3::EvaluateNeumann(Teuchos::ParameterList& params,
                                        DRT::Discretization& discretization,
                                        DRT::Condition& condition,
                                        vector<int>& lm,
                                        Epetra_SerialDenseVector& elevec1,
                                        Epetra_SerialDenseMatrix* elemat1)
{
  // get element displacements
  RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
  if (disp==null) dserror("Cannot get state vector 'displacement'");
  vector<double> mydisp(lm.size());
  DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
  // get element velocities (UNCOMMENT IF NEEDED)
  /*
  RCP<const Epetra_Vector> vel  = discretization.GetState("velocity");
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
    fac = wgt * jacobi_[numgp];

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
 |this function performs an update of the rotation (in quaterion form) at the|
 |numgp-th Gauss point by the incremental rotation deltatheta, by means of a |
 |quaternion product and then computes the respective new triad Tnew at the  |
 | Gauss point							                              (public) cyron02/09|
 *---------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Beam3::updatetriad(const LINALG::Matrix<3,1>& deltatheta, LINALG::Matrix<3,3>& Tnew, const int numgp)
{
  //computing quaternion equivalent to rotation by deltatheta
  LINALG::Matrix<4,1> Qrot;
  LARGEROTATIONS::angletoquaternion(deltatheta,Qrot);

  //computing quaternion Qnew_ for new configuration of Qold_ for old configuration by means of a quaternion product
  LARGEROTATIONS::quaternionproduct(Qold_[numgp],Qrot,Qnew_[numgp]);

  //normalizing quaternion in order to make sure that it keeps unit absolute values through time stepping
  double abs = Qnew_[numgp].Norm2();
  for(int i = 0; i<4; i++)
    Qnew_[numgp](i) = Qnew_[numgp](i) / abs;

  LARGEROTATIONS::quaterniontotriad(Qnew_[numgp],Tnew);

  return;
} //DRT::ELEMENTS::Beam3::updatetriad


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
     *  LARGEROTATIONS::angletoquaternion(deltatheta,q);
     *  LARGEROTATIONS::quaterniontoangle(q,deltatheta);
     *
     * in the very beginning of this method. However, for the above reason we assume that theta lies always in the proper
     * region and thus save the related compuational cost and just throw an error if this prerequesite is unexpectedly not
     * satisfied */
    if(abs_theta > M_PI)
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
  LARGEROTATIONS::quaterniontotriad(Qold_[numgp],Told);

  //compute spin matrix from eq. (17.73)
  LINALG::Matrix<3,3> spin;
  LARGEROTATIONS::computespin(spin, deltatheta);

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
  LARGEROTATIONS::computespin(Sn,stressn);
  LARGEROTATIONS::computespin(Sm,stressm);

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
  LARGEROTATIONS::computespin(Sn,stressn);
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
 | calculation of elastic energy (private)                                                         cyron 12/10|
 *-----------------------------------------------------------------------------------------------------------*/
template<int nnode>
void DRT::ELEMENTS::Beam3::b3_energy( Teuchos::ParameterList& params,
                                      vector<double>& disp,
                                      Epetra_SerialDenseVector* intenergy)
{
  //initialize energies (only one kind of energy computed here
  (*intenergy)(0) = 0.0;

  //constitutive laws from Crisfield, Vol. 2, equation (17.76)
  LINALG::Matrix<3,1> Cm;
  LINALG::Matrix<3,1> Cb;

  //normal/shear strain and bending strain(curvature)
  LINALG::Matrix<3,1> epsilonn;
  LINALG::Matrix<3,1> epsilonm;

  //derivative of x with respect to xi
  LINALG::Matrix<3,1> dxdxi_gp;

  //triad at GP, Crisfiel Vol. 2, equation (17.73)
  LINALG::Matrix<3,3> Tnew;

  //first of all we get the material law
  Teuchos::RCP<const MAT::Material> currmat = Material();
  double ym = 0;
  double sm = 0;

  //assignment of material parameters; only St.Venant material is accepted for this beam
  switch(currmat->MaterialType())
  {
    case INPAR::MAT::m_stvenant:// only linear elastic material supported
    {
      const MAT::StVenantKirchhoff* actmat = static_cast<const MAT::StVenantKirchhoff*>(currmat.get());
      ym = actmat->Youngs();
      sm = ym / (2*(1 + actmat->PoissonRatio()));
    }
    break;
    default:
    dserror("unknown or improper type of material law");
  }

  /*first displacement vector is modified for proper element evaluation in case of periodic boundary conditions; in case that
   *no periodic boundary conditions are to be applied the following code line may be ignored or deleted*/
  NodeShift<nnode,3>(params,disp);

  //Get integrationpoints for underintegration
  DRT::UTILS::IntegrationPoints1D gausspoints(MyGaussRule(nnode,gaussunderintegration));

  //Get DiscretizationType
  const DRT::Element::DiscretizationType distype = Shape();

  //Matrices for h and h,xi
  LINALG::Matrix<1,nnode> deriv;

  //Loop through all GP and calculate their contribution to the forcevector and stiffnessmatrix
  for(int numgp=0; numgp < gausspoints.nquad; numgp++)
  {
    //Get location and weight of GP in parameter space
    const double xi = gausspoints.qxg[numgp][0];
    const double wgt = gausspoints.qwgt[numgp];

    //Get h,xi
    DRT::UTILS::shape_function_1D_deriv1(deriv,xi,distype);

    //set up current dxdxi, theta, and thetaprime at the GP
    dxdxi_gp.Clear();
    for (int dof=0; dof<3; ++dof)//j
      for (int node=0; node<nnode; ++node)
         dxdxi_gp(dof) += (Nodes()[node]->X()[dof]+disp[6*node+dof])*deriv(node);


    //compute current triad at numgp-th Gauss point
    LARGEROTATIONS::quaterniontotriad(Qnew_[numgp],Tnew);

    //setting constitutive parameters , Crisfield, Vol. 2, equation (17.76)
    Cm(0) = ym*crosssec_;
    Cm(1) = sm*crosssecshear_;
    Cm(2) = sm*crosssecshear_;
    Cb(0) = sm*Irr_;
    Cb(1) = ym*Iyy_;
    Cb(2) = ym*Izz_;

    //computing current axial and shear strain epsilon, Crisfield, Vol. 2, equation (17.97)
    epsilonn.Clear();
    epsilonn.MultiplyTN(Tnew,dxdxi_gp);
    epsilonn.Scale(1/jacobi_[numgp]);
    epsilonn(0) -=  1.0;

    epsilonm.Clear();
    epsilonm = curvnew_[numgp];

    //adding elastic energy from epsilonn at this Gauss point
    for(int i=0; i<3; i++)
    {
      (*intenergy)(0) += 0.5*epsilonn(i)*epsilonn(i)*Cm(i)*wgt*jacobi_[numgp];
      (*intenergy)(0) += 0.5*epsilonm(i)*epsilonm(i)*Cb(i)*wgt*jacobi_[numgp];
    }

  }

  // hack in order to just get the contribution of beam3II elements (filaments) -> set contribution to 0
  //(*intenergy)(0) = 0.0;

  return;

} // DRT::ELEMENTS::Beam3::b3_energy


/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   cyron 01/08|
 *-----------------------------------------------------------------------------------------------------------*/
template<int nnode>
void DRT::ELEMENTS::Beam3::b3_nlnstiffmass( Teuchos::ParameterList& params,
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


  //"new" variables have to be adopted to current discplacement

  /*first displacement vector is modified for proper element evaluation in case of periodic boundary conditions; in case that
   *no periodic boundary conditions are to be applied the following code line may be ignored or deleted*/
  NodeShift<nnode,3>(params,disp);

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
				thetaprimenew_gp(dof)     += disp[6*node+3+dof]*deriv(node)/jacobi_[numgp];

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

		//updating local curvature according
		//updatecurvature(Tnew,deltatheta_gp,deltathetaprime_gp,numgp);
		approxupdatecurvature(Tnew,deltatheta_gp,deltathetaprime_gp,numgp);

		epsilonn.Clear();

		//computing current axial and shear strain epsilon, Crisfield, Vol. 2, equation (17.97)
		epsilonn.MultiplyTN(Tnew,dxdxi_gp);
		epsilonn.Scale(1/jacobi_[numgp]);
		epsilonn(0) -=  1.0;

		epsilonm.Clear();

		//computing spin matrix S(dxdxi_gp) according to Crisfield, Vol. 2, equation (17.100)
		LINALG::Matrix<3,3> S_gp;

		S_gp.Clear();

		 LARGEROTATIONS::computespin(S_gp,dxdxi_gp);

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

	    Kstiff_gp.Scale(wgt/jacobi_[numgp]);

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
   			  (*massmatrix)(i,j)+= jacobimass_[numgp]*wgt*massmatrixgp(i,j);


	  }//for(int numgp=0; numgp < gausspointsmass.nquad; numgp++)

  }//if (massmatrix != NULL)

  /*the following function call applied statistical forces and damping matrix according to the fluctuation dissipation theorem;
   * it is dedicated to the application of beam2 elements in the frame of statistical mechanics problems; for these problems a
   * special vector has to be passed to the element packed in the params parameter list; in case that the control routine calling
   * the element does not attach this special vector to params the following method is just doing nothing, which means that for
   * any ordinary problem of structural mechanics it may be ignored*/
  CalcBrownian<nnode,3,6,4>(params,vel,disp,stiffmatrix,force);

  // in statistical mechanics simulations, a deletion influenced by the values of the internal force vector might occur
  if(params.get<string>("internalforces","no")=="yes" && force != NULL)
  	internalforces_ = Teuchos::rcp(new Epetra_SerialDenseVector(*force));

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
void DRT::ELEMENTS::Beam3::EvaluatePTC(Teuchos::ParameterList& params,
                                      Epetra_SerialDenseMatrix& elemat1)
{
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
    LARGEROTATIONS::quaternionproduct(LARGEROTATIONS::inversequaternion(Qconv_[gp]),Qnew_[gp],deltaQ);
    LINALG::Matrix<3,1> deltatheta;
    LARGEROTATIONS::quaterniontoangle(deltaQ,deltatheta);

    //isotropic artificial stiffness
    LINALG::Matrix<3,3> artstiff;
    artstiff = LARGEROTATIONS::Tmatrix(deltatheta);

    //scale artificial damping with crotptc parameter for PTC method
    artstiff.Scale( params.get<double>("crotptc",0.0) );

    for(int i=0; i<nnode; i++)
      for (int j=0; j<nnode; j++)
        for(int k=0; k<3; k++)
          for (int l=0; l<3; l++)
            elemat1(i*6+3+k,j*6+3+l) += artstiff(k,l)*funct(i)*funct(j)*wgt*jacobi_[gp];

    //PTC for translational degrees of freedom
    for(int i=0; i<nnode; i++)
      for (int j=0; j<nnode; j++)
        for(int k=0; k<3; k++)
            elemat1(i*6+k,j*6+k) += params.get<double>("ctransptc",0.0)*funct(i)*funct(j)*wgt*jacobi_[gp];

  }

  return;
} //DRT::ELEMENTS::Beam3::EvaluatePTC

/*-----------------------------------------------------------------------------------------------------------*
 | computes damping coefficients per lengthand stores them in a matrix in the following order: damping of    |
 | translation parallel to filament axis, damping of translation orthogonal to filament axis, damping of     |
 | rotation around filament axis                                             (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Beam3::MyDampingConstants(Teuchos::ParameterList& params,LINALG::Matrix<3,1>& gamma, const INPAR::STATMECH::FrictionModel& frictionmodel)
{
  //translational damping coefficients according to Howard, p. 107, table 6.2;
  gamma(0) = 2*PI*params.get<double>("ETA",0.0);
  gamma(1) = 4*PI*params.get<double>("ETA",0.0);

  /*damping coefficient of rigid straight rod spinning around its own axis according to Howard, p. 107, table 6.2;
   *as this coefficient is very small for thin rods it is increased artificially by a factor for numerical convencience*/
  double rsquare = std::pow((4*Iyy_/PI),0.5);
  double artificial = 4000;//1000;  //1000 not bad for standard Actin3D_10.dat files; for 40 elements also 1 seems to work really well; for large networks 4000 seems good (artificial contribution then still just ~0.1 % of nodal moments)
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
     lrefe += gausspointsdamping.qwgt[gp]*jacobi_[gp];

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
void DRT::ELEMENTS::Beam3::MyBackgroundVelocity(Teuchos::ParameterList& params,  //!<parameter list
                                                const LINALG::Matrix<ndim,1>& evaluationpoint,  //!<point at which background velocity and its gradient has to be computed
                                                LINALG::Matrix<ndim,1>& velbackground,  //!< velocity of background fluid
                                                LINALG::Matrix<ndim,ndim>& velbackgroundgrad) //!<gradient of velocity of background fluid
{

  /*note: this function is not yet a general one, but always assumes a shear flow, where the velocity of the
   * background fluid is always directed in direction params.get<int>("OSCILLDIR",0) and orthogonal to z-axis.
   * In 3D the velocity increases linearly in z and equals zero for z = 0.
   * In 2D the velocity increases linearly in y and equals zero for y = 0. */

  //velocity at upper boundary of domain
  double uppervel = 0.0;

  //default values for background velocity and its gradient
  velbackground.PutScalar(0);
  velbackgroundgrad.PutScalar(0);

  double time = params.get<double>("total time",0.0);
  double starttime = params.get<double>("STARTTIMEACT",0.0);
  double dt = params.get<double>("delta time");
  double shearamplitude = params.get<double> ("SHEARAMPLITUDE", 0.0);
  int curvenumber = params.get<int> ("CURVENUMBER", -1);
  int oscilldir = params.get<int> ("OSCILLDIR", -1);
  INPAR::STATMECH::DBCType dbctype = params.get<INPAR::STATMECH::DBCType>("DBCTYPE", INPAR::STATMECH::dbctype_std);
  Teuchos::RCP<std::vector<double> > defvalues = Teuchos::rcp(new std::vector<double>(3,0.0));
  Teuchos::RCP<std::vector<double> > periodlength = params.get("PERIODLENGTH", defvalues);

  bool shearflow = false;
  if(dbctype==INPAR::STATMECH::dbctype_shearfixed || dbctype==INPAR::STATMECH::dbctype_sheartrans)
    shearflow = true;
  //oscillations start only at params.get<double>("STARTTIMEACT",0.0)
  if(periodlength->at(0) > 0.0)
    if(shearflow && time>starttime && fabs(time-starttime)>dt/1e4 && curvenumber >=  1 && oscilldir >= 0 )
    {
      uppervel = shearamplitude * (DRT::Problem::Instance()->Curve(curvenumber-1).FctDer(time,1))[1];

      //compute background velocity
      velbackground(oscilldir) = (evaluationpoint(ndim-1) / periodlength->at(ndim-1)) * uppervel;

      //compute gradient of background velocity
      velbackgroundgrad(oscilldir,ndim-1) = uppervel / periodlength->at(ndim-1);
    }

}
/*-----------------------------------------------------------------------------------------------------------*
 | computes rotational damping forces and stiffness (public)                                    cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode> //number of nodes
inline void DRT::ELEMENTS::Beam3::MyRotationalDamping(Teuchos::ParameterList& params,  //!<parameter list
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
  INPAR::STATMECH::FrictionModel frictionmodel = DRT::INPUT::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL");

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
    LARGEROTATIONS::quaternionproduct(LARGEROTATIONS::inversequaternion(Qconv_[gp]),Qnew_[gp],deltaQ);

    //rotation between last converged position and current position expressed as a three element rotation vector
    LINALG::Matrix<3,1> deltatheta;
    LARGEROTATIONS::quaterniontoangle(deltaQ,deltatheta);

    //angular velocity at this Gauss point
    LINALG::Matrix<3,1> omega(true);
    omega += deltatheta;
    omega.Scale(1/dt);

    //compute matrix T*W*T^t
    LINALG::Matrix<3,3> Tnew;
    LINALG::Matrix<3,3> TWTt;
    LARGEROTATIONS::quaterniontotriad(Qnew_[gp],Tnew);
    for(int k=0; k<3; k++)
      for(int j = 0; j<3; j++)
        TWTt(k,j) = Tnew(k,0)*Tnew(j,0);

    //compute vector T*W*T^t*\omega
    LINALG::Matrix<3,1> TWTtomega;
    TWTtomega.Multiply(TWTt,omega);

    //compute matrix T*W*T^t*H^(-1)
    LINALG::Matrix<3,3> TWTtHinv;
    TWTtHinv.Multiply(TWTt,LARGEROTATIONS::Tmatrix(deltatheta));

    //compute spin matrix S(\omega)
    LINALG::Matrix<3,3> Sofomega;
    LARGEROTATIONS::computespin(Sofomega,omega);

    //compute matrix T*W*T^t*S(\omega)
    LINALG::Matrix<3,3> TWTtSofomega;
    TWTtSofomega.Multiply(TWTt,Sofomega);

    //compute spin matrix S(T*W*T^t*\omega)
    LINALG::Matrix<3,3> SofTWTtomega;
    LARGEROTATIONS::computespin(SofTWTtomega,TWTtomega);

    //loop over all line nodes
    for(int i=0; i<nnode; i++)
      //loop over three dimensions in line direction
      for(int k=0; k<3; k++)
      {
        if(force != NULL)
          (*force)(i*6+3+k) += gamma(2)*TWTtomega(k)*funct(i)*gausspoints.qwgt[gp]*jacobi_[gp];

        if(stiffmatrix != NULL)
          //loop over all column nodes
          for (int j=0; j<nnode; j++)
            //loop over three dimensions in column direction
            for (int l=0; l<3; l++)
              (*stiffmatrix)(i*6+3+k,j*6+3+l) += gamma(2)*( TWTtHinv(k,l) / dt + TWTtSofomega(k,l) - SofTWTtomega(k,l) )*funct(i)*funct(j)*gausspoints.qwgt[gp]*jacobi_[gp];
      }
  }


  return;
}//DRT::ELEMENTS::Beam3::MyRotationalDamping(.)

/*-----------------------------------------------------------------------------------------------------------*
 | computes translational damping forces and stiffness (public)                                 cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim, int dof> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node
inline void DRT::ELEMENTS::Beam3::MyTranslationalDamping(Teuchos::ParameterList& params,  //!<parameter list
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
  INPAR::STATMECH::FrictionModel frictionmodel = DRT::INPUT::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL");

  //damping coefficients for translational and rotatinal degrees of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma,frictionmodel);

  //get vector jacobi with Jacobi determinants at each integration point (gets by default those values required for consistent damping matrix)
  vector<double> jacobi(jacobimass_);

  //determine type of numerical integration performed (lumped damping matrix via lobatto integration!)
  IntegrationType integrationtype = gaussexactintegration;
  if(frictionmodel == INPAR::STATMECH::frictionmodel_isotropiclumped)
  {
    integrationtype = lobattointegration;
    jacobi = jacobinode_;
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
        tpar(k) += deriv(i)*(Nodes()[i]->X()[k]+disp[dof*i+k]) / jacobi[gp];

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
            (*force)(i*dof+k)+= funct(i)*jacobi[gp]*gausspoints.qwgt[gp]*( (k==l)*gamma(1) + (gamma(0) - gamma(1))*tpar(k)*tpar(l) ) *(velgp(l)- velbackground(l));

          if(stiffmatrix != NULL)
            //loop over all column nodes
            for (int j=0; j<nnode; j++)
            {
              (*stiffmatrix)(i*dof+k,j*dof+l) += gausspoints.qwgt[gp]*funct(i)*funct(j)*jacobi[gp]*(                 (k==l)*gamma(1) + (gamma(0) - gamma(1))*tpar(k)*tpar(l) ) / dt;
              (*stiffmatrix)(i*dof+k,j*dof+l) -= gausspoints.qwgt[gp]*funct(i)*funct(j)*jacobi[gp]*( velbackgroundgrad(k,l)*gamma(1) + (gamma(0) - gamma(1))*tpartparvelbackgroundgrad(k,l) ) ;
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
inline void DRT::ELEMENTS::Beam3::MyStochasticForces(Teuchos::ParameterList& params,  //!<parameter list
                                              const vector<double>&     vel,  //!< element velocity vector
                                              const vector<double>&     disp, //!<element disp vector
                                              Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseVector* force)//!< element internal force vector
{
  //get friction model according to which forces and damping are applied
  INPAR::STATMECH::FrictionModel frictionmodel = DRT::INPUT::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL");

  //damping coefficients for three translational and one rotatinal degree of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma,frictionmodel);


  //get vector jacobi with Jacobi determinants at each integration point (gets by default those values required for consistent damping matrix)
  vector<double> jacobi(jacobimass_);

  //determine type of numerical integration performed (lumped damping matrix via lobatto integration!)
  IntegrationType integrationtype = gaussexactintegration;
  if(frictionmodel == INPAR::STATMECH::frictionmodel_isotropiclumped)
  {
    integrationtype = lobattointegration;
    jacobi = jacobinode_;
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
        tpar(k) += deriv(i)*(Nodes()[i]->X()[k]+disp[dof*i+k]) / jacobi[gp];


    //loop over all line nodes
    for(int i=0; i<nnode; i++)
      //loop dimensions with respect to lines
      for(int k=0; k<ndim; k++)
        //loop dimensions with respect to columns
        for(int l=0; l<ndim; l++)
        {
          if(force != NULL)
            (*force)(i*dof+k) -= funct(i)*(sqrt(gamma(1))*(k==l) + (sqrt(gamma(0)) - sqrt(gamma(1)))*tpar(k)*tpar(l))*(*randomnumbers)[gp*randompergauss+l][LID()]*sqrt(jacobi[gp]*gausspoints.qwgt[gp]);

          if(stiffmatrix != NULL)
            //loop over all column nodes
            for (int j=0; j<nnode; j++)
            {
              (*stiffmatrix)(i*dof+k,j*dof+k) -= funct(i)*deriv(j)*tpar(l)*(*randomnumbers)[gp*randompergauss+l][LID()]*sqrt(gausspoints.qwgt[gp]/ jacobi[gp])*(sqrt(gamma(0)) - sqrt(gamma(1)));
              (*stiffmatrix)(i*dof+k,j*dof+l) -= funct(i)*deriv(j)*tpar(k)*(*randomnumbers)[gp*randompergauss+l][LID()]*sqrt(gausspoints.qwgt[gp]/ jacobi[gp])*(sqrt(gamma(0)) - sqrt(gamma(1)));
            }
        }
  }



  return;
}//DRT::ELEMENTS::Beam3::MyStochasticForces(.)

/*-----------------------------------------------------------------------------------------------------------*
 | computes stochastic moments and (if required) resulting stiffness (public)                   cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int randompergauss> //number of nodes, number of random numbers required per Gauss point, number of random numbers required per Gauss point
inline void DRT::ELEMENTS::Beam3::MyStochasticMoments(Teuchos::ParameterList& params,  //!<parameter list
                                              const vector<double>&     vel,  //!< element velocity vector
                                              const vector<double>&     disp, //!<element disp vector
                                              Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseVector* force)//!< element internal force vector
{

  //get friction model according to which forces and damping are applied
  INPAR::STATMECH::FrictionModel frictionmodel = DRT::INPUT::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL");

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
    LARGEROTATIONS::quaterniontotriad(Qnew_[gp],Tnew);

    //get first column out of Tnew
    LINALG::Matrix<3,1> t1;
    for(int i=0; i<3; i++)
      t1(i) = Tnew(i,0);

    //compute spin matrix from first column of Tnew times random number
    LINALG::Matrix<3,3> S;
    LARGEROTATIONS::computespin(S,t1);
    S.Scale((*randomnumbers)[gp*randompergauss+3][LID()]);


    //loop over all line nodes
    for(int i=0; i<nnode; i++)
      //loop over lines of matrix t_{\par} \otimes t_{\par}
      for(int k=0; k<3; k++)
      {
        if(force != NULL)
          (*force)(i*6+3+k) -= funct(i)*t1(k)*(*randomnumbers)[gp*randompergauss+3][LID()]*sqrt(jacobi_[gp]*gausspoints.qwgt[gp]*gamma(2));

        if(stiffmatrix != NULL)
          //loop over all column nodes
          for (int j=0; j<nnode; j++)
            //loop over three dimensions with respect to columns
            for(int l=0; l<3; l++)
              (*stiffmatrix)(i*6+3+k,j*6+3+l) += funct(i)*funct(j)*S(k,l)*sqrt(jacobi_[gp]*gausspoints.qwgt[gp]*gamma(2));

    }
  }
  return;
}//DRT::ELEMENTS::Beam3::MyStochasticMoments(.)

/*-----------------------------------------------------------------------------------------------------------*
 | Assemble stochastic and viscous forces and respective stiffness according to fluctuation dissipation      |
 | theorem                                                                               (public) cyron 10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
inline void DRT::ELEMENTS::Beam3::CalcBrownian(Teuchos::ParameterList& params,
                                              const vector<double>&           vel,  //!< element velocity vector
                                              const vector<double>&           disp, //!< element displacement vector
                                              Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseVector* force) //!< element internal force vector
{
  //if no random numbers for generation of stochastic forces are passed to the element no Brownian dynamics calculations are conducted
  if( params.get<  RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null) == Teuchos::null)
    return;

  //add stiffness and forces due to translational damping effects
  MyTranslationalDamping<nnode,ndim,dof>(params,vel,disp,stiffmatrix,force);

  //add stiffness and forces (i.e. moments) due to rotational damping effects
  MyRotationalDamping<nnode>(params,vel,disp,stiffmatrix,force);

  //add stochastic forces and (if required) resulting stiffness
  MyStochasticForces<nnode,ndim,dof,randompergauss>(params,vel,disp,stiffmatrix,force);

  //add stochastic moments and resulting stiffness
  MyStochasticMoments<nnode,randompergauss>(params,vel,disp,stiffmatrix,force);


return;

}//DRT::ELEMENTS::Beam3::CalcBrownian(.)

/*-----------------------------------------------------------------------------------------------------------*
 | shifts nodes so that proper evaluation is possible even in case of periodic boundary conditions; if two   |
 | nodes within one element are separated by a periodic boundary, one of them is shifted such that the final |
 | distance in R^3 is the same as the initial distance in the periodic space; the shift affects computation  |
 | on element level within that very iteration step, only (no change in global variables performed)          |                                 |
 |                                                                                       (public) cyron 10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim> //number of nodes, number of dimensions
inline void DRT::ELEMENTS::Beam3::NodeShift(Teuchos::ParameterList& params,  //!<parameter list
                                            vector<double>& disp) //!<element disp vector
{
  /*get number of degrees of freedom per node; note: the following function assumes the same number of degrees
   *of freedom for each element node*/
  int numdof = NumDofPerNode(*(Nodes()[0]));

  double time = params.get<double>("total time",0.0);
	double starttime = params.get<double>("STARTTIMEACT",0.0);
	double dt = params.get<double>("delta time");
  double shearamplitude = params.get<double> ("SHEARAMPLITUDE", 0.0);
  int curvenumber = params.get<int> ("CURVENUMBER", -1);
  int oscilldir = params.get<int> ("OSCILLDIR", -1);
  Teuchos::RCP<std::vector<double> > defvalues = Teuchos::rcp(new std::vector<double>(3,0.0));
  Teuchos::RCP<std::vector<double> > periodlength = params.get("PERIODLENGTH", defvalues);

  INPAR::STATMECH::DBCType dbctype = params.get<INPAR::STATMECH::DBCType>("DBCTYPE", INPAR::STATMECH::dbctype_std);
  bool shearflow = false;
  if(dbctype==INPAR::STATMECH::dbctype_shearfixed || dbctype==INPAR::STATMECH::dbctype_sheartrans)
    shearflow = true;

  /*only if periodic boundary conditions are in use, i.e. params.get<double>("PeriodLength",0.0) > 0.0, this
   * method has to change the displacement variables*/
  if(periodlength->at(0) > 0.0)
    //loop through all nodes except for the first node which remains fixed as reference node
    for(int i=1;i<nnode;i++)
    {
      for(int dof= ndim - 1; dof > -1; dof--)
      {
        /*if the distance in some coordinate direction between some node and the first node becomes smaller by adding or subtracting
         * the period length, the respective node has obviously been shifted due to periodic boundary conditions and should be shifted
         * back for evaluation of element matrices and vectors; this way of detecting shifted nodes works as long as the element length
         * is smaller than half the periodic length*/
        if( fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) + periodlength->at(dof) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) < fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) )
        {
          disp[numdof*i+dof] += periodlength->at(dof);

          /*the upper domain surface orthogonal to the z-direction may be subject to shear Dirichlet boundary condition; the lower surface
           *may be fixed by DBC. To avoid problmes when nodes exit the domain through the upper z-surface and reenter through the lower
           *z-surface, the shear has to be substracted from nodal coordinates in that case */
          if(shearflow && dof == 2 && curvenumber >=  1 && time>starttime && fabs(time-starttime)>dt/1e4)
            disp[numdof*i+oscilldir] += shearamplitude*DRT::Problem::Instance()->Curve(curvenumber-1).f(time);
        }

        if( fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - periodlength->at(dof) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) < fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) )
        {
          disp[numdof*i+dof] -= periodlength->at(dof);

          /*the upper domain surface orthogonal to the z-direction may be subject to shear Dirichlet boundary condition; the lower surface
           *may be fixed by DBC. To avoid problmes when nodes exit the domain through the lower z-surface and reenter through the upper
           *z-surface, the shear has to be added to nodal coordinates in that case */
          if(shearflow && dof == 2 && curvenumber >=  1 && time>starttime && fabs(time-starttime)>dt/1e4)
            disp[numdof*i+oscilldir] -= shearamplitude*DRT::Problem::Instance()->Curve(curvenumber-1).f(time);
        }
      }
    }

return;

}//DRT::ELEMENTS::Beam3::NodeShift



