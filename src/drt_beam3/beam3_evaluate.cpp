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
#include "../drt_lib/linalg_fixedsizematrix.H"
//including random number library of blitz for statistical forces
#include <random/normal.h>
#include "../drt_mat/stvenantkirchhoff.H"

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
  else if (action=="calc_brownian")        act = Beam3::calc_brownian;
  else if (action=="calc_struct_ptcstiff")        act = Beam3::calc_struct_ptcstiff;
  else dserror("Unknown type of action for Beam3");

  switch(act)
  {
    case Beam3::calc_struct_ptcstiff:
    {
      EvaluatePTC(params, elemat1);
    }
    break;
    //action type for evaluating statistical forces
    case Beam3::calc_brownian:
    {

    	/*
      // get element displacements (for use in shear flow fields)
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get state vector 'displacement'");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      */

      /*in case of parallel computing the random forces have to be handled in a special way: normally
       * in the frame of evaluation each processor evaluates forces for its column elements (including
       * ghost elements); later on in the assembly each processor adds the thereby gained forces only
       * for those DOF of whose owner it is. In case of random forces such an assembly would render it
       * impossible to establish certain correlations between forces related to nodes with different ownerss
       * (for each nodes the random forces would be evaluated in an identical process, but due to
       * independent random numbers); as correlation between forces is restricted to the support of at
       * the maximum one element a solution to this problem is, to evaluate all the forces of one element
       * only by means of one specific processor (here we employ the elemnet owner processor); these
       * forces are assembled in a column map vector and later exported to a row map force vector; this
       * export is carried out additively so that it is important not to evaluate any forces at all if
       * this processor is not owner of the element;
       * note: the crucial difference between this assembly and the common one is that for certain nodal
       * forces not the owner of the node is responsible, but the owner of the element*/

      //test whether this processor is row map owner of the element (otherwise no forces added)
     // if(this->Owner() != discretization.Comm().MyPID()) return 0;


      //compute stochastic forces in local frame
      ComputeLocalBrownianForces(params);

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


/*    //the following code block can be used to check quickly whether the nonlinear stiffness matrix is calculated
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
      		vel_aux[6*k + i] += h_rel * params.get<double>("gamma",0.581) / ( params.get<double>("delta time",0.01)*params.get<double>("beta",0.292) );
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


      		for(int u = 0 ; u < 6*nnode ; u++ )
      		{
      			stiff_approx(u,k*6+i)= ( pow(force_aux[u],2) - pow(elevec1(u),2) )/ (h_rel * (force_aux[u] + elevec1(u) ) );// berechnet mit dem Differenzenquotient dFres/du
      		}

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
                                        Epetra_SerialDenseVector& elevec1)
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
  curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);

  // no. of nodes on this element; the following line is only valid for elements with constant number of
  // degrees of freedom per node
  const int numdf = 6;
  const DiscretizationType distype = this->Shape();

  // gaussian points
  const DRT::UTILS::IntegrationPoints1D intpoints(gaussrule_);

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
 | Compute forces for Brownian motion (public)                                                       06/09|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3::ComputeLocalBrownianForces(ParameterList& params)
{
  /*creating a random generator object which creates random numbers with mean = 0 and standard deviation
   * (2kT/dt)^0,5 with thermal energy kT, time step size dt; using Blitz namespace "ranlib" for random number generation*/
  ranlib::Normal<double> normalGen(0,pow(2.0 * params.get<double>("KT",0.0) / params.get<double>("delta time",0.0),0.5));

  //fstoch consists of 4 values, two forces at each node
  LINALG::Matrix<6,1> aux;
  aux(0) = normalGen.random();
  aux(1) = normalGen.random();
  aux(2) = normalGen.random();
  aux(3) = normalGen.random();
  aux(4) = normalGen.random();
  aux(5) = normalGen.random();




  /*calculation of Brownian forces and damping is based on drag coefficient; this coefficient per unit
  * length is approximated by the one of an infinitely long staff for friciton orthogonal to staff axis*/
  double zeta = 4 * PI * lrefe_ * params.get<double>("ETA",0.0);

  int stochasticorder = params.get<int>("STOCH_ORDER",-1);

  switch(stochasticorder)
  {
  case -1:
  {
    //multiply S_loc (cholesky decomposition of C_loc) for C_loc only diagonal
    floc_(0,0) = pow(zeta/2,0.5)*aux(0,0);
    floc_(1,0) = pow(zeta/2,0.5)*aux(1,0);
    floc_(2,0) = pow(zeta/2,0.5)*aux(2,0);
    floc_(3,0) = pow(zeta/2,0.5)*aux(3,0);
    floc_(4,0) = pow(zeta/2,0.5)*aux(4,0);
    floc_(5,0) = pow(zeta/2,0.5)*aux(5,0);


  }
  break;
  case 0:
  {
    //multiply S_loc(cholesky decomposition of C_loc) for gamma_parallel=gamma_perp
  	floc_(0,0) = pow(zeta/3,0.5)*aux(0);
  	floc_(1,0) = pow(zeta/3,0.5)*aux(1);
  	floc_(2,0) = pow(zeta/3,0.5)*aux(2);
  	floc_(3,0) = pow(zeta/12,0.5)*aux(0)+pow(zeta/4,0.5)*aux(3);
  	floc_(4,0) = pow(zeta/12,0.5)*aux(1)+pow(zeta/4,0.5)*aux(4);
  	floc_(5,0) = pow(zeta/12,0.5)*aux(2)+pow(zeta/4,0.5)*aux(5);


  }
  break;
  case 1:
  {
    //multiply S_loc(cholesky decomposition of C_loc) for gamma_parallel=gamma_perp/2
  	floc_(0,0) = pow(zeta/6,0.5)*aux(0);
  	floc_(1,0) = pow(zeta/3,0.5)*aux(1);
  	floc_(2,0) = pow(zeta/3,0.5)*aux(2);
  	floc_(3,0) = pow(zeta/24,0.5)*aux(0)+pow(zeta/8,0.5)*aux(3);
  	floc_(4,0) = pow(zeta/12,0.5)*aux(1)+pow(zeta/4,0.5)*aux(4);
  	floc_(5,0) = pow(zeta/12,0.5)*aux(2)+pow(zeta/4,0.5)*aux(5);
  }
  break;
  default:
      dserror("Unknown type of stochasticorder for Beam3 %d", stochasticorder);
  }

return 0;
}//DRT::ELEMENTS::Beam3::ComputeLocalBrownianForces
/*-----------------------------------------------------------------------------------------------------------*
 | Assemble statistical forces and damping matrix according to fluctuation dissipation theorem (public) 06/09|
 *----------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Beam3::CalcBrownian(ParameterList& params,
                              vector<double>&           vel,  //!< element velocity vector
                              vector<double>&           disp, //!< element displacement vector
                              Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                              Epetra_SerialDenseVector* force) //!< element internal force vector
{
 int stochasticorder = params.get<int>("STOCH_ORDER",-2);// polynomial order for interpolation of stochastic line load

 //if no stochasticorder has been chosen explicitly we exit this function without any action
 if(stochasticorder == -2)
   return;

 double dt = params.get<double>("delta time",0.0);//timestep

 double zeta = 4 * PI * lrefe_ * params.get<double>("ETA",0.0);//damping parameter zeta

 /*the following block adds a rotational damping and related internal forces; this damping goes along with the damping of a
 * rigid straight rod spinning around its own axis; analytically it is given by a damping coefficient
 * gamma_a = 4*pi*eta*r^2*lrefe_ where r is the radius of the crosssection; as this
 * coefficient is very small for thin rods it is increased artificially by a factor for numerical convencience*/
  double rsquare = pow((4*Iyy_/PI),0.5);
  //gamma_a artificially increased by factor artificial; this factor should scale with lrefe_^2 and comprise a heuristically calculated constant (here 60)
  double artificial = 60;
  double gammaa = 4*PI*params.get<double>("ETA",0.0)*(rsquare)*artificial;


  //computing angle increment from current position in comparison with last converged position for damping
  LINALG::Matrix<4,1> deltaQ;
  LINALG::Matrix<3,1> deltatheta;
  quaternionproduct(inversequaternion(Qconv_[0]),Qnew_[0],deltaQ);
  quaterniontoangle(deltaQ,deltatheta);

  //angular velocity
  LINALG::Matrix<3,1> omega = deltatheta;
  omega.Scale(1/ params.get<double>("delta time",0.01));

  //computing special matrix for anisotropic damping
  LINALG::Matrix<3,3> Tconv;
  quaterniontotriad(Qconv_[0],Tconv);
  LINALG::Matrix<3,3> Theta;
  for(int i = 0; i<3; i++)
    for(int j = 0; j<3; j++)
      Theta(i,j) = Tconv(i,0)*Tconv(j,0);

  //inverse exponential map
  LINALG::Matrix<3,3> Hinverse = Hinv(deltatheta);

  //stiffness due to rotational damping
  LINALG::Matrix<3,3> artstiff;
  artstiff.Multiply(Theta,Hinverse);
  artstiff.Scale(gammaa*0.5 / params.get<double>("delta time",0.01));

  //forces due to rotational damping
  LINALG::Matrix<3,1> artforce;
  artforce.Multiply(Theta,omega);
  artforce.Scale(gammaa);

  //adding artificial contributions to stiffness matrix and force vector
  for(int i= 0; i<3; i++)
  {
    for(int j=0;j<3;j++)
    {
      (*stiffmatrix)(3+i, 3+j) += artstiff(i,j);
      (*stiffmatrix)(9+i, 9+j) += artstiff(i,j);
      (*stiffmatrix)(9+i, 3+j) += artstiff(i,j);
      (*stiffmatrix)(3+i, 9+j) += artstiff(i,j);
    }
  }

  for(int i = 0; i<3; i++)
  {
    (*force)[3+i] += artforce(i);
    (*force)[9+i] += artforce(i);
  }


//the following code blocks incalculate damping due to translation of the nodes

 switch(stochasticorder)
 {
 //simple isotropic model of Brownian motion with uncorrelated nodal forces
 case -1:
 {
 	
	//calc brownian damp matrix
	if (stiffmatrix != NULL) // necessary to run stiffmatrix control routine
	{
		(*stiffmatrix)(0,0) += zeta/(2.0*dt);
		(*stiffmatrix)(1,1) += zeta/(2.0*dt);
		(*stiffmatrix)(2,2) += zeta/(2.0*dt);
		(*stiffmatrix)(6,6) += zeta/(2.0*dt);
		(*stiffmatrix)(7,7) += zeta/(2.0*dt);
		(*stiffmatrix)(8,8) += zeta/(2.0*dt);
	}
	
	if (force != NULL)
	{
	 //calc int brownian forces
   (*force)(0) +=zeta/(2.0)*vel[0];
   (*force)(1) +=zeta/(2.0)*vel[1];
   (*force)(2) +=zeta/(2.0)*vel[2];
   (*force)(6) +=zeta/(2.0)*vel[6];
   (*force)(7) +=zeta/(2.0)*vel[7];
   (*force)(8) +=zeta/(2.0)*vel[8];


   //calc ext brownian forces
   (*force)(0) -=floc_(0,0);
   (*force)(1) -=floc_(1,0);
   (*force)(2) -=floc_(2,0);
   (*force)(6) -=floc_(3,0);
   (*force)(7) -=floc_(4,0);
   (*force)(8) -=floc_(5,0);	
	}

 }
 break;

 //isotropic model of Brownian motion with correlated forces
 case 0:
 {

   //calc brownian damp matrix
	 if (stiffmatrix != NULL) // necessary to run stiffmatrix control routine
	 {
   (*stiffmatrix)(0,0) += zeta/(3.0*dt);
   (*stiffmatrix)(1,1) += zeta/(3.0*dt);
   (*stiffmatrix)(2,2) += zeta/(3.0*dt);
   (*stiffmatrix)(6,6) += zeta/(3.0*dt);
   (*stiffmatrix)(7,7) += zeta/(3.0*dt);
   (*stiffmatrix)(8,8) += zeta/(3.0*dt);
   (*stiffmatrix)(0,6) += zeta/(6.0*dt);
   (*stiffmatrix)(1,7) += zeta/(6.0*dt);
   (*stiffmatrix)(2,8) += zeta/(6.0*dt);
   (*stiffmatrix)(6,0) += zeta/(6.0*dt);
   (*stiffmatrix)(7,1) += zeta/(6.0*dt);
   (*stiffmatrix)(8,2) += zeta/(6.0*dt);
	 }

	 if (force != NULL)
	 {
	 //calc int brownian forces
	 (*force)(0) +=zeta/(3.0)*vel[0]+zeta/(6.0)*vel[6];
   (*force)(1) +=zeta/(3.0)*vel[1]+zeta/(6.0)*vel[7];
   (*force)(2) +=zeta/(3.0)*vel[2]+zeta/(6.0)*vel[8];
   (*force)(6) +=zeta/(6.0)*vel[0]+zeta/(3.0)*vel[6];
   (*force)(7) +=zeta/(6.0)*vel[1]+zeta/(3.0)*vel[7];
   (*force)(8) +=zeta/(6.0)*vel[2]+zeta/(3.0)*vel[8];

   //calc ext brownian forces
   (*force)(0) -=floc_(0,0);
   (*force)(1) -=floc_(1,0);
   (*force)(2) -=floc_(2,0);
   (*force)(6) -=floc_(3,0);
   (*force)(7) -=floc_(4,0);
   (*force)(8) -=floc_(5,0);
	 }


  }
  break;
  //anisotropic model of Brownian motion with correlated nodal forces
  case 1:
  {
  	 //local damping matrix
  	 LINALG::Matrix<3,3> dampbasis(true);
  	 dampbasis(0,0)=zeta/2.0;
  	 dampbasis(1,1)=zeta;
  	 dampbasis(2,2)=zeta;
  	
  	 LINALG::Matrix<3,3> Tnew;//Rotation matrix in current iteration step
  	 quaterniontotriad(Qnew_[0],Tnew);
  	
  	 //rotate C to global
  	 LINALG::Matrix<3,3> aux1(true);
  	 aux1.Multiply(Tnew,dampbasis);
  	 dampbasis.MultiplyNT(aux1,Tnew);
  	
  	 //calc internal stiffness
  	 if (stiffmatrix !=NULL)
  	 {
			 //add stiffness due to variation of transl dofs to stiffmatrix
			 for(int i=0; i<3; i++)
			 {
				 for (int j=0; j<3; j++)
				 {
					(*stiffmatrix)(i,j) += dampbasis(i,j)/(3.0*dt);
					(*stiffmatrix)(i,j+6) += dampbasis(i,j)/(6.0*dt);
					(*stiffmatrix)(i+6,j) += dampbasis(i,j)/(6.0*dt);
					(*stiffmatrix)(i+6,j+6) += dampbasis(i,j)/(3.0*dt);
				 }
			 }
  	
			 /* für rotatorische freiheitsgrade noch multiplizieren (hering)
		//to calc the derivation of rot dof, in the 3D case there have to be multiplied Hinv(rotdof)
		//because of the non additivity of rot increments
			 LINALG::Matrix<3,1> theta1(true);//node1
			 LINALG::Matrix<3,1> theta2(true);//node2
			 for (int i=0; i<3; i++)
			 {
				 theta1(i,0)=disp[3+i];
				 theta2(i,0)=disp[9+i];
			 }
			
			 LINALG::Matrix<3,3>H1;
			 LINALG::Matrix<3,3>H2;
			 H1=Hinv(theta1);
			 H2=Hinv(theta2);

			//end Hinv
			 */
			
		
			
			
  	 //define seperate spin matrix for each rot dof
			
		 LINALG::Matrix<3,3> S1(true);//delta alpha1
		 S1(1,2)=-1.0;
		 S1(2,1)=1.0;
		
		 LINALG::Matrix<3,3> S2(true);//delta alpha2
		 S2(0,2)=1.0;
		 S2(2,0)=-1.0;
		
		 LINALG::Matrix<3,3> S3(true);//delta alpha3
		 S3(0,1)=-1.0;
		 S3(1,0)=1.0;
	
		 LINALG::Matrix<3,3> aux2(true);
		
		 //values due to alpha 1
		 aux1.Multiply(S1,dampbasis);//multiply S1 from left
		 aux2.Multiply(dampbasis,S1);//multiply S1 from right
		
		
		 for(int i=0; i<3; i++)//calc difference
   	 {
   		 for (int j=0; j<3; j++)
   		 {
   			aux1(i,j)=aux1(i,j)-aux2(i,j);
   		 }
   	 }
	
		 LINALG::Matrix<6,6> aux3(true); //to store result for two nodes multiply velcoef
		
		 for(int i=0; i<3; i++)//now two nodes
   	 {
   		 for (int j=0; j<3; j++)
   		 {
   			 aux3(i,j)=aux1(i,j)/3.0;
   			 aux3(i+3,j)=aux1(i,j)/6.0;
   			 aux3(i,j+3)=aux1(i,j)/6.0;
   			 aux3(i+3,j+3)=aux1(i,j)/3.0;
   		 }
   	 }
	
		 LINALG::Matrix<6,1> aux4(true);
		
		 for(int i=0; i<6; i++)//multiply velocity
   	 {
   		 for (int j=0; j<3; j++)
   		 {
   			 aux4(i,0)+=aux3(i,j)*vel[j];
   			 aux4(i,0)+=aux3(i,j+3)*vel[j+6];
   		 }
   	 }
	
		 for (int i=0; i<3; i++)//fill stiffness values due to delta alpha1 in stiffmatrix
		 {
			 (*stiffmatrix)(i,3)+=aux4(i,0)*0.5;
			 (*stiffmatrix)(i+6,3)+=aux4(i+3,0)*0.5;
			 (*stiffmatrix)(i,9)+=aux4(i,0)*0.5;
			 (*stiffmatrix)(i+6,9)+=aux4(i+3,0)*0.5;
		 } 		 //end delta alpha1
	
		 //values due to alpha 2
		 aux1.Multiply(S2,dampbasis);//multiply S2 from left
		 aux2.Multiply(dampbasis,S2);//multiply S2 from right
		
		
		 for(int i=0; i<3; i++)//calc difference
	   		 for (int j=0; j<3; j++)
	   			aux1(i,j)=aux1(i,j)-aux2(i,j);

 	
		
		 for(int i=0; i<3; i++)//now two nodes
   	 {
   		 for (int j=0; j<3; j++)
   		 {
   			 aux3(i,j)=aux1(i,j)/3.0;
   			 aux3(i+3,j)=aux1(i,j)/6.0;
   			 aux3(i,j+3)=aux1(i,j)/6.0;
   			 aux3(i+3,j+3)=aux1(i,j)/3.0;
   		 }
   	 }
 	
		
		 for (int i=0; i<6; i++)//set aux4 to zero
			 aux4(i,0)=0.0;
		
		 for(int i=0; i<6; i++)//multiply velocity
   	 {
   		 for (int j=0; j<3; j++)
   		 {
   			 aux4(i,0)+=aux3(i,j)*vel[j];
   			 aux4(i,0)+=aux3(i,j+3)*vel[j+6];
   		 }
   	 }
 	
		for (int i=0; i<3; i++)//fill stiffness values due to delta alpha2 in stiffmatrix
 		{
 			(*stiffmatrix)(i,4)+=aux4(i,0)*0.5;
 			(*stiffmatrix)(i+6,4)+=aux4(i+3,0)*0.5;
 			(*stiffmatrix)(i,10)+=aux4(i,0)*0.5;
 			(*stiffmatrix)(i+6,10)+=aux4(i+3,0)*0.5;
 			
 		}//end delta alpha2
		
		
 		//values due to alpha 3
		aux1.Multiply(S3,dampbasis);//multiply S3 from left
		aux2.Multiply(dampbasis,S3);//multiply S3 from right
		
		
		 for(int i=0; i<3; i++)//calc difference
	   	 {
	   		 for (int j=0; j<3; j++)
	   		 {
	   			aux1(i,j)=aux1(i,j)-aux2(i,j);
	   		 }
	   	 }
 	
		
		 for(int i=0; i<3; i++)//now two nodes
	   	 {
	   		 for (int j=0; j<3; j++)
	   		 {
	   			 aux3(i,j)=aux1(i,j)/3.0;
	   			 aux3(i+3,j)=aux1(i,j)/6.0;
	   			 aux3(i,j+3)=aux1(i,j)/6.0;
	   			 aux3(i+3,j+3)=aux1(i,j)/3.0;
	   		 }
	   	 }
 	
		
		 for (int i=0; i<6; i++)//set aux4 to zero
		 				 aux4(i,0)=0.0;
		
		 for(int i=0; i<6; i++)//multiply velocity
	   	 {
	   		 for (int j=0; j<3; j++)
	   		 {
	   			 aux4(i,0)+=aux3(i,j)*vel[j];
	   			 aux4(i,0)+=aux3(i,j+3)*vel[j+6];
	   		 }
	   	 }
 	
		 for (int i=0; i<3; i++)//fill stiffness values due to delta alpha3 in stiffmatrix
 		 {
 			 (*stiffmatrix)(i,5)+=aux4(i,0)*0.5;
 			 (*stiffmatrix)(i+6,5)+=aux4(i+3,0)*0.5;
 			 (*stiffmatrix)(i,11)+=aux4(i,0)*0.5;
 			 (*stiffmatrix)(i+6,11)+=aux4(i+3,0)*0.5;
 		 } //end delta alpha3
		
		//end internal stiffness
		
		//calc ext stiffness

   		 //begin delta alpha1
 		 aux1.Multiply(S1,Tnew);
		
		
		 for (int i=0; i<3; i++)
		 {
			 for (int j=0; j<3; j++)
			 {
				 aux4(i,0) = aux1(i,j)*floc_(j,0);
				 aux4(i+3,0) = aux1(i,j)*floc_(j+3,0);
			 }
		 }
		
		 for (int i=0; i<3; i++)//fill ext stiffness values due to delta alpha1 in stiffmatrix
 		 {
 			 (*stiffmatrix)(i,3)-=aux4(i,0)*0.5;
 			 (*stiffmatrix)(i+6,3)-=aux4(i+3,0)*0.5;
 			 (*stiffmatrix)(i,9)-=aux4(i,0)*0.5;
 			 (*stiffmatrix)(i+6,9)-=aux4(i+3,0)*0.5;
 		 }//end delta alpha 1
 		
 		
 		 //begin delta alpha2
 		 aux1.Multiply(S2,Tnew);
 		
 		 for (int i=0; i<6; i++)//set aux4 to zero
 			 aux4(i,0)=0.0;
 		  		
 		  		
		 for (int i=0; i<3; i++)
		 {
			 for (int j=0; j<3; j++)
			 {
				 aux4(i,0) += aux1(i,j)*floc_(j,0);
				 aux4(i+3,0) += aux1(i,j)*floc_(j+3,0);
			 }
		 }
		
		 for (int i=0; i<3; i++)//fill extstiffness values due to delta alpha2 in stiffmatrix
		 {
			 (*stiffmatrix)(i,4)-=aux4(i,0)*0.5;
			 (*stiffmatrix)(i+6,4)-=aux4(i+3,0)*0.5;
			 (*stiffmatrix)(i,10)-=aux4(i,0)*0.5;
			 (*stiffmatrix)(i+6,10)-=aux4(i+3,0)*0.5;
		 }//end delta alpha2
		
		 //begin delta alpha3
		 aux1.Multiply(S3,Tnew);
  		  		
		 for (int i=0; i<6; i++)//set aux4 to zero
  			 aux4(i,0)=0.0;
  		
 			 for (int i=0; i<3; i++)
 			 {
 				 for (int j=0; j<3; j++)
 				 {
 					 aux4(i,0) += aux1(i,j)*floc_(j,0);
 					 aux4(i+3,0) += aux1(i,j)*floc_(j+3,0);
 				 }
 			 }

 			 for (int i=0; i<3; i++)//fill ext stiffness values due to delta alpha3 in stiffmatrix
 		 {
 			 (*stiffmatrix)(i,5)-=aux4(i,0)*0.5;
 			 (*stiffmatrix)(i+6,5)-=aux4(i+3,0)*0.5;
 			 (*stiffmatrix)(i,11)-=aux4(i,0)*0.5;
 			 (*stiffmatrix)(i+6,11)-=aux4(i+3,0)*0.5;
 		 } //end delta alpha3
		
		//end calc ext stiffness
	
	 }//end stiffmatrix calc
	
	

	
	 if (force != NULL)
	 {
		
		 //calc internal forces
		 LINALG::Matrix<6,6> aux5(true);
		 for (int i=0; i<3; i++)
		 {
			 for (int j=0; j<3; j++)
			 {
				 aux5(i,j)=dampbasis(i,j)/3.0;
				 aux5(i+3,j)=dampbasis(i,j)/6.0;
				 aux5(i,j+3)=dampbasis(i,j)/6.0;
				 aux5(i+3,j+3)=dampbasis(i,j)/3.0;
			 }
		 }
		

		 for (int i=0; i<3; i++)
		 {
			 for (int j=0; j<3; j++)
			 {
				 (*force)(i) +=aux5(i,j)*vel[j];
				 (*force)(i) +=aux5(i,j+3)*vel[j+6];
				 (*force)(i+6) +=aux5(i+3,j)*vel[j];
				 (*force)(i+6) +=aux5(i+3,j+3)*vel[j+6];
			 }			
		 }
		
		 //calc external forces
		 for (int i=0; i<3; i++)
		 {
			 for (int j=0; j<3; j++)
			 {
				 (*force)(i) -= Tnew(i,j)*floc_(j,0);
				 (*force)(i+6) -= Tnew(i,j)*floc_(j+3,0);
			 }
		 } 			
		
	 }//end calc forces

  	
   } //end case 1
   break;
  }

return;

}//DRT::ELEMENTS::Beam2::CalcBrownian

/*-----------------------------------------------------------------------------------------------------------*
 | Evaluate PTC damping (public)                                                                  cyron 10/08|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3::EvaluatePTC(ParameterList& params,
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
   * artificial PTC damping has to be employed. 
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
   * Cross section A, moment of inertia I: from the above discussed physics one might assume a dependence of the
   * PTC parameter on the ratio of bending and strechting stiffness, i.e. on EI / EA. Such a dependence might
   * considerably exacerbate the application of the PTC algorithm. However, by means of numerical experiments a
   * different rule to deterime the PTC parameter was found: Beyond some ratio EI / EA simulations were found to
   * be unstable
   * 
   * 
   * Length L: reduing
   */



  double basisdamp   = (20e-2)*PI; //2 for network0_5.dat
  double anisofactor = 10;
  


  //computing angle increment from current position in comparison with last converged position for damping
  LINALG::Matrix<4,1> deltaQ;
  LINALG::Matrix<3,1> deltatheta;
  quaternionproduct(inversequaternion(Qconv_[0]),Qnew_[0],deltaQ);
  quaterniontoangle(deltaQ,deltatheta);

  //computing special matrix for anisotropic damping
  LINALG::Matrix<3,3> Tconv;
  quaterniontotriad(Qconv_[0],Tconv);
  LINALG::Matrix<3,3> Theta;
  for(int i = 0; i<3; i++)
    for(int j = 0; j<3; j++)
      Theta(i,j) = Tconv(i,0)*Tconv(j,0);

  //inverse exponential map
  LINALG::Matrix<3,3> Hinverse = Hinv(deltatheta);

  //isotropic artificial stiffness
  LINALG::Matrix<3,3> artstiff = Hinverse;
  artstiff.Scale(basisdamp);

  //anisotropic artificial stiffness
  LINALG::Matrix<3,3> auxstiff;
  auxstiff.Multiply(Theta,Hinverse);
  auxstiff.Scale(anisofactor*basisdamp);
  artstiff += auxstiff;



  //scale artificial damping with dti parameter for PTC method
  artstiff.Scale( params.get<double>("dti",0.0) );


  for(int i= 0; i<3; i++)
  {
    for(int j=0;j<3;j++)
    {
      /*
      //translational damping; seemingly not necessary, rather bad for convergence
      elemat1(  i,   j) += basisdamp;
      elemat1(6+i, 6+j) += basisdamp;
      elemat1(6+i,   j) += basisdamp;
      elemat1(  i, 6+j) += basisdamp;
      */

      //rotational damping
      elemat1(3+i, 3+j) += artstiff(i,j);
      elemat1(9+i, 9+j) += artstiff(i,j);
      elemat1(9+i, 3+j) += artstiff(i,j);
      elemat1(3+i, 9+j) += artstiff(i,j);
    }
  }

  return 0;
} //DRT::ELEMENTS::Beam3::EvaluatePTC




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
      {
        StTCmTt(i,j) += S(k,i)*TCmTt(k,j);
      }
    }
  }

  //calculating StTCmTtS
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      StTCmTtS(i,j) = 0.0;
      for (int k = 0; k < 3; ++k)
      {
        StTCmTtS(i,j) += StTCmTt(i,k)*S(k,j);
      }
    }
  }

  //calculating TCmTtS
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      TCmTtS(i,j) = 0.0;
      for (int k = 0; k < 3; ++k)
      {
        TCmTtS(i,j) += TCmTt(i,k)*S(k,j);
      }
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
     *   LINALG::Matrix<4,1> q;
     *   angletoquaternion(deltatheta,q);
     *   quaterniontoangle(q,deltatheta);
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
  double ym;
  double sm;
  double density;

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

  //Get integrationpoints for the applied gaussrule
  DRT::UTILS::IntegrationPoints1D gausspoints(gaussrule_);

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
				dxdxi_gp(dof)              += (Nodes()[node]->X()[dof]+disp[6*node+dof])*deriv(node);
				
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
	  //The mass matrix has to be integrated completely. Therefore we use more gausspoints
	  gaussrule_ = static_cast<enum DRT::UTILS::GaussRule1D>(nnode);
	
	  //Get the applied integrationpoints for complete integration
	  DRT::UTILS::IntegrationPoints1D gausspointsmass(gaussrule_);	
	
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
	  	
	  //reset gaussrule
	  gaussrule_ = static_cast<enum DRT::UTILS::GaussRule1D>(nnode-1);

  }//if (massmatrix != NULL)

  /*the following function call applied statistical forces and damping matrix according to the fluctuation dissipation theorem;
   * it is dedicated to the application of beam2 elements in the frame of statistical mechanics problems; for these problems a
   * special vector has to be passed to the element packed in the params parameter list; in case that the control routine calling
   * the element does not attach this special vector to params the following method is just doing nothing, which means that for
   * any ordinary problem of structural mechanics it may be ignored*/
   CalcBrownian(params,vel,disp,stiffmatrix,force);


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

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_BEAM3
