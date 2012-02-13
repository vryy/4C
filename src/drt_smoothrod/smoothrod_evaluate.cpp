/*!-----------------------------------------------------------------------------------------------------------
 \file smoothrod_evaluate.cpp
 \brief

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

 *-----------------------------------------------------------------------------------------------------------*/
#ifdef D_SMOOTHROD
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "smoothrod.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../linalg/linalg_fixedsizematrix.H"



/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                                 cyron 01/08|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Smoothrod::Evaluate(ParameterList& params,
    DRT::Discretization& discretization,
    vector<int>& lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::Smoothrod::ActionType act = Smoothrod::calc_none;
  // get the action required
  string action = params.get<string>("action","calc_none");
  if (action == "calc_none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff") act = Smoothrod::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff") act = Smoothrod::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = Smoothrod::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass") act = Smoothrod::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass") act = Smoothrod::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass") act = Smoothrod::calc_struct_nlnstifflmass; //with lumped mass matrix
  else if (action=="calc_struct_stress") act = Smoothrod::calc_struct_stress;
  else if (action=="calc_struct_eleload") act = Smoothrod::calc_struct_eleload;
  else if (action=="calc_struct_fsiload") act = Smoothrod::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep") act = Smoothrod::calc_struct_update_istep;
  else if (action=="calc_struct_update_imrlike") act = Smoothrod::calc_struct_update_imrlike;
  else if (action=="calc_struct_reset_istep") act = Smoothrod::calc_struct_reset_istep;
  else if (action=="calc_struct_ptcstiff")        act = Smoothrod::calc_struct_ptcstiff;
  else dserror("Unknown type of action for Smoothrod");

  switch(act)
  {
    case Smoothrod::calc_struct_ptcstiff:
    {
      //only nonlinear case implemented!
      dserror("no ptc implemented for smoothrod element");
    }
    break;
    /*in case that only linear stiffness matrix is required sr_nlstiffmass is called with zero dispalcement and
     residual values*/
    case Smoothrod::calc_struct_linstiff:
    {
      //only nonlinear case implemented!
      dserror("linear stiffness matrix called, but not implemented");

    }
    break;

    //nonlinear stiffness and mass matrix are calculated even if only nonlinear stiffness matrix is required
    case Smoothrod::calc_struct_nlnstiffmass:
    case Smoothrod::calc_struct_nlnstifflmass:
    case Smoothrod::calc_struct_nlnstiff:
    case Smoothrod::calc_struct_internalforce:
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

      if (act == Smoothrod::calc_struct_nlnstiffmass)
      {
    	  switch(nnode)
    	  {
  	  		case 2:
  	  		{
  	  			sr_nlnstiffmass<2>(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
  	  			break;
  	  		}
  	  		case 3:
  	  		{
  	  			sr_nlnstiffmass<3>(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
  	  			break;
  	  		}
  	  		case 4:
  	  		{
  	  			sr_nlnstiffmass<4>(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
  	  			break;
  	  		}
  	  		case 5:
  	  		{
  	  			sr_nlnstiffmass<5>(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
  	  			break;
  	  		}
  	  		default:
  	  			dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
    	  }
      }
      else if (act == Smoothrod::calc_struct_nlnstifflmass)
      {
    	  switch(nnode)
    	  {
  	  		case 2:
  	  		{
  	  			sr_nlnstiffmass<2>(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
  	  			lumpmass<2>(&elemat2);
  	  			break;
  	  		}
  	  		case 3:
  	  		{
  	  			sr_nlnstiffmass<3>(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
  	  			lumpmass<3>(&elemat2);
  	  			break;
  	  		}
  	  		case 4:
  	  		{
  	  			sr_nlnstiffmass<4>(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
  	  			lumpmass<4>(&elemat2);
  	  			break;
  	  		}
  	  		case 5:
  	  		{
  	  			sr_nlnstiffmass<5>(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
  	  			lumpmass<5>(&elemat2);
  	  			break;
  	  		}
  	  		default:
  	  			dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
    	  }
      }
      else if (act == Smoothrod::calc_struct_nlnstiff)
      {
    	  switch(nnode)
    	  {
  	  		case 2:
  	  		{
  	  			sr_nlnstiffmass<2>(params,myvel,mydisp,&elemat1,NULL,&elevec1);
  	  			break;
  	  		}
  	  		case 3:
  	  		{
  	  			sr_nlnstiffmass<3>(params,myvel,mydisp,&elemat1,NULL,&elevec1);
  	  			break;
  	  		}
  	  		case 4:
  	  		{
  	  			sr_nlnstiffmass<4>(params,myvel,mydisp,&elemat1,NULL,&elevec1);
  	  			break;
  	  		}
  	  		case 5:
  	  		{
  	  			sr_nlnstiffmass<5>(params,myvel,mydisp,&elemat1,NULL,&elevec1);
  	  			break;
  	  		}
  	  		default:
  	  			dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
    	  }
      }

      else if (act == Smoothrod::calc_struct_internalforce)
      {
    	  switch(nnode)
    	  {
  	  		case 2:
  	  		{
  	  			sr_nlnstiffmass<2>(params,myvel,mydisp,NULL,NULL,&elevec1);
  	  			break;
  	  		}
  	  		case 3:
  	  		{
  	  			sr_nlnstiffmass<3>(params,myvel,mydisp,NULL,NULL,&elevec1);
  	  			break;
  	  		}
  	  		case 4:
  	  		{
  	  			sr_nlnstiffmass<4>(params,myvel,mydisp,NULL,NULL,&elevec1);
  	  			break;
  	  		}
  	  		case 5:
  	  		{
  	  			sr_nlnstiffmass<5>(params,myvel,mydisp,NULL,NULL,&elevec1);
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
      		
		  //sr_nlnstiffmass is a templated function. therefore we need to point out the number of nodes in advance
        	  switch(nnode)
        	  {
        	  		case 2:
        	  		{
        	  			sr_nlnstiffmass<2>(params,vel_aux,disp_aux,NULL,NULL,&force_aux);       	  			
        	  			break;
        	  		}
        	  		case 3:
        	  		{
        	  			sr_nlnstiffmass<3>(params,vel_aux,disp_aux,NULL,NULL,&force_aux);
        	  			break;
        	  		}
        	  		case 4:
        	  		{
        	  			sr_nlnstiffmass<4>(params,vel_aux,disp_aux,NULL,NULL,&force_aux);
        	  			break;
        	  		}
        	  		case 5:
        	  		{
        	  			sr_nlnstiffmass<5>(params,vel_aux,disp_aux,NULL,NULL,&force_aux);
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
      dserror("Unknown type of action for Smoothrod %d", act);
  }
  return 0;

}

/*-----------------------------------------------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)                                       cyron 03/08|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Smoothrod::EvaluateNeumann(ParameterList& params,
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



/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   cyron 03/10|
 *-----------------------------------------------------------------------------------------------------------*/
template<int nnode>
void DRT::ELEMENTS::Smoothrod::sr_nlnstiffmass( ParameterList& params,
                                            vector<double>&           vel,
                                            vector<double>&           disp,
                                            Epetra_SerialDenseMatrix* stiffmatrix,
                                            Epetra_SerialDenseMatrix* massmatrix,
                                            Epetra_SerialDenseVector* force)
{




  return;

} // DRT::ELEMENTS::Smoothrod::sr_nlnstiffmass

/*------------------------------------------------------------------------------------------------------------*
 | lump mass matrix					   (private)                                                           cyron 01/08|
 *------------------------------------------------------------------------------------------------------------*/
template<int nnode>
void DRT::ELEMENTS::Smoothrod::lumpmass(Epetra_SerialDenseMatrix* emass)
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
