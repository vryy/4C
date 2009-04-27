/*!-----------------------------------------------------------------------------------------------------------
\file beam2_evaluate.cpp
\brief

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*-----------------------------------------------------------------------------------------------------------*/
#ifdef D_BEAM2R
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "beam2r.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/stvenantkirchhoff.H"

/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                                 cyron 01/08|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam2r::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::Beam2r::ActionType act = Beam2r::calc_none;
  // get the action required
  string action = params.get<string>("action","calc_none");
  if (action == "calc_none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")      act = Beam2r::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")      act = Beam2r::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = Beam2r::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")  act = Beam2r::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")  act = Beam2r::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass") act = Beam2r::calc_struct_nlnstifflmass;
  else if (action=="calc_struct_stress")        act = Beam2r::calc_struct_stress;
  else if (action=="calc_struct_eleload")       act = Beam2r::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")       act = Beam2r::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")  act = Beam2r::calc_struct_update_istep;
  else if (action=="calc_struct_update_imrlike") act = Beam2r::calc_struct_update_imrlike;
  else if (action=="calc_struct_reset_istep")   act = Beam2r::calc_struct_reset_istep;
  else if (action=="calc_brownian_damp")       act = Beam2r::calc_brownian_damp;
  else if (action=="calc_struct_ptcstiff")        act = Beam2r::calc_struct_ptcstiff;
  else dserror("Unknown type of action for Beam2r");

  switch(act)
  {
    case Beam2r::calc_struct_ptcstiff:
    {
    	dserror("Beam2r element does'nt need any special ptc tools to allow stable implicit dynamics with acceptable time step size");
    }
    break;
    //action type for evaluating statistical forces
    case Beam2r::calc_brownian_damp:
    {
    	/*calculation of Brownian forces and damping is based on drag coefficient; this coefficient per unit
    	 * length is approximated by the one of an infinitely long staff for friciton orthogonal to staff axis*/
    	 double zeta = 4 * PI * lrefe_ * params.get<double>("ETA",0.0);
    	          
    	 // get element displacements (for use in shear flow fields)
    	 RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
    	 if (disp==null) dserror("Cannot get state vector 'displacement'");
    	 vector<double> mydisp(lm.size());
    	 DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
    	      
    	      
    	 //first we evaluate the damping matrix
    	 EvaluateBrownianDamp(params,mydisp,zeta,elemat1);
    	      
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
    	 if(this->Owner() != discretization.Comm().MyPID()) return 0;
    	     
    	 //evaluation of statistical forces or displacements
    	 Epetra_SerialDenseVector brownian(lm.size());
    	 EvaluateBrownianForces(params,mydisp,zeta,brownian);

    	 /*all the above evaluated forces or movements are assembled*/
    	 //note carefully: a space between the two subsequal ">" signs is mandatory for the C++ parser in order to avoid confusion with ">>" for streams
    	 RCP<Epetra_Vector>    browniancol = params.get<  RCP<Epetra_Vector> >("statistical vector",Teuchos::null);

    	 for(unsigned int i = 0; i < lm.size(); i++)
    	  {
    		 //note: lm contains the global Ids of the degrees of freedom of this element
    	     //testing whether the browniancol vector has really an element related with the i-th element of brownian by the i-the entry of lm
    	     if (!(browniancol->Map()).MyGID(lm[i])) dserror("Sparse vector browniancol does not have global row %d",lm[i]);

    	     //get local Id of the fstatcol vector related with a certain element of fstat
    	     int lid = (browniancol->Map()).LID(lm[i]);

    	     //add to the related element of fstatcol the contribution of fstat
    	     (*browniancol)[lid] += brownian[i];
    	  }
    }
    break;
    /*in case that only linear stiffness matrix is required b2_nlstiffmass is called with zero displacement and
     residual values*/
    case Beam2r::calc_struct_linstiff:
    {
      //only nonlinear case implemented!
      dserror("linear stiffness matrix called, but not implemented");
    }
    break;

    //nonlinear stiffness and mass matrix are calculated even if only nonlinear stiffness matrix is required
    case Beam2r::calc_struct_nlnstiffmass:
    case Beam2r::calc_struct_nlnstifflmass:
    case Beam2r::calc_struct_nlnstiff:
    case Beam2r::calc_struct_internalforce:
    {
      int lumpedmass = 0;  // 0=consistent, 1=lumped
      if (act==Beam2r::calc_struct_nlnstifflmass) lumpedmass = 1;

      // need current global displacement and residual forces and get them from discretization
      // making use of the local-to-global map lm one can extract current displacement and residual values for each degree of freedom
      //
      // get element displacements
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

      // determine element matrices and forces
      if (act == Beam2r::calc_struct_nlnstiffmass)
        nlnstiffmass(params,myvel,mydisp,&elemat1,&elemat2,&elevec1,lumpedmass);
      else if (act == Beam2r::calc_struct_nlnstifflmass)
      {
        nlnstiffmass(params,myvel,mydisp,&elemat1,&elemat2,&elevec1,lumpedmass);
      }
      else if (act == Beam2r::calc_struct_nlnstiff)
        nlnstiffmass(params,myvel,mydisp,&elemat1,NULL,&elevec1,lumpedmass);
      else if  (act ==  calc_struct_internalforce)
        nlnstiffmass(params,myvel,mydisp,NULL,NULL,&elevec1,lumpedmass);
      

      /*
      //the following code block can be used to check quickly whether the nonlinear stiffness matrix is calculated
      //correctly or not by means of a numerically approximated stiffness matrix
      //if(Id() == 3) //limiting the following tests to certain element numbers
      {
        //variable to store numerically approximated stiffness matrix
        Epetra_SerialDenseMatrix stiff_approx;
        stiff_approx.Shape(6,6);

        //relative error of numerically approximated stiffness matrix
        Epetra_SerialDenseMatrix stiff_relerr;
        stiff_relerr.Shape(6,6);

        //characteristic length for numerical approximation of stiffness
        double h_rel = 1e-8;

        //flag indicating whether approximation lead to significant relative error
        int outputflag = 0;

        //calculating strains in new configuration
        for(int i=0; i<3; i++)
        {
          for(int k=0; k<2; k++)
          {

            Epetra_SerialDenseVector force_aux;
            force_aux.Size(6);

            //create new displacement and velocity vectors in order to store artificially modified displacements
            vector<double> vel_aux(6);
            vector<double> disp_aux(6);
            for(int id = 0;id<6;id++)
            {
                DRT::UTILS::ExtractMyValues(*disp,disp_aux,lm);
                DRT::UTILS::ExtractMyValues(*vel,vel_aux,lm);
            }

            //modifying displacment artificially (for numerical derivative of internal forces):
            disp_aux[i + 3*k] = disp_aux[i + 3*k] + h_rel;
             vel_aux[i + 3*k] =  vel_aux[i + 3*k] + h_rel * params.get<double>("gamma",0.581) / ( params.get<double>("delta time",0.01)*params.get<double>("beta",0.292) );


            nlnstiffmass(params,vel_aux,disp_aux,NULL,NULL,&force_aux,lumpedmass);


            for(int u = 0;u<6;u++)
            {
              stiff_approx(u,i+k*3)= ( pow(force_aux[u],2) - pow(elevec1(u),2) )/ (h_rel * (force_aux[u] + elevec1(u) ) );
            }

          }
        }

       for(int line=0; line<6; line++)
       {
         for(int col=0; col<6; col++)
         {
           stiff_relerr(line,col)= fabs( ( pow(elemat1(line,col),2) - pow(stiff_approx(line,col),2) )/ ( (elemat1(line,col) + stiff_approx(line,col)) * elemat1(line,col) ));

           //suppressing small entries whose effect is only confusing and NaN entires (which arise due to zero entries)
           if ( fabs( stiff_relerr(line,col) ) < h_rel*500 || isnan( stiff_relerr(line,col)) || elemat1(line,col) == 0)
             stiff_relerr(line,col) = 0;

           if ( stiff_relerr(line,col) > 0)
             outputflag = 1;
         }
       }

       if(outputflag ==1)
       {
         std::cout<<"\n\n acutally calculated stiffness matrix"<< elemat1;
         std::cout<<"\n\n approximated stiffness matrix"<< stiff_approx;
         std::cout<<"\n\n rel error stiffness matrix"<< stiff_relerr;
       }
      }
      //end of section in which numerical approximation for stiffness matrix is computed
      */

    }
    break;
    case calc_struct_update_istep:
    case calc_struct_update_imrlike:
    case calc_struct_reset_istep:
    break;
    case calc_struct_stress:
      dserror("No stress output implemented for beam2r elements");
    default:
      dserror("Unknown type of action for beam2r %d", act);
  }
  return 0;

}


/*-----------------------------------------------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)                                       cyron 01/08|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Beam2r::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{
  // element displacements
  RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
  if (disp==null) dserror("Cannot get state vector 'displacement'");
  vector<double> mydisp(lm.size());
  DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);


  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const vector<int>* curve  = condition.Get<vector<int> >("curve");
  int curvenum = -1;
  // number of the load curve related with a specific line Neumann condition called
  if (curve) curvenum = (*curve)[0];
  // amplitude of load curve at current time called
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);

  //jacobian determinant
  double det = lrefe_/2;


  // no. of nodes on this element
  const int iel = NumNode();
  const int numdf = 3;
  const DiscretizationType distype = this->Shape();

  // gaussian points
  const DRT::UTILS::IntegrationPoints1D  intpoints(gaussrule_);


  //declaration of variable in order to store shape function
  Epetra_SerialDenseVector      funct(iel);

  // get values and switches from the condition

  // onoff is related to the first 6 flags of a line Neumann condition in the input file;
  // value 1 for flag i says that condition is active for i-th degree of freedom
  const vector<int>*    onoff = condition.Get<vector<int> >("onoff");
  // val is related to the 6 "val" fields after the onoff flags of the Neumann condition
  // in the input file; val gives the values of the force as a multiple of the prescribed load curve
  const vector<double>* val   = condition.Get<vector<double> >("val");

  //integration loops
  for (int ip=0; ip<intpoints.nquad; ++ip)
  {
    //integration points in parameter space and weights
    const double xi = intpoints.qxg[ip][0];
    const double wgt = intpoints.qwgt[ip];

    //evaluation of shape funcitons at Gauss points
    DRT::UTILS::shape_function_1D(funct,xi,distype);

    double fac=0;
    fac = wgt * det;

    // load vector ar
    double ar[numdf];
    // loop the dofs of a node

    for (int i=0; i<numdf; ++i)
    {
      ar[i] = fac * (*onoff)[i]*(*val)[i]*curvefac;
    }

    //sum up load components
    for (int node=0; node<iel; ++node)
      for (int dof=0; dof<numdf; ++dof)
         elevec1[node*numdf+dof] += funct[node] *ar[dof];

  } // for (int ip=0; ip<intpoints.nquad; ++ip)

  return 0;
}

/*-----------------------------------------------------------------------------------------------------------*
 | Evaluate damping matrix for Brownian motion (public)                                           cyron 03/09|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Beam2r::EvaluateBrownianDamp(ParameterList& params,
                                              vector<double> mydisp,
                                              double zeta,
                                              Epetra_SerialDenseMatrix& elemat1)
{
	// polynomial order for interpolation of stochastic line load (zero corresponds to bead spring model)
	  int stochasticorder = params.get<int>("STOCH_ORDER",0);
	  
	  //simple model with isotropic Brownian motion and uncorrelated nodal forces (like bead-spring-model)
	  if (stochasticorder == 0)
	  {
	    elemat1(0,0) += zeta/2.0;
	    elemat1(1,1) += zeta/2.0;
	    elemat1(3,3) += zeta/2.0;
	    elemat1(4,4) += zeta/2.0;  
	  }

	  //in case of a consistent damping matrix stochastic nodal forces are calculated consistently by methods of weighted integrals
	  else if (stochasticorder == 1)
	  {     
		  dserror("Consistent damp matrix called but not implemented. (stochasticorder=1)");
	  }
}
 //DRT::ELEMENTS::Beam2r::EvaluateBrownianDamp

/*-----------------------------------------------------------------------------------------------------------*
 | Evaluate forces for Brownian motion (public)                                                  cyron 03/09|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Beam2r::EvaluateBrownianForces(ParameterList& params,
                                                vector<double> mydisp,
                                                double zeta,
                                                Epetra_SerialDenseVector& elevec1)
{
	// thermal energy responsible for statistical forces
	  double kT = params.get<double>("KT",0.0);

	  // polynomial order for interpolation of stochastic line load (zero corresponds to bead spring model)
	  int stochasticorder = params.get<int>("STOCH_ORDER",0);
	  
	  //simple isotropic model of Brownian motion with uncorrelated nodal forces
	  if (stochasticorder == 0)
	  {
	    //calculating standard deviation of statistical forces according to fluctuation dissipation theorem
	    double stand_dev_trans = pow(2.0 * kT * (zeta/2.0) / params.get<double>("delta time",0.01),0.5);

	    //creating a random generator object which creates random numbers with mean = 0 and standard deviation
	    //stand_dev; using Blitz namespace "ranlib" for random number generation
	    ranlib::Normal<double> normalGen(0,stand_dev_trans);

	    //adding statistical forces
	    elevec1[0] += normalGen.random();
	    elevec1[1] += normalGen.random();
	    elevec1[3] += normalGen.random();
	    elevec1[4] += normalGen.random(); 
	  }

	  //in case of a consistent damping matrix stochastic nodal forces are calculated consistently by methods of weighted integrals
	  else if (stochasticorder == 1)
	  {     
		  dserror("Consistent damp matrix called but not implemented. (stochasticorder=1)");
	  }

	  return 0;
} //DRT::ELEMENTS::Beam2r::EvaluateBrownianForces


/*-----------------------------------------------------------------------------------------------------------*
 | evaluate auxiliary vectors and matrices for Reissner`s formulation (private)                   cyron 01/08|
 *----------------------------------------------------------------------------------------------------------*/
//notation for this function similar to Crisfield, Volume 1, section 7.4;
inline void DRT::ELEMENTS::Beam2r::local_aux(LINALG::Matrix<3,6>& Bcurr,
                                        LINALG::Matrix<6,1>& rcurr,
                                        LINALG::Matrix<6,1>& zcurr,
                                        LINALG::Matrix<6,1>& s,
                                        const double& lcurr,
                                        const double& thetaav,
                                        const double& c1,
                                        const double& c2)

{
  double cos_thetaav = cos(thetaav);
  double sin_thetaav = sin(thetaav);
  
  //vector r according to Crisfield, Vol. 1, (7.62)
  rcurr(0) = -cos_thetaav;
  rcurr(1) = -sin_thetaav;
  rcurr(2) = 0;
  rcurr(3) = cos_thetaav;
  rcurr(4) = sin_thetaav;
  rcurr(5) = 0;

  //vector z according to Crisfield, Vol. 1, (7.66)
  zcurr(0) = sin_thetaav;
  zcurr(1) = -cos_thetaav;
  zcurr(2) = 0;
  zcurr(3) = -sin_thetaav;
  zcurr(4) = cos_thetaav;
  zcurr(5) = 0;
  
  //assigning values to each element of the Bcurr matrix, Crisfield, Vol. 1, (7.118)
  for(int id_col=0; id_col<6; id_col++)
    {
      Bcurr(0,id_col) = rcurr(id_col) + c2 * s(id_col);
      Bcurr(1,id_col) = 0;
      Bcurr(2,id_col) = zcurr(id_col) + c1 * s(id_col);
    }
    Bcurr(1,2) = -1;
    Bcurr(1,5) =  1;

  return;
} /* DRT::ELEMENTS::Beam2r::local_aux */

/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   cyron 01/08|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam2r::nlnstiffmass( ParameterList& params,
                                            vector<double>&           vel,
                                            vector<double>&           disp,
                                            Epetra_SerialDenseMatrix* stiffmatrix,
                                            Epetra_SerialDenseMatrix* massmatrix,
                                            Epetra_SerialDenseVector* force,
                                            int lumpedmass)
{
  const int numdf = 3;
  const int iel = NumNode();
  //coordinates in current configuration of all the nodes in two dimensions stored in numdf x iel matrices
  LINALG::Matrix<3,2> xcurr;

  //current length of beam in physical space
  double lcurr = 0.0;

  //some geometric auxiliary variables according to Crisfield, Vol. 1 section 7.4
  LINALG::Matrix<6,1> zcurr;
  LINALG::Matrix<6,1> rcurr;
  LINALG::Matrix<6,1> s;
  LINALG::Matrix<3,6> Bcurr;
  //auxiliary matrix storing the product of constitutive matrix C and Bcurr
  LINALG::Matrix<3,6> aux_CB;
  //declaration of local internal forces
  LINALG::Matrix<3,1> force_loc;
  //declaration of material parameters
  double ym; //Young's modulus
  double sm; //shear modulus
  double density; //density

  //Inserting current configuration into xcurr
  for (int k=0; k<iel; ++k)
  {
    xcurr(0,k) = Nodes()[k]->X()[0] + disp[k*numdf+0]; //x for each node in figure 7.9 
    xcurr(1,k) = Nodes()[k]->X()[1] + disp[k*numdf+1]; //z for each node  in figure 7.9

    /*note that xcurr(2,0),xcurr(2,1) are local angles in Crisfield, Vol. 1. They refer to the
     * reference configuration. The exact global angle is never used in the further description
     * ,we will only need the difference. Therefore we do not store it. But we will be able to compute
     * thetaav from thetaav0 and the inkrements xcurr(2,0),xcurr(2,1)  */
    
    xcurr(2,k) = disp[k*numdf+2]; //theta for each node  in figure 7.9
  }

  //current length
  lcurr = pow( pow(xcurr(0,1)-xcurr(0,0),2) + pow(xcurr(1,1)-xcurr(1,0),2) , 0.5 );
  
  //calculate average absolute angle thetaav in current configuration
  double thetaav = thetaav0_ + 0.5 * (xcurr(2,0) + xcurr(2,1));
  
  //calculate sinus and cosinus beta for local_aux from (7.116) and (7.117)
  double sin_beta = ((xcurr(1,1)-xcurr(1,0)) * cos(thetaav) - (xcurr(0,1) - xcurr(0,0)) * sin(thetaav)) / lcurr;
  double cos_beta = ((xcurr(0,1)-xcurr(0,0)) * cos(thetaav) + (xcurr(1,1) - xcurr(1,0)) * sin(thetaav)) / lcurr;
  
  //vector s according to Crisfield, Vol. 1, (7.119)
    s(0) = 0;
    s(1) = 0;
    s(2) = 1;
    s(3) = 0;
    s(4) = 0;
    s(5) = 1;
    
  //calculation of factors needed for Bcurr; (7.120)  and (7.121)
      double c1 = -0.5 * lcurr;
      double c2 = 0.5 * lcurr * sin_beta;

  //calculation of local geometrically important matrices and vectors
  local_aux(Bcurr,rcurr,zcurr,s,lcurr,thetaav,c1,c2);


  // get the material law
  Teuchos::RCP<const MAT::Material> currmat = Material();
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

  //Crisfield, Vol. 1, (7.55) and (7.116)
  force_loc(0) = ym * crosssec_ * ( lcurr * cos_beta - lrefe_ ) / lrefe_;

  //local internal bending moment, Crisfield, Vol. 1, (7.108). Note that for the difference we won`t need the exact thetas
  //we can compute from from our delta_thetas
  force_loc(1) = ym * mominer_ * (xcurr(2,1)-xcurr(2,0)) / lrefe_;

  //local internal shear force, Crisfield, Vol. 1, (7.98) and gamma from (7.117)
  /*gamma is  here given in (7.114)*/
  force_loc(2) = sm * crosssecshear_ * ( lcurr * sin_beta / lrefe_ );

  if (force != NULL)
  {
  //declaration of global internal force
  LINALG::Matrix<6,1> force_glob;
  
  //calculation of global internal forces from Crisfield, Vol. 1, (7.124): q_i = B^T q_{li}
  force_glob.MultiplyTN(Bcurr,force_loc); //führt die o.g. Operation aus. Funktion der class LINALG::Matrix

    for(int k = 0; k<6; k++)
      (*force)(k) = force_glob(k); // force Vektor ist hier fertig berechnet und wurde global übergeben
  }


  //calculating tangential stiffness matrix in global coordinates, Crisfield, Vol. 1, (7.107)
  if (stiffmatrix != NULL)
  {
    //declaration of fixed size matrix for global stiffness
    LINALG::Matrix<6,6> stiff_glob;
    
    //linear elastic part of tangential stiffness matrix including rotation: K_t1 = B^T C_l B / l_0 (7.128)
    for(int id_col=0; id_col<6; id_col++)
    {
      aux_CB(0,id_col) = Bcurr(0,id_col) * (ym*crosssec_/lrefe_);// aux_CB oben definiert. Nur als Ablage für den Zwischenwert da (3x6)
      aux_CB(1,id_col) = Bcurr(1,id_col) * (ym*mominer_/lrefe_);
      aux_CB(2,id_col) = Bcurr(2,id_col) * (sm*crosssecshear_/lrefe_);
    }//im Endeffekt C^T B
    
    stiff_glob.MultiplyTN(aux_CB,Bcurr);
    
    //adding geometric stiffness by axial force: N (s z^T + z s^T) / 2 + N * c1 s s^T/2 (7.128)
    double aux_N_fac = force_loc(0) / 2;
    for(int id_lin=0; id_lin<6; id_lin++)
        for(int id_col=0; id_col<6; id_col++)
        {
          stiff_glob(id_lin,id_col) += aux_N_fac * ( s(id_lin) * zcurr(id_col) + zcurr(id_lin) * s(id_col));
          stiff_glob(id_lin,id_col) += aux_N_fac * c1 * s(id_lin) * s(id_col);
        }

    //adding geometric stiffness by shear force: -Q ( s r^T + r s^T ) / 2 - Q * c2 * s s^T / 2 (7.128)
    double aux_Q_fac = force_loc(2)/2;
    for(int id_lin=0; id_lin<6; id_lin++)
        for(int id_col=0; id_col<6; id_col++)
        {
          stiff_glob(id_lin,id_col) -= aux_Q_fac * ( s(id_lin) * rcurr(id_col) + rcurr(id_lin) * s(id_col));
          stiff_glob(id_lin,id_col) -= aux_Q_fac * c2 * s(id_lin) * s(id_col);
        }
    //shifting values from fixed size matrix to epetra matrix *stiffmatrix
    for(int i = 0; i < 6; i++)
      for(int j = 0; j < 6; j++)
        (*stiffmatrix)(i,j) = stiff_glob(i,j);// Steifigkeitsmatrix ist hier vollständig berechnet und wird global übergeben
    
  }
  


  //calculating mass matrix (local version = global version)
  if (massmatrix != NULL)
  {
      //if lumpedmass == 0 a consistent mass Timoshenko beam mass matrix is applied
      if (lumpedmass == 0)
      {
        //assignment of massmatrix by means of auxiliary diagonal matrix aux_E stored as an array
        double aux_E[3]={density*lrefe_*crosssec_/6.0, density*lrefe_*crosssec_/6.0, density*lrefe_*mominer_/6.0};
        for(int id=0; id<3; id++)
        {
              (*massmatrix)(id,id)     = 2.0*aux_E[id];
              (*massmatrix)(id+3,id+3) = 2.0*aux_E[id];
              (*massmatrix)(id,id+3)   = aux_E[id];
              (*massmatrix)(id+3,id)   = aux_E[id];
        }
      }
      
      /*if lumpedmass == 1 a lumped mass matrix is applied where the cross sectional moment of inertia is
       * assumed to be approximately zero so that the 3,3 and 5,5 element are both zero */
      else if (lumpedmass == 1)
      {
        //note: this is not an exact lumped mass matrix, but it is modified in such a way that it leads
        //to a diagonal mass matrix with constant diagonal entries
        (*massmatrix)(0,0) = density*lrefe_*crosssec_/2.0;
        (*massmatrix)(1,1) = density*lrefe_*crosssec_/2.0;
        (*massmatrix)(2,2) = density*lrefe_*mominer_/2.0;
        (*massmatrix)(3,3) = density*lrefe_*crosssec_/2.0;
        (*massmatrix)(4,4) = density*lrefe_*crosssec_/2.0;
        (*massmatrix)(5,5) = density*lrefe_*mominer_/2.0;
       }
      else
        dserror("improper value of variable lumpedmass");
  }

  return;
} // DRT::ELEMENTS::Beam2r::nlnstiffmass


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_BEAM2R

