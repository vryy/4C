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
#ifdef D_BEAM2
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "beam2.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

//externally defined structure for material data
extern struct _MATERIAL  *mat;



/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                                 cyron 01/08|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam2::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::Beam2::ActionType act = Beam2::calc_none;
  // get the action required
  string action = params.get<string>("action","calc_none");
  if (action == "calc_none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")      act = Beam2::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")      act = Beam2::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = Beam2::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")  act = Beam2::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")  act = Beam2::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass") act = Beam2::calc_struct_nlnstifflmass;
  else if (action=="calc_struct_stress")        act = Beam2::calc_struct_stress;
  else if (action=="calc_struct_eleload")       act = Beam2::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")       act = Beam2::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")  act = Beam2::calc_struct_update_istep;
  else if (action=="calc_struct_update_imrlike") act = Beam2::calc_struct_update_imrlike;
  else if (action=="calc_struct_reset_istep")   act = Beam2::calc_struct_reset_istep;
  else dserror("Unknown type of action for Beam2");
   
  switch(act)
  {
    /*in case that only linear stiffness matrix is required b2_nlstiffmass is called with zero dispalcement and 
     residual values*/ 
     case Beam2::calc_struct_linstiff:
    {   
      //only nonlinear case implemented!
      dserror("linear stiffness matrix called, but not implemented");

    }
    break;
       
    //nonlinear stiffness and mass matrix are calculated even if only nonlinear stiffness matrix is required
    case Beam2::calc_struct_nlnstiffmass:
    case Beam2::calc_struct_nlnstifflmass:
    case Beam2::calc_struct_nlnstiff:
    case Beam2::calc_struct_internalforce:
    {
      int lumpedmass = lumpedflag_;  // 0=consistent, 1=lumped
      if (act==Beam2::calc_struct_nlnstifflmass) lumpedflag_ = 1;
      // need current global displacement and residual forces and get them from discretization
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      
      /*making use of the local-to-global map lm one can extract current displacemnet and residual values for
      / each degree of freedom */
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      
      if (act == Beam2::calc_struct_nlnstiffmass)
        b2_nlnstiffmass(mydisp,&elemat1,&elemat2,&elevec1);
      else if (act == Beam2::calc_struct_nlnstifflmass)
      {
        b2_nlnstiffmass(mydisp,&elemat1,&elemat2,&elevec1);
        lumpedflag_ = lumpedmass;
      }
      else if (act == Beam2::calc_struct_nlnstiff)
        b2_nlnstiffmass(mydisp,&elemat1,NULL,&elevec1);
      else if  (act ==  calc_struct_internalforce)
        b2_nlnstiffmass(mydisp,NULL,NULL,&elevec1);
      
    }
    break;
    case calc_struct_update_istep:
    {
      ;// there is nothing to do here at the moment
    }
    break;
    case calc_struct_update_imrlike:
    {
      ;// there is nothing to do here at the moment
    }
    break;
    case calc_struct_reset_istep:
    {
      ;// there is nothing to do here at the moment
    }
    break;
    default:
      dserror("Unknown type of action for Beam2 %d", act);
  }
  return 0;

}


/*-----------------------------------------------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)                                       cyron 01/08|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Beam2::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{		
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
  const DRT::UTILS::IntegrationPoints1D  intpoints = getIntegrationPoints1D(gaussrule_);
  
  
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
    const double xi = intpoints.qxg[ip];
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
  
  
/*by the following code part stochastic external forces can be applied elementwise; for decoupling 
 * load frequency and time step size stochastical forces should better be applied outside the element*/
  
  if (thermalenergy_ > 0)
  {	  
	  extern struct _MATERIAL  *mat;
	  // get the material law and density
	  MATERIAL* currmat = &(mat[material_-1]);
	  double density;
    
	  //assignment of material parameters; only St.Venant material is accepted for this beam 
	  switch(currmat->mattyp)
	    	{
		    case m_stvenant:// only linear elastic material supported
		      {
		    	  density = currmat->m.stvenant->density;
		      }
		      break;
		      default:
		      dserror("unknown or improper type of material law");
		 }	
	  
	  //calculating diagonal entry of damping matrix  
	  double gamma_trans = params.get<double>("damping factor M",0.0) * crosssec_ * density * lrefe_/2;
	  double gamma_rot   = params.get<double>("damping factor M",0.0) * mominer_  * density * lrefe_/2;
	  
	  //calculating standard deviation of statistical forces according to fluctuation dissipation theorem
	  double stand_dev_trans = pow(2 * thermalenergy_ * gamma_trans / params.get<double>("delta time",0.01),0.5);
	  double stand_dev_rot   = pow(2 * thermalenergy_ * gamma_rot   / params.get<double>("delta time",0.01),0.5);

	  //creating a random generator object which creates random numbers with mean = 0 and standard deviation
	  //stand_dev; using Blitz namespace "ranlib" for random number generation
	  ranlib::Normal<double> normalGen_trans(0,stand_dev_trans);
	  ranlib::Normal<double> normalGen_rot(0,stand_dev_rot);
	  
	  //adding statistical forces 
	  elevec1[0] += normalGen_trans.random();  
	  elevec1[1] += normalGen_trans.random();
	  elevec1[2] += normalGen_rot.random();
  	elevec1[3] += normalGen_trans.random();
  	elevec1[4] += normalGen_trans.random(); 
  	elevec1[5] += normalGen_rot.random();    	   
  }   
  return 0;
}

/*-----------------------------------------------------------------------------------------------------------*
 | evaluate auxiliary vectors and matrices for corotational formulation                           cyron 01/08|							     
 *----------------------------------------------------------------------------------------------------------*/
//notation for this function similar to Crisfield, Volume 1;
void DRT::ELEMENTS::Beam2::b2_local_aux(LINALG::SerialDenseMatrix& Bcurr,
                    			LINALG::SerialDenseVector& rcurr,
                    			LINALG::SerialDenseVector& zcurr,
                                        double& beta,
                    			const LINALG::SerialDenseMatrix& xcurr,
                    			const double& lcurr,
                    			const double& lrefe_)

{
	//this function expects x 3x2 matrix xcurr for 3 DOF of 2 beam2 nodes
	#ifdef DEBUG
	dsassert(xcurr.M()==3,"improper dimension of xcurr");
	dsassert(xcurr.N()==2,"improper dimension of xcurr");
	#endif // #ifdef DEBUG


  // beta is the rotation angle out of x-axis in a x-y-plane 
  double cos_beta = (xcurr(0,1)-xcurr(0,0))/lcurr;
  double sin_beta = (xcurr(1,1)-xcurr(1,0))/lcurr;
  
  /*a crucial point in a corotational frame work is how to calculate the correct rotational angle
   * beta, which is the rotation relative to the reference rotation beta0_; 
   * in case that abs(beta)<PI this can be done easily making use of asin(beta) and acos(beta);
   * however, since in the course of large deformations abs(beta)>PI can happen one needs some special means
   * in order to avoid confusion; here we calculate an angel psi in a local coordinate system rotated by an 
   * angle of halfrotations_*PI; afterwards we know beta = psi * halfrotations_*PI; for this procedure we only
   * have to make sure that abs(psi)>PI i.e. that the true rotation angle lies always in the range +/- PI of the 
   * rotated system. This can be achieved by updating the variable halfrotations_ after each time step in such
   * a way that the current angle leads to abs(psi)<0.5*PI in the coordinate system based on the updated variable
   * halfrotations_; the also in the next time step abs(psi)>PI will be guaranteed as long as no rotations through
   * an angle >0.5*PI take place within one time step throughout the whole simulated structure. This seems imply a 
   * fairly mild restriction to time step size*/
  
  // psi is the rotation angle out of x-axis in a x-y-plane rotatated by halfrotations_*PI
  double cos_psi = pow(-1.0,halfrotations_)*(xcurr(0,1)-xcurr(0,0))/lcurr;
  double sin_psi = pow(-1.0,halfrotations_)*(xcurr(1,1)-xcurr(1,0))/lcurr;
  
  //first we calculate psi in a range between -pi < beta <= pi in the rotated plane
  double psi = 0;
  if (cos_psi >= 0)
  	psi = asin(sin_psi);
  else
  {	if (sin_psi >= 0)
  		psi =  acos(cos_psi);
        else
        	psi = -acos(cos_psi);
   }
  
  beta = psi + PI*halfrotations_ - beta0_;
  
  if (psi > PI/2)
	  halfrotations_++;
  if (psi < -PI/2)
 	  halfrotations_--; 

  rcurr[0] = -cos_beta;
  rcurr[1] = -sin_beta;
  rcurr[2] = 0;
  rcurr[3] = cos_beta;
  rcurr[4] = sin_beta;
  rcurr[5] = 0;
  
  zcurr[0] = sin_beta;
  zcurr[1] = -cos_beta;
  zcurr[2] = 0;
  zcurr[3] = -sin_beta;
  zcurr[4] = cos_beta;
  zcurr[5] = 0;
    
  //assigning values to each element of the Bcurr matrix 
  Bcurr.Shape(3,6);
  
  for(int id_col=0; id_col<6; id_col++)
  	{
	  Bcurr(0,id_col) = rcurr[id_col];
	  Bcurr(2,id_col) = (lrefe_ / lcurr) * zcurr[id_col];
	  if (id_col == 2 || id_col ==5)
	    		Bcurr(2,id_col) = Bcurr(2,id_col) - (lrefe_ / 2);
  	}
    Bcurr(1,2) = 1;
    Bcurr(1,5) = -1;

  return;
} /* DRT::ELEMENTS::Beam2::b2_local_aux */

/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   cyron 01/08|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam2::b2_nlnstiffmass( vector<double>&           disp,
                                            Epetra_SerialDenseMatrix* stiffmatrix,
                                            Epetra_SerialDenseMatrix* massmatrix,
                                            Epetra_SerialDenseVector* force)
{
  const int numdf = 3;
  const int iel = NumNode();
  //coordinates in reference and current configuration of all the nodes in two dimensions stored in 3 x iel matrices
  LINALG::SerialDenseMatrix xrefe;
  LINALG::SerialDenseMatrix xcurr;

  xrefe.Shape(3,2);
  xcurr.Shape(3,2);

  //current length of beam in physical space
  double lcurr = 0;
  //current angle between x-axis and beam in physical space
  double beta;
  //some geometric auxiliary variables according to Crisfield, Vol. 1
  LINALG::SerialDenseVector zcurr;
  LINALG::SerialDenseVector rcurr;
  LINALG::SerialDenseMatrix Bcurr;
  //auxiliary matrix storing the product of constitutive matrix C and Bcurr
  LINALG::SerialDenseMatrix aux_CB;

  zcurr.Size(6);
  rcurr.Size(6);
  Bcurr.Shape(3,6);
  aux_CB.Shape(3,6);

  //calculating refenrence configuration xrefe and current configuration xcurr
  for (int k=0; k<iel; ++k)
  {

    xrefe(0,k) = Nodes()[k]->X()[0];
    xrefe(1,k) = Nodes()[k]->X()[1];
    xrefe(2,k) = Nodes()[k]->X()[2];

    xcurr(0,k) = xrefe(0,k) + disp[k*numdf+0];
    xcurr(1,k) = xrefe(1,k) + disp[k*numdf+1];
    //note: this is actually not the current director angle, but current director angle minus reference director angle
    xcurr(2,k) = xrefe(2,k) + disp[k*numdf+2];

  }
  
  // calculation of local geometrically important matrices and vectors; notation according to Crisfield-------
  //current length
  lcurr = pow( pow(xcurr(0,1)-xcurr(0,0),2) + pow(xcurr(1,1)-xcurr(1,0),2) , 0.5 );
  
  //calculation of local geometrically important matrices and vectors
  b2_local_aux(Bcurr, rcurr, zcurr, beta, xcurr, lcurr, lrefe_);
  
  //calculation of local internal forces
  LINALG::SerialDenseVector force_loc;
  
  force_loc.Size(3);
  
  /* read material parameters using structure _MATERIAL which is defined by inclusion of      /
  / "../drt_lib/drt_timecurve.H"; note: material parameters have to be read in the evaluation /
  / function instead of e.g. beam2_input.cpp or within the Beam2Register class since it is not/
  / sure that structure _MATERIAL is declared within those scopes properly whereas it is within/
  / the evaluation functions */

  // get the material law
  MATERIAL* currmat = &(mat[material_-1]);
  double ym;
  double sm;
  double density;
    
  //assignment of material parameters; only St.Venant material is accepted for this beam 
  switch(currmat->mattyp)
    	{
	    case m_stvenant:// only linear elastic material supported
	      {
	    	  ym = currmat->m.stvenant->youngs;
	    	  sm = ym / (2*(1 + currmat->m.stvenant->possionratio));
	    	  density = currmat->m.stvenant->density;
	      }
	      break;
	      default:
	      dserror("unknown or improper type of material law");
	 }	

  /*in the following internal forces are computed; one has to take into consideration that not global director angles cause
   * internal forces, but local ones. Global and local director angles are different from each other by the reference beam 
   * angle beta0_ so that theta_loc = theta_glob - beta0_; in case of curvature beta0_ cancels out of the calculation, however
   * in case of the internal shear force it has to be accounted for */
   
  //local internal axial force: note: application of the following expression of axial strain leads
  //to numerically significantly better stability in comparison with (lcurr - lrefe_)/lrefe_
  force_loc(0) = ym*crosssec_*(lcurr*lcurr - lrefe_*lrefe_)/(lrefe_*(lcurr + lrefe_));
  
  //local internal bending moment
  force_loc(1) = -ym*mominer_*(xcurr(2,1)-xcurr(2,0))/lrefe_;
  
  //local internal shear force   
  force_loc(2) = -sm*crosssecshear_*( (xcurr(2,1)+xcurr(2,0))/2 - beta - beta0_);  

  //calculating tangential stiffness matrix in global coordinates
  if (stiffmatrix != NULL)
  {
      
      //linear elastic part including rotation
      
      for(int id_col=0; id_col<6; id_col++)
      {
          aux_CB(0,id_col) = Bcurr(0,id_col) * (ym*crosssec_/lrefe_);
          aux_CB(1,id_col) = Bcurr(1,id_col) * (ym*mominer_/lrefe_);
          aux_CB(2,id_col) = Bcurr(2,id_col) * (sm*crosssecshear_/lrefe_);
      }
      
      (*stiffmatrix).Multiply('T','N',1,Bcurr,aux_CB,0);
          
      //adding geometric stiffness by shear force 
      double aux_Q_fac = force_loc(2)*lrefe_ / pow(lcurr,2);
      for(int id_lin=0; id_lin<6; id_lin++)
          for(int id_col=0; id_col<6; id_col++)
          {
              (*stiffmatrix)(id_lin,id_col) -= aux_Q_fac * rcurr(id_lin) * zcurr(id_col);
              (*stiffmatrix)(id_lin,id_col) -= aux_Q_fac * rcurr(id_col) * zcurr(id_lin);
          }
      
      //adding geometric stiffness by axial force 
      double aux_N_fac = force_loc(0)/lcurr; 
      for(int id_lin=0; id_lin<6; id_lin++)
          for(int id_col=0; id_col<6; id_col++)
              (*stiffmatrix)(id_lin,id_col) += aux_N_fac * zcurr(id_lin) * zcurr(id_col);
  }
  
  //calculation of global internal forces from force = B_transposed*force_loc 
  if (force != NULL)
  {
      (*force).Size(6);
      for(int id_col=0; id_col<6; id_col++)
          for(int id_lin=0; id_lin<3; id_lin++)
              (*force)(id_col) += Bcurr(id_lin,id_col)*force_loc(id_lin);
  }
  
  //calculating mass matrix (local version = global version) 
  if (massmatrix != NULL)
  {
      (*massmatrix).Shape(6,6);
      
      //if lumped_flag == 0 a consistent mass Timoshenko beam mass matrix is applied
      if (lumpedflag_ == 0)
      {
              //assignment of massmatrix by means of auxiliary diagonal matrix aux_E stored as an array
              double aux_E[3]={density*lrefe_*crosssec_/6,density*lrefe_*crosssec_/6,density*lrefe_*mominer_/6};
              for(int id=0; id<3; id++)
              {
              	    (*massmatrix)(id,id) = 2*aux_E[id];
                    (*massmatrix)(id+3,id+3) = 2*aux_E[id];
                    (*massmatrix)(id,id+3) = aux_E[id];
                    (*massmatrix)(id+3,id) = aux_E[id];
              }
      }
      /*if lumped_flag == 1 a lumped mass matrix is applied where the cross sectional moment of inertia is
       * assumed to be approximately zero so that the 3,3 and 5,5 element are both zero */
      
      else if (lumpedflag_ == 1)
      {
             (*massmatrix).Shape(6,6);
             //note: this is not an exact lumped mass matrix, but it is modified in such a way that it leads
             //to a diagonal mass matrix with constant diagonal entries
             (*massmatrix)(0,0) = density*lrefe_*crosssec_/2;
             (*massmatrix)(1,1) = density*lrefe_*crosssec_/2;	 
             (*massmatrix)(2,2) = density*lrefe_*mominer_/2; 
             (*massmatrix)(3,3) = density*lrefe_*crosssec_/2;
             (*massmatrix)(4,4) = density*lrefe_*crosssec_/2; 
             (*massmatrix)(5,5) = density*lrefe_*mominer_/2;
       }
      else
              dserror("improper value of variable lumpedflag_");    
  }
  return;
} // DRT::ELEMENTS::Beam2::b2_nlnstiffmass

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_BEAM2

