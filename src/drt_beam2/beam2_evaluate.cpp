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
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"

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
  else if (action=="calc_struct_stress")        act = Beam2::calc_struct_stress;
  else if (action=="calc_struct_eleload")       act = Beam2::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")       act = Beam2::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")  act = Beam2::calc_struct_update_istep;
  else dserror("Unknown type of action for Beam2");
   
  switch(act)
  {
    /*in case that only linear stiffness matrix is required b2_nlstiffmass is called with zero dispalcement and 
     residual values*/ 
     case Beam2::calc_struct_linstiff:
    {   
    
    /*lm is a vector mapping each local degree of freedom to a global one; its length is the overall number of 
    / degrees of freedom of the element */
      vector<double> mydisp(lm.size());        
      for (int i=0; i<(int)mydisp.size(); ++i) mydisp[i] = 0.0;  //setting displacement mydisp to zero
      
      //current residual forces
      vector<double> myres(lm.size());
      for (int i=0; i<(int)myres.size(); ++i) myres[i] = 0.0;    //setting residual myres to zero
          
      b2_nlnstiffmass(mydisp,elemat1,elemat2,elevec1);
    }
    break;
       
    //nonlinear stiffness and mass matrix are calculated even if only nonlinear stiffness matrix is required
    case Beam2::calc_struct_nlnstiffmass:
    case Beam2::calc_struct_nlnstiff:
    {
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
	      
      b2_nlnstiffmass(mydisp,elemat1,elemat2,elevec1);
    }
    break;
    case calc_struct_update_istep:
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
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);
  
  //reference length calculated from nodal data
  double x_sqdist_ref = pow( Nodes()[1]->X()[0] - Nodes()[0]->X()[0],2 );
  double y_sqdist_ref = pow( Nodes()[1]->X()[1] - Nodes()[0]->X()[1],2 );
  double length_refe = pow( x_sqdist_ref + y_sqdist_ref ,0.5 );

  
  //jacobian determinant
  double det = length_refe/2;


  // no. of nodes on this element
  const int iel = NumNode();
  const int numdf = 3;
  const DiscretizationType distype = this->Shape();

  // gaussian points 
  const DRT::UTILS::IntegrationPoints1D  intpoints = getIntegrationPoints1D(gaussrule_);
  
  
  //declaration of variable in order to store shape function
  Epetra_SerialDenseVector      funct(iel);

  // get values and switches from the condition
  const vector<int>*    onoff = condition.Get<vector<int> >("onoff");
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
    double ar[3];
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
	  double gamma = params.get<double>("damping factor M",0.0) * crosssec_ * density * length_refe;  
	  
	  //calculating standard deviation of statistical forces according to fluctuation dissipation theorem
	  double stand_dev = pow(2 * thermalenergy_ * gamma / params.get<double>("delta time",0.01),0.5);
	  
	  //using Blitz namespace for random number generation
	  using namespace ranlib;
	  //creating a random generator object which creates random numbers with mean = 0 and variance = 1
	  Normal<double> normalGen(0,1);

	  //adding statistical forces accounting for connectivity of nodes
	  elevec1[0] += stand_dev * normalGen.random() / sqrt( Nodes()[0]->NumElement() );  
	  elevec1[1] += stand_dev * normalGen.random() / sqrt( Nodes()[0]->NumElement() );
  	  elevec1[3] += stand_dev * normalGen.random() / sqrt( Nodes()[1]->NumElement() );
  	  elevec1[4] += stand_dev * normalGen.random() / sqrt( Nodes()[1]->NumElement() ); 	  	  
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
                    			const double& length_curr,
                    			const double& length_refe)

{
	//this function expects x 3x2 matrix xcurr for 3 DOF of 2 beam2 nodes
	#ifdef DEBUG
	dsassert(xcurr.M()==3,"improper dimension of xcurr");
	dsassert(xcurr.N()==2,"improper dimension of xcurr");
	#endif // #ifdef DEBUG


  // beta is the rotation angle out of x-axis in a x-y-plane
  double cos_beta = (xcurr(0,1)-xcurr(0,0))/length_curr;
  double sin_beta = (xcurr(1,1)-xcurr(1,0))/length_curr;
  
  
  //first we calculate beta in a range between -pi < beta <= pi
  if (cos_beta >= 0)
  	beta = asin(sin_beta);
  else
  {	if (sin_beta >= 0)
  		beta =  acos(cos_beta);
        else
        	beta = -acos(cos_beta);
   }
   /*then we take into consideration that we have to add 2n*pi to beta in order to get a realistic difference/
   / beta - teta1, where teta1 = xcurr(0,3); note that between diff_n1, diff_n2, diff_n3 there is such a great/
   / difference that it's no problem to lose the round off in the following cast operation */
   int teta1 = static_cast<int>( 360*xcurr(2,0)/(2*PI) );
   double n_aux = ( teta1 - (teta1 % 360) )/360;
   double diff_n1 = abs( beta + (n_aux-1) * 2 * PI - xcurr(2,0) );
   double diff_n2 = abs( beta + n_aux * 2 * PI - xcurr(2,0) );
   double diff_n3 = abs( beta + (n_aux+1) * 2 * PI - xcurr(2,0) );
   beta = beta + n_aux * 2 * PI;
   if (diff_n1 < diff_n2 && diff_n1 < diff_n3)
   	beta = beta + (n_aux-1) * 2 * PI;
   if (diff_n3 < diff_n2 && diff_n3 < diff_n1)
   	beta = beta + (n_aux+1) * 2 * PI;
  
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
	  Bcurr(2,id_col) = (length_refe / length_curr) * zcurr[id_col];
	  if (id_col == 2 || id_col ==5)
	    		Bcurr(2,id_col) = Bcurr(2,id_col) - (length_refe / 2);
  	}
    Bcurr(1,2) = 1;
    Bcurr(1,5) = -1;

  return;
} /* DRT::ELEMENTS::Beam2::b2_local_aux */

/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   cyron 01/08|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam2::b2_nlnstiffmass( vector<double>&           disp,
                                            Epetra_SerialDenseMatrix& stiffmatrix,
                                            Epetra_SerialDenseMatrix& massmatrix,
                                            Epetra_SerialDenseVector& force)
{
const int numdf = 3;
const int iel = NumNode();
//coordinates in reference and current configuration of all the nodes in two dimensions stored in 3 x iel matrices
LINALG::SerialDenseMatrix xrefe;
LINALG::SerialDenseMatrix xcurr;

xrefe.Shape(3,2);
xcurr.Shape(3,2);

//current length of beam in physical space
double length_curr = 0;
//length of beam in reference configuration
double length_refe = 0;
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
    xcurr(2,k) = xrefe(2,k) + disp[k*numdf+2];

  }
   
    
  // calculation of local geometrically important matrices and vectors; notation according to Crisfield-------
  //current length
  length_curr = pow( pow(xcurr(0,1)-xcurr(0,0),2) + pow(xcurr(1,1)-xcurr(1,0),2) , 0.5 );
  //length in reference configuration
  length_refe  = pow( pow(xrefe(0,1)-xrefe(0,0),2) + pow(xrefe(1,1)-xrefe(1,0),2) , 0.5 );
  
  //calculation of local geometrically important matrices and vectors
  b2_local_aux(Bcurr, rcurr, zcurr, beta, xcurr, length_curr, length_refe);
  
  //calculation of local internal forces
  LINALG::SerialDenseVector force_loc;
  
  force_loc.Size(3);
  
  /* read material parameters using structure _MATERIAL which is defined by inclusion of      /
  / "../drt_lib/drt_timecurve.H"; note: material parameters have to be read in the evaluation /
  / function instead of e.g. beam2_input.cpp or within the Beam2Register class since it is not/
  / sure that structure _MATERIAL is declared within those scopes properly whereas it is within/
  / the evaluation functions */
  
  extern struct _MATERIAL  *mat;
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

  
  //local internal axial force
  force_loc(0) = ym*crosssec_*(length_curr - length_refe)/length_refe;
  
  //local internal bending moment
  force_loc(1) = -ym*mominer_*(xcurr(2,1)-xcurr(2,0))/length_refe;
  
  //local internal shear force
  force_loc(2) = -sm*crosssecshear_*( (xcurr(2,1)+xcurr(2,0))/2 - beta);
 

  //calculating tangential stiffness matrix in global coordinates---------------------------------------------
  
  //linear elastic part including rotation
  
  for(int id_col=0; id_col<6; id_col++)
  {
	  aux_CB(0,id_col) = Bcurr(0,id_col) * (ym*crosssec_/length_refe);
	  aux_CB(1,id_col) = Bcurr(1,id_col) * (ym*mominer_/length_refe);
	  aux_CB(2,id_col) = Bcurr(2,id_col) * (sm*crosssecshear_/length_refe);
  }
   
  stiffmatrix.Multiply('T','N',1,Bcurr,aux_CB,0);

  //adding geometric stiffness by shear force 
  double aux_Q_fac = force_loc(2)*length_refe / pow(length_curr,2);
  for(int id_lin=0; id_lin<6; id_lin++)
  	for(int id_col=0; id_col<6; id_col++)
  	{
  		stiffmatrix(id_lin,id_col) += aux_Q_fac * rcurr(id_lin) * zcurr(id_col);
  		stiffmatrix(id_lin,id_col) += aux_Q_fac * rcurr(id_col) * zcurr(id_lin);
  	}
  
  //adding geometric stiffness by axial force 
  double aux_N_fac = force_loc(1)/length_curr; 
  for(int id_lin=0; id_lin<6; id_lin++)
  	for(int id_col=0; id_col<6; id_col++)
  		stiffmatrix(id_lin,id_col) += aux_N_fac * zcurr(id_lin) * zcurr(id_col);  
  
  //calculation of global internal forces from force = B_transposed*force_loc---------------------------------- 
  for(int id_col=0; id_col<6; id_col++)
	  for(int id_lin=0; id_lin<3; id_lin++)
    	force(id_col) += Bcurr(id_lin,id_col)*force_loc(id_lin);
  
  //calculating mass matrix (lcoal version = global version)-------------------------------------------------- 
  //the following code lines are based on the assumption that massmatrix is a 6x6 matrix filled with zeros
  #ifdef DEBUG
  dsassert(massmatrix.M()==6,"wrong mass matrix input");
  dsassert(massmatrix.N()==6,"wrong mass matrix input");
  for(int i=0; i<6; i++)
  	for(int j=0; j<6; j++)
  	dsassert(massmatrix(i,j)==0,"wrong mass matrix input in beam2_evaluate, function b2_nlnstiffmass"); 
  #endif // #ifdef DEBUG
  
  
  //if lumped_flag == 0 a consistent mass Timoshenko beam mass matrix is applied
  if (lumpedflag_ == 0)
  {
	  //assignment of massmatrix by means of auxiliary diagonal matrix aux_E stored as an array
	  double aux_E[3]={density*length_refe*crosssec_/6,density*length_refe*crosssec_/6,density*length_refe*mominer_/6};
	  for(int id=0; id<3; id++)
	  {
	  	massmatrix(id,id) = 2*aux_E[id];
	        massmatrix(id+3,id+3) = 2*aux_E[id];
	        massmatrix(id,id+3) = aux_E[id];
	        massmatrix(id+3,id) = aux_E[id];
	  }
  }
  /*if lumped_flag == 1 a lumped mass matrix is applied where the cross sectional moment of inertia is
   * assumed to be approximately zero so that the 3,3 and 5,5 element are both zero */
  
  else if (lumpedflag_ == 1)
  {
 	 massmatrix.Shape(6,6);
 	 massmatrix(0,0) = density*length_refe*crosssec_/2;
 	 massmatrix(1,1) = density*length_refe*crosssec_/2;
 	 massmatrix(3,3) = density*length_refe*crosssec_/2;
 	 massmatrix(4,4) = density*length_refe*crosssec_/2;
   }
  else
	  dserror("improper value of variable lumpedflag_");

  return;
} // DRT::ELEMENTS::Beam2::b2_nlnstiffmass(

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_BEAM2
