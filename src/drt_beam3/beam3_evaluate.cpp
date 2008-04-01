/*!-----------------------------------------------------------------------------------------------------------
\file beam2_evaluate.cpp
\brief

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
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"

//externally defined structure for material data
extern struct _MATERIAL  *mat;



/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                                 cyron 01/08|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
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
  else if (action=="calc_struct_linstiff")      act = Beam3::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")      act = Beam3::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = Beam3::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")  act = Beam3::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")  act = Beam3::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_stress")        act = Beam3::calc_struct_stress;
  else if (action=="calc_struct_eleload")       act = Beam3::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")       act = Beam3::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")  act = Beam3::calc_struct_update_istep;
  else if (action=="calc_struct_update_genalpha_imrlike") act = Beam3::calc_struct_update_genalpha_imrlike;
  else dserror("Unknown type of action for Beam3");
   
  switch(act)
  {
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
    case Beam3::calc_struct_nlnstiff:
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
        
      b3_nlnstiffmass(mydisp,elemat1,elemat2,elevec1);
      

      //the following code block can be used to check quickly whether the nonlinear stiffness matrix is calculated
      //correctly or not: the function b3_nlnstiff_approx(mydisp) calculated the stiffness matrix approximated by
      //finite differences and prints the relative error of every single element in case that it is >1e-5; before
      //activating this part of code also the function b3_nlnstiff_approx(mydisp) has to be activated both in Beam3.H
      //and Beam3_evaluate.cpp
      /*
      
      Epetra_SerialDenseMatrix stiff_approx;
      Epetra_SerialDenseMatrix stiff_relerr;
      stiff_approx.Shape(6,6);
      stiff_relerr.Shape(6,6);      
      double h_rel = 1e-7;
      int outputflag = 0;
      stiff_approx = b3_nlnstiff_approx(mydisp, h_rel);
      
      for(int line=0; line<6; line++)
      {
	      for(int col=0; col<6; col++)
		      {
		      	stiff_relerr(line,col)= abs( (elemat1(line,col) - stiff_approx(line,col))/elemat1(line,col) );
		      	if (stiff_relerr(line,col)<h_rel*100)
		      		//suppressing small entries whose effect is only confusing
		      		stiff_relerr(line,col)=0;
		      	//there is no error if an entry is nan e.g. due to dirichlet boundary conditions
		      	if ( isnan( stiff_relerr(line,col) ) )
		      		stiff_relerr(line,col)=0;
		      	if (stiff_relerr(line,col)>0)
		      	   outputflag = 1;		      	
		      }
      }
	if(outputflag ==1)
	{
	      std::cout<<"\n\n acutally calculated stiffness matrix"<< elemat1;
	      std::cout<<"\n\n approximated stiffness matrix"<< stiff_approx;    
	      std::cout<<"\n\n rel error stiffness matrix"<< stiff_relerr;
	} 
	*/
    
    }
    break;
    case calc_struct_update_istep:
    {
      ;// there is nothing to do here at the moment
    }
    break;
    case calc_struct_update_genalpha_imrlike:
    {
      ;// there is nothing to do here at the moment
    }
    break;
    default:
      dserror("Unknown type of action for Beam3 %d", act);
  }
  return 0;

}


/*-----------------------------------------------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)                                       cyron 03/08|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Beam3::EvaluateNeumann(ParameterList& params,
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
  if (curvenum>=0 && usetime)//notation for this function similar to Crisfield, Volume 1;
    curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);
  
  //jacobian determinant
  double det = lrefe_/2;


  // no. of nodes on this element; the following line is only valid for elements with constant number of 
  // degrees of freedom per node
  const int numdf = Nodes()[0]->NumDofPerNode();
  const DiscretizationType distype = this->Shape();

  // gaussian points 
  const DRT::UTILS::IntegrationPoints1D  intpoints = getIntegrationPoints1D(gaussrule_);
  
  
  //declaration of variable in order to store shape function
  Epetra_SerialDenseVector      funct(NumNode());

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
    for (int node=0; node<NumNode(); ++node)
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
  	  double gamma;
  	  gamma = 4*params.get<double>("damping factor M",0.0) * crosssec_ * density * lrefe_/3;  
  	  
  	  //calculating standard deviation of statistical forces according to fluctuation dissipation theorem
  	  double stand_dev = pow(2 * thermalenergy_ * gamma / params.get<double>("delta time",0.01),0.5);

  	  //creating a random generator object which creates random numbers with mean = 0 and standard deviation
  	  //stand_dev; using Blitz namespace "ranlib" for random number generation
  	  ranlib::Normal<double> normalGen(0,stand_dev);
  	  
  	  //adding statistical forces accounting for connectivity of nodes
  	(int node=0; node<NumNode(); ++node)
  	  for (int idof=0; idof<numdf; ++numdf)  		  
  	  {
  	    for (int inode=0; inode<NumNode(); ++inode)
  		  elevec1[numdf+numdf*inode] += normalGen.random() / sqrt( Nodes()[inode]->NumElement() ); 
  	  }	   
    }  

return 0;
}



/*-----------------------------------------------------------------------------------------------------------*
 | auxiliary functions for dealing with large rotations and nonlinear stiffness                    cyron 04/08|							     
 *----------------------------------------------------------------------------------------------------------*/
//computing spin matrix out of a rotation vector
void computespin(Epetra_SerialDenseMatrix& spin, const Epetra_SerialDenseMatrix rotationangle)
{
  spin.Shape(3,3);
  spin(0,1) = -rotation(2,0);
  spin(1,0) =  rotation(2,0);
  spin(0,2) =  rotation(1,0);
  spin(2,0) = -rotation(1,0);
  spin(1,2) = -rotation(0,0);
  spin(2,1) =  rotation(0,0); 
  return;
} /* DRT::ELEMENTS::Beam3::computespin */

//computing rotation matrix out of a spin matrix
void computerotation(Epetra_SerialDenseMatrix& rotationmatrix, const Epetra_SerialDenseMatrix spin)
{
  rotationmatrix.Shape(3,3);
  rotationmatrix(0,0)=1;
  rotationmatrix(1,1)=1;
  rotationmatrix(2,2)=1;
  rotationmatrix += spin;
  return;
} /* DRT::ELEMENTS::Beam3::computerotation */

//!computing rotation matrix X according to Crisfield, Vol. 2, equation (17.74)
void computeXrot(Epetra_SerialDenseMatrix& Xrot, const Epetra_SerialDenseMatrix T_new, const Epetra_SerialDenseMatrix x21)
{
  x21.Scale(0.5);
  Epetra_SerialDenseMatrix spinx21;
  computespin(spinx21,x21);
  
  Xrot.Shape(12,6);
  for (int j=0; j<6; ++j)
     {
        Xrot(j,j)   = -1;
        Xrot(j+6,j) =  1;
     } 
  
  for (int j=0; j<3; ++j)
     {
        Xrot(j+3,j) = spinx21(j,j);
        Xrot(j+9,j) = spinx21(j,j);
     } 
  
  return;
} /* DRT::ELEMENTS::Beam3::computeXrot*/

//computing stiffens matrix Ksigma1 according to Crisfield, Vol. 2, equation (17.83)
void computeKsig1(Epetra_SerialDenseMatrix& Ksig1, const Epetra_SerialDenseMatrix stressn, const Epetra_SerialDenseMatrix stressm)
{
  Ksig1.Shape(12,12);
  stressn.Scale(0.5);
  stressm.Scale(0.5);
  Epetra_SerialDenseMatrix Sn;
  Epetra_SerialDenseMatrix Sm;
  computespin(Sn,stressn);
  computespin(Sm,stressm);
  
  for (int i=0; i<3; ++i)
    {
      for (int j=0; j<3; ++j)
      {
        Ksig1(i  ,j+3) =  Sn(i,j);
        Ksig1(i  ,j+9) =  Sn(i,j);
        Ksig1(i+6,j+3) = -Sn(i,j);
        Ksig1(i+6,j+9) = -Sn(i,j);
        
        Ksig1(i+3,j+3) =  Sm(i,j);
        Ksig1(i+3,j+9) =  Sm(i,j);
        Ksig1(i+9,j+3) = -Sm(i,j);
        Ksig1(i+9,j+9) = -Sm(i,j);
      }
    }   
  return;
} /* DRT::ELEMENTS::Beam3::computeKsig1*/

//computing stiffens matrix Ksigma1 according to Crisfield, Vol. 2, equation (17.87) and (17.88)
void computeKsig2(Epetra_SerialDenseMatrix& Ksig2, const Epetra_SerialDenseMatrix stressn, const Epetra_SerialDenseMatrix x21)
{
  Ksig2.Shape(12,12);
  
  stressn.Scale(0.5);
  x21.Scale(0.25);
  Epetra_SerialDenseMatrix Sn;
  computespin(Sn,stressn);
  Epetra_SerialDenseMatrix Y;
  Y.Shape(3,3);
  computespin(Y,x21)
  Y.Multiply('N','N',1,Y,Sn,0); 
  
  for (int i=0; i<3; ++i)
    {
      for (int j=0; j<3; ++j)
      {
        Ksig2(i+3,j  ) = -Sn(i,j);
        Ksig2(i+3,j+6) =  Sn(i,j);
        Ksig2(i+9,j  ) = -Sn(i,j);
        Ksig2(i+9,j+6) =  Sn(i,j);
        
        Ksig2(i+3,j+3) =  Y(i,j);
        Ksig2(i+3,j+9) =  Y(i,j);
        Ksig2(i+9,j+3) =  Y(i,j);
        Ksig2(i+9,j+9) =  Y(i,j);
      }
    }   
  return;
} /* DRT::ELEMENTS::Beam3::computeKsig2*/



/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   cyron 01/08|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3::b3_nlnstiffmass( vector<double>&           disp,
                                            Epetra_SerialDenseMatrix& stiffmatrix,
                                            Epetra_SerialDenseMatrix& massmatrix,
                                            Epetra_SerialDenseVector& force)
{
  //normal and shear strain
  Epetra_SerialSendeMatrix espilon;
  epsilon.Shape(3,1);
  //number of nodal degrees of freedom
  const int numdf = Nodes()[0]->NumDofPerNode(); 
  //current position of nodal degrees of freedom
  Epetra_SerialDenseMatrix xcurr;
  xcurr.Shape(numdf,NumNode());
  //rotation matrix XT, Crisfield, Vol. 2, equation (17.74)
  Epetra_SerialDenseMatrix Xrot;
  Xrot.Shape(2*numdf,numdf);
  //stress values n and m, Crisfield, Vol. 2, equation (17.78)
  Epetra_SerialDenseMatrix stressn;
  stressn.Shape(3,1);
  Epetra_SerialDenseMatrix stressm;
  stressm.Shape(3,1);
  //nonlinear parts of stiffness matrix, Crisfiel Vol. 2, equation (17.83) and (17.87)
  Epetra_SerialDenseMatrix Ksig1;
  Epetra_SerialDenseMatrix Ksig2;
  //rotation matrix, Crisfield Vol. 2, equation (17.74)
  Epetra_SerialDenseMatrix XT;
  XT.Shape(12,6);
  
   
  /*update of central nodal triad T; in case that a time step was finished after the last call 
   * of Evaluate formerly Tnew_ becomes Told_ */
  if(params.get("total time", 0) > timenew_ + params.get<double>("delta time",0.01),0.5)/2;
  {                     
    //formerly "new" variables are now the "old" ones
    Told_= Tnew_;
    curvold_ = curvnew_;   
    //new time is current time
    timenew_ = params.get("total time", 0);
    betaplusalphaold_  = betaplusalphanew_;
    betaminusalphaold_ = betaminusalphanew_;
  }
  
  //nodal coordinates in current position
  for (int k=0; k<NumNode(); ++k)   
    {
      for (int j=0; j<numdf; ++j)
      {
        xcurr(j,k) = Nodes()[k]->X()[j] + disp[k*numdf+j]; 
      } 
    }
  
  //first of all "new" variables have to be adopted to dispalcement passed in from BACI driver
  
  //difference between coordinates of both nodes, x21' Crisfield  Vol. 2 equ. (17.66a) and (17.72)
   for (int j=0; j<3; ++j)
   {
     x21(j) = xcurr(j,1) - xcurr(j,0); 
     betaplusalphanew_(j)  = xcurr(j+3,1) + xcurr(j+3,0); 
     betaminusalphanew_(j) = xcurr(j+3,1) - xcurr(j+3,0); 
   } 
  
  //auxiliary matrix for update of Tnew_, Crisfield, Vol 2, equation (17.65)
  computespin(DeltaT,Deltaalpha);
  Tnew_.Multiply('N','N',1,(alphaplusbetanew_-alphaplusbetaold_)/2,Told_,0); 
  
  //computing triad for curvature update, Crisfield, Vol 2, equation (17.73)
  Tmid_.Multiply('N','N',1,computerotation(DeltaT.Scale(0.5)),Told_,0); 
  
  //updating curvature, Crisfield, Vol. 2, equation (17.72)
  curvnew_.Multiply('T','N',1,Tmid_,betaminusalphanew_-betaminusalphaold_,0); 
  curvnew_.Scale(1/lrefe_);
  curvnew_ += curvold_;
  
  //computing current axial and shear strain epsilon, Crisfield, Vol. 2, equation (17.67)
  epsilon.Multiply('T','N',1,Tnew_,x21,0); 
  epsilon.Sacle(1/lrefe_);
  epsilon(1) = epsilon(1) - 1;
  
  //computing rotation matrix XT according to Crisfield, Vol. 2, equation (17.74)
  void computeXrot(Epetra_SerialDenseMatrix& Xrot, const Epetra_SerialDenseMatrix T_new, const Epetra_SerialDenseMatrix x21);
  
  /* read material parameters using structure _MATERIAL which is defined by inclusion of      /
  / "../drt_lib/drt_timecurve.H"; note: material parameters have to be read in the evaluation /
  / function instead of e.g. Beam3_input.cpp or within the Beam3Register class since it is not/
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
    
  //stress values n and m, Crisfield, Vol. 2, equation (17.76) and (17.78)
  stressn = epsilon;
  stressn(0) = stressn(0)*ym*crosssec_;
  stressn(1) = stressn(1)*sm*crosssecshear_;
  stressn(2) = stressn(2)*sm*crosssecshear_;
  stressn.Multiply('N','N',1,Tnew_,stressn,0); 

  stressm = curvnew_;
  stressm(0) = stressm(0)*sm*Irr_;
  stressm(1) = stressm(1)*ym*Iyy_;
  stressm(2) = stressm(2)*ym*Izz_;
  stressm.Multiply('N','N',1,Tnew_,stressm,0); 
  
  //computing global internal forces, Crisfield Vol. 2, equation (17.79)
  force.Size(12);
  for (int i=0; i<12; ++i)
    {
      for (int j=0; j<3; ++j)
      {
        force(i) += Xrot(i,j)*stressn(j)
      }
      for (int j=0; j<3; ++j)
      {
        force(i) += Xrot(i,j+3)*stressm(j)
      }
    } 
   
  //stress dependent nonlinear parts of the stiffness matrix according to Crisfield, Vol. 2 equs. (17.83) and (17.87)
  computeKsig1(Ksig1,stressn,stressm);
  computeKsig2(Ksig2,stressn,x21);
  
  //computing linear stiffness matrix
  stiffmatrix.Shape(6,6);
  stiffmatrix(0,0) = ym*crosssec_/lrefe_; 
  stiffmatrix(1,1) = sm*crosssecshear_/lrefe_; 
  stiffmatrix(2,2) = sm*crosssecshear_/lrefe_;  
  stiffmatrix(3,3) = sm*Irr_/lrefe_; 
  stiffmatrix(4,4) = ym*Iyy_/lrefe_; 
  stiffmatrix(5,5) = ym*Izz_/lrefe_; 
  XT.Multiply('N','N',1,Xrot,Tnew_,0); 
  stiffmatrix.Multiply('N','N',1,XT,stiffmatrix,0); 
  stiffmatrix.Multiply('N','T',1,stiffmatrix,XT,0); 
  
  //adding nonlinear parts to tangent stiffness matrix, Crisfield, Vol. 2, equation (17.89)
  stiffmatrix += Ksig1;
  stiffmatrix += Ksig2;
   
  //calculating mass matrix (lcoal version = global version) 
  massmatrix.Shape(6,6);
  
  //if lumped_flag == 0 a consistent mass Timoshenko beam mass matrix is applied
  if (lumpedflag_ == 0)
  {
	  //assignment of massmatrix by means of auxiliary diagonal matrix aux_E stored as an array
	  double aux_E[3]={density*lrefe_*crosssec_/6,density*lrefe_*crosssec_/6,density*lrefe_*Iyy_/6};
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
 	 //note: this is not an exact lumped mass matrix, but it is modified in such a way that it leads
 	 //to a diagonal mass matrix with constant diagonal entries
 	 massmatrix(0,0) = 4*density*lrefe_*crosssec_/( 3*Nodes()[0]->NumElement() );
 	 massmatrix(1,1) = 4*density*lrefe_*crosssec_/( 3*Nodes()[0]->NumElement() );
 	 
 	 massmatrix(2,2) = 4*density*lrefe_*crosssec_/( 3*Nodes()[0]->NumElement() );
 	 
 	 massmatrix(3,3) = 4*density*lrefe_*crosssec_/( 3*Nodes()[1]->NumElement() );
 	 massmatrix(4,4) = 4*density*lrefe_*crosssec_/( 3*Nodes()[1]->NumElement() );
 	 
 	 massmatrix(5,5) = 4*density*lrefe_*crosssec_/( 3*Nodes()[1]->NumElement() );
   }
  else
	  dserror("improper value of variable lumpedflag_");  
    
  return;
} // DRT::ELEMENTS::Beam3::b3_nlnstiffmass

//the following function can be activated in order to find bugs; it calculates a finite difference
//approximation of the nonlinear stiffness matrix; activate the follwing block for bug fixing only

Epetra_SerialDenseMatrix DRT::ELEMENTS::Beam3::b3_nlnstiff_approx(vector<double>& disp, double h_rel)
{	
	Epetra_SerialDenseMatrix stiff_dummy;
	stiff_dummy.Shape(6,6);
	Epetra_SerialDenseMatrix mass_dummy;
	mass_dummy.Shape(6,6);
	Epetra_SerialDenseVector force_disp;
	force_disp.Size(6);
	Epetra_SerialDenseVector force_disp_delta;
	force_disp_delta.Size(6);
	Epetra_SerialDenseMatrix stiff_approx;
	stiff_approx.Shape(6,6);
	vector<double> disp_delta;
	
	DRT::ELEMENTS::Beam3::b3_nlnstiffmass(disp,stiff_dummy,mass_dummy,force_disp);
	
	for(int col=0; col<6; col++)
	{		
		force_disp_delta.Size(6);
		disp_delta = disp;
		disp_delta[col] = disp_delta[col]+h_rel;
		DRT::ELEMENTS::Beam3::b3_nlnstiffmass(disp_delta,stiff_dummy,mass_dummy,force_disp_delta);
		for(int line=0; line<6; line++)
			stiff_approx(line,col) = (force_disp_delta[line] - force_disp[line])/h_rel;		
	} 
	return stiff_approx;	
}


void DRT::ELEMENTS::Beam3::Arbeit(double& AN,double& AM,double& AQ, double& xv)
  {
  	AN += Arbeit_N;
  	AM += Arbeit_M;
  	AQ += Arbeit_Q;
  	xv = x_verschiebung;
	return;
  } 
void DRT::ELEMENTS::Beam3::Thermik(double& kT)
  {	
	kT = thermalenergy_;
	return;
  } 

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_BEAM3
