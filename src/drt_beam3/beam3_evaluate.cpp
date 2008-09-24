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

//externally defined structure for material data
extern struct _MATERIAL *mat;

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

      if (act == Beam3::calc_struct_nlnstiffmass)
      {
        b3_nlnstiffmass(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
      }
      else if (act == Beam3::calc_struct_nlnstifflmass)
      {
        b3_nlnstiffmass(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
        // lump mass matrix (bborn 07/08)
        // the mass matrix is lumped anyway, cf #b3_nlnstiffmass
        //b3_lumpmass(&elemat2);    
      }
      else if (act == Beam3::calc_struct_nlnstiff)
        b3_nlnstiffmass(params,myvel,mydisp,&elemat1,NULL,&elevec1);     
        
      else if (act == Beam3::calc_struct_internalforce)
      b3_nlnstiffmass(params,myvel,mydisp,NULL,NULL,&elevec1);  
      
      /*at the end of an iteration step the geometric ocnfiguration has to be updated: the starting point for the
       * next iteration step is the configuration at the end of the current step */  
      Qold_ = Qnew_;
      curvold_ = curvnew_;
      betaplusalphaold_ = betaplusalphanew_;
      betaminusalphaold_ = betaminusalphanew_; 
      
      /*
      //the following code block can be used to check quickly whether the nonlinear stiffness matrix is calculated
      //correctly or not by means of a numerically approximated stiffness matrix
      if(Id() == 3) //limiting the following tests to certain element numbers
      {       
        //variable to store numerically approximated stiffness matrix
        Epetra_SerialDenseMatrix stiff_approx;       
        stiff_approx.Shape(12,12);
        
        //relative error of numerically approximated stiffness matrix
        Epetra_SerialDenseMatrix stiff_relerr;         
        stiff_relerr.Shape(12,12);
        
        //characteristic length for numerical approximation of stiffness
        double h_rel = 1e-6;
        
        //flag indicating whether approximation lead to significant relative error
        int outputflag = 0;
        
        //calculating strains in new configuration
        for(int i=0; i<6; i++)
        {
          for(int k=0; k<2; k++)
          {          
            Epetra_SerialDenseVector force_aux;
            force_aux.Size(12);
            
            //create new displacement and velocity vectors in order to store artificially modified displacements
            vector<double> vel_aux(12);
            vector<double> disp_aux(12);
            for(int id = 0;id<12;id++)
            {
                DRT::UTILS::ExtractMyValues(*disp,disp_aux,lm);
                DRT::UTILS::ExtractMyValues(*vel,vel_aux,lm);
            }
            
            //modifying displacment artificially (for numerical derivative of internal forces):
            disp_aux[i + 6*k] = disp_aux[i + 6*k] + h_rel;
             vel_aux[i + 6*k] =  vel_aux[i + 6*k] + h_rel * params.get<double>("gamma",0.581) / ( params.get<double>("delta time",0.01)*params.get<double>("beta",0.292) );
            
            b3_nlnstiffmass(params,vel_aux,disp_aux,NULL,NULL,&force_aux);
                      
            for(int u = 0;u<12;u++)
            {
              stiff_approx(u,i+k*6)= ( pow(force_aux[u],2) - pow(elevec1(u),2) )/ (h_rel * (force_aux[u] + elevec1(u) ) );
            }       
          }
        }

       for(int line=0; line<12; line++)
       {
         for(int col=0; col<12; col++)
         {
           stiff_relerr(line,col)= fabs( ( pow(elemat1(line,col),2) - pow(stiff_approx(line,col),2) )/ ( (elemat1(line,col) + stiff_approx(line,col)) * elemat1(line,col) ));
                    
           //suppressing small entries whose effect is only confusing and NaN entires (which arise due to zero entries)
           if ( fabs( stiff_relerr(line,col) ) < h_rel*100 || isnan( stiff_relerr(line,col) ))
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
      } 
      //end of section in which numerical approximation for stiffness matrix is computed
      */

      
    
    }
    break;
    case calc_struct_update_istep:
    case calc_struct_update_imrlike:
    {
      /*the action calc_struct_update_istep is called in the very end of a time step when the new dynamic
       * equilibrium has finally been found; this is the point where the variable representing the geomatric
       * status of the beam have to be updated; the geometric status is represented by means of the triad Tnew_,
       * the curvature curvnew_ and the angular values betaplusalphanew_ and betaminusalphanew_*/
      Qconv_ = Qnew_;
      curvconv_ = curvnew_;
      betaplusalphaconv_ = betaplusalphanew_;
      betaminusalphaconv_ = betaminusalphanew_; 
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
      betaplusalphaold_ = betaplusalphaconv_;
      betaminusalphaold_ = betaminusalphaconv_; 
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

  //jacobian determinant
  double det = lrefe_/2;

  // no. of nodes on this element; the following line is only valid for elements with constant number of 
  // degrees of freedom per node
  const int numdf = 6;
  const DiscretizationType distype = this->Shape();

  // gaussian points 
  const DRT::UTILS::IntegrationPoints1D intpoints = getIntegrationPoints1D(gaussrule_);

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

  //stochastic forces due to fluctuation-dissipation theorem 
  if (kT_ > 0)
  {     
    //stochastic field of line load is interpolated by zeroth order polynomial functions
    if (stochasticorder_ == 0)
    {
      //calculating standard deviation of statistical forces according to fluctuation dissipation theorem
      double stand_dev_trans = pow(2 * kT_ * (zeta_/2) / params.get<double>("delta time",0.01),0.5);
  
      //creating a random generator object which creates random numbers with mean = 0 and standard deviation
      //stand_dev; using Blitz namespace "ranlib" for random number generation
      ranlib::Normal<double> normalGen(0,stand_dev_trans);
      
      //adding statistical forces 
      elevec1[0] += normalGen.random();  
      elevec1[1] += normalGen.random();
      elevec1[2] += normalGen.random();
      elevec1[6] += normalGen.random();
      elevec1[7] += normalGen.random();  
      elevec1[8] += normalGen.random();   
    }
    //stochastic field of line load is interpolated by first order polynomial functions
    else if (stochasticorder_ == 1)
    {
      //calculating standard deviation of statistical forces according to fluctuation dissipation theorem
      double stand_dev_trans = pow(2 * kT_ * (zeta_/3) / params.get<double>("delta time",0.01),0.5);
  
      //creating a random generator object which creates random numbers with mean = 0 and standard deviation
      //stand_dev; using Blitz namespace "ranlib" for random number generation
      ranlib::Normal<double> normalGen(0,stand_dev_trans);
      
      //adding uncorrelated components of statistical forces 
      elevec1[0] += normalGen.random();  
      elevec1[1] += normalGen.random();
      elevec1[2] += normalGen.random();
      elevec1[6] += normalGen.random();
      elevec1[7] += normalGen.random();  
      elevec1[8] += normalGen.random();   
      
      //adding correlated components of statistical forces 
      double force1 = normalGen.random()/pow(2,0.5);
      double force2 = normalGen.random()/pow(2,0.5);
      double force3 = normalGen.random()/pow(2,0.5);
      elevec1[0] += force1;  
      elevec1[1] += force2;
      elevec1[2] += force3;
      elevec1[6] += force1;
      elevec1[7] += force2;
      elevec1[8] += force3;     
    }       
  }  
  return 0;
}

/*-----------------------------------------------------------------------------------------------------------*
 | auxiliary functions for dealing with large rotations and nonlinear stiffness                    cyron 04/08|							     
 *----------------------------------------------------------------------------------------------------------*/
//computing basis of stiffness matrix of Crisfield, Vol. 2, equation (17.81)
inline void DRT::ELEMENTS::Beam3::computestiffbasis(const BlitzMat3x3& Tnew, const BlitzVec3& Cm, const BlitzVec3& Cb, const BlitzMat3x3& spinx21, Epetra_SerialDenseMatrix& stiffmatrix)
{
  //calculating the first matrix of (17.81) directly involves multiplications of large matrices (e.g. with the 12x6-matrix X)
  //application of the definitions in (17.74) allows blockwise evaluation with multiplication and addition of 3x3-matrices only
  //in the follwoing all the blocks are treated separately; their name is related directly to their content:
  //e.g. TCmTt is the product of the 3 matrices T * C_m * T^t (with T and C_m according to (17.74) and (17.76)
  //for the blockwise calculation on which the following steps are based on the relation S^t = -S for spin matrices was applied
  
  BlitzMat3x3 TCmTt;
  BlitzMat3x3 TCbTt;
  BlitzMat3x3 STCmTt;
  BlitzMat3x3 STCmTtSt;
  
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
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      STCmTt(i,j) = 0.0;
      for (int k = 0; k < 3; ++k)
      {
        STCmTt(i,j) += spinx21(i,k)*TCmTt(k,j);
      }
    }
  }
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      STCmTtSt(i,j) = 0.0;
      for (int k = 0; k < 3; ++k)
      {
        STCmTtSt(i,j) += STCmTt(i,k)*spinx21(j,k);
      }
    }
  }
  //calculating basis of stiffness matrix by means of above blocks
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      stiffmatrix(i  ,j)   =  TCmTt(i,j);
      stiffmatrix(i  ,j+3) =  STCmTt(j,i);
      stiffmatrix(i  ,j+6) = -TCmTt(i,j);
      stiffmatrix(i  ,j+9) =  STCmTt(j,i);
      stiffmatrix(i+3,j)   =  STCmTt(i,j);
      stiffmatrix(i+3,j+3) =  STCmTtSt(i,j) + TCbTt(i,j);
      stiffmatrix(i+3,j+6) = -STCmTt(i,j);
      stiffmatrix(i+3,j+9) =  STCmTtSt(j,i) - TCbTt(i,j);
      stiffmatrix(i+6,j)   = -TCmTt(i,j);
      stiffmatrix(i+6,j+3) = -STCmTt(j,i);
      stiffmatrix(i+6,j+6) =  TCmTt(i,j);
      stiffmatrix(i+6,j+9) = -STCmTt(j,i);
      stiffmatrix(i+9,j)   =  STCmTt(i,j);
      stiffmatrix(i+9,j+3) =  STCmTtSt(i,j)-TCbTt(i,j);
      stiffmatrix(i+9,j+6) = -STCmTt(i,j);
      stiffmatrix(i+9,j+9) =  STCmTtSt(i,j)+TCbTt(i,j);
    }
  }
  
  //Checking variation of strains with respect to displacements
  
  return;
} // DRT::ELEMENTS::Beam3::computestiffbasis

//computing spin matrix out of a rotation vector
inline void DRT::ELEMENTS::Beam3::computespin(BlitzMat3x3& spin, BlitzVec3 rotationangle, const double& spinscale)
{
  BLITZTINY::V_scale<3>(rotationangle,spinscale);
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

/*this function performs an update of the central triad as in principle given in Crisfield, Vol. 2, equation (17.65), but by means of a
 * quaterion product and then calculation of the equivalent rotation matrix according to eq. (16.70*/
inline void DRT::ELEMENTS::Beam3::updatetriad(BlitzVec3 deltabetaplusalpha, BlitzMat3x3& Tnew)
{
  //calculating angle theta by which triad is rotated according to Crisfield, Vol. 2, equation (17.64)
  BLITZTINY::V_scale<3>(deltabetaplusalpha,0.5);
  
  //absolute value of rotation angle theta
  double abs_theta = pow(deltabetaplusalpha(0)*deltabetaplusalpha(0) + deltabetaplusalpha(1)*deltabetaplusalpha(1) + deltabetaplusalpha(2)*deltabetaplusalpha(2) , 0.5);
  
  //computing quaterion for rotation by angle theta
  BlitzVec4 Qrot;
  if (abs_theta > 0)
  {
    Qrot(0) = deltabetaplusalpha(0) * sin(abs_theta / 2) / abs_theta;
    Qrot(1) = deltabetaplusalpha(1) * sin(abs_theta / 2) / abs_theta;
    Qrot(2) = deltabetaplusalpha(2) * sin(abs_theta / 2) / abs_theta;
    Qrot(3) = cos(abs_theta / 2);
  }
  else
  {
    BLITZTINY::PutScalar<4>(Qrot,0);
    Qrot(3) = 1;
  }
  
  //computing quaterion Qnew_ for new configuration of Qold_ for old configuration by means of a quaternion product
  Qnew_(0) = Qrot(3)*Qold_(0) + Qold_(3)*Qrot(0) + Qrot(1)*Qold_(2) - Qold_(1)*Qrot(2);
  Qnew_(1) = Qrot(3)*Qold_(1) + Qold_(3)*Qrot(1) + Qrot(2)*Qold_(0) - Qold_(2)*Qrot(0);
  Qnew_(2) = Qrot(3)*Qold_(2) + Qold_(3)*Qrot(2) + Qrot(0)*Qold_(1) - Qold_(0)*Qrot(1);
  Qnew_(3) = Qrot(3)*Qold_(3) - Qrot(2)*Qold_(2) - Qrot(1)*Qold_(1) - Qrot(0)*Qold_(0);
  
  //separate storage of vector part of Qnew_
  BlitzVec3 Qnewvec;
  for(int i = 0; i<3; i++)
  {
    Qnewvec(i) = Qnew_(i);
  }
  
  //computing the rotation matrix from Crisfield, Vol. 2, equation (17.70) with respect to Qnew_,
  //which is the new center triad Tnew
  computespin(Tnew, Qnewvec, 2*Qnew_(3));
  for(int i = 0; i<3; i++)
  {
    for(int j = 0; j<3; j++)
    {
      Tnew(i,j) += 2*Qnew_(i)*Qnew_(j);
      if(i == j)
        Tnew(i,j) += Qnew_(3)*Qnew_(3) - Qnew_(2)*Qnew_(2) - Qnew_(1)*Qnew_(1) - Qnew_(0)*Qnew_(0);
    }
  }
 
} //DRT::ELEMENTS::Beam3::updatetriad

//updating local curvature according to Crisfield, Vol. 2, pages 209 - 210; not: an exact update of the curvature is computed by
//means of equation (16.148) instead of an approximated one as given by equs. (17.72) and (17.73)
inline void DRT::ELEMENTS::Beam3::updatecurvature(const BlitzMat3x3& Tnew, BlitzVec3 deltabetaplusalpha,BlitzVec3 deltabetaminusalpha)
{
  //-------------------calculating omega-------------------------------------//
  
  //applying proper scaling for rotation angle
  BLITZTINY::V_scale<3>(deltabetaplusalpha,0.5);
 
  //absolute value of rotation vector theta
  double abs_theta = pow(deltabetaplusalpha(0)*deltabetaplusalpha(0) + deltabetaplusalpha(1)*deltabetaplusalpha(1) + deltabetaplusalpha(2)*deltabetaplusalpha(2) , 0.5);  
  
  BlitzVec3 omega = deltabetaplusalpha;
  BlitzVec3 omegaprime = deltabetaplusalpha;
  if (abs_theta > 0)
  {
    BLITZTINY::V_scale<3>(omega,2*tan(0.5*abs_theta) / abs_theta);
    BlitzMat3x3 Aux;
    for(int i = 0; i<3; i++)
    {
      for(int j = 0; j<3; j++)
      {
        Aux(i,j) = 0;
        Aux(i,j) -= (1 - abs_theta / sin(abs_theta) ) * deltabetaplusalpha(i)*deltabetaplusalpha(j) / pow(abs_theta,2);
        if(i==j)
          Aux(i,j) += 1;
          
       Aux(i,j) *= 2*tan(abs_theta / 2) / abs_theta;
      }
    }
    BLITZTINY::V_scale<3>(deltabetaminusalpha,1 / lrefe_);
    BLITZTINY::MV_product<3,3>(Aux,deltabetaminusalpha,omegaprime);
  }
  
  BlitzVec3 curvaux;
  curvaux(0) = 0.5*(omega(1)*omegaprime(2) - omega(2)*omegaprime(1)) ;
  curvaux(1) = 0.5*(omega(2)*omegaprime(0) - omega(0)*omegaprime(2)) ;
  curvaux(2) = 0.5*(omega(0)*omegaprime(1) - omega(1)*omegaprime(0)) ;
  
  curvaux += omegaprime;
  BLITZTINY::V_scale<3>(curvaux, 1/(1 + pow(tan(abs_theta/2),2) ));
  
  BLITZTINY::MtV_product<3,3>(Tnew,curvaux,curvnew_);
  curvnew_ += curvold_;
  
  return;
} //DRT::ELEMENTS::Beam3::updatecurvature

//computing stiffens matrix Ksigma1 according to Crisfield, Vol. 2, equation (17.83)
inline void DRT::ELEMENTS::Beam3::computeKsig1(Epetra_SerialDenseMatrix& Ksig1, const BlitzVec3& stressn, const BlitzVec3& stressm)
{
  Ksig1.Shape(12,12);
  BlitzMat3x3 Sn;
  BlitzMat3x3 Sm;
  computespin(Sn,stressn, 0.5);
  computespin(Sm,stressm, 0.5);

  for (int i=0; i<3; ++i)
  {
    for (int j=0; j<3; ++j)
    {
      Ksig1(i ,j+3) = Sn(i,j);
      Ksig1(i ,j+9) = Sn(i,j);
      Ksig1(i+6,j+3) = -Sn(i,j);
      Ksig1(i+6,j+9) = -Sn(i,j);

      Ksig1(i+3,j+3) = Sm(i,j);
      Ksig1(i+3,j+9) = Sm(i,j);
      Ksig1(i+9,j+3) = -Sm(i,j);
      Ksig1(i+9,j+9) = -Sm(i,j);
    }
  }
  return;
} //DRT::ELEMENTS::Beam3::computeKsig1

//computing stiffens matrix Ksigma1 according to Crisfield, Vol. 2, equation (17.87) and (17.88)
inline void DRT::ELEMENTS::Beam3::computeKsig2(Epetra_SerialDenseMatrix& Ksig2, const BlitzVec3& stressn, const BlitzVec3& x21)
{
  Ksig2.Shape(12,12);
  BlitzMat3x3 Sn;
  BlitzMat3x3 Sx21;
  BlitzMat3x3 Y;
  
  computespin(Sn,stressn, 0.5); 
  computespin(Sx21,x21, 0.5);
  BLITZTINY::MM_product<3,3,3>(Sx21,Sn,Y);

  for (int i=0; i<3; ++i)
  {
    for (int j=0; j<3; ++j)
    {
      Ksig2(i+3,j ) = -Sn(i,j);
      Ksig2(i+3,j+6) = Sn(i,j);
      Ksig2(i+9,j ) = -Sn(i,j);
      Ksig2(i+9,j+6) = Sn(i,j);

      Ksig2(i+3,j+3) = Y(i,j);
      Ksig2(i+3,j+9) = Y(i,j);
      Ksig2(i+9,j+3) = Y(i,j);
      Ksig2(i+9,j+9) = Y(i,j);
    }
  }
  return;
} // DRT::ELEMENTS::Beam3::computeKsig2


/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   cyron 01/08|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam3::b3_nlnstiffmass( ParameterList& params,
                                            vector<double>&           vel,
                                            vector<double>&           disp,
                                            Epetra_SerialDenseMatrix* stiffmatrix,
                                            Epetra_SerialDenseMatrix* massmatrix,
                                            Epetra_SerialDenseVector* force)
{
  //constitutive laws from Crisfield, Vol. 2, equation (17.76)
  BlitzVec3 Cm;
  BlitzVec3 Cb;
  
  //normal/shear strain and bending strain(curvature)
  BlitzVec3 epsilonn;
  BlitzVec3 epsilonm;
  
  //stress values n and m, Crisfield, Vol. 2, equation (17.78)
  BlitzVec3 stressn;
  BlitzVec3 stressm;
  
  //current position of nodal degrees of freedom
  BlitzMat6x2 xcurr;
  
  //difference between coordinates of both nodes in current configuration, x21' Crisfield  Vol. 2 equ. (17.66a) and (17.72)
  BlitzVec3 x21;
  
  //nonlinear parts of stiffness matrix, Crisfiel Vol. 2, equation (17.83) and (17.87)
  Epetra_SerialDenseMatrix Ksig1;
  Epetra_SerialDenseMatrix Ksig2;
  
  //auxiliary variables
  BlitzVec3 deltabetaplusalpha;
  BlitzVec3 deltabetaminusalpha;
  //midpoint triad, Crisfiel Vol. 2, equation (17.73)
  BlitzMat3x3 Tnew;

  //nodal coordinates in current position
  for (int k=0; k<2; ++k) //looping over number of nodes
  {
    for (int j=0; j<3; ++j) 
    {
      xcurr(j,k)   = Nodes()[k]->X()[j] + disp[k*6+j]; //translational DOF
      xcurr(j+3,k) = disp[k*6+j+3]; //rotational DOF
    }
  }

  //first of all "new" variables have to be adopted to dispalcement passed in from BACI driver

  //difference between coordinates of both nodes, x21' Crisfield  Vol. 2 equ. (17.66a) and (17.72)
  for (int j=0; j<3; ++j)
  {
    x21(j) = xcurr(j,1) - xcurr(j,0);
    betaplusalphanew_(j) = xcurr(j+3,1) + xcurr(j+3,0);
    betaminusalphanew_(j) = xcurr(j+3,1) - xcurr(j+3,0);
  }

  deltabetaplusalpha  = betaplusalphanew_;
  deltabetaplusalpha -= betaplusalphaold_;
  
  deltabetaminusalpha  = betaminusalphanew_;
  deltabetaminusalpha -= betaminusalphaold_;
  
  //calculating current central triad like in Crisfield, Vol. 2, equation (17.65), but by a quaternion product
  updatetriad(deltabetaplusalpha,Tnew);
  
  //updating local curvature
  updatecurvature(Tnew, deltabetaplusalpha,deltabetaminusalpha);
   
  //computing current axial and shear strain epsilon, Crisfield, Vol. 2, equation (17.67)
  BLITZTINY::MtV_product<3,3>(Tnew,x21,epsilonn);
  BLITZTINY::V_scale<3>(epsilonn,1/lrefe_);
  epsilonn(0) -=  1;

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
  epsilonn(0) *= ym*crosssec_;
  epsilonn(1) *= sm*crosssecshear_;
  epsilonn(2) *= sm*crosssecshear_;
  BLITZTINY::MV_product<3,3>(Tnew,epsilonn,stressn);  
   
  //computing spin matrix S(x21)/2 according to Crisfield, Vol. 2, equation (17.74)
  BlitzMat3x3 spinx21;
  computespin(spinx21,x21,0.5); 
  
  //lamda is the derivative of current velocity with respect to current displacement
  double lamda = params.get<double>("gamma",0.581) / (params.get<double>("delta time",0.01)*params.get<double>("beta",0.292));
  
  
  //turning bending strain epsilonm into bending stress stressm
  epsilonm = curvnew_;
  epsilonm(0) *= sm*Irr_;
  epsilonm(1) *= ym*Iyy_;
  epsilonm(2) *= ym*Izz_;
  BLITZTINY::MV_product<3,3>(Tnew,epsilonm,stressm); 
  
  //matrix Theta can filter out of a vector the component parallel to beam axis (used for artificial torsional damping):
  BlitzMat3x3 Theta;
  for(int i = 0; i<3; i++)
  {
    for(int j = 0; j<3; j++)
    {
      Theta(i,j) = Tnew(i,0)*Tnew(0,j);
    }
  }
  
  //computing global internal forces, Crisfield Vol. 2, equation (17.79)
  //note: X = [-I 0; -S -I; I 0; -S I] with -S = T^t; and S = S(x21)/2;
  if (force != NULL)
  {
    (*force).Size(12);
    for (int i=0; i<3; ++i)
    {
      (*force)(i)   -= stressn(i);
      (*force)(i+3) -= stressm(i);
      (*force)(i+6) += stressn(i);
      (*force)(i+9) += stressm(i);
      
      for (int j=0; j<3; ++j)
      {      
        (*force)(i+3) -= stressn(j)*spinx21(i,j);      
        (*force)(i+9) -= stressn(j)*spinx21(i,j);    
      }
    }
       
    //for problems of statistical mechanics viscous damping is incalculated
    if(kT_ > 0)
    {
      if (stochasticorder_ == 0)
      {
        //adding internal forces due to viscous damping (by background fluid of thermal bath) with zeroth order stochastic interpolation
        (*force)[0] += (zeta_/2)*vel[0];
        (*force)[1] += (zeta_/2)*vel[1];
        (*force)[2] += (zeta_/2)*vel[2];
        (*force)[6] += (zeta_/2)*vel[6];
        (*force)[7] += (zeta_/2)*vel[7];
        (*force)[8] += (zeta_/2)*vel[8];
      }
      else if (stochasticorder_ == 1)
      {
        //adding internal forces due to viscous damping (by background fluid of thermal bath) with first order stochastic interpolation
        (*force)[0] += zeta_*(vel[0]/3 + vel[6]/6);
        (*force)[1] += zeta_*(vel[1]/3 + vel[7]/6);
        (*force)[2] += zeta_*(vel[2]/3 + vel[8]/6);
        (*force)[6] += zeta_*(vel[6]/3 + vel[0]/6);
        (*force)[7] += zeta_*(vel[7]/3 + vel[1]/6);
        (*force)[8] += zeta_*(vel[8]/3 + vel[2]/6);
      }          
    }
    for(int i = 0; i<3;i++)
    {
      (*force)[i+3] += (zeta_/2)*vel[i+3]*0.01;
      (*force)[i+9] += (zeta_/2)*vel[i+9]*0.01;
    }
    
    /*    
    //adding artificial torsional damping in order to stabilize numerically free fluctuations
    double torsdamp = 0.01;
    for(int i = 0; i<3; i++)
    {
      for(int j = 0; j<3; j++)
      {
        //node 1
        (*force)[3+i] += Theta(i,j)*vel[3+j]*(zeta_/2)*torsdamp;
        //node 2
        (*force)[9+i] += Theta(i,j)*vel[9+j]*(zeta_/2)*torsdamp;
      }
    }
    */
    
    
    
  }
  

  //computing linear stiffness matrix
  if (stiffmatrix != NULL)
  {
    (*stiffmatrix).Shape(12,12);    
    
    //setting constitutive parameters , Crisfield, Vol. 2, equation (17.76)
    Cm(0) = ym*crosssec_/lrefe_;
    Cm(1) = sm*crosssecshear_/lrefe_;
    Cm(2) = sm*crosssecshear_/lrefe_;
    Cb(0) = sm*Irr_/lrefe_;
    Cb(1) = ym*Iyy_/lrefe_;
    Cb(2) = ym*Izz_/lrefe_;
    
    //setting up basis of stiffness matrix according to Crisfield, Vol. 2, equation (17.81)       
    computestiffbasis(Tnew,Cm,Cb,spinx21,(*stiffmatrix));
       
    //adding nonlinear (stress dependent) parts to tangent stiffness matrix, Crisfield, Vol. 2 equs. (17.83), (17.87), (17.89)
    computeKsig1(Ksig1,stressn,stressm);
    computeKsig2(Ksig2,stressn,x21);
    (*stiffmatrix) += Ksig1;
    (*stiffmatrix) += Ksig2;
    
    /*note: the following operations are carried out in order to integrate a semiexplicit viscous damping into the stiffness as necessary
     * for simulations of statistical mechanics; consequently there is no counterpart of the following lines in the textbook of Crisfield*/
    
    //adding diagonal entries for viscous damping "stiffness" (by background fluid of thermal bath) for problems of statistical mechanics
    if(kT_ > 0)
    {
      if (stochasticorder_ == 0)
      {
        //adding internal forces due to viscous damping (by background fluid of thermal bath) with zeroth order stochastic interpolation
        (*stiffmatrix)(0,0) += (zeta_/2)*lamda;
        (*stiffmatrix)(1,1) += (zeta_/2)*lamda;
        (*stiffmatrix)(2,2) += (zeta_/2)*lamda;
        (*stiffmatrix)(6,6) += (zeta_/2)*lamda;
        (*stiffmatrix)(7,7) += (zeta_/2)*lamda;
        (*stiffmatrix)(8,8) += (zeta_/2)*lamda;
      }
      else if (stochasticorder_ == 1)
      {
        //adding internal forces due to viscous damping (by background fluid of thermal bath) with first order stochastic interpolation
        (*stiffmatrix)(0,0) += (zeta_/3)*lamda;
        (*stiffmatrix)(1,1) += (zeta_/3)*lamda;
        (*stiffmatrix)(2,2) += (zeta_/3)*lamda;
        (*stiffmatrix)(6,6) += (zeta_/3)*lamda;
        (*stiffmatrix)(7,7) += (zeta_/3)*lamda;
        (*stiffmatrix)(8,8) += (zeta_/3)*lamda;
        
        (*stiffmatrix)(0,6) += (zeta_/6)*lamda;
        (*stiffmatrix)(6,0) += (zeta_/6)*lamda;
        (*stiffmatrix)(1,7) += (zeta_/6)*lamda;
        (*stiffmatrix)(7,1) += (zeta_/6)*lamda;
        (*stiffmatrix)(2,8) += (zeta_/6)*lamda;
        (*stiffmatrix)(8,2) += (zeta_/6)*lamda;
      }
    }
    
    for(int i = 0; i<3;i++)
    {
      (*stiffmatrix)(i+3,i+3) += (zeta_/2)*lamda*0.01;
      (*stiffmatrix)(i+9,i+9) += (zeta_/2)*lamda*0.01;
    }
    

    /*
    //adding artifical viscosity stiffness with respect to torsional displacement
    for(int i = 0; i<2;i++)
    {
      double torsdamp = 0.01;
      BlitzVec3 aux1;
      BlitzVec3 aux2;
      BlitzMat3x3 Aux1;
      BlitzMat3x3 Aux2;
      for(int j = 0; j<3; j++)
        aux1(j) = (zeta_/2)*vel[j + i*6];
      BLITZTINY::MV_product<3,3>(Theta,aux1,aux2); 
      computespin(Aux1,aux2,-0.5); 
      for(int j= 0; j<3; j++)
      {
        for(int k=0;k<3;k++)
        {
          (*stiffmatrix)(i*6 +j ,3+k ) += Aux1(j,k)*torsdamp;
          (*stiffmatrix)(i*6 +j ,9+k ) += Aux1(j,k)*torsdamp;
        }
      }
      computespin(Aux1,aux1,0.5);
      BLITZTINY::MM_product<3,3,3>(Theta,Aux1,Aux2);
      for(int j= 0; j<3; j++)
      {
        for(int k=0;k<3;k++)
        {
          (*stiffmatrix)(i*6 +j ,3+k ) += Aux2(j,k)*torsdamp;
          (*stiffmatrix)(i*6 +j ,9+k ) += Aux2(j,k)*torsdamp;
        }
      }
      
      for(int j= 0; j<3; j++)
      {
        for(int k=0;k<3;k++)
        {
          (*stiffmatrix)(i*6 +j ,i*6 + k ) += Theta(j,k)*(zeta_/2)*lamda*torsdamp;
        }
      }     
     }
     */
    
    
     
  
   
  }
  
  /*calculating mass matrix; this beam3 element includes only a lumped mass matrix where for torsion and
   * bending the same moments of inertia are assumed; for slender beams the influence of rotational moments
   * of inertia is dilute so that often rotational inertia is just set to zero; since within this code at one
   * point an LU-decomposition of the mass matrix is carried out this was avoided in the here implemented beam3
   * element and instead the above described simplification was assumed; Iyy_ = Izz_ was assumed for similar 
   * reasons */
  if (massmatrix != NULL)
  {
    (*massmatrix).Shape(12,12);
    for (int i=0; i<3; ++i)
    {
      (*massmatrix)(i  ,i  ) = 0.5*density*lrefe_*crosssec_;
      (*massmatrix)(i+3,i+3) = 0.5*density*lrefe_*Iyy_;
      (*massmatrix)(i+6,i+6) = 0.5*density*lrefe_*crosssec_;
      (*massmatrix)(i+9,i+9) = 0.5*density*lrefe_*Iyy_;
    }   
  }
  
  return;
} // DRT::ELEMENTS::Beam3::b3_nlnstiffmass


// lump mass matrix
void DRT::ELEMENTS::Beam3::b3_lumpmass(Epetra_SerialDenseMatrix* emass)
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
