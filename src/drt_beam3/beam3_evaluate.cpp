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
  else if (action=="calc_brownian_damp")        act = Beam3::calc_brownian_damp;
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
    case Beam3::calc_brownian_damp:
    {
      /*in order to understand the way how in the following parallel computing is handled one has to
       * be aware of what usually happens: each processor evaluates forces on each of its elements no
       * matter whether it's a real element or only a ghost element; later on in the assembly of the
       * evaluation method each processor adds the thereby gained DOF forces only for the DOF of
       * which it is the row map owner; if one used this way for random forces the following problem
       * would arise: if two processors evaluate both the same element (one time as a real element, one time
       * as a ghost element) and assemble later only the random forces of their own DOF the assembly for
       * different DOF of the same element might take place by means of random forces evaluated during different
       * element calls (one time the real element was called, one time only the ghost element); as a consequence
       * prescribing any correlation function between random forces of different DOF is in general
       * impossible since the random forces of two different DOF may have been evaluated during different
       * element calls; the only solution to this problem is obviously to allow element evaluation
       * to one processor only (the row map owner processor) and to allow this processor to assemble the
       * thereby gained random forces later on also to DOF of which it is not the row map owner; so first
       * we check in this method whether the current processor is the owner of the element (if not
       * evaluation is cancelled) and finally this processor assembles the evaluated forces for degrees of
       * freedom right no matter whether its the row map owner of these DOF (outside the element routine
       * in the evaluation method of the discretization this would be not possible by default*/


      //test whether current processor is row map owner of the element
      if(this->Owner() != discretization.Comm().MyPID()) return 0;

      // get element displacements (for use in shear flow fields)
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get state vector 'displacement'");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

      //evaluation of statistical forces or movement (both stored in variable "brownian")
      Epetra_SerialDenseVector brownian(lm.size());
      EvaluateBrownian(params,mydisp,brownian,elemat1);

      /*all the above evaluated forces or movements are used in assembly of the in general column map vector no matter if the current processor is
       * the row map owner of the affected DOF*/
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
      //correctly or not by means of a numerically approximated stiffness matrix; remark: due to involved numerics for
      //this element it's normal that numerical derivation entails relative errors ~1e-2,-3,... in the diagonal elements
      //of those blocks of the stiffness matrix which connect angular displacements with normal and shear strains
      if(Id() == 0) //limiting the following tests to certain element numbers
      {
        //variable to store numerically approximated stiffness matrix
        Epetra_SerialDenseMatrix stiff_approx;
        stiff_approx.Shape(12,12);

        //relative error of numerically approximated stiffness matrix
        Epetra_SerialDenseMatrix stiff_relerr;
        stiff_relerr.Shape(12,12);

        //characteristic length for numerical approximation of stiffness
        double h_rel = 1e-7;

        //flag indicating whether approximation lead to significant relative error
        int outputflag = 0;

        //calculating strains in new configuration
        for(int k=0; k<2; k++)
        {
          for(int i=0; i<6; i++)
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
            disp_aux[i + 6*k] += h_rel;
             vel_aux[i + 6*k] += h_rel * params.get<double>("gamma",0.581) / ( params.get<double>("delta time",0.01)*params.get<double>("beta",0.292) );

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

  //jacobian determinant
  double det = lrefe_/2;

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
    for (int node=0; node<NumNode(); ++node)
    for (int dof=0; dof<numdf; ++dof)
    elevec1[node*numdf+dof] += funct[node] *ar[dof];

  } // for (int ip=0; ip<intpoints.nquad; ++ip)

  return 0;
}


/*-----------------------------------------------------------------------------------------------------------*
 | Evaluate Statistical forces and viscous damping (public)                                       cyron 09/08|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Beam3::EvaluateBrownian(ParameterList& params,
                                           vector<double> mydisp,
                                           Epetra_SerialDenseVector& elevec1,
                                           Epetra_SerialDenseMatrix& elemat1)
{

  // thermal energy responsible for statistical forces
  double kT = params.get<double>("KT",0.0);

  // frictional coefficient per unit length (approximated by the one of an infinitely long staff)
  double zeta = 4 * PI * lrefe_ * params.get<double>("ETA",0.0);

  // polynomial order for interpolation of stochastic line load (zero corresponds to bead spring model)
  int stochasticorder = params.get<int>("STOCH_ORDER",0);

  //stochastic field of line load is interpolated by zeroth order polynomial functions
  if (stochasticorder == 0)
  {
    //adding internal forces due to viscous damping (by background fluid of thermal bath) with zeroth order stochastic interpolation
    elemat1(0,0) += zeta/2.0;
    elemat1(1,1) += zeta/2.0;
    elemat1(2,2) += zeta/2.0;
    elemat1(6,6) += zeta/2.0;
    elemat1(7,7) += zeta/2.0;
    elemat1(8,8) += zeta/2.0;


    //calculating standard deviation of statistical forces or Brownian steps according to fluctuation dissipation theorem   
    double stand_dev_trans = 0;
    //if FORCE_OR_DISP == 0 a statistical force is evaluated
    if(params.get<int>("FORCE_OR_DISP",0) == 0)
      stand_dev_trans = pow(2 * kT * (zeta/2) / params.get<double>("delta time",0.01),0.5);
    //if FORCE_OR_DISP == 1 a statistical displacement is evaluated
    else
      stand_dev_trans = pow(2 * kT * params.get<double>("delta time",0.01) / (zeta/2),0.5);
    


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
  else if (stochasticorder == 1)
  {
    //adding internal forces due to viscous damping (by background fluid of thermal bath) with first order stochastic interpolation
    elemat1(0,0) += zeta/3.0;
    elemat1(1,1) += zeta/3.0;
    elemat1(2,2) += zeta/3.0;
    elemat1(6,6) += zeta/3.0;
    elemat1(7,7) += zeta/3.0;
    elemat1(8,8) += zeta/3.0;

    elemat1(0,6) += zeta/6.0;
    elemat1(6,0) += zeta/6.0;
    elemat1(1,7) += zeta/6.0;
    elemat1(7,1) += zeta/6.0;
    elemat1(2,8) += zeta/6.0;
    elemat1(8,2) += zeta/6.0;
    
    //calculating standard deviation of statistical forces or Brownian steps according to fluctuation dissipation theorem   
    double stand_dev_trans = 0;
    //if FORCE_OR_DISP == 0 a statistical force is evaluated
    if(params.get<int>("FORCE_OR_DISP",0) == 0)
      stand_dev_trans = pow(2 * kT * (zeta/6) / params.get<double>("delta time",0.01),0.5);
    //if FORCE_OR_DISP == 1 a statistical displacement is evaluated
    else
      stand_dev_trans = pow(2 * kT * params.get<double>("delta time",0.01) / (zeta/6),0.5);

    
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
    double stat1 = normalGen.random();
    double stat2 = normalGen.random();
    double stat3 = normalGen.random();
    elevec1[0] += stat1;
    elevec1[1] += stat2;
    elevec1[2] += stat3;
    elevec1[6] += stat1;
    elevec1[7] += stat2;
    elevec1[8] += stat3;
  }
  
  return 0;
} //DRT::ELEMENTS::Beam3::EvaluateBrownian



/*-----------------------------------------------------------------------------------------------------------*
 | Evaluate PTC damping (public)                                                                  cyron 10/08|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Beam3::EvaluatePTC(ParameterList& params,
                                      Epetra_SerialDenseMatrix& elemat1)
{
  LINALG::Matrix<3,1> newangle;
  quaterniontoangle(Qnew_, newangle);
  LINALG::Matrix<3,3> Hinverse = Hinv(newangle);
  double dti = params.get<double>("dti",0.0);

  Hinverse.Scale(dti);

  for(int i= 0; i<3; i++)
  {
    for(int j=0;j<3;j++)
    {
      
      //translational damping
      elemat1(  i,   j) += dti*0.5;
      elemat1(6+i, 6+j) += dti*0.5;
      elemat1(6+i,   j) += dti*0.5;
      elemat1(  i, 6+j) += dti*0.5;
      
      
      //rotational damping
      elemat1(3+i, 3+j) += Hinverse(i,j);
      elemat1(9+i, 9+j) += Hinverse(i,j);
      elemat1(9+i, 3+j) += Hinverse(i,j);
      elemat1(3+i, 9+j) += Hinverse(i,j);
    }
  }

  return 0;
} //DRT::ELEMENTS::Beam3::EvaluatePTC

/*-----------------------------------------------------------------------------------------------------------*
 | auxiliary functions for dealing with large rotations and nonlinear stiffness                    cyron 04/08|
 *----------------------------------------------------------------------------------------------------------*/
//computing basis of stiffness matrix of Crisfield, Vol. 2, equation (17.81)
inline void DRT::ELEMENTS::Beam3::computestiffbasis(const LINALG::Matrix<3,3>& Tnew, const LINALG::Matrix<3,1>& Cm, const LINALG::Matrix<3,1>& Cb, const LINALG::Matrix<3,3>& spinx21, Epetra_SerialDenseMatrix& stiffmatrix)
{
  //calculating the first matrix of (17.81) directly involves multiplications of large matrices (e.g. with the 12x6-matrix X)
  //application of the definitions in (17.74) allows blockwise evaluation with multiplication and addition of 3x3-matrices only
  //in the follwoing all the blocks are treated separately; their name is related directly to their content:
  //e.g. TCmTt is the product of the 3 matrices T * C_m * T^t (with T and C_m according to (17.74) and (17.76)
  //for the blockwise calculation on which the following steps are based on the relation S^t = -S for spin matrices was applied

  LINALG::Matrix<3,3> TCmTt;
  LINALG::Matrix<3,3> TCbTt;
  LINALG::Matrix<3,3> STCmTt;
  LINALG::Matrix<3,3> STCmTtSt;

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
      stiffmatrix(i+6,j)   = -TCmTt(i,j);
      stiffmatrix(i+6,j+6) =  TCmTt(i,j);
      stiffmatrix(i  ,j+6) = -TCmTt(i,j);

      stiffmatrix(i  ,j+3) =  STCmTt(j,i);
      stiffmatrix(i  ,j+9) =  STCmTt(j,i);
      stiffmatrix(i+6,j+3) = -STCmTt(j,i);
      stiffmatrix(i+6,j+9) = -STCmTt(j,i);

      stiffmatrix(i+3,j+6) = -STCmTt(i,j);
      stiffmatrix(i+3,j)   =  STCmTt(i,j);
      stiffmatrix(i+9,j)   =  STCmTt(i,j);
      stiffmatrix(i+9,j+6) = -STCmTt(i,j);

      stiffmatrix(i+3,j+3) =  STCmTtSt(i,j) + TCbTt(i,j);
      stiffmatrix(i+3,j+9) =  STCmTtSt(j,i) - TCbTt(i,j);
      stiffmatrix(i+9,j+3) =  STCmTtSt(i,j) - TCbTt(i,j);
      stiffmatrix(i+9,j+9) =  STCmTtSt(i,j) + TCbTt(i,j);
    }
  }


  return;
} // DRT::ELEMENTS::Beam3::computestiffbasis

/*----------------------------------------------------------------------*
 |computes from a quaternion q the related angle theta (public)cyron10/08|
 *----------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Beam3::quaterniontoangle(const LINALG::Matrix<4,1>& q, LINALG::Matrix<3,1>& theta)
{  
  /*the following funciton computes from a quaternion q an angle theta within [-PI; PI]; such an interval is
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

//computing spin matrix out of a rotation vector
inline void DRT::ELEMENTS::Beam3::computespin(LINALG::Matrix<3,3>& spin, LINALG::Matrix<3,1> rotationangle, const double& spinscale)
{
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

//computing a rotation matrix R from a quaternion q, cf. Crisfield, Vol. 2, equation (17.70)
inline void DRT::ELEMENTS::Beam3::quaterniontotriad(const LINALG::Matrix<4,1>& q, LINALG::Matrix<3,3>& R)
{
  //separate storage of vector part of q
  LINALG::Matrix<3,1> qvec;
  for(int i = 0; i<3; i++)
    qvec(i) = q(i);

  //setting R to third summand of equation (17.70)
  computespin(R, qvec, 2*q(3));

  //adding second summand of equation (17.70)
  for(int i = 0; i<3; i++)
  {
    for(int j = 0; j<3; j++)
    {
      R(i,j) += 2*q(i)*q(j);
    }
  }

  //correting diagonal entries according to first summand of equation (17.70)
  R(0,0) = 1 - 2*(q(1)*q(1) + q(2)*q(2));
  R(1,1) = 1 - 2*(q(0)*q(0) + q(2)*q(2));
  R(2,2) = 1 - 2*(q(0)*q(0) + q(1)*q(1));

  return;
} // DRT::ELEMENTS::Beam3::quaterniontotriad



/*computing a quaternion q from a rotation matrix R; all operations are performed according to
* Crisfield, Vol. 2, section 16.10 and the there described Spurrier's algorithm*/
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


/*matrix H^(-1) which turns non-additive spin variables into additive ones according to Crisfield, Vol. 2, equation (16.93)*/
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
    {
      for(int j=0; j<3; j++)
      {
        result(i,j) += theta(i) * theta(j) * (1 - theta_abs/(2*tan(theta_abs/2)) )/pow(theta_abs,2);
      }
    }
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



/*this function performs an update of the central triad as in principle given in Crisfield, Vol. 2, equation (17.65), but by means of a
 * quaterion product and then calculation of the equivalent rotation matrix according to eq. (16.70*/
inline void DRT::ELEMENTS::Beam3::updatetriad(LINALG::Matrix<3,1> deltabetaplusalpha, LINALG::Matrix<3,3>& Tnew)
{
  //calculating angle theta by which triad is rotated according to Crisfield, Vol. 2, equation (17.64)
  deltabetaplusalpha.Scale(0.5);

  //absolute value of rotation angle theta
  double abs_theta = pow(deltabetaplusalpha(0)*deltabetaplusalpha(0) + deltabetaplusalpha(1)*deltabetaplusalpha(1) + deltabetaplusalpha(2)*deltabetaplusalpha(2) , 0.5);

  //computing quaterion for rotation by angle theta
  LINALG::Matrix<4,1> Qrot;
  if (abs_theta > 0)
  {
    Qrot(0) = deltabetaplusalpha(0) * sin(abs_theta / 2) / abs_theta;
    Qrot(1) = deltabetaplusalpha(1) * sin(abs_theta / 2) / abs_theta;
    Qrot(2) = deltabetaplusalpha(2) * sin(abs_theta / 2) / abs_theta;
    Qrot(3) = cos(abs_theta / 2);
  }
  else
  {
    Qrot.PutScalar(0.0);
    Qrot(3) = 1;
  }

  //computing quaterion Qnew_ for new configuration of Qold_ for old configuration by means of a quaternion product
  Qnew_(0) = Qrot(3)*Qold_(0) + Qold_(3)*Qrot(0) + Qrot(1)*Qold_(2) - Qold_(1)*Qrot(2);
  Qnew_(1) = Qrot(3)*Qold_(1) + Qold_(3)*Qrot(1) + Qrot(2)*Qold_(0) - Qold_(2)*Qrot(0);
  Qnew_(2) = Qrot(3)*Qold_(2) + Qold_(3)*Qrot(2) + Qrot(0)*Qold_(1) - Qold_(0)*Qrot(1);
  Qnew_(3) = Qrot(3)*Qold_(3) - Qrot(2)*Qold_(2) - Qrot(1)*Qold_(1) - Qrot(0)*Qold_(0);
  
  //normalizing quaternion in order to make sure that it keeps unit absolute values throught time stepping
  double abs = pow(Qnew_(0)*Qnew_(0) + Qnew_(1)*Qnew_(1) + Qnew_(2)*Qnew_(2) + Qnew_(3)*Qnew_(3),0.5);
  for(int i = 0; i<4; i++)
    Qnew_(i) = Qnew_(i) / abs;

  quaterniontotriad(Qnew_,Tnew);
  

} //DRT::ELEMENTS::Beam3::updatetriad

//updating local curvature according to Crisfield, Vol. 2, pages 209 - 210; not: an exact update of the curvature is computed by
//means of equation (16.148) instead of an approximated one as given by equs. (17.72) and (17.73)
inline void DRT::ELEMENTS::Beam3::updatecurvature(const LINALG::Matrix<3,3>& Tnew, LINALG::Matrix<3,1> deltabetaplusalpha,LINALG::Matrix<3,1> deltabetaminusalpha)
{
  //-------------------calculating omega-------------------------------------//

  //applying proper scaling for rotation angle
  deltabetaplusalpha.Scale(0.5);

  //absolute value of rotation vector theta
  double abs_theta = pow(deltabetaplusalpha(0)*deltabetaplusalpha(0) + deltabetaplusalpha(1)*deltabetaplusalpha(1) + deltabetaplusalpha(2)*deltabetaplusalpha(2) , 0.5);

  LINALG::Matrix<3,1> omega = deltabetaplusalpha;
  LINALG::Matrix<3,1> omegaprime = deltabetaplusalpha;
  if (abs_theta > 0)
  {
    omega.Scale(2*tan(0.5*abs_theta) / abs_theta);
    LINALG::Matrix<3,3> Aux;
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
    deltabetaminusalpha.Scale(1 / lrefe_);
    omegaprime.Multiply(Aux,deltabetaminusalpha);
  }

  LINALG::Matrix<3,1> curvaux;
  curvaux(0) = 0.5*(omega(1)*omegaprime(2) - omega(2)*omegaprime(1)) ;
  curvaux(1) = 0.5*(omega(2)*omegaprime(0) - omega(0)*omegaprime(2)) ;
  curvaux(2) = 0.5*(omega(0)*omegaprime(1) - omega(1)*omegaprime(0)) ;

  curvaux += omegaprime;
  curvaux.Scale( 1/(1 + pow(tan(abs_theta/2),2) ));

  curvnew_.MultiplyTN(Tnew,curvaux);
  curvnew_ += curvold_;

  return;
} //DRT::ELEMENTS::Beam3::updatecurvature


//updating local curvature according approximately to Crisfield, Vol. 2, eqs. (17.72) and (17.73)
inline void DRT::ELEMENTS::Beam3::approxupdatecurvature(const LINALG::Matrix<3,3>& Tnew, LINALG::Matrix<3,1> deltabetaplusalpha,LINALG::Matrix<3,1> deltabetaminusalpha)
{
  //old triad
  LINALG::Matrix<3,3> Told;
  quaterniontotriad(Qold_,Told);
  
  //compute spin matrix from eq. (17.73)
  LINALG::Matrix<3,3> spin;
  computespin(spin, deltabetaplusalpha, 0.5);
  
  //turning spin matrix to left right hand side matrix of eq. (17.73)
  for(int i = 0; i<3; i++)
    spin(i,i) += 1;

  //complete right hand side matrix of eq. (17.73)
  //mid point triad
  LINALG::Matrix<3,3> Tmid;
  Tmid.Multiply(spin,Told);
  
  //eq. (17.72)
  curvnew_.MultiplyTN(Tmid,deltabetaminusalpha);
  curvnew_.Scale(1/lrefe_);

  curvnew_ += curvold_;

  return;
} //DRT::ELEMENTS::Beam3::approxupdatecurvature


//computing stiffens matrix Ksigma1 according to Crisfield, Vol. 2, equation (17.83)
inline void DRT::ELEMENTS::Beam3::computeKsig1(Epetra_SerialDenseMatrix& Ksig1, const LINALG::Matrix<3,1>& stressn, const LINALG::Matrix<3,1>& stressm)
{
  Ksig1.Shape(12,12);
  LINALG::Matrix<3,3> Sn;
  LINALG::Matrix<3,3> Sm;
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
inline void DRT::ELEMENTS::Beam3::computeKsig2(Epetra_SerialDenseMatrix& Ksig2, const LINALG::Matrix<3,1>& stressn, const LINALG::Matrix<3,1>& x21)
{
  Ksig2.Shape(12,12);
  LINALG::Matrix<3,3> Sn;
  LINALG::Matrix<3,3> Sx21;
  LINALG::Matrix<3,3> Y;

  computespin(Sn,stressn, 0.5);
  computespin(Sx21,x21, 0.5);
  Y.Multiply(Sx21,Sn);

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
  LINALG::Matrix<3,1> Cm;
  LINALG::Matrix<3,1> Cb;

  //normal/shear strain and bending strain(curvature)
  LINALG::Matrix<3,1> epsilonn;
  LINALG::Matrix<3,1> epsilonm;

  //stress values n and m, Crisfield, Vol. 2, equation (17.78)
  LINALG::Matrix<3,1> stressn;
  LINALG::Matrix<3,1> stressm;

  //difference between coordinates of both nodes in current configuration, x21' Crisfield  Vol. 2 equ. (17.66a) and (17.72)
  LINALG::Matrix<3,1> x21;

  //nonlinear parts of stiffness matrix, Crisfiel Vol. 2, equation (17.83) and (17.87)
  Epetra_SerialDenseMatrix Ksig1;
  Epetra_SerialDenseMatrix Ksig2;

  //auxiliary variables
  LINALG::Matrix<3,1> deltabetaplusalpha;
  LINALG::Matrix<3,1> deltabetaminusalpha;
  
  //midpoint triad, Crisfiel Vol. 2, equation (17.73)
  LINALG::Matrix<3,3> Tnew;

  //first of all "new" variables have to be adopted to dispalcement passed in from BACI driver

  //difference between coordinates of both nodes, x21' Crisfield  Vol. 2 equ. (17.66a) and (17.72)
  for (int j=0; j<3; ++j)
  {
    x21(j)                = (X_(3+j) - X_(j) ) + ( disp[6+j]  - disp[j] );
    betaplusalphanew_(j)  = disp[9+j] + disp[3+j];
    betaminusalphanew_(j) = disp[9+j] - disp[3+j];
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
  epsilonn.MultiplyTN(Tnew,x21);
  epsilonn.Scale(1/lrefe_);
  epsilonn(0) -=  1;


  /* read material parameters using structure _MATERIAL which is defined by inclusion of      /
   / "../drt_lib/drt_timecurve.H"; note: material parameters have to be read in the evaluation /
   / function instead of e.g. Beam3_input.cpp or within the Beam3Register class since it is not/
   / sure that structure _MATERIAL is declared within those scopes properly whereas it is within/
   / the evaluation functions */

  // get the material law
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

  //computing spin matrix S(x21)/2 according to Crisfield, Vol. 2, equation (17.74)
  LINALG::Matrix<3,3> spinx21;
  computespin(spinx21,x21,0.5);

  //stress values n and m, Crisfield, Vol. 2, equation (17.76) and (17.78)
  epsilonn(0) *= ym*crosssec_;
  epsilonn(1) *= sm*crosssecshear_;
  epsilonn(2) *= sm*crosssecshear_;

  stressn.Multiply(Tnew,epsilonn);
        

  //turning bending strain epsilonm into bending stress stressm
  epsilonm = curvnew_;
  epsilonm(0) *= sm*Irr_;
  epsilonm(1) *= ym*Iyy_;
  epsilonm(2) *= ym*Izz_;
  stressm.Multiply(Tnew,epsilonm);
  
  /*
  std::cout<<"\nstressn"<<stressn;
  std::cout<<"\nstressm"<<stressm;
*/

  //computing global internal forces, Crisfield Vol. 2, equation (17.79)
  //note: X = [-I 0; -S -I; I 0; -S I] with -S = T^t; and S = S(x21)/2;
  if (force != NULL)
  {
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
    

    
    //artificial isotropic rotational damping
    {
      double torsdamp = 0.000005; //0.005 is obviously sufficient for free fluctionations with 10 elements
      LINALG::Matrix<3,1> newangle;
      LINALG::Matrix<3,1> convangle;
      quaterniontoangle(Qnew_, newangle);
      quaterniontoangle(Qconv_, convangle);
      
      LINALG::Matrix<3,1> omega = newangle;
      omega -= convangle;
      omega.Scale( 1.0 / params.get<double>("delta time",0.01) );
      for(int i = 0; i<3; i++)
      {
        //node 1     
        (*force)[3+i] += omega(i)*torsdamp*zeta_;
        //node 2
        (*force)[9+i] += omega(i)*torsdamp*zeta_;      
      }
    }
    
    
    /*
    //artificial isotropic curvature damping
    {
      double torsdamp = 0.05; //0.005 is obviously sufficient for free fluctionations with 10 elements
      
      LINALG::Matrix<3,1> omega = curvnew_;
      omega -= curvconv_;
      omega.Scale( 1.0 / params.get<double>("delta time",0.01) );
      for(int i = 0; i<3; i++)
      {
        //node 1     
        (*force)[3+i] += omega(i)*torsdamp*zeta_;
        //node 2
        (*force)[9+i] += omega(i)*torsdamp*zeta_;      
      }
    }
    */
    

  }
  

  //computing linear stiffness matrix
  if (stiffmatrix != NULL)
  {
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





    
    //artificial isotropic rotational damping stiffness
    {
      double torsdamp = 0.000005; //0.005 is obviously sufficient for free fluctionations with 10 elements
      LINALG::Matrix<3,1> newangle;
      quaterniontoangle(Qnew_, newangle);
      LINALG::Matrix<3,3> Hinverse = Hinv(newangle);


      Hinverse.Scale(zeta_*0.5*torsdamp / params.get<double>("delta time",0.01));

      for(int i= 0; i<3; i++)
      {
        for(int j=0;j<3;j++)
        {
          (*stiffmatrix)(3+i, 3+j) += Hinverse(i,j);
          (*stiffmatrix)(9+i, 9+j) += Hinverse(i,j);
          (*stiffmatrix)(9+i, 3+j) += Hinverse(i,j);
          (*stiffmatrix)(3+i, 9+j) += Hinverse(i,j);
        }
      }
    }
    
    
    /*
    //artificial isotropic curvature damping stiffness
    {
      double torsdamp = 0.05; //0.005 is obviously sufficient for free fluctionations with 10 elements

      LINALG::Matrix<3,3> artstiff = Tnew;


      artstiff.Scale(zeta_*torsdamp / (params.get<double>("delta time",0.01) * lrefe_ ));

      for(int i= 0; i<3; i++)
      {
        for(int j=0;j<3;j++)
        {
          (*stiffmatrix)(3+i, 3+j) += artstiff(j,i);
          (*stiffmatrix)(9+i, 9+j) -= artstiff(j,i);
          (*stiffmatrix)(9+i, 3+j) += artstiff(j,i);
          (*stiffmatrix)(3+i, 9+j) -= artstiff(j,i);
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
