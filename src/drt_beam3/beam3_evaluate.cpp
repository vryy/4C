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
      /*calculation of Brownian forces and damping is based on drag coefficient; this coefficient per unit
       * length is approximated by the one of an infinitely long staff*/
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
      if(Id() == 8) //limiting the following tests to certain element numbers
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
          //std::cout<<"\n\n acutally calculated stiffness matrix"<< elemat1;
          //std::cout<<"\n\n approximated stiffness matrix"<< stiff_approx;
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

int DRT::ELEMENTS::Beam3::EvaluateBrownianDamp(ParameterList& params,
                                           vector<double> mydisp,
                                           double zeta,
                                           Epetra_SerialDenseMatrix& elemat1)
{
  /*
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
  }
  */
  
  
  //local damping matrix
  LINALG::Matrix<3,3> dampbasis(true);
  dampbasis(0,0) = zeta/2;
  dampbasis(1,1) = zeta;
  dampbasis(2,2) = zeta;
  
  LINALG::Matrix<3,3> Tconv;
  quaterniontotriad(Qconv_,Tconv);
  
  //turning local damping matrix into global one
  dampbasis.Multiply(Tconv,dampbasis);
  dampbasis.MultiplyNT(dampbasis,Tconv);


  for(int i = 0; i<3; i++)
  {
    for(int j = 0; j<3; j++)
    {
      elemat1(i,j)     += dampbasis(i,j)/3.0;
      elemat1(i+6,j+6) += dampbasis(i,j)/3.0;
    
      elemat1(i,j+6) += dampbasis(i,j)/6.0;
      elemat1(i+6,j) += dampbasis(i,j)/6.0;
    }
  }
  
  
  return 0;
} //DRT::ELEMENTS::Beam3::EvaluateBrownian


/*-----------------------------------------------------------------------------------------------------------*
 | Evaluate Statistical forces and viscous damping (public)                                       cyron 09/08|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Beam3::EvaluateBrownianForces(ParameterList& params,
                                           vector<double> mydisp,
                                           double zeta,
                                           Epetra_SerialDenseVector& elevec1)
{
  /*
  // thermal energy responsible for statistical forces
  double kT = params.get<double>("KT",0.0);

  // polynomial order for interpolation of stochastic line load (zero corresponds to bead spring model)
  int stochasticorder = params.get<int>("STOCH_ORDER",0);

  //stochastic field of line load is interpolated by zeroth order polynomial functions
  if (stochasticorder == 0)
  {

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
  */
  
  // thermal energy responsible for statistical forces
  double kT = params.get<double>("KT",0.0);

  //local force vectors of first and second node
  LINALG::Matrix<3,1> force1;
  LINALG::Matrix<3,1> force2;
  
  LINALG::Matrix<3,3> Tconv;
  quaterniontotriad(Qconv_,Tconv);
  
 
  double stand_dev_par = pow(2 * kT * (zeta/12) / params.get<double>("delta time",0.01),0.5);
  double stand_dev_ort = pow(2 * kT * (zeta/ 6) / params.get<double>("delta time",0.01),0.5);

  //creating a random generator object which creates random numbers with mean = 0 and standard deviation
  //stand_dev; using Blitz namespace "ranlib" for random number generation
  ranlib::Normal<double> normalGenpar(0,stand_dev_par);
  ranlib::Normal<double> normalGenort(0,stand_dev_ort);

  //adding uncorrelated components of statistical forces
  force1(0) = normalGenpar.random();
  force1(1) = normalGenort.random();
  force1(2) = normalGenort.random();
  force2(0) = normalGenpar.random();
  force2(1) = normalGenort.random();
  force2(2) = normalGenort.random();


  //adding correlated components of statistical forces
  double stat1 = normalGenpar.random();
  double stat2 = normalGenort.random();
  double stat3 = normalGenort.random();
  force1(0) += stat1;
  force1(1) += stat2;
  force1(2) += stat3;
  force2(0) += stat1;
  force2(1) += stat2;
  force2(2) += stat3;
  
  //transform nodal forces into global coordinates
  force1.Multiply(Tconv,force1);
  force2.Multiply(Tconv,force2);
  
  for(int i = 0; i<3; i++)
  {
    elevec1[i]   = force1(i);
    elevec1[i+6] = force2(i);
  }

    

  
  return 0;
} //DRT::ELEMENTS::Beam3::EvaluateBrownianForces



/*-----------------------------------------------------------------------------------------------------------*
 | Evaluate PTC damping (public)                                                                  cyron 10/08|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Beam3::EvaluatePTC(ParameterList& params,
                                      Epetra_SerialDenseMatrix& elemat1)
{
  //get PTC dti parameter
  double dti = params.get<double>("dti",0.0);
  
  double isotorsdamp = 0.0005;
  double anisotorsdamp = 0.0005;  
  
  //computing angle increment from current position in comparison with last converged position for damping
  LINALG::Matrix<4,1> deltaQ;
  LINALG::Matrix<3,1> deltatheta;     
  quaternionproduct(inversequaternion(Qconv_),Qnew_,deltaQ);      
  quaterniontoangle(deltaQ,deltatheta);
           
  //computing special matrix for anisotropic damping
  LINALG::Matrix<3,3> Tconv;
  quaterniontotriad(Qconv_,Tconv);
  LINALG::Matrix<3,3> Theta;
  for(int i = 0; i<3; i++)
    for(int j = 0; j<3; j++)
      Theta(i,j) = Tconv(i,0)*Tconv(j,0);
  
  //inverse exponential map
  LINALG::Matrix<3,3> Hinverse = Hinv(deltatheta);
     
  //isotropic artificial stiffness
  LINALG::Matrix<3,3> artstiff = Hinverse;
  artstiff.Scale(zeta_*0.5*isotorsdamp / params.get<double>("delta time",0.01));    
  
  //anisotropic artificial stiffness      
  LINALG::Matrix<3,3> auxstiff;
  auxstiff.Multiply(Theta,Hinverse);
  auxstiff.Scale(anisotorsdamp*zeta_*0.5 / params.get<double>("delta time",0.01));
  artstiff += auxstiff;
  
  //scale artificial damping with dti parameter for PTC method
  artstiff.Scale(dti);
       

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

  //adding diagonal entries according to first summand of equation (17.70)
  R(0,0) = 1 - 2*(q(1)*q(1) + q(2)*q(2));
  R(1,1) = 1 - 2*(q(0)*q(0) + q(2)*q(2));
  R(2,2) = 1 - 2*(q(0)*q(0) + q(1)*q(1));

  return;
} // DRT::ELEMENTS::Beam3::quaterniontotriad

//!computing a quaternion from an angle vector
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
  
  //computing quaternion equivalent rotation by angle deltabetaplusalpha/2
  LINALG::Matrix<4,1> Qrot;
  angletoquaternion(deltabetaplusalpha,Qrot);

  //computing quaterion Qnew_ for new configuration of Qold_ for old configuration by means of a quaternion product
  quaternionproduct(Qold_,Qrot,Qnew_);
  
  //normalizing quaternion in order to make sure that it keeps unit absolute values throught time stepping
  double abs = Qnew_.Norm2();
  for(int i = 0; i<4; i++)
    Qnew_(i) = Qnew_(i) / abs;

  quaterniontotriad(Qnew_,Tnew);
  

} //DRT::ELEMENTS::Beam3::updatetriad

//computes inverse quaternion q^{-1} for input quaternion q
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

/*quaternion product q12 = q2*q1, Crisfield, Vol. 2, equation (16.71)*/
inline void DRT::ELEMENTS::Beam3::quaternionproduct(const LINALG::Matrix<4,1>& q1,const LINALG::Matrix<4,1>& q2,LINALG::Matrix<4,1>& q12)
{
  q12(0) = q2(3)*q1(0) + q1(3)*q2(0) + q2(1)*q1(2) - q1(1)*q2(2);
  q12(1) = q2(3)*q1(1) + q1(3)*q2(1) + q2(2)*q1(0) - q1(2)*q2(0);
  q12(2) = q2(3)*q1(2) + q1(3)*q2(2) + q2(0)*q1(1) - q1(0)*q2(1);
  q12(3) = q2(3)*q1(3) - q2(2)*q1(2) - q2(1)*q1(1) - q2(0)*q1(0);
} //DRT::ELEMENTS::Beam3::quaternionproduct

//updating local curvature according to Crisfield, Vol. 2, pages 209 - 210; not: an exact update of the curvature is computed by
//means of equation (16.148) instead of an approximated one as given by equs. (17.72) and (17.73)
inline void DRT::ELEMENTS::Beam3::updatecurvature(const LINALG::Matrix<3,3>& Tnew, LINALG::Matrix<3,1> deltabetaplusalpha,LINALG::Matrix<3,1> deltabetaminusalpha)
{
  //-------------------calculating omega-------------------------------------//

  //applying proper scaling for rotation angle
  deltabetaplusalpha.Scale(0.5);
  
  //compute quaternion from angle deltabetaplusalpha/2
  LINALG::Matrix<4,1> q;
  angletoquaternion(deltabetaplusalpha,q);
  quaterniontoangle(q,deltabetaplusalpha);

  //absolute value of rotation vector theta
  double abs_theta = deltabetaplusalpha.Norm2();

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

/*
//updating local curvature according to Crisfield, Vol. 2, pages 209 - 210; not: an exact update of the curvature is computed by
//means of equation (16.148) instead of an approximated one as given by equs. (17.72) and (17.73)
inline void DRT::ELEMENTS::Beam3::updatecurvature(const LINALG::Matrix<3,3>& Tnew, LINALG::Matrix<3,1> deltabetaplusalpha,LINALG::Matrix<3,1> deltabetaminusalpha)
{
  
  //applying proper scaling for rotation angle
  deltabetaplusalpha.Scale(0.5);
  
  //compute quaternion from angle deltabetaplusalpha/2
  LINALG::Matrix<4,1> q;
  angletoquaternion(deltabetaplusalpha,q);
  
  //compute Rodrigues parameters omega for deltabetaplusalpha/2:
  LINALG::Matrix<3,1> omega;
  quaterniontorodrigues(q,omega);
  
  //spatial derivative of Rodrigues parameters according to (16.147)
  LINALG::Matrix<3,1> omegaprime; 
  
  //for use of formula (16.147) we need an angle -PI < theta < PI; an angle in this domain is gained by method quaterniontoangle;
  //note that the resulting angle theta is not necessarily identical to deltabetaplusalpha/2, but that both may differ by a
  //multiple of 2*PI
  LINALG::Matrix<3,1> theta;
  quaterniontoangle(q,theta);
  
  if(theta.Norm2() == 0)
  {
    omegaprime = deltabetaminusalpha;
    omegaprime.Scale(1/lrefe_);
  }
  else
  {
    //compute unit vector e of rotation axis (cf. (16.146) )
    LINALG::Matrix<3,1> e = theta;
    e.Scale(1/theta.Norm2());
      
    //compute matrix e*e^T (cf. (16.147) )
    LINALG::Matrix<3,3> eeT;
    eeT.MultiplyNT(e,e);
    
    //in order to evaluate (16.147) we need to prepare the matrix given there in brackets:
    LINALG::Matrix<3,3> aux = eeT;
    aux.Scale(theta.Norm2() / sin(theta.Norm2()) - 1);
    for(int i = 0; i<3; i++)
      aux(i,i) += 1;
    
    omegaprime.Multiply(aux,deltabetaminusalpha);
    //note: theta cannot be +/- PI/2 since in this case already quaterniontorodrigues would have failed
    omegaprime.Scale( ( 2*tan(theta.Norm2() / 2) )/(lrefe_ * theta.Norm2() ) );   
  }
  
 

  LINALG::Matrix<3,1> curvaux;
  curvaux(0) = 0.5*(omega(1)*omegaprime(2) - omega(2)*omegaprime(1)) ;
  curvaux(1) = 0.5*(omega(2)*omegaprime(0) - omega(0)*omegaprime(2)) ;
  curvaux(2) = 0.5*(omega(0)*omegaprime(1) - omega(1)*omegaprime(0)) ;

  curvaux += omegaprime;
  curvaux.Scale( 1/(1 + 0.25*omega.Norm2()*omega.Norm2() ));

  curvnew_.MultiplyTN(Tnew,curvaux);
  curvnew_ += curvold_;


  return;
} //DRT::ELEMENTS::Beam3::updatecurvature
*/

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
    x21(j)                = (Nodes()[1]->X()[j]  - Nodes()[0]->X()[j] ) + ( disp[6+j]  - disp[j] );
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
  updatecurvature(Tnew,deltabetaplusalpha,deltabetaminusalpha);

  //computing current axial and shear strain epsilon, Crisfield, Vol. 2, equation (17.67)
  epsilonn.MultiplyTN(Tnew,x21);
  epsilonn.Scale(1/lrefe_);
  epsilonn(0) -=  1;


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

    
    //the following block adds artificial rotation damping to force vector and stiffness matrix
    {
      double isotorsdamp = 0; //0.0005 seems a good choice for merely isotropic damping
      double anisotorsdamp = 0.0005;  //0.0005 seems a good choice with isotorsdamp = 0;
      
      //computing angle increment from current position in comparison with last converged position for damping
      LINALG::Matrix<4,1> deltaQ;
      LINALG::Matrix<3,1> deltatheta;     
      quaternionproduct(inversequaternion(Qconv_),Qnew_,deltaQ);      
      quaterniontoangle(deltaQ,deltatheta);
      
      //angular velocity
      LINALG::Matrix<3,1> omega = deltatheta;
      omega.Scale(1/ params.get<double>("delta time",0.01));
              
      //computing special matrix for anisotropic damping
      LINALG::Matrix<3,3> Tconv;
      quaterniontotriad(Qconv_,Tconv);
      LINALG::Matrix<3,3> Theta;
      for(int i = 0; i<3; i++)
        for(int j = 0; j<3; j++)
          Theta(i,j) = Tconv(i,0)*Tconv(j,0);
      
      //inverse exponential map
      LINALG::Matrix<3,3> Hinverse = Hinv(deltatheta);
         
      //isotropic artificial stiffness
      LINALG::Matrix<3,3> artstiff = Hinverse;
      artstiff.Scale(zeta_*0.5*isotorsdamp / params.get<double>("delta time",0.01));    
      
      //anisotropic artificial stiffness      
      LINALG::Matrix<3,3> auxstiff;
      auxstiff.Multiply(Theta,Hinverse);
      auxstiff.Scale(anisotorsdamp*zeta_*0.5 / params.get<double>("delta time",0.01));
      artstiff += auxstiff;
           
      //isotropic artificial forces
      LINALG::Matrix<3,1> artforce = omega;
      artforce.Scale(isotorsdamp*zeta_);
          
      //anisotropic artificial forces  
      LINALG::Matrix<3,1> auxforce;
      auxforce.Multiply(Theta,omega);
      auxforce.Scale(anisotorsdamp*zeta_);
      artforce += auxforce;
      
      //adding artificial contributions to stifness matrix and force vector
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
    }
    
  


  }


  /*calculating mass matrix; this beam3 element includes only a lumped mass matrix for slender beams the 
   * influence of rotational moments of inertia is dilute so that here (as often) rotational inertia is just 
   * set to zero*/
  if (massmatrix != NULL)
  {
    for (int i=0; i<3; ++i)
    {
      (*massmatrix)(i  ,i  ) = 0.5*density*lrefe_*crosssec_;
      (*massmatrix)(i+6,i+6) = 0.5*density*lrefe_*crosssec_;
      
      /*BAIC carries out a LU-decomposition in the very beginning; in order to prevent it from crashing
      *we make sure full rank of the complete stiffenss matrix by adding minor artificial interia entrie
      *for rotational degrees of freedom; note: these are physically not correct, but do not affet much
      * the simulation results*/
      (*massmatrix)(i+3,i+3) = 0.5*density*lrefe_*Iyy_;
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
