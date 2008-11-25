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
  else if (action=="calc_stat_force_damp")      act = Beam2::calc_stat_force_damp;
  else if (action=="calc_struct_ptcstiff")        act = Beam2::calc_struct_ptcstiff;
  else dserror("Unknown type of action for Beam2");

  switch(act)
  {
    case Beam2::calc_struct_ptcstiff:
    {
      //Beam2 element does'nt need any special ptc tools to allow stable implicit dynamics with acceptable time step size
    }
    break;
    //action type for evaluating statistical forces
    case Beam2::calc_stat_force_damp:
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
       * thereby gained random forces later on also to DOF of which it is not the row map owner; so fist
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

      //evaluation of statistical forces with the local statistical forces vector fstat
      Epetra_SerialDenseVector fstat(lm.size());
      EvaluateStatForceDamp(params,mydisp,fstat,elemat1);

      /*all the above evaluated forces are used in assembly of the column map force vector no matter if the current processor if
       * the row map owner of these DOF*/
      //note carefully: a space between the two subsequal ">" signs is mandatory for the C++ parser in order to avoid confusion with ">>" for streams
      RCP<Epetra_Vector>    fstatcol = params.get<  RCP<Epetra_Vector> >("statistical force vector",Teuchos::null);

      for(unsigned int i = 0; i < lm.size(); i++)
      {
        //note: lm contains the global Ids of the degrees of freedom of this element
        //testing whether the fstatcol vector has really an element related with the i-th element of fstat by the i-the entry of lm
        if (!(fstatcol->Map()).MyGID(lm[i])) dserror("Sparse vector fstatcol does not have global row %d",lm[i]);

        //get local Id of the fstatcol vector related with a certain element of fstat
        int lid = (fstatcol->Map()).LID(lm[i]);

        //add to the related element of fstatcol the contribution of fstat
        (*fstatcol)[lid] += fstat[i];
      }

    }
    break;
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
      int lumpedmass = 0;  // 0=consistent, 1=lumped
      if (act==Beam2::calc_struct_nlnstifflmass) lumpedmass = 1;

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

      // determine element matrices and forces
      if (act == Beam2::calc_struct_nlnstiffmass)
        b2_nlnstiffmass(params,myvel,mydisp,&elemat1,&elemat2,&elevec1,lumpedmass);
      else if (act == Beam2::calc_struct_nlnstifflmass)
      {
        b2_nlnstiffmass(params,myvel,mydisp,&elemat1,&elemat2,&elevec1,lumpedmass);
      }
      else if (act == Beam2::calc_struct_nlnstiff)
        b2_nlnstiffmass(params,myvel,mydisp,&elemat1,NULL,&elevec1,lumpedmass);
      else if  (act ==  calc_struct_internalforce)
        b2_nlnstiffmass(params,myvel,mydisp,NULL,NULL,&elevec1,lumpedmass);

      /*
      //the following code block can be used to check quickly whether the nonlinear stiffness matrix is calculated
      //correctly or not by means of a numerically approximated stiffness matrix
      if(Id() == 3) //limiting the following tests to certain element numbers
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
            //save current configuration:
            int hrsave = hrnew_;


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


            b2_nlnstiffmass(params,vel_aux,disp_aux,NULL,NULL,&force_aux,lumpedmass);


            for(int u = 0;u<6;u++)
            {
              stiff_approx(u,i+k*3)= ( pow(force_aux[u],2) - pow(elevec1(u),2) )/ (h_rel * (force_aux[u] + elevec1(u) ) );
            }

            //reset geometric ocnfiguration before  computation of approximated stiffness:
            hrnew_ = hrsave;
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

      hrold_ = hrnew_;
    }
    break;
    case calc_struct_update_istep:
    case calc_struct_update_imrlike:
    {
      hrconv_ = hrnew_;
    }
    break;
    case calc_struct_reset_istep:
    {
      hrold_ = hrconv_;
    }
    break;
    case calc_struct_stress:
      dserror("No stress output implemented for beam2 elements");
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
  // element displacements
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

  return 0;
}

/*-----------------------------------------------------------------------------------------------------------*
 | Evaluate Statistical forces and viscous damping (public)                                       cyron 09/08|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Beam2::EvaluateStatForceDamp(ParameterList& params,
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

  //in case of a lumped damping matrix stochastic forces are applied analogously
  if (stochasticorder == 0)
  {
    //adding entries for lumped viscous damping "stiffness" (by background fluid of thermal bath)
    elemat1(0,0) += zeta/2.0;
    elemat1(1,1) += zeta/2.0;
    elemat1(3,3) += zeta/2.0;
    elemat1(4,4) += zeta/2.0;

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
    //adding entries for consistent viscous damping "stiffness" (by background fluid of thermal bath)
    elemat1(0,0) += zeta/3.0;
    elemat1(1,1) += zeta/3.0;
    elemat1(3,3) += zeta/3.0;
    elemat1(4,4) += zeta/3.0;

    elemat1(0,3) += zeta/6.0;
    elemat1(3,0) += zeta/6.0;
    elemat1(1,4) += zeta/6.0;
    elemat1(4,1) += zeta/6.0;

    //calculating standard deviation of statistical forces according to fluctuation dissipation theorem
    double stand_dev_trans = pow(2.0 * kT * (zeta/6.0) / params.get<double>("delta time",0.01),0.5);

    //creating a random generator object which creates random numbers with mean = 0 and standard deviation
    //stand_dev; using Blitz namespace "ranlib" for random number generation
    ranlib::Normal<double> normalGen(0,stand_dev_trans);

    //adding uncorrelated components of statistical forces
    elevec1[0] += normalGen.random();
    elevec1[1] += normalGen.random();
    elevec1[3] += normalGen.random();
    elevec1[4] += normalGen.random();

    //adding correlated components of statistical forces
    double force1 = normalGen.random();
    double force2 = normalGen.random();
    elevec1[0] += force1;
    elevec1[1] += force2;
    elevec1[3] += force1;
    elevec1[4] += force2;
  }

  return 0;
} //DRT::ELEMENTS::Beam3::EvaluateStatisticalNeumann

/*-----------------------------------------------------------------------------------------------------------*
 | evaluate auxiliary vectors and matrices for corotational formulation                           cyron 01/08|
 *----------------------------------------------------------------------------------------------------------*/
//notation for this function similar to Crisfield, Volume 1;
inline void DRT::ELEMENTS::Beam2::b2_local_aux(LINALG::Matrix<3,6>& Bcurr,
                    			              LINALG::Matrix<6,1>& rcurr,
                    			              LINALG::Matrix<6,1>& zcurr,
                                        double& beta,
                                        const LINALG::Matrix<3,2>& xcurr,
                                        const double& lcurr,
                                        const double& lrefe_)

{

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
  double cos_psi = pow(-1.0,hrold_)*(xcurr(0,1)-xcurr(0,0))/lcurr;
  double sin_psi = pow(-1.0,hrold_)*(xcurr(1,1)-xcurr(1,0))/lcurr;

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

  beta = psi + PI*hrold_ - beta0_;

  if (psi > PI/2)
	  hrnew_ = hrold_ +1;
  if (psi < -PI/2)
    hrnew_ = hrold_ -1;

  rcurr(0) = -cos_beta;
  rcurr(1) = -sin_beta;
  rcurr(2) = 0;
  rcurr(3) = cos_beta;
  rcurr(4) = sin_beta;
  rcurr(5) = 0;

  zcurr(0) = sin_beta;
  zcurr(1) = -cos_beta;
  zcurr(2) = 0;
  zcurr(3) = -sin_beta;
  zcurr(4) = cos_beta;
  zcurr(5) = 0;

  //assigning values to each element of the Bcurr matrix
  for(int id_col=0; id_col<6; id_col++)
  	{
  	  Bcurr(0,id_col) = rcurr(id_col);
  	  Bcurr(1,id_col) = 0;
  	  Bcurr(2,id_col) = (lrefe_ / lcurr) * zcurr(id_col);
  	}
    Bcurr(2,2) -= (lrefe_ / 2);
    Bcurr(2,5) -= (lrefe_ / 2);
    Bcurr(1,2) += 1;
    Bcurr(1,5) -= 1;

  return;
} /* DRT::ELEMENTS::Beam2::b2_local_aux */

/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   cyron 01/08|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam2::b2_nlnstiffmass( ParameterList& params,
                                            vector<double>&           vel,
                                            vector<double>&           disp,
                                            Epetra_SerialDenseMatrix* stiffmatrix,
                                            Epetra_SerialDenseMatrix* massmatrix,
                                            Epetra_SerialDenseVector* force,
                                            int lumpedmass)
{
  const int numdf = 3;
  const int iel = NumNode();
  //coordinates in current configuration of all the nodes in two dimensions stored in 3 x iel matrices
  LINALG::Matrix<3,2> xcurr;

  //current length of beam in physical space
  double lcurr = 0;
  //current angle between x-axis and beam in physical space
  double beta;
  //some geometric auxiliary variables according to Crisfield, Vol. 1
  LINALG::Matrix<6,1> zcurr;
  LINALG::Matrix<6,1> rcurr;
  LINALG::Matrix<3,6> Bcurr;
  //auxiliary matrix storing the product of constitutive matrix C and Bcurr
  LINALG::Matrix<3,6> aux_CB;
  //declaration of local internal forces
  LINALG::Matrix<3,1> force_loc;
  //declaration of material parameters
  double ym; //Young's modulus
  double sm; //shear modulus
  double density; //density

  //calculating refenrence configuration xrefe and current configuration xcurr
  for (int k=0; k<iel; ++k)
  {
    xcurr(0,k) = Nodes()[k]->X()[0] + disp[k*numdf+0];
    xcurr(1,k) = Nodes()[k]->X()[1] + disp[k*numdf+1];
    //note: this is actually not the current director angle, but current director angle minus reference director angle
    xcurr(2,k) = disp[k*numdf+2];
  }

  // calculation of local geometrically important matrices and vectors; notation according to Crisfield-------
  //current length
  lcurr = pow( pow(xcurr(0,1)-xcurr(0,0),2) + pow(xcurr(1,1)-xcurr(1,0),2) , 0.5 );

  //calculation of local geometrically important matrices and vectors
  b2_local_aux(Bcurr, rcurr, zcurr, beta, xcurr, lcurr, lrefe_);


  /* read material parameters using structure _MATERIAL which is defined by inclusion of      /
  / "../drt_lib/drt_timecurve.H"; note: material parameters have to be read in the evaluation /
  / function instead of e.g. beam2_input.cpp or within the Beam2Register class since it is not/
  / sure that structure _MATERIAL is declared within those scopes properly whereas it is within/
  / the evaluation functions */

  // get the material law
  MATERIAL* currmat = &(mat[material_-1]);
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

  if (force != NULL)
  {
  //declaration of global internal force
  LINALG::Matrix<6,1> force_glob;
  //calculation of global internal forces from force = B_transposed*force_loc
  force_glob.MultiplyTN(Bcurr,force_loc);

    for(int k = 0; k<6; k++)
      (*force)(k) = force_glob(k);
  }


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
    BLITZTINY::MtM_product<6,6,3>(aux_CB,Bcurr,*stiffmatrix);

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

  //calculating mass matrix (local version = global version)
  if (massmatrix != NULL)
  {
      //if lumped_flag == 0 a consistent mass Timoshenko beam mass matrix is applied
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
      /*if lumped_flag == 1 a lumped mass matrix is applied where the cross sectional moment of inertia is
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
} // DRT::ELEMENTS::Beam2::b2_nlnstiffmass


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_BEAM2

