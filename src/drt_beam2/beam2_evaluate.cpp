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

#include "beam2.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../drt_inpar/inpar_statmech.H"

/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                                 cyron 01/08|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam2::Evaluate(Teuchos::ParameterList& params,
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
  else if (action=="calc_struct_ptcstiff")        act = Beam2::calc_struct_ptcstiff;
  else dserror("Unknown type of action for Beam2");

  switch(act)
  {
    case Beam2::calc_struct_ptcstiff:
    {
      dserror("Beam2 element does'nt need any special ptc tools to allow stable implicit dynamics with acceptable time step size");
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
           
      //only if random numbers for Brownian dynamics are passed to element, get element velocities
      vector<double> myvel(lm.size());
      if( params.get<  RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null) != Teuchos::null)
      {
        RefCountPtr<const Epetra_Vector> vel  = discretization.GetState("velocity");      
        DRT::UTILS::ExtractMyValues(*vel,myvel,lm);
      }

      // determine element matrices and forces
      if (act == Beam2::calc_struct_nlnstiffmass)
        nlnstiffmass(params,lm,myvel,mydisp,&elemat1,&elemat2,&elevec1,lumpedmass);
      else if (act == Beam2::calc_struct_nlnstifflmass)
      {
        nlnstiffmass(params,lm,myvel,mydisp,&elemat1,&elemat2,&elevec1,lumpedmass);
      }
      else if (act == Beam2::calc_struct_nlnstiff)
        nlnstiffmass(params,lm,myvel,mydisp,&elemat1,NULL,&elevec1,lumpedmass);
      else if  (act ==  calc_struct_internalforce)
        nlnstiffmass(params,lm,myvel,mydisp,NULL,NULL,&elevec1,lumpedmass);

      /*at the end of an iteration step the geometric ocnfiguration has to be updated: the starting point for the
       * next iteration step is the configuration at the end of the current step */
      numperiodsold_ = numperiodsnew_;
      alphaold_ = alphanew_;


      //the following code block can be used to check quickly whether the nonlinear stiffness matrix is calculated
      //correctly or not by means of a numerically approximated stiffness matrix
      /*if(Id() == 2) //limiting the following tests to certain element numbers
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
             vel_aux[i + 3*k] =  vel_aux[i + 3*k] + h_rel * params.get<double>("gamma",1.0) / ( params.get<double>("delta time",0.01)*params.get<double>("beta",1.0) );

            //std::cout << " approximation force calc " <<  i << k;
             //calc forces depend on vel_aux and disp_aux
            nlnstiffmass(params,lm,vel_aux,disp_aux,NULL,NULL,&force_aux,lumpedmass);



            //calc approx stiffmatrix
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
         std::cout<<"\n\n actually calculated stiffness matrix"<< elemat1;
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
       * status of the beam have to be updated;*/
      numperiodsconv_ = numperiodsnew_;
      alphaconv_ = alphanew_;
    }
    case calc_struct_reset_istep:
    {
      /*the action calc_struct_reset_istep is called by the adaptive time step controller; carries out one test
       * step whose purpose is only figuring out a suitabel timestep; thus this step may be a very bad one in order
       * to iterated towards the new dynamic equilibrium and the thereby gained new geometric configuration should
       * not be applied as starting point for any further iteration step; as a consequence the thereby generated change
       * of the geometric configuration should be canceled and the configuration should be reset to the value at the
       * beginning of the time step*/
      numperiodsold_ = numperiodsconv_;
      alphaold_ = alphaconv_;
      numperiodsnew_ = numperiodsconv_;
      alphanew_ = alphaconv_;
    }
    break;
    case calc_struct_stress:
      dserror("No stress output implemented for beam2 elements");
    default:
      dserror("Unknown type of action for Beam2 %d", act);
  }
  return 0;

}//DRT::ELEMENTS::Beam2::Evaluate


/*-----------------------------------------------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)                                       cyron 01/08|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Beam2::EvaluateNeumann(Teuchos::ParameterList&    params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1,
                                           Epetra_SerialDenseMatrix* elemat1)
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
    curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

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
 | compute current rotation absolute rotation angle of element frame out of x-axisrr              cyron 03/09|
 *----------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Beam2::updatealpha(const LINALG::Matrix<3,2>& xcurr,const double& lcurr)
{
  /*befor computing the absolute rotation angle of the element frame we first compute an angle beta \in [-PI;PI[
   * from the current nodal positions; beta denotes a rotation out of the x-axis in a x-y-plane; note that
   * this angle may differ from the absolute rotation angle alpha by a multiple of 2*PI;*/
  double beta;

  // beta is the rotation angle out of x-axis in a x-y-plane
  double cos_beta = (xcurr(0,1)-xcurr(0,0))/lcurr;
  double sin_beta = (xcurr(1,1)-xcurr(1,0))/lcurr;

  //computation of beta according to Crisfield, Vol. 1, (7.60)

  //if coc_beta >= 0 we know -PI/2 <= beta <= PI/2
  if (cos_beta >= 0)
    beta = asin(sin_beta);
  //else we know  beta > PI/2 or beta < -PI/2
  else
  { 
    //if sin_beta >=0 we know beta > PI/2
    if(sin_beta >= 0)
      beta =  acos(cos_beta);
    //elss we know beta > -PI/2
    else
      beta = -acos(cos_beta);
  }

  /* by default we assume that the difference between beta and the absolute rotation angle alpha is the same
   * multiple of 2*PI as in the iteration step before; then beta + numperiodsnew_*2*PI would be the new absolute
   * rotation angle alpha; if the difference between this angle and the absolute angle in the last converged step
   * is smaller than minus PI we assume that beta, which is evaluated in [-PI; PI[ has exceeded the upper limit of
   * this interval in positive direciton from the last to this iteration step; then alpha can be computed from
   * beta by adding (numperiodsnew_ + 1)*2*PI; analogously with a difference greater than +PI we assume that beta
   * has exceeded the lower limit of the interval [-PI; PI[ in negative direction so that alpha can be computed
   * adding (numperiodsnew_ - 1)*2*PI  */
  numperiodsnew_ = numperiodsold_;

  if(beta + numperiodsnew_*2*PI - alphaold_ < -PI)
    numperiodsnew_ ++;
  else if(beta + numperiodsnew_*2*PI - alphaold_ > M_PI)
    numperiodsnew_ --;


  alphanew_ = beta + 2*M_PI*numperiodsnew_;

}

/*-----------------------------------------------------------------------------------------------------------*
 | evaluate auxiliary vectors and matrices for corotational formulation                           cyron 01/08|
 *----------------------------------------------------------------------------------------------------------*/
//notation for this function similar to Crisfield, Volume 1;
inline void DRT::ELEMENTS::Beam2::local_aux(LINALG::Matrix<3,6>& Bcurr,
                                        LINALG::Matrix<6,1>&     rcurr,
                                        LINALG::Matrix<6,1>&     zcurr,
                                        const double&            lcurr,
                                        const double&            lrefe_)
{
  double cos_alpha = cos(alphanew_);
  double sin_alpha = sin(alphanew_);

  //vector r according to Crisfield, Vol. 1, (7.62)
  rcurr(0) = -cos_alpha;
  rcurr(1) = -sin_alpha;
  rcurr(2) = 0;
  rcurr(3) = cos_alpha;
  rcurr(4) = sin_alpha;
  rcurr(5) = 0;

  //vector z according to Crisfield, Vol. 1, (7.66)
  zcurr(0) = sin_alpha;
  zcurr(1) = -cos_alpha;
  zcurr(2) = 0;
  zcurr(3) = -sin_alpha;
  zcurr(4) = cos_alpha;
  zcurr(5) = 0;

  //assigning values to each element of the Bcurr matrix, Crisfield, Vol. 1, (7.99)
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
} /* DRT::ELEMENTS::Beam2::local_aux */

/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   cyron 01/08|
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Beam2::nlnstiffmass(Teuchos::ParameterList&   params,
                                        vector<int>&              lm,
                                        vector<double>&           vel,
                                        vector<double>&           disp,
                                        Epetra_SerialDenseMatrix* stiffmatrix,
                                        Epetra_SerialDenseMatrix* massmatrix,
                                        Epetra_SerialDenseVector* force,
                                        int                       lumpedmass)
{
  const int numdf = 3;
  const int iel = NumNode();
  //coordinates in current configuration of all the nodes in two dimensions stored in 3 x iel matrices
  LINALG::Matrix<3,2> xcurr;

  //current length of beam in physical space
  double lcurr = 0;

  //some geometric auxiliary variables according to Crisfield, Vol. 1
  LINALG::Matrix<6,1> zcurr;
  LINALG::Matrix<6,1> rcurr;
  LINALG::Matrix<3,6> Bcurr;
  //auxiliary matrix storing the product of constitutive matrix C and Bcurr
  LINALG::Matrix<3,6> aux_CB;
  //declaration of local internal forces
  LINALG::Matrix<3,1> force_loc;
  //declaration of material parameters
  double ym = 0; //Young's modulus
  double sm = 0; //shear modulus
  double density = 0; //density

  //calculating refenrence configuration xrefe and current configuration xcurr
  for (int k=0; k<iel; ++k)
  {
    xcurr(0,k) = Nodes()[k]->X()[0] + disp[k*numdf+0];
    xcurr(1,k) = Nodes()[k]->X()[1] + disp[k*numdf+1];

    /*note that xcurr(2,0),xcurr(2,1) are local angles; in Crisfield, Vol. 1, (7.98) they
     *are denoted by theta_{l1},theta_{l2} in contrast to the local  angles theta_{1},theta_{2};
     *the global director angle is not used at all in the present element formulation*/
    xcurr(2,k) = disp[k*numdf+2];
  }

  //current length
  lcurr = std::pow( pow(xcurr(0,1)-xcurr(0,0),2) + pow(xcurr(1,1)-xcurr(1,0),2) , 0.5 );

  //update absolute rotation angle alpha of element frame
  updatealpha(xcurr,lcurr);

  //calculation of local geometrically important matrices and vectors
  local_aux(Bcurr,rcurr,zcurr,lcurr,lrefe_);


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

  //Crisfield, Vol. 1, (7.52 - 7.55)
  force_loc(0) = ym*crosssec_*(lcurr*lcurr - lrefe_*lrefe_)/(lrefe_*(lcurr + lrefe_));

  //local internal bending moment, Crisfield, Vol. 1, (7.97)
  force_loc(1) = -ym*mominer_*(xcurr(2,1)-xcurr(2,0))/lrefe_;

  //local internal shear force, Crisfield, Vol. 1, (7.98)
  /*in the following internal forces are computed; note that xcurr(2,0),xcurr(2,1) are global angles;
   * in Crisfield, Vol. 1, (7.98) they are denoted by theta_{1},theta_{2} in contrast to the local
   * angles theta_{l1},theta_{l2}; note also that these variables do not represent the absolute director
   * angle (since they are zero in the beginning even for an initially rotated beam), but only
   * the absolute director angle minus the reference angle alpha0_; as a consequence the shear force
   * has to be computed by substraction of (alphanew_ - alpha0_)*/
  force_loc(2) = -sm*crosssecshear_*( ( xcurr(2,1) + xcurr(2,0) )/2 - (alphanew_ - alpha0_) );

  if (force != NULL)
  {
  //declaration of global internal force
  LINALG::Matrix<6,1> force_glob;
  //calculation of global internal forces from Crisfield, Vol. 1, (7.102): q_i = B^T q_{li}
  force_glob.MultiplyTN(Bcurr,force_loc);

    for(int k = 0; k<6; k++)
      (*force)(k) = force_glob(k);
  }


  //calculating tangential stiffness matrix in global coordinates, Crisfield, Vol. 1, (7.107)
  if (stiffmatrix != NULL)
  {
    //declaration of fixed size matrix for global tiffness
    LINALG::Matrix<6,6> stiff_glob;

    //linear elastic part including rotation: B^T C_t B / l_0
    for(int id_col=0; id_col<6; id_col++)
    {
      aux_CB(0,id_col) = Bcurr(0,id_col) * (ym*crosssec_/lrefe_);
      aux_CB(1,id_col) = Bcurr(1,id_col) * (ym*mominer_/lrefe_);
      aux_CB(2,id_col) = Bcurr(2,id_col) * (sm*crosssecshear_/lrefe_);
    }

    stiff_glob.MultiplyTN(aux_CB,Bcurr);

    //adding geometric stiffness by shear force: N z z^T / l_n
    double aux_Q_fac = force_loc(2)*lrefe_ / pow(lcurr,2);
    for(int id_lin=0; id_lin<6; id_lin++)
        for(int id_col=0; id_col<6; id_col++)
        {
          stiff_glob(id_lin,id_col) -= aux_Q_fac * rcurr(id_lin) * zcurr(id_col);
          stiff_glob(id_lin,id_col) -= aux_Q_fac * rcurr(id_col) * zcurr(id_lin);
        }

    //adding geometric stiffness by axial force: Q l_0 (r z^T + z r^T) / (l_n)^2
    double aux_N_fac = force_loc(0)/lcurr;
    for(int id_lin=0; id_lin<6; id_lin++)
        for(int id_col=0; id_col<6; id_col++)
          stiff_glob(id_lin,id_col) += aux_N_fac * zcurr(id_lin) * zcurr(id_col);

    //shfting values from fixed size matrix to epetra matrix *stiffmatrix
    for(int i = 0; i < 6; i++)
      for(int j = 0; j < 6; j++)
        (*stiffmatrix)(i,j) = stiff_glob(i,j);

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


 /*the following function call applied statistical forces and damping matrix according to the fluctuation dissipation theorem;
   * it is dedicated to the application of beam2 elements in the frame of statistical mechanics problems; for these problems a
   * special vector has to be passed to the element packed in the params parameter list; in case that the control routine calling
   * the element does not attach this special vector to params the following method is just doing nothing, which means that for
   * any ordinary problem of structural mechanics it may be ignored*/
   CalcBrownian<2,2,3,2>(params,vel,disp,stiffmatrix,force);


  return;
} // DRT::ELEMENTS::Beam2::nlnstiffmass

/*-----------------------------------------------------------------------------------------------------------*
 |computes the number of different random numbers required in each time step for generation of stochastic    |
 |forces;                                                                    (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam2::HowManyRandomNumbersINeed()
{
  /*at each Gauss point one needs for each node as many random numbers as randomly excited degrees of freedom,
   *i.e. two random numbers for the translational degrees of freedom*/
  return (2*NumNode());

}

/*-----------------------------------------------------------------------------------------------------------*
 | computes damping coefficients per lengthand stores them in a matrix in the following order: damping of    |
 | translation parallel to filament axis, damping of translation orthogonal to filament axis, damping of     |
 | rotation around filament axis                                             (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Beam2::MyDampingConstants(Teuchos::ParameterList& params,LINALG::Matrix<3,1>& gamma, const INPAR::STATMECH::FrictionModel& frictionmodel)
{  
  //translational damping coefficients according to Howard, p. 107, table 6.2;
  gamma(0) = 2*PI*params.get<double>("ETA",0.0);
  gamma(1) = 4*PI*params.get<double>("ETA",0.0);
  
  /*no rotation around element axis possible in 2D; now damping coefficient specified for this motion*/
  gamma(2) = 0;
  
  //in case of an isotropic friction model the same damping coefficients are applied parallel to the polymer axis as perpendicular to it
  if(frictionmodel == INPAR::STATMECH::frictionmodel_isotropicconsistent || frictionmodel == INPAR::STATMECH::frictionmodel_isotropiclumped)
    gamma(0) = gamma(1);

  
   /* in the following section damping coefficients are replaced by those suggested in LiTang2004 assuming that actin filament is
   //discretized by one element only and actin diameter 8nm*/
   /*
   double lrefe=0;
   for (int gp=0; gp<nnode-1; gp++)
     lrefe += gausspointsdamping.qwgt[gp]*jacobi_[gp];
     
   double K_r = 2.94382;
   double K_t = 1.921348;
   
   double p=lrefe/(0.008);
   
   gamma(0) = 4.0*PI*lrefe_*params.get<double>("ETA",0.0)/( (1/K_t)*(3*lnp + 0.658) - (1/K_r)*(lnp - 0.447) );
   gamma(1) = K_r*4.0*PI*lrefe_*params.get<double>("ETA",0.0)/(lnp - 0.447);

   */    
}


/*-----------------------------------------------------------------------------------------------------------*
 |computes velocity of background fluid and gradient of that velocity at a certain evaluation point in       |
 |the physical space                                                         (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int ndim> //number of dimensions of embedding space
void DRT::ELEMENTS::Beam2::MyBackgroundVelocity(Teuchos::ParameterList&       params,  //!<parameter list
                                                const LINALG::Matrix<ndim,1>& evaluationpoint,  //!<point at which background velocity and its gradient has to be computed
                                                LINALG::Matrix<ndim,1>&       velbackground,  //!< velocity of background fluid
                                                LINALG::Matrix<ndim,ndim>&    velbackgroundgrad) //!<gradient of velocity of background fluid
{
  
  /*note: this function is not yet a general one, but always assumes a shear flow, where the velocity of the
   * background fluid is always directed in x-direction. In 3D the velocity increases linearly in z and equals zero for z = 0.
   * In 2D the velocity increases linearly in y and equals zero for y = 0. */
  
  velbackground.PutScalar(0);
  velbackground(0) = evaluationpoint(ndim-1) * params.get<double>("CURRENTSHEAR",0.0);
  
  velbackgroundgrad.PutScalar(0);
  velbackgroundgrad(0,ndim-1) = params.get<double>("CURRENTSHEAR",0.0);

}

/*-----------------------------------------------------------------------------------------------------------*
 | computes translational damping forces and stiffness (public)                                 cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim, int dof> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node
inline void DRT::ELEMENTS::Beam2::MyTranslationalDamping(Teuchos::ParameterList& params,  //!<parameter list
                                                  const vector<double>&          vel,  //!< element velocity vector
                                                  const vector<double>&          disp, //!<element disp vector
                                                  Epetra_SerialDenseMatrix*      stiffmatrix,  //!< element stiffness matrix
                                                  Epetra_SerialDenseVector*      force)//!< element internal force vector
{  
  //get time step size
  double dt = params.get<double>("delta time",0.0);
  
  //velocity and gradient of background velocity field
  LINALG::Matrix<ndim,1> velbackground;
  LINALG::Matrix<ndim,ndim> velbackgroundgrad;
  
  //evaluation point in physical space corresponding to a certain Gauss point in parameter space
  LINALG::Matrix<ndim,1> evaluationpoint;
  
  //get friction model according to which forces and damping are applied
  INPAR::STATMECH::FrictionModel frictionmodel = DRT::INPUT::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL");

  //damping coefficients for translational and rotatinal degrees of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma,frictionmodel);
  
  //get vector jacobi with Jacobi determinants at each integration point (gets by default those values required for consistent damping matrix)
  vector<double> jacobi(jacobimass_);
  
  //determine type of numerical integration performed (lumped damping matrix via lobatto integration!)
  IntegrationType integrationtype = gaussexactintegration;
  if(frictionmodel == INPAR::STATMECH::frictionmodel_isotropiclumped)
  {
    integrationtype = lobattointegration;
    jacobi = jacobinode_;
  }
  
  //get Gauss points and weights for evaluation of damping matrix
  DRT::UTILS::IntegrationPoints1D gausspoints(MyGaussRule(nnode,integrationtype));
  
  //matrix to store basis functions and their derivatives evaluated at a certain Gauss point
  LINALG::Matrix<1,nnode> funct;
  LINALG::Matrix<1,nnode> deriv;


  for(int gp=0; gp < gausspoints.nquad; gp++)
  {    
    //evaluate basis functions and their derivatives at current Gauss point
    DRT::UTILS::shape_function_1D(funct,gausspoints.qxg[gp][0],Shape());
    DRT::UTILS::shape_function_1D_deriv1(deriv,gausspoints.qxg[gp][0],Shape());
     
    //compute point in phyiscal space corresponding to Gauss point
    evaluationpoint.PutScalar(0);
    //loop over all line nodes
    for(int i=0; i<nnode; i++)
      //loop over all dimensions
      for(int j=0; j<ndim; j++)
        evaluationpoint(j) += funct(i)*(Nodes()[i]->X()[j]+disp[dof*i+j]);
    
    //compute velocity and gradient of background flow field at evaluationpoint
    MyBackgroundVelocity<ndim>(params,evaluationpoint,velbackground,velbackgroundgrad);

 
    //compute tangent vector t_{\par} at current Gauss point
    LINALG::Matrix<ndim,1> tpar(true);
    for(int i=0; i<nnode; i++)
      for(int k=0; k<ndim; k++)
        tpar(k) += deriv(i)*(Nodes()[i]->X()[k]+disp[dof*i+k]) / jacobi[gp];
    
    //compute velocity vector at this Gauss point
    LINALG::Matrix<ndim,1> velgp(true);
    for(int i=0; i<nnode; i++)
      for(int l=0; l<ndim; l++)
        velgp(l) += funct(i)*vel[dof*i+l]; 
    
    //compute matrix product (t_{\par} \otimes t_{\par}) \cdot velbackgroundgrad
    LINALG::Matrix<ndim,ndim> tpartparvelbackgroundgrad(true);
    for(int i=0; i<ndim; i++)
      for(int j=0; j<ndim; j++)
        for(int k=0; k<ndim; k++)
          tpartparvelbackgroundgrad(i,j) += tpar(i)*tpar(k)*velbackgroundgrad(k,j);
        
    //loop over all line nodes
    for(int i=0; i<nnode; i++)            
      //loop over lines of matrix t_{\par} \otimes t_{\par}
      for(int k=0; k<ndim; k++)
        //loop over columns of matrix t_{\par} \otimes t_{\par}
        for(int l=0; l<ndim; l++)           
        {               
          if(force != NULL)
            (*force)(i*dof+k)+= funct(i)*jacobi[gp]*gausspoints.qwgt[gp]*( (k==l)*gamma(1) + (gamma(0) - gamma(1))*tpar(k)*tpar(l) ) *(velgp(l)- velbackground(l));
          
          if(stiffmatrix != NULL)
            //loop over all column nodes
            for (int j=0; j<nnode; j++) 
            {
              (*stiffmatrix)(i*dof+k,j*dof+l) += gausspoints.qwgt[gp]*funct(i)*funct(j)*jacobi[gp]*(                 (k==l)*gamma(1) + (gamma(0) - gamma(1))*tpar(k)*tpar(l) ) / dt;
              (*stiffmatrix)(i*dof+k,j*dof+l) -= gausspoints.qwgt[gp]*funct(i)*funct(j)*jacobi[gp]*( velbackgroundgrad(k,l)*gamma(1) + (gamma(0) - gamma(1))*tpartparvelbackgroundgrad(k,l) ) ;             
              (*stiffmatrix)(i*dof+k,j*dof+k) += gausspoints.qwgt[gp]*funct(i)*deriv(j)*                                                   (gamma(0) - gamma(1))*tpar(l)*(velgp(l) - velbackground(l));
              (*stiffmatrix)(i*dof+k,j*dof+l) += gausspoints.qwgt[gp]*funct(i)*deriv(j)*                                                   (gamma(0) - gamma(1))*tpar(k)*(velgp(l) - velbackground(l));    
            }    
        }   
  }
 
  return;
}//DRT::ELEMENTS::Beam3::MyTranslationalDamping(.)

/*-----------------------------------------------------------------------------------------------------------*
 | computes stochastic forces and resulting stiffness (public)                                  cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
inline void DRT::ELEMENTS::Beam2::MyStochasticForces(Teuchos::ParameterList& params,  //!<parameter list
                                              const vector<double>&          vel,  //!< element velocity vector
                                              const vector<double>&          disp, //!<element disp vector
                                              Epetra_SerialDenseMatrix*      stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseVector*      force)//!< element internal force vector
{
  //get friction model according to which forces and damping are applied
  INPAR::STATMECH::FrictionModel frictionmodel = DRT::INPUT::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL");
  
  //damping coefficients for three translational and one rotatinal degree of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma,frictionmodel);
  

  //get vector jacobi with Jacobi determinants at each integration point (gets by default those values required for consistent damping matrix)
  vector<double> jacobi(jacobimass_);
  
  //determine type of numerical integration performed (lumped damping matrix via lobatto integration!)
  IntegrationType integrationtype = gaussexactintegration;
  if(frictionmodel == INPAR::STATMECH::frictionmodel_isotropiclumped)
  {
    integrationtype = lobattointegration;
    jacobi = jacobinode_;
  }
  
  //get Gauss points and weights for evaluation of damping matrix
  DRT::UTILS::IntegrationPoints1D gausspoints(MyGaussRule(nnode,integrationtype));
  
  //matrix to store basis functions and their derivatives evaluated at a certain Gauss point
  LINALG::Matrix<1,nnode> funct;
  LINALG::Matrix<1,nnode> deriv;
  
  
  /*get pointer at Epetra multivector in parameter list linking to random numbers for stochastic forces with zero mean
   * and standard deviation (2*kT / dt)^0.5; note carefully: a space between the two subsequal ">" signs is mandatory
   * for the C++ parser in order to avoid confusion with ">>" for streams*/
   RCP<Epetra_MultiVector> randomnumbers = params.get<  RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null);

  for(int gp=0; gp < gausspoints.nquad; gp++)
  {
    //evaluate basis functions and their derivatives at current Gauss point
    DRT::UTILS::shape_function_1D(funct,gausspoints.qxg[gp][0],Shape());
    DRT::UTILS::shape_function_1D_deriv1(deriv,gausspoints.qxg[gp][0],Shape());
    
    //compute tangent vector t_{\par} at current Gauss point
    LINALG::Matrix<ndim,1> tpar(true);
    for(int i=0; i<nnode; i++)
      for(int k=0; k<ndim; k++)
        tpar(k) += deriv(i)*(Nodes()[i]->X()[k]+disp[dof*i+k]) / jacobi[gp];
     
    
    //loop over all line nodes
    for(int i=0; i<nnode; i++)             
      //loop dimensions with respect to lines
      for(int k=0; k<ndim; k++)
        //loop dimensions with respect to columns
        for(int l=0; l<ndim; l++)           
        {
          if(force != NULL)
            (*force)(i*dof+k) -= funct(i)*(sqrt(gamma(1))*(k==l) + (sqrt(gamma(0)) - sqrt(gamma(1)))*tpar(k)*tpar(l))*(*randomnumbers)[gp*randompergauss+l][LID()]*sqrt(jacobi[gp]*gausspoints.qwgt[gp]);          

          if(stiffmatrix != NULL)
            //loop over all column nodes
            for (int j=0; j<nnode; j++) 
            {            
              (*stiffmatrix)(i*dof+k,j*dof+k) -= funct(i)*deriv(j)*tpar(l)*(*randomnumbers)[gp*randompergauss+l][LID()]*sqrt(gausspoints.qwgt[gp]/ jacobi[gp])*(sqrt(gamma(0)) - sqrt(gamma(1)));   
              (*stiffmatrix)(i*dof+k,j*dof+l) -= funct(i)*deriv(j)*tpar(k)*(*randomnumbers)[gp*randompergauss+l][LID()]*sqrt(gausspoints.qwgt[gp]/ jacobi[gp])*(sqrt(gamma(0)) - sqrt(gamma(1)));  
            }
        }  
  }

  return;
}//DRT::ELEMENTS::Beam2r::MyStochasticForces(.)



/*-----------------------------------------------------------------------------------------------------------*
 | Assemble stochastic and viscous forces and respective stiffness according to fluctuation dissipation      |
 | theorem                                                                               (public) cyron 10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
inline void DRT::ELEMENTS::Beam2::CalcBrownian(Teuchos::ParameterList&  params,
                                              const vector<double>&     vel,  //!< element velocity vector
                                              const vector<double>&     disp, //!< element displacement vector
                                              Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseVector* force) //!< element internal force vector
{   
  //if no random numbers for generation of stochastic forces are passed to the element no Brownian dynamics calculations are conducted
  if( params.get<  RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null) == Teuchos::null)
    return;
  
  //add stiffness and forces due to translational damping effects
  MyTranslationalDamping<nnode,ndim,dof>(params,vel,disp,stiffmatrix,force); 
 
  //add stochastic forces and (if required) resulting stiffness
  MyStochasticForces<nnode,ndim,dof,randompergauss>(params,vel,disp,stiffmatrix,force);

return;

}//DRT::ELEMENTS::Beam2::CalcBrownian(.)



