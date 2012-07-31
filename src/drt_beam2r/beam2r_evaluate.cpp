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

#include "beam2r.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/stvenantkirchhoff.H"

/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                                 cyron 01/08|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam2r::Evaluate(Teuchos::ParameterList& params,
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
  else if (action=="calc_struct_ptcstiff")        act = Beam2r::calc_struct_ptcstiff;
  else dserror("Unknown type of action for Beam2r");

  switch(act)
  {
    case Beam2r::calc_struct_ptcstiff:
    {
      dserror("Beam2r element does'nt need any special ptc tools to allow stable implicit dynamics with acceptable time step size");
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
      
      //only if random numbers for Brownian dynamics are passed to element, get element velocities
      vector<double> myvel(lm.size());
      if( params.get<  RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null) != Teuchos::null)
      {
        RefCountPtr<const Epetra_Vector> vel  = discretization.GetState("velocity");      
        DRT::UTILS::ExtractMyValues(*vel,myvel,lm);
      }
      
      
      const int nnode = NumNode();
      // determine element matrices and forces
      // nlinstiffmass is templated. Therefore we need to give the number of nodes to the function
      if (act == Beam2r::calc_struct_nlnstiffmass)
      {
        switch(nnode)
        {
            case 2:
            {
              nlnstiffmass<2>(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
              break;
            }
            case 3:
            {
              nlnstiffmass<3>(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
              break;
            }
            case 4:
            {
              nlnstiffmass<4>(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
              break;
            }
            case 5:
            {
              nlnstiffmass<5>(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
              break;
            }
            default:
              dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
        }
      }
      else if (act == Beam2r::calc_struct_nlnstifflmass)
      {
        switch(nnode)
        {
            case 2:
            {
              nlnstiffmass<2>(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
              lumpedmass<2>(&elemat2);
              break;
            }
            case 3:
            {
              nlnstiffmass<3>(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
              lumpedmass<3>(&elemat2);
              break;
            }
            case 4:
            {
              nlnstiffmass<4>(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
              lumpedmass<4>(&elemat2);
              break;
            }
            case 5:
            {
              nlnstiffmass<5>(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
              lumpedmass<5>(&elemat2);
              break;
            }
            default:
              dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
        }
      }
      else if (act == Beam2r::calc_struct_nlnstiff)
      {
        switch(nnode)
        {
            case 2:
            {
              nlnstiffmass<2>(params,myvel,mydisp,&elemat1,NULL,&elevec1);
              break;
            }
            case 3:
            {
              nlnstiffmass<3>(params,myvel,mydisp,&elemat1,NULL,&elevec1);
              break;
            }
            case 4:
            {
              nlnstiffmass<4>(params,myvel,mydisp,&elemat1,NULL,&elevec1);
              break;
            }
            case 5:
            {
              nlnstiffmass<5>(params,myvel,mydisp,&elemat1,NULL,&elevec1);
              break;
            }
            default:
              dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
        }
      }
      else if  (act ==  calc_struct_internalforce)
      {
        switch(nnode)
        {
            case 2:
            {
              nlnstiffmass<2>(params,myvel,mydisp,NULL,NULL,&elevec1);
              break;
            }
            case 3:
            {
              nlnstiffmass<3>(params,myvel,mydisp,NULL,NULL,&elevec1);
              break;
            }
            case 4:
            {
              nlnstiffmass<4>(params,myvel,mydisp,NULL,NULL,&elevec1);
              break;
            }
            case 5:
            {
              nlnstiffmass<5>(params,myvel,mydisp,NULL,NULL,&elevec1);
              break;
            }
            default:
              dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
        }
      }




      //the following code block can be used to check quickly whether the nonlinear stiffness matrix is calculated
      //correctly or not by means of a numerically approximated stiffness matrix
      //The code block will work for all higher order elements.
     /* if(Id() == 3) //limiting the following tests to certain element numbers
      {
        //variable to store numerically approximated stiffness matrix
        Epetra_SerialDenseMatrix stiff_approx;
        stiff_approx.Shape(3*nnode,3*nnode);

        //relative error of numerically approximated stiffness matrix
        Epetra_SerialDenseMatrix stiff_relerr;
        stiff_relerr.Shape(3*nnode,3*nnode);

        //characteristic length for numerical approximation of stiffness
        double h_rel = 1e-8;

        //flag indicating whether approximation leads to significant relative error
        int outputflag = 0;

        //calculating strains in new configuration
        for(int i=0; i<3; i++) //for all dof
        {
          for(int k=0; k<nnode; k++)//for all nodes
          {

            Epetra_SerialDenseVector force_aux;
            force_aux.Size(3*nnode);

            //create new displacement and velocity vectors in order to store artificially modified displacements
            vector<double> vel_aux(3*nnode);
            vector<double> disp_aux(3*nnode);

              DRT::UTILS::ExtractMyValues(*disp,disp_aux,lm);
              DRT::UTILS::ExtractMyValues(*vel,vel_aux,lm);

            //modifying displacement artificially (for numerical derivative of internal forces):
            disp_aux[3*k + i] += h_rel;
            vel_aux[3*k + i] += h_rel * params.get<double>("gamma",0.581) / ( params.get<double>("delta time",0.01)*params.get<double>("beta",0.292) );
        //nlnstiffmass is a templated function. therefore we need to point out the number of nodes in advance
              switch(nnode)
              {
                  case 2:
                  {
                    nlnstiffmass<2>(params,vel_aux,disp_aux,NULL,NULL,&force_aux);
                    break;
                  }
                  case 3:
                  {
                    nlnstiffmass<3>(params,vel_aux,disp_aux,NULL,NULL,&force_aux);
                    break;
                  }
                  case 4:
                  {
                    nlnstiffmass<4>(params,vel_aux,disp_aux,NULL,NULL,&force_aux);
                    break;
                  }
                  case 5:
                  {
                    nlnstiffmass<5>(params,vel_aux,disp_aux,NULL,NULL,&force_aux);
                    break;
                  }
                  default:
                    dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
              }


            for(int u = 0 ; u < 3*nnode ; u++ )
            {
              stiff_approx(u,i+k*3)= ( pow(force_aux[u],2) - pow(elevec1(u),2) )/ (h_rel * (force_aux[u] + elevec1(u) ) );
            }

          } //for(int k=0; k<nnode; k++)//for all nodes

        } //for(int i=0; i<3; i++) //for all dof

        for(int line=0; line<3*nnode; line++)
        {
          for(int col=0; col<3*nnode; col++)
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

int DRT::ELEMENTS::Beam2r::EvaluateNeumann(Teuchos::ParameterList& params,
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

  // no. of nodes on this element
  const int nnode = NumNode();
  const int numdf = 3;
  const DiscretizationType distype = this->Shape();

  // gaussian points
  const DRT::UTILS::IntegrationPoints1D  intpoints(gaussrule_);


  //declaration of variable in order to store shape function
  Epetra_SerialDenseVector      funct(nnode);

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
    const double det = jacobi_[ip];
    //evaluation of shape functions at Gauss points
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
    for (int node=0; node<nnode; ++node)
      for (int dof=0; dof<numdf; ++dof)
         elevec1[node*numdf+dof] += funct[node] *ar[dof];

  } // for (int ip=0; ip<intpoints.nquad; ++ip)

  return 0;
}

/*-----------------------------------------------------------------------------------------------------------*
 | evaluate auxiliary vectors and matrices for Reissner`s formulation (private)                   cyron 01/08|
 *----------------------------------------------------------------------------------------------------------*/
//notation for this function similar to Crisfield, Volume 1, section 7.4;
template<int nnode>
inline void DRT::ELEMENTS::Beam2r::local_aux(LINALG::Matrix<3,3*nnode>& Bcurr_gp,
                                        LINALG::Matrix<3*nnode,1>& rcurr_gp,
                                        LINALG::Matrix<3*nnode,1>& zcurr_gp,
                                        LINALG::Matrix<3*nnode,1>& scurr_gp,
                                        const double& theta_gp,
                                        const double& c1,
                                        const double& c2,
                                        LINALG::Matrix<1,nnode>& funct,
                                        LINALG::Matrix<1,nnode>& deriv)

{

  //midline helps to set up row 1 of Bcurr
  LINALG::Matrix<3*nnode,1> midline;

  //set up vectors rcurr_gp, zcurr_gp, scurr_gp and midline for current configuration at GP
  for(int id_col=0; id_col< nnode; id_col++)
  {
    //set up rcurrgp according to Crisfield Vol.1 (7.136)
    rcurr_gp(3*id_col)=cos(theta_gp)*deriv(id_col);
    rcurr_gp(3*id_col+1)=sin(theta_gp)*deriv(id_col);
    rcurr_gp(3*id_col+2)=0.0;

    //set up zcurrgp according to Crisfield Vol.1 (7.137)
    zcurr_gp(3*id_col)= -sin(theta_gp)*deriv(id_col);
    zcurr_gp(3*id_col+1)=cos(theta_gp)*deriv(id_col);
    zcurr_gp(3*id_col+2)=0.0;

    //set up s according to Crisfield Vol.1 (7.141)
    scurr_gp(3*id_col)=0.0;
    scurr_gp(3*id_col+1)=0.0;
    scurr_gp(3*id_col+2)=funct(id_col);

    //midline helps to fill up row 1 of Bcurr
    midline(3*id_col)=0;
    midline(3*id_col+1)=0;
    midline(3*id_col+2)=deriv(id_col);

  }

  //assigning values to each element of the Bcurr matrix, Crisfield, Vol. 1, (7.135)
  for(int id_col=0; id_col<3*nnode; id_col++)
  {
      Bcurr_gp(0,id_col) = rcurr_gp(id_col) + c2 * scurr_gp(id_col);
      Bcurr_gp(1,id_col) = midline(id_col);
      Bcurr_gp(2,id_col) = zcurr_gp(id_col) + c1 * scurr_gp(id_col);
  }

  return;
} /* DRT::ELEMENTS::Beam2r::local_aux */

/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   cyron 01/08|
 *-----------------------------------------------------------------------------------------------------------*/
template<int nnode>
void DRT::ELEMENTS::Beam2r::nlnstiffmass( Teuchos::ParameterList& params,
                                            vector<double>&           vel,
                                            vector<double>&           disp,
                                            Epetra_SerialDenseMatrix* stiffmatrix,
                                            Epetra_SerialDenseMatrix* massmatrix,
                                            Epetra_SerialDenseVector* force)
{

  const int numdf = 3;

  /*coordinates in current configuration of all the nodes in two dimensions stored in numdf x nnode matrices
   *
   * [x1 x2 ...]
   * [z1 z2 ...]
   * [O1 O2 ...]
   *
   * */
  LINALG::Matrix<3,nnode> xcurr;

  //some geometric auxiliary variables at GP according to Crisfield, Vol. 1 section 7.4
  LINALG::Matrix<3*nnode,1> zcurr_gp;
  LINALG::Matrix<3*nnode,1> rcurr_gp;
  LINALG::Matrix<3*nnode,1> scurr_gp;
  LINALG::Matrix<3,3*nnode> Bcurr_gp;

  //auxiliary matrix storing the product of constitutive matrix C and Bcurr
  LINALG::Matrix<3,3*nnode> aux_CB;

  //declaration of local internal forces at GP
  LINALG::Matrix<3,1> force_loc_gp; //keeps the same dimension for different order of shapefunctions

  //declaration of material parameters
  double ym = 0; //Young's modulus
  double sm = 0; //shear modulus
  double density = 0; //density

  //Inserting current configuration into xcurr
  for (int k=0; k<nnode; ++k)
  {
    xcurr(0,k) = Nodes()[k]->X()[0] + disp[k*numdf+0]; //x for each node in figure 7.9
    xcurr(1,k) = Nodes()[k]->X()[1] + disp[k*numdf+1]; //z for each node  in figure 7.9

    /*note that xcurr(2,k) are local angles in Crisfield, Vol. 1. They refer to the
     * reference configuration. The exact global angle is only needed for the GP to integrate. For our further calculations
     * we will only calculate theta_gp for each GP to integrate. This is possible because in SetUpReferenceGeometry in
     * beam2r.cpp we stored theta0_ at each gausspoint. The correct value at the GP can be interpolated from xcurr(2,k)
     * Furthermore for a stress-free-reference-configuration we only consider the nodal deflections from ref. conf.
     * */

    xcurr(2,k) = disp[k*numdf+2]; //theta_l for each node  in figure 7.9
  }

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

  //Get integrationpoints for the applied gaussrule
  DRT::UTILS::IntegrationPoints1D gausspoints(gaussrule_);

  //Get DiscretizationType
  const DRT::Element::DiscretizationType distype = Shape();

  //Matrices for h and h,xi
  LINALG::Matrix<1,nnode> funct;
  LINALG::Matrix<1,nnode> deriv;

  //Loop through all GP and calculate their contribution to the forcevector and stiffnessmatrix
  for(int numgp=0; numgp < gausspoints.nquad; numgp++)
  {
      //Get location and weight of GP in parameter space
      const double xi = gausspoints.qxg[numgp][0];
      const double wgt = gausspoints.qwgt[numgp];

      //Get h and h,xi
      DRT::UTILS::shape_function_1D(funct,xi,distype);
      DRT::UTILS::shape_function_1D_deriv1(deriv,xi,distype);

      //Variables to store values for each GP
      double  dxdxi(0.0),
          dzdxi(0.0),
          theta_gp(0.0),
          theta_gp_deriv(0.0);

    //calculate current dx/dxi, dz/dxi and theta as well as theta,xi at the gausspoint
    for(int i=0; i<nnode; i++)
    {
        dxdxi+=deriv(i)*xcurr(0,i);
        dzdxi+=deriv(i)*xcurr(1,i);
        theta_gp+=funct(i)*xcurr(2,i);
        theta_gp_deriv+=deriv(i)*xcurr(2,i);
    }

    //Interpolation of current theta at gp theta = theta0 + SumOverNodes( h(GP) * delta_theta )
    theta_gp += theta0_[numgp];


    // store all strains in one Matrix for better access
    LINALG::Matrix<3,1> strains;

      // epsilon at the GP according to Crisfield Vol.1 (7.132)
    strains(0)=(cos(theta_gp)*dxdxi+sin(theta_gp)*dzdxi)/jacobi_[numgp]-1;
    // chi at the GP according to Crisfield Vol.1 (7.131)
    strains(1)=theta_gp_deriv/jacobi_[numgp];
    // gamma at the GP according to Crisfield Vol.1 (7.133)
      strains(2)=(-sin(theta_gp)*dxdxi+cos(theta_gp)*dzdxi)/jacobi_[numgp];

      //Crisfield Vol.1 (7.139)
      const double c1 = -jacobi_[numgp]*(1 + strains(0));

      //Crisfield Vol.1 (7.140)
      const double c2 = jacobi_[numgp]*strains(2);

      //calculation of local geometrically important matrices and vectors
      local_aux<nnode>(Bcurr_gp,rcurr_gp,zcurr_gp,scurr_gp,theta_gp,c1,c2,funct,deriv);

      //Crisfield, Vol. 1, (7.55) and (7.132)
      force_loc_gp(0) = ym * crosssec_ * strains(0);

      //local internal bending moment, Crisfield, Vol. 1, (7.108)
      force_loc_gp(1) = ym * mominer_ * strains(1);

      //local internal shear force, Crisfield, Vol. 1, (7.98) and gamma from (7.133)
      force_loc_gp(2) = sm * crosssecshear_ * strains(2);

      if (force != NULL)
      {
          //declaration of global internal force
          LINALG::Matrix<3*nnode,1> force_glob_gp;

          //calculation of global internal forces from Crisfield, Vol. 1, (7.124): q_i = B^T q_{li}
          force_glob_gp.MultiplyTN(Bcurr_gp,force_loc_gp);

          for(int k = 0; k<3*nnode; k++)
            (*force)(k) += wgt*force_glob_gp(k); // internal global force-vector is completely calculated and returned
      }

  //calculating tangential stiffness matrix in global coordinates, Crisfield, Vol. 1, (7.107)
  if (stiffmatrix != NULL)
    {
        //declaration of fixed size matrix for global stiffness
        LINALG::Matrix<3*nnode,3*nnode> stiff_glob_gp;

        //linear elastic part of tangential stiffness matrix including rotation: K_t1 = B^T C_l B / jacobi following (7.144)
        for(int id_col=0; id_col<3*nnode; id_col++)
        {
          aux_CB(0,id_col) = Bcurr_gp(0,id_col) * (ym*crosssec_/jacobi_[numgp]);// aux_CB oben definiert. Nur als Ablage für den Zwischenwert da (3x3*nnode)
          aux_CB(1,id_col) = Bcurr_gp(1,id_col) * (ym*mominer_/jacobi_[numgp]);
          aux_CB(2,id_col) = Bcurr_gp(2,id_col) * (sm*crosssecshear_/jacobi_[numgp]);
        }

        stiff_glob_gp.MultiplyTN(aux_CB,Bcurr_gp);

        //adding geometric stiffness by axial force: N (s z^T + z s^T) + N * c1 s s^T following (7.144)
        for(int id_lin=0; id_lin<3*nnode; id_lin++)
          for(int id_col=0; id_col<3*nnode; id_col++)
          {
            stiff_glob_gp(id_lin,id_col) += force_loc_gp(0) * ( scurr_gp(id_lin) * zcurr_gp(id_col) + zcurr_gp(id_lin) * scurr_gp(id_col));
            stiff_glob_gp(id_lin,id_col) += force_loc_gp(0) * c1 * scurr_gp(id_lin) * scurr_gp(id_col);
          }

        //adding geometric stiffness by shear force: -Q ( s r^T + r s^T ) - Q * c2 * s s^T following (7.144)
        for(int id_lin=0; id_lin<3*nnode; id_lin++)
          for(int id_col=0; id_col<3*nnode; id_col++)
          {
            stiff_glob_gp(id_lin,id_col) -= force_loc_gp(2) * ( scurr_gp(id_lin) * rcurr_gp(id_col) + rcurr_gp(id_lin) * scurr_gp(id_col));
            stiff_glob_gp(id_lin,id_col) -= force_loc_gp(2) * c2 * scurr_gp(id_lin) * scurr_gp(id_col);
          }

        //shifting values from fixed size matrix to epetra matrix *stiffmatrix
        for(int i = 0; i < 3*nnode; i++)
          for(int j = 0; j < 3*nnode; j++)
            (*stiffmatrix)(i,j) += wgt * stiff_glob_gp(i,j);
        //stiffnessmatrix is here completely calculated for this GP and added up with the weight

    }//if (stiffmatrix != NULL)

  }//for(int numgp=0; numgp < gausspoints.nquad; numgp++)




  //calculating mass matrix (local version = global version)
  //We use a consistent Timoshenko mass matrix here
 if (massmatrix != NULL)
  {
   //The mass matrix has to be integrated completely. Therefore we use more gausspoints
   gaussrule_ = static_cast<enum DRT::UTILS::GaussRule1D>(nnode);

   //Get the applied integrationpoints for complete integration
   DRT::UTILS::IntegrationPoints1D gausspointsmass(gaussrule_);

   //Matrix to store Shape functions as defined throughout the FE lecture
   LINALG::Matrix<3,3*nnode> N;

   //Matrix to store mass matrix at each GP in
   LINALG::Matrix<3*nnode,3*nnode> massmatrixgp;

   for(int numgp=0; numgp < gausspointsmass.nquad; numgp++)
     {
     //Get location and weight of GP in parameter space
     const double xi = gausspointsmass.qxg[numgp][0];
     const double wgt = gausspointsmass.qwgt[numgp];

     //Get h
     DRT::UTILS::shape_function_1D(funct,xi,distype);

     //Set N to zeros
     N.Clear();

     //Fill up N as defined in the FE lecture
     for (int i=0; i<nnode; i++)
     {
       N(0,3*i)=funct(i);
       N(1,3*i+1)=funct(i);
       N(2,3*i+2)=funct(i);
     }

     //m= density*crosssec* integraloverxi [N_t*N]
     massmatrixgp.MultiplyTN(density*crosssec_,N,N);

     for (int i=0; i<3*nnode; i++)
     {
       for (int j=0; j<nnode; j++)
       {
         massmatrixgp(i,3*j+2)= mominer_/crosssec_ * massmatrixgp(i,3*j+2);
         //This function multiplies all entries associated with
         //the rotation theta. The mominer_ comes from calculations considering
         //the moment created by a rotation theta
       }
     }

     //Sum up the massmatrices calculated at each gp using the lengthfactor and the weight
     for (int i=0; i<3*nnode; i++)
     {
       for (int j=0; j<3*nnode; j++)
       {
         (*massmatrix)(i,j)+= jacobimass_[numgp]*wgt*massmatrixgp(i,j);
       }
     }
    }

   //reset gaussrule
   gaussrule_ = static_cast<enum DRT::UTILS::GaussRule1D>(nnode-1);
  }


   
   /*the following function call applied statistical forces and damping matrix according to the fluctuation dissipation theorem;
    * it is dedicated to the application of beam2 elements in the frame of statistical mechanics problems; for these problems a
    * special vector has to be passed to the element packed in the params parameter list; in case that the control routine calling
    * the element does not attach this special vector to params the following method is just doing nothing, which means that for
    * any ordinary problem of structural mechanics it may be ignored*/
    CalcBrownian<nnode,2,3,2>(params,vel,disp,stiffmatrix,force);
   
     
  return;
} // DRT::ELEMENTS::Beam2r::nlnstiffmass


/*-----------------------------------------------------------------------------------------------------------*
 | Transforms consistent massmatrix into a lumped massmatrix          (private)                   cyron 01/08|
 *-----------------------------------------------------------------------------------------------------------*/

template<int nnode>
inline void DRT::ELEMENTS::Beam2r::lumpedmass(Epetra_SerialDenseMatrix* massmatrix)

{
  //lumped massmatric adds up all translational and rotational parts of the consistent massmatrix and
  //spreads them equally on the diagonal

  if (massmatrix != NULL)
  {
    double translational = 0.0;
    double rotational = 0.0;

    for (int i=0; i<3*nnode; i++)
    {
      for (int j=0; j<nnode; j++)
      {
        //Add up and...
        translational+=(*massmatrix)(i,3*j);
        translational+=(*massmatrix)(i,3*j+1);
        rotational+=(*massmatrix)(i,3*j+2);
        //...set to zero afterwards
        (*massmatrix)(i,3*j)=0;
        (*massmatrix)(i,3*j+1)=0;
        (*massmatrix)(i,3*j+2)=0;
      }
    }
    //The diagonal entries are gained from the sum of all entries
    translational= translational/(2.0*nnode);
    rotational=rotational/nnode;

    for (int i=0; i<nnode; i++)
    {
      //enter diagonal entries
      (*massmatrix)(3*i,3*i)= translational;
      (*massmatrix)(3*i+1,3*i+1)= translational;
      (*massmatrix)(3*i+2,3*i+2)= rotational;
    }

  }

  return;
} /* DRT::ELEMENTS::Beam2r::lumpedmass */

/*-----------------------------------------------------------------------------------------------------------*
 |computes the number of different random numbers required in each time step for generation of stochastic    |
 |forces;                                                                    (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam2r::HowManyRandomNumbersINeed()
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
inline void DRT::ELEMENTS::Beam2r::MyDampingConstants(Teuchos::ParameterList& params,LINALG::Matrix<3,1>& gamma, const INPAR::STATMECH::FrictionModel& frictionmodel,const DRT::UTILS::IntegrationPoints1D gausspoints, const vector<double> jacobi)
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
   double lrefe = 0;
   for (int gp=0; gp<gausspoints.nquad; gp++)
     lrefe += gausspoints.qwgt[gp]*jacobi[gp];
    
   double K_r = 0.2861736;
   double K_t = 0.5204678;
   
   double p = lrefe/(0.008);
   
   gamma(0) = 4.0*PI*params.get<double>("ETA",0.0)/(K_t*(3*log(p) + 0.658) - K_r*(log(p) - 0.447));
   gamma(1) = 4.0*PI*params.get<double>("ETA",0.0)/(K_r*(log(p)-0.447));  
   */ 
}


/*-----------------------------------------------------------------------------------------------------------*
 |computes velocity of background fluid and gradient of that velocity at a certain evaluation point in       |
 |the physical space                                                         (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int ndim> //number of dimensions of embedding space
void DRT::ELEMENTS::Beam2r::MyBackgroundVelocity(Teuchos::ParameterList& params,  //!<parameter list
                                                const LINALG::Matrix<ndim,1>& evaluationpoint,  //!<point at which background velocity and its gradient has to be computed
                                                LINALG::Matrix<ndim,1>& velbackground,  //!< velocity of background fluid
                                                LINALG::Matrix<ndim,ndim>& velbackgroundgrad) //!<gradient of velocity of background fluid
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
inline void DRT::ELEMENTS::Beam2r::MyTranslationalDamping(Teuchos::ParameterList& params,  //!<parameter list
                                                  const vector<double>&     vel,  //!< element velocity vector
                                                  const vector<double>&     disp, //!<element disp vector
                                                  Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                                                  Epetra_SerialDenseVector* force)//!< element internal force vector
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
  
  //damping coefficients for translational and rotatinal degrees of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma,frictionmodel,gausspoints,jacobi);
  
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
              (*stiffmatrix)(i*dof+k,j*dof+k) += gausspoints.qwgt[gp]*funct(i)*deriv(j)*                                               (gamma(0) - gamma(1))*tpar(l)*(velgp(l) - velbackground(l));
              (*stiffmatrix)(i*dof+k,j*dof+l) += gausspoints.qwgt[gp]*funct(i)*deriv(j)*                                               (gamma(0) - gamma(1))*tpar(k)*(velgp(l) - velbackground(l));
            }    
        }
  }
  
 
  return;
}//DRT::ELEMENTS::Beam3::MyTranslationalDamping(.)

/*-----------------------------------------------------------------------------------------------------------*
 | computes stochastic forces and resulting stiffness (public)                                  cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
inline void DRT::ELEMENTS::Beam2r::MyStochasticForces(Teuchos::ParameterList& params,  //!<parameter list
                                              const vector<double>&     vel,  //!< element velocity vector
                                              const vector<double>&     disp, //!<element disp vector
                                              Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseVector* force)//!< element internal force vector
{
  //get friction model according to which forces and damping are applied
  INPAR::STATMECH::FrictionModel frictionmodel = DRT::INPUT::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL");
  
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
  
  //damping coefficients for three translational and one rotatinal degree of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma,frictionmodel,gausspoints,jacobi);
  
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
inline void DRT::ELEMENTS::Beam2r::CalcBrownian(Teuchos::ParameterList& params,
                                              const vector<double>&           vel,  //!< element velocity vector
                                              const vector<double>&           disp, //!< element displacement vector
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


