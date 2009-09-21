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
  else if (action=="calc_brownian")       act = Beam2r::calc_brownian;
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
    case Beam2r::calc_brownian:
    {

      /*calculation of Brownian forces and damping is based on drag coefficient; this coefficient per unit
       * length is approximated by the one of an infinitely long staff for friciton orthogonal to staff axis*/

       LINALG::Matrix<2,2> xrefe;
       for (int k=0; k<2; ++k)
       {
         xrefe(0,k) = Nodes()[k]->X()[0];
         xrefe(1,k) = Nodes()[k]->X()[1];
       }

       //this derivation need to be checked (hering)
       lrefe_  = pow( pow(xrefe(1,1)-xrefe(1,0),2) + pow(xrefe(0,1)-xrefe(0,0),2) , 0.5 );

       //this section refers to previous implementation, where local forces were not stored in class variable
       //but in columnmap/rowmap vector in paramlist
           /*
           // get element displacements (for use in shear flow fields)
           RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
           if (disp==null) dserror("Cannot get state vector 'displacement'");
           vector<double> mydisp(lm.size());
           DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
             */  
          
                
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
           //if(this->Owner() != discretization.Comm().MyPID()) return 0;
           
      
       
       //compute stochastic forces in local frame      
       const int nnode = NumNode();
          switch(nnode)//switch due to order of Ansatzfunctions
          {
              case 2:ComputeLocalBrownianForces<2>(params); break;//linear
              case 3:ComputeLocalBrownianForces<3>(params); break;//quadratic
              case 4:ComputeLocalBrownianForces<4>(params); break;//cubic     
              case 5:ComputeLocalBrownianForces<5>(params); break;//quartic     
              default:dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
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

int DRT::ELEMENTS::Beam2r::EvaluateNeumann(ParameterList& params,
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
    const double det = alpha_[ip];
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
 | Compute forces for Brownian motion (public)                                                       06/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode>
void DRT::ELEMENTS::Beam2r::ComputeLocalBrownianForces(ParameterList& params)
{
  /*creating a random generator object which creates random numbers with mean = 0 and standard deviation
   * (2kT/dt)^0,5 with thermal energy kT, time step size dt; using Blitz namespace "ranlib" for random number generation*/
  ranlib::Normal<double> normalGen(0,pow(2.0 * params.get<double>("KT",0.0) / params.get<double>("delta time",0.0),0.5));
  
  //fstoch consists of nnode*2 values, two forces at each node
  LINALG::Matrix<nnode,2> aux;
  for (int i=0; i<nnode; i++)
  {
    for(int j=0; j<2; j++)
    {
      aux(i,j) = normalGen.random();
    }
  }
     
  /*calculation of Brownian forces and damping is based on drag coefficient; this coefficient per unit
  * length is approximated by the one of an infinitely long staff for friciton orthogonal to staff axis*/
  double zeta = 4 * PI * lrefe_ * params.get<double>("ETA",0.0);
  int stochasticorder = params.get<int>("STOCH_ORDER",-2);
  //if no stochasticorder has been chosen explicitly we exit this function without any action
  if(stochasticorder == -2)
    return;
 
  switch(nnode)
  {
    case 2:     
    { 
      switch(stochasticorder)
      {
        case -1:
        {
          //multiply S_loc (cholesky decomposition of C_loc) for C_loc only diagonal
          for (int i=0; i<2; i++)
          {
            floc_[0](i,0) = pow(zeta/2.0,0.5)*aux(0,i);
            floc_[1](i,0) = pow(zeta/2.0,0.5)*aux(1,i);
          }
        }
        break;
        case 0:
        {
          //multiply S_loc(cholesky decomposition of C_loc) for gamma_parallel=gamma_perp
          for (int i=0; i<2; i++)
          {
            floc_[0](i,0) = pow(zeta/3.0,0.5)*aux(0,i);
            floc_[1](i,0) = pow(zeta/12.0,0.5)*aux(0,i)+pow(zeta/4.0,0.5)*aux(1,i);
          }               
        }
        break;
        case 1:
        {        
          //compute damping coefficients consistently to LiTang2004 assuming that actin filament is
          //discretized by one element only and actin diameter 8nm
          /*
          double K_r = 2.94382;
          double K_t = 1.921348;
          double lnp = log(lrefe_ / 0.008);
          double gamma_perp = K_r*4.0*PI*lrefe_*params.get<double>("ETA",0.0)/(lnp - 0.447);       
          double gamma_par  = 4.0*PI*lrefe_*params.get<double>("ETA",0.0)/( (1/K_t)*(3*lnp + 0.658) - (1/K_r)*(lnp - 0.447) );
          */
                  
          double gamma_perp = zeta;    
          double gamma_par  = zeta/2.0;
                    
          floc_[0](0,0) = pow(gamma_par/3.0,0.5)*aux(0,0);
          floc_[0](1,0) = pow(gamma_perp/3.0,0.5)*aux(0,1);
          floc_[1](0,0) = pow(gamma_par/12.0,0.5)*aux(0,0)+pow(gamma_par/4.0,0.5)*aux(1,0);
          floc_[1](1,0) = pow(gamma_perp/12.0,0.5)*aux(0,1)+pow(gamma_perp/4.0,0.5)*aux(1,1);                  
        }
        break;
      }//end switch stochorder
    }//end nnode=2
    break;
    case 3:
    {
      switch(stochasticorder)
      {
        case -1:
        {
          for (int i=0; i<2; i++)
          {
            floc_[0](i,0) = pow(zeta/6.0,0.5)*aux(0,i);
            floc_[1](i,0) = pow(2.0*zeta/3.0,0.5)*aux(1,i);
            floc_[2](i,0) = pow(zeta/6.0,0.5)*aux(2,i);
          }
          
        }
        break;
        case 0:
        {
          for (int i=0; i<2; i++)
          {
            floc_[0](i,0) = pow(2.0*zeta/15.0,0.5)*aux(0,i);
            floc_[1](i,0) = pow(zeta/30.0,0.5)*aux(0,i)+pow(zeta/2.0,0.5)*aux(1,i);
            floc_[2](i,0) = -pow(zeta/120.0,0.5)*aux(0,i)+pow(zeta/72.0,0.5)*aux(1,i)+pow(zeta/9.0,0.5)*aux(2,i);   
          }
        }
        break;
        case 1:
        {
          dserror("high order not with stochorder 1 implemented");
        }
        break;
      }//end switch stochorder
    }//end nnode=3
    break;
    case 4:
    {
      switch(stochasticorder)
      {
        case -1:
        {
          for (int i=0; i<2; i++)
          {
            floc_[0](i,0) = pow(zeta/8.0,0.5)*aux(0,i);
            floc_[1](i,0) = pow(3.0*zeta/8.0,0.5)*aux(1,i);
            floc_[2](i,0) = pow(3.0*zeta/8.0,0.5)*aux(2,i);
            floc_[3](i,0) = pow(zeta/8.0,0.5)*aux(3,i);
          }         
        }
        break;
        case 0:
        {
          for (int i=0; i<2; i++)
          {
            floc_[0](i,0) = pow(8.0*zeta/105.0,0.5)*aux(0,i);
            floc_[1](i,0) = 99.0/128.0*pow(8.0*zeta/105.0,0.5)*aux(0,i)+pow(3483.0*zeta/10240.0,0.5)*aux(1,i);
            floc_[2](i,0) = -9.0/32.0*pow(8.0*zeta/105.0,0.5)*aux(0,i)-4.0/43.0*pow(3483.0*zeta/10240.0,0.5)*aux(1,i)+pow(81.0*zeta/215.0,0.5)*aux(2,i);
            floc_[3](i,0) = 19.0/118.0*pow(8.0*zeta/105.0,0.5)*aux(0,i)-103.0/1161.0*pow(3483.0*zeta/10240.0,0.5)*aux(1,i)+17.0/108.0*pow(81.0*zeta/215.0,0.5)*aux(2,i)+pow(zeta/16.0,0.5)*aux(3,i);
          }
        }
        break;
        case 1:
        {
          dserror("high order not with stochorder 1 implemented");
        }
        break;
      }//end switch stochorder
    }//end nnode=4    
    break;
    case 5:
    {
      switch(stochasticorder)
      {
        case -1:
        {
          for (int i=0; i<2; i++)
          {
          floc_[0](i,0) = pow(7.0*zeta/90.0,0.5)*aux(0,i);
          floc_[1](i,0) = pow(32.0*zeta/90.0,0.5)*aux(1,i);
          floc_[2](i,0) = pow(4.0*zeta/30.0,0.5)*aux(2,i);
          floc_[3](i,0) = pow(32.0*zeta/90.0,0.5)*aux(3,i);
          floc_[4](i,0) = pow(7.0*zeta/90.0,0.5)*aux(4,i);
          }         
        }
        break;
        case 0:
        {
          for (int i=0; i<2; i++)
          {
            floc_[0](i,0) = 396.0/1745.0*pow(zeta,0.5)*aux(0,i);
            floc_[1](i,0) = 268.0/1165.0*pow(zeta,0.5)*aux(0,i)+653.0/1273.0*pow(zeta,0.5)*aux(1,i);
            floc_[2](i,0) = -278.0/2041.0*pow(zeta,0.5)*aux(0,i)-223.0/3124.0*pow(zeta,0.5)*aux(1,i)+915.0/1652.0*pow(zeta,0.5)*aux(2,i);
            floc_[3](i,0) = 131.0/3010.0*pow(zeta,0.5)*aux(0,i)+259.0/3781.0*pow(zeta,0.5)*aux(1,i)-113.0/1099.0*pow(zeta,0.5)*aux(2,i)+515.0/942.0*pow(zeta,0.5)*aux(3,i);
            floc_[4](i,0) = -46.0/2041.0*pow(zeta,0.5)*aux(0,i)+137.0/4666.0*pow(zeta,0.5)*aux(1,i)-176.0/3081.0*pow(zeta,0.5)*aux(2,i)+409.0/4936.0*pow(zeta,0.5)*aux(3,i)+1.0/5.0*pow(zeta,0.5)*aux(4,i);
          }
        }
        break;
        case 1:
        {
          dserror("high order not with stochorder 1 implemented");
        }
        break;
      }//end switch stochorder
    }
    break;
    default:
      dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
  }//end switch nnode

  
return; 
}//DRT::ELEMENTS::Beam2::ComputeLocalBrownianForces
/*-----------------------------------------------------------------------------------------------------------*
 | Assemble statistical forces and damping matrix according to fluctuation dissipation theorem (public) 06/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode>
inline void DRT::ELEMENTS::Beam2r::CalcBrownian(ParameterList& params,
                              vector<double>&           vel,  //!< element velocity vector
                              vector<double>&           disp,
                              Epetra_SerialDenseMatrix* stiffmatrix,  //!< element stiffness matrix
                              Epetra_SerialDenseVector* force)  //!< element internal force vector
{  
  
  
   //get parameters
   double dt = params.get<double>("delta time",0.0);
   int stochasticorder = params.get<int>("STOCH_ORDER",-2);// polynomial order for interpolation of stochastic line load
   double zeta = 4 * PI * lrefe_ * params.get<double>("ETA",0.0);
   //if no stochasticorder has been chosen explicitly we exit this function without any action
   if(stochasticorder == -2)
     return;
     

   
    switch(nnode)
    {
    case 2:     
    { 
     switch(stochasticorder)
       {
       //simple isotropic model of Brownian motion with uncorrelated nodal forces
       case -1:
       {
         
         //calc brownian damp matrix 
         if (stiffmatrix != NULL) // necessary to run stiffmatrix control routine
         {
           for (int i=0; i<2; i++)
           {           
             (*stiffmatrix)(i,i) += zeta/(2.0*dt);
             (*stiffmatrix)(3+i,3+i) += zeta/(2.0*dt);
           }
         }
         
            
         if (force != NULL)
         {           
           //calc internal brownian forces  
           for (int i=0; i<2; i++)
           {
             (*force)(i) +=zeta/(2.0)*vel[i];
             (*force)(3+i) +=zeta/(2.0)*vel[3+i];
           }
                     
           //calc external brownian forces 
           for (int i=0; i<2; i++)
           {
             (*force)(i) -=floc_[0](i,0);
             (*force)(3+i) -=floc_[1](i,0);
           }
          
         }
            
       }//end case -1
       break;
       
       //isotropic model of Brownian motion with correlated forces
       case 0:
       {           
         //calc brownian damp matrix 
         if (stiffmatrix != NULL) // necessary to run stiffmatrix control routine
         {
           for (int i=0; i<2; i++)
           {
           (*stiffmatrix)(i,i) += zeta/(3.0*dt);
           (*stiffmatrix)(3+i,3+i) += zeta/(3.0*dt);
           (*stiffmatrix)(i,3+i) += zeta/(6.0*dt);
           (*stiffmatrix)(3+i,i) += zeta/(6.0*dt);
           }
         } 
               
         //calc int brownian forces
         if (force !=  NULL)
         {
           for (int i=0; i<2; i++)
           {
           (*force)(i) +=zeta/(3.0)*vel[i]+zeta/(6.0)*vel[3+i];
           (*force)(3+i) +=zeta/(6.0)*vel[i]+zeta/(3.0)*vel[3+i];
           }
    
           //calc ext brownian forces
           for (int i=0; i<2; i++)
           {
             (*force)(i) -=floc_[0](i,0);
             (*force)(3+i) -=floc_[1](i,0);
           }
         }         
       }//end case 0
       break;
       //anisotropic model of Brownian motion with correlated nodal forces
       case 1:
       {           
        
         /* This section is to calc the angle for a geometrical midpoint triad due to the angles stored at nodes
          * in dof vector*/
         
         double theta=(disp[2]+disp[5])*0.5;
         
         //triad to rotate local configuration into global
         LINALG::Matrix<2,2> T(true);
         T(0,0) =  cos(theta);
         T(0,1) = -sin(theta);
         T(1,0) =  sin(theta);
         T(1,1) =  cos(theta);
      
          
         
          //compute damping coefficients consistently to LiTang2004 assuming that actin filament is
          //discretized by one element only and actin diameter 8nm
          /*          
          double K_r = 2.94382;
          double K_t = 1.921348;
          double lnp = log(lrefe_ / 0.008);
          double gamma_perp = K_r*4.0*PI*lrefe_*params.get<double>("ETA",0.0)/(lnp - 0.447);       
          double gamma_par  = 4.0*PI*lrefe_*params.get<double>("ETA",0.0)/( (1/K_t)*(3*lnp + 0.658) - (1/K_r)*(lnp - 0.447) );
          */
          
          double gamma_perp = zeta;    
          double gamma_par  = zeta/2.0;


          //local damping matrix
          LINALG::Matrix<2,2> dampbasis(true);
          dampbasis(0,0) = gamma_par;
          dampbasis(1,1) = gamma_perp;


      
          //turning local into global damping matrix (storing intermediate result in variable "aux1")
          LINALG::Matrix<2,2> aux1;
          aux1.Multiply(T,dampbasis);
          dampbasis.MultiplyNT(aux1,T);
        
            
          //calc brownian damp
          if (stiffmatrix != NULL) // necessary to run stiffmatrix control routine
          {
            
            //complete first term due to variation of velocity     
            for(int i=0; i<2; i++)
            {
              for(int j=0; j<2; j++)
              {
                (*stiffmatrix)(i,j) += dampbasis(i,j)/(3.0*dt);
                (*stiffmatrix)(i+3,j+3) += dampbasis(i,j)/(3.0*dt);
                (*stiffmatrix)(i,j+3) += dampbasis(i,j)/(6.0*dt);
                (*stiffmatrix)(i+3,j) += dampbasis(i,j)/(6.0*dt);
              }
            }
         
            //calc second term due to variation of Triad
            LINALG::Matrix<2,2> Spin(true);//Spinmatrix
            Spin(0,1)=-1.0;
            Spin(1,0)=1.0;
           
            //multiply Spin from right side to C_global
            aux1.Multiply(Spin,dampbasis);
           
            //multiply Spin from left side to C_global
            LINALG::Matrix<2,2> aux2(true);
            aux2.Multiply(dampbasis,Spin);
      
            for (int i=0; i<2; i++)//difference
            {
              for (int j=0; j<2; j++)
              {
                aux1(i,j)=aux1(i,j)-aux2(i,j);
              }
            } 
            
            LINALG::Matrix<4,4> aux3(true);
            for(int i=0; i<2; i++)//multiply velcoef, now two nodes
            {
              for(int j=0; j<2; j++)
              {
                aux3(i,j) += aux1(i,j)/(3.0);
                aux3(i+2,j+2) += aux1(i,j)/(3.0);
                aux3(i,j+2) += aux1(i,j)/(6.0);
                aux3(i+2,j) += aux1(i,j)/(6.0);
              }
            }   
           
           LINALG::Matrix<4,1> aux4(true);
           for (int i=0; i<4; i++)//multiply velocity
           {
             for (int j=0; j<2; j++)
             {
               aux4(i,0)+=aux3(i,j)*vel[j];
               aux4(i,0)+=aux3(i,j+2)*vel[j+3];
             }
           }
         
           //add contribution to stiffmatrix
           for (int i=0; i<2; i++)
           {
             (*stiffmatrix)(i,2)+=aux4(i,0)*0.5;
             (*stiffmatrix)(i+3,2)+=aux4(i+2,0)*0.5;
             (*stiffmatrix)(i,5)+=aux4(i,0)*0.5;
             (*stiffmatrix)(i+3,5)+=aux4(i+2,0)*0.5;
           }
          
           //calc ext. stiffness      
           aux1.Multiply(Spin,T);     
               
           
           LINALG::Matrix<4,1>aux5(true);//multiply fstoch with ST, aux5 vector now for two nodes
           for (int i=0; i<2; i++)
           {
             for (int j=0; j<2; j++)
             {
               aux5(i,0)+=aux1(i,j)*floc_[0](j,0);
               aux5(i+2,0)+=aux1(i,j)*floc_[1](j,0);
             }
           }
           //add contribution to stiffmatrix
           for (int i=0; i<2; i++)
           {
             (*stiffmatrix)(i,2)-=aux5(i,0)*0.5;
             (*stiffmatrix)(i+3,2)-=aux5(i+2,0)*0.5;
             (*stiffmatrix)(i,5)-=aux5(i,0)*0.5;
             (*stiffmatrix)(i+3,5)-=aux5(i+2,0)*0.5;
           }
         
           
          }
          
         if (force != NULL)
          {
            //calc internal brownian forces
            LINALG::Matrix<4,4> aux6(true);
            for(int i=0; i<2; i++)
            {
              for(int j=0; j<2; j++)
              {
                aux6(i,j) = dampbasis(i,j)/(3.0);
                aux6(i+2,j+2) = dampbasis(i,j)/(3.0);
                aux6(i,j+2) = dampbasis(i,j)/(6.0);
                aux6(i+2,j) = dampbasis(i,j)/(6.0);
              }
            }
            
            
            for (int i=0; i<2; i++)
            {
              for (int j=0; j<2; j++)
              {
                (*force)(i)+=aux6(i,j)*vel[j];
                (*force)(i)+=aux6(i,j+2)*vel[j+3];
                (*force)(i+3)+=aux6(i+2,j)*vel[j];
                (*force)(i+3)+=aux6(i+2,j+2)*vel[j+3];
              }
            }
      
            
            //calc ext brownian forces
            for (int i=0; i<2; i++)
            {
              for (int j=0; j<2; j++)
              {
                (*force)(i)-=T(i,j)*floc_[0](j,0);
                (*force)(i+3)-=T(i,j)*floc_[1](j,0);
              }
            } 
      
          }//end calc forces
  
       }//end case 1 
       break;
       
      }//switch       
    }//end nnode=2
    break;
    case 3:
    {
      switch(stochasticorder)
      {
        case -1:
        {
          if (stiffmatrix != NULL) // necessary to run stiffmatrix control routine
           {
             for (int i=0; i<2; i++)
             {              
               (*stiffmatrix)(i,i) += zeta/(6.0*dt);
               (*stiffmatrix)(6+i,6+i) += 2.0*zeta/(3.0*dt);
               (*stiffmatrix)(3+i,3+i) += zeta/(6.0*dt);
             }
           } 
           

           //calc int brownian forces
           if (force !=  NULL)
           {
             for (int i=0; i<2; i++)
             {
               (*force)(i) +=zeta/(6.0)*vel[i];
               (*force)(6+i) +=2.0*zeta/(3.0)*vel[6+i];
               (*force)(3+i) +=zeta/(6.0)*vel[3+i];
             }
      
             //calc ext brownian forces
             for (int i=0; i<2; i++)
             {
               (*force)(i) -=floc_[0](i,0);              
               (*force)(6+i) -=floc_[1](i,0);
               (*force)(3+i) -=floc_[2](i,0);
               
             }
           }
        }
        break;
        case 0:
       {           
         //calc brownian damp matrix 
         if (stiffmatrix != NULL) // necessary to run stiffmatrix control routine
         {
           for (int i=0; i<2; i++)
           {
             (*stiffmatrix)(i,i) += 2.0*zeta/(15.0*dt);
             (*stiffmatrix)(6+i,6+i) += 8.0*zeta/(15.0*dt);
             (*stiffmatrix)(3+i,3+i) += 2.0*zeta/(15.0*dt);

             (*stiffmatrix)(i,6+i) += zeta/(15.0*dt);
             (*stiffmatrix)(6+i,i) += zeta/(15.0*dt);
             (*stiffmatrix)(6+i,3+i) += zeta/(15.0*dt);
             (*stiffmatrix)(3+i,6+i) += zeta/(15.0*dt);
            
             (*stiffmatrix)(i,3+i) -= zeta/(30.0*dt);
             (*stiffmatrix)(3+i,i) -= zeta/(30.0*dt);
           }
         } 
               
         //calc int brownian forces
         if (force !=  NULL)
         {
           for (int i=0; i<2; i++)
           {
             (*force)(i) +=2.0*zeta/(15.0)*vel[i]+zeta/(15.0)*vel[6+i]-zeta/30.0*vel[3+i];
             (*force)(6+i) +=zeta/(15.0)*vel[i]+8.0*zeta/(15.0)*vel[6+i]+zeta/15.0*vel[3+i];
             (*force)(3+i) +=2.0*zeta/(15.0)*vel[3+i]+zeta/(15.0)*vel[6+i]-zeta/30.0*vel[i];
           }
    
           //calc ext brownian forces
           for (int i=0; i<2; i++)
           {
             (*force)(i) -=floc_[0](i,0);
             (*force)(6+i) -=floc_[1](i,0);
             (*force)(3+i) -=floc_[2](i,0);
             
           }
         }         
       }//end case 0
        break;
        case 1:
        {
          dserror("high order not with stochorder 1");
        }
        break;
      }//end switch stochorder
    }//end nnode=3
    break;
    case 4:
    {
      switch(stochasticorder)
      {
        case -1:
        {
          if (stiffmatrix != NULL) // necessary to run stiffmatrix control routine
           {
             for (int i=0; i<2; i++)
             {
               (*stiffmatrix)(i,i) += zeta/(8.0*dt);
               (*stiffmatrix)(6+i,6+i) += 3.0*zeta/(8.0*dt);
               (*stiffmatrix)(9+i,9+i) += 3.0*zeta/(8.0*dt);
               (*stiffmatrix)(3+i,3+i) += zeta/(8.0*dt);
             }
           } 
                 
           //calc int brownian forces
           if (force !=  NULL)
           {
             for (int i=0; i<2; i++)
             {
               (*force)(i) +=zeta/(8.0)*vel[i];
               (*force)(6+i) +=3.0*zeta/(8.0)*vel[6+i];
               (*force)(9+i) +=3.0*zeta/(8.0)*vel[9+i];
               (*force)(3+i) +=zeta/(8.0)*vel[3+i];
             }
      
             //calc ext brownian forces
             for (int i=0; i<2; i++)
             {
               (*force)(i) -=floc_[0](i,0);
               (*force)(6+i) -=floc_[1](i,0);
               (*force)(9+i) -=floc_[2](i,0);
               (*force)(3+i) -=floc_[3](i,0);
             }
           }
        }
        break;
        case 0:
       {           
         //calc brownian damp matrix 
         if (stiffmatrix != NULL) // necessary to run stiffmatrix control routine
         {
           for (int i=0; i<2; i++)
           {
             (*stiffmatrix)(i,i) += 8.0*zeta/(105.0*dt);
             (*stiffmatrix)(6+i,6+i) += 27.0*zeta/(70.0*dt);
             (*stiffmatrix)(9+i,9+i) += 27.0*zeta/(70.0*dt);
             (*stiffmatrix)(3+i,3+i) += 8.0*zeta/(105.0*dt);
             
             (*stiffmatrix)(i,6+i) += 33.0*zeta/(560.0*dt);
             (*stiffmatrix)(6+i,i) += 33.0*zeta/(560.0*dt);
             (*stiffmatrix)(9+i,3+i) += 33.0*zeta/(560.0*dt);
             (*stiffmatrix)(3+i,9+i) += 33.0*zeta/(560.0*dt);

             (*stiffmatrix)(i,9+i) -= 3.0*zeta/(140.0*dt);
             (*stiffmatrix)(9+i,i) -= 3.0*zeta/(140.0*dt);
             (*stiffmatrix)(6+i,3+i) -= 3.0*zeta/(140.0*dt);
             (*stiffmatrix)(3+i,6+i) -= 3.0*zeta/(140.0*dt);

             (*stiffmatrix)(i,3+i) += 19.0*zeta/(1680.0*dt);
             (*stiffmatrix)(3+i,i) += 19.0*zeta/(1680.0*dt);

             (*stiffmatrix)(6+i,9+i) -= 27.0*zeta/(560.0*dt);
             (*stiffmatrix)(9+i,6+i) -= 27.0*zeta/(560.0*dt);
           }
         } 
               
         //calc int brownian forces
         if (force !=  NULL)
         {
           for (int i=0; i<2; i++)
           {
             (*force)(i) +=8.0*zeta/(105.0)*vel[i]+33.0*zeta/(560.0)*vel[6+i]-3.0*zeta/140.0*vel[9+i]+19.0*zeta/(1680.0)*vel[3+i];
             (*force)(6+i) +=33.0*zeta/(560.0)*vel[i]+27.0*zeta/(70.0)*vel[6+i]-27.0*zeta/560.0*vel[9+i]-3.0*zeta/(140.0)*vel[3+i];
             (*force)(9+i) +=33.0*zeta/(560.0)*vel[3+i]+27.0*zeta/(70.0)*vel[9+i]-27.0*zeta/560.0*vel[6+i]-3.0*zeta/(140.0)*vel[i];
             (*force)(3+i) +=8.0*zeta/(105.0)*vel[3+i]+33.0*zeta/(560.0)*vel[9+i]-3.0*zeta/140.0*vel[6+i]+19.0*zeta/(1680.0)*vel[i];
           }
           //calc ext brownian forces
           for (int i=0; i<2; i++)
           {
             (*force)(i) -=floc_[0](i,0);
             (*force)(6+i) -=floc_[1](i,0);
             (*force)(9+i) -=floc_[2](i,0);
             (*force)(3+i) -=floc_[3](i,0);
           }
         }         
       }//end case 0
        break;
        case 1:
        {
          dserror("high order not with stochorder 1");
        }
        break;
      }//end switch stochorder
    }//end case 4     
    break;
    case 5:
    {
      switch(stochasticorder)
      {
        case -1:
        {
          if (stiffmatrix != NULL) // necessary to run stiffmatrix control routine
           {
             for (int i=0; i<2; i++)
             {
               (*stiffmatrix)(i,i) += 7.0*zeta/(90.0*dt);
               (*stiffmatrix)(6+i,6+i) += 32.0*zeta/(90.0*dt);
               (*stiffmatrix)(9+i,9+i) += 4.0*zeta/(30.0*dt);
               (*stiffmatrix)(12+i,12+i) += 32.0*zeta/(90.0*dt);
               (*stiffmatrix)(3+i,3+i) += 7.0*zeta/(90.0*dt);
             }
           } 
                 
           //calc int brownian forces
           if (force !=  NULL)
           {
             for (int i=0; i<2; i++)
             {
               (*force)(i) +=7.0*zeta/(90.0)*vel[i];
               (*force)(6+i) +=32.0*zeta/(90.0)*vel[6+i];
               (*force)(9+i) +=4.0*zeta/(30.0)*vel[9+i];
               (*force)(12+i) +=32.0*zeta/(90.0)*vel[12+i];
               (*force)(3+i) +=7.0*zeta/(90.0)*vel[3+i];
             }
      
             //calc ext brownian forces
             for (int i=0; i<2; i++)
             {
               (*force)(i) -=floc_[0](i,0);
               (*force)(6+i) -=floc_[1](i,0);
               (*force)(9+i) -=floc_[2](i,0);
               (*force)(12+i) -=floc_[3](i,0);
               (*force)(3+i) -=floc_[4](i,0);
             }
           }
        }
        break;
        case 0:
       {           
         //calc brownian damp matrix 
         if (stiffmatrix != NULL) // necessary to run stiffmatrix control routine
         {
           for (int i=0; i<2; i++)
           {
             (*stiffmatrix)(i,i) += 146.0*zeta/(2835.0*dt);
             (*stiffmatrix)(3+i,3+i) += 146.0*zeta/(2835.0*dt);            
             (*stiffmatrix)(6+i,6+i) += 128.0*zeta/(405.0*dt);
             (*stiffmatrix)(12+i,12+i) += 128.0*zeta/(405.0*dt);           
             (*stiffmatrix)(9+i,9+i) += 104.0*zeta/(315.0*dt);

             (*stiffmatrix)(i,6+i) += 148.0*zeta/(2835.0*dt);
             (*stiffmatrix)(6+i,i) += 148.0*zeta/(2835.0*dt);
             (*stiffmatrix)(12+i,3+i) += 148.0*zeta/(2835.0*dt);
             (*stiffmatrix)(3+i,12+i) += 148.0*zeta/(2835.0*dt);
                       
             (*stiffmatrix)(i,9+i) -= 29.0*zeta/(945.0*dt);
             (*stiffmatrix)(9+i,i) -= 29.0*zeta/(945.0*dt);
             (*stiffmatrix)(9+i,3+i) -= 29.0*zeta/(945.0*dt);
             (*stiffmatrix)(3+i,9+i) -= 29.0*zeta/(945.0*dt);
                      
             (*stiffmatrix)(6+i,9+i) -= 64.0*zeta/(945.0*dt);
             (*stiffmatrix)(9+i,6+i) -= 64.0*zeta/(945.0*dt);
             (*stiffmatrix)(9+i,12+i) -= 64.0*zeta/(945.0*dt);
             (*stiffmatrix)(12+i,9+i) -= 64.0*zeta/(945.0*dt);
            
             (*stiffmatrix)(i,12+i) += 4.0*zeta/(405.0*dt);
             (*stiffmatrix)(12+i,i) += 4.0*zeta/(405.0*dt);
             (*stiffmatrix)(6+i,3+i) += 4.0*zeta/(405.0*dt);
             (*stiffmatrix)(3+i,6+i) += 4.0*zeta/(405.0*dt);
                                  
             (*stiffmatrix)(6+i,12+i) += 128.0*zeta/(2835.0*dt);
             (*stiffmatrix)(12+i,6+i) += 128.0*zeta/(2835.0*dt);
            
             (*stiffmatrix)(i,3+i) -= 29.0*zeta/(5670.0*dt);
             (*stiffmatrix)(3+i,i) -= 29.0*zeta/(5670.0*dt);
           }
  
         } 
               
         //calc int brownian forces
         if (force !=  NULL)
         {
           for (int i=0; i<2; i++)
           {
             (*force)(i) +=146.0*zeta/(2835.0)*vel[i]+148.0*zeta/(2835.0)*vel[6+i]-29.0*zeta/945.0*vel[9+i]+4.0*zeta/(405.0)*vel[12+i]-29.0*zeta/(5670.0)*vel[3+i];
             (*force)(6+i) +=148.0*zeta/(2835.0)*vel[i]+128.0*zeta/(405.0)*vel[6+i]-64.0*zeta/945.0*vel[9+i]+128.0*zeta/(2835.0)*vel[12+i]+4.0*zeta/(405.0)*vel[3+i];
             (*force)(9+i) +=-29.0*zeta/(945.0)*vel[i]-64.0*zeta/(945.0)*vel[6+i]+104.0*zeta/315.0*vel[9+i]-64.0*zeta/(945.0)*vel[12+i]-29.0*zeta/(945.0)*vel[3+i];
             (*force)(12+i) +=148.0*zeta/(2835.0)*vel[3+i]+128.0*zeta/(405.0)*vel[12+i]-64.0*zeta/945.0*vel[9+i]+128.0*zeta/(2835.0)*vel[6+i]+4.0*zeta/(405.0)*vel[i];
             (*force)(3+i) +=146.0*zeta/(2835.0)*vel[3+i]+148.0*zeta/(2835.0)*vel[12+i]-29.0*zeta/945.0*vel[9+i]+4.0*zeta/(405.0)*vel[6+i]-29.0*zeta/(5670.0)*vel[i];
           }
           
           //calc ext brownian forces
           for (int i=0; i<2; i++)
           {
             (*force)(i) -=floc_[0](i,0);
             (*force)(6+i) -=floc_[1](i,0);
             (*force)(9+i) -=floc_[2](i,0);
             (*force)(12+i) -=floc_[3](i,0);
             (*force)(3+i) -=floc_[4](i,0);
           }
         }         
       }//end case 0
        break;
        case 1:
        {
          dserror("high order not with stochorder 1");
        }
        break;
      }//end switch stochorder
    }  
    break;
    default:
      dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
    }//end switch nnode
 
  return;

}//DRT::ELEMENTS::Beam2::CalcBrownian



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
void DRT::ELEMENTS::Beam2r::nlnstiffmass( ParameterList& params,
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
  double ym; //Young's modulus
  double sm; //shear modulus
  double density; //density

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
    strains(0)=(cos(theta_gp)*dxdxi+sin(theta_gp)*dzdxi)/alpha_[numgp]-1;
    // chi at the GP according to Crisfield Vol.1 (7.131)
    strains(1)=theta_gp_deriv/alpha_[numgp];
    // gamma at the GP according to Crisfield Vol.1 (7.133)
      strains(2)=(-sin(theta_gp)*dxdxi+cos(theta_gp)*dzdxi)/alpha_[numgp];

      //Crisfield Vol.1 (7.139)
      const double c1 = -alpha_[numgp]*(1 + strains(0));

      //Crisfield Vol.1 (7.140)
      const double c2 = alpha_[numgp]*strains(2);

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

        //linear elastic part of tangential stiffness matrix including rotation: K_t1 = B^T C_l B / alpha following (7.144)
        for(int id_col=0; id_col<3*nnode; id_col++)
        {
          aux_CB(0,id_col) = Bcurr_gp(0,id_col) * (ym*crosssec_/alpha_[numgp]);// aux_CB oben definiert. Nur als Ablage für den Zwischenwert da (3x3*nnode)
          aux_CB(1,id_col) = Bcurr_gp(1,id_col) * (ym*mominer_/alpha_[numgp]);
          aux_CB(2,id_col) = Bcurr_gp(2,id_col) * (sm*crosssecshear_/alpha_[numgp]);
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
         (*massmatrix)(i,j)+= alphamass_[numgp]*wgt*massmatrixgp(i,j);
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
   CalcBrownian<nnode>(params,vel,disp,stiffmatrix,force);
   
     
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


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_BEAM2R

