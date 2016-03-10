/*!-----------------------------------------------------------------------------------------------------------
 \file beam3r_evaluate.cpp
 \brief

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>

 *-----------------------------------------------------------------------------------------------------------*/

#include "beam3r.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_fem_general/largerotations.H"
#include "../drt_inpar/inpar_statmech.H"
#include "../headers/FAD_utils.H"

#include <iostream>
#include <iomanip>

#include <Teuchos_Time.hpp>

/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                                 cyron 01/08|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3r::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization,
    std::vector<int>& lm,
    Epetra_SerialDenseMatrix& elemat1, //nonlinear stiffness matrix
    Epetra_SerialDenseMatrix& elemat2, //nonlinear mass matrix
    Epetra_SerialDenseVector& elevec1, //nonlinear internal (elastic) forces
    Epetra_SerialDenseVector& elevec2, //nonlinear inertia forces
    Epetra_SerialDenseVector& elevec3)
{
  // set flag whether pure Statmech functionalities are needed in this element evaluation call
  // TODO: maybe we could do this once upon initialization of element and already use in SetUpReferenceGeometry
  // so far, no chance to figure out whether Statmech is needed because no parameter list is available upon initialization
  needstatmech_ = (params.get<Teuchos::RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null) != Teuchos::null);

  DRT::ELEMENTS::Beam3r::ActionType act = Beam3r::calc_none;
  // get the action required
  std::string action = params.get<std::string>("action","calc_none");
  if (action == "calc_none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff") act = Beam3r::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff") act = Beam3r::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = Beam3r::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass") act = Beam3r::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass") act = Beam3r::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass") act = Beam3r::calc_struct_nlnstifflmass; //with lumped mass matrix
  else if (action=="calc_struct_stress") act = Beam3r::calc_struct_stress;
  else if (action=="calc_struct_eleload") act = Beam3r::calc_struct_eleload;
  else if (action=="calc_struct_fsiload") act = Beam3r::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep") act = Beam3r::calc_struct_update_istep;
  else if (action=="calc_struct_reset_istep") act = Beam3r::calc_struct_reset_istep;
  else if (action=="calc_struct_ptcstiff")        act = Beam3r::calc_struct_ptcstiff;
  else if (action=="calc_struct_energy")        act = Beam3r::calc_struct_energy;
  else dserror("Unknown type of action for Beam3r");

  switch(act)
  {
    case Beam3r::calc_struct_ptcstiff:
    {
      const int nnode = NumNode();

      switch(nnode)
       {
           case 2:EvaluatePTC<2>(params, elemat1); break;
           case 3:EvaluatePTC<3>(params, elemat1); break;
           case 4:EvaluatePTC<4>(params, elemat1); break;
           case 5:EvaluatePTC<5>(params, elemat1); break;
           default:dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
       }
    }
    break;
    /*in case that only linear stiffness matrix is required b3_nlstiffmass is called with zero dispalcement and
     residual values*/
    case Beam3r::calc_struct_linstiff:
    {
      //only nonlinear case implemented!
      dserror("linear stiffness matrix called, but not implemented");

    }
    break;
    //calculate internal energy
    case Beam3r::calc_struct_energy:
    {
      // need current global displacement and residual forces and get them from discretization
      // making use of the local-to-global map lm one can extract current displacemnet and residual values for each degree of freedom
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

      const int nnode = NumNode();

      switch(nnode)
      {
        case 2:
        {
          b3_energy<2>(params,mydisp,&elevec1);
          break;
        }
        case 3:
        {
          b3_energy<3>(params,mydisp,&elevec1);
          break;
        }
        case 4:
        {
          b3_energy<4>(params,mydisp,&elevec1);
          break;
        }
        case 5:
        {
          b3_energy<5>(params,mydisp,&elevec1);
          break;
        }
        default:
          dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
      }
    }
    break;

    //nonlinear stiffness and mass matrix are calculated even if only nonlinear stiffness matrix is required
    case Beam3r::calc_struct_nlnstiffmass:
    case Beam3r::calc_struct_nlnstifflmass:
    case Beam3r::calc_struct_nlnstiff:
    case Beam3r::calc_struct_internalforce:
    {
      // need current global displacement and residual forces and get them from discretization
      // making use of the local-to-global map lm one can extract current displacemnet and residual values for each degree of freedom
      //
      // get element displcements
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

      // get residual displacements
      Teuchos::RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (res==Teuchos::null) dserror("Cannot get state vectors 'residual displacement'");
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);

      //only if random numbers for Brownian dynamics are passed to element, get element velocities
      std::vector<double> myvel(lm.size());
      std::vector<double> myacc(lm.size());
      if( needstatmech_)
      {
        Teuchos::RCP<const Epetra_Vector> vel  = discretization.GetState("velocity");
        DRT::UTILS::ExtractMyValues(*vel,myvel,lm);
      }

      if(act == Beam3r::calc_struct_nlnstiffmass or act == Beam3r::calc_struct_nlnstifflmass)
      {
        Teuchos::RCP<const Epetra_Vector> vel  = discretization.GetState("velocity");
        if (vel==Teuchos::null) dserror("Cannot get state vectors 'velocity'");
        DRT::UTILS::ExtractMyValues(*vel,myvel,lm);

        Teuchos::RCP<const Epetra_Vector> acc  = discretization.GetState("acceleration");
        if (acc==Teuchos::null) dserror("Cannot get state vectors 'acceleration'");
        DRT::UTILS::ExtractMyValues(*acc,myacc,lm);
      }

      const int nnode = NumNode();

      if (act == Beam3r::calc_struct_nlnstiffmass)
      {
        switch(nnode)
        {
          case 2:
          {
            b3_nlnstiffmass<2>(params,myacc,myvel,mydisp,&elemat1,&elemat2,&elevec1,&elevec2);
            break;
          }
          case 3:
          {
            b3_nlnstiffmass<3>(params,myacc,myvel,mydisp,&elemat1,&elemat2,&elevec1,&elevec2);
            break;
          }
          case 4:
          {
            b3_nlnstiffmass<4>(params,myacc,myvel,mydisp,&elemat1,&elemat2,&elevec1,&elevec2);
            break;
          }
          case 5:
          {
            b3_nlnstiffmass<5>(params,myacc,myvel,mydisp,&elemat1,&elemat2,&elevec1,&elevec2);
            break;
          }
        }
      }
      else if (act == Beam3r::calc_struct_nlnstifflmass)
      {
        dserror("Lumped masss matrix not implemented for beam3r elements so far!");
      }
      else if (act == Beam3r::calc_struct_nlnstiff)
      {
        switch(nnode)
        {
          case 2:
          {
            b3_nlnstiffmass<2>(params,myacc,myvel,mydisp,&elemat1,NULL,&elevec1,NULL);
            break;
          }
          case 3:
          {
            b3_nlnstiffmass<3>(params,myacc,myvel,mydisp,&elemat1,NULL,&elevec1,NULL);
            break;
          }
          case 4:
          {
            b3_nlnstiffmass<4>(params,myacc,myvel,mydisp,&elemat1,NULL,&elevec1,NULL);
            break;
          }
          case 5:
          {
            b3_nlnstiffmass<5>(params,myacc,myvel,mydisp,&elemat1,NULL,&elevec1,NULL);
            break;
          }
          default:
            dserror("Only Line2, Line3, Line4, and Line5 Elements implemented.");
        }
      }
      else if (act == Beam3r::calc_struct_internalforce)
      {
        switch(nnode)
        {
          case 2:
          {
            b3_nlnstiffmass<2>(params,myacc,myvel,mydisp,NULL,NULL,&elevec1,NULL);
            break;
          }
          case 3:
          {
            b3_nlnstiffmass<3>(params,myacc,myvel,mydisp,NULL,NULL,&elevec1,NULL);
            break;
          }
          case 4:
          {
            b3_nlnstiffmass<4>(params,myacc,myvel,mydisp,NULL,NULL,&elevec1,NULL);
            break;
          }
          case 5:
          {
            b3_nlnstiffmass<5>(params,myacc,myvel,mydisp,NULL,NULL,&elevec1,NULL);
            break;
          }
          default:
            dserror("Only Line2, Line3, Line4, and Line5 Elements implemented.");
        }
      }

    /*at the end of an iteration step the geometric configuration has to be updated: the starting point for the
     * next iteration step is the configuration at the end of the current step */
    Qold_ = Qnew_;
    dispthetaold_= dispthetanew_;

//    //the following code block can be used to check quickly whether the nonlinear stiffness matrix is calculated
//    //correctly or not by means of a numerically approximated stiffness matrix
//    //The code block will work for all higher order elements.
//    if(Id() == 0) //limiting the following tests to certain element numbers
//    {
//      //variable to store numerically approximated stiffness matrix
//      Epetra_SerialDenseMatrix stiff_approx;
//      stiff_approx.Shape(6*nnode,6*nnode);
//
//      //relative error of numerically approximated stiffness matrix
//      Epetra_SerialDenseMatrix stiff_relerr;
//      stiff_relerr.Shape(6*nnode,6*nnode);
//
//      //characteristic length for numerical approximation of stiffness
//      double h_rel = 1e-10;
//
//      //flag indicating whether approximation leads to significant relative error
//      int outputflag = 0;
//
//      //calculating strains in new configuration
//      for(int i=0; i<6; i++) //for all dof
//      {
//        for(int k=0; k<nnode; k++)//for all nodes
//        {
//          Epetra_SerialDenseVector force_aux;
//          force_aux.Size(6*nnode);
//
//          //create new displacement and velocity vectors in order to store artificially modified displacements
//          std::vector<double> vel_aux(myvel);
//          std::vector<double> disp_aux(mydisp);
//
//          //modifying displacement artificially (for numerical derivative of internal forces):
//          disp_aux[6*k + i] += h_rel;
//          vel_aux[6*k + i] += h_rel / params.get<double>("delta time",0.01);
//
//          //b3_nlnstiffmass is a templated function. therefore we need to point out the number of nodes in advance
//          switch(nnode)
//          {
//              case 2:
//              {
//                b3_nlnstiffmass<2>(params,myacc,vel_aux,disp_aux,NULL,NULL,&force_aux,NULL);
//                break;
//              }
//              case 3:
//              {
//                b3_nlnstiffmass<3>(params,myacc,vel_aux,disp_aux,NULL,NULL,&force_aux,NULL);
//                break;
//              }
//              case 4:
//              {
//                b3_nlnstiffmass<4>(params,myacc,vel_aux,disp_aux,NULL,NULL,&force_aux,NULL);
//                break;
//              }
//              case 5:
//              {
//                b3_nlnstiffmass<5>(params,myacc,vel_aux,disp_aux,NULL,NULL,&force_aux,NULL);
//                break;
//              }
//              default:
//                dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
//          }
//          //computing derivative d(fint)/du numerically by finite difference
//          for(int u = 0 ; u < 6*nnode ; u++ )
//          {
//            stiff_approx(u,k*6+i)= ( force_aux[u] - elevec1(u) )/ h_rel ;
//            //stiff_approx(u,k*6+i)= ( pow(force_aux[u],2) - pow(elevec1(u),2) )/ (h_rel * (force_aux[u] + elevec1(u) ) );
//          }
//        } //for(int k=0; k<nnode; k++)//for all nodes
//
//      } //for(int i=0; i<3; i++) //for all dof
//
//      for(int line=0; line<6*nnode; line++)
//      {
//        for(int col=0; col<6*nnode; col++)
//        {
//          if (fabs(elemat1(line,col)) > 1.0e-10)
//            stiff_relerr(line,col)= fabs( ( elemat1(line,col) - stiff_approx(line,col) )/ elemat1(line,col) );
//          else if (fabs(stiff_approx(line,col)) < 1.0e-5)
//            stiff_relerr(line,col)=0.0;
//          else
//            stiff_relerr(line,col)=1000.0;
//
//          //suppressing small entries whose effect is only confusing and NaN entires (which arise due to zero entries)
//          if ( fabs( stiff_relerr(line,col) ) < h_rel*500 || isnan( stiff_relerr(line,col))) //isnan = is not a number
//            stiff_relerr(line,col) = 0;
//
//          //if ( stiff_relerr(line,col) > 0)
//            outputflag = 1;
//        } //for(int col=0; col<3*nnode; col++)
//
//      } //for(int line=0; line<3*nnode; line++)
//
//      if(outputflag ==1)
//      {
//
//        std::cout<<"\n\n acutally calculated stiffness matrix\n";
//        for(int line=0; line<6*nnode; line++)
//        {
//          for(int col=0; col<6*nnode; col++)
//          {
//            if(isnan(elemat1(line,col)))
//              std::cout<<"     nan   ";
//            else if(elemat1(line,col) == 0)
//              std::cout<<"     0     ";
//            else if(elemat1(line,col) >= 0)
//              std::cout<<"  "<< std::scientific << std::setprecision(3)<<elemat1(line,col);
//            else
//              std::cout<<" "<< std::scientific << std::setprecision(3)<<elemat1(line,col);
//          }
//          std::cout<<"\n";
//        }
//
//        std::cout<<"\n\n approximated stiffness matrix\n";
//        for(int line=0; line<6*nnode; line++)
//        {
//          for(int col=0; col<6*nnode; col++)
//          {
//            if(isnan(stiff_approx(line,col)))
//              std::cout<<"     nan   ";
//            else if(stiff_approx(line,col) == 0)
//              std::cout<<"     0     ";
//            else if(stiff_approx(line,col) >= 0)
//              std::cout<<"  "<< std::scientific << std::setprecision(3)<<stiff_approx(line,col);
//            else
//              std::cout<<" "<< std::scientific << std::setprecision(3)<<stiff_approx(line,col);
//          }
//          std::cout<<"\n";
//        }
//
//        std::cout<<"\n\n rel error stiffness matrix\n";
//        for(int line=0; line<6*nnode; line++)
//        {
//          for(int col=0; col<6*nnode; col++)
//          {
//            if(isnan(stiff_relerr(line,col)))
//              std::cout<<"     nan   ";
//            else if(stiff_relerr(line,col) == 0)
//              std::cout<<"     0     ";
//            else if(stiff_relerr(line,col) >= 0)
//              std::cout<<"  "<< std::scientific << std::setprecision(3)<<stiff_relerr(line,col);
//            else
//              std::cout<<" "<< std::scientific << std::setprecision(3)<<stiff_relerr(line,col);
//          }
//          std::cout<<"\n";
//        }
//      }
//
//    } //end of section in which numerical approximation for stiffness matrix is computed

    }
    break;
    case calc_struct_update_istep:
    {
      /*the action calc_struct_update_istep is called in the very end of a time step when the new dynamic
       * equilibrium has finally been found; this is the point where the variable representing the geometric
       * status of the beam at the end of the time step has to be stored*/
      Qconv_ = Qnew_;
      Qconvmass_ = Qnewmass_;
      wconvmass_ = wnewmass_;
      aconvmass_ = anewmass_;
      amodconvmass_ = amodnewmass_;
      rttconvmass_ = rttnewmass_;
      rttmodconvmass_ = rttmodnewmass_;
      rtconvmass_ = rtnewmass_;
      dispconvmass_ = dispnewmass_;
      dispthetaconv_ = dispthetanew_;
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
      dispthetaold_ = dispthetaconv_;
      Qnew_ = Qconv_;
      dispthetanew_ = dispthetaconv_;
      Qnewmass_=Qconvmass_;
      wnewmass_=wconvmass_;
      anewmass_=aconvmass_;
      amodnewmass_=amodconvmass_;
      rttnewmass_=rttconvmass_;
      rttmodnewmass_=rttmodconvmass_;
      rtnewmass_=rtconvmass_;
      dispnewmass_=dispconvmass_;
      dispthetanew_=dispthetaconv_;
    }
    break;
    case calc_struct_stress:
      dserror("No stress output implemented for beam3r elements");
    break;
    default:
      dserror("Unknown type of action for Beam3r %d", act);
    break;
  }
  return 0;

}

/*-----------------------------------------------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)                                       cyron 03/08|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3r::EvaluateNeumann(Teuchos::ParameterList& params,
                                        DRT::Discretization& discretization,
                                        DRT::Condition& condition,
                                        std::vector<int>& lm,
                                        Epetra_SerialDenseVector& elevec1,
                                        Epetra_SerialDenseMatrix* elemat1)
{
  // get element displacements
  Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
  if (disp==Teuchos::null) dserror("Cannot get state vector 'displacement'");
  std::vector<double> mydisp(lm.size());
  DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // no. of nodes on this element; the following line is only valid for elements with constant number of
  // degrees of freedom per node
  const int numdf = 6;
  const DiscretizationType distype = this->Shape();

  // find out whether we will use a time curve and get the factor
  const std::vector<int>* curve  = condition.Get<std::vector<int> >("curve");
  // amplitude of load curve at current time called
  std::vector<double> curvefac(numdf,1.0);

  for (int i=0; i<numdf; ++i)
  {
    int curvenum = -1;
    // number of the load curve related with a specific line Neumann condition called
    if (curve) curvenum = (*curve)[i];

    if (curvenum>=0 && usetime)
      curvefac[i] = DRT::Problem::Instance()->Curve(curvenum).f(time);
  }

  // gaussian points
  const DRT::UTILS::IntegrationPoints1D intpoints(MyGaussRule(NumNode(),gaussunderintegration));

  //declaration of variable in order to store shape function
  Epetra_SerialDenseVector funct(NumNode());

  // get values and switches from the condition

  // onoff is related to the first 6 flags of a line Neumann condition in the input file;
  // value 1 for flag i says that condition is active for i-th degree of freedom
  const std::vector<int>* onoff = condition.Get<std::vector<int> >("onoff");
  // val is related to the 6 "val" fields after the onoff flags of the Neumann condition
  // in the input file; val gives the values of the force as a multiple of the prescribed load curve
  const std::vector<double>* val = condition.Get<std::vector<double> >("val");
  // funct is related to the 6 "funct" fields after the val field of the Neumann condition
  // in the input file; funct gives the number of the function defined in the section FUNCT
  const std::vector<int>* functions = condition.Get<std::vector<int> >("funct");

  //integration loops
  for (int numgp=0; numgp<intpoints.nquad; ++numgp)
  {
    //integration points in parameter space and weights
    const double xi = intpoints.qxg[numgp][0];
    const double wgt = intpoints.qwgt[numgp];
    //evaluation of shape funcitons at Gauss points
    DRT::UTILS::shape_function_1D(funct,xi,distype);

    //position vector at the gauss point at reference configuration needed for function evaluation
    std::vector<double> X_ref(3,0.0);
    //calculate coordinates of corresponding Gauss point in reference configuration
    for (int node=0;node<NumNode();node++)
    {
      for (int dof=0;dof<3;dof++)
      {
        X_ref[dof]+=funct[node]*Nodes()[node]->X()[dof];
      }
    }

    double fac=0;
    fac = wgt * jacobi_[numgp];

    // load vector ar
    double ar[numdf];

    // loop the dofs of a node
    for (int dof=0; dof<numdf; ++dof)
      ar[dof] = fac * (*onoff)[dof]*(*val)[dof]*curvefac[dof];
    double functionfac = 1.0;
    int functnum = -1;

    //sum up load components
    for (int dof=0; dof<6; ++dof)
    {
      if (functions) functnum = (*functions)[dof];
      else functnum = -1;

      if (functnum>0)
      {
        // evaluate function at the position of the current GP
        functionfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(dof, &X_ref[0], time, NULL);
      }
      else functionfac = 1.0;

      for (int node=0; node<NumNode(); ++node)
      {
        elevec1[node*numdf+dof] += funct[node] *ar[dof] *functionfac;
      }
    }
  } // for (int numgp=0; numgp<intpoints.nquad; ++numgp)

  return 0;
}


/*----------------------------------------------------------------------------------------------------------------------*
 |compute convected stresses from convected strains and returning also constitutive matrix between both according to    |
 |Jelenic 1999, section 2.4)                                                                                 cyron 04/10|
 *----------------------------------------------------------------------------------------------------------------------*/
template <typename T>
inline void DRT::ELEMENTS::Beam3r::strainstress(const LINALG::TMatrix<T,3,1>& Gamma,
                                                 const LINALG::TMatrix<T,3,1>& K,
                                                 LINALG::TMatrix<T,3,1>& stressN,
                                                 LINALG::TMatrix<T,3,3>& CN,
                                                 LINALG::TMatrix<T,3,1>& stressM,
                                                 LINALG::TMatrix<T,3,3>& CM)
{
  //first of all we get the material law
  Teuchos::RCP<const MAT::Material> currmat = Material();
  double ym = 0;
  double sm = 0;

  //assignment of material parameters; only St.Venant material is accepted for this beam
  switch(currmat->MaterialType())
  {
    case INPAR::MAT::m_stvenant:// only linear elastic material supported
    {
      const MAT::StVenantKirchhoff* actmat = static_cast<const MAT::StVenantKirchhoff*>(currmat.get());
      ym = actmat->Youngs();
      sm = ym / (2*(1 + actmat->PoissonRatio()));
    }
    break;
    default:
      dserror("unknown or improper type of material law");
    break;
  }

  //defining material constitutive matrix CN between Gamma and N according to Jelenic 1999, section 2.4
  CN.PutScalar(0);
  CN(0,0) = ym*crosssec_;
  CN(1,1) = sm*crosssecshear_;
  CN(2,2) = sm*crosssecshear_;

  //defining material constitutive matrix CM between curvature and moment according to Jelenic 1999, section 2.4
  CM.PutScalar(0);
  CM(0,0) = sm*Irr_;
  CM(1,1) = ym*Iyy_;
  CM(2,2) = ym*Izz_;

  //computing stresses by multiplying strains with respective constitutive matrix
  stressN.Multiply(CN,Gamma);
  stressM.Multiply(CM,K);

  return;
} // DRT::ELEMENTS::Beam3r::straintostress

/*----------------------------------------------------------------------------------------------------------------------*
 |push forward stresses and constitutive matrix to their spatial counterparts by rotation matrix Lambda according to    |
 |Romero 2004, eq. (3.10)                                                                                    cyron 04/10|
 *----------------------------------------------------------------------------------------------------------------------*/
template <typename T>
inline void DRT::ELEMENTS::Beam3r::pushforward(const LINALG::TMatrix<T,3,3>& Lambda,
                                                const LINALG::TMatrix<T,3,1>& stressN,
                                                const LINALG::TMatrix<T,3,3>& CN,
                                                const LINALG::TMatrix<T,3,1>& stressM ,
                                                const LINALG::TMatrix<T,3,3>& CM,
                                                LINALG::TMatrix<T,3,1>& stressn,
                                                LINALG::TMatrix<T,3,3>& cn,
                                                LINALG::TMatrix<T,3,1>& stressm,
                                                LINALG::TMatrix<T,3,3>& cm)
{
  //introduce auxiliary variable for pushforward of rotational matrices
  LINALG::TMatrix<T,3,3> temp;

  //push forward translational stresses
  stressn.Multiply(Lambda,stressN);

  //pushforward translational constitutive matrix CN to matrix cn according to Jelenic 1999, paragraph following to (2.22) on page 148
  temp.Multiply(Lambda,CN);
  cn.MultiplyNT(temp,Lambda);

  //push forward rotational stresses stresses
  stressm.Multiply(Lambda,stressM);

  //pushforward translational constitutive matrix CM to matrix cm according to Jelenic 1999, paragraph following to (2.22) on page 148
  temp.Multiply(Lambda,CM);
  cm.MultiplyNT(temp,Lambda);

   return;
} // DRT::ELEMENTS::Beam3r::pushforward

/*------------------------------------------------------------------------------------------------------------*
 | calculation of elastic energy (private)                                                        cyron 12/10|
 *-----------------------------------------------------------------------------------------------------------*/
template<unsigned int nnode>
void DRT::ELEMENTS::Beam3r::b3_energy( Teuchos::ParameterList& params,
                                        std::vector<double>& disp,
                                        Epetra_SerialDenseVector* intenergy)
{
  //initialize energies (only one kind of energy computed here
  intenergy->Scale(0.0);
  Eint_=0.0;

  bool calcenergy = false;
  if(params.isParameter("energyoftype")==false) calcenergy = true;
  else if(params.get<std::string>("energyoftype")=="beam3r") calcenergy =true;

  if(calcenergy)
  {
    //const double t_tot = Teuchos::Time::wallTime();

    // Note: index i refers to the nnode nodes
    //vector whose numgp-th element is a 1xnnode-matrix with all shape functions evaluated at the numgp-th GP
    std::vector<LINALG::Matrix<1,nnode> > I_i(nnode-1);

    //vector whose numgp-th element is a 1xnnode-matrix with the derivatives of all shape functions evaluated at the numgp-th GP
    std::vector<LINALG::Matrix<1,nnode> > I_i_xi(nnode-1);

    //current nodal DOFs in total Lagrangian style, i.e. initial values + displacements
    LINALG::Matrix<6*nnode,1> disp_totlag(true);

    //derivative of initial curve in physical space with respect to curve parameter s
    LINALG::Matrix<3,1> r_s;

    // dummy 3D vector
    LINALG::Matrix<3,1> dummy(true);

    // reference quaternion Q_r corresponding to reference triad Lambda_r
    LINALG::Matrix<4,1> Q_r(true);

    // rotation angles between nodal triads and reference triad according to (3.8), Jelenic 1999
    // use type TMatrix here to enable use of functions Calc_Psi_l and Calc_Psi_l (defined for std::vector<LINALG::TMatrix>)
    std::vector<LINALG::TMatrix<double,3,1> > Psi_li(nnode);

    //interpolated local relative rotation \Psi^l at a certain Gauss point according to (3.11), Jelenic 1999
    LINALG::Matrix<3,1> Psi_l;
    /*derivative of interpolated local relative rotation \Psi^l with respect to arc-length parameter
     * at a certain Gauss point according to (3.11), Jelenic 1999*/
    LINALG::Matrix<3,1> Psi_l_s;
    // triad at GP
    LINALG::Matrix<3,3> Lambda;

    //3D vector related to spin matrix \hat{\kappa} from (2.1), Jelenic 1999
    LINALG::Matrix<3,1> K;
    //3D vector of material axial and shear strains from (2.1), Jelenic 1999
    LINALG::Matrix<3,1> Gamma;

    //material stresses N and M and constitutive matrices C_N and C_M according to section 2.4, Jelenic 1999
    LINALG::Matrix<3,1> stressN;
    LINALG::Matrix<3,1> stressM;
    LINALG::Matrix<3,3> CN;
    LINALG::Matrix<3,3> CM;


    /*first displacement vector is modified for proper element evaluation in case of periodic boundary conditions; in case that
     *no periodic boundary conditions are to be applied the following code line may be ignored or deleted*/
    if(params.isParameter("PERIODLENGTH"))
      NodeShift<nnode,3>(params,disp);

    // get current values of translational nodal DOFs in total Lagrangean manner (initial value + disp)
    for (unsigned int dof=0; dof<3; ++dof)
      for (unsigned int node=0; node<nnode; ++node)
        disp_totlag(6*node+dof) = Nodes()[node]->X()[dof] + disp[6*node+dof];

    // rotational DOFs need different handling and quaternions at nodes Qnew_ should already be updated at this point
    // update depends on choice of multiplicative or additive (FAD) increments and is done in b3_nlnstiffmass, which therefore always needs to be called before b3energy!
    // see comments in b3_nlnstiffmass for further discussion and update routine of Qnew_

    // compute reference triad Lambda_r according to (3.9), Jelenic 1999
    // Qnew_ has already been filled with initial values from input file in ReadElement()
    CalcRefQuaternion<double>(Qnew_[nodeI_],Qnew_[nodeJ_],Q_r,dummy);

    for(unsigned int node=0; node<nnode; node++)
    {
      // compute nodal local rotations according to (3.8), Jelenic 1999
      CalcPsi_li<double>(Qnew_[node],Q_r,Psi_li[node]);
    }

    //Get the applied integrationpoints for underintegration
    DRT::UTILS::IntegrationPoints1D gausspoints_underint(MyGaussRule(nnode,gaussunderintegration));

    // evaluate all shape functions and derivatives with respect to element parameter xi at all specified Gauss points
    EvaluateShapeFunctions<nnode>(gausspoints_underint,I_i,I_i_xi,this->Shape());

    //Loop through all GPs for under-integration and calculate jacobi determinants at the GPs
    for(int numgp=0; numgp < gausspoints_underint.nquad; numgp++)
    {
      //weight of GP in parameter space
      const double wgt = gausspoints_underint.qwgt[numgp];

      Calc_r_s<nnode,double>(disp_totlag, I_i_xi[numgp], jacobi_[numgp], r_s);
      Calc_Psi_l<nnode,double>(Psi_li, I_i[numgp], Psi_l);
      Calc_Psi_l_s<nnode,double>(Psi_li, I_i_xi[numgp], jacobi_[numgp], Psi_l_s);
      Calc_Lambda<double>(Psi_l,Q_r,Lambda);

      //compute material strains gamma and kappa according to Jelenic 1999, eq. (2.12)
      computestrain<double>(Psi_l,Psi_l_s,r_s,Lambda,Gammaref_[numgp],Kref_[numgp],Gamma,K);

      //compute material stress vector from strain vector according to Jelenic 1999, page 147, section 2.4
      strainstress<double>(Gamma,K,stressN,CN,stressM,CM);

      //adding elastic energy at this Gauss point
      if(intenergy->M()==1)
      {
        for(int i=0; i<3; i++)
        {
          (*intenergy)(0) += 0.5*Gamma(i)*stressN(i)*wgt*jacobi_[numgp];
          (*intenergy)(0) += 0.5*K(i)*stressM(i)*wgt*jacobi_[numgp];
        }
        Eint_=(*intenergy)(0);
      }
      else if(intenergy->M()==6)
      {
        for(int i=0; i<3; i++)
        {
          (*intenergy)(i) += 0.5*Gamma(i)*stressN(i)*wgt*jacobi_[numgp];
          (*intenergy)(i+3) += 0.5*K(i)*stressM(i)*wgt*jacobi_[numgp];
        }
      }
      else
        dserror("energy vector of invalid size %i, expected row dimension 1 (total elastic energy of element) or 6 (elastic energy of 3 translational and 3 rotational spatial dimensions)!", intenergy->M());
    }

  }
  return;

} // DRT::ELEMENTS::Beam3r::b3_energy


/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   cyron 01/08|
 *-----------------------------------------------------------------------------------------------------------*/
template<unsigned int nnode>
void DRT::ELEMENTS::Beam3r::b3_nlnstiffmass( Teuchos::ParameterList& params,
                                            std::vector<double>&      acc,
                                            std::vector<double>&      vel,
                                            std::vector<double>&      disp,
                                            Epetra_SerialDenseMatrix* stiffmatrix,
                                            Epetra_SerialDenseMatrix* massmatrix,
                                            Epetra_SerialDenseVector* force,
                                            Epetra_SerialDenseVector* inertia_force)
{
  //internal force vector
  LINALG::TMatrix<FADordouble, 6*nnode, 1> f_int(true);

  // Note: index i refers to the nnode nodes
  //vector whose numgp-th element is a 1xnnode-matrix with all Lagrange polynomial basis functions evaluated at the numgp-th Gauss point
  std::vector<LINALG::Matrix<1,nnode> > I_i(nnode-1);

  //vector whose numgp-th element is a 1xnnode-matrix with the derivatives of all Lagrange polynomial basis functions evaluated at the numgp-th Gauss point
  std::vector<LINALG::Matrix<1,nnode> > I_i_xi(nnode-1);

  //vector whose numgp-th element is a 1xnnode-matrix with all Lagrange polynomial basis functions evaluated at the numgp-th Gauss point
  std::vector<LINALG::Matrix<1,nnode> > Imass_i(nnode);

  //vector with nnode elements, who represent the 3x3-matrix-shaped interpolation function \tilde{I}^nnode at a certain Gauss point for elasticity according to (3.18), Jelenic 1999
  std::vector<LINALG::TMatrix<double,3,3> > Itilde(nnode);

  //vector with nnode elements, who represent the 3x3-matrix-shaped interpolation function \tilde{I'}^nnode at a certain Gauss point for elasticity according to (3.19), Jelenic 1999
  std::vector<LINALG::TMatrix<double,3,3> > Itildeprime(nnode);

  //vector with nnode elements, who represent the 3x3-matrix-shaped interpolation function \tilde{I}^nnode at a certain  Gauss point for mass matrix according to (3.19), Jelenic 1999
  std::vector<LINALG::TMatrix<double,3,3> > Itildemass(nnode);

  //! quaternions of all nodal triads
  std::vector<LINALG::TMatrix<FADordouble,4,1> > Q_i(nnode);

  //reference triad Lambda_r and corresponding quaternion Q_r
  LINALG::TMatrix<FADordouble,3,3> Lambda_r(true);
  LINALG::TMatrix<FADordouble,4,1> Q_r(true);

  //angle of relative rotation between node I and J according to (3.10), Jelenic 1999
  LINALG::TMatrix<FADordouble,3,1> Phi_IJ(true);

  //rotation angles between nodal triads and reference triad according to (3.8), Jelenic 1999
  std::vector<LINALG::TMatrix<FADordouble,3,1> > Psi_li(nnode);

  //rotational displacement at a certain node between this and last iteration step
  LINALG::TMatrix<FADordouble,3,1>  deltatheta;
  //rotational displacement at a certain node between this and last iteration step in quaternion form
  LINALG::TMatrix<FADordouble,4,1>  deltaQ;

  // derivation of beam centerline with respect to arc-length parameter: r'(x) from (2.12), Jelenic 1999
  LINALG::TMatrix<FADordouble,3,1> r_s;
  //spin matrix related to vector r_s
  LINALG::TMatrix<FADordouble,3,3> r_s_hat;
  //interpolated local relative rotation \Psi^l at a certain Gauss point according to (3.11), Jelenic 1999
  LINALG::TMatrix<FADordouble,3,1> Psi_l;
  /*derivative of interpolated local relative rotation \Psi^l with respect to arc-length parameter
   * at a certain Gauss point according to (3.11), Jelenic 1999*/
  LINALG::TMatrix<FADordouble,3,1> Psi_l_s;
  // triad at GP
  LINALG::TMatrix<FADordouble,3,3> Lambda;

  //3D vector related to spin matrix \hat{\kappa} from (2.1), Jelenic 1999
  LINALG::TMatrix<FADordouble,3,1> K;
  //3D vector of material axial and shear strains from (2.1), Jelenic 1999
  LINALG::TMatrix<FADordouble,3,1> Gamma;

  //convected stresses N and M and constitutive matrices C_N and C_M according to section 2.4, Jelenic 1999
  LINALG::TMatrix<FADordouble,3,1> stressN;
  LINALG::TMatrix<FADordouble,3,1> stressM;
  LINALG::TMatrix<FADordouble,3,3> CN;
  LINALG::TMatrix<FADordouble,3,3> CM;

  //spatial stresses n and m according to (3.10), Romero 2004 and spatial constitutive matrices c_n and c_m according to page 148, Jelenic 1999
  LINALG::TMatrix<FADordouble,3,1> stressn;
  LINALG::TMatrix<FADordouble,3,1> stressm;
  LINALG::TMatrix<FADordouble,3,3> cn;
  LINALG::TMatrix<FADordouble,3,3> cm;

  /*first displacement vector is modified for proper element evaluation in case of periodic boundary conditions; in case that
   *no periodic boundary conditions are to be applied the following code line may be ignored or deleted*/
  if(params.isParameter("PERIODLENGTH"))
    NodeShift<nnode,3>(params,disp);

  //current nodal DOFs in total Lagrangian style, i.e. initial values + displacements
  LINALG::TMatrix<FADordouble,6*nnode,1> disp_totlag(true);


  // get current values of translational nodal DOFs in total Lagrangean manner (initial value + disp)
  // rotational DOFs need different handling, depending on whether FAD is used or not (see comment below)
  for (unsigned int dof=0; dof<3; ++dof)
    for (unsigned int node=0; node<nnode; ++node)
      disp_totlag(6*node+dof) = Nodes()[node]->X()[dof] + disp[6*node+dof];


  #ifndef BEAM3RAUTOMATICDIFF
  //Compute current nodal triads
  for (unsigned int node=0; node<nnode; ++node)
  {
    /*rotation increment relative to configuration in last iteration step is difference between current rotation
     *entry in displacement vector minus rotation entry in displacement vector in last iteration step*/

    //This shift is necessary, since our beam formulation and the corresponding linearization is based on multiplicative
    //increments, while in BACI (Newtonfull()) the displacement of the last iteration and the rotation increment
    //of the current iteration are added in an additive manner. This step is reversed in the next two lines in order
    //to recover the multiplicative rotation increment between the last and the current Newton step. This also
    //means that dispthetanew_ has no physical meaning since it is the additive sum of non-additive rotation increments(meier, 03.2014)
    for(int i=0; i<3; i++)
    {
      dispthetanew_[node](i) = disp[6*node+3+i];
    }

    deltatheta  = dispthetanew_[node];
    deltatheta -= dispthetaold_[node];

    //compute quaternion from rotation angle relative to last configuration
    LARGEROTATIONS::angletoquaternion(deltatheta,deltaQ);

    //multiply relative rotation with rotation in last configuration to get rotation in new configuration
    LARGEROTATIONS::quaternionproduct(Qold_[node],deltaQ,Qnew_[node]);

    //renormalize quaternion to keep its absolute value one even in case of long simulations and intricate calculations
    Qnew_[node].Scale(1/Qnew_[node].Norm2());

    // copy quaternions of nodal triads to TMatrix FADordouble
    for (unsigned int node=0; node<nnode; ++node)
      for (int i=0; i<4; ++i)
        Q_i[node](i)= Qnew_[node](i);
  }
  #else


  // set differentiation variables for FAD: translational DOFs
  for (unsigned int node=0; node< nnode; ++node)
    for (int i=0; i<3; ++i)
      disp_totlag(6*node+i).diff(6*node+i,6*nnode);     // TODO remove redundancy !

  #ifndef MULTIPLICATIVEUPDATES
  // note: when using FAD, we do not need to worry about multiplicative rotational increments for rotational DOFs
  // FAD will compute additive increments ("normal" differentiation), just as FullNewton expects, so everything is fine
  // => rotational DOFs are treated just like translational DOFs (initial value + disp)
  for (unsigned int node=0; node<nnode; ++node)
    for (unsigned int dof=0; dof<3; ++dof)
    {
      disp_totlag(6*node+3+dof) = theta0_[node](dof) + disp[6*node+3+dof];
      // set differentiation variables for FAD: rotational DOFs
      disp_totlag(6*node+3+dof).diff(6*node+3+dof,6*nnode);
    }

  // set differentiation variables for FAD
  for (unsigned int i=0; i< 6*nnode; ++i)
      disp_totlag(i).diff(i,6*nnode);

  LINALG::TMatrix<FAD,3,1> theta_totlag(true);

  for (unsigned int node=0; node<nnode; ++node)
  {
    for(int dof=0; dof<3; dof++)
    {
      theta_totlag(dof) = disp_totlag(6*node+3+dof);
    }

    //compute nodal quaternions
    LARGEROTATIONS::angletoquaternion(theta_totlag,Q_i[node]);

    //renormalize quaternion to keep its absolute value one even in case of long simulations and intricate calculations.
    Q_i[node].Scale(1.0/FADUTILS::VectorNorm<4>(Q_i[node]));
  }

  #else


  // compute nodal quaternions based on multiplicative increments of rotational DOFs
  for (unsigned int node=0; node<nnode; ++node)
  {
    for(int i=0; i<3; i++)
    {
      dispthetanew_[node](i) = disp[6*node+3+i];

      deltatheta(i)  = dispthetanew_[node](i);
      deltatheta(i) -= dispthetaold_[node](i);
    }

    // set differentiation variables for FAD: rotational DOFs
    for(int i=0; i<3; i++)
    {
      deltatheta(i).diff(6*node+3+i,6*nnode);
    }

    LINALG::TMatrix<FAD,4,1> Q_i_old(true);
    for (int i=0; i<4; ++i)
      Q_i_old(i) = Qold_[node](i);

    //compute quaternion from rotation angle relative to last configuration
    LARGEROTATIONS::angletoquaternion(deltatheta,deltaQ);

    //multiply relative rotation with rotation in last configuration to get rotation in new configuration
    LARGEROTATIONS::quaternionproduct(Q_i_old,deltaQ,Q_i[node]);

    //renormalize quaternion to keep its absolute value one even in case of long simulations and intricate calculations.
    Q_i[node].Scale(1.0/FADUTILS::VectorNorm<4>(Q_i[node]));
  }

  #endif //#ifndef MULTIPLICATIVEUPDATE

  // copy quaternions of nodal triads to class variable Qnew_
  for (unsigned int node=0; node<nnode; ++node)
    for (int i=0; i<4; ++i)
      Qnew_[node](i)= FADUTILS::CastToDouble(Q_i[node](i));

  #endif //#ifndef BEAM3RAUTOMATICDIFF

  // compute reference triad Lambda_r according to (3.9), Jelenic 1999
  CalcRefQuaternion<FADordouble>(Q_i[nodeI_],Q_i[nodeJ_],Q_r,Phi_IJ);
  LARGEROTATIONS::quaterniontotriad(Q_r,Lambda_r);

  // integration points for elasticity (underintegration)
  DRT::UTILS::IntegrationPoints1D gausspoints(MyGaussRule(nnode,gaussunderintegration));

  // evaluate all shape functions and derivatives with respect to element parameter xi at all specified Gauss points
  EvaluateShapeFunctions<nnode>(gausspoints,I_i,I_i_xi,this->Shape());

  // compute nodal local rotations according to (3.8), Jelenic 1999
  // this is done individually for each node in order to avoid function argument std::vector<LINALG::TMatrix<...> >
  // a function with this argument type cannot be called with an argument of type std::vector<LINALG::Matrix<...> > (needed e.g. in SetUpReferenceGeometry)
  for (unsigned int node=0; node<nnode; ++node)
    CalcPsi_li<FADordouble>(Q_i[node],Q_r,Psi_li[node]);

  // reset norm of maximal bending curvature
  K_max_=0.0;

  //Loop through all GP and calculate their contribution to the forcevector and stiffnessmatrix
  for(int numgp=0; numgp < gausspoints.nquad; numgp++)
  {
    //weight of GP in parameter space
    const double wgt = gausspoints.qwgt[numgp];

    Calc_r_s<nnode,FADordouble>(disp_totlag, I_i_xi[numgp], jacobi_[numgp], r_s);
    Calc_Psi_l<nnode,FADordouble>(Psi_li, I_i[numgp], Psi_l);
    Calc_Psi_l_s<nnode,FADordouble>(Psi_li, I_i_xi[numgp], jacobi_[numgp], Psi_l_s);
    Calc_Lambda<FADordouble>(Psi_l,Q_r,Lambda);

    //compute spin matrix related to vector rprime for later use
    LARGEROTATIONS::computespin<FADordouble>(r_s_hat,r_s);

    //compute material strains gamma and kappa according to Jelenic 1999, eq. (2.12)
    computestrain<FADordouble>(Psi_l,Psi_l_s,r_s,Lambda,Gammaref_[numgp],Kref_[numgp],Gamma,K);

    //determine norm of maximal bending curvature at this GP and store in class variable if needed
    double Kmax = std::sqrt(FADUTILS::CastToDouble(K(1))*FADUTILS::CastToDouble(K(1)) + FADUTILS::CastToDouble(K(2))*FADUTILS::CastToDouble(K(2)) );
    if(Kmax > K_max_)
      K_max_=Kmax;

    //compute material stress vector from strain vector according to Jelenic 1999, page 147, section 2.4
    strainstress<FADordouble>(Gamma,K,stressN,CN,stressM,CM);

    /*compute spatial stresses and constitutive matrices from convected ones according to Jelenic 1999, page 148, paragraph
     *between (2.22) and (2.23) and Romero 2004, (3.10)*/
    pushforward<FADordouble>(Lambda,stressN,CN,stressM,CM,stressn,cn,stressm,cm);

    /*computation of internal forces according to Jelenic 1999, eq. (4.3); computation split up with respect
     *to single blocks of matrix in eq. (4.3); note that Jacobi determinant in diagonal blocks cancels out
     *in implementation, whereas for the lower left block we have to multiply the weight by the jacobi
     *determinant*/
    for (unsigned int node=0; node<nnode; ++node)
    {
      /*upper left block (note: jacobi determinant cancels out as deriv is derivative with respect to
       *parameter in Gauss integration interval and I^{i'} in Jelenic 1999 is derivative with respect to
       *curve length in reference configuration*/
      for (int i=0; i<3; ++i)
        f_int(6*node+i) += I_i_xi[numgp](node)*stressn(i)*wgt;

      //lower left block
      for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
          f_int(6*node+3+i) -= r_s_hat(i,j)*stressn(j)*I_i[numgp](node)*wgt*jacobi_[numgp];

      /*lower right block (note: jacobi determinant cancels out as Iprime is derivative with respect to
       *parameter in Gauss integration interval and I^{i'} in Jelenic 1999 is derivative with respect to
       *curve length in reference configuration*/
      for (int j=0; j<3; ++j)
        f_int(6*node+3+j) += I_i_xi[numgp](node)*stressm(j)*wgt;
    }

    if(force!=NULL)
    {
      for (unsigned int i=0; i< 6*nnode; i++)
      {
        (*force)(i)= FADUTILS::CastToDouble(f_int(i));
      }
    }

    if (stiffmatrix != NULL)
    {
      #ifdef BEAM3RAUTOMATICDIFF
      //Calculating stiffness matrix with FAD
      for (unsigned int i = 0; i < 6*nnode; i++)
      {
        for (unsigned int j = 0; j < 6*nnode; j++)
        {
          (*stiffmatrix)(i,j)=f_int(i).dx(j);
        }
      }

      #ifdef MULTIPLICATIVEUPDATES

      // we need to transform the stiffmatrix because its entries are derivatives with respect to additive rotational increments
      // we want a stiffmatrix containing derivatives with respect to multiplicative rotational increments
      // therefore apply a trafo matrix to all those 3x3 blocks in stiffmatrix which correspond to derivation with respect to rotational DOFs
      // the trafo matrix is simply the T-Matrix (see Jelenic1999, (2.4)): \Delta_{mult} \vec \theta_{inode} = \mat T(\vec \theta_{inode} * \Delta_{addit} \vec \theta_{inode}
      LINALG::TMatrix<FAD,3,3> tempmat(true);
      LINALG::TMatrix<FAD,3,3> newstiffmat(true);
      LINALG::TMatrix<FAD,3,3> Tmat(true);
      LINALG::TMatrix<FAD,3,1> deltatheta_i(true);
      for (unsigned int jnode = 0; jnode < nnode; jnode++)
      {
        for (int i=0; i<3; ++i)
          deltatheta_i(i) = dispthetanew_[jnode](i) - dispthetaold_[jnode](i);

        Tmat = LARGEROTATIONS::Tmatrix(deltatheta_i);

        for (unsigned int inode = 0; inode < nnode; inode++)
        {
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              tempmat(i,j) = (*stiffmatrix)(6*inode+i,6*jnode+3+j);

          newstiffmat.Clear();
          newstiffmat.MultiplyNN(tempmat,Tmat);

          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(6*inode+i,6*jnode+3+j) = newstiffmat(i,j).val();

          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              tempmat(i,j) = (*stiffmatrix)(6*inode+3+i,6*jnode+3+j);

          newstiffmat.Clear();
          newstiffmat.MultiplyNN(tempmat,Tmat);

          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(6*inode+3+i,6*jnode+3+j) = newstiffmat(i,j).val();
        }
      }
      #endif

      #else //#ifdef BEAM3RAUTOMATICDIFF

      computeItilde<nnode>(Psi_l,Itilde,Phi_IJ,Lambda_r,Psi_li,I_i[numgp]);
      computeItildeprime<nnode,double>(Psi_l,Psi_l_s,Itildeprime,Phi_IJ,Lambda_r,Psi_li,I_i[numgp],I_i_xi[numgp],jacobi_[numgp]);

      /*computation of stiffness matrix according to Jelenic 1999, eq. (4.7); computation split up with respect
      *to single blocks of matrix in eq. (4.3). Here, the jacobi determinant jacobi_[numgp] stemming from ds and
      *to the inverse determinant 1/jacobi_[numgp] stemming from I^{i'} have to be considered (these two cancel
      *to out if both appear), while the jacobi determinant of Itildeprime and rprime have already been
      *to considered in the corresponding calculation of these quantities.*/

      //auxiliary variables for storing intermediate matrices in computation of entries of stiffness matrix
      LINALG::Matrix<3,3> auxmatrix1;
      LINALG::Matrix<3,3> auxmatrix2;

      for(unsigned int nodei = 0; nodei < nnode; nodei++)
        for(unsigned int nodej = 0; nodej < nnode; nodej++)
        {
          //upper left block
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(6*nodei+i,6*nodej+j) += I_i_xi[numgp](nodei)*I_i_xi[numgp](nodej)*cn(i,j)*wgt/jacobi_[numgp];

          //upper right block
          auxmatrix2.Multiply(cn,r_s_hat);
          LARGEROTATIONS::computespin(auxmatrix1,stressn);
          auxmatrix2 -= auxmatrix1;
          auxmatrix2.Scale(I_i_xi[numgp](nodei));
          auxmatrix1.Multiply(auxmatrix2,Itilde[nodej]);
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(6*nodei+i,6*nodej+3+j) += auxmatrix1(i,j)*wgt;

          //lower left block; note: error in eq. (4.7), Jelenic 1999: the first factor should be I^i instead of I^j
          auxmatrix2.Multiply(r_s_hat,cn);
          LARGEROTATIONS::computespin(auxmatrix1,stressn);
          auxmatrix1 -= auxmatrix2;
          auxmatrix1.Scale(I_i[numgp](nodei)*I_i_xi[numgp](nodej));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(6*nodei+3+i,6*nodej+j) += auxmatrix1(i,j)*wgt;

          //lower right block
          //first summand
          auxmatrix1.Multiply(cm,Itildeprime[nodej]);
          auxmatrix1.Scale(I_i_xi[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(6*nodei+3+i,6*nodej+3+j) += auxmatrix1(i,j)*wgt;

          //second summand
          LARGEROTATIONS::computespin(auxmatrix2,stressm);
          auxmatrix1.Multiply(auxmatrix2,Itilde[nodej]);
          auxmatrix1.Scale(I_i_xi[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(6*nodei+3+i,6*nodej+3+j) -= auxmatrix1(i,j)*wgt;

          //third summand; note: error in eq. (4.7), Jelenic 1999: the first summand in the parantheses should be \hat{\Lambda N} instead of \Lambda N
          LARGEROTATIONS::computespin(auxmatrix1,stressn);
          auxmatrix2.Multiply(cn,r_s_hat);
          auxmatrix1 -= auxmatrix2;
          auxmatrix2.Multiply(auxmatrix1,Itilde[nodej]);
          auxmatrix1.Multiply(r_s_hat,auxmatrix2);
          auxmatrix1.Scale(I_i[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(6*nodei+3+i,6*nodej+3+j) += auxmatrix1(i,j)*jacobi_[numgp]*wgt;

        }
      #endif

    }//if (stiffmatrix != NULL)
  }//for(int numgp=0; numgp < gausspoints.nquad; numgp++)


  //calculation of mass matrix: According to the paper of Jelenic and Crisfield "Geometrically exact 3D beam theory: implementation of a strain-invariant
  //finite element for statics and dynamics", 1999, page 146, a time integration scheme that delivers angular velocities and angular accelerations as
  //needed for the inertia terms of geometrically exact beams has to be based on multiplicative rotation angle increments between two successive time
  //steps. Since BACI does all displacement updates in an additive manner, the global vector of rotational displacements has no physical meaning and,
  //consequently the global velocity and acceleration vectors resulting from the BACI time integration schemes have no physical meaning, too. Therefore,
  //a mass matrix in combination with this global acceleration vector is meaningless from a physical point of view. For these reasons, we have to apply
  //our own time integration scheme at element level. Up to now, the only implemented integration scheme is the gen-alpha Lie group time integration
  //according to [Arnold, Brüls (2007)], [Brüls, Cardona, 2010] and [Brüls, Cardona, Arnold (2012)] in combination with a constdisvelacc predictor. (Christoph Meier, 04.14)

  const double beta = params.get<double>("rot_beta",1000);

  if (massmatrix != NULL and inertia_force != NULL)
  {

    if(beta < 999)
    {
      #ifdef BEAM3RAUTOMATICDIFF
      #ifndef MULTIPLICATIVEUPDATES
        dserror("beam3r: for dynamic simulations in combination with a FAD stiffness matrix, please define flag MULTIPLICATIVEUPDATES because the implemented GenAlpha Lie group time integration is based on multiplicative updates!");
      #endif
      #endif

      const double gamma = params.get<double>("rot_gamma",1000);
      const double alpha_f = params.get<double>("rot_alphaf",1000);
      const double alpha_m = params.get<double>("rot_alpham",1000);
      const double dt = params.get<double>("delta time",1000);
      const bool materialintegration=true;
      const double diff_factor_vel = gamma/(beta*dt);
      const double diff_factor_acc = (1.0-alpha_m)/(beta*dt*dt*(1.0-alpha_f));

      LINALG::Matrix<3,3> Lambdanewmass(true);
      LINALG::Matrix<3,3> Lambdaconvmass(true);

      //first of all we get the material law
      Teuchos::RCP<const MAT::Material> currmat = Material();
      double rho = 0.0;

      //assignment of material parameters; only St.Venant material is accepted for this beam
      switch(currmat->MaterialType())
      {
        case INPAR::MAT::m_stvenant:// only linear elastic material supported
        {
          const MAT::StVenantKirchhoff* actmat = static_cast<const MAT::StVenantKirchhoff*>(currmat.get());
          rho = actmat->Density();
        }
        break;
        default:
          dserror("unknown or improper type of material law");
        break;
      }

      //tensor of mass moments of inertia and cross-section value. These values are used in order to artificially scale
      //the the translational and rotational inertia terms with given input parameters if necessary:
      LINALG::Matrix<3,3> Jp(true);
      Jp(0,0)=inertscalerot1_*(Iyy_+Izz_);
      Jp(1,1)=inertscalerot2_*Iyy_;
      Jp(2,2)=inertscalerot2_*Izz_;
      Jp.Scale(rho);

      const double scaledcrosssec=inertscaletrans_*crosssec_;

      //integration points for mass matrix (exact integration)
      DRT::UTILS::IntegrationPoints1D gausspoints_exact(MyGaussRule(nnode,gaussexactintegration));

      //evaluate all shape functions at all specified Gauss points
      EvaluateShapeFunctions<nnode>(gausspoints_exact,Imass_i,this->Shape());

      //Calculate current (translational) displacement at gauss points (needed for element intern time integration)
      for (int gp=0; gp<gausspoints_exact.nquad; gp++)//loop through Gauss points
      {
        for (int i=0; i<3; ++i)
        {
          dispnewmass_[gp](i) = 0.0;
          for (unsigned int node=0; node<nnode; ++node)
          {
            dispnewmass_[gp](i) += disp[6*node+i]*Imass_i[gp](node);
          }
        }
      }

      Ekin_=0.0;
      L_=0.0;
      P_=0.0;

      for (int gp=0; gp<gausspoints_exact.nquad; gp++)//loop through Gauss points
      {
        //weight of GP in parameter space
        const double wgtmass = gausspoints_exact.qwgt[gp];

        LINALG::Matrix<3,3> Jp_bar(Jp);
        Jp_bar.Scale(diff_factor_acc);

        LINALG::Matrix<3,1> dL(true);

        // update quaternions at GPs for exact Gauss quadrature
        Calc_Psi_l<nnode,FADordouble>(Psi_li, Imass_i[gp], Psi_l);
        Calc_Qgauss<double>(FADUTILS::CastToDouble<FADordouble,3,1>(Psi_l),FADUTILS::CastToDouble<FADordouble,4,1>(Q_r),Qnewmass_[gp]);
        computeItilde<nnode>(Psi_l,Itildemass,Phi_IJ,Lambda_r,Psi_li,Imass_i[gp]);

        Lambdanewmass.Clear();
        Lambdaconvmass.Clear();
        //compute current and old triad at Gauss point
        LARGEROTATIONS::quaterniontotriad<double>(Qnewmass_[gp],Lambdanewmass);
        LARGEROTATIONS::quaterniontotriad<double>(Qconvmass_[gp],Lambdaconvmass);

        //rotation between last converged position and current position expressend as a quaternion
        LINALG::Matrix<4,1>  deltaQ(true);
        LARGEROTATIONS::quaternionproduct<double>(LARGEROTATIONS::inversequaternion<double>(Qconvmass_[gp]),Qnewmass_[gp],deltaQ);

        //spatial rotation between last converged position and current position expressed as a three element rotation vector
        LINALG::Matrix<3,1> deltatheta(true);
        LARGEROTATIONS::quaterniontoangle<double>(deltaQ,deltatheta);

        //compute material counterparts of spatial vectors
        LINALG::Matrix<3,1> deltaTHETA(true);
        LINALG::Matrix<3,1> Wconvmass(true);
        LINALG::Matrix<3,1> Wnewmass(true);
        LINALG::Matrix<3,1> Aconvmass(true);
        LINALG::Matrix<3,1> Anewmass(true);
        LINALG::Matrix<3,1> Amodconvmass(true);
        LINALG::Matrix<3,1> Amodnewmass(true);
        deltaTHETA.MultiplyTN(Lambdanewmass,deltatheta);
        Wconvmass.MultiplyTN(Lambdaconvmass,wconvmass_[gp]);
        Aconvmass.MultiplyTN(Lambdaconvmass,aconvmass_[gp]);
        Amodconvmass.MultiplyTN(Lambdaconvmass,amodconvmass_[gp]);

        //update angular velocities and accelerations according to Newmark time integration scheme either in
        //material description (see Jelenic, 1999, p. 146, equations (2.8) and (2.9)) or in spatial description
        //(for testing purposes, not recommended by Jelenic). The corresponding equations are adapted according to
        //the gen-alpha Lie group time integration scheme proposed in [Arnold, Brüls (2007)], [Brüls, Cardona, 2010]
        //and [Brüls, Cardona, Arnold (2012)]. In the predictor step of the time integration the following
        //formulas automatically deliver a constant displacement (deltatheta=0), consistent velocity and consistent acceleration
        //predictor. This fact has to be reflected in a consistent manner by the choice of the predictor in the input file:
        if (materialintegration)
        {
          for (int i=0;i<3;i++)
          {
            Anewmass(i)=   (1.0-alpha_m)/(beta*dt*dt*(1.0-alpha_f))*deltaTHETA(i)-(1.0-alpha_m)/(beta*dt*(1.0-alpha_f))*Wconvmass(i)      \
                          -alpha_f/(1.0-alpha_f)*Aconvmass(i)+(alpha_m/(1.0-alpha_f)-(0.5-beta)*(1.0-alpha_m)/(beta*(1.0-alpha_f)))*Amodconvmass(i);

            Wnewmass(i)=gamma/(beta*dt)*deltaTHETA(i)+(1-gamma/beta)*Wconvmass(i)+dt*(1-gamma/(2*beta))*Amodconvmass(i);

            Amodnewmass(i)=1.0/(1.0-alpha_m)*((1.0-alpha_f)*Anewmass(i) + alpha_f*Aconvmass(i) - alpha_m*Amodconvmass(i));
          }
          wnewmass_[gp].Multiply(Lambdanewmass,Wnewmass);
          anewmass_[gp].Multiply(Lambdanewmass,Anewmass);
          amodnewmass_[gp].Multiply(Lambdanewmass,Amodnewmass);
        }
        else
        {
          for (int i=0;i<3;i++)
          {
            wnewmass_[gp](i)=gamma/(beta*dt)*deltatheta(i)+(1-gamma/beta)*wconvmass_[gp](i)+dt*(1-gamma/(2*beta))*amodconvmass_[gp](i);

            anewmass_[gp](i)=   (1.0-alpha_m)/(beta*dt*dt*(1.0-alpha_f))*deltatheta(i)-(1.0-alpha_m)/(beta*dt*(1.0-alpha_f))*wconvmass_[gp](i)      \
                          -alpha_f/(1.0-alpha_f)*aconvmass_[gp](i)+(alpha_m/(1.0-alpha_f)-(0.5-beta)*(1.0-alpha_m)/(beta*(1.0-alpha_f)))*amodconvmass_[gp](i);

            amodnewmass_[gp](i)=1.0/(1.0-alpha_m)*((1.0-alpha_f)*anewmass_[gp](i) + alpha_f*aconvmass_[gp](i) - alpha_m*amodconvmass_[gp](i) );
          }
          Wnewmass.MultiplyTN(Lambdanewmass,wnewmass_[gp]);
          Anewmass.MultiplyTN(Lambdanewmass,anewmass_[gp]);
          Amodnewmass.MultiplyTN(Lambdanewmass,amodnewmass_[gp]);
        }

        LINALG::Matrix<3,1> deltar(true);
        for (int i=0;i<3;i++)
        {
          deltar(i)=dispnewmass_[gp](i)-dispconvmass_[gp](i);
        }
        for (int i=0;i<3;i++)
        {
          rttnewmass_[gp](i)=   (1.0-alpha_m)/(beta*dt*dt*(1.0-alpha_f))*deltar(i)-(1.0-alpha_m)/(beta*dt*(1.0-alpha_f))*rtconvmass_[gp](i)      \
                        -alpha_f/(1.0-alpha_f)*rttconvmass_[gp](i)+(alpha_m/(1.0-alpha_f)-(0.5-beta)*(1.0-alpha_m)/(beta*(1.0-alpha_f)))*rttmodconvmass_[gp](i);

          rtnewmass_[gp](i)=gamma/(beta*dt)*deltar(i)+(1-gamma/beta)*rtconvmass_[gp](i)+dt*(1-gamma/(2*beta))*rttmodconvmass_[gp](i);

          rttmodnewmass_[gp](i)=1.0/(1.0-alpha_m)*((1.0-alpha_f)*rttnewmass_[gp](i) + alpha_f*rttconvmass_[gp](i) - alpha_m*rttmodconvmass_[gp](i) );
        }

        //spin matrix of the material angular velocity, i.e. S(W)
        LINALG::Matrix<3,3> SWnewmass(true);
        LARGEROTATIONS::computespin<double>(SWnewmass,Wnewmass);
        LINALG::Matrix<3,1> Jp_Wnewmass(true);
        LINALG::Matrix<3,1> auxvector1(true);
        LINALG::Matrix<3,1> Pi_t(true);
        Jp_Wnewmass.Multiply(Jp,Wnewmass);
        for (int i=0;i<3;i++)
          for (int j=0;j<3;j++)
            auxvector1(i)+=SWnewmass(i,j)*Jp_Wnewmass(j)+Jp(i,j)*Anewmass(j);

        Pi_t.Multiply(Lambdanewmass,auxvector1);
        LINALG::Matrix<3,1> r_tt(true);
        LINALG::Matrix<3,1> r_t(true);
        LINALG::Matrix<3,1> r(true);
        for (int i=0; i<3; ++i)
          for (unsigned int node=0; node<nnode; ++node)
          {
            r(i) += (Nodes()[node]->X()[i]+disp[6*node+i])*Imass_i[gp](node);
          }

        for (int i=0; i<3; ++i)
        {
          r_tt(i) = rttnewmass_[gp](i);
          r_t(i) =  rtnewmass_[gp](i);
        }

        LINALG::Matrix<3,3> S_r(true);
        LARGEROTATIONS::computespin<double>(S_r,r);
        dL.Multiply(S_r,r_t);
        dL.Scale(rho*scaledcrosssec);
        LINALG::Matrix<3,1> Lambdanewmass_Jp_Wnewmass(true);
        Lambdanewmass_Jp_Wnewmass.Multiply(Lambdanewmass,Jp_Wnewmass);
        dL.Update(1.0,Lambdanewmass_Jp_Wnewmass,1.0);
        for (int i=0;i<3;i++)
        {
          L_(i)+=wgtmass*jacobimass_[gp]*dL(i);
          P_(i)+=wgtmass*jacobimass_[gp]*rho*scaledcrosssec*r_t(i);
        }

        LINALG::Matrix<3,3> S_Pit(true);
        LARGEROTATIONS::computespin<double>(S_Pit,Pi_t);
        LINALG::Matrix<3,3> SJpWnewmass(true);
        LARGEROTATIONS::computespin<double>(SJpWnewmass,Jp_Wnewmass);
        LINALG::Matrix<3,3> SWnewmass_Jp(true);
        SWnewmass_Jp.Multiply(SWnewmass,Jp);
        Jp_bar.Update(diff_factor_vel,SWnewmass_Jp,1.0);
        Jp_bar.Update(-diff_factor_vel,SJpWnewmass,1.0);

        LINALG::Matrix<3,3> Tmatrix(true);
        Tmatrix=LARGEROTATIONS::Tmatrix(deltatheta);

        LINALG::Matrix<3,3> Lambdanewmass_Jpbar(true);
        Lambdanewmass_Jpbar.Multiply(Lambdanewmass, Jp_bar);
        LINALG::Matrix<3,3> LambdaconvmassT_Tmatrix(true);
        LambdaconvmassT_Tmatrix.MultiplyTN(Lambdaconvmass, Tmatrix);
        LINALG::Matrix<3,3> Lambdanewmass_Jpbar_LambdaconvmassT_Tmatrix(true);
        Lambdanewmass_Jpbar_LambdaconvmassT_Tmatrix.Multiply(Lambdanewmass_Jpbar, LambdaconvmassT_Tmatrix);
        LINALG::Matrix<3,3> auxmatrix1(true);
        auxmatrix1.Update(-1.0,S_Pit,1.0);
        auxmatrix1.Update(1.0,Lambdanewmass_Jpbar_LambdaconvmassT_Tmatrix,1.0);

        //inertia forces
        for (unsigned int i=0;i<nnode;i++)
        {
          for (int j=0;j<3;j++)
          {
            //translational contribution
            (*inertia_force)(6*i+j)  += jacobimass_[gp]*wgtmass*rho*scaledcrosssec*Imass_i[gp](i)*r_tt(j);
            //rotational contribution
            (*inertia_force)(6*i+j+3)+= jacobimass_[gp]*wgtmass*Imass_i[gp](i)*Pi_t(j);
          }
        }

        //linearization of inertia forces
        for (unsigned int j=0;j<nnode;j++)
        {
          //translational contribution
          for (unsigned int i=0;i<nnode;i++)
          {
            for (int k=0;k<3;k++)
              (*massmatrix)(6*i+k,6*j+k)+=diff_factor_acc*jacobimass_[gp]*wgtmass*rho*scaledcrosssec*Imass_i[gp](i)*Imass_i[gp](j);
          }

          //rotational contribution
          LINALG::Matrix<3,3> auxmatrix2(true);
          auxmatrix2.Multiply(auxmatrix1,Itildemass[j]);
          for (unsigned int i=0;i<nnode;i++)
          {
            for (int k=0;k<3;k++)
            {
              for (int l=0;l<3;l++)
              {
                (*massmatrix)(6*i+3+k,6*j+3+l)+=jacobimass_[gp]*wgtmass*Imass_i[gp](i)*auxmatrix2(k,l);
              }
            }
          }
        }

        //Calculation of kinetic energy
        LINALG::Matrix<1,1> ekinrot(true);
        LINALG::Matrix<1,1> ekintrans(true);
        ekinrot.MultiplyTN(Wnewmass,Jp_Wnewmass);
        ekintrans.MultiplyTN(r_t,r_t);
        Ekin_+=0.5*(ekinrot.Norm2() + rho*scaledcrosssec*ekintrans.Norm2())*jacobimass_[gp]*wgtmass;
        kintorsionenergy_+=0.5* Wnewmass(0)*Jp_Wnewmass(0)*jacobimass_[gp]*wgtmass;
        kinbendingenergy_+=0.5* Wnewmass(1)*Jp_Wnewmass(1)*jacobimass_[gp]*wgtmass;
        kinbendingenergy_+=0.5* Wnewmass(2)*Jp_Wnewmass(2)*jacobimass_[gp]*wgtmass;
        kintransenergy_+=0.5*rho*scaledcrosssec*ekintrans.Norm2()*jacobimass_[gp]*wgtmass;

        Jp_Wnewmass.Multiply(Jp,Wnewmass);
      }//for (int gp=0; gp<gausspoints_exact.nquad; gp++)
    }
    else
    {
      //This is a dummy mass matrix which is necessary for statmech simulations
      for (unsigned int i=0; i<6*nnode; i++)
        (*massmatrix)(i,i) = 1;
    }
  }//if (massmatrix != NULL and inertia_force != NULL)

  // in statistical mechanics simulations, a deletion influenced by the values of the internal force vector might occur   // TODO what's this ???
  if(gausspoints.nquad==1 && params.get<std::string>("internalforces","no")=="yes" && force != NULL)
  {
    eps_ = FADUTILS::CastToDouble(Gamma(0));
    f_ = *force;
    for (int i=0; i<3; ++i)
      Ngp_(i) = FADUTILS::CastToDouble(stressN(i));
  }

  /*the function CalcBrownian applies statistical forces and damping matrix according to the fluctuation dissipation theorem;
  * it is dedicated to the application of beam elements in the frame of statistical mechanics problems; for these problems a
  * special vector has to be passed to the element packed in the params parameter list; in case that the control routine calling
  * the element does not attach this special vector to params the following method is just doing nothing, which means that for
  * any ordinary problem of structural mechanics it may be ignored*/
  if(needstatmech_)
  {
    // the following is pre-calculated for Statmech specific methods (historic reasons, feel free to adapt it)
    //vector with nnode elements, who represent the 3x3-matrix-shaped interpolation function \tilde{I}^nnode at all Gauss points for mass matrix according to (3.19), Jelenic 1999
    std::vector<std::vector<LINALG::Matrix<3,3> > > Itildemass_allGP(nnode);
    std::vector<LINALG::TMatrix<double,3,3> >  Itildemass(nnode);

    for (unsigned int i=0; i<Itildemass_allGP.size(); ++i)
      Itildemass_allGP[i].resize(nnode);

    //integration points for mass matrix (exact integration)
    DRT::UTILS::IntegrationPoints1D gausspoints_exact(MyGaussRule(nnode,gaussexactintegration));

    //evaluate all shape functions at all specified Gauss points
    EvaluateShapeFunctions<nnode>(gausspoints_exact,Imass_i,this->Shape());

    // TODO do we really need to update the class variable dispnewmass_ in case of pure Statmech simulations?
    //Calculate current (translational) displacement at gauss points
    for (int gp=0; gp<gausspoints_exact.nquad; gp++)//loop through Gauss points
    {
      for (int i=0; i<3; ++i)
      {
        dispnewmass_[gp](i) = 0.0;
        for (unsigned int node=0; node<nnode; ++node)
        {
          dispnewmass_[gp](i) += disp[6*node+i]*Imass_i[gp](node);
        }
      }
    }

    for (int gp=0; gp<gausspoints_exact.nquad; gp++)//loop through Gauss points
    {
      // update quaternions at GPs for exact Gauss quadrature
      Calc_Psi_l<nnode,FADordouble>(Psi_li, Imass_i[gp], Psi_l);
      Calc_Qgauss<double>(FADUTILS::CastToDouble<FADordouble,3,1>(Psi_l),FADUTILS::CastToDouble<FADordouble,4,1>(Q_r),Qnewmass_[gp]);
      computeItilde<nnode>(Psi_l,Itildemass,Phi_IJ,Lambda_r,Psi_li,Imass_i[gp]);

      for (unsigned int inode=0; inode<Itildemass.size(); ++inode)
        for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
            Itildemass_allGP[gp][inode](i,j) = Itildemass[inode](i,j);
    }

  #ifndef BEAM3RCONSTSTOCHFORCE
    CalcBrownian<nnode,3,6,4>(params,vel,disp,stiffmatrix,force,Imass_i,Itildemass_allGP);
  #else
    CalcBrownian<nnode,3,6,3>(params,vel,disp,stiffmatrix,force,Imass_i,Itildemass_allGP);
  #endif
  }

  return;
} // DRT::ELEMENTS::Beam3r::b3_nlnstiffmass

/*------------------------------------------------------------------------------------------------------------*
 | lump mass matrix             (private)                                                   cyron 01/08|
 *------------------------------------------------------------------------------------------------------------*/
template<unsigned int nnode>
void DRT::ELEMENTS::Beam3r::lumpmass(Epetra_SerialDenseMatrix* emass)
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

/*-----------------------------------------------------------------------------------------------------------*
 | Evaluate PTC damping (public)                                                                  cyron 10/08|
 *----------------------------------------------------------------------------------------------------------*/
template<unsigned int nnode>
void DRT::ELEMENTS::Beam3r::EvaluatePTC(Teuchos::ParameterList& params,
                                      Epetra_SerialDenseMatrix& elemat1)
{
  //apply PTC rotation damping term using a Lobatto integration rule; implemented for 2 nodes only
  if(nnode > 2)
    dserror("PTC implemented for 2-noded elements only");

  for (unsigned int node=0; node<nnode; node++)
  {
    //computing angle increment from current position in comparison with last converged position for damping
    LINALG::Matrix<4,1> deltaQ;
    LARGEROTATIONS::quaternionproduct(LARGEROTATIONS::inversequaternion(Qconv_[node]),Qnew_[node],deltaQ);
    LINALG::Matrix<3,1> deltatheta;
    LARGEROTATIONS::quaterniontoangle(deltaQ,deltatheta);

    //isotropic artificial stiffness
    LINALG::Matrix<3,3> artstiff;
    artstiff = LARGEROTATIONS::Tmatrix(deltatheta);

    //scale artificial damping with crotptc parameter for PTC method
    artstiff.Scale( params.get<double>("crotptc",0.0) );

    //each node gets a block diagonal damping term; the Lobatto integration weight is 0.5 for 2-noded elements
    for(int k=0; k<3; k++)
      for (int l=0; l<3; l++)
        elemat1(node*6+3+k,node*6+3+l) += artstiff(k,l)*0.5*jacobinode_[node];

    //PTC for translational degrees of freedom; the Lobatto integration weight is 0.5 for 2-noded elements
    for(int k=0; k<3; k++)
      elemat1(node*6+k,node*6+k) += params.get<double>("ctransptc",0.0)*0.5*jacobinode_[node];
  }

  return;
} //DRT::ELEMENTS::Beam3r::EvaluatePTC

/*-----------------------------------------------------------------------------------------------------------*
 | computes damping coefficients per length and stores them in a matrix in the following order: damping of    |
 | translation parallel to filament axis, damping of translation orthogonal to filament axis, damping of     |
 | rotation around filament axis                                             (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
inline void DRT::ELEMENTS::Beam3r::MyDampingConstants(Teuchos::ParameterList& params,LINALG::Matrix<3,1>& gamma)
{
  //translational damping coefficients according to Howard, p. 107, table 6.2;
  gamma(0) = 2*PI*params.get<double>("ETA",0.0);
  gamma(1) = 4*PI*params.get<double>("ETA",0.0);

  /*damping coefficient of rigid straight rod spinning around its own axis according to Howard, p. 107, table 6.2;
   *as this coefficient is very small for thin rods it is increased artificially by a factor for numerical convencience*/
  double rsquare = std::sqrt(4*Iyy_/PI);
  //TODO: Here the damping constants are artificially enhanced!!!
  double artificial = 4000;//4000;//50;  20000//50 not bad for standard Actin3D_10.dat files; for 40 elements also 1 seems to work really well; for large networks 4000 seems good (artificial contribution then still just ~0.1 % of nodal moments)
  gamma(2) = 4*PI*params.get<double>("ETA",0.0)*rsquare*artificial;

  //in case of an isotropic friction model the same damping coefficients are applied parallel to the polymer axis as perpendicular to it
  if(DRT::INPUT::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL") == INPAR::STATMECH::frictionmodel_isotropicconsistent || DRT::INPUT::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL") == INPAR::STATMECH::frictionmodel_isotropiclumped)
    gamma(0) = gamma(1);

   /* in the following section damping coefficients are replaced by those suggested in Ortega2003, which allows for a
    * comparison of the finite element simulation with the results of that article; note that we assume that the element
    * length is equivalent to the particle length in the following when computing the length to diameter ratio p*/
/*
   double lrefe= 0.3;
   //for (int gp=0; gp<nnode-1; gp++)
     //lrefe += gausspointsdamping.qwgt[gp]*jacobi_[gp];

   double p=lrefe/(pow(crosssec_*4.0/PI,0.5));
   double Ct=0.312+0.565/p-0.100/pow(p,2.0);
   double Cr=-0.662+0.917/p-0.05/pow(p,2.0);
   gamma(0) = 2.0*PI*params.get<double>("ETA",0.0)/(log(p) + 2*Ct - Cr);
   gamma(1) = 4.0*PI*params.get<double>("ETA",0.0)/(log(p)+Cr);
   gamma(2) = 4.0*PI*params.get<double>("ETA",0.0)*rsquare*artificial*(0.96 + 0.64992/p - 0.17568/(p*p));
*/
}

/*-----------------------------------------------------------------------------------------------------------*
 |computes the number of different random numbers required in each time step for generation of stochastic    |
 |forces;                                                                    (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Beam3r::HowManyRandomNumbersINeed()
{
  /*at each Gauss point one needs as many random numbers as randomly excited degrees of freedom, i.e. three
   *random numbers for the translational degrees of freedom and one random number for the rotation around the element axis*/
  #ifndef BEAM3RCONSTSTOCHFORCE
    return (4*NumNode());
  #else
    return (3);
  #endif
}

/*-----------------------------------------------------------------------------------------------------------*
 |computes velocity of background fluid and gradient of that velocity at a certain evaluation point in       |
 |the physical space                                                         (public)           cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<int ndim> //number of dimensions of embedding space
void DRT::ELEMENTS::Beam3r::MyBackgroundVelocity(Teuchos::ParameterList&      params,  //!<parameter list
                                                const LINALG::Matrix<ndim,1>&  evaluationpoint,  //!<point at which background velocity and its gradient has to be computed
                                                LINALG::Matrix<ndim,1>&        velbackground,  //!< velocity of background fluid
                                                LINALG::Matrix<ndim,ndim>&     velbackgroundgrad) //!<gradient of velocity of background fluid
{
  /*note: this function is not yet a general one, but always assumes a shear flow, where the velocity of the
   * background fluid is always directed in direction params.get<int>("DBCDISPDIR",0) and orthogonal to z-axis.
   * In 3D the velocity increases linearly in z and equals zero for z = 0.
   * In 2D the velocity increases linearly in y and equals zero for y = 0. */

  //velocity at upper boundary of domain
  double uppervel = 0.0;

  //default values for background velocity and its gradient
  velbackground.PutScalar(0);
  velbackgroundgrad.PutScalar(0);

  double time = params.get<double>("total time",0.0);
  double starttime = params.get<double>("STARTTIMEACT",0.0);
  double dt = params.get<double>("delta time");
  double shearamplitude = params.get<double> ("SHEARAMPLITUDE", 0.0);
  int curvenumber = params.get<int> ("CURVENUMBER", -1)-1;
  int dbcdispdir = params.get<int> ("DBCDISPDIR", -1)-1;

  Teuchos::RCP<std::vector<double> > defvalues = Teuchos::rcp(new std::vector<double>(3,0.0));
  Teuchos::RCP<std::vector<double> > periodlength = params.get("PERIODLENGTH", defvalues);
  INPAR::STATMECH::DBCType dbctype = params.get<INPAR::STATMECH::DBCType>("DBCTYPE", INPAR::STATMECH::dbctype_std);
  bool shearflow = false;
  if(dbctype==INPAR::STATMECH::dbctype_shearfixed ||
     dbctype==INPAR::STATMECH::dbctype_shearfixeddel ||
     dbctype==INPAR::STATMECH::dbctype_sheartrans ||
     dbctype==INPAR::STATMECH::dbctype_affineshear||
     dbctype==INPAR::STATMECH::dbctype_affinesheardel)
    shearflow = true;

  //oscillations start only at params.get<double>("STARTTIMEACT",0.0)
  if(periodlength->at(0) > 0.0)
    if(shearflow && time>starttime && fabs(time-starttime)>dt/1e4 && curvenumber >=  0 && dbcdispdir >= 0 )
    {
      uppervel = shearamplitude * (DRT::Problem::Instance()->Curve(curvenumber).FctDer(time,1))[1];

      //compute background velocity
      velbackground(dbcdispdir) = (evaluationpoint(ndim-1) / periodlength->at(ndim-1)) * uppervel;

      //compute gradient of background velocity
      velbackgroundgrad(dbcdispdir,ndim-1) = uppervel / periodlength->at(ndim-1);
    }
}

/*-----------------------------------------------------------------------------------------------------------*
 | computes rotational damping forces and stiffness (public)                                    cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<unsigned int nnode> //number of nodes
inline void DRT::ELEMENTS::Beam3r::MyRotationalDamping(Teuchos::ParameterList&                       params,  //!<parameter list
                                              const std::vector<double>&                              vel,  //!< element velocity vector
                                              const std::vector<double>&                              disp, //!<element disp vector
                                              Epetra_SerialDenseMatrix*                               stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseVector*                               force,  //!< element internal force vector
                                              const DRT::UTILS::IntegrationPoints1D&                  gausspointsdamping,
                                              const std::vector<LINALG::Matrix<1,nnode> >&            Idamping,
                                              const std::vector<std::vector<LINALG::Matrix<3,3> > >&  Itildedamping,
                                              const std::vector<LINALG::Matrix<4,1> >&                Qconvdamping,
                                              const std::vector<LINALG::Matrix<4,1> >&                Qnewdamping)
{
  //get time step size
  double dt = params.get<double>("delta time",0.0);

  //auxiliary matrices
  LINALG::Matrix<3,3> sum;
  LINALG::Matrix<3,3> auxmatrix;
  LINALG::Matrix<3,3> Lambdadamping;

  //determine type of numerical integration performed (lumped damping matrix via lobatto integration!)
  std::vector<double> jacobi(jacobimass_);
  if(DRT::INPUT::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL") == INPAR::STATMECH::frictionmodel_isotropiclumped)
    jacobi = jacobinode_;

  //damping coefficients for translational and rotational degrees of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma);

  for (int gp=0; gp<gausspointsdamping.nquad; gp++)//loop through Gauss points
  {
    //compute triad at Gauss point
    LARGEROTATIONS::quaterniontotriad(Qnewdamping[gp],Lambdadamping);

    //rotation between last converged position and current position expressed as a quaternion
    LINALG::Matrix<4,1>  deltaQ;
    LARGEROTATIONS::quaternionproduct(LARGEROTATIONS::inversequaternion(Qconvdamping[gp]),Qnewdamping[gp],deltaQ);

    //rotation between last converged position and current position expressed as a three element rotation vector
    LINALG::Matrix<3,1> deltatheta;
    LARGEROTATIONS::quaterniontoangle(deltaQ,deltatheta);

    //angular velocity at this Gauss point according to backward Euler scheme
    LINALG::Matrix<3,1> omega(true);
    omega += deltatheta;
    omega.Scale(1/dt);

    //compute matrix T*W*T^t
    LINALG::Matrix<3,3> TWTt;
    for(int k=0; k<3; k++)
      for(int j = 0; j<3; j++)
        TWTt(k,j) = (Lambdadamping)(k,0)*(Lambdadamping)(j,0);

    //compute vector T*W*T^t*\omega
    LINALG::Matrix<3,1> TWTtomega;
    TWTtomega.Multiply(TWTt,omega);

    //compute matrix T*W*T^t*H^(-1)
    LINALG::Matrix<3,3> TWTtHinv;
    TWTtHinv.Multiply(TWTt,LARGEROTATIONS::Tmatrix(deltatheta));

    //compute spin matrix S(\omega)
    LINALG::Matrix<3,3> Sofomega;
    LARGEROTATIONS::computespin(Sofomega,omega);

    //compute matrix T*W*T^t*S(\omega)
    LINALG::Matrix<3,3> TWTtSofomega;
    TWTtSofomega.Multiply(TWTt,Sofomega);

    //compute spin matrix S(T*W*T^t*\omega)
    LINALG::Matrix<3,3> SofTWTtomega;
    LARGEROTATIONS::computespin(SofTWTtomega,TWTtomega);

    //loop over all line nodes
    for(unsigned int i=0; i<nnode; i++)
      //loop over three dimensions in line direction
      for(int k=0; k<3; k++)
      {
        if(force != NULL)
          (*force)(i*6+3+k) += gamma(2)*TWTtomega(k)*(Idamping[gp])(i)*gausspointsdamping.qwgt[gp]*jacobi[gp];

        if(stiffmatrix != NULL)
          //loop over all column nodes
          for (unsigned int j=0; j<nnode; j++)
            //loop over three dimensions in column direction
            for (int l=0; l<3; l++)
            {
              sum.PutScalar(0);
              sum += TWTtHinv;
              sum.Scale(1/dt);
              sum += TWTtSofomega;
              sum -= SofTWTtomega;

              auxmatrix.Multiply(sum,(Itildedamping[gp])[j]);

              (*stiffmatrix)(i*6+3+k,j*6+3+l) += gamma(2)*auxmatrix(k,l)*(Idamping[gp])(i)*gausspointsdamping.qwgt[gp]*jacobi[gp];
            }
      }
  }

  return;
}//DRT::ELEMENTS::Beam3r::MyRotationalDamping(.)

/*-----------------------------------------------------------------------------------------------------------*
 | computes translational damping forces and stiffness (public)                                 cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<unsigned int nnode, int ndim, int dof> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node
inline void DRT::ELEMENTS::Beam3r::MyTranslationalDamping(Teuchos::ParameterList& params,  //!<parameter list
                                                  const std::vector<double>&       vel,  //!< element velocity vector
                                                  const std::vector<double>&       disp, //!<element disp vector
                                                  Epetra_SerialDenseMatrix*        stiffmatrix,  //!< element stiffness matrix
                                                  Epetra_SerialDenseVector*        force)//!< element internal force vector
{
  //get time step size
  double dt = params.get<double>("delta time",0.0);

  //velocity and gradient of background velocity field
  LINALG::Matrix<ndim,1> velbackground;
  LINALG::Matrix<ndim,ndim> velbackgroundgrad;

  //evaluation point in physical space corresponding to a certain Gauss point in parameter space
  LINALG::Matrix<ndim,1> evaluationpoint;

  //damping coefficients for translational and rotatinal degrees of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma);

  //get vector jacobi with Jacobi determinants at each integration point (gets by default those values required for consistent damping matrix)
  std::vector<double> jacobi(jacobimass_);

  //determine type of numerical integration performed (lumped damping matrix via lobatto integration!)
  IntegrationType integrationtype = gaussexactintegration;
  if(DRT::INPUT::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL") == INPAR::STATMECH::frictionmodel_isotropiclumped)
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
    for(unsigned int i=0; i<nnode; i++)
      //loop over all dimensions
      for(int j=0; j<ndim; j++)
        evaluationpoint(j) += funct(i)*(Nodes()[i]->X()[j]+disp[dof*i+j]);

    //compute velocity and gradient of background flow field at evaluationpoint
    MyBackgroundVelocity<ndim>(params,evaluationpoint,velbackground,velbackgroundgrad);


    //compute tangent vector t_{\par} at current Gauss point
    LINALG::Matrix<ndim,1> tpar(true);
    for(unsigned int i=0; i<nnode; i++)
      for(int k=0; k<ndim; k++)
        tpar(k) += deriv(i)*(Nodes()[i]->X()[k]+disp[dof*i+k]) / jacobi[gp];

    //compute velocity vector at this Gauss point
    LINALG::Matrix<ndim,1> velgp(true);
    for(unsigned int i=0; i<nnode; i++)
      for(int l=0; l<ndim; l++)
        velgp(l) += funct(i)*vel[dof*i+l];

    //compute matrix product (t_{\par} \otimes t_{\par}) \cdot velbackgroundgrad
    LINALG::Matrix<ndim,ndim> tpartparvelbackgroundgrad(true);
    for(int i=0; i<ndim; i++)
      for(int j=0; j<ndim; j++)
        for(int k=0; k<ndim; k++)
          tpartparvelbackgroundgrad(i,j) += tpar(i)*tpar(k)*velbackgroundgrad(k,j);

    //loop over all line nodes
    for(unsigned int i=0; i<nnode; i++)
      //loop over lines of matrix t_{\par} \otimes t_{\par}
      for(int k=0; k<ndim; k++)
        //loop over columns of matrix t_{\par} \otimes t_{\par}
        for(int l=0; l<ndim; l++)
        {
          if(force != NULL)
            (*force)(i*dof+k)+= funct(i)*jacobi[gp]*gausspoints.qwgt[gp]*( (k==l)*gamma(1) + (gamma(0) - gamma(1))*tpar(k)*tpar(l) ) *(velgp(l)- velbackground(l));

          if(stiffmatrix != NULL)
            //loop over all column nodes
            for (unsigned int j=0; j<nnode; j++)
            {
              (*stiffmatrix)(i*dof+k,j*dof+l) += gausspoints.qwgt[gp]*funct(i)*funct(j)*jacobi[gp]*(                 (k==l)*gamma(1) + (gamma(0) - gamma(1))*tpar(k)*tpar(l) ) / dt;
              (*stiffmatrix)(i*dof+k,j*dof+l) -= gausspoints.qwgt[gp]*funct(i)*funct(j)*jacobi[gp]*( velbackgroundgrad(k,l)*gamma(1) + (gamma(0) - gamma(1))*tpartparvelbackgroundgrad(k,l) ) ;
              (*stiffmatrix)(i*dof+k,j*dof+k) += gausspoints.qwgt[gp]*funct(i)*deriv(j)*                                               (gamma(0) - gamma(1))*tpar(l)*(velgp(l) - velbackground(l));
              (*stiffmatrix)(i*dof+k,j*dof+l) += gausspoints.qwgt[gp]*funct(i)*deriv(j)*                                               (gamma(0) - gamma(1))*tpar(k)*(velgp(l) - velbackground(l));
            }
        }
  }

  return;
}//DRT::ELEMENTS::Beam3r::MyTranslationalDamping(.)

/*-----------------------------------------------------------------------------------------------------------*
 | computes stochastic forces and resulting stiffness (public)                                  cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<unsigned int nnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
inline void DRT::ELEMENTS::Beam3r::MyStochasticForces(Teuchos::ParameterList& params,  //!<parameter list
                                              const std::vector<double>&       vel,  //!< element velocity vector
                                              const std::vector<double>&       disp, //!<element disp vector
                                              Epetra_SerialDenseMatrix*        stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseVector*        force)//!< element internal force vector
{
  //damping coefficients for three translational and one rotational degree of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma);

  //get vector jacobi with Jacobi determinants at each integration point (gets by default those values required for consistent damping matrix)
  std::vector<double> jacobi(jacobimass_);

  //determine type of numerical integration performed (lumped damping matrix via lobatto integration!)
  IntegrationType integrationtype = gaussexactintegration;
  if(DRT::INPUT::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL") == INPAR::STATMECH::frictionmodel_isotropiclumped)
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
  Teuchos::RCP<Epetra_MultiVector> randomnumbers = params.get<  Teuchos::RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null);

  for(int gp=0; gp < gausspoints.nquad; gp++)
  {
    //evaluate basis functions and their derivatives at current Gauss point
    DRT::UTILS::shape_function_1D(funct,gausspoints.qxg[gp][0],Shape());
    DRT::UTILS::shape_function_1D_deriv1(deriv,gausspoints.qxg[gp][0],Shape());

    //compute tangent vector t_{\par} at current Gauss point
    LINALG::Matrix<ndim,1> tpar(true);
    for(unsigned int i=0; i<nnode; i++)
      for(int k=0; k<ndim; k++)
        tpar(k) += deriv(i)*(Nodes()[i]->X()[k]+disp[dof*i+k]) / jacobi[gp];


    //loop over all line nodes
    for(unsigned int i=0; i<nnode; i++)
      //loop dimensions with respect to lines
      for(int k=0; k<ndim; k++)
        //loop dimensions with respect to columns
        for(int l=0; l<ndim; l++)
        {
          if(force != NULL)
          {
            #ifndef BEAM3RCONSTSTOCHFORCE
              (*force)(i*dof+k) -= funct(i)*(sqrt(gamma(1))*(k==l) + (sqrt(gamma(0)) - sqrt(gamma(1)))*tpar(k)*tpar(l))*(*randomnumbers)[gp*randompergauss+l][LID()]*sqrt(jacobi[gp]*gausspoints.qwgt[gp]);
            #else
              (*force)(i*dof+k) -= funct(i)*(sqrt(gamma(1))*(k==l) + (sqrt(gamma(0)) - sqrt(gamma(1)))*tpar(k)*tpar(l))*(*randomnumbers)[l][LID()]*sqrt(jacobi[gp]*gausspoints.qwgt[gp]);
            #endif
          }

          if(stiffmatrix != NULL)
            //loop over all column nodes
            for (unsigned int j=0; j<nnode; j++)
            {
              #ifndef BEAM3RCONSTSTOCHFORCE
                (*stiffmatrix)(i*dof+k,j*dof+k) -= funct(i)*deriv(j)*tpar(l)*(*randomnumbers)[gp*randompergauss+l][LID()]*sqrt(gausspoints.qwgt[gp]/ jacobi[gp])*(sqrt(gamma(0)) - sqrt(gamma(1)));
                (*stiffmatrix)(i*dof+k,j*dof+l) -= funct(i)*deriv(j)*tpar(k)*(*randomnumbers)[gp*randompergauss+l][LID()]*sqrt(gausspoints.qwgt[gp]/ jacobi[gp])*(sqrt(gamma(0)) - sqrt(gamma(1)));
              #else
                (*stiffmatrix)(i*dof+k,j*dof+k) -= funct(i)*deriv(j)*tpar(l)*(*randomnumbers)[l][LID()]*sqrt(gausspoints.qwgt[gp]/ jacobi[gp])*(sqrt(gamma(0)) - sqrt(gamma(1)));
                (*stiffmatrix)(i*dof+k,j*dof+l) -= funct(i)*deriv(j)*tpar(k)*(*randomnumbers)[l][LID()]*sqrt(gausspoints.qwgt[gp]/ jacobi[gp])*(sqrt(gamma(0)) - sqrt(gamma(1)));
              #endif
            }
        }
  }

  return;
}//DRT::ELEMENTS::Beam3r::MyStochasticForces(.)

/*-----------------------------------------------------------------------------------------------------------*
 | computes stochastic moments and (if required) resulting stiffness (public)                   cyron   10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<unsigned int nnode, int randompergauss> //number of nodes, number of random numbers required per Gauss point, number of random numbers required per Gauss point
inline void DRT::ELEMENTS::Beam3r::MyStochasticMoments(Teuchos::ParameterList&                       params,  //!<parameter list
                                              const std::vector<double>&                              vel,  //!< element velocity vector
                                              const std::vector<double>&                              disp, //!<element disp vector
                                              Epetra_SerialDenseMatrix*                               stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseVector*                               force, //!< element internal force vector
                                              const DRT::UTILS::IntegrationPoints1D&                  gausspointsdamping,
                                              const std::vector<LINALG::Matrix<1,nnode> >&            Idamping,
                                              const std::vector<std::vector<LINALG::Matrix<3,3> > >&  Itildedamping,
                                              const std::vector<LINALG::Matrix<4,1> >&                Qconvdamping,
                                              const std::vector<LINALG::Matrix<4,1> >&                Qnewdamping)
{
  //auxiliary matrix
  LINALG::Matrix<3,3> auxmatrix;

  //determine type of numerical integration performed (lumped damping matrix via lobatto integration!)
  std::vector<double> jacobi(jacobimass_);
  if(DRT::INPUT::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL") == INPAR::STATMECH::frictionmodel_isotropiclumped)
    jacobi = jacobinode_;

  //damping coefficients for three translational and one rotational degree of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma);

  /*get pointer at Epetra multivector in parameter list linking to random numbers for stochastic forces with zero mean
   * and standard deviation (2*kT / dt)^0.5; note carefully: a space between the two subsequal ">" signs is mandatory
   * for the C++ parser in order to avoid confusion with ">>" for streams*/
   Teuchos::RCP<Epetra_MultiVector> randomnumbers = params.get<  Teuchos::RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null);

  for(int gp=0; gp < gausspointsdamping.nquad; gp++)
  {
    //get first column out of current triad at Gauss point
    LARGEROTATIONS::quaterniontotriad(Qnewdamping[gp],auxmatrix);
    LINALG::Matrix<3,1> t1;
    for(int i=0; i<3; i++)
      t1(i) = auxmatrix(i,0);

    //compute spin matrix from first column of Tnew times random number
    LINALG::Matrix<3,3> S;
    LARGEROTATIONS::computespin(S,t1);
    S.Scale((*randomnumbers)[gp*randompergauss+3][LID()]);

    //loop over all line nodes
    for(int i=0; i<nnode; i++)
      //loop over lines of matrix t_{\par} \otimes t_{\par}
      for(int k=0; k<3; k++)
      {
        if(force != NULL)
          (*force)(i*6+3+k) -= (Idamping[gp])(i)*t1(k)*(*randomnumbers)[gp*randompergauss+3][LID()]*sqrt(jacobi[gp]*gausspointsdamping.qwgt[gp]*gamma(2));

        if(stiffmatrix != NULL)
          //loop over all column nodes
          for (int j=0; j<nnode; j++)
            //loop over three dimensions with respect to columns
            for(int l=0; l<3; l++)
            {
              auxmatrix.Multiply(S,(Itildedamping[gp])[j]);
              (*stiffmatrix)(i*6+3+k,j*6+3+l) += (Idamping[gp])(i)*auxmatrix(k,l)*sqrt(jacobi[gp]*gausspointsdamping.qwgt[gp]*gamma(2));
            }

    }
  }
  return;
}//DRT::ELEMENTS::Beam3r::MyStochasticMoments(.)

/*-----------------------------------------------------------------------------------------------------------*
 | Assemble stochastic and viscous forces and respective stiffness according to fluctuation dissipation      |
 | theorem                                                                               (public) cyron 10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<unsigned int nnode, int ndim, int dof, int randompergauss> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node, number of random numbers required per Gauss point
inline void DRT::ELEMENTS::Beam3r::CalcBrownian(Teuchos::ParameterList&                       params,
                                              const std::vector<double>&                       vel,  //!< element velocity vector
                                              const std::vector<double>&                       disp, //!< element displacement vector
                                              Epetra_SerialDenseMatrix*                        stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseVector*                        force,
                                              std::vector<LINALG::Matrix<1,nnode> >&           Imass,
                                              std::vector<std::vector<LINALG::Matrix<3,3> > >& Itildemass) //!< element internal force vector
{
  /*for integration of damping matrix always nnode Gauss points required; but in case of Lobatto integration
   *these are identical to the nnode nodes and then the basis functions are no longer the one also required
   *for the mass matrix, but rather their values at the integration points are given by a Kronecker-Delta function*/
  IntegrationType dampingintrule(gaussexactintegration);
  if(DRT::INPUT::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL") == INPAR::STATMECH::frictionmodel_isotropiclumped)
    dampingintrule = lobattointegration;

  DRT::UTILS::IntegrationPoints1D gausspointsdamping(MyGaussRule(nnode,dampingintrule));
  std::vector<LINALG::Matrix<1,nnode> > Idamping(Imass);
  std::vector<std::vector<LINALG::Matrix<3,3> > > Itildedamping(Itildemass);
  std::vector<LINALG::Matrix<4,1> > Qconvdamping(Qconvmass_);
  std::vector<LINALG::Matrix<4,1> > Qnewdamping(Qnewmass_);

  if(DRT::INPUT::get<INPAR::STATMECH::FrictionModel>(params,"FRICTION_MODEL") == INPAR::STATMECH::frictionmodel_isotropiclumped)
  {
    //in case of Lobatto integration nodal triads are triads at Gauss points
    Qconvdamping = Qconv_;
    Qnewdamping  = Qnew_;

    //loop over all Gauss points
    for(int i=0; i < gausspointsdamping.nquad; i++)
      //loop over all nodes
      for(unsigned int j=0; j<nnode; j++)
      {
        if(i == (int)j)
          (Idamping[i])(j) = 1;
        else
          (Idamping[i])(j) = 0;
      }

    //loop through all Gauss points
    for(int i=0; i < gausspointsdamping.nquad; i++)
      //loop through all nodes to calculate respective basis function matrix
      for(unsigned int j=0; j<nnode; j++)
        for(int k=0; k < 3; k++)
          for(int l=0; l < 3; l++)
          {
            if(i == (int)j && k == l)
              ((Itildedamping[i])[j])(k,l) = 1;
            else
              ((Itildedamping[i])[j])(k,l) = 0;
          }
  }

  //now start with evaluation of force vectors and stiffness matrices

  //add stiffness and forces due to translational damping effects
  MyTranslationalDamping<nnode,ndim,dof>(params,vel,disp,stiffmatrix,force);

  //add stiffness and forces (i.e. moments) due to rotational damping effects
  MyRotationalDamping<nnode>(params,vel,disp,stiffmatrix,force,gausspointsdamping,Idamping,Itildedamping,Qconvdamping,Qnewdamping);

  //add stochastic forces and (if required) resulting stiffness
  MyStochasticForces<nnode,ndim,dof,randompergauss>(params,vel,disp,stiffmatrix,force);

  //add stochastic moments and resulting stiffness // TODO what about StochsticMoments???
  //MyStochasticMoments<nnode,randompergauss>(params,vel,disp,stiffmatrix,force,gausspointsdamping,Idamping,Itildedamping,Qconvdamping,Qnewdamping);

return;
}//DRT::ELEMENTS::Beam3r::CalcBrownian(.)

/*-----------------------------------------------------------------------------------------------------------*
 | shifts nodes so that proper evaluation is possible even in case of periodic boundary conditions; if two   |
 | nodes within one element are separated by a periodic boundary, one of them is shifted such that the final |
 | distance in R^3 is the same as the initial distance in the periodic space; the shift affects computation  |
 | on element level within that very iteration step, only (no change in global variables performed)          |                                 |
 |                                                                                       (public) cyron 10/09|
 *----------------------------------------------------------------------------------------------------------*/
template<unsigned int nnode, int ndim> //number of nodes, number of dimensions
inline void DRT::ELEMENTS::Beam3r::NodeShift(Teuchos::ParameterList& params,  //!<parameter list
                                              std::vector<double>&    disp) //!<element disp vector
{
  /*get number of degrees of freedom per node; note: the following function assumes the same number of degrees
   *of freedom for each element node*/
  int numdof = NumDofPerNode(*(Nodes()[0]));

  double time = params.get<double>("total time",0.0);
  double starttime = params.get<double>("STARTTIMEACT",0.0);
  double dt = params.get<double>("delta time");
  double shearamplitude = params.get<double> ("SHEARAMPLITUDE", 0.0);
  int curvenumber = params.get<int> ("CURVENUMBER", -1)-1;
  int dbcdispdir = params.get<int> ("DBCDISPDIR", -1)-1;

  Teuchos::RCP<std::vector<double> > defvalues = Teuchos::rcp(new std::vector<double>(3,0.0));
  Teuchos::RCP<std::vector<double> > periodlength = params.get("PERIODLENGTH", defvalues);
  INPAR::STATMECH::DBCType dbctype = params.get<INPAR::STATMECH::DBCType>("DBCTYPE", INPAR::STATMECH::dbctype_std);
  bool shearflow = false;
  if(dbctype==INPAR::STATMECH::dbctype_shearfixed || dbctype==INPAR::STATMECH::dbctype_sheartrans || dbctype==INPAR::STATMECH::dbctype_affineshear)
    shearflow = true;

  /*only if periodic boundary conditions are in use, i.e. params.get<double>("PeriodLength",0.0) > 0.0, this
   * method has to change the displacement variables*/
  if(periodlength->at(0) > 0.0)
    //loop through all nodes except for the first node which remains fixed as reference node
    for(unsigned int i=1;i<nnode;i++)
    {
      for(int dof= ndim - 1; dof > -1; dof--)
      {
        /*if the distance in some coordinate direction between some node and the first node becomes smaller by adding or subtracting
         * the period length, the respective node has obviously been shifted due to periodic boundary conditions and should be shifted
         * back for evaluation of element matrices and vectors; this way of detecting shifted nodes works as long as the element length
         * is smaller than half the periodic length*/
        if( fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) + periodlength->at(dof) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) < fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) )
        {
          disp[numdof*i+dof] += periodlength->at(dof);

          /*the upper domain surface orthogonal to the z-direction may be subject to shear Dirichlet boundary condition; the lower surface
           *may be fixed by DBC. To avoid problmes when nodes exit the domain through the upper z-surface and reenter through the lower
           *z-surface, the shear has to be substracted from nodal coordinates in that case */
          if(shearflow && dof == 2 && curvenumber >=  0 && time>starttime && fabs(time-starttime)>dt/1e4)
            disp[numdof*i+dbcdispdir] += shearamplitude*DRT::Problem::Instance()->Curve(curvenumber).f(time);
        }

        if( fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - periodlength->at(dof) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) < fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) )
        {
          disp[numdof*i+dof] -= periodlength->at(dof);

          /*the upper domain surface orthogonal to the z-direction may be subject to shear Dirichlet boundary condition; the lower surface
           *may be fixed by DBC. To avoid problmes when nodes exit the domain through the lower z-surface and reenter through the upper
           *z-surface, the shear has to be added to nodal coordinates in that case */
          if(shearflow && dof == 2 && curvenumber >=  0 && time>starttime && fabs(time-starttime)>dt/1e4 )
            disp[numdof*i+dbcdispdir] -= shearamplitude*DRT::Problem::Instance()->Curve(curvenumber).f(time);
        }
      }
    }

return;
}//DRT::ELEMENTS::Beam3r::NodeShift
