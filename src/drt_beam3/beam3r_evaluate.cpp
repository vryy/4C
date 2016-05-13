/*!-----------------------------------------------------------------------------------------------------------
\file beam3r_evaluate.cpp

\brief evaluation methods for 3D nonlinear Reissner beam element

\level 2

\maintainer Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262


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
#include "../drt_structure_new/str_elements_paramsinterface.H"

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
  SetParamsInterfacePtr(params);

  // start with "none"
  ELEMENTS::ActionType act = ELEMENTS::none;

  if (IsParamsInterface())
  {
   act = ParamsInterface().GetActionType();
  }
  else
  {
    // get the action required
    std::string action = params.get<std::string>("action","calc_none");
    if (action == "calc_none")                    dserror("No action supplied");
    else if (action=="calc_struct_linstiff")                               act = ELEMENTS::struct_calc_linstiff;
    else if (action=="calc_struct_nlnstiff")                               act = ELEMENTS::struct_calc_nlnstiff;
    else if (action=="calc_struct_internalforce")                          act = ELEMENTS::struct_calc_internalforce;
    else if (action=="calc_struct_linstiffmass")                           act = ELEMENTS::struct_calc_linstiffmass;
    else if (action=="calc_struct_nlnstiffmass")                           act = ELEMENTS::struct_calc_nlnstiffmass;
    else if (action=="calc_struct_nlnstifflmass")                          act = ELEMENTS::struct_calc_nlnstifflmass; //with lumped mass matrix
    else if (action=="calc_struct_stress")                                 act = ELEMENTS::struct_calc_stress;
    else if (action=="calc_struct_eleload")                                act = ELEMENTS::struct_calc_eleload;
    else if (action=="calc_struct_fsiload")                                act = ELEMENTS::struct_calc_fsiload;
    else if (action=="calc_struct_update_istep")                           act = ELEMENTS::struct_calc_update_istep;
    else if (action=="calc_struct_reset_istep")                            act = ELEMENTS::struct_calc_reset_istep;
    else if (action=="calc_struct_ptcstiff")                               act = ELEMENTS::struct_calc_ptcstiff;
    else if (action=="calc_struct_energy")                                 act = ELEMENTS::struct_calc_energy;
    else dserror("Unknown type of action for Beam3r");
  }

  // nnodetriad: number of nodes used for interpolation of triad field
  const int nnodetriad = NumNode();

  switch(act)
  {
    case ELEMENTS::struct_calc_ptcstiff:
    {
      switch(nnodetriad)
       {
           case 2:EvaluatePTC<2>(params, elemat1); break;
           case 3:EvaluatePTC<3>(params, elemat1); break;
           case 4:EvaluatePTC<4>(params, elemat1); break;
           case 5:EvaluatePTC<5>(params, elemat1); break;
           default:dserror("Only Line2, Line3, Line4 and Line5 Elements implemented.");
       }
    }
    break;

    case ELEMENTS::struct_calc_linstiff:
    {
      // only nonlinear case implemented!
      dserror("linear stiffness matrix called, but not implemented");
    }
    break;

    case ELEMENTS::struct_calc_energy:
    {
      if(elevec1 != Teuchos::null)
      {
        if(elevec1.M()!=1)
          dserror("energy vector of invalid size %i, expected row dimension 1 (total elastic energy of element)!", elevec1.M());
        elevec1(0)=Eint_;
      }

    }
    break;

    // nonlinear stiffness and mass matrix are calculated even if only nonlinear stiffness matrix is required
    case ELEMENTS::struct_calc_nlnstiffmass:
    case ELEMENTS::struct_calc_nlnstifflmass:
    case ELEMENTS::struct_calc_nlnstiff:
    case ELEMENTS::struct_calc_internalforce:
    {
      // need current global displacement and residual forces and get them from discretization
      // making use of the local-to-global map lm one can extract current displacement and residual values for each degree of freedom

      // get element displacements
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==Teuchos::null) dserror("Cannot get state vectors 'displacement'");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

      std::vector<double> myvel(lm.size());
      std::vector<double> myacc(lm.size());       // TODO myacc unused? reuse if Lie-group time integration happens outside the element at some day

      /* only in case of Statmech application (calculation of damping and stochastic forces), get element velocities
       * note: condition if(needstatmech_) alone would call GetState("velocity") also in cases where not needed and
       * therefore not accounted for on discretization level => error "Cannot find state velocity" in drt_discret.H*/
      if (needstatmech_ and (params.get<Teuchos::RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null) != Teuchos::null) )
      {
        Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState("velocity");
        DRT::UTILS::ExtractMyValues(*vel,myvel,lm);
      }

      if (act == ELEMENTS::struct_calc_nlnstiffmass)
      {
        switch(nnodetriad)
        {
          case 2:
          {
            if (!centerline_hermite_)
              b3_nlnstiffmass<2,2,1>(params,myacc,myvel,mydisp,&elemat1,&elemat2,&elevec1,&elevec2);
            else
              b3_nlnstiffmass<2,2,2>(params,myacc,myvel,mydisp,&elemat1,&elemat2,&elevec1,&elevec2);
            break;
          }
          case 3:
          {
            if (!centerline_hermite_)
              b3_nlnstiffmass<3,3,1>(params,myacc,myvel,mydisp,&elemat1,&elemat2,&elevec1,&elevec2);
            else
              b3_nlnstiffmass<3,2,2>(params,myacc,myvel,mydisp,&elemat1,&elemat2,&elevec1,&elevec2);
            break;
          }
          case 4:
          {
            if (!centerline_hermite_)
              b3_nlnstiffmass<4,4,1>(params,myacc,myvel,mydisp,&elemat1,&elemat2,&elevec1,&elevec2);
            else
              b3_nlnstiffmass<4,2,2>(params,myacc,myvel,mydisp,&elemat1,&elemat2,&elevec1,&elevec2);
            break;
          }
          case 5:
          {
            if (!centerline_hermite_)
              b3_nlnstiffmass<5,5,1>(params,myacc,myvel,mydisp,&elemat1,&elemat2,&elevec1,&elevec2);
            else
              b3_nlnstiffmass<5,2,2>(params,myacc,myvel,mydisp,&elemat1,&elemat2,&elevec1,&elevec2);
            break;
          }
        }
      }

      else if (act == ELEMENTS::struct_calc_nlnstifflmass)
      {
        // TODO there is a method 'Beam3r::lumpmass'; check generality and functionality and enable action here
        dserror("Lumped mass matrix not implemented for beam3r elements so far!");
      }

      else if (act == ELEMENTS::struct_calc_nlnstiff)
      {
        switch(nnodetriad)
        {
          case 2:
          {
            if (!centerline_hermite_)
              b3_nlnstiffmass<2,2,1>(params,myacc,myvel,mydisp,&elemat1,NULL,&elevec1,NULL);
            else
              b3_nlnstiffmass<2,2,2>(params,myacc,myvel,mydisp,&elemat1,NULL,&elevec1,NULL);
            break;
          }
          case 3:
          {
            if (!centerline_hermite_)
              b3_nlnstiffmass<3,3,1>(params,myacc,myvel,mydisp,&elemat1,NULL,&elevec1,NULL);
            else
              b3_nlnstiffmass<3,2,2>(params,myacc,myvel,mydisp,&elemat1,NULL,&elevec1,NULL);
            break;
          }
          case 4:
          {
            if (!centerline_hermite_)
              b3_nlnstiffmass<4,4,1>(params,myacc,myvel,mydisp,&elemat1,NULL,&elevec1,NULL);
            else
              b3_nlnstiffmass<4,2,2>(params,myacc,myvel,mydisp,&elemat1,NULL,&elevec1,NULL);
            break;
          }
          case 5:
          {
            if (!centerline_hermite_)
              b3_nlnstiffmass<5,5,1>(params,myacc,myvel,mydisp,&elemat1,NULL,&elevec1,NULL);
            else
              b3_nlnstiffmass<5,2,2>(params,myacc,myvel,mydisp,&elemat1,NULL,&elevec1,NULL);
            break;
          }
          default:
            dserror("Only Line2, Line3, Line4, and Line5 Elements implemented.");
        }

      }
      else if (act == ELEMENTS::struct_calc_internalforce)
      {
        switch(nnodetriad)
        {
          case 2:
          {
            if (!centerline_hermite_)
              b3_nlnstiffmass<2,2,1>(params,myacc,myvel,mydisp,NULL,NULL,&elevec1,NULL);
            else
              b3_nlnstiffmass<2,2,2>(params,myacc,myvel,mydisp,NULL,NULL,&elevec1,NULL);
            break;
          }
          case 3:
          {
            if (!centerline_hermite_)
              b3_nlnstiffmass<3,3,1>(params,myacc,myvel,mydisp,NULL,NULL,&elevec1,NULL);
            else
              b3_nlnstiffmass<3,2,2>(params,myacc,myvel,mydisp,NULL,NULL,&elevec1,NULL);
            break;
          }
          case 4:
          {
            if (!centerline_hermite_)
              b3_nlnstiffmass<4,4,1>(params,myacc,myvel,mydisp,NULL,NULL,&elevec1,NULL);
            else
              b3_nlnstiffmass<4,2,2>(params,myacc,myvel,mydisp,NULL,NULL,&elevec1,NULL);
            break;
          }
          case 5:
          {
            if (!centerline_hermite_)
              b3_nlnstiffmass<5,5,1>(params,myacc,myvel,mydisp,NULL,NULL,&elevec1,NULL);
            else
              b3_nlnstiffmass<5,2,2>(params,myacc,myvel,mydisp,NULL,NULL,&elevec1,NULL);
            break;
          }
          default:
            dserror("Only Line2, Line3, Line4, and Line5 Elements implemented.");
        }
      }

    /* at the end of an iteration step the geometric configuration has to be updated: the starting point for the
     * next iteration step is the configuration at the end of the current step */
    Qoldnode_ = Qnewnode_;
    dispthetaoldnode_= dispthetanewnode_;
    }
    break;

    case ELEMENTS::struct_calc_update_istep:
    {
      /* the action calc_struct_update_istep is called in the very end of a time step when the new dynamic
       * equilibrium has finally been found; this is the point where the variable representing the geometric
       * status of the beam at the end of the time step has to be stored*/
      Qconvnode_ = Qnewnode_;
      QconvGPmass_ = QnewGPmass_;
      wconvGPmass_ = wnewGPmass_;
      aconvGPmass_ = anewGPmass_;
      amodconvGPmass_ = amodnewGPmass_;
      rttconvGPmass_ = rttnewGPmass_;
      rttmodconvGPmass_ = rttmodnewGPmass_;
      rtconvGPmass_ = rtnewGPmass_;
      rconvGPmass_ = rnewGPmass_;
      dispthetaconvnode_ = dispthetanewnode_;
    }
    break;

    case ELEMENTS::struct_calc_reset_istep:
    {
      /* the action calc_struct_reset_istep is called by the adaptive time step controller; carries out one test
       * step whose purpose is only figuring out a suitabel timestep; thus this step may be a very bad one in order
       * to iterated towards the new dynamic equilibrium and the thereby gained new geometric configuration should
       * not be applied as starting point for any further iteration step; as a consequence the thereby generated change
       * of the geometric configuration should be canceled and the configuration should be reset to the value at the
       * beginning of the time step*/
      Qoldnode_ = Qconvnode_;
      dispthetaoldnode_ = dispthetaconvnode_;
      Qnewnode_ = Qconvnode_;
      dispthetanewnode_ = dispthetaconvnode_;
      QnewGPmass_=QconvGPmass_;
      wnewGPmass_=wconvGPmass_;
      anewGPmass_=aconvGPmass_;
      amodnewGPmass_=amodconvGPmass_;
      rttnewGPmass_=rttconvGPmass_;
      rttmodnewGPmass_=rttmodconvGPmass_;
      rtnewGPmass_=rtconvGPmass_;
      rnewGPmass_=rconvGPmass_;
    }
    break;

    case ELEMENTS::struct_calc_stress:
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

  // nnodetriad: number of nodes used for interpolation of triad field
  const int nnodetriad = NumNode();
  // nnodecl: number of nodes used for interpolation of centerline
  // assumptions: nnodecl<=nnodetriad; centerline nodes have local ID 0...nnodecl-1
  int nnodecl = nnodetriad;
  if (centerline_hermite_) nnodecl=2;

  // vpernode: number of interpolated values per node (1: value (i.e. Lagrange), 2: value + derivative of value (i.e. Hermite))
  int vpernode=1;
  if (centerline_hermite_)
    vpernode=2;

  // number of DOFs per node depending on type of node
  const int dofperclnode = 3*vpernode;
  const int dofpertriadnode = 3;
  const int dofpercombinode = dofperclnode+dofpertriadnode;

  const DiscretizationType distype = this->Shape();

  // find out whether we will use a time curve and get the factor
  const std::vector<int>* curve  = condition.Get<std::vector<int> >("curve");
  // amplitude of load curve at current time called, 6 components (3 forces, 3 moments)
  std::vector<double> curvefac(6,1.0);

  for (int i=0; i<6; ++i)
  {
    int curvenum = -1;
    // number of the load curve related with a specific line Neumann condition called
    if (curve) curvenum = (*curve)[i];

    if (curvenum>=0 && usetime)
      curvefac[i] = DRT::Problem::Instance()->Curve(curvenum).f(time);
  }

  // gaussian points
  const DRT::UTILS::IntegrationPoints1D intpoints(MyGaussRule(neumann_lineload));

  // declaration of variables in order to store shape functions
  // used for interpolation of triad field
  Epetra_SerialDenseVector I_i(nnodetriad);
  // used for interpolation of centerline
  Epetra_SerialDenseVector H_i(vpernode*nnodecl);

  // get values and switches from the condition

  // onoff is related to the first numdf flags of a line Neumann condition in the input file;
  // value 1 for flag i says that condition is active for i-th degree of freedom
  const std::vector<int>* onoff = condition.Get<std::vector<int> >("onoff");
  // val is related to the numdf "val" fields after the onoff flags of the Neumann condition
  // in the input file; val gives the values of the force as a multiple of the prescribed load curve
  const std::vector<double>* val = condition.Get<std::vector<double> >("val");
  // funct is related to the numdf "funct" fields after the val field of the Neumann condition
  // in the input file; funct gives the number of the function defined in the section FUNCT
  const std::vector<int>* functions = condition.Get<std::vector<int> >("funct");

  // integration points in parameter space and weights
  double xi=0.0;
  double wgt=0.0;

  // integration loops
  for (int numgp=0; numgp<intpoints.nquad; ++numgp)
  {
    xi = intpoints.qxg[numgp][0];
    wgt = intpoints.qwgt[numgp];

    // evaluation of shape functions at Gauss points
    DRT::UTILS::shape_function_1D(I_i,xi,distype);
    if (centerline_hermite_)
      DRT::UTILS::shape_function_hermite_1D(H_i,xi,reflength_,line2);
    else
      DRT::UTILS::shape_function_1D(H_i,xi,distype);

    // position vector at the gauss point at reference configuration needed for function evaluation
    std::vector<double> X_ref(3,0.0);

    // calculate coordinates of corresponding Gauss point in reference configuration
    for (int node=0;node<nnodecl;node++)
    {
      for (int dim=0;dim<3;dim++)
      {
        X_ref[dim]+=H_i[vpernode*node]*Nodes()[node]->X()[dim];

        if (centerline_hermite_)
          X_ref[dim]+=H_i[vpernode*node+1]*(Trefnode_[node])(dim);
      }
    }

    double fac=0;
    fac = wgt * jacobiGPneumannline_[numgp];

    // load vector ar
    double ar[6];

    // loop the relevant dofs of a node
    for (int dof=0; dof<6; ++dof)
      ar[dof] = fac * (*onoff)[dof] * (*val)[dof] * curvefac[dof];
    double functionfac = 1.0;
    int functnum = -1;

    //sum up load components
    for (int dof=0; dof<6; ++dof)
    {
      if (functions)
        functnum = (*functions)[dof];
      else
        functnum = -1;

      // evaluate function at the position of the current GP
      if (functnum>0)
        functionfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(dof, &X_ref[0], time, NULL);   // TODO X_ref[0] is only x-coordinate. is this done on purpose?
      else
        functionfac = 1.0;

      for (int node=0; node<nnodecl; ++node)
      {
        if (dof<3)
        {
          elevec1[dofpercombinode*node+dof] += H_i[vpernode*node] *ar[dof] *functionfac;

          if (centerline_hermite_)
            elevec1[dofpercombinode*node+6+dof] += H_i[vpernode*node+1] *ar[dof] *functionfac;
        }
        else // dof<6
          elevec1[dofpercombinode*node+dof] += I_i[node] *ar[dof] *functionfac;
      }

      for (int node=nnodecl; node<nnodetriad; ++node)
        if (dof>2 && dof<6)
          elevec1[dofperclnode*nnodecl+dofpertriadnode*node+dof-3] += I_i[node] *ar[dof] *functionfac;
    }
  } // for (int numgp=0; numgp<intpoints.nquad; ++numgp)

  return 0;
}

/*----------------------------------------------------------------------------------------------------------------------*
 | get constitutive matrices from material law                                                               grill 03/16|
 *----------------------------------------------------------------------------------------------------------------------*/
template <typename T>
inline void DRT::ELEMENTS::Beam3r::GetConstitutiveMatrices(LINALG::TMatrix<T,3,3>& CN,
                                                           LINALG::TMatrix<T,3,3>& CM)
{
  // first of all we get the material law
   Teuchos::RCP<const MAT::Material> currmat = Material();
   double ym = 0.0;
   double sm = 0.0;

   // assignment of material parameters; only St.Venant material is accepted for this beam
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

   // defining material constitutive matrix CN between Gamma and N according to Jelenic 1999, section 2.4
   CN.PutScalar(0);
   CN(0,0) = ym*crosssec_;
   CN(1,1) = sm*crosssecshear_;
   CN(2,2) = sm*crosssecshear_;

   // defining material constitutive matrix CM between curvature and moment according to Jelenic 1999, section 2.4
   CM.PutScalar(0);
   CM(0,0) = sm*Irr_;
   CM(1,1) = ym*Iyy_;
   CM(2,2) = ym*Izz_;

  return;
}

/*----------------------------------------------------------------------------------------------------------------------*
 |push forward material stress vector and constitutive matrix to their spatial counterparts by rotation matrix Lambda   |
 |according to Romero 2004, eq. (3.10)                                                                       cyron 04/10|
 *----------------------------------------------------------------------------------------------------------------------*/
template <typename T>
inline void DRT::ELEMENTS::Beam3r::pushforward(const LINALG::TMatrix<T,3,3>& Lambda,
                                                const LINALG::TMatrix<T,3,1>& stress_mat,
                                                const LINALG::TMatrix<T,3,3>& C_mat,
                                                LINALG::TMatrix<T,3,1>& stress_spatial,
                                                LINALG::TMatrix<T,3,3>& c_spatial)
{
  // introduce auxiliary variable for pushforward of rotational matrices
  LINALG::TMatrix<T,3,3> temp;

  // push forward stress vector
  stress_spatial.Multiply(Lambda,stress_mat);

  // push forward constitutive matrix according to Jelenic 1999, paragraph following to (2.22) on page 148
  temp.Multiply(Lambda,C_mat);
  c_spatial.MultiplyNT(temp,Lambda);

   return;
}

/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and mass matrix (private)                                                   cyron 01/08|
 *-----------------------------------------------------------------------------------------------------------*/
template<unsigned int nnodetriad, unsigned int nnodecl, unsigned int vpernode>
void DRT::ELEMENTS::Beam3r::b3_nlnstiffmass(Teuchos::ParameterList&   params,
                                            std::vector<double>&      acc,
                                            std::vector<double>&      vel,
                                            std::vector<double>&      disp,
                                            Epetra_SerialDenseMatrix* stiffmatrix,
                                            Epetra_SerialDenseMatrix* massmatrix,
                                            Epetra_SerialDenseVector* force,
                                            Epetra_SerialDenseVector* inertia_force)
{
  // nnodetriad: number of nodes used for interpolation of triad field
  // nnodecl: number of nodes used for interpolation of centerline
  // assumptions: nnodecl<=nnodetriad; centerline nodes have local ID 0...nnodecl-1
  // vpernode: number of interpolated values per centerline node (1: value (i.e. Lagrange), 2: value + derivative of value (i.e. Hermite))

  /********************************** Initialize/resize variables **************************************
   *****************************************************************************************************/

  //********************************** quantities valid for entire element *****************************
  const int dofperclnode = 3*vpernode;
  const int dofpertriadnode = 3;
  const int dofpercombinode = dofperclnode+dofpertriadnode;

  // internal force vector
  LINALG::TMatrix<FADordouble, dofperclnode*nnodecl+dofpertriadnode*nnodetriad, 1> f_int(true);

  // reference triad Lambda_r and corresponding quaternion Q_r
  LINALG::TMatrix<FADordouble,3,3> Lambda_r(true);
  LINALG::TMatrix<FADordouble,4,1> Q_r(true);

  // angle of relative rotation between node I and J according to (3.10), Jelenic 1999
  LINALG::TMatrix<FADordouble,3,1> Phi_IJ(true);

  // clear internal (elastic) energy
  Eint_=0.0;

  //**************************************** nodal quantities *******************************************

  // current nodal DOFs relevant for centerline interpolation in total Lagrangian style, i.e. initial values + displacements
  LINALG::TMatrix<FADordouble,3*vpernode*nnodecl,1> disp_totlag_centerline(true);

  // quaternions of all nodal triads
  std::vector<LINALG::TMatrix<FADordouble,4,1> > Q_i(nnodetriad);

  // rotation angles between nodal triads and reference triad according to (3.8), Jelenic 1999
  std::vector<LINALG::TMatrix<FADordouble,3,1> > Psi_li(nnodetriad);

  //*************************** physical quantities evaluated at a certain GP ***************************

  // derivation of beam centerline with respect to arc-length parameter: r'(x) from (2.12), Jelenic 1999
  LINALG::TMatrix<FADordouble,3,1> r_s;
  // spin matrix related to vector r_s
  LINALG::TMatrix<FADordouble,3,3> r_s_hat;
  // interpolated local relative rotation \Psi^l at a certain Gauss point according to (3.11), Jelenic 1999
  LINALG::TMatrix<FADordouble,3,1> Psi_l;
  /* derivative of interpolated local relative rotation \Psi^l with respect to arc-length parameter
   * at a certain Gauss point according to (3.11), Jelenic 1999*/
  LINALG::TMatrix<FADordouble,3,1> Psi_l_s;
  // triad at GP
  LINALG::TMatrix<FADordouble,3,3> Lambda;

  // 3D vector related to spin matrix \hat{\kappa} from (2.1), Jelenic 1999
  LINALG::TMatrix<FADordouble,3,1> K;
  // 3D vector of material axial and shear strains from (2.1), Jelenic 1999
  LINALG::TMatrix<FADordouble,3,1> Gamma;

  // convected stresses N and M and constitutive matrices C_N and C_M according to section 2.4, Jelenic 1999
  LINALG::TMatrix<FADordouble,3,1> stressN;
  LINALG::TMatrix<FADordouble,3,1> stressM;
  LINALG::TMatrix<FADordouble,3,3> CN;
  LINALG::TMatrix<FADordouble,3,3> CM;

  // spatial stresses n and m according to (3.10), Romero 2004 and spatial constitutive matrices c_n and c_m according to page 148, Jelenic 1999
  LINALG::TMatrix<FADordouble,3,1> stressn;
  LINALG::TMatrix<FADordouble,3,1> stressm;
  LINALG::TMatrix<FADordouble,3,3> cn;
  LINALG::TMatrix<FADordouble,3,3> cm;

  //********************************** (generalized) shape functions ************************************
  /* Note: index i refers to the i-th shape function (i = 0 ... nnode*vpernode-1)
   * the vectors store individual shape functions, NOT an assembled matrix of shape functions)*/

  /* vector whose numgp-th element is a 1xnnode-matrix with all Lagrange polynomial shape functions evaluated at the numgp-th Gauss point
   * these shape functions are used for the interpolation of the triad field*/
  std::vector<LINALG::Matrix<1,nnodetriad> > I_i;
  // same for the derivatives
  std::vector<LINALG::Matrix<1,nnodetriad> > I_i_xi;

  /* vector whose numgp-th element is a 1x(vpernode*nnode)-matrix with all (Lagrange/Hermite) shape functions evaluated at the numgp-th GP
   * these shape functions are used for the interpolation of the beam centerline*/
  std::vector<LINALG::Matrix<1,vpernode*nnodecl> > H_i;
  // same for the derivatives
  std::vector<LINALG::Matrix<1,vpernode*nnodecl> > H_i_xi;

  // vector with nnode elements, who represent the 3x3-matrix-shaped interpolation function \tilde{I}^nnode at a certain Gauss point according to (3.18), Jelenic 1999
  std::vector<LINALG::TMatrix<double,3,3> > Itilde(nnodetriad);

  // vector with nnode elements, who represent the 3x3-matrix-shaped interpolation function \tilde{I'}^nnode at a certain Gauss point according to (3.19), Jelenic 1999
  std::vector<LINALG::TMatrix<double,3,3> > Itildeprime(nnodetriad);

  /*************************** update/compute quantities valid for entire element **********************
   *****************************************************************************************************/

  /* first displacement vector is modified for proper element evaluation in case of periodic boundary conditions; in case that
   * no periodic boundary conditions are to be applied the following code line may be ignored or deleted*/
  if(params.isParameter("PERIODLENGTH"))
    NodeShift<nnodecl,3>(params,disp);

  // update disp_totlag
  UpdateDispTotLagAndNodalTriads<nnodetriad,nnodecl,vpernode>(disp,disp_totlag_centerline,Q_i);

  // compute reference triad Lambda_r according to (3.9), Jelenic 1999
  CalcRefQuaternion<FADordouble>(Q_i[nodeI_],Q_i[nodeJ_],Q_r,Phi_IJ);
  LARGEROTATIONS::quaterniontotriad(Q_r,Lambda_r);

  // setup constitutive matrices
  GetConstitutiveMatrices<FADordouble>(CN,CM);

  /* compute nodal local rotations according to (3.8), Jelenic 1999
   * this is done individually for each node in order to avoid function argument std::vector<LINALG::TMatrix<...> >
   * a function with this argument type cannot be called with an argument of type std::vector<LINALG::Matrix<...> > (needed e.g. in SetUpReferenceGeometry) */
  for (unsigned int node=0; node<nnodetriad; ++node)
    CalcPsi_li<FADordouble>(Q_i[node],Q_r,Psi_li[node]);

  /******************************* elasticity: compute fint and stiffmatrix ****************************
   *****************************************************************************************************/

  //************************* residual and stiffmatrix contributions from forces ***********************

  // get integration points for elasticity
  DRT::UTILS::IntegrationPoints1D gausspoints_elast_force(MyGaussRule(res_elastic_force));

  // reuse variables for individual shape functions and resize to new numgp
  I_i.resize(gausspoints_elast_force.nquad);
  H_i_xi.resize(gausspoints_elast_force.nquad);

  // evaluate all shape functions and derivatives with respect to element parameter xi at all specified Gauss points
  EvaluateShapeFunctionsAllGPs<nnodetriad,1>(gausspoints_elast_force,I_i,this->Shape());

  EvaluateShapeFunctionDerivsAllGPs<nnodecl,vpernode>(gausspoints_elast_force,H_i_xi,this->Shape());


  // Loop through all GP and calculate their contribution to the forcevector and stiffnessmatrix
  for(int numgp=0; numgp < gausspoints_elast_force.nquad; numgp++)
  {
    // weight of GP in parameter space
    const double wgt = gausspoints_elast_force.qwgt[numgp];

    Calc_r_s<nnodecl,vpernode,FADordouble>(disp_totlag_centerline, H_i_xi[numgp], jacobiGPelastf_[numgp], r_s);

    Calc_Psi_l<nnodetriad,FADordouble>(Psi_li, I_i[numgp], Psi_l);
    Calc_Lambda<FADordouble>(Psi_l,Q_r,Lambda);

    // compute spin matrix related to vector rprime for later use
    LARGEROTATIONS::computespin<FADordouble>(r_s_hat,r_s);

    // compute material strains Gamma and K
    computeGamma<FADordouble>(r_s,Lambda,GammarefGP_[numgp],Gamma);

    // compute material stresses by multiplying strains with constitutive matrix
    stressN.Multiply(CN,Gamma);

    /* compute spatial stresses and constitutive matrices from convected ones according to Jelenic 1999, page 148, paragraph
     * between (2.22) and (2.23) and Romero 2004, (3.10)*/
    pushforward<FADordouble>(Lambda,stressN,CN,stressn,cn);

    /* computation of internal forces according to Jelenic 1999, eq. (4.3); computation split up with respect
     * to single blocks of matrix in eq. (4.3)*/
    for (unsigned int node=0; node<nnodecl; ++node)
    {
      /* upper left block
       * note: jacobi factor cancels out because it is defined by ds=(ds/dxi)*dxi
       *       and I^{i'} in Jelenic1999 is derivative with respect to arc-length parameter in reference configuration s
       *       which can be computed from I_i_xi by multiplication with the inverse determinant: I^{i'}=I_i_s=I_i_xi*(dxi/ds) */
      for (int k=0; k<3; ++k)
      {
        f_int(dofpercombinode*node+k) += H_i_xi[numgp](vpernode*node)*stressn(k)*wgt;
        if (centerline_hermite_)
          f_int(dofpercombinode*node+6+k) += H_i_xi[numgp](vpernode*node+1)*stressn(k)*wgt;
      }

      // lower left block
      for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
          f_int(dofpercombinode*node+3+i) -= r_s_hat(i,j)*stressn(j)*I_i[numgp](node)*wgt*jacobiGPelastf_[numgp];
    }
    for (unsigned int node=nnodecl; node<nnodetriad; ++node)    // this loop is only entered in case of nnodetriad>nnodecl
    {
      // lower left block
      for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
          f_int(dofperclnode*nnodecl+dofpertriadnode*node+i) -= r_s_hat(i,j)*stressn(j)*I_i[numgp](node)*wgt*jacobiGPelastf_[numgp];
    }

#ifndef BEAM3RAUTOMATICDIFF

    if (stiffmatrix != NULL)
    {
      computeItilde<nnodetriad>(Psi_l,Itilde,Phi_IJ,Lambda_r,Psi_li,I_i[numgp]);

      /* computation of stiffness matrix according to Jelenic 1999, eq. (4.7); computation split up with respect
       * to single blocks of matrix in eq. (4.7).
       * note: again, jacobi factor cancels out in terms whith I^{i'}=I_i_s=I_i_xi*(dxi/ds) (see comment above)
       *       but be careful: Itildeprime and rprime are indeed derivatives with respect to arc-length parameter in reference configuration s */

      // auxiliary variables for storing intermediate matrices in computation of entries of stiffness matrix
      LINALG::Matrix<3,3> auxmatrix1;
      LINALG::Matrix<3,3> auxmatrix2;
      LINALG::Matrix<3,3> auxmatrix3;

      for (unsigned int nodei=0; nodei<nnodecl; nodei++)
      {
        for (unsigned int nodej=0; nodej<nnodecl; nodej++)
        {
          // upper left block
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
            {
              (*stiffmatrix)(dofpercombinode*nodei+i,dofpercombinode*nodej+j) += H_i_xi[numgp](vpernode*nodei)*H_i_xi[numgp](vpernode*nodej)*cn(i,j)*wgt/jacobiGPelastf_[numgp];
              if (centerline_hermite_)
              {
                (*stiffmatrix)(dofpercombinode*nodei+6+i,dofpercombinode*nodej+j) += H_i_xi[numgp](vpernode*nodei+1)*H_i_xi[numgp](vpernode*nodej)*cn(i,j)*wgt/jacobiGPelastf_[numgp];
                (*stiffmatrix)(dofpercombinode*nodei+i,dofpercombinode*nodej+6+j) += H_i_xi[numgp](vpernode*nodei)*H_i_xi[numgp](vpernode*nodej+1)*cn(i,j)*wgt/jacobiGPelastf_[numgp];
                (*stiffmatrix)(dofpercombinode*nodei+6+i,dofpercombinode*nodej+6+j) += H_i_xi[numgp](vpernode*nodei+1)*H_i_xi[numgp](vpernode*nodej+1)*cn(i,j)*wgt/jacobiGPelastf_[numgp];
              }
            }

          // lower left block; note: error in eq. (4.7), Jelenic 1999: the first factor should be I^i instead of I^j
          auxmatrix2.Multiply(r_s_hat,cn);
          LARGEROTATIONS::computespin(auxmatrix1,stressn);
          auxmatrix1 -= auxmatrix2;
          auxmatrix1.Scale(I_i[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
            {
              (*stiffmatrix)(dofpercombinode*nodei+3+i,dofpercombinode*nodej+j) += auxmatrix1(i,j)*H_i_xi[numgp](vpernode*nodej)*wgt;
              if (centerline_hermite_)
                (*stiffmatrix)(dofpercombinode*nodei+3+i,dofpercombinode*nodej+6+j) += auxmatrix1(i,j)*H_i_xi[numgp](vpernode*nodej+1)*wgt;
            }

          // upper right block
          auxmatrix2.Multiply(cn,r_s_hat);
          LARGEROTATIONS::computespin(auxmatrix1,stressn);
          auxmatrix2 -= auxmatrix1;       // auxmatrix2: term in parantheses

          auxmatrix3.Multiply(auxmatrix2,Itilde[nodej]);
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
            {
              (*stiffmatrix)(dofpercombinode*nodei+i,dofpercombinode*nodej+3+j) += auxmatrix3(i,j)*H_i_xi[numgp](vpernode*nodei)*wgt;
              if (centerline_hermite_)
                (*stiffmatrix)(dofpercombinode*nodei+6+i,dofpercombinode*nodej+3+j) += auxmatrix3(i,j)*H_i_xi[numgp](vpernode*nodei+1)*wgt;
            }

          // lower right block
          // third summand; note: error in eq. (4.7), Jelenic 1999: the first summand in the parantheses should be \hat{\Lambda N} instead of \Lambda N
          auxmatrix1.Multiply(auxmatrix2,Itilde[nodej]);  // term in parantheses is the same as in upper right block but with opposite sign (note '-=' below)

          auxmatrix3.Multiply(r_s_hat,auxmatrix1);
          auxmatrix3.Scale(I_i[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofpercombinode*nodei+3+i,dofpercombinode*nodej+3+j) -= auxmatrix3(i,j)*jacobiGPelastf_[numgp]*wgt;

        }
        for (unsigned int nodej=nnodecl; nodej<nnodetriad; nodej++)    // this loop is only entered in case of nnodetriad>nnodecl
        {
          // upper right block
          auxmatrix2.Multiply(cn,r_s_hat);
          LARGEROTATIONS::computespin(auxmatrix1,stressn);
          auxmatrix2 -= auxmatrix1;       // auxmatrix2: term in parantheses

          auxmatrix3.Multiply(auxmatrix2,Itilde[nodej]);
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
            {
              (*stiffmatrix)(dofpercombinode*nodei+i,dofperclnode*nnodecl+dofpertriadnode*nodej+j) += auxmatrix3(i,j)*H_i_xi[numgp](vpernode*nodei)*wgt;
              if (centerline_hermite_)
                (*stiffmatrix)(dofpercombinode*nodei+6+i,dofperclnode*nnodecl+dofpertriadnode*nodej+j) += auxmatrix3(i,j)*H_i_xi[numgp](vpernode*nodei+1)*wgt;
            }

          // lower right block
          // third summand; note: error in eq. (4.7), Jelenic 1999: the first summand in the parantheses should be \hat{\Lambda N} instead of \Lambda N
          auxmatrix1.Multiply(auxmatrix2,Itilde[nodej]);  // term in parantheses is the same as in upper right block but with opposite sign (note '-=' below)

          auxmatrix3.Multiply(r_s_hat,auxmatrix1);
          auxmatrix3.Scale(I_i[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofpercombinode*nodei+3+i,dofperclnode*nnodecl+dofpertriadnode*nodej+j) -= auxmatrix3(i,j)*jacobiGPelastf_[numgp]*wgt;
        }
      }
      for (unsigned int nodei=nnodecl; nodei<nnodetriad; nodei++)    // this loop is only entered in case of nnodetriad>nnodecl
      {
        for (unsigned int nodej=0; nodej<nnodecl; nodej++)
        {
          // lower left block; note: error in eq. (4.7), Jelenic 1999: the first factor should be I^i instead of I^j
          auxmatrix2.Multiply(r_s_hat,cn);
          LARGEROTATIONS::computespin(auxmatrix1,stressn);
          auxmatrix1 -= auxmatrix2;
          auxmatrix1.Scale(I_i[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
            {
              (*stiffmatrix)(dofperclnode*nnodecl+dofpertriadnode*nodei+i,dofpercombinode*nodej+j) += auxmatrix1(i,j)*H_i_xi[numgp](vpernode*nodej)*wgt;
              if (centerline_hermite_)
                (*stiffmatrix)(dofperclnode*nnodecl+dofpertriadnode*nodei+i,dofpercombinode*nodej+6+j) += auxmatrix1(i,j)*H_i_xi[numgp](vpernode*nodej+1)*wgt;
            }

          // lower right block
          // third summand; note: error in eq. (4.7), Jelenic 1999: the first summand in the parantheses should be \hat{\Lambda N} instead of \Lambda N
          auxmatrix2.Multiply(cn,r_s_hat);
          LARGEROTATIONS::computespin(auxmatrix1,stressn);
          auxmatrix2 -= auxmatrix1;       // auxmatrix2: term in parantheses

          auxmatrix1.Multiply(auxmatrix2,Itilde[nodej]);  // term in parantheses is the same as in upper right block but with opposite sign (note '-=' below)

          auxmatrix3.Multiply(r_s_hat,auxmatrix1);
          auxmatrix3.Scale(I_i[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofperclnode*nnodecl+dofpertriadnode*nodei+i,dofpercombinode*nodej+3+j) -= auxmatrix3(i,j)*jacobiGPelastf_[numgp]*wgt;

        }
        for (unsigned int nodej=nnodecl; nodej<nnodetriad; nodej++)
        {
          // lower right block
          // third summand; note: error in eq. (4.7), Jelenic 1999: the first summand in the parantheses should be \hat{\Lambda N} instead of \Lambda N
          auxmatrix2.Multiply(cn,r_s_hat);
          LARGEROTATIONS::computespin(auxmatrix1,stressn);
          auxmatrix2 -= auxmatrix1;       // auxmatrix2: term in parantheses

          auxmatrix1.Multiply(auxmatrix2,Itilde[nodej]);  // term in parantheses is the same as in upper right block but with opposite sign (note '-=' below)

          auxmatrix3.Multiply(r_s_hat,auxmatrix1);
          auxmatrix3.Scale(I_i[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofperclnode*nnodecl+dofpertriadnode*nodei+i,dofperclnode*nnodecl+dofpertriadnode*nodej+j) -= auxmatrix3(i,j)*jacobiGPelastf_[numgp]*wgt;
        }
      }

    } // if (stiffmatrix != NULL)
#endif


    // add elastic energy from forces at this GP
    for(int dim=0; dim<3; dim++)
    {
      Eint_ += 0.5 * FADUTILS::CastToDouble(Gamma(dim)) * FADUTILS::CastToDouble(stressN(dim)) *jacobiGPelastf_[numgp]*wgt;
    }

  }//for(int numgp=0; numgp < gausspoints_elast_force.nquad; numgp++)


  //************************* residual and stiffmatrix contributions from moments ***********************

  // get integration points for elasticity
  DRT::UTILS::IntegrationPoints1D gausspoints_elast_moment(MyGaussRule(res_elastic_moment));

  // reuse variables for individual shape functions and resize to new numgp
  I_i.resize(gausspoints_elast_moment.nquad);
  I_i_xi.resize(gausspoints_elast_moment.nquad);

  // evaluate all shape functions and derivatives with respect to element parameter xi at all specified Gauss points
  EvaluateShapeFunctionsAndDerivsAllGPs<nnodetriad,1>(gausspoints_elast_moment,I_i,I_i_xi,this->Shape());

  // reset norm of maximal bending curvature
  Kmax_=0.0;

  // Loop through all GP and calculate their contribution to the forcevector and stiffnessmatrix
  for(int numgp=0; numgp<gausspoints_elast_moment.nquad; numgp++)
  {
    // weight of GP in parameter space
    const double wgt = gausspoints_elast_moment.qwgt[numgp];

    Calc_Psi_l<nnodetriad,FADordouble>(Psi_li, I_i[numgp], Psi_l);
    Calc_Psi_l_s<nnodetriad,FADordouble>(Psi_li, I_i_xi[numgp], jacobiGPelastm_[numgp], Psi_l_s);
    Calc_Lambda<FADordouble>(Psi_l,Q_r,Lambda);

    // compute material curvature K
    computeK<FADordouble>(Psi_l,Psi_l_s,KrefGP_[numgp],K);

    // determine norm of maximal bending curvature at this GP and store in class variable if needed
    double Kmax = std::sqrt(FADUTILS::CastToDouble(K(1))*FADUTILS::CastToDouble(K(1)) + FADUTILS::CastToDouble(K(2))*FADUTILS::CastToDouble(K(2)) );
    if(Kmax > Kmax_)
      Kmax_=Kmax;

    // compute material stresses by multiplying curvature with constitutive matrix
    stressM.Multiply(CM,K);

    /* compute spatial stresses and constitutive matrix from material ones according to Jelenic 1999, page 148, paragraph
     * between (2.22) and (2.23) and Romero 2004, (3.10)*/
    pushforward<FADordouble>(Lambda,stressM,CM,stressm,cm);

    /* computation of internal forces according to Jelenic 1999, eq. (4.3); computation split up with respect
     * to single blocks of matrix in eq. (4.3)*/
    for (unsigned int node=0; node<nnodecl; ++node)
    {
      // lower right block
      for (int i=0; i<3; ++i)
        f_int(dofpercombinode*node+3+i) += I_i_xi[numgp](node)*stressm(i)*wgt;
    }
    for (unsigned int node=nnodecl; node<nnodetriad; ++node)    // this loop is only entered in case of nnodetriad>nnodecl
    {
      // lower right block
      for (int i=0; i<3; ++i)
        f_int(dofperclnode*nnodecl+dofpertriadnode*node+i) += I_i_xi[numgp](node)*stressm(i)*wgt;
    }

#ifndef BEAM3RAUTOMATICDIFF

    if (stiffmatrix != NULL)
    {
      computeItilde<nnodetriad>(Psi_l,Itilde,Phi_IJ,Lambda_r,Psi_li,I_i[numgp]);
      computeItildeprime<nnodetriad,double>(Psi_l,Psi_l_s,Itildeprime,Phi_IJ,Lambda_r,Psi_li,I_i[numgp],I_i_xi[numgp],jacobiGPelastm_[numgp]);

      /* computation of stiffness matrix according to Jelenic 1999, eq. (4.7)*/

      // auxiliary variables for storing intermediate matrices in computation of entries of stiffness matrix
      LINALG::Matrix<3,3> auxmatrix1;
      LINALG::Matrix<3,3> auxmatrix2;

      for(unsigned int nodei=0; nodei<nnodecl; nodei++)
      {
        for(unsigned int nodej=0; nodej<nnodecl; nodej++)
        {
          // lower right block
          // first summand
          auxmatrix1.Multiply(cm,Itildeprime[nodej]);
          auxmatrix1.Scale(I_i_xi[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofpercombinode*nodei+3+i,dofpercombinode*nodej+3+j) += auxmatrix1(i,j)*wgt;

          // second summand
          LARGEROTATIONS::computespin(auxmatrix2,stressm);
          auxmatrix1.Multiply(auxmatrix2,Itilde[nodej]);
          auxmatrix1.Scale(I_i_xi[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofpercombinode*nodei+3+i,dofpercombinode*nodej+3+j) -= auxmatrix1(i,j)*wgt;
        }
        for(unsigned int nodej=nnodecl; nodej<nnodetriad; nodej++)    // this loop is only entered in case of nnodetriad>nnodecl
        {
          // lower right block
          // first summand
          auxmatrix1.Multiply(cm,Itildeprime[nodej]);
          auxmatrix1.Scale(I_i_xi[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofpercombinode*nodei+3+i,dofperclnode*nnodecl+dofpertriadnode*nodej+j) += auxmatrix1(i,j)*wgt;

          // second summand
          LARGEROTATIONS::computespin(auxmatrix2,stressm);
          auxmatrix1.Multiply(auxmatrix2,Itilde[nodej]);
          auxmatrix1.Scale(I_i_xi[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofpercombinode*nodei+3+i,dofperclnode*nnodecl+dofpertriadnode*nodej+j) -= auxmatrix1(i,j)*wgt;
        }
      }

      for(unsigned int nodei=nnodecl; nodei<nnodetriad; nodei++)    // this loop is only entered in case of nnodetriad>nnodecl
      {
        for(unsigned int nodej=0; nodej<nnodecl; nodej++)
        {
          // lower right block
          // first summand
          auxmatrix1.Multiply(cm,Itildeprime[nodej]);
          auxmatrix1.Scale(I_i_xi[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofperclnode*nnodecl+dofpertriadnode*nodei+i,dofpercombinode*nodej+3+j) += auxmatrix1(i,j)*wgt;

          // second summand
          LARGEROTATIONS::computespin(auxmatrix2,stressm);
          auxmatrix1.Multiply(auxmatrix2,Itilde[nodej]);
          auxmatrix1.Scale(I_i_xi[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofperclnode*nnodecl+dofpertriadnode*nodei+i,dofpercombinode*nodej+3+j) -= auxmatrix1(i,j)*wgt;
        }
        for(unsigned int nodej=nnodecl; nodej<nnodetriad; nodej++)
        {
          // lower right block
          // first summand
          auxmatrix1.Multiply(cm,Itildeprime[nodej]);
          auxmatrix1.Scale(I_i_xi[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofperclnode*nnodecl+dofpertriadnode*nodei+i,dofperclnode*nnodecl+dofpertriadnode*nodej+j) += auxmatrix1(i,j)*wgt;

          // second summand
          LARGEROTATIONS::computespin(auxmatrix2,stressm);
          auxmatrix1.Multiply(auxmatrix2,Itilde[nodej]);
          auxmatrix1.Scale(I_i_xi[numgp](nodei));
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofperclnode*nnodecl+dofpertriadnode*nodei+i,dofperclnode*nnodecl+dofpertriadnode*nodej+j) -= auxmatrix1(i,j)*wgt;
        }
      }

    }//if (stiffmatrix != NULL)
#endif

    // add elastic energy from moments at this GP
    for(int dim=0; dim<3; dim++)
    {
      Eint_ += 0.5 * FADUTILS::CastToDouble(K(dim)) * FADUTILS::CastToDouble(stressM(dim)) *jacobiGPelastm_[numgp]*wgt;
    }

  }//for(int numgp=0; numgp < gausspoints_elast_moment.nquad; numgp++)

  if(force!=NULL)
  {
    for (unsigned int i=0; i<dofperclnode*nnodecl+dofpertriadnode*nnodetriad; i++)
    {
      (*force)(i)= FADUTILS::CastToDouble(f_int(i));
    }
  }

#ifdef BEAM3RAUTOMATICDIFF

  if (stiffmatrix != NULL)
  {
    // compute stiffness matrix with FAD
    for (unsigned int i=0; i<dofperclnode*nnodecl+dofpertriadnode*nnodetriad; i++)
    {
      for (unsigned int j=0; j<dofperclnode*nnodecl+dofpertriadnode*nnodetriad; j++)
      {
        (*stiffmatrix)(i,j)=f_int(i).dx(j);
      }
    }

#ifdef MULTIPLICATIVEUPDATES

    /* we need to transform the stiffmatrix because its entries are derivatives with respect to additive rotational increments
     * we want a stiffmatrix containing derivatives with respect to multiplicative rotational increments
     * therefore apply a trafo matrix to all those 3x3 blocks in stiffmatrix which correspond to derivation with respect to rotational DOFs
     * the trafo matrix is simply the T-Matrix (see Jelenic1999, (2.4)): \Delta_{mult} \vec \theta_{inode} = \mat T(\vec \theta_{inode} * \Delta_{addit} \vec \theta_{inode}*/

    LINALG::TMatrix<FAD,3,3> tempmat(true);
    LINALG::TMatrix<FAD,3,3> newstiffmat(true);
    LINALG::TMatrix<FAD,3,3> Tmat(true);
    LINALG::TMatrix<FAD,3,1> deltatheta_j(true);
    LINALG::TMatrix<FAD,3,1> theta_totlag_j(true);
    LINALG::TMatrix<FAD,4,1> Q_i_old(true);
    LINALG::TMatrix<FAD,4,1> deltaQ(true);
    LINALG::TMatrix<FAD,4,1> Q_i(true);

    for (unsigned int jnode=0; jnode<nnodecl; jnode++)
    {
      theta_totlag_j.Clear();

      for (int dim=0; dim<3; ++dim)
        deltatheta_j(dim) = dispthetanewnode_[jnode](dim) - dispthetaoldnode_[jnode](dim);

      Q_i_old.Clear();
      for (int i=0; i<4; ++i)
        Q_i_old(i) = Qoldnode_[jnode](i);

      // compute quaternion from rotation angle relative to last configuration
      LARGEROTATIONS::angletoquaternion(deltatheta_j,deltaQ);

      // multiply relative rotation with rotation in last configuration to get rotation in new configuration
      LARGEROTATIONS::quaternionproduct(Q_i_old,deltaQ,Q_i);

      // compute physical total angle theta_totlag
      LARGEROTATIONS::quaterniontoangle(Q_i,theta_totlag_j);

      // compute Tmatrix of theta_totlag_i
      Tmat = LARGEROTATIONS::Tmatrix(theta_totlag_j);

      for (unsigned int inode=0; inode<nnodecl; inode++)
      {
        // block1: derivative of nodal positions with respect to theta (rotational DOFs)
        for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
            tempmat(i,j) = (*stiffmatrix)(dofpercombinode*inode+i,dofpercombinode*jnode+3+j);

        newstiffmat.Clear();
        newstiffmat.MultiplyNN(tempmat,Tmat);

        for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
            (*stiffmatrix)(dofpercombinode*inode+i,dofpercombinode*jnode+3+j) = newstiffmat(i,j).val();

        // block2: derivative of nodal theta with respect to theta (rotational DOFs)
        for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
            tempmat(i,j) = (*stiffmatrix)(dofpercombinode*inode+3+i,dofpercombinode*jnode+3+j);

        newstiffmat.Clear();
        newstiffmat.MultiplyNN(tempmat,Tmat);

        for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
            (*stiffmatrix)(dofpercombinode*inode+3+i,dofpercombinode*jnode+3+j) = newstiffmat(i,j).val();

        // block3: derivative of nodal tangents with respect to theta (rotational DOFs)
        if(centerline_hermite_)
        {
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              tempmat(i,j) = (*stiffmatrix)(dofpercombinode*inode+6+i,dofpercombinode*jnode+3+j);

          newstiffmat.Clear();
          newstiffmat.MultiplyNN(tempmat,Tmat);

          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofpercombinode*inode+6+i,dofpercombinode*jnode+3+j) = newstiffmat(i,j).val();
        }

      }
      for (unsigned int inode=nnodecl; inode<nnodetriad; inode++)    // this loop is only entered in case of nnodetriad>nnodecl
      {
        // block2: derivative of nodal theta with respect to theta (rotational DOFs)
        for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
            tempmat(i,j) = (*stiffmatrix)(dofperclnode*nnodecl+dofpertriadnode*inode+i,dofpercombinode*jnode+3+j);

        newstiffmat.Clear();
        newstiffmat.MultiplyNN(tempmat,Tmat);

        for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
            (*stiffmatrix)(dofperclnode*nnodecl+dofpertriadnode*inode+i,dofpercombinode*jnode+3+j) = newstiffmat(i,j).val();

      }
    }


    for (unsigned int jnode=nnodecl; jnode<nnodetriad; jnode++)    // this loop is only entered in case of nnodetriad>nnodecl
    {
      theta_totlag_j.Clear();

      for (int dim=0; dim<3; ++dim)
        deltatheta_j(dim) = dispthetanewnode_[jnode](dim) - dispthetaoldnode_[jnode](dim);

      Q_i_old.Clear();
      for (int i=0; i<4; ++i)
        Q_i_old(i) = Qoldnode_[jnode](i);

      // compute quaternion from rotation angle relative to last configuration
      LARGEROTATIONS::angletoquaternion(deltatheta_j,deltaQ);

      // multiply relative rotation with rotation in last configuration to get rotation in new configuration
      LARGEROTATIONS::quaternionproduct(Q_i_old,deltaQ,Q_i);

      // compute physical total angle theta_totlag
      LARGEROTATIONS::quaterniontoangle(Q_i,theta_totlag_j);

      // compute Tmatrix of theta_totlag_i
      Tmat = LARGEROTATIONS::Tmatrix(theta_totlag_j);

      for (unsigned int inode=0; inode<nnodecl; inode++)
      {
        // block1: derivative of nodal positions with respect to theta (rotational DOFs)
        for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
            tempmat(i,j) = (*stiffmatrix)(dofpercombinode*inode+i,dofperclnode*nnodecl+dofpertriadnode*jnode+j);

        newstiffmat.Clear();
        newstiffmat.MultiplyNN(tempmat,Tmat);

        for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
            (*stiffmatrix)(dofpercombinode*inode+i,dofperclnode*nnodecl+dofpertriadnode*jnode+j) = newstiffmat(i,j).val();

        // block2: derivative of nodal theta with respect to theta (rotational DOFs)
        for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
            tempmat(i,j) = (*stiffmatrix)(dofpercombinode*inode+3+i,dofperclnode*nnodecl+dofpertriadnode*jnode+j);

        newstiffmat.Clear();
        newstiffmat.MultiplyNN(tempmat,Tmat);

        for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
            (*stiffmatrix)(dofpercombinode*inode+3+i,dofperclnode*nnodecl+dofpertriadnode*jnode+j) = newstiffmat(i,j).val();

        // block3: derivative of nodal tangents with respect to theta (rotational DOFs)
        if(centerline_hermite_)
        {
          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              tempmat(i,j) = (*stiffmatrix)(dofpercombinode*inode+6+i,dofperclnode*nnodecl+dofpertriadnode*jnode+j);

          newstiffmat.Clear();
          newstiffmat.MultiplyNN(tempmat,Tmat);

          for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
              (*stiffmatrix)(dofpercombinode*inode+6+i,dofperclnode*nnodecl+dofpertriadnode*jnode+j) = newstiffmat(i,j).val();
        }

      }
      for (unsigned int inode=nnodecl; inode<nnodetriad; inode++)
      {
        // block2: derivative of nodal theta with respect to theta (rotational DOFs)
        for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
            tempmat(i,j) = (*stiffmatrix)(dofperclnode*nnodecl+dofpertriadnode*inode+i,dofperclnode*nnodecl+dofpertriadnode*jnode+j);

        newstiffmat.Clear();
        newstiffmat.MultiplyNN(tempmat,Tmat);

        for (int i=0; i<3; ++i)
          for (int j=0; j<3; ++j)
            (*stiffmatrix)(dofperclnode*nnodecl+dofpertriadnode*inode+i,dofperclnode*nnodecl+dofpertriadnode*jnode+j) = newstiffmat(i,j).val();

      }
    }

#endif
  }
#endif


  /******************************* inertia: compute fint and massmatrix ********************************
   *****************************************************************************************************/

  // calculation of inertia forces/moments and massmatrix; in case of Statmech, a dummy massmatrix is computed below
  if (massmatrix != NULL and inertia_force != NULL and (!needstatmech_) )
  {
#ifdef BEAM3RAUTOMATICDIFF
#ifndef MULTIPLICATIVEUPDATES
      dserror("beam3r: for dynamic simulations in combination with a FAD stiffness matrix, please define flag MULTIPLICATIVEUPDATES because the implemented GenAlpha Lie group time integration is based on multiplicative updates!");
#endif
#endif

    /* calculation of mass matrix: According to the paper of Jelenic and Crisfield "Geometrically exact 3D beam theory: implementation of a strain-invariant
     * finite element for statics and dynamics", 1999, page 146, a time integration scheme that delivers angular velocities and angular accelerations as
     * needed for the inertia terms of geometrically exact beams has to be based on multiplicative rotation angle increments between two successive time
     * steps. Since BACI does all displacement updates in an additive manner, the global vector of rotational displacements has no physical meaning and,
     * consequently the global velocity and acceleration vectors resulting from the BACI time integration schemes have no physical meaning, too. Therefore,
     * a mass matrix in combination with this global acceleration vector is meaningless from a physical point of view. For these reasons, we have to apply
     * our own time integration scheme at element level. Up to now, the only implemented integration scheme is the gen-alpha Lie group time integration
     * according to [Arnold, Brüls (2007)], [Brüls, Cardona, 2010] and [Brüls, Cardona, Arnold (2012)] in combination with a constdisvelacc predictor. (Christoph Meier, 04.14)*/

    const double beta = params.get<double>("rot_beta",1000);
    const double gamma = params.get<double>("rot_gamma",1000);
    const double alpha_f = params.get<double>("rot_alphaf",1000);
    const double alpha_m = params.get<double>("rot_alpham",1000);
    const double dt = params.get<double>("delta time",1000);
    const bool materialintegration=true;        // TODO unused?
    const double diff_factor_vel = gamma/(beta*dt);
    const double diff_factor_acc = (1.0-alpha_m)/(beta*dt*dt*(1.0-alpha_f));

    LINALG::Matrix<3,3> Lambdanewmass(true);
    LINALG::Matrix<3,3> Lambdaconvmass(true);

    // get the material law
    Teuchos::RCP<const MAT::Material> currmat = Material();
    double rho = 0.0;

    // assignment of material parameters; only St.Venant material is accepted for this beam
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

    /* tensor of mass moments of inertia and cross-section value. These values are used in order to artificially scale
     * the the translational and rotational inertia terms with given input parameters if necessary*/
    LINALG::Matrix<3,3> Jp(true);
    Jp(0,0)=inertscalerot1_*(Iyy_+Izz_);
    Jp(1,1)=inertscalerot2_*Iyy_;
    Jp(2,2)=inertscalerot2_*Izz_;
    Jp.Scale(rho);

    const double scaledcrosssec=inertscaletrans_*crosssec_;

    // get integration scheme for inertia forces and mass matrix
    DRT::UTILS::IntegrationPoints1D gausspoints_mass(MyGaussRule(res_inertia));
    // reuse variables for individual shape functions and resize to new numgp
    I_i.resize(gausspoints_mass.nquad);
    H_i.resize(gausspoints_mass.nquad);

    // evaluate all shape functions at all specified Gauss points
    EvaluateShapeFunctionsAllGPs<nnodetriad,1>(gausspoints_mass,I_i,this->Shape());
    EvaluateShapeFunctionsAllGPs<nnodecl,vpernode>(gausspoints_mass,H_i,this->Shape());

    // Calculate current centerline position at gauss points (needed for element intern time integration)
    for (int gp=0; gp<gausspoints_mass.nquad; gp++)//loop through Gauss points
      Calc_r<nnodecl,vpernode,double>(FADUTILS::CastToDouble<FADordouble,3*vpernode*nnodecl,1>(disp_totlag_centerline),H_i[gp],rnewGPmass_[gp]);

    Ekin_=0.0;
    L_=0.0;
    P_=0.0;

    for (int gp=0; gp<gausspoints_mass.nquad; gp++)//loop through Gauss points
    {
      // weight of GP in parameter space
      const double wgtmass = gausspoints_mass.qwgt[gp];

      LINALG::Matrix<3,3> Jp_bar(Jp);
      Jp_bar.Scale(diff_factor_acc);

      LINALG::Matrix<3,1> dL(true);

      // update quaternions at GPs for exact Gauss quadrature
      Calc_Psi_l<nnodetriad,FADordouble>(Psi_li, I_i[gp], Psi_l);
      Calc_Qgauss<double>(FADUTILS::CastToDouble<FADordouble,3,1>(Psi_l),FADUTILS::CastToDouble<FADordouble,4,1>(Q_r),QnewGPmass_[gp]);
      computeItilde<nnodetriad>(Psi_l,Itilde,Phi_IJ,Lambda_r,Psi_li,I_i[gp]);

      Lambdanewmass.Clear();
      Lambdaconvmass.Clear();
      // compute current and old triad at Gauss point
      LARGEROTATIONS::quaterniontotriad<double>(QnewGPmass_[gp],Lambdanewmass);
      LARGEROTATIONS::quaterniontotriad<double>(QconvGPmass_[gp],Lambdaconvmass);

      // rotation between last converged position and current position expressed as a quaternion
      LINALG::Matrix<4,1>  deltaQ(true);
      LARGEROTATIONS::quaternionproduct<double>(LARGEROTATIONS::inversequaternion<double>(QconvGPmass_[gp]),QnewGPmass_[gp],deltaQ);

      // spatial rotation between last converged position and current position expressed as a three element rotation vector
      LINALG::Matrix<3,1> deltatheta(true);
      LARGEROTATIONS::quaterniontoangle<double>(deltaQ,deltatheta);

      // compute material counterparts of spatial vectors
      LINALG::Matrix<3,1> deltaTHETA(true);
      LINALG::Matrix<3,1> Wconvmass(true);
      LINALG::Matrix<3,1> Wnewmass(true);
      LINALG::Matrix<3,1> Aconvmass(true);
      LINALG::Matrix<3,1> Anewmass(true);
      LINALG::Matrix<3,1> Amodconvmass(true);
      LINALG::Matrix<3,1> Amodnewmass(true);
      deltaTHETA.MultiplyTN(Lambdanewmass,deltatheta);
      Wconvmass.MultiplyTN(Lambdaconvmass,wconvGPmass_[gp]);
      Aconvmass.MultiplyTN(Lambdaconvmass,aconvGPmass_[gp]);
      Amodconvmass.MultiplyTN(Lambdaconvmass,amodconvGPmass_[gp]);

      /* update angular velocities and accelerations according to Newmark time integration scheme either in
       * material description (see Jelenic, 1999, p. 146, equations (2.8) and (2.9)) or in spatial description
       * (for testing purposes, not recommended by Jelenic). The corresponding equations are adapted according to
       * the gen-alpha Lie group time integration scheme proposed in [Arnold, Brüls (2007)], [Brüls, Cardona, 2010]
       * and [Brüls, Cardona, Arnold (2012)]. In the predictor step of the time integration the following
       * formulas automatically deliver a constant displacement (deltatheta=0), consistent velocity and consistent acceleration
       * predictor. This fact has to be reflected in a consistent manner by the choice of the predictor in the input file*/
      if (materialintegration)
      {
        for (int i=0;i<3;i++)
        {
          Anewmass(i)=   (1.0-alpha_m)/(beta*dt*dt*(1.0-alpha_f))*deltaTHETA(i)-(1.0-alpha_m)/(beta*dt*(1.0-alpha_f))*Wconvmass(i)      \
                        -alpha_f/(1.0-alpha_f)*Aconvmass(i)+(alpha_m/(1.0-alpha_f)-(0.5-beta)*(1.0-alpha_m)/(beta*(1.0-alpha_f)))*Amodconvmass(i);

          Wnewmass(i)=gamma/(beta*dt)*deltaTHETA(i)+(1-gamma/beta)*Wconvmass(i)+dt*(1-gamma/(2*beta))*Amodconvmass(i);

          Amodnewmass(i)=1.0/(1.0-alpha_m)*((1.0-alpha_f)*Anewmass(i) + alpha_f*Aconvmass(i) - alpha_m*Amodconvmass(i));
        }
        wnewGPmass_[gp].Multiply(Lambdanewmass,Wnewmass);
        anewGPmass_[gp].Multiply(Lambdanewmass,Anewmass);
        amodnewGPmass_[gp].Multiply(Lambdanewmass,Amodnewmass);
      }
      else
      {
        for (int i=0;i<3;i++)
        {
          wnewGPmass_[gp](i)=gamma/(beta*dt)*deltatheta(i)+(1-gamma/beta)*wconvGPmass_[gp](i)+dt*(1-gamma/(2*beta))*amodconvGPmass_[gp](i);

          anewGPmass_[gp](i)=   (1.0-alpha_m)/(beta*dt*dt*(1.0-alpha_f))*deltatheta(i)-(1.0-alpha_m)/(beta*dt*(1.0-alpha_f))*wconvGPmass_[gp](i)      \
                        -alpha_f/(1.0-alpha_f)*aconvGPmass_[gp](i)+(alpha_m/(1.0-alpha_f)-(0.5-beta)*(1.0-alpha_m)/(beta*(1.0-alpha_f)))*amodconvGPmass_[gp](i);

          amodnewGPmass_[gp](i)=1.0/(1.0-alpha_m)*((1.0-alpha_f)*anewGPmass_[gp](i) + alpha_f*aconvGPmass_[gp](i) - alpha_m*amodconvGPmass_[gp](i) );
        }
        Wnewmass.MultiplyTN(Lambdanewmass,wnewGPmass_[gp]);
        Anewmass.MultiplyTN(Lambdanewmass,anewGPmass_[gp]);
        Amodnewmass.MultiplyTN(Lambdanewmass,amodnewGPmass_[gp]);
      }

      LINALG::Matrix<3,1> deltar(true);
      for (int i=0;i<3;i++)
      {
        deltar(i)=rnewGPmass_[gp](i)-rconvGPmass_[gp](i);
      }
      for (int i=0;i<3;i++)
      {
        rttnewGPmass_[gp](i)=   (1.0-alpha_m)/(beta*dt*dt*(1.0-alpha_f))*deltar(i)-(1.0-alpha_m)/(beta*dt*(1.0-alpha_f))*rtconvGPmass_[gp](i)      \
                      -alpha_f/(1.0-alpha_f)*rttconvGPmass_[gp](i)+(alpha_m/(1.0-alpha_f)-(0.5-beta)*(1.0-alpha_m)/(beta*(1.0-alpha_f)))*rttmodconvGPmass_[gp](i);

        rtnewGPmass_[gp](i)=gamma/(beta*dt)*deltar(i)+(1-gamma/beta)*rtconvGPmass_[gp](i)+dt*(1-gamma/(2*beta))*rttmodconvGPmass_[gp](i);

        rttmodnewGPmass_[gp](i)=1.0/(1.0-alpha_m)*((1.0-alpha_f)*rttnewGPmass_[gp](i) + alpha_f*rttconvGPmass_[gp](i) - alpha_m*rttmodconvGPmass_[gp](i) );
      }

      // spin matrix of the material angular velocity, i.e. S(W)
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

      r_tt = rttnewGPmass_[gp];
      r_t = rtnewGPmass_[gp];
      r = rnewGPmass_[gp];

      LINALG::Matrix<3,3> S_r(true);
      LARGEROTATIONS::computespin<double>(S_r,r);
      dL.Multiply(S_r,r_t);
      dL.Scale(rho*scaledcrosssec);
      LINALG::Matrix<3,1> Lambdanewmass_Jp_Wnewmass(true);
      Lambdanewmass_Jp_Wnewmass.Multiply(Lambdanewmass,Jp_Wnewmass);
      dL.Update(1.0,Lambdanewmass_Jp_Wnewmass,1.0);
      for (int i=0;i<3;i++)
      {
        L_(i)+=wgtmass*jacobiGPmass_[gp]*dL(i);
        P_(i)+=wgtmass*jacobiGPmass_[gp]*rho*scaledcrosssec*r_t(i);
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

      // inertia forces
      for (int i=0; i<3; i++)
      {
        for (unsigned int node=0; node<nnodecl; node++)
        {
          // translational contribution
          (*inertia_force)(dofpercombinode*node+i) += jacobiGPmass_[gp]*wgtmass*rho*scaledcrosssec*H_i[gp](vpernode*node)*r_tt(i);
          if (centerline_hermite_)
            (*inertia_force)(dofpercombinode*node+6+i) += jacobiGPmass_[gp]*wgtmass*rho*scaledcrosssec*H_i[gp](vpernode*node+1)*r_tt(i);
          //rotational contribution
          (*inertia_force)(dofpercombinode*node+3+i) += jacobiGPmass_[gp]*wgtmass*I_i[gp](node)*Pi_t(i);
        }
        for (unsigned int node=nnodecl; node<nnodetriad; node++)    // this loop is only entered in case of nnodetriad>nnodecl
        {
          // rotational contribution
          (*inertia_force)(dofperclnode*nnodecl+dofpertriadnode*node+i) += jacobiGPmass_[gp]*wgtmass*I_i[gp](node)*Pi_t(i);
        }
      }

      // linearization of inertia forces: massmatrix
      for (unsigned int jnode=0; jnode<nnodecl; jnode++)
      {
        // translational contribution
        for (unsigned int inode=0; inode<nnodecl; inode++)
          for (int k=0;k<3;k++)
          {
            (*massmatrix)(dofpercombinode*inode+k,dofpercombinode*jnode+k) += diff_factor_acc*jacobiGPmass_[gp]*wgtmass*rho*scaledcrosssec*H_i[gp](vpernode*inode)*H_i[gp](vpernode*jnode);
            if (centerline_hermite_)
            {
              (*massmatrix)(dofpercombinode*inode+6+k,dofpercombinode*jnode+6+k) += diff_factor_acc*jacobiGPmass_[gp]*wgtmass*rho*scaledcrosssec*H_i[gp](vpernode*inode+1)*H_i[gp](vpernode*jnode+1);
              (*massmatrix)(dofpercombinode*inode+k,dofpercombinode*jnode+6+k) += diff_factor_acc*jacobiGPmass_[gp]*wgtmass*rho*scaledcrosssec*H_i[gp](vpernode*inode)*H_i[gp](vpernode*jnode+1);
              (*massmatrix)(dofpercombinode*inode+6+k,dofpercombinode*jnode+k) += diff_factor_acc*jacobiGPmass_[gp]*wgtmass*rho*scaledcrosssec*H_i[gp](vpernode*inode+1)*H_i[gp](vpernode*jnode);
            }
          }

        // rotational contribution
        LINALG::Matrix<3,3> auxmatrix2(true);
        auxmatrix2.Multiply(auxmatrix1,Itilde[jnode]);
        for (unsigned int inode=0; inode<nnodecl; inode++)
        {
          for (int i=0; i<3; i++)
            for (int j=0; j<3; j++)
              (*massmatrix)(dofpercombinode*inode+3+i,dofpercombinode*jnode+3+j) += jacobiGPmass_[gp]*wgtmass*I_i[gp](inode)*auxmatrix2(i,j);
        }
        for (unsigned int inode=nnodecl; inode<nnodetriad; inode++)    // this loop is only entered in case of nnodetriad>nnodecl
        {
          for (int i=0; i<3; i++)
            for (int j=0; j<3; j++)
              (*massmatrix)(dofperclnode*nnodecl+dofpertriadnode*inode+i,dofpercombinode*jnode+3+j) += jacobiGPmass_[gp]*wgtmass*I_i[gp](inode)*auxmatrix2(i,j);
        }
      }
      for (unsigned int jnode=nnodecl; jnode<nnodetriad; ++jnode)    // this loop is only entered in case of nnodetriad>nnodecl
      {
        // rotational contribution
        LINALG::Matrix<3,3> auxmatrix2(true);
        auxmatrix2.Multiply(auxmatrix1,Itilde[jnode]);
        for (unsigned int inode=0; inode<nnodecl; inode++)
        {
          for (int i=0; i<3; i++)
            for (int j=0; j<3; j++)
              (*massmatrix)(dofpercombinode*inode+3+i,dofperclnode*nnodecl+dofpertriadnode*jnode+j) += jacobiGPmass_[gp]*wgtmass*I_i[gp](inode)*auxmatrix2(i,j);
        }
        for (unsigned int inode=nnodecl; inode<nnodetriad; inode++)    // this loop is only entered in case of nnodetriad>nnodecl
        {
          for (int i=0; i<3; i++)
            for (int j=0; j<3; j++)
              (*massmatrix)(dofperclnode*nnodecl+dofpertriadnode*inode+i,dofperclnode*nnodecl+dofpertriadnode*jnode+j) += jacobiGPmass_[gp]*wgtmass*I_i[gp](inode)*auxmatrix2(i,j);
        }
      }

      // Calculation of kinetic energy
      LINALG::Matrix<1,1> ekinrot(true);
      LINALG::Matrix<1,1> ekintrans(true);
      ekinrot.MultiplyTN(Wnewmass,Jp_Wnewmass);
      ekintrans.MultiplyTN(r_t,r_t);
      Ekin_+=0.5*(ekinrot.Norm2() + rho*scaledcrosssec*ekintrans.Norm2())*jacobiGPmass_[gp]*wgtmass;
      Ekintorsion_+=0.5* Wnewmass(0)*Jp_Wnewmass(0)*jacobiGPmass_[gp]*wgtmass;
      Ekinbending_+=0.5* Wnewmass(1)*Jp_Wnewmass(1)*jacobiGPmass_[gp]*wgtmass;
      Ekinbending_+=0.5* Wnewmass(2)*Jp_Wnewmass(2)*jacobiGPmass_[gp]*wgtmass;
      Ekintrans_+=0.5*rho*scaledcrosssec*ekintrans.Norm2()*jacobiGPmass_[gp]*wgtmass;

      Jp_Wnewmass.Multiply(Jp,Wnewmass);
    }//for (int gp=0; gp<gausspoints_mass.nquad; gp++)
  }//if (massmatrix != NULL and inertia_force != NULL)


  /**************** damping and stochastic forces: compute fint and stiffmatrix ************************
   *****************************************************************************************************/

  // do the following only for applications of statistical mechanics
  if(needstatmech_)
  {
    // so far, centerline and triad field is interpolated with linear Lagrange polynomials, i.e. only use I_i and nnodetriad (=2) in the following
    if(centerline_hermite_) dserror("Hermite interpolation of centerline not implemented yet for Statmech applications");
    if(nnodetriad!=nnodecl) dserror("you should not be here, different interpolation of centerline and triad field is not implemented yet for Statmech applications");

    // in statistical mechanics simulations, a deletion influenced by the values of the internal force vector might occur
    // TODO check: do we still need this in Statmech?
    if( (gausspoints_elast_force.nquad==1 || gausspoints_elast_moment.nquad==1) && (params.get<std::string>("internalforces","no")=="yes") && force != NULL)
    {
      eps_ = FADUTILS::CastToDouble(Gamma(0));
      f_ = *force;
      for (int i=0; i<3; ++i)
        Ngp_(i) = FADUTILS::CastToDouble(stressN(i));
    }

    //This is a dummy mass matrix which is necessary for statmech simulations
    if(massmatrix != NULL)
    {
      for (unsigned int i=0; i<6*nnodetriad; i++)
        (*massmatrix)(i,i) = 1;
    }


    // check whether random numbers vector exists
    // TODO this check is a leftover from historic code; check: necessary? best location/ best way to do this?
    if(params.get<Teuchos::RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null) != Teuchos::null)
    {
      // get integration scheme for damping and stochastic terms
      DRT::UTILS::GaussRule1D gaussrule_damp_stoch = MyGaussRule(params,res_damp_stoch);
      DRT::UTILS::IntegrationPoints1D gausspoints_damp_stoch(gaussrule_damp_stoch);

      /* the following is pre-calculated for all Statmech specific methods, i.e. damping and stochastic terms
       * vector with nnode elements, who represent the 3x3-matrix-shaped interpolation function \tilde{I}^nnode at all Gauss points according to (3.19), Jelenic 1999*/
      std::vector<std::vector<LINALG::Matrix<3,3> > > Itildedamping_allGP(gausspoints_damp_stoch.nquad, std::vector<LINALG::Matrix<3,3> >(nnodetriad) );

      // if case of frictionmodel_isotropiclumped, the following is done in CalcBrownian() -> skip this block
      if(gaussrule_damp_stoch != DRT::UTILS::intrule_line_lobatto2point)
      {
        // reuse variables for individual shape functions and resize to new numgp
        I_i.resize(gausspoints_damp_stoch.nquad);

        //evaluate all shape functions at all specified Gauss points
        EvaluateShapeFunctionsAllGPs<nnodetriad,1>(gausspoints_damp_stoch,I_i,this->Shape());

        dsassert(gaussrule_damp_stoch == MyGaussRule(res_inertia), "this implementation assumes the same integration scheme for mass and statmech terms because it reuses class variables dispnewGPmass_ and QnewGPmass_");

        for (int gp=0; gp<gausspoints_damp_stoch.nquad; gp++)//loop through Gauss points
        {
          // update quaternions at GPs for exact Gauss quadrature
          Calc_Psi_l<nnodetriad,FADordouble>(Psi_li, I_i[gp], Psi_l);
          Calc_Qgauss<double>(FADUTILS::CastToDouble<FADordouble,3,1>(Psi_l),FADUTILS::CastToDouble<FADordouble,4,1>(Q_r),QnewGPmass_[gp]);
          computeItilde<nnodetriad>(Psi_l,Itilde,Phi_IJ,Lambda_r,Psi_li,I_i[gp]);

          for (unsigned int inode=0; inode<nnodetriad; ++inode)
            Itildedamping_allGP[gp][inode].Update(Itilde[inode]);
        }
      }

      /* the function CalcBrownian applies statistical forces and damping matrix according to the fluctuation dissipation theorem;
       * it is dedicated to the application of beam elements in the frame of statistical mechanics problems*/
#ifndef BEAM3RCONSTSTOCHFORCE
        CalcBrownian<nnodetriad,3,6,4>(params,vel,disp,stiffmatrix,force,I_i,Itildedamping_allGP);
#else
        CalcBrownian<nnode,3,6,3>(params,vel,disp,stiffmatrix,force,I_i,Itildedamping_allGP);
#endif
    }

  }

  return;
} // DRT::ELEMENTS::Beam3r::b3_nlnstiffmass

/*------------------------------------------------------------------------------------------------------------*
 | update (total) displacement vector and set nodal triads (as quaternions)                        grill 03/16|
 *------------------------------------------------------------------------------------------------------------*/
template<unsigned int nnodetriad, unsigned int nnodecl, unsigned int vpernode>
void DRT::ELEMENTS::Beam3r::UpdateDispTotLagAndNodalTriads(const std::vector<double>&                         disp,
                                                           LINALG::TMatrix<FADordouble,3*vpernode*nnodecl,1>& disp_totlag_centerline,
                                                           std::vector<LINALG::TMatrix<FADordouble,4,1> >&    Q_i)
{
  // nnodetriad: number of nodes used for interpolation of triad field
  // nnodecl: number of nodes used for interpolation of centerline
  // assumptions: nnodecl<=nnodetriad; centerline nodes have local ID 0...nnodecl-1
  // vpernode: number of interpolated values per centerline node (1: value (i.e. Lagrange), 2: value + derivative of value (i.e. Hermite))

  const int dofperclnode = 3*vpernode;
  const int dofpertriadnode = 3;
  const int dofpercombinode = dofperclnode+dofpertriadnode;

  // rotational displacement at a certain node between this and last iteration step
  LINALG::TMatrix<FADordouble,3,1> deltatheta;
  // rotational displacement at a certain node between this and last iteration step in quaternion form
  LINALG::TMatrix<FADordouble,4,1> deltaQ;

  // get current values of translational nodal DOFs in total Lagrangean manner (initial value + disp)
  // rotational DOFs need different handling, depending on whether FAD is used or not (see comment below)
  for (unsigned int dim=0; dim<3; ++dim)
  {
    for (unsigned int node=0; node<nnodecl; ++node)
    {
      disp_totlag_centerline(3*vpernode*node+dim) = Nodes()[node]->X()[dim] + disp[dofpercombinode*node+dim];

      // have Hermite interpolation? then update tangent DOFs as well
      if(centerline_hermite_)
        disp_totlag_centerline(3*vpernode*node+3+dim) = Trefnode_[node](dim) + disp[dofpercombinode*node+6+dim];
    }
  }

  /* get current displacement values of rotational DOFs
   * note: these values are unphysical because the last increment has been added although large rotation vectors are non-additive
   *       further processing of these values depends on choice of analytical or automatic differentiation (via FAD) */
  for(int dim=0; dim<3; dim++)
  {
    for (unsigned int node=0; node<nnodecl; ++node)
      dispthetanewnode_[node](dim) = disp[dofpercombinode*node+3+dim];

    for (unsigned int node=nnodecl; node<nnodetriad; ++node)
      dispthetanewnode_[node](dim) = disp[dofperclnode*nnodecl+dofpertriadnode*node+dim];
  }

#ifndef BEAM3RAUTOMATICDIFF
  // Compute current nodal triads
  for (unsigned int node=0; node<nnodetriad; ++node)
  {
    /* rotation increment relative to configuration in last iteration step is difference between current rotation
     * entry in displacement vector minus rotation entry in displacement vector in last iteration step*/

    /* This shift is necessary, since our beam formulation and the corresponding linearization is based on multiplicative
     * increments, while in BACI (Newtonfull()) the displacement of the last iteration and the rotation increment
     * of the current iteration are added in an additive manner. This step is reversed in the next two lines in order
     * to recover the multiplicative rotation increment between the last and the current Newton step. This also
     * means that dispthetanewnode_ has no physical meaning since it is the additive sum of non-additive rotation increments(meier, 03.2014)*/
    deltatheta  = dispthetanewnode_[node];
    deltatheta -= dispthetaoldnode_[node];

    // compute quaternion from rotation angle relative to last configuration
    LARGEROTATIONS::angletoquaternion(deltatheta,deltaQ);

    // multiply relative rotation with rotation in last configuration to get rotation in new configuration
    LARGEROTATIONS::quaternionproduct(Qoldnode_[node],deltaQ,Qnewnode_[node]);

    // renormalize quaternion to keep its absolute value one even in case of long simulations and intricate calculations
    Qnewnode_[node].Scale(1/Qnewnode_[node].Norm2());

    // copy quaternions of nodal triads to TMatrix FADordouble
    Q_i[node] = Qnewnode_[node];
  }

#else
  // set differentiation variables for FAD: translational DOFs
  for (int dim=0; dim<3; ++dim)
  {
    for (unsigned int node=0; node<nnodecl; ++node)
    {
      disp_totlag_centerline(dofperclnode*node+dim).diff(dofpercombinode*node+dim, dofperclnode*nnodecl+dofpertriadnode*nnodetriad);

      // have Hermite interpolation? then set tangent DOFs as well
      if(vpernode==2)
        disp_totlag_centerline(dofperclnode*node+3+dim).diff(dofpercombinode*node+6+dim, dofperclnode*nnodecl+dofpertriadnode*nnodetriad);
    }
  }

  // rotation vector theta at a specific node in a total Lagrangean manner (with respect to global reference coordinate system)
  std::vector<LINALG::TMatrix<FAD,3,1> > theta_totlag_i(nnodetriad);

#ifndef MULTIPLICATIVEUPDATES //use additive increments

  /* FAD will compute linearization with respect to additive increments ("normal" differentiation),
   * just as FullNewton expects, so everything is fine
   * => rotational DOFs are naively treated just like translational DOFs (initial value + disp)*/

  // compute unphysical total angle theta_totlag (large rotation pseudo-vectors are non-additive!)
  for (unsigned int node=0; node<nnodetriad; ++node)
    for (unsigned int dim=0; dim<3; ++dim)
    theta_totlag_i[node](dim) = theta0node_[node](dim) + dispthetanewnode_[node](dim);

#else //use multiplicative increments, just as in case of analytical linearization (see comments above)

  // temporary matrix of type FAD to store quaternion from last configuration
  LINALG::TMatrix<FAD,4,1> Q_i_old;

  // compute nodal quaternions based on multiplicative increments of rotational DOFs
  for (unsigned int node=0; node<nnodetriad; ++node)
  {
    for (unsigned int dim=0; dim<3; ++dim)
    {
      deltatheta(dim)  = dispthetanewnode_[node](dim);
      deltatheta(dim) -= dispthetaoldnode_[node](dim);
    }

    for (int i=0; i<4; ++i)
      Q_i_old(i) = Qoldnode_[node](i);

    // compute quaternion from rotation angle relative to last configuration
    LARGEROTATIONS::angletoquaternion(deltatheta,deltaQ);

    // multiply relative rotation with rotation in last configuration to get rotation in new configuration
    LARGEROTATIONS::quaternionproduct(Q_i_old,deltaQ,Q_i[node]);

    // compute physical total angle theta_totlag
    LARGEROTATIONS::quaterniontoangle(Q_i[node],theta_totlag_i[node]);
  }

#endif //#ifndef MULTIPLICATIVEUPDATE

  // set differentiation variables for FAD: rotational DOFs
  for (unsigned int dim=0; dim<3; ++dim)
  {
    for (unsigned int node=0; node<nnodecl; ++node)
      theta_totlag_i[node](dim).diff(dofpercombinode*node+3+dim, dofperclnode*nnodecl+dofpertriadnode*nnodetriad);

    for (unsigned int node=nnodecl; node<nnodetriad; ++node)
      theta_totlag_i[node](dim).diff(dofperclnode*nnodecl+dofpertriadnode*node+dim, dofperclnode*nnodecl+dofpertriadnode*nnodetriad);
  }

  for (unsigned int node=0; node<nnodetriad; ++node)
  {
    // compute nodal quaternions
    LARGEROTATIONS::angletoquaternion(theta_totlag_i[node],Q_i[node]);

    // renormalize quaternion to keep its absolute value one even in case of long simulations and intricate calculations.
    Q_i[node].Scale(1.0/FADUTILS::VectorNorm<4>(Q_i[node]));
  }

  // copy quaternions of nodal triads to class variable Qnewnode_
  for (unsigned int node=0; node<nnodetriad; ++node)
    for (int i=0; i<4; ++i)
      Qnewnode_[node](i)= FADUTILS::CastToDouble(Q_i[node](i));

#endif //#ifndef BEAM3RAUTOMATICDIFF

  return;
}

/*------------------------------------------------------------------------------------------------------------*
 | lump mass matrix             (private)                                                   cyron 01/08|
 *------------------------------------------------------------------------------------------------------------*/
template<unsigned int nnode>
void DRT::ELEMENTS::Beam3r::lumpmass(Epetra_SerialDenseMatrix* massmatrix)
{
  // lump mass matrix
  if (massmatrix != NULL)
  {
    // we assume #elemat2 is a square matrix
    for (int c=0; c<(*massmatrix).N(); ++c) // parse columns
    {
      double d = 0.0;
      for (int r=0; r<(*massmatrix).M(); ++r) // parse rows
      {
        d += (*massmatrix)(r,c); // accumulate row entries
        (*massmatrix)(r,c) = 0.0;
      }

      (*massmatrix)(c,c) = d; // apply sum of row entries on diagonal
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
  if(nnode>2 or centerline_hermite_)
    dserror("PTC was originally implemented for 2-noded Reissner beam element only. Check functionality for numnodes>2 and/or Hermite interpolation and extend if needed!");

  for (unsigned int node=0; node<nnode; node++)
  {
    //computing angle increment from current position in comparison with last converged position for damping
    LINALG::Matrix<4,1> deltaQ;
    LARGEROTATIONS::quaternionproduct(LARGEROTATIONS::inversequaternion(Qconvnode_[node]),Qnewnode_[node],deltaQ);
    LINALG::Matrix<3,1> deltatheta;
    LARGEROTATIONS::quaterniontoangle(deltaQ,deltatheta);

    //isotropic artificial stiffness
    LINALG::Matrix<3,3> artstiff;
    artstiff = LARGEROTATIONS::Tmatrix(deltatheta);

    //scale artificial damping with crotptc parameter for PTC method
    artstiff.Scale( params.get<double>("crotptc",0.0) );

    //each node gets a block diagonal damping term; the Lobatto integration weight is 0.5 for 2-noded elements
    //jacobi determinant is constant and equals 0.5*refelelength for 2-noded elements
    for(int k=0; k<3; k++)
      for (int l=0; l<3; l++)
        elemat1(node*6+3+k,node*6+3+l) += artstiff(k,l)*0.5*0.5*reflength_;

    //PTC for translational degrees of freedom; the Lobatto integration weight is 0.5 for 2-noded elements
    for(int k=0; k<3; k++)
      elemat1(node*6+k,node*6+k) += params.get<double>("ctransptc",0.0)*0.5*0.5*reflength_;
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
     //lrefe += gausspointsdamping.qwgt[gp]*jacobiGP_[gp];

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

  //get pre-calculated jacobi factors at GPs
  std::vector<double> jacobi(jacobiGPdampstoch_);

  //damping coefficients for translational and rotational degrees of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma);

  for (int gp=0; gp<gausspointsdamping.nquad; gp++)//loop through Gauss points
  {
    //compute triad at Gauss point
    LARGEROTATIONS::quaterniontotriad(Qnewdamping[gp],Lambdadamping);

    //rotation between last converged position and current position expressed as a quaternion
    LINALG::Matrix<4,1> deltaQ;
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

  //damping coefficients for translational and rotational degrees of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma);

  //get vector jacobi with Jacobi determinants at each integration point
  std::vector<double> jacobi(jacobiGPdampstoch_);

  //get Gauss points and weights for evaluation of damping matrix
  DRT::UTILS::GaussRule1D gaussrule = MyGaussRule(params,res_damp_stoch);
  DRT::UTILS::IntegrationPoints1D gausspoints(gaussrule);

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
                                              const std::vector<double>&       vel,  //!< element velocity vector   TODO unused?
                                              const std::vector<double>&       disp, //!<element disp vector
                                              Epetra_SerialDenseMatrix*        stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseVector*        force)//!< element internal force vector
{
  //damping coefficients for three translational and one rotational degree of freedom
  LINALG::Matrix<3,1> gamma(true);
  MyDampingConstants(params,gamma);

  //get vector jacobi with Jacobi determinants at each integration point
  std::vector<double> jacobi(jacobiGPdampstoch_);

  //get Gauss points and weights for evaluation of damping matrix
  DRT::UTILS::GaussRule1D gaussrule = MyGaussRule(params,res_damp_stoch);
  DRT::UTILS::IntegrationPoints1D gausspoints(gaussrule);

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
                                              const std::vector<double>&                              vel,  //!< element velocity vector   TODO unused?
                                              const std::vector<double>&                              disp, //!<element disp vector        TODO unused?
                                              Epetra_SerialDenseMatrix*                               stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseVector*                               force, //!< element internal force vector
                                              const DRT::UTILS::IntegrationPoints1D&                  gausspointsdamping,
                                              const std::vector<LINALG::Matrix<1,nnode> >&            Idamping,
                                              const std::vector<std::vector<LINALG::Matrix<3,3> > >&  Itildedamping,
                                              const std::vector<LINALG::Matrix<4,1> >&                Qconvdamping,    // TODO unused?
                                              const std::vector<LINALG::Matrix<4,1> >&                Qnewdamping)
{
  //auxiliary matrix
  LINALG::Matrix<3,3> auxmatrix;

  //determine type of numerical integration performed (lumped damping matrix via lobatto integration!)
  std::vector<double> jacobi(jacobiGPdampstoch_);

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
inline void DRT::ELEMENTS::Beam3r::CalcBrownian(Teuchos::ParameterList&                        params,
                                              const std::vector<double>&                       vel,  //!< element velocity vector
                                              const std::vector<double>&                       disp, //!< element displacement vector
                                              Epetra_SerialDenseMatrix*                        stiffmatrix,  //!< element stiffness matrix
                                              Epetra_SerialDenseVector*                        force,
                                              std::vector<LINALG::Matrix<1,nnode> >&           Idamping,
                                              std::vector<std::vector<LINALG::Matrix<3,3> > >& Itildedamping) //!< element internal force vector
{
  /* for integration of damping matrix always nnode Gauss points required; but in case of Lobatto integration
   * these are identical to the nnode nodes and then the basis functions are no longer the one also required
   * for the mass matrix, but rather their values at the integration points are given by a Kronecker-Delta function*/
  DRT::UTILS::GaussRule1D gaussrule_damp_stoch = MyGaussRule(params,res_damp_stoch);
  DRT::UTILS::IntegrationPoints1D gausspointsdamping(gaussrule_damp_stoch);

  std::vector<LINALG::Matrix<4,1> > Qconvdamping(QconvGPmass_);
  std::vector<LINALG::Matrix<4,1> > Qnewdamping(QnewGPmass_);

  if(gaussrule_damp_stoch == DRT::UTILS::intrule_line_lobatto2point)
  {
    //in case of Lobatto integration nodal triads are triads at Gauss points
    Qconvdamping = Qconvnode_;
    Qnewdamping  = Qnewnode_;

    Idamping.resize(gausspointsdamping.nquad);

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
  else if(gaussrule_damp_stoch != DRT::UTILS::intrule_line_2point)
    dserror("beam3r::CalcBrownian: this implementation assumes either intrule_line_2point or intrule_line_lobatto2point integration scheme");

  //now start with evaluation of force vectors and stiffness matrices

  //add stiffness and forces due to translational damping effects
  MyTranslationalDamping<nnode,ndim,dof>(params,vel,disp,stiffmatrix,force);

  //add stiffness and forces (i.e. moments) due to rotational damping effects
  MyRotationalDamping<nnode>(params,vel,disp,stiffmatrix,force,gausspointsdamping,Idamping,Itildedamping,Qconvdamping,Qnewdamping);

  //add stochastic forces and (if required) resulting stiffness
  MyStochasticForces<nnode,ndim,dof,randompergauss>(params,vel,disp,stiffmatrix,force);

  //add stochastic moments and resulting stiffness // TODO what about StochasticMoments???
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
                                             std::vector<double>&   disp) //!<element disp vector
{
  if(centerline_hermite_ or nnode>2)
    dserror("method 'Beam3r::NodeShift' was originally implemented for 2-noded Reissner beam element only. Check functionality for numnodes>2 and/or Hermite interpolation and extend if needed!");

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
        if( std::fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) + periodlength->at(dof) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) < std::fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) )
        {
          disp[numdof*i+dof] += periodlength->at(dof);

          /*the upper domain surface orthogonal to the z-direction may be subject to shear Dirichlet boundary condition; the lower surface
           *may be fixed by DBC. To avoid problmes when nodes exit the domain through the upper z-surface and reenter through the lower
           *z-surface, the shear has to be substracted from nodal coordinates in that case */
          if(shearflow && dof == 2 && curvenumber >=  0 && time>starttime && std::fabs(time-starttime)>dt/1e4)
            disp[numdof*i+dbcdispdir] += shearamplitude*DRT::Problem::Instance()->Curve(curvenumber).f(time);
        }

        if( std::fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - periodlength->at(dof) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) < std::fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) )
        {
          disp[numdof*i+dof] -= periodlength->at(dof);

          /*the upper domain surface orthogonal to the z-direction may be subject to shear Dirichlet boundary condition; the lower surface
           *may be fixed by DBC. To avoid problmes when nodes exit the domain through the lower z-surface and reenter through the upper
           *z-surface, the shear has to be added to nodal coordinates in that case */
          if(shearflow && dof == 2 && curvenumber >=  0 && time>starttime && std::fabs(time-starttime)>dt/1e4 )
            disp[numdof*i+dbcdispdir] -= shearamplitude*DRT::Problem::Instance()->Curve(curvenumber).f(time);
        }
      }
    }

return;
}//DRT::ELEMENTS::Beam3r::NodeShift

/*----------------------------------------------------------------------------------------------------------*
 | Get position vector at xi for given nodal displacements                                        popp 02/16|
 *----------------------------------------------------------------------------------------------------------*/
LINALG::Matrix<3,1> DRT::ELEMENTS::Beam3r::GetPos(double& xi, LINALG::Matrix<12,1>& disp_totlag) const
{
  LINALG::Matrix<3,1> r(true);
  LINALG::Matrix<4,1> N_i(true);

  if (!centerline_hermite_)
    dserror("ERROR: GetPos() method only for Hermite center lines");

  DRT::UTILS::shape_function_hermite_1D(N_i,xi,reflength_,line2);

  for (int n=0;n<4;n++)
  {
    for (int i=0;i<3;i++)
    {
      r(i)+=N_i(n)*disp_totlag(3*n+i);
    }
  }

  return (r);
}
