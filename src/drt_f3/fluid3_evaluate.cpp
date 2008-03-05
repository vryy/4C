/*!
\file fluid3_evaluate.cpp
\brief

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
#ifdef D_FLUID3
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif

#include "fluid3.H"
#include "fluid3_impl.H"
#include "fluid3_lin_impl.H"
#include "fluid3_genalpha_resVMM.H"
#include "fluid3_stationary.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_mat/newtonianfluid.H"

#include <blitz/array.h>
#include <Epetra_SerialDenseSolver.h>

using namespace DRT::UTILS;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;


/*
  Depending on the type of the algorithm (the implementation) and the
  element type (tet, hex etc.), the elements allocate common static
  arrays.

  That means that for example all hex8 fluid elements of the stationary
  implementation have a pointer f8 to the same 'implementation class'
  containing all the element arrays for eight noded elements, and all
  wedge15 fluid elements of the same problem have a pointer f15 to
  the 'implementation class' containing all the element arrays for the
  15 noded element.
  
  */

DRT::ELEMENTS::Fluid3Impl* DRT::ELEMENTS::Fluid3::Impl()
{
  switch (NumNode())
  {
  case 8:
  {
    static Fluid3Impl* f8;
    if (f8==NULL)
      f8 = new Fluid3Impl(8);
    return f8;
  }
  case 20:
  {
    static Fluid3Impl* f20;
    if (f20==NULL)
      f20 = new Fluid3Impl(20);
    return f20;
  }
  case 27:
  {
    static Fluid3Impl* f27;
    if (f27==NULL)
      f27 = new Fluid3Impl(27);
    return f27;
  }
  case 4:
  {
    static Fluid3Impl* f4;
    if (f4==NULL)
      f4 = new Fluid3Impl(4);
    return f4;
  }
  case 10:
  {
    static Fluid3Impl* f10;
    if (f10==NULL)
      f10 = new Fluid3Impl(10);
    return f10;
  }
  case 6:
  {
    static Fluid3Impl* f6;
    if (f6==NULL)
      f6 = new Fluid3Impl(6);
    return f6;
  }
  case 15:
  {
    static Fluid3Impl* f15;
    if (f15==NULL)
      f15 = new Fluid3Impl(15);
    return f15;
  }
  case 5:
  {
    static Fluid3Impl* f5;
    if (f5==NULL)
      f5 = new Fluid3Impl(5);
    return f5;
  }

  default:
    dserror("node number %d not supported", NumNode());
  }
  return NULL;
}

DRT::ELEMENTS::Fluid3GenalphaResVMM* DRT::ELEMENTS::Fluid3::GenalphaResVMM()
{
  switch (NumNode())
  {
  case 8:
  {
    static Fluid3GenalphaResVMM* f8;
    if (f8==NULL)
      f8 = new Fluid3GenalphaResVMM(8);
    return f8;
  }
  case 20:
  {
    static Fluid3GenalphaResVMM* f20;
    if (f20==NULL)
      f20 = new Fluid3GenalphaResVMM(20);
    return f20;
  }
  case 27:
  {
    static Fluid3GenalphaResVMM* f27;
    if (f27==NULL)
      f27 = new Fluid3GenalphaResVMM(27);
    return f27;
  }
  case 4:
  {
    static Fluid3GenalphaResVMM* f4;
    if (f4==NULL)
      f4 = new Fluid3GenalphaResVMM(4);
    return f4;
  }
  case 10:
  {
    static Fluid3GenalphaResVMM* f10;
    if (f10==NULL)
      f10 = new Fluid3GenalphaResVMM(10);
    return f10;
  }
  case 6:
  {
    static Fluid3GenalphaResVMM* f6;
    if (f6==NULL)
      f6 = new Fluid3GenalphaResVMM(6);
    return f6;
  }
  case 15:
  {
    static Fluid3GenalphaResVMM* f15;
    if (f15==NULL)
      f15 = new Fluid3GenalphaResVMM(15);
    return f15;
  }
  default:
    dserror("node number %d not supported", NumNode());
  }
  return NULL;
}


DRT::ELEMENTS::Fluid3Stationary* DRT::ELEMENTS::Fluid3::StationaryImpl()
{
  switch (NumNode())
  {
  case 8:
  {
    static Fluid3Stationary* f8;
    if (f8==NULL)
      f8 = new Fluid3Stationary(8);
    return f8;
  }
  case 20:
  {
    static Fluid3Stationary* f20;
    if (f20==NULL)
      f20 = new Fluid3Stationary(20);
    return f20;
  }
  case 27:
  {
    static Fluid3Stationary* f27;
    if (f27==NULL)
      f27 = new Fluid3Stationary(27);
    return f27;
  }
  case 4:
  {
    static Fluid3Stationary* f4;
    if (f4==NULL)
      f4 = new Fluid3Stationary(4);
    return f4;
  }
  case 10:
  {
    static Fluid3Stationary* f10;
    if (f10==NULL)
      f10 = new Fluid3Stationary(10);
    return f10;
  }
  case 6:
  {
    static Fluid3Stationary* f6;
    if (f6==NULL)
      f6 = new Fluid3Stationary(6);
    return f6;
  }
  case 15:
  {
    static Fluid3Stationary* f15;
    if (f15==NULL)
      f15 = new Fluid3Stationary(15);
    return f15;
  }
  case 5:
  {
    static Fluid3Stationary* f5;
    if (f5==NULL)
      f5 = new Fluid3Stationary(5);
    return f5;
  }

  default:
    dserror("node number %d not supported", NumNode());
  }
  return NULL;
}

DRT::ELEMENTS::Fluid3lin_Impl* DRT::ELEMENTS::Fluid3::F3linImpl()
{
  switch (NumNode())
  {
  case 8:
  {
    static Fluid3lin_Impl* f8;
    if (f8==NULL)
      f8 = new Fluid3lin_Impl(8);
    return f8;
  }
  case 20:
  {
    static Fluid3lin_Impl* f20;
    if (f20==NULL)
      f20 = new Fluid3lin_Impl(20);
    return f20;
  }
  case 27:
  {
    static Fluid3lin_Impl* f27;
    if (f27==NULL)
      f27 = new Fluid3lin_Impl(27);
    return f27;
  }
  case 4:
  {
    static Fluid3lin_Impl* f4;
    if (f4==NULL)
      f4 = new Fluid3lin_Impl(4);
    return f4;
  }
  case 10:
  {
    static Fluid3lin_Impl* f10;
    if (f10==NULL)
      f10 = new Fluid3lin_Impl(10);
    return f10;
  }
  case 6:
  {
    static Fluid3lin_Impl* f6;
    if (f6==NULL)
      f6 = new Fluid3lin_Impl(6);
    return f6;
  }
  case 15:
  {
    static Fluid3lin_Impl* f15;
    if (f15==NULL)
      f15 = new Fluid3lin_Impl(15);
    return f15;
  }
  case 5:
  {
    static Fluid3lin_Impl* f5;
    if (f5==NULL)
      f5 = new Fluid3lin_Impl(5);
    return f5;
  }

  default:
    dserror("node number %d not supported", NumNode());
  }
  return NULL;
}

/*---------------------------------------------------------------------*
|  converts a string into an Action for this element                   |
*----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3::ActionType DRT::ELEMENTS::Fluid3::convertStringToActionType(
              const string& action) const
{
  dsassert(action != "none", "No action supplied");

  DRT::ELEMENTS::Fluid3::ActionType act = Fluid3::none;
  if (action == "calc_fluid_systemmat_and_residual")
    act = Fluid3::calc_fluid_systemmat_and_residual;
  else if (action == "calc_linear_fluid")
    act = Fluid3::calc_linear_fluid;
  else if (action == "calc_fluid_genalpha_sysmat_and_residual")
    act = Fluid3::calc_fluid_genalpha_sysmat_and_residual;
  else if (action == "time update for subscales")
    act = Fluid3::calc_fluid_genalpha_update_for_subscales;
  else if (action == "time average for subscales and residual")
    act = Fluid3::calc_fluid_genalpha_average_for_subscales_and_residual;
  else if (action == "calc_fluid_stationary_systemmat_and_residual")
    act = Fluid3::calc_fluid_stationary_systemmat_and_residual;  
  else if (action == "calc_fluid_beltrami_error")
    act = Fluid3::calc_fluid_beltrami_error;
  else if (action == "calc_turbulence_statistics")
    act = Fluid3::calc_turbulence_statistics;
  else if (action == "calc_fluid_box_filter")
    act = Fluid3::calc_fluid_box_filter;
  else if (action == "calc_smagorinsky_const")
    act = Fluid3::calc_smagorinsky_const;
  else
    dserror("Unknown type of action for Fluid3");
  return act;
}

/*----------------------------------------------------------------------*
 // converts a string into an stabilisation action for this element
 //                                                          gammi 08/07
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3::StabilisationAction DRT::ELEMENTS::Fluid3::ConvertStringToStabAction(
  const string& action) const
{
  DRT::ELEMENTS::Fluid3::StabilisationAction act = stabaction_unspecified;

  map<string,StabilisationAction>::const_iterator iter=stabstrtoact_.find(action);

  if (iter != stabstrtoact_.end())
  {
    act = (*iter).second;
  }
  else
  {
    char errorout[200];
    sprintf(errorout,"looking for stab action (%s) not contained in map",action.c_str());
    dserror(errorout);
  }
  return act;
}


 /*----------------------------------------------------------------------*
 |  evaluate the element (public)                            g.bau 03/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Fluid3::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  // get the action required
  const string action = params.get<string>("action","none");
  const DRT::ELEMENTS::Fluid3::ActionType act = convertStringToActionType(action);

  // get the material
  RefCountPtr<MAT::Material> mat = Material();
  if (mat->MaterialType()!=m_fluid)
    dserror("newtonian fluid material expected but got type %d", mat->MaterialType());

  MATERIAL* actmat = static_cast<MAT::NewtonianFluid*>(mat.get())->MaterialData();

  switch(act)
  {
    //--------------------------------------------------
    //--------------------------------------------------
    // the standard one-step-theta implementation
    //--------------------------------------------------
    //--------------------------------------------------
      case calc_fluid_systemmat_and_residual:
      {
        // the number of nodes
        const int numnode = NumNode();

        //--------------------------------------------------
        // get all state vectors
        //--------------------------------------------------
        // 
        // need current velocity and history vector
        RefCountPtr<const Epetra_Vector> velnp = discretization.GetState("velnp");
        RefCountPtr<const Epetra_Vector> hist  = discretization.GetState("hist");
        if (velnp==null || hist==null)
          dserror("Cannot get state vectors 'velnp' and/or 'hist'");

        // extract local values from the global vectors
        vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);
        vector<double> myhist(lm.size());
        DRT::UTILS::ExtractMyValues(*hist,myhist,lm);

        RefCountPtr<const Epetra_Vector> dispnp;
        vector<double> mydispnp;
        RefCountPtr<const Epetra_Vector> gridv;
        vector<double> mygridv;

        if (is_ale_)
        {
          dispnp = discretization.GetState("dispnp");
          if (dispnp==null) dserror("Cannot get state vectors 'dispnp'");
          mydispnp.resize(lm.size());
          DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);

          gridv = discretization.GetState("gridv");
          if (gridv==null) dserror("Cannot get state vectors 'gridv'");
          mygridv.resize(lm.size());
          DRT::UTILS::ExtractMyValues(*gridv,mygridv,lm);
        }

        // create blitz objects for element arrays
        blitz::Array<double, 1> eprenp(numnode);
        blitz::Array<double, 2> evelnp(3,numnode,blitz::ColumnMajorArray<2>());
        blitz::Array<double, 2> evhist(3,numnode,blitz::ColumnMajorArray<2>());
        blitz::Array<double, 2> edispnp(3,numnode,blitz::ColumnMajorArray<2>());
        blitz::Array<double, 2> egridv(3,numnode,blitz::ColumnMajorArray<2>());

        // split velocity and pressure, insert into element arrays
        for (int i=0;i<numnode;++i)
        {
          evelnp(0,i) = myvelnp[0+(i*4)];
          evelnp(1,i) = myvelnp[1+(i*4)];
          evelnp(2,i) = myvelnp[2+(i*4)];

          eprenp(i) = myvelnp[3+(i*4)];

          // the history vector contains the information of time step t_n (mass rhs!)
          evhist(0,i) = myhist[0+(i*4)];
          evhist(1,i) = myhist[1+(i*4)];
          evhist(2,i) = myhist[2+(i*4)];
        }

        if (is_ale_)
        {
          // assign grid velocity and grid displacement to element arrays
          for (int i=0;i<numnode;++i)
          {
            edispnp(0,i) = mydispnp[0+(i*4)];
            edispnp(1,i) = mydispnp[1+(i*4)];
            edispnp(2,i) = mydispnp[2+(i*4)];

            egridv(0,i) = mygridv[0+(i*4)];
            egridv(1,i) = mygridv[1+(i*4)];
            egridv(2,i) = mygridv[2+(i*4)];
          }
        }

        // get fine-scale velocity
        RCP<const Epetra_Vector> fsvelnp;
        blitz::Array<double, 2> fsevelnp(3,numnode,blitz::ColumnMajorArray<2>());

        // get flag for (fine-scale) subgrid viscosity (1=yes, 0=no)
        const int fssgv = params.get<int>("fs subgrid viscosity",0);
        
        if (fssgv > 0)
        {
          fsvelnp = discretization.GetState("fsvelnp");
          if (fsvelnp==null) dserror("Cannot get state vector 'fsvelnp'");
          vector<double> myfsvelnp(lm.size());
          DRT::UTILS::ExtractMyValues(*fsvelnp,myfsvelnp,lm);

          // get fine-scale velocity and insert into element arrays
          for (int i=0;i<numnode;++i)
          {
            fsevelnp(0,i) = myfsvelnp[0+(i*4)];
            fsevelnp(1,i) = myfsvelnp[1+(i*4)];
            fsevelnp(2,i) = myfsvelnp[2+(i*4)];
          }
        }
        else
        {
          for (int i=0;i<numnode;++i)
          {
            fsevelnp(0,i) = 0.0;
            fsevelnp(1,i) = 0.0;
            fsevelnp(2,i) = 0.0;
          }
        }
        
        //--------------------------------------------------
        // get all control parameters for time integration
        // and stabilization
        //--------------------------------------------------
        //
        // get control parameter
        const double time = params.get<double>("total time",-1.0);

        const bool newton = params.get<bool>("include reactive terms for linearisation",false);

        // set parameters for stabilization
        ParameterList& stablist = params.sublist("STABILIZATION");

        // if not available, define map from string to action
        if(stabstrtoact_.empty())
        {
          stabstrtoact_["no_pspg"        ]=pstab_assume_inf_sup_stable;
          stabstrtoact_["yes_pspg"       ]=pstab_use_pspg;
          stabstrtoact_["no_supg"        ]=convective_stab_none;
          stabstrtoact_["yes_supg"       ]=convective_stab_supg;
          stabstrtoact_["no_vstab"       ]=viscous_stab_none;
          stabstrtoact_["vstab_gls"      ]=viscous_stab_gls;
          stabstrtoact_["vstab_gls_rhs"  ]=viscous_stab_gls_only_rhs;
          stabstrtoact_["vstab_usfem"    ]=viscous_stab_usfem;
          stabstrtoact_["vstab_usfem_rhs"]=viscous_stab_usfem_only_rhs;
          stabstrtoact_["no_cstab"       ]=continuity_stab_none;
          stabstrtoact_["cstab_qs"       ]=continuity_stab_yes;
          stabstrtoact_["no_cross"       ]=cross_stress_stab_none;
          stabstrtoact_["cross_complete" ]=cross_stress_stab;
          stabstrtoact_["cross_rhs"      ]=cross_stress_stab_only_rhs;
          stabstrtoact_["no_reynolds"    ]=reynolds_stress_stab_none;
          stabstrtoact_["reynolds_rhs"   ]=reynolds_stress_stab_only_rhs;
        }

        StabilisationAction pspg     = ConvertStringToStabAction(stablist.get<string>("PSPG"));
        StabilisationAction supg     = ConvertStringToStabAction(stablist.get<string>("SUPG"));
        StabilisationAction vstab    = ConvertStringToStabAction(stablist.get<string>("VSTAB"));
        StabilisationAction cstab    = ConvertStringToStabAction(stablist.get<string>("CSTAB"));
        StabilisationAction cross    = ConvertStringToStabAction(stablist.get<string>("CROSS-STRESS"));
        StabilisationAction reynolds = ConvertStringToStabAction(stablist.get<string>("REYNOLDS-STRESS"));

        // One-step-Theta: timefac = theta*dt
        // BDF2:           timefac = 2/3 * dt
        const double timefac = params.get<double>("thsl",-1.0);
        if (timefac < 0.0) dserror("No thsl supplied");

        // --------------------------------------------------
        // init fine scale Smagorinsky 
        // --------------------------------------------------
        
        // get Smagorinsky model parameter for fine-scale subgrid viscosity
        const double Cs_fs = params.get<double>("fs Smagorinsky parameter",0.0);

        // --------------------------------------------------
        // set parameters for classical turbulence models
        // --------------------------------------------------
        ParameterList& turbmodelparams    = params.sublist("TURBULENCE MODEL");

        // initialise the Smagorinsky constant Cs and the viscous length scale l_tau to zero
        double Cs            = 0.0;
        double Cs_delta_sq   = 0.0;
        double l_tau         = 0.0;
        double visceff       = 0.0;

        // the default action is no model        
        TurbModelAction turb_mod_action = no_model;
      
        // remember the layer of averaging for the dynamic Smagorinsky model
        int  nlayer=0;
        
        if (turbmodelparams.get<string>("TURBULENCE_APPROACH", "none") == "CLASSICAL_LES")
        {
          string& physical_turbulence_model = turbmodelparams.get<string>("PHYSICAL_MODEL");
          
          // --------------------------------------------------
          // standard constant coefficient Smagorinsky model
          if (physical_turbulence_model == "Smagorinsky")
          {
            // the classic Smagorinsky model only requires one constant parameter
            turb_mod_action = smagorinsky;
            Cs              = turbmodelparams.get<double>("C_SMAGORINSKY");
          }
          // --------------------------------------------------
          // Smagorinsky model with van Driest damping
          else if (physical_turbulence_model == "Smagorinsky_with_van_Driest_damping")
          {
            // that's only implemented for turbulent channel flow
            if (turbmodelparams.get<string>("CANONICAL_FLOW","no")
                !=
                "channel_flow_of_height_2")
            {
              dserror("van_Driest_damping only for channel_flow_of_height_2\n");
            }

            // for the Smagorinsky model with van Driest damping, we need
            // a viscous length to determine the y+ (heigth in wall units)
            turb_mod_action = smagorinsky_with_wall_damping;

            // get parameters of model
            Cs              = turbmodelparams.get<double>("C_SMAGORINSKY");
            l_tau           = turbmodelparams.get<double>("CHANNEL_L_TAU");

            // this will be the y-coordinate of a point in the element interior
            // we will determine the element layer in which he is contained to
            // be able to do the output of visceff etc.
            double center = 0;
            
            for(int inode=0;inode<numnode;inode++)
            {
              center+=Nodes()[inode]->X()[1];
            }
            center/=numnode;
        
            // node coordinates of plane to the element layer
            RefCountPtr<vector<double> > planecoords
              =
              turbmodelparams.get<RefCountPtr<vector<double> > >("planecoords_");
    
            bool found = false;
            for (nlayer=0;nlayer<(int)(*planecoords).size()-1;)
            {
                if(center<(*planecoords)[nlayer+1])
                {
                  found = true;
                  break;
                }
                nlayer++;
            }
            if (found ==false)
            {
              dserror("could not determine element layer");
            }
          }
          // --------------------------------------------------
          // Smagorinsky model with dynamic Computation of Cs
          else if (physical_turbulence_model == "Dynamic_Smagorinsky")
          {
            turb_mod_action = dynamic_smagorinsky;

            // for turbulent channel flow, use averaged quantities
            if (turbmodelparams.get<string>("CANONICAL_FLOW","no")
                ==
                "channel_flow_of_height_2")
            {
              RCP<vector<double> > averaged_LijMij
                =
                turbmodelparams.get<RCP<vector<double> > >("averaged_LijMij_");
              RCP<vector<double> > averaged_MijMij
                =
                turbmodelparams.get<RCP<vector<double> > >("averaged_MijMij_");
              
              //this will be the y-coordinate of a point in the element interior
              // here, the layer is determined in order to get the correct
              // averaged value from the vector of averaged (M/L)ijMij
              double center = 0;
              for(int inode=0;inode<numnode;inode++)
              {
                center+=Nodes()[inode]->X()[1];
              }
              center/=numnode;

              RCP<vector<double> > planecoords
                =
                turbmodelparams.get<RCP<vector<double> > >("planecoords_");
                
              bool found = false;
              for (nlayer=0;nlayer<(int)(*planecoords).size()-1;)
              {
                if(center<(*planecoords)[nlayer+1])
                {
                  found = true;
                  break;
                }
                nlayer++;
              }
              if (found ==false)
              {
                dserror("could not determine element layer");
              }

              // Cs_delta_sq is set by the averaged quantities              
              Cs_delta_sq = 0.5 * (*averaged_LijMij)[nlayer]/(*averaged_MijMij)[nlayer] ;
              
              // clipping to get algorithm stable
              if (Cs_delta_sq<0)
              {
                Cs_delta_sq=0;
              }
            }
            else
            {
              // when no averaging was done, we just keep the calculated (clipped) value
              Cs_delta_sq = Cs_delta_sq_;
            }
          }
          else
          {
            dserror("Up to now, only Smagorinsky (constant coefficient with and without wall function as well as dynamic) is available");
          }
        }
          
        //--------------------------------------------------
        // wrap epetra serial dense objects in blitz objects
        //--------------------------------------------------
        blitz::Array<double, 2> estif(elemat1.A(),
                                      blitz::shape(elemat1.M(),elemat1.N()),
                                      blitz::neverDeleteData,
                                      blitz::ColumnMajorArray<2>());
        blitz::Array<double, 1> eforce(elevec1.Values(),
                                       blitz::shape(elevec1.Length()),
                                       blitz::neverDeleteData);

        //--------------------------------------------------
        // calculate element coefficient matrix and rhs     
        //--------------------------------------------------
        Impl()->Sysmat(this,
                       evelnp,
                       fsevelnp,
                       eprenp,
                       evhist,
                       edispnp,
                       egridv,
                       estif,
                       eforce,
                       actmat,
                       time,
                       timefac,
                       newton,
                       fssgv,
                       pspg,
                       supg,
                       vstab,
                       cstab,
                       cross,
                       reynolds,
                       turb_mod_action,
                       Cs,
                       Cs_delta_sq,
                       visceff,
                       l_tau,
                       Cs_fs);

        //--------------------------------------------------
        // output values of Cs, visceff and Cs_delta_sq
        //--------------------------------------------------
        if (turbmodelparams.get<string>("TURBULENCE_APPROACH", "none") == "CLASSICAL_LES")
        {
          string& physical_turbulence_model = turbmodelparams.get<string>("PHYSICAL_MODEL");
          
          if (physical_turbulence_model == "Dynamic_Smagorinsky"
              ||
              physical_turbulence_model ==  "Smagorinsky_with_van_Driest_damping"
            )
          {
            if (turbmodelparams.get<string>("CANONICAL_FLOW","no")
                ==
                "channel_flow_of_height_2")
            {
              // Cs was changed in Sysmat (Cs->sqrt(Cs/hk)) to compare it with the standard
              // Smagorinsky Cs
              
              if(this->Owner() == discretization.Comm().MyPID())
              {
                (*(turbmodelparams.get<RCP<vector<double> > >("local_Cs_sum")))         [nlayer]+=Cs;
                (*(turbmodelparams.get<RCP<vector<double> > >("local_Cs_delta_sq_sum")))[nlayer]+=Cs_delta_sq;
                (*(turbmodelparams.get<RCP<vector<double> > >("local_visceff_sum")))    [nlayer]+=visceff;
              }
            }
          }
        }

        
        // This is a very poor way to transport the density to the
        // outside world. Is there a better one?
        params.set("density", actmat->m.fluid->density);

      }
      break;
      case calc_linear_fluid:
      {
        // need current velocity and history vector
        RefCountPtr<const Epetra_Vector> velnp = discretization.GetState("velnp");
        RefCountPtr<const Epetra_Vector> hist  = discretization.GetState("hist");
        if (velnp==null || hist==null)
          dserror("Cannot get state vectors 'velnp' and/or 'hist'");

        // extract local values from the global vectors
        vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);
        vector<double> myhist(lm.size());
        DRT::UTILS::ExtractMyValues(*hist,myhist,lm);

        RefCountPtr<const Epetra_Vector> dispnp;
        vector<double> mydispnp;

        // create blitz objects for element arrays
        const int numnode = NumNode();
        blitz::Array<double, 1> eprenp(numnode);
        blitz::Array<double, 2> evelnp(3,numnode,blitz::ColumnMajorArray<2>());
        blitz::Array<double, 2> evhist(3,numnode,blitz::ColumnMajorArray<2>());

        // split velocity and pressure, insert into element arrays
        for (int i=0;i<numnode;++i)
        {
          evelnp(0,i) = myvelnp[0+(i*4)];
          evelnp(1,i) = myvelnp[1+(i*4)];
          evelnp(2,i) = myvelnp[2+(i*4)];

          eprenp(i) = myvelnp[3+(i*4)];

          // the history vector contains the information of time step t_n (mass rhs!)
          evhist(0,i) = myhist[0+(i*4)];
          evhist(1,i) = myhist[1+(i*4)];
          evhist(2,i) = myhist[2+(i*4)];
        }

        // get control parameter
        const double time = params.get<double>("total time",-1.0);

        // One-step-Theta: timefac = theta*dt
        // BDF2:           timefac = 2/3 * dt
        const double timefac = params.get<double>("thsl",-1.0);
        if (timefac < 0.0) dserror("No thsl supplied");

        // wrap epetra serial dense objects in blitz objects
        blitz::Array<double, 2> estif(elemat1.A(),
                                      blitz::shape(elemat1.M(),elemat1.N()),
                                      blitz::neverDeleteData,
                                      blitz::ColumnMajorArray<2>());
        blitz::Array<double, 1> eforce(elevec1.Values(),
                                       blitz::shape(elevec1.Length()),
                                       blitz::neverDeleteData);

        // calculate element coefficient matrix and rhs     
        F3linImpl()->Sysmat(this,
			    evelnp,
			    eprenp,
			    evhist,
			    estif,
			    eforce,
			    actmat,
			    time,
			    timefac);
      }
      break;
      case calc_fluid_beltrami_error:
      {
        // add error only for elements which are not ghosted
        if(this->Owner() == discretization.Comm().MyPID())
        {

          // need current velocity and history vector
          RefCountPtr<const Epetra_Vector> vel_pre_np = discretization.GetState("u and p at time n+1 (converged)");
          if (vel_pre_np==null) dserror("Cannot get state vectors 'velnp'");

          // extract local values from the global vectors
          vector<double> my_vel_pre_np(lm.size());
          DRT::UTILS::ExtractMyValues(*vel_pre_np,my_vel_pre_np,lm);

          // split "my_vel_pre_np" into velocity part "myvelnp" and pressure part "myprenp"
          const int numnode = NumNode();
          vector<double> myprenp(numnode);
          vector<double> myvelnp(3*numnode);

          for (int i=0;i<numnode;++i)
          {
            myvelnp[0+(i*3)]=my_vel_pre_np[0+(i*4)];
            myvelnp[1+(i*3)]=my_vel_pre_np[1+(i*4)];
            myvelnp[2+(i*3)]=my_vel_pre_np[2+(i*4)];

            myprenp[i]=my_vel_pre_np[3+(i*4)];
          }

          // integrate beltrami error
          f3_int_beltrami_err(myvelnp,myprenp,actmat,params);
        }
      }
      break;
      case calc_turbulence_statistics:
      {
        // do nothing if you do not own this element
        if(this->Owner() == discretization.Comm().MyPID())
        {
          // --------------------------------------------------
          // extract velocities and pressure from the global distributed vectors

          // velocity and pressure values (n+1)
          RefCountPtr<const Epetra_Vector> velnp
            = discretization.GetState("u and p (n+1,converged)");


          if (velnp==null)
          {
            dserror("Cannot get state vectors 'velnp'");
          }

          // extract local values from the global vectors
          vector<double> mysol  (lm.size());
          DRT::UTILS::ExtractMyValues(*velnp,mysol,lm);

          // integrate mean values
          f3_calc_means(mysol,params);

        }
      }
      break;
      case calc_fluid_box_filter:
      {
        // --------------------------------------------------
        // extract velocities from the global distributed vectors

        // velocity and pressure values (most recent
        // intermediate solution, i.e. n+alphaF for genalpha
        // and n+1 for one-step-theta)
        RefCountPtr<const Epetra_Vector> vel =
          discretization.GetState("u and p (trial)");

        if (vel==null)
        {
          dserror("Cannot get state vectors 'vel'");
        }

        // extract local values from the global vectors
        vector<double> myvel(lm.size());
        DRT::UTILS::ExtractMyValues(*vel,myvel,lm);

        // create blitz objects for element arrays
        const int numnode = NumNode();
        blitz::Array<double, 2> evel(3,numnode,blitz::ColumnMajorArray<2>());

        // split velocity and throw away  pressure, insert into element array
        for (int i=0;i<numnode;++i)
        {
          evel(0,i) = myvel[0+(i*4)];
          evel(1,i) = myvel[1+(i*4)];
          evel(2,i) = myvel[2+(i*4)];
        }

        // initialise the contribution of this element to the patch volume to zero
        double volume_contribution = 0;

        // wrap Epetra Objects in Blitz array
        blitz::Array<double, 1> vel_hat(elevec1.Values(),
                                        blitz::shape(elevec1.Length()),
                                        blitz::neverDeleteData);

        blitz::Array<double, 2> reystress_hat(elemat1.A(),
                                              blitz::shape(elemat1.M(),elemat1.N()),
                                              blitz::neverDeleteData,
                                              blitz::ColumnMajorArray<2>());

        blitz::Array<double, 2> modeled_stress_grid_scale_hat(elemat2.A(),
                                                              blitz::shape(elemat2.M(),elemat2.N()),
                                                              blitz::neverDeleteData,
                                                              blitz::ColumnMajorArray<2>());

        // integrate the convolution with the box filter function for this element
        // the results are assembled onto the *_hat arrays
        this->f3_apply_box_filter(evel,
                                  vel_hat,
                                  reystress_hat,
                                  modeled_stress_grid_scale_hat,
                                  volume_contribution);

        // hand down the volume contribution to the time integration algorithm
        params.set<double>("volume_contribution",volume_contribution);

      }
      break;
      case calc_smagorinsky_const:
      {
        RefCountPtr<Epetra_MultiVector> filtered_vel                        =
          params.get<RefCountPtr<Epetra_MultiVector> >("col_filtered_vel");
        RefCountPtr<Epetra_MultiVector> col_filtered_reynoldsstress         =
          params.get<RefCountPtr<Epetra_MultiVector> >("col_filtered_reynoldsstress");
        RefCountPtr<Epetra_MultiVector> col_filtered_modeled_subgrid_stress =
          params.get<RefCountPtr<Epetra_MultiVector> >("col_filtered_modeled_subgrid_stress");

        const int numnode = this->NumNode();
        blitz::Array<double, 2> evelaf_hat                           (3  ,numnode,blitz::ColumnMajorArray<2>());
        blitz::Array<double, 3> ereynoldsstress_hat                  (3,3,numnode,blitz::ColumnMajorArray<3>());
        blitz::Array<double, 3> efiltered_modeled_subgrid_stress_hat (3,3,numnode,blitz::ColumnMajorArray<3>());

        
        for (int nn=0;nn<numnode;++nn)
        {
          DRT::Node* actnode = this->Nodes()[nn];
          
          for (int dimi=0;dimi<3;++dimi)
          {
            evelaf_hat(dimi,nn) = (*((*filtered_vel)(dimi)))[actnode->LID()];
            
            for (int dimj=0;dimj<3;++dimj)
            {
              ereynoldsstress_hat                  (dimi,dimj,nn)
                =
                (*((*col_filtered_reynoldsstress        )(3*dimi+dimj)))[actnode->LID()];
              efiltered_modeled_subgrid_stress_hat (dimi,dimj,nn)
                =
                (*((*col_filtered_modeled_subgrid_stress)(3*dimi+dimj)))[actnode->LID()];
            }
          }
        }

        double LijMij   = 0;
        double MijMij   = 0;
        double center   = 0;

        this->f3_calc_smag_const_LijMij_and_MijMij(evelaf_hat,
                                                   ereynoldsstress_hat,
                                                   efiltered_modeled_subgrid_stress_hat,
                                                   LijMij,
                                                   MijMij,
                                                   center);

        // set Cs_delta_sq without averaging (only clipping)
        Cs_delta_sq_ = 0.5 * LijMij / MijMij;
        if (Cs_delta_sq_<0)
        {
          Cs_delta_sq_= 0;
        }

        params.set<double>("LijMij",LijMij);
        params.set<double>("MijMij",MijMij);
        params.set<double>("center",center);
      }
      break;
      case calc_fluid_genalpha_sysmat_and_residual:
      {
#ifdef PERF
        // --------------------------------------------------
        // extended time measurement
        RefCountPtr<Time> timeelegetvelnp   = params.get<RefCountPtr<Time>  >("get global vectors an set node values");
        RefCountPtr<Time> timeeleinitsmag   = params.get<RefCountPtr<Time>  >("initialise Smagorinsky model");
        RefCountPtr<Time> timeeleinitstab   = params.get<RefCountPtr<Time>  >("initialise stabilization flags");
        RefCountPtr<Time> timeelesysmat     = params.get<RefCountPtr<Time>  >("complete call to sysmat");
#endif     


#ifdef PERF
        RefCountPtr<TimeMonitor> timeelegetvelnp_ref = rcp(new TimeMonitor(*timeelegetvelnp)); 
#endif     
        
        // --------------------------------------------------
        // extract velocities, pressure and accelerations from the
        // global distributed vectors

        // velocity and pressure values (current iterate, n+1)
        RefCountPtr<const Epetra_Vector> velnp = discretization.GetState("u and p (n+1      ,trial)");

        // velocities    (intermediate time step, n+alpha_F)
        RefCountPtr<const Epetra_Vector> velaf = discretization.GetState("u and p (n+alpha_F,trial)");

        // accelerations (intermediate time step, n+alpha_M)
        RefCountPtr<const Epetra_Vector> accam = discretization.GetState("acc     (n+alpha_M,trial)");

        if (velnp==null || velaf==null || accam==null)
        {
          dserror("Cannot get state vectors 'velnp', 'velaf'  and/or 'accam'");
        }

        // extract local values from the global vectors
        vector<double> myvelnp(lm.size());
        DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);

        vector<double> myvelaf(lm.size());
        DRT::UTILS::ExtractMyValues(*velaf,myvelaf,lm);

        vector<double> myaccam(lm.size());
        DRT::UTILS::ExtractMyValues(*accam,myaccam,lm);

        // create blitz matrix objects
        const int numnode = NumNode();
        blitz::Array<double, 1> eprenp    (  numnode);
        blitz::Array<double, 2> evelnp    (3,numnode,blitz::ColumnMajorArray<2>());
        blitz::Array<double, 2> evelaf    (3,numnode,blitz::ColumnMajorArray<2>());
        blitz::Array<double, 2> eaccam    (3,numnode,blitz::ColumnMajorArray<2>());
        blitz::Array<double, 2> edispnp   (3,numnode,blitz::ColumnMajorArray<2>());
        blitz::Array<double, 2> egridvelaf(3,numnode,blitz::ColumnMajorArray<2>());


        // split "my_velnp" into velocity part "myvelnp" and pressure part "myprenp"
        // Additionally only the 'velocity' components of my_velaf
        // and my_accam are important!
        for (int i=0;i<numnode;++i)
        {
          eprenp(i)   = myvelnp[3+(i*4)];

          evelnp(0,i) = myvelnp[0+(i*4)];
          evelnp(1,i) = myvelnp[1+(i*4)];
          evelnp(2,i) = myvelnp[2+(i*4)];

          evelaf(0,i) = myvelaf[0+(i*4)];
          evelaf(1,i) = myvelaf[1+(i*4)];
          evelaf(2,i) = myvelaf[2+(i*4)];

          eaccam(0,i) = myaccam[0+(i*4)];
          eaccam(1,i) = myaccam[1+(i*4)];
          eaccam(2,i) = myaccam[2+(i*4)];
        }


        
        if(is_ale_)
        {
          // get most recent displacements
          RefCountPtr<const Epetra_Vector> dispnp
            =
            discretization.GetState("dispnp");

          // get intermediate grid velocities
          RefCountPtr<const Epetra_Vector> gridvelaf
            =
            discretization.GetState("gridvelaf");

          if (dispnp==null || gridvelaf==null)
          {
            dserror("Cannot get state vectors 'dispnp' and/or 'gridvelaf'");
          }

          vector<double> mydispnp(lm.size());
          DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);

          vector<double> mygridvelaf(lm.size());
          DRT::UTILS::ExtractMyValues(*gridvelaf,mygridvelaf,lm);

          // extract velocity part from "mygridvelaf" and get
          // set element displacements
          for (int i=0;i<numnode;++i)
          {
            egridvelaf(0,i) = mygridvelaf[0+(i*4)];
            egridvelaf(1,i) = mygridvelaf[1+(i*4)];
            egridvelaf(2,i) = mygridvelaf[2+(i*4)];
            
            edispnp(0,i)    = mydispnp   [0+(i*4)];
            edispnp(1,i)    = mydispnp   [1+(i*4)];
            edispnp(2,i)    = mydispnp   [2+(i*4)];
          }
        }
        
#ifdef PERF
        timeelegetvelnp_ref = null;
#endif     
        
        // --------------------------------------------------
        // set parameters for time integration
        ParameterList& timelist = params.sublist("time integration parameters");

        const double alphaM = timelist.get<double>("alpha_M");
        const double alphaF = timelist.get<double>("alpha_F");
        const double gamma  = timelist.get<double>("gamma");
        const double dt     = timelist.get<double>("dt");
        const double time   = timelist.get<double>("time");

        // --------------------------------------------------
        // set parameters for nonlinear treatment

        const bool newton = params.get<bool>("include reactive terms for linearisation");

        // get flag for (fine-scale) subgrid viscosity (1=yes, 0=no)
        const int fssgv = params.get<int>("fs subgrid viscosity",0);

        // get fine-scale velocity
        RCP<const Epetra_Vector> fsvelaf;
        blitz::Array<double, 2> fsevelaf(3,numnode,blitz::ColumnMajorArray<2>());
        if (fssgv > 0)
        {
          fsvelaf = discretization.GetState("fsvelaf");
          if (fsvelaf==null) dserror("Cannot get state vector 'fsvelaf'");
          vector<double> myfsvelaf(lm.size());
          DRT::UTILS::ExtractMyValues(*fsvelaf,myfsvelaf,lm);

          // get fine-scale velocity and insert into element arrays
          for (int i=0;i<numnode;++i)
          {
            fsevelaf(0,i) = myfsvelaf[0+(i*4)];
            fsevelaf(1,i) = myfsvelaf[1+(i*4)];
            fsevelaf(2,i) = myfsvelaf[2+(i*4)];
          }
        }
        else
        {
          for (int i=0;i<numnode;++i)
          {
            fsevelaf(0,i) = 0.0;
            fsevelaf(1,i) = 0.0;
            fsevelaf(2,i) = 0.0;
          }
        }

#ifdef PERF
        RefCountPtr<TimeMonitor> timeeleinitstab_ref = rcp(new TimeMonitor(*timeeleinitstab)); 
#endif     
        
        // --------------------------------------------------
        // set parameters for stabilisation
        ParameterList& stablist = params.sublist("STABILIZATION");

        // if not available, define map from string to action
        if(stabstrtoact_.empty())
        {
          stabstrtoact_["quasistatic"    ]=subscales_quasistatic;
          stabstrtoact_["time_dependent" ]=subscales_time_dependent;
          stabstrtoact_["no_transient"   ]=inertia_stab_drop;
          stabstrtoact_["yes_transient"  ]=inertia_stab_keep;
          stabstrtoact_["no_pspg"        ]=pstab_assume_inf_sup_stable;
          stabstrtoact_["yes_pspg"       ]=pstab_use_pspg;
          stabstrtoact_["no_supg"        ]=convective_stab_none;
          stabstrtoact_["yes_supg"       ]=convective_stab_supg;
          stabstrtoact_["no_vstab"       ]=viscous_stab_none;
          stabstrtoact_["vstab_gls"      ]=viscous_stab_gls;
          stabstrtoact_["vstab_gls_rhs"  ]=viscous_stab_gls_only_rhs;
          stabstrtoact_["vstab_usfem"    ]=viscous_stab_usfem;
          stabstrtoact_["vstab_usfem_rhs"]=viscous_stab_usfem_only_rhs;
          stabstrtoact_["no_ctab"        ]=continuity_stab_none;
          stabstrtoact_["cstab_qs"       ]=continuity_stab_yes;
          stabstrtoact_["cstab_td"       ]=continuity_stab_td;
          stabstrtoact_["no_cross"       ]=cross_stress_stab_none;
          stabstrtoact_["cross_complete" ]=cross_stress_stab;
          stabstrtoact_["cross_rhs"      ]=cross_stress_stab_only_rhs;
          stabstrtoact_["no_reynolds"    ]=reynolds_stress_stab_none;
          stabstrtoact_["reynolds_rhs"   ]=reynolds_stress_stab_only_rhs;
        }

        StabilisationAction tds      = ConvertStringToStabAction(stablist.get<string>("TDS"));
        StabilisationAction inertia  = ConvertStringToStabAction(stablist.get<string>("TRANSIENT"));
        StabilisationAction pspg     = ConvertStringToStabAction(stablist.get<string>("PSPG"));
        StabilisationAction supg     = ConvertStringToStabAction(stablist.get<string>("SUPG"));
        StabilisationAction vstab    = ConvertStringToStabAction(stablist.get<string>("VSTAB"));
        StabilisationAction cstab    = ConvertStringToStabAction(stablist.get<string>("CSTAB"));
        StabilisationAction cross    = ConvertStringToStabAction(stablist.get<string>("CROSS-STRESS"));
        StabilisationAction reynolds = ConvertStringToStabAction(stablist.get<string>("REYNOLDS-STRESS"));

#ifdef PERF
        timeeleinitstab_ref = null; 
#endif

#ifdef PERF
        RefCountPtr<TimeMonitor> timeeleinitsmag_ref = rcp(new TimeMonitor(*timeeleinitsmag)); 
#endif     
                
        // --------------------------------------------------
        // set parameters for turbulence model
        ParameterList& turbmodelparams    = params.sublist("TURBULENCE MODEL");

        // get Smagorinsky model parameter for fine-scale subgrid viscosity
        const double Cs_fs = turbmodelparams.get<double>("C_SMAGORINSKY",0.0);

        // initialise the Smagorinsky constant Cs and the viscous length scale l_tau to zero
        double Cs            = 0.0;
        double Cs_delta_sq   = 0.0;
        double l_tau         = 0.0;
        double visceff       = 0.0;
        
        // the default action is no model        
        TurbModelAction turb_mod_action = no_model;


        int  nlayer=0;
        if (turbmodelparams.get<string>("TURBULENCE_APPROACH", "none") == "CLASSICAL_LES")
        {
          string& physical_turbulence_model = turbmodelparams.get<string>("PHYSICAL_MODEL");
 
          if (physical_turbulence_model == "Smagorinsky")
          {
            // the classic Smagorinsky model only requires one constant parameter
            turb_mod_action = smagorinsky;
            Cs              = turbmodelparams.get<double>("C_SMAGORINSKY");
          }
          else if (physical_turbulence_model == "Smagorinsky_with_van_Driest_damping")
          {
            RefCountPtr<vector<double> > planecoords      = turbmodelparams.get<RefCountPtr<vector<double> > >("planecoords_");

            // for the Smagorinsky model with van Driest damping, we need a viscous length to determine
            // the y+ (heigth in wall units)
            turb_mod_action = smagorinsky_with_wall_damping;
            Cs              = turbmodelparams.get<double>("C_SMAGORINSKY");
            l_tau           = turbmodelparams.get<double>("CHANNEL_L_TAU");

            const int iel = NumNode();
            
            //this will be the y-coordinate of a point in the element interior
            double center = 0;
            for(int inode=0;inode<iel;inode++)
            {
              center+=Nodes()[inode]->X()[1];
            }
            center/=iel;
            
            bool found = false;
            for (nlayer=0;nlayer<(int)(*planecoords).size()-1;)
            {
                if(center<(*planecoords)[nlayer+1])
                {
                  found = true;
                  break;
                }
                nlayer++;
            }
            if (found ==false)
            {
              dserror("could not determine element layer");
            }
          }
          else if (physical_turbulence_model == "Dynamic_Smagorinsky")
          {
            turb_mod_action = dynamic_smagorinsky;
            
            if (turbmodelparams.get<string>("CANONICAL_FLOW","no")
                ==
                "channel_flow_of_height_2")
            {
              RCP<vector<double> > averaged_LijMij
                =
                turbmodelparams.get<RCP<vector<double> > >("averaged_LijMij_");
              RCP<vector<double> > averaged_MijMij
                =
                turbmodelparams.get<RCP<vector<double> > >("averaged_MijMij_");
              
              RCP<vector<double> > planecoords
                =
                turbmodelparams.get<RCP<vector<double> > >("planecoords_");
              
              const int iel = NumNode();
            
              //this will be the y-coordinate of a point in the element interior
              double center = 0;
              for(int inode=0;inode<iel;inode++)
              {
                center+=Nodes()[inode]->X()[1];
              }
              center/=iel;
              
              bool found = false;
              for (nlayer=0;nlayer<(int)(*planecoords).size()-1;)
              {
                if(center<(*planecoords)[nlayer+1])
                {
                  found = true;
                  break;
                }
                nlayer++;
              }
              if (found ==false)
              {
                dserror("could not determine element layer");
              }
              
              Cs_delta_sq = 0.5 * (*averaged_LijMij)[nlayer]/(*averaged_MijMij)[nlayer] ;
              
              // clipping to get algorithm stable
              if (Cs_delta_sq<0)
              {
                Cs_delta_sq=0;
              }
            }
            else
            {
              Cs_delta_sq = Cs_delta_sq_;
            }
          }
          else
          {
            dserror("Up to now, only Smagorinsky (constant coefficient with and without wall function as well as dynamic) is available");
          }
        }


#ifdef PERF
        timeeleinitsmag_ref = null; 
#endif     
                
        
        // --------------------------------------------------
        // specify what to compute
        const bool compute_elemat = params.get<bool>("compute element matrix");

        
#ifdef PERF
        // --------------------------------------------------
        // extended time measurement
        RefCountPtr<Time> timeelederxy2    = params.get<RefCountPtr<Time>  >("time used for second derivatives");
        RefCountPtr<Time> timeelederxy     = params.get<RefCountPtr<Time>  >("time used for secondfirst");
        RefCountPtr<Time> timeeletau       = params.get<RefCountPtr<Time>  >("time used for computation of tau");
        RefCountPtr<Time> timeelegalerkin  = params.get<RefCountPtr<Time>  >("time used for galerkin loops");
        RefCountPtr<Time> timeelepspg      = params.get<RefCountPtr<Time>  >("time used for pspg loop");
        RefCountPtr<Time> timeelesupg      = params.get<RefCountPtr<Time>  >("time used for supg loop");
        RefCountPtr<Time> timeelecstab     = params.get<RefCountPtr<Time>  >("time used for cstab loop");
        RefCountPtr<Time> timeelevstab     = params.get<RefCountPtr<Time>  >("time used for vstab loop");
        RefCountPtr<Time> timeelecrossrey  = params.get<RefCountPtr<Time>  >("time used for cross and reynolds loop");
        RefCountPtr<Time> timeeleintertogp = params.get<RefCountPtr<Time>  >("time used to interpolate to gauss points");
        RefCountPtr<Time> timeeleseteledata= params.get<RefCountPtr<Time>  >("set basic element data");
        RefCountPtr<Time> timeeletdextras  = params.get<RefCountPtr<Time>  >("time dependent fine scale specials");
#endif

   
#ifdef PERF
        RefCountPtr<TimeMonitor> timeelesysmat_ref = rcp(new TimeMonitor(*timeelesysmat)); 
#endif          
        // --------------------------------------------------
        // calculate element coefficient matrix
        GenalphaResVMM()->Sysmat(this,
                                 elemat1,
                                 elevec1,
                                 edispnp,
                                 egridvelaf,
                                 evelnp,
                                 eprenp,
                                 eaccam,
                                 evelaf,
                                 fsevelaf,
                                 actmat,
                                 alphaM,
                                 alphaF,
                                 gamma,
                                 dt,
                                 time,
                                 newton,
                                 fssgv,
                                 Cs_fs,
                                 tds,
                                 inertia,
                                 pspg,
                                 supg,
                                 vstab,
                                 cstab,
                                 cross,
                                 reynolds,
                                 turb_mod_action,
                                 Cs,
                                 Cs_delta_sq,
                                 visceff,
                                 l_tau,
#ifdef PERF
                                 timeelederxy2   ,
                                 timeelederxy    ,
                                 timeeletau      ,
                                 timeelegalerkin ,
                                 timeelepspg     ,
                                 timeelesupg     ,
                                 timeelecstab    ,
                                 timeelevstab    ,
                                 timeelecrossrey ,  
                                 timeeleintertogp,
                                 timeeleseteledata,
                                 timeeletdextras  ,

#endif                                 
                                 compute_elemat
          );
#ifdef PERF
        timeelesysmat_ref = null;
#endif          

        
#ifdef PERF
        timeeleinitsmag_ref = rcp(new TimeMonitor(*timeeleinitsmag)); 
#endif
        
        if (turbmodelparams.get<string>("TURBULENCE_APPROACH", "none") == "CLASSICAL_LES")
        {
          string& physical_turbulence_model = turbmodelparams.get<string>("PHYSICAL_MODEL");
          
          if (physical_turbulence_model == "Dynamic_Smagorinsky"
              ||
              physical_turbulence_model ==  "Smagorinsky_with_van_Driest_damping"
            )
          {
            // Cs was changed in Sysmat (Cs->sqrt(Cs/hk)) to compare it with the standard
            // Smagorinsky Cs

            if(this->Owner() == discretization.Comm().MyPID())
            {
              (*(turbmodelparams.get<RCP<vector<double> > >("local_Cs_sum")))         [nlayer]+=Cs;
              (*(turbmodelparams.get<RCP<vector<double> > >("local_Cs_delta_sq_sum")))[nlayer]+=Cs_delta_sq;
              (*(turbmodelparams.get<RCP<vector<double> > >("local_visceff_sum")))    [nlayer]+=visceff;
            }
          }
        }

#ifdef PERF
        timeeleinitsmag_ref = null; 
#endif        
      }
      break;
      case calc_fluid_genalpha_update_for_subscales:
      {
        // most recent subscale pressure becomes the old subscale pressure
        // for the next timestep
        //
        //  ~n   ~n+1
        //  p <- p
        //
        sub_pre_old_ = sub_pre_;

        // the old subscale acceleration for the next timestep is calculated
        // on the fly, not stored on the element
        /*
                     ~n+1   ~n
             ~ n     u    - u     ~ n   / 1.0-gamma \
            acc  <-  --------- - acc * |  ---------  |
                     gamma*dt           \   gamma   /
        */
        
        const double dt     = params.get<double>("dt");
        const double gamma  = params.get<double>("gamma");

        sub_acc_old_ = (sub_vel_-sub_vel_old_)/(gamma*dt)
                       -
                       sub_acc_old_*(1.0-gamma)/gamma;

        // most recent subscale velocity becomes the old subscale velocity
        // for the next timestep
        //
        //  ~n   ~n+1
        //  u <- u
        //
        sub_vel_old_=sub_vel_;
      }
      break;
      case calc_fluid_genalpha_average_for_subscales_and_residual:
      {
        if(this->Owner() == discretization.Comm().MyPID())
        {
          // --------------------------------------------------
          // extract velocities, pressure and accelerations from the
          // global distributed vectors

          // velocity and pressure values (current iterate, n+1)
          RefCountPtr<const Epetra_Vector> velnp = discretization.GetState("u and p (n+1      ,trial)");

          // velocities    (intermediate time step, n+alpha_F)
          RefCountPtr<const Epetra_Vector> velaf = discretization.GetState("u and p (n+alpha_F,trial)");

          // accelerations (intermediate time step, n+alpha_M)
          RefCountPtr<const Epetra_Vector> accam = discretization.GetState("acc     (n+alpha_M,trial)");
          
          if (velnp==null || velaf==null || accam==null)
          {
            dserror("Cannot get state vectors 'velnp', 'velaf'  and/or 'accam'");
          }
          
          // extract local values from the global vectors
          vector<double> myvelnp(lm.size());
          DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);
          
          vector<double> myvelaf(lm.size());
          DRT::UTILS::ExtractMyValues(*velaf,myvelaf,lm);
          
          vector<double> myaccam(lm.size());
          DRT::UTILS::ExtractMyValues(*accam,myaccam,lm);
          
          // create blitz matrix objects
          
          const int numnode = NumNode();
          blitz::Array<double, 1> eprenp (  numnode);
          blitz::Array<double, 2> evelnp (3,numnode,blitz::ColumnMajorArray<2>());
          blitz::Array<double, 2> evelaf (3,numnode,blitz::ColumnMajorArray<2>());
          blitz::Array<double, 2> eaccam (3,numnode,blitz::ColumnMajorArray<2>());
          
          
          // split "my_velnp" into velocity part "myvelnp" and pressure part "myprenp"
          // Additionally only the 'velocity' components of my_velaf
          // and my_accam are important!
          for (int i=0;i<numnode;++i)
          {
            eprenp(i)   = myvelnp[3+(i*4)];
            
            evelnp(0,i) = myvelnp[0+(i*4)];
            evelnp(1,i) = myvelnp[1+(i*4)];
            evelnp(2,i) = myvelnp[2+(i*4)];
            
            evelaf(0,i) = myvelaf[0+(i*4)];
            evelaf(1,i) = myvelaf[1+(i*4)];
            evelaf(2,i) = myvelaf[2+(i*4)];
            
            eaccam(0,i) = myaccam[0+(i*4)];
            eaccam(1,i) = myaccam[1+(i*4)];
            eaccam(2,i) = myaccam[2+(i*4)];
          }
          
          // --------------------------------------------------
          // set parameters for time integration
          ParameterList& timelist = params.sublist("time integration parameters");
          
          const double alphaM = timelist.get<double>("alpha_M");
          const double alphaF = timelist.get<double>("alpha_F");
          const double gamma  = timelist.get<double>("gamma");
          const double dt     = timelist.get<double>("dt");
          const double time   = timelist.get<double>("time");
          
          // --------------------------------------------------
          // set parameters for stabilisation
          ParameterList& stablist = params.sublist("STABILIZATION");
          
          // if not available, define map from string to action
          if(stabstrtoact_.empty())
          {
            stabstrtoact_["quasistatic"    ]=subscales_quasistatic;
            stabstrtoact_["time_dependent" ]=subscales_time_dependent;
            stabstrtoact_["no_transient"   ]=inertia_stab_drop;
            stabstrtoact_["yes_transient"  ]=inertia_stab_keep;
            stabstrtoact_["no_pspg"        ]=pstab_assume_inf_sup_stable;
            stabstrtoact_["yes_pspg"       ]=pstab_use_pspg;
            stabstrtoact_["no_supg"        ]=convective_stab_none;
            stabstrtoact_["yes_supg"       ]=convective_stab_supg;
            stabstrtoact_["no_vstab"       ]=viscous_stab_none;
            stabstrtoact_["vstab_gls"      ]=viscous_stab_gls;
            stabstrtoact_["vstab_gls_rhs"  ]=viscous_stab_gls_only_rhs;
            stabstrtoact_["vstab_usfem"    ]=viscous_stab_usfem;
            stabstrtoact_["vstab_usfem_rhs"]=viscous_stab_usfem_only_rhs;
            stabstrtoact_["no_cstab"       ]=continuity_stab_none;
            stabstrtoact_["cstab_qs"       ]=continuity_stab_yes;
            stabstrtoact_["cstab_td"       ]=continuity_stab_td;
            stabstrtoact_["no_cross"       ]=cross_stress_stab_none;
            stabstrtoact_["cross_complete" ]=cross_stress_stab;
            stabstrtoact_["cross_rhs"      ]=cross_stress_stab_only_rhs;
            stabstrtoact_["no_reynolds"    ]=reynolds_stress_stab_none;
            stabstrtoact_["reynolds_rhs"   ]=reynolds_stress_stab_only_rhs;
          }
          
          StabilisationAction tds      = ConvertStringToStabAction(stablist.get<string>("TDS"));
          
          blitz::Array<double, 1>  mean_res(3);
          blitz::Array<double, 1>  mean_sacc(3);
          blitz::Array<double, 1>  mean_res_sq(3);
          blitz::Array<double, 1>  mean_sacc_sq(3);

          mean_res    =0;
          mean_sacc   =0;
          mean_res_sq =0;
          mean_sacc_sq=0;

          RefCountPtr<vector<double> > incrres       = params.get<RefCountPtr<vector<double> > >("incrres");
          RefCountPtr<vector<double> > incrres_sq    = params.get<RefCountPtr<vector<double> > >("incrres_sq");
          RefCountPtr<vector<double> > incrsacc      = params.get<RefCountPtr<vector<double> > >("incrsacc");
          RefCountPtr<vector<double> > incrsacc_sq   = params.get<RefCountPtr<vector<double> > >("incrsacc_sq");
          
          RefCountPtr<vector<double> > planecoords   = params.get<RefCountPtr<vector<double> > >("planecoords_");
          
          // --------------------------------------------------
          // calculate element coefficient matrix
          GenalphaResVMM()->CalcRes(this,
                                    evelnp,
                                    eprenp,
                                    eaccam,
                                    evelaf,
                                    actmat,
                                    alphaM,
                                    alphaF,
                                    gamma,
                                    dt,
                                    time,
                                    tds,
                                    mean_res,
                                    mean_sacc,
                                    mean_res_sq,
                                    mean_sacc_sq);
          

          //this will be the y-coordinate of a point in the element interior
          double center = 0;
          
          // get node coordinates of element
          blitz::Array<double,2>  xyze(3,NumNode());
          for(int inode=0;inode<NumNode();inode++)
          {
            xyze(0,inode)=Nodes()[inode]->X()[0];
            xyze(1,inode)=Nodes()[inode]->X()[1];
            xyze(2,inode)=Nodes()[inode]->X()[2];
            
            center+=xyze(1,inode);
          }
          center/=NumNode();

          bool found = false;
          
          int nlayer;
          for (nlayer=0;nlayer<(int)(*planecoords).size()-1;)
          {
            if(center<(*planecoords)[nlayer+1])
            {
              found = true;
              break;
            }
            nlayer++;
          }
          if (found ==false)
          {
              dserror("could not determine element layer");
          }


          for(int mm=0;mm<3;++mm)
          {
            (*incrres)    [3*nlayer+mm] += mean_res    (mm);
            (*incrres_sq) [3*nlayer+mm] += mean_res_sq (mm);
            (*incrsacc)   [3*nlayer+mm] += mean_sacc   (mm);
            (*incrsacc_sq)[3*nlayer+mm] += mean_sacc_sq(mm);
          }
        }
      }
      break;
      case calc_fluid_stationary_systemmat_and_residual:
      {
          // need current velocity/pressure 
          RefCountPtr<const Epetra_Vector> velnp = discretization.GetState("velnp");
          if (velnp==null)
            dserror("Cannot get state vector 'velnp'");

          // extract local values from the global vector
          vector<double> myvelnp(lm.size());
          DRT::UTILS::ExtractMyValues(*velnp,myvelnp,lm);

          if (is_ale_)
          {
        	  dserror("No ALE support within stationary fluid solver.");
          }

          // split velocity and pressure
          // create blitz objects

          const int numnode = NumNode();
          blitz::Array<double, 1> eprenp(numnode);
          blitz::Array<double, 2> evelnp(3,numnode,blitz::ColumnMajorArray<2>());

          for (int i=0;i<numnode;++i)
          {
            evelnp(0,i) = myvelnp[0+(i*4)];
            evelnp(1,i) = myvelnp[1+(i*4)];
            evelnp(2,i) = myvelnp[2+(i*4)];

            eprenp(i) = myvelnp[3+(i*4)];
          }

          // get control parameter
          const double pseudotime = params.get<double>("total time",-1.0);
          if (pseudotime < 0.0)
        	  dserror("no value for total (pseudo-)time in the parameter list");

          const bool newton = params.get<bool>("include reactive terms for linearisation",false);

          // --------------------------------------------------
          // set parameters for stabilisation
          ParameterList& stablist = params.sublist("STABILIZATION");

          // if not available, define map from string to action
          if(stabstrtoact_.empty())
          {
            stabstrtoact_["no_pspg"        ]=pstab_assume_inf_sup_stable;
            stabstrtoact_["yes_pspg"       ]=pstab_use_pspg;
            stabstrtoact_["no_supg"        ]=convective_stab_none;
            stabstrtoact_["yes_supg"       ]=convective_stab_supg;
            stabstrtoact_["no_vstab"       ]=viscous_stab_none;
            stabstrtoact_["vstab_gls"      ]=viscous_stab_gls;
            stabstrtoact_["vstab_gls_rhs"  ]=viscous_stab_gls_only_rhs;
            stabstrtoact_["vstab_usfem"    ]=viscous_stab_usfem;
            stabstrtoact_["vstab_usfem_rhs"]=viscous_stab_usfem_only_rhs;
            stabstrtoact_["no_cstab"       ]=continuity_stab_none;
            stabstrtoact_["cstab_qs"       ]=continuity_stab_yes;
            stabstrtoact_["no_cross"       ]=cross_stress_stab_none;
            stabstrtoact_["cross_complete" ]=cross_stress_stab;
            stabstrtoact_["cross_rhs"      ]=cross_stress_stab_only_rhs;
            stabstrtoact_["no_reynolds"    ]=reynolds_stress_stab_none;
            stabstrtoact_["reynolds_rhs"   ]=reynolds_stress_stab_only_rhs;
          }

          StabilisationAction pspg     = ConvertStringToStabAction(stablist.get<string>("PSPG"));
          StabilisationAction supg     = ConvertStringToStabAction(stablist.get<string>("SUPG"));
          StabilisationAction vstab    = ConvertStringToStabAction(stablist.get<string>("VSTAB"));
          StabilisationAction cstab    = ConvertStringToStabAction(stablist.get<string>("CSTAB"));
          StabilisationAction cross    = ConvertStringToStabAction(stablist.get<string>("CROSS-STRESS"));
          StabilisationAction reynolds = ConvertStringToStabAction(stablist.get<string>("REYNOLDS-STRESS"));

          // get flag for (fine-scale) subgrid viscosity (1=yes, 0=no)
          const int fssgv = params.get<int>("fs subgrid viscosity",0);

          // get fine-scale velocity
          RCP<const Epetra_Vector> fsvelnp;
          blitz::Array<double, 2> fsevelnp(3,numnode,blitz::ColumnMajorArray<2>());
          if (fssgv > 0)
          {
            fsvelnp = discretization.GetState("fsvelnp");
            if (fsvelnp==null) dserror("Cannot get state vector 'fsvelnp'");
            vector<double> myfsvelnp(lm.size());
            DRT::UTILS::ExtractMyValues(*fsvelnp,myfsvelnp,lm);

            // get fine-scale velocity and insert into element arrays
            for (int i=0;i<numnode;++i)
            {
              fsevelnp(0,i) = myfsvelnp[0+(i*4)];
              fsevelnp(1,i) = myfsvelnp[1+(i*4)];
              fsevelnp(2,i) = myfsvelnp[2+(i*4)];
            }
          }
          else
          {
            for (int i=0;i<numnode;++i)
            {
              fsevelnp(0,i) = 0.0;
              fsevelnp(1,i) = 0.0;
              fsevelnp(2,i) = 0.0;
            }
          }

          // get Smagorinsky model parameter for fine-scale subgrid viscosity
          const double Cs_fs = params.get<double>("fs Smagorinsky parameter",0.0);

          // wrap epetra serial dense objects in blitz objects
          blitz::Array<double, 2> estif(elemat1.A(),
                                        blitz::shape(elemat1.M(),elemat1.N()),
                                        blitz::neverDeleteData,
                                        blitz::ColumnMajorArray<2>());
          blitz::Array<double, 1> eforce(elevec1.Values(),
                                         blitz::shape(elevec1.Length()),
                                         blitz::neverDeleteData);

          // calculate element coefficient matrix and rhs         
          StationaryImpl()->Sysmat(this,
                                   evelnp,
                                   fsevelnp,
                                   eprenp,
                                   estif,
                                   eforce,
                                   actmat,
                                   pseudotime,
                                   newton,
                                   fssgv,
                                   pspg,
                                   supg,
                                   vstab,
                                   cstab,
                                   cross,
                                   reynolds,
                                   Cs_fs);

          // This is a very poor way to transport the density to the
          // outside world. Is there a better one?
          params.set("density", actmat->m.fluid->density);	  
      }
      break;
      default:
        dserror("Unknown type of action for Fluid3");
  } // end of switch(act)

  return 0;
} // end of DRT::ELEMENTS::Fluid3::Evaluate


/*----------------------------------------------------------------------*
 |  do nothing (public)                                      gammi 04/07|
 |                                                                      |
 |  The function is just a dummy. For the fluid elements, the           |
 |  integration of the volume neumann loads takes place in the element. |
 |  We need it there for the stabilisation terms!                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Fluid3::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1)
{
  return 0;
}

// get optimal gaussrule for discretization type
GaussRule3D DRT::ELEMENTS::Fluid3::getOptimalGaussrule(const DiscretizationType& distype)
{
    GaussRule3D rule = intrule3D_undefined;
    switch (distype)
    {
    case hex8:
        rule = intrule_hex_8point;
        break;
    case hex20: case hex27:
        rule = intrule_hex_27point;
        break;
    case tet4:
        rule = intrule_tet_4point;
        break;
    case tet10:
        rule = intrule_tet_5point;
        break;
    default:
        dserror("unknown number of nodes for gaussrule initialization");
  }
  return rule;
}


/*---------------------------------------------------------------------*
 |  calculate error for beltrami test problem               gammi 04/07|
 *---------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3::f3_int_beltrami_err(
  vector<double>&           evelnp,
  vector<double>&           eprenp,
  struct _MATERIAL*         material,
  ParameterList& 	    params
  )
{
  const int NSD = 3;

  // add element error to "integrated" error
  double velerr = params.get<double>("L2 integrated velocity error");
  double preerr = params.get<double>("L2 integrated pressure error");

  // set element data
  const int iel = NumNode();
  const DiscretizationType distype = this->Shape();

  Epetra_SerialDenseVector  funct(iel);
  Epetra_SerialDenseMatrix  xjm(3,3);
  Epetra_SerialDenseMatrix  deriv(3,iel);

  // get node coordinates of element
  Epetra_SerialDenseMatrix xyze(3,iel);
  for(int inode=0;inode<iel;inode++)
  {
    xyze(0,inode)=Nodes()[inode]->X()[0];
    xyze(1,inode)=Nodes()[inode]->X()[1];
    xyze(2,inode)=Nodes()[inode]->X()[2];
  }

  // set constants for analytical solution
  const double t = params.get("total time",-1.0);
  dsassert (t >= 0.0, "beltrami: no total time for error calculation");

  const double a      = PI/4.0;
  const double d      = PI/2.0;

  // get viscosity
  const double  visc = material->m.fluid->viscosity;

  double         preint;
  vector<double> velint  (3);
  vector<double> xint    (3);

  vector<double> u       (3);

  double         deltap;
  vector<double> deltavel(3);

  // gaussian points
  const GaussRule3D gaussrule = getOptimalGaussrule(distype);
  const IntegrationPoints3D  intpoints = getIntegrationPoints3D(gaussrule);

  // start loop over integration points
  for (int iquad=0;iquad<intpoints.nquad;iquad++)
  {
    // declaration of gauss point variables
    const double e1 = intpoints.qxg[iquad][0];
    const double e2 = intpoints.qxg[iquad][1];
    const double e3 = intpoints.qxg[iquad][2];
    shape_function_3D(funct,e1,e2,e3,distype);
    shape_function_3D_deriv1(deriv,e1,e2,e3,distype);

    /*----------------------------------------------------------------------*
      | calculate Jacobian matrix and it's determinant (private) gammi  07/07|
      | Well, I think we actually compute its transpose....
      |
      |     +-            -+ T      +-            -+
      |     | dx   dx   dx |        | dx   dy   dz |
      |     | --   --   -- |        | --   --   -- |
      |     | dr   ds   dt |        | dr   dr   dr |
      |     |              |        |              |
      |     | dy   dy   dy |        | dx   dy   dz |
      |     | --   --   -- |   =    | --   --   -- |
      |     | dr   ds   dt |        | ds   ds   ds |
      |     |              |        |              |
      |     | dz   dz   dz |        | dx   dy   dz |
      |     | --   --   -- |        | --   --   -- |
      |     | dr   ds   dt |        | dt   dt   dt |
      |     +-            -+        +-            -+
      |
      *----------------------------------------------------------------------*/
    Epetra_SerialDenseMatrix    xjm(NSD,NSD);

    for (int isd=0; isd<NSD; isd++)
    {
      for (int jsd=0; jsd<NSD; jsd++)
      {
        double dum = 0.0;
        for (int inode=0; inode<iel; inode++)
        {
          dum += deriv(isd,inode)*xyze(jsd,inode);
        }
        xjm(isd,jsd) = dum;
      }
    }

    // determinant of jacobian matrix
    const double det = xjm(0,0)*xjm(1,1)*xjm(2,2)+
                       xjm(0,1)*xjm(1,2)*xjm(2,0)+
                       xjm(0,2)*xjm(1,0)*xjm(2,1)-
                       xjm(0,2)*xjm(1,1)*xjm(2,0)-
                       xjm(0,0)*xjm(1,2)*xjm(2,1)-
                       xjm(0,1)*xjm(1,0)*xjm(2,2);

    if(det < 0.0)
    {
        printf("\n");
        printf("GLOBAL ELEMENT NO.%i\n",Id());
        printf("NEGATIVE JACOBIAN DETERMINANT: %lf\n", det);
        dserror("Stopped not regulary!\n");
    }

    const double fac = intpoints.qwgt[iquad]*det;

    // get velocity sol at integration point
    for (int i=0;i<3;i++)
    {
      velint[i]=0.0;
      for (int j=0;j<iel;j++)
      {
        velint[i] += funct[j]*evelnp[i+(3*j)];
      }
    }

    // get pressure sol at integration point
    preint = 0;
    for (int inode=0;inode<iel;inode++)
    {
      preint += funct[inode]*eprenp[inode];
    }

    // get velocity sol at integration point
    for (int isd=0;isd<3;isd++)
    {
      xint[isd]=0.0;
      for (int inode=0;inode<iel;inode++)
      {
        xint[isd] += funct[inode]*xyze(isd,inode);
      }
    }

    // compute analytical pressure
    const double p = -a*a/2.0 *
        ( exp(2.0*a*xint[0])
        + exp(2.0*a*xint[1])
        + exp(2.0*a*xint[2])
        + 2.0 * sin(a*xint[0] + d*xint[1]) * cos(a*xint[2] + d*xint[0]) * exp(a*(xint[1]+xint[2]))
        + 2.0 * sin(a*xint[1] + d*xint[2]) * cos(a*xint[0] + d*xint[1]) * exp(a*(xint[2]+xint[0]))
        + 2.0 * sin(a*xint[2] + d*xint[0]) * cos(a*xint[1] + d*xint[2]) * exp(a*(xint[0]+xint[1]))
        )* exp(-2.0*visc*d*d*t);

    // compute analytical velocities
    u[0] = -a * ( exp(a*xint[0]) * sin(a*xint[1] + d*xint[2]) +
                  exp(a*xint[2]) * cos(a*xint[0] + d*xint[1]) ) * exp(-visc*d*d*t);
    u[1] = -a * ( exp(a*xint[1]) * sin(a*xint[2] + d*xint[0]) +
                  exp(a*xint[0]) * cos(a*xint[1] + d*xint[2]) ) * exp(-visc*d*d*t);
    u[2] = -a * ( exp(a*xint[2]) * sin(a*xint[0] + d*xint[1]) +
                  exp(a*xint[1]) * cos(a*xint[2] + d*xint[0]) ) * exp(-visc*d*d*t);

    // compute difference between analytical solution and numerical solution
    deltap = preint - p;

    for (int isd=0;isd<NSD;isd++)
    {
      deltavel[isd] = velint[isd]-u[isd];
    }

    // add square to L2 error
    for (int isd=0;isd<NSD;isd++)
    {
      velerr += deltavel[isd]*deltavel[isd]*fac;
    }
    preerr += deltap*deltap*fac;

  } // end of loop over integration points


  // we use the parameterlist as a container to transport the calculated
  // errors from the elements to the dynamic routine

  params.set<double>("L2 integrated velocity error",velerr);
  params.set<double>("L2 integrated pressure error",preerr);

  return;
}

/*---------------------------------------------------------------------*
 | Calculate spatial mean values for channel flow (cartesian mesh)
 |                                                           gammi 07/07
 |
 | The necessary element integration is performed in here. The element
 | is cut into at least two (HEX8) or three (quadratic elements) planes,
 | the spatial functions (velocity, pressure etc.) are integrated over
 | this plane and this element contribution is added to a processor local
 | vector (see formulas below for a exact description of the output).
 | The method assumes, that all elements are of the same rectangular
 | shape in the "inplanedirection". In addition, it is assumed that
 | the sampling planes are distributed equidistant in the element.
 |
 |
 |                      ^ normdirect       integration plane
 |                      |                /
 |                      |               /
 |                      |
 |                +-----|-------------+
 |               /|     |            /|
 |              / |     |           / |
 |             /  |     |          /  |
 |            /   |     |         /   |
 |           /    +-----|--------/----+ ---- additional integration
 |          /    /|     |       /    /|      plane (for quadratic elements)
 |         /    / |     |      /    / |
 |        +-------------------+    /  |
 |        |   /   |     *-----|---+------------>
 |        |  /    +----/------|--/----+         inplanedirect[1]
 |        | /    /    /       | /    /
 |        |/    /    /        |/    /   \
 |        +---------+---------+    /     \
 |        |   /    /          |   /       integration plane
 |        |  /    /           |  /
 |        | /    /            | /
 |        |/    /             |/
 |        +----/--------------+
 |            /
 |           /   inplanedirect[0]
 |
 |
 |  Example for a mean value evaluation:
 |
 |         1.0       /
 |  _               |
 |  u = -------- *  | u(x,y,z) dx dy dz =
 |      +---        |
 |       \         / A
 |       / area
 |      +---
 |
 |
 |        1.0      /
 |                |            area
 |  =  -------- * | u(r,s,t) * ---- dr ds dt
 |     +---       |              4
 |      \        /  [-1:1]^2
 |      / area
 |     +---
 |
 |
 |
 |         1.0      /
 |                 |            1
 |  =   -------- * | u(r,s,t) * - dr ds dt
 |                 |            4
 |       numele   /  [-1:1]^2
 |
 |                |                        |
 |                +------------------------+
 |             this is the integral we compute!
 |
 | The factor 1/4 is necessary since we use a reference element of
 | size 2x2
 |
 | The method computes:
 |                      _             _             _             _
 |             numele * u  , numele * v  , numele * w  , numele * p
 |                      ___           ___           ___           ___
 |                       ^2            ^2            ^2            ^2
 | and         numele * u  , numele * v  , numele * w  , numele * p
 |                      _ _           _ _           _ _  
 | as well as  numele * u*v, numele * u*w, numele * v*w
 | 
 | as well as numele.
 | All results are communicated via the parameter list!
 |
 *---------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3::f3_calc_means(
  vector<double>&           sol  ,
  ParameterList& 	    params
  )
{

  // set element data
  const int iel = NumNode();
  const DiscretizationType distype = this->Shape();


  // the plane normal tells you in which plane the integration takes place
  const int normdirect = params.get<int>("normal direction to homogeneous plane");


  // the vector planes contains the coordinates of the homogeneous planes (in
  // wall normal direction)
  RefCountPtr<vector<double> > planes = params.get<RefCountPtr<vector<double> > >("coordinate vector for hom. planes");

  // get the pointers to the solution vectors
  RefCountPtr<vector<double> > sumu   = params.get<RefCountPtr<vector<double> > >("mean velocity u");
  RefCountPtr<vector<double> > sumv   = params.get<RefCountPtr<vector<double> > >("mean velocity v");
  RefCountPtr<vector<double> > sumw   = params.get<RefCountPtr<vector<double> > >("mean velocity w");
  RefCountPtr<vector<double> > sump   = params.get<RefCountPtr<vector<double> > >("mean pressure p");

  RefCountPtr<vector<double> > sumsqu = params.get<RefCountPtr<vector<double> > >("mean value u^2");
  RefCountPtr<vector<double> > sumsqv = params.get<RefCountPtr<vector<double> > >("mean value v^2");
  RefCountPtr<vector<double> > sumsqw = params.get<RefCountPtr<vector<double> > >("mean value w^2");
  RefCountPtr<vector<double> > sumuv  = params.get<RefCountPtr<vector<double> > >("mean value uv");
  RefCountPtr<vector<double> > sumuw  = params.get<RefCountPtr<vector<double> > >("mean value uw");
  RefCountPtr<vector<double> > sumvw  = params.get<RefCountPtr<vector<double> > >("mean value vw");
  RefCountPtr<vector<double> > sumsqp = params.get<RefCountPtr<vector<double> > >("mean value p^2");


  // get node coordinates of element
  Epetra_SerialDenseMatrix xyze(3,iel);
  for(int inode=0;inode<iel;inode++)
  {
    xyze(0,inode)=Nodes()[inode]->X()[0];
    xyze(1,inode)=Nodes()[inode]->X()[1];
    xyze(2,inode)=Nodes()[inode]->X()[2];
  }

  double min = xyze(normdirect,0);
  double max = xyze(normdirect,0);

  // set maximum and minimum value in wall normal direction
  for(int inode=0;inode<iel;inode++)
  {
    if(min > xyze(normdirect,inode))
    {
      min=xyze(normdirect,inode);
    }
    if(max < xyze(normdirect,inode))
    {
      max=xyze(normdirect,inode);
    }
  }

  // determine the ids of the homogeneous planes intersecting this element
  set<int> planesinele;
  for(unsigned nplane=0;nplane<planes->size();++nplane)
  {
    // get all available wall normal coordinates
    for(int nn=0;nn<iel;++nn)
    {
      if (min-2e-9 < (*planes)[nplane] && max+2e-9 > (*planes)[nplane])
      {
        planesinele.insert(nplane);
      }
    }
  }

  // remove lowest layer from planesinele to avoid double calculations. This is not done
  // for the first level (index 0) --- if deleted, shift the first integration point in
  // wall normal direction
  // the shift depends on the number of sampling planes in the element
  double shift=0;

  // set the number of planes which cut the element
  const int numplanesinele = planesinele.size();

  if(*planesinele.begin() != 0)
  {
    // this is not an element of the lowest element layer
    planesinele.erase(planesinele.begin());

    shift=2.0/((double) numplanesinele - 1.0);
  }
  else
  {
    // this is an element of the lowest element layer. Increase the counter
    // in order to compute the total number of elements in one layer
    int* count = params.get<int*>("count processed elements");

    (*count)++;
  }

  // determine the orientation of the rst system compared to the xyz system
  int elenormdirect=-1;
  bool upsidedown =false;
  // the only thing of interest is how normdirect is oriented in the
  // element coordinate system
  if(xyze(normdirect,4)-xyze(normdirect,0)>2e-9)
  {
    // t aligned
    elenormdirect =2;
    cout << "upsidedown false" <<&endl;
  }
  else if (xyze(normdirect,3)-xyze(normdirect,0)>2e-9)
  {
    // s aligned
    elenormdirect =1;
  }
  else if (xyze(normdirect,1)-xyze(normdirect,0)>2e-9)
  {
    // r aligned
    elenormdirect =0;
  }
  else if(xyze(normdirect,4)-xyze(normdirect,0)<-2e-9)
  {
    cout << xyze(normdirect,4)-xyze(normdirect,0) << &endl;
    // -t aligned
    elenormdirect =2;
    upsidedown =true;
    cout << "upsidedown true" <<&endl;
  }
  else if (xyze(normdirect,3)-xyze(normdirect,0)<-2e-9)
  {
    // -s aligned
    elenormdirect =1;
    upsidedown =true;
  }
  else if (xyze(normdirect,1)-xyze(normdirect,0)<-2e-9)
  {
    // -r aligned
    elenormdirect =0;
    upsidedown =true;
  }
  else
  {
    dserror("cannot determine orientation of plane normal in local coordinate system of element");
  }
  vector<int> inplanedirect;
  {
    set <int> inplanedirectset;
    for(int i=0;i<3;++i)
    {
      inplanedirectset.insert(i);
    }
    inplanedirectset.erase(elenormdirect);

    for(set<int>::iterator id = inplanedirectset.begin();id!=inplanedirectset.end() ;++id)
    {
      inplanedirect.push_back(*id);
    }
  }

  // allocate vector for shapefunctions
  Epetra_SerialDenseVector  funct(iel);

  // get the quad9 gaussrule for the in plane integration
  const IntegrationPoints2D  intpoints = getIntegrationPoints2D(intrule_quad_9point);

  // a hex8 element has two levels, the hex20 and hex27 element have three layers to sample
  // (now we allow even more)
  double layershift=0;

  // loop all levels in element
  for(set<int>::const_iterator id = planesinele.begin();id!=planesinele.end() ;++id)
  {
    // reset temporary values
    double ubar=0;
    double vbar=0;
    double wbar=0;
    double pbar=0;

    double usqbar=0;
    double vsqbar=0;
    double wsqbar=0;
    double uvbar =0;
    double uwbar =0;
    double vwbar =0;
    double psqbar=0;

    // get the intgration point in wall normal direction
    double e[3];

    e[elenormdirect]=-1.0+shift+layershift;
    if(upsidedown)
    {
      e[elenormdirect]*=-1;
    }

    // start loop over integration points in layer
    for (int iquad=0;iquad<intpoints.nquad;iquad++)
    {
      // get the other gauss point coordinates
      for(int i=0;i<2;++i)
      {
        e[inplanedirect[i]]=intpoints.qxg[iquad][i];
      }

      // compute the shape function values
      shape_function_3D(funct,e[0],e[1],e[2],distype);

      // check whether this gausspoint is really inside the desired plane
      {
        double x[3];
        x[0]=0;
        x[1]=0;
        x[2]=0;
        for(int inode=0;inode<iel;inode++)
        {
          x[0]+=funct[inode]*xyze(0,inode);
          x[1]+=funct[inode]*xyze(1,inode);
          x[2]+=funct[inode]*xyze(2,inode);
        }

        if(abs(x[normdirect]-(*planes)[*id])>2e-9)
        {
          dserror("Mixing up element cut planes during integration");
        }
      }

      //interpolated values at gausspoints
      double ugp=0;
      double vgp=0;
      double wgp=0;
      double pgp=0;

      // we assume that every 2d element we are integrating here is of
      // rectangular shape and every element is of the same size.
      // 1/4 is necessary since we use a reference element of size 2x2
      // the factor fac is omitting the element area up to now
      double fac=0.25*intpoints.qwgt[iquad];

      for(int inode=0;inode<iel;inode++)
      {
        ugp += funct[inode]*sol[inode*4  ];
        vgp += funct[inode]*sol[inode*4+1];
        wgp += funct[inode]*sol[inode*4+2];
        pgp += funct[inode]*sol[inode*4+3];
      }

      // add contribution to integral
      ubar   += ugp*fac;
      vbar   += vgp*fac;
      wbar   += wgp*fac;
      pbar   += pgp*fac;

      usqbar += ugp*ugp*fac;
      vsqbar += vgp*vgp*fac;
      wsqbar += wgp*wgp*fac;
      uvbar  += ugp*vgp*fac;
      uwbar  += ugp*wgp*fac;
      vwbar  += vgp*wgp*fac;
      psqbar += pgp*pgp*fac;
    } // end loop integration points

    // add increments from this layer to processor local vectors
    (*sumu  )[*id] += ubar;
    (*sumv  )[*id] += vbar;
    (*sumw  )[*id] += wbar;
    (*sump  )[*id] += pbar;

    (*sumsqu)[*id] += usqbar;
    (*sumsqv)[*id] += vsqbar;
    (*sumsqw)[*id] += wsqbar;
    (*sumuv) [*id] += uvbar;
    (*sumuw) [*id] += uwbar;
    (*sumvw) [*id] += vwbar;
    (*sumsqp)[*id] += psqbar;

    // jump to the next layer in the element.
    // in case of an hex8 element, the two coordinates are -1 and 1(+2)
    // for quadratic elements with three sample planes, we have -1,0(+1),1(+2)

    layershift+=2.0/((double) numplanesinele - 1.0);
  }


  return;
} // DRT::ELEMENTS::Fluid3::f3_calc_means

//----------------------------------------------------------------------
//
//----------------------------------------------------------------------
void DRT::ELEMENTS::Fluid3::f3_apply_box_filter(
    blitz::Array<double, 2>&  evelaf,
    blitz::Array<double, 1>&  vel_hat,
    blitz::Array<double, 2>&  reystr_hat,
    blitz::Array<double, 2>&  modeled_stress_grid_scale_hat,
    double&                   volume
    )
{

  //------------------------------------------------------------------
  //                     BLITZ CONFIGURATION
  //------------------------------------------------------------------
  //
  // We define the variables i,j,k to be indices to blitz arrays.
  // These are used for array expressions, that is matrix-vector
  // products in the following.

  blitz::firstIndex  i;   // Placeholder for the first index
  blitz::secondIndex j;   // Placeholder for the second index
  blitz::thirdIndex  k;   // Placeholder for the third index

  // set element data
  const int iel = NumNode();
  const DiscretizationType distype = this->Shape();

  // allocate arrays for shapefunctions, derivatives and the transposed jacobian
  blitz::Array<double,1>  funct(iel);
  blitz::Array<double,2>  xjm  (3,3);
  blitz::Array<double,2>  deriv(3,iel,blitz::ColumnMajorArray<2>());


  // get node coordinates of element
  blitz::Array<double,2>  xyze(3,iel);
  for(int inode=0;inode<iel;inode++)
  {
    xyze(0,inode)=Nodes()[inode]->X()[0];
    xyze(1,inode)=Nodes()[inode]->X()[1];
    xyze(2,inode)=Nodes()[inode]->X()[2];
  }

  // use one point gauss rule to calculate tau at element center
  DRT::UTILS::GaussRule3D integrationrule_filter=DRT::UTILS::intrule_hex_1point;
  switch (distype)
  {
      case DRT::Element::hex8:
        integrationrule_filter = DRT::UTILS::intrule_hex_1point;
        break;
      case DRT::Element::tet4:
        integrationrule_filter = DRT::UTILS::intrule_tet_1point;
        break;
      case DRT::Element::tet10:
      case DRT::Element::hex20:
      case DRT::Element::hex27:
        dserror("the box filtering operation is only permitted for linear elements\n");
        break;
      default:
        dserror("invalid discretization type for fluid3");
  }

  // gaussian points
  const DRT::UTILS::IntegrationPoints3D intpoints_onepoint(integrationrule_filter);
  
  // shape functions and derivs at element center
  const double e1    = intpoints_onepoint.qxg[0][0];
  const double e2    = intpoints_onepoint.qxg[0][1];
  const double e3    = intpoints_onepoint.qxg[0][2];
  const double wquad = intpoints_onepoint.qwgt[0];
  
  DRT::UTILS::shape_function_3D       (funct,e1,e2,e3,distype);
  DRT::UTILS::shape_function_3D_deriv1(deriv,e1,e2,e3,distype);

  // get Jacobian matrix and determinant
  xjm = blitz::sum(deriv(i,k)*xyze(j,k),k);
  const double det = xjm(0,0)*xjm(1,1)*xjm(2,2)+
                     xjm(0,1)*xjm(1,2)*xjm(2,0)+
                     xjm(0,2)*xjm(1,0)*xjm(2,1)-
                     xjm(0,2)*xjm(1,1)*xjm(2,0)-
                     xjm(0,0)*xjm(1,2)*xjm(2,1)-
                     xjm(0,1)*xjm(1,0)*xjm(2,2);

  //
  //             compute global first derivates
  //
  blitz::Array<double,2>  derxy(3,iel,blitz::ColumnMajorArray<2>());
  /*
    Use the Jacobian and the known derivatives in element coordinate
    directions on the right hand side to compute the derivatives in
    global coordinate directions

          +-                 -+     +-    -+      +-    -+
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |     | ---- |      | ---- |
          |  dr    dr    dr   |     |  dx  |      |  dr  |
          |                   |     |      |      |      |
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |  *  | ---- |   =  | ---- | for all k
          |  ds    ds    ds   |     |  dy  |      |  ds  |
          |                   |     |      |      |      |
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |     | ---- |      | ---- |
          |  dt    dt    dt   |     |  dz  |      |  dt  |
          +-                 -+     +-    -+      +-    -+

          Do one LU factorisation, everything else is backward substitution!

  */
  
  {
    // LAPACK solver
    Epetra_LAPACK          solver;
    
    // this copy of xjm will be used to calculate a in place factorisation
    blitz::Array<double,2> factorU(3,3,blitz::ColumnMajorArray<2>());
    factorU=xjm.copy();
    
    // a vector specifying the pivots (reordering)
    int pivot[3];

    // error code
    int ierr = 0;

    // Perform LU factorisation
    solver.GETRF(3,3,factorU.data(),3,&(pivot[0]),&ierr);
    
    if (ierr!=0)
    {
      dserror("Unable to perform LU factorisation during computation of derxy");
    }
    
    // backward substitution. The copy is required since GETRS replaces
    // the input with the result
    derxy =deriv.copy();
    solver.GETRS('N',3,iel,factorU.data(),3,&(pivot[0]),derxy.data(),3,&ierr);
    
    if (ierr!=0)
    {
      dserror("Unable to perform backward substitution after factorisation of jacobian");
    }
  }
  
  // get velocities (n+alpha_F,i) at integration point
  //
  //                 +-----
  //       n+af       \                  n+af
  //    vel    (x) =   +      N (x) * vel
  //                  /        j         j
  //                 +-----
  //                 node j
  //
  blitz::Array<double,1> velintaf (3);
  velintaf = blitz::sum(funct(j)*evelaf(i,j),j);


  // get velocity (n+alpha_F,i) derivatives at integration point
  //
  //       n+af      +-----  dN (x)
  //   dvel    (x)    \        k         n+af
  //   ----------- =   +     ------ * vel
  //       dx         /        dx        k
  //         j       +-----      j
  //                 node k
  //
  // j : direction of derivative x/y/z
  //
  blitz::Array<double,2>  vderxyaf(3,3);
  vderxyaf = blitz::sum(derxy(j,k)*evelaf(i,k),k);

  /*                            
                            +-     n+af          n+af    -+
          / h \       1.0   |  dvel_i  (x)   dvel_j  (x)  |
     eps | u   |    = --- * |  ----------- + -----------  |
          \   / ij    2.0   |      dx            dx       |
                            +-       j             i     -+
  */
  blitz::Array<double,2> epsilon(3,3,blitz::ColumnMajorArray<2>());
  epsilon = 0.5 * ( vderxyaf(i,j) + vderxyaf(j,i) );
  
  //
  // modeled part of subgrid scale stresses
  // 
  /*    +-                                 -+ 1                 
        |          / h \           / h \    | -         / h \   
        | 2 * eps | u   |   * eps | u   |   | 2  * eps | u   |  
        |          \   / kl        \   / kl |           \   / ij
        +-                                 -+                   
      
        |                                   |
        +-----------------------------------+
             'resolved' rate of strain
  */
  
  double rateofstrain = 0;
  
  for(int rr=0;rr<3;rr++)
  {
    for(int mm=0;mm<3;mm++)
    {
      rateofstrain += epsilon(rr,mm)*epsilon(rr,mm);
    }
  }
  rateofstrain *= 2.0;
  rateofstrain = sqrt(rateofstrain);

  
  //--------------------------------------------------
  // one point integrations

  // determine contribution to patch volume
  volume = wquad*det;

  // add contribution to integral over velocities
  vel_hat += velintaf*wquad*det;
  
  // add contribution to integral over reynolds stresses
  reystr_hat += velintaf(i)*velintaf(j)*wquad*det;
  
  // add contribution to integral over the modeled part of subgrid
  // scale stresses 
  modeled_stress_grid_scale_hat += rateofstrain * epsilon * wquad*det;

  return;
} // DRT::ELEMENTS::Fluid3::f3_apply_box_filter

//----------------------------------------------------------------------
// Calculate the quantities LijMij and MijMij, to compare the influence
// of the modeled and resolved stress tensor --- from this relation, Cs
// will be computed
//----------------------------------------------------------------------


void DRT::ELEMENTS::Fluid3::f3_calc_smag_const_LijMij_and_MijMij(
  blitz::Array<double, 2> evel_hat,
  blitz::Array<double, 3> ereynoldsstress_hat,
  blitz::Array<double, 3> efiltered_modeled_subgrid_stress_hat,
  double&                 LijMij,
  double&                 MijMij,
  double&                 center)
{

  //------------------------------------------------------------------
  //                     BLITZ CONFIGURATION
  //------------------------------------------------------------------
  //
  // We define the variables i,j,k to be indices to blitz arrays.
  // These are used for array expressions, that is matrix-vector
  // products in the following.

  blitz::firstIndex  i;   // Placeholder for the first index
  blitz::secondIndex j;   // Placeholder for the second index
  blitz::thirdIndex  k;   // Placeholder for the third index

  // set element data
  const int iel = NumNode();
  const DiscretizationType distype = this->Shape();

  // allocate arrays for shapefunctions, derivatives and the transposed jacobian
  blitz::Array<double,1>  funct(iel);
  blitz::Array<double,2>  xjm  (3,3);
  blitz::Array<double,2>  deriv(3,iel,blitz::ColumnMajorArray<2>());


  //this will be the y-coordinate of a point in the element interior
  center = 0;
  
  // get node coordinates of element
  blitz::Array<double,2>  xyze(3,iel);
  for(int inode=0;inode<iel;inode++)
  {
    xyze(0,inode)=Nodes()[inode]->X()[0];
    xyze(1,inode)=Nodes()[inode]->X()[1];
    xyze(2,inode)=Nodes()[inode]->X()[2];

    center+=xyze(1,inode);
  }
  center/=iel;
  

  // use one point gauss rule to calculate tau at element center
  DRT::UTILS::GaussRule3D integrationrule_filter=DRT::UTILS::intrule_hex_1point;
  switch (distype)
  {
      case DRT::Element::hex8:
      case DRT::Element::hex20:
      case DRT::Element::hex27:
        integrationrule_filter = DRT::UTILS::intrule_hex_1point;
        break;
      case DRT::Element::tet4:
      case DRT::Element::tet10:
        integrationrule_filter = DRT::UTILS::intrule_tet_1point;
        break;
      default:
        dserror("invalid discretization type for fluid3");
  }

  // gaussian points --- i.e. the midpoint
  const DRT::UTILS::IntegrationPoints3D intpoints_onepoint(integrationrule_filter);
  const double e1    = intpoints_onepoint.qxg[0][0];
  const double e2    = intpoints_onepoint.qxg[0][1];
  const double e3    = intpoints_onepoint.qxg[0][2];
  
  // shape functions and derivs at element center
  DRT::UTILS::shape_function_3D       (funct,e1,e2,e3,distype);
  DRT::UTILS::shape_function_3D_deriv1(deriv,e1,e2,e3,distype);
  
  // get element type constant for tau
  double mk=0.0;
  switch (distype)
  {
      case DRT::Element::tet4:
      case DRT::Element::hex8:
        mk = 0.333333333333333333333;
        break;
      case DRT::Element::hex20:
      case DRT::Element::hex27:
      case DRT::Element::tet10:
        mk = 0.083333333333333333333;
        break;
      default:
        dserror("type unknown!\n");
  }

  // get Jacobian matrix 
  xjm = blitz::sum(deriv(i,k)*xyze(j,k),k);

  //
  //             compute global first derivates
  //
  blitz::Array<double,2>  derxy(3,iel,blitz::ColumnMajorArray<2>());
  /*
    Use the Jacobian and the known derivatives in element coordinate
    directions on the right hand side to compute the derivatives in
    global coordinate directions

          +-                 -+     +-    -+      +-    -+
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |     | ---- |      | ---- |
          |  dr    dr    dr   |     |  dx  |      |  dr  |
          |                   |     |      |      |      |
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |  *  | ---- |   =  | ---- | for all k
          |  ds    ds    ds   |     |  dy  |      |  ds  |
          |                   |     |      |      |      |
          |  dx    dy    dz   |     | dN_k |      | dN_k |
          |  --    --    --   |     | ---- |      | ---- |
          |  dt    dt    dt   |     |  dz  |      |  dt  |
          +-                 -+     +-    -+      +-    -+

          Do one LU factorisation, everything else is backward substitution!

  */
  
  {
    // LAPACK solver
    Epetra_LAPACK          solver;
    
    // this copy of xjm will be used to calculate a in place factorisation
    blitz::Array<double,2> factorU(3,3,blitz::ColumnMajorArray<2>());
    factorU=xjm.copy();
    
    // a vector specifying the pivots (reordering)
    int pivot[3];

    // error code
    int ierr = 0;

    // Perform LU factorisation
    solver.GETRF(3,3,factorU.data(),3,&(pivot[0]),&ierr);
    
    if (ierr!=0)
    {
      dserror("Unable to perform LU factorisation during computation of derxy");
    }
    
    // backward substitution. The copy is required since GETRS replaces
    // the input with the result
    derxy =deriv.copy();
    solver.GETRS('N',3,iel,factorU.data(),3,&(pivot[0]),derxy.data(),3,&ierr);
    
    if (ierr!=0)
    {
      dserror("Unable to perform backward substitution after factorisation of jacobian");
    }
  }
  
  // get velocities (n+alpha_F/1,i) at integration point
  //
  //                   +-----
  //     ^ n+af/1       \                ^ n+af/1
  //    vel      (x) =   +      N (x) * vel
  //                     /        j         j
  //                    +-----
  //                    node j
  //
  blitz::Array<double,1> velint_hat (3);
  velint_hat = blitz::sum(funct(j)*evel_hat(i,j),j);


  // get velocity (n+alpha_F,i) derivatives at integration point
  //
  //     ^ n+af/1      +-----  dN (x)
  //   dvel      (x)    \        k       ^ n+af/1
  //   ------------- =   +     ------ * vel
  //       dx           /        dx        k
  //         j         +-----      j
  //                   node k
  //
  // j : direction of derivative x/y/z
  //
  blitz::Array<double,2>  vderxy_hat(3,3);
  vderxy_hat = blitz::sum(derxy(j,k)*evel_hat(i,k),k);

  // get filtered reynolds stress (n+alpha_F/1,i) at integration point
  //
  //                        +-----
  //        ^   n+af/1       \                   ^   n+af/1
  //    restress      (x) =   +      N (x) * restress
  //            ij           /        k              k, ij
  //                        +-----
  //                        node k
  //
  blitz::Array<double,2>  restress_hat(3,3);
  restress_hat = blitz::sum(funct(k)*ereynoldsstress_hat(i,j,k),k);

  // get filtered modeled subgrid stress (n+alpha_F/1,i) at integration point
  //
  //                                                   
  //                   ^                   n+af/1      
  //    filtered_modeled_subgrid_stress_hat      (x) = 
  //                                       ij          
  //
  //            +-----                                                      
  //             \                              ^                   n+af/1  
  //          =   +      N (x) * filtered_modeled_subgrid_stress_hat        
  //             /        k                                         k, ij   
  //            +-----                                                      
  //            node k                                                                                             
  //                                                   
  //
  blitz::Array<double,2>  filtered_modeled_subgrid_stress_hat(3,3);
  filtered_modeled_subgrid_stress_hat = blitz::sum(funct(k)*efiltered_modeled_subgrid_stress_hat(i,j,k),k);

  
  /*                            
                            +-   ^ n+af/1        ^   n+af/1    -+
      ^   / h \       1.0   |  dvel_i    (x)   dvel_j      (x)  |
     eps | u   |    = --- * |  ------------- + ---------------  |
          \   / ij    2.0   |       dx              dx          |
                            +-        j               i        -+
  */
  blitz::Array<double,2> epsilon_hat(3,3,blitz::ColumnMajorArray<2>());
  epsilon_hat = 0.5 * ( vderxy_hat(i,j) + vderxy_hat(j,i) );
  
  //
  // modeled part of subtestfilter scale stresses
  // 
  /*    +-                                 -+ 1                 
        |      ^   / h \       ^   / h \    | -     ^   / h \   
        | 2 * eps | u   |   * eps | u   |   | 2  * eps | u   |  
        |          \   / kl        \   / kl |           \   / ij
        +-                                 -+                   
      
        |                                   |
        +-----------------------------------+
             'resolved' rate of strain
  */
  
  double rateofstrain_hat = 0;
  
  for(int rr=0;rr<3;rr++)
  {
    for(int mm=0;mm<3;mm++)
    {
      rateofstrain_hat += epsilon_hat(rr,mm)*epsilon_hat(rr,mm);
    }
  }
  rateofstrain_hat *= 2.0;
  rateofstrain_hat = sqrt(rateofstrain_hat);
  
  blitz::Array<double, 2> L_ij (3,3,blitz::ColumnMajorArray<2>());
  blitz::Array<double, 2> M_ij (3,3,blitz::ColumnMajorArray<2>());

  L_ij = restress_hat(i,j) - velint_hat(i)*velint_hat(j);

  // this is sqrt(3)
  double filterwidthratio = 1.73;

  M_ij = filtered_modeled_subgrid_stress_hat(i,j)
         -
         filterwidthratio*filterwidthratio*rateofstrain_hat*epsilon_hat(i,j);

  LijMij =0;
  MijMij =0;
  for(int rr=0;rr<3;rr++)
  {
    for(int mm=0;mm<3;mm++)
    {
      LijMij += L_ij(rr,mm)*M_ij(rr,mm);
      MijMij += M_ij(rr,mm)*M_ij(rr,mm);
    }
  }
  
  return;
} // DRT::ELEMENTS::Fluid3::f3_calc_smag_const_LijMij_and_MijMij

//
// check for higher order derivatives for shape functions
//
bool DRT::ELEMENTS::Fluid3::isHigherOrderElement(
  const DRT::Element::DiscretizationType  distype) const
{
  bool hoel = true;
  switch (distype)
  {
  case hex8: case hex20: case hex27: case tet10: case wedge15:
    hoel = true;
    break;
  case tet4: case wedge6: case pyramid5:
    hoel = false;
    break;
  default:
    dserror("distype unknown!");
  }
  return hoel;
}

//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*----------------------------------------------------------------------*
 |  init the element (public)                                mwgee 12/06|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Fluid3Register::Initialize(DRT::Discretization& dis)
{
  bool dofillcompleteagain = false;
  //-------------------- loop all my column elements and check rewinding
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    // get the actual element
    if (dis.lColElement(i)->Type() != DRT::Element::element_fluid3) continue;
    DRT::ELEMENTS::Fluid3* actele = dynamic_cast<DRT::ELEMENTS::Fluid3*>(dis.lColElement(i));
    if (!actele) dserror("cast to Fluid3* failed");
    
    const DRT::Element::DiscretizationType distype = actele->Shape();
    bool possiblytorewind = false;
    switch(distype)
    {
    case DRT::Element::hex8: case DRT::Element::hex20: case DRT::Element::hex27:
        possiblytorewind = true;
        break;
    case DRT::Element::tet4: case DRT::Element::tet10:
        possiblytorewind = true;
        break;
    case DRT::Element::wedge6: case DRT::Element::wedge15:
        possiblytorewind = true;
        break;
    case DRT::Element::pyramid5:
        possiblytorewind = true;
        break;
    default:
        dserror("invalid discretization type for fluid3");
    }
    
    if ( (possiblytorewind) && (!actele->donerewinding_) ) {
      const bool rewind = checkRewinding3D(actele);

      if (rewind) {
        if (distype==DRT::Element::tet4){
          int iel = actele->NumNode();
          int new_nodeids[iel];
          const int* old_nodeids;
          old_nodeids = actele->NodeIds();
          // rewinding of nodes to arrive at mathematically positive element
          new_nodeids[0] = old_nodeids[0];
          new_nodeids[1] = old_nodeids[2];
          new_nodeids[2] = old_nodeids[1];
          new_nodeids[3] = old_nodeids[3];
          actele->SetNodeIds(iel, new_nodeids);
        }
        else if (distype==DRT::Element::hex8){
          int iel = actele->NumNode();
          int new_nodeids[iel];
          const int* old_nodeids;
          old_nodeids = actele->NodeIds();
          // rewinding of nodes to arrive at mathematically positive element
          new_nodeids[0] = old_nodeids[4];
          new_nodeids[1] = old_nodeids[5];
          new_nodeids[2] = old_nodeids[6];
          new_nodeids[3] = old_nodeids[7];
          new_nodeids[4] = old_nodeids[0];
          new_nodeids[5] = old_nodeids[1];
          new_nodeids[6] = old_nodeids[2];
          new_nodeids[7] = old_nodeids[3];
          actele->SetNodeIds(iel, new_nodeids);
        }
        else if (distype==DRT::Element::wedge6){
          int iel = actele->NumNode();
          int new_nodeids[iel];
          const int* old_nodeids;
          old_nodeids = actele->NodeIds();
          // rewinding of nodes to arrive at mathematically positive element
          new_nodeids[0] = old_nodeids[3];
          new_nodeids[1] = old_nodeids[4];
          new_nodeids[2] = old_nodeids[5];
          new_nodeids[3] = old_nodeids[0];
          new_nodeids[4] = old_nodeids[1];
          new_nodeids[5] = old_nodeids[2];
          actele->SetNodeIds(iel, new_nodeids);
        }
        else if (distype == DRT::Element::pyramid5){
          int iel = actele->NumNode();
          int new_nodeids[iel];
          const int* old_nodeids;
          old_nodeids = actele->NodeIds();
          // rewinding of nodes to arrive at mathematically positive element
          new_nodeids[1] = old_nodeids[3];
          new_nodeids[3] = old_nodeids[1];
          // the other nodes can stay the same
          new_nodeids[0] = old_nodeids[0];
          new_nodeids[2] = old_nodeids[2];
          new_nodeids[4] = old_nodeids[4];
          actele->SetNodeIds(iel, new_nodeids);
        }
        else dserror("no rewinding scheme for this type of fluid3");
      }
      // process of rewinding done
      actele->donerewinding_ = true;
      dofillcompleteagain = true;
    }
  }
  // fill complete again to reconstruct element-node pointers,
  // but without element init, etc.
  if(dofillcompleteagain) dis.FillComplete(false,false,false);
  
  return 0;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
