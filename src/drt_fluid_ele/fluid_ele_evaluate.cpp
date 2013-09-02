/*!
\file fluid_ele_evaluate.cpp
\brief

<pre>
Maintainer: Ursula Rasthofer & Volker Gravemeier
            {rasthofer,vgravem}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236/-245
</pre>
*/


#include "fluid_ele_factory.H"

#include "fluid_ele.H"
#include "fluid_ele_action.H"
#include "fluid_ele_evaluate_utils.H"
#include "fluid_ele_utils.H"
#include "fluid_ele_interface.H"
#include "fluid_ele_parameter.H"
#include "fluid_ele_tds.H"

#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"

#include "../drt_geometry/position_array.H"

#include "../drt_inpar/inpar_fluid.H"

#include "../drt_lib/drt_utils.H"

#include "../drt_mat/arrhenius_pv.H"
#include "../drt_mat/carreauyasuda.H"
#include "../drt_mat/ferech_pv.H"
#include "../drt_mat/herschelbulkley.H"
#include "../drt_mat/mixfrac.H"
#include "../drt_mat/modpowerlaw.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/permeablefluid.H"
#include "../drt_mat/sutherland.H"
#include "../drt_mat/yoghurt.H"

#include "../drt_nurbs_discret/drt_nurbs_discret.H"

#include "../drt_opti/topopt_fluidAdjoint3_interface.H"
#include "../drt_opti/topopt_fluidAdjoint3_impl_parameter.H"

#include "../drt_lib/drt_condition_utils.H"


/*
  Depending on the type of action and the element type (tet, hex etc.),
  the elements allocate common static arrays.

  */

/*---------------------------------------------------------------------*
|  Call the element to set all basic parameter                         |
*----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidType::PreEvaluate(DRT::Discretization&                  dis,
                                            Teuchos::ParameterList&               p,
                                            Teuchos::RCP<LINALG::SparseOperator>  systemmatrix1,
                                            Teuchos::RCP<LINALG::SparseOperator>  systemmatrix2,
                                            Teuchos::RCP<Epetra_Vector>           systemvector1,
                                            Teuchos::RCP<Epetra_Vector>           systemvector2,
                                            Teuchos::RCP<Epetra_Vector>           systemvector3)
{
  const FLD::Action action = DRT::INPUT::get<FLD::Action>(p,"action");

  if (action == FLD::set_general_fluid_parameter)
  {
    Teuchos::RCP<DRT::ELEMENTS::FluidEleParameter> fldpara = DRT::ELEMENTS::FluidEleParameter::Instance();
    fldpara->SetElementGeneralFluidParameter(p,dis.Comm().MyPID());
  }
  else if (action == FLD::set_time_parameter)
  {
    Teuchos::RCP<DRT::ELEMENTS::FluidEleParameter> fldpara = DRT::ELEMENTS::FluidEleParameter::Instance();
    fldpara->SetElementTimeParameter(p);
  }
  else if (action == FLD::set_turbulence_parameter)
  {
    Teuchos::RCP<DRT::ELEMENTS::FluidEleParameter> fldpara = DRT::ELEMENTS::FluidEleParameter::Instance();
    fldpara->SetElementTurbulenceParameter(p);
  }
  else if (action == FLD::set_loma_parameter)
  {
    Teuchos::RCP<DRT::ELEMENTS::FluidEleParameter> fldpara = DRT::ELEMENTS::FluidEleParameter::Instance();
    fldpara->SetElementLomaParameter(p);
  }
  else if (action == FLD::set_poro_parameter)
  {
    Teuchos::RCP<DRT::ELEMENTS::FluidEleParameter> fldpara = DRT::ELEMENTS::FluidEleParameter::Instance();
    fldpara->SetElementPoroParameter(p);
  }
  else if (action == FLD::set_topopt_parameter)
  {
    Teuchos::RCP<DRT::ELEMENTS::FluidEleParameter> fldpara = DRT::ELEMENTS::FluidEleParameter::Instance();
    fldpara->SetElementTopoptParameter(p);
  }
  else if (action == FLD::set_general_adjoint_parameter)
  {
    Teuchos::RCP<DRT::ELEMENTS::FluidAdjoint3ImplParameter> fldpara = DRT::ELEMENTS::FluidAdjoint3ImplParameter::Instance();
    fldpara->SetElementGeneralAdjointParameter(p);
  }
  else if (action == FLD::set_adjoint_time_parameter)
  {
    Teuchos::RCP<DRT::ELEMENTS::FluidAdjoint3ImplParameter> fldpara = DRT::ELEMENTS::FluidAdjoint3ImplParameter::Instance();
    fldpara->SetElementAdjointTimeParameter(p);
  }

  return;
}


 /*----------------------------------------------------------------------*
 |  evaluate the element (public)                            g.bau 03/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Fluid::Evaluate(Teuchos::ParameterList&            params,
                                    DRT::Discretization&      discretization,
                                    std::vector<int>&         lm,
                                    Epetra_SerialDenseMatrix& elemat1,
                                    Epetra_SerialDenseMatrix& elemat2,
                                    Epetra_SerialDenseVector& elevec1,
                                    Epetra_SerialDenseVector& elevec2,
                                    Epetra_SerialDenseVector& elevec3)
{
  // get the action required
  const FLD::Action act = DRT::INPUT::get<FLD::Action>(params,"action");

  // get material
  RCP<MAT::Material> mat = Material();

  // get space dimensions
  const int nsd = DRT::UTILS::getDimension(Shape());

  switch(act)
  {
    //-----------------------------------------------------------------------
    // standard implementation enabling time-integration schemes such as
    // one-step-theta, BDF2, and generalized-alpha (n+alpha_F and n+1)
    //-----------------------------------------------------------------------
    case FLD::calc_fluid_systemmat_and_residual:
    {
      switch(params.get<int>("physical type",INPAR::FLUID::incompressible))
      {
      case INPAR::FLUID::loma:
      {
          return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), "loma")->Evaluate(
              this,
              discretization,
              lm,
              params,
              mat,
              elemat1,
              elemat2,
              elevec1,
              elevec2,
              elevec3 );
          break;
      }
      case INPAR::FLUID::poro:
      {
        return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), "poro")->Evaluate(
            this,
            discretization,
            lm,
            params,
            mat,
            elemat1,
            elemat2,
            elevec1,
            elevec2,
            elevec3 );
        break;
      }
      case INPAR::FLUID::poro_p1:
      {
        return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), "poro_p1")->Evaluate(
            this,
            discretization,
            lm,
            params,
            mat,
            elemat1,
            elemat2,
            elevec1,
            elevec2,
            elevec3 );
        break;
      }
      case INPAR::FLUID::poro_p2:
      {
        return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), "poro_p2")->Evaluate(
            this,
            discretization,
            lm,
            params,
            mat,
            elemat1,
            elemat2,
            elevec1,
            elevec2,
            elevec3 );
        break;
      }
      default:
        return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), "std")->Evaluate(
            this,
            discretization,
            lm,
            params,
            mat,
            elemat1,
            elemat2,
            elevec1,
            elevec2,
            elevec3);
      }

    }
    break;
    //-----------------------------------------------------------------------
    // standard implementation enabling time-integration schemes such as
    // one-step-theta, BDF2, and generalized-alpha (n+alpha_F and n+1)
    // for the particular case of porous flow
    //-----------------------------------------------------------------------
    /***********************************************/
    case FLD::calc_porousflow_fluid_coupling:
    {
      switch(params.get<int>("physical type",INPAR::FLUID::incompressible))
      {
      case INPAR::FLUID::poro:
      {
        return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), "poro")->Evaluate(
            this,
            discretization,
            lm,
            params,
            mat,
            elemat1,
            elemat2,
            elevec1,
            elevec2,
            elevec3,
            true);
        break;
      }
      case INPAR::FLUID::poro_p1:
      {
        return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), "poro_p1")->Evaluate(
            this,
            discretization,
            lm,
            params,
            mat,
            elemat1,
            elemat2,
            elevec1,
            elevec2,
            elevec3,
            true);
        break;
      }
      case INPAR::FLUID::poro_p2:
      {
        return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), "poro_p2")->Evaluate(
            this,
            discretization,
            lm,
            params,
            mat,
            elemat1,
            elemat2,
            elevec1,
            elevec2,
            elevec3,
            true);
        break;
      }
      default:
        dserror("Unknown physical type for poroelasticity\n");
      break;
      }
    }
    break;
    //-----------------------------------------------------------------------
    // standard implementation enabling time-integration schemes such as
    // one-step-theta, BDF2, and generalized-alpha (n+alpha_F and n+1)
    // for evaluation of off-diagonal matrix block for monolithic
    // low-Mach-number solver
    //-----------------------------------------------------------------------
    case FLD::calc_loma_mono_odblock:
    {
      return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), "loma")->Evaluate(
          this,
          discretization,
          lm,
          params,
          mat,
          elemat1,
          elemat2,
          elevec1,
          elevec2,
          elevec3,
          true);
    }
    break;
    //-----------------------------------------------------------------------
    // standard implementation enabling time-integration schemes such as
    // one-step-theta, BDF2, and generalized-alpha (n+alpha_F and n+1)
    // for evaluation of off-diagonal matrix block for monolithic
    // porous flow with scalar transport
    //-----------------------------------------------------------------------
    case FLD::calc_poroscatra_mono_odblock:
    {
      switch(params.get<int>("physical type",INPAR::FLUID::incompressible))
      {
      case INPAR::FLUID::poro_p2:
      {
      return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), "poro_p2")->Evaluate(
          this,
          discretization,
          lm,
          params,
          mat,
          elemat1,
          elemat2,
          elevec1,
          elevec2,
          elevec3,
          true);
      }
      break;
      default:
        dserror("Unknown physical type for poroelasticity\n");
      break;
      }
    }
    break;
    case FLD::calc_div_u:
    {
      // compute divergence of velocity field at the element
      return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), "std")->ComputeDivU(this, discretization, lm, elevec1);
    }
    break;
    case FLD::calc_fluid_error:
    {
      // compute error for a known analytical solution
      return  DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), "std")->ComputeError(this, params, mat, discretization, lm, elevec1);
    }
    break;
    case FLD::calc_turbulence_statistics:
    {
      if (nsd == 3)
      {
        // do nothing if you do not own this element
        if(this->Owner() == discretization.Comm().MyPID())
        {
          // --------------------------------------------------
          // extract velocity and pressure from global
          // distributed vectors
          // --------------------------------------------------
          // velocity and pressure values (n+1)
          RCP<const Epetra_Vector> velnp
          = discretization.GetState("u and p (n+1,converged)");
          if (velnp==Teuchos::null) dserror("Cannot get state vector 'velnp'");

          // extract local values from the global vectors
          std::vector<double> mysol  (lm.size());
          DRT::UTILS::ExtractMyValues(*velnp,mysol,lm);

          std::vector<double> mydisp(lm.size());
          if(is_ale_)
          {
            // get most recent displacements
            RCP<const Epetra_Vector> dispnp
            =
                discretization.GetState("dispnp");

            if (dispnp==Teuchos::null)
            {
              dserror("Cannot get state vectors 'dispnp'");
            }

            DRT::UTILS::ExtractMyValues(*dispnp,mydisp,lm);
          }

          // integrate mean values
          const DiscretizationType distype = this->Shape();

          switch (distype)
          {
          case DRT::Element::hex8:
          {
            FLD::f3_calc_means<8>(this,discretization,mysol,mydisp,params);
            break;
          }
          case DRT::Element::hex20:
          {
            FLD::f3_calc_means<20>(this,discretization,mysol,mydisp,params);
            break;
          }
          case DRT::Element::hex27:
          {
            FLD::f3_calc_means<27>(this,discretization,mysol,mydisp,params);
            break;
          }
          case DRT::Element::nurbs8:
          {
            FLD::f3_calc_means<8>(this,discretization,mysol,mydisp,params);
            break;
          }
          case DRT::Element::nurbs27:
          {
            FLD::f3_calc_means<27>(this,discretization,mysol,mydisp,params);
            break;
          }
          default:
          {
            dserror("Unknown element type for mean value evaluation\n");
          }
          }
        }
      } // end if (nsd == 3)
      else dserror("action 'calc_turbulence_statistics' is a 3D specific action");
    }
    break;
    case FLD::calc_turbscatra_statistics:
    {
      if(nsd == 3)
      {
        // do nothing if you do not own this element
        if(this->Owner() == discretization.Comm().MyPID())
        {
          // --------------------------------------------------
          // extract velocity, pressure, and scalar from global
          // distributed vectors
          // --------------------------------------------------
          // velocity/pressure and scalar values (n+1)
          RCP<const Epetra_Vector> velnp
            = discretization.GetState("u and p (n+1,converged)");
          RCP<const Epetra_Vector> scanp
            = discretization.GetState("scalar (n+1,converged)");
          if (velnp==Teuchos::null || scanp==Teuchos::null)
            dserror("Cannot get state vectors 'velnp' and/or 'scanp'");

          // extract local values from the global vectors
          std::vector<double> myvelpre(lm.size());
          std::vector<double> mysca(lm.size());
          DRT::UTILS::ExtractMyValues(*velnp,myvelpre,lm);
          DRT::UTILS::ExtractMyValues(*scanp,mysca,lm);

          // integrate mean values
          const DiscretizationType distype = this->Shape();

          switch (distype)
          {
          case DRT::Element::hex8:
          {
            FLD::f3_calc_scatra_means<8>(this,discretization,myvelpre,mysca,params);
            break;
          }
          case DRT::Element::hex20:
          {
            FLD::f3_calc_scatra_means<20>(this,discretization,myvelpre,mysca,params);
            break;
          }
          case DRT::Element::hex27:
          {
            FLD::f3_calc_scatra_means<27>(this,discretization,myvelpre,mysca,params);
            break;
          }
          default:
          {
            dserror("Unknown element type for turbulent passive scalar mean value evaluation\n");
          }
          }
        }
      } // end if (nsd == 3)
      else dserror("action 'calc_turbscatra_statistics' is a 3D specific action");
    }
    break;
    case FLD::calc_loma_statistics:
    {
      if(nsd == 3)
      {
        // do nothing if you do not own this element
        if(this->Owner() == discretization.Comm().MyPID())
        {
          // --------------------------------------------------
          // extract velocity, pressure, and temperature from
          // global distributed vectors
          // --------------------------------------------------
          // velocity/pressure and scalar values (n+1)
          RCP<const Epetra_Vector> velnp
          = discretization.GetState("u and p (n+1,converged)");
          RCP<const Epetra_Vector> scanp
          = discretization.GetState("scalar (n+1,converged)");
          if (velnp==Teuchos::null || scanp==Teuchos::null)
            dserror("Cannot get state vectors 'velnp' and/or 'scanp'");

          // extract local values from global vectors
          std::vector<double> myvelpre(lm.size());
          std::vector<double> mysca(lm.size());
          DRT::UTILS::ExtractMyValues(*velnp,myvelpre,lm);
          DRT::UTILS::ExtractMyValues(*scanp,mysca,lm);

          // get factor for equation of state
          const double eosfac = params.get<double>("eos factor",100000.0/287.0);

          // integrate mean values
          const DiscretizationType distype = this->Shape();

          switch (distype)
          {
          case DRT::Element::hex8:
          {
            FLD::f3_calc_loma_means<8>(this,discretization,myvelpre,mysca,params,eosfac);
            break;
          }
          case DRT::Element::hex20:
          {
            FLD::f3_calc_loma_means<20>(this,discretization,myvelpre,mysca,params,eosfac);
            break;
          }
          case DRT::Element::hex27:
          {
            FLD::f3_calc_loma_means<27>(this,discretization,myvelpre,mysca,params,eosfac);
            break;
          }
          default:
          {
            dserror("Unknown element type for low-Mach-number mean value evaluation\n");
          }
          }
        }
      } // end if (nsd == 3)
      else dserror("action 'calc_loma_statistics' is a 3D specific action");
    }
    break;
    case FLD::calc_fluid_box_filter:
    {
      if (nsd == 3)
      {
        const DiscretizationType distype = this->Shape();
        int nen = 0;
        if (distype == DRT::Element::hex8)
            nen = 8;
        else if (distype == DRT::Element::tet4)
            nen = 4;
        else dserror("not supported");

        // --------------------------------------------------
        // extract velocity and pressure from global
        // distributed vectors
        // --------------------------------------------------
        // velocity, pressure and temperature values (most recent
        // intermediate solution, i.e. n+alphaF for genalpha
        // and n+1 for one-step-theta)
        RCP<const Epetra_Vector> vel =
            discretization.GetState("u and p (trial)");
        if (vel==Teuchos::null)
          dserror("Cannot get state vectors 'vel'");
        // extract local values from the global vectors
        std::vector<double> myvel(lm.size());
        DRT::UTILS::ExtractMyValues(*vel,myvel,lm);

        std::vector<double> tmp_temp(lm.size());
        std::vector<double> mytemp(nen);
        double thermpress = 0.0;
        // pointer to class FluidEleParameter (access to the general parameter)
        Teuchos::RCP<DRT::ELEMENTS::FluidEleParameter> fldpara = DRT::ELEMENTS::FluidEleParameter::Instance();
        if (fldpara->PhysicalType()==INPAR::FLUID::loma)
        {
          RCP<const Epetra_Vector> temp = discretization.GetState("T (trial)");
          if (temp==Teuchos::null)
            dserror("Cannot get state vectors 'temp'");
          DRT::UTILS::ExtractMyValues(*temp,tmp_temp,lm);

          for (int i=0;i<nen;i++)
             mytemp[i] = tmp_temp[nsd+(i*(nsd+1))];

          // get thermodynamic pressure
          thermpress = params.get<double>("thermpress");
        }

        // initialize the contribution of this element to the patch volume to zero
        double volume_contribution = 0.0;

        // initialize the contributions of this element to the filtered scalar quantities
        double dens_hat = 0.0;
        double dens_strainrate_hat = 0.0;
        // get pointers for vector quantities
        RCP<std::vector<double> > vel_hat = params.get<RCP<std::vector<double> > >("vel_hat");
        RCP<std::vector<double> > densvel_hat = params.get<RCP<std::vector<double> > >("densvel_hat");
        RCP<std::vector<std::vector<double> > > reynoldsstress_hat = params.get<RCP<std::vector<std::vector<double> > > >("reynoldsstress_hat");
        RCP<std::vector<std::vector<double> > > modeled_subgrid_stress = params.get<RCP<std::vector<std::vector<double> > > >("modeled_subgrid_stress");

        // integrate the convolution with the box filter function for this element
        // the results are assembled onto the *_hat arrays
        switch (distype)
        {
        case DRT::Element::hex8:
        {
          FLD::f3_apply_box_filter<8>(
              this,
              fldpara,
              myvel,
              mytemp,
              thermpress,
              vel_hat,
              densvel_hat,
              reynoldsstress_hat,
              modeled_subgrid_stress,
              volume_contribution,
              dens_hat,
              dens_strainrate_hat);
          break;
        }
        case DRT::Element::tet4:
        {
          FLD::f3_apply_box_filter<4>(
              this,
              fldpara,
              myvel,
              mytemp,
              thermpress,
              vel_hat,
              densvel_hat,
              reynoldsstress_hat,
              modeled_subgrid_stress,
              volume_contribution,
              dens_hat,
              dens_strainrate_hat);
          break;
        }
        default:
        {
          dserror("Unknown element type for box filter application\n");
        }
        }

        // hand down the volume contribution to the time integration algorithm
        params.set<double>("volume_contribution",volume_contribution);
        // as well as the filtered scalar quantities
        params.set<double>("dens_hat",dens_hat);
        params.set<double>("dens_strainrate_hat",dens_strainrate_hat);

      } // end if (nsd == 3)
      else dserror("action 'calc_fluid_box_filter' is 3D specific action");
    }
    break;
    case FLD::calc_smagorinsky_const:
    {
      if(nsd == 3)
      {
        RCP<Epetra_MultiVector> col_filtered_vel                    =
            params.get<RCP<Epetra_MultiVector> >("col_filtered_vel");
        RCP<Epetra_MultiVector> col_filtered_reynoldsstress         =
            params.get<RCP<Epetra_MultiVector> >("col_filtered_reynoldsstress");
        RCP<Epetra_MultiVector> col_filtered_modeled_subgrid_stress =
            params.get<RCP<Epetra_MultiVector> >("col_filtered_modeled_subgrid_stress");

        // pointer to class FluidEleParameter (access to the general parameter)
        Teuchos::RCP<DRT::ELEMENTS::FluidEleParameter> fldpara = DRT::ELEMENTS::FluidEleParameter::Instance();
        // add potential loma specific vectors
        RCP<Epetra_MultiVector> col_filtered_dens_vel = Teuchos::null;
        RCP<Epetra_Vector> col_filtered_dens = Teuchos::null;
        RCP<Epetra_Vector> col_filtered_dens_strainrate = Teuchos::null;
        if (fldpara->PhysicalType()==INPAR::FLUID::loma)
        {
          col_filtered_dens_vel = params.get<RCP<Epetra_MultiVector> >("col_filtered_dens_vel");
          col_filtered_dens = params.get<RCP<Epetra_Vector> >("col_filtered_dens");
          col_filtered_dens_strainrate = params.get<RCP<Epetra_Vector> >("col_filtered_dens_strainrate");
        }

        double LijMij   = 0.0;
        double MijMij   = 0.0;
        double CI_numerator = 0.0;
        double CI_denominator = 0.0;
        double xcenter  = 0.0;
        double ycenter  = 0.0;
        double zcenter  = 0.0;

        const DiscretizationType distype = this->Shape();
        switch (distype)
        {
        case DRT::Element::hex8:
        {
          FLD::f3_calc_smag_const_LijMij_and_MijMij<8>(
              this,
              fldpara,
              col_filtered_vel                   ,
              col_filtered_reynoldsstress        ,
              col_filtered_modeled_subgrid_stress,
              col_filtered_dens_vel              ,
              col_filtered_dens                  ,
              col_filtered_dens_strainrate       ,
              LijMij                             ,
              MijMij                             ,
              CI_numerator                       ,
              CI_denominator                     ,
              xcenter                            ,
              ycenter                            ,
              zcenter                            );
          break;
        }
        case DRT::Element::tet4:
        {
          FLD::f3_calc_smag_const_LijMij_and_MijMij<4>(
              this,
              fldpara,
              col_filtered_vel                   ,
              col_filtered_reynoldsstress        ,
              col_filtered_modeled_subgrid_stress,
              col_filtered_dens_vel              ,
              col_filtered_dens                  ,
              col_filtered_dens_strainrate       ,
              LijMij                             ,
              MijMij                             ,
              CI_numerator                       ,
              CI_denominator                     ,
              xcenter                            ,
              ycenter                            ,
              zcenter                            );
          break;
        }
        default:
        {
          dserror("Unknown element type for box filter application\n");
        }
        }

        double Cs_delta_sq = 0.0;
        double Ci_delta_sq = 0.0;
        // set Cs_delta_sq without averaging (only clipping)
        if (abs(MijMij) < 1E-16) Cs_delta_sq= 0.0;
        else  Cs_delta_sq = 0.5 * LijMij / MijMij;
        if (Cs_delta_sq<0.0)
        {
          Cs_delta_sq= 0.0;
        }

        if (fldpara->PhysicalType()==INPAR::FLUID::loma)
        {
          // and set Ci_delta_sq without averaging (only clipping) for loma
          if (abs(CI_denominator) < 1E-16) Ci_delta_sq= 0.0;
          else Ci_delta_sq = 0.5 * CI_numerator / CI_denominator;
          if (Ci_delta_sq<0.0)
          {
            Ci_delta_sq= 0.0;
          }
        }

        // set all values in parameter list
        params.set<double>("ele_Cs_delta_sq",Cs_delta_sq);
        params.set<double>("ele_Ci_delta_sq",Ci_delta_sq);
        params.set<double>("LijMij",LijMij);
        params.set<double>("MijMij",MijMij);
        params.set<double>("CI_numerator",CI_numerator);
        params.set<double>("CI_denominator",CI_denominator);
        params.set<double>("xcenter",xcenter);
        params.set<double>("ycenter",ycenter);
        params.set<double>("zcenter",zcenter);
      } // end if(nsd == 3)
      else dserror("action 'calc_smagorinsky_const' is a 3D specific action");
    }
    break;
    case FLD::calc_fluid_genalpha_update_for_subscales:
    {
      // time update for time-dependent subgrid-scales
      const double dt     = params.get<double>("dt");
      const double gamma  = params.get<double>("gamma");
      this->TDS()->Update(dt,gamma);
    }
    break;
    case FLD::calc_dissipation:
    {
      if (nsd == 3)
      {
        if (this->Owner() == discretization.Comm().MyPID()) // don't store values of ghosted elements
        {
          return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), "std")->CalcDissipation(
              this,
              params,
              discretization,
              lm,
              mat);
        }
      }
      else dserror("%i D elements does not support calculation of dissipation", nsd);
    }
    break;
    case FLD::calc_model_params_mfsubgr_scales:
    {
      if (nsd == 3)
      {
        // velocity values
        RCP<const Epetra_Vector> velnp = discretization.GetState("velnp");
        // fine-scale velocity values
        RCP<const Epetra_Vector> fsvelnp = discretization.GetState("fsvelnp");
        if (velnp==Teuchos::null or fsvelnp==Teuchos::null)
        {
          dserror("Cannot get state vectors");
        }

        // extract local values from the global vectors
        std::vector<double> myvel(lm.size());
        DRT::UTILS::ExtractMyValues(*velnp,myvel,lm);
        std::vector<double> myfsvel(lm.size());
        DRT::UTILS::ExtractMyValues(*fsvelnp,myfsvel,lm);

        // pointer to class FluidEleParameter (access to the general parameter)
        Teuchos::RCP<DRT::ELEMENTS::FluidEleParameter> fldpara = DRT::ELEMENTS::FluidEleParameter::Instance();

        const DiscretizationType distype = this->Shape();
        switch (distype)
        {
        case DRT::Element::hex8:
        {
          // don't store values of ghosted elements
          if (this->Owner() == discretization.Comm().MyPID())
          {
            FLD::f3_get_mf_params<8,3,DRT::Element::hex8>(this,fldpara,params,mat,myvel,myfsvel);
          }
          break;
        }
        default:
        {
          dserror("Unknown element type\n");
        }
        }
      }
      else dserror("%i D elements does not support calculation of model parameters", nsd);
    }
    break;
    case FLD::calc_mean_Cai:
    {
      if (nsd == 3)
      {
        const DiscretizationType distype = this->Shape();
        int nen = 0;
        if (distype == DRT::Element::hex8)
           nen = 8;
        else dserror("not supported");

        // velocity values
        RCP<const Epetra_Vector> vel = discretization.GetState("velocity");
        // scalar values
        RCP<const Epetra_Vector> sca = discretization.GetState("scalar");
        if (vel==Teuchos::null or sca==Teuchos::null)
        {
          dserror("Cannot get state vectors");
        }
        // extract local values from the global vectors
        std::vector<double> myvel(lm.size());
        DRT::UTILS::ExtractMyValues(*vel,myvel,lm);
        std::vector<double> tmp_sca(lm.size());
        std::vector<double> mysca(nen);
        DRT::UTILS::ExtractMyValues(*sca,tmp_sca,lm);
        for (int i=0;i<nen;i++)
           mysca[i] = tmp_sca[nsd+(i*(nsd+1))];
        // get thermodynamic pressure
        double thermpress = params.get<double>("thermpress");

        // pointer to class FluidEleParameter
        Teuchos::RCP<DRT::ELEMENTS::FluidEleParameter> fldpara = DRT::ELEMENTS::FluidEleParameter::Instance();

        double Cai = 0.0;
        double vol = 0.0;

        bool is_inflow_ele = false;

       std::vector<DRT::Condition*> myinflowcond;

       // check whether all nodes have a unique inflow condition
       DRT::UTILS::FindElementConditions(this, "TurbulentInflowSection", myinflowcond);
       if (myinflowcond.size()>1)
         dserror("More than one inflow condition on one node!");

       if (myinflowcond.size()==1)
         is_inflow_ele = true;

       // exclude elemenets of inflow section
        if (not is_inflow_ele)
        {
          switch (distype)
          {
            case DRT::Element::hex8:
            {
              FLD::f3_get_mf_nwc<8,3,DRT::Element::hex8>(this,fldpara,Cai,vol,myvel,mysca,thermpress);
              break;
            }
            default:
            {
              dserror("Unknown element type\n");
            }
          }
        }
        
        // hand down the Cai and volume contribution to the time integration algorithm
        params.set<double>("Cai_int",Cai);
        params.set<double>("ele_vol",vol);
      }
      else dserror("%i D elements does not support calculation of mean Cai", nsd);
    }
    break;
    case FLD::set_mean_Cai:
    {
      // pointer to class FluidEleParameter
      Teuchos::RCP<DRT::ELEMENTS::FluidEleParameter> fldpara = DRT::ELEMENTS::FluidEleParameter::Instance();
      fldpara->SetCsgsPhi(params.get<double>("meanCai"));
    }
    break;
    case FLD::calc_node_normal:
    {
      if (nsd == 3)
      {
        const DiscretizationType distype = this->Shape();
        switch (distype)
        {
        case DRT::Element::hex27:
        {
          FLD::ElementNodeNormal<DRT::Element::hex27>(this,params,discretization,lm,elevec1);
          break;
        }
        case DRT::Element::hex20:
        {
          FLD::ElementNodeNormal<DRT::Element::hex20>(this,params,discretization,lm,elevec1);
          break;
        }
        case DRT::Element::hex8:
        {
          FLD::ElementNodeNormal<DRT::Element::hex8>(this,params,discretization,lm,elevec1);
          break;
        }
        case DRT::Element::tet4:
        {
          FLD::ElementNodeNormal<DRT::Element::tet4>(this,params,discretization,lm,elevec1);
          break;
        }
        case DRT::Element::tet10:
        {
          FLD::ElementNodeNormal<DRT::Element::tet10>(this,params,discretization,lm,elevec1);
          break;
        }
        default:
        {
          dserror("Unknown element type for shape function integration\n");
        }
        }
      }
      else dserror("action 'calculate node normal' should also work in 2D, but 2D elements are not"
          " added to the template yet. Also it is not tested");
      break;
    }
    case FLD::integrate_shape:
    {
      // integrate shape function for this element
      // (results assembled into element vector)
      return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), "std")->IntegrateShapeFunction(this, discretization, lm, elevec1);
    }
    case FLD::calc_divop:
    {
      // calculate the integrated divergence operator
      return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), "std")->CalcDivOp(this, discretization, lm, elevec1);
    }
    case FLD::calc_mat_deriv_u_and_rot_u:
    {
      // calculate material derivative at specified element coordinates
      return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), "std")->CalcMatDerivAndRotU(
          this,
          params,
          discretization,
          lm,
          elevec1,
          elevec2,
          elevec3);
    }
    case FLD::void_fraction_gaussian_integration:
    {
      // calculate material derivative at specified element coordinates
      return DRT::ELEMENTS::FluidFactory::ProvideImpl(Shape(), "std")->ComputeVoidFraction(this, params, discretization, lm, elevec1);
    }
    case FLD::set_general_fluid_parameter:
    case FLD::set_time_parameter:
    case FLD::set_turbulence_parameter:
    case FLD::set_loma_parameter:
    case FLD::set_poro_parameter:
    case FLD::set_topopt_parameter:
    case FLD::set_general_adjoint_parameter:
    case FLD::set_adjoint_time_parameter:
//    case FLD::calc_adjoint_neumann: // this is done by the surface elements
      break;
    //-----------------------------------------------------------------------
    // adjoint implementation enabling time-integration schemes such as
    // one-step-theta, BDF2, and generalized-alpha (n+alpha_F and n+1)
    //-----------------------------------------------------------------------
    case FLD::calc_adjoint_systemmat_and_residual:
    {
      return DRT::ELEMENTS::FluidAdjoint3ImplInterface::Impl(Shape())->Evaluate(
          this,
          discretization,
          lm,
          params,
          mat,
          elemat1,
          elevec1 );
      break;
    }
    default:
      dserror("Unknown type of action for Fluid");
      break;
  } // end of switch(act)

  return 0;
} // end of DRT::ELEMENTS::Fluid::Evaluate


/*----------------------------------------------------------------------*
 |  do nothing (public)                                      gammi 04/07|
 |                                                                      |
 |  The function is just a dummy. For fluid elements, the integration   |
 |  integration of volume Neumann conditions (body forces) takes place  |
 |  in the element. We need it there for the stabilisation terms!       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Fluid::EvaluateNeumann(Teuchos::ParameterList&    params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           std::vector<int>&         lm,
                                           Epetra_SerialDenseVector& elevec1,
                                           Epetra_SerialDenseMatrix* elemat1)
{
  return 0;
}


