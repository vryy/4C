  /*!----------------------------------------------------------------------
\file fluid_ele_calc_poro_p1_immersed.cpp

 \brief internal implementation of poro immersed fluid element (p1 poro fluid)

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

\level 1

*----------------------------------------------------------------------*/

#include "fluid_ele_poro.H"
#include "fluid_ele_immersed.H"
#include "fluid_ele_parameter_poro.H"
#include "fluid_ele_calc_poro_p1_immersed.H"

#include "../drt_mat/fluidporo.H"
#include "../drt_inpar/inpar_cell.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_immersed_problem/immersed_field_exchange_manager.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcPoroP1Immersed<distype> *
DRT::ELEMENTS::FluidEleCalcPoroP1Immersed<distype>::Instance( bool create )
{
  static FluidEleCalcPoroP1Immersed<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new FluidEleCalcPoroP1Immersed<distype>();
    }
  }
  else
  {
    if ( instance!=NULL )
      delete instance;
    instance = NULL;
  }
  return instance;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP1Immersed<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
    Instance( false );
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcPoroP1Immersed<distype>::FluidEleCalcPoroP1Immersed()
  : DRT::ELEMENTS::FluidEleCalcPoroP1<distype>::FluidEleCalcPoroP1()
{
  // initializations
  exchange_manager_ = DRT::ImmersedFieldExchangeManager::Instance();
  gp_iquad_=-1234.0;
  artifical_flow_mod_type =
      DRT::INPUT::IntegralValue<int>(
          DRT::Problem::Instance()->CellMigrationParams(),
          "ARTIFICIAL_FLUID_MOD");

  // modification type
  // in standard case simplified formulation is active.
  // every modification except for the cancelation of
  // the phi*div(v^f) term is multiplied with 0.0
  mod_type_multiplicator_ = 0.0;

  // in case of full modification the other modifications
  // are included by multiplication with 1.0
  if(artifical_flow_mod_type==INPAR::CELL::mod_full)
    mod_type_multiplicator_ = 1.0;

  immersedele_=NULL;
}


/*----------------------------------------------------------------------*
 * Evaluate                                                 rauch 08/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcPoroP1Immersed<distype>::Evaluate(DRT::ELEMENTS::Fluid*     ele,
                                                 DRT::Discretization &          discretization,
                                                 const std::vector<int> &       lm,
                                                 Teuchos::ParameterList&        params,
                                                 Teuchos::RCP<MAT::Material> &  mat,
                                                 Epetra_SerialDenseMatrix&      elemat1_epetra,
                                                 Epetra_SerialDenseMatrix&      elemat2_epetra,
                                                 Epetra_SerialDenseVector&      elevec1_epetra,
                                                 Epetra_SerialDenseVector&      elevec2_epetra,
                                                 Epetra_SerialDenseVector&      elevec3_epetra,
                                                 bool                           offdiag)
{

  Teuchos::RCP<const MAT::FluidPoro> actmat = Teuchos::rcp_static_cast<const MAT::FluidPoro>(mat);

  my_p::const_permeability_ = (actmat->PermeabilityFunction() == MAT::PAR::const_);

  DRT::ELEMENTS::FluidPoro* poroele = dynamic_cast<DRT::ELEMENTS::FluidPoro*>(ele);

  if(poroele)
    my_p::kintype_ = poroele->KinematicType();


   // get integration rule for fluid elements cut by structural boundary
  int num_gp_fluid_bound = DRT::Problem::Instance()->ImmersedMethodParams().get<int>("NUM_GP_FLUID_BOUND");
  int degree_gp_fluid_bound=3;
  if(num_gp_fluid_bound == 8)
    degree_gp_fluid_bound = 3;
  else if (num_gp_fluid_bound == 64)
    degree_gp_fluid_bound = 7;
  else if (num_gp_fluid_bound == 125)
    degree_gp_fluid_bound = 9;
  else if (num_gp_fluid_bound == 343)
    degree_gp_fluid_bound = 13;
  else if (num_gp_fluid_bound == 729)
    degree_gp_fluid_bound = 17;
  else if (num_gp_fluid_bound == 1000)
    degree_gp_fluid_bound = 19;
  else
    dserror("Invalid value for parameter NUM_GP_FLUID_BOUND (valid parameters are 8, 64, 125, 343, 729 and 1000).");

  const DRT::UTILS::GaussIntegration intpoints_fluid_bound(distype,degree_gp_fluid_bound);

  // different integration rule for fluid elements cut by structural boundary
  immersedele_ = dynamic_cast<DRT::ELEMENTS::FluidImmersedBase* >(ele);

  if(not immersedele_)
    dserror("you shouldn't be here");

  if(immersedele_->IsBoundaryImmersed())
  {
    if (not offdiag) //evaluate diagonal block (pure fluid block)
      return my_p1::Evaluate(    ele,
                                 discretization,
                                 lm,
                                 params,
                                 mat,
                                 elemat1_epetra,
                                 elemat2_epetra,
                                 elevec1_epetra,
                                 elevec2_epetra,
                                 elevec3_epetra,
                                 intpoints_fluid_bound);
      else // evaluate off diagonal block (coupling block)
        return my_p1::EvaluateOD(ele,
                                 discretization,
                                 lm,
                                 params,
                                 mat,
                                 elemat1_epetra,
                                 elemat2_epetra,
                                 elevec1_epetra,
                                 elevec2_epetra,
                                 elevec3_epetra,
                                 intpoints_fluid_bound);
  }
  else
  {
    if (not offdiag) //evaluate diagonal block (pure fluid block)
      return my_p1::Evaluate(    ele,
                                 discretization,
                                 lm,
                                 params,
                                 mat,
                                 elemat1_epetra,
                                 elemat2_epetra,
                                 elevec1_epetra,
                                 elevec2_epetra,
                                 elevec3_epetra,
                                 DRT::ELEMENTS::FluidEleCalc<distype,DRT::ELEMENTS::Fluid::none>::intpoints_);
    else // evaluate off diagonal block (coupling block)
      return my_p1::EvaluateOD(  ele,
                                 discretization,
                                 lm,
                                 params,
                                 mat,
                                 elemat1_epetra,
                                 elemat2_epetra,
                                 elevec1_epetra,
                                 elevec2_epetra,
                                 elevec3_epetra,
                                 DRT::ELEMENTS::FluidEleCalc<distype,DRT::ELEMENTS::Fluid::none>::intpoints_);
  }
} // FluidEleCalcPoroP1Immersed<distype>::Evaluate


/*-----------------------------------------------------------------------*
 *   Computation of porosity                                 rauch 08/15 |
 *-----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP1Immersed<distype>::ComputePorosity(
        Teuchos::ParameterList&             params,
        const double&                       press,
        const double&                       J,
        const int&                          gp,
        const LINALG::Matrix<my::nen_,1>&   shapfct,
        const LINALG::Matrix<my::nen_,1>*   myporosity,
        double&                             porosity,
        double*                             dphi_dp,
        double*                             dphi_dJ,
        double*                             dphi_dJdp,
        double*                             dphi_dJJ,
        double*                             dphi_dpp,
        bool                                save
        )
{
  // set current gp number
  gp_iquad_ = gp;

  return my_p1::ComputePorosity( params,
                                 press,
                                 J,
                                 gp,
                                 shapfct,
                                 myporosity,
                                 porosity,
                                 dphi_dp,
                                 dphi_dJ,
                                 dphi_dJdp,
                                 dphi_dJJ,
                                 dphi_dpp,
                                 save);
} // FluidEleCalcPoroP1Immersed<distype>::ComputePorosity


/*------------------------------------------------------------------------*
 *  Evaluate Pressure Equation                                 rauch 08/15|
 *------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP1Immersed<distype>::EvaluatePressureEquation(
    Teuchos::ParameterList&                       params,
    const double&                                 timefacfacpre,
    const double&                                 rhsfac,
    const double&                                 dphi_dp,
    const double&                                 dphi_dJ,
    const double&                                 dphi_dJdp,
    const double&                                 dphi_dpp,
    const LINALG::Matrix<my::nen_,1>*             eporositydot,
    const LINALG::Matrix<my::nen_,1>*             eporositydotn,
    const LINALG::Matrix<my::nen_,1>&             echist,
    const LINALG::Matrix<my::nsd_,my::nen_>&      dgradphi_dp,
    LINALG::Matrix<my::nen_, my::nen_*my::nsd_>&  estif_q_u,
    LINALG::Matrix<my::nen_,my::nen_>&            ppmat,
    LINALG::Matrix<my::nen_,1>&                   preforce
    )
{
  EvaluatePressureEquationNonTransient( params,
                                        timefacfacpre,
                                        rhsfac,
                                        dphi_dp,
                                        dphi_dJ,
                                        dphi_dJdp,
                                        dphi_dpp,
                                        dgradphi_dp,
                                        estif_q_u,
                                        ppmat,
                                        preforce);

  // now the porosity time derivative (different for standard poro and poro_p1 elements)
// now the porosity time derivative (different for standard poro and poro_p1 elements)
if (my_p::porofldpara_->IsStationaryConti() == false)
{
  if(eporositydot)
  {
    // initialize dphi/dt
    double porositydot = -12345.0;

    // we do not add the dphi/dt term in the artifical domain
    // in case of full modification
    if(immersedele_->IsImmersed() and
       artifical_flow_mod_type == INPAR::CELL::mod_full)
      porositydot = 0.0;
    else if(immersedele_->IsBoundaryImmersed() and
            artifical_flow_mod_type == INPAR::CELL::mod_full)
    {
      // check whether integration point is covered by cell
      bool gp_has_projected_divergence = false;
      if( (immersedele_->GetRCPIntPointHasProjectedDivergence()) != Teuchos::null)
        if(immersedele_->GetRCPIntPointHasProjectedDivergence()->size()>0)
          gp_has_projected_divergence =
              (int)(immersedele_->IntPointHasProjectedDivergence(gp_iquad_));

      if(gp_has_projected_divergence)
        porositydot = my_p::funct_.Dot(*eporositydot);
      else
        porositydot = my_p::funct_.Dot(*eporositydot);
    }
    else // physical domain
      porositydot = my_p::funct_.Dot(*eporositydot);

    // rhs transient term
    for (int vi=0; vi<my_p::nen_; ++vi)
      preforce(vi)-=  rhsfac * porositydot * my_p::funct_(vi) ;

    // just update internal variables, no contribution to rhs
    const double porositydotn = my_p::funct_.Dot(*eporositydotn);

    my_p::histcon_ = my_p::fldparatimint_->OmTheta() * my_p::fldparatimint_->Dt() * porositydotn;

    //rhs from last time step
    my_p::rhscon_ = 1.0/my_p::fldparatimint_->Dt()/my_p::fldparatimint_->Theta() * my_p::histcon_;

    //transient part of continuity equation residual
    my_p::conres_old_ += porositydot - my_p::rhscon_;
  }
  else
    dserror("no porosity time derivative given for poro_p1 element!");
}
else
  dserror("Stationary conti equation is not implemented and not tested for\n"
          "immersed cell-flow interaction.");

  return;
} // FluidEleCalcPoroP1Immersed<distype>::EvaluatePressureEquation


/*-------------------------------------------------------------------------*
 *  evaluation of continuity equation (non transient part)     rauch 08/15 |
 *-------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP1Immersed<distype>::
EvaluatePressureEquationNonTransient(
    Teuchos::ParameterList&                       params,
    const double&                                 timefacfacpre,
    const double&                                 rhsfac,
    const double&                                 dphi_dp,
    const double&                                 dphi_dJ,
    const double&                                 dphi_dJdp,
    const double&                                 dphi_dpp,
    const LINALG::Matrix<my::nsd_,my::nen_>&      dgradphi_dp,
    LINALG::Matrix<my::nen_, my::nen_*my::nsd_>&  estif_q_u,
    LINALG::Matrix<my::nen_,my::nen_>&            ppmat,
    LINALG::Matrix<my::nen_,1>&                   preforce
    )
{
  if (my::fldparatimint_->IsStationary())
    dserror("immersed interstitial flow - cell interaction only tested for transient problems.");

  my_p::EvaluatePressureEquationNonTransient( params,
                                              timefacfacpre,
                                              rhsfac,
                                              dphi_dp,
                                              dphi_dJ,
                                              dphi_dJdp,
                                              dphi_dpp,
                                              dgradphi_dp,
                                              estif_q_u,
                                              ppmat,
                                              preforce);


  // check whether integration point is covered by immersed cell
  bool gp_has_projected_divergence = false;
  if( (immersedele_->GetRCPIntPointHasProjectedDivergence()) != Teuchos::null)
    if(immersedele_->GetRCPIntPointHasProjectedDivergence()->size()>0)
      gp_has_projected_divergence =
          (int)(immersedele_->IntPointHasProjectedDivergence(gp_iquad_));

  // we substract all terms that should not be contained in the flow formulation
  // in the artificial domain.
  if(immersedele_->IsImmersed() and
      exchange_manager_->IsFluidInteraction() and
      artifical_flow_mod_type != INPAR::CELL::mod_none)
  {
    if( my_p::porofldpara_->PoroContiPartInt() == false )
    {
      double vel_grad_porosity = 0.0;
      for (int idim = 0; idim <my::nsd_; ++idim)
        vel_grad_porosity += my_p::grad_porosity_(idim)*my::velint_(idim);

      double    grad_porosity_gridvelint=0.0;
      for (int j =0; j< my::nsd_; j++)
        grad_porosity_gridvelint += my_p::grad_porosity_(j) * my_p::gridvelint_(j);

      // left-hand side
      for (int idim = 0; idim <my::nsd_; ++idim)
      {
        const double grad_porosity_idim = my_p::grad_porosity_(idim);
        for (int ui=0; ui<my::nen_; ++ui)
        {
          const int fui = my::nsd_*ui;
          const double funct_ui = my::funct_(ui);
          const double derxy_idim_ui = my::derxy_(idim,ui);

          for (int vi=0; vi<my::nen_; ++vi)
          {
            estif_q_u(vi,fui+idim) -= timefacfacpre*my::funct_(vi) * (my_p::porosity_ * derxy_idim_ui); // phi*div(v^f) (remove) Dv^f
            estif_q_u(vi,fui+idim) -= timefacfacpre*my::funct_(vi) * grad_porosity_idim * funct_ui * mod_type_multiplicator_; // grad(phi)*v^f (remove) Dv^f
          }
        }
      }

      // right-hand side
      const double rhsfac_vdiv = rhsfac * my::vdiv_;
      for (int vi=0; vi<my::nen_; ++vi)
      {
        preforce(vi) +=   rhsfac_vdiv * my_p::porosity_ * my::funct_(vi); // phi*div(v^f) (remove)
        preforce(vi) +=   rhsfac * vel_grad_porosity * my::funct_(vi) * mod_type_multiplicator_;    // grad(phi)*v^f (remove)
        preforce(vi) +=   rhsfac *my::funct_(vi) * (- grad_porosity_gridvelint ) * mod_type_multiplicator_; // grad(phi)*v^s (remove)
        //preforce(vi) -=   rhsfac *my::funct_(vi) * my_p::porosity_ * my_p::gridvdiv_ * mod_type_multiplicator_; // phi*div(v^s) (add)
      } // vi

    }
    else // my_p::porofldpara_->PoroContiPartInt() == true
      dserror("No adaption of interstitial flow to compressibility of immersed body implemented\n"
              "for partially integrated poro P1 formulation!\n"
              "Set CONTIPARTINT to 'No' in ---POROELASTICITY DYNAMIC section.");
  }
  else if (immersedele_->IsBoundaryImmersed() and
           gp_has_projected_divergence and
           exchange_manager_->IsFluidInteraction() and
           artifical_flow_mod_type != INPAR::CELL::mod_none)
  {
    for (int vi=0; vi<my::nen_; ++vi)
      preforce(vi) += rhsfac * immersedele_->ProjectedIntPointDivergence(gp_iquad_) * my_p::porosity_ * my::funct_(vi);

  } // if 'cut' element
  return;
} // FluidEleCalcPoroP1Immersed<distype>::EvaluatePressureEquationNonTransient


/*----------------------------------------------------------------------*
 |  gp loop for off-diagonal terms                          rauch 10/17 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcPoroP1Immersed<distype>::GaussPointLoopP1OD(
                        Teuchos::ParameterList&                                     params,
                        const LINALG::Matrix<my::nsd_,my::nen_>&                    ebofoaf,
                        const LINALG::Matrix<my::nsd_, my::nen_>&                   evelaf,
                        const LINALG::Matrix<my::nsd_, my::nen_>&                   evelnp,
                        const LINALG::Matrix<my::nsd_, my::nen_>&                   eveln,
                        const LINALG::Matrix<my::nen_, 1>&                          epreaf,
                        const LINALG::Matrix<my::nen_, 1>&                          eprenp,
                        const LINALG::Matrix<my::nen_, 1>&                          epren,
                        const LINALG::Matrix<my::nsd_,my::nen_> &                   emhist,
                        const LINALG::Matrix<my::nen_,1>&                           echist,
                        const LINALG::Matrix<my::nen_, 1> &                         epressnp_timederiv,
                        const LINALG::Matrix<my::nen_, 1> &                         epressam_timederiv,
                        const LINALG::Matrix<my::nen_, 1> &                         epressn_timederiv,
                        const LINALG::Matrix<my::nsd_,my::nen_>&                    eaccam,
                        const LINALG::Matrix<my::nsd_, my::nen_>&                   edispnp,
                        const LINALG::Matrix<my::nsd_, my::nen_>&                   edispn,
                        const LINALG::Matrix<my::nsd_, my::nen_>&                   egridv,
                        const LINALG::Matrix<my::nsd_, my::nen_>&                   egridvn,
                        const LINALG::Matrix<my::nen_,1>&                           escaaf,
                        const LINALG::Matrix<my::nen_,1>*                           eporositynp,
                        LINALG::Matrix<(my::nsd_ + 1) * my::nen_, 1>&               eforce,
                        LINALG::Matrix<my::nen_ * my::nsd_, my::nen_ * my::nsd_>&   ecoupl_u,
                        LINALG::Matrix<my::nen_, my::nen_ * my::nsd_>&              ecoupl_p,
                        LINALG::Matrix<my::nen_ * my::nsd_, my::nen_>&              ecouplp1_u,
                        LINALG::Matrix<my::nen_, my::nen_>&                         ecouplp1_p,
                        Teuchos::RCP<const MAT::Material>                           material,
                        const DRT::UTILS::GaussIntegration &                        intpoints)
{
  my_p1::GaussPointLoopP1OD(
      params,
      ebofoaf,
      evelaf,
      evelnp,
      eveln,
      epreaf,
      eprenp,
      epren,
      emhist,
      echist,
      epressnp_timederiv,
      epressam_timederiv,
      epressn_timederiv,
      eaccam,
      edispnp,
      edispn,
      egridv,
      egridvn,
      escaaf,
      eporositynp,
      eforce,
      ecoupl_u,
      ecoupl_p,
      ecouplp1_u,
      ecouplp1_p,
      material,
      intpoints);

  if(immersedele_->IsImmersed() and
     artifical_flow_mod_type != INPAR::CELL::mod_none)
  {
    for ( DRT::UTILS::GaussIntegration::const_iterator iquad=intpoints.begin(); iquad!=intpoints.end(); ++iquad )
    {
      // evaluate shape functions and derivatives at integration point
      my_p::EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(),iquad.Weight());
      // evaluate shape function derivatives w.r.t. to material coordinates at integration point
      my_p::SetupMaterialDerivatives();

      for (int ui=0; ui<my_p::nen_; ++ui)
        for (int vi=0; vi<my_p::nen_; ++vi)
          ecouplp1_p(vi,ui) -=  my_p::fac_ * my_p::funct_(vi) * my_p::funct_(ui) * mod_type_multiplicator_;

      static LINALG::Matrix<my_p::nen_,1> derxy_fluidvel(false);
      derxy_fluidvel.Clear();

      for (int i =0; i< my_p::nen_; i++)
        for (int j =0; j< my_p::nsd_; j++)
          derxy_fluidvel(i) += my_p::derxy_(j,i) * my::velint_(j);

      static LINALG::Matrix<my_p::nen_,1> derxy_gridvel(false);
      derxy_gridvel.Clear();

      for (int i =0; i< my_p::nen_; i++)
        for (int j =0; j< my_p::nsd_; j++)
          derxy_gridvel(i) += my_p::derxy_(j,i) * (-my_p::gridvelint_(j));

      double timefacfacpre = my_p::fldparatimint_->TimeFacPre() * my_p::fac_;

      for (int ui=0; ui<my_p::nen_; ++ui)
      {
        for (int vi=0; vi<my_p::nen_; ++vi)
        {
          ecouplp1_p(vi,ui) -=
            + timefacfacpre * my_p::vdiv_ * my_p::funct_(vi) * my_p::funct_(ui)
            + timefacfacpre * my_p::funct_(vi) * derxy_fluidvel(ui) * mod_type_multiplicator_
            + timefacfacpre * my_p::funct_(vi) * derxy_gridvel(ui)  * mod_type_multiplicator_
                                 ;
        } // vi
      } // ui
    } // gp loop
  } // if immersed ele
  return;
} // FluidEleCalcPoroP1Immersed<distype>::GaussPointLoopP1OD

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template class DRT::ELEMENTS::FluidEleCalcPoroP1Immersed<DRT::Element::hex8>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1Immersed<DRT::Element::hex20>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1Immersed<DRT::Element::hex27>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1Immersed<DRT::Element::tet4>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1Immersed<DRT::Element::tet10>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1Immersed<DRT::Element::wedge6>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1Immersed<DRT::Element::wedge15>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1Immersed<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1Immersed<DRT::Element::quad4>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1Immersed<DRT::Element::quad8>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1Immersed<DRT::Element::quad9>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1Immersed<DRT::Element::tri3>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1Immersed<DRT::Element::tri6>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1Immersed<DRT::Element::nurbs9>;
template class DRT::ELEMENTS::FluidEleCalcPoroP1Immersed<DRT::Element::nurbs27>;


