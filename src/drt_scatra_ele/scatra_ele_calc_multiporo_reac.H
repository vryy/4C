/*----------------------------------------------------------------------*/
/*!
 \file scatra_ele_calc_multiporo_reac.H

 \brief evaluation class containing routines for calculation of scalar transport
        within multiphase porous medium

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_SCATRA_ELE_SCATRA_ELE_CALC_MULTIPORO_REAC_H_
#define SRC_DRT_SCATRA_ELE_SCATRA_ELE_CALC_MULTIPORO_REAC_H_


#include "scatra_ele_calc_poro_reac.H"
#include "../drt_mat/fluidporo_multiphase.H"
#include "../drt_mat/fluidporo_singlephase.H"
#include "../drt_mat/fluidporo_multiphase_reactions.H"
#include "../drt_porofluidmultiphase_ele/porofluidmultiphase_ele_calc_utils.H"

#include "../drt_porofluidmultiphase_ele/porofluid_phasemanager.H"
#include "../drt_porofluidmultiphase_ele/porofluidmultiphase_ele_action.H"
#include "../drt_porofluidmultiphase_ele/porofluidmultiphase_ele_calc.H"
#include "../drt_porofluidmultiphase_ele/porofluid_variablemanager.H"
#include "../drt_porofluidmultiphase_ele/porofluid_evaluator.H"
#include "../drt_porofluidmultiphase_ele/porofluidmultiphase_ele_parameter.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"

//forward declarations
namespace MAT
{
class ScatraMat;
}

namespace DRT
{
namespace ELEMENTS
{

template<int NSD,int NEN> class ScaTraEleInternalVariableManagerMultiPoro;

template<DRT::Element::DiscretizationType distype>
class ScaTraEleCalcMultiPoroReac : public ScaTraEleCalcPoroReac<distype>
{
private:
  /// private constructor, since we are a Singleton.
  ScaTraEleCalcMultiPoroReac(const int numdofpernode,const int numscal,const std::string& disname);

  typedef ScaTraEleCalc<distype> my;
  typedef ScaTraEleCalcPoroReac<distype> pororeac;
  typedef ScaTraEleCalcPoro<distype> poro;
  typedef ScaTraEleCalcAdvReac<distype> advreac;

public:
  virtual ~ScaTraEleCalcMultiPoroReac() {}

  /// Singleton access method
  static ScaTraEleCalcMultiPoroReac<distype> * Instance(
    const int numdofpernode,
    const int numscal,
    const std::string& disname,
    const ScaTraEleCalcMultiPoroReac* delete_me = NULL );

  /// called upon destruction
  virtual void Done();

  /// Setup element evaluation
  virtual int SetupCalc(
      DRT::Element*               ele,
      DRT::Discretization&        discretization
      );

  protected:

  //! extract element based or nodal values
  //  return extracted values of phinp
  virtual void ExtractElementAndNodeValues(
      DRT::Element*                 ele,
      Teuchos::ParameterList&       params,
      DRT::Discretization&          discretization,
      DRT::Element::LocationArray&  la
  );

  //! extract element based or nodal values --> L2-projection case: called within ExtractElementAndNodeValues
  //  return extracted values of phinp
  virtual void ExtractNodalFlux(
      DRT::Element*                 ele,
      Teuchos::ParameterList&       params,
      DRT::Discretization&          discretization,
      DRT::Element::LocationArray&  la,
      const int                     numfluidphases
  );

  //! set internal variables
  virtual void SetInternalVariablesForMatAndRHS();

  //! evaluate material
  virtual void Materials(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               k,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    double&                                 visc,      //!< fluid viscosity
    const int                               iquad=-1   //!< id of current gauss point (default = -1)
  );

  //! material MatMultiPoro
  virtual void MatMultiPoro(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               k,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    double&                                 visc,      //!< fluid viscosity
    const int                               iquad=-1   //!< id of current gauss point (default = -1)
  );

  //! Set advanced reaction terms and derivatives
  virtual void SetAdvancedReactionTerms(
      const int                                 k,           //!< index of current scalar
      const Teuchos::RCP<MAT::MatListReactions> matreaclist, //!< index of current scalar
      const double* gpcoord                                  //!< current Gauss-point coordinates
  );

  //! compute pore pressure
  virtual double ComputePorePressure(
    );

  //! get internal variable manager for multiporo formulation
  Teuchos::RCP<ScaTraEleInternalVariableManagerMultiPoro<my::nsd_,my::nen_> > VarManager()
      {return Teuchos::rcp_static_cast<ScaTraEleInternalVariableManagerMultiPoro<my::nsd_,my::nen_> >(my::scatravarmanager_);};

  //! calculation of convective element matrix in convective form
  //! the only difference to the base class version is, that there is no scaling with the density
  virtual void CalcMatConv(
    Epetra_SerialDenseMatrix&         emat,        //!< element matrix to be filled
    const int                         k,           //!< index of current scalar
    const double                      timefacfac,  //!< domain-integration factor times time-integration factor
    const double                      densnp,      //!< density at time_(n+1)
    const LINALG::Matrix<my::nen_,1>& sgconv       //!< subgrid-scale convective operator
  );

  //! adaption of convective term for rhs
  //! the only difference to the base class version is, that there is no scaling with the density
  virtual void RecomputeConvPhiForRhs(
    const int                          k,          //!< index of current scalar
    const LINALG::Matrix<my::nsd_,1>&  sgvelint,  //!< subgrid-scale velocity at Gauss point
    const double                       densnp,     //!< density at time_(n+1)
    const double                       densn,      //!< density at time_(n)
    const double                       vdiv       //!< velocity divergence
  );

  //! calculation of convective element matrix: add conservative contributions
  //! the only difference to the base class version is, that there is no scaling with the density
  virtual void CalcMatConvAddCons(
    Epetra_SerialDenseMatrix&     emat,        //!< element matrix to be filled
    const int                     k,           //!< index of current scalar
    const double                  timefacfac,  //!< domain-integration factor times time-integration factor
    const double                  vdiv,        //!< velocity divergence
    const double                  densnp       //!< density at time_(n+1)
  );

  //! calculation of mass element matrix (standard shape functions)
  virtual void CalcMatMass(
    Epetra_SerialDenseMatrix &     emat,        //!< element matrix to be filled
    const int &                    k,           //!< index of current scalar
    const double &                 fac,         //!< domain-integration factor
    const double &                 densam       //!< density at time_(n+am)
  );

  //! calculation of convective element matrix (OD term structure coupling)
  //! difference to base class: linearization of mesh motion + shapederivatives pressure gradient have to be included
  virtual void CalcConvODMesh(
    Epetra_SerialDenseMatrix&       emat,
    const int                       k,
    const int                       ndofpernodemesh,
    const double                    fac,
    const double                    rhsfac,
    const double                    densnp,
    const double                    J,
    const LINALG::Matrix<my::nsd_,1>&   gradphi,
    const LINALG::Matrix<my::nsd_,1>&   convelint
    );

  //! calculation of linearized mass (off diagonal/shapederivative term mesh)
  //! difference to base class: linearization of porosity is included
  virtual void CalcLinMassODMesh(
    Epetra_SerialDenseMatrix&                  emat,              //!< element matrix to be filled
    const int                                  k,                 //!< index of current scalar
    const int                                  ndofpernodemesh,   //!< number of dofs per node of ale element
    const double                               rhsfac,            //!< time-integration factor for rhs times domain-integration factor
    const double                               fac,               //!< domain-integration factor
    const double                               densam,            //!< density at time_(n+am)
    const double                               densnp,            //!< density at time_(n+1)
    const double                               phinp,             //!< scalar at time_(n+1)
    const double                               hist,              //!< history of time integartion
    const double                               J,                 //!< determinant of Jacobian det(dx/ds)
    const LINALG::Matrix<1,my::nsd_*my::nen_>& dJ_dmesh           //!< derivative of det(dx/ds) w.r.t. mesh displacement
  );

  //! standard Galerkin transient, old part of rhs and source term (off diagonal/shapederivative term mesh)
  //! difference to base class: linearization of porosity and advanced reaction terms are included
  virtual void CalcHistAndSourceODMesh(
    Epetra_SerialDenseMatrix&                  emat,             //!< element current to be filled
    const int                                  k,                //!< index of current scalar
    const int                                  ndofpernodemesh,  //!< number of dofs per node of ale element
    const double                               fac,              //!< domain-integration factor
    const double                               rhsint,           //!< rhs at Gauss point
    const double                               J,                //!< determinant of Jacobian det(dx/ds)
    const LINALG::Matrix<1,my::nsd_*my::nen_>& dJ_dmesh,         //!< derivative of det(dx/ds) w.r.t. mesh displacement
    const double                               densnp            //!< density
  );

  //! standard Galerkin diffusive term (off diagonal/shapederivative term mesh)
  //! difference to base class: linearization of porosity and effective diffusivity
  virtual void CalcDiffODMesh(
      Epetra_SerialDenseMatrix&                   emat,             //!< element current to be filled
      const int                                   k,                //!< index of current scalar
      const int                                   ndofpernodemesh,  //!< number of dofs per node of ale element
      const double                                diffcoeff,        //!< diffusion coefficient
      const double                                fac,              //!< domain-integration factor
      const double                                rhsfac,           //!< time-integration factor for rhs times domain-integration factor
      const double                                J,                //!< determinant of Jacobian det(dx/ds)
      const LINALG::Matrix<my::nsd_,1>&           gradphi,          //!< scalar gradient at Gauss point
      const LINALG::Matrix<my::nsd_,1>&           convelint,        //!< convective velocity
      const LINALG::Matrix<1,my::nsd_*my::nen_>&  dJ_dmesh          //!< derivative of det(dx/ds) w.r.t. mesh displacement
  );

  //! reactive terms (standard Galerkin) (off diagonal/shapederivative term mesh)
  //! difference to base class: linearization of porosity
  virtual void CalcReactODMesh(
      Epetra_SerialDenseMatrix&                   emat,             //!< element current to be filled
      const int                                   k,                //!< index of current scalar
      const int                                   ndofpernodemesh,  //!< number of dofs per node of ale element
      const double                                rhsfac,           //!< time-integration factor for rhs times domain-integration factor
      const double                                rea_phi,          //!< reactive term
      const double                                J,                //!< determinant of Jacobian det(dx/ds)
      const LINALG::Matrix<1,my::nsd_*my::nen_>&  dJ_dmesh          //!< derivative of det(dx/ds) w.r.t. mesh displacement
  );

  //! generic linearization of term scaled with porosity for coupling with structure (OD terms)
  void PorosityLinearizationStruct(
      Epetra_SerialDenseMatrix&                   emat,             //!< element current to be filled
      const int                                   k,                //!< index of current scalar
      double                                      prefac            //!< generic factor with which the matrix is multiplied
  );

  //! generic linearization of term scaled with porosity for coupling with fluid (OD terms)
  void PorosityLinearizationFluid(
      Epetra_SerialDenseMatrix&                   emat,                         //!< element current to be filled
      const int                                   k,                            //!< index of current scalar
      const int                                   curphase,                     //!< current phase (k belongs to curphase)
      const int                                   numfluidphases,               //!< number of fluid phases
      const int                                   totalnummultiphasedofpernode, //!< number of fluid phases + number of volumefractions
      double                                      prefac                        //!< generic factor with which the matrix is multiplied
  );

  //! generic linearization of term scaled with saturation for coupling with fluid (OD terms)
  void SaturationLinearization(
      Epetra_SerialDenseMatrix&                   emat,                         //!< element current to be filled
      const int                                   k,                            //!< index of current scalar
      const int                                   curphase,                     //!< current phase (k belongs to curphase)
      const int                                   numfluidphases,               //!< number of fluid phases
      const int                                   totalnummultiphasedofpernode, //!< number of fluid phases + number of volumefractions
      double                                      prefac                        //!< generic factor with which the matrix is multiplied
  );

  //! calculation of convective element matrix in convective form (off diagonal term fluid)
  virtual void CalcMatConvODFluid(
    Epetra_SerialDenseMatrix&         emat,             //!< element matrix to be filled
    const int                         k,                //!< index of current scalar
    const int                         ndofpernodefluid, //!< number of dofs per node of fluid element // only a dummy variable
    const double                      timefacfac,       //!< domain-integration factor times time-integration factor
    const double                      densnp,           //!< density at time_(n+1)
    const LINALG::Matrix<my::nsd_,1>& gradphi           //!< scalar gradient
  );

  //! calculation of convective element matrix in convective form -- additional conservative contributions (off diagonal term fluid)
  //! not yet implemented --> dserror
  virtual void CalcMatConvAddConsODFluid(
    Epetra_SerialDenseMatrix&     emat,             //!< element matrix to be filled
    const int                     k,                //!< index of current scalar
    const int                     ndofpernodefluid, //!< number of dofs per node of fluid element
    const double                  timefacfac,       //!< domain-integration factor times time-integration factor
    const double                  densnp,           //!< density at time_(n+1)
    const double                  phinp             //!< scalar at time_(n+1)
  );

  //! calculation of linearized mass (off diagonal terms fluid)
  //! linearization of porosity*saturation*vtrans
  virtual void CalcLinMassODFluid(
    Epetra_SerialDenseMatrix&          emat,              //!< element matrix to be filled
    const int                          k,                 //!< index of current scalar
    const int                          ndofpernodemesh,   //!< number of dofs per node of fluid element // only a dummy variable
    const double                       rhsfac,            //!< time-integration factor for rhs times domain-integration factor
    const double                       fac,               //!< domain-integration factor
    const double                       densam,            //!< density at time_(n+am)
    const double                       densnp,            //!< density at time_(n+1)
    const double                       phinp,             //!< scalar at time_(n+1)
    const double                       hist               //!< history of time integartion
  );

  //! standard Galerkin transient, old part of rhs and source term (off diagonal terms fluid)
  //! linearization of porosity*saturation*vtrans + advanced reaction terms
  virtual void CalcHistAndSourceODFluid(
    Epetra_SerialDenseMatrix&          emat,             //!< element current to be filled
    const int                          k,                //!< index of current scalar
    const int                          ndofpernodemesh,  //!< number of dofs per node of ale element
    const double                       fac,              //!< domain-integration factor
    const double                       rhsint,           //!< rhs at Gauss point
    const double                       densnp            //!< density
  );

  //! standard Galerkin reactive term (off diagonal terms fluid)
  //! linearization of porosity*saturation*vreact
  virtual void CalcReactODFluid(
    Epetra_SerialDenseMatrix&          emat,             //!< element current to be filled
    const int                          k,                //!< index of current scalar
    const int                          ndofpernodemesh,  //!< number of dofs per node of ale element
    const double                       rhsfac,           //!< time-integration factor for rhs times domain-integration factor
    const double                       rea_phi           //!< rhs at Gauss point
  );

  //! standard Galerkin diffusive term (off diagonal terms fluid)
  //! linearization of porosity*saturation*d_eff
  virtual void CalcDiffODFluid(
    Epetra_SerialDenseMatrix&           emat,             //!< element current to be filled
    const int                           k,                //!< index of current scalar
    const int                           ndofpernodemesh,  //!< number of dofs per node of ale element
    const double                        rhsfac,           //!< time-integration factor for rhs times domain-integration factor
    const LINALG::Matrix<my::nsd_,1>&   gradphi           //!< scalar gradient at Gauss point
  );

  // fill coupling vector and add variables to reaction in order to compute reaction values and derivatives
  void FillCouplingVectorAndAddVariables(
      const int                                       k,
      const Teuchos::RCP<MAT::MatListReactions> matreaclist,
      const Teuchos::RCP<ScaTraEleReaManagerAdvReac>  remanager);

  //! Get right hand side including reaction bodyforce term
  virtual void GetRhsInt(
    double&      rhsint,  //!< rhs containing bodyforce at Gauss point
    const double densnp,  //!< density at t_(n+1)
    const int    k        //!< index of current scalar
  );

  //! calculation of reactive element matrix
  virtual void CalcMatReact(
    Epetra_SerialDenseMatrix&          emat,        //!< element matrix to be filled
    const int                          k,           //!< index of current scalar
    const double                       timefacfac,  //!< domain-integration factor times time-integration factor
    const double                       timetaufac,  //!< domain-integration factor times time-integration factor times tau
    const double                       taufac,      //!< domain-integration factor times tau
    const double                       densnp,      //!< density at time_(n+1)
    const LINALG::Matrix<my::nen_,1>&      sgconv,  //!< subgrid-scale convective operator
    const LINALG::Matrix<my::nen_,1>&      diff     //!< laplace term
  );

  protected:

  //! nodal flux values at t_(n+1)
  std::vector<LINALG::Matrix<my::nsd_,my::nen_> > efluxnp_;

  //! a vector containing all quantities, the equation is coupled with
  //! (i.e. pressures, saturations and porosity)
  std::vector<std::pair<std::string,double> > couplingvalues_;

  //! penalty factor to avoid very small "densities"
  const double penalty_ = 1.0;

  //! do we use L2-projection or evaluation at GP
  bool L2_projection_ = false;

};


template<int NSD, int NEN>
class ScaTraEleInternalVariableManagerMultiPoro : public ScaTraEleInternalVariableManager<NSD, NEN>
{
  typedef ScaTraEleInternalVariableManager<NSD, NEN> my;

public:

  ScaTraEleInternalVariableManagerMultiPoro(int numscal):
    ScaTraEleInternalVariableManager<NSD, NEN>(numscal),
    multiconvelint_(0),
    multiconv_(0),
    pressure_(0),
    saturation_(0),
    density_(0),
    solidpressure_(0.0),
    scalartophaseID_(numscal,-1),
    delta_(numscal,0.0),
    min_sat_(numscal,0.0),
    materialset_(false),
    myaction_(SCATRA::Action::calc_mat_and_rhs)
  {return;}

  virtual ~ScaTraEleInternalVariableManagerMultiPoro(){};

  // compute and set internal variables -- no L2-projection but evaluation at GP
   void SetInternalVariablesMultiPoro(
       const LINALG::Matrix<NEN,1>&                 funct,       //! array for shape functions
       const LINALG::Matrix<NSD,NEN>&               derxy,       //! global derivatives of shape functions w.r.t x,y,z
       const LINALG::Matrix<NSD,NEN>&               deriv,       //! global derivatives of shape functions w.r.t r,s,t
       const LINALG::Matrix<NSD,NSD>&               xjm,
       const LINALG::Matrix<NSD,NEN>&               xyze0,
       const std::vector<LINALG::Matrix<NEN,1> >&   ephinp,      //! scalar at t_(n+1) or t_(n+alpha_F)
       const std::vector<LINALG::Matrix<NEN,1> >&   ephin,       //! scalar at t_(n)
       const std::vector<LINALG::Matrix<NEN,1> >&   ehist        //! history vector of transported scalars
       )
   {
     // call base class (scatra) with dummy variable econvelnp
     const LINALG::Matrix<NSD,NEN> econvelnp(true);
     my::SetInternalVariables(funct,derxy,ephinp,ephin,econvelnp,ehist);

     //------------------------get determinant of Jacobian dX / ds
     // transposed jacobian "dX/ds"
     LINALG::Matrix<NSD,NSD> xjm0;
     xjm0.MultiplyNT(deriv,xyze0);

     // inverse of transposed jacobian "ds/dX"
     const double det0= xjm0.Determinant();

     const double det = xjm.Determinant();

     // determinant of deformationgradient det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds) )^-1
     const double JacobianDefGradient = det/det0;

     // clear current gauss point data for safety
     phasemanager_->ClearGPState();
     variablemanager_->EvaluateGPVariables(funct,derxy);
     // access from outside to the phasemanager: scatra-discretization has fluid-dis on dofset 2
     phasemanager_->EvaluateGPState(JacobianDefGradient,*variablemanager_,ndsscatra_porofluid_);

     const int numfluidphases = phasemanager_->NumFluidPhases();
     const int totalnummultiphasedofpernode = phasemanager_->TotalNumDof();

     // resize all phase related vectors
     multiconvelint_.resize(numfluidphases);
     multiconv_.resize(numfluidphases);
     pressure_.resize(numfluidphases);
     saturation_.resize(numfluidphases);
     density_.resize(numfluidphases);
     pressuregrad_.resize(numfluidphases);
     abspressuregrad_.resize(numfluidphases);
     volfrac_.resize(totalnummultiphasedofpernode-numfluidphases);

     const std::vector<LINALG::Matrix<NSD,1> >& fluidgradphi = *(variablemanager_->GradPhinp());

     volfrac_ = phasemanager_->VolFrac();

     for (int k = 0; k < numfluidphases; ++k)
     {
       // current pressure gradient
       pressuregrad_[k].Clear();

       // compute the pressure gradient from the phi gradients
       for (int idof=0; idof<numfluidphases; ++idof)
         pressuregrad_[k].Update(phasemanager_->PressureDeriv(k,idof),fluidgradphi[idof],1.0);

       // compute the absolute value of the pressure gradient from the phi gradients
       abspressuregrad_[k]=0.0;
       for (int i=0; i<NSD; i++)
         abspressuregrad_[k]+=pressuregrad_[k](i)*pressuregrad_[k](i);
       abspressuregrad_[k]=sqrt(abspressuregrad_[k]);

       // diffusion tensor
       LINALG::Matrix<NSD,NSD> difftensor(true);
       phasemanager_->PermeabilityTensor(k,difftensor);
       difftensor.Scale(phasemanager_->RelPermeability(k)/phasemanager_->DynViscosity(k, abspressuregrad_[k],ndsscatra_porofluid_));

       // Insert Darcy's law: porosity*S^\pi*(v^\pi - v_s) = - k/\mu * grad p
       // additional scaling with density since
       // \rho*porosity*S^\pi*(v^\pi - v_s) = - rho*k/\mu * grad p
       multiconvelint_[k].Multiply(-phasemanager_->Density(k),difftensor,pressuregrad_[k]);

       // phase pressure
       pressure_[k] = phasemanager_->Pressure(k);
       // phase saturation
       saturation_[k] = phasemanager_->Saturation(k);
       // phase density
       density_[k] = phasemanager_->Density(k);
       // convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
       multiconv_[k].MultiplyTN(derxy,multiconvelint_[k]);
     }

     // solid pressure
     solidpressure_ = phasemanager_->SolidPressure();

     // check if the IDs make sense
     for (std::vector<int>::iterator it = scalartophaseID_.begin() ; it != scalartophaseID_.end(); ++it)
     {
        if( *it<0 or *it>=numfluidphases )
          dserror("Invalid phase ID %i",*it);
     }

     // overwrite convective term
     // - rho * k/\mu*grad p * grad phi
     for (int k = 0; k < my::numscal_; ++k)
     {
       my::conv_phi_[k] = multiconvelint_[scalartophaseID_[k]].Dot(my::gradphi_[k]);
     }

     return;
   };

   // adapt convective term in case of L2-projection
    void AdaptConvectiveTermForL2(
        const LINALG::Matrix<NEN,1>&                 funct,       //! array for shape functions
        const LINALG::Matrix<NSD,NEN>&               derxy,       //! global derivatives of shape functions w.r.t x,y,z
        const std::vector<LINALG::Matrix<NSD,NEN> >& econvelnp    //! nodal convective velocity values at t_(n+1) or t_(n+alpha_F)
    )
    {

      const int numfluidphases = econvelnp.size();

      // resize all phase related vectors
      multiconvelint_.resize(numfluidphases);
      multiconv_.resize(numfluidphases);

      for (int k = 0; k < numfluidphases; ++k)
      {
        // convective velocity
        multiconvelint_[k].Multiply(phasemanager_->Density(k),econvelnp[k],funct);
        // convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
        multiconv_[k].MultiplyTN(derxy,multiconvelint_[k]);
      }

      // check if the IDs make sense
      for (std::vector<int>::iterator it = scalartophaseID_.begin() ; it != scalartophaseID_.end(); ++it)
      {
         if( *it<0 or *it>=numfluidphases )
           dserror("Invalid phase ID %i",*it);
      }

      // overwrite convective term
      for (int k = 0; k < my::numscal_; ++k)
      {
        my::conv_phi_[k] = multiconvelint_[scalartophaseID_[k]].Dot(my::gradphi_[k]);
      }

      return;
    };

     // Set the fluid-material in the scatra-Varmanager
    void SetFluidPoromultiphaseMaterial(
        DRT::Element*               ele
        )
    {
      // check if we actually have three materials
      if(ele->NumMaterial()<3)
        dserror("no third material available");

      //here we rely that the PoroMultiPhase material has been added as third material
      multiphasemat_ = Teuchos::rcp_dynamic_cast<MAT::FluidPoroMultiPhase>(ele->Material(ndsscatra_porofluid_));
      if(multiphasemat_==Teuchos::null)
        dserror("cast to MAT::FluidPoroMultiPhase failed!");

      materialset_ = true;
    }

  /*========================================================================*/
  //! @name return methods for internal variables
  /*========================================================================*/

  //! return convective velocity
  virtual const LINALG::Matrix<NSD,1>& ConVel(const int k) const
  {
    return multiconvelint_[scalartophaseID_[k]];
  };

  //! return convective part in convective form
  virtual const LINALG::Matrix<NEN,1>& Conv(const int k) const
  {
    return multiconv_[scalartophaseID_[k]];
  };

  //! return pressure associated with scalar k
  double Pressure(const int k) const
  {
    return pressure_[scalartophaseID_[k]];
  };

  //! return pressures
  const std::vector<double>& Pressure() const
  {
    return pressure_;
  };

  //! return saturation associated with scalar k
  double Saturation(const int k) const
  {
    return saturation_[scalartophaseID_[k]];
  };

  //! return saturation associated with scalar k
  const std::vector<double>& Saturation() const
  {
    return saturation_;
  };

  //! return density associated with scalar k
  double Density(const int k) const
  {
    return density_[scalartophaseID_[k]];
  };

  //! return density associated with scalar k
  const std::vector<double>& Density() const
  {
    return density_;
  };

  //! return volfrac vector
  const std::vector<double>& VolFrac() const
  {
    return volfrac_;
  };

  //! return solid pressure
  double SolidPressure() const
  {
    return solidpressure_;
  };

  //! set scalar ID to phase ID mapping
  void SetPhaseID(const int scalarID, const int phaseID)
  {
    scalartophaseID_[scalarID] = phaseID;
  };

  //! get phase ID from scalar ID
  int GetPhaseID(const int scalarID)
  {
    return scalartophaseID_[scalarID];
  };

  //! set delta for evaluation of effective diffusivity
  void SetDelta(const double delta, const int k)
  {
    // set
    delta_[k] = delta;
  };

  //! set minimum saturation
  void SetMinSat(const double minsat, const int k)
  {
    min_sat_[k] = minsat;
  };

  //! set action
  void SetAction(const SCATRA::Action action)
  {
    // set
    myaction_ = action;
  };

  //! get delta
  double GetDelta(const int k)
  {
    return delta_[k];
  };

  //! get minimum saturation
  double GetMinSat(const int k)
  {
    return min_sat_[k];
  };

  //! get action
  SCATRA::Action GetAction()
  {
    // return
    return myaction_;
  };

  //! return pressure gradient of current phase
  const LINALG::Matrix<NSD,1>& PressureGradient(const int curphase) const
  {
    return pressuregrad_[curphase];
  };

  //! return absolute value of pressure gradient of current phase
  double AbsPressureGradient(const int curphase) const
  {
    return abspressuregrad_[curphase];
  };

  //! return fluid material
  Teuchos::RCP<MAT::FluidPoroMultiPhase> MultiphaseMat()
  {
    if(!materialset_)
      dserror("Fluid-Multiphase Material has not yet been set in Variablemanager");

    return multiphasemat_;
  }

  //! Setup phasemanager and variablemanager of fluid
  void SetupPoroFluidManagers(
      DRT::Element*                 ele,
      Teuchos::ParameterList&       params,
      DRT::Discretization&          discretization,
      DRT::Element::LocationArray&  la,
      const int                     numfluidphases,
      const int                     totalnummultiphasedofpernode
      )
  {

    // dummy parameter list
    DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter* para =
        DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter::Instance(discretization.Name());

    // create phase-manager
    phasemanager_ =
        DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerInterface::CreatePhaseManager(
            *para,
            NSD,
            MultiphaseMat()->MaterialType(),
            POROFLUIDMULTIPHASE::Action::get_access_from_scatra,
            totalnummultiphasedofpernode,
            numfluidphases
            );

    // access from outside to the phasemanager: scatra-discretization has fluid-dis on dofset 2
    phasemanager_->Setup(ele,ndsscatra_porofluid_);

    // create variablemanager
    variablemanager_ =
        DRT::ELEMENTS::POROFLUIDMANAGER::VariableManagerInterface<NSD,NEN>::CreateVariableManager(
            *para,
            POROFLUIDMULTIPHASE::Action::get_access_from_scatra,
            totalnummultiphasedofpernode
            );

    return;

  }

  // extract the element and node values of the poro-fluid --> extract them from its variablemanager
  void ExtractElementAndNodeValuesOfPoroFluid(
      DRT::Element*                 ele,
      DRT::Discretization&          discretization,
      DRT::Element::LocationArray&  la,
      LINALG::Matrix<NSD, NEN>& xyze
      )
  {

    // access from outside to the variablemananger: scatra-discretization has fluid-dis on dofset 2
    variablemanager_->ExtractElementAndNodeValues(*ele,discretization,la,xyze,ndsscatra_porofluid_);

    return;

  }

  // get the phasemanager of the fluid
  Teuchos::RCP< DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerInterface > FluidPhaseManager()
  {
    return phasemanager_;
  }

  // get the variablemanager of the fluid
  Teuchos::RCP< DRT::ELEMENTS::POROFLUIDMANAGER::VariableManagerInterface<NSD,NEN> > FluidVarManager()
  {
    return variablemanager_;
  }

private:

  //! convective velocity
  std::vector<LINALG::Matrix<NSD,1> > multiconvelint_;
  //! convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
  std::vector<LINALG::Matrix<NEN,1> > multiconv_;

  //! phase pressure
  std::vector<double> pressure_;
  //! phase saturation
  std::vector<double> saturation_;
  //! phase density
  std::vector<double> density_;
  //! solid pressure
  double solidpressure_;
  //! pressure gradient
  std::vector<LINALG::Matrix<NSD,1> > pressuregrad_;
  //! norm of pressure-gradient
  std::vector<double> abspressuregrad_;
  //! phase pressure
  std::vector<double> volfrac_;

  //! fluid-poro multiphase material
  Teuchos::RCP<MAT::FluidPoroMultiPhase> multiphasemat_;

  //! mapping from scalar ID to phase ID
  std::vector<int > scalartophaseID_;
  //! delta for effective diffusivity
  std::vector<double> delta_;
  //! minimum saturation under which the corresponding mass fraction is assumed to be zero
  std::vector<double> min_sat_;
  //! check if multiphase material has been set
  bool materialset_;

  SCATRA::Action myaction_;

  //! phase manager of the fluid
  Teuchos::RCP< DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerInterface > phasemanager_;

  //! variable manager of the fluid
  Teuchos::RCP< DRT::ELEMENTS::POROFLUIDMANAGER::VariableManagerInterface<NSD,NEN> > variablemanager_;

  //! dofset of fluid field on scatra dis
  // TODO: find a better way to do this
  const int ndsscatra_porofluid_ = 2;

};

}
}

#endif /* SRC_DRT_SCATRA_ELE_SCATRA_ELE_CALC_MULTIPORO_REAC_H_ */
