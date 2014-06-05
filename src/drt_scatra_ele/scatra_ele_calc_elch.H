/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_elch.H

\brief evalution of ScaTra elements for ion-transport equation

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15252
</pre>
*/
/*--------------------------------------------------------------------------*/

#ifndef SCATRA_ELE_CALC_ELCH_H
#define SCATRA_ELE_CALC_ELCH_H

#include "scatra_ele_calc.H"
// due to class ScaTraEleInternalVariableManagerElch()
#include "scatra_ele_parameter_elch.H"

#include "../drt_inpar/inpar_elch.H"

namespace DRT
{
namespace ELEMENTS
{
  // forward declarations
  class ScaTraEleParameterElch;
  class ScaTraEleDiffManagerElch;
  template<int NSD, int NEN>
  class ScaTraEleInternalVariableManagerElch;

template<DRT::Element::DiscretizationType distype>
class ScaTraEleCalcElch : public ScaTraEleCalc<distype>
{
protected:

  /// (private) protected constructor, since we are a Singleton.
  /// this constructor is called from a derived class
  /// -> therefore, it has to be protected instead of private
  ScaTraEleCalcElch(const int numdofpernode,const int numscal);

  typedef ScaTraEleCalc<distype> my;

public:

  virtual ~ScaTraEleCalcElch() {}

  /// called upon destruction
  virtual void Done() = 0;

  /// In this class we do not define a static ScaTraEle...* Instance
  /// since only derived child classes are free to be allocated!!

  //! evaluate the element
  int Evaluate(DRT::ELEMENTS::Transport*     ele,
               Teuchos::ParameterList&       params,
               DRT::Discretization &         discretization,
               const std::vector<int> &      lm,
               Epetra_SerialDenseMatrix&     elemat1_epetra,
               Epetra_SerialDenseMatrix&     elemat2_epetra,
               Epetra_SerialDenseVector&     elevec1_epetra,
               Epetra_SerialDenseVector&     elevec2_epetra,
               Epetra_SerialDenseVector&     elevec3_epetra);

  //! evaluate service routines
  int EvaluateService(DRT::ELEMENTS::Transport*     ele,
                      Teuchos::ParameterList&       params,
                      DRT::Discretization &         discretization,
                      const std::vector<int> &      lm,
                      Epetra_SerialDenseMatrix&     elemat1_epetra,
                      Epetra_SerialDenseMatrix&     elemat2_epetra,
                      Epetra_SerialDenseVector&     elevec1_epetra,
                      Epetra_SerialDenseVector&     elevec2_epetra,
                      Epetra_SerialDenseVector&     elevec3_epetra);

protected:

  /*========================================================================*/
  //! @name general framework
  /*========================================================================*/

  //! Prepare everything what is needed in CallMatAndRhs() to calculate the sysmat and the rhs
  void Sysmat(
    DRT::Element*                         ele,       //!< the element we are dealing with
    Epetra_SerialDenseMatrix&             emat,      //!< element matrix to calculate
    Epetra_SerialDenseVector&             erhs,      //!< element rhs to calculate
    Epetra_SerialDenseVector&             subgrdiff  //!< subgrid-diff.-scaling vector
  );

  //! calculate contributions to matrix and rhs (inside of loop over all scalars)
  virtual void CalMatAndRhs(
    Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> >& vm,
    Epetra_SerialDenseMatrix&                emat,         //!< element matrix to calculate
    Epetra_SerialDenseVector&                erhs,         //!< element rhs to calculate+
    const int                                k,            //!< index of current scalar
    const double                             fac,          //!< domain-integration factor
    const double                             timefacfac,   //!< domain-integration factor times time-integration factor
    const double                             rhsfac,       //!< time-integration factor for rhs times domain-integration factor
    Teuchos::RCP<ScaTraEleDiffManagerElch>&  dme,          //!< diffusion manager
    double&                                  rhsint,       //!< rhs at Gauss point
    const double                             hist          //!< history
  ){return;};

  //! calculate contributions to matrix and rhs (outside of loop over all scalars)
  virtual void CalMatAndRhsOutsideScalarLoop(
    Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> >& vm,
    Epetra_SerialDenseMatrix&                emat,         //!< element matrix to calculate
    Epetra_SerialDenseVector&                erhs,         //!< element rhs to calculate
    const double                             fac,          //!< domain-integration factor
    const double                             timefacfac,   //!< domain-integration factor times time-integration factor
    const double                             rhsfac,       //!< time-integration factor for rhs times domain-integration factor
    Teuchos::RCP<ScaTraEleDiffManagerElch>&  dme           //!< diffusion manager
  ){return;};

  //! Correction for additional flux terms / currents across Dirichlet boundaries
  virtual void CorrectionForFluxAcrossDC(
    DRT::Discretization&        discretization,   //!< discretization
    const std::vector<int>&     lm,               //!< location vector
    Epetra_SerialDenseMatrix&   emat,             //!< element matrix to calculate
    Epetra_SerialDenseVector&   erhs              //!< element rhs to calculate
  ){return;};

  //! Formulation-specific internal variables are set in the variable manager
  virtual void SetFormulationSpecificInternalVariables(
    Teuchos::RCP<ScaTraEleDiffManagerElch>&                                   dme,  //!< rcp to diffusion manager Elch
    Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> >&  vm    //!< rcp to variable manager Elch
  ){return;};

  /*========================================================================*/
  //! @name material and related and related functions
  /*========================================================================*/

  //! material ion
  void MatIon(
    const Teuchos::RCP<const MAT::Material> material,     //!< pointer to current material
    const int                               k,            //!< id of current scalar
    double&                                 densn,        //!< density at t_(n)
    double&                                 densnp,       //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,       //!< density at t_(n+alpha_M)
    Teuchos::RCP<ScaTraEleDiffManagerElch>  diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
    Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
    double&                                 visc,         //!< fluid viscosity
    const int                               iquad         //!< id of current gauss point
  );

  /*========================================================================*/
  //! @name methods for evaluation of individual terms
  /*========================================================================*/

  //! Potential equation ENC
  void CalcMatPotEquENC(
    Epetra_SerialDenseMatrix&                 emat,     //!< element matrix to be filled
    const int                                 k,        //!< index of current scalar
    const double                              fac,      //!< domain-integration factor
    const double                              alphaf,   //!< time factor for ENC
    Teuchos::RCP<ScaTraEleDiffManagerElch>&   dme       //!< diffusion manager
  );

  //! CalcRhs: Potential equation ENC
  void CalcRhsPotEquENC(
    Epetra_SerialDenseVector&                 erhs,     //!< element vector to be filled
    const int                                 k,        //!< index of current scalar
    const double                              fac,      //!< domain-integration factor
    Teuchos::RCP<ScaTraEleDiffManagerElch>&   dme,   //!< diffusion manager
    const double                              conint    //!< concentration at GP
  );

  /*========================================================================*/
  //! @name additional service routines
  /*========================================================================*/

  //! validity check for all elements with respect to
  //! formulation-specific parameter, degree's of freedom, number of scalars, ...
  void CheckElchElementParameter(
    DRT::ELEMENTS::Transport*  ele  //!< the element we are dealing with
    );

  //!  calculate weighted mass flux (no reactive flux so far) -> elch-specific implementation
  virtual void CalculateFlux(
    Epetra_SerialDenseVector&       fluxx,      //!< flux in x-direction to be computed
    Epetra_SerialDenseVector&       fluxy,      //!< flux in y-direction to be computed
    Epetra_SerialDenseVector&       fluxz,      //!< flux in z-direction to be computed
    const INPAR::SCATRA::FluxType   fluxtype,   //!< type fo flux
    const int                       k,          //!< index of current scalar
    const double                    fac,         //!< integration factor
    Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> >& vm,  //!< variable manager
    Teuchos::RCP<ScaTraEleDiffManagerElch>&   dme                                 //!< diffusion manager
  ){return;};

  //!  calculate weighted current flux (no reactive flux so far) -> elch-specific implementation
  virtual void CalculateCurrent(
    Epetra_SerialDenseVector&       fluxx,      //!< flux in x-direction to be computed
    Epetra_SerialDenseVector&       fluxy,      //!< flux in y-direction to be computed
    Epetra_SerialDenseVector&       fluxz,      //!< flux in z-direction to be computed
    const INPAR::SCATRA::FluxType   fluxtype,   //!< type fo flux
    const double                    fac,        //!< integration factor
    Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> >& vm,  //!< variable manager
    Teuchos::RCP<ScaTraEleDiffManagerElch>&                                  dme  //!< diffusion manager
  ){return;};

  //! calculate error of numerical solution with respect to analytical solution
  void CalErrorComparedToAnalytSolution(
    const DRT::Element*                   ele,        //!< the element we are dealing with
    Teuchos::ParameterList&               params,     //!< parameter list
    Epetra_SerialDenseVector&             errors      //!< vector containing L2-error norm
  );

  //! calculate conductivity of electrolyte solution
  void CalculateConductivity(
    const DRT::Element*               ele,        //!< the element we are dealing with
    const enum INPAR::ELCH::EquPot    equpot,     //!< type of closing equation for electric potential
    Epetra_SerialDenseVector&         sigma       //!< conductivity of all single ions + overall electrolyte solution
  );

  //! Solving for electrical potenial field: Evaluate system matrix and rhs
  void CalculateElectricPotentialField(
    const DRT::Element*               ele,        //!< the element we are dealing with
    const enum INPAR::ELCH::EquPot    equpot,     //!< type of closing equation for electric potential
    Epetra_SerialDenseMatrix&         emat,       //!< element matrix
    Epetra_SerialDenseVector&         erhs        //!< element rhs
  );

  //! Solving for intial time derivative: finish preparation of matrix and rhs
  virtual void PrepMatAndRhsInitialTimeDerivative(
    Epetra_SerialDenseMatrix&  elemat1_epetra,  //!< element matrix
    Epetra_SerialDenseVector&  elevec1_epetra   //!< element rhs
  ){return;};

  // Get conductivity from material
  virtual void GetConductivity(
    const enum INPAR::ELCH::EquPot    equpot,     //!< type of closing equation for electric potential
    double&                           sigma_all,  //!< conductivity of electrolyte solution
    Epetra_SerialDenseVector&         sigma       //!< conductivity of all single ions + overall electrolyte solution
  ){return;};

  // Solving for electric potential field : Calculate Matrix and Rhs
  virtual void CalMatAndRhsElectricPotentialField(
    Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> >& vm,
    const enum INPAR::ELCH::EquPot    equpot,     //!< type of closing equation for electric potential
    Epetra_SerialDenseMatrix&         emat,       //!< element matrix
    Epetra_SerialDenseVector&         erhs,       //!< element rhs
    const double                      fac,        //!< integration factor
    Teuchos::RCP<ScaTraEleDiffManagerElch>& dme   //!< diffusion manager
  ){return;};

  /*========================================================================*/
  //! @name flags and enums
  /*========================================================================*/

  //! pointer to parameter class elch
  DRT::ELEMENTS::ScaTraEleParameterElch* elchpara_;

  Teuchos::RCP<ScaTraEleInternalVariableManagerElch <my::nsd_,my::nen_> > varmanager_;

  /*========================================================================*/
  //! @name scalar degrees of freedom and related
  /*========================================================================*/

  //! electric potential
  LINALG::Matrix<my::nen_,1> epotnp_;

};


/// ScaTraEleDiffManagerElch implementation
/*!
  This class keeps all elch-specific transport parameter needed for the evaluation of an element. The
  ScaTraEleDiffManagerElch is derived from the standard ScaTraEleDiffManager.
*/
class ScaTraEleDiffManagerElch : public ScaTraEleDiffManager
{
  public:

    ScaTraEleDiffManagerElch(int numscal):
    ScaTraEleDiffManager(numscal),
    valence_(numscal,0.0)
    {return;}

    ~ScaTraEleDiffManagerElch(){};

    //! Set valence of the single ionic species
    void SetValence(const double valence,const int k)
    {
      valence_[k]=valence;
      return;
    };

    //! Access routine for valence of the single ionic species
    const double GetValence(const int k){return valence_[k];};

    // TODO: ELCH: Has to be moved to scatra_ele_calc_elch_NP.cpp
    unsigned int LengthVector(){return diff_.size();};

    void IncreaseLengthVector(const int k, const int numscal)
    {
      if(diff_.size()== (unsigned) numscal)
      {
        valence_.push_back(0.0);
        diff_.push_back(0.0);
      }

      return;
    };

  protected:

    //! valence of the single ionic species
    std::vector<double> valence_;
};

/// ScaTraEleInternalVariableManagerElch implementation
/*!
  This class manages all internal variables needed for the evaluation of an element.
  The internal variables stored in this class are used by the Nernst-Planck formulation
  as well as the diffusion-conduction formulation.
  All formulation-specific internal variables are stored and managed by a class derived from this class
  (class ScaTraEleInternalVariableManagerElchNP and class ScaTraEleInternalVariableManagerElchDiffCond).
*/
template<int NSD, int NEN>
class ScaTraEleInternalVariableManagerElch
{
  public:

    ScaTraEleInternalVariableManagerElch(int numscal, int nsd,DRT::ELEMENTS::ScaTraEleParameterElch* elchpara):
      numscal_(numscal),
      nsd_(nsd),
      // constant internal variables
      faraday_(INPAR::ELCH::faraday_const),
      epsilon_(INPAR::ELCH::epsilon_const),
      frt_(elchpara->FRT()),
      // internal variables evaluated at the Gauss point
      gradphi_(numscal),
      conint_(numscal),
      conintinv_(numscal),
      gradpot_(true),
      convelint_(true),
      conv_(true),
      convphi_(numscal),
      sgconv_(true)
    {return;}

    virtual ~ScaTraEleInternalVariableManagerElch(){};

    // compute and set internal variables used by Nerst-Planck formulation and the Diffusion-Conduction formulation
    void SetInternalVariablesElch(
      LINALG::Matrix<NEN,1>&                   funct,      //!< array for shape functions
      LINALG::Matrix<NSD,NEN>&                 derxy,      //!< global derivatives of shape functions w.r.t x,y,z
      std::vector<LINALG::Matrix<NEN,1> >&     ephinp,     //!< nodal concentration values at t_(n+1) or t_(n+alpha_F)
      LINALG::Matrix<NEN,1>&                   epotnp,     //!< nodal electric potential values at t_(n+1) or t_(n+alpha_F)
//      LINALG::Matrix<NSD,NEN>&                 evelnp,     //!< nodal velocity values at t_(n+1) or t_(n+alpha_F)
      LINALG::Matrix<NSD,NEN>&                 econvelnp   //!< nodal convective velocity values at t_(n+1) or t_(n+alpha_F)
    )
    {
      // loop over all transported scalars
      // get concentration of transported scalar k at integration point
      // evaluation of all concentrations is necessary at this point since
      // -> homogeneous reactions of scalar k may depend on all concentrations
      // -> concentration depending material parameters for the diffusion-convection formulation
      // -> avoiding of possible errors (concentration was always defined as a vector where only one
      //    entry was filled)
      for (int k = 0; k < numscal_;++k)
      {
        // calculate concentration at GP
        conint_[k] = funct.Dot(ephinp[k]);
        // calculate 1/concentration at GP
        conintinv_[k] = 1/conint_[k];
        // calculatd gradient of concentration at GP
        gradphi_[k].Multiply(derxy,ephinp[k]);
      }
      // calculate gradient of electric potential at GP
      gradpot_.Multiply(derxy,epotnp);
      // calculate velocity at GP
//      velint_.Multiply(evelnp,funct);
      // calculate convective velocity at GP
      convelint_.Multiply(econvelnp,funct);
      // convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
      conv_.MultiplyTN(derxy,convelint_);

      // convective term using current scalar value
      for (int k = 0; k < numscal_;++k)
      {
        //convphi_[k] = convelint_.Dot(gradphi_[k]);
        convphi_[k] = conv_.Dot(ephinp[k]);
      }

      return;
    };

  /*========================================================================*/
  //! @name return constant internal variables
  /*========================================================================*/

    const double Faraday(){return faraday_;};
    const double Epsilon(){return epsilon_;};
    const double FRT(){return frt_;};

  /*========================================================================*/
  //! @name return methods for GP values
  /*========================================================================*/

    //! return gradient of phi for species k at GP
    const LINALG::Matrix<NSD,1> GradPhi(const int k){return gradphi_[k];};
    //! return gradient of phi for all species in a vector at GP
    const std::vector<LINALG::Matrix<NSD,1> > GradPhi(){return gradphi_;};
    //! return concentration of species k at GP
    const double ConInt(const int k){return conint_[k];};
    //! return 1/concentration of species k at GP
    const double ConIntInv(const int k){return conintinv_[k];};
    //! return 1/concentration of all species in a vector at GP
    const std::vector<double> ConIntInv(){return conintinv_;};
    //! return gradient of electric potential at GP
    const LINALG::Matrix<NSD,1> GradPot(){return gradpot_;};
    //! return velocity at GP
//    const LINALG::Matrix<NSD,1> VelInt(){return velint_;};
    //! return convective velocity at GP
    const LINALG::Matrix<NSD,1> ConVelInt(){return convelint_;};
    //! return convective part in convective form
    const LINALG::Matrix<NEN,1> Conv(){return conv_;};
    //! return convective contribution
    const double ConvPhi(const int k ){return convphi_[k];};
    //! return subgrid velocity
    const LINALG::Matrix<NEN,1> SGConv(){return sgconv_;};

  protected:

  /*========================================================================*/
  //! @name number of scalars and number of dof's
  /*========================================================================*/

    //! number of scalars
    int numscal_;
    //! number of space dimensions
    int nsd_;

  /*========================================================================*/
  //! @name constant internal variables
  /*========================================================================*/

    //! Faraday constant
    double faraday_;
    //! dielectric constant
    double epsilon_;
    //! pre-calculation of divisions which are regularly used in following
    //! (a division is much more expensive than a multiplication)
    //! constant parameter F/RT
    double frt_;

  /*========================================================================*/
  //! @name internal variables evaluated at the Gauss point
  /*========================================================================*/

    //! gradient of concentration at GP
    std::vector<LINALG::Matrix<NSD,1> > gradphi_;
    //! concentration at GP
    std::vector<double> conint_;
    //! 1/concentration at GP
    std::vector<double> conintinv_;
    //! gradient of electric potential at integration point
    LINALG::Matrix<NSD,1> gradpot_;
    //! velocity at GP
//    LINALG::Matrix<NSD,1> velint_;
    //! convective velocity at GP
    LINALG::Matrix<NSD,1> convelint_;
    //! convective part in convective form: u_x*N,x+ u_y*N,y
    LINALG::Matrix<NEN,1> conv_;
    //! convective contribution
    std::vector<double> convphi_;
    // subgrid velocity
    LINALG::Matrix<NEN,1> sgconv_;
};

}
}

#endif

