/*!----------------------------------------------------------------------
\file scatra_ele_calc.H

\brief main file containing routines for calculation of scatra element

\level 1

\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251

*----------------------------------------------------------------------*/
#ifndef SCATRA_ELE_CALC_H
#define SCATRA_ELE_CALC_H

#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

#include <Epetra_IntSerialDenseVector.h>

#include "scatra_ele_action.H"
#include "scatra_ele_calc_utils.H"
#include "scatra_ele_interface.H"

#include "../drt_fluid_ele/fluid_ele.H"

namespace FLD
{
  template<DRT::Element::DiscretizationType distype,int numdofpernode, DRT::ELEMENTS::Fluid::EnrichmentType enrtype> class RotationallySymmetricPeriodicBC;
}

namespace DRT
{
namespace ELEMENTS
{
  // forward declarations
  class ScaTraEleParameterStd;
  class ScaTraEleParameterTimInt;
  class ScaTraEleParameterTurbulence;

  class ScaTraEleDiffManager;
  class ScaTraEleReaManager;

  template<int NSD,int NEN> class ScaTraEleInternalVariableManager;

  /// Scatra element implementation
  /*!
    This internal class keeps all the working arrays needed to
    calculate the scatra element. Additionally, the method Sysmat()
    provides a clean and fast element implementation.

    <h3>Purpose</h3>

    The idea is to separate the element maintenance (class Fluid) from the
    mathematical contents (this class). There are different
    implementations of the scatra element, this is just one such
    implementation.

    The scatra element will allocate exactly one object of this class for all
    scatra elements with the same number of nodes in the mesh. This
    allows us to use exactly matching working arrays (and keep them
    around.)

    The code is meant to be as clean as possible. This is the only way
    to keep it fast. The number of working arrays has to be reduced to
    a minimum so that the element fits into the cache. (There might be
    room for improvements.)

    <h3>Usability</h3>

    The calculations are done by the Evaluate() method. There are two
    version. The virtual method that is inherited from ScatraEleInterface
    (and called from Scatra) and the non-virtual one that does the actual
    work. The non-virtual Evaluate() method must be callable without an actual
    Scatra object.
  */

  template<DRT::Element::DiscretizationType distype,int probdim=DRT::UTILS::DisTypeToDim<distype>::dim>
  class ScaTraEleCalc : public ScaTraEleInterface
  {
  protected:

    /// (private) protected constructor, since we are a Singleton.
    /// this constructor is called from a derived class
    /// -> therefore, it has to be protected instead of private
    ScaTraEleCalc(const int numdofpernode,const int numscal, const std::string& disname);

  public:

    virtual ~ScaTraEleCalc() {}

    /// called upon destruction
    virtual void Done() = 0;

    /// In this class we do not define a static ScaTraEle...* Instance
    /// since only derived child classes are free to be allocated!!

    /// Setup element evaluation
    virtual int SetupCalc(
        DRT::Element*               ele,
        DRT::Discretization&        discretization
        );

    /// Evaluate the element
    /*!
      Generic virtual interface function. Called via base pointer.
     */
    virtual int Evaluate(DRT::Element*                 ele,
                         Teuchos::ParameterList&       params,
                         DRT::Discretization&          discretization,
                         DRT::Element::LocationArray&  la,
                         Epetra_SerialDenseMatrix&     elemat1_epetra,
                         Epetra_SerialDenseMatrix&     elemat2_epetra,
                         Epetra_SerialDenseVector&     elevec1_epetra,
                         Epetra_SerialDenseVector&     elevec2_epetra,
                         Epetra_SerialDenseVector&     elevec3_epetra);

    //! evaluate action
    virtual int EvaluateAction(
        DRT::Element*                 ele,
        Teuchos::ParameterList&       params,
        DRT::Discretization&          discretization,
        const SCATRA::Action&         action,
        DRT::Element::LocationArray&  la,
        Epetra_SerialDenseMatrix&     elemat1_epetra,
        Epetra_SerialDenseMatrix&     elemat2_epetra,
        Epetra_SerialDenseVector&     elevec1_epetra,
        Epetra_SerialDenseVector&     elevec2_epetra,
        Epetra_SerialDenseVector&     elevec3_epetra
        );

    //! evaluate action for off-diagonal system matrix block
    virtual int EvaluateActionOD(
        DRT::Element*                 ele,
        Teuchos::ParameterList&       params,
        DRT::Discretization&          discretization,
        const SCATRA::Action&         action,
        DRT::Element::LocationArray&  la,
        Epetra_SerialDenseMatrix&     elemat1_epetra,
        Epetra_SerialDenseMatrix&     elemat2_epetra,
        Epetra_SerialDenseVector&     elevec1_epetra,
        Epetra_SerialDenseVector&     elevec2_epetra,
        Epetra_SerialDenseVector&     elevec3_epetra
        );

    //! evaluate service routine
    int EvaluateService(
        DRT::Element*                 ele,
        Teuchos::ParameterList&       params,
        DRT::Discretization &         discretization,
        DRT::Element::LocationArray&  la,
        Epetra_SerialDenseMatrix&     elemat1_epetra,
        Epetra_SerialDenseMatrix&     elemat2_epetra,
        Epetra_SerialDenseVector&     elevec1_epetra,
        Epetra_SerialDenseVector&     elevec2_epetra,
        Epetra_SerialDenseVector&     elevec3_epetra
        );

    int EvaluateOD(
        DRT::Element*                  ele,
        Teuchos::ParameterList&        params,
        DRT::Discretization &          discretization,
        DRT::Element::LocationArray&   la,
        Epetra_SerialDenseMatrix&      elemat1_epetra,
        Epetra_SerialDenseMatrix&      elemat2_epetra,
        Epetra_SerialDenseVector&      elevec1_epetra,
        Epetra_SerialDenseVector&      elevec2_epetra,
        Epetra_SerialDenseVector&      elevec3_epetra
        );

     /*========================================================================*/
     //! @name static member variables
     /*========================================================================*/

     //! number of element nodes (nomenclature: T. Hughes, The finite element method)
     static const unsigned nen_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

     //! number of space dimensions
     static const unsigned nsd_ = probdim;
     static const unsigned nsd_ele_ = DRT::UTILS::DisTypeToDim<distype>::dim;

     //! element-type specific flag if second derivatives are needed
     static const bool use2ndderiv_ = SCATRA::Use2ndDerivs<distype>::use;

     //! number of components necessary to store second derivatives
     // 1 component  for nsd=1:  (N,xx)
     // 3 components for nsd=2:  (N,xx ; N,yy ; N,xy)
     // 6 components for nsd=3:  (N,xx ; N,yy ; N,zz ; N,xy ; N,xz ; N,yz)
     static const unsigned numderiv2_ = DRT::UTILS::DisTypeToNumDeriv2<distype>::numderiv2;

     typedef ScaTraEleInternalVariableManager<nsd_,nen_> varmanager;

  protected:

    /*========================================================================*/
    //! @name general framework
    /*========================================================================*/

    //! extract element based or nodal values
    //  return extracted values of phinp
    virtual void ExtractElementAndNodeValues(
      DRT::Element*                 ele,
      Teuchos::ParameterList&       params,
      DRT::Discretization&          discretization,
      DRT::Element::LocationArray&  la
    );

    //! extract turbulence approach
    void ExtractTurbulenceApproach(
      DRT::Element*                  ele,
      Teuchos::ParameterList&        params,
      DRT::Discretization&           discretization,
      DRT::Element::LocationArray&   la,
      int&                           nlayer
    );

    //! calculate matrix and rhs. Here the whole thing is hidden.
    virtual void Sysmat(
      DRT::Element*                         ele,       //!< the element we are dealing with
      Epetra_SerialDenseMatrix&             emat,      //!< element matrix to calculate
      Epetra_SerialDenseVector&             erhs,      //!< element rhs to calculate
      Epetra_SerialDenseVector&             subgrdiff  //!< subgrid-diff.-scaling vector
    );

    //! calculate matrix. Here the whole thing is hidden.
    virtual void SysmatODMesh(
      DRT::Element*                         ele,                  //!< the element we are dealing with
      Epetra_SerialDenseMatrix&             emat,                 //!< element matrix to calculate
      const int                             numdofpernode_mesh    //!< number of DOF of mesh displacement field
    );

    //! calculate matrix . Here the whole thing is hidden.
    virtual void SysmatODFluid(
      DRT::Element*                         ele,                  //!< the element we are dealing with
      Epetra_SerialDenseMatrix&             emat,                 //!< element matrix to calculate
      const int                             numdofpernode_fluid   //!< number of DOF of fluid field
    );

    //! get the body force
    virtual void BodyForce(
      const DRT::Element*  ele //!< the element we are dealing with
    );

    //! further node-based source terms not given via Neumann volume condition
    void OtherNodeBasedSourceTerms(
      const std::vector<int>&    lm,                //!< location vector
      DRT::Discretization&       discretization,    //!< discretization
      Teuchos::ParameterList&    params             //!< parameterlist
    );

    //! read element coordinates
    virtual void ReadElementCoordinates(
        const DRT::Element*     ele
    );

    //! evaluate shape functions and their derivatives at element center
    virtual double EvalShapeFuncAndDerivsAtEleCenter();

    //! evaluate shape functions and their derivatives at current integration point
    double EvalShapeFuncAndDerivsAtIntPoint(
      const DRT::UTILS::IntPointsAndWeights<nsd_ele_>& intpoints, //!< integration points
      const int                                        iquad      //!< id of current Gauss point
    );

    //! evaluate shape functions and their derivatives at current integration point
    double EvalShapeFuncAndDerivsInParameterSpace();

    //! set internal variables
    virtual void SetInternalVariablesForMatAndRHS();

    // add nodal displacements to point coordinates
    virtual void UpdateNodeCoordinates()
    {
      xyze_ += edispnp_;
      return;
    };

    //! finite difference check for debugging purposes
    virtual void FDCheck(
      DRT::Element*                              ele,          //!< the element we are dealing with
      Epetra_SerialDenseMatrix&                  emat,         //!< element matrix to calculate
      Epetra_SerialDenseVector&                  erhs,         //!< element rhs to calculate
      Epetra_SerialDenseVector&                  subgrdiff     //!< subgrid-diff.-scaling vector
    );

    /*========================================================================*/
    //! @name routines for additional element evaluations (called from EvaluateAction)
    /*========================================================================*/

    //! calculate mass matrix and rhs for determining initial time derivative
    virtual void CalcInitialTimeDerivative(
      DRT::Element*                 ele,              //!< current element
      Epetra_SerialDenseMatrix&     emat,             //!< element matrix
      Epetra_SerialDenseVector&     erhs,             //!< element residual
      Teuchos::ParameterList&       params,           //!< parameter list
      DRT::Discretization&          discretization,   //!< discretization
      DRT::Element::LocationArray&  la                //!< location array
    );

    //! Correct RHS calculated from CalcRHSLinMass() for the linearized mass term
    virtual void CorrectRHSFromCalcRHSLinMass(
      Epetra_SerialDenseVector&     erhs,
      const int                     k,
      const double                  fac,
      const double                  densnp,
      const double                  phinp
    );

    //!  integrate shape functions over domain
    void IntegrateShapeFunctions(
      const DRT::Element*             ele,       //!< the element we are dealing with
      Epetra_SerialDenseVector&       elevec1,   //!< rhs vector
      const Epetra_IntSerialDenseVector& dofids  //!< index of current scalar
    );

    //!  calculate weighted mass flux (no reactive flux so far)
    virtual void CalculateFlux(
      LINALG::Matrix<3,nen_>&         flux,      //!< flux to be computed
      const DRT::Element*             ele,       //!< the element we are dealing with
      const INPAR::SCATRA::FluxType   fluxtype,  //!< type fo flux
      const int                       k          //!< index of current scalar
    );

    //! calculate domain integral, i.e., surface area or volume of domain element
    void CalcDomainIntegral(
        const DRT::Element*         ele,     //!< the element we are dealing with
        Epetra_SerialDenseVector&   scalar   //!< result vector for scalar integral to be computed
        );

    //! calculate scalar(s) and domain integral
    virtual void CalculateScalars(
      const DRT::Element*             ele,       //!< the element we are dealing with
      Epetra_SerialDenseVector&       scalars,   //!< scalar to be computed
      const bool                      inverting  //!< bool indicating inversion
    );

    void CalculateMomentumAndVolume(
    const DRT::Element*             ele,                //!< the element we are dealing with
    Epetra_SerialDenseVector&       momandvol,          //!< element (volume) momentum vector and minus domain volume
    const double                    interface_thickness //!< interface thickness of smoothing function.
    );

    //! calculate/reconstruct gradient(s) at nodes using L2-projection.
    virtual void CalcGradientAtNodes(
      const DRT::Element*             ele,       //!< the element we are dealing with
      Epetra_SerialDenseMatrix&       elemat1,   //!< element matrix, to be filled
      Epetra_SerialDenseVector&       elevec1,   //!< element RHS vector (x-direction)
      Epetra_SerialDenseVector&       elevec2,   //!< element RHS vector (y-direction)
      Epetra_SerialDenseVector&       elevec3    //!< element RHS vector (z-direction)
      );

    //! calculate/reconstruct gradient(s) at nodes using L2-projection.
    virtual void CalcGradientAtNodes(
      const DRT::Element*             ele,       //!< the element we are dealing with
      Epetra_SerialDenseMatrix&       elemat1,   //!< element matrix, to be filled
      Epetra_SerialDenseMatrix&       elemat2,   //!< element RHS vector (x-direction)
      const INPAR::SCATRA::L2ProjectionSystemType &  systemtype  //!< system type
      );

    //! calculate/reconstruct gradient at element center
    virtual void CalcGradientEleCenter(
      const DRT::Element*             ele,       //!< the element we are dealing with
      Epetra_SerialDenseVector&       elevec1,   //!< gradient at element center
      Epetra_SerialDenseVector&       elevec2    //!< gauss point at element center
      );

    //! calculate/reconstruct gradient(s) at nodes using L2-projection.
    virtual void CalcCurvatureAtNodes(
      const DRT::Element*                             ele,          //!< the element we are dealing with
      Epetra_SerialDenseMatrix&                       elemat1,      //!< element matrix, to be filled
      Epetra_SerialDenseVector&                       elevec1,      //!< element RHS vector (x-direction)
      const std::vector<LINALG::Matrix<nen_,nsd_> >&   egradphinp   //!< the gradient of phi in the x-direction
      );

    //! calculate filtered fields for turbulent Prandtl number
    void CalcBoxFilter(
        DRT::Element*                  ele,
        Teuchos::ParameterList&        params,
        DRT::Discretization&           discretization,
        DRT::Element::LocationArray&   la
        );

    //! calculate error of numerical solution with respect to analytical solution
    void CalErrorComparedToAnalytSolution(
      const DRT::Element*                   ele,        //!< the element we are dealing with
      Teuchos::ParameterList&               params,     //!< parameter list
      Epetra_SerialDenseVector&             errors      //!< vector containing L2-error norm
    );

    //! calculate matrix and rhs. Here the whole thing is hidden.
    virtual void CalcHeteroReacMatAndRHS(
      DRT::Element*                         ele,       //!< the element we are dealing with
      Epetra_SerialDenseMatrix&             emat,      //!< element matrix to calculate
      Epetra_SerialDenseVector&             erhs      //!< element rhs to calculate
    );


    /*========================================================================*/
    //! @name material and related functions
    /*========================================================================*/

    //! get the material parameters
    virtual void GetMaterialParams(
      const DRT::Element*  ele,       //!< the element we are dealing with
      std::vector<double>& densn,     //!< density at t_(n)
      std::vector<double>& densnp,    //!< density at t_(n+1) or t_(n+alpha_F)
      std::vector<double>& densam,    //!< density at t_(n+alpha_M)
      double&              visc,       //!< fluid viscosity
      const int            iquad=-1    //!< id of current gauss point (default = -1)
    );

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

    //! material ScaTra
    virtual void MatScaTra(
      const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
      const int                               k,        //!< id of current scalar
      double&                                 densn,    //!< density at t_(n)
      double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
      double&                                 densam,   //!< density at t_(n+alpha_M)
      double&                                 visc,     //!< fluid viscosity
      const int                               iquad=-1  //!< id of current gauss point (default = -1)
    );

    //! evaluate multi-scale scalar transport material
    virtual void MatScaTraMultiScale(
        const Teuchos::RCP<const MAT::Material>   material,   //!< multi-scale scalar transport material
        double&                                   densn,      //!< density at time t_(n)
        double&                                   densnp,     //!< density at time t_(n+1) or t_(n+alpha_f)
        double&                                   densam      //!< density at time t_(n+alpha_m)
        ) const;

    //! evaluate electrode material
    void MatElectrode(
        const Teuchos::RCP<const MAT::Material>   material   //!< electrode material
        );

    /*========================================================================*/
    //! @name stabilization and related functions
    /*========================================================================*/

    //! calculate stabilization parameter
    virtual void CalcTau(
      double&                         tau,         //!< the stabilisation parameters (one per transported scalar)
      const double                    diffus,      //!< diffusivity or viscosity
      const double                    reacoeff,    //!< reaction coefficient
      const double                    densnp,      //!< density at t_(n+1)
      const LINALG::Matrix<nsd_,1>&   convelint,   //!< convective velocity at integration point
      const double                    vol          //!< element volume
    );

    //! calculate characteristic element length
    virtual double CalcCharEleLength(
      const double                  vol,        //!< element volume
      const double                  vel_norm,   //!< norm of velocity
      const LINALG::Matrix<nsd_,1>& convelint   //!< convective velocity at integration point
    );

    //! calculation of tau according to Taylor, Hughes and Zarins
    virtual void CalcTauTaylorHughesZarins(
      double&                      tau,       //!< the stabilisation parameters (one per transported scalar)
      const double                 diffus,    //!< diffusivity or viscosity
      const double                 reacoeff,  //!< reaction coefficient
      const double                 densnp,    //!< density at t_(n+1)
      const LINALG::Matrix<nsd_,1>& convelint  //!< convective velocity at integration point
    );

    //! calculation of tau according to Franca and Valentin
    virtual void CalcTauFrancaValentin(
      double&                      tau,       //!< the stabilisation parameters (one per transported scalar)
      const double                 diffus,    //!< diffusivity or viscosity
      const double                 reacoeff,  //!< reaction coefficient
      const double                 densnp,    //!< density at t_(n+1)
      const LINALG::Matrix<nsd_,1>& convelint, //!< convective velocity at integration point
      const double                 vol        //!< element volume
    );

    //! calculation of tau according to Franca, Shakib and Codina
    virtual void CalcTauFrancaShakibCodina(
      double&                      tau,       //!< the stabilisation parameters (one per transported scalar)
      const double                 diffus,    //!< diffusivity or viscosity
      const double                 reacoeff,  //!< reaction coefficient
      const double                 densnp,    //!< density at t_(n+1)
      const LINALG::Matrix<nsd_,1>& convelint, //!< convective velocity at integration point
      const double                 vol        //!< element volume
    );

    //! calculation of tau according to Codina
    virtual void CalcTauCodina(
      double&                      tau,       //!< the stabilisation parameters (one per transported scalar)
      const double                 diffus,    //!< diffusivity or viscosity
      const double                 reacoeff,  //!< reaction coefficient
      const double                 densnp,    //!< density at t_(n+1)
      const LINALG::Matrix<nsd_,1>& convelint, //!< convective velocity at integration point
      const double                 vol        //!< element volume
    );

    //! calculation of tau according to Franca, Madureira and Valentin
    virtual void CalcTauFrancaMadureiraValentin(
      double&                      tau,       //!< the stabilisation parameters (one per transported scalar)
      const double                 diffus,    //!< diffusivity or viscosity
      const double                 reacoeff,  //!< reaction coefficient
      const double                 densnp,    //!< density at t_(n+1)
      const double                 vol        //!< element volume
    );

    //! exact calculation of tau for 1D
    virtual void CalcTau1DExact(
      double&                      tau,       //!< the stabilisation parameters (one per transported scalar)
      const double                 diffus,    //!< diffusivity or viscosity
      const double                 reacoeff,  //!< reaction coefficient
      const double                 densnp,    //!< density at t_(n+1)
      const LINALG::Matrix<nsd_,1>& convelint, //!< convective velocity at integration point
      const double                 vol        //!< element volume
    );

    //! calculation of strong residual for stabilization
    virtual void CalcStrongResidual(
      const int                   k,          //!< index of current scalar
      double&                     scatrares,  //!< residual of convection-diffusion-reaction eq
      const double                densam,     //!< density at t_(n+am)
      const double                densnp,     //!< density at t_(n+1)
      const double                rea_phi,    //!< reactive contribution
      const double                rhsint,     //!< rhs at integration point
      const double                tau         //!< the stabilisation parameter
    );

    //! calculate artificial diffusivity
    void CalcArtificialDiff(
      const double                        vol,        //!< element volume
      const int                           k,          //!< id of current scalar
      const double                        densnp,     //!< density at t_(n+1)
      const LINALG::Matrix<nsd_,1>&       convelint,  //!< convective velocity at integration point
      const LINALG::Matrix<nsd_,1>&       gradphi,    //!< scalar gradient
      const double                        conv_phi,   //!< convective contribution
      const double                        scatrares,  //!< residual of convection-diffusion-reaction eq
      const double                        tau         //!< the stabilisation parameter
    );

    //! calculate subgrid-scale velocity
    void CalcSubgrVelocity(
      const DRT::Element*           ele,        //!< the element we are dealing with
      LINALG::Matrix<nsd_,1>&       sgvelint,   //!< subgrid velocity at integration point
      const double                  densam,     //!< density at t_(n+am)
      const double                  densnp,     //!< density at t_(n+1)
      const double                  visc,       //!< fluid viscosity
      const LINALG::Matrix<nsd_,1>& convelint,  //!< convective velocity at integration point
      const double                  tau         //!< the stabilisation parameter
    );

    /*========================================================================*/
    //! @name turbulence and related functions
    /*========================================================================*/

    //! calculate mean turbulent Prandtl number
    void GetMeanPrtOfHomogenousDirection(
      Teuchos::ParameterList&    turbmodelparams,  //!< turbulence parameter list
      int&                       nlayer            //!< layer of homogeneous plane
    );

    //! output of model parameters
    void StoreModelParametersForOutput(
      const DRT::Element*                   ele,             //!< the element we are dealing with
      const bool                            isowned,         //!< owner
      Teuchos::ParameterList&               turbulencelist,  //!< turbulence parameter list
      const int                             nlayer           //!< layer of homogeneous plane
    );

    //! calculate all-scale art. subgrid diffusivity
    void CalcSubgrDiff(
      double&                               visc,      //!< fluid viscosity
      const double                          vol,       //!< element volume
      const int                             k,         //!< index of current scalar
      const double                          densnp     //!< density at t_(n+1)
   );

    //!  calculate fine-scale art. subgrid diffusivity
    void CalcFineScaleSubgrDiff(
      double&                               sgdiff,    //!< subgrid-scale diffusion
      Epetra_SerialDenseVector&             subgrdiff, //!< subgrid-scale diffusion vector
      DRT::Element*                         ele,       //!< the element we are dealing with
      const double                          vol,       //!< element volume
      const int                             k,         //!< index of current scalar
      const double                          densnp,    //!< density at t_(n+1)
      const double                          diffus,    //!< diffusion
      const LINALG::Matrix<nsd_,1>          convelint  //!< convective velocity
    );

    //! calculation of coefficients B and D for multifractal subgrid-scales
    void CalcBAndDForMultifracSubgridScales(
      LINALG::Matrix<nsd_,1>&                     B_mfs,      //!< coefficient for fine-scale velocity (will be filled)
      double &                                    D_mfs,      //!< coefficient for fine-scale scalar (will be filled)
      const double                                vol,        //!< volume of element
      const int                                   k,          //!< index of current scalar
      const double                                densnp,     //!< density at t_(n+1)
      const double                                diffus,     //!< diffusivity
      const double                                visc,       //!< fluid vicosity
      const LINALG::Matrix<nsd_,1>                convelint,  //!< convective velocity
      const LINALG::Matrix<nsd_,1>                fsvelint    //!< fine-scale velocity
    );

    // calculate reference length for multifractal subgrid-scales
    double CalcRefLength(
      const double                       vol,        //!< volume of element
      const LINALG::Matrix<nsd_,1>       convelint  //!< convective velocity
    );

    //! calculate filtered quantities for dynamic Smagorinsky model
    void scatra_apply_box_filter(
      double&                                            dens_hat,
      double&                                            temp_hat,
      double&                                            dens_temp_hat,
      double&                                            phi2_hat,
      double&                                            phiexpression_hat,
      Teuchos::RCP<std::vector<double> >                 vel_hat,
      Teuchos::RCP<std::vector<double> >                 densvel_hat,
      Teuchos::RCP<std::vector<double> >                 densveltemp_hat,
      Teuchos::RCP<std::vector<double> >                 densstraintemp_hat,
      Teuchos::RCP<std::vector<double> >                 phi_hat,
      Teuchos::RCP<std::vector<std::vector<double> > >   alphaijsc_hat,
      double&                                            volume,
      const DRT::Element*                                ele,
      Teuchos::ParameterList&                            params
    );

    //! get density at integration point
    virtual double GetDensity(
        const DRT::Element*                 ele,
        Teuchos::RCP<const MAT::Material>   material,
        Teuchos::ParameterList&             params,
        const double                        tempnp
        );

    //! calculate viscous part of subgrid-scale velocity
    virtual void CalcSubgrVelocityVisc(
        LINALG::Matrix<nsd_,1>&   epsilonvel
        );

    //! calculate turbulent Prandtl number for dynamic Smagorinsky model
    void scatra_calc_smag_const_LkMk_and_MkMk(
      Teuchos::RCP<Epetra_MultiVector>&  col_filtered_vel,
      Teuchos::RCP<Epetra_MultiVector>&  col_filtered_dens_vel,
      Teuchos::RCP<Epetra_MultiVector>&  col_filtered_dens_vel_temp,
      Teuchos::RCP<Epetra_MultiVector>&  col_filtered_dens_rateofstrain_temp,
      Teuchos::RCP<Epetra_Vector>&       col_filtered_temp,
      Teuchos::RCP<Epetra_Vector>&       col_filtered_dens,
      Teuchos::RCP<Epetra_Vector>&       col_filtered_dens_temp,
      double&                   LkMk,
      double&                   MkMk,
      double&                   xcenter,
      double&                   ycenter,
      double&                   zcenter,
      const DRT::Element*       ele
    );

    void scatra_calc_vreman_dt(
      Teuchos::RCP<Epetra_MultiVector>& col_filtered_phi,
      Teuchos::RCP<Epetra_Vector>& col_filtered_phi2,
      Teuchos::RCP<Epetra_Vector>&   col_filtered_phiexpression,
      Teuchos::RCP<Epetra_MultiVector>& col_filtered_alphaijsc,
      double& dt_numerator,
      double& dt_denominator,
      const DRT::Element* ele
    );

    //! calculate normalized subgrid-diffusivity matrix
    virtual void CalcSubgrDiffMatrix(
      const DRT::Element*        ele,    //!< the element we are dealing with
      Epetra_SerialDenseMatrix&  emat   //!< element matrix to calculate
    );

    //! calculate dissipation introduced by stabilization and turbulence models
    void CalcDissipation(
      Teuchos::ParameterList&               params,          //!< parameter list
      DRT::Element*                         ele,             //!< pointer to element
      DRT::Discretization&                  discretization,  //!< scatra discretization
      DRT::Element::LocationArray&          la               //!< location array
    );

    /*========================================================================*/
    //! @name methods for evaluation of individual terms
    /*========================================================================*/

    //! calculate the Laplacian in strong form for all shape functions
    void GetLaplacianStrongForm(
      LINALG::Matrix<nen_,1>&                diff  //!< laplace term to be computed
    );

    //! calculate divergence of vector field (e.g., velocity)
    void GetDivergence(
      double&                          vdiv,
      const LINALG::Matrix<nsd_,nen_>& evel
    );

    //! calculate rate of strain of (fine-scale) velocity
    inline double GetStrainRate(
      const LINALG::Matrix<nsd_,nen_>& evel
    )
    {
      // evel is tranferred here since the evaluation of the strain rate can be performed
      // for various velocities such as velint_, fsvel_, ...

      double rateofstrain=0;

      // get velocity derivatives at integration point
      //
      //              +-----  dN (x)
      //   dvel (x)    \        k
      //   -------- =   +     ------ * vel
      //      dx       /        dx        k
      //        j     +-----      j
      //              node k
      //
      // j : direction of derivative x/y/z
      //
      LINALG::Matrix<nsd_,nsd_> velderxy;
      velderxy.MultiplyNT(evel,derxy_);

      // compute (resolved) rate of strain
      //
      //          +-                                 -+ 1
      //          |          /   \           /   \    | -
      //          | 2 * eps | vel |   * eps | vel |   | 2
      //          |          \   / ij        \   / ij |
      //          +-                                 -+
      //
      LINALG::Matrix<nsd_,nsd_> two_epsilon;
      for (unsigned rr=0;rr<nsd_;++rr)
      {
        for (unsigned mm=0;mm<nsd_;++mm)
        {
          two_epsilon(rr,mm) = velderxy(rr,mm) + velderxy(mm,rr);
        }
      }

      for (unsigned rr=0;rr<nsd_;rr++)
      {
        for (unsigned mm=0;mm<nsd_;mm++)
        {
          rateofstrain += two_epsilon(rr,mm)*two_epsilon(mm,rr);
        }
      }

      return(sqrt(rateofstrain/2.0));
    };

    //! calculate the Laplacian (weak form)
    void GetLaplacianWeakForm(
      double&                          val,   //!< ?
      const int                        vi,    //!< ?
      const int                        ui     //!< ?
    )
    {
      val = 0.0;
      for (unsigned j = 0; j<nsd_; j++)
      {
        val += derxy_(j, vi)*derxy_(j, ui);
      }
      return;
    };

    //! calculate the Laplacian (weak form)
    void GetLaplacianWeakForm(
      double&                          val,   //!< ?
      const LINALG::Matrix<nsd_,nsd_>& diffus3, //!< ?
      const int                        vi,    //!< ?
      const int                        ui     //!< ?
    )
    {
      val = 0.0;
      for (unsigned j = 0; j<nsd_; j++)
      {
        for (unsigned i = 0; i<nsd_; i++)
        {
          val += derxy_(j, vi)*diffus3(j,i)*derxy_(i, ui);
        }
      }
      return;
    };

    //! calculate the Laplacian (weak form)
    void GetLaplacianWeakFormRHS(
      double& val,                              //!< ?
      const LINALG::Matrix<nsd_,1>&    gradphi, //!< ?
      const int vi                              //!< ?
    )
    {
      val = 0.0;
      for (unsigned j = 0; j<nsd_; j++)
      {
        val += derxy_(j,vi)*gradphi(j);
      }
      return;
    };

    //! compute rhs containing bodyforce
    virtual void GetRhsInt(
      double&      rhsint,  //!< rhs containing bodyforce at integration point
      const double densnp,  //!< density at t_(n+1)
      const int    k        //!< index of current scalar
    );

    //! calculation of convective element matrix in convective form
    virtual void CalcMatConv(
      Epetra_SerialDenseMatrix&     emat,        //!< element matrix to be filled
      const int                     k,           //!< index of current scalar
      const double                  timefacfac,  //!< domain-integration factor times time-integration factor
      const double                  densnp,      //!< density at time_(n+1)
      const LINALG::Matrix<nen_,1>& sgconv       //!< subgrid-scale convective operator
    );

    //! calculation of convective element matrix: add conservative contributions
    virtual void CalcMatConvAddCons(
      Epetra_SerialDenseMatrix&     emat,        //!< element matrix to be filled
      const int                     k,           //!< index of current scalar
      const double                  timefacfac,  //!< domain-integration factor times time-integration factor
      const double                  vdiv,        //!< velocity divergence
      const double                  densnp       //!< density at time_(n+1)
    );

    //! calculation of diffusive element matrix
    virtual void CalcMatDiff(
      Epetra_SerialDenseMatrix&     emat,        //!< element matrix to be filled
      const int                     k,           //!< index of current scalar
      const double                  timefacfac   //!< domain-integration factor times time-integration factor
    );

    //!calculation of stabilization element matrix
    void CalcMatTransConvDiffStab(
      Epetra_SerialDenseMatrix&     emat,        //!< element matrix to be filled
      const int                     k,           //!< index of current scalar
      const double                  timetaufac,  //!< domain-integration factor times time-integration factor times tau
      const double                  densnp,      //!< density  at time_(n+1)
      const LINALG::Matrix<nen_,1>& sgconv,      //!< subgrid-scale convective operator
      const LINALG::Matrix<nen_,1>& diff         //!< laplace term
    );

    //! calculation of mass element matrix (standard shape functions)
    virtual void CalcMatMass(
      Epetra_SerialDenseMatrix &     emat,        //!< element matrix to be filled
      const int &                    k,           //!< index of current scalar
      const double &                 fac,         //!< domain-integration factor
      const double &                 densam       //!< density at time_(n+am)
    );

    //! calculation of mass element matrix (dual shape functions)
    void CalcMatMass(
      Epetra_SerialDenseMatrix &     emat,        //!< element matrix to be filled
      const int &                    k,           //!< index of current scalar
      const double &                 fac,         //!< domain-integration factor
      const double &                 densam,      //!< density at time_(n+am)
      const LINALG::Matrix<nen_,1> & sfunct,      //!< solution function values
      const LINALG::Matrix<nen_,1> & tfunct       //!< test function values
      ) const;

    //! calculation of stabilization mass element matrix
    void CalcMatMassStab(
      Epetra_SerialDenseMatrix&     emat,        //!< element matrix to be filled
      const int                     k,           //!< index of current scalar
      const double                  taufac,      //!< domain-integration factor times tau
      const double                  densam,      //!< density at time_(n+am)
      const double                  densnp,      //!< density at time_(n+1)
      const LINALG::Matrix<nen_,1>& sgconv,      //!< subgrid-scale convective operator
      const LINALG::Matrix<nen_,1>& diff         //!< laplace term
    );

    //! calculation of reactive element matrix
    virtual void CalcMatReact(
      Epetra_SerialDenseMatrix&          emat,        //!< element matrix to be filled
      const int                          k,           //!< index of current scalar
      const double                       timefacfac,  //!< domain-integration factor times time-integration factor
      const double                       timetaufac,  //!< domain-integration factor times time-integration factor times tau
      const double                       taufac,      //!< domain-integration factor times tau
      const double                       densnp,      //!< density at time_(n+1)
      const LINALG::Matrix<nen_,1>&      sgconv,      //!< subgrid-scale convective operator
      const LINALG::Matrix<nen_,1>&      diff         //!< laplace term
    );

    virtual void CalcMatChemo(
      Epetra_SerialDenseMatrix&           emat,
      const int                           k,
      const double                        timefacfac,
      const double                        timetaufac,
      const double                        densnp,
      const double                        scatrares,
      const LINALG::Matrix<nen_,1>&       sgconv,
      const LINALG::Matrix<nen_,1>&       diff
    ){return;};

    //! calculation of linearized mass rhs vector
    void CalcRHSLinMass(
      Epetra_SerialDenseVector&     erhs,        //!< element vector to be filled
      const int                     k,           //!< index of current scalar
      const double                  rhsfac,      //!< time-integration factor for rhs times domain-integration factor
      const double                  fac,         //!< domain-integration factor
      const double                  densam,      //!< density at time_(n+am)
      const double                  densnp       //!< density at time_(n+1)
    );

    //! adaption of rhs with respect to time integration
    void ComputeRhsInt(
      double&                       rhsint,      //!< rhs at Gauss point
      const double                  densam,      //!< density at time_(n+am)
      const double                  densnp,      //!< density at time_(n+1)
      const double                  hist         //!< history of time integartion
    );

    //! adaption of residual with respect to time integration
    void RecomputeScatraResForRhs(
      double&                       scatrares,  //!< residual of convection-diffusion-reaction eq
      const int                     k,          //!< index of current scalar
      const LINALG::Matrix<nen_,1>&  diff,      //!< laplace term
      const double                  densn,      //!< density at time_(n)
      const double                  densnp,     //!< density at time_(n+1)
      double&                       rea_phi,    //!< reactive contribution
      const double                  rhsint      //!< rhs at Gauss point
    );

    //! adaption of convective term for rhs
    virtual void RecomputeConvPhiForRhs(
      const int                     k,          //!< index of current scalar
      const LINALG::Matrix<nsd_,1>&  sgvelint,  //!< subgrid-scale velocity at Gauss point
      const double                  densnp,     //!< density at time_(n+1)
      const double                  densn,      //!< density at time_(n)
      const double                  vdiv       //!< velocity divergence
    );

    //! standard Galerkin transient, old part of rhs and source term
    void CalcRHSHistAndSource(
      Epetra_SerialDenseVector&     erhs,        //!< element vector to be filled
      const int                     k,           //!< index of current scalar
      const double                  fac,         //!< domain-integration factor
      const double                  rhsint       //!< rhs at Gauss point
    );

    //! standard Galerkin convective term on right hand side
    void CalcRHSConv(
      Epetra_SerialDenseVector&     erhs,        //!< element vector to be filled
      const int                     k,           //!< index of current scalar
      const double                  rhsfac       //!< time-integration factor for rhs times domain-integration factor
    );

    //! standard Galerkin diffusive term on right hand side
    virtual void CalcRHSDiff(
      Epetra_SerialDenseVector&     erhs,        //!< element vector to be filled
      const int                     k,           //!< index of current scalar
      const double                  rhsfac       //!< time-integration factor for rhs times domain-integration factor
    );

    //! transient, convective and diffusive stabilization terms on right hand side
    void CalcRHSTransConvDiffStab(
      Epetra_SerialDenseVector&     erhs,        //!< element vector to be filled
      const int                     k,           //!< index of current scalar
      const double                  rhstaufac,   //!< time-integration factor for rhs times domain-integration factor times tau
      const double                  densnp,     //!< density at time_(n+1)
      const double                  scatrares,   //!< residual of convection-diffusion-reaction eq
      const LINALG::Matrix<nen_,1>& sgconv,      //!< subgrid-scale convective operator
      const LINALG::Matrix<nen_,1>& diff         //!< laplace term
    );

    //! reactive terms (standard Galerkin and stabilization) on rhs
    void CalcRHSReact(
      Epetra_SerialDenseVector&     erhs,        //!< element vector to be filled
      const int                     k,           //!< index of current scalar
      const double                  rhsfac,      //!< time-integration factor for rhs times domain-integration factor
      const double                  rhstaufac,   //!< time-integration factor for rhs times domain-integration factor times tau
      const double                  rea_phi,     //!< reactive term
      const double                  densnp,      //!< density at time_(n+1)
      const double                  scatrares    //!< residual of convection-diffusion-reaction eq
    );

    virtual void CalcRHSChemo(
        Epetra_SerialDenseVector&     erhs,
        const int                     k,
        const double                  rhsfac,
        const double                  rhstaufac,
        const double                  scatrares,
        const double                  densnp
    ){return;};

    //! fine-scale subgrid-diffusivity term on right hand side
    void CalcRHSFSSGD(
      Epetra_SerialDenseVector&     erhs,        //!< element vector to be filled
      const int                     k,           //!< index of current scalar
      const double                  rhsfac,      //!< time-integration factor for rhs times domain-integration factor
      const double                  sgdiff,      //!< subgrid-scale diffusivity
      const LINALG::Matrix<nsd_,1>  fsgradphi    //!< gardient of fine-scale velocity
    );

    //! multifractal subgrid-scale modeling on right hand side only rasthofer 11/13  |
    void CalcRHSMFS(
      Epetra_SerialDenseVector&     erhs,        //!< element vector to be filled
      const int                     k,           //!< index of current scalar
      const double                  rhsfac,      //!< time-integration factor for rhs times domain-integration factor
      const double                  densnp,      //!< density at time_(n+1)
      const LINALG::Matrix<nsd_,1>  mfsggradphi, //!< gradient of multifractal fine-scale scalar
      const LINALG::Matrix<nsd_,1>  mfsgvelint,  //!< multifractal fine-scale velocity
      const double                  mfssgphi,    //!< multifractal fine-scale scalar
      const double                  mfsvdiv      //!< divergence of fine-scale velocity
    );

    //! macro-scale matrix and vector contributions arising from macro-micro coupling in multi-scale simulations
    virtual void CalcMatAndRhsMultiScale(
        const DRT::Element* const   ele,          //!< element
        Epetra_SerialDenseMatrix&   emat,         //!< element matrix
        Epetra_SerialDenseVector&   erhs,         //!< element right-hand side vector
        const int                   k,            //!< species index
        const int                   iquad,        //!< Gauss point index
        const double                timefacfac,   //!< domain integration factor times time integration factor
        const double                rhsfac        //!< domain integration factor times time integration factor for right-hand side vector
    );

    //! calculation of convective element matrix in convective form (off diagonal term fluid)
    virtual void CalcMatConvODFluid(
      Epetra_SerialDenseMatrix&     emat,             //!< element matrix to be filled
      const int                     k,                //!< index of current scalar
      const int                     ndofpernodefluid, //!< number of dofs per node of fluid element
      const double                  timefacfac,       //!< domain-integration factor times time-integration factor
      const double                  densnp,           //!< density at time_(n+1)
      const LINALG::Matrix<nsd_,1>& gradphi           //!< scalar gradient
    );

    //! calculation of convective element matrix: add conservative contributions (off diagonal term fluid)
    virtual void CalcMatConvAddConsODFluid(
      Epetra_SerialDenseMatrix&     emat,             //!< element matrix to be filled
      const int                     k,                //!< index of current scalar
      const int                     ndofpernodefluid, //!< number of dofs per node of fluid element
      const double                  timefacfac,       //!< domain-integration factor times time-integration factor
      const double                  densnp,           //!< density at time_(n+1)
      const double                  phinp             //!< scalar at time_(n+1)
    );

    //! calculation of linearized mass (off diagonal/shapederivative term mesh)
    virtual void CalcLinMassODMesh(
      Epetra_SerialDenseMatrix&          emat,              //!< element matrix to be filled
      const int                          k,                 //!< index of current scalar
      const int                          ndofpernodemesh,   //!< number of dofs per node of ale element
      const double                       rhsfac,            //!< time-integration factor for rhs times domain-integration factor
      const double                       fac,               //!< domain-integration factor
      const double                       densam,            //!< density at time_(n+am)
      const double                       densnp,            //!< density at time_(n+1)
      const double                       phinp,             //!< scalar at time_(n+1)
      const double                       hist,              //!< history of time integartion
      const double                       J,                 //!< determinant of Jacobian det(dx/ds)
      const LINALG::Matrix<1,nsd_*nen_>& dJ_dmesh           //!< derivative of det(dx/ds) w.r.t. mesh displacement
    );

    //! standard Galerkin transient, old part of rhs and source term (off diagonal/shapederivative term mesh)
    virtual void CalcHistAndSourceODMesh(
      Epetra_SerialDenseMatrix&          emat,             //!< element current to be filled
      const int                          k,                //!< index of current scalar
      const int                          ndofpernodemesh,  //!< number of dofs per node of ale element
      const double                       fac,              //!< domain-integration factor
      const double                       rhsint,           //!< rhs at Gauss point
      const double                       J,                //!< determinant of Jacobian det(dx/ds)
      const LINALG::Matrix<1,nsd_*nen_>& dJ_dmesh,         //!< derivative of det(dx/ds) w.r.t. mesh displacement
      const double                       densnp
    );

    //! standard Galerkin convective term (off diagonal)
    virtual void CalcConvODMesh(
        Epetra_SerialDenseMatrix&       emat,             //!< element current to be filled
        const int                       k,                //!< index of current scalar
        const int                       ndofpernodemesh,  //!< number of dofs per node of ale element
        const double                    fac,              //!< domain-integration factor
        const double                    rhsfac,           //!< time-integration factor for rhs times domain-integration factor
        const double                    densnp,           //!< density at time_(n+1)
        const double                    J,                //!< determinant of Jacobian det(dx/ds)
        const LINALG::Matrix<nsd_,1>&   gradphi,          //!< scalar gradient at Gauss point
        const LINALG::Matrix<nsd_,1>&   convelint         //!< convective velocity
    );

    //! standard Galerkin convective term (only shapederivative term mesh)
    void ApplyShapeDerivsConv(
        Epetra_SerialDenseMatrix&       emat,
        const int                       k,
        const double                    rhsfac,
        const double                    densnp,
        const double                    J,
        const LINALG::Matrix<nsd_,1>&   gradphi,
        const LINALG::Matrix<nsd_,1>&   convelint
    );

    //! standard Galerkin diffusive term (off diagonal/shapederivative term mesh)
    virtual void CalcDiffODMesh(
        Epetra_SerialDenseMatrix&           emat,             //!< element current to be filled
        const int                           k,                //!< index of current scalar
        const int                           ndofpernodemesh,  //!< number of dofs per node of ale element
        const double                        fac,              //!< domain-integration factor
        const double                        rhsfac,           //!< time-integration factor for rhs times domain-integration factor
        const double                        J,                //!< determinant of Jacobian det(dx/ds)
        const LINALG::Matrix<nsd_,1>&       gradphi,          //!< scalar gradient at Gauss point
        const LINALG::Matrix<nsd_,1>&       convelint,        //!< convective velocity
        const LINALG::Matrix<1,nsd_*nen_>&  dJ_dmesh          //!< derivative of det(dx/ds) w.r.t. mesh displacement
    );

    //! reactive terms (standard Galerkin) (off diagonal/shapederivative term mesh)
    virtual void CalcReactODMesh(
        Epetra_SerialDenseMatrix&           emat,             //!< element current to be filled
        const int                           k,                //!< index of current scalar
        const int                           ndofpernodemesh,  //!< number of dofs per node of ale element
        const double                        rhsfac,           //!< time-integration factor for rhs times domain-integration factor
        const double                        rea_phi,          //!< reactive term
        const double                        J,                //!< determinant of Jacobian det(dx/ds)
        const LINALG::Matrix<1,nsd_*nen_>&  dJ_dmesh          //!< derivative of det(dx/ds) w.r.t. mesh displacement
    );

    //! calculation of linearized mass (off diagonal terms fluid)
    virtual void CalcLinMassODFluid(
      Epetra_SerialDenseMatrix&          emat,              //!< element matrix to be filled
      const int                          k,                 //!< index of current scalar
      const int                          ndofpernodemesh,   //!< number of dofs per node of ale element
      const double                       rhsfac,            //!< time-integration factor for rhs times domain-integration factor
      const double                       fac,               //!< domain-integration factor
      const double                       densam,            //!< density at time_(n+am)
      const double                       densnp,            //!< density at time_(n+1)
      const double                       phinp,             //!< scalar at time_(n+1)
      const double                       hist               //!< history of time integartion
    );

    //! standard Galerkin transient, old part of rhs and source term (off diagonal terms fluid)
    virtual void CalcHistAndSourceODFluid(
      Epetra_SerialDenseMatrix&          emat,             //!< element current to be filled
      const int                          k,                //!< index of current scalar
      const int                          ndofpernodemesh,  //!< number of dofs per node of ale element
      const double                       fac,              //!< domain-integration factor
      const double                       rhsint,           //!< rhs at Gauss point
      const double                       densnp            //!< density
    );

    //! standard Galerkin reactive term (off diagonal terms fluid)
    virtual void CalcReactODFluid(
      Epetra_SerialDenseMatrix&          emat,             //!< element current to be filled
      const int                          k,                //!< index of current scalar
      const int                          ndofpernodemesh,  //!< number of dofs per node of ale element
      const double                       rhsfac,           //!< time-integration factor for rhs times domain-integration factor
      const double                       rea_phi           //!< reactive term
    );

    //! standard Galerkin diffusive term (off diagonal terms fluid)
    virtual void CalcDiffODFluid(
        Epetra_SerialDenseMatrix&           emat,             //!< element current to be filled
        const int                           k,                //!< index of current scalar
        const int                           ndofpernodemesh,  //!< number of dofs per node of ale element
        const double                        rhsfac,           //!< time-integration factor for rhs times domain-integration factor
        const LINALG::Matrix<nsd_,1>&       gradphi           //!< scalar gradient at Gauss point
    );

    /*========================================================================*/
    //! @name dofs and nodes
    /*========================================================================*/

    //! number of dof per node
    const int numdofpernode_;
    //! number of transported scalars (numscal_ <= numdofpernode_)
    const int numscal_;

    /*========================================================================*/
    //! @name parameter lists
    /*========================================================================*/

    //! pointer to general scalar transport parameter class
    DRT::ELEMENTS::ScaTraEleParameterStd* scatrapara_;
    //! pointer to turbulence parameter class
    DRT::ELEMENTS::ScaTraEleParameterTurbulence* turbparams_;
    //! pointer to time integration parameter class
    DRT::ELEMENTS::ScaTraEleParameterTimInt* scatraparatimint_;

    /*========================================================================*/
    //! @name manager classes for efficient application to various problems
    /*========================================================================*/

    //! manager for diffusion
    Teuchos::RCP<ScaTraEleDiffManager> diffmanager_;
    //! manager for reaction
    Teuchos::RCP<ScaTraEleReaManager> reamanager_;

    /*========================================================================*/
    //! @name scalar degrees of freedom and related
    /*========================================================================*/

    //! state variables at t_(n)
    std::vector<LINALG::Matrix<nen_,1> > ephin_;
    //! state variables at t_(n+1) or t_(n+alpha_F)
    std::vector<LINALG::Matrix<nen_,1> > ephinp_;
    //! history vector of transported scalars
    std::vector<LINALG::Matrix<nen_,1> > ehist_;
    //! fine-scale solution ?
    std::vector<LINALG::Matrix<nen_,1> > fsphinp_;

    /*========================================================================*/
    //! @name velocity, pressure, and related
    /*========================================================================*/

    //! for the handling of rotationally symmetric periodic boundary conditions
    Teuchos::RCP<FLD::RotationallySymmetricPeriodicBC<distype,nsd_+1,DRT::ELEMENTS::Fluid::none> > rotsymmpbc_;
    //! nodal velocity values at t_(n+1) or t_(n+alpha_F)
    LINALG::Matrix<nsd_,nen_> evelnp_;
    //! nodal convective velocity values at t_(n+1) or t_(n+alpha_F)
    LINALG::Matrix<nsd_,nen_> econvelnp_;
    //! nodal fine-scale velocity values at t_(n+1) or t_(n+alpha_F)
    //! required for fine-scale subgrid diffusivity of type smagorinsky_small and multifractal subgrid scales
    LINALG::Matrix<nsd_,nen_> efsvel_;
    //! nodal acceleration values at t_(n+1) or t_(n+alpha_F)
    LINALG::Matrix<nsd_,nen_> eaccnp_;
    //! nodal displacement values for ALE
    LINALG::Matrix<nsd_,nen_> edispnp_;
    //! nodal pressure values at t_(n+1) or t_(n+alpha_F)
    LINALG::Matrix<nen_,1> eprenp_;


    /*========================================================================*/
    //! @name element coefficients and related extracted in evaluate
    /*========================================================================*/

    //! turbulent Prandtl number
    double tpn_;

    /*========================================================================*/
    //! @name Galerkin approximation and related
    /*========================================================================*/

    //! coordinates of current integration point in reference coordinates
    LINALG::Matrix<nsd_ele_,1> xsi_;
    //! node coordinates
    LINALG::Matrix<nsd_,nen_> xyze_;
    //! array for shape functions
    LINALG::Matrix<nen_,1> funct_;
    //! array for dual shape functions
    LINALG::Matrix<nen_,1> dual_funct_;
    //! array for shape function derivatives w.r.t r,s,t
    LINALG::Matrix<nsd_,nen_> deriv_;
    //! array for second derivatives of shape function w.r.t r,s,t
    LINALG::Matrix<numderiv2_,nen_> deriv2_;
    //! global derivatives of shape functions w.r.t x,y,z
    LINALG::Matrix<nsd_,nen_> derxy_;
    //! global second derivatives of shape functions w.r.t x,y,z
    LINALG::Matrix<numderiv2_,nen_> derxy2_;

    //! transposed jacobian "dx/ds"
    LINALG::Matrix<nsd_,nsd_> xjm_;
    //! inverse of transposed jacobian "ds/dx"
    LINALG::Matrix<nsd_,nsd_> xij_;
    //! 2nd derivatives of coord.-functions w.r.t r,s,t
    LINALG::Matrix<numderiv2_,nsd_> xder2_;

    //! bodyforce in element nodes
    std::vector<LINALG::Matrix<nen_,1> > bodyforce_;
//
    //! weights for nurbs elements
    LINALG::Matrix<nen_,1> weights_;
    //! knot vector for nurbs elements
    std::vector<Epetra_SerialDenseVector> myknots_;

    /*========================================================================*/
    //! @name can be very useful
    /*========================================================================*/

    int eid_;
    DRT::Element* ele_;

    //! variable manager for Gauss point values
    Teuchos::RCP<ScaTraEleInternalVariableManager<nsd_,nen_> > scatravarmanager_;
  };


  /// Scatra diffusion manager
  /*!
      This is a basic class to handle diffusion. It exclusively contains
      the isotropic diffusion coefficient. For anisotropic diffusion or
      more advanced diffusion laws, e.g., nonlinear ones, a derived class
      has to be constructed in the problem-dependent subclass for element
      evaluation.
  */
  class ScaTraEleDiffManager
  {
  public:

    ScaTraEleDiffManager(int numscal):
      diff_(numscal,0.0),
      sgdiff_(numscal,0.0)
      {return;}

    virtual ~ScaTraEleDiffManager(){};

    //! Set the isotropic diffusion coefficient
    virtual void SetIsotropicDiff(const double & diff, const int & k)
    {
//      if (diff < 0.0) dserror("negative (physical) diffusivity: %f",0,diff);

      diff_[k] = diff;
      return;
    }

    //! Set the isotropic subgrid diffusion coefficient
    virtual void SetIsotropicSubGridDiff(const double sgdiff, const int k)
    {
      diff_[k] += sgdiff;
      sgdiff_[k] += sgdiff;
      return;
    }

    //! Return the stored isotropic diffusion coefficients
    virtual std::vector<double> GetIsotropicDiff(){return diff_;}
    virtual double GetIsotropicDiff(const int k){return diff_[k];}

    //! Return the stored sub-grid diffusion coefficient
    virtual double GetSubGrDiff(const int k){return sgdiff_[k];}

  protected:

    //! scalar diffusion coefficient
    std::vector<double> diff_;

    //! subgrid diffusion coefficient
    std::vector<double> sgdiff_;
  };


  /// Scatra reaction manager
  /*!
      This is a basic class to handle constant reaction terms. For even more advanced reaction laws,
      a derived class has to be constructed in the problem-dependent subclass for element evaluation.
  */
  class ScaTraEleReaManager
    {
    public:

      ScaTraEleReaManager(int numscal):
        include_me_(false), // is reaction manager active?
        reacoeff_(numscal,0.0) // size of vector + initialized to zero
        {return;}

      virtual ~ScaTraEleReaManager(){};

      //! @name set routines

      //! Clear everything and resize to length numscal
      virtual void Clear(int numscal)
      {
        //clear
        reacoeff_.resize(0);
        //resize
        reacoeff_.resize(numscal,0.0);
        include_me_=false;
        return;
      }

      //! Set the reaction coefficient
      void SetReaCoeff(const double reacoeff, const int k)
      {
        // NOTE: it is important that this reaction coefficient set here does not depend on ANY concentration.
        // If so, use the advanced reaction framework to get a proper linearisation.
        reacoeff_[k] = reacoeff;
        if( reacoeff!=0.0 )
          include_me_=true;

        return;
      }

      //! @name access routines

      //! Return the reaction coefficient
      double GetReaCoeff(const int k) const {return reacoeff_[k];}

      //! Return the stabilization coefficient
      virtual double GetStabilizationCoeff(const int k, const double phinp_k) const {return reacoeff_[k];}

      //! return flag: reaction activated
      bool Active() const {return include_me_;}

    protected:
      //! flag for reaction
      bool include_me_;

    private:

      //NOTE: it is important that this reaction coefficient set here does not depend on ANY concentration.
      //! scalar reaction coefficient
      std::vector<double> reacoeff_;
    };

  /// ScaTraEleInternalVariableManager implementation
  /*!
    This class manages all internal variables needed for the evaluation of an element.

    All formulation-specific internal variables are stored and managed by a class derived from this class.
  */
  template<int NSD, int NEN>
  class ScaTraEleInternalVariableManager
  {
    public:

    ScaTraEleInternalVariableManager(int numscal):
      phinp_(numscal,0.0),
      phin_(numscal,0.0),
      convelint_(true),
      conv_(true),
      gradphi_(numscal),
      conv_phi_(numscal,0.0),
      hist_(numscal,0.0),
      numscal_(numscal)
      {return;}

    virtual ~ScaTraEleInternalVariableManager(){};

    // compute and set internal variables
    void SetInternalVariables(
        const LINALG::Matrix<NEN,1>&                 funct,       //! array for shape functions
        const LINALG::Matrix<NSD,NEN>&               derxy,       //! global derivatives of shape functions w.r.t x,y,z
        const std::vector<LINALG::Matrix<NEN,1> >&   ephinp,      //! scalar at t_(n+1) or t_(n+alpha_F)
        const std::vector<LINALG::Matrix<NEN,1> >&   ephin,       //! scalar at t_(n)
        const LINALG::Matrix<NSD,NEN>&               econvelnp,   //! nodal convective velocity values at t_(n+1) or t_(n+alpha_F)
        const std::vector<LINALG::Matrix<NEN,1> >&   ehist        //! history vector of transported scalars
        )
    {
      // convective velocity
      convelint_.Multiply(econvelnp,funct);
      // convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
      conv_.MultiplyTN(derxy,convelint_);

      for (int k = 0; k < numscal_; ++k)
      {
        // calculate scalar at t_(n+1) or t_(n+alpha_F)
        phinp_[k] = funct.Dot(ephinp[k]);
        // calculate scalar at t_(n)
        phin_[k] = funct.Dot(ephin[k]);
        // spatial gradient of current scalar value
        gradphi_[k].Multiply(derxy,ephinp[k]);
        // convective term
        conv_phi_[k] = convelint_.Dot(gradphi_[k]);
        // history data (or acceleration)
        hist_[k] = funct.Dot(ehist[k]);
      }

      return;
    };

  /*========================================================================*/
  //! @name return methods for internal variables
  /*========================================================================*/

    //! return scalar values at t_(n+1) or t_(n+alpha_F)
    virtual const std::vector<double>& Phinp() const {return phinp_;};
    //! return scalar value at t_(n+1) or t_(n+alpha_F)
    virtual const double& Phinp(const int k) const {return phinp_[k];};
    //! return scalar values at t_(n)
    virtual const std::vector<double>& Phin() const {return phin_;};
    //! return scalar value at t_(n)
    virtual const double& Phin(const int k) const {return phin_[k];};
    //! return convective velocity
    virtual const LINALG::Matrix<NSD,1>& ConVel(const int k) const {return convelint_;};
    //! return convective part in convective form
    virtual const LINALG::Matrix<NEN,1>& Conv(const int k) const {return conv_;};
    //! return spatial gradient of all scalar values
    virtual const std::vector<LINALG::Matrix<NSD,1> >& GradPhi() const {return gradphi_;};
    //! return spatial gradient of current scalar value
    virtual const LINALG::Matrix<NSD,1>& GradPhi(const int k) const {return gradphi_[k];};
    //! return convective term of current scalar value
    virtual const double& ConvPhi(const int k) const {return conv_phi_[k];};
    //! return history term of current scalar value
    virtual const double& Hist(const int k) const {return hist_[k];};

    /*========================================================================*/
    //! @name set methods for internal variables
    /*========================================================================*/
    // in a perfect world, all values should have been set in SetInternalVariables() and
    // related methods. Single(!) variables should not be manipulated, because this
    // easily leads to inconsistencies. However, this needs to be done in certain places
    // at the moment. Therefore the following set methods. Try to circumvent
    // using them and do not add new ones, if possible.
    //TODO: reduce number of set methods

    //! set spatial gradient of current scalar value
    virtual void SetGradPhi(const int k,LINALG::Matrix<NSD,1>& gradphi) {gradphi_[k] = gradphi;};
    //! set convective term of current scalar value
    virtual void SetConvPhi(const int k,double conv_phi) {conv_phi_[k] = conv_phi;};
    //! set convective term of current scalar value
    virtual void AddToConvPhi(const int k,double conv_phi) {conv_phi_[k] += conv_phi;};
    //! set convective term of current scalar value
    virtual void ScaleConvPhi(const int k,double scale) {conv_phi_[k] *= scale;};

  protected:

  /*========================================================================*/
  //! @name internal variables evaluated at element center or Gauss point
  /*========================================================================*/

    //! scalar at t_(n+1) or t_(n+alpha_F)
    std::vector<double> phinp_;
    //! scalar at t_(n)
    std::vector<double> phin_;
    //! convective velocity
    LINALG::Matrix<NSD,1> convelint_;
    //! convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
    LINALG::Matrix<NEN,1> conv_;
    //! spatial gradient of current scalar value
    std::vector<LINALG::Matrix<NSD,1> > gradphi_;
    //! convective term
    std::vector<double> conv_phi_;
    //! history data (or acceleration)
    std::vector<double> hist_;

  /*========================================================================*/
  //! @name number of scalars
  /*========================================================================*/

    //! number of scalars
    const int numscal_;

  };

  }
}


#endif
