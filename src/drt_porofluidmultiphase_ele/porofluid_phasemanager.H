/*----------------------------------------------------------------------*/
/*!
 \file porofluid_phasemanager.H

 \brief manager class for handling the phases and their dofs on element level

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROFLUIDMULTPHASE_ELE_POROFLUID_PHASEMANAGER_H_
#define SRC_DRT_POROFLUIDMULTPHASE_ELE_POROFLUID_PHASEMANAGER_H_

#include "porofluidmultiphase_ele_action.H"

#include "../drt_inpar/inpar_material.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "Teuchos_RCP.hpp"
#include <vector>


//forward declarations
class Epetra_SerialDenseMatrix;
class Epetra_SerialDenseVector;

namespace MAT
{
  class Material;
  class StructPoro;
  class FluidPoroSinglePhase;
  class FluidPoroMultiPhase;
}

namespace DRT
{
class Element;

namespace ELEMENTS
{
class PoroFluidMultiPhaseEleParameter;

namespace POROFLUIDMANAGER
{
class VariableManagerMinAccess;

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/

/*!
\brief interface to phase manager classes

These classes manage all accesses to variables at gauss points for the different
phases. The phase manager is responsible for the choice of primary variables.
As we can choose pressures or saturation as primary variables this class
will make sure everything is accessed correctly. Therefore it differentiates
between the phases (as each phase can have a different physical primary variable)

(in contrast to variable manager, which manages the 'generic' primary variables. Note that
these two manager classes are closely related).

The idea is, that there are the methods Setup(...) and EvaluateGPState(..), which
need to be called before evaluation.
All other methods are (more or less) constant access methods.

All implementations are derived from this class.

Two factory methods (in fact, only one as one calls the other), provide
the specialization depending on the evaluation action.

\author vuong
*/

class PhaseManagerInterface
  {
  public:
    //! constructor
    PhaseManagerInterface(){};

    //! destructor
    virtual ~PhaseManagerInterface(){};

    //! factory method
    static Teuchos::RCP< DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerInterface > CreatePhaseManager(
        const DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter& para,
        int nsd,
        INPAR::MAT::MaterialType mattype,
        const POROFLUIDMULTIPHASE::Action&     action,
        int totalnumdofpernode,
        int numfluidphases
        );

    //! factory method
    static Teuchos::RCP< DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerInterface > WrapPhaseManager(
        const DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter& para,
        int nsd,
        INPAR::MAT::MaterialType mattype,
        const POROFLUIDMULTIPHASE::Action&     action,
        Teuchos::RCP< PhaseManagerInterface > corephasemanager
        );

    //! setup (matnum is the material number of the porofluid-material on the current element)
    //! default is set to zero, if called from a porofluidmultiphase-element
    //! otherwise it has to be explicitly passed from the caller
    virtual void Setup(const DRT::Element* ele, const int matnum = 0)=0;

    //! evaluate pressures, saturations and derivatives at GP (matnum is the material number of the porofluid-material on the current element)
    //! default is set to zero, if called from a porofluidmultiphase-element
    //! otherwise it has to be explicitly passed from the caller
    virtual void EvaluateGPState(
        double                 J,
        const VariableManagerMinAccess& varmanager,
        const int matnum = 0)=0;

    //! clear the states
    virtual void ClearGPState()=0;

    //! check for reactions
    virtual bool IsReactive(int phasenum) const=0;

    //! get scalar to phase ID
    virtual int ScalarToPhaseID(int iscal) const=0;

    //! check if EvaluateGPState() was called
    virtual void CheckIsEvaluated() const=0;

    //! check if Setup() was called
    virtual void CheckIsSetup() const=0;

    //! @name Access methods

    //! get the number of phases
    virtual int NumFluidPhases() const=0;

    //! get the total number of dofs (number of phases + number of volume fractions)
    virtual int TotalNumDof() const=0;

    //! get solid pressure
    virtual double SolidPressure() const=0;

    //! recalculate solid pressure
    virtual void RecalculateSolidPressure(const double porosity, const std::vector<double>& volfracpressure)=0;

    //! get saturation of phase 'phasenum'
    virtual double Saturation(int phasenum) const=0;

    //! get pressure of phase 'phasenum'
    virtual double Pressure(int phasenum) const=0;

    //! get saturation of all phases
    virtual const std::vector<double>& Saturation() const =0;

    //! get volfracs of all phases
    virtual const std::vector<double>& VolFrac() const =0;

    //! get volfrac of volfrac 'volfracnum'
    virtual double VolFrac(int volfracnum) const =0;

    //! get sum of additional volume fractions
    virtual double SumAddVolFrac() const=0;

    //! get pressure of all phases
    virtual const std::vector<double>& Pressure() const =0;

    //! get bulk modulus of phase 'phasenum'
    virtual double InvBulkmodulus(
        int phasenum) const=0;

    //! check if fluid phase 'phasenum' is incompressible (very low compressibility < 1e-14)
    virtual bool IncompressibleFluidPhase(
        int phasenum) const=0;

    //! get inverse bulk modulus of solid phase
    virtual double InvBulkmodulusSolid() const=0;

    //! check if solid is incompressible (either very low compressibility < 1e-14 or MAT_PoroLawConstant)
    virtual bool IncompressibleSolid() const = 0;

    //! get density of phase 'phasenum'
    virtual double Density(
        int phasenum) const=0;

    //! get the current element the manager was set up with
    virtual const DRT::Element* Element() const = 0;

    //! get porosity
    virtual double Porosity() const=0;

    //! get Jacobian of deformation gradient
    virtual double JacobianDefGrad() const=0;

    //! get derivative of porosity wrt JacobianDefGrad
    virtual double PorosityDerivWrtJacobianDefGrad() const=0;

    //! get derivative of porosity w.r.t. DOF 'doftoderive'
    virtual double PorosityDeriv(
            int doftoderive) const=0;

    //! check if porosity depends on fluid (pressure)
    virtual bool PorosityDependsOnFluid() const=0;

    //! check if porosity depends on structure (basically Jacobian of def. gradient)
    virtual bool PorosityDependsOnStruct() const=0;

    //! get derivative of saturation of phase 'phasenum' w.r.t. DOF 'doftoderive'
    virtual double SaturationDeriv(
            int phasenum,
            int doftoderive) const=0;

    //! get 2nd derivative of saturation of phase 'phasenum' w.r.t. DOF 'doftoderive'
    virtual double SaturationDerivDeriv(
            int phasenum,
            int firstdoftoderive,
            int seconddoftoderive) const=0;

    //! get derivative of pressure of phase 'phasenum' w.r.t. DOF 'doftoderive'
    virtual double PressureDeriv(
            int phasenum,
            int doftoderive) const=0;

    //! get derivative of solid pressure  w.r.t. DOF 'doftoderive'
    virtual double SolidPressureDeriv(
            int doftoderive) const=0;

    //! get derivative of pressure of phase 'phasenum'
    //! w.r.t. DOF 'doftoderive' (first derivative)
    //! and w.r.t. DOF 'doftoderive2' (second derivative)
    virtual double SolidPressureDerivDeriv(
            int doftoderive,
            int doftoderive2) const=0;

    //! get the reaction term
    virtual double ReacTerm(int phasenum) const=0;

    //! get the derivative of the reaction term
    virtual double ReacDeriv(int phasenum, int doftoderive) const=0;

    //! get total number of scalars in system
    virtual int NumScal() const = 0;

    //! get the derivative of the reaction term w.r.t. scalar 'scaltoderive'
    virtual double ReacDerivScalar(int phasenum, int scaltoderive) const=0;

    //! get the derivative of the reaction term w.r.t. porosity
    virtual double ReacDerivPorosity(int phasenum) const=0;

    //! get the diffusion tensor
    virtual void PermeabilityTensor(int phasenum,  LINALG::Matrix<3,3>& permeabilitytensor) const=0;
    //! get the diffusion tensor
    virtual void PermeabilityTensor(int phasenum,  LINALG::Matrix<2,2>& permeabilitytensor) const=0;
    //! get the diffusion tensor
    virtual void PermeabilityTensor(int phasenum,  LINALG::Matrix<1,1>& permeabilitytensor) const=0;

    //! check for constant relpermeability
    virtual bool HasConstantRelPermeability(int phasenum) const=0;
    //! get relative diffusivity of phase
    virtual double RelPermeability(int phasenum) const=0;
    //! get derivative of relative permeability of phase
    virtual double RelPermeabilityDeriv(int phasenum) const=0;

    //! check for constant dynamic viscosity
    virtual bool HasConstantDynViscosity(int phasenum) const=0;

    //! get dynamic viscosity of phase (matnum is the material number of the porofluid-material on the current element)
    //! default is set to zero, if called from a porofluidmultiphase-element
    //! otherwise it has to be explicitly passed from the caller
    virtual double DynViscosity(int phasenum, double abspressgrad, int matnum = 0) const=0;
    //! get dynamic viscosity of phase
    virtual double DynViscosity(const MAT::Material& material, int phasenum, double abspressgrad) const=0;
    //! get derivative of dynamic viscosity of phase
    virtual double DynViscosityDeriv(int phasenum, double abspressgrad) const=0;
    //! get derivative dynamic viscosity of phase
    virtual double DynViscosityDeriv(const MAT::Material& material, int phasenum, double abspressgrad) const=0;

    //! get the diffusion tensor
    virtual void DiffTensorVolFrac(int volfracnum,  LINALG::Matrix<3,3>& difftensorvolfrac) const=0;
    //! get the diffusion tensor
    virtual void DiffTensorVolFrac(int volfracnum,  LINALG::Matrix<2,2>& difftensorvolfrac) const=0;
    //! get the diffusion tensor
    virtual void DiffTensorVolFrac(int volfracnum,  LINALG::Matrix<1,1>& difftensorvolfrac) const=0;

    //! get density of volume frac 'volfracnum'
    virtual double VolFracDensity(int volfracnum) const=0;

    //! get pressure of volume frac 'volfracnum'
    virtual double VolFracPressure(int volfracnum) const=0;

    //! check if volume frac 'volfracnum' has additional scalar dependent flux
    virtual bool HasAddScalarDependentFlux(int volfracnum) const=0;

    //! check if volume frac 'volfracnum' has additional scalar dependent flux of scalar 'iscal'
    virtual bool HasAddScalarDependentFlux(int volfracnum, int iscal) const=0;

    //! check if volume frac 'volfracnum' has receptor kinetic-law of scalar 'iscal'
    //! see: Anderson, A. R. A. & Chaplain, M. A. J.
    //       Continuous and discrete mathematical models of tumor-induced angiogenesis
    virtual bool HasReceptorKineticLaw(int volfracnum, int iscal) const=0;

    //! return scalar diffusivity of of scalar 'iscal' of volume fraction 'volfracnum'
    virtual double ScalarDiff(int volfracnum, int iscal) const=0;

    //! return omega half of scalar 'iscal' of volume fraction 'volfracnum' for receptor kinetic law
    //! see: Anderson, A. R. A. & Chaplain, M. A. J.
    //       Continuous and discrete mathematical models of tumor-induced angiogenesis
    virtual double OmegaHalf(int volfracnum, int iscal) const=0;

    //@}

  };

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/

/*!
\brief standard phase manager, holding pressures and saturations

This class is a minimal version of a phase manager. It basically provides
access to the primary variables (saturation, pressure, etc.) and quantities
that can be accessed via the material in almost all cases (bulkmodulus, ...).

\author vuong
*/
class PhaseManagerCore : public PhaseManagerInterface
  {
  public:
    //! constructor
    explicit PhaseManagerCore(int totalnumdofpernode, int numphase);

    //! copy constructor
    PhaseManagerCore(const PhaseManagerCore& old);

    //! destructor
    virtual ~PhaseManagerCore(){};

    //! setup (matnum is the material number of the porofluid-material on the current element)
    //! default is set to zero, if called from a porofluidmultiphase-element
    //! otherwise it has to be explicitly passed from the caller
    void Setup(const DRT::Element* ele, const int matnum = 0);

    //! evaluate pressures, saturations and derivatives at GP (matnum is the material number of the porofluid-material on the current element)
    //! default is set to zero, if called from a porofluidmultiphase-element
    //! otherwise it has to be explicitly passed from the caller
    void EvaluateGPState(
        double                          J,
        const VariableManagerMinAccess& varmanager,
        const int matnum = 0);

    //! clear the states
    void ClearGPState();

    //! check for reactions
    bool IsReactive(int phasenum) const {return false;};

    //! @name Access methods

    //! get the number of phases
    int NumFluidPhases() const {return numfluidphases_;};

    //! get the total number of dofs (number of phases + number of volume fractions)
    int TotalNumDof() const {return totalnumdofpernode_;};

    //! get solid pressure
    double SolidPressure() const;

    //! recalculate solid pressure
    void RecalculateSolidPressure(const double porosity, const std::vector<double>& volfracpressure);

    //! get saturation of phase 'phasenum'
    double Saturation(int phasenum) const;

    //! get pressure of phase 'phasenum'
    double Pressure(int phasenum) const;

    //! get saturation of all phases
    const std::vector<double>& Saturation() const;

    //! get volfracs of all phases
    const std::vector<double>& VolFrac() const;

    //! get volfrac of volfrac 'volfracnum'
    double VolFrac(int volfracnum) const;

    //! get sum of additional volume fractions
    double SumAddVolFrac() const;

    //! get pressure of all phases
    const std::vector<double>& Pressure() const;

    //! get bulk modulus of phase 'phasenum'
    double InvBulkmodulus(
        int phasenum) const;

    //! check if fluid phase 'phasenum' is incompressible (very low compressibility < 1e-14)
    bool IncompressibleFluidPhase(
        int phasenum) const;

    //! get inverse bulk modulus of solid phase
    double InvBulkmodulusSolid() const;

    //! check if solid is incompressible (either very low compressibility < 1e-14 or MAT_PoroLawConstant)
    bool IncompressibleSolid() const;

    //! get density of phase 'phasenum'
    double Density(
        int phasenum) const;

    //! get the current element the manager was set up with
    virtual const DRT::Element* Element() const {return ele_;};

    //@}

    //! check if EvaluateGPState() was called
    void CheckIsEvaluated() const
    {
      if(not isevaluated_)
        dserror("Gauss point states have not been set!");
    }

    //! check if EvaluateGPState() was called
    void CheckIsSetup() const
    {
      if(not issetup_)
        dserror("Setup() was not called!");
    }

    //! get porosity
    virtual double Porosity() const
    {dserror("Porosity not available for this phase manager!");return 0.0;};

    //! get porosity
    virtual double JacobianDefGrad() const
    {dserror("JacobianDefGrad not available for this phase manager!");return 0.0;};

    //! get derivative of porosity wrt JacobianDefGrad
    virtual double PorosityDerivWrtJacobianDefGrad() const
    {dserror("Derivative of Porosity w.r.t. JacobianDefGrad not available for this phase manager!");return 0.0;};

    //! get derivative of porosity w.r.t. DOF 'doftoderive'
    virtual double PorosityDeriv(
            int doftoderive) const
    {dserror("Derivative of porosity not available for this phase manager!");return 0.0;};

    //! check if porosity depends on fluid (pressure)
    virtual bool PorosityDependsOnFluid() const
    {dserror("PorosityDependsOnFluid() not available for this phase manager!");return false;};

    //! check if porosity depends on structure (basically Jacobian of def. gradient)
    virtual bool PorosityDependsOnStruct() const
    {dserror("PorosityDependsOnStruct() not available for this phase manager!");return false;};

    //! get derivative of saturation of phase 'phasenum' w.r.t. DOF 'doftoderive'
    virtual double SaturationDeriv(
            int phasenum,
            int doftoderive) const
    {dserror("Derivative of saturation not available for this phase manager!");return 0.0;};

    //! get 2nd derivative of saturation of phase 'phasenum' w.r.t. DOF 'doftoderive'
    virtual double SaturationDerivDeriv(
            int phasenum,
            int firstdoftoderive,
            int seconddoftoderive) const
    {dserror("2nd Derivative of saturation not available for this phase manager!");return 0.0;};

    //! get derivative of pressure of phase 'phasenum' w.r.t. DOF 'doftoderive'
    virtual double PressureDeriv(
            int phasenum,
            int doftoderive) const
    {dserror("Derivative of pressure not available for this phase manager!");return 0.0;};

    //! get derivative of solid pressure  w.r.t. DOF 'doftoderive'
    virtual double SolidPressureDeriv(
            int doftoderive) const
    {dserror("Derivative of solid pressure not available for this phase manager!");return 0.0;};

    //! get derivative of pressure of phase 'phasenum'
    //! w.r.t. DOF 'doftoderive' (first derivative)
    //! and w.r.t. DOF 'doftoderive2' (second derivative)
    virtual double SolidPressureDerivDeriv(
            int doftoderive,
            int doftoderive2) const
    {dserror("Second derivative of solid pressure not available for this phase manager!");return 0.0;};

    //! get the reaction term
    virtual double ReacTerm(int phasenum) const
    {dserror("Reaction term not available for this phase manager!");return 0.0;};

    //! get total number of scalars in system
    virtual int NumScal() const
    {dserror("Number of scalars not available for this phase manager"); return 0;};

    //! get the derivative of the reaction term
    virtual double ReacDeriv(int phasenum, int doftoderive) const
    {dserror("Reaction term derivative not available for this phase manager!");return 0.0;};

    //! get the derivative of the reaction term w.r.t. scalar 'scaltoderive'
    virtual double ReacDerivScalar(int phasenum, int scaltoderive) const
    {dserror("Reaction term derivative (scalar) not available for this phase manager!");return 0.0;};

    //! get the derivative of the reaction term w.r.t. porosity
    virtual double ReacDerivPorosity(int phasenum) const
    {dserror("Reaction term derivative (porosity) not available for this phase manager!");return 0.0;};

    //! get the diffusion tensor
    virtual void PermeabilityTensor(int phasenum,  LINALG::Matrix<3,3>& permeabilitytensor) const
        {dserror("Diffusion tensor (3D) not available for this phase manager!");};
    //! get the diffusion tensor
    virtual void PermeabilityTensor(int phasenum,  LINALG::Matrix<2,2>& permeabilitytensor) const
        {dserror("Diffusion tensor (2D) not available for this phase manager!");};
    //! get the diffusion tensor
    virtual void PermeabilityTensor(int phasenum,  LINALG::Matrix<1,1>& permeabilitytensor) const
        {dserror("Diffusion tensor (1D) not available for this phase manager!");};

    //! check for constant relpermeability
    virtual bool HasConstantRelPermeability(int phasenum) const
        {dserror("Check for Constant Relative Permeability not available for this phase manager!");return false;};
    //! get relative diffusivity of phase
    virtual double RelPermeability(int phasenum) const
        {dserror("Relative Diffusivity not available for this phase manager!");return 0.0;};
    //! get derivative of relative permeability of phase
    virtual double RelPermeabilityDeriv(int phasenum) const
        {dserror("Derivativ of relativ permeability not available for this phase manager!");return 0.0;};

    //! check for constant dynamic viscosity
    virtual bool HasConstantDynViscosity(int phasenum) const
        {dserror("Check for Constant Dynamic Viscosity not available for this phase manager!");return false;};
    //! get relative diffusivity of phase
    virtual double DynViscosity(int phasenum, double abspressgrad, int matnum = 0) const
        {dserror("Dynamic Viscosity not available for this phase manager!");return 0.0;};
    //! get dynamic viscosity of phase
    virtual double DynViscosity(const MAT::Material& material, int phasenum, double abspressgrad) const
        {dserror("Dynamic Viscosity not available for this phase manager!");return 0.0;};
    //! get derivative of dynamic viscosity of phase
    virtual double DynViscosityDeriv(int phasenum, double abspressgrad) const
        {dserror("Derivative of dynamic Viscosity not available for this phase manager!");return 0.0;};
    //! get derivative dynamic viscosity of phase
    double DynViscosityDeriv(const MAT::Material& material, int phasenum, double abspressgrad) const
        {dserror("Derivative of dynamic Viscosity not available for this phase manager!");return 0.0;};

    //! get the diffusion tensor
    virtual void DiffTensorVolFrac(int volfracnum,  LINALG::Matrix<3,3>& difftensorvolfrac) const
    {dserror("Diffusion tensor for volume fractions (3D) not available for this phase manager!");};
    //! get the diffusion tensor
    virtual void DiffTensorVolFrac(int volfracnum,  LINALG::Matrix<2,2>& difftensorvolfrac) const
    {dserror("Diffusion tensor for volume fractions (2D) not available for this phase manager!");};
    //! get the diffusion tensor
    virtual void DiffTensorVolFrac(int volfracnum,  LINALG::Matrix<1,1>& difftensorvolfrac) const
    {dserror("Diffusion tensor for volume fractions (1D) not available for this phase manager!");};

    //! get density of volume frac 'volfracnum'
    virtual double VolFracDensity(int volfracnum) const
    {dserror("VolFracDensity not available for this phase manager!");return 0.0;};

    //! get pressure of volume frac 'volfracnum'
    virtual double VolFracPressure(int volfracnum) const
    {dserror("VolFracPressure not available for this phase manager!");return 0.0;};

    //! check if volume frac 'volfracnum' has additional scalar dependent flux
    virtual bool HasAddScalarDependentFlux(int volfracnum) const
    {dserror("HasAddScalarDependentFlux not available for this phase manager!");return false;};

    //! check if volume frac 'volfracnum' has additional scalar dependent flux of scalar 'iscal'
    virtual bool HasAddScalarDependentFlux(int volfracnum, int iscal) const
    {dserror("HasAddScalarDependentFlux not available for this phase manager!");return false;};

    //! check if volume frac 'volfracnum' has receptor kinetic-law of scalar 'iscal'
    //! see: Anderson, A. R. A. & Chaplain, M. A. J.
    //       Continuous and discrete mathematical models of tumor-induced angiogenesis
    virtual bool HasReceptorKineticLaw(int volfracnum, int iscal) const
    {dserror("HasReceptorKineticLaw not available for this phase manager!");return false;};

    //! return scalar diffusivities of scalar 'iscal' of volume fraction 'volfracnum'
    virtual double ScalarDiff(int volfracnum, int iscal) const
    {dserror("ScalarDiff not available for this phase manager!");return 0.0;};

    //! return omega half of scalar 'iscal' of volume fraction 'volfracnum' for receptor kinetic law
    //! see: Anderson, A. R. A. & Chaplain, M. A. J.
    //       Continuous and discrete mathematical models of tumor-induced angiogenesis
    virtual double OmegaHalf(int volfracnum, int iscal) const
    {dserror("OmegaHalf not available for this phase manager!");return 0.0;};

    //! get scalar to phase ID
    virtual int ScalarToPhaseID(int iscal) const
    {dserror("ScalarToPhaseID not available for this phase manager!"); return 0;};

  private:
    //! total number of dofs per node (numfluidphases + numvolfrac)
    const int totalnumdofpernode_;
    //! number of phases
    const int numfluidphases_;

    //! generalized pressure
    std::vector<double> genpressure_;
    //! additional volume fraction
    std::vector<double> volfrac_;
    //! sum of additional volume fractions
    double sumaddvolfrac_;
    //! true pressure
    std::vector<double> pressure_;
    //! saturation
    std::vector<double> saturation_;
    //! densities
    std::vector<double> density_;
    //! solid pressure
    double solidpressure_;
    //! inverse bulk moduli of the fluid phases
    std::vector<double> invbulkmodulifluid_;
    //! inverse solid bulk modulus
    double invbulkmodulussolid_;

    //! the current element
    const DRT::Element* ele_;

    //! flag indicating of gauss point state has been set and evaluated
    bool isevaluated_;

    //! flag of Setup was called
    bool issetup_;

  };

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/

/*!
\brief wrapper class, base class for extensions to phase manager

The idea is to use a core phase manager (class PhaseManagerCore) and extend
it when the evaluation to be performed demands it. For this a decorator
 pattern is used, i.e. this class wraps another phase manager, extending it if necessary.

This is the base class for all decorators.

\author vuong
*/
class PhaseManagerDecorator : public PhaseManagerInterface
  {
  public:
    //! constructor
    explicit PhaseManagerDecorator(Teuchos::RCP< POROFLUIDMANAGER::PhaseManagerInterface > phasemanager)
    : phasemanager_(phasemanager)
    {};

    //! destructor
    virtual ~PhaseManagerDecorator(){};

    //! setup (matnum is the material number of the porofluid-material on the current element)
    //! default is set to zero, if called from a porofluidmultiphase-element
    //! otherwise it has to be explicitly passed from the caller
    virtual void Setup(const DRT::Element* ele, const int matnum = 0)
    {phasemanager_->Setup(ele, matnum);};

    //! check if EvaluateGPState() was called
    virtual void CheckIsEvaluated() const
    {phasemanager_->CheckIsEvaluated();};

    //! check if Setup() was called
    virtual void CheckIsSetup() const
    {phasemanager_->CheckIsSetup();};

    //! @name Access methods

    //! get derivative of porosity w.r.t. DOF 'doftoderive'
    virtual double PorosityDeriv(
            int doftoderive) const
    {return phasemanager_->PorosityDeriv(doftoderive);};

    //! check if porosity depends on fluid (pressure)
    virtual bool PorosityDependsOnFluid() const
    {return phasemanager_->PorosityDependsOnFluid();};

    //! check if porosity depends on structure (basically Jacobian of def. gradient)
    virtual bool PorosityDependsOnStruct() const
    {return phasemanager_->PorosityDependsOnStruct();};

    //! get derivative of saturation of phase 'phasenum' w.r.t. DOF 'doftoderive'
    virtual double SaturationDeriv(
            int phasenum,
            int doftoderive) const
    {return phasemanager_->SaturationDeriv(phasenum,doftoderive);};

    //! get 2nd derivative of saturation of phase 'phasenum' w.r.t. DOF 'doftoderive'
    virtual double SaturationDerivDeriv(
            int phasenum,
            int firstdoftoderive,
            int seconddoftoderive) const
    {return phasemanager_->SaturationDerivDeriv(phasenum,firstdoftoderive,seconddoftoderive);};

    //! get derivative of pressure of phase 'phasenum' w.r.t. DOF 'doftoderive'
    virtual double PressureDeriv(
            int phasenum,
            int doftoderive) const
    {return phasemanager_->PressureDeriv(phasenum,doftoderive);};

    //! get derivative of solid pressure  w.r.t. DOF 'doftoderive'
    virtual double SolidPressureDeriv(
            int doftoderive) const
    {return phasemanager_->SolidPressureDeriv(doftoderive);};

    //! get derivative of pressure of phase 'phasenum'
    //! w.r.t. DOF 'doftoderive' (first derivative)
    //! and w.r.t. DOF 'doftoderive2' (second derivative)
    virtual double SolidPressureDerivDeriv(
            int doftoderive,
            int doftoderive2) const
    {return phasemanager_->SolidPressureDerivDeriv(doftoderive,doftoderive2);};

    //! check if the current phase is involved in a reaction
    virtual bool IsReactive(int phasenum) const
    {return phasemanager_->IsReactive(phasenum);};

    //! get scalar to phase ID
    virtual int ScalarToPhaseID(int iscal) const
    {return phasemanager_->ScalarToPhaseID(iscal);};

    //! get the number of phases
    virtual int NumFluidPhases() const
    {return phasemanager_->NumFluidPhases();};

    //! get the number of dofs (number of phases + number of volfracs)
    virtual int TotalNumDof() const
    {return phasemanager_->TotalNumDof();};

    //! get solid pressure
    virtual double SolidPressure() const
    {return phasemanager_->SolidPressure();};

    //! recalculate solid pressure
    virtual void RecalculateSolidPressure(const double porosity, const std::vector<double>& volfracpressure)
    {return phasemanager_->RecalculateSolidPressure(porosity,volfracpressure);};

    //! get saturation of phase 'phasenum'
    virtual double Saturation(int phasenum) const
    {return phasemanager_->Saturation(phasenum);};

    //! get pressure of phase 'phasenum'
    double Pressure(int phasenum) const
    {return phasemanager_->Pressure(phasenum);};

    //! get saturation of all phases
    virtual const std::vector<double>& Saturation() const
    {return phasemanager_->Saturation();};

    //! get volfracs of all phases
    virtual const std::vector<double>& VolFrac() const
    {return phasemanager_->VolFrac();};

    //! get volfrac of volfrac 'volfracnum'
    virtual double VolFrac(int volfracnum) const
    {return phasemanager_->VolFrac(volfracnum);};

    //! get solid pressure
    virtual double SumAddVolFrac() const
    {return phasemanager_->SumAddVolFrac();};

    //! get pressure of all phases
    virtual const std::vector<double>& Pressure() const
    {return phasemanager_->Pressure();};

    //! get bulk modulus of phase 'phasenum'
    virtual double InvBulkmodulus(
        int phasenum) const
    {return phasemanager_->InvBulkmodulus(phasenum);};

    //! check if fluid phase 'phasenum' is incompressible (very low compressibility < 1e-14)
    virtual bool IncompressibleFluidPhase(
        int phasenum) const
    {return phasemanager_->IncompressibleFluidPhase(phasenum);};

    //! get inverse bulk modulus of solid phase
    double InvBulkmodulusSolid() const
    {return phasemanager_->InvBulkmodulusSolid();};

    //! check if solid is incompressible (either very low compressibility < 1e-14 or MAT_PoroLawConstant)
    bool IncompressibleSolid() const
    {return phasemanager_->IncompressibleSolid();};

    //! get porosity
    virtual double Porosity() const
    {return phasemanager_->Porosity();};

    //! get JacobianDefGrad
    virtual double JacobianDefGrad() const
    {return phasemanager_->JacobianDefGrad();};

    //! get derivative of porosity wrt JacobianDefGrad
    virtual double PorosityDerivWrtJacobianDefGrad() const
    {return phasemanager_->PorosityDerivWrtJacobianDefGrad();};

    //! get density of phase 'phasenum'
    virtual double Density(
        int phasenum) const
    {return phasemanager_->Density(phasenum);};

    //! get the current element the manager was set up with
    virtual const DRT::Element* Element() const
    {return phasemanager_->Element();};

    //! get the reaction term
    virtual double ReacTerm(int phasenum) const
    {return phasemanager_->ReacTerm(phasenum);};

    //! get total number of scalars in system
    virtual int NumScal() const
    {return phasemanager_->NumScal();};

    //! get the derivative of the reaction term
    virtual double ReacDeriv(int phasenum, int doftoderive) const
    {return phasemanager_->ReacDeriv(phasenum,doftoderive);};

    //! get the derivative of the reaction term w.r.t. scalar 'scaltoderive'
    virtual double ReacDerivScalar(int phasenum, int scaltoderive) const
    {return phasemanager_->ReacDerivScalar(phasenum,scaltoderive);};

    //! get the derivative of the reaction term w.r.t. porosity
    virtual double ReacDerivPorosity(int phasenum) const
    {return phasemanager_->ReacDerivPorosity(phasenum);};

    //! get the diffusion tensor
    virtual void PermeabilityTensor(int phasenum,  LINALG::Matrix<3,3>& permeabilitytensor) const
    {phasemanager_->PermeabilityTensor(phasenum,permeabilitytensor);};
    //! get the diffusion tensor
    virtual void PermeabilityTensor(int phasenum,  LINALG::Matrix<2,2>& permeabilitytensor) const
    {phasemanager_->PermeabilityTensor(phasenum,permeabilitytensor);};
    //! get the diffusion tensor
    virtual void PermeabilityTensor(int phasenum,  LINALG::Matrix<1,1>& permeabilitytensor) const
    {phasemanager_->PermeabilityTensor(phasenum,permeabilitytensor);};

    //! check for constant relpermeability
    virtual bool HasConstantRelPermeability(int phasenum) const
    {return phasemanager_->HasConstantRelPermeability(phasenum);};

    //! get relative diffusivity of phase
    virtual double RelPermeability(int phasenum) const
    {return phasemanager_->RelPermeability(phasenum);};

    //! get derivative of relative permeability of phase
    virtual double RelPermeabilityDeriv(int phasenum) const
    {return phasemanager_->RelPermeabilityDeriv(phasenum);};

    //! check for constant dynamic visosity
    virtual bool HasConstantDynViscosity(int phasenum) const
    {return phasemanager_->HasConstantDynViscosity(phasenum);};
    //! get dynamic viscosity of phase
    virtual double DynViscosity(int phasenum, double abspressgrad, int matnum = 0) const
    {return phasemanager_->DynViscosity(phasenum, abspressgrad, matnum = 0);};
    //! get dynamic viscosity of phase
    virtual double DynViscosity(const MAT::Material& material, int phasenum, double abspressgrad) const
    {return phasemanager_->DynViscosity(material, phasenum, abspressgrad);};
    //! get derivative of dynamic viscosity of phase
    virtual double DynViscosityDeriv(int phasenum, double abspressgrad) const
    {return phasemanager_->DynViscosityDeriv(phasenum, abspressgrad);};
    //! get derivative dynamic viscosity of phase
    double DynViscosityDeriv(const MAT::Material& material, int phasenum, double abspressgrad) const
    {return phasemanager_->DynViscosityDeriv(material, phasenum, abspressgrad);};

    //! get the diffusion tensor
    virtual void DiffTensorVolFrac(int volfracnum,  LINALG::Matrix<3,3>& difftensorvolfrac) const
    {phasemanager_->DiffTensorVolFrac(volfracnum,difftensorvolfrac);};
    //! get the diffusion tensor
    virtual void DiffTensorVolFrac(int volfracnum,  LINALG::Matrix<2,2>& difftensorvolfrac) const
    {phasemanager_->DiffTensorVolFrac(volfracnum,difftensorvolfrac);};
    //! get the diffusion tensor
    virtual void DiffTensorVolFrac(int volfracnum,  LINALG::Matrix<1,1>& difftensorvolfrac) const
    {phasemanager_->DiffTensorVolFrac(volfracnum,difftensorvolfrac);};


    //! get density of volume frac 'volfracnum'
    virtual double VolFracDensity(int volfracnum) const
    {return phasemanager_->VolFracDensity(volfracnum);};

    //! get pressure of volume frac 'volfracnum'
    virtual double VolFracPressure(int volfracnum) const
    {return phasemanager_->VolFracPressure(volfracnum);};

    //! check if volume frac 'volfracnum' has additional scalar dependent flux
    virtual bool HasAddScalarDependentFlux(int volfracnum) const
    {return phasemanager_->HasAddScalarDependentFlux(volfracnum);};

    //! check if volume frac 'volfracnum' has additional scalar dependent flux of scalar 'iscal'
    virtual bool HasAddScalarDependentFlux(int volfracnum, int iscal) const
    {return phasemanager_->HasAddScalarDependentFlux(volfracnum, iscal);};

    //! check if volume frac 'volfracnum' has receptor kinetic-law of scalar 'iscal'
    //! see: Anderson, A. R. A. & Chaplain, M. A. J.
    //       Continuous and discrete mathematical models of tumor-induced angiogenesis
    virtual bool HasReceptorKineticLaw(int volfracnum, int iscal) const
    {return phasemanager_->HasReceptorKineticLaw(volfracnum, iscal);};

    //! return scalar diffusivities of scalar 'iscal' of volume fraction 'volfracnum'
    virtual double ScalarDiff(int volfracnum, int iscal) const
    {return phasemanager_->ScalarDiff(volfracnum, iscal);};

    //! return omega half of scalar 'iscal' of volume fraction 'volfracnum' for receptor kinetic law
    //! see: Anderson, A. R. A. & Chaplain, M. A. J.
    //       Continuous and discrete mathematical models of tumor-induced angiogenesis
    virtual double OmegaHalf(int volfracnum, int iscal) const
    {return phasemanager_->OmegaHalf(volfracnum, iscal);};

    //@}

  protected:

    //! wrapped phase manager
    Teuchos::RCP< POROFLUIDMANAGER::PhaseManagerInterface > phasemanager_;

  };

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/

/*!
\brief wrapper class, extensions for derivatives

This class is a decorator for a phase manager, including the derivatives
of the pressures and saturations w.r.t. to the primary variables

\author vuong
*/
class PhaseManagerDeriv : public PhaseManagerDecorator
  {
  public:
    //! constructor
    explicit PhaseManagerDeriv(Teuchos::RCP< POROFLUIDMANAGER::PhaseManagerInterface > phasemanager);

    //! destructor
    virtual ~PhaseManagerDeriv(){};

    //! evaluate pressures, saturations and derivatives at GP (matnum is the material number of the porofluid-material on the current element)
    //! default is set to zero, if called from a porofluidmultiphase-element
    //! otherwise it has to be explicitly passed from the caller
    void EvaluateGPState(
        double                 J,
        const VariableManagerMinAccess& varmanager,
        const int matnum = 0);

    //! clear the states
    void ClearGPState();

    //! @name Access methods

    //! get derivative of saturation of phase 'phasenum' w.r.t. DOF 'doftoderive'
    double SaturationDeriv(
            int phasenum,
            int doftoderive) const;

    //! get 2nd derivative of saturation of phase 'phasenum' w.r.t. DOF 'doftoderive' (dS/dphi_{ii})
    double SaturationDerivDeriv(
            int phasenum,
            int firstdoftoderive,
            int seconddoftoderive) const;

    //! get derivative of pressure of phase 'phasenum' w.r.t. DOF 'doftoderive'
    double PressureDeriv(
            int phasenum,
            int doftoderive) const;

    //! get derivative of solid pressure  w.r.t. DOF 'doftoderive'
    double SolidPressureDeriv(
            int doftoderive) const;

    //! get derivative of pressure of phase 'phasenum'
    //! w.r.t. DOF 'doftoderive' (first derivative)
    //! and w.r.t. DOF 'doftoderive2' (second derivative)
    double SolidPressureDerivDeriv(
            int doftoderive,
            int doftoderive2) const;

    //@}

  private:

    //! derivative of true pressure w.r.t. degrees of freedom
    // first index: pressure, second index: dof
    Teuchos::RCP<Epetra_SerialDenseMatrix> pressurederiv_;
    //! derivative of saturations w.r.t. degrees of freedom
    // first index: saturation, second index: dof
    Teuchos::RCP<Epetra_SerialDenseMatrix> saturationderiv_;
    //! second derivative of saturations w.r.t. degrees of freedom
    Teuchos::RCP<std::vector<Epetra_SerialDenseMatrix> > saturationderivderiv_;

    //! derivative of solid pressure w.r.t. degrees of freedom
    Teuchos::RCP<Epetra_SerialDenseVector> solidpressurederiv_;
    //! second derivative of solid pressure w.r.t. degrees of freedom;
    Teuchos::RCP<Epetra_SerialDenseMatrix> solidpressurederivderiv_;


  };

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/

/*!
\brief wrapper class, extensions for derivatives and porosity

This class is a decorator for a phase manager, including evaluation of
the pressures and saturations w.r.t. to the primary variables and the porosity.

\author vuong
*/
class PhaseManagerDerivAndPorosity : public PhaseManagerDeriv
  {
  public:
    //! constructor
    explicit PhaseManagerDerivAndPorosity(Teuchos::RCP< POROFLUIDMANAGER::PhaseManagerInterface > phasemanager);

    //! destructor
    virtual ~PhaseManagerDerivAndPorosity(){};

    //! evaluate pressures, saturations and derivatives at GP (matnum is the material number of the porofluid-material on the current element)
    //! default is set to zero, if called from a porofluidmultiphase-element
    //! otherwise it has to be explicitly passed from the caller
    void EvaluateGPState(
        double                 J,
        const VariableManagerMinAccess& varmanager,
        const int matnum = 0);

    //! clear the states
    void ClearGPState();

    //! @name Access methods

    //! get porosity
    double Porosity() const;

    //! get porosity
    double JacobianDefGrad() const;

    //! get derivative of porosity wrt JacobianDefGrad
    double PorosityDerivWrtJacobianDefGrad() const;

    //! get derivative of porosity w.r.t. DOF 'doftoderive'
    double PorosityDeriv(
            int doftoderive) const;

    //! check if porosity depends on fluid (pressure)
    bool PorosityDependsOnFluid() const;

    //! check if porosity depends on structure (basically Jacobian of def. gradient)
    bool PorosityDependsOnStruct() const;

    //@}

  private:

    //! porosity
    double porosity_;

    //! Jacobian of def gradient
    double J_;

    //! derivative of porosity w.r.t. Jacobian of defgradient
    double dporosity_dJ_;

    //! derivative of porosity w.r.t. solid pressure
    double dporosity_dp_;

    //! derivative of porosity w.r.t. degrees of freedom
    Teuchos::RCP<Epetra_SerialDenseVector> porosityderiv_;

  };

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/

/*!
\brief wrapper class, extensions for reaction/mass exchange terms

This class is a decorator for a phase manager, including evaluation of
reaction/mass exchange terms.

\author vuong
*/
class PhaseManagerReaction : public PhaseManagerDecorator
  {
  public:
      //! constructor
    PhaseManagerReaction(Teuchos::RCP< POROFLUIDMANAGER::PhaseManagerInterface > phasemanager);

    //! destructor
    virtual ~PhaseManagerReaction(){};

    //! setup (matnum is the material number of the porofluid-material on the current element)
    //! default is set to zero, if called from a porofluidmultiphase-element
    //! otherwise it has to be explicitly passed from the caller
    void Setup(const DRT::Element* ele, const int matnum = 0);

    //! evaluate pressures, saturations and derivatives at GP (matnum is the material number of the porofluid-material on the current element)
    //! default is set to zero, if called from a porofluidmultiphase-element
    //! otherwise it has to be explicitly passed from the caller
    void EvaluateGPState(
        double                 J,
        const VariableManagerMinAccess& varmanager,
        const int matnum = 0);

    //! clear the states
    void ClearGPState();

    //! @name Access methods

    //! get the reaction term
    double ReacTerm(int phasenum) const;

    //! get the derivative of the reaction term
    double ReacDeriv(int phasenum, int doftoderive) const;

    //! check if the current phase is involved in a reaction
    bool IsReactive(int phasenum) const {CheckIsSetup();return isreactive_[phasenum];};

    //! get scalar to phase ID
    int ScalarToPhaseID(int iscal) const;

    //! get total number of scalars in system
    int NumScal() const;

    //! get the derivative of the reaction term w.r.t. scalar 'scaltoderive'
    double ReacDerivScalar(int phasenum, int scaltoderive) const;

    //! get the derivative of the reaction term w.r.t. porosity
    double ReacDerivPorosity(int phasenum) const;

    //@}

  private:
    //! reaction terms
    std::vector<double> reacterms_;
    //! derivatives of reaction terms w.r.t. fluid primary dofs
    std::vector<std::vector<double> > reactermsderivs_;
    //! derivatives of reaction terms w.r.t. (true) pressures
    std::vector<std::vector<double> > reactermsderivspressure_;
    //! derivatives of reaction terms w.r.t. saturations
    std::vector<std::vector<double> > reactermsderivssaturation_;
    //! derivatives of reaction terms w.r.t. porosity
    std::vector<double> reactermsderivsporosity_;
    //! derivatives of reaction terms w.r.t. scalars --> needed for off-diagonal matrices
    std::vector<std::vector<double> > reactermsderivsscalar_;
    //! derivatives of reaction terms w.r.t. volume fraction
    std::vector<std::vector<double> > reactermsderivsvolfrac_;

    //! flags indicating whether the phase is involved in a reaction
    std::vector<bool> isreactive_;
    //! scalar to phase ids
    std::vector<int> scalartophaseid_;
    //! number of scalars
    int numscal_;

  };

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/

/*!
\brief wrapper class, extensions for diffusion

This class is a decorator for a phase manager, including evaluation of
the diffusion tensor (inverse permeability). As the tensor is saved
as a LINALG::Matrix, it is templated by the number of space dimensions

\author vuong
*/
template<int nsd>
class PhaseManagerDiffusion : public PhaseManagerDecorator
  {
  public:
     //! constructor
    PhaseManagerDiffusion(Teuchos::RCP< POROFLUIDMANAGER::PhaseManagerInterface > phasemanager);

    //! destructor
    virtual ~PhaseManagerDiffusion(){};

    //! setup (matnum is the material number of the porofluid-material on the current element)
    //! default is set to zero, if called from a porofluidmultiphase-element
    //! otherwise it has to be explicitly passed from the caller
    void Setup(const DRT::Element* ele, const int matnum = 0);

    //! evaluate pressures, saturations and derivatives at GP (matnum is the material number of the porofluid-material on the current element)
    //! default is set to zero, if called from a porofluidmultiphase-element
    //! otherwise it has to be explicitly passed from the caller
    void EvaluateGPState(
        double                 J,
        const VariableManagerMinAccess& varmanager,
        const int matnum = 0);

    //! clear the states
    void ClearGPState();

    //! @name Access methods

    //! get the diffusion tensor
    void PermeabilityTensor(int phasenum, LINALG::Matrix<nsd,nsd>& permeabilitytensor) const;

    //! check for constant relpermeability
    bool HasConstantRelPermeability(int phasenum) const;

    //! get relative diffusivity of phase
    double RelPermeability(int phasenum) const;

    //! get derivative of relative permeability of phase
    double RelPermeabilityDeriv(int phasenum) const;

    //! check for constant dynamic viscosity
    bool HasConstantDynViscosity(int phasenum) const;
    //! get dynamic viscosity of phase
    double DynViscosity(int phasenum, double abspressgrad, int matnum = 0) const;
    //! get dynamic viscosity of phase
    double DynViscosity(const MAT::Material& material, int phasenum, double abspressgrad) const;
    //! get derivative of dynamic viscosity of phase
    double DynViscosityDeriv(int phasenum, double abspressgrad) const;
    //! get derivative dynamic viscosity of phase
    double DynViscosityDeriv(const MAT::Material& material, int phasenum, double abspressgrad) const;

    //@}

  private:
    //! diffusion tensor
    std::vector<LINALG::Matrix<nsd,nsd> > permeabilitytensors_;
    //! relative diffusivities
    std::vector<double> relpermeabilities_;
    //! derivative of relative permeabilities w.r.t. saturation
    std::vector<double> derrelpermeabilities_;
    //! check for constant relative permeabilities of phase
    std::vector<bool> constrelpermeability_;
    //! check for constant dynamic viscosities of phase
    std::vector<bool> constdynviscosity_;

  };

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/

/*!
\brief wrapper class, extensions for additional volume fractions

This class is a decorator for a phase manager, including evaluation of
the additional volume fractions, e.g. for endothelial cells

\author vuong
*/
template<int nsd>
class PhaseManagerVolFrac : public PhaseManagerDecorator
  {
  public:
     //! constructor
  PhaseManagerVolFrac(Teuchos::RCP< POROFLUIDMANAGER::PhaseManagerInterface > phasemanager);

    //! destructor
    virtual ~PhaseManagerVolFrac(){};

    //! setup (matnum is the material number of the porofluid-material on the current element)
    //! default is set to zero, if called from a porofluidmultiphase-element
    //! otherwise it has to be explicitly passed from the caller
    void Setup(const DRT::Element* ele, const int matnum = 0);

    //! evaluate pressures, saturations and derivatives at GP (matnum is the material number of the porofluid-material on the current element)
    //! default is set to zero, if called from a porofluidmultiphase-element
    //! otherwise it has to be explicitly passed from the caller
    void EvaluateGPState(
        double                 J,
        const VariableManagerMinAccess& varmanager,
        const int matnum = 0);

    //! clear the states
    void ClearGPState();

    //! @name Access methods

    //! get the diffusion tensor
    void DiffTensorVolFrac(int volfracnum, LINALG::Matrix<nsd,nsd>& difftensorvolfrac) const;

    //! get density of volume frac 'volfracnum'
    virtual double VolFracDensity(int volfracnum) const;

    //! get pressure of volume frac 'volfracnum'
    virtual double VolFracPressure(int volfracnum) const;

    //! check if volume frac 'volfracnum' has additional scalar dependent flux
    virtual bool HasAddScalarDependentFlux(int volfracnum) const;

    //! return scalar diffusivity of scalar 'iscal' of volume fraction 'volfracnum'
    virtual double ScalarDiff(int volfracnum, int iscal) const;

    //! return omega half of scalar 'iscal' of volume fraction 'volfracnum' for receptor kinetic law
    //! see: Anderson, A. R. A. & Chaplain, M. A. J.
    //       Continuous and discrete mathematical models of tumor-induced angiogenesis
    virtual double OmegaHalf(int volfracnum, int iscal) const;

    //! check if volume frac 'volfracnum' has additional scalar dependent flux of scalar 'iscal'
    virtual bool HasAddScalarDependentFlux(int volfracnum, int iscal) const;

    //! check if volume frac 'volfracnum' has receptor kinetic-law of scalar 'iscal'
    //! see: Anderson, A. R. A. & Chaplain, M. A. J.
    //       Continuous and discrete mathematical models of tumor-induced angiogenesis
    virtual bool HasReceptorKineticLaw(int volfracnum, int iscal) const;

    //@}

  private:
    //! diffusion tensors
    std::vector<LINALG::Matrix<nsd,nsd> > difftensorsvolfrac_;
    //! densities of volume fractions
    std::vector<double> volfracdensity_;
    //! pressure of volume fractions
    std::vector<double> volfracpressure_;
    //! does the material have additional scalar dependent flux
    std::vector<bool> hasaddscalardpendentflux_;
    //! matrix for scalar diffusivities
    std::vector < std::vector<double> > scalardiffs_;
    //! matrix for omega half values of receptor-kinetic law
    std::vector < std::vector<double> > omega_half_;

  };


/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
} // namespace POROFLUIDMANAGER
} // namespace ELEMENT
} // namespace DRT


#endif /* SRC_DRT_POROFLUIDMULTPHASE_ELE_POROFLUID_PHASEMANAGER_H_ */
