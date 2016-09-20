/*----------------------------------------------------------------------*/
/*!
 \file porofluid_phasemanager.H

 \brief manager class for handling the phases and their dofs on element level

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROFLUIDMULTPHASE_ELE_POROFLUID_PHASEMANAGER_H_
#define SRC_DRT_POROFLUIDMULTPHASE_ELE_POROFLUID_PHASEMANAGER_H_

#include "porofluidmultiphase_ele_action.H"

#include "../drt_inpar/inpar_material.H"
#include "../drt_lib/drt_dserror.H"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RCP.hpp"
#include <vector>

//forward declarations
namespace MAT
{
  class Material;
  class StructPoro;
  class FluidPoroSinglePhase;
  class FluidPoroMultiPhase;
}

namespace DRT
{
namespace ELEMENTS
{
class PoroFluidMultiPhaseEleParameter;

namespace POROFLUIDMANAGER
{
class VariableManagerMinAccess;

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/

//! interface to phase manager classes
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
        INPAR::MAT::MaterialType mattype,
        const POROFLUIDMULTIPHASE::Action&     action,
        int numphases
        );

    //! factory method
    static Teuchos::RCP< DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerInterface > WrapPhaseManager(
        const DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter& para,
        INPAR::MAT::MaterialType mattype,
        const POROFLUIDMULTIPHASE::Action&     action,
        Teuchos::RCP< PhaseManagerInterface > corephasemanager
        );

    //! setup
    virtual void Setup(const MAT::Material& material)=0;

    //! evaluate pressures, saturations and derivatives at GP
    virtual void EvaluateGPState(
        const MAT::Material&   material,
        MAT::StructPoro&       matstructporo,
        double                 J,
        const VariableManagerMinAccess& varmanager)=0;

    //! clear the states
    virtual void ClearGPState()=0;

    //! check for reactions
    virtual bool IsReactive(int phasenum) const=0;

    //! check if EvaluateGPState() was called
    virtual void CheckIsEvaluated() const=0;

    //! unwrap the manager
    virtual Teuchos::RCP< DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerInterface >
    UnWrap()=0;

    //! @name Access methods

    //! get the number of phases
    virtual int NumPhases() const=0;

    //! evaluate relative diffusivity of phase 'phasenum'
    virtual double EvaluateRelDiffusivity(
        const MAT::Material& material,
        int phasenum) const=0;

    //! get solid pressure
    virtual double SolidPressure() const=0;

    //! get saturation of phase 'phasenum'
    virtual double Saturation(int phasenum) const=0;

    //! get pressure of phase 'phasenum'
    virtual double Pressure(int phasenum) const=0;

    //! get saturation of all phases
    virtual const std::vector<double>& Saturation() const =0;

    //! get pressure of all phases
    virtual const std::vector<double>& Pressure() const =0;

    //! get porosity
    virtual double Porosity() const =0;

    //! get bulk modulus of phase 'phasenum'
    virtual double Bulkmodulus(
        const MAT::Material& material,
        int phasenum) const=0;

    //! get density of phase 'phasenum'
    virtual double Density(
        const MAT::Material& material,
        int phasenum) const=0;

    //! get derivative of saturation of phase 'phasenum' w.r.t. DOF 'doftoderive'
    virtual double SaturationDeriv(
            int phasenum,
            int doftoderive) const
    {dserror("Derivative of saturation not available for this phase manager!");return 0.0;};

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

    //! get the derivative of the reaction term
    virtual double ReacDeriv(int phasenum, int doftoderive) const
    {dserror("Reaction term not available for this phase manager!");return 0.0;};

    //@}

  };

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/

//! standard phase manager, holding pressures and saturations
class PhaseManagerCore : public PhaseManagerInterface
  {
  public:
    //! constructor
    explicit PhaseManagerCore(int numphase);

    //! copy constructor
    PhaseManagerCore(const PhaseManagerCore& old);

    //! destructor
    virtual ~PhaseManagerCore(){};

    //! setup
    void Setup(const MAT::Material& material);

    //! evaluate pressures, saturations and derivatives at GP
    void EvaluateGPState(
        const MAT::Material&   material,
        MAT::StructPoro&       matstructporo,
        double                 J,
        const VariableManagerMinAccess& varmanager);

    //! clear the states
    void ClearGPState();

    //! check for reactions
    bool IsReactive(int phasenum) const {return false;};

    //! unwrap the manager
    virtual Teuchos::RCP< DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerInterface >
    UnWrap() {return Teuchos::rcp(new PhaseManagerCore(*this));};

    //! @name Access methods

    //! get the number of phases
    int NumPhases() const {return numphases_;};

    //! evaluate relative diffusivity of phase 'phasenum'
    double EvaluateRelDiffusivity(
        const MAT::Material& material,
        int phasenum) const;

    //! get solid pressure
    double SolidPressure() const;

    //! get saturation of phase 'phasenum'
    double Saturation(int phasenum) const;

    //! get pressure of phase 'phasenum'
    double Pressure(int phasenum) const;

    //! get saturation of all phases
    const std::vector<double>& Saturation() const;

    //! get pressure of all phases
    const std::vector<double>& Pressure() const;

    //! get porosity
    double Porosity() const;

    //! get bulk modulus of phase 'phasenum'
    double Bulkmodulus(
        const MAT::Material& material,
        int phasenum) const;

    //! get density of phase 'phasenum'
    double Density(
        const MAT::Material& material,
        int phasenum) const;

    //@}

    //! check if EvaluateGPState() was called
    void CheckIsEvaluated() const
    {
      if(not isevaluated_)
        dserror("Gauss point states have not been set!");
    }

  protected:
    //! transform generalized pressures to true pressures
    void TransformGenPresToTruePres(
        const std::vector<double>& phi,
        std::vector<double>& phi_transformed);

    //! compute porosity
    double ComputePorosity(
        MAT::StructPoro &  structmat,    //!< material of the structure
        double             J,             //!< determinant of deformation gradient
        double pres
      );

    //! transformation of degrees of freedom to true pressures
    Teuchos::RCP<Epetra_SerialDenseMatrix> dof2pres_;

    //! number of phases
    const int numphases_;

    //! generalized pressure
    std::vector<double> genpressure_;
    //! true pressure
    std::vector<double> pressure_;
    //! saturation
    std::vector<double> saturation_;
    //! solid pressure
    double solidpressure_;
    //! porosity
    double porosity_;

    //! flag indicating of gauss point state has been set and evaluated
    bool isevaluated_;

    //! flag of Setup was called
    bool issetup_;

  private:
    //! no standard constructor
    PhaseManagerCore();
  };

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/

//! wrapper class, extensions for derivatives
class PhaseManagerDeriv : public PhaseManagerInterface
  {
  public:
    //! constructor
    explicit PhaseManagerDeriv(Teuchos::RCP< POROFLUIDMANAGER::PhaseManagerInterface > phasemanager);

    //! destructor
    virtual ~PhaseManagerDeriv(){};

    //! setup
    void Setup(const MAT::Material& material)
    {phasemanager_->Setup(material);};

    //! evaluate pressures, saturations and derivatives at GP
    void EvaluateGPState(
        const MAT::Material&   material,
        MAT::StructPoro&       matstructporo,
        double                 J,
        const VariableManagerMinAccess& varmanager);

    //! clear the states
    void ClearGPState();

    //! check if EvaluateGPState() was called
    void CheckIsEvaluated() const
    {phasemanager_->CheckIsEvaluated();};

    //! unwrap the manager
    virtual Teuchos::RCP< DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerInterface >
    UnWrap() {return phasemanager_->UnWrap();};

    //! @name Access methods

    //! get derivative of saturation of phase 'phasenum' w.r.t. DOF 'doftoderive'
    double SaturationDeriv(
            int phasenum,
            int doftoderive) const;

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

    //! check if the current phase is involved in a reaction
    bool IsReactive(int phasenum) const
    {return phasemanager_->IsReactive(phasenum);};

    //! get the number of phases
    int NumPhases() const
    {return phasemanager_->NumPhases();};

    //! evaluate relative diffusivity of phase 'phasenum'
    double EvaluateRelDiffusivity(
        const MAT::Material& material,
        int phasenum) const
    {return phasemanager_->EvaluateRelDiffusivity(material,phasenum);};

    //! get solid pressure
    double SolidPressure() const
    {return phasemanager_->SolidPressure();};

    //! get saturation of phase 'phasenum'
    double Saturation(int phasenum) const
    {return phasemanager_->Saturation(phasenum);};

    //! get pressure of phase 'phasenum'
    double Pressure(int phasenum) const
    {return phasemanager_->Pressure(phasenum);};

    //! get saturation of all phases
    const std::vector<double>& Saturation() const
    {return phasemanager_->Saturation();};

    //! get pressure of all phases
    const std::vector<double>& Pressure() const
    {return phasemanager_->Pressure();};

    //! get bulk modulus of phase 'phasenum'
    double Bulkmodulus(
        const MAT::Material& material,
        int phasenum) const
    {return phasemanager_->Bulkmodulus(material,phasenum);};

    //! get porosity
    double Porosity() const
    {return phasemanager_->Porosity();};

    //! get density of phase 'phasenum'
    double Density(
        const MAT::Material& material,
        int phasenum) const
    {return phasemanager_->Density(material,phasenum);};

    //@}

  private:

    //! evaluate derivative of saturation with respect to pressure
    double EvaluateDerivOfSaturationWrtPressure(
        const MAT::Material& material,
        int phasenum,
        int doftoderive,
        const std::vector<double>& state) const;

    //! evaluate derivative of degree of freedom with respect to pressure
    double EvaluateDerivOfDofWrtPressure(
        const MAT::Material& material,
        int phasenum,
        int doftoderive,
        const std::vector<double>& state) const;

    //! wrapped phase manager
    Teuchos::RCP< POROFLUIDMANAGER::PhaseManagerInterface > phasemanager_;

    //! derivative of true pressure w.r.t. degrees of freedom
    // first index: pressure, second index: dof
    Teuchos::RCP<Epetra_SerialDenseMatrix> pressurederiv_;
    //! derivative of saturations w.r.t. degrees of freedom
    // first index: saturation, second index: dof
    Teuchos::RCP<Epetra_SerialDenseMatrix> saturationderiv_;
    //! derivative of solid pressure w.r.t. degrees of freedom
    Teuchos::RCP<Epetra_SerialDenseVector> solidpressurederiv_;
    //! second derivative of solid pressure w.r.t. degrees of freedom;
    Teuchos::RCP<Epetra_SerialDenseMatrix> solidpressurederivderiv_;

    //! no standard constructor
    PhaseManagerDeriv();

  };

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/

//! wrapper class, extensions for reaction/mass exchange terms
class PhaseManagerReaction : public PhaseManagerInterface
  {
  public:
      //! constructor
    PhaseManagerReaction(Teuchos::RCP< POROFLUIDMANAGER::PhaseManagerInterface > phasemanager);

    //! destructor
    virtual ~PhaseManagerReaction(){};

    //! setup
    void Setup(const MAT::Material& material);

    //! evaluate pressures, saturations and derivatives at GP
    void EvaluateGPState(
        const MAT::Material&   material,
        MAT::StructPoro&       matstructporo,
        double                 J,
        const VariableManagerMinAccess& varmanager);

    //! clear the states
    void ClearGPState();

    //! unwrap the manager
    virtual Teuchos::RCP< DRT::ELEMENTS::POROFLUIDMANAGER::PhaseManagerInterface >
    UnWrap() {return phasemanager_->UnWrap();};

    //! @name Access methods

    //! get the reaction term
    double ReacTerm(int phasenum) const;

    //! get the derivative of the reaction term
    double ReacDeriv(int phasenum, int doftoderive) const;

    //! check if the current phase is involved in a reaction
    bool IsReactive(int phasenum) const {CheckIsEvaluated();return isreactive_[phasenum];};

    //! get the number of phases
    int NumPhases() const
    {return phasemanager_->NumPhases();};

    //! get derivative of saturation of phase 'phasenum' w.r.t. DOF 'doftoderive'
    double SaturationDeriv(
            int phasenum,
            int doftoderive) const
    {return phasemanager_->SaturationDeriv(phasenum,doftoderive);};

    //! get derivative of pressure of phase 'phasenum' w.r.t. DOF 'doftoderive'
    double PressureDeriv(
            int phasenum,
            int doftoderive) const
    {return phasemanager_->PressureDeriv(phasenum,doftoderive);};

    //! get derivative of solid pressure  w.r.t. DOF 'doftoderive'
    double SolidPressureDeriv(
            int doftoderive) const
    {return phasemanager_->SolidPressureDeriv(doftoderive);};

    //! get derivative of pressure of phase 'phasenum'
    //! w.r.t. DOF 'doftoderive' (first derivative)
    //! and w.r.t. DOF 'doftoderive2' (second derivative)
    double SolidPressureDerivDeriv(
            int doftoderive,
            int doftoderive2) const
    {return phasemanager_->SolidPressureDerivDeriv(doftoderive,doftoderive2);};

    //! evaluate relative diffusivity of phase 'phasenum'
    double EvaluateRelDiffusivity(
        const MAT::Material& material,
        int phasenum) const
    {return phasemanager_->EvaluateRelDiffusivity(material,phasenum);};

    //! get solid pressure
    double SolidPressure() const
    {return phasemanager_->SolidPressure();};

    //! get saturation of phase 'phasenum'
    double Saturation(int phasenum) const
    {return phasemanager_->Saturation(phasenum);};

    //! get pressure of phase 'phasenum'
    double Pressure(int phasenum) const
    {return phasemanager_->Pressure(phasenum);};

    //! get saturation of all phases
    const std::vector<double>& Saturation() const
    {return phasemanager_->Saturation();};

    //! get pressure of all phases
    const std::vector<double>& Pressure() const
    {return phasemanager_->Pressure();};

    //! get bulk modulus of phase 'phasenum'
    double Bulkmodulus(
        const MAT::Material& material,
        int phasenum) const
    {return phasemanager_->Bulkmodulus(material,phasenum);};

    //! get porosity
    double Porosity() const
    {return phasemanager_->Porosity();};

    //! get density of phase 'phasenum'
    double Density(
        const MAT::Material& material,
        int phasenum) const
    {return phasemanager_->Density(material,phasenum);};

    //@}

    //! check if EvaluateGPState() was called
    void CheckIsEvaluated() const
    {phasemanager_->CheckIsEvaluated();};

  private:
    //! reaction terms
    std::vector<double> reacterms_;
    //! derivatives of reaction terms w.r.t. dofs
    std::vector<std::vector<double> > reactermsderivs_;
    //! derivatives of reaction terms w.r.t. (true) pressures
    std::vector<std::vector<double> > reactermsderivspressure_;
    //! derivatives of reaction terms w.r.t. saturations
    std::vector<std::vector<double> > reactermsderivssaturation_;
    //! derivatives of reaction terms w.r.t. porosity
    std::vector<double> reactermsderivsporosity_;

    //! flags indicating whether the phase is involved in a reaction
    std::vector<bool> isreactive_;

    //! wrapped phase manager
    Teuchos::RCP< POROFLUIDMANAGER::PhaseManagerInterface > phasemanager_;

    //! no standard constructor
    PhaseManagerReaction();
  };
}
}
}


#endif /* SRC_DRT_POROFLUIDMULTPHASE_ELE_POROFLUID_PHASEMANAGER_H_ */
