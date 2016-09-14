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

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_RCP.hpp"
#include <vector>

//forward declarations
namespace MAT
{
  class Material;
  class FluidPoroSinglePhase;
  class FluidPoroMultiPhase;
}

namespace DRT
{
namespace ELEMENTS
{

namespace POROFLUIDEVALUATOR
{
class ReactionEvaluator;
}


class PoroFluidPhaseManager
  {
  public:
    //! constructor
    explicit PoroFluidPhaseManager(int numphase);

    //! destructor
    virtual ~PoroFluidPhaseManager(){};

    //! setup
    void Setup(const MAT::Material& material);

    //! evaluate pressures, saturations and derivatives at GP
    void EvaluateGPState(
        const MAT::Material& material,
        const std::vector<double>& phinp);

    //! evaluate reactionss at GP
    void EvaluateGPReactions(
        const MAT::Material& material,
        double porosity,
        const std::vector<double>& scalars);

    //! clear the states
    void ClearGPState();

    //! check for reactions
    bool IsReactive(int phasenum);

    //! @name Access methods

    //! get the number of phases
    int NumPhases() const {return numphases_;};

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
    const std::vector<double>& Saturation() const {return saturation_;};

    //! get pressure of all phases
    const std::vector<double>& Pressure() const {return pressure_;};

    //! get bulk modulus of phase 'phasenum'
    double Bulkmodulus(
        const MAT::Material& material,
        int phasenum) const;

    //! get density of phase 'phasenum'
    double Density(
        const MAT::Material& material,
        int phasenum) const;

    //! get the reaction term
    double ReacTerm(int phasenum) const;

    //! get the derivative of the reaction term
    double ReacDeriv(int phasenum, int doftoderive) const;

    //@}

  private:
    //! transform generalized pressures to true pressures
    void TransformGenPresToTruePres(
        const std::vector<double>& phi,
        std::vector<double>& phi_transformed);

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

    //! flag indicating of gauss point state has been set and evaluated
    bool statesset_;

    //! helper class for evaluation of reaction terms
    Teuchos::RCP< POROFLUIDEVALUATOR::ReactionEvaluator > reactionevaluator_;
  };

}
}


#endif /* SRC_DRT_POROFLUIDMULTPHASE_ELE_POROFLUID_PHASEMANAGER_H_ */
