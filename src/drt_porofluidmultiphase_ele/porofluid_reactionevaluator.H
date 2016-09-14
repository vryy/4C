/*----------------------------------------------------------------------*/
/*!
 \file porofluid_reactionevaluator.H

 \brief evaluation class for reactive terms (mass sources and sinks)

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUID_REACTIONEVALUATOR_H_
#define SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUID_REACTIONEVALUATOR_H_

#include <vector>

namespace MAT
{
  class Material;
}

namespace DRT
{
namespace ELEMENTS
{
class PoroFluidPhaseManager;

namespace POROFLUIDEVALUATOR
{
class ReactionEvaluator
  {
  public:
    //! constructor
    ReactionEvaluator();

    //! destructor
    virtual ~ReactionEvaluator(){};

    //! setup
    void Setup(
        const PoroFluidPhaseManager&       phasemanager,
        const MAT::Material& material);

    //! evaluate reaction terms at GP
    void EvaluateGPState(
        const PoroFluidPhaseManager&       phasemanager,      //!< phase manager
        const MAT::Material&               material,
        const double&                      porosity,
        const std::vector<double> &        scalar
        );

    //! clear the states
    void ClearGPState();

    //! get the reaction term
    double GetReacTerm(int phasenum) const;

    //! get the derivative of the reaction term
    double GetReacDeriv(int phasenum, int doftoderive) const;

    //! check if the current phase is involved in a reaction
    bool IsReactive(int phasenum) const {return isreactive_[phasenum];};

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

    //! flag indicating gauss point was evaluated
    bool isevaluated_;
  };
}
}
}


#endif /* SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUID_REACTIONEVALUATOR_H_ */
