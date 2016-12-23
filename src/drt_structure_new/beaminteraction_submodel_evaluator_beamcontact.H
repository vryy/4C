/*-----------------------------------------------------------*/
/*!
\file beaminteraction_submodel_evaluator_beamcontact.H

\brief class for submodel beam contact

\maintainer Jonas Eichinger, Maximilian Grill

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef BEAMINTERACTION_SUBMODEL_EVALUATOR_BEAMCONTACT_H_
#define BEAMINTERACTION_SUBMODEL_EVALUATOR_BEAMCONTACT_H_

#include <Teuchos_RCP.hpp>
#include "beaminteraction_submodel_evaluator_generic.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_inpar/inpar_binstrategy.H"
#include "../drt_inpar/inpar_beaminteraction.H"

//forward declaration ...
class Epetra_Vector;
class Epetra_Map;
namespace DRT
{
  class Exporter;
  class Element;
  class Node;
  class PackBuffer;
}
namespace LINALG {
  class SparseMatrix;
  class SerialDenseVector;
  class SerialDenseMatrix;
}
namespace BEAMINTERACTION {
class BeamContactParams;
class BeamToBeamInteraction;

namespace SUBMODELEVALUATOR {

  class BeamContact : public Generic
  {
  public:
    //! constructor
    BeamContact();

    //! destructor
    virtual ~BeamContact() {};

    //! setup class variables
    virtual void Setup();

    //! Returns the type of the current model evaluator
    INPAR::BEAMINTERACTION::SubModelType Type() const
    { return INPAR::BEAMINTERACTION::submodel_beamcontact; }

    //! @name Derived public STR::MODELEVALUATOR::Generic methods
    //! @{
    //! \brief reset model specific variables (without jacobian)
    virtual void Reset();

    //! \brief Evaluate the current right-hand-side at \f$t_{n+1}\f$
    virtual bool EvaluateForce();

    //! \brief Evaluate the current tangential stiffness matrix at \f$t_{n+1}\f$
    virtual bool EvaluateStiff();

    //! \brief Evaluate the current right-hand-side vector and tangential stiffness matrix at \f$t_{n+1}\f$
    virtual bool EvaluateForceStiff();

    virtual void UpdateStepState(const double& timefac_n);

    virtual void PreUpdateStepElement();

    virtual void UpdateStepElement();

    virtual void PostUpdateStepElement();

    //! output routine for model evlaluator
    virtual void OutputStepState(IO::DiscretizationWriter& iowriter) const;

    //! reset routine for model evlaluator
    virtual void ResetStepState();

    //! @}


    //!@name routines that are not derived and handle beam to beam contact
    //! @{
    /// get neighbouring eles in discret
    virtual void FindAndStoreNeighboringElements();

    /// exclude certein neighbors from interaction evaluation
    virtual void SelectElesToBeConsideredForContactEvaluation(
        DRT::Element* currele,
        std::set<DRT::Element*>& neighbors
    ) const;

    /// create instances of class BeamToBeamInteraction that will be evaluated
    //  to get force and stiffness contributions from beam interactions
    virtual void CreateBeamToBeamContactElementPairs();

    /// print
    void PrintAllBeamToBeamContactElementPairs(std::ostream& out) const;

    /// print
    void PrintActiveBeamToBeamContactSet(std::ostream& out) const;
    //! @}

  private:
    //! data container holding all beam contact related parameters
    Teuchos::RCP<BEAMINTERACTION::BeamContactParams> beam_contact_params_ptr_;

    //! type of eles in bins
    INPAR::BINSTRATEGY::BinContent bin_beamcontent_;

     //! interacting pairs of beam elements that might exert forces on each other
    std::vector<Teuchos::RCP<BEAMINTERACTION::BeamToBeamInteraction> > BTB_contact_elepairs_;

    //! mapping beam ele (elegid) to set of spatially proximal eles (pointer to elements)
    std::map<int, std::set<DRT::Element*> > nearby_elements_map_;

  };

}
}

#endif
