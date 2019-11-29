/*----------------------------------------------------------------------*/
/*! \file
\file  ad_fbi_constraintbridge.H

\brief Abstract class to be overloaded by different adapter implementations connecting a constraint
enforcement technique with a discretization approach for Fluid-beam interaction.

\level 3

\maintainer Nora Hagmeyer
*----------------------------------------------------------------------*/

#ifndef AD_FBI_CONSTRAINTBRIDGE_H_
#define AD_FBI_CONSTRAINTBRIDGE_H_

#include <Teuchos_RCP.hpp>
#include "../linalg/linalg_utils.H"
#include <Epetra_Map.h>

// Forward decleration

namespace DRT
{
  class Discretization;
  class Element;
}  // namespace DRT
namespace BEAMINTERACTION
{
  class BeamContactPair;
}
namespace GEOMETRYPAIR
{
  class LineTo3DEvaluationData;
}
namespace FBI
{
  class BeamToFluidMeshtyingParams;
}
namespace ADAPTER
{
  class ConstraintEnforcerFactory;
  class FBIConstraintBridge
  {
    friend ConstraintEnforcerFactory;

   public:
    /// empty destructor
    virtual ~FBIConstraintBridge(){};

    /**
     * \brief Initializes all members of the class     *
     */
    virtual void Setup(const Epetra_Map* beam_map, const Epetra_Map* fluid_map);
    /**
     * \brief Computes the coupling matrices
     *
     * This is where the magic happens. The global meshtying contributions are integrated using
     * information of the beam elements, the fluid elements and their position relative to each
     * other.
     *
     */
    virtual void Evaluate(
        Teuchos::RCP<const std::vector<Teuchos::RCP<DRT::Discretization>>> discretizations) = 0;

    /**
     * \brief Wraps the ResetState function of the pair
     *
     * Here, the current values of the single fields are handed to/updated in the pair
     */
    virtual void ResetPair(Teuchos::RCP<const std::vector<double>> beam_centerline_dofvec,
        const Teuchos::RCP<const std::vector<double>> fluid_nodal_dofvec,
        Teuchos::RCP<BEAMINTERACTION::BeamContactPair> interactionpair);

    /// Creates a volume_beam_meshtying pair
    virtual void CreatePair(std::vector<DRT::Element const*> elements,
        Teuchos::RCP<std::vector<double>> beam_centerline_dofvec,
        Teuchos::RCP<std::vector<double>> fluid_nodal_dofvec);

    virtual const Teuchos::RCP<std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>>>
    GetPairs() const final
    {
      return meshtying_pairs_;
    };
    /// returns data container holding all beam interaction related parameters
    virtual Teuchos::RCP<FBI::BeamToFluidMeshtyingParams> GetParams() const final
    {
      return beam_interaction_params_;
    };

    /// returns data container holding all geometry related evaluation data
    virtual const Teuchos::RCP<GEOMETRYPAIR::LineTo3DEvaluationData> GetGeometryData() const final
    {
      return geometry_evaluation_data_;
    };

    /// Clears the pair vector and builds the pairs anew
    virtual void ResetBridge();

   protected:
    /// Standard constructor
    FBIConstraintBridge();

    /// Provides access to data container holding all beam interaction related parameters
    virtual void SetParams(Teuchos::RCP<FBI::BeamToFluidMeshtyingParams> params) final
    {
      beam_interaction_params_ = params;
    };

    /// Provides access to data container holding all geometry related evaluation data
    virtual void SetGeometryData(Teuchos::RCP<GEOMETRYPAIR::LineTo3DEvaluationData> data) final
    {
      geometry_evaluation_data_ = data;
    };

    Teuchos::RCP<std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>>>
        meshtying_pairs_;  // todo make this private?
   private:
    /// data container holding all beam interaction related parameters
    Teuchos::RCP<FBI::BeamToFluidMeshtyingParams> beam_interaction_params_;

    /// data container holding all geometry related evaluation data
    Teuchos::RCP<GEOMETRYPAIR::LineTo3DEvaluationData> geometry_evaluation_data_;
  };
}  // namespace ADAPTER

#endif /* AD_FBI_CONSTRAINTBRIDGE_H_ */
