// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_STR_PASIWRAPPER_HPP
#define FOUR_C_ADAPTER_STR_PASIWRAPPER_HPP

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_adapter_str_wrapper.hpp"
#include "4C_pasi_str_model_evaluator_partitioned.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | forward declarations                                                 |
 *----------------------------------------------------------------------*/
namespace Solid
{
  class MapExtractor;

  namespace ModelEvaluator
  {
    class PartitionedPASI;
  }  // namespace ModelEvaluator
}  // namespace Solid

/*----------------------------------------------------------------------*
 | class declarations                                                   |
 *----------------------------------------------------------------------*/
namespace Adapter
{
  class PASIStructureWrapper : public StructureWrapper
  {
   public:
    //! constructor
    explicit PASIStructureWrapper(std::shared_ptr<Structure> structure);

    //! set pointer to model evaluator
    void set_model_evaluator_ptr(std::shared_ptr<Solid::ModelEvaluator::PartitionedPASI> me)
    {
      pasi_model_evaluator_ = me;
      return;
    }

    //! communication object at the interface
    virtual std::shared_ptr<const Solid::MapExtractor> interface() const { return interface_; }

    //! apply interface force to structure interface
    void apply_interface_force(std::shared_ptr<const Core::LinAlg::Vector<double>> intfforce);

   protected:
    //! The structural model evaluator object.
    //! Your PASI algorithm calls methods in this adapter.
    //! If this method is related to the structural field,
    //! a corresponding method in the model evaluator may be
    //! called, if necessary.
    std::shared_ptr<Solid::ModelEvaluator::PartitionedPASI> pasi_model_evaluator_;

    //! access the pasi model evaluator
    std::shared_ptr<Solid::ModelEvaluator::PartitionedPASI> pasi_model_evaluator()
    {
      return pasi_model_evaluator_;
    };

    //! the interface map setup for interface <-> full translation
    std::shared_ptr<Solid::MapExtractor> interface_;
  };

}  // namespace Adapter

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
