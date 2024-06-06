/*---------------------------------------------------------------------*/
/*! \file
\brief Contact parameter interface. Necessary for the communication
       between the structural time integration framework and the
       contact strategies.

\level 3


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_CONTACT_PARAMSINTERFACE_HPP
#define FOUR_C_CONTACT_PARAMSINTERFACE_HPP

#include "4C_config.hpp"

#include "4C_mortar_paramsinterface.hpp"

FOUR_C_NAMESPACE_OPEN

namespace NOX
{
  namespace Nln
  {
    enum class CorrectionType : int;
  }  // namespace Nln
}  // namespace NOX
namespace Inpar
{
  namespace STR
  {
    enum PredEnum : int;
  }  // namespace STR
  namespace CONTACT
  {
    enum VariationalApproach : int;
    enum class CouplingScheme : int;
  }  // namespace CONTACT
}  // namespace Inpar

namespace STR
{
  namespace MODELEVALUATOR
  {
    class Generic;
  }  // namespace MODELEVALUATOR
}  // namespace STR

namespace CONTACT
{
  class ParamsInterface : public Mortar::ParamsInterface
  {
   public:
    //! get the predictor status
    virtual bool is_predictor() const = 0;

    //! Is the current state coming from a predictor step?
    virtual bool is_predictor_state() const = 0;

    //! \brief get the currently active predictor type
    /** \note If the execution of the predictor is finished, this
     *  function will return Inpar::STR::pred_vague. \author hiermeier */
    virtual enum Inpar::STR::PredEnum get_predictor_type() const = 0;

    //! get the current step length
    virtual double get_step_length() const = 0;

    //! get number of linear system corrections (modified Newton approach)
    virtual int get_number_of_modified_newton_corrections() const = 0;

    //! get the is_default_step indicator
    virtual bool is_default_step() const = 0;

    //! get correction type
    virtual NOX::Nln::CorrectionType get_correction_type() const = 0;

    //! get the current time step [derived]
    virtual double get_delta_time() const = 0;

    //! get a pointer to the contact model evaluator
    virtual const STR::MODELEVALUATOR::Generic& get_model_evaluator() const = 0;

    //! get the output file path
    virtual std::string get_output_file_path() const = 0;

    //! get the variational approach type
    virtual enum Inpar::CONTACT::VariationalApproach get_variational_approach_type() const = 0;

    //! set the variational approach type
    virtual void set_variational_approach_type(
        const enum Inpar::CONTACT::VariationalApproach var_type) = 0;

    //! set the coupling approach mode
    virtual enum Inpar::CONTACT::CouplingScheme get_coupling_scheme() const = 0;

    //! set the coupling scheme
    virtual void set_coupling_scheme(const enum Inpar::CONTACT::CouplingScheme scheme) = 0;
  };
}  // namespace CONTACT


FOUR_C_NAMESPACE_CLOSE

#endif
