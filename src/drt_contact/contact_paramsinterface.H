/*---------------------------------------------------------------------*/
/*!
\file contact_paramsinterface.H

\brief Contact parameter interface. Necessary for the communication
       between the structural time integration framework and the
       contact strategies.

\level 3

\maintainer Matthias Mayr

*/
/*---------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_CONTACT_PARAMSINTERFACE_H_
#define SRC_DRT_CONTACT_CONTACT_PARAMSINTERFACE_H_

#include "../drt_mortar/mortar_paramsinterface.H"

namespace NOX
{
  namespace NLN
  {
    enum class CorrectionType : int;
  }  // namespace NLN
}  // namespace NOX
namespace INPAR
{
  namespace STR
  {
    enum PredEnum : int;
  }  // namespace STR
  namespace CONTACT
  {
    enum VariationalApproach : int;
  }  // namespace CONTACT
}  // namespace INPAR

namespace STR
{
  namespace MODELEVALUATOR
  {
    class Generic;
  }  // namespace MODELEVALUATOR
}  // namespace STR

namespace CONTACT
{
  class ParamsInterface : public MORTAR::ParamsInterface
  {
   public:
    //! destructor
    virtual ~ParamsInterface(){};

    //! get the predictor status
    virtual bool IsPredictor() const = 0;

    //! Is the current state coming from a predictor step?
    virtual bool IsPredictorState() const = 0;

    //! \brief get the currently active predictor type
    /** \note If the execution of the predictor is finished, this
     *  function will return INPAR::STR::pred_vague. \author hiermeier */
    virtual enum INPAR::STR::PredEnum GetPredictorType() const = 0;

    //! get the current step length
    virtual double GetStepLength() const = 0;

    //! get number of linear system corrections (modified Newton approach)
    virtual int GetNumberOfModifiedNewtonCorrections() const = 0;

    //! get the is_default_step indicator
    virtual bool IsDefaultStep() const = 0;

    //! get correction type
    virtual NOX::NLN::CorrectionType GetCorrectionType() const = 0;

    //! get the current time step [derived]
    virtual double GetDeltaTime() const = 0;

    //! get a pointer to the contact model evaluator
    virtual const STR::MODELEVALUATOR::Generic& GetModelEvaluator() const = 0;

    //! get the output file path
    virtual std::string GetOutputFilePath() const = 0;

    //! get the variational approach type
    virtual enum INPAR::CONTACT::VariationalApproach GetVariationalApproachType() const = 0;

    //! get the variational approach type
    virtual void SetVariationalApproachType(
        const enum INPAR::CONTACT::VariationalApproach var_type) = 0;
  };
}  // namespace CONTACT


#endif /* SRC_DRT_CONTACT_CONTACT_PARAMSINTERFACE_H_ */
