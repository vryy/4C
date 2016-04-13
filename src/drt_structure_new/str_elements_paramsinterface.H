/*-----------------------------------------------------------*/
/*!
\file str_elements_paramsinterface.H

\maintainer Michael Hiermeier

\date Mar 24, 2016

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_DRT_STRUCTURE_NEW_STR_ELEMENTS_PARAMSINTERFACE_H_
#define SRC_DRT_STRUCTURE_NEW_STR_ELEMENTS_PARAMSINTERFACE_H_

#include "../drt_inpar/inpar_structure.H"         // enums
#include "../drt_lib/drt_elements_paramsinterface.H"    // base class
#include "../solver_nonlin_nox/nox_nln_statustest_factory.H" // enums

namespace STR {
namespace ELEMENTS {
  //! evaluation error flags
  enum EvalErrorFlag
  {
    ele_error_none = 0,                     //!< no error occurred (default)
    ele_error_negative_def_gradient = 1     //!< a negative deformation gradient was detected
  }; // EvalErrorFlag

  //! Map evaluation error flag to a std::string
  static inline const std::string EvalErrorFlag2String(const enum EvalErrorFlag& errorflag)
  {
    switch (errorflag)
    {
      case ele_error_none:
        return "ele_error_none";
        break;
      case ele_error_negative_def_gradient:
        return "ele_error_negative_def_gradient";
        break;
      default:
        return "unknown";
        break;
    }
    return "";
  }; // EvalErrorFlag2String

  /*! \brief Parameter interface for the structural elements and the STR::Integrator data exchange
   *
   *  This class is a special case of the DRT::ELEMENTS::Interface class and gives you all the
   *  basic function definitions which you can use to get access to the STR::Integrator and many
   *  more objects.
   *  Please consider to derive a special interface class, if you need special parameters inside
   *  of your element. Keep the Evaluate call untouched and cast the interface object to the
   *  desired specification.
   *
   *  ToDo Currently we set the interface in the elements via the Teuchos::ParameterList.
   *  Theoretically, the Teuchos::ParameterList can be replaced by the interface itself!
   *
   *  \date 03/2016
   *  \author hiermeier */
  class ParamsInterface : public DRT::ELEMENTS::ParamsInterface
  {
  public:
    //! return the damping type
    virtual enum INPAR::STR::DampKind GetDampingType() const = 0;

    virtual bool IsTolerateErrors() const = 0;

    //! @name General time integration parameters
    //! @{
    virtual double GetTimIntFactorDisp() const = 0;

    virtual double GetTimIntFactorVel() const = 0;
    //! @}

    //! @name Access control parameters for the handling of element internal variables (e.g. EAS)
    //! @{
    //! get the current step length
    virtual double GetStepLength() const = 0;

    //! Is the current step a default step, or e.g. a line search step?
    virtual bool IsDefaultStep() const = 0;
    //! @}

    //! @name Set functions
    //! @{
    /*! \brief set evaluation error flag
     *
     *  See the EvalErrorFlag enumerators for more information. */
    virtual void SetEleEvalErrorFlag(const enum EvalErrorFlag& error_flag) = 0;
    //! @}

    //! @name output related functions
    //! @{
    virtual Teuchos::RCP<std::vector<char> >& MutableStressDataPtr() = 0;

    virtual Teuchos::RCP<std::vector<char> >& MutableStrainDataPtr() = 0;

    virtual Teuchos::RCP<std::vector<char> >& MutablePlasticStrainDataPtr() = 0;

    //! get the current stress type
    virtual enum INPAR::STR::StressType GetStressOutputType() const = 0;

    //! get the current strain type
    virtual enum INPAR::STR::StrainType GetStrainOutputType() const = 0;

    //! get the current plastic strain type
    virtual enum INPAR::STR::StrainType GetPlasticStrainOutputType() const = 0;

    //! add the current partial update norm of the given quantity
    virtual void SumIntoMyUpdateNorm(
        const enum    NOX::NLN::StatusTest::QuantityType& qtype,
        const int&    numentries,
        const double* my_update_values,
        const double* my_new_sol_values,
        const double& step_length,
        const int&    owner) = 0;

    //! collects and calculates the solution norm of the previous accepted Newton step on the current proc
    virtual void SumIntoMyPreviousSolNorm(
        const enum    NOX::NLN::StatusTest::QuantityType& qtype,
        const int&    numentries,
        const double* my_old_values,
        const int&    owner) = 0;
    //! @}
  }; // class ParamsInterface
} // namespace ELEMENTS
} // namespace STR



#endif /* SRC_DRT_STRUCTURE_NEW_STR_ELEMENTS_PARAMSINTERFACE_H_ */
