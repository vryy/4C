/*-----------------------------------------------------------*/
/*!
\file str_predict_generic.H

\brief Generic class for all predictors.

\maintainer Michael Hiermeier

\date Aug 13, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STR_PREDICT_GENERIC_H_
#define STR_PREDICT_GENERIC_H_

#include <Teuchos_RCP.hpp>
// necessary due to enums
#include "../drt_inpar/inpar_structure.H"

// forward declaration ...
class Epetra_Vector;
namespace NOX {
namespace Abstract {
  class Group;
} // namespace Abstract
} // namespace NOX
namespace STR {
  class Dbc;
namespace IMPLICIT {
  class Generic;
} // namespace IMPLICIT
namespace TIMINT {
  class BaseDataGlobalState;
  class BaseDataIO;
} // namespace TIMINT
namespace PREDICT {
  class Generic
  {
  public:
    //! constructor
    Generic();

    //! destructor
    virtual ~Generic() {};

    //! initialize the base class variables
    virtual void Init(const enum INPAR::STR::PredEnum& type,
        const Teuchos::RCP<STR::IMPLICIT::Generic>& implint_ptr,
        const Teuchos::RCP<STR::Dbc>& dbc_ptr,
        const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
        const Teuchos::RCP<STR::TIMINT::BaseDataIO>& iodata_ptr,
        const Teuchos::RCP<Teuchos::ParameterList>& noxparams_ptr);

    //! setup of the specific predictor
    virtual void Setup() = 0;

    //! Get the predictor type enum
    const INPAR::STR::PredEnum& GetType() const { return type_; };

    //! returns the name of the used predictor
    virtual const std::string Name() const;

    //! Preprocess the predictor step
    virtual void PrePredict(NOX::Abstract::Group& grp);

    //! Pre-/Postprocess the specific predictor step
    void Predict(NOX::Abstract::Group& grp);

    //! Calculate the specific predictor step
    virtual void Compute(NOX::Abstract::Group& grp) = 0;

    //! Postprocess the predictor step
    virtual void PostPredict(NOX::Abstract::Group& grp);

    //! return a constant reference to the global state object (read only)
    const STR::TIMINT::BaseDataGlobalState& GlobalState() const;

    //! print the result of the predictor step
    void Print() const;

    //! Run before the external force are computed and assembled
    virtual bool PreApplyForceExternal( Epetra_Vector& fextnp ) const;

  protected:
    //! returns init state
    const bool& IsInit() const { return isinit_; };

    //! returns setup state
    const bool& IsSetup() const { return issetup_; };

    void CheckInit() const;

    void CheckInitSetup() const;

    Teuchos::RCP<STR::IMPLICIT::Generic>& ImplIntPtr();
    STR::IMPLICIT::Generic& ImplInt();

    Teuchos::RCP<STR::Dbc>& DbcPtr();
    STR::Dbc& Dbc();

    Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& GlobalStatePtr();
    STR::TIMINT::BaseDataGlobalState& GlobalState();

    Teuchos::RCP<STR::TIMINT::BaseDataIO>& IODataPtr();
    STR::TIMINT::BaseDataIO& IOData();

    Teuchos::RCP<Teuchos::ParameterList>& NoxParamsPtr();
    Teuchos::ParameterList& NoxParams();

  protected:
    //! indicates if the Init() function has been called
    bool isinit_;

    //! indicates if the Setup() function has been called
    bool issetup_;

  private:
    //! predictor type
    enum INPAR::STR::PredEnum type_;

    //! pointer to the implicit integrator
    Teuchos::RCP<STR::IMPLICIT::Generic> implint_ptr_;

    //! pointer to the dirichlet boundary condition object
    Teuchos::RCP<STR::Dbc> dbc_ptr_;

    //! global state pointer
    Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> gstate_ptr_;

    //! input/output data pointer
    Teuchos::RCP<STR::TIMINT::BaseDataIO> iodata_ptr_;

    Teuchos::RCP<Teuchos::ParameterList> noxparams_ptr_;
  };  // class  Generic
} // namespace PREDICT
} // namespace STR


#endif /* STR_PREDICT_GENERIC_H_ */
