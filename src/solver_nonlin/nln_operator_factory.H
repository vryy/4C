/*----------------------------------------------------------------------------*/
/*!
\file nln_operator_factory.H

\brief Factory for nonlinear operators

\level 3

\maintainer Matthias Mayr
*/

/*----------------------------------------------------------------------------*/

#ifndef NLN_OPERATOR_FACTORY_H
#define NLN_OPERATOR_FACTORY_H

/*----------------------------------------------------------------------------*/
/* headers */

// Teuchos
#include <Teuchos_RCP.hpp>

/*----------------------------------------------------------------------------*/
/* forward declarations */

namespace NLNSOL
{
  class NlnOperatorBase;

  namespace UTILS
  {
    class NlnConfig;
  }  // namespace UTILS
}  // namespace NLNSOL

namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
  /*! \class NlnOperatorFactory
   *  \brief Factory for nonlinear operators
   *
   *  Factory class to create any kind of nonlinear operator. This class just
   *  implements a Create()-routine that creates the wanted object and returns
   *  it as a NLNSOL::NlnOperatorBase. The choice which object has to be created
   *  is based on an entry in an Teuchos::ParameterList.
   *
   *  Add any future nonlinear operator to this factory class. Use only the
   *  NLNSOL::NlnOperatorFactory to create nonlinear operators, never call there
   *  constructors explicitly.
   *
   *  \author mayr.mt \date 06/2014
   */
  class NlnOperatorFactory
  {
   public:
    //! @name Construction / Destruction
    //@{

    //! Constructor (empty)
    NlnOperatorFactory();

    //! Destructor
    virtual ~NlnOperatorFactory(){};

    //@}

    //! @name Creation on nonlinear operators
    //@{

    /*! \brief Create a user-given nonlinear operator
     *
     *  The nonlinear operator to be created is configured via the parameter
     *  "nonlinear operator: type" in the parameter list \c list in the global configuration
     *  list \c config.
     *
     *  \author mayr.mt \date 10/2015
     */
    virtual Teuchos::RCP<NLNSOL::NlnOperatorBase> Create(
        Teuchos::RCP<const NLNSOL::UTILS::NlnConfig> config,  ///< configuration of nonlinear solver
        const std::string listname  ///< name of sublist with configuration of nonlinear operator
    );

    //@}

   protected:
   private:
  };  // class NlnOperatorFactory
}  // namespace NLNSOL

#endif /* NLN_OPERATOR_FACTORY_H */
