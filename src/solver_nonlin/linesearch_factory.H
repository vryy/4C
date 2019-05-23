/*----------------------------------------------------------------------------*/
/*!

\brief Line search factory class

\level 3

\maintainer Matthias Mayr
*/

/*----------------------------------------------------------------------------*/

#ifndef LINESEARCH_FACTORY_H
#define LINESEARCH_FACTORY_H

/*----------------------------------------------------------------------------*/
/* headers */

// Teuchos
#include <Teuchos_RCP.hpp>

/*----------------------------------------------------------------------------*/
/* forward declarations */

namespace NLNSOL
{
  class LineSearchBase;
}  // namespace NLNSOL

namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
  /*! \class LineSearchFactory
   *  \brief Line search factory class
   *
   *  Factory class to create any kind of line search algorithm. This class just
   *  implements a Create()-routine that creates the wanted object and returns
   *  it as a NLNSOL::LineSearchBase. The choice which object has to be created
   *  is based on an entry in an Teuchos::ParameterList.
   *
   *  Add any future line search algorithms to this factory class. Use only the
   *  NLNSOL::LineSearchFactory to create line search objects, never call there
   *  constructors explicitly.
   *
   * \author mayr.mt \date 05/2014
   */
  class LineSearchFactory
  {
   public:
    //! @name Construction / Destruction
    //@{

    //! constructor (empty)
    LineSearchFactory();

    //! Destructor
    virtual ~LineSearchFactory(){};

    //@}

    //! @name Creation of line search algorithms
    //@{

    /*! \brief Create a user-given line search algorithm
     *
     *  The line search to be created is configured via the parameter "line
     *  search" in the parameter list \c listname in the global configuration
     *  list.
     *
     *  \author mayr.mt \date 10/2015
     */
    virtual Teuchos::RCP<NLNSOL::LineSearchBase> Create(
        Teuchos::RCP<const NLNSOL::UTILS::NlnConfig> config,  ///< configuration of nonlinear solver
        const std::string listname  ///< name of sublist with configuration of line search algorithm
    );

    //@}

   protected:
   private:
  };  // class LineSearchFactory
}  // namespace NLNSOL

#endif /* LINESEARCH_FACTORY_H */
