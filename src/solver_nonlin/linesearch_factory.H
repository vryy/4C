/*----------------------------------------------------------------------*/
/*!
\file linesearch_factory.H

\brief Line search factory class

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------*/

#ifndef LINESEARCH_FACTORY_H_
#define LINESEARCH_FACTORY_H_

/*----------------------------------------------------------------------*/
/* headers */

// Teuchos
#include <Teuchos_RCP.hpp>

// baci
#include "linesearch_base.H"

/*----------------------------------------------------------------------*/
/* forward declarations */

namespace Teuchos
{
  class ParameterList;
}

/*----------------------------------------------------------------------*/
/* definition of classes */
namespace NLNSOL
{
  /*!
   * \brief Line search factory class
   *
   * Choose from parameter list, which line search object is needed. Create the
   * line search algorithmic object and deliver it.
   *
   * \author mayr.mt
   * \date 05/2014
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

    virtual Teuchos::RCP<NLNSOL::LineSearchBase> Create(
        const Teuchos::ParameterList& params ///< line search parameter list
        );

    //@}

  protected:

  private:

  };
}

#endif /* LINESEARCH_FACTORY_H_ */
