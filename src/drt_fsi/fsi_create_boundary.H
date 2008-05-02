/*----------------------------------------------------------------------*/
/*!
\file fsi_create_boundary.H

\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef FSI_CREATE_BOUNDARY_H
#define FSI_CREATE_BOUNDARY_H

#include <string>
#include <Teuchos_RCP.hpp>

namespace DRT
{
  class Discretization;
}

//! create discretization from a given condition
Teuchos::RCP<DRT::Discretization> CreateDiscretizationFromCondition(
  Teuchos::RCP<DRT::Discretization> cutterdis,    ///< condition with discretization
  const std::string&                condname,     ///< name of the condition, by which the derived discretization is identified
  const std::string&                discret_name, ///< name of the new condition
  const std::string&                element_name  ///< name/type of the elements to be created
  );

#endif
#endif
