/*----------------------------------------------------------------------*/
/*!
\file scatra_utils.H

\brief utility functions for scalar transport problems

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef SCATRA_UTILS_H
#define SCATRA_UTILS_H

namespace SCATRA
{
/// automatically create a matching scalar transport discretization 
/// from an already existing fluid discretization
void CreateScaTraDiscretization(
    Teuchos::RCP<DRT::Discretization>& fluiddis,  ///< RCP to fluid discretization
    Teuchos::RCP<DRT::Discretization>& scatradis, ///< RCP to empty scalar transport discretization
    std::map<string,string>& conditions_to_copy,  ///< condition names to copy (and new names)
    bool makequadratic                            ///< increase polynomial order
    );

}

#endif
#endif
