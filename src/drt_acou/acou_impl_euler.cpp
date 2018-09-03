/*!----------------------------------------------------------------------
\file acou_impl_euler.cpp
\brief specific acoustic implicit euler implementation

<pre>
\level 2

\maintainer Luca Berardocco
            berardoccoo@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
</pre>

*----------------------------------------------------------------------*/

#include "acou_impl_euler.H"

/*----------------------------------------------------------------------*
 |  Constructor (public)                                 schoeder 01/14 |
 *----------------------------------------------------------------------*/
ACOU::TimIntImplEuler::TimIntImplEuler(const Teuchos::RCP<DRT::DiscretizationHDG>& actdis,
    const Teuchos::RCP<LINALG::Solver>& solver, const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output)
    : AcouImplicitTimeInt(actdis, solver, params, output)
{
}
