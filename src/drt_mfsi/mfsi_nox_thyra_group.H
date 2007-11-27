/*----------------------------------------------------------------------*/
/*!
\file mfsi_nox_thyra_group.H

\brief NOX Thyra group enhancement

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef MFSI_NOX_THYRA_GROUP_H
#define MFSI_NOX_THYRA_GROUP_H

#include <NOX_Thyra_Group.H>

namespace MFSI {

  // forward declaration
  class Algorithm;

  /// slightly enhanced NOX group
  /*!

    This class fakes the behaviour of the original NOX group to
    account for our finite element implementation that always
    calculates the field matrix and the right hand side at the same
    time. If NOX asks for one of them, we set both and flag them as
    current.

    \author u.kue
    \date 11/07
   */
  class NOX_Thyra_Group : public NOX::Thyra::Group
  {
  public:
    NOX_Thyra_Group(const NOX::Thyra::Vector& initial_guess,
                    const Algorithm* mfsi);

    //! @name compute functions
    /// These functions always calculate both, RHS and Jacobian.

    NOX::Abstract::Group::ReturnType computeF();
    NOX::Abstract::Group::ReturnType computeJacobian();

    //@}

  private:
    const Algorithm* mfsi_;
  };
}

#endif
#endif
