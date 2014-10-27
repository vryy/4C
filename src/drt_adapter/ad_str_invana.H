/*----------------------------------------------------------------------*/
/*!
\file ad_str_invana.H

\brief Wrapper for the structural time integration to be used
in the inverse analysis

<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef ADAPTER_STRUCTURE_INVANA_H
#define ADAPTER_STRUCTURE_INVANA_H

#include "ad_str_timeloop.H"

namespace ADAPTER
{
  //! documentation
  class StructureInvana : public StructureTimeLoop
  {
  public:

    /// constructor
    explicit StructureInvana(Teuchos::RCP<Structure> structure);

    virtual void PostPredict();

    virtual void PostUpdate();

    Teuchos::RCP<DRT::UTILS::TimIntMStep<Epetra_Vector> > DispSteps(){return sdis_;}

    Teuchos::RCP<DRT::UTILS::TimIntMStep<double> > TimeSteps(){return stime_;}


  private:

    /// resize storage of quantities to be passed out of the time integrator
    void ResizeStorage();

    //* Storage of quantities:
    //! time_
    Teuchos::RCP<DRT::UTILS::TimIntMStep<double> > stime_;

    //! displacements
    Teuchos::RCP<DRT::UTILS::TimIntMStep<Epetra_Vector> > sdis_;

    //! veclocities
    Teuchos::RCP<DRT::UTILS::TimIntMStep<Epetra_Vector> > svel_;

    //! number of steps stored
    int sizestorage_;

  };

}

#endif
