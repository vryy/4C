/*----------------------------------------------------------------------*/
/*!
\brief Solve FSI problems using a Robin partitioning approach

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef FSI_ROBIN_H
#define FSI_ROBIN_H

#include "fsi_partitioned.H"

namespace FSI {

/// Robin coupled interface system
class Robin : public Partitioned
{
public:

  explicit Robin(Epetra_Comm& comm);

protected:

  /// composed FSI operator
  virtual void FSIOp(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag);

  /// interface fluid operator
  virtual Teuchos::RCP<Epetra_Vector> FluidOp(Teuchos::RCP<Epetra_Vector> idisp,
                                              const FillType fillFlag);

  /// interface structural operator
  virtual Teuchos::RCP<Epetra_Vector> StructOp(Teuchos::RCP<Epetra_Vector> iforce,
                                               const FillType fillFlag);

private:

  /// do robin coupling on fluid field
  bool fluidrobin_;

  /// do robin coupling on structural field
  bool structrobin_;
};

}

#endif
#endif
