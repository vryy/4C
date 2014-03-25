  /*!----------------------------------------------------------------------
\file immersed_partitioned_fsi_neumannneumann.H

\brief partitioned immersed fsi algorithm for neumann-neumann coupling (volume force coupling)

<pre>
Maintainers: Andreas Rauch & Anh-Tu Vuong
             {rauch,vuong}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289--15240 / 15264
</pre>
*----------------------------------------------------------------------*/
#ifndef IMMERSED_PARTITIONED_FSI_NEUMANNNEUMANN_H
#define IMMERSED_PARTITIONED_FSI_NEUMANNNEUMANN_H

#include "immersed_partitioned_fsi.H"
#include "../drt_fsi/fsi_partitioned.H"

namespace IMMERSED {

class ImmersedPartitionedFSINeumannNeumann : public ImmersedPartitionedFSI, public FSI::Partitioned
{

public:

  /// setup partitioned ImmersedFSI algorithm
  explicit ImmersedPartitionedFSINeumannNeumann(const Epetra_Comm& comm);

protected:

  /// composed FSI operator
  virtual void FSIOp(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag);

  /// interface fluid operator
  virtual Teuchos::RCP<Epetra_Vector> FluidOp(Teuchos::RCP<Epetra_Vector> bforce,
                                              const FillType fillFlag);

  /// interface structural operator
  virtual Teuchos::RCP<Epetra_Vector> StructOp(Teuchos::RCP<Epetra_Vector> iforce,
                                               const FillType fillFlag);
  /// initial guess
  virtual Teuchos::RCP<Epetra_Vector> InitialGuess();

  //! @name Various global forces
  //@{
  Teuchos::RCP<Epetra_Vector> fintn_;   //!< internal structural force at time n+1
  Teuchos::RCP<Epetra_Vector> ffsi_;    //!< body force vector being applied to fluid field
  //@}


};// class ImmersedFSI
} // namspace IMMERSED

#endif
