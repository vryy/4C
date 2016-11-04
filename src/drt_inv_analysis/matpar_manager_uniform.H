/*----------------------------------------------------------------------*/
/*!
\file matpar_manager_uniform.H

\brief Patch-wise uniform distribution of material parameters

<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/
/* definitions */
#ifndef MATPAR_MANAGER_UNIFORM_H
#define MATPAR_MANAGER_UNIFORM_H

/*----------------------------------------------------------------------*/
/* headers */
#include "matpar_manager.H"


// forward declarations
namespace DRT
{
  class Discretization;
}

namespace INVANA

{
/*! \class MatParManagerUniform
 *  \brief Class to have blockwise constant material parameters
 *
 *  This class handles optimization parameters which represent blockwise
 *  patches of material. The patches must be given by different input file
 *  materials. The final unique vector of parameters is kept on proc 0.
 *
 *  For the parallel computation of derivatives, a completely redundant version
 *  of the parameters is used which is then later 'contracted' by INVANA::Finalize.
 *
 *  \author kehl \date somewhen2013
 */
class MatParManagerUniform  : public MatParManager
{
public:

  //! Constructor
  MatParManagerUniform(Teuchos::RCP<DRT::Discretization> discret);

  //! Destructor
  virtual ~MatParManagerUniform(){}

  //! Setup specific parametrization layout
  virtual void Setup();

  //! Account for distributed chain-rule application
  virtual void Finalize(Teuchos::RCP<Epetra_MultiVector> source,
      Teuchos::RCP<Epetra_MultiVector> target
      );

protected:

  //! Get current material parameters
  virtual void FillParameters(Teuchos::RCP<Epetra_MultiVector> params);

  //! Initialize optimization parameters
  virtual void InitParameters(int parapos, double val);

  //! Chain rule application
  virtual void ContractGradient(Teuchos::RCP<Epetra_MultiVector> dfint,
      double val, int elepos, int paraposglobal, int paraposlocal);

private:

  //! build a graph of elements in the optimization parameter vector
  virtual void FillAdjacencyMatrix(const Epetra_Map& paramrowmap, Teuchos::RCP<Epetra_CrsMatrix> graph);

};
}  // namespace INVANA

/*----------------------------------------------------------------------*/
#endif
