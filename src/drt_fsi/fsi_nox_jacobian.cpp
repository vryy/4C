/*----------------------------------------------------------------------*/
/*!

\brief Matrix-free Newton-Krylov for FSI

\level 1

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------*/

#include "fsi_nox_jacobian.H"

#include <iostream>

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_RowMatrix.h>
#include <NOX_Abstract_Group.H>
#include <NOX_Epetra_Vector.H>
#include <NOX_Epetra_VectorSpace.H>
#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Utils.H>


NOX::FSI::FSIMatrixFree::FSIMatrixFree(Teuchos::ParameterList& printParams,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& i, const NOX::Epetra::Vector& x)
    : label("FSI-Matrix-Free"),
      interface(i),
      currentX(x),
      perturbX(x),
      perturbY(x),
      useGroupForComputeF(false),
      utils(printParams)
{
  perturbX.init(0.0);
  perturbY.init(0.0);

  // Epetra_Operators require Epetra_Maps, so anyone using block maps
  // (Epetra_BlockMap) won't be able to directly use the AztecOO solver.
  // We get around this by creating an Epetra_Map from the Epetra_BlockMap.
  const Epetra_Map* testMap = 0;
  testMap = dynamic_cast<const Epetra_Map*>(&currentX.getEpetraVector().Map());
  if (testMap != 0)
  {
    epetraMap = Teuchos::rcp(new Epetra_Map(*testMap));
  }
  else
  {
    int size = currentX.getEpetraVector().Map().NumGlobalPoints();
    int mySize = currentX.getEpetraVector().Map().NumMyPoints();
    int indexBase = currentX.getEpetraVector().Map().IndexBase();
    const Epetra_Comm& comm = currentX.getEpetraVector().Map().Comm();
    epetraMap = Teuchos::rcp(new Epetra_Map(size, mySize, indexBase, comm));
  }
}


NOX::FSI::FSIMatrixFree::~FSIMatrixFree() {}


int NOX::FSI::FSIMatrixFree::SetUseTranspose(bool UseTranspose)
{
  if (UseTranspose == true)
  {
    utils.out()
        << "ERROR: FSIMatrixFree::SetUseTranspose() - Transpose is unavailable in Matrix-Free mode!"
        << std::endl;
    throw "NOX Error";
  }
  return (-1);
}


int NOX::FSI::FSIMatrixFree::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  // Calculate the matrix-vector product:
  //
  // y = R' x = S'(F(d)) F'(d) x - x
  //
  // that comes down to a FSI residuum call with linear field solvers.
  //
  // We make use of the special structure of the FSI Residuum (this
  // approach is not general purpose) and neglect the dependence of
  // the fluid field on the interface displacements.

  // Convert X and Y from an Epetra_MultiVector to a Epetra_Vectors
  // and NOX::epetra::Vectors.  This is done so we use a consistent
  // vector space for norms and inner products.
  Teuchos::RCP<Epetra_Vector> wrappedX = Teuchos::rcp(new Epetra_Vector(View, X, 0));
  Teuchos::RCP<Epetra_Vector> wrappedY = Teuchos::rcp(new Epetra_Vector(View, Y, 0));
  NOX::Epetra::Vector nevX(wrappedX, NOX::Epetra::Vector::CreateView);
  NOX::Epetra::Vector nevY(wrappedY, NOX::Epetra::Vector::CreateView);

  // The trial vector x is not guaranteed to be a suitable interface
  // displacement. It might be much too large to fit the ALE
  // algorithm. But we know our residual to be linear, so we can
  // easily scale x.

  double xscale = 1e4 * nevX.norm();
  // double xscale = nevX.norm();
  if (xscale == 0)
  {
    // In the first call is x=0. No need to calculate the
    // residuum. y=0 in that case.
    nevY.init(0.);
    return 0;
  }

  // For some strange reason currentX.Map()!=X.Map() and we are bound
  // to call computeF with the right map.
  perturbX = currentX;
  // perturbX.update(1./xscale,nevX,0.0);
  perturbX.update(1., nevX, 0.0);

  if (!useGroupForComputeF)
  {
    interface->computeF(perturbX.getEpetraVector(), perturbY.getEpetraVector(),
        NOX::Epetra::Interface::Required::User);
  }
  else
  {
    groupPtr->setX(perturbX);
    groupPtr->computeF();
    perturbY = groupPtr->getF();
  }

  // scale back
  // nevY.update(xscale, perturbY, 0.0);
  nevY.update(1., perturbY, 0.0);

  return 0;
}


int NOX::FSI::FSIMatrixFree::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  utils.out() << "ERROR: FSIMatrixFree::ApplyInverse - Not available for Matrix Free!" << std::endl;
  throw "NOX Error";
  return (-1);
}


double NOX::FSI::FSIMatrixFree::NormInf() const
{
  utils.out() << "ERROR: FSIMatrixFree::NormInf() - Not Available for Matrix-Free mode!"
              << std::endl;
  throw "NOX Error";
  return 1.0;
}


const char* NOX::FSI::FSIMatrixFree::Label() const { return label.c_str(); }


bool NOX::FSI::FSIMatrixFree::UseTranspose() const { return false; }


bool NOX::FSI::FSIMatrixFree::HasNormInf() const { return false; }


const Epetra_Comm& NOX::FSI::FSIMatrixFree::Comm() const
{
  return currentX.getEpetraVector().Map().Comm();
}


const Epetra_Map& NOX::FSI::FSIMatrixFree::OperatorDomainMap() const { return *epetraMap; }


const Epetra_Map& NOX::FSI::FSIMatrixFree::OperatorRangeMap() const { return *epetraMap; }


bool NOX::FSI::FSIMatrixFree::computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac)
{
  // Remember the current interface displacements.
  currentX = x;

  // Nothing to do here. The work is done when we apply a vector to
  // the Jacobian.
  bool ok = true;
  return ok;
}


void NOX::FSI::FSIMatrixFree::setGroupForComputeF(const NOX::Abstract::Group& group)
{
  useGroupForComputeF = true;
  groupPtr = group.clone();
  return;
}
