
#include "fsi_nox_jacobian.H"

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_RowMatrix.h>
#include <NOX_Abstract_Group.H>
#include <NOX_Epetra_Vector.H>
#include <NOX_Epetra_VectorSpace.H>
#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Utils.H>


FSIMatrixFree::FSIMatrixFree(Teuchos::ParameterList& printParams,
                             const Teuchos::RefCountPtr<NOX::Epetra::Interface::Required>& i,
                             const NOX::Epetra::Vector& x) :
  label("FSI-Matrix-Free"),
  interface(i),
  currentX(x),
  useGroupForComputeF(false),
  utils(printParams)
{
  // Epetra_Operators require Epetra_Maps, so anyone using block maps
  // (Epetra_BlockMap) won't be able to directly use the AztecOO solver.
  // We get around this by creating an Epetra_Map from the Epetra_BlockMap.
  const Epetra_Map* testMap = 0;
  testMap = dynamic_cast<const Epetra_Map*>(&currentX.getEpetraVector().Map());
  if (testMap != 0) {
    epetraMap = Teuchos::rcp(new Epetra_Map(*testMap));
  }
  else {
    int size = currentX.getEpetraVector().Map().NumGlobalPoints();
    int mySize = currentX.getEpetraVector().Map().NumMyPoints();
    int indexBase = currentX.getEpetraVector().Map().IndexBase();
    const Epetra_Comm& comm = currentX.getEpetraVector().Map().Comm();
    epetraMap = Teuchos::rcp(new Epetra_Map(size, mySize, indexBase, comm));
  }

}

FSIMatrixFree::~FSIMatrixFree()
{
}


int FSIMatrixFree::SetUseTranspose(bool UseTranspose)
{
  if (UseTranspose == true) {
    utils.out() << "ERROR: FSIMatrixFree::SetUseTranspose() - Transpose is unavailable in Matrix-Free mode!"
                << endl;
    throw "NOX Error";
  }
  return (-1);
}

int FSIMatrixFree::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  // Calculate the matrix-vector product:
  //
  // y = R' x = x - S'(F(d)) F'(d) x
  //
  // that comes down to a FSI residuum call with linear field solvers.
  //
  // We make use of the special structure of the FSI Residuum (this
  // approach is not general purpose) and neglect the dependence of
  // the fluid field on the interface displacements.

  // Convert X and Y from an Epetra_MultiVector to a Epetra_Vectors
  // and NOX::epetra::Vectors.  This is done so we use a consistent
  // vector space for norms and inner products.
  Teuchos::RefCountPtr<Epetra_Vector> wrappedX = Teuchos::rcp(new Epetra_Vector(View, X, 0));
  Teuchos::RefCountPtr<Epetra_Vector> wrappedY = Teuchos::rcp(new Epetra_Vector(View, Y, 0));
  NOX::Epetra::Vector nevX(wrappedX, NOX::Epetra::Vector::CreateView);
  NOX::Epetra::Vector nevY(wrappedY, NOX::Epetra::Vector::CreateView);

  if (!useGroupForComputeF)
    interface->computeF(nevX.getEpetraVector(), nevY.getEpetraVector(),
			NOX::Epetra::Interface::Required::User);
  else
  {
    groupPtr->setX(nevX);
    groupPtr->computeF();
    nevY = groupPtr->getF();
  }

  //nevY.update(-1.0, fp, 1.0);

#if 0
  // Use a directional derivative to compute y = Jx
  /*
   * eta = scalar perturbation
   * u = solution vector used to evaluate f
   * f = function evaluation (RHS)
   * x = vector that J is applied to
   *
   *        f(u+eta*x) - f(u)
   * Jx =   -----------------
   *               eta
   */


  // Compute perturbation constant, eta
  // Taken from LOCA v1.0 manual SAND2002-0396 p. 28 eqn. 2.43
  // eta = lambda*(lambda + 2norm(u)/2norm(x))
  double solutionNorm = 1.0;
  double vectorNorm = 1.0;

  solutionNorm = currentX.norm();
  vectorNorm = currentX.getVectorSpace()->norm(*wrappedX);

  // Make sure the norm is not zero, otherwise we can get an inf perturbation
  if (vectorNorm == 0.0) {
    //utils.out(Utils::Warning) << "Warning: NOX::Epetra::FSIMatrixFree::Apply() - vectorNorm is zero" << endl;
    vectorNorm = 1.0;
  }

#if 0
  // Create an extra perturbed residual vector pointer if needed
  if ( diffType == Centered )
    if ( Teuchos::is_null(fmPtr) )
      fmPtr = Teuchos::rcp(new NOX::Epetra::Vector(fo));
#endif

  double scaleFactor = 1.0;
#if 0
  if ( diffType == Backward )
  scaleFactor = -1.0;
#endif

  if (computeEta) {
    if (useNewPerturbation) {
      double dotprod = currentX.getVectorSpace()->
	innerProduct(currentX.getEpetraVector(), *wrappedX);
      if (dotprod==0.0)
	dotprod = 1.0e-12;
      eta = lambda*(1.0e-12/lambda + fabs(dotprod)/(vectorNorm * vectorNorm))
	* dotprod/fabs(dotprod);
    }
    else
      eta = lambda*(lambda + solutionNorm/vectorNorm);
  }
  else
    eta = userEta;

  // Compute the perturbed RHS
  perturbX = currentX;
  Y = X;
  Y.Scale(eta);
  perturbX.update(1.0,nevY,1.0);

  if (!useGroupForComputeF)
    interface->computeF(perturbX.getEpetraVector(), fp.getEpetraVector(),
			NOX::Epetra::Interface::Required::MF_Res);
  else{
    groupPtr->setX(perturbX);
    groupPtr->computeF();
    fp = dynamic_cast<const NOX::Epetra::Vector&>
      (groupPtr->getF());
  }

#if 0
  if ( diffType == Centered ) {
    Y.Scale(-2.0);
    perturbX.update(scaleFactor,nevY,1.0);
    if (!useGroupForComputeF)
      interface->computeF(perturbX.getEpetraVector(), fmPtr->getEpetraVector(),
			  NOX::Epetra::Interface::Required::MF_Res);
    else{
      groupPtr->setX(perturbX);
      groupPtr->computeF();
      *fmPtr = dynamic_cast<const NOX::Epetra::Vector&>
        (groupPtr->getF());
    }
  }
#endif

  // Compute the directional derivative
  if ( diffType != Centered ) {
    nevY.update(1.0, fp, -1.0, fo, 0.0);
    nevY.scale( 1.0/(scaleFactor * eta) );
  }
  else {
    nevY.update(1.0, fp, -1.0, *fmPtr, 0.0);
    nevY.scale( 1.0/(2.0 * eta) );
  }

#endif
  return 0;
}

int FSIMatrixFree::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  utils.out() << "ERROR: FSIMatrixFree::ApplyInverse - Not available for Matrix Free!"
              << endl;
  throw "NOX Error";
  return (-1);
}

double FSIMatrixFree::NormInf() const
{
  utils.out() << "ERROR: FSIMatrixFree::NormInf() - Not Available for Matrix-Free mode!" << endl;
  throw "NOX Error";
  return 1.0;
}


const char* FSIMatrixFree::Label () const
{
  return label.c_str();
}

bool FSIMatrixFree::UseTranspose() const
{
  return false;
}

bool FSIMatrixFree::HasNormInf() const
{
  return false;
}

const Epetra_Comm & FSIMatrixFree::Comm() const
{
  return currentX.getEpetraVector().Map().Comm();
}

const Epetra_Map& FSIMatrixFree::OperatorDomainMap() const
{
  return *epetraMap;
}

const Epetra_Map& FSIMatrixFree::OperatorRangeMap() const
{
  return *epetraMap;
}

bool FSIMatrixFree::computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac)
{
#if 0
  // Since we have no explicit Jacobian we set our currentX to the
  // incoming value and evaluate the RHS.  When the Jacobian is applied,
  // we compute the perturbed residuals and the directional
  // derivative.  Ignore Jac.
  currentX = x;

  bool ok = false;
  if (!useGroupForComputeF)
    ok = interface->computeF(x, fo.getEpetraVector(),
			     NOX::Epetra::Interface::Required::MF_Jac);
  else {
    groupPtr->setX(currentX);
    groupPtr->computeF();
    fo = dynamic_cast<const NOX::Epetra::Vector&>
      (groupPtr->getF());
    ok = true;
  }
#endif

  // Remember the current interface displacements.
  currentX = x;

  // Nothing to do here. The work is done when we apply a vector to
  // the Jacobian.
  bool ok = true;
  return ok;
}

void FSIMatrixFree::setGroupForComputeF(const NOX::Abstract::Group& group)
{
  useGroupForComputeF = true;
  groupPtr = group.clone();
  return;
}
