
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include <vector>
#include <algorithm>

#include "fsi_utils.H"
#include "../drt_lib/drt_utils.H"

#include <Epetra_CrsMatrix.h>
#include <EpetraExt_RowMatrixOut.h>
#include <NOX.H>
#include <NOX_Epetra.H>
#include <Epetra_SerialDenseMatrix.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

extern "C" /* stuff which is c and is accessed from c++ */
{
#include "../headers/standardtypes.h"
}

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;


void FSI::Utils::DumpJacobian(NOX::Epetra::Interface::Required& interface,
                              double alpha,
                              double beta,
                              Teuchos::RefCountPtr<Epetra_Vector> soln,
                              std::string filename)
{
  // that's really stupid again
  const Epetra_BlockMap& bmap = soln->Map();
  Epetra_Map map(bmap.NumGlobalElements(),
                 bmap.NumMyElements(),
                 bmap.MyGlobalElements(),
                 0,
                 bmap.Comm());

  RefCountPtr<Epetra_CrsMatrix> jacobian = rcp(new Epetra_CrsMatrix(Copy, map, map.NumGlobalElements()));

  int nummyelements = map.NumMyElements();
  int mypos = DRT::Utils::FindMyPos(nummyelements, map.Comm());
  double eta = 0.0;

  Epetra_Vector fo(*soln);
  Epetra_Vector fp(*soln);
  Epetra_Vector Jc(*soln);

  // Compute the RHS at the initial solution
  interface.computeF(*soln, fo, NOX::Epetra::Interface::Required::FD_Res);

  Epetra_Vector x_perturb = *soln;

  for (int i=0; i<map.NumGlobalElements(); ++i)
  {
    if (map.Comm().MyPID()==0)
      cout << "calculate column " << i << "\n";

    int proc = 0;
    int idx = 0;
    if (i>=mypos and i<mypos+nummyelements)
    {
      eta = alpha*(*soln)[i-mypos] + beta;
      x_perturb[i-mypos] += eta;
      idx = map.GID(i-mypos);
      proc = map.Comm().MyPID();
    }

    // Find what proc eta is on
    int broadcastProc = 0;
    map.Comm().SumAll(&proc, &broadcastProc, 1);

    // Send the perturbation variable, eta, to all processors
    map.Comm().Broadcast(&eta, 1, broadcastProc);

    map.Comm().Broadcast(&idx, 1, broadcastProc);

    // Compute the perturbed RHS
    interface.computeF(x_perturb, fp, NOX::Epetra::Interface::Required::FD_Res);

    // Compute the column k of the Jacobian
    Jc.Update(1.0, fp, -1.0, fo, 0.0);
    Jc.Scale( 1.0/eta );

    // Insert nonzero column entries into the jacobian
    for (int j = 0; j < map.NumMyElements(); ++j)
    {
      int gid = map.GID(j);
      if (Jc[j] != 0.0)
      {
        int err = jacobian->SumIntoGlobalValues(gid,1,&Jc[j],&idx);
        if (err>0)
        {
          err = jacobian->InsertGlobalValues(gid,1,&Jc[j],&idx);
        }
        if (err != 0)
          dserror("Assembly failed");
      }
    }

    // Unperturb the solution vector
    x_perturb = *soln;
  }

  jacobian->FillComplete();

  EpetraExt::RowMatrixToMatlabFile(filename.c_str(),*jacobian);
}


void FSI::Utils::MGS(NOX::Epetra::Interface::Required& interface,
                     Teuchos::RefCountPtr<Epetra_Vector> soln,
                     int m)
{

  // Kann man denn den orthonormalen Vektor ins computeF hineingeben?
  // Sollte man statt dessen nicht besser die Jacobian angewendet auf
  // einen Vektor berechnen?
  std::vector<Teuchos::RefCountPtr<Epetra_Vector> > q;

  Epetra_SerialDenseMatrix r(m+1,m);

  Teuchos::RefCountPtr<Epetra_Vector> v = rcp(new Epetra_Vector(*soln));
  q.push_back(v);

  for (int k=0; k<m; ++k)
  {
    Teuchos::RefCountPtr<Epetra_Vector> w = rcp(new Epetra_Vector(*soln));
    interface.computeF(*v, *w, NOX::Epetra::Interface::Required::FD_Res);
    for (int i=0; i<=k; ++i)
    {
      if (w->Dot(*q[i], &r(i,k))) dserror("Dot failed");
      if (w->Update(-r(i,k), *q[i], 1.)) dserror("Update failed");
    }
    if (w->Norm2(&r(k+1,k))) dserror("Norm2 failed");
    if (w->Scale(1./r(k+1,k))) dserror("Scale failed");
    q.push_back(w);
  }

  // Ich habe A nicht vorliegen.

  //[m, n] = size(A);
  //q = zeros(m, n);
  //r = zeros(n, n);
  //for k = 1:n
  //  r(k,k) = norm(A(1:m, k));
  //  if r(k,k) == 0
  //    break;
  //  end
  //  q(1:m, k) = A(1:m, k) / r(k,k);
  //  for j = k+1:n
  //    r(k, j) = dot(q(1:m, k), A(1:m, j));
  //    A(1:m, j) = A(1:m, j) - r(k, j) * q(1:m, k);
  //  end
  //end
  //return [q, r]
}


void FSI::Utils::MGS(const NOX::Abstract::Group& grp,
                     const NOX::Abstract::Vector& dir,
                     int m,
                     std::ostream &stream)
{
  std::vector<Teuchos::RefCountPtr<NOX::Abstract::Vector> > q;
  Epetra_SerialDenseMatrix r(m+1,m);

  Teuchos::RefCountPtr<NOX::Abstract::Vector> v = dir.clone();
  q.push_back(v);

  for (int k=0; k<m; ++k)
  {
    Teuchos::RefCountPtr<NOX::Abstract::Vector> w = computeDirectionalDerivative(*v, grp);

    for (int i=0; i<=k; ++i)
    {
      r(i,k) = w->innerProduct(*q[i]);
      w->update(-r(i,k), *q[i], 1.);
    }
    r(k+1,k) = w->norm();
    if (r(k+1,k)==0)
      break;
    w->scale(1./r(k+1,k));
    q.push_back(w);
  }

  for (unsigned i=0; i<q.size(); ++i)
  {
    q[i]->print(stream);
  }
  stream << r;
}


Teuchos::RefCountPtr<NOX::Abstract::Vector>
FSI::Utils::computeDirectionalDerivative(const NOX::Abstract::Vector& dir,
                                         const NOX::Abstract::Group& grp)
{
  // Allocate space for vecPtr and grpPtr if necessary
  Teuchos::RefCountPtr<NOX::Abstract::Vector> vecPtr = dir.clone(NOX::ShapeCopy);

  // We make a copy of the group and therefore we make sure that the
  // current F in the original group is not lost.
  Teuchos::RefCountPtr<NOX::Abstract::Group> grpPtr = grp.clone(NOX::ShapeCopy);

  // Check that F exists
  if (!grp.isF())
  {
    dserror("Invalid F");
  }

  // Compute the perturbation parameter
  double lambda = 1.0e-6;
  double denominator = dir.norm();

  // Don't divide by zero
  if (denominator == 0.0)
    denominator = 1.0;

  double eta = lambda * (lambda + grp.getX().norm() / denominator);

  // Don't divide by zero
  if (eta == 0.0)
    eta = 1.0e-6;

  // Perturb the solution vector
  vecPtr->update(eta, dir, 1.0, grp.getX(), 0.0);

  // Compute the new F --> F(x + eta * dir)
  grpPtr->setX(*vecPtr);
  grpPtr->computeF();

  // Compute Js = (F(x + eta * dir) - F(x))/eta
  vecPtr->update(-1.0/eta, grp.getF(), 1.0/eta, grpPtr->getF(), 0.0);

  return vecPtr;
}

#endif
#endif
