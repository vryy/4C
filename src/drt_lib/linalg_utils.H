/*!----------------------------------------------------------------------
\file linalg_utils.H
\brief A collection of helper methods for namespace LINALG

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef LINALG_UTILS_H
#define LINALG_UTILS_H


#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_LAPACK.h"
#include "Teuchos_RefCountPtr.hpp"
#include "drt_exporter.H"

using namespace std;
using namespace Teuchos;

/*!
\brief LINALG: namespace of the ccarat linear algebra module

*/
namespace LINALG
{
  /*!
  \brief Create a new Epetra_CrsMatrix and return RefcountPtr to it

  \param rowmap (in): row map of matrix
  \param npr (in): estimated number of entries per row.
                   (need not be exact, better should be too big rather then too small)
  */
  RefCountPtr<Epetra_CrsMatrix> CreateMatrix(const Epetra_Map& rowmap, const int npr);

  /*!
  \brief Create a new Epetra_Vector and return RefcountPtr to it

  \param rowmap (in): row map of vector
  \param init (in): initializa vector to zero upon construction
  */
  RefCountPtr<Epetra_Vector> CreateVector(const Epetra_Map& rowmap, const bool init = true);

  /*!
  \brief Export a vector to a different map

  Values of source are copied to target where maps don't have to match.
  Prerequisite: Either the map of source OR the map of target has to be unique
                (will be tested)
  \warning When source is overlapping (and therefore target is unique), values
           in the overlapping region are inserted into the target on a first come
           first serve basis, meaning they should be equal in the source to
           be deterministic
  \param source (in) : source vector values are taken from
  \param target (out): target vector values will be inserted in
  */
  void Export(const Epetra_Vector& source, Epetra_Vector& target);

  /*!
  \brief Compute all eigenvalues and, optionally, eigenvectors
   of a real symmetric matrix A

  \param A (in/out): Matrix to be analysed, if eigv=true stores eigenvectors
  \param L (in/out): Vector of eigenvalues in ascending order
  \param dim (in) :  Dimension of matrix
  \param eigv (in):  flag to evaluate also eigenvectors
  */
  void SymmetricEigen(Epetra_SerialDenseMatrix& A,
                      Epetra_SerialDenseVector& L,const int dim, const char eigv='N');

  /*!
  \brief Invert a symmetric dim*dim square matrix

  \param A (in/out): Matrix to be inverted
  \param dim (in) :  Dimension of matrix
  */
  void SymmetricInverse(Epetra_SerialDenseMatrix& A, const int dim);

  /*!
  \brief Invert a nonsymmetric dim*dim square matrix

  \param A (in/out): Matrix to be inverted
  \param dim (in) :  Dimension of matrix
  */
  void NonSymmetricInverse(Epetra_SerialDenseMatrix& A, const int dim);

  /*!
  \brief Assemble a Epetra_SerialDenseMatrix into a Epetra_CrsMatrix

  This is an individual call.
  Will only assemble locally and will never do any commmunication.
  All values that can not be assembled locally will be ignored.
  Will use the communicator and rowmap from matrix A to determine ownerships.
  Local matrix Aele has to be square.

  If matrix A is Filled(), it stays so and you can only assemble to
  places already masked. An attempt to assemble into a non-existing
  place is a grave mistake.

  If matrix A is not Filled(), the matrix is enlarged as required.

  \note Assembling to a non-Filled() matrix is much more expensive
        than to a Filled() matrix. If the sparse mask
        does not change it pays to keep the matrix around and assemble
        into the Filled() matrix.

  \param A (out)   : Sparse matrix to be assembled on
  \param Aele (in) : dense matrix to be assembled
  \param lm (in) : vector with gids
  \param lmowner (in) : vector with owner procs of gids
  */
  void Assemble(Epetra_CrsMatrix& A, const Epetra_SerialDenseMatrix& Aele,
                const vector<int>& lm, const vector<int>& lmowner);

  /*!
  \brief Assemble a Epetra_SerialDenseVector into a Epetra_Vector

  This is an individual call.
  Will only assemble locally and will never do any commmunication.
  All values that can not be assembled locally will be ignored.
  Will use the communicator from vector V to determine ownerships.

  \param V (out)   : Sparse matrix to be assembled on
  \param Vele (in) : dense matrix to be assembled
  \param lm (in) : vector with gids
  \param lmowner (in) : vector with owner procs of gids
  */
  void Assemble(Epetra_Vector& V, const Epetra_SerialDenseVector& Vele,
                const vector<int>& lm, const vector<int>& lmowner);


  /*!
  \brief Call FillComplete on a Epetra_CrsMatrix
  */
  void Complete(Epetra_CrsMatrix& A);

  /*!
  \brief Add a (transposed) Epetra_CrsMatrix to another: B = B*scalarB + A(^T)*scalarA

  Add one matrix to another. the matrix B to be added to must not be
  completed. Sparsity patterns of A and B need not match and A and B can be
  nonsymmetric in value and pattern.
  Range, Row and Domain maps of A(^T) and B have to match.

  Note that this is a true parallel add, even in the transposed case!

  \param A          (in)     : Matrix to add to B (must have Filled()==true)
  \param transposeA (in)     : flag indicating whether transposed of A should be used
  \param scalarA    (in)     : scaling factor for A
  \param B          (in/out) : Matrix to be added to (must have Filled()==false)
  \param scalarB    (in)     : scaling factor for B
  */
  void Add(const Epetra_CrsMatrix& A, const bool transposeA, const double scalarA,
           Epetra_CrsMatrix& B, const double scalarB);

  /*!
  \brief Apply dirichlet boundary condition to a linear system of equations

  Modifies a system of equations such that dirichlet boundary conditions are enforced.
  Prescribed dirichlet BC values are supplied in dbcval and dbctoggle, where
  a prescribed value is dbcval[i] and dbctoggle[i] = 1.0. No BC is enforced in
  all places where dbctoggle[i] != 1.0.<br>
  Let us denote the \f$ A_{2 \times 2} \f$ blocks of \f$A\f$ by
  \f$A_{ff}, A_{fD}, A_{Df}, A_{DD}\f$, where \f$f\f$ stands for 'free' and
  \f$D\f$ stands for 'Dirichlet BC'. Then, after a call to this method<br>

  \f$ A_{ff} = A_{ff}, \f$<br>
  \f$ A_{fD} = A_{fD}, \f$<br>
  \f$ A_{Df} = 0_{Df}, \f$<br>
  \f$ A_{DD} = I_{DD}, \f$<br>
  \f$ x_{D} = dbcval_{D}, \f$<br>
  \f$ b_{D} = dbcval_{D} \f$<br>

  and

  \f$ A_{ff} x_f + A_{fD} x_D = b_f \f$<br>
  \f$ 0 x_f + I_{DD} x_D = x_D \f$.<br>

  \note The matrix is then nonsymmetric. When using iterative methods on this
        linear system of equations that depend on the symmetry of the matrix (such as e.g. CG),
        the initial guess supplied to the solver has to be exact at the
        Dirichlet BCs. This should be easy, as the values at the Dirichlet BCs
        are known.

  \note The mask of matrix \f$A\f$ is not modified. That is the
        entries in \f$A_{Df}\f$ and \f$A_{DD}\f$ are set to zero, not
        removed. This way the matrix can be reused in the next step.

  \param A         (in/out) : Matrix of Ax=b
  \param x         (in/out) : initial guess vector x of Ax=b
  \param b         (in/out) : rhs vector b of Ax=b
  \param dbcval    (in)     : vector holding prescribed dirichlet values
  \param dbctoggle (in)     : vector holding 1.0 where dirichlet should be applied
                              and 0.0 everywhere else
  */
  void ApplyDirichlettoSystem(RefCountPtr<Epetra_CrsMatrix>&   A,
                              RefCountPtr<Epetra_Vector>&      x,
                              RefCountPtr<Epetra_Vector>&      b,
                              const RefCountPtr<Epetra_Vector> dbcval,
                              const RefCountPtr<Epetra_Vector> dbctoggle);

  /*!
  \brief Apply dirichlet boundary condition to a linear system of equations


  \param A         (in/out) : Matrix to of Ax=b
  \param dbctoggle (in)     : vector holding 1.0 where dirichlet should be applied
                              and 0.0 everywhere else
  */
  void ApplyDirichlettoSystem(RefCountPtr<Epetra_CrsMatrix>&   A,
                              const RefCountPtr<Epetra_Vector> dbctoggle);

  /*!
  \brief Apply dirichlet boundary condition to a linear system of equations


  \param x (in/out)         : vector x of Ax=b
  \param b (in/out)         : vector b of Ax=b
  \param dbcval (in)        : vector holding values that are supposed to be prescribed
  \param dbctoggle (in)     : vector holding 1.0 where dirichlet should be applied
                              and 0.0 everywhere else
  */
  void ApplyDirichlettoSystem(RefCountPtr<Epetra_Vector>&      x,
                              RefCountPtr<Epetra_Vector>&      b,
                              const RefCountPtr<Epetra_Vector> dbcval,
                              const RefCountPtr<Epetra_Vector> dbctoggle);




  /*!
  \brief Gather information of type vector<T> on a subset of processors

  This template gathers information provided in sdata on a subset of
  processors tprocs, where the length of the array tprocs is ntargetprocs.
  The redistributed data is returned in rdata which has appropiate size
  on output (size of rdata on input is arbitrary). ntargetprocs can be
  one to reduce data to one proc, it also can be equal to the total number
  of processors to make sdata redundant on all procs.

  \note Functionality of this method is equal to that of Epetra_Comm::GatherAll
        except for that the Epetra version demands the data to be of constant
        size over all procs which this method does not require!

  \param sdata (in) : Information to be gathered on tprocs.
                      Length of sdata can be different on every proc.
  \param rdata (out): Information from sdata gathered on a subset of procs.
                      size of rdata on input is arbitrary, it is exact on output.
  \param ntargetprocs (in): length of tprocs
  \param tprocs (in): vector of procs ids the information in sdata shall be
                      gathered on.
  \param comm (in):   communicator to be used.


  */
  template<typename T> void Gather(vector<T>&         sdata,
                                   vector<T>&         rdata,
                                   const int          ntargetprocs,
                                   const int*         tprocs,
                                   const Epetra_Comm& comm)
  {
    const int myrank  = comm.MyPID();
    const int numproc = comm.NumProc();
    if (numproc==1)
    {
      rdata = sdata;
      return; // nothing to do in serial
    }
    // build a map of data
    map<int,vector<T> > datamap;
    datamap[myrank] = sdata;
    // build a source map
    Epetra_Map source(numproc,1,&myrank,0,comm);
    // build a target map which is redundant on all target procs and zero everywhere else
    bool iamtarget = false;
    for (int i=0; i<ntargetprocs; ++i)
      if (tprocs[i]==myrank)
      {
        iamtarget = true;
        break;
      }
    vector<int> targetvec(0);
    if (iamtarget)
    {
      targetvec.resize(numproc);
      for (int i=0; i<numproc; ++i) targetvec[i] = i;
    }
    const int tnummyelements = (int)targetvec.size();
    Epetra_Map target(-1,tnummyelements,&targetvec[0],0,comm);
    // build an exporter and export data
    DRT::Exporter exporter(source,target,comm);
    exporter.Export(datamap);
    // put data from map in rdata
    rdata.clear();
    int count=0;
    map<int,vector<int> >::iterator curr;
    for (curr=datamap.begin(); curr != datamap.end(); ++curr)
    {
      vector<T>& current = curr->second;
      const int size = (int)current.size();
      rdata.resize((int)rdata.size()+size);
      for (int i=0; i<size; ++i)
      rdata[count+i] = current[i];
      count += size;
    }
    return;
  }


} // namespace LINALG


#endif  // #ifndef LINALG_UTILS_H
#endif  // #ifdef CCADISCRET
