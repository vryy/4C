/*!----------------------------------------------------------------------
\file linalg_downwindmatrix.H

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef LINALG_DOWNWINDMATRIX_H
#define LINALG_DOWNWINDMATRIX_H

// Trilinos includes
#include <Epetra_CrsMatrix.h>
#include <EpetraExt_Reindex_CrsMatrix.h>
#include <EpetraExt_Reindex_LinearProblem.h>
#include <EpetraExt_Reindex_MultiVector.h>
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Time.h"
#include "Epetra_IntVector.h"

#include "drt_dserror.H"
#include "linalg_mapextractor.H"
#include "linalg_systemmatrix.H"
#include "linalg_ana.H"

using namespace std;
using namespace Teuchos;

namespace LINALG
{

/*!
\author gee (gee@lnm.mw.tum.de)

\brief Creates a convection downwind numbering of a linear problem
*/
class  DownwindMatrix
{
public:

  /*!
  \brief Standard Constructor
  
  \param A (in) : Matrix to create downwind numbering for
  \param nv (in): number of velocity degrees of freedom per node, may not be zero
  \param np (in): number of pressure degrees of freedom per node, may be zero
  \param tau (in): threshold parameter, somewhere between 1.0 and about 4.0
  \param outlevel (in): level of output to screen (0 no output, >0 output)
  */
  explicit DownwindMatrix(RCP<Epetra_CrsMatrix> A, const int nv, 
                          const int np, const double tau, const int outlevel);

  /*!
  \brief Destructor
  */
  virtual ~DownwindMatrix() {}

  /*!
  \brief return reindexed linear problem
  
  Takes a linear problem and reindexes it. 
  The lp provided must contain a Epetra_CrsMatrix, a rhs and lhs MultiVector.
  Returns a ref to the reindexed problem.
  Note that this class remains owner of the reindexed problem and destroys it
  in its ctor.
  
  \param lp (in): original linear problem
  \return reference to reindexed linear problem
  */
  Epetra_LinearProblem& ReindexedProblem(Epetra_LinearProblem& lp)
  {
    // shallow redinex the problem, always create new reindexer
    lpreindex_  = rcp(new EpetraExt::LinearProblem_Reindex(ndofrowmap_.get()));
    Epetra_LinearProblem& reproblem = (*lpreindex_)(lp);
    
    // Allocate sorted reindex linear problem
    sreproblem_ = rcp(new Epetra_LinearProblem());
    
    // get matrix, lhs and rhs
    Epetra_CrsMatrix* matrix = dynamic_cast<Epetra_CrsMatrix*>(reproblem.GetMatrix());
    Epetra_MultiVector* lhs = reproblem.GetLHS();
    Epetra_MultiVector* rhs = reproblem.GetRHS();
    if (!matrix || !lhs || !rhs) dserror("linear problem broken");
    
    Asre_ = rcp(new Epetra_CrsMatrix(Copy,*sndofrowmap_,matrix->MaxNumEntries()));
    Asre_->Export(*matrix,*sexporter_,Insert);
    Asre_->FillComplete();
        
    // export lhs,rhs to new layout
    slhs_ = rcp(new Epetra_Vector(*sndofrowmap_));
    slhs_->Export(*lhs,*sexporter_,Insert);
    srhs_ = rcp(new Epetra_Vector(*sndofrowmap_));
    srhs_->Export(*rhs,*sexporter_,Insert);
    
    // put stuff in linear problem
    sreproblem_->SetOperator(Asre_.get());
    sreproblem_->SetLHS(slhs_.get());
    sreproblem_->SetRHS(srhs_.get());
    
    return *sreproblem_;
  }
  
  /*!
  \brief return reindexed matrix
  
  Returns a ptr to the reindexed matrix
  Note that this class remains owner of the reindexed matrix and destroys it
  in its ctor. The calling method must NOT delete the matrix.
  
  \return ptr to reindexed matrix
  */
  inline Epetra_CrsMatrix* ReindexedMatrix()
  {
    if (Asre_ != null) return Asre_.get();
    else dserror("Do not have downwind matrix");
    return NULL;
  }
  
  /*!
  \brief undo reindexing of lhs vector
  
  Returns values from the reindexed lhs vector into a supplied lhs vector with
  original ordering.
  
  \param lhs (out): lhs vector in original ordering to be filled with values from
                    the reindexed lhs
  */
  inline void UndoDownWindLHS(Epetra_MultiVector* lhs)
  {
    // shallow reindex the vector
    EpetraExt::MultiVector_Reindex reindex(*ndofrowmap_);
    Epetra_MultiVector& relhs = reindex(*lhs);
    // export back to the original vector
    relhs.Export(*(sreproblem_->GetLHS()),*rexporter_,Insert);
  }
  
  /*!
  \brief undo reindexing of rhs vector
  
  Returns values from the reindexed rhs vector into a supplied rhs vector with
  original ordering.
  
  \param rhs (out): rhs vector in original ordering to be filled with values from
                    the reindexed rhs
  */
  inline void UndoDownWindRHS(Epetra_MultiVector* rhs)
  {
    // shallow reindex the vector
    EpetraExt::MultiVector_Reindex reindex(*ndofrowmap_);
    Epetra_MultiVector& rerhs = reindex(*rhs);
    // export back to the original vector
    rerhs.Export(*(sreproblem_->GetRHS()),*rexporter_,Insert);
  }
  
  /*!
  \brief clear temporary data
  
  */
  inline void ClearData()
  {
    A_ = null;
    lpreindex_ = null;
    sreproblem_ = null;
    Asre_ = null;
    slhs_ = null;
    srhs_ = null;
  }

private:

  // don't want copy-ctor and = operator
  DownwindMatrix(DownwindMatrix& old);
  DownwindMatrix operator = (const DownwindMatrix& old);

  /// setup phase of downwinding
  void Setup();
  
  /// sets an index in the Bey&Wittum method
  void SetF(const int i, int& nf, Epetra_IntVector& index, 
            const Epetra_CrsMatrix& graph, int rec);

  /// sets an index in the Hackbusch method
  void SetL(const int i, int& nf, Epetra_IntVector& index, 
            const Epetra_CrsMatrix& graph,int rec);

  /// Do downwinding according to Bey & Wittum
  RCP<Epetra_Map> DownwindBeyWittum(const Epetra_CrsMatrix& nnodegraph, 
                                    const Epetra_IntVector& oninflow);
  
  /// Do downwinding according to Hackbusch
  RCP<Epetra_Map> DownwindHackbusch(const Epetra_CrsMatrix& nnodegraph, 
                                    const Epetra_IntVector& oninflow);

  /// test whether k is successor of i
  inline bool IsSuccessor(const int k,const int gi, 
                          const Epetra_CrsMatrix& graph)
  {
    int numentries;
    int* indices;
    double* values;
    graph.ExtractMyRowView(k,numentries,values,indices);
    for (int j=0; j<numentries; ++j)
      if (graph.ColMap().GID(indices[j])==gi) return true;
    return false;
  }

  /// get communicator (A_ not present all the time!)
  const Epetra_Comm& Comm() const {return(A_->Comm());}

  RCP<Epetra_CrsMatrix>                 A_;           // original input matrix
  const int                             outlevel_;    // level of output
  const int                             nv_;          // number velocity dofs
  const int                             np_;          // number pressure dofs
  const int                             bs_;          // total nodal block size (nv_ + np_)
  const double                          tau_;         // strong connections cutoff ratio

  RCP<EpetraExt::LinearProblem_Reindex> lpreindex_;   // linear problem shallow reindexer

  RCP<Epetra_LinearProblem>             sreproblem_;  // sorted reindexed linear problem

  RCP<Epetra_Map>                       ndofrowmap_;  // result new row map of system
  RCP<Epetra_Map>                       sndofrowmap_; // result new row map of system with dofs sorted in ascending order

  RCP<Epetra_Export>                    sexporter_;   // export from ndofrowmap_ to sndofrowmap_
  RCP<Epetra_Export>                    rexporter_;   // export the other direction

  RCP<Epetra_CrsMatrix>                 Asre_;        // deeply reindexed and sorted matrix
  RCP<Epetra_Vector>                    slhs_;        // deeply reindexed vectors
  RCP<Epetra_Vector>                    srhs_;

  
}; // class  DownwindMatrix
}  // namespace LINALG


#endif // LINALG_DOWNWINDMATRIX_H
#endif  // #ifdef CCADISCRET
