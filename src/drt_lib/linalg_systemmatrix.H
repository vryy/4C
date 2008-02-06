#ifdef CCADISCRET

#ifndef LINALG_SYSTEMMATRIX_H
#define LINALG_SYSTEMMATRIX_H

#include <vector>

#include <Epetra_CrsMatrix.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_Map.h>
#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>


/*

- Koennten wir Epetra_Operator implementieren?
  Was ist, wenn jemand sowas wie matrixfree operator an die Loeser geben will?

- Platin-style-doxygen Dokumentation notwendig (inkl aller Parameter und 
  Beschreibung, wie das mit der Maske und den Dirichlet-RBs en detail funktioniert)
  
- Umbenennungen: SparseMatrix BaseSparseMatrix BlockSparseMatrix etc.
                 mask_ ->graph_
                 realdirichlet -> irgendwas mit Aussagekraft
                 
- support          A*x
                   A*B
                   A*B*C
                   Transpose
                   Split (the old way)
                   MergeIntoOne (fuer die Loeser)?
                   
- move LINALG functions to class wherever reasonable, remove
  respective LINALG functions.
                   
- change remaining LINALG functions to directly support these new classes

- be much more explicit on treatment of mask, how do I get rid of it if I don't want it?
  How can I explicitly flag the behavior I want
  Be very explicit on whether to use old mask or not
  Do not mix notion of how to handle Dirichlet with notion of how to handle graph
  (the first is a consequence of the second, isn't it?)
  
- ColBlock method is spooky

- For keeping learning curves up and tribute to (important) simplicity:
  Absolutely no templates here
  
- Absolutely no mixing of DRT objects and LINALG algebra objects, 
  which means if you want to have a condition-style system matrix, you can have 
  it in your own sandbox, but not in LINALG (and can not make it be specially
  supported by DRT::Discretization etc.)
  
- give excessive seminar on it so everybody understands it, answer all
  questions by young chair-members in the following.
  (I'm drowning in questions about all sorts of stuff)  

*/


namespace LINALG
{

  ///
  class SystemMatrix
  {
  public:
    virtual ~SystemMatrix();

    virtual void Zero() = 0;

    virtual void Assemble(const Epetra_SerialDenseMatrix& Aele,
                          const std::vector<int>& lm,
                          const std::vector<int>& lmowner)
    {
      Assemble(Aele,lm,lmowner,lm);
    }

    virtual void Assemble(const Epetra_SerialDenseMatrix& Aele,
                          const std::vector<int>& lmrow,
                          const std::vector<int>& lmrowowner,
                          const std::vector<int>& lmcol) = 0;

    virtual void Complete() = 0;
    virtual void Complete(const Epetra_Map& domainmap, const Epetra_Map& rangemap) = 0;

    virtual void ApplyDirichlet(const Teuchos::RCP<Epetra_Vector> dbctoggle) = 0;
  };


  ///
  class SingleSystemMatrix : public SystemMatrix
  {
  public:
    explicit SingleSystemMatrix(bool realdirichlet=true);
    ~SingleSystemMatrix();

    void Setup(const Epetra_Map& rowmap, const int npr);

    Teuchos::RCP<Epetra_CrsMatrix>& Matrix() { return sysmat_; }

    bool Filled() const { return sysmat_->Filled(); }
    const Epetra_Map& RowMap() const { return sysmat_->RowMap(); }
    const Epetra_Map& ColMap() const { return sysmat_->ColMap(); }

    void Zero();
    void Assemble(const Epetra_SerialDenseMatrix& Aele,
                  const std::vector<int>& lmrow,
                  const std::vector<int>& lmrowowner,
                  const std::vector<int>& lmcol);

    void Assemble(double val, int rgid, int cgid);

    void Complete();
    void Complete(const Epetra_Map& domainmap, const Epetra_Map& rangemap);

    void ApplyDirichlet(const Teuchos::RCP<Epetra_Vector> dbctoggle);

  private:
    Teuchos::RCP<Epetra_CrsMatrix> sysmat_;
    bool realdirichlet_;
    Teuchos::RCP<Epetra_CrsGraph> mask_;
  };


  ///
  class BlockSystemMatrixBase : public SystemMatrix
  {
  public:
    BlockSystemMatrixBase(const Epetra_Map& fullrangemap,
                          const Epetra_Map& fulldomainmap,
                          int rows,
                          int cols,
                          std::vector<Epetra_Map> rangemaps,
                          std::vector<Epetra_Map> domainmaps);

    SingleSystemMatrix& Matrix(int r, int c) { return blocks_[r*cols_+c]; }

    void Zero();
    void Complete();
    void Complete(const Epetra_Map& domainmap, const Epetra_Map& rangemap);

    bool Filled() const;

    int Rows() const { return rows_; }
    int Cols() const { return cols_; }

    Epetra_Map& RangeMap(int r) { return rangemaps_[r]; }
    Epetra_Map& DomainMap(int r) { return domainmaps_[r]; }

    Epetra_Map& FullRangeMap() { return fullrangemap_; }
    Epetra_Map& FullDomainMap() { return fulldomainmap_; }

  private:
    Epetra_Map fullrangemap_;
    Epetra_Map fulldomainmap_;
    int rows_;
    int cols_;
    std::vector<Epetra_Map> rangemaps_;
    std::vector<Epetra_Map> domainmaps_;
    std::vector<SingleSystemMatrix> blocks_;
  };


  ///
  template <class Condition>
  class BlockSystemMatrix : public BlockSystemMatrixBase,
                            public Condition
  {
  public:
    BlockSystemMatrix(int rows,
                      int cols,
                      std::vector<Epetra_Map> rowmaps,
                      std::vector<Epetra_Map> colmaps);

    void Assemble(const Epetra_SerialDenseMatrix& Aele,
                  const std::vector<int>& lmrow,
                  const std::vector<int>& lmrowowner,
                  const std::vector<int>& lmcol);

    void ApplyDirichlet(const Teuchos::RCP<Epetra_Vector> dbctoggle);

    void Complete();
  };


  ///
  class DefaultBlockMatrixCondition
  {
  public:
    explicit DefaultBlockMatrixCondition(BlockSystemMatrixBase* mat);

    int RowBlock(int lrow, int rgid);
    int ColBlock(int rblock, int lcol, int cgid);

    void Assemble(double val,
                  int lrow, int rgid, int rblock,
                  int lcol, int cgid, int cblock);

    void Complete();

  private:
    BlockSystemMatrixBase* mat_;
    std::map<int,std::map<int,double> > ghost_;
  };
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <class Condition>
LINALG::BlockSystemMatrix<Condition>::BlockSystemMatrix(int rows,
                                                        int cols,
                                                        std::vector<Epetra_Map> rowmaps,
                                                        std::vector<Epetra_Map> colmaps)
  : BlockSystemMatrixBase(rows,cols,rowmaps,colmaps),
    Condition(this)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <class Condition>
void LINALG::BlockSystemMatrix<Condition>::Assemble(const Epetra_SerialDenseMatrix& Aele,
                                                    const std::vector<int>& lmrow,
                                                    const std::vector<int>& lmrowowner,
                                                    const std::vector<int>& lmcol)
{
  const int lrowdim = (int)lmrow.size();
  const int lcoldim = (int)lmcol.size();

  const int myrank = fullrangemap_.Comm().MyPID();

  // loop rows of local matrix
  for (int lrow=0; lrow<lrowdim; ++lrow)
  {

    // check ownership of row
    if (lmrowowner[lrow]!=myrank)
      continue;

    int rgid = lmrow[lrow];
    int rblock = this->RowBlock(lrow,rgid);

    for (int lcol=0; lcol<lcoldim; ++lcol)
    {
      double val = Aele(lrow,lcol);
      int cgid = lmcol[lcol];
      int cblock = this->ColBlock(rblock,lcol,cgid);

      this->Assemble(val,lrow,rgid,rblock,lcol,cgid,cblock);
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <class Condition>
void LINALG::BlockSystemMatrix<Condition>::Complete()
{
  Condition::Complete();
  BlockSystemMatrixBase::Complete();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <class Condition>
void LINALG::BlockSystemMatrix<Condition>::ApplyDirichlet(const Teuchos::RCP<Epetra_Vector> dbctoggle)
{
}


#endif
#endif
