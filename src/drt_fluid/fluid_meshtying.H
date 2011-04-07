/*!----------------------------------------------------------------------
\file fluidrmeshtying.H

\brief Methods to apply meshtying to fluid and scatra systems

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

 *----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef FLUID_MESHTYING_H_
#define FLUID_MESHTYING_H_

#include <Teuchos_RCP.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>

// get rid of it ???
#include "fluid_utils_mapextractor.H"
#include "../drt_adapter/adapter_coupling_mortar.H"
#include "../drt_inpar/inpar_fluid.H"

using namespace std;
using namespace Teuchos;

//forward declarations
namespace DRT
{
  class Discretization;
}
namespace Adapter
{
  class CouplingMortar;
}

namespace LINALG
{
  class Solver;
  class SparseMatrix;
  class MultiMapExtractor;
  class MapExtractor;
  class BlockSparseMatrixBase;
  class SparseOperator;
}

/*!
\brief  This class manages local transformations(rotation) of velocity fields

        o used for rotationally symmetric boundary conditions
        o used for rotationally symmetric boundary conditions

\author ehrl

 */

namespace FLD
{
class Meshtying
{
  public:

    Meshtying(RCP<DRT::Discretization>      dis,
              int       meshtyingoption);           ///> actual discretisation


    virtual ~Meshtying(){return;};

    RCP<Epetra_Map> MergedMap(){return mergedmap_;};

    RCP<LINALG::SparseOperator> Setup();

    //! matrix condensation
    void Condensation(RCP<LINALG::SparseOperator>&  sysmat,
                      RCP<Epetra_Vector>&           residual);

    void SolveMeshtying(LINALG::Solver&                 solver,
                        RCP<LINALG::SparseOperator>     sysmat,
                        RCP<Epetra_Vector>&              incvel,
                        RCP<Epetra_Vector>              residual,
                        int                             itnum,
                        RCP<Epetra_MultiVector>         w,
                        RCP<Epetra_MultiVector>         c,
                        bool                            project
                        );

    void KrylovProjection(RCP<Epetra_Vector>   c);

   void ResidualSaddlePointSystem(
       RCP<Epetra_Vector>      feff);

   void UpdateLag();


    void OutputSetUp();

    void OutputSparseMatrixSplit(
        RCP<LINALG::SparseOperator>                conmat
        );

    void OutputBlockMatrix(RCP<LINALG::BlockSparseMatrixBase>      blockmatrix,
                           RCP<Epetra_Vector>                      residual);

    void OutputVectorSplit(
        RCP<Epetra_Vector>                conres,
        std::vector<RCP<Epetra_Vector> >&  vectorsplit);


  private:

    void CondensationBlockMatrix(RCP<LINALG::SparseOperator>&   sysmat,
                                 RCP<Epetra_Vector>&          residual);

    void CondensationSparseMatrix(
        RCP<LINALG::SparseOperator>&   sysmat,
        RCP<Epetra_Vector>&            residual);

    void SplitMatrix(RCP<LINALG::SparseOperator>                 matrix,
                     std::vector<RCP<LINALG::SparseMatrix> >&   splitmatrix);

    void SplitVector(RCP<Epetra_Vector>                 vector,
                     std::vector<RCP<Epetra_Vector> >&   splitvector);

    void CondensationOperationSparseMatrix(
        RCP<LINALG::SparseOperator>&                 sysmat,
        RCP<Epetra_Vector>&                        residual,
        std::vector<RCP<LINALG::SparseMatrix> > & splitmatrix,
        std::vector<RCP<Epetra_Vector> >&         splitvector);

    void CondensationOperationBlockMatrix(
        RCP<LINALG::SparseOperator>&      sysmat,
        RCP<Epetra_Vector>&                      residual,
        std::vector<RCP<Epetra_Vector> >&        splitvector);

    void UpdateSlaveDOF(RCP<Epetra_Vector>&   inc);


    void PrepareSaddlePointSystem(
        RCP<LINALG::SparseOperator>    sysmat,
        RCP<LINALG::SparseMatrix>      mergedsysmat,
        RCP<Epetra_Vector>             residual,
        RCP<Epetra_Vector>             mergedresidual);

    void PrepareSaddlePointSystemPC(
        RCP<LINALG::SparseOperator>        sysmat,
        RCP<LINALG::BlockSparseMatrixBase> blocksysmat,
        RCP<Epetra_Vector>                 residual,
        RCP<Epetra_Vector>                 mergedresidual);

    void UpdateSaddlePointSystem(
        RCP<Epetra_Vector>      inc,
        RCP<Epetra_Vector>      mergedinc);

  private:

    //
    RCP<DRT::Discretization> discret_;

    //
    int msht_;

    //! the processor ID from the communicator
    int                           myrank_;

    //
    ADAPTER::CouplingMortar adaptermeshtying_;

    RCP<Epetra_Map> problemrowmap_;

    const Epetra_Map* dofrowmap_;

    // slave dof rowmap
    RCP<Epetra_Map> gndofrowmap_;

    // slave & master dof rowmap
    RCP<Epetra_Map> gsmdofrowmap_;

    // slave dof rowmap
    RCP<Epetra_Map> gsdofrowmap_;

    // master dof rowmap
    RCP<Epetra_Map> gmdofrowmap_;

    // master dof rowmap
    RCP<Epetra_Map> glmdofrowmap_;

    // map extractor mesh tying
    LINALG::MultiMapExtractor extractor_;

    // lagrange multiplier
    RCP<Epetra_Vector> lag_;

    // lagrange multiplier
    RCP<Epetra_Vector> lagold_;

    double theta_;

    RCP<Epetra_Map> mergedmap_;

};  // end  class Meshtying
} // end namespace FLD

#endif   /* FLUID_MESHTYING_H_ */
#endif  // #ifdef CCADISCRET
