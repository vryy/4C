/*----------------------------------------------------------------------*/
/*!
 \file poro_monolithicstructuresplit.H

 \brief  monolithic structuresplit poroelasticity algorithms

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef PORO_MONOLITHICFLUIDSPLIT_H_
#define PORO_MONOLITHICFLUIDSPLIT_H_

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/

#include "poroelast_monolithic.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/

namespace FSI
{
  //class OverlappingBlockMatrix;

  namespace UTILS
  {
    class MatrixRowTransform;
    class MatrixColTransform;
    class MatrixRowColTransform;
  }
}

namespace POROELAST
{
  class MonolithicFluidSplit: public Monolithic
  {

    public:
      //! create using a Epetra_Comm
      explicit MonolithicFluidSplit(const Epetra_Comm& comm,
                                    const Teuchos::ParameterList& timeparams);

      /*! do the setup for the monolithic system

       <pre>

       1.) setup coupling
       2.) get maps for all blocks in the system (and for the whole system as well)
       create combined map
       3.) create system matrix

       </pre>

       \note We want to do this setup after reading the restart information, not
       directly in the constructor. This is necessary since during restart (if
       ReadMesh is called), the dofmaps for the blocks might get invalid.
       */
      // Setup the monolithic Poroelasticity system
      virtual void SetupSystem();

      /// setup composed right hand side from field solvers
      virtual void SetupRHS(bool firstcall=false);

      /// setup composed system matrix from field solvers
      virtual void SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat
          );

      //! start a new time step
      virtual void PrepareTimeStep();

    private:

      /// build block vector from field vectors
      void SetupVector(Epetra_Vector &f,
                       Teuchos::RCP<const Epetra_Vector> sv,
                       Teuchos::RCP<const Epetra_Vector> fv,
                       double fluidscale);

      //! extract the field vectors from a given composed vector
      /*!
       \param x  (i) composed vector that contains all field vectors
       \param sx (o) structural vector (e.g. displacements)
       \param fx (o) fluid vector (e.g. velocities and pressure)
       */
      virtual void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
                                       Teuchos::RCP<const Epetra_Vector>& sx,
                                       Teuchos::RCP<const Epetra_Vector>& fx,
                                       bool firstcall=false);

      //! combined DBC map
      //! unique map of all dofs that should be constrained with DBC
      virtual Teuchos::RCP<Epetra_Map> CombinedDBCMap();

      /// recover Lagrange multiplier \f$\lambda_\Gamma\f$ at the interface at the end of each time step
      /// (i.e. condensed forces onto the structure) needed for rhs in next time step
      virtual void RecoverLagrangeMultiplier();

      Teuchos::RCP<FSI::UTILS::MatrixRowColTransform> fggtransform_;
      Teuchos::RCP<FSI::UTILS::MatrixRowTransform> fgitransform_;
      Teuchos::RCP<FSI::UTILS::MatrixColTransform> figtransform_;
      Teuchos::RCP<FSI::UTILS::MatrixRowTransform> cfggtransform_;
      Teuchos::RCP<FSI::UTILS::MatrixColTransform> csggtransform_;
      Teuchos::RCP<FSI::UTILS::MatrixRowTransform> cfgitransform_;
      Teuchos::RCP<FSI::UTILS::MatrixColTransform> csigtransform_;

      //! flag indicating whether there are no slip conditions to be evaluated at the interface
      bool evaluateinterface_;

      //! map containing DOFs with fsi- and DBC conditions
      Teuchos::RCP<Epetra_Map> fsibcmap_;

      /// @name Some quantities to recover the Langrange multiplier at the end of each time step

      //! Lagrange multiplier \f$\lambda_\Gamma^n\f$ at the interface (ie condensed forces onto the structure)
      //! evalutated at old time step \f$t_n\f$ but needed for next time step \f$t_{n+1}\f$
      Teuchos::RCP<Epetra_Vector> lambda_;

      //! interface force \f$f_{\Gamma,i+1}^{S,n+1}\f$ onto the fluid at current iteration \f$i+1\f$
      Teuchos::RCP<const Epetra_Vector> fgcur_;

      //! interface force \f$f_{\Gamma,i}^{S,n+1}\f$ onto the fluid at previous iteration \f$i\f$
      //Teuchos::RCP<const Epetra_Vector> fgpre_;

      //! inner structural displacement increment \f$\Delta(\Delta d_{I,i+1}^{n+1})\f$ at current iteration \f$i+1\f$
      Teuchos::RCP<Epetra_Vector> ddiinc_;

      //! inner fluid velocity increment \f$\Delta(\Delta u_{I,i+1}^{n+1})\f$ at current iteration \f$i+1\f$
      Teuchos::RCP<Epetra_Vector> duiinc_;

      //! inner displacement solution of the structure at previous iteration
      Teuchos::RCP<const Epetra_Vector> solipre_;

      //! inner velocity/pressure solution of the fluid at previous iteration
      Teuchos::RCP<const Epetra_Vector> solivelpre_;

      //! structural interface displacement increment \f$\Delta(\Delta d_{\Gamma,i+1}^{n+1})\f$ at current iteration \f$i+1\f$
      Teuchos::RCP<Epetra_Vector> ddginc_;

      //! fluid interface velocity increment \f$\Delta(\Delta d_{\Gamma,i+1}^{n+1})\f$ at current iteration \f$i+1\f$
      Teuchos::RCP<Epetra_Vector> duginc_;

      //! interface displacement solution of the structure at previous iteration
      Teuchos::RCP<const Epetra_Vector> solgpre_;

      //! interface velocity solution of the fluid at previous iteration
      Teuchos::RCP<const Epetra_Vector> solgvelpre_;

      //! block \f$F_{\Gamma I,i+1}\f$ of fluid matrix at current iteration \f$i+1\f$
      Teuchos::RCP<const LINALG::SparseOperator> fgicur_;

      //! block \f$F_{\Gamma I,i}\f$ of fluid matrix at previous iteration \f$i\f$
     // Teuchos::RCP<const LINALG::SparseOperator> fgipre_;

      //! block \f$F_{\Gamma\Gamma,i+1}\f$ of fluid matrix at current iteration \f$i+1\f$
      Teuchos::RCP<const LINALG::SparseOperator> fggcur_;

      //! block \f$F_{\Gamma\Gamma,i}\f$ of fluid matrix at previous iteration \f$i\f$
      //Teuchos::RCP<const LINALG::SparseOperator> fggpre_;

      //! block \f$C_{\Gamma\Gamma,i+1}\f$ of fluid matrix at current iteration \f$i+1\f$
      Teuchos::RCP<const LINALG::SparseOperator> cgicur_;

      //! block \f$C_{\Gamma\Gamma,i}\f$ of fluid matrix at previous iteration \f$i\f$
      //Teuchos::RCP<const LINALG::SparseOperator> cgipre_;

      //! block \f$C_{\Gamma\Gamma,i+1}\f$ of fluid matrix at current iteration \f$i+1\f$
      Teuchos::RCP<const LINALG::SparseOperator> cggcur_;

      //! block \f$C_{\Gamma\Gamma,i}\f$ of fluid matrix at previous iteration \f$i\f$
      //Teuchos::RCP<const LINALG::SparseOperator> cggpre_;

      //@}

  }; // Monolithic

} // namespace POROELAST

#endif /* PORO_MONOLITHICFLUIDSPLIT_H_ */
