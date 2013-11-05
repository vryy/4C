/*----------------------------------------------------------------------*/
/*!
\file fpsi_monolithic.h

\brief General framework for monolithic fpsi solution schemes

<pre>
Maintainer: Andreas Rauch
            rauch.andreas@mytum.de
</pre>
*/

/*----------------------------------------------------------------------*/

#ifndef FPSI_MONOLITHIC_H
#define FPSI_MONOLITHIC_H

// FPSI includes
#include "fpsi.H"

// FSI includes
#include "../drt_fsi/fsi_matrixtransform.H"
#include "../drt_fsi/fsi_monolithic.H"

//POROELAST includes
#include "../drt_poroelast/poro_base.H"

//ALE includes
#include "../drt_ale/ale.H"
#include "../drt_ale/ale_lin.H"

//ADAPTER includes
#include "../drt_adapter/ad_fld_base_algorithm.H"

//LINALG includes
#include "../linalg/linalg_mapextractor.H"

//INPAR includes
#include "../drt_inpar/inpar_fpsi.H"

/*----------------------------------------------------------------------*
 | forward declarations                                                  |
 *----------------------------------------------------------------------*/
namespace LINALG
{
  class SparseMatrix;
  class MapExtractor;
}

namespace POROELAST
{
  class Monolithic;
}

namespace ADAPTER
{
  class Coupling;
}

namespace FSI
{
  namespace UTILS
  {
    class DebugWriter;
    class MonolithicDebugWriter;
    class MatrixColTransform;
  }
}
/*----------------------------------------------------------------------*/

namespace FPSI
{

class MonolithicBase : public FPSI_Base
  {
  public:

    // ctor
    explicit MonolithicBase(const Epetra_Comm& comm,
                            const Teuchos::ParameterList& fpsidynparams,
                            const Teuchos::ParameterList& poroelastdynparams);

    // dtor
    virtual ~MonolithicBase();

    // read restart data
    virtual void ReadRestart(int step);

    // start a new time step
    virtual void PrepareTimeStep();

    // take current results for converged and save for next time step
    virtual void Update();

    // calculate stresses, strains, energies
    virtual void PrepareOutput();

    // Output routine accounting for Lagrange multiplier at the interface
    virtual void Output();


    const Teuchos::RCP<POROELAST::Monolithic>& PoroField(){return poroelast_subproblem_;}
    const Teuchos::RCP< ::ADAPTER::Fluid>&     FluidField(){return fluid_subproblem_;}
    const Teuchos::RCP<ALE::Ale>&              AleField(){ return ale_; }

    Teuchos::RCP<DRT::DofSet> porostructdofset;
    Teuchos::RCP<DRT::DofSet> porofluiddofset;
    Teuchos::RCP<DRT::DofSet> fluiddofset;
    Teuchos::RCP<DRT::DofSet> aledofset;

    Teuchos::RCP<std::map<int,int> > Fluid_PoroFluid_InterfaceMap;
    Teuchos::RCP<std::map<int,int> > PoroFluid_Fluid_InterfaceMap;

  protected:

    // underlying poroelast problem
    Teuchos::RCP<POROELAST::Monolithic>             poroelast_subproblem_;
    // underlying fluid of the FSI problem
    Teuchos::RCP< ::ADAPTER::Fluid>                 fluid_subproblem_;
    // underlying ale of the FSI problem
    Teuchos::RCP<ALE::Ale>                          ale_;

    ADAPTER::Coupling& PoroFluidFluidMatchingCoupling() { return *couple_porofluid_fluid_matching_; }
    ADAPTER::Coupling& PoroFluidFluidCoupling()     { return *couple_porofluid_fluid_; }
    ADAPTER::Coupling& StructureFluidCoupling()     { return *coupsf_; }
    ADAPTER::Coupling& SmallStructureFluidCoupling(){ return *smallcoupsf_; }
    ADAPTER::Coupling& FluidAleCoupling()           { return *coupfa_; }
    ADAPTER::Coupling& StructureAleCoupling()       { return *coupsa_; }
    ADAPTER::Coupling& InnerStructureFluidCoupling(){ return *innercoupsf_;}

    const ADAPTER::Coupling& PoroFluidFluidMatchingCoupling()     const { return *couple_porofluid_fluid_matching_; }
    const ADAPTER::Coupling& PoroFluidFluidCoupling()     const { return *couple_porofluid_fluid_; }
    const ADAPTER::Coupling& StructureFluidCoupling()     const { return *coupsf_; }
    const ADAPTER::Coupling& SmallStructureFluidCoupling()const { return *smallcoupsf_; }
    const ADAPTER::Coupling& FluidAleCoupling()           const { return *coupfa_; }
    const ADAPTER::Coupling& StructureAleCoupling()       const { return *coupsa_; }
    const ADAPTER::Coupling& InnerStructureFluidCoupling()const { return *innercoupsf_;}

    virtual Teuchos::RCP<Epetra_Vector> StructToFluid   (Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> FluidToStruct   (Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> SmallFluidToStruct   (Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> SmallStructToFluid   (Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> StructToAle     (Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> AleToStruct     (Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> FluidToAle      (Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> AleToFluid      (Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> FluidToPorofluid(Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> PorofluidToFluid(Teuchos::RCP<const Epetra_Vector> iv) const;

    Teuchos::RCP<ADAPTER::Coupling>                  couple_porofluid_fluid_matching_;
    Teuchos::RCP<ADAPTER::Coupling>                  couple_porofluid_fluid_;
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform>     couplingrowtransform_;
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform>     couplingrowtransform2_;
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform>     couplingrowtransform3_;
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform>     couplingrowtransform4_;
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform>     couplingrowtransform5_;
    Teuchos::RCP<FSI::UTILS::MatrixColTransform>     couplingcoltransform_;
    Teuchos::RCP<FSI::UTILS::MatrixColTransform>     couplingcoltransform2_;
    Teuchos::RCP<FSI::UTILS::MatrixRowColTransform>  couplingrowcoltransform_;
    Teuchos::RCP<FSI::UTILS::MatrixRowColTransform>  couplingrowcoltransform2_;


  private:

    // coupling of structure and fluid at the interface
    Teuchos::RCP<ADAPTER::Coupling>                coupsf_;
    // coupling of structure and fluid at the interface
    Teuchos::RCP<ADAPTER::Coupling>                smallcoupsf_;
    // coupling of structure and ale at the interface
    Teuchos::RCP<ADAPTER::Coupling>                coupsa_;
    // coupling of fluid and ale in the entire fluid volume
    Teuchos::RCP<ADAPTER::Coupling>                coupfa_;
    // coupling of all structure and fluid dofs
    Teuchos::RCP<ADAPTER::Coupling>               innercoupsf_;

  }; // MonolithicBase


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<  MonolithicBase -> Monolithic  >>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


class Monolithic : public MonolithicBase
{

  friend class FSI::UTILS::MonolithicDebugWriter;

public:

  // ctor
  explicit Monolithic(const Epetra_Comm& comm,
                      const Teuchos::ParameterList& fpsidynparams,
                      const Teuchos::ParameterList& poroelastdynparams);

  virtual void SetupSystem();

  void Timeloop();

  void PrepareTimeloop();

  virtual void TimeStep();

  void TestResults(const Epetra_Comm& comm);

  Teuchos::RCP<LINALG::SparseMatrix> FluidPoroCouplingMatrix();

  Teuchos::RCP<LINALG::SparseMatrix> PoroFluidCouplingMatrix();

  // Evaluate mechanical-fluid system matrix
  void ApplyCouplingTerms(  Teuchos::RCP<LINALG::SparseOperator>  k_pf,
                            Teuchos::RCP<LINALG::SparseOperator>  k_fp,
                            Teuchos::RCP<LINALG::SparseMatrix>    p,
                            Teuchos::RCP<LINALG::SparseMatrix>    f,
                            Teuchos::RCP<LINALG::SparseOperator>  a,
                            Teuchos::RCP<LINALG::SparseOperator>  k_fa
                          );

  // Evaluate fluid-mechanical system matrix
  void ApplyFluidCouplMatrix(Teuchos::RCP<LINALG::SparseOperator> k_fp);

  //coupling matrices
  Teuchos::RCP<LINALG::SparseOperator> k_pf_;
  Teuchos::RCP<LINALG::SparseOperator> k_fp_;
  //sub matrices
  Teuchos::RCP<LINALG::SparseOperator> k_pf_struct_;
  Teuchos::RCP<LINALG::SparseOperator> k_pf_porofluid_;
  Teuchos::RCP<LINALG::SparseOperator> k_fp_porofluid_;
  Teuchos::RCP<LINALG::SparseOperator> k_pfs_;
  //ale matrices
  Teuchos::RCP<LINALG::SparseOperator> k_pa_;
  Teuchos::RCP<LINALG::SparseOperator> k_ap_;
  Teuchos::RCP<LINALG::SparseOperator> k_aa_;
  Teuchos::RCP<LINALG::SparseOperator> k_fa_;

  virtual void SetupRHS(bool firstcall=false) = 0;

  virtual void SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat) = 0;
  virtual void SetupSystemMatrix() {SetupSystemMatrix(*SystemMatrix()); }
  virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> SystemMatrix() const = 0;

  // Evaluate all fields at x^n+1 with x^n+1 = x_n + stepinc
  virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> stepinc); ///< increment between time step n and n+1

  void SetupNewton();

  void FPSIFDCheck();

  void LinearSolve();

  bool SetupSolver();

  void BuildConvergenceNorms();

  void PrintNewtonIter();

  void PrintNewtonIterHeader(FILE* ofile);

  void PrintNewtonIterText(FILE* ofile);

  bool Converged();

  // full monolithic dof row map
  Teuchos::RCP<const Epetra_Map> DofRowMap() const { return blockrowdofmap_.FullMap(); }

  // map of all dofs on Dirichlet-Boundary
  virtual Teuchos::RCP<Epetra_Map> CombinedDBCMap();

  // extractor to communicate between full monolithic map and block maps
  const LINALG::MultiMapExtractor& Extractor() const { return blockrowdofmap_; }


protected:

  // block systemmatrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> systemmatrix_;
  // dof row map splitted in (field) blocks
  LINALG::MultiMapExtractor     blockrowdofmap_;
  // dof row map (not splitted)
  Teuchos::RCP<Epetra_Map>      fullmap_;
  // increment between Newton steps k and k+1
  Teuchos::RCP<Epetra_Vector>   iterinc_;
  Teuchos::RCP<Epetra_Vector>   iterincold_;
  // zero vector of full length
  Teuchos::RCP<Epetra_Vector>   zeros_;
  // linear algebraic solver
  Teuchos::RCP<LINALG::Solver>  solver_;
  // rhs of FPSI system
  Teuchos::RCP<Epetra_Vector>   rhs_;
  Teuchos::RCP<Epetra_Vector>   rhsold_;

  Teuchos::RCP<const Epetra_Vector>   meshdispold_;

  Teuchos::RCP<Epetra_Vector>   porointerfacedisplacementsold_;

  // adapt solver tolerancePoroField()->SystemSparseMatrix()
  bool   solveradapttol_;
  int    linesearch_;
  double linesearch_counter;
  double solveradaptolbetter_;

  /// extract the three field vectors from a given composed vector
      /*
        \param x  (i) composed vector that contains all field vectors
        \param sx (o) poroelast dofs
        \param fx (o) free fluid velocities and pressure
        \param ax (o) ale displacements
       */
      virtual void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
                                       Teuchos::RCP<const Epetra_Vector>& sx,
                                       Teuchos::RCP<const Epetra_Vector>& fx,
                                       Teuchos::RCP<const Epetra_Vector>& ax) = 0;

      /// setup list with default parameters
      void SetDefaultParameters(const Teuchos::ParameterList& fpsidynparams);

private:

	enum INPAR::FPSI::ConvergenceNorm normtypeinc_;
	enum INPAR::FPSI::ConvergenceNorm normtypefres_;
	enum INPAR::FPSI::BinaryOp combinedconvergence_;

	double toleranceiterinc_;
	double toleranceresidualforces_;

	int maximumiterations_;
	int minimumiterations_;
	double normofrhs_;
	double normofrhsold_;
	double normofiterinc_;
	double normofiterincold_;

	double normrhsfluidvelocity_;
	double normrhsfluidpressure_;
	double normrhsporofluidvelocity_;
	double normrhsporofluidpressure_;
	double normrhsporointerface_;
	double normrhsfluidinterface_;
	double normrhsporostruct_;
	double normrhsfluid_;
	double normrhsale_;

	 double normofiterincporostruct_;
	 double normofiterincporofluid_;
	 double normofiterincfluid_;
	 double normofiterincporofluidvelocity_;
	 double normofiterincporofluidpressure_;
	 double normofiterincfluidvelocity_;
	 double normofiterincfluidpressure_;
	 double normofiterincale_;
	 double normofiterincfluidinterface_;
	 double normofiterincporointerface_;
	 double normofiterincporo_;

   // iteration step
   int iter_;

   int printscreen_;     // print infos to standard out every printscreen_ steps
   bool printiter_;      // print intermediate iterations during solution
   bool printerrfile_;   // print infos to error file #errfile_
   FILE* errfile_;       // error file handle
   // timer for solution technique
    Epetra_Time timer_;

   bool isfirsttimestep_;

protected:
   bool   islinesearch_;
   // flag is true if this is the first Newton iteration, false otherwise
   bool firstcall_;


}; // class Monolithic

} // namespace FPSI

#endif // FPSI_MONOLITHIC_H
