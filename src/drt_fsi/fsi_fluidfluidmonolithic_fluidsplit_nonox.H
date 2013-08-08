#ifndef FSI_FLUIDFLUIDMONOLITHIC_FLUIDSPLIT_NONOX_H_
#define FSI_FLUIDFLUIDMONOLITHIC_FLUIDSPLIT_NONOX_H_



#include "../drt_inpar/inpar_xfem.H"
//Parent class MonolithicNoNOX
#include "fsi_monolithic_nonox.H"



/********************************************************************/
namespace ADAPTER
{
	class Coupling;
	class FluidFluidFSI;
}

/********************************************************************/

//UTILS for matrix transformations
namespace FSI
{
	namespace UTILS
	{
		class MatrixRowTransform;
		class MatrixColTransform;
		class MatrixRowColTransform;
	}

}

/********************************************************************/


namespace FSI
{
//rk
///Monolithic fluid-fluid FSI with Fluid-Split
/*!
  Background fluid, embedded fluid & structure.
  Split of the (embedded) fluid matrix, while the structure matrix remains
  unmodified.
 */

	class FluidFluidMonolithicFluidSplitNoNOX : public MonolithicNoNOX
	{

	friend class FSI::FSIResultTest;

	public:

	    explicit FluidFluidMonolithicFluidSplitNoNOX(const Epetra_Comm& comm,const Teuchos::ParameterList& timeparams);


	    /// Initialize, read parameters and create global DOF row map
		virtual void SetupSystem();

	    /// Setup composed system block matrix
		virtual void SetupSystemMatrix();

		/// Repeat setup, if maps have changed
		void SetupNewSystem();

		/// Setup composed RHS-vector
		virtual void SetupRHS(Epetra_Vector& f, bool firstcall=false);

		/// Call IntialGuess() on all fields
		virtual void InitialGuess(Teuchos::RCP<Epetra_Vector> ig );

		/// Create map containing the DOFs with Dirichlet-BCs from all fields!
	    virtual Teuchos::RCP<Epetra_Map> CombinedDBCMap();

		/// Newton-Loop
		virtual void Newton();


	protected:
	    /// extract the three field vectors from a given composed vector
	    /*
	      \param x  (i) composed vector that contains all field vectors
	      \param sx (o) structural displacements
	      \param fx (o) fluid velocities and pressure
	      \param ax (o) ale displacements
	     */
	    virtual void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
	                                     Teuchos::RCP<const Epetra_Vector>& sx,
	                                     Teuchos::RCP<const Epetra_Vector>& fx,
	                                     Teuchos::RCP<const Epetra_Vector>& ax);

	    /// Time Update: recover the Lagrange multiplier and relax ALE (where required)
	    virtual void Update();

	    virtual void PrepareTimeStep();

	    /// Compute the Lagrange multiplier (FSI stresses) for the current time step
	    virtual void RecoverLagrangeMultiplier();

	    // Compute the residual and incremental norms required for convergence check
	    virtual void BuildCovergenceNorms();

	    virtual void ReadRestart(int step);

	    virtual void Output();


	private:

	    /// build block vector from field vectors
	    void SetupVector(Epetra_Vector &f,
	                     Teuchos::RCP<const Epetra_Vector> sv,
	                     Teuchos::RCP<const Epetra_Vector> fv,
	                     Teuchos::RCP<const Epetra_Vector> av,
	                     double fluidscale);


	    //F_\Gamma_\Gamma
	    Teuchos::RCP<UTILS::MatrixRowColTransform> fggtransform_;

	    Teuchos::RCP<UTILS::MatrixRowTransform> fmggtransform_;

	    //F_\Gamma_I
	    Teuchos::RCP<UTILS::MatrixRowTransform> fgitransform_;
	    //F_I_\Gamma
	    Teuchos::RCP<UTILS::MatrixColTransform> figtransform_;

	    //F^G_\Gamma_I
	    Teuchos::RCP<UTILS::MatrixRowColTransform> fmgitransform_;
	    //F^G_I_I
	    Teuchos::RCP<UTILS::MatrixColTransform> fmiitransform_;

	    //A_I_\Gamma
	    Teuchos::RCP<UTILS::MatrixColTransform> aigtransform_;


	    //additional ALE residual to avoid incremental ALE errors
	    Teuchos::RCP<Epetra_Vector> aleresidual_;

	    INPAR::FSI::LinearBlockSolver linearsolverstrategy_;

	    Teuchos::RCP<FLD::UTILS::FluidXFluidMapExtractor> xfluidfluidsplitter_;


	    enum INPAR::XFEM::Monolithic_xffsi_Approach monolithic_approach_;

	    //ALE-relaxation algorithm applied?
	    bool relaxing_ale_;

	    //In case of ALE-relaxation, the time step interval
	    //between relaxation has to be specified
	    int relaxing_ale_every_;

	    int currentstep_;

	    //Lagrange multiplier
	    Teuchos::RCP<Epetra_Vector> lambda_;


	    //Fluid matrices: f
	    //Moving fluid matrices: fm
	    //g: \Gamma (FSI interface)
	    //i: I (internal field)
	    //curr: current

	    Teuchos::RCP<const LINALG::SparseMatrix> fggcurr_;
	    Teuchos::RCP<const LINALG::SparseMatrix> fgicurr_;
	    Teuchos::RCP<const LINALG::SparseMatrix> fmggcurr_;
	    Teuchos::RCP<const LINALG::SparseMatrix> fmgicurr_;

	    //Current RHS @ the interface for LM-recovery
	   // Teuchos::RCP<Epetra_Vector> fgcurr_;

	   // not needed here, it is needed for the implementation with nox
	   // to get the iteration increment in ExtractFliedVectors
	   //Inner fluid increment sum
	   //Teuchos::RCP<const Epetra_Vector>  solipre_;
	   //Inner ALE increment sum
	   //Teuchos::RCP<const Epetra_Vector>  solialepre_;
	   //Structural interface displacement increment sum
       //Teuchos::RCP<const Epetra_Vector>  solgpre_;


       //Inner ale displacement increment
	   Teuchos::RCP<Epetra_Vector>  ddialeinc_;
	   //Interface displacement increment
	   Teuchos::RCP<Epetra_Vector>  ddginc_;
	   //Inner fluid velocity increment
	   Teuchos::RCP<Epetra_Vector>  duiinc_;




	};

}


#endif /* FSI_FLUIDFLUIDMONOLITHIC_FLUIDSPLIT_NONOX_H_ */
