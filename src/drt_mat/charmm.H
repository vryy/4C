/*!----------------------------------------------------------------------
\file charmm.H
\brief CHARMm Interface to compute the mechanical properties of integrins

<pre>
Maintainer: Robert Metzke
	    metzke@lnm.mw.tum.de
	    http://www.lnm.mw.tum.de
	    089 - 289-15244
</pre>
 *----------------------------------------------------------------------*/
#ifndef CHARMM_H_
#define CHARMM_H_

#include "material.H"
#include "matpar_parameter.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"


namespace MAT {

    namespace PAR {

	/*----------------------------------------------------------------------*/
	/// material parameters for Integrin

	class CHARMM : public Parameter {
	public:

	    /// standard constructor
	    CHARMM(
		    Teuchos::RCP<MAT::PAR::Material> matdata
		    );

	    /// destructor

	    virtual ~CHARMM() {
		;
	    }

	    /// @name material parameters
	    //@{

	    /// \b ORIGIN: boolean (1/0) to decide if F is computed at the element origin
	    const int origin_;
	    /// \b FCL: first characteristic (fc) length (double) [AKMA]
	    const double fcl_;
	    /// \b FCD_Type type of the fc direction (string) (principal/vector/none)
	    const string* fcd_type_;
	    /// \b FCD: direction of fc direction if fcd_type_ is vector (unit vector)
	    const vector<double>* fcd_;
	    /// \b FCD_Space: directional space of fc (unit vector)
	    const vector<double>* fcds_;

	    /// \b SCL: second characteristic (sc) length [AKMA]
	    const double scl_;
	    /// \b SCD_Type: type of the sc direction (string) (principal/vector/none)
	    const string* scd_type_;
	    /// \b SCD: direction of sc direction if scd_type_ is vector (unit vector)
	    const vector<double>* scd_;
	    /// \b SCD_Space: directional space of sc (unit vector)
	    const vector<double>* scds_;

	    /// \b FCD_Acceleration: boolean (1/0) if F'' in fcd is computed and as force given to MD
	    const int fcdacc_;
	    /// \b AtomicMass: Atomic mass of the moving part of the protein [AKMA]
	    const double atomicmass_;
	    /// \b Facc_Scale: Scale factor from whatever force unit used in FE to pN
	    const double facc_scale_;
	    /// \b Time_AKMA: Scale factor from FE delta t time to AKMA time
	    const double timeakma_;
	    /// \b Time_Scale: Linear scaling between FE and MD delta t time
	    const double timescale_;

	    /// \b HARD: Use hard coded results e.g. for testing or external run MD
	    const int hard_;
	    /// \b c_Scale: Scale factor for material parameter c1 (Neo-Hookean) depending on used FE units
	    const double c_scale_;

	    /// \b PATH: Path to the folder structre for CHARMm problem type
	    const string* path_;
            /// \b USE_OLD_RSULTS: Use previously computed results from CHARMm
	    const int use_old_results_;
            /// \b SERPAR: Serial or parallel CHARMm execution
	    const string* serpar_;
            /// \b CHARMM: CHARMm executable
	    const string* charmm_;
            /// \b INPUT: CHARMm input file
            const string* input_;

            /// \b NUE: poisson ratio
	    const double nue_;
	    /// \b DENS: mass density
	    const double density_;

	    //@}

            /// create material instance of matching type with my parameters
            virtual Teuchos::RCP<MAT::Material> CreateMaterial();

	}; // class CHARMM

    } // namespace PAR

class CHARMMType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "CHARMMType"; }

  static CHARMMType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static CHARMMType instance_;
};

    /*----------------------------------------------------------------------*/

    /*!
    \brief CHARMm Interface to compute the mechanical properties of integrins

     The material is a interface to the molecular dynamics code CHARMm.
     The idea of ACE (Atomistic to Continuum Energy based) coupling is that,
     FE is giving strains to MD, which is using this strains to compute an
     energy difference between the both configurations. This energy difference
     is then given back to FE and used to compute the material properties of
     the element.

     References:

     CHARMM: A Program for Macromolecular Energy, Minimization, and Dynamics
     Calculations, J. Comp. Chem. 4, 187-217 (1983), by B. R. Brooks,
     R. E. Bruccoleri, B. D. Olafson, D. J. States, S. Swaminathan, and
     M. Karplus.

     CHARMM: The Energy Function and Its Parameterization with an Overview of
     the Program, in The Encyclopedia of Computational Chemistry, 1, 271-277,
     P. v. R. Schleyer et al., editors (John Wiley & Sons: Chichester, 1998),
     by A. D. MacKerell, Jr., B. Brooks,C. L. Brooks, III, L. Nilsson, B. Roux,
     Y. Won, and M. Karplus.

    */

    class CHARMM : public Material {
    public:

	/// construct empty material object
	CHARMM();

	/// construct the material object given material parameters
	explicit CHARMM(MAT::PAR::CHARMM* params);

	//! @name Packing and Unpacking

	/*!
	  \brief Return unique ParObject id

	  every class implementing ParObject needs a unique id defined at the
	  top of drt_parobject.H (this file) and should return it in this method.
	 */
	virtual int UniqueParObjectId() const {
          return CHARMMType::Instance().UniqueParObjectId();
	}

	/*!
	  \brief Pack this class so it can be communicated

	  Resizes the vector data and stores all information of a class in it.
	  The first information to be stored in data has to be the
	  unique parobject id delivered by UniqueParObjectId() which will then
	  identify the exact class on the receiving processor.

	  \param data (in/out): char vector to store class information
	 */
	virtual void Pack(DRT::PackBuffer& data) const;

	/*!
	  \brief Unpack data from a char vector into this class

	  The vector data contains all information to rebuild the
	  exact copy of an instance of a class on a different processor.
	  The first entry in data has to be an integer which is the unique
	  parobject id defined at the top of this file and delivered by
	  UniqueParObjectId().

	  \param data (in) : vector storing all data to be unpacked into this
	  instance.
	 */
	virtual void Unpack(const vector<char>& data);

	//@}

	/// material type
	virtual INPAR::MAT::MaterialType MaterialType() const {
	    return INPAR::MAT::m_charmm;
	}

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new CHARMM( *this ) ); }

	/// Initialize history variables
	void Setup(
		DRT::Container& data_ ///< Container for History
	);


	/*!
	 \brief The main routine to evaluate stress and constitutive tensor

	 This is main routine to evaluate the stress and constitutive tensor.
	 The calling of CHARMm is initiated here as well. All preparations for
	 the characteristic strains and the computation of the material
	 properties is done here as well.

	 \param glstrain (in) : Green-Lagrange strains
	 \param cmat (out) : constitutive tensor in Voigt notation
	 \param stress (out) : Second Piola-Kirchoff stress tensor
	 \param ele_ID (in) : unique element id
	 \param gp (in) : Gaussian point
	 \param data_ (in/out) : Container for history data
	 \param time (in) : current time
	 \param xrefe (in) : reference coordinates
	 \param xcurr (in) : current coordinates
	 */
	void Evaluate(
		const LINALG::Matrix<NUM_STRESS_3D, 1 > * glstrain,
		LINALG::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>* cmat,
		LINALG::Matrix<NUM_STRESS_3D, 1 > * stress,
		const int ele_ID,
		const int gp,
		DRT::Container& data_,
		const double time,
		const LINALG::SerialDenseMatrix& xrefe,
		const LINALG::SerialDenseMatrix& xcurr
	);

	/*!
	 \brief Compute strains in the charateristic directions

	 This routine computes the principal strains (lambda) and the unit
	 vector of the direction in the charateristic directions. Depending on
	 the wish of the user, it might be done on the origin of the element.

	 \param origin (in) : Strains at the origin (true/false)
	 \param xrefe (in) : reference coordinates for F_0
	 \param xcurr (in) : current coordinates for F_0
	 \param C (in/out) : Right Cauchy Green tensor (in); Vector (out)
	 \param lambda (out) : principal strains
	 */
	void EvalStrain(
		const bool& origin,
		const LINALG::SerialDenseMatrix& xrefe,
		const LINALG::SerialDenseMatrix& xcurr,
		LINALG::Matrix < 3, 3 > & C,
		LINALG::SerialDenseVector& lambda
	);

	/*!
	 \brief Compute second time derivative of F in diretcion of the first characteristic direction

	 Here the acceleration in the first characteritic direction is computed.
	 This force may be later used for pulling on the atom or as force in various
	 constraints in MD. This addes a dynamic information to the coupling.

	 \param STARTD (in) : length of the last timestep
	 \param ENDD (in) : length at the current timestep
	 \param time_STARTD (in) : time at last timestep
	 \param time_ENDD (in) : time at current timestep
	 \param v_his (in) : v at last timestep
	 \param atomic_mass (in) : Atomic mass of moving part
	 \param Facc_scale (in) : Scale your force down to pN
	 \param FCD_v (out) : current v
	 \param a (out) : acceleration of the moving molecule
	 \param Force (out) : Force on the moving part
	 */
	void EvalAccForce(
		const double& STARTD,
		const double& ENDD,
		const double& time_STARTD,
		const double& time_ENDD,
		const double& v_his,
		const double& atomic_mass,
		const double& Facc_scale,
		double& FCD_v,
		double& a,
		double& Force
	);

	/*!
	 \brief Call external Charmm binary

	 Here CHARMm gets actually called. The business of MPI for CHARMm is
	 handeled here as well. In- and Output files are defined.

	 \param CHARMmPar (in) : Map with all necessary CHARMm parameter
	 \param charmm_result (out) : The CHARMm result vector
	*/
	void CHARMmfileapi(
		map<string, double>& CHARMmPar,
		LINALG::SerialDenseVector& charmm_result);

        /*!
	 \brief Read CHARMm Results

	 After CHARMm finished the results must be read in and the status of
         the run must be checked. All this business is done here.

	 \param CHARMmPar (in) : Map with all necessary CHARMm parameter
	 \param charmm_result (out) : The CHARMm result vector
	*/
	void Readresults(
		map<string, double>& CHARMmPar,
		LINALG::SerialDenseVector& charmm_result);

        /*!
	 \brief The actual routine wich opens and reads the results and status

	 Reader for status and result files written by CHARMm

	 \param file (in) : output string stream of the file
	 \param charmm_result (out) : content of the file as map
	*/
	void Reader(
		const ostringstream& file,
		map<string, double>& content);


	/*!
	 \brief Hardcoded CHARMm results for testing

	 This functions allows for testing purpose and external run dynamics
	 runs the storing of CHARMm results.

	 \param STARTD (in) : FCD beginning
	 \param ENDD (in) : FCD end
	 \param charmm_result (out) : CHARMm result vector
	 */
	void CHARMmfakeapi(
		const double STARTD,
		const double ENDD,
		LINALG::SerialDenseVector& charmm_result
	);

	/// \b ORIGIN: boolean (1/0) to decide if F is computed at the element origin
	inline bool Origin() const {
	    if (params_->origin_ == 1) return true;
	    return false;
	}
	/// \b FCL: first characteristic (fc) length (double) [AKMA]
	inline double FCL() const { return params_->fcl_; }
	/// \b FCD_Type: type of the fc direction (string) (principal/vector/none)
	inline string FCDType() { return (*params_->fcd_type_); }
	/// \b FCD: direction of fc direction if fcd_type_ is vector (unit vector)
	inline vector<double> FCD() { return (*params_->fcd_); }
	/// \b FCD_Space: directional space of fc (unit vector)
	inline vector<double> FCDS() { return (*params_->fcds_); }

	/// second characteristic (sc) length [AKMA]
	inline double SCL() const { return params_->scl_; }
	/// \b SCD_Type: type of the sc direction (string) (principal/vector/none)
	inline string SCDType() { return (*params_->scd_type_); }
	/// \b SCD: direction of sc direction if scd_type_ is vector (unit vector)
	inline vector<double> SCD() { return (*params_->scd_); }
	/// \b SCD_Space: directional space of sc (unit vector)
	inline vector<double> SCDS() { return (*params_->scds_); }

	/// \b FCD_Acceleration: boolean (1/0) if F'' in fcd is computed and as force given to MD
	inline bool FCDAcc() const {
		if (params_->fcdacc_ == 1) return true;
		return false;
	}
	/// \b AtomicMass: Atomic mass of the moving part of the protein [AKMA]
	inline double AtomicMass() const { return params_->atomicmass_; }
	/// \b Facc_Scale: Scale factor from whatever force unit used in FE to pN
	inline double FaccScale() const { return params_->facc_scale_; }
	/// \b Time_AKMA: Scale factor from FE delta t time to AKMA time
	inline double TimeAKMA() const { return params_->timeakma_; }
	/// \b Time_Scale: Linear scaling between FE and MD delta t time
	inline double TimeScale() const { return params_->timescale_; }

	// \b HARD: Use hard coded results e.g. for testing or external run MD
	inline bool HARD() const {
		if (params_->hard_ == 1) return true;
		return false;
	}
	/// \b c_Scale: Scale factor for material parameter c1 (Neo-Hookean) depending on used FE units
	inline double CScale() const { return params_->c_scale_; }
	/// \b PATH: Path to the folder structre for CHARMm problem type
	inline string Path() { return (*params_->path_); }
	/// \b USE_OLD_RESULTS: boolean (1/0) if previously computed results from CHARMm should be used
        inline bool Use_old_Results() const {
	    if (params_->use_old_results_ == 1) return true;
	    return false;
	}
	/// \b SERPAR: serial or parallel execution of CHARMm
	inline string Serpar() { return (*params_->serpar_); }
	/// \b CHARMM: CHARMm executable
	inline string CHARMMEXE() { return (*params_->charmm_); }
	/// \b INPUT: CHARMm input file
	inline string INPUT() { return (*params_->input_); }
	/// \b NUE: Poisson ratio
	inline double NUE() const { return params_->nue_; }
	/// \b DENS: density
	virtual double Density() const {
	    return params_->density_;
	}

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }


    private:

	/// my material parameters
	MAT::PAR::CHARMM* params_;

    };
}

#endif // CHARMM_H_

