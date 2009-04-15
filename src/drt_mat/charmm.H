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
#ifdef CCADISCRET

#include "material.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"


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

	    /// mass density
	    const double density_;

	    //@}

	}; // class CHARMM

    } // namespace PAR

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
	    return ParObject_CHARMM;
	}

	/*!
	  \brief Pack this class so it can be communicated

	  Resizes the vector data and stores all information of a class in it.
	  The first information to be stored in data has to be the
	  unique parobject id delivered by UniqueParObjectId() which will then
	  identify the exact class on the receiving processor.

	  \param data (in/out): char vector to store class information
	 */
	virtual void Pack(vector<char>& data) const;

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
	 \brief Call external Charmm binary

	 Here CHARMm gets actually called. The business of MPI for CHARMm is
	 handeled here as well. In- and Output files are defined.

	 \param FCD_STARTD (in) : First characteristic strain beginning
	 \param FCD_ENDD (in) : First characteristic strain end
	 \param FCD_direction (in) : First charteristic directional vector
	 \param SCD_STARTD (in) : Second characteristic strain beginning
	 \param SCD_ENDD (in) : Second characteristic strain end
	 \param SCD_direction (in) : Second charteristic directional vector
	 \param charmm_result (out) : The CHARMm result vector
	*/
	void CHARMmfileapi(
		const double FCD_STARTD,
		const double FCD_ENDD,
		const LINALG::SerialDenseVector FCD_direction,
		const double SCD_STARTD,
		const double SCD_ENDD,
		const LINALG::SerialDenseVector SCD_direction,
		LINALG::SerialDenseVector& charmm_result);

	/*!
	 \brief Function to read results from the minimized CHARMm output

	 Here the actual output from CHARMm is read and processed. The success
	 of the run is checked as well. All read results are from a cold
	 (minimized only) run from CHARMm.

	 Example CHARMm input file: 1dzi_fem_min.inp

	 \param outputfile (in) : Name of the CHARMm outputfile
	 \param energyfile (in) : Name of the energy file
	 \param volumefile (in) : Name of the volume file
	 \param debug (in) : printout level
	 \param charmm_result (out) : CHARMm result vector
	 */
	void Readcoldresults(
		const ostringstream& outputfile,
		const ostringstream& energyfile,
		const ostringstream& volumefile,
		const int debug,
		LINALG::SerialDenseVector& charmm_result
	);

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

	/// density
	inline double Density() const {
	    return params_->density_;
	}


    private:

	/// my material parameters
	MAT::PAR::CHARMM* params_;

    };
}

#endif // CCADISCRET
#endif // CHARMM_H_

