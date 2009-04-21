/*!----------------------------------------------------------------------
\file protein.H
\brief CHARMm Interface to compute the mechanical properties of protein

<pre>
Maintainer: Robert Metzke
	    metzke@lnm.mw.tum.de
	    http://www.lnm.mw.tum.de
	    089 - 289-15244
</pre>
 *----------------------------------------------------------------------*/
#ifndef PROTEIN_H_
#define PROTEIN_H_
#ifdef CCADISCRET

#include "material.H"

namespace MAT {

    namespace PAR {

	/*----------------------------------------------------------------------*/
	/// material parameters for Integrin

	class PROTEIN : public Parameter {
	public:

	    /// standard constructor
	    PROTEIN(
		    Teuchos::RCP<MAT::PAR::Material> matdata
		    );

	    /// destructor

	    virtual ~PROTEIN() {
		;
	    }

	    /// @name material parameters
	    //@{

	    /// mass density
	    const double density_;

	    //@}

	}; // class PROTEIN

    } // namespace PAR

    /*----------------------------------------------------------------------*/

    /*!
    \brief CHARMm Interface to compute the mechanical properties of protein

    */

    class PROTEIN : public Material {
    public:

	/// construct empty material object
	PROTEIN();

	/// construct the material object given material parameters
	explicit PROTEIN(MAT::PAR::PROTEIN* params);

	//! @name Packing and Unpacking

	/*!
	  \brief Return unique ParObject id

	  every class implementing ParObject needs a unique id defined at the
	  top of drt_parobject.H (this file) and should return it in this method.
	 */
	virtual int UniqueParObjectId() const {
	    return ParObject_Protein;
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
	    return INPAR::MAT::m_protein;
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

	/// density
	inline double Density() const {
	    return params_->density_;
	}


    private:

	/// my material parameters
	MAT::PAR::PROTEIN* params_;

    };
}

#endif // CCADISCRET
#endif // PROTEIN_H_
