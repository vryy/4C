/*----------------------------------------------------------------------*/
/*!
\file myocard.H
\brief

<pre>
Maintainer: Lasse Jagschies
            lasse.jagschies@tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
</pre>
*/

/*----------------------------------------------------------------------*
 |  definitions                                              dano 09/09 |
 *----------------------------------------------------------------------*/
#ifndef MYOCARD_H
#define MYOCARD_H

/*----------------------------------------------------------------------*
 |  headers                                                  dano 09/09 |
 *----------------------------------------------------------------------*/
#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>

/*----------------------------------------------------------------------*
 |                                                           dano 09/09 |
 *----------------------------------------------------------------------*/
namespace MAT
{

    namespace PAR
	{

  /*----------------------------------------------------------------------*/
	/// parameters for scalar transport material
	class Myocard : public Parameter
	    {
	    public:

	    /// standard constructor
	    Myocard( Teuchos::RCP<MAT::PAR::Material> matdata );

	    /// destructor
	    virtual ~Myocard() { ; }

	    /// @name material parameters
	    //@{

	    /// Diffusivity
	    const double maindirdiffusivity;
	    const double offdirdiffusivity;

            /// Pertubation for calculation of derivative
	    const double dt_deriv;

	    /// Model parameter
	    const double u_o; //
	    const double u_u;
	    const double Theta_v;
	    const double Theta_w;
	    const double Theta_vm;
	    const double Theta_o;
	    const double Tau_v1m;
	    const double Tau_v2m;
	    const double Tau_vp;
	    const double Tau_w1m;
	    const double Tau_w2m;
	    const double k_wm;
	    const double u_wm;
	    const double Tau_wp;
	    const double Tau_fi;
	    const double Tau_o1;
	    const double Tau_o2;
	    const double Tau_so1;
	    const double Tau_so2;
	    const double k_so;
	    const double u_so;
	    const double Tau_s1;
	    const double Tau_s2;
	    const double k_s;
	    const double u_s;
	    const double Tau_si;
	    const double Tau_winf;
	    const double w_infs;

	    //@}

	    /// create material instance of matching type with my parameters
	    virtual Teuchos::RCP<MAT::Material> CreateMaterial();

	    }; // class myocard

	} // namespace PAR

    class MyocardType : public DRT::ParObjectType
	{
	public:

	std::string Name() const { return "MyocardType"; }

	static MyocardType & Instance() { return instance_; };

	virtual DRT::ParObject* Create( const std::vector<char> & data );

	private:

	static MyocardType instance_;
	};

    /*----------------------------------------------------------------------*/
    /// Myocard material according to [1], [2]
    ///
    /// This is a reaction-diffusion law of anisotropic, instationary electric conductivity in cardiac muscle tissue
    ///
    /// <h3>References</h3>
    /// <ul>
    /// <li> [1] A Bueno-Orovio et. al., "Minimal model for human ventricular action potentials in tissue", Journal of Theoretical Biology 253 (2008) 544-560
    /// <li> [2] KHWJ ten Tusscher et. al., "A model for human ventricular tissue", American Journal of Physiology Heart and Circulatory Physiology 286 (2004) H1573-H1589
    /// <li> [3] S Rush and H Larsen, "A practical algorithm for solving dynamic membrane equations", IEEE Trans on Biomedical Engineering 25 (1978) 389-392
    /// </ul>
    ///
    /// \author ljag
    /// \date 06/12
    class Myocard : public Material
	{
	public:

	/// construct empty material object
	Myocard();

	/// constructor with given material parameters
	Myocard(MAT::PAR::Myocard* params);

	/// @name Packing and Unpacking
	//@{

	/// Return unique ParObject id
	///
	///  every class implementing ParObject needs a unique id defined at the
	///  top of drt_parobject.H (this file) and should return it in this method.
	virtual int UniqueParObjectId() const { return MyocardType::Instance().UniqueParObjectId(); }

	/// Pack this class so it can be communicated
	///
	/// Resizes the vector data and stores all information of a class in it.
	/// The first information to be stored in data has to be the
	/// unique parobject id delivered by UniqueParObjectId() which will then
	/// identify the exact class on the receiving processor.
	virtual void Pack(DRT::PackBuffer& data) const;  ///< (in/out): char vector to store class information

	/// \brief Unpack data from a char vector into this class
	///
	/// The vector data contains all information to rebuild the
	/// exact copy of an instance of a class on a different processor.
	/// The first entry in data has to be an integer which is the unique
	/// parobject id defined at the top of this file and delivered by
	/// UniqueParObjectId().
	///
	virtual void Unpack(const std::vector<char>& data);  ///< vector storing all data to be unpacked into this

	//@}

	/// material type
	virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_myocard; }

	/// return copy of this material object
	virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new Myocard( *this ) ); }


	/// material call from DRT::ELEMENTS::Transport::ReadElement function
        /// to setup conductivity tensor for each element
	void Setup(DRT::INPUT::LineDefinition* linedef);
  void Setup(const std::vector<double> &fiber1);

	/// diffusivity
	void ComputeDiffusivity(LINALG::Matrix<1,1>& diffus3) const;
	void ComputeDiffusivity(LINALG::Matrix<2,2>& diffus3) const;
	void ComputeDiffusivity(LINALG::Matrix<3,3>& diffus3) const;
	void SetupDiffusionTensor(const std::vector<double> &fiber1);

	/// compute reaction coefficient
	double ComputeReactionCoeff(const double phi, const double dt) const;

	/// compute reaction coefficient derivative
	double ComputeReactionCoeffDeriv(const double csnp, const double dt) const;

	/// compute Heaviside step function
	double GatingFunction(const double Gate1, const double Gate2, const double p, const double var, const double thresh) const;

	/// time update for this material
	void Update(const double phi, const double dt);

	/// @name Access material constants
	//@{

	//@}

	/// Return quick accessible material parameter data
	virtual MAT::PAR::Parameter* Parameter() const { return params_; }

	private:

	/// my material parameters
	MAT::PAR::Myocard* params_;

	/// conductivity tensor
	LINALG::Matrix<3,3> difftensor_;

	/// initial gating variables
	double v0_; /// fast inward current
	double w0_; /// slow inward current
	double s0_; /// slow inward current


	}; // Myocard

} // namespace MAT

#endif
