/*----------------------------------------------------------------------*/
/*!
\file viscoelasthyper.H

<pre>
Maintainer: Anna Birzle
	    anna.birzle@tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef VISCOELASTHYPER_H
#define VISCOELASTHYPER_H



#include "matpar_parameter.H"
#include "so3_material.H"
#include "../drt_lib/drt_parobjectfactory.H"

// forward declaration due to avoid header definition
namespace MAT
{
  namespace ELASTIC
  {
    class Summand;
  }

  // forward declaration
  class ViscoElastHyper;

  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// Collection of viscohyperelastic materials
    ///
    /// Storage map of hyperelastic summands.


    class ViscoElastHyper : public Parameter
    {

      friend class MAT::ViscoElastHyper;

    public:

      /// standard constructor
      ///
      /// This constructor recursively calls the constructors of the
      /// parameter sets of the hyperelastic summands.
      ViscoElastHyper(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~ViscoElastHyper() { ; }

      /// @name material parameters
      //@{

//       /// provide access to material/summand by its ID
//       Teuchos::RCP<const MAT::ELASTIC::Summand> MaterialById(
//         const int id  ///< ID to look for in collection of summands
//         ) const;

      /// length of material list
      const int nummat_;

      /// the list of material IDs
      const std::vector<int>* matids_;

      /// material mass density
      const double density_;

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      //@}

    };  // class ViscoElastHyper

  }  // namespace PAR

class ViscoElastHyperType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "ViscoElastHyperType"; }

  static ViscoElastHyperType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static ViscoElastHyperType instance_;
};


  /*----------------------------------------------------------------------*/
  /// Collection of hyperelastic materials
  ///
  /// This collection offers to additively compose a stress response
  /// based on summands defined separately.  This is possible, because
  /// we deal with hyperelastic materials, which are composed
  /// of (Helmholtz free energy density) potentials.  Effectively, we want
  ///\f[
  ///  \Psi(\boldsymbol{C}) = \sum_i \Psi_i(\boldsymbol{C})
  ///\f]
  /// in which the individual \f$\Psi_i\f$ is implemented as #MAT::ELASTIC::Summand.
  ///
  /// Quite often the right Cauchy-Green 2-tensor \f$\boldsymbol{C}\f$
  /// is replaced by its various invariant forms as argument.
  ///
  /// The task of ElastHyper is the evaluation of the
  /// potential energies and their derivatives to obtain the actual
  /// stress response and the elasticity tensor. The storage is located
  /// at the associated member #params_.
  ///
  /// <h3>References</h3>
  /// <ul>
  /// <li> [1] GA Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
  /// </ul>
  ///
  /// \author rausch,tk,bborn
  /// \date 05/09
  class Material;

  class ViscoElastHyper : public So3Material
  {
  public:

    /// construct empty material object
    ViscoElastHyper();

    /// construct the material object given material parameters
    explicit ViscoElastHyper(MAT::PAR::ViscoElastHyper* params);

    ///@name Packing and Unpacking
    //@{

    /// \brief Return unique ParObject id
    ///
    /// every class implementing ParObject needs a unique id defined at the
    /// top of drt_parobject.H (this file) and should return it in this method.
    virtual int UniqueParObjectId() const { return ViscoElastHyperType::Instance().UniqueParObjectId(); }

    /// \brief Pack this class so it can be communicated
    ///
    /// Resizes the vector data and stores all information of a class in it.
    /// The first information to be stored in data has to be the
    /// unique parobject id delivered by UniqueParObjectId() which will then
    /// identify the exact class on the receiving processor.
    ///
    /// \param data (in/out): char vector to store class information
    virtual void Pack(DRT::PackBuffer& data) const;

    /// \brief Unpack data from a char vector into this class
    ///
    /// The vector data contains all information to rebuild the
    /// exact copy of an instance of a class on a different processor.
    /// The first entry in data has to be an integer which is the unique
    /// parobject id defined at the top of this file and delivered by
    /// UniqueParObjectId().
    ///
    /// \param data (in) : vector storing all data to be unpacked into this
    ///                    instance.
    virtual void Unpack(const std::vector<char>& data);

    //@}

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_viscoelasthyper; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new ViscoElastHyper( *this ) ); }

    /// Check if history variables are already initialized
    virtual bool Initialized() const
    {
    	return isinitvis_&&(histmodrcgcurr_!=Teuchos::null);
    }

    /// number of materials
    virtual int NumMat() const { return params_->nummat_; }

    /// deliver material ID of index i'th potential summand in collection
    virtual int MatID(
      const unsigned index  ///< index i
      ) const;

    /// provide access to material by its ID
//     virtual Teuchos::RCP<const MAT::ELASTIC::Summand> MaterialById(const int id) const { return params_->MaterialById(id); }

    /// material mass density
    virtual double Density() const { return params_->density_; }

    /// a shear modulus equivalent
    virtual double ShearMod() const;

    /// evaluate strain energy function
    virtual void StrainEnergy(const LINALG::Matrix<6,1>& glstrain, ///< Green-Lagrange strain
                              double& psi);                        ///< Strain energy function

    /// hyperelastic stress response plus elasticity tensor
    virtual void Evaluate(const LINALG::Matrix<3,3>* defgrd,  ///< Deformation gradient
                          const LINALG::Matrix<6,1>* glstrain,///< Green-Lagrange strain
                          Teuchos::ParameterList& params,     ///< Container for additional information
                          LINALG::Matrix<6,1>* stress,        ///< 2nd Piola-Kirchhoff stresses
                          LINALG::Matrix<6,6>* cmat);         ///< Constitutive matrix

    virtual void Setup(int numgp, DRT::INPUT::LineDefinition* linedef);

    virtual void ResetAll(const int numgp);
    virtual void SetupTimeVariables();

    virtual void Update();

    /// setup patient-specific AAA stuff
    virtual void SetupAAA(Teuchos::ParameterList& params);

    /// return if isotropic not splitted formulation
    virtual bool IsotropicPrincipal() const { return isoprinc_;}

    /// return if isotropic and splitted formulation
    virtual bool IsotropicModified() const { return isomod_;}

    /// return if anisotropic not splitted formulation
    virtual bool AnisotropicPrincipal() const { return anisoprinc_;}

    /// return if anisotropic and splitted formulation
    virtual bool AnisotropicModified() const { return anisomod_;}

    /// return if isotropic, splitted and viscos formulation
    virtual bool IsotropicModifiedVisco() const { return isomodvisco_;}

    /// get fiber vectors
    virtual void GetFiberVecs(std::vector<LINALG::Matrix<3,1> >& fibervecs);

    /// evaluate fiber directions from locsys and alignment angle, pull back
    virtual void EvaluateFiberVecs
    (
      const double newgamma,   ///< new angle
      const LINALG::Matrix<3,3>& locsys,   ///< local coordinate system
      const LINALG::Matrix<3,3>& defgrd   ///< deformation gradient
    );

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    /// Return names of visualization data
    virtual void VisNames(std::map<std::string,int>& names);

    /// Return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID);

  protected:
    /// calculates the kinematic quantities and tensors used afterwards
    virtual void EvaluateKinQuant(
      const LINALG::Matrix<6,1>& glstrain,
      LINALG::Matrix<6,1>& id2,
      LINALG::Matrix<6,1>& scg,
      LINALG::Matrix<6,1>& rcg,
      LINALG::Matrix<6,1>& icg,
      LINALG::Matrix<6,6>& id4,
      LINALG::Matrix<6,6>& id4sharp,
      LINALG::Matrix<3,1>& prinv,
      LINALG::Matrix<3,1>& modinv);

    /// calculates the kinematic quantities and tensors used afterwards for viscos part
    virtual void EvaluateKinQuantVis(
      LINALG::Matrix<6,1> rcg,
      LINALG::Matrix<6,1>& modrcg,
      LINALG::Matrix<6,1> icg,
      Teuchos::ParameterList& params,
      LINALG::Matrix<3,1> prinv,
      LINALG::Matrix<6,1>& modrcgrate,
      LINALG::Matrix<7,1>& modrateinv
      );

    /// calculates the factors associated to the hyperelastic laws
    virtual void EvaluateGammaDelta(
        const LINALG::Matrix<6,1> rcg,
        LINALG::Matrix<3,1> prinv,
        LINALG::Matrix<3,1> modinv,
        LINALG::Matrix<3,1>& gamma,
        LINALG::Matrix<8,1>& delta,
        LINALG::Matrix<3,1>& modgamma,
        LINALG::Matrix<5,1>& moddelta
        );

    /// calculates the factors associated to the viscos laws
    virtual void EvaluateMyXi(
        LINALG::Matrix<3,1> modinv,
        LINALG::Matrix<8,1>& modmy,
        LINALG::Matrix<33,1>& modxi,
        LINALG::Matrix<7,1>& modrateinv,
        Teuchos::ParameterList& params
        );

    /// calculates the isotropic stress and elasticity tensor for coupled configuration
    virtual void EvaluateIsotropicPrinc(
        LINALG::Matrix<6,1>& stressisoprinc,
        LINALG::Matrix<6,6>& cmatisoprinc,
        LINALG::Matrix<6,1> scg,
        LINALG::Matrix<6,1> id2,
        LINALG::Matrix<6,1> icg,
        LINALG::Matrix<6,6> id4sharp,
        LINALG::Matrix<3,1> gamma,
        LINALG::Matrix<8,1> delta
        );

    /// calculates the isotropic stress and elasticity tensor for decoupled configuration
    virtual void EvaluateIsotropicMod(
        LINALG::Matrix<6,1>& stressisomodiso,
        LINALG::Matrix<6,1>& stressisomodvol,
        LINALG::Matrix<6,6>& cmatisomodiso,
        LINALG::Matrix<6,6>& cmatisomodvol,
        LINALG::Matrix<6,1> rcg,
        LINALG::Matrix<6,1> id2,
        LINALG::Matrix<6,1> icg,
        LINALG::Matrix<6,6> id4,
        LINALG::Matrix<6,6> id4sharp,
        LINALG::Matrix<3,1> modinv,
        LINALG::Matrix<3,1> prinv,
        LINALG::Matrix<3,1> modgamma,
        LINALG::Matrix<5,1> moddelta
        );

    /// calculates the isotropic stress and elasticity tensor for viscos decoupled configuration
    virtual void EvaluateIsoModVisco(
        LINALG::Matrix<6,1>& stressisomodisovisco,
        LINALG::Matrix<6,1>& stressisomodvolvisco,
        LINALG::Matrix<6,6>& cmatisomodisovisco,
        LINALG::Matrix<6,6>& cmatisomodvolvisco,
        LINALG::Matrix<3,1> prinv,
        LINALG::Matrix<3,1> modinv,
        LINALG::Matrix<8,1> modmy,
        LINALG::Matrix<33,1> modxi,
        LINALG::Matrix<6,1> rcg,
        LINALG::Matrix<6,1> id2,
        LINALG::Matrix<6,1> icg,
        LINALG::Matrix<6,6> id4,
        LINALG::Matrix<6,1> modrcgrate
        );

    /// calculates the anisotropic stress and elasticity tensor for coupled configuration
    virtual void EvaluateAnisotropicPrinc(
        LINALG::Matrix<6,1>& stressanisoprinc,
        LINALG::Matrix<6,6>& cmatanisoprinc,
        LINALG::Matrix<6,1> rcg,
        Teuchos::ParameterList& params
        );

    /// calculates the anisotropic stress and elasticity tensor for decoupled configuration
    virtual void EvaluateAnisotropicMod(
        LINALG::Matrix<6,1>& stressanisoprinc,
        LINALG::Matrix<6,6>& cmatanisoprinc,
        LINALG::Matrix<6,1> rcg,
        LINALG::Matrix<6,1> icg,
        LINALG::Matrix<3,1> prinv
        );

    /// calculate principal invariants
    virtual void InvariantsPrincipal(
      LINALG::Matrix<3,1>& prinv,  ///< principal invariants
      const LINALG::Matrix<6,1>& rcg  ///< symmetric Cartesian 2-tensor in strain-like 6-Voigt notation
      );

    /// calculate modified invariants
    virtual void InvariantsModified(
      LINALG::Matrix<3,1>& modinv,  ///< modified invariants
      const LINALG::Matrix<3,1>& prinv  ///< principal invariants
      );

    /// calculate principal stretches \f$\lambda_\alpha\f$ and material principal directions
    void StretchesPrincipal(
     LINALG::Matrix<3,1>& prstr,  ///< principal stretches
     LINALG::Matrix<3,3>& prdir,  ///< principal directions, stored col-by-col
     const LINALG::Matrix<6,1>& rcg);  ///< symmetric Cartesian 2-tensor in strain-like 6-Voigt notation


    /// calculate modified principal stretches \f$\bar{\lambda}_\alpha\f$
    void StretchesModified(
     LINALG::Matrix<3,1>& modstr,  ///< modified principal stretches
     const LINALG::Matrix<3,1>& prstr  ///< principal stretches
     );

    /// determine if at least one summand uses principal stretches
    bool HaveCoefficientsStretchesPrincipal();

    /// determine if at least one summand uses modified principal stretches
    bool HaveCoefficientsStretchesModified();

    /// Determine PK2 stress response and material elasticity tensor
    /// due to energy densities (MAT::ELASTIC::Summand) described
    /// in (modified) principal stretches.
    ///
    /// The stress response is achieved by collecting the coefficients
    /// \f$\gamma_\alpha\f$ and \f$\delta_{\alpha\beta}\f$ due to
    /// MAT::ELASTIC::Summand::AddCoefficientsStretchesPrincipal()
    /// (and/or \f$\bar{\gamma}_\alpha\f$ and \f$\bar{\delta}_{\alpha\beta}\f$
    /// due to MAT::ELASTIC::Summand::AddCoefficientsStretchesModified()).
    /// The collected coefficients build the principal 2nd Piola--Kirchhoff
    /// stresses which are transformed into ordinary Cartesian co-ordinates
    /// applying the principal directions. Similarly, the elasticity
    /// 4-tensor is obtained.
    ///
    /// Please note, unlike suggested by Holzapfel, p 263-264, the modified
    /// coefficients are transformed to unmodified coefficients and than added
    /// onto the respective quantities. The transformation goes along the following
    /// lines. The first derivatives or \f$\gamma_\alpha\f$ coefficients are related
    /// to the modified coefficients \f$\bar{\gamma}_\alpha\f$ via the chain rule
    ///\f[
    ///   \gamma_\alpha
    ///   = \frac{\partial \Psi(\boldsymbol{\lambda})}{\partial\lambda_\alpha}
    ///   = \frac{\partial \Psi(\bar{\boldsymbol{\lambda}})}{\partial \bar{\lambda}_\eta}
    ///   \, \frac{\partial \bar{\lambda}_\eta}{\partial\lambda_\alpha}
    ///   = \bar{\gamma}_\eta \, \frac{\partial \bar{\lambda}_\eta}{\partial\lambda_\alpha}
    ///\f]
    /// utilising Holzapfel Eq (6.142):
    ///\f[
    ///   \frac{\partial\bar{\lambda}_\alpha}{\partial\lambda_\beta}
    ///   = J^{-1/3} \Big( 1_{\alpha\beta} - \frac{1}{3} \lambda_\alpha \lambda_\beta^{-1} \Big)
    ///\f]
    /// in which Kronecker's delta is denoted \f$1_{\alpha\beta}\f$.
    /// Likewise (and once without Holzapfel), the second derivatives \f$\delta_{\alpha\beta}\f$ can
    /// be recovered directly by knowledge of the modified coefficients
    /// \f$\bar{\gamma}_\alpha\f$ and \f$\bar{\delta}_{\alpha\beta}\f$, i.e.
    ///\f[
    ///   \delta_{\alpha\beta}
    ///   = \frac{\partial^2\Psi(\boldsymbol{\lambda})}{\partial\lambda_\alpha \partial\lambda_\beta}
    ///   = \frac{\partial^2\Psi(\bar{\boldsymbol{\lambda}})}{\partial\bar{\lambda}_\eta\partial\bar{\lambda}_\epsilon}
    ///     \, \frac{\partial \bar{\lambda}_\eta}{\partial\lambda_\alpha}
    ///     \, \frac{\partial \bar{\lambda}_\epsilon}{\partial\lambda_\beta}
    ///   + \frac{\partial \Psi(\bar{\boldsymbol{\lambda}})}{\partial \bar{\lambda}_\eta}
    ///     \, \frac{\partial^2 \bar{\lambda}_\eta}{\partial\lambda_\alpha\partial\lambda_\beta}
    ///   = \bar{\delta}_{\eta\epsilon}
    ///     \, \frac{\partial \bar{\lambda}_\eta}{\partial\lambda_\alpha}
    ///     \, \frac{\partial \bar{\lambda}_\epsilon}{\partial\lambda_\beta}
    ///   + \bar{\gamma}_\eta
    ///     \, \frac{\partial^2 \bar{\lambda}_\eta}{\partial\lambda_\alpha\partial\lambda_\beta}
    ///\f]
    /// with
    ///\f[
    ///   \frac{\partial^2 \bar{\lambda}_\alpha}{\partial\lambda_\beta\partial\lambda_\eta}
    ///   = -\frac{1}{3} J^{-1/3} \Big(
    ///     1_{\alpha\beta} \lambda_\eta^{-1}
    ///     + 1_{\alpha\eta}\lambda_\beta^{-1}
    ///     - \big( 1_{\beta\eta} + \frac{1}{3} \big) \lambda_\alpha \lambda_\beta^{-1} \lambda_\eta^{-1}
    ///   \Big)
    ///\f]
    ///
    /// <h3>References</h3>
    /// See Holzapfel, p 245-246, p 257-259, p 263-264
    ///
    /// \author bborn
    /// \date 08/09

    virtual void ResponseStretches(
    LINALG::Matrix<6,6>& cmat,  ///< material elasticity matrix
    LINALG::Matrix<6,1>& stress,  ///< 2nd Piola-Kirchhoff stress vector
    const LINALG::Matrix<6,1> rcg,  ///< Right Cauchy-Green tensor
    const bool& havecoeffstrpr,  ///< true, if at least one summand is based on pr.stretches
    const bool& havecoeffstrmod  ///< true, if at least one summand is based on mod.pr.stretches
    );

    /// Flags to specify the formulations
    bool isoprinc_;
    bool isomod_;
    bool anisoprinc_ ;
    bool anisomod_ ;
    bool isomodvisco_;

    /// @name Voigt index arrays
    //@{
    static const int VOIGT6ROW_[];  ///< 6-Voigt row index of corresponding 2-tensor
    static const int VOIGT6COL_[];  ///< 6-Voigt column index of corresponding 2-tensor
    static const int VOIGT3X3SYM_[]; ///< go from symmtric 2-tensor index pair to 6-Voigt index
                                     ///< by [3*i+j] for any i,j=0,1,2
    //@}

    /// my material parameters
    MAT::PAR::ViscoElastHyper* params_;

    /// map to materials/potential summands
    std::vector<Teuchos::RCP<MAT::ELASTIC::Summand> > potsum_;

  private:

    /// visco history Cauchy-Green-Tensor
    Teuchos::RCP<std::vector<LINALG::Matrix<6,1> > > histmodrcgcurr_; ///< current Cauchy-Green-Tensor
    Teuchos::RCP<std::vector<LINALG::Matrix<6,1> > > histmodrcglast_; ///< Cauchy-Green-Tensor of last converged state

    bool isinitvis_;  ///< indicates if #Initialized routine has been called
  };

}

#endif

