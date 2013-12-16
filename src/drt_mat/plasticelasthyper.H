/*----------------------------------------------------------------------*/
/*!
\file plasticelasthyper.H

<pre>
Maintainer: Alexander Seitz
            seitz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef PLASTICELASTHYPER_H
#define PLASTICELASTHYPER_H



#include "matpar_parameter.H"
#include "so3_material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "elasthyper.H"

// forward declaration due to avoid header definition
namespace MAT
{
  namespace ELASTIC
  {
    class Summand;
  }

  // forward declaration
  class PlasticElastHyper;

  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// Collection of hyperelastic materials
    ///
    /// Storage map of hyperelastic summands.
    class PlasticElastHyper : public Parameter
    {

      friend class MAT::PlasticElastHyper;

    public:

      /// standard constructor
      ///
      /// This constructor recursively calls the constructors of the
      /// parameter sets of the hyperelastic summands.
      PlasticElastHyper(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~PlasticElastHyper() { ; }

      /// @name material parameters
      //@{

       /// provide access to material/summand by its ID
       Teuchos::RCP<const MAT::ELASTIC::Summand> MaterialById(
         const int id  ///< ID to look for in collection of summands
         ) const;


      /// length of material list
      const int nummat_;

      /// the list of material IDs
      const std::vector<int>* matids_;

      /// material mass density
      const double density_;

      // plastic material parameters ************************************
      // initial yield stress
      const double inityield_;
      // linear isotropic hardening parameter
      const double isohard_;
      // exponent for nonlinear isotropic hardening
      const double expisohard_;
      // saturation yield stress for nonlinear isotropic hardening
      const double infyield_;
      // linear kinematic hardening parameter
      const double kinhard_;
      // parameters F,G,H,L,M,N of Hill yield criterion
      const double rY_11_;
      const double rY_22_;
      const double rY_33_;
      const double rY_12_;
      const double rY_23_;
      const double rY_13_;
      const double plspin_eta_;
      // plastic material parameters ************************************

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();
      //@}

    };  // class PlasticElastHyper

  }  // namespace PAR

class PlasticElastHyperType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "PlasticElastHyperType"; }

  static PlasticElastHyperType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const std::vector<char> & data );

private:

  static PlasticElastHyperType instance_;
};


  /*----------------------------------------------------------------------*/
  /// Collection of hyperelastic materials
  class Material;

  class PlasticElastHyper : public So3Material
  {
  public:

    /// construct empty material object
    PlasticElastHyper();

    /// construct the material object given material parameters
    explicit PlasticElastHyper(MAT::PAR::PlasticElastHyper* params);

    ///@name Packing and Unpacking
    //@{

    /// \brief Return unique ParObject id
    ///
    /// every class implementing ParObject needs a unique id defined at the
    /// top of drt_parobject.H (this file) and should return it in this method.
    virtual int UniqueParObjectId() const { return PlasticElastHyperType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_plelasthyper; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new PlasticElastHyper( *this ) ); }

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

    /// initial yield stress
    virtual double Inityield() const { return params_->inityield_;}

    /// linear isotropic hardening modulus
    virtual double Isohard() const { return params_->isohard_;}

    /// exponent for nonlinear isotropic hardening
    virtual double Expisohard() const {return params_->expisohard_;}

    /// saturation yield stress for nonlinear isotropic hardening
    virtual double Infyield() const { return params_->infyield_;}

    /// linear kinematic hardening modulus
    virtual double Kinhard() const { return params_->kinhard_;}

    /// a shear modulus equivalent
    virtual double ShearMod() const;

    /// a shear modulus equivalent
    virtual bool HaveHillPlasticity() const {return (params_->rY_11_!=0.); }

    /// plastic anisotropy tensor
    virtual LINALG::Matrix<5,5> PlAniso() const {return *PlAniso_; }

    /// inverse plastic anisotropy tensor
    virtual LINALG::Matrix<5,5> InvPlAniso() const {return *InvPlAniso_; }

    /// plastic spin eta
    virtual double PlSpinEta() const {return params_->plspin_eta_; }

    /// hyperelastic stress response plus elasticity tensor
    virtual void Evaluate(const LINALG::Matrix<3,3>* defgrd,  ///< Deformation gradient
                          const LINALG::Matrix<3,3>* invpldefgrd, ///< inverse plastic deformation gradient
                          Teuchos::ParameterList& params,   ///< Container for additional information
                          LINALG::Matrix<6,1>* PK2stress,   ///< 2nd Piola-Kirchhoff stresses
                          LINALG::Matrix<6,6>* cmat_ref,    ///< Constitutive matrix
                          LINALG::Matrix<6,9>* dPK2dFpinv,   ///< derivative of PK2 w.r.t. F_p^{-1}
                          LINALG::Matrix<3,3>* MandelStress, ///< Mandel stress tensor
                          LINALG::Matrix<6,6>* dMdC,         ///< derivative of Mandel stress w.r.t. C
                          LINALG::Matrix<6,9>* dMdFpinv      ///< derivative of Mandel stress w.r.t. F_p^{-1}
                          );

    /// hyperelastic stress response plus elasticity tensor
    virtual void Evaluate(const LINALG::Matrix<3,3>* defgrd,  ///< Deformation gradient
                          const LINALG::Matrix<6,1>* glstrain,///< Green-Lagrange strain
                          Teuchos::ParameterList& params,     ///< Container for additional information
                          LINALG::Matrix<6,1>* stress,        ///< 2nd Piola-Kirchhoff stresses
                          LINALG::Matrix<6,6>* cmat)          ///< Constitutive matrix
    {dserror("PlasticHyperElast cannot be called by standard evaluate routine"); return;}

    virtual void Setup(int numgp, DRT::INPUT::LineDefinition* linedef);

    virtual void SetupHillPlasticity(DRT::INPUT::LineDefinition* linedef);

    virtual void Update();

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
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp);

  protected:
    /// calculates the kinematic quantities and tensors used afterwards
    virtual void EvaluatePlastKinQuant(
        const LINALG::Matrix<3,3>* defgrd,
        const LINALG::Matrix<3,3>* invpldefgrd,
        LINALG::Matrix<6,1>& Cpi,
        LINALG::Matrix<6,1>& CpiCCpi,
        LINALG::Matrix<6,1>& ircg,
        LINALG::Matrix<6,1>& Ce,
        LINALG::Matrix<3,3>& CeM,
        LINALG::Matrix<6,1>& Ce2,
        LINALG::Matrix<6,1>& id2V,
        LINALG::Matrix<3,3>& id2,
        LINALG::Matrix<3,3>& CpiC,
        LINALG::Matrix<3,3>& FpiCe,
        LINALG::Matrix<9,1>& CFpiCei,
        LINALG::Matrix<9,1>& CFpi,
        LINALG::Matrix<3,3>& FpiTC,
        LINALG::Matrix<9,1>& CFpiCe,
        LINALG::Matrix<3,3>& CeFpiTC,
        LINALG::Matrix<3,1>& prinv
    );

    /// convert non-symmetric 2-tensor to 9x1 vector
    void Matrix3x3to9x1(LINALG::Matrix<3,3> A, LINALG::Matrix<9,1>& Out);

    /// calculates the factors associated to the hyperelastic laws
    virtual void EvaluateGammaDelta(
        LINALG::Matrix<3,1> prinv,
        LINALG::Matrix<3,1>& gamma,
        LINALG::Matrix<8,1>& delta
        );

    /// add a contribution to the derivative of the PK2 stress with respect to the plastic deformation gradient
    /// in index notation this method does:
    /// dPK2dFpinv_IJKL += fac*( A_IK * B_JL + A_JK * B_IL )
    /// Note that the result is symmetric within the first two indices
    /// Thus the first two indices are stored in stress-like Voigt notation.
    /// The 2nd index pair has no symmetries and is therefore stored as 9-vector
    /// Be aware that this is NOT symmetric w.r.t. A and B !!!
    void AddtodPK2dFpinv(LINALG::Matrix<6,9>& dPK2dFpinv,
        LINALG::Matrix<3,3> A,
        LINALG::Matrix<3,3> B,
        double fac);

    /// in index notation this does:
    /// dMdC_abCD += fac*( A_Ca * A_Db + A_Cb * A_Da )
    /// result has minor symmetries but no major symmetry
    void AddtodMdC_gamma1(LINALG::Matrix<6,6>& dMdC,
        LINALG::Matrix<3,3> A, double fac);

    /// in index notation this does:
    /// dMdC_abCD += fac*( A_Ca*B_Db + A_Da*B_Cb + A_Db*B_Ca + A_Cb*B_Da )
    /// result has minor symmetries but no major symmetry
    // Symmetric w.r.t. A,B
    void AddtodMdC_gamma2(LINALG::Matrix<6,6>& dMdC,
        LINALG::Matrix<3,3> A, LINALG::Matrix<3,3> B, double fac);

    /// calculates the isotropic stress and elasticity tensor for coupled configuration
    virtual void EvaluateIsotropicPrincPlast(
        LINALG::Matrix<6,1>& stressisoprinc,
        LINALG::Matrix<6,6>& cmatisoprinc,
        LINALG::Matrix<6,9>& dPK2dFpinvIsoprinc,
        LINALG::Matrix<3,3>& MandelStressIsoprinc,
        LINALG::Matrix<6,6>& dMdCisoprinc,
        LINALG::Matrix<6,9>& dMdFpinvIsoprinc,
        LINALG::Matrix<6,1> Cpi,
        LINALG::Matrix<6,1> CpiCCpi,
        LINALG::Matrix<6,1> ircg,
        LINALG::Matrix<6,1> Ce,
        LINALG::Matrix<3,3> CeM,
        LINALG::Matrix<6,1> Ce2,
        LINALG::Matrix<6,1> id2V,
        LINALG::Matrix<3,3> id2,
        LINALG::Matrix<3,3> CpiC,
        LINALG::Matrix<3,3> FpiCe,
        LINALG::Matrix<3,3> Fpi,
        LINALG::Matrix<9,1> CFpiCei,
        LINALG::Matrix<9,1> CFpi,
        LINALG::Matrix<3,3> FpiTC,
        LINALG::Matrix<9,1> CFpiCe,
        LINALG::Matrix<3,3> CeFpiTC,
        LINALG::Matrix<3,1> gamma,
        LINALG::Matrix<8,1> delta
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
    bool isomodvisco_ ;

    /// @name Voigt index arrays
    //@{
    static const int VOIGT6ROW_[];  ///< 6-Voigt row index of corresponding 2-tensor
    static const int VOIGT6COL_[];  ///< 6-Voigt column index of corresponding 2-tensor
    static const int VOIGT3X3SYM_[]; ///< go from symmtric 2-tensor index pair to 6-Voigt index
                                     ///< by [3*i+j] for any i,j=0,1,2
    static const int VOIGT3X3_[3][3];

    //@}

    /// my material parameters
    MAT::PAR::PlasticElastHyper* params_;

    /// map to materials/potential summands
    std::vector<Teuchos::RCP<MAT::ELASTIC::Summand> > potsum_;

    // plastic anisotropy tensor for Hill-plasticity
    // as there is only a dependency on the deviatoric stress, we use the reduced matrix
    // notation S=[S_11 S_22 S_12 S_23 S_13]
    Teuchos::RCP<LINALG::Matrix<5,5> > PlAniso_;
    Teuchos::RCP<LINALG::Matrix<5,5> > InvPlAniso_;

  };

}

#endif

