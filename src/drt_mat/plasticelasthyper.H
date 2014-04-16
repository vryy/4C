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
                          LINALG::Matrix<6,9>* dMdFpinv,      ///< derivative of Mandel stress w.r.t. F_p^{-1}
                          const int eleGID);                ///< Element GID

    /// hyperelastic stress response plus elasticity tensor
    virtual void Evaluate(const LINALG::Matrix<3,3>* defgrd,  ///< Deformation gradient
                          const LINALG::Matrix<6,1>* glstrain,///< Green-Lagrange strain
                          Teuchos::ParameterList& params,     ///< Container for additional information
                          LINALG::Matrix<6,1>* stress,        ///< 2nd Piola-Kirchhoff stresses
                          LINALG::Matrix<6,6>* cmat,          ///< Constitutive matrix
                          const int eleGID)                    ///< Element GID
    {dserror("PlasticHyperElast cannot be called by standard evaluate routine"); return;}

    virtual void Setup(int numgp, DRT::INPUT::LineDefinition* linedef);

    virtual void SetupHillPlasticity(DRT::INPUT::LineDefinition* linedef);

    virtual void Update()
    {dserror("Elastic summands in PlasticHyperElast are not allowed to have internal variables"
        "that would need an Update-routine!");};

    /// return if isotropic not splitted formulation
    virtual bool IsotropicPrincipal() const { return isoprinc_;}

    /// return if isotropic and splitted formulation
    virtual bool IsotropicModified() const { return isomod_;}

    /// return if anisotropic not splitted formulation
    virtual bool AnisotropicPrincipal() const { return anisoprinc_;}

    /// return if anisotropic and splitted formulation
    virtual bool AnisotropicModified() const { return anisomod_;}

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    /// Return whether the material requires the deformation gradient for its evaluation
    virtual bool NeedsDefgrd() { return true; };

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

