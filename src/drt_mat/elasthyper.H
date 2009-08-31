/*----------------------------------------------------------------------*/
/*!
\file elasthyper.H

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef ELASTHYPER_H
#define ELASTHYPER_H

#ifdef CCADISCRET

#include "matpar_parameter.H"
#include "material.H"
#include "../drt_matelast/elast_summand.H"

namespace MAT
{

  // forward declaration
  class ElastHyper;

  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// Collection of hyperelastic materials
    ///
    /// Storage map of hyperelastic summands.
    ///
    /// \author rausch,tk,bborn
    /// \date 05/09
    class ElastHyper : public Parameter
    {

      friend class MAT::ElastHyper;

    public:

      /// standard constructor
      ///
      /// This constructor recursively calls the constructors of the
      /// parameter sets of the hyperelastic summands.
      ElastHyper(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~ElastHyper() { ; }

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

      double gamma_;

   private:

      /// map to materials/potential summands
      std::map<int,Teuchos::RCP<MAT::ELASTIC::Summand> > potsum_;

      //@}

    };  // class ElastHyper

  }  // namespace PAR


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
  class ElastHyper : public Material
  {
  public:

    /// construct empty material object
    ElastHyper();

    /// construct the material object given material parameters
    explicit ElastHyper(MAT::PAR::ElastHyper* params);

    ///@name Packing and Unpacking
    //@{

    /// \brief Return unique ParObject id
    ///
    /// every class implementing ParObject needs a unique id defined at the
    /// top of drt_parobject.H (this file) and should return it in this method.
    virtual int UniqueParObjectId() const { return ParObject_ElastHyper; }

    /// \brief Pack this class so it can be communicated
    ///
    /// Resizes the vector data and stores all information of a class in it.
    /// The first information to be stored in data has to be the
    /// unique parobject id delivered by UniqueParObjectId() which will then
    /// identify the exact class on the receiving processor.
    ///
    /// \param data (in/out): char vector to store class information
    virtual void Pack(vector<char>& data) const;

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
    virtual void Unpack(const vector<char>& data);

    //@}

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_elasthyper; }

    /// number of materials
    int NumMat() const { return params_->nummat_; }

    /// deliver material ID of index i'th potential summand in collection
    int MatID(
      const unsigned index  ///< index i
      ) const;

    /// provide access to material by its ID
    Teuchos::RCP<const MAT::ELASTIC::Summand> MaterialById(const int id) const { return params_->MaterialById(id); }

    /// material mass density
    double Density() const { return params_->density_; }

    /// a shear modulus equivalent
    double ShearMod() const;

    /// hyperelastic stress response plus elasticity tensor
    void Evaluate(
      const LINALG::Matrix<6,1>& glstrain,  ///< Green-Lagrange strain
      LINALG::Matrix<6,6>& cmat,  ///< elasticity matrix
      LINALG::Matrix<6,1>& stress  ///< 2nd Piola--Kirchhoff stress
      );

    void Setup(DRT::INPUT::LineDefinition* linedef);

    /// return if anisotropic
    bool Anisotropic() const { return anisotropic_;}

    /// return a1s
    LINALG::Matrix<3,1> Geta1() const {return a1_;}

    /// return a2s
    LINALG::Matrix<3,1> Geta2() const {return a2_;}

  private:

    /// calculate principal invariants
    void InvariantsPrincipal(
      LINALG::Matrix<3,1>& prinv,  ///< principal invariants
      const LINALG::Matrix<6,1>& rcg  ///< symmetric Cartesian 2-tensor in strain-like 6-Voigt notation
      );

    /// calculate principal invariants
    void InvariantsPrincipalAniso(
      LINALG::Matrix<6,1>& prinv,  ///< principal invariants
      const LINALG::Matrix<6,1>& rcg  ///< symmetric Cartesian 2-tensor in strain-like 6-Voigt notation
      );

    /// calculate modified invariants
    void InvariantsModified(
      LINALG::Matrix<3,1>& modinv,  ///< modified invariants
      const LINALG::Matrix<3,1>& prinv  ///< principal invariants
      );

    /// calculate principal stretches \f$\lambda_\alpha\f$ and material principal directions
    void StretchesPrincipal(
      LINALG::Matrix<3,1>& prstr,  ///< principal stretches
      LINALG::Matrix<3,3>& prdir,  ///< principal directions, stored col-by-col
      const LINALG::Matrix<6,1>& rcg  ///< symmetric Cartesian 2-tensor in strain-like 6-Voigt notation
      );

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
    void ResponseStretches(
      LINALG::Matrix<6,6>& cmat,  ///< material elasticity matrix
      LINALG::Matrix<6,1>& stress,  ///< 2nd Piola-Kirchhoff stress vector
      const LINALG::Matrix<6,1>& rcg,  ///< right Cauchy-Green strain vector
      const bool& havecoeffstrpr,  ///< true, if at least one summand is based on pr.stretches
      const bool& havecoeffstrmod  ///< true, if at least one summand is based on mod.pr.stretches
      );

    /// my material parameters
    MAT::PAR::ElastHyper* params_;

    /// Am I anisotropic?
    bool anisotropic_;

    /// structural tensors in voigt notation for anisotropy
    LINALG::Matrix<6,1>  A1_;
    LINALG::Matrix<6,1>  A2_;
    /// vectors of local coordinate system
    LINALG::Matrix<3,1>  a1_;
    LINALG::Matrix<3,1>  a2_;
    /// mixed structural tensor (asymmetric) \f$a1 \otimes a2\f$
    LINALG::Matrix<3,3>  A1A2_;

    /// @name Voigt index arrays
    //@{
    static const int VOIGT6ROW_[];  ///< 6-Voigt row index of corresponding 2-tensor
    static const int VOIGT6COL_[];  ///< 6-Voigt column index of corresponding 2-tensor
    static const int VOIGT3X3SYM_[]; ///< go from symmtric 2-tensor index pair to 6-Voigt index
                                     ///< by [3*i+j] for any i,j=0,1,2
    //@}

  };

}

#endif

#endif
