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
      Teuchos::RCP<const MAT::ELAST::Summand> MaterialById(
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
      std::map<int,Teuchos::RCP<MAT::ELAST::Summand> > potsum_;

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
  /// in which the individual \f$\Psi_i\f$ is implemented as #MAT::ELAST::Summand.
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
    Teuchos::RCP<const MAT::ELAST::Summand> MaterialById(const int id) const { return params_->MaterialById(id); }

    /// material mass density
    double Density() const { return params_->density_; }

    /// hyperelastic stress response plus elasticity tensor
    void Evaluate(
      const LINALG::Matrix<6,1>& glstrain,  ///< Green-Lagrange strain
      LINALG::Matrix<6,6>& cmat,  ///< elasticity matrix
      LINALG::Matrix<6,1>& stress  ///< 2nd Piola--Kirchhoff stress
      );

    void Setup();
    
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
    
    /// my material parameters
    MAT::PAR::ElastHyper* params_;
    
    /// Am I anisotropic?
    bool anisotropic_;

    /// structural tensors in voigt notation for anisotropy
    LINALG::Matrix<6,1>  A1_;
    LINALG::Matrix<6,1>  A2_;
    /// mixed structural tensor (asymmetric) \f$a1 \otimes a2\f$
    LINALG::Matrix<3,3>  A1A2_;
  };

}

#endif

#endif
