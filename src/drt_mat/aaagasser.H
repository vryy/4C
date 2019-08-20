/*----------------------------------------------------------------------*/
/*! \file
\brief
Be sure to include 'AAA' in your BACI configurations file!!!
This file contains the routines required for OGDEN-like material model
according to GASSER: "Failure properties of intraluminal thrombus in abdominal
aortic aneurysm under static and pulsating mechanical load", Journal of Vascular
Surgery Volume 48, Number 1, July 2008.

the input line should read:
  MAT 1 MAT_Struct_AAAGasser DENS 0.0001 VOL OgSiMi NUE 0.49 BETA -2.0 CLUM 2.62e-3 CMED 2.13e-3
CABLUM 1.98e-3

\level 3

\maintainer Christoph Schmidt


*/
/*----------------------------------------------------------------------*/
#ifndef AAAGASSER_H_
#define AAAGASSER_H_

#include "matpar_parameter.H"
#include "so3_material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for thrombus material acc. to Gasser [2008]
    class AAAgasser : public Parameter
    {
     public:
      /// standard constructor
      AAAgasser(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~AAAgasser() { ; }

      /// @name material parameters
      //@{

      /// mass density
      const double density_;
      /// Type of volumetric Strain Energy Density ('OgSiMi' for Ogden-Simo-Miehe, 'SuBa' for
      /// Sussman-Bathe, 'SiTa' for Simo-Taylor)
      const std::string* vol_;
      /// poisson's ratio
      const double nue_;
      /// parameter from Holzapfel
      const double beta_;
      /// stiffness parameter (luminal)
      const double Clum_;
      /// stiffness parameter (medial)
      const double Cmed_;
      /// stiffness parameter (abluminal)
      const double Cablum_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class AAAgasser

  }  // namespace PAR

  class AAAgasserType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "AAAgasserType"; }

    static AAAgasserType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static AAAgasserType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Be sure to include 'AAA' in your BACI configurations file!!! <br>
  ///
  /// <h3>Isochoric Material Description:</h3>
  ///
  /// OGDEN-like material according to GASSER [1].
  ///
  /// isochoric strain energy function is given by
  ///
  ///\f[
  ///   \Psi_{iso} = c \sum\limits_{\alpha = 1}^{3} (\bar{\lambda}_{\alpha}^{4} - 1) = \ldots =
  ///          c \left[ (\bar{I}_{\boldsymbol C})^2 - 2 \, \bar{II}_{\boldsymbol C} - 3 \right] =
  ///          c \left[ (I_{\boldsymbol C} \cdot III_{\boldsymbol C}^{-\frac{1}{3}})^2 - 2 \,
  ///          II_{\boldsymbol C} \cdot III_{\boldsymbol C}^{-\frac{2}{3}} - 3 \right]
  ///\f]
  ///
  /// ... and modified to slight compressibility !
  ///
  /// with: \f$ \bar{I}_{\boldsymbol C} = I_{\boldsymbol C} \cdot III_{\boldsymbol
  /// C}^{-\frac{1}{3}}, \, \bar{II}_{\boldsymbol C} = II_{\boldsymbol C} \cdot III_{\boldsymbol
  /// C}^{-\frac{2}{3}}\f$
  ///
  /// symbols:
  ///   /// \f$c\f$ .. position dependend stiffness parameter <br>
  /// \f$\Psi\f$  .. strain energy function <br>
  /// \f$\mathbf{C}\f$  .. right CAUCHY-GREEN tensor <br>
  /// \f$\bar{\lambda}_{\alpha}\f$  .. modified principle streches \f$(\alpha = 1, \ldots ,3)\f$
  /// <br> \f$I_{\mathbf{C}}, II_{\mathbf{C}}, III_{\mathbf{C}}\f$ .. principle invariants of the
  /// right CG tensor <br> \f$\bar{I}_{\mathbf{C}}, \bar{II}_{\mathbf{C}}\f$ .. modified invariants
  /// <br>
  ///   ///
  /// <h3>Volumetric Material Description:</h3>
  ///
  /// material according to OGDEN / SIMO & MIEHE [2].
  ///
  /// volumetric strain energy function is given by
  ///
  ///\f[
  ///   \Psi_{vol}^{OSM} = \frac{\kappa}{\beta^2}(\beta \mbox{ln}J + J^{-\beta} - 1)
  ///\f]
  ///
  /// material according to SUSSMAN / BATHE [3].
  ///
  /// volumetric strain energy function is given by
  ///
  ///\f[
  ///   \Psi_{vol}^{SB} = \frac{\kappa}{2}(J - 1)^2
  ///\f]
  ///
  /// material according to SIMO / TAYLOR [3].
  ///
  /// volumetric strain energy function is given by
  ///
  ///\f[
  ///   \Psi_{vol}^{ST} = \frac{\kappa}{4}\left[(J-1)^2 + (\mbox{ln}J)^2\right]
  ///\f]
  ///
  ///
  /// symbols:
  ///   /// \f$J\f$ .. volume ratio (JACOBIAN determinant) (\f$= III_{\mathbf{C}}^{-\frac{1}{2}}\f$)
  ///   <br>
  /// \f$\kappa\f$  .. dilatational modulus acc. to \f$ \kappa = 24c / (3 - 6\nu) \f$ [4] <br>
  /// \f$\beta\f$ .. material parameter acc. to OGDEN <br>
  ///   ///
  ///
  /// <h3>Reference:</h3>
  /// <ul>
  /// <li>[1] Gasser, T.C. and G&ouml;rg&uuml;l&uuml;, G. and Folkesson, M. and Swedenborg, J. :
  ///     Failure properties of intraluminal thrombus in abdominal aortic
  ///     aneurysm under static and pulsating mechanical load, Journal of
  ///     Vascular Surgery, Volume 48, Number 1, July 2008 </li>
  /// <li>[2] Holzapfel G. A.: Nonlinear Solid Mechanics - A Continuum Approach For Engineering,
  /// (reprint) 2007</li> <li>[3] Doll, S. ; Schweizerhof, K.: On the Development of Volumetric
  /// Strain Energy Functions. In: Journal of Applied Mechanics 67 (2000), March</li> <li>[4] Bonet,
  /// Javier ; Wood, Richard D.: Nonlinear Continuum Mechanics for Finite Element Analysis. 2.
  /// Cambridge University Press, 2008</li>
  /// </ul>
  ///
  /// \author goe
  /// \date 11/09
  class AAAgasser : public So3Material
  {
   public:
    /// empty constructor
    AAAgasser();

    /// constructor with given material parameters
    AAAgasser(MAT::PAR::AAAgasser* params);

    /// @name Packing and Unpacking

    ///  \brief Return unique ParObject id
    ///
    ///  every class implementing ParObject needs a unique id defined at the
    ///  top of drt_parobject.H (this file) and should return it in this method.
    virtual int UniqueParObjectId() const { return AAAgasserType::Instance().UniqueParObjectId(); }

    ///  \brief Pack this class so it can be communicated
    ///
    ///  Resizes the vector data and stores all information of a class in it.
    ///  The first information to be stored in data has to be the
    ///  unique parobject id delivered by UniqueParObjectId() which will then
    ///  identify the exact class on the receiving processor.
    ///
    ///  \param data (in/out): char vector to store class information
    virtual void Pack(DRT::PackBuffer& data) const;

    ///  \brief Unpack data from a char vector into this class
    ///
    ///  The vector data contains all information to rebuild the
    ///  exact copy of an instance of a class on a different processor.
    ///  The first entry in data has to be an integer which is the unique
    ///  parobject id defined at the top of this file and delivered by
    ///  UniqueParObjectId().
    ///
    ///  \param data (in) : vector storing all data to be unpacked into this
    ///  instance.
    virtual void Unpack(const std::vector<char>& data);

    //@}


    /// material mass density
    virtual double Density() const { return params_->density_; }

    /// shear modulus
    double ShearMod() const
    {
      dserror("Cannot provide shear modulus equivalent");
      return 0.0;
    }

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_aaagasser; }

    /// check if element kinematics and material kinematics are compatible
    virtual void ValidKinematics(INPAR::STR::KinemType kinem)
    {
      if (!(kinem == INPAR::STR::kinem_nonlinearTotLag))
        dserror("element and material kinematics are not compatible");
    }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new AAAgasser(*this)); }

    /// <h2>Calculation of location dependent stiffness parameter c in GASSER material
    /// formulation:</h2>
    ///
    /// follows...
    ///
    ///
    /// <h2>the material routine - hyperelastic stress response and elasticity tensor:</h2>
    ///
    /// symbols:
    ///     /// \f$c\f$ .. position dependend stiffness parameter acc. to GASSER <br>
    /// \f$\Psi\f$ .. strain energy function <br>
    /// \f$\mathbf{C}\f$ .. right CAUCHY-GREEN tensor <br>
    /// \f$I_{\mathbf{C}}, II_{\mathbf{C}}, III_{\mathbf{C}}\f$ .. invariants of the right CG tensor
    /// <br> \f$\mathbf{I}\f$ .. identity second order <br> \f$\mathbf{S}\f$ .. second
    /// PIOLA-KIRCHHOFF stress tensor <br> \f$\mathbb{C}\f$ .. fourth order elasticity tensor <br>
    /// \f$\mathbb{S}\f$ .. fourth order identity (see Holzapfel p. 261) <br>
    /// \f$J\f$ .. volume ratio (JACOBIAN determinant) (\f$= III_{\mathbf{C}}^{-\frac{1}{2}}\f$)
    /// <br> \f$\kappa\f$ .. dilatational modulus <br> \f$\beta\f$ .. material parameter acc. to
    /// OGDEN <br>
    ///     ///
    /// <h3>Isochoric Part:</h3>
    ///
    /// Determine isochoric PK2 stress response \f$\mathbf{S}_{iso}\f$ and isochoric material
    /// elasticity tensor \f$\mathbb{C}_{iso}\f$ due to strain energy function according to GASSER
    /// [2008] described in principal invariants.
    ///
    /// The stress response is achieved by collecting the coefficients
    /// \f$\gamma_i \, (i=1, \ldots ,3)\f$ (see Holzapfel p.248):
    ///
    ///\f[
    ///   \gamma_1
    ///   = 2 \, \left( \frac{\partial \Psi_{iso}}{\partial I_{\mathbf{C}}} + I_{\mathbf{C}} \,
    ///   \frac{\partial \Psi_{iso}}{\partial II_{\mathbf{C}}} \right) = \ldots = 0
    ///\f]
    ///
    ///\f[
    ///   \gamma_2
    ///   = -2 \, \frac{\partial \Psi_{iso}}{\partial II_{\mathbf{C}}} = \ldots
    ///   = 4c \, III_{\mathbf{C}}^{-\frac{2}{3}}
    ///\f]
    ///
    ///\f[
    ///   \gamma_3
    ///   = 2 \, III_{\mathbf{C}} \, \frac{\partial \Psi_{iso}}{\partial III_{\mathbf{C}}} = \ldots
    ///   = \frac{4}{3}c \, III_{\mathbf{C}}^{-\frac{2}{3}} \, \left(-I_{\mathbf{C}}^2 + 2 \,
    ///   II_{\mathbf{C}}\right)
    ///\f]
    ///
    /// entering Holzapfel Eq (6.32) results:
    ///\f[
    ///   \mathbf{S}_{iso} = 2 \, \frac{\partial
    ///   \Psi_{iso}(I_{\mathbf{C}},II_{\mathbf{C}},III_{\mathbf{C}})}{\partial \mathbf{C}} =
    ///   \gamma_1\mathbf{I} + \gamma_2\mathbf{C} + \gamma_3\mathbf{C}^{-1}
    ///\f]
    ///
    /// The isochoric elasticity tensor is reached in the same fashion. Firstly, the coefficients
    /// \f$\delta_j \, (j=1, \ldots ,8)\f$ are calculated:
    ///
    ///\f[
    ///   \delta_{1}
    ///   = 4 \, \left( \frac{\partial^2 \Psi_{iso}}{\partial I_{\mathbf{C}}\partial I_{\mathbf{C}}}
    ///           + 2 \, I_{\mathbf{C}} \, \frac{\partial^2 \Psi_{iso}}{\partial
    ///           I_{\mathbf{C}}\partial II_{\mathbf{C}}}
    ///           + \frac{\partial \Psi_{iso}}{\partial II_{\mathbf{C}}}
    ///           + I_{\mathbf{C}}^2 \, \frac{\partial^2 \Psi_{iso}}{\partial
    ///           II_{\mathbf{C}}\partial II_{\mathbf{C}}} \right) = \ldots
    ///   = 0
    ///\f]
    ///
    ///\f[
    ///   \delta_{2}
    ///   = -4 \, \left( \frac{\partial^2 \Psi_{iso}}{\partial I_{\mathbf{C}}\partial
    ///   II_{\mathbf{C}}}
    ///           + I_{\mathbf{C}} \, \frac{\partial^2 \Psi_{iso}}{\partial II_{\mathbf{C}}\partial
    ///           II_{\mathbf{C}}} \right) = \ldots
    ///   = 0
    ///\f]
    ///
    ///\f[
    ///   \delta_{3}
    ///   = 4 \, \left( III_{\mathbf{C}}\frac{\partial^2 \Psi_{iso}}{\partial I_{\mathbf{C}}\partial
    ///   III_{\mathbf{C}}}
    ///           + I_{\mathbf{C}}III_{\mathbf{C}} \, \frac{\partial^2 \Psi_{iso}}{\partial
    ///           II_{\mathbf{C}}\partial III_{\mathbf{C}}} \right) = \ldots
    ///   = 0
    ///\f]
    ///
    ///\f[
    ///   \delta_{4}
    ///   = 4 \, \frac{\partial^2 \Psi_{iso}}{\partial II_{\mathbf{C}}\partial II_{\mathbf{C}}} =
    ///   \ldots = 0
    ///\f]
    ///
    ///\f[
    ///   \delta_{5}
    ///   = -4 III_{\mathbf{C}} \, \frac{\partial^2 \Psi_{iso}}{\partial II_{\mathbf{C}}\partial
    ///   III_{\mathbf{C}}} = \ldots = -\frac{16}{3}c \, III_{\mathbf{C}}^{-\frac{2}{3}}
    ///\f]
    ///
    ///\f[
    ///   \delta_{6}
    ///   = 4 \, \left( III_{\mathbf{C}}\frac{\partial \Psi_{iso}}{\partial III_{\mathbf{C}}}
    ///           + III_{\mathbf{C}}^2 \, \frac{\partial^2 \Psi_{iso}}{\partial
    ///           III_{\mathbf{C}}\partial III_{\mathbf{C}}} \right) = \ldots
    ///   =  \frac{16}{9}c \, III_{\mathbf{C}}^{-\frac{2}{3}} \, \left(I_{\mathbf{C}}^2 - 2 \,
    ///   II_{\mathbf{C}}\right)
    ///\f]
    ///
    ///\f[
    ///   \delta_{7}
    ///   = -4 III_{\mathbf{C}} \, \frac{\partial \Psi_{iso}}{\partial III_{\mathbf{C}}} = \ldots
    ///   = \frac{8}{3}c \, III_{\mathbf{C}}^{-\frac{2}{3}} \, \left(I_{\mathbf{C}}^2 - 2 \,
    ///   II_{\mathbf{C}}\right)
    ///\f]
    ///
    ///\f[
    ///   \delta_{8}
    ///   = -4 \, \frac{\partial \Psi_{iso}}{\partial II_{\mathbf{C}}} = \ldots
    ///   = 8 c \, III_{\mathbf{C}}^{-\frac{2}{3}}
    ///\f]
    ///
    /// then Holzapfel Eq (6.193) is evaluated:
    ///
    ///\f[
    ///   \mathbb{C}_{iso} = 2 \, \frac{\partial \mathbf{S}_{iso}}{\partial \mathbf{C}} = 4 \,
    ///   \frac{\partial^2 \Psi_{iso}}{\partial \mathbf{C}\partial \mathbf{C}} =
    ///\f]
    ///\f[
    ///   = \delta_1\mathbf{I}\otimes\mathbf{I} +
    ///   \delta_2\left(\mathbf{I}\otimes\mathbf{C}+\mathbf{C}\otimes\mathbf{I}\right)
    ///   + \delta_3\left(\mathbf{I}\otimes\mathbf{C}^{-1}+\mathbf{C}^{-1}\otimes\mathbf{I}\right) +
    ///   \delta_4\mathbf{C}\otimes\mathbf{C}
    ///   + \delta_5\left(\mathbf{C}\otimes\mathbf{C}^{-1}+\mathbf{C}^{-1}\otimes\mathbf{C}\right) +
    ///   \delta_6\mathbf{C}^{-1}\otimes\mathbf{C}^{-1}
    ///   + \delta_7\mathbf{C}^{-1}\odot\mathbf{C}^{-1} + \delta_8\mathbb{S}
    ///\f]
    ///
    ///
    /// <h3>Dilatational Part:</h3>
    ///
    /// Determine volumetric PK2 stress response \f$\mathbf{S}_{vol}\f$ and volumetric material
    /// elasticity tensor \f$\mathbb{C}_{vol}\f$ due to selected strain energy function.
    ///
    /// The stress response is achieved by deriving the hydrostatic pressure
    /// \f$p\f$ (see Holzapfel p.230, 6.91):
    ///\f[
    ///   p^{OSM}
    ///   = \frac{\partial \Psi_{vol}(J)}{\partial J}  = \ldots
    ///   = \frac{\kappa}{\beta \, J}(1 - J^{-\beta})
    ///\f]
    ///
    ///\f[
    ///   p^{SB}
    ///   = \frac{\partial \Psi_{vol}(J)}{\partial J}  = \ldots
    ///   = \kappa(J - 1)
    ///\f]
    ///
    ///\f[
    ///   p^{ST}
    ///   = \frac{\partial \Psi_{vol}(J)}{\partial J}  = \ldots
    ///   = \frac{1}{2}\kappa(J - \frac{1}{J}\mbox{ln}J - 1)
    ///\f]
    ///
    /// entering Holzapfel Eq (6.89) yields:
    ///\f[
    ///   \mathbf{S}_{vol} = Jp\mathbf{C}^{-1}
    ///\f]
    ///
    /// For the volumetric elasticity tensor the fictious pressure \f$\tilde{p}\f$
    /// has to be calculated:
    ///
    ///\f[
    ///   \tilde{p}^{OSM}
    ///   = \kappa \, J^{-\beta-1}
    ///\f]
    ///
    ///\f[
    ///   \tilde{p}^{SB}
    ///   = \kappa (2J - 1)
    ///\f]
    ///
    ///\f[
    ///   \tilde{p}^{ST}
    ///   = \frac{1}{2}\kappa (2J + \frac{1}{J} - 1)
    ///\f]
    ///
    /// from Holzapfel Eq (6.166) follows:
    ///\f[
    ///   \mathbb{C}_{vol} = 2 \, \frac{\partial \mathbf{S}_{vol}}{\partial \mathbf{C}} = 2 \,
    ///   \frac{\partial (Jp\mathbf{C}^{-1})}{\partial \mathbf{C}} =
    ///   J\tilde{p}\mathbf{C}^{-1}\otimes\mathbf{C}^{-1} - 2Jp\mathbf{C}^{-1}\odot\mathbf{C}^{-1}
    ///\f]
    ///
    /// <h3>Summation:</h3>
    ///
    /// Addition yields the stress and elasticity tensors respectively:
    ///
    ///\f[
    ///   \mathbf{S} = \mathbf{S}_{iso} + \mathbf{S}_{vol}
    ///\f]
    ///
    ///\f[
    ///   \mathbb{C} = \mathbb{C}_{iso} + \mathbb{C}_{vol}
    ///\f]
    ///
    /// <h3>References</h3>
    /// Holzapfel G. A., Nonlinear Solid Mechanics - A Continuum Approach For Engineering, (reprint)
    /// [2007]
    ///
    /// \author goe
    /// \date 11/09
    virtual void Evaluate(const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<6, 1>* glstrain,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat,
        const int eleGID);


    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

   private:
    /// my material parameters
    MAT::PAR::AAAgasser* params_;
  };
}  // namespace MAT

#endif  // AAAGASSER_H_
