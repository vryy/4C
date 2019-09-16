/*----------------------------------------------------------------------*/
/*! \file
\brief This file contains the routines required to calculate the isochoric contribution
of the aaagasser material and the corresponding volumetric contribution.

MAT 20 ELAST_IsoVolAAAGasser CLUM 2.62E3 CMED 1.98E3 CABLUM 1.73E3 NUE 0.49 BETA -2.0

\level 1

\maintainer Amadeus Gebauer

*----------------------------------------------------------------------*/
#ifndef ELAST_ISOVOLAAAGASSER_H_
#define ELAST_ISOVOLAAAGASSER_H_

#include "elast_summand.H"
#include "../drt_mat/matpar_parameter.H"

namespace MAT
{
  namespace ELASTIC
  {
    namespace PAR
    {
      /*----------------------------------------------------------------------*/
      /// This material provides the isochoric and(!) the volumetric
      /// contribution for a AAA thromubs. The input line should read
      /// like this:
      ///     MAT 21 ELAST_IsoVolAAAGasser CLUM 2.62E3 CMED 1.98E3
      ///     CABLUM 1.73E3 NUE 0.49 BETA -2.0
      ///
      /// An Ogden type material
      /// \f$\Psi=c\underset{i=1}{\overset{3}{\sum}}(\lambda_{i}^{4}-1)\f$ is chosen for the
      /// isochoric contribution with spatially varying stiffness parameters governed by luminal
      /// thrombus layer stiffness CLUM, the medial layer stiffness CMED and the abluminal thrombus
      /// layer stiffness CABLUM (see: Gasser et al., Failure properties of intraluminal thrombus in
      /// abdominal aortic aneurysm under static and pulsating
      /// mechanical loads, J Vasc Surg, 2008). The volumetric
      /// contribution is modeled by an Ogden-Simo_Miehe type SEF:
      /// \f$\Psi=\frac {\kappa}{\beta^2}(\beta lnJ + J^{-\beta}-1)\f$ (see  Doll and Schweizerhof,
      /// On the Development of Volumetric Strain Energy Functions, Journal of Applied Mechanics,
      /// 2000).
      ///
      /// Strain energy function is given by
      ///\f[
      ///   \Psi = c (\overline{I}_1^2 - 2 \overline{I}_2 -3)
      ///          + \frac{8c}{(1-2 \nu) \beta^2} \big(\beta log(J) + J^{-\beta} -1 \big)
      ///\f]


      class IsoVolAAAGasser : public MAT::PAR::Parameter
      {
       public:
        /// standard constructor
        IsoVolAAAGasser(Teuchos::RCP<MAT::PAR::Material> matdata);

        /// destructor
        virtual ~IsoVolAAAGasser() { ; }

        /// @name material parameters
        //@{

        // !brief enum for mapping between material parameter and entry in the matparams_ vector
        enum matparamnames_
        {
          clum,      /// stiffness parameter (luminal)
          cmed,      /// stiffness parameter (medial)
          cablum,    /// stiffness parameter (abluminal)
          nue,       /// poisson ration
          beta,      /// numerical parameter (should be -2.0 (Doll and Schwizerhof)
                     /// or 9.0(Holzapfel and Ogden))
          normdist,  /// normalized ILT thickness
          cele,
          mu_lum,
          mu_med,
          mu_ablum,
          sigma_lum,
          sigma_med,
          sigma_ablum,
          xi,
          first = clum,
          last = xi
        };

        virtual void OptParams(std::map<std::string, int>* pnames);

        //@}

        /// create material instance of matching type with my parameters
        ///
        /// Remark:                                                  AN 03/2014
        /// Function has to be implemented, since it is a pure virtual
        /// function in the MAT::PAR::Material class.
        /// Function stays empty, since the call within the material factory is done only
        /// respective the elasthyper class
        ///
        virtual Teuchos::RCP<MAT::Material> CreateMaterial() { return Teuchos::null; };

        void SetInitToTrue() { isinit_ = true; };

        bool IsInit() { return isinit_; };

       private:
        // has this material been properly initialized?
        double isinit_;

      };  // class IsoVolAAAGasser

    }  // namespace PAR

    class IsoVolAAAGasser : public Summand
    {
     public:
      /// constructor with given material parameters
      IsoVolAAAGasser(MAT::ELASTIC::PAR::IsoVolAAAGasser* params);

      ///@name Packing and Unpacking
      //@{

      virtual void PackSummand(DRT::PackBuffer& data) const;

      virtual void UnpackSummand(
          const std::vector<char>& data, std::vector<char>::size_type& position);

      //@}

      /// @name Access material constants
      //@{

      /// material type
      virtual INPAR::MAT::MaterialType MaterialType() const
      {
        return INPAR::MAT::mes_isovolaaagasser;
      }

      //@}

      void CalcCele(const int eleID);

      /// Setup of patient-specific stuff
      virtual void SetupAAA(Teuchos::ParameterList& params, const int eleGID);

      // add strain energy
      virtual void AddStrainEnergy(double& psi,  ///< strain energy function
          const LINALG::Matrix<3, 1>& prinv,  ///< principal invariants of right Cauchy-Green tensor
          const LINALG::Matrix<3, 1>& modinv,  ///< modified invariants of right Cauchy-Green tensor
          const LINALG::Matrix<6, 1>& glstrain,  ///< Green-Lagrange strain
          const int eleGID                       ///< element GID
      );


      // birzle 12/2014
      // Add derivatives with respect to modified invariants.
      virtual void AddDerivativesModified(
          LINALG::Matrix<3, 1>& dPmodI,  ///< first derivative with respect to modified invariants
          LINALG::Matrix<6, 1>&
              ddmodPII,  ///< second derivative with respect to modified invariants
          const LINALG::Matrix<3, 1>& modinv,  ///< modified invariants of right Cauchy-Green tensor
          const int eleGID                     ///< element GID
      );

      /// Indicator for formulation
      virtual void SpecifyFormulation(
          bool& isoprinc,     ///< global indicator for isotropic principal formulation
          bool& isomod,       ///< global indicator for isotropic splitted formulation
          bool& anisoprinc,   ///< global indicator for anisotropic principal formulation
          bool& anisomod,     ///< global indicator for anisotropic splitted formulation
          bool& viscogeneral  ///< global indicator, if one viscoelastic formulation is used
      )
      {
        isomod = true;
        return;
      };

      //! @name Visualization methods

      /// Return names of visualization data
      virtual void VisNames(std::map<std::string, int>& names);

      /// Return visualization data
      virtual bool VisData(
          const std::string& name, std::vector<double>& data, int numgp, int eleId);

      //@}

     private:
      /// my material parameters
      MAT::ELASTIC::PAR::IsoVolAAAGasser* params_;

    };  // class IsoVolAAAGasser

  }  // namespace ELASTIC
}  // namespace MAT

#endif
