/*----------------------------------------------------------------------*/
/*!
\brief a common base class for all solid materials

\level 1

\maintainer Fabian Braeu
*/
/*----------------------------------------------------------------------*/

#ifndef SO3MATERIAL_H
#define SO3MATERIAL_H

#include "material.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_inpar/inpar_structure.H"

namespace MAT
{
  class So3Material : public Material
  {
   public:
    //! @name Packing and Unpacking

    /// Return unique ParObject ID
    virtual int UniqueParObjectId() const = 0;

    /// Pack this class so it can be communicated
    virtual void Pack(DRT::PackBuffer& data) const = 0;

    /// Unpack data from a char vector into this class
    virtual void Unpack(const std::vector<char>& data) = 0;

    //@}

    //! @name Evaluation methods

    /// Evaluate material law
    virtual void Evaluate(const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<6, 1>* glstrain,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat,
        const int eleGID) = 0;

    /// Evaluate material law
    virtual void EvaluateNonLinMass(const LINALG::Matrix<3, 3>* defgrd,
        const LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
        LINALG::Matrix<6, 1>* linmass_disp, LINALG::Matrix<6, 1>* linmass_vel, const int eleGID)
    {
      dserror("Material does not support evaluation of nonlinear mass matrix");
    };

    /// Evaluate strain energy function (for hyperelastic materials only)
    virtual void StrainEnergy(const LINALG::Matrix<6, 1>& glstrain, double& psi, const int eleGID)
    {
      dserror("Material does not support calculation of strain energy");
    };

    /// Evaluate material law for GEMM
    virtual void EvaluateGEMM(LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* stress,
        LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>* cmat, double* density,
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain_m,
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain_new,
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain_old, LINALG::Matrix<3, 3>* rcg_new,
        LINALG::Matrix<3, 3>* rcg_old, const int eleGID)
    {
      dserror("Material does not support evaluation for GEMM");
    };

    /**
     * \brief Evaluate the Cauchy stress and linearizations with given deformation gradient.
     *
     * Cauchy stress is evaluated within this function call. If requested, the required
     * linearizations are calculated. A potential thermal dependency is handled if the temperature
     * is handed in.
     *
     * \param[in]     defgrd    deformation gradient (\f[\bold{F}\f])
     * \param[in]     n         vector n (\f[\bold{n}\f])
     * \param[in]     t         vector t (\f[\bold{t}\f])
     * \param[out]    snt       cauchy stress tensor contracted using the vectors n and t (\f[
                                \bold{\sigma} \cdot \bold{n} \cdot \bold{t} \f])
     * \param[out]    DsntDn    derivative of snt w.r.t. vector n (\f[ \frac{\mathrm{d}
                                \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}{\mathrm{d} \bold{n}}
                                \f])
     * \param[out]    DsntDt    derivative of snt w.r.t. vector t (\f[ \frac{\mathrm{d}
                                \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}{\mathrm{d} \bold{t}}
                                \f])
     * \param[out]    DsntDF    derivative of snt w.r.t. deformation gradient (\f[\frac{\mathrm{d}
                                \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}{\mathrm{d}\bold{F}}\f])
     * \param[out]    D2sntDF2  second derivative of snt w.r.t. deformation gradient (\f[ \frac{
                                \mathrm{d}^2 \bold{\sigma} \cdot \bold{n} \cdot \bold{t}}
                                {\mathrm{d} \bold{F}^2} \f])
     * \param[out]    D2sntDFDn second derivative of snt w.r.t. deformation gradient and vector n
                                (\f[ \frac{ \mathrm{d}^2 \bold{\sigma} \cdot \bold{n} \cdot
                                \bold{t}} {\mathrm{d} \bold{F} \mathrm{d} \bold{n} } \f])
     * \param[out]    D2sntDFDt second derivative of snt w.r.t. deformation gradient and vector n
                                (\f[ \frac{ \mathrm{d}^2 \bold{\sigma} \cdot \bold{n} \cdot
                                \bold{t}} {\mathrm{d} \bold{F} \mathrm{d} \bold{n} } \f])
     * \param[in]     eleGID    global ID of element
     * \param[in]     temp      temperature (if thermal dependency shall be handeled)
     * \param[out]    DsntDT    derivative of snt w.r.t. temperature (\f[ \frac{ \mathrm{d}
                                \bold{\sigma} \cdot \bold{n} \cdot \bold{t} }{\mathrm{d} T} \f])
     * \param[out]    D2sntDFDT second derivative of snt w.r.t. deformation gradient and temperature
                                (\f[ \frac{ \mathrm{d}^2 \bold{\sigma} \cdot \bold{n} \cdot \bold{t}
                                } {\mathrm{d} \bold{F} \mathrm{d} T } \f])
     */
    virtual void EvaluateCauchy(const LINALG::Matrix<3, 3>& defgrd, const LINALG::Matrix<3, 1>& n,
        const LINALG::Matrix<3, 1>& t, double& snt, LINALG::Matrix<3, 1>* DsntDn,
        LINALG::Matrix<3, 1>* DsntDt, LINALG::Matrix<9, 1>* DsntDF, LINALG::Matrix<9, 9>* D2sntDF2,
        LINALG::Matrix<9, 3>* D2sntDFDn, LINALG::Matrix<9, 3>* D2sntDFDt, const int eleGID,
        const double* temp = NULL, double* DsntDT = NULL, LINALG::Matrix<9, 1>* D2sntDFDT = NULL)
    {
      dserror("EvaluateCauchy not implementend for this material");
      return;
    }

    //@}

    /// Return whether material includes a varying material density
    virtual bool VaryingDensity() const { return false; }

    //! @name Handling of Gauss point data

    /// check if element kinematics and material kinematics are compatible
    virtual void ValidKinematics(INPAR::STR::KinemType kinem) = 0;

    /// Set up for materials with GP data (e.g., history variables)
    virtual void Setup(int numgp, DRT::INPUT::LineDefinition* linedef){
        /* do nothing for simple material models */};

    /// Update of GP data (e.g., history variables)
    virtual void Update(){/* do nothing for simple material models */};

    //! \brief Indicator, whether the extended update call is used
    //!
    //! Return true, if the material needs the Update(defgrd, gp, params, eleGID) call
    virtual bool UsesExtendedUpdate() { return false; }

    //! \brief Update of GP data (e.g., history variables)
    //!
    //! This method is currently only called
    //! from specific element types. If you need the additional functionality compared to Update()
    //! with any other element than the adapted ones, you need to implement it yourself. Currently
    //! only HEX8 and HEX8FBAR elements are supported
    //!
    //! Materials that use this method need to return true in UsesExtendedUpdate()
    //!
    //! \param defgrd deformation gradient
    //! \param gp Gauß point
    //! \param params Container for additional information
    //! \param eleGID Element id
    virtual void Update(LINALG::Matrix<3, 3> const& defgrd, int const gp,
        Teuchos::ParameterList& params,
        int const eleGID){/* Do nothing for simple material models*/};

    /// Reset time step (for time adaptivity)
    virtual void ResetStep(){/* do nothing for simple material models */};

    /// Reset internal variables to state in the beginning of the computation (needed for inverse
    /// analysis)
    virtual void ResetAll(int numgp){/* do nothing for simple material models */};

    /// Store internal history variables to be eventually reset at some point
    virtual void StoreHistory(int timestep)
    { /* do nothing for simple material models */
    }

    /// Set history variables from time point given as input
    virtual void SetHistory(int timestep)
    { /* do nothing for simple material models */
    }

    //@}

    //! @name Visualization methods

    /// Return names of visualization data
    virtual void VisNames(std::map<std::string, int>& names){
        /* do nothing for simple material models */};

    // Return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp)
    { /* do nothing for simple material models */
      return false;
    };

    /// Return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleId)
    { /* do nothing for simple material models */
      return false;
    };

    //@}

    //! @name Query methods

    /// Return whether the material requires the deformation gradient for its evaluation
    virtual bool NeedsDefgrd()
    { /* usual materials are based on GL-strains */
      return false;
    };

    //@}
  };
}  // namespace MAT

#endif
