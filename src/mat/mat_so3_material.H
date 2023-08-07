/*----------------------------------------------------------------------*/
/*! \file
\brief a common base class for all solid materials

\level 1

*/
/*----------------------------------------------------------------------*/

#ifndef MAT_SO3_MATERIAL_H
#define MAT_SO3_MATERIAL_H

#include "mat_material.H"

#include <unordered_map>

#include "inpar_structure.H"

#include "lib_linedefinition.H"

namespace MAT
{
  class So3Material : public Material
  {
   public:
    //! @name Packing and Unpacking

    int UniqueParObjectId() const override = 0;

    void Pack(DRT::PackBuffer& data) const override = 0;

    void Unpack(const std::vector<char>& data) override = 0;

    //@}

    //! @name Evaluation methods

    /*!
     * @brief Evaluate the material law, i.e. the stress tensor and the constitutive tensor
     *
     * @param[in]  defgrad   Deformation gradient
     * @param[in]  glstrain  Green-Lagrange strain
     * @param[in]  params    Container for additional information
     * @param[out] stress    2nd Piola-Kirchhoff stresses
     * @param[out] cmat      Constitutive matrix
     * @param[in]  gp        current gauss point
     * @param[in]  eleGID    Element ID
     */
    virtual void Evaluate(const LINALG::Matrix<3, 3>* defgrad, const LINALG::Matrix<6, 1>* glstrain,
        Teuchos::ParameterList& params, LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat,
        int gp, int eleGID) = 0;

    /// Evaluate material law
    virtual void EvaluateNonLinMass(const LINALG::Matrix<3, 3>* defgrd,
        const LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
        LINALG::Matrix<6, 1>* linmass_disp, LINALG::Matrix<6, 1>* linmass_vel, int gp, int eleGID)
    {
      dserror("Material does not support evaluation of nonlinear mass matrix");
    };

    /// Evaluate strain energy function (for hyperelastic materials only)
    virtual void StrainEnergy(const LINALG::Matrix<6, 1>& glstrain, double& psi, int gp, int eleGID)
    {
      dserror("Material does not support calculation of strain energy");
    };

    /// Evaluate material law for GEMM
    virtual void EvaluateGEMM(LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* stress,
        LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>* cmat,
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain_m,
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain_new,
        LINALG::Matrix<MAT::NUM_STRESS_3D, 1>* glstrain_old, LINALG::Matrix<3, 3>* rcg_new,
        LINALG::Matrix<3, 3>* rcg_old, const int gp, const int eleGID)
    {
      dserror("Material does not support evaluation for GEMM");
    };

    /*!
     * @brief  Evaluate the Cauchy stress contracted with normal and direction vector and
     * its linearizations with given deformation gradient.
     *
     * Cauchy stress is evaluated within this function call. If requested, the required
     * linearizations are calculated. A potential thermal dependency is handled if the temperature
     * is handed in.
     *
     * @param[in] defgrd  deformation gradient (\f[\mathbf{F}\f])
     * @param[in] n       normal vector (\f[\mathbf{n}\f])
     * @param[in] dir     direction vector (\f[\mathbf{v}\f]), can be either normal or tangential
     *                    vector
     * @param[out] cauchy_n_dir  cauchy stress tensor contracted using the vectors n and dir
     *                           (\f[ \mathbf{\sigma} \cdot \mathbf{n} \cdot \mathbf{v} \f])
     * @param[out] d_cauchyndir_dn   derivative of cauchy_n_dir w.r.t. vector n
     *                        (\f[ \frac{ \mathrm{d} \mathbf{\sigma} \cdot \mathbf{n} \cdot
     * \mathbf{v}} { \mathrm{d} \mathbf{n}} \f])
     * @param[out] d_cauchyndir_ddir  derivative of cauchy_n_dir w.r.t. direction vector v
     *                        (\f[ \frac{ \mathrm{d} \mathbf{\sigma} \cdot \mathbf{n} \cdot
     * \mathbf{v}} { \mathrm{d} \mathbf{v}} \f])
     * @param[out] d_cauchyndir_dF    derivative of cauchy_n_dir w.r.t. deformation gradient
     *                        (\f[ \frac{ \mathrm{d} \mathbf{\sigma} \cdot \mathbf{n} \cdot
     * \mathbf{v}} { \mathrm{d} \mathbf{F}} \f])
     * @param[out] d2_cauchyndir_dF2  second derivative of cauchy_n_dir w.r.t. deformation gradient
     *                       (\f[ \frac{\mathrm{d}^2 \mathbf{\sigma} \cdot \mathbf{n} \cdot
     * \mathbf{v}}
     *                                 {\mathrm{d} \mathbf{F} \mathrm{d} \mathbf{F} } \f])
     * @param[out] d2_cauchyndir_dF_dn  second derivative of cauchy_n_dir w.r.t. deformation
     *                                  gradient and normal vector
     *                       (\f[ \frac{\mathrm{d}^2 \mathbf{\sigma} \cdot \mathbf{n} \cdot
     * \mathbf{v}}
     *                                 {\mathrm{d} \mathbf{F} \mathrm{d} \mathbf{n} } \f])
     * @param[out] d2_cauchyndir_dF_ddir  second derivative of cauchy_n_dir w.r.t. deformation
     *                                  gradient and direction vector
     *                       (\f[ \frac{\mathrm{d}^2 \mathbf{\sigma} \cdot \mathbf{n} \cdot
     * \mathbf{v}}
     *                                 {\mathrm{d} \mathbf{F} \mathrm{d} \mathbf{v} } \f])
     * @param[in] gp             Gauss point id
     * @param[in] eleGID         global id of element
     * @param[in] concentration  concentration
     * @param[in] temp           temperature
     * @param[out] d_cauchyndir_dT    derivative of cauchy_n_dir w.r.t. temperature
     *                        (\f[ \frac{ \mathrm{d} \mathbf{\sigma} \cdot \mathbf{n} \cdot
     * \mathbf{v}} { \mathrm{d} T} \f])
     * @param[out] d2_cauchyndir_dF_dT  second derivative of cauchy_n_dir w.r.t. deformation
     *                                  gradient and temperature
     *                       (\f[ \frac{\mathrm{d}^2 \mathbf{\sigma} \cdot \mathbf{n} \cdot
     * \mathbf{v}}
     *                                 {\mathrm{d} \mathbf{F} \mathrm{d} T } \f])
     */
    virtual void EvaluateCauchyNDirAndDerivatives(const LINALG::Matrix<3, 3>& defgrd,
        const LINALG::Matrix<3, 1>& n, const LINALG::Matrix<3, 1>& dir, double& cauchy_n_dir,
        LINALG::Matrix<3, 1>* d_cauchyndir_dn, LINALG::Matrix<3, 1>* d_cauchyndir_ddir,
        LINALG::Matrix<9, 1>* d_cauchyndir_dF, LINALG::Matrix<9, 9>* d2_cauchyndir_dF2,
        LINALG::Matrix<9, 3>* d2_cauchyndir_dF_dn, LINALG::Matrix<9, 3>* d2_cauchyndir_dF_ddir,
        int gp, int eleGID, const double* concentration, const double* temp,
        double* d_cauchyndir_dT, LINALG::Matrix<9, 1>* d2_cauchyndir_dF_dT)
    {
      dserror("EvaluateCauchyNDirAndDerivatives not implemented for this material");
    }

    /*!
     * @brief  Evaluate the derivative of the deformation gradient w.r.t. degree of freedom x
     *
     * @param[in] defgrd         deformation gradient
     * @param[in] concentration  concentration at gauss point
     * @param[out] d_F_dx        derivative of deformation gradient w.r.t. degree of freedom x
     */
    virtual void EvaluateLinearizationOD(
        const LINALG::Matrix<3, 3>& defgrd, double concentration, LINALG::Matrix<9, 1>* d_F_dx)
    {
      dserror("EvaluateLinearizationOD not implemented for this material");
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

    /*!
     * Post setup routine which will be called after all elements were read and set up
     *
     * \brief This method will be called after the input phase to setup the material with
     * input data that has not yet been read during the Setup(int,DRT::INPUT::LineDefinition*) call.
     *
     * @param params Container for additional information passed from the element
     * @param eleGID Global element id
     */
    virtual void PostSetup(Teuchos::ParameterList& params, const int eleGID)
    {
      // Do nothing for simple material models
    }

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
    //! \param gp Gauss point
    //! \param params Container for additional information
    //! \param eleGID Element id
    virtual void Update(LINALG::Matrix<3, 3> const& defgrd, int const gp,
        Teuchos::ParameterList& params,
        int const eleGID){/* Do nothing for simple material models*/};

    virtual void UpdatePrestress(LINALG::Matrix<3, 3> const& defgrd, int const gp,
        Teuchos::ParameterList& params, int const eleGID){
        // only have to be implemented by those materials that have their own prestressing
        // implementation
    };

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
    //@{

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

    /*!
     * \brief Register names of the internal data that should be saved during runtime vtk output
     *
     * \param name_and_size [out] : unordered map of names of the data with the respective vector
     * size
     */
    virtual void RegisterVtkOutputDataNames(
        std::unordered_map<std::string, int>& names_and_size) const
    {
      // do nothing if material does not have vtk runtime output
    }

    /*!
     * \brief Evaluate internal data for every Gauss point saved for output during runtime vtk
     * output
     *
     * \param name [in] : Name of the data to export
     * \param data [out] : NUMGPxNUMDATA Matrix holding the data
     *
     * \return true if data is set by the material, otherwise false
     */
    virtual bool EvaluateVtkOutputData(
        const std::string& name, Epetra_SerialDenseMatrix& data) const
    {
      return false;
    }

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
