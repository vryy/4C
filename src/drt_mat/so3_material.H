/*----------------------------------------------------------------------*/
/*!
\file so3_material.H

<pre>
Maintainer: Lena Yoshihara
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef SO3MATERIAL_H
#define SO3MATERIAL_H

#include "material.H"
#include "../drt_lib/drt_linedefinition.H"

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
    virtual void Evaluate(const LINALG::Matrix<3,3>* defgrd,
                          const LINALG::Matrix<6,1>* glstrain,
                          Teuchos::ParameterList& params,
                          LINALG::Matrix<6,1>* stress,
                          LINALG::Matrix<6,6>* cmat,
                          const int eleGID) = 0;

    /// Evaluate material law
    virtual void EvaluateNonLinMass(const LINALG::Matrix<3,3>* defgrd,
                          const LINALG::Matrix<6,1>* glstrain,
                          Teuchos::ParameterList& params,
                          LINALG::Matrix<6,1>* linmass_disp,
                          LINALG::Matrix<6,1>* linmass_vel,
                          const int eleGID)
    { dserror("Material does not support evaluation of nonlinear mass matrix"); };

    /// Evaluate strain energy function (for hyperelastic materials only)
    virtual void StrainEnergy(const LINALG::Matrix<6,1>& glstrain,
                              double& psi)
    { dserror("Material does not support calculation of strain energy"); };

    /// Evaluate material law for GEMM
    virtual void EvaluateGEMM(LINALG::Matrix<MAT::NUM_STRESS_3D,1>* stress,
                              LINALG::Matrix<MAT::NUM_STRESS_3D,MAT::NUM_STRESS_3D>* cmat,
                              double* density,
                              LINALG::Matrix<MAT::NUM_STRESS_3D,1>* glstrain_m,
                              LINALG::Matrix<MAT::NUM_STRESS_3D,1>* glstrain_new,
                              LINALG::Matrix<MAT::NUM_STRESS_3D,1>* glstrain_old,
                              LINALG::Matrix<3,3>* rcg_new,
                              LINALG::Matrix<3,3>* rcg_old,
                              const int eleGID)
    { dserror("Material does not support evaluation for GEMM"); };

    //@}

    /// Return whether material includes a varying material density
    virtual bool VaryingDensity() const
    {
      return false;
    }

    //! @name Handling of Gauss point data

    /// Set up for materials with GP data (e.g., history variables)
    virtual void Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
    { /* do nothing for simple material models */ };

    /// Update of GP data (e.g., history variables)
    virtual void Update()
    { /* do nothing for simple material models */ };

    /// Reset time step (for time adaptivity)
    virtual void ResetStep()
    { /* do nothing for simple material models */ };

    /// Reset internal variables to state in the beginning of the computation (needed for inverse analysis)
    virtual void ResetAll(int numgp)
    { /* do nothing for simple material models */ };

    //@}

    //! @name Visualization methods

    /// Return names of visualization data
    virtual void VisNames(std::map<std::string,int>& names)
    { /* do nothing for simple material models */ };

    // Return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp)
    { /* do nothing for simple material models */ return false; };

    /// Return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleId)
    { /* do nothing for simple material models */ return false; };

    //@}
  };
}

#endif
