/*----------------------------------------------------------------------*/
/*!
\file fouriervar.H
\brief material for heat conduction according to fourier's law with variable conductivity

\level 1
<pre>
\maintainer Sebastian Pröll
</pre>
*/

/*----------------------------------------------------------------------*
 |  definitions                                                         |
 *----------------------------------------------------------------------*/
#ifndef FOURIERVAR_H_
#define FOURIERVAR_H_

/*----------------------------------------------------------------------*
 |  headers                                                             |
 *----------------------------------------------------------------------*/
#include "matpar_parameter.H"
#include "material.H"
#include "consolidation.H"

#include "../drt_lib/drt_parobjectfactory.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for FourierVar material
    ///
    /// <h3>Input line</h3>
    /// MAT 2 THERM_FourierVar CAPA 1737210 CONDUCTNUM 2 CONDUCTFUNC 3 4
    class FourierVar : public Parameter
    {
     public:
      /// standard constructor
      FourierVar(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~FourierVar() { ; }

      /// @name material parameters
      //@{

      /// functions for volumetric heat capacity
      const std::vector<int> capafunct_;
      /// functions for heat conductivity evaluation
      const std::vector<int> conductfunct_;
      // material specifying the consolidation implementation
      const int consolmat_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class FourierVar

  }  // namespace PAR

  class FourierVarType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "FourierVarType"; }

    static FourierVarType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static FourierVarType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// FourierVar material
  ///
  /// This is a an extension of Fourier's isotropic law of heat conduction to
  /// multiple phases and varying heat capacity and conductivity.
  ///
  /// the input line should read something like:
  /// MAT 2 THERM_FourierVar CAPANUM § CAPAFUNCT 5 6 6 CONDUCTNUM 3 CONDUCTFUNCT 4 3 4 CONSOLMAT 3
  ///
  /// \author proell
  /// \date 04/18
  class FourierVar : public Material
  {
   public:
    /// empty constructor
    FourierVar();

    /// constructor with given material parameters
    FourierVar(MAT::PAR::FourierVar* params);

    /// @name Packing and Unpacking
    //@{

    /// Return unique ParObject id
    ///
    ///  every class implementing ParObject needs a unique id defined at the
    ///  top of drt_parobject.H (this file) and should return it in this method.
    virtual int UniqueParObjectId() const { return FourierVarType::Instance().UniqueParObjectId(); }

    /// Pack this class so it can be communicated
    ///
    /// Resizes the vector data and stores all information of a class in it.
    /// The first information to be stored in data has to be the
    /// unique parobject id delivered by UniqueParObjectId() which will then
    /// identify the exact class on the receiving processor.
    virtual void Pack(DRT::PackBuffer& data  ///< (in/out): char vector to store class information
        ) const;

    /// \brief Unpack data from a char vector into this class
    ///
    /// The vector data contains all information to rebuild the
    /// exact copy of an instance of a class on a different processor.
    /// The first entry in data has to be an integer which is the unique
    /// parobject id defined at the top of this file and delivered by
    /// UniqueParObjectId().
    ///
    virtual void Unpack(
        const std::vector<char>& data  ///< vector storing all data to be unpacked into this
    );

    //@}

    /// @name Access material constants
    //@{

    /// temperature-dependent volumetric heat capacity
    double Capacity(const double temperature,  ///< temperature at GP
        const int gp                           ///< phase at GP
    );

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_th_fourier_var; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new FourierVar(*this)); }

    //@}

    // setup rotine to create history variables based on number of Gauss points
    void Setup(const int numgp);

    // update history variables
    void Update();

    /// Main material call to determine heat flux and constitutive tensor
    template <unsigned nsd>
    void Evaluate(const LINALG::Matrix<nsd, 1>& gradtemp, LINALG::Matrix<nsd, nsd>& cmat,
        LINALG::Matrix<nsd, 1>& heatflux, const double temp, const int gp)
    {
      // conductivity tensor
      cmat.Clear();
      double cond = Conductivity(temp, gp);
      for (unsigned i = 0; i < nsd; ++i) cmat(i, i) = cond;

      // heatflux
      heatflux.MultiplyNN(cmat, gradtemp);

      return;
    }
    /// Evaluate time derivative of cmat
    template <unsigned nsd>
    void GetCmat_T(const LINALG::Matrix<nsd, 1>& gradtemp, LINALG::Matrix<nsd, nsd>& cmat,
        const double temp, const int gp)
    {
      // conductivity tensor
      cmat.Clear();
      double cond_deriv = Conductivity_T(temp, gp);
      for (unsigned int i = 0; i < gradtemp.Rows(); ++i) cmat(i, i) = cond_deriv;

      return;
    }

    /// Evaluate thermal derivative of conductivity
    double Capacity_T(const double temperature,  ///< temperature at GP
        const int gp                             ///< current GP
    );

    /// Get a capacity value useable for heat integration method
    double HeatIntegrationCapacity();

    /// get a pointer to the (converged!) consolidation fraction, for postproc
    template <unsigned int nen>
    void NodalPhase(LINALG::Matrix<nen, 1>& nodaltemp,  ///< temperatures at element nodes
        LINALG::Matrix<nen, 1>& nodalphase              ///< (return) fraction of solid material
    )
    {
      consol_->NodalPhase(nodaltemp, nodalphase);
    }


    /// evaluate conductivity and capacity at nodes
    ///
    /// only for post-processing purposes
    template <unsigned int nen>
    void NodalConductivityCapacity(
        LINALG::Matrix<nen, 1>& nodaltemp,  ///< temperatures at element nodes
        LINALG::Matrix<nen, 1>& nodalcond,  ///< (return) conductivity at nodes
        LINALG::Matrix<nen, 1>& nodalcapa   ///< (return) capacity at nodes
    )
    {
      consol_->EvaluateTempDependentFunctionAtNodes(nodaltemp, nodalcond, params_->conductfunct_);
      consol_->EvaluateTempDependentFunctionAtNodes(nodaltemp, nodalcapa, params_->capafunct_);
    }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    /// Get direct access to the Consolidation object
    virtual Teuchos::RCP<MAT::Consolidation> Consolidation() { return consol_; }

   private:
    /// temperature-dependent conductivity
    double Conductivity(const double temperature,  ///< temperature at GP
        const int gp                               ///< current GP
    );

    /// Evaluate thermal derivative of conductivity
    double Conductivity_T(const double temperature,  ///< temperature at GP
        const int gp                                 ///< current GP
    );

    /// my material parameters
    MAT::PAR::FourierVar* params_;

    /// pointer to the consolidation material managing the actual evaluation of function
    Teuchos::RCP<MAT::Consolidation> consol_;

  };  // FourierVar

}  // namespace MAT
/*----------------------------------------------------------------------*/
#endif  // FOURIERVAR_H_
