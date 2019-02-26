/*!-----------------------------------------------------------------------*
 \file fluidporo.H

\brief  fluid material for poroelasticity problems

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249

\level 2
 *-----------------------------------------------------------------------*/

#ifndef FLUIDPORO_H_
#define FLUIDPORO_H_

/*---------------------------------------------------------------------*
 | headers                                                              |
 *---------------------------------------------------------------------*/
#include "material.H"
#include "matpar_parameter.H"
#include "../drt_lib/drt_parobjectfactory.H"

/*---------------------------------------------------------------------*
 *---------------------------------------------------------------------*/
namespace MAT
{
  namespace PAR
  {
    enum PoroFlowType
    {
      undefined,
      darcy,
      darcy_brinkman
    };

    enum PoroFlowPermeabilityFunction
    {
      pf_undefined,
      const_,
      kozeny_karman
    };

    /*----------------------------------------------------------------------*/
    /// material parameters for Permeable fluid
    ///
    /// This object exists only once for each read Newton fluid.
    class FluidPoro : public Parameter
    {
     public:
      /// standard constructor
      FluidPoro(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~FluidPoro() { ; }

      /// set inital porosity from structural material and calculate permeabilitycorrectionfactor_
      void SetInitialPorosity(double initialporosity);

      /// @name material parameters
      //@{

      /// kinematic or dynamic viscosity
      const double viscosity_;
      /// density
      const double density_;
      /// permeability
      const double permeability_;
      /// flow type: Darcy or Darcy-Brinkman
      PoroFlowType type_;
      /// flag indicating varying permeablity
      const bool varyingpermeability_;
      /// permeability function type: Const or Kozeny_Carman
      PoroFlowPermeabilityFunction permeabilityfunc_;
      /// a correction factor to ensure the permeability set in the input file
      double permeabilitycorrectionfactor_;
      /// initial porosity
      double initialporosity_;

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class FluidPoro

  }  // namespace PAR

  class FluidPoroType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "FluidPoroType"; }

    static FluidPoroType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static FluidPoroType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for poro fluid flow material
  ///
  /// This object exists (several times) at every element
  class FluidPoro : public Material
  {
   public:
    /// construct empty material object
    FluidPoro();

    /// construct the material object given material parameters
    explicit FluidPoro(MAT::PAR::FluidPoro* params);

    //! @name Packing and Unpacking

    /*!
     \brief Return unique ParObject id

     every class implementing ParObject needs a unique id defined at the
     top of drt_parobject.H (this file) and should return it in this method.
     */
    virtual int UniqueParObjectId() const { return FluidPoroType::Instance().UniqueParObjectId(); }

    /*!
     \brief Pack this class so it can be communicated

     Resizes the vector data and stores all information of a class in it.
     The first information to be stored in data has to be the
     unique parobject id delivered by UniqueParObjectId() which will then
     identify the exact class on the receiving processor.

     \param data (in/out): char vector to store class information
     */
    virtual void Pack(DRT::PackBuffer& data) const;

    /*!
     \brief Unpack data from a char vector into this class

     The vector data contains all information to rebuild the
     exact copy of an instance of a class on a different processor.
     The first entry in data has to be an integer which is the unique
     parobject id defined at the top of this file and delivered by
     UniqueParObjectId().

     \param data (in) : vector storing all data to be unpacked into this
     instance.
     */
    virtual void Unpack(const std::vector<char>& data);

    //@}

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_fluidporo; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp(new FluidPoro(*this)); }

    /// compute reaction coefficient
    double ComputeReactionCoeff() const;

    /// compute reaction coefficient
    void ComputeReactionTensor(
        LINALG::Matrix<2, 2>& reactiontensor, double J, double porosity) const;

    /// compute reaction coefficient
    void ComputeReactionTensor(
        LINALG::Matrix<3, 3>& reactiontensor, double J, double porosity) const;

    void ComputeLinMatReactionTensor(LINALG::Matrix<2, 2>& linreac_dphi,
        LINALG::Matrix<2, 2>& linreac_dJ, double J, double porosity) const;

    void ComputeLinMatReactionTensor(LINALG::Matrix<3, 3>& linreac_dphi,
        LINALG::Matrix<3, 3>& linreac_dJ, double J, double porosity) const;

    /// correction factor to ensure the permeability set in the input file
    double PermeabilityCorrectionFactor() const { return params_->permeabilitycorrectionfactor_; }

    /// effective viscosity (zero for Darcy and greater than zero for Darcy-Brinkman)
    double EffectiveViscosity() const;

    /// return type
    PAR::PoroFlowType Type() const { return params_->type_; }

    /// return viscosity
    double Viscosity() const { return params_->viscosity_; }

    /// return density
    double Density() const { return params_->density_; }

    /// return permeability
    double Permeability() const { return params_->permeability_; }

    /// return permeabilityfunction
    PAR::PoroFlowPermeabilityFunction PermeabilityFunction() const
    {
      return params_->permeabilityfunc_;
    }

    /// evaluate material law
    void EvaluateViscStress(LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat,
        const LINALG::Matrix<6, 1>* glstrain, const int gp, Teuchos::ParameterList& params);

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    // flag indicating a varying permeability
    bool VaryingPermeablity() const { return params_->varyingpermeability_; }

   private:
    /// my material parameters
    MAT::PAR::FluidPoro* params_;
  };

}  // namespace MAT

#endif  // FluidPoro_H_
