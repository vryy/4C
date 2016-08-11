/*----------------------------------------------------------------------*/
/*!
 \file fluidporo_singlephase.H

 \brief single phase material for multiphase porous flow

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_MAT_FLUIDPORO_SINGLEPHASE_H_
#define SRC_DRT_MAT_FLUIDPORO_SINGLEPHASE_H_



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

  //forward declaration
  class FluidPoroPhaseDof;

    /*----------------------------------------------------------------------*/
    /// material parameters for a single phase of porous multiphase fluid
    ///
    /// This object exists only once for each read fluid.
    class FluidPoroSinglePhase: public Parameter
    {

      public:

        /// standard constructor
        FluidPoroSinglePhase(Teuchos::RCP<MAT::PAR::Material> matdata);

        /// destructor
        virtual ~FluidPoroSinglePhase()
        {
          ;
        }

        /// create material instance of matching type with my parameters
        virtual Teuchos::RCP<MAT::Material> CreateMaterial();

        /// @name material parameters
        //@{

        /// kinematic or dynamic viscosity
        const double viscosity_;
        /// density
        const double density_;
        /// permeability
        const double permeability_;
        /// bulk modulus
        const double bulkmodulus_;
        /// ID of type of degree of freedom
        const int phasedofID_;
        /// type of degree of freedom
        FluidPoroPhaseDof* phasedof_;

        //@}

    }; // class FluidPoroSinglePhase



  } // namespace PAR

  class FluidPoroSinglePhaseType: public DRT::ParObjectType
  {
    public:
      std::string Name() const
      {
        return "FluidPoroSinglePhaseType";
      }

      static FluidPoroSinglePhaseType & Instance()
      {
        return instance_;
      }
      ;

      virtual DRT::ParObject* Create(const std::vector<char> & data);

    private:

      static FluidPoroSinglePhaseType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for a single porous fluid phase within multiphase porous flow
  ///
  /// This object exists (several times) at every element
  class FluidPoroSinglePhase: public Material
  {
    public:

      /// construct empty material object
      FluidPoroSinglePhase();

      /// construct the material object given material parameters
      explicit FluidPoroSinglePhase(MAT::PAR::FluidPoroSinglePhase* params);

      //! @name Packing and Unpacking

      /*!
       \brief Return unique ParObject id

       every class implementing ParObject needs a unique id defined at the
       top of drt_parobject.H (this file) and should return it in this method.
       */
      virtual int UniqueParObjectId() const
      {
        return FluidPoroSinglePhaseType::Instance().UniqueParObjectId();
      }

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
      virtual INPAR::MAT::MaterialType MaterialType() const
      {
        return INPAR::MAT::m_fluidporo_singlephase;
      }

      /// return copy of this material object
      virtual Teuchos::RCP<Material> Clone() const
      {
        return Teuchos::rcp(new FluidPoroSinglePhase(*this));
      }

      /// return viscosity
      double Viscosity() const
      {
        return params_->viscosity_;
      }

      /// return density
      double Density() const
      {
        return params_->density_;
      }

      /// return permeability
      double Permeability() const
      {
        return params_->permeability_;
      }

      /// return bulkmodulus
      double Bulkmodulus() const
      {
        return params_->bulkmodulus_;
      }

      /// return type of degree of freedom
      INPAR::MAT::MaterialType PoroDofType() const;

      /// mark dofs associated with this phase in a given row (=numphase) in a matrix
      void FillDoFMatrix(
          Epetra_SerialDenseMatrix& dofmat,
          int numphase) const;

      /// evaluate saturation of the phase
      double EvaluateSaturation(
          int phasenum,
          const std::vector<double>& state,
          const std::vector<double>& pressure) const;

      /// evaluate the generalized(!) pressure of this phase
      double EvaluateGenPressure(
              int phasenum,
              const std::vector<double>& state) const;

      //! evaluate derivative of saturation with respect to pressure
      double EvaluateDerivOfSaturationWrtPressure(
          int phasenum,
          int doftoderive,
          const std::vector<double>& state) const;

      //! evaluate derivative of degree of freedom with respect to pressure
      double EvaluateDerivOfDofWrtPressure(
              int phasenum,
              int doftoderive,
              const std::vector<double>& state) const;

      /// Return quick accessible material parameter data
      virtual MAT::PAR::Parameter* Parameter() const
      {
        return params_;
      }

    private:

      /// my material parameters
      MAT::PAR::FluidPoroSinglePhase* params_;
  };

}


#endif /* SRC_DRT_MAT_FLUIDPORO_SINGLEPHASE_H_ */
