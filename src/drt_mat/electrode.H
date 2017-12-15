/*!----------------------------------------------------------------------*/
/*!
\file electrode.H

\brief electrode material carrying concentration and electric potential as degrees of freedom

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef ELECTRODE_H
#define ELECTRODE_H

#include "elchsinglemat.H"

#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  namespace PAR
  {
    //! models for half cell open circuit potential of electrode
    enum OCPModels
    {
      ocp_undefined,
      ocp_redlichkister,
      ocp_taralov,
      ocp_polynomial
    };

    //! parameters for electrode material
    class Electrode : public ElchSingleMat
    {
      public:

        //! constructor
        Electrode(
            Teuchos::RCP<MAT::PAR::Material> matdata
            );

        //! destructor
        virtual ~Electrode(){return;};

        //! create instance of electrode material
        virtual Teuchos::RCP<MAT::Material> CreateMaterial();

        //! @name parameters for electrode material
        //! @{
        //! saturation value of intercalated Lithium concentration
        const double cmax_;

        //! model for half cell open circuit potential
        const OCPModels ocpmodel_;

        //! number of parameters underlying half cell open circuit potential model
        const int ocpparanum_;

        //! parameters underlying half cell open circuit potential model
        const std::vector<double> ocppara_;
        //! @}

      private:

        //! convert string to model for half cell open circuit potential
        OCPModels StringToOCPModel(const std::string& ocpmodelstring) const;
    };  // class MAT::PAR::Electrode
  } // namespace PAR


  /*----------------------------------------------------------------------*/
  class ElectrodeType : public DRT::ParObjectType
  {
    public:

      std::string Name() const {return "ElectrodeType";};

      static ElectrodeType& Instance(){return instance_;};

      virtual DRT::ParObject* Create(const std::vector<char>& data);

    private:

      static ElectrodeType instance_;
  };


  /*----------------------------------------------------------------------*/
  //! wrapper for electrode material
  class Electrode : public ElchSingleMat
  {
    public:

      //! construct empty electrode material
      Electrode();

      //! construct electrode material with specific material parameters
      explicit Electrode(MAT::PAR::Electrode* params);

      //! @name packing and unpacking
      /*!
        \brief Return unique ParObject id

        Every class implementing ParObject needs a unique id defined at the
        top of drt_parobject.H (this file) and should return it in this method.
      */
      virtual int UniqueParObjectId() const {return ElectrodeType::Instance().UniqueParObjectId();};

      /*!
        \brief Pack this class so it can be communicated

        Resizes the vector data and stores all information of a class in it.
        The first information to be stored in data has to be the
        unique ParObject ID delivered by UniqueParObjectId() which will then
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

        \param data (in) : vector storing all data to be unpacked into this instance.
      */
      virtual void Unpack(const std::vector<char>& data);
      //@}

      //! return material type
      virtual INPAR::MAT::MaterialType MaterialType() const {return INPAR::MAT::m_electrode;};

      //! clone electrode material
      virtual Teuchos::RCP<Material> Clone() const {return Teuchos::rcp(new Electrode(*this));};

      //! return saturation value of intercalated Lithium concentration
      double CMax() const {return params_->cmax_;};

      //! compute half cell open circuit potential
      double ComputeOpenCircuitPotential(
          const double concentration,   //!< concentration
          const double faraday,         //!< Faraday constant
          const double frt              //!< factor F/RT
          ) const;

      //! compute first derivative of half cell open circuit potential with respect to concentration
      double ComputeFirstDerivOpenCircuitPotential(
          const double concentration,   //!< concentration
          const double faraday,         //!< Faraday constant
          const double frt              //!< factor F/RT
          ) const;

      //! compute second derivative of half cell open circuit potential with respect to concentration
      double ComputeSecondDerivOpenCircuitPotential(
          const double concentration,   //!< concentration
          const double faraday,         //!< Faraday constant
          const double frt              //!< factor F/RT
          ) const;

    protected:

      //! return material parameters
      virtual MAT::PAR::Parameter* Parameter() const {return params_;}

    private:

      //! my material parameters
      MAT::PAR::Electrode* params_;
  };
}
#endif
