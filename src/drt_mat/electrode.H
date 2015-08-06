/*!----------------------------------------------------------------------*/
/*!
\file electrode.H

\brief electrode material carrying concentration and electric potential as degrees of freedom

<pre>
Maintainer: Rui Fang
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

    protected:

      //! return material parameters
      virtual MAT::PAR::Parameter* Parameter() const {return params_;}

    private:

      //! my material parameters
      MAT::PAR::Electrode* params_;
  };
}
#endif
