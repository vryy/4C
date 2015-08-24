/*----------------------------------------------------------------------*/
/*!
 \file structporo_reaction_ecm.H

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/


#ifndef STRUCTPORO_REACTION_ECM_H_
#define STRUCTPORO_REACTION_ECM_H_

#include "structporo_reaction.H"

namespace MAT
{

  // forward declaration
  class StructPoroReactionECM;

  namespace PAR
  {

    class StructPoroReactionECM: public PAR::StructPoroReaction
    {

      friend class MAT::StructPoroReactionECM;

      public:

        /// standard constructor
      StructPoroReactionECM(Teuchos::RCP<MAT::PAR::Material> matdata);

        /// destructor
        virtual ~StructPoroReactionECM()
        {
          ;
        }

        /// create material instance of matching type with my parameters
        virtual Teuchos::RCP<MAT::Material> CreateMaterial();

        /// @name material parameters
        //@{

        //@}

    }; // class StructPoroReactionECM

  } // namespace PAR

  class StructPoroReactionECMType: public DRT::ParObjectType
  {
    public:

      std::string Name() const
      {
        return "StructPoroReactionECMType";
      }

      static StructPoroReactionECMType & Instance()
      {
        return instance_;
      }
      ;

      virtual DRT::ParObject* Create(const std::vector<char> & data);

    private:

      static StructPoroReactionECMType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for StructPoroReactionECM material
  ///
  /// This object exists (several times) at every element
  class StructPoroReactionECM: public StructPoroReaction
  {
    public:

      /// construct empty material object
    StructPoroReactionECM();

      /// construct the material object given material parameters
      explicit StructPoroReactionECM(MAT::PAR::StructPoroReactionECM* params);

      //! @name Packing and Unpacking

      /*!
       \brief Return unique ParObject id

       every class implementing ParObject needs a unique id defined at the
       top of drt_parobject.H (this file) and should return it in this method.
       */
      virtual int UniqueParObjectId() const
      {
        return StructPoroReactionECMType::Instance().UniqueParObjectId();
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
        return INPAR::MAT::m_structpororeactionECM;
      }

      /// return copy of this material object
      virtual Teuchos::RCP<Material> Clone() const
      {
        return Teuchos::rcp(new StructPoroReactionECM(*this));
      }

      /// Initialize internal variables
      void Setup
      (
        int numgp,  ///< number of Gauss points
        DRT::INPUT::LineDefinition* linedef
      );

      /// Return quick accessible material parameter data
      virtual MAT::PAR::Parameter* Parameter() const
      {
        return params_;
      }

      /// Update of GP data (e.g., history variables)
      virtual void Update();

      double BodyForceTerm(double porosity) const;

      //! @name Visualization methods

      // Return names of visualization data
      virtual void VisNames(std::map<std::string,int>& names);

      // Return visualization data
      virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID);

      //@}

    protected:

      virtual void Reaction(Teuchos::RCP<std::vector<double> >& scalars,
                            double porosity,
                            Teuchos::ParameterList& params);

      double refporosity_old_;

      double refporositydot_old_;

      ///reference porosity
      double conc_m2c1_;

      double conc_m2c1_old_;

      double conc_m2c1dot_;
      double conc_m2c1dot_old_;

      /// my material parameters
      MAT::PAR::StructPoroReactionECM* params_;

  };

}

#endif /* STRUCTPORO_REACTION_ECM_H_ */
