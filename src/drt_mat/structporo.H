/*!-----------------------------------------------------------------------*
 \file structporo.H

 <pre>
     Maintainer: Anh-Tu Vuong
                 vuong@lnm.mw.tum.de
                 http://www.lnm.mw.tum.de
                 089 - 289-15264
 </pre>
 *-----------------------------------------------------------------------*/

#ifndef STRUCTPORO_H_
#define STRUCTPORO_H_


#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{

  // forward declaration
  class StructPoro;

  namespace PAR
  {

    class StructPoro: public Parameter
    {

      friend class MAT::StructPoro;

      public:

        /// standard constructor
        StructPoro(Teuchos::RCP<MAT::PAR::Material> matdata);

        /// destructor
        virtual ~StructPoro()
        {
          ;
        }

        /// create material instance of matching type with my parameters
        virtual Teuchos::RCP<MAT::Material> CreateMaterial();

        /// @name material parameters
        //@{

        int matid_;

        double bulkmodulus_;

        double penaltyparameter_;

        double initporosity_;

        //@}

    }; // class StructPoro

  } // namespace PAR

  class StructPoroType: public DRT::ParObjectType
  {
    public:

      std::string Name() const
      {
        return "StructPoroType";
      }

      static StructPoroType & Instance()
      {
        return instance_;
      }
      ;

      virtual DRT::ParObject* Create(const vector<char> & data);

    private:

      static StructPoroType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for StructPoro material
  ///
  /// This object exists (several times) at every element
  class StructPoro: public Material
  {
    public:

      /// construct empty material object
      StructPoro();

      /// construct the material object given material parameters
      explicit StructPoro(MAT::PAR::StructPoro* params);

      //! @name Packing and Unpacking

      /*!
       \brief Return unique ParObject id

       every class implementing ParObject needs a unique id defined at the
       top of drt_parobject.H (this file) and should return it in this method.
       */
      virtual int UniqueParObjectId() const
      {
        return StructPoroType::Instance().UniqueParObjectId();
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
      virtual void Unpack(const vector<char>& data);

      //@}

      /// material type
      virtual INPAR::MAT::MaterialType MaterialType() const
      {
        return INPAR::MAT::m_structporo;
      }

      /// return material
      Teuchos::RCP<MAT::Material> GetMaterial() const
      {
        return mat_;
      }

      /// return material ID
      int MatID() const
      {
        return params_->matid_;
      }

      /// return porosity at gausspoints
      RCP< vector<double> > Porosity() const
      {
        return porosity_;
      }

      //void SetPorosityAtGP(std::vector<double> porosity_gp);

      /// return porosity average
      double PorosityAv() const;

      /// return bulk modulus
      double Bulkmodulus() const
      {
        return params_->bulkmodulus_;
      }

      /// return prenalty parameter
      double Penaltyparameter() const
      {
        return params_->penaltyparameter_;
      }

      /// return initial porosity
      double Initporosity() const
      {
        return params_->initporosity_;
      }

      /// compute current porosity (not used)
      void ComputePorosity(double press, double J,int gp,
          double& porosity,double& dphi_dp,double& dphi_dJ,double& dphi_dJdp,double& dphi_dJJ,double& dphi_dpp) const;

      /// return copy of this material object
      virtual Teuchos::RCP<Material> Clone() const
      {
        return Teuchos::rcp(new StructPoro(*this));
      }

      /// Initialize internal variables
      void Setup
      (
        const int numgp  ///< number of Gauss points
      );

      /// Return quick accessible material parameter data
      virtual MAT::PAR::Parameter* Parameter() const
      {
        return params_;
      }

      /// Return porosity evaluated at some gauss point.           %urrecha 04/12
      double GetPorosityAtGP(int gp) const
      {
        return porosity_->at(gp);
      }

      void   SetDPoroDtAtGP(std::vector<double> dporodt_gp);

      /// Return derivative of porosity wrt time evaluated at some gauss point.           %urrecha 04/12
      double GetDPoroDtAtGP(int gp) const
      {
        return dporodt_->at(gp);
      }

      void   SetGradPorosityAtGP(std::vector< LINALG::Matrix<3,1> > gradporosity_gp);

      /// Return porosity gradient evaluated at some gauss point.  %urrecha 04/12
      LINALG::Matrix<3,1> GetGradPorosityAtGP(int gp) const
      {
        return gradporosity_->at(gp);
      }

    private:

      /// my material parameters
      MAT::PAR::StructPoro* params_;

      /// actual material
      Teuchos::RCP<MAT::Material> mat_;


      RCP<vector<double> > porosity_;
      RCP<vector<double> > dporodt_;
      Teuchos::RCP<vector<LINALG::Matrix<3,1> > > gradporosity_;

      //void SetupPorosity(int numgp);
  };

}

#endif /* STRUCTPORO_H_ */
