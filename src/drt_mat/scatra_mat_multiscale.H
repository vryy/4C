/*----------------------------------------------------------------------*/
/*!
\file scatra_mat_multiscale.H

\brief material for macro-scale elements in multi-scale simulations of scalar transport problems

\level 3

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef SCATRA_MAT_MULTISCALE_H
#define SCATRA_MAT_MULTISCALE_H

#include "scatra_mat.H"

namespace MAT
{
  // forward declaration
  class ScatraMatMultiScaleGP;

  namespace PAR
  {
    //! material parameters
    class ScatraMatMultiScale : public ScatraMat
    {
      public:

        //! constructor
        ScatraMatMultiScale(
            Teuchos::RCP<MAT::PAR::Material> matdata
            );

        //! destructor
        virtual ~ScatraMatMultiScale() { return; };

        //! create material
        virtual Teuchos::RCP<MAT::Material> CreateMaterial();

        //! return name of micro-scale input file
        const std::string MicroInputFileName() const { return microfile_; };

        //! return number of micro-scale discretization
        int MicroDisNum() const { return microdisnum_; };

        //! return porosity
        double Porosity() const { return porosity_; };

        //! return tortuosity
        double Tortuosity() const { return tortuosity_; };

        //! return specific micro-scale surface area
        double A_s() const { return A_s_; };

      private:

        //! @name material parameters
        //@{
        //! name of micro-scale input file
        const std::string microfile_;

        //! number of micro-scale discretization
        const int microdisnum_;

        //! porosity
        const double porosity_;

        //! tortuosity
        const double tortuosity_;

        //! specific micro-scale surface area
        const double A_s_;
        //@}
    };  // class MAT::PAR::ScatraMatMultiScale
  } // namespace PAR


  /*----------------------------------------------------------------------*/
  class ScatraMatMultiScaleType : public DRT::ParObjectType
  {
    public:

      std::string Name() const { return "ScatraMatMultiScaleType"; };

      static ScatraMatMultiScaleType& Instance() { return instance_; };

      DRT::ParObject* Create(const std::vector<char>& data);

    private:

      static ScatraMatMultiScaleType instance_;
  };


  /*----------------------------------------------------------------------*/
  //! material wrapper
  class ScatraMatMultiScale : public ScatraMat
  {
    public:

      //! construct empty material
      ScatraMatMultiScale();

      //! construct material with specific material parameters
      explicit ScatraMatMultiScale(MAT::PAR::ScatraMatMultiScale* params);

      //! @name packing and unpacking
      /*!
        \brief Return unique ParObject id

        Every class implementing ParObject needs a unique id defined at the
        top of drt_parobject.H (this file) and should return it in this method.
      */
      int UniqueParObjectId() const { return ScatraMatMultiScaleType::Instance().UniqueParObjectId(); };

      /*!
        \brief Pack this class so it can be communicated

        Resizes the vector data and stores all information of a class in it.
        The first information to be stored in data has to be the
        unique ParObject ID delivered by UniqueParObjectId() which will then
        identify the exact class on the receiving processor.

        \param data (in/out): char vector to store class information
      */
      void Pack(DRT::PackBuffer& data) const;

      /*!
        \brief Unpack data from a char vector into this class

        The vector data contains all information to rebuild the
        exact copy of an instance of a class on a different processor.
        The first entry in data has to be an integer which is the unique
        parobject id defined at the top of this file and delivered by
        UniqueParObjectId().

        \param data (in) : vector storing all data to be unpacked into this instance.
      */
      void Unpack(const std::vector<char>& data);
      //@}

      //! return material type
      INPAR::MAT::MaterialType MaterialType() const {return INPAR::MAT::m_scatra_multiscale;};

      //! clone material
      Teuchos::RCP<Material> Clone() const {return Teuchos::rcp(new ScatraMatMultiScale(*this));};

      //! initialize multi-scale scalar transport material
      void Initialize(
          const int   ele_id,   //!< macro-scale element ID
          const int   gp_id     //!< macro-scale Gauss point ID
          );

      //! prepare time step on micro scale
      void PrepareTimeStep(
          const int      gp_id,        //!< macro-scale Gauss point ID
          const double   phinp_macro   //!< macro-scale state variable
          ) const;

      //! evaluate multi-scale scalar transport material
      void Evaluate(
          const int      gp_id,          //!< macro-scale Gauss point ID
          const double   phinp_macro,    //!< macro-scale state variable
          double&        q_micro,        //!< micro-scale flux
          double&        dq_dphi_micro   //!< derivative of micro-scale flux w.r.t. macro-scale state variable
          ) const;

      //! update multi-scale scalar transport material
      void Update(
          const int   gp_id   //!< macro-scale Gauss point ID
          ) const;

      //! create output on micro scale
      void Output(
          const int   gp_id   //!< macro-scale Gauss point ID
          ) const;

      //! read restart on micro scale
      void ReadRestart(
          const int   gp_id   //!< macro-scale Gauss point ID
          ) const;

      //! return name of micro-scale input file
      const std::string MicroInputFileName() const { return params_->MicroInputFileName(); };

      //! return number of micro-scale discretization
      int MicroDisNum() const { return params_->MicroDisNum(); };

      //! return porosity
      double Porosity() const { return params_->Porosity(); };

      //! return tortuosity
      double Tortuosity() const { return params_->Tortuosity(); };

      //! return specific micro-scale surface area
      double A_s() const { return params_->A_s(); };

    private:

      //! material parameters
      MAT::PAR::ScatraMatMultiScale* params_;

      //! map between Gauss point ID and Gauss point submaterial
      std::map<int,Teuchos::RCP<ScatraMatMultiScaleGP> > matgp_;
  };
}
#endif
