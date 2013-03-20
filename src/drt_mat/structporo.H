/*!-----------------------------------------------------------------------*
 \file structporo.H

 \brief wrapper for structure material of porous media

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
#include "so3_material.H"
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

        /// material ID of sub-material
        int matid_;

        /// bulk modulus of skeleton phase
        double bulkmodulus_;

        /// penalty parameter for porosity
        double penaltyparameter_;

        /// initial porosity
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

      virtual DRT::ParObject* Create(const std::vector<char> & data);

    private:

      static StructPoroType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for StructPoro material
  ///
  /// This object exists (several times) at every element
  class StructPoro: public So3Material
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
      virtual void Unpack(const std::vector<char>& data);

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
      Teuchos::RCP<std::vector<double> > Porosity() const
      {
        return porosity_;
      }

      /// return porosity at gausspoints of surface element
      std::vector<double> SurfPorosity(int surfnum) const
      {
        dsassert(surfporosity_!=Teuchos::null,"surface porosity has not been evaluated");
        return (*surfporosity_)[surfnum];
      }

      /// return porosity average
      double PorosityAv() const;

      /// return porosity gradient average
      //void PorosityGradientAv( LINALG::Matrix<3,1>& porosityav) const;
      //void PorosityGradientAv( LINALG::Matrix<2,1>& porosityav) const;

      /// return porosity time derivative average
      //double DPorosityDtAv() const;

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

      /// return time deriviative of reference porosity (only nonzero with reaction)
      virtual double RefPorosityTimeDeriv() const
      {
        return 0.0;
      }

      /// compute current porosity and save it
      virtual void ComputePorosity(  Teuchos::ParameterList& params, ///< (i) element parameter list
                             double press,                  ///< (i) pressure at gauss point
                             double J,                      ///< (i) determinant of jacobian at gauss point
                             int gp,                        ///< (i) number of current gauss point
                             double& porosity,              ///< (o) porosity at gauss point
                             double* dphi_dp,               ///< (o) first derivative of porosity w.r.t. pressure at gauss point
                             double* dphi_dJ,               ///< (o) first derivative of porosity w.r.t. jacobian at gauss point
                             double* dphi_dJdp,             ///< (o) derivative of porosity w.r.t. pressure and jacobian at gauss point
                             double* dphi_dJJ,              ///< (o) second derivative of porosity w.r.t. jacobian at gauss point
                             double* dphi_dpp,               ///< (o) second derivative of porosity w.r.t. pressure at gauss point
                             bool save = true
                             ) ;

      /// compute current porosity and save it
      void ComputePorosity(  Teuchos::ParameterList& params, ///< (i) element parameter list
                             double press,       ///< (i) pressure at gauss point
                             double J,           ///< (i) determinant of jacobian at gauss point
                             int gp,             ///< (i) number of current gauss point
                             double& porosity,    ///< (o) porosity at gauss point
                             bool save = true
                             ) ;

      /// compute current surface porosity and save it
      void ComputeSurfPorosity( Teuchos::ParameterList& params, ///< (i) element parameter list
                                double press,       ///< (i) pressure at gauss point
                                double J,           ///< (i) determinant of jacobian at gauss point
                                const int surfnum,  ///< (i) number of surface
                                int gp,             ///< (i) number of current gauss point
                                double& porosity,   ///< (o) porosity at gauss point
                                double* dphi_dp,    ///< (o) first derivative of porosity w.r.t. pressure at gauss point
                                double* dphi_dJ,    ///< (o) first derivative of porosity w.r.t. jacobian at gauss point
                                double* dphi_dJdp,  ///< (o) derivative of porosity w.r.t. pressure and jacobian at gauss point
                                double* dphi_dJJ,   ///< (o) second derivative of porosity w.r.t. jacobian at gauss point
                                double* dphi_dpp,    ///< (o) second derivative of porosity w.r.t. pressure at gauss point
                                bool save = true
                                ) ;

      /// compute current surface porosity and save it
      void ComputeSurfPorosity( Teuchos::ParameterList& params, ///< (i) element parameter list
                                double press,       ///< (i) pressure at gauss point
                                double J,           ///< (i) determinant of jacobian at gauss point
                                const int surfnum,  ///< (i) number of surface
                                int gp,             ///< (i) number of current gauss point
                                double& porosity,    ///< (o) porosity at gauss point
                                bool save = true
                                ) ;

      /// return copy of this material object
      virtual Teuchos::RCP<Material> Clone() const
      {
        return Teuchos::rcp(new StructPoro(*this));
      }

      /// Initialize internal variables
      virtual void PoroSetup
      (
        int numgp,  ///< number of Gauss points
        DRT::INPUT::LineDefinition* linedef
      );

      //Calculate coupling part of homogenized 2 Piola-Kirchhoff stress
      //3D
      void CouplStress( const LINALG::Matrix<3,3>& defgrad,   ///< (i) deformation gradient
                        const LINALG::Matrix<3,1>& fluidvel,  ///< (i) fluid velocity in gausspoint
                        const double& press,                  ///< (i) pressure in gauss point
                        LINALG::Matrix<6,1>& couplstress      ///< (o) coupling stress in gauss point
                        ) const;
      //2D
      void CouplStress( const LINALG::Matrix<2,2>& defgrad,   ///< (i) deformation gradient
                        const LINALG::Matrix<2,1>& fluidvel,  ///< (i) fluid velocity in gausspoint
                        const double& press,                  ///< (i) pressure in gauss point
                        LINALG::Matrix<3,1>& couplstress      ///< (o) coupling stress in gauss point
                        ) const;

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

      /*
      void SetDPoroDtAtGP(std::vector<double> dporodt_gp);

      /// Return derivative of porosity wrt time evaluated at some gauss point.           %urrecha 04/12
      double GetDPoroDtAtGP(int gp) const
      {
        return dporodt_->at(gp);
      }

      /// Set porosity gradient evaluated at gauss points
      void   SetGradPorosityAtGP(std::vector< LINALG::Matrix<3,1> > gradporosity_gp);
      void   SetGradPorosityAtGP(std::vector< LINALG::Matrix<2,1> > gradporosity_gp);

      /// Set gradient of jacobian evaluated at gauss points
      void   SetGradJAtGP(std::vector< LINALG::Matrix<1,3> > gradporosity_gp);
      void   SetGradJAtGP(std::vector< LINALG::Matrix<1,2> > gradporosity_gp);

      /// Return porosity gradient evaluated at some gauss point.  %urrecha 04/12
      void GetGradPorosityAtGP(LINALG::Matrix<3,1>& gradporosity, int gp) const;
      void GetGradPorosityAtGP(LINALG::Matrix<2,1>& gradporosity, int gp) const;

      /// Return gradient of jacobian evaluated at some gauss point
      void GetGradJAtGP(LINALG::Matrix<1,3>& gradJ, int gp) const;
      void GetGradJAtGP(LINALG::Matrix<1,2>& gradJ, int gp) const;
      */

      //! @name Evaluation methods

      /// evaluate material law
      virtual void Evaluate(const LINALG::Matrix<3,3>* defgrd,
                            const LINALG::Matrix<6,1>* glstrain,
                            Teuchos::ParameterList& params,
                            LINALG::Matrix<6,1>* stress,
                            LINALG::Matrix<6,6>* cmat)
      {mat_->Evaluate(defgrd,glstrain,params,stress,cmat);};

      /// Evaluate strain energy function (for hyperelastic materials only)
      virtual void StrainEnergy(const LINALG::Matrix<6,1>& glstrain,
                                double& psi)
      {mat_->StrainEnergy(glstrain,psi);};

      /// Evaluate material law for GEMM
      virtual void  EvaluateGEMM(LINALG::Matrix<MAT::NUM_STRESS_3D,1>* stress,
                                 LINALG::Matrix<MAT::NUM_STRESS_3D,MAT::NUM_STRESS_3D>* cmat,
                                 double* density,
                                 LINALG::Matrix<MAT::NUM_STRESS_3D,1>* glstrain_m,
                                 LINALG::Matrix<MAT::NUM_STRESS_3D,1>* glstrain_new,
                                 LINALG::Matrix<MAT::NUM_STRESS_3D,1>* glstrain_old,
                                 LINALG::Matrix<3,3>* rcg_new,
                                 LINALG::Matrix<3,3>* rcg_old)
      {mat_->EvaluateGEMM(stress,cmat,density,glstrain_m,glstrain_new,glstrain_old,rcg_new,rcg_old);};

      //@}

      /// Return material density (if provided by the specific material)
      virtual double Density() const {return mat_->Density();};

      //! @name Handling of Gauss point data. Here, the poro material just calls the underlying material

      /// Set up for materials with GP data (e.g., history variables)
      virtual void Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
      {
        //setup the underlying material
        // Note: poro material itself is setup when calling PoroSetup()
        mat_->Setup(numgp,linedef);
      };

      /// Update of GP data (e.g., history variables)
      virtual void Update(){mat_->Update();};

      /// Reset time step (for time adaptivity)
      virtual void ResetStep(){mat_->ResetStep();};

      /// Reset internal variables to state in the beginning of the computation (needed for inverse analysis)
      virtual void ResetAll(int numgp){ mat_->ResetAll(numgp);};

      //@}

      //! @name Visualization methods

      /// Return names of visualization data
      virtual void VisNames(std::map<string,int>& names);

      /// Return visualization data
      virtual bool VisData(const string& name, std::vector<double>& data, int numgp);

      //@}

      protected:

      /// compute current porosity and save it
      void ComputePorosity(  const double& initporosity,    ///< (i) initial/reference porosity at gauss point
                             const double& press,           ///< (i) pressure at gauss point
                             const double& J,               ///< (i) determinant of jacobian at gauss point
                             const int& gp,                 ///< (i) number of current gauss point
                             double& porosity,              ///< (o) porosity at gauss point
                             double* dphi_dp,               ///< (o) first derivative of porosity w.r.t. pressure at gauss point
                             double* dphi_dJ,               ///< (o) first derivative of porosity w.r.t. jacobian at gauss point
                             double* dphi_dJdp,             ///< (o) derivative of porosity w.r.t. pressure and jacobian at gauss point
                             double* dphi_dJJ,              ///< (o) second derivative of porosity w.r.t. jacobian at gauss point
                             double* dphi_dpp,              ///< (o) second derivative of porosity w.r.t. pressure at gauss point
                             double* dphi_dphiref,          ///< (o) derivative of porosity w.r.t. reference porosity (only nonzero with reaction)
                             bool save = true
                             ) ;

      /// my material parameters
      MAT::PAR::StructPoro* params_;

      /// actual material
      Teuchos::RCP<MAT::So3Material> mat_;

      /// porosity at gauss points
      Teuchos::RCP<std::vector<double> > porosity_;

      /// porosity at gauss points of surface element
      Teuchos::RCP<std::map<int, std::vector< double > > > surfporosity_;

      /// porosity time derivative at gauss points
      //Teuchos::RCP<std::vector<double> > dporodt_;

      /// porosity gradient at gauss points
     // Teuchos::RCP<std::vector<LINALG::Matrix<3,1> > > gradporosity_;

      /// gradient of jacobian at gauss points
      //Teuchos::RCP<std::vector<LINALG::Matrix<1,3> > > gradJ_;

      /// flag indicating initialization of attributes
      bool isinitialized_;

  };

}

#endif /* STRUCTPORO_H_ */
