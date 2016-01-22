/*!----------------------------------------------------------------------
\file biochemo_mechano_cell_activefiber.H

\brief Implementation of Biochemo-Mechano Coupled Stress Fiber Formation and Dissociation.
       Micro-Macro Approach (Volume Averaging) for Angular Stress Fiber Distribution.
       Active Stress is Additive to Passive Stress Modeled by Arbitrary Law.
       Interacts with Transport and Reaction of Scalar Quantities.

 Exemplary Input Line:

 MAT 1 MAT_CELLCONTRACTION DENS 1.0 IDMATPASSIVE 2 KFOR 10.0 KBACK 10.0 NMAX 1.0 KSTRESS 1.0

<pre>
Maintainer: Andreas Rauch
            rauch@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*----------------------------------------------------------------------*/

#ifndef BIOCHEMO_MECHANO_CELL_ACTIVEFIBER_H_
#define BIOCHEMO_MECHANO_CELL_ACTIVEFIBER_H_

#include "so3_material.H"
#include "matpar_parameter.H"
#include "../drt_lib/drt_parobjectfactory.H"

#define option3
#define nphi 20 // total (!) number of integration points, i.e. in 3D we have nphi/ntheta integration points in phi direction
#define ntheta 4 // in case of 3D integration: number of integration points in theta direction


namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters
    class BioChemoMechanoCellActiveFiber : public Parameter
    {
    public:

      /// standard constructor
      BioChemoMechanoCellActiveFiber(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~BioChemoMechanoCellActiveFiber() { ; }

      /// @name material parameters
      //@{
      /// density
      const double density_;
      /// decay constant of activation signal
//      const double decayconst_;
//      /// elastic material number
      const int idmatpassive_;
      /// non-dimensional parameter governing the rate of formation of stress fibers
      const double kforwards_;
      /// non-dimensional parameter governing the rate of dissociation of stress fibers
      const double kbackwards_;
      const double kRockEta_;
      const double kActin_;
      //const double DissCoupling_;
      const double ratemax_;
      const double nmax_;
      const double kstress_;
      const double sourceconst_;
      const std::string myintmethod_;

//      /// build analytical constitutive matrix
      bool analyticalmaterialtangent_;
      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class ActiveFiber

  }  // namespace PAR

  class BioChemoMechanoCellActiveFiberType : public DRT::ParObjectType
  {
  public:

    std::string Name() const { return "BioChemoMechanoCellActiveFiberType"; }

    static BioChemoMechanoCellActiveFiberType & Instance() { return instance_; };

    virtual DRT::ParObject* Create( const std::vector<char> & data );

  private:

    static BioChemoMechanoCellActiveFiberType instance_;
  }; //class ActiveFiberType

  /*----------------------------------------------------------------------*/
  /// Wrapper for active fiber material
  class BioChemoMechanoCellActiveFiber : public So3Material
  {
  public:

    /// construct empty material object
    BioChemoMechanoCellActiveFiber();

    /// construct the material object given material parameters
    explicit BioChemoMechanoCellActiveFiber(MAT::PAR::BioChemoMechanoCellActiveFiber* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return BioChemoMechanoCellActiveFiberType::Instance().UniqueParObjectId(); }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.
      This material contains history variables, which are packed for restart purposes.

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
      History data is unpacked in restart.

      \param data (in) : vector storing all data to be unpacked into this
      instance.
    */
    virtual void Unpack(const std::vector<char>& data);

    //@}

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_biochemomechano; }

    /// check if element kinematics and material kinematics are compatible
    virtual void ValidKinematics(INPAR::STR::KinemType kinem)
    { if (!(kinem == INPAR::STR::kinem_nonlinearTotLag))
    dserror("element and material kinematics are not compatible");}

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new BioChemoMechanoCellActiveFiber( *this ) ); }

    //! check if history variables are already initialized
    bool Initialized() const
    {
      return ( isinit_ and (histdefgrdcurr_ != Teuchos::null) );
    }

    /// Setup
    virtual void Setup(int numgp, DRT::INPUT::LineDefinition* linedef);

    /// Reset internal variables to state in the beginning of the computation (needed for inverse analysis)
    virtual void ResetAll(int numgp);

    /// Update
    virtual void Update();

    /// Reset time step
    virtual void ResetStep();

    /// Evaluate material
    virtual void Evaluate(const LINALG::Matrix<3,3>* defgrd,
                          const LINALG::Matrix<6,1>* glstrain,  ///<green lagrange strain
                          Teuchos::ParameterList& params,
                          LINALG::Matrix<6,1>* stress,  ///< 2nd PK-stress
                          LINALG::Matrix<6,6>* cmat,  ///< material stiffness matrix
                          const int eleGID);

    /// Return density
    virtual double Density() const {return params_->density_;}

    /// Return elastic material
    //Teuchos::RCP<MAT::Material> Matpassive() const { return matpassive_; }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    /// Return names of visualization data
    virtual void VisNames(std::map<std::string,int>& names);

    /// Return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID);

    /// Return parameter of this material
    MAT::PAR::BioChemoMechanoCellActiveFiber* GetMaterialParams(){return params_;};

  private:

    /// my material parameters
    MAT::PAR::BioChemoMechanoCellActiveFiber* params_;

    /// (current) fiber activation level
    Teuchos::RCP<std::vector<LINALG::Matrix<3,3> > > etacurr_;
    /// fiber activation level of old timestep (last converged state)
    Teuchos::RCP<std::vector<LINALG::Matrix<3,3> > > etalast_;
//#ifdef option3
    Teuchos::RCP<std::vector<LINALG::Matrix<nphi,1> > > etalastphi_;
    Teuchos::RCP<std::vector<LINALG::Matrix<nphi,1> > > etalastphisurf_;
    Teuchos::RCP<std::vector<LINALG::Matrix<nphi,1> > > etacurrphi_;
    Teuchos::RCP<std::vector<LINALG::Matrix<nphi,1> > > etacurrphisurf_;
//#endif
    /// average intensity level at every point of the cytoplasm
    Teuchos::RCP<std::vector<double> > Nfil_;
    Teuchos::RCP<std::vector<double> > Nfillast_;
    Teuchos::RCP<std::vector<double> > rate_;
    Teuchos::RCP<std::vector<double> > D_;
    Teuchos::RCP<std::vector<double> > ccurr_;
    Teuchos::RCP<std::vector<double> > clast_;
    Teuchos::RCP<std::vector<double> > nsfhor_;
    Teuchos::RCP<std::vector<double> > nsfver_;
    Teuchos::RCP<std::vector<double> > nsfdiagup_;
    Teuchos::RCP<std::vector<double> > nsfdiagdown_;

    /// (current) deformation gradient
    Teuchos::RCP<std::vector<LINALG::Matrix<3,3> > > histdefgrdcurr_; ///< active history deformation gradient
    Teuchos::RCP<std::vector<LINALG::Matrix<6,1> > >strainratecurr_;
    Teuchos::RCP<std::vector<LINALG::Matrix<6,1> > >strainratecurrsurf_;
    Teuchos::RCP<std::vector<LINALG::Matrix<3,3> > > histdefgrdcurrsurf_; ///< active history deformation gradient
    /// deformation gradient of old timestep (last converged state)
    Teuchos::RCP<std::vector<LINALG::Matrix<3,3> > > histdefgrdlast_; ///< active history of deformation gradient
    // active history of deformation gradient at surface for scatra boundary calculation
    Teuchos::RCP<std::vector<LINALG::Matrix<3,3> > > histdefgrdlastsurf_;
    Teuchos::RCP<std::vector<LINALG::Matrix<6,1> > >strainratelast_;
    Teuchos::RCP<std::vector<LINALG::Matrix<6,1> > >strainratelastsurf_;
    /// passive material
    Teuchos::RCP<MAT::So3Material> matpassive_;
    /// indicates if material is initialized
    bool isinit_;


    /// Setup defgrd rate, rotation tensor, strain rate and rotation rate
    void SetupRates(
        LINALG::Matrix<3,3> defgrd,
        LINALG::Matrix<3,3> invdefgrd,
        Teuchos::ParameterList& params,
        LINALG::Matrix<3,3>& defgrdrate,
        LINALG::Matrix<3,3>& R,
        LINALG::Matrix<6,1>& strainrate,
        LINALG::Matrix<3,3>& rotationrate,
        const int& gp,
        const double& dt
    );



    void Dissociation(
      double rate,
      double& D,
      double k_B,
      double ratemax
      );

    void GPtoNodes(
        LINALG::Matrix<8,1> invalue,
        LINALG::Matrix<8,1>& outvalue
        );

    void EtaGPtoNodes(
        LINALG::Matrix<nphi,1>& node1,
        LINALG::Matrix<nphi,1>& node2,
        LINALG::Matrix<nphi,1>& node3,
        LINALG::Matrix<nphi,1>& node4,
        LINALG::Matrix<nphi,1>& node5,
        LINALG::Matrix<nphi,1>& node6,
        LINALG::Matrix<nphi,1>& node7,
        LINALG::Matrix<nphi,1>& node8
        );

    void EtaAtBrdyGP(
            LINALG::Matrix<nphi,1>& outvalue,
            Teuchos::ParameterList& params
            );
    void cROCKlastAtBrdyGP(
            double& cROCKused,
             Teuchos::ParameterList& params
            );
    void cROCKusedAtBrdyGP(
            double& cROCKused,
             Teuchos::ParameterList& params
            );

    void DefGradAtBrdyGP(
            LINALG::Matrix<6,1>& outvalue,
            Teuchos::ParameterList& params
            );

    void DefGradGPtoNodes(
            LINALG::Matrix<6,1>& node1,
            LINALG::Matrix<6,1>& node2,
            LINALG::Matrix<6,1>& node3,
            LINALG::Matrix<6,1>& node4,
            LINALG::Matrix<6,1>& node5,
            LINALG::Matrix<6,1>& node6,
            LINALG::Matrix<6,1>& node7,
            LINALG::Matrix<6,1>& node8
            );

    /// Convert spatial stress to material stress
    void CauchytoPK2(
        LINALG::Matrix<6,1>& Sactive, ///< active material stress
        LINALG::Matrix<3,3>& cauchystress, ///< cauchy stress in matrix notation
        LINALG::Matrix<3,3> defgrd,   ///< deformation gradient
        LINALG::Matrix<3,3> invdefgrd,///< inverse deformation gradient
        LINALG::Matrix<6,1> sigma     ///< active spatial stress
     );

    /// Convert material stress to spatial stress
    void PK2toCauchy(
        LINALG::Matrix<6,1> Spassive, ///< passive material stress
        LINALG::Matrix<3,3>& PK2stress, ///< PK2 stress in matrix notation
        LINALG::Matrix<3,3> defgrd,   ///< deformation gradient
        LINALG::Matrix<6,1>& cauchystess     ///< passive spatial stress
     );

    /// Convert material strain to spatial strain
    void GLtoEA(
        LINALG::Matrix<6,1> glstrain,
        LINALG::Matrix<3,3> invdefgrd,
        LINALG::Matrix<3,3>& eastrain
     );

    /// Transformation of strain rate tensor in (omega,phi) direction
    void TransformStrainRate(
        LINALG::Matrix<6,1> strainrate, ///< strain rate tensor epsilonpunkt_ij at every point in cartesian coordinates
        double * epsomegaphi ///< strain of fibers in (omega,phi) direction
    );


  };  //class BioChemoMechanoCellActiveFiber
  }  // namespace MAT


#endif /*BIOCHEMO_MECHANO_CELL_ACTIVEFIBER_H_*/
