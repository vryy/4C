/*!----------------------------------------------------------------------
\file biochemo_mechano_cell_passivefiber.H

\brief Implementation of Biochemo-Mechano Coupled passive, viscoelastic material model for the cell.

\level 3

\maintainer Andreas Rauch

*----------------------------------------------------------------------*/

#ifndef BIOCHEMO_MECHANO_CELL_PASSIVEFIBER_H_
#define BIOCHEMO_MECHANO_CELL_PASSIVEFIBER_H_

#include "so3_material.H"
#include "matpar_parameter.H"
#include "../drt_lib/drt_parobjectfactory.H"


namespace MAT
{

  namespace PAR
  {

  /*! \brief Material parameters for biochemo-mechano cell cytoskeleton model.
   *
   * \date 01/2017
   * \author rauch (rauch@lnm.mw.tum.de)
  */
    class BioChemoMechanoCellPassiveFiber : public Parameter
    {
    public:

      /// standard constructor
      BioChemoMechanoCellPassiveFiber(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~BioChemoMechanoCellPassiveFiber() { ; }

      /// @name material parameters
      //@{
      /// elastic material number
      const int idmatelast_;
      /// viscosity of isotropic viscous part 2*mu_*I
      const double mu_;
      /// build analytical constitutive matrix
      bool analyticalmaterialtangent_;
      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class BioChemoMechanoCellPassiveFiber
  }  // namespace PAR


  /*! \brief Type BioChemoMechanoCellPassiveFiberType.
   *
   * \date 01/2017
   * \author rauch (rauch@lnm.mw.tum.de)
  */
  class BioChemoMechanoCellPassiveFiberType : public DRT::ParObjectType
  {
  public:

    /// return material name
    std::string Name() const { return "BioChemoMechanoCellPassiveFiberType"; }

    /// return the instance
    static BioChemoMechanoCellPassiveFiberType & Instance() { return instance_; };

    /// create ParObject
    virtual DRT::ParObject* Create( const std::vector<char> & data );

  private:

    /// static instance of passive fiber type
    static BioChemoMechanoCellPassiveFiberType instance_;
  }; //class PassiveFiberType


  /*! \brief Implementation of Biochemo-Mechano Coupled Cytoskeleton Formation.
   *
   *
   *   Exemplary Input Line:
   *
   *    MAT 1 MAT_BIOCHEMOMECHANO_PASSIVE VISC 10.0
   *
   *
   * \date 01/2017
   * \author rauch (rauch@lnm.mw.tum.de)
  */
  class BioChemoMechanoCellPassiveFiber : public So3Material
  {
  public:

    /// construct empty material object
    BioChemoMechanoCellPassiveFiber();

    /// construct the material object given material parameters
    explicit BioChemoMechanoCellPassiveFiber(MAT::PAR::BioChemoMechanoCellPassiveFiber* params);

    //! @name Packing and Unpacking
    //@{
    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return BioChemoMechanoCellPassiveFiberType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_biochemomechano_passive; }

    /// check if element kinematics and material kinematics are compatible
    virtual void ValidKinematics(INPAR::STR::KinemType kinem)
    { if (!(kinem == INPAR::STR::kinem_nonlinearTotLag))
    dserror("element and material kinematics are not compatible");}

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new BioChemoMechanoCellPassiveFiber( *this ) ); }

    /// check if history variables are already initialized
    bool Initialized() const
    {
      return ( isinit_ );
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
    virtual void Evaluate(const LINALG::Matrix<3,3>* defgrd,   ///< deformation gradient
                          const LINALG::Matrix<6,1>* glstrain, ///< green lagrange strain
                          Teuchos::ParameterList& params,      ///< parameter list
                          LINALG::Matrix<6,1>* stress,         ///< 2nd PK-stress
                          LINALG::Matrix<6,6>* cmat,           ///< material stiffness matrix
                          const int eleGID                     ///< global id of element
                          );

    /// Return Viscosity
    virtual double Viscosity() const {return params_->mu_;}

    /// Return Density
    virtual double Density() const {return matelast_->Density();};

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    /// Return names of visualization data
    virtual void VisNames(std::map<std::string,int>& names);

    /// Return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID);

    /// Return parameter of this material
    MAT::PAR::BioChemoMechanoCellPassiveFiber* GetMaterialParams(){return params_;};

  private:

    /// (current) deformation gradient
     Teuchos::RCP<std::vector<LINALG::Matrix<3,3> > > histdefgrdcurr_;
    /// deformation gradient of old timestep (last converged state)
    Teuchos::RCP<std::vector<LINALG::Matrix<3,3> > > histdefgrdlast_;

    /// elastic material
    Teuchos::RCP<MAT::So3Material> matelast_;
    /// my material parameters
    MAT::PAR::BioChemoMechanoCellPassiveFiber* params_;
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

    /// Convert spatial stress to material stress
    void CauchytoPK2(
        LINALG::Matrix<6,1>& Sactive, ///< active material stress
        LINALG::Matrix<3,3>& cauchystress, ///< cauchy stress in matrix notation
        LINALG::Matrix<3,3> defgrd,   ///< deformation gradient
        LINALG::Matrix<3,3> invdefgrd,///< inverse deformation gradient
        LINALG::Matrix<6,1> sigma     ///< active spatial stress
     );

  }; // class BioChemoMechanoCellPassiveFiber
  } // namespace MAT


#endif /*BIOCHEMO_MECHANO_CELL_PASSIVEFIBER_H_*/
