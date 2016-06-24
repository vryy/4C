/*----------------------------------------------------------------------*/
/*!
 \file growth_ip.H

 \brief This file contains routines for an integration point based isotropic volumetric growth law.

\level 3

 <pre>
   \maintainer Moritz Thon
               thon@lnm.mw.tum.de
               http://www.mhpc.mw.tum.de
               089 - 289-10364
 </pre>
 *----------------------------------------------------------------------*/
#ifndef GROWTH_IP_H
#define GROWTH_IP_H


#include "so3_material.H"
#include "matpar_parameter.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_lib/drt_container.H"

/*------------------------------------------------------------------------*/
/* forward declarations */

namespace MAT
{

class GrowthLaw;

  namespace PAR
  {

  /*----------------------------------------------------------------------*/
  /*! \class Growth
   *  \brief Common parameters for volumetric growth
   *
   *  \author kehl
   *  \date 6/2015
   */
  class Growth : public Parameter
  {
  public:

    /// standard constructor
    Growth(
        Teuchos::RCP<MAT::PAR::Material> matdata
    );

    /// destructor
    virtual ~Growth() { ; }

    /// @name material parameters
    //@{
    /// elastic material number
    const int idmatelastic_;
    /// growthlaw material number
    const int idgrowthlaw_;
    /// start growth after starttime
    const double starttime_;
    /// stop growth after endtime
    const double endtime_;
    /// growth law
    Teuchos::RCP<MAT::GrowthLaw> growthlaw_;
    //@}

    /// create material instance of matching type with my parameters
    virtual Teuchos::RCP<MAT::Material> CreateMaterial();

  };  // class Growth

  }  // namespace PAR


  /*----------------------------------------------------------------------*/
  /*! \class Growth
       \brief Wrapper material for volumetric growth

       Base class for isovolumetric growth materials.

       It keeps a pointer to the elastic material and to a specific
       growth law (via its Parameters) and some variables common to
       all growth materials.

       Further it basically passes through the interfaces functions
       Evaluate(...) and EvaluateNonLinMass(...).

       \author kehl
       \date 6/2015
  */
  class Growth : public So3Material
  {
  public:

    /// construct empty material object
    Growth();

    /// construct the material object given material parameters
    explicit Growth(MAT::PAR::Growth* params);

    //! @name Packing and Unpacking
    //@{

    /*!
      \brief Return unique ParObject id

      \sa DRT::ParObject
    */
    virtual int UniqueParObjectId() const = 0;

    /*!
      \brief Pack this class so it can be communicated

      \sa DRT::ParObject
      \param data (in/out): char vector to store class information
    */
    virtual void Pack(DRT::PackBuffer& data) const;

    /*!
      \brief Unpack data from a char vector into this class

      \sa DRT::ParObject
      \param data (in) : vector storing all data to be unpacked into this
      instance.
    */
    virtual void Unpack(const std::vector<char>& data);

    //@}

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_growth_volumetric; }

    /// check if element kinematics and material kinematics are compatible
    virtual void ValidKinematics(INPAR::STR::KinemType kinem)
    { if (!(kinem == INPAR::STR::kinem_nonlinearTotLag))
      dserror("element and material kinematics are not compatible");}

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const = 0;

    /// Setup
    virtual void Setup(int numgp, DRT::INPUT::LineDefinition* linedef);

    /// Reset internal variables to state in the beginning of the computation (needed for inverse analysis)
    virtual void ResetAll(int numgp);

    /// Update
    void Update();

    /// Reset time step
    void ResetStep();

    /// Store history/internal variables
    virtual void StoreHistory(int timestep);

    /// Set history/internal variables
    virtual void SetHistory(int timestep);


    //! @name Evaluation methods
    //@{

    /*! \brief Evaluate material stresses and cmat
     *
     * \param In
     * defgrd - deformation gradient
     * \param In
     * glstrain - green lagrange strain
     * \param In
     * params - a parameter list as handed in from the element
     * \param Out
     * stress - stresses from material evaluation
     * \param Out
     * cmat - material stiffness matrix
     * \param In
     * eleGID - the calling element's GID
     *
     *  \author kehl
     * \date 06/2015
     */
    virtual void Evaluate(const LINALG::Matrix<3,3>* defgrd,
                          const LINALG::Matrix<6,1>* glstrain,
                          Teuchos::ParameterList& params,
                          LINALG::Matrix<6,1>* stress,
                          LINALG::Matrix<6,6>* cmat,
                          const int eleGID ) = 0;

    /*! \brief Evaluate mass change
     *
     * \param In
     * defgrd - deformation gradient
     * \param In
     * glstrain - green lagrange strain
     * \param In
     * params - a parameter list as handed in from the element
     * \param Out
     * linmass_disp - linearization wrt to displacements
     * \param Out
     * linmass_vel - linearization wrt to velocities
     * \param In
     * eleGID - the calling element's GID
     *
     *  \author kehl
     * \date 06/2015
     */
    virtual void EvaluateNonLinMass(const LINALG::Matrix<3,3>* defgrd,
                          const LINALG::Matrix<6,1>* glstrain,
                          Teuchos::ParameterList& params,
                          LINALG::Matrix<6,1>* linmass_disp,
                          LINALG::Matrix<6,1>* linmass_vel,
                          const int eleGID) = 0;

    /*! \brief Evaluate elastic material stresses and cmat
     *
     * \param In
     * defgrd - deformation gradient
     * \param In
     * glstrain - green lagrange strain
     * \param Out
     * stress - stresses from material evaluation
     * \param Out
     * cmat - material stiffness matrix
     * \param In
     * params - a parameter list as handed in from the element
     * \param In
     * eleGID - the calling element's GID
     *
     *  \author kehl
     * \date 06/2015
     */
    void EvaluateElastic(const LINALG::Matrix<3,3>* defgrd,
                                 const LINALG::Matrix<6,1>* glstrain,
                                 LINALG::Matrix<6,1>* stress,
                                 LINALG::Matrix<6,6>* cmat,
                                 Teuchos::ParameterList& params,
                                 const int eleGID );
    //@}


    //! @name Query  methods
    //@{

    /// Return density
    double Density() const {dserror("growth material needs gauss point data for density!"); return -1.0;}

    virtual double Density(int gp) const
    {
      const double density = matelastic_->Density();
      const double theta_gp = (*theta_)[gp];
      return density * theta_gp * theta_gp * theta_gp;
    };

    /// Return whether material has a varying material density
    bool VaryingDensity() const { return true;}

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Growth* Parameter() const { return params_; }

    /// access to the elastic material
    const Teuchos::RCP<MAT::So3Material> Matelastic() const {return matelastic_;}

    // read access to thetaold_
    Teuchos::RCP<const std::vector<double> > ThetaOld() const {return thetaold_;}

    //@}

  protected:

    /// growth stretch
    Teuchos::RCP<std::vector<double> > theta_;

    /// indicates if material is initialized
    bool isinit_;

  private:

    /// my material parameters
    MAT::PAR::Growth* params_;

    /// elastic material
    Teuchos::RCP<MAT::So3Material> matelastic_;

    /// growth stretch at the time step before
    Teuchos::RCP<std::vector<double> > thetaold_;

    /// data container to store and recall history data
    DRT::Container histdata_;

  }; //class Growth


  /*! \class GrowthMandelType
       \brief ParObjectType instance

       \sa DRT::ParObjectType

    \author kehl
    \date 6/2015
  */
  class GrowthMandelType : public DRT::ParObjectType
  {
  public:

    std::string Name() const { return "GrowthMandelType"; }

    static GrowthMandelType & Instance() { return instance_; };

    virtual DRT::ParObject* Create( const std::vector<char> & data );

  private:

    static GrowthMandelType instance_;
  };


  /*----------------------------------------------------------------------*/
  /*! \class GrowthMandel
      \brief Split of growth and elastic part of the deformation

      Implementation of the multiplicative split into isovolumetric growth
      and elastic part of the deformation.

      Evaluation of the growth factor \f$\theta\f$ and its derivative
      \f$\frac{\partial \vartheta}{\partial C}\f$ is passed to a generic growth law.

      \author kehl
      \date 6/2015
  */
  class GrowthMandel : public Growth
  {
  public:

    /// construct empty material object
    GrowthMandel();

    /// construct the material object given material parameters
    explicit GrowthMandel(MAT::PAR::Growth* params);

    //! @name Packing and Unpacking
    //@{

    /*!
      \brief Return unique ParObject id

      \sa DRT::ParObject
    */
    virtual int UniqueParObjectId() const { return GrowthMandelType::Instance().UniqueParObjectId(); }

    /*!
      \brief Pack this class so it can be communicated

      \sa DRT::ParObject
      \param data (in/out): char vector to store class information
    */
    virtual void Pack(DRT::PackBuffer& data) const;

    /*!
      \brief Unpack data from a char vector into this class

      \sa DRT::ParObject
      \param data (in) : vector storing all data to be unpacked into this
      instance.
    */
    virtual void Unpack(const std::vector<char>& data);

    //@}

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new GrowthMandel( *this ) ); }

    /// Setup
    virtual void Setup(int numgp, DRT::INPUT::LineDefinition* linedef);

    /// Reset internal variables to state in the beginning of the computation (needed for inverse analysis)
    virtual void ResetAll(int numgp);

    /// Return names of visualization data
    virtual void VisNames(std::map<std::string,int>& names);

    /// Return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID);

    //! @name Evaluation methods
    //@{

    /*! \brief Evaluate material stresses and cmat
     *
     * \param In
     * defgrd - deformation gradient
     * \param In
     * glstrain - green lagrange strain
     * \param In
     * params - a parameter list as handed in from the element
     * \param Out
     * stress - stresses from material evaluation
     * \param Out
     * cmat - material stiffness matrix
     * \param In
     * eleGID - the calling element's GID
     *
     *  \author kehl
     * \date 06/2015
     */
    virtual void Evaluate(const LINALG::Matrix<3,3>* defgrd,
                          const LINALG::Matrix<6,1>* glstrain,
                          Teuchos::ParameterList& params,
                          LINALG::Matrix<6,1>* stress,
                          LINALG::Matrix<6,6>* cmat,
                          const int eleGID );

    /*! \brief Evaluate mass change
     *
     * \param In
     * defgrd - deformation gradient
     * \param In
     * glstrain - green lagrange strain
     * \param In
     * params - a parameter list as handed in from the element
     * \param Out
     * linmass_disp - linearization wrt to displacements
     * \param Out
     * linmass_vel - linearization wrt to velocities
     * \param In
     * eleGID - the calling element's GID
     *
     *  \author kehl
     * \date 06/2015
     */
    virtual void EvaluateNonLinMass(const LINALG::Matrix<3,3>* defgrd,
                          const LINALG::Matrix<6,1>* glstrain,
                          Teuchos::ParameterList& params,
                          LINALG::Matrix<6,1>* linmass_disp,
                          LINALG::Matrix<6,1>* linmass_vel,
                          const int eleGID);

    /*! \brief evaluate the volumetric growth factor
     *
     * This function just hands over the comuptation of the growth factor
     * theta to the attached growthlaw.
     *
     * \param Out
     * theta - the volumetric growth factor
     * \param Out
     * dthetadC - derivative of theta wrt the chauchy green tensor
     * \param In
     * matelastic - an elastic material
     * \param In
     * defgrad - the deformation gradient
     * \param In
     * glstrain - the green lagrange strains
     * \param In
     * params - a parameter list as handed in from the element
     * \param In
     * eleGID - the calling element's GID
     *
     *  \author kehl
     * \date 06/2015
     */
    virtual void EvaluateGrowth(double* theta,
                                LINALG::Matrix<6,1>* dthetadC,
                                const LINALG::Matrix<3,3>* defgrd,
                                const LINALG::Matrix<6,1>* glstrain,
                                Teuchos::ParameterList& params,
                                const int eleGID );

    //@}

    //! @name Query  methods
    //@{

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Growth* Parameter() const { return paramsMandel_; }

    //@}

  protected:

    /// mandel stress
    Teuchos::RCP<std::vector<double> > mandel_;

  private:

    /// my material parameters
    MAT::PAR::Growth* paramsMandel_;

  }; //class GrowthMandel

}  // namespace MAT

#endif /*GROWTH_IP_H_*/
