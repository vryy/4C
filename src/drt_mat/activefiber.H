/*!----------------------------------------------------------------------
\file activefiber.H

\brief active fiber material for cell contractibility

\level 3

 This file contains routines for a local material law that contains active
 fiber formation and orientation for the modeling of living cells

 example input line
 MAT 1 MAT_ACTIVEFIBER DENS 1.0 DECAY 720.0 IDMATPASSIVE 2 KFOR 10.0 KBACK 1.0 KVAR 10.0 SIGMAX 3.9E-03 EPSNULL 2.8E-04

 For a detailed description of the model see:

 - Deshpande, V., McMeeking, R.M., Evans, A.G., 2007. A model for the
   contractility of the cytoskeleton including the effects of stress-fibre
   formation and dissociation, Proceedings of the Royal Society A:
   Mathematical, Physical and Engineering Sciences 463, 787-815.

   ____________________________________________________________________________________
  | !!! ATTENTION !!! Many major mistakes in literature that were corrected here !!!!  |
  |____________________________________________________________________________________|

<pre>
\maintainer Andreas Rauch
            rauch@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*----------------------------------------------------------------------*/

#ifndef ACTIVEFIBER_H_
#define ACTIVEFIBER_H_


#include "so3_material.H"
#include "matpar_parameter.H"
#include "../drt_lib/drt_parobjectfactory.H"

# define numbgp 10
# define twice 20


namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters
    class ActiveFiber : public Parameter
    {
    public:

      /// standard constructor
      ActiveFiber(
        Teuchos::RCP<MAT::PAR::Material> matdata
        );

      /// destructor
      virtual ~ActiveFiber() { ; }

      /// @name material parameters
      //@{
      /// density
      const double density_;
      /// decay constant of activation signal
      const double decayconst_;
      /// elastic material number
      const int idmatpassive_;
      /// non-dimensional parameter governing the rate of formation of stress fibers
      const double kforwards_;
      /// non-dimensional parameter governing the rate of dissociation of stress fibers
      const double kbackwards_;
      /// non-dimensional fiber rate sensitivity
      const double kvariance_;
      /// maximum tension exerted by stress fibers
      const double sigmamax_;
      /// reference strain rate of cross-bridge dynamics law
      const double epsilonnull_;
      /// build analytical constitutive matrix
      bool analyticalmaterialtangent_;
      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class ActiveFiber

  }  // namespace PAR

  class ActiveFiberType : public DRT::ParObjectType
  {
  public:

    std::string Name() const { return "ActiveFiberType"; }

    static ActiveFiberType & Instance() { return instance_; };

    virtual DRT::ParObject* Create( const std::vector<char> & data );

  private:

    static ActiveFiberType instance_;
  }; //class ActiveFiberType

  /*----------------------------------------------------------------------*/
  /// Wrapper for active fiber material
  class ActiveFiber : public So3Material
  {
  public:

    /// construct empty material object
    ActiveFiber();

    /// construct the material object given material parameters
    explicit ActiveFiber(MAT::PAR::ActiveFiber* params);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ActiveFiberType::Instance().UniqueParObjectId(); }

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
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_activefiber; }

    /// check if element kinematics and material kinematics are compatible
    virtual void ValidKinematics(INPAR::STR::KinemType kinem)
    { if (!(kinem == INPAR::STR::kinem_nonlinearTotLag))
    dserror("element and material kinematics are not compatible");}

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new ActiveFiber( *this ) ); }

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
    MAT::PAR::ActiveFiber* GetMaterialParams(){return params_;};

  private:

    /// my material parameters
    MAT::PAR::ActiveFiber* params_;

    /// (current) fiber activation level
    Teuchos::RCP<std::vector<LINALG::Matrix<numbgp,twice> > > etacurr_;
    /// fiber activation level of old timestep (last converged state)
    Teuchos::RCP<std::vector<LINALG::Matrix<numbgp,twice> > > etalast_;
    /// average intensity level at every point of the cytoplasm
    Teuchos::RCP<std::vector<double> > etahat_;
    /// average intensity levels at three different directions
    Teuchos::RCP<std::vector<double> > etahor_;
    Teuchos::RCP<std::vector<double> > etaver_;
    Teuchos::RCP<std::vector<double> > etadiag_;
//    Teuchos::RCP<std::vector<double> > dxx_;
//    Teuchos::RCP<std::vector<double> > dyy_;
//    Teuchos::RCP<std::vector<double> > dzz_;
//    Teuchos::RCP<std::vector<double> > dxy_;
//    Teuchos::RCP<std::vector<double> > dyz_;
//    Teuchos::RCP<std::vector<double> > dxz_;
    /// (current) stress in (omega,phi) direction
    Teuchos::RCP<std::vector<LINALG::Matrix<numbgp,twice> > > sigmaomegaphicurr_;
    /// stress in (omega,phi) direction of old timestep (last converged state)
    Teuchos::RCP<std::vector<LINALG::Matrix<numbgp,twice> > > sigmaomegaphilast_;
    /// (current) deformation gradient
    Teuchos::RCP<std::vector<LINALG::Matrix<3,3> > > histdefgrdcurr_; ///< active history deformation gradient
    /// deformation gradient of old timestep (last converged state)
    Teuchos::RCP<std::vector<LINALG::Matrix<3,3> > > histdefgrdlast_; ///< active history of deformation gradient
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

    /// Calculate activation signal of current and last timestep
    void CalcActivationSignal(
        double * Csignal, ///< activation signal at current timestep
        Teuchos::ParameterList& params,
        double * Csignalold ///< activation signal at last timestep
    );

    /// Evaluate fiber formation law
    void EvaluateFormationLaw(
        Teuchos::ParameterList& params,
        double * etacurr, ///< fiber activation level in (omega,phi) direction of iteration step i+1
        double etalast, ///< fiber activation level in (omega,phi) direction of old time step (last converged step)
        double Csignal, ///< value of activation signal of new time step
        double Csignalold, ///< value of activation signal of old time step (last converged step)
        double sigmaomegaphicurr, ///< stress in (omega,phi) direction of iteration step i (i.e. from last iteration step)
        double sigmaomegaphilast ///< stress in (omega,phi) direction of old time step (last converged state)
    );

    /// Evaluate Hill stress law
    void EvaluateHillStressLaw(
        double * sigmaomegaphinew, ///< value of stress in (omega,phi) direction of iteration step i+1 (i.e. newest iteration value)
        double etacurr, ///< fiber activation level of (omega,phi) direction of iteration step i+1
        double epsomegaphi ///< strain rate in (omega,phi) direction of iteration step i (i.e. from last iteration step)
    );

    /// Numerical volume averaging of active stress by considering a finite number of stress fibers
    void NumCalcInteg(
        double sigmaomegaphinew,     ///< value of stress in (omega,phi) direction of iteration step i+1 (i.e. newest iteration value)
        LINALG::Matrix<6,1>& sigma  ///< sigma stress in cartesian coordinates in iteration step i+1
    );

    /// Convert spatial stress to material stress
    void CauchytoPK2(
        LINALG::Matrix<6,1>& Sactive, ///< active material stress
        LINALG::Matrix<3,3>& cauchystress, ///< cauchy stress in matrix notation
        LINALG::Matrix<3,3> defgrd,   ///< deformation gradient
        LINALG::Matrix<3,3> invdefgrd,///< inverse deformation gradient
        LINALG::Matrix<6,1> sigma     ///< active spatial stress
     );

    /// Convert material strain to spatial strain
    void GLtoEA(
        LINALG::Matrix<6,1> glstrain,
        LINALG::Matrix<3,3> invdefgrd,
        LINALG::Matrix<3,3>& eastrain
     );


    /// Setup 4-Tensor from 6x6 Voigt notation
    void Setup4Tensor(
        double FourTensor[3][3][3][3],
        LINALG::Matrix<6,6> VoigtMatrix
    );

    /// Setup 6x6 Voigt matrix from 4-tensor with minor symmetries
    void Setup6x6VoigtMatrix(
        LINALG::Matrix<6,6>& VoigtMatrix,
        double FourTensor[3][3][3][3]
    );

    /// Multiply 4-Tensor*Matrix
    void MultFourTensorMatrix(
        double FourTensorResult[3][3][3][3],
        double FourTensor[3][3][3][3],
        LINALG::Matrix<3,3> Matrix,
        bool clearresulttensor = true
    );

    /// Multiply Matrix*4-Tensor
    void MultMatrixFourTensor(
        double FourTensorResult[3][3][3][3],
        LINALG::Matrix<3,3> Matrix,
        double FourTensor[3][3][3][3],
        bool clearresulttensor = true
    );

    /// Set every tensor value to zero
    void Clear4Tensor(double FourTensor[3][3][3][3]);

    /// Transformation of strain rate tensor in (omega,phi) direction
    void TransformStrainRate(
        LINALG::Matrix<6,1> strainrate, ///< strain rate tensor epsilonpunkt_ij at every point in cartesian coordinates
        double * epsomegaphi ///< strain of fibers in (omega,phi) direction
    );

    /// computes elasticity tensor for the active stress part in matrix notion for 3d
    void SetupCmatActive(
        LINALG::Matrix<6,6>& cmatactive,
        LINALG::Matrix<3,3> rotationrate,
        LINALG::Matrix<6,1> strainrate,
        LINALG::Matrix<3,3> defgrd,
        LINALG::Matrix<3,3> defgrdrate,
        LINALG::Matrix<3,3> R,
        LINALG::Matrix<3,3> invdefgrd,
        LINALG::Matrix<numbgp,twice> etanew,
        LINALG::Matrix<numbgp,twice> sigmaomegaphicurr,
        LINALG::Matrix<3,3> cauchystress,
        Teuchos::ParameterList& params,
        double theta,
        double Csignal
    );

    /// Calculate root of a 3x3 matrix (symmetric or non-symmetric)
    void MatrixRoot3x3(LINALG::Matrix<3,3>& MatrixInOut);

    /// Calculate derivative of the root of a symmetric 3x3 matrix
    void MatrixRootDerivativeSym3x3(
        const LINALG::Matrix<3,3> MatrixIn,
        LINALG::Matrix<6,6>& MatrixRootDeriv);

    /// Transpose four tensor wrt to basis vectors 1 and 2
    void TransposeFourTensor12(
                                double ResultTensor[3][3][3][3],
                                double FourTensor[3][3][3][3]
                               );

    /// print four tensor
    void PrintFourTensor(double FourTensor[3][3][3][3]);

  };  //class ActiveFiber

  /// Debug output to gmsh-file
  /* this needs to be copied to STR::TimInt::OutputStep() to enable debug output
  {
    discret_->SetState("displacement",Dis());
    MAT::ActiveFiberOutputToGmsh(discret_, GetStep(), 1);
  }
  don't forget to include activefiber.H */
  void ActiveFiberOutputToGmsh
  (
    const Teuchos::RCP<DRT::Discretization> dis, ///< discretization with displacements
    const int timestep, ///< index of timestep
    const int iter ///< iteration index of newton iteration
  );

}  // namespace MAT


#endif /*ACTIVEFIBER_*/


