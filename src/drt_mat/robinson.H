/*----------------------------------------------------------------------*/
/*!
\file robinson.H
\brief Robinson's visco-plastic material

      example input line:
      MAT 1 MAT_Struct_Robinson  KIND Arya_NarloyZ  YOUNG POLY 2 1.47e9 -7.05e5
        NUE 0.34  DENS 8.89e-3  THEXPANS 0.0  HRDN_FACT 3.847e-12  HRDN_EXPO 4.0
        SHRTHRSHLD POLY 2 69.88e8 -0.067e8   RCVRY 6.083e-3  ACTV_ERGY 40000.0
        ACTV_TMPR 811.0  G0 0.04  M_EXPO 4.365  BETA POLY 3 0.8 0.0 0.533e-6
        H_FACT 1.67e16
      12.01.12 as first step only implement geometric non-linear case, i.e., total Lagrange


<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef VISCOPLASTICROBINSON_H_
#define VISCOPLASTICROBINSON_H_

#ifdef CCADISCRET

#include "matpar_parameter.H"
#include "material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ParameterList.hpp>

using namespace Teuchos;

namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// material parameters for visco-plastic Robinson's material
    class Robinson : public Parameter
    {
    public:

      /// standard constructor
      Robinson(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~Robinson() { ; }

      /// @name material parameters
      //@{

      /// kind of Robinson material (slight differences:vague,butler,arya,arya_narloyz,arya_crmosteel)
      const string* kind_;
      /// number of young's modulus data
      const int youngsnum_;
      /// Young's modulus (temperature dependent --> polynomial expression)
      const vector<double> youngs_;
      /// Possion's ratio
      const double poissonratio_;
      /// mass density
      const double density_;
      /// linear coefficient of thermal expansion
      const double thermexpans_;
      /// hardening factor 'A' (needed for flow law)
      const double hrdn_fact_;
      /// hardening power 'n'  (exponential of F in the flow law)
      const double hrdn_expo_;
      /// number of 'K^2' data
      const int shrthrshldnum_;
      /// Bingam-Prager shear stress threshold \kappa^2 'K^2=K^2(K_0)'
      const vector<double> shrthrshld_;
      /// recovery factor 'R_0'
      const double rcvry_;
      /// activation energy 'Q_0'
      const double actv_ergy_;
      /// activation temperature 'T_0'
      const double actv_tmpr_;
      /// 'G_0' (temperature independent, minimum value attainable by G )
      const double g0_;
      /// 'm' (temperature independent, exponential in evolution law of back stress)
      const double m_;
      /// number of 'beta' data
      const int betanum_;
      /// 'beta' (temperature independent)
      /// Arya_NarloyZ: beta = 0.533e-6 T^2 + 0.8
      const vector<double> beta_;
      /// h
      /// Arya_NarloyZ: H = 1.67e4 . (6.895)^(beta - 1) / (3 . K_0^2)
      const double h_;

      //@}

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class Robinson

  }  // namespace PAR


class RobinsonType : public DRT::ParObjectType
{
public:

  std::string Name() const { return "RobinsonType"; }

  static RobinsonType & Instance() { return instance_; };

  virtual DRT::ParObject* Create( const vector<char> & data );

private:

  static RobinsonType instance_;

};  // RobinsonType


  /*----------------------------------------------------------------------*/
  /// wrapper for visco-plastic Robinson's material
  class Robinson : public Material
  {
  public:

    /// construct empty material object
    Robinson();

    /// construct the material object given material parameters
    explicit Robinson(MAT::PAR::Robinson* params);

    //! @name Packing and Unpacking

    //!  \brief return unique ParObject id
    //!
    //!  every class implementing ParObject needs a unique id defined at the
    //!  top of drt_parobject.H (this file) and should return it in this method.
    virtual int UniqueParObjectId() const { return RobinsonType::Instance().UniqueParObjectId(); }

    //!  \brief Pack this class so it can be communicated
    //!
    //!  Resizes the vector data and stores all information of a class in it.
    //!  The first information to be stored in data has to be the
    //!  unique parobject id delivered by UniqueParObjectId() which will then
    //!  identify the exact class on the receiving processor.
    //!
    //!  \param data (in/out): char vector to store class information
    virtual void Pack(DRT::PackBuffer& data) const;

    //!  \brief Unpack data from a char vector into this class
    //!
    //!  The vector data contains all information to rebuild the
    //!  exact copy of an instance of a class on a different processor.
    //!  The first entry in data has to be an integer which is the unique
    //!  parobject id defined at the top of this file and delivered by
    //!  UniqueParObjectId().
    //!
    //!  \param data (in) : vector storing all data to be unpacked into this
    //!  instance.
    virtual void Unpack(const vector<char>& data);

    //@}

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const { return INPAR::MAT::m_vp_robinson; }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const { return Teuchos::rcp( new Robinson( *this ) ); }

    /// initialize internal stress variables
    void Setup(
      const int numgp, ///< number of Gauss points
//      KinematicType kintype
      DRT::INPUT::LineDefinition* linedef
      );

    /// update internal stress variables
    void Update();

    /// reset internal stress variables
    void Reset();

    /// evaluate material
    void Evaluate(
      const LINALG::Matrix<NUM_STRESS_3D,1>& strain,
      LINALG::Matrix<NUM_STRESS_3D,1>& plstrain, ///< return the plastic strains to the element
      const int gp,  ///< current Gaus point
      Teuchos::ParameterList& params,  ///< parameter list for communication
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& cmat,  ///< material stiffness matrix
      LINALG::Matrix<NUM_STRESS_3D,1>& stress  ///< 2nd PK-stress
      );

    // computes stress
    void Stress(
      const double p,  //!< volumetric stress tensor
      const LINALG::Matrix<NUM_STRESS_3D,1>& devstress,  //!< deviatoric stress tensor
      LINALG::Matrix<NUM_STRESS_3D,1>& stress //!< 2nd PK-stress
      );

    // computes relative stress eta = stress - back stress
    void RelDevStress(
      const LINALG::Matrix<NUM_STRESS_3D,1>& devstrain,  //!< deviatoric stress tensor
      const LINALG::Matrix<NUM_STRESS_3D,1>& beta,  //!< back stress tensor
      LINALG::Matrix<NUM_STRESS_3D,1>& eta //!< relative stress
      );

    /// computes isotropic elasticity tensor in matrix notion for 3d
    void SetupCmat(
      double temp,  //!< current temperature
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& cmat  //!< material tangent
      );

    void CalcBEViscousStrainRate(
    //  \param   ele          ELEMENT*      (i)   element
    //  \param   mat_robin    VP_ROBINSON*  (i)   element Robinson material
      const double dt,  // (i) time step size
      double tempnp, //         DOUBLE        (i)   temperature
      const LINALG::Matrix<NUM_STRESS_3D,1>&  strain_p,  //      DOUBLE[]      (i)   viscous strain at t_n
      const LINALG::Matrix<NUM_STRESS_3D,1>&  strain_pn,  //      DOUBLE[]      (i)   viscous strain at t_{n+1}^<i>
      const LINALG::Matrix<NUM_STRESS_3D,1>&  devstress,  //    DOUBLE[]      (i)   stress deviator at t_{n+1}^<i>
      const LINALG::Matrix<NUM_STRESS_3D,1>&  eta,  //    DOUBLE[]      (i)   over stress at t_{n+1}^<i>
      LINALG::Matrix<NUM_STRESS_3D,1>& strainrate_p,  //    DOUBLE[]      (o)   viscous strain residual
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> kve,  // (o)
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> kvv,  // (o)
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> kva  // (o)
      );

    // residual of BE-discretised back stress rate flow rule at Gauss point
    void CalcBEBackStressFlow(
      const double dt,
      const double tempnp,
      const LINALG::Matrix<NUM_STRESS_3D,1>&  strain_p,  // double vscstn[NUMSTR_SOLID3],
      const LINALG::Matrix<NUM_STRESS_3D,1>&  strain_pn,  // double vscstnn[NUMSTR_SOLID3],
      const LINALG::Matrix<NUM_STRESS_3D,1>&  devstress,  // double devstsn[NUMSTR_SOLID3],
      const LINALG::Matrix<NUM_STRESS_3D,1>& backstresslast, //  const double bacsts[NUMSTR_SOLID3],
      const LINALG::Matrix<NUM_STRESS_3D,1>& backstresscurr, //const double bacstsn[NUMSTR_SOLID3],
      //TODO 22.11.11INT* bckstss,
      LINALG::Matrix<NUM_STRESS_3D,1>& backstress_res, // double bckstsr[NUMSTR_SOLID3],
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& kae, // double kae[NUMSTR_SOLID3][NUMSTR_SOLID3],
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& kav, // double kav[NUMSTR_SOLID3][NUMSTR_SOLID3],
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& kaa // double kaa[NUMSTR_SOLID3][NUMSTR_SOLID3])
      );

    // Reduce (statically condense) system in eps,eps^v,al to purely eps
    /*!
    The linearised stress and internal residuals are

          [ sig   ]         [ sig    ]^i
      Lin [ res^v ]       = [ res^v  ]
          [ res^al]_{n+1}   [ res^al ]_{n+1}

                               [ kee  kev  kea ]^i  [ iinc eps   ]^i
                            +  [ kve  kvv  kva ]    [ iinc eps^v ]
                               [ kae  kav  kaa ]    [ iinc al    ]_{n+1}

                            [ sig ]
                          = [  0  ]  on every element (e)
                            [  0  ]  and at each Gauss point <gp>

    with - strain increment         iinc eps   -->  epsii
         - viscous strain increment iinc eps^v -->  strain_pn TODO???
         - back stress increment    iinc al    -->  beta
         -                          kee        -->  cmat

    Due to the fact that the internal residuals (the BE-discretised evolution
    laws of the viscous strain and the back stress) are C^{-1}-continuous
    across element boundaries. We can statically condense this system.
    The iterative increments inc eps^v and inc al are expressed in inc eps.
    We achieve

      [ iinc eps^v ]   [ kvv  kva ]^{-1} (   [ res^v  ]   [ kve ]                )
      [            ] = [          ]      ( - [        ] - [     ] . [ iinc eps ] )
      [ iinc al    ]   [ kav  kaa ]      (   [ res^al ]   [ kae ]                )

    thus

                                         [ kvv  kva ]^{-1} [ res^v  ]^i
      sig_red^i = sig^i - [ kev  kea ]^i [          ]      [        ]
                                         [ kav  kaa ]      [ res^al ]

    and
                                         [ kvv  kva ]^{-1} [ kve ]^i
      kee_red^i = kee^i - [ kev  kea ]^i [          ]      [     ]
                                         [ kav  kaa ]      [ kae ]

      ==> condensed system: kee_red^i . iinc eps = sig_red^i

      TODO 17.01.12 what happens with sigma on the rhs???

    \param  stress          (6x1)  (io)   stress vector
    \param  cmat == kee     (6x6)  (io)   material stiffness matrix, constitutive tensor
    \param  kev             (6x6)  (i)    \frac{\pd \sig}{\pd \eps^v}
    \param  kea             (6x6)  (i)    \frac{\pd \sig}{\pd \al}
    \param  strain_pres     (6x1)  (i)    viscous strain residual (CCARAT: vscstnr)
    \param  kve             (6x6)  (i)    \frac{\pd res^v}{\pd \eps}
    \param  kvv             (6x6)  (i)    \frac{\pd res^v}{\pd \eps^v}
    \param  kva             (6x6)  (i)    \frac{\pd res^v}{\pd \al}
    \param  backstress_pres (6x1)  (i)    viscous strain residual (CCARAT: bckstsr)
    \param  kae             (6x6)  (i)    \frac{\pd res^al}{\pd \eps}
    \param  kav             (6x6)  (i)    \frac{\pd res^al}{\pd \eps^v}
    \param  kaa             (6x6)  (i)    \frac{\pd res^al}{\pd \al}
    \param  kvarva          (12x1) (o)    condensed matrix of residual
    \param  kvakvae         (12x6) (o)    condensed matrix of tangent
    */
    void CalculateCondensedSystem(
      LINALG::Matrix<NUM_STRESS_3D,1>& stress,
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& cmat,
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& kev,
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& kea,
      const LINALG::Matrix<NUM_STRESS_3D,1>& strain_pres,
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& kve,
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& kvv,
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& kva,
      const LINALG::Matrix<NUM_STRESS_3D,1>& backstressrate,
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& kae,
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& kav,
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& kaa,
      LINALG::Matrix<(2*NUM_STRESS_3D),1>& kvarva,
      LINALG::Matrix<(2*NUM_STRESS_3D),NUM_STRESS_3D>& kvakvae
      );

    //! \brief iterative update of material internal variables
    //!
    //! material internal variables (viscous strain and back stress) are updated by
    //! their iterative increments.
    //! Their iterative increments are expressed in terms of the iterative increment
    //! of the total strain.
    //! Here the reduction matrices (kvarvam,kvakvae) stored at previous call of
    //! CalculateCondensedSystem() care used.
    //!
    //! CCARAT: so3_mat_robinson_be_mivupditer()
    //!
    //! \params epsii (i)
    //! \params Delta strain_p (o)
    //! \params Delta backstress (o)
    void IterativeUpdateOfInternalVariables(
      const int numgp,
      const LINALG::Matrix<NUM_STRESS_3D,1> epsii
      );

    // incremental update of material internal variables
    void IncrementalUpdateOfInternalVariables(
      const int gp  //!< total number of GPs in domain
      );

    //! \brief select Robinson's material and return stress
    //!
    //! \param   gp          int       (i)   Gauss point
    //! \param   strain      (6x1)     (i)   total strain
    //! \param   stress      (6x1)     (o)   stress
    //! \param   cmat        (6x6)     (o)   elasticity tensor
    void SelectRobinsonsMaterialBEReturnStress(
      const int gp, //!< number of Gauss points
      const LINALG::Matrix<NUM_STRESS_3D,1>& strain,  //!< total strain vector
      LINALG::Matrix<NUM_STRESS_3D,1>& stress,  //!< stress
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& cmat //!< material stiffness matrix
      );

    /// return density
    virtual double Density() const { return params_->density_; }

    /// check if history variables are already initialized
    bool Initialized() const
    {
      return ( isinit_ && (strainplcurr_!=Teuchos::null) );
    }

    /// return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    //! flag plastic step was called
    bool plastic_step;

    /// @name temperature specific methods
    //@{

    //! calculate temperature dependent material parameter
    void GetMatParameterAtTempnp(
      const vector<double>* paramvector, //  (i) given parameter is a vector
      double tempnp,  // current temperature
      double parambytempn // parameter at current temperature
      );

    //! calculate temperature dependent material parameter
    void GetMatParameterAtTempnp(
      const double paramconst, // (i) given parameter is a constant
      double tempnp,  // current temperature
      double parambytempn // parameter at current temperature
      );

    //@}

    /// @name specific methods for TSI and plastic material
    //@{

    /// additive split of strain rates needed for TSI
    //   \param stepsize      (in): stepsize
    //   \param strainlinrate (in): total strain rate ( B d')
    //   \param strainelrate_ (out): elastic strain rate ( B d^e' = B d - Dgama N)
    void StrainRateSplit(
      int gp,
      const double stepsize,
      LINALG::Matrix<NUM_STRESS_3D,1>& strainlinrate
      );

    /// return current plastic strain vector \f${\varepsilon}^p_{n+1}\f$
    LINALG::Matrix<NUM_STRESS_3D,1> PlasticStrain(int gp) const
    {
      return strainplcurr_->at(gp);
    }

    //@}


  private:

    /// my material parameters
    MAT::PAR::Robinson* params_;

    /// indicator if #Initialize routine has been called
    bool isinit_;

    //! visco-plastic strain vector Ev^<gp> at t_{n} for every Gauss point gp
    //!    Ev^<g>T = [ E_11  E_22  E_33  2*E_12  2*E_23  2*E_31 ]^<g> */
    //!< \f${\varepsilon}^p_{n}\f$
    // CCARAT: ARRAY vicstn;
    Teuchos::RCP< vector< LINALG::Matrix<NUM_STRESS_3D,1> > > strainpllast_;
    //! current visco-plastic strain vector Ev^<gp> at t_{n+1} for every Gauss point gp
    //!    Ev^<g>T = [ E_11  E_22  E_33  2*E_12  2*E_23  2*E_31 ]^<g> */
    // CCARAT: ARRAY vicstnn;
    Teuchos::RCP< vector< LINALG::Matrix<NUM_STRESS_3D,1> > > strainplcurr_;  //!< \f${\varepsilon}^p_{n+1}\f$
    //! old back stress vector Av^<gp> at t_n for every Gauss point gp
    //!    Av^<gp>T = [ A_11  A_22  A_33  A_12  A_23  A_31 ]^<gp>
    // CCARAT: ARRAY bacsts;
    Teuchos::RCP< vector< LINALG::Matrix<NUM_STRESS_3D,1> > > backstresslast_;  //!< \f${\beta}_{n}\f$
    //! current back stress vector Av^<gp> at t_{n+1} for every Gauss point gp
    //!< \f${\beta}_{n+1}\f$
    //!    Av^<gp>T = [ A_11  A_22  A_33  A_12  A_23  A_31 ]^<gp> */
    // CCARAT: ARRAY bacstsn;
    Teuchos::RCP< vector< LINALG::Matrix<NUM_STRESS_3D,1> > > backstresscurr_;
    //! update vector for MIV iterative increments
    //!          [ kvv  kva ]^{-1}   [ res^v  ]
    //! kvarva = [          ]      . [        ]
    //!          [ kav  kaa ]      . [ res^al ]
    Teuchos::RCP< vector< LINALG::Matrix<(2*NUM_STRESS_3D),1> > > kvarva_;
    //! update matrix for MIV iterative increments
    //!              [ kvv  kva ]^{-1}   [ kve ]
    //!    kvakvae = [          ]      . [     ]
    //!              [ kav  kaa ]      . [ kae ]
    Teuchos::RCP< vector< LINALG::Matrix<(2*NUM_STRESS_3D),NUM_STRESS_3D> > > kvakvae_;

    //! kinematic types of so_hex8 (default: total-lagrangian)
    // TODO 12.01.12 implement kinematic type that is passed from the structural
    // element to the material
    enum KinematicType
    {
      soh8_geolin,
      soh8_totlag
    };
    //! kinematic type
    KinematicType kintype_;

  };  // class Robinson : public Material
}  // namespace MAT


/*----------------------------------------------------------------------*/
#endif  // CCADISCRET
#endif  // ROBINSON_H_


