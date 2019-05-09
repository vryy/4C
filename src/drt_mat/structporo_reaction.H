/*----------------------------------------------------------------------*/
/*!
 \brief wrapper for structure material of porous media including reactive reference porosity

\maintainer Christoph Schmidt

\level 3
 *----------------------------------------------------------------------*/


#ifndef STRUCTPORO_REACTION_H_
#define STRUCTPORO_REACTION_H_

#include "structporo.H"

namespace MAT
{
  // forward declaration
  class StructPoroReaction;

  namespace PAR
  {
    class StructPoroReaction : public PAR::StructPoro
    {
      friend class MAT::StructPoroReaction;

     public:
      /// standard constructor
      StructPoroReaction(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~StructPoroReaction() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

      /// @name material parameters
      //@{

      int dofIDReacScalar_;
      //@}

    };  // class StructPoro

  }  // namespace PAR

  class StructPoroReactionType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "StructPoroReactionType"; }

    static StructPoroReactionType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static StructPoroReactionType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for StructPoro material
  ///
  /// This object exists (several times) at every element
  class StructPoroReaction : public StructPoro
  {
   public:
    /// construct empty material object
    StructPoroReaction();

    /// construct the material object given material parameters
    explicit StructPoroReaction(MAT::PAR::StructPoroReaction* params);

    //! @name Packing and Unpacking

    /*!
     \brief Return unique ParObject id

     every class implementing ParObject needs a unique id defined at the
     top of drt_parobject.H (this file) and should return it in this method.
     */
    virtual int UniqueParObjectId() const
    {
      return StructPoroReactionType::Instance().UniqueParObjectId();
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
      return INPAR::MAT::m_structpororeaction;
    }

    /// return initial porosity
    double RefPorosityTimeDeriv() const { return dphiDphiref_ * refporositydot_; }

    /// compute current porosity and save it
    void ComputePorosity(Teuchos::ParameterList& params,  ///< (i) element parameter list
        double press,                                     ///< (i) pressure at gauss point
        double J,          ///< (i) determinant of jacobian at gauss point
        int gp,            ///< (i) number of current gauss point
        double& porosity,  ///< (o) porosity at gauss point
        double* dphi_dp,   ///< (o) first derivative of porosity w.r.t. pressure at gauss point
        double* dphi_dJ,   ///< (o) first derivative of porosity w.r.t. jacobian at gauss point
        double*
            dphi_dJdp,  ///< (o) derivative of porosity w.r.t. pressure and jacobian at gauss point
        double* dphi_dJJ,  ///< (o) second derivative of porosity w.r.t. jacobian at gauss point
        double* dphi_dpp,  ///< (o) second derivative of porosity w.r.t. pressure at gauss point
        bool save = true);

    //! evaluate constitutive relation for porosity and compute derivatives
    virtual void ConstitutiveDerivatives(Teuchos::ParameterList& params,  ///< (i) parameter list
        double press,        ///< (i) fluid pressure at gauss point
        double J,            ///< (i) Jacobian determinant at gauss point
        double porosity,     ///< (i) porosity at gauss point
        double* dW_dp,       ///< (o) derivative of potential w.r.t. pressure
        double* dW_dphi,     ///< (o) derivative of potential w.r.t. porosity
        double* dW_dJ,       ///< (o) derivative of potential w.r.t. jacobian
        double* dW_dphiref,  ///< (o) derivative of potential w.r.t. jacobian
        double* W            ///< (o) inner potential
    );

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new StructPoroReaction(*this));
    }

    /// Initialize internal variables
    void Setup(int numgp,  ///< number of Gauss points
        DRT::INPUT::LineDefinition* linedef);

    /// Return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const { return params_; }

    /// return reference porosity average
    double RefPorosityAv() const;

    //! @name Evaluation methods

    /// evaluate material law
    virtual void Evaluate(const LINALG::Matrix<3, 3>* defgrd,  ///< (i) deformation gradient
        const LINALG::Matrix<6, 1>* glstrain,                  ///< (i) green lagrange strain
        Teuchos::ParameterList& params,                        ///< (i) parameter list
        LINALG::Matrix<6, 1>* stress,  ///< (o) second piola kirchhoff stress
        LINALG::Matrix<6, 6>* cmat,    ///< (o) constitutive matrix
        const int eleGID);

    //@}

    //! @name Visualization methods

    /// Return names of visualization data
    virtual void VisNames(std::map<std::string, int>& names);

    /// Return visualization data
    virtual bool VisData(const std::string& name, std::vector<double>& data, int numgp, int eleID);

    //@}

   protected:
    virtual void Reaction(const double porosity, const double J,
        Teuchos::RCP<std::vector<double>> scalars, Teuchos::ParameterList& params);

    /// my material parameters
    MAT::PAR::StructPoroReaction* params_;

    /// reference porosity
    double refporosity_;

    /// derivative of porosity w.r.t. reference porosity
    double dphiDphiref_;

    /// time derivative of reference porosity
    double refporositydot_;
  };

}  // namespace MAT


#endif /* STRUCTPORO_REACTION_H_ */
