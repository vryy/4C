/*----------------------------------------------------------------------*/
/*!
 \file structporo_reaction_ecm.H

 \brief porous material for dissolution reaction in ECM model

 <pre>
 \level 3

 \maintainer Anh-Tu Vuong
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
  virtual ~StructPoroReactionECM(){};

  /// create material instance of matching type with my parameters
  virtual Teuchos::RCP<MAT::Material> CreateMaterial();

  /// @name material parameters
  //@{

  /// density of collagen
  double densCollagen_;
  //@}

};
// class StructPoroReactionECM

}// namespace PAR

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
  void Setup(int numgp, ///< number of Gauss points
      DRT::INPUT::LineDefinition* linedef);

  /// Return quick accessible material parameter data
  virtual MAT::PAR::Parameter* Parameter() const
  {
    return params_;
  }

  /// evaluate chemical potential
  virtual void ChemPotential(
      const LINALG::Matrix<6, 1>& glstrain, ///< (i) green lagrange strain
      const double porosity, ///< (i) porosity
      const double press, ///< (i) pressure at gauss point
      const double J, ///< (i) determinant of jacobian at gauss point
      int EleID, ///< (i) element GID
      double& psi, ///< (o) chemical potential
      const int gp
      );

  /// Update of GP data (e.g., history variables)
  virtual void Update();

  //! @name Visualization methods

  /// Return names of visualization data
  virtual void VisNames(std::map<std::string, int>& names);

  /// Return visualization data
  virtual bool VisData(const std::string& name, std::vector<double>& data,
      int numgp, int eleID);

protected:

  virtual void Reaction(
      const double porosity,
      const double J,
      Teuchos::RCP<std::vector<double> > scalars,
      Teuchos::ParameterList& params);

  /// reference porosity at time step n
  double refporosity_old_;

  /// time derivative of reference porosity at time step n
  double refporositydot_old_;

  /// chemical potential
  std::vector <double> chempot_;

  /// initial chemical potential
  std::vector <double> chempot_init_;

  /// my material parameters
  MAT::PAR::StructPoroReactionECM* params_;

};
}

#endif /* STRUCTPORO_REACTION_ECM_H_ */
