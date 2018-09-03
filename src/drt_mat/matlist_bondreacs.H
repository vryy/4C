/*----------------------------------------------------------------------*/
/*!
 \file matlist_bondreacs.H

 \brief material list for bond reactions.

 \level 2

 <pre>
   \maintainer Andreas Rauch
               rauch@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289 - 15240
 </pre>
 *----------------------------------------------------------------------*/

#ifndef MATLIST_BONDREACS_H_
#define MATLIST_BONDREACS_H_



#include "matpar_parameter.H"
#include "material.H"
#include "matlist.H"
#include "matlist_reactions.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace MAT
{
  namespace PAR
  {
    /*----------------------------------------------------------------------*/
    /// material parameters for list of materials
    class MatListBondReacs : public MatListReactions
    {
     public:
      /// standard constructor
      MatListBondReacs(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~MatListBondReacs() { ; }

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial();

    };  // class MatListBondReacs

  }  // namespace PAR

  class MatListBondReacsType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "MatListBondReacsType"; }

    static MatListBondReacsType& Instance() { return instance_; };

    virtual DRT::ParObject* Create(const std::vector<char>& data);

   private:
    static MatListBondReacsType instance_;
  };

  /*----------------------------------------------------------------------*/
  /// Wrapper for a list of materials
  class MatListBondReacs : public MatListReactions
  {
   public:
    /// construct empty material object
    MatListBondReacs();

    /// construct the material object given material parameters
    explicit MatListBondReacs(MAT::PAR::MatListBondReacs* params);

    //! @name Packing and Unpacking

    /*!
        \brief Return unique ParObject id

        every class implementing ParObject needs a unique id defined at the
        top of drt_parobject.H (this file) and should return it in this method.
     */
    virtual int UniqueParObjectId() const
    {
      return MatListBondReacsType::Instance().UniqueParObjectId();
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

    /// initialize
    virtual void Initialize();

    /// material type
    virtual INPAR::MAT::MaterialType MaterialType() const
    {
      return INPAR::MAT::m_matlist_bondreacs;
    }

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const
    {
      return Teuchos::rcp(new MatListBondReacs(*this));
    }

    /// Return quick accessible material parameter data
    virtual MAT::PAR::MatListBondReacs* Parameter() const { return paramsbondreac_; }

    /// advanced reaction terms
    double CalcReaBodyForceTerm(const int k, const std::vector<double>& phinp,
        const std::vector<double>& phin, const double violation, const double porosity,
        const double* gpcoord,  //!< Gauss-point coordinates
        const double scale = 1.0) const;

    /// calculate advanced reaction term derivatives
    void CalcReaBodyForceDerivMatrix(const int k, std::vector<double>& derivs,
        const std::vector<double>& phinp, const std::vector<double>& phin, const double violation,
        const double porosity,
        const double* gpcoord,  //!< Gauss-point coordinates
        const double scale = 1.0) const;

    /// advanced reaction terms
    double CalcReaBodyForceTerm(const int k, const std::vector<double>& phinp,
        const std::vector<double>& phin,
        const std::vector<std::pair<std::string, double>>& constants, const double violation,
        const double porosity,
        const double* gpcoord,  //!< Gauss-point coordinates
        const double scale = 1.0) const;

    /// calculate advanced reaction term derivatives
    void CalcReaBodyForceDerivMatrix(const int k, std::vector<double>& derivs,
        const std::vector<double>& phinp, const std::vector<double>& phin,
        const std::vector<std::pair<std::string, double>>& constants, const double violation,
        const double porosity,
        const double* gpcoord,  //!< Gauss-point coordinates
        const double scale = 1.0) const;


   private:
    /// clear everything
    void Clear();

    /// my material parameters
    MAT::PAR::MatListBondReacs* paramsbondreac_;

  };  // class MatListBondReacs

}  // namespace MAT


#endif /* MATLIST_BONDREACS_H_ */
