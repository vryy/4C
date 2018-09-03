/*----------------------------------------------------------------------*/
/*!
\file beam_link.H

\brief One beam-to-beam pair (two beam elements) connected by a mechanical link

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------*/
#ifndef BeamLink_H
#define BeamLink_H

#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "../drt_inpar/inpar_beaminteraction.H"

//#include <Teuchos_RCP.hpp>

// forward declarations
// class Epetra_Vector;
// class Epetra_SerialDenseMatrix;
namespace DRT
{
  class Element;
  class PackBuffer;
  namespace ELEMENTS
  {
    class Beam3Base;
    class Beam3r;
  }  // namespace ELEMENTS
}  // namespace DRT
namespace LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace LINALG

namespace BEAMINTERACTION
{
  class BeamLinkType : public DRT::ParObjectType
  {
   public:
    std::string Name() const { return "BeamLinkType"; };

    static BeamLinkType& Instance() { return instance_; };

   private:
    static BeamLinkType instance_;
  };


  /*!
   \brief element for interaction of two 3D beam elements via a mechanical linkage
   */
  class BeamLink : public DRT::ParObject
  {
   public:
    //! @name Constructors and destructors and related methods

    //! Constructor
    BeamLink();

    /*!
    \brief Copy Constructor

    Makes a deep copy of a Element

    */
    BeamLink(const BeamLink& old);

    //! Destructor
    virtual ~BeamLink() {}

    //! Initialization
    virtual void Init(const int id, const std::vector<std::pair<int, int>>& eleids,
        const std::vector<LINALG::Matrix<3, 1>>& initpos,
        const std::vector<LINALG::Matrix<3, 3>>& inittriad,
        INPAR::BEAMINTERACTION::CrosslinkerType linkertype, double timelinkwasset);

    //! Setup
    virtual void Setup(const int matnum);

    /*!
    \brief Return unique ParObject id

    Every class implementing ParObject needs a unique id defined at the
    top of drt_parobject.H
    */
    virtual int UniqueParObjectId() const = 0;

    //! return copy of this linking object object
    virtual Teuchos::RCP<BeamLink> Clone() const = 0;

    /*!
    \brief Pack this class so it can be communicated

    \ref Pack and \ref Unpack are used to communicate this element

    */
    virtual void Pack(DRT::PackBuffer& data) const;

    /*!
    \brief Unpack data from a char vector into this class

    \ref Pack and \ref Unpack are used to communicate this element

    */
    virtual void Unpack(const std::vector<char>& data);

    //@}

    //! @name Access methods

    /*!
    \brief Return global id
    */
    virtual inline int Id() const { return id_; }

    /*!
    \brief Return gid of first/second element (specified via given local element number 0/1)
    */
    virtual inline const int& GetEleGid(const int& elenum) const { return bspotIds_[elenum].first; }

    /*!
    \brief Return element-local binding spot number of first/second element
           (specified via given local element number 0/1)
    */
    virtual inline const int& GetLocBSpotNum(const int& elenum) const
    {
      return bspotIds_[elenum].second;
    }

    //! return position of first connection site
    inline const LINALG::TMatrix<double, 3, 1>& GetBindSpotPos1() const { return bspotpos1_; }

    //! return position of second connection site
    inline const LINALG::TMatrix<double, 3, 1>& GetBindSpotPos2() const { return bspotpos2_; }

    inline INPAR::BEAMINTERACTION::CrosslinkerType GetLinkerType() const { return linkertype_; }

    //! return time at which linker was set
    inline double GetTimeLinkWasSet() const { return timelinkwasset_; }

    //! return time at which linker was set
    inline double GetReferenceLength() const { return reflength_; }

    //! get force in first or second binding spot
    virtual void GetBindingSpotForce(int bspotid, LINALG::SerialDenseVector& bspotforce) const = 0;

    //! get internal linker energy
    virtual double GetInternalEnergy() const = 0;

    //! get kinetic linker energy
    virtual double GetKineticEnergy() const = 0;

    //! scale linker element reference length
    virtual void ScaleLinkerReferenceLength(double scalefac) = 0;

    //@}

    //! @name Public evaluation methods

    /*
    \brief Update position and triad of both connection sites (a.k.a. binding spots)
    */
    virtual void ResetState(
        std::vector<LINALG::Matrix<3, 1>>& bspotpos, std::vector<LINALG::Matrix<3, 3>>& bspottriad);

    /*!
    \brief Evaluate forces
    */
    virtual bool EvaluateForce(
        LINALG::SerialDenseVector& forcevec1, LINALG::SerialDenseVector& forcevec2) = 0;

    /*!
    \brief Evaluate stiffness contribution
    */
    virtual bool EvaluateStiff(LINALG::SerialDenseMatrix& stiffmat11,
        LINALG::SerialDenseMatrix& stiffmat12, LINALG::SerialDenseMatrix& stiffmat21,
        LINALG::SerialDenseMatrix& stiffmat22) = 0;

    /*!
    \brief Evaluate forces and stiffness contribution
    */
    virtual bool EvaluateForceStiff(LINALG::SerialDenseVector& forcevec1,
        LINALG::SerialDenseVector& forcevec2, LINALG::SerialDenseMatrix& stiffmat11,
        LINALG::SerialDenseMatrix& stiffmat12, LINALG::SerialDenseMatrix& stiffmat21,
        LINALG::SerialDenseMatrix& stiffmat22) = 0;


    void Print(std::ostream& out) const;

    //@}

   protected:
    //! returns init state
    inline const bool& IsInit() const { return isinit_; };

    //! returns setup state
    inline const bool& IsSetup() const { return issetup_; };

    //! Check the init state
    inline void CheckInit() const
    {
      if (not IsInit()) dserror("Call Init() first!");
    }

    //! Check the init and setup state
    inline void CheckInitSetup() const
    {
      if (not IsInit() or not IsSetup()) dserror("Call Init() and Setup() first!");
    }



   protected:
    //! @name member variables

    //! indicates if the Init() function has been called
    bool isinit_;

    //! indicates if the Setup() function has been called
    bool issetup_;

   private:
    //! a unique global id
    int id_;

    //! unique identifiers for first [0] and second [1] binding spot:
    // each is a pair of element GID and local binding spot number
    std::vector<std::pair<int, int>> bspotIds_;

    //! current position of the two connection sites (a.k.a. binding spots)
    LINALG::Matrix<3, 1> bspotpos1_;
    LINALG::Matrix<3, 1> bspotpos2_;

    //! type of filament element belongs to
    INPAR::BEAMINTERACTION::CrosslinkerType linkertype_;

    //! stores the the time the link was set (can e.g. be used to calculate
    //  lifetime of a link or check if link is new in certain time step)
    double timelinkwasset_;

    //! linker reference length
    double reflength_;

    //@}
  };

}  // namespace BEAMINTERACTION

#endif
