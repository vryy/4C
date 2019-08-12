/*----------------------------------------------------------------------*/
/*!
\brief Classes for mortar coupling in 2D.

\level 2

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------*/
#ifndef MORTAR_COUPLING2D_H
#define MORTAR_COUPLING2D_H

#include <vector>
#include "Epetra_Comm.h"
#include "../drt_inpar/inpar_mortar.H"

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace MORTAR
{
  // forward declarations
  class MortarElement;
  class ParamsInterface;

  /*!
   \brief A class representing the framework for mortar coupling of ONE
   slave element and ONE master element of a mortar interface in
   2D. Concretely, this class controls projection, overlap detection
   and finally integration of the mortar coupling matrices D and M
   and possibly of the weighted gap vector g~.

   \author popp (popp@lnm.mw.tum.de)
   */

  class Coupling2d
  {
   public:
    /*!
     \brief Constructor with shape function specification

     Constructs an instance of this class and enables custom shape function types.<br>
     Note that this is \b not a collective call as coupling is
     performed in parallel by individual processes.

     */
    Coupling2d(DRT::Discretization& idiscret, int dim, bool quad, Teuchos::ParameterList& params,
        MORTAR::MortarElement& sele, MORTAR::MortarElement& mele);

    /*!
     \brief Destructor

     */
    virtual ~Coupling2d() {}

    //! @name Access methods

    /*!
     \brief Get interface discretization

     */
    virtual DRT::Discretization& Discret() const { return idiscret_; }

    /*!
     \brief Get coupling slave element

     */
    virtual MORTAR::MortarElement& SlaveElement() const { return sele_; }

    /*!
     \brief Get coupling master element

     */
    virtual MORTAR::MortarElement& MasterElement() const { return mele_; }

    /*!
     \brief Get problem dimension (here: 2D)

     */
    virtual const int& Dim() { return dim_; };

    /*!
     \brief Get coupling / FE ansatz type (true = quadratic)

     */
    virtual const bool& Quad() { return quad_; };

    /*!
     \brief Return the LM interpolation / testing type for quadratic FE

     */
    INPAR::MORTAR::LagMultQuad LagMultQuad()
    {
      return DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(imortar_, "LM_QUAD");
    }

    /*!
     \brief Return the LM shape fcn type

     */
    INPAR::MORTAR::ShapeFcn ShapeFcn()
    {
      return DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(imortar_, "LM_SHAPEFCN");
    }

    /*!
     \brief Get interface contact parameter list

     */
    virtual Teuchos::ParameterList& InterfaceParams() { return imortar_; };

    /*!
     \brief Get projection status of the four end nodes

     */
    virtual const std::vector<bool>& HasProj() { return hasproj_; }

    /*!
     \brief Get overlap regions in parameter spaces

     */
    virtual const std::vector<double>& XiProj() { return xiproj_; }

    /*!
     \brief Get overlap status

     */
    virtual const bool& Overlap() { return overlap_; };

    //@}

    //! @name Evaluation methods

    /*!
     \brief Projection of slave / master pair

     This method projects the nodes of the slave MortarElement sele_ onto
     the master MortarElement mele_ and vice versa. The variable hasproj
     stores a boolean variable for each of the 4 end nodes, indicating
     whether a feasible projection was found or not. The local element
     coordinates of the 4 projection points are stored in xiproj.

     */
    virtual bool Project();

    /*!
     \brief Detect overlap of slave / master pair

     This method evaluates the overlap of the current MortarElement pair
     sele_ / mele_ based on the projection status of the 4 end nodes
     (hasproj) and the coordinates of the projection points (xiproj).
     According to the detected overlap case, the integration limits
     are determined and written into xiproj and the overlap status
     is stored in the boolean variable overlap.

     */
    virtual bool DetectOverlap();

    /*!
     \brief Integrate overlap of slave / master pair

     This method integrates the overlap of the current MortarElement
     pair sele_ / mele_ based on the integration limits (xiproj). The
     integration always includes the Mortar matrices D/M and the gap g.

     */
    virtual bool IntegrateOverlap(const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);

    //@}
   private:
    /*!
     \brief Checks roughly whether the two elements are near

     This methods checks the orientation of slave and master element.
     Projection and overlap detection only make sense if the scalar
     product of the two center normals is negative (i.e. if the two
     elements form an angle smaller than 90°).

     */
    virtual bool RoughCheckOrient();

   protected:
    /*!
     \brief Get communicator

     */
    virtual const Epetra_Comm& Comm() const;

    // don't want = operator and cctor
    Coupling2d operator=(const Coupling2d& old);
    Coupling2d(const Coupling2d& old);

    DRT::Discretization& idiscret_;    // discretization of the contact interface
    int dim_;                          // problem dimension (here: 2D)
    bool quad_;                        // flag indicating coupling type (true = quadratic)
    Teuchos::ParameterList& imortar_;  // containing contact input parameters
    MORTAR::MortarElement& sele_;      // slave element to perform coupling for
    MORTAR::MortarElement& mele_;      // master element to perform coupling for

    std::vector<bool> hasproj_;   // projection status of the four end nodes
    std::vector<double> xiproj_;  // overlap regions in parameter spaces
    bool overlap_;                // overlap status
  };
  // class Coupling2d

  /*!
   \brief A class representing the framework for mortar coupling of ONE
   slave element and SEVERAL master elements of a mortar interface in
   2D. This class simply stores and manages several Coupling2d objects.

   \author popp (popp@lnm.mw.tum.de)
   */

  class Coupling2dManager
  {
   public:
    /*!
     \brief Constructor with shape function specification

     Constructs an instance of this class and enables custom shape function types.<br>
     Note that this is \b not a collective call as coupling is
     performed in parallel by individual processes.

     */
    Coupling2dManager(DRT::Discretization& idiscret, int dim, bool quad,
        Teuchos::ParameterList& params, MORTAR::MortarElement* sele,
        std::vector<MORTAR::MortarElement*> mele);

    /*!
     \brief Destructor

     */
    virtual ~Coupling2dManager() {}

    /*!
     \brief Get coupling slave element

     */
    virtual MORTAR::MortarElement& SlaveElement() const { return *sele_; }

    /*!
     \brief Get one specific coupling master element

     */
    virtual MORTAR::MortarElement& MasterElement(int k) const { return *(mele_[k]); }

    /*!
     \brief Get all coupling master elements

     */
    virtual std::vector<MORTAR::MortarElement*> MasterElements() const { return mele_; }

    /*!
     \brief Get coupling pairs

     */
    virtual std::vector<Teuchos::RCP<MORTAR::Coupling2d>>& Coupling() { return coup_; }

    /*!
     \brief Get type of quadratic LM interpolation

     */
    virtual INPAR::MORTAR::LagMultQuad LagMultQuad()
    {
      return DRT::INPUT::IntegralValue<INPAR::MORTAR::LagMultQuad>(imortar_, "LM_QUAD");
    }

    /*!
     \brief Get integration type

     */
    virtual INPAR::MORTAR::IntType IntType()
    {
      return DRT::INPUT::IntegralValue<INPAR::MORTAR::IntType>(imortar_, "INTTYPE");
    }

    /*!
     \brief Evaluate coupling pairs

     */
    virtual bool EvaluateCoupling(const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);

    /*!
     \brief Get coupling type

     */
    virtual const bool& Quad() { return quad_; };

    /*!
     \brief Return the LM shape fcn type

     */
    INPAR::MORTAR::ShapeFcn ShapeFcn()
    {
      return DRT::INPUT::IntegralValue<INPAR::MORTAR::ShapeFcn>(imortar_, "LM_SHAPEFCN");
    }

    //@}
   protected:
    /*!
     \brief Evaluate mortar coupling pairs

     */
    virtual void IntegrateCoupling(const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);

    /*!
     \brief Calculate consistent dual shape functions in boundary elements

     */
    virtual void ConsistDualShape();

   protected:
    // don't want = operator and cctor
    Coupling2dManager operator=(const Coupling2dManager& old);
    Coupling2dManager(const Coupling2dManager& old);

    DRT::Discretization& idiscret_;             // discretization of the contact interface
    int dim_;                                   // problem dimension (here: 2D)
    bool quad_;                                 // flag indicating coupling type (true = quadratic)
    Teuchos::ParameterList& imortar_;           // containing contact input parameters
    MORTAR::MortarElement* sele_;               // slave element
    std::vector<MORTAR::MortarElement*> mele_;  // master elements
    std::vector<Teuchos::RCP<Coupling2d>> coup_;  // coupling pairs
  };
  // class Coupling2dManager
}  // namespace MORTAR

#endif  // #ifndef MORTAR_COUPLING2D_H
