/*----------------------------------------------------------------------------*/
/** \file
\brief compute the intersection of the zero level-set iso-contour on the
  SCATRA surface discretization and evaluate important quantities


\level 3
*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_XCONTACT_XCONTACT_LEVELSET_INTERSECTION_H_
#define SRC_DRT_CONTACT_XCONTACT_XCONTACT_LEVELSET_INTERSECTION_H_

#include "../drt_levelset/levelset_intersection_utils.H"

namespace XCONTACT
{
  namespace LEVELSET
  {
    class Intersection : public SCATRA::LEVELSET::Intersection
    {
     public:
      /// constructor
      Intersection();

      /// destructor
      virtual ~Intersection(){};

      /** \brief capture the zero iso line of the levelset function without calculation
       *  of the plus and minus volumes and boundary surface
       *
       *  \param phi                     (in)  : levelset values
       *  \param scatradis               (in)  : scalar transport discretization
       *  \param elementBoundaryIntCells (out) : map containg the element boundary cell pointers
       *                                         (key is the ID of the parent element ( coming from
       *                                         the ScaTra discretization)
       *
       *  \author hiermeier \date 11/16 */
      void CaptureZeroLevelSet(const Epetra_Vector& phi, const DRT::Discretization& scatradis,
          std::map<int, GEO::BoundaryIntCellPtrs>& elementBoundaryIntCells);

      /** \brief capture the zero iso line of the levelset function and calculate
       *  the plus and minus volumes and boundary surface
       *
       *  \param phi                     (in)  : levelset values
       *  \param scatradis               (in)  : scalar transport discretization
       *  \param volumedomainminus       (out) : global volume of the minus domain
       *  \param volumedomainplus        (out) : global volume of the plus domain
       *  \param zerosurface             (out) : global surface value of the boundary cells
       *  \param elementBoundaryIntCells (out) : map containg the element boundary cell pointers
       *                                         (key is the ID of the parent element ( coming from
       *                                         the ScaTra discretization)
       *
       *  \author hiermeier \date 11/16 */
      void CaptureZeroLevelSet(const Epetra_Vector& phi, const DRT::Discretization& scatradis,
          double& volumedomainminus, double& volumedomainplus, double& zerosurface,
          std::map<int, GEO::BoundaryIntCellPtrs>& elementBoundaryIntCells);

     protected:
      /// lists the supported boundary cell types
      virtual void CheckBoundaryCellType(DRT::Element::DiscretizationType distype_bc) const;


      virtual void AddToBoundaryIntCellsPerEle(const LINALG::SerialDenseMatrix& xyze,
          const GEO::CUT::BoundaryCell& bcell, DRT::Element::DiscretizationType distype_ele);
    };  // class Intersection

  }  // namespace LEVELSET
}  // namespace XCONTACT


#endif /* SRC_DRT_CONTACT_XCONTACT_XCONTACT_LEVELSET_INTERSECTION_H_ */
