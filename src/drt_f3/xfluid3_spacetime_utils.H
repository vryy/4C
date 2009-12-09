/*!
\file xfluid3_utils.H

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET
#ifndef XFLUID3_SPACETIME_UTILS_H
#define XFLUID3_SPACETIME_UTILS_H

#include "../drt_lib/drt_timecurve.H"
#include "../drt_xfem/xfem_enums.H"
#include "../drt_geometry/intersection_service.H"
#include <Epetra_SerialDenseSolver.h>

namespace XFLUID
{

  /// generate old velocity/acceleration values, if integration point was in a void during the last time step
  inline LINALG::Matrix<3,1> getIntegrationPointValueFromBoundary(
      const Teuchos::RCP<DRT::Discretization>&   boundarydis,   ///< connection to the boundary dis
      const Epetra_Vector&                       ivalcol,       ///< global vector with boundary dofs
      const DRT::Element*                        boundaryele,   ///< the element of which element tangent stiffness matrix is calculated
      const LINALG::Matrix<8,1>&                 funct_ST,      ///< shape function of the space time slab
      const std::size_t                          numnode_slab   ///< number of nodes of space time slab
  )
  {
    const std::size_t nsd = 3;
    const std::size_t numnodefix_slab = 8;

    const std::size_t numnode_boundary = boundaryele->NumNode();

    // get interface velocities at the boundary element nodes
    LINALG::Matrix<nsd,numnodefix_slab> val_boundary;

    const DRT::Node*const* nodes = boundaryele->Nodes();

    std::vector<double> myval(nsd);
    for (std::size_t inode = 0; inode < numnode_boundary; ++inode)
    {
      const DRT::Node* node = nodes[inode];
      const std::vector<int> lm = boundarydis->Dof(node);

      DRT::UTILS::ExtractMyValues(ivalcol,myval,lm);
      for (std::size_t isd = 0; isd < nsd; ++isd)
        val_boundary(isd,inode  ) = myval[isd];

      DRT::UTILS::ExtractMyValues(ivalcol,myval,lm);
      for (std::size_t isd = 0; isd < nsd; ++isd)
        val_boundary(isd,inode+4) = myval[isd];
    }

    return XFLUID::interpolateVectorFieldToIntPoint(val_boundary , funct_ST, numnode_slab);
  }

  /// generate old velocity/acceleration values, if integration point was in a void during the last time step
  template<DRT::Element::DiscretizationType DISTYPE, class M>
  bool modifyOldTimeStepsValues(
      const DRT::Element*                        ele,       ///< the element of which element tangent stiffness matrix is calculated
      const Teuchos::RCP<XFEM::InterfaceHandleXFSI>&  ih,   ///< connection to the interface handler
      const M&                                   xyze,
      const LINALG::Matrix<3,1>&                 posXiDomain,
      const bool                                 is_in_fluid,
      const bool                                 was_in_fluid,
      double                                     dt,
      const Epetra_Vector&                       ivelcolnp,
      const Epetra_Vector&                       ivelcoln,
      const Epetra_Vector&                       ivelcolnm,
      const Epetra_Vector&                       iacccoln,
      LINALG::Matrix<3,1>&                       gpveln,
      LINALG::Matrix<3,1>&                       gpvelnm,
      LINALG::Matrix<3,1>&                       gpaccn
  )
  {
    const std::size_t nsd = 3;
    const std::size_t numnodefix_slab = 8;

    LINALG::Matrix<nsd,1> posx_gp;
    GEO::elementToCurrentCoordinatesT<DISTYPE>(xyze, posXiDomain, posx_gp);

    if (not is_in_fluid)
    {
      std::cout << "should I arrive here?" << std::endl;
      exit(1);
//      return false;
    }

    const bool in_space_time_slab_area = (is_in_fluid and (not was_in_fluid));

    if (in_space_time_slab_area)
    {
      XFEM::SpaceTimeBoundaryCell slab;
      LINALG::Matrix<nsd,1> rst(true);

      const bool found_cell = ih->FindSpaceTimeLayerCell(posx_gp,slab,rst);

      if (found_cell)
      {
        const std::size_t numnode_slab = DRT::UTILS::getNumberOfElementNodes(slab.Shape());

        if (numnode_slab > numnodefix_slab)
          dserror("increase numnodefix_slab!");

        const double delta_slab = -(rst(2)-1.0)*0.5;

        if (not ((-1.0e-7) < delta_slab and delta_slab < (1.0+1.0e-7)))
        {
          cout << rst(2) <<  "  " << delta_slab << endl;
          cout << slab.toString() << endl << endl;
          dserror("wrong value of delta_slab");
        }

        const DRT::Element* boundaryele = ih->cutterdis()->gElement(slab.getBeleId());
        const size_t numnode_boundary = boundaryele->NumNode();
        if (numnode_boundary != 4)
          dserror("needs more generalizashun!");
        if (ivelcoln.GlobalLength() > (int)nsd)
        {
          //LINALG::SerialDenseVector funct_ST(numnode_boundary*2);
          LINALG::Matrix<numnodefix_slab,1> funct_ST;
          DRT::UTILS::shape_function_3D(funct_ST,rst(0),rst(1),rst(2),slab.Shape());

          gpveln  = getIntegrationPointValueFromBoundary(ih->cutterdis(), ivelcoln,  boundaryele, funct_ST, numnode_slab);
          gpvelnm = getIntegrationPointValueFromBoundary(ih->cutterdis(), ivelcolnm, boundaryele, funct_ST, numnode_slab);
          gpaccn  = getIntegrationPointValueFromBoundary(ih->cutterdis(), iacccoln,  boundaryele, funct_ST, numnode_slab);

//          const std::size_t nsd = 3;
//          const std::size_t numnodefix_slab = 8;
//
//          const std::size_t numnode_boundary = boundaryele->NumNode();
//
//          // get interface velocities at the boundary element nodes
//          LINALG::Matrix<nsd,numnodefix_slab> val_boundary;
//
//          const DRT::Node*const* nodes = boundaryele->Nodes();
//
//          std::vector<double> myval(nsd);
//          for (std::size_t inode = 0; inode < numnode_boundary; ++inode)
//          {
//            const DRT::Node* node = nodes[inode];
//            const std::vector<int> lm = ih->cutterdis()->Dof(node);
//
//            DRT::UTILS::ExtractMyValues(ivelcoln,myval,lm);
//            for (std::size_t isd = 0; isd < nsd; ++isd)
//              val_boundary(isd,inode  ) = myval[isd];
//
//            DRT::UTILS::ExtractMyValues(ivelcoln,myval,lm);
//            for (std::size_t isd = 0; isd < nsd; ++isd)
//              val_boundary(isd,inode+numnode_boundary) = myval[isd];
//          }
//
//          gpveln = XFLUID::interpolateVectorFieldToIntPoint(val_boundary , funct_ST, numnode_slab);

        }
        else
        {
          gpveln.Clear();
          gpvelnm.Clear();
          gpaccn.Clear();
        }
      }
      else
      {
//        cout << "did not find cell" << endl;
      }
    }
    else
    {
    }
    return true;
  }


  /// generate old velocity/acceleration values, if integration point was in a void during the last time step
    template<DRT::Element::DiscretizationType DISTYPE, class M>
    bool ProjectSpaceTimeValuesToNewMesh(
        const DRT::Element*                        ele,       ///< the element of which element tangent stiffness matrix is calculated
        const Teuchos::RCP<XFEM::InterfaceHandleXFSI>&  ih,   ///< connection to the interface handler
        const M&                                   xyze,
        const LINALG::Matrix<3,1>&                 posXiDomain,
        const int                                  labelnp,
        const Epetra_Vector&                       ivelcoln,
        const Epetra_Vector&                       iacccoln,
        LINALG::Matrix<3,1>&                       gpveln,
        LINALG::Matrix<3,1>&                       gpaccn
    )
    {
      const std::size_t nsd = 3;
      const std::size_t numnodefix_slab = 8;

      LINALG::Matrix<nsd,1> posx_gp;
      GEO::elementToCurrentCoordinatesT<DISTYPE>(xyze, posXiDomain, posx_gp);

      const bool is_in_fluid = (0 == labelnp);

      if (not is_in_fluid)
      {
        std::cout << "should I arrive here?" << std::endl;
        exit(1);
  //      return false;
      }

      const bool was_in_fluid = (ih->PositionWithinConditionN(posx_gp) == 0);

      const bool in_space_time_slab_area = (is_in_fluid and (not was_in_fluid));

      if (in_space_time_slab_area)
      {
        XFEM::SpaceTimeBoundaryCell slab;
        LINALG::Matrix<nsd,1> rst(true);

        const bool found_cell = ih->FindSpaceTimeLayerCell(posx_gp,slab,rst);

        if (found_cell)
        {
          const std::size_t numnode_slab = DRT::UTILS::getNumberOfElementNodes(slab.Shape());

          if (numnode_slab > numnodefix_slab)
            dserror("increase numnodefix_slab!");

          const double delta_slab = -(rst(2)-1.0)*0.5;

          if (not ((-1.0e-7) < delta_slab and delta_slab < (1.0+1.0e-7)))
          {
            cout << rst(2) <<  "  " << delta_slab << endl;
            cout << slab.toString() << endl << endl;
            dserror("wrong value of delta_slab");
          }

          const DRT::Element* boundaryele = ih->cutterdis()->gElement(slab.getBeleId());
          const size_t numnode_boundary = boundaryele->NumNode();
          if (numnode_boundary != 4)
            dserror("needs more generalizashun!");
          if (ivelcoln.GlobalLength() > (int)nsd)
          {
            //LINALG::SerialDenseVector funct_ST(numnode_boundary*2);
            LINALG::Matrix<numnodefix_slab,1> funct_ST;
            DRT::UTILS::shape_function_3D(funct_ST,rst(0),rst(1),rst(2),slab.Shape());

            gpveln  = getIntegrationPointValueFromBoundary(ih->cutterdis(), ivelcoln,  boundaryele, funct_ST, numnode_slab);
            gpaccn  = getIntegrationPointValueFromBoundary(ih->cutterdis(), iacccoln,  boundaryele, funct_ST, numnode_slab);

          }
          else
          {
            gpveln.Clear();
            gpaccn.Clear();
          }
        }
        else
        {
  //        cout << "did not find cell" << endl;
        }
      }

      return true;
    }

} // namespace XFLUID

#endif
#endif  // #ifdef CCADISCRET
