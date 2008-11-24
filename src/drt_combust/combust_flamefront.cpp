/*!-----------------------------------------------------------------------------------------------*
 \file combust_flamefront.cpp

 \brief 

  detailed description in header file combust_interface.H

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
 *------------------------------------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "combust_flamefront.H"
#include "../drt_lib/standardtypes_cpp.H"
// #include "../drt_lib/drt_utils.H"
// #include "../drt_geometry/integrationcell.H"

// extern struct _FILES  allfiles;

/*------------------------------------------------------------------------------------------------*
 | constructor                                                                         henke 10/08 | 
 *------------------------------------------------------------------------------------------------*/
COMBUST::FlameFront::FlameFront(
    Teuchos::RCP<const DRT::Discretization>  fluiddis,
    Teuchos::RCP<const DRT::Discretization>  gfuncdis,
    map< int, GEO::DomainIntCells >&            fluidintcells,            ///< domainintegrationcells for each intersected element
    map< int, GEO::BoundaryIntCells >&          flamefrontintcells        ///< boundaryintegrationcells for each intersected element
    )
{
  if (fluiddis->Comm().MyPID() == 0)
    std::cout << "Constructing flame front" << std::endl;
}
/*------------------------------------------------------------------------------------------------*
 | destructor                                                                         henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
COMBUST::FlameFront::~FlameFront()
{
  return;
}

/*------------------------------------------------------------------------------------------------*
 | flame front control routine                                                        henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::ProcessFlameFront()
{
  /* This function is accessible from outside the FlameFront class. It can be called e.g. by the
   * interface handle constructor. If the FlameFront class will always do the same thing, this
   * function could be spared. Then the content of this function could be added to the constructor
   * of this class. henke/30.10.08
   * 
   * order of calls:
   * 
   * RefineFlameFront()
   * CaptureFlameFront()
   *
   */
  return;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*
 * class section: refinement                                                                      *
 *------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------------------*
 | refine the region around flame front                                               henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::RefineFlameFront() // input. searchtree from interfacehandle
{
  return;
  /*
     initialize the searchtree
     start a recursive loop to refine the region around the flame front until the maximal number 
       of refinements is reached
       {
         identify the intersection points (zeros of the level set function (G-equation) on the
           refinement cell edges) -> FindFlameFront()
         if a refinement cell is intersected, it will be refined once more -> RefineCell()
       }
   */
}

/*------------------------------------------------------------------------------------------------*
 | find zeros of level set on refinement cell edges                                   henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::FindFlameFront() // LocateFlameFront() // FindIntersectionPoints()
{
  return;
}

/*------------------------------------------------------------------------------------------------*
 | refine a cell embedded in the octree at the level above                            henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::RefineCell()
{
  return;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*
 * class section: capture of flame front                                                          *
 *------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------------------*
 | capture flame front to provide integration cells                                   henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::CaptureFlameFront()
{
  return;
  /*
     triangulate flame front surface (spanned by the intersection points) in every refinement cell
     create a piecewise linear complex (PLC) in Tetgen format from every refinement cell
     call the Constraint Delaunay Tetrahedrization (CDT) to produce burnt and unburnt integration cells
     transform integration cells from Tetgen format to baci format
     
     order function calls:
     
     TriangulateFlameFront()
     buildPLC()
     CreateIntegrationCells()
     TransformIntegrationCells()
  */
}

/*------------------------------------------------------------------------------------------------*
 | triangulate the intersection surface of a refinement cell                          henke 10/08 | 
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::TriangulateFlameFront() //input intersection points for a refinement cell
{
  return;
}

/*------------------------------------------------------------------------------------------------*
 | build piecewise linear complex (PLC) in Tetgen format,based on Ursulas preparePLC() henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::buildPLC()
{
  return;
}

/*------------------------------------------------------------------------------------------------*
 | calls the CDT to create burnt and unburnt integration cells                        henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::CreateIntegrationCells() // output: integration cells in Tetgen format
{
  return;
}

/*------------------------------------------------------------------------------------------------*
 | this could easlily be integrated into CreateIntegrationCells()                     henke 10/08 |
 *------------------------------------------------------------------------------------------------*/
void COMBUST::FlameFront::TransformIntegrationCells() // output: integration cells in baci format
{
  return;
}

#endif // #ifdef CCADISCRET
