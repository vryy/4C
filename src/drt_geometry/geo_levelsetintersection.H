/*!-----------------------------------------------------------------------------------------------*
\file geo_levelsetintersection.H

\brief Class that provides the wizard for the specific functionality for a mesh cut based on a level set field

<pre>
Maintainer: Benedikt Schott and Magnus Winter
            schott@lnm.mw.tum.de, winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef GEO_LEVELSETINTERSECTION_H
#define GEO_LEVELSETINTERSECTION_H

#include "geo_intersection.H"


namespace DRT
{
  class Discretization;
  class Element;
}

namespace GEO
{
  namespace CUT
  {
    class LevelSetIntersection;
  }

  /// contains the cut
  class CutWizardLevelSet : public CutWizard
  {
  public:

    /*!
    \brief Constructor
    */
    CutWizardLevelSet( const DRT::Discretization & dis );

    /*========================================================================*/
    //! @name Add functionality for elements
    /*========================================================================*/

    /*!
    \brief Add this background mesh element to the cut libraries
    */
    void AddElement( DRT::Element * ele, std::vector<double> myphinp, bool include_inner = true );

    /*========================================================================*/
    //! @name Cut functionality routines
    /*========================================================================*/

    /*!
    \brief Cut routine for standard non-parallel framework which is used only for cut_test (level set specific)
    */
    void Cut(
        bool include_inner,
        INPAR::CUT::VCellGaussPts VCellgausstype,  //!< Gauss point generation method for Volumecell
        INPAR::CUT::BCellGaussPts BCellgausstype,  //!< Gauss point generation method for Boundarycell
        bool screenoutput=true
    );

    /*========================================================================*/
    //! @name Get routines
    /*========================================================================*/

    /*!
    \brief Get the level set intersection object
     */
    GEO::CUT::LevelSetIntersection & Mesh() {return *levelsetintersection_; }

  protected:

    void FindNodePositions();


  private:


    /*========================================================================*/
    //! @name private class variables
    /*========================================================================*/

    Teuchos::RCP<GEO::CUT::LevelSetIntersection> levelsetintersection_;  ///< level set intersection object
  };
}

#endif
