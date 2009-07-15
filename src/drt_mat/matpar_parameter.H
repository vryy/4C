/*----------------------------------------------------------------------*/
/*!
\file matpar_parameter.H

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifndef MATPAR_PARAMETER_H
#define MATPAR_PARAMETER_H

#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "Teuchos_RefCountPtr.hpp"
#include "../drt_lib/standardtypes_cpp.H"
#include "matpar_material.H"

/*----------------------------------------------------------------------*/
/* forward declarations */
namespace MAT
{
  namespace PAR
  {
    class Material;
  }
}

/*----------------------------------------------------------------------*/
/* declarations */
namespace MAT
{

  /// Namespace of objects holding the material input parameters.
  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// Base object to hold 'quick' access material parameters
    ///
    /// MAT::PAR::Parameters is derived for the various implemented
    /// materials. These provide the 'quick' access to the read-in
    /// material parameters.
    ///
    /// For every read-in material will exist a single instance (of
    /// a derived class) of this object.
    ///
    /// \author bborn
    /// \date 02/09
    class Parameter
    {

    public:

      /// construct the material object given material parameters
      Parameter(
        Teuchos::RCP<const MAT::PAR::Material> matdata  ///< read and validated material data (of 'slow' access)
        );

      /// destructor
      virtual ~Parameter() { ; }

      /// (unique) material ID
      int Id() const { return id_; }

      /// material type
      INPAR::MAT::MaterialType Type() const { return type_; }

      /// material name
      std::string Name() const { return name_; }

    private:

      /// material ID, as defined in input file
      int id_;

      /// material type
      INPAR::MAT::MaterialType type_;

      /// material name
      std::string name_;

    };  // class Parameter

  }  // namespace PAR

}  // namespace DRT

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
#endif  // #ifndef MATPAR_PARAMETER_H
