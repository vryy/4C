/*----------------------------------------------------------------------*/
/*!
\file material.H

\brief Interface class for complex materials at Gauss points

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef MATERIAL_H
#define MATERIAL_H

#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

/// MAT: complex materials
namespace MAT
{

/// Interface class for complex materials at Gauss points
/*!
   The Material class provides a general way to store history variables
   in Gauss points arising from complex material behavior. It is an
   empty interface and only used to derive classes depending on the
   material type. Therefore it contains nothing except the virtual
   destructor.


   The idea is to have a vector of pointers to this class in each element
   that has to support complex materials.

   <code>
   vector<RefCountPtr<Material> >    mat_;
   </code>

   Different material classes will (in general) not share a common set
   of functions, so the element needs to dynamic_cast to the right
   material implementation.

   <code>
   MicroMaterial* micromat =
        dynamic_cast<MicroMaterial*>(mat_[gp].get());
   </code>

   <h3>Storage</h3>

   In order to support storage of material data (for restart, post
   processing, ...) this class will have to implement DRT::ParObject. It
   is the elements responsibility to pack and unpack its materials along
   with its other data.

   <h3>Alternatives</h3>

   Simple materials that need to rememder nothing or just very few
   values at the Gauss points can be implemented inside the elements
   themselves. Any vectors or matrices the materials needs can be put
   to the elements data_ container and therefore handled automatically
   when in comes to pack/unpack. However, for complex materials that
   could be used with many elements (element variants) an extra class
   might be more appropiate.

   \author u. kue and l.w.
   \date 06/07
*/
class Material
{
public:

  virtual ~Material();

};

}

#endif
#endif
#endif
