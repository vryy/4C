/*----------------------------------------------------------------------*/
/*!
\file sutherland_condif.H

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef SUTHERLAND_CONDIF_H
#define SUTHERLAND_CONDIF_H

#ifdef CCADISCRET

#include "material.H"

namespace MAT
{

  /// temperature-dependent viscosity according to Sutherland law
  class SutherlandCondif : public Material
  {
  public:

    /// construct empty material object
    SutherlandCondif();

    /// construct the material object given material parameters
    explicit SutherlandCondif(MATERIAL* matdata);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */

    virtual int UniqueParObjectId() const { return ParObject_SutherlandCondif; }

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique parobject id delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.

      \param data (in/out): char vector to store class information
    */
    virtual void Pack(vector<char>& data) const;

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
    virtual void Unpack(const vector<char>& data);

    //@}

    /// material type 
    virtual MATERIAL_TYP MaterialType() const { return m_sutherland_condif; }

    /// return material parameters for element calculation

    /*!
      To smooth transition of fluid elements to the unified class
      based material structure we simply return the globally defined
      material parameters. The actual stress calculations are still
      done inside the elements. This will change once we have to
      support different materials in fluid as well.
     */
    MATERIAL* MaterialData() const { return matdata_; }


  private:

    /// my material parameters
    MATERIAL* matdata_;
  };

}

#endif

#endif
