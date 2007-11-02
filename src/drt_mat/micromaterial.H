/*----------------------------------------------------------------------*/
/*!
\file micromaterial.H

\brief class for handling of micro-macro transitions

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef MICROMATERIAL_H
#define MICROMATERIAL_H

#ifdef CCADISCRET

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <list>
#include <string>
#include <algorithm>

#include <Epetra_Comm.h>
#include <Teuchos_RefCountPtr.hpp>

#include "../io/io_drt.H"
#include "../io/io_drt_micro.H"

#include "../drt_lib/drt_discret.H"
#include "material.H"

#include "../drt_stru_multi/microstrugenalpha.H"

using namespace std;
using namespace Teuchos;
using namespace IO;


namespace MAT {

  // forward
  class MicroMaterialGP;

class MicroMaterial : public Material
{
public:

    /// construct empty material object
    MicroMaterial();

    /// construct the material object given material parameters
    explicit MicroMaterial(MATERIAL* matdata);

    //! @name Packing and Unpacking

    /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
    */
    virtual int UniqueParObjectId() const { return ParObject_MicroMaterial; }

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
    virtual MATERIAL_TYP MaterialType() const { return m_struct_multiscale; }

    void Evaluate(const Epetra_SerialDenseMatrix* defgrd,
                  Epetra_SerialDenseMatrix* cmat,
                  Epetra_SerialDenseVector* stress,
                  double* density,
                  const int gp,
                  const int ele_ID,
                  const double time,
                  const string action);

private:

  vector<RefCountPtr<MicroMaterialGP> > matgp_;
  /// my material parameters
  MATERIAL* matdata_;
};

}

#endif
#endif

