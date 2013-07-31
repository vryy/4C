/*======================================================================*/
/*!
\file matpar_material.H

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifndef MATPAR_MATERIAL_H
#define MATPAR_MATERIAL_H

/*----------------------------------------------------------------------*/
/* headers */
#include "Epetra_Comm.h"
#include "Teuchos_RCP.hpp"

#include "../drt_lib/drt_container.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_inpar/inpar_material.H"
#include "../drt_lib/drt_parobjectfactory.H"
#include "matpar_parameter.H"

/*----------------------------------------------------------------------*/
/* foward declarations */
namespace MAT
{
  namespace PAR
  {
    class Parameter;
  }
}

/*----------------------------------------------------------------------*/
/* declarations */
namespace MAT
{

  namespace PAR
  {
    class ParMaterialType : public DRT::ParObjectType
    {
    public:

      std::string Name() const { return "ParMaterialType"; }

      static ParMaterialType & Instance() { return instance_; };

    private:

      static ParMaterialType instance_;
};

    /// Container for read-in materials
    ///
    /// This object stores the validated material parameters as
    /// DRT::Container.
    ///
    /// \author bborn
    /// \date 02/09
    class Material : public DRT::Container
    {
    public:

      /// @name life span
      //@{

      /// standard constructor
      Material(
        const int id,  ///< unique material ID
        const INPAR::MAT::MaterialType type,  ///< type of material
        const std::string name  ///< name of material
        );


      /// Empty Constructor with type condition_none
      Material();


      /// Copy Constructor
      ///
      /// Makes a deep copy of a condition
      Material(const MAT::PAR::Material& old);

      /// Destructor
      virtual ~Material();

      /// Return unique ParObject id
      ///
      /// every class implementing ParObject needs a unique id defined at the
      /// top of this file.
      virtual int UniqueParObjectId() const { return ParMaterialType::Instance().UniqueParObjectId(); }

      /// Pack this class so it can be communicated
      ///
      /// \ref Pack and \ref Unpack are used to communicate this class
      virtual void Pack(DRT::PackBuffer& data) const;

      /// Unpack data from a char vector into this class
      ///
      /// \ref Pack and \ref Unpack are used to communicate this class
      virtual void Unpack(const std::vector<char>& data);

      /// Set pointer to readily allocated 'quick access' material parameters
      ///
      /// This function is called by the material factory MAT::Material::Factory.
      /// To circumvent more than this single major switch of material type to
      /// object, #params_ are allocated externally.
      inline void SetParameter(
        MAT::PAR::Parameter* matparam  ///< the pointer
        )
      {
        params_ = Teuchos::rcp(matparam);
      }

      //@}

      /// @name Query methods
      //@{

      /// Return material id
      inline virtual int Id() const { return id_; }

      /// Return material name
      inline virtual std::string Name() const { return name_; }

      /// Print this Condition (std::ostream << is also implemented for DRT::Condition)
      virtual void Print(std::ostream& os) const;

      /// Return type of condition
      inline virtual INPAR::MAT::MaterialType Type() const { return type_; }

      /// Return communicator
      inline const Teuchos::RCP<Epetra_Comm> Comm() const { return comm_; }

      /// Return quick accessible material parameter data
      ///
      /// These quick access parameters are stored in separate member #params_;
      /// whereas the originally read ones are stored in DRT::Container base
      inline MAT::PAR::Parameter* Parameter() const { return params_.get();}

      //@}

    protected:

      /// don't want = operator
      Material operator = (const Material& old);

      /// Unique ID of this material, no second material of same ID may exist
      int id_;

      /// Type of this condition
      INPAR::MAT::MaterialType type_;

      /// Name
      std::string name_;

      /// A communicator
      Teuchos::RCP<Epetra_Comm> comm_;

      /// Unwrapped material data for 'quick' access
      Teuchos::RCP<MAT::PAR::Parameter> params_;

    };  // class Material

  }  // namespec PAR

}  // namespace MAT


/// out stream operator
std::ostream& operator << (std::ostream& os, const MAT::PAR::Material& node);


#endif  // #ifndef MATPAR_MATERIAL_H
