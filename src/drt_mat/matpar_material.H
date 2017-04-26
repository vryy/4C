/*!-----------------------------------------------------------------------
\file matpar_material.H

\level 1

<pre>
\maintainer Anna Birzle
            birzle@lnm.mw.tum.de
            089/289 15255
</pre>

*-----------------------------------------------------------------------*/

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
