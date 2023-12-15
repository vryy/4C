/*----------------------------------------------------------------------------*/
/*! \file

\brief Converter to use ADAPTER::Coupling type objects in both coupling directions

\level 1

*/
/*----------------------------------------------------------------------------*/

#ifndef BACI_COUPLING_ADAPTER_CONVERTER_H
#define BACI_COUPLING_ADAPTER_CONVERTER_H

/*----------------------------------------------------------------------------*/
/* headers */
#include <map>

class Epetra_Map;
class Epetra_Vector;

namespace Teuchos
{
  template <typename T>
  class RCP;
}

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/* forward declarations */
namespace CORE::ADAPTER
{
  class Coupling;
}

/*----------------------------------------------------------------------------*/
/* definition of classes */
namespace CORE::ADAPTER
{
  /*! \class CouplingConverter
   *  \brief Abstract converter base for master/slave conversion of data
   *
   *  The point is that many generic coupling algorithms that transfer data
   *  between master and slave might be used in both directions. These
   *  algorithms can utilize a Converter to enable use in both directions.
   */
  class CouplingConverter
  {
   public:
    virtual ~CouplingConverter() = default;
    virtual Teuchos::RCP<Epetra_Vector> SrcToDst(Teuchos::RCP<const Epetra_Vector> s) const = 0;

    virtual Teuchos::RCP<Epetra_Vector> DstToSrc(Teuchos::RCP<const Epetra_Vector> d) const = 0;

    virtual Teuchos::RCP<const Epetra_Map> SrcMap() const = 0;

    virtual Teuchos::RCP<const Epetra_Map> DstMap() const = 0;

    virtual Teuchos::RCP<const Epetra_Map> PermSrcMap() const = 0;

    virtual Teuchos::RCP<const Epetra_Map> PermDstMap() const = 0;

    virtual void FillSrcToDstMap(std::map<int, int>& rowmap) const = 0;
  };

  /// master to slave converter
  class CouplingMasterConverter : public CouplingConverter
  {
   public:
    explicit CouplingMasterConverter(const Coupling& coup) : coup_(coup) {}

    Teuchos::RCP<Epetra_Vector> SrcToDst(Teuchos::RCP<const Epetra_Vector> s) const override;

    Teuchos::RCP<Epetra_Vector> DstToSrc(Teuchos::RCP<const Epetra_Vector> d) const override;

    Teuchos::RCP<const Epetra_Map> SrcMap() const override;

    Teuchos::RCP<const Epetra_Map> DstMap() const override;

    Teuchos::RCP<const Epetra_Map> PermSrcMap() const override;

    Teuchos::RCP<const Epetra_Map> PermDstMap() const override;

    void FillSrcToDstMap(std::map<int, int>& rowmap) const override;

   private:
    const Coupling& coup_;
  };

  /// slave to master converter
  class CouplingSlaveConverter : public CouplingConverter
  {
   public:
    explicit CouplingSlaveConverter(const Coupling& coup) : coup_(coup) {}

    Teuchos::RCP<Epetra_Vector> SrcToDst(Teuchos::RCP<const Epetra_Vector> s) const override;

    Teuchos::RCP<Epetra_Vector> DstToSrc(Teuchos::RCP<const Epetra_Vector> d) const override;

    Teuchos::RCP<const Epetra_Map> SrcMap() const override;

    Teuchos::RCP<const Epetra_Map> DstMap() const override;

    Teuchos::RCP<const Epetra_Map> PermSrcMap() const override;

    Teuchos::RCP<const Epetra_Map> PermDstMap() const override;

    void FillSrcToDstMap(std::map<int, int>& rowmap) const override;

   private:
    const Coupling& coup_;
  };

}  // namespace CORE::ADAPTER

BACI_NAMESPACE_CLOSE

#endif  // COUPLING_ADAPTER_CONVERTER_H
