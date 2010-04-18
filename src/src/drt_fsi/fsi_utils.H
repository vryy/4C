
#ifdef CCADISCRET

#ifndef FSI_UTILS_H
#define FSI_UTILS_H

#include <set>
#include <string>
#include <functional>
#include <iostream>

#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Abstract_Vector.H>
#include <NOX_Abstract_Group.H>
#include <Teuchos_RefCountPtr.hpp>

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_mapextractor.H"


namespace FSI
{
  /// FSI::UTILS: Random stuff that might be helpful when dealing with FSI problems
  namespace UTILS
  {

    /// write interface jacobian to file
    void DumpJacobian(NOX::Epetra::Interface::Required& interface,
                      double alpha,
                      double beta,
                      Teuchos::RefCountPtr<Epetra_Vector> soln,
                      std::string filename);

    /// shift given map to have no overlap with maps in vector
    Teuchos::RCP<const Epetra_Map> ShiftMap(Teuchos::RCP<const Epetra_Map> emap,
                                            const std::vector<Teuchos::RCP<const Epetra_Map> >& vecSpaces);

    /// create ale discretization parallel to the fluid one
    //void CreateAleDiscretization();



    /*!
    \brief implementation of special clone strategy for automatic generation
           of an ALE-mesh from a given fluid discretization
     */
    class AleFluidCloneStrategy
    {
    public:

      /// constructor
      explicit AleFluidCloneStrategy(){}
      /// destructor
      virtual ~AleFluidCloneStrategy() {}

    protected:

      /// determine element type string and whether element is copied or not
      bool DetermineEleType(
          DRT::Element* actele,
          const bool ismyele,
          vector<string>& eletype);

      /// set element-specific data (material etc.)
      void SetElementData(
          RCP<DRT::Element> newele,
          DRT::Element* oldele,
          const int matid,
          const bool nurbsdis);

      /// returns conditions names to be copied (source and target name)
      map<string,string> ConditionsToCopy();

      /// check for correct material
      void CheckMaterialType(const int matid);

    private:

    }; // class AleFluidCloneStrategy

  }
}

#endif
#endif
