/*======================================================================*/
/*!
\file matpar_bundle.H
\brief

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifndef MATPAR_BUNDLE_H
#define MATPAR_BUNDLE_H
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "../drt_inpar/inpar_material.H"
#include "matpar_material.H"
#include "material.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Teuchos_RCP.hpp>


namespace MAT
{
  class Material;

  namespace PAR
  {
    class Material;
  }
}


namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// bundle holds all read-in materials of a #DRT::Problem
    class Bundle
    {
    public:

      /// construct
      Bundle();

      /// destructor
      ~Bundle();

      /// insert new par of material ID and its data
      void Insert(
        int matid,
        Teuchos::RCP<MAT::PAR::Material> mat
        );

      /// check if a material exists to provided ID
      int Find(
        const int id
        ) const;

      /// provide access to material map (a li'l dirty)
      const std::map<int,Teuchos::RCP<MAT::PAR::Material> >* Map() const { return &(matmap_); }

      /// make quick access parameters
      void MakeParameters();

      /// return number of defined materials
      int Num() const { return matmap_.size(); }

      /// return materials by ID
      Teuchos::RCP<MAT::PAR::Material> ById(
        const int num
        ) const;

      /// return material parameters
      MAT::PAR::Parameter* ParameterById(
        const int num
        ) const
      {
        return ById(num)->Parameter();
      }

      /// return (first) ID by material type
      /// 
      /// \return The ID of seached for material type.
      ///         If the search is unsuccessful -1 is returned
      int FirstIdByType(
        const INPAR::MAT::MaterialType type
        ) const;

      /// return problem index to read from
      int GetReadFromProblem() const { return materialreadfromproblem_; }

      /// set problem index to read from
      void SetReadFromProblem(const int p) { materialreadfromproblem_ = p; }

      /// reset problem index to read from, i.e. to index 0
      void ResetReadFromProblem() { materialreadfromproblem_ = 0; }
      
    private:

      /// the map linking material IDs to input materials
      std::map<int,Teuchos::RCP<MAT::PAR::Material> > matmap_;

      /// the problem instance of which read-in is performed
      int materialreadfromproblem_;

    };

  }  // namespace PAR

}  // namespace MAT


#endif  // #ifdef CCADISCRET
#endif  // #ifndef MATPAR_BUNDLE_H
