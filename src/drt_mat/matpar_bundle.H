/*======================================================================*/
/*!
\file matpar_bundle.H
\brief

\level 1

<pre>
\maintainer Anna Birzle
            birzle@lnm.mw.tum.de
            089/289 15255
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifndef MATPAR_BUNDLE_H
#define MATPAR_BUNDLE_H

/*----------------------------------------------------------------------*/
/* headers */
#include "../drt_inpar/inpar_material.H"
#include "matpar_material.H"
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Teuchos_RCP.hpp>

/*----------------------------------------------------------------------*/
/* forward declarations */
namespace MAT
{
  class Material;

  namespace PAR
  {
    class Material;
  }
}

/*----------------------------------------------------------------------*/
/* declarations */
namespace MAT
{

  namespace PAR
  {

    /*----------------------------------------------------------------------*/
    /// bundle holds all read-in materials of a #DRT::Problem
    ///
    /// <h4>About</h4>
    /// The bundle provides an interface between unique material IDs and
    /// associated material parameters/data. The material ID is set via
    /// input file. It has to be unique, larger than zero or equal zero.
    /// Material ID and data are hold in #matmap_.
    ///
    /// <h4>Special issues for multi-problem-instance applications</h4>
    /// We have for each DRT::Problem instance an individual material bundle.
    /// However, this fact is not transparanet at the read time of the elements
    /// which are only aware of their material ID. The variable #materialreadfromproblem_
    /// of the material bundle of the 0th DRT::Problem instance make it possible to switch
    /// among different DRT::Problem. (The variable #materialreadfromproblem_ is redundant
    /// in the material bundles of all non-0th DRT::Problem instances.)
    ///
    /// \author bborn
    /// \date 02/09
    class Bundle
    {
    public:

      /// construct
      Bundle();

      /// destructor
      ~Bundle();

      /// insert new par of material ID and its data
      void Insert(
        int matid,  ///< material ID
        Teuchos::RCP<MAT::PAR::Material> mat  ///< (validated) material parameters
        );

      /// check if a material exists to provided ID
      ///
      /// \return Upon failure -1 is returned, otherwise >=0
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
        const int num  ///< request is made for this material ID
        ) const;

      /// return material parameters
      MAT::PAR::Parameter* ParameterById(
        const int num  ///< request is made for this material ID
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
      void SetReadFromProblem(
        const int p  ///< index of DRT::Problem instance to read for
        )
      { materialreadfromproblem_ = p; }

      /// reset problem index to read from, i.e. to index 0
      void ResetReadFromProblem() { materialreadfromproblem_ = 0; }

    private:

      /// the map linking material IDs to input materials
      std::map<int,Teuchos::RCP<MAT::PAR::Material> > matmap_;

      /// the index of problem instance of which material read-in shall be performed
      int materialreadfromproblem_;

    };

  }  // namespace PAR

}  // namespace MAT


#endif  // #ifndef MATPAR_BUNDLE_H
