/*----------------------------------------------------------------------*/
/*!
\brief structure-specific utils and auxiliary functions

\level 1

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------*/


#ifndef STRUCTURE_UTILS_MAPEXTRACTOR_H
#define STRUCTURE_UTILS_MAPEXTRACTOR_H

#include "../linalg/linalg_mapextractor.H"
#include "../drt_inpar/inpar_structure.H"

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace STR
{
  namespace AUX
  {
    /// Determine norm of force residual
    double CalculateVectorNorm(const enum INPAR::STR::VectorNorm norm,  ///< type of norm to use
        const Teuchos::RCP<Epetra_Vector> vect,                         ///< the vector of interest
        const int numneglect =
            0  ///< number of DOFs that have to be neglected for possible length scaling
    );

    /// specific MultiMapExtractor to handle the structure field
    class MapExtractor : public LINALG::MultiMapExtractor
    {
     public:
      enum
      {
        cond_other = 0,
        cond_fsi = 1,
        cond_lung_asi = 2,
        cond_bio_gr = 3,
        cond_ale_wear = 4,
        cond_fpsi = 5,
        cond_immersed = 6,
        cond_pasi = 7
      };

      /// setup the whole thing
      void Setup(
          const DRT::Discretization& dis, const Epetra_Map& fullmap, bool overlapping = false);

      /// get all element gids those nodes are touched by any condition
      Teuchos::RCP<std::set<int>> ConditionedElementMap(const DRT::Discretization& dis) const;

      MAP_EXTRACTOR_VECTOR_METHODS(Other, cond_other)
      MAP_EXTRACTOR_VECTOR_METHODS(FSICond, cond_fsi)
      MAP_EXTRACTOR_VECTOR_METHODS(LungASICond, cond_lung_asi)
      MAP_EXTRACTOR_VECTOR_METHODS(BioGrCond, cond_bio_gr)
      MAP_EXTRACTOR_VECTOR_METHODS(AleWearCond, cond_ale_wear)
      MAP_EXTRACTOR_VECTOR_METHODS(FPSICond, cond_fpsi)
      MAP_EXTRACTOR_VECTOR_METHODS(IMMERSEDCond, cond_immersed)
      MAP_EXTRACTOR_VECTOR_METHODS(PASICond, cond_pasi)
    };
  }  // namespace AUX
}  // namespace STR

#endif
