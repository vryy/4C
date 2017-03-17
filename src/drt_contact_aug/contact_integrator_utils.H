/*---------------------------------------------------------------------*/
/*!
\file contact_integrator_utils.H

\brief Utility methods for the contact integration.

\level 2

\maintainer Michael Hiermeier

\date Mar 8, 2017

*/
/*---------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_AUG_CONTACT_INTEGRATOR_UTILS_H_
#define SRC_DRT_CONTACT_AUG_CONTACT_INTEGRATOR_UTILS_H_

#include "../headers/pairedvector.H"
#include "../linalg/linalg_fixedsizematrix.H"

namespace MORTAR {
  class MortarElement;
}

namespace CONTACT {
  class CoIntegrator;
namespace INTEGRATOR {
  class UniqueProjInfo;
  typedef GEN::pairedvector<MORTAR::MortarElement*, CONTACT::INTEGRATOR::UniqueProjInfo > UniqueProjInfoPair;

  /** \brief (Re)set a paired-vector
   *
   *  \param reserve_size (in) : maximal necessary capacity of the paired vector
   *  \param vec          (out): paired-vector which will be reseted
   *
   *  \author hiermeier \date 03/17 */
  template <class T>
  inline void ResetPairedVector(
      const unsigned reserve_size,
      GEN::pairedvector<int,T>& paired_vec )
  {
    paired_vec.clear();
    if ( paired_vec.capacity() < reserve_size )
      paired_vec.resize( reserve_size );
  }

  /** \brief (Re)set a std::vector of paired-vectors
   *
   *  \param reserve_size (in) : maximal necessary capacity of the paired vector
   *  \param vec          (out): vector of paired-vectors which will be reseted
   *
   *  \author hiermeier \date 03/17 */
  template <unsigned num_vec, class T >
  inline void ResetPairedVector(
      const unsigned reserve_size,
      std::vector<GEN::pairedvector<int,T> >& vec )
  {
    if ( vec.size() != num_vec )
    {
      vec.clear();
      vec.resize( num_vec, GEN::pairedvector<int,T>( reserve_size ) );
      return;
    }

    for ( typename std::vector<GEN::pairedvector<int,T> >::iterator it =
        vec.begin(); it != vec.end(); ++it )
    {
      ResetPairedVector<T>( reserve_size, *it );
    }
  }

  /** \brief Find a feasible master element in a given set of master elements
   *
   *  The master element as well as the projection information are stored in
   *  the projInfo paired vector.
   *
   *  \param sele (in)         : slave element
   *  \param meles (in)        : set of master elements
   *  \param boundary_ele (in) : Is the given slave element a boundary element?
   *  \param wrapper (in)      : call-back to the calling integrator class
   *  \param projInfo (out)    : found feasible projection info and the feasible
   *                             master element
   *
   *  \author hiermeier \date 03/17 */
  bool FindFeasibleMasterElement3D(
      MORTAR::MortarElement& sele,
      const std::vector<MORTAR::MortarElement*> & meles,
      bool boundary_ele,
      CoIntegrator& wrapper,
      UniqueProjInfoPair& projInfo );

  /** \brief Find a feasible master element in a given set of master elements
   *
   *  The master element as well as the projection information are stored in
   *  the projInfo paired vector.
   *
   *  \param sele (in)         : slave element
   *  \param meles (in)        : set of master elements
   *  \param wrapper (in)      : call-back to the calling integrator class
   *  \param projInfo (out)    : found feasible projection info and the feasible
   *                             master element
   *
   *  \author hiermeier \date 03/17 */
  bool FindFeasibleMasterElement2D(
      MORTAR::MortarElement& sele,
      const std::vector<MORTAR::MortarElement*> & meles,
      CoIntegrator& wrapper,
      UniqueProjInfoPair& projInfo);

  class UniqueProjInfo
  {
  public:
    /// empty constructor
    UniqueProjInfo()
        : gaussPoints_(0),
          uniqueProjAlpha_(0),
          uniqueMxi_(0),
          reserve_size_(0)
    {

    }

    /// standard constructor
    UniqueProjInfo( int numGP, unsigned msize )
        : gaussPoints_(0),
          uniqueProjAlpha_(0),
          uniqueMxi_(0),
          reserve_size_( std::ceil( numGP / msize ) )
    {

    }

    void Insert( int gp, double uniqueProjAlpha, double uniqueMxi[] )
    {
      ReserveSize( true );

      gaussPoints_.push_back( gp );
      uniqueProjAlpha_.push_back( uniqueProjAlpha );
      uniqueMxi_.push_back( LINALG::Matrix<2,1>( uniqueMxi, false ) );
    }

    void Insert( int gp, double uniqueMxi[] )
    {
      ReserveSize( false );

      gaussPoints_.push_back( gp );
      uniqueMxi_.push_back( LINALG::Matrix<2,1>( uniqueMxi, false ) );
    }

    void Print( std::ostream& stream ) const
    {
      stream << "--- UniqueProjInfo ---\n";
      stream << "#GuassPoints: " << gaussPoints_.size() << " ( capacity: "
          << gaussPoints_.capacity() << " )\n{ ";
      for ( std::vector<int>::const_iterator cit = gaussPoints_.begin();
            cit != gaussPoints_.end(); ++cit )
        stream << *cit << " ";
      stream << "}\n\n";

      stream << "#ProjAlpha: " << uniqueProjAlpha_.size() << " ( capacity: "
             << uniqueProjAlpha_.capacity() << " )\n{ ";
      for ( std::vector<double>::const_iterator cit = uniqueProjAlpha_.begin();
            cit != uniqueProjAlpha_.end(); ++cit )
        stream << *cit << " ";
      stream << "}\n\n";

      stream << "#Mxi: " << uniqueMxi_.size() << " ( capacity: "
          << uniqueMxi_.capacity() << " )\n{ ";
      for ( std::vector<LINALG::Matrix<2,1> >::const_iterator cit =
            uniqueMxi_.begin(); cit != uniqueMxi_.end(); ++cit )
        stream << " [ " << (*cit)(0) << ", " << (*cit)(1) << " ] ";
      stream << " }\n\n";
      stream << std::flush;
    }

    std::vector<int> gaussPoints_;
    std::vector<double> uniqueProjAlpha_;
    std::vector<LINALG::Matrix<2,1> > uniqueMxi_;

  private:
    inline void ReserveSize( bool reserve_alpha )
    {
      if ( gaussPoints_.capacity() > 0 )
        return;

      gaussPoints_.reserve( reserve_size_ );

      if ( reserve_alpha )
        uniqueProjAlpha_.reserve( reserve_size_ );

      uniqueMxi_.reserve( reserve_size_ );
    }

    int reserve_size_;

  };

} // namespace INTEGRATOR
} // namespace CONTACT


#endif /* SRC_DRT_CONTACT_AUG_CONTACT_INTEGRATOR_UTILS_H_ */
