/*---------------------------------------------------------------------*/
/*!
\brief Utility methods for the contact integration.

\level 2

\maintainer Matthias Mayr
*/
/*---------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_AUG_CONTACT_INTEGRATOR_UTILS_H_
#define SRC_DRT_CONTACT_AUG_CONTACT_INTEGRATOR_UTILS_H_

#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "contact_aug_utils.H"
#include <unordered_map>

namespace MORTAR
{
  class MortarElement;
  class MortarNode;
}  // namespace MORTAR

namespace CONTACT
{
  class CoIntegrator;
  class CoNode;
  namespace INTEGRATOR
  {
    struct ElementNormal;
    class UniqueProjInfo;

    /// type definitions
    typedef GEN_DATA::pairedvector<MORTAR::MortarElement*, CONTACT::INTEGRATOR::UniqueProjInfo>
        UniqueProjInfoPair;
    typedef CONTACT::AUG::Deriv1stMap Deriv1stMap;
    typedef CONTACT::AUG::Deriv1stVecMap Deriv1stVecMap;
    typedef CONTACT::AUG::Deriv2ndMap Deriv2ndMap;
    typedef CONTACT::AUG::Deriv2ndVecMap Deriv2ndVecMap;

    /** return the size of the current stored entries and not the size of
     *  the potential size (a.k.a. capacity)
     *
     *  \author hiermeier */
    template <typename U, typename T0, typename... Ts>
    double SizeOf(const GEN_DATA::pairedvector<U, T0, Ts...>& p_vec)
    {
      return p_vec.size() * (sizeof(U) + sizeof(T0));
    }

    /** return the size of the current stored entries and not the size of
     *  the potential size (a.k.a. capacity)
     *
     *  \author hiermeier */
    template <typename U, typename... Ts>
    double SizeOf(const GEN_DATA::pairedmatrix<U, Ts...>& pp_mat)
    {
      double sum_size = 0.0;
      for (auto& p_vec : pp_mat) sum_size += sizeof(p_vec.first) + SizeOf(p_vec.second);

      return sum_size;
    }

    /** return the size of the current stored entries and not the size of
     *  the potential size (a.k.a. capacity)
     *
     *  \author hiermeier */
    template <typename T>
    double SizeOf(const std::vector<T>& vec)
    {
      double sum_size = 0.0;
      for (auto& v : vec) sum_size += SizeOf(v);

      return sum_size;
    }

    /** return the size of the current stored entries and not the size of
     *  the potential size (a.k.a. capacity)
     *
     *  \author hiermeier */
    template <typename T, typename U>
    double SizeOf(const std::map<T, U>& map)
    {
      double sum_size = 0.0;
      for (auto& m : map) sum_size += sizeof(m.first) + SizeOf(m.second);

      return sum_size;
    }

    /** perform a symmetry test of the given 2-nd order derivative map
     *
     *  \author hiermeier */
    bool CheckSymmetry(const Deriv2ndMap& deriv2nd, const bool output = false,
        const double rel_tol = 1.0e-12, const double abs_tol = 1.0e-12);

    /** perform a symmetry test of the given second order derivatives
     *
     *  \author hiermeier */
    template <typename T>
    bool TestForSymmetry(const T& deriv2nd, const bool output = false)
    {
      std::cout << "Check for symmetry... ";
      const bool symm_check = CheckSymmetry(deriv2nd, output);
      std::cout << (symm_check ? "succeeded!\n" : "failed!\n");

      return symm_check;
    }

    /** \brief perform a symmetry test of the given vector containing second order
     *  derivatives
     *
     *  \author hiermeier */
    template <typename... Ts>
    bool TestForSymmetry(
        const std::vector<GEN_DATA::pairedmatrix<int, Ts...>>& deriv2nd, const bool output = false)
    {
      for (unsigned i = 0; i < deriv2nd.size(); ++i)
      {
        std::cout << "Check vector component " << i << " for symmetry... ";
        const bool symm_check = CheckSymmetry(deriv2nd[i], output);
        std::cout << (symm_check ? "succeeded!\n" : "failed!\n");

        if (not symm_check) return false;
      }
      return true;
    }

    double BuildAveragedNormalAtSlaveNode(
        std::vector<ElementNormal>& adj_ele_normals, MORTAR::MortarNode& slavenode);

    double UnitSlaveElementNormal(const MORTAR::MortarElement& sele,
        const LINALG::Matrix<3, 2>& tau, LINALG::Matrix<3, 1>& unit_normal);

    void Deriv1st_AveragedSlaveNormal(CONTACT::CoNode& cnode,
        const std::vector<ElementNormal>& adj_ele_normals, const double avg_normal_length,
        Deriv1stVecMap& d_nodal_avg_normal);

    void Deriv1st_NonUnitSlaveNormal(
        const double* xi, MORTAR::MortarElement& sele, Deriv1stVecMap& d_non_unit_normal);

    template <DRT::Element::DiscretizationType slavetype>
    void Deriv1st_NonUnitSlaveNormal(
        MORTAR::MortarElement& sele, const double* xi, Deriv1stVecMap& d_non_unit_normal);

    void Deriv1st_UnitSlaveNormal(const DRT::Element::DiscretizationType slavetype,
        const LINALG::Matrix<3, 1>& unit_normal, const double length_n_inv,
        const Deriv1stVecMap& d_non_unit_normal, Deriv1stVecMap& d_unit_normal, const bool reset);

    void Deriv2nd_AveragedSlaveNormal(CONTACT::CoNode& cnode,
        const std::vector<ElementNormal>& adj_ele_normals, const double avg_normal_length,
        const Deriv1stVecMap& d_nodal_avg_normal);

    void Deriv2nd_NonUnitSlaveNormal(
        const double* xi, MORTAR::MortarElement& sele, Deriv2ndVecMap& dd_non_unit_normal);

    template <DRT::Element::DiscretizationType slavetype>
    void Deriv2nd_NonUnitSlaveNormal(
        MORTAR::MortarElement& sele, const double* xi, Deriv2ndVecMap& dd_non_unit_normal);

    void Deriv2nd_UnitSlaveNormal(const DRT::Element::DiscretizationType slavetype,
        const LINALG::Matrix<3, 1>& unit_normal, const double length_n_inv,
        const Deriv1stVecMap& d_non_unit_normal, const Deriv1stVecMap& d_unit_normal,
        const Deriv2ndVecMap& dd_non_unit_normal, Deriv2ndVecMap& dd_unit_normal);

    /// collect all element nodal dofs and store them in one matrix
    template <unsigned probdim, unsigned numnode>
    void GetElementNodalDofs(
        const MORTAR::MortarElement& ele, LINALG::TMatrix<int, probdim, numnode>& nodal_dofs);

    /// evaluate the Levi Civita pseudo tensor
    double LeviCivitaSymbol(const int i_1, const int i_2, const int i_3);

    /** \brief Find a feasible master element in a given set of master elements
     *
     *  The master element as well as the projection information are stored in
     *  the projInfo paired vector.
     *
     *  \param sele (in)         : slave element
     *  \param meles (in)        : set of master elements
     *  \param boundary_ele (in) : Is the given slave element a boundary element?
     *  \param wrapper (in)      : call-back to the calling integrator object
     *  \param projInfo (out)    : found feasible projection info and the feasible
     *                             master element
     *
     *  \author hiermeier \date 03/17 */
    bool FindFeasibleMasterElements(MORTAR::MortarElement& sele,
        const std::vector<MORTAR::MortarElement*>& meles, bool boundary_ele, CoIntegrator& wrapper,
        UniqueProjInfoPair& projInfo);

    /** \brief Find a feasible master element in a given set of master elements
     *
     *  The master element as well as the projection information are stored in
     *  the projInfo paired vector. In contrast to the alternative call, this
     *  function does not throw a warning if a non-boundary slave element
     *  has non-projectable GPs.
     *
     *  \param sele (in)         : slave element
     *  \param meles (in)        : set of master elements
     *  \param wrapper (in)      : call-back to the calling integrator object
     *  \param projInfo (out)    : found feasible projection info and the feasible
     *                             master element
     *
     *  \author hiermeier \date 03/17 */
    inline bool FindFeasibleMasterElements(MORTAR::MortarElement& sele,
        const std::vector<MORTAR::MortarElement*>& meles, CoIntegrator& wrapper,
        UniqueProjInfoPair& projInfo)
    {
      return FindFeasibleMasterElements(sele, meles, true, wrapper, projInfo);
    }

    /** \brief Is the given GP inside the bounds of the element of the given type?
     *
     *  \param mxi  (in) : projected GP (solution of the local Newton)
     *  \param type (in) : type of the target element
     *  \param tol  (in) : optional (positive) tolerance for the check
     *
     *  \author hiermeier \date 07/17 */
    bool WithinBounds(
        const double* mxi, const DRT::Element::DiscretizationType type, const double tol = 0.0);

    /*--------------------------------------------------------------------------*/
    /** Container class for the unique projection information
     *
     *  \author hiermeier \date 03/17 */
    class UniqueProjInfo
    {
     public:
      /// empty constructor
      UniqueProjInfo()
          : gaussPoints_(0), uniqueProjAlpha_(0), uniqueMxi_(0), scaling_(0), reserve_size_(0)
      { /* empty */
      }

      /// standard constructor
      UniqueProjInfo(unsigned reserve_size)
          : gaussPoints_(0),
            uniqueProjAlpha_(0),
            uniqueMxi_(0),
            scaling_(0),
            reserve_size_(reserve_size)
      { /* empty */
      }

      /// default copy constructor
      UniqueProjInfo(const UniqueProjInfo& projinfo) = default;

      /** \brief Fill the container with new information
       *
       *  \param[in] gp               Gauss point id number
       *  \param[in] uniqueProjAlpha  auxiliary distance factor
       *  \param[in] uniqueMxi        projected gauss point parametric coordinates
       *  \param[in] scaling          scaling factor if the gp has more than
       *                              one feasible target element
       *
       *  \author hiermeier \date 03/17 */
      void Insert(const int gp, const double uniqueProjAlpha, const double uniqueMxi[],
          const double scaling)
      {
        ReserveSize();

        gaussPoints_.push_back(gp);
        uniqueProjAlpha_.push_back(uniqueProjAlpha);
        uniqueMxi_.push_back(LINALG::Matrix<2, 1>(uniqueMxi, false));
        scaling_.push_back(scaling);
      }

      /** \brief Print the unique projection info class
       *
       *  \param[in/out] stream  output stream object
       *
       *  \author hiermeier \date 03/17 */
      void Print(std::ostream& stream) const
      {
        stream << "--- UniqueProjInfo ---\n";
        stream << "#GuassPoints: " << gaussPoints_.size()
               << " ( capacity: " << gaussPoints_.capacity() << " )\n{ ";
        for (std::vector<int>::const_iterator cit = gaussPoints_.begin(); cit != gaussPoints_.end();
             ++cit)
          stream << *cit << " ";
        stream << "}\n\n";

        stream << "#ProjAlpha: " << uniqueProjAlpha_.size()
               << " ( capacity: " << uniqueProjAlpha_.capacity() << " )\n{ ";
        for (std::vector<double>::const_iterator cit = uniqueProjAlpha_.begin();
             cit != uniqueProjAlpha_.end(); ++cit)
          stream << *cit << " ";
        stream << "}\n\n";

        stream << "#Mxi: " << uniqueMxi_.size() << " ( capacity: " << uniqueMxi_.capacity()
               << " )\n{ ";
        for (std::vector<LINALG::Matrix<2, 1>>::const_iterator cit = uniqueMxi_.begin();
             cit != uniqueMxi_.end(); ++cit)
          stream << " [ " << (*cit)(0) << ", " << (*cit)(1) << " ] ";
        stream << " }\n\n";

        stream << "#GP-Scaling: " << scaling_.size() << " ( capacity: " << scaling_.capacity()
               << " )\n{ ";
        for (std::vector<double>::const_iterator cit = scaling_.begin(); cit != scaling_.end();
             ++cit)
          stream << *cit << " ";
        stream << " }\n\n";
        stream << std::flush;
      }

      /// Gauss point id numbers
      std::vector<int> gaussPoints_;

      /// auxiliary distance values
      std::vector<double> uniqueProjAlpha_;

      /// projected master parametric coordinates
      std::vector<LINALG::Matrix<2, 1>> uniqueMxi_;

      /** Scale the GP weight if the slave gp projects onto more than one
       *  master element, otherwise the scaling factor is equal to 1.0 */
      std::vector<double> scaling_;

     private:
      /// reserve capacity for the member variables
      inline void ReserveSize()
      {
        if (gaussPoints_.capacity() > 0) return;

        gaussPoints_.reserve(reserve_size_);
        uniqueProjAlpha_.reserve(reserve_size_);
        uniqueMxi_.reserve(reserve_size_);
        scaling_.reserve(reserve_size_);
      }

      int reserve_size_;
    };

    /*--------------------------------------------------------------------------*/
    /// \brief container for an element normal at the parametric position xi
    struct ElementNormal
    {
      /// constructor
      ElementNormal() : xi_{0.0, 0.0}, unit_n_(true), ele_(NULL), length_n_inv_(-1.0){/* empty */};

      /// destructor
      ~ElementNormal(){/* empty */};

      /// access the element unit normal components (read-only)
      double operator()(unsigned i) const { return unit_n_(i, 0); }

      /// access the element unit normal components
      double& operator()(unsigned i) { return unit_n_(i, 0); }

      /// corresponding parametric coordinates
      double xi_[2];

      /// unit normal at position xi_
      LINALG::Matrix<3, 1> unit_n_;

      /// pointer to the corresponding element
      MORTAR::MortarElement* ele_;

      /// inverse length of the non-unit normal
      double length_n_inv_;
    };

  }  // namespace INTEGRATOR
}  // namespace CONTACT


#endif /* SRC_DRT_CONTACT_AUG_CONTACT_INTEGRATOR_UTILS_H_ */
