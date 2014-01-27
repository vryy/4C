/*!---------------------------------------------------------------------------
\file drt_meshfree_cell_utils.H

\brief utils for meshfree cells in meshfree discretisations

<pre>
Maintainer: Keijo Nissen
            nissen@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>

*--------------------------------------------------------------------------*/
#ifndef DRT_MESHFREE_CELL_UTILS_H
#define DRT_MESHFREE_CELL_UTILS_H

#include "../drt_lib/drt_element.H"  // for distype enum 'DRT::Element::DiscretizationType'
#include "../drt_lib/drt_dserror.H"  // for constructor in class GetCellShapeFunktionsDerivsWeightsAtGaussPointXi
#include "../linalg/linalg_serialdensematrix.H" // for member in class GetCellShapeFunktionsDerivsWeightsAtGaussPointXi
#include "../linalg/linalg_serialdensevector.H" // for member in class GetCellShapeFunktionsDerivsWeightsAtGaussPointXi
#include "../drt_fem_general/drt_utils_integration.H" // for Gauss rules
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H" // DRT::UTILS::DisTypeToDim

namespace DRT
{

namespace MESHFREE
{
/*==========================================================================*/
//forward declarations:
/*==========================================================================*/
class MeshfreeNode;

/*--------------------------------------------------------------------------*/
/*!
 * \brief Standard Gauss rule for integrating meshfree discretisations over an
 * backgroundmesh. (Paradox, I know...) Don't forget to update struct
 * DisTypeToNumGaussPoints as well, when making any changes!
 *///                                                   (public) nis Feb12
/*--------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
struct DisTypeToOptGaussRule {};
template<> struct DisTypeToOptGaussRule<DRT::Element::hex8 > {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_64point; };
template<> struct DisTypeToOptGaussRule<DRT::Element::tet4 > {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_tet_4point;  };
template<> struct DisTypeToOptGaussRule<DRT::Element::quad4> {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_16point;};
template<> struct DisTypeToOptGaussRule<DRT::Element::tri3 > {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_tri_7point;  };
template<> struct DisTypeToOptGaussRule<DRT::Element::line2> {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule_line_4point; };

/*--------------------------------------------------------------------------*/
/*!
 * \brief Number of Gauss points for standard Gauss rule for integrating
 * meshfree discretisations over an backgroundmesh. To be kept up to date with
 * struct DisTypeToOptGaussRule
 *///                                                   (public) nis Feb12
/*--------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
struct DisTypeToNumGaussPoints {};
template<> struct DisTypeToNumGaussPoints<DRT::Element::hex8 > {static const int ngp = 64;};
template<> struct DisTypeToNumGaussPoints<DRT::Element::tet4 > {static const int ngp =  4;};
template<> struct DisTypeToNumGaussPoints<DRT::Element::quad4> {static const int ngp = 16;};
template<> struct DisTypeToNumGaussPoints<DRT::Element::tri3 > {static const int ngp =  7;};
template<> struct DisTypeToNumGaussPoints<DRT::Element::line2> {static const int ngp =  4;};

//template<DRT::Element::DiscretizationType distype>
//struct DisTypeMinGaussRule {};  //dserror("Distype not implemented for meshfree discretisations.");};
//template<> struct DisTypeMinGaussRule<DRT::Element::hex8>  {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_1point; };
//template<> struct DisTypeMinGaussRule<DRT::Element::tet4>  {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_tet_1point; };
//template<> struct DisTypeMinGaussRule<DRT::Element::quad4> {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_1point;};
//template<> struct DisTypeMinGaussRule<DRT::Element::tri3>  {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_tri_1point; };
//template<> struct DisTypeMinGaussRule<DRT::Element::line2> {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule_line_1point;};

/*--------------------------------------------------------------------------*/
/*!
 * \brief Provides interface to CellGaussPoints
 *///                                                   (public) nis Feb12
/*--------------------------------------------------------------------------*/
class CellGaussPointInterface : public DRT::SingletonDestruction
{
public:
  /*------------------------------------------------------------------------*/
  /*! \brief provides access to derived class
   *///                                                   (public) nis Mar12
  /*------------------------------------------------------------------------*/
  static CellGaussPointInterface* Impl(DRT::Element::DiscretizationType distype);

  /*------------------------------------------------------------------------*/
  /*! \brief map local Gauss point to global real space - dummy
   *///                                                   (public) nis Mar12
  /*------------------------------------------------------------------------*/
  virtual int GetCellGaussPointsAtX(
    LINALG::SerialDenseMatrix &  Xa_sdm, //!< (in)  position of points
    LINALG::SerialDenseMatrix &  X_sdm,  //!< (out) xyz-values of Gauss points
    LINALG::SerialDenseVector &  w_sdm   //!< (out) weights of Gauss points in xyz-system
    ) const = 0;
};


/*--------------------------------------------------------------------------*/
/*!
 * \brief Gives global positions of Gauss points for standard Gauss rule for
 * integrating meshfree discretisations over an backgroundmesh.
 *
 *///                                                   (public) nis Feb12
/*--------------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
class CellGaussPoints : public CellGaussPointInterface
{
public:

  /*------------------------------------------------------------------------*/
  /*! \brief Singleton access method
   *///                                                   (public) nis Mar12
  /*------------------------------------------------------------------------*/
  static CellGaussPoints<distype> * Instance(bool create);

  /*------------------------------------------------------------------------*/
  /*! \brief Called upon singleton destruction
   *///                                                   (public) nis Mar12
  /*------------------------------------------------------------------------*/
  virtual void Done();

  /*------------------------------------------------------------------------*/
  /*! \brief map local Gauss point to global real space
   *///                                                   (public) nis Mar12
  /*------------------------------------------------------------------------*/
  virtual int GetCellGaussPointsAtX(
    LINALG::SerialDenseMatrix &  Xa_sdm, //!< (in)  position of points
    LINALG::SerialDenseMatrix &  X_sdm,  //!< (out) xyz-values of Gauss points
    LINALG::SerialDenseVector &  w_sdm   //!< (out) weights of Gauss points in xyz-system
    ) const;

private:
  /*------------------------------------------------------------------------*/
  /*!
   * \brief private constructor
   *
   * to be called only once and only by Instance(). Add more distypes ad lib
   * but with care!
   *///                                                  (public) nis Feb12
  /*------------------------------------------------------------------------*/
  CellGaussPoints();

  //! number of space dimensions
  static const int nsd_ = DRT::UTILS::DisTypeToDim<distype>::dim;
  //! number of element points = number of shape functions
  static const int nek_ = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
  //! number of Gauss points
  static const int ngp_ = DisTypeToNumGaussPoints<distype>::ngp;

  //! matrix of cell shape function values at Gauss points
  static LINALG::SerialDenseMatrix N_;
  //! vector of weights of Gauss points
  static LINALG::SerialDenseVector w_;
  //! vector of matrices of values of cell shape function derivatives at Gauss points
  static std::vector<LINALG::SerialDenseMatrix> dN_;
  //! flag indicating instantiation status of CellGaussPoints<distype> - initilized to false
  static bool is_instantiated_;
};

} // end namespace MESHFREE
} // end namespace DRT

#endif // DRT_MESHFREE_CELL_UTILS_H
