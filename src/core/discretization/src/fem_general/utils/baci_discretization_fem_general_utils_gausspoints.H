/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of gauss formulas

\level 1


*----------------------------------------------------------------------*/
#ifndef BACI_DISCRETIZATION_FEM_GENERAL_UTILS_GAUSSPOINTS_H
#define BACI_DISCRETIZATION_FEM_GENERAL_UTILS_GAUSSPOINTS_H

#include "baci_config.H"

#include "baci_discretization_fem_general_utils_fem_shapefunctions.H"
#include "baci_discretization_geometry_element_coordtrafo.H"

BACI_NAMESPACE_OPEN

namespace CORE::FE
{
  /// base class of gauss point collection
  class GaussPoints
  {
   public:
    virtual ~GaussPoints() = default;
    /// number of gauss points
    virtual int NumPoints() const = 0;

    /// spacial dimension
    virtual int NumDimension() const = 0;

    /// gauss point coordinates
    virtual const double* Point(int point) const = 0;

    /// gauss weight
    virtual double Weight(int point) const = 0;

    /// debug print
    virtual void Print() const = 0;
  };

  /// specific collected gauss points for xfem usage
  class CollectedGaussPoints : public GaussPoints
  {
   public:
    CollectedGaussPoints(int size = 0) { gp_.reserve(size); }

    void Append(double x, double y, double z, double w)
    {
      // if ( w > 1e-15 )
      gp_.push_back(Entry(x, y, z, w));
    }

    void Append(const CORE::LINALG::Matrix<2, 1>& xi, double w)
    {
      gp_.push_back(Entry(xi(0), xi(1), 0., w));
    }

    void Append(const CORE::LINALG::Matrix<3, 1>& xi, double w)
    {
      gp_.push_back(Entry(xi(0), xi(1), xi(2), w));
    }

    void IncreaseReserved(int size) { gp_.reserve(gp_.size() + size); }

    int NumPoints() const override { return gp_.size(); }

    int NumDimension() const override { return 3; }

    const double* Point(int point) const override { return &gp_[point].data[0]; }

    double Weight(int point) const override { return gp_[point].data[3]; }

    void Print() const override
    {
      std::cout << " collected gauss points:\n";
      for (int i = 0; i < NumPoints(); ++i)
      {
        std::cout << "    ";
        for (int j = 0; j < NumDimension(); ++j) std::cout << gp_[i].data[j] << " ";
        std::cout << gp_[i].data[3] << "\n";
      }
    }

   private:
    struct Entry
    {
      Entry(double x, double y, double z, double w)
      {
        data[0] = x;
        data[1] = y;
        data[2] = z;
        data[3] = w;
      }
      double data[4];
    };
    std::vector<Entry> gp_;
  };

  /// specific collected gauss points for xfem usage
  ///
  /// This is meant to be an inverted collection, where the first set of gauss
  /// points integrates a whole element and all following sets substract from
  /// the element.
  class GaussPointsComposite : public GaussPoints
  {
   public:
    explicit GaussPointsComposite(int size) { gp_.reserve(size); }

    void Append(Teuchos::RCP<GaussPoints> gp) { gp_.push_back(gp); }

    /// number of gauss points
    int NumPoints() const override
    {
      int numpoints = 0;
      for (std::vector<Teuchos::RCP<GaussPoints>>::const_iterator i = gp_.begin(); i != gp_.end();
           ++i)
      {
        Teuchos::RCP<GaussPoints> gp = *i;
        numpoints += gp->NumPoints();
      }
      return numpoints;
    }

    /// spacial dimension
    int NumDimension() const override { return gp_[0]->NumDimension(); }

    /// gauss point coordinates
    const double* Point(int point) const override
    {
      Teuchos::RCP<GaussPoints> gp = Find(point);
      return gp->Point(point);
    }

    /// gauss weight
    double Weight(int point) const override
    {
      Teuchos::RCP<GaussPoints> gp = Find(point);
      return gp->Weight(point);
      //     if ( gp_[0]==gp )
      //     {
      //       return gp->Weight( point );
      //     }
      //     else
      //     {
      //       return -1. * gp->Weight( point );
      //     }
    }

    /// debug print
    void Print() const override
    {
      for (std::vector<Teuchos::RCP<GaussPoints>>::const_iterator i = gp_.begin(); i != gp_.end();
           ++i)
      {
        Teuchos::RCP<GaussPoints> gp = *i;
        gp->Print();
      }
    }

   private:
    Teuchos::RCP<GaussPoints> Find(int& point) const
    {
      for (std::vector<Teuchos::RCP<GaussPoints>>::const_iterator i = gp_.begin(); i != gp_.end();
           ++i)
      {
        Teuchos::RCP<GaussPoints> gp = *i;
        int numpoints = gp->NumPoints();
        if (numpoints > point)
        {
          return gp;
        }
        point -= numpoints;
      }
      throw std::runtime_error("gauss point not available");
    }

    std::vector<Teuchos::RCP<GaussPoints>> gp_;
  };

  /// remember calculated gauss points so we do not need to calculate again
  class GaussPointCache
  {
   public:
    static GaussPointCache& Instance();

    Teuchos::RCP<GaussPoints> Create(CORE::FE::CellType distype, int degree);

   private:
    /// cache of already created gauss rules
    std::map<std::pair<CORE::FE::CellType, int>, Teuchos::RCP<GaussPoints>> gp_cache_;
  };

  /// gauss integration interface
  class GaussIntegration
  {
   public:
    /// Very simple internal gauss point iterator.
    ///
    /// With this iterator our gauss point loop looks familiar. Furthermore, we
    /// store a raw pointer and avoid the Teuchos::rcp indirection. This is the gauss
    /// loop, after all!
    class GaussPointIterator
    {
     public:
      GaussPointIterator(GaussPoints* gp, int point) : gp_(gp), point_(point) {}

      /// increment iterator
      void operator++() { point_ += 1; }

      /// point coordinates
      const double* Point() const { return gp_->Point(point_); }

      /// gauss weight at point
      double Weight() const { return gp_->Weight(point_); }

      /// actual point we are at
      int operator*() const { return point_; }

      /// actual point we are at
      int operator->() const { return point_; }

      /// compare
      bool operator==(const GaussPointIterator& other) const
      {
        return gp_ == other.gp_ and point_ == other.point_;
      }

      /// compare
      bool operator!=(const GaussPointIterator& other) const
      {
        return gp_ != other.gp_ or point_ != other.point_;
      }

     private:
      GaussPoints* gp_;
      int point_;
    };

    typedef GaussPointIterator iterator;
    typedef GaussPointIterator const_iterator;

    /// construct the optimal (normal) rule for a given element shape
    GaussIntegration(CORE::FE::CellType distype);

    /// construct rule for a given element shape
    GaussIntegration(CORE::FE::CellType distype, int degree);

    /// construct with a known set of gauss points
    explicit GaussIntegration(Teuchos::RCP<GaussPoints> gp) : gp_(gp) {}

    void Clear() { gp_ = Teuchos::null; }

    iterator begin() { return GaussPointIterator(&*gp_, 0); }

    const_iterator begin() const { return GaussPointIterator(&*gp_, 0); }

    iterator end() { return GaussPointIterator(&*gp_, gp_->NumPoints()); }

    const_iterator end() const { return GaussPointIterator(&*gp_, gp_->NumPoints()); }

    /// number of gauss points
    int NumPoints() const { return gp_->NumPoints(); }

    /// spacial dimension
    int NumDimension() const { return gp_->NumDimension(); }

    /// gauss point coordinates
    const double* Point(int point) const { return gp_->Point(point); }

    /// gauss weight
    double Weight(int point) const { return gp_->Weight(point); }

    /// debug print
    void Print() const { gp_->Print(); }

    Teuchos::RCP<GaussPoints> Points() const { return gp_; }

    void SetPoints(Teuchos::RCP<GaussPoints> gp) { gp_ = gp; }

    /// Create Gauss integration rule of given degree
    template <CORE::FE::CellType distype>
    static Teuchos::RCP<GaussPoints> CreateProjected(
        const CORE::LINALG::Matrix<CORE::FE::dim<distype>, CORE::FE::num_nodes<distype>>& xie,
        int degree)
    {
      Teuchos::RCP<GaussPoints> gp = GaussPointCache::Instance().Create(distype, degree);
      Teuchos::RCP<CollectedGaussPoints> cgp =
          Teuchos::rcp(new CollectedGaussPoints(gp->NumPoints()));

      GaussIntegration intpoints(gp);

      ProjectGaussPointsLocalToGlobal<distype>(xie, intpoints, cgp);

      return cgp;
    }

    /// Project the given Gauss points in local coordinate system of the element to its global
    /// coordinate system
    template <CORE::FE::CellType distype>
    static void ProjectGaussPointsLocalToGlobal(
        const CORE::LINALG::Matrix<CORE::FE::dim<distype>, CORE::FE::num_nodes<distype>>& xie,
        GaussIntegration& intpoints, Teuchos::RCP<CollectedGaussPoints>& cgp)
    {
      const int nsd = CORE::FE::dim<distype>;
      const int nen = CORE::FE::num_nodes<distype>;

      CORE::LINALG::Matrix<nen, 1> funct;
      CORE::LINALG::Matrix<nsd, nen> deriv;

      CORE::LINALG::Matrix<nsd, nsd> xjm;
      CORE::LINALG::Matrix<nsd, 1> xi;

      for (CORE::FE::GaussIntegration::iterator iquad = intpoints.begin(); iquad != intpoints.end();
           ++iquad)
      {
        CORE::LINALG::Matrix<nsd, 1> eta(iquad.Point());

        // cell shape functions and their first derivatives
        CORE::FE::shape_function<distype>(eta, funct);
        CORE::FE::shape_function_deriv1<distype>(eta, deriv);

        // local coordinates of gauss point w.r.to background element
        xi.Multiply(xie, funct);

        // get transposed of the jacobian matrix d x / d \xi
        // xjm(i,j) = deriv(i,k)*xyze(j,k)
        xjm.MultiplyNT(deriv, xie);

        double det = xjm.Determinant();

        cgp->Append(xi, iquad.Weight() * det);
      }
    }

    /// Project the given Gauss points in global coordinate system of the element to its local
    /// coordinate system
    template <CORE::FE::CellType distype>
    static Teuchos::RCP<GaussPoints> ProjectGaussPointsGlobalToLocal(
        const CORE::LINALG::Matrix<CORE::FE::dim<distype>, CORE::FE::num_nodes<distype>>& xie,
        GaussIntegration& intpoints, const bool& throw_error = true)
    {
      const int nsd = CORE::FE::dim<distype>;
      const int nen = CORE::FE::num_nodes<distype>;

      CORE::LINALG::Matrix<nen, 1> funct;
      CORE::LINALG::Matrix<nsd, nen> deriv;

      CORE::LINALG::Matrix<nsd, nsd> xjm;
      CORE::LINALG::Matrix<nsd, 1> xi;

      Teuchos::RCP<CORE::FE::CollectedGaussPoints> cgp =
          Teuchos::rcp(new CORE::FE::CollectedGaussPoints(intpoints.NumPoints()));

      for (CORE::FE::GaussIntegration::iterator iquad = intpoints.begin(); iquad != intpoints.end();
           ++iquad)
      {
        CORE::LINALG::Matrix<nsd, 1> glo(iquad.Point());

        bool insideele = CORE::GEO::currentToVolumeElementCoordinates(distype, xie, glo, xi);
        if (not insideele && throw_error)
        {
          dserror("Given Gauss points not inside the element?");
        }

        // cell shape functions and their first derivatives
        CORE::FE::shape_function<distype>(xi, funct);
        CORE::FE::shape_function_deriv1<distype>(xi, deriv);

        // get transposed of the jacobian matrix d x / d \xi
        // xjm(i,j) = deriv(i,k)*xyze(j,k)
        xjm.MultiplyNT(deriv, xie);

        double det = xjm.Determinant();

        cgp->Append(xi, iquad.Weight() / det);
      }
      return cgp;
    }

   private:
    /// internal collection
    Teuchos::RCP<GaussPoints> gp_;
  };

}  // namespace CORE::FE

BACI_NAMESPACE_CLOSE

#endif
