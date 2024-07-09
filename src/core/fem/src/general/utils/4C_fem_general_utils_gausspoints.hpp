/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of gauss formulas

\level 1


*----------------------------------------------------------------------*/
#ifndef FOUR_C_FEM_GENERAL_UTILS_GAUSSPOINTS_HPP
#define FOUR_C_FEM_GENERAL_UTILS_GAUSSPOINTS_HPP

#include "4C_config.hpp"

#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_fem_geometry_element_coordtrafo.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  /// base class of gauss point collection
  class GaussPoints
  {
   public:
    virtual ~GaussPoints() = default;
    /// number of gauss points
    virtual int num_points() const = 0;

    /// spacial dimension
    virtual int num_dimension() const = 0;

    /// gauss point coordinates
    virtual const double* point(int point) const = 0;

    /// gauss weight
    virtual double weight(int point) const = 0;

    /// debug print
    virtual void print() const = 0;
  };

  /// specific collected gauss points for xfem usage
  class CollectedGaussPoints : public GaussPoints
  {
   public:
    CollectedGaussPoints(int size = 0) { gp_.reserve(size); }

    void append(double x, double y, double z, double w)
    {
      // if ( w > 1e-15 )
      gp_.push_back(Entry(x, y, z, w));
    }

    void append(const Core::LinAlg::Matrix<2, 1>& xi, double w)
    {
      gp_.push_back(Entry(xi(0), xi(1), 0., w));
    }

    void append(const Core::LinAlg::Matrix<3, 1>& xi, double w)
    {
      gp_.push_back(Entry(xi(0), xi(1), xi(2), w));
    }

    void increase_reserved(int size) { gp_.reserve(gp_.size() + size); }

    int num_points() const override { return gp_.size(); }

    int num_dimension() const override { return 3; }

    const double* point(int point) const override { return &gp_[point].data[0]; }

    double weight(int point) const override { return gp_[point].data[3]; }

    void print() const override
    {
      std::cout << " collected gauss points:\n";
      for (int i = 0; i < num_points(); ++i)
      {
        std::cout << "    ";
        for (int j = 0; j < num_dimension(); ++j) std::cout << gp_[i].data[j] << " ";
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

    void append(Teuchos::RCP<GaussPoints> gp) { gp_.push_back(gp); }

    /// number of gauss points
    int num_points() const override
    {
      int numpoints = 0;
      for (std::vector<Teuchos::RCP<GaussPoints>>::const_iterator i = gp_.begin(); i != gp_.end();
           ++i)
      {
        Teuchos::RCP<GaussPoints> gp = *i;
        numpoints += gp->num_points();
      }
      return numpoints;
    }

    /// spacial dimension
    int num_dimension() const override { return gp_[0]->num_dimension(); }

    /// gauss point coordinates
    const double* point(int point) const override
    {
      Teuchos::RCP<GaussPoints> gp = find(point);
      return gp->point(point);
    }

    /// gauss weight
    double weight(int point) const override
    {
      Teuchos::RCP<GaussPoints> gp = find(point);
      return gp->weight(point);
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
    void print() const override
    {
      for (std::vector<Teuchos::RCP<GaussPoints>>::const_iterator i = gp_.begin(); i != gp_.end();
           ++i)
      {
        Teuchos::RCP<GaussPoints> gp = *i;
        gp->print();
      }
    }

   private:
    Teuchos::RCP<GaussPoints> find(int& point) const
    {
      for (std::vector<Teuchos::RCP<GaussPoints>>::const_iterator i = gp_.begin(); i != gp_.end();
           ++i)
      {
        Teuchos::RCP<GaussPoints> gp = *i;
        int numpoints = gp->num_points();
        if (numpoints > point)
        {
          return gp;
        }
        point -= numpoints;
      }
      FOUR_C_THROW("gauss point not available");
    }

    std::vector<Teuchos::RCP<GaussPoints>> gp_;
  };

  /// remember calculated gauss points so we do not need to calculate again
  class GaussPointCache
  {
   public:
    static GaussPointCache& instance();

    Teuchos::RCP<GaussPoints> create(Core::FE::CellType distype, int degree);

   private:
    /// cache of already created gauss rules
    std::map<std::pair<Core::FE::CellType, int>, Teuchos::RCP<GaussPoints>> gp_cache_;
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
      const double* point() const { return gp_->point(point_); }

      /// gauss weight at point
      double weight() const { return gp_->weight(point_); }

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
    GaussIntegration(Core::FE::CellType distype);

    /// construct rule for a given element shape
    GaussIntegration(Core::FE::CellType distype, int degree);

    /// construct with a known set of gauss points
    explicit GaussIntegration(Teuchos::RCP<GaussPoints> gp) : gp_(gp) {}

    void clear() { gp_ = Teuchos::null; }

    iterator begin() { return GaussPointIterator(&*gp_, 0); }

    const_iterator begin() const { return GaussPointIterator(&*gp_, 0); }

    iterator end() { return GaussPointIterator(&*gp_, gp_->num_points()); }

    const_iterator end() const { return GaussPointIterator(&*gp_, gp_->num_points()); }

    /// number of gauss points
    int num_points() const { return gp_->num_points(); }

    /// spacial dimension
    int num_dimension() const { return gp_->num_dimension(); }

    /// gauss point coordinates
    const double* point(int point) const { return gp_->point(point); }

    /// gauss weight
    double weight(int point) const { return gp_->weight(point); }

    /// debug print
    void print() const { gp_->print(); }

    Teuchos::RCP<GaussPoints> points() const { return gp_; }

    void set_points(Teuchos::RCP<GaussPoints> gp) { gp_ = gp; }

    /// Create Gauss integration rule of given degree
    template <Core::FE::CellType distype>
    static Teuchos::RCP<GaussPoints> create_projected(
        const Core::LinAlg::Matrix<Core::FE::dim<distype>, Core::FE::num_nodes<distype>>& xie,
        int degree)
    {
      Teuchos::RCP<GaussPoints> gp = GaussPointCache::instance().create(distype, degree);
      Teuchos::RCP<CollectedGaussPoints> cgp =
          Teuchos::rcp(new CollectedGaussPoints(gp->num_points()));

      GaussIntegration intpoints(gp);

      project_gauss_points_local_to_global<distype>(xie, intpoints, cgp);

      return cgp;
    }

    /// Project the given Gauss points in local coordinate system of the element to its global
    /// coordinate system
    template <Core::FE::CellType distype>
    static void project_gauss_points_local_to_global(
        const Core::LinAlg::Matrix<Core::FE::dim<distype>, Core::FE::num_nodes<distype>>& xie,
        GaussIntegration& intpoints, Teuchos::RCP<CollectedGaussPoints>& cgp)
    {
      const int nsd = Core::FE::dim<distype>;
      const int nen = Core::FE::num_nodes<distype>;

      Core::LinAlg::Matrix<nen, 1> funct;
      Core::LinAlg::Matrix<nsd, nen> deriv;

      Core::LinAlg::Matrix<nsd, nsd> xjm;
      Core::LinAlg::Matrix<nsd, 1> xi;

      for (Core::FE::GaussIntegration::iterator iquad = intpoints.begin(); iquad != intpoints.end();
           ++iquad)
      {
        Core::LinAlg::Matrix<nsd, 1> eta(iquad.point());

        // cell shape functions and their first derivatives
        Core::FE::shape_function<distype>(eta, funct);
        Core::FE::shape_function_deriv1<distype>(eta, deriv);

        // local coordinates of gauss point w.r.to background element
        xi.multiply(xie, funct);

        // get transposed of the jacobian matrix d x / d \xi
        // xjm(i,j) = deriv(i,k)*xyze(j,k)
        xjm.multiply_nt(deriv, xie);

        double det = xjm.determinant();

        cgp->append(xi, iquad.weight() * det);
      }
    }

    /// Project the given Gauss points in global coordinate system of the element to its local
    /// coordinate system
    template <Core::FE::CellType distype>
    static Teuchos::RCP<GaussPoints> project_gauss_points_global_to_local(
        const Core::LinAlg::Matrix<Core::FE::dim<distype>, Core::FE::num_nodes<distype>>& xie,
        GaussIntegration& intpoints, const bool& throw_error = true)
    {
      const int nsd = Core::FE::dim<distype>;
      const int nen = Core::FE::num_nodes<distype>;

      Core::LinAlg::Matrix<nen, 1> funct;
      Core::LinAlg::Matrix<nsd, nen> deriv;

      Core::LinAlg::Matrix<nsd, nsd> xjm;
      Core::LinAlg::Matrix<nsd, 1> xi;

      Teuchos::RCP<Core::FE::CollectedGaussPoints> cgp =
          Teuchos::rcp(new Core::FE::CollectedGaussPoints(intpoints.num_points()));

      for (Core::FE::GaussIntegration::iterator iquad = intpoints.begin(); iquad != intpoints.end();
           ++iquad)
      {
        Core::LinAlg::Matrix<nsd, 1> glo(iquad.point());

        bool insideele = Core::Geo::currentToVolumeElementCoordinates(distype, xie, glo, xi);
        if (not insideele && throw_error)
        {
          FOUR_C_THROW("Given Gauss points not inside the element?");
        }

        // cell shape functions and their first derivatives
        Core::FE::shape_function<distype>(xi, funct);
        Core::FE::shape_function_deriv1<distype>(xi, deriv);

        // get transposed of the jacobian matrix d x / d \xi
        // xjm(i,j) = deriv(i,k)*xyze(j,k)
        xjm.multiply_nt(deriv, xie);

        double det = xjm.determinant();

        cgp->append(xi, iquad.weight() / det);
      }
      return cgp;
    }

   private:
    /// internal collection
    Teuchos::RCP<GaussPoints> gp_;
  };

}  // namespace Core::FE

FOUR_C_NAMESPACE_CLOSE

#endif
