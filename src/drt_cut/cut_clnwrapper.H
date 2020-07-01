/*---------------------------------------------------------------------*/
/*! \file

\brief Class which stores floating points numbers of arbitrary precision (equivalent to double )

\level 3


*----------------------------------------------------------------------*/

#ifndef CUT_CLNWRAPPER_H
#define CUT_CLNWRAPPER_H
#include <stddef.h>
#include "../drt_lib/drt_dserror.H"
#include <boost/unordered_map.hpp>
#include "cut_memory_manager.H"

// whether we want to use CLN at all or not
#define CLN_CALC
#ifdef CLN_CALC
#include <cln/cln.h>
#endif

// starting precision
#define CLN_START_PRECISION 17

// number of decimal points, that after (numerical)  failure of previous calculation
#define CLN_INCREMENT_STEP 10

// limiting maximum cln precision
#define CLN_LIMIT_PREC 50

// limiting number of increasing precison in the cut_kernel.H
#define CLN_LIMIT_ITER 7

// maximum achievable CLN precison value, error is computed with respect to it
#define CLN_REFERENCE_PREC (CLN_START_PRECISION + CLN_LIMIT_ITER * CLN_INCREMENT_STEP)

#ifndef CLN_CALC
namespace cln
{
  // dummy class to make compiler happy
  struct cl_F
  {
    cl_F(const char* istring);
    cl_F() { dserror("Should not be called!"); };
  };
  struct float_format_t
  {
    float_format_t();
  };
}  // namespace cln
#endif


void* mallocwrap(size_t size);
void deallocwrap(size_t size);

namespace GEO
{
  namespace CUT
  {
    /// Wrapper around CLN long floating point type, that gives better conversion operators,
    /// maintains precision across instances, caches converted double values and
    /// supports running in a custom memory manager
    class ClnWrapper
    {
     public:
      ClnWrapper(const cln::cl_F& a) : value_(a)
      {
#ifndef CLN_CALC
        dserror(
            "In order to call this function BACI should be compiled with CLN library and CLN_CALC "
            "flag on!");
#endif
      }
      /// Default constructor
      ClnWrapper() : value_(CachedConvert(0.0, precision_))
      {
#ifndef CLN_CALC
        dserror(
            "In order to call this function BACI should be compiled with CLN library and CLN_CALC "
            "flag on!");
#endif
      }
      /// initialization from the string.
      /// E.g 0.271828182845904523536028747135266249775724709369996e+1_40
      /// to construct 'e' with precision fo 40 decimal points
      ClnWrapper(const char* istring) : value_(istring)
      {
#ifndef CLN_CALC
        dserror(
            "In order to call this function BACI should be compiled with CLN library and CLN_CALC "
            "flag on!");
#endif
      }
      /// Initialization from the constant double
      ClnWrapper(double a) : value_(CachedConvert(a, precision_))
      {
#ifndef CLN_CALC
        dserror(
            "In order to call this function BACI should be compiled with CLN library and CLN_CALC "
            "flag on!");
#endif
      };
      ClnWrapper(double& a)
      {
        dserror("Constructor for not compile time double to cln::cl_F is not allowed");
      }


#ifdef CLN_CALC
      inline ClnWrapper& operator=(const ClnWrapper& other)
      {
        value_ = other.value_;
        return *this;
      }
      inline ClnWrapper& operator+=(const ClnWrapper& other)
      {
        value_ += other.value_;
        return *this;
      }
      inline ClnWrapper& operator-=(const ClnWrapper& other)
      {
        value_ -= other.value_;
        return *this;
      }
      inline ClnWrapper& operator/=(const ClnWrapper& other)
      {
        value_ /= other.value_;
        return *this;
      }
      inline ClnWrapper& operator*=(const ClnWrapper& other)
      {
        value_ *= other.value_;
        return *this;
      }
      inline ClnWrapper& operator=(double other)
      {
        value_ = CachedConvert(other, precision_);
        return *this;
      }
      // note: we MUST not have implicit convertion, since it will lead to unpredictable behaviour
      // lose of precision and potential memory leaks
      explicit inline operator double() const { return cln::double_approx(value_); }
      inline ClnWrapper& operator+=(double other)
      {
        value_ += CachedConvert(other, value_);
        return *this;
      }
      inline ClnWrapper& operator-=(double other)
      {
        value_ -= CachedConvert(other, value_);
        return *this;
      }
      inline ClnWrapper& operator/=(double other)
      {
        value_ /= CachedConvert(other, value_);
        return *this;
      }
      inline ClnWrapper& operator*=(double other)
      {
        value_ *= CachedConvert(other, value_);
        return *this;
      }
      inline ClnWrapper& operator=(double& other)
      {
        dserror("Unexpected convertion between not-constant double and cln::cl_F");
        return *this;
      }
      inline ClnWrapper& operator+=(double& other)
      {
        dserror("Unexpected convertion between not-constant double and cln::cl_F");
        return *this;
      }
      inline ClnWrapper& operator-=(double& other)
      {
        dserror("Unexpected convertion between not-constant double and cln::cl_F");
        return *this;
      }
      inline ClnWrapper& operator/=(double& other)
      {
        dserror("Unexpected convertion between not-constant double and cln::cl_F");
        return *this;
      }
      inline ClnWrapper& operator*=(double& other)
      {
        dserror("Unexpected convertion between not-constant double and cln::cl_F");
        return *this;
      }
      inline ClnWrapper operator-() const { return (-value_); }
#else
      inline ClnWrapper& operator=(const ClnWrapper& other)
      {
        dserror(
            "In order to call this function BACI should be compiled with CLN library and CLN_CALC "
            "flag on!");
        return *this;
      }
      inline ClnWrapper& operator+=(const ClnWrapper& other)
      {
        dserror(
            "In order to call this function BACI should be compiled with CLN library and CLN_CALC "
            "flag on!");
        return *this;
      }
      inline ClnWrapper& operator-=(const ClnWrapper& other)
      {
        dserror(
            "In order to call this function BACI should be compiled with CLN library and CLN_CALC "
            "flag on!");
        return *this;
      }
      inline ClnWrapper& operator/=(const ClnWrapper& other)
      {
        dserror(
            "In order to call this function BACI should be compiled with CLN library and CLN_CALC "
            "flag on!");
        return *this;
      }
      inline ClnWrapper& operator*=(const ClnWrapper& other)
      {
        dserror(
            "In order to call this function BACI should be compiled with CLN library and CLN_CALC "
            "flag on!");
        return *this;
      }
      inline ClnWrapper& operator=(double other)
      {
        dserror(
            "In order to call this function BACI should be compiled with CLN library and CLN_CALC "
            "flag on!");
        return *this;
      }
      inline ClnWrapper& operator+=(double other)
      {
        dserror(
            "In order to call this function BACI should be compiled with CLN library and CLN_CALC "
            "flag on!");
        return *this;
      }
      inline ClnWrapper& operator-=(double other)
      {
        dserror(
            "In order to call this function BACI should be compiled with CLN library and CLN_CALC "
            "flag on!");
        return *this;
      }
      inline ClnWrapper& operator/=(double other)
      {
        dserror(
            "In order to call this function BACI should be compiled with CLN library and CLN_CALC "
            "flag on!");
        return *this;
      }
      inline ClnWrapper& operator*=(double other)
      {
        dserror(
            "In order to call this function BACI should be compiled with CLN library and CLN_CALC "
            "flag on!");
        return *this;
      }
      inline ClnWrapper& operator=(double& other)
      {
        dserror(
            "In order to call this function BACI should be compiled with CLN library and CLN_CALC "
            "flag on!");
        return *this;
      }
      inline ClnWrapper& operator+=(double& other)
      {
        dserror(
            "In order to call this function BACI should be compiled with CLN library and CLN_CALC "
            "flag on!");
        return *this;
      }
      inline ClnWrapper& operator-=(double& other)
      {
        dserror(
            "In order to call this function BACI should be compiled with CLN library and CLN_CALC "
            "flag on!");
        return *this;
      }
      inline ClnWrapper& operator/=(double& other)
      {
        dserror(
            "In order to call this function BACI should be compiled with CLN library and CLN_CALC "
            "flag on!");
        return *this;
      }
      inline ClnWrapper& operator*=(double& other)
      {
        dserror(
            "In order to call this function BACI should be compiled with CLN library and CLN_CALC "
            "flag on!");
        return *this;
      }
      inline ClnWrapper operator-() const { return *this; }
      // needed to trick to compiler, see macros for no CLN version below
      operator bool() const
      {
        dserror(
            "In order to call this function BACI should be compiled with CLN library and CLN_CALC "
            "flag on!");
        return true;
      }
#endif
      const cln::cl_F& Value() const { return value_; }

      static void SetPrecision(int precision)
      {
        if (precision <= 0) dserror("Invalid preciso of %d", precision);
        precision_ = precision;
      }

      static void ResetPrecision() { precision_ = CLN_START_PRECISION; }

      static unsigned int GetPrecision() { return precision_; }

      //      /// Double to cln conversion
      //      template <class Reference>
      //      static cln::cl_F& CachedConvert(double a, Reference ref);

      template <class ReferenceVal>
      static cln::cl_F& CachedConvert(double a, ReferenceVal ref)
      {
#ifdef CLN_CALC
        // If this is true, floating point underflow returns zero instead of throwing an exception.
        cln::cl_inhibit_floating_point_underflow = true;
        static std::ios_base::Init StreamInitializer;
        // id in the cache containers
        int prec_id = (GEO::CUT::ClnWrapper::precision_ - CLN_START_PRECISION) / CLN_INCREMENT_STEP;
        // we create special vector for zeroes, for faster lookup in the table, since they are
        // needed more often
        static std::vector<std::pair<cln::cl_F, bool>> zeros_cache(
            static_cast<int>(double(CLN_REFERENCE_PREC) / double(CLN_INCREMENT_STEP)));
        // for some reason, native conversion from double 0.0 to CLN loses precision, so we convert
        // explicitly
        if (a == 0.0)
        {
          if (zeros_cache[prec_id].second == false)
          {
            // convert from the string representation
            std::stringstream string_buffer;
            int nsize = GEO::CUT::ClnWrapper::precision_;
            string_buffer << nsize;
            std::string clnumstr = "0e+1_" + string_buffer.str();
            zeros_cache[prec_id].first = clnumstr.c_str();
            zeros_cache[prec_id].second = true;
          }
          return zeros_cache[prec_id].first;
        }
        typedef boost::unordered_map<double, cln::cl_F> maptype;
        // initialize look-up vector  with maximum number of different precisions
        static std::vector<maptype> clnval_cache(
            static_cast<int>(double(CLN_REFERENCE_PREC) / double(CLN_INCREMENT_STEP)));
        maptype::iterator it = clnval_cache[prec_id].find(a);
        if (it != clnval_cache[prec_id].end())
        {
          return it->second;
        }
        else
        {
          // check if this happenning during the const memory container, otherwise it will be never
          // freed (until destruction of static variables, in the end of the program ), but probably
          // that is fine
#ifdef CUSTOM_MEMORY_ALLOCATOR

          bool reset_allocator = false;
          if (GEO::CUT::MemorySingleton::getInstance().IsMemoryPool())
          {
#if EXTENDED_CUT_DEBUG_OUTPUT
            std::cout
                << "NOTICE: Creating a static variable during run on the pointpool! This should "
                   "not happen too often otherwise container will be exhausted!"
                << std::endl;
#endif
            // let it be allocated in the global memory
            GEO::CUT::MemorySingleton::getInstance().SwitchState();
            reset_allocator = true;
          }

#endif

          cln::cl_F newval;
          if (a == 0.0)
          {
            dserror("Should not happen at this point");
          }
          else
          {
            newval = cln::cl_float(a, cln::float_format(ref));
          }
          std::pair<maptype::iterator, bool> in_ins =
              clnval_cache[prec_id].insert(std::make_pair(a, newval));

          // continue running on previous memory allocator state
#ifdef CUSTOM_MEMORY_ALLOCATOR

          if (reset_allocator) GEO::CUT::MemorySingleton::getInstance().SwitchState();

#endif

          return (in_ins.first->second);
        }
#else
        dserror(
            "In order to call this function BACI should be compiled with CLN library and CLN_CALC "
            "flag "
            "on!");
#endif
      }
      //
      //
      //
      ////      // initial value of precision_
      ////      precision_ = CLN_START_PRECISION;

      friend ::std::ostream& operator<<(std::ostream& stream, const GEO::CUT::ClnWrapper& a)
      {
        stream << a.Value();
        return stream;
      }

     private:
      // real CLN value
      cln::cl_F value_;
      // precision of value_
      static unsigned int precision_;
    };


#define MAKE_OPERATOR(_ret_type, _operator)                                              \
  inline _ret_type operator _operator(const ClnWrapper& first, const ClnWrapper& second) \
  {                                                                                      \
    return first.Value() _operator second.Value();                                       \
  }                                                                                      \
  inline _ret_type operator _operator(const ClnWrapper& first, double& second)           \
  {                                                                                      \
    dserror("Unexpected convertion between not-constant double and cln::cl_F");          \
    return first.Value() _operator ClnWrapper::CachedConvert(second, first.Value());     \
  }                                                                                      \
  inline _ret_type operator _operator(const ClnWrapper& first, double second)            \
  {                                                                                      \
    return first.Value() _operator ClnWrapper::CachedConvert(second, first.Value());     \
  }                                                                                      \
  inline _ret_type operator _operator(double first, const ClnWrapper& second)            \
  {                                                                                      \
    return ClnWrapper::CachedConvert(first, second.Value()) _operator second.Value();    \
  }                                                                                      \
  inline _ret_type operator _operator(double& first, const ClnWrapper& second)           \
  {                                                                                      \
    dserror("Unexpected convertion between not-constant double and cln::cl_F");          \
    return ClnWrapper::CachedConvert(first, second.Value()) _operator second.Value();    \
  }

    // Generating all the necessary operators

#ifdef CLN_CALC
    MAKE_OPERATOR(ClnWrapper, +);
    MAKE_OPERATOR(ClnWrapper, *);
    MAKE_OPERATOR(ClnWrapper, /);
    MAKE_OPERATOR(ClnWrapper, -);
    MAKE_OPERATOR(bool, ==);
    MAKE_OPERATOR(bool, !=);
    MAKE_OPERATOR(bool, <);
    MAKE_OPERATOR(bool, >);
    MAKE_OPERATOR(bool, <=);
    MAKE_OPERATOR(bool, >=);
#else

// dummies to trick the compile believe, that such operators are implemented
#define MAKE_DUMMY_OPERATOR(_ret_type, _operator)                                               \
  inline _ret_type operator _operator(const ClnWrapper& first, const ClnWrapper& second)        \
  {                                                                                             \
    dserror(                                                                                    \
        "In order to call this function BACI should be compiled with CLN library and CLN_CALC " \
        "flag on!");                                                                            \
    return first;                                                                               \
  }                                                                                             \
  inline _ret_type operator _operator(const ClnWrapper& first, double& second)                  \
  {                                                                                             \
    dserror(                                                                                    \
        "In order to call this function BACI should be compiled with CLN library and CLN_CALC " \
        "flag on!");                                                                            \
    return first;                                                                               \
  }                                                                                             \
  inline _ret_type operator _operator(const ClnWrapper& first, double second)                   \
  {                                                                                             \
    dserror(                                                                                    \
        "In order to call this function BACI should be compiled with CLN library and CLN_CALC " \
        "flag on!");                                                                            \
    return first;                                                                               \
  }                                                                                             \
  inline _ret_type operator _operator(double first, const ClnWrapper& second)                   \
  {                                                                                             \
    dserror(                                                                                    \
        "In order to call this function BACI should be compiled with CLN library and CLN_CALC " \
        "flag on!");                                                                            \
    return second;                                                                              \
  }                                                                                             \
  inline _ret_type operator _operator(double& first, const ClnWrapper& second)                  \
  {                                                                                             \
    dserror(                                                                                    \
        "In order to call this function BACI should be compiled with CLN library and CLN_CALC " \
        "flag on!");                                                                            \
    return second;                                                                              \
  }

    MAKE_DUMMY_OPERATOR(ClnWrapper, +);
    MAKE_DUMMY_OPERATOR(ClnWrapper, *);
    MAKE_DUMMY_OPERATOR(ClnWrapper, /);
    MAKE_DUMMY_OPERATOR(ClnWrapper, -);
    MAKE_DUMMY_OPERATOR(bool, ==);
    MAKE_DUMMY_OPERATOR(bool, !=);
    MAKE_DUMMY_OPERATOR(bool, <);
    MAKE_DUMMY_OPERATOR(bool, >);
    MAKE_DUMMY_OPERATOR(bool, <=);
    MAKE_DUMMY_OPERATOR(bool, >=);

#endif

  }  // namespace CUT

}  // namespace GEO

// Additional functions in the cln namespace to callback from ClnWrapper
namespace cln
{
#ifdef CLN_CALC

  inline GEO::CUT::ClnWrapper fabs(const GEO::CUT::ClnWrapper& a) { return cln::abs(a.Value()); };
  inline cln::float_format_t float_format(const GEO::CUT::ClnWrapper& a)
  {
    return cln::float_format(a.Value());
  }
  inline GEO::CUT::ClnWrapper abs(const GEO::CUT::ClnWrapper& a) { return cln::abs(a.Value()); };

  inline GEO::CUT::ClnWrapper sqrt(const GEO::CUT::ClnWrapper& a) { return cln::sqrt(a.Value()); };
  inline double double_approx(const GEO::CUT::ClnWrapper& a)
  {
    return cln::double_approx(a.Value());
  };
  inline cln::cl_F cl_float(const GEO::CUT::ClnWrapper& a, const cln::float_format_t& format)
  {
    return cln::cl_float(a.Value(), format);
  }
  inline GEO::CUT::ClnWrapper pow(const GEO::CUT::ClnWrapper& base, int exponent)
  {
    dserror("Not implemented!");
    return base;
  }

#else

  inline GEO::CUT::ClnWrapper fabs(const GEO::CUT::ClnWrapper& a)
  {
    dserror(
        "In order to call this function BACI should be compiled with CLN library and CLN_CALC flag "
        "on!");
    return a;
  };
  inline cln::float_format_t float_format(const GEO::CUT::ClnWrapper& a)
  {
    dserror(
        "In order to call this function BACI should be compiled with CLN library and CLN_CALC flag "
        "on!");
    return cln::float_format_t();
  }
  inline GEO::CUT::ClnWrapper abs(const GEO::CUT::ClnWrapper& a)
  {
    dserror(
        "In order to call this function BACI should be compiled with CLN library and CLN_CALC flag "
        "on!");
    return a;
  };
  inline GEO::CUT::ClnWrapper sqrt(const GEO::CUT::ClnWrapper& a)
  {
    dserror(
        "In order to call this function BACI should be compiled with CLN library and CLN_CALC flag "
        "on!");
    return a;
  };
  inline double double_approx(const GEO::CUT::ClnWrapper& a)
  {
    dserror(
        "In order to call this function BACI should be compiled with CLN library and CLN_CALC flag "
        "on!");
    return a;
  };
  inline cln::cl_F cl_float(const GEO::CUT::ClnWrapper& a, const cln::float_format_t& format)
  {
    dserror(
        "In order to call this function BACI should be compiled with CLN library and CLN_CALC flag "
        "on!");
    return a.Value();
  }

  inline cln::cl_F least_negative_float(const cln::float_format_t& format)
  {
    dserror(
        "In order to call this function BACI should be compiled with CLN library and CLN_CALC flag "
        "on!");
    return cln::cl_F();
  }

  inline cln::cl_F least_positive_float(const cln::float_format_t& format)
  {
    dserror(
        "In order to call this function BACI should be compiled with CLN library and CLN_CALC flag "
        "on!");
    return cln::cl_F();
  }

  inline GEO::CUT::ClnWrapper pow(const GEO::CUT::ClnWrapper& base, int exponent)
  {
    dserror(
        "In order to call this function BACI should be compiled with CLN library and CLN_CALC flag "
        "on!");
    return base;
  }

#endif

}  // namespace cln


namespace GEO
{
  namespace CUT
  {
    // the only function that can perform conversion of
    template <class Matrix_cln, class Matrix_double>
    void ConvClnDouble(const Matrix_cln& in, Matrix_double& out)
    {
      for (unsigned int idx = 0; idx < in.Rows() * in.Columns(); ++idx)
      {
        out.A()[idx] = cln::double_approx(in.A()[idx]);
      }
    }

    template <class Matrix_cln, class Matrix_double>
    void ConvDoubleCLN(const Matrix_double& in, Matrix_cln& out, int precision = 20)
    {
      for (unsigned int idx = 0; idx < in.Rows() * in.Columns(); ++idx)
      {
        ClnWrapper clnnum;
        // zeros do not convert properly to CLN (lose of precision)
        if ((in.A()[idx] == 0.0))
        {
          // returning the cached value from the ClnWrapper cln table
          clnnum = 0.0;
        }
        else
          clnnum = cln::cl_float(in.A()[idx], cln::float_format(precision));
        out.A()[idx] = clnnum;
      }
    }

#ifdef CLN_CALC
    // instantation of the static variables
    template cln::cl_F& ClnWrapper::CachedConvert<unsigned int>(double a, unsigned int ref);
    template cln::cl_F& ClnWrapper::CachedConvert<const cln::cl_F&>(double a, const cln::cl_F& ref);
    template cln::cl_F& ClnWrapper::CachedConvert<cln::cl_F>(double a, cln::cl_F ref);
// initial value of precision_
// unsigned int ClnWrapper::precision_ = CLN_START_PRECISION;
#endif  // endif CLN_CALC

  }  // namespace CUT
}  // namespace GEO


// providing overloads for the arithmetical functions
// (actually dissallowed by the standart), but compiler will most
// probably support it
namespace std
{
#ifdef CLN_CALC

  inline GEO::CUT::ClnWrapper abs(const GEO::CUT::ClnWrapper& a) { return cln::abs(a.Value()); };

  inline GEO::CUT::ClnWrapper sqrt(const GEO::CUT::ClnWrapper& a) { return cln::sqrt(a.Value()); };

//  inline GEO::CUT::ClnWrapper pow(const GEO::CUT::ClnWrapper& base, int exponent)
//  {
//    dserror("Not implemented!");
//    return base;
//  }
#else

  inline GEO::CUT::ClnWrapper abs(const GEO::CUT::ClnWrapper& a)
  {
    dserror(
        "In order to call this function BACI should be compiled with CLN library and CLN_CALC flag "
        "on!");
    return a;
  };

  inline GEO::CUT::ClnWrapper sqrt(const GEO::CUT::ClnWrapper& a)
  {
    dserror(
        "In order to call this function BACI should be compiled with CLN library and CLN_CALC flag "
        "on!");
    return a;
  };

    //  inline GEO::CUT::ClnWrapper pow(const GEO::CUT::ClnWrapper& base, int exponent)
    //  {
    //    dserror(
    //        "In order to call this function BACI should be compiled with CLN library and CLN_CALC
    //        flag " "on!");
    //    return base;
    //  }

#endif

}  // namespace std

#ifdef CLN_CALC
inline GEO::CUT::ClnWrapper fabs(const GEO::CUT::ClnWrapper& a) { return cln::abs(a.Value()); };

inline GEO::CUT::ClnWrapper pow(const GEO::CUT::ClnWrapper& base, int exponent)
{
  dserror("Not implemented!");
  return base;
}
#else
inline GEO::CUT::ClnWrapper fabs(const GEO::CUT::ClnWrapper& a)
{
  dserror(
      "In order to call this function BACI should be compiled with CLN library and CLN_CALC flag "
      "on!");
  return a;
};

inline GEO::CUT::ClnWrapper pow(const GEO::CUT::ClnWrapper& base, int exponent)
{
  dserror(
      "In order to call this function BACI should be compiled with CLN library and CLN_CALC flag "
      "on!");
  return base;
}
#endif

#endif /* CUT_CLNWRAPPER_H */
