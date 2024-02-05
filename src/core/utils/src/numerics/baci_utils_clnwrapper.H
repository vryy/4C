/*---------------------------------------------------------------------*/
/*! \file

\brief Class which stores floating points numbers of arbitrary precision (equivalent to double )

\level 3


*----------------------------------------------------------------------*/

#ifndef BACI_UTILS_CLNWRAPPER_H
#define BACI_UTILS_CLNWRAPPER_H
#include "baci_config.H"

#include "baci_utils_exceptions.H"

#include <cln/cln.h>

#include <cstddef>
#include <sstream>
#include <unordered_map>
#include <vector>

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

BACI_NAMESPACE_OPEN

namespace CORE::CLN
{
  /// Wrapper around CLN long floating point type, that gives better conversion operators,
  /// maintains precision across instances, caches converted double values and
  /// supports running in a custom memory manager
  class ClnWrapper
  {
   public:
    ClnWrapper(const cln::cl_F& a) : value_(a) {}
    /// Default constructor
    ClnWrapper() : value_(CachedConvert(0.0, precision_)) {}
    /// initialization from the string.
    /// E.g 0.271828182845904523536028747135266249775724709369996e+1_40
    /// to construct 'e' with precision fo 40 decimal points
    ClnWrapper(const char* istring) : value_(istring) {}
    /// Initialization from the constant double
    ClnWrapper(double a) : value_(CachedConvert(a, precision_)){};
    ClnWrapper(double& a)
    {
      dserror("Constructor for not compile time double to cln::cl_F is not allowed");
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

    const cln::cl_F& Value() const { return value_; }

    static void SetPrecision(int precision)
    {
      if (precision <= 0) dserror("Invalid preciso of %d", precision);
      precision_ = precision;
    }

    static void ResetPrecision() { precision_ = CLN_START_PRECISION; }

    static unsigned int GetPrecision() { return precision_; }

    template <class ReferenceVal>
    static cln::cl_F& CachedConvert(double a, ReferenceVal ref)
    {
      // If this is true, floating point underflow returns zero instead of throwing an exception.
      cln::cl_inhibit_floating_point_underflow = true;
      static std::ios_base::Init StreamInitializer;
      // id in the cache containers
      int prec_id = (CORE::CLN::ClnWrapper::precision_ - CLN_START_PRECISION) / CLN_INCREMENT_STEP;
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
          int nsize = CORE::CLN::ClnWrapper::precision_;
          string_buffer << nsize;
          std::string clnumstr = "0e+1_" + string_buffer.str();
          zeros_cache[prec_id].first = clnumstr.c_str();
          zeros_cache[prec_id].second = true;
        }
        return zeros_cache[prec_id].first;
      }
      using maptype = std::unordered_map<double, cln::cl_F>;
      // initialize look-up vector  with maximum number of different precisions
      static std::vector<maptype> clnval_cache(
          static_cast<int>(double(CLN_REFERENCE_PREC) / double(CLN_INCREMENT_STEP)));
      auto it = clnval_cache[prec_id].find(a);
      if (it != clnval_cache[prec_id].end())
      {
        return it->second;
      }
      else
      {
        // check if this happenning during the const memory container, otherwise it will be never
        // freed (until destruction of static variables, in the end of the program ), but probably
        // that is fine
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

        return (in_ins.first->second);
      }
    }

    friend ::std::ostream& operator<<(std::ostream& stream, const CORE::CLN::ClnWrapper& a)
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

}  // namespace CORE::CLN

BACI_NAMESPACE_CLOSE

#endif
