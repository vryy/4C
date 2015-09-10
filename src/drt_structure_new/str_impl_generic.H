/*
 * str_impl_generic.H
 *
 *  Created on: Aug 17, 2015
 *      Author: farah
 */

#ifndef STR_IMPL_GENERIC_H_
#define STR_IMPL_GENERIC_H_

namespace STR {
namespace IMPLICIT {
  class Generic
  {
  public:
    //! constructor
    Generic() {};

    //! destructor
    virtual ~Generic() {};

    void Init();

    void Setup();

  protected:
    //! returns init state
    const bool& IsInit() const { return isInit_; };

    //! returns setup state
    const bool& IsSetup() const { return isSetup_; };

protected:
    //! indicates if the Init() function has been called
    bool isInit_;

    //! indicates if the Setup() function has been called
    bool isSetup_;

  };
} // namespace IMPLICIT
} // namespace STR


#endif /* STR_IMPL_GENERIC_H_ */
