/*!----------------------------------------------------------------------
\file dstrc.H
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifndef DSTRC_H
#define DSTRC_H
/*----------------------------------------------------------------------*/
/*!
  \brief Simplify dstrc_enter and dstrc_exit calls in c++     uk 10/06
 */
/*----------------------------------------------------------------------*/
#ifdef __cplusplus
// Depreciated! Do not use any more!
class DSTraceHelper
{
public:
  explicit DSTraceHelper(char* functname)
    {}
};
#endif

#endif  /* #ifndef DSTRC_H */
