/*----------------------------------------------------------------------*/
/*! \file

\brief Color escape sequences for those that prefer colorful output.

\level 2


*/
/*----------------------------------------------------------------------*/

#ifndef DRT_COLORS_H
#define DRT_COLORS_H

#ifdef COLOROUTPUT
#define BLACK "[30m"
#define BLACK_LIGHT "[30;1m"
#define RED "[31m"
#define RED_LIGHT "[31;1m"
#define GREEN "[32m"
#define GREEN_LIGHT "[32;1m"
#define YELLOW "[33m"
#define YELLOW_LIGHT "[33;1m"
#define BLUE "[34m"
#define BLUE_LIGHT "[34;1m"
#define MAGENTA "[35m"
#define MAGENTA_LIGHT "[35;1m"
#define BLUE2 "[36m"
#define BLUE2_LIGHT "[36;1m"
#define GRAY "[37m"
#define GRAY_LIGHT "[37;1m"
#define END_COLOR "[m"
#else
#define BLACK ""
#define BLACK_LIGHT ""
#define RED ""
#define RED_LIGHT ""
#define GREEN ""
#define GREEN_LIGHT ""
#define YELLOW ""
#define YELLOW_LIGHT ""
#define BLUE ""
#define BLUE_LIGHT ""
#define MAGENTA ""
#define MAGENTA_LIGHT ""
#define BLUE2 ""
#define BLUE2_LIGHT ""
#define GRAY ""
#define GRAY_LIGHT ""
#define END_COLOR ""
#endif

#endif
