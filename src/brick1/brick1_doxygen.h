/*!
\addtogroup BRICK1
\brief This module 'BRICK1' contains all routines and structures necessary
for the 3D hexahedral and tetrahedral-elements.

<pre>
Maintainer: Andreas Lipka
            lipka@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/lipka/
            0711 - 685-6575
</pre>

This module 'BRICK1' contains all routines and structures necessary
for the 3D hexahedral and tetrahedral-elements.

*  Implementation:
*
*  - general:
*           -# hex with 8/20 nodes
*           -# tet with 4/10 nodes
*
*  - material laws:
*           -# St.Venant (geom.nl)
*           -# open/closed cell linear foam material
*           -# porous versions of St.Venant, foam material
*           -# Mises plasticity (large strain - Roehl)
*           -# Drucker Prager plasticity
*
*  - integration:
*           -# 8,27 point gauss -> Hex
*           -# 4,5  point gauss -> Tet
*
*  - Element input:
*
*      -#  1 BRICK1 TET4       5      8      2     12 MAT       1 GP 1 1 1 GP_TET 4
*      -#  1 BRICK1 TET10      27     43     42     10     31     36     30     17     21     19 MAT       1 GP 1 1 1 GP_TET 5
*
*      -#  1 BRICK1 HEX8      15     25     23     14     12     22     20     10 MAT       1 GP 2 2 2 GP_TET 0  HYB 0 FORM 2
*      -#  1 BRICK1 HEX20      92    104    102     90     87     99     97     85     96    103     95     91     89    101    100     88     94     98     93     86 MAT       1 GP 3 3 3 GP_TET 0 HYB=0 FORM=2



*/
