C--------------------------------------------------------------------------
C
C \file
C \brief shape functions and their natural derivatives for fluid3 element
C
C
C <pre>
C Maintainer: Malte Neumann
C             neumann@statik.uni-stuttgart.de
C             http://www.uni-stuttgart.de/ibs/members/neumann/
C             0711 - 685-6121
C </pre>
C
C--------------------------------------------------------------------------



C--------------------------------------------------------------------------
C
C \brief shape functions and their natural derivatives for hexaeder
C
C In this routine the shape functions and their natural first and second
C derivatives with respect to r/s/t are evaluated for H E X A H E D E R
C
C \param funct()         real*8   (o) shape functions
C \param deriv()()       real*8   (o) deriv. of shape funcs
C \param deriv2()()      real*8   (o) 2nd deriv. of sh. funcs
C \param r               real*8   (i) coords of gauss points
C \param s               real*8   (i) coords of gauss points
C \param t               real*8   (i) coords of gauss points
C \param typ             integer  (i) element typ
C \param icode           integer  (i) evaluation flag
C \param sizevec(6)      integer  (i) some sizes
C
C \return void
C
C \author mn
C \date   10/04
C
C--------------------------------------------------------------------------
      subroutine f3fhex( funct, deriv, deriv2, r, s, t, typ,
     &         icode, sizevec)

      implicit none

      integer sizevec(6)
      real*8 funct(sizevec(1))
      real*8 deriv(sizevec(1),3)
      real*8 deriv2(sizevec(1),6)
      real*8 r,s,t
      integer typ
      integer icode

      real*8 zero,one,five,TWO
      real*8 rp,rm,sp,sm,tp,tm, q18 ,Q12,Q14
      real*8 rrm,ssm,ttm
      real*8 drm1,dr00,drp1,dsm1,ds00,dsp1,dtm1,dt00,dtp1
      real*8 rm1,r00,rp1,sm1,s00,sp1,tm1,t00,tp1

      TWO=2.0
      zero=0.0
      one=1.0
      five=5.0
      q18 =one/8.0
      Q12=one/2.0
      Q14=one/4.0

      if (typ.EQ.8) then
        rp=one+r
        rm=one-r
        sp=one+s
        sm=one-s
        tp=one+t
        tm=one-t

        funct(1)=q18*rp*sm*tm
        funct(2)=q18*rp*sp*tm
        funct(3)=q18*rm*sp*tm
        funct(4)=q18*rm*sm*tm
        funct(5)=q18*rp*sm*tp
        funct(6)=q18*rp*sp*tp
        funct(7)=q18*rm*sp*tp
        funct(8)=q18*rm*sm*tp

        if(icode.GT.1) then

          deriv(1,1)= q18*sm*tm
          deriv(2,1)= q18*sp*tm
          deriv(3,1)=-deriv(2,1)
          deriv(4,1)=-deriv(1,1)
          deriv(5,1)= q18*sm*tp
          deriv(6,1)= q18*sp*tp
          deriv(7,1)=-deriv(6,1)
          deriv(8,1)=-deriv(5,1)

          deriv(1,2)=-q18*tm*rp
          deriv(2,2)=-deriv(1,2)
          deriv(3,2)= q18*tm*rm
          deriv(4,2)=-deriv(3,2)
          deriv(5,2)=-q18*tp*rp
          deriv(6,2)=-deriv(5,2)
          deriv(7,2)= q18*tp*rm
          deriv(8,2)=-deriv(7,2)

          deriv(1,3)=-q18*rp*sm
          deriv(2,3)=-q18*rp*sp
          deriv(3,3)=-q18*rm*sp
          deriv(4,3)=-q18*rm*sm
          deriv(5,3)=-deriv(1,3)
          deriv(6,3)=-deriv(2,3)
          deriv(7,3)=-deriv(3,3)
          deriv(8,3)=-deriv(4,3)

        endif
        if (icode.EQ.3) then

          deriv2(1,1) =  zero
          deriv2(1,2) =  zero
          deriv2(1,3) =  zero
          deriv2(1,4) = -q18*tm
          deriv2(1,5) = -q18*sm
          deriv2(1,6) =  q18*rp

          deriv2(2,1) =  zero
          deriv2(2,2) =  zero
          deriv2(2,3) =  zero
          deriv2(2,4) = -deriv2(1,4)
          deriv2(2,5) = -q18*sp
          deriv2(2,6) = -deriv2(1,6)

          deriv2(3,1) =  zero
          deriv2(3,2) =  zero
          deriv2(3,3) =  zero
          deriv2(3,4) =  deriv2(1,4)
          deriv2(3,5) = -deriv2(2,5)
          deriv2(3,6) = -q18*rm

          deriv2(4,1) =  zero
          deriv2(4,2) =  zero
          deriv2(4,3) =  zero
          deriv2(4,4) = -deriv2(1,4)
          deriv2(4,5) = -deriv2(1,5)
          deriv2(4,6) = -deriv2(3,6)

          deriv2(5,1) =  zero
          deriv2(5,2) =  zero
          deriv2(5,3) =  zero
          deriv2(5,4) = -q18*tp
          deriv2(5,5) = -deriv2(1,5)
          deriv2(5,6) = -deriv2(1,6)

          deriv2(6,1) =  zero
          deriv2(6,2) =  zero
          deriv2(6,3) =  zero
          deriv2(6,4) = -deriv2(5,4)
          deriv2(6,5) = -deriv2(2,5)
          deriv2(6,6) =  deriv2(1,6)

          deriv2(7,1) =  zero
          deriv2(7,2) =  zero
          deriv2(7,3) =  zero
          deriv2(7,4) =  deriv2(5,4)
          deriv2(7,5) =  deriv2(2,5)
          deriv2(7,6) = -deriv2(3,6)

          deriv2(8,1) =  zero
          deriv2(8,2) =  zero
          deriv2(8,3) =  zero
          deriv2(8,4) = -deriv2(5,4)
          deriv2(8,5) =  deriv2(1,5)
          deriv2(8,6) =  deriv2(3,6)

        endif

        elseif (typ.EQ.20) then

c       form basic values
        rp=one+r
        rm=one-r
        sp=one+s
        sm=one-s
        tp=one+t
        tm=one-t
        rrm=one-r*r
        ssm=one-s*s
        ttm=one-t*t

        funct(1 ) =q18*rp*sm*tm*(rp+sm+tm-five)
        funct(2 ) =q18*rp*sp*tm*(rp+sp+tm-five)
        funct(3 ) =q18*rm*sp*tm*(rm+sp+tm-five)
        funct(4 ) =q18*rm*sm*tm*(rm+sm+tm-five)
        funct(5 ) =q18*rp*sm*tp*(rp+sm+tp-five)
        funct(6 ) =q18*rp*sp*tp*(rp+sp+tp-five)
        funct(7 ) =q18*rm*sp*tp*(rm+sp+tp-five)
        funct(8 ) =q18*rm*sm*tp*(rm+sm+tp-five)
        funct(9 ) =q14*rp*ssm*tm
        funct(10) =q14*rrm*sp*tm
        funct(11)=q14*rm*ssm*tm
        funct(12)=q14*rrm*sm*tm
        funct(13)=q14*rp*ssm*tp
        funct(14)=q14*rrm*sp*tp
        funct(15)=q14*rm*ssm*tp
        funct(16)=q14*rrm*sm*tp
        funct(17)=q14*rp*sm*ttm
        funct(18)=q14*rp*sp*ttm
        funct(19)=q14*rm*sp*ttm
        funct(20)=q14*rm*sm*ttm

        if(icode.GT.1) then

          deriv(1,1) = q18*sm*tm*(TWO*rp+sm+tm-five)
          deriv(1,2) =-q18*tm*rp*(TWO*sm+tm+rp-five)
          deriv(1,3) =-q18*rp*sm*(TWO*tm+rp+sm-five)

          deriv(2,1) = q18*sp*tm*(TWO*rp+sp+tm-five)
          deriv(2,2) = q18*tm*rp*(TWO*sp+tm+rp-five)
          deriv(2,3) =-q18*rp*sp*(TWO*tm+rp+sp-five)

          deriv(3,1) =-q18*sp*tm*(TWO*rm+sp+tm-five)
          deriv(3,2) = q18*tm*rm*(TWO*sp+tm+rm-five)
          deriv(3,3) =-q18*rm*sp*(TWO*tm+rm+sp-five)

          deriv(4,1) =-q18*sm*tm*(TWO*rm+sm+tm-five)
          deriv(4,2) =-q18*tm*rm*(TWO*sm+tm+rm-five)
          deriv(4,3) =-q18*rm*sm*(TWO*tm+rm+sm-five)

          deriv(5,1) = q18*sm*tp*(TWO*rp+sm+tp-five)
          deriv(5,2) =-q18*tp*rp*(TWO*sm+tp+rp-five)
          deriv(5,3) = q18*rp*sm*(TWO*tp+rp+sm-five)

          deriv(6,1) = q18*sp*tp*(TWO*rp+sp+tp-five)
          deriv(6,2) = q18*tp*rp*(TWO*sp+tp+rp-five)
          deriv(6,3) = q18*rp*sp*(TWO*tp+rp+sp-five)

          deriv(7,1) =-q18*sp*tp*(TWO*rm+sp+tp-five)
          deriv(7,2) = q18*tp*rm*(TWO*sp+tp+rm-five)
          deriv(7,3) = q18*rm*sp*(TWO*tp+rm+sp-five)

          deriv(8,1) =-q18*sm*tp*(TWO*rm+sm+tp-five)
          deriv(8,2) =-q18*tp*rm*(TWO*sm+tp+rm-five)
          deriv(8,3) = q18*rm*sm*(TWO*tp+rm+sm-five)

          deriv(9,1) = q14*ssm*tm
          deriv(9,2) =-q12*s*tm*rp
          deriv(9,3) =-q14*ssm*rp

          deriv(10,1) =-q12*r*sp*tm
          deriv(10,2) = q14*rrm*tm
          deriv(10,3) =-q14*rrm*sp

          deriv(11,1)=-deriv(9,1)
          deriv(11,2)=-q12*s*tm*rm
          deriv(11,3)=-q14*ssm*rm

          deriv(12,1)=-q12*r*sm*tm
          deriv(12,2)=-deriv(10,2)
          deriv(12,3)=-q14*rrm*sm

          deriv(13,1)= q14*ssm*tp
          deriv(13,2)=-q12*s*tp*rp
          deriv(13,3)=-deriv(9,3)

          deriv(14,1)=-q12*r*sp*tp
          deriv(14,2)= q14*rrm*tp
          deriv(14,3)=-deriv(9,3)

          deriv(15,1)=-deriv(13,1)
          deriv(15,2)=-q12*s*tp*rm
          deriv(15,3)=-deriv(11,3)

          deriv(16,1)=-q12*r*sm*tp
          deriv(16,2)=-deriv(14,2)
          deriv(16,3)=-deriv(12,3)

          deriv(17,1)= q14*sm*ttm
          deriv(17,2)=-q14*ttm*rp
          deriv(17,3)=-q12*t*rp*sm

          deriv(18,1)= q14*sp*ttm
          deriv(18,2)=-deriv(17,2)
          deriv(18,3)=-q12*t*rp*sp

          deriv(19,1)=-deriv(18,1)
          deriv(19,2)= q14*ttm*rm
          deriv(19,3)=-q12*t*rm*sp

          deriv(20,1)=-deriv(17,1)
          deriv(20,2)=-deriv(19,2)
          deriv(20,3)=-q12*t*rm*sm

        endif
        if(icode.EQ.3) then

          deriv2(1,1) = q14*sm*tm
          deriv2(1,2) = q14*tm*rp
          deriv2(1,3) = q14*rp*sm
          deriv2(1,4) =-q18*(tm*(2*rp+sm+tm-five+sm*tm))
          deriv2(1,5) =-q18*(sm*(2*rp+sm+tm-five+sm*tm))
          deriv2(1,6) = q18*(rp*(2*sm+tm+rp-five+tm*rp))

          deriv2(2,1) = q14*sp*tm
          deriv2(2,2) = deriv2(2,3)
          deriv2(2,3) = q14*rp*sp
          deriv2(2,4) =-q18*(tm*(2*rp+sp+tm-five+sp*tm))
          deriv2(2,5) =-q18*(sp*(2*rp+sp+tm-five+sp*tm))
          deriv2(2,6) =-q18*(rp*(2*sp+tm+rp-five+tm*rp))

          deriv2(3,1) =-deriv2(3,2)
          deriv2(3,2) = q14*tm*rm
          deriv2(3,3) = q14*rm*sp
          deriv2(3,4) =-q18*(tm*(2*rm+sp+tm-five+sp*tm))
          deriv2(3,5) = q18*(sp*(2*rm+sp+tm-five+sp*tm))
          deriv2(3,6) =-q18*(rm*(2*sp+tm+rm-five+tm*rm))

          deriv2(4,1) =-deriv2(2,1)
          deriv2(4,2) = deriv2(4,3)
          deriv2(4,3) = q14*rm*sm
          deriv2(4,4) =-q18*(tm*(2*rm+sm+tm-five+sm*tm))
          deriv2(4,5) = q18*(sm*(2*rm+sm+tm-five+sm*tm))
          deriv2(4,6) = q18*(rm*(2*sm+tm+rm-five+tm*rm))

          deriv2(5,1) = q14*sm*tp
          deriv2(5,2) = q14*tp*rp
          deriv2(5,3) = deriv2(2,4)
          deriv2(5,4) =-q18*(tp*(2*rp+sm+tp-five+sm*tp))
          deriv2(5,5) = q18*(sm*(2*rp+sm+tp-five+sm*tp))
          deriv2(5,6) =-q18*(rp*(2*sm+tp+rp-five+tp*rp))

          deriv2(6,1) = q14*sp*tp
          deriv2(6,2) = deriv2(6,3)
          deriv2(6,3) = deriv2(3,4)
          deriv2(6,4) =-q18*(tp*(2*rp+sp+tp-five+sp*tp))
          deriv2(6,5) = q18*(sp*(2*rp+sp+tp-five+sp*tp))
          deriv2(6,6) = q18*(rp*(2*sp+tp+rp-five+tp*rp))

          deriv2(7,1) =-deriv2(7,2)
          deriv2(7,2) = q14*tp*rm
          deriv2(7,3) = deriv2(4,4)
          deriv2(7,4) =-q18*(tp*(2*rm+sp+tp-five+sp*tp))
          deriv2(7,5) =-q18*(sp*(2*rm+sp+tp-five+sp*tp))
          deriv2(7,6) = q18*(rm*(2*sp+tp+rm-five+tp*rm))

          deriv2(8,1) =-deriv2(6,2)
          deriv2(8,2) = deriv2(8,3)
          deriv2(8,3) = deriv2(5,4)
          deriv2(8,4) =-q18*(tp*(2*rm+sm+tp-five+sm*tp))
          deriv2(8,5) =-q18*(sm*(2*rm+sm+tp-five+sm*tp))
          deriv2(8,6) =-q18*(rm*(2*sm+tp+rm-five+tp*rm))

          deriv2(9,1) = zero
          deriv2(9,2) = -q12*tm*rp
          deriv2(9,3) = zero
          deriv2(9,4) =-q12*s*tm
          deriv2(9,5) =-q14*ssm
          deriv2(9,6) = q12*s*rp

          deriv2(10,1)=-q12*sp*tm
          deriv2(10,2)= zero
          deriv2(10,3)= zero
          deriv2(10,4)=-q12*r*tm
          deriv2(10,5)= q12*r*sp
          deriv2(10,6)=-q14*rrm

          deriv2(11,1)= zero
          deriv2(11,2)= -q12*tm*rm
          deriv2(11,3)= zero
          deriv2(11,4)= q12*s*tm
          deriv2(11,5)=-deriv2(9,5)
          deriv2(11,6)= q12*s*rm

          deriv2(12,1)=-q12*sm*tm
          deriv2(12,2)= zero
          deriv2(12,3)= zero
          deriv2(12,4)= q12*r*tm
          deriv2(12,5)= q12*r*sm
          deriv2(12,6)=-deriv2(10,6)

          deriv2(13,1)= zero
          deriv2(13,2)= -q12*tp*rp
          deriv2(13,3)= zero
          deriv2(13,4)=-q12*s*tp
          deriv2(13,5)=-deriv2(9,5)
          deriv2(13,6)=-deriv2(9,6)

          deriv2(14,1)=-q12*sp*tp
          deriv2(14,2)= zero
          deriv2(14,3)= zero
          deriv2(14,4)=-q12*r*tp
          deriv2(14,5)=-deriv2(10,5)
          deriv2(14,6)=-deriv2(10,6)

          deriv2(15,1)= zero
          deriv2(15,2)= -q12*tp*rm
          deriv2(15,3)= zero
          deriv2(15,4)= q12*s*tp
          deriv2(15,5)= deriv2(9,5)
          deriv2(15,6)=-deriv2(11,6)

          deriv2(16,1)=-q12*sm*tp
          deriv2(16,2)= zero
          deriv2(16,3)= zero
          deriv2(16,4)= q12*r*tp
          deriv2(16,5)=-deriv2(12,5)
          deriv2(16,6)= deriv2(10,6)

          deriv2(17,1)= zero
          deriv2(17,2)= zero
          deriv2(17,3)= zero
          deriv2(17,4)=-q14*ttm
          deriv2(17,5)=-q12*t*sm
          deriv2(17,6)= q12*t*rp

          deriv2(18,1)= zero
          deriv2(18,2)= zero
          deriv2(18,3)= zero
          deriv2(18,4)= q14*ttm
          deriv2(18,5)=-q12*t*sp
          deriv2(18,6)=-deriv2(17,6)

          deriv2(19,1)= zero
          deriv2(19,2)= zero
          deriv2(19,3)= zero
          deriv2(19,4)= deriv2(17,4)
          deriv2(19,5)= q12*t*sp
          deriv2(19,6)= q12*t*rm

          deriv2(20,1)= zero
          deriv2(20,2)= zero
          deriv2(20,3)= zero
          deriv2(20,4)= deriv2(18,4)
          deriv2(20,5)= q12*t*sm
          deriv2(20,6)=-deriv2(19,6)

        endif

        elseif(typ.EQ.27)  then
c       quadratic shape functions and their natural derivatives with central nodes

c       form basic values
        rm1=q12*r*(r - one)
        r00=(one - r*r)
        rp1=q12*r*(r + one)
        sm1=q12*s*(s - one)
        s00=(one - s*s)
        sp1=q12*s*(s + one)
        tm1=q12*t*(t - one)
        t00=(one - t*t)
        tp1=q12*t*(t + one)

        funct(1 )= rp1 * sm1 * tm1
        funct(2 )= rp1 * sp1 * tm1
        funct(3 )= rm1 * sp1 * tm1
        funct(4 )= rm1 * sm1 * tm1
        funct(5 )= rp1 * sm1 * tp1
        funct(6 )= rp1 * sp1 * tp1
        funct(7 )= rm1 * sp1 * tp1
        funct(8 )= rm1 * sm1 * tp1
        funct(9 )= rp1 * s00 * tm1
        funct(10)= r00 * sp1 * tm1
        funct(11)= rm1 * s00 * tm1
        funct(12)= r00 * sm1 * tm1
        funct(13)= rp1 * s00 * tp1
        funct(14)= r00 * sp1 * tp1
        funct(15)= rm1 * s00 * tp1
        funct(16)= r00 * sm1 * tp1
        funct(17)= rp1 * sm1 * t00
        funct(18)= rp1 * sp1 * t00
        funct(19)= rm1 * sp1 * t00
        funct(20)= rm1 * sm1 * t00
        funct(21)= rp1 * s00 * t00
        funct(22)= r00 * sp1 * t00
        funct(23)= rm1 * s00 * t00
        funct(24)= r00 * sm1 * t00
        funct(25)= r00 * s00 * tp1
        funct(26)= r00 * s00 * tm1
        funct(27)= r00 * s00 * t00

        if(icode.GT.1) then

          drm1 = r - q12
          dr00 = -TWO * r
          drp1 = r + q12
          dsm1 = s - q12
          ds00 = -TWO * s
          dsp1 = s + q12
          dtm1 = t - q12
          dt00 = -TWO * t
          dtp1 = t + q12

          deriv(1 ,1)= drp1 * sm1 * tm1
          deriv(2 ,1)= drp1 * sp1 * tm1
          deriv(3 ,1)= drm1 * sp1 * tm1
          deriv(4 ,1)= drm1 * sm1 * tm1
          deriv(5 ,1)= drp1 * sm1 * tp1
          deriv(6 ,1)= drp1 * sp1 * tp1
          deriv(7 ,1)= drm1 * sp1 * tp1
          deriv(8 ,1)= drm1 * sm1 * tp1
          deriv(9 ,1)= drp1 * s00 * tm1
          deriv(10,1)= dr00 * sp1 * tm1
          deriv(11,1)= drm1 * s00 * tm1
          deriv(12,1)= dr00 * sm1 * tm1
          deriv(13,1)= drp1 * s00 * tp1
          deriv(14,1)= dr00 * sp1 * tp1
          deriv(15,1)= drm1 * s00 * tp1
          deriv(16,1)= dr00 * sm1 * tp1
          deriv(17,1)= drp1 * sm1 * t00
          deriv(18,1)= drp1 * sp1 * t00
          deriv(19,1)= drm1 * sp1 * t00
          deriv(20,1)= drm1 * sm1 * t00
          deriv(21,1)= drp1 * s00 * t00
          deriv(22,1)= dr00 * sp1 * t00
          deriv(23,1)= drm1 * s00 * t00
          deriv(24,1)= dr00 * sm1 * t00
          deriv(25,1)= dr00 * s00 * tp1
          deriv(26,1)= dr00 * s00 * tm1
          deriv(27,1)= dr00 * s00 * t00

          deriv(1,2)= rp1 * dsm1 * tm1
          deriv(2,2)= rp1 * dsp1 * tm1
          deriv(3,2)= rm1 * dsp1 * tm1
          deriv(4,2)= rm1 * dsm1 * tm1
          deriv(5,2)= rp1 * dsm1 * tp1
          deriv(6,2)= rp1 * dsp1 * tp1
          deriv(7,2)= rm1 * dsp1 * tp1
          deriv(8,2)= rm1 * dsm1 * tp1
          deriv(9,2)= rp1 * ds00 * tm1
          deriv(10,2)= r00 * dsp1 * tm1
          deriv(11,2)= rm1 * ds00 * tm1
          deriv(12,2)= r00 * dsm1 * tm1
          deriv(13,2)= rp1 * ds00 * tp1
          deriv(14,2)= r00 * dsp1 * tp1
          deriv(15,2)= rm1 * ds00 * tp1
          deriv(16,2)= r00 * dsm1 * tp1
          deriv(17,2)= rp1 * dsm1 * t00
          deriv(18,2)= rp1 * dsp1 * t00
          deriv(19,2)= rm1 * dsp1 * t00
          deriv(20,2)= rm1 * dsm1 * t00
          deriv(21,2)= rp1 * ds00 * t00
          deriv(22,2)= r00 * dsp1 * t00
          deriv(23,2)= rm1 * ds00 * t00
          deriv(24,2)= r00 * dsm1 * t00
          deriv(25,2)= r00 * ds00 * tp1
          deriv(26,2)= r00 * ds00 * tm1
          deriv(27,2)= r00 * ds00 * t00

          deriv(1,3)= rp1 * sm1 * dtm1
          deriv(2,3)= rp1 * sp1 * dtm1
          deriv(3,3)= rm1 * sp1 * dtm1
          deriv(4,3)= rm1 * sm1 * dtm1
          deriv(5,3)= rp1 * sm1 * dtp1
          deriv(6,3)= rp1 * sp1 * dtp1
          deriv(7,3)= rm1 * sp1 * dtp1
          deriv(8,3)= rm1 * sm1 * dtp1
          deriv(9,3)= rp1 * s00 * dtm1
          deriv(10,3)= r00 * sp1 * dtm1
          deriv(11,3)= rm1 * s00 * dtm1
          deriv(12,3)= r00 * sm1 * dtm1
          deriv(13,3)= rp1 * s00 * dtp1
          deriv(14,3)= r00 * sp1 * dtp1
          deriv(15,3)= rm1 * s00 * dtp1
          deriv(16,3)= r00 * sm1 * dtp1
          deriv(17,3)= rp1 * sm1 * dt00
          deriv(18,3)= rp1 * sp1 * dt00
          deriv(19,3)= rm1 * sp1 * dt00
          deriv(20,3)= rm1 * sm1 * dt00
          deriv(21,3)= rp1 * s00 * dt00
          deriv(22,3)= r00 * sp1 * dt00
          deriv(23,3)= rm1 * s00 * dt00
          deriv(24,3)= r00 * sm1 * dt00
          deriv(25,3)= r00 * s00 * dtp1
          deriv(26,3)= r00 * s00 * dtm1
          deriv(27,3)= r00 * s00 * dt00
        endif

      endif

      return
      end  subroutine


C--------------------------------------------------------------------------
C
C \brief shape functions and their natural derivatives for tetraeder
C
C In this routine the shape functions and their natural first and second
C derivatives with respect to r/s/t are evaluated for T E T R A E D E R
C
C \param funct()         real*8   (o) shape functions
C \param deriv()()       real*8   (o) deriv. of shape funcs
C \param deriv2()()      real*8   (o) 2nd deriv. of sh. funcs
C \param r               real*8   (i) coords of gauss points
C \param s               real*8   (i) coords of gauss points
C \param t               real*8   (i) coords of gauss points
C \param typ             integer  (i) element typ
C \param icode           integer  (i) evaluation flag
C \param sizevec(6)      integer  (i) some sizes
C
C \return void
C
C \author mn
C \date   10/04
C
C--------------------------------------------------------------------------
      subroutine f3ftet(funct, deriv, deriv2, r, s, t, typ,
     &        icode, sizevec)
      implicit none

      integer sizevec(6)
      real*8 funct(sizevec(1))
      real*8 deriv(sizevec(1),3)
      real*8 deriv2(sizevec(1),6)
      real*8 r,s,t
      integer typ
      integer icode

      real*8 zero,one,five
      real*8 t1,t2,t3,t4

c     form basic values
      zero = 0.0
      one  = 1.0
      five = 5.0

      if(typ.EQ.4) then
        t1=r
        t2=s
        t3=t
        t4=one-r-s-t

        funct(1)= t1
        funct(2)= t2
        funct(3)= t3
        funct(4)= t4

        if(icode.GT.1) then

          deriv(1,1)= one
          deriv(2,1)= zero
          deriv(3,1)= zero
          deriv(4,1)=-one

          deriv(1,2)= zero
          deriv(2,2)= one
          deriv(3,2)= zero
          deriv(4,2)=-one

          deriv(1,3)= zero
          deriv(2,3)= zero
          deriv(3,3)= one
          deriv(4,3)=-one
        endif

      endif

      return
      end subroutine





C--------------------------------------------------------------------------
C
C \brief jacobian matrix
C
C In this routine the jacobian matrix and its determinant is calculated.
C
C \param funct()         real*8   (i) shape functions
C \param deriv()()       real*8   (i) deriv. of shape funcs
C \param xjm()()()       real*8   (o) jacobian matrix
C \param det()           real*8   (o) determinant of the jacobian
C \param elecord()()()   real*8   (i) vector containing nodal coordinates
C \param sizevec(6)      integer  (i) some sizes
C
C \return void
C
C \author mn
C \date   10/04
C
C--------------------------------------------------------------------------
      subroutine f3fjaco(funct, deriv, xjm, det, elecord, sizevec)

      implicit none

      integer sizevec(6)
      real*8 funct(sizevec(1))
      real*8 deriv(sizevec(1),3)
      real*8 xjm(sizevec(4),3,3)
      real*8 det(sizevec(4))
      real*8 elecord(sizevec(4),sizevec(1),3)

      integer i,j,l,k
      real*8 dum

c     determine jacobian at point r,s,t
      do i=1,3
        do j=1,3

          do k=1,sizevec(5)
            dum=0.0
            do l=1,sizevec(2)
              dum =dum+ deriv(l,i)*elecord(k,l,j)
            enddo
            xjm(k,j,i)=dum
          enddo

        enddo
      enddo

c     determinant of jacobian
      do k=1,sizevec(5)

        det(k) =
     &  xjm(k,1,1)*xjm(k,2,2)*xjm(k,3,3)+
     &  xjm(k,2,1)*xjm(k,3,2)*xjm(k,1,3)+
     &  xjm(k,3,1)*xjm(k,1,2)*xjm(k,2,3)-
     &  xjm(k,3,1)*xjm(k,2,2)*xjm(k,1,3)-
     &  xjm(k,1,1)*xjm(k,3,2)*xjm(k,2,3)-
     &  xjm(k,2,1)*xjm(k,1,2)*xjm(k,3,3)

      enddo

      return
      end subroutine


C--------------------------------------------------------------------------
C
C \brief global derivates
C
C In this routine the global derivatives w.r.t. x,y,z at point r,s,t are
C calculated.
C
C \param derxy()()()     real*8   (o) global derivatives
C \param deriv()()       real*8   (i) deriv. of shape funcs
C \param xjm()()()       real*8   (i) jacobian matrix
C \param xjm()()()       real*8   (i) inverse of jacobian matrix
C \param det()           real*8   (i) determinant of the jacobian
C \param sizevec(6)      integer  (i) some sizes
C
C \return void
C
C \author mn
C \date   10/04
C
C--------------------------------------------------------------------------
      subroutine f3fgder(derxy, deriv, xjm, xji, det, sizevec)

      implicit none

      integer sizevec(6)
      real*8 derxy(sizevec(4),sizevec(1),3)
      real*8 deriv(sizevec(1),3)
      real*8 xjm(sizevec(4),3,3)
      real*8 xji(sizevec(4),32,32)
      real*8 det(sizevec(4))

      integer  k,j
      real*8 zero
      real*8 invdet(sizevec(4))


c     initialistion
      zero=0.0

      do k=1,sizevec(2)
        do j=1,sizevec(5)
          derxy(j,k,1)=zero
          derxy(j,k,2)=zero
          derxy(j,k,3)=zero
        enddo
      enddo

c     inverse of jacobian
      do j=1,sizevec(5)
        invdet(j) = 1.0 / det(j)

        xji(j,1,1) = (  xjm(j,2,2)*xjm(j,3,3) - xjm(j,2,3)*xjm(j,3,2))
     &             *invdet(j)
        xji(j,1,2) = (- xjm(j,1,2)*xjm(j,3,3) + xjm(j,1,3)*xjm(j,3,2))
     &            *invdet(j)
        xji(j,1,3) = (  xjm(j,1,2)*xjm(j,2,3) - xjm(j,1,3)*xjm(j,2,2))
     &            *invdet(j)
        xji(j,2,1) = (- xjm(j,2,1)*xjm(j,3,3) + xjm(j,2,3)*xjm(j,3,1))
     &            *invdet(j)
        xji(j,2,2) = (  xjm(j,1,1)*xjm(j,3,3) - xjm(j,1,3)*xjm(j,3,1))
     &            *invdet(j)
        xji(j,2,3) = (- xjm(j,1,1)*xjm(j,2,3) + xjm(j,1,3)*xjm(j,2,1))
     &            *invdet(j)
        xji(j,3,1) = (  xjm(j,2,1)*xjm(j,3,2) - xjm(j,2,2)*xjm(j,3,1))
     &            *invdet(j)
        xji(j,3,2) = (- xjm(j,1,1)*xjm(j,3,2) + xjm(j,1,2)*xjm(j,3,1))
     &            *invdet(j)
        xji(j,3,3) = (  xjm(j,1,1)*xjm(j,2,2) - xjm(j,1,2)*xjm(j,2,1))
     &            *invdet(j)
      enddo


c     calculate global derivatives
      do k=1,sizevec(2)
        do j=1,sizevec(5)
          derxy(j,k,1) = derxy(j,k,1) + xji(j,1,1) * deriv(k,1)
     &  + xji(j,2,1) * deriv(k,2) + xji(j,3,1) * deriv(k,3)
          derxy(j,k,2) = derxy(j,k,2) + xji(j,1,2) * deriv(k,1)
     &  + xji(j,2,2) * deriv(k,2) + xji(j,3,2) * deriv(k,3)
          derxy(j,k,3) = derxy(j,k,3) + xji(j,1,3) * deriv(k,1)
     &  + xji(j,2,3) * deriv(k,2) + xji(j,3,3) * deriv(k,3)
        enddo
      enddo


      return
      end subroutine






C--------------------------------------------------------------------------
C
C \brief second global derivatives (looping version)
C
C In this routine the second global derivatives w.r.t x/y/z at point r,s,t
C are calculated.
C
C \param elecord()()()   real*8   (i) vector containing nodal coordinates
C \param xjm()()()       real*8   (i) jacobian matrix
C \param bm_in()()()     real*8   (i) working array
C \param xder2_in()()()  real*8   (i) working array
C \param derxy()()()     real*8   (i) global derivatives
C \param derxy2()()()    real*8   (o) 2nd global derivatives
C \param sizevec(6)      integer  (i) some sizes
C
C \return void
C
C \author mn
C \date   10/04
C
C--------------------------------------------------------------------------
      subroutine f3fgder2loop(elecord_f, xjm, bm_in, xder2_in, derxy,
     &          derxy2, deriv2, sizevec)

      implicit none

      integer sizevec(6)
      real*8 elecord_f(sizevec(4),sizevec(1),3)
      real*8 xjm(sizevec(4),3,3)
      real*8 bm_in(sizevec(4),32,32)
      real*8 xder2_in(sizevec(4),32,32)
      real*8 derxy(sizevec(4),sizevec(1),3)
      real*8 derxy2(sizevec(4),sizevec(1),6)
      real*8 deriv2(sizevec(1),6)

      integer i,j,value,k
      real*8 r0,r1,r2,r3,r4,r5
      real*8 TWO


c     calculate elements of jacobian_bar matrix
      TWO=2.0

      do k=1,sizevec(5)

        bm_in(k,1,1) = xjm(k,1,1)*xjm(k,1,1)
        bm_in(k,1,2) = xjm(k,1,2)*xjm(k,1,2)
        bm_in(k,1,3) = xjm(k,1,3)*xjm(k,1,3)
        bm_in(k,1,4) = xjm(k,1,1)*xjm(k,1,2)
        bm_in(k,1,5) = xjm(k,1,1)*xjm(k,1,3)
        bm_in(k,1,6) = xjm(k,1,2)*xjm(k,1,3)

        bm_in(k,2,1) = xjm(k,2,1)*xjm(k,2,1)
        bm_in(k,2,2) = xjm(k,2,2)*xjm(k,2,2)
        bm_in(k,2,3) = xjm(k,2,3)*xjm(k,2,3)
        bm_in(k,2,4) = xjm(k,2,1)*xjm(k,2,2)
        bm_in(k,2,5) = xjm(k,2,1)*xjm(k,2,3)
        bm_in(k,2,6) = xjm(k,2,2)*xjm(k,2,3)

        bm_in(k,3,1) = xjm(k,3,1)*xjm(k,3,1)
        bm_in(k,3,2) = xjm(k,3,2)*xjm(k,3,2)
        bm_in(k,3,3) = xjm(k,3,3)*xjm(k,3,3)
        bm_in(k,3,4) = xjm(k,3,1)*xjm(k,3,2)
        bm_in(k,3,5) = xjm(k,3,1)*xjm(k,3,3)
        bm_in(k,3,6) = xjm(k,3,2)*xjm(k,3,3)

        bm_in(k,4,1) = TWO*xjm(k,1,1)*xjm(k,2,1)
        bm_in(k,4,2) = TWO*xjm(k,1,2)*xjm(k,2,2)
        bm_in(k,4,3) = TWO*xjm(k,1,3)*xjm(k,2,3)
        bm_in(k,4,4) =     xjm(k,1,1)*xjm(k,2,2)+xjm(k,1,2)*xjm(k,2,1)
        bm_in(k,4,5) =     xjm(k,1,1)*xjm(k,2,3)+xjm(k,1,3)*xjm(k,2,1)
        bm_in(k,4,6) =     xjm(k,1,2)*xjm(k,2,3)+xjm(k,1,3)*xjm(k,2,2)

        bm_in(k,5,1) = TWO*xjm(k,1,1)*xjm(k,3,1)
        bm_in(k,5,2) = TWO*xjm(k,1,2)*xjm(k,3,2)
        bm_in(k,5,3) = TWO*xjm(k,1,3)*xjm(k,3,3)
        bm_in(k,5,4) =     xjm(k,1,1)*xjm(k,3,2)+xjm(k,1,2)*xjm(k,3,1)
        bm_in(k,5,5) =     xjm(k,1,1)*xjm(k,3,3)+xjm(k,1,3)*xjm(k,3,1)
        bm_in(k,5,6) =     xjm(k,1,2)*xjm(k,3,3)+xjm(k,1,3)*xjm(k,3,2)

        bm_in(k,6,1) = TWO*xjm(k,2,1)*xjm(k,3,1)
        bm_in(k,6,2) = TWO*xjm(k,2,2)*xjm(k,3,2)
        bm_in(k,6,3) = TWO*xjm(k,2,3)*xjm(k,3,3)
        bm_in(k,6,4) =     xjm(k,2,1)*xjm(k,3,2)+xjm(k,2,2)*xjm(k,3,1)
        bm_in(k,6,5) =     xjm(k,2,1)*xjm(k,3,3)+xjm(k,2,3)*xjm(k,3,1)
        bm_in(k,6,6) =     xjm(k,2,2)*xjm(k,3,3)+xjm(k,2,3)*xjm(k,3,2)

      enddo


      call f3finv62(bm_in, value,sizevec(5),sizevec(4))


c     initialise
      do i=1,3
        do j=1,6
          do k=1,sizevec(5)
            xder2_in(k,i,j)=0.0
          enddo
        enddo
      enddo

      do i=1,sizevec(2)
        do j=1,6
          do k=1,sizevec(5)
            derxy2(k,i,j)=0.0
          enddo
        enddo
      enddo

c     determine 2nd derivatives of coord.-functions
      do i=1,sizevec(2)
        do k=1,sizevec(5)
          xder2_in(k,1,1)=xder2_in(k,1,1)+deriv2(i,1) * elecord_f(k,i,1)
          xder2_in(k,1,2)=xder2_in(k,1,2)+deriv2(i,2) * elecord_f(k,i,1)
          xder2_in(k,1,3)=xder2_in(k,1,3)+deriv2(i,3) * elecord_f(k,i,1)
          xder2_in(k,1,4)=xder2_in(k,1,4)+deriv2(i,4) * elecord_f(k,i,1)
          xder2_in(k,1,5)=xder2_in(k,1,5)+deriv2(i,5) * elecord_f(k,i,1)
          xder2_in(k,1,6)=xder2_in(k,1,6)+deriv2(i,6) * elecord_f(k,i,1)

          xder2_in(k,2,1)=xder2_in(k,2,1)+deriv2(i,1) * elecord_f(k,i,2)
          xder2_in(k,2,2)=xder2_in(k,2,2)+deriv2(i,2) * elecord_f(k,i,2)
          xder2_in(k,2,3)=xder2_in(k,2,3)+deriv2(i,3) * elecord_f(k,i,2)
          xder2_in(k,2,4)=xder2_in(k,2,4)+deriv2(i,4) * elecord_f(k,i,2)
          xder2_in(k,2,5)=xder2_in(k,2,5)+deriv2(i,5) * elecord_f(k,i,2)
          xder2_in(k,2,6)=xder2_in(k,2,6)+deriv2(i,6) * elecord_f(k,i,2)

          xder2_in(k,3,1)=xder2_in(k,3,1)+deriv2(i,1) * elecord_f(k,i,3)
          xder2_in(k,3,2)=xder2_in(k,3,2)+deriv2(i,2) * elecord_f(k,i,3)
          xder2_in(k,3,3)=xder2_in(k,3,3)+deriv2(i,3) * elecord_f(k,i,3)
          xder2_in(k,3,4)=xder2_in(k,3,4)+deriv2(i,4) * elecord_f(k,i,3)
          xder2_in(k,3,5)=xder2_in(k,3,5)+deriv2(i,5) * elecord_f(k,i,3)
          xder2_in(k,3,6)=xder2_in(k,3,6)+deriv2(i,6) * elecord_f(k,i,3)
        enddo
      enddo

c     calculate second global derivatives
      do i=1,sizevec(2)
        do k=1,sizevec(5)
          r0 = deriv2(i,1) - xder2_in(k,1,1)*derxy(k,i,1)
     &                   - xder2_in(k,2,1)*derxy(k,i,2)
     &                   - xder2_in(k,3,1)*derxy(k,i,3)
          r1 = deriv2(i,2) - xder2_in(k,1,2)*derxy(k,i,1)
     &                   - xder2_in(k,2,2)*derxy(k,i,2)
     &                   - xder2_in(k,3,2)*derxy(k,i,3)
          r2 = deriv2(i,3) - xder2_in(k,1,3)*derxy(k,i,1)
     &                   - xder2_in(k,2,3)*derxy(k,i,2)
     &                   - xder2_in(k,3,3)*derxy(k,i,3)
          r3 = deriv2(i,4) - xder2_in(k,1,4)*derxy(k,i,1)
     &                   - xder2_in(k,2,4)*derxy(k,i,2)
     &                   - xder2_in(k,3,4)*derxy(k,i,3)
          r4 = deriv2(i,5) - xder2_in(k,1,5)*derxy(k,i,1)
     &                   - xder2_in(k,2,5)*derxy(k,i,2)
     &                   - xder2_in(k,3,5)*derxy(k,i,3)
          r5 = deriv2(i,6) - xder2_in(k,1,6)*derxy(k,i,1)
     &                   - xder2_in(k,2,6)*derxy(k,i,2)
     &                   - xder2_in(k,3,6)*derxy(k,i,3)

          derxy2(k,i,1)=derxy2(k,i,1)+bm_in(k,1,1)*r0 + bm_in(k,2,1)*r1
     &     + bm_in(k,3,1)*r2+bm_in(k,4,1)*r3 
     &     + bm_in(k,5,1)*r4 + bm_in(k,6,1)*r5
          derxy2(k,i,2)=derxy2(k,i,2)+bm_in(k,1,2)*r0 + bm_in(k,2,2)*r1
     &     + bm_in(k,3,2)*r2+bm_in(k,4,2)*r3
     &     + bm_in(k,5,2)*r4 + bm_in(k,6,2)*r5
          derxy2(k,i,3)=derxy2(k,i,3)+bm_in(k,1,3)*r0 + bm_in(k,2,3)*r1
     &     + bm_in(k,3,3)*r2+bm_in(k,4,3)*r3
     &     + bm_in(k,5,3)*r4 + bm_in(k,6,3)*r5
          derxy2(k,i,4)=derxy2(k,i,4)+bm_in(k,1,4)*r0 + bm_in(k,2,4)*r1
     &     + bm_in(k,3,4)*r2+bm_in(k,4,4)*r3
     &     + bm_in(k,5,4)*r4 + bm_in(k,6,4)*r5
          derxy2(k,i,5)=derxy2(k,i,5)+bm_in(k,1,5)*r0 + bm_in(k,2,5)*r1
     &     + bm_in(k,3,5)*r2+bm_in(k,4,5)*r3
     &     + bm_in(k,5,5)*r4 + bm_in(k,6,5)*r5
          derxy2(k,i,6)=derxy2(k,i,6)+bm_in(k,1,6)*r0 + bm_in(k,2,6)*r1
     &     + bm_in(k,3,6)*r2+bm_in(k,4,6)*r3
     &     + bm_in(k,5,6)*r4 + bm_in(k,6,6)*r5
        enddo
      enddo


      return
      end subroutine





C--------------------------------------------------------------------------
C
C \brief second global derivatives (NO looping version)
C
C In this routine the second global derivatives w.r.t x/y/z at point r,s,t
C are calculated.
C
C \param elecord()()()   real*8   (i) vector containing nodal coordinates
C \param xjm()()()       real*8   (i) jacobian matrix
C \param bm_in()()()     real*8   (i) working array
C \param xder2_in()()()  real*8   (i) working array
C \param derxy()()()     real*8   (i) global derivatives
C \param derxy2()()()    real*8   (o) 2nd global derivatives
C \param sizevec(6)      integer  (i) some sizes
C
C \return void
C
C \author mn
C \date   10/04
C
C--------------------------------------------------------------------------
      subroutine f3fgder2(elecord_f, xjm, bm_in, xder2_in, derxy,
     &          derxy2, deriv2, sizevec)

      implicit none

      integer sizevec(6)
      real*8 elecord_f(sizevec(4),sizevec(1),3)
      real*8 xjm(sizevec(4),3,3)
      real*8 bm_in(sizevec(4),32,32)
      real*8 xder2_in(sizevec(4),32,32)
      real*8 derxy(sizevec(4),sizevec(1),3)
      real*8 derxy2(sizevec(4),sizevec(1),6)
      real*8 deriv2(sizevec(1),6)

      integer i,j,value,k
      real*8 r0,r1,r2,r3,r4,r5
      real*8 TWO
      real*8 bm(6,6), xder2(3,6)


c     calculate elements of jacobian_bar matrix
      TWO=2.0

      do k=1,sizevec(5)

        bm(1,1) = xjm(k,1,1)*xjm(k,1,1)
        bm(1,2) = xjm(k,1,2)*xjm(k,1,2)
        bm(1,3) = xjm(k,1,3)*xjm(k,1,3)
        bm(1,4) = xjm(k,1,1)*xjm(k,1,2)
        bm(1,5) = xjm(k,1,1)*xjm(k,1,3)
        bm(1,6) = xjm(k,1,2)*xjm(k,1,3)

        bm(2,1) = xjm(k,2,1)*xjm(k,2,1)
        bm(2,2) = xjm(k,2,2)*xjm(k,2,2)
        bm(2,3) = xjm(k,2,3)*xjm(k,2,3)
        bm(2,4) = xjm(k,2,1)*xjm(k,2,2)
        bm(2,5) = xjm(k,2,1)*xjm(k,2,3)
        bm(2,6) = xjm(k,2,2)*xjm(k,2,3)

        bm(3,1) = xjm(k,3,1)*xjm(k,3,1)
        bm(3,2) = xjm(k,3,2)*xjm(k,3,2)
        bm(3,3) = xjm(k,3,3)*xjm(k,3,3)
        bm(3,4) = xjm(k,3,1)*xjm(k,3,2)
        bm(3,5) = xjm(k,3,1)*xjm(k,3,3)
        bm(3,6) = xjm(k,3,2)*xjm(k,3,3)

        bm(4,1) = TWO*xjm(k,1,1)*xjm(k,2,1)
        bm(4,2) = TWO*xjm(k,1,2)*xjm(k,2,2)
        bm(4,3) = TWO*xjm(k,1,3)*xjm(k,2,3)
        bm(4,4) =     xjm(k,1,1)*xjm(k,2,2)+xjm(k,1,2)*xjm(k,2,1)
        bm(4,5) =     xjm(k,1,1)*xjm(k,2,3)+xjm(k,1,3)*xjm(k,2,1)
        bm(4,6) =     xjm(k,1,2)*xjm(k,2,3)+xjm(k,1,3)*xjm(k,2,2)

        bm(5,1) = TWO*xjm(k,1,1)*xjm(k,3,1)
        bm(5,2) = TWO*xjm(k,1,2)*xjm(k,3,2)
        bm(5,3) = TWO*xjm(k,1,3)*xjm(k,3,3)
        bm(5,4) =     xjm(k,1,1)*xjm(k,3,2)+xjm(k,1,2)*xjm(k,3,1)
        bm(5,5) =     xjm(k,1,1)*xjm(k,3,3)+xjm(k,1,3)*xjm(k,3,1)
        bm(5,6) =     xjm(k,1,2)*xjm(k,3,3)+xjm(k,1,3)*xjm(k,3,2)

        bm(6,1) = TWO*xjm(k,2,1)*xjm(k,3,1)
        bm(6,2) = TWO*xjm(k,2,2)*xjm(k,3,2)
        bm(6,3) = TWO*xjm(k,2,3)*xjm(k,3,3)
        bm(6,4) =     xjm(k,2,1)*xjm(k,3,2)+xjm(k,2,2)*xjm(k,3,1)
        bm(6,5) =     xjm(k,2,1)*xjm(k,3,3)+xjm(k,2,3)*xjm(k,3,1)
        bm(6,6) =     xjm(k,2,2)*xjm(k,3,3)+xjm(k,2,3)*xjm(k,3,2)


c       inverse of jacobian_bar matrix
        call f3finv6(bm, value)


c       initialise
        do i=1,3
          do j=1,6
            xder2(i,j)=0.0
          enddo
        enddo

        do i=1,sizevec(2)
          do j=1,6
            derxy2(k,i,j)=0.0
          enddo
        enddo

c       determine 2nd derivatives of coord.-functions
        do i=1,sizevec(2)
          xder2(1,1) =xder2(1,1)+deriv2(i,1) * elecord_f(k,i,1)
          xder2(1,2) =xder2(1,2)+deriv2(i,2) * elecord_f(k,i,1)
          xder2(1,3) =xder2(1,3)+deriv2(i,3) * elecord_f(k,i,1)
          xder2(1,4) =xder2(1,4)+deriv2(i,4) * elecord_f(k,i,1)
          xder2(1,5) =xder2(1,5)+deriv2(i,5) * elecord_f(k,i,1)
          xder2(1,6) =xder2(1,6)+deriv2(i,6) * elecord_f(k,i,1)

          xder2(2,1) =xder2(2,1)+deriv2(i,1) * elecord_f(k,i,2)
          xder2(2,2) =xder2(2,2)+deriv2(i,2) * elecord_f(k,i,2)
          xder2(2,3) =xder2(2,3)+deriv2(i,3) * elecord_f(k,i,2)
          xder2(2,4) =xder2(2,4)+deriv2(i,4) * elecord_f(k,i,2)
          xder2(2,5) =xder2(2,5)+deriv2(i,5) * elecord_f(k,i,2)
          xder2(2,6) =xder2(2,6)+deriv2(i,6) * elecord_f(k,i,2)

          xder2(3,1) =xder2(3,1)+deriv2(i,1) * elecord_f(k,i,3)
          xder2(3,2) =xder2(3,2)+deriv2(i,2) * elecord_f(k,i,3)
          xder2(3,3) =xder2(3,3)+deriv2(i,3) * elecord_f(k,i,3)
          xder2(3,4) =xder2(3,4)+deriv2(i,4) * elecord_f(k,i,3)
          xder2(3,5) =xder2(3,5)+deriv2(i,5) * elecord_f(k,i,3)
          xder2(3,6) =xder2(3,6)+deriv2(i,6) * elecord_f(k,i,3)
        enddo

c       calculate second global derivatives
        do i=1,sizevec(2)
          r0 = deriv2(i,1) - xder2(1,1)*derxy(k,i,1)
     &                   - xder2(2,1)*derxy(k,i,2)
     &                   - xder2(3,1)*derxy(k,i,3)
          r1 = deriv2(i,2) - xder2(1,2)*derxy(k,i,1)
     &                   - xder2(2,2)*derxy(k,i,2)
     &                   - xder2(3,2)*derxy(k,i,3)
          r2 = deriv2(i,3) - xder2(1,3)*derxy(k,i,1)
     &                   - xder2(2,3)*derxy(k,i,2)
     &                   - xder2(3,3)*derxy(k,i,3)
          r3 = deriv2(i,4) - xder2(1,4)*derxy(k,i,1)
     &                   - xder2(2,4)*derxy(k,i,2)
     &                   - xder2(3,4)*derxy(k,i,3)
          r4 = deriv2(i,5) - xder2(1,5)*derxy(k,i,1)
     &                   - xder2(2,5)*derxy(k,i,2)
     &                   - xder2(3,5)*derxy(k,i,3)
          r5 = deriv2(i,6) - xder2(1,6)*derxy(k,i,1)
     &                   - xder2(2,6)*derxy(k,i,2)
     &                   - xder2(3,6)*derxy(k,i,3)

          derxy2(k,i,1) = derxy2(k,i,1)+bm(1,1)*r0 + bm(2,1)*r1
     &          + bm(3,1)*r2 +bm(4,1)*r3 + bm(5,1)*r4 + bm(6,1)*r5
          derxy2(k,i,2) = derxy2(k,i,2)+bm(1,2)*r0 + bm(2,2)*r1
     &          + bm(3,2)*r2 +bm(4,2)*r3 + bm(5,2)*r4 + bm(6,2)*r5
          derxy2(k,i,3) = derxy2(k,i,3)+bm(1,3)*r0 + bm(2,3)*r1
     &          + bm(3,3)*r2 +bm(4,3)*r3 + bm(5,3)*r4 + bm(6,3)*r5
          derxy2(k,i,4) = derxy2(k,i,4)+bm(1,4)*r0 + bm(2,4)*r1
     &          + bm(3,4)*r2 +bm(4,4)*r3 + bm(5,4)*r4 + bm(6,4)*r5
          derxy2(k,i,5) = derxy2(k,i,5)+bm(1,5)*r0 + bm(2,5)*r1
     &          + bm(3,5)*r2 +bm(4,5)*r3 + bm(5,5)*r4 + bm(6,5)*r5
          derxy2(k,i,6) = derxy2(k,i,6)+bm(1,6)*r0 + bm(2,6)*r1
     &          + bm(3,6)*r2 +bm(4,6)*r3 + bm(5,6)*r4 + bm(6,6)*r5
        enddo


      enddo

      return
      end subroutine



