      subroutine hjetmass_boxints_init(i1,i2,i3,i4,boxints,ord,s)
        use mod_qcdloop_c
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'scale.f'
          include 'masses.f'

          integer i1,i2,i3,i4
          integer ord
          complex*32 boxints(16)
          real*16 t, s(mxpart,mxpart)
          real*16 mt2, mhsq, qmusq

        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        mt2 = mt**2
        mhsq = s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)
        qmusq = musq

        boxints(:) = (0q0,0q0)

!      (   0,   0,   0,s123, s12, s23;mts,mts,mts,mts)
        call qlI4q(0q0, 0q0, 0q0, t(i1,i2,i3), s(i1,i2), s(i2,i3),
     &          mt2, mt2, mt2, mt2, qmusq, ord, boxints(1))

!      (   0,   0,   0,s124, s14, s12;mts,mts,mts,mts)
        call qlI4q(0q0, 0q0, 0q0, t(i1,i2,i4), s(i1,i4), s(i1,i2),
     &          mt2, mt2, mt2, mt2, qmusq, ord, boxints(2))

!      (   0,   0,   0,s134, s34, s14;mts,mts,mts,mts)
        call qlI4q(0q0, 0q0, 0q0, t(i1,i3,i4), s(i3,i4), s(i1,i4),
     &          mt2, mt2, mt2, mt2, qmusq, ord, boxints(3))

!      (   0,   0,   0,s234, s23, s34;mts,mts,mts,mts)
        call qlI4q(0q0, 0q0, 0q0, t(i2,i3,i4), s(i2,i3), s(i3,i4),
     &          mt2, mt2, mt2, mt2, qmusq, ord, boxints(4))

!      (   0,   0, s12,mhsq, s34,s123;mts,mts,mts,mts)
        call qlI4q(0q0, 0q0, s(i1,i2), mhsq, s(i3,i4), t(i1,i2,i3),
     &          mt2, mt2, mt2, mt2, qmusq, ord, boxints(5))

!      (   0,   0, s12,mhsq, s34,s124;mts,mts,mts,mts)
        call qlI4q(0q0, 0q0, s(i1,i2), mhsq, s(i3,i4), t(i1,i2,i4),
     &          mt2, mt2, mt2, mt2, qmusq, ord, boxints(6))

!      (   0,   0, s14,mhsq, s23,s124;mts,mts,mts,mts)
        call qlI4q(0q0, 0q0, s(i1,i4), mhsq, s(i2,i3), t(i1,i2,i4),
     &          mt2, mt2, mt2, mt2, qmusq, ord, boxints(7))

!      (   0,   0, s14,mhsq, s23,s134;mts,mts,mts,mts)
        call qlI4q(0q0, 0q0, s(i1,i4), mhsq, s(i2,i3), t(i1,i3,i4),
     &          mt2, mt2, mt2, mt2, qmusq, ord, boxints(8))

!      (   0,   0, s23,mhsq, s14,s123;mts,mts,mts,mts)
        call qlI4q(0q0, 0q0, s(i2,i3), mhsq, s(i1,i4), t(i1,i2,i3),
     &          mt2, mt2, mt2, mt2, qmusq, ord, boxints(9))

!      (   0,   0, s23,mhsq, s14,s234;mts,mts,mts,mts)
        call qlI4q(0q0, 0q0, s(i2,i3), mhsq, s(i1,i4), t(i2,i3,i4),
     &          mt2, mt2, mt2, mt2, qmusq, ord, boxints(10))

!      (   0,   0, s34,mhsq, s12,s134;mts,mts,mts,mts)
        call qlI4q(0q0, 0q0, s(i3,i4), mhsq, s(i1,i2), t(i1,i3,i4),
     &          mt2, mt2, mt2, mt2, qmusq, ord, boxints(11))

!      (   0,   0, s34,mhsq, s12,s234;mts,mts,mts,mts)
        call qlI4q(0q0, 0q0, s(i3,i4), mhsq, s(i1,i2), t(i2,i3,i4),
     &          mt2, mt2, mt2, mt2, qmusq, ord, boxints(12))

!      (   0, s12,   0,mhsq,s123,s124;mts,mts,mts,mts)
        call qlI4q(0q0, s(i1,i2), 0q0, mhsq, t(i1,i2,i3), t(i1,i2,i4),
     &          mt2, mt2, mt2, mt2, qmusq, ord, boxints(13))

!      (   0, s14,   0,mhsq,s124,s134;mts,mts,mts,mts)
        call qlI4q(0q0, s(i1,i4), 0q0, mhsq, t(i1,i2,i4), t(i1,i3,i4),
     &          mt2, mt2, mt2, mt2, qmusq, ord, boxints(14))

!      (   0, s23,   0,mhsq,s123,s234;mts,mts,mts,mts)
        call qlI4q(0q0, s(i2,i3), 0q0, mhsq, t(i1,i2,i3), t(i2,i3,i4),
     &          mt2, mt2, mt2, mt2, qmusq, ord, boxints(15))

!      (   0, s34,   0,mhsq,s134,s234;mts,mts,mts,mts)
        call qlI4q(0q0, s(i3,i4), 0q0, mhsq, t(i1,i3,i4), t(i2,i3,i4),
     &          mt2, mt2, mt2, mt2, qmusq, ord, boxints(16))

      end subroutine
