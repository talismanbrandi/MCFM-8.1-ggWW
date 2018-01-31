      subroutine hjetmass_bubints_init(i1,i2,i3,i4,bubints,ord,s)
          use mod_qcdloop_c
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'scale.f'
          include 'masses.f'
          integer i1,i2,i3,i4
          integer ord
          complex*32 bubints(9)
          real*16 mhsq_q, qmt2
          real*16 t, s(mxpart,mxpart)
          real*16 mt2, mhsq, qmusq

        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        mt2 = mt**2
        mhsq = s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)
        qmusq = musq

        bubints(:) = (0q0,0q0)

!           (mhsq;mtsq,mtsq)
        call qlI2q(mhsq, mt2, mt2, qmusq, ord, bubints(1))

!           (s123;mtsq,mtsq)
        call qlI2q(t(i1,i2,i3), mt2, mt2, qmusq, ord, bubints(2))

!           (s124;mtsq,mtsq)
        call qlI2q(t(i1,i2,i4), mt2, mt2, qmusq, ord, bubints(3))

!           ( s12;mtsq,mtsq)
        call qlI2q(s(i1,i2), mt2, mt2, qmusq, ord, bubints(4))

!           (s134;mtsq,mtsq)
        call qlI2q(t(i1,i3,i4), mt2, mt2, qmusq, ord, bubints(5))

!           ( s14;mtsq,mtsq)
        call qlI2q(s(i1,i4), mt2, mt2, qmusq, ord, bubints(6))

!           (s234;mtsq,mtsq)
        call qlI2q(t(i2,i3,i4), mt2, mt2, qmusq, ord, bubints(7))

!           ( s23;mtsq,mtsq)
        call qlI2q(s(i2,i3), mt2, mt2, qmusq, ord, bubints(8))

!           ( s34;mtsq,mtsq)
        call qlI2q(s(i3,i4), mt2, mt2, qmusq, ord, bubints(9))

      end subroutine
