      subroutine hjetmass_triints_init(i1,i2,i3,i4,triints,ord,s)
        use mod_qcdloop_c
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'scale.f'
          include 'masses.f'
          integer i1,i2,i3,i4
          integer ord
          complex*32 triints(18)
          real*16 t, s(mxpart,mxpart)
          real*16 mt2, mhsq, qmusq

        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        mt2 = mt**2
        mhsq = s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)
        qmusq = musq

        triints(:) = (0q0,0q0)

!         ( s12,   0,   0;mtsq,mtsq,mtsq)
        call qlI3q(s(i1,i2), 0q0, 0q0, mt2, mt2, mt2, qmusq, ord,
     &  triints(1))

!         (s123,   0,mhsq;mtsq,mtsq,mtsq)
        call qlI3q(t(i1,i2,i3), 0q0, mhsq, mt2, mt2, mt2, qmusq, ord,
     &  triints(2))

!         (s124,   0,mhsq;mtsq,mtsq,mtsq)
        call qlI3q(t(i1,i2,i4), 0q0, mhsq, mt2, mt2, mt2, qmusq, ord,
     &  triints(3))

!         ( s12,s123,   0;mtsq,mtsq,mtsq)
        call qlI3q(s(i1,i2), t(i1,i2,i3), 0q0, mt2, mt2, mt2, qmusq,
     &  ord, triints(4))

!         ( s12,s124,   0;mtsq,mtsq,mtsq)
        call qlI3q(s(i1,i2), t(i1,i2,i4), 0q0, mt2, mt2, mt2, qmusq,
     &  ord, triints(5))

!         (s134,   0,mhsq;mtsq,mtsq,mtsq)
        call qlI3q(t(i1,i3,i4), 0q0, mhsq, mt2, mt2, mt2, qmusq, ord,
     &  triints(6))

!         ( s14,   0,   0;mtsq,mtsq,mtsq)
        call qlI3q(s(i1,i4), 0q0, 0q0, mt2, mt2, mt2, qmusq, ord,
     &  triints(7))

!         ( s14,s124,   0;mtsq,mtsq,mtsq)
        call qlI3q(s(i1,i4), t(i1,i2,i4), 0q0, mt2, mt2, mt2, qmusq,
     &  ord, triints(8))

!         ( s14,s134,   0;mtsq,mtsq,mtsq)
        call qlI3q(s(i1,i4), t(i1,i3,i4), 0q0, mt2, mt2, mt2, qmusq,
     &  ord, triints(9))

!         ( s23,   0,   0;mtsq,mtsq,mtsq)
        call qlI3q(s(i2,i3), 0q0, 0q0, mt2, mt2, mt2, qmusq, ord,
     &  triints(10))

!         (s234,   0,mhsq;mtsq,mtsq,mtsq)
        call qlI3q(t(i2,i3,i4), 0q0, mhsq, mt2, mt2, mt2, qmusq, ord,
     &  triints(11))

!         ( s23,mhsq, s14;mtsq,mtsq,mtsq)
        call qlI3q(s(i2,i3), mhsq, s(i1,i4), mt2, mt2, mt2, qmusq, ord,
     &  triints(12))

!         ( s23,s123,   0;mtsq,mtsq,mtsq)
        call qlI3q(s(i2,i3), t(i1,i2,i3), 0q0, mt2, mt2, mt2, qmusq,
     &  ord, triints(13))

!         ( s23,s234,   0;mtsq,mtsq,mtsq)
        call qlI3q(s(i2,i3), t(i2,i3,i4), 0q0, mt2, mt2, mt2, qmusq,
     &  ord, triints(14))

!         ( s34,   0,   0;mtsq,mtsq,mtsq)
        call qlI3q(s(i3,i4), 0q0, 0q0, mt2, mt2, mt2, qmusq, ord,
     &  triints(15))

!         ( s34,mhsq, s12;mtsq,mtsq,mtsq)
        call qlI3q(s(i3,i4), mhsq, s(i1,i2), mt2, mt2, mt2, qmusq, ord,
     &  triints(16))

!         ( s34,s134,   0;mtsq,mtsq,mtsq)
        call qlI3q(s(i3,i4), t(i1,i3,i4), 0q0, mt2, mt2, mt2, qmusq,
     &  ord, triints(17))

!         ( s34,s234,   0;mtsq,mtsq,mtsq)
        call qlI3q(s(i3,i4), t(i2,i3,i4), 0q0, mt2, mt2, mt2, qmusq,
     &  ord, triints(18))

      end subroutine
