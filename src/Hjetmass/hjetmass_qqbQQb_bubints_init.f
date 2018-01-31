      subroutine hjetmass_qqbQQb_bubints_init
     &   (i1,i2,i3,i4,bubints,ord,s)
        use mod_qcdloop_c
          implicit none
          include 'types.f'
          include 'mxpart.f'
          include 'HqqbQQblabels.f'
          include 'constants.f'
          real*16 s(mxpart,mxpart)
          include 'scale.f'
          include 'masses.f'
          integer i1,i2,i3,i4
          integer ord
          complex*32 bubints(3)
          real*16 t
          real*16 mt2, mhsq, qmusq

        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        mt2 = mt**2
        mhsq = s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)
        qmusq = musq

        bubints(:) = (0d0,0d0)

c   "(mhsq;mtsq,mtsq)",
        call qlI2q(mhsq, mt2, mt2, qmusq, ord, bubints(1))

c   "( s34;mtsq,mtsq)"]
        call qlI2q(s(i3,i4), mt2, mt2, qmusq, ord, bubints(2))

c   "( s12;mtsq,mtsq)",
        call qlI2q(s(i1,i2), mt2, mt2, qmusq, ord, bubints(3))

      end subroutine

