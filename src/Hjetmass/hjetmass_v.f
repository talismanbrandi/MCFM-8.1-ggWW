      subroutine hjetmass_v(p,msq)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'constants.f'
      include 'masses.f'
      include 'hdecaymode.f'
      include 'sprods_com.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'scheme.f'
      include 'epinv.f'
      include 'scale.f'
      include 'asymptotic.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),gg,qg,gq,hdecay
      double precision gg_new, gq_old, gg_mix
      double precision gq_new, gq_mix
      double precision s34, sman, tman, uman
      double precision msqgamgam
      double precision virtme_gg, virtme_gq, qa, aq
      double precision virtme_gq_test
      double precision virtme_gq_new, virtme_gg_multi
      double precision mesq_gg, mesq_gq
      double precision msqlo(-nf:nf, -nf:nf)
      double precision msqv(-nf:nf, - nf:nf)
      double precision virtgg, virtqa, virtaq, virtqg, virtgq, Asq, fac
      double precision gg_insert, gq_insert, test, dummy
      double precision x1000, x100, x10, ep2, ep1, fin, epsave
      double precision ep2_mcfm, ep1_mcfm, fin_mcfm
      double precision insarray(-nf:nf,-nf:nf)
      double precision insarray_new(-nf:nf,-nf:nf)
      s34=(p(3,4)+p(4,4))**2
     & -(p(3,1)+p(4,1))**2-(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2

      scheme = 'tH-V'

C   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          call htautaudecay(p,3,4,hdecay)
      elseif (hdecaymode == 'bqba') then
          call hbbdecay(p,3,4,hdecay)
      elseif (hdecaymode == 'gaga') then
          hdecay=msqgamgam(hmass)
      else
      write(6,*) 'Unimplemented process in hjetmass'
      stop
      endif

      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)

      sman = s(1,2)
      tman = s(1,5)
      uman = s(2,5)

      gg = virtme_gg_multi(sman,tman,uman,mtex) * hdecay
      gq = virtme_gq_new(uman,tman,sman,mtex) * hdecay
      qg = virtme_gq_new(tman,uman,sman,mtex) * hdecay
      aq = -8d0/3d0 * virtme_gq_new(sman,uman,tman,mtex) * hdecay
      qa = -8d0/3d0 * virtme_gq_new(sman,tman,uman,mtex) * hdecay

      call hjetmass(p,msqlo)
      call insert_exact_new(p,insarray_new)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0 

      if ((j.eq.0).and.(k.eq.0)) then
        msq(j,k) = (gg - insarray_new(0,0))
      endif
      if ((j.gt.0).and.(k.eq.-j)) then
        msq(j,k) = (qa - insarray_new(j,k))
      endif
      if ((j.lt.0).and.(k.eq.-j)) then
        msq(j,k) = (aq - insarray_new(j,k))
      endif
      if ((j.eq.0).and.(k.ne.0)) then
        msq(j,k) = (gq - insarray_new(j,k))
      endif
      if ((j.ne.0).and.(k.eq.0)) then
        msq(j,k) = (qg - insarray_new(j,k))
      endif

      enddo
      enddo

      return
      end

      subroutine insert_exact_new(p, msqlo)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        include 'nf.f'
        include 'constants.f'
        include 'masses.f'
        include 'scale.f'
        include 'epinv.f'
        include 'qcdcouple.f'
        include 'ewcouple.f'
        include 'sprods_com.f'
        double precision, intent(in) :: p(mxpart,4)
        double precision, intent(out) :: msqlo(-nf:nf, -nf:nf)
        double precision sman, tman, uman
        double complex insert_gg
        double complex cdlogwrap
        double precision kg, mh, nl, pcut, insgq_new, bf
        double complex LogAbsMus, LogAbsMut, LogAbsMuu
        integer j,k

        sman = s(1,2)
        tman = s(1,5)
        uman = s(2,5)

        mh = hmass

        call hjetmass(p,msqlo)

        nl = 5d0
        pcut = 1d0
        bf = (11d0/6d0*ca - 2d0/3d0*tr*nl)
        kg = ca*(67d0/18d0 - pisqo6) - 10d0/9d0 * tr*nl
     &          - ca*dlog(pcut)**2
     &          + (11d0/6d0*ca - 2d0/3d0*tr*nl)*(pcut - 1 - dlog(pcut))

        LogAbsMus = cdlogwrap(-musq/sman)
        LogAbsMut = cdlogwrap(-musq/tman)
        LogAbsMuu = cdlogwrap(-musq/uman)

      insert_gg =
     & 3.D0*Kg + 3.D0*bf+ 3.D0*bf*epinv + bf*LogAbsMus + 
     & bf*LogAbsMut + bf*LogAbsMuu - cA*pi**2 + 3.D0*cA*epinv**2
     &  + cA*LogAbsMus*epinv + 1.D0/2.D0*cA*LogAbsMus**2 + cA*LogAbsMut
     & *epinv + 1.D0/2.D0*cA*LogAbsMut**2 + cA*LogAbsMuu*epinv + 1.D0/2.
     & D0*cA*LogAbsMuu**2

      do j=-nf,nf
      do k=-nf,nf

      if ((j.eq.0).and.(k.eq.0)) then
        msqlo(j,k) = msqlo(j,k)*dreal(insert_gg)
      endif
      if ((j.gt.0).and.(k.eq.-j)) then
        msqlo(j,k) = msqlo(j,k) * insgq_new(sman,tman,uman)
      endif
      if ((j.lt.0).and.(k.eq.-j)) then
        msqlo(j,k) = msqlo(j,k) * insgq_new(sman,uman,tman)
      endif
      if ((j.eq.0).and.(k.ne.0)) then
        msqlo(j,k) = msqlo(j,k) * insgq_new(uman,tman,sman)
      endif
      if ((j.ne.0).and.(k.eq.0)) then
        msqlo(j,k) = msqlo(j,k) * insgq_new(tman,uman,sman)
      endif

      ! convert prefactor from lo to virt
      msqlo(j,k) = msqlo(j,k) * as/2/pi

      enddo
      enddo
    
      end subroutine insert_exact_new

      double precision function insgq_new(sman, tman, uman)
        implicit none
        include 'types.f'
        include 'constants.f'
        include 'masses.f'
        include 'scale.f'
        include 'epinv.f'
        include 'qcdcouple.f'
        include 'ewcouple.f'
        double precision kg, kq
        double precision, intent(in) :: sman, tman, uman
        double precision kquark, kgluon, pcut, nl

        double complex LogMums, LogMumt, LogMumu
        double complex insop
        double complex cdlogwrap


        nl = 5d0
        pcut = 1d0
        kg = ca*(67d0/18d0 - pisqo6) - 10d0/9d0 * tr*nl
     &          - ca*dlog(pcut)**2
     &          + (11d0/6d0*ca - 2d0/3d0*tr*nl)*(pcut - 1 - dlog(pcut))
        kq = cf*(7d0/2d0 - pisqo6)
     &          - cf*dlog(pcut)**2
     &          + 3d0/2d0*cf*(pcut - 1 - dlog(pcut))

        kquark = kq
        kgluon = kg

        LogMumt = cdlogwrap(-musq/tman)
        LogMumu = cdlogwrap(-musq/uman)
        LogMums = cdlogwrap(-musq/sman)

      insop =
     &  + epinv * ( 3.D0*CF + 2.D0*CF*LogMums - 2.D0/3.D0*nl*tr + 11.D0/
     &    6.D0*cA + cA*LogMumt + cA*LogMumu - cA*LogMums )
      insop = insop + epinv**2 * ( 2.D0*CF + cA )
      insop = insop + 2.D0*Kquark + Kgluon + 3.D0*CF + 3.D0*CF*LogMums
     &     + CF*LogMums**2 - 2.D0/3.D0*CF*pi**2 - 2.D0/3.D0*nl*tr - 1.D0
     &    /3.D0*nl*tr*LogMumt - 1.D0/3.D0*nl*tr*LogMumu + 11.D0/6.D0*cA
     &     + 5.D0/3.D0*cA*LogMumt + 1.D0/2.D0*cA*LogMumt**2 + 5.D0/3.D0
     &    *cA*LogMumu + 1.D0/2.D0*cA*LogMumu**2 - 3.D0/2.D0*cA*LogMums
     &     - 1.D0/2.D0*cA*LogMums**2 - 1.D0/3.D0*cA*pi**2

      insgq_new = dreal(insop)

      end function insgq_new

