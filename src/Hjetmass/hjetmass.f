      subroutine hjetmass(p,msq)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'constants.f'
      include 'masses.f'
      include 'hdecaymode.f'
      double precision s(mxpart,mxpart)
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'asymptotic.f'
      include 'first.f'
      integer j,k
      double precision msq(-nf:nf,-nf:nf),p(mxpart,4),gg,qg,gq,qq,hdecay
      double precision msq_tmp(-nf:nf, -nf:nf)
      double precision ehsvm3,ehsvm4,s34
      double precision msqgamgam
      double precision origmbsq
      double precision virtme_gg
      double precision gg_htl, gq_htl, qg_htl, qq_htl
      double precision dot

      double precision sman, tman, uman

      ! sushi amplitude defs
      double precision cex_prefactor, aex_prefactor
      double precision c1ex, c2ex, c3ex, c4ex
      double precision c1exSq, c2exSq, c3exSq, c4exSq
      double precision aexSq, fac, Asq
      double precision asmh_real, asmh_abs
      double precision c1smh_real, c2smh_real
      double precision c1smh_abs, c2smh_abs
      double complex c1smh, c2smh, asmh
      double precision sushi_gg, sushi_gq, sushi_qg, sushi_qq
      ! end sushi amplitude defs
      integer, parameter :: iglue = 5

      call dotem(iglue,p,s)

      s34=(p(3,4)+p(4,4))**2
     & -(p(3,1)+p(4,1))**2-(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2

C   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          call htautaudecay(p,3,4,hdecay)
      elseif (hdecaymode == 'bqba') then
          call hbbdecay(p,3,4,hdecay)
      elseif (hdecaymode == 'gaga') then
          hdecay=msqgamgam(sqrt(s34))
      else
      write(6,*) 'Unimplemented process in hjetmass'
      stop
      endif

      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)


      sman = s(1,2)
      tman = s(1,5)
      uman = s(2,5)

      !! amplitude from SusHi
      cex_prefactor = sqrt(gsq)**3 / 4d0 / pi**2 / sqrt(vevsq)

      if (first .eqv. .true.) then
        call sushi_bernini(18)
      else
        first =  .false. 
      endif

      c1exSq = (c1smh_abs(sman,tman,uman,mt**2) *
     & cex_prefactor * mt**2 / 16d0)**2

      c2exSq = (c2smh_abs(sman,tman,uman,mt**2) *
     & cex_prefactor * mt**2 / 16d0)**2

      c3exSq = (c2smh_abs(tman,sman,uman,mt**2) *
     & cex_prefactor * mt**2 / 16d0)**2

      c4exSq = (c2smh_abs(uman,sman,tman,mt**2) *
     & cex_prefactor * mt**2 / 16d0)**2

      sushi_gg = (c1exSq + c2exSq + c3exSq + c4exSq)*
     &              24*16/sman/tman/uman
      gg =  sushi_gg*hdecay/256d0

      aex_prefactor = sqrt(gsq)**3 / 32d0 / pi**2 / sqrt(vevsq)
      aexSq = (asmh_abs(tman,sman,uman,mt**2) / tman * aex_prefactor *
     &             mt**2 * 8d0/3d0)**2
      sushi_qg = aexSq*8d0/2d0 * (-uman**2 * tman - sman**2 * tman)
      qg =  sushi_qg*hdecay/96d0

      aexSq = (asmh_abs(uman,sman,tman,mt**2) / uman * aex_prefactor *
     &             mt**2 * 8d0/3d0)**2
      sushi_gq = aexSq*8d0/2d0 * (-tman**2 * uman - sman**2 * uman)
      gq = sushi_gq*hdecay/96d0


      aexSq = (asmh_abs(sman,uman,tman,mt**2) / sman * aex_prefactor *
     &             mt**2 * 8d0/3d0)**2
      sushi_qq = aexSq*8d0/2d0 * (-tman**2 * sman - uman**2 * sman)
      qq = sushi_qq*hdecay/(-36d0)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0

      if ((j.eq. 0) .or. (k.eq.0)) then
           if ((j.eq. 0) .and. (k.eq.0)) then
                msq(j,k)=gg
           elseif ((j.eq.0).and.(k.ne.0)) then
                msq(j,k)=gq
           elseif ((j.ne.0).and.(k.eq.0)) then
                msq(j,k)=qg
           endif
      elseif ((j.eq.-k).and. (j.ne.0)) then
           msq(j,k)=qq
      endif

      enddo
      enddo

      return
      end

