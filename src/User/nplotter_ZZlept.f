      subroutine nplotter_ZZlept(p,wt,wt2,switch)
      implicit none
      include 'types.f'
c--- Variable passed in to this routine:
c
c---      p:  4-momenta of particles in the format p(i,4)
c---          with the particles numbered according to the input file
c---          and components labelled by (px,py,pz,E)
c
c---     wt:  weight of this event
c
c---    wt2:  weight^2 of this event
c
c--- switch:  an integer:: equal to 0 or 1, depending on the type of event
c---                0  --> lowest order, virtual or real radiation
c---                1  --> counterterm for real radiation
      
      include 'vegas_common.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'histo.f'
      include 'outputflags.f'
      include 'interference.f'
      real(dp):: p(mxpart,4),wt,wt2,m3456,pt34,pttwo
c--- Modifications for adding a discriminator for ggZZ -- MS
      real(dp):: D_MELA
c--- End Modification -- MS
      integer:: switch,n,nplotmax
      integer tag
      logical, save::first=.true.
      common/nplotmax/nplotmax
ccccc!$omp threadprivate(first,/nplotmax/)

************************************************************************
*                                                                      *
*     INITIAL BOOKKEEPING                                              *
*                                                                      *
************************************************************************
      
    
      
      if (first) then
c--- Initialize histograms, without computing any quantities; instead
c--- set them to dummy values
        tag=tagbook
        goto 99
      else
c--- Add event in histograms
        tag=tagplot
      endif

************************************************************************
*                                                                      *
*     DEFINITIONS OF QUANTITIES TO PLOT                                *
*                                                                      *
************************************************************************

      m3456=(p(3,4)+p(4,4)+p(5,4)+p(6,4))**2
     &     -(p(3,1)+p(4,1)+p(5,1)+p(6,1))**2
     &     -(p(3,2)+p(4,2)+p(5,2)+p(6,2))**2
     &     -(p(3,3)+p(4,3)+p(5,3)+p(6,3))**2
      m3456=sqrt(max(m3456,zip))
      
      pt34=pttwo(3,4,p)
       
************************************************************************
*                                                                      *
*     FILL HISTOGRAMS                                                  *
*                                                                      *
************************************************************************

c--- Call histogram routines
   99 continue

c--- Book and fill ntuple if that option is set, remembering to divide
c--- by # of iterations now that is handled at end for regular histograms
      if (creatent .eqv. .true.) then
        call bookfill(tag,p,wt/real(itmx,dp))  
      endif

c--- "n" will count the number of histograms
      n=nextnplot

c--- Syntax of "bookplot" routine is:
c
c---   call bookplot(n,tag,titlex,var,wt,wt2,xmin,xmax,dx,llplot)
c
c---        n:  internal number of histogram
c---      tag:  "book" to initialize histogram, "plot" to fill
c---   titlex:  title of histogram
c---      var:  value of quantity being plotted
c---       wt:  weight of this event (passed in)
c---      wt2:  weight of this event (passed in)
c---     xmin:  lowest value to bin
c---     xmax:  highest value to bin
c---       dx:  bin width
c---   llplot:  equal to "lin"/"log" for linear/log scale

c--- Plots of m(3456) in specific regions
c--- Modifications for adding a few relevant extended histograms for ggZZ -- AP
      call bookplot(n,tag,'60 < m(3456) < 2010',
     & m3456,wt,wt2,60._dp,2010._dp,10._dp,'log')
      n=n+1

      call bookplot(n,tag,'250 < m(3456) < 4000',
     & m3456,wt,wt2,250._dp,4000._dp,20._dp,'log')
      n=n+1

      call bookplot(n,tag,'250 < m(3456) < 10000',
     & m3456,wt,wt2,250._dp,10000._dp,50._dp,'log')
      n=n+1
c--- End Modification -- AP
      call bookplot(n,tag,'10 < m(3456) < 2010',
     & m3456,wt,wt2,10._dp,2010._dp,20._dp,'log')
      n=n+1
      
      call bookplot(n,tag,'130 < m(3456) < 2010',
     & m3456,wt,wt2,130._dp,2010._dp,20._dp,'log')
      n=n+1
           
      call bookplot(n,tag,'300 < m(3456) < 2020',
     & m3456,wt,wt2,300._dp,2020._dp,20._dp,'log')
      n=n+1
      
      call bookplot(n,tag,'10 < m(3456) < 130',
     & m3456,wt,wt2,10._dp,130._dp,5._dp,'lin')
      n=n+1
      
      call bookplot(n,tag,'pt(Z)',
     & pt34,wt,wt2,0._dp,2._dp,0.02_dp,'lin')
      n=n+1
      
      call bookplot(n,tag,'+INTEGRAL+ pt(Z)',
     & pt34,wt,wt2,0._dp,10._dp,0.1_dp,'lin')
      n=n+1
      
      call bookplot(n,tag,'50 < m(3456) < 250',
     & m3456,wt,wt2,50._dp,250._dp,2._dp,'log')
      n=n+1

c--- Modifications for adding a few relevant extended histograms for ggZZ Discriminator -- MS & AP      
      call get_MELA_Discr_ggZZ(p,D_MELA)
      call bookplot(n,tag,'D_GG_BKG(Sgg vs. Bgg)',
     & D_MELA,wt,wt2,0._dp,1._dp,0.005_dp,'lin')
      n=n+1
      
      call get_MELA_Discr_ggZZ_BSM_1(p,D_MELA)
      call bookplot(n,tag,'D_BSM_1(S_BSM vs. S_SM v1)',
     & D_MELA,wt,wt2,0._dp,1._dp,0.005_dp,'lin')
      n=n+1
      
      call get_MELA_Discr_ggZZ_BSM_2(p,D_MELA)
      call bookplot(n,tag,'D_BSM_2(S_BSM vs. S_SM v2)',
     & D_MELA,wt,wt2,-1._dp,1._dp,0.01_dp,'lin')
      n=n+1
            
      call get_MELA_Discr_ppZZ(p,D_MELA)
      call bookplot(n,tag,'D_PP_BKG(Spp vs. Bpp)',
     & D_MELA,wt,wt2,0._dp,1._dp,0.005_dp,'lin')
      n=n+1
c--- End Modification -- MS & AP
        
        
        
        
      
c--- usual plots for 3+4
      call autoplot2(p,34,3,4,tag,wt,wt2,n)

c--- usual plots for 5+6
      call autoplot2(p,56,5,6,tag,wt,wt2,n)

      if (interference) then
c--- usual plots for 3+6
        call autoplot2(p,36,3,6,tag,wt,wt2,n)

c--- usual plots for 4+5
        call autoplot2(p,45,4,5,tag,wt,wt2,n)
      endif

c--- usual plots for 3+4+5+6
      call autoplot4(p,3456,3,4,5,6,tag,wt,wt2,n)

c--- additional plots that may be present at NLO       
      if (abs(p(7,4)) > 1.e-8_dp) then
        call autoplot1(p,7,tag,wt,wt2,n)
      else
        n=n+2
      endif

************************************************************************
*                                                                      *
*     FINAL BOOKKEEPING                                                *
*                                                                      *
************************************************************************

c--- We have over-counted the number of histograms by 1 at this point
      n=n-1

c--- Ensure the built-in maximum number of histograms is not exceeded    
      call checkmaxhisto(n)

c--- Set the maximum number of plots, on the first call
      if (first) then
        first=.false.
        nplotmax=n
      endif
      
      return
      end

      
c--- NOTE: The following has been moved to User/mela.f      
c--- Modifications for MELA discriminators -- MS & AP  
C      subroutine get_MELA_Discr_ggZZ(p,D_MELA)
C      implicit none
C      include 'types.f'
C      include 'mxpart.f'      
C      real(dp):: p(mxpart,4),D_MELA      
C      real(dp) :: msq_SIGINT(-5:5,-5:5),msq_SIG(-5:5,-5:5),msq_BKG
C      logical:: SM
C      
C      SM=.false.
C      
C        msq_SIG(:,:)=0._dp    
C        msq_SIGINT(:,:)=0._dp      
C        msq_BKG=0._dp   
C        
C               
C        call gg_hZZ_tb(p,msq_SIG,SM)      ! this is the signal ME (gg->H->ZZ)
C        call gg_zz(p,msq_BKG,SM)          ! this is the gg-bkg ME (gg->ZZ)
!         call gg_zz_Hpi(p,msq_SIGINT) ! this includes the interf. term (hence not positive definit)
C
C        D_MELA = msq_SIG(0,0)/(msq_SIG(0,0)+msq_BKG)      
!         D_MELA = msq_SIGINT(0,0)/(msq_SIGINT(0,0)+msq_BKG)
C
C      
C      return
C      end
C      
C      
C      subroutine get_MELA_Discr_ggZZ_BSM_1(p,D_MELA)
C      implicit none
C      include 'types.f'
C      include 'mxpart.f'      
C      real(dp):: p(mxpart,4),D_MELA      
C      real(dp):: msq_SM(-5:5,-5:5),msq_BSM(-5:5,-5:5)
C      
C        msq_SM(:,:)=0._dp    
C        msq_BSM(:,:)=0._dp      
C        
C               
C        call gg_hZZ_tb(p,msq_SM,.true.)  ! this is the BSM signal ME (gg->H->ZZ)
C        call gg_hZZ_tb(p,msq_BSM,.false.)  ! this is the SM signal ME (gg->H->ZZ)
C
C        D_MELA = msq_BSM(0,0)/(msq_BSM(0,0)+msq_SM(0,0))      
C      
C      return
C      end
C      
C      subroutine get_MELA_Discr_ggZZ_BSM_2(p,D_MELA)
C      implicit none
C      include 'types.f'
C      include 'mxpart.f'      
C      real(dp):: p(mxpart,4),D_MELA      
C      real(dp):: msq_SM(-5:5,-5:5),msq_BSM(-5:5,-5:5)
C      
C        msq_SM(:,:)=0._dp    
C        msq_BSM(:,:)=0._dp      
C        
C               
C        call gg_hZZ_tb(p,msq_SM,.true.)  ! this is the BSM signal ME (gg->H->ZZ)
C        call gg_hZZ_tb(p,msq_BSM,.false.)  ! this is the SM signal ME (gg->H->ZZ)
C
C        D_MELA = (msq_BSM(0,0)-msq_SM(0,0))/(msq_BSM(0,0)+msq_SM(0,0))      
C      
C      return
C      end
C      
C
C      subroutine get_MELA_Discr_PPZZ(p,D_MELA)
C      implicit none
C      include 'types.f'
C      include 'energy.f'      
C      include 'mxpart.f'  
C      include 'nf.f'
C      include 'facscale.f'            
C      include 'nflav.f'      
C      real(dp):: p(mxpart,4),D_MELA      
C      real(dp) :: msq_SIG(-5:5,-5:5),msq_BKGg,msq_BKGq(-5:5,-5:5)
C      real(dp) :: fx1(-nf:nf),fx2(-nf:nf),xmsqjk,xx(1:2)
C      real(dp) :: Collider_Energy,Etot,Pztot
C      integer:: ih1,ih2,j,k
C      common/density/ih1,ih2
C      logical:: SM
C      
C      SM=.false.
C      
C        msq_SIG(:,:)=0._dp    
C        msq_BKGg=0._dp
C        msq_BKGq(:,:) = 0d0
C        
C               
C        call gg_hZZ_tb(p,msq_SIG,SM)      ! this is the signal ME (gg->H->ZZ)
C        call gg_zz(p,msq_BKGg,SM)         ! this is the gg-bkg ME (gg->ZZ)
C        call qqb_zz(p,msq_BKGq)        ! this is the qq-bkg ME (qqb->ZZ)
C
C        
C        Collider_Energy = sqrts
C        Etot = p(3,4)+p(4,4)+p(5,4)+p(6,4)
C        Pztot= p(3,3)+p(4,3)+p(5,3)+p(6,3)
C        xx(1) = (Etot+Pztot)/Collider_Energy
C        xx(2) = (Etot-Pztot)/Collider_Energy
C        
C        call fdist(ih1,xx(1),facscale,fx1)
C        call fdist(ih2,xx(2),facscale,fx2)
C
C        msq_SIG(0,0) = fx1(0)*fx2(0) * msq_SIG(0,0)
C        msq_BKGg     = fx1(0)*fx2(0) * msq_BKGg
C
C        xmsqjk = 0d0
C        do j=-nflav,nflav
C        do k=-nflav,nflav                
C          xmsqjk = xmsqjk + fx1(j)*fx2(k) * msq_BKGq(j,k)
C        enddo 
C        enddo
C       
C        
C        D_MELA = msq_SIG(0,0)/(msq_SIG(0,0)+msq_BKGg+xmsqjk)      
C
C      
C      return
C      end
c--- End Modification -- MS & AP     
      
      
