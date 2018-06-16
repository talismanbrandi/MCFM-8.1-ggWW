c--- Modifications for MELA discriminators -- MS & AP  
      subroutine get_MELA_Discr_ggZZ(p,D_MELA)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'      
      real(dp):: p(mxpart,4),D_MELA      
      real(dp) :: msq_SIGINT(-5:5,-5:5),msq_SIG(-5:5,-5:5),msq_BKG
      
        msq_SIG(:,:)=0._dp    
        msq_SIGINT(:,:)=0._dp      
        msq_BKG=0._dp   
        
               
        call gg_hZZ_tb(p,msq_SIG,.true.,one,zip)      ! this is the signal ME (gg->H->ZZ)
        call gg_zz(p,msq_BKG,.true.,one,zip,one,one)    ! this is the gg-bkg ME (gg->ZZ)
!         call gg_zz_Hpi(p,msq_SIGINT) ! this includes the interf. term (hence not positive definit)

        D_MELA = msq_SIG(0,0)/(msq_SIG(0,0)+msq_BKG)      
!         D_MELA = msq_SIGINT(0,0)/(msq_SIGINT(0,0)+msq_BKG)

      
      return
      end
      

      subroutine get_MELA_Discr_PPZZ(p,D_MELA)
      implicit none
      include 'types.f'
      include 'energy.f'      
      include 'mxpart.f'  
      include 'nf.f'
      include 'facscale.f'            
      include 'nflav.f'
      include 'constants.f'      
      real(dp):: p(mxpart,4),D_MELA      
      real(dp) :: msq_SIG(-5:5,-5:5),msq_BKGg,msq_BKGq(-5:5,-5:5)
      real(dp) :: fx1(-nf:nf),fx2(-nf:nf),xmsqjk,xx(1:2)
      real(dp) :: Collider_Energy,Etot,Pztot
      integer:: ih1,ih2,j,k
      common/density/ih1,ih2
      
        msq_SIG(:,:)=0._dp    
        msq_BKGg=0._dp
        msq_BKGq(:,:) = 0._dp
        
               
        call gg_hZZ_tb(p,msq_SIG,.true.,one,zip)      ! this is the signal ME (gg->H->ZZ)
        call gg_zz(p,msq_BKGg,.true.,one,zip,one,one)   ! this is the gg-bkg ME (gg->ZZ)
        call qqb_zz(p,msq_BKGq)                      ! this is the qq-bkg ME (qqb->ZZ)

        
        Collider_Energy = sqrts
        Etot = p(3,4)+p(4,4)+p(5,4)+p(6,4)
        Pztot= p(3,3)+p(4,3)+p(5,3)+p(6,3)
        xx(1) = (Etot+Pztot)/Collider_Energy
        xx(2) = (Etot-Pztot)/Collider_Energy
        
        call fdist(ih1,xx(1),facscale,fx1)
        call fdist(ih2,xx(2),facscale,fx2)

        msq_SIG(0,0) = fx1(0)*fx2(0) * msq_SIG(0,0)
        msq_BKGg     = fx1(0)*fx2(0) * msq_BKGg

        xmsqjk = 0d0
        do j=-nflav,nflav
        do k=-nflav,nflav                
          xmsqjk = xmsqjk + fx1(j)*fx2(k) * msq_BKGq(j,k)
        enddo 
        enddo
       
        
        D_MELA = msq_SIG(0,0)/(msq_SIG(0,0)+msq_BKGg+xmsqjk)      

      
      return
      end

      subroutine get_MELA_Discr_PPZZ_CMS(p,D_MELA)
      implicit none
      include 'types.f'
      include 'energy.f'      
      include 'mxpart.f'  
      include 'nf.f'
      include 'facscale.f'            
      include 'nflav.f'
      include 'constants.f'      
      real(dp) :: p(mxpart,4),D_MELA      
      real(dp) :: msq_gg(-5:5,-5:5)
      real(dp) :: msq_qq(-5:5,-5:5)
      real(dp) :: fx1(-nf:nf),fx2(-nf:nf),xmsqjk,xx(1:2)
      real(dp) :: Collider_Energy,Etot,Pztot
      integer:: ih1,ih2,j,k
      common/density/ih1,ih2
      
        msq_gg(:,:)=0._dp
        msq_qq(:,:)=0._dp
        
        call gg_zz_all(p,msq_gg,.true.,sqrt(10._dp),zip,one,one)   ! this is the signal ME (gg->H->ZZ + gg->ZZ)
        call qqb_zz(p,msq_qq)                                ! this is the qq-bkg ME (qqb->ZZ)

        
        Collider_Energy = sqrts
        Etot = p(3,4)+p(4,4)+p(5,4)+p(6,4)
        Pztot= p(3,3)+p(4,3)+p(5,3)+p(6,3)
        xx(1) = (Etot+Pztot)/Collider_Energy
        xx(2) = (Etot-Pztot)/Collider_Energy
        
        call fdist(ih1,xx(1),facscale,fx1)
        call fdist(ih2,xx(2),facscale,fx2)

        msq_gg(0,0)  = fx1(0)*fx2(0) * msq_gg(0,0)

        xmsqjk = 0d0
        do j=-nflav,nflav
        do k=-nflav,nflav                
          xmsqjk = xmsqjk + fx1(j)*fx2(k) * msq_qq(j,k)
        enddo 
        enddo
        
        D_MELA = msq_gg(0,0)/(msq_gg(0,0)+xmsqjk)      

      return
      end
      
      subroutine get_MELA_Discr_PPZZ_BSM(p,D_MELA,cti,cgi)
      implicit none
      include 'types.f'
      include 'energy.f'      
      include 'mxpart.f'  
      include 'nf.f'
      include 'facscale.f'            
      include 'nflav.f'
      include 'constants.f'      
      real(dp) :: p(mxpart,4),D_MELA      
      real(dp) :: msq_gg(-5:5,-5:5)
      real(dp) :: msq_qq(-5:5,-5:5)
      real(dp) :: fx1(-nf:nf),fx2(-nf:nf),xmsqjk,xx(1:2)
      real(dp) :: Collider_Energy,Etot,Pztot
      real(dp) :: cti,cgi
      integer:: ih1,ih2,j,k
      common/density/ih1,ih2
      
        msq_gg(:,:)=0._dp
        msq_qq(:,:)=0._dp
        
        call gg_zz_all(p,msq_gg,.true.,cti,cgi,one,one)      ! this is the signal ME (gg->H->ZZ + gg->ZZ)
        call qqb_zz(p,msq_qq)                                ! this is the qq-bkg ME (qqb->ZZ)

        
        Collider_Energy = sqrts
        Etot = p(3,4)+p(4,4)+p(5,4)+p(6,4)
        Pztot= p(3,3)+p(4,3)+p(5,3)+p(6,3)
        xx(1) = (Etot+Pztot)/Collider_Energy
        xx(2) = (Etot-Pztot)/Collider_Energy
        
        call fdist(ih1,xx(1),facscale,fx1)
        call fdist(ih2,xx(2),facscale,fx2)

        msq_gg(0,0)  = fx1(0)*fx2(0) * msq_gg(0,0)

        xmsqjk = 0d0
        do j=-nflav,nflav
        do k=-nflav,nflav                
          xmsqjk = xmsqjk + fx1(j)*fx2(k) * msq_qq(j,k)
        enddo 
        enddo
        
        D_MELA = msq_gg(0,0)/(msq_gg(0,0)+xmsqjk)      

      return
      end  
        
      subroutine get_MELA_Discr_ggZZ_BSM(p,D_MELA,cti,cgi)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'      
      real(dp):: p(mxpart,4),D_MELA      
      real(dp) :: msq_BSM(-5:5,-5:5),msq_SM(-5:5,-5:5)
      real(dp) :: cti,cgi
      
        msq_BSM(:,:)=0._dp    
        msq_SM(:,:)=0._dp      
        
               
        call gg_ZZ_all(p,msq_BSM,.true.,cti,cgi,one,one)      ! this is the signal ME (gg->H->ZZ + gg->ZZ) BSM
        call gg_ZZ_all(p,msq_SM,.true.,one,zip,one,one)       ! this is the signal ME (gg->H->ZZ + gg->ZZ) SM

        D_MELA = msq_SM(0,0)/(msq_SM(0,0)+msq_BSM(0,0))

      
      return
      end
c--- End Modification -- MS & AP
