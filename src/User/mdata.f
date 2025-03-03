      subroutine mdata
      implicit none
      include 'types.f'
************************************************************************
*     Calculational scheme for EW couplings                            *
************************************************************************
c
c     ewscheme=-1  : Old MCFM default 
c                    input values = Gf,alpha(m_Z),m_W,m_Z
c                    output values = sin^2(theta_W),mtop
c
c     ewscheme=0   : Old MadEvent default (= AlpGen with iewopt=2)
c                    input values = sin^2(theta_W),alpha(m_Z),m_Z
c                    output values = m_W,Gf.
c
c     ewscheme=1   : New MCFM default, also Madevent default, "G_mu scheme"
c                    = LUSIFER and AlpGen (iewopt=3) defaults
c                    input values = G_F,m_Z,m_W
c                    output values = sin^2(theta_W),alpha(m_Z).
c
c     ewscheme=2   : input  values = G_F,sin^2(theta_W),alpha(m_Z)
c                    output values = m_W,m_Z.
c
c     ewscheme=3   : User choice. All parameters are left as they are
c                    input here. You have to know what you're doing.
c
      
      include 'ewinput.f'
      data ewscheme  / +1                  /   ! Chooses EW scheme
      data Gf_inp    / 1.16639e-5_dp          /   ! G_F
      data aemmz_inp / 7.7585538055706e-03_dp /   ! alpha_EM(m_Z)=1/128.89
      data xw_inp    / 0.2223_dp            /   ! sin^2(theta_W)
      data wmass_inp / 80.385_dp            /   ! W mass
      data zmass_inp / 91.1876_dp           /   ! Z mass
************************************************************************


************************************************************************
*     Definition to use for computing transverse quantities            *
*         useEt=.false.    transverse momentum [previous default]      *
*         useEt=.true.     transverse energy                           *
************************************************************************
      include 'useet.f'
      data useEt/.false./
************************************************************************


************************************************************************
*     Masses, widths and initial-state flavour information             *
************************************************************************
      include 'masses.f'
      include 'nflav.f'
      include 'nores.f'
c--- if true, nores removes all of the gg contribution
      data nores/.false./
c--- Masses: note that "mtausq" is typically used throughout the
c--- program to calculate couplings that depend on the mass, while
c--- "mtau" is the mass that appears in the rest of the matrix
c--- elements and phase space (and may be set to zero in the program,
c--- depending on the process number) 

      data mtau,mtausq/1.777_dp,3.157729_dp/
c----   Note: after v5.6, the masses for top, bottom and charm quarks
c----         are set in the input file

c---  Widths: note that the top width is calculated in the program
c---  The W width of 2.1054 is derived using the measured BR of
c---    10.86 +/- 0.09 % (PDG 2015) and the LO partial width calculation
c---    for Mw=80.385 GeV
      data wwidth,zwidth/2.093_dp,2.4952_dp/
      data tauwidth/2.269e-12_dp/
c--- Number of active flavours in the initial state: this parameter
c--- may be changed in the program for some processes
      data nflav/5/
c--- Masses below here are currently unused      
      data md,mu,ms/66.e-3_dp,66.e-3_dp,0.15_dp/
      data mel,mmu/0.51099892e-3_dp,0.105658369_dp/

************************************************************************


************************************************************************
*     CKM matrix entries                                               *
************************************************************************
      real(dp):: Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/cabib/Vud,Vus,Vub,Vcd,Vcs,Vcb
      data  Vud  ,  Vus  ,  Vub  ,
     &      Vcd  ,  Vcs  ,  Vcb
     &   /0.975_dp,0.222_dp,0.000_dp,
     &    0.222_dp,0.975_dp,0.000_dp/
!      data  Vud  ,  Vus  ,  Vub  ,
!     &      Vcd  ,  Vcs  ,  Vcb
!     &   /1._dp,0._dp,0.000_dp,
!     &    0._dp,1._dp,0.000_dp/

************************************************************************


************************************************************************
*     MS-bar quark massses, mQ_MS(mQ_MS)                               *
*     computed using 1-loop running from pole masses of                *
*     mb=4.75, mt=173.5                                                *
************************************************************************
      include 'msbarmasses.f'
      data mc_msbar/1.275_dp/
      data mb_msbar/4.37_dp/
      data mt_msbar/166._dp/
************************************************************************


************************************************************************
*     Relevant for the H+b process only :                              *
*       susycoup: the deviation of the Higgs coupling from the         *
*                 Standard Model value (S.M. = 1._dp)                    *
************************************************************************
      include 'susycoup.f'
      data susycoup/1._dp/
************************************************************************


************************************************************************
*     Dim. Reg. parameter epsilon, used for checking the proper        *
*      operation of the NLO code in the program                        *
************************************************************************
      include 'epinv.f'
      include 'epinv2.f'
      data epinv/ 1d3/
      data epinv2/1d3/
************************************************************************
      return
      end
