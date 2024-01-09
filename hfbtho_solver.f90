!***********************************************************************
!
!    Copyright (c) 2020, Lawrence Livermore National Security, LLC.
!                        Produced at the Lawrence Livermore National
!                        Laboratory.
!                        Lead developer: Nicolas Schunck, schunck1@llnl.gov
!    HFODD
!    -----
!      LLNL-CODE-710577 All rights reserved.
!      LLNL-CODE-470611 All rights reserved.
!
!      Copyright 2017, N. Schunck, J. Dobaczewski, W. Satula, P. Baczyk,
!                      J. Dudek, Y. Gao, M. Konieczka, K. Sato, Y. Shi,
!                      X.B. Wang, and T.R. Werner
!      Copyright 2012, N. Schunck, J. Dobaczewski, J. McDonnell,
!                      W. Satula, J.A. Sheikh, A. Staszczak,
!                      M. Stoitsov, P. Toivanen
!      Copyright 2009, J. Dobaczewski, W. Satula, B.G. Carlsson, J. Engel,
!                      P. Olbratowski, P. Powalowski, M. Sadziak,
!                      J. Sarich, N. Schunck, A. Staszczak, M. Stoitsov,
!                      M. Zalewski, H. Zdunczuk
!      Copyright 2004, 2005, J. Dobaczewski, P. Olbratowski
!      Copyright 1997, 2000, J. Dobaczewski, J. Dudek
!
!    HFBTHO
!    -----
!      LLNL-CODE-728299 All rights reserved.
!      LLNL-CODE-573953 All rights reserved.
!
!      Copyright 2017, R. Navarro Perez, N. Schunck, R. Lasseri, C. Zhang,
!                      J. Sarich
!      Copyright 2012, M.V. Stoitsov, N. Schunck, M. Kortelainen, H.A. Nam,
!                      N. Michel, J. Sarich, S. Wild
!      Copyright 2005, M.V. Stoitsov, J. Dobaczewski, W. Nazarewicz, P.Ring
!
!
!    This file is part of DFTNESS.
!
!    DFTNESS is free software: you can redistribute it and/or modify it
!    under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    DFTNESS is distributed in the hope that it will be useful, but
!    WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with DFTNESS. If not, see <http://www.gnu.org/licenses/>.
!
!    OUR NOTICE AND TERMS AND CONDITIONS OF THE GNU GENERAL PUBLIC
!    LICENSE
!
!    Our Preamble Notice
!
!      A. This notice is required to be provided under our contract
!         with the U.S. Department of Energy (DOE). This work was
!         produced at the Lawrence Livermore National Laboratory under
!         Contract No. DE-AC52-07NA27344 with the DOE.
!      B. Neither the United States Government nor Lawrence Livermore
!         National Security, LLC nor any of their employees, makes any
!         warranty, express or implied, or assumes any liability or
!         responsibility for the accuracy, completeness, or usefulness
!         of any information, apparatus, product, or process disclosed,
!         or represents that its use would not infringe privately-owned
!         rights.
!      C. Also, reference herein to any specific commercial products,
!         process, or services by trade name, trademark, manufacturer
!         or otherwise does not necessarily constitute or imply its
!         endorsement, recommendation, or favoring by the United States
!         Government or Lawrence Livermore National Security, LLC. The
!         views and opinions of authors expressed herein do not
!         necessarily state or reflect those of the United States
!         Government or Lawrence Livermore National Security, LLC, and
!         shall not be used for advertising or product endorsement
!         purposes.
!
!    The precise terms and conditions for copying, distribution and
!    modification are contained in the file COPYING.
!
!***********************************************************************

! ==================================================================== !
!                                                                      !
!                        DFT SOLVER PACKAGE                            !
!                                                                      !
! ==================================================================== !

!----------------------------------------------------------------------
!> This module provides the main \theCode DFT solver. It includes
!> routines for the calculation and diagonalization of the HFB matrix;
!> the definition of densities on the quadrature mesh and in configuration
!> space; the self-consistent loop the calculation of expectation values of
!> observables. It also includes a package of a few routines to perform
!> particle number projection.
!>
!>  @author
!>    Mario Stoitsov, Nicolas Schunck, Markus kortelainen, Rodrigo Navarro Perez
!----------------------------------------------------------------------
!  Subroutines: - HFBTHO_DFT_SOLVER
!               - heading
!               - thodefh
!               - preparer
!               - coordinateLST
!               - iter
!               - hfbdiag
!               - ALambda
!               - Canonical
!               - resu
!               - initialize_HFBTHO_NAMELIST
!               - read_HFBTHO_NAMELIST
!               - check_consistency
!               - initialize_HFBTHO_SOLVER
!               - base0
!               - base
!               - gaupol
!               - start
!               - printRHO
!               - nucleus
!               - stab
!               - coulom1
!               - coulom
!               - coulom_test
!               - optHFBTHO
!               - DENSIT
!               - field
!               - gamdel
!               - gamdel_gogny
!               - broyden_min
!               - expect
!               - Constraint_or_not
!               - getLagrange
!               - requested_blocked_level
!----------------------------------------------------------------------
Module HFBTHO_solver

  Use HFBTHO_utilities
  Use HFBTHO
  Use HFBTHO_multipole_moments
  Use HFBTHO_collective
  Use HFBTHO_fission_fragments
  Use HFBTHO_gogny
  Use HFBTHO_Lipkin
  Use HFBTHO_PNP
  Use HFBTHO_io
#if(USE_QRPA==1)
  Use HFBTHO_storage
#endif
  Use math

  implicit none

  !===============================QQst
  ! lpfam variables
  !===============================
  !=================================================================================================================================
  ! time-even fields
  Complex(pr), Public, Allocatable         :: vhbn_c(:),vn_c(:),vrn_c(:),vzn_c(:),vdn_c(:),vsn_c(:),davn_c(:),dbvn_c(:)  
  Complex(pr), Public, Allocatable         :: vhbp_c(:),vp_c(:),vrp_c(:),vzp_c(:),vdp_c(:),vsp_c(:),davp_c(:),dbvp_c(:)    
  Complex(pr), Public, Allocatable, Target :: vJRRn_c(:),vJRFIn_c(:),vJRZn_c(:),vJFIRn_c(:),vJFIFIn_c(:)
  Complex(pr), Public, Allocatable, Target :: vJFIZn_c(:),vJZRn_c(:),vJZFIn_c(:),vJZZn_c(:)
  Complex(pr), Public, Allocatable, Target :: vJRRp_c(:),vJRFIp_c(:),vJRZp_c(:),vJFIRp_c(:),vJFIFIp_c(:)
  Complex(pr), Public, Allocatable, Target :: vJFIZp_c(:),vJZRp_c(:),vJZFIp_c(:),vJZZp_c(:)
  Complex(pr), Public, Allocatable         :: cou_c(:), vc_c(:,:) 
  ! time-even densities
  Complex(pr), Public, Allocatable, Target :: aka_c(:,:),bka_c(:,:),ro_c(:,:),tau_c(:,:),dro_c(:,:),dj_c(:,:) 
  Complex(pr), Public, Allocatable, Target :: JRR_c(:,:),JRFI_c(:,:),JRZ_c(:,:),JFIR_c(:,:),JFIFI_c(:,:)
  Complex(pr), Public, Allocatable, Target :: JFIZ_c(:,:),JZR_c(:,:),JZFI_c(:,:),JZZ_c(:,:)
  Complex(pr), Public, Allocatable, Target :: NABLAR_c(:,:),NABLAZ_c(:,:)
  Complex(pr), Public, Allocatable, Target :: rk_c(:,:),ak_c(:,:)   !! rk_c for LN. Not used
  ! HFB densitieds for density dependent part
  Complex(pr), Public, Allocatable, Target :: ro_c_hfb(:,:),tau_c_hfb(:,:),kap_c_hfb(:,:)
  ! time-odd coupling constants
  Real(pr)                                 :: Css(0:1),CDss(0:1),CDeltass(0:1),Cjcjc(0:1),Csnablajc(0:1),CsT(0:1),CsF(0:1)
  ! time-odd fields
  Complex(pr), Public, Allocatable         :: vSIGRn_c(:),vSIGZn_c(:),vSIGFIn_c(:)     !! s, neutrons
  Complex(pr), Public, Allocatable         :: vSIGRp_c(:),vSIGZp_c(:),vSIGFIp_c(:)     !! s, protons
  Complex(pr), Public, Allocatable         :: vSIDRn_c(:),vSIDZn_c(:),vSIDFIn_c(:)     !! \delta s, neutrons
  Complex(pr), Public, Allocatable         :: vSIDRp_c(:),vSIDZp_c(:),vSIDFIp_c(:)     !! \delta s, protons
  Complex(pr), Public, Allocatable         :: vjRn_c(:),vjZn_c(:),vjFIn_c(:)           !! j, neutrons
  Complex(pr), Public, Allocatable         :: vjRp_c(:),vjZp_c(:),vjFIp_c(:)           !! j, protons  
  Complex(pr), Public, Allocatable         :: vTRn_c(:),vTZn_c(:),vTFIn_c(:)           !! T, neutrons
  Complex(pr), Public, Allocatable         :: vTRp_c(:),vTZp_c(:),vTFIp_c(:)           !! T, protons  
  ! time-odd densities.
  Complex(pr), Public, Allocatable, Target :: sor_c(:,:),soz_c(:,:),sofi_c(:,:)        !! spin
  Complex(pr), Public, Allocatable, Target :: dsor_c(:,:),dsoz_c(:,:),dsofi_c(:,:)     !! laplacian spin
  Complex(pr), Public, Allocatable, Target :: rjr_c(:,:),rjz_c(:,:),rjfi_c(:,:)        !! current (= j)
  Complex(pr), Public, Allocatable, Target :: tsor_c(:,:),tsoz_c(:,:),tsofi_c(:,:)     !! spin kinetic (= T)
  Complex(pr), Public, Allocatable, Target :: curjr_c(:,:),curjz_c(:,:),curjfi_c(:,:)  !! curl current (= curl j)
  Complex(pr), Public, Allocatable, Target :: cursr_c(:,:),cursz_c(:,:),cursfi_c(:,:)  !! curl spin
  ! gamdel matrix elements
  Complex(pr), Public, Allocatable         :: hN_c(:),hP_c(:),dN_c(:),dP_c(:)
  Complex(pr), Public, Allocatable         :: BhN_c(:),BhP_c(:),BdN_c(:),BdP_c(:)
  Complex(pr), Public, Allocatable         :: h11N_c(:),h11P_c(:),Bh11N_c(:),Bh11P_c(:)
  Complex(pr), Public, Allocatable         :: FhNab_c(:),FhPab_c(:),BFhNab_c(:),BFhPab_c(:) !! operator F matrix elements 
  Complex(pr), Public, Allocatable         :: FhN_c(:),FhP_c(:),BFhN_c(:),BFhP_c(:)         !! F20 and F02
  Complex(pr), Public, Allocatable         :: F11N_c(:),F11P_c(:),BF11N_c(:),BF11P_c(:)     !! F11 and F11b
  ! operators to remove spurious states
  Complex(pr), Public, Allocatable         :: JyNab_c(:),JyPab_c(:),BJyNab_c(:),BJyPab_c(:) !! J_y operator
  Complex(pr), Public, Allocatable         :: JyN_c(:),JyP_c(:),BJyN_c(:),BJyP_c(:)         !! J_y^20 and J_y^02
  ! uv amplitudes
  Complex(pr), Public, Allocatable         :: GreenNplus(:),GreenPplus(:),GreenNminus(:),GreenPminus(:)  
  Complex(pr), Public, Allocatable         :: GreenNefa(:),GreenPefa(:)
  Complex(pr), Public, Allocatable, Target :: VqpN(:),UqpN(:),VqpP(:),UqpP(:) 
  ! eta density matrices for DENSIT
  Complex(pr), Public, Allocatable :: rhoetadm(:),kapetadm(:),kapetabdm(:)
  ! x,y amplitudes
  Complex(pr), Public, Allocatable, Target :: XN_c(:),YN_c(:),XP_c(:),YP_c(:)
  ! p amplitudes
  Complex(pr), Public, Allocatable, Target :: P1N_c(:),P2N_c(:),P1P_c(:),P2P_c(:)
  !
  ! f-factors for efa
  Complex(pr), Public, Allocatable         :: efa_fN_factor1(:),efa_fN_factor2(:)
  Complex(pr), Public, Allocatable         :: efa_fP_factor1(:),efa_fP_factor2(:)
  !
  ! temporary arrays
  Complex(pr), Public, Allocatable         :: qrpa_AUX(:),qrpa_BUX(:),qrpa_AUX20(:),qrpa_BUX20(:)
  Complex(pr), Public, Allocatable         :: hfb_UA(:),hfb_VA(:)
  !
  ! indexing of blocks
  integer(ipr), Public                     :: iHFB_NB, iQRPA_NB                      !! total number of HFB  and QRPA blocks
  integer(ipr), Public                     :: iHFB_NQP                               !! total number of HFB quasiparticles
  integer(ipr), Public                     :: iHFB_LEN, iQRPA_LEN                    !! total length of HFB and QRPA data arrays     
  integer(ipr), Public, Allocatable        :: iHFB_ID(:), iQRPA_IDI(:), iQRPA_IDF(:) !! dimensions of HFB and QRPA blocks
  integer(ipr), Public, Allocatable        :: iHFB_LP(:), iQRPA_LPI(:), iQRPA_LPF(:) !! Lambda+ of HFB, and Lambda+_i and Lambda+_f of QRPA blocks
  integer(ipr), Public, Allocatable        :: iHFB_OM(:), iQRPA_OMI(:), iQRPA_OMF(:) !! Omega of HFB, and Omega_i and Omega_f of QRPA blocks. Number is 2*real value. E.g. Omega=-3/2 is stored as -3
  integer(ipr), Public, Allocatable        :: iHFB_PR(:)                             !! Parity of the HFB block
  integer(ipr), Public, Allocatable        :: iHFB_DS(:), iQRPA_DS(:)                !! index indicating where the block data starts for given block
  integer(ipr), Public, Allocatable        :: iQRPA_BLI(:), iQRPA_BLF(:)             !! HFB block indeces of the qrpa block referring to the HFBTHO basis wave functions.
  integer(ipr), Public, Allocatable        :: iQRPA_BHI(:), iQRPA_BHF(:)             !! HFB block indeces of the qrpa block referring to the iHFB_??() indexing
  integer(ipr), Public, Allocatable        :: iQRPA_TRANSPOSE(:)                     !! index of QRPA block which corresponds transpose of the given qrpa block index
  Integer(ipr), Public, Allocatable        :: iHFB_QPBS(:)                           !! index of first q.p. state strating a block. Dimension = iHFB_NB
  Integer(ipr), Public, Allocatable        :: ipwi_hfb(:,:)                          !! states inside pairing window
  !
  ! qrpa_Broyden
  Character(1), Public                     :: qrpa_bbroyden
  Real(pr), Public, Allocatable            :: qrpa_broin(:),qrpa_broout(:)
  ! qrpa common variables
  Complex(pr), Public                      :: snz_qrpa(2),dsnz_qrpa(2),qrpa_omega
  Integer(ipr), Public                     :: gfpoles_imN,gfpoles_imP,iter_qrpa
  Integer(ipr)                             :: icacou_qrpa
  Real(pr), Public                         :: ImStrengthN,ImStrengthP,ReStrengthN,ReStrengthP
  Real(pr), Public                         :: ImTransStrN,ImTransStrP,ReTransStrN,ReTransStrP
  Real(pr), Public                         :: eeffN, eeffP, gleffN, gleffP, gseffN, gseffP  !! effective charges
  ! qrpa variables taken from the qrpa.inp file and deduced from those
  Integer(ipr), Public                     :: qrpa_calc_mode,qrpa_file_io                                !! calculation mode and file io
  Integer(ipr), Public                     :: T_qrpa_responce,L_qrpa_responce,K_qrpa_responce,P_qrpa_response,QRPA_Operator_type !! isospin-type, L, K, parity, and opr type
  Integer(ipr), Public                     :: max_iter_qrpa,max_qrpa_points
  Integer(ipr), Public                     :: qrpa_nbroyden
  Real(pr), Public                         :: qrpa_alphamix
  Real(pr), Public                         :: qrpa_eps
  Complex(pr), Public                      :: qrpa_eta !! unused
  Real(pr), Public                         :: qrpa_Romega,qrpa_Iomega,qrpa_Romega_step,qrpa_Iomega_step  !! line params
  Real(pr), Public                         :: qrpa_circRe,qrpa_circIm,qrpa_circRad                       !! circle params
  Real(pr), Public                         :: qrpa_HcircRe,qrpa_HcircRad                                 !! half-circle params
  !
  ! various setups
  Logical                                  :: USE_HFB_PWI                            !! If true, densities are calculated with HFB pairing window
  Logical                                  :: USE_ONLY_F20 = .false.                 !! If true, X and Y directly from F20,F02,F11. H20, H02, H11 ignored. _use for debug only_
  Logical                                  :: USE_ONLY_F11EFA = .false.               !! If true, only F11 with efa used. _use for debug only_
  Logical                                  :: USE_NORMALHO_K0 = .false.              !! If true, normal HO basis used for K=0 case
  Logical                                  :: ENABLE_N_EFA = .false.                 !! If true, activates neutron EFA. Set true in initialization, if needed
  Logical                                  :: ENABLE_P_EFA = .false.                 !! If true, activates proton EFA. Set true in initialization, if needed
  !
  !=================================================================================================================================
  !=================================================================================================================================
  ! variable for constructing AB matrix
  !-----------------------------------
  ! Four indices for A and B supermatrix
  !-----------------------------------
  integer :: qqalpha, qqbeta, qqgamma, qqsigma
  !-----------------------------------
  ! indices for dH20, X and Y
  integer :: qqmu, qqnu, qqmup, qqnup
  !----------------------------------
  integer :: qqi, qqj
  integer :: WX_sum_nd

  integer :: numberofXN_c
  integer, allocatable :: WX_sum_nd_colm1(:)          !!!Q: dimension is iHFB_NB, record sum of nd of first column minus 1 blocks.
  integer, allocatable, public :: QQ_select(:,:)  !!!Q: qrpaINI variables, made public here os it can be used elsewhere.
  integer, allocatable, public :: AB_row(:), AB_col(:)
  complex, allocatable, public :: AN_qrpa(:,:), BN_qrpa(:,:)  !!: neutron AB matrix
  complex, allocatable, public :: AP_qrpa(:,:), BP_qrpa(:,:)  !!: proton AB matrix
  real:: omega_N0, omega_NN0, omega_PN0, v_N
    !-------------------
    !AB matrix.
    !--------------------
    integer:: ct1, ct2, ct3, inbrr, inbrc, inbcc, qqk





    !===============================QQend






Contains

  !-----------------------------------------------------------------------------------------------
  !
  !   Axially-deformed configurational constrained and/or unconstrained Hartree-Fock-Bogoliubov
  !   calculations with Skyrme-like functionals and delta pairing using the Harmonic-Oscillator
  !   (HO), and/or Transformed HO (THO) basis with or without reflection symmetry imposed, with
  !   or without the Lipkin-Nogami procedure. The solver can handle all Skyrme-like functionals,
  !   DME-functionals, Fayans-functionals, calculate infinite nuclear matter properties, finite
  !   nuclei (even-even, odd-even, odd-odd), and neutron drops, isoscalar and isovector monopo-
  !   le FAM QRPA calculations for spherical and deformed nuclei.
  !
  !   All necessary input variables contain the suffix _INI. Below, the complete list of these
  !   variables with some example values:
  !
  !   ======== hfbtho_NAMELIST.dat
  !   n00_INI=20; npr1_INI=70; npr2_INI=50;  kindhfb_INI=-1; inin_INI=-1
  !   b0_INI=2.234776; q_INI=0.0; cdef_INI=0.0; cqad_INI=0.5; skyrme_INI='SLY4'; nkblo_INI=0
  !   ILST_INI=0; keypj_INI=1; iproj_INI=0; npr1pj_INI=0;
  !   icou_INI=2; IDEBUG_INI=0; npr2pj_INI=0;
  !   Parity_INI=.False.; epsi_INI=0.00001_pr; MAX_ITER_INI=101
  !   Add_Pairing_INI=.False.; DO_FITT_INI=.False.; Print_PTHO_Namelist_INI=.True.
  !
  !   ======== from read_UNEDF_NAMELIST
  !   DMEORDER=-1; DMELDA=0; use_TMR_pairing=0
  !   HBZERO=20.73553000000000;    E2CHARG=1.4399784085965135; CRHO(0)=-933.3423749999999;  CRHO(1)=830.0524855000001;
  !   CDRHO(0)=861.0625000000000;  CDRHO(1)=-1064.2732500000;  CTAU(0)=57.12868750000000;   CTAU(1)=24.65673650000000;
  !   CRDR(0)=-76.99620312499999;  CRDR(1)=15.65713512500000;  CRDJ(0)=-92.25000000000000;  CRDJ(1)=-30.7500000000000;
  !   CJ(0)=17.20961150000000;     CJ(1)=64.57581250000000;    CPV0(0)=-258.2000000000000;  CPV0(1)=-258.2000000000000;
  !   CPV1(0)=0.5000000000000000;  CPV1(1)=0.500000000000000;  SIGMA=0.1666666666666667;    CEXPAR=1.000000000000000;
  !   E_NM=-15.97214914144462;     K_NM=229.9009644826037;     SMASS_NM =1.439546988976078; RHO_NM =0.1595387567117334;
  !   ASS_NM =32.00430281505202;   LASS_NM=45.96175148046161;  VMASS_NM =1.249838547196253;
  !   MPI=0.6995945261023822;      GA=1.290000000000000;       FPI=0.4683223517486062;      C1=-0.1598130000000000;
  !   C3 =-0.6708200000000;        C4 =0.6708200000000000;     CD =-2.062000000000000;      CE=-0.6250000000000;
  !   LAMBDAX =3.547896604156107;  USE_INM=.false.;            USE_CM_COR =.true.;          USE_DME3N_TERMS=.true.;
  !   USE_J2TERMS =.true.;         USE_CHARGE_DENSITY=.false.; PRINT_NAMELIST=.true.;
  !
  ! Memo:
  !  -  inin_INI switches scratch unconstrained (inin=1,2,3) or constrained (inin=100,200,300)
  !     calculations. Unconstrained mode begins with a small number of constrained iterations.
  !  -  inin_INI switches unconstrained (inin=-1,-2,-3) or constrained (inin=-100,-200,-300)
  !     calculations from a previous solution if the latter exists. If not, the solver sets
  !     inin=Abs(inin) and resumes from scratch.
  !  -  The same holds for odd nuclei. If even-even solution for the odd nucleus does not exists
  !     it is calculated first.
  !  -  Print_Screen=T/F for n00_INI=+/-: output is generated and written in thoout.dat file
  !	only if n00_INI>0. For n00_INI<0, the number of shells is set to abs(n00_INI) but all
  !     output is supressed.
  !
  !  -  At the end of the solution, the solver provides all results in the arrays
  !
  !                  nucname,ereslbl(1:2),eres(1:ierest)
  !
  !     which contain:
  !
  !        LBL,BLKN,BLKZ,Jsi,JININ,A,N',Z,Efn,Efp,JEtot,Jbett,Jbetn,Jbetp,JQt,JQn,JQp
  !        JpEn,JpEp,JpDn,JpDp,JAsn,JAsp,Jrt,Jrn,Jrp,Jrc,Jht,Jhn,Jhp,Jqht,Jqhn,Jqhp,
  !        JKINt,JKINn,JKINp,JSO,JCDIR,JCEX,JDisn,JDisp,JV2Mn,JV2Mp,JILST,JKIND,JL,
  !        JECMPAV1,JECMPAV2,JECMPAV3,JA,JN,JZ,ITER,UEtot,Ubett,Ubetn,Ubetp,UQt,UQn,
  !        UQp,Uln,Ulp,UpEn,UpEp,UpDn,UpDp,UAsn,UAsp,Urt,Urn,Urp,Urc,Uht,Uhn,Uhp,
  !        Uqht,Uqhn,Uqhp,UKINT,UKINN,UKINP,USO,UCDIR,UCEX,UDisn,UDisp,UV2Mn,UV2Mp,
  !        UECMT,UECMN,UECMP,UROTT,UROTN,UROTP,USQUJT,USQUJN,USQUJP,UCRANT,UCRANN,
  !        UCRANP,UERIGT,UERIGN,UERIGP,EHFBLN,EHFB,LNbet,LNben,LNbep,LNQt,LNQn,LNQp,
  !        LNpEn,LNpEp,LNpDn,LNpDp,LNrt,LNrn,LNrp,LNrC,LNam2n,LNam2p,LNe2n,LNe2p,
  !        BlEqpN,BlDEqpN,BlOvrN,BlEqpZ,BlDEqpZ,BlOvrZ
  !
  !  -  Standard inputs from               'hfbtho_NAMELIST.dat'
  !  -  USer-defined functional read from  'UNEDF_NAMELIST.DAT'
  !  -  Final solution stored to files     '*.hel' and/or '*.tel'
  !  -  Output is written to files         'thoout.dat', 'thodef.dat' and 'hodef.dat'
  !  -  Output files *.dat may exist as    'thoout','thores','hodenp',
  !                                      'thodenp','thodef','thoene',
  !                                      'thoprc','dat0.1.2.3.4'
  !  -  External accuracy pr/ipr always set in UNEDF module
  !
  !  -  n00    Number of oscillator shells
  !             n00>0 prints to thoout.dat & screen
  !             n00<0 no print at all
  !             n00=0 program stops (NB!)
  !  -  b0     Oscillator Basis parameter b0>0 (If b0<0 it takes a default value)
  !  -  beta0  Value of Basis deformation parameter
  !  -  AN     Number of neutrons N
  !  -  AZ     Number of protons Z
  !  -  FTST'  Fayance forces label
  !  -  kind   Kind of calculations 1: noLN, -1:LN
  !  -  inin   Unconstraint Iterations from peviouse solution
  !             -1: (from spherical *.hel or *.tel file)
  !             -2: (from prolate   *.hel or *.tel file)
  !             -3: (from oblate    *.hel or *.tel file)
  !            Unconstraint Iterations fom scratch with
  !            preliminary constraint at deformation Cbeta, (i):
  !              1: Spherical scratch
  !              2: Prolate scratch
  !              3: Oblate scratch
  !            Constrained calculation (icstr) at Cbeta, see (i)
  !              100, 200, 300 fom scratch
  !             -100,-200,-300 fom previouse solution
  !  -  blNeutrons: a group responsible for blocking a particular neutron level
  !     The group consists of 5 numbers, e.g., for 7-[ 3, 0, 3]: 7 -1 3 0 3
  !      k1  2 \times \Omega
  !         =0: the whole group (k) is disregarded (n0 blocking)
  !         >0: blocking in N+1 nucleus
  !         <0: blocking in N-1 nucleus
  !      k2  parity (+1 or -1); NB! when k2=0, the ground state walker is applied
  !      k3,k4,k5  Nilson quantum numbers
  !  -  blProtons: exactly the same as (k) but for protons
  !
  !-----------------------------------------------------------------------------------------------
  Subroutine HFBTHO_DFT_SOLVER()
    Use HFBTHO_THO
    Implicit None
    Integer(ipr) :: iw,ib,j,i,it,l,maxi0,icstr0,iterMax,icons,il,kickoff,iexit
    Real(pr)     :: epsi0,qq,f,f1,f2,f3,r,g,g1
    Real(pr)     :: ehfb_tmp = 999.999


    !---------------------------QQst
   !Logical      :: lpr
    Integer(ipr) :: iblocase(2),nkblocase(2,5)
    Integer(ipr) :: iodd,ios,envstat,envlen
    integer(ipr) :: qrpa_points_thiscore,qrpa_start_thiscore,idivm
    integer(ipr) :: arrayTaskID,arrayTaskCount,arrayTaskMin,arrayTaskMax
    CHARACTER    :: envRead*20
    real(pr),allocatable :: qrpa_omega_points_re(:),qrpa_omega_points_im(:)
    !----------------------------QQend




    !-------------------------------------------------------------
    ! Initializing all according to *_INI values
    !-------------------------------------------------------------
    Call initialize_HFBTHO_SOLVER
    If(ierror_flag.Ne.0) Return
    !-------------------------------------------------------------
    ! Opening files and deciding when to print based on MPI structure
    ! (team_rank = 0 by default even if no MPI is requested)
    !-------------------------------------------------------------
#if(DO_PES==1)
                       do_print = 0
    If(team_rank == 0) do_print = 1
    If(lout <= lfile .And. do_print == 1) Open(lfile,file='thoout'//row_string//'.dat',status='unknown')
#endif
    ! Mass tables or driplines modes
#if(DO_MASSTABLE==1 || DRIP_LINES==1)
    do_print = 1
    If(lout <= lfile) Open(lfile,file='thoout'//row_string//'.dat',status='unknown')
#endif
    ! All other cases
#if(DO_MASSTABLE==0 && DRIP_LINES==0 && DO_PES==0)
    do_print = 1
    !Open(lout,file='hfbtho.out',status='unknown')
    If(lout < lfile) Open(lfile,file='thoout.dat',status='unknown')
#endif
    Call Constraint_or_not(inin_INI,inin,icstr)
    If(ierror_flag.Ne.0) Return
    !-------------------------------------------------------------------------
    ! Main blocking loop while irestart .Ne. 0 - Tries to read from wel file
    ! if none, scratch calc of even core, then read from wel for blocking
    !-------------------------------------------------------------------------
    irestart=0; iexit=0
    matrix_elements_calculated = .false.
    Do
       n00=Abs(n00_INI);  b0=b0_INI;           q=q_INI;           iLST=iLST_INI;
       maxi=MAX_ITER_INI; npr(1)=npr_INI(1);   npr(2)=npr_INI(2); npr(3)=npr(1)+npr(2);
       skyrme=skyrme_INI; kindhfb=kindhfb_INI
       keypjn=keypjn_INI; keypjp=keypjp_INI;   iproj=iproj_INI;   npr1pj=npr1pj_INI; npr2pj=npr2pj_INI;
       nkblo=nkblo_INI;
       basis_HFODD=basis_HFODD_INI
       mini=1; iiter=0
       !-------------------------------------------------------------
       ! Define the set of constraints
       !-------------------------------------------------------------
       numberCons=0; kickoff=0
       Do l=1,lambdaMax
          If(lambda_active(l).Gt.0) numberCons = numberCons + 1
          If(lambda_active(l).Lt.0) kickoff = kickoff + 1
       End Do
       ! Add constraint on the neck
       If(neck_constraints) Then
          numberCons = numberCons + 1
          neckLag = zero
          !kickoff = kickoff + 1
       End If
       !
       multLambda=0; multLag = zero; multRequested=zero
       icons=0
       Do l=1,lambdaMax
          If(lambda_active(l).Gt.0) Then
             icons=icons+1
             multLambda(icons)=lambda_values(l)
          End If
          multRequested(l) = expectation_values(l)
       End Do
       If(neck_constraints) Then
          icons=icons+1
          multLambda(icons)=0
       End If
       !-------------------------------------------------------------
       ! PARTICLE NUMBER
       !-------------------------------------------------------------
       Do it=1,2
          If(nkblo(it,1).Ne.0) Then
             If(npr_INI(it).Eq.2*(npr_INI(it)/2)) Then
                If(nkblo(it,1).Gt.0) Then
                   ! particle state
                   npr(it)=npr(it)+1
                   iparenti(it)=-1
                Else
                   ! hole state
                   npr(it)=npr(it)-1
                   iparenti(it)=+1
                End If
             Else
                 iparenti(it)=999 ! Odd noblock
             End If
             nkblo(it,1)=Abs(nkblo(it,1))
          End If
       End Do
       npr(3) = npr(1) + npr(2)
       !-------------------------------------------------------------
       ! HFB+HO CALCULATIONS
       !-------------------------------------------------------------
       If(ILST.Le.0) Then
          icacou=0; icahartree=0
          Call preparer(.True.)
          If(ierror_flag.Ne.0) Return
          ! Reading input file: if not possible, start from scratch
          Call inout(1,iexit)
          ! When skyrme='READ', we overwrite everything we can with what is contained in the binary file,
          ! among others, the oscillator lengths b0, bp and bz. Since these are used in the routine
          ! preparer() to set up the characteristics of the HO basis, we need to repeat a number of
          ! operations here to account for the fact that (b0,bp,bz) as read from the binary file may
          ! be different from (b0,bp,bz) as specified by the user. This case should only occur when the
          ! binary file was produced by HFODD. This is because HFODD uses a different prescription to
          ! calculate said bp and bz, which is not equivalent to the HFBTHO prescription. The HFODD
          ! bp and bz are passed as arguments to the HFBTHO kernel.
          If(iexit>0) Then
             Call start
          Else
             If(Trim(functional)=='READ') Then
                Call base0(.True.)  ! basis space (calculate configurational space)
                Call realloc        ! global allocation
                Call gausspoints    ! Gauss mesh points
                Call base(.True.)   ! oscillator configurations (set up quantum numbers)
                Call gaupol(.True.) ! basis wf at gauss mesh points
                Call print_functional_parameters(skyrme) ! print new coupling constants (read from file)
             End If
             Do it=1,2
                tz(it) = Real(npr(it),Kind=pr)
             End Do
             If (VERSION_READ == 5 .And. no_HFB_mat) Call gamdel(.false.,.false.)
          End If
          If(ierror_flag.Ne.0) Return
          !-------------------------------------------------------------
          ! BLOCKING
          !-------------------------------------------------------------
          Call handle_blocking(irestart,iexit)
          If(irestart.Eq.-1) Exit
          If(ierror_flag.Ne.0) Return
          !-------------------------------------------------------------
          ! Preliminary constrained calculations
          !-------------------------------------------------------------
          If(kickoff.Gt.0) Then
             icstr0=icstr; epsi0=epsi; ! remember accuracy
             icstr=1                   ! constraint true
             epsi=1.0_pr               ! small accuracy
             mini = 1                  ! initial iteration number
             if(is_NEDF) then
                iterMax = maxi; maxi = 25
             else
                iterMax = maxi; maxi = 10
             endif
             numberCons=0
             Do l=1,lambdaMax
                If(Abs(lambda_active(l)).Gt.0) Then
                   numberCons=numberCons+1
                   multLambda(numberCons)=lambda_values(l)
                End If
             End Do
             If(do_print == 1) Then
                Do iw=lout,lfile
                   If(Parity) Then
                      Write(iw,'(/,a,i3,a,i2,a,/)') '  ### INITIAL STAGE(constrained calculations, reflection symmetry used)'
                   Else
                      Write(iw,'(/,a,i3,a,i2,a,/)') '  ### INITIAL STAGE(constrained calculations, no reflection symmetry used)'
                   End If
                End Do
             End If
             Call iter(.True., .True.) ! small constraint iterations, keep history
             If(ierror_flag.Ne.0) Return
             ! For the next phase, use only true constraints
             icstr=icstr0; epsi=epsi0
             maxi = iterMax
             numberCons=0
             Do l=1,lambdaMax
                If(lambda_active(l).Gt.0) Then
                   numberCons=numberCons+1
                   multLambda(numberCons)=lambda_values(l)
                End If
             End Do
          End If
          !-------------------------------------------------------------
          ! REGULAR HFB+HO ITERATIONS
          !-------------------------------------------------------------
          If(do_print == 1) Then
             Do iw=lout,lfile
                If(Parity) Then
                   Write(iw,'(/,a,i3,a,i2,a,/)')    '  ### REGULAR STAGE (reflection symmetry imposed)'
                Else
                   Write(iw,'(/,a,i3,a,i2,a,/)')    '  ### REGULAR STAGE (no reflection symmetry imposed)'
                End If
             End Do
          End If
          Call iter(.True., .False.) ! Remove history
          If(ierror_flag.Ne.0) Return
          Call resu(1)
          If(ierror_flag.Ne.0) Return
       End If
       !-------------------------------------------------------------
       ! HFB+THO CALCULATIONS FROM HFB+HO
       !-------------------------------------------------------------
       If(ILST.Lt.0) Then
          ILST1=1; icacou=0; icahartree=0
          Call coordinateLST(.False.) ! THO basis
          If(ierror_flag.Ne.0) Return
          Call densit                 ! THO densities
          If(ierror_flag.Ne.0) Return
          Call field                  ! Nuclear fields
          If(ierror_flag.Ne.0) Return
          Call iter(.True., .False.)  ! HFB+THO iterations, remove history
          If(ierror_flag.Ne.0) Return
          Call resu(1)                ! print/record results
          If(ierror_flag.Ne.0) Return
       End If
       !-------------------------------------------------------------
       ! HFB+HO ITERATIONS
       !-------------------------------------------------------------
       If(ILST.Gt.0) Then
          If(inin.Gt.0) Then
             ierror_flag=ierror_flag+1
             ierror_info(ierror_flag)= ' Forbidden iLST>0, inin>0 '
             Return
          End If
          icacou=0; icahartree=0
          Call preparer(.True.)
          If(ierror_flag.Ne.0) Return
          ! Reading input file: if not possible, start from scratch
          Call inout(1,iexit)
          ! See comment above
          If(iexit>0) Then
             Call start
          Else
             If(Trim(functional)=='READ') Then
                Call base0(.True.)  ! basis space (calculate configurational space)
                Call gausspoints    ! Gauss mesh points
                Call base(.True.)   ! oscillator configurations (set up quantum numbers)
                Call gaupol(.True.) ! basis wf at gauss mesh points
                Call print_functional_parameters(skyrme) ! print new coupling constants (read from file)
             End If
             Do it=1,2
                tz(it) = Real(npr(it),Kind=pr)
             End Do
             If (VERSION_READ == 5 .And. no_HFB_mat) Call gamdel(.false.,.false.)
          End If
          If(ierror_flag.Ne.0) Return
          !-------------------------------------------------------------
          ! BLOCKING
          !-------------------------------------------------------------
          Call handle_blocking(irestart,iexit)
          If(irestart.Eq.-1) Exit
          If(ierror_flag.Ne.0) Return
          Call iter(.True., .False.)  ! HFB+THO iterations, remove history
          If(ierror_flag.Ne.0) Return
          Call resu(1)                ! print/record results
          If(ierror_flag.Ne.0) Return
       End If
       ! ---------------------------------------------------------------
       ! Save only the most bound blocking solution to pnFAM file
       ! (not enabled for large scale modes)
       ! ---------------------------------------------------------------
#if(DO_MASSTABLE==0 && DRIP_LINES==0 && DO_PES==0)
#if(USE_QRPA==1)
       ! Avoid writing solution for even-even core if blocking is on
       If(nkblo(1,1).Eq.abs(nkblo_INI(1,1)) .And. nkblo(2,1).Eq.abs(nkblo_INI(2,1))) Then
          If(ehfb.Lt.ehfb_tmp) Then
             ehfb_tmp = ehfb
             call save_HFBTHO_solution
          End If
       End If
#endif
#endif
       !-------------------------------------------------------------
       ! Go for the requested blocking state in a case of odd nuclei
       ! if restarted due to corrupted/missing previous solution
       !-------------------------------------------------------------
       inin=-Abs(inin)
       If(irestart.Eq.0) Exit
    End Do

    !=======================================QQst
    ! lpfam after hfbtho converge.
    !=======================================
      ! READ qrpa.inp FILE
      !Open(112,file=trim(adjustl(qq_io_path))//'qrpa.inp')
       Open(112,file='qrpa.inp')
       Read(112,*) 
       Read(112,*) qrpa_calc_mode,qrpa_file_io
       Read(112,*) QRPA_Operator_type, T_qrpa_responce,L_qrpa_responce,K_qrpa_responce       
       Read(112,*) qrpa_eps;         
       Read(112,*) max_iter_qrpa;    
       Read(112,*) qrpa_nbroyden
       Read(112,*) qrpa_alphamix;    
       Read(112,*) qrpa_Romega,qrpa_Iomega,qrpa_Romega_step,qrpa_Iomega_step  !! line params
       Read(112,*) qrpa_circRe,qrpa_circIm,qrpa_circRad;                      !! circle params
       Read(112,*) qrpa_HcircRe,qrpa_HcircRad;                                !! half-circle params
       Read(112,*) max_qrpa_points
       Close(112)   
       !
       ! allocate table for all points to be calculated
       if(allocated(qrpa_omega_points_re)) deallocate(qrpa_omega_points_re,qrpa_omega_points_im)
       allocate(qrpa_omega_points_re(max_qrpa_points),qrpa_omega_points_im(max_qrpa_points))
       !
       ! default settings without mpi
       qrpa_start_thiscore = 1                   !! start by default from first item on the table
       qrpa_points_thiscore = max_qrpa_points    !! by default this core calculates all points
       !
       ! set all omega points in tables. all mpi threads will have the same table
       select case(qrpa_calc_mode)
       case(1) !! line
         Do i=1,max_qrpa_points 
           qrpa_omega_points_re(i) = qrpa_Romega + real(i-1,kind=pr)*qrpa_Romega_step
           qrpa_omega_points_im(i) = qrpa_Iomega + real(i-1,kind=pr)*qrpa_Iomega_step
         end do

       case(2)  !! circle
         Do i=1,max_qrpa_points 
           qrpa_omega_points_re(i) = qrpa_circRe + qrpa_circRad * cos( 2.0_pr*pi*real(i-1,kind=pr)/real(max_qrpa_points,kind=pr) )
           qrpa_omega_points_im(i) = qrpa_circIm + qrpa_circRad * sin( 2.0_pr*pi*real(i-1,kind=pr)/real(max_qrpa_points,kind=pr) )
         end do

       case(3)  !! half-circle
         Do i=1,max_qrpa_points 
           qrpa_omega_points_re(i) = qrpa_HcircRe + qrpa_HcircRad * sin( pi*real(i-1,kind=pr)/real(max_qrpa_points-1,kind=pr) )
           qrpa_omega_points_im(i) = -qrpa_HcircRad * cos( pi*real(i-1,kind=pr)/real(max_qrpa_points-1,kind=pr) )
         end do
       case default
         STOP "unknown qrpa calculation mode"
       end select
       !
       ! check file io mode. 0: noting is read or written, 1: write amplitudes, 2: read operators
       if((qrpa_file_io.gt.2) .or. (qrpa_file_io.lt.0)) stop "wrong qrpa_file_io"
       !
       Do i=qrpa_start_thiscore,qrpa_start_thiscore+qrpa_points_thiscore-1
        qrpa_omega=Cmplx(qrpa_omega_points_re(i),qrpa_omega_points_im(i),kind=pr) ! set omega as complex number 
        Call qrpa_HFBTHO_solver  !!!: Use hfbtho_version                          ! solve qrpa for this omega
       End Do
       deallocate(qrpa_omega_points_re,qrpa_omega_points_im)
 

 
    call lptest 
    !================================QQend


    !=======================================

    !
  End Subroutine HFBTHO_DFT_SOLVER
  !=======================================================================
  !> Print heading to screen 'lout' and to tape thoout.dat 'lfile'
  !=======================================================================
  Subroutine heading
#if(USE_OPENMP==1)
    Use omp_lib
#endif
#if(USE_MPI==1)
    Use mpi
#endif
    Implicit None
#if(USE_MPI==1 || USE_MPI==2)
    Integer(ipr) :: rank,ierr
#endif
    Integer(ipr) :: iw,idt(8),numThreads,idThread
    Character(len=12) rcl(3)
    Character(len=50) today
    !
#if(USE_OPENMP==1)
!$OMP PARALLEL SHARED(iw,lout,lfile,do_print) PRIVATE(numThreads,idThread)
    numThreads = omp_get_num_threads()
    idThread = omp_get_thread_num()
#if(DO_MASSTABLE==1)
    If (idThread .Eq. 0.and.irow.le.1) Then
#elif(DRIP_LINES==1)
    If (idThread .Eq. 0.and.irow.le.1.and.calc_counter.eq.1) Then
#elif(USE_MPI==1)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    If(idThread .Eq. 0 .and. rank .eq. 0) Then
#else
    If(idThread .Eq. 0) Then
#endif
       If(do_print == 1) Then
          Do iw=lout,lfile
             Write(iw,'("Multi-threading framework with OpenMP:",i2," threads/task")') numThreads
          End Do
       End If
    End If
!$OMP END PARALLEL
#endif
    Call Date_and_time(rcl(1),rcl(2),rcl(3),idt)
    Write(today,'(a,i2,a,i2,a,i4,a,i2,a,i2,a)')'(',idt(2),'/',idt(3),'/',idt(1),', ',idt(5),':',idt(6),')'
    If(do_print == 1) Then
       Do iw=lout,lfile
          Write(iw,*)
          Write(iw,'("  ==========================================")')
          Write(iw,'("           FORTRAN 95 CODE (KIND=",i2,") ")') pr
          Write(iw,'("               Version: ",a)') Version
          Write(iw,'("  ==========================================")')
          Write(iw,'("       AXIALLY DEFORMED CONFIGURATIONAL     ")')
          Write(iw,'("     HARTREE-FOCK-BOGOLIUBOV CALCULATIONS   ")')
          Write(iw,'("                     WITH                   ")')
          Write(iw,'("      SKYRME+DELTA PAIRING OR GOGNY EDFs    ")')
          Write(iw,'("                  USING THE                 ")')
          Write(iw,'("             HARMONIC-OSCILLATOR            ")')
          Write(iw,'("                    AND/OR                  ")')
          Write(iw,'("        TRANSFORMED HARMONIC-OSCILLATOR     ")')
          Write(iw,'("                    BASIS                   ")')
          Write(iw,'("                     ---                    ")')
          Write(iw,'("                v1.66  (2005)               ")')
          Write(iw,'("    Stoitsov, Dobaczewski, Nazarewicz, Ring ")')
          Write(iw,'("                v2.00d (2012)               ")')
          Write(iw,'("        Stoitsov, Schunck, Kortelainen      ")')
          Write(iw,'("                v3.00  (2016)               ")')
          Write(iw,'("            Schunck, Navarro Perez          ")')
          Write(iw,'("  ==========================================")')
          Write(iw,'("    Nucleus: ",a," (A=",i4,", N=",i3,", Z=",i3,")")') nucname,npr(1)+npr(2),npr(1),npr(2)
          If(Parity) Then
             Write(iw,'("       Reflection Symmetry Imposed       ")')
          Else
             Write(iw,'("      No Reflection Symmetry Imposed     ")')
          End If
          Write(iw,'("            ",a)') today
          Write(iw,'("  ==========================================")')
          Write(iw,*)
          Write(iw,*)
       End Do
    End If
  End Subroutine heading
  !=======================================================================
  !> Print labels to hodef.dat or/and thodef.dat files
  !=======================================================================
  Subroutine thodefh(iw1)
    Implicit None
    Integer(ipr) :: iw1
    hlabels(1)='LBL';       hlabels(11)='JEtot';    hlabels(21)='JpDp';
    hlabels(2)='BLKN';      hlabels(12)='Jbett';    hlabels(22)='JAsn';
    hlabels(3)='BLKZ';      hlabels(13)='Jbetn';    hlabels(23)='JAsp';
    hlabels(4)='Jsi';       hlabels(14)='Jbetp';    hlabels(24)='Jrt';
    hlabels(5)='JININ';     hlabels(15)='JQt';      hlabels(25)='Jrn';
    hlabels(6)='A';         hlabels(16)='JQn';      hlabels(26)='Jrp';
    hlabels(7)='N';         hlabels(17)='JQp';      hlabels(27)='Jrc';
    hlabels(8)='Z';         hlabels(18)='JpEn';     hlabels(28)='Jht';
    hlabels(9)='Efn';       hlabels(19)='JpEp';     hlabels(29)='Jhn';
    hlabels(10)='Efp';      hlabels(20)='JpDn';     hlabels(30)='Jhp';
    !
    hlabels(31)='Jqht';    hlabels(41)='JDisp';    hlabels(51)='JN';
    hlabels(32)='Jqhn';    hlabels(42)='JV2Mn';    hlabels(52)='JZ';
    hlabels(33)='Jqhp';    hlabels(43)='JV2Mp';    hlabels(53)='ITER';
    hlabels(34)='JKINt';   hlabels(44)='JILST';    hlabels(54)='UEtot';
    hlabels(35)='JKINn';   hlabels(45)='JKIND';    hlabels(55)='Ubett';
    hlabels(36)='JKINp';   hlabels(46)='JL';       hlabels(56)='Ubetn';
    hlabels(37)='JSO';     hlabels(47)='JECMPAV1'; hlabels(57)='Ubetp';
    hlabels(38)='JCDIR';   hlabels(48)='JECMPAV2'; hlabels(58)='UQt';
    hlabels(39)='JCEX';    hlabels(49)='JECMPAV3'; hlabels(59)='UQn';
    hlabels(40)='JDisn';   hlabels(50)='JA';       hlabels(60)='UQp';
    !
    hlabels(61)='Uln';     hlabels(71)='Urp';      hlabels(81)='UKINP';
    hlabels(62)='Ulp';     hlabels(72)='Urc';      hlabels(82)='USO';
    hlabels(63)='UpEn';    hlabels(73)='Uht';      hlabels(83)='UCDIR';
    hlabels(64)='UpEp';    hlabels(74)='Uhn';      hlabels(84)='UCEX';
    hlabels(65)='UpDn';    hlabels(75)='Uhp';      hlabels(85)='UDisn';
    hlabels(66)='UpDp';    hlabels(76)='Uqht';     hlabels(86)='UDisp';
    hlabels(67)='UAsn';    hlabels(77)='Uqhn';     hlabels(87)='UV2Mn';
    hlabels(68)='UAsp';    hlabels(78)='Uqhp';     hlabels(88)='UV2Mp';
    hlabels(69)='Urt';     hlabels(79)='UKINT';    hlabels(89)='UECMT';
    hlabels(70)='Urn';     hlabels(80)='UKINN';    hlabels(90)='UECMN';
    !
    hlabels(91)='UECMP';   hlabels(101)='UERIGT';  hlabels(111)='LNQp';
    hlabels(92)='UROTT';   hlabels(102)='UERIGN';  hlabels(112)='LNpEn';
    hlabels(93)='UROTN';   hlabels(103)='UERIGP';  hlabels(113)='LNpEp';
    hlabels(94)='UROTP';   hlabels(104)='EHFBLN';  hlabels(114)='LNpDn';
    hlabels(95)='USQUJT';  hlabels(105)='EHFB';    hlabels(115)='LNpDp';
    hlabels(96)='USQUJN';  hlabels(106)='LNbet';   hlabels(116)='LNrt';
    hlabels(97)='USQUJP';  hlabels(107)='LNben';   hlabels(117)='LNrn';
    hlabels(98)='UCRANT';  hlabels(108)='LNbep';   hlabels(118)='LNrp';
    hlabels(99)='UCRANN';  hlabels(109)='LNQt';    hlabels(119)='LNrC';
    hlabels(100)='UCRANP'; hlabels(110)='LNQn';    hlabels(120)='LNam2n';
    !
    hlabels(121)='LNam2p';
    hlabels(122)='LNe2n';
    hlabels(123)='LNe2p';
    hlabels(124)='BlEqpN';
    hlabels(125)='BlDEqpN';
    hlabels(126)='BlOvrN';
    hlabels(127)='BlEqpZ';
    hlabels(128)='BlDEqpZ';
    hlabels(129)='BlOvrZ';
    !
    Write(iw1,'((1x,a,2x),6x,660(a,2x))') hlabels
    !
    ! HELP
    !Do i=1,129
    ! Write(iw1,'(1x,i3,a,a)',advance='NO') i,':',trim(hlabels(i))
    !End Do
    ! 1:LBL  2:BLKN  3:BLKZ  4:Jsi  5:JININ  6:A  7:N  8:Z  9:Efn 10:Efp
    ! 11:JEtot 12:Jbett 13:Jbetn 14:Jbetp 15:JQt 16:JQn 17:JQp 18:JpEn 19:JpEp 20:JpDn
    ! 21:JpDp 22:JAsn 23:JAsp 24:Jrt 25:Jrn 26:Jrp 27:Jrc 28:Jht 29:Jhn 30:Jhp
    ! 31:Jqht 32:Jqhn 33:Jqhp 34:JKINt 35:JKINn 36:JKINp 37:JSO 38:JCDIR 39:JCEX 40:JDisn
    ! 41:JDisp 42:JV2Mn 43:JV2Mp 44:JILST 45:JKIND 46:JL 47:JECMPAV1 48:JECMPAV2 49:JECMPAV3 50:JA
    ! 51:JN 52:JZ 53:ITER 54:UEtot 55:Ubett 56:Ubetn 57:Ubetp 58:UQt 59:UQn 60:UQp
    ! 61:Uln 62:Ulp 63:UpEn 64:UpEp 65:UpDn 66:UpDp 67:UAsn 68:UAsp 69:Urt 70:Urn
    ! 71:Urp 72:Urc 73:Uht 74:Uhn 75:Uhp 76:Uqht 77:Uqhn 78:Uqhp 79:UKINT 80:UKINN
    ! 81:UKINP 82:USO 83:UCDIR 84:UCEX 85:UDisn 86:UDisp 87:UV2Mn 88:UV2Mp 89:UECMT 90:UECMN
    ! 91:UECMP 92:UROTT 93:UROTN 94:UROTP 95:USQUJT 96:USQUJN 97:USQUJP 98:UCRANT 99:UCRANN 100:UCRANP
    ! 101:UERIGT 102:UERIGN 103:UERIGP 104:EHFBLN 105:EHFB 106:LNbet 107:LNben 108:LNbep 109:LNQt 110:LNQn
    ! 111:LNQp 112:LNpEn 113:LNpEp 114:LNpDn 115:LNpDp 116:LNrt 117:LNrn 118:LNrp 119:LNrC 120:LNam2n
    ! 121:LNam2p 122:LNe2n 123:LNe2p 124:BlEqpN 125:BlDEqpN 126:BlOvrN 127:BlEqpZ 128:BlDEqpZ 129:BlOvrZ
  End Subroutine thodefh
  !=======================================================================
  !> Allocates arrays depending on the nunber of shells/states and on
  !> the number of Gauss points
  !=======================================================================
  Subroutine thoalloc
    Implicit None
    Integer :: ier,ib,ND
    !
    ! number of int.points
    If(Parity) Then
       ngh=ngh_INI; ngl=ngl_INI; nleg=nleg_INI     !Yesp
    Else
       ngh=2*ngh_INI; ngl=ngl_INI; nleg=nleg_INI   !Nop
    End If
    !
    !nbx=2*n00+1                   ! maximal number of k-blocks
    !ntx=(n00+1)*(n00+2)*(n00+3)/6 ! max.num. p/n levels
    !nzx=n00                       ! maximal nz-quantum number
    !nrx=n00/2+1                   ! maximal nr-quantum number
    !nlx=n00                       ! maximal ml-quantum number
    !ndx=(n00+2)*(n00+2)/4         ! maximal dim. of one k-block
    !nhhdim=number of nonzero HH matrix elements
    !
    nzrlx=(nzx+1)*(nrx+1)*(nlx+1)   ! phy(:,:,nzrlx)
    nghl=ngh*ngl                    ! nghl=ngh*ngl
    nqx=ndx*ndx; nb2x=nbx+nbx; ndx2=ndx+ndx
    ilnqx=ilpjmax*nqx; ilnghl=ilpjmax*nghl
    nhfbx=ndx+ndx; nhfbqx=nhfbx*nhfbx; nkx=ntx; ndxs=ndx*(ndx+1)/2
    !-----------------------------------------
    !Arrays depending on gauss points
    !-----------------------------------------
    If(Allocated(xleg)) Deallocate(xleg,wleg)
    If(nleg.Gt.0) Allocate(xleg(nleg),wleg(nleg))
    If(Allocated(xh)) Deallocate(xh,wh,xl,sxl,wl,vc ,vhbn,vn,vrn,vzn,vdn,vsn,dvn,vhbp,vp,vrp,vzp,vdp,vsp,dvp,  &
                                 vSZFIn,vSFIZn,vSRFIn,vSFIRn,vSZFIp,vSFIZp,vSRFIp,vSFIRp, &
                                 fl,fli,fh,fd,fp1,fp2,fp3,fp4,fp5,fp6, fs1,fs2,fs3,fs4,fs5,fs6, &
                                 wdcor,wdcori,cou,vDHartree,vhart00,vhart01,vhart11)
    Allocate(xh(ngh),wh(ngh),xl(ngl),sxl(ngl),wl(ngl),vc(nghl,nghl))
    Allocate(vhbn(nghl),vn(nghl),vrn(nghl),vzn(nghl),vdn(nghl),vsn(nghl),dvn(nghl),  &
             vhbp(nghl),vp(nghl),vrp(nghl),vzp(nghl),vdp(nghl),vsp(nghl),dvp(nghl),  &
             vSZFIn(nghl),vSFIZn(nghl),vSRFIn(nghl),vSFIRn(nghl),  &
             vSZFIp(nghl),vSFIZp(nghl),vSRFIp(nghl),vSFIRp(nghl))
    Allocate(fl(nghl),fli(nghl),fh(nghl),fd(nghl),fp1(nghl),fp2(nghl),fp3(nghl),  &
             fp4(nghl),fp5(nghl),fp6(nghl),fs1(nghl),fs2(nghl),fs3(nghl),fs4(nghl),  &
             fs5(nghl),fs6(nghl),wdcor(nghl),wdcori(nghl),cou(nghl),vDHartree(nghl,2), &
             vhart00(nghl,nghl),vhart01(nghl,nghl),vhart11(nghl,nghl))
    If(Allocated(aka)) Deallocate(aka,ro,tau,dro,dj,NABLAR,NABLAZ,SZFI,SFIZ,SRFI,SFIR)
    Allocate(aka(nghl,2),ro(nghl,2),tau(nghl,2),dro(nghl,2),dj(nghl,2),  &
             SZFI(nghl,2),SFIZ(nghl,2),SRFI(nghl,2),SFIR(nghl,2),NABLAR(nghl,2),NABLAZ(nghl,2))
    If(Allocated(qfield)) Deallocate(qfield)
    Allocate(qfield(nghl,lambdaMax+1)) ! constraining fields: lambdaMax multipoles + neck
    If(Allocated(MEFFn)) Deallocate(MEFFn,MEFFp)
    Allocate(MEFFn(nghl),MEFFp(nghl))
    If(Allocated(geff_inv)) Deallocate(geff_inv)
    Allocate(geff_inv(nghl,2))
    !-----------------------------------------
    ! Arrays depending on configurations
    !-----------------------------------------
    If(Allocated(rk)) Deallocate(rk,ak,qh,qh1,ql,ql1,nz,nr,nl,ns,npar,id,  &
                                 ia,ikb,ipb,ka,kd,tb,txb,numax,ek,dk,vk,vk1,uk,vkmax,ddc,ddc1,hfb1,lcanon)
    Allocate(rk(nqx,nb2x),ak(nqx,nb2x),qh(0:nzx,1:ngh+1),  &
             qh1(0:nzx,1:ngh+1),ql(0:nrx,0:nlx,1:ngl+1),ql1(0:nrx,0:nlx,1:ngl+1),  &
             nz(ntx),nr(ntx),nl(ntx),ns(ntx),npar(ntx),id(nbx),ia(nbx),ikb(nbx),lcanon(0:nbx,2),  &
             ipb(nbx),ka(nbx,2),kd(nbx,2),tb(ntx),txb(nbx),numax(0:nkx,2),  &
             ek(nkx,2),dk(nkx,2),vk(nkx,2),vk1(nkx,2),uk(nkx,2),vkmax(nkx,2),  &
             ddc(ndx,nkx,2),ddc1(ndx,nkx,2),hfb1(nhfbx,2))
    If(finite_range.or.coulomb_gaussian) Then
       If(Allocated(nrr)) Deallocate(nrr,nll,nss,noo,nzz,nzzx)
       Allocate(nrr(ntx),nll(ntx),nss(ntx),noo(ntx),nzz(ntx,ntx),nzzx(ntx))
    End If
    !-----------------------------------------
    ! HFB Arrays
    !-----------------------------------------
    If(Allocated(erhfb)) Deallocate(erhfb,drhfb,erhfb1,drhfb1)
    Allocate(erhfb(nkx),drhfb(nkx),erhfb1(nkx),drhfb1(nkx))
    If(Allocated(hfb)) Deallocate(hfb,zhfb,evvk,hfbcan,evvkcan)
    Allocate(hfb(ndx2,ndx2),zhfb(ndx2),evvk(ndx2),hfbcan(ndx,ndx),evvkcan(ndx))
    If(Allocated(AN)) Deallocate(AN,ANk,PFIU,PFID,FIU,FID,FIUR,FIDR,FIUD2N,FIDD2N,FIUZ,FIDZ)
    Allocate(AN(nqx),ANk(nqx),PFIU(ndx),PFID(ndx),FIU(ndx),FID(ndx),  &
             FIUR(ndx),FIDR(ndx),FIUD2N(ndx),FIDD2N(ndx),FIUZ(ndx),FIDZ(ndx))
    !-----------------------------------------
    ! Optimal LAPACK storage
    !-----------------------------------------
    If(Allocated(alwork)) Deallocate(alwork,lwork)
#if(USE_ESSL==0)
    ialwork=1+6*ndx+2*ndx**2; ilwork=3+5*ndx;
    Allocate(alwork(ialwork),lwork(ilwork));alwork = 0.0; lwork = 1
#else
    ialwork=0; ilwork=5*ndx;
    Allocate(alwork(1),lwork(ilwork));alwork = 0.0; lwork = 0
#endif
    !-----------------------------------------
    ! Eqp, U,V
    !-----------------------------------------
    If(Allocated(RVqpN)) Deallocate(RVqpN,RVqpP,RUqpN,RUqpP,REqpN,REqpP)
    Allocate(RVqpN(nuv),RVqpP(nuv),RUqpN(nuv),RUqpP(nuv),REqpN(nqp),REqpP(nqp))
    If(Allocated(KpwiP)) Deallocate(KpwiP,KpwiN,KqpN,KqpP)
    Allocate(KpwiN(nqp),KpwiP(nqp),KqpN(nqp),KqpP(nqp))
    If(Allocated(fn_T)) Deallocate(fn_T,fp_T)
    Allocate(fn_T(nqp),fp_T(nqp))
    !-----------------------------------------
    ! PNP ARRAYS: CONF. AND GAUGE ANGLE
    !-----------------------------------------
    If(Allocated(exp1iphy))Deallocate(ropj,taupj,dropj,djpj,akapj,coupj,pjk,  &
         SZFIpj,SFIZpj,SRFIpj,SFIRpj,epj,cpj,ypj,rpj,ddepj,phypj,sinphy,  &
         exp1iphy,exp2iphy,exp1iphym,exp2iphym,xg_norm,xg_phase)
    Allocate(ropj(nghl,ilpjmax,2),taupj(nghl,ilpjmax,2),dropj(nghl,ilpjmax,2),djpj(nghl,ilpjmax,2),akapj(nghl,ilpjmax,2),  &
             SZFIpj(nghl,ilpjmax,2),SFIZpj(nghl,ilpjmax,2),SRFIpj(nghl,ilpjmax,2),SFIRpj(nghl,ilpjmax,2), &
             coupj(nghl,ilpjmax))
    Allocate(epj(ilpjmax,2),pjk(ilpjmax,2),xg_norm(ilpjmax,2),xg_phase(ilpjmax,2))
    Allocate(cpj(nkx,ilpjmax,2),ypj(nkx,ilpjmax,2),rpj(nkx,ilpjmax,2),ddepj(nqx,ilpjmax,nb2x))
    Allocate(phypj(ilpjmax),sinphy(ilpjmax),exp1iphy(ilpjmax),exp2iphy(ilpjmax),exp1iphym(ilpjmax),exp2iphym(ilpjmax))
    !-----------------------------------------
    ! FIELDS INITIALIZATION (NB! optimize)
    !-----------------------------------------
    ro=zero;     tau=zero;    dro=zero;    dj=zero;  aka=zero; rk=zero
    vn=zero;     vsn=zero;    vhbn=zero;   vrn=zero; vzn=zero; vdn=zero;
    vp=zero;     vsp=zero;    vhbp=zero;   vrp=zero; vzp=zero; vdp=zero;
    dvn=zero;    dvp=zero;
    vSFIZn=zero; vSZFIn=zero; vSFIRn=zero; vSRFIn=zero;  vDHartree=zero;
    vSFIZp=zero; vSZFIp=zero; vSFIRp=zero; vSRFIp=zero;
    ! Jason
    If(Allocated(allhfb)) Then
       Do ib=1,oldnb
          Deallocate(allhfb(ib)%arr,allevvk(ib)%arr,allalwork(ib)%arr,alllwork(ib)%arr)
       End Do
       Deallocate (allhfb,allevvk,allalwork,alllwork)
       Deallocate (allIALWORK,allILWORK,allISUPPZ)
    End If
    If (Allocated(allibro)) Deallocate(allibro)
    !
  End Subroutine thoalloc
  !=======================================================================
  !> Allocates arrays depending on the nunber of shells/states and on
  !> the number of Gauss points
  !=======================================================================
  Subroutine realloc
    Implicit None
    Integer :: ier,ib,ND
    !
    !nbx=2*n00+1                   ! maximal number of k-blocks
    !ntx=(n00+1)*(n00+2)*(n00+3)/6 ! max.num. p/n levels
    !nzx=n00                       ! maximal nz-quantum number
    !nrx=n00/2+1                   ! maximal nr-quantum number
    !nlx=n00                       ! maximal ml-quantum number
    !ndx=(n00+2)*(n00+2)/4         ! maximal dim. of one k-block
    !nhhdim=number of nonzero HH matrix elements
    !
    nzrlx=(nzx+1)*(nrx+1)*(nlx+1)   ! phy(:,:,nzrlx)
    nghl=ngh*ngl                    ! nghl=ngh*ngl
    nqx=ndx*ndx; nb2x=nbx+nbx; ndx2=ndx+ndx
    ilnqx=ilpjmax*nqx; ilnghl=ilpjmax*nghl
    nhfbx=ndx+ndx; nhfbqx=nhfbx*nhfbx; nkx=ntx; ndxs=ndx*(ndx+1)/2
    !-----------------------------------------
    ! Arrays depending on configurations
    !-----------------------------------------
    If(Allocated(rk)) Deallocate(rk,ak,qh,qh1,ql,ql1,nz,nr,nl,ns,npar,id,  &
                                 ia,ikb,ipb,ka,kd,tb,txb,numax,ek,dk,vk,vk1,uk,vkmax,ddc,ddc1,hfb1,lcanon)
    Allocate(rk(nqx,nb2x),ak(nqx,nb2x),qh(0:nzx,1:ngh+1),  &
             qh1(0:nzx,1:ngh+1),ql(0:nrx,0:nlx,1:ngl+1),ql1(0:nrx,0:nlx,1:ngl+1),  &
             nz(ntx),nr(ntx),nl(ntx),ns(ntx),npar(ntx),id(nbx),ia(nbx),ikb(nbx),lcanon(0:nbx,2),  &
             ipb(nbx),ka(nbx,2),kd(nbx,2),tb(ntx),txb(nbx),numax(0:nkx,2),  &
             ek(nkx,2),dk(nkx,2),vk(nkx,2),vk1(nkx,2),uk(nkx,2),vkmax(nkx,2),  &
             ddc(ndx,nkx,2),ddc1(ndx,nkx,2),hfb1(nhfbx,2))
    If(finite_range.or.coulomb_gaussian) Then
       If(Allocated(nrr)) Deallocate(nrr,nll,nss,noo,nzz,nzzx)
       Allocate(nrr(ntx),nll(ntx),nss(ntx),noo(ntx),nzz(ntx,ntx),nzzx(ntx))
    End If
    !-----------------------------------------
    ! HFB Arrays
    !-----------------------------------------
    If(Allocated(erhfb)) Deallocate(erhfb,drhfb,erhfb1,drhfb1)
    Allocate(erhfb(nkx),drhfb(nkx),erhfb1(nkx),drhfb1(nkx))
    If(Allocated(hfb)) Deallocate(hfb,zhfb,evvk,hfbcan,evvkcan)
    Allocate(hfb(ndx2,ndx2),zhfb(ndx2),evvk(ndx2),hfbcan(ndx,ndx),evvkcan(ndx))
    If(Allocated(AN)) Deallocate(AN,ANk,PFIU,PFID,FIU,FID,FIUR,FIDR,FIUD2N,FIDD2N,FIUZ,FIDZ)
    Allocate(AN(nqx),ANk(nqx),PFIU(ndx),PFID(ndx),FIU(ndx),FID(ndx),  &
             FIUR(ndx),FIDR(ndx),FIUD2N(ndx),FIDD2N(ndx),FIUZ(ndx),FIDZ(ndx))
    !-----------------------------------------
    ! Optimal LAPACK storage
    !-----------------------------------------
    If(Allocated(alwork)) Deallocate(alwork,lwork)
#if(USE_ESSL==0)
    ialwork=1+6*ndx+2*ndx**2; ilwork=3+5*ndx;
    Allocate(alwork(ialwork),lwork(ilwork));alwork = 0.0; lwork = 1
#else
    ialwork=0; ilwork=5*ndx;
    Allocate(alwork(1),lwork(ilwork));alwork = 0.0; lwork = 0
#endif
    !-----------------------------------------
    ! Eqp, U,V
    !-----------------------------------------
    If(Allocated(RVqpN)) Deallocate(RVqpN,RVqpP,RUqpN,RUqpP,REqpN,REqpP)
    Allocate(RVqpN(nuv),RVqpP(nuv),RUqpN(nuv),RUqpP(nuv),REqpN(nqp),REqpP(nqp))
    If(Allocated(KpwiP)) Deallocate(KpwiP,KpwiN,KqpN,KqpP)
    Allocate(KpwiN(nqp),KpwiP(nqp),KqpN(nqp),KqpP(nqp))
    If(Allocated(fn_T)) Deallocate(fn_T,fp_T)
    Allocate(fn_T(nqp),fp_T(nqp))
    ! Jason
    If(Allocated(allhfb)) Then
       Do ib=1,oldnb
          Deallocate(allhfb(ib)%arr,allevvk(ib)%arr,allalwork(ib)%arr,alllwork(ib)%arr)
       End Do
       Deallocate (allhfb,allevvk,allalwork,alllwork)
       Deallocate (allIALWORK,allILWORK,allISUPPZ)
    End If
    If (Allocated(allibro)) Deallocate(allibro)
    !
  End Subroutine realloc
  !=======================================================================
  !> Setup routine: initializes most variables - including NAMELISTS; set
  !> up the basis and the quadrature grid; computes matrix elements of the
  !> Gogny force (if finite range); prints summary information of the run
  !=======================================================================
  Subroutine preparer(lpr)
    Use HFBTHO_gauss
    Implicit None
    Logical :: lpr
    Integer(ipr) :: iw,l,icount
    Real(pr) :: amas_base
    !
    If(n00.Eq.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)=' STOP: No more nuclei pass to the solver'
       Return
    End If
    !-----------------------------------------
    ! select the symbol of the nucleus
    !-----------------------------------------
    Call nucleus(1,npr(2),nucname)
    If(ierror_flag.Ne.0) Return
    !-----------------------------------------
    ! print headings to screen/'thoout.dat'
    !-----------------------------------------
    If(lpr .And. do_print == 1) Then
       Call heading
       Call print_functional_parameters(skyrme)
       Do iw=lout,lfile
          If(ierror_flag.Ne.0) Return
          If(Print_HFBTHO_Namelist) Then
             Write(iw,'(100(2x,a,f15.8))')
             Write(iw,'(100(2x,a,f15.8))') 'NAMELIST CONTENT (copy/past to hfbtho_NAMELIST.dat and modify)'
             Write(iw,'(100(2x,a,f15.8))') '-------------------------------------------------------------'
             Write(iw,HFBTHO_GENERAL)
             Write(iw,HFBTHO_INITIAL)
             Write(iw,HFBTHO_ITERATIONS)
             Write(iw,HFBTHO_FUNCTIONAL)
             Write(iw,HFBTHO_PAIRING)
             Write(iw,HFBTHO_CONSTRAINTS)
             Write(iw,HFBTHO_BLOCKING)
             Write(iw,HFBTHO_PROJECTION)
             Write(iw,HFBTHO_FEATURES)
             Write(iw,HFBTHO_TDDFT)
             Write(iw,HFBTHO_NECK)
             Write(iw,HFBTHO_TEMPERATURE)
             Write(iw,HFBTHO_DEBUG)
             Write(iw,HFBTHO_RESTORATION)
        End If
       End Do
    End If
    !-----------------------------------------
    ! particle number as real variable
    !-----------------------------------------
    tz(1)=Real(npr(1),Kind=pr); tz(2)=Real(npr(2),Kind=pr)
    If(fragment_properties) Then
       tz=tz_fragments(1:2)
    End If
    amas=tz(1)+tz(2)
    ! In case of blocking on top of an even-even core, force the same oscillator length
    ! for the odd nucleus as for the core. At this point, particle number npr() have already
    ! been increased/decreased by 1 unit to acocunt for the odd nucleus, so after this trick
    ! amas_base becomes equal to the mass number of the even-even core.
    amas_base = amas
    If(.Not.odd_noBlock) Then
       If(nkblo_INI(1,1).Gt.0.And.npr(1).Ne.2*(npr(1)/2)) amas_base = amas_base-one
       If(nkblo_INI(1,1).lt.0.And.npr(1).Ne.2*(npr(1)/2)) amas_base = amas_base+one
    End If
    If(.Not.odd_noBlock) Then
       If(nkblo_INI(2,1).Gt.0.And.npr(2).Ne.2*(npr(2)/2)) amas_base = amas_base-one
       If(nkblo_INI(2,1).lt.0.And.npr(2).Ne.2*(npr(2)/2)) amas_base = amas_base+one
    End If
    drhoi=zero
    !-----------------------------------------
    ! default combinations
    !-----------------------------------------
    chargee2=e2charg
    coex=-chargee2*(three/pi)**p13; cex=-0.750_pr*coex
    !-----------------------------------------
    ! hbzero from forces [hqc**2/(two*amu)]
    !-----------------------------------------
    if(.not.hb0_charge_dependent) then
       hbzeron = hbzero
       hbzerop = hbzero
    endif
    hb0=hbzero;   If (use_cm_cor) hb0=hb0*(one-one/amas)
    hb0n=hbzeron; If (use_cm_cor) hb0n=hb0n*(one-one/amas)
    hb0p=hbzerop; If (use_cm_cor) hb0p=hb0p*(one-one/amas)
    !-----------------------------------------
    ! basis parameter q
    !-----------------------------------------
    beta0=q; q=Exp((3.0_pr*Sqrt(5.0_pr/(16.0_pr*pi)))*beta0)
    !-----------------------------------------
    ! basis parameters b0,bp,bz
    !-----------------------------------------
    If(b0.Le.zero) Then
       ! define oscillator frequency from default with empirical factor 1.2,
       ! and set length accordingly
       r00=r0*amas_base**p13; r02=r00**2; r04=r02**2
       hom=41.0_pr*amas_base**(-p13)*r0
       b0=Sqrt(two*hbzero/hom)
    Else
       ! define oscillator frequency from user-defined length, and set default
       ! empirical factor accordingly
       hom=hqc**2/(amn*b0**2)
       r0=(hom/41.0_pr)*amas_base**(p13)
       r00=r0*amas_base**p13; r02=r00**2; r04=r02**2
    End If
    ! Provides a mechanism to overwrite HO lengths with input, say from HFODD
    If (bz_INI > zero .And. bp_INI > zero) Then
        bp=bp_INI;bz=bz_INI;b0=b0_INI
    Else
        bp=b0*q**(-one/6.0_pr); bz=b0*q**(one/3.0_pr)
    End If
    bpp=bp*bp
    !-----------------------------------------
    ! constraint in terms of beta
    !-----------------------------------------
    ty20=Sqrt(5.0_pr/pi)*hom/b0**2/two
    !-----------------------------------------
    ! projection: number of grid points
    !-----------------------------------------
    keypjn=Max(1,keypjn); keypjp=Max(1,keypjp); ilpjn=keypjn; ilpjp=keypjp; ilpj2=ilpjn*ilpjp;
    ilpjmax=Max(ilpjn,ilpjp)
    If(iproj.Eq.0) Then
       npr1pj=npr(1); npr2pj=npr(2)
    Else
       npr1pj=npr(1)+npr1pj; npr2pj=npr(2)+npr2pj
    End If
    !-----------------------------------------
    ! blocking window
    !-----------------------------------------
    !pwiblo=Min(Max(25.0_pr/Sqrt(Real(npr(1)+npr(2),Kind=pr)),2.0_pr),8.0_pr)
    pwiblo=1.0_pr
    !-----------------------------------------
    ! THO
    !-----------------------------------------
    ass=zero; iasswrong=0
    !-----------------------------------------
    ! iterations
    !-----------------------------------------
    etot=zero; varmas=zero; rms=zero; ept=-two; del=one; alast=-seven; siold=one
    varmasNZ=zero; pjmassNZ=zero; ass=zero; skass=zero
    !---------------------------------------------------------
    ! statistics to screen('lout')/file('lfile')
    !---------------------------------------------------------
    If(lpr) Then
       If(iLST.Eq.0)  Then         ! HFB+HO case only
          iLST1=0
       Else                        ! HFB+THO case
          If(iLST.Gt.0) Then       ! HFB+THO  only
             iLST1=1
          Else                     ! HFB+THO after HFB+HO
             iLST1=0
          End If
       End If
    End If
    If(lpr .And. do_print == 1) Then
       Do iw=lout,lfile
          Write(iw,*)
          Write(iw,'(a)')             '  ---------------------------------------'
          Write(iw,'(a)')             '        Characteristics of the run       '
          Write(iw,'(a)')             '  ---------------------------------------'
          Write(iw,'(a,i5)')          '  Output file ................: ',lfile
          Write(iw,'(a,2x,a2,i4)')    '  Nucleus ....................: ',nucname,npr(1)+npr(2)
          Write(iw,'(a,i5)')          '  Number of HO shells ........: ',n00
          Write(iw,'(a,f20.14)')      '  HO length b0 (fm) ..........: ',b0
          Write(iw,'(a,f8.3,a,f8.3)') '  Basis deformation ..........:  beta0=',beta0,' q=',q
          Write(iw,'(a,5(1x,e15.8))') '  HO: b0,1/b0,bp,bz,q ........: ',b0,one/b0,bp,bz,q
          Write(iw,'(a,3(1x,e15.8))') '  h**2/(2m_n), cmc, e**2 .....: ',hbzeron,hb0n,chargee2
          Write(iw,'(a,3(1x,e15.8))') '  h**2/(2m_p), cmc, e**2 .....: ',hbzerop,hb0p,chargee2
          Write(iw,'(a,2(1X,e15.8))') '  hom=f*41.0_pr*A^{-1/3}, f...: ',hom,r0
          If(automatic_basis) Then
             Write(iw,'(a)')       '    Adjusted basis is ........:   ON'
          End If
          If(iLST.Eq.0)  Then         ! HFB+HO case only
             Write(iw,'(a)')          '  THO basis is ...............:  OFF'
          Else                        ! HFB+THO case
             Write(iw,'(a)')          '  THO basis is ...............:   ON'
             If(iLST.Gt.0) Then       ! HFB+THO  only
                If(inin.Gt.0) Then
                   ierror_flag=ierror_flag+1
                   ierror_info(ierror_flag)=' Stop: Forbidden iLST>0, inin>0 combination.'
                   Return
                End If
                Write(iw,'(a)')       '    THO parameters from tholst.wel'
             Else                     ! HFB+THO after HFB+HO
                Write(iw,'(a)')       '    HFB+THO after a HFB+HO run    '
             End If
          End If
          Write(iw,'(a,i5)')          '  Maximal number of iterations: ',maxi
          Write(iw,'(a,f6.3)')        '  Initial mixing parameter ...: ',xmix
          If(inin.Eq.1)  Then
             Write(iw,'(a)')          '  Initial w.f. ...............:  from spherical scratch'
          End If
          If(inin.Eq.2)  Then
             Write(iw,'(a)')          '  Initial w.f. ...............:  from prolate scratch'
          End If
          If(inin.Eq.3)  Then
             Write(iw,'(a)')          '  Initial w.f. ...............:  from oblate scratch'
          End If
          If(inin.Lt.0) Then
             Write(iw,'(a)')          '  Initial wave functions from :  tape'
          End If
          Write(iw,'(a,3x,a)')        '  Energy functional ..........: ',skyrme
          If(finite_range) Then
             Write(iw,'(a,f6.2,a)')   '    with a finite-range central force'
          End If
          If(icou.Eq.-4) Write(iw,'(a)')    '    direct Coulomb by sum of Gaussians, exchange Coulomb with Slater approximation'
          If(icou.Eq.-3) Write(iw,'(a)')    '    direct, exchange and pairing Coulomb by sum of Gaussians'
          If(icou.Eq.-2) Write(iw,'(a)')    '    direct and exchange Coulomb by sum of Gaussians'
          If(icou.Eq.-1) Write(iw,'(a)')    '    direct Coulomb force only by sum of Gaussians'
          If(icou.Eq. 0) Write(iw,'(a)')    '    without Coulomb forces'
          If(icou.Eq. 1) Write(iw,'(a)')    '    direct Coulomb by substitution method'
          If(icou.Eq. 2) Write(iw,'(a)')    '    direct Coulomb by substitution method, exchange Coulomb with Slater approximation'
          Write(iw,'(a,f6.2,a)')          '  Pairing cutoff energy ......: ',pwi, ' MeV'
          If(kindhfb.Lt.0) Then
             Write(iw,'(a)')          '  Lipkin-Nogami procedure is .:   ON'
          Else
             Write(iw,'(a)')          '  Lipkin-Nogami procedure is .:  OFF'
          End If
          If((ilpjn-1.Eq.0.And.ilpjp-1.Eq.0).Or.iproj.Eq.0) Then
             Write(iw,'(a)')          '  PAV procedure is ...........:  OFF'
          Else
             Write(iw,'(a)')          '  PAV procedure is ...........:   ON'
             Write(iw,'(a,2i5)')       '    Number of gauge points....: ',keypjn,keypjp
          End If
          If(icstr.Eq.0) Then
             Write(iw,'(a)')          '  Constraint calculation is ..:  OFF'
          Else
             Write(iw,'(a)')          '  Constraint calculation is ..:   ON'
             icount=0
             Do l=1,8
                If(Abs(lambda_active(l)).Gt.0) Then
                   icount=icount+1
                   Write(iw,'(a,i1,a,i1,a,f8.3)') '    Constraint ',icount,' .............: lambda=',l, &
                                                  ' Ql=',multRequested(l)
                End If
             End Do
             If(neck_constraints) Then
                icount=icount+1
                Write(iw,'(a,i1,a,a,f8.3)') '    Neck       ',icount,' .............: lambda=0', &
                                                  ' Ql=',neckRequested
             End If
          End If
          If(iparenti(1).Ne.0) Then
             Write(iw,'(a)')          '  Neutron blocking is ........:   ON'
          End If
          If(iparenti(2).Ne.0) Then
             Write(iw,'(a)')          '  Proton blocking is .........:   ON'
          End If
          If(switch_on_temperature) Then
             Write(iw,'(a,f6.2,a)')   '  Temperature T ..............: ',temper,' MeV'
          Else
             Write(iw,'(a,f6.2)')     '  Temperature T ..............:   0.00 MeV'
          End If
          If(pairing_regularization) Then
             Write(iw,'(a,f6.2,a)')   '  Pairing regularization is ..:   ON'
          Else
             Write(iw,'(a,f6.2,a)')   '  Pairing regularization is ..:  OFF'
          End If
          If(collective_inertia) Then
             Write(iw,'(a,f6.2,a)')   '  Collective inertia is ......:   ON'
          Else
             Write(iw,'(a,f6.2,a)')   '  Collective inertia is ......:  OFF'
          End If
          If(fission_fragments) Then
             Write(iw,'(a,f6.2,a)')   '  Fission fragments are ......:   ON'
          Else
             Write(iw,'(a,f6.2,a)')   '  Fission fragments are ......:  OFF'
          End If
          Write(iw,'(a,i3)')          '  Restart indicator ..........: ',inin
          If(nbroyden.Eq.0) Then
             Write(iw,'(a,i3)')       '  Linear mixing ..............: ',nbroyden
          Else
             Write(iw,'(a,i3)')       '  Broyden mixing (#iterations): ',nbroyden
          End If
       End Do
    End If
    !-----------------------------------------
    ! BASIS, GAUSS POINTS, HOWF
    !-----------------------------------------
    Call gfv                      ! factorials
    If(ierror_flag.Ne.0) Return
    Call base0(lpr)               ! basis space (calculate configurational space)
    If(ierror_flag.Ne.0) Return
    Call thoalloc                 ! global allocation
    If(ierror_flag.Ne.0) Return
    Call gausspoints              ! GAUSS mesh points
    If(ierror_flag.Ne.0) Return
    Call recompute_coulomb_expansion ! expansion of the Coulomb potential on Gaussians
    If(ierror_flag.Ne.0) Return
    Call base(lpr)                ! oscillator configurations (set up quantum numbers)
    If(ierror_flag.Ne.0) Return
    Call gaupol(lpr)              ! basis wf at gauss mesh points
    If(ierror_flag.Ne.0) Return
    If(finite_range.or.coulomb_gaussian) Then
       Call gogny_matrix_elements
    endif
    If(ierror_flag.Ne.0) Return
    !
  End Subroutine preparer
  !====================================================================
  !> Defines and stores the weights and nodes of integration in the
  !> dimensionless coordinates \f$ \xi\f$ and \f$ \eta \f$. Defines the
  !> scaling transformation for the THO basis.
  !====================================================================
  Subroutine coordinateLST(lpr)
    Use HFBTHO_THO, Only: f01234
    Implicit None
    Logical :: lpr
    Integer(ipr) :: i,il,ih
    If(iLST1.Eq.0) Then
       ! HO-basis
       Do il=1,ngl
          Do ih=1,ngh
             i=ih+(il-1)*ngh
             fh(i)=bz*xh(ih)
             fl(i)=bp*Sqrt(xl(il))
             wdcor(i)=pi*wh(ih)*wl(il)*bz*bp*bp
             wdcori(i)=one/wdcor(i)
          End Do
       End Do
    Else
       ! THO basis
       Call f01234(.False.)
       If(ierror_flag.Ne.0) Return
    End If
    !
    Call optHFBTHO ! optimal HO/THO combinations
    If(ierror_flag.Ne.0) Return
    !
  End Subroutine coordinateLST
  !====================================================================
  !> Self-consistent loop
  !====================================================================
  Subroutine iter(lpr, history)
    Implicit None
    Logical, intent(in) :: lpr, history
    Real(pr)       :: time5,assprn,delln(2)
    Real(pr), Save :: time
    Integer(ipr)   :: iw,it,ite
    integer(ipr)   :: wct1,wct2,countrate,countmax
    Logical :: iter_conv, do_broyden, exit_loop
    !---------------------------------------------------
    ! print to screen('lout')/thoout.dat('lfile')
    !---------------------------------------------------
    If(do_print == 1) Then
       Do iw=lout,lfile
          If(iLST.Eq.0) Then
             Write(iw,'(a,f7.3,4(a,i3),a)')  &
                  '  |HFB+HO> iterations(b0=',b0,', Nsh=',n00,  &
                  ', inin=',inin,', N=',npr(1),', Z=',npr(2),')...'
          Else
             If(iLST1.Eq.0.Or.iasswrong(3).Ne.0) Then
                If(iasswrong(3).Ne.0) Then
                   Write(iw,'(a,f7.3,a,i3,a)')  &
                        '  |HFB+THO substituted by HFB+HO> iterations (b0=',  &
                        b0,', Nsh=',n00,')...'
                Else
                   Write(iw,'(a,f7.3,a)')'  towards |hfb+tho> iterations...'
                   Write(iw,'(a,f7.3,a)')
                   Write(iw,'(a,f7.3,a,i3,a)')  &
                        '  |Preliminary HFB+HO> iterations (b0=',b0,', Nsh=',n00,')...'
                End If
             Else
                If(itass.Eq.1) Then
                   Write(iw,'(2(a,f7.3),a,i3,a)')  &
                        '  |HFB+THO> iterations(b0=',b0,', neutron density decay=',  &
                        decay,', Nsh=',n00,')...'
                Else
                   Write(iw,'(2(a,f7.3),a,i3,a)')  &
                        '  |HFB+THO> iterations(b0=',b0,', proton density decay=',  &
                        decay,', Nsh=',n00,')...'
                End If
             End If
          End If
          Write(iw,'(2x,130("-"))')
          Write(iw,'(20(a))') '  i','          si ','    mix ','  beta',    &
         &    '      Etot ','      A ','      rn','      rp ','        En', &
         &    '      Dn','      Ep','      Dp','        Ln  ','    Lp ',    &
         &    '    time  time(Gog.)'
          Write(iw,'(2x,130("-"))')
       End Do
    End If
    !---------------------------------------------------------------------
    ! main HFB iteration loop
    !---------------------------------------------------------------------
    iError_in_HO=0; iError_in_THO=0; time=0.0_pr; time5=0.0_pr; iter_fail=0
    iter_conv=.false.; do_broyden=.true.; exit_loop=.false.
    If(maxi==0) exit_loop=.true.
    Do ite=mini,maxi+2 ! +1 for no broyden, +1 for final printout
       !-------------------------------------------------
       ! Print current iteration (iiter=0 for initial hfbmatrix)
       !-------------------------------------------------
       If(lpr .and. iiter.Ge.mini .or. iiter.Eq.0) Then
          assprn=ass(1); If(assprn.Gt.ass(2)) assprn=-ass(2) ! protons come with '-'
          delLN=del+frdel; If(kindhfb.Lt.0) delLN=del+frdel+ala2 ! LN case
          If(do_print == 1) Then
             Do iw=lout,lfile
                Write(iw,'(i4,a1,1x,f12.8,f5.2,f7.3,f13.6,1x,f6.1,2(f8.3)," | ",4(f8.3)," | ",4(f8.3))') &
                      iiter,bbroyden,si,xmix,bet,etot,varmas,rms(1),rms(2),ept(1)+frept(1),delLN(1), &
                      ept(2)+frept(2),delLN(2),alast(1),alast(2),time,wct_gogny
             End Do
          End If
       End If
       !-------------------------------------------------
       ! Exit loop (1 extra iteration without broyden for storage)
       !-------------------------------------------------
       ! Exit here to set correct variables in printout above
       If(exit_loop) Exit
       ! Checking if we have had too many failures in determining the Fermi level (a bad sign)
       If(iter_fail > iter_fail_max) Then
          ierror_flag=ierror_flag+1
          Write(ierror_info(ierror_flag),'(" STOP: in iter(), iter_fail=",i3," exceeds iter_fail_max=",i3)')  iter_fail,iter_fail_max
          Return
       End If
       ! Converged
       If(iiter.Ge.mini+1.And.si.Lt.epsi) Then
          iter_conv = .true.
          exit_loop = .true.
       End If
       ! Slow to converge and lambda>0
       If(iiter.Ge.1000 .And. (alast(1).Gt.zero.Or.alast(2).Gt.zero)) Then
          iter_conv = .false.
          exit_loop = .true.
       End If
       ! Iteration limit reached
       If(iiter+1==maxi+1) exit_loop = .true.
       ! Final iteration without broyden, for storage
       If(exit_loop) Then
          If(history) Then
             Exit
          Else
             do_broyden = .false.
          End If
       End If
       !-------------------------------------------------
       ! Begin the iteration
       !-------------------------------------------------
       Call system_clock(wct1,countrate,countmax)
       If(do_broyden) iiter=ite
       !-------------------------------------------------
       ! HFBDIAG
       !-------------------------------------------------
       Do it=itmin,itmax
          Call hfbdiag(it,0) ! HFB diagonalization with minimal canonical
          If(ierror_flag.Ne.0) Return
       End Do
       !-------------------------------------------------
       ! EXPECT (->DENSIT, COULOMB), FIELD, GAMDEL
       !-------------------------------------------------
       Call expect(.False.) ! expectation values
       If (numberCons.Gt.0) Call getLagrange(ite) ! new Lagrange parameters for constraints
       If(ierror_flag.Ne.0) Return
       Call field ! new fields
       If(ierror_flag.Ne.0) Return
       Call gamdel(.false.,do_broyden) ! new HFB matrix
       If(ierror_flag.Ne.0) Return
       !-------------------------------------------------
       ! Dumping control (old linear mixing)
       !-------------------------------------------------
       If(do_broyden) Then
          If(si.Lt.siold) Then
             xmix = Min(xmax,xmix * 1.10_pr);  ! old value: 1.13
          Else
             xmix = 0.5_pr ! Mario: xmix0
          End If
       End If
       siold=si
       !-------------------------------------------------
       ! Time per iteration
       !-------------------------------------------------
       call system_clock(wct2,countrate,countmax)
       time=(wct2-wct1)/real(countrate,kind=pr)
       time5=time5+time
    End Do
    ! Increase start for subsequent calls to iter
    mini = iiter + 1
    !
    If(iter_conv) Then
       !-------------------------------------------------
       ! Solution converged
       !-------------------------------------------------
       If(iLST1.Eq.0) Then
          iError_in_HO=0
       Else
          iError_in_THO=0
       End If
       If(.Not.lpr) Then ! (deprecated?)
          If(do_print == 1) Then
             Do iw=lout,lfile
                Write(iw,'(2x,130("-"))')
                Write(iw,'("  *   iteration interrupted after",i4," steps   si=",f17.10," ho=",i3," tho=",i3)') &
                      iiter,si,iError_in_HO,iError_in_THO
                Write(iw,'(2x,130("-"))')
                Write(iw,'(a,f8.3,a)') '  Total CPU time=',time5/60.0_pr,' minutes'
             End Do
          End If
       End If
       If(do_print == 1) Then
          Do iw=lout,lfile
             Write(iw,'(2x,130("-"))')
             Write(iw,'("  *   iteration converged   after",i4," steps   si=",f17.10," ho=",i3," tho=",i3)') &
                 iiter,si,iError_in_HO,iError_in_THO
             Write(iw,'(2x,130("-"))')
             Write(iw,'(a,f8.3,a)') '  Total CPU time=',time5/60.0_pr,' minutes'
          End Do
       End If
    Else
       !-------------------------------------------------
       ! Solution interrupted due to iterations limit
       !-------------------------------------------------
       If(iLST1.Eq.0) Then
          iError_in_HO=-1
       Else
          iError_in_THO=-1
       End If
       If(do_print == 1) Then
          Do iw=lout,lfile
             Write(iw,'(2x,130("-"))')
             Write(iw,'("  *   iterations limit interrupt after",i4," steps   si=",f17.10," ho=",i3," tho=",i3)') &
                 iiter,si,iError_in_HO,iError_in_THO
             Write(iw,'(2x,130("-"))')
             Write(iw,'("  Total CPU time=",f8.3," minutes")') time5/60.0_pr
          End Do
       End If
    End If
  End Subroutine iter
  !====================================================================
  !> Block-diagonalization of the HFB matrix, loop over the Fermi level and
  !> definition of the density matrix and pairing tensor in configuration space.
  !====================================================================
  Subroutine hfbdiag(it,icanon)
#if(USE_OPENMP==1)
    Use omp_lib
#endif
    Implicit None
    Logical :: lpr_pwi,norm_to_improve
    Character(Len=1) :: char1,char2,char3
    Integer(ipr) :: iw,it,i0,icanon,ibiblo,ier,i,j,k,k0,kl,lc,ib,nd,  &
                    nhfb,n1,n2,kaib,m,ndk,nd1,nd2,kdib,k1,k2,id1,id2, &
                    n12,n21,ntz,nhhph,nhhpp,ibro,ibroib,i_uv,i_eqp,jj,&
                    tid,IL,IU,NUMFOU,jlwork,jalwork,ldw,ldi,ii,counter
    Real(pr) :: al,al2,emin,hla,dla,pn,eqpe,ela,enb,enb1,ekb, &
                s1,s2,s3,alnorm,sitest,fac1,fac2,fT,exponent, &
                VL,VU,ABSTOL,buffer
    Integer(ipr), Pointer :: KpwiPo(:),KqpPo(:)
    Real(pr), Pointer     :: EqpPo(:),VqpPo(:),UqpPo(:),f_T(:)
    Integer(ipr), Allocatable :: ISUPPZ(:),lwork_p(:)
    Real(pr), Allocatable :: alwork_p(:),eigenv(:),eigenf(:,:),hfbmat(:,:)
    Real(pr), Allocatable :: H_diag(:,:),Wmat(:,:),bufmat(:,:),diag_elems(:),Ub(:),Vb(:)
    Real(pr), External :: DLAMCH
    !
    If (debug_solver.Ge.1) Call get_CPU_time('hfbdiag',0)
    !
    If(it.Eq.1) Then
       EqpPo=>REqpN; VqpPo=>RVqpN; UqpPo=>RUqpN; KpwiPo=>KpwiN; KqpPo=>KqpN; f_T=>fn_T
    Else
       EqpPo=>REqpP; VqpPo=>RVqpP; UqpPo=>RUqpP; KpwiPo=>KpwiP; KqpPo=>KqpP; f_T=>fp_T
    End If
    KpwiPo=0; KqpPo=0; f_T=zero
    !
    nhhph=(it-1)*nhhdim; nhhpp=(it+1)*nhhdim
    If(.Not. Allocated(allhfb)) Then
       oldnb = nb ! for destroying data structures in next computation
       Allocate(allhfb((nb)),allevvk((nb)),allalwork((nb)),alllwork((nb)),allIALWORK(nb),allILWORK(nb))
       Allocate(allISUPPZ(nb))
       Do ib=1,nb
          nhfb=2*id(ib)
          Allocate(allhfb(ib)%arr(1:nhfb,1:nhfb))
          Allocate(allevvk(ib)%arr(1:nhfb))
          !jalwork=1+6*nhfb+2*nhfb**2; allIALWORK(ib)=jalwork  ! DSYEVD
          !jlwork=3+5*nhfb; allILWORK(ib)=jlwork               ! DSYEVD
#if(USE_ESSL==0)
          jalwork=26*nhfb; allIALWORK(ib)=jalwork  ! DSYEVR
          jlwork=10*nhfb; allILWORK(ib)=jlwork     ! DSYEVR
#else
          jalwork=8*nhfb; allIALWORK(ib)=jalwork   ! DSYEVX
          jlwork=5*nhfb; allILWORK(ib)=jlwork      ! DSYEVX
#endif
          Allocate(allalwork(ib)%arr(1:jalwork))
          Allocate(alllwork(ib)%arr(1:jlwork))
          Allocate(allISUPPZ(ib)%arr(1:2*nhfb)) ! DSYEVR
       End Do
    End If
    !
    ABSTOL=2.0_pr*DLAMCH('S')
    !
    If (.Not. Allocated(allibro)) Then
       Allocate(allibro(1:NB))
       allibro(1)=0
       Do ib=2,NB
          allibro(ib) = allibro(ib-1) + (ID(ib-1)*(ID(ib-1)+1)/2)
       End Do
    End If
    !
    !------------------------------------------------------------------
    ! Loop the internal normalization
    !------------------------------------------------------------------
    !sitest=Max(Min(0.10_pr,si*0.010_pr),0.000010_pr)
    sitest=Min(0.10_pr,si*0.010_pr)
    norm_to_improve=.True.; inner(it)=-1; sumnz(it)=one
    Do While(norm_to_improve)
       !
       inner(it)=inner(it)+1
       !
       If(Abs(sumnz(it)).Lt.sitest.Or.inner(it).Eq.20) norm_to_improve=.False.
       !
       sumnz(it)=zero; entropy(it)=zero; v2min(it)=one; Dispersion(it)=zero
       !
       kl=0; emin=1000.0_pr; al=ala(it)
       !
       !------------------------------------------------------------------
       !                          BLOCKING
       !------------------------------------------------------------------
       If(iparenti(it).Eq.0) blomax(it)=0
       blo123d(it)=0; blok1k2d(it)=0; blocanon(it)=0;
       ibiblo=bloblo(keyblo(it),it)
       !------------------------------------------------------------------
       ! Runs over Omega blocks
       !------------------------------------------------------------------
       i_uv=0; i_eqp=0
       lc=0; lcanon(0,it)=0; klmax=0
!$OMP Parallel Default(None) &
!$OMP& SHARED(nb,id,ia,it,nbx,allibro,brin,allhfb,allevvk,lout, &
!$OMP&        allALWORK,allLWORK,allIALWORK,allILWORK,nhhph,nhhpp,al, &
!$OMP&        zhfb,ndx2,allISUPPZ,ABSTOL) &
!$OMP& PRIVATE(ib,nd,nhfb,i0,m,ibro,n1,nd1,n2,nd2,hla,dla,ier,tid,char1,char2, &
!$OMP&         NUMFOU,IL,IU,VL,VU,eigenf,eigenv,hfbmat,ISUPPZ,alwork_p,lwork_p, &
!$OMP&         ldw,ldi,char3,jalwork,jlwork)
#if(USE_OPENMP==1)
       tid = OMP_GET_THREAD_NUM()
#endif
!$OMP DO SCHEDULE(DYNAMIC)
       Do ib=1,nb
          nd=id(ib); nhfb=nd+nd; i0=ia(ib); m=ib+(it-1)*nbx; ibro=allibro(ib)
#if(USE_ESSL==0)
          jalwork=26*nhfb; jlwork=10*nhfb
#else
          jalwork=8*nhfb; jlwork=5*nhfb
#endif
          allhfb(ib)%arr(1:nhfb,1:nhfb)=0.0_pr; allevvk(ib)%arr(1:nhfb)=0.0_pr
          allALWORK(ib)%arr(1:jalwork)=0.0_pr; allLWORK(ib)%arr(1:jlwork)=0;
          allISUPPZ(ib)%arr(1:2*nhfb)=0
          !------------------------------------------------------------------
          !  hfb-matrix
          !------------------------------------------------------------------
          Allocate(hfbmat(nhfb,nhfb))
          Do n1=1,nd
             nd1=n1+nd
             Do n2=1,n1
                nd2=n2+nd; ibro=ibro+1
                hla=brin(nhhph+ibro); dla=brin(nhhpp+ibro)
                hfbmat(n1,n2)=hla;    hfbmat(nd2,n1)=dla
                hfbmat(nd1,n2)=dla;   hfbmat(nd1,nd2)=-hla
             End Do
             hfbmat(n1,n1)  =hfbmat(n1,n1)  -al
             hfbmat(nd1,nd1)=hfbmat(nd1,nd1)+al
          End Do
          char1='V'; char2='I'; char3='L'; NUMFOU=0
          VL=0.0_pr ;VU=0.0_pr; IL=1; IU=nhfb; ldw=allIALWORK(ib); ldi=allILWORK(ib)
          Allocate(eigenv(nhfb)); eigenv(1:nhfb)=0.0_pr
          Allocate(eigenf(nhfb,nhfb)); eigenf(1:nhfb,1:nhfb)=0.0_pr
          Allocate(alwork_p(ldw)); Allocate(lwork_p(ldi))
#if(USE_ESSL==0)
          ier=0; Allocate(ISUPPZ(2*nhfb))
          Call DSYEVR(char1,char2,char3,nhfb,hfbmat,nhfb,VL,VU,IL,IU,ABSTOL,NUMFOU, &
                      eigenv,eigenf,nhfb,ISUPPZ,alwork_p,ldw,lwork_p,ldi,ier)
#else
          ier=0; Allocate(ISUPPZ(nhfb))
          Call DSYEVX('V','A','L',nhfb,hfbmat,nhfb,VL,VU,IL,IU,ABSTOL,NUMFOU, &
                      eigenv,eigenf,nhfb,alwork_p,ldw,lwork_p,ISUPPZ,ier)
#endif
          allevvk(ib)%arr(1:nhfb) = eigenv(1:nhfb)
          allhfb(ib)%arr(1:nhfb,1:nhfb) = eigenf(1:nhfb,1:nhfb)
          If(ier.NE.0) Then
#if(USE_OPENMP==1)
             If(tid == 0) Then
                Write(lout,'("The algorithm failed to compute eigenvalues")')
                Write(lout,'("Fermi level: ",f20.14," block size: ",i4," vector sizes:",2i8)') al,nd,nhhph,nhhpp
             End If
             Write(lout,'("Thread ",i4," block ib=",i4," - Error code: ",i4)') tid,ib,ier
#else
             Write(lout,'("The algorithm failed to compute eigenvalues")')
             Write(lout,'("Block ib=",i4," - Error code: ",i4)') tid,ib,ier
#endif
          End If
          Deallocate(eigenf,eigenv,hfbmat,ISUPPZ,alwork_p,lwork_p)
       End Do ! ib
!$OMP End Do
!$OMP End Parallel
       Do ib=1,NB
          nd=id(ib); nhfb=nd+nd; i0=ia(ib); m=ib+(it-1)*nbx; ibro=allibro(ib)
          !------------------------------------------------------------------
          !                          BLOCKING
          !------------------------------------------------------------------
          ! Specific blocking:
          !   keyblo==0 still (not set in handle_blocking), set blocking config here
          !   NOTE!:
          !   requested_blocked_level overwrites keyblo=1 of even-even candidate
          !   list with blocked QP. The result is what is stored in odd binaries.
          ! Automatic blocking:
          !   keyblo/=0 here (set in handle_blocking), so this is skipped
          !------------------------------------------------------------------
          If(iiter.Le.1.And.inner(it).Eq.0) Then ! run if iiter=0,1
             ! Even-even binary: okay to lookup request level from candidates
             If(blocking_never_done(it)) Then
                If(iparenti(it).Ne.0.And.keyblo(it).Eq.0) Then
                   Call requested_blocked_level(ib,it)
                   If(ierror_flag.Ne.0) Return
                   ibiblo=bloblo(keyblo(it),it)
                End If
             ! Odd binary: use previously requested level (CAN'T change level)
             Else
                If(iparenti(it).Ne.0.And.keyblo(it).Eq.0) Then
                   keyblo(it)=1
                   ibiblo=bloblo(keyblo(it),it)
                End If
             End If
          End If
          ! Finding out actual qp to block
          k0=0
          If(ibiblo.Eq.ib) Then
              If(iiter.Le.1.And.inner(it).Eq.0) Then ! run if iiter=0,1
                ! blocked level as in the even-even nucleus
                k0=blo123(keyblo(it),it); ndk=k0+nd
                Do n2=1,nd
                   nd2=n2+nd
                   hfb1(n2,it)=allhfb(ib)%arr(n2,ndk)    !U
                   hfb1(nd2,it)=allhfb(ib)%arr(nd2,ndk)  !V
                End Do
                ! number of states in the block to be tested
                blocross(it)=Min(blomax(it)+10,nd)
             End If
             ! overlap between new and old blocked levels
             s3=zero
             Do n1=1,blocross(it)
                ndk=n1+nd; s1=zero
                Do n2=1,nd
                   nd2=n2+nd
                   s1=s1+Abs(hfb1(nd2,it)*allhfb(ib)%arr(nd2,ndk)) !VV
                   s1=s1+Abs(hfb1(n2,it)*allhfb(ib)%arr(n2,ndk))   !UU
                End Do
                ! Maximum overlap
                If(s1.Gt.s3) Then
                   s3=s1; k0=n1
                End If
             End Do
             ! Update for possibly new state being blocked
             blo123d(it)=k0
             blo123(keyblo(it),it)=k0
             If(.Not.norm_to_improve) Then
                ! find maximal HO component
                ndk=k0+nd
                s1=zero
                Do n1=1,nd
                   nd1=n1+nd
                   hfb1(n1,it)=allhfb(ib)%arr(n1,ndk); hfb1(nd1,it)=allhfb(ib)%arr(nd1,ndk)
                   s2=Max(s1,Abs(allhfb(ib)%arr(n1,ndk)),Abs(allhfb(ib)%arr(nd1,ndk)))
                   ! Maximum overlap
                   If(s2.Gt.s1) Then
                      s1=s2; i=n1+i0  ! labels in k[k1,k2] numbering
                   End If
                End Do
                ! print blocked state
                If(do_print == 1) Then
                   Do iw=lout,lfile
                      Write(iw,'(4x,a,2(a,i3),2x,3(a,1x,f12.8,1x),(i3,a,i3,1x),a)')  &
                           protn(it),' Blocking: block=',ib,  &
                           ' state=',k0,  &
                           ' Eqp=',allevvk(ib)%arr(k0+nd),  &
                           ' Dqpe=',allevvk(ib)%arr(k0+nd)-eqpmin(it),  &
                           ' Ovlp=',s3, keyblo(it),'/',blomax(it), tb(i)
                   End Do
                End If
                ! ieresbl=6, 'BLKN','BLKZ'
                ereslbl(it)=tb(i)
                If(it.Eq.1) Then
                   ! 'BlEqpN','BlDEqpN','BlOvrN'
                   eresbl(1)=allevvk(ib)%arr(k0+nd); eresbl(2)=allevvk(ib)%arr(k0+nd)-eqpmin(it); eresbl(3)=s1
                Else
                   ! 'BlEqpZ','BlDEqpZ','BlOvrZ'
                   eresbl(4)=allevvk(ib)%arr(k0+nd); eresbl(5)=allevvk(ib)%arr(k0+nd)-eqpmin(it); eresbl(6)=s1
                End If
             End If
          End If
          !------------------------------------------------------------------
          ! Run over all qp states k in the block
          !------------------------------------------------------------------
          kaib=kl
          Do k=1,nd
             ndk=k+nd
             ! referent spectra
             pn=zero
             Do i=1,nd
                hla=allhfb(ib)%arr(i+nd,ndk)**2; pn=pn+hla
             End Do
             ! Blocking
             If(k.Eq.k0) Then
                If(nkblo_INI(it,1).Ne.0) blocking_never_done(it) = .False.
                n1=k0+nd
                Do i=1,nd
                   hla=allhfb(ib)%arr(i+nd,n1)**2; dla=allhfb(ib)%arr(i,n1)**2; pn=pn-half*(hla-dla)
                End Do
             End If
             eqpe=allevvk(ib)%arr(nd+k); ela=eqpe*(one-two*pn)
             enb=ela+al;                 ekb=Sqrt(Abs(eqpe**2-ela**2))
             i_eqp=i_eqp+1
             !------------------------------------------------------------------
             ! cut-off condition: energy pwi + Fermi cut-off function
             !------------------------------------------------------------------
                                                              exponent=Huge(1.0_pr)
             If(Abs(100.0_pr*(enb-pwi)).Lt.Log(Huge(1.0_pr))) exponent=Exp(100.0_pr*(enb-pwi))
             If(basis_HFODD) Then
                lpr_pwi=enb.Le.pwi !jacek sharp cut off for hfodd
                !lpr_pwi=enb.Le.pwi.Or.Abs(one/(one+exponent)).Gt.cutoff_tol
             Else
                lpr_pwi=enb.Le.pwi.Or.Abs(one/(one+exponent)).Gt.cutoff_tol
             End If
             if(finite_range) then
                lpr_pwi = .true.
             endif
             !------------------------------------------------------------------
             ! Remember the whole qp solution
             !------------------------------------------------------------------
             If(.Not.norm_to_improve) Then
                EqpPo(i_eqp)=eqpe              ! Eqp_k
                If(lpr_pwi) KqpPo(kl+1)=i_eqp  ! below pwi otherwise zero
                If(lpr_pwi) KpwiPo(kl+1)=i_uv  ! below pwi otherwise zero
                ! Attention: The q.p. spinors are filled out even if Eqp > Ecut
                !            In other words, UqpPo and VqpPO (to be assigned to
                !            RVqpN and RUqpN contain *all* spinors irrespective
                !            cut-off.
                Do n2=1,nd
                   nd2=n2+nd; i_uv=i_uv+1
                   UqpPo(i_uv)=allhfb(ib)%arr(n2,ndk)   ! U_ak
                   VqpPo(i_uv)=allhfb(ib)%arr(nd2,ndk)  ! V_ak
                End Do
             End If
             !------------------------------------------------------------------
             ! Define Fermi-Dirac occupations
             !------------------------------------------------------------------
             fT=zero
             If(switch_on_temperature.And.temper.Gt.1.D-12) Then
                fT = half*(one-Tanh(half*eqpe/temper))
                ! factor two comes from K>0 states only
                buffer = zero
                If(fT.Gt.zero.And.fT.Lt.one) Then
                   buffer = two*fT*Log(fT) + two*(one-fT)*Log(one-fT)
                End If
                entropy(it) = entropy(it) - buffer
                f_T(i_eqp) = fT
             End If
             !------------------------------------------------------------------
             ! Pairing window
             !------------------------------------------------------------------
             If(lpr_pwi) Then
                kl=kl+1                                        !number of active states
                ! Update for possibly new state being blocked
                If(k0.Eq.k) blok1k2d(it)=kl                    !blocking: dynamic #: k[k1,k2] numbering
                If(k0.Eq.k) blok1k2(keyblo(it),it)=kl
                If((eqpe.Le.emin).And.(pn.Gt.0.0001_pr)) Then  !to avoid unocc at magic numbers
                   emin=eqpe; alnorm=pn                        !min qpe and its occupation
                End If
                erhfb(kl)=enb; drhfb(kl)=ekb; uk(kl,it)=pn     !ref.s.p. energies, deltas, occupancies
                sumnz(it)=sumnz(it)+two*pn+two*(one-two*pn)*fT !internal normalization
             End If
          End Do ! End k
          !
          If(.Not.norm_to_improve) Then
             ! Debugging: For blocking, store the Bogoliubov matrix before it is overwritten,
             !            swap U and V to do the blocking and check that the resulting matrix
             !            is still unitary
             If(debug_solver.Ge.2 .And. ibiblo.Eq.ib .And. Abs(si).Le.1.e-9_pr) Then
                If(Allocated(Wmat)) Deallocate(Wmat)
                Allocate(Wmat(nhfb,nhfb)); Wmat=zero
                Do n1=1,nhfb
                   Do n2=1,nhfb
                      Wmat(n1,n2) = allhfb(ib)%arr(n1,n2)
                   End Do
                End Do
                ! blocking
                Allocate(Ub(nd),Vb(nd))
                i = nd+k0
                Vb = Wmat(1:nd,i); Ub = Wmat(nd+1:nhfb,i)
                Wmat(1:nd,i) = Ub; Wmat(nd+1:nhfb,i) = Vb
                i = nd-k0+1
                Ub = Wmat(1:nd,i); Vb = Wmat(nd+1:nhfb,i)
                Wmat(1:nd,i) = Vb; Wmat(nd+1:nhfb,i) = Ub
                Deallocate(Ub,Vb)
                ! Orthonormalization of W after blocking prescription
                Allocate(H_diag(1:nhfb,1:nhfb)); H_diag=zero
                Call dgemm('t','n',nhfb,nhfb,nhfb,one,Wmat,nhfb,Wmat,nhfb,zero,H_diag,nhfb)
                counter=0
                Do k=1,nhfb
                   Do n1=1,nhfb
                      If(n1.Ne.k.And.Abs(H_diag(n1,k)).Gt.1.e-10_pr) Then
                         Write(6,'("n1=",i4," n2=",i4," H_diag=",f20.14)') n1,k,H_diag(n1,k)
                         counter=counter+1
                      End If
                   End Do
                End Do
                ! Print diagonal elements
                Write(6,'("ib=",i3," counter=",i12)') ib,counter
                Do k=1,nhfb
                   Write(6,'("ib=",i2," k=",i4," H_kk=",f20.14," Eqp=",f12.6)') ib,k,H_diag(k,k),allevvk(ib)%arr(k)
                End Do
                Deallocate(H_diag)
             End If
             !------------------------------------------------------------------
             !  Density matrices (computed only when norm_to_improve = False,
             !  e.g., when the particle number is conserved within 10^-5)
             !------------------------------------------------------------------
             kdib=kl-kaib; ka(ib,it)=kaib; kd(ib,it)=kdib
             k1=kaib+1; k2=kaib+kdib
             eqpe=0.0_pr
             Do n2=1,nd
                Do n1=n2,nd
                   s1=zero; s2=zero
                   If(k1.Le.k2) Then
                      Do k=k1,k2
                         ! temperature
                         fac1 = one; fac2 = zero
                         If(switch_on_temperature) Then
                            ii=KqpPo(k); fac1=one-f_T(ii); fac2=f_T(ii)
                         End If
                         nd1=KpwiPo(k)+n1; nd2=KpwiPo(k)+n2
                         s1=s1+VqpPo(nd1)*fac1*VqpPo(nd2)+UqpPo(nd1)*fac2*UqpPo(nd2)
                         s2=s2+UqpPo(nd1)*fac1*VqpPo(nd2)+VqpPo(nd1)*fac2*UqpPo(nd2) &
                              +VqpPo(nd2)*fac1*UqpPo(nd1)+UqpPo(nd2)*fac2*VqpPo(nd1)
                      End Do
                      s1=two*s1; s2=half*s2 ! two:due to m-projection, half:due to symmetrization
                      ! blocking
                      If(ibiblo.Eq.ib) Then
                         If(blok1k2d(it).Eq.0) Then
                            ierror_flag=ierror_flag+1
                            ierror_info(ierror_flag)=' STOP: No blocking candidate found (1)!!!'
                            Return
                         End If
                         i=blok1k2d(it); id1=KpwiPo(i)+n1; id2=KpwiPo(i)+n2
                         s1=s1-VqpPo(id1)*VqpPo(id2)+UqpPo(id1)*UqpPo(id2)
                         s2=s2-half*(UqpPo(id1)*VqpPo(id2)+VqpPo(id1)*UqpPo(id2))
                      End If
                   End If
                   n12=n1+(n2-1)*nd; n21=n2+(n1-1)*nd
                   rk(n12,m)=s1; rk(n21,m)=s1     !  V V'
                   ak(n12,m)=-s2; ak(n21,m)=-s2   !- U V', ak=half*(pairing density)
                   ! Store the density matrix \rho
                   hfbcan(n1,n2)=s1; allhfb(ib)%arr(n1,n2)=s1
                End Do !n1
             End Do !n2
             ! Debugging: For blocking, test the structure of the HFB matrix in the
             !            quasiparticle basis after blocking
             If(debug_solver.Ge.2 .And. ibiblo.Eq.ib .And. Abs(si).Le.1.e-9_pr) Then
                ! Reconstruct the HFB matrix in HO basis
                Allocate(hfbmat(nhfb,nhfb))
                Do n1=1,nd
                   nd1=n1+nd
                   Do n2=1,n1
                      nd2=n2+nd; ibro=ibro+1
                      hla=brin(nhhph+ibro); dla=brin(nhhpp+ibro)
                      hfbmat(n1,n2)=hla;    hfbmat(nd2,n1)=dla
                      hfbmat(nd1,n2)=dla;   hfbmat(nd1,nd2)=-hla
                   End Do
                   hfbmat(n1,n1)  =hfbmat(n1,n1)  -al
                   hfbmat(nd1,nd1)=hfbmat(nd1,nd1)+al
                End Do
                ! Symmetrize
                Do n1=1,nhfb
                   Do n2=n1+1,nhfb
                      hfbmat(n1,n2)=hfbmat(n2,n1)
                  End Do
                End Do
                ! Compute the HFB matrix in the qp basis by multiplying by the Bogoliubov matrix
                Allocate(H_diag(1:nhfb,1:nhfb),bufmat(1:nhfb,1:nhfb)); H_diag=zero
                Call dgemm('n','n',nhfb,nhfb,nhfb,one,hfbmat,nhfb,Wmat,nhfb,zero,bufmat,nhfb)
                Call dgemm('t','n',nhfb,nhfb,nhfb,one,Wmat,nhfb,bufmat,nhfb,zero,H_diag,nhfb)
                ! Sort diagonal elements
                Allocate(diag_elems(1:nhfb))
                counter=0
                Do k=1,nhfb
                   diag_elems(k) = H_diag(k,k)
                   Do n1=1,nhfb
                      If(n1.Ne.k.And.Abs(H_diag(n1,k)).Gt.1.e-10_pr) Then
                         Write(6,'("n1=",i4," n2=",i4," H_diag=",f20.14)') n1,k,H_diag(n1,k)
                         counter=counter+1
                      End If
                   End Do
                End Do
                !Call ord(nhfb,diag_elems)
                ! Print diagonal elements
                Write(6,'("ib=",i3," counter=",i12)') ib,counter
                Do k=1,nhfb
                   Write(6,'(i2,i4,2f20.14)') ib,k,diag_elems(k),allevvk(ib)%arr(k)
                End Do
                Deallocate(hfbmat,H_diag,bufmat,diag_elems)
                If(Allocated(Wmat)) Deallocate(Wmat)
             End If
             !------------------------------------------------------------------
             ! Canonical basis
             !------------------------------------------------------------------
             If(k1.Le.k2) Then
                Call Canonical(it,icanon,k2,k1,nd,i0,lc,ib,ibiblo,m,ibro)
                If(ierror_flag.Ne.0) Return
             End If
             lcanon(ib,it)=lc
             !
          End If
          !
       End Do !ib
       !
       If(kl.Eq.0) Then
          ierror_flag=ierror_flag+1
          ierror_info(ierror_flag)=' STOP: kl=zero, no states below pwi!!!'
          Return
       End If
       If(iparenti(it).Ne.0.And.ibiblo.Eq.0) Then
          ierror_flag=ierror_flag+1
          ierror_info(ierror_flag)='STOP: No blocking candidate found (2)!!!'
          Return
       End If
       eqpmin(it)=emin; klmax(it)=kl; sumnz(it)=sumnz(it)-tz(it)
       !------------------------------------------------------------------
       ! Lambda search
       !------------------------------------------------------------------
       Call ALambda(al,it,kl)
       If(ierror_flag.Ne.0) Return
       If(keyblo(it).Eq.0) Then
          ala(it)=al
       Else
          ala(it)=ala(it)+0.50_pr*(al-ala(it))
       End If
       ! NB! 'alast' instead of 'al' at small pairing
       alast(it)=al
       If(Abs(ept(it)+frept(it)).Lt.0.0001_pr.And.(.Not.switch_on_temperature)) Then
          ntz=Int(tz(it)+0.1_pr); ntz=ntz/2
          Do k=1,kl
             drhfb(k)=erhfb(k)
          End Do
          Call ord(kl,drhfb)
          alast(it)=drhfb(ntz)  !last bound s.p. energy
       End If
       !------------------------------------------------------------------
       ! THO asymptotic decay
       !------------------------------------------------------------------
       ! density asymptotic decay \rho(r)->Exp(-ass(it)*r)
       ! ass(it)=2*Sqrt((E_min-\lambda)/((A-1)/A)*hbar**2/(2*m)))
       al2=zero
       If(kindhfb.Lt.0) Then
          al2=al+two*ala2(it)*(one-two*alnorm) ! al=al+two*ala2(it)
       End If
       if(it.eq.1) then
          al2=(emin-al2)/hb0n
       else
          al2=(emin-al2)/hb0p
       endif
       ! wrong asymptotic
       iasswrong(it)=0; If(al2.Le.zero) iasswrong(it)=1; ass(it)=two*Sqrt(Abs(al2))
       !
    End Do ! While(norm_to_improve)
    !
    If (debug_solver.Ge.1) Call get_CPU_time('hfbdiag',1)
    !
  End Subroutine hfbdiag
  !=======================================================================
  !> Determination of the Fermi energy based on the value of the particle
  !> number and the current density matrix
  !=======================================================================
  Subroutine ALambda(al,it,kl)
    Implicit None
    Integer(ipr) :: it,i,k,kl,icze,lit,ntz,iw
    Real(pr), Save :: fm7=1.0e-7_pr,fm10=1.0e-10_pr
    Real(pr) :: al,vh,xinf,xsup,esup,ez,dez,dfz,dvh,y,a,b,einf,absez,sn
    Real(pr) :: fT,dfT
    Real(pr), Pointer :: f_T(:)
    !-------------------------------------------------
    ! Fermi-Dirac occupations
    !-------------------------------------------------
    If(switch_on_temperature) Then
       If(it.Eq.1) Then
          f_T=>fn_T
       Else
          f_T=>fp_T
       End If
    End If
    !-------------------------------------------------
    ! Chemical potential without pairing
    !-------------------------------------------------
    If(CpV0(it-1).Eq.zero.and.trim(skyrme).ne.'D1S') Then
       ntz=Int(tz(it)+0.1_pr); ntz=ntz/2
       Do k=1,kl
          drhfb(k)=erhfb(k)
       End Do
       Call ord(kl,drhfb)
       If(ntz.Lt.kl) Then
          al=half*(drhfb(ntz)+drhfb(ntz+1))
       Else
          al=drhfb(ntz)+0.001_pr
       End If
       Return
    End If
    !-------------------------------------------------
    ! Chemical potential with pairing
    !-------------------------------------------------
    xinf=-1000.0_pr; xsup=1000.0_pr; esup=one; icze=0
    Do lit=1,500
       sn=zero;dez=zero;dfz=zero
       Do i=1,kl
          vh=zero; dvh=zero; fT=zero; dfT=zero
          y=erhfb(i)-al; a=y*y+drhfb(i)**2; b=Sqrt(a)
          !
          If(switch_on_temperature.And.temper.Gt.1.e-12_pr) Then
             fT =half*(one-Tanh(half*b/temper))
             dfT=y/b/temper*fT*(one-fT)
             f_T(i)=fT
          Else
             fT =zero
             dfT=zero
          End If
          !
          If(b.Gt.zero)  vh=half*(one-y/b)
          !
          !If(b.Lt.fm7.And.icze.Eq.1) vh=-einf/(esup-einf) !no pairing
          If(vh.Lt.1.e-12_pr) vh = zero
          If((vh-one).Gt.1.e-12_pr)  vh = one
          If(b.Gt.zero) dvh=half*drhfb(i)**2/(a*b)         ! D[ez,al](i)
          ! blocking
          If(i.Eq.blok1k2d(it)) Then
             vh=half; dvh=zero
          End If
          sn=sn+two*vh+two*(one-two*vh)*fT
          dez=dez+two*(one-two*fT)*dvh
          dfz=dfz+two*(one-two*vh)*dfT   ! D[ez,al]
       End Do
       ez=sn-tz(it); absez=Abs(ez)/tz(it)
       dez=dez+dfz
       !-------------------------------------------------
       ! Correcting bounds
       !-------------------------------------------------
       If(ez.Lt.zero) Then
          xinf=Max(xinf,al); einf=ez
       Else
          xsup=Min(xsup,al); esup=ez
       End If
       If(lit.Eq.1) Then
          If(absez.Le.0.10_pr) Then
             al=al-ez
          Else
             al=al-0.10_pr*Sign(one,ez)
          End If
       Else
          al=al-ez/(dez+1.e-20_pr)                      ! Newton method
       End If
       If(xsup-xinf.Lt.fm7) icze=1                      ! low/upp close
       If(al.Lt.xinf.Or.al.Gt.xsup) al=half*(xinf+xsup) ! mean upp/low
       If(absez.Le.fm10) Return
    End Do
    !-------------------------------------------------
    ! Low accuracy warning
    !-------------------------------------------------
    If(do_print == 1) Then
       Do iw=lout,lfile
         Write(iw,'(a,2(e12.5,2x),a,2(2x,f8.4),a,i2)') ' Low accuracy=',sn,ez,' for N,Z=',tz,' it=',it
       End Do
    End If
    iter_fail=iter_fail+1 ! Incrementing iteration failure counter
  End Subroutine ALambda
  !=========================================================================================
  !> Determination of the canonical basis by diagonalization of the density matrix
  !=========================================================================================
  Subroutine Canonical(it,icanon,k2,k1,nd,i0,lc,ib,ibiblo,m,ibroib)
    Use HFBTHO_utilities
    Use HFBTHO
    Implicit None
    Integer(ipr) :: it,i0,icanon,ibiblo,i,iw,k,kk,lc,ib,nd,n1,n2,m,nd1,k1,k2,n12,ier
    Integer(ipr) :: nhhph,nhhpp,ibro,ibroib,il,iu,NUMFOU
    Integer(ipr), Allocatable :: ifail(:)
    Real(pr) :: s1,s2,vx,h1,d1,h2,d2,ddn1,ddn2,vl,vu,ABSTOL
    Real(pr), Allocatable :: hh(:,:),de(:,:),ewavef(:,:)
    Real(pr), Pointer     :: EqpPo(:),VqpPo(:),UqpPo(:)
    Integer(ipr), Pointer :: KpwiPo(:),KqpPo(:)
    Integer(ipr), Allocatable :: ISUPPZ(:)
    Real(pr), Allocatable :: eigenv(:),eigenf(:,:)
    Real(pr), External :: DLAMCH
    !
    If(it.Eq.1) Then
       EqpPo=>REqpN; VqpPo=>RVqpN; UqpPo=>RUqpN; KpwiPo=>KpwiN; KqpPo=>KqpN
    Else
       EqpPo=>REqpP; VqpPo=>RVqpP; UqpPo=>RUqpP; KpwiPo=>KpwiP; KqpPo=>KqpP
    End If
    !
#if(USE_ESSL==0)
    If(Allocated(alwork)) Deallocate(alwork,lwork)
    ialwork=26*nd; ilwork=10*nd
    Allocate(ALWORK(ialwork),LWORK(ilwork)); ALWORK = 0.0_pr; LWORK = 1
#else
    If(Allocated(alwork)) Deallocate(alwork,lwork)
    ialwork=8*nd; ilwork=5*nd;
    Allocate(ALWORK(ialwork),LWORK(ilwork)); ALWORK = 0.0_pr; LWORK = 0
#endif
    !
    ABSTOL=2.0_pr*DLAMCH('S')
    !
    !------------------------------------------------------
    ! Vanishing pairing => just taking the HF states
    !------------------------------------------------------
    If(Abs(ept(it)).Lt.0.0001_pr.And.(.Not.switch_on_temperature)) Then
       Do k=1,nd
          kk=k1+k-1; lc=lc+1                              ! total number of the canonical states
          ddc(1:nd,lc,it)=zero; vk(lc,it)=zero            ! zeros: nd could be larger then k2-k1+1
          If(kk.Gt.k2) Cycle
          vx=zero
          Do i=1,nd
             h1=VqpPo(KpwiPo(kk)+i)**2; vx=vx+h1
          End Do
          If (vx.Le.zero) vx=zero                         ! roundoff errors
          If (vx.Ge.one ) vx=one
          Do i=1,nd
             If(vx.Ge.half) Then
                ddc(i,lc,it)=VqpPo(KpwiPo(kk)+i)          ! (ph) s.p. orbitals in conf.space
             Else
                ddc(i,lc,it)=UqpPo(KpwiPo(kk)+i)          ! (ph) s.p. orbitals in conf.space
             End If
          End Do
          Dispersion(it)=Dispersion(it)+four*vx*(one-vx)  ! internal P/N Dispersion
          If(Abs(vx-half).Le.v2min(it)) Then
             v2min(it)=Abs(vx-half); v2minv(it)=vx        ! divergent condition
             lcc=lc
          End If
          vk(lc,it)=vx                                    ! (ph) s.p. occupations v^2
          !------------------------------------------------------
          ! RESU only
          !------------------------------------------------------
          If(icanon.Ne.0) Then
             ek(lc,it)=EqpPo(KqpPo(kk))*(one-two*vx)+ala(it)      ! (ph) s.p. energies
             dk(lc,it)=zero                                       ! (ph) s.p. deltas
          End If
       End Do !k
    Else
       !------------------------------------------------------
       ! Pairing => calculate canonical basis
       !------------------------------------------------------
#if(USE_ESSL==0)
       VL=0.0_pr; VU=0.0_pr; IL=1; IU=nd; NUMFOU=0
       Allocate(ISUPPZ(2*nd))
       Allocate(eigenv(nd)); eigenv(1:nd)=0.0_pr
       Allocate(eigenf(nd,nd)); eigenf(1:nd,1:nd)=0.0_pr
       ier=0; Call DSYEVR('V','A','L',nd,hfbcan,ndx,VL,VU,IL,IU,ABSTOL,NUMFOU, &
                          eigenv,eigenf,nd,ISUPPZ,ALWORK,ialwork,LWORK,ilwork,ier)
       evvkcan(1:nd) = eigenv(1:nd)
       hfbcan(1:nd,1:nd) = eigenf(1:nd,1:nd)
       Deallocate(eigenv,eigenf,ISUPPZ)
       !ier=0; Call DSYEVD('V','L',nd,hfbcan,ndx,evvkcan,ALWORK,ialwork,LWORK,ilwork,ier)
#else
       vl=0.0_pr; vu=0.0_pr; il=1; iu=1; m=0; abstol=0.0_pr
       Allocate(ifail(nd),ewavef(ndx,ndx))
       ier=0; Call DSYEVX('V','A','L',nd,hfbcan,ndx,vl,vu,il,iu,abstol,m, &
                          evvkcan,ewavef,ndx,ALWORK,ialwork,LWORK,ifail,ier)
       hfbcan(1:ndx,1:ndx) = ewavef(1:ndx,1:ndx)
       Deallocate(ifail,ewavef)
#endif
       ! bug in LAPACK
       If(ier.Gt.0) Then
          Do iw=lout,lfile
             Write(iw,*) 'FATAL ERROR CONDITION IN CANONICAL DSYEVR, ier=',ier,'(RECOVERED)'
          End Do
          Do n2=1,nd
             Do n1=n2,nd
                vx=allhfb(ib)%arr(n1,n2)
                hfbcan(n2,n1)=vx; hfbcan(n1,n2)=vx
             End Do
          End Do
          Call sdiag(ndx,nd,hfbcan,evvkcan,hfbcan,zhfb,+1)
       End If
       !------------------------------------------------------
       ! Eigenvalues and wavefunctions
       !------------------------------------------------------
       Do k=1,nd
          lc=lc+1 ! total number of canonical states
          Do i=1,nd
             ddc(i,lc,it)=hfbcan(i,k) ! (ph) canon orbitals in conf.space
          End Do
          vx=evvkcan(k)*half
          If (vx.Le.zero) vx=zero ! roundoff errors
          If (vx.Ge.one ) vx=one
          ! blocking
          If(ibiblo.Eq.ib.And.vx.Gt.0.49_pr.And.vx.Le.0.51_pr) blocanon(it)=lc
          Dispersion(it)=Dispersion(it)+four*vx*(one-vx) ! internal P/N Dispersion
          If(Abs(vx-half).Le.v2min(it)) Then
             v2min(it)=Abs(vx-half); v2minv(it)=vx ! divergent condition
             lcc=lc
          End If
          vk(lc,it)=vx ! (ph) canon occupations v^2
          !------------------------------------------------------
          ! RESU only
          !------------------------------------------------------
          If(icanon.Ne.0) Then
             ! canonical energies and deltas (no physical meaning in PNP)
             nhhph=(it-1)*nhhdim; nhhpp=(it+1)*nhhdim
             Allocate(hh(nd,nd),de(nd,nd))
             ibro=ibroib
             Do n1=1,nd
                Do n2=1,n1
                   ibro=ibro+1
                   vx=brin(nhhph+ibro); hh(n2,n1)=vx; hh(n1,n2)=vx
                   vx=brin(nhhpp+ibro); de(n2,n1)=vx; de(n1,n2)=vx
                End Do
             End Do
             h1=zero; d1=zero
             Do n2=1,nd
                h2=zero; d2=zero
                Do n1=1,nd
                   ddn1=hfbcan(n1,k)
                   h2=h2+ddn1*hh(n1,n2)
                   d2=d2+ddn1*de(n1,n2)
                End Do
                ddn2=hfbcan(n2,k)
                h1=h1+h2*ddn2
                d1=d1+d2*ddn2
             End Do
             ek(lc,it)=h1 ! (ph) canon s.p. energies
             dk(lc,it)=d1 ! (ph) canon s.p. deltas
             Deallocate(hh,de)
          End If
          !
       End Do !k
    End If
    !------------------------------------------------------
    ! RESU only
    !------------------------------------------------------
    If(icanon.Ne.0) Then
       !------------------------------------------------------
       ! Find maximal HO components of all qp states
       !------------------------------------------------------
       Do k=k1,k2
          s1=zero
          Do n1=1,nd
             nd1=nd+n1
             s2=Max(s1,Abs(VqpPo(KpwiPo(k)+n1)),Abs(UqpPo(KpwiPo(k)+n1)))
             If(s2.Gt.s1) Then
                s1=s2
                vkmax(k,it)=s1    ! maximal overlap
                numax(k,it)=n1+i0 ! its number in k[k1,k2] numbering
             End If
          End Do
       End Do
       !------------------------------------------------------
       ! Searching for possible blocking candidates
       !------------------------------------------------------
       If(iparenti(1).Eq.0 .and. iparenti(2).Eq.0) Then
          n1=0
          Do k=k1,k2
             n1=n1+1
              ! Search within |(1-2*N)*Eqpe| lover than 'pwiblo'
              ! The levels number n1 is 1,2,3,... for the given block ([123] numbering)
             If(Abs(EqpPo(KqpPo(k))-eqpmin(it)).Le.pwiblo) Then
                blomax(it)=blomax(it)+1                         ! blocked state #, maximel # of block candidates
                If(blomax(it).Gt.bloall) Then
                   ierror_flag=ierror_flag+1
                   ierror_info(ierror_flag)='Too many blocking candidates! Increase bloall and run again'
                   Return
                End If
                bloblo(blomax(it),it)=ib ! Omega block where to block
                blo123(blomax(it),it)=n1 ! index of the blocking candidate within the block
                blok1k2(blomax(it),it)=k ! index of the blocking candidate amongst all quasiparticles
                bloqpdif(blomax(it),it)=Abs(EqpPo(KqpPo(k))-eqpmin(it)) ! difference in energy between current qp and one with lowest energy
             End If
          End Do
       End If
    End If
    !
  End Subroutine Canonical
  !=======================================================================
  !> Displays results at convergence: single particle energies, densities, fields
  !=======================================================================
  Subroutine resu(irecord)
    Implicit None
    Integer(ipr) :: it,iw,ib,im,m,nd,k,k0,k1,k2,j,n,imax,nhfb,irecord,iexit
    Real(pr)     :: sum,eqpe,pn,ela,enb,ek0,vk0,ekk,delb,ovmax,s,uuvv,  &
         dk0,skk,summ(4),vvs,vvc,enjacek
    Real(pr), Pointer     :: EqpPo(:),VqpPo(:),UqpPo(:)
    Integer(ipr), Pointer :: KpwiPo(:),KqpPo(:)
    Logical :: l_write_wel
    !
    iexit=0
    !--------------------------------------------
    ! Write binary file
    !--------------------------------------------
    ! Write if: nucleus is even-even, odd+specific blocking, odd+blocking deactivated
    l_write_wel = ( (write_hel .And. fragment_properties) .Or.&
        (write_hel .And. .Not.fragment_properties .And. .Not.odd_noBlock .And. &
        (blocking_never_done(1) .And. blocking_never_done(2) .Or. (blocking_mode(1).Eq.1 .Or. blocking_mode(2).Eq.1) ) ) .Or. &
        (write_hel .And. .Not.fragment_properties .And. odd_noBlock) )
    !
    ! Write pre-diagonalization data (fields, etc., without broyden)
    ! This ensures that, reading from binary, the 1st iteration (or here,
    ! if maxi=0) exactly reproduces what is written to thoout.dat
    If(l_write_wel) Then
       Call inout(2,iexit)
       If(ierror_flag.Ne.0) Return
    End If
    !
    ! Get the SOLUTION for fields in binary (blocking cands, QP energies, etc.)
    Do it=itmin,itmax
       Call hfbdiag(it,1) ! Maximal canonical
       If(ierror_flag.Ne.0) Return
    End Do

    !
    ! Write post-diagonalization data. This is mainly for blocking candidates,
    ! since handling of odd particle number makes reproducing this after reading
    ! in a binary difficult (c.f. preparer), but also applies to e.g. FT occupations.
    If(l_write_wel) Then
       Call inout(3,iexit)
       If(ierror_flag.Ne.0) Return
    End If

    !--------------------------------------------
    ! Printing quasiparticle states
    !--------------------------------------------
    Do it=itmin,itmax
       If(it.Eq.1) Then
          EqpPo=>REqpN; VqpPo=>RVqpN; UqpPo=>RUqpN; KpwiPo=>KpwiN; KqpPo=>KqpN
       Else
          EqpPo=>REqpP; VqpPo=>RVqpP; UqpPo=>RUqpP; KpwiPo=>KpwiP; KqpPo=>KqpP
       End If
       !
       If(Print_Screen .And. do_print == 1) Then
          iw=lfile
          Write(iw,200) tit(it)
          Write(iw,*) ' eqp(k) -> q.p. energy '
          Write(iw,*) ' e(k)   -> referent s.p. energy '
          Write(iw,*) ' p(k)   -> occ.probability '
          Write(iw,*) ' del(k) -> referent s.p. gap '
          Write(iw,*) ' fermi energy alast=',alast(it)
          Write(iw,'(a,a)')  &
               '  #k  block#    eqp(k)     e(k)       (1-2N)E      decay        p(k)',  &
               '        del(k)    overl      labels'
200     Format(//,' #quasiparticle energies ',a,/,1x,32('-'))
       End If
       sum=zero
       Do ib=1,nb
          nd=id(ib); im=ia(ib); m=ib+(it-1)*nbx; nhfb=nd+nd
          k1=ka(ib,it)+1
          k2=ka(ib,it)+kd(ib,it)
          If(k1.Le.k2) Then
             Do k=k1,k2                                 ! print active states only
                pn=uk(k,it)                             ! qp probabilities
                j=k
                If(pn.Gt.-1.d-14) Then                   ! print If signIficant pn
                   ! main oscillator component
                   ovmax=vkmax(k,it)                    ! maximal overlap
                   imax=numax(k,it)                     ! its number
                   ! printing
                   eqpe=EqpPo(KqpPo(k))                 ! qp energies
                   if(it.eq.1) then
                      skk=two*Sqrt(Abs(eqpe-ala(it))/hb0n)  ! qp decay
                   else
                      skk=two*Sqrt(Abs(eqpe-ala(it))/hb0p)  ! qp decay
                   endif
                   ela=eqpe*(one-two*pn)
                   enb=ela+ala(it)                      ! ref. s.p. energies
                   delb=Sqrt(Abs(eqpe**2-ela**2))       ! ref. s.p. delta
                   sum=sum+two*pn                       ! particle number
                   If(Print_Screen .And. do_print == 1) Then
                      iw=lfile
                      Write(iw,201) k,ib,eqpe,enb,(one-two*pn)*eqpe,skk,pn,delb,ovmax,tb(imax)
201                   Format(i4,2x,i3,1x,f12.6,f12.6,f12.6,f12.6,2x,f12.8,  &
                             2(2x,f7.4),' ',a13)
                   End If
                End If
             End Do
          End If
       End Do !ib
       !--------------------------------------------
       ! Printing canonical single particle states
       !--------------------------------------------
       If(Print_Screen .And. do_print == 1) Then
          iw=lfile
          Write(iw,'(a,i4,a,i4)')  &
               '#all active are ',j,' q.p. states out of ',nt
          Write(iw,'(a,f6.1)') '#since the cut off is pwi=',pwi
          Write(iw,'(3a,f6.1)')'#check: number of ',tit(it),'=',sum
          Write(iw,100) tit(it)
          Write(iw,*) ' labels -> {2*omega}{parity}[nn=nz+2*nr+nl,nz,nl]'
          Write(iw,*) ' cqpe   -> canonical q.p. energies'
          Write(iw,*) ' ce     -> canonical s.p. energies'
          Write(iw,*) ' fermi energy=',alast(it)
          Write(iw,*) ' average cdelt=',del(it)
          Write(iw,'(a,a)')'  k0      ceqp        ce         v*v',  &
               '       u*v        cdel     overl      labels'
100       Format(//,' #canonical s.p. energies ',a,/,1x,33('-'),//)
       End If
       k0=0
       summ=zero; enjacek=zero
       Do ib=1,nb
          nd=id(ib); im=ia(ib)
          k1=ka(ib,it)+1; k2=ka(ib,it)+kd(ib,it)
          If(k1.Le.k2) Then
             Do k=1,nd
                k0=k0+1
                ! for Lipkin Nogami
                vvs=two*Sqrt(vk(k0,it))*Sqrt(one-vk(k0,it))    !2vu
                vvc=two*vk(k0,it)-one                          !2v^2-1
                summ(1)=summ(1)+vvs**2
                summ(2)=summ(2)+vvs**2*vvc
                summ(3)=summ(3)+vvs**4
                summ(4)=summ(4)+(vvs*vvc)**2
                ! search for main oscillator component
                ovmax=zero
                Do n=1,nd
                   s=Abs(ddc(n,k0,it))                         !canon orbitals in conf.space
                   If (s.Ge.ovmax) Then
                      ovmax=s; imax=n
                   End If
                End Do
                ! printing
                ek0=ek(k0,it)                                  !canon s.p. energies
                enjacek=enjacek+ek0*vk(k0,it)
                If(ek0.Lt.pwi) Then                            !print up to 'pwi'
                   vk0=vk(k0,it)                               !canon occupations v^2
                   If(vk0.Gt.-1.d-4) Then                       !print If signIficant v^2
                      dk0=-dk(k0,it)                           !canon s.p. deltas
                      ekk=Sqrt((ek0-ala(it))**2+dk(k0,it)**2)  !resulting cqpe
                      uuvv=Sqrt(Abs(vk0*(one-vk0)))            !resulting u*v
                      If(Print_Screen .And. do_print == 1) Then
                         iw=lfile
                         Write(iw,101) k0,ekk+ala(it),ek0,vk0,uuvv,dk0,ovmax,tb(im+imax)
101                      Format(i4,2f12.6,2(1x,f12.8),2(2x,f7.4),' ',a13)
                      End If
                   End If
                End If
             End Do !k0
          End If
       End Do !ib
       !--------------------------------------------
       ! Lipkin-Nogami
       !--------------------------------------------
       ssln(1,it)=summ(1)
       ssln(2,it)=summ(2)
       ssln(3,it)=summ(4)*summ(1)-summ(2)**2+summ(1)**3/4.0_pr-half*summ(3)*summ(1)
       If(Print_Screen .And. do_print == 1) Then
          iw=lfile
          Write(iw,*) ' Sum canonical e_v*V^2_k=',two*enjacek
       End If
    End Do !it
    !--------------------------------------------
    ! To thoout.dat, thodef.dat and hodef.dat
    !--------------------------------------------
    If(irecord.Ne.0) Then
       iappend=1
       Call expect(.True.)    !print  & record HFB+PAV results
       If(ierror_flag.Ne.0) Return
       iappend=0
    Else
       Call expect(.True.)    !print HFB+PAV results
       If(ierror_flag.Ne.0) Return
    End If
    !
  End Subroutine resu
  !=======================================================================
  !> Initializes all NAMELISTS
  !=======================================================================
  Subroutine initialize_HFBTHO_NAMELIST
    Implicit None
    ! HFBTHO_GENERAL
    number_of_shells    = 10
    oscillator_length   = -one
    basis_deformation   = zero
    proton_number       = 24
    neutron_number      = 26
    type_of_calculation = 1
    ! HFBTHO_INITIAL
    beta2_deformation = zero
    beta3_deformation = zero
    beta4_deformation = zero
    ! HFBTHO_ITERATIONS
    number_iterations = 100
    accuracy          = 1.D-5
    restart_file      = -1
    ! HFBTHO_FUNCTIONAL
    functional          = 'SLY4'
    add_initial_pairing = .False.
    type_of_coulomb     = 2
    include_3N_force    = .False.
    ! HFBTHO_PAIRING
    user_pairing    = .False.
    vpair_n         = -300.0_pr
    vpair_p         = -300.0_pr
    pairing_cutoff  =   60.0_pr
    pairing_feature =    0.5_pr
    ! HFBTHO_CONSTRAINTS
    lambda_values       = (/ 0, 0, 0, 0, 0, 0, 0, 0 /)
    lambda_active       = (/ 0, 0, 0, 0, 0, 0, 0, 0 /)
    expectation_values  = (/ 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr, 0.0_pr /)
    ! HFBTHO_BLOCKING
    proton_blocking  = (/ 0, 0, 0, 0, 0 /)
    neutron_blocking = (/ 0, 0, 0, 0, 0 /)
    ! HFBTHO_PROJECTION
    switch_to_THO    = 0
    projection_is_on = 0
    gauge_points     = 1
    delta_Z          = 0
    delta_N          = 0
    ! HFBTHO_TEMPERATURE
    set_temperature = .False.
    temperature     = zero
    ! HFBTHO_FEATURES
    collective_inertia     = .False.
    fission_fragments      = .False.
    pairing_regularization = .False.
    automatic_basis        = .False.
    ! HFBTHO_TDDFT
    filter                 = .False.
    fragment_properties    = .False.
    real_Z                 = 24.0_pr
    real_N                 = 26.0_pr
    ! HFBTHO_NECK
    set_neck_constrain     = .False.
    neck_value             = 0.5_pr
    ! HFBTHO_DEBUG
    number_Gauss        =  40
    number_Laguerre     =  40
    number_Legendre     =  80
    compatibility_HFODD = .False.
    number_states       = 500
    force_parity        = .True.
    print_time          = 0
    ! HFBTHO_RESTORATION
    PNP_is_on                   = 0
    number_of_gauge_points      = 1
    delta_neutrons              = 0
    delta_protons               = 0
    AMP_is_on                   = 0
    number_of_rotational_angles = 1
    maximal_angular_momentum    = 0
    !
  End Subroutine initialize_HFBTHO_NAMELIST
  !=======================================================================
  !> Read all NAMELISTS from hfbtho_NAMELIST.dat
  !=======================================================================
  Subroutine read_HFBTHO_NAMELIST(filename_nml)
    Implicit None
    Character(*), Intent(IN), Optional :: filename_nml
    ! TODO: Make this length a program-wide constant
    Character(256) :: filename
    Integer(ipr) :: ios,lnamelist=16
    If (present(filename_nml)) Then
        filename = filename_nml
    Else
        filename = 'hfbtho_NAMELIST.dat'
    End If
    Open(lnamelist,file=filename,DELIM='APOSTROPHE') ! 'QUOTE'
    !
    ierror_flag = 0
    !
    ! General input data
    Read(UNIT=lnamelist,NML=HFBTHO_GENERAL,iostat=ios)
    If (ios.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='Error in HFBTHO_GENERAL read'
       Return
    End If
    !
    ! Deformations of the initial WS solution
    Read(UNIT=lnamelist,NML=HFBTHO_INITIAL,iostat=ios)
    If (ios.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='Error in HFBTHO_INITIAL read'
       Return
    End If
    !
    ! Iterations
    Read(UNIT=lnamelist,NML=HFBTHO_ITERATIONS,iostat=ios)
    If (ios.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='Error in HFBTHO_ITERATIONS read'
       Return
    End If
    !
    ! Type of functional
    Read(UNIT=lnamelist,NML=HFBTHO_FUNCTIONAL,iostat=ios)
    If (ios.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='Error in HFBTHO_FUNCTIONAL read'
       Return
    End If
    !
    ! Characteristics of pairing
    Read(UNIT=lnamelist,NML=HFBTHO_PAIRING,iostat=ios)
    If (ios.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='Error in HFBTHO_PAIRING read'
       Return
    End If
    !
    ! Constraints
    Read(UNIT=lnamelist,NML=HFBTHO_CONSTRAINTS,iostat=ios)
    If (ios.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='Error in HFBTHO_CONSTRAINTS read'
       Return
    End If
    !
    ! Blocking
    Read(UNIT=lnamelist,NML=HFBTHO_BLOCKING,iostat=ios)
    If (ios.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='Error in HFBTHO_BLOCKING read'
       Return
    End If
    !
    ! Particle number projection
    Read(UNIT=lnamelist,NML=HFBTHO_PROJECTION,iostat=ios)
    If (ios.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='Error in HFBTHO_PROJECTION read'
       Return
    End If
    !
    ! Finite temperature
    Read(UNIT=lnamelist,NML=HFBTHO_TEMPERATURE,iostat=ios)
    If (ios.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='Error in HFBTHO_TEMPERATURE read'
       Return
    End If
    !
    ! Various features of the calculation
    Read(UNIT=lnamelist,NML=HFBTHO_FEATURES,iostat=ios)
    If (ios.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='Error in HFBTHO_FEATURES read'
       Return
    End If
    !
    ! Interface with TDDFT codes
    Read(UNIT=lnamelist,NML=HFBTHO_TDDFT,iostat=ios)
    If (ios.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='Error in HFBTHO_TDDFT read'
       Return
    End If
    !
    ! Constraint on the neck
    Read(UNIT=lnamelist,NML=HFBTHO_NECK,iostat=ios)
    If (ios.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='Error in HFBTHO_NECK read'
       Return
    End If
    !
    ! Debug
    Read(UNIT=lnamelist,NML=HFBTHO_DEBUG,iostat=ios)
    If (ios.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='Error in HFBTHO_DEBUG read'
       Return
    End If
    !
    ! Restoration of symmetries
    Read(UNIT=lnamelist,NML=HFBTHO_RESTORATION,iostat=ios)
    If (ios.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='Error in HFBTHO_RESTORATION read'
       Return
    End If
    !
    Close(lnamelist)
    !
  End Subroutine read_HFBTHO_NAMELIST
  !=======================================================================
  !> Performs some basis consistency checks of input data
  !=======================================================================
  Subroutine check_consistency
    Implicit None
    Integer(ipr) :: counter, i
    Real(pr) :: A, preset_inin(3)
    Character(30), Dimension(:) :: preset_forces(45)
    !
    If((abs(n00_INI).Lt.1).Or.(abs(n00_INI).GT.50)) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("number_of_shells = ",i6," out-of-bounds: [1,50]")') &
             n00_INI
       Return
    End If
    !
    If((npr_INI(1).Lt.1).Or.(npr_INI(2).Lt.1)) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("Z = ",i6," N = ",i6," out-of-bounds: (Z,N)>1")') &
             npr_INI(2),npr_INI(1)
       Return
    End If
    !
    If(Abs(kindhfb_INI).Ne.1) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("type_of_calculation = ",i6," unrecognized: (-1,1)")') &
             kindhfb_INI
       Return
    End If
    !
    If(epsi_INI.Lt.0) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("accuracy = ",e24.12," out-of-bounds: >0")') &
             epsi_INI
       Return
    End If
    !
    preset_inin( 1) = 1
    preset_inin( 2) = 2
    preset_inin( 3) = 3
    !
    counter=0
    Do i=1, 3
       If(Abs(inin_INI).Eq.preset_inin(i)) Then
          counter=1
          Exit
       End If
    End Do
    !
    If(counter.Eq.0) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("restart_file = ",i6," unrecognized: see list in publi")') &
             inin_INI
       Return
    End If
    !
    preset_forces( 1) = 'SIII'
    preset_forces( 2) = 'SKM*'
    preset_forces( 3) = 'SKP'
    preset_forces( 4) = 'SLY4'
    preset_forces( 5) = 'SLY5'
    preset_forces( 6) = 'SLY6'
    preset_forces( 7) = 'SLY7'
    preset_forces( 8) = 'SKI3'
    preset_forces( 9) = 'SKO'
    preset_forces(10) = 'SKX'
    preset_forces(11) = 'UNE0'
    preset_forces(12) = 'UNE1'
    preset_forces(13) = 'UNE2'
    preset_forces(14) = 'N0LO'
    preset_forces(15) = 'N1LO'
    preset_forces(16) = 'N2LO'
    preset_forces(17) = 'FITS'
    preset_forces(18) = 'D1'
    preset_forces(19) = 'D1p'
    preset_forces(20) = 'D1S'
    preset_forces(21) = 'D1N'
    preset_forces(22) = 'T0X0'
    preset_forces(23) = 'DME_LO'
    preset_forces(24) = 'DME_NLO'
    preset_forces(25) = 'DME_N2LO'
    preset_forces(26) = 'DME_NLOD'
    preset_forces(27) = 'DME_N2LOD'
    preset_forces(28) = 'REG_LO'
    preset_forces(29) = 'REG_NLO'
    preset_forces(30) = 'REG_N2LO'
    preset_forces(31) = 'REG_NLOD'
    preset_forces(32) = 'REG_N2LOD'
    preset_forces(33) = 'NEDF'
    preset_forces(34) = 'SeaLL1'
    preset_forces(35) = 'NEDF1'
    preset_forces(36) = 'NEDF2'
    preset_forces(37) = 'NEDF3'
    preset_forces(38) = 'NEDF4'
    preset_forces(39) = 'NEDF5'
    preset_forces(40) = 'HFB1'
    preset_forces(41) = 'SKM*mod'
    preset_forces(42) = 'SKOP'
    preset_forces(43) = 'SLY4mod'
    preset_forces(44) = 'READ'
    preset_forces(45) = 'RAND'
    !
    counter=0
    Do i=1, 45
       If(Trim(functional).Eq.Trim(preset_forces(i))) Then
          counter=1
          Exit
       End If
    End Do
    ! Functional must be in preset list
    If(counter.Eq.0) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("functional = ",a30," unrecognized: see list in publi")') &
             skyrme_INI
       Return
    End If
    ! Pairing cut-off must be positive
    If(pwi_INI.Lt.0) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("pairing_cutoff = ",i4," out-of-bounds: >=0")') &
             pwi_INI
       Return
    End If
    ! Pairing cut-off must be positive
    If(cpv1_INI.Lt.0.0.Or.cpv1_INI.Gt.1.0) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("pairing_feature = ",i4," out-of-bounds: [0.0,1.0]")') &
             cpv1_INI
       Return
    End If
    ! Options for Coulomb: -3, -2, -1, 0, 1, 2
    If(icou_INI.Lt.-4.Or.icou_INI.Gt.2) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("type_of_coulomb = ",i4," unrecognized: (-4,-3,-2,-1,0,1,2)")') &
             icou_INI
       Return
    End If
    ! Choices of basis (HO or THO): -1, 0, 1
    If(Abs(iLST_INI).Gt.1) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("switch_to_THO = ",i4," unrecognized: (-1,0,1)")') &
             iLST_INI
       Return
    End If
    ! At least one gauge point if projection is required
    If((keypjn_INI.Le.0.Or.keypjp_INI.Le.0).And.iproj_INI.Ne.0) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("gauge_points = ",2i4," out-of-bounds: >=0")') &
             keypjn_INI,keypjp_INI
       Return
    End If
    ! Number of protons must be greater than 0 for projection
    If((npr_INI(1)+npr1pj_INI).Lt.1.And.iproj_INI.Ne.0) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("delta_N = ",i4," out-of-bounds: N+dN>=1")') &
             npr1pj_INI
       Return
    End If
    ! Number of neutrons must be greater than 0 for projection
    If((npr_INI(2)+npr2pj_INI).Lt.1.And.iproj_INI.Ne.0) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("delta_Z = ",i4," out-of-bounds: Z+dZ>=1")') &
             npr2pj_INI
       Return
    End If
    ! Temperature must be positive
    If(temper.Lt.zero) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("temperature = ",i4," out-of-bounds: T>=0")') &
             temper
       Return
    End If
    ! Number of Gauss-Laguerre integration points between 0 and 100
    If(ngh_INI.Lt.1.Or.ngh_INI.Gt.100) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("number_Gauss = ",i4," out-of-bounds: [1,100]")') &
             ngh_INI
       Return
    End If
    ! Number of Gauss-Hermite integration points between 0 and 100
    If(ngl_INI.Lt.1.Or.ngl_INI.Gt.100) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("number_Laguerre = ",i4," out-of-bounds: [1,100]")') &
             ngl_INI
       Return
    End If
    ! Number of Gauss-Legendre integration points lower than 100
    If(nleg_INI.Gt.100) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("number_Legendre = ",i4," out-of-bounds: [-infty,100]")') &
             nleg_INI
       Return
    End If
    ! Number of Gauss-Legendre integration points between 1 and 100 for PNP
    If((nleg_INI.Lt.1.Or.nleg_INI.Gt.100).And.iproj_INI.Ne.0) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("number_Legendre = ",i4," out-of-bounds: [1,100]")') &
             nleg_INI
       Return
    End If
    ! Number of basis states must be greater than 0
    If(nstate_INI.Lt.1.And.basis_HFODD_INI) Then
       ierror_flag=ierror_flag+1
       Write(ierror_info(ierror_flag),'("number_states = ",i4," out-of-bounds: >0")') &
             nstate_INI
       Return
    End If
    !
  End Subroutine check_consistency
  !=======================================================================
  !> Initializes the solver based on user-defined input data. In particular
  !> defines the parameters of the functional
  !=======================================================================
  Subroutine initialize_HFBTHO_SOLVER
    !---------------------------------------------------------------------
    ! default parameters
    !---------------------------------------------------------------------
    Implicit None
    Real(pr) :: A
    Integer(ipr) :: it
    !------------------------------------
    ! tapes
    !------------------------------------
    lwin=41; lwou=42;  lwel=52; lres=57; lin=3
    !------------------------------------
    ! From Namelist or default values
    !------------------------------------
    nstate                = nstate_INI
    epsi                  = epsi_INI                  ! stop criteria
    Add_Pairing           = Add_Pairing_INI           ! add pairing starting from file
    icou                  = icou_INI                  ! Coulomb flag, see routine
    DO_FITT               = DO_FITT_INI               ! calculates quantities for reg.optimization
    debug_solver          = IDEBUG_INI                ! debug
    Parity                = Parity_INI                ! reflection symmetry
    Print_HFBTHO_Namelist = Print_HFBTHO_Namelist_INI ! Print Namelist
    !---------------------------------------------------------------------
    ! Pairing set by user
    !---------------------------------------------------------------------
    rho_c=0.160_pr
    If(set_pairing) Then
       CpV0(0)=V0n_INI
       CpV0(1)=V0p_INI
       CpV1(0)=cpv1_INI
       CpV1(1)=cpv1_INI
    End If
    pwi=pwi_INI
    !------------------------------------
    ! output control
    !------------------------------------
    If(n00_INI.Gt.0) Then
       Print_Screen = .True.
#if(DO_MASSTABLE==1 || DRIP_LINES==1 || DO_PES==1)
       if(mpi_size.eq.1) then
          lfile = lout+1 ! condensed output to screen & full thoout.dat file
       else
          lout = lfile   ! full output to thoout.dat file only
       endif
#else
       lfile = lout+1    ! condensed output to screen & full thoout.dat file
#endif
    Else
       Print_Screen=.False.
#if(DO_MASSTABLE==1 || DRIP_LINES==1 || DO_PES==1)
       lout = lfile ! condensed output to thoout.dat file only
#else
       lfile = lout-1 ! no output to screen & thoout.dat
#endif
    End If
    !------------------------------------
    ! Pi
    !------------------------------------
    PI=four*Atan(one)
    !------------------------------------
    ! blocking
    !------------------------------------
    blocking_never_done = .True.
    bloblo=0; blo123=0; blok1k2=0;  keyblo=0
    blomax=0; nkblo=0;  iparenti=0
    blocanon=0;         eqpmin=zero
    blocking_mode = 0 ! 0 = Off, 1 = Specific, -1 = Auto
    Do it=1,2
       If(nkblo_INI(it,1).Ne.0) Then
          If(nkblo_INI(it,2).Ne.0) Then
             blocking_mode(it) = 1
          Else
             blocking_mode(it) = -1
          End If
       End If
    End Do
    !------------------------------------
    ! buffers
    !------------------------------------
    eres=zero;  eresu=zero;  eresl=zero;
    eresj=zero; eresbl=zero; ereslbl=' 00[00,00,00]'
    !------------------------------------
    ! def parameters
    !------------------------------------
    ffdef3=Sqrt(five/(four*pi))/two
    ffdef4=Sqrt(117.0_pr)/(four*pi)
    ffdef5=Sqrt(nine/(four*pi))/eight
    ffdef6=Sqrt(five*pi)/three
    ffdef7=Sqrt(pi)/four
    !------------------------------------
    ! former linear mixing
    !------------------------------------
    xmix0 = 0.1_pr ! lowest mixing parameter  (redefined later)
    xmix = 0.3_pr  ! initial mixing parameter (changes every iteration)
    xmax = 0.9_pr  ! Mario: 1.0, Nicolas: 0.9
    bbroyden = 'L'
    !------------------------------------
    ! misc (redefined later)
    !------------------------------------
    rehfbcan=0.0_pr; depnp=0.0_pr; ala2=0.00_pr
    ept=-2.0_pr; del=1.0_pr; ala=-7.0_pr; frept=0.0_pr; coept=0.0_pr
    ala1(1)=-14.6851_pr; ala1(2)=-3.7522_pr; si=1.0_pr
    iqrpa=0; icacou=0;  icacoupj=0; icahartree=0; iasswrong=0
    iError_in_HO=0;  iError_in_THO=0
    ECMHFB=0.0_pr; ECMPAV=0.0_pr
    If(use_full_cm_cor) Then
       tz(1:2) = npr_INI(1:2)
       A = tz(1) + tz(2)
       facECM = A/(A-1.0_pr)
    End If
    entropy(:)=zero
    !------------------------------------
    ! Saxon-Woods: von koepf und ring, z.phys. (1991)
    !------------------------------------
    v0ws=-71.28_pr; akv=0.4616_pr; r0v=1.2334_pr; av=0.6150_pr
    vso=11.1175_pr; rso=1.1443_pr; aso=0.6476_pr
    !------------------------------------
    ! fixed text
    !------------------------------------
    tp(1)='+'; tp(2)='-'; tis(1)='n'; tis(2)='p';
    tit(1)='neutrons'; tit(2)='protons '
    tl(0)='s'; tl(1)='p'; tl(2)='d'; tl(3)='f'; tl(4)='g'
    tl(5)='h'; tl(6)='i'; tl(7)='j'; tl(8)='k'; tl(9)='l'
    tl(10)='m'; tl(11)='n'; tl(12)='o'; tl(13)='p'; tl(14)='q'
    tl(15)='r'; tl(16)='s'; tl(17)='t'; tl(18)='u'; tl(19)='v'; tl(20)='w'
    !------------------------------------
    ! fixed parity sign
    !------------------------------------
    tpar(1)=+1; tpar(2)=-1;
    !------------------------------------
    ! physical constants
    !------------------------------------
    amn=938.90590_pr
    amu=931.4940130_pr; r0=1.20_pr
    alphi=137.036020_pr; hqc=197.328910_pr
    !------------------------------------
    ! Coulomb
    !------------------------------------
    !chargee2=hqc/alphi
    !chargee2=1.43997841_pr
    !-----------------------------------
    ! set the loops over particle types
    !-----------------------------------
    itmin=1 ; itmax = 2;
    If(npr_INI(1).Eq.0) itmin = 2
    If(npr_INI(2).Eq.0) itmax = 1
    !-----------------------------------
    ! error flag and info
    !-----------------------------------
    ierror_flag=0
    ierror_info(ierror_flag)='No errors in the solver!'
    !
    Call set_functional_parameters(skyrme_INI,.False.)
    !-----------------------------------
    ! set multipole moments units
    !-----------------------------------
    Call moments_setUnits
    !
  End Subroutine initialize_HFBTHO_SOLVER
  !=======================================================================
  !> Adjust deformation beta2 and oscillator frequency of the basis and
  !> deformations beta2,beta4 of the initial Woods-Saxon potential
  !=======================================================================
  Subroutine adjust_basis(q2val,initial_deformations,zz,nn)
    Implicit None
    Logical, Intent(IN) :: initial_deformations
    Integer, Intent(IN) :: zz,nn
    Real(pr), Intent(IN) :: q2val
    Real(pr) :: OMEGA0,hbarc,mass_neut,SCALEM
    ! Loose fit based on fission of 240Pu, see PRC 90, 054305 (2014)
    If(q2val.Gt.30.0_pr) Then
       OMEGA0=0.1_pr*EXP(-0.02_pr*q2val)*q2val + 6.5_pr
    Else
       OMEGA0=8.1464_pr
    End If
    hbarc=197.32891_pr; mass_neut=938.9059_pr
    ! Rescaling so that the fit makes sense for mass A, not just 240Pu
    oscillator_length=hbarc/Sqrt(OMEGA0*mass_neut)
    basis_deformation=0.05_pr*Sqrt(abs(q2val))
    If(initial_deformations) Then
       beta2_deformation=basis_deformation
       beta4_deformation=beta2_deformation/100.0_pr
    End If
  End Subroutine adjust_basis
  !=======================================================================
  !> Counts and orders basis states in cylindrical coordinates
  !=======================================================================
  Subroutine base0(lpr)
    Implicit None
    Logical :: lpr
    Integer(ipr) :: iw,k,nre,nze,ke,la,le,ip,ir,iz,il,is,Iall,ilauf,jlauf,ib,nd
    Integer(ipr) :: NOSCIL
    Real(pr), Allocatable :: e(:)
    Real(pr) :: hbz,hbp,ee
    !
    If(n00.Gt.n00max) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='STOP: too large n00 versus n00max'
       Return
    End If
    !-----------------------------------------------
    ! MAXIMUM NUMBER OF THE HO SHELLS (n00,NOSCIL)
    ! (7,120),(8,165),(9,220),(10,286),(11,364)
    ! (12,455),(14,680),(16,969),(18,1330),(20,1771)
    !-----------------------------------------------
    NOSCIL=(n00+1)*(n00+2)*(n00+3)/6
    !-----------------------------------------------
    ! count all states for n00max
    !-----------------------------------------------
    nze=n00max; nre=n00max/2; ke=n00max
    If(basis_HFODD) Then
       nze=n00; nre=n00/2; ke=n00
    End If
    Iall=0;
    Do k=1,ke+1
       la=k-1; le=min0(ke,k)
       Do ip=1,2
          Do ir=0,nre
             Do iz=0,nze
                Do il=la,le
                   Do is=+1,-1,-2
                      If (iz+2*ir+il.Gt.n00max)    Cycle
                      If (il+(is+1)/2.Ne.k)        Cycle
                      If (Mod(iz+il,2).Ne.ip-1)    Cycle
                      Iall=Iall+1
                   End Do
                End Do
             End Do
          End Do
       End Do
    End Do
    !-----------------------------------------------
    ! charge all energies for n00max
    !-----------------------------------------------
    Allocate(e(Iall))
    hbz=two*hbzero/bz**2; hbp=two*hbzero/bp**2;
    Iall=0;
    Do k=1,ke+1
       la=k-1; le=min0(ke,k)
       Do ip=1,2
          Do ir=0,nre
             Do iz=0,nze
                Do il=la,le
                   Do is=+1,-1,-2
                      If (iz+2*ir+il.Gt.n00max)    Cycle
                      If (il+(is+1)/2.Ne.k)        Cycle
                      If (Mod(iz+il,2).Ne.ip-1)    Cycle
                      Iall=Iall+1
                      e(Iall)=hbz*(Real(iz,Kind=pr)+half) &
                         +hbp*(two*Real(ir,Kind=pr)+Real(il,Kind=pr)+one)
                   End Do
                End Do
             End Do
          End Do
       End Do
    End Do
    !-----------------------------------------------
    ! sort energies and derive base cut-off energy
    !-----------------------------------------------
    Call ord(Iall,e);
    If(Iall.Gt.NOSCIL) Then
       EBASECUT=E(NOSCIL)+1.0D-5
    Else
       EBASECUT=E(Iall)+1.0D-5
    End If
    If(basis_HFODD) EBASECUT=E(nstate)+1.0D-5
    Deallocate(e)
    !-----------------------------------------------
    ! calculate the actual states
    !-----------------------------------------------
    nze=n00max; nre=n00max/2; ke=n00max
    If(basis_HFODD) Then
       nze=n00; nre=n00/2; ke=n00
    EndIf
    ib=0; ilauf=0; ndx=0; nzx=0; nrx=0; nlx=0; nqp=0; nuv=0; nnx=0; nrlx=0; nox=0
    ! loop over k-quantum number
    Do k=1,ke+1
       la=k-1; le=min0(ke,k)
       ! loop over parity
       If(.Not.Parity) jlauf=ilauf !Nop
       Do ip=1,2
          If(Parity) jlauf=ilauf !Yesp
          Do ir=0,nre
             Do iz=0,nze
                Do il=la,le
                   Do is=+1,-1,-2
                      If (iz+2*ir+il.Gt.n00max)    Cycle
                      If (il+(is+1)/2.Ne.k)        Cycle
                      If (Mod(iz+il,2).Ne.(ip-1))  Cycle
                      ee=hbz*(Real(iz,Kind=pr)+half)&
                    +hbp*(two*Real(ir,Kind=pr)+Real(il,Kind=pr)+one)
                      If(ee.Lt.EBASECUT) Then
                         ilauf=ilauf+1
                         nzx=Max(nzx,iz); nrx=Max(nrx,ir); nlx=Max(nlx,il)
                         nnx=Max(nnx,iz+2*ir+il); nrlx=Max(nrlx,2*ir+il)
                         nox=max(nox,k-1)
                      End If
                   End Do
                End Do
             End Do
          End Do
          If(Parity) Then                !Yesp
             If (ilauf.Gt.jlauf) Then
                ib=ib+1
                nd=ilauf-jlauf
                ndx=Max(ndx,nd)
                nqp=nqp+nd; nuv=nuv+nd*nd
             End If
          End If
       End Do
       If(.Not.Parity) Then              !Nop
          If(ilauf.Gt.jlauf) Then
             ib=ib+1
             nd=ilauf-jlauf
             ndx=Max(ndx,nd)
             nqp=nqp+nd; nuv=nuv+nd*nd
          End If
       End If
    End Do
    nbx=ib; ntx=ilauf
    !-----------------------------------------------
    ! print statistics
    !-----------------------------------------------
    If(lpr .And. do_print == 1) Then
       Do iw=lout,lfile
          Write(iw,*)
          Write(iw,'(a)')  '  ---------------------------------------'
          Write(iw,'(a)')  '        Harmonic Oscillator Basis        '
          Write(iw,'(a)')  '  ---------------------------------------'
          Write(iw,'(a,2(i6,2x),a)') '  NUV, NQP:                      ',nuv,nqp
          Write(iw,'(a,2(i6,2x),a)') '  Comparison with bookkeeping spherical basis:'
          Write(iw,'(a,2(i6,2x),a)') '  n00:                           ',n00,n00,  &
               'Maximal number of shells'
          Write(iw,'(a,2(i6,2x),a)') '  nbx, 2*n00+1:                  ',nbx,2*n00+1,  &
               'Maximal number of K-blocks'
          Write(iw,'(a,2(i6,2x),a)') '  ntx, (n00+1)*(n00+2)*(n00+3)/6 ',ntx,(n00+1)*(n00+2)*(n00+3)/6,  &
               'Max.num. p/n levels'
          Write(iw,'(a,2(i6,2x),a)') '  nzx, n00:                      ',nzx,n00,  &
               'Maximal nz-quantum number'
          Write(iw,'(a,2(i6,2x),a)') '  nrx, n00/2  :                  ',nrx,n00/2,  &
               'Maximal nr-quantum number'
          Write(iw,'(a,2(i6,2x),a)') '  nlx, n00:                      ',nlx,n00,  &
               'Maximal ml-quantum number'
          Write(iw,'(a,2(i6,2x),a)') '  ndx, (n00+2)*(n00+2)/4:        ',ndx,(n00+2)*(n00+2)/4,  &
               'Maximal dim. of one k-block'
          Write(iw,*)
       End Do
    End If
    !
  End Subroutine base0
  !=======================================================================
  !> Define quantum numbers for basis states and computes matrix sizes
  !=======================================================================
  Subroutine base(lpr)
    Implicit None
    Logical :: lpr
    Integer(ipr) :: nze,nre,ke,ib,ilauf,jlauf,nom,nnm,  &
         k,la,le,ip,ir,iz,il,is,nn,ND,IBX,N1,N2,iw, &
         klauf
    Real(pr) :: hbz,hbp,ee
    !
    hbz=two*hbzero/bz**2; hbp=two*hbzero/bp**2;
    !
    nze=n00max; nre=n00max/2; ke=n00max
    If(basis_HFODD) Then
       nze=n00; nre=n00/2; ke=n00
    End If
    ib=0; ilauf=0; nzm=0; nrm=0; nlm=0; nom=0; nnm=0; jlauf = 0
    If(finite_range.or.coulomb_gaussian) Then
       If(Allocated(ib_zrls)) Deallocate(i_zrls,ib_zrls)
       If(.not.allocated(ib_zrls)) Then
          Allocate( i_zrls(0:nze,0:nre,0:n00max,0:1))
          Allocate(ib_zrls(0:nze,0:nre,0:n00max,0:1))
       End If
       i_zrls = 0
       ib_zrls = 0
       !-----------------------------------------------
       ! loop over k-quantum number with z=0
       !-----------------------------------------------
       Do k=1,n00max+1
          la=k-1; le=min0(n00max,k)
          klauf = ilauf
          Do ir=0,nre
             Do il=la,le
                Do is=+1,-1,-2
                   jlauf = 0
                   do ip = 1,2
                      do iz = 0,nze
                         If (iz+2*ir+il.Gt.n00max) Cycle
                         If (il+(is+1)/2.Ne.k)     Cycle
                         If (Mod(iz+il,2).Ne.ip-1) Cycle
                         ee=hbz*(Real(iz,Kind=pr)+half) &
                           +hbp*(two*Real(ir,Kind=pr)+Real(il,Kind=pr)+one)
                         If(ee.Lt.EBASECUT) Then
                            jlauf=jlauf+1
                            nzz(ilauf+1,jlauf) = iz
                         End If
                      enddo
                   enddo
                   if(jlauf.ne.0) then
                      ilauf = ilauf + 1
                      nrr(ilauf)=ir; nll(ilauf)=il; nss(ilauf)=is;
                      noo(ilauf) = k; nzzx(ilauf) = jlauf
                   endif
                End Do
             End Do
          End Do
       End Do ! end k
       nttx=ilauf
    End If
    !
    nze=n00max; nre=n00max/2; ke=n00max
    If(basis_HFODD) Then
       nze=n00; nre=n00/2; ke=n00
    End If
    ib=0; ilauf=0; nzm=0; nrm=0; nlm=0; nom=0; nnm=0
    !-----------------------------------------------
    ! loop over k-quantum number
    !-----------------------------------------------
    Do k=1,ke+1
       la=k-1; le=min0(ke,k)
       ! loop over parity
       If(.Not.Parity) jlauf=ilauf !Nop
       Do ip=1,2
          If(Parity) jlauf=ilauf   !Yesp
          Do ir=0,nre
             Do iz=0,nze
                Do il=la,le
                   Do is=+1,-1,-2
                      If (iz+2*ir+il.Gt.n00max) Cycle
                      If (il+(is+1)/2.Ne.k)     Cycle
                      If (Mod(iz+il,2).Ne.ip-1) Cycle
                      ee=hbz*(Real(iz,Kind=pr)+half)&
                    +hbp*(two*Real(ir,Kind=pr)+Real(il,Kind=pr)+one)
                      If(ee.Lt.EBASECUT) Then
                         ilauf=ilauf+1
                         If (ilauf.Gt.ntx) Then
                            ierror_flag=ierror_flag+1
                            ierror_info(ierror_flag)='STOP: in base: ntx too small'
                            Return
                         End If
                         nz(ilauf)=iz; nr(ilauf)=ir; nl(ilauf)=il; ns(ilauf)=is; npar(ilauf)=ip
                         If(finite_range.or.coulomb_gaussian) Then
                            i_zrls(iz,ir,il,(is+1)/2) = ilauf
                            ib_zrls(iz,ir,il,(is+1)/2) = ib+1
                         End If
                         nn =iz+2*ir+il
                         Write(tb(ilauf),100) 2*k-1,tp(ip),nn,iz,il
100                      Format(i2,a1,'[',i2,',',i2,',',i2,']')
                         Do iw=lout,lfile
                            If(lpr.And.debug_solver.Ge.3 .And. do_print == 1) &
                               Write(iw,'(i4,a,i2,a,i2,a,i2,a,i2,a,i2,a,2x,a,1x,a,f14.8)')  &
                                 ilauf,'   nn=',nn,'   nz=',iz,'   nr=',ir,  &
                                 '   ml=',il,'  ms=',is,' /2',tb(ilauf),'e=',ee
                         End Do
                         nzm=Max(nzm,iz); nrm=Max(nrm,ir); nlm=Max(nlm,il)
                         nom=Max(nom,2*k-1); nnm=Max(nnm,iz+2*ir+il)
                      End If
                   End Do
                End Do
             End Do
          End Do
          !-----------------------------------------------
          ! Block memory
          !-----------------------------------------------
          If(Parity) Then                !Yesp
             If (ilauf.Gt.jlauf) Then
                ib=ib+1
                ia(ib)=jlauf; id(ib)=ilauf-jlauf
                ikb(ib)=k; ipb(ib)=ip
                Write(txb(ib),'(i3,a,i2,a,a1)') ib,'. block:  k=',k+k-1,'/2',tp(ip)
                !ir=(ib+1)/2
                !Write(*,*)  ib,2*k-1,'2*Omega=',2*ir - 1
                Do iw=lout,lfile
                   If(lpr.And.debug_solver.Ge.3 .And. do_print == 1) &
                      Write(iw,'(/,a,i3,a,a1)')'  For the above block:  k=',k+k-1,'/2',tp(ip)
                End Do
             End If
             If(id(ib).Eq.0) Then
                ierror_flag=ierror_flag+1
                ierror_info(ierror_flag)='STOP: in base Block Memory(1)'
                Return
             End If
          End If
       End Do ! end of ip
       !-----------------------------------------------
       ! Block memory
       !-----------------------------------------------
       If(.Not.Parity) Then               !Nop
          If (ilauf.Gt.jlauf) Then
             ib=ib+1
             ia(ib)=jlauf; id(ib)=ilauf-jlauf
             nn = nz(ilauf)+2*nr(ilauf)+nl(ilauf); ip = 2 - Mod(nn,2)
             ikb(ib)=k; ipb(ib)=ip
             Write(txb(ib),'(i3,a,i2,a,a1)') ib,'. block:  k=',k+k-1,'/2',tp(ip)
             Do iw=lout,lfile
                If(lpr.And.debug_solver.Ge.3 .And. do_print == 1) &
                   Write(iw,'(/,a,i3,a,a1)')'  For the above block:  k=',k+k-1,'/2',tp(ip)
             End Do
          End If
          If(id(ib).Eq.0) Then
             ierror_flag=ierror_flag+1
             ierror_info(ierror_flag)='STOP: in base Block Memory(2)'
             Return
          End If
       End If
    End Do ! end k
    nb=ib;  nt=ilauf
    !-----------------------------------------------
    ! broyden/linear mixing (storage)
    !-----------------------------------------------
    nhhdim=0
    Do ib=1,NB
       ND=ID(ib)
       Do N1=1,ND
          Do N2=1,N1
             nhhdim=nhhdim+1
          End Do
       End Do
    End Do
    nhhdim2=2*nhhdim; nhhdim3=3*nhhdim; nhhdim4=4*nhhdim
    If(Allocated(brin)) Deallocate(brin,brout)
    If(neck_constraints) Then
       Allocate(brin(nhhdim4+lambdaMax+1),brout(nhhdim4+lambdaMax+1))
    Else
       Allocate(brin(nhhdim4+lambdaMax),brout(nhhdim4+lambdaMax))
    End If
    brin=zero; brout=zero
    !-----------------------------------------------
    ! Print statistics
    !-----------------------------------------------
    If(lpr .And. do_print == 1) Then
       Do iw=lout,lfile
          Write(iw,'(a,i4)')   '  Actual basis used'
          Write(iw,'(a,i4)')   '  Number of blocks: nb .......: ',nb
          Write(iw,'(a,i4)')   '  Number of levels: nt .......: ',nt
          Write(iw,'(a,i4)')   '  Maximal 2*omega : nom ......: ',nom
          Write(iw,'(a,i4)')   '  Maximal nz:       nzm ......: ',nzm
          Write(iw,'(a,i4)')   '  Maximal nr:       nrm ......: ',nrm
          Write(iw,'(a,i4)')   '  Maximal ml:       nlm ......: ',nlm
          Write(iw,'(a,i4)')   '  Maximal N=nz+2*nr+nl .......: ',nnm
          Write(iw,'(a,i4)')   '  2 x biggest block dim. .....: ',ndx2
          Write(iw,'(a,i8)')   '  Non-zero elements of h .....: ',nhhdim
          Write(iw,'(a,i8)')   '  Number of Broyden elements .: ',nhhdim4
          Write(iw,'(a,i4)')
       End Do
    End If
    If(nzm.Ge.n00max.Or.(nom-1)/2.Eq.n00max) Then
       Write(*,*) 'nzm=',nzm,'  (nom-1)/2=',(nom-1)/2,'  n00max=',n00max
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='STOP: Please increase n00max to have correct basis'
    End If
  End Subroutine base
  !=======================================================================
  !> Subroutine gaupol() computes the harmonic oscillator wave functions
  !> and their first derivatives in cylindrical coordinates \f$ (\rho,\theta,z )\f$. The full HO
  !> wave function reads
  !>  \f[
  !>     \psi_{n_{r}\Lambda n_{z}}(\rho, \theta, z) =
  !>     \psi_{n_{r}}^{|\Lambda|}(\eta) \psi_{n_{z}}(\xi) \frac{e^{i\Lambda\theta}}{\sqrt{2\pi}}
  !>                  \chi_{\Sigma}(\sigma),
  !>  \f]
  !>  where the dimensionless coordinates are \f$ \eta = \beta_{\perp}^{2}\rho^2 \f$,
  !>  \f$ \xi = \beta_{z}z \f$, \f$ \Lambda \f$ is the eigenvalue of \f$ \hat{\ell}_{z} \f$,
  !>  \f$ \Sigma \f$ the eigenvalue of \f$ \hat{s}_{z} \f$ and \f$ \sigma \f$ is a number
  !>  that labels the spin eigenstate in coordinate space, by opposition to \f$ \chi_{\Sigma} \f$,
  !>  which refers to a vector in spin space. Under the action of time-reversal, the HO
  !>  basis states turn into \f$ \psi_{n_{r}-\Lambda n_{z}}(\rho, \theta, z) \f$, where
  !>  the spatial part proportional to \f$ \rho, z \f$ is unchanged.
  !=======================================================================
  Subroutine gaupol(lpr)
    Implicit None
    Logical :: lpr
    Real(pr) :: w0,z,x,s,s0,s1,w00,w4pii,dsq,d1,d2,d3,d4,hs0,hs1
    Integer(ipr) :: ih,il,iw,ix,n,l,n1,n2
    !-----------------------------------------------
    !> - For the z-coordinate, the HO wave function reads
    !>     \f[
    !>        \psi_{n_{z}}(z) = \mathcal{N}_{n_{z}}\beta_{z}^{1/2} e^{-\xi^{2}/2} H_{n_{z}}(\xi),
    !>     \f]
    !>   with the normalization factor
    !>     \f[
    !>        \mathcal{N}_{n_{z}}
    !>      = \left( \frac{1}{\sqrt{\pi}2^{n_{z}}n_{z}!} \right)^{1/2}.
    !>     \f]
    !>   The first derivative of the HO wave function is
    !>     \f[
    !>        \frac{\partial\psi_{n_{z}}}{\partial z} = \beta_{z}
    !>            \left[ \sqrt{2n_{z}} \psi_{n_{z}-1}(\xi) - \xi \psi_{n_{z}}(\xi) \right].
    !>     \f]
    !>   In the code, variable \f$ \mathtt{qh} \f$ contains \f$ (1/\beta^{1/2}_{z})\psi_{n_{z}} \f$
    !>   on the nodes of the Gauss-Hermite integration mesh; variable \f$ \mathtt{qh1} \f$
    !>   contains the value of \f$ (1/\beta^{3/2}_{z})\partial\psi_{n_{z}}/\partial z \f$ on that same
    !>   mesh. Note that this is consistent with a dimension analysis: \f$ \psi_{n_{z}}(z) \f$
    !>   has dimension \f$ [L]^{-1/2} \f$ so that the integral over space of the square
    !>   \f$ |\psi_{n_{z}}(z)|^2 \f$ is dimensionless and can be interpreted as a probability.
    !>   Similarly, the dimension of the derivative \f$ \partial\psi_{n_{z}}/\partial z \f$ is
    !>   \f$ [L]^{-3/2} \f$. In the code, objects \f$ \mathtt{qh} \f$ and \f$ \mathtt{qh1} \f$
    !>   are thus dimensionless.
    !-----------------------------------------------
    w4pii=pi**(-0.250_pr)
    Do ih=1,ngh
       z=xh(ih); w0=w4pii*Exp(-half*z*z)
       w0 =w0*Sqrt(wh(ih))
       qh(0,ih)=w0;       qh(1,ih)=sq(2)*w0*z
       qh1(0,ih)=-w0*z;   qh1(1,ih)=sq(2)*w0*(one-z*z)
       Do n=2,nzm
          qh(n,ih)=sqi(n)*(sq(2)*z*qh(n-1,ih)-sq(n-1)*qh(n-2,ih))
          qh1(n,ih)=sq(n+n)*qh(n-1,ih)-z*qh(n,ih)
       End Do
    End Do
    !-----------------------------------------------
    !> - For the \f$\rho\f$-coordinate, the HO wave function reads
    !>     \f[
    !>        \psi_{n_{r}}^{\Lambda}(\eta) =
    !>        \mathcal{N}_{n_{r}}^{\Lambda} \beta_{\perp}\sqrt{2}
    !>              \eta^{|\Lambda|/2} e^{-\eta/2} L_{n_{r}}^{|\Lambda|}(\eta),
    !>     \f]
    !>   with the normalization factor
    !>     \f[
    !>        \mathcal{N}_{n_{r}}^{\Lambda}
    !>      = \left( \frac{n_{r}}{ (n_{r} + |\Lambda|)!} \right)^{1/2}.
    !>     \f]
    !>   The first derivative of the HO wave function is
    !>     \f[
    !>        \frac{\partial\psi_{n_{r}}^{\Lambda}}{\partial\rho} = \frac{\beta_{\perp}}{\sqrt{\eta}}
    !>            \left[ (2n_{r} + |\Lambda| - \eta)\psi_{n_{r}}^{\Lambda}(\eta)
    !>                   -2\sqrt{n_{r}(n_{r}+|\Lambda|)}\psi_{n_{r}-1}^{\Lambda}(\eta) \right].
    !>     \f]
    !>   In the code, variable \f$ \mathtt{ql} \f$ contains the value of
    !>   \f$ (1/\beta_{\perp}\sqrt{2})\psi_{n_{r}}^{\Lambda} \f$ on the nodes of the Gauss-Laguerre
    !>   integration mesh; variable \f$ \mathtt{ql1} \f$ contains the value of
    !>   \f$ (\sqrt{\eta}/\beta_{\perp}^{2}\sqrt{2}) \partial\psi_{n_{r}}^{\Lambda}/\partial \rho \f$
    !>   on that same mesh.
    !-----------------------------------------------
    Do il=1,ngl
       x=xl(il); w00=sq(2)*Exp(-half*x)
       Do l=0,nlm
          w0=w00*Sqrt(half*wl(il)*x**l)
          ql(0,l,il)=wfi(l)*w0;         ql(1,l,il)=(l+1-x)*wfi(l+1)*w0
          ql1(0,l,il)=(l-x)*wfi(l)*w0;  ql1(1,l,il)=(Real(l*l+l,Kind=pr) &
                                                  -x*Real(l+l+3,Kind=pr)+x*x)*wfi(l+1)*w0
          Do n=2,nrm
             dsq=sq(n)*sq(n+l); d1=Real(n+n+l-1,Kind=pr)-x
             d2=sq(n-1)*sq(n-1+l); d3=n+n+l-x; d4=two*dsq
             ql(n,l,il)=(d1*ql(n-1,l,il)-d2*ql(n-2,l,il))/dsq
             ql1(n,l,il)=d3*ql(n,l,il)-d4*ql(n-1,l,il)
          End Do
       End Do
    End Do
    !-----------------------------------------------
    ! Test accuracy for Hermite orthonormalization
    !-----------------------------------------------
    hs0=zero; hs1=two
    Do n1=0,nzm
       Do n2=0,n1
          If (Mod(n1-n2,2).Eq.0) Then
             s=zero
             Do ih=1,ngh
                s=s+qh(n1,ih)*qh(n2,ih)
             End Do
             If(n1.Ne.n2) Then
                hs0=Max(s,hs0)
             Else
                hs1=Min(s,hs1)
             End If
          End If
       End Do
    End Do
    !-----------------------------------------------
    ! Test accuracy for Laguerre orthonormalization
    !-----------------------------------------------
    s0=zero; s1=two
    Do l=0,nlm
       Do n1=0,nrm
          Do n2=0,n1
             s=zero
             Do il=1,ngl
                s=s+ql(n1,l,il)*ql(n2,l,il)
             End Do
             If(n1.Ne.n2) Then
                s0=Max(s,s0)
             Else
                s1=Min(s,s1)
             End If
          End Do
       End Do
    End Do
    !-----------------------------------------------
    ! print accuracy
    !-----------------------------------------------
    If(lpr .And. do_print == 1) Then
       Do iw=lout,lfile
          Write(iw,'(a)')  '  ---------------------------------------'
          Write(iw,'(a)')  '            Integration Meshes           '
          Write(iw,'(a)')  '  ---------------------------------------'
          Write(iw,'(a,i3)')  &
               '  Number of Gauss-Hermite mesh points ngh ....: ',ngh
          Write(iw,'(a,i3)')  &
               '  Number of Gauss-Laguerre mesh points ngl ...: ',ngl
          Write(iw,'(a,i3)')  &
               '  Number of Gauss-Legendre mesh points nleg ..: ',nleg
          Write(iw,'(a)') &
               '  Integration boundaries'
          Write(iw,'(2(a,f12.8))')  &
               '    Hermite  - from xh(1)  =',xh(1),  ' to xh(ngh)   =',xh(ngh)
          Write(iw,'(2(a,f12.8))')  &
               '    Laguerre - From xl(1)  =',xl(1),  ' to xl(ngl)   =',xl(ngl)
          If(nleg.Gt.0) Then
             Write(iw,'(2(a,f12.8))')  &
                  '    Legendre - From xleg(1)=',xleg(1),' to xleg(nleg)=',xleg(nleg)
          End If
          Write(iw,*)  &
               ' Max.dev.in:     Orthogonality            Normalization'
          Write(iw,*) ' Hermite  ',hs0,Abs(one-hs1)
          Write(iw,*) ' Laguerre ',s0,Abs(one-s1)
       End Do
    End If
    !-----------------------------------------------
    ! debug
    !-----------------------------------------------
    If (lpr.And.debug_solver.Ge.4 .And. do_print == 1) Then
       ix=3
       Do iw=lout,lfile
          Write(iw,*) ' nz    qh(nz,ih=1,...)'
          Do n=0,nzm
             Write(iw,'(i4,3f15.8)') n,(qh(n,ih),ih=1,ix)
             Write(iw,'(i4,3f15.8)') n,(qh1(n,ih),ih=1,ix)
             Write(iw,*) ' '
          End Do
          Do l=0,nlm
             Write(iw,*) ' nr ml    ql(nr,l,il=1,...)'
             Do n=0,nrm
                Write(iw,'(i4,i3,3f15.8)') n,l,(ql(n,l,il),il=1,ix)
                Write(iw,'(i4,i3,3f15.8)') n,l,(ql1(n,l,il),il=1,ix)
                Write(iw,*) ' '
             End Do
          End Do
       End Do
       !-----------------------------------------------
       ! Test for Hermite polynomials normalization
       !-----------------------------------------------
       Do n1=0,nzm
          Do n2=0,n1
             If (Mod(n1-n2,2).Eq.0) Then
                s=zero
                Do ih=1,ngh
                   s=s+qh(n1,ih)*qh(n2,ih)
                End Do
                Do iw=lout,lfile
                   Write(iw,100) n1,n2,s
                End Do
100             Format(' Gauss-Hermite: n1=',i3,'  n2=',i3,f20.8)
             End If
          End Do
       End Do
       !-----------------------------------------------
       ! Test for Laguerre polynomials normalization
       !-----------------------------------------------
       Do l=0,nlm
          Do n1=0,nrm
             Do n2=0,n1
                s=zero
                Do il=1,ngl
                   s=s+ql(n1,l,il)*ql(n2,l,il)
                End Do
                Do iw=lout,lfile
                   Write(iw,101) l,n1,n2,s
101                Format(' Gauss Laguerre: l='  &
                        ,i2,' n1=',i3,'  n2=',i3,f20.8)
                End Do
             End Do
          End Do
       End Do
    End If
    !
    Call coordinateLST(.False.)  ! coordinate LST
    !
  End Subroutine gaupol
  !=======================================================================
  !>
  !> Scratch initialization of the self-consistent loop by diagonalization
  !> of the deformed Woods-Saxon Hamiltonian. We assume the nuclear surface
  !> is parametrized as follows
  !>
  !>   \f[
  !>      R(\theta,\varphi) = R_0\, c(\beta) \left[ 1 + \sum_{\lambda=1}^{\lambda_{\rm max}}
  !>                         \beta_{\lambda}Y_{\lambda 0}(\theta,\varphi) \right]
  !>                        = R_0\, c(\beta) f(\theta,\varphi)
  !>   \f]
  !>
  !>  with \f$ R_0 = 1.2 A^{1/3} \f$, \f$ c(\beta) \f$ the volume conservation
  !>  coefficient, \f$ \beta_{\lambda} \f$ the deformation parameters and
  !>  \f$ Y_{\lambda\mu}(\theta,\varphi) \f$ the usual spherical harmonics. Recall
  !>  that \f$ Y_{\lambda 0}(\theta,\varphi) = \sqrt{(2\lambda+1)/4\pi}
  !>  P_{\lambda}(\cos\theta)\f$ with \f$ P_{\lambda}(x)\f$ the Legendre polynomial
  !>  of order \f$ \lambda \f$. In this routine, we include deformations of the
  !>  nuclear surface up to \f$ \lambda_{\rm max} =4\f$ (not to be confused with
  !>  the maximum multipolarity of multipole moments used as constraints, which is
  !>  set in @ref hfbtho_variables.f90 to \f$ \lambda_{\rm max} = 8\f$). The volume
  !>  conservation coefficient ensures that the volume is constant as we deform the
  !>  nuclear surface. It is computed based on the condition
  !>
  !>    \f[
  !>       \frac{4}{3}\pi R_{0}^{3} = \int_{0}^{2\pi}d\varphi \int_{0}^{\pi} \sin\theta d\theta
  !>                                  \int_{0}^{R(\theta,\varphi)} r^{2}dr ,
  !>    \f]
  !>
  !>  which leads to
  !>
  !>    \f[
  !>       c(\beta) = \left[ \frac{4\pi}{\displaystyle\int_{0}^{2\pi}d\varphi
  !>                  \int_{0}^{\pi} f^{3}(\theta,\varphi)\sin\theta d\theta} \right]^{1/3} .
  !>    \f]
  !>
  !>  In the case of an octupole deformation, we first compute the position of the
  !>  center of mass according to the definition
  !>
  !>    \f[
  !>       z_{\rm cm} = \frac{1}{M} \int_{0}^{2\pi}d\varphi \int_{0}^{\pi} \sin\theta d\theta
  !>                                 \int_{0}^{R(\theta,\varphi)} z \rho(\boldsymbol{r})r^{2}dr .
  !>    \f]
  !>
  !>  We then assume a constant density inside the nucleus, \f$ \rho(\boldsymbol{r}) = \rho_0 \f$,
  !>  and this leads to the analytical formula
  !>
  !>    \f[
  !>       z_{\rm cm} = \frac{3R_0}{16\pi} c^{4}(\beta)
  !>                    \int_{0}^{2\pi}d\varphi \int_{0}^{\pi} f^{4}(\theta,\varphi)\cos\theta\sin\theta d\theta .
  !>    \f]
  !>
  !>  Calculations of the distance from the current point to the point on the
  !>  surface simply require the substitution \f$ z \rightarrow z - z_{\rm cm} \f$.
  !>  The MATHEMATICA code used to obtain the expression for the variables *fac* and
  !>  *zcm* is given below.
  !>
  !>~~~~~~~~~~~~~~~~~~~~~
  !>      ClearAll[P1, P2, P3, P4, b1, b2, b3, b4, f, c, x, zcm, R];
  !>      P1 := LegendreP[1, x];
  !>      P2 := LegendreP[2, x];
  !>      P3 := LegendreP[3, x];
  !>      P4 := LegendreP[4, x];
  !>      f := 1 + b1*P1 + b2*P2 + b3*P3 + b4*P4;
  !>      c := f^3;
  !>      zzz := x*f^4;
  !>      calpha := Integrate[c, {x, -1, +1}];
  !>      FortranForm[calpha]
  !>      zcm := 3*R0*calpha^4/16/Pi * Integrate[zzz, {x, -1, +1}];
  !>      FortranForm[zcm]
  !>~~~~~~~~~~~~~~~~~~~~~
  !>
  !=======================================================================
  Subroutine start
    !---------------------------------------------------------------------
    ! initializes scratch Saxon-Woods potentials
    !---------------------------------------------------------------------
    Implicit None
    Integer(ipr) :: iw,i,ih,il,ihl,it,ita
    Real(pr) :: zb(ngh),rrb(ngl),rb(ngl),rav,rao,vpws,vls,betas,gamma,fac, &
                facb,zz,rr,r,ctet,cphi,p2,p20,p22,s,u,w,f,rc,c,beta00,     &
                b1_ws,b2_ws,b3_ws,b4_ws,pleg1,pleg2,pleg3,pleg4,zcm
    !----------------------------------------------------------------------------
    ! Re-initializing all again since scratch calculation
    !----------------------------------------------------------------------------
    Call initialize_HFBTHO_SOLVER
    If(ierror_flag.Ne.0) Return
    Call Constraint_or_not(inin_INI,inin,icstr)
    If(ierror_flag.Ne.0) Return
    Do it=itmin,itmax
       ! If blocking active, reset npr to namelist values for initial run
       If(blocking_mode(it).Ne.0 .And. .Not.fragment_properties) Then
          npr(it)=npr_INI(it)
       End If
    End Do
    npr(3)=npr(1)+npr(2)
    If(do_print == 1) Then
       If(blocking_mode(1).Ne.0 .Or. blocking_mode(2).Ne.0) Then
          ! odd nucleus requested but no even-even solution, recalculate the even-even nucleus from scratch
          If(.Not.odd_noBlock) Then
             Do iw=lout,lfile
                Write(iw,'(1x,a,2i4)')
                Write(iw,'(1x,a,2i4)') ' Initialization for the even-even core (N,Z)=: ',npr(1:2)
             End Do
          Else
             Do iw=lout,lfile
                Write(iw,'(1x,a,2i4)')
                Write(iw,'(1x,a,2i4)') ' Initialization without blocking for (N,Z)=: ',npr(1:2)
             End Do
          End If
       Else
          ! scratch for the even-even nucleus requested
          Do iw=lout,lfile
             Write(iw,'(1x,a,2i4)')
             Write(iw,'(a,a,3i4)')    '  Scratch initialization for the nucleus: ',nucname,npr(1:2)
             Write(iw,'(1x,a,2i4)')
          End Do
       End If
    End If
    n00=Abs(n00_INI);  b0=b0_INI;           q=q_INI; iLST=iLST_INI
    maxi=MAX_ITER_INI; inin=inin_INI;
    skyrme=skyrme_INI; kindhfb=kindhfb_INI
    iproj=iproj_INI;   npr1pj=npr1pj_INI;   npr2pj=npr2pj_INI;
    icacou=0; icahartree=0
    !
    Call preparer(.False.)
    !
    If(ierror_flag.Ne.0) Return
    inin=Abs(inin)          ! positive even if inin_INI is not
    If(Abs(b2_0).Gt.3.0_pr) b2_0=3.0_pr ! Avoid crazy initial points (quadrupole deformation)
    If(Abs(b3_0).Gt.1.2_pr) b3_0=1.2_pr ! Avoid crazy initial points (octupole deformation)
    If(Abs(b4_0).Gt.1.0_pr) b4_0=1.0_pr ! Avoid crazy initial points (hexadecapole deformation)
    !-----------------------------------
    ! Saxon-Woods potentials
    !-----------------------------------
    beta00=bet     ! wf to requested deformation
    If(do_print == 1) Then
       Do iw=lout,lfile
          Write(iw,'(/,a)') '  Initial potentials of Saxon-Woods shape '
          Write(iw,'(a,2f14.8)') '  v0ws   =',v0ws
          Write(iw,'(a,2f14.8)') '  kappa  =',akv
          Write(iw,'(a,2f14.8)') '  vs0    =',vso
          Write(iw,'(a,2f14.8)') '  r0     =',r0v
          Write(iw,'(a,2f14.8)') '  a      =',av
          Write(iw,'(a,2f14.8)') '  r0-so  =',rso
          Write(iw,'(a,2f14.8)') '  a-so   =',aso
          Write(iw,'(a,f14.8)')  '  b2_ws  =',b2_0
          Write(iw,'(a,f14.8)')  '  b3_ws  =',b3_0
          Write(iw,'(a,f14.8)')  '  b4_ws  =',b4_0
       End Do
    End If
    !-----------------------------------
    ! Densities
    !-----------------------------------
    Do it=itmin,itmax
       ita=3-it; rav=r0v(it)*amas**p13; rao=rso(it)*amas**p13
       vpws=v0ws*(one-akv*(npr(it)-npr(ita))/amas)
       vls=half*(hqc/amu)**2*vpws*vso(it)
       ! Deformations of the surface
       b1_0 = zero
       b1_ws = b1_0 * Sqrt(3.0_pr/(4.0_pr*pi))
       b2_ws = b2_0 * Sqrt(5.0_pr/(4.0_pr*pi))
       b3_ws = b3_0 * Sqrt(7.0_pr/(4.0_pr*pi))
       b4_ws = b4_0 * Sqrt(9.0_pr/(4.0_pr*pi))
       ! Volume conservation condition
       fac = two+ (two*b4_ws**2)/three + (40.0_pr*b2_ws*b4_ws**2)/231.0_pr + (36.0_pr*b4_ws**3)/1001.0_pr  &
                + (six*b2_ws**2)/five  + (four*b2_ws**3)/35.0_pr + (12.0_pr*b2_ws**2*b4_ws)/35.0_pr &
                + (six*b3_ws**2)/seven + (eight*b2_ws*b3_ws**2)/35.0_pr + (12.0_pr*b3_ws**2*b4_ws)/77.0_pr
       ! z-coordinate of the center of mass, vanishes when b3=0.
       zcm = ((two + (six*b2_ws**2)/five + (four*b2_ws**3)/35.0_pr + (six*b3_ws**2)/seven + (eight*b2_ws*b3_ws**2)/35.0_pr &
            + (12.0_pr*b2_ws**2*b4_ws)/35.0_pr + (12.0_pr*b3_ws**2*b4_ws)/77.0_pr &
            + (two*b4_ws**2)/3.0_pr + (40.0_pr*b2_ws*b4_ws**2)/231.0_pr + (36.0_pr*b4_ws**3)/1001.0_pr)**4 &
            * (three*b3_ws*(13923.0_pr*b2_ws**3 + 68.0_pr*b3_ws**2*(117.0_pr + 122.0_pr*b4_ws) &
             + 34.0_pr*b2_ws**2*(1287.0_pr + 812.0_pr*b4_ws) + 20.0_pr*b4_ws*(2431.0_pr + 867.0_pr*b4_ws + 336.0_pr*b4_ws**2) &
             + 17.0_pr*b2_ws*(3861.0_pr + 621.0_pr*b3_ws**2 + 3692.0_pr*b4_ws + 1281.0_pr*b4_ws**2)))*rav)/(510510.0_pr*pi)
       zcm = (three*(one/(two + (six*b2_ws**2)/five + (four*b2_ws**3)/35.0_pr + (six*b3_ws**2)/seven + (eight*b2_ws*b3_ws**2)/35.0_pr + (12.0_pr*b2_ws**2*b4_ws)/35.0_pr + &
              (12.0_pr*b3_ws**2*b4_ws)/77.0_pr   + (two*b4_ws**2)/three                      + (40.0_pr*b2_ws*b4_ws**2)/231.0_pr          + (36.0_pr*b4_ws**3)/1001.0_pr))**(four/three) &
           * ((72.0_pr*b2_ws*b3_ws)/35.0_pr      + (  48.0_pr*b2_ws**2*b3_ws)/35.0_pr        + (  24.0_pr*b2_ws**3*b3_ws)/55.0_pr         + (  96.0_pr*b3_ws**3)/385.0_pr + (1656.0_pr*b2_ws*b3_ws**3)/5005.0_pr + &
              (32.0_pr*b3_ws*b4_ws)/21.0_pr      + (2272.0_pr*b2_ws*b3_ws*b4_ws)/1155.0_pr   + (1856.0_pr*b2_ws**2*b3_ws*b4_ws)/2145.0_pr + (3904.0_pr*b3_ws**3*b4_ws)/15015.0_pr + &
             (544.0_pr*b3_ws*b4_ws**2)/1001.0_pr + ( 488.0_pr*b2_ws*b3_ws*b4_ws**2)/715.0_pr + ( 512.0_pr*b3_ws*b4_ws**3)/2431.0_pr)*R0)/(two*two**(one/three))
       ! Volume conservation condition
       fac=(two/fac)**(p13)
       ! z,r-coordinates in fm
       zb=xh*bz; rrb=xl*bp*bp; rb=Sqrt(rrb)
       Do ih=1,ngh
          zz=(zb(ih)-zcm)**2
          Do il=1,ngl
             rr=rrb(il)+zz; r=Sqrt(rr); ctet=zz/rr
             ! Deformed surface
             pleg1 = ctet
             pleg2 = half*(three*ctet - one)
             pleg3 = half*(five*ctet**(1.5) - three*Sqrt(ctet))
             pleg4 = (35.0_pr*ctet**2 - 30.0_pr*ctet + three)/eight
             facb=fac*(one + b2_ws*pleg2 + b3_ws*pleg3 + b4_ws*pleg4)
             ! Woods-Saxon potential
             u= vpws/( one+Exp( (r-rav*facb) / av(it) ))
             w=-vls /( one+Exp( (r-rao*facb) / aso(it)))
             ihl=ih+(il-1)*ngh
             If(it.Eq.1) Then
                vhbn(ihl)=hb0n; vn(ihl)=u; vsn(ihl)=w;
                vrn(ihl)=zero; vzn(ihl)=zero; vdn(ihl)=zero;
                vSFIZn(ihl)=zero; vSZFIn(ihl)=zero;
                vSFIRn(ihl)=zero; vSRFIn(ihl)=zero;
             Else
                vhbp(ihl)=hb0p; vp(ihl)=u; vsp(ihl)=w;
                vrp(ihl)=zero; vzp(ihl)=zero; vdp(ihl)=zero;
                vSFIZp(ihl)=zero; vSZFIp(ihl)=zero;
                vSFIRp(ihl)=zero; vSRFIp(ihl)=zero;
             End If
             ro(ihl,it)=u
             aka(ihl,it)=5.0d-3*Exp((r-rav*facb)/2.0_pr)
          End Do
       End Do
       s=tz(it)/Sum(ro(:,it))
       Do il=1,ngl
          Do ih=1,ngh
             ihl=ih+(il-1)*ngh
             f=s/(pi*wh(ih)*wl(il)* bz*bp*bp); ro(ihl,it)=f*ro(ihl,it)
          End Do
       End Do
       !-----------------------------------
       ! pairing
       !-----------------------------------
       Do il=1,nghl
          If(it.Eq.1) Then
             dvn(il)=-100.0_pr*aka(il,it)
          Else
             dvp(il)=-100.0_pr*aka(il,it)
          End If
       End Do
    End Do
    !-----------------------------------
    ! coulomb
    !-----------------------------------
    If(icou.le.0) Then
       cou=zero
    Else
       rc=r0v(2)*amas**p13
       Do il=1,ngl
          Do ih=1,ngh
             r=Sqrt((zb(ih)-zcm)**2+rrb(il))
             If (r.Lt.rc) Then
                c=half*(3/rc-r*r/(rc**3))
             Else
                c=one/r
             End If
             cou(ih+(il-1)*ngh)=c*tz(2)/alphi
          End Do
       End Do
    End If
    !-----------------------------------
    ! initial ph+pp matrix elements
    !-----------------------------------
    ak=0.1_pr; rk=0.1_pr ! initial density matrix elements (improve later)
    brin=zero             ! initial matrix elements to zero
    iiter=0               ! iteration number iiter to zero
    Call gamdel(.true.,.true.)
    !
  End Subroutine start
  !=======================================================================
  !> Display the density matrix \f$ \rho \f$ and pairing tensor \f$ \kappa \f$
  !> as a function of the (dimensionless) radius \f$ r = \sqrt{\eta^2 + \xi^2}\f$
  !=======================================================================
  Subroutine printRHO
    Implicit None
    Integer(ipr), Save :: ihli,Ifle
    Ifle=76+iLST1
#if(DO_MASSTABLE==1 || DRIP_LINES==1 || DO_PES==1)
    If(Ifle.Eq.76) Open(Ifle,file='ho_den'//row_string//'.dat',status='unknown')
    If(Ifle.Eq.77) Open(Ifle,file='tho_den'//row_string//'.dat',status='unknown')
#else
    If(Ifle.Eq.76) Open(Ifle,file='ho_den.dat',status='unknown')
    If(Ifle.Eq.77) Open(Ifle,file='tho_den.dat',status='unknown')
#endif
    Write(Ifle,*) 'r  denN  denP  akaN  akaP '
    Do ihli=1,nghl
       Write(Ifle,'(12(1x,e16.8))') Sqrt(fh(ihli)**2+fl(ihli)**2)  &
            ,ro(ihli,1),ro(ihli,2),aka(ihli,1),aka(ihli,2)
    End Do
    Close(Ifle)
  End Subroutine printRHO
  !=======================================================================
  !> Determines the symbol of the current element (is=1) or the proton
  !> number of a given element (is=2)
  !=======================================================================
  Subroutine nucleus(is,npr2,te)
    Use HFBTHO_utilities, Only: pr,ipr
    Use HFBTHO, Only: ierror_flag,ierror_info
    Implicit None
    Integer(ipr) :: is,npr2,np
    Integer(ipr) :: maxz
    Parameter (maxz=133)
    Character(2) te
    Character(2*maxz+2) t
    T(  1: 40)=' n HHeLiBe B C N O FNeNaMgAlSi P SClAr K'
    T( 41: 80)='CaScTi VCrMnFeCoNiCuZnGaGeASSeBrKrRbSr Y'
    T( 81:120)='ZrNbMoTcRoRhPdAgCdInSnSbTe IXeCsBaLaCePr'
    T(121:160)='NdPmSmEuGdTbDyHoErTmYbLuHfTa WReOsIrPtAu'
    T(161:200)='HgTlPbBiPoAtRnFrRaAcThPa UNpPuAmCmBkCfEs'
    T(201:240)='FmMdNoLrRfDbSgBhHsMtDsRgCnNhFlMcLvTsOgXx'
    T(241:265)='XxXxXxXxXxXxXxXxXxXxXxXxX'
    If (is.Eq.1) Then
       If (npr2.Lt.0.Or.npr2.Gt.maxz) Then
          ierror_flag=ierror_flag+1
          ierror_info(ierror_flag)='STOP: in nucleus npr2 is wrong:'
          Return
       End If
       te=t(2*npr2+1:2*npr2+2)
       Return
    Else
       Do np=0,maxz
          If (te.Eq.t(2*np+1:2*np+2)) Then
             npr2=np
             Return
          End If
       End Do
    End If
    ierror_flag=ierror_flag+1
    ierror_info(ierror_flag)='STOP: in nucleus the nucleus is unknown!'
  End Subroutine nucleus
  !=======================================================================
  !> Given 'Z' returns mass number 'A' on the stability line
  !=======================================================================
  Subroutine stab(npr2,npr3)
    Implicit None
    Integer(ipr) :: npr2,npr3
    Real(pr), Save :: sn,sz,dsn,c,c5
    c=0.0060_pr; c5=5.0_pr*c/3.0_pr; sz=npr2; sn=npr2
    Do While(Abs(dsn).Lt.1.0d-5)
       dsn=sz-sn+c*(sn+sz)**(5.0_pr/3.0_pr)
       dsn=dsn/(-1.0_pr+c5*(sn+sz)**(2.0_pr/3.0_pr))
       sn=sn-dsn
    End Do
    npr3=Int(sn)
  End Subroutine stab
  !=======================================================================
  !> Coulomb field (direct part) Vautherin prescription
  !> Ref.: Phys. Rev. C 7, 296 (1973)
  !=======================================================================
  Subroutine coulom1
    Use EllipticIntegral
    Implicit None
    Integer(ipr), Save :: i,k
    Real(pr) :: zd2,rhl,y1,y2,xx1,xx2,s1,s2,e1,e2,vik,f,r,r1,r4,  &
                rr2,z,z1,zd1,x1,x2,fac1,fac2
    !
    If (debug_solver.Ge.1) Call get_CPU_time('coulom1',0)
    !
    If(icacou.Eq.0) Then
       !
       icacou=1
       !
       ! For parity-breaking shapes, the Coulomb potential was incorrectly
       ! calculated by assuming the two intervals [0,+\infty[ and ]-infty,0]
       ! were equivalent (see also routine coulom() below). This bug was
       ! corrected in version 200d
       If(Parity) Then
          fac1 = one;  fac2 = one
       Else
          fac1 = zero; fac2 = two
       End If
       !
       f=half*chargee2/pi
       ! See notes in subroutine coulom for explanations about some numerical
       ! factors apparently missing here.
!$OMP PARALLEL DO        &
!$OMP& DEFAULT(NONE)     &
!$OMP& SCHEDULE(DYNAMIC) &
!$OMP& SHARED(nghl,fl,fh,fac1,fac2,wdcor,vc,f) &
!$OMP& PRIVATE(i,r,z,r4,k,r1,z1,rr2,rhl,zd1,y1,xx1,s1,zd2,y2,xx2,s2,vik)
       Do i=1,nghl
          r=fl(i); z=fh(i)
          r4=four*r
          Do k=1,i
             r1=fl(k); z1=fh(k)
             rhl=r4*r1     ! 4 r r'
             rr2=(r+r1)**2 ! (r+r')^2
             ! z>0 part
             zd1=(z-z1)**2 ! (z-z')^2
             y1=zd1+rr2    ! d(r,z) = (r+r')^2 + (z-z')^2
             xx1=rhl/y1    ! 4 r r' / d(r,z)
             s1=Sqrt(y1)   ! sqrt(d(r,z))
             ! z<0 part
             zd2=(z+z1)**2
             y2=zd2+rr2
             xx2=rhl/y2
             s2=Sqrt(y2)
             !
             vik = f*fac2*(s1*CompleteEllipticFunction_2nd(xx1) &
                          +s2*CompleteEllipticFunction_2nd(xx2)*fac1)
             !
             vc(i,k)=vik*wdcor(k)  !wdcor=pi*wh*wl*bz*bp*bp
             vc(k,i)=vik*wdcor(i)  !wdcor=pi*wh*wl*bz*bp*bp
             !
          End Do  !k
       End Do  !i
!$OMP End Parallel Do
    End If
    ! Calculation of the coulomb field (each iteration)
    cou=zero
    Call dgemm('n','n',nghl,1,nghl,1.0_pr,vc,nghl,dro(:,2),nghl,0.0_pr,cou,nghl)
    !
    If (debug_solver.Ge.1) Call get_CPU_time('coulom1',1)
    !
  End Subroutine coulom1
  !=======================================================================
  !> Coulomb field (direct part), Gogny prescription
  !> Ref.: Phys. Rev. C 27, 2317 (1983)
  !=======================================================================
  Subroutine coulom
    Use bessik
    Implicit None
    Integer(ipr), Save :: i,j,k
    Real(pr), Save :: zd2,y1,y2,xx1,s1,vik,f,r,r1,fac1,fac2,rr2,z,z1,zd1,t,  &
                      bb,r2,r12,rrr,rz1,rz2,rrz1,rrz2,xx,rk1,rip1,rkp1,alpha,&
                      beta,xxx
    !
    If (debug_solver.Ge.1) Call get_CPU_time('coulom',0)
    !
    If(icacou.Eq.0) Then
       !
       icacou=1
       !
       ! For parity-breaking shapes, the Coulomb potential was incorrectly
       ! calculated by assuming the two intervals [0,+\infty[ and ]-infty,0]
       ! were equivalent (see also below). This bug was corrected in version
       ! 139a
       If(Parity) Then
          fac1 = one;  fac2 = one
       Else
          fac1 = zero; fac2 = two
       End If
       ! Notes:
       !   - Missing factor 2 compared to Eq. (58) CPC paper because the density
       !     ro(:,it) already contains it (see routine DENSIT) due to T-invariance
       !   - Missing factor 1/2 when applying Gauss-Legendre quadrature (from [0,1]
       !     to the proper [-1,1] interval because it will be put back in subroutine
       !     expect() and is cancelled by a factor 2 in the HF field
       !   - For conserved parity, Gauss-Hermite points are all positive, the full
       !     integral over z' is split in z'<0 and z'>0, values of z and z1 below
       !     refer to the absolute values of z' (=-z' if z'<0)
       !
       bb=50.0_pr          ! Length scale L
       beta=2.00_pr
       alpha=one/beta
       f=chargee2/Sqrt(pi) ! e^2/Sqrt(pi)
       !
!$OMP PARALLEL DO        &
!$OMP& DEFAULT(NONE)     &
!$OMP& SCHEDULE(DYNAMIC) &
!$OMP& SHARED(nghl,fl,fh,nleg,xleg,bb,fac1,fac2,wleg,wdcor,vc,f,alpha,beta) &
!$OMP& PRIVATE(i,r,z,k,r1,z1,rrr,rr2,zd1,zd2,rz1,rz2,rrz1,rrz2, &
!$OMP&         xx1,j,xx,y1,s1,t,y2,vik,xxx)
       Do i=1,nghl
          r = fl(i); z = fh(i)
          Do k=1,i
             !
             r1 = fl(k); z1 = fh(k)
             rrr = two*r*r1; rr2 = (r - r1)**2
             ! z>0 part
             zd1 = (z - z1)**2
             rz1 = rr2 + zd1
             ! z<0 part
             zd2 = (z + z1)**2
             rz2 = rr2 + zd2
             ! Gauss-Legendre integration over u from 0 to D
             xx1=zero
             Do j=1,nleg
                xx=(one-xleg(j)**beta)**alpha ! change of variable to 0 <= u <= 1
                xxx=(one-xleg(j)**beta)**(alpha+one)
                y1=(xleg(j)/(bb*xx))**2 ! u^2
                s1=y1*rrr               ! 2 u^2 r r'
                y2=besei0(s1)           ! I0( 2 u^2 r r' ) * exp(-2 u^2 r r')
                xx1=xx1+fac2*wleg(j)*y2*(Exp(-rz1*y1) + fac1*Exp(-rz2*y1)) / xxx
             End Do
             vik=f*xx1/bb
             !
             vc(i,k)=vik*wdcor(k)  !wdcor=pi*wh*wl*bz*bp*bp
             vc(k,i)=vik*wdcor(i)  !wdcor=pi*wh*wl*bz*bp*bp
             !
          End Do  !k
       End Do  !i
!$OMP End Parallel Do
       !
    End If

    ! Calculation of the Coulomb field
    cou=zero
    Call dgemm('n','n',nghl,1,nghl,1.0_pr,vc,nghl,ro(:,2),nghl,0.0_pr,cou,nghl)
    !
    If (debug_solver.Ge.1) Call get_CPU_time('coulom',1)
    !
  End Subroutine coulom
  !=======================================================================
  ! Coulomb field (direct part), Gogny prescription
  ! Ref.: Phys. Rev. C 27, 2317 (1983)
  !=======================================================================
  Subroutine coulom_test
    Use bessik
    Implicit None
    Integer(ipr), Save :: i,j,k
    Real(pr), Save :: zd2,y1,y2,xx1,s1,vik,f,r,r1,fac1,fac2,rr2,z,z1,zd1,t,  &
                      bb,r2,r12,rrr,rz1,rz2,rrz1,rrz2,xx,rk1,rip1,rkp1,alpha,&
                      beta,xxx,func
    !
    If (debug_solver.Ge.1) Call get_CPU_time('coulom_test',0)
    !
    !
    ! For parity-breaking shapes, the Coulomb potential was incorrectly
    ! calculated by assuming the two intervals [0,+\infty[ and ]-infty,0]
    ! were equivalent (see also below). This bug was corrected in version
    ! 139a
    If(Parity) Then
       fac1 = one;  fac2 = one
    Else
       fac1 = zero; fac2 = two
    End If
    !
    bb=5.0_pr          ! Length scale L
    beta=2.00_pr
    alpha=one/beta
    !f=chargee2/Sqrt(pi) ! e^2/Sqrt(pi)
    f=one/Sqrt(pi)       ! 1/Sqrt(pi)
    !
    Do j=1,nleg
       ! Gauss-Legendre integration over u from 0 to D
       xx=(one-xleg(j)**beta)**alpha ! change of variable to 0 <= u <= 1
       xxx=(one-xleg(j)**beta)**(alpha+one)
       !
       func=zero
       Do i=1,nghl
          r = fl(i); z = fh(i)
          Do k=1,i
             !
             r1 = fl(k); z1 = fh(k)
             rrr = two*r*r1; rr2 = (r - r1)**2
             ! z>0 part
             zd1 = (z - z1)**2
             rz1 = rr2 + zd1
             ! z<0 part
             zd2 = (z + z1)**2
             rz2 = rr2 + zd2
             y1=(xleg(j)/(bb*xx))**2 ! u^2
             s1=y1*rrr               ! 2 u^2 r r'
             y2=besei0(s1)           ! I0( 2 u^2 r r' ) * exp(-2 u^2 r r')
             xx1=fac2*wleg(j)*y2*(Exp(-rz1*y1) + fac1*Exp(-rz2*y1)) / xxx
             vik=f*xx1/bb
             !
             func=func+vik*wdcor(k)*ro(k,2)*wdcor(i)*ro(i,2)  !wdcor=pi*wh*wl*bz*bp*bp
             !
          End Do ! k
       End Do  ! i
       Write(6,'(2f30.14)') xleg(j),func
       !
    End Do  !j
    !
    If (debug_solver.Ge.1) Call get_CPU_time('coulom_test',1)
    !
  End Subroutine coulom_test
  !=======================================================================
  !> Optimization arrays
  !>    FI2D_opt(JA,ihil) == Laplacian(r,z) HOwf
  !>    FID2D-xlamy2*FID  == Laplacian(r,z,phy) FID
  !>    FIU2D-xlapy2*FIU  == Laplacian(r,z,phy) FIU
  !=======================================================================
  Subroutine optHFBTHO
    Implicit None
    Integer(ipr) :: i,ih,il,ib,ibx,nd,nza,nra,nla,nsa
    Integer(ipr) :: ihil,laplus,im,JA,N1,N2,ndnd,n12,n21
    Real(pr)    :: qla,v2,v4,yi,y,y2,qha,qhla,xmi,u,u2,un,up,xxx
    Real(pr)    :: sml2,cnzaa,cnraa,a,b
    Real(pr)    :: FITW1,FITW2,FITW3,FITW4
    Real(pr)    :: fi1r,fi1z,fi2d,QHL1A,QH1LA,vh,vdh,vsh,hbh
    Real(pr)    :: SRFIh,SFIRh,SFIZh,SZFIh,SNABLARh,SNABLAZh
    Real(pr)    :: xlam,xlam2,xlamy,xlamy2,xlap,xlap2,xlapy,xlapy2,XLAMPY
    Real(pr)    :: bpi,bpi2,bzi,bzi2,xh2
    !
    bpi=one/bp; bpi2=bpi*bpi; bzi=one/bz; bzi2=bzi*bzi
    !
    !-----------------------------------------
    ! Allocate the optimization arrays
    !-----------------------------------------
    If(Allocated(QHLA_opt)) Deallocate(QHLA_opt,FI1R_opt,FI1Z_opt,FI2D_opt,y_opt)
    Allocate(QHLA_opt(ntx,nghl),FI1R_opt(ntx,nghl),FI1Z_opt(ntx,nghl),FI2D_opt(ntx,nghl),y_opt(nghl))
    !----------------------------------------------
    ! START BLOCKS
    !----------------------------------------------
    Do ib=1,NB
       ND=ID(ib); IM=ia(ib)
       If(Parity) Then
          LAPLUS=(ib+1)/2 !Yesp
       Else
          LAPLUS=ib       !Nop
       End If
       XLAP=LAPLUS; XLAM=XLAP-ONE; xlap2=xlap*xlap; xlam2=xlam*xlam
       !----------------------------------------------
       ! SUM OVER GAUSS INTEGRATION POINTS
       !----------------------------------------------
       Do IL=1,ngl
          v2=half/xl(il); v4=v2*v2
          Do IH=1,ngh
             ihil=ih+(il-1)*ngh; xh2=xh(ih)**2
             If(iLST1.Eq.0) Then
                ! HO-basis
                yi=Sqrt(xl(il))*bp; y=one/yi; y2=y*y
                xlamy=xlam*y; xlamy2=xlam2*y2; xlapy=xlap*y; xlapy2=xlap2*y2; XLAMPY=XLAMY+XLAPY
             Else
                ! THO-basis
                y=fli(ihil); y2=y*y; xlamy=xlam*y; u=xh(ih); u2=u*u;
                xlamy2=xlam2*y2; xlapy=xlap*y; xlapy2=xlap2*y2; XLAMPY=XLAMY+XLAPY
             End If
             y_opt(ihil)=y
             !----------------------------------------------
             ! SCAN OVER BASIS STATES
             !----------------------------------------------
             Do N1=1,ND
                JA=N1+IM; NLA=NL(JA); NRA=NR(JA); NZA=NZ(JA); NSA=NS(JA)
                SML2=NLA*NLA; CNZAA=NZA+NZA+1; CNRAA=NRA+NRA+NLA+1
                QHA=QH(NZA,IH); QLA=QL(NRA,NLA,IL); QHLA=QHA*QLA
                QHL1A=QHA*QL1(NRA,NLA,IL)*V2; QH1LA=QH1(NZA,IH)*QLA
                If(iLST1.Eq.0) Then
                   ! HO-basis
                   FI1R=(two*Sqrt(xl(il))*bpi)*QHL1A
                   FI1Z=bzi*QH1LA
                   FI2D=((xh2-CNZAA)*bzi2+four*(p14-CNRAA*V2+SML2*V4)*xl(il)*bpi2 )*QHLA
                Else
                   ! THO-basis
                   u=xh(ih); u2=u*u;
                   FI1R=FP4(IHIL)*QHLA+FP5(IHIL)*QH1LA+FP6(IHIL)*QHL1A
                   FI1Z=FP1(IHIL)*QHLA+FP2(IHIL)*QH1LA+FP3(IHIL)*QHL1A
                   FI2D=(FS1(IHIL)*QH1LA*QH1LA+FS2(IHIL)*QHL1A*QHL1A  &
                        +FOUR*FS4(IHIL)*QH1LA*QHL1A  &
                        +TWO*(FS5(IHIL)*QH1LA+FS6(IHIL)*QHL1A)*QHLA  &
                        +((U2-CNZAA)*FS1(IHIL)+(p14-CNRAA*V2+SML2*V4)*FS2(IHIL)  &
                        +FS3(IHIL))*QHLA*QHLA-TWO*(FI1R*FI1R+FI1Z*FI1Z))/(TWO*QHLA)
                End If
                QHLA_opt(JA,ihil)=QHLA; FI2D_opt(JA,ihil)=FI2D; FI1R_opt(JA,ihil)=FI1R; FI1Z_opt(JA,ihil)=FI1Z
             End Do !N1
             !
          End Do !IH
       End Do !IL
    End Do !IB
  End Subroutine optHFBTHO
  !=======================================================================
  !> Computes all Skyrme densities at Gauss quadrature points
  !=========================================================================
  Subroutine DENSIT
    Implicit None
    Integer(ipr) :: iw,nsa,nza,nra,nla,k,i,nd,il,ih,ihil,laplus,ii
    Integer(ipr) :: imen,ib,im,it,J,JJ,JA,JN,k0,k1,k2,ibiblo
    Integer(ipr) :: bb,size,ndxmax
    Parameter(ndxmax=(n00max+2)*(n00max+2)/4)
    Real(pr)     :: s,ss,sd,yi,y,y2,sml2,cnzaa,cnraa,u,u2,v2,v4
    Real(pr)     :: anik,pnik,qhla,qh1la,qhl1a,qla,qha,fi1r,fi1z,fi2d,fidd
    Real(pr)     :: xlam,xlam2,xlamy,xlamy2,xlap,xlap2,xlapy,xlapy2,XLAMPY
    Real(pr)     :: TFIU,TFID,TFIUR,TFIDR,TFIUZ,TFIDZ,TFIUD2,TFIDD2
    Real(pr)     :: TPFIU,TPFID,TPFIUR,TPFIDR,TPFIUZ,TPFIDZ,TPFIUD2,TPFIDD2
    Real(pr)     :: PIU,PIUZ,PIUR,PIUD2,PID,PIDZ,PIDR,PIDD2
    Real(pr)     :: TEMP1,TEMP2,TEMP3,TEMP4,TEMP5,TEMP6,TEMP7,TEMP8,TEMP9,TEMP10,TEMP11,TW_T,PW_T,WGT(nghl)
    Real(pr)     :: Takaihil,Troihil,Tdjihil,Ttauihil,Tdroihil,TSRFIihil
    Real(pr)     :: TSFIRihil,TSFIZihil,TSZFIihil,TNABLARIHIL,TNABLAZIHIL
    Real(pr), Pointer :: TAKA(:),TRO(:),TDJ(:),TTAU(:),TDRO(:)
    Real(pr), Pointer :: TSRFI(:),TSFIR(:),TSFIZ(:),TSZFI(:),TNABLAR(:),TNABLAZ(:)
    Real(pr)     :: time1,time2,fk,f1k
    Real(pr), Pointer     :: EqpPo(:),VqpPo(:),UqpPo(:)
    Integer(ipr), Pointer :: KpwiPo(:),KqpPo(:)
#if(USE_OPENMP==1)
    Real(pr), Allocatable :: OMPTAKA(:,:,:),OMPTRO(:,:,:),OMPTDJ(:,:,:),OMPTTAU(:,:,:),OMPTDRO(:,:,:)
    Real(pr), Allocatable :: OMPTSRFI(:,:,:),OMPTSFIR(:,:,:),OMPTSFIZ(:,:,:),OMPTSZFI(:,:,:), &
                             OMPTSZIF(:,:,:),OMPTNABLAR(:,:,:),OMPTNABLAZ(:,:,:)
#else
    Real(pr) :: OMPTAKA(nghl,2),OMPTRO(nghl,2),OMPTDJ(nghl,2),OMPTTAU(nghl,2),OMPTDRO(nghl,2)
    Real(pr) :: OMPTSRFI(nghl,2),OMPTSFIR(nghl,2),OMPTSFIZ(nghl,2),OMPTSZFI(nghl,2), &
                OMPTSZIF(nghl,2),OMPTNABLAR(nghl,2),OMPTNABLAZ(nghl,2)
#endif
    !
    Real(pr) :: OMPFIU(ndxmax),OMPFID(ndxmax),OMPFIUR(ndxmax),OMPFIDR(ndxmax),OMPFIUZ(ndxmax)
    Real(pr) :: OMPFIDZ(ndxmax),OMPFIUD2N(ndxmax),OMPFIDD2N(ndxmax)
    !
    Real(pr) :: OMPPFIU(ndxmax),OMPPFID(ndxmax),OMPPFIUR(ndxmax),OMPPFIDR(ndxmax),OMPPFIUZ(ndxmax)
    Real(pr) :: OMPPFIDZ(ndxmax),OMPPFIUD2N(ndxmax),OMPPFIDD2N(ndxmax)
    !
    Real(pr) :: OMPAN(ndxmax*ndxmax),OMPANK(ndxmax*ndxmax)
    Real(pr) :: f_T(ndxmax),f1_T(ndxmax)
    Real(pr) :: dnrm2
    external dnrm2
    !
    If (debug_solver.Ge.1) Call get_CPU_time('densit',0)
    If( .Not.Allocated(ro_normalization)) Allocate(ro_normalization(nghl,2))
    !
    !-----------------------------------------------
    ! ZERO N & P DENSITIES
    !-----------------------------------------------
    ro=zero; tau=zero; dj=zero; dro=zero; aka=zero; szfi=zero; sfiz=zero
    srfi=zero; sfir=zero; nablar=zero; nablaz=zero; varmas=zero
    !
#if(USE_OPENMP==1)
    Allocate(OMPTAKA(1:NB,1:nghl,2),OMPTRO(1:NB,1:nghl,2),OMPTDJ(1:NB,1:nghl,2),&
             OMPTTAU(1:NB,1:nghl,2),OMPTDRO(1:NB,1:nghl,2))
    Allocate(OMPTSRFI(1:NB,1:nghl,2),OMPTSFIR(1:NB,1:nghl,2),OMPTSFIZ(1:NB,1:nghl,2), &
             OMPTSZFI(1:NB,1:nghl,2),OMPTSZIF(1:NB,1:nghl,2),OMPTNABLAR(1:NB,1:nghl,2),&
             OMPTNABLAZ(1:NB,1:nghl,2))
#endif
    OMPTAKA   = ZERO; OMPTRO     = ZERO; OMPTDJ     = ZERO; OMPTTAU  = ZERO
    OMPTDRO   = ZERO; OMPTSRFI   = ZERO; OMPTSFIR   = ZERO; OMPTSFIZ = ZERO
    OMPTSZFI  = ZERO; OMPTNABLAR = ZERO; OMPTNABLAZ = ZERO;
    !
#if(USE_OPENMP==1)
!$OMP PARALLEL DO        &
!$OMP& DEFAULT(NONE)     &
!$OMP& SCHEDULE(DYNAMIC) &
!$OMP& SHARED(NB,npr_INI,ID,ia,nghl,NS,Parity,bloblo,keyblo,blo123d,ka,kd,KpwiN,RVqpN,RUqpN,&
!$OMP&        fn_T,KpwiP,RVqpP,RUqpP,fp_T,y_opt,QHLA_opt,FI2D_opt,FI1Z_opt,FI1R_opt,&
!$OMP&        switch_on_temperature,OMPTaka,OMPTro,OMPTdj,OMPTtau,OMPTdro,OMPTSRFI,&
!$OMP&        OMPTSFIR,OMPTSFIZ,OMPTSZFI,OMPTNABLAR,OMPTNABLAZ) &
!$OMP& PRIVATE(bb,it,ib,ND,IM,LAPLUS,XLAP,XLAM,xlap2,xlam2,ibiblo,K0,k1,k2,imen,&
!$OMP&         ompan,ompank,f_T,f1_T,J,JJ,K,ihil,y,xlamy,xlapy,XLAMPY,y2,xlamy2,&
!$OMP&         xlapy2,OMPFIU,OMPFIUZ,OMPFIUR,OMPFID,OMPFIDZ,OMPFIDR,OMPFIUD2N,OMPFIDD2N,&
!$OMP&         OMPPFIU,OMPPFIUZ,OMPPFIUR,OMPPFID,OMPPFIDZ,OMPPFIDR,OMPPFIUD2N,OMPPFIDD2N,&
!$OMP&         PIU,PIUZ,PIUR,PIUD2,PID,PIDZ,PIDR,PIDD2,JN,JA,NSA,QHLA,PNIK,FI2D,FI1Z,FI1R,&
!$OMP&         Takaihil,Troihil,Tdjihil,Ttauihil,Tdroihil,TSRFIihil,&
!$OMP&         TFIU,TFID,TFIUR,TFIDR,TFIUZ,TFIDZ,TFIUD2,TFIDD2,fk,f1k,&
!$OMP&         TPFIU,TPFID,TPFIUR,TPFIDR,TPFIUZ,TPFIDZ,TPFIUD2,TPFIDD2,&
!$OMP&         TEMP1,TEMP2,TEMP3,TEMP4,TEMP5,TEMP6,TEMP7,TEMP8,TEMP9,TEMP10,TEMP11,TW_T,PW_T,&
!$OMP&         TSFIRihil,TSFIZihil,TSZFIihil,TNABLARIHIL,TNABLAZIHIL)
#endif
    Do bb=0,2*NB-1
       it = bb/NB + 1
       ib = Mod(bb,NB)+1
       !
       ! case of zero particle number, only flush densities
       If((npr_INI(1).Eq.0).And.(it.Eq.1)) Cycle
       If((npr_INI(2).Eq.0).And.(it.Eq.2)) Cycle
       !-----------------------------------------------
       ! SCAN OVER BLOCKS
       !-----------------------------------------------
       ND=ID(ib); IM=ia(ib)
       If(Parity) Then
          LAPLUS=(ib+1)/2 !Yesp
       Else
          LAPLUS=ib       !Nop
       End If
       XLAP=LAPLUS; XLAM=XLAP-ONE; xlap2=xlap*xlap; xlam2=xlam*xlam
       !
       ! blocking
       ibiblo=bloblo(keyblo(it),it)
       K0=0; If(ibiblo.Eq.ib) K0=blo123d(it)
       !
       !----------------------------------------------
       ! PAIRING WINDOW QP WAVE FUNCTIONS
       !----------------------------------------------
       k1=ka(ib,it)+1; k2=ka(ib,it)+kd(ib,it); imen=k2-k1+1
       If(IMEN.Gt.0) Then
          ompan=ZERO; ompank=ZERO; f_T=ZERO; f1_T=ZERO
          J=0
          If(it.Eq.1) then
             Do JJ=1,nd ! basis
                Do K=K1,K2 ! qp
                   J=J+1; I=KpwiN(K)+JJ; ompan(J)=RVqpN(I); ompank(J)=RUqpN(I)
                End Do
             End Do
             J=0
             Do K=K1,K2
                J=J+1;JJ=K !KpwiN(K)
                f_T(J)=one-fn_T(JJ);f1_T(J)=fn_T(JJ)
             End Do
          Else
             Do JJ=1,nd ! basis
                Do K=K1,K2 ! qp
                   J=J+1; I=KpwiP(K)+JJ; ompan(J)=RVqpP(I); ompank(J)=RUqpP(I)
                End Do
             End Do
             J=0
             Do K=K1,K2
                J=J+1;JJ=K !KpwiP(K)
                f_T(J)=one-fp_T(JJ);f1_T(J)=fp_T(JJ)
             End Do
          End If
          !-----------------------------------------------
          ! SCAN OVER GAUSS INTEGRATION POINTS
          !-----------------------------------------------
          Do ihil=1,nghl
             y=y_opt(ihil); xlamy =xlam*y;    xlapy =xlap*y;   XLAMPY=XLAMY+XLAPY
             y2=y*y;        xlamy2=xlam2*y2;  xlapy2=xlap2*y2
             Do K=1,IMEN
                ! V_k components
                OMPFIU(K)    = ZERO; OMPFIUZ(K)   = ZERO; OMPFIUR(K) = ZERO
                OMPFID(K)    = ZERO; OMPFIDZ(K)   = ZERO; OMPFIDR(K) = ZERO
                OMPFIUD2N(K) = ZERO; OMPFIDD2N(K) = ZERO;
                ! U_k components
                OMPPFIU(K)    = ZERO; OMPPFIUZ(K)   = ZERO; OMPPFIUR(K) = ZERO
                OMPPFID(K)    = ZERO; OMPPFIDZ(K)   = ZERO; OMPPFIDR(K) = ZERO
                OMPPFIUD2N(K) = ZERO; OMPPFIDD2N(K) = ZERO;
             End Do
             If(K0.Ne.0) Then
                PIU=ZERO;  PIUZ=ZERO; PIUR=ZERO; PIUD2=ZERO
                PID=ZERO;  PIDZ=ZERO; PIDR=ZERO; PIDD2=ZERO
             End If
             !-----------------------------------------------
             ! SUM OVER BASIS STATES
             !-----------------------------------------------
             JN=0
             Do I=1,ND
                JA=IM+I; NSA=NS(JA); JN=(I-1)*imen
                QHLA=QHLA_opt(JA,ihil); FI2D=FI2D_opt(JA,ihil)
                FI1Z=FI1Z_opt(JA,ihil); FI1R=FI1R_opt(JA,ihil)
                !-----------------------------------------------
                ! QUASIPARTICLE WF IN COORDINATE SPACE
                !-----------------------------------------------
                If (NSA.Gt.0) Then
                   ! SPIN Up
                   Call DAXPY(IMEN,-QHLA,OMPANK(JN+1),1,OMPPFIU,1)
                   ! temperature
                   If(switch_on_temperature) Then
                      Call DAXPY(IMEN,-FI2D,OMPANK(JN+1),1,OMPPFIUD2N,1)
                      Call DAXPY(IMEN,-FI1R,OMPANK(JN+1),1,OMPPFIUR,1)
                      Call DAXPY(IMEN,-FI1Z,OMPANK(JN+1),1,OMPPFIUZ,1)
                   End If
                   Call DAXPY(IMEN, QHLA,OMPAN(JN+1) ,1,OMPFIU,1)
                   Call DAXPY(IMEN, FI2D,OMPAN(JN+1) ,1,OMPFIUD2N,1)
                   Call DAXPY(IMEN, FI1R,OMPAN(JN+1) ,1,OMPFIUR,1)
                   Call DAXPY(IMEN, FI1Z,OMPAN(JN+1) ,1,OMPFIUZ,1)
                   ! blocking
                   If(K0.Ne.0) Then
                      PNIK  = OMPANK(JN+K0)
                      PIU   = PIU   + PNIK*QHLA
                      PIUD2 = PIUD2 + PNIK*FI2D
                      PIUR  = PIUR  + PNIK*FI1R
                      PIUZ  = PIUZ  + PNIK*FI1Z
                   End If
                Else
                   ! SPIN Down
                   Call DAXPY(IMEN,-QHLA,OMPANK(JN+1),1,OMPPFID,1)
                   ! temperature
                   If(switch_on_temperature) Then
                      Call DAXPY(IMEN,-FI2D,OMPANK(JN+1),1,OMPPFIDD2N,1)
                      Call DAXPY(IMEN,-FI1R,OMPANK(JN+1),1,OMPPFIDR,1)
                      Call DAXPY(IMEN,-FI1Z,OMPANK(JN+1),1,OMPPFIDZ,1)
                   End If
                   Call DAXPY(IMEN, QHLA,OMPAN(JN+1) ,1,OMPFID,1)
                   Call DAXPY(IMEN, FI2D,OMPAN(JN+1) ,1,OMPFIDD2N,1)
                   Call DAXPY(IMEN, FI1R,OMPAN(JN+1) ,1,OMPFIDR,1)
                   Call DAXPY(IMEN, FI1Z,OMPAN(JN+1) ,1,OMPFIDZ,1)
                   ! blocking
                   If(K0.Ne.0) Then
                      PNIK  = OMPANK(JN+K0)
                      PID   = PID   + PNIK*QHLA
                      PIDD2 = PIDD2 + PNIK*FI2D
                      PIDR  = PIDR  + PNIK*FI1R
                      PIDZ  = PIDZ  + PNIK*FI1Z
                   End If
                End If
             End Do ! I=1,ND
             !-----------------------------------------------
             ! DENSITIES IN COORDINATE SPACE
             !-----------------------------------------------
             Takaihil=zero;    Troihil=zero;    Tdjihil=zero;   Ttauihil=zero;  Tdroihil=zero
             TSRFIihil=zero;   TSFIRihil=zero;  TSFIZihil=zero; TSZFIihil=zero; TNABLARIHIL=zero; TNABLAZIHIL=zero
             !
             Do K=1,IMEN
                TFIU=OMPFIU(K); TFIUZ=OMPFIUZ(K); TFIUR=OMPFIUR(K); TFIUD2=OMPFIUD2N(K); TPFIU=OMPPFIU(K)
                TFID=OMPFID(K); TFIDZ=OMPFIDZ(K); TFIDR=OMPFIDR(K); TFIDD2=OMPFIDD2N(K); TPFID=OMPPFID(K)
                !
                If(switch_on_temperature) Then
                   !
                   fk=f_T(K); f1k=f1_T(K)
                   !
                   TPFIUZ=OMPPFIUZ(K); TPFIUR=OMPPFIUR(K); TPFIUD2=OMPPFIUD2N(K)
                   TPFIDZ=OMPPFIDZ(K); TPFIDR=OMPPFIDR(K); TPFIDD2=OMPPFIDD2N(K)
                   !
                   TEMP1  = (TPFIU*TFIU+TPFID*TFID)*fk-(TFIU*TPFIU+TFID*TPFID)*f1k
                            TAKAIHIL = TAKAIHIL + TEMP1
                   TEMP2  = (TFIU*TFIU+TFID*TFID)*fk+(TPFIU*TPFIU+TPFID*TPFID)*f1k
                            TROIHIL = TROIHIL + TEMP2
                   TEMP3  = (TFIUR *TFIDZ -TFIDR *TFIUZ +XLAMY*TFIU *(TFIUR -TFIDZ) -XLAPY*TFID *(TFIDR +TFIUZ)) *fk &
                          + (TPFIUR*TPFIDZ-TPFIDR*TPFIUZ+XLAMY*TPFIU*(TPFIUR-TPFIDZ)-XLAPY*TPFID*(TPFIDR+TPFIUZ))*f1k
                            TDJIHIL = TDJIHIL + TEMP3
                   !
                   TW_T=(TFIUR *TFIUR +TFIDR *TFIDR +TFIUZ *TFIUZ +TFIDZ *TFIDZ)*fk&
                       +(TPFIUR*TPFIUR+TPFIDR*TPFIDR+TPFIUZ*TPFIUZ+TPFIDZ*TPFIDZ)*f1k
                   !
                   TEMP4  = (XLAMY2*TFIU *TFIU +XLAPY2*TFID *TFID) *fk &
                          + (XLAMY2*TPFIU*TPFIU+XLAPY2*TPFID*TPFID)*f1k + TW_T
                            TTAUIHIL = TTAUIHIL + TEMP4
                   TEMP5  = (TFIU*TFIUD2+TFID*TFIDD2)*fk + (TPFIU*TPFIUD2+TPFID*TPFIDD2)*f1k + TW_T
                            TDROIHIL = TDROIHIL + TEMP5
                   TEMP6  = (TFIUR*TFID-TFIDR*TFIU)*fk + (TPFIUR*TPFID-TPFIDR*TPFIU)*f1k
                            TSRFIIHIL = TSRFIIHIL + TEMP6
                   TEMP7  = (TFIU*TFID*XLAMPY)*fk + (TPFIU*TPFID*XLAMPY)*f1k
                            TSFIRIHIL = TSFIRIHIL + TEMP7
                   TEMP8  = (XLAMY*TFIU*TFIU-XLAPY*TFID*TFID)*fk + (XLAMY*TPFIU*TPFIU-XLAPY*TPFID*TPFID)*f1k
                            TSFIZIHIL = TSFIZIHIL + TEMP8
                   TEMP9  = (TFIUZ*TFID-TFIDZ*TFIU)*fk + (TPFIUZ*TPFID-TPFIDZ*TPFIU)*f1k
                            TSZFIIHIL = TSZFIIHIL + TEMP9
                   TEMP10 = (TFIUR*TFIU+TFIDR*TFID)*fk + (TPFIUR*TPFIU+TPFIDR*TPFID)*f1k
                            TNABLARIHIL = TNABLARIHIL + TEMP10
                   TEMP11 = (TFIUZ*TFIU+TFIDZ*TFID)*fk + (TPFIUZ*TPFIU+TPFIDZ*TPFID)*f1k
                            TNABLAZIHIL = TNABLAZIHIL + TEMP11
                   !
                Else
                   !
                   TEMP1  = TPFIU*TFIU+TPFID*TFID;                  TAKAIHIL    = TAKAIHIL   + TEMP1
                   TEMP2  = TFIU*TFIU+TFID*TFID;                    TROIHIL     = TROIHIL    + TEMP2
                   TEMP3  = TFIUR*TFIDZ-TFIDR*TFIUZ  &
                           +XLAMY*TFIU*(TFIUR-TFIDZ) &
                           -XLAPY*TFID*(TFIDR+TFIUZ) ;              TDJIHIL     = TDJIHIL    + TEMP3
                   !
                   TW_T=TFIUR*TFIUR+TFIDR*TFIDR+TFIUZ*TFIUZ+TFIDZ*TFIDZ
                   !
                   TEMP4  = XLAMY2*TFIU*TFIU+XLAPY2*TFID*TFID+TW_T; TTAUIHIL    = TTAUIHIL   + TEMP4
                   TEMP5  = TFIU*TFIUD2+TFID*TFIDD2          +TW_T; TDROIHIL    = TDROIHIL   + TEMP5
                   TEMP6  = TFIUR*TFID-TFIDR*TFIU;                  TSRFIIHIL   = TSRFIIHIL  + TEMP6
                   TEMP7  = TFIU*TFID*XLAMPY;                       TSFIRIHIL   = TSFIRIHIL  + TEMP7
                   TEMP8  = XLAMY*TFIU*TFIU-XLAPY*TFID*TFID;        TSFIZIHIL   = TSFIZIHIL  + TEMP8
                   TEMP9  = TFIUZ*TFID-TFIDZ*TFIU;                  TSZFIIHIL   = TSZFIIHIL  + TEMP9
                   TEMP10 = TFIUR*TFIU+TFIDR*TFID;                  TNABLARIHIL = TNABLARIHIL+ TEMP10
                   TEMP11 = TFIUZ*TFIU+TFIDZ*TFID;                  TNABLAZIHIL = TNABLAZIHIL+ TEMP11
                   !
                End If
                !
                If(K.Ne.K0) Cycle
                !
                ! blocking
                TAKAIHIL    = TAKAIHIL    - TEMP1;                 TEMP1  = PIU*PIU+PID*PID
                TROIHIL     = TROIHIL     - HALF*(TEMP2 - TEMP1);  TEMP2  = PIUR*PIDZ-PIDR*PIUZ+XLAMY*PIU*(PIUR-PIDZ) &
                                                                                               -XLAPY*PID*(PIDR+PIUZ)
                !
                PW_T=PIUR*PIUR+PIDR*PIDR+PIUZ*PIUZ+PIDZ*PIDZ
                TDJIHIL     = TDJIHIL     - HALF*(TEMP3 - TEMP2);  TEMP3  = PW_T+XLAMY2*PIU*PIU+XLAPY2*PID*PID
                TTAUIHIL    = TTAUIHIL    - HALF*(TEMP4 - TEMP3);  TEMP4  = PW_T+PIU*PIUD2+PID*PIDD2;
                TDROIHIL    = TDROIHIL    - HALF*(TEMP5 - TEMP4);  TEMP5  = PIUR*PID-PIDR*PIU;
                TSRFIIHIL   = TSRFIIHIL   - HALF*(TEMP6 - TEMP5);  TEMP6  = PIU*PID*XLAMPY;
                TSFIRIHIL   = TSFIRIHIL   - HALF*(TEMP7 - TEMP6);  TEMP7  = XLAMY*PIU*PIU-XLAPY*PID*PID;
                TSFIZIHIL   = TSFIZIHIL   - HALF*(TEMP8 - TEMP7);  TEMP8  = PIUZ*PID-PIDZ*PIU;
                TSZFIIHIL   = TSZFIIHIL   - HALF*(TEMP9 - TEMP8);  TEMP9  = PIUR*PIU+PIDR*PID;
                TNABLARIHIL = TNABLARIHIL - HALF*(TEMP10- TEMP9);  TEMP10 = PIUZ*PIU+PIDZ*PID;
                TNABLAZIHIL = TNABLAZIHIL - HALF*(TEMP11- TEMP10)
             End Do !K
             ! For the definition of the densities below, see Eqs.(25)-(26) in CPC 167, 43 (2005)
#if(USE_OPENMP==1)
             OMPTaka(ib,ihil,it)    = OMPTaka(ib,ihil,it)    + TAKAIHIL
             OMPTro(ib,ihil,it)     = OMPTro(ib,ihil,it)     + TROIHIL
             OMPTdj(ib,ihil,it)     = OMPTdj(ib,ihil,it)     + TDJIHIL
             OMPTtau(ib,ihil,it)    = OMPTtau(ib,ihil,it)    + TTAUIHIL
             OMPTdro(ib,ihil,it)    = OMPTdro(ib,ihil,it)    + TDROIHIL
             OMPTSRFI(ib,ihil,it)   = OMPTSRFI(ib,ihil,it)   + TSRFIIHIL
             OMPTSFIR(ib,ihil,it)   = OMPTSFIR(ib,ihil,it)   + TSFIRIHIL
             OMPTSFIZ(ib,ihil,it)   = OMPTSFIZ(ib,ihil,it)   + TSFIZIHIL
             OMPTSZFI(ib,ihil,it)   = OMPTSZFI(ib,ihil,it)   + TSZFIIHIL
             OMPTNABLAR(ib,IHIL,IT) = OMPTNABLAR(ib,IHIL,IT) + TNABLARIHIL
             OMPTNABLAZ(ib,IHIL,IT) = OMPTNABLAZ(ib,IHIL,IT) + TNABLAZIHIL
#else
             OMPTaka(ihil,it)    = OMPTaka(ihil,it)    + TAKAIHIL  ! \kappa
             OMPTro(ihil,it)     = OMPTro(ihil,it)     + TROIHIL   ! \rho
             OMPTdj(ihil,it)     = OMPTdj(ihil,it)     + TDJIHIL   ! \div J
             OMPTtau(ihil,it)    = OMPTtau(ihil,it)    + TTAUIHIL  ! \tau
             OMPTdro(ihil,it)    = OMPTdro(ihil,it)    + TDROIHIL  ! \Delta\rho
             OMPTSRFI(ihil,it)   = OMPTSRFI(ihil,it)   + TSRFIIHIL ! J_{r\varphi}
             OMPTSFIR(ihil,it)   = OMPTSFIR(ihil,it)   + TSFIRIHIL ! J_{\varphi r}
             OMPTSFIZ(ihil,it)   = OMPTSFIZ(ihil,it)   + TSFIZIHIL ! J_{\varphi z}
             OMPTSZFI(ihil,it)   = OMPTSZFI(ihil,it)   + TSZFIIHIL ! J_{z\varphi}
             OMPTNABLAR(IHIL,IT) = OMPTNABLAR(IHIL,IT) + TNABLARIHIL
             OMPTNABLAZ(IHIL,IT) = OMPTNABLAZ(IHIL,IT) + TNABLAZIHIL
#endif
          End Do !ihil
       End If
    End Do !bb
#if(USE_OPENMP==1)
!$OMP End Parallel Do
!
!$OMP PARALLEL DO        &
!$OMP& DEFAULT(NONE)     &
!$OMP& SCHEDULE(DYNAMIC) &
!$OMP& SHARED(nghl,nb,AKA,RO,DJ,TAU,DRO,SRFI,SFIR,SFIZ,SZFI,NABLAR,NABLAZ, &
!$OMP&        OMPTaka,OMPTro,OMPTdj,OMPTtau,OMPTdro,OMPTSRFI,&
!$OMP&        OMPTSFIR,OMPTSFIZ,OMPTSZFI,OMPTNABLAR,OMPTNABLAZ)
    Do ihil = 1,nghl
       Do ib=1,nb
          AKA(ihil,1)    = AKA(ihil,1)    + OMPTaka(ib,ihil,1)
          RO(ihil,1)     = RO(ihil,1)     + OMPTro(ib,ihil,1)
          DJ(ihil,1)     = DJ(ihil,1)     + OMPTdj(ib,ihil,1)
          TAU(ihil,1)    = TAU(ihil,1)    + OMPTtau(ib,ihil,1)
          DRO(ihil,1)    = DRO(ihil,1)    + OMPTdro(ib,ihil,1)
          SRFI(ihil,1)   = SRFI(ihil,1)   + OMPTSRFI(ib,ihil,1)
          SFIR(ihil,1)   = SFIR(ihil,1)   + OMPTSFIR(ib,ihil,1)
          SFIZ(ihil,1)   = SFIZ(ihil,1)   + OMPTSFIZ(ib,ihil,1)
          SZFI(ihil,1)   = SZFI(ihil,1)   + OMPTSZFI(ib,ihil,1)
          NABLAR(ihil,1) = NABLAR(ihil,1) + OMPTNABLAR(ib,ihil,1)
          NABLAZ(ihil,1) = NABLAZ(ihil,1) + OMPTNABLAZ(ib,ihil,1)
          !
          AKA(ihil,2)    = AKA(ihil,2)    + OMPTaka(ib,ihil,2)
          RO(ihil,2)     = RO(ihil,2)     + OMPTro(ib,ihil,2)
          DJ(ihil,2)     = DJ(ihil,2)     + OMPTdj(ib,ihil,2)
          TAU(ihil,2)    = TAU(ihil,2)    + OMPTtau(ib,ihil,2)
          DRO(ihil,2)    = DRO(ihil,2)    + OMPTdro(ib,ihil,2)
          SRFI(ihil,2)   = SRFI(ihil,2)   + OMPTSRFI(ib,ihil,2)
          SFIR(ihil,2)   = SFIR(ihil,2)   + OMPTSFIR(ib,ihil,2)
          SFIZ(ihil,2)   = SFIZ(ihil,2)   + OMPTSFIZ(ib,ihil,2)
          SZFI(ihil,2)   = SZFI(ihil,2)   + OMPTSZFI(ib,ihil,2)
          NABLAR(ihil,2) = NABLAR(ihil,2) + OMPTNABLAR(ib,ihil,2)
          NABLAZ(ihil,2) = NABLAZ(ihil,2) + OMPTNABLAZ(ib,ihil,2)
       End Do
    End Do
!$OMP End Parallel Do
    Deallocate(OMPTaka,OMPTro,OMPTdj,OMPTtau,OMPTdro,OMPTSRFI,OMPTSFIR,OMPTSFIZ, &
               OMPTSZFI,OMPTNABLAR,OMPTNABLAZ)
#else
    Do it=1,2
       Do ihil = 1,nghl
          AKA(ihil,it)    = OMPTaka(ihil,it)
          RO(ihil,it)     = OMPTro(ihil,it)
          DJ(ihil,it)     = OMPTdj(ihil,it)
          TAU(ihil,it)    = OMPTtau(ihil,it)
          DRO(ihil,it)    = OMPTdro(ihil,it)
          SRFI(ihil,it)   = OMPTSRFI(ihil,it)
          SFIR(ihil,it)   = OMPTSFIR(ihil,it)
          SFIZ(ihil,it)   = OMPTSFIZ(ihil,it)
          SZFI(ihil,it)   = OMPTSZFI(ihil,it)
          NABLAR(ihil,it) = OMPTNABLAR(ihil,it)
          NABLAZ(ihil,it) = OMPTNABLAZ(ihil,it)
       End Do
    End Do
#endif
    If(iiter.Ge.2.And.si.Lt.epsi.and.iError_in_HO.eq.0) ro_normalization=ro ! Save density for later normalization
    Do it = 1,2
       TRO=>ro(:,it);         TTAU=>tau(:,it);       TDJ=>dj(:,it);     TDRO=>dro(:,it)
       TSZFI=>SZFI(:,it);     TSFIZ=>SFIZ(:,it);     TSRFI=>SRFI(:,it); TSFIR=>SFIR(:,it)
       TNABLAR=>NABLAR(:,it); TNABLAZ=>NABLAZ(:,it); TAKA=>aka(:,it)
       s=two*Sum(tro); sd=four*Sum(tdro); drhoi(it)=sd
       Sumnz(it)=Abs(s-tz(it))
       varmas=varmas+s
       DNFactor(it)=tz(it)/s
       !----------------------------------------------------
       ! REMOVES INT.WEIGHTS AND MULTIPLIES BY THE JACOBIAN
       !----------------------------------------------------
       piu=two*tz(it)/s
       WGT=wdcori
       Call dscal(NGHL,piu,WGT,1)
       Tro=Tro*WGT; Ttau=Ttau*WGT; Taka=Half*Taka*WGT;
       TSRFI=TSRFI*WGT; TSFIR=TSFIR*WGT;
       TSFIZ=TSFIZ*WGT; TSZFI=TSZFI*WGT;
       Call dscal(NGHL,two,WGT,1)
       Tdro=Tdro*WGT; Tdj=Tdj*WGT
       TNABLAR=TNABLAR*WGT; TNABLAZ=TNABLAZ*WGT
       !
    End Do !it
    DNFactor(3)=DNFactor(1)+DNFactor(2)
    !
    If (debug_solver.Ge.1) Call get_CPU_time('densit',1)
    !----------------------------------------------------
    ! COULOMB AND HARTREE FIELDS
    !----------------------------------------------------
    If(icou.Gt.0) Then
       If(nleg.Lt.0) Then
          Call coulom1
       Else
          Call coulom
       End If
    Else
       cou=zero
    End If
    !Call HartreeDir
  End Subroutine DENSIT
  !=======================================================================
  !> When expressed in terms of isoscalar and isovector densities, the total energy takes the most general form:
  !>
  !>  \f{eqnarray}{
  !>      E & = & \frac{1}{2}hb0(\tau_0 + \tau_1) + \frac{1}{2}hb0(\tau_0 - \tau_1) \\
  !>
  !>      & + & U^{\rho\tau}_{0,0}\rho_0 \tau_0 + U^{\rho\tau}_{1,0}\rho_1 \tau_1
  !>          + U^{\rho\tau}_{2,0}\rho_0 \tau_1 + U^{\rho\tau}_{3,0}\rho_1 \tau_0 \\
  !>
  !>      & + & U^{\rho\rho}_{0,0}\rho_0^2 + U^{\rho\rho}_{1,0}\rho_1^2 + \big[ U^{\rho\rho}_{2,0} + U^{\rho\rho}_{3,0} \big]\rho_0 \rho_1 \\
  !>
  !>      & + & U^{\rho\Delta\rho}_{0,0}\rho_0\Delta\rho_0 + U^{\rho\Delta\rho}_{1,0}\rho_1\Delta\rho_1
  !>          + U^{\rho\Delta\rho}_{2,0}\rho_0\Delta\rho_1 + U^{\rho\Delta\rho}_{3,0}\rho_1\Delta\rho_0 \\
  !>
  !>      & + & U^{\rho\nabla J}_{0,0}\rho_0\nabla\cdot J_0 + U^{\rho\nabla J}_{1,0}\rho_1\nabla\cdot J_1
  !>          + U^{\rho\nabla J}_{2,0}\rho_0\nabla\cdot J_1 + U^{\rho\nabla J}_{3,0}\rho_1\nabla\cdot J_0 \\
  !>
  !>      & + & \color{red}{
  !>            U^{\nabla\rho\nabla\rho}_{0,0}\big[ \nabla_{\rho}\rho_{0}\nabla_{\rho}\rho_{0} + \nabla_{z}\rho_{0}\nabla_{z}\rho_{0} \big]
  !>          + U^{\nabla\rho\nabla\rho}_{1,0}\big[ \nabla_{\rho}\rho_{1}\nabla_{\rho}\rho_{1} + \nabla_{z}\rho_{1}\nabla_{z}\rho_{1} \big]
  !>          +\big[ U^{\nabla\rho\nabla\rho}_{3,0}+U^{\nabla\rho\nabla\rho}_{2,0} \big]
  !>           \big[ \nabla_{\rho}\rho_{0}\nabla_{\rho}\rho_{1} + \nabla_{z}\rho_{0}\nabla_{z}\rho_{1} \big] } \\
  !>
  !>      & + & \color{red}{
  !>            U^{J\nabla\rho}_{0,0}\big[ \nabla_{\rho}\rho_{0}(J_{\varphi z,0}-J_{z\varphi,0}) - \nabla_{z}\rho_{0}(J_{\varphi\rho,0}-J_{\rho\varphi,0}) \big]
  !>          + U^{J\nabla\rho}_{1,0}\big[ \nabla_{\rho}\rho_{1}(J_{\varphi z,1}-J_{z\varphi,1}) - \nabla_{z}\rho_{1}(J_{\varphi\rho,1}-J_{\rho\varphi,1}) \big]
  !>          + U^{J\nabla\rho}_{2,0}\big[ \nabla_{\rho}\rho_{1}(J_{\varphi z,0}-J_{z\varphi,0}) - \nabla_{z}\rho_{1}(J_{\varphi\rho,0}-J_{\rho\varphi,0}) \big]
  !>          + U^{J\nabla\rho}_{3,0}\big[ \nabla_{\rho}\rho_{0}(J_{\varphi z,1}-J_{z\varphi,1}) - \nabla_{z}\rho_{0}(J_{\varphi\rho,1}-J_{\rho\varphi,1}) \big] } \\
  !>
  !>      & + & U^{JJ}_{0,0}J^2_0 + U^{JJ}_{1,0}J^2_1
  !>         + \big[ U^{JJ}_{3,0} + U^{JJ}_{2,0} \big]
  !>           \big[ J_{\varphi z,0}J_{\varphi z,1} + J_{\varphi\rho,0}J_{\varphi\rho,1} + J_{z\varphi, 0}J_{z\varphi, 1} + J_{\rho\varphi,0}J_{\rho\varphi,1} \big] \\
  !>     & + & 2U^{JJ}_{0,0}\big[ J_{\varphi z,0}J_{z\varphi, 0} + J_{\varphi\rho,0}J_{\rho\varphi,0} \big]
  !>         + 2U^{JJ}_{1,0}\big[ J_{\varphi z,1}J_{z\varphi, 1} + J_{\varphi\rho,1}J_{\rho\varphi,1} \big]
  !>       + \big[ U^{JJ}_{3,0} + U^{JJ}_{2,0} \big] \big[ J_{\varphi z,0}J_{z\varphi, 1} + J_{\varphi\rho,0}J_{\rho\varphi,1}
  !>                                                     + J_{z\varphi, 0}J_{\varphi z,1} + J_{\rho\varphi,0}J_{\varphi\rho,1} \big]
  !> \f}
  !>
  !> The terms marked in red vanish identically when the coupling constants \f$ U_{t,0}^{uu'} \f$  are not density-dependent, since in
  !> this case, integration by part of these terms simply give an extra contribution to the \f$ \rho\Delta\rho \f$ and \f$ \rho\nabla J\f$
  !> terms.
  !=======================================================================
  Subroutine field
    USE UNEDF, Only: is_NEDF
    Implicit None
    Integer(ipr) :: iw,it,ita,ihli,lambda,icons
    Real(pr) :: ra,ra2,rs,rsa,rsa0,z,rrr,rear_pair
    Real(pr) :: rt,rt1,ds,da,dt,dt1,tts,tta,tt,tt1,djs,dja,djt,djt1
    Real(pr) :: rsa0A,rsa0A1,V0V1,v01a,rns,rps,rsa1,rsa12,rsa10
    Real(pr) :: rsa0An,rsa0An1,rsa0As,rsa0As1
    Real(pr) :: RHO_0,RHO_1,TAU_0,TAU_1,DRHO_0,DRHO_1,DJ_0,DJ_1
    Real(pr) :: SZFIN,SFIZN,SRFIN,SFIRN,SZFIP,SFIZP,SRFIP,SFIRP
    Real(pr) :: SZFI_0,SFIZ_0,SRFI_0,SFIR_0,SZFI_1,SFIZ_1,SRFI_1,SFIR_1
    Real(pr) :: SNABLARN,SNABLAZN,SNABLARP,SNABLAZP
    Real(pr) :: SNABLAR_0,SNABLAZ_0,SNABLAR_1,SNABLAZ_1
    Real(pr) :: J2_0,J2_1,JabJba_0,JabJba_1
    Real(pr) :: cx,x,fac_n,fac_p,gr,Pi,ec,kc,kf,lc,a,b
    Real(pr), Dimension(0:8) :: Qval
    Real(pr),Dimension(2) :: pUr,pUt,pUNr,pUNz,pUDr,pUDj,pUFIZ,pUZFI,pUFIR,pURFI
    Real(pr),Dimension(2) :: tUr,tUt,tUNr,tUNz,tUDr,tUDj,tUFIZ,tUZFI,tUFIR,tURFI
    Real(pr), Save :: ALAMBDA=0.0_pr,AEPSI=1.0_pr,CSPR=1.0_pr
    !
    If (debug_solver.Ge.1) Call get_CPU_time('field',0)
    !
    Pi = 4.0_pr*Atan(1.0_pr)
    !
    ! fields
    Do ihli=1,nghl
       ! Scalar densities and components of div J
       RHO_0  =   ro(ihli,1) + ro(ihli,2) ; RHO_1  =   ro(ihli,1)  -ro(ihli,2)
       TAU_0  =  tau(ihli,1) +tau(ihli,2) ; TAU_1  =  tau(ihli,1) -tau(ihli,2)
       DRHO_0 =  dro(ihli,1) +dro(ihli,2) ; DRHO_1 =  dro(ihli,1) -dro(ihli,2)
       DJ_0   =   dj(ihli,1)  +dj(ihli,2) ; DJ_1   =   dj(ihli,1)  -dj(ihli,2)
       SFIZ_0 = SFIZ(ihli,1)+SFIZ(ihli,2) ; SFIZ_1 = SFIZ(ihli,1)-SFIZ(ihli,2)
       SFIR_0 = SFIR(ihli,1)+SFIR(ihli,2) ; SFIR_1 = SFIR(ihli,1)-SFIR(ihli,2)
       SZFI_0 = SZFI(ihli,1)+SZFI(ihli,2) ; SZFI_1 = SZFI(ihli,1)-SZFI(ihli,2)
       SRFI_0 = SRFI(ihli,1)+SRFI(ihli,2) ; SRFI_1 = SRFI(ihli,1)-SRFI(ihli,2)
       ! Components r and z of  \nabla\rho
       SNABLAR_0 = NABLAR(ihli,1) + NABLAR(ihli,2)
       SNABLAR_1 = NABLAR(ihli,1) - NABLAR(ihli,2)
       SNABLAZ_0 = NABLAZ(ihli,1) + NABLAZ(ihli,2)
       SNABLAZ_1 = NABLAZ(ihli,1) - NABLAZ(ihli,2)
       ! Various components of the tensor density
       J2_0 = SFIZ_0**2 + SFIR_0**2 + SZFI_0**2 + SRFI_0**2
       J2_1 = SFIZ_1**2 + SFIR_1**2 + SZFI_1**2 + SRFI_1**2
       JabJba_0 = 2*(SFIZ_0*SZFI_0 + SFIR_0*SRFI_0)
       JabJba_1 = 2*(SFIZ_1*SZFI_1 + SFIR_1*SRFI_1)
       !
       tUr=zero ; tUDr=zero ; tUNr=zero ; tUNz=zero
       tUt=zero ; tUDj=zero ; tUFIZ=zero ; tUZFI=zero
       tUFIR=zero ; tURFI=zero ;
       !
       Call calculate_U_parameters(RHO_0,RHO_1,TAU_0,TAU_1,DRHO_0,DRHO_1,(SNABLAR_0**2+SNABLAZ_0**2),(SNABLAR_1**2+SNABLAZ_1**2) )
       !
       ! tUr(1)=dE/d RHO_0;       tUr(2)=dE/d RHO_1
       ! tUt(1)=dE/d TAU_0;       tUt(2)=dE/d TAU_1
       ! tUDr(1)=dE/d DeltaRHO_0; tUDr(2)=dE/d DeltaRHO_1
       ! and so on ...
       !
       ! Contributions in the case 'u' depends on RHO_0
       ! \rho^2
       tUr(1)  =  tUr(1) + two*Urhorho(0,0)*RHO_0 + Urhorho(0,1)*RHO_0*RHO_0 + Urhorho(1,1)*RHO_1*RHO_1  &
                             +(Urhorho(3,0)+Urhorho(2,0))*RHO_1 + (Urhorho(3,1)+ Urhorho(2,1))*RHO_0*RHO_1
       tUr(2)  =  tUr(2) + two*Urhorho(1,0)*RHO_1+Urhorho(0,2)*RHO_0*RHO_0 + Urhorho(1,2)*RHO_1*RHO_1  &
                             +(Urhorho(3,0)+Urhorho(2,0))*RHO_0 + (Urhorho(3,2)+Urhorho(2,2))*RHO_0*RHO_1
       tUr(1)  =  tUr(1) + vDHartree(ihli,1)
       tUr(2)  =  tUr(2) + vDHartree(ihli,2)
       ! \rho\tau
       tUr(1)  =  tUr(1) + Urhotau(0,0)*TAU_0 + Urhotau(0,1)*TAU_0*RHO_0 + Urhotau(1,1)*TAU_1*RHO_1  &
                         + Urhotau(2,0)*TAU_1 + Urhotau(2,1)*RHO_0*TAU_1 + Urhotau(3,1)*RHO_1*TAU_0
       tUt(1)  =  tUt(1) + Urhotau(0,0)*RHO_0 + Urhotau(3,0)*RHO_1
       tUr(2)  =  tUr(2) + Urhotau(1,0)*TAU_1 + Urhotau(1,2)*TAU_1*RHO_1 + Urhotau(0,2)*TAU_0*RHO_0  &
                         + Urhotau(3,0)*TAU_0 + Urhotau(3,2)*RHO_1*TAU_0 + Urhotau(2,2)*RHO_0*TAU_1
       tUt(2)  =  tUt(2) + Urhotau(1,0)*RHO_1 + Urhotau(2,0)*RHO_0
       ! \rho\Delta\rho
       tUr(1)  =  tUr(1) + UrhoDrho(0,0)*DRHO_0 + UrhoDrho(0,1)*RHO_0*DRHO_0 + UrhoDrho(1,1)*RHO_1*DRHO_1  &
                         + UrhoDrho(2,0)*DRHO_1 + UrhoDrho(2,1)*RHO_0*DRHO_1 + UrhoDrho(3,1)*RHO_1*DRHO_0
       tUDr(1) = tUDr(1) + UrhoDrho(0,0)*RHO_0  + UrhoDrho(3,0)*RHO_1
       tUr(2 ) = tUr(2)  + UrhoDrho(1,0)*DRHO_1 + UrhoDrho(1,2)*RHO_1*DRHO_1 + UrhoDrho(0,2)*RHO_0*DRHO_0  &
                         + UrhoDrho(3,0)*DRHO_0 + UrhoDrho(3,2)*RHO_1*DRHO_0 + UrhoDrho(2,2)*RHO_0*DRHO_1
       tUDr(2) = tUDr(2) + UrhoDrho(1,0)*RHO_1  + UrhoDrho(2,0)*RHO_0
       ! (\nabla\rho)^2 (cancels for density-independent coupling constants)
       tUr(1)  = tUr(1)  +  Unablarho(0,1)*(SNABLAR_0**2+SNABLAZ_0**2) + Unablarho(1,1)*(SNABLAR_1**2+SNABLAZ_1**2) &
                         + (Unablarho(2,1)+Unablarho(3,1))*(SNABLAR_0*SNABLAR_1+SNABLAZ_0*SNABLAZ_1)
       tUNr(1) = tUNr(1) + two*Unablarho(0,0)*SNABLAR_0 + (Unablarho(2,0) + Unablarho(3,0))*SNABLAR_1
       tUNz(1) = tUNz(1) + two*Unablarho(0,0)*SNABLAZ_0 + (Unablarho(2,0) + Unablarho(3,0))*SNABLAZ_1
       tUr(2)  = tUr(2)  +  Unablarho(0,2)*(SNABLAR_0**2+SNABLAZ_0**2) + Unablarho(1,2)*(SNABLAR_1**2+SNABLAZ_1**2) &
                         + (Unablarho(2,2)+Unablarho(3,2))*(SNABLAR_0*SNABLAR_1+SNABLAZ_0*SNABLAZ_1)
       tUNr(2) = tUNr(2) + two*Unablarho(1,0)*SNABLAR_1 + (Unablarho(2,0) + Unablarho(3,0))*SNABLAR_0
       tUNz(2) = tUNz(2) + two*Unablarho(1,0)*SNABLAZ_1 + (Unablarho(2,0) + Unablarho(3,0))*SNABLAZ_0
       ! \rho\nabla J (spin-orbit)
       tUr(1)  = tUr(1)  + UrhonablaJ(0,0)*DJ_0+UrhonablaJ(0,1)*DJ_0*RHO_0+UrhonablaJ(1,1)*DJ_1*RHO_1  &
                         + UrhonablaJ(2,0)*DJ_1+UrhonablaJ(2,1)*RHO_0*DJ_1+UrhonablaJ(3,1)*RHO_1*DJ_0
       tUDj(1) = tUDj(1) + UrhonablaJ(0,0)*RHO_0+UrhonablaJ(3,0)*RHO_1
       tUr(2)  = tUr(2)  + UrhonablaJ(1,0)*DJ_1+UrhonablaJ(1,2)*DJ_1*RHO_1+UrhonablaJ(0,2)*DJ_0*RHO_0  &
                         + UrhonablaJ(3,0)*DJ_0+UrhonablaJ(3,2)*RHO_1*DJ_0+UrhonablaJ(2,2)*RHO_0*DJ_1
       tUDj(2) = tUDj(2) + UrhonablaJ(1,0)*RHO_1+UrhonablaJ(2,0)*RHO_0
       ! J\nabla\rho (cancels for density-independent coupling constants)
       tUr(1)   = tUr(1)   + UJnablarho(0,1)*(SNABLAR_0*(SFIZ_0-SZFI_0)-SNABLAZ_0*(SFIR_0-SRFI_0))
       tUr(1)   = tUr(1)   + UJnablarho(1,1)*(SNABLAR_1*(SFIZ_1-SZFI_1)-SNABLAZ_1*(SFIR_1-SRFI_1))
       tUr(1)   = tUr(1)   + UJnablarho(2,1)*(SNABLAR_1*(SFIZ_0-SZFI_0)-SNABLAZ_1*(SFIR_0-SRFI_0))
       tUr(1)   = tUr(1)   + UJnablarho(3,1)*(SNABLAR_0*(SFIZ_1-SZFI_1)-SNABLAZ_0*(SFIR_1-SRFI_1))
       tUr(2)   = tUr(2)   + UJnablarho(0,2)*(SNABLAR_0*(SFIZ_0-SZFI_0)-SNABLAZ_0*(SFIR_0-SRFI_0))
       tUr(2)   = tUr(2)   + UJnablarho(1,2)*(SNABLAR_1*(SFIZ_1-SZFI_1)-SNABLAZ_1*(SFIR_1-SRFI_1))
       tUr(2)   = tUr(2)   + UJnablarho(2,2)*(SNABLAR_1*(SFIZ_0-SZFI_0)-SNABLAZ_1*(SFIR_0-SRFI_0))
       tUr(2)   = tUr(2)   + UJnablarho(3,2)*(SNABLAR_0*(SFIZ_1-SZFI_1)-SNABLAZ_0*(SFIR_1-SRFI_1))
       tUNr(1)  = tUNr(1)  + UJnablarho(0,0)*(SFIZ_0-SZFI_0)
       tUNr(2)  = tUNr(2)  + UJnablarho(1,0)*(SFIZ_1-SZFI_1)
       tUNz(1)  = tUNz(1)  - UJnablarho(0,0)*(SFIR_0-SRFI_0)
       tUNz(2)  = tUNz(2)  - UJnablarho(1,0)*(SFIR_1-SRFI_1)
       tUFIZ(1) = tUFIZ(1) + UJnablarho(0,0)*SNABLAR_0*half
       tUFIZ(2) = tUFIZ(2) + UJnablarho(1,0)*SNABLAR_1*half
       tUZFI(1) = tUZFI(1) - UJnablarho(0,0)*SNABLAR_0*half
       tUZFI(2) = tUZFI(2) - UJnablarho(1,0)*SNABLAR_1*half
       tURFI(1) = tURFI(1) + UJnablarho(0,0)*SNABLAZ_0*half
       tURFI(2) = tURFI(2) + UJnablarho(1,0)*SNABLAZ_1*half
       tUFIR(1) = tUFIR(1) - UJnablarho(0,0)*SNABLAZ_0*half
       tUFIR(2) = tUFIR(2) - UJnablarho(1,0)*SNABLAZ_1*half
       ! Tensor terms: Diagonal terms J^2
       tUr(1) = tUr(1) + UJJ(0,1)*J2_0 + UJJ(1,1)*J2_1 &
                       +(UJJ(3,1)+UJJ(2,1))*(SFIZ_0*SFIZ_1+SFIR_0*SFIR_1+SZFI_0*SZFI_1+SRFI_0*SRFI_1)
       tUr(2) = tUr(2) + UJJ(0,2)*J2_0 + UJJ(1,2)*J2_1 &
                       +(UJJ(3,2)+UJJ(2,2))*(SFIZ_0*SFIZ_1+SFIR_0*SFIR_1+SZFI_0*SZFI_1+SRFI_0*SRFI_1)
       tUFIZ(1) = tUFIZ(1) + UJJ(0,0)*SFIZ_0 + half*(UJJ(3,0)+UJJ(2,0))*SFIZ_1
       tUFIR(1) = tUFIR(1) + UJJ(0,0)*SFIR_0 + half*(UJJ(3,0)+UJJ(2,0))*SFIR_1
       tUZFI(1) = tUZFI(1) + UJJ(0,0)*SZFI_0 + half*(UJJ(3,0)+UJJ(2,0))*SZFI_1
       tURFI(1) = tURFI(1) + UJJ(0,0)*SRFI_0 + half*(UJJ(3,0)+UJJ(2,0))*SRFI_1
       tUFIZ(2) = tUFIZ(2) + UJJ(1,0)*SFIZ_1 + half*(UJJ(3,0)+UJJ(2,0))*SFIZ_0
       tUFIR(2) = tUFIR(2) + UJJ(1,0)*SFIR_1 + half*(UJJ(3,0)+UJJ(2,0))*SFIR_0
       tUZFI(2) = tUZFI(2) + UJJ(1,0)*SZFI_1 + half*(UJJ(3,0)+UJJ(2,0))*SZFI_0
       tURFI(2) = tURFI(2) + UJJ(1,0)*SRFI_1 + half*(UJJ(3,0)+UJJ(2,0))*SRFI_0
       ! Tensor terms: Cross terms Jab.Jba
       tUr(1) = tUr(1) + UJabJba(0,1)*JabJba_0 + UJabJba(1,1)*JabJba_1 &
                       +(UJabJba(3,1)+UJabJba(2,1))*(SFIZ_0*SZFI_1+SFIR_0*SRFI_1+SZFI_0*SFIZ_1+SRFI_0*SFIR_1)
       tUr(2) = tUr(2) + UJabJba(0,2)*JabJba_0 + UJabJba(1,2)*JabJba_1 &
                       +(UJabJba(3,2)+UJabJba(2,2))*(SFIZ_0*SZFI_1+SFIR_0*SRFI_1+SZFI_0*SFIZ_1+SRFI_0*SFIR_1)
       tUFIZ(1) = tUFIZ(1) + UJabJba(0,0)*SZFI_0 + half*(UJabJba(3,0)+UJabJba(2,0))*SZFI_1
       tUFIR(1) = tUFIR(1) + UJabJba(0,0)*SRFI_0 + half*(UJabJba(3,0)+UJabJba(2,0))*SRFI_1
       tUZFI(1) = tUZFI(1) + UJabJba(0,0)*SFIZ_0 + half*(UJabJba(3,0)+UJabJba(2,0))*SFIZ_1
       tURFI(1) = tURFI(1) + UJabJba(0,0)*SFIR_0 + half*(UJabJba(3,0)+UJabJba(2,0))*SFIR_1
       tUFIZ(2) = tUFIZ(2) + UJabJba(1,0)*SZFI_1 + half*(UJabJba(3,0)+UJabJba(2,0))*SZFI_0
       tUFIR(2) = tUFIR(2) + UJabJba(1,0)*SRFI_1 + half*(UJabJba(3,0)+UJabJba(2,0))*SRFI_0
       tUZFI(2) = tUZFI(2) + UJabJba(1,0)*SFIZ_1 + half*(UJabJba(3,0)+UJabJba(2,0))*SFIZ_0
       tURFI(2) = tURFI(2) + UJabJba(1,0)*SFIR_1 + half*(UJabJba(3,0)+UJabJba(2,0))*SFIR_0

       tUr(1)=tUr(1)+UFnonstdr(0)                                        !! other amplitudes
       tUr(2)=tUr(2)+UFnonstdr(1)

       !!  External Field
       !!
       !!tUr(1)=tUr(1)+Vexternal(0,zero,fl(ihli),fh(ihli))
       !!tUr(2)=tUr(2)+Vexternal(1,zero,fl(ihli),fh(ihli))

       ! Contributions in the case 'u' depends on TAU_0
       tUt(1)=tUt(1)+Urhotau(0,6)*RHO_0*TAU_0  &
                    +Urhotau(1,6)*RHO_1*TAU_1+Urhotau(2,6)*RHO_0*TAU_1  &
                    +Urhotau(3,6)*RHO_1*TAU_0+Urhorho(0,6)*RHO_0**2  &
                    +Urhorho(1,6)*RHO_1**2+(Urhorho(2,6)+Urhorho(3,6))*RHO_0*RHO_1  &
                    +UrhoDrho(0,6)*RHO_0*DRHO_0+UrhoDrho(1,6)*RHO_1*DRHO_1  &
                    +UrhoDrho(2,6)*RHO_0*DRHO_1+UrhoDrho(3,6)*RHO_1*DRHO_0  &
                    +Unablarho(0,6)*(SNABLAR_0*SNABLAR_0+SNABLAZ_0*SNABLAZ_0)  &
                    +Unablarho(1,6)*(SNABLAR_1*SNABLAR_1+SNABLAZ_1*SNABLAZ_1)  &
                    +(Unablarho(2,6)+Unablarho(3,6))*(SNABLAR_0*SNABLAR_1+SNABLAZ_0*SNABLAZ_1)  &
                    +UrhonablaJ(0,6)*RHO_0*DJ_0+UrhonablaJ(1,6)*RHO_1*DJ_1  &
                    +UrhonablaJ(2,6)*RHO_0*DJ_1+UrhonablaJ(3,6)*RHO_1*DJ_0  &
                    +UJnablarho(0,6)*(SNABLAR_0*(SFIZ_0-SZFI_0)-SNABLAZ_0*(SFIR_0-SRFI_0))  &
                    +UJnablarho(1,6)*(SNABLAR_1*(SFIZ_1-SZFI_1)-SNABLAZ_1*(SFIR_1-SRFI_1))  &
                    +UJnablarho(2,6)*(SNABLAR_1*(SFIZ_0-SZFI_0)-SNABLAZ_1*(SFIR_0-SRFI_0))  &
                    +UJnablarho(3,6)*(SNABLAR_0*(SFIZ_1-SZFI_1)-SNABLAZ_0*(SFIR_1-SRFI_1))
       tUt(1)=tUt(1)+UJJ(0,6)*J2_0+UJJ(1,6)*J2_1  &
                   +(UJJ(2,6)+UJJ(3,6))*(SFIZ_0*SFIZ_1+SFIR_0*SFIR_1+SZFI_0*SZFI_1+SRFI_0*SRFI_1)
       tUt(1)=tUt(1)+UJabJba(0,6)*JabJba_0+UJJ(1,6)*JabJba_1  &
                   +(UJabJba(2,6)+UJabJba(3,6))*(SFIZ_0*SZFI_1+SFIR_0*SRFI_1+SZFI_0*SFIZ_1+SRFI_0*SFIR_1)
       !
       ! Contributions in the case 'u' depends on DeltaRHO_0
       tUDr(1)=tUDr(1)+Urhotau(0,7)*RHO_0*TAU_0  &
                      +Urhotau(1,7)*RHO_1*TAU_1+Urhotau(2,7)*RHO_0*TAU_1  &
                      +Urhotau(3,7)*RHO_1*TAU_0+Urhorho(0,7)*RHO_0**2  &
                      +Urhorho(1,7)*RHO_1**2+(Urhorho(2,7)+Urhorho(3,7))*RHO_0*RHO_1  &
                      +UrhoDrho(0,7)*RHO_0*DRHO_0+UrhoDrho(1,7)*RHO_1*DRHO_1  &
                      +UrhoDrho(2,7)*RHO_0*DRHO_1+UrhoDrho(3,7)*RHO_1*DRHO_0  &
                      +Unablarho(0,7)*(SNABLAR_0*SNABLAR_0+SNABLAZ_0*SNABLAZ_0)  &
                      +Unablarho(1,7)*(SNABLAR_1*SNABLAR_1+SNABLAZ_1*SNABLAZ_1)  &
                      +(Unablarho(2,7)+Unablarho(3,7))*(SNABLAR_0*SNABLAR_1+SNABLAZ_0*SNABLAZ_1)  &
                      +UrhonablaJ(0,7)*RHO_0*DJ_0+UrhonablaJ(1,7)*RHO_1*DJ_1  &
                      +UrhonablaJ(2,7)*RHO_0*DJ_1+UrhonablaJ(3,7)*RHO_1*DJ_0  &
                      +UJnablarho(0,7)*(SNABLAR_0*(SFIZ_0-SZFI_0)-SNABLAZ_0*(SFIR_0-SRFI_0))  &
                      +UJnablarho(1,7)*(SNABLAR_1*(SFIZ_1-SZFI_1)-SNABLAZ_1*(SFIR_1-SRFI_1))  &
                      +UJnablarho(2,7)*(SNABLAR_1*(SFIZ_0-SZFI_0)-SNABLAZ_1*(SFIR_0-SRFI_0))  &
                      +UJnablarho(3,7)*(SNABLAR_0*(SFIZ_1-SZFI_1)-SNABLAZ_0*(SFIR_1-SRFI_1))
       tUDr(1)=tUDr(1)+UJJ(0,7)*J2_0+UJJ(1,7)*J2_1  &
                     +(UJJ(2,7)+UJJ(3,7))*(SFIZ_0*SFIZ_1+SFIR_0*SFIR_1+SZFI_0*SZFI_1+SRFI_0*SRFI_1)
       tUDr(1)=tUDr(1)+UJabJba(0,7)*JabJba_0+UJJ(1,7)*JabJba_1  &
                     +(UJabJba(2,7)+UJabJba(3,7))*(SFIZ_0*SZFI_1+SFIR_0*SRFI_1+SZFI_0*SFIZ_1+SRFI_0*SFIR_1)

       ! Back to proton-neutron representation
       pUr(1)  =tUr(1)+tUr(2);            pUr(2)  =tUr(1)  -tUr(2)
       pUt(1)  =tUt(1)+tUt(2)+hb0n*facECM;pUt(2)  =tUt(1)  -tUt(2)+hb0p*facECM
       pUDr(1) =tUDr(1)+tUDr(2);          pUDr(2) =tUDr(1) -tUDr(2)
       pUNr(1) =tUNr(1)+tUNr(2);          pUNr(2) =tUNr(1) -tUNr(2)
       pUNz(1) =tUNz(1)+tUNz(2);          pUNz(2) =tUNz(1) -tUNz(2)
       pUDj(1) =tUDj(1)+tUDj(2);          pUDj(2) =tUDj(1) -tUDj(2)
       pUFIZ(1)=tUFIZ(1)+tUFIZ(2);        pUFIZ(2)=tUFIZ(1)-tUFIZ(2)
       pUZFI(1)=tUZFI(1)+tUZFI(2);        pUZFI(2)=tUZFI(1)-tUZFI(2)
       pUFIR(1)=tUFIR(1)+tUFIR(2);        pUFIR(2)=tUFIR(1)-tUFIR(2)
       pURFI(1)=tURFI(1)+tURFI(2);        pURFI(2)=tURFI(1)-tURFI(2)
       !
       Do it=itmin,itmax   !! loop over n  & p
          ita=3-it
          ! constraining potential
          If (numberCons.Gt.0) Then
              z=fh(ihli); rrr=fl(ihli)**2
              Call moments_valueMesh(z,rrr,Qval)
              do icons=1,numberCons
                 lambda=multLambda(icons)
                 If(lambda.Ge.1) Then
                    pUr(it)= pUr(it) - multLag(lambda)*Qval(lambda)
                 End If
                 If(lambda.Eq.0) Then
                    pUr(it)= pUr(it) - neckLag*Exp(-((z-Z_NECK*bz)/AN_VAL)**2)
                 End If
              end do
          End If
          ! Coulomb
          If(it.Eq.2) Then
             If(icou.Ge.1) pUr(it)=pUr(it)+cou(ihli)
             If(icou.Eq.2.Or.icou.Eq.-4) pUr(it)=pUr(it)+CExPar*coex*ro(ihli,it)**p13 ! Slater approximation
          End If
          rsa0=(ro(ihli,it)+ro(ihli,ita))/rho_c
          ! Pairing regularization
          If(pairing_regularization) Then
             !
             ec = pwi - ala(it) ! cut-off energy for s.p. (for comparisons with HFBRAD)
             gr = (CpV0(it-1)*(ONE-rsa0*CpV1(it-1))) ! original pairing strength
             !
             If(it.Eq.1) Then
                If(is_NEDF .Or. skyrme.Eq.'SLY4mod') Then
                   MEFFn(ihli) = 0.5_pr*(hb0n+hb0p)
                Else
                   MEFFn(ihli) = hb0n + (Ctau(0)-Ctau(1))*RHO_0 + 2.0_pr*Ctau(1)*ro(ihli,it)
                End If
                fac_n=Sqrt(one/MEFFn(ihli))
                If(pUr(1) + ec - ala(it) .lt. 0.0_pr) Then
                   kc = fac_n*Sqrt(ala(it)+ec-pUr(1))
                   lc = fac_n*Sqrt(ala(it)-ec-pUr(1))
                   kf = fac_n*Sqrt(ala(it)-pUr(1))
                   a = 0.25_pr*(kc/MEFFn(ihli)/Pi**2)*(one-0.5_pr*kf/kc*Log((kc+kf)/(kc-kf)))
                   b = 0.25_pr*(lc/MEFFn(ihli)/Pi**2)*(one-0.5_pr*kf/lc*Log((lc+kf)/(kf-lc)))
                   geff_inv(ihli,it) = one/gr - a -b
                Else
                   If(ala(it) .Gt. pUr(1)) Then
                      kc = fac_n*Sqrt(ala(it)+ec-pUr(1))
                      kf = fac_n*Sqrt(ala(it)-pUr(1))
                      geff_inv(ihli,it) = one/gr - 0.25_pr*(kc/MEFFn(ihli)/Pi**2)*(one-0.5_pr*kf/kc   &
                                                                                *Log((kc+kf)/(kc-kf)))
                   Else
                      If(pUr(1) - ec - ala(it) .lt. 0.0_pr) Then
                         kc = fac_n*Sqrt(ala(it)+ec-pUr(1))
                         kf = fac_n*Sqrt(pUr(1)-ala(it))
                         geff_inv(ihli,it) = one/gr - 0.25_pr*(kc/MEFFn(ihli)/Pi**2)*(one+kf/kc &
                                                                                    *Atan(kf/kc))
                      Else
                         geff_inv(ihli,it) = one/gr
                      End If
                   End If
                End If
                ! contribution to Delta
                dvn(ihli)=aka(ihli,it)/geff_inv(ihli,it)
                ! contribution to rearrangement
                !rear_pair =
             Else
                If(is_NEDF .Or. skyrme.Eq.'SLY4mod') Then
                   MEFFp(ihli) = 0.5_pr*(hb0n+hb0p)
                Else
                   MEFFp(ihli) = hb0p + (Ctau(0)-Ctau(1))*RHO_0 + 2.0_pr*Ctau(1)*ro(ihli,it)
                End If
                fac_p=Sqrt(one/MEFFp(ihli))
                If(pUr(2) + ec - ala(it) .lt. 0.0_pr) Then
                   kc = fac_p*Sqrt(ala(it)+ec-pUr(2))
                   lc = fac_p*Sqrt(ala(it)-ec-pUr(2))
                   kf = fac_p*Sqrt(ala(it)-pUr(2))
                   a = 0.25_pr*(kc/MEFFp(ihli)/Pi**2)*(one-0.5_pr*kf/kc*Log((kc+kf)/(kc-kf)))
                   b = 0.25_pr*(lc/MEFFp(ihli)/Pi**2)*(one-0.5_pr*kf/lc*Log((lc+kf)/(kf-lc)))
                   geff_inv(ihli,it) = one/gr - a -b
                Else
                   If(ala(it).Gt.pUr(2)) Then
                      kc = fac_p*Sqrt(ala(it)+ec-pUr(2))
                      kf = fac_p*Sqrt(ala(it)-pUr(2))
                      geff_inv(ihli,it) = one/gr - 0.25_pr*(kc/MEFFp(ihli)/Pi**2)*(one-0.5_pr*kf/kc &
                                                                                *Log((kc+kf)/(kc-kf)))
                   Else
                      If(pUr(2) - ec - ala(it) .lt. 0.0_pr) Then
                         kc = fac_p*Sqrt(ala(it)+ec-pUr(2))
                         kf = fac_p*Sqrt(pUr(2)-ala(it))
                         geff_inv(ihli,it) = one/gr - 0.25_pr*(kc/MEFFp(ihli)/Pi**2)*(one+kf/kc &
                                                                                    *Atan(kf/kc))
                      Else
                         geff_inv(ihli,it) = one/gr
                      End If
                   End If
                End If
                ! contribution to Delta
                dvp(ihli)=aka(ihli,it)/geff_inv(ihli,it)
                ! contribution to rearrangement
                !rear_pair =
             End If
          Else
              ! pairing contribution to delta dv(ihli,it)
              If(it.Eq.1) Then
                 dvn(ihli)=(CpV0(it-1)*(ONE-rsa0*CpV1(it-1)))*aka(ihli,it)
              Else
                 dvp(ihli)=(CpV0(it-1)*(ONE-rsa0*CpV1(it-1)))*aka(ihli,it)
              End If
             ! pairing contribution to rearrangement term
              If(use_TMR_pairing.eq.0) then
                 rear_pair = CpV0(it-1) *CpV1(it-1) /rho_c*aka(ihli,it)**2 &
                           + CpV0(ita-1)*CpV1(ita-1)/rho_c*aka(ihli,ita)**2
                 pUr(it)=pUr(it)-rear_pair
              endif
          End If
       End Do !it
       !
       vn(ihli)=pUr(1)       ; vp(ihli)=pUr(2)        !* RHO_ij
       vhbn(ihli)=pUt(1)     ; vhbp(ihli)=pUt(2)      !* TAU_ij
       vrn(ihli)=pUNr(1)     ; vrp(ihli)=pUNr(2)      !* NABLAr RHO__ij
       vzn(ihli)=pUNz(1)     ; vzp(ihli)=pUNz(2)      !* NABLAz RHO__ij
       vdn(ihli)=pUDr(1)     ; vdp(ihli)=pUDr(2)      !* DELTA RHO_ij
       vsn(ihli)=pUDj(1)     ; vsp(ihli)=pUDj(2)      !* NABLA . J__ij
       vSFIZn(ihli)=pUFIZ(1) ; vSFIZp(ihli)=pUFIZ(2)  !* JFIZ_ij
       vSZFIn(ihli)=pUZFI(1) ; vSZFIp(ihli)=pUZFI(2)  !* JZFI_ij
       vSFIRn(ihli)=pUFIR(1) ; vSFIRp(ihli)=pUFIR(2)  !* JFIR_ij
       vSRFIn(ihli)=pURFI(1) ; vSRFIp(ihli)=pURFI(2)  !* JRFI_ij
       !
    End Do !ihli
    !
    If (debug_solver.Ge.1) Call get_CPU_time('field',1)
    !
  End Subroutine field
  !=======================================================================
  !> Computes the constraining fields at the Gauss quadrature points
  !=======================================================================
  Subroutine constraining_field()
    Implicit None
    Integer(ipr) :: ihli,lambda,icons
    Real(pr) :: z,rrr
    Real(pr), Dimension(0:8) :: Qval
    ! fields
    Do icons=1,numberCons
       lambda=multLambda(icons)
       ! Regular multipole moments
       If(lambda.Ge.1) Then
          Do ihli=1,nghl
             z=fh(ihli); rrr=fl(ihli)**2
             Call moments_valueMesh(z,rrr,Qval)
             qfield(ihli,lambda) = Qval(lambda)
          End do
       End If
       ! Gaussian neck operator
       If(lambda.Eq.0) Then
          Do ihli=1,nghl
             z=fh(ihli)
             qfield(ihli,lambda) = Exp(-((z-Z_NECK*bz)/AN_VAL)**2)
          End Do
       End If
    End Do
  End Subroutine constraining_field
  !=======================================================================
  !> Computes the mean field and pairing field of the HFB matrix in
  !> configuration space
  !=======================================================================
  Subroutine gamdel(WoodsSaxon,DoMixing)
    Implicit None
    logical, intent(in) :: WoodsSaxon,DoMixing
    Integer(ipr) :: i,ih,il,ib,ibx,nd,nd2,nza,nra,nla,nsa,nsb,nsab,icons,lambda
    Integer(ipr) :: ihil,laplus,im,JA,N1,N2,ndnd,n12,n21
    Integer(ipr) :: i1,i2,i3
    Real(pr) :: qla,yi,y,y2,qha,qhla,xmi,u2,un,up,xxx
    Real(pr) :: sml2,cnzaa,cnraa,SSU,SSD
    Real(pr) :: FITW1,FITW2,FITW3,FITW4
    Real(pr) :: fi1r,fi1z,fi2d,QHL1A,QH1LA
    Real(pr) :: vh,vdh,vsh,hbh,vsum
    Real(pr) :: SRFIh,SFIRh,SFIZh,SZFIh,SNABLARh,SNABLAZh
    Real(pr) :: xlam,xlam2,xlamy,xlamy2,xlap,xlap2,xlapy,xlapy2,XLAMPY
    Real(pr) :: FIUN1,FIDN1,FIURN1,FIDRN1,FIUZN1,FIDZN1,FIUD2N1,FIDD2N1
    Real(pr) :: FIUN2,FIDN2,FIURN2,FIDRN2,FIUZN2,FIDZN2,FIUD2N2,FIDD2N2
    Real(pr) :: FIUN12,FIDN12,FIURN12,FIDRN12,FIUZN12,FIDZN12
    Real(pr) :: vnhl,vrnhl,vznhl,vdnhl,vsnhl,vhbnhl,vSRFInhl,vSFIRnhl
    Real(pr) :: vSFIZnhl,vSZFInhl,vphl,vrphl,vzphl,vdphl,vsphl,vhbphl
    Real(pr) :: vSRFIphl,vSFIRphl,vSFIZphl,vSZFIphl,dvnhl,dvphl
    Integer(ipr) :: ibro
    Integer(ipr) :: ndxmax
    Parameter(ndxmax=(n00max+2)*(n00max+2)/4)
    Real(pr) :: OMPFIU(ndxmax),OMPFID(ndxmax),OMPFIUR(ndxmax),OMPFIDR(ndxmax),OMPFIUZ(ndxmax), &
                OMPFIDZ(ndxmax),OMPFIUD2N(ndxmax),OMPFIDD2N(ndxmax)
    Integer(ipr) :: t1,t2,countrate,countmax
    !
    If (debug_solver.Ge.1) Call get_CPU_time('gamdel',0)
    !
    If (neck_constraints) Call QNFIND()
    !
    !----------------------------------------------
    ! START BLOCKS
    !----------------------------------------------
    brout=zero; ibro=0
    If (.Not. Allocated(allibro)) Then
       Allocate(allibro(1:NB))
       allibro(1)=0
       Do ib=2,NB
          allibro(ib) = allibro(ib-1) + (ID(ib-1)*(ID(ib-1)+1)/2)
       End Do
    End If
!$OMP PARALLEL DO        &
!$OMP& DEFAULT(NONE)     &
!$OMP& SCHEDULE(DYNAMIC) &
!$OMP& SHARED(NB,ID,IA,NBX,NS,nghl, &
!$OMP&        NHHDIM2,NHHDIM3,NHHDIM4,allibro, &
!$OMP&        vSRFIn,vSFIRn,vSFIZn,vSZFIn, &
!$OMP&        vSRFIp,vSFIRp,vSFIZp,vSZFIp, &
!$OMP&        vn,vrn,vzn,vdn,vsn,vhbn,dvn, &
!$OMP&        vp,vrp,vzp,vdp,vsp,vhbp,dvp, &
!$OMP&        QHLA_opt,FI1R_opt, FI1Z_opt, FI2D_opt, y_opt, &
!$OMP&        nhhdim,kindhfb,ALA2,RK,brout,Parity) &
!$OMP& PRIVATE(I,ND,IB,IM,IBX,LAPLUS,XLAM,XLAP,XLAM2,IL,IH,IHIL,Y,Y2, &
!$OMP&         XLAMY,XLAMY2,XLAP2,XLAPY,XLAPY2,XLAMPY,N1,JA,NSA,SSU,SSD, &
!$OMP&         vnhl,vrnhl,vznhl,vdnhl,vsnhl,vhbnhl,dvnhl, &
!$OMP&         vphl,vrphl,vzphl,vdphl,vsphl,vhbphl,dvphl, &
!$OMP&         vSRFInhl,vSFIRnhl,vSFIZnhl,vSZFInhl,&
!$OMP&         vSRFIphl,vSFIRphl,vSFIZphl,vSZFIphl,&
!$OMP&         FI2D,i1,i2,i3,NSB,NSAB,SNABLARh, SNABLAZh,FI1R,FI1Z, &
!$OMP&         FIUD2N1,FIDD2N1,FIUD2N2,FIDD2N2,FITW3,FITW4,&
!$OMP&         OMPFIUD2N,OMPFIDD2N,OMPFIU,OMPFIUR,OMPFIUZ,OMPFID,OMPFIDR,OMPFIDZ, &
!$OMP&         FIUN1,FIDN1,FIURN1,FIDRN1,FIUZN1,FIDZN1,N2,FIUN2,FIDN2,FIURN2,  &
!$OMP&         FIDRN2,FIUZN2,FIDZN2,FIUN12,FIDN12,FIURN12,FIDRN12,FIUZN12,FIDZN12,VH,&
!$OMP&         HBH,VDH,VSH,SRFIH,SFIRH,SFIZH,SZFIH,UN,UP,N12,QHLA)
    Do ib=1,NB
       ND=ID(ib); IM=ia(ib); ibx=ib+nbx
       If(Parity) Then
          LAPLUS=(ib+1)/2 !Yesp
       Else
          LAPLUS=ib       !Nop
       End If
       XLAP=LAPLUS; XLAM=XLAP-ONE; xlap2=xlap*xlap; xlam2=xlam*xlam
       !----------------------------------------------
       ! SUM OVER GAUSS INTEGRATION POINTS
       !----------------------------------------------
       Do ihil=1,nghl
          y=y_opt(ihil); xlamy=xlam*y;     xlapy=xlap*y;   XLAMPY=XLAMY+XLAPY
          y2=y*y;        xlamy2=xlam2*y2;  xlapy2=xlap2*y2
          ! Neutron mean field
          vnhl=vn(ihil);         vrnhl=vrn(ihil);       vznhl=vzn(ihil);       vdnhl=vdn(ihil)
          vsnhl=vsn(ihil);       vhbnhl=vhbn(ihil);     vSRFInhl=vSRFIn(IHIL); vSFIRnhl=vSFIRn(IHIL)
          vSFIZnhl=vSFIZn(IHIL); vSZFInhl=vSZFIn(IHIL)
          ! Proton mean field
          vphl=vp(ihil);         vrphl=vrp(ihil);       vzphl=vzp(ihil);       vdphl=vdp(ihil)
          vsphl=vsp(ihil);       vhbphl=vhbp(ihil);     vSRFIphl=vSRFIp(IHIL); vSFIRphl=vSFIRp(IHIL)
          vSFIZphl=vSFIZp(IHIL); vSZFIphl=vSZFIp(IHIL)
          ! Pairing field (neutrons and protons)
          dvnhl=dvn(ihil);       dvphl=dvp(ihil)
          !
          Do N1=1,ND
             JA=IM+N1;               NSA=NS(JA);             SSU=Max(NSA,0);         SSD=Max(-NSA,0)
             QHLA=QHLA_opt(JA,ihil); FI1R=FI1R_opt(JA,ihil); FI1Z=FI1Z_opt(JA,ihil); FI2D=FI2D_opt(JA,ihil)
             OMPFIU(N1)=QHLA*SSU;    OMPFIUR(N1)=fi1r*SSU
             OMPFIUZ(N1)=fi1z*SSU;   OMPFIUD2N(N1)=(FI2D-XLAMY2*QHLA)*SSU
             OMPFID(N1)=QHLA*SSD;    OMPFIDR(N1)=fi1r*SSD
             OMPFIDZ(N1)=fi1z*SSD;   OMPFIDD2N(N1)=(FI2D-XLAPY2*QHLA)*SSD
          End Do
          !
          i=allibro(ib)
          Do N1=1,ND
             JA=IM+N1;               NSA=NS(JA)
             FIUN1=OMPFIU(N1);       FIURN1=OMPFIUR(N1);
             FIUZN1=OMPFIUZ(N1);     FIUD2N1=OMPFIUD2N(N1)
             FIDN1=OMPFID(N1);       FIDRN1=OMPFIDR(N1);
             FIDZN1=OMPFIDZ(N1);     FIDD2N1=OMPFIDD2N(N1)
             Do N2=1,N1
                I=I+1; i1=i+nhhdim; i2=i+nhhdim2; i3=i+nhhdim3; NSB=NS(N2+IM); NSAB=NSA+NSB
                If (NSAB.Ne.0) Then
                   If (NSB.Gt.0) Then                                    !spin:UpUp
                      FIUN2    = OMPFIU(N2);    FIURN2 = OMPFIUR(N2)
                      FIUD2N2  = OMPFIUD2N(N2); FIUZN2 = OMPFIUZ(N2)
                      vh       = FIUN1*FIUN2
                      hbh      = vh*XLAMY2+FIURN1*FIURN2+FIUZN1*FIUZN2
                      vdh      = hbh+hbh+FIUN1*FIUD2N2+FIUN2*FIUD2N1
                      SNABLARh = FIURN1*FIUN2+FIURN2*FIUN1
                      SNABLAZh = FIUZN1*FIUN2+FIUZN2*FIUN1
                      vsh      = SNABLARh*XLAMY
                      SFIZh    = (vh+vh)*XLAMY ! =SFIZh (v103)
                   Else                                                  !spin:DoDo
                      FIDN2    = OMPFID(N2);  FIDRN2  = OMPFIDR(N2);
                      FIDZN2   = OMPFIDZ(N2); FIDD2N2 = OMPFIDD2N(N2)
                      vh       = FIDN1*FIDN2
                      hbh      = vh*XLAPY2+FIDRN1*FIDRN2+FIDZN1*FIDZN2
                      vdh      = hbh+hbh+FIDN1*FIDD2N2+FIDN2*FIDD2N1;
                      SNABLARh = FIDRN1*FIDN2+FIDRN2*FIDN1
                      SNABLAZh = FIDZN1*FIDN2+FIDZN2*FIDN1
                      vsh      =-SNABLARh*XLAPY
                      SFIZh    =-(vh+vh)*XLAPY ! =SFIZh (v103)
                   End If
                   brout(i )=brout(i )+vSFIZnhl*SFIZh+vh*vnhl+SNABLARh*vrnhl+SNABLAZh*vznhl+vdh*vdnhl+vsh*vsnhl+hbh*vhbnhl
                   brout(i1)=brout(i1)+vSFIZphl*SFIZh+vh*vphl+SNABLARh*vrphl+SNABLAZh*vzphl+vdh*vdphl+vsh*vsphl+hbh*vhbphl
                   brout(i2)=brout(i2)+vh*dvnhl
                   brout(i3)=brout(i3)+vh*dvphl
                Else
                   If (NSB.Gt.0) Then                                                                !spin:DoUp
                      !vh=ZERO; hbh=ZERO; vdh=ZERO; SNABLARh=ZERO; SNABLAZh=ZERO; SFIZh=ZERO
                      FIUN2   = OMPFIU(N2);    FIURN2 = OMPFIUR(N2);
                      FIUD2N2 = OMPFIUD2N(N2); FIUZN2 = OMPFIUZ(N2)
                      FITW3   =-FIDZN1*FIUN2; FITW4=FIUZN2*FIDN1
                      vsh     =-FIDRN1*FIUZN2+FIURN2*FIDZN1+FITW3*XLAMY-FITW4*XLAPY
                      SRFIh   =-FIDRN1*FIUN2+FIURN2*FIDN1
                      SFIRh   = FIDN1*FIUN2*XLAMPY
                      SZFIh   = FITW3+FITW4
                   Else                                                                             !spin:UpDo
                      !vh=ZERO; hbh=ZERO; vdh=ZERO; SNABLARh=ZERO; SNABLAZh=ZERO; SFIZh=ZERO
                      FIDN2  = OMPFID(N2);   FIDRN2  = OMPFIDR(N2);
                      FIDZN2 = OMPFIDZ(N2);  FIDD2N2 = OMPFIDD2N(N2)
                      FITW3  =-FIDZN2*FIUN1; FITW4=FIUZN1*FIDN2
                      vsh    = FIURN1*FIDZN2-FIDRN2*FIUZN1-FITW4*XLAPY+FITW3*XLAMY ! -vsh (v103)
                      SRFIh  = FIURN1*FIDN2-FIDRN2*FIUN1 !=SRFIh (v103)
                      SFIRh  = FIUN1*FIDN2*XLAMPY !=SFIRh(v103)
                      SZFIh  = FITW3+FITW4 !=SZFIh(v103)
                   End If
                   brout(i )=brout(i )+vsh*vsnhl+vSRFInhl*SRFIh+vSFIRnhl*SFIRh+vSZFInhl*SZFIh
                   brout(i1)=brout(i1)+vsh*vsphl+vSRFIphl*SRFIh+vSFIRphl*SFIRh+vSZFIphl*SZFIh
                End If
                !----------------------------------------------
                ! LN PH PART
                !----------------------------------------------
                If(kindhfb.Lt.0) Then
                   If(ihil.Eq.1) Then
                      un=zero; up=zero;
                      If(N1.Eq.N2) Then
                         un=-ala2(1); up=-ala2(2)
                      End If
                      n12=N1+(N2-1)*ND
                      brout(i )=brout(i )+two*(ala2(1)*rk(n12,ib )+un)
                      brout(i1)=brout(i1)+two*(ala2(2)*rk(n12,ibx)+up)
                   End If
                End If
             End Do !N2
          End Do !N1
       End Do !ihil
    End Do !IB
!$OMP End Parallel Do
    If (debug_solver.Ge.1) Call get_CPU_time('gamdel',1)
    !
    If(finite_range.or.coulomb_gaussian) Then
       Call allocate_fr()
       If(.not.WoodsSaxon) then
          Call system_clock(t1,countrate,countmax)
          Call gamdel_gogny
          Call system_clock(t2,countrate,countmax)
       wct_gogny = (t2-t1)/real(countrate,kind=pr)
       End If
    Else
       wct_gogny = 0._pr
    End If
    ! Lagrange parameters for the constraints
    Do lambda=1,lambdaMax
       brout(nhhdim4+lambda)=multLag(lambda)
    End Do
    If(neck_constraints) Then
       brout(nhhdim4+lambdaMax+1)=neckLag
    End If
    !----------------------------------------------
    ! BROYDEN/LINEAR MIXING
    !----------------------------------------------
    !
    If(debug_solver.Ge.1) Call get_CPU_time('broyden',0)
    !
    If(DoMixing) Then
       If(neck_constraints) Then
          Call broyden_min(nhhdim4+lambdaMax+1,brout,brin,alphamix,si,iiter,nbroyden,bbroyden)
       Else
          Call broyden_min(nhhdim4+lambdaMax,brout,brin,alphamix,si,iiter,nbroyden,bbroyden)
       End If
    Else
       bbroyden='N'
       brin=brout
    End If
    !
    If(debug_solver.Ge.1) Call get_CPU_time('broyden',1)
    !
    Do lambda=1,lambdaMax
       multLag(lambda)=brin(nhhdim4+lambda)
    End Do
    If(neck_constraints) Then
       neckLag=brin(nhhdim4+lambdaMax+1)
    End If
    !
  End Subroutine gamdel
  !=====================================================================
  !>  Calculates the HFB fields Gamma and Delta for a finite range
  !>  Gogny type functional using the axial and radial components of
  !>  two body potential matrix elements calculated with the subroutines
  !>  in the \ref hfbtho_gogny module.
  !>
  !>  If the preprocessor variable GOGNY_SYMMETRIES is set to 1 the
  !>  calculation takes advantage of the separability of the matrix
  !>  elements and first 'folds' the axial matrix elements with the
  !>  density and pairing matrices and creates the objects Zrho and
  !>  Zkappa. Subsequently this objects are 'folded' with the radial
  !>  matrix elements to obtain the Gamma and Delta fields. This
  !>  significantly reduces the number of operations and speeds up the
  !>  calculation. In this mode several other symmetries of the axial
  !>  matrix element are used to reduced the computation time.
  !>
  !>  For any other value of GOGNY_SYMMETRIES the loops to calculate
  !>  the HFB fields goes over every combination of quantum numbers.
  !>  This mode of course requires a larger computation time and is only
  !>  used when debugging the calculation of the fields.
  !======================================================================
  Subroutine gamdel_gogny
    Use UNEDF, Only: e2charg
    !--------------------------------------------------------------------
    ! ph- and pp- matrices in configurational space
    !--------------------------------------------------------------------
    Implicit None
    integer(ipr) :: ib,jb,ind,jnd,im,jm,ibx,jbx,i,ig
    integer(ipr) :: i1,i2,i3
    integer(ipr) :: n1,n2,n3,n4,j1,j2,j3,j4,n43,n12,n21,n34
    integer(ipr) :: nlb,nld,nzac
    integer(ipr) :: iz_abcd,iz_abdc,iz_acbd
    integer(ipr) :: iz_adcb,iz_adbc,iz_acdb
    integer(ipr) :: ir_abcd,ir_abdc,ir_acbd
    integer(ipr) :: ir_adcb,ir_adbc,ir_acdb
    integer(ipr) :: jr_abcd,jr_abdc,jr_acbd
    integer(ipr) :: jr_adcb,jr_adbc,jr_acdb
    real(pr) :: rho_n,rho_p,kap_n,kap_p
    real(pr) :: gamma_n_dir,gamma_n_exc,gamma_p_dir,gamma_p_exc
    real(pr) :: delta_n_dir,delta_n_exc,delta_p_dir,delta_p_exc
    real(pr) :: gamma_c_dir,gamma_c_exc,delta_c_dir
    real(pr) :: Vz_abcd, Vz_abdc
    real(pr) :: V_abcd,V_abdc,V_acbd,V_acdb
    real(pr) :: U_abcd,U_abdc,U_acbd,U_acdb
    real(pr) :: Vdir,Udir,Vexc,Uexc,Vdel,Udel
    real(pr) :: W,B,H,M
    integer(ipr) :: oa,ob,oc,od,la,lb,lc,ld,iba,ibb,ibc,ibd,fact
    real(pr) :: Zrho_n_dir(1:n_g_all,0:nrx,0:nrx,0:1,0:1,1:nox+1)
    real(pr) :: Zrho_p_dir(1:n_g_all,0:nrx,0:nrx,0:1,0:1,1:nox+1)
    real(pr) :: Zrho_n_exc(1:n_g_all,0:nrx,0:nrx,0:1,0:1,1:nox+1)
    real(pr) :: Zrho_p_exc(1:n_g_all,0:nrx,0:nrx,0:1,0:1,1:nox+1)
    real(pr) :: Zkap_n_dir(1:n_g_all,0:nrx,0:nrx,0:1,0:1,1:nox+1)
    real(pr) :: Zkap_p_dir(1:n_g_all,0:nrx,0:nrx,0:1,0:1,1:nox+1)
    real(pr) :: zrho_nd,zrho_pd,zrho_ne,zrho_pe,zkap_n,zkap_p
    integer(ipr) :: ka,pa,ra,za,lmina,lmaxa,sa,ja,na,ida,ima,ibax
    integer(ipr) :: kb,pb,rb,zb,lminb,lmaxb,sb,nb,idb,imb,ibbx
    integer(ipr) :: kc,pc,rc,zc,lminc,lmaxc,sc,jc,nc,idc,imc,ibcx
    integer(ipr) :: kd,pd,rd,zd,lmind,lmaxd,sd,jd,nd,idd,imd,ibdx
    integer(ipr) :: nac,ndb,sac,sdb,n_rlx
    integer(ipr) :: ita,itb,itc,itd,npa,npb,npc,npd,izb,izd
    gamma_g_dir = zero
    gamma_g_exc = zero
    delta_g_dir = zero
    coulf_g_dir = zero
    coulf_g_exc = zero
    coulf_d_dir = zero
#if(GOGNY_SYMMETRIES==1)
    n_rlx = max(2*nrx,nlx)
    if(.not.force_is_dme) then
!$OMP Parallel Default(None)&
!$OMP& SHARED(nzx,nttx,nrr,nll,nss,noo,nzzx,nzz,ib_zrls,id,ia,nbx,&
!$OMP&        i_zrls,rk,ak,VzGogny,ntx,nz,nr,nl,ns,nox,nrx,nrlx,n_rlx,&
!$OMP&        W_g,B_g,H_g,M_g,VrGogny,allibro,nhhdim,nhhdim2,nhhdim3,&
!$OMP&        gamma_g_dir,gamma_g_exc,delta_g_dir,brout,n_g_all,n_g, &
!$OMP&        V_g_coul,coulf_g_dir,coulf_g_exc,coulf_d_dir,icou,e2charg) &
!$OMP& PRIVATE(nzac,Zrho_n_dir,Zrho_p_dir,Zrho_n_exc,Zrho_p_exc,&
!$OMP&         Zkap_n_dir,Zkap_p_dir,zc,za,itb,rb,lb,sb,npb,nlb,itd,npd,&
!$OMP&         rd,ld,sd,nld,izb,zb,ibb,idb,imb,ibbx,jb,nb,izd,zd,ibd,jd,&
!$OMP&         imd,nd,ndb,rho_n,rho_p,kap_n,kap_p,iz_abdc,iz_abcd,&
!$OMP&         iz_acbd,ig,V_abcd,V_abdc,V_acbd,ita,ra,la,sa,iba,ida,ima,&
!$OMP&         ibax,na,itc,rc,lc,sc,ibc,imc,nc,sac,gamma_n_dir,&
!$OMP&         gamma_p_dir,gamma_n_exc,gamma_p_exc,delta_n_dir,&
!$OMP&         delta_p_dir,ob,sdb,ir_abdc,ir_abcd,ir_acbd,jr_abdc,&
!$OMP&         jr_abcd,jr_acbd,W,B,H,M,zrho_nd,zrho_pd,zrho_ne,zrho_pe,&
!$OMP&         zkap_n,zkap_p,U_abcd,U_abdc,U_acbd,Vdir,Udir,Vexc,Uexc,&
!$OMP&         Vdel,Udel,i,i1,i2,i3,nac,&
!$OMP&         gamma_c_dir,gamma_c_exc,delta_c_dir)
!$OMP DO SCHEDULE(DYNAMIC)
    do nzac = 1,(nzx+1)**2
       Zrho_n_dir = 0
       Zrho_p_dir = 0
       Zrho_n_exc = 0
       Zrho_p_exc = 0
       Zkap_n_dir = 0
       Zkap_p_dir = 0
       zc = mod(nzac-1,nzx+1)
       za = (nzac-zc-1)/(nzx+1)
       do itb = 1,nttx
          rb = nrr(itb); lb = nll(itb); sb = nss(itb); npb = noo(itb)
          nlb = mod(npb-lb+1,2)
          do itd = 1,nttx
             npd = noo(itd)
             if(npb.ne.npd) cycle
             rd = nrr(itd); ld = nll(itd); sd = nss(itd)
             nld = mod(npd-ld+1,2)
             do izb = 1,nzzx(itb)
                zb = nzz(itb,izb)
                ibb = ib_zrls(zb,rb,lb,(sb+1)/2)
                idb = id(ibb)
                imb=ia(ibb); ibbx=ibb+nbx
                jb = i_zrls(zb,rb,lb,(sb+1)/2)
                nb = jb-imb
                do izd = 1,nzzx(itd)
                   zd = nzz(itd,izd)
                   if(mod(za+zb+zc+zd,2).ne.0) cycle
                   ibd = ib_zrls(zd,rd,ld,(sd+1)/2)
                   if(ibb.ne.ibd) cycle
                   jd = i_zrls(zd,rd,ld,(sd+1)/2)
                   imd=ia(ibd)
                   nd = jd - imd
                   ndb = nd + (nb-1)*idb
                   rho_n = rk(ndb,ibb )*0.5_pr !rho_db
                   rho_p = rk(ndb,ibbx)*0.5_pr
                   kap_n = ak(ndb,ibb )        !kappa_db
                   kap_p = ak(ndb,ibbx)
                   iz_abdc=zindex(za,zb,zd,zc)
                   if(zc.eq.zd) then
                      iz_abcd=iz_abdc
                   else
                      iz_abcd=zindex(za,zb,zc,zd)
                   endif
                   if(zd.eq.zb) then
                      iz_acbd=iz_abdc
                   elseif(zc.eq.zb) then
                      iz_acbd=iz_abcd
                   else
                      iz_acbd=zindex(za,zc,zb,zd)
                   endif
                   do ig = 1,n_g_all
                      V_abcd = VzGogny(ig,iz_abcd)
                      V_abdc = VzGogny(ig,iz_abdc)
                      V_acbd = VzGogny(ig,iz_acbd)
                      Zrho_n_dir(ig,rd,rb,nld,nlb,npb) = &
                           Zrho_n_dir(ig,rd,rb,nld,nlb,npb)+V_abcd*rho_n
                      Zrho_p_dir(ig,rd,rb,nld,nlb,npb) = &
                           Zrho_p_dir(ig,rd,rb,nld,nlb,npb)+V_abcd*rho_p
                      Zrho_n_exc(ig,rd,rb,nld,nlb,npb) = &
                           Zrho_n_exc(ig,rd,rb,nld,nlb,npb)+V_abdc*rho_n
                      Zrho_p_exc(ig,rd,rb,nld,nlb,npb) = &
                           Zrho_p_exc(ig,rd,rb,nld,nlb,npb)+V_abdc*rho_p
                      Zkap_n_dir(ig,rd,rb,nld,nlb,npb) = &
                           Zkap_n_dir(ig,rd,rb,nld,nlb,npb)+V_acbd*kap_n
                      Zkap_p_dir(ig,rd,rb,nld,nlb,npb) = &
                           Zkap_p_dir(ig,rd,rb,nld,nlb,npb)+V_acbd*kap_p
                   enddo !ig
                enddo ! izd
             enddo !izb
          enddo !itd
       enddo !itb
       do ita = 1,ntx
          if(za.ne.nz(ita)) cycle
          ra = nr(ita); la = nl(ita); sa = ns(ita)
          iba = ib_zrls(za,ra,la,(sa+1)/2)
          ida=id(iba); ima=ia(iba); ibax=iba+nbx
          na = ita-ima
          do itc = 1, ita
             if(zc.ne.nz(itc)) cycle
             rc = nr(itc); lc = nl(itc); sc = ns(itc)
             ibc = ib_zrls(zc,rc,lc,(sc+1)/2)
             if(ibc.ne.iba) cycle
             imc=ia(ibc)
             nc = itc-imc
             sac = sa + sc
             gamma_n_dir = zero
             gamma_p_dir = zero
             gamma_n_exc = zero
             gamma_p_exc = zero
             delta_n_dir = zero
             delta_p_dir = zero
             gamma_c_dir = zero
             gamma_c_exc = zero
             delta_c_dir = zero
             do ob = 0,nox
                ibb = ob + 1
                do nlb = 0,1
                   lb = ob+nlb
                   sb = -2*nlb+1
                   do nld = 0,1
                      ld = ob+nld
                      if(la+lb.ne.lc+ld.and.la-lb.ne.lc-ld) cycle
                      sd = -2*nld+1
                      sdb = sd + sb
                      do rb = 0,nrx
                         if(2*rb+lb.gt.nrlx) cycle
                         do rd = 0,nrx
                            if(2*rd+ld.gt.nrlx) cycle
                            if(sac.ne.0) then
                               ir_abdc=rindex(&
                                    ra,rb,rd,rc,la, lb, ld, lc,n_rlx)
                               if(rb.eq.rd.and.lb.eq.ld) then
                                  ir_acbd=ir_abdc
                               else
                                  ir_acbd=rindex(&
                                       ra,rc,rb,rd,la,-lc,lb,-ld,n_rlx)
                               endif
                               if(rc.eq.rd.and.lc.eq.ld) then
                                  ir_abcd = ir_abdc
                               elseif(rc.eq.rb.and.lc.eq.lb) then
                                  ir_abcd = ir_acbd
                               else
                                  ir_abcd=rindex(&
                                       ra,rb,rc,rd,la,lb,lc,ld,n_rlx)
                               endif
                               if(lb.eq.0.and.ld.eq.0) then
                                  jr_abdc = ir_abdc
                                  jr_acbd = ir_acbd
                                  jr_abcd = ir_abcd
                               else
                                  jr_abdc=rindex(&
                                       ra,rb,rd,rc,la,-lb,-ld,lc,n_rlx)
                                  if(rb.eq.rd .and.lb.eq.ld) then
                                     jr_acbd = jr_abdc
                                  else
                                     jr_acbd=rindex(&
                                          ra,rc,rb,rd,la,-lc,-lb,ld,n_rlx)
                                  endif
                                  if(rc.eq.rb.and.lb.eq.lc) then
                                     jr_abcd = ir_acbd
                                  elseif(rc.eq.rd.and.lc.eq.ld) then
                                     jr_abcd = ir_abdc
                                  else
                                     jr_abcd=rindex(&
                                          ra,rb,rc,rd,la,-lb,lc,-ld,n_rlx)
                                  endif
                               endif
                            else
                               if(sa.eq.sb) then
                                  jr_abdc=rindex(&
                                       ra,rb,rd,rc,la,-lb,-ld, lc,n_rlx)
                                  ir_acbd=rindex(&
                                       ra,rc,rb,rd,la,-lc, lb,-ld,n_rlx)
                               else
                                  ir_abdc=rindex(&
                                       ra,rb,rd,rc,la, lb, ld,lc,n_rlx)
                                  jr_acbd=rindex(&
                                       ra,rc,rb,rd,la,-lc,-lb,ld,n_rlx)
                               endif
                            endif
                            do ig = 1,n_g
                               W = W_g(ig)
                               B = B_g(ig)
                               H = H_g(ig)
                               M = M_g(ig)
                               zrho_nd = Zrho_n_dir(ig,rd,rb,nld,nlb,ibb)
                               zrho_pd = Zrho_p_dir(ig,rd,rb,nld,nlb,ibb)
                               zrho_ne = Zrho_n_exc(ig,rd,rb,nld,nlb,ibb)
                               zrho_pe = Zrho_p_exc(ig,rd,rb,nld,nlb,ibb)
                               zkap_n  = Zkap_n_dir(ig,rd,rb,nld,nlb,ibb)
                               zkap_p  = Zkap_p_dir(ig,rd,rb,nld,nlb,ibb)
                               if(sac.ne.0.and.sdb.ne.0) then
                                  V_abdc=VrGogny(ig,ir_abdc)
                                  U_abdc=VrGogny(ig,jr_abdc)
                                  V_abcd=VrGogny(ig,ir_abcd)
                                  U_abcd=VrGogny(ig,jr_abcd)
                                  V_acbd=VrGogny(ig,ir_acbd)
                                  U_acbd=VrGogny(ig,jr_acbd)
                                  if(sa.eq.sb) then
                                     !up up up up|or|down down down down
                                     Vdir = V_abcd
                                     Udir = U_abcd
                                     Vexc = V_abdc
                                     Uexc = U_abdc
                                     Vdel = V_acbd
                                     Udel = U_acbd
                                  else
                                     !up up down down|or|down down up up
                                     Vdir = U_abcd
                                     Udir = V_abcd
                                     Vexc = U_abdc
                                     Uexc = V_abdc
                                     Vdel = U_acbd
                                     Udel = V_acbd
                                  endif
                                  gamma_n_dir = gamma_n_dir &
                                       + Vdir*((W+B-M-H)*zrho_nd +&
                                       (W+B)*zrho_pd)&
                                       + Udir*(    (W-H)*zrho_nd +&
                                       W*zrho_pd)
                                  gamma_p_dir = gamma_p_dir &
                                       + Vdir*((W+B-M-H)*zrho_pd +&
                                       (W+B)*zrho_nd)&
                                       + Udir*(    (W-H)*zrho_pd +&
                                       W*zrho_nd)
                                  gamma_n_exc = gamma_n_exc &
                                       + Vexc*((M+H-W-B)*zrho_ne +&
                                       (M+H)*zrho_pe)&
                                       + Uexc*(    (M-B)*zrho_ne +&
                                       M*zrho_pe)
                                  gamma_p_exc = gamma_p_exc &
                                       + Vexc*((M+H-W-B)*zrho_pe +&
                                       (M+H)*zrho_ne)&
                                       + Uexc*(    (M-B)*zrho_pe +&
                                       M*zrho_ne)
                                  delta_n_dir = delta_n_dir &
                                       + (Vdel*(W-H)+Udel*(M-B))*zkap_n
                                  delta_p_dir = delta_p_dir &
                                       + (Vdel*(W-H)+Udel*(M-B))*zkap_p
                               elseif(sac.eq.0.and.sdb.eq.0) then
                                  if(sa.eq.sb) then
                                     !up down down up|or|down up up down
                                     Vexc = -VrGogny(ig,jr_abdc)
                                     Vdel =  VrGogny(ig,ir_acbd)
                                  else
                                     !up down up down|or|down up down up
                                     Vexc =  VrGogny(ig,ir_abdc)
                                     Vdel = -VrGogny(ig,jr_acbd)
                                  endif
                                  gamma_n_exc = gamma_n_exc &
                                       + Vexc*((H-W)*zrho_ne + H*zrho_pe)
                                  gamma_p_exc = gamma_p_exc &
                                       + Vexc*((H-W)*zrho_pe + H*zrho_ne)
                                  delta_n_dir = delta_n_dir &
                                       + Vdel*zkap_n*(W-H+B-M)
                                  delta_p_dir = delta_p_dir &
                                       + Vdel*zkap_p*(W-H+B-M)
                               endif
                            enddo !ig
                            do ig = n_g+1,n_g_all
                               W = V_g_coul(ig-n_g) * e2charg
                               zrho_pd = Zrho_p_dir(ig,rd,rb,nld,nlb,ibb)
                               if(icou.le.-4) then
                                  zrho_pe = 0._pr
                                  zkap_p  = 0._pr
                               elseif(icou.le.-3) then
                                  zrho_pe = Zrho_p_exc(ig,rd,rb,nld,nlb,ibb)
                                  zkap_p  = Zkap_p_dir(ig,rd,rb,nld,nlb,ibb)
                               elseif(icou.le.-2) then
                                  zrho_pe = Zrho_p_exc(ig,rd,rb,nld,nlb,ibb)
                                  zkap_p  = 0._pr
                               else
                                  zrho_pe = 0._pr
                                  zkap_p  = 0._pr
                               endif
                               if(sac.ne.0.and.sdb.ne.0) then
                                  if(sa.eq.sb) then
                                     !up up up up|or|down down down down
                                     Vdir = VrGogny(ig,ir_abcd)
                                     Udir = VrGogny(ig,jr_abcd)
                                     if(icou.le.-4) then
                                        Vexc = 0._pr
                                        Vdel = 0._pr
                                     elseif(icou.le.-3) then
                                        Vexc = VrGogny(ig,ir_abdc)
                                        Vdel = VrGogny(ig,ir_acbd)
                                     elseif(icou.le.-2) then
                                        Vexc = VrGogny(ig,ir_abdc)
                                        Vdel = 0._pr
                                     else
                                        Vexc = 0._pr
                                        Vdel = 0._pr
                                     endif
                                  else
                                     !up up down down|or|down down up up
                                     Vdir = VrGogny(ig,jr_abcd)
                                     Udir = VrGogny(ig,ir_abcd)
                                     if(icou.le.-4) then
                                        Vexc = 0._pr
                                        Vdel = 0._pr
                                     elseif(icou.le.-3) then
                                        Vexc = VrGogny(ig,jr_abdc)
                                        Vdel = VrGogny(ig,jr_acbd)
                                     elseif(icou.le.-2) then
                                        Vexc = VrGogny(ig,jr_abdc)
                                        Vdel = 0._pr
                                     else
                                        Vexc = 0._pr
                                        Vdel = 0._pr
                                     endif
                                  endif
                                  gamma_c_dir = gamma_c_dir &
                                       + (Vdir+Udir)*W*zrho_pd
                                  gamma_c_exc = gamma_c_exc &
                                       - Vexc*W*zrho_pe
                                  delta_c_dir = delta_c_dir &
                                       + Vdel*W*zkap_p
                               elseif(sac.eq.0.and.sdb.eq.0.and.icou.le.-2) then
                                  if(sa.eq.sb) then
                                     !up down down up|or|down up up down
                                     Vexc = -VrGogny(ig,jr_abdc)
                                     if(icou.eq.-3) then
                                        Vdel =  VrGogny(ig,ir_acbd)
                                     else
                                        Vdel = 0._pr
                                     endif
                                  else
                                     !up down up down|or|down up down up
                                     Vexc =  VrGogny(ig,ir_abdc)
                                     if(icou.eq.-3) then
                                        Vdel = -VrGogny(ig,jr_acbd)
                                     else
                                        Vdel = 0._pr
                                     endif
                                  endif
                                  gamma_c_exc = gamma_c_exc &
                                       - Vexc*W*zrho_pe
                                  delta_c_dir = delta_c_dir &
                                       + Vdel*zkap_p*W
                               endif
                            enddo !ig
                         enddo !rd
                      enddo !rb
                   enddo !nld
                enddo !nlb
             enddo !ob
             i=allibro(iba)+na*(na-1)/2+nc
             i1=i+nhhdim
             i2=i+nhhdim2
             i3=i+nhhdim3
             nac = na+(nc-1)*ida;
             gamma_g_dir(nac,iba )=gamma_n_dir
             gamma_g_dir(nac,ibax)=gamma_p_dir
             gamma_g_exc(nac,iba )=gamma_n_exc
             gamma_g_exc(nac,ibax)=gamma_p_exc
             delta_g_dir(nac,iba )=delta_n_dir
             delta_g_dir(nac,ibax)=delta_p_dir
             coulf_g_dir(nac,ibax)=gamma_c_dir
             coulf_g_exc(nac,ibax)=gamma_c_exc
             coulf_d_dir(nac,ibax)=delta_c_dir
             brout(i )=brout(i )+gamma_n_dir+gamma_n_exc
             brout(i1)=brout(i1)+gamma_p_dir+gamma_p_exc+gamma_c_dir+gamma_c_exc
             brout(i2)=brout(i2)+delta_n_dir
             brout(i3)=brout(i3)+delta_p_dir+delta_c_dir
          enddo !itc
       enddo !ita
    enddo !nzac
!$OMP End Do
!$OMP End Parallel
 else
!$OMP Parallel Default(None)&
!$OMP& SHARED(nzx,nttx,nrr,nll,nss,noo,nzzx,nzz,ib_zrls,id,ia,nbx,&
!$OMP&        i_zrls,rk,ak,VzGogny,ntx,nz,nr,nl,ns,nox,nrx,nrlx,n_rlx,&
!$OMP&        W_g,B_g,H_g,M_g,VrGogny,allibro,nhhdim,nhhdim2,nhhdim3,&
!$OMP&        gamma_g_dir,gamma_g_exc,delta_g_dir,brout,n_g,force_is_dme) &
!$OMP& PRIVATE(nzac,Zrho_n_dir,Zrho_p_dir,Zrho_n_exc,Zrho_p_exc,&
!$OMP&         Zkap_n_dir,Zkap_p_dir,zc,za,itb,rb,lb,sb,npb,nlb,itd,npd,&
!$OMP&         rd,ld,sd,nld,izb,zb,ibb,idb,imb,ibbx,jb,nb,izd,zd,ibd,jd,&
!$OMP&         imd,nd,ndb,rho_n,rho_p,kap_n,kap_p,iz_abdc,iz_abcd,&
!$OMP&         iz_acbd,ig,V_abcd,V_abdc,V_acbd,ita,ra,la,sa,iba,ida,ima,&
!$OMP&         ibax,na,itc,rc,lc,sc,ibc,imc,nc,sac,gamma_n_dir,&
!$OMP&         gamma_p_dir,gamma_n_exc,gamma_p_exc,delta_n_dir,&
!$OMP&         delta_p_dir,ob,sdb,ir_abdc,ir_abcd,ir_acbd,jr_abdc,&
!$OMP&         jr_abcd,jr_acbd,W,B,H,M,zrho_nd,zrho_pd,zrho_ne,zrho_pe,&
!$OMP&         zkap_n,zkap_p,U_abcd,U_abdc,U_acbd,Vdir,Udir,Vexc,Uexc,&
!$OMP&         Vdel,Udel,i,i1,i2,i3,nac)
!$OMP DO SCHEDULE(DYNAMIC)
    do nzac = 1,(nzx+1)**2
       Zrho_n_dir = 0
       Zrho_p_dir = 0
       zc = mod(nzac-1,nzx+1)
       za = (nzac-zc-1)/(nzx+1)
       do itb = 1,nttx
          rb = nrr(itb); lb = nll(itb); sb = nss(itb); npb = noo(itb)
          nlb = mod(npb-lb+1,2)
          do itd = 1,nttx
             npd = noo(itd)
             if(npb.ne.npd) cycle
             rd = nrr(itd); ld = nll(itd); sd = nss(itd)
             nld = mod(npd-ld+1,2)
             do izb = 1,nzzx(itb)
                zb = nzz(itb,izb)
                ibb = ib_zrls(zb,rb,lb,(sb+1)/2)
                idb = id(ibb)
                imb=ia(ibb); ibbx=ibb+nbx
                jb = i_zrls(zb,rb,lb,(sb+1)/2)
                nb = jb-imb
                do izd = 1,nzzx(itd)
                   zd = nzz(itd,izd)
                   if(mod(za+zb+zc+zd,2).ne.0) cycle
                   ibd = ib_zrls(zd,rd,ld,(sd+1)/2)
                   if(ibb.ne.ibd) cycle
                   jd = i_zrls(zd,rd,ld,(sd+1)/2)
                   imd=ia(ibd)
                   nd = jd - imd
                   ndb = nd + (nb-1)*idb
                   rho_n = rk(ndb,ibb )*0.5_pr !rho_db
                   rho_p = rk(ndb,ibbx)*0.5_pr
                   iz_abcd=zindex(za,zb,zc,zd)
                   do ig = 1,n_g
                      V_abcd = VzGogny(ig,iz_abcd)
                      Zrho_n_dir(ig,rd,rb,nld,nlb,npb) = &
                           Zrho_n_dir(ig,rd,rb,nld,nlb,npb)+V_abcd*rho_n
                      Zrho_p_dir(ig,rd,rb,nld,nlb,npb) = &
                           Zrho_p_dir(ig,rd,rb,nld,nlb,npb)+V_abcd*rho_p
                   enddo !ig
                enddo ! izd
             enddo !izb
          enddo !itd
       enddo !itb
       do ita = 1,ntx
          if(za.ne.nz(ita)) cycle
          ra = nr(ita); la = nl(ita); sa = ns(ita)
          iba = ib_zrls(za,ra,la,(sa+1)/2)
          ida=id(iba); ima=ia(iba); ibax=iba+nbx
          na = ita-ima
          do itc = 1, ita
             if(zc.ne.nz(itc)) cycle
             rc = nr(itc); lc = nl(itc); sc = ns(itc)
             ibc = ib_zrls(zc,rc,lc,(sc+1)/2)
             if(ibc.ne.iba) cycle
             imc=ia(ibc)
             nc = itc-imc
             sac = sa + sc
             if(sac.eq.0) cycle
             gamma_n_dir = zero
             gamma_p_dir = zero
             do ob = 0,nox
                ibb = ob + 1
                do nlb = 0,1
                   lb = ob+nlb
                   sb = -2*nlb+1
                   do nld = 0,1
                      ld = ob+nld
                      if(la+lb.ne.lc+ld.and.la-lb.ne.lc-ld) cycle
                      sd = -2*nld+1
                      sdb = sd + sb
                      if(sdb.eq.0) cycle
                      do rb = 0,nrx
                         if(2*rb+lb.gt.nrlx) cycle
                         do rd = 0,nrx
                            if(2*rd+ld.gt.nrlx) cycle
                            ir_abcd=rindex(&
                                 ra,rb,rc,rd,la,lb,lc,ld,n_rlx)
                            if(lb.eq.0.and.ld.eq.0) then
                               jr_abcd = ir_abcd
                            else
                               jr_abcd=rindex(&
                                    ra,rb,rc,rd,la,-lb,lc,-ld,n_rlx)
                            endif
                            do ig = 1,n_g
                               W = W_g(ig)
                               B = B_g(ig)
                               H = H_g(ig)
                               M = M_g(ig)
                               zrho_nd = Zrho_n_dir(ig,rd,rb,nld,nlb,ibb)
                               zrho_pd = Zrho_p_dir(ig,rd,rb,nld,nlb,ibb)
                               V_abcd=VrGogny(ig,ir_abcd)
                               U_abcd=VrGogny(ig,jr_abcd)
                               if(sa.eq.sb) then
                                  !up up up up|or|down down down down
                                  Vdir = V_abcd
                                  Udir = U_abcd
                               else
                                  !up up down down|or|down down up up
                                  Vdir = U_abcd
                                  Udir = V_abcd
                               endif
                               gamma_n_dir = gamma_n_dir &
                                    + Vdir*((W+B-M-H)*zrho_nd +&
                                    (W+B)*zrho_pd)&
                                    + Udir*(    (W-H)*zrho_nd +&
                                    W*zrho_pd)
                               gamma_p_dir = gamma_p_dir &
                                    + Vdir*((W+B-M-H)*zrho_pd +&
                                    (W+B)*zrho_nd)&
                                    + Udir*(    (W-H)*zrho_pd +&
                                    W*zrho_nd)
                            enddo !ig
                         enddo !rd
                      enddo !rb
                   enddo !nld
                enddo !nlb
             enddo !ob
             i=allibro(iba)+na*(na-1)/2+nc
             i1=i+nhhdim
             i2=i+nhhdim2
             i3=i+nhhdim3
             nac = na+(nc-1)*ida;
             gamma_g_dir(nac,iba )=gamma_n_dir
             gamma_g_dir(nac,ibax)=gamma_p_dir
             brout(i )=brout(i )+gamma_n_dir
             brout(i1)=brout(i1)+gamma_p_dir
          enddo !itc
       enddo !ita
    enddo !nzac
!$OMP End Do
!$OMP End Parallel
 endif
#else
    do ita = 1,ntx
       ra = nr(ita); za = nz(ita); la = nl(ita); sa = ns(ita)
       iba = ib_zrls(za,ra,la,(sa+1)/2)
       ida=id(iba); ima=ia(iba); ibax=iba+nbx
       na = ita-ima
       do itc = 1, ita
          rc = nr(itc); zc = nz(itc); lc = nl(itc); sc = ns(itc)
          ibc = ib_zrls(zc,rc,lc,(sc+1)/2)
          idc=id(ibc); imc=ia(ibc); ibcx=ibc+nbx
          nc = itc-imc
          sac = sa + sc
          if(ibc.ne.iba) cycle
          gamma_n_dir = zero
          gamma_p_dir = zero
          gamma_n_exc = zero
          gamma_p_exc = zero
          delta_n_dir = zero
          delta_p_dir = zero
          do itb = 1,ntx
             rb = nr(itb); zb = nz(itb); lb = nl(itb); sb = ns(itb)
             ibb = ib_zrls(zb,rb,lb,(sb+1)/2)
             idb=id(ibb); imb=ia(ibb); ibbx=ibb+nbx
             nb = itb-imb
             do itd = 1,ntx
                rd = nr(itd); zd = nz(itd); ld = nl(itd); sd = ns(itd)
                ibd = ib_zrls(zd,rd,ld,(sd+1)/2)
                idd=id(ibd); imd=ia(ibd); ibdx=ibd+nbx
                nd = itd-imd
                sdb = sd + sb
                if(ibb.ne.ibd) cycle
                ndb = nd + (nb-1)*idd
                rho_n = rk(ndb,ibb )*0.5_pr !rho_db
                rho_p = rk(ndb,ibbx)*0.5_pr
                kap_n = ak(ndb,ibb )        !kappa_db
                kap_p = ak(ndb,ibbx)
                do ig = 1,n_g
                   W = W_g(ig)
                   B = B_g(ig)
                   H = H_g(ig)
                   M = M_g(ig)
                   if(sac.ne.0.and.sdb.ne.0) then
                      V_abdc=Vr_Gogny(ig,ra,la,rb, lb,rd, ld,rc, lc)
                      U_abdc=Vr_Gogny(ig,ra,la,rb,-lb,rd,-ld,rc, lc)
                      V_abcd=Vr_Gogny(ig,ra,la,rb, lb,rc, lc,rd, ld)
                      U_abcd=Vr_Gogny(ig,ra,la,rb,-lb,rc, lc,rd,-ld)
                      V_acbd=Vr_Gogny(ig,ra,la,rc,-lc,rb, lb,rd,-ld)
                      U_acbd=Vr_Gogny(ig,ra,la,rc,-lc,rb,-lb,rd, ld)
                      if(sa.eq.sb) then
                         ! up up   up up   |or| down down down down
                         Vdir = V_abcd*Vz_Gogny(ig,za,zb,zc,zd)
                         Udir = U_abcd*Vz_Gogny(ig,za,zb,zc,zd)
                         Vexc = V_abdc*Vz_Gogny(ig,za,zb,zd,zc)
                         Uexc = U_abdc*Vz_Gogny(ig,za,zb,zd,zc)
                         Vdel = V_acbd*Vz_Gogny(ig,za,zc,zb,zd)
                         Udel = U_acbd*Vz_Gogny(ig,za,zc,zb,zd)
                      else
                         !up up  down down |or| down down up up
                         Vdir = U_abcd*Vz_Gogny(ig,za,zb,zc,zd)
                         Udir = V_abcd*Vz_Gogny(ig,za,zb,zc,zd)
                         Vexc = U_abdc*Vz_Gogny(ig,za,zb,zd,zc)
                         Uexc = V_abdc*Vz_Gogny(ig,za,zb,zd,zc)
                         Vdel = U_acbd*Vz_Gogny(ig,za,zc,zb,zd)
                         Udel = V_acbd*Vz_Gogny(ig,za,zc,zb,zd)
                      endif
                      gamma_n_dir = gamma_n_dir &
                           + Vdir*((W+B-M-H)*rho_n + (W+B)*rho_p)&
                           + Udir*(    (W-H)*rho_n +     W*rho_p)
                      gamma_p_dir = gamma_p_dir &
                           + Vdir*((W+B-M-H)*rho_p + (W+B)*rho_n)&
                           + Udir*(    (W-H)*rho_p +     W*rho_n)
                      gamma_n_exc = gamma_n_exc &
                           + Vexc*((M+H-W-B)*rho_n + (M+H)*rho_p)&
                           + Uexc*(    (M-B)*rho_n +     M*rho_p)
                      gamma_p_exc = gamma_p_exc &
                           + Vexc*((M+H-W-B)*rho_p + (M+H)*rho_n)&
                           + Uexc*(    (M-B)*rho_p +     M*rho_n)
                      delta_n_dir = delta_n_dir &
                           + (Vdel*(W-H)+Udel*(M-B))*kap_n
                      delta_p_dir = delta_p_dir &
                           + (Vdel*(W-H)+Udel*(M-B))*kap_p
                   elseif(sac.eq.0.and.sdb.eq.0) then
                      if(sa.eq.sb) then
                         !up down down up|or| down up up down
                         Vexc = -Vr_Gogny(ig,ra,la,rb,-lb,rd,-ld,rc,lc)&
                              *Vz_Gogny(ig,za,zb,zd,zc)
                         Vdel =  Vr_Gogny(ig,ra,la,rc,-lc,rb,lb,rd,-ld)&
                              *Vz_Gogny(ig,za,zc,zb,zd)
                      else
                         !up down up down |or| down up down up
                         Vexc = Vr_Gogny(ig,ra,la,rb,lb,rd,ld,rc,lc)&
                              *Vz_Gogny(ig,za,zb,zd,zc)
                         Vdel =-Vr_Gogny(ig,ra,la,rc,-lc,rb,-lb,rd,ld)&
                              *Vz_Gogny(ig,za,zc,zb,zd)
                      endif
                      gamma_n_exc = gamma_n_exc &
                           + Vexc*((H-W)*rho_n + H*rho_p)
                      gamma_p_exc = gamma_p_exc &
                           + Vexc*((H-W)*rho_p + H*rho_n)
                      delta_n_dir = delta_n_dir + Vdel*kap_n*(W-H+B-M)
                      delta_p_dir = delta_p_dir + Vdel*kap_p*(W-H+B-M)
                   endif
                enddo !ig
             enddo !itd
          enddo !itb
          if(force_is_dme) then
             !only the Hartree Part
             gamma_n_exc = zero
             gamma_p_exc = zero
             delta_n_dir = zero
             delta_p_dir = zero
          endif
          i=allibro(iba)+na*(na-1)/2+nc
          i1=i+nhhdim
          i2=i+nhhdim2
          i3=i+nhhdim3
          nac = na+(nc-1)*ida;
          gamma_g_dir(nac,iba )=gamma_n_dir
          gamma_g_dir(nac,ibax)=gamma_p_dir
          gamma_g_exc(nac,iba )=gamma_n_exc
          gamma_g_exc(nac,ibax)=gamma_p_exc
          delta_g_dir(nac,iba )=delta_n_dir
          delta_g_dir(nac,ibax)=delta_p_dir
          brout(i )=brout(i )+gamma_n_dir+gamma_n_exc
          brout(i1)=brout(i1)+gamma_p_dir+gamma_p_exc
          brout(i2)=brout(i2)+delta_n_dir
          brout(i3)=brout(i3)+delta_p_dir
       enddo !itc
    enddo !ita
#endif
  End Subroutine gamdel_gogny

  !=======================================================================
  !               PARTICLE NUMBER PROJECTION MODULE
  !=======================================================================
  !
  !=======================================================================
  !> Records on disk the occupations both of canonical states and of
  !> quasiparticles for use in Marc Verriere statistical projection module
  !=======================================================================
  Subroutine sample_occupations()
    Implicit None
    Character(Len=1) :: type_part
    Integer(ipr) :: it,ib,k,k1,k2,ierr,k0,nd,nstate
    Real(pr) :: pn,ek0,vk0
    !
    Open(55,file='occup_canonical.dat',status='unknown',form='formatted',IOSTAT=ierr)
    Write(55,'(2i4)') npr(1),npr(2)
    Write(55,'("Aneck   0.0000000000000E+0.000")')
    Write(55,'("Sneck   0.0000000000000E+0.000")')
    Write(55,'("Qneck 0 0.0000000000000E+0.000")')
    Write(55,'("Qneck 1 0.0000000000000E+0.000")')
    Write(55,'("n_gh   001")')
    Do it=itmin,itmax
       If(it.Eq.1) type_part='0'
       If(it.Eq.2) type_part='1'
       k0=0; nstate=0
       Do ib=1,nb
          nd=id(ib)
          k1=ka(ib,it)+1
          k2=ka(ib,it)+kd(ib,it)
          If(k1.Le.k2) Then
             Do k=1,nd                   ! print active states only
                k0=k0+1
                ek0=ek(k0,it)            ! canonical s.p. energies
                vk0=vk(k0,it)            ! canonical occupations v^2
                If(vk0.Gt.-1.0e-4_pr) Then   ! print If signIficant v^2
                   nstate=nstate+1
                End If
             End Do
          End If
       End Do
       Write(55,'("basis_",a1,i5)') type_part,nstate
     End Do
     Do it=itmin,itmax
       If(it.Eq.1) type_part='0'
       If(it.Eq.2) type_part='1'
       k0=0
       Do ib=1,nb
          nd=id(ib)
          k1=ka(ib,it)+1
          k2=ka(ib,it)+kd(ib,it)
          If(k1.Le.k2) Then
             Do k=1,nd                   ! print active states only
                k0=k0+1
                ek0=ek(k0,it)            ! canonical s.p. energies
                vk0=vk(k0,it)            ! canonical occupations v^2
                If(vk0.Gt.-1.0e-4_pr) Then   ! print If signIficant v^2
                   pn=vk(k0,it)          ! qp probabilities
                   Write(55,'(a1,1x,e20.13)') type_part,pn
                End If
             End Do
          End If
       End Do !ib
    End Do !it
    Close(55)
    ! Write qp energies and occupations
    Open(55,file='occup_qp.dat',status='unknown',form='formatted',IOSTAT=ierr)
    Write(55,'(2i4)') npr(1),npr(2)
    Write(55,'("Aneck   0.0000000000000E+0.000")')
    Write(55,'("Sneck   0.0000000000000E+0.000")')
    Write(55,'("Qneck 0 0.0000000000000E+0.000")')
    Write(55,'("Qneck 1 0.0000000000000E+0.000")')
    Write(55,'("n_gh   001")')
    Do it=itmin,itmax
       If(it.Eq.1) type_part='0'
       If(it.Eq.2) type_part='1'
       nstate=0
       Do ib=1,nb
          k1=ka(ib,it)+1
          k2=ka(ib,it)+kd(ib,it)
          If(k1.Le.k2) Then
             Do k=k1,k2                               ! print active states only
                pn=uk(k,it)                             ! qp probabilities
                nstate=nstate+1
             End Do
          End If
       End Do
       Write(55,'("basis_",a1,i5)') type_part,nstate
     End Do
     Do it=itmin,itmax
       If(it.Eq.1) type_part='0'
       If(it.Eq.2) type_part='1'
       Do ib=1,nb
          k1=ka(ib,it)+1
          k2=ka(ib,it)+kd(ib,it)
          If(k1.Le.k2) Then
             Do k=k1,k2                                 ! print active states only
                pn=uk(k,it)                             ! qp probabilities
                Write(55,'(a1,1x,e20.13)') type_part,pn
             End Do
          End If
       End Do !ib
    End Do !it
    Close(55)
  End Subroutine sample_occupations
  !=======================================================================
  !> Modified Broyden's method: D.D.Johnson, PRB 38, 12807 (1988)
  !> Adopted from: (C) 2001 PWSCF group
  !=======================================================================
  Subroutine broyden_min(N,vout,vin,alpha,si,iter,M,bbroyden)
    !---------------------------------------------------------------------
    ! Input :
    !  N      dimension of arrays vin,vout
    !  vin    outpu at previous iteration
    !  vout   output at current iteration
    !  alpha  mixing factor (0 < alpha <= 1)
    !  iter   current iteration number
    !  M      number of iterations in Broyden history
    !  M=0    Linear mixing
    ! Output:
    !  si     MaxVal(|vout-vin|)
    !  vin    Broyden/Linear mixing result
    !  vout   vout-vin
    !  bbroyden='B' Broyden mixing, curvature>0
    !  bbroyden='L' Linear mixing,  curvature<0
    !---------------------------------------------------------------------
    Use HFBTHO_utilities, Only: pr,ipr
    Use HFBTHO, Only: ierror_flag,ierror_info
    Implicit None
    Integer(ipr),     Intent(In)    :: N,iter,M
    Real(pr),        Intent(In)     :: alpha
    Real(pr),        Intent(Out)    :: si
    Character(1),      Intent(Out)  :: bbroyden
    Real(pr),        Intent(InOut)  :: vout(N),vin(N)
    Integer(ipr)                    :: i,j,iter_used,ipos,inext,info
    Integer(ipr), Allocatable, Save :: iwork(:)
    Real(pr),    Allocatable, Save  :: beta(:,:),work(:)
    Real(pr),    Allocatable, Save  :: df(:,:),dv(:,:),curv(:)
    Real(pr),                 Save  :: w0
    Real(pr)                        :: DDOT,DNRM2,normi,gamma,curvature,sf
    !
    sf=-1.0_pr; Call DAXPY(N,sf,vin,1,vout,1)
    si=Maxval(Abs(vout))
    ! Linear mixing
    If(M.Le.0) Then
       bbroyden='L'; Call DAXPY(N,xmix,vout,1,vin,1)
       If(iter.Eq.0) Write(6,*) '  Linear mixing (alpha) : ',xmix
       Return
    End If
    If(iter.Eq.0) Then
       bbroyden='L'; Call DAXPY(N,alpha,vout,1,vin,1)
       !If(iter.Eq.0) Write(6,*) '  Linear mixing (alpha) : ',alpha
       Return
    End If
    ! Broyden mixing
    iter_used=Min(iter-1,M)
    ipos=iter-1-((iter-2)/M)*M
    inext=iter-((iter-1)/M)*M
    If (iter.Eq.1) Then
       w0=0.010_pr
       If(Allocated(df)) Deallocate(curv,df,dv,beta,work,iwork)
       Allocate(curv(N),df(N,M),dv(N,M),beta(M,M),work(M),iwork(M))
    Else
       df(:,ipos)=vout(:)-df(:,ipos); dv(:,ipos)=vin(:)-dv(:,ipos)
       Normi=1.0_pr/Sqrt((DNRM2(N,df(1,ipos),1))**2)
       Call DSCAL(N,Normi,df(1,ipos),1)
       Call DSCAL(N,Normi,dv(1,ipos),1)
    End If
    Do i=1,iter_used
       Do j=i+1,iter_used
          beta(i,j)=DDOT(N,df(1, j),1,df(1,i),1)
       End Do
       beta(i,i)=1.0_pr+w0*w0
    End Do
#if(USE_ESSL==0)
    Call DSYTRF('U', iter_used, beta, M, iwork, work, M, info)
#else
    Call DPOTRF('U', iter_used, beta, M, info)
#endif
    If(info.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='STOP: In Broyden: info at DSYTRF '
       Return
    End If
#if(USE_ESSL==0)
    Call DSYTRI('U', iter_used, beta, M, iwork, work, info)
#else
    Call DPOTRI('U', iter_used, beta, M, info)
#endif
    If(info.Ne.0) Then
       ierror_flag=ierror_flag+1
       ierror_info(ierror_flag)='STOP: In Broyden: info at DSYTRI '
       Return
    End If
    Do i=1,iter_used
       Do j=i+1,iter_used
          beta(j,i)=beta(i,j)
       End Do
       work(i)=DDOT(N,df(1,i),1,vout,1)
    End Do
    curv=alpha*vout
    Do i=1,iter_used
       gamma=0.0_pr
       Do j=1,iter_used
          gamma=gamma+beta(j,i)*work(j)
       End Do
       curv=curv-gamma*(dv(:,i)+alpha*df(:,i))
    End Do
    Call DCOPY(N,vout,1,df(1,inext),1)
    Call DCOPY(N,vin,1,dv(1,inext),1)
    curvature=DDOT(N,vout,1,curv,1)
    If(curvature.Gt.-1.0_pr) Then
       bbroyden='B'; sf=+1.0_pr; Call DAXPY(N,sf,curv,1,vin,1)
    Else
       bbroyden='L'; sf=alpha*0.50_pr; Call DAXPY(N,sf,vout,1,vin,1)
    End If
  End Subroutine broyden_min
  !=======================================================================
  !> Calculates expectation values (xn is the particle number)
  !> for lpr=.true. also calculates PAV corrections
  !=======================================================================
  Subroutine expect(lpr)
    Implicit None
    Logical :: lpr
    Integer(ipr) :: i,it,ihli,iw,LAMACT,I_TYPE
    Real(pr) :: ekt(3),xn(3),q4(3),def(3),bet2(3),het4(3),econst
    Real(pr) :: z,zz,rrr,p2,p3,p4,row,r212,r222,rc
    Real(pr) :: eso,ecodi,ecoex,rn,rp,rnp1,rnp2,rt,whl,tnt,tpt,tt
    Real(pr) :: dn,dp,dt,akn,akp,akn2,akp2,adn,adp,evol,esurf,ecoul
    Real(pr) :: etens,dd1n,dd1p,rt1,tt1,dt1,djn,djp,djt,djt1
    Real(pr) :: RHO_0,RHO_1,TAU_0,TAU_1,DRHO_0,DRHO_1,DJ_0,DJ_1,J2_0,J2_1,JabJba_0,JabJba_1
    Real(pr) :: SZFIN,SFIZN,SRFIN,SFIRN,SZFIP,SFIZP,SRFIP,SFIRP
    Real(pr) :: SZFI_0,SFIZ_0,SRFI_0,SFIR_0,SZFI_1,SFIZ_1,SRFI_1,SFIR_1
    Real(pr) :: SNABLARN,SNABLAZN,SNABLARP,SNABLAZP
    Real(pr) :: SNABLAR_0,SNABLAZ_0,SNABLAR_1,SNABLAZ_1
    Real(pr) :: xn1,xn2,rms1,rms2,q21,q22,q41,q42,EKIN_N,EKIN_P,ept1,ept2,del1,del2,rsa0
    Real(pr) :: ESURF_rho_DELTA_rho,ESURF_NABLA_rho_NABLA_rho,ESO_rho_NABLA_J,ESO_NABLA_rho_J
    Real(pr) :: EVOL_rho_tau,EVOL_rho_rho,EExtra,E_HARTREE_DIR,tempE_Crho0,tempREARR
    Real(pr) :: E_EXT_FIELD
    Real(pr) :: QLMLEF,QLMRIG
    Real(pr) :: RpRMSsq,RnRMSsq,DarwinFoldy
    Real(pr) :: e_gamma_fr_dir,e_gamma_fr_exc
    Real(pr) :: delrho1,delrho2,delkap1,delkap2
    Real(pr) :: cdelrho1,cdelrho2,cdelkap1,cdelkap2
    Real(pr), Pointer     :: EqpPo(:),VqpPo(:),UqpPo(:)
    Integer(ipr), Pointer :: KpwiPo(:),KqpPo(:)
    !------------------------------------------------
    ! Part called during iterations (lpr=F)
    !------------------------------------------------
      !
    Call DENSIT
    If(ierror_flag.Ne.0) Return
    !
    If (debug_solver.Eq.1) Call get_CPU_time('expect',0)
    !
    !------------------------------------------------
    ! zero energy variables
    !------------------------------------------------
    EKIN_N=zero;              EKIN_P=zero;
    EVOL_rho_tau=zero;        EVOL_rho_rho=zero;
    ESURF_rho_DELTA_rho=zero; ESURF_NABLA_rho_NABLA_rho=zero;
    ESO_rho_NABLA_J=zero;     ESO_NABLA_rho_J=zero; E_HARTREE_DIR=zero
    ept1=zero; ept2=zero;     del1=zero; del2=zero;
    ecodi=zero; ecoex=zero; etens=zero
    EExtra=zero ; E_EXT_FIELD = zero ;
    xn1=zero; xn2=zero; rms1=zero; rms2=zero
    q21=zero; q22=zero; q41=zero; q42=zero
    tempE_Crho0=zero; tempREARR=zero
    DEROT=zero; SQUJ=zero; CRAN=zero; ERIGHFB=zero
    !------------------------------------------------
    ! zero optimization variables
    !------------------------------------------------
    If(DO_FITT) Then
       efit_0=zero; efitV0=zero; dfitV0=zero
       efit_rhorho=zero; efit_rhorhoD=zero;
       efit_rhotau=zero; efit_rhoDrho=zero;
       efit_rhonablaJ=zero; efit_JJ=zero;
    End If
    !------------------------------------------------
    ! Integration in coordinate space
    !------------------------------------------------
    Do ihli=1,nghl
       whl=wdcor(ihli)
       !------------------------------------------------
       ! np-representation
       !------------------------------------------------
       rn=ro(ihli,1);      rp=ro(ihli,2); rnp2=rn**2+rp**2; rnp1=rn-rp
       tnt=tau(ihli,1);    tpt=tau(ihli,2);
       dn=dro(ihli,1);     dp=dro(ihli,2);
       djn=dj(ihli,1);     djp=dj(ihli,2);
       akn=aka(ihli,1);    akp=aka(ihli,2)
       akn2=akn*akn;       akp2=akp*akp
       adn=akn*rn;         adp=akp*rp
       SFIZN=SFIZ(IHLI,1); SFIZP=SFIZ(IHLI,2);
       SFIRN=SFIR(IHLI,1); SFIRP=SFIR(IHLI,2);
       SZFIN=SZFI(IHLI,1); SZFIP=SZFI(IHLI,2);
       SRFIN=SRFI(IHLI,1); SRFIP=SRFI(IHLI,2);
       SNABLARN=NABLAR(IHLI,1); SNABLARP=NABLAR(IHLI,2);
       SNABLAZN=NABLAZ(IHLI,1); SNABLAZP=NABLAZ(IHLI,2);
       !------------------------------------------------
       ! t-representation
       !------------------------------------------------
       RHO_0=rn+rp;        RHO_1=rn-rp;
       TAU_0=tnt+tpt;      TAU_1=tnt-tpt;
       DRHO_0=dn+dp;       DRHO_1=dn-dp;
       DJ_0=djn+djp;       DJ_1=djn-djp;
       SFIZ_0=SFIZN+SFIZP; SFIZ_1=SFIZN-SFIZP;
       SFIR_0=SFIRN+SFIRP; SFIR_1=SFIRN-SFIRP;
       SZFI_0=SZFIN+SZFIP; SZFI_1=SZFIN-SZFIP;
       SRFI_0=SRFIN+SRFIP; SRFI_1=SRFIN-SRFIP;
       SNABLAR_0=SNABLARN+SNABLARP; SNABLAR_1=SNABLARN-SNABLARP;
       SNABLAZ_0=SNABLAZN+SNABLAZP; SNABLAZ_1=SNABLAZN-SNABLAZP;
       J2_0=SFIZ_0**2+SFIR_0**2+SZFI_0**2+SRFI_0**2
       J2_1=SFIZ_1**2+SFIR_1**2+SZFI_1**2+SRFI_1**2
       JabJba_0=2*(SFIZ_0*SZFI_0+SFIR_0*SRFI_0)
       JabJba_1=2*(SFIZ_1*SZFI_1+SFIR_1*SRFI_1)
       !
       Call calculate_U_parameters(RHO_0,RHO_1,TAU_0,TAU_1,DRHO_0,DRHO_1,  &
            (SNABLAR_0**2+SNABLAZ_0**2),(SNABLAR_1**2+SNABLAZ_1**2))
       !------------------------------------------------
       ! rms and deformations
       !------------------------------------------------
       z=fh(ihli); zz=z*z; rrr=zz+fl(ihli)**2
       p2=p32*zz-half*rrr; p3=p53*z*p2-p23*rrr*z; p4=p74*z*p3-p34*rrr*p2
       row=whl*rn; xn1=xn1+row; rms1=rms1+row*rrr; q21=q21+row*p2; q41=q41+row*p4
       row=whl*rp; xn2=xn2+row; rms2=rms2+row*rrr; q22=q22+row*p2; q42=q42+row*p4
       !------------------------------------------------
       ! PH energies
       !------------------------------------------------
       EKIN_N=EKIN_N+hb0n*(TAU_0+TAU_1)*HALF*whl*facECM                         ! kinetic, n
       EKIN_P=EKIN_P+hb0p*(TAU_0-TAU_1)*HALF*whl*facECM                         ! kinetic, p
       EVOL_rho_tau=EVOL_rho_tau+(Urhotau(0,0)*RHO_0*TAU_0  &                   ! volume rho tau
            +Urhotau(1,0)*RHO_1*TAU_1+Urhotau(2,0)*RHO_0*TAU_1  &
            +Urhotau(3,0)*RHO_1*TAU_0 )*whl
       EVOL_rho_rho=EVOL_rho_rho+(Urhorho(0,0)*RHO_0**2  &                      ! volume density dependent
            +Urhorho(1,0)*RHO_1**2+(Urhorho(3,0)+Urhorho(2,0))*RHO_0*RHO_1)*whl
       ESURF_rho_DELTA_rho =ESURF_rho_DELTA_rho+(UrhoDrho(0,0)*RHO_0*DRHO_0  &  ! surface: rho delta rho
            +UrhoDrho(1,0)*RHO_1*DRHO_1+UrhoDrho(2,0)*RHO_0*DRHO_1  &
            +UrhoDrho(3,0)*RHO_1*DRHO_0 )*whl
       ESURF_NABLA_rho_NABLA_rho=ESURF_NABLA_rho_NABLA_rho  &                   ! surface: (nabla rho)**2
            +(Unablarho(0,0)*(SNABLAR_0*SNABLAR_0+SNABLAZ_0*SNABLAZ_0)  &
            +Unablarho(1,0)*(SNABLAR_1*SNABLAR_1+SNABLAZ_1*SNABLAZ_1)  &
            +(Unablarho(3,0)+Unablarho(2,0))*(SNABLAR_0*SNABLAR_1+SNABLAZ_0*SNABLAZ_1) )*whl
       ESO_rho_NABLA_J=ESO_rho_NABLA_J+(UrhonablaJ(0,0)*RHO_0*DJ_0  &           ! spin-orbit rho Nabla . J
            +UrhonablaJ(1,0)*RHO_1*DJ_1+UrhonablaJ(2,0)*RHO_0*DJ_1  &
            +UrhonablaJ(3,0)*RHO_1*DJ_0 )*whl
       ESO_NABLA_rho_J=ESO_NABLA_rho_J  &
            +(UJnablarho(0,0)*(SNABLAR_0*(SFIZ_0-SZFI_0)-SNABLAZ_0*(SFIR_0-SRFI_0))  &  ! spin-orbit J . Nabla rho
            + UJnablarho(1,0)*(SNABLAR_1*(SFIZ_1-SZFI_1)-SNABLAZ_1*(SFIR_1-SRFI_1))  &
            + UJnablarho(2,0)*(SNABLAR_1*(SFIZ_0-SZFI_0)-SNABLAZ_1*(SFIR_0-SRFI_0))  &
            + UJnablarho(3,0)*(SNABLAR_0*(SFIZ_1-SZFI_1)-SNABLAZ_0*(SFIR_1-SRFI_1)) )*whl
       ETENS=ETENS+(UJJ(0,0)*J2_0+UJJ(1,0)*J2_1  &                              ! tensor J^2
            +(UJJ(3,0)+UJJ(2,0))*(SFIZ_0*SFIZ_1+SFIR_0*SFIR_1+SZFI_0*SZFI_1+SRFI_0*SRFI_1)&
            +UJabJba(0,0)*JabJba_0 + UJabJba(1,0)*JabJba_1 &
            +(UJabJba(3,0)+UJabJba(2,0))*(SFIZ_0*SZFI_1+SFIR_0*SRFI_1+SZFI_0*SFIZ_1+SRFI_0*SFIR_1))*whl
       EExtra=EExtra+(UEnonstdr(0)+UEnonstdr(1))*whl                            ! extra field if needed
       E_EXT_FIELD=E_EXT_FIELD + ( Vexternal(0,zero,fl(ihli),z)*RHO_0 &         ! external field
            +Vexternal(1,zero,fl(ihli),z)*RHO_1 )*whl
       !------------------------------------------------
       ! Coulomb & Hartree
       !------------------------------------------------
       If (icou.Ge.1) ecodi=ecodi+half*cou(ihli)*rp*whl
       If (icou.Eq.2.Or.icou.Eq.-4) ecoex=ecoex-CExPar*cex*rp**p43*whl ! Slater approximation
       E_HARTREE_DIR=E_HARTREE_DIR +half*vDHartree(ihli,1)*RHO_0*whl+half*vDHartree(ihli,2)*RHO_1*whl
       ! just for printing
       tempE_Crho0=tempE_Crho0+RHO_0**2*whl
       tempREARR=tempREARR+(Cdrho(0)*RHO_0**2+Cdrho(1)*RHO_1**2)*RHO_0**sigma*whl
       !------------------------------------------------
       ! pairing energy and delta
       !------------------------------------------------
       rsa0=(RHO_0/rho_c)
       If(pairing_regularization .and. All(geff_inv .Ne. 0.0_pr)) Then
          ept1=ept1+akn2*whl/geff_inv(ihli,1); del1=del1-adn*whl/geff_inv(ihli,1)
          ept2=ept2+akp2*whl/geff_inv(ihli,2); del2=del2-adp*whl/geff_inv(ihli,2)
       Else
          dd1n=CpV0(0)*(ONE-rsa0*CpV1(0))*whl
          dd1p=CpV0(1)*(ONE-rsa0*CpV1(1))*whl
          ept1=ept1+dd1n*akn2; del1=del1-dd1n*adn
          ept2=ept2+dd1p*akp2; del2=del2-dd1p*adp
       End If
       !------------------------------------------------
       ! optimization quantities
       !------------------------------------------------
       If(DO_FITT) Then
          efitV0(0)=efitV0(0)+(ONE-rsa0*CpV1(0))*akn2*whl
          efitV0(1)=efitV0(1)+(ONE-rsa0*CpV1(1))*akp2*whl
          dfitV0(0)=dfitV0(0)-(ONE-rsa0*CpV1(0))*adn*whl
          dfitV0(1)=dfitV0(1)-(ONE-rsa0*CpV1(1))*adp*whl
          !
          efit_rhotau(0)=efit_rhotau(0)+RHO_0*TAU_0*whl              ! rho tau
          efit_rhotau(1)=efit_rhotau(1)+RHO_1*TAU_1*whl              ! rho tau
          efit_rhorho(0)=efit_rhorho(0)+RHO_0**2*whl                 ! rho^2
          efit_rhorho(1)=efit_rhorho(1)+RHO_1**2*whl                 ! rho^2
          efit_rhorhoD(0)=efit_rhorhoD(0)+RHO_0**sigma*RHO_0**2*whl  ! rho^2
          efit_rhorhoD(1)=efit_rhorhoD(1)+RHO_0**sigma*RHO_1**2*whl  ! rho^2
          efit_rhoDrho(0)=efit_rhoDrho(0)+RHO_0*DRHO_0*whl           ! rho Delta rho
          efit_rhoDrho(1)=efit_rhoDrho(1)+RHO_1*DRHO_1*whl           ! rho Delta rho
          efit_rhonablaJ(0)=efit_rhonablaJ(0)+RHO_0*DJ_0*whl         ! rho nabla J J
          efit_rhonablaJ(1)=efit_rhonablaJ(1)+RHO_1*DJ_1*whl         ! rho nabla J J
          efit_JJ(0)=efit_JJ(0)+J2_0*whl                             ! J.J
          efit_JJ(1)=efit_JJ(1)+J2_1*whl                             ! J.J
       End If
    End Do !ihli
    !------------------------------------------------
    ! Finite range traces
    !------------------------------------------------
    if(finite_range) then
       e_gamma_fr_dir = trace_product(gamma_g_dir,rk)*0.5_pr
       e_gamma_fr_exc = trace_product(gamma_g_exc,rk)*0.5_pr
       call trace_product_2(delta_g_dir,rk,delrho1,delrho2)
       call trace_product_2(delta_g_dir,ak,delkap1,delkap2)
       delrho1 = -delrho1; delrho2 = -delrho2
    else
       e_gamma_fr_dir = zero; e_gamma_fr_exc = zero
       delrho1 = zero; delrho2 = zero
       delkap1 = zero; delkap2 = zero
    endif
    !------------------------------------------------
    ! Coulomb Gaussians traces
    !------------------------------------------------
    if(coulomb_gaussian) then
       if(icou.le.-1) ecodi = trace_product(coulf_g_dir,rk)*0.5_pr
       if(icou.eq.-2.or.icou.eq.-3) ecoex = trace_product(coulf_g_exc,rk)*0.5_pr
       call trace_product_2(coulf_d_dir,rk,cdelrho1,cdelrho2)
       call trace_product_2(coulf_d_dir,ak,cdelkap1,cdelkap2)
       cdelrho1 = -cdelrho1; cdelrho2 = -cdelrho2
    else
       cdelrho1 = zero; cdelrho2 = zero
       cdelkap1 = zero; cdelkap2 = zero
    endif
    !------------------------------------------------
    ! after the integration
    !------------------------------------------------
    xn(1)=xn1;              xn(2)=xn2;         xn(3)=xn1+xn2;
    rms(1)=rms1;            rms(2)=rms2
    q2(1)=q21;              q2(2)=q22;
    q4(1)=q41;              q4(2)=q42
    ekt(1)=EKIN_N;          ekt(2)=EKIN_P;     ekt(3)=ekt(1)+ekt(2)
    ept(1)=ept1;            ept(2)=ept2;       ept(3)=ept(1)+ept(2)
    frept(1)=delkap1;       frept(2)=delkap2;  frept(3)=frept(1)+frept(2)
    coept(1)=cdelkap1;      coept(2)=cdelkap2; coept(3)=coept(1)+coept(2)
    del(1)=del1/xn(1);      del(2)=del2/xn(2);
    frdel(1)=delrho1/xn(1); frdel(2)=delrho2/xn(2);
    codel(1)=cdelrho1/xn(1);codel(2)=cdelrho2/xn(2);
    !
    ept = ept + frept + coept
    del = del + frdel + codel
    !
    EVOL=EVOL_rho_tau+EVOL_rho_rho+E_HARTREE_DIR
    esurf=ESURF_rho_DELTA_rho+ESURF_NABLA_rho_NABLA_rho
    ESO=ESO_rho_NABLA_J+ESO_NABLA_rho_J
    ecoul=ecodi+ecoex
    etot=ekt(3)+evol+esurf+eso+ecoul+ept(3)+ETENS+EExtra+E_EXT_FIELD&
         + e_gamma_fr_dir+e_gamma_fr_exc
    ehfb=etot
    entropy(3)=entropy(1)+entropy(2)
    !------------------------------------------------
    ! rms and deformations
    !------------------------------------------------
    Do it=itmin,itmax
       rms(it)=Sqrt(rms(it)/xn(it))
       q2(it)=two*q2(it)       !Qnp=<2r^2P_2(teta)>=<2z^2-x^2-y^2>
       q4(it)=ffdef4*q4(it)    !Hn=8r^4P_4(teta)=8z^4-24z^2(x^2+y^2)+3(x^2+y^2)^2
       def(it)=Sqrt(pi/5.0_pr)*q2(it)/(rms(it)**2*xn(it))
    End Do
    r212=rms(1)**2; r222=rms(2)**2
    rms(3)=Sqrt((xn(1)*r212+xn(2)*r222)/amas)
    q2(3)=q2(1)+q2(2)          ! quadrupole moment
    q4(3)=q4(1)+q4(2)          ! hexadecapole moment
    def(3)=Sqrt(pi/5.0_pr)*q2(3)/(rms(3)**2*amas) !deformation
    bet=def(3)
    !bet=ffdef6*q2(3)/(amas*r02)  ! bet=Q2*Sqrt(5 Pi)/(3A x^2);  x=r0 A^(1/3)
    !------------------------------------------------
    ! Lipkin-Nogami energy
    !------------------------------------------------
    If(kindhfb.Lt.0) Then
       Call tracesln
       If(ierror_flag.Ne.0) Return
       etot=etot+etr(3)
    End If
    !------------------------------------------------
    ! optimization quantities
    !------------------------------------------------
    If(DO_FITT) Then
       efV_0=0.0_pr
       If(kindhfb.Lt.0) Then
          efV_0(0)=ala2(1)
          efV_0(1)=ala2(2)
       End If
       dfitV0(0)=dfitV0(0)/xn(1)
       dfitV0(1)=dfitV0(1)/xn(2)
       efit_0=etot-efitV0(0)*CpV0(0)-efitV0(1)*CpV0(1)  &
            -efit_rhotau(0)*Ctau(0)-efit_rhotau(1)*Ctau(1)  &
            -efit_rhorho(0)*Crho(0)-efit_rhorho(1)*Crho(1)  &
            -efit_rhorhoD(0)*Cdrho(0)-efit_rhorhoD(1)*Cdrho(1)  &
            -efit_rhoDrho(0)*CrDr(0)-efit_rhoDrho(1)*CrDr(1)  &
            -efit_rhonablaJ(0)*CrdJ(0)-efit_rhonablaJ(1)*CrdJ(1)  &
            -efit_JJ(0)*CJ(0)-efit_JJ(1)*CJ(1)
    End If
    !------------------------------------------------
    ! expectation values of multipole moments
    !------------------------------------------------
    Call moments_computeValue()
    !
    ! ATDHFB and GCM+GOA collective mass at the perturbative cranking approximation
    If (collective_inertia .And. lpr) Then
        Call calculate_collective_mass()
    End If
    !
    ! Fission fragment characteristics at convergence
    If (fission_fragments .And. lpr) Then
        ! Number of particles in the neck
        Call QNFIND()
        ! Position of the fragment centers of mass
        Call center_of_mass(Z_NECK,CENLEF,CENRIG)
        ! Occupation of quasiparticles
        If(qp_occupation) Call wave_localization(Z_NECK)
        ! Mass multipole moments in the fragment intrinsic frame
        If(.Not.Allocated(QLMTOT)) Allocate(QLMTOT(0:lambdaMax,0:1),QLMPRO(0:lambdaMax,0:1))
        I_TYPE=1
        Do LAMACT=0,lambdaMax
           Call QLMFRA(Z_NECK,LAMACT,QLMLEF,QLMRIG,CENLEF,CENRIG,I_TYPE)
           QLMTOT(LAMACT,0) = QLMLEF
           QLMTOT(LAMACT,1) = QLMRIG
        End Do
        ! Charge of the fission fragments
        I_TYPE=2
        Do LAMACT=0,lambdaMax
           Call QLMFRA(Z_NECK,LAMACT,QLMLEF,QLMRIG,CENLEF,CENRIG,I_TYPE)
           QLMPRO(LAMACT,0) = QLMLEF
           QLMPRO(LAMACT,1) = QLMRIG
        End Do
        ! Record occupation probabilities for canonical states and quasiparticles
        ! This information is used to
        If(statistical_PNP) Call sample_occupations()
    End If
    !------------------------------------------------
    ! debug
    !------------------------------------------------
    If(Print_Screen.And.debug_solver.Ge.2) Then
       If(do_print == 1) Then
          Write(lout,'(4(a12,g13.6))')  &
               ' Tn=     ',ekt(1),           ' Tp=     ',ekt(2), &
               ' EPn=    ',ept(1),           ' EPp=    ',ept(2),  &
               ' EVOL=   ',EVOL,             ' Esurf=  ',esurf,  &
               ' NrNr=   ',ESURF_NABLA_rho_NABLA_rho,' rDr=    ',ESURF_rho_DELTA_rho,  &
               ' Etens=  ',ETENS,            ' Eso=   ',eso,  &
               ' rNJ=    ',ESO_rho_NABLA_J,  ' NrJ=   ',ESO_NABLA_rho_J,  &
               ' ECd=    ',ecodi,            ' ECex=  ',ecoex, &
               ' EHd=    ',E_HARTREE_DIR,    ' Ir0^2= ',tempE_Crho0, &
               ' Eextra= ',EExtra,           ' Ext.Fl= ',E_EXT_FIELD, &
               ' Etot=  ',etot
          If(DO_FITT) Then
             Write(*,'(4(a12,g13.6))')
             Write(*,'(4(a12,g13.6))')  &
                  ' efrr0= ',efit_rhorho(0),     ' efrr1= ',efit_rhorho(1), &
                  ' efrrD0=   ',efit_rhorhoD(0),          ' efrr1D=  ',efit_rhorhoD(1),  &
                  ' efrt0= ',efit_rhotau(0),     ' efrt1= ',efit_rhotau(1), &
                  ' efrDr0=   ',efit_rhoDrho(0),          ' efrDr1=  ',efit_rhoDrho(1),  &
                  ' efrDj0=',efit_rhonablaJ(0),  ' efrDj1=',efit_rhonablaJ(1), &
                  ' efjj0=    ',efit_JJ(0),               ' efjj1=   ',efit_JJ(1),  &
                  ' efV0_0=',efitV0(0),          ' efV0_1=',efitV0(1), &
                  ' dfV0_0=   ',dfitV0(0),                ' dfV0_1=  ',dfitV0(1),  &
                  ' efV0=  ',efV_0(0),           ' efV_1= ',efV_0(1), &
                  ' ef0=      ',efit_0,                   ' etot=    ',etot
          End If
       End If
    End If
    !------------------------------------------------
    ! Part called at the very end only (lpr=T)
    !------------------------------------------------
    If (lpr) Then
       !------------------------------------------------
       ! other definitions of deformations  (ffdef6=Sqrt(5.0_pr*pi)/3.0_pr)
       !------------------------------------------------
       bet2(1)=ffdef6*q2(1)/(xn(1)*r02) ! beta_n=Qn*Sqrt(5 Pi)/(3N x^2)
       bet2(2)=ffdef6*q2(2)/(xn(2)*r02) ! x=r0 A^(1/3)
       bet2(3)=ffdef6*q2(3)/(amas*r02)
       het4(1)=ffdef7*q4(1)/(xn(1)*r04)
       het4(2)=ffdef7*q4(2)/(xn(2)*r04)
       het4(3)=ffdef7*q4(3)/(amas*r04)
       ! Charge radius, from Adv. Nucl. Phys. 8, 219 (1975)
       RpRMSsq=0.769_pr
       RnRMSsq=-0.1161_pr   ! J. Phys. G 33, 1 (2006)
       DarwinFoldy=0.033_pr ! Phys. Rev. A 56, 4579 (1997)
       rc = Sqrt(r222 + RpRMSsq + (xn(1)/xn(2))*RnRMSsq + DarwinFoldy)
       ! transitions to barn,barn^2,barn^4
       Do i=1,3
          q2(i)=q2(i)/100.0_pr; q4(i)=q4(i)/10000.0_pr
       End Do
       !------------------------------------------------
       ! STORE to unprojected buffer 'eresu'
       !------------------------------------------------
       ! ieresu=50 from module definitions
       ! ,'UEtot','Ubett','Ubetn','Ubetp',' UQt ',' UQn ',' UQp '  &
       eresu(1)=etot; eresu(2)=def(3); eresu(3)=def(1); eresu(4)=def(2);
       eresu(5)=q2(3); eresu(6)=q2(1); eresu(7)=q2(2);
       ! ,' Uln ',' Ulp ',' UpEn',' UpEp',' UpDn',' UpDp',' UAsn',' UAsp'  &
       eresu(8)=alast(1); eresu(9)=alast(2); eresu(10)=ept(1); eresu(11)=ept(2);
       eresu(12)=del(1); eresu(13)=del(2); eresu(14)=ass(1); eresu(15)=ass(2);
       ! ,' Urt ',' Urn ',' Urp ',' Urc ',' Uht ',' Uhn ',' Uhp '  &
       eresu(16)=rms(3); eresu(17)=rms(1); eresu(18)=rms(2); eresu(19)=rc;
       eresu(20)=het4(3); eresu(21)=het4(1); eresu(22)=het4(2);
       ! ,' Uqht',' Uqhn',' Uqhp'  &
       eresu(23)=q4(3); eresu(24)=q4(1); eresu(25)=q4(2);
       ! ,'UKINT','UKINN','UKINP',' USO ','UCDIR',' UCEX','UDisn','UDisp'  &
       eresu(26)=ekt(3); eresu(27)=ekt(1); eresu(28)=ekt(2); eresu(29)=eso;
       eresu(30)=ecodi; eresu(31)=ecoex; eresu(32)=Dispersion(1); eresu(33)=Dispersion(2);
       ! ,'UV2Mn','UV2Mp'
       eresu(34)=v2min(1); eresu(35)=v2min(2);
       !  ,'UECMT','UECMN','UECMP'
       eresu(36)=ECMHFB(3); eresu(37)=ECMHFB(1); eresu(38)=ECMHFB(2);
       !  ,'UROTT','UROTN','UROTP'
       eresu(39)=DEROT(3); eresu(40)=DEROT(1); eresu(41)=DEROT(2);
       !  ,'USQUJT','USQUJN','USQUJP'
       eresu(42)=SQUJ(3); eresu(43)=SQUJ(1); eresu(44)=SQUJ(2);
       !  ,'UCRANT','UCRANN','UCRANP'
       eresu(45)=CRAN(3); eresu(46)=CRAN(1); eresu(47)=CRAN(2);
       !  ,'UERIGT','UERIGN','UERIGP'
       eresu(48)=ERIGHFB(3); eresu(49)=ERIGHFB(1); eresu(50)=ERIGHFB(2);
       !
       ! nucleus with wrong assymptotic
       If(iasswrong(3).Ne.0) eresu(16)=-eresu(16)
       !------------------------------------------------
       ! WRITE UNPROJECTED OUTPUT
       !------------------------------------------------
       If(do_print == 1) Then
          Do iw=lout,lfile
             Write(iw,*)
             Write(iw,'(a,9x,a)')            '  NB! From expect (UNPROJECTED RESULTS)'
             Write(iw,*)
             If(iLST1.Ne.0)  &
                  Write(iw,'(a,3f15.6)') '  hfb decay const. ass ',ass
             Write(iw,'(a,5f15.6)') '  pairing: CpV0,CpV1,...    ',CpV0,CpV1
             Write(iw,'(a,a)')      '  forces:   ',skyrme
             If(keyblo(1).Ne.0) Write(iw,'(a,i4,a,f10.3)')  '  Blocked neutron block    ', bloblo(keyblo(1),1)
             If(keyblo(2).Ne.0) Write(iw,'(a,i4,a,f10.3)')  '  Blocked proton  block    ', bloblo(keyblo(2),2)
             Write(iw,*)
             Write(iw,'(/,28x,a,8x,a,9x,a)') ' neutrons ','protons','total'
             Write(iw,'(a,6f15.6)')          '  Requested part.numbs.',tz,Sum(tz)
             Write(iw,'(a,6f15.6)')          '  UnPj(av) part.numbs .',xn
             Write(iw,'(a,3f15.6)')          '  b0, bz, bp ..........',b0,bz,bp
             Write(iw,*)
             Write(iw,'(a,3f15.6)') '  lambda (ala) ........',ala
             Write(iw,'(a,3f15.6)') '  Lambda (alast) ......',alast
             Write(iw,'(a,3f15.6)') '  delta(n,p), pwi .....',del,pwi
             If(finite_range) then
                Write(iw,'(a,3f15.6)') '  pairing energy ......',frept(3)
             Else
                Write(iw,'(a,3f15.6)') '  pairing energy ......',ept
             End If
             If(kindhfb.Lt.0) Then
                Write(iw,'(a,3f15.6)') '  LN lambda_2 ... ala2 ',ala2
                Write(iw,'(a,3f15.6)') '  LN energies .........',etr
                Write(iw,'(a,3f15.6)') '  delta(n,p)+ala2 .....',del+ala2
                Write(iw,'(a,3f15.6)') '  Geff(n,p) ...........',Geff
             End If
             Write(iw,*)
             Write(iw,'(a,3f15.6)') '  rms-radius ..........',rms
             Write(iw,'(a,15x,2f15.6)') '  charge-radius, r0 ...',rc,r00
             Write(iw,'(a,3f15.6)') '  deformation beta2....',def
             Write(iw,'(a,3f15.6)') '  dipole moment[fm] ...',(qmoment(1,it),it=1,3)
             Write(iw,'(a,3f15.6)') '  quadrupole moment[b] ',(qmoment(2,it),it=1,3)
             Write(iw,'(a,3f15.6)') '  octupole moment .....',(qmoment(3,it),it=1,3)
             Write(iw,'(a,3f15.6)') '  hexadecapole moment .',(qmoment(4,it),it=1,3)
             Write(iw,'(a,3f15.6)') '  q5 ..................',(qmoment(5,it),it=1,3)
             Write(iw,'(a,3f15.6)') '  q6 ..................',(qmoment(6,it),it=1,3)
             Write(iw,'(a,3f15.6)') '  q7 ..................',(qmoment(7,it),it=1,3)
             Write(iw,'(a,3f15.6)') '  q8 ..................',(qmoment(8,it),it=1,3)
             Write(iw,*)
             Write(iw,'(a,3f15.6)')    '  kinetic energy ......',ekt
             Write(iw,'(a,30x,f15.6)') '  volume energy .......',evol
             Write(iw,'(a,30x,f15.6)') '        rho_tau .......',EVOL_rho_tau
             Write(iw,'(a,30x,f15.6)') '        rho_rho .......',EVOL_rho_rho
             Write(iw,'(a,30x,f15.6)') '  surface energy ......',esurf
             Write(iw,'(a,30x,f15.6)') '   rho_DELTA_rho ......',ESURF_rho_DELTA_rho
             Write(iw,'(a,30x,f15.6)') '   (NABLA_rho)^2 ......',ESURF_NABLA_rho_NABLA_rho
             Write(iw,'(a,30x,f15.6)') '  spin-orbit energy ...',eso
             Write(iw,'(a,30x,f15.6)') '        rho_NABLA_J ...',ESO_rho_NABLA_J
             Write(iw,'(a,30x,f15.6)') '        NABLA_rho_J ...',ESO_NABLA_rho_J
             if(finite_range) then
                Write(iw,'(a,30x,f15.6)') ' finite range direct...',e_gamma_fr_dir
                Write(iw,'(a,30x,f15.6)') ' finite range exchange.',e_gamma_fr_exc
                Write(iw,'(a,30x,f15.6)') ' finite range total....',e_gamma_fr_dir+e_gamma_fr_exc
                Write(iw,'(a,30x,f15.6)') ' finite range pairing..',frept(3)
             endif
             Write(iw,'(a,30x,f15.6)') '  coulomb energy ......',ecodi+ecoex+coept(3)
             Write(iw,'(a,30x,f15.6)') '          direct ......',ecodi
             Write(iw,'(a,30x,f15.6)') '          exchange ....',ecoex
             Write(iw,'(a,30x,f15.6)') '          pairing .....',coept(3)
             Write(iw,'(a,30x,f15.6)') '  tensor energy .......',etens
             Write(iw,'(a,30x,f15.6)') '  direct Hartree E  ...',E_HARTREE_DIR
             Write(iw,'(a,30x,f15.6)') '  Extra E .............',EEXTRA
             Write(iw,'(a,30x,f15.6)') '  External field E ....',E_EXT_FIELD
             Write(iw,'(a,3f15.6)')    '  Entropy .............',entropy
             Write(iw,*)
             Write(iw,'(a,30x,f15.6)')    '  tEnergy: ehfb (qp)...',ehfb
             If(kindhfb.Lt.0) Then
                Write(iw,'(a,30x,f15.6)') '  tEnergy: ehfb(qp)+LN ',etot
             End If
             Write(iw,*)
             Write(iw,'(a,6f15.6)')    '  Calculated but not added corrections'
             Write(iw,'(a,6f15.6)')    '  ===================================='
             Write(iw,'(a,6f15.6)')    '  cmc-diagonal part ...',ekt(1)/hb0n*hbzeron-ekt(1),&
                  ekt(2)/hb0p*hbzerop-ekt(2),ekt(1)/hb0n*hbzeron-ekt(1)+ekt(2)/hb0p*hbzerop-ekt(2)
             Write(iw,'(a,6f15.6)')    '  cmc-hfb .............',ECMHFB
             Write(iw,'(a,6f15.6)')    '  cranking rot corr ...',DEROT
             Write(iw,*)
             Write(iw,'(a,6f15.6)')    '  SQUJ ................',SQUJ
             Write(iw,'(a,6f15.6)')    '  CRAN x 4 ............',4.0_pr*CRAN
             Write(iw,'(a,6f15.6)')    '  Rigit Body ..........',ERIGHFB
             Write(iw,*)
             ! Print collective inertia mass tensor
             If(collective_inertia) Then
                Write(iw,'(a,6f15.6)')    '  ZPE (ATDHFB) ........',E0_ATD
                Write(iw,'(a,6f15.6)')    '  ZPE (GCM) ...........',E0_GCM
                Write(iw,*)
                Call print_collective(iw)
             End If
             ! Printing fision fragment characteristics
             If(fission_fragments) Call print_moments(iw)
          End Do
       End If
       !------------------------------------------------
       ! START corrected Lipkin-Nogami characteristics
       !------------------------------------------------
       If(kindhfb.Lt.0) Then
          Call densitln  !LN corrections to the densities
          If(ierror_flag.Ne.0) Return
          Do it=itmin,itmax
             xn(it)=zero
             rms(it)=zero; q2(it)=zero; q4(it)=zero
          End Do
          !
          Do ihli=1,nghl
             whl=wdcor(ihli)
             rn=ro(ihli,1); rp=ro(ihli,2); rnp2=rn**2+rp**2
             ! rms and deformations
             z=fh(ihli); zz=z*z; rrr=zz+fl(ihli)**2
             p2=p32*zz   -half*rrr
             p3=p53*z*p2 -p23*rrr*z
             p4=p74*z*p3 -p34*rrr*p2
             row=whl*rn
             xn(1)=xn(1)+row
             rms(1)=rms(1)+row*rrr
             q2(1)=q2(1)+row*p2
             q4(1)=q4(1)+row*p4
             row=whl*rp
             xn(2)=xn(2)+row
             rms(2)=rms(2)+row*rrr
             q2(2)=q2(2)+row*p2
             q4(2)=q4(2)+row*p4
          End Do !ihli
          !------------------------------------------------
          ! rms and deformations
          !------------------------------------------------
          Do it=itmin,itmax
             rms(it)=Sqrt(rms(it)/xn(it))
             q2(it)=two*q2(it)     !Qnp=<2r^2P_2(teta)>=<2z^2-x^2-y^2>
             q4(it)=ffdef4*q4(it)  !Hn=<8r^4P_4(teta)>=<8z^4-24z^2(x^2+y^2)+3(x^2+y^2)^2>
             def(it)=Sqrt(pi/5.0_pr)*q2(it)/(rms(it)**2*xn(it))
          End Do
          r212=rms(1)**2; r222=rms(2)**2
          rms(3)=Sqrt((xn(1)*r212+xn(2)*r222)/amas)
          q2(3)=q2(1)+q2(2)  ! quadrupole moment
          q4(3)=q4(1)+q4(2)  ! hexadecapole moment
          def(3)=Sqrt(pi/5.0_pr)*q2(3)/(rms(3)**2*amas) !deformation
          ! other definitions of the same quantitsies
          bet2(1)=ffdef6*q2(1)/(xn(1)*r02) !beta_n=Qn*Sqrt(5Pi)/(3N x^2)
          bet2(2)=ffdef6*q2(2)/(xn(2)*r02) !x=r0=1.2A^(1/3)
          bet2(3)=ffdef6*q2(3)/(amas*r02)
          het4(1)=ffdef7*q4(1)/(xn(1)*r04)
          het4(2)=ffdef7*q4(2)/(xn(2)*r04)
          het4(3)=ffdef7*q4(3)/(amas*r04)
          xn(3)=xn(1)+xn(2)
          bet=def(3)
          RpRMSsq=0.769_pr
          RnRMSsq=-0.1161_pr   ! J. Phys. G 33, 1 (2006)
          DarwinFoldy=0.033_pr ! Phys. Rev. A 56, 4579 (1997)
          rc = Sqrt(r222 + RpRMSsq + (xn(1)/xn(2))*RnRMSsq + DarwinFoldy)
          ! transitions to barn,barn^2,barn^4
          Do i=1,3
             q2(i)=q2(i)/100.0_pr; q4(i)=q4(i)/10000.0_pr
          End Do
          !------------------------------------------------
          ! STORE to unprojected LN buffer 'eresl'
          !------------------------------------------------
          ! ieresl=20 from module definitions
          ! ,' EHFBLN',' EHFB',' LNEt','LNbet','LNben','LNbep',' LNQt',' LNQn',' LNQp'  &
          eresl(1)=etot; eresl(2)=etot-etr(3);
          eresl(3)=def(3); eresl(4)=def(1); eresl(5)=def(2)
          eresl(6)=q2(3); eresl(7)=q2(1); eresl(8)=q2(2);
          ! ,'LNpEn','LNpEp','LNpDn','LNpDp',' LNrt',' LNrn',' LNrC'  &
          eresl(9)=ept(1); eresl(10)=ept(2); eresl(11)=del(1)+ala2(1); eresl(12)=del(2)+ala2(2);
          eresl(13)=rms(3); eresl(14)=rms(1); eresl(15)=rms(2); eresl(16)=rc;
          ! ,' LNam2n',' LNam2p',' LNe2n',' LNe2p'
          eresl(17)=ala2(1); eresl(18)=ala2(2); eresl(19)=etr(1); eresl(20)=etr(2)
          !------------------------------------------------
          ! WRITE UNPROJECTED LN OUTPUT
          !------------------------------------------------
          If(do_print == 1) Then
             Do iw=lout,lfile
                Write(iw,'(a,3f15.6)')
                Write(iw,'(a,3f15.6)') '  With Lipkin-Nogami Corrections'
                Write(iw,'(a,3f15.6)') '================================'
                Write(iw,'(a,3f15.6)') '  rms-radius ..........',rms
                Write(iw,'(a,15x,2f15.6)') '  charge-radius, r0 ...',rc,r00
                Write(iw,'(a,3f15.6)') '  deformation beta ....',def
                Write(iw,'(a,3f15.6)') '  quadrupole moment[b] ',q2
                Write(iw,'(a,3f15.6)') '  hexadecapole moment .',q4
                Write(iw,'(a,3f15.6)') '================================'
                Write(iw,'(a,3f15.6)')
             End Do
          End If
       End If
       !------------------------------------------------
       ! WRITE all blocking candidates
       !------------------------------------------------
       If(do_print == 1) Then
          Do iw=lout,lfile
             Write(iw,*)
             Do it=itmin,itmax
                If(it.Eq.1) Then
                   EqpPo=>REqpN; VqpPo=>RVqpN; UqpPo=>RUqpN; KpwiPo=>KpwiN; KqpPo=>KqpN
                Else
                   EqpPo=>REqpP; VqpPo=>RVqpP; UqpPo=>RUqpP; KpwiPo=>KpwiP; KqpPo=>KqpP
                End If
                !
                Write(iw,*) ' ',' Blocking candidates are:'
                Write(iw,*) '  ',protn(it),' eqpmin=',eqpmin(it),' pwiblo=',pwiblo
                Do i=1,blomax(it)
                   If(keyblo(3).Ne.0.and.i.Ne.keyblo(it)) cycle
                   Write(iw,'(a,i4,a,i4,a,i4,2x,i4,3(a,1x,f12.8,1x),a)') '    num=',i,  &
                        ' block=',bloblo(i,it),  &
                        ' state=',blo123(i,it),blok1k2(i,it),  &
                        ' Eqp=',EqpPo(KqpPo(blok1k2(i,it))),  &
                        ' (1-2N)E=',(one-two*uk(blok1k2(i,it),it))*EqpPo(KqpPo(blok1k2(i,it))),  &
                        ' overl=',vkmax(blok1k2(i,it),it),  &
                        tb(numax(blok1k2(i,it),it))
                End Do
                Write(iw,*)
             End Do
             Write(iw,*)
          End Do
       End If
       !
       !------------------------------------------------
       ! PAV
       !------------------------------------------------
       ! Projecting on different nucleus
       If(iproj.Ne.0) Then
          npr(1)=Int(npr1pj); npr(2)=Int(npr2pj)
          !tz(1)=Real(npr(1),Kind=pr); tz(2)=Real(npr(2),Kind=pr)
          Call expectpj(.True.)
       End If
    End If
    !
    If (debug_solver.Eq.1) Call get_CPU_time('expect',1)
    !
  End Subroutine expect
  !=======================================================================
  !> Determines if calculations include constraints on multipole moments
  !=======================================================================
  Subroutine Constraint_or_not(inin_INI0,inin0,icstr0)
    Implicit None
    Integer(ipr), Intent(in) :: inin_INI0
    Integer(ipr), Intent(inout) :: inin0,icstr0
    Integer(ipr) :: icount,l
    icount=0
    Do l=1,lambdaMax
       If(lambda_active(l).Gt.0) icount=icount+1
    End Do
    If(icount.Gt.0) Then
       icstr0=1; inin0=inin_INI0
    Else
       icstr0=0; inin0=inin_INI0
    End If
  End Subroutine Constraint_or_not
  !=======================================================================
  !> This routine updates the Lagrange multipliers of the multi-dimensional
  !> linear constraints based on the variation of the generalized density
  !> matrix and the QRPA matrix at the cranking approximation; see
  !> \cite berger1980 \cite younes2009-a \cite schunck2012
  !=======================================================================
  Subroutine getLagrange(ite)
    Implicit None
    Character(Len=1) :: trans
    Integer(ipr) :: ite,icons,lambda,icount,it,i,j,l,ierror,bb
    Integer(ipr) :: ib,nd,nd2,nhfb,i0,m,k1,k2,n1,n2,nd1,k,kk,kkk,ll
    Integer(ipr) :: i_uvN,i_uvP,incx,incy
    Integer(ipr), allocatable :: ipivot(:),iftN(:),iftP(:)
    Real(pr) :: minu,hla,t_term,temp_k,temp_l,result,brakev,epsilo
    Real(pr), Allocatable :: multMatElems(:)
    Real(pr), allocatable :: EqpN(:),EqpP(:)
    Real(pr), allocatable :: vecold(:),qmultt(:),veclam(:),veccns(:)
    Real(pr), allocatable :: cnsorg(:,:),cnsmat(:,:,:),cnsvec(:)
    Real(pr), allocatable :: fn12pl(:,:,:),fp12pl(:,:,:)
    Real(pr), allocatable :: fn11pl(:,:,:),fp11pl(:,:,:),fn11mi(:,:,:),fp11mi(:,:,:)
    Real(pr), allocatable :: doubln(:,:),doublp(:,:),dsum_n(:,:),dsum_p(:,:)
    Real(pr), allocatable :: workcn(:),dblmul(:,:),Umatr(:,:),Vmatr(:,:)
    !
    If (debug_solver.Eq.1) Call get_CPU_time('getLagrange',0)

    minu=-one
    epsilo=1.E-14
    !
    ! initializing the multipole moment template array
    Allocate(qmultt(0:lambdaMax));qmultt=zero
    Do lambda=0,lambdaMax
       qmultt(lambda)=qmoment(lambda,3)
    End Do
    !
    ! constructing the vector of the deviations of the current constraint
    ! from the requested values
    !
    Allocate(vecold(1:numberCons));vecold=zero
    Allocate(cnsvec(1:numberCons));cnsvec=zero
    Allocate(veclam(1:numberCons));veclam=zero
    !
    Do icons=1,numberCons
       lambda=multLambda(icons)
       ! regular multipole
       If(lambda.Ge.1) Then
          cnsvec(icons)=multRequested(lambda)-qmultt(lambda)
          If (nbroyden.lt.1) Then
              vecold(icons)=multLag(lambda)
          Else
              vecold(icons)=brin(nhhdim4+lambda)
          End If
       End If
       ! Gaussian neck
       If(lambda.Eq.0) Then
          cnsvec(icons)=neckRequested-neckValue
          If (nbroyden.lt.1) Then
              vecold(icons)=neckLag
          Else
              vecold(icons)=brin(nhhdim4+lambdaMax+1)
          End If
       End If
       veclam(icons)=vecold(icons)
    End Do
    !
    ! proceeding to determine the matrix of the constraint operators
    ! in the q.p. basis
    ! loop over the K blocks
    Allocate(cnsmat(numberCons,numberCons,0:2*nb-1));cnsmat=zero
    !
    i_uvN=0 ! new index referring to all q.p. vectors
    i_uvP=0 ! new index referring to all q.p. vectors
    !
!$OMP Parallel Do &
!$OMP& Default(None) &
!$OMP& SCHEDULE(Dynamic) &
!$OMP& SHARED(nb,id,ia,nbx,kd,ka,fn_T,KqpN,REqpN,RVqpN,RUqpN,KpwiN,fp_T,KqpP,REqpP,RVqpP,RUqpP,KpwiP, &
!$OMP&        numberCons,epsilo,temper,multLambda,minu,cnsmat,switch_on_temperature) &
!$OMP& PRIVATE(bb,ib,it,nd,nd2,nhfb,i0,m,Umatr,Vmatr,k,kk,kkk,n1,icons,lambda, &
!$OMP&         multMatElems,dblmul,j,n2,hla,gaussian_neck,l,ll,temp_k,temp_l,t_term,result, &
!$OMP&         doubln,fn12pl,fn11pl,EqpN,iftN,i_uvN,dsum_n,doublp,fp12pl,fp11pl,EqpP,iftP,i_uvP,dsum_p)
    Do bb=0,2*NB-1
       it = bb/NB + 1
       ib = Mod(bb,NB)+1
       !
       !------------------------------------------------------
       ! neutron sector
       !------------------------------------------------------
       !
       If(it.Eq.1) Then
          !------------------------------------------------------
          ! matrix of the constraint in q.p. basis
          !------------------------------------------------------
          nd=id(ib); nd2=nd*nd; nhfb=nd+nd; i0=ia(ib); m=ib+(it-1)*nbx
          !
          If(kd(ib,it).Gt.0) Then
             Allocate(doubln(nd,kd(ib,it))); doubln=zero
             Allocate(fn12pl(kd(ib,it),kd(ib,it),numberCons)); fn12pl=zero
             Allocate(Umatr(nd,kd(ib,it))); Umatr=zero
             Allocate(Vmatr(nd,kd(ib,it))); Vmatr=zero
             Allocate(EqpN(kd(ib,it))); EqpN=zero
             Allocate(ifTN(kd(ib,it))); ifTN=1
             ! temperature
             If(switch_on_temperature) Then
                Allocate(fn11pl(kd(ib,it),kd(ib,it),numberCons)); fn11pl=zero
             End If
             !
             ! U and V for this block (v. 101)
             Do k=1,kd(ib,it)
                ifTN(k)=ka(ib,it)+k; kk=KqpN(ka(ib,it)+k); EqpN(k)=REqpN(kk); kkk=KpwiN(ka(ib,it)+k)
                Do n1=1,nd
                   !i_uvN=i_uvN+1
                   i_uvN=kkk+n1
                   Vmatr(n1,k)=RVqpN(i_uvN)
                   Umatr(n1,k)=RUqpN(i_uvN)
                End Do
             End Do
             !
             Do icons=1,numberCons
                !
                lambda=multLambda(icons)
                If(lambda.Ge.1) Then
                   Allocate(multMatElems(1:nd2)); multMatElems=zero
                   Call moments_computeField(lambda,ib,multMatElems)
                   ! matrix of the constraints in HO basis (size nd x nd)
                   Allocate(dblmul(nd,nd));dblmul=zero
                   j=0
                   Do n1=1,nd
                      Do n2=1,n1
                         j=j+1;hla=multMatElems(j)
                         dblmul(n1,n2)=hla;dblmul(n2,n1)=hla
                      End Do
                   End Do
                End If
                If(lambda.Eq.0) Then
                   Allocate(gaussian_neck(1:nd2)); gaussian_neck=zero
                   Call neck_computeField(ib)
                   ! matrix of the constraints in HO basis (size nd x nd)
                   Allocate(dblmul(nd,nd));dblmul=zero
                   j=0
                   Do n1=1,nd
                      Do n2=1,n1
                         j=j+1;hla=gaussian_neck(j)
                         dblmul(n1,n2)=hla;dblmul(n2,n1)=hla
                      End Do
                   End Do
                End If
                !
                ! matrix of the constraint operator in the qp basis. due to
                ! the q.p. cut-off the actual size of the q.p. basis is not
                ! the same as the s.p. (ho) basis, and it is  not the  same
                ! for protons and neutrons.  the formulas implemented below
                ! differ from the 2 references for 3 reasons:
                !  - different  phase convention for the bogoliubov  matrix
                !  - block structure of the bogoliubov matrix in hfodd
                !  - storage in a() and b() arrays correspond to complex
                !    conjugate of the actual matrices
                !
                ! second term: v^{+} f^{*} u^{*} = v^{T} f u
                Call dgemm('n','n',nd,kd(ib,it),nd,one,dblmul,nd,Umatr,nd,zero,doubln,nd)
                Call dgemm('t','n',kd(ib,it),kd(ib,it),nd,one,Vmatr,nd,doubln,nd,zero,fn12pl(1,1,icons),kd(ib,it))
                !
                ! first term:  u^{+} f v^{*} = u^{T} f v
                Call dgemm('n','n',nd,kd(ib,it),nd,one,dblmul,nd,Vmatr,nd,zero,doubln,nd)
                Call dgemm('t','n',kd(ib,it),kd(ib,it),nd,one,Umatr,nd,doubln,nd,minu,fn12pl(1,1,icons),kd(ib,it))
                !
                ! temperature - computing \tilde{f}^{11}
                If(switch_on_temperature) Then
                  !
                  ! second term: v^{+} f^{*} v = v^{T} f v
                  Call dgemm('n','n',nd,kd(ib,it),nd,one,dblmul,nd,Vmatr,nd,zero,doubln,nd)
                  Call dgemm('t','n',kd(ib,it),kd(ib,it),nd,one,Vmatr,nd,doubln,nd,zero,fn11pl(1,1,icons),kd(ib,it))
                  ! first term:  u^{+} f u = u^{T} f u
                  Call dgemm('n','n',nd,kd(ib,it),nd,one,dblmul,nd,Umatr,nd,zero,doubln,nd)
                  Call dgemm('t','n',kd(ib,it),kd(ib,it),nd,one,Umatr,nd,doubln,nd,minu,fn11pl(1,1,icons),kd(ib,it))
                  !
                End If
                !
                If(Allocated(multMatElems)) Deallocate(multMatElems)
                If(Allocated(gaussian_neck)) Deallocate(gaussian_neck)
                Deallocate(dblmul)
                !
             End Do ! end icons (neutrons)
             Deallocate(doubln,Umatr,Vmatr)
             !------------------------------------------------------
             ! constraint correlation matrix
             !------------------------------------------------------
              Allocate(doubln(kd(ib,it),kd(ib,it))); doubln=zero
              Allocate(dsum_n(kd(ib,it),kd(ib,it))); dsum_n=zero
              !
              Do i=1,numberCons
                 Do j=1,numberCons
                    !
                    ! temperature
                    If((.Not.switch_on_temperature)) Then
                       !
                       Do l=1,kd(ib,it)
                          Do k=1,kd(ib,it)
                             If(Abs(EqpN(k)+EqpN(l)).Gt.Epsilo) Then
                                doubln(k,l)=fn12pl(k,l,i)/(EqpN(k)+EqpN(l))
                             Else
                                doubln(k,l)=zero
                             End If
                          End Do
                       End do
                       !
                       Call dgemm('t','n',kd(ib,it),kd(ib,it),kd(ib,it),one,fn12pl(1,1,j),kd(ib,it),&
                                                             doubln,kd(ib,it),zero,dsum_n,kd(ib,it))
                    Else
                       !
                       ! term corresponding to f^12
                       Do l=1,kd(ib,it)
                          Do k=1,kd(ib,it)
                             kk=iftN(k);ll=iftN(l)
                             temp_k=fn_T(kk)
                             temp_l=fn_T(ll)
                             If(Abs(EqpN(k)+EqpN(l)).Gt.Epsilo) Then
                                doubln(k,l)=fn12pl(k,l,i)*(one+temp_k+temp_l)/(EqpN(k)+EqpN(l))
                             Else
                                doubln(k,l)=zero
                             End If
                          End Do
                       End Do
                       !
                       Call dgemm('t','n',kd(ib,it),kd(ib,it),kd(ib,it),one,fn12pl(1,1,j),kd(ib,it),&
                                                             doubln,kd(ib,it),zero,dsum_n,kd(ib,it))
                       !
                       ! first term: positive simplex
                       Do l=1,kd(ib,it)
                          Do k=1,kd(ib,it)
                             kk=iftN(k);ll=iftN(l)
                             temp_k=fn_T(kk)
                             temp_l=fn_T(ll)
                             If(k.ne.l.And.(Abs(EqpN(k)-EqpN(l)).Gt.Epsilo)) Then
                                t_term=-(temp_k-temp_l)/(EqpN(k)-EqpN(l))
                             Else
                                t_term=-temp_k*(temp_k-one)/temper
                             End If
                             doubln(k,l)=half*t_term*fn11pl(k,l,i)
                          End Do
                       End Do
                       !
                       Call dgemm('t','n',kd(ib,it),kd(ib,it),kd(ib,it),one,fn11pl(1,1,j),kd(ib,it),&
                                                              doubln,kd(ib,it),one,dsum_n,kd(ib,it))
                       !
                    End If
                    !
                    ! taking the trace of the resulting matrix
                    !
                    result=zero
                    Do l=1,kd(ib,it)
                       result=result+dsum_n(l,l)
                    End Do
                    !
                    cnsmat(i,j,bb)=cnsmat(i,j,bb)+result
                    !
                 End Do ! end of loop over j constraint
              End Do ! end of loop over i constraint
              Deallocate(doubln,dsum_n,fn12pl,EqpN,ifTN)
              If(switch_on_temperature) Deallocate(fn11pl)
           End If
       !------------------------------------------------------
       ! Proton sector
       !------------------------------------------------------
       Else
          !------------------------------------------------------
          ! matrix of the constraint in q.p. basis
          !------------------------------------------------------
          nd=id(ib); nd2=nd*nd; nhfb=nd+nd; i0=ia(ib); m=ib+(it-1)*nbx
          !
          If(kd(ib,it).Gt.0) Then
             Allocate(doublp(nd,kd(ib,it))); doublp=zero
             Allocate(fp12pl(kd(ib,it),kd(ib,it),numberCons)); fp12pl=zero
             Allocate(Umatr(nd,kd(ib,it))); Umatr=zero
             Allocate(Vmatr(nd,kd(ib,it))); Vmatr=zero
             Allocate(EqpP(kd(ib,it))); EqpP=zero
             Allocate(ifTP(kd(ib,it))); ifTP=1
             ! temperature
             If(switch_on_temperature) Then
                Allocate(fp11pl(kd(ib,it),kd(ib,it),numberCons)); fp11pl=zero
             End If
             !
             ! U and V for this block
             Do k=1,kd(ib,it)
                ifTP(k)=ka(ib,it)+k; kk=KqpP(ka(ib,it)+k); EqpP(k)=REqpP(kk); kkk=KpwiP(ka(ib,it)+k)
                Do n1=1,nd
                   !i_uvP=i_uvP+1
                   i_uvP=kkk+n1
                   Vmatr(n1,k)=RVqpP(i_uvP)
                   Umatr(n1,k)=RUqpP(i_uvP)
                End Do
             End Do
             !
             Do icons=1,numberCons
                !
                lambda=multLambda(icons)
                If(lambda.Ge.1) Then
                   Allocate(multMatElems(1:nd2)); multMatElems=zero
                   Call moments_computeField(lambda,ib,multMatElems)
                   Allocate(dblmul(nd,nd));dblmul=zero
                   j=0
                   Do n1=1,nd
                      Do n2=1,n1
                         j=j+1;hla=multMatElems(j)
                         dblmul(n1,n2)=hla;dblmul(n2,n1)=hla
                      End Do
                   End Do
                End If
                If(lambda.Eq.0) Then
                   Allocate(gaussian_neck(1:nd2)); gaussian_neck=zero
                   Call neck_computeField(ib)
                   Allocate(dblmul(nd,nd));dblmul=zero
                   j=0
                   Do n1=1,nd
                      Do n2=1,n1
                         j=j+1;hla=gaussian_neck(j)
                         dblmul(n1,n2)=hla;dblmul(n2,n1)=hla
                      End Do
                   End Do
                End If
                !
                ! matrix of the constraint operator in the qp basis. due to
                ! the q.p. cut-off the actual size of the q.p. basis is not
                ! the same as the s.p. (ho) basis, and it is  not the  same
                ! for protons and neutrons.  the formulas implemented below
                ! differ from the 2 references for 3 reasons:
                !  - different  phase convention for the bogoliubov  matrix
                !  - block structure of the bogoliubov matrix in hfodd
                !  - storage in a() and b() arrays correspond to complex
                !    conjugate of the actual matrices
                !
                ! second term: v^{+} f^{*} u^{*} = v^{t} f u
                Call dgemm('n','n',nd,kd(ib,it),nd,one,dblmul,nd,Umatr,nd,zero,doublp,nd)
                Call dgemm('t','n',kd(ib,it),kd(ib,it),nd,one,Vmatr,nd,doublp,nd,zero,fp12pl(1,1,icons),kd(ib,it))
                !
                ! first term:  u^{+} f v^{*} = u^{t} f v
                Call dgemm('n','n',nd,kd(ib,it),nd,one,dblmul,nd,Vmatr,nd,zero,doublp,nd)
                Call dgemm('t','n',kd(ib,it),kd(ib,it),nd,one,Umatr,nd,doublp,nd,minu,fp12pl(1,1,icons),kd(ib,it))
                !
                ! temperature - computing \tilde{f}^{11}
                If(switch_on_temperature) Then
                  !
                  ! second term: v f^{*} v = v^{T} f v
                  Call dgemm('n','n',nd,kd(ib,it),nd,one,dblmul,nd,Vmatr,nd,zero,doublp,nd)
                  Call dgemm('t','n',kd(ib,it),kd(ib,it),nd,one,Vmatr,nd,doublp,nd,zero,fp11pl(1,1,icons),kd(ib,it))
                  ! first term:  u f u = u^{T} f u
                  Call dgemm('n','n',nd,kd(ib,it),nd,one,dblmul,nd,Umatr,nd,zero,doublp,nd)
                  Call dgemm('t','n',kd(ib,it),kd(ib,it),nd,one,Umatr,nd,doublp,nd,minu,fp11pl(1,1,icons),kd(ib,it))
                  !
                End If
                !
                If(Allocated(multMatElems)) Deallocate(multMatElems)
                If(Allocated(gaussian_neck)) Deallocate(gaussian_neck)
                Deallocate(dblmul)
                !
             End Do ! end icons (protons)
             Deallocate(doublp,Umatr,Vmatr)
             !------------------------------------------------------
             ! constraint correlation matrix
             !------------------------------------------------------
             Allocate(doublp(kd(ib,it),kd(ib,it))); doublp=zero
             Allocate(dsum_p(kd(ib,it),kd(ib,it))); dsum_p=zero
             !
             Do i=1,numberCons
                Do j=1,numberCons
                   !
                   ! temperature
                   If((.Not.switch_on_temperature)) Then
                      !
                      Do l=1,kd(ib,it)
                         Do k=1,kd(ib,it)
                            If(Abs(EqpP(k)+EqpP(l)).Gt.Epsilo) Then
                               doublp(k,l)=fp12pl(k,l,i)/(EqpP(k)+EqpP(l))
                            Else
                               doublp(k,l)=zero
                            End If
                         End Do
                      End do
                      !
                      Call dgemm('t','n',kd(ib,it),kd(ib,it),kd(ib,it),one,fp12pl(1,1,j),kd(ib,it),&
                                                            doublp,kd(ib,it),zero,dsum_p,kd(ib,it))
                   Else
                      !
                      ! term corresponding to f^12
                      Do l=1,kd(ib,it)
                         Do k=1,kd(ib,it)
                            kk=iftP(k);ll=iftP(l)
                            temp_k=fp_T(kk)
                            temp_l=fp_T(ll)
                            If(Abs(EqpP(k)+EqpP(l)).Gt.Epsilo) Then
                               doublp(k,l)=fp12pl(k,l,i)*(one+temp_k+temp_l)/(EqpP(k)+EqpP(l))
                            Else
                               doublp(k,l)=zero
                            End If
                         End Do
                      End Do
                      !
                      Call dgemm('t','n',kd(ib,it),kd(ib,it),kd(ib,it),one,fp12pl(1,1,j),kd(ib,it),&
                                                            doublp,kd(ib,it),zero,dsum_p,kd(ib,it))
                      !
                      ! first term: positive simplex
                      Do l=1,kd(ib,it)
                         Do k=1,kd(ib,it)
                            kk=iftP(k);ll=iftP(l)
                            temp_k=fp_T(kk)
                            temp_l=fp_T(ll)
                            If(k.ne.l.And.(Abs(EqpP(k)-EqpP(l)).Gt.Epsilo)) Then
                               t_term=-(temp_k-temp_l)/(EqpP(k)-EqpP(l))
                            Else
                               t_term=-temp_k*(temp_k-one)/temper
                            End If
                            doublp(k,l)=half*t_term*fp11pl(k,l,i)
                         End Do
                      End Do
                      !
                      Call dgemm('t','n',kd(ib,it),kd(ib,it),kd(ib,it),one,fp11pl(1,1,j),kd(ib,it),&
                                                             doublp,kd(ib,it),one,dsum_p,kd(ib,it))
                      !
                   End If
                   !
                   ! taking the trace of the resulting matrix
                   !
                   result=zero
                   Do l=1,kd(ib,it)
                      result=result+dsum_p(l,l)
                   End Do
                   !
                   cnsmat(i,j,bb)=cnsmat(i,j,bb)+result
                   !
                End Do ! end of loop over j constraint
             End Do ! end of loop over i constraint
             Deallocate(doublp,dsum_p,fp12pl,EqpP,ifTP)
             If(switch_on_temperature) Deallocate(fp11pl)
          End If
          !
       End If ! it
       !
    End Do ! end of loop over blocks ib
!$OMP End Parallel Do
    !
    ! computing the inverse of the correlation matrix
    Allocate(cnsorg(numberCons,numberCons));cnsorg=zero
    Do l=1,numberCons
       Do k=1,numberCons
          result = zero
          Do bb=0,2*nb-1
             result = result + cnsmat(k,l,bb)
          End Do
          cnsorg(k,l)=result
       End Do
    End Do
    !
    ierror=0
    Allocate(ipivot(numberCons))
    Call dgetrf(numberCons,numberCons,cnsorg,numberCons,ipivot,ierror)
    !
    ierror=0
    Allocate(workcn(numberCons))
    Call dgetri(numberCons,cnsorg,numberCons,ipivot,workcn,numberCons,ierror)
    Deallocate(ipivot)
    !
    ! constructing the vector of variations of the linear constraints
    trans='N'; incx=1; incy=1
    Call dgemv(trans,numberCons,numberCons,one,cnsorg,numberCons,cnsvec,incx,zero,workcn,incy)
    !
    ! updating the linear constraint vector (mixing has to be done simultaneously).
    If (ite.Eq.0) Then
        brakev=zero
    Else
        brakev=xmix
    End If
    !
    Allocate(veccns(numberCons))
    Do i=1,numberCons
       veccns(i)=veclam(i)+workcn(i)
       lambda=multLambda(i)
       ! regular multipole
       If(lambda.Ge.1) Then
          If(nbroyden.lt.1) Then
             multLag(lambda)=brakev*vecold(i)+(1.0-brakev)*veccns(i)
          Else
             multLag(lambda)=veccns(i)
             brout(nhhdim4+lambda)=multLag(lambda)
          End If
       End If
       If(lambda.Eq.0) Then
          If(nbroyden.lt.1) Then
             neckLag=brakev*vecold(i)+(1.0-brakev)*veccns(i)
          Else
             neckLag=veccns(i)
             brout(nhhdim4+lambdaMax+1)=neckLag
          End If
       End If
    End Do
    !
    Deallocate(veccns,vecold,workcn)
    Deallocate(cnsmat,cnsorg)
    Deallocate(qmultt,cnsvec,veclam)
    !
    If (debug_solver.Eq.1) Call get_CPU_time('getLagrange',1)
    !
    Return
  End Subroutine getLagrange
  !=======================================================================
  !> Search for the requested state to block
  !=======================================================================
  Subroutine requested_blocked_level(ib,it)
    Implicit None
    Integer(ipr), Intent(in)  :: ib,it
    Integer(ipr) :: nd,im,k,ndk,na2,nad2,iqn,k0,LAPLUS,OMEGA,n1,n2,n3
    Real(pr) :: s1,s2,UUk,VVk
    k0=0
    If(nkblo(it,2).Eq.0) Return
    If(Parity) Then
       LAPLUS=(ib+1)/2 !Yesp
    Else
       LAPLUS=ib       !Nop
    End If
    OMEGA=2*LAPLUS-1
    If(nkblo(it,1).Ne.OMEGA) Return
    nd=ID(ib); im=ia(ib);
    Do k=1,nd
       ndk=k+nd; s1=zero
       Do na2=1,nd
          nad2=na2+nd
          UUk=allhfb(ib)%arr(na2,ndk)
          VVk=allhfb(ib)%arr(nad2,ndk)
          s2=Max(s1,Abs(UUk),Abs(VVk))
          If(s2.Gt.s1) Then
             s1=s2
             iqn=na2+im  ! the position in [123] numbering
          End If
       End Do
       ! quantum numbers: Omega,P[n1,n2,n3]=>OMEGA,tpar(npar(iqn))[nz(iqn)+2*nr(iqn)+nl(iqn),nz(iqn),nl(iqn)]
       If(nkblo(it,2).Ne.tpar(npar(iqn))) Cycle
       n3=nl(iqn);          If(nkblo(it,5).Ne.n3) Cycle
       n2=nz(iqn);          If(nkblo(it,4).Ne.n2) Cycle
       n1=n2+2*nr(iqn)+n3;  If(nkblo(it,3).Ne.n1) Cycle
       k0=iqn
       keyblo(it)=1
       bloblo(keyblo(it),it)=ib
       blo123(keyblo(it),it)=k
       Exit
    End Do
  End Subroutine requested_blocked_level
  !=======================================================================
  !
  !=======================================================================
  Subroutine handle_blocking(irestart, iexit)
     Integer(ipr), Intent(Inout) :: irestart
     Integer(ipr), Intent(In)    :: iexit
     Integer(ipr) :: it, prim, sec
     Integer(ipr), Dimension(2) :: nloops=1

     ! Check the blocking mode
     !! If blocking active, make sure we initialize from an even-even core
     !! (for now, otherwise particle number handling will be wrong)
     !! If (npr_INI(1).Ne.2*(npr_INI(1)/2) .Or. npr_INI(2).Ne.2*(npr_INI(2)/2)) Then
     !!   ierror_flag=ierror_flag+1
     !!   ierror_info(ierror_flag)='STOP: Must request an even-even core for blocking.'
     !! End If
     Do it=1,2
        If(blocking_mode(it).Ne.0) Then
           ! Don't allow specific blocking and auto blocking at the same time (for now)
           If (blocking_mode(1).Eq. -blocking_mode(2)) Then
              ierror_flag=ierror_flag+1
              ierror_info(ierror_flag)='STOP: Invalid blocking mode (Auto+Specific).'
              Return
           End If
        End If
     End Do

     ! If even-even nucleus (no block or 1st iter of auto block), enforce
     ! blocking_never_done = True, no matter what was read from binary file
     If ((blocking_mode(1).Ne.1 .And. blocking_mode(2).Ne.1)) blocking_never_done = .True.
     ! If blocking and read from file on 1st iter, skip the "even core" loop (via below)
     If(irestart.Eq.0 .And. iexit.Eq.0) Then
        Do it=1,2
           If(blocking_mode(it).Ne.0) Then
              irestart=1
           End If
        End Do
     End If

     ! Primary/secondary index for looping over P or N or both
     prim=1; sec=1
     Do it=1,2
        If(blocking_mode(it).Ne.0) Then
           prim = it; sec = it+1
        End If
     End Do
     If(sec .Gt. 2) sec=1

     ! Even core, keep keyblo=0 if blocking, but loop
     If(irestart.Eq.0) Then
        Do it=1,2
           If(blocking_mode(it).Ne.0) Then
              irestart=1
           End If
        End Do
     ! 1st blocking iteration, initialize loop variables
     Else If(irestart.Eq.1) Then
        ! Don't print blocking candidates for odd solutions.
        ! 1) Theres no reason to block again on an odd number parity state
        ! 2) It prints the new (blocked) info for block+state of old candidates
        !    from even core, which is meaningless
        keyblo(3)=1
        Do it=1,2
           ! Auto block: keyblo=1 indicates start blocking, irestart/=0 to loop (unless nloops=1)
           ! Specific block: nloops=1, keyblo=0 to be set in requested_blocked_level
           If(blocking_mode(it).Eq.-1) Then
              keyblo(it)=1
              irestart=2
              nloops(it)=blomax(it)
           End If
        End Do
        If(nloops(1)*nloops(2).Eq.1) irestart = 0
        ! Sanity check
        If(nloops(1)*nloops(2).Le.0) Then
           ierror_flag=ierror_flag+1
           ierror_info(ierror_flag)='STOP: Invalid number of blocking candidates.'
           Return
        End If
     ! Rest of iterations for automatic blocking
     Else
        irestart=irestart+1
        keyblo(prim)=keyblo(prim)+1
        ! Odd-Odd
        If(blocking_mode(sec).Ne.0) then
           If(keyblo(prim).Gt.blomax(prim)) Then
              keyblo(prim) = 1
              keyblo(sec)  = keyblo(sec)+1
           End If
           If(keyblo(1).Eq.blomax(1).And.keyblo(2).Eq.blomax(2)) irestart =0
        ! Odd-Even or Even-Odd
        Else
           If(keyblo(prim).Eq.blomax(prim)) Then
              irestart =0
           End If
        End If
     End if
  End Subroutine handle_blocking

  !==============================================QQst
  ! lpfam subroutines
  !==============================================QQ

  subroutine lptest
     write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
     write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
     write(*,*) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
     write(*,*) pi, pi, pi
  end subroutine lptest

  !=================================================================================================================================
  Subroutine qrpa_calculate_timeodd_CC
    Use HFBTHO, only : Crho,Cdrho,Ctau,CrDr,CrdJ,CJ  ! basic coupling constants
    ! from skyrme force
    Cjcjc(0)     =  -Ctau(0) ;       Cjcjc(1) =  -Ctau(1)
    Csnablajc(0) =   Crdj(0) ;   Csnablajc(1) =   Crdj(1) 
    Css(0)       = -2.0_pr*Crho(0)/3.0_pr-Crho(1)
    Css(1)       = -Crho(0)/3.0_pr
    CDss(0)      = -2.0_pr*CDrho(0)/3.0_pr-CDrho(1)
    CDss(1)      = -CDrho(0)/3.0_pr
    CDeltass(0)  = (Ctau(0)+3.0_pr*Ctau(1)  -4.0_pr*(Crdr(0)+Crdr(1)))/8.0_pr   
    CDeltass(1)  = (3.0_pr*(Ctau(0)-Ctau(1))-4.0_pr*(Crdr(0)+Crdr(1)))/24.0_pr
    CsT(0)       = -CJ(0)
    CsT(1)       = -CJ(1)
    !
  End Subroutine qrpa_calculate_timeodd_CC
  !=================================================================================================================================
  !
  !   Modified Bessel I_n(x)
  ! Returns I_n(x) * Exp(-x)
  ! Uses bessel I_0(x) and I_1(x) from HFBTHO
  !=================================================================================================================================
  function modbesselI(n,x)
    use Bessik
    implicit none
    integer(ipr), intent(in) :: n
    Real(pr), intent(in) :: x
    real(pr) :: modbesselI 
    !
    Real(pr) :: I0,I1,I2,I3, eps
    Real(pr) :: I4, I5, I6
    !
    eps=Spacing(1.0_pr);
    !
    I0 = besei0(x)
    I1 = besei1(x)
   !I2 = I0 - 2.0_pr*I1/(x+eps)
   !I3 = I1 - 4.0_pr*I2/(x+eps)
    I2 = I0 - 2.0_pr*I1/(x)
    I3 = I1 - 4.0_pr*I2/(x)
    I4 = I2 - 6.0_pr*I3/(x)
    I5 = I3 - 8.0_pr*I4/(x)
    I6 = I4 -10.0_pr*I5/(x)


    !
    modbesselI = 0.0_pr
    select case(n)
      case(0)
        modbesselI = I0
      case(1)
        modbesselI = I1
      case(2) 
        modbesselI = I2
      case(3)
        modbesselI = I3
      case(4)
        modbesselI = I4
      case(5)
        modbesselI = I5
      case(6)
        modbesselI = I6
      case default
        STOP "modBesselI: n not in [0,6]!"
   end select
   return


  end function modbesselI
  !=================================================================================================================================
  !
  ! Direct Coulomb part. This is now done in a similar manner as in HFBTHO v200d
  !
  !=================================================================================================================================
  Subroutine qrpa_coulom
    !-------------------------------------------------------------------------------------------------------------------------------
    ! Coulom-field (direct part), Gogny prescription, complex numbers
    !-------------------------------------------------------------------------------------------------------------------------------
    Use HFBTHO
    Use bessik
    Implicit None
    Integer(ipr), Save :: i,j,k
    Real(pr), Save :: zd2,y1,y2,xx1,s1,vik,f,r,r1,fac1,fac2,rr2,z,z1,zd1,t,  &
                      bb,r2,r12,rrr,rz1,rz2,rrz1,rrz2,xx,rk1,rip1,rkp1,alpha,&
                      beta,xxx
   !Complex(pr) :: cone, czero
   !cone = 1.0_pr;  czero = 0.0_pr
    !
    Call get_CPU_time('qrpa_coulom',0)
    !
    If(icacou_qrpa.Eq.0 .AND. K_QRPA_RESPONCE .EQ. 0) Then

      icacou_qrpa = 1
      vc_c = vc       ! vc_c can be taken from the HFB part, since there too K=0

    else if (icacou_qrpa .EQ. 0 .AND. K_QRPA_RESPONCE .NE. 0) Then

      icacou_qrpa = 1
      !
      ! For parity-breaking shapes, the Coulomb potential was incorrectly
      ! calculated by assuming the two intervals [0,+\infty[ and ]-infty,0]
      ! were equivalent (see also below). This bug was corrected in version
      ! 139a
      If(Parity) Then
        fac1 = (-1)**K_QRPA_RESPONCE;  fac2 = one 
      Else
        fac1 = zero; fac2 = two
      End If
      bb=50.0_pr          ! Length scale L
      beta=2.00_pr
      alpha=one/beta
      f=chargee2/Sqrt(pi) ! e^2/Sqrt(pi)

      Do i=1,nghl
        r = fl(i); z = fh(i)
        Do k=1,i
           !
           r1 = fl(k); z1 = fh(k)
           rrr = two*r*r1; rr2 = (r - r1)**2
           ! z>0 part
           zd1 = (z - z1)**2
           rz1 = rr2 + zd1
           ! z<0 part
           zd2 = (z + z1)**2
           rz2 = rr2 + zd2
           ! Gauss-Legendre integration over u from 0 to D
           xx1=zero
           Do j=1,nleg
              xx=(one-xleg(j)**beta)**alpha ! change of variable to 0 <= u <= 1
              xxx=(one-xleg(j)**beta)**(alpha+one)
              y1=(xleg(j)/(bb*xx))**2 ! u^2
              s1=y1*rrr               ! 2 u^2 r r'
              !!y2=besei0(s1)           ! I0( 2 u^2 r r' ) * exp(-2 u^2 r r')
              y2=modbesselI(K_QRPA_RESPONCE,s1)
              xx1=xx1+fac2*wleg(j)*y2*(Exp(-rz1*y1) + fac1*Exp(-rz2*y1)) / xxx
           End Do
           vik=f*xx1/bb
           !
           vc_c(i,k)=vik*wdcor(k)  !wdcor=pi*wh*wl*bz*bp*bp
           vc_c(k,i)=vik*wdcor(i)  !wdcor=pi*wh*wl*bz*bp*bp
           !
        End Do  !k
      End Do  !i
    End If  !! else if (icacou_qrpa .EQ. 0 .AND. K_QRPA_RESPONCE .NE. 0)
    ! Calculation of the Coulomb field
    cou_c=zero
    Call zgemm('n','n',nghl,1,nghl,cone,vc_c,nghl,ro_c(:,2),nghl,czero,cou_c,nghl)
    !
   !If (IDEBUG.Eq.1) Call get_CPU_time('coulom',1)
  End Subroutine qrpa_coulom
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_DENSIT
    !---------------------------------------------------------------------
    ! local densities in coordinate space. 
    !---------------------------------------------------------------------
    Use HFBTHO, only : kindhfb
    Use HFBTHO, only : npr_ini
    Use HFBTHO, only : ia, ns
    Use HFBTHO, only : QHLA_opt,FI1R_opt,FI1Z_opt,FI2D_opt,y_opt   
    Use HFBTHO, only : cj
    Use HFBTHO, only : wdcori
    Use hfbtho_utilities, only : get_CPU_time
    Implicit None
    Real(pr)     :: y
    Complex(pr)  :: ss
    Complex(pr)  :: TEMP1_c,TEMP2_c,TEMP3_c,TEMP4_c,TEMP5_c,TEMP6_c,TEMP7_c,TEMP8_c,TEMP9_c,TEMP1_c0,TEMP1_c1
    Complex(pr)  :: TEMP10_c,TEMP11_c,TEMP12_c,TEMP13_c,TEMP14_c,TEMP15_c,TEMP16_c,TEMP17_c,TEMP18_c
    Complex(pr)  :: TEMP_rr,TEMP_rf,TEMP_rz,TEMP_fr,TEMP_ff,TEMP_fz,TEMP_zr,TEMP_zf,TEMP_zz
    Complex(pr)  :: FIA1FIB2
    !
    Integer(ipr) :: nsa,nla,ihil,n1,n2
    Integer(ipr) :: imena,imenb,im,ib,it,J,JA
    Integer(ipr), Allocatable :: LAM1(:), LAM2(:), SS1(:), SS2(:)
    integer(ipr) :: ovlap
    complex(pr)  :: ovlapr,ovlapphi,ovlapz
    !
    ! pointers
    Complex(pr),  Pointer :: TAKA_c(:),TBKA_c(:),TRO_c(:),TDJ_c(:),TTAU_c(:),TDRO_c(:)
    Complex(pr),  Pointer :: TJRR_c(:),TJRFI_c(:),TJRZ_c(:),TJFIR_c(:),TJFIFI_c(:),TJFIZ_c(:)
    Complex(pr),  Pointer :: TJZR_c(:), TJZFI_c(:),TJZZ_c(:)
    Complex(pr),  Pointer :: TNABLAR_c(:),TNABLAZ_c(:)
    !
    Complex(pr),  Pointer :: Psor_c(:),Psoz_c(:),Psofi_c(:)       !! spin
    Complex(pr),  Pointer :: Pdsor_c(:),Pdsoz_c(:),Pdsofi_c(:)    !! laplacian spin
    Complex(pr),  Pointer :: Prjr_c(:),Prjz_c(:),Prjfi_c(:)       !! current (=j)
    Complex(pr),  Pointer :: Ptsor_c(:),Ptsoz_c(:),Ptsofi_c(:)    !! spin kinetic
    Complex(pr),  Pointer :: Pcurjr_c(:),Pcurjz_c(:),Pcurjfi_c(:) !! curl current (= curl j)
    Complex(pr),  Pointer :: Pcursr_c(:),Pcursz_c(:),Pcursfi_c(:) !! curl s (= curl s)
    Complex(pr),  Pointer :: VqpPoA(:),UqpPoA(:),VqpPoB(:),UqpPoB(:)
    Integer(ipr), Pointer :: Kpwi_qrpa(:)
    Complex(pr),  Pointer :: Xqrpa(:),Yqrpa(:)
    Complex(pr),  Pointer :: P1qrpa(:),P2qrpa(:)
    Complex(pr),  Pointer :: Vqp(:),Uqp(:)
    ! helpers
    Integer(ipr) :: ipa,ipq,ipt,ibHFBi,ibHFBf,ibTHOi,ibTHOf,lda
    Real(pr)     :: SSU,SSD
    ! density matrix elements
    Complex(pr) :: rhme,kpme,kbme
    !
    Complex(pr)               :: iunit
    Complex(pr), Allocatable  :: FI_cA(:),FI_cB(:), FIR_cA(:),FIR_cB(:), FIZ_cA(:),FIZ_cB(:)
    Complex(pr), Allocatable  :: FID2_cA(:),FID2_cB(:)
    !
    Complex(pr), Allocatable  :: FIU_cA(:),FID_cA(:),FIUR_cA(:),FIDR_cA(:)
    Complex(pr), Allocatable  :: FIUD2N_cA(:),FIDD2N_cA(:),FIUZ_cA(:),FIDZ_cA(:)    
    Complex(pr), Allocatable  :: FIU_cB(:),FID_cB(:),FIUR_cB(:),FIDR_cB(:)
    Complex(pr), Allocatable  :: FIUD2N_cB(:),FIDD2N_cB(:),FIUZ_cB(:),FIDZ_cB(:)    
    !
    Parameter(iunit = Cmplx(0.0_pr,1.0_pr,kind=pr))    
    !
    Call get_CPU_time('qrpa_DENSIT',0)
    !
    Do it=1,2 !! flush also densities in the case of zero paricle number
       ! pointers
       If(it.Eq.1) Then
         Vqp=>VqpN; Uqp=>UqpN; Xqrpa=>XN_c; Yqrpa=>YN_c; P1qrpa=>P1N_c(:); P2qrpa=>P2N_c(:);
       Else
         Vqp=>VqpP; Uqp=>UqpP; Xqrpa=>XP_c; Yqrpa=>YP_c; P1qrpa=>P1P_c(:); P2qrpa=>P2P_c(:);
       Endif
       !
       !------------------------------------------------
       ! Calculation of density matrices.
       ! The use of HFB pairing window here is set to be 
       ! equivalent to the old K=0 FAM module. Therefore,
       ! in matrix product, only the other HFB matrix is
       ! modified
       !------------------------------------------------
       if(USE_HFB_PWI) then
         ! set aux arrays
         If(it.Eq.1) Then
           hfb_UA = UqpN ; hfb_VA = VqpN
         Else
           hfb_UA = UqpP ; hfb_VA = VqpP
         Endif
         ! remove states outside pairing window
         do ib=1,iHFB_NB
           im=0
           do ja=1,iHFB_ID(ib)
             do j=1,iHFB_ID(ib) 
               if(ja.gt.ipwi_hfb(ib,it)) then
                hfb_UA( iHFB_DS(ib)+im ) = zero;
                hfb_VA( iHFB_DS(ib)+im ) = zero;
               end if
               im=im+1
             end do !! components
           end do !! q.p. states
         end do !! blocks
         !
         ! density matrices with pairing window
         call qrpa_TriProd('n','n','c',Uqp,Xqrpa,hfb_Va,qrpa_aux)
         call qrpa_TriProd('n','t','c',hfb_Va,Yqrpa,Uqp,rhoetadm)
         rhoetadm = -rhoetadm -qrpa_aux
         !
         call qrpa_TriProd('n','t','c',Vqp,Yqrpa,hfb_Va,kapetadm)
         call qrpa_TriProd('n','n','c',hfb_Ua,Xqrpa,Uqp,qrpa_aux)
         kapetadm = -kapetadm +qrpa_aux
         call qrpa_TriProd('n','c','c',Vqp,Xqrpa,hfb_Va,kapetabdm)
         call qrpa_TriProd('n','*','c',hfb_Ua,Yqrpa,Uqp,qrpa_aux)
         kapetabdm = -kapetabdm +qrpa_aux
         !
         ! add efa contibutions
         if((it.eq.1 .and. ENABLE_N_EFA) .or. ((it.eq.2 .and. ENABLE_P_EFA)) ) then

           !! ivestigate here proper use of pairing window! Now pairing window is set for both matrices.
           !! with efa, this should not have impact since blocked states are always withing pairing window (?)
           !! with finite temperature, this may have impact...
           call qrpa_TriProd('n','n','c',hfb_Ua,P1qrpa,hfb_Ua,qrpa_aux)
           call qrpa_TriProd('n','t','c',hfb_Va,P2qrpa,hfb_Va,qrpa_bux)
           rhoetadm = rhoetadm +qrpa_aux -qrpa_bux

           call qrpa_TriProd('n','t','c',hfb_Va,P2qrpa,hfb_Ua,qrpa_aux)
           call qrpa_TriProd('n','n','c',hfb_Ua,P1qrpa,hfb_Va,qrpa_bux)
           kapetadm = kapetadm +qrpa_aux +qrpa_bux

           call qrpa_TriProd('n','*','c',hfb_Va,P2qrpa,hfb_Ua,qrpa_aux)
           call qrpa_TriProd('n','c','c',hfb_Ua,P1qrpa,hfb_Va,qrpa_bux)
           kapetabdm = kapetabdm +qrpa_aux +qrpa_bux

         end if

       else
         !
         ! density matrices without pairing window
         call qrpa_TriProd('n','n','c',Uqp,Xqrpa,Vqp,qrpa_aux)
         call qrpa_TriProd('n','t','c',Vqp,Yqrpa,Uqp,rhoetadm)
         rhoetadm = -rhoetadm -qrpa_aux
         !
         call qrpa_TriProd('n','t','c',Vqp,Yqrpa,Vqp,kapetadm)
         call qrpa_TriProd('n','n','c',Uqp,Xqrpa,Uqp,qrpa_aux)
         kapetadm = -kapetadm +qrpa_aux
         call qrpa_TriProd('n','c','c',Vqp,Xqrpa,Vqp,kapetabdm)
         call qrpa_TriProd('n','*','c',Uqp,Yqrpa,Uqp,qrpa_aux)
         kapetabdm = -kapetabdm +qrpa_aux

         if((it.eq.1 .and. ENABLE_N_EFA) .or. ((it.eq.2 .and. ENABLE_P_EFA)) ) then
           call qrpa_TriProd('n','n','c',Uqp,P1qrpa,Uqp,qrpa_aux)
           call qrpa_TriProd('n','t','c',Vqp,P2qrpa,Vqp,qrpa_bux)
           rhoetadm = rhoetadm +qrpa_aux -qrpa_bux

           call qrpa_TriProd('n','t','c',Vqp,P2qrpa,Uqp,qrpa_aux)
           call qrpa_TriProd('n','n','c',Uqp,P1qrpa,Vqp,qrpa_bux)
           kapetadm = kapetadm +qrpa_aux +qrpa_bux

           call qrpa_TriProd('n','*','c',Vqp,P2qrpa,Uqp,qrpa_aux)
           call qrpa_TriProd('n','c','c',Uqp,P1qrpa,Vqp,qrpa_bux)
           kapetabdm = kapetabdm +qrpa_aux +qrpa_bux
         end if

       end if !! USE_HFB_PWI
       !
       ! pointers to densities, time even
       TRO_c=>ro_c(:,it);         TTAU_c=>tau_c(:,it);       TDJ_c=>dj_c(:,it);     TDRO_c=>dro_c(:,it);  
       TJRR_c=>JRR_c(:,it);       TJRFI_c=>JRFI_c(:,it);     TJRZ_c=>JRZ_c(:,it);   TJFIR_c=>JFIR_c(:,it);
       TJFIFI_c=>JFIFI_c(:,it);   TJFIZ_c=>JFIZ_c(:,it);     TJZR_c=>JZR_c(:,it);   TJZFI_c=>JZFI_c(:,it);
       TJZZ_c=>JZZ_c(:,it);
       TNABLAR_c=>NABLAR_c(:,it); TNABLAZ_c=>NABLAZ_c(:,it); TAKA_c=>aka_c(:,it);   TBKA_c=>bka_c(:,it);     
       !
       ! pointers to densities, time odd
       Psor_c=>sor_c(:,it);       Psoz_c=>soz_c(:,it);       Psofi_c=>sofi_c(:,it);
       Pdsor_c=>dsor_c(:,it);     Pdsoz_c=>dsoz_c(:,it);     Pdsofi_c=>dsofi_c(:,it);
       Prjr_c=>rjr_c(:,it);       Prjz_c=>rjz_c(:,it);       Prjfi_c=>rjfi_c(:,it);
       Ptsor_c=>tsor_c(:,it);     Ptsoz_c=>tsoz_c(:,it);     Ptsofi_c=>tsofi_c(:,it);
       Pcurjr_c=>curjr_c(:,it);   Pcurjz_c=>curjz_c(:,it);   Pcurjfi_c=>curjfi_c(:,it);
       Pcursr_c=>cursr_c(:,it);   Pcursz_c=>cursz_c(:,it);   Pcursfi_c=>cursfi_c(:,it);	   
       ! ZERO THE DENSITIES
       Do ihil=1,nghl
          TRO_c(ihil)=ZERO;     TTAU_c(ihil)=ZERO;    TDJ_c(ihil)=ZERO;   TDRO_c(ihil)=ZERO;
          TJRR_c(ihil)=ZERO;    TJRFI_c(ihil)=ZERO;   TJRZ_c(ihil)=ZERO;  TJFIR_c(ihil)=ZERO;
          TJFIFI_c(ihil)=ZERO;  TJFIZ_c(ihil)=ZERO;   TJZR_c(ihil)=ZERO;  TJZFI_c(ihil)=ZERO;
          TJZZ_c(ihil)=ZERO;
          TNABLAR_c(ihil)=ZERO; TNABLAZ_c(ihil)=ZERO; TAKA_c(ihil)=ZERO;  TBKA_c(ihil)=ZERO; 
          Psor_c(ihil)=ZERO;    Psoz_c(ihil)=ZERO;    Psofi_c(ihil)=ZERO; 
          Pdsor_c(ihil)=ZERO;   Pdsoz_c(ihil)=ZERO;   Pdsofi_c(ihil)=ZERO; 
          Prjr_c(ihil)=ZERO;    Prjz_c(ihil)=ZERO;    Prjfi_c(ihil)=ZERO; 
          Ptsor_c(ihil)=ZERO;   Ptsoz_c(ihil)=ZERO;   Ptsofi_c(ihil)=ZERO;
          Pcurjr_c(ihil)=ZERO;  Pcurjz_c(ihil)=ZERO;  Pcurjfi_c(ihil)=ZERO;
          Pcursr_c(ihil)=ZERO;  Pcursz_c(ihil)=ZERO;  Pcursfi_c(ihil)=ZERO;
       End Do
       ! in the case of zero particle number, only flush densities
       If((npr_INI(1).Eq.0).And.(it.Eq.1)) Cycle
       If((npr_INI(2).Eq.0).And.(it.Eq.2)) Cycle
       !-----------------------------------------------
       ! SCAN OVER GAUSS INTEGRATION POINTS
       !-----------------------------------------------
       Do ihil=1,nghl
         !-----------------------------------------------
         ! SCAN OVER QRPA BLOCKS 
         !-----------------------------------------------
         Do ib=1,iQRPA_NB
           If(kindhfb.Lt.0) Then
             STOP "LN does not work yet!"
           End if
           ibHFBi = iQRPA_BHI(ib) ; ibHFBf = iQRPA_BHF(ib)
           ibTHOi = iQRPA_BLI(ib) ; ibTHOf = iQRPA_BLF(ib)
           ipa = iHFB_DS(ibHFBi) ; ipq = iQRPA_DS(ib) ; ipt = iQRPA_DS(iQRPA_TRANSPOSE(ib))
           !
           ! allocate arrays
           IMENA=iQRPA_IDF(ib) ; IMENB=iQRPA_IDI(ib)
           If(Allocated(FI_cA)) Deallocate(FI_cA,FIR_cA,FIZ_cA,FID2_cA)
           Allocate(FI_cA(IMENA),FIR_cA(IMENA),FIZ_cA(IMENA),FID2_cA(IMENA))
           If(Allocated(FI_cB)) Deallocate(FI_cB,FIR_cB,FIZ_cB,FID2_cB)
           Allocate(FI_cB(IMENB),FIR_cB(IMENB),FIZ_cB(IMENB),FID2_cB(IMENB))
           if(allocated(lam1)) Deallocate(lam1,ss1,lam2,ss2)
           allocate(LAM1(IMENA),SS1(IMENA),LAM2(IMENB),SS2(IMENB))
           !
           ! run trough the left basis vectors and fill the array
           Do N1=1,iQRPA_IDF(ib)
             IM=ia(iQRPA_BLF(ib)); JA=IM+N1; NSA=NS(JA); SSU=Max(NSA,0); SSD=Max(-NSA,0)
             FI_cA(N1)   = QHLA_opt(JA,ihil)
             FIR_cA(N1)  = FI1R_opt(JA,ihil)
             FIZ_cA(N1)  = FI1Z_opt(JA,ihil)
             FID2_cA(N1) = FI2D_opt(JA,ihil)
             SS1(N1) = NSA
             if(NSA.gt.0) then
               LAM1(N1) = (iQRPA_OMF(ib) -1)/2
             else
               LAM1(N1) = -(iQRPA_OMF(ib) +1)/2
             end if
             If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) LAM1(N1) = abs(LAM1(N1))    !! positive omega states have positive lambda
           End Do
           !
           ! run trough the right basis vectors and fill the array
           Do N1=1,iQRPA_IDI(ib)
             IM=ia(iQRPA_BLI(ib)); JA=IM+N1; NSA=NS(JA); SSU=Max(NSA,0); SSD=Max(-NSA,0)
             FI_cB(N1)   = QHLA_opt(JA,ihil)
             FIR_cB(N1)  = FI1R_opt(JA,ihil)
             FIZ_cB(N1)  = FI1Z_opt(JA,ihil)
             FID2_cB(N1) = FI2D_opt(JA,ihil)
             SS2(n1) = NSA
             if(NSA.gt.0) then
               LAM2(N1) = (iQRPA_OMI(ib) -1)/2
             else
               LAM2(N1) = -(iQRPA_OMI(ib) +1)/2
             end if
             If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) LAM2(N1) = abs(LAM2(N1))
           End Do

           y=y_opt(ihil);
           j=iQRPA_DS(ib)-1
           !
           Do N2=1,iQRPA_IDI(ib)
             Do N1=1,iQRPA_IDF(ib)
                !-----------------------------------------------
                ! go one step forward in density matrix arrays
                !-----------------------------------------------
                j=j+1  
                !-----------------------------------------------
                ! induced density matrix and kappa matrices
                !-----------------------------------------------
                rhme = rhoetadm(j) ; kpme = kapetadm(j) ; kbme = kapetabdm(j)  
                !-----------------------------------------------
                ! TIME-EVEN PART
                !-----------------------------------------------
                !
                ovlap = 0 ; 
                if( abs(LAM1(n1) + LAM2(n2)) .eq. abs(K_qrpa_responce) ) ovlap = 1
                if( abs(LAM1(n1) - LAM2(n2)) .eq. abs(K_qrpa_responce) ) ovlap = 1
                !
                ovlapphi = 0; ovlapr = 0; 
                ovlapz = -ovlap;
                ovlapr = (1-ovlap)*(-IUnit)
                ovlapphi = (1-ovlap)*(-1)
                FIA1FIB2 = FI_cA(N1)*FI_cB(N2)  !! \Phi_a(r,z) \Phi_b(r,z). Used to calculate many densities
                !
                ! define overlaps for normal HO basis. 
                If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) then
                   if(LAM1(n1) .eq. LAM2(n2)) then ; ovlap = 1 ; else ; ovlap = 0 ; end if
                   ovlapz = ovlap
                   ovlapr = (1-ovlap)
                   ovlapphi = (1-ovlap)*Iunit
                   if( SS1(n1).lt.0 ) ovlapz = ovlapz * (-1.0_pr)
                   if( SS2(n2).lt.0 ) ovlapphi = ovlapphi * (-1.0_pr)
                end if
                !
                !  ****************
                !   time-even part
                !  ****************
                !
                TEMP1_c = FIA1FIB2*ovlap
                !
                TEMP2_c = ( FIR_cA(N1)*FIR_cB(N2) +FIZ_cA(N1)*FIZ_cB(N2) &
                           +y*y*FIA1FIB2*LAM1(N1)*LAM2(N2) )*ovlap
                !
                TEMP3_c = ( FID2_cA(N1)*FI_cB(N2) +FI_cA(N1)*FID2_cB(N2) +two*FIR_cA(N1)*FIR_cB(N2) &
                           +two*FIZ_cA(N1)*FIZ_cB(N2) -y*y*FIA1FIB2*K_qrpa_responce**2 )*ovlap 
                !
                TEMP4_c = (-( IUnit*LAM2(N2)*y*FIR_cA(N1)*FI_cB(N2)  + IUnit*LAM1(N1)*y*FI_cA(N1)*FIR_cB(N2) )*ovlapz & 
                           +( -IUnit*LAM1(N1)*y*FI_cA(N1)*FIZ_cB(N2) + IUnit*LAM2(N2)*y*FIZ_cA(N1)*FI_cB(N2) )*ovlapr &
                           +( FIZ_cA(N1)*FIR_cB(N2)-FIR_cA(N1)*FIZ_cB(N2) )*ovlapphi )/IUnit
                !
                ! K=0 and !=0 tensor part
                TEMP_rf = -iunit*( FI_cA(N1)*FIR_cB(N2) -FIR_cA(N1)*FI_cB(N2) )*half*ovlapphi
                TEMP_fr = -iunit*( -iunit*(LAM1(N1)-LAM2(N2))*y*FIA1FIB2 )*half*ovlapr
                TEMP_fz = -iunit*( -iunit*(LAM1(N1)+LAM2(N2))*y*FIA1FIB2 )*half*ovlapz
                TEMP_zf = -iunit*( FI_cA(N1)*FIZ_cB(N2) -FIZ_cA(N1)*FI_cB(N2) )*half*ovlapphi
                !
                ! K.ne.0 tensor part
                TEMP_rr = -iunit*( FI_cA(N1)*FIR_cB(N2) -FIR_cA(N1)*FI_cB(N2) )*half*ovlapr
                TEMP_rz = -iunit*( FI_cA(N1)*FIR_cB(N2) -FIR_cA(N1)*FI_cB(N2) )*half*ovlapz
                TEMP_ff = -iunit*( -iunit*(LAM1(N1)-LAM2(N2))*y*FIA1FIB2 )*half*ovlapphi
                TEMP_zr = -iunit*( FI_cA(N1)*FIZ_cB(N2) -FIZ_cA(N1)*FI_cB(N2) )*half*ovlapr
                TEMP_zz = -iunit*( FI_cA(N1)*FIZ_cB(N2) -FIZ_cA(N1)*FI_cB(N2) )*half*ovlapz
                !
                ! these densities have (e^iKphi - e^-iKphi) structure
                ! here e^+iKphi part defines the density and sign
                if(LAM1(n1) - LAM2(n2) .gt. 0) then
                  TEMP_rz = -TEMP_rz; TEMP_zz = -TEMP_zz;
                end if
                if(LAM1(n1) + LAM2(n2) +1 .eq. -abs(K_qrpa_responce)) then
                  TEMP_rr = -TEMP_rr; TEMP_ff = -TEMP_ff; TEMP_zr = -TEMP_zr;
                end if
                ! put zero on those which don't contribute when K=0
                if(K_qrpa_responce.eq.0) then
                  TEMP_rr = zero; TEMP_rz = zero; TEMP_ff = zero; TEMP_zr = zero; TEMP_zz = zero; 
                end if
                !
                ! pairing
                TAKA_c(IHIL)=TAKA_c(IHIL)+TEMP1_c*kpme
                TBKA_c(ihil)=TBKA_c(ihil)+TEMP1_c*kbme
                !
                ! correct densities when using normal HO basis
                If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) then
                  TEMP4_c = (+( IUnit*LAM2(N2)*y*FIR_cA(N1)*FI_cB(N2) + IUnit*LAM1(N1)*y*FI_cA(N1)*FIR_cB(N2) )*ovlapz & 
                             -( IUnit*LAM1(N1)*y*FI_cA(N1)*FIZ_cB(N2) + IUnit*LAM2(N2)*y*FIZ_cA(N1)*FI_cB(N2) )*ovlapr &
                             +( FIZ_cA(N1)*FIR_cB(N2)-FIR_cA(N1)*FIZ_cB(N2) )*ovlapphi )/IUnit
                  !! correct here tensor terms too
                  if(abs(CJ(0)) .ge. spacing(1.0_pr) .or. abs(CJ(1)) .ge. spacing(1.0_pr)) stop "tensor term support missing for this case"
                end if
                !
                ! time-even matter densities
                TRO_c(IHIL)   = TRO_c(IHIL)   + TEMP1_c*rhme
                TTAU_c(IHIL)  = TTAU_c(IHIL)  + TEMP2_c*rhme
                TDRO_c(IHIL)  = TDRO_c(IHIL)  + TEMP3_c*rhme
                TDJ_c(IHIL)   = TDJ_c(IHIL)   + TEMP4_c*rhme
                TJRFI_c(ihil) = TJRFI_c(ihil) + TEMP_rf*rhme
                TJFIR_c(ihil) = TJFIR_c(ihil) + TEMP_fr*rhme
                TJFIZ_c(ihil) = TJFIZ_c(ihil) + TEMP_fz*rhme
                TJZFI_c(ihil) = TJZFI_c(ihil) + TEMP_zf*rhme
                !
                TJRR_c(ihil)  = TJRR_c(ihil)  + TEMP_rr*rhme
                TJRZ_c(ihil)  = TJRZ_c(ihil)  + TEMP_rz*rhme
                TJFIFI_c(ihil)= TJFIFI_c(ihil)+ TEMP_ff*rhme
                TJZR_c(ihil)  = TJZR_c(ihil)  + TEMP_zr*rhme
                TJZZ_c(ihil)  = TJZZ_c(ihil)  + TEMP_zz*rhme
                !  ***************
                !   time-odd part
                !  ***************
                !
                TEMP3_c = 0.0_pr; TEMP4_c = 0.0_pr; TEMP6_c = 0.0_pr; TEMP7_c = 0.0_pr; TEMP9_c = 0.0_pr; TEMP10_c= 0.0_pr; 
                TEMP11_c= 0.0_pr; TEMP12_c= 0.0_pr; TEMP13_c= 0.0_pr; TEMP14_c= 0.0_pr; TEMP15_c= 0.0_pr; TEMP16_c= 0.0_pr;
                TEMP17_c= 0.0_pr; TEMP18_c= 0.0_pr;  
                !
                TEMP1_c = -iunit*( FIR_cA(N1)*FI_cB(N2) - FI_cA(N1)*FIR_cB(N2) )*half*ovlap; !! j_r
                TEMP2_c = -iunit*( FIZ_cA(N1)*FI_cB(N2) - FI_cA(N1)*FIZ_cB(N2) )*half*ovlap; !! j_z
                !
                TEMP5_c = FIA1FIB2*ovlapphi                                                  !! s_phi
                !
                TEMP8_c = (FID2_cA(N1)*FI_cB(N2)+FI_cA(N1)*FID2_cB(N2)+two*FIR_cA(N1)*FIR_cB(N2)+two*FIZ_cA(N1)*FIZ_cB(N2) &
                          -y*y*FIA1FIB2*K_qrpa_responce**2 -y*y*FIA1FIB2 )*ovlapphi                  !! (Delta s)_phi
                !
                TEMP11_c = -iunit*( FIR_cA(N1)*FIZ_cB(N2) - FIZ_cA(N1)*FIR_cB(N2) )*two*half*ovlap   !! (curl j)_phi
                !
                TEMP13_c = -(FIZ_cA(N1)*FI_cB(N2) +FI_cA(N1)*FIZ_cB(N2))*ovlapphi                    !! (curl s)_r
                TEMP15_c = (y*FIA1FIB2 +FIR_cA(N1)*FI_cB(N2) +FI_cA(N1)*FIR_cB(N2))*ovlapphi         !! (curl s)_z
                !
                TEMP17_c = ( FIR_cA(N1)*FIR_cB(N2) +FIZ_cA(N1)*FIZ_cB(N2) &
                            -y*y*FIA1FIB2*LAM1(N1)*LAM2(N2) )*ovlapphi                               !! T_fi
                !
                ! these densities have (e^iKphi - e^-iKphi) structure
                ! here e^+iKphi part defines the density and sign
                if(LAM1(n1) + LAM2(n2) +1 .eq. abs(K_qrpa_responce)) then                
                  TEMP4_c = FIA1FIB2*ovlapr                                                 !! s_r
                  !
                  TEMP7_c = (FID2_cA(N1)*FI_cB(N2)+FI_cA(N1)*FID2_cB(N2)+two*FIR_cA(N1)*FIR_cB(N2)+two*FIZ_cA(N1)*FIZ_cB(N2) &
                            -y*y*FIA1FIB2*K_qrpa_responce**2 -y*y*FIA1FIB2 )*ovlapr         !! (Delta s)_r
                  TEMP8_c = TEMP8_c -two*IUnit*K_qrpa_responce*y*y*FIA1FIB2*ovlapr          !! (Delta s)_phi term: (2/r^2)(dv_r/dphi)
                  !
                  TEMP14_c = TEMP14_c + (FIZ_cA(N1)*FI_cB(N2)+FI_cA(N1)*FIZ_cB(N2))*ovlapr  !! (curl s)_phi
                  TEMP15_c = TEMP15_c + IUnit*K_qrpa_responce*y*FIA1FIB2*ovlapr             !! (curl s)_z
                  TEMP16_c = ( FIR_cA(N1)*FIR_cB(N2) +FIZ_cA(N1)*FIZ_cB(N2)     &
                              -y*y*FIA1FIB2*LAM1(N1)*LAM2(N2) )*ovlapr                      !! T_r
                end if
                if(LAM1(n1) + LAM2(n2) +1 .eq. -abs(K_qrpa_responce)) then                
                  TEMP4_c = -FIA1FIB2*ovlapr                                                !! s_r
                  !
                  TEMP7_c = -(FID2_cA(N1)*FI_cB(N2)+FI_cA(N1)*FID2_cB(N2)+two*FIR_cA(N1)*FIR_cB(N2)+two*FIZ_cA(N1)*FIZ_cB(N2) &
                            -y*y*FIA1FIB2*K_qrpa_responce**2 -y*y*FIA1FIB2 )*ovlapr         !! (Delta s)_r
                  TEMP8_c = TEMP8_c +two*IUnit*K_qrpa_responce*y*y*FIA1FIB2*ovlapr          !! (Delta s)_phi term: (2/r^2)(dv_r/dphi)
                  !
                  TEMP14_c = TEMP14_c - (FIZ_cA(N1)*FI_cB(N2)+FI_cA(N1)*FIZ_cB(N2))*ovlapr  !! (curl s)_phi
                  TEMP15_c = TEMP15_c - IUnit*K_qrpa_responce*y*FIA1FIB2*ovlapr             !! (curl s)_z 
                  TEMP16_c = -( FIR_cA(N1)*FIR_cB(N2) +FIZ_cA(N1)*FIZ_cB(N2)     &
                               -y*y*FIA1FIB2*LAM1(N1)*LAM2(N2) )*ovlapr                     !! T_r
                end if
                !
                TEMP7_c = TEMP7_c +two*IUnit*K_qrpa_responce*y*y*FIA1FIB2*ovlapphi !! (Delta s)_r term: (-2/r^2)(dv_phi/dphi)
                !
                ! these densities have (e^iKphi - e^-iKphi) structure
                ! here e^+iKphi part defines the density and sign
                if(LAM1(n1) - LAM2(n2) .gt. 0) then
                  TEMP3_c  = -half*y*FIA1FIB2*(LAM1(N1)+LAM2(N2))*ovlap !! j_phi
                  TEMP6_c  = -FIA1FIB2*ovlapz                           !! s_z
                  TEMP9_c  = -(FID2_cA(N1)*FI_cB(N2)+FI_cA(N1)*FID2_cB(N2)+two*FIR_cA(N1)*FIR_cB(N2)+two*FIZ_cA(N1)*FIZ_cB(N2) &
                             -y*y*FIA1FIB2*K_qrpa_responce**2 )*ovlapz  !! (Delta s)_z
                  !
                  TEMP10_c = TEMP10_c +half*y*(FIZ_cA(N1)*FI_cB(N2)+FI_cA(N1)*FIZ_cB(N2))*(LAM1(N1)+LAM2(N2))*ovlap     !! (curl j)_r
                  TEMP12_c = TEMP12_c -half*y*( +FIR_cA(N1)*FI_cB(N2) +FI_cA(N1)*FIR_cB(N2))*(LAM1(N1)+LAM2(N2))*ovlap  !! (curl j)_z
                  !
                  TEMP13_c = TEMP13_c + IUnit*K_qrpa_responce*y*FIA1FIB2*ovlapz             !! (curl s)_r
                  TEMP14_c = TEMP14_c + (FIR_cA(N1)*FI_cB(N2) +FI_cA(N1)*FIR_cB(N2))*ovlapz !! (curl s)_phi
                  TEMP18_c = -( FIR_cA(N1)*FIR_cB(N2) +FIZ_cA(N1)*FIZ_cB(N2)     &
                               +y*y*FIA1FIB2*LAM1(N1)*LAM2(N2) )*ovlapz                     !! T_z
                else 
                  TEMP3_c  = half*y*FIA1FIB2*(LAM1(N1)+LAM2(N2))*ovlap !! j_phi
                  TEMP6_c  = FIA1FIB2*ovlapz                           !! s_z
                  TEMP9_c  = (FID2_cA(N1)*FI_cB(N2)+FI_cA(N1)*FID2_cB(N2)+two*FIR_cA(N1)*FIR_cB(N2)+two*FIZ_cA(N1)*FIZ_cB(N2) &
                             -y*y*FIA1FIB2*K_qrpa_responce**2 )*ovlapz !! (Delta s)_z
                  !
                  TEMP10_c = TEMP10_c -half*y*(FIZ_cA(N1)*FI_cB(N2)+FI_cA(N1)*FIZ_cB(N2))*(LAM1(N1)+LAM2(N2))*ovlap    !! (curl j)_r
                  TEMP12_c = TEMP12_c +half*y*( +FIR_cA(N1)*FI_cB(N2) +FI_cA(N1)*FIR_cB(N2))*(LAM1(N1)+LAM2(N2))*ovlap !! (curl j)_z
                  !
                  TEMP13_c = TEMP13_c - IUnit*K_qrpa_responce*y*FIA1FIB2*ovlapz             !! (curl s)_r
                  TEMP14_c = TEMP14_c - (FIR_cA(N1)*FI_cB(N2) +FI_cA(N1)*FIR_cB(N2))*ovlapz !! (curl s)_phi
                  TEMP18_c = ( FIR_cA(N1)*FIR_cB(N2) +FIZ_cA(N1)*FIZ_cB(N2)     &
                              +y*y*FIA1FIB2*LAM1(N1)*LAM2(N2) )*ovlapz                      !! T_z
                end if
                !
                ! these densities have (e^iKphi - e^-iKphi) structure
                ! here e^+iKphi part defines the density and sign
                TEMP10_c = TEMP10_c -IUnit*K_qrpa_responce*y*(-iunit*( FIZ_cA(N1)*FI_cB(N2) - FI_cA(N1)*FIZ_cB(N2) ))*half*ovlap !! (curl j)_r
                TEMP12_c = TEMP12_c +IUnit*K_qrpa_responce*y*(-iunit*( FIR_cA(N1)*FI_cB(N2) - FI_cA(N1)*FIR_cB(N2) ))*half*ovlap !! (curl j)_z
                !
                ! correct time-odd densities when using normal HO basis
                If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) then
                 !! time-odd densities are already correct
                end if
                !
                ! time-odd matter densities
                Prjr_c(ihil)    = Prjr_c(ihil) +TEMP1_c*rhme         !! j_r
                Prjz_c(ihil)    = Prjz_c(ihil) +TEMP2_c*rhme         !! j_z
                Psofi_c(ihil)   = Psofi_c(ihil) +TEMP5_c*rhme        !! s_fi
                Pdsofi_c(ihil)  = Pdsofi_c(ihil) +TEMP8_c*rhme       !! (Delta s)_fi
                Pcurjfi_c(ihil) = Pcurjfi_c(ihil) -TEMP11_c*rhme     !! (curl j)_fi
                Pcursr_c(ihil)  = Pcursr_c(ihil) -TEMP13_c*rhme      !! (curl s)_r
                Pcursz_c(ihil)  = Pcursz_c(ihil) -TEMP15_c*rhme      !! (curl s)_z
                Ptsofi_c(ihil)  = Ptsofi_c(ihil) +TEMP17_c*rhme      !! T_fi
                if(K_qrpa_responce.ne.0) then 
                  Prjfi_c(ihil)   = Prjfi_c(ihil) +TEMP3_c*rhme      !! j_fi
                  Psor_c(ihil)    = Psor_c(ihil) +TEMP4_c*rhme       !! s_r
                  Psoz_c(ihil)    = Psoz_c(ihil) +TEMP6_c*rhme       !! s_z
                  Pdsor_c(ihil)   = Pdsor_c(ihil) +TEMP7_c*rhme      !! (Delta s)_r
                  Pdsoz_c(ihil)   = Pdsoz_c(ihil) +TEMP9_c*rhme      !! (Delta s)_z
                  Pcurjr_c(ihil)  = Pcurjr_c(ihil) -TEMP10_c*rhme    !! (curl j)_r
                  Pcurjz_c(ihil)  = Pcurjz_c(ihil) -TEMP12_c*rhme    !! (curl j)_z
                  Pcursfi_c(ihil) = Pcursfi_c(ihil) -TEMP14_c*rhme   !! (curl s)_fi
                  Ptsor_c(ihil)   = Ptsor_c(ihil) +TEMP16_c*rhme     !! T_r
                  Ptsoz_c(ihil)   = Ptsoz_c(ihil) +TEMP18_c*rhme     !! T_z
                end if
                !
             End Do
           End Do

         End Do !! qrpa blocks
         !
       End DO  !! ihil
       !----------------------------------------------------
       ! REMOVES INT.WEIGHTS AND MULTIPLIES BY THE JACOBIAN
       !----------------------------------------------------
       snz_qrpa(it)=Two*Sum(TRO_c); dsnz_qrpa(it)=Two*Sum(TDRO_c)
       Do ihil=1,nghl      
          ss=wdcori(ihil)
          if(K_QRPA_RESPONCE .ne. 0) ss = ss*half !! normalization
          !----------------------------------
          ! TIME-EVEN PART
          !----------------------------------
          TAKA_c(ihil)    = TAKA_c(ihil)*ss
          TBKA_c(ihil)    = TBKA_c(ihil)*ss
          !
          ss=ss*Two
          TJRR_c(ihil)     = TJRR_c(ihil)*ss     ! J_r,r
          TJRFI_c(ihil)    = TJRFI_c(ihil)*ss    ! J_r,fi
          TJRZ_c(ihil)     = TJRZ_c(ihil)*ss     ! J_r,z
          TJFIR_c(ihil)    = TJFIR_c(ihil)*ss    ! J_fi,r
          TJFIFI_c(ihil)   = TJFIFI_c(ihil)*ss   ! J_fi,fi
          TJFIZ_c(ihil)    = TJFIZ_c(ihil)*ss    ! J_fi,z
          TJZR_c(ihil)     = TJZR_c(ihil)*ss     ! J_z,r
          TJZFI_c(ihil)    = TJZFI_c(ihil)*ss    ! J_z,fi
          TJZZ_c(ihil)     = TJZZ_c(ihil)*ss     ! J_z,z
          !
          TRO_c(ihil)     = TRO_c(ihil)*ss
          TDRO_c(ihil)    = TDRO_c(ihil)*ss
          TTAU_c(ihil)    = TTAU_c(ihil)*ss
          TDJ_c(ihil)     = TDJ_c(ihil)*ss
          TNABLAR_c(IHIL) = TNABLAR_c(IHIL)*ss
          TNABLAZ_c(IHIL) = TNABLAZ_c(IHIL)*ss
          !----------------------------------
          ! TIME-ODD PART
          !----------------------------------		  
          Prjr_c(ihil)     = Prjr_c(ihil)*ss     ! j_r
          Prjfi_c(ihil)    = Prjfi_c(ihil)*ss    ! j_fi
          Prjz_c(ihil)     = Prjz_c(ihil)*ss     ! j_z
          Psor_c(ihil)     = Psor_c(ihil)*ss     ! s_r
          Psofi_c(ihil)    = Psofi_c(ihil)*ss    ! s_fi          
          Psoz_c(ihil)     = Psoz_c(ihil)*ss     ! s_z
          Pcursr_c(ihil)   = Pcursr_c(ihil)*ss   ! cur s
          Pcursfi_c(ihil)  = Pcursfi_c(ihil)*ss  ! cur s
          Pcursz_c(ihil)   = Pcursz_c(ihil)*ss   ! cur s
          Pdsor_c(ihil)    = Pdsor_c(ihil)*ss    ! (laplacian s)_r
          Pdsofi_c(ihil)   = Pdsofi_c(ihil)*ss   ! (laplacian s)_fi
          Pdsoz_c(ihil)    = Pdsoz_c(ihil)*ss    ! (laplacian s)_z
          Pcurjr_c(ihil)   = Pcurjr_c(ihil)*ss   ! (cur j)_r
          Pcurjfi_c(ihil)  = Pcurjfi_c(ihil)*ss  ! (cur j)_fi
          Pcurjz_c(ihil)   = Pcurjz_c(ihil)*ss   ! (cur j)_z
          Ptsor_c(ihil)    = Ptsor_c(ihil)*ss    ! T_r
          Ptsofi_c(ihil)   = Ptsofi_c(ihil)*ss   ! T_fi
          Ptsoz_c(ihil)    = Ptsoz_c(ihil)*ss    ! T_z
       Enddo
    End Do !it
    !
    Call get_CPU_time('qrpa_DENSIT',1)
    !
  End Subroutine qrpa_DENSIT
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_field
    !---------------------------------------------------------------------
    ! calculates fields in r-space form axially symmetric densities
    ! NB! includes the external field
    !---------------------------------------------------------------------
    Use HFBTHO, only : CExPar, coex
    Use HFBTHO, only : icou, itmin, itmax 
    Use HFBTHO, only : p23 
    Use HFBTHO, only : rho_c, sigma
    Use HFBTHO, only : Cdrho, Ctau, Crdr, Crdj, Cj, CpV0, CpV1, Crho 
    Use hfbtho_utilities, only : get_CPU_time
    Implicit None
    Integer(ipr)   :: iw,it,ita,ihli
    Real(pr)       :: eps
    !
    Complex(pr)    :: RHO_0,RHO_1,TAU_0,TAU_1,DRHO_0,DRHO_1,DJ_0,DJ_1,davn0,dbvn0
    Complex(pr)    :: JRR_0,JRR_1,JRFI_0,JRFI_1,JRZ_0,JRZ_1,JFIR_0,JFIR_1,JFIFI_0,JFIFI_1
    Complex(pr)    :: JFIZ_0,JFIZ_1,JZR_0,JZR_1,JZFI_0,JZFI_1,JZZ_0,JZZ_1
    Complex(pr)    :: SNABLARN,SNABLAZN,SNABLARP,SNABLAZP
    Complex(pr)    :: SNABLAR_0,SNABLAZ_0,SNABLAR_1,SNABLAZ_1
    Complex(pr)    :: J2_0,J2_1,rsa0
    ! time-even fields
    Complex(pr)    :: pUr(2),pUt(2),pUNr(2),pUNz(2),pUDr(2),pUDj(2)
    Complex(pr)    :: tUr(2),tUt(2),tUNr(2),tUNz(2),tUDr(2),tUDj(2)
    Complex(pr)    :: pUJrr(2),pUJrfi(2),pUJrz(2),pUJfir(2),pUJfifi(2),pUJfiz(2),pUJzr(2),pUJzfi(2),pUJzz(2)
    Complex(pr)    :: tUJrr(2),tUJrfi(2),tUJrz(2),tUJfir(2),tUJfifi(2),tUJfiz(2),tUJzr(2),tUJzfi(2),tUJzz(2)
    !
    Complex(pr)    :: RHO_0_HFB, RHO_1_HFB, TAU_0_HFB, TAU_1_HFB, KAP_P_HFB, KAP_N_HFB
    !
    Complex(pr)    :: SPR_0,SPR_1,SPZ_0,SPZ_1,SPFI_0,SPFI_1
    Complex(pr)    :: DSPR_0,DSPR_1,DSPZ_0,DSPZ_1,DSPFI_0,DSPFI_1
    Complex(pr)    :: CJCr_0,CJCr_1,CJCz_0,CJCz_1,CJCfi_0,CJCfi_1
    Complex(pr)    :: Csr_0,Csr_1,Csz_0,Csz_1,Csfi_0,Csfi_1
    Complex(pr)    :: JCR_0,JCR_1,JCZ_0,JCZ_1,JCFI_0,JCFI_1
    Complex(pr)    :: Tr_0,Tr_1,Tfi_0,Tfi_1,Tz_0,Tz_1
    ! time-odd fields
    Complex(pr)    :: pSir(2),pSiz(2),pSifi(2),pSidr(2),pSidz(2),pSidfi(2),pjcr(2),pjcz(2),pjcfi(2),crpj(2),czpj(2),cfipj(2)
    Complex(pr)    :: tSir(2),tSiz(2),tSifi(2),tSidr(2),tSidz(2),tSidfi(2),tjcr(2),tjcz(2),tjcfi(2),crtj(2),cztj(2),cfitj(2)
    Complex(pr)    :: crts(2),czts(2),cfits(2),crps(2),czps(2),cfips(2)
    complex(pr)    :: tTTr(2),tTTz(2),tTTfi(2)
    complex(pr)    :: pTTr(2),pTTz(2),pTTfi(2)
    !
    Call get_CPU_time('qrpa_field',0)
    !
    ! fields
    eps=Spacing(1.0_pr);
    Do ihli=1,nghl
       !------------------------------------------------------------------------------------
       ! TIME-EVEN PART
       !------------------------------------------------------------------------------------
       !
       ! form the isoscalar and -vector densities from the arrays containin p-n densities
       RHO_0     = ro_c(ihli,1)     + ro_c(ihli,2);         RHO_1     = ro_c(ihli,1)     - ro_c(ihli,2)
       TAU_0     = tau_c(ihli,1)    + tau_c(ihli,2);        TAU_1     = tau_c(ihli,1)    - tau_c(ihli,2)
       DRHO_0    = dro_c(ihli,1)    + dro_c(ihli,2);        DRHO_1    = dro_c(ihli,1)    - dro_c(ihli,2)
       DJ_0      = dj_c(ihli,1)     + dj_c(ihli,2);         DJ_1      = dj_c(ihli,1)     - dj_c(ihli,2)
       JRR_0     = Jrr_c(ihli,1)    + Jrr_c(ihli,2);        JRR_1     = Jrr_c(ihli,1)    - Jrr_c(ihli,2)
       JRFI_0    = Jrfi_c(ihli,1)   + Jrfi_c(ihli,2);       JRFI_1    = Jrfi_c(ihli,1)   - Jrfi_c(ihli,2)
       JRZ_0     = Jrz_c(ihli,1)    + Jrz_c(ihli,2);        JRZ_1     = Jrz_c(ihli,1)    - Jrz_c(ihli,2)
       JFIR_0    = Jfir_c(ihli,1)   + Jfir_c(ihli,2);       JFIR_1    = Jfir_c(ihli,1)   - Jfir_c(ihli,2)
       JFIFI_0   = Jfifi_c(ihli,1)  + Jfifi_c(ihli,2);      JFIFI_1   = Jfifi_c(ihli,1)  - Jfifi_c(ihli,2)
       JFIZ_0    = Jfiz_c(ihli,1)   + Jfiz_c(ihli,2);       JFIZ_1    = Jfiz_c(ihli,1)   - Jfiz_c(ihli,2)
       JZR_0     = Jzr_c(ihli,1)    + Jzr_c(ihli,2);        JZR_1     = Jzr_c(ihli,1)    - Jzr_c(ihli,2)
       JZFI_0    = Jzfi_c(ihli,1)   + Jzfi_c(ihli,2);       JZFI_1    = Jzfi_c(ihli,1)   - Jzfi_c(ihli,2)
       JZZ_0     = Jzz_c(ihli,1)    + Jzz_c(ihli,2);        JZZ_1     = Jzz_c(ihli,1)    - Jzz_c(ihli,2)
       SNABLAR_0 = NABLAR_c(ihli,1) + NABLAR_c(ihli,2);     SNABLAR_1 = NABLAR_c(ihli,1) - NABLAR_c(ihli,2)
       SNABLAZ_0 = NABLAZ_c(ihli,1) + NABLAZ_c(ihli,2);     SNABLAZ_1 = NABLAZ_c(ihli,1) - NABLAZ_c(ihli,2)
       !
       RHO_0_HFB = ro_c_hfb(ihli,1)  + ro_c_hfb(ihli,2)  ;  RHO_1_HFB = ro_c_hfb(ihli,1)  - ro_c_hfb(ihli,2)     ! HFB rho_0 and rho_1 density
       TAU_0_HFB = tau_c_hfb(ihli,1) + tau_c_hfb(ihli,2) ;  TAU_1_HFB = tau_c_hfb(ihli,1) - tau_c_hfb(ihli,2)    ! HFB tau_0 and tau_1 density
       KAP_N_HFB = kap_c_hfb(ihli,1)                     ;  KAP_P_HFB = kap_c_hfb(ihli,2)                        ! HFB pairing densities
       !    
       ! At the moment basic Skyrme EDF is assumed with no extra density dependence besides the usual ones
       !
       ! initialize to zero these isospin indexed fields
       tUr=zero; tUDr=zero; tUNr=zero; tUNz=zero; tUt=zero; tUDj=zero; 
       tUJrr=zero; tUJrfi=zero; tUJrz=zero; tUJfir=zero; tUJfifi=zero;
       tUJfiz=zero; tUJzr=zero; tUJzfi=zero; tUJzz=zero;
       !
       ! rho^2 terms. Here the fields have been explicitly linearized for density dependent part
       tUr(1)=tUr(1)+two*Crho(0)*RHO_0 +(2.0_pr+sigma)*(1.0_pr+sigma)*CDrho(0)*(RHO_0_HFB**sigma)*RHO_0 &
                    +CDrho(1)*sigma*two*RHO_1*RHO_1_HFB*((RHO_0_HFB+eps)**(sigma-1.0_pr)) &
                    +CDrho(1)*sigma*(sigma-1.0_pr)*RHO_0*(RHO_1_HFB**2)*((RHO_0_HFB+eps)**(sigma-2.0_pr)) 
       tUr(2)=tUr(2)+two*Crho(1)*RHO_1 +two*CDrho(1)*(RHO_0_HFB**sigma)*RHO_1 &
                    +two*CDrho(1)*sigma*RHO_0*RHO_1_HFB*((RHO_0_HFB+eps)**(sigma-1.0_pr))
       !
       ! rho tau
       tUr(1)=tUr(1)+Ctau(0)*TAU_0
       tUr(2)=tUr(2)+Ctau(1)*TAU_1
       tUt(1)=tUt(1)+Ctau(0)*RHO_0
       tUt(2)=tUt(2)+Ctau(1)*RHO_1
       !
       ! rho Delta rho
       tUr(1)=tUr(1)+Crdr(0)*DRHO_0
       tUr(2)=tUr(2)+Crdr(1)*DRHO_1
       tUDr(1)=tUDr(1)+Crdr(0)*RHO_0
       tUDr(2)=tUDr(2)+Crdr(1)*RHO_1
       !
       ! rho nabla J
       tUr(1)=tUr(1)+Crdj(0)*DJ_0
       tUr(2)=tUr(2)+Crdj(1)*DJ_1
       tUdj(1)=tUdj(1)+Crdj(0)*RHO_0
       tUdj(2)=tUdj(2)+Crdj(1)*RHO_1
       !
       ! J^2
       tUJrr(1)   = tUJrr(1)   + CJ(0)*two*Jrr_0;     tUJrr(2) = tUJrr(2)   + CJ(1)*two*Jrr_1;
       tUJrfi(1)  = tUJrfi(1)  + CJ(0)*two*Jrfi_0;   tUJrfi(2) = tUJrfi(2)  + CJ(1)*two*Jrfi_1;
       tUJrz(1)   = tUJrz(1)   + CJ(0)*two*Jrz_0;     tUJrz(2) = tUJrz(2)   + CJ(1)*two*Jrz_1;
       tUJfir(1)  = tUJfir(1)  + CJ(0)*two*Jfir_0;   tUJfir(2) = tUJfir(2)  + CJ(1)*two*Jfir_1;
       tUJfifi(1) = tUJfifi(1) + CJ(0)*two*Jfifi_0; tUJfifi(2) = tUJfifi(2) + CJ(1)*two*Jfifi_1;
       tUJfiz(1)  = tUJfiz(1)  + CJ(0)*two*Jfiz_0;   tUJfiz(2) = tUJfiz(2)  + CJ(1)*two*Jfiz_1;
       tUJzr(1)   = tUJzr(1)   + CJ(0)*two*Jzr_0;     tUJzr(2) = tUJzr(2)   + CJ(1)*two*Jzr_1;
       tUJzfi(1)  = tUJzfi(1)  + CJ(0)*two*Jzfi_0;   tUJzfi(2) = tUJzfi(2)  + CJ(1)*two*Jzfi_1;
       tUJzz(1)   = tUJzz(1)   + CJ(0)*two*Jzz_0;     tUJzz(2) = tUJzz(2)   + CJ(1)*two*Jzz_1;
       !
       ! Proton-neutron representation
       pUr(1)  =tUr(1)   + tUr(2);       pUr(2)  =tUr(1)   - tUr(2)
       pUt(1)  =tUt(1)   + tUt(2);       pUt(2)  =tUt(1)   - tUt(2)
       pUDr(1) =tUDr(1)  + tUDr(2);      pUDr(2) =tUDr(1)  - tUDr(2)
       pUdj(1) =tUdj(1)  + tUdj(2);      pUdj(2) =tUdj(1)  - tUdj(2)
       pUNr(1) =tUNr(1)  + tUNr(2);      pUNr(2) =tUNr(1)  - tUNr(2)
       pUNz(1) =tUNz(1)  + tUNz(2);      pUNz(2) =tUNz(1)  - tUNz(2)
       !
       pUJrr(1)   = tUJrr(1)   + tUJrr(2);    pUJrr(2)   = tUJrr(1)   - tUJrr(2);
       pUJrfi(1)  = tUJrfi(1)  + tUJrfi(2);   pUJrfi(2)  = tUJrfi(1)  - tUJrfi(2);
       pUJrz(1)   = tUJrz(1)   + tUJrz(2);    pUJrz(2)   = tUJrz(1)   - tUJrz(2);
       pUJfir(1)  = tUJfir(1)  + tUJfir(2);   pUJfir(2)  = tUJfir(1)  - tUJfir(2);
       pUJfifi(1) = tUJfifi(1) + tUJfifi(2);  pUJfifi(2) = tUJfifi(1) - tUJfifi(2);
       pUJfiz(1)  = tUJfiz(1)  + tUJfiz(2);   pUJfiz(2)  = tUJfiz(1)  - tUJfiz(2);
       pUJzr(1)   = tUJzr(1)   + tUJzr(2);    pUJzr(2)   = tUJzr(1)   - tUJzr(2);
       pUJzfi(1)  = tUJzfi(1)  + tUJzfi(2);   pUJzfi(2)  = tUJzfi(1)  - tUJzfi(2);
       pUJzz(1)   = tUJzz(1)   + tUJzz(2);    pUJzz(2)   = tUJzz(1)   - tUJzz(2);
       !
       Do it=itmin,itmax   !! loop over n  & p
          ita=3-it
          ! coulomb
          If(it.Eq.2) Then 
             If(icou.Ge.1) pUr(it)=pUr(it)+cou_c(ihli) 
             If(icou.Eq.2) pUr(it)=pUr(it)+CExPar*coex*ro_c(ihli,it) / ( 3.0_pr*(eps+ro_c_hfb(ihli,2)**p23))
          Endif
          !------------------------------------------------------------------------------------
          ! PAIRING 		  
          !------------------------------------------------------------------------------------
          !
          ! density independent part
          davn0 = aka_c(ihli,it)*CpV0(it-1)
          dbvn0 = bka_c(ihli,it)*CpV0(it-1)
          !
          ! density dependent part
          davn0 = davn0 - aka_c(ihli,it)*CpV0(it-1)*CpV1(it-1)*RHO_0_HFB/rho_c
          dbvn0 = dbvn0 - bka_c(ihli,it)*CpV0(it-1)*CpV1(it-1)*RHO_0_HFB/rho_c
          !
          ! rearrangement term. density dep. proportional to rho_n+rho_p -> the same term for protons and neutrons
          pUr(it)=pUr(it)-CpV0(0)*CpV1(0)*conjg(KAP_N_HFB*aka_c(ihli,1) +KAP_N_HFB*Conjg(bka_c(ihli,1)))/rho_c &
                         -CpV0(1)*CpV1(1)*conjg(KAP_P_HFB*aka_c(ihli,2) +KAP_P_HFB*Conjg(bka_c(ihli,2)))/rho_c
          !
          If(it.Eq.1) Then 
             davn_c(ihli)=davn0 -CpV0(it-1)*CpV1(it-1)*KAP_N_HFB*RHO_0/rho_c        !! neutron fields
             dbvn_c(ihli)=dbvn0 -CpV0(it-1)*CpV1(it-1)*KAP_N_HFB*Conjg(RHO_0)/rho_c !! neutron fields
          Else
             davp_c(ihli)=davn0 -CpV0(it-1)*CpV1(it-1)*KAP_P_HFB*RHO_0/rho_c        !! proton fields
             dbvp_c(ihli)=dbvn0 -CpV0(it-1)*CpV1(it-1)*KAP_P_HFB*Conjg(RHO_0)/rho_c !! proton fields
          End If
          !------------------------------------------------------------------------------------
       End Do !it
       vn_c(ihli)     = pUr(1);   vp_c(ihli)     = pUr(2)    !* RHO_ij  
       vhbn_c(ihli)   = pUt(1);   vhbp_c(ihli)   = pUt(2)    !* TAU_ij  
       vrn_c(ihli)    = pUNr(1);  vrp_c(ihli)    = pUNr(2)   !* NABLAr RHO__ij  
       vzn_c(ihli)    = pUNz(1);  vzp_c(ihli)    = pUNz(2)   !* NABLAz RHO__ij  
       vdn_c(ihli)    = pUDr(1);  vdp_c(ihli)    = pUDr(2)   !* DELTA RHO_ij  
       vsn_c(ihli)    = pUdj(1);  vsp_c(ihli)    = pUdj(2)   !* NABLA . J__ij  
       !
       vJrrn_c(ihli)   = pUJrr(1);   vJrrp_c(ihli)   = pUJrr(2);  
       vJrfin_c(ihli)  = pUJrfi(1);  vJrfip_c(ihli)  = pUJrfi(2);  
       vJrzn_c(ihli)   = pUJrz(1);   vJrzp_c(ihli)   = pUJrz(2);  
       vJfirn_c(ihli)  = pUJfir(1);  vJfirp_c(ihli)  = pUJfir(2);  
       vJfifin_c(ihli) = pUJfifi(1); vJfifip_c(ihli) = pUJfifi(2);  
       vJfizn_c(ihli)  = pUJfiz(1);  vJfizp_c(ihli)  = pUJfiz(2);  
       vJzrn_c(ihli)   = pUJzr(1);   vJzrp_c(ihli)   = pUJzr(2);  
       vJzfin_c(ihli)  = pUJzfi(1);  vJzfip_c(ihli)  = pUJzfi(2);  
       vJzzn_c(ihli)   = pUJzz(1);   vJzzp_c(ihli)   = pUJzz(2);  
       !
       !------------------------------------------------------------------------------------
       ! * TIME-ODD PART *
       !------------------------------------------------------------------------------------
       !
       ! T-odd densities
       SPR_0   = sor_c(ihli,1)    + sor_c(ihli,2);        SPR_1 = sor_c(ihli,1)    - sor_c(ihli,2)
       SPZ_0   = soz_c(ihli,1)    + soz_c(ihli,2);        SPZ_1 = soz_c(ihli,1)    - soz_c(ihli,2)
       SPFI_0  = sofi_c(ihli,1)   + sofi_c(ihli,2);      SPFI_1 = sofi_c(ihli,1)   - sofi_c(ihli,2)
       DSPR_0  = dsor_c(ihli,1)   + dsor_c(ihli,2);      DSPR_1 = dsor_c(ihli,1)   - dsor_c(ihli,2)
       DSPZ_0  = dsoz_c(ihli,1)   + dsoz_c(ihli,2);      DSPZ_1 = dsoz_c(ihli,1)   - dsoz_c(ihli,2)
       DSPFI_0 = dsofi_c(ihli,1)  + dsofi_c(ihli,2);    DSPFI_1 = dsofi_c(ihli,1)  - dsofi_c(ihli,2)
       JCR_0   = rjr_c(ihli,1)    + rjr_c(ihli,2);        JCR_1 = rjr_c(ihli,1)    - rjr_c(ihli,2)
       JCZ_0   = rjz_c(ihli,1)    + rjz_c(ihli,2);        JCZ_1 = rjz_c(ihli,1)    - rjz_c(ihli,2)
       JCFI_0  = rjfi_c(ihli,1)   + rjfi_c(ihli,2);      JCFI_1 = rjfi_c(ihli,1)   - rjfi_c(ihli,2)
       CJCr_0  = curjr_c(ihli,1)  + curjr_c(ihli,2);     CJCr_1 = curjr_c(ihli,1)  - curjr_c(ihli,2);
       CJCz_0  = curjz_c(ihli,1)  + curjz_c(ihli,2);     CJCz_1 = curjz_c(ihli,1)  - curjz_c(ihli,2);
       CJCfi_0 = curjfi_c(ihli,1) + curjfi_c(ihli,2);   CJCfi_1 = curjfi_c(ihli,1) - curjfi_c(ihli,2);
       Csr_0   = cursr_c(ihli,1)  + cursr_c(ihli,2);      Csr_1 = cursr_c(ihli,1)  - cursr_c(ihli,2);
       Csz_0   = cursz_c(ihli,1)  + cursz_c(ihli,2);      Csz_1 = cursz_c(ihli,1)  - cursz_c(ihli,2);
       Csfi_0  = cursfi_c(ihli,1) + cursfi_c(ihli,2);    Csfi_1 = cursfi_c(ihli,1) - cursfi_c(ihli,2);
       Tr_0    = tsor_c(ihli,1)   + tsor_c(ihli,2);        Tr_1 = tsor_c(ihli,1)   - tsor_c(ihli,2);
       Tfi_0   = tsofi_c(ihli,1)  + tsofi_c(ihli,2);      Tfi_1 = tsofi_c(ihli,1)  - tsofi_c(ihli,2);
       Tz_0    = tsoz_c(ihli,1)   + tsoz_c(ihli,2);        Tz_1 = tsoz_c(ihli,1)   - tsoz_c(ihli,2);
       !
       ! put zero on this and use it for density dependent parts
       tUr=zero 
       ! 
       ! initialize isospin indexed fields to zero
       tSir=zero; tSiz=zero;  tSifi=zero; tSidr=zero; tSidz=zero; tSidfi=zero 
       tjcr=zero; tjcz=zero;  tjcfi=zero; crtj=zero;  cztj=zero;  cfitj=zero  
       crts=zero; czts=zero;  cfits=zero;
       tTTr=zero; tTTfi=zero; tTTz=zero;
       !
       ! j^2 term  (will be integreted with j_ij)
       tjcr(1)=tjcr(1)+two*Cjcjc(0)*JCR_0
       tjcr(2)=tjcr(2)+two*Cjcjc(1)*JCR_1
       tjcfi(1)=tjcfi(1)+two*Cjcjc(0)*JCFI_0
       tjcfi(2)=tjcfi(2)+two*Cjcjc(1)*JCFI_1
       tjcz(1)=tjcz(1)+two*Cjcjc(0)*JCZ_0
       tjcz(2)=tjcz(2)+two*Cjcjc(1)*JCZ_1
       !
       ! C^s s^2 term (will be integreted with s_ij)
       tSir(1) =tSir(1) +two*Css(0)*SPR_0
       tSir(2) =tSir(2) +two*Css(1)*SPR_1
       tSifi(1)=tSifi(1)+two*Css(0)*SPFI_0
       tSifi(2)=tSifi(2)+two*Css(1)*SPFI_1
       tSiz(1) =tSiz(1) +two*Css(0)*SPZ_0
       tSiz(2) =tSiz(2) +two*Css(1)*SPZ_1
       !
       ! CD^s s^2 density dependent term (will be integreted with s_ij)
       tSir(1) =tSir(1) +two*CDss(0)*SPR_0*(RHO_0_HFB**sigma)
       tSir(2) =tSir(2) +two*CDss(1)*SPR_1*(RHO_0_HFB**sigma)
       tSifi(1)=tSifi(1)+two*CDss(0)*SPFI_0*(RHO_0_HFB**sigma)
       tSifi(2)=tSifi(2)+two*CDss(1)*SPFI_1*(RHO_0_HFB**sigma)
       tSiz(1) =tSiz(1) +two*CDss(0)*SPZ_0*(RHO_0_HFB**sigma)
       tSiz(2) =tSiz(2) +two*CDss(1)*SPZ_1*(RHO_0_HFB**sigma)
       !
       ! CD^s s^2 -rearrangment term not needed (s-density zero for HFB)
       !
       ! s.\delta s term (will be integreted with s_ij) 
       tSir(1) =tSir(1) +two*CDeltass(0)*DSPR_0
       tSir(2) =tSir(2) +two*CDeltass(1)*DSPR_1
       tSifi(1)=tSifi(1)+two*CDeltass(0)*DSPFI_0
       tSifi(2)=tSifi(2)+two*CDeltass(1)*DSPFI_1
       tSiz(1) =tSiz(1) +two*CDeltass(0)*DSPZ_0
       tSiz(2) =tSiz(2) +two*CDeltass(1)*DSPZ_1
       !
       ! s dot curl j (will be integreted with j_ij and s_ij) 
       tSir(1) =tSir(1) +Csnablajc(0)*CJCr_0
       tSir(2) =tSir(2) +Csnablajc(1)*CJCr_1
       tSifi(1)=tSifi(1)+Csnablajc(0)*CJCfi_0
       tSifi(2)=tSifi(2)+Csnablajc(1)*CJCfi_1
       tSiz(1) =tSiz(1) +Csnablajc(0)*CJCz_0
       tSiz(2) =tSiz(2) +Csnablajc(1)*CJCz_1
       !
       tjcr(1) =tjcr(1) +Csnablajc(0)*Csr_0
       tjcr(2) =tjcr(2) +Csnablajc(1)*Csr_1
       tjcfi(1)=tjcfi(1)+Csnablajc(0)*Csfi_0
       tjcfi(2)=tjcfi(2)+Csnablajc(1)*Csfi_1
       tjcz(1) =tjcz(1) +Csnablajc(0)*Csz_0
       tjcz(2) =tjcz(2) +Csnablajc(1)*Csz_1
       !       
       ! s dot T term (will be integreted with s_ij and T_ij) 
       tSir(1) =tSir(1) +CsT(0)*Tr_0
       tSir(2) =tSir(2) +CsT(1)*Tr_1
       tSifi(1)=tSifi(1)+CsT(0)*Tfi_0
       tSifi(2)=tSifi(2)+CsT(1)*Tfi_1
       tSiz(1) =tSiz(1) +CsT(0)*Tz_0
       tSiz(2) =tSiz(2) +CsT(1)*Tz_1
       tTTr(1) =tTTr(1) +CsT(0)*SPR_0
       tTTr(2) =tTTr(2) +CsT(1)*SPR_1
       tTTfi(1)=tTTfi(1)+CsT(0)*SPFI_0
       tTTfi(2)=tTTfi(2)+CsT(1)*SPFI_1
       tTTz(1) =tTTz(1) +CsT(0)*SPZ_0
       tTTz(2) =tTTz(2) +CsT(1)*SPZ_1
       !
       ! add here s dot F later
       !
       ! Proton-neutron representation. 
       ! From density dependent part
       pUr(1)  =tUr(1)+tUr(2);        pUr(2)  =tUr(1)  -tUr(2)
       !
       ! delta s
       pSir(1) =tSir(1)+tSir(2);      pSir(2) =tSir(1)-tSir(2)
       pSifi(1) =tSifi(1)+tSifi(2);   pSifi(2) =tSifi(1)-tSifi(2)
       pSiz(1) =tSiz(1)+tSiz(2);      pSiz(2) =tSiz(1)-tSiz(2)
       !
       ! delta j
       pjcr(1) =tjcr(1)+tjcr(2);      pjcr(2) =tjcr(1)-tjcr(2)
       pjcz(1) =tjcz(1)+tjcz(2);      pjcz(2) =tjcz(1)-tjcz(2)
       pjcfi(1)=tjcfi(1)+tjcfi(2);    pjcfi(2)=tjcfi(1)-tjcfi(2)  
       !	   
       ! delta T
       pTTr(1)=tTTr(1)+tTTr(2);       pTTr(2)=tTTr(1)-tTTr(2);
       pTTfi(1)=tTTfi(1)+tTTfi(2);    pTTfi(2)=tTTfi(1)-tTTfi(2);
       pTTz(1)=tTTz(1)+tTTz(2);       pTTz(2)=tTTz(1)-tTTz(2);
       !
       ! Add to the time-even field arrays
       vn_c(ihli)=vn_c(ihli)+pUr(1); vp_c(ihli)=vp_c(ihli)+pUr(2)        !* RHO_ij  
       !
       ! Set to time-odd field arrays
       vSIGRn_c(ihli)  = pSir(1)  ;  vSIGRp_c(ihli)  = pSir(2)           !* s_r_ij 
       vSIGZn_c(ihli)  = pSiz(1)  ;  vSIGZp_c(ihli)  = pSiz(2)           !* s_z_ij
       vSIGFIn_c(ihli) = pSifi(1) ;  vSIGFIp_c(ihli) = pSifi(2)          !* s_fi_ij
       !
       vjRn_c(ihli)    = pjcr(1)  ;  vjRp_c(ihli)    = pjcr(2)           !* j_r_ij
       vjZn_c(ihli)    = pjcz(1)  ;  vjZp_c(ihli)    = pjcz(2)           !* j_z_ij
       vjFIn_c(ihli)   = pjcfi(1) ;  vjFIp_c(ihli)   = pjcfi(2)          !* j_fi_ij
       !
       vTRn_c(ihli)    = pTTr(1)  ;  vTRp_c(ihli)    = pTTr(2);          !* T_r_ij
       vTZn_c(ihli)    = pTTz(1)  ;  vTZp_c(ihli)    = pTTz(2);          !* T_z_ij
       vTFIn_c(ihli)   = pTTfi(1) ;  vTFIp_c(ihli)   = pTTfi(2);         !* T_fi_ij
       !
       ! add here later F -terms
       !
    End Do !ihli
    !
    Call get_CPU_time('qrpa_field',1)
    !
  End Subroutine qrpa_field
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_gamdel
    !---------------------------------------------------------------------
    ! Complex qrpa matrices BhN/P_c and BdN/P_c in configurational space 
    !---------------------------------------------------------------------
  ! Use HFBTHO, only : Print_Screen, idebug
    Use HFBTHO, only : ia, ns
    Use HFBTHO, only : QHLA_opt,FI1R_opt,FI1Z_opt,FI2D_opt,y_opt  
    Use hfbtho_utilities, only : get_CPU_time
    Implicit None
    Integer(ipr)    :: i,j,ihil,nd,im,nsa,JA,N1,N2,ibq
    Integer(ipr)    :: dim1,dim2,ovlap
    complex(pr)     :: vh,vhb,vdel,vjr,vjz,vso
    complex(pr)     :: vnabr,vnabz,vsfiz,vszfi,vsfir,vsrfi,vspinfi,vspinr
    complex(pr)     :: vsob,vspinfib,vspinrb,vspinz,vspinzb,vjfi
    complex(pr)     :: vjrf,vjrfb,vjfr,vjfrb,vjfz,vjfzb,vjzf,vjzfb
    complex(pr)     :: vjrr,vjrrb,vjrz,vjrzb,vjff,vjffb,vjzr,vjzrb,vjzz,vjzzb
    complex(pr)     :: vttbase,vttr,vttfi,vttz, vttrb,vttfib,vttzb
    complex(pr)     :: ovlapr,ovlapphi,ovlapz
    !
    Integer(ipr), allocatable :: lam1(:),lam2(:),ss1(:),ss2(:)
    Real(pr), allocatable     :: FI1(:),FI2(:),FIR1(:),FIR2(:),FIZ1(:),FIZ2(:),FIDE1(:),FIDE2(:)
    !
    Real(pr) :: y !!,XLAPI,XLAMI,XLAPF,XLAMF, XLAPIY,XLAMIY,XLAPFY,XLAMFY
    !
    Complex(pr)  :: iunit,iunit2
    !
    Parameter(iunit = Cmplx(0.0_pr,1.0_pr,kind=pr), iunit2 = iunit/two)
    !
    Call get_CPU_time('qrpa_gamdel',0)
    !
    hN_c=zero; hP_c=zero; BhN_c=zero; BhP_c=zero; dN_c=zero; dP_c=zero; BdN_c=zero; BdP_c=zero
    !
    !----------------------------------------------
    ! START BLOCKS
    !----------------------------------------------
    Do ibq=1,iQRPA_NB    !! blocks
      !
      dim1 = iQRPA_IDF(ibq) ; dim2 = iQRPA_IDI(ibq)
      !
      if(allocated(FI1)) Deallocate(FI1,FI2,FIR1,FIR2,FIZ1,FIZ2,FIDE1,FIDE2)
      allocate(FI1(dim1),FI2(dim2),FIR1(dim1),FIR2(dim2),FIZ1(dim1),FIZ2(dim2),FIDE1(dim1),FIDE2(dim2))
      !
      if(allocated(lam1)) Deallocate(lam1,ss1,lam2,ss2)
      allocate(lam1(dim1),lam2(dim2),ss1(dim1),ss2(dim2))
      !
      Do ihil=1,nghl     !! mesh
        !
        ! run trough the left basis vectors and fill the array
        Do N1=1,iQRPA_IDF(ibq)
          IM=ia(iQRPA_BLF(ibq)); JA=IM+N1; NSA=NS(JA);
          if(NSA.gt.0) then
            LAM1(N1) = (iQRPA_OMF(ibq) -1)/2
          else
            LAM1(N1) = -(iQRPA_OMF(ibq) +1)/2
          end if
          If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) LAM1(N1) = abs(LAM1(N1))  !! normal HO, positive omega states have positive lambda
          FI1(N1)  = QHLA_opt(JA,ihil)
          FIR1(N1) = FI1R_opt(JA,ihil)
          FIZ1(N1) = FI1Z_opt(JA,ihil)
          FIDE1(N1)= FI2D_opt(JA,ihil) -FI1(N1)*((y_opt(ihil)*LAM1(N1))**2)
          SS1(N1) = NSA
        End Do
        !
        ! run trough the right basis vectors and fill the array
        Do N1=1,iQRPA_IDI(ibq)
          IM=ia(iQRPA_BLI(ibq)); JA=IM+N1; NSA=NS(JA);
          if(NSA.gt.0) then
            LAM2(N1) = (iQRPA_OMI(ibq) -1)/2
          else
            LAM2(N1) = -(iQRPA_OMI(ibq) +1)/2
          end if
          If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) LAM2(N1) = abs(LAM2(N1))  !! normal HO, positive omega states have positive lambda
          FI2(N1)  = QHLA_opt(JA,ihil)
          FIR2(N1) = FI1R_opt(JA,ihil)
          FIZ2(N1) = FI1Z_opt(JA,ihil)
          FIDE2(N1)=FI2D_opt(JA,ihil) -FI2(N1)*((y_opt(ihil)*LAM2(N1))**2)
          SS2(N1) = NSA
        End Do
        !
        y=y_opt(ihil);
        !
        j=iQRPA_DS(ibq)
        !
        Do N2=1,iQRPA_IDI(ibq) !! states in the block
          Do N1=1,iQRPA_IDF(ibq)   !! states in the block
             !
             ovlap = 0;
             if( abs(LAM1(n1) + LAM2(n2)) .eq. abs(K_qrpa_responce) ) ovlap = 1  !! matching lambdas
             if( abs(LAM1(n1) - LAM2(n2)) .eq. abs(K_qrpa_responce) ) ovlap = 1  !! 
             !
             ovlapphi = 0; ovlapr = 0; ovlapz = 0;
             ovlapz = -ovlap;
             ovlapr = (1-ovlap)*(-IUnit)
             ovlapphi = (1-ovlap)*(-1)
             !
             ! define overlaps for normal HO basis. 
             If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) then
                if(LAM1(n1) .eq. LAM2(n2)) then ; ovlap = 1 ; else ; ovlap = 0 ; end if
                ovlapz = ovlap
                ovlapr = (1-ovlap)
                ovlapphi = (-1)*(1-ovlap)*Iunit                    !! opposite sign here compared to DENSIT
                if( SS1(n1).lt.0 ) ovlapz = ovlapz * (-1.0_pr)
                if( SS2(n2).lt.0 ) ovlapphi = ovlapphi * (-1.0_pr)
             end if

             ! some aux. variables, time even
             vh = FI1(N1)*FI2(N2) *ovlap                                                                   !! RHO_ij
             vhb = ( FIR1(N1)*FIR2(N2) +FIZ1(N1)*FIZ2(N2) +LAM1(N1)*LAM2(N2)*y*y*FI1(N1)*FI2(N2) )*ovlap   !! TAU_ij
             vdel = ( FIDE1(N1)*FI2(N2) + FI1(N1)*FIDE2(N2))*ovlap +two*vhb                                !! Laplacian Rho_ij
             vso = (-( IUnit*LAM2(N2)*y*FIR1(N1)*FI2(N2)  + IUnit*LAM1(N1)*y*FI1(N1)*FIR2(N2) ) *ovlapz &  !! nabla.J_ij
                    -( -IUnit*LAM1(N1)*y*FI1(N1)*FIZ2(N2) + IUnit*LAM2(n2)*y*FIZ1(N1)*FI2(N2) ) *ovlapr  &
                    -( FIZ1(N1)*FIR2(N2) - FIR1(N1)*FIZ2(N2) )*ovlapphi )/IUnit
             vsob= (-( IUnit*LAM1(N1)*y*FIR2(N2)*FI1(N1)  + IUnit*LAM2(N2)*y*FI2(N2)*FIR1(N1) ) *ovlapz &   !! nabla.J_ij
                    -( -IUnit*LAM2(N2)*y*FI2(N2)*FIZ1(N1) + IUnit*LAM1(N1)*y*FIZ2(N2)*FI1(N1) ) *ovlapr  &  !! for simplex -i states
                    -( FIZ2(N2)*FIR1(N1) - FIR2(N2)*FIZ1(N1) )*ovlapphi )/IUnit
             vjrf=-iunit*( FIR1(N1)*FI2(N2) -FI1(N1)*FIR2(N2) )*half*ovlapphi             !! J_r,phi_ij
             vjfr=-iunit*( iunit*(LAM1(N1)-LAM2(N2))*y*FI1(N1)*FI2(N2) )*half*ovlapr      !! J_phi,r_ij
             vjfz=-iunit*( -iunit*(LAM1(N1)+LAM2(N2))*y*FI1(N1)*FI2(N2) )*half*ovlapz     !! J_phi,z_ij
             vjzf=-iunit*( FIZ1(N1)*FI2(N2) -FI1(N1)*FIZ2(N2) )*half*ovlapphi             !! J_z,phi_ij
             !
             ! tensor components with only K.ne.0 case
             vjrr=-iunit*( FIR1(N1)*FI2(N2) -FI1(N1)*FIR2(N2) )*half*ovlapr               !! J_r,r_ij
             vjrz=-iunit*( FIR1(N1)*FI2(N2) -FI1(N1)*FIR2(N2) )*half*ovlapz               !! J_r,r_ij
             vjff=-iunit*( iunit*(LAM1(N1)-LAM2(N2))*y*FI1(N1)*FI2(N2))*half*ovlapphi     !! J_phi,phi_ij
             vjzr=-iunit*( FIZ1(N1)*FI2(N2) -FI1(N1)*FIZ2(N2) )*half*ovlapr               !! J_z,r_ij
             vjzz=-iunit*( FIZ1(N1)*FI2(N2) -FI1(N1)*FIZ2(N2) )*half*ovlapz               !! J_z,z_ij
             !
             ! these fields have (e^+iKphi - e^-iKphi) structure, hence the sign
             if(LAM1(n1) + LAM2(n2) +1 .eq. abs(K_qrpa_responce)) then
               vjrr = -vjrr; vjff = -vjff; vjzr = -vjzr;
             end if
             !
             ! these fields have (e^+iKphi - e^-iKphi) structure, hence the sign
             if(LAM1(n1) - LAM2(n2) .gt. 0) then
               vjrz = -vjrz; vjzz = -vjzz
             end if
             !
             ! zero those that don't contribute in K=0 case
             if(K_qrpa_responce.eq.0) then
               vjrr = zero; vjff = zero; vjzr = zero; vjrz = zero; vjzz = zero;
             end if
             !
             ! matrix elements from conjugate states (simplex -i states)
             vjrfb = -vjrf ; vjfrb = -vjfr ; vjfzb =  vjfz ; vjzfb = -vjzf;
             vjrrb = -vjrr ; vjrzb =  vjrz ; vjffb = -vjff ; vjzrb = -vjzr;
             vjzzb =  vjzz ;
             !
             ! correct matrix element terms when using normal HO basis
             If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) then
               vso = (+( IUnit*LAM2(N2)*y*FIR1(N1)*FI2(N2) + IUnit*LAM1(N1)*y*FI1(N1)*FIR2(N2) ) *ovlapz &  !! nabla.J_ij
                      -( IUnit*LAM1(N1)*y*FI1(N1)*FIZ2(N2) + IUnit*LAM2(n2)*y*FIZ1(N1)*FI2(N2) ) *ovlapr &
                      -( FIZ1(N1)*FIR2(N2) - FIR1(N1)*FIZ2(N2) )*ovlapphi )/IUnit
               vsob=vso
               !! add here correction for tensor terms
             end if   
             ! * * time odd part * *
             !
             ! variation w.r.t. j
             vjr = ( FIR1(N1)*FI2(N2) - FI1(N1)*FIR2(N2) )*ovlap  !! imag. unit added later
             vjz = ( FIZ1(N1)*FI2(N2) - FI1(N1)*FIZ2(N2) )*ovlap  !! imag. unit added later
             vjfi = half*y*FI1(N1)*FI2(N2)*(LAM1(N1)+LAM2(N2))*ovlap
             !
             ! variation w.r.t. s
             vspinfi  = FI1(N1)*FI2(N2)*ovlapphi
             vspinr   = FI1(N1)*FI2(N2)*ovlapr
             vspinz   = FI1(N1)*FI2(N2)*ovlapz
             !
             ! variation w.r.t. T
             vttbase = FIR1(N1)*FIR2(N2) +FIZ1(N1)*FIZ2(N2)
             vttfi = ( vttbase -LAM1(N1)*LAM2(N2)*y*y*FI1(N1)*FI2(N2) )*ovlapphi
             vttr  = ( vttbase -LAM1(N1)*LAM2(N2)*y*y*FI1(N1)*FI2(N2) )*ovlapr
             vttz  = ( vttbase +LAM1(N1)*LAM2(N2)*y*y*FI1(N1)*FI2(N2) )*ovlapz
             !
             ! these fields have (e^+iKphi - e^-iKphi) structure, hence the sign
             if(LAM1(n1) + LAM2(n2) +1 .eq. abs(K_qrpa_responce)) then
               vspinr = -vspinr
               vttr  = -vttr
             end if
             !
             ! these fields have (e^+iKphi - e^-iKphi) structure, hence the sign
             if(LAM1(n1) - LAM2(n2) .gt. 0) then
               vspinz = -vspinz
               vjfi = -vjfi
               vttz = -vttz
             end if
             !
             ! matrix elements from conjugate states (simplex -i states)
             vspinrb  = -vspinr  ; vspinfib = -vspinfi ;  vspinzb  = vspinz
             vttrb = -vttr       ; vttfib = -vttfi     ;  vttzb = vttz
             !
             ! zero those that don't contribute in K=0 case
             if(K_qrpa_responce.eq.0) then
               vspinr = zero; vspinz = zero; vspinrb = zero; vspinzb = zero; vjfi = zero;
               vttr   = zero; vttz   = zero; vttrb   = zero; vttzb   = zero;
             end if
             !
             ! correct matrix element terms when using normal HO basis
             If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) then
               vspinfib = vspinfi
             end if
             !
             ! matrix elements with respect of +i states
             hN_c(j) = hN_c(j) +vh*vn_c(ihil)     &                    !! rho_n
                               +vhb*vhbn_c(ihil)  &                    !! tau_n
                               +vjr*vjRN_c(ihil)*iunit2 &              !! jr_n
                               +vjz*vjZN_c(ihil)*iunit2 &              !! jz_n
                               +vjfi*vjFIN_c(ihil) &                   !! jfi_n
                               +vdel*vdn_c(ihil) &                     !! Delta_n
                               +vso*vsn_c(ihil)  &                     !! nabla J
                               +vjrf*vJrfin_c(ihil) &                  !! J_r,fi_n
                               +vjfr*vJfirn_c(ihil) &                  !! J_fi,r_n
                               +vjfz*vJfizn_c(ihil) &                  !! J_fi,z_n
                               +vjzf*vJzfin_c(ihil) &                  !! J_z,fi_n
                               +vjrr*vJrrn_c(ihil) &                   !! J_r,r_n
                               +vjrz*vJrzn_c(ihil) &                   !! J_r,z_n
                               +vjff*vJfifin_c(ihil) &                 !! J_fi,fi_n
                               +vjzr*vJzrn_c(ihil) &                   !! J_z,r_n
                               +vjzz*vJzzn_c(ihil) &                   !! J_z,z_n
                               +vspinr*vSIGRn_c(ihil) &                !! spin_r
                               +vspinfi*vSIGFIn_c(ihil) &              !! spin_fi 
                               +vspinz*vSIGZn_c(ihil) &                !! spin_z
                               +vttr*vTRn_c(ihil) &                    !! T_r
                               +vttfi*vTFIn_c(ihil) &                  !! T_fi
                               +vttz*vTZn_c(ihil)                      !! T_z
             !
             hP_c(j) = hP_c(j) +vh*vp_c(ihil)     &                    !! rho_p
                               +vhb*vhbp_c(ihil)  &                    !! tau_n
                               +vjr*vjRP_c(ihil)*iunit2 &              !! jr_p
                               +vjz*vjZP_c(ihil)*iunit2 &              !! jz_p
                               +vjfi*vjFIP_c(ihil) &                   !! jfi_p
                               +vdel*vdp_c(ihil) &                     !! Delta_p
                               +vso*vsp_c(ihil) &                      !! nabla J
                               +vjrf*vJrfip_c(ihil) &                  !! J_r,fi_p
                               +vjfr*vJfirp_c(ihil) &                  !! J_fi,r_p
                               +vjfz*vJfizp_c(ihil) &                  !! J_fi,z_p
                               +vjzf*vJzfip_c(ihil) &                  !! J_z,fi_p
                               +vjrr*vJrrp_c(ihil) &                   !! J_r,r_p
                               +vjrz*vJrzp_c(ihil) &                   !! J_r,z_p
                               +vjff*vJfifip_c(ihil) &                 !! J_fi,fi_p
                               +vjzr*vJzrp_c(ihil) &                   !! J_z,r_p
                               +vjzz*vJzzp_c(ihil) &                   !! J_z,z_p
                               +vspinr*vSIGRp_c(ihil) &                !! spin_r
                               +vspinfi*vSIGFIp_c(ihil) &              !! spin_fi 
                               +vspinz*vSIGZp_c(ihil) &                !! spin_z
                               +vttr*vTRp_c(ihil) &                    !! T_r
                               +vttfi*vTFIp_c(ihil) &                  !! T_fi
                               +vttz*vTZp_c(ihil)                      !! T_z
             !
             ! matrix elements with respect of conjugate (-i) simplex-y states (note that spin-independent part is the same)
             BhN_c(j)=BhN_c(j) +vh*vn_c(ihil)     &                    !! rho_n
                               +vhb*vhbn_c(ihil)  &                    !! tau_n
                               +vjr*vjRN_c(ihil)*iunit2 &              !! jr_n
                               +vjz*vjZN_c(ihil)*iunit2 &              !! jz_n
                               +vjfi*vjFIN_c(ihil) &                   !! jfi_n
                               +vdel*vdn_c(ihil) &                     !! Delta_n
                               +vsob*vsn_c(ihil)  &                    !! nabla J
                               +vjrfb*vJrfin_c(ihil) &                 !! J_r,fi_n
                               +vjfrb*vJfirn_c(ihil) &                 !! J_fi,r_n
                               +vjfzb*vJfizn_c(ihil) &                 !! J_fi,z_n
                               +vjzfb*vJzfin_c(ihil) &                 !! J_z,fi_n   
                               +vjrrb*vJrrn_c(ihil) &                  !! J_r,r_n
                               +vjrzb*vJrzn_c(ihil) &                  !! J_r,z_n
                               +vjffb*vJfifin_c(ihil) &                !! J_fi,fi_n
                               +vjzrb*vJzrn_c(ihil) &                  !! J_z,r_n
                               +vjzzb*vJzzn_c(ihil) &                  !! J_z,z_n
                               +vspinrb*vSIGRn_c(ihil) &               !! spin_r
                               +vspinfib*vSIGFIn_c(ihil) &             !! spin_fi 
                               +vspinzb*vSIGZn_c(ihil) &               !! spin_z
                               +vttrb*vTRn_c(ihil) &                   !! T_r
                               +vttfib*vTFIn_c(ihil) &                 !! T_fi
                               +vttzb*vTZn_c(ihil)                     !! T_z

             !
             BhP_c(j)=BhP_c(j) +vh*vp_c(ihil)     &                    !! rho_p
                               +vhb*vhbp_c(ihil)  &                    !! tau_n
                               +vjr*vjRP_c(ihil)*iunit2 &              !! jr_p
                               +vjz*vjZP_c(ihil)*iunit2 &              !! jz_p
                               +vjfi*vjFIP_c(ihil) &                   !! jfi_p
                               +vdel*vdp_c(ihil) &                     !! Delta_p
                               +vsob*vsp_c(ihil) &                     !! nabla J
                               +vjrfb*vJrfip_c(ihil) &                 !! J_r,fi_p
                               +vjfrb*vJfirp_c(ihil) &                 !! J_fi,r_p
                               +vjfzb*vJfizp_c(ihil) &                 !! J_fi,z_p
                               +vjzfb*vJzfip_c(ihil) &                 !! J_z,fi_p
                               +vjrrb*vJrrp_c(ihil) &                  !! J_r,r_p
                               +vjrzb*vJrzp_c(ihil) &                  !! J_r,z_p
                               +vjffb*vJfifip_c(ihil) &                !! J_fi,fi_p
                               +vjzrb*vJzrp_c(ihil) &                  !! J_z,r_p
                               +vjzzb*vJzzp_c(ihil) &                  !! J_z,z_p
                               +vspinrb*vSIGRp_c(ihil) &               !! spin_r
                               +vspinfib*vSIGFIp_c(ihil) &             !! spin_fi 
                               +vspinzb*vSIGZp_c(ihil) &               !! spin_z
                               +vttrb*vTRp_c(ihil) &                   !! T_r
                               +vttfib*vTFIp_c(ihil) &                 !! T_fi
                               +vttzb*vTZp_c(ihil)                     !! T_z
             !
             dN_c(j)=dN_c(j)+ vh*davn_c(ihil)
             BdN_c(j)=BdN_c(j)+ vh*dbvn_c(ihil)
             dP_c(j)=dP_c(j)+ vh*davp_c(ihil)
             BdP_c(j)=BdP_c(j)+ vh*dbvp_c(ihil)
             !
             j=j+1  !! go one step forward in the matrix array
             !
          End Do !! Do N2=1,iQRPA_IDI(ib) !! states in the block
        End Do !! Do N1=1,iQRPA_IDF(ib)   !! states in the block
      End Do !! Do ihil=1,nghl
    End Do ! ibq=1,iQRPA_NB
    !
    if(USE_ONLY_F20) then
      hN_c = zero; hP_c = zero; BhN_c = zero; BhP_c = zero; 
      dN_c = zero; dP_c = zero; BdN_c = zero; BdP_c = zero; 
    end if
    !
!WRITE(*,*) "gamdel:"
!call qrpa_printProp("hN_c",hN_c)
!call qrpa_printProp("BhN_c",BhN_c)

   !
    Call get_CPU_time('qrpa_gamdel',1)
    !
  End Subroutine qrpa_gamdel
  !=================================================================================================================================
  !  
  !=================================================================================================================================
  Subroutine qrpa_f_gamdel
    !---------------------------------------------------------------------
    ! Fh_ab  in configurational space
    ! Fh     in qp space, both referenced in full
    !
    ! note: ovlapsy may not work properly, if operator is
    !       a product of sigma_y and phi dependent quantity
    !---------------------------------------------------------------------
    Use HFBTHO, only : nqx
    Use HFBTHO, only : fl, fh, ia, ns, nr, nz, npr
    Use HFBTHO, only : QHLA_opt,FI1R_opt,FI1Z_opt,FI2D_opt,y_opt   
    Implicit None
    Integer(ipr) :: i,j,ihil,im,nsa,JA,N1,N2,ibq,ovlap
    Real(pr)     :: vfn,vfp,zcoord,rcoord,y
    complex(pr)  :: vh,vly1,vly2, vspiny, vspinyb, vspinz, vlz, vspinzb, vlzb
    complex(pr)  :: vpz,vpx1,vpx2,vcoz,vcox
    Real(pr), allocatable :: FI1(:),FI2(:),FIR1(:),FIR2(:),FIZ1(:),FIZ2(:)
    Integer (ipr), allocatable :: LAM1(:),LAM2(:),SS1(:),SS2(:),NZ1(:),NZ2(:),NR1(:),NR2(:)
    Real(pr)     :: vfna(nghl),vfpa(nghl),vfnb(nghl),vfpb(nghl)
    complex(pr)  :: iunit
    complex(pr)  :: ovlaplz,ovlaply    !! overlap when l_z or l_y operator
    complex(pr)  :: ovlapsz,ovlapsy    !! overlap when sigma_z or sigma_y operator
    complex(pr)  :: ovlappz,ovlappx    !! overlap for momentum P_z and P_x operator
    complex(pr)  :: ovlapcoz,ovlapcox  !! overlap for coordinate Q_z operator 
    !
    Parameter(iunit = Cmplx(0.0_pr,1.0_pr,kind=pr))
    !-----------------------------------------------------------------
    ! F_ab in configurational space: FhNab_c,FhPab_c,BFhNab_c,BFhPab_c
    !-----------------------------------------------------------------
    FhNab_c=zero; FhPab_c=zero; BFhNab_c=zero; BFhPab_c=zero;
    ! fields
    select case(QRPA_Operator_type)
    case(1)
      vfna = qrpa_electric_field(0,zero,fl,fh)
      vfpa = qrpa_electric_field(1,zero,fl,fh) 
    case(2)
      vfna = qrpa_spin_field(0,zero,fl,fh)
      vfpa = qrpa_spin_field(1,zero,fl,fh) 
    case(3)
      vfna = qrpa_spin_field(0,zero,fl,fh)
      vfpa = qrpa_spin_field(1,zero,fl,fh) 
      vfnb = qrpa_l_field(0,zero,fl,fh)
      vfpb = qrpa_l_field(1,zero,fl,fh) 
    case(4)
    case(5)
    case(6)
    case(7)
    case default
      STOP "Unknown operator type"
    end select
    !
    if(allocated(FI1)) Deallocate(FI1,FI2,FIR1,FIR2,FIZ1,FIZ2)
    allocate(FI1(nqx),FI2(nqx),FIR1(nqx),FIR2(nqx),FIZ1(nqx),FIZ2(nqx))
    !
    if(allocated(lam1)) Deallocate(lam1,ss1,lam2,ss2,nz1,nz2,nr1,nr2)
    allocate(lam1(nqx),lam2(nqx),ss1(nqx),ss2(nqx),nz1(nqx),nz2(nqx),nr1(nqx),nr2(nqx))
    !
    Do ibq=1,iQRPA_NB
       Do ihil=1,nghl
         !
         ! run trough the left basis vectors and fill the array
         Do N1=1,iQRPA_IDF(ibq)
           IM=ia(iQRPA_BLF(ibq)); JA=IM+N1; NSA=NS(JA);
           FI1(N1)  = QHLA_opt(JA,ihil);
           FIR1(N1) = FI1R_opt(JA,ihil)
           FIZ1(N1) = FI1Z_opt(JA,ihil)
           SS1(N1) = NSA
           NZ1(N1) = NZ(JA)
           NR1(N1) = NR(JA)
           if(NSA.gt.0) then
             LAM1(N1) = (iQRPA_OMF(ibq) -1)/2
           else
             LAM1(N1) = -(iQRPA_OMF(ibq) +1)/2
           end if
           If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) LAM1(N1) = abs(LAM1(N1))  !! normal HO, positive omega states have positive lambda
         End Do
         !
         ! run trough the right basis vectors and fill the array
         Do N1=1,iQRPA_IDI(ibq)
           IM=ia(iQRPA_BLI(ibq)); JA=IM+N1; NSA=NS(JA);
           FI2(N1)  = QHLA_opt(JA,ihil);
           FIR2(N1) = FI1R_opt(JA,ihil)
           FIZ2(N1) = FI1Z_opt(JA,ihil)
           SS2(N1) = NSA
           NZ2(N1) = NZ(JA)
           NR2(N1) = NR(JA)
           if(NSA.gt.0) then
             LAM2(N1) = (iQRPA_OMI(ibq) -1)/2
           else
             LAM2(N1) = -(iQRPA_OMI(ibq) +1)/2
           end if
           If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) LAM2(N1) = abs(LAM2(N1))  !! normal HO, positive omega states have positive lambda
         End Do
         !
         vfn=vfna(ihil); vfp=vfpa(ihil);
         j=iQRPA_DS(ibq)
         zcoord=fh(ihil)
         rcoord=fl(ihil)
         y=y_opt(ihil)   !! y=1/r
         !
         Do N2=1,iQRPA_IDI(ibq)
           Do N1=1,iQRPA_IDF(ibq)
              !
              ovlap = 0; ovlapsz = zero; ovlapsy = zero; ovlaply = zero; ovlaplz = zero; ovlapsz = zero;
              ovlappz = zero; ovlappx = zero ; ovlapcoz = zero; ovlapcox = zero
              if( abs(LAM1(n1) + LAM2(n2)) .eq. abs(K_qrpa_responce) ) ovlap = 1
              if( abs(LAM1(n1) - LAM2(n2)) .eq. abs(K_qrpa_responce) ) ovlap = 1
              !
              if( (LAM1(n1) .eq. -LAM2(n2)) .and. (abs(K_qrpa_responce).eq.1) ) ovlapsy = -1.0_pr; !! this is also used for M1. s_y is K=1 operator
              !
              if(abs(LAM1(n1)-LAM2(n2)).eq.1) ovlaply = 1.0_pr
              if(abs(LAM1(n1)+LAM2(n2)).eq.1) ovlaply = 1.0_pr
              if(abs(K_qrpa_responce).ne.1)   ovlaply = zero;  !! l_y is K=1 operator

              ! full conditions from the angular integrals of exponentials for the LINEAR MOMENTUM and COORDINATE operator
              if( (abs(LAM1(n1) - LAM2(n2)) .eq. abs(K_qrpa_responce)) .and. (abs(K_qrpa_responce).eq.0) )   ovlappz  = 1.0_pr;              ! P_z
              if( (abs(LAM1(n1) - LAM2(n2)) .eq. abs(K_qrpa_responce)) .and. (abs(K_qrpa_responce).eq.1) )   ovlappx  = 1.0_pr;              ! P_x
              if( (abs(LAM1(n1) - LAM2(n2)) .eq. abs(K_qrpa_responce)) .and. (abs(K_qrpa_responce).eq.0) )   ovlapcoz = 1.0_pr;              ! Q_z
              if( (abs(LAM1(n1) - LAM2(n2)) .eq. abs(K_qrpa_responce)) .and. (abs(K_qrpa_responce).eq.1) )   ovlapcox = 1.0_pr;              ! Q_x
              !
              !
              !! hormal HO and s_z seems not to work properly... fix this part in future version
              If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) then
                 if( LAM1(n1).eq.LAM2(n2) ) then ; ovlaplz = 1.0_pr ; ovlapsz = 1.0 ; end if 
                 if( SS1(n1).lt.0 ) ovlapsz = ovlapsz * (-1.0_pr)
              end if

              vh=FI1(N1)*FI2(N2)*ovlap
              vspiny=FI1(N1)*FI2(N2)*half*ovlapsy  ! s = (1/2) sigma
              vspinz=FI1(N1)*FI2(N2)*half*ovlapsz  ! s = (1/2) sigma
              !
              ! terms for l_y
              vly1 =-0.25_pr*iunit*rcoord*(FI1(N1)*FIZ2(N2)-FIZ1(N1)*FI2(N2))*ovlaply &
                    +0.25_pr*iunit*zcoord*(FI1(N1)*FIR2(N2)-FIR1(N1)*FI2(N2))*ovlaply

              vly2 = 0.25_pr*iunit*( (LAM2(N2)+LAM1(N1))*zcoord*y*FI1(N1)*FI2(N2))*ovlaply
              !
              vlz = LAM1(N1)*FI1(N1)*FI2(N2)*ovlaplz
              !
              vpz  = 0.50_pr*iunit*(FI1(N1)*FIZ2(N2)-FIZ1(N1)*FI2(N2))*ovlappz
              vpx1 = 0.25_pr*iunit*(FI1(N1)*FIR2(N2)-FIR1(N1)*FI2(N2))*ovlappx
              vpx2 = 0.25_pr*iunit*FI1(N1)*y*FI2(N2)*(LAM1(N1)+LAM2(N2))*ovlappx
              !
              ! CONDITION FOR THE SIGN
              if(LAM1(N1) - LAM2(N2) .gt. 0) then
                vpx2 = -vpx2                                                     !krisa! vpx2 changes sign depending on the chosen condition for lambda and K numbers
              end if
              !
              !Q_z and Q_x momentum operator (element of Q_y vanishes)
              vcoz = (1/Dble(npr(1)+npr(2)))*FI1(N1)*zcoord*FI2(N2)*ovlapcoz
              vcox = (1/(2*Dble(npr(1)+npr(2))))*FI1(N1)*rcoord*FI2(N2)*ovlapcox
              !
              ! this term has (e^+iKphi - e^-iKphi) structure, hence the sign
              if(LAM1(n1) - LAM2(n2) .gt. 0) then
                 vly2 = -vly2
              end if
              !
              ! simplex -i states
              vspinyb = -vspiny 
              !
              ! negative Omega states
              vspinzb = -vspinz; vlzb = -vlz
              !
              select case(QRPA_Operator_type)
              case(1)
                ! electric operator
                ! simplex +i states and simplex -i state are same for electric operator
                FhNab_c(j) = FhNab_c(j)  + vh *vfn 
                FhPab_c(j) = FhPab_c(j)  + vh *vfp 
                BFhNab_c(j)= BFhNab_c(j) + vh *vfn 
                BFhPab_c(j)= BFhPab_c(j) + vh *vfp 
              case(2)
                ! spin operator
                ! simplex +i states and simplex -i state are different
                ! note: operator s_11 = i*s_y
                 FhNab_c(j)  = FhNab_c(j)  + (iunit*vspiny +vspinz) *vfna(ihil)
                 FhPab_c(j)  = FhPab_c(j)  + (iunit*vspiny +vspinz) *vfpa(ihil)
                 BFhNab_c(j) = BFhNab_c(j) + (iunit*vspinyb+vspinzb)*vfna(ihil)
                 BFhPab_c(j) = BFhPab_c(j) + (iunit*vspinyb+vspinzb)*vfpa(ihil)
              case(3)
                ! M1 operator
                ! simplex +i states and simplex -i state are different for spin part
                ! note: operator l_11 = i*l_y
                FhNab_c(j)  = FhNab_c(j)  + (iunit*vspiny +vspinz) *vfna(ihil) + (iunit*(vly1+vly2) +vlz)*vfnb(ihil)
                FhPab_c(j)  = FhPab_c(j)  + (iunit*vspiny +vspinz) *vfpa(ihil) + (iunit*(vly1+vly2) +vlz)*vfpb(ihil)
                BFhNab_c(j) = BFhNab_c(j) + (iunit*vspinyb+vspinzb)*vfna(ihil) + (iunit*(vly1+vly2) +vlzb)*vfnb(ihil)
                BFhPab_c(j) = BFhPab_c(j) + (iunit*vspinyb+vspinzb)*vfpa(ihil) + (iunit*(vly1+vly2) +vlzb)*vfpb(ihil)
              case(4)
                ! LINEAR MOMENTUM (P10,P11) OPERATOR
                ! Since P10 and P11 are not present simultaneusly, they can be added here
                FhNab_c(j)  = FhNab_c(j)  + vpz - (vpx1+vpx2) 
                FhPab_c(j)  = FhPab_c(j)  + vpz - (vpx1+vpx2)
                BFhNab_c(j) = BFhNab_c(j) + vpz - (vpx1+vpx2)
                BFhPab_c(j) = BFhPab_c(j) + vpz - (vpx1+vpx2)
              case(5)
                ! COORDINATE (Q10,Q11) OPERATOR ------------------------------------------------
                FhNab_c(j)  = FhNab_c(j)  + vcoz - vcox
                FhPab_c(j)  = FhPab_c(j)  + vcoz - vcox
                BFhNab_c(j) = BFhNab_c(j) + vcoz - vcox
                BFhPab_c(j) = BFhPab_c(j) + vcoz - vcox
              case(6)
                ! TOTAL ANGULAR MOMENTUM (J11) OPERATOR ----------------------------------------
                ! simplex +i states and simplex -i state are different for spin part !
                ! note: (tilde)L_11 = 1/Sqrt[2] * (L_(11)+L_(1-1)) = -i*L_y is the K = 1 operator (plus in the bracket is because L (axial vector) is parity=+1 operator)
                FhNab_c(j)  = FhNab_c(j)  + iunit*vspiny  + iunit*(vly1+vly2)
                FhPab_c(j)  = FhPab_c(j)  + iunit*vspiny  + iunit*(vly1+vly2)
                BFhNab_c(j) = BFhNab_c(j) + iunit*vspinyb + iunit*(vly1+vly2)
                BFhPab_c(j) = BFhPab_c(j) + iunit*vspinyb + iunit*(vly1+vly2)
              case(7)
                ! TOTAL ANGULAR MOMENTUM (Jy) OPERATOR ----------------------------------------
                ! simplex +i states and simplex -i state are different for spin part !
                ! note: we use simply J_y = S_y + L_y  (J_y is K = 1 operator)
                FhNab_c(j)  = FhNab_c(j)  + vspiny  + (vly1+vly2)
                FhPab_c(j)  = FhPab_c(j)  + vspiny  + (vly1+vly2)
                BFhNab_c(j) = BFhNab_c(j) + vspinyb + (vly1+vly2)
                BFhPab_c(j) = BFhPab_c(j) + vspinyb + (vly1+vly2)
              case default
                STOP "Unknown operator type"
              end select
              !
              j=j+1
              !
           End Do   ! N1=1,iQRPA_IDF(ibq)
         End Do     ! N2=1,iQRPA_IDI(ibq)
      End Do  !! ihil
    End Do !! ibq
    !
    !----------------------------------------------
    ! F20 and F02 in qp space: FhN_c, FhP_c, BFhN_c, BFhP_c
    ! F11 calculated if required by EFA
    !----------------------------------------------
    Call qrpa_F20(VqpN,UqpN,FhNab_c,BFhNab_c,FhN_c,BFhN_c)
    Call qrpa_F20(VqpP,UqpP,FhPab_c,BFhPab_c,FhP_c,BFhP_c)
    if(ENABLE_N_EFA) Call qrpa_F11(VqpN,UqpN,FhNab_c,BFhNab_c,F11N_c,BF11N_c)
    if(ENABLE_P_EFA) Call qrpa_F11(VqpP,UqpP,FhPab_c,BFhPab_c,F11P_c,BF11P_c)
!
!if(ENABLE_N_EFA) then
!write(*,*) "f11"
!call qrpa_printprop("F11",F11N_c)
!call qrpa_printprop("F11b",BF11N_c)
!end if

if(.false.) then
call qrpa_printprop("FhNab_c",FhNab_c)
!do ibq=99,222; write(*,'(F10.6,SP,F10.6,"i  ")',advance='NO') FhNab_c(ibq) ; end do;
call qrpa_printprop("BFhNab_c",BFhNab_c)
!do ibq=99,222; write(*,'(F10.6,SP,F10.6,"i  ")',advance='NO') BFhNab_c(ibq) ; end do;

qrpa_aux = FhNab_c-BFhNab_c ; call qrpa_printProp("FhNab_c-BFhNab_c",qrpa_aux)
!do ibq=99,222; write(*,'(F10.6,SP,F10.6,"i  ")',advance='NO') qrpa_aux(ibq) ; end do;

qrpa_aux = FhNab_c+BFhNab_c ; call qrpa_printProp("FhNab_c+BFhNab_c",qrpa_aux)
!do ibq=99,222; write(*,'(F10.6,SP,F10.6,"i  ")',advance='NO') qrpa_aux(ibq) ; end do;


!qrpa_aux = F11N_c -BF11N_c ; call qrpa_printProp("F11-BF11",qrpa_aux)
!qrpa_aux = F11N_c +BF11N_c ; call qrpa_printProp("F11+BF11",qrpa_aux)
!stop
end if

    !
    Deallocate(FI1,FI2,FIR1,FIR2,FIZ1,FIZ2,LAM1,LAM2,ss1,ss2,nz1,nz2,nr1,nr2)
    !
  End Subroutine qrpa_f_gamdel
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpaINI_CUT_ALLOCATE
    !-----------------------------------------------------------------------------------------------------
    ! Alocate qrpa arrays and impose cut-off conditions
    !-----------------------------------------------------------------------------------------------------
    Use HFBTHO, only : ia
    Use HFBTHO, only : brout, brin
    Use HFBTHO, only : RVqpN, RVqpP, RUqpN, RUqpP
    Use HFBTHO, only : ro, tau, aka
    Use HFBTHO, only : alphamix, iiter
    Use HFBTHO, only : nuv, nqp, nb
    Use HFBTHO, only : Parity
    Use HFBTHO, only : id, nb, kd, nl,ns, npar
    Use HFBTHO, only : keyblo 
   !Use hfbtho_solver, only : densit, field, gamdel
 
    Implicit None
    Integer(ipr) :: i,j,k,iblk,iline,ip,iq,ib,ibb,nd,ndb,iw,jprot,jneut,n1,n2,jtr,nsa,lam
    Real(pr)     :: EqpCUT,PqpCUT,pn,eqpe,ela,enb
    complex(pr)  :: phase
    logical :: isactive

    !--------------------------------QQst       modify11-01
    integer(ipr)  :: qqj
    integer(ipr)  :: NumberofH11
   !integer, allocatable :: QQ_select(:,:)
   !open(43,file=trim(adjustl(qq_io_path))//'QQ_binary_output_binfo',status='unknown',form='unformatted')         
    open(43,file='QQ_binary_output_binfo',status='unknown',form='unformatted')         
    open(44,file='44_select_matrix')         
    !--------------------------------QQend      modify11-01
 
    !--------------------------------------------------------------------------
    ! test HFBTHO U,V for orthonormality
    !--------------------------------------------------------------------------
    Call qrpa_test_HFB_UV(RVqpN,RUqpN,1)  
    Call qrpa_test_HFB_UV(RVqpP,RUqpP,2)
    !
    ! Check here if EFA activated
    !
    !--------------------------------------------------------------------------
    ! QRPA: Allocate Arrays 
    !--------------------------------------------------------------------------
    ! Arrays depending on gauss points for time-even densities and fields
    If(Allocated(vn_c)) Deallocate(aka_c,bka_c,vhbn_c,vn_c,vrn_c,vzn_c,vdn_c,vsn_c,davn_c,dbvn_c &
         ,vhbp_c,vp_c,vrp_c,vzp_c,vdp_c,vsp_c,davp_c,dbvp_c  &
         ,ro_c,tau_c,dro_c,dj_c,NABLAR_c,NABLAZ_c,cou_c,vc_c )
    Allocate(aka_c(nghl,2),bka_c(nghl,2),ro_c(nghl,2),tau_c(nghl,2),dro_c(nghl,2),dj_c(nghl,2) & 
         ,vhbn_c(nghl),vn_c(nghl),vrn_c(nghl),vzn_c(nghl),vdn_c(nghl),vsn_c(nghl),davn_c(nghl),dbvn_c(nghl) &
         ,vhbp_c(nghl),vp_c(nghl),vrp_c(nghl),vzp_c(nghl),vdp_c(nghl),vsp_c(nghl),davp_c(nghl),dbvp_c(nghl)  &         
         ,NABLAR_c(nghl,2),NABLAZ_c(nghl,2),cou_c(nghl),vc_c(nghl,nghl) )
    ! Arrays depending on gauss points for time-even densities and fields: tensor terms
    if(Allocated(JRR_c)) Deallocate(JRR_c,JRFI_c,JRZ_c,JFIR_c,JFIFI_c,JFIZ_c,JZR_c,JZFI_c,JZZ_c &
         ,vJRRn_c,vJRFIn_c,vJRZn_c,vJFIRn_c,vJFIFIn_c,vJFIZn_c,vJZRn_c,vJZFIn_c,vJZZn_c &
         ,vJRRp_c,vJRFIp_c,vJRZp_c,vJFIRp_c,vJFIFIp_c,vJFIZp_c,vJZRp_c,vJZFIp_c,vJZZp_c )
    Allocate(JRR_c(nghl,2),JRFI_c(nghl,2),JRZ_c(nghl,2),JFIR_c(nghl,2),JFIFI_c(nghl,2),JFIZ_c(nghl,2) &
         ,JZR_c(nghl,2),JZFI_c(nghl,2),JZZ_c(nghl,2) &
         ,vJRRn_c(nghl),vJRFIn_c(nghl),vJRZn_c(nghl),vJFIRn_c(nghl),vJFIFIn_c(nghl) &
         ,vJFIZn_c(nghl),vJZRn_c(nghl),vJZFIn_c(nghl),vJZZn_c(nghl) &
         ,vJRRp_c(nghl),vJRFIp_c(nghl),vJRZp_c(nghl),vJFIRp_c(nghl),vJFIFIp_c(nghl) &
         ,vJFIZp_c(nghl),vJZRp_c(nghl),vJZFIp_c(nghl),vJZZp_c(nghl)  )
    !
    ! Arrays depending on gauss points for time-odd densities and fields
    If(Allocated(sor_c)) Deallocate(sor_c,soz_c,sofi_c,dsor_c,dsoz_c,dsofi_c, &
         rjr_c,rjz_c,rjfi_c,tsor_c,tsoz_c,tsofi_c,curjr_c,curjz_c,curjfi_c,cursr_c,cursz_c,cursfi_c, &
         vSIGRn_c,vSIGZn_c,vSIGFIn_c,vSIGRp_c,vSIGZp_c,vSIGFIp_c, &
         vSIDRn_c,vSIDZn_c,vSIDFIn_c,vSIDRp_c,vSIDZp_c,vSIDFIp_c, &
         vTRn_c,vTZn_c,vTFIn_c,vTRp_c,vTZp_c,vTFIp_c)

    Allocate(sor_c(nghl,2),soz_c(nghl,2),sofi_c(nghl,2),dsor_c(nghl,2), &
         dsoz_c(nghl,2),dsofi_c(nghl,2),rjr_c(nghl,2),rjz_c(nghl,2), &
         rjfi_c(nghl,2),tsor_c(nghl,2),tsoz_c(nghl,2),tsofi_c(nghl,2), &
         curjr_c(nghl,2),curjz_c(nghl,2),curjfi_c(nghl,2),cursr_c(nghl,2),cursz_c(nghl,2),cursfi_c(nghl,2), &
         vSIGRn_c(nghl),vSIGZn_c(nghl),vSIGFIn_c(nghl),vSIGRp_c(nghl),vSIGZp_c(nghl),vSIGFIp_c(nghl), &
         vSIDRn_c(nghl),vSIDZn_c(nghl),vSIDFIn_c(nghl),vSIDRp_c(nghl),vSIDZp_c(nghl),vSIDFIp_c(nghl), &
         vjRn_c(nghl),vjZn_c(nghl),vjFIn_c(nghl),vjRp_c(nghl),vjZp_c(nghl),vjFIp_c(nghl),             &
         vTRn_c(nghl),vTZn_c(nghl),vTFIn_c(nghl), &
         vTRp_c(nghl),vTZp_c(nghl),vTFIp_c(nghl) )

    !
    ! for density dependent case, use hfb densities
    If(Allocated(ro_c_hfb)) Deallocate(ro_c_hfb,tau_c_hfb,kap_c_hfb)
    Allocate(ro_c_hfb(nghl,2),tau_c_hfb(nghl,2),kap_c_hfb(nghl,2))
    ro_c_hfb = ro ; tau_c_hfb = tau ; kap_c_hfb = aka ;
    ! 
    ! use less cryptic variable names
    iHFB_LEN = nuv
    iHFB_NQP = nqp
    iHFB_NB = NB
    !
    ! Arrays depending on HFB block structure
    If(Allocated(iHFB_ID)) Deallocate(iHFB_ID,iHFB_DS,iHFB_LP,iHFB_OM,iHFB_QPBS,iHFB_PR)
    Allocate(iHFB_ID(iHFB_NB),iHFB_DS(iHFB_NB),iHFB_LP(iHFB_NB),iHFB_OM(iHFB_NB), &
             iHFB_QPBS(iHFB_NB),iHFB_PR(iHFB_NB))
    If(Allocated(VqpN)) Deallocate(VqpN,VqpP,UqpN,UqpP)
    Allocate(VqpN(iHFB_LEN),VqpP(iHFB_LEN),UqpN(iHFB_LEN),UqpP(iHFB_LEN)) 
    If(Allocated(hfb_UA)) Deallocate(hfb_UA,hfb_VA)
    Allocate(hfb_UA(iHFB_LEN),hfb_VA(iHFB_LEN))
    If(Allocated(ipwi_hfb)) Deallocate(ipwi_hfb)
    Allocate(ipwi_hfb(iHFB_NB,2))
    !
    ! calculate block sizes and structure for HFB
    i=1 ; j=1
    Do ib=1,NB
       ND=ID(ib)
       iHFB_DS(ib)=i
       iHFB_ID(ib)=ND
       iHFB_QPBS(ib)=j
       ipwi_hfb(ib,1)=kd(ib,1)
       ipwi_hfb(ib,2)=kd(ib,2)
       If(Parity) Then
         iHFB_LP(ib)=(ib+1)/2
         iHFB_OM(ib)=2*nl(j)+ns(j)
         iHFB_PR(ib) = 3-2*npar(j)
       else
         iHFB_LP(ib)=ib
         iHFB_OM(ib)=2*ib-1
         iHFB_PR(ib) = 0  ! parity not conserved         
       end if
       i = i+(ND*ND)
       j = j+ND
    End Do
    !
    ! parity of the operator 
    select case(QRPA_Operator_type)
    case(1) !! electric
      if( mod(abs(L_qrpa_responce),2).eq.0 ) then
        P_qrpa_response = +1
      else
        P_qrpa_response = -1
      end if
    case(2) !! spin
      if( abs(L_qrpa_responce).ne.1 ) STOP "Wrong L for spin operator"
      P_qrpa_response = +1
    case(3) !! magnetic
      if( abs(L_qrpa_responce).ne.1 ) STOP "Only L=1 case supported for magnetic operator so far..."
      P_qrpa_response = +1
    case(4) !! linear momentum
      if( abs(L_qrpa_responce).ne.1 ) STOP "Wrong L for linear momentum opr."
      P_qrpa_response = -1
    case(5) !! coordinate
      if( abs(L_qrpa_responce).ne.1 ) STOP "Wrong L for coordinate opr."
      P_qrpa_response = -1
    case(6) !! J_11
      if( abs(L_qrpa_responce).ne.1 ) STOP "Wrong L for J_11 opr."
      P_qrpa_response = +1
    case(7) !! J_11
      if( abs(L_qrpa_responce).ne.1 ) STOP "Wrong L for J_y opr."
      P_qrpa_response = +1
    case default
    end select

    If(Parity.and.(P_qrpa_response.eq.-1)) STOP "Parity conserved mode incompatible with requested operator!"  !! check if this condition can be relaxed later
    !
    ! calculate block sizes and structure for QRPA
    ! first the number of QRPA blocks
    !-------------------------------------QQst  modify11-02
    NumberofH11 = 0
    !
    if(allocated(QQ_select)) deallocate(QQ_select)
    Allocate(QQ_select(iHFB_NB,iHFB_NB))
    QQ_select = 0
    !-------------------------------------QQend modify11-02
    !
    iQRPA_NB = 0
    Do ib=1,iHFB_NB
      Do ibb=1,iHFB_NB
        isactive = .false.
        ! check Omega quantum number
        If( iHFB_OM(ib)+iHFB_OM(ibb)+2*K_qrpa_responce .eq. 0 ) isactive = .true.
        If( iHFB_OM(ib)-iHFB_OM(ibb)+2*K_qrpa_responce .eq. 0 ) isactive = .true.
        If( iHFB_OM(ib)+iHFB_OM(ibb)-2*K_qrpa_responce .eq. 0 ) isactive = .true.
        If( iHFB_OM(ib)-iHFB_OM(ibb)-2*K_qrpa_responce .eq. 0 ) isactive = .true.
        ! check parity
        if(parity .and. ((iHFB_PR(ib)*iHFB_PR(ibb)).ne.P_qrpa_response)) isactive = .false.
        ! add block
        if(isactive) iQRPA_NB = iQRPA_NB+1
        !-----------------------------------------------QQst            modify11-03
        ! assign selection rule matrix, nozero: change to 1, zero : keep to 0
        !----------------------
        if(isactive)  then 
           QQ_select(ib, ibb) = 1
           NumberofH11 = NumberofH11 + iHFB_id(ib)*iHFB_id(ibb)         !!!Q: iHFB_id, defined by Markus
        end if
        !-----------------------------------------------QQend           modify11-03
        !
      End Do
    End Do
    !
    !------------------------------------------QQst     modify11-04
    ! print out assigned selection rule matrix to check
    !---------------------------
    do ib=1, iHFB_NB
       do ibb=1, iHFB_NB
          if(QQ_select(ib, ibb)==1) write(44, *) ib, ibb, QQ_select(ib, ibb)
       end do
    end do
    !--------------------------
    ! write binfo to binary file
    !--------------------------
    write(43) real(qrpa_omega)
    write(43) iHFB_NB 
    write(43) NumberofH11 
    write(43) id              !!!Q: id is an array, id(ib) is the number of basis in block ib. Temporary also as nd
    write(43) QQ_select       !!!Q: id is variable defined by Markus.
    !
   !deallocate(QQ_select)
    !------------------------------------------QQend    modify11-04


    If(iQRPA_NB .eq. 0) STOP "iQRPA_NB .eq. 0"
    If(Allocated(iQRPA_IDI)) Deallocate(iQRPA_IDI,iQRPA_IDF,iQRPA_LPI,iQRPA_LPF,iQRPA_DS, &
                                       iQRPA_OMI,iQRPA_OMF,iQRPA_BLI,iQRPA_BLF,iQRPA_BHI, &
                                       iQRPA_BHF,iQRPA_TRANSPOSE)
    Allocate(iQRPA_IDI(iQRPA_NB),iQRPA_IDF(iQRPA_NB),iQRPA_LPI(iQRPA_NB),iQRPA_LPF(iQRPA_NB), &
             iQRPA_DS(iQRPA_NB),iQRPA_OMI(iQRPA_NB),iQRPA_OMF(iQRPA_NB),iQRPA_BLI(iQRPA_NB),  &
             iQRPA_BLF(iQRPA_NB),iQRPA_BHI(iQRPA_NB),iQRPA_BHF(iQRPA_NB),iQRPA_TRANSPOSE(iQRPA_NB))
    !
    ! QRPA block sizes etc. 
    i=1 ; j=1; iQRPA_LEN=0;
    Do ib=1,iHFB_NB
      Do ibb=1,iHFB_NB
        isactive = .false.
        ! check Omega quantum number
        If( iHFB_OM(ib)+iHFB_OM(ibb)+2*K_qrpa_responce .eq. 0 ) isactive = .true.
        If( iHFB_OM(ib)-iHFB_OM(ibb)+2*K_qrpa_responce .eq. 0 ) isactive = .true.
        If( iHFB_OM(ib)+iHFB_OM(ibb)-2*K_qrpa_responce .eq. 0 ) isactive = .true.
        If( iHFB_OM(ib)-iHFB_OM(ibb)-2*K_qrpa_responce .eq. 0 ) isactive = .true.
        ! check parity
        if(parity .and. ((iHFB_PR(ib)*iHFB_PR(ibb)).ne.P_qrpa_response)) isactive = .false.
        If(isactive) then
          iQRPA_IDI(i) = iHFB_ID(ib)      ;  iQRPA_IDF(i) = iHFB_ID(ibb)
          iQRPA_OMI(i) = iHFB_OM(ib)      ;  iQRPA_OMF(i) = iHFB_OM(ibb)
          iQRPA_LPI(i) = iHFB_LP(ib)      ;  iQRPA_LPF(i) = iHFB_LP(ibb)
          iQRPA_BLI(i) = ib               ;  iQRPA_BLF(i) = ibb
          iQRPA_BHI(i) = ib               ;  iQRPA_BHF(i) = ibb
          iQRPA_DS(i) = j
          iQRPA_LEN = iQRPA_LEN +iHFB_ID(ib)*iHFB_ID(ibb)
          j=j+iHFB_ID(ib)*iHFB_ID(ibb)
          i=i+1
        end if
      End Do
    End Do
    !
    ! figure out blocks corresponding transpose
    if((K_qrpa_responce .eq. 0).and.(P_qrpa_response.eq.+1)) then
      Do ib = 1,iQRPA_NB
        iQRPA_TRANSPOSE(ib)=ib
      end do
    else
      Do ib = 1,iQRPA_NB
        Do ibb = 1,iQRPA_NB
          if( (iQRPA_BHI(ib) .eq. iQRPA_BHF(ibb)) .and. (iQRPA_BHF(ib) .eq. iQRPA_BHI(ibb)) ) iQRPA_TRANSPOSE(ib)=ibb
        end do
      end do
    end if
    !
    ! print block structure
    if(.false.) then 
      Write(*,*) "iHFB_NB,iQRPA_NB,nuv,iHFB_LEN,iQRPA_LEN",iHFB_NB,iQRPA_NB,nuv,iHFB_LEN,iQRPA_LEN
      WRITE(*,*) "HFB blocks"
      WRITE(*,*) "iHFB_ID(ib),iHFB_DS(ib),iHFB_OM(ib),iHFB_LP(ib),iHFB_QPBS(ib),Parity"
      Do ib=1,iHFB_NB
        WRITE(*,*) iHFB_ID(ib),iHFB_DS(ib),iHFB_OM(ib),iHFB_LP(ib),iHFB_QPBS(ib),iHFB_PR(ib)
      End do
      WRITE(*,*) " "
      WRITE(*,*) "Operator parity:", P_qrpa_response
      WRITE(*,*) "QRPA blocks"
      WRITE(*,*) "ib,iQRPA_IDI(ib),iQRPA_IDF(ib),iQRPA_DS(ib),iQRPA_OMI(ib),iQRPA_OMF(ib),iQRPA_LPI(ib),iQRPA_LPF(ib),iQRPA_BLI(ib),iQRPA_BLF(ib),transpose,iQRPA_BHI(ib),iQRPA_BHF(ib)"
      Do ib=1,iQRPA_NB
        WRITE(*,*) ib,iQRPA_IDI(ib),iQRPA_IDF(ib),iQRPA_DS(ib),iQRPA_OMI(ib),iQRPA_OMF(ib),iQRPA_LPI(ib),iQRPA_LPF(ib),iQRPA_BLI(ib), &
                 iQRPA_BLF(ib),iQRPA_TRANSPOSE(ib),iQRPA_BHI(ib),iQRPA_BHF(ib)
      End do
      WRITE(*,*) "IA",IA
    end if
    !
    ! Arrays depending on QRPA block structure
    If(Allocated(XN_c)) Deallocate(XN_c,YN_c,XP_c,YP_c)
    Allocate(XN_c(iQRPA_LEN),YN_c(iQRPA_LEN),XP_c(iQRPA_LEN),YP_c(iQRPA_LEN))         
    If(Allocated(hN_c)) Deallocate(hN_c,hP_c,dN_c,dP_c)
    Allocate(hN_c(iQRPA_LEN),hP_c(iQRPA_LEN),dN_c(iQRPA_LEN),dP_c(iQRPA_LEN)) 
    If(Allocated(BhN_c)) Deallocate(BhN_c,BhP_c,BdN_c,BdP_c)
    Allocate(BhN_c(iQRPA_LEN),BhP_c(iQRPA_LEN),BdN_c(iQRPA_LEN),BdP_c(iQRPA_LEN)) 
    If(Allocated(FhN_c)) Deallocate(FhN_c,FhP_c,BFhN_c,BFhP_c,FhNab_c,FhPab_c,BFhNab_c,BFhPab_c)
    Allocate(FhN_c(iQRPA_LEN),FhP_c(iQRPA_LEN),BFhN_c(iQRPA_LEN),BFhP_c(iQRPA_LEN), &
             FhNab_c(iQRPA_LEN),FhPab_c(iQRPA_LEN),BFhNab_c(iQRPA_LEN),BFhPab_c(iQRPA_LEN)) 
    If(Allocated(GreenNplus)) Deallocate(GreenNplus,GreenPplus,GreenNminus,GreenPminus)
    Allocate(GreenNplus(iQRPA_LEN),GreenPplus(iQRPA_LEN),GreenNminus(iQRPA_LEN),GreenPminus(iQRPA_LEN))  
    If(Allocated(qrpa_AUX)) Deallocate(qrpa_AUX,qrpa_BUX,qrpa_AUX20,qrpa_BUX20)
    Allocate(qrpa_AUX(iQRPA_LEN),qrpa_BUX(iQRPA_LEN),qrpa_AUX20(iQRPA_LEN),qrpa_BUX20(iQRPA_LEN))
    If(Allocated(rhoetadm)) Deallocate(rhoetadm,kapetadm,kapetabdm)
    Allocate(rhoetadm(iQRPA_LEN),kapetadm(iQRPA_LEN),kapetabdm(iQRPA_LEN))
    If(Allocated(JyNab_c)) Deallocate(JyNab_c,JyPab_c,BJyNab_c,BJyPab_c,JyN_c,JyP_c,BJyN_c,BJyP_c)
    Allocate(JyNab_c(iQRPA_LEN),JyPab_c(iQRPA_LEN),BJyNab_c(iQRPA_LEN),BJyPab_c(iQRPA_LEN), &
             JyN_c(iQRPA_LEN),JyP_c(iQRPA_LEN),BJyN_c(iQRPA_LEN),BJyP_c(iQRPA_LEN))
    !
    ! Check blocking (later also finite temperature)
    ENABLE_N_EFA = .false. ; ENABLE_P_EFA = .false. ; !! default value
    If(keyblo(1).Ne.0) ENABLE_N_EFA = .true.
    If(keyblo(2).Ne.0) ENABLE_P_EFA = .true.
    !
    ! EFA arrays
    if(Allocated(F11N_c)) Deallocate(F11N_c,BF11N_c,h11N_c,Bh11N_c,P1N_c,P2N_c,GreenNefa, &
                                     efa_fN_factor1,efa_fN_factor2)
    if(ENABLE_N_EFA) Allocate(F11N_c(iQRPA_LEN),BF11N_c(iQRPA_LEN),h11N_c(iQRPA_LEN),Bh11N_c(iQRPA_LEN), & 
                              P1N_c(iQRPA_LEN),P2N_c(iQRPA_LEN),GreenNefa(iQRPA_LEN),efa_fN_factor1(iQRPA_LEN), &
                              efa_fN_factor2(iQRPA_LEN))
    if(Allocated(F11P_c)) Deallocate(F11P_c,BF11P_c,h11P_c,Bh11P_c,P1P_c,P2P_c,GreenPefa, &
                                     efa_fP_factor1,efa_fP_factor2)
    if(ENABLE_P_EFA) Allocate(F11P_c(iQRPA_LEN),BF11P_c(iQRPA_LEN),h11P_c(iQRPA_LEN),Bh11P_c(iQRPA_LEN), & 
                              P1P_c(iQRPA_LEN),P2P_c(iQRPA_LEN),GreenPefa(iQRPA_LEN),efa_fP_factor1(iQRPA_LEN), &
                              efa_fP_factor2(iQRPA_LEN))
    ! set efa arrays to zero, except energy denominators
    if(ENABLE_N_EFA) then; F11N_c = zero; BF11N_c = zero; h11N_c = zero; Bh11N_c = zero; P1N_c = zero; 
                            P2N_c = zero; GreenNefa = one; efa_fN_factor1 = one; efa_fN_factor2 = one; end if
    if(ENABLE_P_EFA) then; F11P_c = zero; BF11P_c = zero; h11P_c = zero; Bh11P_c = zero; P1P_c = zero; 
                            P2P_c = zero; GreenPefa = one; efa_fP_factor1 = one; efa_fP_factor2 = one; end if
    !---------------------------------------------------------------------
    ! Charging HFBTHO brin
    !---------------------------------------------------------------------
    alphamix=1.0_pr; iiter=0
    Call densit
    Call field	
    Call gamdel(.false.,.true.)
    !---------------------------------------------------------------------
    ! make HFB matrices complex 
    !---------------------------------------------------------------------
    VqpN=Zero; UqpN=Zero; VqpP=Zero; UqpP=Zero
    i=0;
    Do ib=1,NB
       ND=ID(ib)
       Do N1=1,ND
         Do N2=1,ND
           i = i+1
           VqpN(i)=RVqpN(i)
           VqpP(i)=RVqpP(i)
           UqpN(i)=RUqpN(i)
           UqpP(i)=RUqpP(i)
         End Do
       End Do
    End Do
    !---------------------------------------------------------------------
    ! set U and V to simplex-y phase.
    !---------------------------------------------------------------------
    i=0;
    Do ib=1,NB
       ND=ID(ib)
       Do N1=1,ND
         Do N2=1,ND
           NSA=NS(ia(ib)+N2);
           if(NSA.gt.0) then
             LAM = (iHFB_OM(ib) -1)/2
           else
             LAM = -(iHFB_OM(ib) +1)/2
           end if
           phase = 1.0_pr
           if( mod(abs(Lam),2) .eq. 1 ) then
             phase = cmplx(0.0_pr,1.0_pr)
             if(lam.lt.0) phase = -phase
           end if
           i = i+1
           If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) phase = cmplx(1.0_pr,0.0_pr)  !! don't change phase for normal HO basis
           VqpN(i)=VqpN(i)*phase
           VqpP(i)=VqpP(i)*phase
           UqpN(i)=UqpN(i)*phase
           UqpP(i)=UqpP(i)*phase
         End Do
       End Do
    End Do
    !
    !---------------------------------------------------------------------
    ! use HFB pairing window for calculation of densities
    !---------------------------------------------------------------------
    USE_HFB_PWI = .true.
    !
    icacou_qrpa=0 !! recalculate coulomb
    !
    call qrpa_calculate_timeodd_CC
    !---------------------------------------------------------------------
    ! deallocate some HFBTHO arrrays ...
    !---------------------------------------------------------------------
    Deallocate(brout,brin,RVqpN,RUqpN,RVqpP,RUqpP)
    !---------------------------------------------------------------------
    ! calculate the external field
    !---------------------------------------------------------------------
    Call qrpa_f_gamdel
    !
    close(43)
    close(44)
    !
  End Subroutine qrpaINI_CUT_ALLOCATE
 

  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_N00(V,U,AN00,AN11)
    !---------------------------------------------------------------------
    ! The diagonal of N00 =  V'* V to AN00 
    ! The diagonal of N11 =  U'* U - V'* V  to AN11 
    !   AN00 =>  0 for particles,  1 for holes
    !   AN11 => -1 for particles, +1 for holes
    !
    ! Use this for the HFB matrices imported from HFBTHO
    !---------------------------------------------------------------------
   !Use HFBTHO, only : nuv, nqp
   !Use HFBTHO, only : nb, id 
    Implicit None
    Real(pr), Intent(in)  :: U(NUV),V(NUV)
    Real(pr), Intent(out) :: AN00(nqp),AN11(nqp)
    Real(pr)     :: REAUX(NUV),REBUX(NUV)
    Real(pr)     :: one,zero,onem
    Integer(ipr) :: i,ib,nd,ip,iq
    !
    one=1.0_pr; zero=0.0_pr; onem=-one
    !
    ip = 1; iq=1
    Do ib=1,NB
       ND=ID(ib)
       Call dsyrk('L','t',nd,nd,one, V(ip),nd,zero,REAUX(ip),nd) !  V'* V 
       Call dsyrk('L','t',nd,nd,one ,U(ip),nd,zero,REBUX(ip),nd) !  U'* U 
       Call dsyrk('L','t',nd,nd,onem,V(ip),nd,one ,REBUX(ip),nd) !- V'* V
       Do i=1,nd
          AN00(i+iq-1)=REAUX(ip+(i-1)*nd+i-1)
          AN11(i+iq-1)=REBUX(ip+(i-1)*nd+i-1)
       End Do
       ip = ip + nd*nd
       iq = iq + nd
    End Do
    !
  End Subroutine qrpa_N00
  !=================================================================================================================================
  !  
  !=================================================================================================================================
  Subroutine qrpa_test_HFB_UV(V,U,it)
    !---------------------------------------------------------------------
    ! calculates V' U - U' V and V' V + U' U; UV-real matrices. 
    !---------------------------------------------------------------------
   !Use HFBTHO, only : nuv
   !Use HFBTHO, only : nb, id
   !Use HFBTHO, only : txb 
    Implicit None
    Integer(ipr), Intent(in) :: it
    Real(pr),     Intent(in) :: U(NUV),V(NUV)
    Real(pr)     :: REAUX(NUV)
    Real(pr)     :: onem,sum1,qf1,qf1t,ef1,sum2,qf2,qf2t,ef2 
    Integer(ipr) :: ib,nd,ip,i,j,k 
    !
    onem=-one 
    !
    ip = 1; k=0; qf1t = zero; qf2t = zero
    Do ib=1,NB
       ND=ID(ib)       
       Call dgemm('t','n',nd,nd,nd,onem,U(ip),nd,V(ip),nd,zero,REAUX(ip),nd) !   V' U       
       Call dgemm('t','n',nd,nd,nd,one ,V(ip),nd,U(ip),nd,one ,REAUX(ip),nd) ! - U' V
       sum1 = zero
       Do i=1,ND
          Do j=1,ND
             sum1 = sum1 + REAUX(ip+(j-1)*ND+i-1)**2
          End Do
       End Do
       Call dgemm('t','n',nd,nd,nd,one ,U(ip),nd,U(ip),nd,zero,REAUX(ip),nd) !  V' V 
       Call dgemm('t','n',nd,nd,nd,one ,V(ip),nd,V(ip),nd,one ,REAUX(ip),nd) !+ U' U
       sum2 = zero
       Do i=1,ND
          Do j=1,ND
             sum2 = sum2 + REAUX(ip+(j-1)*ND+i-1)**2
          End Do
       End Do
       qf1 = sum1/Dfloat(Nd*Nd); qf1 = Sqrt(qf1)
       qf2 = sum2/Dfloat(ND);    qf2 = Sqrt(qf2)
       ef1 = dabs(qf1);          ef2 = Abs(qf2-one)
       If ((ef1.Gt.1.0d-11).Or.(ef2.Gt.1.0d-11)) Then
          Write(*,'(2x,a,2(a,d12.3,2x))') txb(ib),'Max VU-UV:',ef1,'VV+UU-1:',ef2  
          k=k+1    
       End If
       ip = ip + nd*nd
    End Do
    If(k.Eq.0.And.it.Eq.1) Write(*,'(2x,a,2(a,d12.3,2x))') ' NEUTRON HFB UV SOLUTION IS OK'
    If(k.Eq.0.And.it.Eq.2) Write(*,'(2x,a,2(a,d12.3,2x))') ' PROTON  HFB UV SOLUTION IS OK'
  End Subroutine qrpa_test_HFB_UV
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_green_function(omega) 
    !---------------------------------------------------------------------
    ! Calculates energy denomitors. With EFA calculates also the f-factors
    !
    ! GreenN_k1k2=1/(Eqp(k1)+Eqp(k2) +/- omega), all complex
    !---------------------------------------------------------------------
   !Use HFBTHO, only : nb, id
   !Use HFBTHO, only : ReqpN, ReqpP, keyblo
   !Use HFBTHO, only : keyblo, blo123d, bloblo
    Implicit None
    Complex(pr), Intent(in) :: omega
    Integer(ipr) :: i,ib,idi,idf,k1,k2
    Integer(ipr), allocatable :: iHFB_EDS(:)
    Real(pr)     :: gf_AeNm,gf_AePm,gf_ReN,gf_ReP,Romega
    Complex(pr)  :: onec,gf_eN,gf_eP,ffn1,ffn2,ffp1,ffp2
    integer(ipr) :: ibln,iblp,k0n1,k0n2,k0p1,k0p2
    Complex(pr)  :: epsn,epsp
    !    
    Romega=Real(omega,kind=pr);  onec=1.0_pr
    gf_AeNm=1000.0_pr; gf_AePm=1000.0_pr
    !
    if(.not.allocated(iHFB_EDS)) allocate(iHFB_EDS(iHFB_NB))
    !
    i=0;
    do ib=1,NB
      IDI=ID(ib)
      iHFB_EDS(ib) = i
      i=i+IDI
    end do
    !
    i=0;
    !
    ! avoid division by zero with efa (in Eqp_i - Eqp_j part)
    epsn = zero; epsp = zero;
    If(ENABLE_N_EFA) epsn = spacing(1.0_pr)
    If(ENABLE_P_EFA) epsp = spacing(1.0_pr)  
    !
    Do ib=1,iQRPA_NB
       Do k1=1,iQRPA_IDI(ib)
         Do k2=1,iQRPA_IDF(ib)
             i=i+1
             ! neutrons 
             gf_eN=REqpN(iHFB_EDS(iQRPA_BHI(ib))+k1)+REqpN(iHFB_EDS(iQRPA_BHF(ib))+k2); GreenNminus(i)=onec/(gf_eN-omega); GreenNplus(i) =onec/(gf_eN+omega)
             ! protons 
             gf_eP=REqpP(iHFB_EDS(iQRPA_BHI(ib))+k1)+REqpP(iHFB_EDS(iQRPA_BHF(ib))+k2); GreenPminus(i)=onec/(gf_eP-omega); GreenPplus(i) =onec/(gf_eP+omega)
             !
             ! blocking
             ibln=bloblo(keyblo(1),1)  !! hfb block of blocked state
             iblp=bloblo(keyblo(2),2) 
             ffn1 = onec; ffn2 = zero; ffp1 = onec; ffp2 = zero;
             If(ENABLE_N_EFA .and. ibln.Eq. iQRPA_BHI(ib) ) then; K0n1=blo123d(1); else; k0n1 = 0; end if
             If(ENABLE_N_EFA .and. ibln.Eq. iQRPA_BHF(ib) ) then; K0n2=blo123d(1); else; k0n2 = 0; end if
             If(ENABLE_P_EFA .and. iblp.Eq. iQRPA_BHI(ib) ) then; K0p1=blo123d(2); else; k0p1 = 0; end if
             If(ENABLE_P_EFA .and. iblp.Eq. iQRPA_BHF(ib) ) then; K0p2=blo123d(2); else; k0p2 = 0; end if
             !
             If(ENABLE_N_EFA .and. k1.eq.k0n1) then; ffn1 = ffn1 - half; ffn2 = ffn2 -half; end if;
             If(ENABLE_N_EFA .and. k2.eq.k0n2) then; ffn1 = ffn1 - half; ffn2 = ffn2 +half; end if;
             If(ENABLE_P_EFA .and. k1.eq.k0p1) then; ffp1 = ffp1 - half; ffp2 = ffp2 -half; end if;
             If(ENABLE_P_EFA .and. k2.eq.k0p2) then; ffp1 = ffp1 - half; ffp2 = ffp2 +half; end if;
             !              
             ! proton poles
             gf_ReP=gf_eP 
             If(gf_AePm.Ge.Abs(gf_ReP-Romega)) Then
                gf_AePm=Abs(gf_ReP-Romega); gfpoles_imP=i
             Endif
             ! neutron poles
             gf_ReN=gf_eN
             If(gf_AeNm.Ge.Abs(gf_ReN-Romega)) Then
                gf_AeNm=Abs(gf_ReN-Romega); gfpoles_imN=i
             Endif
             ! 
             ! efa arrays CHECK THE SIGN OF ENERGY W.R.T F-FACTORS
             gf_eN=REqpN(iHFB_EDS(iQRPA_BHI(ib))+k1)-REqpN(iHFB_EDS(iQRPA_BHF(ib))+k2)+epsn;
             gf_eP=REqpP(iHFB_EDS(iQRPA_BHI(ib))+k1)-REqpP(iHFB_EDS(iQRPA_BHF(ib))+k2)+epsp;
             if(ENABLE_N_EFA) then
               GreenNefa(i) = onec/(gf_eN-omega);
               efa_fN_factor1(i) = ffn1;  
               efa_fN_factor2(i) = ffn2;
             end if
             if(ENABLE_P_EFA) then 
               GreenPefa(i) = onec/(gf_eP-omega);
               efa_fP_factor1(i) = ffp1;  
               efa_fP_factor2(i) = ffp2;
             end if
          Enddo
       Enddo
    Enddo

!WRITE(*,*) "green:"
!write(*,*) "efa1:"
!WRITE(*,'(20(1x,f5.2))') REAL(efa_fN_factor1)
!write(*,*) "efa2:"
!WRITE(*,'(20(1x,f5.2))') REAL(efa_fN_factor2)
!call qrpa_printprop("efa_fN_factor1",efa_fN_factor1)
!call qrpa_printprop("efa_fN_factor2",efa_fN_factor2)
!call qrpa_printprop("GreenNefa",GreenNefa)

    deallocate(iHFB_EDS)
  End Subroutine qrpa_green_function
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_xy 
    !---------------------------------------------------------------------
    ! without EFA:
    ! XN_c=-hN_c*GreenNminus; YN_c=-BhN_c*GreenNplus         
    ! XP_c=-hP_c*GreenPminus; YP_c=-BhP_c*GreenPplus
    !
    ! with EFA: X,Y,P1,P2 amplitudes, see separate documentation
    !---------------------------------------------------------------------
    Implicit None
    Integer(ipr) :: i
    Do i=1,iQRPA_LEN
      !
      ! neutrons
      if(.not.ENABLE_N_EFA) then
        XN_c(i)=-hN_c(i)*GreenNminus(i); YN_c(i)=-BhN_c(i)*GreenNplus(i)    !! normal even-even qrpa
      else
        XN_c(i)  = -hN_c(i)*GreenNminus(i)*efa_fN_factor1(i)
        YN_c(i)  = -BhN_c(i)*GreenNplus(i)*efa_fN_factor1(i)
        P1N_c(i) = -h11N_c(i)*GreenNefa(i)*efa_fN_factor2(i)
        P2N_c(i) = -Bh11N_c(i)*GreenNefa(i)*efa_fN_factor2(i)
      end if
      !
      ! protons
      if(.not.ENABLE_P_EFA) then
        XP_c(i)=-hP_c(i)*GreenPminus(i); YP_c(i)=-BhP_c(i)*GreenPplus(i)    !! normal even-even qrpa
      else
        XP_c(i)  = -hP_c(i)*GreenPminus(i)*efa_fP_factor1(i)
        YP_c(i)  = -BhP_c(i)*GreenPplus(i)*efa_fP_factor1(i)
        P1P_c(i) = -h11P_c(i)*GreenPefa(i)*efa_fP_factor2(i)
        P2P_c(i) = -Bh11P_c(i)*GreenPefa(i)*efa_fP_factor2(i)
      end if
    End Do

!if(ENABLE_N_EFA) then
!WRITE(*,*) "xy:"
!call qrpa_printprop("XN",XN_c)
!call qrpa_printprop("YN",YN_c)
!call qrpa_printprop("P1N",P1N_c)
!call qrpa_printprop("P2N",P2N_c)
!call qrpa_printprop("GreenNefa",GreenNefa)
!call qrpa_printprop("efa_fN_factor2",efa_fN_factor2)
!call qrpa_printprop("efa_fN_factor2*GreenNefa",efa_fN_factor2*GreenNefa)
!call qrpa_printprop("h11n",h11n_c)
!call qrpa_printprop("bh11n",bh11n_c)

!write(*,*) "f11:"
!WRITE(*,'(20(1x,f5.2))') REAL(F11N_c)
!write(*,*) "f11*efa2:"
!WRITE(*,'(20(1x,f5.2))') REAL(efa_fN_factor2*F11N_c)


!WRITE(*,*) "efa_fN_factor2*GreenNefa:"
!WRITE(*,*) REAL(efa_fN_factor2*GreenNefa)
!WRITE(*,*) "Re h11n"
!WRITE(*,*) REAL(h11n_c)
!WRITE(*,*) "Im h11n"
!WRITE(*,*) Imag(h11n_c)
!end if

  End Subroutine  qrpa_xy
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_TriProd(op1,op2,op3,U,X,V,Fout)
    !---------------------------------------------------------------------
    ! Computes: Fout = op1(U) . op2(X) . op3(V)
    ! Use lowercase characters for op's!
    ! U and V blocks are assumed to be square matrices and X to be general block matrix
    ! Op's as in ZGEMM, but only lowercase. New op '*' added for complex conjugation
    !---------------------------------------------------------------------
    Implicit None
    Complex(pr), Intent(in)  :: U(iHFB_LEN),V(iHFB_LEN),X(iQRPA_LEN)
    Complex(pr), Intent(out) :: Fout(iQRPA_LEN)
    Character, Intent(in)    :: op1,op2,op3
    Character                :: op1l,op2l,op3l
    Complex(pr)              :: czero,cone,conem
    Complex(pr), allocatable :: aux1(:),aux2(:),aux3(:),aux(:)
    Integer(ipr) :: ib,nd,ipa,ipq,ipt,ibHFBi,ibHFBf,ibTHOi,ibTHOf
    !
    czero=0.0_pr; cone=1.0_pr; conem=-cone
    !
    if ((op2 .eq. 'T') .or. (op2 .eq. 'C')) then
      STOP "Use lowercase chars"
    end if
    !
    If(allocated(aux1)) Deallocate(aux1,aux2,aux3)
    Allocate(aux1(iHFB_LEN),aux2(iQRPA_LEN),aux3(iHFB_LEN),aux(iQRPA_LEN))

    op1l = op1; op2l = op2; op3l = op3; 
    aux = czero;

    if(op1 .eq. '*') then ; aux1=Conjg(U); op1l = 'n'; else ; aux1=U; end if
    if(op2 .eq. '*') then ; aux2=Conjg(X); op2l = 'n'; else ; aux2=X; end if
    if(op3 .eq. '*') then ; aux3=Conjg(V); op3l = 'n'; else ; aux3=V; end if

    ! first X.V
    Do ib=1,iQRPA_NB
      ibHFBi = iQRPA_BHI(ib) ; ibHFBf = iQRPA_BHF(ib)
      ibTHOi = iQRPA_BLI(ib) ; ibTHOf = iQRPA_BLF(ib)
      ipa = iHFB_DS(ibHFBi) ; ipq = iQRPA_DS(ib) ; ipt = iQRPA_DS(iQRPA_TRANSPOSE(ib))
      !
      if ((op2l .eq. 't') .or. (op2l .eq. 'c')) then
        Call ZGEMM(op2l,op3l,iHFB_ID(ibHFBf),iHFB_ID(ibHFBi),iHFB_ID(ibHFBi),cone,aux2(ipt),iHFB_ID(ibHFBi),aux3(ipa),iHFB_ID(ibHFBi),czero,aux(ipq),iHFB_ID(ibHFBf)) !! XT * op3(V) or XH * op3(V)
      else
        Call ZGEMM(op2l,op3l,iHFB_ID(ibHFBf),iHFB_ID(ibHFBi),iHFB_ID(ibHFBi),cone,aux2(ipq),iHFB_ID(ibHFBf),aux3(ipa),iHFB_ID(ibHFBi),czero,aux(ipq),iHFB_ID(ibHFBf)) !! X * op3(V)
      end if
    End Do
    ! 
    ! then U.(XV)
    Do ib=1,iQRPA_NB
      ibHFBi = iQRPA_BHI(ib) ; ibHFBf = iQRPA_BHF(ib)
      ibTHOi = iQRPA_BLI(ib) ; ibTHOf = iQRPA_BLF(ib)
      ipa = iHFB_DS(ibHFBf) ; ipq = iQRPA_DS(ib) ; ipt = iQRPA_DS(iQRPA_TRANSPOSE(ib))
      Call ZGEMM(op1l,'n',iHFB_ID(ibHFBf),iHFB_ID(ibHFBi),iHFB_ID(ibHFBf),cone,aux1(ipa),iHFB_ID(ibHFBf),aux(ipq),iHFB_ID(ibHFBf),czero,Fout(ipq),iHFB_ID(ibHFBf)) !! op1(U) * (XV) to Fout
    End Do
    Deallocate(aux1,aux2,aux3,aux)
  End Subroutine qrpa_TriProd
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_DiProdDiag(op1,A,B,Fout)
    !---------------------------------------------------------------------
    ! Computes diagonal blocks of Fout = op1(A) . B, where A and B have
    ! the same block structure as QRPA matrices.
    ! Op as in ZGEMM, but only lowercase. New op '*' added for complex conjugation
    ! NOTE: Off-diagonal blocks are NOT calculated!
    !---------------------------------------------------------------------
    Implicit None
    Complex(pr), Intent(in)  :: A(iQRPA_LEN),B(iQRPA_LEN)
    Complex(pr), Intent(out) :: Fout(iHFB_LEN)
    Character, Intent(in)    :: op1
    Complex(pr)              :: czero,cone,conem
    Complex(pr), allocatable :: aux1(:),aux2(:)
    Character                :: op1l
    Integer(ipr) :: ib,ic,ict,ihdata,iadata,ibdata
    !
    czero=0.0_pr; cone=1.0_pr; conem=-cone
    !
    If(allocated(aux1)) Deallocate(aux1,aux2)
    Allocate(aux1(iHFB_LEN),aux2(iQRPA_LEN))
    aux1=czero;
    !
    op1l = op1
    if(op1 .eq. '*') then ; aux2=Conjg(A); op1l = 'n'; else ; aux2=A; end if
    !
    Do ib=1,iHFB_NB
      ihdata = iHFB_DS(ib)
      Do ic=1,iQRPA_NB
        if(iQRPA_BHF(ic).eq.ib) then
          ict = iQRPA_TRANSPOSE(ic)
          iadata = iQRPA_DS(ic)
          ibdata = iQRPA_DS(ict)
          !
          if ((op1l .eq. 't') .or. (op1l .eq. 'c')) then
            Call ZGEMM(op1l,'n',iHFB_ID(ib),iQRPA_IDF(ic),iQRPA_IDI(ic),cone,aux2(ibdata),iQRPA_IDI(ic),B(ibdata),iQRPA_IDI(ic),cone,aux1(ihdata),iHFB_ID(ib)) 
          else
            Call ZGEMM(op1l,'n',iHFB_ID(ib),iQRPA_IDF(ic),iQRPA_IDI(ic),cone,aux2(iadata),iHFB_ID(ib),B(ibdata),iQRPA_IDI(ic),cone,aux1(ihdata),iHFB_ID(ib)) 
          end if
          !
        End If
      End Do
    End Do
    Fout = aux1
    Deallocate(aux1,aux2)
  End Subroutine qrpa_DiProdDiag
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_GiveTransposed(Fin,Fout)
    !
    ! Gives transpose of QRPA-block type of matrix. Output overwritten in Fout.
    ! Not used, but could be usedful for debug.
    !
    Complex(pr), Intent(in)  :: Fin(iQRPA_LEN)
    Complex(pr), Intent(out) :: Fout(iQRPA_LEN)
    Integer(ipr) :: i,j,ib,ibt,iadata,ibdata
    !
    Do ib=1,iQRPA_NB
      ibt = iQRPA_TRANSPOSE(ib)
      iadata = iQRPA_DS(ib)
      ibdata = iQRPA_DS(ibt)
      Do i=1,iQRPA_IDI(ib)
        Do j=1,iQRPA_IDF(ib)
          Fout(ibdata-1 +i +(j-1)*iQRPA_IDI(IB) ) = Fin(iadata-1 +j +(i-1)*iQRPA_IDF(IB) )
        End Do
      End Do
    End Do
  End Subroutine qrpa_GiveTransposed
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Function qrpa_HFBtrace(Fin)
    !
    ! gives trace of a HFB-block type of matrix
    !
    Complex(pr), Intent(in)  :: Fin(iHFB_LEN)
    Complex(pr) :: qrpa_HFBtrace
    Complex(pr) :: trtmp,ciptrace
    integer(ipr) :: ib,nd,ip
    !
    trtmp = 0.0_pr
    Do ib = 1,iHFB_NB
      nd = iHFB_ID(ib)
      ip = iHFB_DS(ib)
      Call qrpa_ctrace(Fin(ip),nd,ciptrace)
      trtmp = trtmp + ciptrace
    End Do
    qrpa_HFBtrace = trtmp
    !
  End Function qrpa_HFBtrace
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_H20(V,U,H,D,BH,BD)
    !
    ! Calculates H20 and H02 matrix elements. Output to H and BH, which are overwritten
    !
    Implicit None
    Complex(pr), Intent(in)     :: U(iHFB_LEN),V(iHFB_LEN),D(iQRPA_LEN),BD(iQRPA_LEN)
    Complex(pr), Intent(inout)  :: H(iQRPA_LEN),BH(iQRPA_LEN)
    Complex(pr), Allocatable    :: temp(:)
    !
    if(allocated(temp)) Deallocate(temp)
    Allocate(temp(iQRPA_LEN))
    !
    ! V and U are complex
    call qrpa_TriProd('c','n','n',U,H,V,qrpa_AUX)
    call qrpa_TriProd('c','t','n',V,BH,U,qrpa_AUX20)
    call qrpa_TriProd('c','c','n',V,BD,V,qrpa_BUX)
    call qrpa_TriProd('c','n','n',U,D,U,qrpa_BUX20)

    temp = -qrpa_AUX -qrpa_AUX20 -qrpa_BUX +qrpa_BUX20

    call qrpa_TriProd('t','t','*',U,H,V,qrpa_AUX)
    call qrpa_TriProd('t','n','*',V,BH,U,qrpa_AUX20)
    call qrpa_TriProd('t','t','*',V,D,V,qrpa_BUX)
    call qrpa_TriProd('t','*','*',U,BD,U,qrpa_BUX20)

    H = temp
    BH = -qrpa_AUX -qrpa_AUX20 -qrpa_BUX +qrpa_BUX20
    Deallocate(temp)    
  End Subroutine qrpa_H20
  !=================================================================================================================================
  !   
  !=================================================================================================================================
  subroutine qrpa_H11(V,U,H,D,BH,BD)
    !
    ! Calculates H11_1 and H11_2 matrix elements.
    !
    Implicit None
    Complex(pr), Intent(in)     :: U(iHFB_LEN),V(iHFB_LEN),D(iQRPA_LEN),BD(iQRPA_LEN),H(iQRPA_LEN),BH(iQRPA_LEN)
    Complex(pr)                 :: H11(iQRPA_LEN),BH11(iQRPA_LEN)
    !
    ! V and U are complex
    call qrpa_TriProd('c','n','n',U,H,U,qrpa_AUX)
    call qrpa_TriProd('c','t','n',V,BH,V,qrpa_AUX20)
    call qrpa_TriProd('c','n','n',U,D,V,qrpa_BUX)
   !call qrpa_TriProd('c','t','n',V,BD,U,qrpa_BUX20)
    call qrpa_TriProd('c','c','n',V,BD,U,qrpa_BUX20)

   !H11 = qrpa_AUX -qrpa_AUX20 +qrpa_BUX -qrpa_BUX20
    H11 = qrpa_AUX -qrpa_AUX20 +qrpa_BUX +qrpa_BUX20

    call qrpa_TriProd('t','n','*',U,BH,U,qrpa_AUX)
    call qrpa_TriProd('t','t','*',V,H,V,qrpa_AUX20)
    call qrpa_TriProd('t','t','*',U,D,V,qrpa_BUX)
   !call qrpa_TriProd('t','n','*',V,BD,U,qrpa_BUX20)
    call qrpa_TriProd('t','*','*',V,BD,U,qrpa_BUX20)

   !BH11 = qrpa_AUX -qrpa_AUX20 +qrpa_BUX -qrpa_BUX20
    BH11 = qrpa_AUX -qrpa_AUX20 +qrpa_BUX +qrpa_BUX20
    !
  end subroutine qrpa_H11
  !=================================================================================================================================
  !   
  !=================================================================================================================================
  Subroutine qrpa_strength(X,Y,P1,P2,F20,F02,F11,F11b,efa,ImStrength,ReStrength)
    !
    ! Computes imag and real parts of qrpa strength function: sum (F20*X +F02*Y)
    ! Note: transition strength can be obtained from imaginary part by multiplying with -1/pi
    !
    Implicit None
    Complex(pr), Intent(in)  :: X(iQRPA_LEN),Y(iQRPA_LEN),F20(iQRPA_LEN),F02(iQRPA_LEN)
    Complex(pr), Intent(in)  :: P1(iQRPA_LEN),P2(iQRPA_LEN),F11(iQRPA_LEN),F11b(iQRPA_LEN)
    logical, Intent(in)      :: efa
    Real(pr), Intent(out)    :: ImStrength,ReStrength
    Complex(pr)              :: ctrace
    Complex(pr), Allocatable :: temp1(:)
    !
    if(Allocated(temp1)) Deallocate(temp1)
    Allocate(temp1(iQRPA_LEN))
    !
    temp1 = 0.0_pr ; ImStrength = 0.0_pr; ReStrength = 0.0_pr
    temp1 = conjg(F20)*X +conjg(F02)*Y
!write(*,*) "F20 part terms: ",sum(conjg(F20)*X),sum(conjg(F02)*Y)

    if(efa) temp1 = temp1 + conjg(F11)*P1 + conjg(F11b)*P2
    !
!if(efa) write(*,*) "F11 part terms: ",sum(conjg(F11)*P1),sum(conjg(F11b)*P2)
!if(efa) call qrpa_printProp("P1",P1)
!if(efa) call qrpa_printProp("P2",P2)
!if(efa) then ; qrpa_aux = P1-P2 ; call qrpa_printProp("P1-P2",qrpa_aux) ; end if
!if(efa) then ; qrpa_aux = P1+P2 ; call qrpa_printProp("P1+P2",qrpa_aux) ; end if
!if(efa) then ; qrpa_aux = F11-F11b ; call qrpa_printProp("F11-F11b",qrpa_aux) ; end if
!if(efa) then ; qrpa_aux = F11+F11b ; call qrpa_printProp("F11+F11b",qrpa_aux) ; end if



    ImStrength = Imag( sum(temp1) )
    ReStrength = Real( sum(temp1) ,kind=pr)
    Deallocate(temp1)
  End Subroutine qrpa_strength
  !=================================================================================================================================
  ! 
  !=================================================================================================================================
  Subroutine qrpa_F20(V,U,F1,F2,FF1,FF2)
    !---------------------------------------------------------------------
    ! Computes: 
    ! FF1 = F20 = -u^\dag f1 V - v^\dag f2^T u
    ! FF2 = F02 = -u^T f1^T V^* - v^T f2 u^*
    !---------------------------------------------------------------------
    Implicit None
    Complex(pr), Intent(in)  :: U(iHFB_LEN),V(iHFB_LEN),F1(iQRPA_LEN),F2(iQRPA_LEN)
    Complex(pr), Intent(out) :: FF1(iQRPA_LEN),FF2(iQRPA_LEN)
    Complex(pr), allocatable :: aux1(:),aux2(:)
    Integer(ipr) :: ib,nd,ipa,ipq,ipt,ibHFBi,ibHFBf,ibTHOi,ibTHOf
    !
    If(allocated(aux1)) Deallocate(aux1,aux2)
    Allocate(aux1(iQRPA_LEN),aux2(iQRPA_LEN))
    !
    call qrpa_TriProd('c','n','n',U,F1,V,aux1)
    call qrpa_TriProd('c','t','n',V,F2,U,aux2)
    FF1 = -aux1-aux2
    !
    call qrpa_TriProd('t','*','*',U,F1,V,aux1)
    call qrpa_TriProd('t','c','*',V,F2,U,aux2)
    FF2 = -aux1-aux2
    !
    Deallocate(aux1,aux2)
  End Subroutine qrpa_F20
  !=================================================================================================================================
  ! 
  !=================================================================================================================================
  Subroutine qrpa_F11(V,U,F1,F2,FF1,FF2)
    !---------------------------------------------------------------------
    ! Computes: 
    ! FF1 = F11  = u^\dag f1 u - v^\dag f2^\dag v
    ! FF2 = F11b = u^T f2^* u^*  - v^T f1^T v^*
    !---------------------------------------------------------------------
    Implicit None
    Complex(pr), Intent(in)  :: U(iHFB_LEN),V(iHFB_LEN),F1(iQRPA_LEN),F2(iQRPA_LEN)
    Complex(pr), Intent(out) :: FF1(iQRPA_LEN),FF2(iQRPA_LEN)
    Complex(pr), allocatable :: aux1(:),aux2(:)
    Integer(ipr) :: ib,nd,ipa,ipq,ipt,ibHFBi,ibHFBf,ibTHOi,ibTHOf
    !
    If(allocated(aux1)) Deallocate(aux1,aux2)
    Allocate(aux1(iQRPA_LEN),aux2(iQRPA_LEN))
    !
    call qrpa_TriProd('c','n','n',U,F1,U,aux1)
    call qrpa_TriProd('c','t','n',V,F2,V,aux2)
    FF1 = aux1-aux2
    !
    call qrpa_TriProd('t','n','*',U,F2,U,aux1)
    call qrpa_TriProd('t','t','*',V,F1,V,aux2)
    FF2 = aux1-aux2
!
! MOD HERE !
!call qrpa_TriProd('c','n','n',U,F1,U,aux1)
!call qrpa_TriProd('c','t','n',V,F2,V,aux2)
!FF1=aux1
!FF2=aux2


    !
    Deallocate(aux1,aux2)
  End Subroutine qrpa_F11
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine Calculate_QRPA_H20
    !---------------------------------------------------------------------
    ! calculates H20 using:               Output: 
    !       - HFB matrices                 - H20, H02
    !       - QRPA amplitudes              - H11, H11b (if needed)
    !
    !---------------------------------------------------------------------
   !Use HFBTHO , only : Print_Screen , idebug
    Implicit None
    integer(ipr) :: i !! remove later
    !---------------------------------------------------------------------
    ! calculates h,d,Bd for neutrons and protons
    !---------------------------------------------------------------------
    !
    Call qrpa_DENSIT
    !
    Call qrpa_coulom
    !
    Call qrpa_field 
    !
    ! matrices h1,h2,d,Bd w.r.t basis state indices
    Call qrpa_gamdel
    !
    ! H11, w.r.t quasiparticel indices
   !if(ENABLE_N_EFA) Call qrpa_H11(VqpN,UqpN,hN_c,dN_c,BhN_c,BdN_c,h11N_c,Bh11N_c)
   !if(ENABLE_P_EFA) Call qrpa_H11(VqpP,UqpP,hP_c,dP_c,BhP_c,BdP_c,h11P_c,Bh11P_c)


    !---------------------------------
    ! phonon correction procedure
    !---------------------------------
    if(si.Lt.qrpa_eps) Call qrpa_NH11(VqpN,UqpN,hN_c,dN_c,BhN_c,BdN_c)
    if(si.Lt.qrpa_eps) Call qrpa_PH11(VqpP,UqpP,hP_c,dP_c,BhP_c,BdP_c)


    !---------------------------------------------------------------------
    ! H20,BH20 as h,Bh. Note: subroutines below overwrite h_c,d_c,Bh_c,Bd_c
    !---------------------------------------------------------------------
    Call qrpa_H20(VqpN,UqpN,hN_c,dN_c,BhN_c,BdN_c)
    Call qrpa_H20(VqpP,UqpP,hP_c,dP_c,BhP_c,BdP_c)
    !
    !---------------------------------------------------------------------
    !---------------------------------------------------------------------
    ! add external fields
    !---------------------------------------------------------------------
    !
    ! neutrons
    hN_c=hN_c+FhN_c; BhN_c=BhN_c+BFhN_c
    if(ENABLE_N_EFA) then ; h11N_c = h11N_c+F11N_c;  Bh11N_c = Bh11N_c+BF11N_c; end if
    if(ENABLE_N_EFA.and.USE_ONLY_F11EFA) then ; h11N_c = F11N_c;  Bh11N_c = BF11N_c; end if
    !
    ! protons
    hP_c=hP_c+FhP_c; BhP_c=BhP_c+BFhP_c
    if(ENABLE_P_EFA) then ; h11P_c = h11P_c+F11P_c;  Bh11P_c = Bh11P_c+BF11P_c; end if
    if(ENABLE_P_EFA.and.USE_ONLY_F11EFA) then ; h11P_c = F11P_c;  Bh11P_c = BF11P_c; end if
    !
  End Subroutine Calculate_QRPA_H20
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_ctrace(A,N,ctrace)
    !---------------------------------------------------------------------
    ! block trace - expand A(ip) to A(N,N) 
    !---------------------------------------------------------------------
    Implicit None
    Complex(pr), Intent(in)  :: A(N,N)
    Complex(pr), Intent(out) :: ctrace
    Integer(ipr) :: N,I 
    ctrace = 0.0_pr
    Do I=1,N
       ctrace = ctrace + A(I,I)
    End Do
  End Subroutine qrpa_ctrace
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_broyden(niter)
    !-------------------------------------
    ! Broyden mixing, could be optimized
    !-------------------------------------
   !Use HFBTHO, only : si
    Implicit None
    Integer(ipr), Intent(in)  :: niter
    Integer(ipr) :: id,nd,i,j,k1,k2 

    ! determine broyden array dimension
    id=8*iQRPA_LEN
    if(ENABLE_N_EFA) id = id +4*iQRPA_LEN
    if(ENABLE_P_EFA) id = id +4*iQRPA_LEN
    !
    If(niter.Eq.0) Then
       ! zero iteration
       qrpa_bbroyden='i'; si=1.0_pr;
       If(Allocated(qrpa_broin)) Deallocate(qrpa_broin,qrpa_broout)    
       Allocate(qrpa_broin(id),qrpa_broout(id)); qrpa_broin=0.0_pr
    End If
    !
    ! to real variable
    Do i=1,iQRPA_LEN
      qrpa_broout(i+0*iQRPA_LEN)=Real(XN_c(i),kind=pr)
      qrpa_broout(i+1*iQRPA_LEN)=Imag(XN_c(i))
      qrpa_broout(i+2*iQRPA_LEN)=Real(XP_c(i),kind=pr)
      qrpa_broout(i+3*iQRPA_LEN)=Imag(XP_c(i)) 
      qrpa_broout(i+4*iQRPA_LEN)=Real(YN_c(i),kind=pr)
      qrpa_broout(i+5*iQRPA_LEN)=Imag(YN_c(i))
      qrpa_broout(i+6*iQRPA_LEN)=Real(YP_c(i),kind=pr)
      qrpa_broout(i+7*iQRPA_LEN)=Imag(YP_c(i))
      if(ENABLE_N_EFA) then
        qrpa_broout(i+8*iQRPA_LEN)=Real(P1N_c(i),kind=pr)
        qrpa_broout(i+9*iQRPA_LEN)=Imag(P1N_c(i))
        qrpa_broout(i+10*iQRPA_LEN)=Real(P2N_c(i),kind=pr)
        qrpa_broout(i+11*iQRPA_LEN)=Imag(P2N_c(i))
      end if
      if(ENABLE_P_EFA .and. .not.ENABLE_N_EFA) then
        qrpa_broout(i+8*iQRPA_LEN)=Real(P1P_c(i),kind=pr)
        qrpa_broout(i+9*iQRPA_LEN)=Imag(P1P_c(i))
        qrpa_broout(i+10*iQRPA_LEN)=Real(P2P_c(i),kind=pr)
        qrpa_broout(i+11*iQRPA_LEN)=Imag(P2P_c(i))
      end if
      if(ENABLE_P_EFA .and. ENABLE_N_EFA) then
        qrpa_broout(i+12*iQRPA_LEN)=Real(P1P_c(i),kind=pr)
        qrpa_broout(i+13*iQRPA_LEN)=Imag(P1P_c(i))
        qrpa_broout(i+14*iQRPA_LEN)=Real(P2P_c(i),kind=pr)
        qrpa_broout(i+15*iQRPA_LEN)=Imag(P2P_c(i))
      end if
    End Do
    !
    ! Broyden
    Call qrpa_broyden_method(id,qrpa_broout,qrpa_broin,qrpa_alphamix,si,niter,qrpa_nbroyden,qrpa_bbroyden)
    If(niter.Eq.0) si=1.0_pr
    !
    Do i=1,iQRPA_LEN
      XN_c(i)=Cmplx(qrpa_broin(i+0*iQRPA_LEN),qrpa_broin(i+1*iQRPA_LEN),kind=pr)
      XP_c(i)=Cmplx(qrpa_broin(i+2*iQRPA_LEN),qrpa_broin(i+3*iQRPA_LEN),kind=pr)
      YN_c(i)=Cmplx(qrpa_broin(i+4*iQRPA_LEN),qrpa_broin(i+5*iQRPA_LEN),kind=pr)
      YP_c(i)=Cmplx(qrpa_broin(i+6*iQRPA_LEN),qrpa_broin(i+7*iQRPA_LEN),kind=pr)
      if(ENABLE_N_EFA) then
        P1N_c(i)=Cmplx(qrpa_broin(i+8*iQRPA_LEN),qrpa_broin(i+9*iQRPA_LEN),kind=pr)
        P2N_c(i)=Cmplx(qrpa_broin(i+10*iQRPA_LEN),qrpa_broin(i+11*iQRPA_LEN),kind=pr)
      end if
      if(ENABLE_P_EFA .and. .not.ENABLE_N_EFA) then
        P1P_c(i)=Cmplx(qrpa_broin(i+8*iQRPA_LEN),qrpa_broin(i+9*iQRPA_LEN),kind=pr)
        P2P_c(i)=Cmplx(qrpa_broin(i+10*iQRPA_LEN),qrpa_broin(i+11*iQRPA_LEN),kind=pr)
      end if
      if(ENABLE_P_EFA .and. ENABLE_N_EFA) then
        P1P_c(i)=Cmplx(qrpa_broin(i+12*iQRPA_LEN),qrpa_broin(i+13*iQRPA_LEN),kind=pr)
        P2P_c(i)=Cmplx(qrpa_broin(i+14*iQRPA_LEN),qrpa_broin(i+15*iQRPA_LEN),kind=pr)
      end if
    End Do
  End Subroutine qrpa_broyden
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_broyden_method(N,vout,vin,alpha,si,iter,M,bbroyden)
    !---------------------------------------------------------------------
    ! Modified Broyden's method: D.D.Johnson, PRB 38, 12807 (1988)
    ! Adopted from: (C) 2001 PWSCF group
    ! Input :
    !  N      dimension of arrays vin,vout
    !  vin    outpu at previous iteration
    !  vout   output at current iteration
    !  alpha  mixing factor (0 < alpha <= 1)
    !  iter   current iteration number
    !  M      number of iterations in Broyden history
    !         (M=0 Linear mixing)
    ! Output:
    !  si     MaxVal(Abs(vout-vin))
    !  vin    Broyden/Linear mixing result
    !  vout   vout-vin
    !  bbroyden='B' Broyden mixing, ='L' Linear mixing
    !---------------------------------------------------------------------
   !Use HFBTHO, only : Print_Screen, idebug
    Implicit None
    Integer(ipr),    Intent(In)     :: N,iter,M
    Real(pr),        Intent(In)     :: alpha
    Real(pr),        Intent(Out)    :: si
    Character(1),    Intent(Out)    :: bbroyden
    Real(pr),        Intent(InOut)  :: vout(N),vin(N)  
    Integer(ipr)                    :: i,j,iter_used,ipos,inext
    Integer(ipr), Allocatable,Save  :: iwork(:)
    Real(pr),    Allocatable, Save  :: beta(:,:),work(:)
    Real(pr),    Allocatable, Save  :: df(:,:),dv(:,:),curv(:)
    Real(pr),                 Save  :: w0
    Real(pr)                        :: DDOT,DNRM2,normi,gamma,curvature,sf
    !
    sf=-1.0_pr; Call DAXPY(N,sf,vin,1,vout,1)         ! vout = vout - vin
    si=Maxval(Abs(vout))    
    !---------------------------------------------------------------------
    ! Linear mixing
    !---------------------------------------------------------------------
    If(M.Eq.0.Or.iter.Eq.0) Then
       bbroyden='L'; Call DAXPY(N,alpha,vout,1,vin,1) ! vin = vin + alpha*vout         
      !If(Print_Screen.And.IDEBUG.Gt.10) Then 
      !   If(iter.Eq.0) Write(6,*) '  Linear mixing (alpha) : ',alpha
      !Endif
       Return
    End If
    !---------------------------------------------------------------------
    ! Broyden mixing
    !---------------------------------------------------------------------
    bbroyden='B'
    iter_used=Min(iter-1,M); ipos=iter-1-((iter-2)/M)*M; inext=iter-((iter-1)/M)*M
    If (iter.Eq.1) Then       
       w0=0.010_pr
       If(Allocated(curv)) Deallocate(curv,df,dv,beta,work,iwork)
       Allocate(curv(N),df(N,M),dv(N,M),beta(M,M),work(M),iwork(M))
     ! If(Print_Screen.And.IDEBUG.Gt.10) Then 	   
     !    Write(6,'(a,i3,3(2x,f18.8),a)') '   Broyden mixing (M,alpha,w0,mem) : '  &
     !         ,M,alpha,w0,(2*N*M+N)*8._pr/1.e6,' MB'
     ! End If
    Else
       df(:,ipos)=vout(:)-df(:,ipos) 
       dv(:,ipos)= vin(:)-dv(:,ipos)
       Normi=1.0_pr/DNRM2(N,df(1,ipos),1)
       Call dscal(N,Normi,df(1,ipos),1)
       Call dscal(N,Normi,dv(1,ipos),1)
    Endif
    Do i=1,iter_used
       Do j=i+1,iter_used
          beta(i,j)=DDOT(N,df(1, j),1,df(1,i),1)
       Enddo
       beta(i,i)= w0*w0  + 1.0_pr 
    Enddo
    Call DSYTRF('U',iter_used,beta,M,iwork,work,M,i)
    If(i.Ne.0) Stop '  In Broyden: info at DSYTRF '
    Call DSYTRI('U',iter_used,beta,M,iwork,work,i)
    If(i.Ne.0) Stop '  In Broyden: info at DSYTRI '
    Do i=1,iter_used
       Do j=i+1,iter_used
          beta(j,i)=beta(i,j)
       Enddo
       work(i)=DDOT(N,df(1,i),1,vout,1)      
    Enddo
    curv=alpha*vout
    Do i=1,iter_used
       gamma=0.0_pr
       Do j=1,iter_used
          gamma=gamma+beta(j,i)*work(j)
       Enddo
       curv=curv-gamma*(dv(:,i)+alpha*df(:,i))
    Enddo
    Call dcopy(N,vout,1,df(1,inext),1)
    Call dcopy(N,vin ,1,dv(1,inext),1)
    sf=+1.0_pr; Call DAXPY(N,sf,curv,1,vin,1)      
  End Subroutine qrpa_broyden_method
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Elemental Function qrpa_electric_field(t,x,y,z)
    !------------------------------------------------------------------
    ! For electric transitions operator
    ! called in qrpa_f_gamdel: qrpa_f_gamdel(t,fl(ihli),fh(ihli)) 
    ! t=0:neutrons, 1:protons
    ! Memo: fl(ihli)^2+fh(ihli)^2=x^2+y^2+z^2
    ! Note: Here in HFBTHO implementation y is the coordinate to radial direction and x = 0
    !------------------------------------------------------------------
   !Use HFBTHO, only : npr
    Implicit None
    Integer(ipr), Intent(in) :: t
    Real(pr),     Intent(in) :: x,y,z
    Real(pr)                 :: t_phase,qrpa_electric_field,qrpa_operator_strenght,qrpa_operator_field
    !
    qrpa_operator_strenght = 0.0_pr
    qrpa_operator_field = 0.0_pr
    !
    select case(L_qrpa_responce)
    case(0)
    !----------------------------------------------------    
    ! Monopole operator F_00
    !----------------------------------------------------    
      qrpa_operator_field    = (x**2+y**2+z**2) 
    !
    case(1)
    !----------------------------------------------------    
    ! Dipole operator F_1K
    !----------------------------------------------------    
      select case(K_qrpa_responce)
        case(-1)
          qrpa_operator_field    = 0.5_pr*Sqrt(1.5_pr/Pi)*y  !! here y <--> x
        case(0)
          qrpa_operator_field    = 0.5_pr*Sqrt(3.0_pr/Pi)*z
        case(1)
          qrpa_operator_field    = -0.5_pr*Sqrt(1.5_pr/Pi)*y  !! here y <--> x
        case default
      end select
    !
    case(2)
    !----------------------------------------------------    
    ! Quadrupole operator F_2K
    !----------------------------------------------------    
      select case(K_qrpa_responce)
        case(0)
          qrpa_operator_field    =  0.25_pr*Sqrt(5.0_pr/Pi)*(2.0_pr*z**2-(x**2+y**2))
        case(1)
          qrpa_operator_field    =  0.50_pr*Sqrt(7.5_pr/Pi)*z*y
        case(-1)
          qrpa_operator_field    =  0.50_pr*Sqrt(7.5_pr/Pi)*z*y
        case(2)
          qrpa_operator_field    = -0.25_pr*Sqrt(7.5_pr/Pi)*y*y
        case(-2)
          qrpa_operator_field    = -0.25_pr*Sqrt(7.5_pr/Pi)*y*y
        case default
      end select
    !----------------------------------------------------    
    ! Octupole operator F_3K
    !----------------------------------------------------    
    case(3)
       select case(K_qrpa_responce)
        case(0)
           qrpa_operator_field    = 0.25_pr*Sqrt(7.0_pr/Pi)*(2.0_pr*z**2 -3.0_pr*(x**2+y**2))*z
        case(1)
           qrpa_operator_field    = -0.125_pr*Sqrt(21.0_pr/Pi)*(4.0_pr*z**2 - y**2 )*y
        case(2)
           qrpa_operator_field    = 0.25_pr*Sqrt(105.0_pr/(2.0_pr*Pi))*(x**2+y**2)*z
        case(3)
           qrpa_operator_field    = -0.125_pr*Sqrt(35.0_pr/Pi)*y*y*y
         case(-3)
           qrpa_operator_field    = -0.125_pr*Sqrt(35.0_pr/Pi)*y*y*y
         case default
       end select

    case(4)
       select case(K_qrpa_responce)
         case(0)
           qrpa_operator_field = 0.1875_pr*Sqrt(1.0_pr/Pi)*(35*z*z*z*z-30*z*z*(x**2+y**2+z**2)+3*(x**2+y**2+z**2)**2)
         case(1)
           qrpa_operator_field = 0.375_pr*Sqrt(5.0_pr/Pi)*y*(7.0_pr*z*z*z - 3.0_pr*z*(x**2+y**2+z**2))
         case(2)
           qrpa_operator_field = 0.375_pr*Sqrt(2.5_pr/Pi)*y*y*(7.0_pr*z*z - (x**2+y**2+z**2))
         case(3)
           qrpa_operator_field = -0.375_pr*Sqrt(35.0_pr/Pi)*y*y*y*z
         case(4)
           qrpa_operator_field = 0.1875_pr*Sqrt(17.5_pr/Pi)*y*y*y*y
         case default
       end select

    case(5)
       select case(K_qrpa_responce)
         case(0)
            qrpa_operator_field = (1.0_pr/16.0_pr)*Sqrt(11.0_pr/Pi)*(63.0_pr*z**5 & 
                 - 70.0_pr*z**3*(x**2+y**2+z**2) + 15.0_pr*z*(x**2+y**2+z**2)**2)
         case(1)
            qrpa_operator_field = -(1.0_pr/16.0_pr)*Sqrt(82.5_pr/Pi)*y*(21.0_pr*z**4 - 14.0_pr*z**2*(x**2+y**2+z**2) & 
                 +  (x**2+y**2+z**2)**2)
         case(2)
            qrpa_operator_field = (1.0_pr/8.0_pr)*Sqrt(577.5_pr/Pi)*y**2*(3.0_pr*z**3 - z*(x**2+y**2+z**2))
         case(3)
            qrpa_operator_field = -(1.0_pr/32.0_pr)*Sqrt(385.0_pr/Pi)*y**3*(9.0_pr*z**2 - (x**2+y**2+z**2))
         case(4)
            qrpa_operator_field = (3.0_pr/16.0_pr)*Sqrt(192.5_pr/Pi)*y**4*z
         case(5)
            qrpa_operator_field = (3.0_pr/32.0_pr)*Sqrt(77.0_pr/Pi)*y**5
         case default
       end select
    case(6)
       select case(K_qrpa_responce)
         case(0)
            qrpa_operator_field = (1.0_pr/32.0_pr)*Sqrt(13.0_pr/Pi)*(231.0_pr*z**6 - 315.0_pr*z**4*(x**2+y**2+z**2) &
                   +105.0_pr*z*z*(x**2+y**2+z**2)**2 -5.0_pr*(x**2+y**2+z**2)**3)
         case(1)
            qrpa_operator_field = +(1.0_pr/16.0_pr)*Sqrt(136.5_pr/Pi)*y*(33.0_pr*z**5 - 30.0_pr*z**3*(x**2+y**2+z**2) &
                   + 5.0_pr*z*(x**2+y**2+z**2)**2)
         case(2)
            qrpa_operator_field = (1.0_pr/64.0_pr)*Sqrt(1365/Pi)*y*y*(33.0_pr*z*z*z*z - 18.0_pr*z*z*(x**2+y**2+z**2)  &
                   + (x**2+y**2+z**2)**2)
         case(3) 
            qrpa_operator_field = -(1.0_pr/32.0_pr)*Sqrt(1365.0_pr/Pi)*y**3*(11.0_pr*z**3 - 3.0_pr*z*(x**2+y**2+z**2))
         case(4)
            qrpa_operator_field = (3.0_pr/32.0_pr)*Sqrt(45.5_pr/Pi)*y**4*(11.0_pr*z**2 - (x**2+y**2+z**2))
         case(5)
            qrpa_operator_field = (3.0_pr/32.0_pr)*Sqrt(1001.0_pr/Pi)*y**5*z
         case(6)
            qrpa_operator_field = -(1.0_pr/64.0_pr)*Sqrt(3003.0_pr/Pi)*y**6
         case default
       end select


     end select
    !----------------------------------------------------    
    ! T of the operator. Note: convention to the old FAM
    ! module has been changed. The operator is no longer
    ! multiplied by charge e
    ! Note: with EFA, npr(3) is not the total A, but 
    ! npr(1)+npr(2) is correct
    !----------------------------------------------------    
    select case(T_qrpa_responce)
    case(0)
      ! isoscalar F = e (Z/A) * opr_i
      If(t.Eq.0) Then
      ! neutrons
        qrpa_operator_strenght = + (Dble(npr(2))/Dble(npr(1)+npr(2)))
      Else
      ! Protons
        qrpa_operator_strenght = + (Dble(npr(2))/Dble(npr(1)+npr(2)))
      Endif 
    case(1)
      ! isovector
      If(t.Eq.0) Then
      ! neutrons
        qrpa_operator_strenght = + (Dble(npr(2))/Dble(npr(1)+npr(2)))
      Else
      ! Protons
        qrpa_operator_strenght = - (Dble(npr(1))/Dble(npr(1)+npr(2)))
      Endif 
    case(2)
      ! bare charge
      If(t.Eq.0) Then
      ! neutrons
        qrpa_operator_strenght = 0.0_pr
      Else
      ! Protons
        qrpa_operator_strenght = +1.0_pr
      Endif 
    case default
    end select
    !
    if( K_qrpa_responce .ne. 0 ) qrpa_operator_field = qrpa_operator_field *(1.0_pr/Sqrt(2.0_pr))  !! simplex-y
    qrpa_electric_field    = qrpa_operator_strenght*qrpa_operator_field
  End Function qrpa_electric_field
  !=================================================================================================================================
  ! 
  !=================================================================================================================================
  Elemental Function qrpa_spin_field(t,x,y,z)
    !
    ! constant external field for spin operator. Depends only on 
    ! particle type and T_qrpa_responce
   !Use HFBTHO, only : npr
    Implicit None
    Integer(ipr), Intent(in) :: t
    Real(pr),     Intent(in) :: x,y,z
    Real(pr)                 :: qrpa_spin_field
    !
    qrpa_spin_field = 0.0_pr
    select case(T_qrpa_responce)
    case(0)
      ! isoscalar type
      If(t.Eq.0) Then
      ! neutrons
        qrpa_spin_field = + (Dble(npr(2))/Dble(npr(1)+npr(2)))
      Else
      ! Protons
        qrpa_spin_field = + (Dble(npr(2))/Dble(npr(1)+npr(2)))
      Endif 
    case(1)
      ! isovector type
      If(t.Eq.0) Then
      ! neutrons
        qrpa_spin_field = + (Dble(npr(2))/Dble(npr(1)+npr(2)))
      Else
      ! Protons
        qrpa_spin_field = - (Dble(npr(1))/Dble(npr(1)+npr(2)))
      Endif
    case(2)
      ! with bare g-factors
      If(t.Eq.0) Then
      ! neutrons
        qrpa_spin_field = -3.826_pr
      Else
      ! Protons
        qrpa_spin_field = +5.586_pr
      Endif  
    case default
    end select
    !
    qrpa_spin_field = qrpa_spin_field*Sqrt(0.75_pr/Pi)
  End Function qrpa_spin_field
  !=================================================================================================================================
  ! 
  !=================================================================================================================================
  Elemental Function qrpa_l_field(t,x,y,z)
    !
    ! constant external field for l operator. Depends only on 
    ! particle type and T_qrpa_responce
   !Use HFBTHO, only : npr 
    Implicit None
    Integer(ipr), Intent(in) :: t
    Real(pr),     Intent(in) :: x,y,z
    Real(pr)                 :: qrpa_l_field
    !
    qrpa_l_field = 0.0_pr
    select case(T_qrpa_responce)
    case(0)
      ! isoscalar type
      If(t.Eq.0) Then
      ! neutrons
        qrpa_l_field = + (Dble(npr(2))/Dble(npr(1)+npr(2)))
      Else
      ! Protons
        qrpa_l_field = + (Dble(npr(2))/Dble(npr(1)+npr(2)))
      Endif 
    case(1)
      ! isovector type
      If(t.Eq.0) Then
      ! neutrons
        qrpa_l_field = + (Dble(npr(2))/Dble(npr(1)+npr(2)))
      Else
      ! Protons
        qrpa_l_field = - (Dble(npr(1))/Dble(npr(1)+npr(2)))
      Endif
    case(2)
      ! with bare g-factors
      If(t.Eq.0) Then
      ! neutrons
        qrpa_l_field = 0.0_pr
      Else
      ! Protons
        qrpa_l_field = +1.0_pr
      Endif  
    case default
    end select
    !
    qrpa_l_field = qrpa_l_field*Sqrt(0.75_pr/Pi)*2.0_pr/(real(abs(L_qrpa_responce),kind=pr)+1.0_pr)
    !
  End Function qrpa_l_field
  !=================================================================================================================================
  ! 
  !=================================================================================================================================
  Subroutine qrpa_HFBTHO_solver
    !----------------------------------------------------------------------------------------------------
    ! qrpa solver. Input parameters:
    !   integer:   max_iter_qrpa,qrpa_nbroyden
    !   real:      qrpa_eps,qrpa_alphamix
    !   complex:   qrpa_eta,qrpa_omega
    !----------------------------------------------------------------------------------------------------
    Use HFBTHO, only : iqrpa
    Implicit None
    Logical      :: file_exists
    Integer(ipr) :: i,j,k1,k2,kk1,kk2,imen,ih,il,ib,ibx,nd,nza,nra,nla,nsa,iw       
    Real(pr)     :: snorm_xN,snorm_yN,snorm_N,snorm_xP,snorm_yP,snorm_P
    Real(pr)     :: xymix,DZNRM2,xynorm_N,xynorm_P
    Real(pr)     :: Rbox,eFermiN
    Real(pr)     :: time1,time2
    Real(pr)     :: strN,strP
    Character(LEN=100) :: fname
    Character(LEN=6)   :: tname,kname,mname,opname
    !----------------------------------------------------------------------------------------------------
    ! QRPA INITIALIZATION
    !----------------------------------------------------------------------------------------------------
    If(iqrpa.Eq.0) Then
       iqrpa=1
       Call qrpaINI_CUT_ALLOCATE      

       Do iw=lout,lfile
         Write(iw,*) '          < Starting FAM calculation >'
        !Write(iw,*) '  FAM QRPA version: ',qrpa_Version
         Write(iw,*) '  '
         Write(iw,*) '  Time-odd coupling constants:'
         Write(iw,*) '  -----------------------------------'
         Write(iw,'("   Css(0)=  ",g26.18,"; Css(1)=  ",g26.18)') Css
         Write(iw,'("   CDss(0)= ",g26.18,"; CDss(1)= ",g26.18)') CDss
         Write(iw,'("   Cjj(0)=  ",g26.18,"; Cjj(1)=  ",g26.18)') Cjcjc
         Write(iw,'("   CsNj(0)= ",g26.18,"; CsNj(1)= ",g26.18)') Csnablajc
         Write(iw,'("   CsDs(0)= ",g26.18,"; CsDs(1)= ",g26.18)') CDeltass
         Write(iw,'("   CsT(0)=  ",g26.18,"; CsT(1)=  ",g26.18)') CsT
         Write(iw,*) '  '
         select case(QRPA_Operator_type)
         case(1)
           Write(iw,*) '  Electric transition operator used'
         case(2)
           Write(iw,*) '  Spin transition operator used'
         case(3)
           Write(iw,*) '  Magnetic transition operator used'
         case(4)
           Write(iw,*) '  Linear momentum operator used'
         case(5)
           Write(iw,*) '  Coordinate operator used'
         case(6)
           Write(iw,*) '  Total angular momentum J11 operator used'
         case(7)
           Write(iw,*) '  Total angular momentum J_y operator used'
         case default
           STOP "Unknown operator type"
         end select

         If(USE_NORMALHO_K0 .and. K_qrpa_responce.eq.0) Write(iw,*) '  Using normal harmonic oscillator basis.'
         if(USE_ONLY_F20) then
          Write(iw,*) '  '
          Write(iw,*) '  >>  WARNING! WARNING!  H20, H02 and all induced fields are set to zero!  WARNING! WARNING!  <<'
          Write(iw,*) '  '
         end if 

         if(USE_ONLY_F11EFA.and.(ENABLE_N_EFA.or.ENABLE_P_EFA)) then
          Write(iw,*) '  '
          Write(iw,*) '  >>  WARNING! WARNING!  Using only F11 part with EFA!  WARNING! WARNING!  <<'
          Write(iw,*) '  '
         end if 

         Write(iw,*) '  '

         if(ENABLE_N_EFA) Write(iw,*) '  Enabled neutron equal filling approximation due to q.p. blocking '
         if(ENABLE_P_EFA) Write(iw,*) '  Enabled proton equal filling approximation due to q.p. blocking '
         if(ENABLE_N_EFA.or.ENABLE_P_EFA) Write(iw,*) '  '

         Write(iw,*) '  QRPA CALCULATIONS '
         Write(iw,*) '  DIMENSION OF U,V iHFB_LEN:',iHFB_LEN
         Write(iw,*) '  DIMENSION OF Eqp iHFB_NQP:',iHFB_NQP
         Write(iw,*) '  MAX ITERATIONS           :',max_iter_qrpa   
         Write(iw,*) '  BROYDEN HISTORY          :',qrpa_nbroyden
         Write(iw,*) '  MIXING PARAMETER         :',qrpa_alphamix
         Write(iw,*) '  QRPA STOP CRITERIA       :',qrpa_eps   
       End Do
    Endif
    !----------------------------------------------------------------------------------------------------
    ! Set the file name for results
    !----------------------------------------------------------------------------------------------------
    select case(T_qrpa_responce)
    case(0)
      tname = 'is'
    case(1)
      tname = 'iv'
    case(2)
      tname = 'bare'
    case default
      STOP "Unknown T_qrpa_responce"
    end select
    !
    select case(QRPA_Operator_type)
    case(1)
      opname = 'E'
    case(2)
      opname = 'S'
    case(3)
      opname = 'M'
    case(4)
      opname = 'P'
    case(5)
      opname = 'C'
    case(6)
      opname = 'J11'
    case(7)
      opname = 'J_y'
    case default
      STOP "Unknown operator type"
    end select

 !  write (kname, "(I1,I1)") L_qrpa_responce,K_qrpa_responce
 !  write (mname, "(I6)") iam_mpi
 !  fname = trim('qrpa_results_' // trim(tname) // '_' //  trim(opname) //  trim(kname) // '_' // trim(adjustl(mname)) // '.dat')
    !----------------------------------------------------------------------------------------------------
    ! QRPA ITERATIONS 
    !----------------------------------------------------------------------------------------------------
    Do iter_qrpa=0,max_iter_qrpa
       Call Cpu_time(time1)                                          ! qrpa time per iteration
       If(iter_qrpa.Eq.0) Then
          Call qrpa_green_function(qrpa_omega)                      ! calculate green functions		
          hN_c=zero;  hP_c=zero;  BhN_c=zero;  BhP_c=zero           ! start with X=Y=0
          si=1.0_pr                                                 ! defect on zero iteration 
          Do iw=lout,lfile                                          ! print header information
             Write(iw,*) 
             Write(iw,*) '  OMEGA               :',qrpa_omega
             Write(iw,*) '  NEUTRON POLE(-)     :',one/GreenNminus(gfpoles_imN)
             Write(iw,*) '  PROTON  POLE(-)     :',one/GreenPminus(gfpoles_imP)
             Write(iw,*) '  T,L,K RESP.OPERATOR :',T_qrpa_responce,L_qrpa_responce,K_qrpa_responce             
             select case(QRPA_Operator_type)
             case(1)
               if(L_qrpa_responce.eq.0) then
                 Write(iw,*) '  Transition strength in units of [e^2 MeV^{-1} fm^{4}]'
               else
                 Write(iw,*) '  Transition strength in units of [e^2 MeV^{-1} fm^{2L}]'
               end if
             case(2,3)
                 if(L_qrpa_responce.eq.1) Write(iw,*) '  Transition strength in units of [(\mu_N/c)^2 MeV^{-1} ]'
             case(4)
                 if(L_qrpa_responce.eq.1) Write(iw,*) '  Transition strength in units of [ fm^{-2} MeV^{-1} ]'        !krisa! 4 is P operator
             case(5)
                 if(L_qrpa_responce.eq.1) Write(iw,*) '  Transition strength in units of [ fm^2 MeV^{-1} ]'           !krisa! 5 is Q operator
             case(6,7)
                 if(L_qrpa_responce.eq.1) Write(iw,*) '  Transition strength in units of [ MeV^{-1} ]'                !krisa! 6 is J11 and 7 is J_y
             end select
             Write(iw,'(a,f7.3,3(a,i3),(a,g14.8),(a,2g14.8),a)') &
                  '  |QRPA> iterations... b0=',b0,', Nsh=',n00,', N=',npr(1),', Z=',npr(2), &
                  ', Eta=',Real(qrpa_eta,kind=pr),', Omega=',qrpa_omega
             Write(iw,'(2x,154(''-''))')
             Write(iw,'(a,a)') '    i        si    mix        reSn         reSp         reSt         imSn         imSp         imSt', &
                               '           TrStrn       TrStrp       TrStrt          Time'
             Write(iw,'(2x,154(''-''))')
          End Do
       Endif
       Call qrpa_xy                                                        ! calculate X,Y from h,Bh
       Call qrpa_broyden(iter_qrpa)                                        ! use Broyden
       Call Calculate_QRPA_H20                                             ! calculate h,Bh from X,Y 

       Call qrpa_strength(XN_c,YN_c,P1N_c,P2N_c,FhN_c,BFhN_c,F11N_c,BF11N_c, &
                          ENABLE_N_EFA,ImStrengthN,ReStrengthN)                   ! qrpa strenght function neutrons, imag and real parts
       Call qrpa_strength(XP_c,YP_c,P1P_c,P2P_c,FhP_c,BFhP_c,F11P_c,BF11P_c, &
                          ENABLE_P_EFA,ImStrengthP,ReStrengthP)                   ! qrpa strenght function protons, imag and real parts
       !
       ! transition strength
       ReTransStrN = -ReStrengthN/pi
       ReTransStrP = -ReStrengthP/pi
       ImTransStrN = -ImStrengthN/pi
       ImTransStrP = -ImStrengthP/pi
       !		
       Call Cpu_time(time2); time2=time2-time1                             ! time per iteration 
       If(iter_qrpa.Lt.500.Or.si.Lt.qrpa_eps*10.0_pr) Then
          Do iw=lout,lfile
             Write(iw,'(i4,a,1x,f12.8,f5.2," | ",6(1x,g12.6)," | ",3(1x,g12.6)," | ",f7.2)') &
                  iter_qrpa,qrpa_bbroyden,si,qrpa_alphamix,&
                  ReStrengthN,ReStrengthP,ReStrengthN+ReStrengthP,ImStrengthN,ImStrengthP,ImStrengthN+ImStrengthP,&
                  ImTransStrN,ImTransStrP,ImTransStrN+ImTransStrP,time2     !krisa! new variables introduced, snz_qrpa(1:2) variables were removed
          End Do
       Endif
       If(si.Lt.qrpa_eps) Exit   ! Exit if converges
    End Do ! iter_qrpa 
    !----------------------------------------------------------------------------------------------------
    ! APPEND RESULTS TO qrpa_results.dat FILE 
    !----------------------------------------------------------------------------------------------------
   !Inquire(FILE=fname, EXIST=file_exists)   
   !Open (101,file=fname,position='append')
   !If(.Not.file_exists) Write(101,'(20a)') '!  N 1 Z 2 ITER 3  ', &                                !krisa! new strength variable names introduced
   !     ' si 4      ReOmega 5      ImOmega 6      reStrN 7       reStrP  8      ', &
   !     'reStrT 9       imStrN 10      imStrP 11      imStrT 12      TrnStrN 13     TrnStrP 14     TrnStrT 15'
   !Write(101,'(1x,3(1x,I3),1x,f9.6,200(f15.6))') npr(1),npr(2),iter_qrpa,si,qrpa_omega,               &
   !     ReStrengthN,ReStrengthP,(ReStrengthN+ReStrengthP),ImStrengthN,ImStrengthP,(ImStrengthN+ImStrengthP),    &
   !     ImTransStrN,ImTransStrP,(ImTransStrN+ImTransStrP)
   !Close(101)
    !
    !----------------------------------------------------------------------------------------------------
    ! Write amplitudes on file
    !----------------------------------------------------------------------------------------------------
    if(qrpa_file_io.eq.1) call write_amplitudes
    call write_amplitudes
    !
   
  End Subroutine qrpa_HFBTHO_solver
 !=================================================================================================================================
  !
  ! Subroutine which writes the FAM-QRPA amplitudes on disk in
  ! binary format. The file content is following:
  !
  ! T_qrpa_responce,L_qrpa_responce,K_qrpa_responce,P_qrpa_response,QRPA_Operator_type,iQRPA_LEN  : 6 integers
  ! qrpa_omega                                                                                    : 1 complex
  ! ImStrengthN,ImStrengthP,ReStrengthN,ReStrengthP                                               : 4 reals
  ! XN_c(i),YN_c(i),XP_c(i),YP_c(i)                                                               : 4*iQRPA_LEN complex (i=1 to iQRPA_LEN)
  !
  !=================================================================================================================================

  Subroutine write_amplitudes
    Implicit None
    Integer(ipr) :: i
    !--------------------------------
    integer(kind = 4) :: qqi,qqj  
    Real :: imag_qrpa_omega, real_qrpa_omega
    Complex(pr), target :: CQRPA_XN_c(iQRPA_LEN), CQRPA_YN_c(iQRPA_LEN), CQRPA_XP_c(iQRPA_LEN), CQRPA_YP_c(iQRPA_LEN)
    Real(pr) :: Norm_N, Norm_P,Normfactor2               !!!Q: Loca0  
   !open(46,file=trim(adjustl(qq_io_path))//'46_amplitude',access='append')
    open(46,file='46_amplitude',access='append')
    open(47,file='47_optimize',access='append')
    !---------------------------------

    open(20,file="amplitudes.dat",form='unformatted',status='replace',action='write')

    Write(20) T_qrpa_responce,L_qrpa_responce,K_qrpa_responce,P_qrpa_response,QRPA_Operator_type,iQRPA_LEN
    Write(20) ENABLE_N_EFA,ENABLE_P_EFA
    Write(20) qrpa_omega
    Write(20) ImStrengthN,ImStrengthP,ReStrengthN,ReStrengthP
    Do i = 1,iQRPA_LEN
      Write(20) XN_c(i),YN_c(i),XP_c(i),YP_c(i)
    End Do

    !--------------------------------
      Norm_N = 0
      Norm_P = 0
      !
      imag_qrpa_omega = aimag(qrpa_omega)
      real_qrpa_omega = real(qrpa_omega)
      ! 
      !     
      Normfactor2 = imag_qrpa_omega/( Pi*ImTransStrN + Pi*ImTransStrP )
      Do qqi = 1, iQRPA_LEN                             !!!Q: make QRPA X and Y complex
         !-------------------------------------------QQst    modify17
         CQRPA_XN_c(QQi)= sqrt(Normfactor2)*XN_c(QQi)                            
         CQRPA_YN_c(QQi)= sqrt(Normfactor2)*YN_c(QQi)
         CQRPA_XP_c(QQi)= sqrt(Normfactor2)*XP_c(QQi)
         CQRPA_YP_c(QQi)= sqrt(Normfactor2)*YP_c(QQi)
         Norm_N = Norm_N   & 
                   +  real(CQRPA_XN_c(QQi))* real(CQRPA_XN_c(QQi)) + aimag(CQRPA_XN_c(QQi))* aimag(CQRPA_XN_c(QQi)) &
                   -  real(CQRPA_YN_c(QQi))* real(CQRPA_YN_c(QQi)) - aimag(CQRPA_YN_c(QQi))* aimag(CQRPA_YN_c(QQi))
         Norm_P = Norm_P   & 
                   +  real(CQRPA_XP_c(QQi))* real(CQRPA_XP_c(QQi)) + aimag(CQRPA_XP_c(QQi))* aimag(CQRPA_XP_c(QQi)) &
                   -  real(CQRPA_YP_c(QQi))* real(CQRPA_YP_c(QQi)) - aimag(CQRPA_YP_c(QQi))* aimag(CQRPA_YP_c(QQi))
      end do
      write(46,*) real_qrpa_omega,',', Norm_N+Norm_P,',', Normfactor2,',', ImTransStrN + ImTransStrP 
 

!! ADD HERE ALSO EFA AMPLITUDES LATER

       !------------------------------------
       ! AB matrix 
       !--------------
       if(allocated(WX_sum_nd_colm1)) deallocate(WX_sum_nd_colm1)
       allocate(WX_sum_nd_colm1(iHFB_NB))

          !
       WX_sum_nd=0
       WX_sum_nd_colm1=0
       do qqj=1,iHFB_NB
          WX_sum_nd=WX_sum_nd + id(qqj)
          WX_sum_nd_colm1(qqj)=WX_sum_nd - id(qqj)
       end do

       !-----------------------------------------------
       !translate index of X and Y to /mu and /nu, a table stored in AB_row and AB_col.
       !----------------------------------
       if(allocated(AB_row)) deallocate(AB_row, AB_col)
       allocate(AB_row(size(XN_c)), AB_col(size(XN_c)))
          !
          qqk = 0

          ct1=1
          ct3=0
          inbcc=0
          inbrc=0
          inbrr=0
          do qqj=1,iHFB_NB
             do qqi=1,IHFB_NB
                if(QQ_select(qqi,qqj)==1)  then
                   inbcc=WX_sum_nd_colm1(qqj)+1
                   inbrc=WX_sum_nd_colm1(qqi)
                   ct3=ct3+id(qqj)*id(qqi)
                      do ct2=ct1, ct3 
                          inbrc=inbrc+1 !row index
                          inbrr=inbrr+1 !col index
                          qqk = qqk +1
                          AB_row(qqk) = inbrc
                          AB_col(qqk) = inbcc
                         !write(*,*) 'qqk', qqk, 'sizeX', size(XN_c), AB_row(qqk), AB_col(qqk),'id', qqi, id(qqi), qqj, id(qqj)
                        if(inbrr.eq. id(qqi)) then
                           inbcc=inbcc+1        !!!Q: col increase by 1
                           inbrc=WX_sum_nd_colm1(qqi)   !!!Q: reset inbrc   
                           inbrr=0              !!!Q: reset inbrr
                        end if
                      end do
                      ct1=ct1+id(qqj)*id(qqi)
                end if
                !rite(*,*) 'assign to matX, matY', qqj, qqi
             end do
          end do

       !=================================
       ! v_N = TV_N/omega_N
       ! TV_N = omega_N - omega_N0
       ! omega_N0 = sum(E_mu + E_nu)(|X_mu,nu|^2 + |Y_mu, nu|^2)
       !==================================
       omega_N0 = 0
       omega_NN0 = 0
       omega_PN0 = 0
       do qqi = 1, iQRPA_LEN
          omega_NN0 = omega_NN0 + (REqpN(AB_row(qqi)) + REqpN(AB_col(qqi))) * &
                   (  real(CQRPA_XN_c(QQi))* real(CQRPA_XN_c(QQi)) + aimag(CQRPA_XN_c(QQi))* aimag(CQRPA_XN_c(QQi)) &
                    + real(CQRPA_YN_c(QQi))* real(CQRPA_YN_c(QQi)) + aimag(CQRPA_YN_c(QQi))* aimag(CQRPA_YN_c(QQi)) )
 
          omega_PN0 = omega_PN0 + (REqpP(AB_row(qqi)) + REqpP(AB_col(qqi))) * &
                   (  real(CQRPA_XP_c(QQi))* real(CQRPA_XP_c(QQi)) + aimag(CQRPA_XP_c(QQi))* aimag(CQRPA_XP_c(QQi)) &
                    + real(CQRPA_YP_c(QQi))* real(CQRPA_YP_c(QQi)) + aimag(CQRPA_YP_c(QQi))* aimag(CQRPA_YP_c(QQi)) )
 
       end do

       omega_N0 = omega_NN0 + omega_PN0
       v_N = (real_qrpa_omega - omega_N0)/real_qrpa_omega

      write(47,*) real_qrpa_omega, ',', ImTransStrN + ImTransStrP, ',', abs(v_N) 
 
    close(20)
    close(46)
    close(47)
  End subroutine
 
  !=================================================================================================================================

   subroutine qrpa_NH11(V,U,H,D,BH,BD)
    !
    ! Calculates H11_1 and H11_2 matrix elements.
    !
    Implicit None
    Complex(pr), Intent(in)     :: U(iHFB_LEN),V(iHFB_LEN),D(iQRPA_LEN),BD(iQRPA_LEN),H(iQRPA_LEN),BH(iQRPA_LEN)
    Complex(pr)                 :: NH11(iQRPA_LEN),BNH11(iQRPA_LEN)
    integer :: qqj


    Real :: imag_qrpa_omega, real_qrpa_omega
    Real(pr) :: Normfactor2               !!!Q: Loca0  
    Complex(pr)              :: iunit
    Parameter(iunit = Cmplx(0.0_pr,1.0_pr,kind=pr))




   !open(200,file=trim(adjustl(qq_io_path))//'QQ_binary_NH11_X', status='unknown',form='unformatted')
    open(200,file='QQ_binary_NH11_X', status='unknown',form='unformatted')
   !open(30,file=trim(adjustl(qq_io_path))//'NH11_print',access='append')
    open(30,file='NH11_print',access='append')


    imag_qrpa_omega = aimag(qrpa_omega)
    real_qrpa_omega = real(qrpa_omega)
    Normfactor2 = imag_qrpa_omega/( Pi*ImTransStrN + Pi*ImTransStrP )
    !
    ! V and U are complex
    call qrpa_TriProd('c','n','n',U,H,U,qrpa_AUX)
    call qrpa_TriProd('c','t','n',V,BH,V,qrpa_AUX20)
    call qrpa_TriProd('c','n','n',U,D,V,qrpa_BUX)
   !call qrpa_TriProd('c','t','n',V,BD,U,qrpa_BUX20)
    call qrpa_TriProd('c','c','n',V,BD,U,qrpa_BUX20)

   !H11 = qrpa_AUX -qrpa_AUX20 +qrpa_BUX -qrpa_BUX20
    NH11 = qrpa_AUX -qrpa_AUX20 +qrpa_BUX +qrpa_BUX20



    call qrpa_TriProd('t','n','*',U,BH,U,qrpa_AUX)
    call qrpa_TriProd('t','t','*',V,H,V,qrpa_AUX20)
    call qrpa_TriProd('t','t','*',U,D,V,qrpa_BUX)
   !call qrpa_TriProd('t','n','*',V,BD,U,qrpa_BUX20)
    call qrpa_TriProd('t','*','*',V,BD,U,qrpa_BUX20)

   !BH11 = qrpa_AUX -qrpa_AUX20 +qrpa_BUX -qrpa_BUX20
    BNH11 = qrpa_AUX -qrpa_AUX20 +qrpa_BUX +qrpa_BUX20
    !

    NH11 = iunit*sqrt(Normfactor2)*NH11
    BNH11 = iunit*sqrt(Normfactor2)*BNH11

    do qqj=1,size(NH11)    
       write(30,*) qqj,NH11(qqj), BNH11(qqj)
    end do

    write(200) NH11
    write(200) BNH11
 
    !
    close(30)
    close(200)
  end subroutine qrpa_NH11
 
  !
  ! Subroutine which reads the canonical conjugate operator Q
  ! not tested yet.
  !
  !=================================================================================================================================

  !=================================================================================================================================
  subroutine qrpa_PH11(V,U,H,D,BH,BD)
    !
    ! Calculates H11_1 and H11_2 matrix elements.
    !
    Implicit None
    Complex(pr), Intent(in)     :: U(iHFB_LEN),V(iHFB_LEN),D(iQRPA_LEN),BD(iQRPA_LEN),H(iQRPA_LEN),BH(iQRPA_LEN)
    Complex(pr)                 :: PH11(iQRPA_LEN),BPH11(iQRPA_LEN)
    integer :: qqj


    Real :: imag_qrpa_omega, real_qrpa_omega
    Real(pr) :: Normfactor2               !!!Q: Loca0  
    Complex(pr)              :: iunit
    Parameter(iunit = Cmplx(0.0_pr,1.0_pr,kind=pr))



   !open(201,file=trim(adjustl(qq_io_path))//'QQ_binary_PH11_X', status='unknown',form='unformatted')
    open(201,file='QQ_binary_PH11_X', status='unknown',form='unformatted')
   !open(31,file=trim(adjustl(qq_io_path))//'PH11_print',access='append')
    open(31,file='PH11_print',access='append')

    imag_qrpa_omega = aimag(qrpa_omega)
    real_qrpa_omega = real(qrpa_omega)
    Normfactor2 = imag_qrpa_omega/( Pi*ImTransStrN + Pi*ImTransStrP )
 
    !
    ! V and U are complex
    call qrpa_TriProd('c','n','n',U,H,U,qrpa_AUX)
    call qrpa_TriProd('c','t','n',V,BH,V,qrpa_AUX20)
    call qrpa_TriProd('c','n','n',U,D,V,qrpa_BUX)
   !call qrpa_TriProd('c','t','n',V,BD,U,qrpa_BUX20)
    call qrpa_TriProd('c','c','n',V,BD,U,qrpa_BUX20)

   !H11 = qrpa_AUX -qrpa_AUX20 +qrpa_BUX -qrpa_BUX20
    PH11 = qrpa_AUX -qrpa_AUX20 +qrpa_BUX +qrpa_BUX20



    call qrpa_TriProd('t','n','*',U,BH,U,qrpa_AUX)
    call qrpa_TriProd('t','t','*',V,H,V,qrpa_AUX20)
    call qrpa_TriProd('t','t','*',U,D,V,qrpa_BUX)
   !call qrpa_TriProd('t','n','*',V,BD,U,qrpa_BUX20)
    call qrpa_TriProd('t','*','*',V,BD,U,qrpa_BUX20)

   !BH11 = qrpa_AUX -qrpa_AUX20 +qrpa_BUX -qrpa_BUX20
    BPH11 = qrpa_AUX -qrpa_AUX20 +qrpa_BUX +qrpa_BUX20
    !

    PH11 = iunit*sqrt(Normfactor2)*PH11
    BPH11 = iunit*sqrt(Normfactor2)*BPH11


    do qqj=1,size(PH11)    
       write(31,*) qqj,PH11(qqj),BPH11(qqj)
    end do

    write(201) PH11
    write(201) BPH11
    !
 
    close(31)
    close(201)
  end subroutine qrpa_PH11
 




  !=================================================================================================================================
  !
  ! Subroutine which reads the canonical conjugate operator Q
  ! not tested yet.
  !
  !=================================================================================================================================
  Subroutine read_operator(QN,QP)
    Implicit None
    Complex(pr), Intent(out) :: QN(iQRPA_LEN),QP(iQRPA_LEN)
    Complex(pr)              :: omegaRead,XNread,YNread,XPread,YPread
    Real(pr)                 :: ImStN,ImStP,ReStN,ReStP
    Complex(pr)              :: iunit
    Integer(ipr)             :: i,Tread,Lread,Kread,Pread,OPRread,LenRead
    logical                  :: ENABLE_N_EFAin,ENABLE_P_EFAin
    Parameter(iunit = Cmplx(0.0_pr,1.0_pr,kind=pr))

    open(20,file="amplitudes.dat",form='unformatted',status='old',action='read')

    Read(20) Tread,Lread,Kread,Pread,OPRread,LenRead
    Read(20) ENABLE_N_EFAin,ENABLE_P_EFAin
    Read(20) omegaRead
    if(ENABLE_N_EFAin.or.ENABLE_P_EFAin) STOP "EFA NOT YET SUPPORTED WITH READ AMPLITUDES"
    if(Kread.ne.K_qrpa_responce) stop "Wrong K in the read operator!"
    if(Pread.ne.P_qrpa_response) stop "Wrong parity in the read operator!"
    if(LenRead.ne.iQRPA_LEN)     stop "Incorrect length in read operator!"
    if(abs(omegaRead).gt.spacing(1.0_pr)) stop "Solution on the file has omega > 0!"

    Read(20) ImStN,ImStP,ReStN,ReStP
    Do i = 1,iQRPA_LEN
      Read(20) XNread,YNread,XPread,YPread
      QN(i) = -0.5_pr*iunit*(XNread+conjg(YNread))/cmplx(ReStN+ReStP,ImStN+ImStP,kind=pr)
      QP(i) = -0.5_pr*iunit*(XPread+conjg(YPread))/cmplx(ReStN+ReStP,ImStN+ImStP,kind=pr)
    End Do

    close(20)
  End subroutine
  !=================================================================================================================================
  !  
  ! Helper subroutine to print qrpa densitiens to file
  !
  !=================================================================================================================================
  Subroutine print_density(densi,fnamepart)
   !Use HFBTHO, only : fl, fh
    Implicit None
    Complex(pr), Intent(in)       :: densi(nghl)
    Character(LEN=*), intent(in) :: fnamepart
    Integer(ipr)      :: ihli
    integer(ipr)      :: outformat
    Character(LEN=100) :: fnamer,fnamei
    !
    outformat = 1 !! Table
    fnamer = trim('qrpa_density_' // trim(fnamepart) // '_re.dat')
    fnamei = trim('qrpa_density_' // trim(fnamepart) // '_im.dat')

    Open (123,file=fnamer,status='replace')
    Write(123,*) "! Real part of density"

    if(outformat.eq.1) write(123,*) "! list of r_p, z, Re(density)"
    Do ihli=1,nghl
      if(outformat.eq.1) write(123,'(f11.6, " ", f11.6, " ", f18.12)') fl(ihli),fh(ihli),Real(densi(ihli))
    End do
    Close(123)

    Open (123,file=fnamei,status='replace')
    Write(123,*) "! Imaginary part of density"

    if(outformat.eq.1) write(123,*) "! list of r_p, z, Im(density)"
    Do ihli=1,nghl
      if(outformat.eq.1) write(123,'(f11.6, " ", f11.6, " ", f18.12)') fl(ihli),fh(ihli),Imag(densi(ihli))
    End do
    Close(123)
    !
  End Subroutine print_density
  !=================================================================================================================================
  !  
  ! Helper subroutine to print hfb densitiens to file
  !
  !=================================================================================================================================
  Subroutine print_density_hfb(densi,fnamepart)
   !Use HFBTHO, only : fl, fh
    Implicit None
    real(pr), Intent(in)          :: densi(nghl)
    Character(LEN=*), intent(in) :: fnamepart
    Integer(ipr)      :: ihli
    integer(ipr)      :: outformat
    Character(LEN=100) :: fname
    !
    outformat = 1 !! Table
    fname = trim('hfb_density_' // trim(fnamepart) // '_re.dat')

    Open (123,file=fname,status='replace')
    Write(123,*) "! Real part of density"

    if(outformat.eq.1) write(123,*) "! list of r_p, z, density"
    Do ihli=1,nghl
      if(outformat.eq.1) write(123,'(f11.6, " ", f11.6, " ", f18.12)') fl(ihli),fh(ihli),densi(ihli)
    End do
    Close(123)
  End Subroutine print_density_hfb
  !=================================================================================================================================
  !
  ! Helper routines for QRPA matrix handling
  !
  !=================================================================================================================================
  function qrpa_isHermitian(A,epsin)
    !-----------------------------------------------------------------------
    ! Returns true if A is Hermitian. Tolerance set by spacing, if not
    ! given. Used for debug. Note: If A=zero, returns true
    !-----------------------------------------------------------------------
    Implicit None
    Complex(pr), Intent(in)        :: A(:)
    Real(pr), Intent(in), optional :: epsin
    Complex(pr), Allocatable    :: temp(:)
    Real(pr) :: eps
    logical  :: qrpa_isHermitian
    integer  :: i
 
    eps = spacing(1.0_pr)
    if(present(epsin)) eps = epsin
    qrpa_isHermitian = .true.
    Allocate(temp(iQRPA_LEN))
    call qrpa_GiveTransposed(A,temp)
    temp = (A - conjg(temp))
    Do i=1,iQRPA_LEN
      if( abs(temp(i)).gt.eps ) then
         qrpa_isHermitian = .false.
         return
      end if
    End do
    Deallocate(temp)
  End function qrpa_isHermitian
  !=================================================================================================================================
  !
  !=================================================================================================================================
  function qrpa_isSymmetric(A,epsin)
    !-----------------------------------------------------------------------
    ! Returns true if A is symmetric. Tolerance set by spacing, if not
    ! given. Used for debug. Note: If A=zero, returns true
    !-----------------------------------------------------------------------
    Implicit None
    Complex(pr), Intent(in)        :: A(:)
    Real(pr), Intent(in), optional :: epsin
    Complex(pr), Allocatable    :: temp(:)
    Real(pr) :: eps
    logical  :: qrpa_isSymmetric
    integer  :: i
 
    eps = spacing(1.0_pr)
    if(present(epsin)) eps = epsin
    qrpa_isSymmetric = .true.
    Allocate(temp(iQRPA_LEN))
    call qrpa_GiveTransposed(A,temp)
    temp = (A - temp)
    Do i=1,iQRPA_LEN
      if( abs(temp(i)).gt.eps ) then
         qrpa_isSymmetric = .false.
         Deallocate(temp)
         return
      end if
    End do
    Deallocate(temp)
  End function qrpa_isSymmetric
  !=================================================================================================================================
  !
  !=================================================================================================================================
  function qrpa_isAntiSymmetric(A,epsin)
    !-----------------------------------------------------------------------
    ! Returns true if A is antisymmetric. Tolerance set by spacing, if not
    ! given. Used for debug. Note: If A=zero, returns true
    !-----------------------------------------------------------------------
    Implicit None
    Complex(pr), Intent(in)        :: A(:)
    Real(pr), Intent(in), optional :: epsin
    Complex(pr), Allocatable    :: temp(:)
    Real(pr) :: eps
    logical  :: qrpa_isAntiSymmetric
    integer  :: i
 
    eps = spacing(1.0_pr)
    if(present(epsin)) eps = epsin
    qrpa_isAntiSymmetric = .true.
    Allocate(temp(iQRPA_LEN))
    call qrpa_GiveTransposed(A,temp)
    temp = (A + temp)
    Do i=1,iQRPA_LEN
      if( abs(temp(i)).gt.eps ) then
         qrpa_isAntiSymmetric = .false.
         Deallocate(temp)
         return
      end if
    End do
    Deallocate(temp)
  End function qrpa_isAntiSymmetric
  !=================================================================================================================================
  !
  !=================================================================================================================================
  function qrpa_isReal(A,epsin)
    !-----------------------------------------------------------------------
    ! Returns true if A is real. Tolerance set by spacing, if not
    ! given. Used for debug. Note: If A=zero, returns true
    !-----------------------------------------------------------------------
    Implicit None
    Complex(pr), Intent(in)        :: A(:)
    Real(pr), Intent(in), optional :: epsin
    Complex(pr), Allocatable    :: temp(:)
    Real(pr) :: eps
    logical  :: qrpa_isReal
    integer  :: i
 
    eps = spacing(1.0_pr)
    if(present(epsin)) eps = epsin
    qrpa_isReal = .true.
    Allocate(temp(iQRPA_LEN))
    temp = Imag(A)
    Do i=1,iQRPA_LEN
      if( abs(temp(i)).gt.eps ) then
         qrpa_isReal = .false.
         Deallocate(temp)
         return
      end if
    End do
    Deallocate(temp)
  End function qrpa_isReal
  !=================================================================================================================================
  !
  !=================================================================================================================================
  function qrpa_isImag(A,epsin)
    !-----------------------------------------------------------------------
    ! Returns true if A is purely imaginary. Tolerance set by spacing, if not
    ! given. Used for debug. Note: If A=zero, returns true
    !-----------------------------------------------------------------------
    Implicit None
    Complex(pr), Intent(in)        :: A(:)
    Real(pr), Intent(in), optional :: epsin
    Complex(pr), Allocatable    :: temp(:)
    Real(pr) :: eps
    logical  :: qrpa_isImag
    integer  :: i
 
    eps = spacing(1.0_pr)
    if(present(epsin)) eps = epsin
    qrpa_isImag = .true.
    Allocate(temp(iQRPA_LEN))
    temp = Real(A)
    Do i=1,iQRPA_LEN
      if( abs(temp(i)).gt.eps ) then
         qrpa_isImag = .false.
         Deallocate(temp)
         return
      end if
    End do
    Deallocate(temp)
  End function qrpa_isImag
  !=================================================================================================================================
  !
  !=================================================================================================================================
  function qrpa_isZero(A,epsin)
    !-----------------------------------------------------------------------
    ! Returns true if A is real. Tolerance set by spacing, if not
    ! given. Used for debug
    !-----------------------------------------------------------------------
    Implicit None
    Complex(pr), Intent(in)        :: A(:)
    Real(pr), Intent(in), optional :: epsin
    Complex(pr), Allocatable    :: temp(:)
    Real(pr) :: eps
    logical  :: qrpa_isZero
    integer  :: i
 
    eps = spacing(1.0_pr)
    if(present(epsin)) eps = epsin
    qrpa_isZero = .true.
    Do i=1,iQRPA_LEN
      if( abs(A(i)).gt.eps ) then
        qrpa_isZero = .false.
        return
      end if
    End do
  End function qrpa_isZero
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_printProp(s,A)
    !
    ! prints qrpa matrix A properties, for debug. Put variable name, or such, for s
    !
    Implicit None
    Character(LEN=*), intent(in) :: s
    Complex(pr), Intent(in)      :: A(:)
    write(*,*) trim(s)," is Hermitian: ", qrpa_isHermitian(A,0.0001_pr), &
                       " symmetric: "   , qrpa_isSymmetric(A,0.0001_pr), &
                       " antisymmetric: ",qrpa_isAntiSymmetric(A,0.0001_pr), &
                       " real: "        , qrpa_isReal(A,0.0001_pr), &
                       " imaginary: "   , qrpa_isImag(A,0.0001_pr), &
                       " zero: "        , qrpa_isZero(A,0.0001_pr)
  end subroutine qrpa_printProp
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_SYM(A)
    !-----------------------------------------------------------------------
    ! symmetrize QRPA matrix so that A^T = A. This routine is only for debug
    !-----------------------------------------------------------------------
    Implicit None
    Complex(pr), Intent(inout)  :: A(:)
    Complex(pr), Allocatable    :: temp(:)

    Allocate(temp(iQRPA_LEN))
    call qrpa_GiveTransposed(A,temp)
    A = half * (A + temp)    
    Deallocate(temp)
  End Subroutine qrpa_SYM
  !=================================================================================================================================
  !
  !=================================================================================================================================
  Subroutine qrpa_ADJ(A)
    !-------------------------------------------------------------------------
    ! Hermitice QRPA matrix so that A^\dag = A. This routine is only for debug
    !-------------------------------------------------------------------------
    Implicit None
    Complex(pr), Intent(inout)  :: A(:)
    Complex(pr), Allocatable    :: temp(:)

    Allocate(temp(iQRPA_LEN))
    call qrpa_GiveTransposed(A,temp)
    A = half * (A + conjg(temp))
    Deallocate(temp)
  End Subroutine qrpa_ADJ
  !  


  !==============================================QQend

  !=======================================================================
  !
  !=======================================================================
End Module HFBTHO_solver

