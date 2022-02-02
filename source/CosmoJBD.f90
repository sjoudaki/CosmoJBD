! Likelihood code for joint analysis of measurements from overlapping spectroscopic and tomographic lensing surveys.
! Data vector includes 5 statistics: the 2pt shear correlation functions \xi_+ and \xi_-, 
! the galaxy-galaxy-lensing angular cross-correlation \gamma_t, 
! and redshift-space galaxy clustering spectra in the form of the monopole P_0 and quadrupole P_2. 
! Datasets: KiDS, CFHTLenS, 2dFLenS, BOSS, WiggleZ
! Original version for cosmic shear in 2015 (S. Joudaki), updated for 3x2pt in 2017 (S. Joudaki, A. Johnson),
! and updated for JBD gravity in 2020 (S. Joudaki, N. A. Lima).
!--------------------------------------------------------------------------------------------------------------------


    !Extra Module for dealing with window functions, i.e. importing and starting splines
    !note: dimensions of the interpolation grids are k (wavenumber), L (L indicates the multipole moment 0,2,4 - > 1,2,3)
    module multipoleinfo
    use settings
    use CosmologyTypes
    use Interpolation
    implicit none

    Type(TCubicSpline),  allocatable :: Win1D_NGC(:)
    Type(TCubicSpline),  allocatable :: Win1D_SGC(:)
    Type(TInterpGrid2D), allocatable :: Win2D_NGC(:,:)
    Type(TInterpGrid2D), allocatable :: Win2D_SGC(:,:)
    REAL(mcp) :: W02_norm_IC_NGC,W02_norm_IC_SGC,temp_k
    INTEGER :: NUM_k1D_win =  800,NUM_kp2D_win = 800,NUM_k2D_win = 80,num_multipoles = 3
    real(mcp) :: delta_kp = 0.0005, k_spacing = 0.005
    real(mcp) :: kp_ini = 0.00025, k_initial = 0.0025, piaj = 3.14159265359d0

    contains

    subroutine multipoleinfo_init(Datasets_Mp)
        type(TSettingIni) :: Datasets_Mp
        Type(TTextFile) :: F
        real(mcp), allocatable, dimension(:) :: k_values_1D,k_values_2D_k,kp_values2D
        real(mcp), allocatable, dimension(:) :: Win1D_L0_NGC,Win1D_L2_NGC,Win1D_L4_NGC
        real(mcp), allocatable, dimension(:) :: Win1D_L0_SGC,Win1D_L2_SGC,Win1D_L4_SGC
        real(mcp), dimension(:,:), pointer :: Win2D_NGC00,Win2D_NGC02,Win2D_NGC04
        real(mcp), dimension(:,:), pointer :: Win2D_SGC00,Win2D_SGC02,Win2D_SGC04
        real(mcp), dimension(:,:), pointer :: Win2D_NGC20,Win2D_NGC22,Win2D_NGC24
        real(mcp), dimension(:,:), pointer :: Win2D_SGC20,Win2D_SGC22,Win2D_SGC24
        real(mcp), dimension(:,:), pointer :: Win2D_NGC40,Win2D_NGC42,Win2D_NGC44
        real(mcp), dimension(:,:), pointer :: Win2D_SGC40,Win2D_SGC42,Win2D_SGC44
        character(LEN=:), allocatable :: window1D_cmass_NGC,window1D_cmass_SGC
        character(LEN=:), allocatable :: window2D_cmass_NGC,window2D_cmass_SGC
        integer   :: ik,inum,inum2,it,jt
        real(mcp) :: kgrid,dummy

        allocate(Win1D_NGC(num_multipoles))
        allocate(Win1D_SGC(num_multipoles))
        allocate(Win2D_NGC(num_multipoles,num_multipoles))
        allocate(Win2D_SGC(num_multipoles,num_multipoles))

        allocate(k_values_1D(NUM_k1D_win))
        allocate(Win1D_L0_NGC(NUM_k1D_win))
        allocate(Win1D_L2_NGC(NUM_k1D_win))
        allocate(Win1D_L4_NGC(NUM_k1D_win))
        allocate(Win1D_L0_SGC(NUM_k1D_win))
        allocate(Win1D_L2_SGC(NUM_k1D_win))
        allocate(Win1D_L4_SGC(NUM_k1D_win))

        allocate(k_values_2D_k(NUM_k2D_win))
        allocate(kp_values2D(NUM_kp2D_win))

        allocate(Win2D_NGC00(NUM_k2D_win,NUM_kp2D_win))
        allocate(Win2D_NGC02(NUM_k2D_win,NUM_kp2D_win))
        allocate(Win2D_NGC04(NUM_k2D_win,NUM_kp2D_win))
        allocate(Win2D_SGC00(NUM_k2D_win,NUM_kp2D_win))
        allocate(Win2D_SGC02(NUM_k2D_win,NUM_kp2D_win))
        allocate(Win2D_SGC04(NUM_k2D_win,NUM_kp2D_win))

        allocate(Win2D_NGC20(NUM_k2D_win,NUM_kp2D_win))
        allocate(Win2D_NGC22(NUM_k2D_win,NUM_kp2D_win))
        allocate(Win2D_NGC24(NUM_k2D_win,NUM_kp2D_win))
        allocate(Win2D_SGC20(NUM_k2D_win,NUM_kp2D_win))
        allocate(Win2D_SGC22(NUM_k2D_win,NUM_kp2D_win))
        allocate(Win2D_SGC24(NUM_k2D_win,NUM_kp2D_win))

        allocate(Win2D_NGC40(NUM_k2D_win,NUM_kp2D_win))
        allocate(Win2D_NGC42(NUM_k2D_win,NUM_kp2D_win))
        allocate(Win2D_NGC44(NUM_k2D_win,NUM_kp2D_win))
        allocate(Win2D_SGC40(NUM_k2D_win,NUM_kp2D_win))
        allocate(Win2D_SGC42(NUM_k2D_win,NUM_kp2D_win))
        allocate(Win2D_SGC44(NUM_k2D_win,NUM_kp2D_win))

        !-----------------------------------------------------------------------------------------------
    ! Import 1D and 2D window function values for both NGC and SGC
    !-----------------------------------------------------------------------------------------------

        window1D_cmass_NGC = Datasets_Mp%ReadFileName('window1D_NGC')
        window1D_cmass_SGC = Datasets_Mp%ReadFileName('window1D_SGC')
        window2D_cmass_NGC = Datasets_Mp%ReadFileName('window2D_NGC')
        window2D_cmass_SGC = Datasets_Mp%ReadFileName('window2D_SGC')

        call F%Open(window1D_cmass_NGC)
        do ik = 1,NUM_k1D_win
                read (F%unit,*) k_values_1D(ik), Win1D_L0_NGC(ik),Win1D_L2_NGC(ik),Win1D_L4_NGC(ik),W02_norm_IC_NGC
        end do

        call F%Open(window1D_cmass_SGC)
        do ik = 1,NUM_k1D_win
                read (F%unit,*) k_values_1D(ik), Win1D_L0_SGC(ik),Win1D_L2_SGC(ik),Win1D_L4_SGC(ik),W02_norm_IC_SGC
        end do

        call F%Open(window2D_cmass_NGC)
        do it = 1,NUM_k2D_win
                do jt = 1,NUM_kp2D_win
                        read (F%unit,*) dummy,dummy, Win2D_NGC00(it,jt), Win2D_NGC02(it,jt), &
                        Win2D_NGC04(it,jt),Win2D_NGC20(it,jt), Win2D_NGC22(it,jt), &
                        Win2D_NGC24(it,jt),Win2D_NGC40(it,jt), Win2D_NGC42(it,jt),&
                        Win2D_NGC44(it,jt)
                end do
        end do

        call F%Open(window2D_cmass_SGC)
        do it = 1,NUM_k2D_win
                do jt = 1,NUM_kp2D_win
                        read (F%unit,*) dummy,dummy, Win2D_SGC00(it,jt), &
                        Win2D_SGC02(it,jt),Win2D_SGC04(it,jt),Win2D_SGC20(it,jt), &
                        Win2D_SGC22(it,jt),Win2D_SGC24(it,jt),Win2D_SGC40(it,jt), &
                        Win2D_SGC42(it,jt),Win2D_SGC44(it,jt)
                end do
        end do

        !-----------------------------------------------------------------------------------------------
    ! Setup splines for all window functions (1D + 2D) and for SGC and NGC
    !-----------------------------------------------------------------------------------------------

        do inum = 1, NUM_k2D_win
                k_values_2D_k(inum) = k_initial + k_spacing*(inum - 1.0)
        end do
	do inum =1, NUM_kp2D_win
                kp_values2D(inum) = kp_ini + delta_kp*(inum - 1.0)
        end do

        call Win1D_NGC(1)%Init(k_values_1D,Win1D_L0_NGC,n=NUM_k1D_win)
        call Win1D_NGC(2)%Init(k_values_1D,Win1D_L2_NGC,n=NUM_k1D_win)
        call Win1D_NGC(3)%Init(k_values_1D,Win1D_L4_NGC,n=NUM_k1D_win)

        call Win1D_SGC(1)%Init(k_values_1D,Win1D_L0_SGC,n=NUM_k1D_win)
        call Win1D_SGC(2)%Init(k_values_1D,Win1D_L2_SGC,n=NUM_k1D_win)
        call Win1D_SGC(3)%Init(k_values_1D,Win1D_L4_SGC,n=NUM_k1D_win)

        call Win2D_NGC(1,1)%Init(k_values_2D_k,kp_values2D,Win2D_NGC00)
        call Win2D_NGC(1,2)%Init(k_values_2D_k,kp_values2D,Win2D_NGC02)
        call Win2D_NGC(1,3)%Init(k_values_2D_k,kp_values2D,Win2D_NGC04)

        call Win2D_NGC(2,1)%Init(k_values_2D_k,kp_values2D,Win2D_NGC20)
        call Win2D_NGC(2,2)%Init(k_values_2D_k,kp_values2D,Win2D_NGC22)
        call Win2D_NGC(2,3)%Init(k_values_2D_k,kp_values2D,Win2D_NGC24)

        call Win2D_NGC(3,1)%Init(k_values_2D_k,kp_values2D,Win2D_NGC40)
        call Win2D_NGC(3,2)%Init(k_values_2D_k,kp_values2D,Win2D_NGC42)
        call Win2D_NGC(3,3)%Init(k_values_2D_k,kp_values2D,Win2D_NGC44)

        call Win2D_SGC(1,1)%Init(k_values_2D_k,kp_values2D,Win2D_SGC00)
        call Win2D_SGC(1,2)%Init(k_values_2D_k,kp_values2D,Win2D_SGC02)
        call Win2D_SGC(1,3)%Init(k_values_2D_k,kp_values2D,Win2D_SGC04)

        call Win2D_SGC(2,1)%Init(k_values_2D_k,kp_values2D,Win2D_SGC20)
        call Win2D_SGC(2,2)%Init(k_values_2D_k,kp_values2D,Win2D_SGC22)
        call Win2D_SGC(2,3)%Init(k_values_2D_k,kp_values2D,Win2D_SGC24)

        call Win2D_SGC(3,1)%Init(k_values_2D_k,kp_values2D,Win2D_SGC40)
        call Win2D_SGC(3,2)%Init(k_values_2D_k,kp_values2D,Win2D_SGC42)
        call Win2D_SGC(3,3)%Init(k_values_2D_k,kp_values2D,Win2D_SGC44)

        DEALLOCATE(k_values_1D,Win1D_L0_NGC,Win1D_L2_NGC,Win1D_L4_NGC)
        DEALLOCATE(Win1D_L0_SGC,Win1D_L2_SGC,Win1D_L4_SGC)
        DEALLOCATE(k_values_2D_k,kp_values2D)

        DEALLOCATE(Win2D_NGC00,Win2D_NGC02,Win2D_NGC04)
        DEALLOCATE(Win2D_SGC00,Win2D_SGC02,Win2D_SGC04)
        DEALLOCATE(Win2D_NGC20,Win2D_NGC22,Win2D_NGC24)
        DEALLOCATE(Win2D_SGC20,Win2D_SGC22,Win2D_SGC24)
        DEALLOCATE(Win2D_NGC40,Win2D_NGC42,Win2D_NGC44)
        DEALLOCATE(Win2D_SGC40,Win2D_SGC42,Win2D_SGC44)

    end subroutine multipoleinfo_init

    end module multipoleinfo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    module CosmoJBD 
    use CosmologyTypes
    use CAMB, only : ComovingRadialDistance, AngularDiameterDistance, AngularDiameterDistance2, f_K, Hofz  !distance also in Mpc no h units
    use constants
    use Precision
    use likelihood
    use settings
    use Interpolation
    use MatrixUtils
    use omp_lib
    use Likelihood_Cosmology
    use Calculator_Cosmology
    use CosmoTheory
    use ModelParams
    use multipoleinfo
    implicit none  
!    private

    TYPE WiggleZ_Multipole_data
        real(mcp), allocatable, dimension(:,:)    :: P0_data_lowZ,P2_data_lowZ,P4_data_lowZ
        real(mcp), allocatable, dimension(:,:)    :: P0_data_highZ,P2_data_highZ,P4_data_highZ
        real(mcp), allocatable, dimension(:,:,:)  :: INV_Cov_matrix_highZ,INV_Cov_matrix_lowZ
        real(mcp), allocatable, dimension(:,:,:)  :: Covolution_matrix_highZ,Covolution_matrix_lowZ
        character(LEN=:), allocatable :: dataset_dir
        real(mcp) :: z_eff_low,z_eff_high,dk_theory,k_min_theory
        INTEGER :: k_num_obs,size_cov,size_window,size_convolution
        INTEGER :: Num_regions,k_num_conv,total_num_k_bins
        LOGICAL :: fit_k_030,fit_k_019,fit_k_015,fit_k_017,use_AP
    end type WiggleZ_Multipole_data

    TYPE CMASS_Multipole_data
        real(mcp), allocatable, dimension(:) :: P0_multipole_NGC,P2_multipole_NGC,P4_multipole_NGC
        real(mcp), allocatable, dimension(:) :: P0_multipole_SGC,P2_multipole_SGC,P4_multipole_SGC
        real(mcp) :: z_eff, k_spacing,k_min,k_max,kmin_theory,kmax_theory,dk_theory,k_min_conv_matrix,dk_conv_matrix
        integer :: k_num, size_cov,size_cov_015,numk_values_conv,numk_values_full_conv
        character(LEN=:), allocatable :: dataset_dir
        logical :: fit_k20_cmass,fit_k15_cmass,fit_k10_cmass,Use_Conv_Matrix,use_full_conv
    end type CMASS_Multipole_data

    TYPE CMASS_Multipole_overlap
        real(mcp), allocatable, dimension(:)    :: P0_data,P2_data
        real(mcp), allocatable, dimension(:,:)  :: Covolution_matrix
        character(LEN=:), allocatable :: dataset_dir
        real(mcp) :: z_eff,dk_theory,k_min_theory
        INTEGER :: k_num_obs,size_cov,size_window,size_convolution
        INTEGER :: Num_regions,k_num_conv,total_num_k_bins
        real(mcp) :: k_spacing_obs,k_min_obs
        LOGICAL :: fit_k_0175,fit_k_0125,fit_k_0075,use_AP
    end type CMASS_Multipole_overlap

    TYPE LOWZ_Multipole_overlap
        real(mcp), allocatable, dimension(:)    :: P0_data,P2_data
        real(mcp), allocatable, dimension(:,:)  :: Covolution_matrix
        character(LEN=:), allocatable :: dataset_dir
        real(mcp) :: z_eff,dk_theory,k_min_theory
        INTEGER :: k_num_obs,size_cov,size_window,size_convolution
        INTEGER :: Num_regions,k_num_conv,total_num_k_bins
        real(mcp) :: k_spacing_obs,k_min_obs
        LOGICAL :: fit_k_0175,fit_k_0125,fit_k_0075,use_AP
    end type LOWZ_Multipole_overlap

    TYPE twodfloz_Multipole_overlap
        real(mcp), allocatable, dimension(:)    :: P0_data,P2_data
        real(mcp), allocatable, dimension(:,:)  :: Covolution_matrix
        character(LEN=:), allocatable :: dataset_dir
        real(mcp) :: z_eff,dk_theory,k_min_theory
        INTEGER :: k_num_obs,size_cov,size_window,size_convolution
        INTEGER :: Num_regions,k_num_conv,total_num_k_bins
        real(mcp) :: k_spacing_obs,k_min_obs
        LOGICAL :: fit_k_0175,fit_k_0125,fit_k_0075,use_AP
    end type twodfloz_Multipole_overlap

    TYPE twodfhiz_Multipole_overlap
        real(mcp), allocatable, dimension(:)    :: P0_data,P2_data
        real(mcp), allocatable, dimension(:,:)  :: Covolution_matrix
        character(LEN=:), allocatable :: dataset_dir
        real(mcp) :: z_eff,dk_theory,k_min_theory
        INTEGER :: k_num_obs,size_cov,size_window,size_convolution
        INTEGER :: Num_regions,k_num_conv,total_num_k_bins
        real(mcp) :: k_spacing_obs,k_min_obs
        LOGICAL :: fit_k_0175,fit_k_0125,fit_k_0075,use_AP
    end type twodfhiz_Multipole_overlap


    TYPE, extends(TCosmoCalcLikelihood) :: CosmoJBDLikelihood  
        real(mcp), dimension(2,70) :: arraysjorig1,arraysjorig2,arraysjorig3,arraysjorig4
        real(mcp), dimension(2,70,4) :: arraysjfull
        real(mcp), dimension(2,29) :: arraysjlenslowz,arraysjlens2dfloz
        real(mcp), dimension(2,28) :: arraysjlenscmass,arraysjlens2dfhiz
        real(mcp), allocatable, dimension(:) :: xipm !7 tom-bins, 7 ang-bins, 2 for +/- gives 28*7*2 = 392
        real(mcp), allocatable, dimension(:,:) :: covxipm, invcovxipm,covxipminv !nzbins*(nzbins+1)*nangbins
        real(mcp), dimension(9) :: thetacfhtini, thetaradcfhtini !nangbins
        real(mcp), dimension(58999) :: ellgentestarrini
        real(mcp), allocatable,dimension(:) :: maskelements
        real(mcp), dimension(9,58999) :: bes0arr,bes4arr,bes2arr
        integer :: size_cov,size_covmask,size_covmaskplanck,sizcov,sizcovpremask,klinesum,set_scenario,size_covgtmask,size_covallmask
        logical :: use_morell, use_rombint
        type (CMASS_Multipole_data) cmass_mp
        type (WiggleZ_Multipole_data) wigglez_mp
        type (CMASS_Multipole_overlap) cmass_mp_overlap
        type (LOWZ_Multipole_overlap) lowz_mp_overlap
        type (twodfloz_Multipole_overlap) twodfloz_mp_overlap
        type (twodfhiz_Multipole_overlap) twodfhiz_mp_overlap
        logical :: use_wigglez_full, use_cmass_full, use_cmass_overlap, use_lowz_overlap, use_2dfloz_overlap, use_2dfhiz_overlap

    contains

    procedure :: LogLike => CosmoJBD_LnLike
    procedure :: ReadIni => CosmoJBD_ReadIni

    END TYPE CosmoJBDLikelihood

    real(mcp), allocatable, dimension(:,:) ::  INV_Cov_matrix_NGC,INV_Cov_matrix_SGC

    logical :: use_CosmoJBD = .false.
!    public CosmoJBDLikelihood, CosmoJBDLikelihood_Add

    contains


    !--------------------------------------------------------------------------
    !Called in datalikelihoods.f90 to add likelihood
    !---------------------------------------------------------------------------
    subroutine CosmoJBDLikelihood_Add(LikeList, Ini)
    use settings
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: Ini
    Type(CosmoJBDLikelihood), pointer :: this
    Type(TSettingIni) :: DataSets
    Integer :: i= 1
    
	if(.not. Ini%Read_Logical('use_CosmoJBD',.false.)) return
    call Ini%TagValuesForName('CosmoJBD_dataset', DataSets)
    allocate(this)
    this%LikelihoodType = 'CosmoJBD'
    this%needs_nonlinear_pk = .true.
    this%needs_background_functions = .true.
    this%needs_powerspectra = .true.
    this%needs_exact_z = .true.
    this%speed = -1
    this%num_z = 37
    this%size_cov = 392 !nzbins*(nzbins+1)*nangbins

    CosmoSettings%use_growth_kz = .true.
    this%needs_weylpower = .true.

    this%use_wigglez_full = Ini%Read_Logical('use_WiggleZ_full',.false.)
    this%use_cmass_full = Ini%Read_Logical('use_CMASS_full',.false.)
    this%use_cmass_overlap = Ini%Read_Logical('use_CMASS_overlap',.false.)
    this%use_lowz_overlap = Ini%Read_Logical('use_LOWZ_overlap',.false.)
    this%use_2dfloz_overlap = Ini%Read_Logical('use_2dfloz_overlap',.false.)
    this%use_2dfhiz_overlap = Ini%Read_Logical('use_2dfhiz_overlap',.false.)

    if(this%use_cmass_full .AND. this%use_cmass_overlap) stop 'Change settings .ini file'
	
    this%cmass_mp%dataset_dir = DataSets%Value(i)
    this%cmass_mp%fit_k10_cmass = Ini%Read_Logical('fit_k10_cmass',.false.)
    this%cmass_mp%fit_k15_cmass = Ini%Read_Logical('fit_k15_cmass',.false.)
    this%cmass_mp%fit_k20_cmass = Ini%Read_Logical('fit_k20_cmass',.false.)
    this%cmass_mp%Use_Conv_Matrix = .false.
    this%cmass_mp%use_full_conv = .true.
	
    this%wigglez_mp%dataset_dir = DataSets%Value(i)
    this%wigglez_mp%fit_k_030 = Ini%Read_Logical('fit_k_030_wigZ',.false.)
    this%wigglez_mp%fit_k_019 = Ini%Read_Logical('fit_k_019_wigZ',.false.)
    this%wigglez_mp%fit_k_017 = Ini%Read_Logical('fit_k_017_wigZ',.false.)
    this%wigglez_mp%fit_k_015 = Ini%Read_Logical('fit_k_015_wigZ',.false.)
	
    this%cmass_mp_overlap%dataset_dir = DataSets%Value(i)
    this%cmass_mp_overlap%fit_k_0175 = Ini%Read_Logical('fit_k_0175_cmass_overlap',.false.)
    this%cmass_mp_overlap%fit_k_0125 = Ini%Read_Logical('fit_k_0125_cmass_overlap',.false.)
    this%cmass_mp_overlap%fit_k_0075 = Ini%Read_Logical('fit_k_0075_cmass_overlap',.false.)
	
    this%lowz_mp_overlap%dataset_dir = DataSets%Value(i)
    this%lowz_mp_overlap%fit_k_0175 = Ini%Read_Logical('fit_k_0175_lowz_overlap',.false.)
    this%lowz_mp_overlap%fit_k_0125 = Ini%Read_Logical('fit_k_0125_lowz_overlap',.false.)
    this%lowz_mp_overlap%fit_k_0075 = Ini%Read_Logical('fit_k_0075_lowz_overlap',.false.)

    this%twodfloz_mp_overlap%dataset_dir = DataSets%Value(i)
    this%twodfloz_mp_overlap%fit_k_0175 = Ini%Read_Logical('fit_k_0175_2dfloz_overlap',.false.)
    this%twodfloz_mp_overlap%fit_k_0125 = Ini%Read_Logical('fit_k_0125_2dfloz_overlap',.false.)
    this%twodfloz_mp_overlap%fit_k_0075 = Ini%Read_Logical('fit_k_0075_2dfloz_overlap',.false.)

    this%twodfhiz_mp_overlap%dataset_dir = DataSets%Value(i)
    this%twodfhiz_mp_overlap%fit_k_0175 = Ini%Read_Logical('fit_k_0175_2dfhiz_overlap',.false.)
    this%twodfhiz_mp_overlap%fit_k_0125 = Ini%Read_Logical('fit_k_0125_2dfhiz_overlap',.false.)
    this%twodfhiz_mp_overlap%fit_k_0075 = Ini%Read_Logical('fit_k_0075_2dfhiz_overlap',.false.)

    call this%ReadDatasetFile(DataSets%Value(1))
!    this%LikelihoodType = 'CosmoJBD'
    this%tag = DataSets%Name(1)
    call this%loadParamNames(trim(DataDir)//'CosmoJBD.paramnames')
    call LikeList%Add(this)
    if (Feedback>1) write(*,*) 'Imported CosmoJBD and multipole data'

    end subroutine CosmoJBDLikelihood_Add



    !-----------------------------------------------------------------
    !Read in data
    !-----------------------------------------------------------------
    subroutine CosmoJBD_ReadIni(this, Ini)
    use MatrixUtils
    use settings
    use omp_lib
    external DGETRF, DGETRI 
    class(CosmoJBDLikelihood) this
    class(TSettingIni) :: Ini
    integer, allocatable, dimension(:) :: ipiv
    real(mcp), allocatable, dimension(:) :: work_inverse
    real(mcp), allocatable, dimension(:,:) :: xipmtemp,covxipmtemp,covxipmtempinv,masktemp
    integer :: info,lwork,zinit,hay,iii,jjj,it,jt,kt,setscenario,sizcovish,sizcovishsq,sizcovishpremask

    this%set_scenario = Ini%Read_Int('set_scenario')
    this%use_morell = Ini%Read_Logical('use_morell',.false.)
    this%use_rombint = Ini%Read_Logical('use_rombint',.false.)

    setscenario = this%set_scenario

    this%size_covmask = 130
    this%size_covmaskplanck = 56
    this%size_covgtmask = 210
    this%size_covallmask = 170 ! 210 + 6 x 4 (lenses)

    if(setscenario == 0) then !no masking
        sizcovish = this%size_cov !180 elements post-masking
        sizcovishpremask = 180 !elements pre-masking, same number as post-masking, as this is the non-masking xi_+/- scenario
    end if
    if(setscenario == 1) then !fiducial masking advocated in Joudaki et al (2016)
        sizcovish = this%size_covmask !130 elements post-masking
        sizcovishpremask = 180 !elements pre-masking
    end if
    if(setscenario == 2) then
        sizcovish = this%size_covgtmask !210 elements post-masking
        sizcovishpremask = 324 !elements pre-masking, 2*9*4*(4+1)/2 + 9*4*4 (9 angular bins, 4 tomographic bins, 4 lenses)
    end if
    if(setscenario == 3) then
        sizcovish = this%size_covallmask !234 elements post-masking
        sizcovishpremask = 340 !elements pre-masking
    end if

    sizcovishsq = sizcovish**2
    this%sizcov = sizcovish
    this%sizcovpremask = sizcovishpremask

    allocate(this%xipm(sizcovish))
    allocate(this%covxipm(sizcovish,sizcovish))
    allocate(this%invcovxipm(sizcovish,sizcovish))
    allocate(this%covxipminv(sizcovish,sizcovish))
    allocate(this%maskelements(sizcovishpremask))
    allocate(masktemp(2,sizcovishpremask))
    allocate(xipmtemp(4,sizcovish))
    allocate(covxipmtemp(3,sizcovishsq)) !!32400 = 180^2

    !!!!!Lens redshifts!!!!!!!
    open(7,file='data/lensingrsdfiles/nz_lowz_modelsj.dat') !USED
    READ (7,*) this%arraysjlenslowz
    close(7)
    open(7,file='data/lensingrsdfiles/nz_cmass_modelsj.dat') !USED
    READ (7,*) this%arraysjlenscmass
    close(7)
    open(7,file='data/lensingrsdfiles/nz_2dfloz_modelsj.dat') !USED
    READ (7,*) this%arraysjlens2dfloz
    close(7)
    open(7,file='data/lensingrsdfiles/nz_2dfhiz_modelsj.dat') !USED
    READ (7,*) this%arraysjlens2dfhiz
    close(7)

    !!!Reading in source distributions
    !open(7,file='data/lensingrsdfiles/nz_z1_kids_binned.dat') !USED
    open(7,file='data/lensingrsdfiles/nz_z1_kids_binned_bootstrap.dat') !USED
    READ (7,*) this%arraysjorig1
    close(7)
    !open(7,file='data/lensingrsdfiles/nz_z2_kids_binned.dat') !USED
    open(7,file='data/lensingrsdfiles/nz_z2_kids_binned_bootstrap.dat') !USED
    READ (7,*) this%arraysjorig2
    close(7)
    !open(7,file='data/lensingrsdfiles/nz_z3_kids_binned.dat') !USED
    open(7,file='data/lensingrsdfiles/nz_z3_kids_binned_bootstrap.dat') !USED
    READ (7,*) this%arraysjorig3
    close(7)
    !open(7,file='data/lensingrsdfiles/nz_z4_kids_binned.dat') !USED
    open(7,file='data/lensingrsdfiles/nz_z4_kids_binned_bootstrap.dat') !USED
    READ (7,*) this%arraysjorig4
    close(7)

    this%arraysjorig1(2,:) = this%arraysjorig1(2,:)/(sum(this%arraysjorig1(2,:))*0.05d0)
    this%arraysjorig2(2,:) = this%arraysjorig2(2,:)/(sum(this%arraysjorig2(2,:))*0.05d0)
    this%arraysjorig3(2,:) = this%arraysjorig3(2,:)/(sum(this%arraysjorig3(2,:))*0.05d0)
    this%arraysjorig4(2,:) = this%arraysjorig4(2,:)/(sum(this%arraysjorig4(2,:))*0.05d0)
    this%arraysjfull(:,:,1) = this%arraysjorig1
    this%arraysjfull(:,:,2) = this%arraysjorig2
    this%arraysjfull(:,:,3) = this%arraysjorig3
    this%arraysjfull(:,:,4) = this%arraysjorig4

    this%arraysjlenslowz(2,:) = this%arraysjlenslowz(2,:)/(sum(this%arraysjlenslowz(2,:))*0.01d0)
    this%arraysjlenscmass(2,:) = this%arraysjlenscmass(2,:)/(sum(this%arraysjlenscmass(2,:))*0.01d0)
    this%arraysjlens2dfloz(2,:) = this%arraysjlens2dfloz(2,:)/(sum(this%arraysjlens2dfloz(2,:))*0.01d0)
    this%arraysjlens2dfhiz(2,:) = this%arraysjlens2dfhiz(2,:)/(sum(this%arraysjlens2dfhiz(2,:))*0.01d0)

    !!!Reading in measurements
    if(setscenario == 0) then
        open(7,file='data/lensingrsdfiles/xipm_kids4tom_regcomb_blind1_edinburghsj.dat') 
    end if
    if(setscenario == 1) then
        open(7,file='data/lensingrsdfiles/xipmcut_kids_regcomb_blind2_swinburnesj.dat')
    end if
    if(setscenario == 2) then
        open(7,file='data/lensingrsdfiles/xipmgtcut_kids_regcomb_blind1sj.dat') !USED
    end if
    if(setscenario == 3) then
        open(7,file='data/lensingrsdfiles/pythonmeasurements170.dat') !USED
    end if
    READ (7,*) xipmtemp
    close(7)
    this%xipm = xipmtemp(2,:)


    !!!Masking file
    if(setscenario == 0) then
        open(7,file='data/lensingrsdfiles/xipm_kids4tom_selectsjnocut.dat') !USED
    end if
    if(setscenario == 1) then
        open(7,file='data/lensingrsdfiles/xipm_kids4tom_selectsj.dat') !USED
    end if
    if(setscenario == 2) then
        open(7,file='data/lensingrsdfiles/xipmgt_kids4tom_selectsj.dat') !USED
    end if
    if(setscenario == 3) then
        open(7,file='data/lensingrsdfiles/xipmgtpllarge4_kids4tom_selectsj_withoutboss.dat') !USED
    end if
    READ (7,*) masktemp
    close(7)
    this%maskelements = masktemp(2,:)

    if(setscenario == 1) then
        open(7,file='data/lensingrsdfiles/xipmcutcov_kids_regcomb_blind2_swinburnesj.dat') !USED
    end if
    if(setscenario == 2) then
        open(7,file='data/lensingrsdfiles/xipmgtcutcov_kids_regcomb_blind1sj.dat') !USED
    end if
    if(setscenario == 3) then
        open(7,file='data/lensingrsdfiles/pythoncov170.dat') !USED
    end if
    READ (7,*) covxipmtemp
    close(7)

    kt = 0
    do it=1,sizcovish
        do jt=1,sizcovish
            kt = kt+1
            this%covxipm(it,jt) = covxipmtemp(3,kt)
        end do
    end do

    !converted from arcmin to degrees
    this%thetacfhtini = (/ 0.71336d0, 1.45210d0, 2.95582d0, 6.01675d0, 12.24745d0, 24.93039d0, 50.74726d0, 103.29898d0, 210.27107d0 /)/60.0d0 !athena angular scales (deg) --- 9 ang bins --- slightly diff from athena
    this%thetaradcfhtini = this%thetacfhtini*3.14159265359d0/180.0d0 !converted to rad
    this%ellgentestarrini = (/ (hay,hay=2,59000,1) /)

    !compute Bessel functions
    do iii=1,9
        do jjj=1,58999
            this%bes0arr(iii,jjj) = bessel_j0(this%ellgentestarrini(jjj)*this%thetaradcfhtini(iii))
            this%bes4arr(iii,jjj) = bessel_jn(4,this%ellgentestarrini(jjj)*this%thetaradcfhtini(iii))
            this%bes2arr(iii,jjj) = bessel_jn(2,this%ellgentestarrini(jjj)*this%thetaradcfhtini(iii))
        end do 
    end do

    allocate(this%exact_z(this%num_z)) !allocate array for matter power spectrum redshifts
    !assign redshifts for matter power spectrum
    this%exact_z(1) = 0.0d0
    this%exact_z(2) = 0.025d0
    do zinit=3,37
        this%exact_z(zinit) = this%exact_z(zinit-1) + 0.154d0*this%exact_z(zinit-1)
    end do
    this%exact_z(37) = 3.474999d0


    if(this%use_cmass_full)    call Read_In_CMASS_data(this,Ini)
    if(this%use_wigglez_full)  call Read_In_Wigglez_data(this, Ini)
    if(this%use_cmass_overlap) call Read_In_CMASS_data_Overlap(this, Ini)
    if(this%use_lowz_overlap)  call Read_In_LOWZ_data_Overlap(this, Ini)
    if(this%use_2dfloz_overlap)  call Read_In_2dfloz_data_Overlap(this, Ini)
    if(this%use_2dfhiz_overlap)  call Read_In_2dfhiz_data_Overlap(this, Ini)

    this%cmass_mp%z_eff = 0.57
    this%wigglez_mp%z_eff_low = 0.44
    this%wigglez_mp%z_eff_high = 0.73
    this%lowz_mp_overlap%z_eff = 0.32
    this%twodfloz_mp_overlap%z_eff = 0.31
    this%twodfhiz_mp_overlap%z_eff = 0.56

    !calculate inverse of covariance
    allocate(work_inverse(sizcovish*sizcovish))         !for inverse calculation
    allocate(ipiv(sizcovish))                           !for inverse calculation
    this%invcovxipm = this%covxipm
    LWORK = sizcovish**2
    call DGETRF(sizcovish,sizcovish,this%invcovxipm,sizcovish,ipiv,info)         !LU decomposition
    call DGETRI(sizcovish,this%invcovxipm,sizcovish,ipiv,work_inverse,lwork,info)    !inverse from LU decompostion
    ipiv(:) =  0
    work_inverse(:) = 0    
    if (info .ne. 0) stop 'Problem with the matrix inverse calculation'

    if(setscenario == 0) then
        this%invcovxipm = (this%invcovxipm)*(930.0-180.0-2.0)/(930.0-1.0) !Hartlap correction factor
    end if
    if(setscenario == 1) then
        this%invcovxipm = this%invcovxipm !Sellentin-Heavens or Analytic covariance instead
    end if
    if(setscenario == 2) then
        this%invcovxipm = (this%invcovxipm)*(930.0-210.0-2.0)/(930.0-1.0) !Hartlap correction factor
    end if
    if(setscenario == 3) then
        this%invcovxipm = this%invcovxipm !use Sellentin-Heavens instead
    end if

    DEALLOCATE(work_inverse) 
    DEALLOCATE(ipiv) 
    DEALLOCATE(xipmtemp) 
    DEALLOCATE(covxipmtemp) 

    this%klinesum = 0

    end subroutine CosmoJBD_ReadIni


    !-------------------------------------------------------------------------------------------
    ! Read in CMASS Data
    !-------------------------------------------------------------------------------------------
    subroutine Read_In_CMASS_data(This,Ini)
    use MatrixUtils
    use multipoleinfo
    class(CosmoJBDLikelihood) this
    type(TSettingIni) :: Ini
    external DGETRF, DGETRI
    type(TSettingIni) :: Datasets_Mp
    Type(TTextFile) :: F
    REAL(mcp) :: N_s,n_b,Scaling_Factor,k_value,dummy1
    REAL(mcp), allocatable, dimension(:,:) :: Cov_matrix_NGC,Cov_matrix_SGC
    INTEGER, allocatable, dimension(:) :: IPIV
    REAL(mcp), allocatable, dimension(:) :: WORK_INVERSE,k_grid
    INTEGER :: inum,inum2,iopb, ios, temp_num, INFO,LWORK,it,jt
        character(LEN=:), allocatable :: CMASS_multipoles_NGC, CMASS_multipoles_SGC, CMASS_NGC_Cov,CMASS_SGC_Cov

        this%cmass_mp%z_eff = 0.57
        this%cmass_mp%k_spacing = 0.005
        this%cmass_mp%k_min = 0.0125

    if(this%cmass_mp%fit_k20_cmass) then
        this%cmass_mp%k_num = 38
        this%cmass_mp%size_cov = 76
    endif

    if(this%cmass_mp%fit_k15_cmass) then
        this%cmass_mp%k_num = 28
        this%cmass_mp%size_cov = 56
    endif

    if(this%cmass_mp%fit_k10_cmass) then
        this%cmass_mp%k_num = 18
        this%cmass_mp%size_cov = 36
    endif

!----------------------------------------------------------------------------------------------------------
!  k range to compute theory multipoles
!----------------------------------------------------------------------------------------------------------

    this%cmass_mp%kmin_theory = 0.000750
    this%cmass_mp%dk_theory = 0.0005
    this%cmass_mp%kmax_theory = 0.4100

!----------------------------------------------------------------------------------------------------------
!  Parameters required for convolution matrix
!----------------------------------------------------------------------------------------------------------

    this%cmass_mp%numk_values_conv= 150
    this%cmass_mp%numk_values_full_conv= 300
    this%cmass_mp%k_min_conv_matrix = 0.001
    this%cmass_mp%dk_conv_matrix = 0.001

    allocate(this%cmass_mp%P0_multipole_NGC(this%cmass_mp%k_num))
    allocate(this%cmass_mp%P2_multipole_NGC(this%cmass_mp%k_num))
    allocate(this%cmass_mp%P4_multipole_NGC(this%cmass_mp%k_num))
    allocate(this%cmass_mp%P0_multipole_SGC(this%cmass_mp%k_num))
    allocate(this%cmass_mp%P2_multipole_SGC(this%cmass_mp%k_num))
    allocate(this%cmass_mp%P4_multipole_SGC(this%cmass_mp%k_num))
    allocate(INV_Cov_matrix_NGC(this%cmass_mp%size_cov,this%cmass_mp%size_cov))
    allocate(INV_Cov_matrix_SGC(this%cmass_mp%size_cov,this%cmass_mp%size_cov))
    allocate(k_grid(this%cmass_mp%k_num))
    allocate(Cov_matrix_NGC(this%cmass_mp%size_cov,this%cmass_mp%size_cov))
    allocate(Cov_matrix_SGC(this%cmass_mp%size_cov,this%cmass_mp%size_cov))
    allocate(WORK_INVERSE(this%cmass_mp%size_cov*this%cmass_mp%size_cov))
    allocate(IPIV(this%cmass_mp%size_cov + 1))

        !----------------------------------------------------------------------------------------------------------
    ! Import interpolation grids for window functions + setup 1D and 2D splines + compute convolution matrices
    !----------------------------------------------------------------------------------------------------------

        call Datasets_Mp%Open(this%cmass_mp%dataset_dir)
	
    call multipoleinfo_init(Datasets_Mp)

        if(this%cmass_mp%Use_Conv_Matrix) stop 'not currently a feature in this code -- see old code'

        !----------------------------------------------------------------------------------------------------------
    ! Import P0,P2,P4 + covariance matrices
    !----------------------------------------------------------------------------------------------------------
	
    CMASS_multipoles_NGC = Datasets_Mp%ReadFileName('CMASS_NGC')
    CMASS_multipoles_SGC = Datasets_Mp%ReadFileName('CMASS_SGC')

    if(this%cmass_mp%fit_k20_cmass) then
        CMASS_NGC_Cov = Datasets_Mp%ReadFileName('NGC_Cov')
        CMASS_SGC_Cov = Datasets_Mp%ReadFileName('SGC_Cov')
    endif
    if(this%cmass_mp%fit_k15_cmass) then
        CMASS_NGC_Cov = Datasets_Mp%ReadFileName('NGC_Cov_015')
        CMASS_SGC_Cov = Datasets_Mp%ReadFileName('SGC_Cov_015')
    endif
    if(this%cmass_mp%fit_k10_cmass) then
        CMASS_NGC_Cov = Datasets_Mp%ReadFileName('NGC_Cov_010')
        CMASS_SGC_Cov = Datasets_Mp%ReadFileName('SGC_Cov_010')
    endif
	
	call Datasets_Mp%Close()
	
        call F%Open(CMASS_multipoles_NGC)
        do it = 1,this%cmass_mp%k_num
            read (F%unit,*) k_grid(it),this%cmass_mp%P0_multipole_NGC(it),this%cmass_mp%P2_multipole_NGC(it), &
                this%cmass_mp%P4_multipole_NGC(it),dummy1,dummy1
        end do

        call F%Open(CMASS_multipoles_SGC)
        do it = 1,this%cmass_mp%k_num
            read (F%unit,*) k_grid(it),this%cmass_mp%P0_multipole_SGC(it),this%cmass_mp%P2_multipole_SGC(it), &
                this%cmass_mp%P4_multipole_SGC(it),dummy1,dummy1
        end do

    call F%Open(CMASS_NGC_Cov)
    do it=1,this%cmass_mp%size_cov
        do jt=1,this%cmass_mp%size_cov
            read (F%unit,*) dummy1,dummy1,Cov_matrix_NGC(it,jt)
        end do
    end do

    call F%Open(CMASS_SGC_Cov)
    do it=1,this%cmass_mp%size_cov
        do jt=1,this%cmass_mp%size_cov
           read (F%unit,*) dummy1,dummy1,Cov_matrix_SGC(it,jt)
        end do
    end do
	
	!----------------------------------------------------------------------------------------------------------
    ! Find Inverse of NHC and SGC covariance matrices
    !----------------------------------------------------------------------------------------------------------

    LWORK = this%cmass_mp%size_cov**2

    do it = 1,this%cmass_mp%size_cov
        do jt = 1,this%cmass_mp%size_cov
            INV_Cov_matrix_NGC(it,jt) = Cov_matrix_NGC(it,jt)
            INV_Cov_matrix_SGC(it,jt) = Cov_matrix_SGC(it,jt)
        end do
    end do

     call DGETRF(this%cmass_mp%size_cov,this%cmass_mp%size_cov,INV_Cov_matrix_NGC,this%cmass_mp%size_cov,IPIV,INFO)         !! LU decomposition
     call DGETRI(this%cmass_mp%size_cov,INV_Cov_matrix_NGC,this%cmass_mp%size_cov,IPIV,WORK_INVERSE,LWORK,INFO)    !! Inverse from LU decompostion
     if (INFO .ne. 0) stop 'Problem with matrix inverse calculation 1'

     call DGETRF(this%cmass_mp%size_cov,this%cmass_mp%size_cov,INV_Cov_matrix_SGC,this%cmass_mp%size_cov,IPIV,INFO)
     call DGETRI(this%cmass_mp%size_cov,INV_Cov_matrix_SGC,this%cmass_mp%size_cov,IPIV,WORK_INVERSE,LWORK,INFO)
     if (INFO .ne. 0) stop 'Problem with matrix inverse calculation 2'

        !----------------------------------------------------------------------------------------------------------
    ! Re-scale inverse covariance matrix (given the skewed nature of the Wishart distribution see (F. Beulter et at.))
    !
    ! N_s  Number of mocks
    ! n_b Number of k-bins
    !----------------------------------------------------------------------------------------------------------

    N_s = 999.0
    n_b = this%cmass_mp%size_cov

    Scaling_Factor = (N_s - n_b -2.0)/(N_s - 1.0)

    INV_Cov_matrix_NGC =  Scaling_Factor*INV_Cov_matrix_NGC
    INV_Cov_matrix_SGC =  Scaling_Factor*INV_Cov_matrix_SGC
	
    DEALLOCATE(Cov_matrix_NGC)
    DEALLOCATE(Cov_matrix_SGC)
    DEALLOCATE(WORK_INVERSE)
    DEALLOCATE(IPIV)
    DEALLOCATE(k_grid)
	
    end subroutine Read_In_CMASS_data

    !-----------------------------------------------------------------------------------------------
    !Read in WiggleZ Data
    !-----------------------------------------------------------------------------------------------
    subroutine Read_In_Wigglez_data(this, Ini)
    use MatrixUtils
    external DGETRF, DGETRI
    class(CosmoJBDLikelihood) this
    type(TSettingIni) :: Ini
    type(TSettingIni) :: Datasets_Wig
    Type(TTextFile) :: F
    REAL(mcp) :: N_s,n_b,Scaling_Factor
    INTEGER, allocatable, dimension(:) :: IPIV
    REAL(mcp), allocatable, dimension(:) :: WORK_INVERSE,k_grid
    real(mcp), allocatable, dimension(:,:) :: temp_matrix
    REAL(mcp), allocatable, dimension(:)   :: Error_P0,Error_P2,Error_P4
        INTEGER :: inum,jnum,inumnew,jnumnew,iopb, ios, dummy,INFO,LWORK,region_num,file_num
        character(LEN=:), allocatable :: data_lowz,data_highz,temp1,temp2
        character(LEN=:), allocatable :: cov_lowz,cov_highz,temp_file
        character(LEN=:), allocatable :: Window_highZ,Window_lowZ

        call Datasets_Wig%Open(this%wigglez_mp%dataset_dir)

        !-----------------------------------------------------------------------------------------------
        ! True AP effect on
        !-----------------------------------------------------------------------------------------------

        this%wigglez_mp%use_AP = .true.

        !-----------------------------------------------------------------------------------------------
        ! Sets the fitting range from .ini file
        !-----------------------------------------------------------------------------------------------

        if(this%wigglez_mp%fit_k_030) then
                !print*,'using kmax = 0.3'
                this%wigglez_mp%k_num_obs = 14
                this%wigglez_mp%size_cov = 42
        endif
	if(this%wigglez_mp%fit_k_019) then
                !print*,'using kmax = 0.19'
                this%wigglez_mp%k_num_obs = 9
                this%wigglez_mp%size_cov = 27
        endif
	if(this%wigglez_mp%fit_k_017) then
                !print*,'using kmax = 0.17'
                this%wigglez_mp%k_num_obs = 8
                this%wigglez_mp%size_cov = 24
        endif
	if(this%wigglez_mp%fit_k_015) then
                !print*,'using kmax = 0.15'
                this%wigglez_mp%k_num_obs = 7
                this%wigglez_mp%size_cov = 21
        endif

	!---------------------------------------------------------------------------------------------
        !  Set relevant parameters for wiggleZ -- i.e. fitting range, region num, z_eff, ect.
        !---------------------------------------------------------------------------------------------

        this%wigglez_mp%k_min_theory = 0.01
        this%wigglez_mp%dk_theory = 0.02
        this%wigglez_mp%Num_regions = 6
        this%wigglez_mp%z_eff_low = 0.44
        this%wigglez_mp%z_eff_high = 0.73
        this%wigglez_mp%total_num_k_bins = 45
        this%wigglez_mp%size_convolution = 75
        this%wigglez_mp%k_num_conv = 25

        allocate(Error_P0(this%wigglez_mp%k_num_obs))
        allocate(Error_P2(this%wigglez_mp%k_num_obs))
        allocate(Error_P4(this%wigglez_mp%k_num_obs))
        allocate(this%wigglez_mp%P0_data_lowZ(this%wigglez_mp%k_num_obs,this%wigglez_mp%num_regions))
        allocate(this%wigglez_mp%P2_data_lowZ(this%wigglez_mp%k_num_obs,this%wigglez_mp%num_regions))
        allocate(this%wigglez_mp%P4_data_lowZ(this%wigglez_mp%k_num_obs,this%wigglez_mp%num_regions))
        allocate(this%wigglez_mp%P0_data_highZ(this%wigglez_mp%k_num_obs,this%wigglez_mp%num_regions))
        allocate(this%wigglez_mp%P2_data_highZ(this%wigglez_mp%k_num_obs,this%wigglez_mp%num_regions))
        allocate(this%wigglez_mp%P4_data_highZ(this%wigglez_mp%k_num_obs,this%wigglez_mp%num_regions))
        allocate(this%wigglez_mp%INV_Cov_matrix_lowZ(this%wigglez_mp%size_cov,this%wigglez_mp%size_cov,this%wigglez_mp%num_regions))
        allocate(this%wigglez_mp%INV_Cov_matrix_highZ(this%wigglez_mp%size_cov,this%wigglez_mp%size_cov,this%wigglez_mp%num_regions))
        allocate(this%wigglez_mp%Covolution_matrix_lowZ(this%wigglez_mp%size_convolution,this%wigglez_mp%size_convolution,this%wigglez_mp%num_regions))
        allocate(this%wigglez_mp%Covolution_matrix_highZ(this%wigglez_mp%size_convolution,this%wigglez_mp%size_convolution,this%wigglez_mp%num_regions))
        allocate(k_grid(this%wigglez_mp%k_num_obs))
        allocate(WORK_INVERSE(this%wigglez_mp%size_cov*this%wigglez_mp%size_cov))
        allocate(IPIV(this%wigglez_mp%size_cov + 1))


        !--------------------------------------------------------------------------------------
        !  import P0,P2,P4 measurement -- for each of the 6 regions
        !---------------------------------------------------------------------------------------

        do region_num = 1,this%wigglez_mp%num_regions

                if(region_num .EQ. 1) then
                   data_lowz =  Datasets_Wig%ReadFileName('region1_z0pt4_data')
                   data_highz = Datasets_Wig%ReadFileName('region1_z0pt8_data')
                endif
                if(region_num .EQ. 2) then
                   data_lowz = Datasets_Wig%ReadFileName('region2_z0pt4_data')
                   data_highz = Datasets_Wig%ReadFileName('region2_z0pt8_data')
                endif
                if(region_num .EQ. 3) then
                   data_lowz = Datasets_Wig%ReadFileName('region3_z0pt4_data')
                   data_highz = Datasets_Wig%ReadFileName('region3_z0pt8_data')
                endif
                if(region_num .EQ. 4) then
                   data_lowz =  Datasets_Wig%ReadFileName('region4_z0pt4_data')
                   data_highz = Datasets_Wig%ReadFileName('region4_z0pt8_data')
                endif
                if(region_num .EQ. 5) then
                   data_lowz = Datasets_Wig%ReadFileName('region5_z0pt4_data')
                   data_highz = Datasets_Wig%ReadFileName('region5_z0pt8_data')
                endif
                if(region_num .EQ. 6) then
                   data_lowz = Datasets_Wig%ReadFileName('region6_z0pt4_data')
                   data_highz = Datasets_Wig%ReadFileName('region6_z0pt8_data')
                endif

                call F%Open(data_lowz)
                do inum=1,this%wigglez_mp%k_num_obs
                        read (F%unit,*) k_grid(inum),this%wigglez_mp%P0_data_lowZ(inum,region_num),Error_P0(inum), &
                        this%wigglez_mp%P2_data_lowZ(inum,region_num), Error_P2(inum),this%wigglez_mp%P4_data_lowZ(inum,region_num),Error_P4(inum)
                end do

                call F%Open(data_highz)
                do inum=1,this%wigglez_mp%k_num_obs
                        read (F%unit,*) k_grid(inum),this%wigglez_mp%P0_data_highZ(inum,region_num),Error_P0(inum), &
                        this%wigglez_mp%P2_data_highZ(inum,region_num),Error_P2(inum),this%wigglez_mp%P4_data_highZ(inum,region_num),Error_P4(inum)
                end do

        end do

        !--------------------------------------------------------------------------------
        !  import Covariance matrix -- 6 regions * 2 redshifts
        !--------------------------------------------------------------------------------

        IF(this%wigglez_mp%fit_k_030) file_num=1
        IF(this%wigglez_mp%fit_k_019) file_num=2
        IF(this%wigglez_mp%fit_k_017) file_num=3
        IF(this%wigglez_mp%fit_k_015) file_num=4
	
        DO region_num = 1,this%wigglez_mp%num_regions
                if(region_num .EQ. 1) then
                        cov_lowz =  Datasets_Wig%ReadFileName(numcat_local('region1_z0pt4_Cov',file_num))
                        cov_highz = Datasets_Wig%ReadFileName(numcat_local('region1_z0pt8_Cov',file_num))
                endif
                if(region_num .EQ. 2) then
                  cov_lowz =  Datasets_Wig%ReadFileName(numcat_local('region2_z0pt4_Cov',file_num))
                  cov_highz = Datasets_Wig%ReadFileName(numcat_local('region2_z0pt8_Cov',file_num))
                endif
                if(region_num .EQ. 3) then
                        cov_lowz =  Datasets_Wig%ReadFileName(numcat_local('region3_z0pt4_Cov',file_num))
                        cov_highz = Datasets_Wig%ReadFileName(numcat_local('region3_z0pt8_Cov',file_num))
                endif
                if(region_num .EQ. 4) then
                        cov_lowz =  Datasets_Wig%ReadFileName(numcat_local('region4_z0pt4_Cov',file_num))
                        cov_highz = Datasets_Wig%ReadFileName(numcat_local('region4_z0pt8_Cov',file_num))
                endif
                if(region_num .EQ. 5) then
                        cov_lowz = Datasets_Wig%ReadFileName(numcat_local('region5_z0pt4_Cov',file_num))
                        cov_highz = Datasets_Wig%ReadFileName(numcat_local('region5_z0pt8_Cov',file_num))
                endif
                if(region_num .EQ. 6) then
                          cov_lowz = Datasets_Wig%ReadFileName(numcat_local('region6_z0pt4_Cov',file_num))
                          cov_highz = Datasets_Wig%ReadFileName(numcat_local('region6_z0pt8_Cov',file_num))
                endif

                call F%Open(cov_lowz)
                do inum=1,this%wigglez_mp%size_cov
                        do jnum=1,this%wigglez_mp%size_cov
                                read (F%unit,*) dummy,dummy,this%wigglez_mp%INV_Cov_matrix_lowZ(inum,jnum,region_num)
                        end do
                end do

                call F%Open(cov_highz)
                do inum=1,this%wigglez_mp%size_cov
                        do jnum=1,this%wigglez_mp%size_cov
                                read (F%unit,*) dummy,dummy,this%wigglez_mp%INV_Cov_matrix_highZ(inum,jnum,region_num)
                        end do
                end do
        end do
	
        !-------------------------------------------------------------------------------------
        !  import convolution matrix -- 6 regions * 2 redshifts
        !-------------------------------------------------------------------------------------

        do region_num = 1,this%wigglez_mp%num_regions

                if(region_num .EQ. 1) then
                   Window_lowZ =  Datasets_Wig%ReadFileName('region1_z0pt4_WindowFn')
                   Window_highZ = Datasets_Wig%ReadFileName('region1_z0pt8_WindowFn')
                endif
                if(region_num .EQ. 2) then
                   Window_lowZ = Datasets_Wig%ReadFileName('region2_z0pt4_WindowFn')
                   Window_highZ = Datasets_Wig%ReadFileName('region2_z0pt8_WindowFn')
                endif
                if(region_num .EQ. 3) then
                   Window_lowZ = Datasets_Wig%ReadFileName('region3_z0pt4_WindowFn')
                   Window_highZ = Datasets_Wig%ReadFileName('region3_z0pt8_WindowFn')
                endif
                if(region_num .EQ. 4) then
                   Window_lowZ = Datasets_Wig%ReadFileName('region4_z0pt4_WindowFn')
                   Window_highZ = Datasets_Wig%ReadFileName('region4_z0pt8_WindowFn')
                endif
                if(region_num .EQ. 5) then
                   Window_lowZ = Datasets_Wig%ReadFileName('region5_z0pt4_WindowFn')
                   Window_highZ = Datasets_Wig%ReadFileName('region5_z0pt8_WindowFn')
                endif
                if(region_num .EQ. 6) then
                   Window_lowZ = Datasets_Wig%ReadFileName('region6_z0pt4_WindowFn')
                   Window_highZ = Datasets_Wig%ReadFileName('region6_z0pt8_WindowFn')
                endif

                call F%Open(Window_lowZ)
                do inum=1,this%wigglez_mp%size_convolution
                        do jnum=1,this%wigglez_mp%size_convolution
                                read (F%unit,*) dummy,dummy,this%wigglez_mp%Covolution_matrix_lowZ(inum,jnum,region_num)
                        end do
                end do

                call F%Open(Window_highZ)
                do inum=1,this%wigglez_mp%size_convolution
                        do jnum=1,this%wigglez_mp%size_convolution
                                read (F%unit,*) dummy,dummy,this%wigglez_mp%Covolution_matrix_highZ(inum,jnum,region_num)
                        end do
                end do

        enddo

	call Datasets_Wig%Close()

        !---------------------------------------------------------------------------------------------------------
        !  Invert all covariance matrices using LU decomposition
        !---------------------------------------------------------------------------------------------------------

        DO region_num = 1,this%wigglez_mp%num_regions
                LWORK = this%wigglez_mp%size_cov**2
                call DGETRF(this%wigglez_mp%size_cov,this%wigglez_mp%size_cov,this%wigglez_mp%INV_Cov_matrix_lowZ(:,:,region_num),this%wigglez_mp%size_cov,IPIV,INFO)
                call DGETRI(this%wigglez_mp%size_cov,this%wigglez_mp%INV_Cov_matrix_lowZ(:,:,region_num),this%wigglez_mp%size_cov,IPIV,WORK_INVERSE,LWORK,INFO)
                IPIV(:) =  0
                WORK_INVERSE(:) = 0
                if (INFO .ne. 0) stop 'Problem with matrix inverse calculation low-Z'
        ENDDO

	DO region_num = 1,this%wigglez_mp%num_regions
                LWORK = this%wigglez_mp%size_cov**2
                call DGETRF(this%wigglez_mp%size_cov,this%wigglez_mp%size_cov,this%wigglez_mp%INV_Cov_matrix_highZ(:,:,region_num),this%wigglez_mp%size_cov,IPIV,INFO)
                call DGETRI(this%wigglez_mp%size_cov,this%wigglez_mp%INV_Cov_matrix_highZ(:,:,region_num),this%wigglez_mp%size_cov,IPIV,WORK_INVERSE,LWORK,INFO)
                IPIV(:) =  0
                WORK_INVERSE(:) = 0
                if (INFO .ne. 0) stop 'Problem with matrix inverse calculation high Z'
        ENDDO

	!----------------------------------------------------------------------------------------------------------
        !  Re-scale inverse covariance matrix (given the skewed nature of the Wishart distribution see (F. Beulter et at.))
        !----------------------------------------------------------------------------------------------------------

        N_s = 600.0
        n_b = this%wigglez_mp%size_cov
        Scaling_Factor = (N_s - n_b -2.0)/(N_s - 1.0)

        DO region_num = 1,this%wigglez_mp%num_regions
                this%wigglez_mp%INV_Cov_matrix_lowZ(:,:,region_num) =  Scaling_Factor*this%wigglez_mp%INV_Cov_matrix_lowZ(:,:,region_num)
                this%wigglez_mp%INV_Cov_matrix_highZ(:,:,region_num) =  Scaling_Factor*this%wigglez_mp%INV_Cov_matrix_highZ(:,:,region_num)
        ENDDO

	DEALLOCATE(Error_P0)
        DEALLOCATE(Error_P2)
        DEALLOCATE(Error_P4)
        DEALLOCATE(WORK_INVERSE)
        DEALLOCATE(IPIV)

    end subroutine Read_In_Wigglez_data

    !-----------------------------------------------------------------------------------------------
    ! Read in data for the overlapping region of CMASS -- overlap with KiDS
    !-----------------------------------------------------------------------------------------------
    subroutine Read_In_CMASS_data_Overlap(this, Ini)
    use MatrixUtils
    external DGETRF, DGETRI
    class(CosmoJBDLikelihood) this
    type(TSettingIni) :: Ini
    type(TSettingIni) :: Datasets_CMASS_overlap
    Type(TTextFile) :: F
    REAL(mcp) :: temp_error, temp_error2
    REAL(mcp), allocatable, dimension(:) :: k_grid
    REAL(mcp), allocatable, dimension(:)   :: Error_P0,Error_P2
        INTEGER :: inum,jnum,inumnew,jnumnew,iopb, ios, dummy,INFO,LWORK,region_num,file_num
        character(LEN=:), allocatable :: data_main,Window_fn
	
        call Datasets_CMASS_overlap%Open(this%cmass_mp_overlap%dataset_dir)
	
        !print*,'Reading in CMASS data for the overlap region'
	
        !-----------------------------------------------------------------------------------------------
        ! Turn AP effect on
        !-----------------------------------------------------------------------------------------------
	
        this%cmass_mp_overlap%use_AP = .true.
	
        !-----------------------------------------------------------------------------------------------
        ! Sets the fitting range from .ini file
        !-----------------------------------------------------------------------------------------------

        if(this%cmass_mp_overlap%fit_k_0175) then
                !print*,'using kmax = 0.175 for CMASS overlap'
                this%cmass_mp_overlap%k_num_obs = 3
                this%cmass_mp_overlap%size_cov = 9
        endif
	if(this%cmass_mp_overlap%fit_k_0125) then
                !print*,'using kmax = 0.125 for CMASS overlap'
                this%cmass_mp_overlap%k_num_obs = 2
                this%cmass_mp_overlap%size_cov = 4
        endif
	if(this%cmass_mp_overlap%fit_k_0075) then
                !print*,'using kmax = 0.075 for CMASS overlap'
                this%cmass_mp_overlap%k_num_obs = 1
                this%cmass_mp_overlap%size_cov = 1
        endif

	!---------------------------------------------------------------------------------------------
        !  Set relevant parameters for CMASS overlap -- i.e. fitting range, region num, z_eff, ect.
        !---------------------------------------------------------------------------------------------

        this%cmass_mp_overlap%k_min_theory = 0.00
        this%cmass_mp_overlap%dk_theory = 0.05
        this%cmass_mp_overlap%z_eff = 0.57
        this%cmass_mp_overlap%size_convolution = 30
        this%cmass_mp_overlap%k_num_conv = 10
        this%cmass_mp_overlap%k_spacing_obs = 0.05
        this%cmass_mp_overlap%k_min_obs = 0.075
	
        allocate(Error_P0(this%cmass_mp_overlap%k_num_obs))
        allocate(Error_P2(this%cmass_mp_overlap%k_num_obs))
        allocate(this%cmass_mp_overlap%Covolution_matrix(this%cmass_mp_overlap%size_convolution,this%cmass_mp_overlap%size_convolution))
        allocate(k_grid(this%cmass_mp_overlap%k_num_obs))

        !-------------------------------------------------------------------------------------
        !  import convolution matrix -- 1 region, 1 redshift
        !-------------------------------------------------------------------------------------

        Window_fn =  Datasets_CMASS_overlap%ReadFileName('CMASS_overlap_conv')

        call F%Open(Window_fn)
        do inum=1,this%cmass_mp_overlap%size_convolution
                do jnum=1,this%cmass_mp_overlap%size_convolution
                        read (F%unit,*) dummy,dummy,this%cmass_mp_overlap%Covolution_matrix(inum,jnum)
                end do
        end do

        call Datasets_CMASS_overlap%Close()
    end subroutine Read_In_CMASS_data_Overlap

    !-----------------------------------------------------------------------------------------------
    ! Read in data for the overlapping region of LOWZ -- overlap with KiDS
    !-----------------------------------------------------------------------------------------------
    subroutine Read_In_LOWZ_data_Overlap(this, Ini)
    use MatrixUtils
    external DGETRF, DGETRI
    class(CosmoJBDLikelihood) this
    type(TSettingIni) :: Ini
    type(TSettingIni) :: Datasets_LOWZ_overlap
    Type(TTextFile) :: F
    REAL(mcp) :: temp_error,temp_error2
    REAL(mcp), allocatable, dimension(:) :: k_grid
    REAL(mcp), allocatable, dimension(:)   :: Error_P0,Error_P2
        INTEGER :: inum,jnum,inumnew,jnumnew,iopb, ios, dummy,INFO,LWORK,region_num,file_num
        character(LEN=:), allocatable :: data_main,Window_fn
	
        call Datasets_LOWZ_overlap%Open(this%lowz_mp_overlap%dataset_dir)

        !print*,'Reading in LOWZ data for the overlap region'
	
        !-----------------------------------------------------------------------------------------------
        ! Turn AP effect on
        !-----------------------------------------------------------------------------------------------
	
        this%lowz_mp_overlap%use_AP = .true.

        !-----------------------------------------------------------------------------------------------
        ! Sets the fitting range from .ini file
        !-----------------------------------------------------------------------------------------------

        if(this%lowz_mp_overlap%fit_k_0175) then
                !print*,'using kmax = 0.175 for LOWZ overlap'
                this%lowz_mp_overlap%k_num_obs = 3
                this%lowz_mp_overlap%size_cov = 9
        endif
	if(this%lowz_mp_overlap%fit_k_0125) then
                !print*,'using kmax = 0.125 for LOWZ overlap'
                this%lowz_mp_overlap%k_num_obs = 2
                this%lowz_mp_overlap%size_cov = 4
        endif
	if(this%lowz_mp_overlap%fit_k_0075) then
                !print*,'using kmax = 0.075 for LOWZ overlap'
                this%lowz_mp_overlap%k_num_obs = 1
                this%lowz_mp_overlap%size_cov = 1
        endif

	!---------------------------------------------------------------------------------------------
        !  Set relevant parameters for lowz overlap -- i.e. fitting range, region num, z_eff, ect.
        !---------------------------------------------------------------------------------------------

        this%lowz_mp_overlap%k_min_theory = 0.00 ! Should be 0, numerical issues?
        this%lowz_mp_overlap%dk_theory = 0.05
        this%lowz_mp_overlap%z_eff = 0.32
        this%lowz_mp_overlap%size_convolution = 30
        this%lowz_mp_overlap%k_num_conv = 10
        this%lowz_mp_overlap%k_spacing_obs = 0.05
        this%lowz_mp_overlap%k_min_obs = 0.075

        allocate(Error_P0(this%lowz_mp_overlap%k_num_obs))
        allocate(Error_P2(this%lowz_mp_overlap%k_num_obs))
        allocate(this%lowz_mp_overlap%Covolution_matrix(this%lowz_mp_overlap%size_convolution,this%lowz_mp_overlap%size_convolution))
        allocate(k_grid(this%lowz_mp_overlap%k_num_obs))

        !-------------------------------------------------------------------------------------
        !  import convolution matrix -- 1 region, 1 redshift
        !-------------------------------------------------------------------------------------

        Window_fn =  Datasets_LOWZ_overlap%ReadFileName('LOWZ_overlap_conv')

        call F%Open(Window_fn)
        do inum=1,this%lowz_mp_overlap%size_convolution
                do jnum=1,this%lowz_mp_overlap%size_convolution
                        read (F%unit,*) dummy,dummy,this%lowz_mp_overlap%Covolution_matrix(inum,jnum)
                end do
        end do

        call Datasets_LOWZ_overlap%Close()
	
    end subroutine Read_In_LOWZ_data_Overlap


    !-----------------------------------------------------------------------------------------------
    ! Read in data for the overlapping region of 2dfloz -- overlap with KiDS
    !-----------------------------------------------------------------------------------------------
    subroutine Read_In_2dfloz_data_Overlap(this, Ini)
    use MatrixUtils
    external DGETRF, DGETRI
    class(CosmoJBDLikelihood) this
    type(TSettingIni) :: Ini
    type(TSettingIni) :: Datasets_2dfloz_overlap
    Type(TTextFile) :: F
    REAL(mcp) :: temp_error,temp_error2
    REAL(mcp), allocatable, dimension(:) :: k_grid
    REAL(mcp), allocatable, dimension(:)   :: Error_P0,Error_P2
        INTEGER :: inum,jnum,inumnew,jnumnew,iopb, ios, dummy,INFO,LWORK,region_num,file_num
        character(LEN=:), allocatable :: data_main,Window_fn
	
        call Datasets_2dfloz_overlap%Open(this%twodfloz_mp_overlap%dataset_dir)

        !-----------------------------------------------------------------------------------------------
        ! Turn AP effect on
        !-----------------------------------------------------------------------------------------------
	
        this%twodfloz_mp_overlap%use_AP = .true.

        !-----------------------------------------------------------------------------------------------
        ! Sets the fitting range from .ini file
        !-----------------------------------------------------------------------------------------------

        if(this%twodfloz_mp_overlap%fit_k_0175) then
                !print*,'using kmax = 0.175 for 2dfloz overlap'
                this%twodfloz_mp_overlap%k_num_obs = 3
                this%twodfloz_mp_overlap%size_cov = 9
        endif
	if(this%twodfloz_mp_overlap%fit_k_0125) then
                !print*,'using kmax = 0.125 for 2dfloz overlap'
                this%twodfloz_mp_overlap%k_num_obs = 2
                this%twodfloz_mp_overlap%size_cov = 4
        endif
	if(this%twodfloz_mp_overlap%fit_k_0075) then
                !print*,'using kmax = 0.075 for 2dfloz overlap'
                this%twodfloz_mp_overlap%k_num_obs = 1
                this%twodfloz_mp_overlap%size_cov = 1
        endif

	!---------------------------------------------------------------------------------------------
        !  Set relevant parameters for 2dfloz overlap -- i.e. fitting range, region num, z_eff, ect.
        !---------------------------------------------------------------------------------------------

        this%twodfloz_mp_overlap%k_min_theory = 0.00 ! Should be 0, numerical issues?
        this%twodfloz_mp_overlap%dk_theory = 0.05
        this%twodfloz_mp_overlap%z_eff = 0.31
        this%twodfloz_mp_overlap%size_convolution = 30
        this%twodfloz_mp_overlap%k_num_conv = 10
        this%twodfloz_mp_overlap%k_spacing_obs = 0.05
        this%twodfloz_mp_overlap%k_min_obs = 0.075

        allocate(Error_P0(this%twodfloz_mp_overlap%k_num_obs))
        allocate(Error_P2(this%twodfloz_mp_overlap%k_num_obs))
        allocate(this%twodfloz_mp_overlap%Covolution_matrix(this%twodfloz_mp_overlap%size_convolution,this%twodfloz_mp_overlap%size_convolution))
        allocate(k_grid(this%twodfloz_mp_overlap%k_num_obs))

        !-------------------------------------------------------------------------------------
        !  import convolution matrix -- 1 region, 1 redshift
        !-------------------------------------------------------------------------------------

        Window_fn =  Datasets_2dfloz_overlap%ReadFileName('2dfloz_overlap_conv')

        call F%Open(Window_fn)
        do inum=1,this%twodfloz_mp_overlap%size_convolution
                do jnum=1,this%twodfloz_mp_overlap%size_convolution
                        read (F%unit,*) dummy,dummy,this%twodfloz_mp_overlap%Covolution_matrix(inum,jnum)
                end do
        end do

        call Datasets_2dfloz_overlap%Close()
	
    end subroutine Read_In_2dfloz_data_Overlap


    !-----------------------------------------------------------------------------------------------
    ! Read in data for the overlapping region of 2dfhiz -- overlap with KiDS
    !-----------------------------------------------------------------------------------------------
    subroutine Read_In_2dfhiz_data_Overlap(this, Ini)
    use MatrixUtils
    external DGETRF, DGETRI
    class(CosmoJBDLikelihood) this
    type(TSettingIni) :: Ini
    type(TSettingIni) :: Datasets_2dfhiz_overlap
    Type(TTextFile) :: F
    REAL(mcp) :: temp_error,temp_error2
    REAL(mcp), allocatable, dimension(:) :: k_grid
    REAL(mcp), allocatable, dimension(:)   :: Error_P0,Error_P2
        INTEGER :: inum,jnum,inumnew,jnumnew,iopb, ios, dummy,INFO,LWORK,region_num,file_num
        character(LEN=:), allocatable :: data_main,Window_fn
	
        call Datasets_2dfhiz_overlap%Open(this%twodfhiz_mp_overlap%dataset_dir)

        !print*,'Reading in 2dfhiz data for the overlap region'
	
        !-----------------------------------------------------------------------------------------------
        ! Turn AP effect on
        !-----------------------------------------------------------------------------------------------
	
        this%twodfhiz_mp_overlap%use_AP = .true.

        !-----------------------------------------------------------------------------------------------
        ! Sets the fitting range from .ini file
        !-----------------------------------------------------------------------------------------------

        if(this%twodfhiz_mp_overlap%fit_k_0175) then
                !print*,'using kmax = 0.175 for 2dfhiz overlap'
                this%twodfhiz_mp_overlap%k_num_obs = 3
                this%twodfhiz_mp_overlap%size_cov = 9
        endif
	if(this%twodfhiz_mp_overlap%fit_k_0125) then
                !print*,'using kmax = 0.125 for 2dfhiz overlap'
                this%twodfhiz_mp_overlap%k_num_obs = 2
                this%twodfhiz_mp_overlap%size_cov = 4
        endif
	if(this%twodfhiz_mp_overlap%fit_k_0075) then
                !print*,'using kmax = 0.075 for 2dfhiz overlap'
                this%twodfhiz_mp_overlap%k_num_obs = 1
                this%twodfhiz_mp_overlap%size_cov = 1
        endif

	!---------------------------------------------------------------------------------------------
        !  Set relevant parameters for 2dfhiz overlap -- i.e. fitting range, region num, z_eff, ect.
        !---------------------------------------------------------------------------------------------

        this%twodfhiz_mp_overlap%k_min_theory = 0.00 ! Should be 0, numerical issues?
        this%twodfhiz_mp_overlap%dk_theory = 0.05
        this%twodfhiz_mp_overlap%z_eff = 0.56
        this%twodfhiz_mp_overlap%size_convolution = 30
        this%twodfhiz_mp_overlap%k_num_conv = 10
        this%twodfhiz_mp_overlap%k_spacing_obs = 0.05
        this%twodfhiz_mp_overlap%k_min_obs = 0.075

        allocate(Error_P0(this%twodfhiz_mp_overlap%k_num_obs))
        allocate(Error_P2(this%twodfhiz_mp_overlap%k_num_obs))
        allocate(this%twodfhiz_mp_overlap%Covolution_matrix(this%twodfhiz_mp_overlap%size_convolution,this%twodfhiz_mp_overlap%size_convolution))
        allocate(k_grid(this%twodfhiz_mp_overlap%k_num_obs))

        !-------------------------------------------------------------------------------------
        !  import convolution matrix -- 1 region, 1 redshift
        !-------------------------------------------------------------------------------------

        Window_fn =  Datasets_2dfhiz_overlap%ReadFileName('2dfhiz_overlap_conv')

        call F%Open(Window_fn)
        do inum=1,this%twodfhiz_mp_overlap%size_convolution
                do jnum=1,this%twodfhiz_mp_overlap%size_convolution
                        read (F%unit,*) dummy,dummy,this%twodfhiz_mp_overlap%Covolution_matrix(inum,jnum)
                end do
        end do

        call Datasets_2dfhiz_overlap%Close()
	
    end subroutine Read_In_2dfhiz_data_Overlap


    !----------------------------------------------------------------------------------------------------------
    ! Main likelihood calculation
    !----------------------------------------------------------------------------------------------------------
    function CosmoJBD_LnLike(this,CMB,Theory,DataParams)
    implicit none
    class(CosmoJBDLikelihood) :: this
    class(CMBParams) :: CMB
    class(TCosmoTheoryPredictions), target :: Theory
    type my_type
        type(TCosmoTheoryPredictions), allocatable :: theoryomp
        real(mcp) :: pisjomp,h0omp,homp,omdmomp,ombomp,omkomp,distlensflatomp,mumin2vomp,mumaxomp,tol_erroromp,ampiayesomp,redziayesomp,lumiayesomp,nellbinsomp,bbiaz0yesomp,bbiaz1yesomp,bbiaz2yesomp,bbiaz3yesomp,bbiazconstomp
        real(mcp), allocatable, dimension(:) :: ellarromp, aphotarromp,exact_zmodomp,lumarromp,exact_zmodomplenslowz,exact_zmodomplenscmass
        integer, allocatable, dimension(:) :: exact_z_index,momp,m1arromp,m2arromp
        integer :: wchooseomp,ellgenomp,bnumomp,bnumomplens,uselensomp,kline,m1omp,m2omp,binwomp,wtrapmaxomp,wtrapmaxlensomp,setscenarioomp,eftflag,designereftmodel
        real(mcp) :: mumax1omp,mumax2omp,mumaxlensomp,mumaxwomp,mumaxcrossomplowz,mumin2vcrossomplowz,mumaxcrossompcmass,mumin2vcrossompcmass
        real(mcp) :: meadnuis1omp,meadnuis2omp
        real(mcp), allocatable, dimension(:,:) :: weightarromp,weightpsarromp,weightnoarromp,weightarromplowz,weightarrompcmass,weightpsarromplowz,weightpsarrompcmass,weightnoarromplenslowz,weightnoarromplenscmass,weightarromp2dfloz,weightarromp2dfhiz,weightpsarromp2dfloz,weightpsarromp2dfhiz,weightnoarromplens2dfloz,weightnoarromplens2dfhiz
        real(mcp), allocatable, dimension(:) :: weightarromplenslowz,weightarromplenscmass,weightarromplens2dfloz,weightarromplens2dfhiz
        type(SCubicSpline),allocatable, dimension(:) :: psourcetypeomp
        type(SCubicSpline), allocatable :: plenstypeomplowz,plenstypeompcmass,plenstypeomp2dfloz,plenstypeomp2dfhiz
        type(SCubicSpline), allocatable :: gcubicomp
        Type(TCubicSpline), allocatable :: clinterptypeomp,clinterptypeiiomp,clinterptypegiomp,clinterptypecrosslowzomp,clinterptypecrosscmassomp,clinterptypecross2dflozomp,clinterptypecross2dfhizomp,clinterptypecrosslowzgiomp,clinterptypecrosscmassgiomp,clinterptypecross2dflozgiomp,clinterptypecross2dfhizgiomp !keep this as tcubicspline
    end type my_type
    real(mcp) rombint
    real(mcp) rombint_obj
    external rombint
    external rombint_obj
    type(my_type) obj
    real(mcp) :: CosmoJBD_LnLike,ampiayes,redziayes,lumiayes
    real(mcp) :: tol_error
    real(mcp) :: sigma_v1,sigma_v2,b1_bias_lowZ,b1_bias_highZ,bbiaz0yes,bbiaz1yes,bbiaz2yes,bbiaz3yes
    real(mcp) :: a1phot,a2phot,a3phot,a4phot,a5phot,a6phot,a7phot
    real(mcp) :: meadnuis1,meadnuis2
    real(mcp) :: mumax,mumin,mumin2v, pisj,likesj,chisqsj,chisqsjsellentin,start,finish,mumax1,mumax2,mumaxw,morellfac,ckmz,fitta1,fitta2,fittdiff,fittdist,mumaxcrosslowz,mumin2vcrosslowz,mumaxcrosscmass,mumin2vcrosscmass,startrsd
    INTEGER   :: zim,bbint,bbintigen,bbintmera,bbintmeraigen
    real(mcp), allocatable, dimension(:) :: xiplus,xiplusii,xiplusgi,xicrosslowz,xicrosscmass,xicrosslowzgi,xicrosscmassgi,xicross2dfloz,xicross2dfhiz,xicross2dflozgi,xicross2dfhizgi
    real(mcp), allocatable, dimension(:) :: ximinus,ximinusii,ximinusgi
    real(mcp), allocatable, dimension(:) :: xiplusminus, xiplusminusgg,xiplusminusii,xiplusminusgi,dataminustheory,dataminustheoryhihi,dataminustheoryhoho,finaltheory
    real(mcp), allocatable, dimension(:) :: xicrossarrlowz,xicrossarrlowzgi,xicrossarrcmass,xicrossarrcmassgi,xicrossarr2dfloz,xicrossarr2dflozgi,xicrossarr2dfhiz,xicrossarr2dfhizgi
    real(mcp), allocatable, dimension(:,:) :: weightarr,weightpsarr,weightnoarr,weightnoarrlenslowz,weightnoarrlenscmass,weightarrlowz,weightarrcmass,weightpsarrlowz,weightpsarrcmass,weightnoarrlens2dfloz,weightnoarrlens2dfhiz,weightarr2dfloz,weightarr2dfhiz,weightpsarr2dfloz,weightpsarr2dfhiz
    real(mcp), allocatable, dimension(:) :: weightarrlenslowz,weightarrlenscmass,weightarrlens2dfloz,weightarrlens2dfhiz
    real(mcp), allocatable, dimension(:) :: m1arr,m2arr,binwarr,exact_zmodlenslowz,exact_zmodlenscmass,exact_zmodlens2dfloz,exact_zmodlens2dfhiz
    real(mcp), allocatable, dimension(:) :: ellarr,trapzarr,trapzarrii,trapzarrgi,trapzarrcrosslowz,trapzarrcrosscmass,trapzarrcross2dfloz,trapzarrcross2dfhiz,trapzarrcrosslowzgi,trapzarrcrosscmassgi,trapzarrcross2dflozgi,trapzarrcross2dfhizgi
    real(mcp), allocatable, dimension(:) :: mumax1arr,mumax2arr, aphotarr,garr,rcslum,cfhtlum,lumarr,mnomp
    real(mcp), allocatable, dimension(:,:) :: clarr,clarrii,clarrgi,clarrcrosslowz,clarrcrosscmass,clarrcross2dfloz,clarrcross2dfhiz,clarrcrosslowzgi,clarrcrosscmassgi,clarrcross2dflozgi,clarrcross2dfhizgi
    real(mcp), allocatable, dimension(:) :: clarrbig,ellgentestarr,clarrbigii,clarrbiggi,clarrbigcrosslowz,clarrbigcrosscmass,clarrbigcross2dfloz,clarrbigcross2dfhiz,clarrbigcrosslowzgi,clarrbigcrosscmassgi,clarrbigcross2dflozgi,clarrbigcross2dfhizgi
    integer :: hoj,bnum,gigi,how,howin
    integer :: loopee,yeye,ellgen,intnum,bin,biin,wchoose, ho, m1, m2, mint,ttt,wtrapmax,wtrapmaxlens,wttt,binw
    integer :: nangbins, nzbins, nellbins,  nell, nzp
    real(mcp) :: ellgentest
    logical :: use_cfht = .true.
    logical :: use_rombintz = .true. !if true use Romberg integration instead of trapezoidal for outer integral
    logical :: use_morellz = .false. !if true use 101 ell values, otherwise 31 ell values
    !!CMASS, LOWZ, 2dfloz, 2dfhiz interpolation objects.
    Type(TCubicSpline) ::  P0_theory_spline_cmass,P2_theory_spline_cmass
    Type(TCubicSpline) ::  P0_theory_spline_cmass_overlap,P2_theory_spline_cmass_overlap
    Type(TCubicSpline) ::  P0_theory_spline_lowz_overlap,P2_theory_spline_lowz_overlap
    Type(TCubicSpline) ::  P0_theory_spline_2dfloz_overlap,P2_theory_spline_2dfloz_overlap
    Type(TCubicSpline) ::  P0_theory_spline_2dfhiz_overlap,P2_theory_spline_2dfhiz_overlap
    !! CMASS full variables
    REAL(mcp), allocatable, dimension(:) :: P0_theory_win_NGC,P2_theory_win_NGC
    REAL(mcp), allocatable, dimension(:) :: P0_theory_win_SGC,P2_theory_win_SGC
    REAL(mcp), allocatable, dimension(:) :: diff_NGC,diff_SGC
    REAL(mcp) :: a_perp_cmass, a_par_cmass,z_cmass,k_val_cmass !Need to be global variables
    REAL(mcp) :: scale_factor_L0,scale_factor_L2,scale_factor_L4 !Need to be global variables
    REAL(mcp) :: like_NGC,like_SGC,sigma_v_cmass,b1_cmass,N_shot_cmass,temp_num
    INTEGER   :: inum,k_grid_num_cmass
    !! WiggleZ full variables
    REAL(mcp), allocatable, dimension(:,:) :: P0P2P4_lowz_final_wigz,P0P2P4_highz_final_wigz
    REAL(mcp), allocatable, dimension(:,:) :: P0P2P4_lowz_obs_wigz,P0P2P4_highz_obs_wigz
    REAL(mcp) :: a_perp_lowz_wigz, a_par_lowz_wigz,a_perp_highz_wigz, a_par_highz_wigz !Need to be global variables.
    REAL(mcp) :: Mp_WigZ_LnLike,sig_v1_wigz,sig_v2_wigz,b_lowz_wigz,b_highz_wigz  !Need to be global variables.
    REAL(mcp) :: z_lowz_wigz,z_highz_wigz, k_val_wigz  !global
    REAL(mcp), allocatable, dimension(:)   :: like_lowZ,like_highZ
    REAL(mcp), allocatable, dimension(:,:) :: diff_lowz_wigz,diff_highz_wigz
    INTEGER   :: region_num
    !! CMASS overlap variables
    REAL(mcp), allocatable, dimension(:) :: P0P2_final_cmass_overlap,P0P2_obs_cmass_overlap,diff_cmass_overlap
    REAL(mcp) :: a_perp_cmass_overlap,a_par_cmass_overlap,k_val_overlap_cmass
    !! LOWZ overlap variables
    REAL(mcp), allocatable, dimension(:) :: P0P2_final_lowz_overlap,P0P2_obs_lowz_overlap,diff_lowz_overlap
    REAL(mcp) :: a_perp_lowz_overlap,a_par_lowz_overlap,k_val_overlap_lowz,z_lowz,sigv_lowz,b1_lowz,N_shot_lowz
    !! 2dfloz overlap variables
    REAL(mcp), allocatable, dimension(:) :: P0P2_final_2dfloz_overlap,P0P2_obs_2dfloz_overlap,diff_2dfloz_overlap
    REAL(mcp) :: a_perp_2dfloz_overlap,a_par_2dfloz_overlap,k_val_overlap_2dfloz,z_2dfloz,sigv_2dfloz,b1_2dfloz,N_shot_2dfloz
    !! 2dfhiz overlap variables
    REAL(mcp), allocatable, dimension(:) :: P0P2_final_2dfhiz_overlap,P0P2_obs_2dfhiz_overlap,diff_2dfhiz_overlap
    REAL(mcp) :: a_perp_2dfhiz_overlap,a_par_2dfhiz_overlap,k_val_overlap_2dfhiz,z_2dfhiz,sigv_2dfhiz,b1_2dfhiz,N_shot_2dfhiz
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real(mcp) :: DataParams(:)
    sigma_v_cmass = DataParams(1)
    b1_cmass = DataParams(2)
    N_shot_cmass = DataParams(3)
    sig_v1_wigz = DataParams(4)
    sig_v2_wigz = DataParams(5)
    b_lowz_wigz = DataParams(6)
    b_highz_wigz = DataParams(7)
    sigv_lowz = DataParams(8)
    b1_lowz = DataParams(9)
    N_shot_lowz = DataParams(10)
    b1_2dfloz = DataParams(11)
    b1_2dfhiz = DataParams(12)
    sigv_2dfloz = DataParams(13)
    sigv_2dfhiz = DataParams(14)
    N_shot_2dfloz = DataParams(15)
    N_shot_2dfhiz = DataParams(16)
    ampiayes = DataParams(17)
    redziayes = DataParams(18)
    lumiayes = DataParams(19)
    a1phot = DataParams(20)
    a2phot = DataParams(21)
    a3phot = DataParams(22)
    a4phot = DataParams(23)
    a5phot = DataParams(24)
    a6phot = DataParams(25)
    a7phot = DataParams(26)
    meadnuis1 = DataParams(27)
    meadnuis2 = DataParams(28)


    CosmoJBD_LnLike = 0.0d0
    startrsd  = OMP_get_wtime()

    bbiaz0yes = b1_lowz
    bbiaz1yes = b1_cmass
    bbiaz2yes = b1_2dfloz
    bbiaz3yes = b1_2dfhiz

    obj%setscenarioomp = this%set_scenario

    k_grid_num_cmass  = this%cmass_mp%k_num

    if(this%use_cmass_full) then
        allocate(P0_theory_win_NGC(k_grid_num_cmass))
        allocate(P2_theory_win_NGC(k_grid_num_cmass))
        allocate(P0_theory_win_SGC(k_grid_num_cmass))
        allocate(P2_theory_win_SGC(k_grid_num_cmass))
        allocate(diff_NGC(this%cmass_mp%size_cov))
        allocate(diff_SGC(this%cmass_mp%size_cov))
    end if

    if(this%use_wigglez_full) then
        allocate(P0P2P4_lowz_final_wigz(this%wigglez_mp%size_cov,this%wigglez_mp%num_regions))
        allocate(P0P2P4_highz_final_wigz(this%wigglez_mp%size_cov,this%wigglez_mp%num_regions))
        allocate(P0P2P4_lowz_obs_wigz(this%wigglez_mp%size_cov,this%wigglez_mp%num_regions))
        allocate(P0P2P4_highz_obs_wigz(this%wigglez_mp%size_cov,this%wigglez_mp%num_regions))
        allocate(diff_lowz_wigz(this%wigglez_mp%size_cov,this%wigglez_mp%num_regions))
        allocate(diff_highz_wigz(this%wigglez_mp%size_cov,this%wigglez_mp%num_regions))
        allocate(like_lowZ(this%wigglez_mp%num_regions))
        allocate(like_highZ(this%wigglez_mp%num_regions))
    end if

    if(this%use_cmass_overlap) then
        allocate(P0P2_final_cmass_overlap(2*this%cmass_mp_overlap%k_num_obs)) !! Only want P0,P2 no measurement of P4
        allocate(diff_cmass_overlap(2*this%cmass_mp_overlap%k_num_obs))
        end if

        if(this%use_lowz_overlap) then
        allocate(P0P2_final_lowz_overlap(2*this%lowz_mp_overlap%k_num_obs)) !! Only want P0,P2 no measurement of P4
        allocate(diff_lowz_overlap(2*this%lowz_mp_overlap%k_num_obs))
        end if

        if(this%use_2dfloz_overlap) then
        allocate(P0P2_final_2dfloz_overlap(2*this%twodfloz_mp_overlap%k_num_obs)) !! Only want P0,P2 no measurement of P4
        allocate(diff_2dfloz_overlap(2*this%twodfloz_mp_overlap%k_num_obs))
        end if

        if(this%use_2dfhiz_overlap) then
        allocate(P0P2_final_2dfhiz_overlap(2*this%twodfhiz_mp_overlap%k_num_obs)) !! Only want P0,P2 no measurement of P4
        allocate(diff_2dfhiz_overlap(2*this%twodfhiz_mp_overlap%k_num_obs))
        end if

    !-----------------------------------------------------------------------------------------
    ! Setup nuisances parameters -- this will need to be changed..
    !-----------------------------------------------------------------------------------------

    like_NGC = 0.0d0
    like_SGC = 0.0d0
    Mp_WigZ_LnLike = 0.0d0

    scale_factor_L0 = 2.00
    scale_factor_L2 = 2.0/5.0
    scale_factor_L4 = 2.0/9.0

    !-----------------------------------------------------------------------------------------
    ! Compute theory -- Only calculated if needed.
    !-----------------------------------------------------------------------------------------

        z_cmass = this%cmass_mp%z_eff
        z_lowz_wigz = this%wigglez_mp%z_eff_low
        z_highz_wigz = this%wigglez_mp%z_eff_high
        z_lowz = this%lowz_mp_overlap%z_eff
        z_2dfloz = this%twodfloz_mp_overlap%z_eff
        z_2dfhiz = this%twodfhiz_mp_overlap%z_eff


        if(this%use_cmass_full) call Get_CMASS_TheoryMPK(P0_theory_win_NGC,P2_theory_win_NGC,P0_theory_win_SGC,P2_theory_win_SGC)

        if(this%use_wigglez_full) call Get_WiggleZ_TheoryMPK(P0P2P4_lowz_final_wigz,P0P2P4_highz_final_wigz)

        if(this%use_cmass_overlap) call Get_CMASS_overlap_TheoryMPK(P0P2_final_cmass_overlap)
        if(this%use_lowz_overlap) call Get_LOWZ_overlap_TheoryMPK(P0P2_final_lowz_overlap)
        if(this%use_2dfloz_overlap) call Get_2dfloz_overlap_TheoryMPK(P0P2_final_2dfloz_overlap)
        if(this%use_2dfhiz_overlap) call Get_2dfhiz_overlap_TheoryMPK(P0P2_final_2dfhiz_overlap)

        !----------------------------------------------------------------------------------------------------------
        ! Compute likelihood -- CMASS full
        !----------------------------------------------------------------------------------------------------------

        if(this%use_cmass_full) then
                do inum=1,this%cmass_mp%size_cov
                        if(inum .LE. this%cmass_mp%k_num) then

                                diff_NGC(inum) = this%cmass_mp%P0_multipole_NGC(inum) - P0_theory_win_NGC(inum)
                                diff_SGC(inum) = this%cmass_mp%P0_multipole_SGC(inum) - P0_theory_win_SGC(inum)
                        else
                            	temp_num = inum - this%cmass_mp%k_num

                                diff_NGC(inum) = this%cmass_mp%P2_multipole_NGC(temp_num) - P2_theory_win_NGC(temp_num)
                                diff_SGC(inum) = this%cmass_mp%P2_multipole_SGC(temp_num) - P2_theory_win_SGC(temp_num)

                        endif
                end do
                like_NGC = Matrix_QuadForm(INV_Cov_matrix_NGC,diff_NGC)/2.0
                like_SGC = Matrix_QuadForm(INV_Cov_matrix_SGC,diff_SGC)/2.0
        end if

        !----------------------------------------------------------------------------------------------------------
        ! Compute likelihood -- WiggleZ  (sum of 6 regions..)
        !----------------------------------------------------------------------------------------------------------

        if(this%use_wigglez_full) then
        DO region_num = 1,this%wigglez_mp%num_regions
                diff_lowz_wigz(:,region_num) =  P0P2P4_lowz_obs_wigz(:,region_num) - P0P2P4_lowz_final_wigz(:,region_num)
                diff_highz_wigz(:,region_num) = P0P2P4_highz_obs_wigz(:,region_num) - P0P2P4_highz_final_wigz(:,region_num)
        ENDDO

	DO region_num = 1,this%wigglez_mp%num_regions
                like_lowZ(region_num) = (Matrix_QuadForm(this%wigglez_mp%INV_Cov_matrix_lowZ(:,:,region_num),diff_lowz_wigz(:,region_num)))/2.0
                like_highZ(region_num) = (Matrix_QuadForm(this%wigglez_mp%INV_Cov_matrix_highZ(:,:,region_num),diff_highz_wigz(:,region_num)))/2.0
        ENDDO

	Mp_WigZ_LnLike = SUM(like_lowZ) + SUM(like_highZ)

        end if
	!----------------------------------------------------------------------------------------------------------
        ! Compute Delta = Obs - Theory data vector for CMASS overlap, LOWZ overlap, 2dfloz/2dfhiz overlap
        !----------------------------------------------------------------------------------------------------------

        if((obj%setscenarioomp == 3) .and. this%use_cmass_overlap) then
!               diff_cmass_overlap(:) = P0P2_obs_cmass_overlap(:) - P0P2_final_cmass_overlap(:)
                diff_cmass_overlap(:) = P0P2_final_cmass_overlap(:)
        end if
	if((obj%setscenarioomp == 3) .and. this%use_lowz_overlap) then
!               diff_lowz_overlap(:) = P0P2_obs_lowz_overlap(:) - P0P2_final_lowz_overlap(:)
                diff_lowz_overlap(:) = P0P2_final_lowz_overlap(:)
        end if
	if((obj%setscenarioomp == 3) .and. this%use_2dfloz_overlap) then
                diff_2dfloz_overlap(:) = P0P2_final_2dfloz_overlap(:)
        end if
	if((obj%setscenarioomp == 3) .and. this%use_2dfhiz_overlap) then
                diff_2dfhiz_overlap(:) = P0P2_final_2dfhiz_overlap(:)
        end if

        if(obj%setscenarioomp < 3) then
                CosmoJBD_LnLike = like_NGC + like_SGC + Mp_WigZ_LnLike
        end if


    use_morellz = this%use_morell
    use_rombintz = this%use_rombint
    nangbins = 9
    nzbins = 4
    if(use_morellz == .false.) nellbins = 31
    if(use_morellz == .true.) nellbins = 101
    obj%nellbinsomp = nellbins
    nell = 58999
    nzp = 70
    wtrapmax = 37
    wtrapmaxlens = 17 !assuming wtrapmaxlens < wtrapmax
    obj%wtrapmaxomp = wtrapmax
    obj%wtrapmaxlensomp = wtrapmaxlens

    allocate(obj%theoryomp)
    allocate(obj%ellarromp(nellbins))
    allocate(obj%exact_z_index(this%num_z))
    allocate(obj%momp(2))
    allocate(obj%psourcetypeomp(nzbins))
    allocate(obj%plenstypeomplowz)
    allocate(obj%plenstypeompcmass)
    allocate(obj%plenstypeomp2dfloz)
    allocate(obj%plenstypeomp2dfhiz)
    allocate(obj%aphotarromp(nzbins))
    allocate(obj%lumarromp(nzbins))
    allocate(obj%exact_zmodomp(wtrapmax))
    allocate(obj%exact_zmodomplenslowz(wtrapmaxlens))
    allocate(obj%exact_zmodomplenscmass(wtrapmaxlens))
    allocate(obj%weightarromp(nzbins,wtrapmax-1))
    allocate(obj%weightarromplowz(nzbins,wtrapmaxlens-1)) !same as weightarromp just evaluated at lowz redshifts
    allocate(obj%weightarrompcmass(nzbins,wtrapmaxlens-1)) !same as weightarromp just evaluated at cmass redshifts
    allocate(obj%weightarromp2dfloz(nzbins,wtrapmaxlens-1)) !same as weightarromp just evaluated at lowz redshifts
    allocate(obj%weightarromp2dfhiz(nzbins,wtrapmaxlens-1)) !same as weightarromp just evaluated at cmass redshifts
    allocate(obj%weightpsarromplowz(nzbins,wtrapmaxlens-1)) !same as weightarromp just evaluated at lowz redshifts
    allocate(obj%weightpsarrompcmass(nzbins,wtrapmaxlens-1)) !same as weightarromp just evaluated at cmass redshifts
    allocate(obj%weightpsarromp2dfloz(nzbins,wtrapmaxlens-1)) !same as weightarromp just evaluated at lowz redshifts
    allocate(obj%weightpsarromp2dfhiz(nzbins,wtrapmaxlens-1)) !same as weightarromp just evaluated at cmass redshifts
    allocate(obj%weightpsarromp(nzbins,wtrapmax-1))
    allocate(obj%weightarromplenslowz(wtrapmaxlens-1))
    allocate(obj%weightarromplenscmass(wtrapmaxlens-1))
    allocate(obj%weightarromplens2dfloz(wtrapmaxlens-1))
    allocate(obj%weightarromplens2dfhiz(wtrapmaxlens-1))
    allocate(obj%m1arromp(nzbins*(nzbins+1)/2))
    allocate(obj%m2arromp(nzbins*(nzbins+1)/2))
    allocate(obj%weightnoarromp(nellbins,wtrapmax-1))
    allocate(obj%weightnoarromplenslowz(nellbins,wtrapmaxlens-1))
    allocate(obj%weightnoarromplenscmass(nellbins,wtrapmaxlens-1))
    allocate(obj%weightnoarromplens2dfloz(nellbins,wtrapmaxlens-1))
    allocate(obj%weightnoarromplens2dfhiz(nellbins,wtrapmaxlens-1))
    allocate(obj%clinterptypeomp)
    allocate(obj%clinterptypeiiomp)
    allocate(obj%clinterptypegiomp)
    allocate(obj%clinterptypecrosslowzomp)
    allocate(obj%clinterptypecrosscmassomp)
    allocate(obj%clinterptypecross2dflozomp)
    allocate(obj%clinterptypecross2dfhizomp)
    allocate(obj%clinterptypecrosslowzgiomp)
    allocate(obj%clinterptypecrosscmassgiomp)
    allocate(obj%clinterptypecross2dflozgiomp)
    allocate(obj%clinterptypecross2dfhizgiomp)
    allocate(obj%gcubicomp)

    allocate(clarrbig(nell))
    allocate(clarrbiggi(nell))
    allocate(clarrbigii(nell))
    allocate(clarrbigcrosslowz(nell))
    allocate(clarrbigcrosscmass(nell))
    allocate(clarrbigcross2dfloz(nell))
    allocate(clarrbigcross2dfhiz(nell))
    allocate(clarrbigcrosslowzgi(nell))
    allocate(clarrbigcrosscmassgi(nell))
    allocate(clarrbigcross2dflozgi(nell))
    allocate(clarrbigcross2dfhizgi(nell))
    allocate(ellgentestarr(nell))
    allocate(xiplusminus(nzbins*(nzbins+1)/2*nangbins*2))
    allocate(xiplusminusgg(nzbins*(nzbins+1)/2*nangbins*2))
    allocate(xiplusminusii(nzbins*(nzbins+1)/2*nangbins*2))
    allocate(xiplusminusgi(nzbins*(nzbins+1)/2*nangbins*2))
    allocate(dataminustheory(this%sizcov))
    allocate(dataminustheoryhihi(this%sizcov))
    allocate(dataminustheoryhoho(this%sizcov))
    allocate(finaltheory(this%sizcovpremask))
    allocate(xicrossarrlowz(nzbins*nangbins))
    allocate(xicrossarrlowzgi(nzbins*nangbins))
    allocate(xicrossarrcmass(nzbins*nangbins))
    allocate(xicrossarrcmassgi(nzbins*nangbins))
    allocate(xicrossarr2dfloz(nzbins*nangbins))
    allocate(xicrossarr2dflozgi(nzbins*nangbins))
    allocate(xicrossarr2dfhiz(nzbins*nangbins))
    allocate(xicrossarr2dfhizgi(nzbins*nangbins))
    allocate(clarr(nellbins,nzbins*(nzbins+1)/2))
    allocate(clarrii(nellbins,nzbins*(nzbins+1)/2))
    allocate(clarrgi(nellbins,nzbins*(nzbins+1)/2))
    allocate(clarrcrosslowz(nellbins,nzbins))
    allocate(clarrcrosscmass(nellbins,nzbins))
    allocate(clarrcross2dfloz(nellbins,nzbins))
    allocate(clarrcross2dfhiz(nellbins,nzbins))
    allocate(clarrcrosslowzgi(nellbins,nzbins))
    allocate(clarrcrosscmassgi(nellbins,nzbins))
    allocate(clarrcross2dflozgi(nellbins,nzbins))
    allocate(clarrcross2dfhizgi(nellbins,nzbins))
    allocate(ellarr(nellbins))
    allocate(xiplus(nangbins))
    allocate(xiplusii(nangbins))
    allocate(xiplusgi(nangbins))
    allocate(ximinus(nangbins))
    allocate(ximinusii(nangbins))
    allocate(ximinusgi(nangbins))
    allocate(xicrosslowz(nangbins))
    allocate(xicrosslowzgi(nangbins))
    allocate(xicrosscmass(nangbins))
    allocate(xicrosscmassgi(nangbins))
    allocate(xicross2dfloz(nangbins))
    allocate(xicross2dflozgi(nangbins))
    allocate(xicross2dfhiz(nangbins))
    allocate(xicross2dfhizgi(nangbins))
    allocate(mumax1arr(nzbins*(nzbins+1)/2))
    allocate(mumax2arr(nzbins*(nzbins+1)/2))
    allocate(aphotarr(nzbins))
    allocate(rcslum(nzbins))
    allocate(cfhtlum(nzbins))
    allocate(lumarr(nzbins))
    allocate(m1arr(nzbins*(nzbins+1)/2))
    allocate(m2arr(nzbins*(nzbins+1)/2))
    allocate(mnomp(2))
    allocate(weightarr(nzbins,wtrapmax-1))
    allocate(weightarrlowz(nzbins,wtrapmaxlens-1)) !same as weightarr just evaluated at lowz redshifts
    allocate(weightarrcmass(nzbins,wtrapmaxlens-1)) !same as weightarr just evaluated at cmass redshifts
    allocate(weightarr2dfloz(nzbins,wtrapmaxlens-1)) !same as weightarr just evaluated at lowz redshifts
    allocate(weightarr2dfhiz(nzbins,wtrapmaxlens-1)) !same as weightarr just evaluated at cmass redshifts
    allocate(weightpsarrlowz(nzbins,wtrapmaxlens-1)) !same as weightarr just evaluated at lowz redshifts
    allocate(weightpsarrcmass(nzbins,wtrapmaxlens-1)) !same as weightarr just evaluated at cmass redshifts
    allocate(weightpsarr2dfloz(nzbins,wtrapmaxlens-1)) !same as weightarr just evaluated at lowz redshifts
    allocate(weightpsarr2dfhiz(nzbins,wtrapmaxlens-1)) !same as weightarr just evaluated at cmass redshifts
    allocate(weightpsarr(nzbins,wtrapmax-1))
    allocate(weightarrlenslowz(wtrapmaxlens-1))
    allocate(weightarrlenscmass(wtrapmaxlens-1))
    allocate(weightarrlens2dfloz(wtrapmaxlens-1))
    allocate(weightarrlens2dfhiz(wtrapmaxlens-1))
    allocate(binwarr(nzbins))
    allocate(trapzarr(wtrapmax-2))
    allocate(trapzarrii(wtrapmax-2))
    allocate(trapzarrgi(wtrapmax-2))
    allocate(trapzarrcrosslowz(wtrapmaxlens-2))
    allocate(trapzarrcrosscmass(wtrapmaxlens-2))
    allocate(trapzarrcross2dfloz(wtrapmaxlens-2))
    allocate(trapzarrcross2dfhiz(wtrapmaxlens-2))
    allocate(trapzarrcrosslowzgi(wtrapmaxlens-2))
    allocate(trapzarrcrosscmassgi(wtrapmaxlens-2))
    allocate(trapzarrcross2dflozgi(wtrapmaxlens-2))
    allocate(trapzarrcross2dfhizgi(wtrapmaxlens-2))
    allocate(weightnoarr(nellbins,wtrapmax-1))
    allocate(weightnoarrlenslowz(nellbins,wtrapmaxlens-1))
    allocate(weightnoarrlenscmass(nellbins,wtrapmaxlens-1))
    allocate(weightnoarrlens2dfloz(nellbins,wtrapmaxlens-1))
    allocate(weightnoarrlens2dfhiz(nellbins,wtrapmaxlens-1))
    allocate(garr(Theory%MPK%ny))

    aphotarr(1) = a1phot
    aphotarr(2) = a2phot
    aphotarr(3) = a3phot
    aphotarr(4) = a4phot
    aphotarr(5) = a5phot
    obj%aphotarromp(1) = a1phot
    obj%aphotarromp(2) = a2phot
    obj%aphotarromp(3) = a3phot
    obj%aphotarromp(4) = a4phot
    obj%aphotarromp(5) = a5phot

    start  = OMP_get_wtime()
    cfhtlum = (/ 0.0174d0, 0.0694d0, 0.1518d0, 0.2182d0 /) !not used for KiDS
    if(use_cfht == .true.) lumarr = cfhtlum
    obj%lumarromp = lumarr

    !EFT
    obj%eftflag = CosmoSettings%EFTFlag
    obj%designereftmodel = CosmoSettings%DesignerEFTmodel

    obj%exact_zmodomp = this%exact_z
    tol_error = 0.005d0 !Numerical integration relative error, for Romberg
    pisj = 3.14159265359d0
    obj%pisjomp = 3.14159265359d0
    obj%h0omp = CMB%H0
    obj%homp = CMB%h
    !Modified for EFT start
    if (obj%eftflag == 2 .and. (obj%designereftmodel == 4 .or. obj%designereftmodel == 5 )) then
    obj%omdmomp = CMB%omdm/(CMB%h**2.0)
    obj%ombomp = CMB%omb/(CMB%h**2.0)
    obj%omkomp = CMB%omk/(CMB%h**2.0)    
    else
    obj%omdmomp = CMB%omdm
    obj%ombomp = CMB%omb
    obj%omkomp = CMB%omk
    end if
    !Modified for EFT end
    obj%theoryomp = Theory




!trapezoid z-steps for lowz
obj%exact_zmodomplenslowz(1) = 0.0d0
obj%exact_zmodomplenslowz(2) = 0.15d0
    do zim=3,17
     obj%exact_zmodomplenslowz(zim) = 0.15d0 + (0.43d0-0.15d0)/15.0d0*(zim-2)
    end do
obj%exact_zmodomplenslowz(17) = 0.43d0 !0.429999d0 !3.49d0

!trapezoid z-steps for cmass
obj%exact_zmodomplenscmass(1) = 0.0d0
obj%exact_zmodomplenscmass(2) = 0.43d0
    do zim=3,17
     obj%exact_zmodomplenscmass(zim) = 0.43 + (0.7d0-0.43d0)/15.0d0*(zim-2)
    end do
obj%exact_zmodomplenscmass(17) = 0.7d0 !0.699999d0 !3.49d0

    mumin=0.0d0
    mumin2v=0.025d0 !min integration redshift
    mumax=3.474999d0 !max integration redshift
    mumin2vcrosslowz=0.15 !0.150001d0 !min for lowz and 2dfloz
    mumaxcrosslowz=0.43d0 !0.4299999d0 !max for lowz and 2dfloz
    mumin2vcrosscmass=0.430001d0 !min for lowz and 2dfloz
    mumaxcrosscmass=0.7d0 !0.6999999d0 !max for lowz and 2dfloz

    !Generate array containing l-values where C(l) is evaluated
    if(use_morellz == .false.) morellfac = 0.5
    if(use_morellz == .true.) morellfac = 0.1
    do ellgen=1,nellbins
        if(ellgen < 10) ellarr(ellgen) = ellgen+1
        if(ellgen > 10 .or. ellgen == 10) ellarr(ellgen) = ellarr(ellgen-1)+morellfac*ellarr(ellgen-1)
    end do
    ellarr = int(ellarr)
    obj%ellarromp = ellarr
    obj%tol_erroromp = tol_error
    obj%mumin2vomp = mumin2v
    obj%mumaxomp = mumax
    obj%mumin2vcrossomplowz = mumin2vcrosslowz
    obj%mumaxcrossomplowz = mumaxcrosslowz
    obj%mumin2vcrossompcmass = mumin2vcrosscmass
    obj%mumaxcrossompcmass = mumaxcrosscmass

    obj%meadnuis1omp = meadnuis1 
    obj%meadnuis2omp = meadnuis2 
    obj%theoryomp%NL_MPK%z = Theory%NL_MPK%z

    obj%ampiayesomp = ampiayes
    obj%redziayesomp = redziayes
    obj%lumiayesomp = lumiayes
    obj%bbiaz0yesomp = bbiaz0yes
    obj%bbiaz1yesomp = bbiaz1yes
    obj%bbiaz2yesomp = bbiaz2yes
    obj%bbiaz3yesomp = bbiaz3yes
    obj%exact_z_index = this%exact_z_index
    wchoose = 1
    obj%wchooseomp = 1

    !interpolation for source redshift distributions
    do ho=1,nzbins
        call obj%psourcetypeomp(ho)%Init(this%arraysjfull(1,:,ho),this%arraysjfull(2,:,ho),nzp)
    end do

   !interpolation for lens redshift distribution
   call obj%plenstypeomplowz%Init(this%arraysjlenslowz(1,:),this%arraysjlenslowz(2,:),29)
   call obj%plenstypeompcmass%Init(this%arraysjlenscmass(1,:),this%arraysjlenscmass(2,:),28)
   call obj%plenstypeomp2dfloz%Init(this%arraysjlens2dfloz(1,:),this%arraysjlens2dfloz(2,:),29)
   call obj%plenstypeomp2dfhiz%Init(this%arraysjlens2dfhiz(1,:),this%arraysjlens2dfhiz(2,:),28)

    obj%kline = 0

    !tomographic ordering
    m1arr = (/ 1, 1, 1, 1, 2, 2, 2, 3, 3, 4 /)
    m2arr = (/ 1, 2, 3, 4, 2, 3, 4, 3, 4, 4 /)

    !compute integration range accounting for photo-z error
    mint = 0
    do m1=1,nzbins
        do m2=1,nzbins
            if(m2 >= m1) then
                mint = mint + 1
                mumax1arr(mint) = mumax + aphotarr(m1arr(mint))
                mumax2arr(mint) = mumax + aphotarr(m2arr(mint))
            end if
        end do
    end do

    !requires the redshift in CosmologyTypes.f90 file to be set to z>3.5 in case of photo-z varying
    do gigi=1,Theory%MPK%ny
        garr(gigi) = exp(rombint(sjgrowtha,1.0d0/(1.0d0+Theory%MPK%y(gigi)),1.0d0,tol_error))
    end do
    call obj%gcubicomp%Init(Theory%MPK%y,garr,Theory%MPK%ny)
    bbint = 0
    bbintigen = 0
    bbintmera = 0
    bbintmeraigen = 0


    if(use_rombintz == .true.) then !Romberg Integration for the outer integral

        !for nz tomographic bins, compute nz*(nz+1)/2 C(l)-combinations
        do bin = 1,nzbins*(nzbins+1)/2
            mint = 0
            do m1=1,nzbins
                do m2=1,nzbins
                    if(m2 >= m1) then 
                        mint = mint + 1
                        if(mint == bin) then
                            obj%m1omp=m1arr(mint)
                            obj%m2omp=m2arr(mint)
                            obj%momp(1)=m1arr(mint)
                            obj%momp(2)=m2arr(mint)
                        end if
                    end if
                end do
            end do
            mnomp = obj%momp
            bnum = bin
            obj%bnumomp = bin
            mumax1 = mumax1arr(bin)
            mumax2 = mumax2arr(bin)
            obj%mumax1omp = mumax1arr(bin)
            obj%mumax2omp = mumax2arr(bin)
            obj%mumaxlensomp = min(mumax1,mumax2)
    
            !!!##############COMPUTE LENSING POWER SPECTRA (GG, GI, II, Gg) AT DISTINCT L-VALUES. PARALLELIZATION EMPLOYED############
            !$OMP PARALLEL DO FIRSTPRIVATE(ellgen,obj)
                do ellgen=1,nellbins
                    obj%ellgenomp = ellgen
                    clarr(ellgen,obj%bnumomp) = rombint_obj(obj,sjclsobjsf,1.0d0/(1.0d0+obj%mumaxlensomp),1.0d0/(1.0d0+obj%mumin2vomp),obj%tol_erroromp) !!good!!
                    clarrii(ellgen,obj%bnumomp) = rombint_obj(obj,sjclsiiobjsf,1.0d0/(1.0d0+obj%mumaxlensomp),1.0d0/(1.0d0+obj%mumin2vomp),obj%tol_erroromp) !good!!
                    clarrgi(ellgen,obj%bnumomp) = rombint_obj(obj,sjclsgiobjsf,1.0d0/(1.0d0+obj%mumaxlensomp),1.0d0/(1.0d0+obj%mumin2vomp),obj%tol_erroromp) !!good!!
                end do
            !$OMP END PARALLEL DO

            if((m1arr(bin) == m2arr(bin)) .and. (obj%setscenarioomp > 1)) then
                bbint = bbint+1
                obj%bnumomplens = bbint
                mumaxw = mumax1arr(bin)
                obj%mumaxwomp = mumaxw

                !$OMP PARALLEL DO FIRSTPRIVATE(ellgen,obj)
                    do ellgen=1,nellbins
                        obj%ellgenomp = ellgen
                        obj%uselensomp = 0
                        obj%bbiazconstomp = obj%bbiaz0yesomp
                        clarrcrosslowz(ellgen,obj%bnumomplens) = rombint_obj(obj,sjclscrossobjsf,1.0d0/(1.0d0+obj%mumaxcrossomplowz),1.0d0/(1.0d0+obj%mumin2vcrossomplowz),obj%tol_erroromp)
                        clarrcrosslowzgi(ellgen,obj%bnumomplens) = rombint_obj(obj,sjclscrossgiobjsf,1.0d0/(1.0d0+obj%mumaxcrossomplowz),1.0d0/(1.0d0+obj%mumin2vcrossomplowz),obj%tol_erroromp)
                        !!!clarrcrosslowz(ellgen,obj%bnumomplens) = rombint_obj(obj,sjclscrossobj,obj%mumin2vcrossomplowz,obj%mumaxcrossomplowz,obj%tol_erroromp) !with z instead of sf
                        obj%uselensomp = 1
                        obj%bbiazconstomp = obj%bbiaz1yesomp
                        clarrcrosscmass(ellgen,obj%bnumomplens) = rombint_obj(obj,sjclscrossobjsf,1.0d0/(1.0d0+obj%mumaxcrossompcmass),1.0d0/(1.0d0+obj%mumin2vcrossompcmass),obj%tol_erroromp)
                        clarrcrosscmassgi(ellgen,obj%bnumomplens) = rombint_obj(obj,sjclscrossgiobjsf,1.0d0/(1.0d0+obj%mumaxcrossompcmass),1.0d0/(1.0d0+obj%mumin2vcrossompcmass),obj%tol_erroromp)
                        !!!clarrcrosscmass(ellgen,obj%bnumomplens) = rombint_obj(obj,sjclscrossobj,obj%mumin2vcrossompcmass,obj%mumaxcrossompcmass,obj%tol_erroromp) !with z instead of sf
                        obj%uselensomp = 2
                        obj%bbiazconstomp = obj%bbiaz2yesomp
                        clarrcross2dfloz(ellgen,obj%bnumomplens) = rombint_obj(obj,sjclscrossobjsf,1.0d0/(1.0d0+obj%mumaxcrossomplowz),1.0d0/(1.0d0+obj%mumin2vcrossomplowz),obj%tol_erroromp)
                        clarrcross2dflozgi(ellgen,obj%bnumomplens) = rombint_obj(obj,sjclscrossgiobjsf,1.0d0/(1.0d0+obj%mumaxcrossomplowz),1.0d0/(1.0d0+obj%mumin2vcrossomplowz),obj%tol_erroromp)
                        obj%uselensomp = 3
                        obj%bbiazconstomp = obj%bbiaz3yesomp
                        clarrcross2dfhiz(ellgen,obj%bnumomplens) = rombint_obj(obj,sjclscrossobjsf,1.0d0/(1.0d0+obj%mumaxcrossompcmass),1.0d0/(1.0d0+obj%mumin2vcrossompcmass),obj%tol_erroromp)
                        clarrcross2dfhizgi(ellgen,obj%bnumomplens) = rombint_obj(obj,sjclscrossgiobjsf,1.0d0/(1.0d0+obj%mumaxcrossompcmass),1.0d0/(1.0d0+obj%mumin2vcrossompcmass),obj%tol_erroromp)
                    end do
                !$OMP END PARALLEL DO
            end if
        end do
    end if 

    !Trapezoid Integration for the outer integral!!! for the case use_rombintz=F !!! ------ do not use this option when varying photo-z params, in that case set use_rombintz = T
    if(use_rombintz == .false.) then

        binwarr = (/ 1, 5, 8, 10 /) !4 tomographic bins, m1arr(i) = m2arr(i) for these indices --- need to change this array for different number of tomographic bins
        weightarr(:,:) = 0.0d0
        weightarrlowz(:,:) = 0.0d0
        weightarrcmass(:,:) = 0.0d0
        weightpsarrlowz(:,:) = 0.0d0
        weightpsarrcmass(:,:) = 0.0d0
        weightarrlenslowz(:) = 0.0d0
        weightarrlenscmass(:) = 0.0d0
        weightarr2dfloz(:,:) = 0.0d0
        weightarr2dfhiz(:,:) = 0.0d0
        weightpsarr2dfloz(:,:) = 0.0d0
        weightpsarr2dfhiz(:,:) = 0.0d0
        weightarrlens2dfloz(:) = 0.0d0
        weightarrlens2dfhiz(:) = 0.0d0
        weightpsarr(:,:) = 0.0d0
        weightnoarr(:,:) = 0.0d0
        weightnoarrlenslowz(:,:) = 0.0d0
        weightnoarrlenscmass(:,:) = 0.0d0
        weightnoarrlens2dfloz(:,:) = 0.0d0
        weightnoarrlens2dfhiz(:,:) = 0.0d0

        do binw = 1,nzbins
            mumaxw = mumax1arr(binwarr(binw))
            obj%mumaxwomp = mumaxw
            bnum = binwarr(binw)
            obj%bnumomp = bnum
            obj%m1omp=m1arr(binwarr(binw))
            obj%m2omp=m2arr(binwarr(binw))
            obj%momp(1)=m1arr(binwarr(binw))
            obj%momp(2)=m2arr(binwarr(binw))
            mnomp = obj%momp
            obj%binwomp = binw

            !$OMP PARALLEL DO FIRSTPRIVATE(wttt,obj)
                do wttt=2,obj%wtrapmaxomp
                    weightarr(obj%binwomp,wttt-1) = sjclsobjonlyweight(obj,obj%exact_zmodomp(wttt))   !have to make sure mumaxw > obj%exact_zmodomp(wttt)
                    weightpsarr(obj%binwomp,wttt-1) = sjclsiiobjonlyweight(obj,obj%exact_zmodomp(wttt))   !have to make sure mumaxw > obj%exact_zmodomp(wttt)
                    if((wttt < (obj%wtrapmaxlensomp + 1)) .and. (obj%setscenarioomp > 1)) then
                        weightarrlowz(obj%binwomp,wttt-1) = sjclsobjonlyweight(obj,obj%exact_zmodomplenslowz(wttt))   !have to make sure mumaxw > obj%exact_zmodomp(wttt)
                        weightarrcmass(obj%binwomp,wttt-1) = sjclsobjonlyweight(obj,obj%exact_zmodomplenscmass(wttt))   !have to make sure mumaxw > obj%exact_zmodomp(wttt)
                        weightarr2dfloz(obj%binwomp,wttt-1) = weightarrlowz(obj%binwomp,wttt-1) !since same redshifts
                        weightarr2dfhiz(obj%binwomp,wttt-1) = weightarrcmass(obj%binwomp,wttt-1) !since same redshifts
                        weightpsarrlowz(obj%binwomp,wttt-1) = sjclsiiobjonlyweight(obj,obj%exact_zmodomplenslowz(wttt))   !have to make sure mumaxw > obj%exact_zmodomp(wttt)
                        weightpsarrcmass(obj%binwomp,wttt-1) = sjclsiiobjonlyweight(obj,obj%exact_zmodomplenscmass(wttt))   !have to make sure mumaxw > obj%exact_zmodomp(wttt)
                        weightpsarr2dfloz(obj%binwomp,wttt-1) = weightpsarrlowz(obj%binwomp,wttt-1) !since same redshifts
                        weightpsarr2dfhiz(obj%binwomp,wttt-1) = weightpsarrcmass(obj%binwomp,wttt-1) !since same redshifts
                    end if
                    if((obj%binwomp == 1) .and. (wttt < (obj%wtrapmaxlensomp+1)) .and. (obj%setscenarioomp > 1)) then
                        obj%uselensomp = 0
                        obj%bbiazconstomp = obj%bbiaz0yesomp
                        weightarrlenslowz(wttt-1) = sjclscrossobjonlyweight(obj,obj%exact_zmodomplenslowz(wttt))   !have to make sure mumaxw > obj%exact_zmodomp(wttt)
                        obj%uselensomp = 1
                        obj%bbiazconstomp = obj%bbiaz1yesomp
                        weightarrlenscmass(wttt-1) = sjclscrossobjonlyweight(obj,obj%exact_zmodomplenscmass(wttt))   !have to make sure mumaxw > obj%exact_zmodomp(wttt)
                        obj%uselensomp = 2
                        obj%bbiazconstomp = obj%bbiaz2yesomp
                        weightarrlens2dfloz(wttt-1) = sjclscrossobjonlyweight(obj,obj%exact_zmodomplenslowz(wttt))   !have to make sure mumaxw > obj%exact_zmodomp(wttt)
                        obj%uselensomp = 3
                        obj%bbiazconstomp = obj%bbiaz3yesomp
                        weightarrlens2dfhiz(wttt-1) = sjclscrossobjonlyweight(obj,obj%exact_zmodomplenscmass(wttt))   !have to make sure mumaxw > obj%exact_zmodomp(wttt)
                    end if
                end do
            !$OMP END PARALLEL DO
        end do

        do ellgen=1,nellbins
            obj%ellgenomp = ellgen
            !$OMP PARALLEL DO FIRSTPRIVATE(wttt,obj)
                do wttt=2,obj%wtrapmaxomp
                    weightnoarr(obj%ellgenomp,wttt-1) = sjclsobjnoweight(obj,obj%exact_zmodomp(wttt))   !have to make sure mumaxw > this%exact_z(wttt)
                    if((wttt < (obj%wtrapmaxlensomp+1)) .and. (obj%setscenarioomp > 1)) then
                        weightnoarrlenslowz(obj%ellgenomp,wttt-1) = sjclscrossobjnoweight(obj,obj%exact_zmodomplenslowz(wttt))   !have to make sure mumaxw > this%exact_z(wttt)
                        weightnoarrlenscmass(obj%ellgenomp,wttt-1) = sjclscrossobjnoweight(obj,obj%exact_zmodomplenscmass(wttt))   !have to make sure mumaxw > this%exact_z(wttt)
                        weightnoarrlens2dfloz(obj%ellgenomp,wttt-1) = weightnoarrlenslowz(obj%ellgenomp,wttt-1)
                        weightnoarrlens2dfhiz(obj%ellgenomp,wttt-1) = weightnoarrlenscmass(obj%ellgenomp,wttt-1)
                    end if
                end do
            !$OMP END PARALLEL DO
        end do

        !for nz tomographic bins, compute nz*(nz+1)/2 C(l)-combinations
        do bin = 1,nzbins*(nzbins+1)/2
            if(m1arr(bin) == m2arr(bin)) then
                bbintigen = bbintigen+1
                obj%bnumomplens = bbintigen
            end if
            mint = 0
            do m1=1,nzbins
                do m2=1,nzbins
                    if(m2 >= m1) then 
                        mint = mint + 1
                        if(mint == bin) then
                            obj%m1omp=m1arr(mint)
                            obj%m2omp=m2arr(mint)
                            obj%momp(1)=m1arr(mint)
                            obj%momp(2)=m2arr(mint)
                        end if
                    end if
                end do
            end do
            mnomp = obj%momp
            bnum = bin
            obj%bnumomp = bin
            mumax1 = mumax1arr(bin)
            mumax2 = mumax2arr(bin)
            obj%mumax1omp = mumax1arr(bin)
            obj%mumax2omp = mumax2arr(bin)
            obj%mumaxlensomp = min(mumax1,mumax2)
            obj%weightarromp = weightarr
            obj%weightarromplowz = weightarrlowz
            obj%weightarrompcmass = weightarrcmass
            obj%weightpsarromplowz = weightpsarrlowz
            obj%weightpsarrompcmass = weightpsarrcmass
            obj%weightarromplenslowz = weightarrlenslowz
            obj%weightarromplenscmass = weightarrlenscmass
            obj%weightarromp2dfloz = weightarr2dfloz
            obj%weightarromp2dfhiz = weightarr2dfhiz
            obj%weightarromplens2dfloz = weightarrlens2dfloz
            obj%weightarromplens2dfhiz = weightarrlens2dfhiz
            obj%weightpsarromp2dfloz = weightpsarr2dfloz
            obj%weightpsarromp2dfhiz = weightpsarr2dfhiz
            obj%weightpsarromp = weightpsarr
            obj%weightnoarromp = weightnoarr
            obj%weightnoarromplenslowz = weightnoarrlenslowz
            obj%weightnoarromplenscmass = weightnoarrlenscmass
            obj%weightnoarromplens2dfloz = weightnoarrlens2dfloz
            obj%weightnoarromplens2dfhiz = weightnoarrlens2dfhiz
            obj%m1arromp = m1arr
            obj%m2arromp = m2arr

            !!COMPUTE LENSING POWER SPECTRA (GG, GI, II, Gg) AT DISTINCT L-VALUES. PARALLELIZATION EMPLOYED
            do ellgen=1,obj%nellbinsomp
                obj%ellgenomp = ellgen
                !$OMP PARALLEL DO FIRSTPRIVATE(ttt,obj)
                    do ttt=2,36
                        trapzarr(ttt-1) = 0.5d0*(obj%weightnoarromp(obj%ellgenomp,ttt-1)*obj%weightarromp(obj%m1arromp(obj%bnumomp),ttt-1)*obj%weightarromp(obj%m2arromp(obj%bnumomp),ttt-1) + obj%weightnoarromp(obj%ellgenomp,ttt)*obj%weightarromp(obj%m1arromp(obj%bnumomp),ttt)*obj%weightarromp(obj%m2arromp(obj%bnumomp),ttt))*(obj%exact_zmodomp(ttt+1)-obj%exact_zmodomp(ttt))
                        trapzarrii(ttt-1) = 0.5d0*(obj%weightnoarromp(obj%ellgenomp,ttt-1)*obj%weightpsarromp(obj%m1arromp(obj%bnumomp),ttt-1)*obj%weightpsarromp(obj%m2arromp(obj%bnumomp),ttt-1) + obj%weightnoarromp(obj%ellgenomp,ttt)*obj%weightpsarromp(obj%m1arromp(obj%bnumomp),ttt)*obj%weightpsarromp(obj%m2arromp(obj%bnumomp),ttt))*(obj%exact_zmodomp(ttt+1)-obj%exact_zmodomp(ttt))
                        trapzarrgi(ttt-1) = 0.5d0*(obj%weightnoarromp(obj%ellgenomp,ttt-1)*(obj%weightarromp(obj%m1arromp(obj%bnumomp),ttt-1)*obj%weightpsarromp(obj%m2arromp(obj%bnumomp),ttt-1)+obj%weightarromp(obj%m2arromp(obj%bnumomp),ttt-1)*obj%weightpsarromp(obj%m1arromp(obj%bnumomp),ttt-1)) + obj%weightnoarromp(obj%ellgenomp,ttt)*(obj%weightarromp(obj%m1arromp(obj%bnumomp),ttt)*obj%weightpsarromp(obj%m2arromp(obj%bnumomp),ttt)+obj%weightarromp(obj%m2arromp(obj%bnumomp),ttt)*obj%weightpsarromp(obj%m1arromp(obj%bnumomp),ttt)))*(obj%exact_zmodomp(ttt+1)-obj%exact_zmodomp(ttt))
                    end do
                !$OMP END PARALLEL DO
                clarr(ellgen,obj%bnumomp)=sum(trapzarr)
                clarrii(ellgen,obj%bnumomp)=sum(trapzarrii)
                clarrgi(ellgen,obj%bnumomp)=sum(trapzarrgi)
                if((m1arr(bin) == m2arr(bin)) .and. (obj%setscenarioomp > 1)) then
                    !$OMP PARALLEL DO FIRSTPRIVATE(ttt,obj)
                        do ttt=2,16
                            trapzarrcrosslowz(ttt-1) = 0.5d0*(obj%weightnoarromplenslowz(obj%ellgenomp,ttt-1)*obj%weightarromplowz(obj%m1arromp(obj%bnumomp),ttt-1)*obj%weightarromplenslowz(ttt-1) + obj%weightnoarromplenslowz(obj%ellgenomp,ttt)*obj%weightarromplowz(obj%m1arromp(obj%bnumomp),ttt)*obj%weightarromplenslowz(ttt))*(obj%exact_zmodomplenslowz(ttt+1)-obj%exact_zmodomplenslowz(ttt))
                            trapzarrcrosscmass(ttt-1) = 0.5d0*(obj%weightnoarromplenscmass(obj%ellgenomp,ttt-1)*obj%weightarrompcmass(obj%m1arromp(obj%bnumomp),ttt-1)*obj%weightarromplenscmass(ttt-1) + obj%weightnoarromplenscmass(obj%ellgenomp,ttt)*obj%weightarrompcmass(obj%m1arromp(obj%bnumomp),ttt)*obj%weightarromplenscmass(ttt))*(obj%exact_zmodomplenscmass(ttt+1)-obj%exact_zmodomplenscmass(ttt))
                            trapzarrcross2dfloz(ttt-1) = 0.5d0*(obj%weightnoarromplens2dfloz(obj%ellgenomp,ttt-1)*obj%weightarromp2dfloz(obj%m1arromp(obj%bnumomp),ttt-1)*obj%weightarromplens2dfloz(ttt-1) + obj%weightnoarromplens2dfloz(obj%ellgenomp,ttt)*obj%weightarromp2dfloz(obj%m1arromp(obj%bnumomp),ttt)*obj%weightarromplens2dfloz(ttt))*(obj%exact_zmodomplenslowz(ttt+1)-obj%exact_zmodomplenslowz(ttt))
                            trapzarrcross2dfhiz(ttt-1) = 0.5d0*(obj%weightnoarromplens2dfhiz(obj%ellgenomp,ttt-1)*obj%weightarromp2dfhiz(obj%m1arromp(obj%bnumomp),ttt-1)*obj%weightarromplens2dfhiz(ttt-1) + obj%weightnoarromplens2dfhiz(obj%ellgenomp,ttt)*obj%weightarromp2dfhiz(obj%m1arromp(obj%bnumomp),ttt)*obj%weightarromplens2dfhiz(ttt))*(obj%exact_zmodomplenscmass(ttt+1)-obj%exact_zmodomplenscmass(ttt))
                            trapzarrcrosslowzgi(ttt-1) = 0.5d0*(obj%weightnoarromplenslowz(obj%ellgenomp,ttt-1)*obj%weightpsarromplowz(obj%m1arromp(obj%bnumomp),ttt-1)*obj%weightarromplenslowz(ttt-1) + obj%weightnoarromplenslowz(obj%ellgenomp,ttt)*obj%weightpsarromplowz(obj%m1arromp(obj%bnumomp),ttt)*obj%weightarromplenslowz(ttt))*(obj%exact_zmodomplenslowz(ttt+1)-obj%exact_zmodomplenslowz(ttt))
                            trapzarrcrosscmassgi(ttt-1) = 0.5d0*(obj%weightnoarromplenscmass(obj%ellgenomp,ttt-1)*obj%weightpsarrompcmass(obj%m1arromp(obj%bnumomp),ttt-1)*obj%weightarromplenscmass(ttt-1) + obj%weightnoarromplenscmass(obj%ellgenomp,ttt)*obj%weightpsarrompcmass(obj%m1arromp(obj%bnumomp),ttt)*obj%weightarromplenscmass(ttt))*(obj%exact_zmodomplenscmass(ttt+1)-obj%exact_zmodomplenscmass(ttt))
                            trapzarrcross2dflozgi(ttt-1) = 0.5d0*(obj%weightnoarromplens2dfloz(obj%ellgenomp,ttt-1)*obj%weightpsarromp2dfloz(obj%m1arromp(obj%bnumomp),ttt-1)*obj%weightarromplens2dfloz(ttt-1) + obj%weightnoarromplens2dfloz(obj%ellgenomp,ttt)*obj%weightpsarromp2dfloz(obj%m1arromp(obj%bnumomp),ttt)*obj%weightarromplens2dfloz(ttt))*(obj%exact_zmodomplenslowz(ttt+1)-obj%exact_zmodomplenslowz(ttt))
                            trapzarrcross2dfhizgi(ttt-1) = 0.5d0*(obj%weightnoarromplens2dfhiz(obj%ellgenomp,ttt-1)*obj%weightpsarromp2dfhiz(obj%m1arromp(obj%bnumomp),ttt-1)*obj%weightarromplens2dfhiz(ttt-1) + obj%weightnoarromplens2dfhiz(obj%ellgenomp,ttt)*obj%weightpsarromp2dfhiz(obj%m1arromp(obj%bnumomp),ttt)*obj%weightarromplens2dfhiz(ttt))*(obj%exact_zmodomplenscmass(ttt+1)-obj%exact_zmodomplenscmass(ttt))
                        end do
                    !$OMP END PARALLEL DO
                    clarrcrosslowz(ellgen,obj%bnumomplens)=sum(trapzarrcrosslowz)
                    clarrcrosscmass(ellgen,obj%bnumomplens)=sum(trapzarrcrosscmass)
                    clarrcross2dfloz(ellgen,obj%bnumomplens)=sum(trapzarrcross2dfloz)
                    clarrcross2dfhiz(ellgen,obj%bnumomplens)=sum(trapzarrcross2dfhiz)
                    clarrcrosslowzgi(ellgen,obj%bnumomplens)=sum(trapzarrcrosslowzgi)
                    clarrcrosscmassgi(ellgen,obj%bnumomplens)=sum(trapzarrcrosscmassgi)
                    clarrcross2dflozgi(ellgen,obj%bnumomplens)=sum(trapzarrcross2dflozgi)
                    clarrcross2dfhizgi(ellgen,obj%bnumomplens)=sum(trapzarrcross2dfhizgi)
                end if
            end do
        end do
    end if

    !Interpolation for C(l)
    ellgentestarr = this%ellgentestarrini
    intnum=nellbins
    do biin = 1,nzbins*(nzbins+1)/2
        call obj%clinterptypeomp%Init(ellarr,clarr(:,biin),intnum)
        call obj%clinterptypeiiomp%Init(ellarr,clarrii(:,biin),intnum)
        call obj%clinterptypegiomp%Init(ellarr,clarrgi(:,biin),intnum)
        if((m1arr(biin) == m2arr(biin)) .and. (obj%setscenarioomp > 1)) then
            bbintmera = bbintmera+1
            call obj%clinterptypecrosslowzomp%Init(ellarr,clarrcrosslowz(:,bbintmera),intnum)
            call obj%clinterptypecrosscmassomp%Init(ellarr,clarrcrosscmass(:,bbintmera),intnum)
            call obj%clinterptypecross2dflozomp%Init(ellarr,clarrcross2dfloz(:,bbintmera),intnum)
            call obj%clinterptypecross2dfhizomp%Init(ellarr,clarrcross2dfhiz(:,bbintmera),intnum)
            call obj%clinterptypecrosslowzgiomp%Init(ellarr,clarrcrosslowzgi(:,bbintmera),intnum)
            call obj%clinterptypecrosscmassgiomp%Init(ellarr,clarrcrosscmassgi(:,bbintmera),intnum)
            call obj%clinterptypecross2dflozgiomp%Init(ellarr,clarrcross2dflozgi(:,bbintmera),intnum)
            call obj%clinterptypecross2dfhizgiomp%Init(ellarr,clarrcross2dfhizgi(:,bbintmera),intnum)
        end if
        ellgentest = 1.0d0
        !$OMP PARALLEL DO FIRSTPRIVATE(obj,ellgen)
            do ellgen=1,58999
                clarrbig(ellgen) = obj%clinterptypeomp%Value(ellgen+1.0d0)
                clarrbigii(ellgen) = obj%clinterptypeiiomp%Value(ellgen+1.0d0)
                clarrbiggi(ellgen) = obj%clinterptypegiomp%Value(ellgen+1.0d0)
            end do
        !$OMP END PARALLEL DO
        if((m1arr(biin) == m2arr(biin)) .and. (obj%setscenarioomp > 1)) then
            ellgentest = 1.0d0
            !$OMP PARALLEL DO FIRSTPRIVATE(obj,ellgen)
                do ellgen=1,58999
                    clarrbigcrosslowz(ellgen) = obj%clinterptypecrosslowzomp%Value(ellgen+1.0d0)
                    clarrbigcrosscmass(ellgen) = obj%clinterptypecrosscmassomp%Value(ellgen+1.0d0)
                    clarrbigcross2dfloz(ellgen) = obj%clinterptypecross2dflozomp%Value(ellgen+1.0d0)
                    clarrbigcross2dfhiz(ellgen) = obj%clinterptypecross2dfhizomp%Value(ellgen+1.0d0)
                    clarrbigcrosslowzgi(ellgen) = obj%clinterptypecrosslowzgiomp%Value(ellgen+1.0d0)
                    clarrbigcrosscmassgi(ellgen) = obj%clinterptypecrosscmassgiomp%Value(ellgen+1.0d0)
                    clarrbigcross2dflozgi(ellgen) = obj%clinterptypecross2dflozgiomp%Value(ellgen+1.0d0)
                    clarrbigcross2dfhizgi(ellgen) = obj%clinterptypecross2dfhizgiomp%Value(ellgen+1.0d0)
                end do
            !$OMP END PARALLEL DO
        end if

        xiplus(:) = 0.0d0
        ximinus(:) = 0.0d0
        xiplusii(:) = 0.0d0
        ximinusii(:) = 0.0d0
        xiplusgi(:) = 0.0d0
        ximinusgi(:) = 0.0d0
        if((m1arr(biin) == m2arr(biin)) .and. (obj%setscenarioomp > 1)) then
            xicrosslowz(:) = 0.0d0
            xicrosscmass(:) = 0.0d0
            xicross2dfloz(:) = 0.0d0
            xicross2dfhiz(:) = 0.0d0
             xicrosslowzgi(:) = 0.0d0
             xicrosscmassgi(:) = 0.0d0
             xicross2dflozgi(:) = 0.0d0
             xicross2dfhizgi(:) = 0.0d0
        end if
        hoj=0
        do yeye=1,nangbins
            do loopee=1,nell
                xiplus(yeye) = xiplus(yeye) + 1.0d0/pisj/2.0d0*clarrbig(loopee)*ellgentestarr(loopee)*this%bes0arr(yeye,loopee)
                ximinus(yeye) = ximinus(yeye) + 1.0d0/pisj/2.0d0*clarrbig(loopee)*ellgentestarr(loopee)*this%bes4arr(yeye,loopee)
                xiplusii(yeye) = xiplusii(yeye) + 1.0d0/pisj/2.0d0*clarrbigii(loopee)*ellgentestarr(loopee)*this%bes0arr(yeye,loopee)
                ximinusii(yeye) = ximinusii(yeye) + 1.0d0/pisj/2.0d0*clarrbigii(loopee)*ellgentestarr(loopee)*this%bes4arr(yeye,loopee)
                xiplusgi(yeye) = xiplusgi(yeye) + 1.0d0/pisj/2.0d0*clarrbiggi(loopee)*ellgentestarr(loopee)*this%bes0arr(yeye,loopee)
                ximinusgi(yeye) = ximinusgi(yeye) + 1.0d0/pisj/2.0d0*clarrbiggi(loopee)*ellgentestarr(loopee)*this%bes4arr(yeye,loopee)
                if((m1arr(biin) == m2arr(biin)) .and. (obj%setscenarioomp > 1)) then
                    xicrosslowz(yeye) = xicrosslowz(yeye) + 1.0d0/pisj/2.0d0*clarrbigcrosslowz(loopee)*ellgentestarr(loopee)*this%bes2arr(yeye,loopee)
                    xicrosscmass(yeye) = xicrosscmass(yeye) + 1.0d0/pisj/2.0d0*clarrbigcrosscmass(loopee)*ellgentestarr(loopee)*this%bes2arr(yeye,loopee)
                    xicross2dfloz(yeye) = xicross2dfloz(yeye) + 1.0d0/pisj/2.0d0*clarrbigcross2dfloz(loopee)*ellgentestarr(loopee)*this%bes2arr(yeye,loopee)
                    xicross2dfhiz(yeye) = xicross2dfhiz(yeye) + 1.0d0/pisj/2.0d0*clarrbigcross2dfhiz(loopee)*ellgentestarr(loopee)*this%bes2arr(yeye,loopee)
                    xicrosslowzgi(yeye) = xicrosslowzgi(yeye) + 1.0d0/pisj/2.0d0*clarrbigcrosslowzgi(loopee)*ellgentestarr(loopee)*this%bes2arr(yeye,loopee)
                    xicrosscmassgi(yeye) = xicrosscmassgi(yeye) + 1.0d0/pisj/2.0d0*clarrbigcrosscmassgi(loopee)*ellgentestarr(loopee)*this%bes2arr(yeye,loopee)
                    xicross2dflozgi(yeye) = xicross2dflozgi(yeye) + 1.0d0/pisj/2.0d0*clarrbigcross2dflozgi(loopee)*ellgentestarr(loopee)*this%bes2arr(yeye,loopee)
                    xicross2dfhizgi(yeye) = xicross2dfhizgi(yeye) + 1.0d0/pisj/2.0d0*clarrbigcross2dfhizgi(loopee)*ellgentestarr(loopee)*this%bes2arr(yeye,loopee)
                end if
            end do
        end do

        !old Catherine's scheme
        xiplusminus(1+(biin-1)*nangbins*2:nangbins+(biin-1)*nangbins*2) = xiplus + xiplusii + xiplusgi !full xi+ accounting for intrinsic alignments
        xiplusminus((nangbins+1)+(biin-1)*nangbins*2:nangbins*2+(biin-1)*nangbins*2) = ximinus + ximinusii + ximinusgi !full xi- accounting for intrinsic alignments
        xiplusminusgg(1+(biin-1)*nangbins*2:nangbins+(biin-1)*nangbins*2) = xiplus
        xiplusminusgg((nangbins+1)+(biin-1)*nangbins*2:nangbins*2+(biin-1)*nangbins*2) = ximinus
        xiplusminusii(1+(biin-1)*nangbins*2:nangbins+(biin-1)*nangbins*2) = xiplusii
        xiplusminusii((nangbins+1)+(biin-1)*nangbins*2:nangbins*2+(biin-1)*nangbins*2) = ximinusii
        xiplusminusgi(1+(biin-1)*nangbins*2:nangbins+(biin-1)*nangbins*2) = xiplusgi
        xiplusminusgi((nangbins+1)+(biin-1)*nangbins*2:nangbins*2+(biin-1)*nangbins*2) = ximinusgi
        if((m1arr(biin) == m2arr(biin)) .and. (obj%setscenarioomp > 1)) then
            bbintmeraigen = bbintmeraigen+1
             xicrossarrlowz(1+(bbintmeraigen-1)*nangbins:nangbins+(bbintmeraigen-1)*nangbins) = xicrosslowz + xicrosslowzgi
             xicrossarrcmass(1+(bbintmeraigen-1)*nangbins:nangbins+(bbintmeraigen-1)*nangbins) = xicrosscmass + xicrosscmassgi
             xicrossarr2dfloz(1+(bbintmeraigen-1)*nangbins:nangbins+(bbintmeraigen-1)*nangbins) = xicross2dfloz + xicross2dflozgi
             xicrossarr2dfhiz(1+(bbintmeraigen-1)*nangbins:nangbins+(bbintmeraigen-1)*nangbins) = xicross2dfhiz + xicross2dfhizgi
        end if

    end do

    this%klinesum = this%klinesum + obj%kline


    finaltheory(1:size(xiplusminus)) = xiplusminus
    if(obj%setscenarioomp > 1) then !ordering is 2dfloz, 2dfhiz, cmass, lowz
        finaltheory(1+size(xiplusminus):size(xiplusminus)+size(xicrossarr2dfloz)) = xicrossarr2dfloz
        finaltheory(1+size(xiplusminus)+size(xicrossarr2dfloz):size(xiplusminus)+size(xicrossarr2dfloz)+size(xicrossarr2dfhiz)) = xicrossarr2dfhiz
        finaltheory(1+size(xiplusminus)+size(xicrossarr2dfloz)+size(xicrossarr2dfhiz):size(xiplusminus)+size(xicrossarr2dfloz)+size(xicrossarr2dfhiz)+size(xicrossarrcmass)) = xicrossarrcmass
        finaltheory(1+size(xiplusminus)+size(xicrossarr2dfloz)+size(xicrossarr2dfhiz)+size(xicrossarrcmass):size(xiplusminus)+size(xicrossarr2dfloz)+size(xicrossarr2dfhiz)+size(xicrossarrcmass)+size(xicrossarrlowz)) = xicrossarrlowz
    end if
    if(obj%setscenarioomp == 3) then
        finaltheory(1+size(xiplusminus)+size(xicrossarr2dfloz)+size(xicrossarr2dfhiz)+size(xicrossarrcmass)+size(xicrossarrlowz):size(xiplusminus)+size(xicrossarr2dfloz)+size(xicrossarr2dfhiz)+size(xicrossarrcmass)+size(xicrossarrlowz)+size(P0P2_final_2dfloz_overlap)) = P0P2_final_2dfloz_overlap
        finaltheory(1+size(xiplusminus)+size(xicrossarr2dfloz)+size(xicrossarr2dfhiz)+size(xicrossarrcmass)+size(xicrossarrlowz)+size(P0P2_final_2dfloz_overlap):size(xiplusminus)+size(xicrossarr2dfloz)+size(xicrossarr2dfhiz)+size(xicrossarrcmass)+size(xicrossarrlowz)+size(P0P2_final_2dfloz_overlap)+size(P0P2_final_2dfhiz_overlap)) = P0P2_final_2dfhiz_overlap
        finaltheory(1+size(xiplusminus)+size(xicrossarr2dfloz)+size(xicrossarr2dfhiz)+size(xicrossarrcmass)+size(xicrossarrlowz)+size(P0P2_final_2dfloz_overlap)+size(P0P2_final_2dfhiz_overlap):size(xiplusminus)+size(xicrossarr2dfloz)+size(xicrossarr2dfhiz)+size(xicrossarrcmass)+size(xicrossarrlowz)+size(P0P2_final_2dfloz_overlap)+size(P0P2_final_2dfhiz_overlap)+size(P0P2_final_cmass_overlap)) = P0P2_final_cmass_overlap
        finaltheory(1+size(xiplusminus)+size(xicrossarr2dfloz)+size(xicrossarr2dfhiz)+size(xicrossarrcmass)+size(xicrossarrlowz)+size(P0P2_final_2dfloz_overlap)+size(P0P2_final_2dfhiz_overlap)+size(P0P2_final_cmass_overlap):size(xiplusminus)+size(xicrossarr2dfloz)+size(xicrossarr2dfhiz)+size(xicrossarrcmass)+size(xicrossarrlowz)+size(P0P2_final_2dfloz_overlap)+size(P0P2_final_2dfhiz_overlap)+size(P0P2_final_cmass_overlap)+size(P0P2_final_lowz_overlap)) = P0P2_final_lowz_overlap
    end if


    howin = 1
    do how=1,this%sizcovpremask
        if(this%maskelements(how) == 1) then
            dataminustheory(howin) = this%xipm(howin) - finaltheory(how)
            dataminustheoryhoho(howin) = finaltheory(how)
            dataminustheoryhihi(howin) = this%xipm(howin)
            howin = howin + 1
        end if
    end do

    if(obj%setscenarioomp < 3) then
        likesj = (Matrix_QuadForm(this%invcovxipm,dataminustheory))/2.0d0 !compute likelihood
        chisqsj = likesj*2.0
        chisqsjsellentin = 930.0d0*dlog(1.0d0 + likesj*2.0d0/(930.0d0-1.0d0)) !KiDS
        CosmoJBD_LnLike = 930.0d0/2.0d0*dlog(1.0d0 + likesj*2.0d0/(930.0d0-1.0d0)) !!KiDS: Sellentin-Heavens
    end if

    if(obj%setscenarioomp == 3) then
        likesj = (Matrix_QuadForm(this%invcovxipm,dataminustheory))/2.0d0 !compute likelihood
        chisqsj = likesj*2.0
        chisqsjsellentin = 930.0d0*dlog(1.0d0 + likesj*2.0d0/(930.0d0-1.0d0)) !KiDS
        CosmoJBD_LnLike = 930.0d0/2.0d0*dlog(1.0d0 + likesj*2.0d0/(930.0d0-1.0d0)) !!KiDS: Sellentin-Heavens
    end if

    finish  = OMP_get_wtime()
    print *, 'finish-startrsd, total 2*CosmoJBD_LnLike,tol_error', finish-startrsd, CosmoJBD_LnLike*2.0d0,tol_error

    if(this%use_cmass_full) then
        DEALLOCATE(P0_theory_win_NGC)
        DEALLOCATE(P2_theory_win_NGC)
        DEALLOCATE(P0_theory_win_SGC)
        DEALLOCATE(P2_theory_win_SGC)
        DEALLOCATE(diff_NGC)
        DEALLOCATE(diff_SGC)
    end if
    if(this%use_wigglez_full) then
        DEALLOCATE(P0P2P4_lowz_final_wigz)
        DEALLOCATE(P0P2P4_highz_final_wigz)
        DEALLOCATE(P0P2P4_lowz_obs_wigz)
        DEALLOCATE(P0P2P4_highz_obs_wigz)
        DEALLOCATE(diff_lowz_wigz)
        DEALLOCATE(diff_highz_wigz)
        DEALLOCATE(like_lowZ)
        DEALLOCATE(like_highZ)
    end if
    if(this%use_cmass_overlap) then
        DEALLOCATE(P0P2_final_cmass_overlap)
        DEALLOCATE(diff_cmass_overlap)
    end if
    if(this%use_lowz_overlap) then
        DEALLOCATE(P0P2_final_lowz_overlap)
        DEALLOCATE(diff_lowz_overlap)
    end if
    if(this%use_2dfloz_overlap) then
        DEALLOCATE(P0P2_final_2dfloz_overlap)
        DEALLOCATE(diff_2dfloz_overlap)
    end if
    if(this%use_2dfhiz_overlap) then
        DEALLOCATE(P0P2_final_2dfhiz_overlap)
        DEALLOCATE(diff_2dfhiz_overlap)
    end if
!!!!!!!!!!!!!!!!!!!
    deallocate(obj%theoryomp)
    deallocate(obj%ellarromp)
    deallocate(obj%exact_z_index)
    deallocate(obj%psourcetypeomp)
    deallocate(obj%plenstypeomplowz)
    deallocate(obj%plenstypeompcmass)
    deallocate(obj%plenstypeomp2dfloz)
    deallocate(obj%plenstypeomp2dfhiz)
    deallocate(obj%momp)
    deallocate(obj%aphotarromp)
    deallocate(obj%weightarromp)
    deallocate(obj%weightarromplowz)
    deallocate(obj%weightarrompcmass)
    deallocate(obj%weightarromplenslowz)
    deallocate(obj%weightarromplenscmass)
    deallocate(obj%weightarromp2dfloz)
    deallocate(obj%weightarromp2dfhiz)
    deallocate(obj%weightarromplens2dfloz)
    deallocate(obj%weightarromplens2dfhiz)
    deallocate(obj%weightpsarromp)
    deallocate(obj%weightnoarromp)
    deallocate(obj%weightnoarromplenslowz)
    deallocate(obj%weightnoarromplenscmass)
    deallocate(obj%weightnoarromplens2dfloz)
    deallocate(obj%weightnoarromplens2dfhiz)
    deallocate(obj%m1arromp)
    deallocate(obj%m2arromp)
    deallocate(obj%clinterptypeomp)
    deallocate(obj%clinterptypeiiomp)
    deallocate(obj%clinterptypegiomp)
    deallocate(obj%clinterptypecrosslowzomp)
    deallocate(obj%clinterptypecrosscmassomp)
    deallocate(obj%clinterptypecross2dflozomp)
    deallocate(obj%clinterptypecross2dfhizomp)
    deallocate(obj%exact_zmodomp)
    deallocate(obj%exact_zmodomplenslowz)
    deallocate(obj%exact_zmodomplenscmass)
    deallocate(obj%gcubicomp)
    deallocate(obj%lumarromp)

    deallocate(clarrbig)
    deallocate(clarrbiggi)
    deallocate(clarrbigii)
    deallocate(clarrbigcrosslowz)
    deallocate(clarrbigcrosscmass)
    deallocate(clarrbigcross2dfloz)
    deallocate(clarrbigcross2dfhiz)
    deallocate(ellgentestarr)
    deallocate(xiplusminus)
    deallocate(xiplusminusgg)
    deallocate(xiplusminusii)
    deallocate(xiplusminusgi)
    deallocate(dataminustheory)
    deallocate(dataminustheoryhihi)
    deallocate(dataminustheoryhoho)
    deallocate(finaltheory)
    deallocate(xicrossarrlowz)
    deallocate(xicrossarrlowzgi)
    deallocate(xicrossarrcmass)
    deallocate(xicrossarrcmassgi)
    deallocate(xicrossarr2dfloz)
    deallocate(xicrossarr2dflozgi)
    deallocate(xicrossarr2dfhiz)
    deallocate(xicrossarr2dfhizgi)
    deallocate(clarr)
    deallocate(clarrii)
    deallocate(clarrgi)
    deallocate(clarrcrosslowz)
    deallocate(clarrcrosscmass)
    deallocate(clarrcross2dfloz)
    deallocate(clarrcross2dfhiz)
    deallocate(ellarr)
    deallocate(xiplus)
    deallocate(xiplusii)
    deallocate(xiplusgi)
    deallocate(xicrosslowz)
    deallocate(xicrosslowzgi)
    deallocate(xicrosscmass)
    deallocate(xicrosscmassgi)
    deallocate(xicross2dfloz)
    deallocate(xicross2dflozgi)
    deallocate(xicross2dfhiz)
    deallocate(xicross2dfhizgi)
    deallocate(ximinus)
    deallocate(ximinusii)
    deallocate(ximinusgi)
    deallocate(mumax1arr)
    deallocate(mumax2arr)
    deallocate(aphotarr)
    deallocate(m1arr)
    deallocate(m2arr)
    deallocate(mnomp)
    deallocate(rcslum)
    deallocate(cfhtlum)
    deallocate(lumarr)
    deallocate(weightarr)
    deallocate(weightarrlowz)
    deallocate(weightarrcmass)
    deallocate(weightarrlenslowz)
    deallocate(weightarrlenscmass)
    deallocate(weightarr2dfloz)
    deallocate(weightarr2dfhiz)
    deallocate(weightarrlens2dfloz)
    deallocate(weightarrlens2dfhiz)
    deallocate(weightpsarr)
    deallocate(weightnoarr)
    deallocate(weightnoarrlenslowz)
    deallocate(weightnoarrlenscmass)
    deallocate(weightnoarrlens2dfloz)
    deallocate(weightnoarrlens2dfhiz)
    deallocate(binwarr)
    deallocate(trapzarr)
    deallocate(trapzarrii)
    deallocate(trapzarrgi)
    deallocate(trapzarrcrosslowz)
    deallocate(trapzarrcrosscmass)
    deallocate(trapzarrcross2dfloz)
    deallocate(trapzarrcross2dfhiz)
    deallocate(garr)

    deallocate(obj%weightpsarromplowz)
    deallocate(obj%weightpsarrompcmass)
    deallocate(obj%weightpsarromp2dfloz)
    deallocate(obj%weightpsarromp2dfhiz)
    deallocate(obj%clinterptypecrosslowzgiomp)
    deallocate(obj%clinterptypecrosscmassgiomp)
    deallocate(obj%clinterptypecross2dflozgiomp)
    deallocate(obj%clinterptypecross2dfhizgiomp)
    deallocate(clarrbigcrosslowzgi)
    deallocate(clarrbigcrosscmassgi)
    deallocate(clarrbigcross2dflozgi)
    deallocate(clarrbigcross2dfhizgi)
    deallocate(clarrcrosslowzgi)
    deallocate(clarrcrosscmassgi)
    deallocate(clarrcross2dflozgi)
    deallocate(clarrcross2dfhizgi)
    deallocate(weightpsarrlowz)
    deallocate(weightpsarrcmass)
    deallocate(weightpsarr2dfloz)
    deallocate(weightpsarr2dfhiz)
    deallocate(trapzarrcrosslowzgi)
    deallocate(trapzarrcrosscmassgi)
    deallocate(trapzarrcross2dflozgi)
    deallocate(trapzarrcross2dfhizgi)

    contains 


!----------------------------------------------------------------------------------------------------------
! Theory calc CMASS
!----------------------------------------------------------------------------------------------------------

        subroutine Get_CMASS_TheoryMPK(P0_theory_win_NGC,P2_theory_win_NGC,P0_theory_win_SGC,P2_theory_win_SGC)
        REAL(mcp) rombint
    external rombint
        REAL(mcp), allocatable, dimension(:) :: P0_theory_win_NGC,P2_theory_win_NGC
    REAL(mcp), allocatable, dimension(:) :: P0_theory_win_SGC,P2_theory_win_SGC
        REAL(mcp), allocatable, dimension(:) :: P0_multipole_theory,P2_multipole_theory
        REAL(mcp), allocatable, dimension(:) :: k_grid_multipole,k_grid_array
        REAL(mcp) :: mu_max,mu_min,nu_total
    REAL(mcp) :: k_grid_int1, dk_grid1,k_grid_int2, dk_grid2,k_grid_max2
    REAL(mcp) :: k_prime,tol_error
    REAL(mcp) :: k_min_window_int,k_max_window_int,dk_window_int,k_temp1
        REAL(mcp) :: scale_multipole_L0,scale_multipole_L2,AP_scaling,scale_multipole_L4
    REAL(mcp), parameter :: H_fid_no_h=95.49, DA_fid_no_h=1345.87,h_fiduical=0.7    !! Fiducial parameters at z= 0.57 units [Mpc]
    REAL(mcp) :: H_fid, DA_fid
        INTEGER   :: k_index,k_num_window_int,k_num2

    H_fid =  H_fid_no_h/h_fiduical           !! Units [...h] Needed  for AP -- should be independent on H0
    DA_fid = DA_fid_no_h*h_fiduical          !! Units  [Mpc h^-1]

        !----------------------------------------------------------------------------------------------------------
    ! Setup Alcock-Pansiki parameters
    !----------------------------------------------------------------------------------------------------------

        a_perp_cmass    =  D_AhUnit(this%cmass_mp%z_eff,CMB)/DA_fid
        a_par_cmass     =  H_fid/HofzhUnit(this%cmass_mp%z_eff,CMB,c)

        !----------------------------------------------------------------------------------------------------------
    ! Setup Numerical Integration variables i.e., width,range ect..
    !----------------------------------------------------------------------------------------------------------

    ! k values to evaluate multipoles convolved with window function (plus integral constraint term)
    ! this needs to match the values of the observed multipoles
    k_grid_int1 = this%cmass_mp%k_min
    dk_grid1    = this%cmass_mp%k_spacing

    ! Integration over mu to find multipoles (over 2D power spectrum)
    mu_max = 1.0
    mu_min = -1.0

    ! k values to evaluate multipoles at (this is the theory prediction)
    k_grid_int2 = this%cmass_mp%kmin_theory
    dk_grid2 = this%cmass_mp%dk_theory
    k_grid_max2 = this%cmass_mp%kmax_theory
    k_num2 =  INT((k_grid_max2 - k_grid_int2)/(dk_grid2))


    ! Integration over k_prime, to convolve multipoles with window function (needs larger range
    ! than above given window functions are not delta functions)
    dk_window_int = 0.00010
    k_max_window_int = 0.25 !!0.399750  !! Fine as we only conisder k_max = 0.10
    k_min_window_int = 0.001 !! fine as we use k_min of 0.015

    k_num_window_int = INT((k_max_window_int - k_min_window_int)/(dk_window_int))

    nu_total = 1000

    allocate(k_grid_multipole(k_num2))
    allocate(P0_multipole_theory(k_num2))
    allocate(P2_multipole_theory(k_num2))
    allocate(k_grid_array(k_grid_num_cmass))

    !----------------------------------------------------------------------------------------------------------
    ! Setup wavenumber grids where multipole predictions will be calculated
    !----------------------------------------------------------------------------------------------------------

    do k_index = 1,k_grid_num_cmass
        k_val_cmass = k_grid_int1 + dk_grid1*(k_index - 1.0)
        k_grid_array(k_index) = k_val_cmass
    end do

    do k_index = 1,k_num2
        k_val_cmass = k_grid_int2 + dk_grid2*(k_index - 1.0)
        k_grid_multipole(k_index) = k_val_cmass
    end do

    !----------------------------------------------------------------------------------------------------------
    ! Calculate Theory predictions for Monopole and quadrupole
    !----------------------------------------------------------------------------------------------------------

        scale_multipole_L0 = 1.0/2.0
    scale_multipole_L2 = 5.0/2.0

    tol_error = 0.0001 ! Numerical integration relative error..

    AP_scaling = 1.0/(a_par_cmass*(a_perp_cmass**2))

    do k_index = 1,k_num2
        k_val_cmass = k_grid_multipole(k_index)
        P0_multipole_theory(k_index) = AP_scaling*scale_multipole_L0*rombint(P0_multipole_FN_BOSS,mu_min,mu_max,tol_error)
        P2_multipole_theory(k_index) = AP_scaling*scale_multipole_L2*rombint(P2_multipole_FN_BOSS,mu_min,mu_max,tol_error)
    end do

    !----------------------------------------------------------------------------------------------------------
    ! Setup interpolation for Monopole and quadrupole theory predictions
    !----------------------------------------------------------------------------------------------------------

    call P0_theory_spline_cmass%Init(k_grid_multipole,P0_multipole_theory,k_num2)
    call P2_theory_spline_cmass%Init(k_grid_multipole,P2_multipole_theory,k_num2)

        !----------------------------------------------------------------------------------------------------------
        ! Full convolution for NGC -- uses local functions
        !----------------------------------------------------------------------------------------------------------

        do k_index = 1,k_grid_num_cmass
                k_val_cmass = k_grid_array(k_index)

                P0_theory_win_NGC(k_index) = &
                2*piaj*rombint(P0_multipole_FN_NGC_WINDOW_2D,k_min_window_int,k_max_window_int,tol_error) - &
                2*piaj*(Win1D_NGC(1)%value(k_val_cmass)/W02_norm_IC_NGC)*rombint(P0_multipole_FN_NGC_WINDOW_1D_IC,k_min_window_int,k_max_window_int,tol_error)

                P2_theory_win_NGC(k_index) = &
                2*piaj*rombint(P2_multipole_FN_NGC_WINDOW_2D,k_min_window_int,k_max_window_int,tol_error) - &
                2*piaj*(Win1D_NGC(2)%value(k_val_cmass)/W02_norm_IC_NGC)*rombint(P2_multipole_FN_NGC_WINDOW_1D_IC,k_min_window_int,k_max_window_int,tol_error)

        end do

        !----------------------------------------------------------------------------------------------------------
        ! Full convolution for SGC -- uses local functions
        !----------------------------------------------------------------------------------------------------------

        do k_index = 1,k_grid_num_cmass
                k_val_cmass = k_grid_array(k_index)

                P0_theory_win_SGC(k_index) = &
                2*piaj*rombint(P0_multipole_FN_SGC_WINDOW_2D,k_min_window_int,k_max_window_int,tol_error) - &
                2*piaj*(Win1D_SGC(1)%value(k_val_cmass)/W02_norm_IC_SGC)*rombint(P0_multipole_FN_SGC_WINDOW_1D_IC,k_min_window_int,k_max_window_int,tol_error)

                P2_theory_win_SGC(k_index) = &
                2*piaj*rombint(P2_multipole_FN_SGC_WINDOW_2D,k_min_window_int,k_max_window_int,tol_error) - &
                2*piaj*(Win1D_SGC(2)%value(k_val_cmass)/W02_norm_IC_SGC)*rombint(P2_multipole_FN_SGC_WINDOW_1D_IC,k_min_window_int,k_max_window_int,tol_error)

        end do

    DEALLOCATE(k_grid_multipole)
    DEALLOCATE(P0_multipole_theory)
    DEALLOCATE(P2_multipole_theory)
    DEALLOCATE(k_grid_array)

        end subroutine

!----------------------------------------------------------------------------------------------------------
! Theory calculation WiggleZ
!----------------------------------------------------------------------------------------------------------

    subroutine Get_WiggleZ_TheoryMPK(P0P2P4_lowz_final_wigz,P0P2P4_highz_final_wigz)
    REAL(mcp), allocatable, dimension(:,:) :: P0P2P4_lowz_final_wigz,P0P2P4_highz_final_wigz
    REAL(mcp), allocatable, dimension(:)   :: P0_multipole_THEORY_lowZ,P2_multipole_THEORY_lowZ,P4_multipole_THEORY_lowZ
    REAL(mcp), allocatable, dimension(:)   :: P0_multipole_THEORY_highZ,P2_multipole_THEORY_highZ,P4_multipole_THEORY_highZ
    REAL(mcp), allocatable, dimension(:,:) :: P0P2_Conv_lowZ,P0P2_Conv_highZ
    REAL(mcp), allocatable, dimension(:)   :: Theory_vec_lowZ,Theory_vec_highZ
    REAL(mcp) :: mu_max,mu_min,mu_value,d_mu,nu_total,tol_error
    REAL(mcp) :: qnum,wnum
    REAL(mcp) :: scale_multipole_L0,scale_multipole_L2,scale_multipole_L4,AP_scaling_highZ,AP_scaling_lowZ
    REAL(mcp), parameter :: H_fid_lowZ_no_h = 87.5602, DA_fid_lowZ_no_h = 1174.2433, h_fiducial = 0.720     !! Fiducial parameters at z= 0.44
    REAL(mcp), parameter :: H_fid_highZ_no_h = 103.163, DA_fid_highZ_no_h = 1507.0329                       !! Fiducial parameters at z= 0.73
    REAL(mcp) :: H_fid_highZ,DA_fid_highZ,H_fid_lowZ,DA_fid_lowZ
    INTEGER   :: inum,k_index
    REAL(mcp) rombint
    external rombint

    allocate(P0_multipole_THEORY_lowZ(this%wigglez_mp%k_num_conv))
    allocate(P2_multipole_THEORY_lowZ(this%wigglez_mp%k_num_conv))
    allocate(P4_multipole_THEORY_lowZ(this%wigglez_mp%k_num_conv))
    allocate(P0_multipole_THEORY_highZ(this%wigglez_mp%k_num_conv))
    allocate(P2_multipole_THEORY_highZ(this%wigglez_mp%k_num_conv))
    allocate(P4_multipole_THEORY_highZ(this%wigglez_mp%k_num_conv))
    allocate(Theory_vec_lowZ(this%wigglez_mp%size_convolution))
    allocate(Theory_vec_highZ(this%wigglez_mp%size_convolution))
    allocate(P0P2_Conv_lowZ(this%wigglez_mp%size_convolution,this%wigglez_mp%num_regions))
    allocate(P0P2_Conv_highZ(this%wigglez_mp%size_convolution,this%wigglez_mp%num_regions))

    !----------------------------------------------------------------------------------------------------------
    ! Setup AP parameters -- can turn on or off
    !----------------------------------------------------------------------------------------------------------

    H_fid_highZ = H_fid_highZ_no_h/h_fiducial
    H_fid_lowZ = H_fid_lowZ_no_h/h_fiducial

    DA_fid_highZ = DA_fid_highZ_no_h*h_fiducial
    DA_fid_lowZ = DA_fid_lowZ_no_h*h_fiducial

        a_perp_lowz_wigz  = D_AhUnit(this%wigglez_mp%z_eff_low,CMB)/DA_fid_lowZ
        a_par_lowz_wigz   = H_fid_lowZ/HofzhUnit(this%wigglez_mp%z_eff_low,CMB,c)

    a_perp_highz_wigz  = D_AhUnit(this%wigglez_mp%z_eff_high,CMB)/DA_fid_highZ
    a_par_highz_wigz   = H_fid_highZ/HofzhUnit(this%wigglez_mp%z_eff_high,CMB,c)

    !----------------------------------------------------------------------------------------------------------
    ! start theory P0,P2,P4 calculation
    !----------------------------------------------------------------------------------------------------------

    tol_error = 0.000001
    mu_max = 1.0
    mu_min = -1.0

    scale_multipole_L0 = 1.0/2.0
    scale_multipole_L2 = 5.0/2.0
    scale_multipole_L4 = 9.0/2.0

    AP_scaling_highZ = 1.0/(a_par_highz_wigz*(a_perp_highz_wigz**2))
    AP_scaling_lowZ = 1.0/(a_par_lowz_wigz*(a_perp_lowz_wigz**2))

    do k_index = 1,this%wigglez_mp%k_num_conv

        k_val_wigz = this%wigglez_mp%k_min_theory + this%wigglez_mp%dk_theory*(k_index - 1.0)

        P0_multipole_THEORY_highZ(k_index) = AP_scaling_highZ*scale_multipole_L0*rombint(P0_multipole_FN_highZ,mu_min,mu_max,tol_error)
        P2_multipole_THEORY_highZ(k_index) = AP_scaling_highZ*scale_multipole_L2*rombint(P2_multipole_FN_highZ,mu_min,mu_max,tol_error)
        P4_multipole_THEORY_highZ(k_index) = AP_scaling_highZ*scale_multipole_L4*rombint(P4_multipole_FN_highZ,mu_min,mu_max,tol_error)

    end do

    do k_index = 1,this%wigglez_mp%k_num_conv

        k_val_wigz = this%wigglez_mp%k_min_theory + this%wigglez_mp%dk_theory*(k_index - 1.0)

        P0_multipole_THEORY_lowZ(k_index) = AP_scaling_lowZ*scale_multipole_L0*rombint(P0_multipole_FN_lowZ,mu_min,mu_max,tol_error)
        P2_multipole_THEORY_lowZ(k_index) = AP_scaling_lowZ*scale_multipole_L2*rombint(P2_multipole_FN_lowZ,mu_min,mu_max,tol_error)
        P4_multipole_THEORY_lowZ(k_index) = AP_scaling_lowZ*scale_multipole_L4*rombint(P4_multipole_FN_lowZ,mu_min,mu_max,tol_error)

    end do

    !----------------------------------------------------------------------------------------------------------
    ! Change data structure and Compute convolved multipoles -- done using convolution matrix
    !----------------------------------------------------------------------------------------------------------

    Theory_vec_lowZ(1:25)   = P0_multipole_THEORY_lowZ
    Theory_vec_lowZ(26:50)  = P2_multipole_THEORY_lowZ
    Theory_vec_lowZ(51:75)  = P4_multipole_THEORY_lowZ
    Theory_vec_highZ(1:25)  = P0_multipole_THEORY_highZ
    Theory_vec_highZ(26:50) = P2_multipole_THEORY_highZ
    Theory_vec_highZ(51:75) = P4_multipole_THEORY_highZ

    DO region_num = 1,this%wigglez_mp%num_regions
        CALL Matrix_MulVec(this%wigglez_mp%Covolution_matrix_lowZ(:,:,region_num),Theory_vec_lowZ(:),P0P2_Conv_lowZ(:,region_num))
        CALL Matrix_MulVec(this%wigglez_mp%Covolution_matrix_highZ(:,:,region_num),Theory_vec_highZ(:),P0P2_Conv_highZ(:,region_num))
    ENDDO

    !----------------------------------------------------------------------------------------------------------
    ! Again, changing format of theory and observations, stack into single vector P = [P0,P2,P4]
    !----------------------------------------------------------------------------------------------------------

    DO region_num = 1,this%wigglez_mp%num_regions
        DO inum=1,this%wigglez_mp%k_num_obs
            qnum = inum  + this%wigglez_mp%k_num_obs
            wnum = inum  + this%wigglez_mp%k_num_obs*2
            P0P2P4_lowz_obs_wigz(inum,region_num) =  this%wigglez_mp%P0_data_lowZ(inum,region_num)
            P0P2P4_highz_obs_wigz(inum,region_num) = this%wigglez_mp%P0_data_highZ(inum,region_num)
            P0P2P4_lowz_obs_wigz(qnum,region_num) = this%wigglez_mp%P2_data_lowZ(inum,region_num)
            P0P2P4_highz_obs_wigz(qnum,region_num) = this%wigglez_mp%P2_data_highZ(inum,region_num)
            P0P2P4_lowz_obs_wigz(wnum,region_num) = this%wigglez_mp%P4_data_lowZ(inum,region_num)
            P0P2P4_highz_obs_wigz(wnum,region_num) = this%wigglez_mp%P4_data_highZ(inum,region_num)
        ENDDO
    ENDDO

    DO region_num = 1,this%wigglez_mp%num_regions
        DO inum=1,this%wigglez_mp%k_num_obs

            P0P2P4_lowz_final_wigz(inum,region_num) = P0P2_Conv_lowZ(1 + inum,region_num)
            P0P2P4_lowz_final_wigz(inum + this%wigglez_mp%k_num_obs,region_num) = P0P2_Conv_lowZ(1 + inum + this%wigglez_mp%k_num_conv,region_num)
            P0P2P4_lowz_final_wigz(inum + this%wigglez_mp%k_num_obs*2,region_num) = P0P2_Conv_lowZ(1 + inum + 2*this%wigglez_mp%k_num_conv,region_num)

            P0P2P4_highz_final_wigz(inum,region_num) = P0P2_Conv_highZ(1 + inum,region_num)
            P0P2P4_highz_final_wigz(inum + this%wigglez_mp%k_num_obs,region_num) = P0P2_Conv_highZ(1 + inum + this%wigglez_mp%k_num_conv,region_num)
            P0P2P4_highz_final_wigz(inum + this%wigglez_mp%k_num_obs*2,region_num) = P0P2_Conv_highZ(1 + inum + 2*this%wigglez_mp%k_num_conv,region_num)
        ENDDO
    ENDDO

    DEALLOCATE(P0_multipole_THEORY_lowZ)
    DEALLOCATE(P2_multipole_THEORY_lowZ)
    DEALLOCATE(P4_multipole_THEORY_lowZ)
    DEALLOCATE(P0_multipole_THEORY_highZ)
    DEALLOCATE(P2_multipole_THEORY_highZ)
    DEALLOCATE(P4_multipole_THEORY_highZ)
    DEALLOCATE(Theory_vec_lowZ)
    DEALLOCATE(Theory_vec_highZ)
    DEALLOCATE(P0P2_Conv_lowZ)
    DEALLOCATE(P0P2_Conv_highZ)

    end subroutine

!----------------------------------------------------------------------------------------------------------
! Theory calculation CMASS - overlap region
!----------------------------------------------------------------------------------------------------------

    subroutine Get_CMASS_overlap_TheoryMPK(P0P2_final_cmass_overlap)
    REAL(mcp), allocatable, dimension(:) :: P0P2_final_cmass_overlap
    REAL(mcp), allocatable, dimension(:)   :: P0_multipole_theory,P2_multipole_theory,P4_multipole_theory
    REAL(mcp), allocatable, dimension(:)   :: P0P2P4_Conv
    REAL(mcp), allocatable, dimension(:)   :: theory_vec,k_values_obs,k_values_conv
    REAL(mcp) :: mu_max,mu_min,mu_value,d_mu,nu_total,tol_error
    REAL(mcp) :: qnum,wnum
    REAL(mcp) :: scale_multipole_L0,scale_multipole_L2,scale_multipole_L4,AP_scaling
    REAL(mcp), parameter :: H_fid_no_h=95.49, DA_fid_no_h=1345.9,h_fiduical=0.7    !! Fiducial parameters at z= 0.57 units [Mpc]
    REAL(mcp) :: H_fid,DA_fid,k_val,k_step
        INTEGER   :: inum,k_index,num_k
        REAL(mcp) rombint
    external rombint

        allocate(P0_multipole_theory(this%cmass_mp_overlap%k_num_conv))
    allocate(P2_multipole_theory(this%cmass_mp_overlap%k_num_conv))
    allocate(P4_multipole_theory(this%cmass_mp_overlap%k_num_conv))

    allocate(theory_vec(this%cmass_mp_overlap%size_convolution))
    allocate(P0P2P4_Conv(this%cmass_mp_overlap%size_convolution))
        allocate(k_values_obs(this%cmass_mp_overlap%k_num_obs))
        allocate(k_values_conv(this%cmass_mp_overlap%k_num_conv))

    !----------------------------------------------------------------------------------------------------------
    ! Compute k vectors (observed values and values for convolution matrix)
    !----------------------------------------------------------------------------------------------------------

        do k_index = 1,this%cmass_mp_overlap%k_num_obs
        k_val = this%cmass_mp_overlap%k_min_obs + this%cmass_mp_overlap%k_spacing_obs*(k_index - 1.0)
        k_values_obs(k_index) = k_val
    end do

	do k_index = 1,this%cmass_mp_overlap%k_num_conv
        k_val = this%cmass_mp_overlap%k_min_theory + this%cmass_mp_overlap%dk_theory*(k_index - 1.0)
        k_values_conv(k_index) = k_val
    end do

    k_values_conv(1) = exp(Theory%MPK%x(1))

    !----------------------------------------------------------------------------------------------------------
    ! Setup AP parameters -- can turn on or off
    !----------------------------------------------------------------------------------------------------------

        H_fid = H_fid_no_h/h_fiduical
        DA_fid = DA_fid_no_h*h_fiduical

        a_perp_cmass_overlap  = D_AhUnit(this%cmass_mp_overlap%z_eff,CMB)/DA_fid
        a_par_cmass_overlap   = H_fid/HofzhUnit(this%cmass_mp_overlap%z_eff,CMB,c)

    !----------------------------------------------------------------------------------------------------------
    ! start theory P0,P2,P4 calculation
    !----------------------------------------------------------------------------------------------------------

        tol_error = 0.000001
    mu_max = 1.0
    mu_min = -1.0

    scale_multipole_L0 = 1.0/2.0
    scale_multipole_L2 = 5.0/2.0
    scale_multipole_L4 = 9.0/2.0

    AP_scaling = 1.0/(a_par_cmass_overlap*(a_perp_cmass_overlap**2))

    do k_index = 1,this%cmass_mp_overlap%k_num_conv

!        k_val_overlap_cmass = this%cmass_mp_overlap%k_min_theory + this%cmass_mp_overlap%dk_theory*(k_index - 1.0)
        k_val_overlap_cmass = k_values_conv(k_index)

        P0_multipole_theory(k_index) = AP_scaling*scale_multipole_L0*rombint(P0_multipole_FN_cmass_overlap,mu_min,mu_max,tol_error)
        P2_multipole_theory(k_index) = AP_scaling*scale_multipole_L2*rombint(P2_multipole_FN_cmass_overlap,mu_min,mu_max,tol_error)
        P4_multipole_theory(k_index) = AP_scaling*scale_multipole_L4*rombint(P4_multipole_FN_cmass_overlap,mu_min,mu_max,tol_error)

    end do

    do k_index = 1,this%cmass_mp_overlap%k_num_conv
                k_val_overlap_cmass = this%cmass_mp_overlap%k_min_theory + this%cmass_mp_overlap%dk_theory*(k_index - 1.0)
    !            print*,'k_values, P0',k_val_overlap_cmass,P0_multipole_theory(k_index)
    end do

    do k_index = 1,this%cmass_mp_overlap%k_num_conv
                k_val_overlap_cmass = this%cmass_mp_overlap%k_min_theory + this%cmass_mp_overlap%dk_theory*(k_index - 1.0)
    !            print*,'k_values, P2',k_val_overlap_cmass,P2_multipole_theory(k_index)
    end do

    do k_index = 1,this%cmass_mp_overlap%k_num_conv
                k_val_overlap_cmass = this%cmass_mp_overlap%k_min_theory + this%cmass_mp_overlap%dk_theory*(k_index - 1.0)
    !            print*,'k_values, P4',k_val_overlap_cmass,P4_multipole_theory(k_index)
    end do

    !------------------------------------------------------------------------------------------------------------------------
    ! Stack into single vector and Compute convolved multipoles -- done using convolution matrix,
    ! note,     P0P2P4_Conv(:) is the convolved vector!
    !-------------------------------------------------------------------------------------------------------------------------

    theory_vec(1:10)   = P0_multipole_theory
    theory_vec(11:20)  = P2_multipole_theory
    theory_vec(21:30)  = P4_multipole_theory

    CALL Matrix_MulVec(this%cmass_mp_overlap%Covolution_matrix(:,:),theory_vec(:),P0P2P4_Conv(:))

        !print*,'convolved theory vector is',P0P2P4_Conv

    !----------------------------------------------------------------------------------------------------------
    ! Stack observations into single vector P = [P0,P2,P4]
    !----------------------------------------------------------------------------------------------------------

    num_k = this%cmass_mp_overlap%k_num_obs

!    P0P2_obs_cmass_overlap(1:num_k) = this%cmass_mp_overlap%P0_data(:)
!    P0P2_obs_cmass_overlap(num_k+1:2*num_k) = this%cmass_mp_overlap%P2_data(:)

        !----------------------------------------------------------------------------------------------------------
    ! Setup interpolation for (convolved) Monopole and quadrupole theory predictions and evaluation at require values.
    !----------------------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------------

    call P0_theory_spline_cmass_overlap%Init(k_values_conv,P0P2P4_Conv(1:10),this%cmass_mp_overlap%k_num_conv)
    call P2_theory_spline_cmass_overlap%Init(k_values_conv,P0P2P4_Conv(11:20),this%cmass_mp_overlap%k_num_conv)

        !print*,'check input range',P0P2P4_Conv(11:20)
        !print*,'check input k range',k_values_conv(:)
        !print*,'check num k values input',this%cmass_mp_overlap%k_num_conv

    DO inum=1,this%cmass_mp_overlap%k_num_obs
                k_val = k_values_obs(inum)
                k_step = this%cmass_mp_overlap%k_num_obs
            P0P2_final_cmass_overlap(inum) = P0_theory_spline_cmass_overlap%value(k_val)
            P0P2_final_cmass_overlap(inum + k_step) = P2_theory_spline_cmass_overlap%value(k_val)
       !     print*,'k value is',k_val
       !     print*,'computed P0 is',P0_theory_spline_cmass_overlap%value(k_val)
       !     print*,'computed P2 is',P2_theory_spline_cmass_overlap%value(k_val)
    ENDDO

	!print*,'Computed data vector is',P0P2_final_cmass_overlap(:)

    DEALLOCATE(P0_multipole_theory)
    DEALLOCATE(P2_multipole_theory)
    DEALLOCATE(P4_multipole_theory)
    DEALLOCATE(theory_vec)
    DEALLOCATE(P0P2P4_Conv)
    DEALLOCATE(k_values_obs)
    DEALLOCATE(k_values_conv)

    end subroutine

!----------------------------------------------------------------------------------------------------------
! Theory calculation LOWZ - overlap region
!----------------------------------------------------------------------------------------------------------

    subroutine Get_LOWZ_overlap_TheoryMPK(P0P2_final_lowz_overlap)
    REAL(mcp), allocatable, dimension(:) :: P0P2_final_lowz_overlap
    REAL(mcp), allocatable, dimension(:)   :: P0_multipole_theory,P2_multipole_theory,P4_multipole_theory
    REAL(mcp), allocatable, dimension(:)   :: P0P2P4_Conv
    REAL(mcp), allocatable, dimension(:)   :: theory_vec,k_values_obs,k_values_conv
    REAL(mcp) :: mu_max,mu_min,mu_value,d_mu,nu_total,tol_error
    REAL(mcp) :: qnum,wnum
    REAL(mcp) :: scale_multipole_L0,scale_multipole_L2,scale_multipole_L4,AP_scaling
    REAL(mcp), parameter :: H_fid_no_h=82.53, DA_fid_no_h=960.2, h_fiduical=0.7  !! Fiducial parameters at z= 0.32 units [Mpc]
    REAL(mcp) :: H_fid,DA_fid,k_val,k_step
    INTEGER   :: inum,k_index,num_k
    REAL(mcp) rombint
    external rombint

    allocate(P0_multipole_theory(this%lowz_mp_overlap%k_num_conv))
    allocate(P2_multipole_theory(this%lowz_mp_overlap%k_num_conv))
    allocate(P4_multipole_theory(this%lowz_mp_overlap%k_num_conv))

    allocate(theory_vec(this%lowz_mp_overlap%size_convolution))
    allocate(P0P2P4_Conv(this%lowz_mp_overlap%size_convolution))

        allocate(k_values_obs(this%lowz_mp_overlap%k_num_obs))
        allocate(k_values_conv(this%lowz_mp_overlap%k_num_conv))

        !----------------------------------------------------------------------------------------------------------
    ! Compute k vectors (observed values and values for convolution matrix)
    !----------------------------------------------------------------------------------------------------------

        do k_index = 1,this%lowz_mp_overlap%k_num_obs
        k_val = this%lowz_mp_overlap%k_min_obs + this%lowz_mp_overlap%k_spacing_obs*(k_index - 1.0)
        k_values_obs(k_index) = k_val
    end do

	do k_index = 1,this%lowz_mp_overlap%k_num_conv
        k_val = this%lowz_mp_overlap%k_min_theory + this%lowz_mp_overlap%dk_theory*(k_index - 1.0)
        k_values_conv(k_index) = k_val
    end do

    k_values_conv(1) = exp(Theory%MPK%x(1))

    !----------------------------------------------------------------------------------------------------------
    ! Setup AP parameters -- can turn on or off
    !----------------------------------------------------------------------------------------------------------

    H_fid = H_fid_no_h/h_fiduical
    DA_fid = DA_fid_no_h*h_fiduical

    a_perp_lowz_overlap  = D_AhUnit(this%lowz_mp_overlap%z_eff,CMB)/DA_fid
    a_par_lowz_overlap   = H_fid/HofzhUnit(this%lowz_mp_overlap%z_eff,CMB,c)

    !----------------------------------------------------------------------------------------------------------
    ! start theory P0,P2,P4 calculation
    !----------------------------------------------------------------------------------------------------------

    tol_error = 0.000001
    mu_max = 1.0
    mu_min = -1.0

    scale_multipole_L0 = 1.0/2.0
    scale_multipole_L2 = 5.0/2.0
    scale_multipole_L4 = 9.0/2.0
    scale_multipole_L4 = 9.0/2.0

    AP_scaling = 1.0/(a_par_lowz_overlap*(a_perp_lowz_overlap**2))

    do k_index = 1,this%lowz_mp_overlap%k_num_conv

        k_val_overlap_lowz = k_values_conv(k_index)

        P0_multipole_theory(k_index) = AP_scaling*scale_multipole_L0*rombint(P0_multipole_FN_lowz_overlap,mu_min,mu_max,tol_error)
        P2_multipole_theory(k_index) = AP_scaling*scale_multipole_L2*rombint(P2_multipole_FN_lowz_overlap,mu_min,mu_max,tol_error)
        P4_multipole_theory(k_index) = AP_scaling*scale_multipole_L4*rombint(P4_multipole_FN_lowz_overlap,mu_min,mu_max,tol_error)

    end do

    do k_index = 1,this%lowz_mp_overlap%k_num_conv
                k_val_overlap_lowz = this%lowz_mp_overlap%k_min_theory + this%lowz_mp_overlap%dk_theory*(k_index - 1.0)
    end do

    do k_index = 1,this%lowz_mp_overlap%k_num_conv
                k_val_overlap_lowz = this%lowz_mp_overlap%k_min_theory + this%lowz_mp_overlap%dk_theory*(k_index - 1.0)
    end do

    do k_index = 1,this%lowz_mp_overlap%k_num_conv
                k_val_overlap_lowz = this%lowz_mp_overlap%k_min_theory + this%lowz_mp_overlap%dk_theory*(k_index - 1.0)
    end do

    !------------------------------------------------------------------------------------------------------------------------
    ! Stack into single vector and Compute convolved multipoles -- done using convolution matrix,
    ! note,     P0P2P4_Conv(:) is the convolved vector!
    !-------------------------------------------------------------------------------------------------------------------------

    theory_vec(1:10)   = P0_multipole_theory
    theory_vec(11:20)  = P2_multipole_theory
    theory_vec(21:30)  = P4_multipole_theory

    CALL Matrix_MulVec(this%lowz_mp_overlap%Covolution_matrix(:,:),theory_vec(:),P0P2P4_Conv(:))

    !----------------------------------------------------------------------------------------------------------
    ! Stack observations into single vector P = [P0,P2,P4]
    !----------------------------------------------------------------------------------------------------------

    num_k = this%lowz_mp_overlap%k_num_obs

    !----------------------------------------------------------------------------------------------------------
    ! Setup interpolation for (convolved) Monopole and quadrupole theory predictions and evaluation at require values.
    !----------------------------------------------------------------------------------------------------------

    call P0_theory_spline_lowz_overlap%Init(k_values_conv,P0P2P4_Conv(1:10),this%lowz_mp_overlap%k_num_conv)
    call P2_theory_spline_lowz_overlap%Init(k_values_conv,P0P2P4_Conv(11:20),this%lowz_mp_overlap%k_num_conv)

    DO inum=1,this%lowz_mp_overlap%k_num_obs
                k_val = k_values_obs(inum)
                k_step = this%lowz_mp_overlap%k_num_obs
            P0P2_final_lowz_overlap(inum) = P0_theory_spline_lowz_overlap%value(k_val)
            P0P2_final_lowz_overlap(inum + k_step) = P2_theory_spline_lowz_overlap%value(k_val)
    ENDDO

    DEALLOCATE(P0_multipole_theory)
    DEALLOCATE(P2_multipole_theory)
    DEALLOCATE(P4_multipole_theory)
    DEALLOCATE(theory_vec)
    DEALLOCATE(P0P2P4_Conv)
    DEALLOCATE(k_values_obs)
    DEALLOCATE(k_values_conv)

    end subroutine



!----------------------------------------------------------------------------------------------------------
! Theory calculation 2dfloz - overlap region
!----------------------------------------------------------------------------------------------------------

    subroutine Get_2dfloz_overlap_TheoryMPK(P0P2_final_2dfloz_overlap)
    REAL(mcp), allocatable, dimension(:) :: P0P2_final_2dfloz_overlap
    REAL(mcp), allocatable, dimension(:)   :: P0_multipole_theory,P2_multipole_theory,P4_multipole_theory
    REAL(mcp), allocatable, dimension(:)   :: P0P2P4_Conv
    REAL(mcp), allocatable, dimension(:)   :: theory_vec,k_values_obs,k_values_conv
    REAL(mcp) :: mu_max,mu_min,mu_value,d_mu,nu_total,tol_error
    REAL(mcp) :: qnum,wnum
    REAL(mcp) :: scale_multipole_L0,scale_multipole_L2,scale_multipole_L4,AP_scaling
    REAL(mcp), parameter :: H_fid_no_h=82.07, DA_fid_no_h=939.7, h_fiduical=0.7  !! Fiducial parameters at z= 0.31 units [Mpc]
    REAL(mcp) :: H_fid,DA_fid,k_val,k_step
    INTEGER   :: inum,k_index,num_k
    REAL(mcp) rombint
    external rombint

    allocate(P0_multipole_theory(this%twodfloz_mp_overlap%k_num_conv))
    allocate(P2_multipole_theory(this%twodfloz_mp_overlap%k_num_conv))
    allocate(P4_multipole_theory(this%twodfloz_mp_overlap%k_num_conv))

    allocate(theory_vec(this%twodfloz_mp_overlap%size_convolution))
    allocate(P0P2P4_Conv(this%twodfloz_mp_overlap%size_convolution))

    allocate(k_values_obs(this%twodfloz_mp_overlap%k_num_obs))
    allocate(k_values_conv(this%twodfloz_mp_overlap%k_num_conv))

    !----------------------------------------------------------------------------------------------------------
    ! Compute k vectors (observed values and values for convolution matrix)
    !----------------------------------------------------------------------------------------------------------

    do k_index = 1,this%twodfloz_mp_overlap%k_num_obs
        k_val = this%twodfloz_mp_overlap%k_min_obs + this%twodfloz_mp_overlap%k_spacing_obs*(k_index - 1.0)
        k_values_obs(k_index) = k_val
    end do

    do k_index = 1,this%twodfloz_mp_overlap%k_num_conv
        k_val = this%twodfloz_mp_overlap%k_min_theory + this%twodfloz_mp_overlap%dk_theory*(k_index - 1.0)
        k_values_conv(k_index) = k_val
    end do

    k_values_conv(1) = exp(Theory%MPK%x(1))

    !----------------------------------------------------------------------------------------------------------
    ! Setup AP parameters -- can turn on or off
    !----------------------------------------------------------------------------------------------------------

    H_fid = H_fid_no_h/h_fiduical
    DA_fid = DA_fid_no_h*h_fiduical

    a_perp_2dfloz_overlap  = D_AhUnit(this%twodfloz_mp_overlap%z_eff,CMB)/DA_fid
    a_par_2dfloz_overlap   = H_fid/HofzhUnit(this%twodfloz_mp_overlap%z_eff,CMB,c)

    !----------------------------------------------------------------------------------------------------------
    ! start theory P0,P2,P4 calculation
    !----------------------------------------------------------------------------------------------------------

    tol_error = 0.000001
    mu_max = 1.0
    mu_min = -1.0

    scale_multipole_L0 = 1.0/2.0
    scale_multipole_L2 = 5.0/2.0
    scale_multipole_L4 = 9.0/2.0
    scale_multipole_L4 = 9.0/2.0

    AP_scaling = 1.0/(a_par_2dfloz_overlap*(a_perp_2dfloz_overlap**2))

    do k_index = 1,this%twodfloz_mp_overlap%k_num_conv

        k_val_overlap_2dfloz = k_values_conv(k_index)

        P0_multipole_theory(k_index) = AP_scaling*scale_multipole_L0*rombint(P0_multipole_FN_2dfloz_overlap,mu_min,mu_max,tol_error)
        P2_multipole_theory(k_index) = AP_scaling*scale_multipole_L2*rombint(P2_multipole_FN_2dfloz_overlap,mu_min,mu_max,tol_error)
        P4_multipole_theory(k_index) = AP_scaling*scale_multipole_L4*rombint(P4_multipole_FN_2dfloz_overlap,mu_min,mu_max,tol_error)

    end do

    do k_index = 1,this%twodfloz_mp_overlap%k_num_conv
                k_val_overlap_2dfloz = this%twodfloz_mp_overlap%k_min_theory + this%twodfloz_mp_overlap%dk_theory*(k_index - 1.0)
    end do

    do k_index = 1,this%twodfloz_mp_overlap%k_num_conv
                k_val_overlap_2dfloz = this%twodfloz_mp_overlap%k_min_theory + this%twodfloz_mp_overlap%dk_theory*(k_index - 1.0)
    end do

    do k_index = 1,this%twodfloz_mp_overlap%k_num_conv
                k_val_overlap_2dfloz = this%twodfloz_mp_overlap%k_min_theory + this%twodfloz_mp_overlap%dk_theory*(k_index - 1.0)
    end do

    !------------------------------------------------------------------------------------------------------------------------
    ! Stack into single vector and Compute convolved multipoles -- done using convolution matrix,
    ! note,     P0P2P4_Conv(:) is the convolved vector!
    !-------------------------------------------------------------------------------------------------------------------------

    theory_vec(1:10)   = P0_multipole_theory
    theory_vec(11:20)  = P2_multipole_theory
    theory_vec(21:30)  = P4_multipole_theory

    CALL Matrix_MulVec(this%twodfloz_mp_overlap%Covolution_matrix(:,:),theory_vec(:),P0P2P4_Conv(:))

    !----------------------------------------------------------------------------------------------------------
    ! Stack observations into single vector P = [P0,P2,P4]
    !----------------------------------------------------------------------------------------------------------

    num_k = this%twodfloz_mp_overlap%k_num_obs

    !P0P2_obs_2dfloz_overlap(1:num_k) = this%twodfloz_mp_overlap%P0_data(:)
    !P0P2_obs_2dfloz_overlap(num_k+1:2*num_k) = this%twodfloz_mp_overlap%P2_data(:)

    !----------------------------------------------------------------------------------------------------------
    ! Setup interpolation for (convolved) Monopole and quadrupole theory predictions and evaluation at require values.
    !----------------------------------------------------------------------------------------------------------

    call P0_theory_spline_2dfloz_overlap%Init(k_values_conv,P0P2P4_Conv(1:10),this%twodfloz_mp_overlap%k_num_conv)
    call P2_theory_spline_2dfloz_overlap%Init(k_values_conv,P0P2P4_Conv(11:20),this%twodfloz_mp_overlap%k_num_conv)

    DO inum=1,this%twodfloz_mp_overlap%k_num_obs
                k_val = k_values_obs(inum)
                k_step = this%twodfloz_mp_overlap%k_num_obs
            P0P2_final_2dfloz_overlap(inum) = P0_theory_spline_2dfloz_overlap%value(k_val)
            P0P2_final_2dfloz_overlap(inum + k_step) = P2_theory_spline_2dfloz_overlap%value(k_val)
    ENDDO

    DEALLOCATE(P0_multipole_theory)
    DEALLOCATE(P2_multipole_theory)
    DEALLOCATE(P4_multipole_theory)
    DEALLOCATE(theory_vec)
    DEALLOCATE(P0P2P4_Conv)
    DEALLOCATE(k_values_obs)
    DEALLOCATE(k_values_conv)

    end subroutine



!----------------------------------------------------------------------------------------------------------
! Theory calculation 2dfhiz - overlap region
!----------------------------------------------------------------------------------------------------------

    subroutine Get_2dfhiz_overlap_TheoryMPK(P0P2_final_2dfhiz_overlap)
    REAL(mcp), allocatable, dimension(:) :: P0P2_final_2dfhiz_overlap
    REAL(mcp), allocatable, dimension(:)   :: P0_multipole_theory,P2_multipole_theory,P4_multipole_theory
    REAL(mcp), allocatable, dimension(:)   :: P0P2P4_Conv
    REAL(mcp), allocatable, dimension(:)   :: theory_vec,k_values_obs,k_values_conv
    REAL(mcp) :: mu_max,mu_min,mu_value,d_mu,nu_total,tol_error
    REAL(mcp) :: qnum,wnum
    REAL(mcp) :: scale_multipole_L0,scale_multipole_L2,scale_multipole_L4,AP_scaling
    REAL(mcp), parameter :: H_fid_no_h=94.92, DA_fid_no_h=1334.3,h_fiduical=0.7    !! Fiducial parameters at z= 0.56 units [Mpc]
    REAL(mcp) :: H_fid,DA_fid,k_val,k_step
    INTEGER   :: inum,k_index,num_k
    REAL(mcp) rombint
    external rombint

    allocate(P0_multipole_theory(this%twodfhiz_mp_overlap%k_num_conv))
    allocate(P2_multipole_theory(this%twodfhiz_mp_overlap%k_num_conv))
    allocate(P4_multipole_theory(this%twodfhiz_mp_overlap%k_num_conv))

    allocate(theory_vec(this%twodfhiz_mp_overlap%size_convolution))
    allocate(P0P2P4_Conv(this%twodfhiz_mp_overlap%size_convolution))

    allocate(k_values_obs(this%twodfhiz_mp_overlap%k_num_obs))
    allocate(k_values_conv(this%twodfhiz_mp_overlap%k_num_conv))

    !----------------------------------------------------------------------------------------------------------
    ! Compute k vectors (observed values and values for convolution matrix)
    !----------------------------------------------------------------------------------------------------------

    do k_index = 1,this%twodfhiz_mp_overlap%k_num_obs
        k_val = this%twodfhiz_mp_overlap%k_min_obs + this%twodfhiz_mp_overlap%k_spacing_obs*(k_index - 1.0)
        k_values_obs(k_index) = k_val
    end do

    do k_index = 1,this%twodfhiz_mp_overlap%k_num_conv
        k_val = this%twodfhiz_mp_overlap%k_min_theory + this%twodfhiz_mp_overlap%dk_theory*(k_index - 1.0)
        k_values_conv(k_index) = k_val
    end do

    k_values_conv(1) = exp(Theory%MPK%x(1))

    !----------------------------------------------------------------------------------------------------------
    ! Setup AP parameters -- can turn on or off
    !----------------------------------------------------------------------------------------------------------

    H_fid = H_fid_no_h/h_fiduical
    DA_fid = DA_fid_no_h*h_fiduical

    a_perp_2dfhiz_overlap  = D_AhUnit(this%twodfhiz_mp_overlap%z_eff,CMB)/DA_fid
    a_par_2dfhiz_overlap   = H_fid/HofzhUnit(this%twodfhiz_mp_overlap%z_eff,CMB,c)

    !----------------------------------------------------------------------------------------------------------
    ! start theory P0,P2,P4 calculation
    !----------------------------------------------------------------------------------------------------------

    tol_error = 0.000001
    mu_max = 1.0
    mu_min = -1.0

    scale_multipole_L0 = 1.0/2.0
    scale_multipole_L2 = 5.0/2.0
    scale_multipole_L4 = 9.0/2.0
    scale_multipole_L4 = 9.0/2.0

    AP_scaling = 1.0/(a_par_2dfhiz_overlap*(a_perp_2dfhiz_overlap**2))

    do k_index = 1,this%twodfhiz_mp_overlap%k_num_conv

        k_val_overlap_2dfhiz = k_values_conv(k_index)

        P0_multipole_theory(k_index) = AP_scaling*scale_multipole_L0*rombint(P0_multipole_FN_2dfhiz_overlap,mu_min,mu_max,tol_error)
        P2_multipole_theory(k_index) = AP_scaling*scale_multipole_L2*rombint(P2_multipole_FN_2dfhiz_overlap,mu_min,mu_max,tol_error)
        P4_multipole_theory(k_index) = AP_scaling*scale_multipole_L4*rombint(P4_multipole_FN_2dfhiz_overlap,mu_min,mu_max,tol_error)

    end do

    do k_index = 1,this%twodfhiz_mp_overlap%k_num_conv
                k_val_overlap_2dfhiz = this%twodfhiz_mp_overlap%k_min_theory + this%twodfhiz_mp_overlap%dk_theory*(k_index - 1.0)
    end do

    do k_index = 1,this%twodfhiz_mp_overlap%k_num_conv
                k_val_overlap_2dfhiz = this%twodfhiz_mp_overlap%k_min_theory + this%twodfhiz_mp_overlap%dk_theory*(k_index - 1.0)
    end do

    do k_index = 1,this%twodfhiz_mp_overlap%k_num_conv
                k_val_overlap_2dfhiz = this%twodfhiz_mp_overlap%k_min_theory + this%twodfhiz_mp_overlap%dk_theory*(k_index - 1.0)
    end do

    !------------------------------------------------------------------------------------------------------------------------
    ! Stack into single vector and Compute convolved multipoles -- done using convolution matrix,
    ! note,     P0P2P4_Conv(:) is the convolved vector!
    !-------------------------------------------------------------------------------------------------------------------------

    theory_vec(1:10)   = P0_multipole_theory
    theory_vec(11:20)  = P2_multipole_theory
    theory_vec(21:30)  = P4_multipole_theory

    CALL Matrix_MulVec(this%twodfhiz_mp_overlap%Covolution_matrix(:,:),theory_vec(:),P0P2P4_Conv(:))

    !----------------------------------------------------------------------------------------------------------
    ! Stack observations into single vector P = [P0,P2,P4]
    !----------------------------------------------------------------------------------------------------------

    num_k = this%twodfhiz_mp_overlap%k_num_obs

    !P0P2_obs_2dfhiz_overlap(1:num_k) = this%twodfhiz_mp_overlap%P0_data(:)
    !P0P2_obs_2dfhiz_overlap(num_k+1:2*num_k) = this%twodfhiz_mp_overlap%P2_data(:)

    !----------------------------------------------------------------------------------------------------------
    ! Setup interpolation for (convolved) Monopole and quadrupole theory predictions and evaluation at require values.
    !----------------------------------------------------------------------------------------------------------

    call P0_theory_spline_2dfhiz_overlap%Init(k_values_conv,P0P2P4_Conv(1:10),this%twodfhiz_mp_overlap%k_num_conv)
    call P2_theory_spline_2dfhiz_overlap%Init(k_values_conv,P0P2P4_Conv(11:20),this%twodfhiz_mp_overlap%k_num_conv)

    DO inum=1,this%twodfhiz_mp_overlap%k_num_obs
            k_val = k_values_obs(inum)
            k_step = this%twodfhiz_mp_overlap%k_num_obs
            P0P2_final_2dfhiz_overlap(inum) = P0_theory_spline_2dfhiz_overlap%value(k_val)
            P0P2_final_2dfhiz_overlap(inum + k_step) = P2_theory_spline_2dfhiz_overlap%value(k_val)
    ENDDO

    DEALLOCATE(P0_multipole_theory)
    DEALLOCATE(P2_multipole_theory)
    DEALLOCATE(P4_multipole_theory)
    DEALLOCATE(theory_vec)
    DEALLOCATE(P0P2P4_Conv)
    DEALLOCATE(k_values_obs)
    DEALLOCATE(k_values_conv)

    end subroutine


    !----------------------------------------------------------------------------------------------------------
    ! functions for calculating theory multipoles for WiggleZ -- using 2D power spectrum (local functions)
    !----------------------------------------------------------------------------------------------------------

    function P0_multipole_FN_highZ(mu)
        REAL(mcp) :: mu,P0_multipole_FN_highZ
        P0_multipole_FN_highZ = &
        LinearPower_2D_s_wigz(Theory,k_val_wigz,mu,z_highz_wigz,sig_v2_wigz,a_perp_highz_wigz,a_par_highz_wigz,b_highz_wigz)*Legendre_Pn(mu,0)
    END function P0_multipole_FN_highZ

    function P2_multipole_FN_highZ(mu)
        REAL(mcp) :: mu,P2_multipole_FN_highZ

        P2_multipole_FN_highZ = &
        LinearPower_2D_s_wigz(Theory,k_val_wigz,mu,z_highz_wigz,sig_v2_wigz,a_perp_highz_wigz,a_par_highz_wigz,b_highz_wigz)*Legendre_Pn(mu,2)
    end function P2_multipole_FN_highZ

    function P4_multipole_FN_highZ(mu)
        REAL(mcp) :: mu,P4_multipole_FN_highZ

        P4_multipole_FN_highZ = &
        LinearPower_2D_s_wigz(Theory,k_val_wigz,mu,z_highz_wigz,sig_v2_wigz,a_perp_highz_wigz,a_par_highz_wigz,b_highz_wigz)*Legendre_Pn(mu,4)
    end function P4_multipole_FN_highZ

    function P0_multipole_FN_lowZ(mu)
        REAL(mcp) :: mu,P0_multipole_FN_lowZ

        P0_multipole_FN_lowZ = &
        LinearPower_2D_s_wigz(Theory,k_val_wigz,mu,z_lowz_wigz,sig_v1_wigz,a_perp_lowz_wigz,a_par_lowz_wigz,b_lowz_wigz)*Legendre_Pn(mu,0)
    END function P0_multipole_FN_lowZ

    function P2_multipole_FN_lowZ(mu)
        REAL(mcp) :: mu,P2_multipole_FN_lowZ

        P2_multipole_FN_lowZ = &
        LinearPower_2D_s_wigz(Theory,k_val_wigz,mu,z_lowz_wigz,sig_v1_wigz,a_perp_lowz_wigz,a_par_lowz_wigz,b_lowz_wigz)*Legendre_Pn(mu,2)
    end function P2_multipole_FN_lowZ

    function P4_multipole_FN_lowZ(mu)
        REAL(mcp) :: mu,P4_multipole_FN_lowZ

        P4_multipole_FN_lowZ = &
        LinearPower_2D_s_wigz(Theory,k_val_wigz,mu,z_lowz_wigz,sig_v1_wigz,a_perp_lowz_wigz,a_par_lowz_wigz,b_lowz_wigz)*Legendre_Pn(mu,4)
    end function P4_multipole_FN_lowZ

        !----------------------------------------------------------------------------------------------------------
        ! functions for calculating theory multipoles for CMASS -- using 2D power spectrum (local functions)
        !----------------------------------------------------------------------------------------------------------

        function P0_multipole_FN_BOSS(mu)
        REAL(mcp) :: mu,P0_multipole_FN_BOSS

        P0_multipole_FN_BOSS = &
        LinearPower_2D_s_cmass(Theory,k_val_cmass,mu,z_cmass,sigma_v_cmass,a_perp_cmass,a_par_cmass,b1_cmass,N_shot_cmass)*Legendre_Pn(mu,0)
        END function P0_multipole_FN_BOSS

        function P2_multipole_FN_BOSS(mu)
        REAL(mcp) :: mu,P2_multipole_FN_BOSS

        P2_multipole_FN_BOSS = &
        LinearPower_2D_s_cmass(Theory,k_val_cmass,mu,z_cmass,sigma_v_cmass,a_perp_cmass,a_par_cmass,b1_cmass,N_shot_cmass)*Legendre_Pn(mu,2)
        end function P2_multipole_FN_BOSS

        function P4_multipole_FN_BOSS(mu)
        REAL(mcp) :: mu,P4_multipole_FN_BOSS

        P4_multipole_FN_BOSS = &
        LinearPower_2D_s_cmass(Theory,k_val_cmass,mu,z_cmass,sigma_v_cmass,a_perp_cmass,a_par_cmass,b1_cmass,N_shot_cmass)*Legendre_Pn(mu,4)
        end function P4_multipole_FN_BOSS

        !----------------------------------------------------------------------------------------------------------
        ! functions for calculating theory multipoles for CMASS overlap region-- using 2D power spectrum (local functions)
        !----------------------------------------------------------------------------------------------------------

        function P0_multipole_FN_cmass_overlap(mu)
        REAL(mcp) :: mu,P0_multipole_FN_cmass_overlap

        P0_multipole_FN_cmass_overlap = &
        LinearPower_2D_s_cmass(Theory,k_val_overlap_cmass,mu,z_cmass,sigma_v_cmass,a_perp_cmass_overlap,a_par_cmass_overlap,b1_cmass,N_shot_cmass)*Legendre_Pn(mu,0)
        END function P0_multipole_FN_cmass_overlap

        function P2_multipole_FN_cmass_overlap(mu)
        REAL(mcp) :: mu,P2_multipole_FN_cmass_overlap

        P2_multipole_FN_cmass_overlap = &
        LinearPower_2D_s_cmass(Theory,k_val_overlap_cmass,mu,z_cmass,sigma_v_cmass,a_perp_cmass_overlap,a_par_cmass_overlap,b1_cmass,N_shot_cmass)*Legendre_Pn(mu,2)
        end function P2_multipole_FN_cmass_overlap

        function P4_multipole_FN_cmass_overlap(mu)
        REAL(mcp) :: mu,P4_multipole_FN_cmass_overlap

        P4_multipole_FN_cmass_overlap = &
        LinearPower_2D_s_cmass(Theory,k_val_overlap_cmass,mu,z_cmass,sigma_v_cmass,a_perp_cmass_overlap,a_par_cmass_overlap,b1_cmass,N_shot_cmass)*Legendre_Pn(mu,4)
        end function P4_multipole_FN_cmass_overlap

        !----------------------------------------------------------------------------------------------------------
        ! functions for calculating theory multipoles for LOWZ overlap region-- using 2D power spectrum (local functions)
        !----------------------------------------------------------------------------------------------------------

        function P0_multipole_FN_lowz_overlap(mu)
        REAL(mcp) :: mu,P0_multipole_FN_lowz_overlap

        P0_multipole_FN_lowz_overlap = &
        LinearPower_2D_s_cmass(Theory,k_val_overlap_lowz,mu,z_lowz,sigv_lowz,a_perp_lowz_overlap,a_par_lowz_overlap,b1_lowz,N_shot_lowz)*Legendre_Pn(mu,0)
        END function P0_multipole_FN_lowz_overlap

        function P2_multipole_FN_lowz_overlap(mu)
        REAL(mcp) :: mu,P2_multipole_FN_lowz_overlap

        P2_multipole_FN_lowz_overlap = &
        LinearPower_2D_s_cmass(Theory,k_val_overlap_lowz,mu,z_lowz,sigv_lowz,a_perp_lowz_overlap,a_par_lowz_overlap,b1_lowz,N_shot_lowz)*Legendre_Pn(mu,2)
        end function P2_multipole_FN_lowz_overlap

        function P4_multipole_FN_lowz_overlap(mu)
        REAL(mcp) :: mu,P4_multipole_FN_lowz_overlap

        P4_multipole_FN_lowz_overlap = &
        LinearPower_2D_s_cmass(Theory,k_val_overlap_lowz,mu,z_lowz,sigv_lowz,a_perp_lowz_overlap,a_par_lowz_overlap,b1_lowz,N_shot_lowz)*Legendre_Pn(mu,4)
        end function P4_multipole_FN_lowz_overlap


        !----------------------------------------------------------------------------------------------------------
        ! functions for calculating theory multipoles for 2dfloz overlap region-- using 2D power spectrum (local functions)
        !----------------------------------------------------------------------------------------------------------

        function P0_multipole_FN_2dfloz_overlap(mu)
        REAL(mcp) :: mu,P0_multipole_FN_2dfloz_overlap

        P0_multipole_FN_2dfloz_overlap = &
        LinearPower_2D_s_cmass(Theory,k_val_overlap_2dfloz,mu,z_2dfloz,sigv_2dfloz,a_perp_2dfloz_overlap,a_par_2dfloz_overlap,b1_2dfloz,N_shot_2dfloz)*Legendre_Pn(mu,0)
        END function P0_multipole_FN_2dfloz_overlap

        function P2_multipole_FN_2dfloz_overlap(mu)
        REAL(mcp) :: mu,P2_multipole_FN_2dfloz_overlap

        P2_multipole_FN_2dfloz_overlap = &
        LinearPower_2D_s_cmass(Theory,k_val_overlap_2dfloz,mu,z_2dfloz,sigv_2dfloz,a_perp_2dfloz_overlap,a_par_2dfloz_overlap,b1_2dfloz,N_shot_2dfloz)*Legendre_Pn(mu,2)
        end function P2_multipole_FN_2dfloz_overlap

        function P4_multipole_FN_2dfloz_overlap(mu)
        REAL(mcp) :: mu,P4_multipole_FN_2dfloz_overlap

        P4_multipole_FN_2dfloz_overlap = &
        LinearPower_2D_s_cmass(Theory,k_val_overlap_2dfloz,mu,z_2dfloz,sigv_2dfloz,a_perp_2dfloz_overlap,a_par_2dfloz_overlap,b1_2dfloz,N_shot_2dfloz)*Legendre_Pn(mu,4)
        end function P4_multipole_FN_2dfloz_overlap


        !----------------------------------------------------------------------------------------------------------
        ! functions for calculating theory multipoles for 2dfhiz overlap region-- using 2D power spectrum (local functions)
        !----------------------------------------------------------------------------------------------------------

        function P0_multipole_FN_2dfhiz_overlap(mu)
        REAL(mcp) :: mu,P0_multipole_FN_2dfhiz_overlap

        P0_multipole_FN_2dfhiz_overlap = &
        LinearPower_2D_s_cmass(Theory,k_val_overlap_2dfhiz,mu,z_2dfhiz,sigv_2dfhiz,a_perp_2dfhiz_overlap,a_par_2dfhiz_overlap,b1_2dfhiz,N_shot_2dfhiz)*Legendre_Pn(mu,0)
        END function P0_multipole_FN_2dfhiz_overlap

        function P2_multipole_FN_2dfhiz_overlap(mu)
        REAL(mcp) :: mu,P2_multipole_FN_2dfhiz_overlap

        P2_multipole_FN_2dfhiz_overlap = &
        LinearPower_2D_s_cmass(Theory,k_val_overlap_2dfhiz,mu,z_2dfhiz,sigv_2dfhiz,a_perp_2dfhiz_overlap,a_par_2dfhiz_overlap,b1_2dfhiz,N_shot_2dfhiz)*Legendre_Pn(mu,2)
        end function P2_multipole_FN_2dfhiz_overlap

        function P4_multipole_FN_2dfhiz_overlap(mu)
        REAL(mcp) :: mu,P4_multipole_FN_2dfhiz_overlap

        P4_multipole_FN_2dfhiz_overlap = &
        LinearPower_2D_s_cmass(Theory,k_val_overlap_2dfhiz,mu,z_2dfhiz,sigv_2dfhiz,a_perp_2dfhiz_overlap,a_par_2dfhiz_overlap,b1_2dfhiz,N_shot_2dfhiz)*Legendre_Pn(mu,4)
        end function P4_multipole_FN_2dfhiz_overlap

        !----------------------------------------------------------------------------------------------------------
        ! functions for calculating convolved multipoles NGC
        !----------------------------------------------------------------------------------------------------------

        function P0_multipole_FN_NGC_WINDOW_2D(k_prime)
        REAL(mcp) :: P0_multipole_FN_NGC_WINDOW_2D,k_prime

        P0_multipole_FN_NGC_WINDOW_2D = &
        k_prime**2*(P0_theory_spline_cmass%value(k_prime)*Win2D_NGC(1,1)%value(k_val_cmass,k_prime) + &
        P2_theory_spline_cmass%value(k_prime)*Win2D_NGC(1,2)%value(k_val_cmass,k_prime))

        end function P0_multipole_FN_NGC_WINDOW_2D

        function P0_multipole_FN_NGC_WINDOW_1D_IC(k_prime)
        REAL(mcp) :: P0_multipole_FN_NGC_WINDOW_1D_IC,k_prime

        P0_multipole_FN_NGC_WINDOW_1D_IC = (k_prime**2)*( &
        scale_factor_L0*P0_theory_spline_cmass%value(k_prime)*Win1D_NGC(1)%value(k_prime) + &
        scale_factor_L2*P2_theory_spline_cmass%value(k_prime)*Win1D_NGC(2)%value(k_prime))

        end function P0_multipole_FN_NGC_WINDOW_1D_IC

        function P2_multipole_FN_NGC_WINDOW_2D(k_prime)
        REAL(mcp) :: P2_multipole_FN_NGC_WINDOW_2D,k_prime

        P2_multipole_FN_NGC_WINDOW_2D = &
        k_prime**2*(P0_theory_spline_cmass%value(k_prime)*Win2D_NGC(2,1)%value(k_val_cmass,k_prime) + &
        P2_theory_spline_cmass%value(k_prime)*Win2D_NGC(2,2)%value(k_val_cmass,k_prime))

        end function P2_multipole_FN_NGC_WINDOW_2D

        function P2_multipole_FN_NGC_WINDOW_1D_IC(k_prime)
        REAL(mcp) :: k_prime,P2_multipole_FN_NGC_WINDOW_1D_IC

        P2_multipole_FN_NGC_WINDOW_1D_IC =(k_prime**2)*( &
        scale_factor_L0*P0_theory_spline_cmass%value(k_prime)*Win1D_NGC(1)%value(k_prime) + &
        scale_factor_L2*P2_theory_spline_cmass%value(k_prime)*Win1D_NGC(2)%value(k_prime))

        end function P2_multipole_FN_NGC_WINDOW_1D_IC

        !----------------------------------------------------------------------------------------------------------
        ! functions for calculating convolved multipoles SGC
        !----------------------------------------------------------------------------------------------------------

        function P0_multipole_FN_SGC_WINDOW_2D(k_prime)
        REAL(mcp) :: P0_multipole_FN_SGC_WINDOW_2D,k_prime
        P0_multipole_FN_SGC_WINDOW_2D = &
        k_prime**2*(P0_theory_spline_cmass%value(k_prime)*Win2D_SGC(1,1)%value(k_val_cmass,k_prime) + &
        P2_theory_spline_cmass%value(k_prime)*Win2D_SGC(1,2)%value(k_val_cmass,k_prime))
        end function P0_multipole_FN_SGC_WINDOW_2D

        function P0_multipole_FN_SGC_WINDOW_1D_IC(k_prime)
        REAL(mcp) :: P0_multipole_FN_SGC_WINDOW_1D_IC,k_prime

        P0_multipole_FN_SGC_WINDOW_1D_IC = (k_prime**2)*( &
        scale_factor_L0*P0_theory_spline_cmass%value(k_prime)*Win1D_SGC(1)%value(k_prime) + &
        scale_factor_L2*P2_theory_spline_cmass%value(k_prime)*Win1D_SGC(2)%value(k_prime))

        end function P0_multipole_FN_SGC_WINDOW_1D_IC

        function P2_multipole_FN_SGC_WINDOW_2D(k_prime)
        REAL(mcp) :: P2_multipole_FN_SGC_WINDOW_2D,k_prime

        P2_multipole_FN_SGC_WINDOW_2D = &
        k_prime**2*(P0_theory_spline_cmass%value(k_prime)*Win2D_SGC(2,1)%value(k_val_cmass,k_prime) + &
        P2_theory_spline_cmass%value(k_prime)*Win2D_SGC(2,2)%value(k_val_cmass,k_prime))

        end function P2_multipole_FN_SGC_WINDOW_2D

        function P2_multipole_FN_SGC_WINDOW_1D_IC(k_prime)
        REAL(mcp) :: k_prime,P2_multipole_FN_SGC_WINDOW_1D_IC

        P2_multipole_FN_SGC_WINDOW_1D_IC =(k_prime**2)*( &
        scale_factor_L0*P0_theory_spline_cmass%value(k_prime)*Win1D_SGC(1)%value(k_prime) + &
        scale_factor_L2*P2_theory_spline_cmass%value(k_prime)*Win1D_SGC(2)%value(k_prime))

        end function P2_multipole_FN_SGC_WINDOW_1D_IC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function sjgrowtha(reda) !integrate to obtain growth function D(0)/D(z)
        REAL(mcp) :: reda,sjgrowtha

        sjgrowtha = (1.0d0/reda)*Theory%growth_z%Value(1.0d0/reda-1.0d0)/Theory%sigma8_z%Value(1.0d0/reda-1.0d0)
    END function sjgrowtha

    !non-weight component for trapezoid integration
    function sjclsobjnoweight(obj,zlens)
        type(my_type) :: obj
        REAL(mcp) :: zlens,sjclsobjnoweight,weightpart,weight1,weight2,distz,kval,deltafinal,hubblez,ckms
        real(mcp) :: kminsj,kmaxsj

        ckms = 299792.458d0
        distz = ComovingRadialDistance(zlens)
        obj%distlensflatomp = distz
        distz = f_K(distz) !with curvature        
        kval = (obj%ellarromp(obj%ellgenomp)+1.0d0/2.0d0)/((obj%h0omp/100.0d0)*distz)

        kminsj = 1.0d-5
        kmaxsj = 100.0d0
        deltafinal = 0.0d0
        if(kval < kminsj) obj%kline = obj%kline + 1
        if (obj%eftflag == 2 .and. (obj%designereftmodel == 4 .or. obj%designereftmodel == 5)) then
           if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = obj%theoryomp%NL_MPK_WEYL%PowerAt(kval,zlens)
        else
           if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = obj%theoryomp%NL_MPK%PowerAt(kval,zlens)
        end if 
        
        hubblez = Hofz(zlens)
        weightpart = 3.0d0/2.0d0*(obj%omdmomp+obj%ombomp)*(obj%h0omp/ckms)**2.0d0*distz*(1.0d0+zlens)
        if (obj%eftflag == 2 .and. (obj%designereftmodel == 4 .or. obj%designereftmodel == 5)) then
           sjclsobjnoweight = 1.0d0/(distz**2.0d0)/hubblez*deltafinal/(obj%homp)**3.0d0/(9.0d0/4.0d0*(obj%h0omp/ckms)**4.0d0*(obj%omdmomp+obj%ombomp)**2.0d0/(obj%homp)**3.0d0*(1.0d0+zlens)**2.0d0)
        else 
           sjclsobjnoweight = 1.0d0/(distz**2.0d0)/hubblez*deltafinal/(obj%homp)**3.0d0
        end if
    END function sjclsobjnoweight


    !weight for trapezoid integration
    function sjclsobjonlyweight(obj,zlens)
        type(my_type) :: obj
        REAL(mcp) :: zlens,sjclsobjonlyweight,weightpart,weight1,weight2,distz,kval,deltafinal,hubblez,ckms
        real(mcp) :: kminsj,kmaxsj

        ckms = 299792.458d0
        distz = ComovingRadialDistance(zlens)
        obj%distlensflatomp = distz !assuming flatness
        distz = f_K(distz) !with curvature        
        weightpart = 3.0d0/2.0d0*(obj%omdmomp+obj%ombomp)*(obj%h0omp/ckms)**2.0d0*distz*(1.0d0+zlens)
        obj%wchooseomp = 1 !doesn't matter here because imposed momp(1) = momp(2)
        sjclsobjonlyweight = weightpart*rombint_obj(obj,weightobjcubic,zlens,obj%mumaxwomp,obj%tol_erroromp) !note changed mumax to mumax1

    END function sjclsobjonlyweight

    !IA-weight for trapezoid integration
    function sjclsiiobjonlyweight(obj,zlens)
        type(my_type) :: obj
        REAL(mcp) :: zlens,sjclsiiobjonlyweight,kval,deltafinal,hubblez,growthsf,growthnorm
        real(mcp) :: kminsj,kmaxsj

        growthnorm = obj%gcubicomp%Value(zlens)	!this is D(0)/D(z)
	hubblez = Hofz(zlens)
        obj%wchooseomp = 1 !doesn't matter here
        sjclsiiobjonlyweight = psourceobjcubic(obj,zlens)*hubblez*(-obj%ampiayesomp*5.0d-14*2.77536627d11*(obj%omdmomp+obj%ombomp)*growthnorm*((1.0d0+zlens)/(1.0d0+0.3d0))**obj%redziayesomp)*(obj%lumarromp(obj%momp(1)))**obj%lumiayesomp !enforced momp(1) = momp(2), so either is fine

    END function sjclsiiobjonlyweight


    !WEAK LENSING INTEGRAND (GG) ---- wrt scale factor instead
    function sjclsobjsf(obj,alens)
        type(my_type) :: obj
        REAL(mcp) :: alens,zlens,sjclsobjsf,weightpart,weight1,weight2,distz,kval,deltafinal,hubblez,ckms
        real(mcp) :: kminsj,kmaxsj

        zlens = 1.0d0/alens - 1.0d0
        ckms = 299792.458d0
        distz = ComovingRadialDistance(zlens)
        obj%distlensflatomp = distz !assuming flatness
        distz = f_K(distz) !with curvature        
        kval = (obj%ellarromp(obj%ellgenomp)+1.0d0/2.0d0)/((obj%h0omp/100.0d0)*distz)
        kminsj = 1.0d-5
        kmaxsj = 100.0d0
        deltafinal = 0.0d0
        if(kval < kminsj) obj%kline = obj%kline + 1

        if (obj%eftflag == 2 .and. (obj%designereftmodel == 4 .or. obj%designereftmodel == 5)) then
           if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = obj%theoryomp%NL_MPK_WEYL%PowerAt(kval,zlens)
        else
           if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = obj%theoryomp%NL_MPK%PowerAt(kval,zlens)
	end if

        hubblez = Hofz(zlens)
        weightpart = 3.0d0/2.0d0*(obj%omdmomp+obj%ombomp)*(obj%h0omp/ckms)**2.0d0*distz*(1.0d0+zlens)
        obj%wchooseomp = 1
        weight1 = weightpart*rombint_obj(obj,weightobjcubic,zlens,obj%mumax1omp,obj%tol_erroromp) !note changed mumax to mumax1
        if(obj%bnumomp == 1 .or. obj%bnumomp == 5 .or. obj%bnumomp == 8 .or. obj%bnumomp == 10) then !generalize later
            weight2 = weight1
        else
            obj%wchooseomp = 2
            weight2 = weightpart*rombint_obj(obj,weightobjcubic,zlens,obj%mumax2omp,obj%tol_erroromp)
        end if

        if (obj%eftflag == 2 .and. (obj%designereftmodel == 4 .or. obj%designereftmodel == 5)) then
           sjclsobjsf = ((1.0d0+zlens)**2.0d0)*weight1*weight2/(distz**2.0d0)/hubblez*deltafinal/(obj%homp)**3.0d0/(9.0d0/4.0d0*(obj%h0omp/ckms )**4.0d0*(obj%omdmomp+obj%ombomp)**2.0d0/(obj%homp)**3.0d0*( 1.0d0+zlens)**2.0d0)
        else
           sjclsobjsf = ((1.0d0+zlens)**2.0d0)*weight1*weight2/(distz**2.0d0)/hubblez*deltafinal/(obj%homp)**3.0d0
        end if

    END function sjclsobjsf


    !INTRINSIC ALIGNMENT INTEGRAND (II) ---- wrt scale factor instead
    function sjclsiiobjsf(obj,alens)
        type(my_type) :: obj
        REAL(mcp) :: alens,zlens,sjclsiiobjsf,growthnorm,weightpart,weight1,weight2,distz,kval,deltafinal,hubblez,ckms,growthsf
        real(mcp) :: kminsj,kmaxsj

        zlens = 1.0d0/alens - 1.0d0
        ckms = 299792.458d0
        distz = ComovingRadialDistance(zlens)
        obj%distlensflatomp = distz !assuming flatness ---- need to fix this
        distz = f_K(distz) !with curvature        
        kval = (obj%ellarromp(obj%ellgenomp)+1.0d0/2.0d0)/((obj%h0omp/100.0d0)*distz)
        kminsj = 1.0d-5
        kmaxsj = 100.0d0
        deltafinal = 0.0d0

        if (obj%eftflag == 2 .and. (obj%designereftmodel == 4 .or. obj%designereftmodel == 5)) then
        if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = obj%theoryomp%NL_MPK_WEYL%PowerAt(kval,zlens)
        else
        if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = obj%theoryomp%NL_MPK%PowerAt(kval,zlens)
        end if 

        growthnorm = obj%gcubicomp%Value(zlens)
        deltafinal = deltafinal*((obj%lumarromp(obj%momp(1))*obj%lumarromp(obj%momp(2)))**obj%lumiayesomp)*(-obj%ampiayesomp*5.0d-14*2.77536627d11*(obj%omdmomp+obj%ombomp)*growthnorm*((1.0d0+zlens)/(1.0d0+0.3d0))**obj%redziayesomp)**2.0d0
	hubblez = Hofz(zlens)
        obj%wchooseomp = 1
        weight1 = psourceobjcubic(obj,zlens)*hubblez
        if(obj%bnumomp == 1 .or. obj%bnumomp == 5 .or. obj%bnumomp == 8 .or. obj%bnumomp == 10) then !generalize later
            weight2 = weight1
        else
            obj%wchooseomp = 2
            weight2 = psourceobjcubic(obj,zlens)*hubblez
        end if

        if (obj%eftflag == 2 .and. (obj%designereftmodel == 4 .or. obj%designereftmodel == 5)) then
        sjclsiiobjsf = ((1.0d0+zlens)**2.0d0)*weight1*weight2/(distz**2.0d0)/hubblez*deltafinal/(obj%homp)**3.0d0/(9.0d0/4.0d0*(obj%h0omp/ckms )**4.0d0*(obj%omdmomp+obj%ombomp)**2.0d0/(obj%homp)**3.0d0*( 1.0d0+zlens)**2.0d0)
        else
        sjclsiiobjsf = ((1.0d0+zlens)**2.0d0)*weight1*weight2/(distz**2.0d0)/hubblez*deltafinal/(obj%homp)**3.0d0
        end if

    END function sjclsiiobjsf


    !LENSING - INTRINSIC ALIGNMENT INTEGRAND (GI) --- wrt scale factor instead
    function sjclsgiobjsf(obj,alens)
        type(my_type) :: obj
        REAL(mcp) :: alens,zlens,sjclsgiobjsf,weightpart,weight1,weight2,distz,kval,deltafinal,hubblez,ckms,weight3,weight4,growthsf,growthnorm
        real(mcp) :: kminsj,kmaxsj

        zlens = 1.0d0/alens - 1.0d0
        ckms = 299792.458d0
        distz = ComovingRadialDistance(zlens)
        obj%distlensflatomp = distz !assuming flatness
        distz = f_K(distz) !with curvature        
        kval = (obj%ellarromp(obj%ellgenomp)+1.0d0/2.0d0)/((obj%h0omp/100.0d0)*distz)
        kminsj = 1.0d-5
        kmaxsj = 100.0d0
        deltafinal = 0.0d0
        if (obj%eftflag == 2 .and. (obj%designereftmodel == 4 .or. obj%designereftmodel == 5)) then
        if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = obj%theoryomp%NL_MPK_WEYL%PowerAt(kval,zlens)
        else
        if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = obj%theoryomp%NL_MPK%PowerAt(kval,zlens)
        end if
        growthnorm = obj%gcubicomp%Value(zlens)	!this is D(0)/D(z)
        deltafinal = deltafinal*(-obj%ampiayesomp*5.0d-14*2.77536627d11*(obj%omdmomp+obj%ombomp)*growthnorm*((1.0d0+zlens)/(1.0d0+0.3d0))**obj%redziayesomp)
	hubblez = Hofz(zlens)
        weightpart = 3.0d0/2.0d0*(obj%omdmomp+obj%ombomp)*(obj%h0omp/ckms)**2.0d0*distz*(1.0d0+zlens)
        obj%wchooseomp = 1
        weight1 = weightpart*rombint_obj(obj,weightobjcubic,zlens,obj%mumax1omp,obj%tol_erroromp) !note changed mumax to mumax1
        weight3 = psourceobjcubic(obj,zlens)*hubblez
        if(obj%bnumomp == 1 .or. obj%bnumomp == 5 .or. obj%bnumomp == 8 .or. obj%bnumomp == 10) then !generalize later
            weight2 = weight1
            weight4 = weight3
        else
            obj%wchooseomp = 2
            weight2 = weightpart*rombint_obj(obj,weightobjcubic,zlens,obj%mumax2omp,obj%tol_erroromp) !note changed mumax to mumax2
            weight4 = psourceobjcubic(obj,zlens)*hubblez
        end if
        if (obj%eftflag == 2 .and. (obj%designereftmodel == 4 .or. obj%designereftmodel == 5)) then
        sjclsgiobjsf = ((1.0d0+zlens)**2.0d0)*(weight1*weight4*(obj%lumarromp(obj%momp(2)))**obj%lumiayesomp + weight2*weight3*(obj%lumarromp(obj%momp(1)))**obj%lumiayesomp)/(distz**2.0d0)/hubblez*deltafinal/(obj%homp)**3.0d0/(9.0d0/4.0d0*(obj%h0omp/ckms )**4.0d0*(obj%omdmomp+obj%ombomp)**2.0d0/(obj%homp)**3.0d0*( 1.0d0+zlens)**2.0d0)
        else
        sjclsgiobjsf = ((1.0d0+zlens)**2.0d0)*(weight1*weight4*(obj%lumarromp(obj%momp(2)))**obj%lumiayesomp + weight2*weight3*(obj%lumarromp(obj%momp(1)))**obj%lumiayesomp)/(distz**2.0d0)/hubblez*deltafinal/(obj%homp)**3.0d0
        end if
    END function sjclsgiobjsf


    !LENSING-GALAXY INTEGRAND (Gg)
    function sjclscrossobjnoweight(obj,zlens)
        type(my_type) :: obj
        REAL(mcp) :: zlens,sjclscrossobjnoweight,weightpart,weight1,weight2,distz,kval,deltafinal,hubblez
        real(mcp) :: kminsj,kmaxsj,ckms
        
        ckms = 299792.458d0
        hubblez = Hofz(zlens)
        distz = ComovingRadialDistance(zlens)
        obj%distlensflatomp = distz
        distz = f_K(distz) !with curvature
        kval = (obj%ellarromp(obj%ellgenomp)+1.0d0/2.0d0)/((obj%h0omp/100.0d0)*distz)

        kminsj = 1.0d-5
        kmaxsj = 100.0d0
        deltafinal = 0.0d0

        if (obj%eftflag == 2 .and. (obj%designereftmodel == 4 .or. obj%designereftmodel == 5)) then
        if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = obj%theoryomp%NL_MPK_WEYL%PowerAt(kval,zlens)
        else
        if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = obj%theoryomp%NL_MPK%PowerAt(kval,zlens)
        end if

        if (obj%eftflag == 2 .and. (obj%designereftmodel == 4 .or. obj%designereftmodel == 5)) then
        sjclscrossobjnoweight = 1.0d0/(distz**2.0d0)/hubblez*deltafinal/(obj%homp**3.0d0)/(9.0d0/4.0d0*(obj%h0omp/ckms)**4.0d0*(obj%omdmomp+obj%ombomp)**2.0d0/(obj%homp)**3.0d0*(1.0d0+zlens)**2.0d0)
        else
        sjclscrossobjnoweight = 1.0d0/(distz**2.0d0)/hubblez*deltafinal/(obj%homp**3.0d0)
        end if

    END function sjclscrossobjnoweight


    !LENSING-GALAXY INTEGRAND (Gg)
    function sjclscrossobjonlyweight(obj,zlens)
        type(my_type) :: obj
        REAL(mcp) :: zlens,sjclscrossobjonlyweight,weight2,hubblez
        real(mcp) :: kminsj,kmaxsj

        hubblez = Hofz(zlens)

        if(obj%uselensomp == 0) then
            weight2 = plensobjcubiclowz(obj,zlens) !single spec bin
        end if
        if(obj%uselensomp == 1) then
            weight2 = plensobjcubiccmass(obj,zlens) !single spec bin
        end if
        if(obj%uselensomp == 2) then
            weight2 = plensobjcubic2dfloz(obj,zlens) !single spec bin
        end if
        if(obj%uselensomp == 3) then
            weight2 = plensobjcubic2dfhiz(obj,zlens) !single spec bin
        end if
      	sjclscrossobjonlyweight = weight2*obj%bbiazconstomp*hubblez !constant bias within lens-bin

    END function sjclscrossobjonlyweight


    !LENSING-GALAXY INTEGRAND (Gg)
    function sjclscrossobj(obj,zlens)
        type(my_type) :: obj
        REAL(mcp) :: zlens,sjclscrossobj,weightpart,weight1,weight2,distz,kval,deltafinal,hubblez,ckms
        real(mcp) :: kminsj,kmaxsj

        ckms = 299792.458d0
        distz = ComovingRadialDistance(zlens)
        obj%distlensflatomp = distz
        distz = f_K(distz) !with curvature

        kval = (obj%ellarromp(obj%ellgenomp)+1.0d0/2.0d0)/((obj%h0omp/100.0d0)*distz)

        kminsj = 1.0d-5
        kmaxsj = 100.0d0
        deltafinal = 0.0d0
        
        if (obj%eftflag == 2 .and. (obj%designereftmodel == 4 .or. obj%designereftmodel == 5)) then  
        if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = obj%theoryomp%NL_MPK_WEYL%PowerAt(kval,zlens)
        else
        if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = obj%theoryomp%NL_MPK%PowerAt(kval,zlens)
        end if
        
        weightpart = 3.0d0/2.0d0*(obj%omdmomp+obj%ombomp)*(obj%h0omp/ckms)**2.0d0*distz*(1.0d0+zlens)
        obj%wchooseomp = 1 !doesn't matter here because m1 = m2
        weight1 = weightpart*rombint_obj(obj,weightobjcubic,zlens,obj%mumaxwomp,obj%tol_erroromp) !note mumaxw

        if(obj%uselensomp == 0) then
            weight2 = plensobjcubiclowz(obj,zlens) !single spec bin
        end if
        if(obj%uselensomp == 1) then
            weight2 = plensobjcubiccmass(obj,zlens) !single spec bin
        end if
        if(obj%uselensomp == 2) then
            weight2 = plensobjcubic2dfloz(obj,zlens) !single spec bin
        end if
        if(obj%uselensomp == 3) then
            weight2 = plensobjcubic2dfhiz(obj,zlens) !single spec bin
        end if

        !the hubble term that should be multiplying plensobjcubic cancels the 1/hubble term that should be below
        if (obj%eftflag == 2 .and. (obj%designereftmodel == 4 .or. obj%designereftmodel == 5 )) then 
        sjclscrossobj = weight1*weight2/(distz**2.0d0)*deltafinal/(obj%homp**3.0d0)*obj%bbiazconstomp/(9.0d0/4.0d0*(obj%h0omp/ckms)**4.0d0*(obj%omdmomp+obj%ombomp)**2.0d0/(obj%homp)**3.0d0*(1.0d0+zlens)**2.0d0) !constant bias within lens-bin
        else
        sjclscrossobj = weight1*weight2/(distz**2.0d0)*deltafinal/(obj%homp**3.0d0)*obj%bbiazconstomp !constant bias within lens-bin
        end if
    END function sjclscrossobj


    !LENSING-GALAXY INTEGRAND (Gg) ---- wrt scale factor instead
    function sjclscrossobjsf(obj,alens)
        type(my_type) :: obj
        REAL(mcp) :: alens,zlens,sjclscrossobjsf,weightpart,weight1,weight2,distz,kval,deltafinal,hubblez,ckms
        real(mcp) :: kminsj,kmaxsj

        zlens = 1.0d0/alens - 1.0d0
        ckms = 299792.458d0
        distz = ComovingRadialDistance(zlens)
        obj%distlensflatomp = distz
        distz = f_K(distz) !with curvature

        kval = (obj%ellarromp(obj%ellgenomp)+1.0d0/2.0d0)/((obj%h0omp/100.0d0)*distz)

        kminsj = 1.0d-5
        kmaxsj = 100.0d0
        deltafinal = 0.0d0
        if (obj%eftflag == 2 .and. (obj%designereftmodel == 4 .or. obj%designereftmodel == 5)) then
        if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = obj%theoryomp%NL_MPK_WEYL%PowerAt(kval,zlens)
        else 
        if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = obj%theoryomp%NL_MPK%PowerAt(kval,zlens)
        end if

        weightpart = 3.0d0/2.0d0*(obj%omdmomp+obj%ombomp)*(obj%h0omp/ckms)**2.0d0*distz*(1.0d0+zlens)
        obj%wchooseomp = 1 !doesn't matter here because m1 = m2
        weight1 = weightpart*rombint_obj(obj,weightobjcubic,zlens,obj%mumaxwomp,obj%tol_erroromp) !note mumaxw

        if(obj%uselensomp == 0) then
            weight2 = plensobjcubiclowz(obj,zlens) !single spec bin
        end if
        if(obj%uselensomp == 1) then
            weight2 = plensobjcubiccmass(obj,zlens) !single spec bin
        end if
        if(obj%uselensomp == 2) then
            weight2 = plensobjcubic2dfloz(obj,zlens) !single spec bin
        end if
        if(obj%uselensomp == 3) then
            weight2 = plensobjcubic2dfhiz(obj,zlens) !single spec bin
        end if

        if (obj%eftflag == 2 .and. (obj%designereftmodel == 4 .or. obj%designereftmodel == 5)) then
        sjclscrossobjsf = ((1.0d0+zlens)**2.0d0)*weight1*weight2/(distz**2.0d0)*deltafinal/(obj%homp**3.0d0)*obj%bbiazconstomp/(9.0d0/4.0d0*(obj%h0omp/ckms)**4.0d0*(obj%omdmomp+obj%ombomp)**2.0d0/(obj%homp)**3.0d0*(1.0d0+zlens)**2.0d0)
        else
        !the hubble term that should be multiplying plensobjcubic cancels the 1/hubble term that should be below
        sjclscrossobjsf = ((1.0d0+zlens)**2.0d0)*weight1*weight2/(distz**2.0d0)*deltafinal/(obj%homp**3.0d0)*obj%bbiazconstomp !constant bias within lens bin
        end if
    END function sjclscrossobjsf


     !INTRINSIC-GALAXY INTEGRAND (Ig) ---- wrt scale factor instead
     function sjclscrossgiobjsf(obj,alens)
         type(my_type) :: obj
         REAL(mcp) :: alens,zlens,sjclscrossgiobjsf,weightpart,weight1,weight2,distz,kval,deltafinal,hubblez,ckms,growthnorm
         real(mcp) :: kminsj,kmaxsj
 
         zlens = 1.0d0/alens - 1.0d0
         ckms = 299792.458d0
         distz = ComovingRadialDistance(zlens)
         obj%distlensflatomp = distz
         distz = f_K(distz) !with curvature
 
         kval = (obj%ellarromp(obj%ellgenomp)+1.0d0/2.0d0)/((obj%h0omp/100.0d0)*distz)
 
         kminsj = 1.0d-5
         kmaxsj = 100.0d0
         deltafinal = 0.0d0
         
         if (obj%eftflag == 2 .and. (obj%designereftmodel == 4 .or. obj%designereftmodel == 5)) then
         if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = obj%theoryomp%NL_MPK_WEYL%PowerAt(kval,zlens)
         else
         if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = obj%theoryomp%NL_MPK%PowerAt(kval,zlens)
         end if
 
         growthnorm = obj%gcubicomp%Value(zlens) !this is D(0)/D(z)
         deltafinal = deltafinal*(-obj%ampiayesomp*5.0d-14*2.77536627d11*(obj%omdmomp+obj%ombomp)*growthnorm*((1.0d0+zlens)/(1.0d0+0.3d0))**obj%redziayesomp)*(obj%lumarromp(obj%momp(1)))**obj%lumiayesomp
         hubblez = Hofz(zlens)
 
         obj%wchooseomp = 1 !doesn't matter here because m1 = m2
         weight1 = psourceobjcubic(obj,zlens)*hubblez
 
         if(obj%uselensomp == 0) then
             weight2 = plensobjcubiclowz(obj,zlens) !single spec bin
         end if
         if(obj%uselensomp == 1) then
             weight2 = plensobjcubiccmass(obj,zlens) !single spec bin
         end if
         if(obj%uselensomp == 2) then
             weight2 = plensobjcubic2dfloz(obj,zlens) !single spec bin
         end if
         if(obj%uselensomp == 3) then
             weight2 = plensobjcubic2dfhiz(obj,zlens) !single spec bin
         end if

        if (obj%eftflag == 2 .and. (obj%designereftmodel == 4 .or. obj%designereftmodel == 5 )) then
        sjclscrossgiobjsf = ((1.0d0+zlens)**2.0d0)*weight1*weight2/(distz**2.0d0)*deltafinal/(obj%homp**3.0d0)*obj%bbiazconstomp/(9.0d0/4.0d0*(obj%h0omp/ckms)**4.0d0*(obj%omdmomp+obj%ombomp)**2.0d0/(obj%homp)**3.0d0*(1.0d0+zlens)**2.0d0)
        else 
        !the hubble term that should be multiplying plensobjcubic cancels the 1/hubble term that should be below
         sjclscrossgiobjsf = ((1.0d0+zlens)**2.0d0)*weight1*weight2/(distz**2.0d0)*deltafinal/(obj%homp**3.0d0)*obj%bbiazconstomp !constant bias within lens bin
        end if
     END function sjclscrossgiobjsf


    !lensing weight, cubic spline of source distribution
    function weightobjcubic(obj,zs)
        type(my_type) :: obj
        REAL(mcp) :: zs,weightobjcubic,chis,diff,ckms

        chis = ComovingRadialDistance(zs)
        diff = f_K(chis-obj%distlensflatomp) !with curvature

        weightobjcubic=0.0d0
        if((zs-obj%aphotarromp(obj%momp(obj%wchooseomp))) >= obj%mumin2vomp)  weightobjcubic = obj%psourcetypeomp(obj%momp(obj%wchooseomp))%Value(zs-obj%aphotarromp(obj%momp(obj%wchooseomp)))
        weightobjcubic = weightobjcubic*diff/f_K(chis)

    END function weightobjcubic


    !cubic spline of source distribution
    function psourceobjcubic(obj,zs)
        type(my_type) :: obj
        REAL(mcp) :: zs,psourceobjcubic

        psourceobjcubic=0.0d0
        if((zs-obj%aphotarromp(obj%momp(obj%wchooseomp))) >= obj%mumin2vomp)  psourceobjcubic = obj%psourcetypeomp(obj%momp(obj%wchooseomp))%Value(zs-obj%aphotarromp(obj%momp(obj%wchooseomp)))

    END function psourceobjcubic


    !cubic spline of lowz lens distribution
    function plensobjcubiclowz(obj,zs)
        type(my_type) :: obj
        REAL(mcp) :: zs,plensobjcubiclowz

        plensobjcubiclowz = obj%plenstypeomplowz%Value(zs)

    END function plensobjcubiclowz

    !cubic spline of cmass lens distribution
    function plensobjcubiccmass(obj,zs)
        type(my_type) :: obj
        REAL(mcp) :: zs,plensobjcubiccmass

        plensobjcubiccmass = obj%plenstypeompcmass%Value(zs)

    END function plensobjcubiccmass

    !cubic spline of 2dfloz lens distribution
    function plensobjcubic2dfloz(obj,zs)
        type(my_type) :: obj
        REAL(mcp) :: zs,plensobjcubic2dfloz

        plensobjcubic2dfloz = obj%plenstypeomp2dfloz%Value(zs)

    END function plensobjcubic2dfloz

    !cubic spline of cmass lens distribution
    function plensobjcubic2dfhiz(obj,zs)
        type(my_type) :: obj
        REAL(mcp) :: zs,plensobjcubic2dfhiz

        plensobjcubic2dfhiz = obj%plenstypeomp2dfhiz%Value(zs)

    END function plensobjcubic2dfhiz


    END function CosmoJBD_LnLike


    !----------------------------------------------------------------------------------------------------------
    ! functions for calculating  2D (linear) power spectrum in redshift space with AP distortions:
    !
    ! 	 see W. E. Ballinger 1996 (http://arxiv.org/pdf/astro-ph/9605017v1.pdf) for derivation or Beulter 2014 for formula used
    ! We perform the scaling from observed to true k and mu, before calcualting the 2D power spectrum.
    !
    !----------------------------------------------------------------------------------------------------------

    function LinearPower_2D_s_cmass(Theory,k_obs,mu_obs,z,sigma_v_cmass,a_perp_cmass,a_par_cmass,b1_cmass,N_shot_cmass)
        implicit none
        class(TCosmoTheoryPredictions) :: Theory
        real(mcp),   intent(in) :: z
        REAL(mcp), intent(in) :: k_obs,mu_obs,b1_cmass
        REAL(mcp), intent(in) :: a_par_cmass,N_shot_cmass,sigma_v_cmass,a_perp_cmass
        REAL(mcp) :: LinearPower_2D_s_cmass,F_AP,beta,growth_rate
        REAL(mcp) :: mu_TRUE,k_TRUE,scaling_factor,FoG_term

        F_AP = a_par_cmass/a_perp_cmass
        scaling_factor = 1.0 + (mu_obs**2)*(1.0/F_AP**2 - 1.0)
        mu_TRUE = (mu_obs/F_AP)*(scaling_factor)**(-0.5)
        k_TRUE =  (k_obs/a_perp_cmass)*(scaling_factor)**(0.5)
    
        !Modified for EFT start
        if (CosmoSettings%eftflag == 2 .and. (CosmoSettings%designereftmodel == 4 .or. CosmoSettings%designereftmodel == 5)) then
           growth_rate = Theory%growth_k_z%PowerAt(k_TRUE,z)
        else
           growth_rate = Theory%growth_z%Value(z)/Theory%sigma8_z%Value(z)
           ! AL uses growth = f\sigma8, not f
        end if
        !Modified for EFT end
        beta =  growth_rate/b1_cmass
        FoG_term = exp(-(k_TRUE*mu_TRUE*sigma_v_cmass)**2)

        LinearPower_2D_s_cmass = &
        b1_cmass**2*(Theory%MPK%PowerAt(k_TRUE,z) + N_shot_cmass)*(1.0 + beta*mu_TRUE**2)**2*FoG_term

    END function LinearPower_2D_s_cmass

    function LinearPower_2D_s_wigz(Theory,k_obs,mu_obs,z_wig,sig_v,a_perp,a_par,b1_bias)
        implicit none
        class(TCosmoTheoryPredictions) :: Theory
        REAL(mcp), intent(in) :: k_obs,mu_obs,b1_bias,z_wig
        REAL(mcp), intent(in) :: a_par,sig_v,a_perp
        REAL(mcp) :: LinearPower_2D_s_wigz,F_AP,beta,growth_rate
        REAL(mcp) :: mu_TRUE,k_TRUE,scaling_factor,FoG_term

        F_AP = a_par/a_perp
        scaling_factor = (1.0 + (mu_obs**2)*(1.0/(F_AP**2.0) - 1.0))
        mu_TRUE = (mu_obs/F_AP)*(scaling_factor)**(-0.5)
        k_TRUE =  (k_obs/a_perp)*(scaling_factor)**(0.5)
        !Modified for EFT start
        if (CosmoSettings%eftflag == 2 .and. (CosmoSettings%designereftmodel == 4 .or. CosmoSettings%designereftmodel == 5)) then
           growth_rate = Theory%growth_k_z%PowerAt(k_TRUE,z_wig)
        else
           growth_rate = Theory%growth_z%Value(z_wig)/Theory%sigma8_z%Value(z_wig)
           !AL uses growth = f\sigma8, not f
        end if
        !Modified for EFT end
        beta = growth_rate/b1_bias
        FoG_term = exp(-(k_TRUE*mu_TRUE*sig_v*growth_rate)**2)

        LinearPower_2D_s_wigz = (b1_bias**2)*Theory%MPK%PowerAt(k_TRUE,z_wig)*((1.0 + beta*mu_TRUE**2)**2)*FoG_term

    END function LinearPower_2D_s_wigz

    !----------------------------------------------------------------------------------------------------------
    ! Legendre polynomials needed
    !----------------------------------------------------------------------------------------------------------

    function Legendre_Pn(mu,n)
    implicit none
    integer :: n
    real(mcp) :: Legendre_Pn,mu

    if(n .EQ. 0) Legendre_Pn = 1.0
    if(n .EQ. 2) Legendre_Pn = (1.0/2.0)*(3.0*mu**2 - 1.0)
    if(n .EQ. 4) Legendre_Pn = (1.0/8.0)*(35.0*mu**4 - 30.0*mu**2 + 3.0)

    end function Legendre_Pn

    !----------------------------------------------------------------------------------------------------------
    ! Functions needed for AP calculation
    !----------------------------------------------------------------------------------------------------------

    function D_AhUnit(z,CMB)
        class(CMBParams) :: CMB
        real(mcp) :: z,D_AhUnit

        D_AhUnit = AngularDiameterDistance(z)*(CMB%h0/100)

    end function D_AhUnit

    function HofzhUnit(z,CMB,c)
        class(CMBParams) :: CMB
        real(mcp) :: z,HofzhUnit,c

        HofzhUnit = (c*Hofz(z)/1.d3)/(CMB%h0/100)

    end function HofzhUnit

    function numcat_local(S, num)
        character(LEN=*) S
        character(LEN=1024) numcat_local, numstr
        integer num

        write (numstr, *) num
        numcat_local = trim(S) // trim(adjustl(numstr))

    end function numcat_local

    end module CosmoJBD
