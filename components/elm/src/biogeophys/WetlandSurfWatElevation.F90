module WetlandSurfWatElevation

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! CLM Satelitte Phenology model (SP) ecosystem dynamics (phenology, vegetation).
  ! Allow some subroutines to be used by the CLM Carbon Nitrogen model (CLMCN)
  ! so that DryDeposition code can get estimates of LAI differences between months.
  !
  ! !USES:
  use shr_strdata_mod , only : shr_strdata_type, shr_strdata_create
  use shr_strdata_mod , only : shr_strdata_print, shr_strdata_advance
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_kind_mod    , only : CL => shr_kind_CL
  use shr_kind_mod    , only : CX => shr_kind_CXX
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun
  use elm_varctl      , only : iulog, scmlat,scmlon,single_column
  use elm_varcon      , only : grlnd
  use controlMod      , only : NLFilename
  use decompMod       , only : gsmap_lnd_gdc2glo
  use domainMod       , only : ldomain
  use fileutils       , only : getavu, relavu
  use VegetationType       , only : veg_pp
  use WaterstateType  , only : waterstate_type
  use perf_mod        , only : t_startf, t_stopf
  use spmdMod         , only : masterproc
  use spmdMod         , only : mpicom, comp_id
  use mct_mod
  use ncdio_pio   
  use topounit_varcon , only : max_topounits
  use GridcellType    , only : grc_pp
  use ColumnDataType    , only : column_water_state
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: interpDailyWetSurfWatElev       ! interpolate monthly vegetation data
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  private :: readDailyWetSurfWatElev  ! read monthly vegetation data for two months
  !
  ! !PRIVATE TYPES:

  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine interpDailyWetSurfWatElev (bounds, col_ws)
    !
    ! !DESCRIPTION:
    ! Determine if 2 new months of data are to be read.
    !
    ! !USES:
    use elm_varctl      , only : fsurdat
    use clm_time_manager, only : get_curr_date, get_step_size, get_nstep,get_curr_calday
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    type(column_water_state), intent(inout) :: col_ws
    !
    ! !LOCAL VARIABLES:
    !real(r8):: dtime       ! land model time step (sec)

    !-----------------------------------------------------------------------
    
    !dtime = get_step_size()

    call readDailyWetSurfWatElev (bounds, fsurdat, col_ws)

  end subroutine interpDailyWetSurfWatElev
  
  !-----------------------------------------------------------------------
  
  subroutine readDailyWetSurfWatElev (bounds, &
       fveg, col_ws)
    !
    ! !DESCRIPTION:
    ! Read monthly vegetation data for two consec. months.
    !
    ! !USES:
    use elm_varpar       , only : numpft
    use pftvarcon        , only : noveg
    use fileutils        , only : getfil
    use spmdMod          , only : masterproc, mpicom, MPI_REAL8, MPI_INTEGER
    use shr_scam_mod     , only : shr_scam_getCloseLatLon
    use clm_time_manager , only : get_nstep,get_curr_calday,get_curr_date
    use LandunitType   , only: lun_pp
    use ColumnType   , only: col_pp
    use landunit_varcon   , only: istwet
    use netcdf
    
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in) :: bounds
    character(len=*)  , intent(in) :: fveg      ! file with monthly vegetation data
    type(column_water_state), intent(inout) :: col_ws
    !
    ! !LOCAL VARIABLES:
    character(len=256) :: locfn           ! local file name
    type(file_desc_t)  :: ncid            ! netcdf id
    integer :: g,p,c     ! indices
    integer :: dimid,varid                ! input netCDF id's
    integer :: ntim                       ! number of input data time samples
    integer :: nlon_i                     ! number of input data longitudes
    integer :: nlat_i                     ! number of input data latitudes
    integer :: npft_i                     ! number of input data pft types
    integer :: ier                        ! error code
    integer :: closelatidx,closelonidx
    real(r8):: closelat,closelon
    logical :: readvar
    real(r8), pointer :: SurfWatElev(:,:,:)        ! lai read from input files
    character(len=32) :: subname = 'readDailyWetSurfWatElev'
    integer :: jday        ! julian day of the year
    integer :: year,mon,day,sec    !Date components for end of current time step
    
    !-----------------------------------------------------------------------

    ! Determine necessary indices

    jday =  get_curr_calday()
    call get_curr_date(year,mon,day,sec)
    
    allocate(&
         SurfWatElev(bounds%begg:bounds%endg,0:numpft,0:364), &
         stat=ier)
    if (ier /= 0) then
       write(iulog,*)subname, 'allocation big error '
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end if

    ! ----------------------------------------------------------------------
    ! Open Surface Water Elevation file
    ! Read data and convert from gridcell to pft data
    ! ----------------------------------------------------------------------

    call getfil(fveg, locfn, 0)
    call ncd_pio_openfile (ncid, trim(locfn), 0)

    if (single_column) then
       call shr_scam_getCloseLatLon (ncid, scmlat, scmlon, closelat, closelon,&
            closelatidx, closelonidx)
    endif

    if (year == 2012) then
       call ncd_io(ncid=ncid, varname='SurfWatElev_Wet_2012', flag='read', data=SurfWatElev, dim1name=grlnd, &
            readvar=readvar)
       if (.not. readvar) call endrun(msg=' ERROR: SurfWatElev_Wet NOT on fveg file'//errMsg(__FILE__, __LINE__))
    elseif (year == 2013) then
       call ncd_io(ncid=ncid, varname='SurfWatElev_Wet_2013', flag='read', data=SurfWatElev, dim1name=grlnd, &
            readvar=readvar)
       if (.not. readvar) call endrun(msg=' ERROR: SurfWatElev_Wet NOT on fveg file'//errMsg(__FILE__, __LINE__))
    elseif (year == 2014) then
       call ncd_io(ncid=ncid, varname='SurfWatElev_Wet_2014', flag='read', data=SurfWatElev, dim1name=grlnd, &
            readvar=readvar)
       if (.not. readvar) call endrun(msg=' ERROR: SurfWatElev_Wet NOT on fveg file'//errMsg(__FILE__, __LINE__))
    elseif (year == 2015) then
       call ncd_io(ncid=ncid, varname='SurfWatElev_Wet_2015', flag='read', data=SurfWatElev, dim1name=grlnd, &
            readvar=readvar)
       if (.not. readvar) call endrun(msg=' ERROR: SurfWatElev_Wet NOT on fveg file'//errMsg(__FILE__, __LINE__))
    elseif (year == 2016) then
       call ncd_io(ncid=ncid, varname='SurfWatElev_Wet_2016', flag='read', data=SurfWatElev, dim1name=grlnd, &
            readvar=readvar)
       if (.not. readvar) call endrun(msg=' ERROR: SurfWatElev_Wet NOT on fveg file'//errMsg(__FILE__, __LINE__))
    elseif (year == 2017) then
       call ncd_io(ncid=ncid, varname='SurfWatElev_Wet_2017', flag='read', data=SurfWatElev, dim1name=grlnd, &
            readvar=readvar)
       if (.not. readvar) call endrun(msg=' ERROR: SurfWatElev_Wet NOT on fveg file'//errMsg(__FILE__, __LINE__))
    elseif (year == 2018) then
       call ncd_io(ncid=ncid, varname='SurfWatElev_Wet_2018', flag='read', data=SurfWatElev, dim1name=grlnd, &
            readvar=readvar)
       if (.not. readvar) call endrun(msg=' ERROR: SurfWatElev_Wet NOT on fveg file'//errMsg(__FILE__, __LINE__))
    elseif (year == 2019) then
       call ncd_io(ncid=ncid, varname='SurfWatElev_Wet_2019', flag='read', data=SurfWatElev, dim1name=grlnd, &
            readvar=readvar)
       if (.not. readvar) call endrun(msg=' ERROR: SurfWatElev_Wet NOT on fveg file'//errMsg(__FILE__, __LINE__))
    elseif (year == 2020) then
       call ncd_io(ncid=ncid, varname='SurfWatElev_Wet_2020', flag='read', data=SurfWatElev, dim1name=grlnd, &
            readvar=readvar)
       if (.not. readvar) call endrun(msg=' ERROR: SurfWatElev_Wet NOT on fveg file'//errMsg(__FILE__, __LINE__))
    elseif (year == 2021) then
       call ncd_io(ncid=ncid, varname='SurfWatElev_Wet_2021', flag='read', data=SurfWatElev, dim1name=grlnd, &
            readvar=readvar)
       if (.not. readvar) call endrun(msg=' ERROR: SurfWatElev_Wet NOT on fveg file'//errMsg(__FILE__, __LINE__))
    elseif (year == 2022) then
       call ncd_io(ncid=ncid, varname='SurfWatElev_Wet_2022', flag='read', data=SurfWatElev, dim1name=grlnd, &
            readvar=readvar)
       if (.not. readvar) call endrun(msg=' ERROR: SurfWatElev_Wet NOT on fveg file'//errMsg(__FILE__, __LINE__))
    else
       call ncd_io(ncid=ncid, varname='SurfWatElev_Wet', flag='read', data=SurfWatElev, dim1name=grlnd, &
            readvar=readvar)
       if (.not. readvar) call endrun(msg=' ERROR: SurfWatElev_Wet NOT on fveg file'//errMsg(__FILE__, __LINE__))
    endif

       ! Assign surface water depth to each pft
       do p = bounds%begp,bounds%endp
          g =veg_pp%gridcell(p)
          c = veg_pp%column(p)

          if (lun_pp%itype(col_pp%landunit(c)) == istwet) then
             if (SurfWatElev(g,p,jday-1) <= 1200) then
                col_ws%h2osfc_wet(c) = SurfWatElev(g,13,jday-1)
             else
                col_ws%h2osfc_wet(c) = 1200_r8
             endif
          else
             col_ws%h2osfc_wet(c) = 0_r8
          endif
                       
       end do   ! end of loop over patches


    call ncd_pio_closefile(ncid)

    deallocate(SurfWatElev)

  end subroutine readDailyWetSurfWatElev

end module WetlandSurfWatElevation
