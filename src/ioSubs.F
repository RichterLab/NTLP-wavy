

      subroutine read_input_file
      use pars
      use particles
      use con_data
      use tank_pars
      implicit none

      character(48) :: label
      character(180) :: params_dir
      namelist /step_params/ iti,itmax,imean,ihst,itape,
     +                       itstr,it_his,it_viz,iviz,itn

      namelist /grid_params/ ncpu_s, Uo, Ttop, Tbot,
     +         qstar, tsfcc, ugcont, vgcont,
     +         zi, zl, xl, yl, zw1,
     +         zo, ustar_tank,waveslope,c_speed,wave_l,minwavez,
     +         dt_new,cfl

      namelist /path_names/ path_seed,path_part,path_res,
     +         path_ran

      namelist /flags/ iradup,ifilt,itcut,isubs,ibrcl,iocean,idebug,
     +         ivis0,new_vis,ismlt,ifree,iupwnd,iskew,ibuoy,ifix_dt,
     +         i_move_step,method,iz_space,iflat,i3d_surf,i_pm,lpbc,
     +         k_stab,no_pgrad,ispray,iDNS,ifields,ievap,isfc,i_tank,
     +         i_zeta,inetcdf,i_pinertia

      namelist /constants/ grav, t00, fcor, tnumpart,
     +         radius_init,part_grav,Vpmax,nuf,Pra,Sc,
     +         vp_init,Tp_init,rhoa,rhow,
     +         Cpp,Cpa,Mw,Ru,Ms,Sal,Gam,Ion,Os


      !params.in contains namelists to read
      !open(12, file="./params.in", status="old")
      call get_command_argument(1,params_dir)
      open(12,file=params_dir,status="old")

      read(12,nml=step_params)
      if (myid==0) print step_params

      read(12,nml=flags)
      if (myid==0) print flags

      read(12,nml=grid_params)
      if (myid==0) print grid_params

      read(12,nml=path_names)
      if (myid==0) print path_names

      read(12,nml=constants)
      if (myid==0) print constants
      CpaCpp = Cpa/Cpp

      dz_U = (zl-minwavez)/maxnz


      end subroutine read_input_file

c       =========== WRITE HISTOGRAMS ==========
        subroutine write_histograms(istep)
        use pars
        use fields
        use con_data
        use con_stats
        use particles
        implicit none

        include 'mpif.h'
        integer :: i,j
        integer :: ierr


        integer,intent(in) :: istep
        integer :: iblnk
        integer :: num_entries

        real :: sumbuf_restime(histbins+2,histbins+2)
        real :: sumbuf_locations(histbins+2,histbins+2)
        real :: sumbuf_height(histbins+2,histbins+2)
        real :: sumbuf_momentum(histbins+2,histbins+2)


        character*80 path_restime,path_locations,path_height
        character*80 path_momentum

        character cgrid*10, num*3,iblks*16, fformat*10



c -------------- Collect all histogram from different processors
        num_entries = (histbins+2)*(histbins+2)
        call mpi_reduce(hist_restime,sumbuf_restime,num_entries,
     +   mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
        call mpi_reduce(hist_locations,sumbuf_locations,num_entries,
     +   mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
        call mpi_reduce(hist_height,sumbuf_height,num_entries,
     +   mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
        call mpi_reduce(hist_momentum,sumbuf_momentum,num_entries,
     +   mpi_real8,mpi_sum,0,mpi_comm_world,ierr)


        if (myid==0) then

c
c --------------- build character strings for file name
c
        cgrid = '.le.'
        write(num,'(i3.3)') itn

        iblnk = index(path_histog,' ')
        call blnk(iblks)
        write(iblks(1:7),'(i7.7)') istep
        iblks(8:8) = '_'
        write(iblks(9:15),'(i7.7)') (istep+itape)
        call blnk(path_sav_hist)


        path_restime = path_histog(1:iblnk-1)//'/restime'//
     +          cgrid(1:4)//case(1:3)//'.'//iblks(1:15)//'.txt'
        path_locations = path_histog(1:iblnk-1)//'/locations'//
     +          cgrid(1:4)//case(1:3)//'.'//iblks(1:15)//'.txt'
        path_height = path_histog(1:iblnk-1)//'/height'//
     +          cgrid(1:4)//case(1:3)//'.'//iblks(1:15)//'.txt'
        path_momentum = path_histog(1:iblnk-1)//'/momentum'//
     +          cgrid(1:4)//case(1:3)//'.'//iblks(1:15)//'.txt'

c ---------------- save data

        write(num,'(i3.3)') histbins+2
        fformat = '('//num//'e15.6)'
        open (unit=900,file=path_restime,
     +  form='formatted',access='sequential')

        write(900,fformat),bins(:,1)
        write(900,fformat),bins(:,2)

        do i =1,histbins+2
           write(900,fformat),sumbuf_restime(1:histbins+2,i)
        end do


        close(900)

        open (unit=901,file=path_locations,
     +  form='formatted',access='sequential')

        write(901,fformat),bins(:,2)
        write(901,fformat),bins(:,3)

        do i =1,histbins+2
           write(901,fformat),sumbuf_locations(1:histbins+2,i)
        end do


        close(901)


        open (unit=902,file=path_height,
     +  form='formatted',access='sequential')

        write(902,fformat),bins(:,2)
        write(902,fformat),bins(:,4)

        do i =1,histbins+2
           write(902,fformat),sumbuf_height(1:histbins+2,i)
        end do


        close(902)


        open (unit=903,file=path_momentum,
     +  form='formatted',access='sequential')

        write(903,fformat),bins(:,2)
        write(903,fformat),bins(:,5)

        do i =1,histbins+2
           write(903,fformat),sumbuf_momentum(1:histbins+2,i)
        end do


        close(903)
        end if

        end subroutine write_histograms
        subroutine set_binsdata(binsdata,sizea,numvars,lmin,lmax,idx)
        use pars
        use fields
        use con_data
        use con_stats
        use particles
        implicit none

        integer :: i,nbin,ibin,nbinnew,idx
        real :: dhl

        integer,intent(in) :: sizea,numvars
        real,intent(inout) :: binsdata(sizea,numvars)

        real :: lmin,lmax

        nbin = histbins !From Module Particle


c       Calculate size of interval
        dhl = (lmax-lmin)/nbin

c       ===== update x-axis for each bin =====
        binsdata(1,idx) = lmin-dhl
        do i = 1,histbins+1
          binsdata(i+1,idx)= dhl+binsdata(i,idx)
        end do

        end subroutine set_binsdata
        subroutine add_histogram(binsdata1,binsdata2,
     +            histdata,sizea,val1,val2)

        use pars
        use fields
        use con_data
        use con_stats
        use particles
        implicit none

        integer :: i,nbin,nbinnew
        integer :: ibin(2)
        real :: dht,dhl
        real,intent(in) :: val1,val2
c       real,intent(inout),dimension(:,:) :: histdata

        integer,intent(in) :: sizea
        real,intent(in) :: binsdata1(sizea),binsdata2(sizea)
        real,intent(inout) :: histdata(sizea,sizea)

        real :: tmin,lmin,tmax,lmax

        tmin = binsdata1(2)
        lmin = binsdata2(2)
        nbin = histbins !From Module Particle
        tmax = binsdata1(nbin+1)
        lmax = binsdata2(nbin+1)

c       Calculate size of interval
        dht = (tmax-tmin)/nbin
        dhl = (lmax-lmin)/nbin
c       Find corresponding bin of input data


        if (val1 .gt. tmax) then
                ibin(1) = nbin + 2
        else if (val1 .lt. tmin) then
                ibin(1) = 1
        else
                ibin(1) = (floor((val1-tmin)/dht)+1)+1
        end if

        if (val2 .gt. lmax) then
                ibin(2) = nbin + 2
        else if (val2 .lt. lmin) then
                ibin(2) = 1
        else
                ibin(2) = (floor((val2-lmin)/dhl)+1)+1
        end if


c       Add the current event to the histogram
        histdata(ibin(1),ibin(2)) = histdata(ibin(1),ibin(2)) + 1.0

        end subroutine add_histogram
        subroutine clean_histograms
        use pars
        use fields
        use con_data
        use con_stats
        use particles
        hist_restime = 0
        hist_locations = 0
        hist_height = 0
        hist_momentum = 0

        end subroutine clean_histograms
