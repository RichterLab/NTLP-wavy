

      subroutine read_input_file
        use pars
        use particles
        use con_data
        use tank_pars
        implicit none

        character(48) :: label
        character(180) :: params_dir
        namelist /step_params/ iti,itmax,imean,ihst,itape,
     +      itstr,it_his,it_viz,iviz,itn

        namelist /io_params/ inetcdf,ifields,itrajout,
     +      path_seed,path_part,
     +      path_res,path_ran,wave_file

        namelist /computational_params/ iradup,ifilt,itcut,isubs,ibrcl,
     +      iocean,idebug,ivis0,new_vis,ismlt,ifree,iupwnd,iskew,ibuoy,
     +      ifix_dt,move_step,method,iz_space,iflat,i3d_surf,i_pm,
     +      lpbc,k_stab,no_pgrad,ispray,iDNS,ievap,isfc,
     +      i_zeta,pinertia,verb,ncpu_s,minwavez,zi,zw1,dt_new,cfl

        namelist /physical_params/ nuf,Pra,Sc,Cpa,ustar_tank,ugcont,
     +      vgcont,t00,fcor,grav,zl,xl,yl,zo,Uo,qstar,tsfcc,Ttop,Tbot,
     +      tnumpart,rhoa,rhow,part_grav,Vpmax,Cpp,Mw,Ru,Ms,Sal,Gam,Ion,
     +      Os,radius_init,Tp_init,vp_init,waveslope,c_speed,wave_l


        !params.in contains namelists to read
        call get_command_argument(1,params_dir)
        open(12,file=params_dir,status="old")

        read(12,nml=step_params)
        if (myid==0) print step_params

        read(12,nml=io_params)
        if (myid==0) print io_params

        read(12,nml=computational_params)
        if (myid==0) print computational_params

        read(12,nml=physical_params)
        if (myid==0) print physical_params
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
        real :: sumbuf_slopes(histbins+2,histbins+2)
        real :: sumbuf_initV(histbins+2,histbins+2)


        character*80 path_restime,path_locations,path_height
        character*80 path_momentum, path_slopes, path_initV

        character cgrid*10, num*3,iblks*16, fformat*10



c -------------- Collect all histogram from different processors
        num_entries = (histbins+2)*(histbins+2)
        call mpi_reduce(hist_restime,sumbuf_restime,num_entries,
     +      mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
        call mpi_reduce(hist_locations,sumbuf_locations,num_entries,
     +      mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
        call mpi_reduce(hist_height,sumbuf_height,num_entries,
     +      mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
        call mpi_reduce(hist_momentum,sumbuf_momentum,num_entries,
     +      mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
        call mpi_reduce(hist_slopes,sumbuf_slopes,num_entries,
     +      mpi_real8,mpi_sum,0,mpi_comm_world,ierr)
        call mpi_reduce(hist_initV,sumbuf_initV,num_entries,
     +      mpi_real8,mpi_sum,0,mpi_comm_world,ierr)


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
     +        cgrid(1:4)//case(1:3)//'.'//iblks(1:15)//'.txt'
          path_locations = path_histog(1:iblnk-1)//'/locations'//
     +        cgrid(1:4)//case(1:3)//'.'//iblks(1:15)//'.txt'
          path_height = path_histog(1:iblnk-1)//'/height'//
     +        cgrid(1:4)//case(1:3)//'.'//iblks(1:15)//'.txt'
          path_momentum = path_histog(1:iblnk-1)//'/momentum'//
     +        cgrid(1:4)//case(1:3)//'.'//iblks(1:15)//'.txt'
          path_slopes = path_histog(1:iblnk-1)//'/slopes'//
     +        cgrid(1:4)//case(1:3)//'.'//iblks(1:15)//'.txt'
          path_initV = path_histog(1:iblnk-1)//'/initV'//
     +        cgrid(1:4)//case(1:3)//'.'//iblks(1:15)//'.txt'

c ---------------- save data

          write(num,'(i3.3)') histbins+2
          fformat = '('//num//'e15.6)'
          open (unit=900,file=path_restime,
     +        form='formatted',access='sequential')

          write(900,fformat),bins(:,1)
          write(900,fformat),bins(:,2)

          do i =1,histbins+2
            write(900,fformat),sumbuf_restime(1:histbins+2,i)
          end do


          close(900)

          open (unit=901,file=path_locations,
     +        form='formatted',access='sequential')

          write(901,fformat),bins(:,2)
          write(901,fformat),bins(:,3)

          do i =1,histbins+2
            write(901,fformat),sumbuf_locations(1:histbins+2,i)
          end do


          close(901)


          open (unit=902,file=path_height,
     +        form='formatted',access='sequential')

          write(902,fformat),bins(:,2)
          write(902,fformat),bins(:,4)

          do i =1,histbins+2
            write(902,fformat),sumbuf_height(1:histbins+2,i)
          end do


          close(902)


          open (unit=903,file=path_momentum,
     +        form='formatted',access='sequential')

          write(903,fformat),bins(:,2)
          write(903,fformat),bins(:,5)

          do i =1,histbins+2
            write(903,fformat),sumbuf_momentum(1:histbins+2,i)
          end do


          close(903)


          open (unit=904,file=path_slopes,
     +        form='formatted',access='sequential')

          write(904,fformat),bins(:,1)
          write(904,fformat),bins(:,6)

          do i =1,histbins+2
            write(904,fformat),sumbuf_slopes(1:histbins+2,i)
          end do


          close(904)


          open (unit=905,file=path_initV,
     +        form='formatted',access='sequential')

          write(905,fformat),bins(:,4)
          write(905,fformat),bins(:,7)

          do i =1,histbins+2
            write(905,fformat),sumbuf_initV(1:histbins+2,i)
          end do


          close(905)
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
     +      histdata,sizea,val1,val2)

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


        if (val1 > tmax) then
          ibin(1) = nbin + 2
        else if (val1 < tmin) then
          ibin(1) = 1
        else
          ibin(1) = (floor((val1-tmin)/dht)+1)+1
        end if

        if (val2 > lmax) then
          ibin(2) = nbin + 2
        else if (val2 < lmin) then
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
        hist_slopes = 0
        hist_initV = 0

      end subroutine clean_histograms
      subroutine save_particles
        use particles
        use pars
        implicit none
        include 'mpif.h'

        integer :: istatus(mpi_status_size), ierr, fh
        integer(kind=mpi_offset_kind) :: zoffset,offset
        integer :: pnum_vec(numprocs)
        integer :: iproc,i
        type(particle) :: writebuf(numpart),tmp

      !Do this with mpi_write_at_all
      !Need to figure out the displacements - need numpart from each proc
        call mpi_allgather(numpart,1,mpi_integer,pnum_vec,1,mpi_integer,
     +      mpi_comm_world,ierr)

      !Package all the particles into writebuf:
        i = 1
        part => first_particle
        do while (associated(part))
          writebuf(i) = part
      !write(*,'a5,3e15.6') 'xp:',part%xp(1:3)
          part => part%next
          i = i + 1
        end do

      !Now only write to the file if you actually have particles
      !EXCEPTION: proc 0, which needs to write tnumpart regardless
        call mpi_file_open(mpi_comm_world, path_sav_part,
     +      mpi_mode_create+mpi_mode_rdwr,
     +      mpi_info_null,fh,ierr)

        zoffset = 0
      !Write tnumpart first:
        if (myid==0) then
          call mpi_file_write_at(fh,zoffset,tnumpart,1,mpi_integer,
     +        istatus,ierr)
          write(*,*) 'wrote tnumpart = ',tnumpart
        end if

        zoffset = zoffset + 4

      !Now compute the offset (in bytes!):
        offset = zoffset
        do iproc = 0,myid-1
          offset = offset + pnum_vec(iproc+1)*(sizeof(tmp)-pad_diff)
        end do

      !Now everyone else write, ONLY if numpart > 0
        if (numpart > 0) then
          call mpi_file_write_at(fh,offset,writebuf,numpart,
     +        particletype,istatus,ierr)
        end if

        call mpi_file_close(fh,ierr)

        write(*,*) 'proc',myid,'wrote numpart = ',numpart

        if (myid==0) write(*,7000) path_sav_part
 7000   format(' PARTICLE DATA IS WRITTEN IN FILE  ',a80)

      end subroutine save_particles
      subroutine read_part_res
        use pars
        use particles
        implicit none
        include 'mpif.h'

        integer :: istatus(mpi_status_size), ierr, fh
        integer(kind=mpi_offset_kind) :: zoffset,offset
        integer :: myp,totalp
        integer :: iproc,i,pidxmax
        type(particle), allocatable :: readbuf(:)

        if (myid==0) write(*,7000) path_part
 7000   format(' READING PARTICLE DATA FROM  ',a80)


        call mpi_file_open(mpi_comm_world,path_part,
     +      mpi_mode_rdonly,
     +      mpi_info_null,fh,ierr)


      !Read in the total number of particles:
        offset = 0
        call mpi_file_read_at_all(fh,offset,tnumpart,1,
     +      mpi_integer,istatus,ierr)
        if (myid==0) write(*,*) 'read tnumpart = ',tnumpart

        offset = 4
        allocate(readbuf(tnumpart))
        call mpi_file_read_at_all(fh,offset,readbuf,tnumpart,
     +      particletype,istatus,ierr)

        do i = 1,tnumpart
        !Now - does it lie within this proc's bounds?
          if (readbuf(i)%xp(2) > ymin .AND.
     +        readbuf(i)%xp(2) < ymax .AND.
     +        readbuf(i)%xp(1) > xmin .AND.
     +        readbuf(i)%xp(1) < xmax) then
            if (.NOT. associated(first_particle)) then
              allocate(first_particle)
              first_particle = readbuf(i)
              nullify(first_particle%prev,first_particle%next)
              part => first_particle
            else
              allocate(part%next)
              part%next = readbuf(i)
              part%next%prev => part
              part => part%next
              nullify(part%next)
            end if

          end if
        end do

        deallocate(readbuf)

        call mpi_file_close(fh,ierr)

      !Now just check how many each processor obtained:
      !At the same time, figure out max(pidx) and set ngidx
      !to one plus this value:
        pidxmax = 0
        part => first_particle
        myp = 0
        do while (associated(part))
          myp = myp+1
          if (part%pidx > pidxmax) pidxmax = part%pidx
          part => part%next
        end do

      !Set ngidx (the index for creating new particles) to 1+pidmax:
        ngidx = pidxmax + 1

        numpart = myp

        call mpi_allreduce(myp,totalp,1,mpi_integer,mpi_sum,
     +      mpi_comm_world,ierr)

        write(*,*) 'proc',myid,'read in numpart:',myp
        if (myid==0) write(*,*) 'total number of particles read:',totalp

      end subroutine read_part_res

      subroutine set_paths
c
c ------------- set file path for RESTART, and
c               directories for saving, history, and viz files
c
c     path_res  --- the path and file name of the velocity restart file
c     path_sav  --- the path where the new 3d volumes are to be saved
c     path_his  --- the path where the new history files are to be saved
c     path_viz  --- the path where xy, xz, or yz planes of data will be stored
c     path_stuf --- the path where fun facts about the viz planes of
c                   data will be stored
c     path_bndy --- the path where the boundary shape is stored
c     path_seed --- seed name
c
        use pars
c
c
c --------- initialization
c
        path_sav    = path_seed
        path_his    = path_seed
        path_viz_xy = path_seed
        path_viz_xz = path_seed
        path_viz_yz = path_seed
        path_stuf   = path_seed
        path_bndy   = path_seed
        path_histog = path_seed
c
        return
      end

      subroutine get_fields
c ----------- special routine to read just 3d fields
c             as an initial guess, easy to customize
c             if missing data, etc..

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
#if defined(SWAP)
        use module_byteswap
#endif
        include 'mpif.h'

        integer status(mpi_status_size), ierr
        integer(kind=mpi_offset_kind) :: offset, disp
        integer(kind=k8)              :: nsize, nsize2
        real, allocatable, dimension(:,:,:) :: temp
        logical there

        allocate(temp(nvar+4,nnx,iys:iye))

c ---------- input file to read from

c --------------------- get restart file from local directory
c                       reuse unit number

        close(nvel)

        inquire(file=path_ran,exist=there)
        if(there) then
          if(l_root) write(6,6001) path_ran
        else
          if(l_root) write(6,6005) path_ran
          stop
        end if

c ---- open file

        call mpi_file_open(mpi_comm_world, path_ran,
     +      mpi_mode_create+mpi_mode_rdwr,
     +      mpi_info_null, nvel, ierr)

c ---- set file view

        disp = 0
        call mpi_file_set_view(nvel,disp,mpi_real8,mpi_real8,
     +      'native',mpi_info_null,ierr)

c ------------ read 3d fields, make rhs*8
c              grid speed is in file but don't need it for this restart

        nsize  = int(nvar+4,k8)*nnx*nny
        nsize2 = int(nvar+4,k8)*nnx*(iys-1)
        n_read = (nvar+4)*nnx*(iye+1-iys)

        do k=izs,ize
          offset = int((k-1),k8)*nsize + nsize2
          call mpi_file_read_at_all(nvel,offset,temp,n_read,
     +        mpi_real8,status,ierr)
          if (ierr /= 0) goto 9992
#if defined(SWAP)
          call byteswap(temp)
#endif
          u_off_set = 0.0
          do j=iys,iye
            do i=1,nnx
              u(i,j,k)  = temp(1,i,j) + u_off_set
              v(i,j,k)  = temp(2,i,j)
              w(i,j,k)  = temp(3,i,j)
              uf(i,j,k) = temp(4,i,j) + u_off_set
              vf(i,j,k) = temp(5,i,j)
              wf(i,j,k) = temp(6,i,j)
              e(i,j,k)  = temp(nvar+3,i,j)
            end do
          end do
          do is = 1,nscl
            do j = iys,iye
              do i = 1,nnx
                t(i,j,is,k) = temp(6+is,i,j)
              end do
            end do
          end do

        end do

c ---- close file

        call mpi_file_close(nvel, ierr)

        deallocate(temp)

        do k=1,nnz
          ug(k) = ugcont
          vg(k) = vgcont
        end do

        do iz=izs,ize

c ---------------- initial guess for pressure

          do iy=iys,iye
            do ix=1,nnx
              p(ix,iy,iz) = 0.0
            end do
          end do
        end do

        return
c ---------------------------- process errors
  100   continue
        write(6,9000) path_ran, nvel
        call mpi_finalize(ierr)
        stop

 9992   continue
        write(6,6100) nvel,iz
        call mpi_finalize(ierr)
        stop
c ---------------------
 6001   format(' SR. GET_FIELDS: FILE READ FOR INITIALIZATION = ',a80)
 6005   format(' 6005, SR. GET_FIELDS: cannot find restart file = ',a80)
 6100   format(' SR. GET_FIELDS: file read error on unit number = ',i2,/
     +      ,'               at iz = ',i4)
 9000   format(' 9000, SR. GET_FIELDS: cannot open file =',a80,/,
     +      ' to unit number = ',i2)
      end
      subroutine read_res

c -------------- read restart file including constant file
c                changed for iys:iye

        use pars
        use fields
        use con_data
        use con_stats
#if defined(SWAP)
        use module_byteswap
#endif
        include 'mpif.h'

        integer status(mpi_status_size), ierr
        integer(kind=mpi_offset_kind) :: offset, disp
        integer(kind=k8)              :: nsize, nsize2
        real, allocatable, dimension(:,:,:) :: temp
        allocate(temp(nvar+4,nnx,iys:iye))

c ---- open file

        call mpi_file_open(mpi_comm_world, path_res,
     +      mpi_mode_create+mpi_mode_rdwr,
     +      mpi_info_null, nvel, ierr)

c ---- set file view

        disp = 0
        call mpi_file_set_view(nvel,disp,mpi_real8,mpi_real8,
     +      'native',mpi_info_null,ierr)

c ------------ read 3d fields

        nsize  = int(nvar+4,k8)*nnx*nny
        nsize2 = int(nvar+4,k8)*nnx*(iys-1)
        n_read = (nvar+4)*nnx*(iye+1-iys)

        do k=izs,ize
          offset = int((k-1),k8)*nsize + nsize2
          call mpi_file_read_at_all(nvel,offset,temp,n_read,
     +        mpi_real8,status,ierr)
          if (ierr /= 0) goto 9992
#if defined(SWAP)
          call byteswap(temp)
#endif
          do j=iys,iye
            do i=1,nnx
              u(i,j,k)    = temp(1,i,j)
              v(i,j,k)    = temp(2,i,j)
              w(i,j,k)    = temp(3,i,j)
              uf(i,j,k)   = temp(4,i,j)
              vf(i,j,k)   = temp(5,i,j)
              wf(i,j,k)   = temp(6,i,j)
              e(i,j,k)    = temp(nvar+3,i,j)
              gm(i,j,k,5) = temp(nvar+4,i,j)
            end do
          end do
          do is = 1,nscl
            do j = iys,iye
              do i = 1,nnx
                t(i,j,is,k) = temp(6+is,i,j)
              end do
            end do
          end do

        end do

c ---- close file

        call mpi_file_close(nvel, ierr)

        deallocate(temp)

c ------------ every mpi process reads constant file

        rewind(nvelc)
        read(nvelc,err=9993) c_c, c_s, case
        close(nvelc)

        if(l_root) write(6,4001) case
 4001   format(' 4001, SR. RESTART: case from restart = ',a3)

c ----- special restart conditions -------------------------------------

c -------- set case name to case input

        case   = case_inp
        if(l_root) write(6,4002) case_inp, utau, utausv
 4002   format(' 4002, SR. RESTART:',/,
     +      ' files will be saved with case name = ',a3,/,
     +      ' utau = ',e15.6,' utausv = ',e15.6)

c ------------------- if new vis model set match point for
c                     outer grid
        nmatch = 48
        utau = utausv

c -------- hand coded changes to restart if needed

        qstars = 0.000
        wtsfcs = 0.000
        zos    = 250.0e-06

c --------- flaw in t_zero ... fix in later runs

        if(iti < move_step) then
          t_zero = time
        end if

        if(l_root) write(6,4012) time, t_zero, dt, dt1, dt_new
        if(l_root) write(6,4013) qstar(1) , nmatch, case, zo

        call get_dz

        return
 9992   continue
        write(6,6100) nvel,iz
 6100   format(' SR. READ_RES: error reading file on unit number = ',i2,
     +      /,'               at iz = ',i4)
        call mpi_finalize(ierr)
        stop
c ---------------------
 9993   continue
        write(6,6200) nvelc
 6200   format(' SR. READ_RES:',/,
     +      '    error reading constant file on unit number = ',i2)
        call mpi_finalize(ierr)
        stop
c ---------------------
 4012   format(' SR. RESTART: ',/,
     +      ' T = ',e15.6,5x,' T_zero = ',e15.6,/,
     +      ' dt = ',e15.6,' dt1 = ',e15.6,' dt_new = ',e15.6)
 4013   format(' after restart qstar = ',e15.6,' nmatch = ',i5,
     +      ' case = ',a3,' zo = ',e15.6)
      end
      subroutine save_viz(it)

c --------------- save multiple (x-y), (x-z), (y-z), planes of data .
c                 modify recl in all open statements for more or less
c                 variables.
c                 Constant - x, implies yz planes
c                 Constant - y, implies xz planes
c                 Constant - z, implies xy planes

c ------------- routine uses send/recv to get information in y-z planes

        use pars
        use fields
        use con_data
        use con_stats
        use fftwk
#if defined(SWAP)
        use module_byteswap
#endif
        include 'mpif.h'

c ------ for xz and yz planes output z grid also

        parameter(nvar_o = 6, nvar_o1 = nvar_o + 1)

        integer ix_pick(maxnx),  iy_pick(maxny),  iz_pick(maxnz),
     +      ix_order(maxnx), iy_order(maxny), iz_order(maxnz)

        integer istatus(mpi_status_size), ierr
        integer(kind=mpi_offset_kind) :: offset, disp
        integer(kind=k8)              :: nsize, nsize2

        real(kind=4), dimension(nvar_o1,nny,izs:ize) :: temp_x
        real(kind=4), dimension(nvar_o1,nnx,izs:ize) :: temp_y
        real(kind=4), dimension(nvar_o,nnx,iys:iye)  :: temp_z
        real, dimension(nvar_o1,iys:iye,izs:ize)     :: buf_send
        real(kind=8), dimension(nnx,0:nnz)           :: z_xz
        real(kind=8), dimension(nny,0:nnz)           :: z_yz

c ------------- don't touch

        data iviz_x,  iviz_y,  iviz_z  /0, 0, 0/
        data ionce_x, ionce_y, ionce_z, istuff /0, 0, 0, 0/
        data ix_pick, iy_pick, iz_pick /maxnx*0, maxny*0, maxnz*0/
        data ix_order, iy_order, iz_order /maxnx*0, maxny*0, maxnz*0/
        save iviz_x,  iviz_y,  iviz_z,
     +      ix_pick, iy_pick, iz_pick,
     +      ix_order, iy_order, iz_order,
     +      ionce_x, ionce_y, ionce_z, istuff,
     +      npln_x, npln_y, npln_z

c ----------- turn on z levels to save. Customize for your own use.
c             Set iz_pick(iz) = iz, ix_pick(ix) = ix, iy_pick(iy) = iy
c             Data is round-robin alternated in the data file for more than
c             1 plane for any particular view.

        iz_pick(3)   = 3
        iz_pick(6)   = 6
        iz_pick(10)  = 10
        iz_pick(20)  = 20
        iz_pick(40)  = 40

c -------------- pick an x-z plane of data (can add more)

        iy_pick(1)         = 1
        iy_pick(nny/4)     = nny/4
        iy_pick(nny/2)     = nny/2
        iy_pick((3*nny)/4) = (3*nny)/4

c -------------- pick a y-z plane of data (can add more)

        ix_pick(nnx/4) = nnx/4
        ix_pick(nnx/2) = nnx/2
        ix_pick(nnx)   = nnx

c ------ find total number of z's turned on and open file once

        if(ionce_z == 0) then
          npln_z = 0
          do k=1,nnz
            if(iz_pick(k) == k) then
              npln_z = npln_z + 1
              iz_order(k) = npln_z
            end if
          end do
          ionce_z = 1
          iviz_z =  -npln_z
          if(npln_z /= 0) then
            call mpi_file_open(mpi_comm_world, path_viz_xy,
     +          mpi_mode_create+mpi_mode_rdwr,
     +          mpi_info_null, nviz_z, ierr)
            disp = 0
            call mpi_file_set_view(nviz_z,disp,mpi_real4,mpi_real4,
     +          'native',mpi_info_null,ierr)
          end if
        end if

c ------ find total number of y's turned on and open file once

        if (ionce_y == 0) then
          npln_y = 0
          do j=1,nny
            if(iy_pick(j) == j) then
              npln_y = npln_y + 1
              iy_order(j) = npln_y
            end if
          end do
          ionce_y = 1
          iviz_y  = -npln_y
          if(npln_y /= 0) then
            call mpi_file_open(mpi_comm_world, path_viz_xz,
     +          mpi_mode_create+mpi_mode_rdwr,
     +          mpi_info_null, nviz_y, ierr)
            disp = 0
            call mpi_file_set_view(nviz_y,disp,mpi_real4,mpi_real4,
     +          'native',mpi_info_null,ierr)
          end if
        end if

c ------ find total number of x's turned on and open file once

        if(ionce_x == 0) then
          npln_x = 0
          do i=1,nnx
            if(ix_pick(i) == i) then
              npln_x = npln_x + 1
              ix_order(i) = npln_x
            end if
          end do
          ionce_x = 1
          iviz_x  = -npln_x
          if(npln_x /= 0) then
            call mpi_file_open(mpi_comm_world, path_viz_yz,
     +          mpi_mode_create+mpi_mode_rdwr,
     +          mpi_info_null, nviz_x, ierr)
            disp = 0
            call mpi_file_set_view(nviz_x,disp,mpi_real4,mpi_real4,
     +          'native',mpi_info_null,ierr)
          end if
        end if

        if(istuff == 0 .and. l_root) then
          open(nviz_s,file=path_stuf)
          istuff = 1
        end if

c --------- write data, subtract t_ref to increase
c           resolution on 32 bit machines

c ---------- xy planes of data

        iviz_z  = iviz_z + npln_z
        nsize   = int(nvar_o,k8)*nnx*nny
        nsize2  = int(nvar_o,k8)*nnx*(iys-1)
        n_write = nvar_o*nnx*(iye+1-iys)
        do k=izs,ize
          if(iz_pick(k) == k) then
            km1 = k - 1
            do j=iys,iye
              do i=1,nnx
                temp_z(1,i,j) = u(i,j,k)
                temp_z(2,i,j) = v(i,j,k)
                temp_z(3,i,j) = w(i,j,k)
                temp_z(4,i,j) = (t(i,j,1,k) - t_ref)

c ---------- get the pressure field assumes flux form for advection

                temp_z(5,i,j) = p(i,j,k) - c23*e(i,j,k) - pxym(k)

c ---------- get the flux velocity minus the grid speed

                temp_z(6,i,j) = wf(i,j,k) - gm(i,j,k,5)
              end do
            end do
#if defined(SWAP)
            call byteswap(temp_z)
#endif
            offset = int((iviz_z + iz_order(k) - 1),k8)*nsize + nsize2
            call mpi_file_write_at(nviz_z,offset,temp_z,n_write,
     +          mpi_real4,istatus,ierr)
            if (ierr /= 0) go to 9991
          end if
        end do

c ---------- xz planes of data

        iviz_y = iviz_y + npln_y
        nsize  = int(nvar_o1,k8)*nnx*nnz
        nsize2 = int(nvar_o1,k8)*nnx*(izs-1)
        nwrite = (nvar_o1)*nnx*(ize+1-izs)
        do j=iys,iye
          if(iy_pick(j) == j) then
            do k=izs,ize
              km1 = k - 1
              do i=1,nnx
                temp_y(1,i,k) = u(i,j,k)
                temp_y(2,i,k) = v(i,j,k)
                temp_y(3,i,k) = w(i,j,k)
                temp_y(4,i,k) = (t(i,j,1,k) - t_ref)

c ---------- get the pressure field

                temp_y(5,i,k) =  p(i,j,k) - pxym(k) - c23*e(i,j,k)

c ---------- get the flux velocity minus the grid speed

                temp_y(6,i,k) = wf(i,j,k) - gm(i,j,k,5)

c ---------- the shape of the gridlines

                temp_y(7,i,k) =  0.5*(zgrid_w(i,j,km1,1) +
     +              zgrid_w(i,j,k,1))
              end do
            end do
#if defined(SWAP)
            call byteswap(temp_y)
#endif
            offset = int((iviz_y + iy_order(j) - 1),k8)*nsize + nsize2
            call mpi_file_write_at(nviz_y,offset,temp_y,nwrite,
     +          mpi_real4,istatus,ierr)
            if (ierr /= 0) goto 9992
          end if
        end do

c ---------- yz planes that cut across all processors
c            just have root node on that slab write data

        iviz_x  = iviz_x + npln_x
        n_write = nvar_o1*nny*(ize+1-izs)
        nsize   = int(nvar_o1,k8)*nny*nnz
        nsize2  = int(nvar_o1,k8)*nny*(izs-1)
        n_send  = nvar_o1*(ize+1-izs)*(iye+1-iys)
        do i=1,nnx
          if(ix_pick(i) == i) then

c ----------- build send buffer

            do k=izs,ize
              km1 = k - 1
              do j=iys,iye
                buf_send(1,j,k) = u(i,j,k)
                buf_send(2,j,k) = v(i,j,k)
                buf_send(3,j,k) = w(i,j,k)
                buf_send(4,j,k) = (t(i,j,1,k) - t_ref)

c ---------- get the pressure field

                buf_send(5,j,k) = p(i,j,k) - pxym(k) - c23*e(i,j,k)

c ---------- get the flux velocity minus the grid speed

                buf_send(6,j,k) = wf(i,j,k) - gm(i,j,k,5)

c ---------- shape of the gridlines

                buf_send(7,j,k) = 0.5*(zgrid_w(i,j,km1,1) +
     +              zgrid_w(i,j,k,1))
              end do
            end do
            if(myid /= iss) then
              call mpi_send(buf_send(1,iys,izs),n_send,
     +            mpi_real8,iss,1,
     +            mpi_comm_world,ierr)
            else
              do k=izs,ize
                do j=iys,iye
                  do ii=1,nvar_o1
                    temp_x(ii,j,k) = buf_send(ii,j,k)
                  end do
                end do
              end do
              do l=iss+1,ise
                call recv_yz_var(temp_x,nvar_o1,nny,
     +              iy_s(l),iy_e(l),izs,ize,l)
              end do
#if defined(SWAP)
              call byteswap(temp_x)
#endif
              offset = int((iviz_x + ix_order(i) - 1),k8)*nsize + nsize2
              call mpi_file_write_at(nviz_x,offset,temp_x,n_write,
     +            mpi_real4,istatus,ierr)
              if (ierr /= 0) goto 9993
            end if
          end if
        end do

c ------------- ascii file with facts in it that goes
c               with visualization

        if(l_root) then
          write(nviz_s,5000) time, amonin, zi, utau
 5000     format(4e20.8)
        end if

c ---- last time step or hit mtape close the files

        if (it == itmax .or. mtape) then
          if(npln_z /= 0) then
            call mpi_file_close(nviz_z, ierr)
            ionce_z = 0
          end if
          if(npln_y /= 0) then
            call mpi_file_close(nviz_y, ierr)
            ionce_y = 0
          end if
          if(npln_x /= 0) then
            call mpi_file_close(nviz_x, ierr)
            ionce_x = 0
          end if
          if(l_root) then
            close(nviz_s)
            istuff = 0
          end if
        end if

        return
c --------------------------  errors in writing viz file
 9991   continue
        write(6,6000) nviz_z, iz
 6000   format(' SR. SAVE_VIS:',/,
     +      '    trouble cannot write xy viz file on unit = ',i2,/,
     +      '             at iz = ',i4)
        call mpi_finalize(ierr)
        stop
c --------------------------  errors in writing viz file
 9992   continue
        write(6,6100) nviz_y, iz, iviz_y
 6100   format(' SR. SAVE_VIS:',/,
     +      '    trouble cannot write xz viz file on unit = ',i2,/,
     +      '             at iz = ',i4,/,
     +      '            iviz_y = ',i8)
c --------------------------  errors in writing viz file
 9993   continue
        write(6,6200) nviz_x, iz, iviz_x
 6200   format(' SR. SAVE_VIS:',/,
     +      '    trouble cannot write yz viz file on unit = ',i2,/,
     +      '             at iz = ',i4,/,
     +      '            iviz_x = ',i8)
        call mpi_finalize(ierr)
        stop
      end
      subroutine save_bndy(it)

c --------------- save boundary shape and boundary speed
c                 as function of time
c                 at same time step as viz files

        use pars
        use fields
        use con_data
        use con_stats
        use fftwk
#if defined(SWAP)
        use module_byteswap
#endif
        include 'mpif.h'

        parameter(nvar_o = 2)

        integer istatus(mpi_status_size), ierr
        integer(kind=mpi_offset_kind) :: offset, disp
        integer(kind=k8)              :: nsize, nsize2

        real(kind=4), dimension(2,nnx,iys:iye) :: temp_z

c ------------- don't touch

        data ionce_z, iviz_z /0, 0/
        save ionce_z, iviz_z

c ------ open file once

c     write(nprt,4341) nbndy, path_bndy
 4341   format(' 4341 nbndy = ',i5,' path_bndy = ',a80)

        if(ionce_z == 0) then
          ionce_z = 1
          call mpi_file_open(mpi_comm_world, path_bndy,
     +        mpi_mode_create+mpi_mode_rdwr,
     +        mpi_info_null, nbndy, ierr)
          disp = 0
          call mpi_file_set_view(nbndy,disp,mpi_real4,mpi_real4,
     +        'native',mpi_info_null,ierr)
        end if

c ---------- xy plane of data

        iviz_z  = iviz_z + 1
        nsize   = int(nvar_o,k8)*nnx*nny
        nsize2  = int(nvar_o,k8)*nnx*(iys-1)
        n_write = nvar_o*nnx*(iye+1-iys)

        if(iss /= 0) go to 999

        k = 0

        do j=iys,iye
          do i=1,nnx
            temp_z(1,i,j) = bndy_t(i,j,1)
            temp_z(2,i,j) = gm(i,j,k,5)
          end do
        end do
#if defined(SWAP)
        call byteswap(temp_z)
#endif
        offset = int((iviz_z - 1),k8)*nsize + nsize2
        call mpi_file_write_at(nbndy,offset,temp_z,n_write,
     +      mpi_real4,istatus,ierr)
        if (ierr /= 0) go to 9991

  999   continue

c ---- last time step or hit mtape close the files

        if(it == itmax .or. mtape) then
          call mpi_file_close(nbndy, ierr)
          ionce_z = 0
          iviz_z  = 0
        end if

        return
c --------------------------  errors in writing bndy file
 9991   continue
        write(6,6000) nbndy
 6000   format(' SR. SAVE_BNDY:',/,
     +      '    trouble cannot write boundary file on unit = ',i2)
        call mpi_finalize(ierr)
        stop
      end
      subroutine save_v(it)

c --------------- save 3d fields

        use pars
        use fields
#if defined(SWAP)
        use module_byteswap
#endif
        include 'mpif.h'
        logical there

        integer status(mpi_status_size), ierr
        integer(kind=mpi_offset_kind) :: offset, disp
        integer(kind=k8)                 nsize, nsize2

        real, allocatable, dimension(:,:,:) :: temp
        allocate(temp(nvar+4,nnx,iys:iye))

c ---- open file

        call mpi_file_open(mpi_comm_world, path_sav_v,
     +      mpi_mode_create+mpi_mode_rdwr,
     +      mpi_info_null, nvel, ierr)

c ---- set file view

        disp = 0
        call mpi_file_set_view(nvel,disp,mpi_real8,mpi_real8,
     +      'native',mpi_info_null,ierr)

c ---- write data

        nsize   = int(nvar+4,k8)*nnx*nny
        nsize2  = int(nvar+4,k8)*nnx*(iys-1)
        n_write = (nvar+4)*nnx*(iye+1-iys)

        do k=izs,ize
          do j = iys,iye
            do i = 1,nnx
              temp(1,i,j)      = u(i,j,k)
              temp(2,i,j)      = v(i,j,k)
              temp(3,i,j)      = w(i,j,k)
              temp(4,i,j)      = uf(i,j,k)
              temp(5,i,j)      = vf(i,j,k)
              temp(6,i,j)      = wf(i,j,k)
              temp(nvar+3,i,j) = e(i,j,k)
              temp(nvar+4,i,j) = gm(i,j,k,5)
            end do
          end do
          do is = 1,nscl
            do j = iys,iye
              do i = 1,nnx
                temp(6+is,i,j) = t(i,j,is,k)
              end do
            end do
          end do
#if defined(SWAP)
          call byteswap(temp)
#endif

          offset = int((k-1),k8)*nsize + nsize2
          call mpi_file_write_at(nvel,offset,temp,n_write,
     +        mpi_real8,status,ierr)
          if (ierr /= 0) goto 9991

        end do

c ---- close file

        call mpi_file_close(nvel, ierr)

c ---- check file

        if (l_root) then
          inquire(file=path_sav_v,exist=there)
          if(.not.there) then
            write(6,8000) nvel,myid
            call mpi_finalize(ierr)
            stop
          end if
          write(6,7000) it,path_sav_v
        end if

        deallocate(temp)

        return
c --------------------------  errors in writing restart file
 9991   continue
        write(6,6000) nvel, iz
 6000   format(' SR. SAVE_V:',/,
     +      '    trouble cannot write restart file on unit = ',i2,/,
     +      '             at iz = ',i4)
        call mpi_finalize(ierr)
        stop
c --------------------
 7000   format(' **** DATA SET AT IT = ',I6,/,
     +      '      VELOCITY DATA IS WRITTEN IN FILE  ',a80)
 8000   format(' in SAVE_V: trouble writing file ',i5,'  myid = ',i5,
     +      ' at iz = ',i5)
      end
      subroutine save_c(it)

c --------------- root process writes constant file
c                 sequential fortan binary

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats

        logical there
        character options*8, passwd*1

c ---- open file

        open(nvelc,err=9992,file=path_sav_c,form='unformatted',
     +      status='unknown')
        write(nvelc,err=9992) c_c, c_s, case
        close(nvelc)

        inquire(file=path_sav_c,exist=there)
        if(.not.there) then
          write(6,8001) path_sav_c
          call mpi_finalize(ierr)
          stop
        end if
c -----------------------------  output ok message
        write(6,7001) path_sav_c

        return
c --------------------------  errors in writing constant file
 9992   continue
        write(6,6100) nvelc
 6100   format(' SR. SAVE_V:',/,
     +      '    trouble cannot open/write constant file on unit = ',i2)
        call mpi_finalize(ierr)
        stop
c ---------------------
 7001   format('      CONSTANT DATA IS WRITTEN IN FILE  ',a80)
 8001   format(' SR. SAVE_C: Trouble constant file not in path =',a80)
      end
      subroutine save_p

c -------------- save pressure file

        use pars
        use fields
#if defined(SWAP)
        use module_byteswap
#endif
        include 'mpif.h'
        logical there

        integer status(mpi_status_size), ierr
        integer(kind=mpi_offset_kind) :: offset, disp
        integer(kind=k8)              :: nsize, nsize2

        real temp(nnx,iys:iye)

c ---- open file

        call mpi_file_open(mpi_comm_world, path_sav_p,
     +      mpi_mode_create+mpi_mode_rdwr,
     +      mpi_info_null, npre, ierr)

c ---- set file view

        disp = 0
        call mpi_file_set_view(npre,disp,mpi_real8,mpi_real8,
     +      'native',mpi_info_null,ierr)

c ---- write data

        nsize   = int(nnx,k8)*nny
        nsize2  = int(nnx,k8)*(iys -1)
        n_write = nnx*(iye+1-iys)
        do k=izs,ize
          do j=iys,iye
            do i=1,nnx
              temp(i,j) = p(i,j,k)
            end do
          end do
#if defined(SWAP)
          call byteswap(temp)
#endif
          offset = int((k-1),k8)*nsize + nsize2
          call mpi_file_write_at(npre,offset,temp,n_write,
     +        mpi_real8,status,ierr)
        end do

c ---- close file

        call mpi_file_close(npre, ierr)

c ---- check file

        if (l_root) then
          inquire(file=path_sav_p,exist=there)
          if(.not.there) then
            write(6,8000) path_sav_p
            call mpi_finalize(ierr)
            stop
          end if
          write(6,7000) path_sav_p
        end if

        return
c -------------------- process write errors
 9991   continue
        write(6,6000) npre, iz
 6000   format(' SR. SAVE_P:',/,
     +      '    trouble cannot write pressure file on unit = ',i2,/,
     +      '             at iz = ',i4)
        call mpi_finalize(ierr)
        stop
c -----------------------
 7000   format('      PRESSURE DATA IS WRITTEN IN FILE  ',a80)
 8000   format(' SR. SAVE_P: Trouble pressure file not in path =',a80)
      end
      subroutine get_units

        use pars

c -------------- unit numbers for files

        nvel   = 20
        npre   = 30
        nhis1  = 40
        nvelc  = 50
        nhisp  = 60
        nviz_z = 80
        nviz_y = 82
        nviz_x = 84
        nviz_s = 90
        nbndy  = 91
        ngrd   = 25

c ------------- unit number for standard print out
c               for each mpi task

        nprt = 1

c ------------- open unit for standard printout

        return
      end
      subroutine get_output_filenames

c ----------- build file names for velocity, pressure, and constants

        use pars
        include 'mpif.h'
        character cgrid*10, num*3

c --------------- build character strings for file name

        cgrid = '.le.'
        if (itn > 999) then
          num = "OVR"
        else
          write(num,'(i3.3)') itn
        end if
        call blnk(path_sav_v)
        call blnk(path_sav_p)
        call blnk(path_sav_c)
        call blnk(path_sav_part)
        path_sav_v = trim(path_sav)//'/u'//
     +      cgrid(1:4)//case(1:3)//num(1:3)
        path_sav_p = trim(path_sav)//'/p'//
     +      cgrid(1:4)//case(1:3)//num(1:3)
        path_sav_c = trim(path_sav)//'/u'//
     +      cgrid(1:4)//case(1:3)//num(1:3)//'.con'
        path_sav_part = trim(path_sav)//'/part'//
     +      cgrid(1:4)//case(1:3)//num(1:3)
c     path_sav_g = trim(path_sav)//'/'//case(1:3)//'.grid'

        return
      end
      subroutine open_his(istep)

c ------------------- open history files by root
c                     isize determined in sr. fill_cs

        use pars
        include 'mpif.h'
        character cgrid*4, iblks*16
        logical there

c --------------- build character strings for ascii history file name

        cgrid = '.le.'
        call blnk(iblks)
        write(iblks(1:7),'(i7.7)') istep
        iblks(8:8) = '_'
        write(iblks(9:15),'(i7.7)') (istep + itape)
        iblnk = index(path_his,' ')
        call blnk(path_sav_h)
        path_sav_h = path_his(1:iblnk-1)//'/his'//
     +      cgrid(1:4)//case(1:3)//'.'//iblks(1:15)

c --------------- build character strings for ieee profile history file
c                 set record counter for direct access file = 0

        krec = 0
        cgrid = '.le.'
        call blnk(iblks)
        write(iblks(1:7),'(i7.7)') istep
        iblks(8:8) = '_'
        write(iblks(9:15),'(i7.7)') (istep + itape)
        iblnk = index(path_his,' ')
        call blnk(path_sav_hp)
        path_sav_hp = path_his(1:iblnk-1)//'/his'//
     +      cgrid(1:4)//case(1:3)//'.'//iblks(1:15)//'.ieee'

c ----------------- save data in directory

        if(l_root) then

          close(nhis1)
          open(nhis1,err=3000,file=path_sav_h,form='formatted')

          close(nhisp)
          open(nhisp,err=4000,file=path_sav_hp,
     +        form='unformatted',access='direct',recl=isize*j_recl,
     +        status='unknown')
        end if

        return
c ------------------- process errors
 3000   continue
        write(6,6301) nhis1, path_sav_h
 6301   format(' 6301, SR. OPEN_HIS:',/,
     +      '    cannot open history1 file on unit = ',i2,/,
     +      '    path = ',a80)
        stop
c-------------------
 4000   continue
        write(6,6302) nhisp, path_sav_hp
 6302   format(' 6302, SR. OPEN_HIS:',/,
     +      '    cannot open history profile file on unit = ',i2,/,
     +      '    path = ',a80)
        stop
      end
      subroutine viz_output_filename(istep)

c ------------------- set visualization files,
c                     leaves files in scratch directory

        use pars
        include 'mpif.h'
        character iblks*16

c --------------- build character strings for file names
c                 with time step

        call blnk(iblks)
        iblks(1:1) = '.'
        write(iblks(2:8),'(i7.7)') istep
        iblks(9:9) = '_'
        write(iblks(10:16),'(i7.7)') (istep + itape)

        iloc = index(path_seed,' ')
        path_viz_xy = path_seed(1:iloc-1)
     +      //'/viz.'//case(1:3)//iblks(1:16)//'.xy.data'

        path_viz_xz = path_seed(1:iloc-1)
     +      //'/viz.'//case(1:3)//iblks(1:16)//'.xz.data'

        path_viz_yz = path_seed(1:iloc-1)
     +      //'/viz.'//case(1:3)//iblks(1:16)//'.yz.data'

        path_stuf = path_seed(1:iloc-1)
     +      //'/stuff.'//case(1:3)//iblks(1:16)//'.data'

        path_bndy = path_seed(1:iloc-1)
     +      //'/bndy.'//case(1:3)//iblks(1:16)//'.data'

        return
      end
