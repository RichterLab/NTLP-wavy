!!!!! ----- D. Richter, 3/11/14 ---
!!!!! ----- Adding particles to Peter's wavy code

      subroutine particle_init
      use particles
      use pars
      use con_data
      use tank_pars
      implicit none
      include 'mpif.h'
      integer :: values(8)
      integer :: idx,ix
      real :: xv,yv,zv,ran2,xp_init(3)
      real :: yp1,ypn,xpart(nnx),y2part(nnx),ddx

      !Create the seed for the random number generator:
      call date_and_time(VALUES=values)
      iseed = -(myid+values(8)+values(7)+values(6))


      !Initialize ngidx, the particle global index
      ngidx = 1

      !For the channel case, set the total number of particles:
      numpart = tnumpart/numprocs
      if (myid==0) then
         numpart = numpart + mod(tnumpart,numprocs)
      end if

      !Initialize the linked list of particles:
      nullify(part,first_particle)

      yp1 = 2.0e+30
      ypn = yp1
      !Set up xpart for interpolation of h:
      ddx = xl/float(nnx)
      do ix=1,nnx
         xpart(ix) = float(ix-1)*ddx
      end do
      ! call spline(xpart,z_tank_i,nnx,yp1,ypn,y2part)

      !Now initialize all particles with a random location on that processor
      do idx=1,numpart
         xv = ran2(iseed)*(xmax-xmin) + xmin
         yv = ran2(iseed)*(ymax-ymin) + ymin
         zv = ran2(iseed)*zl
         xp_init = (/xv,yv,zv/)
         call create_particle(xp_init,vp_init,Tp_init,radius_init,idx,
     +       z_tank_i(1:nnx),xpart(1:nnx),y2part(1:nnx),0.0)
      end do

      partTsrc = 0.0
      partTsrc_t = 0.0

      end subroutine particle_init
      subroutine particle_setup

      use particles
      use pars
      use tank_pars
      use con_data
      implicit none
      include 'mpif.h'

      integer :: blcts(3),types(3)
      integer :: ierr
      integer(kind=MPI_ADDRESS_KIND) :: extent,lb
      integer(kind=MPI_ADDRESS_KIND) :: extent2,lb2,displs(3)
      real :: ddx
      integer :: ix
      integer :: num_reals

      !First set up the neighbors for the interpolation stage:
      call assign_nbrs

      !Also assign the x,y,z max and mins to track particles leaving
      call set_bounds

      num_reintro = 0
      tot_reintro = 0
      his_reintro = 0

      !Recall that these are the bin centers
      !Setting the dimensional maximum time:
!      call set_binsdata(bins,histbins+2,numhistvars,
!     +                  log10(0.001),log10(100.0),1)  !Residence time
      !Setting the residence time based on NONDIMENSIONAL time
      call set_binsdata(bins,histbins+2,numhistvars,
     +                  log10(0.001),log10(100.0*zl/ustar_tank),1)  !Residence time
      call set_binsdata(bins,histbins+2,numhistvars,0.0,1.0,2)  !wave-relative origin
      call set_binsdata(bins,histbins+2,numhistvars,0.0,1.0,3)  !wave-relative death
      call set_binsdata(bins,histbins+2,numhistvars,0.0,zl,4)  !max height
      call set_binsdata(bins,histbins+2,numhistvars,-1.0,20.0,5)  !ufinal - u0

      !Initialize the linked list of particles:
      nullify(part,first_particle)

      !Set up MPI datatypes for sending particle information
      !MUST UPDATE IF THINGS ARE ADDED/REMOVED FROM PARTICLE STRUCTURE

      num_reals = 3*10 + 9

      blcts(1:3) = (/2,num_reals,2/)
      displs(1) = 0
      types(1) = mpi_integer
      call mpi_type_get_extent(mpi_integer,lb,extent,ierr)

      !Displace 2*size of mpi_integer (2 integer: pidx,procidx)
      displs(2) = extent*2
      types(2) = mpi_real8
      call mpi_type_get_extent(mpi_real8,lb,extent,ierr)
      !Displace num_reals*size of mpi_real8 (6 3-vectors plus Tp,Tprhs,Tf)
      displs(3) = displs(2) + extent*num_reals
      types(3) = mpi_integer8

      !Now define the type:
      call mpi_type_create_struct(3,blcts,
     +            displs,types,particletype,ierr)


       call mpi_type_get_true_extent(particletype,lb2,extent2,ierr)
       call mpi_type_get_extent(particletype,lb2,extent,ierr)
       if (extent .NE. sizeof(part) ) then
          if (myid==0) then
          write(*,*) 'WARNING: extent of particletype not equal
     +                  to sizeof(part):'
          write(*,*) 'sizeof(part) = ', sizeof(part)
!          write(*,*) 'sizeof(part%pidx) = ', sizeof(part%pidx)
          write(*,*) 'mpi_type_get_true_extent(particletype) = ',extent2
          write(*,*) 'mpi_type_get_extent(particletype) = ',extent
          end if
       end if

      !Need to compute any padding which may exist in particle struct:
      pad_diff = extent-extent2
      if (myid==0) then
      write(*,*) 'mpi_get_extent = ',extent
      write(*,*) 'mpi_get_true_extent = ',extent2
      write(*,*) 'sizeof(part) = ',sizeof(part)
      write(*,*) 'DIFF = ',pad_diff
      end if
      if (pad_diff .LT. 0) then
        write(*,*) 'WARNING: mpi_get_extent - mpi_get_true_extent LT 0!'
        call mpi_finalize(ierr)
        stop
      end if


      call mpi_type_commit(particletype,ierr)

      end subroutine particle_setup
      subroutine assign_nbrs
        use pars
        use particles
        include 'mpif.h'
      !Figure out which processors lie to all sides:
      !For the case where particles lie in columns not
      !aligning with the velocity, there will be no MPI_PROC_NULL since
      !x and y are BOTH periodic

      !On right boundary:
      if ( mod(myid+1,ncpu_s) == 0 ) then
         !On the top:
         if ( myid .GE. ncpu_s*(ncpu_z-1) ) then
            rproc = myid-ncpu_s+1
            trproc = 0
            tproc = ncpu_s-1
            tlproc = ncpu_s-2
            lproc = myid-1
            blproc = myid-ncpu_s-1
            bproc = myid-ncpu_s
            brproc = myid-ncpu_s - ncpu_s+1
         !On the bottom:
         elseif ( myid .LT. ncpu_s ) then
            rproc = myid-ncpu_s+1
            trproc = myid+1
            tproc = myid+ncpu_s
            tlproc = myid+ncpu_s-1
            lproc = myid-1
            blproc = myid+ncpu_s*(ncpu_z-1)-1
            bproc = myid+ncpu_s*(ncpu_z-1)
            brproc = ncpu_s*(ncpu_z-1)
         !In the middle of right side:
         else
            rproc = myid-ncpu_s+1
            trproc = myid+1
            tproc = myid+ncpu_s
            tlproc = myid+ncpu_s-1
            lproc = myid-1
            blproc = myid-ncpu_s-1
            bproc = myid-ncpu_s
            brproc = myid-ncpu_s - ncpu_s+1
         end if

      !On the left boundary:
      elseif ( mod(myid,ncpu_s) == 0) then
         !On the top:
         if ( myid .GE. ncpu_s*(ncpu_z-1) ) then
            rproc = myid+1
            trproc = 1
            tproc = 0
            tlproc = ncpu_s-1
            lproc = myid+ncpu_s-1
            blproc = myid-1
            bproc = myid-ncpu_s
            brproc = myid-ncpu_s+1
         !On the bottom:
         elseif ( myid .LT. ncpu_s ) then
            rproc = myid+1
            trproc = myid+ncpu_s+1
            tproc = myid+ncpu_s
            tlproc = myid+ncpu_s+ncpu_s-1
            lproc = myid+ncpu_s-1
            blproc = numprocs-1
            bproc = ncpu_s*(ncpu_z-1)
            brproc = ncpu_s*(ncpu_z-1)+1
         !In the middle of left side:
         else
            rproc = myid+1
            trproc = myid+ncpu_s+1
            tproc = myid+ncpu_s
            tlproc = myid+ncpu_s + ncpu_s-1
            lproc = myid+ncpu_s-1
            blproc = myid-1
            bproc = myid-ncpu_s
            brproc = myid-ncpu_s+1
         end if
      !On the top boundary
      elseif ( myid .GE. ncpu_s*(ncpu_z-1) ) then
         !Only check if in the middle:
         if ( .NOT. ( mod(myid,ncpu_s) == 0) ) then
            if ( .NOT. (mod(myid+1,ncpu_s) == 0) ) then
               rproc = myid+1
               trproc = myid-(ncpu_s*(ncpu_z-1))+1
               tproc = myid-(ncpu_s*(ncpu_z-1))
               tlproc = myid-(ncpu_s*(ncpu_z-1))-1
               lproc = myid-1
               blproc = myid-ncpu_s-1
               bproc = myid-ncpu_s
               brproc = myid-ncpu_s+1
            end if
         end if
      !On the bottom boundary
      elseif ( myid .LT. ncpu_s) then
         if ( .NOT. ( mod(myid,ncpu_s) == 0) ) then
            if ( .NOT. (mod(myid+1,ncpu_s) == 0) ) then
               rproc = myid+1
               trproc = myid+ncpu_s+1
               tproc = myid+ncpu_s
               tlproc = myid+ncpu_s-1
               lproc = myid-1
               blproc = myid+ncpu_s*(ncpu_z-1)-1
               bproc = myid+ncpu_s*(ncpu_z-1)
               brproc = myid+ncpu_s*(ncpu_z-1)+1
            end if
         end if
      !Everywhere else:
      else
         rproc = myid+1
         trproc = myid+ncpu_s+1
         tproc = myid+ncpu_s
         tlproc = myid+ncpu_s-1
         lproc = myid-1
         blproc = myid-ncpu_s-1
         bproc = myid-ncpu_s
         brproc = myid-ncpu_s+1
      end if

      return
      end subroutine assign_nbrs
      subroutine set_bounds
      use pars
      use particles
      use con_data
      use con_stats
      implicit none
      include 'mpif.h'

      !Assign the bounds of each processor
      !(this probably already exists in code, but this is for my sake)
      ymin = dy*(iys-1)
      ymax = dy*(iye)
      zmin = z(izs-1)
      zmax = z(ize)
      xmin = dx*(mxs-1)
      xmax = dx*(mxe)

      end subroutine set_bounds
      subroutine create_particle(xp,vp,Tp,rad,idx,ztank,xles,y2,t_born)
      use particles
      use pars
      implicit none

      real :: xp(3),vp(3),Tp,hpart,rad,pi
      real :: ztank(nnx),xles(nnx),y2(nnx),t_born
      integer :: idx,ix

      if (.NOT. associated(first_particle)) then
         allocate(first_particle)
         part => first_particle
         nullify(part%next,part%prev)
      else
         !Add to beginning of list since it's more convenient
         part => first_particle
         allocate(part%prev)
         first_particle => part%prev
         part%prev%next => part
         part => first_particle
         nullify(part%prev)
      end if

      pi   = 4.0*atan(1.0)

      part%xip(1:3) = xp(1:3)
      ! call splint(xles,ztank,y2,nnx,part%xip(1),hpart)
      call linInterp(xles,ztank,nnx,part%xip(1),hpart)
      part%xp(1:2) = part%xip(1:2)
      part%xp(3) = part%xip(3) +
     +             hpart*(1.0-part%xip(3)/zl)**a_decay

      part%xip_i(1:3) = part%xip(1:3)
      part%xp_i(1:3) = part%xp(1:3)

      part%vp(1:3) = vp(1:3)
      part%u0 = vp(1)
      part%Tp = Tp
      part%Tf = Tp
      part%uf(1:3) = 0.0
      part%xrhs(1:3) = 0.0
      part%vrhs(1:3) = 0.0
      part%Tprhs_s = 0.0
      part%Tprhs_L = 0.0
      part%radrhs = 0.0
      part%qinf = 0.0
      part%pidx = idx
      part%procidx = myid
      part%radius = rad
      part%m_s = Sal/1000.0*4.0/3.0*pi*rad**3*rhow
      part%xp_0(1:3) = part%xp(1:3)
      part%t0 = t_born
      part%zpmax = 0.0

      end subroutine create_particle
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
     +                   mpi_comm_world,ierr)

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
     +                   mpi_mode_create+mpi_mode_rdwr,
     +                   mpi_info_null,fh,ierr)

      zoffset = 0
      !Write tnumpart first:
      if (myid==0) then
      call mpi_file_write_at(fh,zoffset,tnumpart,1,mpi_integer,
     +                       istatus,ierr)
      write(*,*) 'wrote tnumpart = ',tnumpart
      end if

      zoffset = zoffset + 4

      !Now compute the offset (in bytes!):
      offset = zoffset
      do iproc = 0,myid-1
         offset = offset + pnum_vec(iproc+1)*(sizeof(tmp)-pad_diff)
      end do

      !Now everyone else write, ONLY if numpart > 0
      if (numpart .GT. 0) then
      call mpi_file_write_at(fh,offset,writebuf,numpart,
     +                       particletype,istatus,ierr)
!      write(*,*) 'wrote writebuf(1)%pidx = ',writebuf(1)%pidx
!      write(*,'a30,3e15.6') 'wrote writebuf(1)%vp =',writebuf(1)%vp(1:3)
!      write(*,'a30,3e15.6') 'wrote writebuf(1)%xp =',writebuf(1)%xp(1:3)
!      write(*,*) 'wrote writebuf(2)%pidx = ',writebuf(2)%pidx
!      write(*,'a30,3e15.6') 'wrote writebuf(2)%vp =',writebuf(2)%vp(1:3)
!      write(*,'a30,3e15.6') 'wrote writebuf(2)%xp =',writebuf(2)%xp(1:3)
      end if

      call mpi_file_close(fh,ierr)

      write(*,*) 'proc',myid,'wrote numpart = ',numpart

      if (myid==0) write(*,7000) path_sav_part
 7000 format(' PARTICLE DATA IS WRITTEN IN FILE  ',a80)

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
 7000 format(' READING PARTICLE DATA FROM  ',a80)


      call mpi_file_open(mpi_comm_world,path_part,
     +                   mpi_mode_rdonly,
     +                   mpi_info_null,fh,ierr)


      !Read in the total number of particles:
      offset = 0
      call mpi_file_read_at_all(fh,offset,tnumpart,1,
     +                      mpi_integer,istatus,ierr)
      if (myid==0) write(*,*) 'read tnumpart = ',tnumpart

      offset = 4
      allocate(readbuf(tnumpart))
      call mpi_file_read_at_all(fh,offset,readbuf,tnumpart,
     +                       particletype,istatus,ierr)

      do i = 1,tnumpart
        !Now - does it lie within this proc's bounds?
        if (readbuf(i)%xp(2) .GT. ymin .AND.
     +       readbuf(i)%xp(2) .LT. ymax .AND.
     +       readbuf(i)%xp(1) .GT. xmin .AND.
     +       readbuf(i)%xp(1) .LT. xmax) then
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
         if (part%pidx .gt. pidxmax) pidxmax = part%pidx
         part => part%next
      end do

      !Set ngidx (the index for creating new particles) to 1+pidmax:
      ngidx = pidxmax + 1

      numpart = myp

      call mpi_allreduce(myp,totalp,1,mpi_integer,mpi_sum,
     +                   mpi_comm_world,ierr)

      write(*,*) 'proc',myid,'read in numpart:',myp
      if (myid==0) write(*,*) 'total number of particles read:',totalp

      end subroutine read_part_res
      subroutine particle_update_rk3(it,istage)
      use pars
      use particles
      use con_data
      use con_stats
      use tank_pars
      implicit none
      include 'mpif.h'

      integer :: istage,ierr,it
      real :: denom,dtl,sigma,ttest
      integer :: ix,iy,iz
      real :: Rep,diff(3),diffnorm,corrfac,myRep_avg,xtmp(3),vtmp(3)
      real :: Nup,Tptmp,taup_i,rhop,Volp,radiustmp,Shp
      real :: TfC,einf,TpC,Lv,Eff_C,Eff_S,estar
      real :: yp1,ypn,hpart,xpart(nnx),y2part(nnx),ddx

      real :: t_s,t_f,t_s1,t_f1



      !First fill extended velocity field for interpolation
      !t_s = mpi_wtime()
      call fill_ext
      !t_f = mpi_wtime()
      !call mpi_barrier(mpi_comm_world,ierr)
      !if (myid==5) write(*,*) 'time fill_ext:',t_f-t_s

      !Get spline fit of surface boundary to interpolate h(xp,yp):
      yp1 = 2.0e+30
      ypn = yp1
      ddx = xl/float(nnx)
      do ix=1,nnx
         xpart(ix) = float(ix-1)*ddx
      end do
      ! call spline(xpart,z_tank_i,nnx,yp1,ypn,y2part)

      !update xip since waves have moved since last particle
      part => first_particle
      do while (associated(part))
            !Update part%xip
            ! call splint(xpart,z_tank_i,y2part,nnx,part%xp(1),hpart)
            call linInterp(xpart,z_tank_i,nnx,part%xp(1),hpart)
            part%xip(1:3) = part%xp(1:3)
            if (iz_space == 5) then    !Algebraic stretching to smooth out waves
                call get_zeta_algeb(a_decay,zl,hpart,part%xip(1:3))
            end if
            part => part%next
      end do

      !Apply bottom boundary conditions here since mesh will have moved from previous:
      num_reintro = 0
      call particle_bcs_nonperiodic

      partcount_t = 0.0
      partcountU_t = 0.0
      if (istage .eq. 1) then
         partflux_t = 0.0
         partfluxup_t = 0.0
         partfluxdn_t = 0.0
      end if
      vpsum_t = 0.0
      upwp_t = 0.0
      vpsqrsum_t = 0.0
      Tpsum_t = 0.0
      Tpsqrsum_t = 0.0
      wpTpsum_t = 0.0
      myRep_avg = 0.0
      Tfsum_t = 0.0
      qfsum_t = 0.0
      radsum_t = 0.0
      rad2sum_t = 0.0
      qstarsum_t = 0.0

      !t_s = mpi_wtime()


      !Loop over the linked list of particles:
      part => first_particle
      do while (associated(part))

         !First, interpolate to get the fluid velocity part%uf(1:3):
         !For wavy code, must pass in the wave height at particle location
         !Get this from spline:
         ! call splint(xpart,z_tank_i,y2part,nnx,part%xp(1),hpart)
         call linInterp(xpart,z_tank_i,nnx,part%xp(1),hpart)

         call uf_interp(hpart)

         if (it .le. 1) then
            part%vp(1:3) = part%uf
            part%Tp = part%Tf
         end if

         diff(1:3) = part%vp - part%uf
         diffnorm = sqrt(diff(1)**2 + diff(2)**2 + diff(3)**2)
         Rep = 2.0*part%radius*diffnorm/nuf
         myRep_avg = myRep_avg + Rep
         corrfac = (1.0 + 0.15*Rep**(0.687))

         Volp = pi2*2.0/3.0*part%radius**3
         rhop = (part%m_s+Volp*rhow)/Volp
         taup_i = 18.0*rhoa*nuf/rhop/(2.0*part%radius)**2

         !Compute Nusselt number for particle:
         !Ranz-Marshall relation
         Nup = 2.0 + 0.6*Rep**(1.0/2.0)*Pra**(1.0/3.0)
         Shp = 2.0 + 0.6*Rep**(1.0/2.0)*Sc**(1.0/3.0)

         !Mass Transfer calculations
         TfC = part%Tf-273.15
         einf = 610.94*EXP(17.6257*TfC/(TfC+243.04))
         TpC = part%Tp-273.15
         Lv =(25.0 - 0.02274*26.0)*10.0**5

         Eff_C = 2.0*Mw*Gam/(Ru*rhow*part%radius*part%Tp)
         Eff_S = Ion*Os*part%m_s*Mw/Ms/(Volp*rhop-part%m_s)
         estar = einf*exp(Mw*Lv/Ru*(1.0/part%Tf-1.0/part%Tp)
     +                                           +Eff_C-Eff_S)
         part%qstar = Mw/Ru*estar/part%Tp/rhoa

         !Temp vectors for RK3
         xtmp(1:3) = part%xp(1:3) + dt*zetas(istage)*part%xrhs(1:3)
         vtmp(1:3) = part%vp(1:3) + dt*zetas(istage)*part%vrhs(1:3)
         Tptmp = part%Tp + dt*zetas(istage)*part%Tprhs_s
         Tptmp = Tptmp + dt*zetas(istage)*part%Tprhs_L
         radiustmp = part%radius + dt*zetas(istage)*part%radrhs


         !Assemble RHS
         part%xrhs(1:3) = part%vp(1:3)
         part%vrhs(1:3) = corrfac*taup_i*(part%uf(1:3)-part%vp(1:3))
     +                                  - part_grav(1:3)

         if (ievap .EQ. 1) then
            part%radrhs = Shp/9.0/Sc*rhop/rhow*part%radius*
     +                       taup_i*(part%qinf-part%qstar) !assumes qinf=rhov/rhoa rather than rhov/rhom
         else
            part%radrhs = 0.0
         end if

         part%Tprhs_s = -Nup/3.0/Pra*CpaCpp*
     +                     rhop/rhow*taup_i*(part%Tp-part%Tf)
         part%Tprhs_L = 3.0*Lv/Cpp/part%radius*part%radrhs


         !Update the particle properties
         if (istage==1) then
             part%xip_i(1:3) = part%xip(1:3)
             part%xp_i(1:3) = part%xp(1:3)
         end if

         part%xp(1:3) = xtmp(1:3) + dt*gama(istage)*part%xrhs(1:3)

         if (i_pinertia==1) then
           part%vp(1:3) = vtmp(1:3) + dt*gama(istage)*part%vrhs(1:3)
         else
           part%vp = part%uf - part_grav/taup_i
         end if
         part%Tp = Tptmp + dt*gama(istage)*part%Tprhs_s
         part%Tp = part%Tp + dt*gama(istage)*part%Tprhs_L
         part%radius = radiustmp + dt*gama(istage)*part%radrhs


         !Update part%xip as well:
         ! call splint(xpart,z_tank_i,y2part,nnx,part%xp(1),hpart)
         call linInterp(xpart,z_tank_i,nnx,part%xp(1),hpart)
         part%xip(1:3) = part%xp(1:3)
         if (iz_space == 5) then    !Algebraic stretching to smooth out waves
             call get_zeta_algeb(a_decay,zl,hpart,part%xip(1:3))
         end if

        !Keep track of max height particle reaches during lifetime
        if (part%xip(3) .gt. part%zpmax) part%zpmax = part%xip(3)


        part => part%next
      end do
      !t_f1 = mpi_wtime()
      !write(*,*) 'proc,loop time: ',myid,t_f1-t_s
      !call mpi_barrier(mpi_comm_world,ierr)
      !t_f = mpi_wtime()
      !if (myid==5) write(*,*) 'time loop:', t_f-t_s

      !Calculate particle fluxes before removing particles (to get lowest grid point correct)
      if (istage .eq. 3) call update_particle_flux

      !Enforce nonperiodic bcs (either elastic or destroying particles)
      !t_s = mpi_wtime()
      call particle_bcs_nonperiodic
      !call mpi_barrier(mpi_comm_world,ierr)
      !t_f = mpi_wtime()
      !if (myid==5) write(*,*) 'time bc_non:', t_f - t_s

      if (istage == 1) tot_reintro = 0
      call particle_reintro
      if (istage == 3) his_reintro = tot_reintro

      !Check to see if particles left processor
      !If they did, remove from one list and add to another
      !t_s = mpi_wtime()
      call particle_exchange
      !call mpi_barrier(mpi_comm_world,ierr)
      !t_f = mpi_wtime()
      !if (myid==5) write(*,*) 'time exchg:', t_f - t_s

      !Now enforce periodic bcs
      !just updates x,y locations if over xl,yl or under 0
      !t_s = mpi_wtime()
      call particle_bcs_periodic
      !call mpi_barrier(mpi_comm_world,ierr)
      !t_f = mpi_wtime()
      !if (myid==5) write(*,*) 'time bc_per:', t_f - t_s

      !Now that particles are in their updated position,
      !compute their contribution to the momentum coupling:
      !t_s = mpi_wtime()
      call update_particle_stats
      !call mpi_barrier(mpi_comm_world,ierr)
      !t_f = mpi_wtime()
      !if (myid==5) write(*,*) 'time cpl: ', t_f - t_s

      !Finally, now that coupling and statistics arrays are filled,
      !Transpose them back to align with the velocities:
      call ztox_trans(partsrc_t(0:nnz+1,iys:iye,mxs:mxe,1),
     +                partsrc(1:nnx,iys:iye,izs-1:ize+1,1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(partsrc_t(0:nnz+1,iys:iye,mxs:mxe,2),
     +                partsrc(1:nnx,iys:iye,izs-1:ize+1,2),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(partsrc_t(0:nnz+1,iys:iye,mxs:mxe,3),
     +                partsrc(1:nnx,iys:iye,izs-1:ize+1,3),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(partTsrc_t(0:nnz+1,iys:iye,mxs:mxe),
     +                partTsrc(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      !Try only calling these when the history data is being written:
      if(mtrans  .and. istage .eq. 3) then
      call ztox_trans(upwp_t(0:nnz+1,iys:iye,mxs:mxe),
     +                upwp(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(vpsum_t(0:nnz+1,iys:iye,mxs:mxe,1),
     +                vpsum(1:nnx,iys:iye,izs-1:ize+1,1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(vpsum_t(0:nnz+1,iys:iye,mxs:mxe,2),
     +                vpsum(1:nnx,iys:iye,izs-1:ize+1,2),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(vpsum_t(0:nnz+1,iys:iye,mxs:mxe,3),
     +                vpsum(1:nnx,iys:iye,izs-1:ize+1,3),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(vpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe,1),
     +                vpsqrsum(1:nnx,iys:iye,izs-1:ize+1,1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(vpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe,2),
     +                vpsqrsum(1:nnx,iys:iye,izs-1:ize+1,2),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(vpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe,3),
     +                vpsqrsum(1:nnx,iys:iye,izs-1:ize+1,3),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(Tpsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +                Tpsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(Tpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +                Tpsqrsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(wpTpsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +                wpTpsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(partcount_t(0:nnz+1,iys:iye,mxs:mxe),
     +                partcount(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(partcountU_t(0:nnz+1,iys:iye,mxs:mxe),
     +                partcountU(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(partflux_t(0:nnz+1,iys:iye,mxs:mxe),
     +                partflux(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(partfluxup_t(0:nnz+1,iys:iye,mxs:mxe),
     +                partfluxup(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(partfluxdn_t(0:nnz+1,iys:iye,mxs:mxe),
     +                partfluxdn(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)

      call ztox_trans(Tfsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +                Tfsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(qfsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +                qfsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(radsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +                radsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(rad2sum_t(0:nnz+1,iys:iye,mxs:mxe),
     +                rad2sum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      call ztox_trans(qstarsum_t(0:nnz+1,iys:iye,mxs:mxe),
     +                qstarsum(1:nnx,iys:iye,izs-1:ize+1),nnx,nnz,mxs,
     +                mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,myid,
     +                ncpu_s,numprocs)
      end if

      !t_s = mpi_wtime
      !Get particle count:
      numpart = 0
      part => first_particle
      do while (associated(part))
      numpart = numpart + 1
      part => part%next
      end do
      !call mpi_barrier(mpi_comm_world,ierr)
      !t_f = mpi_wtime()
      !if (myid==5) write(*,*) 'time numpart: ', t_f - t_s

      !t_s = mpi_wtime()
      !Compute total number of particles
      call mpi_allreduce(numpart,tnumpart,1,mpi_integer,mpi_sum,
     +                   mpi_comm_world,ierr)
      !Compute average particle Reynolds number
      call mpi_allreduce(myRep_avg,Rep_avg,1,mpi_real8,mpi_sum,
     +                   mpi_comm_world,ierr)
      Rep_avg = Rep_avg/tnumpart
      !call mpi_barrier(mpi_comm_world,ierr)
      !t_f = mpi_wtime()
      !if (myid==5) write(*,*) 'time mpi_allreduce: ', t_f - t_s

      end subroutine particle_update_rk3
      subroutine particle_bcs_nonperiodic
      use particles
      use con_data
      use pars
      use tank_pars
      implicit none
      real :: top,bot
      real :: restime,beg_phase,end_phase
      real :: t_bndy,wave_local


      t_bndy = t_stage - t_zero

      part => first_particle
      do while (associated(part))

      !Remove particles which fall under the radius
      !Rebound particles which hit top

      top = zl - radius_init
      bot = 0.0 + radius_init

      if (part%xip(3) .GT. top) then
         part%xp(3) = top - (part%xp(3)-top)
         part%vp(3) = -part%vp(3)
         part => part%next
      elseif (part%xip(3) .LT. bot) then

         restime = t_stage - part%t0

         !Need to know where the wave is relative to periodic box:
         wave_local = mod(c_speed*(part%t0-t_zero),xl)
         if (part%xp_0(1) .gt. wave_local) then
         !Particle ahead of wave phase inside domain
         beg_phase = mod(part%xp_0(1)-wave_local,wave_l)/wave_l
         else
         !Particle behind wave phase inside domain
         beg_phase = mod((part%xp_0(1)+xl)-wave_local,wave_l)/wave_l
         end if


         wave_local = mod(c_speed*t_bndy,xl)
         if (part%xp(1) .gt. wave_local) then
         end_phase = mod(part%xp(1)-wave_local,wave_l)/wave_l
         else
         end_phase = mod((part%xp(1)+xl)-wave_local,wave_l)/wave_l
         end if

         !write(*,'(a5,2e15.6)') 'loc:',beg_phase,end_phase
         !write(*,'(a5,2e15.6)') 'res:',hist_restime,beg_phase


         !JPDF of residence time and wave-relative origin
         call add_histogram(bins(:,1),bins(:,2),hist_restime,
     +    histbins+2,log10(restime),beg_phase)

         !JPDF of wave-relative origin and wave-relative end
         call add_histogram(bins(:,2),bins(:,3),hist_locations,
     +    histbins+2,beg_phase,end_phase)

         !JPDF of wave-relative origin and max height
         call add_histogram(bins(:,2),bins(:,4),hist_height,
     +    histbins+2,beg_phase,part%zpmax)

         !JPDF of wave-relative origin and max height
         call add_histogram(bins(:,2),bins(:,5),hist_momentum,
     +    histbins+2,beg_phase,part%vp(1)-part%u0)

         call destroy_particle
         num_reintro = num_reintro + 1

      else
         part => part%next
      end if

      end do

      end subroutine particle_bcs_nonperiodic
      subroutine particle_reintro(istage)
      use particles
      use pars
      use con_data
      use tank_pars
      implicit none
      include 'mpif.h'

      integer :: ierr,randproc,np,tmp_reintro,istage,ix
      real :: xv,yv,zv,ran2,vinj,t_bndy,uorb,worb
      real :: xpart(nnx),y2part(nnx),yp1,ypn,ddx
      real :: u2part(nnx),w2part(nnx)

      t_bndy = t_stage - t_zero

      call mpi_allreduce(num_reintro,tmp_reintro,1,mpi_integer,mpi_sum,
     +       mpi_comm_world,ierr)

      if (myid==0) write(*,*) 'time,tot_reintro: ',time,tmp_reintro

      tot_reintro = tot_reintro + tmp_reintro

      !Set up xpart for interpolation of orbital speeds:
      yp1 = 2.0e+30
      ypn = yp1
      ddx = xl/float(nnx)
      do ix=1,nnx
         xpart(ix) = float(ix-1)*ddx
      end do
      ! call spline(xpart,u_tank_i,nnx,yp1,ypn,u2part)
      ! call spline(xpart,w_tank_i,nnx,yp1,ypn,w2part)

      !Now everyone knows how many total particles need to be reintroduced
      !For each one to be reintroduced, get a random processor who creates it:

      do np=1,tmp_reintro

      !Proc 0 gets a random proc ID, broadcasts it out:
      if (myid==0) randproc = floor(ran2(iseed)*numprocs)

      call mpi_bcast(randproc,1,mpi_integer,0,mpi_comm_world,ierr)

       if (myid==randproc) then

        xv = ran2(iseed)*(xmax-xmin) + xmin
        yv = ran2(iseed)*(ymax-ymin) + ymin

        !Inject the particles upwards:
        !Include local orbital velocity as well as a random component for ejection
        ! uorb = waveslope*c_speed*cos(pi2/wave_l*(xv-c_speed*t_bndy))
        ! worb = -waveslope*c_speed*sin(pi2/wave_l*(xv-c_speed*t_bndy))
        call linInterp(xpart,u_tank_i,nnx,xv,uorb)
        call linInterp(xpart,w_tank_i,nnx,xv,worb)

        vinj = worb + ran2(iseed)*Vpmax
        !zv = radius_init
        zv = zw1

        call create_particle((/xv,yv,zv/),(/uorb,0.0,vinj/),300.0,
     +   radius_init,2,z_tank_i(1:nnx),xpart(1:nnx),
     +   y2part(1:nnx),t_stage)

       end if
      end do

      end subroutine particle_reintro
      subroutine particle_bcs_periodic
      use particles
      use con_data
      use pars
      implicit none

      !Assumes domain goes from [0,xl),[0,yl),[0,zl]
      !Also maintain the number of particles on each proc

      part => first_particle
      do while (associated(part))

      !x,y periodic

      if (part%xp(1) .GT. xl) then
         part%xp(1) = part%xp(1)-xl
      elseif (part%xp(1) .LT. 0) then
         part%xp(1) = xl + part%xp(1)
      end if

      if (part%xp(2) .GT. yl) then
         part%xp(2) = part%xp(2)-yl
      elseif (part%xp(2) .LT. 0) then
         part%xp(2) = yl + part%xp(2)
      end if

      part => part%next

      end do


      end subroutine particle_bcs_periodic
      subroutine particle_exchange
      use pars
      use particles
      use con_data
      use con_stats
      implicit none
      include 'mpif.h'

      type(particle), pointer :: tmp
      integer :: idx,psum,csum
      integer :: ir,itr,itop,itl,il,ibl,ib,ibr
      integer :: istatus(mpi_status_size),ierr
      integer :: status_array(mpi_status_size,16),req(16)
      type(particle), allocatable :: rbuf_s(:),trbuf_s(:)
      type(particle), allocatable :: tbuf_s(:),tlbuf_s(:)
      type(particle), allocatable :: lbuf_s(:),blbuf_s(:)
      type(particle), allocatable :: bbuf_s(:),brbuf_s(:)
      type(particle), allocatable :: rbuf_r(:),trbuf_r(:)
      type(particle), allocatable :: tbuf_r(:),tlbuf_r(:)
      type(particle), allocatable :: lbuf_r(:),blbuf_r(:)
      type(particle), allocatable :: bbuf_r(:),brbuf_r(:)
      type(particle), allocatable :: totalbuf(:)

      !Zero out the counters for how many particles to send each dir.
      pr_s=0;ptr_s=0;pt_s=0;ptl_s=0;pl_s=0;pbl_s=0;pb_s=0;pbr_s=0

      !As soon as the location is updated, must check to see if it left the proc:
      !May be a better way of doing this, but it seems most reasonable:
      part => first_particle
      do while (associated(part))

         !First get numbers being sent to all sides:
         if (part%xp(2) .GT. ymax) then
            if (part%xp(1) .GT. xmax) then !top right
               ptr_s = ptr_s + 1
            elseif (part%xp(1) .LT. xmin) then !bottom right
               pbr_s = pbr_s + 1
            else  !right
               pr_s = pr_s + 1
            end if
         elseif (part%xp(2) .LT. ymin) then
            if (part%xp(1) .GT. xmax) then !top left
               ptl_s = ptl_s + 1
            else if (part%xp(1) .LT. xmin) then !bottom left
               pbl_s = pbl_s + 1
            else  !left
               pl_s = pl_s + 1
            end if
         elseif ( (part%xp(1) .GT. xmax) .AND.
     +           (part%xp(2) .LT. ymax) .AND.
     +           (part%xp(2) .GT. ymin) ) then !top
            pt_s = pt_s + 1
         elseif ( (part%xp(1) .LT. xmin) .AND.
     +           (part%xp(2) .LT. ymax) .AND.
     +           (part%xp(2) .GT. ymin) ) then !bottom
            pb_s = pb_s + 1
         end if

         part => part%next
      end do

      !Now allocate the send buffers based on these counts:
      allocate(rbuf_s(pr_s),trbuf_s(ptr_s),tbuf_s(pt_s),tlbuf_s(ptl_s))
      allocate(lbuf_s(pl_s),blbuf_s(pbl_s),bbuf_s(pb_s),brbuf_s(pbr_s))

      !Now loop back through the particles and fill the buffers:
      !NOTE: If it finds one, add it to buffer and REMOVE from list
      ir=1;itr=1;itop=1;itl=1;il=1;ibl=1;ib=1;ibr=1

      part => first_particle
      do while (associated(part))

         if (part%xp(2) .GT. ymax) then
            if (part%xp(1) .GT. xmax) then !top right
!               write(*,*) 'Proc',myid,'about to send part to
!     +proc',trproc,'TOP RIGHT'
               trbuf_s(itr) = part
               call destroy_particle
               itr = itr + 1
            elseif (part%xp(1) .LT. xmin) then !bottom right
!               write(*,*) 'Proc',myid,'about to send part to
!     +proc',brproc,'BOTTOM RIGHT'
               brbuf_s(ibr) = part
               call destroy_particle
               ibr = ibr + 1
            else   !right
!               write(*,*) 'Proc',myid,'about to send part to
!     +proc',rproc,'RIGHT'
               rbuf_s(ir) = part
               call destroy_particle
               ir = ir + 1
            end if
         elseif (part%xp(2) .LT. ymin) then
            if (part%xp(1) .GT. xmax) then !top left
!               write(*,*) 'Proc',myid,'about to send part to
!     +proc',tlproc,'TOP LEFT'
               tlbuf_s(itl) = part
               call destroy_particle
               itl = itl + 1
            else if (part%xp(1) .LT. xmin) then !bottom left
!               write(*,*) 'Proc',myid,'about to send part to
!     +proc',blproc,'BOTTOM LEFT'
               blbuf_s(ibl) = part
               call destroy_particle
               ibl = ibl + 1
            else  !left
!               write(*,*) 'Proc',myid,'about to send part to
!     +proc',lproc,'LEFT'
               lbuf_s(il) = part
               call destroy_particle
               il = il + 1
            end if
         elseif ( (part%xp(1) .GT. xmax) .AND.
     +           (part%xp(2) .LT. ymax) .AND.
     +           (part%xp(2) .GT. ymin) ) then !top
!               write(*,*) 'Proc',myid,'about to send part to
!     +proc',tproc,'TOP'
            tbuf_s(itop) = part
            call destroy_particle
            itop = itop + 1
         elseif ( (part%xp(1) .LT. xmin) .AND.
     +           (part%xp(2) .LT. ymax) .AND.
     +           (part%xp(2) .GT. ymin) ) then !bottom
!               write(*,*) 'Proc',myid,'about to send part to
!     +proc',bproc,'BOTTOM'
            bbuf_s(ib) = part
            call destroy_particle
            ib = ib + 1
         else
         part => part%next
         end if

      end do

      !Now everyone exchanges the counts with all neighbors:
      !Left/right:
      call MPI_Sendrecv(pr_s,1,mpi_integer,rproc,3,
     +        pl_r,1,mpi_integer,lproc,3,mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(pl_s,1,mpi_integer,lproc,4,
     +        pr_r,1,mpi_integer,rproc,4,mpi_comm_world,istatus,ierr)

      !Top/bottom:
      call MPI_Sendrecv(pt_s,1,mpi_integer,tproc,5,
     +        pb_r,1,mpi_integer,bproc,5,mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(pb_s,1,mpi_integer,bproc,6,
     +        pt_r,1,mpi_integer,tproc,6,mpi_comm_world,istatus,ierr)

      !Top right/bottom left:
      call MPI_Sendrecv(ptr_s,1,mpi_integer,trproc,7,
     +        pbl_r,1,mpi_integer,blproc,7,
     +        mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(pbl_s,1,mpi_integer,blproc,8,
     +        ptr_r,1,mpi_integer,trproc,8,
     +        mpi_comm_world,istatus,ierr)

       !Top left/bottom right:
      call MPI_Sendrecv(ptl_s,1,mpi_integer,tlproc,9,
     +        pbr_r,1,mpi_integer,brproc,9,
     +        mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(pbr_s,1,mpi_integer,brproc,10,
     +         ptl_r,1,mpi_integer,tlproc,10,
     +         mpi_comm_world,istatus,ierr)

      !Now everyone has the number of particles arriving from every neighbor
      !If the count is greater than zero, exchange:

      !Allocate room to receive from each side
      allocate(rbuf_r(pr_r),trbuf_r(ptr_r),tbuf_r(pt_r),tlbuf_r(ptl_r))
      allocate(lbuf_r(pl_r),blbuf_r(pbl_r),bbuf_r(pb_r),brbuf_r(pbr_r))

      !Send to right:
      if (pr_s .GT. 0) then
      call mpi_isend(rbuf_s,pr_s,particletype,rproc,11,
     +               mpi_comm_world,req(1),ierr)
      else
      req(1) = mpi_request_null
      end if

      !Receive from left:
      if (pl_r .GT. 0) then
      call mpi_irecv(lbuf_r,pl_r,particletype,lproc,11,
     +               mpi_comm_world,req(2),ierr)
      else
      req(2) = mpi_request_null
      end if

      !Send to left:
      if (pl_s .GT. 0) then
      call mpi_isend(lbuf_s,pl_s,particletype,lproc,12,
     +               mpi_comm_world,req(3),ierr)
      else
      req(3) = mpi_request_null
      end if

      !Receive from right:
      if (pr_r .GT. 0) then
      call mpi_irecv(rbuf_r,pr_r,particletype,rproc,12,
     +               mpi_comm_world,req(4),ierr)
      else
      req(4) = mpi_request_null
      end if

      !Send to top:
      if (pt_s .GT. 0) then
      call mpi_isend(tbuf_s,pt_s,particletype,tproc,13,
     +                mpi_comm_world,req(5),ierr)
      else
      req(5) = mpi_request_null
      end if

      !Receive from bottom:
      if (pb_r .GT. 0) then
      call mpi_irecv(bbuf_r,pb_r,particletype,bproc,13,
     +                mpi_comm_world,req(6),ierr)
      else
      req(6) = mpi_request_null
      end if

      !Send to bottom:
      if (pb_s .GT. 0) then
      call mpi_isend(bbuf_s,pb_s,particletype,bproc,14,
     +                mpi_comm_world,req(7),ierr)
      else
      req(7) = mpi_request_null
      end if

      !Recieve from top:
      if (pt_r .GT. 0) then
      call mpi_irecv(tbuf_r,pt_r,particletype,tproc,14,
     +                mpi_comm_world,req(8),ierr)
      else
      req(8) = mpi_request_null
      end if

      !Send to top right:
      if (ptr_s .GT. 0) then
      call mpi_isend(trbuf_s,ptr_s,particletype,trproc,15,
     +                mpi_comm_world,req(9),ierr)
      else
      req(9) = mpi_request_null
      end if

      !Receive from bottom left:
      if (pbl_r .GT. 0) then
      call mpi_irecv(blbuf_r,pbl_r,particletype,blproc,15,
     +                mpi_comm_world,req(10),ierr)
      else
      req(10) = mpi_request_null
      end if

      !Send to bottom left:
      if (pbl_s .GT. 0) then
      call mpi_isend(blbuf_s,pbl_s,particletype,blproc,16,
     +                mpi_comm_world,req(11),ierr)
      else
      req(11) = mpi_request_null
      end if

      !Receive from top right:
      if (ptr_r .GT. 0) then
      call mpi_irecv(trbuf_r,ptr_r,particletype,trproc,16,
     +                mpi_comm_world,req(12),ierr)
      else
      req(12) = mpi_request_null
      end if

      !Send to top left:
      if (ptl_s .GT. 0) then
      call mpi_isend(tlbuf_s,ptl_s,particletype,tlproc,17,
     +                mpi_comm_world,req(13),ierr)
      else
      req(13) = mpi_request_null
      end if

      !Receive from bottom right:
      if (pbr_r .GT. 0) then
      call mpi_irecv(brbuf_r,pbr_r,particletype,brproc,17,
     +                mpi_comm_world,req(14),ierr)
      else
      req(14) = mpi_request_null
      end if

      !Send to bottom right:
      if (pbr_s .GT. 0) then
      call mpi_isend(brbuf_s,pbr_s,particletype,brproc,18,
     +                mpi_comm_world,req(15),ierr)
      else
      req(15) = mpi_request_null
      end if

      !Receive from top left:
      if (ptl_r .GT. 0) then
      call mpi_irecv(tlbuf_r,ptl_r,particletype,tlproc,18,
     +                mpi_comm_world,req(16),ierr)
      else
      req(16) = mpi_request_null
      end if

      call mpi_waitall(16,req,status_array,ierr)

      !Now add incoming particles to linked list:
      !NOTE: add them to beginning since it's easiest to access (first_particle)

      !Form one large buffer to loop through and add:
      psum = pr_r+ptr_r+pt_r+ptl_r+pl_r+pbl_r+pb_r+pbr_r
      csum = 0
      allocate(totalbuf(psum))
      if (pr_r .GT. 0) then
         totalbuf(1:pr_r) = rbuf_r(1:pr_r)
         csum = csum + pr_r
      end if
      if (ptr_r .GT. 0) then
         totalbuf(csum+1:csum+ptr_r) = trbuf_r(1:ptr_r)
         csum = csum + ptr_r
      end if
      if (pt_r .GT. 0) then
         totalbuf(csum+1:csum+pt_r) = tbuf_r(1:pt_r)
         csum = csum + pt_r
      end if
      if (ptl_r .GT. 0) then
         totalbuf(csum+1:csum+ptl_r) = tlbuf_r(1:ptl_r)
         csum = csum + ptl_r
      end if
      if (pl_r .GT. 0) then
         totalbuf(csum+1:csum+pl_r) = lbuf_r(1:pl_r)
         csum = csum + pl_r
      end if
      if (pbl_r .GT. 0) then
         totalbuf(csum+1:csum+pbl_r) = blbuf_r(1:pbl_r)
         csum = csum + pbl_r
      end if
      if (pb_r .GT. 0) then
         totalbuf(csum+1:csum+pb_r) = bbuf_r(1:pb_r)
         csum = csum + pb_r
      end if
      if (pbr_r .GT. 0) then
         totalbuf(csum+1:csum+pbr_r) = brbuf_r(1:pbr_r)
         csum = csum + pbr_r
      end if

      do idx = 1,psum
        if (.NOT. associated(first_particle)) then
           allocate(first_particle)
           first_particle = totalbuf(idx)
           nullify(first_particle%next,first_particle%prev)
        else
           allocate(first_particle%prev)
           tmp => first_particle%prev
           tmp = totalbuf(idx)
           tmp%next => first_particle
           nullify(tmp%prev)
           first_particle => tmp
           nullify(tmp)
        end if
      end do

      deallocate(rbuf_s,trbuf_s,tbuf_s,tlbuf_s)
      deallocate(lbuf_s,blbuf_s,bbuf_s,brbuf_s)
      deallocate(rbuf_r,trbuf_r,tbuf_r,tlbuf_r)
      deallocate(lbuf_r,blbuf_r,bbuf_r,brbuf_r)
      deallocate(totalbuf)

      end subroutine particle_exchange
      subroutine destroy_particle
      use particles
      implicit none

      type(particle), pointer :: tmp

      !Is it the first and last in the list?
      if (associated(part,first_particle) .AND.
     +    (.NOT. associated(part%next)) ) then
          nullify(first_particle)
          deallocate(part)
      else
        if (associated(part,first_particle)) then !Is it the first particle?
           first_particle => part%next
           part => first_particle
           deallocate(part%prev)
        elseif (.NOT. associated(part%next)) then !Is it the last particle?
           nullify(part%prev%next)
           deallocate(part)
        else
           tmp => part
           part => part%next
           tmp%prev%next => tmp%next
           tmp%next%prev => tmp%prev
           deallocate(tmp)
        end if
      end if

      end subroutine destroy_particle

      subroutine fill_ext
      use pars
      use particles
      use fields
      use con_stats
      use con_data
      implicit none
      include 'mpif.h'

      integer :: istatus(mpi_status_size),ierr
      integer :: ix,iy,iz
      !preceding letter: r=right,l=left,t=top,b=bot.
      !_s: buf of things to send TO r,l,t,b
      !_r: buf of things to recv FROM r,l,t,b
      real :: tbuf_s(nnz+2,iye-iys+1,2,5),tbuf_r(nnz+2,iye-iys+1,3,5)
      real :: bbuf_s(nnz+2,iye-iys+1,3,5),bbuf_r(nnz+2,iye-iys+1,2,5)
      real :: rbuf_s(nnz+2,2,mxe-mxs+1,5),rbuf_r(nnz+2,3,mxe-mxs+1,5)
      real :: lbuf_s(nnz+2,3,mxe-mxs+1,5),lbuf_r(nnz+2,2,mxe-mxs+1,5)

      !Corners:
      real :: trbuf_s(nnz+2,2,2,5),trbuf_r(nnz+2,3,3,5)
      real :: brbuf_s(nnz+2,2,3,5),brbuf_r(nnz+2,3,2,5)
      real :: blbuf_s(nnz+2,3,3,5),blbuf_r(nnz+2,2,2,5)
      real :: tlbuf_s(nnz+2,3,2,5),tlbuf_r(nnz+2,2,3,5)
      !MPI send counts:
      integer :: rc_s,rc_r,trc_s,trc_r,tc_s,tc_r,tlc_s,tlc_r
      integer :: lc_s,lc_r,blc_s,blc_r,bc_s,bc_r,brc_s,brc_r

      !Debugging:
      real :: xv,yv,zv

      !To update the particle ODE in time, need the interpolated
      !velocity field
      !This requires filling uext,vext,wext from nearby procs
      uext = 0.0
      vext = 0.0
      wext = 0.0
      Text = 0.0
      T2ext = 0.0

      !First fill the center, since this is just u,v,w on that proc:

      !In the column setup, need to tranpose u,v,w first into u_t,v_t,w_t:
      call xtoz_trans(u(1:nnx,iys:iye,izs-1:ize+1),u_t,nnx,nnz,
     +                mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,
     +                myid,ncpu_s,numprocs)
      call xtoz_trans(v(1:nnx,iys:iye,izs-1:ize+1),v_t,nnx,nnz,
     +                mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,
     +                myid,ncpu_s,numprocs)
      call xtoz_trans(w(1:nnx,iys:iye,izs-1:ize+1),w_t,nnx,nnz,
     +                mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,
     +                myid,ncpu_s,numprocs)
      call xtoz_trans(t(1:nnx,iys:iye,1,izs-1:ize+1),T_t,nnx,nnz,
     +                mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,
     +                myid,ncpu_s,numprocs)
      call xtoz_trans(t(1:nnx,iys:iye,2,izs-1:ize+1),T2_t,nnx,nnz,
     +                mxs,mxe,mx_s,mx_e,iys,iye,izs,ize,iz_s,iz_e,
     +                myid,ncpu_s,numprocs)

      uext(0:nnz+1,iys:iye,mxs:mxe) = u_t(0:nnz+1,iys:iye,mxs:mxe)
      vext(0:nnz+1,iys:iye,mxs:mxe) = v_t(0:nnz+1,iys:iye,mxs:mxe)
      wext(0:nnz+1,iys:iye,mxs:mxe) = w_t(0:nnz+1,iys:iye,mxs:mxe)
      Text(0:nnz+1,iys:iye,mxs:mxe) = T_t(0:nnz+1,iys:iye,mxs:mxe)
      T2ext(0:nnz+1,iys:iye,mxs:mxe) = T2_t(0:nnz+1,iys:iye,mxs:mxe)

      !Recall that SR assign_nbrs assigned rproc,lproc, etc.

      !Going to call 6 sendrecv calls - one for each proc. nbr.:

      !Fill the send buffers:

      !I know these are redundant, but so I can keep them straight...
      tc_s = 5*(nnz+2)*2*(iye-iys+1)
      tc_r = 5*(nnz+2)*3*(iye-iys+1)
      trc_s = 5*(nnz+2)*2*2
      trc_r = 5*(nnz+2)*3*3
      rc_s = 5*(nnz+2)*(mxe-mxs+1)*2
      rc_r = 5*(nnz+2)*(mxe-mxs+1)*3
      tlc_s = 5*(nnz+2)*3*2
      tlc_r = 5*(nnz+2)*2*3
      bc_s = 5*(nnz+2)*3*(iye-iys+1)
      bc_r = 5*(nnz+2)*2*(iye-iys+1)
      blc_s = 5*(nnz+2)*3*3
      blc_r = 5*(nnz+2)*2*2
      lc_s = 5*(nnz+2)*(mxe-mxs+1)*3
      lc_r = 5*(nnz+2)*(mxe-mxs+1)*2
      brc_s = 5*(nnz+2)*2*3
      brc_r = 5*(nnz+2)*3*2

      !First u:
      tbuf_s(1:nnz+2,1:iye-iys+1,1:2,1) = u_t(0:nnz+1,iys:iye,mxe-1:mxe)
      trbuf_s(1:nnz+2,1:2,1:2,1) = u_t(0:nnz+1,iye-1:iye,mxe-1:mxe)
      rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,1) = u_t(0:nnz+1,iye-1:iye,mxs:mxe)
      brbuf_s(1:nnz+2,1:2,1:3,1) = u_t(0:nnz+1,iye-1:iye,mxs:mxs+2)
      bbuf_s(1:nnz+2,1:iye-iys+1,1:3,1) = u_t(0:nnz+1,iys:iye,mxs:mxs+2)
      blbuf_s(1:nnz+2,1:3,1:3,1) = u_t(0:nnz+1,iys:iys+2,mxs:mxs+2)
      lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,1) = u_t(0:nnz+1,iys:iys+2,mxs:mxe)
      tlbuf_s(1:nnz+2,1:3,1:2,1) = u_t(0:nnz+1,iys:iys+2,mxe-1:mxe)

      !v:
      tbuf_s(1:nnz+2,1:iye-iys+1,1:2,2) = v_t(0:nnz+1,iys:iye,mxe-1:mxe)
      trbuf_s(1:nnz+2,1:2,1:2,2) = v_t(0:nnz+1,iye-1:iye,mxe-1:mxe)
      rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,2) = v_t(0:nnz+1,iye-1:iye,mxs:mxe)
      brbuf_s(1:nnz+2,1:2,1:3,2) = v_t(0:nnz+1,iye-1:iye,mxs:mxs+2)
      bbuf_s(1:nnz+2,1:iye-iys+1,1:3,2) = v_t(0:nnz+1,iys:iye,mxs:mxs+2)
      blbuf_s(1:nnz+2,1:3,1:3,2) = v_t(0:nnz+1,iys:iys+2,mxs:mxs+2)
      lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,2) = v_t(0:nnz+1,iys:iys+2,mxs:mxe)
      tlbuf_s(1:nnz+2,1:3,1:2,2) = v_t(0:nnz+1,iys:iys+2,mxe-1:mxe)

      !w:
      tbuf_s(1:nnz+2,1:iye-iys+1,1:2,3) = w_t(0:nnz+1,iys:iye,mxe-1:mxe)
      trbuf_s(1:nnz+2,1:2,1:2,3) = w_t(0:nnz+1,iye-1:iye,mxe-1:mxe)
      rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,3) = w_t(0:nnz+1,iye-1:iye,mxs:mxe)
      brbuf_s(1:nnz+2,1:2,1:3,3) = w_t(0:nnz+1,iye-1:iye,mxs:mxs+2)
      bbuf_s(1:nnz+2,1:iye-iys+1,1:3,3) = w_t(0:nnz+1,iys:iye,mxs:mxs+2)
      blbuf_s(1:nnz+2,1:3,1:3,3) = w_t(0:nnz+1,iys:iys+2,mxs:mxs+2)
      lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,3) = w_t(0:nnz+1,iys:iys+2,mxs:mxe)
      tlbuf_s(1:nnz+2,1:3,1:2,3) = w_t(0:nnz+1,iys:iys+2,mxe-1:mxe)

      !T:
      tbuf_s(1:nnz+2,1:iye-iys+1,1:2,4) = T_t(0:nnz+1,iys:iye,mxe-1:mxe)
      trbuf_s(1:nnz+2,1:2,1:2,4) = T_t(0:nnz+1,iye-1:iye,mxe-1:mxe)
      rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,4) = T_t(0:nnz+1,iye-1:iye,mxs:mxe)
      brbuf_s(1:nnz+2,1:2,1:3,4) = T_t(0:nnz+1,iye-1:iye,mxs:mxs+2)
      bbuf_s(1:nnz+2,1:iye-iys+1,1:3,4) = T_t(0:nnz+1,iys:iye,mxs:mxs+2)
      blbuf_s(1:nnz+2,1:3,1:3,4) = T_t(0:nnz+1,iys:iys+2,mxs:mxs+2)
      lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,4) = T_t(0:nnz+1,iys:iys+2,mxs:mxe)
      tlbuf_s(1:nnz+2,1:3,1:2,4) = T_t(0:nnz+1,iys:iys+2,mxe-1:mxe)

      !T2:
      tbuf_s(1:nnz+2,1:iye-iys+1,1:2,5) =T2_t(0:nnz+1,iys:iye,mxe-1:mxe)
      trbuf_s(1:nnz+2,1:2,1:2,5) = T2_t(0:nnz+1,iye-1:iye,mxe-1:mxe)
      rbuf_s(1:nnz+2,1:2,1:mxe-mxs+1,5) =T2_t(0:nnz+1,iye-1:iye,mxs:mxe)
      brbuf_s(1:nnz+2,1:2,1:3,5) = T2_t(0:nnz+1,iye-1:iye,mxs:mxs+2)
      bbuf_s(1:nnz+2,1:iye-iys+1,1:3,5) =T2_t(0:nnz+1,iys:iye,mxs:mxs+2)
      blbuf_s(1:nnz+2,1:3,1:3,5) = T2_t(0:nnz+1,iys:iys+2,mxs:mxs+2)
      lbuf_s(1:nnz+2,1:3,1:mxe-mxs+1,5) =T2_t(0:nnz+1,iys:iys+2,mxs:mxe)
      tlbuf_s(1:nnz+2,1:3,1:2,5) = T2_t(0:nnz+1,iys:iys+2,mxe-1:mxe)

      !Zero out recieve buffers
      rbuf_r=0.0;trbuf_r=0.0;tbuf_r=0.0;tlbuf_r=0.0;lbuf_r=0.0
      blbuf_r=0.0;bbuf_r=0.0;brbuf_r=0.0

      !Left/right:
      call MPI_Sendrecv(rbuf_s,rc_s,mpi_real8,rproc,3,
     +        lbuf_r,lc_r,mpi_real8,lproc,3,mpi_comm_world,istatus,ierr)

      call mpi_barrier(mpi_comm_world,ierr)
      call MPI_Sendrecv(lbuf_s,lc_s,mpi_real8,lproc,4,
     +        rbuf_r,rc_r,mpi_real8,rproc,4,mpi_comm_world,istatus,ierr)

      !Top/bottom:
      call MPI_Sendrecv(tbuf_s,tc_s,mpi_real8,tproc,5,
     +        bbuf_r,bc_r,mpi_real8,bproc,5,mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(bbuf_s,bc_s,mpi_real8,bproc,6,
     +        tbuf_r,tc_r,mpi_real8,tproc,6,mpi_comm_world,istatus,ierr)

      !Top right/bottom left:
      call MPI_Sendrecv(trbuf_s,trc_s,mpi_real8,trproc,7,
     +        blbuf_r,blc_r,mpi_real8,blproc,7,
     +        mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(blbuf_s,blc_s,mpi_real8,blproc,8,
     +        trbuf_r,trc_r,mpi_real8,trproc,8,
     +        mpi_comm_world,istatus,ierr)

       !Top left/bottom right:
      call MPI_Sendrecv(tlbuf_s,tlc_s,mpi_real8,tlproc,9,
     +        brbuf_r,brc_r,mpi_real8,brproc,9,
     +        mpi_comm_world,istatus,ierr)

      call MPI_Sendrecv(brbuf_s,brc_s,mpi_real8,brproc,10,
     +         tlbuf_r,tlc_r,mpi_real8,tlproc,10,
     +         mpi_comm_world,istatus,ierr)

      !Now fill the ext arrays with the recieved buffers:
      uext(0:nnz+1,iys:iye,mxe+1:mxe+3) =
     +     tbuf_r(1:nnz+2,1:iye-iys+1,1:3,1)
      uext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = trbuf_r(1:nnz+2,1:3,1:3,1)
      uext(0:nnz+1,iye+1:iye+3,mxs:mxe) =
     +     rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,1)
      uext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = brbuf_r(1:nnz+2,1:3,1:2,1)
      uext(0:nnz+1,iys:iye,mxs-2:mxs-1) =
     +     bbuf_r(1:nnz+2,1:iye-iys+1,1:2,1)
      uext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = blbuf_r(1:nnz+2,1:2,1:2,1)
      uext(0:nnz+1,iys-2:iys-1,mxs:mxe) =
     +     lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,1)
      uext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = tlbuf_r(1:nnz+2,1:2,1:3,1)

      vext(0:nnz+1,iys:iye,mxe+1:mxe+3) =
     +     tbuf_r(1:nnz+2,1:iye-iys+1,1:3,2)
      vext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = trbuf_r(1:nnz+2,1:3,1:3,2)
      vext(0:nnz+1,iye+1:iye+3,mxs:mxe) =
     +     rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,2)
      vext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = brbuf_r(1:nnz+2,1:3,1:2,2)
      vext(0:nnz+1,iys:iye,mxs-2:mxs-1) =
     +     bbuf_r(1:nnz+2,1:iye-iys+1,1:2,2)
      vext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = blbuf_r(1:nnz+2,1:2,1:2,2)
      vext(0:nnz+1,iys-2:iys-1,mxs:mxe) =
     +     lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,2)
      vext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = tlbuf_r(1:nnz+2,1:2,1:3,2)

      wext(0:nnz+1,iys:iye,mxe+1:mxe+3) =
     +     tbuf_r(1:nnz+2,1:iye-iys+1,1:3,3)
      wext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = trbuf_r(1:nnz+2,1:3,1:3,3)
      wext(0:nnz+1,iye+1:iye+3,mxs:mxe) =
     +     rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,3)
      wext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = brbuf_r(1:nnz+2,1:3,1:2,3)
      wext(0:nnz+1,iys:iye,mxs-2:mxs-1) =
     +     bbuf_r(1:nnz+2,1:iye-iys+1,1:2,3)
      wext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = blbuf_r(1:nnz+2,1:2,1:2,3)
      wext(0:nnz+1,iys-2:iys-1,mxs:mxe) =
     +     lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,3)
      wext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = tlbuf_r(1:nnz+2,1:2,1:3,3)

      Text(0:nnz+1,iys:iye,mxe+1:mxe+3) =
     +     tbuf_r(1:nnz+2,1:iye-iys+1,1:3,4)
      Text(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) = trbuf_r(1:nnz+2,1:3,1:3,4)
      Text(0:nnz+1,iye+1:iye+3,mxs:mxe) =
     +     rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,4)
      Text(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) = brbuf_r(1:nnz+2,1:3,1:2,4)
      Text(0:nnz+1,iys:iye,mxs-2:mxs-1) =
     +     bbuf_r(1:nnz+2,1:iye-iys+1,1:2,4)
      Text(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) = blbuf_r(1:nnz+2,1:2,1:2,4)
      Text(0:nnz+1,iys-2:iys-1,mxs:mxe) =
     +     lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,4)
      Text(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) = tlbuf_r(1:nnz+2,1:2,1:3,4)

      T2ext(0:nnz+1,iys:iye,mxe+1:mxe+3) =
     +     tbuf_r(1:nnz+2,1:iye-iys+1,1:3,5)
      T2ext(0:nnz+1,iye+1:iye+3,mxe+1:mxe+3) =trbuf_r(1:nnz+2,1:3,1:3,5)
      T2ext(0:nnz+1,iye+1:iye+3,mxs:mxe) =
     +     rbuf_r(1:nnz+2,1:3,1:mxe-mxs+1,5)
      T2ext(0:nnz+1,iye+1:iye+3,mxs-2:mxs-1) =brbuf_r(1:nnz+2,1:3,1:2,5)
      T2ext(0:nnz+1,iys:iye,mxs-2:mxs-1) =
     +     bbuf_r(1:nnz+2,1:iye-iys+1,1:2,5)
      T2ext(0:nnz+1,iys-2:iys-1,mxs-2:mxs-1) =blbuf_r(1:nnz+2,1:2,1:2,5)
      T2ext(0:nnz+1,iys-2:iys-1,mxs:mxe) =
     +     lbuf_r(1:nnz+2,1:2,1:mxe-mxs+1,5)
      T2ext(0:nnz+1,iys-2:iys-1,mxe+1:mxe+3) =tlbuf_r(1:nnz+2,1:2,1:3,5)


      end subroutine fill_ext
      subroutine get_zeta_algeb(alpha,L,h,xip)
      implicit none

      real :: f,fp,h,alpha,tol,err,zeta,zetanew,L,z
      real, intent(inout) :: xip(1:3)

      err = 1.0
      z = xip(3)
      zeta = z
      do while (err .gt. 1e-10)
        f = zeta + h*(1.0 - zeta/L)**alpha - z
        fp = 1 - alpha*h/L*(1.0 - zeta/L)**(alpha-1.0)

        zetanew = zeta - f/fp
        err = abs(zeta-zetanew)/zeta
        zeta = zetanew
      end do
      xip(3) = zeta

      end subroutine get_zeta_algeb
      subroutine uf_interp(hpart)
      use pars
      use particles
      use fields
      use con_stats
      use con_data
      implicit none

      integer :: ix,iy,izuv,izw,iz,i,k,j
      integer :: first,last
      real :: xkval,xjval,pj,dxvec(2)
      integer :: ijpts(2,6),kuvpts(6)
      real :: wt(3,6)
      real :: ran2,hpart

      !Interpolates to get the fluid velocity
      !NOTE: interpolating in COMPUTATIONAL space!
      !part%xp is in physical space - first must figure out where in
      !computational space that lies and put into part%xip!
      !then interpolate velocity to that position:

      !Newton's method to back out the computational location (xi,eta,zeta)
      !into part%xip given h at that location:

      part%xip(1:3) = part%xp(1:3)
      if (iz_space == 5) then    !Algebraic stretching to smooth out waves
          call get_zeta_algeb(a_decay,zl,hpart,part%xip(1:3))
      end if

      !Following Orlando's code, get the "leftmost" node
      !This is just the minimum (i,j,k) on the volume

      ijpts(1,3) = floor(part%xip(1)/dx) + 1
      ijpts(2,3) = floor(part%xip(2)/dy) + 1

      if (ijpts(1,3) .le. 0)  then
      write(*,*) 'PROB:',part%xip(1),part%xip(2),part%xip(3)
      end if

      !Fill in the neighbors:
      ijpts(1,2) = ijpts(1,3)-1
      ijpts(1,1) = ijpts(1,2)-1
      ijpts(1,4) = ijpts(1,3)+1
      ijpts(1,5) = ijpts(1,4)+1
      ijpts(1,6) = ijpts(1,5)+1

      ijpts(2,2) = ijpts(2,3)-1
      ijpts(2,1) = ijpts(2,2)-1
      ijpts(2,4) = ijpts(2,3)+1
      ijpts(2,5) = ijpts(2,4)+1
      ijpts(2,6) = ijpts(2,5)+1

      !Finding the k-lhnode is different since grid may be stretched

      !Do (u,v) loop first:
      do iz=0,nnz+1
         if (zz(iz) .GT. part%xip(3)) go to 1199
      end do
1199  continue

      kuvpts(3) = iz-1
      !Then fill in the rest:
      kuvpts(4) = kuvpts(3)+1
      kuvpts(5) = kuvpts(4)+1
      kuvpts(6) = kuvpts(5)+1
      kuvpts(2) = kuvpts(3)-1
      kuvpts(1) = kuvpts(2)-1


!---------
!      !As an aside, use ijpts,kuvpts to update the particle numbers at
!      each cell:
!
!      !call inputs  are the nodes associated with the volume that each
!      part. lies in
!      call particle_stats(ijpts(1,3),ijpts(2,3),kuvpts(3))
!
!---------


      !Fill in the weights:
      !First for x and y since they are periodic:
      wt(1:3,1:6) = 0.0
      dxvec(1) = dx
      dxvec(2) = dy
      do iz = 1,2
      do j = 1,6
         xjval = dxvec(iz)*(ijpts(iz,j)-1)
         pj = 1.0
         do k = 1,6
            xkval = dxvec(iz)*(ijpts(iz,k)-1)
            if (j .NE. k) then
                  pj = pj*(part%xip(iz)-xkval)/(xjval-xkval)
            end if
         end do
         wt(iz,j) = pj
       end do
       end do

       !Now compute weights in z-dir
       !There are 2 sections: weights at (u,v) nodes (kuvpts)
       !And weights computed at w nodes (kwpts)

       !Compute weights at kuvpts
       !Must check to see how close we are to a top/bot boundary
       if (kuvpts(3) == 1) then
          first = 3
          last = 4
          !Set these equal to 1 so uext(-1) won't be accessed
          !Note: the value doesn't matter since weight will be 0
          kuvpts(1) = 1
          kuvpts(2) = 1
       elseif (kuvpts(3) == 0) then
          first = 4
          last = 5
          kuvpts(1) = 1
          kuvpts(2) = 1
          kuvpts(3) = 1
       elseif (kuvpts(3) .LT. 0) then
          first = 0
          last = 0
       elseif (kuvpts(3) == 2) then
          first = 2
          last = 5
       !Between top cell center and the domain boundary
       elseif (kuvpts(3) == nnz) then
          first = 2
          last = 3
          kuvpts(4) = nnz
          kuvpts(5) = nnz
          kuvpts(6) = nnz
       elseif (kuvpts(3) .GT. nnz) then
          first = 0
          last = 0
       !Between 2nd to last and last cell center at top
       elseif (kuvpts(3) == nnz-1) then
          first = 3
          last = 4
          kuvpts(5) = nnz
          kuvpts(6) = nnz
       elseif (kuvpts(3) == nnz-2) then
          first = 2
          last = 5
       else
          first = 1
          last = 6
       end if

       !Recall that wt has been set to zero, so
       !weights will be zero if (first,last) isn't (1,6)
       do j = first,last
           xjval = zz(kuvpts(j))
           pj = 1.0
           do k = first,last
              xkval = zz(kuvpts(k))
              if (j .NE. k) then
                 pj = pj*(part%xip(3)-xkval)/(xjval-xkval)
              end if
           end do
           wt(3,j) = pj
      end do

      !Now we have the weights - compute the velocity at xp:
        part%uf(1:3) = 0.0
        part%Tf = 0.0
        part%qinf = 0.0
        do k = 1,6
        do j = 1,6
        do i = 1,6
            ix = ijpts(1,i)
            iy = ijpts(2,j)
            izuv = kuvpts(k)
            part%uf(1) = part%uf(1)+uext(izuv,iy,ix)*wt(1,i)*
     +                   wt(2,j)*wt(3,k)
            part%uf(2) = part%uf(2)+vext(izuv,iy,ix)*wt(1,i)*
     +                   wt(2,j)*wt(3,k)
            part%uf(3) = part%uf(3)+wext(izuv,iy,ix)*wt(1,i)*
     +                   wt(2,j)*wt(3,k)
            part%Tf = part%Tf+Text(izuv,iy,ix)*wt(1,i)*
     +                   wt(2,j)*wt(3,k)
            part%qinf = part%qinf+T2ext(izuv,iy,ix)*wt(1,i)*
     +                   wt(2,j)*wt(3,k)
         end do
         end do
         end do


      end subroutine uf_interp

      subroutine particle_stats(ipt,jpt,kpt)
      use particles
      use pars
      use con_data
      implicit none
      integer :: i,ipt,jpt,kpt
      real :: partmass,Volp,rhop

      !Takes in ipt,jpt,kpt as the node to the "bottom left" of the
      !particle
      !(i.e. the node in the negative direction for x,y,z)
      !and computes quantities needed to get particle statistics

      !partcount_t(kpt,jpt,ipt) = partcount_t(kpt,jpt,ipt) + 1.0

      Volp = pi2*2.0/3.0*part%radius**3
      rhop = (part%m_s+Volp*rhow)/Volp
      partmass = rhop*Volp
      partcount_t(kpt,jpt,ipt) = partcount_t(kpt,jpt,ipt) + partmass  !PARTICLE AVERAGES ARE MASS-WEIGHTED

!      !Get su mean, mean-squared of particle velocities at each level
!      upwp_t(kpt,jpt,ipt) = upwp_t(kpt,jpt,ipt) +
!     +                      partmass*part%vp(1)*part%vp(3)
!      do i = 1,3
!      vpsum_t(kpt,jpt,ipt,i) = vpsum_t(kpt,jpt,ipt,i) +
!     +                         part%vp(i)*partmass
!      vpsqrsum_t(kpt,jpt,ipt,i)=vpsqrsum_t(kpt,jpt,ipt,i)+
!     +                          partmass*part%vp(i)**2
!      end do
!
!      Tpsum_t(kpt,jpt,ipt) = Tpsum_t(kpt,jpt,ipt) +
!     +                       partmass*part%Tp
!      Tpsqrsum_t(kpt,jpt,ipt) = Tpsqrsum_t(kpt,jpt,ipt) +
!     +                          partmass*part%Tp**2
!      wpTpsum_t(kpt,jpt,ipt) = wpTpsum_t(kpt,jpt,ipt) +
!     +                          partmass*part%Tp*part%vp(3)
!      Tfsum_t(kpt,jpt,ipt) = Tfsum_t(kpt,jpt,ipt) +
!     +                       partmass*part%Tf
!      qfsum_t(kpt,jpt,ipt) = qfsum_t(kpt,jpt,ipt) +
!     +                       partmass*part%qinf
!      radsum_t(kpt,jpt,ipt) = radsum_t(kpt,jpt,ipt) +
!     +                       partmass*part%radius
!      rad2sum_t(kpt,jpt,ipt) = rad2sum_t(kpt,jpt,ipt) +
!     +                       partmass*part%radius**2
!      qstarsum_t(kpt,jpt,ipt) = qstarsum_t(kpt,jpt,ipt) +
!     +                       partmass*part%qstar

      end subroutine particle_stats

      subroutine particle_stats_U(ipt,jpt,kpt)
      use particles
      use pars
      use con_data
      implicit none
      integer :: i,ipt,jpt,kpt
      real :: partmass,Volp,rhop

      !Takes in ipt,jpt,kpt as the node to the "bottom left" of the
      !particle
      !(i.e. the node in the negative direction for x,y,z)
      !and computes quantities needed to get particle statistics

      !partcount_t(kpt,jpt,ipt) = partcount_t(kpt,jpt,ipt) + 1.0

      Volp = pi2*2.0/3.0*part%radius**3
      rhop = (part%m_s+Volp*rhow)/Volp
      partmass = rhop*Volp
      partcountU_t(kpt,jpt,ipt) = partcountU_t(kpt,jpt,ipt) + partmass  !PARTICLE AVERAGES ARE MASS-WEIGHTED

      !Get su mean, mean-squared of particle velocities at each level
      upwp_t(kpt,jpt,ipt) = upwp_t(kpt,jpt,ipt) +
     +                      partmass*part%vp(1)*part%vp(3)
      do i = 1,3
      vpsum_t(kpt,jpt,ipt,i) = vpsum_t(kpt,jpt,ipt,i) +
     +                         part%vp(i)*partmass
      vpsqrsum_t(kpt,jpt,ipt,i)=vpsqrsum_t(kpt,jpt,ipt,i)+
     +                          partmass*part%vp(i)**2
      end do

      Tpsum_t(kpt,jpt,ipt) = Tpsum_t(kpt,jpt,ipt) +
     +                       partmass*part%Tp
      Tpsqrsum_t(kpt,jpt,ipt) = Tpsqrsum_t(kpt,jpt,ipt) +
     +                          partmass*part%Tp**2
      wpTpsum_t(kpt,jpt,ipt) = wpTpsum_t(kpt,jpt,ipt) +
     +                          partmass*part%Tp*part%vp(3)
      Tfsum_t(kpt,jpt,ipt) = Tfsum_t(kpt,jpt,ipt) +
     +                       partmass*part%Tf
      qfsum_t(kpt,jpt,ipt) = qfsum_t(kpt,jpt,ipt) +
     +                       partmass*part%qinf
      radsum_t(kpt,jpt,ipt) = radsum_t(kpt,jpt,ipt) +
     +                       partmass*part%radius
      rad2sum_t(kpt,jpt,ipt) = rad2sum_t(kpt,jpt,ipt) +
     +                       partmass*part%radius**2
      qstarsum_t(kpt,jpt,ipt) = qstarsum_t(kpt,jpt,ipt) +
     +                       partmass*part%qstar

      end subroutine particle_stats_U

      subroutine update_particle_flux

      use particles
      use pars
      use con_data
      use con_stats
      use tank_pars
      implicit none
      include 'mpif.h'
      real :: rhop,taup_i,partmass
      real :: xv,yv,zv
      integer :: i,j,k,ncount,ipt,jpt,kpt,kwpt,kUpt
      integer :: ix,iy,iz
      real :: lref,Cvv,Cva,Cpv
      real :: zu_1,zu_nnz,ztocheck,ztocheck_i
      integer :: iz_beg,iz_end


      part => first_particle
      do while (associated(part))

      ipt = floor(part%xp(1)/dx) + 1
      jpt = floor(part%xp(2)/dy) + 1


      if (i_zeta .eq. 1) then  !Check fluxes according to zeta levels
          ztocheck = part%xip(3)
          ztocheck_i = part%xip_i(3)
      elseif (i_zeta .eq. 0) then
          ztocheck = part%xp(3)  !Check fluxes according to z-levels
          ztocheck_i = part%xp_i(3)
      end if

      kUpt = floor((ztocheck-minwavez)/dz_U) + 1

      !Also compute particle flux statistics, which is slightly more
      !complicated:
      zu_1 = minwavez + 0.5*dz_U
      zu_nnz = zl - 0.5*dz_U

      if (ztocheck .gt. minwavez) then
      !Compute particle flux statistics based on the move:
      if (ztocheck_i .lt. ztocheck) then  !Particle moved up

        iz_beg = floor((ztocheck_i-zu_1)/dz_U) + 2
        iz_end = floor((ztocheck-zu_1)/dz_U) + 1

        !Calling this before particle_bcs_nonperiodic, so must account for when particles are below
        if (iz_beg .lt. 1) then
           iz_beg = 1
        end if

         do iz=iz_beg,iz_end

            if (iz .lt. 1) then
            write(*,*) 'up iz<1:',iz,ztocheck,ztocheck_i
            end if
            if (iz .gt. nnz) then
            write(*,*) 'up iz>nnz:',iz,ztocheck,ztocheck_i
            end if

            !partflux_t(iz,jpt,ipt) = partflux_t(iz,jpt,ipt) + 1.0
            !partfluxup_t(iz,jpt,ipt) = partfluxup_t(iz,jpt,ipt) + 1.0
         end do

      elseif (ztocheck_i .gt. ztocheck) then  !Particle moved down

        iz_beg = floor((ztocheck-zu_1)/dz_U) + 2
        iz_end = floor((ztocheck_i-zu_1)/dz_U) + 1

        if (iz_beg .lt. 1) then
           iz_beg = 1
        end if

         do iz=iz_beg,iz_end

            if (iz .lt. 1) then
            write(*,*) 'dn iz<1:',iz,ztocheck,ztocheck_i
            end if
            if (iz .gt. nnz) then
            write(*,*) 'dn iz>nnz:',iz,ztocheck,ztocheck_i
            end if

            !partflux_t(iz,jpt,ipt) = partflux_t(iz,jpt,ipt) - 1.0
            !partfluxdn_t(iz,jpt,ipt) = partfluxdn_t(iz,jpt,ipt) - 1.0
         end do

      end if
      end if


      part => part%next
      end do


      end subroutine update_particle_flux

      subroutine update_particle_stats

      use particles
      use pars
      use con_data
      use con_stats
      use tank_pars
      implicit none
      include 'mpif.h'
      real :: rhop,taup_i,partmass
      real :: xv,yv,zv
      integer :: i,j,k,ncount,ipt,jpt,kpt,kwpt,kUpt
      integer :: ix,iy,iz
      real :: lref,Cvv,Cva,Cpv
      real :: zu_1,zu_nnz,ztocheck
      integer :: iz_beg,iz_end


      part => first_particle
      do while (associated(part))

      ipt = floor(part%xp(1)/dx) + 1
      jpt = floor(part%xp(2)/dy) + 1
      kUpt = floor((part%xip(3)-0.0)/dz_U) + 1

      if (i_zeta .eq. 1) then  !particle stats according to zeta levels
          ztocheck = part%xip(3)
      elseif (i_zeta .eq. 0) then
          ztocheck = part%xp(3)  !particle stats according to z-levels
      end if

      kUpt = floor((ztocheck-minwavez)/dz_U) + 1

      call particle_stats_U(ipt,jpt,kUpt)



      !First, as done in uf_interp, must find the "leftmost" node
      !of volume where particle belongs:
      !(must repeat since now particle locations have been updated)

      do iz=0,nnz+1
         if (zz(iz) .GT. part%xip(3)) go to 1299
      end do
1299  continue
      kpt = iz-1

      !Do again for w:
      do iz = 0,nnz+1
         if (z(iz) .GT. part%xip(3)) go to 1120
      end do
1120  continue

      kwpt = iz


!---------
      !As an aside, use ipt,jpt,kwpt to update the particle numbers at
      !each cell:

      !call inputs are the nodes associated with the volume that each
      !part. lies in
      call particle_stats(ipt,jpt,kwpt)


!---------

      part => part%next
      end do


      end subroutine update_particle_stats
