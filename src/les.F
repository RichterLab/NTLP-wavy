      program les_mpi_2d
c
        use pars
        use fields
        use particles
        use con_data
        use con_stats
        use netcdf_io
        include 'mpif.h'
c
c ------------- definition of internal flags
c
c       iDNS    = 0; call the subgrid computation of vis_m and vis_s
c               = 1; call the molecular viscosity and diffusivity
c
c       igrdr   =  3; data comes from restart file
c               =  2; data comes from initialization (random)
c               =  1; data comes from coarser grid (or otherwise)
c
c       ibcu    = 1 ; upper boundary condition set by radiation bc
c               = 0 ; fixed value = 0.
c               = -1; wind tunnel setup, gradient = 0 bcs
c               = 2; Sullivan et al. Couette flow, sets u = Uo and w = v = 0 (NO SCALAR CHANGE YET)
c
c       ibcl    = 0 ; lower boundary condition set by similarity theory (sr. setup)
c               = -1; value defined by coarser mesh for all variables
c
c       ifix_dt = 0 ; variable time step with fixed cfl number in setcon
c               = 1 ; fixed time step set in sr. get_dt
c
c       ifree   = 0 ; use spatially averaged surface conditions for MO (call lower)
c               = 1 ; use point-by-point conditions for MO free convection (call lower_free)
c               = 2 ; use point-by-point conditions for MO with curvy boundaries (call lower_curvy)
c               = 3 ; DNS no-slip conditions on the bottom boundary (call lower_curvy_dns)
c
c       iflat   = 0 ; flat lower bottom, no pressure iteration
c               = 1 ; curvy lower boundary terrain following grid
c                     iterate for pressure field
c               = 2 ; time dependent lower boundary shape
c
c       lpbc    = 0 ; dpdz = 0 lower boundary condition for (u,v,w,uf,vf) updates
c               = 1 ; linear extrapolation to get ghost point for p(iz=0)
c               = 2 ; solve for p(iz=0) from wf equation at boundary
c
c       imesh   = 0 ; stationary bottom
c               = 1 ; moving lower bottom
c
c       move_step  < 0  ; never move the bottom
c                    = it ; step where the moving grid starts
c
c       i_pm    = 1 ; use Peirson-Moskowitz wave spectrum
c              /= 1 ; use Donelan spectrum with wave age and wind dependence
c
c       i3d_surf   = 0  ; call simple l_bndy surface
c                  = 1  ; call moving 3d surface (water waves)
c                  = 2  ; call moving 3d surface with time varying amplitude
c                  = 3  ; call moving 3d surface with time varying amplitude (tank)
c
c       ihst    = nn; frequency at which global variables are output in history file
c               < 0 ; no history files
c
c       iviz    = nn; frequency at which global variables are output in viz files
c               < 0 ; no viz files
c
c       it_his  = time step where history files start, incremented by itape
c
c       it_viz  = time step where viz files start, incremented by itape
c
c       ismlt   = 1 ; use businger formulas in MO
c                 0 ; use large and everyone elses formulas in MO
c
c       iupwnd  = 0;  use skew symmetric formulas for all derivatives
c                     in scalar equations
c               = 1;  use hybrid upwind scheme for all derivatives
c                     in scalar equations
c
c       iskew   = 0;  use skew symmetric form for advection in (u,v,w) equations
c               = 1;  use flux form for advection in momentum equations
c
c       ivis0   = 0; old eddy viscosity model
c               = 1; new eddy viscosity model
c
c       iz_space = 1; non-uniform spacing with first grid point = zw1
c                = 2; non-uniform grid with adjusted first point off
c                     the surface to give constant stretching
c                = 3; non-uniform exponential curvy grid, good for
c                     moving boundaries
c                = 4; non-uniform exponential curvy grid with constant
c                     grid space at boundary and first grid point, good for
c                     moving boundaries
c                = 5; algebraic mapping with control over decay with height
c                     good for wavy boundaries. see a_decay for grid decay
c                = 0; uniform spacing
c
c       k_stab   = an iz (or k) index where stability corrected length scale is turned on
c                = -1 correction is always on
c                > nnz correction is never on
c                = value correction is on at iz = value
c
c       new_vis = step; the iteration step for which the new model
c                       is turned on when ivis0=1
c               < 0; new model is on at all steps for ivis0=1
c
c       nscl  .ge. 1   number of scalars to be followed set in parameter statements
c                      change entries in sr. init, and sr. suft for surface bc's
c
c -------------------------------------------------------------------------------
c
c ---------- initialize MPI, get myid, numprocs,
c            test if on root process
c
        call mpi_init(ierr)
        call mpi_comm_rank(mpi_comm_world,myid,ierr)
        call mpi_comm_size(mpi_comm_world,numprocs,ierr)
        if (myid==0.and.verb>0) write(*,*) 'Finished mpi initialization'
c
        i_root = 0
        l_root = .false.
        if(myid == i_root) l_root = .true.
c
        l_debug = .false.
        if(idebug == 1) l_debug = .true.
c
        ts_mpi = mpi_wtime()

      !Read the input file for all parameters
        call read_input_file

c
c -------- set number of x-y slab cpus
c
c
        case_inp = 'wav'
c
        call get_units
        call gridd
        call setcon
        call set_paths

        call change_RH_bcs_to_q
c
c -------------- scratch run
c
        ! if this is a fresh run
        if (iti==0)  then
          igrdr = 2
          case = case_inp
          call init
          call setup(it)
          call mesh
          if (inetcdf == 1) then
            call netcdf_init
          end if

          if (ispray==1) then
            call particle_setup
            call particle_init
          end if

c
c ---------- choose routine for getting initial guess
c
          if (myid==0.and.verb>0) write(*,*)
     +        'Choosing routine for getting initial guess'
          if(iocean == 1) then
            call randoc
          else
            if (ifields == 1) then
              call get_fields
            else
              call random
               !call random_f
               !call random_wt
            end if
          end if
          call dns_vis
          call get_max

        ! if this is a restart run
        else ! (iti\=0)
          igrdr = 3
          call restart

          call get_max
          call setup(it)
          call mesh

          if (ispray==1) then
            call particle_setup
            call read_part_res
          end if

          if (inetcdf .eq. 1) then
             call netcdf_res
          endif

        end if



        if (myid == 0) then
          idum = 1
          write(*,*) 'ran1(1)',  ran1(idum)
          ! write(*,*) 'ran1(1)',  ran1(1)
          ! write(*,*) 'ran1(-1)', ran1(-1)
          ! write(*,*) 'ran1(-1)', ran1(-1)
          ! write(*,*) 'ran1(1)',  ran1(1)
          ! write(*,*) 'ran1(1)',  ran1(1)
          ! write(*,*) 'ran1(-1)', ran1(-1)
          ! write(*,*) 'ran1(-1)', ran1(-1)
          !
          ! write(*,*) 'ran2(1)',  ran2(1)
          ! write(*,*) 'ran2(1)',  ran2(1)
          ! write(*,*) 'ran2(-1)', ran2(-1)
          ! write(*,*) 'ran2(-1)', ran2(-1)
          ! write(*,*) 'ran2(1)',  ran2(1)
          ! write(*,*) 'ran2(1)',  ran2(1)
          ! write(*,*) 'ran2(-1)', ran2(-1)
          ! write(*,*) 'ran2(-1)', ran2(-1)


        end if
c
c --------------- time loop ------------
c
        tzero = time
        call get_dt(it)
 9000   continue
        call set_sav(it,iti)

        if (myid==0) then
          write(*,*) 'Starting time loop'
          write(*,*) 'it,time = ',it,time
        end if
        if (ispray==1) then
          if (it == 1) numpart = 0
          part => first_particle
          do while (associated(part))
            if (part%pidx == 1 .and. part%procidx == 0) then
              write(*,'(a4,4e15.6)') 'xp1:',time,part%xp(1:3)
              write(*,'(a5,4e15.6)') 'xip1:',time,part%xip(1:3)
              write(*,'(a4,4e15.6)') 'vp1:',time,part%vp(1:3)
              write(*,'(a4,4e15.6)') 'uf1:',time,part%uf(1:3)
              write(*,'(a4,2e15.6)') 'Tp1:',time,part%Tp
              write(*,'(a4,2e15.6)') 'Tf1:',time,part%Tf
              write(*,'(a4,2e15.6)') 'qf1:',time,part%qinf
              write(*,'(a5,2e15.6)') 'rad1:',time,part%radius
            end if
            part => part%next
            if (it == 1) numpart = numpart + 1
          end do
          if (myid==0) write(*,*) 'time,tnumpart:',time,tnumpart
        end if

c
c --------- specially designed surface cooling routine
c           for gabls run
c
c     call forcing
c
        if(it >= new_vis .and. ivis0 == 1) then
          ivis = 1
        else
          ivis = 0
        end if
c
c ---------------- 3 stage runge-kutta time stepping
c
        do  8999 istage=1,3
c
          if (myid==0.and.verb>0) write(*,*) 'rk3 stage:', istage
          dtzeta   = dt*zetas(istage)
          dtgama   = dt*gama(istage)
          dtalp    = dtzeta + dtgama
          t_stage  = time + dt*etas(istage)
          t_stage2 = time + dt*etas(istage+1)
          if(istage <= 2) then
            t_stage3 = time + dt*etas(istage+2)
          else
            t_stage3 = time + dt + dt1*etas(istage-1)
          end if
c
c ---------- compute derivatives of (u,v,w)
c
          call exchange
          call get_derv
c
c --------- check for moving grid
c
          if(imesh == 1) then
            call new_mesh(it,istage)
            call grid_speed(it,istage)
          end if

c
c --------- new eddy viscosity, and bcs
c
          if (myid==0.and.verb>0) write(*,*) 'Handling lower bc'
          if(iss == 0 .and. ifree == 0) then
            call lower(it)
          else if(ifree == 1) then
            call lower_free(it)
          else if(ifree == 2) then
            call lower_curvy(it,istage)
          else if(ifree == 3) then
            call lower_curvy_dns(it,istage)
          end if
c
          if (myid==0.and.verb>0) write(*,*) 'Handling upper bc'
          if(ise == numprocs-1) then
            if (iDNS == 1) then
              call upper_dns
            else
              call upper
            end if
          end if
          call bcast_pbc
          call bcast_surf
          call get_means(istage)
          if(ivis == 1) then
            call iso(it)
            call surfvis(it)
          end if
          if(istage == 1)then
            call xy_stats
            call tke_budget
            call pbltop(itop)
          end if
c
c ------------ save velocity aand pressure fields
c
          if(msave .and. istage == 1) then
            if(myid==0.and.verb>0)write(*,*)'Saving fields and parts'
            call save_v(it)
            call save_p
            if (ispray==1) call save_particles
          end if
          if(msave_v .and. istage == 1) then
            call save_bndy(it)
!         call save_viz(it)
          end if

          if (msave .and. istage == 1) then
            call write_histograms(it-itape+1)
            call clean_histograms
          end if

c
c --------- get rhs/J for all equations
c
          call comp1(istage,it)
c
          if(istage == 1) then
            if(msave .and. l_root) call save_c(it)
          end if
c
c --------- solve for pressure
c
          call psolver(istage)
c     call chk_div
c
c --------- add pressure gradient and dealias
c
          call comp2
c
          if(micut) then
            call dealias
            call dealias_grid_speed
          end if
c     call chk_div

c
c -------- update particles
c
          if (myid==0.and.verb>0) write(*,*) 'Updating particles'
          if (ispray==1) then
         !t_s = mpi_wtime()
            call particle_update_rk3(it,istage)
         !t_f = mpi_wtime()
         !call mpi_barrier(mpi_comm_world,ierr)
         !if (myid==5) write(*,*) 'time part: ',t_f-t_s
          end if

          if(mnout .and. istage == 1)  then
            if(l_debug) call print(nprt,it,izs,ize)
            if(l_root) call print(6,it,1,nnz)
          end if
          if(l_root) then
            if (inetcdf == 1) then
              if(mhis .and. istage == 1) then
                call write_his_netcdf
                if(mtape) call close_his_netcdf
              end if
            else
              if(mhis .and. istage == 1) then
                call write_his(itop)
                if(mtape) call close_his
              end if
            end if !(inetcdf == 1)
          end if !(l_root)
 8999   continue
        call get_max
        call get_dt(it)

        if (mhis .and. itrajout) then
           call particle_write_traj(it)
        end if
c
        if (it>=itmax) go to 99000
        go to 9000
c
99000   continue
        te_mpi = mpi_wtime()
        write(6,9997) (te_mpi - ts_mpi)
 9997   format(' Job Execution Time = ',e15.6)
c
 9998   continue
        call mpi_finalize(ierr)
c
        stop
      end
