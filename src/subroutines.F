
      subroutine get_max
c
c --------- routine computes max velocities as sweep through
c           the flux velocity field
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
c
      real u_send(6), u_recv(6)
c
      dx_i = 1.0/dx
      dy_i = 1.0/dy
c
      u_temp   = 0.0
      v_temp   = 0.0
      w_temp   = 0.0
      e_temp   = 0.0
      vis_temp = 0.0
      do iz=izs,ize
c
        u_xy = 0.0
        v_xy = 0.0
        w_xy = 0.0
        e_xy = 0.0
        g_xy = gm(1,iys,iz,1)
        do iy=iys,iye
        do ix=1,nnx
c
c          u_xy = amax1(u_xy,abs(u(ix,iy,iz)+stokes(iz)))
c          v_xy = amax1(v_xy,abs(v(ix,iy,iz)))
c          w_xy = amax1(w_xy,abs(w(ix,iy,iz)))
           u_xy = amax1(u_xy,abs(uf(ix,iy,iz)))
           v_xy = amax1(v_xy,abs(vf(ix,iy,iz)))
c          w_xy = amax1(w_xy,abs(wf(ix,iy,iz) - gm(ix,iy,iz,5)))
           w_xy = amax1(w_xy,abs(wf(ix,iy,iz)),abs(gm(ix,iy,iz,5)))
           e_xy = amax1(e_xy,e(ix,iy,iz))
           g_xy = amax1(g_xy,gm(ix,iy,iz,1))
        enddo
        enddo
        u_xy   = u_xy*dx_i
        v_xy   = v_xy*dy_i
        wsav   = w_xy
        w_xy   = w_xy/abs(dzw(iz))
c
c ------- get max allowable diffusivity for scalars
c
        vis_xy = dfac(iz)*3.0*ck*dsl_z(iz)*(g_xy**pow3)*
     +           sqrt(e_xy)/amin1(dx,dy,dzw(iz))**2
c
        u_temp   = amax1(u_xy,u_temp)
        v_temp   = amax1(v_xy,v_temp)
        w_temp   = amax1(w_xy,w_temp)
        e_temp   = amax1(e_xy,e_temp)
        vis_temp = amax1(vis_xy,vis_temp)
c
c       if(iz .le. 15) then
c         write(6,6000) iz, wmax
c6000     format(' in get_dt iz = ',i3,' wmax = ',e15.6)
c       endif
c
      enddo
      u_send(1) = u_temp
      u_send(2) = v_temp
      u_send(3) = w_temp
      u_send(4) = wsav
      u_send(5) = e_temp
      u_send(6) = vis_temp
c
      call mpi_allreduce(u_send,u_recv,6,mpi_real8,
     +     mpi_max,mpi_comm_world,ierror)
c
      umax   = u_recv(1)
      vmax   = u_recv(2)
      wmax   = u_recv(3)
      wabs   = u_recv(4)
      emax   = u_recv(5)
      vismax = u_recv(6)
c
      return
      end
      subroutine get_dt(it)
c
c ---------- routine computes max time step for given cfl number
c            from max's found previously
c
      use pars
      use con_data
      use con_stats
c
      ucfl = umax
      vcfl = vmax
      wcfl = wmax
      vel_max = amax1(ucfl, vcfl, wcfl)
      if(vel_max .le. 0.0) then
          write(6,6000) ucfl, vcfl, wcfl, vel_max
 6000     format('6000, sr. get_dt bad news, umax = ',e15.6,/,
     +           ' vmax = ',e15.6,' wmax = ',e15.6,/,
     +           ' vel_max = ',e15.5,/,
     +           ' infinite time step !!!')
          stop
      endif
c
c ---------------- choose fixed or variable time step
c
      if(ifix_dt .ne. 0) then
c
c ------------- if used, change to fit your problem
c
         dt_new = 0.000370099
      else
c
c ------------------- new estimate of best time step
c                     from cfl constraint
c
      dt_new  = cfl/vel_max
      dt_new = amin1(dt_new, 5.0)
c     dt_new = amin1(dt_new, 1.0)
c     dt_new = amin1(dt_new, 20.0)
c     dt_new = amin1(dt_new, 10.0)
      endif
c
c ---------------- compare against viscous stability limit
c
      if(vismax*dt_new .gt. 0.5) then
         dt_cfl = dt_new
         dt_new = 0.5/vismax
         if(l_root) then
            write(6,6200) dt_new, dt_cfl, vismax
 6200       format(' 6200 get_dt: cfl time step too large',/,
     +      '   viscous time step = ',e15.6,
     +      ' cfl time step = ',e15.6,' vismax = ',e15.6)
         endif
      endif
c
c -------- for safety if restart set timestep = saved timestep in
c          restart file
c
CHECK ... for moving grids?
c
      if(it .eq. iti .and. iti .ne. 0) then
        dt_new = dt1
        dt1    = dt
        write(nprt,4001) it, iti, dt_new, dt1, dt
 4001   format(' 4001: it = ',i5,' iti = ',i5,
     +         ' dt_new = ',e15.6,' dt1 = ',e15.6,
     +         ' dt = ',e15.6)
      endif
c
      return
      end
      subroutine lterp(n,zary,zpt,i,ip1,ratio)
c
c ---- linear interpolation for zpt in zary, where zary is
c      monotonic increasing or decreasing function
c
      dimension zary(*)
      nm1 = n-1
      if(n.le.1) then
          i = 1
          ip1 = 1
          ratio = 0.0
          go to 999
      endif
      if(zary(1) .lt. zary(2)) go to 1
                               go to 101
    1 continue
c
c **** monotonic increasing array
c
        if(zpt .lt. zary(1)) then
          i = 1
          ip1 = 1
          ratio = 0.0
          go to 999
        else if(zpt .gt. zary(n)) then
          i = n
          ip1 = n
          ratio = 1.0
          go to 999
        endif
        do j=1,nm1
              if(zpt .ge. zary(j) .and.
     $           zpt .le. zary(j+1)) then
                 i = j
                 ip1 = j+1
                 ratio = (zpt - zary(i))/(zary(ip1) - zary(i))
                 go to 999
              endif
        enddo
c
c **** decreasing array
c
  101 continue
        if(zpt .gt. zary(1)) then
          i = 1
          ip1 = 1
          ratio = 0.0
          go to 999
        else if(zpt .lt. zary(n)) then
          i = n
          ip1 = n
          ratio = 1.0
          go to 999
        endif
        do j=1,nm1
              if(zpt .le. zary(j) .and.
     $           zpt .ge. zary(j+1)) then
                 i = j
                 ip1 = j+1
                 ratio = (zpt - zary(i))/(zary(ip1) - zary(I))
                 go to 999
              endif
        enddo
  999 continue
      return
      end
      subroutine setcon
c
      use pars
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
      external get_zi
c
c ----------------- get machine type, can read datadir also
c
      open(99,file='mach.file')
      read(99,9000) imach
 9000 format(i1)
c
      close(99)
c
c ----------- initialize fft
c
      call rffti(nnx,trigx(1,1))
      call rffti(nny,trigx(1,2))
      call cffti(nny,trigc(1))
      call cffti(nnx,trigcx(1))
c
c ----------- start step for history and viz files
c
      it_his_nxt = it_his
      it_viz_nxt = it_viz
c
c ---------------- set min value of e
c
      if(iocean .eq. 1) then
c
         smal_e = 0.0
         smal_e = 1.0e-12
c        smal_e = 6.0e-03
      else
         smal_e = 1.0e-09
c        smal_e = 0.0
      endif
c
c
c ---------------------- set constants in eddy viscosity model
c
      ck       = 0.1
      ceps     = 0.93
      csmag    = sqrt(ck*sqrt(ck/ceps))
      stab_c   = 0.76
      pow3     = 1.0/3.0
c
c ----------------- set stability constant
c
      stabmin = 1.0e-12
c
c ---------------- minimum dsl length constant
c
      almin_c = 0.0001
c
c -------------------- set viscosity model parameters
c                      to usual values, overridden by iso
c                      dfac is in constant file for calculation
c                      of max eddy viscosity
c
      viscon      = 0.0
      xksurf      = 0.0
      nmatch      = -1
      myid_newvis = 0
      do iz=1,nnz
         dfac(iz) = 1.0
      enddo
c
c -------------- initialize grid restart flag
c
      igrdr = 1
c
c ---------- set constants for pressure iteration
c
      p_surf = 0.0
c
      if(iflat .eq. 0) then
         ip_max    = 2
         ptol      = 0.0
         accel_max = 1.0
         accel_min = 1.0
         do i=1,ip_accel
            accel(i) = 1.0
         enddo
      else
c
c ------- with effectively flat just use 2 iterations
c
         ip_max  = 2
         ptol    = 1.0e-07
c        ptol    = 1.0e-20
c
         if(i3d_surf == 3) then   !  strategy for breaking waves in tank
            ptol      = 1.0e-10
            ip_max    = 50
            ip_max    = 35
            accel_max = 1.0
            accel_min = 0.1
            ip_ramp   = 30
c
c           accel_max = 1.0
c           accel_min = 1.0
c           ip_max    = 15
c
c
c ---------- for initialization with flat
c
c           ip_max    = 2
c           ptol      = 1.0e-07
c           accel_max = 1.0
c           accel_min = 1.0
c           ip_ramp   = 30
         else                     !  strategy for pbl with surface waves or hills
            ptol      = 1.0e-10
            ip_max    = 30
c           ip_max    = 2
c
c           ip_max    = 1 ! with wave amplitude = 0
c
            accel_max = 1.0
            accel_min = 1.0
            ip_ramp   = 10
         endif
         d_accel   = (accel_max - accel_min)/float(ip_ramp - 1)
         do i=1,ip_accel
            accel(i) = accel_max
         enddo
         do i=1,ip_ramp
            accel(i) = float(i-1)*d_accel + accel_min
         enddo
         write(nprt,1017) (i,accel(i),i=1,ip_accel)
 1017    format(' acceleration parameters',/,
     +          ' i ',5x,' accel',/,(i5,e15.6))
      endif
c
c -------------- set grid generation factor for type of simulation
c                in routines exp_grid, exp_grid2
c
         if(i3d_surf == 0) then
            a_exp_grd   = 0.020
            a_exp_grd2  = 0.00010
            a_decay     = 3.0
         elseif(i3d_surf == 1) then
            a_exp_grd   = 0.020
            a_exp_grd2  = 0.00010
            a_decay     = 3.0
         elseif(i3d_surf == 2) then
            a_exp_grd   = 0.020
            a_exp_grd2  = 0.00010
            a_decay     = 3.0
         elseif(i3d_surf == 3) then
c           a_exp_grd   = 50.00
            a_exp_grd   = 11.65
            a_exp_grd2  = -1000.0 ! fix for tank simulations
            a_decay     = 3.0
         endif
c
c -------------- create mpi operation to find max and location
c                using local gradient method
c
      call mpi_op_create(get_zi,.true.,ziloc,ierror)
c
c ---------- define coefficients for 3-order runge-kutta
c            time stepping scheme, borrowed from Spalart,
c            Moser and Rogers, J. Comp. Physics 3/21/90
c            Note this is a simplier version since all terms
c            are lumped in the non-linear terms.
c            cfl number is for an entire runge-kutta step
c            in this case three stages. cfl = max(u)*dt/dx
c
c            see jcp article by berkooz for time steps
c
      zetas(1) = 0.0
      zetas(2) = -17.0/60.0
      zetas(3) = -5.0/12.0
      gama(1)  = 8.0/15.0
      gama(2)  = 5.0/12.0
      gama(3)  = 3.0/4.0
c
      etas(1)  = 0.0
      etas(2)  = 8.0/15.0
      etas(3)  = 2.0/3.0
c
c ----------- a full step, at the new time
c
      etas(4)  = 1.0
c
c     cfl = 0.63
c
c     cfl = 0.50
c     cfl = 0.30
c     cfl = 0.10
c
c --------- set stationary flag on a restart
c           see sr. set_sav for tests
c
      imesh = 0
c
      return
c
      end
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
      subroutine setup(it)
c
      use particles
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      it = iti
c
c ------------ turn on new sgs model at a particular step
c
      if(it .ge. new_vis .and. ivis0 .eq. 1) then
          ivis = 1
      else
          ivis = 0
      endif
c
      if(igrdr . eq. 3) then
         if(l_root) then
            write(6,6)iti,utau,tsfcc(1) ,qstar(1)
            write(6,510)
            write(6,520) wtsfc(1),utau,amonin,dtdzf(1),zody,zo
     +            ,cdbtm,ugcont
            call print(6,it,1,nnz)
         endif
          if(l_debug) then
            write(nprt,6)iti,utau,tsfcc(1) ,qstar(1)
            write(nprt,510)
            write(nprt,520) wtsfc(1),utau,amonin,dtdzf(1),zody,zo
     +            ,cdbtm,ugcont
            call print(nprt,it,izs,ize)
          endif
      endif
c     if(ifilt.eq.1)call filter
      if(l_root) then
         write(6,1) nnx,nny,nnz,ismlt,ifilt,iti,itmax,
     +             iupwnd,ibuoy,itcut,
     +             dt,zo,tsfcc(1),isubs,ibrcl,
     +             method, iocean, ivis
      endif
      if(l_debug) then
         write(nprt,1) nnx,nny,nnz,ismlt,ifilt,iti,itmax,
     +             iupwnd,ibuoy,itcut,
     +             dt,zo,tsfcc(1),isubs,ibrcl,
     +             method, iocean, ivis
      endif
c
c -------------- boundary condition flags
c
      !ibcu = 2
      ibcu = -1
      ibcl = 0
c
c -------------------- wavenumbers, introduce a normalized
c                      set of wavenumbers to eliminate computation
c                      in derivatives , xderiv, yderiv
c
      do i=1,nnx
         xkn(i) = float(i-1)*pi2/xl
         if(i.gt.ncx)xkn(i) = -float(nnx-i+1)*pi2/xl
      enddo
      fn = 1.0/float(nnx)
      do i=1,nnx
         xk(i)  = xkn(i)*fn
         xkk(i) = xkn(i)*xkn(i)*fn
      enddo
      do i=1,nny
         ykn(i) = float(i-1)*pi2/yl
         if(i.gt.ncy)ykn(i) = -float(nny-i+1)*pi2/yl
      enddo
      fn = 1.0/float(nny)
      do i=1,nny
         yk(i)  = ykn(i)*fn
         ykk(i) = ykn(i)*ykn(i)*fn
      enddo
      ii = -1
      do i=1,ncx
         ii = ii + 2
         temp = xkn(i)**2
         do j=1,nny
            temp1       = temp + ykn(j)**2
            xks(ii,j)   = temp1
            xks(ii+1,j) = temp1
         enddo
      enddo
      xnn = abs(batag*dtdzf(1))
c
c ----------- choose correct sign so gravity waves
c             propagate out of the domain
c
      sgn = 1.0
      if(iocean .eq. 1) sgn = -1.0
      if(ibcu.eq.1) then
         do iy=1,nny
         do ix=1,nnxp2
            if(xks(ix,iy) .le. 0.) then
              wavexy(ix,iy) = 0.0
            else
              wavexy(ix,iy) = sgn*sqrt(xnn/xks(ix,iy))
            endif
         enddo
         enddo
      endif
c
c -------------------- set length scale for SGS model
c
      if(iz_space .eq. 0) then
c
c ------------- uniform vertical spacing
c
      dx32 = dx*3./2.
      dy32 = dy*3./2.
      dsl  = (abs(dx32*dy32*dzw(1)))**(1./3.)
      dslg = dsl
c
c -------
c
      do iz=0,nnzp1
         dsl_z(iz) = dslg
      enddo
      if(l_root)  write(6,2000) dsl
      if(l_debug) write(nprt,2000) dsl
      if(l_root) write(6, 4033) (iz,dsl_z(iz),iz=0,nnz+1)
c
c ------------- variable vertical spacing
c
      else
c
c ----------- just estimate dsl for average spacing
c
         dx32 = dx*3./2.
         dy32 = dy*3./2.
c
         dsl_max = (abs(dx32*dy32*dzw(0)))**pow3
         do iz=0,nnzp1
            dsl_z(iz) = (abs(dx32*dy32*dzw(iz)))**pow3
            if(dsl_z(iz) .gt. dsl_max) dsl_max = dsl_z(iz)
         enddo
c        do iz=0,nnzp1
c           dsl_z(iz) = dsl_max
c        enddo
         dsl  = dsl_max
         dslg = dsl
         if(l_root) write(6, 4033) (iz,dsl_z(iz),iz=0,nnz+1)
      endif
c
c -------- SPECIAL for grid refinement study
c
c     nnx_c = 128
c     nny_c = 128
c     nnz_c = 128
c     dx32  = float(nnx)/float(nnx_c)*dx*1.5
c     dy32  = float(nny)/float(nny_c)*dy*1.5
c     dzw_c = float(nnz)/float(nnz_c)*dzw(1)
c
c     dsl  = (abs(dx32*dy32*dzw_c))**(1./3.)
c     dslg = dsl
c     if(l_root)  write(6,2100) dsl
c     if(l_debug) write(nprt,2100) dsl
c
c --------------------- create dsl array for easy indexing in comp1
c
c     do iz=0,nnzp1
c        dsl_z(iz) = dslg
c     enddo
      if(l_root) write(6, 4033) (iz,dsl_z(iz),iz=0,nnz+1)
c
      gridr = 1.0
      sml_eg = smal_e*gridr
c
c -------------------- set viscosity model parameters
c                      now in setcon
c
c     if(ivis .ne. 1) then
c       viscon = 0.0
c       xksurf = 0.0
c       nmatch = -1
c       myid_newvis = 0
c       do iz=1,nnz
c          dfac(iz) = 1.0
c       enddo
c     endif
c ------------------- set stokes velocity for atmos/oceanic flow
c
      call stokesv
c
c --------- special for moving 3d surface
c
      if(i3d_surf .eq. 1) then
         call surf_setup
      elseif(i3d_surf .eq. 2) then
         call surf_setup_gabls
      elseif(i3d_surf .eq. 3) then
         call surf_setup_tank
      endif
c
c --------- can add a time factor so as to skip into any part of
c           the specified geostrophic arrays. time factor in seconds
c
      t_factor = 7200.0
c
c ---------- for print out to get more digits
c
c     t_ref = t00
c
c -------------- for wind tunnel simulations
c                account for box height
c
      t_ref      = 0.0
c
       dpdx_con_z = rhoa*ustar_tank**2/zl
      !dpdx_con_z = dpdx_con/zl
      !dpdx_con_z = 0.0
c
c -------------------- specify cooling rate and initial
c                      temperature even for restarts
c                      from runs with DM
c
      c_rate   = 0.375/3600.0
      t_surf_i = 287.99999
c
c -------------------- do not look for zi below zi_min
c
      zi_min = 30.0
      if(iocean .eq. 1) zi_min = -5.0
      iz_min = 1
      do iz=1,nnz-1
         if(zz(iz) .lt. zi_min .and.
     +      zz(iz+1) .ge. zi_min) iz_min = iz
      enddo
      if(l_root) then
         write(6,9000) zi_min, iz_min
      endif
c
c -------- for curvy bottom set the distribution
c          of specified heat flux. customize for
c          particular use
c
      call set_lower_bc
c
 9998 continue
      return
c --------------------------- format statements
    6 format(///,' DATA FROM RESTART FILE AT STEP =',I5,
     +       ' U_* = ',e15.6,' TS = ',e15.6,' Q_* = ',e15.6,///)
  510 format(' RESTART ***** CASE WITH : ******',/)
  520 format(' WT = ',e12.4,', U_* = ',e12.4,', L = ',e12.4,
     +       ', DTDZ FREE = ',e12.4,', ZODY = ',e12.4,/,10x,
     +       '  ZO(BTM) = ',e12.4,', CDBTM = ',e12.4,
     +       ', UG = ',e12.4)
    1 format(10x,' NNX = ',i5,',  NNY = ',i5,
     + ',  NNZ = ',i5,/,10x,' SFC SMLT = ',i1,
     + ',  FILTER = ',i1,
     + ',  ITI = ',i6,',  ITMAX = ',i6,/,10x,
     + ' IUPWIND = ',i1,',  BUYNCY = ',i1,
     + ',  ITCUT = ',i1,/,10x,
     + ' DT = ',e15.6,',  ZO = ',e15.6,',  TS = ',e15.6,
     + ',  SUBSD = ',i1,/,
     + 10x,' BRCLICITY = ',i1,',  METHOD = ',i1,',  IOCEAN = ',i1,
     + ',  IVIS = ',i1)
 4033 format(' Eddy viscosity length scale for variable z meshes',/,
     +          ' iz ',5x, ' dsl ',/,(i3,5x,e15.6))
 2000 format(10x,' DSL = ',e15.6)
 2100 format(10x,' DSL coarse = ',e15.6)
 9000 format(' Search for zi above the height = ',e15.6,/,
     +       ' iz_min = ',i5)
      end
      subroutine set_lower_bc
c
c -------- for curvy lower boundary set
c          distribution of heat flux or other
c          scalar sources and sinks
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      do iscl = 1,nscl

      if(isfc(iscl) .eq. 0) then ! given qstar compute surface temp
c
      do j=iys,iye
      do i=1,nnx
         t3tbc(i,j,iscl)  = qstar(iscl)
         t_grnd(i,j,iscl) = tsfcc(iscl)
      enddo
      enddo
c
      else ! given surface temp compute qstar
c
      do j=iys,iye
      do i=1,nnx
         t3tbc(i,j,iscl)  = -99999.0
         t_grnd(i,j,iscl) = tsfcc(iscl)
      enddo
      enddo
c
      endif

      end do
c
      return
      end
      subroutine nblnk(word)
      parameter (nmax=304)
      character wordt*304, word*(*)
      nchar = len(word)
      if(nchar .gt. nmax) then
         write(6,6000) nchar,nmax
 6000    format(' TROUBLE, IN SR. NBLNK : NCHAR = ',i6,
     +          ' EXCEEDS NMAX = ',i6)
         stop
      endif
      jj = 0
      do j=1,nchar
         if(word(j:j) .ne. ' ') then
            jj = jj + 1
            wordt(jj:jj) = word(j:j)
         endif
         word(j:j) = ' '
      enddo
      do j=1,jj
         word(j:j) = wordt(j:j)
      enddo
c
      return
      end
      subroutine blnk(word)
      character word*(*)
      nchar = len(word)
      do j=1,nchar
         word(j:j) = ' '
      enddo
c
      return
      end
      subroutine iso(it)
c
CHECK new eddy viscosity ivis0 = 0 for testing
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
c
      real sfk(1:nnz)
c
c ---- get isotropy factor and scale it to match at the matching
c      height. uses boundary conditions from lower and upper.
c
      do iz=1,nnz
c        dfac(iz) = 1.0
         dfac(iz) = 0.0
         sfk(iz)  = 0.0
      enddo
      do iz=izs,ize
         dfac(iz) = 1.0
      enddo
c
c ------ set nmatch equal to fraction of initial zi in sr. random
c
c     nmatch = izi/2
c     nmatch = 16
      nmatch = 48
      do i=0,numprocs-1,ncpu_s
         if(nmatch .ge. iz_s(i) .and.
     +      nmatch .le. iz_e(i)) myid_newvis = i
      enddo
c
      do iz=izs,min(ize,nmatch)
         izp1 = iz + 1
         izm1 = iz - 1
c        weit = dzw(iz)/(dzw(iz) + dzw(izp1))
c        weit1 = 1.0 - weit
c
c ---- get fluctuating strains
c
         do j=iys,iye
         do i=1,nnx
            s11 = weitm(iz)*ux(i,j,iz)**2 + weit(iz)*ux(i,j,izp1)**2
            s22 = weitm(iz)*vy(i,j,iz)**2 + weit(iz)*vy(i,j,izp1)**2
            wz  = (w(i,j,iz)-w(i,j,izm1))*dzw_i(iz)
            wzp = (w(i,j,izp1)-w(i,j,iz))*dzw_i(izp1)
            s33 = weit(iz)*wzp**2 + weitm(iz)*wz**2
            s12 = weitm(iz)*(uy(i,j,iz) + vx(i,j,iz))**2 +
     +            weit(iz)*(uy(i,j,izp1) + vx(i,j,izp1))**2
            s13 = (((u(i,j,izp1) - u(i,j,iz) +
     +            u_mn(iz) - u_mn(izp1))*dzu_i(izp1) +
     +            wx(i,j,iz)))**2
            s23 = (((v(i,j,izp1) - v(i,j,iz) +
     +          v_mn(iz) - v_mn(izp1))*dzu_i(izp1) +
     +          wy(i,j,iz)))**2
            sfk(iz) = sfk(iz) + 2.0*(s11 + s22 + s33) +
     +                       s12 + s13 + s23
         enddo
         enddo
         sfk(iz) = sfk(iz)*fnxy
      enddo
      call mpi_sum_z(sfk,i_root,myid,nnz,1)
c
      do iz=izs,min(ize,nmatch)
         izp1 = iz + 1
         izm1 = iz - 1
c
         sfk(iz) = sqrt(sfk(iz))
         smk = sqrt((u_mn(izp1)-u_mn(iz))**2 +
     +              (v_mn(izp1)-v_mn(iz))**2)*abs(dzu_i(izp1))
         if(sfk(iz) .le. 0. .and. smk .le. 0.) then
           dfac(iz) = 1.0
         else
           dfac(iz) = sfk(iz)/(sfk(iz) + smk)
         endif
c     if(l_root) write(6,6001) iz,sfk(iz),smk,dfac(iz)
 6001 format(' iz = ',i3,' sfk = ',e15.6,
     +       ' smk = ',e15.6,' dfac = ',e15.6)
      enddo
c
c
c ---- rescale ratio to give unity at match height
c      and if nested grid match value at upper boundary
c      of coarser grid
c
      if(myid .eq. myid_newvis) then
         dfacm = dfac(nmatch)
      endif
c
      call mpi_bcast(dfacm,1,mpi_real8,
     +              myid_newvis,mpi_comm_world,ierr)
c
      do iz=izs,min(ize,nmatch)
         dfac(iz) = dfac(iz)/dfacm
         dfac(iz) = amax1(dfac(iz), 0.1)
         dfac(iz) = amin1(dfac(iz), 1.0)
      enddo
c
c --------- gather dfac on all processes for printing and use in tke_vis
c           use reduce and divide by number of slab cpus
c
      call mpi_sum_z(dfac,i_root,myid,nnz,1)
      fncpu_s = 1.0/float(ncpu_s)
      do iz=1,nnz
         dfac(iz) = dfac(iz)*fncpu_s
      enddo
c
c     call mpi_gath_root(dfac(izs),dfac(1),iz_s,iz_e,izs,ize,nnz,
c    +                   myid,numprocs,ncpu_s)
c
c     if(l_root) write(6,6000) nmatch,ivis,(iz,dfac(iz),iz=1,nnz)
 6000 format(' in sr. iso, nmatch = ',i3,/,
     +       ' ivis = ',i3,'iz',5x,'dfac',/,(i3,1x,e15.6))
c     write(nprt,3001) (iz,dfac(iz),iz=1,nnz)
 3001 format(' iz ',5x,' dfac ',/,(i5,e15.6))
      return
      end
      subroutine surfvis(it)
c
c FIX for transformed grids
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
      real xkvis(nnx,iys:iye), alwk(nnx,iys:iye)
c
      real send(3), buf(3)
c
      xksurf = 0.0
      viscon = 0.0
      vise   = 0.0
c
c ----------- only root process(es) compute
c
      if(iss .eq. 0) then

c     ck = 0.1
c     csmag = 0.18
c     xkmax  = dzdz/dt/5.
      iz   = 1
      izm1 = iz - 1
      izp1 = iz + 1
c     xkmax  = dzu(izp1)*dzu(izp1)/(5.0*dt)
      dz_i = dzu_i(izp1)
      if(iocean .eq. 1) then
         call sufto(it)
      else
         call suft(it)
      endif
      if(qstar(1) .eq. 0.) then
         zeta = 0.0
      else
         zeta = abs(z(1))/amonin
      endif
      if(ismlt .eq. 1) then
          call busngr(zeta,phim,phis,psim,psis)
      else
          call fzol(zeta,phim,phis,psim,psis)
      endif
      viscon = vk*abs(z(1))/(utau*phim)
      vise   = utau*vk*abs(z(1))/phim
c
c ---- get special value at z1 to match with surface layer
c
      uws = 0.0
      vws = 0.0
      do iy=iys,iye
      do ix=1,nnx
         uws = uws + 0.5*(u(ix,iy,iz)-u_mn(iz) +
     +         u(ix,iy,izp1) - u_mn(izp1))*(w(ix,iy,iz)-w_mn(iz))
         vws = vws + 0.5*(v(ix,iy,iz)-v_mn(iz) +
     +         v(ix,iy,izp1) - v_mn(izp1))*(w(ix,iy,iz)-w_mn(iz))
      enddo
      enddo
      uws = uws*fnxy
      vws = vws*fnxy
c
c ---- get average fluctuating eddy viscsoity
c
      do iy=iys,iye
      do ix=1,nnx
         e(ix,iy,iz)=amax1(e(ix,iy,iz),sml_eg)
      enddo
      enddo
      dslk = amin1(dsl,vk*abs(z(iz))/csmag)
c     stabmin = 1.e-12
c     almin = 0.0001*dsl
      almin = almin_c*dsl_z(iz)
      do iy=iys,iye
      do ix=1,nnx
         alwk(ix,iy)=dslk
c
c --------no stability corrected length scales when
c         new eddy viscosity is on
c
c         stab=batag*(t(ix,iy,1,izp1)-t(ix,iy,1,iz))*dz_i
c         if(stab.gt.stabmin) then
c           als = stab_c*sqrt(e(ix,iy,iz)/stab)
c           alwk(ix,iy) = amin1(dslk,als)
c         endif
c         alwk(ix,iy)=amax1(almin,alwk(ix,iy))
         xkvis(ix,iy)=ck*alwk(ix,iy)*sqrt(e(ix,iy,iz))*dfac(1)
c        xkvis(ix,iy)=amin1(xkvis(ix,iy),xkmax)
      enddo
      enddo
c
c ---- get average viscosity
c
      xkavg = 0.0
      do iy=iys,iye
      do ix=1,nnx
         xkavg = xkavg + xkvis(ix,iy)
      enddo
      enddo
      xkavg = xkavg*fnxy
c
      buf(1) = uws
      buf(2) = vws
      buf(3) = xkavg
      call mpi_sum_xy(buf,myid,iss,ise,3)
c
      uws   = buf(1)
      vws   = buf(2)
      xkavg = buf(3)
c
      xkz1 = vise - sqrt(uws**2 + vws**2)*viscon
      xksurf =  xkz1 - xkavg
      xksurf = amax1(xksurf,0.0)
      xksurf = amin1(xksurf,vise)
c     if(l_root) write(6,6000) dfac(1), xkavg, xkz1, vise, xksurf
 6000 format(' dfac = ',e12.4,' xkavg = ',e12.4,' xkz1 = ',e12.4,/,
     +       ' vise = ',e12.4,' xksurf = ',e12.4)
c
      endif
c
c ---------- broadcast values to other processes
c
      send(1) = xksurf
      send(2) = viscon
      send(3) = vise
c
      call mpi_bcast(send,3,mpi_real8,
     +              i_root,mpi_comm_world,ierr)
c
      xksurf = send(1)
      viscon = send(2)
      vise   = send(3)
c
      return
      end
      subroutine comp1(istep,it)
c
c ----- 3-order runge-kutta time stepping and monotone scalar fluxes in x,y,z.
c       designed to use mpi in x & y directions.
c       variables advanced are (u,v,w,e,theta)/J
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      parameter(js = 7, ns = 3, nstat = js + ns*nscl)
      real stat(1:nnz,nstat)
c
c ------ temp arrays to hold rhs from step n-1 and
c        field variables from step n
c
      real urhs(nnx,iys:iye,izs:ize),
     +     vrhs(nnx,iys:iye,izs:ize),
     +     wrhs(nnx,iys:iye,izs:ize),
     +     erhs(nnx,iys:iye,izs:ize),
     +     trhs(nnx,iys:iye,nscl,izs:ize)
c
c ------------ build rhs with u/J from last time step
c              when the routine is entered (u,v,w,e,theta) are
c              free of (1/J)
c
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            urhs(i,j,k) = u(i,j,k)*gm(i,j,k,1) + dtzeta*r1(i,j,k)
            vrhs(i,j,k) = v(i,j,k)*gm(i,j,k,1) + dtzeta*r2(i,j,k)
            wrhs(i,j,k) = w(i,j,k)*gm(i,j,k,1) + dtzeta*r3(i,j,k)
            erhs(i,j,k) = e(i,j,k)*gm(i,j,k,1) + dtzeta*r5(i,j,k)
         enddo
         enddo
      enddo
c
      do k=izs,ize
         do l=1,nscl
         do j=iys,iye
         do i=1,nnx
            trhs(i,j,l,k) = t(i,j,l,k)*gm(i,j,k,1) + dtzeta*r4(i,j,l,k)
         enddo
         enddo
         enddo
      enddo
c
c --------- get viscosity and rhs of (e,u,v,w)-equations
c           at next step
c
      if (iDNS .eq. 1) then
         call dns_vis
         call rhs_uvw_dns(istep)
      else
         call tke_vis(istep)
         call rhs_uvw(istep)
      endif
c
c -------- evaluate rhs of scalar equations
c
      do l=1,nscl
         if (iDNS .eq. 1) then
            call rhs_scl_dns(istep,l)
         else
            call rhs_scl(istep,l)
         endif
      enddo
c
c ---------- gather stat sums on root processor
c            using mpi_reduction over all processors
c
      if(istep .eq. 1) then
c
        do j=1,nstat
        do k=1,nnz
           stat(k,j) = 0.0
        enddo
        enddo
        do k=izs,ize
           stat(k,1) = uwsb(k)
           stat(k,2) = vwsb(k)
           stat(k,3) = wwsb(k)
           stat(k,4) = tr_tau(k)
           stat(k,5) = triz(k)
           stat(k,6) = shrz(k)
           stat(k,7) = t_diss(k)
        enddo
        m1 = js
        m2 = js + nscl
        m3 = js + 2*nscl
        do l=1,nscl
           do k=izs,ize
              stat(k,m1+l) = utsb(k,l)
              stat(k,m2+l) = vtsb(k,l)
              stat(k,m3+l) = wtsb(k,l)
           enddo
        enddo
        call mpi_sum_z(stat(1,1),i_root,myid,nstat*nnz,1)
        do k=1,nnz
           uwsb(k)   = stat(k,1)
           vwsb(k)   = stat(k,2)
           wwsb(k)   = stat(k,3)
           tr_tau(k) = stat(k,4)
           triz(k)   = stat(k,5)
           shrz(k)   = stat(k,6)
           t_diss(k) = stat(k,7)
        enddo
        do l=1,nscl
           do k=1,nnz
              utsb(k,l) = stat(k,m1+l)
              vtsb(k,l) = stat(k,m2+l)
              wtsb(k,l) = stat(k,m3+l)
           enddo
        enddo
        do k=1,nnz
           buyz(k) = batag*wtsb(k,1)
        enddo
c
c -------- end if block
c
      endif
c
c ------- save old rhs in field variables for RK-advancement
c         add new rhs for computation of flux velocities
c         at this stage field variables are (u,v,w,e,theta)/J
c
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            u(i,j,k) = urhs(i,j,k) + dtgama*r1(i,j,k)
            v(i,j,k) = vrhs(i,j,k) + dtgama*r2(i,j,k)
            w(i,j,k) = wrhs(i,j,k) + dtgama*r3(i,j,k)
            e(i,j,k) = erhs(i,j,k) + dtgama*r5(i,j,k)
         enddo
         enddo
      enddo
      do k=izs,ize
         do l=1,nscl
         do j=iys,iye
         do i=1,nnx
            t(i,j,l,k) = trhs(i,j,l,k) + dtgama*r4(i,j,l,k)
         enddo
         enddo
         enddo
      enddo
c
c --------- generate rhs of contravariant flux velocities equations
c           at cell faces
c
      call rhs_flux_vel(istep)
c
      return
      end
      subroutine rhs_flux_vel(istep)
c
c ------ get rhs of intermediate velocity (uf,vf,wf)
c        by interpolating from cell center values.
c        extrapolate from interior to get surface value for wf.
c        assumes that rhs of the (u,v,w) equations are
c        multiplied by (1/J) from comp1, ie,
c        the right sides used in the interpolation here
c        are (u/J,v/J,w/J).
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
      real fs(nnx,iys:iye,3), fr(nnx,iys:iye,3)
      integer istatus(mpi_status_size)
      real wfi_sum(0:nnz)
c
c --------- Send updated (u,v,w) (from comp1)
c           to processor below the current myid to
c           get the contravariant flux velocities.
c
      nb = myid - ncpu_s
      nt = myid + ncpu_s
      if(iss .eq. 0) then
         nb = mpi_proc_null
      endif
      if(ise .eq. numprocs-1) then
         nt = mpi_proc_null
      endif
      nsend = 3*nnx*(iye + 1 - iys)
      nrecv = nsend
      do j=iys,iye
      do i=1,nnx
         fs(i,j,1) = u(i,j,izs)
         fs(i,j,2) = v(i,j,izs)
         fs(i,j,3) = w(i,j,izs)
      enddo
      enddo
      call mpi_sendrecv(
     +     fs(1,iys,1),nsend,mpi_real8,nb,2,
     +     fr(1,iys,1),nrecv,mpi_real8,nt,2,
     +     mpi_comm_world,istatus,ierr)
      if(ise .ne. numprocs-1) then
         do j=iys,iye
         do i=1,nnx
            u(i,j,ize+1)  = fr(i,j,1)
            v(i,j,ize+1)  = fr(i,j,2)
            w(i,j,ize+1)  = fr(i,j,3)
         enddo
         enddo
      endif

c
c ---------- get mean of wfi
c
      do k=0,nnz
         wfi_sum(k) = 0.0
      enddo
c
      do k=izs,ize
         km1 = k - 1
         kp1 = k + 1
         if(k .eq. nnz) then
c
c ------------ don't worry about upper bc for wf,
c              correct value is used in psolver
c
           do j=iys,iye
           do i=1,nnx
              wfi(i,j,k) = wfbc(i,j,1)
              ufi(i,j,k) = u(i,j,k)
              vfi(i,j,k) = v(i,j,k)
           enddo
           enddo
         else
           do j=iys,iye
           do i=1,nnx
              ufi(i,j,k) = u(i,j,k)
              vfi(i,j,k) = v(i,j,k)
              wf_k       = w(i,j,k)*gm(i,j,k,4)  +
     +                     u(i,j,k)*gm(i,j,k,2)  +
     +                     v(i,j,k)*gm(i,j,k,3)
              wf_kp1     = w(i,j,kp1)*gm(i,j,kp1,4)  +
     +                     u(i,j,kp1)*gm(i,j,kp1,2)  +
     +                     v(i,j,kp1)*gm(i,j,kp1,3)
              wfi(i,j,k) = wf_k*weitm(k) + wf_kp1*weit(k)
           enddo
           enddo
         endif
c
         do j=iys,iye
         do i=1,nnx
            wfi_sum(k) = wfi_sum(k) + wfi(i,j,k)
         enddo
         enddo
         wfi_sum(k) = wfi_sum(k)*fnxy
c
c --------- linear extrapolation to get wfi surface value
c
         if(k == 1) then
            w1 = -dzw(k)/(dzw(k) + dzw(kp1))
            w2 = 1.0 - w1
            do j=iys,iye
            do i=1,nnx
               ufi(i,j,km1) = u(i,j,k)
               vfi(i,j,km1) = v(i,j,k)
               wf_k         = w(i,j,k)*gm(i,j,k,4)  +
     +                        u(i,j,k)*gm(i,j,k,2)  +
     +                        v(i,j,k)*gm(i,j,k,3)
               wf_kp1       = w(i,j,kp1)*gm(i,j,kp1,4)  +
     +                        u(i,j,kp1)*gm(i,j,kp1,2)  +
     +                        v(i,j,kp1)*gm(i,j,kp1,3)
               wfi(i,j,km1) = wf_k*w2 + wf_kp1*w1
               wfi_sum(km1) = wfi_sum(km1) + wfi(i,j,km1)
            enddo
            enddo
            wfi_sum(km1) = wfi_sum(km1)*fnxy
         endif
c
      enddo
c
      call mpi_sum_z(wfi_sum(0),i_root,myid,nnz+1,1)
c
c -------  take out mean here
c
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            wfi(i,j,k) = wfi(i,j,k) - wfi_sum(k)
         enddo
         enddo
      enddo
c
      if(izs == 1) then
         k = 0
         do j=iys,iye
         do i=1,nnx
            wfi(i,j,k) = wfi(i,j,k) - wfi_sum(k)
         enddo
         enddo
      endif
c
c     if(izs .eq. 1) then
c       write(nprt,7002) istep
c7002   format(' 7002 rhs flux istage = ',i3)
c       k = 1
c       km1 = k - 1
c       kp1 = k + 1
c       j   = iys
c       write(nprt,7001) (i,wfi(i,j,km1), wfi(i,j,k), wfi(i,j,k+1),
c    +                      ufi(i,j,k), vfi(i,j,k),i=1,nnx)
c7001   format(' 7001 rhs_flux',/,
c    +  ' i ',5x,' wfi 0',5x,' wfi 1',5x,' wfi 2',5x,
c    +  ' ufi 1 ',5x,' vfi ',5x,/,(i5,5e15.6))
c     endif
c
      return
      end
      subroutine rhs_uvw(istep)
c
c ---------- get right hand sides of (u,v,w) equations
c            for pencil size (nnx, iys:iye, izs:ize).
c            rhs are (r1,r2,r3)/J, but the field
c            variables used to build rhs are (u,v,w,e,theta).
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      real fntd(nnx,iys:iye,izs:ize)
      real fnt1(nnx,iys:iye), fnt2(nnx,iys:iye), fnt3(nnx,iys:iye)
      real fnt4(nnx,iys:iye), fnt5(nnx,iys:iye)
      real tauij(nnx,iys:iye,6,2), sij(nnx,iys:iye,6)
      real pm_z(1:nnz)
c
      if(iskew .eq. 0) then
         upfac = 0.5
      else
         upfac = 1.0
      endif
c
c --------- special first grid point arrays, to avoid repeating
c
      if(izs .eq. 1) then
         k   = 1
         kp1 = k + 1
         km1 = k - 1
         do j=iys,iye
         do i=1,nnx
            fnt4(i,j) = weit(k)*gm(i,j,kp1,1) + weitm(k)*gm(i,j,k,1)
            fnt4(i,j) = 2.0*(weit(k)*vis_m(i,j,kp1) +
     +                       weitm(k)*vis_m(i,j,k))*fnt4(i,j)
            fnt5(i,j) = weit(km1)*gm(i,j,k,1) + weitm(km1)*gm(i,j,km1,1)
         enddo
         enddo
      endif
c
      do k=izs,ize
c
      km1 = k - 1
      kp1 = k + 1
c
c ------------ geostrophic wind
c
c     u_geo = -fcor*vg(k)
c     v_geo =  fcor*(ug(k)-ugal)
c
c ----------- for wind tunnel type simulations, account for box height
c
      u_geo = dpdx_con_z
      v_geo = 0.0
c
      do j=iys,iye
      do i=1,nnx
c
c ----------- boussinesq approximation for buoyancy (no hydro static piece)
c
         w_buy = batag*(t(i,j,1,k) - t00)
c
c ------------ coriolis, vertical and horizontal components
c
         u_cor     =  fcor*v(i,j,k) - fcor_h*w(i,j,k)
         v_cor     = -fcor*u(i,j,k)
         w_cor     =  fcor_h*u(i,j,k)
c
c ------------ build rhs times 1/J
c
         r1(i,j,k) = (u_cor + u_geo)*gm(i,j,k,1)
         r2(i,j,k) = (v_cor + v_geo)*gm(i,j,k,1)
         r3(i,j,k) = (w_cor + w_buy)*gm(i,j,k,1)
      enddo
      enddo
c
      enddo ! end z-loop
c
c ---------- advection for skew symmetric form
c
      if(iskew .eq. 0) then    ! skew symmetric form no grid speed
c
      do k=izs,ize
      km1 = k - 1
      kp1 = k + 1
      do j=iys,iye
      do i=1,nnx
         uzm     = (u(i,j,k) - u(i,j,km1))*dzu_i(k)
         vzm     = (v(i,j,k) - v(i,j,km1))*dzu_i(k)
         wzm     = (w(i,j,k) - w(i,j,km1))*dzu_i(k)
         uz      = (u(i,j,kp1) - u(i,j,k))*dzu_i(kp1)
         vz      = (v(i,j,kp1) - v(i,j,k))*dzu_i(kp1)
         wz      = (w(i,j,kp1) - w(i,j,k))*dzu_i(kp1)
c
         u_avg_u = u(i,j,k)*weitm(k) + u(i,j,kp1)*weit(k)
         v_avg_u = v(i,j,k)*weitm(k) + v(i,j,kp1)*weit(k)
         w_avg_u = w(i,j,k)*weitm(k) + w(i,j,kp1)*weit(k)
         u_avg_l = u(i,j,k)*weit(km1) + u(i,j,km1)*weitm(km1)
         v_avg_l = v(i,j,k)*weit(km1) + v(i,j,km1)*weitm(km1)
         w_avg_l = w(i,j,k)*weit(km1) + w(i,j,km1)*weitm(km1)
c
         u_adv = 0.5*(uf(i,j,k)*ux(i,j,k) + vf(i,j,k)*uy(i,j,k) +
     +           (wf(i,j,k)*u_avg_u - wf(i,j,km1)*u_avg_l)*dzw_i(k)) +
     +           0.25*(wf(i,j,k)*uz + wf(i,j,km1)*uzm)
         v_adv = 0.5*(uf(i,j,k)*vx(i,j,k) + vf(i,j,k)*vy(i,j,k) +
     +           (wf(i,j,k)*v_avg_u - wf(i,j,km1)*v_avg_l)*dzw_i(k)) +
     +           0.25*(wf(i,j,k)*vz + wf(i,j,km1)*vzm)
         w_adv = 0.5*(uf(i,j,k)*wx(i,j,k) + vf(i,j,k)*wy(i,j,k) +
     +           (wf(i,j,k)*w_avg_u - wf(i,j,km1)*w_avg_l)*dzw_i(k)) +
     +           0.25*(wf(i,j,k)*wz + wf(i,j,km1)*wzm)
c
c ------------ totals
c
         r1(i,j,k) = r1(i,j,k) - u_adv
         r2(i,j,k) = r2(i,j,k) - v_adv
         r3(i,j,k) = r3(i,j,k) - w_adv
      enddo
      enddo
      enddo ! end z-loop
c
      else   ! flux form with grid speed
c
      do k=izs,ize
      km1 = k - 1
      kp1 = k + 1
      do j=iys,iye
      do i=1,nnx
c
         u_avg_u = u(i,j,k)*weitm(k) + u(i,j,kp1)*weit(k)
         v_avg_u = v(i,j,k)*weitm(k) + v(i,j,kp1)*weit(k)
         w_avg_u = w(i,j,k)*weitm(k) + w(i,j,kp1)*weit(k)
         u_avg_l = u(i,j,k)*weit(km1) + u(i,j,km1)*weitm(km1)
         v_avg_l = v(i,j,k)*weit(km1) + v(i,j,km1)*weitm(km1)
         w_avg_l = w(i,j,k)*weit(km1) + w(i,j,km1)*weitm(km1)
c
         wflx_u  = wf(i,j,k) - gm(i,j,k,5)
         wflx_l  = wf(i,j,km1) - gm(i,j,km1,5)
         u_adv   = (wflx_u*u_avg_u - wflx_l*u_avg_l)*dzw_i(k)
         v_adv   = (wflx_u*v_avg_u - wflx_l*v_avg_l)*dzw_i(k)
         w_adv   = (wflx_u*w_avg_u - wflx_l*w_avg_l)*dzw_i(k)
c
c ------------ totals
c
         r1(i,j,k) = r1(i,j,k) - u_adv
         r2(i,j,k) = r2(i,j,k) - v_adv
         r3(i,j,k) = r3(i,j,k) - w_adv
      enddo
      enddo
      enddo ! end z-loop
c
      endif ! end of skew check
c
c ----------- SGS cell center x-contributions + resolved fluxes
c             advective terms can be written in skew symmetric
c             or flux form. treat first point special if at lower boundary,
c             average boundary and upper cell tau flux
c
      do k=izs,ize
         if(k .eq. 1) then
            call strain_u(sij,k)
            do j=iys,iye
            do i=1,nnx
               tauij(i,j,1,2) = 0.5*(fnt4(i,j)*sij(i,j,1) -
     +                               t11bc(i,j)*fnt5(i,j))
               tauij(i,j,2,2) = 0.5*(fnt4(i,j)*sij(i,j,2) -
     +                               t12bc(i,j)*fnt5(i,j))
               tauij(i,j,3,2) = 0.5*(fnt4(i,j)*sij(i,j,3) -
     +                               t13bc(i,j)*fnt5(i,j))
            enddo
            enddo
         else
            call strain_c(sij,k)
            do j=iys,iye
            do i=1,nnx
               tauij(i,j,1,2) = 2.0*vis_m(i,j,k)*sij(i,j,1)*gm(i,j,k,1)
               tauij(i,j,2,2) = 2.0*vis_m(i,j,k)*sij(i,j,2)*gm(i,j,k,1)
               tauij(i,j,3,2) = 2.0*vis_m(i,j,k)*sij(i,j,3)*gm(i,j,k,1)
            enddo
            enddo
         endif
         do j=iys,iye
         do i=1,nnx
            fnt1(i,j) = tauij(i,j,1,2) - upfac*uf(i,j,k)*u(i,j,k)
            fnt2(i,j) = tauij(i,j,2,2) - upfac*uf(i,j,k)*v(i,j,k)
            fnt3(i,j) = tauij(i,j,3,2) - upfac*uf(i,j,k)*w(i,j,k)
         enddo
         enddo
         call xderivp(fnt1(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call xderivp(fnt2(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call xderivp(fnt3(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         do j=iys,iye
         do i=1,nnx
            r1(i,j,k) = r1(i,j,k) + fnt1(i,j)
            r2(i,j,k) = r2(i,j,k) + fnt2(i,j)
            r3(i,j,k) = r3(i,j,k) + fnt3(i,j)
         enddo
         enddo
      enddo
c
c ----------- SGS z-contributions
c
      if(izs .eq. 1) then
         do j=iys,iye
         do i=1,nnx
            tauij(i,j,1,2) = -t11bc(i,j)
            tauij(i,j,2,2) = -t12bc(i,j)
            tauij(i,j,3,2) = -t13bc(i,j)
            tauij(i,j,4,2) = -t22bc(i,j)
            tauij(i,j,5,2) = -t23bc(i,j)
            tauij(i,j,6,2) = -t33bc(i,j)
         enddo
         enddo
      else
         call strain_u(sij,izs-1)
         do n=1,6
         do j=iys,iye
         do i=1,nnx
            fac   = 2.0*(weit(izs-1)*vis_m(i,j,izs) +
     +              weitm(izs-1)*vis_m(i,j,izs-1))
            tauij(i,j,n,2) = fac*sij(i,j,n)
         enddo
         enddo
         enddo
      endif
      do k=izs,ize
         kp1 = k + 1
         km1 = k - 1
         do n=1,6
         do j=iys,iye
         do i=1,nnx
            tauij(i,j,n,1) = tauij(i,j,n,2)
         enddo
         enddo
         enddo
         call strain_u(sij,k)
         do n=1,6
         do j=iys,iye
         do i=1,nnx
            fac   = 2.0*(weit(k)*vis_m(i,j,kp1) +
     +                   weitm(k)*vis_m(i,j,k))
            tauij(i,j,n,2) = fac*sij(i,j,n)
         enddo
         enddo
         enddo
         uw_sum = 0.0
         vw_sum = 0.0
         ww_sum = 0.0
         do j=iys,iye
         do i=1,nnx
            gm_xu     = weit(k)*gm(i,j,kp1,1)*gm(i,j,kp1,2) +
     +                  weitm(k)*gm(i,j,k,1)*gm(i,j,k,2)
            gm_xl     = weit(km1)*gm(i,j,k,1)*gm(i,j,k,2) +
     +                  weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,2)
            gm_yu     = weit(k)*gm(i,j,kp1,1)*gm(i,j,kp1,3) +
     +                  weitm(k)*gm(i,j,k,1)*gm(i,j,k,3)
            gm_yl     = weit(km1)*gm(i,j,k,1)*gm(i,j,k,3) +
     +                  weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,3)
            tau_up    = tauij(i,j,1,2)*gm_xu + tauij(i,j,2,2)*gm_yu +
     +                  tauij(i,j,3,2)
            tau_lo    = tauij(i,j,1,1)*gm_xl + tauij(i,j,2,1)*gm_yl +
     +                  tauij(i,j,3,1)
            uw_sum    = uw_sum + tau_up
            r1(i,j,k) = r1(i,j,k) + (tau_up - tau_lo)*dzw_i(k)
            tau_up    = tauij(i,j,2,2)*gm_xu + tauij(i,j,4,2)*gm_yu +
     +                  tauij(i,j,5,2)
            tau_lo    = tauij(i,j,2,1)*gm_xl + tauij(i,j,4,1)*gm_yl +
     +                  tauij(i,j,5,1)
            vw_sum    = vw_sum + tau_up
            r2(i,j,k) = r2(i,j,k) + (tau_up - tau_lo)*dzw_i(k)
            tau_up    = tauij(i,j,3,2)*gm_xu + tauij(i,j,5,2)*gm_yu +
     +                  tauij(i,j,6,2)
            tau_lo    = tauij(i,j,3,1)*gm_xl + tauij(i,j,5,1)*gm_yl +
     +                  tauij(i,j,6,1)
            ww_sum    = ww_sum + tau_up
            r3(i,j,k) = r3(i,j,k) + (tau_up - tau_lo)*dzw_i(k)
         enddo
         enddo
c
c -------- save SGS fluxes for printout
c
      if(istep .eq. 1) then
         tr_tau(k) = 0.0
         do j=iys,iye
         do i=1,nnx
            ufluc    = (u(i,j,kp1) - uxym(kp1))*weit(k) +
     +                 (u(i,j,k) - uxym(k))*weitm(k)
            vfluc    = (v(i,j,kp1) - vxym(kp1))*weit(k) +
     +                 (v(i,j,k) - vxym(k))*weitm(k)
            tr_tau(k) = tr_tau(k) -
     +                 tauij(i,j,3,2)*ufluc - tauij(i,j,5,2)*vfluc
         enddo
         enddo
         uwsb(k)   = -uw_sum*fnxy
         vwsb(k)   = -vw_sum*fnxy
         wwsb(k)   = -ww_sum*fnxy
         tr_tau(k) = tr_tau(k)*fnxy
      endif
c
c ---------- end z loop
c
      enddo
c
c ---------- SGS fluxes tau_12, tau_22, tau_23 that depend on
c            y-derivatives at cell centers. special treatment
c            of first point at lower boundary. assumes strain
c            at k = 1 is at upper face, see strain_12,_22,_23
c
      do k=izs,ize
         call strain_12(sij(1,iys,1),k)
         if(k .eq. 1) then
            do j=iys,iye
            do i=1,nnx
               tau_c       = 0.5*(fnt4(i,j)*sij(i,j,1) -
     +                            t12bc(i,j)*fnt5(i,j))
               fntd(i,j,k) = tau_c  - upfac*vf(i,j,k)*u(i,j,k)
            enddo
            enddo
         else
            do j=iys,iye
            do i=1,nnx
               tau_c       = 2.0*gm(i,j,k,1)*vis_m(i,j,k)*sij(i,j,1)
               fntd(i,j,k) = tau_c - upfac*vf(i,j,k)*u(i,j,k)
            enddo
            enddo
         endif
      enddo
      call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            r1(i,j,k)   = r1(i,j,k) + fntd(i,j,k)
         enddo
         enddo
      enddo
c
      do k=izs,ize
         call strain_22(sij(1,iys,1),k)
         if(k .eq. 1) then
            do j=iys,iye
            do i=1,nnx
               tau_c       = 0.5*(fnt4(i,j)*sij(i,j,1) -
     +                            t22bc(i,j)*fnt5(i,j))
               fntd(i,j,k) = tau_c  - upfac*vf(i,j,k)*v(i,j,k)
            enddo
            enddo
         else
            do j=iys,iye
            do i=1,nnx
               tau_c       = 2.0*gm(i,j,k,1)*vis_m(i,j,k)*sij(i,j,1)
               fntd(i,j,k) = tau_c - upfac*vf(i,j,k)*v(i,j,k)
            enddo
            enddo
         endif
      enddo
      call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            r2(i,j,k)   = r2(i,j,k) + fntd(i,j,k)
         enddo
         enddo
      enddo
c
      do k=izs,ize
         call strain_23(sij(1,iys,1),k)
         if(k .eq. 1) then
            do j=iys,iye
            do i=1,nnx
               tau_c       = 0.5*(fnt4(i,j)*sij(i,j,1) -
     +                            t23bc(i,j)*fnt5(i,j))
               fntd(i,j,k) = tau_c - upfac*vf(i,j,k)*w(i,j,k)
            enddo
            enddo
         else
            do j=iys,iye
            do i=1,nnx
               tau_c       = 2.0*gm(i,j,k,1)*vis_m(i,j,k)*sij(i,j,1)
               fntd(i,j,k) = tau_c - upfac*vf(i,j,k)*w(i,j,k)
            enddo
            enddo
         endif
      enddo
      call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            r3(i,j,k) = r3(i,j,k) + fntd(i,j,k)
         enddo
         enddo
      enddo
c
c --------- get average of r3 as estimate of mean pressure gradient
c           from turbulence and add it to right hand sides of all equations
c
      do k=1,nnz
         pm_z(k) = 0.0
      enddo
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            pm_z(k) = pm_z(k) + r3(i,j,k)
         enddo
         enddo
         pm_z(k) = pm_z(k)*fnxy
      enddo
      call mpi_sum_z(pm_z,i_root,myid,nnz,1)
c
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            r3(i,j,k) = r3(i,j,k) - pm_z(k)
         enddo
         enddo
      enddo
c
c ------- integrate p-gradient to get values at wf points
c
      pm(0) = p_surf
      do k=1,nnz
         pm(k) = pm(k-1) + pm_z(k)*dzw(k)
      enddo
c
      if(no_pgrad == 1) then
c
c --------- add mean p-gradients to right hand sides
c           of u and v equations using non-conservative form
c           to avoid specifying the surface pressure
c
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            r1(i,j,k) = r1(i,j,k) - pm_z(k)*gm(i,j,k,1)*gm(i,j,k,2)
            r2(i,j,k) = r2(i,j,k) - pm_z(k)*gm(i,j,k,1)*gm(i,j,k,3)
         enddo
         enddo
      enddo
c
      endif ! check for no mean dpdx in (u,v)
c
c --------- leave with (r1,r2,r3)/J
c
      return
      end
      subroutine rhs_uvw_dns(istep)
c
c ---------- get right hand sides of (u,v,w) equations
c            for pencil size (nnx, iys:iye, izs:ize).
c            rhs are (r1,r2,r3)/J, but the field
c            variables used to build rhs are (u,v,w,e,theta).
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      real fntd(nnx,iys:iye,izs:ize)
      real fnt1(nnx,iys:iye), fnt2(nnx,iys:iye), fnt3(nnx,iys:iye)
      real fnt4(nnx,iys:iye), fnt5(nnx,iys:iye)
      real tauij(nnx,iys:iye,6,2), sij(nnx,iys:iye,6)
      real pm_z(1:nnz)
      real sfc_flx(2)
c
      if(iskew .eq. 0) then
         upfac = 0.5
      else
         upfac = 1.0
      endif
      sfc_flx = 0.0
c
c --------- special first grid point arrays, to avoid repeating
c
      if(izs .eq. 1) then
         k   = 1
         kp1 = k + 1
         km1 = k - 1
         do j=iys,iye
         do i=1,nnx
            fnt4(i,j) = weit(k)*gm(i,j,kp1,1) + weitm(k)*gm(i,j,k,1)
            fnt4(i,j) = 2.0*(weit(k)*vis_m(i,j,kp1) +
     +                       weitm(k)*vis_m(i,j,k))*fnt4(i,j)
            fnt5(i,j) = weit(km1)*gm(i,j,k,1) + weitm(km1)*gm(i,j,km1,1)
         enddo
         enddo
      endif
c
      do k=izs,ize
c
      km1 = k - 1
      kp1 = k + 1
c
c ------------ geostrophic wind
c
c     u_geo = -fcor*vg(k)
c     v_geo =  fcor*(ug(k)-ugal)
c
c ----------- for wind tunnel type simulations, account for box height
c
      u_geo = dpdx_con_z
      v_geo = 0.0
c
      do j=iys,iye
      do i=1,nnx
c
c ----------- boussinesq approximation for buoyancy (no hydro static piece)
c
         w_buy = batag*(t(i,j,1,k) - t00)
c
c ------------ coriolis, vertical and horizontal components
c
         u_cor     =  fcor*v(i,j,k) - fcor_h*w(i,j,k)
         v_cor     = -fcor*u(i,j,k)
         w_cor     =  fcor_h*u(i,j,k)
c
c ------------ build rhs times 1/J
c
         r1(i,j,k) = (u_cor + u_geo)*gm(i,j,k,1)
         r2(i,j,k) = (v_cor + v_geo)*gm(i,j,k,1)
         r3(i,j,k) = (w_cor + w_buy)*gm(i,j,k,1)
      enddo
      enddo
c
      enddo ! end z-loop
c
c ---------- advection for skew symmetric form
c
      if(iskew .eq. 0) then    ! skew symmetric form no grid speed
c
      do k=izs,ize
      km1 = k - 1
      kp1 = k + 1
      do j=iys,iye
      do i=1,nnx
         uzm     = (u(i,j,k) - u(i,j,km1))*dzu_i(k)
         vzm     = (v(i,j,k) - v(i,j,km1))*dzu_i(k)
         wzm     = (w(i,j,k) - w(i,j,km1))*dzu_i(k)
         uz      = (u(i,j,kp1) - u(i,j,k))*dzu_i(kp1)
         vz      = (v(i,j,kp1) - v(i,j,k))*dzu_i(kp1)
         wz      = (w(i,j,kp1) - w(i,j,k))*dzu_i(kp1)
c
         u_avg_u = u(i,j,k)*weitm(k) + u(i,j,kp1)*weit(k)
         v_avg_u = v(i,j,k)*weitm(k) + v(i,j,kp1)*weit(k)
         w_avg_u = w(i,j,k)*weitm(k) + w(i,j,kp1)*weit(k)
         u_avg_l = u(i,j,k)*weit(km1) + u(i,j,km1)*weitm(km1)
         v_avg_l = v(i,j,k)*weit(km1) + v(i,j,km1)*weitm(km1)
         w_avg_l = w(i,j,k)*weit(km1) + w(i,j,km1)*weitm(km1)
c
         u_adv = 0.5*(uf(i,j,k)*ux(i,j,k) + vf(i,j,k)*uy(i,j,k) +
     +           (wf(i,j,k)*u_avg_u - wf(i,j,km1)*u_avg_l)*dzw_i(k)) +
     +           0.25*(wf(i,j,k)*uz + wf(i,j,km1)*uzm)
         v_adv = 0.5*(uf(i,j,k)*vx(i,j,k) + vf(i,j,k)*vy(i,j,k) +
     +           (wf(i,j,k)*v_avg_u - wf(i,j,km1)*v_avg_l)*dzw_i(k)) +
     +           0.25*(wf(i,j,k)*vz + wf(i,j,km1)*vzm)
         w_adv = 0.5*(uf(i,j,k)*wx(i,j,k) + vf(i,j,k)*wy(i,j,k) +
     +           (wf(i,j,k)*w_avg_u - wf(i,j,km1)*w_avg_l)*dzw_i(k)) +
     +           0.25*(wf(i,j,k)*wz + wf(i,j,km1)*wzm)
c
c ------------ totals
c
         r1(i,j,k) = r1(i,j,k) - u_adv
         r2(i,j,k) = r2(i,j,k) - v_adv
         r3(i,j,k) = r3(i,j,k) - w_adv
      enddo
      enddo
      enddo ! end z-loop
c
      else   ! flux form with grid speed
c
      do k=izs,ize
      km1 = k - 1
      kp1 = k + 1
      do j=iys,iye
      do i=1,nnx
c
         u_avg_u = u(i,j,k)*weitm(k) + u(i,j,kp1)*weit(k)
         v_avg_u = v(i,j,k)*weitm(k) + v(i,j,kp1)*weit(k)
         w_avg_u = w(i,j,k)*weitm(k) + w(i,j,kp1)*weit(k)
         u_avg_l = u(i,j,k)*weit(km1) + u(i,j,km1)*weitm(km1)
         v_avg_l = v(i,j,k)*weit(km1) + v(i,j,km1)*weitm(km1)
         w_avg_l = w(i,j,k)*weit(km1) + w(i,j,km1)*weitm(km1)
c
         wflx_u  = wf(i,j,k) - gm(i,j,k,5)
         wflx_l  = wf(i,j,km1) - gm(i,j,km1,5)
         u_adv   = (wflx_u*u_avg_u - wflx_l*u_avg_l)*dzw_i(k)
         v_adv   = (wflx_u*v_avg_u - wflx_l*v_avg_l)*dzw_i(k)
         w_adv   = (wflx_u*w_avg_u - wflx_l*w_avg_l)*dzw_i(k)
c
c ------------ totals
c
         r1(i,j,k) = r1(i,j,k) - u_adv
         r2(i,j,k) = r2(i,j,k) - v_adv
         r3(i,j,k) = r3(i,j,k) - w_adv
      enddo
      enddo
      enddo ! end z-loop
c
      endif ! end of skew check
c
c ----------- SGS cell center x-contributions + resolved fluxes
c             advective terms can be written in skew symmetric
c             or flux form. treat first point special if at lower boundary,
c             average boundary and upper cell tau flux
c
      do k=izs,ize
!         if(k .eq. 1) then
!            call strain_u(sij,k)
!            do j=iys,iye
!            do i=1,nnx
!               tauij(i,j,1,2) = 0.5*(fnt4(i,j)*sij(i,j,1) -
!     +                               t11bc(i,j)*fnt5(i,j))
!               tauij(i,j,2,2) = 0.5*(fnt4(i,j)*sij(i,j,2) -
!     +                               t12bc(i,j)*fnt5(i,j))
!               tauij(i,j,3,2) = 0.5*(fnt4(i,j)*sij(i,j,3) -
!     +                               t13bc(i,j)*fnt5(i,j))
!            enddo
!            enddo
!         else
            call strain_c_dns(sij,k)
            do j=iys,iye
            do i=1,nnx
               tauij(i,j,1,2) = 2.0*vis_m(i,j,k)*sij(i,j,1)*gm(i,j,k,1)
               tauij(i,j,2,2) = 2.0*vis_m(i,j,k)*sij(i,j,2)*gm(i,j,k,1)
               tauij(i,j,3,2) = 2.0*vis_m(i,j,k)*sij(i,j,3)*gm(i,j,k,1)
            enddo
            enddo
!         endif
         do j=iys,iye
         do i=1,nnx
            fnt1(i,j) = tauij(i,j,1,2) - upfac*uf(i,j,k)*u(i,j,k)
            fnt2(i,j) = tauij(i,j,2,2) - upfac*uf(i,j,k)*v(i,j,k)
            fnt3(i,j) = tauij(i,j,3,2) - upfac*uf(i,j,k)*w(i,j,k)
         enddo
         enddo
         call xderivp(fnt1(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call xderivp(fnt2(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call xderivp(fnt3(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         do j=iys,iye
         do i=1,nnx
            r1(i,j,k) = r1(i,j,k) + fnt1(i,j)
            r2(i,j,k) = r2(i,j,k) + fnt2(i,j)
            r3(i,j,k) = r3(i,j,k) + fnt3(i,j)
         enddo
         enddo
      enddo
c
c ----------- SGS z-contributions
c
!      if(izs .eq. 1) then
!         do j=iys,iye
!         do i=1,nnx
!            tauij(i,j,1,2) = -t11bc(i,j)
!            tauij(i,j,2,2) = -t12bc(i,j)
!            tauij(i,j,3,2) = -t13bc(i,j)
!            tauij(i,j,4,2) = -t22bc(i,j)
!            tauij(i,j,5,2) = -t23bc(i,j)
!            tauij(i,j,6,2) = -t33bc(i,j)
!         enddo
!         enddo
!      else
         call strain_u(sij,izs-1)  !NOTE: strain_u does not need dns revision -- doesn't treat k=1 specially
         do n=1,6
         do j=iys,iye
         do i=1,nnx
            fac   = 2.0*(weit(izs-1)*vis_m(i,j,izs) +
     +              weitm(izs-1)*vis_m(i,j,izs-1))
            tauij(i,j,n,2) = fac*sij(i,j,n)
         enddo
         enddo
         enddo
!      endif
      do k=izs,ize
         kp1 = k + 1
         km1 = k - 1
         do n=1,6
         do j=iys,iye
         do i=1,nnx
            tauij(i,j,n,1) = tauij(i,j,n,2)
         enddo
         enddo
         enddo
         call strain_u(sij,k)
         do n=1,6
         do j=iys,iye
         do i=1,nnx
            fac   = 2.0*(weit(k)*vis_m(i,j,kp1) +
     +                   weitm(k)*vis_m(i,j,k))
            tauij(i,j,n,2) = fac*sij(i,j,n)
         enddo
         enddo
         enddo
         uw_sum = 0.0
         vw_sum = 0.0
         ww_sum = 0.0
         do j=iys,iye
         do i=1,nnx
            gm_xu     = weit(k)*gm(i,j,kp1,1)*gm(i,j,kp1,2) +
     +                  weitm(k)*gm(i,j,k,1)*gm(i,j,k,2)
            gm_xl     = weit(km1)*gm(i,j,k,1)*gm(i,j,k,2) +
     +                  weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,2)
            gm_yu     = weit(k)*gm(i,j,kp1,1)*gm(i,j,kp1,3) +
     +                  weitm(k)*gm(i,j,k,1)*gm(i,j,k,3)
            gm_yl     = weit(km1)*gm(i,j,k,1)*gm(i,j,k,3) +
     +                  weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,3)
            tau_up    = tauij(i,j,1,2)*gm_xu + tauij(i,j,2,2)*gm_yu +
     +                  tauij(i,j,3,2)
            tau_lo    = tauij(i,j,1,1)*gm_xl + tauij(i,j,2,1)*gm_yl +
     +                  tauij(i,j,3,1)
            uw_sum    = uw_sum + tau_up
            r1(i,j,k) = r1(i,j,k) + (tau_up - tau_lo)*dzw_i(k)

            if (k==1) then
            sfc_flx(1) = sfc_flx(1) + tau_lo
            end if

            tau_up    = tauij(i,j,2,2)*gm_xu + tauij(i,j,4,2)*gm_yu +
     +                  tauij(i,j,5,2)
            tau_lo    = tauij(i,j,2,1)*gm_xl + tauij(i,j,4,1)*gm_yl +
     +                  tauij(i,j,5,1)
            vw_sum    = vw_sum + tau_up

            if (k==1) then
            sfc_flx(2) = sfc_flx(2) + tau_lo
            end if

            r2(i,j,k) = r2(i,j,k) + (tau_up - tau_lo)*dzw_i(k)
            tau_up    = tauij(i,j,3,2)*gm_xu + tauij(i,j,5,2)*gm_yu +
     +                  tauij(i,j,6,2)
            tau_lo    = tauij(i,j,3,1)*gm_xl + tauij(i,j,5,1)*gm_yl +
     +                  tauij(i,j,6,1)
            ww_sum    = ww_sum + tau_up
            r3(i,j,k) = r3(i,j,k) + (tau_up - tau_lo)*dzw_i(k)

         enddo
         enddo
c
c -------- save SGS fluxes for printout
c
      if(istep .eq. 1) then
         tr_tau(k) = 0.0
         do j=iys,iye
         do i=1,nnx
            ufluc    = (u(i,j,kp1) - uxym(kp1))*weit(k) +
     +                 (u(i,j,k) - uxym(k))*weitm(k)
            vfluc    = (v(i,j,kp1) - vxym(kp1))*weit(k) +
     +                 (v(i,j,k) - vxym(k))*weitm(k)
            tr_tau(k) = tr_tau(k) -
     +                 tauij(i,j,3,2)*ufluc - tauij(i,j,5,2)*vfluc
         enddo
         enddo
         uwsb(k)   = -uw_sum*fnxy
         vwsb(k)   = -vw_sum*fnxy
         wwsb(k)   = -ww_sum*fnxy
         tr_tau(k) = tr_tau(k)*fnxy
      endif

      !As an aside, compute the values uwsfc and vwsfc for output:
      if (k==1) then
         call mpi_sum_xy(sfc_flx,myid,iss,ise,2)
         uwsfc = sfc_flx(1)*fnxy
         vwsfc = sfc_flx(2)*fnxy
      end if
c
c ---------- end z loop
c
      enddo
c
c ---------- SGS fluxes tau_12, tau_22, tau_23 that depend on
c            y-derivatives at cell centers. special treatment
c            of first point at lower boundary. assumes strain
c            at k = 1 is at upper face, see strain_12,_22,_23
c
      do k=izs,ize
         !call strain_12(sij(1,iys,1),k)
         call strain_12_dns(sij(1,iys,1),k)
!         if(k .eq. 1) then
!            do j=iys,iye
!            do i=1,nnx
!               tau_c       = 0.5*(fnt4(i,j)*sij(i,j,1) -
!     +                            t12bc(i,j)*fnt5(i,j))
!               fntd(i,j,k) = tau_c  - upfac*vf(i,j,k)*u(i,j,k)
!            enddo
!            enddo
!         else
            do j=iys,iye
            do i=1,nnx
               tau_c       = 2.0*gm(i,j,k,1)*vis_m(i,j,k)*sij(i,j,1)
               fntd(i,j,k) = tau_c - upfac*vf(i,j,k)*u(i,j,k)
            enddo
            enddo
!         endif
      enddo
      call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            r1(i,j,k)   = r1(i,j,k) + fntd(i,j,k)
         enddo
         enddo
      enddo
c
      do k=izs,ize
         !call strain_22(sij(1,iys,1),k)
         call strain_22_dns(sij(1,iys,1),k)
!         if(k .eq. 1) then
!            do j=iys,iye
!            do i=1,nnx
!               tau_c       = 0.5*(fnt4(i,j)*sij(i,j,1) -
!     +                            t22bc(i,j)*fnt5(i,j))
!               fntd(i,j,k) = tau_c  - upfac*vf(i,j,k)*v(i,j,k)
!            enddo
!            enddo
!         else
            do j=iys,iye
            do i=1,nnx
               tau_c       = 2.0*gm(i,j,k,1)*vis_m(i,j,k)*sij(i,j,1)
               fntd(i,j,k) = tau_c - upfac*vf(i,j,k)*v(i,j,k)
            enddo
            enddo
!         endif
      enddo
      call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            r2(i,j,k)   = r2(i,j,k) + fntd(i,j,k)
         enddo
         enddo
      enddo
c
      do k=izs,ize
         !call strain_23(sij(1,iys,1),k)
         call strain_23_dns(sij(1,iys,1),k)
!         if(k .eq. 1) then
!            do j=iys,iye
!            do i=1,nnx
!               tau_c       = 0.5*(fnt4(i,j)*sij(i,j,1) -
!     +                            t23bc(i,j)*fnt5(i,j))
!               fntd(i,j,k) = tau_c - upfac*vf(i,j,k)*w(i,j,k)
!            enddo
!            enddo
!         else
            do j=iys,iye
            do i=1,nnx
               tau_c       = 2.0*gm(i,j,k,1)*vis_m(i,j,k)*sij(i,j,1)
               fntd(i,j,k) = tau_c - upfac*vf(i,j,k)*w(i,j,k)
            enddo
            enddo
!         endif
      enddo
      call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            r3(i,j,k) = r3(i,j,k) + fntd(i,j,k)
         enddo
         enddo
      enddo
c
c --------- get average of r3 as estimate of mean pressure gradient
c           from turbulence and add it to right hand sides of all equations
c
      do k=1,nnz
         pm_z(k) = 0.0
      enddo
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            pm_z(k) = pm_z(k) + r3(i,j,k)
         enddo
         enddo
         pm_z(k) = pm_z(k)*fnxy
      enddo
      call mpi_sum_z(pm_z,i_root,myid,nnz,1)
c
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            r3(i,j,k) = r3(i,j,k) - pm_z(k)
         enddo
         enddo
      enddo
c
c ------- integrate p-gradient to get values at wf points
c
      pm(0) = p_surf
      do k=1,nnz
         pm(k) = pm(k-1) + pm_z(k)*dzw(k)
      enddo
c
      if(no_pgrad == 1) then
c
c --------- add mean p-gradients to right hand sides
c           of u and v equations using non-conservative form
c           to avoid specifying the surface pressure
c
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            r1(i,j,k) = r1(i,j,k) - pm_z(k)*gm(i,j,k,1)*gm(i,j,k,2)
            r2(i,j,k) = r2(i,j,k) - pm_z(k)*gm(i,j,k,1)*gm(i,j,k,3)
         enddo
         enddo
      enddo
c
      endif ! check for no mean dpdx in (u,v)
c
c --------- leave with (r1,r2,r3)/J
c
      return
      end
      subroutine rhs_uvw_old(istep)
c
c ---------- get right hand sides of (u,v,w) equations
c            for pencil size (nnx, iys:iye, izs:ize).
c            rhs are (r1,r2,r3)/J, but the field
c            variables used to build rhs are (u,v,w,e,theta).
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      real fntd(nnx,iys:iye,izs:ize)
      real fnt1(nnx,iys:iye), fnt2(nnx,iys:iye), fnt3(nnx,iys:iye)
      real tauij(nnx,iys:iye,6,2), sij(nnx,iys:iye,6)
      real r3_sum(1:nnz)
c
      if(iskew .eq. 0) then
         upfac = 0.5
      else
         upfac = 1.0
      endif
c
      do k=izs,ize
c
      km1 = k - 1
      kp1 = k + 1
c
c ------------ geostrophic wind
c
      u_geo = -fcor*vg(k)
      v_geo =  fcor*(ug(k)-ugal)
c
c ----------- for wind tunnel type simulations
c
c     u_geo = dpdx_con
c     v_geo = 0.0
c
      do j=iys,iye
      do i=1,nnx
c
c ------------ buoyancy (with hydrostatic part)
c              will correct at the bottom
c              of the routine
c
c        w_buy = batag*(t(i,j,1,k) - t_mn(k,1))
         w_buy = batag*t(i,j,1,k)
c
c ------------ coriolis, vertical and horizontal components
c
         u_cor     =  fcor*v(i,j,k) - fcor_h*w(i,j,k)
         v_cor     = -fcor*u(i,j,k)
         w_cor     =  fcor_h*u(i,j,k)
c
c ------------ build rhs times 1/J
c
         r1(i,j,k) = (u_cor + u_geo)*gm(i,j,k,1)
         r2(i,j,k) = (v_cor + v_geo)*gm(i,j,k,1)
         r3(i,j,k) = (w_cor + w_buy)*gm(i,j,k,1)
      enddo
      enddo
c
      enddo ! end z-loop
c
c ---------- advection for skew symmetric form, fix for moving grids
c
      if(iskew .eq. 0) then    ! skew symmetric form no grid speed
c
      do k=izs,ize
      km1 = k - 1
      kp1 = k + 1
      do j=iys,iye
      do i=1,nnx
         uzm     = (u(i,j,k) - u(i,j,km1))*dzu_i(k)
         vzm     = (v(i,j,k) - v(i,j,km1))*dzu_i(k)
         wzm     = (w(i,j,k) - w(i,j,km1))*dzu_i(k)
         uz      = (u(i,j,kp1) - u(i,j,k))*dzu_i(kp1)
         vz      = (v(i,j,kp1) - v(i,j,k))*dzu_i(kp1)
         wz      = (w(i,j,kp1) - w(i,j,k))*dzu_i(kp1)
c
         u_avg_u = u(i,j,k)*weitm(k) + u(i,j,kp1)*weit(k)
         v_avg_u = v(i,j,k)*weitm(k) + v(i,j,kp1)*weit(k)
         w_avg_u = w(i,j,k)*weitm(k) + w(i,j,kp1)*weit(k)
         u_avg_l = u(i,j,k)*weit(km1) + u(i,j,km1)*weitm(km1)
         v_avg_l = v(i,j,k)*weit(km1) + v(i,j,km1)*weitm(km1)
         w_avg_l = w(i,j,k)*weit(km1) + w(i,j,km1)*weitm(km1)
c
         u_adv = 0.5*(uf(i,j,k)*ux(i,j,k) + vf(i,j,k)*uy(i,j,k) +
     +           (wf(i,j,k)*u_avg_u - wf(i,j,km1)*u_avg_l)*dzw_i(k)) +
     +           0.25*(wf(i,j,k)*uz + wf(i,j,km1)*uzm)
         v_adv = 0.5*(uf(i,j,k)*vx(i,j,k) + vf(i,j,k)*vy(i,j,k) +
     +           (wf(i,j,k)*v_avg_u - wf(i,j,km1)*v_avg_l)*dzw_i(k)) +
     +           0.25*(wf(i,j,k)*vz + wf(i,j,km1)*vzm)
         w_adv = 0.5*(uf(i,j,k)*wx(i,j,k) + vf(i,j,k)*wy(i,j,k) +
     +           (wf(i,j,k)*w_avg_u - wf(i,j,km1)*w_avg_l)*dzw_i(k)) +
     +           0.25*(wf(i,j,k)*wz + wf(i,j,km1)*wzm)
c
c ------------ totals
c
         r1(i,j,k) = r1(i,j,k) - u_adv
         r2(i,j,k) = r2(i,j,k) - v_adv
         r3(i,j,k) = r3(i,j,k) - w_adv
      enddo
      enddo
      enddo ! end z-loop
c
      else   ! flux form with grid speed
c
      do k=izs,ize
      km1 = k - 1
      kp1 = k + 1
      do j=iys,iye
      do i=1,nnx
c
         u_avg_u = u(i,j,k)*weitm(k) + u(i,j,kp1)*weit(k)
         v_avg_u = v(i,j,k)*weitm(k) + v(i,j,kp1)*weit(k)
         w_avg_u = w(i,j,k)*weitm(k) + w(i,j,kp1)*weit(k)
         u_avg_l = u(i,j,k)*weit(km1) + u(i,j,km1)*weitm(km1)
         v_avg_l = v(i,j,k)*weit(km1) + v(i,j,km1)*weitm(km1)
         w_avg_l = w(i,j,k)*weit(km1) + w(i,j,km1)*weitm(km1)
c
         wflx_u  = wf(i,j,k) - gm(i,j,k,5)
         wflx_l  = wf(i,j,km1) - gm(i,j,km1,5)
         u_adv   = (wflx_u*u_avg_u - wflx_l*u_avg_l)*dzw_i(k)
         v_adv   = (wflx_u*v_avg_u - wflx_l*v_avg_l)*dzw_i(k)
         w_adv   = (wflx_u*w_avg_u - wflx_l*w_avg_l)*dzw_i(k)
c
c ------------ totals
c
         r1(i,j,k) = r1(i,j,k) - u_adv
         r2(i,j,k) = r2(i,j,k) - v_adv
         r3(i,j,k) = r3(i,j,k) - w_adv
      enddo
      enddo
      enddo ! end z-loop
c
      endif ! end of skew check
c
c ----------- SGS cell center x-contributions + resolved fluxes
c             advective terms written in skew symmetric form
c
      do k=izs,ize
         call strain_c(sij,k)
         do j=iys,iye
         do i=1,nnx
            fac       = 2.0*gm(i,j,k,1)*vis_m(i,j,k)
            fnt1(i,j) = fac*sij(i,j,1) - upfac*uf(i,j,k)*u(i,j,k)
            fnt2(i,j) = fac*sij(i,j,2) - upfac*uf(i,j,k)*v(i,j,k)
            fnt3(i,j) = fac*sij(i,j,3) - upfac*uf(i,j,k)*w(i,j,k)
         enddo
         enddo
         call xderivp(fnt1(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call xderivp(fnt2(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call xderivp(fnt3(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         do j=iys,iye
         do i=1,nnx
            r1(i,j,k) = r1(i,j,k) + fnt1(i,j)
            r2(i,j,k) = r2(i,j,k) + fnt2(i,j)
            r3(i,j,k) = r3(i,j,k) + fnt3(i,j)
         enddo
         enddo
      enddo
c
c ----------- SGS z-contributions
c
      if(izs .eq. 1) then
         do j=iys,iye
         do i=1,nnx
            tauij(i,j,1,2) = -t11bc(i,j)
            tauij(i,j,2,2) = -t12bc(i,j)
            tauij(i,j,3,2) = -t13bc(i,j)
            tauij(i,j,4,2) = -t22bc(i,j)
            tauij(i,j,5,2) = -t23bc(i,j)
            tauij(i,j,6,2) = -t33bc(i,j)
         enddo
         enddo
      else
         call strain_u(sij,izs-1)
         do n=1,6
         do j=iys,iye
         do i=1,nnx
            fac   = 2.0*(weit(izs-1)*vis_m(i,j,izs) +
     +              weitm(izs-1)*vis_m(i,j,izs-1))
            tauij(i,j,n,2) = fac*sij(i,j,n)
         enddo
         enddo
         enddo
      endif
      do k=izs,ize
         kp1 = k + 1
         km1 = k - 1
         do n=1,6
         do j=iys,iye
         do i=1,nnx
            tauij(i,j,n,1) = tauij(i,j,n,2)
         enddo
         enddo
         enddo
         call strain_u(sij,k)
         do n=1,6
         do j=iys,iye
         do i=1,nnx
            fac   = 2.0*(weit(k)*vis_m(i,j,kp1) +
     +                   weitm(k)*vis_m(i,j,k))
            tauij(i,j,n,2) = fac*sij(i,j,n)
         enddo
         enddo
         enddo
         uw_sum = 0.0
         vw_sum = 0.0
         ww_sum = 0.0
         do j=iys,iye
         do i=1,nnx
            gm_xu     = weit(k)*gm(i,j,kp1,1)*gm(i,j,kp1,2) +
     +                  weitm(k)*gm(i,j,k,1)*gm(i,j,k,2)
            gm_xl     = weit(km1)*gm(i,j,k,1)*gm(i,j,k,2) +
     +                  weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,2)
            gm_yu     = weit(k)*gm(i,j,kp1,1)*gm(i,j,kp1,3) +
     +                  weitm(k)*gm(i,j,k,1)*gm(i,j,k,3)
            gm_yl     = weit(km1)*gm(i,j,k,1)*gm(i,j,k,3) +
     +                  weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,3)
            tau_up    = tauij(i,j,1,2)*gm_xu + tauij(i,j,2,2)*gm_yu +
     +                  tauij(i,j,3,2)
            tau_lo    = tauij(i,j,1,1)*gm_xl + tauij(i,j,2,1)*gm_yl +
     +                  tauij(i,j,3,1)
            uw_sum    = uw_sum + tau_up
            r1(i,j,k) = r1(i,j,k) + (tau_up - tau_lo)*dzw_i(k)
            tau_up    = tauij(i,j,2,2)*gm_xu + tauij(i,j,4,2)*gm_yu +
     +                  tauij(i,j,5,2)
            tau_lo    = tauij(i,j,2,1)*gm_xl + tauij(i,j,4,1)*gm_yl +
     +                  tauij(i,j,5,1)
            vw_sum    = vw_sum + tau_up
            r2(i,j,k) = r2(i,j,k) + (tau_up - tau_lo)*dzw_i(k)
            tau_up    = tauij(i,j,3,2)*gm_xu + tauij(i,j,5,2)*gm_yu +
     +                  tauij(i,j,6,2)
            tau_lo    = tauij(i,j,3,1)*gm_xl + tauij(i,j,5,1)*gm_yl +
     +                  tauij(i,j,6,1)
            ww_sum    = ww_sum + tau_up
            r3(i,j,k) = r3(i,j,k) + (tau_up - tau_lo)*dzw_i(k)
         enddo
         enddo
c
c -------- save SGS fluxes for printout
c
      if(istep .eq. 1) then
         tr_tau(k) = 0.0
         do j=iys,iye
         do i=1,nnx
            ufluc    = (u(i,j,kp1) - uxym(kp1))*weit(k) +
     +                 (u(i,j,k) - uxym(k))*weitm(k)
            vfluc    = (v(i,j,kp1) - vxym(kp1))*weit(k) +
     +                 (v(i,j,k) - vxym(k))*weitm(k)
            tr_tau(k) = tr_tau(k) -
     +                 tauij(i,j,3,2)*ufluc - tauij(i,j,5,2)*vfluc
         enddo
         enddo
         uwsb(k)   = -uw_sum*fnxy
         vwsb(k)   = -vw_sum*fnxy
         wwsb(k)   = -ww_sum*fnxy
         tr_tau(k) = tr_tau(k)*fnxy
      endif
c
c ---------- end z loop
c
      enddo
c
c ---------- SGS fluxes tau_12, tau_22, tau_23 that depend on
c            y-derivatives at cell centers
c
      do k=izs,ize
         call strain_12(sij(1,iys,1),k)
         do j=iys,iye
         do i=1,nnx
            fac         = 2.0*gm(i,j,k,1)*vis_m(i,j,k)
c           fntd(i,j,k) = fac*sij(i,j,1) - upfac*vf(i,j,k)*u(i,j,k)
            fntd(i,j,k) = tauij(i,j,2,2) - upfac*vf(i,j,k)*u(i,j,k)
         enddo
         enddo
      enddo
      call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            r1(i,j,k)   = r1(i,j,k) + fntd(i,j,k)
         enddo
         enddo
      enddo
      do k=izs,ize
         call strain_22(sij(1,iys,1),k)
         do j=iys,iye
         do i=1,nnx
            fac         = 2.0*gm(i,j,k,1)*vis_m(i,j,k)
c           fntd(i,j,k) = fac*sij(i,j,1) - upfac*vf(i,j,k)*v(i,j,k)
            fntd(i,j,k) = tauij(i,j,4,2) - upfac*vf(i,j,k)*u(i,j,k)
         enddo
         enddo
      enddo
      call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            r2(i,j,k)   = r2(i,j,k) + fntd(i,j,k)
         enddo
         enddo
      enddo
      do k=izs,ize
         call strain_23(sij(1,iys,1),k)
         do j=iys,iye
         do i=1,nnx
            fac       = 2.0*gm(i,j,k,1)*vis_m(i,j,k)
c           fntd(i,j,k) = fac*sij(i,j,1) - upfac*vf(i,j,k)*w(i,j,k)
            fntd(i,j,k) = tauij(i,j,5,2) - upfac*vf(i,j,k)*u(i,j,k)
         enddo
         enddo
      enddo
      call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
c
c --------- add a constant to r3 to account for hydrostatic pressure
c           and for periodic sidewalls. we have r3/J therefore the
c           constant to be added  = - <J*r3/J>/J
c
      do k=1,nnz
         r3_sum(k) = 0.0
      enddo
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            r3(i,j,k) = r3(i,j,k) + fntd(i,j,k)
            r3_sum(k) = r3_sum(k) + r3(i,j,k)*gm(i,j,k,4)
         enddo
         enddo
         r3_sum(k) = r3_sum(k)*fnxy
      enddo
c
      call mpi_sum_z(r3_sum,i_root,myid,nnz,1)
c
c ------- make sure <r3/J> = 0 and set r3 = 0 at top
c
      do k=izs,ize
         if(k .eq. nnz) then
            do j=iys,iye
            do i=1,nnx
               r3(i,j,k) = 0.0
            enddo
            enddo
         else
            do j=iys,iye
            do i=1,nnx
               r3(i,j,k) = r3(i,j,k) - r3_sum(k)*gm(i,j,k,1)
            enddo
            enddo
         endif
      enddo
c
c --------- leave with (r1,r2,r3)/J
c
      return
      end
      subroutine strain_u(sij_u,k)
c
c ---------- get strain matrix at the upper face at the k-level
c
c        strain definitions
c
c        ux   = uxi  + uzet*gm_x
c        uy   = ueta + uzet*gm_y
c        uz   = uzet*gm_z
c        vx   = vxi  + vzet*gm_x
c        vy   = veta + vzet*gm_y
c        vz   = vzet*gm_z
c        wx   = wxi  + wzet*gm_x
c        wy   = weta + wzet*gm_y
c        wz   = wzet*gm_z
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real sij_u(nnx,iys:iye,6)
c
      km1  = k - 1
      kp1  = k + 1
      w1   = weit(k)
      w1m  = weitm(k)
      w2   = weit(km1)
      w2m  = weitm(km1)
c
c -------- at upper face
c
      do j=iys,iye
      do i=1,nnx
         uzet  = (u(i,j,kp1) - u(i,j,k))*dzu_i(kp1)
         vzet  = (v(i,j,kp1) - v(i,j,k))*dzu_i(kp1)
         wzet  = (w(i,j,kp1) - w(i,j,k))*dzu_i(kp1)
         uxi   = w1*ux(i,j,kp1) + w1m*ux(i,j,k)
         vxi   = w1*vx(i,j,kp1) + w1m*vx(i,j,k)
         wxi   = w1*wx(i,j,kp1) + w1m*wx(i,j,k)
         ueta  = w1*uy(i,j,kp1) + w1m*uy(i,j,k)
         veta  = w1*vy(i,j,kp1) + w1m*vy(i,j,k)
         weta  = w1*wy(i,j,kp1) + w1m*wy(i,j,k)
         gm_x  = w1*gm(i,j,kp1,2) + w1m*gm(i,j,k,2)
         gm_y  = w1*gm(i,j,kp1,3) + w1m*gm(i,j,k,3)
         gm_z  = w1*gm(i,j,kp1,4) +  w1m*gm(i,j,k,4)
         sij_u(i,j,1) = uxi  + uzet*gm_x
         sij_u(i,j,2) = 0.5*(ueta + uzet*gm_y + vxi  + vzet*gm_x)
         sij_u(i,j,3) = 0.5*(uzet*gm_z + wxi  + wzet*gm_x)
         sij_u(i,j,4) = veta + vzet*gm_y
         sij_u(i,j,5) = 0.5*(vzet*gm_z + weta + wzet*gm_y)
         sij_u(i,j,6) = wzet*gm_z
      enddo
      enddo
c
      return
      end
      subroutine strain_c_dns(sij_c,k)
c
c ---------- get strain matrix at the cell center
c ---------- for DNS just don't treak k = 1 any different
c
c        strain definitions
c
c        ux   = uxi  + uzet*gm_x
c        uy   = ueta + uzet*gm_y
c        uz   = uzet*gm_z
c        vx   = vxi  + vzet*gm_x
c        vy   = veta + vzet*gm_y
c        vz   = vzet*gm_z
c        wx   = wxi  + wzet*gm_x
c        wy   = weta + wzet*gm_y
c        wz   = wzet*gm_z
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real sij_c(nnx,iys:iye,6)
c
      km1  = k - 1
      kp1  = k + 1
      w1   = weit(k)
      w1m  = weitm(k)
      w2   = weit(km1)
      w2m  = weitm(km1)
c
!      if(k == 1) then
!c
!c ------- vertical gradients at upper face
!Cnew      now computed at cell centers using surface bcs
!c
!      do j=iys,iye
!      do i=1,nnx
!         uzet  = (u(i,j,kp1) - u(i,j,k))*dzu_i(kp1)
!         vzet  = (v(i,j,kp1) - v(i,j,k))*dzu_i(kp1)
!         wzet  = (w(i,j,kp1) - w(i,j,k))*dzu_i(kp1)
!         gm_x  = gm(i,j,k,2)
!         gm_y  = gm(i,j,k,3)
!         gm_z  = gm(i,j,k,4)
!         sij_c(i,j,1) = ux(i,j,k)  + uzet*gm_x
!         sij_c(i,j,2) = 0.5*(uy(i,j,k) + uzet*gm_y +
!     +                       vx(i,j,k)  + vzet*gm_x)
!         sij_c(i,j,3) = 0.5*(uzet*gm_z + wx(i,j,k)  + wzet*gm_x)
!         sij_c(i,j,4) = vy(i,j,k) + vzet*gm_y
!         sij_c(i,j,5) = 0.5*(vzet*gm_z + wy(i,j,k) + wzet*gm_y)
!         sij_c(i,j,6) = wzet*gm_z
!      enddo
!      enddo
!
!      else
c
c -------- at cell center
c
      do j=iys,iye
      do i=1,nnx
         uup   = w1*u(i,j,kp1) + w1m*u(i,j,k)
         ulo   = w2*u(i,j,k) + w2m*u(i,j,km1)
         uzet  = (uup - ulo)*dzw_i(k)
         vup   = w1*v(i,j,kp1) + w1m*v(i,j,k)
         vlo   = w2*v(i,j,k) + w2m*v(i,j,km1)
         vzet  = (vup - vlo)*dzw_i(k)
         wup   = w1*w(i,j,kp1) + w1m*w(i,j,k)
         wlo   = w2*w(i,j,k) + w2m*w(i,j,km1)
         wzet  = (wup - wlo)*dzw_i(k)
         gm_x  = gm(i,j,k,2)
         gm_y  = gm(i,j,k,3)
         gm_z  = gm(i,j,k,4)
         sij_c(i,j,1) = ux(i,j,k)  + uzet*gm_x
         sij_c(i,j,2) = 0.5*(uy(i,j,k) + uzet*gm_y +
     +                       vx(i,j,k)  + vzet*gm_x)
         sij_c(i,j,3) = 0.5*(uzet*gm_z + wx(i,j,k)  + wzet*gm_x)
         sij_c(i,j,4) = vy(i,j,k) + vzet*gm_y
         sij_c(i,j,5) = 0.5*(vzet*gm_z + wy(i,j,k) + wzet*gm_y)
         sij_c(i,j,6) = wzet*gm_z
      enddo
      enddo
c
!      endif
c
      return
      end
      subroutine strain_c(sij_c,k)
c
c ---------- get strain matrix at the cell center
c            special if at first level off the wall
c
c        strain definitions
c
c        ux   = uxi  + uzet*gm_x
c        uy   = ueta + uzet*gm_y
c        uz   = uzet*gm_z
c        vx   = vxi  + vzet*gm_x
c        vy   = veta + vzet*gm_y
c        vz   = vzet*gm_z
c        wx   = wxi  + wzet*gm_x
c        wy   = weta + wzet*gm_y
c        wz   = wzet*gm_z
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real sij_c(nnx,iys:iye,6)
c
      km1  = k - 1
      kp1  = k + 1
      w1   = weit(k)
      w1m  = weitm(k)
      w2   = weit(km1)
      w2m  = weitm(km1)
c
      if(k == 1) then
c
c ------- vertical gradients at upper face
Cnew      now computed at cell centers using surface bcs
c
      do j=iys,iye
      do i=1,nnx
         uzet  = (u(i,j,kp1) - u(i,j,k))*dzu_i(kp1)
         vzet  = (v(i,j,kp1) - v(i,j,k))*dzu_i(kp1)
         wzet  = (w(i,j,kp1) - w(i,j,k))*dzu_i(kp1)
         gm_x  = gm(i,j,k,2)
         gm_y  = gm(i,j,k,3)
         gm_z  = gm(i,j,k,4)
         sij_c(i,j,1) = ux(i,j,k)  + uzet*gm_x
         sij_c(i,j,2) = 0.5*(uy(i,j,k) + uzet*gm_y +
     +                       vx(i,j,k)  + vzet*gm_x)
         sij_c(i,j,3) = 0.5*(uzet*gm_z + wx(i,j,k)  + wzet*gm_x)
         sij_c(i,j,4) = vy(i,j,k) + vzet*gm_y
         sij_c(i,j,5) = 0.5*(vzet*gm_z + wy(i,j,k) + wzet*gm_y)
         sij_c(i,j,6) = wzet*gm_z
      enddo
      enddo

      else
c
c -------- at cell center
c
      do j=iys,iye
      do i=1,nnx
         uup   = w1*u(i,j,kp1) + w1m*u(i,j,k)
         ulo   = w2*u(i,j,k) + w2m*u(i,j,km1)
         uzet  = (uup - ulo)*dzw_i(k)
         vup   = w1*v(i,j,kp1) + w1m*v(i,j,k)
         vlo   = w2*v(i,j,k) + w2m*v(i,j,km1)
         vzet  = (vup - vlo)*dzw_i(k)
         wup   = w1*w(i,j,kp1) + w1m*w(i,j,k)
         wlo   = w2*w(i,j,k) + w2m*w(i,j,km1)
         wzet  = (wup - wlo)*dzw_i(k)
         gm_x  = gm(i,j,k,2)
         gm_y  = gm(i,j,k,3)
         gm_z  = gm(i,j,k,4)
         sij_c(i,j,1) = ux(i,j,k)  + uzet*gm_x
         sij_c(i,j,2) = 0.5*(uy(i,j,k) + uzet*gm_y +
     +                       vx(i,j,k)  + vzet*gm_x)
         sij_c(i,j,3) = 0.5*(uzet*gm_z + wx(i,j,k)  + wzet*gm_x)
         sij_c(i,j,4) = vy(i,j,k) + vzet*gm_y
         sij_c(i,j,5) = 0.5*(vzet*gm_z + wy(i,j,k) + wzet*gm_y)
         sij_c(i,j,6) = wzet*gm_z
      enddo
      enddo
c
      endif
c
      return
      end
      subroutine strain_12_dns(sij_y,k)
c
c ---------- get strain s12
c            for special evaluation of eta derivatives
c            using matrix transpose
c
c        strain definitions
c
c        ux   = uxi  + uzet*gm_x
c        uy   = ueta + uzet*gm_y
c        uz   = uzet*gm_z
c        vx   = vxi  + vzet*gm_x
c        vy   = veta + vzet*gm_y
c        vz   = vzet*gm_z
c        wx   = wxi  + wzet*gm_x
c        wy   = weta + wzet*gm_y
c        wz   = wzet*gm_z
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real sij_y(nnx,iys:iye)
c
      km1  = k - 1
      kp1  = k + 1
      w1   = weit(k)
      w1m  = weitm(k)
      w2   = weit(km1)
      w2m  = weitm(km1)
c
!      if(k == 1) then
!c
!c ------ at upper face
!c
!      do j=iys,iye
!      do i=1,nnx
!         uzet       = (u(i,j,kp1) - u(i,j,k))*dzu_i(kp1)
!         vzet       = (v(i,j,kp1) - v(i,j,k))*dzu_i(kp1)
!         vxi        = w1*vx(i,j,kp1) + w1m*vx(i,j,k)
!         ueta       = w1*uy(i,j,kp1) + w1m*uy(i,j,k)
!         gm_x       = w1*gm(i,j,kp1,2) + w1m*gm(i,j,k,2)
!         gm_y       = w1*gm(i,j,kp1,3) + w1m*gm(i,j,k,3)
!         sij_y(i,j) = 0.5*(ueta + uzet*gm_y + vxi  + vzet*gm_x)
!      enddo
!      enddo
!c
!      else
c
c -------- at cell center
c
      do j=iys,iye
      do i=1,nnx
         uup   = w1*u(i,j,kp1) + w1m*u(i,j,k)
         ulo   = w2*u(i,j,k) + w2m*u(i,j,km1)
         uzet  = (uup - ulo)*dzw_i(k)
         vup   = w1*v(i,j,kp1) + w1m*v(i,j,k)
         vlo   = w2*v(i,j,k) + w2m*v(i,j,km1)
         vzet  = (vup - vlo)*dzw_i(k)
         sij_y(i,j) = 0.5*(uy(i,j,k) + uzet*gm(i,j,k,3) +
     +                     vx(i,j,k) + vzet*gm(i,j,k,2))
      enddo
      enddo
c
!      endif
c
      return
      end
      subroutine strain_12(sij_y,k)
c
c ---------- get strain s12
c            for special evaluation of eta derivatives
c            using matrix transpose
c
c        strain definitions
c
c        ux   = uxi  + uzet*gm_x
c        uy   = ueta + uzet*gm_y
c        uz   = uzet*gm_z
c        vx   = vxi  + vzet*gm_x
c        vy   = veta + vzet*gm_y
c        vz   = vzet*gm_z
c        wx   = wxi  + wzet*gm_x
c        wy   = weta + wzet*gm_y
c        wz   = wzet*gm_z
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real sij_y(nnx,iys:iye)
c
      km1  = k - 1
      kp1  = k + 1
      w1   = weit(k)
      w1m  = weitm(k)
      w2   = weit(km1)
      w2m  = weitm(km1)
c
      if(k == 1) then
c
c ------ at upper face
c
      do j=iys,iye
      do i=1,nnx
         uzet       = (u(i,j,kp1) - u(i,j,k))*dzu_i(kp1)
         vzet       = (v(i,j,kp1) - v(i,j,k))*dzu_i(kp1)
         vxi        = w1*vx(i,j,kp1) + w1m*vx(i,j,k)
         ueta       = w1*uy(i,j,kp1) + w1m*uy(i,j,k)
         gm_x       = w1*gm(i,j,kp1,2) + w1m*gm(i,j,k,2)
         gm_y       = w1*gm(i,j,kp1,3) + w1m*gm(i,j,k,3)
         sij_y(i,j) = 0.5*(ueta + uzet*gm_y + vxi  + vzet*gm_x)
      enddo
      enddo
c
      else
c
c -------- at cell center
c
      do j=iys,iye
      do i=1,nnx
         uup   = w1*u(i,j,kp1) + w1m*u(i,j,k)
         ulo   = w2*u(i,j,k) + w2m*u(i,j,km1)
         uzet  = (uup - ulo)*dzw_i(k)
         vup   = w1*v(i,j,kp1) + w1m*v(i,j,k)
         vlo   = w2*v(i,j,k) + w2m*v(i,j,km1)
         vzet  = (vup - vlo)*dzw_i(k)
         sij_y(i,j) = 0.5*(uy(i,j,k) + uzet*gm(i,j,k,3) +
     +                     vx(i,j,k) + vzet*gm(i,j,k,2))
      enddo
      enddo
c
      endif
c
      return
      end
      subroutine strain_22_dns(sij_y,k)
c
c ---------- get strain s22 at cell center
c            for special evaluation of eta derivatives
c            using matrix transpose
c
c        strain definitions
c
c        ux   = uxi  + uzet*gm_x
c        uy   = ueta + uzet*gm_y
c        uz   = uzet*gm_z
c        vx   = vxi  + vzet*gm_x
c        vy   = veta + vzet*gm_y
c        vz   = vzet*gm_z
c        wx   = wxi  + wzet*gm_x
c        wy   = weta + wzet*gm_y
c        wz   = wzet*gm_z
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real sij_y(nnx,iys:iye)
c
      km1  = k - 1
      kp1  = k + 1
      w1   = weit(k)
      w1m  = weitm(k)
      w2   = weit(km1)
      w2m  = weitm(km1)
c
!      if(k == 1) then
!c
!c -------at upper face
!c
!      do j=iys,iye
!      do i=1,nnx
!         vzet       = (v(i,j,kp1) - v(i,j,k))*dzu_i(kp1)
!         veta       = w1*vy(i,j,kp1) + w1m*vy(i,j,k)
!         gm_y       = w1*gm(i,j,kp1,3) + w1m*gm(i,j,k,3)
!         sij_y(i,j) = veta + vzet*gm_y
!      enddo
!      enddo
!      else
c
c -------- at cell center
c
      do j=iys,iye
      do i=1,nnx
         vup   = w1*v(i,j,kp1) + w1m*v(i,j,k)
         vlo   = w2*v(i,j,k) + w2m*v(i,j,km1)
         vzet  = (vup - vlo)*dzw_i(k)
         sij_y(i,j) = vy(i,j,k) + vzet*gm(i,j,k,3)
      enddo
      enddo
c
!      endif
c
      return
      end
      subroutine strain_22(sij_y,k)
c
c ---------- get strain s22 at cell center
c            for special evaluation of eta derivatives
c            using matrix transpose
c
c        strain definitions
c
c        ux   = uxi  + uzet*gm_x
c        uy   = ueta + uzet*gm_y
c        uz   = uzet*gm_z
c        vx   = vxi  + vzet*gm_x
c        vy   = veta + vzet*gm_y
c        vz   = vzet*gm_z
c        wx   = wxi  + wzet*gm_x
c        wy   = weta + wzet*gm_y
c        wz   = wzet*gm_z
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real sij_y(nnx,iys:iye)
c
      km1  = k - 1
      kp1  = k + 1
      w1   = weit(k)
      w1m  = weitm(k)
      w2   = weit(km1)
      w2m  = weitm(km1)
c
      if(k == 1) then
c
c -------at upper face
c
      do j=iys,iye
      do i=1,nnx
         vzet       = (v(i,j,kp1) - v(i,j,k))*dzu_i(kp1)
         veta       = w1*vy(i,j,kp1) + w1m*vy(i,j,k)
         gm_y       = w1*gm(i,j,kp1,3) + w1m*gm(i,j,k,3)
         sij_y(i,j) = veta + vzet*gm_y
      enddo
      enddo
      else
c
c -------- at cell center
c
      do j=iys,iye
      do i=1,nnx
         vup   = w1*v(i,j,kp1) + w1m*v(i,j,k)
         vlo   = w2*v(i,j,k) + w2m*v(i,j,km1)
         vzet  = (vup - vlo)*dzw_i(k)
         sij_y(i,j) = vy(i,j,k) + vzet*gm(i,j,k,3)
      enddo
      enddo
c
      endif
c
      return
      end
      subroutine strain_23_dns(sij_y,k)
c
c ---------- get strain s23 at cell center
c            for special evaluation of eta derivatives
c            using matrix transpose
c
c        strain definitions
c
c        ux   = uxi  + uzet*gm_x
c        uy   = ueta + uzet*gm_y
c        uz   = uzet*gm_z
c        vx   = vxi  + vzet*gm_x
c        vy   = veta + vzet*gm_y
c        vz   = vzet*gm_z
c        wx   = wxi  + wzet*gm_x
c        wy   = weta + wzet*gm_y
c        wz   = wzet*gm_z
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real sij_y(nnx,iys:iye)
c
      km1  = k - 1
      kp1  = k + 1
      w1   = weit(k)
      w1m  = weitm(k)
      w2   = weit(km1)
      w2m  = weitm(km1)
c
!      if(k == 1) then
!c
!c -------- at upper face
!c
!      do j=iys,iye
!      do i=1,nnx
!         vzet       = (v(i,j,kp1) - v(i,j,k))*dzu_i(kp1)
!         wzet       = (w(i,j,kp1) - w(i,j,k))*dzu_i(kp1)
!         weta       = w1*wy(i,j,kp1) + w1m*wy(i,j,k)
!         gm_y       = w1*gm(i,j,kp1,3) + w1m*gm(i,j,k,3)
!         gm_z       = w1*gm(i,j,kp1,4) +  w1m*gm(i,j,k,4)
!         sij_y(i,j) = 0.5*(vzet*gm_z + weta + wzet*gm_y)
!      enddo
!      enddo
!      else
c
c -------- at cell center
c
      do j=iys,iye
      do i=1,nnx
         vup   = w1*v(i,j,kp1) + w1m*v(i,j,k)
         vlo   = w2*v(i,j,k) + w2m*v(i,j,km1)
         vzet  = (vup - vlo)*dzw_i(k)
         wup   = w1*w(i,j,kp1) + w1m*w(i,j,k)
         wlo   = w2*w(i,j,k) + w2m*w(i,j,km1)
         wzet  = (wup - wlo)*dzw_i(k)
         sij_y(i,j) = 0.5*(vzet*gm(i,j,k,4) +
     +                     wy(i,j,k) + wzet*gm(i,j,k,3))
      enddo
      enddo
c
!      endif
c
      return
      end
      subroutine strain_23(sij_y,k)
c
c ---------- get strain s23 at cell center
c            for special evaluation of eta derivatives
c            using matrix transpose
c
c        strain definitions
c
c        ux   = uxi  + uzet*gm_x
c        uy   = ueta + uzet*gm_y
c        uz   = uzet*gm_z
c        vx   = vxi  + vzet*gm_x
c        vy   = veta + vzet*gm_y
c        vz   = vzet*gm_z
c        wx   = wxi  + wzet*gm_x
c        wy   = weta + wzet*gm_y
c        wz   = wzet*gm_z
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real sij_y(nnx,iys:iye)
c
      km1  = k - 1
      kp1  = k + 1
      w1   = weit(k)
      w1m  = weitm(k)
      w2   = weit(km1)
      w2m  = weitm(km1)
c
      if(k == 1) then
c
c -------- at upper face
c
      do j=iys,iye
      do i=1,nnx
         vzet       = (v(i,j,kp1) - v(i,j,k))*dzu_i(kp1)
         wzet       = (w(i,j,kp1) - w(i,j,k))*dzu_i(kp1)
         weta       = w1*wy(i,j,kp1) + w1m*wy(i,j,k)
         gm_y       = w1*gm(i,j,kp1,3) + w1m*gm(i,j,k,3)
         gm_z       = w1*gm(i,j,kp1,4) +  w1m*gm(i,j,k,4)
         sij_y(i,j) = 0.5*(vzet*gm_z + weta + wzet*gm_y)
      enddo
      enddo
      else
c
c -------- at cell center
c
      do j=iys,iye
      do i=1,nnx
         vup   = w1*v(i,j,kp1) + w1m*v(i,j,k)
         vlo   = w2*v(i,j,k) + w2m*v(i,j,km1)
         vzet  = (vup - vlo)*dzw_i(k)
         wup   = w1*w(i,j,kp1) + w1m*w(i,j,k)
         wlo   = w2*w(i,j,k) + w2m*w(i,j,km1)
         wzet  = (wup - wlo)*dzw_i(k)
         sij_y(i,j) = 0.5*(vzet*gm(i,j,k,4) +
     +                     wy(i,j,k) + wzet*gm(i,j,k,3))
      enddo
      enddo
c
      endif
c
      return
      end
      subroutine rhs_scl_dns(istep,iscl)
c
c ------ get right hand side of scalar equation (iscl)
c        monotone scalar fluxes only in z
c        for pencil size (nnx, iys:iye, izs:ize)
c        care is taken so that if monotone is on then
c        conservative horizontal flux form is used!
c        uses flux conservation form for advection always
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      real fnt1(nnx,iys:iye), fntd(nnx,iys:iye,izs:ize)
      real tx(nnx,iys:iye,izs-1:ize+1), ty(nnx,iys:iye,izs-1:ize+1)
      real flux_u(nnx,iys:iye), flux_l(nnx,iys:iye)
      real tflx(nnx,iys:iye,3,2), scl_g(nnx,iys:iye,3)
      real sfc_flx
c
c --------- set sign for ocean simulations that use monotone
c
      sgn = 1.0
      if(iocean .eq. 1) sgn = -1.0
      if(iupwnd .eq. 1) then
         upfac = 1.0
      else
         upfac = 0.5
      endif
      sfc_flx = 0.0
c
c --------- get x & y derivatives for [izs-1:ize+1]
c
      do k=izs-1,ize+1
         do j=iys,iye
         do i=1,nnx
            tx(i,j,k) = t(i,j,iscl,k)
            ty(i,j,k) = t(i,j,iscl,k)
         enddo
         enddo
         call xderivp(tx(1,iys,k),trigx(1,1),xk(1),nnx,iys,iye)
      enddo
      call yd_mpi(ty(1,iys,izs-1),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
c
c -------- SGS term tflx(.,.,1) + resolved flux at cell center
c
      do k=izs,ize
         kp1 = k + 1
         km1 = k - 1
!         if(k .eq. 1) then
!            call grads_u(scl_g,k,tx,ty,iscl)
!            do j=iys,iye
!            do i=1,nnx
!               gm_u  = weit(k)*gm(i,j,kp1,1) + weitm(k)*gm(i,j,k,1)
!               gm_l  = weit(km1)*gm(i,j,k,1) + weitm(km1)*gm(i,j,km1,1)
!               fac   = weit(k)*vis_s(i,j,kp1) + weitm(k)*vis_s(i,j,k)
!               tflx(i,j,1,1) = 0.5*(fac*scl_g(i,j,1)*gm_u -
!     +                               t1tbc(i,j,iscl)*gm_l)
!               tflx(i,j,2,1) = 0.5*(fac*scl_g(i,j,2)*gm_u -
!     +                               t2tbc(i,j,iscl)*gm_l)
!               fnt1(i,j)     = tflx(i,j,1,1) -
!     +                         upfac*uf(i,j,k)*t(i,j,iscl,k)
!               fntd(i,j,k)   = tflx(i,j,2,1) -
!     +                         upfac*vf(i,j,k)*t(i,j,iscl,k)
!            enddo
!            enddo
!         else
            !call grads_c(scl_g,k,tx,ty,iscl)
            call grads_c_dns(scl_g,k,tx,ty,iscl)
            do j=iys,iye
            do i=1,nnx
               tflx(i,j,1,1) = vis_s(i,j,k)*scl_g(i,j,1)*gm(i,j,k,1)
               tflx(i,j,2,1) = vis_s(i,j,k)*scl_g(i,j,2)*gm(i,j,k,1)
               fnt1(i,j)     = tflx(i,j,1,1) -
     +                         upfac*uf(i,j,k)*t(i,j,iscl,k)
               fntd(i,j,k)   = tflx(i,j,2,1) -
     +                         upfac*vf(i,j,k)*t(i,j,iscl,k)
            enddo
            enddo
!         endif
         call xderivp(fnt1(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         do j=iys,iye
         do i=1,nnx
            r4(i,j,iscl,k) = fnt1(i,j)
         enddo
         enddo
c
c -------- save SGS fluxes for printout
c
         if(istep .eq. 1) then
            utsb(k,iscl) = 0.0
            vtsb(k,iscl) = 0.0
            do j=iys,iye
            do i=1,nnx
               utsb(k,iscl) = utsb(k,iscl) - tflx(i,j,1,1)
               vtsb(k,iscl) = vtsb(k,iscl) - tflx(i,j,2,1)
            enddo
            enddo
            utsb(k,iscl) = utsb(k,iscl)*fnxy
            vtsb(k,iscl) = vtsb(k,iscl)*fnxy
         endif
      enddo
c
c ---------- SGS fluxes tflx(.,.,2) that depend on
c            y-derivative at cell centers
c
      call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            r4(i,j,iscl,k) = r4(i,j,iscl,k) + fntd(i,j,k)
         enddo
         enddo
      enddo
c
c ------- resolved advection for skew symmetric form
c         for no monotone
c
      if(iupwnd .ne. 1) then
c
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            r4(i,j,iscl,k) = r4(i,j,iscl,k) -
     +            upfac*(uf(i,j,k)*tx(i,j,k) + vf(i,j,k)*ty(i,j,k))
         enddo
         enddo
      enddo
      endif ! end of iupwnd check for horizontal advection
c
c ----------- SGS z-contributions
c
!      if(izs .eq. 1) then
!         do j=iys,iye
!         do i=1,nnx
!            tflx(i,j,1,2) = -t1tbc(i,j,iscl)
!            tflx(i,j,2,2) = -t2tbc(i,j,iscl)
!            tflx(i,j,3,2) = -t3tbc(i,j,iscl)
!         enddo
!         enddo
!      else
         call grads_u(scl_g,izs-1,tx,ty,iscl)
         do n=1,3
         do j=iys,iye
         do i=1,nnx
            fac   = weit(izs-1)*vis_s(i,j,izs) +
     +              weitm(izs-1)*vis_s(i,j,izs-1)
            tflx(i,j,n,2) = fac*scl_g(i,j,n)
         enddo
         enddo
         enddo
!      endif
      do k=izs,ize
         kp1 = k + 1
         km1 = k - 1
         do n=1,3
         do j=iys,iye
         do i=1,nnx
            tflx(i,j,n,1) = tflx(i,j,n,2)
         enddo
         enddo
         enddo
         call grads_u(scl_g,k,tx,ty,iscl)
         do n=1,3
         do j=iys,iye
         do i=1,nnx
            fac   = weit(k)*vis_s(i,j,kp1) + weitm(k)*vis_s(i,j,k)
            tflx(i,j,n,2) = fac*scl_g(i,j,n)
         enddo
         enddo
         enddo
         do j=iys,iye
         do i=1,nnx
            gm_xu      = weit(k)*gm(i,j,kp1,1)*gm(i,j,kp1,2) +
     +                   weitm(k)*gm(i,j,k,1)*gm(i,j,k,2)
            gm_xl      = weit(km1)*gm(i,j,k,1)*gm(i,j,k,2) +
     +                   weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,2)
            gm_yu      = weit(k)*gm(i,j,kp1,1)*gm(i,j,kp1,3) +
     +                   weitm(k)*gm(i,j,k,1)*gm(i,j,k,3)
            gm_yl      = weit(km1)*gm(i,j,k,1)*gm(i,j,k,3) +
     +                   weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,3)
            tflx_up    = tflx(i,j,1,2)*gm_xu + tflx(i,j,2,2)*gm_yu +
     +                   tflx(i,j,3,2)
            tflx_lo    = tflx(i,j,1,1)*gm_xl + tflx(i,j,2,1)*gm_yl +
     +                   tflx(i,j,3,1)
            r4(i,j,iscl,k) = r4(i,j,iscl,k) +
     +                       (tflx_up - tflx_lo)*dzw_i(k)
            if (k==1) then
            sfc_flx = sfc_flx + tflx_lo
            end if
         enddo
         enddo
c
c -------- save SGS fluxes for printout, gather sums on exit
c
         if(istep .eq. 1) then
            wtsb(k,iscl) = 0.0
            do j=iys,iye
            do i=1,nnx
               wtsb(k,iscl) = wtsb(k,iscl) - tflx(i,j,3,2)
            enddo
            enddo
            wtsb(k,iscl) = wtsb(k,iscl)*fnxy
         endif

      !As an aside, compute the values uwsfc and vwsfc for output:
      if (k==1 .and. isfc(1)==1) then
         call mpi_sum_xy(sfc_flx,myid,iss,ise,1)
         wtsfc(iscl) = sfc_flx
      end if

c
c ---------- end z loop
c
      enddo
c
c ---------- resolved vertical advection
c
      if(iupwnd .ne. 1) then
c
c --------- skew symmetric advective form for
c           vertical flux = 0.5(wf*dt/dz + d/dz(wf*t))
c           no grid speed
c
      do k=izs,ize
         km1 = k - 1
         kp1 = k + 1
         do j=iys,iye
         do i=1,nnx
            tzm     = (t(i,j,iscl,k) - t(i,j,iscl,km1))*dzu_i(k)
            tz      = (t(i,j,iscl,kp1) - t(i,j,iscl,k))*dzu_i(kp1)
            theta_u = weitm(k)*t(i,j,iscl,k) +
     +                weit(k)*t(i,j,iscl,kp1)
            theta_l = weit(km1)*t(i,j,iscl,k) +
     +                weitm(km1)*t(i,j,iscl,km1)
            t_adv   = 0.5*(wf(i,j,k)*theta_u -
     +                     wf(i,j,km1)*theta_l)*dzw_i(k) +
     +                0.25*(wf(i,j,k)*tz + wf(i,j,km1)*tzm)
            r4(i,j,iscl,k) = r4(i,j,iscl,k) - t_adv
         enddo
         enddo
      enddo
c
      else
c
c ----------- z-direction special with monotone
c
      do k=izs,ize
         km2 = k - 2
         km1 = k - 1
         kp1 = k + 1
         kp2 = k + 2
         if(k .eq. 1) then
              do j=iys,iye
              do i=1,nnx
                 wflx_l      = wf(i,j,km1) - gm(i,j,km1,5)
                 wflx_u      = wf(i,j,k)   - gm(i,j,k,5)
                 flux_l(i,j) = sgn*0.5*wflx_l*
     +                         (t(i,j,iscl,km1)+t(i,j,iscl,k))
                 flux_u(i,j) =
     +           amax1(sgn*wflx_u,0.)*(t(i,j,iscl,k) +
     +           rlim(t(i,j,iscl,kp1),t(i,j,iscl,k),
     +                t(i,j,iscl,km1))) +
     +           amin1(sgn*wflx_u,0.)*(t(i,j,iscl,kp1) +
     +           rlim(t(i,j,iscl,k),t(i,j,iscl,kp1),
     +                t(i,j,iscl,kp2)))
              enddo
              enddo
         else if(k .eq. nnz) then
              do j=iys,iye
              do i=1,nnx
                 wflx_l      = wf(i,j,km1) - gm(i,j,km1,5)
                 wflx_u      = wf(i,j,k) - gm(i,j,k,5)
                 flux_u(i,j) = sgn*0.5*wflx_u*
     +                        (t(i,j,iscl,kp1)+t(i,j,iscl,k))
                 flux_l(i,j) =
     +           amax1(sgn*wflx_l,0.)*(t(i,j,iscl,km1) +
     +           rlim(t(i,j,iscl,k),t(i,j,iscl,km1),
     +                t(i,j,iscl,km2))) +
     +           amin1(sgn*wflx_l,0.)*(t(i,j,iscl,k) +
     +           rlim(t(i,j,iscl,km1),t(i,j,iscl,k),
     +                t(i,j,iscl,kp1)))
              enddo
              enddo
         else
              do j=iys,iye
              do i=1,nnx
                 wflx_l      = wf(i,j,km1) - gm(i,j,km1,5)
                 wflx_u      = wf(i,j,k) - gm(i,j,k,5)
                 flux_u(i,j) =
     +           amax1(sgn*wflx_u,0.)*(t(i,j,iscl,k) +
     +           rlim(t(i,j,iscl,kp1),t(i,j,iscl,k),
     +                t(i,j,iscl,km1))) +
     +           amin1(sgn*wflx_u,0.)*(t(i,j,iscl,kp1) +
     +           rlim(t(i,j,iscl,k),t(i,j,iscl,kp1),
     +                t(i,j,iscl,kp2)))
                 flux_l(i,j) =
     +           amax1(sgn*wflx_l,0.)*(t(i,j,iscl,km1) +
     +           rlim(t(i,j,iscl,k),t(i,j,iscl,km1),
     +                t(i,j,iscl,km2))) +
     +           amin1(sgn*wflx_l,0.)*(t(i,j,iscl,k) +
     +           rlim(t(i,j,iscl,km1),t(i,j,iscl,k),
     +                t(i,j,iscl,kp1)))
              enddo
              enddo
         endif
c
c ---------- sum vertical monotone flux
c
         do j=iys,iye
         do i=1,nnx
            r4(i,j,iscl,k) = r4(i,j,iscl,k)
     +          - sgn*(flux_u(i,j) - flux_l(i,j))*dzw_i(k)
         enddo
         enddo
c
c ---------- end z loop
c
      enddo
c
c -------- end monotone if block
c
      endif
c
      return
      end
      subroutine rhs_scl(istep,iscl)
c
c ------ get right hand side of scalar equation (iscl)
c        monotone scalar fluxes only in z
c        for pencil size (nnx, iys:iye, izs:ize)
c        care is taken so that if monotone is on then
c        conservative horizontal flux form is used!
c        uses flux conservation form for advection always
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      real fnt1(nnx,iys:iye), fntd(nnx,iys:iye,izs:ize)
      real tx(nnx,iys:iye,izs-1:ize+1), ty(nnx,iys:iye,izs-1:ize+1)
      real flux_u(nnx,iys:iye), flux_l(nnx,iys:iye)
      real tflx(nnx,iys:iye,3,2), scl_g(nnx,iys:iye,3)
c
c --------- set sign for ocean simulations that use monotone
c
      sgn = 1.0
      if(iocean .eq. 1) sgn = -1.0
      if(iupwnd .eq. 1) then
         upfac = 1.0
      else
         upfac = 0.5
      endif
c
c --------- get x & y derivatives for [izs-1:ize+1]
c
      do k=izs-1,ize+1
         do j=iys,iye
         do i=1,nnx
            tx(i,j,k) = t(i,j,iscl,k)
            ty(i,j,k) = t(i,j,iscl,k)
         enddo
         enddo
         call xderivp(tx(1,iys,k),trigx(1,1),xk(1),nnx,iys,iye)
      enddo
      call yd_mpi(ty(1,iys,izs-1),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
c
c -------- SGS term tflx(.,.,1) + resolved flux at cell center
c          special treatment of first grid point off the wall,
c          average upper cell and boundary value to get cell
c          center value.
c
      do k=izs,ize
         kp1 = k + 1
         km1 = k - 1
         if(k .eq. 1) then
            call grads_u(scl_g,k,tx,ty,iscl)
            do j=iys,iye
            do i=1,nnx
               gm_u  = weit(k)*gm(i,j,kp1,1) + weitm(k)*gm(i,j,k,1)
               gm_l  = weit(km1)*gm(i,j,k,1) + weitm(km1)*gm(i,j,km1,1)
               fac   = weit(k)*vis_s(i,j,kp1) + weitm(k)*vis_s(i,j,k)
               tflx(i,j,1,1) = 0.5*(fac*scl_g(i,j,1)*gm_u -
     +                               t1tbc(i,j,iscl)*gm_l)
               tflx(i,j,2,1) = 0.5*(fac*scl_g(i,j,2)*gm_u -
     +                               t2tbc(i,j,iscl)*gm_l)
               fnt1(i,j)     = tflx(i,j,1,1) -
     +                         upfac*uf(i,j,k)*t(i,j,iscl,k)
               fntd(i,j,k)   = tflx(i,j,2,1) -
     +                         upfac*vf(i,j,k)*t(i,j,iscl,k)
            enddo
            enddo
         else
            call grads_c(scl_g,k,tx,ty,iscl)
            do j=iys,iye
            do i=1,nnx
               tflx(i,j,1,1) = vis_s(i,j,k)*scl_g(i,j,1)*gm(i,j,k,1)
               tflx(i,j,2,1) = vis_s(i,j,k)*scl_g(i,j,2)*gm(i,j,k,1)
               fnt1(i,j)     = tflx(i,j,1,1) -
     +                         upfac*uf(i,j,k)*t(i,j,iscl,k)
               fntd(i,j,k)   = tflx(i,j,2,1) -
     +                         upfac*vf(i,j,k)*t(i,j,iscl,k)
            enddo
            enddo
         endif
         call xderivp(fnt1(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         do j=iys,iye
         do i=1,nnx
            r4(i,j,iscl,k) = fnt1(i,j)
         enddo
         enddo
c
c -------- save SGS fluxes for printout
c
         if(istep .eq. 1) then
            utsb(k,iscl) = 0.0
            vtsb(k,iscl) = 0.0
            do j=iys,iye
            do i=1,nnx
               utsb(k,iscl) = utsb(k,iscl) - tflx(i,j,1,1)
               vtsb(k,iscl) = vtsb(k,iscl) - tflx(i,j,2,1)
            enddo
            enddo
            utsb(k,iscl) = utsb(k,iscl)*fnxy
            vtsb(k,iscl) = vtsb(k,iscl)*fnxy
         endif
      enddo
c
c ---------- SGS fluxes tflx(.,.,2) that depend on
c            y-derivative at cell centers
c
      call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            r4(i,j,iscl,k) = r4(i,j,iscl,k) + fntd(i,j,k)
         enddo
         enddo
      enddo
c
c ------- resolved advection for skew symmetric form
c         for no monotone
c
      if(iupwnd .ne. 1) then
c
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            r4(i,j,iscl,k) = r4(i,j,iscl,k) -
     +            upfac*(uf(i,j,k)*tx(i,j,k) + vf(i,j,k)*ty(i,j,k))
         enddo
         enddo
      enddo
      endif ! end of iupwnd check for horizontal advection
c
c ----------- SGS z-contributions
c
      if(izs .eq. 1) then
         do j=iys,iye
         do i=1,nnx
            tflx(i,j,1,2) = -t1tbc(i,j,iscl)
            tflx(i,j,2,2) = -t2tbc(i,j,iscl)
            tflx(i,j,3,2) = -t3tbc(i,j,iscl)
         enddo
         enddo
      else
         call grads_u(scl_g,izs-1,tx,ty,iscl)
         do n=1,3
         do j=iys,iye
         do i=1,nnx
            fac   = weit(izs-1)*vis_s(i,j,izs) +
     +              weitm(izs-1)*vis_s(i,j,izs-1)
            tflx(i,j,n,2) = fac*scl_g(i,j,n)
         enddo
         enddo
         enddo
      endif
      do k=izs,ize
         kp1 = k + 1
         km1 = k - 1
         do n=1,3
         do j=iys,iye
         do i=1,nnx
            tflx(i,j,n,1) = tflx(i,j,n,2)
         enddo
         enddo
         enddo
         call grads_u(scl_g,k,tx,ty,iscl)
         do n=1,3
         do j=iys,iye
         do i=1,nnx
            fac   = weit(k)*vis_s(i,j,kp1) + weitm(k)*vis_s(i,j,k)
            tflx(i,j,n,2) = fac*scl_g(i,j,n)
         enddo
         enddo
         enddo
         do j=iys,iye
         do i=1,nnx
            gm_xu      = weit(k)*gm(i,j,kp1,1)*gm(i,j,kp1,2) +
     +                   weitm(k)*gm(i,j,k,1)*gm(i,j,k,2)
            gm_xl      = weit(km1)*gm(i,j,k,1)*gm(i,j,k,2) +
     +                   weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,2)
            gm_yu      = weit(k)*gm(i,j,kp1,1)*gm(i,j,kp1,3) +
     +                   weitm(k)*gm(i,j,k,1)*gm(i,j,k,3)
            gm_yl      = weit(km1)*gm(i,j,k,1)*gm(i,j,k,3) +
     +                   weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,3)
            tflx_up    = tflx(i,j,1,2)*gm_xu + tflx(i,j,2,2)*gm_yu +
     +                   tflx(i,j,3,2)
            tflx_lo    = tflx(i,j,1,1)*gm_xl + tflx(i,j,2,1)*gm_yl +
     +                   tflx(i,j,3,1)
            r4(i,j,iscl,k) = r4(i,j,iscl,k) +
     +                       (tflx_up - tflx_lo)*dzw_i(k)
         enddo
         enddo
c
c -------- save SGS fluxes for printout, gather sums on exit
c
         if(istep .eq. 1) then
            wtsb(k,iscl) = 0.0
            do j=iys,iye
            do i=1,nnx
               wtsb(k,iscl) = wtsb(k,iscl) - tflx(i,j,3,2)
            enddo
            enddo
            wtsb(k,iscl) = wtsb(k,iscl)*fnxy
         endif
c
c ---------- end z loop
c
      enddo
c
c ---------- resolved vertical advection
c
      if(iupwnd .ne. 1) then
c
c --------- skew symmetric advective form for
c           vertical flux = 0.5(wf*dt/dz + d/dz(wf*t))
c           no grid speed
c
      do k=izs,ize
         km1 = k - 1
         kp1 = k + 1
         do j=iys,iye
         do i=1,nnx
            tzm     = (t(i,j,iscl,k) - t(i,j,iscl,km1))*dzu_i(k)
            tz      = (t(i,j,iscl,kp1) - t(i,j,iscl,k))*dzu_i(kp1)
            theta_u = weitm(k)*t(i,j,iscl,k) +
     +                weit(k)*t(i,j,iscl,kp1)
            theta_l = weit(km1)*t(i,j,iscl,k) +
     +                weitm(km1)*t(i,j,iscl,km1)
            t_adv   = 0.5*(wf(i,j,k)*theta_u -
     +                     wf(i,j,km1)*theta_l)*dzw_i(k) +
     +                0.25*(wf(i,j,k)*tz + wf(i,j,km1)*tzm)
            r4(i,j,iscl,k) = r4(i,j,iscl,k) - t_adv
         enddo
         enddo
      enddo
c
      else
c
c ----------- z-direction special with monotone
c
      do k=izs,ize
         km2 = k - 2
         km1 = k - 1
         kp1 = k + 1
         kp2 = k + 2
         if(k .eq. 1) then
              do j=iys,iye
              do i=1,nnx
                 wflx_l      = wf(i,j,km1) - gm(i,j,km1,5)
                 wflx_u      = wf(i,j,k)   - gm(i,j,k,5)
                 flux_l(i,j) = sgn*0.5*wflx_l*
     +                         (t(i,j,iscl,km1)+t(i,j,iscl,k))
                 flux_u(i,j) =
     +           amax1(sgn*wflx_u,0.)*(t(i,j,iscl,k) +
     +           rlim(t(i,j,iscl,kp1),t(i,j,iscl,k),
     +                t(i,j,iscl,km1))) +
     +           amin1(sgn*wflx_u,0.)*(t(i,j,iscl,kp1) +
     +           rlim(t(i,j,iscl,k),t(i,j,iscl,kp1),
     +                t(i,j,iscl,kp2)))
              enddo
              enddo
         else if(k .eq. nnz) then
              do j=iys,iye
              do i=1,nnx
                 wflx_l      = wf(i,j,km1) - gm(i,j,km1,5)
                 wflx_u      = wf(i,j,k) - gm(i,j,k,5)
                 flux_u(i,j) = sgn*0.5*wflx_u*
     +                        (t(i,j,iscl,kp1)+t(i,j,iscl,k))
                 flux_l(i,j) =
     +           amax1(sgn*wflx_l,0.)*(t(i,j,iscl,km1) +
     +           rlim(t(i,j,iscl,k),t(i,j,iscl,km1),
     +                t(i,j,iscl,km2))) +
     +           amin1(sgn*wflx_l,0.)*(t(i,j,iscl,k) +
     +           rlim(t(i,j,iscl,km1),t(i,j,iscl,k),
     +                t(i,j,iscl,kp1)))
              enddo
              enddo
         else
              do j=iys,iye
              do i=1,nnx
                 wflx_l      = wf(i,j,km1) - gm(i,j,km1,5)
                 wflx_u      = wf(i,j,k) - gm(i,j,k,5)
                 flux_u(i,j) =
     +           amax1(sgn*wflx_u,0.)*(t(i,j,iscl,k) +
     +           rlim(t(i,j,iscl,kp1),t(i,j,iscl,k),
     +                t(i,j,iscl,km1))) +
     +           amin1(sgn*wflx_u,0.)*(t(i,j,iscl,kp1) +
     +           rlim(t(i,j,iscl,k),t(i,j,iscl,kp1),
     +                t(i,j,iscl,kp2)))
                 flux_l(i,j) =
     +           amax1(sgn*wflx_l,0.)*(t(i,j,iscl,km1) +
     +           rlim(t(i,j,iscl,k),t(i,j,iscl,km1),
     +                t(i,j,iscl,km2))) +
     +           amin1(sgn*wflx_l,0.)*(t(i,j,iscl,k) +
     +           rlim(t(i,j,iscl,km1),t(i,j,iscl,k),
     +                t(i,j,iscl,kp1)))
              enddo
              enddo
         endif
c
c ---------- sum vertical monotone flux
c
         do j=iys,iye
         do i=1,nnx
            r4(i,j,iscl,k) = r4(i,j,iscl,k)
     +          - sgn*(flux_u(i,j) - flux_l(i,j))*dzw_i(k)
         enddo
         enddo
c
c ---------- end z loop
c
      enddo
c
c -------- end monotone if block
c
      endif
c
      return
      end
      subroutine rhs_scl_old(istep,iscl)
c
c ------ get right hand side of scalar equation (iscl)
c        monotone scalar fluxes only in z
c        for pencil size (nnx, iys:iye, izs:ize)
c        care is taken so that if monotone is on then
c        conservative horizontal flux form is used!
c        uses flux conservation form for advection always
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      real fnt1(nnx,iys:iye), fntd(nnx,iys:iye,izs:ize)
      real tx(nnx,iys:iye,izs-1:ize+1), ty(nnx,iys:iye,izs-1:ize+1)
      real flux_u(nnx,iys:iye), flux_l(nnx,iys:iye)
      real tflx(nnx,iys:iye,3,2), scl_g(nnx,iys:iye,3)
c
c --------- set sign for ocean simulations that use monotone
c
      sgn = 1.0
      if(iocean .eq. 1) sgn = -1.0
      if(iupwnd .eq. 1) then
         upfac = 1.0
      else
         upfac = 0.5
      endif
c
c --------- get x & y derivatives for [izs-1:ize+1]
c
      do k=izs-1,ize+1
         do j=iys,iye
         do i=1,nnx
            tx(i,j,k) = t(i,j,iscl,k)
            ty(i,j,k) = t(i,j,iscl,k)
         enddo
         enddo
         call xderivp(tx(1,iys,k),trigx(1,1),xk(1),nnx,iys,iye)
      enddo
      call yd_mpi(ty(1,iys,izs-1),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
c
c -------- SGS term tflx(.,.,1) + resolved flux at cell center
c
      do k=izs,ize
         call grads_c(scl_g,k,tx,ty,iscl)
         do j=iys,iye
         do i=1,nnx
            tflx(i,j,1,1) = vis_s(i,j,k)*scl_g(i,j,1)
            tflx(i,j,2,1) = vis_s(i,j,k)*scl_g(i,j,2)
            fnt1(i,j)     = gm(i,j,k,1)*tflx(i,j,1,1) -
     +                      upfac*uf(i,j,k)*t(i,j,iscl,k)
            fntd(i,j,k)   = gm(i,j,k,1)*tflx(i,j,2,1) -
     +                      upfac*vf(i,j,k)*t(i,j,iscl,k)
         enddo
         enddo
         call xderivp(fnt1(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         do j=iys,iye
         do i=1,nnx
            r4(i,j,iscl,k) = fnt1(i,j)
         enddo
         enddo
c
c -------- save SGS fluxes for printout
c
         if(istep .eq. 1) then
            utsb(k,iscl) = 0.0
            vtsb(k,iscl) = 0.0
            do j=iys,iye
            do i=1,nnx
               utsb(k,iscl) = utsb(k,iscl) - tflx(i,j,1,1)
               vtsb(k,iscl) = vtsb(k,iscl) - tflx(i,j,2,1)
            enddo
            enddo
            utsb(k,iscl) = utsb(k,iscl)*fnxy
            vtsb(k,iscl) = vtsb(k,iscl)*fnxy
         endif
      enddo
c
c ---------- SGS fluxes tflx(.,.,2) that depend on
c            y-derivative at cell centers
c
      call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            r4(i,j,iscl,k) = r4(i,j,iscl,k) + fntd(i,j,k)
         enddo
         enddo
      enddo
c
c ------- resolved advection for skew symmetric form
c         for no monotone
c
      if(iupwnd .ne. 1) then
c
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            r4(i,j,iscl,k) = r4(i,j,iscl,k) -
     +            upfac*(uf(i,j,k)*tx(i,j,k) + vf(i,j,k)*ty(i,j,k))
         enddo
         enddo
      enddo
      endif
c
c ----------- SGS z-contributions
c
      if(izs .eq. 1) then
         do j=iys,iye
         do i=1,nnx
            tflx(i,j,1,2) = -t1tbc(i,j,iscl)
            tflx(i,j,2,2) = -t2tbc(i,j,iscl)
            tflx(i,j,3,2) = -t3tbc(i,j,iscl)
         enddo
         enddo
      else
         call grads_u(scl_g,izs-1,tx,ty,iscl)
         do n=1,3
         do j=iys,iye
         do i=1,nnx
            fac   = weit(izs-1)*vis_s(i,j,izs) +
     +              weitm(izs-1)*vis_s(i,j,izs-1)
            tflx(i,j,n,2) = fac*scl_g(i,j,n)
         enddo
         enddo
         enddo
      endif
      do k=izs,ize
         kp1 = k + 1
         km1 = k - 1
         do n=1,3
         do j=iys,iye
         do i=1,nnx
            tflx(i,j,n,1) = tflx(i,j,n,2)
         enddo
         enddo
         enddo
         call grads_u(scl_g,k,tx,ty,iscl)
         do n=1,3
         do j=iys,iye
         do i=1,nnx
            fac   = weit(k)*vis_s(i,j,kp1) + weitm(k)*vis_s(i,j,k)
            tflx(i,j,n,2) = fac*scl_g(i,j,n)
         enddo
         enddo
         enddo
         do j=iys,iye
         do i=1,nnx
            gm_xu      = weit(k)*gm(i,j,kp1,1)*gm(i,j,kp1,2) +
     +                   weitm(k)*gm(i,j,k,1)*gm(i,j,k,2)
            gm_xl      = weit(km1)*gm(i,j,k,1)*gm(i,j,k,2) +
     +                   weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,2)
            gm_yu      = weit(k)*gm(i,j,kp1,1)*gm(i,j,kp1,3) +
     +                   weitm(k)*gm(i,j,k,1)*gm(i,j,k,3)
            gm_yl      = weit(km1)*gm(i,j,k,1)*gm(i,j,k,3) +
     +                   weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,3)
            tflx_up    = tflx(i,j,1,2)*gm_xu + tflx(i,j,2,2)*gm_yu +
     +                   tflx(i,j,3,2)
            tflx_lo    = tflx(i,j,1,1)*gm_xl + tflx(i,j,2,1)*gm_yl +
     +                   tflx(i,j,3,1)
            r4(i,j,iscl,k) = r4(i,j,iscl,k) +
     +                       (tflx_up - tflx_lo)*dzw_i(k)
         enddo
         enddo
c
c -------- save SGS fluxes for printout, gather sums on exit
c
         if(istep .eq. 1) then
            wtsb(k,iscl) = 0.0
            do j=iys,iye
            do i=1,nnx
               wtsb(k,iscl) = wtsb(k,iscl) - tflx(i,j,3,2)
            enddo
            enddo
            wtsb(k,iscl) = wtsb(k,iscl)*fnxy
         endif
c
c ---------- end z loop
c
      enddo
c
c ---------- resolved vertical advection
c
      if(iupwnd .ne. 1) then
c
c --------- skew symmetric advective form for
c           vertical flux = 0.5(wf*dt/dz + d/dz(wf*t))
c           no grid speed
c
      do k=izs,ize
         km1 = k - 1
         kp1 = k + 1
         do j=iys,iye
         do i=1,nnx
            tzm     = (t(i,j,iscl,k) - t(i,j,iscl,km1))*dzu_i(k)
            tz      = (t(i,j,iscl,kp1) - t(i,j,iscl,k))*dzu_i(kp1)
            theta_u = weitm(k)*t(i,j,iscl,k) +
     +                weit(k)*t(i,j,iscl,kp1)
            theta_l = weit(km1)*t(i,j,iscl,k) +
     +                weitm(km1)*t(i,j,iscl,km1)
            t_adv   = 0.5*(wf(i,j,k)*theta_u -
     +                     wf(i,j,km1)*theta_l)*dzw_i(k) +
     +                0.25*(wf(i,j,k)*tz + wf(i,j,km1)*tzm)
            r4(i,j,iscl,k) = r4(i,j,iscl,k) - t_adv
         enddo
         enddo
      enddo
c
      else
c
c ----------- z-direction special with monotone
c
      do k=izs,ize
         km2 = k - 2
         km1 = k - 1
         kp1 = k + 1
         kp2 = k + 2
         if(k .eq. 1) then
              do j=iys,iye
              do i=1,nnx
                 wflx_l      = wf(i,j,km1) - gm(i,j,km1,5)
                 wflx_u      = wf(i,j,k)   - gm(i,j,k,5)
                 flux_l(i,j) = sgn*0.5*wflx_l*
     +                         (t(i,j,iscl,km1)+t(i,j,iscl,k))
                 flux_u(i,j) =
     +           amax1(sgn*wflx_u,0.)*(t(i,j,iscl,k) +
     +           rlim(t(i,j,iscl,kp1),t(i,j,iscl,k),
     +                t(i,j,iscl,km1))) +
     +           amin1(sgn*wflx_u,0.)*(t(i,j,iscl,kp1) +
     +           rlim(t(i,j,iscl,k),t(i,j,iscl,kp1),
     +                t(i,j,iscl,kp2)))
              enddo
              enddo
         else if(k .eq. nnz) then
              do j=iys,iye
              do i=1,nnx
                 wflx_l      = wf(i,j,km1) - gm(i,j,km1,5)
                 wflx_u      = wf(i,j,k) - gm(i,j,k,5)
                 flux_u(i,j) = sgn*0.5*wflx_u*
     +                        (t(i,j,iscl,kp1)+t(i,j,iscl,k))
                 flux_l(i,j) =
     +           amax1(sgn*wflx_l,0.)*(t(i,j,iscl,km1) +
     +           rlim(t(i,j,iscl,k),t(i,j,iscl,km1),
     +                t(i,j,iscl,km2))) +
     +           amin1(sgn*wflx_l,0.)*(t(i,j,iscl,k) +
     +           rlim(t(i,j,iscl,km1),t(i,j,iscl,k),
     +                t(i,j,iscl,kp1)))
              enddo
              enddo
         else
              do j=iys,iye
              do i=1,nnx
                 wflx_l      = wf(i,j,km1) - gm(i,j,km1,5)
                 wflx_u      = wf(i,j,k) - gm(i,j,k,5)
                 flux_u(i,j) =
     +           amax1(sgn*wflx_u,0.)*(t(i,j,iscl,k) +
     +           rlim(t(i,j,iscl,kp1),t(i,j,iscl,k),
     +                t(i,j,iscl,km1))) +
     +           amin1(sgn*wflx_u,0.)*(t(i,j,iscl,kp1) +
     +           rlim(t(i,j,iscl,k),t(i,j,iscl,kp1),
     +                t(i,j,iscl,kp2)))
                 flux_l(i,j) =
     +           amax1(sgn*wflx_l,0.)*(t(i,j,iscl,km1) +
     +           rlim(t(i,j,iscl,k),t(i,j,iscl,km1),
     +                t(i,j,iscl,km2))) +
     +           amin1(sgn*wflx_l,0.)*(t(i,j,iscl,k) +
     +           rlim(t(i,j,iscl,km1),t(i,j,iscl,k),
     +                t(i,j,iscl,kp1)))
              enddo
              enddo
         endif
c
c ---------- sum vertical monotone flux
c
         do j=iys,iye
         do i=1,nnx
            r4(i,j,iscl,k) = r4(i,j,iscl,k)
     +          - sgn*(flux_u(i,j) - flux_l(i,j))*dzw_i(k)
         enddo
         enddo
c
c ---------- end z loop
c
      enddo
c
c -------- end monotone if block
c
      endif
c
      return
      end
      subroutine grads_u(scl_u,k,tx,ty,iscl)
c
c ---------- get scalar gradient at upper face at the k-level
c
c        gradient definitions
c
c        tx   = txi  + tzet*gm_x
c        ty   = teta + tzet*gm_y
c        tz   = tzet*gm_z
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real scl_u(nnx,iys:iye,3)
      real tx(nnx,iys:iye,izs-1:ize+1), ty(nnx,iys:iye,izs-1:ize+1)
c
      km1  = k - 1
      kp1  = k + 1
      w1   = weit(k)
      w1m  = weitm(k)
      w2   = weit(km1)
      w2m  = weitm(km1)
c
c -------- at upper face
c
      do j=iys,iye
      do i=1,nnx
         tzet  = (t(i,j,iscl,kp1) - t(i,j,iscl,k))*dzu_i(kp1)
         txi   = w1*tx(i,j,kp1) + w1m*tx(i,j,k)
         teta  = w1*ty(i,j,kp1) + w1m*ty(i,j,k)
         gm_x  = w1*gm(i,j,kp1,2) + w1m*gm(i,j,k,2)
         gm_y  = w1*gm(i,j,kp1,3) + w1m*gm(i,j,k,3)
         gm_z  = w1*gm(i,j,kp1,4) +  w1m*gm(i,j,k,4)
         scl_u(i,j,1) = txi  + tzet*gm_x
         scl_u(i,j,2) = teta + tzet*gm_y
         scl_u(i,j,3) = tzet*gm_z
      enddo
      enddo
c
      return
      end
      subroutine grads_c_dns(scl_c,k,tx,ty,iscl)
c
c ---------- get scalar gradient at the cell center
c
c        gradient definitions
c
c        tx   = txi  + tzet*gm_x
c        ty   = teta + tzet*gm_y
c        tz   = tzet*gm_z
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real scl_c(nnx,iys:iye,3)
      real tx(nnx,iys:iye,izs-1:ize+1), ty(nnx,iys:iye,izs-1:ize+1)
c
      km1  = k - 1
      kp1  = k + 1
      w1   = weit(k)
      w1m  = weitm(k)
      w2   = weit(km1)
      w2m  = weitm(km1)
c
      if(k == 1) then
c
      do j=iys,iye
      do i=1,nnx
         tzet  = (t(i,j,iscl,kp1) - t(i,j,iscl,k))*dzu_i(kp1)
         scl_c(i,j,1) = tx(i,j,k) + tzet*gm(i,j,k,2)
         scl_c(i,j,2) = ty(i,j,k) + tzet*gm(i,j,k,3)
         scl_c(i,j,3) = tzet*gm(i,j,k,4)
      enddo
      enddo
      else
c
c -------- at cell center
c
      do j=iys,iye
      do i=1,nnx
         tup   = w1*t(i,j,iscl,kp1) + w1m*t(i,j,iscl,k)
         tlo   = w2*t(i,j,iscl,k) + w2m*t(i,j,iscl,km1)
         tzet  = (tup - tlo)*dzw_i(k)
         scl_c(i,j,1) = tx(i,j,k) + tzet*gm(i,j,k,2)
         scl_c(i,j,2) = ty(i,j,k) + tzet*gm(i,j,k,3)
         scl_c(i,j,3) = tzet*gm(i,j,k,4)
      enddo
      enddo
c
      endif
c
      return
      end
      subroutine grads_c(scl_c,k,tx,ty,iscl)
c
c ---------- get scalar gradient at the cell center
c
c        gradient definitions
c
c        tx   = txi  + tzet*gm_x
c        ty   = teta + tzet*gm_y
c        tz   = tzet*gm_z
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real scl_c(nnx,iys:iye,3)
      real tx(nnx,iys:iye,izs-1:ize+1), ty(nnx,iys:iye,izs-1:ize+1)
c
      km1  = k - 1
      kp1  = k + 1
      w1   = weit(k)
      w1m  = weitm(k)
      w2   = weit(km1)
      w2m  = weitm(km1)
c
      if(k == 1) then
c
      do j=iys,iye
      do i=1,nnx
         tzet  = (t(i,j,iscl,kp1) - t(i,j,iscl,k))*dzu_i(kp1)
         scl_c(i,j,1) = tx(i,j,k) + tzet*gm(i,j,k,2)
         scl_c(i,j,2) = ty(i,j,k) + tzet*gm(i,j,k,3)
         scl_c(i,j,3) = tzet*gm(i,j,k,4)
      enddo
      enddo
      else
c
c -------- at cell center
c
      do j=iys,iye
      do i=1,nnx
         tup   = w1*t(i,j,iscl,kp1) + w1m*t(i,j,iscl,k)
         tlo   = w2*t(i,j,iscl,k) + w2m*t(i,j,iscl,km1)
         tzet  = (tup - tlo)*dzw_i(k)
         scl_c(i,j,1) = tx(i,j,k) + tzet*gm(i,j,k,2)
         scl_c(i,j,2) = ty(i,j,k) + tzet*gm(i,j,k,3)
         scl_c(i,j,3) = tzet*gm(i,j,k,4)
      enddo
      enddo
c
      endif
c
      return
      end
      subroutine dns_vis
      use particles
      use pars
      use fields
      implicit none

!     In DNS mode, just set the molecular viscosity (and scalar diffusivities)
!     Also, to make the rest of code work, set the rhs of e equation to 0

      !Both for air at the moment:
      vis_m = nuf
      !Use Prantdl number for thermal diffusivity:
      !NOTE: ASSUMING Sc = Pr since vis_s only has 1 field!
      vis_s(1:nnx,iys:iye,izs-1:ize+1) = nuf/Pra

      r5 = 0.0
      e = 0.0

      end
      subroutine tke_vis(istep)
c
c ------ get viscosity using deardorff tke model with
c        stability correction. fixes for surface layer.
c        get rhs of e-equation/J. field variables are
c        (u,v,w,e,theta) all at cell center. get viscosity
c        at cell center also.
c        "length scale" based on the volume in physical space
c
c        volume = gm(.,.,.,1)*d_zeta*d_eta*d_xi
c        (fletcher, vol2 1988, p. 51)
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      real ex(nnx,iys:iye,izs-1:ize+1), ey(nnx,iys:iye,izs-1:ize+1)
      real fnt1(nnx,iys:iye), fntd(nnx,iys:iye,izs:ize)
      real dissp(nnx,iys:iye), alk(nnx,iys:iye), sij(nnx,iys:iye,6)
      real eflx(nnx,iys:iye,3,2), e_g(nnx,iys:iye,3)
      real rj3(nnx,iys:iye)
c
c -------- see par statments for setting k_stab
c
      do k=max(izs-1,1),min(nnz,ize+1)
c
      kp1  = k + 1
      km1  = k - 1
      w1   = weit(k)
      w1m  = weitm(k)
      w2   = weit(km1)
      w2m  = weitm(km1)
      dslk = dsl_z(k)
c
      do j=iys,iye
      do i=1,nnx
         rj3(i,j) = gm(i,j,k,1)**pow3
      enddo
      enddo
c
      if(k .gt. 0) dslk  = amin1(dsl_z(k),vk*abs(z(k))/csmag)
      almin = almin_c*dsl_z(k)
      if(k .eq. 0 .or. k .eq. nnz+1) then
         dfack = 1.0
      else
         dfack = dfac(k)
      endif
c     if(ivis .eq. 1 .and. k .le. nmatch) then
c
      if(k .le. k_stab) then
c
c --------------- no stability corrected length scales
c
         do j=iys,iye
         do i=1,nnx
            alk(i,j) = dslk*rj3(i,j)
         end do
         end do
      else
         do j=iys,iye
         do i=1,nnx
            alk(i,j) = dslk*rj3(i,j)
            tup  = w1*t(i,j,1,kp1) + w1m*t(i,j,1,k)
            tlo  = w2*t(i,j,1,k) + w2m*t(i,j,1,km1)
            stab = gm(i,j,k,4)*batag*(tup - tlo)*dzw_i(k)
            if(stab .gt. stabmin) then
              als = stab_c*sqrt(e(i,j,k)/stab)
              alk(i,j) = amin1(dslk,als)*rj3(i,j)
            endif
            alk(i,j)  = amax1(almin*rj3(i,j),alk(i,j))
         enddo
         enddo
      endif
      do j=iys,iye
      do i=1,nnx
         vis_m(i,j,k) = ck*alk(i,j)*sqrt(e(i,j,k))*dfack
         vis_s(i,j,k) = (1.0 + 2.0*alk(i,j)/(dslk*rj3(i,j)))*
     +                  vis_m(i,j,k)
      enddo
      enddo
c
c -------------- special case at bottom and top
c
      if(k == 1) then
         do j=iys,iye
         do i=1,nnx
            vis_m(i,j,km1) = vis_m(i,j,k)
            vis_s(i,j,km1) = vis_s(i,j,k)
         enddo
         enddo
      elseif(k == nnz) then
         do j=iys,iye
         do i=1,nnx
            vis_m(i,j,kp1) = vis_m(i,j,k)
            vis_s(i,j,kp1) = vis_s(i,j,k)
         enddo
         enddo
      endif
c
      if(k .ge. izs .and. k .le. ize) then
c
c -------- to avoid storing length scale get dissipation
c          and shear production with cell centered strains
c          special at lower boundary
c
      call strain_c(sij,k)
      if(k .eq. 1) then
         do j=iys,iye
         do i=1,nnx
            fac       = 2.0*(weit(k)*vis_m(i,j,k+1) +
     +                       weitm(k)*vis_m(i,j,k))
            tau11_s11 = 0.5*(fac*sij(i,j,1) - t11bc(i,j))*sij(i,j,1)
            tau12_s12 = 0.5*(fac*sij(i,j,2) - t12bc(i,j))*sij(i,j,2)
            tau13_s13 = 0.5*(fac*sij(i,j,3) - t13bc(i,j))*sij(i,j,3)
            tau22_s22 = 0.5*(fac*sij(i,j,4) - t22bc(i,j))*sij(i,j,4)
            tau23_s23 = 0.5*(fac*sij(i,j,5) - t23bc(i,j))*sij(i,j,5)
            tau33_s33 = 0.5*(fac*sij(i,j,6) - t33bc(i,j))*sij(i,j,6)
            fnt1(i,j) = tau11_s11 + tau22_s22 + tau33_s33 +
     +                  2.0*(tau12_s12 + tau13_s13 + tau23_s23)
         enddo
         enddo
      else
         do j=iys,iye
         do i=1,nnx
            fnt1(i,j)  = 2.0*vis_m(i,j,k)*(
     +            sij(i,j,1)**2 + sij(i,j,4)**2 + sij(i,j,6)**2 +
     +            2.0*(sij(i,j,2)**2 + sij(i,j,3)**2 + sij(i,j,5)**2))
         enddo
         enddo
      endif
c
      do j=iys,iye
      do i=1,nnx
         dissp(i,j) = (0.19 + 0.74*alk(i,j)/(dslk*rj3(i,j)))*
     +                 e(i,j,k)*sqrt(e(i,j,k))/alk(i,j)
         r5(i,j,k) = (fnt1(i,j) - dissp(i,j))*gm(i,j,k,1)
      enddo
      enddo
c
      if(istep .eq. 1) then
         t_diss(k) = 0.0
         shrz(k)   = 0.0
         triz(k)   = 0.0
         do j=iys,iye
         do i=1,nnx
            t_diss(k) = t_diss(k) + dissp(i,j)
            shrz(k)   = shrz(k) + fnt1(i,j)
         enddo
         enddo
         t_diss(k) = t_diss(k)*fnxy
         shrz(k)   = shrz(k)*fnxy
      endif
c
      endif   ! endif for k = [izs,ize]
c
c -------------- end z loop
c
      enddo
c
c -------------- if special 2 part surface layer model is on
c                get "mean" viscosity
CHECK 2-part for curvy?
c
c     do k=izs-1,ize
c        km1         = k - 1
c        kp1         = k + 1
c        vis_mean(k) = 0.0
c        if(ivis .eq. 1 .and. k .le. nmatch) then
c           if(k .le. 1) then
c             vis_mean(k) = xksurf
c           else
c             stravg = sqrt((u_mn(kp1)-u_mn(k))**2 +
c    +              (v_mn(kp1)-v_mn(k))**2)*abs(dzu_i(kp1))
c             vis_mean(k) = xksurf*viscon*stravg
c           endif
c        endif
c     enddo
c
c --------- get x & y derivatives for [izs-1:ize+1]
c
      do k=izs-1,ize+1
         do j=iys,iye
         do i=1,nnx
            ex(i,j,k) = e(i,j,k)
            ey(i,j,k) = e(i,j,k)
         enddo
         enddo
         call xderivp(ex(1,iys,k),trigx(1,1),xk(1),nnx,iys,iye)
      enddo
      call yd_mpi(ey(1,iys,izs-1),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
c
c -------- diffusion term eflx(.,.,1) + resolved flux at cell center
c
      do k=izs,ize
         call grd_e_c(e_g,k,ex,ey)
         do j=iys,iye
         do i=1,nnx
            eflx(i,j,1,1) = 2.0*vis_m(i,j,k)*e_g(i,j,1)
            eflx(i,j,2,1) = 2.0*vis_m(i,j,k)*e_g(i,j,2)
            fnt1(i,j)     = gm(i,j,k,1)*eflx(i,j,1,1) -
     +                      uf(i,j,k)*e(i,j,k)
            fntd(i,j,k)   = gm(i,j,k,1)*eflx(i,j,2,1) -
     +                      vf(i,j,k)*e(i,j,k)
         enddo
         enddo
         call xderivp(fnt1(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         do j=iys,iye
         do i=1,nnx
            r5(i,j,k) = r5(i,j,k) + fnt1(i,j)
         enddo
         enddo
      enddo
c
c -------- diffusion term eflx(.,.,2) + resolved flux at cell center
c
      call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            r5(i,j,k) = r5(i,j,k) + fntd(i,j,k)
         enddo
         enddo
      enddo
c
c --------- z-diffusion term + resolved vertical flux
c           no divergence of flux at first gridpoint
c
      k   = izs
      km1 = izs - 1
      kp1 = izs + 1
      if(izs .eq. 1) then
         do j=iys,iye
         do i=1,nnx
            e(i,j,km1)  = e(i,j,k)
            ex(i,j,km1) = ex(i,j,k)
            ey(i,j,km1) = ey(i,j,k)
         enddo
         enddo
      endif
      call grd_e_u(e_g,k,ex,ey)
      do n=1,3
      do j=iys,iye
      do i=1,nnx
         fac   = weit(k)*vis_m(i,j,kp1) +
     +              weitm(k)*vis_m(i,j,k)
         eflx(i,j,n,2) = 2.0*fac*e_g(i,j,n)
      enddo
      enddo
      enddo
c
      do k=izs,ize
         kp1 = k + 1
         km1 = k - 1
         do n=1,3
         do j=iys,iye
         do i=1,nnx
            eflx(i,j,n,1) = eflx(i,j,n,2)
         enddo
         enddo
         enddo
         call grd_e_u(e_g,k,ex,ey)
         do n=1,3
         do j=iys,iye
         do i=1,nnx
            fac   = weit(k)*vis_m(i,j,kp1) + weitm(k)*vis_m(i,j,k)
            eflx(i,j,n,2) = 2.0*fac*e_g(i,j,n)
         enddo
         enddo
         enddo
         do j=iys,iye
         do i=1,nnx
            gm_xu     = weit(k)*gm(i,j,kp1,1)*gm(i,j,kp1,2) +
     +                  weitm(k)*gm(i,j,k,1)*gm(i,j,k,2)
            gm_xl     = weit(km1)*gm(i,j,k,1)*gm(i,j,k,2) +
     +                  weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,2)
            gm_yu     = weit(k)*gm(i,j,kp1,1)*gm(i,j,kp1,3) +
     +                  weitm(k)*gm(i,j,k,1)*gm(i,j,k,3)
            gm_yl     = weit(km1)*gm(i,j,k,1)*gm(i,j,k,3) +
     +                  weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,3)
            eflx_up   = eflx(i,j,1,2)*gm_xu + eflx(i,j,2,2)*gm_yu +
     +                  eflx(i,j,3,2)
            eflx_lo   = eflx(i,j,1,1)*gm_xl + eflx(i,j,2,1)*gm_yl +
     +                  eflx(i,j,3,1)
            r5(i,j,k) = r5(i,j,k) + (eflx_up - eflx_lo)*dzw_i(k)
         enddo
         enddo
c
c ---------- resolved vertical advection, flux form
c            plus grid speed
c
         do j=iys,iye
         do i=1,nnx
            e_u       = weitm(k)*e(i,j,k) + weit(k)*e(i,j,kp1)
            e_l       = weit(km1)*e(i,j,k) + weitm(km1)*e(i,j,km1)
            wflx_u    = wf(i,j,k) - gm(i,j,k,5)
            wflx_l    = wf(i,j,km1) - gm(i,j,km1,5)
            e_adv     = (wflx_u*e_u - wflx_l*e_l)*dzw_i(k)
            r5(i,j,k) = r5(i,j,k) - e_adv
         enddo
         enddo
c
c --------- get buoyancy term/J in sgs e equation
c
      if(k .eq. 1) then
         do j=iys,iye
         do i=1,nnx
            fac   = weit(k)*vis_s(i,j,kp1) + weitm(k)*vis_s(i,j,k)
            tzet  = (t(i,j,1,kp1) - t(i,j,1,k))*dzu_i(kp1)
            r5(i,j,k) = r5(i,j,k) - fac*tzet*batag*gm(i,j,k,1)
         enddo
         enddo
      else
         do j=iys,iye
         do i=1,nnx
            tup   = weit(k)*t(i,j,1,kp1) + weitm(k)*t(i,j,1,k)
            tlo   = weit(km1)*t(i,j,1,k) + weitm(km1)*t(i,j,1,km1)
            tzet  = (tup - tlo)*dzw_i(k)
            r5(i,j,k) = r5(i,j,k) - vis_s(i,j,k)*tzet*batag*gm(i,j,k,1)
         enddo
         enddo
      endif
c
c -------- end z loop
c
      enddo
c
c ------------- 9/1989 add ihflt=1 option--mean shear does not generate sgs tke
c
c     uxymm=0.
c     uxymp=0.
c     vxymm=0.
c     vxymp=0.
c     if(ivis .eq. 1 .and. k .le. nmatch) then
c        uxymm = u_mn(k)
c        uxymp = u_mn(kp1)
c        vxymm = v_mn(k)
c        vxymp = v_mn(kp1)
c     endif
c
c -------------- end z loop
c
c     enddo
c
c -------------- leave with r5/J
c
      return
      end
      subroutine tke_vis_old(istep)
c
c ------ get viscosity using deardorff tke model with
c        stability correction. fixes for surface layer.
c        get rhs of e-equation/J. field variables are
c        (u,v,w,e,theta) all at cell center. get viscosity
c        at cell center also.
c        "length scale" based on the volume in physical space
c
c        volume = gm(.,.,.,1)*d_zeta*d_eta*d_xi
c        (fletcher, vol2 1988, p. 51)
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      real ex(nnx,iys:iye,izs-1:ize+1), ey(nnx,iys:iye,izs-1:ize+1)
      real fnt1(nnx,iys:iye), fntd(nnx,iys:iye,izs:ize)
      real dissp(nnx,iys:iye), alk(nnx,iys:iye), sij(nnx,iys:iye,6)
      real eflx(nnx,iys:iye,3,2), e_g(nnx,iys:iye,3)
c
      do k=max(izs-1,1),min(nnz,ize+1)
c
      kp1  = k + 1
      km1  = k - 1
      w1   = weit(k)
      w1m  = weitm(k)
      w2   = weit(km1)
      w2m  = weitm(km1)
      dslk = dsl_z(k)
      if(k .gt. 0) dslk  = amin1(dsl_z(k),vk*abs(z(k))/csmag)
      almin = almin_c*dsl_z(k)
      if(k .eq. 0 .or. k .eq. nnz+1) then
         dfack = 1.0
      else
         dfack = dfac(k)
      endif
c     if(ivis .eq. 1 .and. k .le. nmatch) then
c
CHECK FIX
c
      if(k .le. ((3*nnz)/6)) then
c
c --------------- no stability corrected length scales
c
         do j=iys,iye
         do i=1,nnx
            alk(i,j) = dslk*gm(i,j,k,1)
         end do
         end do
      else
         do j=iys,iye
         do i=1,nnx
            alk(i,j) = dslk*gm(i,j,k,1)
            tup  = w1*t(i,j,1,kp1) + w1m*t(i,j,1,k)
            tlo  = w2*t(i,j,1,k) + w2m*t(i,j,1,km1)
            stab = gm(i,j,k,4)*batag*(tup - tlo)*dzw_i(k)
            if(stab .gt. stabmin) then
              als = stab_c*sqrt(e(i,j,k)/stab)
              alk(i,j) = amin1(dslk,als)*gm(i,j,k,1)
            endif
            alk(i,j)  = amax1(almin,alk(i,j))
         enddo
         enddo
      endif
      do j=iys,iye
      do i=1,nnx
         vis_m(i,j,k) = ck*alk(i,j)*sqrt(e(i,j,k))*dfack
         vis_s(i,j,k) = (1.0 + 2.0*alk(i,j)/(dslk*gm(i,j,k,1)))*
     +                  vis_m(i,j,k)
      enddo
      enddo
c
c -------------- special case at bottom and top
c
      if(k == 1) then
         do j=iys,iye
         do i=1,nnx
            vis_m(i,j,km1) = vis_m(i,j,k)
            vis_s(i,j,km1) = vis_s(i,j,k)
         enddo
         enddo
      elseif(k == nnz) then
         do j=iys,iye
         do i=1,nnx
            vis_m(i,j,kp1) = vis_m(i,j,k)
            vis_s(i,j,kp1) = vis_s(i,j,k)
         enddo
         enddo
      endif
c
      if(k .ge. izs .and. k .le. ize) then
c
c -------- to avoid storing length scale get dissipation
c          and shear production with cell centered strains
c
      call strain_c(sij,k)
      do j=iys,iye
      do i=1,nnx
         dissp(i,j) = (0.19 + 0.74*alk(i,j)/(dslk*gm(i,j,k,1)))*
     +                 e(i,j,k)*sqrt(e(i,j,k))/alk(i,j)
         fnt1(i,j)  = 2.0*vis_m(i,j,k)*(
     +            sij(i,j,1)**2 + sij(i,j,4)**2 + sij(i,j,6)**2 +
     +            2.0*(sij(i,j,2)**2 + sij(i,j,3)**2 + sij(i,j,5)**2))
         r5(i,j,k) = (fnt1(i,j) - dissp(i,j))*gm(i,j,k,1)
      enddo
      enddo
c
      if(istep .eq. 1) then
         t_diss(k) = 0.0
         shrz(k)   = 0.0
         triz(k)   = 0.0
         do j=iys,iye
         do i=1,nnx
            t_diss(k) = t_diss(k) + dissp(i,j)
            shrz(k)   = shrz(k) + fnt1(i,j)
         enddo
         enddo
         t_diss(k) = t_diss(k)*fnxy
         shrz(k)   = shrz(k)*fnxy
      endif
c
      endif   ! endif for k = [izs,ize]
c
c -------------- end z loop
c
      enddo
c
c -------------- if special 2 part surface layer model is on
c                get "mean" viscosity
CHECK 2-part for curvy?
c
c     do k=izs-1,ize
c        km1         = k - 1
c        kp1         = k + 1
c        vis_mean(k) = 0.0
c        if(ivis .eq. 1 .and. k .le. nmatch) then
c           if(k .le. 1) then
c             vis_mean(k) = xksurf
c           else
c             stravg = sqrt((u_mn(kp1)-u_mn(k))**2 +
c    +              (v_mn(kp1)-v_mn(k))**2)*abs(dzu_i(kp1))
c             vis_mean(k) = xksurf*viscon*stravg
c           endif
c        endif
c     enddo
c
c --------- get x & y derivatives for [izs-1:ize+1]
c
      do k=izs-1,ize+1
         do j=iys,iye
         do i=1,nnx
            ex(i,j,k) = e(i,j,k)
            ey(i,j,k) = e(i,j,k)
         enddo
         enddo
         call xderivp(ex(1,iys,k),trigx(1,1),xk(1),nnx,iys,iye)
      enddo
      call yd_mpi(ey(1,iys,izs-1),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
c
c -------- diffusion term eflx(.,.,1) + resolved flux at cell center
c
      do k=izs,ize
         call grd_e_c(e_g,k,ex,ey)
         do j=iys,iye
         do i=1,nnx
            eflx(i,j,1,1) = 2.0*vis_m(i,j,k)*e_g(i,j,1)
            eflx(i,j,2,1) = 2.0*vis_m(i,j,k)*e_g(i,j,2)
            fnt1(i,j)     = gm(i,j,k,1)*eflx(i,j,1,1) -
     +                      uf(i,j,k)*e(i,j,k)
            fntd(i,j,k)   = gm(i,j,k,1)*eflx(i,j,2,1) -
     +                      vf(i,j,k)*e(i,j,k)
         enddo
         enddo
         call xderivp(fnt1(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         do j=iys,iye
         do i=1,nnx
            r5(i,j,k) = r5(i,j,k) + fnt1(i,j)
         enddo
         enddo
      enddo
c
c -------- diffusion term eflx(.,.,2) + resolved flux at cell center
c
      call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            r5(i,j,k) = r5(i,j,k) + fntd(i,j,k)
         enddo
         enddo
      enddo
c
c --------- z-diffusion term + resolved vertical flux
c           no divergence of flux at first gridpoint
c
      k   = izs
      km1 = izs - 1
      kp1 = izs + 1
      if(izs .eq. 1) then
         do j=iys,iye
         do i=1,nnx
            e(i,j,km1)  = e(i,j,k)
            ex(i,j,km1) = ex(i,j,k)
            ey(i,j,km1) = ey(i,j,k)
         enddo
         enddo
      endif
      call grd_e_u(e_g,k,ex,ey)
      do n=1,3
      do j=iys,iye
      do i=1,nnx
         fac   = weit(k)*vis_m(i,j,kp1) +
     +              weitm(k)*vis_m(i,j,k)
         eflx(i,j,n,2) = 2.0*fac*e_g(i,j,n)
      enddo
      enddo
      enddo
c
      do k=izs,ize
         kp1 = k + 1
         km1 = k - 1
         do n=1,3
         do j=iys,iye
         do i=1,nnx
            eflx(i,j,n,1) = eflx(i,j,n,2)
         enddo
         enddo
         enddo
         call grd_e_u(e_g,k,ex,ey)
         do n=1,3
         do j=iys,iye
         do i=1,nnx
            fac   = weit(k)*vis_m(i,j,kp1) + weitm(k)*vis_m(i,j,k)
            eflx(i,j,n,2) = 2.0*fac*e_g(i,j,n)
         enddo
         enddo
         enddo
         do j=iys,iye
         do i=1,nnx
            gm_xu     = weit(k)*gm(i,j,kp1,1)*gm(i,j,kp1,2) +
     +                  weitm(k)*gm(i,j,k,1)*gm(i,j,k,2)
            gm_xl     = weit(km1)*gm(i,j,k,1)*gm(i,j,k,2) +
     +                  weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,2)
            gm_yu     = weit(k)*gm(i,j,kp1,1)*gm(i,j,kp1,3) +
     +                  weitm(k)*gm(i,j,k,1)*gm(i,j,k,3)
            gm_yl     = weit(km1)*gm(i,j,k,1)*gm(i,j,k,3) +
     +                  weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,3)
            eflx_up   = eflx(i,j,1,2)*gm_xu + eflx(i,j,2,2)*gm_yu +
     +                  eflx(i,j,3,2)
            eflx_lo   = eflx(i,j,1,1)*gm_xl + eflx(i,j,2,1)*gm_yl +
     +                  eflx(i,j,3,1)
            r5(i,j,k) = r5(i,j,k) + (eflx_up - eflx_lo)*dzw_i(k)
         enddo
         enddo
c
c ---------- resolved vertical advection, flux form
c            plus grid speed
c
         do j=iys,iye
         do i=1,nnx
            e_u       = weitm(k)*e(i,j,k) + weit(k)*e(i,j,kp1)
            e_l       = weit(km1)*e(i,j,k) + weitm(km1)*e(i,j,km1)
            wflx_u    = wf(i,j,k) - gm(i,j,k,5)
            wflx_l    = wf(i,j,km1) - gm(i,j,km1,5)
            e_adv     = (wflx_u*e_u - wflx_l*e_l)*dzw_i(k)
            r5(i,j,k) = r5(i,j,k) - e_adv
         enddo
         enddo
c
c --------- get buoyancy term/J in sgs e equation
c
      if(k .eq. 1) then
         do j=iys,iye
         do i=1,nnx
            fac   = weit(k)*vis_s(i,j,kp1) + weitm(k)*vis_s(i,j,k)
            tzet  = (t(i,j,1,kp1) - t(i,j,1,k))*dzu_i(kp1)
            r5(i,j,k) = r5(i,j,k) - fac*tzet*batag*gm(i,j,k,1)
         enddo
         enddo
      else
         do j=iys,iye
         do i=1,nnx
            tup   = weit(k)*t(i,j,1,kp1) + weitm(k)*t(i,j,1,k)
            tlo   = weit(km1)*t(i,j,1,k) + weitm(km1)*t(i,j,1,km1)
            tzet  = (tup - tlo)*dzw_i(k)
            r5(i,j,k) = r5(i,j,k) - vis_s(i,j,k)*tzet*batag*gm(i,j,k,1)
         enddo
         enddo
      endif
c
c -------- end z loop
c
      enddo
c
c ------------- 9/1989 add ihflt=1 option--mean shear does not generate sgs tke
c
c     uxymm=0.
c     uxymp=0.
c     vxymm=0.
c     vxymp=0.
c     if(ivis .eq. 1 .and. k .le. nmatch) then
c        uxymm = u_mn(k)
c        uxymp = u_mn(kp1)
c        vxymm = v_mn(k)
c        vxymp = v_mn(kp1)
c     endif
c
c -------------- end z loop
c
c     enddo
c
c -------------- leave with r5/J
c
      return
      end
      subroutine grd_e_u(scl_u,k,ex,ey)
c
c ---------- get scalar gradient at upper face at the k-level
c
c        gradient definitions
c
c        ex   = exi  + ezet*gm_x
c        ey   = eeta + ezet*gm_y
c        ez   = tzet*gm_z
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real scl_u(nnx,iys:iye,3)
      real ex(nnx,iys:iye,izs-1:ize+1), ey(nnx,iys:iye,izs-1:ize+1)
c
      km1  = k - 1
      kp1  = k + 1
      w1   = weit(k)
      w1m  = weitm(k)
      w2   = weit(km1)
      w2m  = weitm(km1)
c
c -------- at upper face
c
      do j=iys,iye
      do i=1,nnx
         ezet  = (e(i,j,kp1) - e(i,j,k))*dzu_i(kp1)
         exi   = w1*ex(i,j,kp1) + w1m*ex(i,j,k)
         eeta  = w1*ey(i,j,kp1) + w1m*ey(i,j,k)
         gm_x  = w1*gm(i,j,kp1,2) + w1m*gm(i,j,k,2)
         gm_y  = w1*gm(i,j,kp1,3) + w1m*gm(i,j,k,3)
         gm_z  = w1*gm(i,j,kp1,4) +  w1m*gm(i,j,k,4)
         scl_u(i,j,1) = exi  + ezet*gm_x
         scl_u(i,j,2) = eeta + ezet*gm_y
         scl_u(i,j,3) = ezet*gm_z
      enddo
      enddo
c
      return
      end
      subroutine grd_e_c(scl_c,k,ex,ey)
c
c ---------- get scalar gradient at the cell center
c
c        gradient definitions
c
c        ex   = exi  + ezet*gm_x
c        ey   = eeta + ezet*gm_y
c        ez   = ezet*gm_z
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real scl_c(nnx,iys:iye,3)
      real ex(nnx,iys:iye,izs-1:ize+1), ey(nnx,iys:iye,izs-1:ize+1)
c
      km1  = k - 1
      kp1  = k + 1
      w1   = weit(k)
      w1m  = weitm(k)
      w2   = weit(km1)
      w2m  = weitm(km1)
c
      if(k == 1) then
c
      do j=iys,iye
      do i=1,nnx
         ezet  = (e(i,j,kp1) - e(i,j,k))*dzu_i(kp1)
         scl_c(i,j,1) = ex(i,j,k) + ezet*gm(i,j,k,2)
         scl_c(i,j,2) = ey(i,j,k) + ezet*gm(i,j,k,3)
         scl_c(i,j,3) = ezet*gm(i,j,k,4)
      enddo
      enddo
      else
c
c -------- at cell center
c
      do j=iys,iye
      do i=1,nnx
         eup   = w1*e(i,j,kp1) + w1m*e(i,j,k)
         elo   = w2*e(i,j,k) + w2m*e(i,j,km1)
         ezet  = (eup - elo)*dzw_i(k)
         scl_c(i,j,1) = ex(i,j,k) + ezet*gm(i,j,k,2)
         scl_c(i,j,2) = ey(i,j,k) + ezet*gm(i,j,k,3)
         scl_c(i,j,3) = ezet*gm(i,j,k,4)
      enddo
      enddo
c
      endif
c
      return
      end
      subroutine lower(it)
c
c ------ setup lower boundary condition for entire plane at (k = 1)
c        using either businger or large formulas with wind.
c        index f(.,.,2)  indicates lower.
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real sfc_flx(2+nscl)
c
      k    = 1
      km1  = k - 1
      dz_i = dzu_i(k)
c
      do j=iys,iye
      do i=1,nnx
         ebc(i,j,2)  = amax1(e(i,j,k),sml_eg)
         wbc(i,j,2)  = 0.0
         pbc(i,j,2)  = 0.0
         pbc2(i,j,2) = 0.0
         wfbc(i,j,2) = 0.0
         ufbc(i,j,2) = 0.0
         vfbc(i,j,2) = 0.0
      enddo
      enddo
c
      if(iocean .eq. 1) then ! water surface layer
         call sufto(it)
         do j=iys,iye
         do i=1,nnx
            tau13m(i,j) = -au13m
            tau23m(i,j) = -au23m
         enddo
         enddo
         do iscl=1,nscl
           do j=iys,iye
           do i=1,nnx
              taut3m(i,j,iscl) = aut3m(iscl)
           enddo
           enddo
         enddo
c
      else  !   land surface layer
c
         call suft(it)
         fac = -utau**2/(windm*sqrt(u1xy**2 + v1xy**2))
         do j=iys,iye
         do i=1,nnx
            tau13m(i,j) = fac*(windm*(u(i,j,k)+ugal-u1xy)+
     +                     wind(i,j)*u1xy)
            tau23m(i,j) = fac*(windm*(v(i,j,k)-v1xy)+
     +                     wind(i,j)*v1xy)
         enddo
         enddo
         do iscl=1,nscl
            dnom3=t10xy(iscl)*windm
            if(dnom3 .ne. 0.) then
               dnom_i = 1.0/dnom3
               do j=iys,iye
               do i=1,nnx
                  taut3m(i,j,iscl) = aut3m(iscl)*
     +                 (windm*(t(i,j,iscl,k) - t1xy(iscl))+
     +                  wind(i,j)*(t1xy(iscl) - tsfcc(iscl)))*dnom_i
               enddo
               enddo
            else
               do j=iys,iye
               do i=1,nnx
                  taut3m(i,j,iscl) = aut3m(iscl)
               enddo
               enddo
            endif
         enddo
c
      endif ! end of surface layer type
c
      do j=iys,iye
      do i=1,nnx
         t11bc(i,j) = 0.0
         t12bc(i,j) = 0.0
         t13bc(i,j) = tau13m(i,j)
         t22bc(i,j) = 0.0
         t23bc(i,j) = tau23m(i,j)
         t33bc(i,j) = 0.0
      enddo
      enddo
c
c -------- partial sums of surface fluxes and mean scalar
c
      sfc_flx(1) = 0.0
      sfc_flx(2) = 0.0
      do j=iys,iye
      do i=1,nnx
         sfc_flx(1) = sfc_flx(1) + tau13m(i,j)
         sfc_flx(2) = sfc_flx(2) + tau23m(i,j)
      enddo
      enddo
      do iscl=1,nscl
         sfc_flx(2+iscl) = 0.0
         do j=iys,iye
         do i=1,nnx
            sfc_flx(2+iscl) = sfc_flx(2+iscl) + taut3m(i,j,iscl)
         enddo
         enddo
      enddo
c
      call mpi_sum_xy(sfc_flx,myid,iss,ise,(2+nscl))
      uwsfc = sfc_flx(1)*fnxy
      vwsfc = sfc_flx(2)*fnxy
      uusfc = 0.0
      uvsfc = 0.0
      vvsfc = 0.0
      wwsfc = 0.0
      do iscl=1,nscl
         wtsfc(iscl) = sfc_flx(2+iscl)*fnxy
      enddo
c     write(nprt,2345) uwsfc, vwsfc, wtsfc(nscl), tsfcc(nscl)
 2345 format(' in lower 2345 uwsfc = ',e15.6,' vwsfc = ',e15.6,
     +       ' wtsfc = ',e15.6,' tsfcc = ',e15.6)
c
      do j=iys,iye
      do i=1,nnx
         dudz       = 2.*(u(i,j,k) + ugal)*dz_i
         dvdz       = 2.*v(i,j,k)*dz_i
         ubc(i,j,2) = u(i,j,k) - dudz*dzu(k)
         vbc(i,j,2) = v(i,j,k) - dvdz*dzu(k)
      enddo
      enddo
      do iscl=1,nscl
         do j=iys,iye
         do i=1,nnx
            dtdz     = 2.*(t(i,j,iscl,k)-tsfcc(iscl))*dz_i
            tbc(i,j,iscl,2) = t(i,j,iscl,k) - dtdz*dzu(k)
         enddo
         enddo
      enddo
c
c ------------ initialize u, v, w, t derivatives are
c              computed in get_derv at km1
c              no need to call derivatives here since
c              wbc = 0, change for more general lower bc
c
      do j=iys,iye
      do i=1,nnx
         u(i,j,km1)  = ubc(i,j,2)
         v(i,j,km1)  = vbc(i,j,2)
         w(i,j,km1)  = wbc(i,j,2)
         r3(i,j,km1) =  0.0
         e(i,j,km1)  = ebc(i,j,2)
         wf(i,j,km1) = wfbc(i,j,2)
         uf(i,j,km1) = ufbc(i,j,2)
         vf(i,j,km1) = vfbc(i,j,2)
      enddo
      enddo
c
      do iscl=1,nscl
         do j=iys,iye
         do i=1,nnx
            t1tbc(i,j,iscl) = 0.0
            t2tbc(i,j,iscl) = 0.0
            t3tbc(i,j,iscl) = taut3m(i,j,iscl)
            t(i,j,iscl,km1) = tbc(i,j,iscl,2)
         enddo
         enddo
      enddo
c
      return
      end
      subroutine lower_free(it)
c
c --------------- setup lower boundary condition for free
c                 convection where each processor applies
c                 log-law at several (ix,iy) for k = 1.
c                 index f(.,.,2)  indicates lower
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
c
      real u_level1(nnx,iys:iye,2+nscl), buf(2+2*nscl)
      real sbuf(2+2*nscl,mxs:mxe,iys:iye)
      real rbuf((2+2*nscl)*nnx*(iye+1-iys))
      integer istatus(mpi_status_size)
c
c -------------- broadcast level 1 data everywhere
c
      if(iss .eq. 0) then
         do iy=iys,iye
         do ix=1,nnx
            u_level1(ix,iy,1) = u(ix,iy,1)
            u_level1(ix,iy,2) = v(ix,iy,1)
         enddo
         enddo
         do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            u_level1(ix,iy,2+iscl) = t(ix,iy,iscl,1)
         enddo
         enddo
         enddo
      endif
      num = nnx*(iye + 1 - iys)*(2+nscl)
c
c ------ send all of root data to other processors
c
      call mpi_send_root(u_level1(1,iys,1),
     +             num,myid,numprocs,ncpu_s)
c
c --------- every task gets their own fluxes and surface scalars
c
      call suft2(u_level1,it)
c
c --------- send surface scalars and momentum fluxes
c           back to root(s)
c
      if(numprocs .eq. 1) go to 999
c
      do iy=iys,iye
      do ix=mxs,mxe
         sbuf(1,ix,iy)  = tau13m(ix,iy)
         sbuf(2,ix,iy)  = tau23m(ix,iy)
      enddo
      enddo
      do iscl=1,nscl
      do iy=iys,iye
      do ix=mxs,mxe
         sbuf(2+iscl,ix,iy)      = taut3m(ix,iy,iscl)
         sbuf(2+nscl+iscl,ix,iy) = t_grnd(ix,iy,iscl)
      enddo
      enddo
      enddo
c
      irow_r = mod(myid,ncpu_s)
      if(myid .ge. ncpu_s) then
        num = (2+2*nscl)*(mxe+1-mxs)*(iye+1-iys)
        call mpi_send(sbuf(1,mxs,iys),num,mpi_real8,irow_r,1,
     +       mpi_comm_world,ierr)
      else
        do l=irow_r+ncpu_s,numprocs-1,ncpu_s
           num = (2+2*nscl)*(mx_e(l)+1-mx_s(l))*(iye+1-iys)
           call mpi_recv(rbuf(1),num,mpi_real8,l,1,
     +          mpi_comm_world,istatus,ierr)
c          call f_suft2(rbuf,maxnx,maxny,mx_s(l),mx_e(l),iys,iye,nscl,
           call f_suft2(rbuf,nnx,mx_s(l),mx_e(l),iys,iye,nscl,
     +                  tau13m,tau23m,taut3m,t_grnd)
        enddo
      endif
c
  999 continue
c
c ------------ only for root row = 0
c              get sums of surface conditions
c              and set surface boundary conditions
c
      if(iss .eq. 0) then
c
         buf(1) = 0.0
         buf(2) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            buf(1) = buf(1) + tau13m(ix,iy)
            buf(2) = buf(2) + tau23m(ix,iy)
         enddo
         enddo
         do iscl=1,nscl
            buf(2+iscl)      = 0.
            buf(2+nscl+iscl) = 0.
            do iy=iys,iye
            do ix=1,nnx
               buf(2+iscl)      = buf(2+iscl) + taut3m(ix,iy,iscl)
               buf(2+nscl+iscl) = buf(2+nscl+iscl) + t_grnd(ix,iy,iscl)
            enddo
            enddo
         enddo
c
         call mpi_sum_xy(buf,myid,iss,ise,2+2*nscl)
         uwsfc = buf(1)*fnxy
         vwsfc = buf(2)*fnxy
         do iscl=1,nscl
            wtsfc(iscl) = buf(2+iscl)*fnxy
            tsfcc(iscl) = buf(2+nscl+iscl)*fnxy
         enddo
c
         do j=iys,iye
         do i=1,nnx
            t11bc(i,j) = 0.0
            t12bc(i,j) = 0.0
            t13bc(i,j) = tau13m(i,j)
            t22bc(i,j) = 0.0
            t23bc(i,j) = tau23m(i,j)
            t33bc(i,j) = 0.0
            t1tbc(i,j,1) = 0.0
            t2tbc(i,j,1) = 0.0
            t3tbc(i,j,1) = taut3m(i,j,1)
         enddo
         enddo
c
         iz   = 1
         izm1 = iz - 1
         dz_i = dzu_i(iz)
c
         do iy=iys,iye
         do ix=1,nnx
            ebc(ix,iy,2)=amax1(e(ix,iy,iz),sml_eg)
            wbc(ix,iy,2)= 0.0
            pbc(ix,iy,2) = 0.0
            pbc2(ix,iy,2) = 0.0
c
            wfbc(i,j,2)  = 0.0
            ufbc(i,j,2)  = 0.0
            vfbc(i,j,2)  = 0.0
         enddo
         enddo
c
         do iy=iys,iye
         do ix=1,nnx
            dudz     = 2.*u(ix,iy,iz)*dz_i
            dvdz     = 2.*v(ix,iy,iz)*dz_i
            ubc(ix,iy,2) = u(ix,iy,iz) - dudz*dzu(iz)
            vbc(ix,iy,2) = v(ix,iy,iz) - dvdz*dzu(iz)
         enddo
         enddo
         do iscl=1,nscl
            do iy=iys,iye
            do ix=1,nnx
               dtdz     = 2.*(t(ix,iy,iscl,iz)-tsfcc(iscl))*dz_i
               tbc(ix,iy,iscl,2) = t(ix,iy,iscl,iz) - dtdz*dzu(iz)
            enddo
            enddo
         enddo
c
c ------------ initialize u, v, w, t and derivatives at izm1
c
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,izm1)  = ubc(ix,iy,2)
            v(ix,iy,izm1)  = vbc(ix,iy,2)
            w(ix,iy,izm1)  = wbc(ix,iy,2)
            r3(ix,iy,izm1) =  0.0
            e(ix,iy,izm1)  = ebc(ix,iy,2)
            ux(ix,iy,izm1) = 0.0
            uy(ix,iy,izm1) = 0.0
            vx(ix,iy,izm1) = 0.0
            vy(ix,iy,izm1) = 0.0
            wx(ix,iy,izm1) = wbc(ix,iy,2)
            wy(ix,iy,izm1) = wbc(ix,iy,2)
         enddo
         enddo
c
         do iscl=1,nscl
            do iy=iys,iye
            do ix=1,nnx
               t(ix,iy,iscl,izm1) = tbc(ix,iy,iscl,2)
            enddo
            enddo
         enddo
c
c ----- end of if block for root row
c
      endif
c
 7999 continue
c
      return
      end
      subroutine f_suft2(rbuf,nnx,mxs,mxe,iys,iye,nscl,
     +                  tau13m,tau23m,taut3m,t_grnd)
c
c ------ fill surface arrays on root processors
c
      real rbuf(2+2*nscl,mxs:mxe,iys:iye)
      real tau13m(nnx,iys:iye), tau23m(nnx,iys:iye),
     +     taut3m(nnx,iys:iye,nscl), t_grnd(nnx,iys:iye,nscl)
c
      do iy=iys,iye
      do ix=mxs,mxe
         tau13m(ix,iy) = rbuf(1,ix,iy)
         tau23m(ix,iy) = rbuf(2,ix,iy)
      enddo
      enddo
      do iscl=1,nscl
         do iy=iys,iye
         do ix=mxs,mxe
            taut3m(ix,iy,iscl) = rbuf(2+iscl,ix,iy)
            t_grnd(ix,iy,iscl) = rbuf(2+nscl+iscl,ix,iy)
         enddo
         enddo
      enddo
c
      return
      end
      subroutine lower_curvy(it,istage)
c
c --------- setup lower boundary condition for case
c           with curvy lower boundary. each processor
c           applies log-law at several (i,j) for k = 1
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
c
c ------- set the number of variables to pass and
c         get back
c
      parameter (mvar= 9+3*nscl, jvar = 8+2*nscl)
      real u_level1(nnx,iys:iye,mvar), buf(jvar)
      real sbuf(jvar,mxs:mxe,iys:iye)
      real rbuf(jvar*nnx*(iye+1-iys))
      integer istatus(mpi_status_size)
c
      t_bndy = t_stage - t_zero
c
c        write(nprt, 1010) t_bndy, t_stage, t_zero
c1010    format(' 1010: lower curvy',/,
c    +          ' t_bndy = ',e15.6,
c    +          ' t_stage = ',e15.6,' t_zero = ',e15.6)
c
      if(i3d_surf .eq. 0) then
         call l_bndy(t_bndy)
      elseif(i3d_surf .eq. 1) then
         call get_surf(t_bndy)
      elseif(i3d_surf .eq. 2) then
         call get_surf_gabls(t_bndy)
      elseif(i3d_surf .eq. 3) then
         call get_surf_tank(t_bndy)
      endif
c
c     call dealias_grid_speed
c
c ----- broadcast level 1 data to processors above root
c
      k   = 1
      km1 = k - 1
c
      if(iss .eq. 0) then
         do j=iys,iye
         do i=1,nnx
c
c -------- subtract off surface movement when applying the drag law
c          plus the orbital velocity of the wave.
c          w at the boundary is close to grid speed.
c          Note for iz_space = 4 grid speed is constant across
c          the first grid cell and thus z_t at w point = z_t at
c          the boundary. Can improve the surface estimate FIX
c          and couple it to the underlying water using the
c          fuller definition
c
c           u_level1(i,j,1) = u(i,j,k)
c           u_level1(i,j,2) = v(i,j,k)
c           u_level1(i,j,3) = w(i,j,k)
c
            u_level1(i,j,1) = u(i,j,k) - u_orbit(i,j)
            u_level1(i,j,2) = v(i,j,k) - v_orbit(i,j)
c           u_level1(i,j,3) = w(i,j,k) - gm(i,j,km1,5)
            u_level1(i,j,3) = w(i,j,k) - w_orbit(i,j)
            u_level1(i,j,4) = gm(i,j,k,1)
            u_level1(i,j,5) = gm(i,j,k,2)
            u_level1(i,j,6) = gm(i,j,k,3)
            u_level1(i,j,7) = u_orbit(i,j)
            u_level1(i,j,8) = v_orbit(i,j)
            u_level1(i,j,9) = gm(i,j,km1,5)
         enddo
         enddo
c           j = iys
c           write(nprt,3002) t_bndy, dt, istage
c3002       format(' 3002 debug ',5x,/,
c    +             ' t bndy = ',e20.10,' dt = ',e15.6,' istage = ',i5)
c           write(nprt,3001) (i,bndy(i,j),gm(i,j,k,1),gm(i,j,km1,5),
c    +                        gm(i,j,k,5), w(i,j,k),u(i,j,k),
c    +                        u_orbit(i,j),u_level1(i,j,1),
c    +                        u_level1(i,j,3),i=1,nnx)
c3001  format(' i ',10x,' h ',10x,' 1/J',10x,' gm5 ',10x,' gm5',10x,
c    +        ' w ',10x,' u ',10x,' u orb ',10x,' du ',5x,' dw ',/,
c    +        (i5,9e15.6))
c
c          write(nprt,3003) (i,gm(i,j,km1,5),wf(i,j,km1), w(i,j,km1),
c    +                          gm(i,j,k,5), wf(i,j,k), w(i,j,k),
c    +                          gm(i,j,k+1,5), wf(i,j,k+1),w(i,j,k+1),
c    +                          i=1,nnx)
c3003  format(' 3003 debug',/,
c    +        ' i ',10x,' gm 5 0',10x,' wf 0 ',10x,' w 0',8x,
c    +                  ' gm5 1',10x,' wf 1 ',8x,' w 1 ',8x,
c    +        ' gm 5 2',10x,' wf 2 ',10x,' w 2',/,
c    +        (i5,9e15.6))
         do iscl=1,nscl
         do j=iys,iye
         do i=1,nnx
            u_level1(i,j,9+iscl)        = t3tbc(i,j,iscl)
            u_level1(i,j,9+nscl+iscl)   = t(i,j,iscl,k)
            u_level1(i,j,9+2*nscl+iscl) = t_grnd(i,j,iscl)
         enddo
         enddo
         enddo
      endif
c
      num = nnx*(iye + 1 - iys)*mvar
c
c ------ send all of root data to other processors
c
      call mpi_send_root(u_level1(1,iys,1),
     +             num,myid,numprocs,ncpu_s)
c
c --------- every task gets their own fluxes and surface scalars
c
      do j=iys,iye
      do i=1,nnx
         ustar(i,j)    = -999999.99
c        t_grnd(i,j,1) = -999999.99
      enddo
      enddo
c

      call drglaw_l(u_level1,mvar,it)
c
c     j = iys
c     write(nprt,8201) (i,t11bc(i,j),t12bc(i,j),t13bc(i,j),
c    +                  t22bc(i,j),t23bc(i,j),t33bc(i,j),
c    +                  ustar(i,j),i=1,nnx)
c8201 format(' 8201 debug tau',/,
c    +       ' i ',5x,' t11 ',5x,' t12 ',5x,' t13 ',5x,' t22',5x,
c    +       ' t23 ',5x,' t33',5x,' ustar ',/,(i5,7e15.6))
c
c --------- send surface scalars and momentum fluxes
c           back to root(s)
c
      if(numprocs .eq. 1) go to 999
c
      do j=iys,iye
      do i=mxs,mxe
         sbuf(1,i,j)  = t11bc(i,j)
         sbuf(2,i,j)  = t12bc(i,j)
         sbuf(3,i,j)  = t13bc(i,j)
         sbuf(4,i,j)  = t22bc(i,j)
         sbuf(5,i,j)  = t23bc(i,j)
         sbuf(6,i,j)  = t33bc(i,j)
         sbuf(7,i,j)  = ustar(i,j)
         sbuf(8,i,j)  = amol(i,j)
      enddo
      enddo
      do iscl=1,nscl
      do j=iys,iye
      do i=mxs,mxe
         sbuf(8+iscl,i,j)      = t3tbc(i,j,iscl)
         sbuf(8+nscl+iscl,i,j) = t_grnd(i,j,iscl)
      enddo
      enddo
      enddo
c
      irow_r = mod(myid,ncpu_s)
      if(myid .ge. ncpu_s) then
        num = jvar*(mxe+1-mxs)*(iye+1-iys)
        call mpi_send(sbuf(1,mxs,iys),num,mpi_real8,irow_r,1,
     +       mpi_comm_world,ierr)
      else
        do l=irow_r+ncpu_s,numprocs-1,ncpu_s
           num = jvar*(mx_e(l)+1-mx_s(l))*(iye+1-iys)
           call mpi_recv(rbuf(1),num,mpi_real8,l,1,
     +          mpi_comm_world,istatus,ierr)
           call f_suft_curvy(rbuf,nnx,mx_s(l),mx_e(l),
     +          iys,iye,nscl,jvar,t11bc,t12bc,t13bc,t22bc,
     +          t23bc,t33bc,ustar,t3tbc,t_grnd,amol)
        enddo
      endif
c
  999 continue
c
      if(iss .eq. 0) then
c
c     j = iys
c     write(nprt,7201) (i,t11bc(i,j),t12bc(i,j),t13bc(i,j),
c    +                  t22bc(i,j),t23bc(i,j),t33bc(i,j),
c    +                  ustar(i,j),i=1,nnx)
c7201 format(' 7201 debug tau',/,
c    +   ' i ',5x,' t11 ',5x,' t12 ',5x,' t13 ',5x,' t22',5x,
c    +   ' t23 ',5x,' t33',5x,' ustar ',/,(i5,7e15.6))
c     write(nprt,8202) (i, -gm(i,j,k,2)*gm(i,j,k,1),
c    +                  u(i,j,k), u_orbit(i,j),
c    +                  w(i,j,k), w_orbit(i,j),
c    +                  (u(i,j,k) - u_orbit(i,j)),
c    +                  (w(i,j,k) - w_orbit(i,j)),
c    +                  t13bc(i,j), ustar(i,j)**2,z1b(i,j),
c    +                  i=1,nnx)
c8202 format(' 8202 debug surf',/,
c    +       ' i ',5x,' gm(2)',5x, 'u ',5x,' uorbit',
c    +       ' w ',5x,' worbit',5x,' du ',5x,' dw',5x,
c    +       ' t13bc ',5x,' u**2 ',5x,' z1b',/,
c    +       (i5,10e15.6))
c
c ------------ only for root row = 0
c              get sums of surface conditions
c              and set surface boundary conditions
c
c
         do i=1,jvar
            buf(i) = 0.0
         enddo
c
c ------- add SGS pieces that appear on rhs of [u,v,w] equations
c
         do j=iys,iye
         do i=1,nnx
            fac1   = gm(i,j,k,1)*gm(i,j,k,2)
            fac2   = gm(i,j,k,1)*gm(i,j,k,3)
            buf(1) = buf(1) + t11bc(i,j)
            buf(2) = buf(2) + t12bc(i,j)
            buf(3) = buf(3) + t13bc(i,j) +
     +               fac1*t11bc(i,j) + fac2*t12bc(i,j)
            buf(4) = buf(4) + t22bc(i,j)
            buf(5) = buf(5) + t23bc(i,j) +
     +               fac1*t12bc(i,j) + fac2*t22bc(i,j)
            buf(6) = buf(6) + t33bc(i,j) +
     +               fac1*t13bc(i,j) + fac2*t23bc(i,j)
            buf(7) = buf(7) + ustar(i,j)**2
            buf(8) = buf(8) + amol(i,j)
         enddo
         enddo
         do iscl=1,nscl
            do j=iys,iye
            do i=1,nnx
               buf(8+iscl)      = buf(8+iscl) + t3tbc(i,j,iscl)
               buf(8+nscl+iscl) = buf(8+nscl+iscl) +
     +                            t_grnd(i,j,iscl)
            enddo
            enddo
         enddo
c
         call mpi_sum_xy(buf,myid,iss,ise,jvar)
         uusfc   =  buf(1)*fnxy
         uvsfc   =  buf(2)*fnxy
         uwsfc   =  buf(3)*fnxy
         vvsfc   =  buf(4)*fnxy
         vwsfc   =  buf(5)*fnxy
         wwsfc   =  buf(6)*fnxy
         utau    =  sqrt(buf(7)*fnxy)
         amonin  =  buf(8)*fnxy
         utausv  = utau
         do iscl=1,nscl
            wtsfc(iscl) = buf(8+iscl)*fnxy
            tsfcc(iscl) = buf(8+nscl+iscl)*fnxy
         enddo
         if(ibuoy == 0 .or. wtsfc(1) == 0.) then
            zol = 0.0
         else
            zol = zz(1)/amonin
         endif
c
         dz_i = dzu_i(k)
c
         do j=iys,iye
         do i=1,nnx
            !ebc(i,j,2)  = amax1(e(i,j,k),sml_eg)
            ebc(i,j,2) = 0.0
            pbc(i,j,2)  = 0.0
            pbc2(i,j,2) = 0.0
c
            ubc(i,j,2)  = 2.0*u_orbit(i,j) - u(i,j,k)
            vbc(i,j,2)  = 2.0*v_orbit(i,j) - v(i,j,k)
c
Cnew ------ bc for vertical velocity
c
            wbc(i,j,2)  = 2.0*gm(i,j,km1,5) - w(i,j,k)
c
            ufbc(i,j,2)  = -uf(i,j,k)
            vfbc(i,j,2)  = -vf(i,j,k)
c
c ------- no flow through the surface at t = t_stage
c
            wfbc(i,j,2)  = gm(i,j,km1,5)
         enddo
         enddo
c
         do iscl=1,nscl
            do j=iys,iye
            do i=1,nnx
               tbc(i,j,iscl,2) = 2.0*t_grnd(i,j,iscl) - t(i,j,iscl,k)
            enddo
            enddo
         enddo
c
c ------------ initialize u, v, w, t and derivatives at km1
c
         do j=iys,iye
         do i=1,nnx
            u(i,j,km1)  = ubc(i,j,2)
            v(i,j,km1)  = vbc(i,j,2)
            w(i,j,km1)  = wbc(i,j,2)
            r3(i,j,km1) =  0.0
            e(i,j,km1)  = ebc(i,j,2)
            ux(i,j,km1) = 0.0
            uy(i,j,km1) = 0.0
            vx(i,j,km1) = 0.0
            vy(i,j,km1) = 0.0
            wx(i,j,km1) = wbc(i,j,2)
            wy(i,j,km1) = wbc(i,j,2)
c
            wf(i,j,km1) = wfbc(i,j,2)
            uf(i,j,km1) = ufbc(i,j,2)
            vf(i,j,km1) = vfbc(i,j,2)
         enddo
         enddo
c
         do iscl=1,nscl
            do j=iys,iye
            do i=1,nnx
               t(i,j,iscl,km1) = tbc(i,j,iscl,2)
            enddo
            enddo
         enddo
c
      endif ! end of if block for root row
c
      return
      end
      subroutine lower_curvy_dns(it,istage)
c
c --------- setup lower boundary condition for DNS
c           with curvy lower boundary.
c           This is no-slip for a moving wavy surface
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
c
c ------- set the number of variables to pass and
c         get back
c
      parameter (mvar= 9+3*nscl, jvar = 8+2*nscl)
      real u_level1(nnx,iys:iye,mvar), buf(jvar)
      real sbuf(jvar,mxs:mxe,iys:iye)
      real rbuf(jvar*nnx*(iye+1-iys))
      integer istatus(mpi_status_size)
c
      t_bndy = t_stage - t_zero
c
c
      if(i3d_surf .eq. 0) then
         call l_bndy(t_bndy)
      elseif(i3d_surf .eq. 1) then
         call get_surf(t_bndy)
      elseif(i3d_surf .eq. 2) then
         call get_surf_gabls(t_bndy)
      elseif(i3d_surf .eq. 3) then
         call get_surf_tank(t_bndy)
      endif
c
c ----- broadcast level 1 data to processors above root
c
      k   = 1
      km1 = k - 1
c
c
      if(iss .eq. 0) then
c
c ------------ only for root row = 0
c              get sums of surface conditions
c              and set surface boundary conditions
c
c
         do i=1,jvar
            buf(i) = 0.0
         enddo
c
c ------- add SGS pieces that appear on rhs of [u,v,w] equations
c ------- DHR (3/17/15) I'm not sure if I should leave this, but I am
c
!         do j=iys,iye
!         do i=1,nnx
!            fac1   = gm(i,j,k,1)*gm(i,j,k,2)
!            fac2   = gm(i,j,k,1)*gm(i,j,k,3)
!            buf(1) = buf(1) + t11bc(i,j)
!            buf(2) = buf(2) + t12bc(i,j)
!            buf(3) = buf(3) + t13bc(i,j) +
!     +               fac1*t11bc(i,j) + fac2*t12bc(i,j)
!            buf(4) = buf(4) + t22bc(i,j)
!            buf(5) = buf(5) + t23bc(i,j) +
!     +               fac1*t12bc(i,j) + fac2*t22bc(i,j)
!            buf(6) = buf(6) + t33bc(i,j) +
!     +               fac1*t13bc(i,j) + fac2*t23bc(i,j)
!            buf(7) = buf(7) + ustar(i,j)**2
!            buf(8) = buf(8) + amol(i,j)
!         enddo
!         enddo
!         do iscl=1,nscl
!            do j=iys,iye
!            do i=1,nnx
!               buf(8+iscl)      = buf(8+iscl) + t3tbc(i,j,iscl)
!               buf(8+nscl+iscl) = buf(8+nscl+iscl) +
!     +                            t_grnd(i,j,iscl)
!            enddo
!            enddo
!         enddo
!c
!         call mpi_sum_xy(buf,myid,iss,ise,jvar)
!         uusfc   =  buf(1)*fnxy
!         uvsfc   =  buf(2)*fnxy
!         uwsfc   =  buf(3)*fnxy
!         vvsfc   =  buf(4)*fnxy
!         vwsfc   =  buf(5)*fnxy
!         wwsfc   =  buf(6)*fnxy
!         utau    =  sqrt(buf(7)*fnxy)
!         amonin  =  buf(8)*fnxy
!         utausv  = utau
!         do iscl=1,nscl
!            wtsfc(iscl) = buf(8+iscl)*fnxy
!            tsfcc(iscl) = buf(8+nscl+iscl)*fnxy
!         enddo
!         if(ibuoy == 0 .or. wtsfc(1) == 0.) then
!            zol = 0.0
!         else
!            zol = zz(1)/amonin
!         endif
c
         dz_i = dzu_i(k)
c
         do j=iys,iye
         do i=1,nnx
            ebc(i,j,2)  = amax1(e(i,j,k),sml_eg)
            pbc(i,j,2)  = 0.0
            pbc2(i,j,2) = 0.0
c
            ubc(i,j,2)  = 2.0*u_orbit(i,j) - u(i,j,k)
            vbc(i,j,2)  = 2.0*v_orbit(i,j) - v(i,j,k)
c
Cnew ------ bc for vertical velocity
c
            wbc(i,j,2)  = 2.0*gm(i,j,km1,5) - w(i,j,k)
c           wbc(i,j,2)  = -w(i,j,k)
c
            ufbc(i,j,2)  = -uf(i,j,k)
            vfbc(i,j,2)  = -vf(i,j,k)
c
c ------- no flow through the surface at t = t_stage
c
            wfbc(i,j,2)  = gm(i,j,km1,5)
         enddo
         enddo
c
         do iscl=1,nscl
            do j=iys,iye
            do i=1,nnx
               if (isfc(iscl)==1) then
                  tbc(i,j,iscl,2) = 2.0*Tbot(iscl) - t(i,j,iscl,k)
               end if
               if (isfc(iscl)==0) then
                  tbc(i,j,iscl,2) = t(i,j,iscl,k) +
     +           dzu(0)*wtsfc(1)/vis_s(i,j,k)
               end if
               !tbc(i,j,iscl,2) = 2.0*t_grnd(i,j,iscl) - t(i,j,iscl,k)
            enddo
            enddo
         enddo
c
c ------------ initialize u, v, w, t and derivatives at km1
c
         do j=iys,iye
         do i=1,nnx
            u(i,j,km1)  = ubc(i,j,2)
            v(i,j,km1)  = vbc(i,j,2)
            w(i,j,km1)  = wbc(i,j,2)
            r3(i,j,km1) =  0.0
            e(i,j,km1)  = ebc(i,j,2)
            ux(i,j,km1) = 0.0
            uy(i,j,km1) = 0.0
            vx(i,j,km1) = 0.0
            vy(i,j,km1) = 0.0
            wx(i,j,km1) = wbc(i,j,2)
            wy(i,j,km1) = wbc(i,j,2)
c
            wf(i,j,km1) = wfbc(i,j,2)
            uf(i,j,km1) = ufbc(i,j,2)
            vf(i,j,km1) = vfbc(i,j,2)
         enddo
         enddo
c
         do iscl=1,nscl
            do j=iys,iye
            do i=1,nnx
               t(i,j,iscl,km1) = tbc(i,j,iscl,2)
            enddo
            enddo
         enddo
c
      endif ! end of if block for root row
c
      return
      end
      subroutine f_suft_curvy(rbuf,nnx,mxs,mxe,iys,iye,nscl,jvar,
     +     t11bc,t12bc,t13bc,t22bc,t23bc,t33bc,ustar,t3tbc,t_grnd,
     +     amol)
c
c ------ fill surface arrays on root processors
c
      real rbuf(jvar,mxs:mxe,iys:iye)
      real t11bc(nnx,iys:iye), t12bc(nnx,iys:iye),
     +     t13bc(nnx,iys:iye), t22bc(nnx,iys:iye),
     +     t23bc(nnx,iys:iye), t33bc(nnx,iys:iye),
     +     ustar(nnx,iys:iye), amol(nnx,iys:iye),
     +     t3tbc(nnx,iys:iye,nscl), t_grnd(nnx,iys:iye,nscl)
c
      do j=iys,iye
      do i=mxs,mxe
         t11bc(i,j)  =  rbuf(1,i,j)
         t12bc(i,j)  =  rbuf(2,i,j)
         t13bc(i,j)  =  rbuf(3,i,j)
         t22bc(i,j)  =  rbuf(4,i,j)
         t23bc(i,j)  =  rbuf(5,i,j)
         t33bc(i,j)  =  rbuf(6,i,j)
         ustar(i,j)  =  rbuf(7,i,j)
         amol(i,j)   = rbuf(8,i,j)
      enddo
      enddo
      do iscl=1,nscl
         do j=iys,iye
         do i=mxs,mxe
            t3tbc(i,j,iscl)  = rbuf(8+iscl,i,j)
            t_grnd(i,j,iscl) = rbuf(8+nscl+iscl,i,j)
         enddo
         enddo
      enddo
c
      return
      end
      subroutine drglaw_l(u_level1,mvar,it)
c
c --------- apply drag law across a set of gridpoints
c           with curvy lower boundary. each processor
c           applies log-law at several (i,j) for k = 1
c
c     u_level1(.,.,1)            =  u
c     u_level1(.,.,2)            =  v
c     u_level1(.,.,3)            =  w
c     u_level1(.,.,4)            =  gm(.,.,1,1)
c     u_level1(.,.,5)            =  gm(.,.,1,2)
c     u_level1(.,.,6)            =  gm(.,.,1,3)
c     u_level1(.,.,7)            =  u_orbit
c     u_level1(.,.,8)            =  v_orbit
c     u_level1(.,.,9)            =  surface grid speed
c     u_level1(.,.,9+iscl)       =  t3tbc(.,.,iscl)
c     u_level1(.,.,9+nscl+iscl)  =  t(.,.,iscl,1)
c     u_level1(.,.,9+2nscl+iscl) =  t_grnd(.,.,iscl)
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      real u_level1(nnx,iys:iye,mvar)
      real aij(3,3), us(3), usn(3)
      real tlv1(nscl), qflux(nscl), tsurf(nscl)
c
c ------- get matrix of direction cosines for surface coordinate system
c         and compute surface fluxes for a set of gridpoints [mxs:mxe,iys:iye]
c
      do j=iys,iye
      do i=mxs,mxe
c
c -------- vector aligned with \xi (t1)
c
         e1   = 1.0
         e3   = -u_level1(i,j,4)*u_level1(i,j,5)
         emag = sqrt(e1*e1 + e3*e3)
         e1   = e1/emag
         e3   = e3/emag
c
c -------- vector aligned with \eta (t2)
c
         f2   = 1.0
         f3   = -u_level1(i,j,4)*u_level1(i,j,6)
         fmag = sqrt(f2*f2 + f3*f3)
         f2   = f2/fmag
         f3   = f3/fmag
c
c -------- surface normal n
c
         g1   = -e3*f2
         g2   = -e1*f3
         g3   = e1*f2
         gmag = sqrt(g1*g1 + g2*g2 + g3*g3)
         g1   = g1/gmag
         g2   = g2/gmag
         g3   = g3/gmag
c
c -------- surface wind (2 vectors, us_1 and us_2)
c
         us1a      = u_level1(i,j,1)*e1 + u_level1(i,j,3)*e3
         us2a      = u_level1(i,j,2)*f2 + u_level1(i,j,3)*f3
         us(1)     = us1a*e1
         us(2)     = us2a*f2
         us(3)     = us1a*e3 + us2a*f3
         wind(i,j) = sqrt(us(1)**2 + us(2)**2 + us(3)**2)
         usn(1)    = us(1)/wind(i,j)
         usn(2)    = us(2)/wind(i,j)
         usn(3)    = us(3)/wind(i,j)
c        d_wind    = sqrt((us(1) - u_level1(i,j,7))**2 +
c    +                    (us(2) - u_level1(i,j,8))**2 +
c    +                    (us(3) - u_level1(i,j,9))**2)
c
c -------- third dimension r = n cross us
c
         rb1   = g2*usn(3) - g3*usn(2)
         rb2   = g3*usn(1) - g1*usn(3)
         rb3   = g1*usn(2) - g2*usn(1)
         rbmag = sqrt(rb1**2 + rb2**2 + rb3**2)
         rb1   = rb1/rbmag
         rb2   = rb2/rbmag
         rb3   = rb3/rbmag
c
c -------- matrix of direction cosines
c
         aij(1,1) = usn(1)
         aij(1,2) = usn(2)
         aij(1,3) = usn(3)
         aij(2,1) = rb1
         aij(2,2) = rb2
         aij(2,3) = rb3
         aij(3,1) = g1
         aij(3,2) = g2
         aij(3,3) = g3
c
c --------- set simple scalar arrays
c           fix here if specified surface value
c
         do l=1,nscl
            qflux(l) = u_level1(i,j,9+l)
            tlv1(l)  = u_level1(i,j,9+nscl+l)
            tsurf(l) = u_level1(i,j,9+2*nscl+l)
         enddo
c
         zody = alog(z1b(i,j)/zo + 1.0)
         call mo(ibuoy,isfc,ismlt,vk,batag,zody,z1b(i,j),
     +           wind(i,j), ustar(i,j), nscl,
     +           tlv1(1), tsurf(1), qflux(1),
     +           zol, cd, ch, amol(i,j),RHbot)
c
c ------------ build the momentum fluxes
c
         utau2      = -ustar(i,j)**2
         t11bc(i,j) = utau2*2.0*aij(1,1)*aij(3,1)
         t12bc(i,j) = utau2*(aij(1,1)*aij(3,2) +
     +                        aij(3,1)*aij(1,2))
         t13bc(i,j) = utau2*(aij(1,1)*aij(3,3) +
     +                        aij(1,3)*aij(3,1))
         t22bc(i,j) = utau2*2.0*aij(1,2)*aij(3,2)
         t23bc(i,j) = utau2*(aij(1,2)*aij(3,3) +
     +                        aij(1,3)*aij(3,2))
         t33bc(i,j) = utau2*2.0*aij(1,3)*aij(3,3)
c
         do iscl=1,nscl
            t1tbc(i,j,iscl)  = 0.0
            t2tbc(i,j,iscl)  = 0.0
            t3tbc(i,j,iscl)  = qflux(iscl)
            t_grnd(i,j,iscl) = tsurf(iscl)
         enddo
      enddo
      enddo
c
      return
      end
      subroutine mo(ibuoy,isfc,ismlt,vk,batag,zody,z1,
     +              wind,utau,nscl,t1xy,tsfcc,qstar,
     +              zeta,cd,ch,amonin,RHbot)
      use particles
c
c ---------- routine to find surface fluxes and values
c            using MO theory. Iterate for z/L
c
c            Gridpoints [mxs:mxe,iys:iye].
c
c    [assume:  z1 is constant, zody is constant, minimum wind ]
c
c            isfc = 1, then compute tstar and qstar given tsfcc
c                 = 0,   "     "    tsfcc given qstar
c
c    can add more scalars using definition below
c
      parameter (iter_mo = 30, zeta_min = -6.0, zeta_max = 3.0)
      parameter (wind_min = 0.00001)
      real t1xy(nscl), tsfcc(nscl), qstar(nscl)
      integer isfc(nscl)
      real botC,RHB,RHbot
c

      if(ismlt .eq. 1) then
          vk74in = 0.74/vk
      else
          vk74in = 1.0/vk
      endif
      batagk = batag*vk
c
      windm = amax1(wind_min,wind)
c
c ---------- limits for zeta
c
      zeta_mn = zeta_min
      zeta_mx = zeta_max
      if(isfc(1) .eq. 0) then
         f_con = z1*batagk*qstar(1)/((windm*vk)**3)
      else
         d_theta = vk74in*(tsfcc(1) - t1xy(1))
         f_con   = z1*batagk*d_theta/((windm*vk)**2)
      endif
c
c --------- iteration for zeta
c
      do iter=1,iter_mo
         zeta_a = 0.5*(zeta_mn + zeta_mx)
         if(ismlt .eq. 1) then
             call busngr(zeta_a,phim,phis,psim,psis)
         else
             call fzol(zeta_a,phim,phis,psim,psis)
         endif
         u_fac = (zody - psim)
         if(isfc(1) .eq. 0) then
            f_new =  zeta_a + f_con*u_fac**3
         else
            t_fac = 1.0/amax1((zody - psis),0.01)
            f_new =  zeta_a + f_con*u_fac*u_fac*t_fac
         endif
         if(f_new .lt. 0.0) then
            zeta_mn = zeta_a
         else
            zeta_mx = zeta_a
         endif
      enddo
c
c --------- check if neutral surface layer
c
      if(ibuoy == 0 .or. qstar(1) == 0.) then
          amonin   = 1000.
          zeta     = 0.0
          utau     = windm*vk/zody
          tsfcc(1) = t1xy(1)
          dnom = zody*vk74in
      else
         dnom      = amax1((zody-psis),0.01)*vk74in
         utau      = windm*vk/(zody-psim)
         if(isfc(1) == 0) then
            thstar = -qstar(1)/utau
            tsfcc(1)  = t1xy(1) - thstar*dnom
         else
            thstar = (t1xy(1) - tsfcc(1))/dnom
            qstar(1)  = -utau*thstar
         endif
         amonin = -utau**3/(batagk*qstar(1))
         zeta   = z1/amonin
      endif
      cd = utau**2/windm**2
      ch = 0.0
      dnom_ch = windm*(tsfcc(1) - t1xy(1))
      if(dnom_ch .ne. 0.) then
        ch = qstar(1)/dnom_ch
      endif

      !Do the lower BC for RH:

      if (isfc(2) == 0) then  !Set surface flux, solve for tsfcc
         thstar = -qstar(2)/utau
         tsfcc(2) = t1xy(2) - thstar*dnom
      elseif (isfc(2) == 1) then !Set tsfcc, solve for qstar
         botC=(tsfcc(1)-273.15)
         RHB = RHbot/100.0
         tsfcc(2) =  RHB*Mw/Ru/tsfcc(1)*610.94*
     +               EXP(17.6257*botC/(botC+243.04))/rhoa

         thstar = (t1xy(2) - tsfcc(2))/dnom
         qstar(2) = -utau*thstar
      end if

c
c ------ example of other scalars (specified surface value and source)
c
c     qstar(2) = qstar(1)*(t1xy(2) - tsfcc(2))/(t1xy(1) - tsfcc(1))
c
c     tsfcc(2) = t1xy(2) - qstar(2)*(t1xy(1) - tsfcc(1))/qstar(1)

c
      return
      end
      subroutine upper
c
c ---- set boundary condition on upper boundary iz=nnz
c      option for special radiation boundary condition
c                 index f(.,.,1)  indicates upper.
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      iz   = nnz
      izm1 = iz - 1
      izm2 = iz - 2
      izp1 = iz + 1
      izp2 = iz + 2
c
c     call dealias_grid_speed
c
      if(ibcu .eq. -1) then ! channel flow type simulation
c
c --------- boundary conditions are gradient conditions
c
c            dudzbc = 0.0
c            dvdzbc = 0.0
c            dtdzbc = 0.0
c            wbc    = 0.0
c            dedzbc = 0.0
c
        do iy=iys,iye
        do ix=1,nnx
c
Cnew ------ bc for vertical velocity
c
           wbc(ix,iy,1) = 2.0*gm(ix,iy,iz,5) - w(ix,iy,iz)
c          wbc(ix,iy,1) = w(ix,iy,iz)
c          ebc(ix,iy,1) = 0.0
           ebc(ix,iy,1) = e(ix,iy,iz)
           ubc(ix,iy,1) = u(ix,iy,iz)
           vbc(ix,iy,1) = v(ix,iy,iz)
           pbc(ix,iy,1) = 0.0
           pbc2(ix,iy,1)= 0.0
           ufbc(ix,iy,1) = uf(ix,iy,iz)
           vfbc(ix,iy,1) = vf(ix,iy,iz)
c
c --------- upper boundary = grid speed
c
           wfbc(ix,iy,1) = gm(ix,iy,iz,5)
        enddo
        enddo
        do iscl=1,nscl
c
c ---------- average scalar gradient
c
           dtdzf(iscl) = 0.0
        enddo
c
        call mpi_sum_xy(dtdzf,myid,iss,ise,nscl)
c
        do iscl=1,nscl
           do iy=iys,iye
           do ix=1,nnx
              tbc(ix,iy,iscl,1) = t(ix,iy,iscl,iz) +
     +                            dtdzf(iscl)*dzu(nnzp1)
           enddo
           enddo
        enddo
c
      else if(ibcu .eq. 0) then
c
c --------- boundary conditions are gradient conditions
c
c            dudzbc = 0.0
c            dvdzbc = 0.0
c            dtdzbc = dtdzf
c            wbc    = 0.0
c            ebc    = 0.0
c
        do iy=iys,iye
        do ix=1,nnx
           wbc(ix,iy,1) = 2.0*gm(ix,iy,iz,5) - w(ix,iy,iz)
c          wbc(ix,iy,1) = w(ix,iy,iz)
           ebc(ix,iy,1) = 0.0
c          ebc(ix,iy,1) = e(ix,iy,iz)
           ubc(ix,iy,1) = u(ix,iy,iz)
           vbc(ix,iy,1) = v(ix,iy,iz)
           pbc(ix,iy,1) = 0.0
           pbc2(ix,iy,1)= 0.0
           ufbc(ix,iy,1) = uf(ix,iy,iz)
           vfbc(ix,iy,1) = vf(ix,iy,iz)
c
c --------- upper boundary = grid speed
c
           wfbc(ix,iy,1) = gm(ix,iy,iz,5)
        enddo
        enddo
        do iscl=1,nscl
c
c ---------- first get average scalar gradient
c
           dtdzf(iscl) = 0.0
           do iy=iys,iye
           do ix=1,nnx
              dtdzf(iscl) = dtdzf(iscl) + (t(ix,iy,iscl,nnz) -
     +                      t(ix,iy,iscl,nnz-1))*dzu_i(nnz)
           enddo
           enddo
           dtdzf(iscl) = dtdzf(iscl)*fnxy
        enddo
c
        call mpi_sum_xy(dtdzf,myid,iss,ise,nscl)
c
        do iscl=1,nscl
           do iy=iys,iye
           do ix=1,nnx
              tbc(ix,iy,iscl,1) = t(ix,iy,iscl,iz) +
     +                            dtdzf(iscl)*dzu(nnzp1)
           enddo
           enddo
        enddo
      else if(ibcu .eq. 2) then
c
c --------- boundary conditions are Sullivan et al. (2000) Couette flow
c
        do iy=iys,iye
        do ix=1,nnx
           wbc(ix,iy,1) = 0.0
           ebc(ix,iy,1) = 0.0
           ubc(ix,iy,1) = 2.0*Uo-u(ix,iy,iz)
           vbc(ix,iy,1) = -v(ix,iy,iz)
           pbc(ix,iy,1) = 0.0
           pbc2(ix,iy,1)= 0.0
           ufbc(ix,iy,1) = Uo
           vfbc(ix,iy,1) = 0.0
c
c --------- upper boundary = grid speed
c
           wfbc(ix,iy,1) = gm(ix,iy,iz,5)
        enddo
        enddo
        do iscl=1,nscl
c
c ---------- first get average scalar gradient
c
           dtdzf(iscl) = 0.0
           do iy=iys,iye
           do ix=1,nnx
              dtdzf(iscl) = dtdzf(iscl) + (t(ix,iy,iscl,nnz) -
     +                      t(ix,iy,iscl,nnz-1))*dzu_i(nnz)
           enddo
           enddo
           dtdzf(iscl) = dtdzf(iscl)*fnxy
        enddo
c
        call mpi_sum_xy(dtdzf,myid,iss,ise,nscl)
c
        do iscl=1,nscl
           do iy=iys,iye
           do ix=1,nnx
              tbc(ix,iy,iscl,1) = t(ix,iy,iscl,iz) +
     +                            dtdzf(iscl)*dzu(nnzp1)
           enddo
           enddo
        enddo
      else if(ibcu .eq. 1) then
c
c ------------- special if iradup boundary condition
c               get estimate of w from continuity and
c               linearized relation for pressure
c
      xmeanp = 0.0
      grad_ug = ug(nnz) - ug((nnz-1))
      do iy=iys,iye
      do ix=1,nnx
         wfbc(ix,iy,1) = wf(ix,iy,izm1)-
     +                  (ux(ix,iy,iz) + vy(ix,iy,iz))*dzw(iz)
         pbc(ix,iy,1) = .5*(wf(ix,iy,izm1) + wfbc(ix,iy,1))
         ebc(ix,iy,1) = 0.0
         ubc(ix,iy,1) = u(ix,iy,iz) + grad_ug
         vbc(ix,iy,1) = v(ix,iy,iz)
         pbc2(ix,iy,1)=0.5*(u(ix,iy,iz)**2 + v(ix,iy,iz)**2) +
     +              0.25*(wf(ix,iy,izm1)**2 + wfbc(ix,iy,1)**2)
         xmeanp = xmeanp + pbc2(ix,iy,1)
c
c -------- FIX upper radiation boundary conditions
c
         ufbc(ix,iy,1) = ubc(ix,iy,1)
         vfbc(ix,iy,1) = vbc(ix,iy,1)
c        wfbc(ix,iy,1) = wbc(ix,iy,1)
         wbc(ix,iy,1)  = w(ix,iy,iz)
      enddo
      enddo
      call mpi_sum_xy(xmeanp,myid,iss,ise,1)
      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            tbc(ix,iy,iscl,1) = t(ix,iy,iscl,iz) +
     +                          dtdzf(iscl)*dzu(nnzp1)
         enddo
         enddo
      enddo
      xmeanp = xmeanp*fnxy
      do iy=iys,iye
      do ix=1,nnx
         pbc2(ix,iy,1) = pbc2(ix,iy,1) - xmeanp
      enddo
      enddo
c
c ---------- end if block
c
      endif
c
      do iy=iys,iye
      do ix=1,nnx
         w(ix,iy,izp1)   = wbc(ix,iy,1)
         e(ix,iy,izp1)   = ebc(ix,iy,1)
c
         wf(ix,iy,iz)   = wfbc(ix,iy,1)
         uf(ix,iy,izp1) = ufbc(ix,iy,1)
         vf(ix,iy,izp1) = vfbc(ix,iy,1)
c
         u(ix,iy,izp1) = ubc(ix,iy,1)
         v(ix,iy,izp1) = vbc(ix,iy,1)
c ------------- note w and e nnz+1 values are not needed
c        w(ix,iy,izp1) = wbc(ix,iy,1)
c        e(ix,iy,izp1) = ebc(ix,iy,1)
         r3(ix,iy,izp1)= 0.0
         r5(ix,iy,izp1)= 0.0
c
c ---------- set derivatives at top of box (wx,wy not needed)
c            ux,uy,vx,vy are used in e production, but neglect
c            at top of box becuase of bc
c
         wx(ix,iy,izp1) = 0.0
         wy(ix,iy,izp1) = 0.0
         ux(ix,iy,izp1) = 0.0
         uy(ix,iy,izp1) = 0.0
         vx(ix,iy,izp1) = 0.0
         vy(ix,iy,izp1) = 0.0
      enddo
      enddo
      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            t(ix,iy,iscl,izp1) = tbc(ix,iy,iscl,1)
            t(ix,iy,iscl,izp2) = tbc(ix,iy,iscl,1)
         enddo
         enddo
      enddo
c
      return
      end
      subroutine upper_dns
c
c ---- set boundary condition on upper boundary iz=nnz
c ---- for DNS conditions (no-slip, molecular diffusion)
c ---- assumes the grid transformation has completely decayed
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      iz   = nnz
      izm1 = iz - 1
      izm2 = iz - 2
      izp1 = iz + 1
      izp2 = iz + 2


c
c --------- boundary conditions are Sullivan et al. (2000) Couette flow
c
        do iy=iys,iye
        do ix=1,nnx
           wbc(ix,iy,1) = 0.0
           ebc(ix,iy,1) = 0.0
           ubc(ix,iy,1) = 2.0*Uo-u(ix,iy,iz)
           vbc(ix,iy,1) = -v(ix,iy,iz)
           pbc(ix,iy,1) = 0.0
           pbc2(ix,iy,1)= 0.0
           ufbc(ix,iy,1) = Uo
           vfbc(ix,iy,1) = 0.0
c
c --------- upper boundary = grid speed
c
           wfbc(ix,iy,1) = gm(ix,iy,iz,5)
        enddo
        enddo

!Set the scalar boundary condition based on isfc:
! isfc = 0, specified surface heat flux (through qstar)
! isfc = 1, specified surface temperature (through tsfcc)

        do iscl=1,nscl
           do iy=iys,iye
           do ix=1,nnx
              if (isfc(iscl)==1) then
              tbc(ix,iy,iscl,1) = 2.0*Ttop(iscl)-t(ix,iy,iscl,iz)
              end if
              if (isfc(iscl)==0) then
              tbc(ix,iy,iscl,1) = t(ix,iy,iscl,iz) +
     +           dzu(izp1)*wtsfc(iscl)/vis_s(ix,iy,izp1)
              end if
           enddo
           enddo
        enddo


      do iy=iys,iye
      do ix=1,nnx
         w(ix,iy,izp1)   = wbc(ix,iy,1)
         e(ix,iy,izp1)   = ebc(ix,iy,1)
c
         wf(ix,iy,iz)   = wfbc(ix,iy,1)
         uf(ix,iy,izp1) = ufbc(ix,iy,1)
         vf(ix,iy,izp1) = vfbc(ix,iy,1)
c
         u(ix,iy,izp1) = ubc(ix,iy,1)
         v(ix,iy,izp1) = vbc(ix,iy,1)
c ------------- note w and e nnz+1 values are not needed
c        w(ix,iy,izp1) = wbc(ix,iy,1)
c        e(ix,iy,izp1) = ebc(ix,iy,1)
         r3(ix,iy,izp1)= 0.0
         r5(ix,iy,izp1)= 0.0
c
c ---------- set derivatives at top of box (wx,wy not needed)
c            ux,uy,vx,vy are used in e production, but neglect
c            at top of box becuase of bc
c
         wx(ix,iy,izp1) = 0.0
         wy(ix,iy,izp1) = 0.0
         ux(ix,iy,izp1) = 0.0
         uy(ix,iy,izp1) = 0.0
         vx(ix,iy,izp1) = 0.0
         vy(ix,iy,izp1) = 0.0
      enddo
      enddo
      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            t(ix,iy,iscl,izp1) = tbc(ix,iy,iscl,1)
            t(ix,iy,iscl,izp2) = tbc(ix,iy,iscl,1)
         enddo
         enddo
      enddo
c
      return
      end
      subroutine psolver(istage)
c
c -------- solver poisson equation for the pressure
c
c          iflat == 0 ; no iteration flat bottom
c          iflat /= 0 ; iterate using diagonal preconditioner
c                       set max number of iterations ip_max
c                       and convergence tolerance p_tol
c          solve for fluctuating pressure then mean value.
c          array p(.,.,.) is p' array pm(.) is the mean value
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      call get_avg_j
c
      do k=izs-1,ize+1
         do j=iys,iye
         do i=1,nnx
            p(i,j,k) = 0.0
         enddo
         enddo
      enddo
c
c -------- special sequence for flat bottom
c
      if(iflat .eq. 0) then
         ip_iter = 0
         call build_rhs(ip_iter)
         call div_flux_vel
         call pressure
         call build_rhs(ip_iter)
      else
c
c ------- iteration sequence
c
         ip_iter = 0
    1    continue
            ip_iter = ip_iter + 1
            call build_rhs(ip_iter)
            call div_flux_vel
            call precon_p(ip_iter)
            call pressure
            call update_p(istage,ip_iter,ip_con)
         if(ip_iter .lt. ip_max .and. ip_con .eq. 0) go to 1
c
c -------- use last pressure iteration in flux velocity
c
         call build_rhs(ip_iter)
      endif
c
      return
      end
      subroutine build_rhs(ip_iter)
c
c ---------- build the rhs of the poisson equation
c            and update the flux velocities
c            use new time metrics. use wf equation to
c            get dp/dz at k = 0
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real fnt1(nnx,iys:iye,izs:ize+1), fnt2(nnx,iys:iye)
      real wf_sum(nnz)
      include 'mpif.h'
      integer istatus(mpi_status_size)
      real fs(nnx,iys:iye), fr(nnx,iys:iye)
c
c ------------ Send p to processor above the current myid.
c              to get p(.,.,izs-1)
c
      nb = myid - ncpu_s
      nt = myid + ncpu_s
      if(iss .eq. 0) then
         nb = mpi_proc_null
      endif
      if(ise .eq. numprocs-1) then
         nt = mpi_proc_null
      endif
      nsend = nnx*(iye + 1 - iys)
      nrecv = nsend
      do j=iys,iye
      do i=1,nnx
         fs(i,j) = p(i,j,ize)
      enddo
      enddo
c
      call mpi_sendrecv(
     +     fs(1,iys),nsend,mpi_real8,nt,2,
     +     fr(1,iys),nrecv,mpi_real8,nb,2,
     +     mpi_comm_world,istatus,ierr)
      if(iss .ne. 0) then
         do j=iys,iye
         do i=1,nnx
            p(i,j,izs-1) = fr(i,j)
         enddo
         enddo
      endif
c
      if(iss .eq. 0) then
c
c -------- flux velocity is set = grid speed at t = t_stage2
c          updated boundary speed at the next time step
c          picked to match future grid speed
c
         do j=iys,iye
         do i=1,nnx
            wf(i,j,izs-1) = wfbc(i,j,3)
         enddo
         enddo
      endif
c
c --------- get bottom ghost point for fluctuating pressure
c
      if(izs == 1) then
c
      if(lpbc == 2) then
        k    = izs
        kp1  = k+1
        km1  = k-1
        gami = 1.0/dtgama
        if(ip_iter <= 1) then
          do j=iys,iye
          do i=1,nnx
             fnt1(i,j,k) = p(i,j,k)
             fnt2(i,j)   = fnt1(i,j,k)
          enddo
          enddo
        else
          do j=iys,iye
          do i=1,nnx
             fnt1(i,j,k) = 0.5*(p(i,j,k) + p_old(i,j,km1))
             fnt2(i,j)   = fnt1(i,j,k)
          enddo
          enddo
        endif
        call xderivp(fnt2(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
        call yd_mpi(fnt1(1,iys,izs),trigx(1,2),yk(1),
     +            nnx,nny,ixs,ixe,ix_s,ix_e,
     +            iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
        do j=iys,iye
        do i=1,nnx
            fac   = gm2(i,j,k,1)*(gm2(i,j,k,2)**2 +
     +              gm2(i,j,k,3)**2 + gm2(i,j,k,4)**2)
            wbndy = (wfi(i,j,km1) - wf(i,j,km1))*gami
            facx  = gm2(i,j,k,1)*gm2(i,j,k,2)*fnt2(i,j)
            facy  = gm2(i,j,k,1)*gm2(i,j,k,3)*fnt1(i,j,k)
            p(i,j,km1) = p(i,j,k) - (wbndy - facx - facy)*dzu(k)/fac
        enddo
        enddo
c
      elseif(lpbc == 1) then
c
        w1 = -2.0*dzw(izs)/(dzw(izs) + dzw(izs+1))
        w2 = 1.0 - w1
        do j=iys,iye
        do i=1,nnx
           p(i,j,izs-1) = p(i,j,izs)*w2 + w1*p(i,j,izs+1)
        enddo
        enddo
c
      else
c
        do j=iys,iye
        do i=1,nnx
           p(i,j,izs-1) = p(i,j,izs)
        enddo
        enddo
      endif ! end of check for lpbc
c
      endif ! end of check for izs
c
      if(ize .eq. nnz) then
        do j=iys,iye
        do i=1,nnx
           p(i,j,ize+1) = p(i,j,ize)
        enddo
        enddo
      endif
c
c     if(izs .eq. 1)  then
c       k = 1
c       km1 = k - 1
c       kp1 = k + 1
c       j   = iys
c       write(nprt,8200) (i,p_old(i,j,km1), p(i,j,km1),
c    +                    p_old(i,j,k), p(i,j,k),
c    +                    i = 1,5)
c8200   format(' build rhs',/,
c    +         ' i ',5x,' p old 0 ',5x,' p 0',5x,
c    +         ' p old 1 ',5x,' p 1',/,
c    +         (i5,4e15.6))
c     endif
c
c --------- save the last pressure iteration
c
      do k=izs-1,ize+1
         do j=iys,iye
         do i=1,nnx
            p_old(i,j,k) = p(i,j,k)
         enddo
         enddo
      enddo
c
      do k=izs,ize+1
         do j=iys,iye
         do i=1,nnx
            fnt1(i,j,k) = p_old(i,j,k)
         enddo
         enddo
         call xderivp(fnt1(1,iys,k),trigx(1,1),xk(1),nnx,iys,iye)
      enddo
c
      do k=izs,ize
         kp1   = k + 1
         km1   = k - 1
         do j=iys,iye
         do i=1,nnx
            pu_zeta   = weit(k)*p_old(i,j,kp1) + weitm(k)*p_old(i,j,k)
            pl_zeta   = weit(km1)*p_old(i,j,k) +
     +                  weitm(km1)*p_old(i,j,km1)
            dpdz_c    = (pu_zeta - pl_zeta)*dzw_i(k)
            uf(i,j,k) = ufi(i,j,k) - dtgama*dpdz_c*
     +                  gm2(i,j,k,2)*gm2(i,j,k,1)
            vf(i,j,k) = vfi(i,j,k) - dtgama*dpdz_c*
     +                  gm2(i,j,k,3)*gm2(i,j,k,1)
            dpdz      = (p_old(i,j,kp1) - p_old(i,j,k))*dzu_i(kp1)
c
            facu      = gm2(i,j,kp1,1)*(gm2(i,j,kp1,2)**2 +
     +                  gm2(i,j,kp1,3)**2 + gm2(i,j,kp1,4)**2)
            facl      = gm2(i,j,k,1)*(gm2(i,j,k,2)**2 +
     +                  gm2(i,j,k,3)**2 + gm2(i,j,k,4)**2)
            wf(i,j,k) = wfi(i,j,k) - dtgama*(weit(k)*facu +
     +                     weitm(k)*facl)*dpdz
         enddo
         enddo
         do j=iys,iye
         do i=1,nnx
            uf(i,j,k) = uf(i,j,k) - dtgama*fnt1(i,j,k)*gm2(i,j,k,1)
            facu      = gm2(i,j,kp1,1)*gm2(i,j,kp1,2)*fnt1(i,j,kp1)
            facl      = gm2(i,j,k,1)*gm2(i,j,k,2)*fnt1(i,j,k)
            wf(i,j,k) = wf(i,j,k) - dtgama*
     +                 (weit(k)*facu + weitm(k)*facl)
         enddo
         enddo
      enddo
c
c -------- get the y contributions
c
      do k=izs,ize+1
         do j=iys,iye
         do i=1,nnx
            fnt1(i,j,k) = p_old(i,j,k)
         enddo
         enddo
      enddo
      call yd_mpi(fnt1(1,iys,izs),trigx(1,2),yk(1),
     +            nnx,nny,ixs,ixe,ix_s,ix_e,
     +            iys,iye,iy_s,iy_e,izs,ize+1,myid,ncpu_s,numprocs)
c
      do k=izs,ize
         kp1 = k + 1
         km1 = k - 1
         do j=iys,iye
         do i=1,nnx
            vf(i,j,k) = vf(i,j,k) - dtgama*fnt1(i,j,k)*gm2(i,j,k,1)
            facu      = gm2(i,j,kp1,1)*gm2(i,j,kp1,3)*fnt1(i,j,kp1)
            facl      = gm2(i,j,k,1)*gm2(i,j,k,3)*fnt1(i,j,k)
            wf(i,j,k) = wf(i,j,k) - dtgama*
     +                  (weit(k)*facu + weitm(k)*facl)
         enddo
         enddo
      enddo
c
      if(ise .eq. numprocs-1) then
c
c -------- upper boundary speed at t = t_stage2
c
         do j=iys,iye
         do i=1,nnx
            wf(i,j,ize) = wfbc(i,j,4)
         enddo
         enddo
      endif
c
c ------------ make sure continuity is satisfied
c
      do k=1,nnz
         wf_sum(k) = 0.0
      enddo
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            wf_sum(k) = wf_sum(k) + wf(i,j,k)
         enddo
         enddo
         wf_sum(k) = wf_sum(k)*fnxy
      enddo
c
c ---------- gather partial sums
c
      call mpi_sum_z(wf_sum(1),i_root,myid,nnz,1)
c
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            wf(i,j,k) = wf(i,j,k) - wf_sum(k)
         enddo
         enddo
      enddo
c
      return
      end
      subroutine div_flux_vel
c
c ----------- get the divergence of the flux velocity
c             and impose the boundary conditions on wf
c             new time metrics and weight for pressure
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      real fnt1(nnx,iys:iye,izs:ize)
      real fs(nnx,iys:iye), fr(nnx,iys:iye)
      real sum_div(1:nnz)
c
      gami = 1.0/dtgama
c
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            fnt1(i,j,k) = uf(i,j,k)
         enddo
         enddo
         call xderivp(fnt1(1,iys,k),trigx(1,1),xk(1),nnx,iys,iye)
         do j=iys,iye
         do i=1,nnx
            p(i,j,k) = fnt1(i,j,k)
         enddo
         enddo
      enddo
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            fnt1(i,j,k) = vf(i,j,k)
         enddo
         enddo
      enddo
      call yd_mpi(fnt1(1,iys,izs),trigx(1,2),yk(1),
     +            nnx,nny,ixs,ixe,ix_s,ix_e,
     +            iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
c
c ------------ Send updated wf (from psolver)
c              to processor above the current myid.
c
      nb = myid - ncpu_s
      nt = myid + ncpu_s
      if(iss .eq. 0) then
         nb = mpi_proc_null
      endif
      if(ise .eq. numprocs-1) then
         nt = mpi_proc_null
      endif
      nsend = nnx*(iye + 1 - iys)
      nrecv = nsend
      do j=iys,iye
      do i=1,nnx
         fs(i,j) = wf(i,j,ize)
      enddo
      enddo
c
      call mpi_sendrecv(
     +     fs(1,iys),nsend,mpi_real8,nt,2,
     +     fr(1,iys),nrecv,mpi_real8,nb,2,
     +     mpi_comm_world,istatus,ierr)
      if(iss .ne. 0) then
         do j=iys,iye
         do i=1,nnx
            wf(i,j,izs-1) = fr(i,j)
         enddo
         enddo
      endif
c
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            div_wf   = (wf(i,j,k) - wf(i,j,k-1))*dzw_i(k)
            p(i,j,k) = (p(i,j,k) + fnt1(i,j,k) + div_wf)*gami
         enddo
         enddo
      enddo
c
c debug -------
c
      do k=1,nnz
         sum_div(k) = 0.0
      enddo
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            sum_div(k) = sum_div(k) + (p(i,j,k)/gami)**2
         enddo
         enddo
         sum_div(k) = sum_div(k)*fnxy
      enddo
c
c ---------- gather partial sums for rj_dif
c
      call mpi_sum_z(sum_div,i_root,myid,nnz,1)
c
c     write(nprt,9010) (k, sum_div(k),k=izs,ize)
c9010 format(' 9010 debug div flux ',/,
c    +       ' k ',5x,' div ',/,(i5,e20.10))
c
      return
      end
      subroutine get_avg_j
c
c ------- compute the average of (1/J,J) for use in preconditioner
c         metrics are at the new time
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
c --------- 1/J values stored at cell centers
c
      do k=0,nnz+1
         avg_j(k) = 0.0
      enddo
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            avg_j(k) = avg_j(k) + gm2(i,j,k,1)
         enddo
         enddo
         avg_j(k) = avg_j(k)*fnxy
      enddo
c
c ---------- gather partial sums for avg_j
c
      call mpi_sum_z(avg_j(1),i_root,myid,nnz,1)
c
c ---------- J values stored at cell faces
c
      avg_j(0)     = avg_j(1)
      avg_j(nnz+1) = avg_j(nnz)
      do k=1,nnz
         fac        = avg_j(k+1)*weit(k) + avg_j(k)*weitm(k)
         avg_j_i(k) = 1.0/fac
      enddo
      avg_j_i(0)     = 1.0/avg_j(1)
      avg_j_i(nnz+1) = 1.0/avg_j(nnz)
c
c     write(nprt,1001) (k,avg_j(k),avg_j_i(k),k=0,nnz+1)
c1001 format(' 1000 get_avg_j ',/,
c    +       ' k ',5x,' avg j',5x,' avg j_i ',/,
c    +       (i5,2e15.6))
c
      return
      end
      subroutine precon_p(ip_iter)
c
c -------- evaluate the preconditioner at the last iteration
c          for p, p_old has last iteration
c          new weight for pressure
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      real fnt1(nnx,iys:iye,izs:ize)
c
      gami   = 1.0/dtgama
      relax  = accel(ip_iter)
      relax1 = 1.0 - relax
      relax1 = 1.0
c     write(nprt, 1017) ip_iter, relax
c1017 format(' 1017: precon ip_iter = ',i5,
c    +       ' relax = ',e15.6)
c
c ------- get x contribution
c
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            fnt1(i,j,k) = p_old(i,j,k)
         enddo
         enddo
         call xdd(fnt1(1,iys,k),trigx(1,1),xkk(1),nnx,iys,iye)
         do j=iys,iye
         do i=1,nnx
            p(i,j,k) = fnt1(i,j,k)*avg_j(k)*relax1 + p(i,j,k)*relax
         enddo
         enddo
      enddo
c
c ------- get y contribution
c
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            fnt1(i,j,k) = p_old(i,j,k)
         enddo
         enddo
      enddo
      call ydd_mpi(fnt1(1,iys,izs),trigx(1,2),ykk(1),
     +            nnx,nny,ixs,ixe,ix_s,ix_e,
     +            iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            p(i,j,k) = p(i,j,k) + fnt1(i,j,k)*avg_j(k)*relax1
         enddo
         enddo
      enddo
c
c ------- get z contribution
c         always use dp/dz = 0 for pre-conditioner
c         at izs=1
c
c
      do k=max(izs,2),ize
         kp1    = k + 1
         km1    = k - 1
         facu_i = avg_j_i(k)
         facl_i = avg_j_i(km1)
         do j=iys,iye
         do i=1,nnx
            dpdz_u = facu_i*(p_old(i,j,kp1) - p_old(i,j,k))*dzu_i(kp1)
            dpdz_l = facl_i*(p_old(i,j,k) - p_old(i,j,km1))*dzu_i(k)
            p(i,j,k) = p(i,j,k) + (dpdz_u - dpdz_l)*dzw_i(k)*relax1
         enddo
         enddo
      enddo
      if(izs == 1) then
         k      = izs
         kp1    = izs + 1
         facu_i = avg_j_i(k)
         do j=iys,iye
         do i=1,nnx
            dpdz_u   = facu_i*(p_old(i,j,kp1) - p_old(i,j,k))*dzu_i(kp1)
            p(i,j,k) = p(i,j,k) + dpdz_u*dzw_i(k)*relax1
         enddo
         enddo
      endif
c
c     if(izs .eq. 1)  then
c       k = 1
c       km1 = k - 1
c       kp1 = k + 1
c       j   = iys
c       write(nprt,8020) (i, p(i,j,km1), p(i,j,k), p(i,j,k+1),
c    +                    i = 1,5)
c8020   format(' 8020 precon',/,
c    +         ' i ',5x,' p source  0 ',5x,' p 1',5x,
c    +         ' p 2',/,
c    +         (i5,3e15.6))
c     endif
c
      return
      end
      subroutine xdd(ax,trigx,xkk,nnx,iys,iye)
c
c -------- get multiple x second derivatives using fftpack routines
c          use fftpack storage a0, (a1,b1), (a2,b2),...,an
c          assumes that first wavenumber xk(1) = 0.0
c
c          assumes that wavenumbers are normalized by number of points
c
      real xkk(nnx), trigx(2*nnx+15), ax(nnx,iys:iye)
c
c     fn = 1.0/float(nnx)
      do iy=iys,iye
         call rfftf(nnx,ax(1,iy),trigx)
         ii = 1
         ax(1,iy) = 0.0
         ax(nnx,iy) = 0.0
         do ix=2,nnx-1,2
            ii          = ii + 1
            ax(ix,iy)   = -xkk(ii)*ax(ix,iy)
            ax(ix+1,iy) = -xkk(ii)*ax(ix+1,iy)
         enddo
         call rfftb(nnx,ax(1,iy),trigx)
      enddo
c
      return
      end
      subroutine ydd_mpi(ay,trigy,ykk,
     +           nx,ny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
c
c -------- get multiple y second derivatives using fftpack routines and mpi
c          use fftpack storage a_0, (a1,b1), (a2,b2), ...,
c          assumes that first wavenumber yk(1) = 0.0
c          wavenumbers are normalized by number of points, ny
c
      real ykk(ny), trigy(2*ny+15), ay(nx,iys:iye,iz1:iz2)
      real ayt(ny,ixs:ixe,iz1:iz2)
c
      integer ix_s(0:np-1), ix_e(0:np-1),
     +        iy_s(0:np-1), iy_e(0:np-1)
c
      call xtoy_trans(ay,ayt,nx,ny,ixs,ixe,ix_s,ix_e,
     +         iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
c
c     fn = 1.0/float(nny)
      do iz=iz1,iz2
         do ix=ixs,ixe
            call rfftf(ny,ayt(1,ix,iz),trigy)
            ii = 1
            ayt(1,ix,iz)  = 0.0
            ayt(ny,ix,iz) = 0.0
            do iy=2,ny-1,2
               ii              = ii + 1
               ayt(iy,ix,iz)   = -ykk(ii)*ayt(iy,ix,iz)
               ayt(iy+1,ix,iz) = -ykk(ii)*ayt(iy+1,ix,iz)
            enddo
            call rfftb(ny,ayt(1,ix,iz),trigy)
         enddo
      enddo
      call ytox_trans(ayt,ay,nx,ny,ixs,ixe,ix_s,ix_e,
     +         iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
c
      return
      end
      subroutine update_p(istage,ip_iter,ip_con)
c
c -------- check convergence and update p guess
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
      integer istatus(mpi_status_size)
      data small_p /1.0e-10/
      save small_p
c
      pmax  = 0.0
      pdiff = 0.0
      do k=izs,ize
         p_temp  = abs(p(1,iys,k))
         pd_temp = abs(p(1,iys,k) - p_old(1,iys,k))
         do j=iys,iye
         do i=1,nnx
            p_temp  = amax1(abs(p(i,j,k)),p_temp)
            pd_temp = amax1(abs(p(i,j,k) - p_old(i,j,k)),pd_temp)
         enddo
         enddo
         pmax  = amax1(p_temp,pmax)
         pdiff = amax1(pd_temp,pdiff)
      enddo
      pmax      = amax1(pmax,small_p)
      ptest_loc = pdiff/pmax
c
      call mpi_allreduce(ptest_loc,ptest,1,mpi_real8,
     +     mpi_max,mpi_comm_world,ierror)
c
c     write(nprt,6000) ip_iter, pmax, ptest, ptol
c6000 format(' ip_iter = ',i5,' pmax = ',e15.6,
c    +       ' ptest = ',e15.6,' ptol = ',e15.6)
c
      ip_con = 0
      if(ptest .lt. ptol) ip_con = 1
      if(istage .eq. 1) then
         p_iter = float(ip_iter)
         p_err  = ptest
      endif
c
      return
      end
      subroutine comp2
c
c ------- add p gradients to rhs of (u,v,w). Use already defined p
c         at izs-1 and ize+1 to get (u,v,w). see build_rhs
c
c         multiplied by gm2 to get (u,v,w,e,theta) at the new time
c         pressure gradient is at new time and is weighted by dtgama
c
c         included mean pressure in vertical gradients in sr. rhs_uvw
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      real fnt1(nnx,iys:iye,izs:ize)
c
c ---------  ghost point pressure comes from build_rhs
c
c     if(izs .eq. 1) then
c       j = iys
c       k = 1
c       write(nprt,4001) (i,p(i,j,k-1), p(i,j,k), p(i,j,k+1),
c    +                    p(i,j,k+2),i=1,nnx)
c4001   format(' 4001 comp2 ',/,
c    +   ' i ',5x,' p km1 ',10x,' p k ',10x,' p kp1',10x,' p kp2',/,
c    +   (i5,4e15.6))
c     endif
c
      do k=izs,ize
         km1 = k - 1
         kp1 = k + 1
         do j=iys,iye
         do i=1,nnx
            fnt1(i,j,k) = p(i,j,k)*gm2(i,j,k,1)
         end do
         end do
         call xderivp(fnt1(1,iys,k),trigx(1,1),xk(1),nnx,iys,iye)
c
         do j=iys,iye
         do i=1,nnx
            p_kp1     = p(i,j,kp1)
            p_k       = p(i,j,k)
            p_km1     = p(i,j,km1)
c
c ------- u equation
c
            pu_zeta   = weit(k)*p_kp1*gm2(i,j,kp1,2)*
     +                  gm2(i,j,kp1,1) +
     +                  weitm(k)*p_k*gm2(i,j,k,2)*gm2(i,j,k,1)
            pl_zeta   = weit(km1)*p_k*gm2(i,j,k,2)*gm2(i,j,k,1) +
     +                  weitm(km1)*p_km1*gm2(i,j,km1,2)*
     +                  gm2(i,j,km1,1)
c
            dpdz_c    = (pu_zeta - pl_zeta)*dzw_i(k)
c
c --------  add pressure gradient to rhs
c
            r1(i,j,k) = r1(i,j,k) - (fnt1(i,j,k) + dpdz_c)
            u(i,j,k)  = (u(i,j,k) - dtgama*(fnt1(i,j,k) + dpdz_c))*
     +                  gm2(i,j,k,4)
c
c ------- w equation (use zeta_z/J = 1)
c
            pu_zeta   = weit(k)*p_kp1 + weitm(k)*p_k
            pl_zeta   = weit(km1)*p_k + weitm(km1)*p_km1
            dpdz_c    = (pu_zeta - pl_zeta)*dzw_i(k)
c
c -------- add pressure gradient to rhs
c
            r3(i,j,k) = r3(i,j,k) - dpdz_c
            w(i,j,k)  = (w(i,j,k) - dtgama*dpdz_c)*gm2(i,j,k,4)
c
c ------- e equation correct for Jacobian, result is e (not e/J)
c
            e(i,j,k) = e(i,j,k)*gm2(i,j,k,4)
         enddo
         enddo
      enddo
c
c ------ v equation
c
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            fnt1(i,j,k) = p(i,j,k)*gm2(i,j,k,1)
         enddo
         enddo
      enddo
      call yd_mpi(fnt1(1,iys,izs),trigx(1,2),yk(1),
     +            nnx,nny,ixs,ixe,ix_s,ix_e,
     +            iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
      do k=izs,ize
         kp1   = k + 1
         km1   = k - 1
         do j=iys,iye
         do i=1,nnx
            p_kp1     = p(i,j,kp1)
            p_k       = p(i,j,k)
            p_km1     = p(i,j,km1)
            pu_zeta   = weit(k)*p_kp1*gm2(i,j,kp1,3)*
     +                  gm2(i,j,kp1,1) +
     +                  weitm(k)*p_k*gm2(i,j,k,3)*
     +                  gm2(i,j,k,1)
            pl_zeta   = weit(km1)*p_k*gm2(i,j,k,3)*
     +                  gm2(i,j,k,1) +
     +                  weitm(km1)*p_km1*gm2(i,j,km1,3)*
     +                  gm2(i,j,km1,1)
            dpdz_c    = (pu_zeta - pl_zeta)*dzw_i(k)
c
c --------  add pressure gradient to rhs now
c
            r2(i,j,k) = r2(i,j,k) - (fnt1(i,j,k) + dpdz_c)
            v(i,j,k)  = (v(i,j,k) - dtgama*(fnt1(i,j,k) + dpdz_c))*
     +                  gm2(i,j,k,4)
         enddo
         enddo
      enddo
c
c ------------ scalar equations correct for J, result is theta (not theta/J)
c
      do k=izs,ize
         do l=1,nscl
         do j=iys,iye
         do i=1,nnx
            t(i,j,l,k) = t(i,j,l,k)*gm2(i,j,k,4)
         enddo
         enddo
         enddo
      enddo
c
      return
      end
      subroutine pressure
c
c -------- solve for fluctuating pressure using a matrix transpose
c          across mpi tasks and tridiagonal solver.
c          The transposed array
c          is dimensioned (0:nnz+1). Values
c          (0 & nnz+1) are not needed but are useful in the
c          matrix transpose when we return (see send_ztox).
c          On exit p is defined at all [izs-1:ize+1].
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real pfft(nny,jxs:jxe,izs-1:ize+1)
      real pt(0:nnz+1,jxs:jxe,jys:jye)
      real ptopfft(nny,jxs:jxe,1:2)
      real psum(1:nnz)
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
c ------------ Fourier analyze the right hand side
c              at all iz = izs,ize. results are in pfft
c
c
      call fft2d_mpi(p(1,iys,izs),pfft(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)
c
c ------------ Fourier analyze the radiation bc arrays
c
      if(ibcu .eq. 1) then
        call fft2d_mpi(ptop(1,iys,1),ptopfft(1,jxs,1),
     +           trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           1,2,myid,ncpu_s,numprocs,-2)
      endif
c
c ---------- transpose first and last index of array
c            the order of pfft is (y,x,z)
c
      call xtoz_trans(pfft,pt,nny,nnz,jys,jye,jy_s,jy_e,
     +                jxs,jxe,izs,ize,iz_s,iz_e,myid,ncpu_s,
     +                numprocs)
      call solve_trid(pt, ptopfft)
c
c ------------- transpose back
c
      call ztox_trans(pt,pfft,nny,nnz,jys,jye,jy_s,jy_e,
     +                jxs,jxe,izs,ize,iz_s,iz_e,myid,ncpu_s,
     +                numprocs)
c
      iz_ee = ize+1
      if(ise .eq. numprocs-1) then
         iz_ee = ize
      endif
c
c --------- inverse fft at all iz=izs,iz_ee to get p
c           see z indices
c
      call fft2d_mpi(p(1,iys,izs),pfft(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,iz_ee,myid,ncpu_s,numprocs,2)
c
c -------- partial sums for pressure
c
      do iz=1,nnz
         psum(iz) = 0.0
      enddo
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            psum(iz) = psum(iz) + p(ix,iy,iz)
         enddo
         enddo
         psum(iz) = psum(iz)*fnxy
      enddo
      call mpi_sum_z(psum,i_root,myid,nnz,1)
c
      do iz=izs,iz_ee
c        psum(iz) = -psum(iz) + engz(iz) + c23*engsbz(iz)
         do iy=iys,iye
         do ix=1,nnx
            p(ix,iy,iz) = p(ix,iy,iz) - psum(iz)
         enddo
         enddo
      enddo
c
      return
      end
      subroutine solve_trid(pt, ptop)
c
c --------- tridiagonal solver. odd order for ptop, ptop2
c           because of 2d-fft
c
      use pars
      use con_data
      use con_stats
c
      real ptop(nny,jxs:jxe,1:2)
      real pt(0:nnz+1,jxs:jxe,jys:jye)
      real aa(nnz,jxs:jxe),bb(nnz,jxs:jxe),
     +     dd(nnz,jxs:jxe),rh(nnz,jxs:jxe)
      real fac_u(nnz), fac_l(nnz), fac_a(nnz)
c
      do iz=1,nnz
         fac_u(iz) = avg_j_i(iz)/(dzw(iz)*dzu(iz+1))
         fac_l(iz) = avg_j_i(iz-1)/(dzw(iz)*dzu(iz))
         fac_a(iz) = fac_l(iz) + fac_u(iz)
      enddo
c
      do kp=jys,jye
         do lp=jxs,jxe
         do iz=2,nnz-1
            bb(iz,lp)  = fac_l(iz)
            aa(iz,lp)  = fac_u(iz)
            dd(iz,lp)  = -xks(lp,kp)*avg_j(iz) - fac_a(iz)
            rh(iz,lp)  = pt(iz,lp,kp)
         enddo
         enddo
c
c --------------- lower boundary, fill exterior pressure (not used)
c
         do lp=jxs,jxe
            bb(1,lp)  = 1.0
            aa(1,lp)  = fac_u(1)
            dd(1,lp)  = -xks(lp,kp)*avg_j(1) - fac_u(1)
            rh(1,lp)  = pt(1,lp,kp)
            pt(0,lp,kp) = 0.0
         enddo
c
c --------------- upper boundary, fill exterior pressure (not used)
c
         if(ibcu .eq. 1) then
            do lp=jxs,jxe
              bb(nnz,lp) = 0.0
              aa(nnz,lp) = 0.0
              dd(nnz,lp) = 1.0
              rh(nnz,lp) = ptop(kp,lp,1)*wavexy(lp,kp) + ptop(kp,lp,2)
              pt(nnz+1,lp,kp) = 0.0
            enddo
         else
            do lp=jxs,jxe
               bb(nnz,lp) = fac_l(nnz)
               aa(nnz,lp) = 1.0
               dd(nnz,lp) = -xks(lp,kp)*avg_j(nnz) - fac_l(nnz)
               rh(nnz,lp) = pt(nnz,lp,kp)
               pt(nnz+1,lp,kp) = 0.0
            enddo
         endif
c
c ---------------- special situation for zeroth mode
c                  makes mean pressure = 0
c
         if(kp .eq. 1 .and. jxs .eq. 1) then
           do iz=1,nnz
              dd(iz,1) = 1.0
              rh(iz,1) = 0.0
              aa(iz,1) = 0.0
              bb(iz,1) = 0.0
              dd(iz,2) = 1.0
              rh(iz,2) = 0.0
              aa(iz,2) = 0.0
              bb(iz,2) = 0.0
           enddo
         endif
c
c --------------- solve system
c
         call tridv(bb,dd,aa,rh,nnz,jxs,jxe)
         do lp=jxs,jxe
         do iz=1,nnz
            pt(iz,lp,kp) = rh(iz,lp)
         enddo
         enddo
      enddo
c
      return
      end
      subroutine solve_mean_p(wf_sum,sum_fz)
c
c --------- tridiagonal solver for mean pressure from wf equation
c           set ghost values of pm also
c
      use pars
      use con_data
      use con_stats
c
      real wf_sum(0:nnz), sum_fz(nnz)
      real aa(nnz),bb(nnz),dd(nnz),rh(nnz)
c
      gami = 1.0/dtgama
c     tiny = 1.0e-20
      tiny = 0.0
c
      do k=2,nnz-1
         fac_l  = sum_fz(k-1)*dzu_i(k)
         fac_u  = sum_fz(k)*dzu_i(k+1)
         bb(k)  = fac_l
         aa(k)  = fac_u
         dd(k)  = -(fac_u + fac_l + tiny)
         rh(k)  = gami*(wf_sum(k) - wf_sum(k-1))
      enddo
c
c --------------- lower boundary, assume average wf at boundary = 0
c
      k       = 1
      fac_u   = sum_fz(k)*dzu_i(k+1)
      bb(k)   = 0.0
      aa(k)   = fac_u
      dd(k)   = -(fac_u + tiny)
      rh(k)   = gami*(wf_sum(k) - wf_sum(k-1))
c
c --------------- upper boundary, average wf at boundary = 0
c
      k       = nnz
      fac_l   = sum_fz(k-1)*dzu_i(k)
      bb(k)   = fac_l
      aa(k)   = 0.0
      dd(k)   = -(fac_l + tiny)
      rh(k)   = gami*(wf_sum(k) - wf_sum(k-1))
c
c --------------- solve system
c
      call trids(bb,dd,aa,rh,nnz)
      do k=1,nnz
         pm(k) = rh(k)
      enddo
c
      if(lpbc == 1) then
        k       = 1
        w1      = -2.0*dzw(k)/(dzw(k) + dzw(k+1))
        w2      = 1.0 - w1
        pm(k-1) = pm(k)*w2 + w1*pm(k+1)
      elseif(lpbc == 0) then
        k       = 1
        pm(k-1) = pm(k)
      endif
c
      pm(nnz+1) = pm(nnz)
c
      return
      end
      subroutine solve_mean_p_old(wf_sum,sum_fz)
c
c --------- tridiagonal solver for mean pressure from wf equation
c           set ghost values of pm also
c
      use pars
      use con_data
      use con_stats
c
      real wf_sum(nnz), sum_fz(nnz)
      real aa(nnz),bb(nnz),dd(nnz),rh(nnz)
c
      gami = 1.0/dtgama
c     tiny = 1.0e-20
      tiny = 0.0
c
      do k=2,nnz-1
         fac_l  = sum_fz(k-1)/(dzw(k)*dzu(k))
         fac_u  = sum_fz(k)/(dzw(k)*dzu(k+1))
         bb(k)  = fac_l
         aa(k)  = fac_u
         dd(k)  = -(fac_u + fac_l + tiny)
         rh(k)  = gami*(wf_sum(k) - wf_sum(k-1))*dzw_i(k)
      enddo
c
c --------------- lower boundary, assume average wf at boundary = 0
c
      k       = 1
      fac_u   = sum_fz(k)/(dzw(k)*dzu(k+1))
      bb(k)   = 0.0
      aa(k)   = fac_u
      dd(k)   = -(fac_u + tiny)
      rh(k)   = gami*wf_sum(k)*dzw_i(k)
c
c --------------- upper boundary, average wf at boundary = 0
c
      k       = nnz
      fac_l   = sum_fz(k-1)/(dzw(k)*dzu(k))
      bb(k)   = fac_l
      aa(k)   = 0.0
      dd(k)   = -(fac_l + tiny)
      rh(k)   = -gami*wf_sum(k-1)*dzw_i(k)
c
c --------------- solve system
c
      call trids(bb,dd,aa,rh,nnz)
      do k=1,nnz
         pm(k) = rh(k)
      enddo
c
      if(lpbc == 1) then
        k       = 1
        w1      = -2.0*dzw(k)/(dzw(k) + dzw(k+1))
        w2      = 1.0 - w1
        pm(k-1) = pm(k)*w2 + w1*pm(k+1)
      elseif(lpbc == 0) then
        k       = 1
        pm(k-1) = pm(k)
      endif
c
      pm(nnz+1) = pm(nnz)
c
      return
      end
      subroutine trids(b,d,a,r,n)
c
c --- tridiagonal matrix solver with single vector
c
c --- input:   n   size of a,b,d and r
c              b   below diagonal elements (b(1) not used)
c              d   diagonal elements
c              a   above diagonal elements (a(n) not used)
c              r   right hand side
c
c --- output:  r   solution vector
c
      real b(n), d(n), a(n), r(n)
c
      if(n .le. 1 ) then
        r(1) = r(1)/d(1)
         go to 999
      endif
      d(1) = 1.0/d(1)
      do i=2,n
         fac = b(i)*d(i-1)
         d(i) = 1.0/(d(i) - fac*a(i-1))
         r(i) = r(i) - fac*r(i-1)
      enddo
      r(n) = r(n)*d(n)
      do i=n-1,1,-1
         r(i) = d(i)*(r(i) - a(i)*r(i+1))
      enddo
  999 continue
c
      return
      end
      subroutine tridv(b,d,a,r,n,j1,j2)
c
c --- tridiagonal matrix solver with multiple vectors
c     (note j and i loops are reversed from cray version)
c
c --- input:   n   size of a,b,d and r
c              b   below diagonal elements (b(1) not used)
c              d   diagonal elements
c              a   above diagonal elements (a(n) not used)
c              r   right hand side
c              j1:j2  range of input vectors
c
c --- output:  r   solution vector
c
      real b(n,j1:j2), d(n,j1:j2), a(n,j1:j2), r(n,j1:j2)
c
      if(n .le. 1 ) then
         do j=j1,j2
            r(1,j) = r(1,j)/d(1,j)
         enddo
         go to 999
      endif
      do j=j1,j2
         d(1,j) = 1.0/d(1,j)
      enddo
      do j=j1,j2
      do i=2,n
         fac = b(i,j)*d(i-1,j)
         d(i,j) = 1.0/(d(i,j) - fac*a(i-1,j))
         r(i,j) = r(i,j) - fac*r(i-1,j)
      enddo
      enddo
      do j=j1,j2
         r(n,j) = r(n,j)*d(n,j)
      enddo
      do j=j1,j2
      do i=n-1,1,-1
         r(i,j) = d(i,j)*(r(i,j) - a(i,j)*r(i+1,j))
      enddo
      enddo
  999 continue
c
      return
      end
      subroutine get_derv
c
c ------- get ux,uy,vx,vy at all z for this node
c         using parallel fft. can be improved (?)
c         by using exchange to send derivatives
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      iz_ss = izs-1
      iz_ee = ize+1
      if(iss .eq. 0) then
         iz_ss = izs
      endif
      if(ise .eq. numprocs-1) then
         iz_ee = ize
      endif
c
c ------- make sure <w> = 0
c
      do iz=izs-1,ize+1
         w_sum = 0.0
         do iy=iys,iye
         do ix=1,nnx
            w_sum = w_sum + w(ix,iy,iz)
         enddo
         enddo
         w_sum = w_sum*fnxy
         call mpi_sum_xy(w_sum,myid,iss,ise,1)
         do iy=iys,iye
         do ix=1,nnx
            w(ix,iy,iz) = w(ix,iy,iz) - w_sum
         enddo
         enddo
      enddo
c
      do iz=izs-1,ize+1
         do iy=iys,iye
         do ix=1,nnx
            ux(ix,iy,iz) = u(ix,iy,iz)
            vx(ix,iy,iz) = v(ix,iy,iz)
            wx(ix,iy,iz) = w(ix,iy,iz)
            uy(ix,iy,iz) = u(ix,iy,iz)
            vy(ix,iy,iz) = v(ix,iy,iz)
            wy(ix,iy,iz) = w(ix,iy,iz)
         enddo
         enddo
c        endif
         call xderivp(ux(1,iys,iz),trigx(1,1),xk(1),
     +                 nnx,iys,iye)
         call xderivp(vx(1,iys,iz),trigx(1,1),xk(1),
     +                 nnx,iys,iye)
         call xderivp(wx(1,iys,iz),trigx(1,1),xk(1),
     +                 nnx,iys,iye)
      enddo
c
c ---------- get y derivatives for (u,v,w)
c
      call yd_mpi(uy(1,iys,izs-1),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
      call yd_mpi(vy(1,iys,izs-1),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
      call yd_mpi(wy(1,iys,izs-1),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
c
      return
      end
      subroutine chk_div
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real fnt1(nnx,iys:iye,izs:ize+1), fnt2(nnx,iys:iye)
      real div_c(nnz)
c
      do iz=1,nnz
         div_c(iz) = 0.0
      enddo
        do k=izs,ize
           do j=iys,iye
           do i=1,nnx
              fnt1(i,j,k) = vf(i,j,k)
           enddo
           enddo
        enddo
        call yd_mpi(fnt1(1,iys,izs),trigx(1,2),yk(1),
     +              nnx,nny,ixs,ixe,ix_s,ix_e,
     +              iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
        do k=izs,ize
           do j=iys,iye
           do i=1,nnx
              fnt2(i,j) = uf(i,j,k)
           enddo
           enddo
           call xderivp(fnt2(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
           km1 = k - 1
           do j=iys,iye
           do i=1,nnx
              div_c(k) = div_c(k) +
     +                  (fnt2(i,j) + fnt1(i,j,k) +
     +                  (wf(i,j,k) - wf(i,j,km1))*dzw_i(k))**2
           enddo
           enddo
           div_c(k) = div_c(k)*fnxy
        enddo
        call mpi_sum_z(div_c(1),i_root,myid,nnz,1)
c
      write(nprt,9010) (k, div_c(k),k=izs,ize)
 9010 format(' 9010 debug chk div',/,
     +       ' k ',5x,' div ',/,(i5,e20.10))
c
      return
      end
      subroutine get_means(istage)
c
c ------------ get means for all variables
c              for use in iso, surfvis, comp1, tke_stats.
c              simple average along coordinate lines
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real fnt1(nnx,iys:iye,izs:ize+1), fnt2(nnx,iys:iye)
c
      do iz=0,nnz+1
         u_mn(iz)   = 0.0
         v_mn(iz)   = 0.0
         w_mn(iz)   = 0.0
         engz(iz)   = 0.0
         engsbz(iz) = 0.0
         divz(iz)   = 0.0
         pxym(iz)   = 0.0
      enddo
      do iscl=1,nscl
         do iz=0,nnz+1
            t_mn(iz,iscl) = 0.0
         enddo
      enddo
      iz_ee = ize
      if(ize .eq. nnz) iz_ee = nnzp1
      do iz=izs,iz_ee
         do iy=iys,iye
         do ix=1,nnx
            u_mn(iz) = u_mn(iz) + u(ix,iy,iz)
            v_mn(iz) = v_mn(iz) + v(ix,iy,iz)
            w_mn(iz) = w_mn(iz) + w(ix,iy,iz)
         enddo
         enddo
         u_mn(iz) = u_mn(iz)*fnxy
         v_mn(iz) = v_mn(iz)*fnxy
         w_mn(iz) = w_mn(iz)*fnxy
         do iscl=1,nscl
            t_mn(iz,iscl) = 0.0
            do iy=iys,iye
            do ix=1,nnx
               t_mn(iz,iscl) = t_mn(iz,iscl) + t(ix,iy,iscl,iz)
            enddo
            enddo
            t_mn(iz,iscl) = t_mn(iz,iscl)*fnxy
         enddo
      enddo
      call mpi_sum_z(u_mn(1),i_root,myid,nnzp1,1)
      call mpi_sum_z(v_mn(1),i_root,myid,nnzp1,1)
      call mpi_sum_z(w_mn(1),i_root,myid,nnzp1,1)
      do iscl=1,nscl
         call mpi_sum_z(t_mn(1,iscl),i_root,myid,nnzp1,1)
      enddo
c
c -------- set e to minimum value
c
      do iz=izs-1,ize+1
         do iy=iys,iye
         do ix=1,nnx
            e(ix,iy,iz) = amax1(e(ix,iy,iz ),sml_eg)
         enddo
         enddo
      enddo
c
c ------- get terms which contribute to mean pressure
c         careful with the sum, get the mean p_star pressure
c         CHECK if rotational form of the equations correct pxym see flat code
c
      do iz=izs,ize
         izm1 = iz - 1
         do iy=iys,iye
         do ix=1,nnx
            q_temp     =  0.5*(u(ix,iy,iz)**2 +
     +                         v(ix,iy,iz)**2 +
     +                         w(ix,iy,iz)**2)
            engz(iz)   = engz(iz) + q_temp
            engsbz(iz) = engsbz(iz) + e(ix,iy,iz)
            pxym(iz)   = pxym(iz) + (p(ix,iy,iz) - c23*e(ix,iy,iz))
         enddo
         enddo
         engz(iz)   = engz(iz)*fnxy
         engsbz(iz) = engsbz(iz)*fnxy
         pxym(iz)   = pxym(iz)*fnxy
      enddo
      call mpi_sum_z(engz(1),i_root,myid,nnzp1,1)
      call mpi_sum_z(engsbz(1),i_root,myid,nnzp1,1)
      call mpi_sum_z(pxym(1),i_root,myid,nnz,1)
c
c ------------ save means and divergence for printout and compmn
c              all cpus have means over all z
c
      if(istage .eq. 1) then
        do k=izs,ize
           do j=iys,iye
           do i=1,nnx
              fnt1(i,j,k) = vf(i,j,k)
           enddo
           enddo
        enddo
        call yd_mpi(fnt1(1,iys,izs),trigx(1,2),yk(1),
     +              nnx,nny,ixs,ixe,ix_s,ix_e,
     +              iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
        do k=izs,ize
           do j=iys,iye
           do i=1,nnx
              fnt2(i,j) = uf(i,j,k)
           enddo
           enddo
           call xderivp(fnt2(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
           km1 = k - 1
           do j=iys,iye
           do i=1,nnx
              divz(k) = divz(k) +
     +                  (fnt2(i,j) + fnt1(i,j,k) +
     +                  (wf(i,j,k) - wf(i,j,km1))*dzw_i(k))**2
           enddo
           enddo
           divz(k) = divz(k)*fnxy
        enddo
        call mpi_sum_z(divz(1),i_root,myid,nnz,1)
c
        do iz=1,nnz+1
           uxym(iz) = u_mn(iz)
           vxym(iz) = v_mn(iz)
           wxym(iz) = w_mn(iz)
        enddo
        do iscl=1,nscl
           do iz=1,nnz
              txym(iz,iscl) = t_mn(iz,iscl)
           enddo
        enddo
      endif
c
      return
      end
      subroutine spline(x,y,n,yp1,ypn,y2)
      integer n, nmax
      real yp1, ypn, x(n), y(n), y2(n)
      parameter (nmax=5000)
      integer i, k
      real p, qn, sig, un, u(nmax)
      if(yp1 .gt. .99e30) then
         y2(1) = 0.0
         u(1)  = 0.0
      else
         y2(1) = -0.5
         u(1) = (3./(x(2) - x(1)))*((y(2) - y(1))/(x(2) - x(1)) - yp1)
      endif
      do i=2,n-1
         sig = (x(i) - x(i-1))/(x(i+1) - x(i-1))
         p = sig*y2(i-1) + 2.0
         y2(i) = (sig - 1.0)/p
         u(i) = (6.0*((y(i+1) - y(i))/(x(i+1) - x(i)) - (y(i) - y(i-1))
     +          /(x(i) - x(i-1)))/(x(i+1) - x(i-1)) - sig*u(i-1))/p
      enddo
      if(ypn .gt. .99e+30) then
         qn = 0.0
         un = 0.0
      else
         qn = 0.5
         un = (3.0/(x(n) - x(n-1)))*
     +        (ypn - (y(n) - y(n-1))/(x(n) - x(n-1)))
      endif
      y2(n) = (un - qn*u(n-1))/(qn*y2(n-1) + 1.0)
      do k=n-1,1,-1
         y2(k) = y2(k)*y2(k+1) + u(k)
      enddo
c
      return
      end
      subroutine splint(xa,ya,y2a,n,x,y)
      integer n
      real x,y, xa(n), y2a(n), ya(n)
      integer k,khi,klo
      real a,b,h
      klo = 1
      khi = n
    1 continue
        if(khi - klo .gt. 1) then
           k = (khi + klo)/2
           if(xa(k) .gt. x) then
              khi = k
           else
              klo = k
           endif
           go to 1
      endif
      h = xa(khi) - xa(klo)
      a = (xa(khi) - x)/h
      b = (x - xa(klo))/h
      y = a*ya(klo) + b*ya(khi) +
     +    ((a**3 - a)*y2a(klo) + (b**3 - b)*y2a(khi))*(h**2)/6.0
c
      return
      end
      subroutine hx_deriv(ax_s,trigx,xk,nnx,ax)
c
c -------- get an x derivative saving the input. Use fftpack routines
c          with fftpack storage a0, (a1,b1), (a2,b2),...,an
c          assumes that first wavenumber xk(1) = 0.0
c
c          assumes that wavenumbers are normalized by number of points
c
      real xk(nnx), trigx(2*nnx+15), ax_s(nnx), ax(nnx)
c
      fn = 1.0/float(nnx)
      do i=1,nnx
         ax(i) = ax_s(i)*fn
      enddo
c
      call rfftf(nnx,ax(1),trigx)
c
      ii      = 1
      ax(1)   = 0.0
      ax(nnx) = 0.0
      do ix=2,nnx-1,2
         ii       = ii + 1
         temp     = ax(ix)
         ax(ix)   = -xk(ii)*ax(ix+1)
         ax(ix+1) = xk(ii)*temp
      enddo
      call rfftb(nnx,ax(1),trigx)
c
      return
      end
      subroutine gal_interp(f,trigx,nnx,xk,tb,grav,fi)

c
c ---------- propagate array f forward in space by distance dx
c            input is fourier coefficients in fftpack order
c
      real f(nnx), trigx(15+4*nnx), fi(nnx), xk(nnx), c_wav
c
      ncx = nnx/2 + 1
c
c --------- generate new coefficients
c
      cc    = 1.0
      ss    = 0.0
      fi(1) = f(1)*cc
      kk    = 1
      do i=2,nnx-1,2
         kk = kk + 1
         c_wav  = sqrt(grav/xk(i))
         if(i.gt.ncx) c_wav = sqrt(-grav/xk(i))
         cc = cos(xk(kk)*-c_wav*tb)
         ss = sin(xk(kk)*-c_wav*tb)
         aa = f(i)
         bb = f(i+1)
         fi(i)   = aa*cc - bb*ss
         fi(i+1) = bb*cc + aa*ss
      enddo
      cc      = cos(xk(ncx)*-c_wav*tb)
      ss      = 0.0
      aa      = f(nnx)
      fi(nnx) = aa*cc
c
c --------- transform back to get function at correct x points
c
      call rfftb(nnx,fi(1),trigx)
c
      return
      end
      subroutine hdot_interp(f,trigx,nnx,xk,tb,grav,fi)
c
c ----------- get hdot assuming packet moves with speed cs
c
      real f(nnx), trigx(15+4*nnx), fi(nnx), xk(nnx), c_wav
c
      ncx = nnx/2 + 1
c
c --------- generate modified fourier coefficients
c
      fi(1) = 0.0
      kk    = 1
      do i=2,nnx-1,2
         kk = kk + 1
         c_wav  = sqrt(grav/xk(i))
         if(i.gt.ncx) c_wav = sqrt(-grav/xk(i))
         cc = cos(xk(kk)*-c_wav*tb)
         ss = sin(xk(kk)*-c_wav*tb)
         aa = f(i)
         bb = f(i+1)
         uk = xk(kk)*c_wav
         fi(i)     = uk*(aa*ss + bb*cc)
         fi(i+1)   = uk*(bb*ss - aa*cc)
      enddo
      fi(nnx) = 0.0
c
c --------- transform back to get function at correct x points
c
      call rfftb(nnx,fi(1),trigx)
c
      return
      end
      subroutine xderivp(ax,trigx,xk,nnx,iys,iye)
c
c -------- get multiple x derivatives using fftpack routines
c          use fftpack storage a0, (a1,b1), (a2,b2),...,an
c          assumes that first wavenumber xk(1) = 0.0
c
c          assumes that wavenumbers are normalized by number of points
c
      real xk(nnx), trigx(2*nnx+15), ax(nnx,iys:iye)
c
c     fn = 1.0/float(nnx)
      do iy=iys,iye
         call rfftf(nnx,ax(1,iy),trigx)
         ii = 1
         ax(1,iy) = 0.0
         ax(nnx,iy) = 0.0
         do ix=2,nnx-1,2
            ii          = ii + 1
            temp        = ax(ix,iy)
            ax(ix,iy)   = -xk(ii)*ax(ix+1,iy)
            ax(ix+1,iy) = xk(ii)*temp
         enddo
         call rfftb(nnx,ax(1,iy),trigx)
      enddo
c
      return
      end
      subroutine fft2d_mpi(ax,at,trigx,trigc,nx,ny,
     +           jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           iz1,iz2,myid,ncpu,np,isgn)
c
c -------- get 2d fft using fftpack routines and parallel mpi
c          use fftpack storage a0, (a1,b1), (a2,b2),...,
c
c         isgn = -1 do forward transform, get coefficients
c                   incoming array is ax(nx+2,iys:iye,iz1:iz2)
c                   outgoing array is ax(nx+2,iys:iye,iz1:iz2)
c
c         isgn = -2 do forward transform, get coefficients
c                   incoming array is ax(nx+2,iys:iye,iz1:iz2)
c                   outgoing array is at(ny,jxs:jxe,iz1:iz2)
c
c         isgn =  1 do inverse transform, move to physical space
c                   incoming array is ax(nx+2,iys:iye,iz1:iz2)
c                   outgoing array is ax(nx+2,iys:iye,iz1:iz2)
c
c         isgn =  2 do inverse transform, move to physical space
c                   incoming array is at(ny,jxs:jxe,iz1:iz2)
c                   outgoing array is ax(nx+2,iys:iye,iz1:iz2)
c
      real ax(nx+2,iys:iye,iz1:iz2), at(ny,jxs:jxe,iz1:iz2),
     +     trigx(2*nx+15), trigc(4*ny+15),
     +     a2d(2,ny), a_wrk(nx)
      integer jx_s(0:np-1), jx_e(0:np-1),
     +        iy_s(0:np-1), iy_e(0:np-1)
c
      nxp2 = nx + 2
      if(isgn .lt. 0) then
         fn   = 1.0/(float(nx)*float(ny))
c
c ------ 1d fft in x over [iys,iye] for all z
c
         do iz=iz1,iz2
            do iy=iys,iye
               do ix=1,nx
                  a_wrk(ix) = ax(ix,iy,iz)*fn
               enddo
               call rfftf(nx,a_wrk(1),trigx(1))
               ax(1,iy,iz) = a_wrk(1)
               ax(2,iy,iz) = 0.0
               do ix=2,nx
                  ax(ix+1,iy,iz) = a_wrk(ix)
               enddo
               ax(nx+2,iy,iz) = 0.0
            enddo
         enddo
         call xtoy_trans(ax,at,nxp2,ny,jxs,jxe,jx_s,jx_e,
     +        iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
c
c ------ 1d fft in y over [jxs,jxe] for all z
c
         do iz=iz1,iz2
            do ix=jxs,jxe,2
               do iy=1,ny
                  a2d(1,iy) = at(iy,ix,iz)
                  a2d(2,iy) = at(iy,ix+1,iz)
               enddo
               call cfftf(ny,a2d(1,1),trigc(1))
               do iy=1,ny
                  at(iy,ix,iz)   = a2d(1,iy)
                  at(iy,ix+1,iz) = a2d(2,iy)
               enddo
            enddo
         enddo
c
c ---- decide whether to transpose back or leave as is
c
         if(isgn .eq. -1) then
            call ytox_trans(at,ax,nxp2,ny,jxs,jxe,jx_s,jx_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
         endif
c
      else
c
c ---- decide whether to first transpose or leave as is
c
         if(isgn .eq. 1) then
            call xtoy_trans(ax,at,nxp2,ny,jxs,jxe,jx_s,jx_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
         endif
c
c ------ 1d fft in y over [jxs,jxe] for all z
c
         do iz=iz1,iz2
            do ix=jxs,jxe,2
               do iy=1,ny
                  a2d(1,iy) = at(iy,ix,iz)
                  a2d(2,iy) = at(iy,ix+1,iz)
               enddo
               call cfftb(ny,a2d(1,1),trigc(1))
               do iy=1,ny
                  at(iy,ix,iz)   = a2d(1,iy)
                  at(iy,ix+1,iz) = a2d(2,iy)
               enddo
            enddo
         enddo
         call ytox_trans(at,ax,nxp2,ny,jxs,jxe,jx_s,jx_e,
     +        iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
c
c ------  1d fft in x over [iys,iye] for all z
c
         do iz=iz1,iz2
            do iy=iys,iye
               a_wrk(1) = ax(1,iy,iz)
               do ix=2,nx
                  a_wrk(ix) = ax(ix+1,iy,iz)
               enddo
               call rfftb(nx,a_wrk(1),trigx(1))
               do ix=1,nx
                  ax(ix,iy,iz) = a_wrk(ix)
               enddo
            enddo
         enddo
      endif
c
      return
      end
      subroutine fft2d_cmplx_mpi(ax,at,trigx,trigy,nx,nx2,ny,
     +           i2xs,i2xe,i2x_s,i2x_e,iys,iye,iy_s,iy_e,
     +           iz1,iz2,myid,ncpu,np)
c
c -------- get 2d inverse fft using complex fftpack routines
c          move from fourier space to physical space. specifically
c          designed for flow over waves with full complex wave spectrum
c
      real trigx(4*nx+15), trigy(4*ny+15), ax(nx2,iys:iye,iz1:iz2)
      real a2d(2,ny), at(ny,i2xs:i2xe,iz1:iz2)
      integer i2x_s(0:np-1), i2x_e(0:np-1), iy_s(0:np-1), iy_e(0:np-1)
c
c --------- transform in x direction
c
      do k=iz1,iz2
      do j=iys,iye
         call cfftb(nx,ax(1,j,k),trigx(1))
      enddo
      enddo
c
c --------- transform in y direction
c
      call xtoy_trans(ax,at,nx2,ny,i2xs,i2xe,i2x_s,i2x_e,
     +         iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
c
      do k=iz1,iz2
      do i=i2xs,i2xe,2
         do j=1,ny
            a2d(1,j) = at(j,i,k)
            a2d(2,j) = at(j,i+1,k)
         enddo
         call cfftb(ny,a2d(1,1),trigy(1))
         do j=1,ny
            at(j,i,k)   = a2d(1,j)
            at(j,i+1,k) = a2d(2,j)
         enddo
      enddo
      enddo
c
      call ytox_trans(at,ax,nx2,ny,i2xs,i2xe,i2x_s,i2x_e,
     +         iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
c
c     do iy=1,ny
c        call cfftb(nx,ax(1,1,iy),trigx(1))
c     enddo
c
      return
      end
      subroutine yderiv(ay,trigy,yk,nnx,nny)
c
c -------- get multiple y derivatives using fftpack routines
c          use fftpack storage a_0, (a1,b1), (a2,b2), ...,
c          assumes that first wavenumber yk(1) = 0.0
c
c          assumes that wavenumbers are normalized by number of points
c
      real yk(nny), trigy(2*nny+15), ay(nnx,nny)
      real a_trans(nny)
c
c     fn = 1.0/float(nny)
      do ix=1,nnx
         do iy=1,nny
            a_trans(iy) = ay(ix,iy)
         enddo
         call rfftf(nny,a_trans(1),trigy)
         ii = 1
         a_trans(1)   = 0.0
         a_trans(nny) = 0.0
         do iy=2,nny-1,2
            ii            = ii + 1
            temp          = a_trans(iy)
            a_trans(iy)   = -yk(ii)*a_trans(iy+1)
            a_trans(iy+1) = yk(ii)*temp
         enddo
         call rfftb(nny,a_trans(1),trigy)
         do iy=1,nny
            ay(ix,iy) = a_trans(iy)
         enddo
      enddo
c
      return
      end
      subroutine yd_mpi(ay,trigy,yk,
     +           nx,ny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
c
c -------- get multiple y derivatives using fftpack routines and mpi
c          use fftpack storage a_0, (a1,b1), (a2,b2), ...,
c          assumes that first wavenumber yk(1) = 0.0
c          wavenumbers are normalized by number of points, ny
c
      real yk(ny), trigy(2*ny+15), ay(nx,iys:iye,iz1:iz2)
      real ayt(ny,ixs:ixe,iz1:iz2)
c
      integer ix_s(0:np-1), ix_e(0:np-1),
     +        iy_s(0:np-1), iy_e(0:np-1)
c
      call xtoy_trans(ay,ayt,nx,ny,ixs,ixe,ix_s,ix_e,
     +         iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
c
c     fn = 1.0/float(nny)
      do iz=iz1,iz2
         do ix=ixs,ixe
            call rfftf(ny,ayt(1,ix,iz),trigy)
            ii = 1
            ayt(1,ix,iz)  = 0.0
            ayt(ny,ix,iz) = 0.0
            do iy=2,ny-1,2
               ii              = ii + 1
               temp            = ayt(iy,ix,iz)
               ayt(iy,ix,iz)   = -yk(ii)*ayt(iy+1,ix,iz)
               ayt(iy+1,ix,iz) = yk(ii)*temp
            enddo
            call rfftb(ny,ayt(1,ix,iz),trigy)
         enddo
      enddo
      call ytox_trans(ayt,ay,nx,ny,ixs,ixe,ix_s,ix_e,
     +         iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
c
      return
      end
      function rlim(d1,d2,d3)
c
c ------------- Cees's kappa=1/3 scheme
c
      r = (d1-d2+1.e-100)/(d2-d3-1.e-100)
      rlim = (d2-d3)*amax1(0.,amin1(r,amin1(1./6.+1./3.*r,1.)))
c
c ------------- Cees's kappa=-1 scheme
c
c     r = (d1-d2+1.e-100)/(d2-d3-1.e-100)
c     rlim = (d2-d3)*amin1(abs(r),0.5)
c
c ------------- first order upwind
c
c     rlim = 0.0
c
c ------------- QUICK scheme
c
c     rlim = -0.25*d2 - 0.125*d3 + 0.375*d1
c
      return
      end
      function ran1(idum)
c
c ----------- stolen from numerical recipes,p. 271
c
      integer idum, ia, im, iq, ir, ntab, ndiv
      real ran1, am, eps, rnmx
      parameter (ia=16807,im=2147483647,am=1.0/im,iq=127773,ir=2836.0,
     +           ntab=32,ndiv=1+(im-1)/ntab,eps=1.2e-07,rnmx=1.0-eps)
      integer j, k, iv(ntab), iy
      save iv, iy
      data iv /ntab*0/, iy /0/
      if(idum .le. 0 .or. iy .eq. 0) then
         idum = max(-idum,1)
         do j=ntab+8,1,-1
            k = idum/iq
            idum = ia*(idum - k*iq) - ir*k
            if(idum .lt. 0) idum = idum + im
            if(j .le. ntab) iv(j) = idum
         enddo
         iy = iv(1)
      endif
      k     = idum/iq
      idum  = ia*(idum - k*iq) - ir*k
      if(idum .lt. 0) idum = idum + im
      j     = 1 + iy/ndiv
      iy    = iv(j)
      iv(j) = idum
      ran1  = min(am*iy, rnmx)
c
      return
      end
      function ranf()
      data inc /1/
      save inc, ix, ia, m, fm
      if(inc.eq.1) then
        inc = 2
        m = 2**20
        fm = float(m)
        ix = 566387
        ia = 2**10 + 3
      endif
      ix = mod(ia*ix,m)
      fx = float(ix)
      ranf = fx/fm
      return
      end
      subroutine stokesv
c
c ----------- get stokes drift velocity for assumed wavelength stokesw
c             and wave amplitude stokesa. Changed sign for z.
c
c
      use pars
      use con_data
      use con_stats
      include 'mpif.h'
c
      if(iocean .eq. 1) then
c
c ----------- compute stokes velocity for ocean pbls
c
c        stokesw = pi2/20.0
         stokesw = pi2/76.5
c        ak      = 0.04
         ak      = 0.00
c        stokesa = 1.0
         stokesa = ak/stokesw
         sigma = sqrt(abs(grav)*stokesw)
         stokess = sigma*stokesw*stokesa**2
         do iz=1,nnzp1
            stokes(iz) = stokess*exp(2.0*stokesw*zz(iz))
         enddo
         if(l_root) then
            write(6,6000) (iz,zz(iz),stokes(iz),iz=1,nnz)
 6000       format(' iz ',10x,' zz',10x,' stokes',/,(1x,i3,2e12.4))
         endif
c
      else
c
c ----------------- set stokes velocity = 0 for atmos. pbls
c
         do iz=1,nnzp1
            stokes(iz) = 0.0
         enddo
         stokess = 0.0
         udrift = 0.0
         vdrift = 0.0
      endif
c
      return
      end
      subroutine busngr(zeta,phim,phis,psim,psis)
c
c ---- Businger's version of similarity theory
c
      data pih /1.57079633/
      save pih
c
      if(zeta .lt. 0.) then
         x=(1.0 - 15.0*zeta)**0.25
         phim = 1.0/x
         psim = 2.0*alog((1.0+x)/2.0) + alog((1.0+x*x)/2.0) -
     +          2.0*atan(x)+pih
         if(psim.gt.2.0)psim=2.0
         y = sqrt(1.0-9.0*zeta)
         phis = 0.74/y
         psis = alog((1.0+y)/2.0)*2.0
      else if(zeta .gt. 0) then
         phim = 1.0 + 4.7*zeta
         phis = 0.74 + 4.7*zeta
         psim = -4.7*zeta
         psis = -4.7*zeta
      else
         phim = 1.0
         phis = 0.74
         psim = 0.0
         psis = 0.0
      endif
      return
      end
      subroutine fzol(zeta,phim,phis,psim,psis)
c        estimate the stability functions for momentum, m
c                                         and scalars,  c
c        from input of the stability parameter zeta = z/L

      data c1/5./
      data a3,b3,a4,b4/1.258,8.382,-28.862,98.9545/
      data zetam,zetas/-0.2,-1.0/
      save c1, a3, b3, a4, b4, zetam, zetas
c
      psimu(Y)  = 1.571 + 2.0*(alog(0.5*(1.0 + Y)) - atan(Y)) +
     +            alog(0.5 + 0.5*Y**2)
      psisu(Y)  = 2.0*alog(0.5 + 0.5*Y)
      psicu(Y,G)= (1.0 - G)*alog(abs(Y - 1.0))
     +          + 0.5*(G + 2.0)*alog(abs(Y**2 + Y + 1.0))
     +          - (2.0*G + 1.0) / sqrt(3.0) *
     +            atan((Y + 0.5)*2.0/sqrt(3.0))
      Xm(zol)   = (1.0 - 16.0*zol)**0.25
      Xs(zol)   = sqrt(1.0 - 16.0*zol)
      Xc(zol,f) =  abs(1.0 - f*zol)**(4.0/3.0)/(1.0 - f*zol)
c
      if(zeta.ge.0.0)       then
c                                          STABLE
      if(zeta.le.1.0) then
        phim = 1.0 + c1 * zeta
        psim = - c1 * zeta
        phis = phim
        psis = psim
                      else
c                                   use limiting form
        phim = c1 + zeta
        psim = (1.0 - c1)*(1.0 + alog(zeta) ) - zeta
        phis = phim
        psis = psim
                      endif

                            else
c                                         UNSTABLE
c                                                  momentum
       if(zeta.ge.zetam) then
         phim = 1.0 / Xm(zeta)
         psim = psimu(Xm(zeta))
                         else
c                            use convective limit for momentum
         X = (1.0 - b3/a3 * zeta)**(1.0/3.0)

         fm = a3**(-1.0/3.0)
         phim = fm / Xc(zeta,b3/a3)
         psim = psimu(Xm(zetam))
     *        + psicu(Xc(zeta,b3/a3),fm)
     *        - psicu(Xc(zetam,b3/a3),fm)
                         endif

c                                         UNSTABLE scalars
       if(zeta.ge.zetas) then
         phis = 1.0/Xs(zeta)
         psis = psisu(Xs(zeta))
                         else
c                              use convective limit for scalars
         fs =   abs(a4)**(-1.0/3.0)*abs(a4)/a4
         phis = (a4 - b4*zeta)**(-1.0/3.0)
         psis = psisu(Xs(zetas))
     *        + psicu(Xc(zeta,b4/a4),fs)
     *        - psicu(Xc(zetas,b4/a4),fs)
                         endif

                            endif
       return
       end
      subroutine suft(it)
c
c ---------- iterate for zeta = z/L using bisection method
c            either businger or large functions can be specified
c
c            isfc = 0, specified surface heat flux
c                 = 1, specified surface temperature
c
      use pars
      use fields
      use con_data
      use con_stats
      real buf(3+nscl)
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      parameter (iter_mo = 30, zeta_min = -6.0, zeta_max = 3.0)
c
c ---------- limiting value for wind
c
      ufree = 0.07*(abs(batag*qstar(1)*dzw(1)))**(1./3.)
c
c ---- save old utau
c
      utausv = utau
      utau2  = utau*utau
c
      iz   = 1
      izp1 = iz + 1
      izm1 = iz - 1
c
      buf(1)  = 0.0
      buf(2)  = 0.0
      buf(3)  = 0.0
      tol     = 0.01
      do iy=iys,iye
      do ix=1,nnx
         buf(1) = buf(1) + u(ix,iy,iz)
         buf(2) = buf(2) + v(ix,iy,iz)
         wind(ix,iy) = sqrt((u(ix,iy,iz)+ugal)**2
     +                    +v(ix,iy,iz)*v(ix,iy,iz))
         buf(3) = buf(3) + wind(ix,iy)
      enddo
      enddo
      do iscl=1,nscl
         buf(3+iscl) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            buf(3+iscl) = buf(3+iscl) + t(ix,iy,iscl,iz)
         enddo
         enddo
      enddo
c
c -------- get x-y slab sums
c
      call mpi_sum_xy(buf,myid,iss,ise,(3+nscl))
      u1xy  = buf(1)*fnxy + ugal
      v1xy  = buf(2)*fnxy
      windm = buf(3)*fnxy
      do iscl=1,nscl
         t1xy(iscl) = buf(3+iscl)*fnxy
      enddo
      vsfc  = sqrt(u1xy*u1xy+v1xy*v1xy)
      windm = amax1(windm,ufree)
      vsfc  = amax1(vsfc,ufree)
c
c ---------- limits for zeta
c
      zeta_mn = zeta_min
      zeta_mx = zeta_max
      if(isfc(1) .eq. 0) then
         f_con = z1*batag*vk*qstar(1)/((windm*vk)**3)
      else
         d_theta = vk74in*(tsfcc(1) - t1xy(1))
         f_con   = z1*batag*vk*d_theta/((windm*vk)**2)
      endif
c
c --------- iteration for zeta
c
      do iter=1,iter_mo
         zeta_a = 0.5*(zeta_mn + zeta_mx)
         if(ismlt .eq. 1) then
             call busngr(zeta_a,phim,phis,psim,psis)
         else
             call fzol(zeta_a,phim,phis,psim,psis)
         endif
         u_fac = (zody - psim)
         if(isfc(1) .eq. 0) then
            f_new =  zeta_a + f_con*u_fac**3
         else
            t_fac = 1.0/(zody - psis)
            f_new =  zeta_a + f_con*u_fac*u_fac*t_fac
         endif
         if(f_new .lt. 0.0) then
            zeta_mn = zeta_a
         else
            zeta_mx = zeta_a
         endif
c
c ----------- iteration details
c
c        utau      = windm*vk/(zody-psim)
c        write(nprt,1000) iter, zeta_a, utau, phim, psim
c1000    format(' 1000 iter = ',i5,' zeta = ',e15.6,' u_* = ',e15.6,
c    +          ' phim = ',e15.6,' psim = ',e15.6)
      enddo
c
c --------- check if neutral surface layer
c
      if (ibuoy.eq.0 .or. qstar(1) .eq. 0.) then
          amonin    = 1000.
          zeta      = 0.0
          utau      = windm*vk/zody
          thstar(1) = 0.0
          t10xy(1)  = 0.0
          tsfcc(1)  = t1xy(1)
      else
         utau = windm*vk/(zody-psim)
         dnom = (zody-psis)*vk74in
         if(isfc(1) .eq. 0) then
            thstar(1) = -qstar(1)/utau
            tsfcc(1)  = t1xy(1)-thstar(1)*dnom
            t10xy(1)  = thstar(1)*dnom
         else
            thstar(1) = (t1xy(1) - tsfcc(1))/dnom
            t10xy(1)  = thstar(1)*dnom
            qstar(1)  = -utau*thstar(1)
         endif
         amonin = -utau**3/(batagk*qstar(1))
         zeta   = z1/amonin
      endif
c
c ------- surface details, for debug
c
c     write(nprt,2000) windm, utau, qstar(1), tsfcc(1), amonin, zeta,
c    +              z1, batag, vk, batagk, zo
c2000 format(' 2000 suft ',/,
c    +       '    windm = ',e15.6,' utau = ',e15.6,' qstar = ',e15.6,/,
c    +       '    tsfcc = ',e15.6,' MO L = ',e15.6,' z1/L = ',e15.6,/,
c    +       '    z1 = ',e15.6,' batag = ',e15.6,' vk = ',e15.6,/,
c    +       '    batagk = ',e15.6,' zo = ',e15.6)
c
      if (utau.gt.10.0) then
         write(6,9000)
         write(6,9200) utau,windm
         go to 9999
      endif
      if (t10xy(1).gt.0. .and. qstar(1) .gt. 0.) then
         write(6,9000)
         write(6,9300) u1xy,v1xy,t1xy(1),
     +                 tsfcc(1),amonin,utau,it
         go to 9999
      endif

      !Second scalar:
      if(isfc(2) .eq. 0) then
         thstar(2) = -qstar(2)/utau
         tsfcc(2)  = t1xy(2)-thstar(2)*dnom
         t10xy(2)  = thstar(2)*dnom
      else
         thstar(2) = (t1xy(2) - tsfcc(2))/dnom
         t10xy(2)  = thstar(2)*dnom
         qstar(2)  = -utau*thstar(2)
      endif

c
c ---------- examples of two other scalars
c
c     c
c     c **** get flux of b scalar, specified surface value
c     c
c           dnom      = (zody-psis)*vk74in
c           thstar(2) = (t1xy(2)-tsfcc(2))/dnom
c           qstar(2)  = -thstar(2)*utau
c           t10xy(2)  = thstar(2)*dnom
c           aut3m(2)  =  qstar(2)
c
c **** get surface value of c scalar, specified surface flux
c
c     dnom      = (zody-psis)*vk74in
c     thstar(2) = -qstar(2)/utau
c     tsfcc(2)  = t1xy(2) - dnom*thstar(2)
c     t10xy(2)  = thstar(2)*dnom
c     aut3m(2)  = qstar(2)
c
      zol = zeta
      hol = zol*zi/z1
c
c ---- note roundoff problem in angles if close to multiples of pi
c
      tep = u1xy/vsfc
      if(tep.gt.1.)  tep = 1.0
      if(tep.lt.-1.) tep = -1.0
      thta      = acos(tep)
      utau2     = utau*utau
      au13m     = -utau2*cos(thta)
      au23m     = -utau2*sin(thta)*sign(1.,v1xy)
      aut3m(1)  =  qstar(1)
c
      return
c
c -------- iteration did not converge
c
 9999 continue
 9000 format(' Trouble in SR. suft')
 9200 format(' Stop because utau = ',e15.6,' windm = ',e15.6)
 9300 format(' ** CHECK SFC U = ',e15.6,' V=',e15.6,' T,TS = ',2e15.6,
     +       ' L =',e15.6,' U_* = ',e15.6,' AT IT = ',i5)
      call mpi_abort(mpi_comm_world,icode,ierr)
      end
      subroutine sufto(it)
c
      use pars
      use fields
      use con_data
      use con_stats
      real buf(3+nscl)
c
c ------- version of similarity theory adpated for ocean flows
c      option to use businger or large version of similarity theory
c
      iz    = 1
      izm1  = iz - 1
      izp1  = iz + 1
      z1_a  = abs(z1)
      buf(1)  = 0.0
      buf(2)  = 0.0
      buf(3)  = 0.0
      tol     = 0.01
      do iy=iys,iye
      do ix=1,nnx
         buf(1) = buf(1) + u(ix,iy,iz)
         buf(2) = buf(2) + v(ix,iy,iz)
         wind(ix,iy) = sqrt((u(ix,iy,iz)+ugal)**2
     +                    +v(ix,iy,iz)*v(ix,iy,iz))
         buf(3) = buf(3) + wind(ix,iy)
      enddo
      enddo
      do iscl=1,nscl
         buf(3+iscl) = 0.0
         do iy=iys,iye
         do ix=1,nnx
            buf(3+iscl) = buf(3+iscl) + t(ix,iy,iscl,iz)
         enddo
         enddo
      enddo
c
c -------- get x-y slab sums
c
      call mpi_sum_xy(buf,myid,iss,ise,(3+nscl))
      u1xy  = buf(1)*fnxy + ugal
      v1xy  = buf(2)*fnxy
      windm = buf(3)*fnxy
      do iscl=1,nscl
         t1xy(iscl) = buf(3+iscl)*fnxy
      enddo
      vsfc  = sqrt(u1xy*u1xy+v1xy*v1xy)
      windm = amax1(windm,ufree)
      vsfc  = amax1(vsfc,ufree)
c
      t10xy(1)=-qstar(1)/utau*zody*vk74in
c
c ---- check for temperature boundary condition
c
      if(isfc(1) .eq. 0 ) then
         tsfcc(1)=t1xy(1)-t10xy(1)
      endif
c     vsfc=sqrt(u1xy*u1xy+v1xy*v1xy)
c     if(windm.le.0.01)windm=0.01
c     if(vsfc .le.0.01)vsfc =0.01
c
c ----------- input surface wind stress (tau = 0.0184n/m*m)
c             density rho = 1000kg/m^3
c
c     utau = 4.29e-03
c     utau = 6.10e-03
      utau = 7.00e-03
c
c **** save old utau
      utausv = utau
      utau2  = utau*utau
      if (ibuoy.eq.0 .or. qstar(1) .eq. 0.) then
          amonin    = 1000.
          zeta      = 0.
          thstar(1) = 0.0
          t10xy(1)  = 0.0
      else
          amonin = -utau2*utau/(batagk*qstar(1))
          zeta   = z1_a/amonin
      endif
      if (t10xy(1).lt.0. .and. qstar(1) .lt. 0.) then
         write(6,1234)u1xy,v1xy,t1xy(1),tsfcc(1),amonin,utau,it
 1234    format(' ** check sfc u=',e12.3,' v=',e12.3,' t,ts=',2f10.3,
     +     ' l=',e12.3,' u*=',e12.3,' at it=',i5)
         go to 9999
      endif
c
c -------- for stable,neutral and unstable pbl get drift velocity
c
      if(ismlt .eq. 1) then
          call busngr(zeta,phim,phis,psim,psis)
      else
          call fzol(zeta,phim,phis,psim,psis)
      endif
      udrift = windm + stokes(1) - stokess + utau*(zody-psim)*vkin
      vdrift = 0.0
      dnom      = (zody-psis)*vk74in
      if (isfc(1).eq.1) then
         thstar(1) = (t1xy(1) - tsfcc(1))/dnom
         t10xy(1)  = thstar(1)*dnom
         qstar(1)  = - utau*thstar(1)
      else
         thstar(1)  = -qstar(1)/utau
         tsfcc(1)   = t1xy(1)-thstar(1)*dnom
         t10xy(1)   = thstar(1)*dnom
      endif
      zol = zeta
      hol = zol*zi/z1
c
c ---------- examples of two other scalars
c
c     c
c     c **** get flux of b scalar, specified surface value
c     c
c           dnom      = (zody-psis)*vk74in
c           thstar(2) = (t1xy(2)-tsfcc(2))/dnom
c           qstar(2)  = -thstar(2)*utau
c           t10xy(2)  = thstar(2)*dnom
c           aut3m(2)  =  qstar(2)
c     c
c     c **** get surface value of c scalar, specified surface flux
c     c
c           dnom      = (zody-psis)*vk74in
c           thstar(3) = -qstar(3)/utau
c           tsfcc(3)  = t1xy(3) - dnom*thstar(3)
c           t10xy(3)  = thstar(3)*dnom
c           aut3m(3)  = qstar(3)
c
c **** note roundoff problem in angles are close to multiples of pi
c     tep=u1xy/vsfc
c     if(tep.gt.1.)tep=1.
c     if(tep.lt.-1.)tep=-1.
c     thta=acos(tep)
      utau2 = utau*utau
c     au13m = -utau2*cos(thta)
c     au23m = -utau2*sin(thta)*sign(1.,v1xy)
      au13m = utau2
      au23m = 0.0
      aut3m(1)= qstar(1)
c
      return
c
c --------- trouble in sl routine
c
 9999 continue
c
      write(nprt,9000)
 9000 format(' Trouble in SR. sufto')
      call mpi_finalize(ierr)
      stop
      end
      subroutine suft2(u_level1,it)
c
      use pars
      use fields
      use con_data
      use con_stats
c
      real u_level1(nnx,iys:iye,2+nscl)
c
c     u_level1(.,.,1) = u
c     u_level1(.,.,2) = v
c     u_level1(.,.,3) = theta
c     u_level1(.,.,4) = more scalars
c
      tol = 0.01
      ufree=0.07*(abs(batag*qstar(1)*dzw(1)))**(1./3.)
      zeta_mn = -6.0
      zeta_mn_i = 1.0/zeta_mn
      iz   = 1
c     izm1 = iz - 1
c     izp1 = iz + 1
c
c      write(nprt,3131) myid, utau, zody, vk74in, batagk,
c    +               u_level1(jxs,1,3), u_level1(jxe,1,3),
c    +               u_level1(jxs,1,1), u_level1(jxe,1,1),
c    +               u_level1(jxs,1,2), u_level1(jxe,1,2)
c3131  format(' in suft2 myid = ',i4,/,
c    +        ' utau = ',e15.6,' zody = ',e15.6,/,
c    +        ' vk74in = ',e15.6,' batagk = ',e15.6,/,
c    +        ' t(jxs) = ',e15.6,' t(jxe) = ',e15.6,/,
c    +        ' u(jxs) = ',e15.6,' u(jxe) = ',e15.6,/,
c    +        ' v(jxs) = ',e15.6,' v(jxe) = ',e15.6)
c
      do iy=iys,iye
      do ix=mxs,mxe
c
c ----------------- first guess for utau
c
      utau = .001
c
      t10xy(1) = -qstar(1)/utau*zody*vk74in
      tsfcc(1) = u_level1(ix,iy,3) - t10xy(1)
      vsfc2    = u_level1(ix,iy,1)**2 + u_level1(ix,iy,2)**2
      vsfc     = sqrt(vsfc2)
      windm    = ufree+vsfc
      utausv   = utau
      utau2    = utau*utau
      amonin   = -utau2*utau/(batagk*qstar(1))
      if(amonin.eq.0.) then
            write(6,5050) ix,iy,it,utau,amonin
 5050       format(' 5050, sr. suft2, trouble at ',/,
     +             ' ix = ',i6,'iy = ',i6,' it = ',i6,' utau = ',e15.6,
     +             ' amonin = ',e15.6)
            stop
      endif
c
c ---- for unstable, free convection pbl
c
      iter = 0
 100  continue
c
c ----------------- limit the min (-l/z) change to accmmodate stable flow
c
      zeta_i = amin1(amonin/z1,zeta_mn_i)
      zeta_a = 1.0/zeta_i
c
      if(ismlt .eq. 1) then
          call busngr(zeta_a,phim,phis,psim,psis)
      else
          call fzol(zeta_a,phim,phis,psim,psis)
      endif
      utau     = windm*vk/(zody-psim)
      thstar(1)=-qstar(1)/utau
      amonold  = amonin
      amonin   = utau*utau/(batagk*thstar(1))
      diff     = abs(amonin - amonold)
c      write(nprt,5656)iter,psim,utau,zeta,amonin,dmonin,diff
c 5656 format(' iter=',i4,' phm=',e10.3,' utau=',e10.3,
c     1      ' zeta=',e10.3,' l=',e10.3,' diff = ',e12.4)
      iter = iter+1
      if(iter.gt.10)go to 1000
      if(diff.gt.abs(tol*amonin)) go to 100
 1000 continue
c
 2000 continue
c
      if (utau.gt.10.) then
         write(6,232)utau,windm
  232    format(' stop because utau=',e15.6,' windm=',e15.6)
         stop 9999
      endif
      t10xy(1) = -qstar(1)/utau*vk74in*(zody-psis)
      t_grnd(ix,iy,1) = u_level1(ix,iy,3) - t10xy(1)
c
      zol = zeta_a
      hol = zol*zi/z1
      tep = u_level1(ix,iy,1)/windm
      if(tep.gt.1.)  tep = 1.0
      if(tep.lt.-1.) tep = -1.0
      thta  = acos(tep)
      utau2 = utau*utau
c     au13m=-utau2*cos(thta)
c     au23m=-utau2*sin(thta)*sign(1.,u_level1(ix,iy,2))
c     aut3m(1)= qstar(1)
c
      tau13m(ix,iy)   = -utau2*cos(thta)
      tau23m(ix,iy)   = -utau2*sin(thta)*sign(1.,u_level1(ix,iy,2))
      taut3m(ix,iy,1) = qstar(1)
c
c **** get surface value of c scalar, specified surface flux
c
c     dnom      = (zody-psis)*vk74in
c     thstar(2) = -qstar(2)/utau
c     tsfcc(2)  = u_level1(ix,iy,4) - dnom*thstar(2)
c     t_grnd(ix,iy,2)  = u_level1(ix,iy,4) - dnom*thstar(2)
c     t10xy(2)  = thstar(2)*dnom
c     taut3m(ix,iy,2)  = qstar(2)
c
c
c ------- end of x-y loops
c
      enddo
      enddo
c
      return
      end
      subroutine init
c
      use pars
      use particles
      use fields
      use con_data
      use con_stats
c

      real :: RHT,RHB,topC,botC

      pi   = 4.0*atan(1.0)
      pi2  = 2.0*pi
      bfac = 1.0
      if(ibuoy.eq.0) bfac = 0.
c
c -------------------- case specific data
c
      if(iocean .eq. 1) then
         t00b    = 5000.0
         cp      = 4.20e03
         gcp     = grav/cp
         batag   = bfac*grav/t00b
         fcor_h  = 0.0
         vgcont  = 0.
c        wtsfc(1)=0.00
c        wtsfc(1)=4.96e-07
         wtsfc(1)=1.190476e-06
         qstar(1)=wtsfc(1)
c        dtdzf(1)=0.000
         dtdzf(1)=0.2548
         dtjump  = 0.
         divgls  = 0.
         zo      = 0.0001
         zi      = -5.
c        izi     = (55*nnz)/100
c        izi     = nnz
         izi     = 55
         xl      = 50.
         yl      = 50.
         zl      = -20.
c
c ---------- if stretched grid specify location of first point
c
      else
         cp      = 1.e3
         gcp     = grav/cp
         batag   = bfac*grav/t00
         fcor_h  = 0.0


c
c ----------- wind tunnel simulation
c
c

         !Convert RH at bottom (given by scalar (2)) to specific humidity:
         RHbot = tsfcc(2)  !Store the RH for use later
         botC=(tsfcc(1)-273.15)
         RHB = tsfcc(2)/100.0

         tsfcc(2) =  RHB*Mw/Ru/tsfcc(1)*610.94*
     +               EXP(17.6257*botC/(botC+243.04))/rhoa



         wtsfc(1)=qstar(1)
         wtsfc(2)=qstar(2)
c
         dtdzf(1)=0.000
         dtdzf(2)=0.000
c
c
c ----------
c
         dtjump  = 0.0
         divgls  = 0.0

c
c        izi     = int(float(nnz)*zi/zl)
c        izi     = (60*nnz)/100
         izi     = 102
c
c ------------ check type of grid generation, adjust box size
c              so that you get the desired size in physical space
c
         if(iz_space .eq. 4) then
           zl = zl + zw1
         endif
      endif
c
      time     = 0.0
      t_stage  = 0.0
      t_stage2 = 0.0
      t_zero   = 0.0
c
c ---------- set mesh stationary flag
c
      imesh = 0
c
c ---------- outermost coarse grid  indicies are bounds of grid
c
      izlow = 1
      izup  = nnz
      dz    = zl/nnz
      dzg   = abs(dz)
      if(l_root) write(6,4040) zl,nnz,dzg
c
c --------------- generate z grids for particular mesh from
c                 iz = 0,1,...,nnz+1; this allows indexing
c                 to array elements z(0), etc.
c
      zwstrt = 0.0
c
c ------------ if uniform vertical spacing then
c
      if(iz_space .eq. 0) then
c
c ------------ build z grid for w points
c
         do iz=0,nnz+1
            z(iz) = dz*float(iz) + zwstrt
         enddo
      else
        z_gues = 1.1
        !call vgrid(zw1,zi,zl,z_gues,nnz,z(0),l_root,l_debug)
        call vgrid_channel(zw1,zi,zl,nnz,z(0),l_root,l_debug)
      endif
c
      call get_dz
c
      if(l_root) then
         write(6,8002) zwstrt
         write(6,8003) (iz,z(iz),zz(iz),iz=0,nnz+1)
      endif
c
      nnzm1 = nnz-1
      dx    = xl/nnx
      dy    = yl/nny
      fnxy  = 1./float(nxy)
      dzdz  = dzw(1)*dzw(1)
      z1    = zz(1)
c
      c23  = 2.0/3.0
      dsl  = (dx*1.5*dy*1.5*abs(dzw(1)))**(1./3.)
      dslg = dsl
      cs   = 0.2
c
      vk     = 0.4
      batagk = batag*vk
      vkin   = 1./vk
      ttmean = 0.
      zody   = alog(abs(z1/zo))
      write(nprt, 9901) z1,zo,zody
 9901 format(' 9901 z1 = ',e15.6,' zo = ',e15.6,/,
     +       ' zody = ',e15.6)
      zodyin = 1./zody
      wstar  = abs(batag*zi*wtsfc(1))**(1./3.)
      if(ismlt .eq. 1) then
c
c ---- set constants for businger similarity functions
c
         vk74   = vk*0.74
         vk74in = 0.74/vk
         zody74 = zody*0.74
      else
c
c ---- set constants for large similarity functions
c
        vk74    = vk
        vk74in  = 1.0/vk
        zody74  = zody
      endif
      ugal   = 0.0
c     ugal   = ugcont*0.5
c     ugcont = ugcont - ugal
      cdbtm  = vk*vk/zody/zody
      if(iocean .eq. 1) then
c ----------- set surface friction velocity here and in sr. sufto
c        utau = 4.29e-03
         utau = 7.00e-03
      else
         ufree = 0.07*(abs(batag*qstar(1)*dzw(1)))**(1./3.)
c
c ---- note : new estimate for utau !!!
c
         utau  = vk*(ufree+ugcont)/zody
c        utau  = vk*(ufree)/zody
      endif
      utau2    = utau*utau
      if(ibuoy .eq. 0 .or. qstar(1) .eq. 0.) then
        amonin = 1000.0
      else
        amonin = -utau2*utau/(batagk*qstar(1))
      endif
      hol   = abs(zi)/amonin
      zol   = abs(z1)/amonin
      uwsfc = -utau*utau
      vwsfc = -utau*utau
c
c ------- make sure tsfcc is gt than t00 for both isfc=0 or 1
c
      if(l_root) then
         write(6,80)
         write(6,2)wtsfc(1),utau,amonin,dtdzf(1),zody,zo
     +         ,cdbtm,ugcont
      endif
c
      if(l_debug) then
         write(nprt,80)
         write(nprt,2)wtsfc(1),utau,amonin,dtdzf(1),zody,zo
     +         ,cdbtm,ugcont
      endif
c
      return
c ------------------------
   2  format(10x,' WT =',e12.4,',  U* =',e12.4,',  L =',e12.4,/,
     +       10x,' DTDZ FREE =',e12.4,',  ZODY=',e12.4,/,10x,
     +       ' ZO(BTM) =',e12.4,',  CDBTM=',e12.4,
     +       ',  UG = ',e12.4)
  80  format(///,' ***** SCRATCH RUN ***** ',//)
 4040 format(' zl = ',e15.6,' nnz = ',i5,' dzg = ',e15.6)
 4043 format(' znest = ',e15.6,' nnz = ',i5,' dzg = ',e15.6)
 8002 format(' zwstrt = ',e12.4)
 8003 format(' iz ',5x,' zw',5x,' zu ',5x,/,(i3,2e12.4))
      end
      subroutine cgrid(zl,nnz,z)
c
c --------- uniform vertical grid
c
      real z(0:nnz+1)
c
      ddz = zl/float(nnz)
      do k=1,nnz
         z(k) = float(k)*ddz
      enddo
      z(nnz) = zl
      z(0)   = 0.0
      z(nnz+1) = z(nnz) + (z(nnz) - z(nnz-1))
c
      return
      end
      subroutine surf_setup_gabls
c
c ------- set basic parameters of gabls surface and evaluate pieces that are
c         independent of time
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use wavy_pars
c
c --------- allocate surface arrays
c
      allocate(eta0(2,nnx,nny))
c
      wave = xl/4.0
      fac  = pi2/wave
      ak   = 0.3
      ampl = ak*wave/pi2
      dx   = xl/float(nnx)
      do j=iys,iye
      do i=1,nnx
         eta0(1,i,j)  = ampl*cos(fac*dx*float(i-1))
         eta0(2,i,j)  = 0.0
         u_orbit(i,j) = 0.0
         v_orbit(i,j) = 0.0
      enddo
      enddo
c
      return
      end
      subroutine get_surf_gabls(t_bndy)
c
c ------- generate new gabls surface at t = t_bndy. smoothly add or subtract
c         bumps by varying the amplitude
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use wavy_pars
c
      t1_st  = 100.0
      t1_en  = 400.0
      b_rate = 1.0/(t1_en - t1_st)
c
      if(t_bndy .lt. t1_st) then
         ampl_fac  = 1.0
         ampl_rate = 0.0
      elseif(t_bndy .ge. t1_st .and. t_bndy .lt. t1_en) then
         ampl_fac  = amax1((1.0 - b_rate*(t_bndy - t1_st)),0.0)
         ampl_rate = -b_rate
      elseif(t_bndy .ge. t1_en) then
         ampl_fac  = 0.0
         ampl_rate = 0.0
      endif
c
      do j=iys,iye
      do i=1,nnx
         bndy(i,j)    = ampl_fac*eta0(1,i,j)
         hdot(i,j)    = ampl_rate*eta0(1,i,j)
         w_orbit(i,j) = hdot(i,j)
      enddo
      enddo
c
      return
      end
      subroutine surf_setup
c
c ------- set basic parameters of surface and evaluate pieces that are
c         independent of time
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use wavy_pars
c
      real amp(nnx,nny), ran_amp(2,nnx,nny)
c
c --------- allocate surface arrays
c
      allocate(omega(nnx,nny),
     +         eta0(2,nnx,nny),
     +         s_wrk(2,nnx,iys:iye),
     +         v_wrk(2,nnx,iys:iye),
     +         etax(nnxp2,iys:iye),
     +         etay(nnxp2,iys:iye),
     +         eta_wrk(nny,i2xs:i2xe))
c
      u10    = 15.0
      age    = 1.2
      angle  = atan(0.0)
      wind_x = cos(angle)
      wind_y = sin(angle)
      dkx    = pi2/xl
      dky    = pi2/yl
c
c ------------- generate an array of random phases or amplitudes
c               and radial frequencies
c
      do j=1,nny
      do i=1,nnx
         omega(i,j) = sqrt(sqrt(xkn(i)**2 + ykn(j)**2)*grav)
      enddo
      enddo
      idum = -1
      if(i_ranp .eq. 0) then
         do j=1,nny
         do i=1,nnx
            aa_phase       = pi2*(ran1(idum) - 0.5)
            ran_amp(1,i,j) = cos(aa_phase)
            ran_amp(2,i,j) = sin(aa_phase)
         enddo
         enddo
      else
c
c ------- normalize amplitudes so that their squared sum = 1
c
         fac = 1.0/sqrt(2.0)
         do j=1,nny
         do i=1,nnx
            ran_amp(1,i,j) =  gasdev(idum)*fac
            ran_amp(2,i,j) =  gasdev(idum)*fac
         enddo
         enddo
         sum1 = 0.0
         sum2 = 0.0
         do j=1,nny
         do i=1,nnx
            sum1 = sum1 + ran_amp(1,i,j)
            sum2 = sum2 + ran_amp(2,i,j)
         enddo
         enddo
         sum1 = sum1*fnxy
         sum2 = sum2*fnxy
         var1 = 0.0
         var2 = 0.0
         do j=1,nny
         do i=1,nnx
            var1 = var1 + (ran_amp(1,i,j) - sum1)**2
            var2 = var2 + (ran_amp(2,i,j) - sum2)**2
         enddo
         enddo
         var1 = var1*fnxy
         var2 = var2*fnxy
         write(nprt,1345) sum1,var1, sum2,var2
 1345    format(' 1345 get surf:',/,
     +          '      sum1 = ',e15.6,' var1 = ',e15.6,/,
     +          '      sum2 = ',e15.6,' var2 = ',e15.6)
      endif
c
      sum1  = 0.0
      sum2  = 0.0
      p_max = 0.0
      p_min = 0.0
      do j=1,nny
      do i=1,nnx
         sum1 = sum1 + ran_amp(1,i,j)
         sum2 = sum2 + ran_amp(2,i,j)
         if(ran_amp(1,i,j) .gt. p_max) p_max = ran_amp(1,i,j)
         if(ran_amp(1,i,j) .lt. p_min) p_max = ran_amp(1,i,j)
      enddo
      enddo
      sum1 = sum1*fnxy
      sum2 = sum2*fnxy
      write(nprt,5000) p_max, p_min, sum1, sum2
 5000 format(' phase max = ',e15.6,/,
     +       ' phase min = ',e15.6,/,
     +       ' sum r     = ',e15.6,/,
     +       ' sum c     = ',e15.6)
c
c --------- first construct the pieces independent of t
c
      fac_a = sqrt(2.0)
      if(i_pm .eq. 1) then
c
c ---------- pierson-moskowitz spectrum
c
         call pm64(dkx,dky,xkn,ykn,nnx,nny,u10,amp)
c
      else
c
c ---------- donelan spectrum
c
         call donelan(dkx,dky,xkn,ykn,nnx,nny,u10,age,angle,amp)
c
      endif
c
c ------- use a gaussian filter to smooth the wave shape
c
      fac     = 3.0
      f_width = sqrt(fac*dx*fac*dy)
      cons    = f_width*f_width/24.0
      sum_spec = 0.0
      do j=1,nny
      do i=1,nnx
         gauss       = exp(-cons*(xkn(i)**2 + ykn(j)**2))
         amp_root    = sqrt(amp(i,j)*dkx*dky)*gauss
         eta0(1,i,j) = fac_a*amp_root*ran_amp(1,i,j)
         eta0(2,i,j) = fac_a*amp_root*ran_amp(2,i,j)
         sum_spec    = sum_spec + amp(i,j)*dkx*dky
      enddo
      enddo
c
      write(6,8001) sum_spec
 8001 format(' 8001 sum_spec = ',e15.6)
c
c --------- build coefficients at t = 0
c
      do j=1,nny
         eta0(1,ncx,j) = 0.0
         eta0(2,ncx,j) = 0.0
      enddo
      do i=1,nnx
         eta0(1,i,ncy) = 0.0
         eta0(2,i,ncy) = 0.0
      enddo
c
c --------- make surface smooth
c
c     mx = nnx/4
c     mx = 3*(nnx/8)
c     my = 3*(nny/8)
c
c     mcx   = mx/2 + 1
c
      m_cut_x = nnx/3
      m_cut_y = nny/3
c
      i_cut_l = m_cut_x
      i_cut_u = nnx + 2 - i_cut_l
      do j=1,nny
      do i=i_cut_l,i_cut_u
         eta0(1,i,j) = 0.0
         eta0(2,i,j) = 0.0
      enddo
      enddo
      j_cut_l = m_cut_y
      j_cut_u = nny + 2 - j_cut_l
      do j=j_cut_l,j_cut_u
      do i=1,nnx
         eta0(1,i,j) = 0.0
         eta0(2,i,j) = 0.0
      enddo
      enddo
c
c ------------ set amplitudes = 0 for iflat case
c
c     do j=1,nny
c     do i=1,nnx
c        eta0(1,i,j) = 0.0
c        eta0(2,i,j) = 0.0
c     enddo
c     enddo
c
      return
      end
      subroutine surf_setup_tank
c
c ------- basic parameters of wavy surface from lab experiments
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use wavy_pars
      use tank_pars
c
      real xx(nx_tank), zz_tank(nx_tank), us(nx_tank), ws(nx_tank),
     +     y2(nx_tank), u2(nx_tank), w2(nx_tank), xles(nnx)

c
c --------- allocate surface arrays
c
      allocate(z_tank(nnx), u_tank(nnx), w_tank(nnx),
     +         z_tank_f(nnx), u_tank_f(nnx), w_tank_f(nnx),
     +         z_tank_i(nnx), u_tank_i(nnx), w_tank_i(nnx),
     +         hdot_i(nnx), hx_i(nnx))
c
      lut    = 89
      close(lut)
      if(i_tank == 1) then     ! no breaking case
         open(lut,file='./les_tank_nobreak.dat')
         npts_t  = 1485
         wave_l  = 0.232780E+02*0.01
         c_speed = 0.602861E+00
         wave_k  = pi2/wave_l
         ustar_tank = 0.381 ! from table 2 of banner 1990
      elseif(i_tank == 2) then ! breaking case
         open(lut,file='./les_tank_break.dat')
         npts_t  = 1485
         wave_l  = 0.235351E+02*0.01
         c_speed = 0.606335E+00
         wave_k  = pi2/wave_l
         ustar_tank = 0.4928 ! from table 2 of banner 1990
      elseif(i_tank == 3) then     ! mono wave case
         open(lut,file='./les_tank_mono.dat')
         npts_t  = 1500
         wave_l  = 0.232780
         c_speed = 0.602861E+00
         wave_k  = pi2/wave_l
         ustar_tank = 0.381 ! from table 2 of banner 1990
      elseif(i_tank == 6) then     ! mono wave case with no drift
         open(lut,file='./les_tank_mono_3.dat')
         npts_t  = 1500
         wave_l  = 0.232780E+02*0.01
         c_speed = 0.602861E+00
         wave_k  = pi2/wave_l
         ustar_tank = 0.381 ! from table 2 of banner 1990
      elseif(i_tank == 7) then     ! no breaking case with drift = 0
         open(lut,file='./les_tank_nobreak_orb.dat')
         npts_t  = 1485
         wave_l  = 0.232780E+02*0.01
         c_speed = 0.602861E+00
         wave_k  = pi2/wave_l
         ustar_tank = 0.381 ! from table 2 of banner 1990
      elseif(i_tank == 8) then     ! mono wave case ak = 0.1 with no drift
         open(lut,file='./les_tank_mono_4.dat')
         npts_t  = 1500
         wave_l  = 0.232780E+02*0.01
         c_speed = 0.602861E+00
         wave_k  = pi2/wave_l
         ustar_tank = 0.381 ! from table 2 of banner 1990
      elseif(i_tank == 9) then ! breaking case with no drift
         open(lut,file='./les_tank_break_nodrift.dat')
         npts_t  = 1485
         wave_l  = 0.235351E+02*0.01
         c_speed = 0.606335E+00
         wave_k  = pi2/wave_l
         ustar_tank = 0.4928 ! from table 2 of banner 1990
      elseif(i_tank == 10) then !7/12/18: Particle-laden monochromatic turbulence
         open(lut,file='./monochromatic.dat')
         npts_t  = 1000
         wave_k  = pi2/wave_l
      elseif(i_tank == 11) then !7/12/18: Particle-laden monochromatic turbulence with ak=0.2
         open(lut,file='./monochromatic_ak2.dat')
         npts_t  = 1000
         wave_k  = pi2/wave_l
      elseif(i_tank == 12) then !7/12/18: Particle-laden monochromatic turbulence with ak=0.2
         open(lut,file='./flat.dat')
         npts_t  = 1000
         wave_k  = pi2/wave_l
      elseif(i_tank == 13) then !7/12/18: Particle-laden real wave
         open(lut,file='./polychromatic.dat')
         npts_t  = maxnx
         wave_k  = pi2/wave_l
      endif
c
      do i=1,4
         read(lut,'(a1)') idum
      enddo
      do i=1,npts_t
         read(lut,*) xx(i), zz_tank(i), us(i), ws(i)
      enddo
      if (myid==0) write(*,*) 'Finished reading tank data'
c
c ---------- change lengths to m from cm
c
!Commenting out -- BANNER DATA IN CM, MONOCHROMATIC IN M!
!      do i=1,npts_t
!        xx(i)      = xx(i)*0.01
!        zz_tank(i) = zz_tank(i)*0.01
!        us(i)      = us(i)
!        ws(i)      = ws(i)
!     enddo
c
!Commenting out the nondimensionalizations:
!      c_speed = c_speed/ustar_tank
!      wave_l  = 1.0
!      wave_k  = pi2
c
      write(nprt,1001) (i,xx(i), zz_tank(i),us(i),ws(i),
     +                  i=1,nnx)
 1001 format(' 1001: tank ',/,
     +       ' i ',5x,' xx ',5x,' zz ',5x,' us ',5x,
     +       ' ws ',/,(i5,4e15.6))
c
c ---------- fit a spline to the data
c
      yp1 = 2.0e+30
      ypn = yp1
      ! call spline(xx,zz_tank,npts_t,yp1,ypn,y2)
      ! call spline(xx,us,npts_t,yp1,ypn,u2)
      ! call spline(xx,ws,npts_t,yp1,ypn,w2)
c
c ---------- interpolate into the spline function at les points
c
      ddx = xl/float(nnx)
      do i=1,nnx
         xles(i) = float(i-1)*ddx
         ! call splint(xx,zz_tank,y2,npts_t,xles(i),z_tank(i))
         ! call splint(xx,us,u2,npts_t,xles(i),u_tank(i))
         ! call splint(xx,ws,w2,npts_t,xles(i),w_tank(i))
         z_tank(i) = zz_tank(i)
         u_tank(i) = us(i)
         w_tank(i) = ws(i)
      enddo
c
c     write(nprt,4002) (xx(i),zz_tank(i),i=1,npts_t)
c4002 format(' xx ',5x,' zz tank ',/,(2e15.6))
c     write(nprt,4004) (xx(i),ws(i),i=1,npts_t)
c4004 format(' xx ',5x,' ws tank ',/,(2e15.6))
c     write(nprt,4003) (xles(i), z_tank(i),i=1,nnx)
c4003 format(' xles ',5x,' zz les ',/,(2e15.6))
c     write(nprt,4005) (xles(i), w_tank(i),i=1,nnx)
c4005 format(' xles ',5x,' w les ',/,(2e15.6))
c
      sum = 0.0
      do i=1,nnx
         sum = sum + z_tank(i)
      enddo
      sum = sum/float(nnx)
      do i=1,nnx
         z_tank(i) = z_tank(i) - sum
      enddo
c
c --------- force mean w = 0 from the data
c
      sum_w = 0.0
      do i=1,nnx
         sum_w = sum_w + w_tank(i)
      enddo
      sum_w = sum_w/float(nnx)
      do i=1,nnx
         w_tank(i) = w_tank(i) - sum_w
      enddo
c
      write(nprt,4001) sum, sum_w
 4001 format(' Mean wave height = ',e15.6,' sum w = ',e15.6)
c
c --------- get fourier coefficients at t = 0
c           to use in building shape at t > 0
c
      fnx = 1.0/float(nnx)
      do i=1,nnx
         z_tank_f(i) = z_tank(i)*fnx
         u_tank_f(i) = u_tank(i)*fnx
         w_tank_f(i) = w_tank(i)*fnx
      enddo
      call rfftf(nnx,z_tank_f(1),trigx(1,1))
      call rfftf(nnx,u_tank_f(1),trigx(1,1))
      call rfftf(nnx,w_tank_f(1),trigx(1,1))
c
c ------ dealias the inputs to match the grid resolution
c        careful with the indices, match fftpack order
c
      ix_cut   = 2*int(float(nnx)/3.) + 2
c     ix_cut   = 44
c
      do i=ix_cut,nnx
         z_tank_f(i) = 0.0
         u_tank_f(i) = 0.0
         w_tank_f(i) = 0.0
      enddo
c
c ------ DEBUG
c
c     do i=1,nnx
c        z_tank_i(i) = z_tank(i)
c        u_tank_i(i) = u_tank(i)
c        w_tank_i(i) = w_tank(i)
c     enddo
c
      return
      end
      subroutine get_surf(t_bndy)
c
c ------- generate new surface at t = t_bndy
c         assuming its the real part of the complex fft
c         allow for time varying amplitude to smoothly add wave
c         spectrum from a restart
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use wavy_pars
      real bndyx(nnx,iys:iye), bndyy(nnx,iys:iye,izs:izs)
c
c ------------ for doing restart with gradually growing bumps
c              specify initial and end times to blend over
c
      t1_st  = 0.59781897E+04 + 0.5 ! time for volume ja1062
      t1_en  = t1_st + 400.0
      b_rate = 1.0/(t1_en - t1_st)
c
c --------- for debug turn on waves right away
c
      if(t_bndy .lt. t1_st) then
         ampl_fac  = 0.0
         ampl_rate = 0.0
      elseif(t_bndy .ge. t1_st .and. t_bndy .lt. t1_en) then
         ampl_fac  = amax1((b_rate*(t_bndy - t1_st)),0.0)
         ampl_rate = b_rate
      elseif(t_bndy .ge. t1_en) then
         ampl_fac  = 1.0
         ampl_rate = 0.0
      endif
c
c -------------- make surface flat, for initialization
c
c     ampl_fac  = 0.0
c     ampl_rate = 0.0
c
      nx2 = 2*nnx
      do j=iys,iye
      do i=1,nnx
         aa           = eta0(1,i,j)
         bb           = eta0(2,i,j)
         comg         = cos(omega(i,j)*t_bndy)
         somg         = sin(omega(i,j)*t_bndy)
         s_wrk(1,i,j) = comg*aa + somg*bb
         s_wrk(2,i,j) = comg*bb - somg*aa
      enddo
      enddo
c
c     j = iys
c     write(nprt,7001) (i, bndy(i,j),eta0(i,j),i=1,nnx)
c7001 format(' 7001 debug ',/,' i',5x,' bndy ',5x,' eta0',/,(i5,2e15.6))
c     write(nprt,7010) (i, omega(i,j),i=1,ncx)
c7010 format(' 7010 debug ',/,' i',5x,' omega ',/,(i5,e15.6))
c
c --------- get the u_orbital velocity of the wave
c
      do j=iys,iye
      do i=2,nnx
         fac_o        = (grav/omega(i,j))*ampl_fac
         v_wrk(1,i,j) = fac_o*xkn(i)*s_wrk(1,i,j)
         v_wrk(2,i,j) = fac_o*xkn(i)*s_wrk(2,i,j)
      enddo
      enddo
      i  = 1
      do j=max(2,iys),iye
         fac_o        = (grav/omega(i,j))*ampl_fac
         v_wrk(1,i,j) = fac_o*xkn(i)*s_wrk(1,i,j)
         v_wrk(2,i,j) = fac_o*xkn(i)*s_wrk(2,i,j)
      enddo
      if(iys .eq. 1) then
         v_wrk(1,i,iys) = 0.0
         v_wrk(2,i,iys) = 0.0
      endif
      call fft2d_cmplx_mpi(v_wrk(1,1,iys),eta_wrk(1,i2xs),trigcx,trigc,
     +           nnx,nx2,nny,i2xs,i2xe,i2x_s,i2x_e,iys,iye,iy_s,iy_e,
     +           1,1,myid,ncpu_s,numprocs)
      do j=iys,iye
      do i=1,nnx
         u_orbit(i,j) = v_wrk(1,i,j)
      enddo
      enddo
c
c --------- get the v_orbital velocity of the wave
c
      do j=iys,iye
      do i=2,nnx
         fac_o        = (grav/omega(i,j))*ampl_fac
         v_wrk(1,i,j) = fac_o*ykn(j)*s_wrk(1,i,j)
         v_wrk(2,i,j) = fac_o*ykn(j)*s_wrk(2,i,j)
      enddo
      enddo
      i  = 1
      do j=max(2,iys),iye
         fac_o        = (grav/omega(i,j))*ampl_fac
         v_wrk(1,i,j) = fac_o*ykn(j)*s_wrk(1,i,j)
         v_wrk(2,i,j) = fac_o*ykn(j)*s_wrk(2,i,j)
      enddo
      if(iys .eq. 1) then
         v_wrk(1,i,iys) = 0.0
         v_wrk(2,i,iys) = 0.0
      endif
      call fft2d_cmplx_mpi(v_wrk(1,1,iys),eta_wrk(1,i2xs),trigcx,trigc,
     +           nnx,nx2,nny,i2xs,i2xe,i2x_s,i2x_e,iys,iye,iy_s,iy_e,
     +           1,1,myid,ncpu_s,numprocs)
      do j=iys,iye
      do i=1,nnx
         v_orbit(i,j) = v_wrk(1,i,j)
      enddo
      enddo
c
c ------------- get hdot
c
      do j=iys,iye
      do i=2,nnx
         v_wrk(1,i,j) = omega(i,j)*s_wrk(2,i,j)*ampl_fac +
     +                  ampl_rate*s_wrk(1,i,j)
         v_wrk(2,i,j) = -omega(i,j)*s_wrk(1,i,j)*ampl_fac +
     +                  ampl_rate*s_wrk(2,i,j)
      enddo
      enddo
      i  = 1
      do j=max(2,iys),iye
         v_wrk(1,i,j) = omega(i,j)*s_wrk(2,i,j)*ampl_fac +
     +                  ampl_rate*s_wrk(1,i,j)
         v_wrk(2,i,j) = -omega(i,j)*s_wrk(1,i,j)*ampl_fac +
     +                  ampl_rate*s_wrk(2,i,j)
      enddo
      if(iys .eq. 1) then
         v_wrk(1,i,iys) = 0.0
         v_wrk(2,i,iys) = 0.0
      endif
      call fft2d_cmplx_mpi(v_wrk(1,1,iys),eta_wrk(1,i2xs),trigcx,trigc,
     +           nnx,nx2,nny,i2xs,i2xe,i2x_s,i2x_e,iys,iye,iy_s,iy_e,
     +           1,1,myid,ncpu_s,numprocs)
      do j=iys,iye
      do i=1,nnx
         hdot(i,j) = v_wrk(1,i,j)
      enddo
      enddo
c
c -------- get the surface
c
      call fft2d_cmplx_mpi(s_wrk(1,1,iys),eta_wrk(1,i2xs),trigcx,trigc,
     +           nnx,nx2,nny,i2xs,i2xe,i2x_s,i2x_e,iys,iye,iy_s,iy_e,
     +           1,1,myid,ncpu_s,numprocs)
      do j=iys,iye
      do i=1,nnx
         bndy(i,j)          = s_wrk(1,i,j)*ampl_fac
         bndyx(i,j)         = bndy(i,j)
         bndyy(i,j,izs:izs) = bndy(i,j)
      enddo
      enddo
c
c -------- get variance of surface
c
      sumb = 0.0
      do j=iys,iye
      do i=1,nnx
         sumb = sumb + bndy(i,j)
      enddo
      enddo
      sumb = sumb*fnxy
      call mpi_sum_xy(sumb,myid,iss,ise,1)
c
      varb = 0.0
      do j=iys,iye
      do i=1,nnx
         varb = varb + (bndy(i,j) - sumb)**2
      enddo
      enddo
      varb = varb*fnxy
      call mpi_sum_xy(varb,myid,iss,ise,1)
c     write(6,8102) sumb, varb
c8102 format(' 8102: sum bndy = ',e15.6,' variance bndy = ',e15.6)
c
c --------  w_o = hdot +  u_o*h_x + v_o*h_y
c
      call xderivp(bndyx(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
      call yd_mpi(bndyy(1,iys,izs),trigx(1,2),yk(1),
     +            nnx,nny,ixs,ixe,ix_s,ix_e,
     +            iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
      do j=iys,iye
      do i=1,nnx
         w_orbit(i,j) = hdot(i,j) + u_orbit(i,j)*bndyx(i,j) +
     +                  v_orbit(i,j)*bndyy(i,j,izs)
      enddo
      enddo
c
c ----------- make sure mean of w_orbit = 0
c
      sum_w = 0.0
      do j=iys,iye
      do i=1,nnx
         sum_w = sum_w + w_orbit(i,j)
      enddo
      enddo
      sum_w = sum_w*fnxy
      call mpi_sum_xy(sum_w,myid,iss,ise,1)
c
      do j=iys,iye
      do i=1,nnx
         w_orbit(i,j) = w_orbit(i,j) - sum_w
      enddo
      enddo
c
c     j = iys
c     write(nprt,8002) t_bndy
c8002 format(' t bndy = ',e15.6)
c     write(nprt,8001) (i,bndy(i,j),u_orbit(i,j),w_orbit(i,j),i=1,nnx)
c8001 format(' 8001 get surf',/,
c    +       ' i ',5x,' bndy ',5x,' u orb',5x,' w orb ',/,
c    +       (i5,3e15.6))
c
      return
      end
      subroutine get_surf_tank(t_bndy)
c
c ------- generate new surface at t = t_bndy
c         for the laboratory wave. Build w_o and hdot so
c         that they satisfy the inputs for u_o and h
c         with fix sr. hx_deriv
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use wavy_pars
      use tank_pars
c
c ------------ for doing restart with gradually growing bumps
c              specify initial and end times to blend over
c              more aggressive strategy to trigger turbulence
c
c     t1_st  = 0.85340788E+01 + 0.01 ! time from volume ub1020
c
      t1_st  =  0.0 ! time relative to t_zero set by i_move_step
c
c     t1_st  = 100000.0 ! to make surface always flat for initialization
c
      t1_en  = t1_st + 0.25
      b_rate = 1.0/(t1_en - t1_st)
c
c --------- for debug turn on waves right away
c
      if(t_bndy .lt. t1_st) then
         ampl_fac  = 0.0
         ampl_rate = 0.0
      elseif(t_bndy .ge. t1_st .and. t_bndy .lt. t1_en) then
         ampl_fac  = amax1((b_rate*(t_bndy - t1_st)),0.0)
         ampl_rate = b_rate
      elseif(t_bndy .ge. t1_en) then
         ampl_fac  = 1.0
         ampl_rate = 0.0
      endif
c
c ------------ get (u,w) components of surface velocity and wave
c              height. Make h, hdot, uo, wo match
c
c     dist = -c_speed*t_bndy
c     dist = 0
      call gal_interp(u_tank_f,trigx(1,1),nnx,xkn,t_bndy,grav,u_tank_i)
      !call gal_interp(w_tank_f,trigx(1,1),nnx,xkn,t_bndy,grav,w_tank_i)
      call gal_interp(z_tank_f,trigx(1,1),nnx,xkn,t_bndy,grav,z_tank_i)
      call hdot_interp(z_tank_f,trigx(1,1),nnx,xkn,t_bndy,grav,
     +    hdot_i)
      call hx_deriv(z_tank_i,trigx(1,1),xkn,nnx,hx_i)
c
c     write(nprt,1012) c_speed, t_bndy
c1012 format(' 1012: tank s_speed = ',e15.6,' t_bndy = ',e15.6)
c
      sum   = 0.0
      sum_w = 0.0
      do i=1,nnx
         sum = sum + z_tank_i(i)
      enddo
      sum = sum/float(nnx)
c     write(nprt,1013) sum
c1013 format(' 1013: tank sum z = ',e15.6)
c
      ampl_fac = 1.0
      do j=iys,iye
      do i=1,nnx
         u_orbit(i,j) = u_tank_i(i)*ampl_fac
         v_orbit(i,j) = 0.0
         bndy(i,j)    = z_tank_i(i)*ampl_fac
      enddo
      enddo
c
c     j = iys
c     write(nprt,7001) (i, bndy(i,j),u_orbit(i,j),w_orbit(i,j),i=1,nnx)
c7001 format(' 7001 debug ',/,' i',5x,' bndy ',5x,' u orbit',
c    +       5x,' w orbit',/,(i5,3e15.6))
c
c -------- get variance of surface
c
      sumb = 0.0
      do j=iys,iye
      do i=1,nnx
         sumb = sumb + bndy(i,j)
      enddo
      enddo
      sumb = sumb*fnxy
      call mpi_sum_xy(sumb,myid,iss,ise,1)
c
      varb = 0.0
      do j=iys,iye
      do i=1,nnx
         varb = varb + (bndy(i,j) - sumb)**2
      enddo
      enddo
      varb = varb*fnxy
      call mpi_sum_xy(varb,myid,iss,ise,1)
c     write(6,8102) sumb, varb
c8102 format(' 8102: sum bndy = ',e15.6,' variance bndy = ',e15.6)
c
c -------- use w_o = hdot + u_o*h_x + v_o*h_y
c
      do j=iys,iye
      do i=1,nnx
         hdot(i,j)    = hdot_i(i)*ampl_fac
         w_orbit(i,j) = hdot(i,j) + u_orbit(i,j)*hx_i(i)*ampl_fac
      enddo
      enddo
c
c ----------- make sure mean of w_orbit = 0
c
      sum_w = 0.0
      do j=iys,iye
      do i=1,nnx
         sum_w = sum_w + w_orbit(i,j)
      enddo
      enddo
      sum_w = sum_w*fnxy
      call mpi_sum_xy(sum_w,myid,iss,ise,1)
c
      do j=iys,iye
      do i=1,nnx
         w_orbit(i,j)    = w_orbit(i,j) - sum_w
      enddo
      enddo
c
c     j = iys
c     write(nprt,8002) t_bndy
c8002 format(' t bndy = ',e15.6)
c     write(nprt,8001) (i,bndy(i,j),u_orbit(i,j),w_orbit(i,j),i=1,nnx)
c8001 format(' 8001 get surf',/,
c    +       ' i ',5x,' bndy ',5x,' u orb',5x,' w orb ',/,
c    +       (i5,3e15.6))
c
      return
      end
      subroutine pm64(dkx,dky,xk,yk,nx,ny,u10,amp)
c
c --------- evaluate PM wave spectrum in terms of wavenumber
c
      real, parameter ::
     +      grav  = 9.81, f_pm4 = 0.13, a_pm4 = 0.0081, b_pm4 = 1.25
      real kmag
      real amp(nx,ny), xk(nx), yk(ny)
c
      pi2 = 8.0*atan(1.0)
c
      sum_dir  = 0.0
      do j=1,ny
      do i=1,nx
c
c ----------- convert wave number into frequency
c
         kmag = sqrt(xk(i)**2 + yk(j)**2)
         f    = sqrt(grav*kmag)/pi2
c
         if(f .le. 0.0) then
            amp(i,j) = 0.0
         else
            t1        = a_pm4*grav*grav/(f*(pi2*f)**4)
            fpm       = grav*f_pm4/u10
            t4        = b_pm4*(fpm/f)**4
            amp(i,j)  = t1*exp(-t4)
c
c ----------- convert into wavenumber
c
            amp(i,j)  = amp(i,j)*0.5*sqrt(grav/kmag)/(kmag*pi2)
c
c ----------- simple directional spectrum
c
            if(xk(i) .ge. 0.0 .and. yk(j) .ge. 0.0) then
              angle = acos(xk(i)/kmag)
            elseif(xk(i) .ge. 0.0 .and. yk(j) .le. 0.0) then
              angle =  -acos(xk(i)/kmag)
            elseif(xk(i) .le. 0.0 .and. yk(j) .ge. 0.0) then
              angle = pi2*0.5 - acos(abs(xk(i))/kmag)
            elseif(xk(i) .le. 0.0 .and. yk(j) .le. 0.0) then
              angle = -pi2*0.5 + acos(abs(xk(i))/kmag)
            endif
            if(angle .lt. -pi2/4.0 .or. angle .gt. pi2/4.0) then
              fac = 0.0
            else
              fac1     = cos(angle)
              fac      = fac1**8
            endif
            amp(i,j) = amp(i,j)*fac
            sum_dir  = sum_dir + fac*dkx*dky
         endif
c
      enddo
      enddo
c
      sum_dir_i = 1.0/sum_dir
      sum_spec  = 0.0
      do j=1,ny
      do i=1,nx
         amp(i,j) = amp(i,j)*sum_dir_i
         sum_spec = sum_spec + amp(i,j)*dkx*dky
      enddo
      enddo
c
      write(6,3000) sum_spec, sum_dir
 3000 format(' 3000: PM 64 sum_spec = ',e15.6,' sum_dir = ',e15.6)
c
      return
      end
      subroutine donelan(dkx,dky,xk,yk,nx,ny,u10,age,dir_mean,amp)
c
c -------- get amplitude of wavenumber spectrum according to donelan
c          prescription in Komen etal, p. 187 that includes
c          directional dependence. The routine converts the
c          frequency-directional spectrum F(omega,phi)
c          into a wavenumber-directional spectrum G(k,phi), see notes.
c
c     Inputs are :
c
c     (xk,yk) --- (kx,ky) wavenumbers
c     u10     --- wind speed at 10 m
c     age     --- wave age = cp/u10
c     dir_mean -- mean direction of wave field in radians from x-axis
c
c     Output :
c
c     amp --- amplitude of spectrum at wavenumber (xk,yk)
c
      real xk(nx), yk(ny), amp(nx,ny)
c
      data ionce /0/
      save grav, pi2, ionce
      real kmag
c
      if(ionce .eq. 0) then
        grav  = 9.81
        pi2   = 8.0*atan(1.0)
        ionce = 1
      endif
c
c     write(6,5601) u10, age, dir_mean
 5601 format(' 5601: u10 = ',e15.6,' age = ',e15.6,/,
     +       ' dir mean = ',e15.6)

c
      cp      = age*u10
      omega_p = grav/cp
      age_i   = 1.0/age
c
      aa      = amin1(age_i,5.0)
      aa      = amax1(aa, 1.0)
      sigma_d = 0.08*(1.0 + 4.0/(aa**3))
c
      aa      = amin1(age_i,5.0)
      aa      = amax1(aa, 0.83)
      alpha_d = 0.006*(aa**0.55)
c
      if(age_i .le. 1.0) then
        gamma_d = 1.7
      else
        aa = amin1(age_i,5.0)
        gamma_d = 1.7 + 6.0*alog10(age_i)
      endif
c
      sum_spec = 0.0
      sum_dir  = 0.0
c
      do j=1,ny
      do i=1,nx
c
         kmag    = sqrt(xk(i)**2 + yk(j)**2)
         omega   = sqrt(grav*kmag)
         if(kmag .eq. 0.0) then
           amp(i,j) = 0.0
           go to 999
         endif
c
         if(xk(i) .ge. 0.0 .and. yk(j) .ge. 0.0) then
           angle = acos(xk(i)/kmag)
         elseif(xk(i) .ge. 0.0 .and. yk(j) .le. 0.0) then
           angle =  -acos(xk(i)/kmag)
         elseif(xk(i) .le. 0.0 .and. yk(j) .ge. 0.0) then
           angle = pi2*0.5 - acos(abs(xk(i))/kmag)
         elseif(xk(i) .le. 0.0 .and. yk(j) .le. 0.0) then
           angle = -pi2*0.5 + acos(abs(xk(i))/kmag)
         endif
c
         gamma = exp(-(omega - omega_p)**2/
     +           (2.0*(sigma_d**2)*omega_p**2))
         fac   = exp(-(omega_p/omega)**4)
         fac1  = gamma_d**gamma
         spec  = alpha_d*(grav**2)*(omega/omega_p)*
     +           fac*fac1/(omega**5)
c
c --------- directional spreading
c
         ratio = omega/omega_p
         if(ratio .gt. 0.56 .and. ratio .lt. 0.95) then
           beta = 2.61*ratio**1.3
         elseif(ratio .ge. 0.95 .and. ratio .lt. 1.6) then
           beta = 2.28/ratio**1.3
         else
           beta = 1.24
         endif
c
         fac      = (angle - dir_mean)*beta
         sech     = 2.0/(exp(fac) + exp(-fac))
         spec_dir = beta*0.5*sech**2
c
         amp(i,j)  = spec*spec_dir
c
c --------- weight amp by the conversion factor to put it
c           into wavenumber space ... see notes
c
         amp(i,j) = amp(i,j)*0.5*sqrt(grav/kmag)/kmag
c
c --------- get integrals in wavenumber space
c
         sum_spec  = sum_spec + amp(i,j)*dkx*dky
         sum_dir   = sum_dir + spec_dir*dkx*dky
c
  999 continue
c
      enddo
      enddo
c
      write(6,6000) sum_spec, sum_dir
 6000 format(' 6000: in Donelan ',/,
     +       ' Sum spectrum = ',e15.6,/,
     +       ' Sum directional spectrum = ',e15.6)
c
c ---------- rescale the directional
c            spectrum to make sure its integral = 1
c
      sum_dir_i = 1.0/sum_dir
c
c     do j=1,ny
c     do i=1,nx
c        amp(i,j) = amp(i,j)*sum_dir_i
c     enddo
c     enddo
c
      return
      end

      real function gasdev(idum)
      integer idum, iset
      real fac, gset, rsq, v1, v2, ran1
      save iset, gset
      data iset /0/
      if(iset .eq. 0) then
    1 continue
        v1  = 2.0*ran1(idum) - 1.0
        v2  = 2.0*ran1(idum) - 1.0
        rsq = v1**2 + v2**2
        if(rsq .ge. 1.0 .or. rsq .eq. 0.0) go to 1
        fac    = sqrt(-2.0*log(rsq)/rsq)
        gset   = v1*fac
        gasdev = v2*fac
        iset   = 1
      else
        gasdev = gset
        iset   = 0
      endif
c
      return
      end
      subroutine l_bndy(t_bndy)
c
c ---------- build the shape of the lower boundary
c            and compute its speed of movement
c            at time t_bndy
c
      use pars
      use fields
      use con_data
      use con_stats
c
c -------- pick boundary shape
c
      ikase = 1 ! wavy bottom with 2 waves
      ikase = 4 ! 3-D hill
      ikase = 2 ! monochromatic wave
      ikase = 6 ! propagating monochromatic wave plus cross swell
      ikase = 7 ! stationary monochromatic wave , ak = 0.5
      ikase = 3 ! 2-D hill
      ikase = 5 ! propagating monochromatic wave
c
      if(ikase .eq. 1) then ! multiple waves
c
         wave = xl/4.0
         fac  = pi2/wave
         ak   = 0.1
         ampl = ak*wave/pi2
         wave2 = xl/8.0
         fac2  = pi2/wave2
         ak2   = 0.01
         ampl2 = ak2*wave2/pi2
         phase = 0.0
         do j=iys,iye
         do i=1,nnx
            bndy(i,j) = ampl*cos(fac*xgrid(i) + phase) +
     +                  ampl2*sin(fac2*xgrid(i))
            hdot(i,j) = 0.0
            u_orbit(i,j) = 0.0
            v_orbit(i,j) = 0.0
         enddo
         enddo
c
      else if(ikase .eq. 2) then ! monochromatic wave
c
         wave = xl/4.0
         fac  = pi2/wave
         ak   = 0.1
         ampl = ak*wave/pi2
         do j=iys,iye
         do i=1,nnx
            bndy(i,j) = ampl*cos(fac*xgrid(i))
            hdot(i,j) = 0.0
            u_orbit(i,j) = 0.0
            v_orbit(i,j) = 0.0
         enddo
         enddo
c
      elseif(ikase .eq. 3) then ! 2-D hill
c
         a      = 25.0
         xcntr  = 0.5*xl
         hill_l = 67.0
c
         a      = 25.0
         hill_l = 67.0
c
c --------- flat hill for generating wt flow
c
c        a      = 0.0
c
         do j=iys,iye
         do i=1,nnx
            if(xgrid(i) .ge. (xcntr - 2.0*hill_l) .and.
     +         xgrid(i) .le. (xcntr + 2.0*hill_l)) then
              x_temp    = xgrid(i) - xcntr
              bndy(i,j) = a*(1.0 + cos(pi2*x_temp/(4.0*hill_l)))
            else
              bndy(i,j) = 0.0
            endif
            hdot(i,j) = 0.0
            u_orbit(i,j) = 0.0
            v_orbit(i,j) = 0.0
         enddo
         enddo
c
      elseif(ikase .eq. 4) then ! 3-D hill
c
         xo    = xl/2.0
         yo    = yl/2.0
         ampl  = 60.000
         rlen  = 1000.0
         do j=iys,iye
         do i=1,nnx
            rad = sqrt((xgrid(i)-xo)**2 + (ygrid(j)-yo)**2)
            if(rad .le. rlen*0.5) then
               bndy(i,j) = ampl*(1.0 + cos(pi2*rad/rlen))
            else
               bndy(i,j) = 0.0
            endif
            hdot(i,j) = 0.0
            u_orbit(i,j) = 0.0
            v_orbit(i,j) = 0.0
         enddo
         enddo
c
      else if(ikase .eq. 5) then ! propagating monochromatic wave
c
         wave_l  = xl/12.0
         wave_k  = pi2/wave_l
         ak      = 0.1
         ampl    = ak*wave_l/pi2
         c_speed = sqrt(grav/wave_k)
c
c ---------- for start set ampl = 0
c
c        ampl = 0.0
c
c ---------- get Cartesian orbital velocities
c
         f_rad = sqrt(grav*sqrt(wave_k**2))
         con_u = grav*wave_k/f_rad
         do j=iys,iye
         do i=1,nnx
            bndy(i,j) = ampl*cos(wave_k*(xgrid(i) - c_speed*t_bndy))
            hdot(i,j) = ampl*wave_k*c_speed*
     +                    sin(wave_k*(xgrid(i) - c_speed*t_bndy))
            u_orbit(i,j) = bndy(i,j)*con_u
            v_orbit(i,j) = 0.0
         enddo
         enddo
c
      else if(ikase .eq. 6) then ! propagating monochromatic wave
                                 ! long wave moving vertically
c
         wave_l  = xl/12.0
         wave_k  = pi2/wave_l
         ak      = 0.1
         ampl    = ak*wave_l/pi2
         c_speed = sqrt(grav/wave_k)
c
         wave_l2 = xl/6.0
         wave_k2 = pi2/wave_l2
         ak2     = 0.025
         ampl2   = ak2*wave_l2/pi2
         c_spd2  = sqrt(grav/wave_k2)
         do j=iys,iye
         do i=1,nnx
            bndy(i,j) = ampl*cos(wave_k*(xgrid(i) - c_speed*t_bndy))
            hdot(i,j) = ampl*wave_k*c_speed*
     +                    sin(wave_k*(xgrid(i) - c_speed*t_bndy))
            bndy2     = ampl2*cos(wave_k2*(ygrid(j) - c_spd2*t_bndy))
            bndy(i,j) = bndy2 + bndy(i,j)
         enddo
         enddo
c
c ---------- get Cartesian orbital velocities
c
         f_rad = sqrt(grav*sqrt(wave_k**2 + wave_k2**2))
         con_u = grav*wave_k/f_rad
         con_v = grav*wave_k2/f_rad
         do j=iys,iye
         do i=1,nnx
            u_orbit(i,j) = bndy(i,j)*con_u
            v_orbit(i,j) = bndy(i,j)*con_v
         enddo
         enddo
c
c        i = nnx/4
c        j = (iys + iye)/2
c        write(nprt,3101) wave_l, c_speed, bndy(i,j),hdot(i,j)
c3101    format(' 3101 bndy ',/,
c    +   ' wave l = ',e15.6,' c = ',e15.6,
c    +   ' bndy = ',e15.6,' hdot = ',e15.6)
c
      else if(ikase .eq. 7) then ! stationary monochromatic wave
                                 ! long wave moving vertically
c
         wave_l  = xl/2.0
         wave_k  = pi2/wave_l
         ak      = 0.5
c
c -------- make it flat
c
         ak      = 0.0
         ampl    = ak*wave_l/pi2
         c_speed = sqrt(grav/wave_k)
         c_speed = 0.0
c
         do j=iys,iye
         do i=1,nnx
            bndy(i,j) = ampl*cos(wave_k*(xgrid(i)))
            hdot(i,j) = 0.0
         enddo
         enddo
c
c ---------- get Cartesian orbital velocities
c
         f_rad = sqrt(grav*sqrt(wave_k**2))
         con_u = grav*wave_k/f_rad
         con_u = 0.0
         do j=iys,iye
         do i=1,nnx
            u_orbit(i,j) = 0.0
            v_orbit(i,j) = 0.0
         enddo
         enddo
c
c        i = nnx/4
c        j = (iys + iye)/2
c        write(nprt,3101) wave_l, c_speed, bndy(i,j),hdot(i,j)
c3101    format(' 3101 bndy ',/,
c    +   ' wave l = ',e15.6,' c = ',e15.6,
c    +   ' bndy = ',e15.6,' hdot = ',e15.6)
c
c
      endif
c
c     j = iys
c     write(nprt,5001) (i,bndy(i,j),i=1,nnx)
c5001 format(' i ',5x,' bndy ',/,(i5,e15.6))
c
      return
      end
      subroutine mesh
c
c ---------- get the mesh of the terrain following or flat grid
c            metrics stored at cell centers
c            should be made more general to avoid grid generation
c            by every processor
c
      use pars
      use fields
#if defined(SWAP)
      use module_byteswap
#endif
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
      logical there
      real fs(nnx,iys:iye), fr(nnx,iys:iye)
      integer istatus(mpi_status_size)
c
      integer(kind=mpi_offset_kind) :: offset, disp
      integer(kind=k8)                 nsize, nsize2
      parameter(nvar_g = 4)
c
      real ztemp(0:nnz+2), fnt1(nnx,iys:iye,izs-1:ize+1)
      real, allocatable, dimension(:,:,:) :: temp
      allocate(temp(nvar_g,nnx,iys:iye))
c
c ------------ metric definitions
c
c     gm(.,.,.,1)  = z_zeta = 1/J
c     gm(.,.,.,2)  = zeta_x = - z_xi*J
c     gm(.,.,.,3)  = zeta_y = - z_eta*J
c     gm(.,.,.,4)  = zeta_z = 1/z_zeta = J
c     gm(.,.,.,5)  = z_t = grid speed; zeta_t = - z_t*J
c
      dx = xl/float(nnx)
      do i=1,nnx
         xgrid(i) = float(i-1)*dx
      enddo
      dy = yl/float(nny)
      do j=1,nny
         ygrid(j) = float(j-1)*dy
      enddo
c
      if(iflat .eq. 0) then
         do k=izs-1,ize+1
            do j=iys,iye
            do i=1,nnx
               zgrid_w(i,j,k,1) = z(k)
               zgrid_u(i,j,k)   = zz(k)
               z_zeta           = 1.0
               gm(i,j,k,1)      = z_zeta
               gm(i,j,k,4)      = 1.0/z_zeta
               gm(i,j,k,2)      = 0.0
               gm(i,j,k,3)      = 0.0
               gm(i,j,k,5)      = 0.0
            enddo
            enddo
         enddo
         do j=iys,iye
         do i=1,nnx
            z1b(i,j)      = zz(1)
            bndy_t(i,j,1) = zz(0)
            hdot_t(i,j,1) = 0.0
         enddo
         enddo
c
         if(l_debug) then
            ix = nnx/2
            iy = iys
            write(nprt,1020) (k,xgrid(ix),zgrid_u(ix,iy,k),
     +                  zgrid_w(ix,iy,k,1), k=izs-1,ize+1)
         endif
c
      else
c
c --------- the z mesh
c
c        if(imesh .eq. 0) then
c          t_bndy = 0.0
c        else
c          t_bndy = amax1(0.,time - t_zero)
c        endif
c
         t_bndy = amax1(0.,time - t_zero)
c
         write(nprt,2001) t_bndy, time, t_zero
 2001    format(' 2001: mesh ',/,
     +          ' t_bndy = ',e15.6,' time = ',e15.6,
     +          ' t_zero = ',e15.6)
c
         if(i3d_surf .eq. 0) then
            call l_bndy(t_bndy)
         elseif(i3d_surf .eq. 1) then
            call get_surf(t_bndy)
         elseif(i3d_surf .eq. 2) then
            call get_surf_gabls(t_bndy)
         elseif(i3d_surf .eq. 3) then
            call get_surf_tank(t_bndy)
         endif
c
c -------- initial guess for stretching factor
c
         z_gues = 1.1
c
         do j=iys,iye
         do i=1,nnx
            if(iz_space == 1) then
               call vgrid(zw1,zi,zl-bndy(i,j),z_gues,nnz,
     +                    ztemp(0),l_root,l_debug)
            elseif(iz_space == 2) then
               z_first = zw1*(1.0 - bndy(i,j)/zl)
               call vgrid(z_first,zi,zl-bndy(i,j),z_gues,nnz,
     +                    ztemp(0),l_root,l_debug)
            elseif(iz_space == 3) then
               call exp_grid(bndy(i,j),zl,z(0),ztemp(0),nnz,a_exp_grd)
            elseif(iz_space == 4) then
               call exp_grid2(bndy(i,j),zl,z(0),ztemp(0),nnz,a_exp_grd2)
            elseif(iz_space == 5) then
               call algeb_grid(bndy(i,j),zl,z(0),ztemp(0),nnz,a_decay)
            else
               call cgrid(zl-bndy(i,j),nnz,ztemp(0))
            endif
            do k=izs-1,ize+1
               zgrid_w(i,j,k,1) = ztemp(k) + bndy(i,j)
            enddo
            do k=izs,ize+1
               zgrid_u(i,j,k) = bndy(i,j) + 0.5*(ztemp(k) + ztemp(k-1))
            enddo
            if(izs == 1) then
               zgrid_u(i,j,izs-1) = 2.0*zgrid_w(i,j,izs-1,1) -
     +                              zgrid_u(i,j,izs)
            else
               zgrid_u(i,j,izs-1) = bndy(i,j) +
     +                              0.5*(ztemp(izs-1) + ztemp(izs-2))
            endif
c
c --------- distance off the surface (every processor)
c
            z1b(i,j)      = 0.5*(ztemp(1) + ztemp(0))
            bndy_t(i,j,1) = bndy(i,j)
            hdot_t(i,j,1) = hdot(i,j)
         enddo
         enddo
c
         if(l_debug) then
            ix = nnx/2
            iy = iys
            write(nprt,1020) (k,xgrid(ix),zgrid_u(ix,iy,k),
     +                     zgrid_w(ix,iy,k,1),
     +                     k=izs-1,ize+1)
         endif
c
c ------ generate metric components at cell centers
c        set metric components outside domain = to boundary values
c
        do k=izs,ize+1
           do j=iys,iye
           do i=1,nnx
              gm(i,j,k,1) = (zgrid_w(i,j,k,1) - zgrid_w(i,j,k-1,1))*
     +                       dzw_i(k)
           enddo
           enddo
        enddo
c
c ---------- pass the vertical gradient at ize
c            to the processor above at level izs-1
c
      nb = myid - ncpu_s
      nt = myid + ncpu_s
      if(iss .eq. 0) then
         nb = mpi_proc_null
      endif
      if(ise .eq. numprocs-1) then
         nt = mpi_proc_null
      endif
      nsend = nnx*(iye + 1 - iys)
      nrecv = nsend
      do j=iys,iye
      do i=1,nnx
         fs(i,j) = gm(i,j,ize,1)
      enddo
      enddo
      call mpi_sendrecv(
     +     fs(1,iys),nsend,mpi_real8,nt,2,
     +     fr(1,iys),nrecv,mpi_real8,nb,2,
     +     mpi_comm_world,istatus,ierr)
      if(iss .ne. 0) then
         do j=iys,iye
         do i=1,nnx
            gm(i,j,izs-1,1) = fr(i,j)
         enddo
         enddo
      else
         do j=iys,iye
         do i=1,nnx
            gm(i,j,izs-1,1) = gm(i,j,izs,1)
         enddo
         enddo
      endif
c
        do k=izs-1,ize+1
           do j=iys,iye
           do i=1,nnx
              gm(i,j,k,4) = 1.0/gm(i,j,k,1)
              fnt1(i,j,k) = zgrid_u(i,j,k)
              gm(i,j,k,5) = 0.0
           enddo
           enddo
           call xderivp(fnt1(1,iys,k),trigx(1,1),xk(1),nnx,iys,iye)
           do j=iys,iye
           do i=1,nnx
              gm(i,j,k,2) = -fnt1(i,j,k)*gm(i,j,k,4)
           enddo
           enddo
        enddo
        if(ize .eq. nnz) then
         k = ize + 1
         j = (iys + iye )/2
         write(nprt,3131) (i,zgrid_u(i,j,k),gm(i,j,k,4),
     +       gm(i,j,k,2),i=1,nnx)
 3131 format('3131 mesh',/,
     +       ' i ',5x,' zu ',5x,' gm 4 ',5x,' gm 2',/,
     +      (i5,3e15.6))
        endif
c
        do k=izs-1,ize+1
           do j=iys,iye
           do i=1,nnx
              fnt1(i,j,k) = zgrid_u(i,j,k)
           enddo
           enddo
        enddo
        call yd_mpi(fnt1(1,iys,izs-1),trigx(1,2),yk(1),
     +              nnx,nny,ixs,ixe,ix_s,ix_e,
     +              iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
        do k=izs-1,ize+1
           do j=iys,iye
           do i=1,nnx
              gm(i,j,k,3) = -fnt1(i,j,k)*gm(i,j,k,4)
           enddo
           enddo
        enddo
c
      endif
      if(l_debug) then
         ix = nnx/2
         iy = iys
         write(nprt,1010) (k,gm(ix,iy,k,1),gm(ix,iy,k,2),
     +                  gm(ix,iy,k,3),gm(ix,iy,k,4),
     +                  gm(ix,iy,k,5),k=izs-1,ize+1)
 1010    format(' ind ',5x,' 1/J',5x,' -z_xi*J',5x,' -z_eta*J',5x,
     +          ' J',5x,' -z_t*J',/,
     +          (i5,5e15.6))
 1020    format(' k',5x,' xgrid',5x,' zgrid_u',5x,'zgrid_w',/,
     +          (i5,3e15.6))
         k  = izs-1
         iy = iys
         write(nprt,1010) (ix,gm(ix,iy,k,1),gm(ix,iy,k,2),
     +                  gm(ix,iy,k,3),gm(ix,iy,k,4),
     +                  gm(ix,iy,k,5),ix=1,nnx)
         k  = izs-1
         ix = nnx/2
         write(nprt,1010) (iy,gm(ix,iy,k,1),gm(ix,iy,k,2),
     +                  gm(ix,iy,k,3),gm(ix,iy,k,4),
     +                  gm(ix,iy,k,5),iy=iys,iye)
      endif
c
c ----------- save the grid file
c
      path_sav_g = trim(path_sav)//'/'//case(1:3)//'.grid'
c     write(nprt,9001) path_sav_g
c9001 format(' 9001: path = ',a80)
c
      call mpi_file_open(mpi_comm_world, path_sav_g,
     +                   mpi_mode_create+mpi_mode_rdwr,
     +                   mpi_info_null, ngrd, ierr)
c
c ---- set file view
c
      disp = 0
      call mpi_file_set_view(ngrd,disp,mpi_real8,mpi_real8,
     +                      'native',mpi_info_null,ierr)
c
c ---- write data
c
      nsize   = int(nvar_g,k8)*nnx*nny
      nsize2  = int(nvar_g,k8)*nnx*(iys-1)
      n_write = nvar_g*nnx*(iye+1-iys)
c
      if(ize .eq. nnz) then
         izen = ize + 1
      else
         izen = ize
      endif
c     do k=izs,ize
      do k=izs,izen
         do j=iys,iye
         do i=1,nnx
            temp(1,i,j) = zgrid_w(i,j,k-1,1)
            temp(2,i,j) = gm(i,j,k-1,1)
            temp(3,i,j) = gm(i,j,k-1,2)
            temp(4,i,j) = gm(i,j,k-1,3)
         enddo
         enddo
#if defined(SWAP)
      call byteswap(temp)
#endif
c
         offset = int((k-1),k8)*nsize + nsize2
c
c -------- mpi i/o with non-uniform blocks of data use _at
c
         call mpi_file_write_at(ngrd,offset,temp,n_write,
     +                              mpi_real8,status,ierr)
         if (ierr /= 0) goto 9991

      enddo
c
c ---- close file
c
      call mpi_file_close(ngrd, ierr)
c
c ---- check file
c
      if (l_root) then
         inquire(file=path_sav_g,exist=there)
         if(.not.there) then
            write(6,8000) ngrd,myid
            call mpi_abort(mpi_comm_world,icode,ierr)
         endif
         write(6,7000) path_sav_g
      endif
c
      deallocate(temp)
c
c ----------- store metrics for possible time variation
c
      do l=1,5
         do k=izs-1,ize+1
         do j=iys,iye
         do i=1,nnx
            gm2(i,j,k,l) = gm(i,j,k,l)
         enddo
         enddo
         enddo
      enddo
c
      do k=izs-1,ize+1
      do j=iys,iye
      do i=1,nnx
         zgrid_w(i,j,k,2) = zgrid_w(i,j,k,1)
      enddo
      enddo
      enddo
c
      do j=iys,iye
      do i=1,nnx
         z1b_2(i,j)    = z1b(i,j)
         bndy_t(i,j,2) = bndy_t(i,j,1)
         hdot_t(i,j,2) = hdot_t(i,j,1)
      enddo
      enddo
c
      return
c --------------------------  errors in writing grid file
 9991 continue
      write(6,6000) ngrd, k
 6000 format(' SR. MESH:',/,
     +       '    trouble cannot write grid file on unit = ',i2,/,
     +       '             at iz = ',i4)
      call mpi_abort(mpi_comm_world,icode,ierr)
c --------------------
 7000 format(' **** GRID IS WRITTEN IN FILE  ',a80)
 8000 format(' in MESH: trouble writing file ',i5,'  myid = ',i5)
      end
      subroutine new_mesh(it,istage)
c
c ---------- get a new mesh at time step t_stage3
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
      logical there
      real fs(nnx,iys:iye), fr(nnx,iys:iye)
      integer istatus(mpi_status_size)
c
      real ztemp(0:nnz+2), fnt1(nnx,iys:iye,izs-1:ize+1)
c
c ------------ metric definitions
c
c     gm2(.,.,.,1)  = z_zeta = 1/J
c     gm2(.,.,.,2)  = zeta_x = - z_xi*J
c     gm2(.,.,.,3)  = zeta_y = - z_eta*J
c     gm2(.,.,.,4)  = zeta_z = 1/z_zeta = J
c     gm2(.,.,.,5)  = z_t = grid speed; zeta_t = - z_t*J
c
c --------- update the metrics [note gm2(.,.,.,5) = 0]
c
      do l=1,5
         do k=izs-1,ize+1
         do j=iys,iye
         do i=1,nnx
            gm(i,j,k,l)  = gm2(i,j,k,l)
            gm2(i,j,k,l) = gm3(i,j,k,l)
         enddo
         enddo
         enddo
      enddo
c
      do k=izs-1,ize+1
      do j=iys,iye
      do i=1,nnx
         zgrid_w(i,j,k,1) = zgrid_w(i,j,k,2)
         zgrid_w(i,j,k,2) = zgrid_w(i,j,k,3)
      enddo
      enddo
      enddo
c
c --------- distance off the surface
c
         do j=iys,iye
         do i=1,nnx
            z1b(i,j)      = z1b_2(i,j)
            z1b_2(i,j)    = z1b_3(i,j)
            bndy_t(i,j,1) = bndy_t(i,j,2)
            bndy_t(i,j,2) = bndy_t(i,j,3)
            hdot_t(i,j,1) = hdot_t(i,j,2)
            hdot_t(i,j,2) = hdot_t(i,j,3)
         enddo
         enddo
c
c --------- get the new boundary shape
c           at t = t_stage3
c
         t_bndy = t_stage3 - t_zero
c
c        write(nprt, 1010) t_bndy, t_stage3, t_zero
c1010    format(' 1010: new mesh',/,
c    +          ' t_bndy = ',e15.6,
c    +          ' t_stage3 = ',e15.6,' t_zero = ',e15.6)
c
         if(i3d_surf .eq. 0) then
            call l_bndy(t_bndy)
         elseif(i3d_surf .eq. 1) then
            call get_surf(t_bndy)
         elseif(i3d_surf .eq. 2) then
            call get_surf_gabls(t_bndy)
         elseif(i3d_surf .eq. 3) then
            call get_surf_tank(t_bndy)
         endif
c
c -------- initial guess for stretching factor
c
         z_gues = 1.1
c
         do j=iys,iye
         do i=1,nnx
            if(iz_space == 1) then
               call vgrid(zw1,zi,zl-bndy(i,j),z_gues,nnz,
     +                    ztemp(0),l_root,l_debug)
            elseif(iz_space == 2) then
               z_first = zw1*(1.0 - bndy(i,j)/zl)
               call vgrid(z_first,zi,zl-bndy(i,j),z_gues,nnz,
     +                    ztemp(0),l_root,l_debug)
            elseif(iz_space == 3) then
               call exp_grid(bndy(i,j),zl,z(0),ztemp(0),nnz,a_exp_grd)
            elseif(iz_space == 4) then
               call exp_grid2(bndy(i,j),zl,z(0),ztemp(0),nnz,a_exp_grd2)
            elseif(iz_space == 5) then
               call algeb_grid(bndy(i,j),zl,z(0),ztemp(0),nnz,a_decay)
            else
               call cgrid(zl-bndy(i,j),nnz,ztemp(0))
            endif
            do k=izs-1,ize+1
               zgrid_w(i,j,k,3) = ztemp(k) + bndy(i,j)
            enddo
            do k=izs,ize+1
               zgrid_u(i,j,k) = bndy(i,j) + 0.5*(ztemp(k) + ztemp(k-1))
            enddo
            if(izs == 1) then
               zgrid_u(i,j,izs-1) = 2.0*zgrid_w(i,j,izs-1,3) -
     +                              zgrid_u(i,j,izs)
            else
               zgrid_u(i,j,izs-1) = bndy(i,j) +
     +                              0.5*(ztemp(izs-1) + ztemp(izs-2))
            endif
c
c --------- distance off the surface (every processor)
c
            z1b_3(i,j)    = 0.5*(ztemp(1) + ztemp(0))
            bndy_t(i,j,3) = bndy(i,j)
            hdot_t(i,j,3) = hdot(i,j)
         enddo
         enddo
c
c        j = iys
c        k = 1
c        write(nprt,6501) (i, bndy(i,j), zgrid_w(i,j,k,3),
c    +                     z1b(i,j),z1b_2(i,j),
c    +                     z1b_3(i,j),i=1,nnx)
c6501    format(' 6501 new mesh:',/,
c    + ' i ',5x,' bndy ',5x,' zw(1)',5x,' z1b',5x,
c    + ' z2b ',5x,' z3b',/,(i5,5e15.6))
c
c ------ generate metric components at cell centers
c        set metric components outside domain = to boundary values
c
        do k=izs,ize+1
           do j=iys,iye
           do i=1,nnx
              gm3(i,j,k,1) = (zgrid_w(i,j,k,3) - zgrid_w(i,j,k-1,3))*
     +                       dzw_i(k)
           enddo
           enddo
        enddo
c
c ---------- pass the 1/J at ize
c            to the processor above at level izs-1
c
      nb = myid - ncpu_s
      nt = myid + ncpu_s
      if(iss .eq. 0) then
         nb = mpi_proc_null
      endif
      if(ise .eq. numprocs-1) then
         nt = mpi_proc_null
      endif
      nsend = nnx*(iye + 1 - iys)
      nrecv = nsend
      do j=iys,iye
      do i=1,nnx
         fs(i,j) = gm3(i,j,ize,1)
      enddo
      enddo
      call mpi_sendrecv(
     +     fs(1,iys),nsend,mpi_real8,nt,2,
     +     fr(1,iys),nrecv,mpi_real8,nb,2,
     +     mpi_comm_world,istatus,ierr)
      if(iss .ne. 0) then
         do j=iys,iye
         do i=1,nnx
            gm3(i,j,izs-1,1) = fr(i,j)
         enddo
         enddo
      else
         do j=iys,iye
         do i=1,nnx
            gm3(i,j,izs-1,1) = gm3(i,j,izs,1)
         enddo
         enddo
      endif
c
c ----- compute grid speed gm2(.,.,.,5) in sr. grid_speed
c
      do k=izs-1,ize+1
         do j=iys,iye
         do i=1,nnx
            gm3(i,j,k,4) = 1.0/gm3(i,j,k,1)
            fnt1(i,j,k) = zgrid_u(i,j,k)
         enddo
         enddo
         call xderivp(fnt1(1,iys,k),trigx(1,1),xk(1),nnx,iys,iye)
         do j=iys,iye
         do i=1,nnx
            gm3(i,j,k,2) = -fnt1(i,j,k)*gm3(i,j,k,4)
         enddo
         enddo
      enddo
c
      do k=izs-1,ize+1
         do j=iys,iye
         do i=1,nnx
            fnt1(i,j,k) = zgrid_u(i,j,k)
         enddo
         enddo
      enddo
      call yd_mpi(fnt1(1,iys,izs-1),trigx(1,2),yk(1),
     +            nnx,nny,ixs,ixe,ix_s,ix_e,
     +            iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
      do k=izs-1,ize+1
         do j=iys,iye
         do i=1,nnx
            gm3(i,j,k,3) = -fnt1(i,j,k)*gm3(i,j,k,4)
         enddo
         enddo
      enddo
c
      return
      end
      subroutine start_mesh(it,istage)
c
c ---------- case with moving mesh
c            initialize results at time step t_stage2
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
      logical there
      real fs(nnx,iys:iye), fr(nnx,iys:iye)
      integer istatus(mpi_status_size)
c
      real ztemp(0:nnz+2), fnt1(nnx,iys:iye,izs-1:ize+1)
c
c ------------ metric definitions
c
c     gm2(.,.,.,1)  = z_zeta = 1/J
c     gm2(.,.,.,2)  = zeta_x = - z_xi*J
c     gm2(.,.,.,3)  = zeta_y = - z_eta*J
c     gm2(.,.,.,4)  = zeta_z = 1/z_zeta = J
c     gm2(.,.,.,5)  = z_t = grid speed; zeta_t = - z_t*J
c
c --------- update the metrics [note gm2(.,.,.,5) = 0]
c
      do l=1,5
         do k=izs-1,ize+1
         do j=iys,iye
         do i=1,nnx
            gm2(i,j,k,l)  = gm(i,j,k,l)
         enddo
         enddo
         enddo
      enddo
c
      do j=iys,iye
      do i=1,nnx
         z1b_2(i,j)    = z1b(i,j)
         bndy_t(i,j,2) = bndy_t(i,j,1)
         hdot_t(i,j,2) = hdot_t(i,j,1)
      enddo
      enddo
c
c --------- get the new boundary shape
c           at t = t_stage2
c
         t_bndy = t_stage2 - t_zero
c
         if(i3d_surf .eq. 0) then
            call l_bndy(t_bndy)
         elseif(i3d_surf .eq. 1) then
            call get_surf(t_bndy)
         elseif(i3d_surf .eq. 2) then
            call get_surf_gabls(t_bndy)
         elseif(i3d_surf .eq. 3) then
            call get_surf_tank(t_bndy)
         endif
c
c -------- initial guess for stretching factor
c
         z_gues = 1.1
c
         do j=iys,iye
         do i=1,nnx
            if(iz_space == 1) then
               call vgrid(zw1,zi,zl-bndy(i,j),z_gues,nnz,
     +                    ztemp(0),l_root,l_debug)
            elseif(iz_space == 2) then
               z_first = zw1*(1.0 - bndy(i,j)/zl)
               call vgrid(z_first,zi,zl-bndy(i,j),z_gues,nnz,
     +                    ztemp(0),l_root,l_debug)
            elseif(iz_space == 3) then
               call exp_grid(bndy(i,j),zl,z(0),ztemp(0),nnz,a_exp_grd)
            elseif(iz_space == 4) then
               call exp_grid2(bndy(i,j),zl,z(0),ztemp(0),nnz,a_exp_grd2)
            elseif(iz_space == 5) then
               call algeb_grid(bndy(i,j),zl,z(0),ztemp(0),nnz,a_decay)
            else
               call cgrid(zl-bndy(i,j),nnz,ztemp(0))
            endif
            do k=izs-1,ize+1
               zgrid_w(i,j,k,3) = ztemp(k) + bndy(i,j)
            enddo
            do k=izs,ize+1
               zgrid_u(i,j,k) = bndy(i,j) + 0.5*(ztemp(k) + ztemp(k-1))
            enddo
            if(izs == 1) then
               zgrid_u(i,j,izs-1) = 2.0*zgrid_w(i,j,izs-1,3) -
     +                              zgrid_u(i,j,izs)
            else
               zgrid_u(i,j,izs-1) = bndy(i,j) +
     +                              0.5*(ztemp(izs-1) + ztemp(izs-2))
            endif
c
c --------- distance off the surface (every processor)
c
            z1b_3(i,j)    =  0.5*(ztemp(1) + ztemp(0))
            bndy_t(i,j,3) =  bndy(i,j)
            hdot_t(i,j,3) =  hdot(i,j)
         enddo
         enddo
c
c ------ generate metric components at cell centers
c        set metric components outside domain = to boundary values
c
        do k=izs,ize+1
           do j=iys,iye
           do i=1,nnx
              gm3(i,j,k,1) = (zgrid_w(i,j,k,3) - zgrid_w(i,j,k-1,3))*
     +                       dzw_i(k)
           enddo
           enddo
        enddo
c
c ---------- pass the 1/J at ize
c            to the processor above at level izs-1
c
      nb = myid - ncpu_s
      nt = myid + ncpu_s
      if(iss .eq. 0) then
         nb = mpi_proc_null
      endif
      if(ise .eq. numprocs-1) then
         nt = mpi_proc_null
      endif
      nsend = nnx*(iye + 1 - iys)
      nrecv = nsend
      do j=iys,iye
      do i=1,nnx
         fs(i,j) = gm3(i,j,ize,1)
      enddo
      enddo
      call mpi_sendrecv(
     +     fs(1,iys),nsend,mpi_real8,nt,2,
     +     fr(1,iys),nrecv,mpi_real8,nb,2,
     +     mpi_comm_world,istatus,ierr)
      if(iss .ne. 0) then
         do j=iys,iye
         do i=1,nnx
            gm3(i,j,izs-1,1) = fr(i,j)
         enddo
         enddo
      else
         do j=iys,iye
         do i=1,nnx
            gm3(i,j,izs-1,1) = gm3(i,j,izs,1)
         enddo
         enddo
      endif
c
c ----- compute grid speed gm2(.,.,.,5) in sr. grid_speed
c
      do k=izs-1,ize+1
         do j=iys,iye
         do i=1,nnx
            gm3(i,j,k,4) = 1.0/gm3(i,j,k,1)
            fnt1(i,j,k) = zgrid_u(i,j,k)
         enddo
         enddo
         call xderivp(fnt1(1,iys,k),trigx(1,1),xk(1),nnx,iys,iye)
         do j=iys,iye
         do i=1,nnx
            gm3(i,j,k,2) = -fnt1(i,j,k)*gm3(i,j,k,4)
         enddo
         enddo
      enddo
c
      do k=izs-1,ize+1
         do j=iys,iye
         do i=1,nnx
            fnt1(i,j,k) = zgrid_u(i,j,k)
         enddo
         enddo
      enddo
      call yd_mpi(fnt1(1,iys,izs-1),trigx(1,2),yk(1),
     +            nnx,nny,ixs,ixe,ix_s,ix_e,
     +            iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
      do k=izs-1,ize+1
         do j=iys,iye
         do i=1,nnx
            gm3(i,j,k,3) = -fnt1(i,j,k)*gm3(i,j,k,4)
         enddo
         enddo
      enddo
c
      return
      end
      subroutine grid_speed(it,istage)
c
c ---------- get grid speeds for the mesh at t = t_stage
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real fs(nnx,iys:iye), fr(nnx,iys:iye), rh(nnz)
      real gt(0:nnz+1,iys:iye,mxs:mxe)
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      if(istage .le. 2) then
         dtgama_2 = dt*gama(istage+1)
         dtzeta_2 = dt*zetas(istage+1)
      else
         dtgama_2 = dt1*gama(istage-2)
         dtzeta_2 = dt1*zetas(istage-2)
      endif
c
c ------- get the grid speed at t = t_stage
c
      gami = 1.0/dtgama_2
c
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            gm2(i,j,k,5) = dzw(k)*gami*(gm3(i,j,k,1) - gm2(i,j,k,1) -
     +                       dtzeta_2*rj(i,j,k))
         enddo
         enddo
      enddo
      call xtoz_trans(gm2(1,iys,izs-1,5),gt,nnx,nnz,mxs,mxe,mx_s,mx_e,
     +                iys,iye,izs,ize,iz_s,iz_e,myid,ncpu_s,
     +                numprocs)
c
c ---------- integrate from bottom up to find grid speeds
c            using surface value as boundary condition
c            hdot comes from boundary routine
c
      do i=mxs,mxe
      do j=iys,iye
         do k=1,nnz
            rh(k) = gt(k,j,i)
         enddo
         gt(0,j,i) = hdot_t(i,j,2)
         do k=1,nnz
            gt(k,j,i) = rh(k) + gt(k-1,j,i)
         enddo
         gt(nnz+1,j,i) = gt(nnz,j,i)
      enddo
      enddo
c
c --------- transpose back defining ghost point values [izs-1,ize+1]
c           surface value [izs-1] /= value at [izs]
c           with non-uniform spacing of first gridline.
c
      call ztox_trans(gt,gm2(1,iys,izs-1,5),nnx,nnz,mxs,mxe,mx_s,mx_e,
     +                iys,iye,izs,ize,iz_s,iz_e,myid,ncpu_s,
     +                numprocs)
c
c ------- get the new rhs of the 1/J equation
c
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            rj(i,j,k) = (gm2(i,j,k,5) - gm2(i,j,k-1,5))*dzw_i(k)
         enddo
         enddo
      enddo
c
c ---------- determine boundary speeds at t = t_stage
c            for use in psolver
c
      if(iss .eq. 0) then
         k = 0
         do j=iys,iye
         do i=1,nnx
            wfbc(i,j,3) = gm2(i,j,k,5)
         enddo
         enddo
      endif
c
      if(ise .eq. numprocs-1) then
         k = nnz
         do j=iys,iye
         do i=1,nnx
            wfbc(i,j,4) = gm2(i,j,k,5)
         enddo
         enddo
      endif
c
      return
      end
      subroutine algeb_grid(h,zl,zeta,zw,nz,alp)
c
c ----------- algebraic mapping with grid decay.
c             nicely maps boundary shape as function of(x,y)
c             on to computational coordinates. Use a_decay (alp)
c             to control how fast bumps decay with height.
c             Large values cause non-uniform spacing of
c             first grid point. see hres test code
c
      real zeta(0:nz+1), zw(0:nz+1)
c
      do i=0,nz
         zw(i) = zeta(i) + h*(1.0 - zeta(i)/zl)**alp
      enddo
      zw(nz+1) = zw(nz) + (zw(nz) - zw(nz-1))
c
c --------- fix to match rest of code
c
      do i=0,nz+1
         zw(i) = zw(i) - h
      enddo
c
      return
      end
      subroutine exp_grid(h,zl,zeta,zw,nz,alp)
c
c ----------- simple decaying exponential grid
c
      real zeta(0:nz+1), zw(0:nz+1)
c
c --------- choose alp to match physical problem
c
c     alp = 0.017
c     alp = 0.020
c
c ----------- tank simulations
c
c     alp = 66.67
c     alp = 50.00
c
      do i=0,nz
         zw(i) = zeta(i) + h*exp(-alp*zeta(i))
      enddo
c     zw(nz)   = zeta(nz)
      zw(nz+1) = zw(nz) + (zw(nz) - zw(nz-1))
c
c --------- fix to match rest of code
c
      do i=0,nz+1
         zw(i) = zw(i) - h
      enddo
c
      return
      end
      subroutine exp_grid2(h,zl,zeta,zw,nz,alp)
c
c ----------- simple decaying exponential grid
c             but with first grid z(1) = zeta(1)
c             works best to pick box height
c             zl(code) = zl(desired) + zeta(1)
c             otherwise top of box in physical space
c             will be shorter by first grid point off the
c             surface
c
      real zeta(0:nz+1), zw(0:nz+1)
c
c     alp = 0.020
c     alp = 0.030
c     alp = 0.017
c     alp = 0.00005
c     alp = 0.00010
c
      zw(0) = h
      do i=1,nz
         d_zeta = zeta(i) - zeta(1)
         zw(i) = d_zeta + (h + zeta(1))*exp(-alp*d_zeta**2)
      enddo
c     zw(nz)   = zeta(nz)
      zw(nz+1) = zw(nz) + (zw(nz) - zw(nz-1))
c
c --------- fix to match rest of code
c
      do i=0,nz+1
         zw(i) = zw(i) - h
      enddo
c
      return
      end
      subroutine vgrid(z1,zi,zl,z_gues,nnz,z,l_root,ldebug)
c
c ------ iterate sufficiently so that z(nnz) = zl (use r8)
c
      real z(0:nnz+1)
      logical l_root, l_debug
      data ig_max /500/
      save ig_max
c
c ----------------- build grid up to zi first
c
      z_frst = z1
      z_cntr = zl
      n_pbl  = nnz
      z_fac1 = z_cntr/z_frst
      z_fac2 = 1.0/float(nnz)
      z_fac  = z_gues
      knt    = 0
      tol    = 1.0e-10
   10 continue
        knt = knt + 1
        z_facn = (z_fac1*(z_fac - 1.0) + 1.0)**z_fac2
        test   = abs(1.0 - z_facn/z_fac)
        if(knt .gt. ig_max) then
            if(l_root) write(6,9000) z_fac, z_facn, knt
 9000       format(' Cannot find stretching factor',/,
     +             ' z_fac = ',e15.6,' z_facn = ',e15.6,' knt = ',i3)
            stop
        endif
        z_fac = z_facn
        if(test .gt. tol) go to 10
c
c     if(l_root) write(6,9100) z_fac, z_cntr, z1, knt
c9100 format(' Stretching factor = ',e15.6,/,
c    +       ' Match point       = ',e15.6,/,
c    +       ' First z           = ',e15.6,/,
c    +       ' Number of iters   = ',i4)
      z(1) = z_frst
      do iz=2,nnz
         z(iz) = z_frst*(z_fac**(float(iz)) - 1.0)/(z_fac - 1.0)
      enddo
      z(nnz)   = zl
      z(0)     = 0.0
      z(nnz+1) = z(nnz) + (z(nnz) - z(nnz-1))
c
      return
      end
      subroutine get_dz
c
c --------------- compute spacing for given vertical
c                 point distribution
c
      use pars
      use fields
      use con_data
      use con_stats
      include 'mpif.h'
c
      do iz=1,nnz+1
         dzw(iz) = z(iz) - z(iz-1)
      enddo
      dzw(0)     = dzw(1)
      dzw(nnz+2) = dzw(nnz+1)
      do iz=0,nnz+2
         dzw_i(iz) = 1.0/dzw(iz)
      enddo
c
c ------------ build weights for vertical interpolation
c              of a centered variable to a wf-level
c
      do iz=0,nnz+1
         weit(iz)  = dzw(iz)/(dzw(iz) + dzw(iz+1))
         weitm(iz) = 1.0 - weit(iz)
      enddo
c
c ------------ build z grid for u points
c
      dzovr2 = dz*0.5
      do iz=1,nnz+1
         zz(iz) = 0.5*(z(iz) + z(iz-1))
      enddo
      zz(0) = - zz(1)
      do iz=1,nnz+1
         dzu(iz) = zz(iz) - zz(iz-1)
      enddo
      dzu(0)     = dzu(1)
      dzu(nnz+2) = dzu(nnz+1)
      do iz=0,nnz+2
         dzu_i(iz) = 1.0/dzu(iz)
      enddo
c
      return
      end
      subroutine vgrid_channel(z1,zi,zl,nnz,z,l_root,ldebug)
c
      real z(0:nnz+1)
      integer :: zidx
      logical l_root, l_debug
c
c ----------------- build grid up to zi first
c
      z_frst = z1
      z_cntr = zi*0.5
      n_pbl  = nnz
c     n_pbl  = (5*nnz)/8
      z_fac1 = z_cntr/z_frst
      z_fac2 = 1.0/float(n_pbl/2)
      z_fac  = 1.1
      knt = 0
      tol = 0.00001
   10 continue
        knt = knt + 1
        z_facn = (z_fac1*(z_fac - 1.0) + 1.0)**z_fac2
        test   = abs(1.0 - z_facn/z_fac)
        if(knt .gt. 50) then
            if(l_root) write(6,9000) z_fac, z_facn, knt
 9000       format(' Cannot find stretching factor',/,
     +             ' z_fac = ',e15.6,' z_facn = ',e15.6,' knt = ',i3)
            stop
        endif
        z_fac = z_facn
        if(test .gt. tol) go to 10
      if(l_root) write(6,9100) z_fac, z_cntr, z1, knt
 9100 format(' Stretching factor = ',e15.6,/,
     +       ' Match point       = ',e15.6,/,
     +       ' First z           = ',e15.6,/,
     +       ' Number of iters   = ',i4)
      z(1) = z_frst
      do iz=2,n_pbl/2-1
         z(iz) = z_frst*(z_fac**(float(iz)) - 1.0)/(z_fac - 1.0)
      enddo
      z(n_pbl/2) = z_cntr
      do iz=1,n_pbl/2 - 1
         z(n_pbl-iz) = zi - z(iz)
      enddo
      z(n_pbl) = zi
      z(0)   = 0.0
c
      if(l_root) write(6,5300) n_pbl
 5300 format(' n_pbl = ',i4)
c
c -------------- build grid from zi on up
!     For the channel, zi represents the channel centerline
!     Want the mesh to be a mirror image across this:
c
!      zidx = 1
!      do iz=n_pbl+1,nnz
!         z(iz) = zi + (zi - z(n_pbl-zidx))
!         zidx = zidx + 1
!      enddo
      z(nnz+1) = z(nnz) + (z(nnz) - z(nnz-1))
c
      return
      end subroutine vgrid_channel
      subroutine random
c
c ----------- geostrophic winds designed for comparison case
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use particles
      real psi(nnx,iys:iye), psix(nnx,iys:iye),
     +     psiy(nnx,iys:iye,izs:izs), uxx(nnx,iys:iye),
     +     vyy(nnx,iys:iye,izs:izs)

      !Initialize partcount to 0:
      partcount = 0.0
      partcount_t = 0.0
      partcountU = 0.0
      partcountU_t = 0.0
      partflux = 0.0
      partflux_t = 0.0
      partfluxup = 0.0
      partfluxup_t = 0.0
      partfluxdn = 0.0
      partfluxdn_t = 0.0
      vpsum = 0.0
      vpsum_t = 0.0
      vpsqrsum = 0.0
      vpsqrsum_t = 0.0
      upwp_t = 0.0
      upwp = 0.0
      Tpsum = 0.0
      Tpsum_t = 0.0
      Tpsqrsum = 0.0
      Tpsqrsum_t = 0.0
      wpTpsum = 0.0
      wpTpsum_t = 0.0
      partsrc = 0.0
      partsrc_t = 0.0
      partTsrc = 0.0
      partTsrc_t = 0.0
      qfsum = 0.0
      qfsum_t = 0.0
      radsum = 0.0
      radsum_t = 0.0
      rad2sum = 0.0
      rad2sum_t = 0.0
      qstarsum = 0.0
      qstarsum_t = 0.0


c
c ------------ note set nmatch in sr. iso so that
c              it is compatible with conditions here
c
      do iz=1,nnz
c        ug(iz)   = ugcont*(zz(iz)/zl)
         ug(iz)   = ugcont
         vg(iz)   = vgcont
         divz(iz) = 0.0
      enddo
c
CHECK
c
c     izi = (50*nnz)/100
c
      izi = 92
      zi  = z(izi)
c
c     z_lower = zi - 50.0
c     t_lower = 300.0
c     z_upper = zi + 50.0
c     t_upper = 308.0
c     slope   = (t_upper - t_lower)/(z_upper - z_lower)
c
c -------- alternate strategy for specifying theta
c
      z_lower = 300.0
      t_lower = 290.0
c
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz) = ugcont-ugal
c           u(ix,iy,iz) = ug(iz) - ugal
            v(ix,iy,iz) = vgcont
            w(ix,iy,iz) = 0.0
            e(ix,iy,iz) = 0.0
         enddo
         enddo
         do iy=iys,iye
         do ix=1,nnx
!            if(zgrid_w(ix,iy,iz,1) .le. z_lower) then
!              t(ix,iy,1,iz) = t_lower
!            else
!              t(ix,iy,1,iz) = t_lower +
!     +                        (zgrid_u(ix,iy,iz+1) - z_lower)*dtdzf(1)
!            endif
             t(ix,iy,1,iz) = 300.0
             t(ix,iy,2,iz) = 0.01
         enddo
         enddo
         do iy=iys,iye
         do ix=1,nnx
            w(ix,iy,iz)   = 0.
            r1(ix,iy,iz)  = 0.
            r2(ix,iy,iz)  = 0.
            r3(ix,iy,iz)  = 0.
            r4(ix,iy,1,iz)= 0.
            r5(ix,iy,iz)  = 0.
            p(ix,iy,iz)   = 0.
            rj(ix,iy,iz)  = 0.
c
c ---------- flux velocities simple guess
c
c           uf(ix,iy,iz)  = u(ix,iy,iz)*gm(ix,iy,iz,1)
c           vf(ix,iy,iz)  = v(ix,iy,iz)*gm(ix,iy,iz,1)
c           wf(ix,iy,iz)  = u(ix,iy,iz)*gm(ix,iy,iz,2) +
c    +                      v(ix,iy,iz)*gm(ix,iy,iz,3)
            uf(ix,iy,iz)  = u(ix,iy,iz)
            vf(ix,iy,iz)  = v(ix,iy,iz)
            wf(ix,iy,iz)  = 0.0
         enddo
         enddo
      enddo
c
c ------------- set initial random field to be
c               divergence free
c
      idum = -1 - myid
      do iz=izs,ize
c
c ----------- ampv and ampt are max amplitudes of random
c             velocity and temperature fields
c             make sure ampv is set if free convection so
c             that we have motions at first time step
c
         ampv = 0.0
         ampt = 0.10
c
c ------- simple random field scaled between -0.5 and 0.5
c
         sum_psi = 0.0
         do iy=iys,iye
         do ix=1,nnx
            psi(ix,iy) = ran1(idum)
            sum_psi = sum_psi + psi(ix,iy)
         enddo
         enddo
         sum_psi = sum_psi*fnxy
         call mpi_sum_xy(sum_psi,myid,iss,ise,1)
         do iy=iys,iye
         do ix=1,nnx
            psi(ix,iy) = psi(ix,iy) - sum_psi
            psix(ix,iy)     = psi(ix,iy)
            psiy(ix,iy,izs) = psi(ix,iy)
         enddo
         enddo
         call xderivp(psix(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call yd_mpi(psiy(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
         vmaxx = 0.0
         do iy=iys,iye
         do ix=1,nnx
            vmag = sqrt(psix(ix,iy)**2 + psiy(ix,iy,izs)**2)
            if(vmag .gt. vmaxx) vmaxx = vmag
         enddo
         enddo
         facv = ampv/vmaxx
c
         if (z(iz) .le. 50.0) then
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz)   = u(ix,iy,iz) - psiy(ix,iy,izs)*facv
            v(ix,iy,iz)   = v(ix,iy,iz) + psix(ix,iy)*facv
            t(ix,iy,1,iz) = t(ix,iy,1,iz) + psi(ix,iy)*ampt
         enddo
         enddo
         endif
c
         if(z(iz) .le. 250.0) then
         do iy=iys,iye
         do ix=1,nnx
            e(ix,iy,iz) = 0.4*(1.0 - z(iz)/250.0)**3
         enddo
         enddo
         endif
c
c ---------- check divergence of initial field
c
         do iy=iys,iye
         do ix=1,nnx
            uxx(ix,iy)     = uf(ix,iy,iz)
            vyy(ix,iy,izs) = vf(ix,iy,iz)
         enddo
         enddo
         call xderivp(uxx(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call yd_mpi(vyy(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
         do iy=iys,iye
         do ix=1,nnx
            divz(iz) = divz(iz) + (uxx(ix,iy) + vyy(ix,iy,izs))**2
         enddo
         enddo
         divz(iz) = divz(iz)*fnxy
c
c -------- end z loop
c
      enddo
c
      call mpi_sum_z(divz(1),i_root,myid,nnz,1)
c
      write(nprt,6000)
 6000 format(' check of divergence for initial state',/,
     +         ' iz ',5x,' divergence')
      write(nprt,6100) (iz,divz(iz),iz=izs,ize)
 6100 format(i5,e15.6)
c
c ------------ fix for baroclinic and subsidence effects !!
c
c     do iz=izs,ize
c        ug(iz)=ugcont
c        vg(iz)=vgcont
c        if (.not.(ibrcl.eq.1)) go to 19988
c        if (.not.(iz.le.izi)) go to 19987
c        ug(iz)=0.
c        vg(iz)=0.
c 19987    continue
c 19988    continue
c        zz2=zz(iz)
c        wls(iz)=-divgls*zz2
c        if (.not.(iz.eq.1)) go to 19986
c        do ix=1,nnx
c        uls(ix)=divgls*(dx*float(ix-1)-xl*.5)
c        enddo
c     enddo
c     write(nprt,9)(uls(ix),ix=1,nnx)
c  9  format(1x,8e12.3)
c 19986 continue
c
      return
      end
      subroutine random_wt
c
c ----------- random conditions for a wind tunnel flow
c             made dimensionless by u_* and wavelength
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real psi(nnx,iys:iye), psix(nnx,iys:iye),
     +     psiy(nnx,iys:iye,izs:izs), uxx(nnx,iys:iye),
     +     vyy(nnx,iys:iye,izs:izs)
c
c ------------ note set nmatch in sr. iso so that
c              it is compatible with conditions here
c
      do iz=1,nnz
c        ug(iz)   = ugcont*(zz(iz)/zl)
         ug(iz)   = ugcont
         vg(iz)   = vgcont
         divz(iz) = 0.0
      enddo
c
CHECK
c
c     izi = (50*nnz)/100
c
      izi = 100
      zi  = z(izi)
c
      ustar_gues = 1.0
      ustar_vk   = ustar_gues/vk
      do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
!            u(ix,iy,iz) = ustar_vk*
!     +                    alog((zgrid_u(ix,iy,iz) - bndy(ix,iy))/zo)
            u(ix,iy,iz) = ugcont
            v(ix,iy,iz) = 0.0
            w(ix,iy,iz) = 0.0
            e(ix,iy,iz) = 0.0
         enddo
         enddo
         do iy=iys,iye
         do ix=1,nnx
            t(ix,iy,1,iz) = tsfcc(1)
         enddo
         enddo
         do iy=iys,iye
         do ix=1,nnx
            w(ix,iy,iz)   = 0.
            r1(ix,iy,iz)  = 0.
            r2(ix,iy,iz)  = 0.
            r3(ix,iy,iz)  = 0.
            r4(ix,iy,1,iz)= 0.
            r5(ix,iy,iz)  = 0.
            p(ix,iy,iz)   = 0.
            rj(ix,iy,iz)  = 0.
c
c ---------- flux velocities simple guess
c
c           uf(ix,iy,iz)  = u(ix,iy,iz)*gm(ix,iy,iz,1)
c           vf(ix,iy,iz)  = v(ix,iy,iz)*gm(ix,iy,iz,1)
c           wf(ix,iy,iz)  = u(ix,iy,iz)*gm(ix,iy,iz,2) +
c    +                      v(ix,iy,iz)*gm(ix,iy,iz,3)
            uf(ix,iy,iz)  = u(ix,iy,iz)
            vf(ix,iy,iz)  = v(ix,iy,iz)
            wf(ix,iy,iz)  = 0.0
         enddo
         enddo
      enddo
c
c ------------- set initial random field to be
c               divergence free
c
      idum = -1 - myid
      do iz=izs,ize
c
c ----------- ampv and ampt are max amplitudes of random
c             velocity and temperature fields
c             make sure ampv is set if free convection so
c             that we have motions at first time step
c             made dimensionless by ustar
c
         ampv = 0.01
         ampt = 0.01
c
c ------- simple random field scaled between -0.5 and 0.5
c
         sum_psi = 0.0
         do iy=iys,iye
         do ix=1,nnx
            psi(ix,iy) = ran1(idum)
            sum_psi = sum_psi + psi(ix,iy)
         enddo
         enddo
         sum_psi = sum_psi*fnxy
         call mpi_sum_xy(sum_psi,myid,iss,ise,1)
         do iy=iys,iye
         do ix=1,nnx
            psi(ix,iy) = psi(ix,iy) - sum_psi
            psix(ix,iy)     = psi(ix,iy)
            psiy(ix,iy,izs) = psi(ix,iy)
         enddo
         enddo
         call xderivp(psix(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call yd_mpi(psiy(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
         vmaxx = 0.0
         do iy=iys,iye
         do ix=1,nnx
            vmag = sqrt(psix(ix,iy)**2 + psiy(ix,iy,izs)**2)
            if(vmag .gt. vmaxx) vmaxx = vmag
         enddo
         enddo
         facv = ampv/vmaxx
c
         if (z(iz) .le. 0.2) then
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz)   = u(ix,iy,iz) - psiy(ix,iy,izs)*facv
            v(ix,iy,iz)   = v(ix,iy,iz) + psix(ix,iy)*facv
            t(ix,iy,1,iz) = t(ix,iy,1,iz) + psi(ix,iy)*ampt
         enddo
         enddo
         endif
c
         if(z(iz) .le. 0.2) then
         do iy=iys,iye
         do ix=1,nnx
c           e(ix,iy,iz) = 0.4*(1.0 - z(iz)/zl)**3
            e(ix,iy,iz) = 0.05*(1.0 - z(iz)/zl)**3
         enddo
         enddo
         endif
c
c ---------- check divergence of initial field
c
         do iy=iys,iye
         do ix=1,nnx
            uxx(ix,iy)     = uf(ix,iy,iz)
            vyy(ix,iy,izs) = vf(ix,iy,iz)
         enddo
         enddo
         call xderivp(uxx(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call yd_mpi(vyy(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
         do iy=iys,iye
         do ix=1,nnx
            divz(iz) = divz(iz) + (uxx(ix,iy) + vyy(ix,iy,izs))**2
         enddo
         enddo
         divz(iz) = divz(iz)*fnxy
c
c -------- end z loop
c
      enddo
c
      call mpi_sum_z(divz(1),i_root,myid,nnz,1)
c
      write(nprt,6000)
 6000 format(' check of divergence for initial state',/,
     +         ' iz ',5x,' divergence')
      write(nprt,6100) (iz,divz(iz),iz=izs,ize)
 6100 format(i5,e15.6)
c
      return
      end
      subroutine random_f
c
c ---------- example of using given (sparse) initial
c            sounding profiles (FIX for ncpu_s).
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real psi(nnx,iys:iye), psix(nnx,iys:iye),
     +     psiy(nnx,iys:iye,izs:izs), uxx(nnx,iys:iye),
     +     vyy(nnx,iys:iye,izs:izs)
c
      parameter (nt=12, nz=11)
      real zg(nz), u_i(nz,nt), v_i(nz,nt), theta_i(nz,nt)
      real ui_temp(nz), vi_temp(nz), ti_temp(nz)
      real time_g(nt)
c
      data time_g /
     +  0.0000E+00,  0.3600E+04,  0.7200E+04,  0.1080E+05,  0.1440E+05,
     +  0.1800E+05,  0.2160E+05,  0.2520E+05,  0.2880E+05,  0.3240E+05,
     +  0.3600E+05,  0.3960E+05
     +/
      data zg /
     +  0.1000E+02,  0.3000E+02,  0.5500E+02,  0.9000E+02,  0.1400E+03,
     +  0.2150E+03,  0.3300E+03,  0.5000E+03,  0.7500E+03,  0.1100E+04,
     +  0.1600E+04
     +/
      data u_i /
     + -0.1510E+01, -0.1560E+01, -0.1580E+01, -0.1580E+01, -0.1560E+01,
     + -0.1530E+01, -0.1510E+01, -0.9000E+00, -0.1390E+01, -0.1220E+01,
     + -0.5100E+00,
     + -0.1090E+01, -0.1110E+01, -0.1120E+01, -0.1120E+01, -0.1030E+01,
     + -0.9900E+00, -0.9500E+00, -0.6200E+00, -0.1230E+01, -0.9400E+00,
     +  0.2800E+00,
     + -0.9100E+00, -0.9200E+00, -0.9100E+00, -0.9000E+00, -0.8800E+00,
     + -0.8400E+00, -0.8000E+00, -0.6500E+00, -0.1510E+01, -0.1070E+01,
     +  0.2400E+00,
     + -0.8900E+00, -0.8900E+00, -0.8900E+00, -0.8800E+00, -0.8700E+00,
     + -0.8500E+00, -0.8100E+00, -0.7000E+00, -0.1830E+01, -0.8400E+00,
     +  0.3500E+00,
     + -0.1250E+01, -0.1260E+01, -0.1260E+01, -0.1250E+01, -0.1240E+01,
     + -0.1220E+01, -0.1160E+01, -0.8800E+00, -0.1980E+01, -0.1900E+00,
     +  0.7500E+00,
     + -0.1800E+01, -0.1810E+01, -0.1820E+01, -0.1820E+01, -0.1800E+01,
     + -0.1780E+01, -0.1710E+01, -0.1150E+01, -0.1960E+01,  0.3900E+00,
     +  0.9200E+00,
     + -0.2110E+01, -0.2130E+01, -0.2140E+01, -0.2140E+01, -0.2130E+01,
     + -0.2110E+01, -0.2050E+01, -0.9300E+00, -0.1400E+01,  0.8800E+00,
     +  0.9600E+00,
     + -0.2250E+01, -0.2280E+01, -0.2290E+01, -0.2300E+01, -0.2290E+01,
     + -0.2260E+01, -0.2070E+01, -0.4000E-01, -0.1600E+00,  0.1440E+01,
     +  0.1190E+01,
     + -0.2160E+01, -0.2200E+01, -0.2220E+01, -0.2220E+01, -0.2220E+01,
     + -0.2190E+01, -0.1610E+01,  0.1470E+01,  0.1420E+01,  0.2050E+01,
     +  0.1610E+01,
     + -0.2230E+01, -0.2270E+01, -0.2290E+01, -0.2300E+01, -0.2300E+01,
     + -0.2260E+01, -0.1350E+01,  0.2480E+01,  0.2380E+01,  0.2320E+01,
     +  0.1740E+01,
     + -0.1890E+01, -0.1930E+01, -0.1950E+01, -0.1950E+01, -0.1940E+01,
     + -0.1890E+01, -0.1120E+01,  0.3010E+01,  0.3030E+01,  0.2800E+01,
     +  0.2000E+01,
     + -0.1210E+01, -0.1230E+01, -0.1240E+01, -0.1230E+01, -0.1210E+01,
     + -0.1140E+01, -0.4600E+00,  0.3320E+01,  0.3510E+01,  0.3420E+01,
     +  0.2340E+01
     +/
      data v_i /
     +  0.4800E+00,  0.5100E+00,  0.5300E+00,  0.5700E+00,  0.6900E+00,
     +  0.7300E+00,  0.7600E+00,  0.1410E+01, -0.4200E+00, -0.3060E+01,
     + -0.3500E+01,
     +  0.7800E+00,  0.8100E+00,  0.8400E+00,  0.8900E+00,  0.1060E+01,
     +  0.1110E+01,  0.1130E+01,  0.1190E+01, -0.1040E+01, -0.2900E+01,
     + -0.3440E+01,
     +  0.3000E+00,  0.3200E+00,  0.3400E+00,  0.3800E+00,  0.4800E+00,
     +  0.5300E+00,  0.5800E+00,  0.5300E+00, -0.1330E+01, -0.2040E+01,
     + -0.2830E+01,
     + -0.2700E+00, -0.2600E+00, -0.2400E+00, -0.2200E+00, -0.1800E+00,
     + -0.1300E+00, -0.5000E-01,  0.1000E+00, -0.1170E+01, -0.1100E+01,
     + -0.2370E+01,
     + -0.5500E+00, -0.5400E+00, -0.5300E+00, -0.5100E+00, -0.4800E+00,
     + -0.4100E+00, -0.2600E+00,  0.1700E+00, -0.4200E+00, -0.2200E+00,
     + -0.2080E+01,
     + -0.2700E+00, -0.2600E+00, -0.2500E+00, -0.2400E+00, -0.2100E+00,
     + -0.1600E+00, -0.1000E-01,  0.8500E+00,  0.9700E+00,  0.3500E+00,
     + -0.2250E+01,
     +  0.5300E+00,  0.5400E+00,  0.5600E+00,  0.5700E+00,  0.6000E+00,
     +  0.6500E+00,  0.7600E+00,  0.1960E+01,  0.2280E+01,  0.3600E+00,
     + -0.2590E+01,
     +  0.1590E+01,  0.1630E+01,  0.1650E+01,  0.1680E+01,  0.1720E+01,
     +  0.1780E+01,  0.2010E+01,  0.3260E+01,  0.3110E+01,  0.1600E+00,
     + -0.2580E+01,
     +  0.2560E+01,  0.2620E+01,  0.2660E+01,  0.2690E+01,  0.2740E+01,
     +  0.2830E+01,  0.3400E+01,  0.4030E+01,  0.3030E+01, -0.7000E-01,
     + -0.2320E+01,
     +  0.3500E+01,  0.3600E+01,  0.3650E+01,  0.3700E+01,  0.3750E+01,
     +  0.3860E+01,  0.4580E+01,  0.4100E+01,  0.2450E+01,  0.6000E-01,
     + -0.1770E+01,
     +  0.4500E+01,  0.4640E+01,  0.4700E+01,  0.4760E+01,  0.4830E+01,
     +  0.4930E+01,  0.5420E+01,  0.3960E+01,  0.2000E+01,  0.5000E+00,
     + -0.1150E+01,
     +  0.5290E+01,  0.5470E+01,  0.5550E+01,  0.5620E+01,  0.5690E+01,
     +  0.5790E+01,  0.6070E+01,  0.4000E+01,  0.1910E+01,  0.9700E+00,
     + -0.5600E+00
     +/
      data theta_i /
     +  0.2936E+03,  0.2936E+03,  0.2937E+03,  0.2937E+03,  0.2938E+03,
     +  0.2942E+03,  0.2948E+03,  0.2980E+03,  0.3027E+03,  0.3092E+03,
     +  0.3186E+03,
     +  0.2937E+03,  0.2937E+03,  0.2937E+03,  0.2938E+03,  0.2939E+03,
     +  0.2942E+03,  0.2946E+03,  0.2978E+03,  0.3024E+03,  0.3090E+03,
     +  0.3184E+03,
     +  0.2939E+03,  0.2939E+03,  0.2939E+03,  0.2939E+03,  0.2939E+03,
     +  0.2941E+03,  0.2944E+03,  0.2976E+03,  0.3023E+03,  0.3089E+03,
     +  0.3182E+03,
     +  0.2940E+03,  0.2939E+03,  0.2939E+03,  0.2939E+03,  0.2939E+03,
     +  0.2941E+03,  0.2943E+03,  0.2975E+03,  0.3022E+03,  0.3087E+03,
     +  0.3181E+03,
     +  0.2940E+03,  0.2940E+03,  0.2939E+03,  0.2939E+03,  0.2939E+03,
     +  0.2940E+03,  0.2942E+03,  0.2974E+03,  0.3021E+03,  0.3086E+03,
     +  0.3180E+03,
     +  0.2941E+03,  0.2940E+03,  0.2940E+03,  0.2940E+03,  0.2941E+03,
     +  0.2941E+03,  0.2941E+03,  0.2973E+03,  0.3019E+03,  0.3085E+03,
     +  0.3179E+03,
     +  0.2942E+03,  0.2942E+03,  0.2942E+03,  0.2942E+03,  0.2941E+03,
     +  0.2941E+03,  0.2941E+03,  0.2973E+03,  0.3020E+03,  0.3086E+03,
     +  0.3179E+03,
     +  0.2943E+03,  0.2943E+03,  0.2943E+03,  0.2943E+03,  0.2943E+03,
     +  0.2943E+03,  0.2943E+03,  0.2975E+03,  0.3022E+03,  0.3087E+03,
     +  0.3181E+03,
     +  0.2945E+03,  0.2945E+03,  0.2945E+03,  0.2945E+03,  0.2945E+03,
     +  0.2944E+03,  0.2946E+03,  0.2978E+03,  0.3025E+03,  0.3090E+03,
     +  0.3184E+03,
     +  0.2947E+03,  0.2947E+03,  0.2947E+03,  0.2947E+03,  0.2946E+03,
     +  0.2946E+03,  0.2949E+03,  0.2980E+03,  0.3027E+03,  0.3093E+03,
     +  0.3187E+03,
     +  0.2949E+03,  0.2949E+03,  0.2949E+03,  0.2948E+03,  0.2948E+03,
     +  0.2948E+03,  0.2950E+03,  0.2982E+03,  0.3028E+03,  0.3094E+03,
     +  0.3188E+03,
     +  0.2950E+03,  0.2950E+03,  0.2950E+03,  0.2950E+03,  0.2950E+03,
     +  0.2950E+03,  0.2950E+03,  0.2982E+03,  0.3029E+03,  0.3095E+03,
     +  0.3188E+03
     +/
c
      save time_g, zg, u_i, v_i, theta_i
c
c --------- find time location of initial profiles
c
      call lterp(nt,time_g,t_factor,jt,jtp1,t_weit)
c
      do iz=1,nz
         ui_temp(iz) = u_i(iz,jt)*(1.0 - t_weit) +
     +                 u_i(iz,jtp1)*t_weit
         vi_temp(iz) = v_i(iz,jt)*(1.0 - t_weit) +
     +                 v_i(iz,jtp1)*t_weit
         ti_temp(iz) = theta_i(iz,jt)*(1.0 - t_weit) +
     +                 theta_i(iz,jtp1)*t_weit
      enddo
c
c ----------- interpolate vertically
c
      do iz=izs,ize
         call lterp(nz,zg,zz(iz),kk,kkp1,v_weit)
         u_temp = ui_temp(kk)*(1.0 - v_weit) +
     +            ui_temp(kkp1)*v_weit
         v_temp = vi_temp(kk)*(1.0 - v_weit) +
     +            vi_temp(kkp1)*v_weit
         theta_temp = ti_temp(kk)*(1.0 - v_weit) +
     +            ti_temp(kkp1)*v_weit
c
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz)   = u_temp
            v(ix,iy,iz)   = v_temp
            t(ix,iy,1,iz) = theta_temp
            w(ix,iy,iz)   = 0.
            r1(ix,iy,iz)  = 0.
            r2(ix,iy,iz)  = 0.
            r3(ix,iy,iz)  = 0.
            r4(ix,iy,1,iz)= 0.
            r5(ix,iy,iz)  = 0.
         enddo
         enddo
      enddo
c
c ------------- set initial random field to be
c               divergence free
c
      idum = -1
      do iz=izs,ize
         if (iz.le.8) then
c
c ----------- ampv and ampt are max amplitudes of random
c             velocity and temperature fields
c
         ampv = 0.5
         ampt = 0.1
c
c ------- simple random field scaled between 0 and 1
c
         do iy=iys,iye
         do ix=1,nnx
            psi(ix,iy) = ran1(idum)
         enddo
         enddo
c
         do iy=iys,iye
         do ix=1,nnx
            psix(ix,iy) = psi(ix,iy)
            psiy(ix,iy,izs) = psi(ix,iy)
         enddo
         enddo
         call xderivp(psix(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call yd_mpi(psiy(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
c
         vmaxx = 0.0
         do iy=iys,iye
         do ix=1,nnx
            vmag = sqrt(psix(ix,iy)**2 + psiy(ix,iy,izs)**2)
            if(vmag .gt. vmaxx) vmaxx = vmag
         enddo
         enddo
         facv = ampv/vmaxx
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz)   = u(ix,iy,iz) - psiy(ix,iy,izs)*facv
            v(ix,iy,iz)   = v(ix,iy,iz) + psix(ix,iy)*facv
            t(ix,iy,1,iz) = t(ix,iy,1,iz) + psi(ix,iy)*ampt
            e(ix,iy,iz)   = 1.0
         enddo
         enddo
         endif
c
c ---------- check divergence of initial field
c
         do iy=iys,iye
         do ix=1,nnx
            uxx(ix,iy)     = u(ix,iy,iz)
            vyy(ix,iy,izs) = v(ix,iy,iz)
         enddo
         enddo
         call xderivp(uxx(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call yd_mpi(vyy(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
         do iy=iys,iye
         do ix=1,nnx
            divz(iz) = divz(iz) + (uxx(ix,iy) + vyy(ix,iy,izs))**2
         enddo
         enddo
         divz(iz) = divz(iz)*fnxy
c
c -------- end z loop
c
      enddo
c
      call mpi_sum_z(divz(1),i_root,myid,nnz,1)
c
      write(nprt,6000)
 6000 format(' check of divergence for initial state',/,
     +         ' iz ',5x,' divergence')
      write(nprt,6100) (iz,divz(iz),iz=izs,ize)
 6100 format(i5,e15.6)
c
c ------------ fix for baroclinic and subsidence effects !!
c
c     do iz=izs,ize
c        ug(iz)=ugcont
c        vg(iz)=vgcont
c        if (.not.(ibrcl.eq.1)) go to 19988
c        if (.not.(iz.le.izi)) go to 19987
c        ug(iz)=0.
c        vg(iz)=0.
c 19987    continue
c 19988    continue
c        zz2=zz(iz)
c        wls(iz)=-divgls*zz2
c        if (.not.(iz.eq.1)) go to 19986
c        do ix=1,nnx
c        uls(ix)=divgls*(dx*float(ix-1)-xl*.5)
c        enddo
c     enddo
c     write(nprt,9)(uls(ix),ix=1,nnx)
c  9  format(1x,8e12.3)
c 19986 continue
c
      return
      end
      subroutine randoc
c
c -------- random initial conditions for an
c          ocean simulation
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real psi(nnx,iys:iye), psix(nnx,iys:iye),
     +     psiy(nnx,iys:iye,izs:izs), uxx(nnx,iys:iye),
     +     vyy(nnx,iys:iye,izs:izs)
c
      izi=(5*nnz)/20
      zi=z(izi)
      tmixed = 283.0
      do iz=izs,ize
         if (iz.le.izi) then
            do iy=iys,iye
            do ix=1,nnx
               u(ix,iy,iz)   = ugcont-ugal
               v(ix,iy,iz)   = vgcont
               w(ix,iy,iz)   = 0.0
               t(ix,iy,1,iz) = tmixed
               e(ix,iy,iz)   = 0.0
            enddo
            enddo
         endif
         if (iz.gt.izi) then
            do iy=iys,iye
            do ix=1,nnx
               u(ix,iy,iz)   = ugcont-ugal
               v(ix,iy,iz)   = vgcont
               w(ix,iy,iz)   = 0.0
               t(ix,iy,1,iz) = tmixed + dtdzf(1)*(zz(iz)-zi)
               e(ix,iy,iz)   = 0.0
            enddo
            enddo
         endif
         do iy=iys,iye
         do ix=1,nnx
            w(ix,iy,iz)    = 0.0
            r1(ix,iy,iz)   = 0.0
            r2(ix,iy,iz)   = 0.0
            r3(ix,iy,iz)   = 0.0
            r4(ix,iy,1,iz) = 0.0
            r5(ix,iy,iz)   = 0.0
         enddo
         enddo
      enddo
c
c ------------- set initial random field to be
c               divergence free
c
      idum = -1
      do iz=izs,ize
      if (iz.le.4) then
c
c ----------- ampv and ampt are max amplitudes of random
c             velocity and temperature fields
c
         ampv = 0.01
c        ampt = 0.00
         ampt = 0.0001
c
c ------- simple random field scaled between 0 and 1
c
         do iy=iys,iye
         do ix=1,nnx
            psi(ix,iy) = ran1(idum)
         enddo
         enddo
c
         do iy=iys,iye
         do ix=1,nnx
            psix(ix,iy) = psi(ix,iy)
            psiy(ix,iy,izs) = psi(ix,iy)
         enddo
         enddo
         call xderivp(psix(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call yd_mpi(psiy(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
c
         vmaxx = 0.0
         do iy=iys,iye
         do ix=1,nnx
            vmag = sqrt(psix(ix,iy)**2 + psiy(ix,iy,izs)**2)
            if(vmag .gt. vmaxx) vmaxx = vmag
         enddo
         enddo
         facv = ampv/vmaxx
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,iz) = u(ix,iy,iz) - psiy(ix,iy,izs)*facv
            v(ix,iy,iz) = v(ix,iy,iz) + psix(ix,iy)*facv
            t(ix,iy,1,iz) = t(ix,iy,1,iz) + psi(ix,iy)*ampt
            e(ix,iy,iz) = 0.0001
         enddo
         enddo
      endif
c
c ---------- check divergence of initial field
c
         do iy=iys,iye
         do ix=1,nnx
            uxx(ix,iy) = u(ix,iy,iz)
            vyy(ix,iy,izs) = v(ix,iy,iz)
         enddo
         enddo
         call xderivp(uxx(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
         call yd_mpi(vyy(1,iys,izs),trigx(1,2),yk(1),
     +           nnx,nny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
         do iy=iys,iye
         do ix=1,nnx
            divz(iz) = divz(iz) + (uxx(ix,iy) + vyy(ix,iy,izs))**2
         enddo
         enddo
         divz(iz) = divz(iz)*fnxy
c
c -------- end z loop
c
      enddo
c
      call mpi_sum_z(divz(1),i_root,myid,nnz,1)
c
      write(nprt,6000)
 6000 format(' check of divergence for initial state',/,
     +         ' iz ',5x,' divergence')
      write(nprt,6100) (iz,divz(iz),iz=izs,ize)
 6100 format(i5,e15.6)
c
      do iz=izs,ize
         ug(iz)=ugcont
         vg(iz)=vgcont
      enddo
c
      return
      end
      subroutine get_fields
c
c ----------- special routine to read just 3d fields
c             as an initial guess, easy to customize
c             if missing data, etc..
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
#if defined(SWAP)
      use module_byteswap
#endif
      include 'mpif.h'
c
      integer status(mpi_status_size), ierr
      integer(kind=mpi_offset_kind) :: offset, disp
      integer(kind=k8)              :: nsize, nsize2
      real, allocatable, dimension(:,:,:) :: temp
      logical there
c
      allocate(temp(nvar+4,nnx,iys:iye))
c
c ---------- input file to read from
c
c --------------------- get restart file from local directory
c                       reuse unit number
c
      close(nvel)
c
      inquire(file=path_ran,exist=there)
      if(there) then
         if(l_root) write(6,6001) path_ran
      else
         if(l_root) write(6,6005) path_ran
         stop
      endif
c
c ---- open file
c
      call mpi_file_open(mpi_comm_world, path_ran,
     +                   mpi_mode_create+mpi_mode_rdwr,
     +                   mpi_info_null, nvel, ierr)
c
c ---- set file view
c
      disp = 0
      call mpi_file_set_view(nvel,disp,mpi_real8,mpi_real8,
     +                      'native',mpi_info_null,ierr)
c
c ------------ read 3d fields, make rhs*8
c              grid speed is in file but don't need it for this restart
c
      nsize  = int(nvar+4,k8)*nnx*nny
      nsize2 = int(nvar+4,k8)*nnx*(iys-1)
      n_read = (nvar+4)*nnx*(iye+1-iys)
c
      do k=izs,ize
         offset = int((k-1),k8)*nsize + nsize2
         call mpi_file_read_at_all(nvel,offset,temp,n_read,
     +                              mpi_real8,status,ierr)
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
         enddo
         enddo
         do is = 1,nscl
            do j = iys,iye
            do i = 1,nnx
              t(i,j,is,k) = temp(6+is,i,j)
            enddo
            enddo
         enddo
c
      enddo
c
c ---- close file
c
      call mpi_file_close(nvel, ierr)
c
      deallocate(temp)
c
      do k=1,nnz
         ug(k) = ugcont
         vg(k) = vgcont
      enddo
c
      do iz=izs,ize
c
c ---------------- initial guess for pressure
c
         do iy=iys,iye
         do ix=1,nnx
            p(ix,iy,iz) = 0.0
         enddo
         enddo
      enddo
c
      return
c ---------------------------- process errors
  100 continue
      write(6,9000) path_ran, nvel
      call mpi_finalize(ierr)
      stop
c
 9992 continue
      write(6,6100) nvel,iz
      call mpi_finalize(ierr)
      stop
c ---------------------
 6001 format(' SR. GET_FIELDS: FILE READ FOR INITIALIZATION = ',a80)
 6005 format(' 6005, SR. GET_FIELDS: cannot find restart file = ',a80)
 6100 format(' SR. GET_FIELDS: file read error on unit number = ',i2,/,
     +       '               at iz = ',i4)
 9000 format(' 9000, SR. GET_FIELDS: cannot open file =',a80,/,
     +       ' to unit number = ',i2)
      end
      subroutine forcing
c
c ----------- update surface temperature based on a
c             constant cooling rate
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      include 'mpif.h'
c
      t_surface = t_surf_i - c_rate*time
c
      do j=iys,iye
      do i=1,nnx
         t_grnd(i,j,1) = t_surface
      enddo
      enddo
c
      return
      end
      subroutine set_sav(it,istart)
c
c --------- time stepping layout
c
c ----------------------------------------------------------------------
c      it-1          it             it+1
c       0 ----------> 1               2
c       *             *               *
c             dt             dt1             dt_new
c     time
c ----------------------------------------------------------------------
c
c                it-1 = 1 ---------> it=2           3
c                     *               *             *
c                         dt = dt1       dt1 = dt_new
c                                                update dt_new from it = 1
c                    time
c ----------------------------------------------------------------------
c
c --------- dt  = current increment between [it-1,it]
c           dt1 = next increment between [it,it+1]
c           dt_new from flow field at end of current it+1 step
c
      use pars
      use fields
      use con_data
      use con_stats
      use netcdf_io
c
      data ionce /0/
      save ionce
c
      if(it .ne. istart) then
c
c ---------- increment time if not first time through
c
         time=time+dt
      endif
c
      if(istart .eq. 0) then
         dt1 = dt_new
      endif
c
      it=it+1
c
      dt    = dt1
      dt1   = dt_new
      mnout = (mod(it,imean).eq.0).or. (it.eq.1)
      mtape = (mod(it,itape).eq.0)
      micut = (mod(it,itcut).eq.0)
      mviz  = (mod(it,iviz).eq.0)
c
c ---------- decide frequency history fields are saved
c
      if(ihst .lt. 0) then
         mhis = .false.
      else
         mhis = (mod(it,ihst).eq.0 .and. it .ge. it_his)
      endif
      mtrans = (mod(it+1,ihst).eq.0)
c
c ---------- decide frequency viz fields are saved
c
      if(iviz .lt. 0) then
         msave_v = .false.
      else
         msave_v = (mod(it,iviz).eq.0 .and. it .ge. it_viz)
      endif
c
c ---------- decide whether velocity fields are saved
c
      msave = .false.
      if(it .ge. itstr .and. mtape) then
         itn=itn+1
         msave = .true.
         call get_output_filenames
      endif
c
c --------- decide whether history files are to be saved
c
      if((ihst .gt. 0) .and. (it .ge. it_his_nxt)) then
         if (inetcdf .eq. 1) then
            call open_his_netcdf
         else
            call open_his(it)
         end if
         it_his_nxt = it_his_nxt + itape
      endif
c
c --------- decide whether viz files are to be saved
c
      if((iviz .gt. 0) .and. (it .ge. it_viz_nxt)) then
c     if(l_root) then
c       write(6,4501) iviz, it_viz_nxt
c4501   format(' 4501: iviz = ',i5,' it_viz_nxt = ',i5)
c     endif
         call viz_output_filename(it)
         it_viz_nxt = it_viz_nxt + itape
      endif
c
c ---------- if first time thru and moving mesh
c            setup the grid and set t_zero = reference time
c
      if(i_move_step .lt. 0 .or.
     +   it .lt. i_move_step) then
        imesh  = 0
        t_zero = time
        go to 999
      endif
c
      if(it .eq. i_move_step .and. imesh .eq. 0) then
c
c -------- first time thru set W = grid speed everywhere
c
        imesh    = 1
        t_zero   = time
        t_stage  = time + dt*etas(1)
        t_stage2 = time + dt*etas(2)
        call start_mesh(it,1)
        call grid_speed(it,0)
        do k=izs-1,ize+1
           do j=iys,iye
           do i=1,nnx
              wf(i,j,k) = gm2(i,j,k,5) + wf(i,j,k)
c             wf(i,j,k) = gm2(i,j,k,5)
           enddo
           enddo
        enddo
c
      elseif(it .gt. i_move_step .and. imesh .eq. 0) then
c
c ---------- restart with moving grid
c
        imesh    = 1
        t_stage  = time + dt*etas(1)
        t_stage2 = time + dt*etas(2)
        call start_mesh(it,1)
        call grid_speed(it,0)
      endif
c
  999 continue
c
c     write(nprt,9901) it, i_move_step, t_zero
c9901 format(' 9901 set_sav',/,
c    +       ' it = ',i5,' i_move_step = ',i8,' t_zero = ',e15.6)
c
      return
      end
      subroutine pbltop(itop)
c
c ---------- get estimate of pbl top
c
c            method = 0, min of wt flux
c                        (good for buoyancy cases)
c            method = 1, uw flux less than critical value
c                        (good for ekman cases)
c            method = 2, running t average exceeds criterion
c                        (good for neutral cases with capping
c                         inversions)
c            method = 3, maximum gradient in temperature field
c                        (good for finding local zi see jas paper)
c                        with minimum search height (sr. setup)
c
c ------------ if method uses average statistics then only root
c              process need find zi
c
      use pars
      use fields
      use con_data
      use con_stats
      real trun(maxnz)
      include 'mpif.h'
      real gradloc(2,nnx,nny), gradmax(2,nnx,nny)
      external get_zi
c
      if(method .le. 2 .and. l_root) then
c
      sgn = 1.0
      if(iocean .eq. 1) sgn = 1.0
      if (method .le. 0 .or. method .gt. 2) then
         itop=1
         wttot=wtle(1,1)+wtsb(1,1)
         wtmin=wttot*sgn
         do iz=2,nnz
            wttot=(wtle(iz,1)+wtsb(iz,1))*sgn
            if (wttot.le.wtmin) then
               itop=iz
               wtmin=wttot
            endif
         enddo
         zi=z(itop)
      else if (method .eq. 1) then
         itop = 1
         crit = 0.05
         uwsf = utau*utau
         do iz=1,nnzm1
               uwtot = (uwle(iz) + uwsb(iz))**2 +
     $                 (vwle(iz) + vwsb(iz))**2
               uwtot = sqrt(uwtot)
               if(uwtot/uwsf .gt. crit) then
                  itop=iz
               endif
         enddo
         zi=z(itop)
      else if (method .eq. 2) then
         trun(1) = txym(1,1)
         do iz=2,nnz
             weight = z(iz-1)/z(iz)
             trun(iz) = trun(iz-1)*weight + (1.0-weight)*txym(iz,1)
         enddo
         itop = 1
         tcrit = 0.25
         if(iocean .eq. 1) tcrit = 0.1
         do iz=2,nnz
                if(txym(iz,1) .gt. (trun(iz) + tcrit)) then
                  itop = iz
                  go to 320
                endif
         enddo
  320    continue
         zi=z(itop)
      endif
      do iy=1,nny
      do ix=1,nnx
         gradmax(2,ix,iy) = zi
      enddo
      enddo
c
c ----------- use gradient method, every process computes
c
      elseif(method .eq. 3) then
c
c ---------------- get local zi from gradient in temperaure field
c
c     dz_i = dzu_i(izs+1)
c     do iy=1,nny
c     do ix=1,nnx
c        gradloc(1,ix,iy) = (t(ix,iy,1,izs+1) - t(ix,iy,1,izs))*dz_i
c        gradloc(2,ix,iy) = z(izs)
c     enddo
c     enddo
c
c ------- similar to zeroing the stat array in sr. mean_stat
c
      do iy=1,nny
      do ix=1,nnx
         gradloc(1,ix,iy) = 0.0
         gradloc(2,ix,iy) = z(iz_min)
      enddo
      enddo
c
c ------------- now all z in this process
c
      if(iz_min .le. ize) then
      do iz=max(izs,iz_min),ize
         izp1 = iz + 1
         do iy=iys,iye
         do ix=1,nnx
            grad = (t(ix,iy,1,izp1) - t(ix,iy,1,iz))*dzu_i(izp1)
            if(grad .gt. gradloc(1,ix,iy)) then
               gradloc(1,ix,iy) = grad
               gradloc(2,ix,iy) = z(iz)
            endif
         enddo
         enddo
      enddo
      endif
c
c     call mpi_reduce(gradloc,gradmax,2*nnx*nny,mpi_real8,ziloc,
c    +                i_root,mpi_comm_world,ierror)
c
c ----------- alternate version using already defined function in mpi
c             passes 2 real8 variables
c
      call mpi_reduce(gradloc,gradmax,nnx*nny,mpi_2double_precision,
     +                mpi_maxloc,i_root,mpi_comm_world,ierror)
c
c ------------ get average on root process
c
      if(l_root) then
         zi_avg = 0.0
         do iy=1,nny
         do ix=1,nnx
            zi_avg = zi_avg + gradmax(2,ix,iy)
         enddo
         enddo
         zi = zi_avg*fnxy
         itop = 1
         do k=1,nnz-1
            if(zi .gt. z(k) .and. zi .le. z(k+1)) then
              itop = k
            endif
         enddo
      endif
c
      endif
c
c -------- send average zi everywhere
c
      call mpi_bcast(zi,1,mpi_real8,
     +              i_root,mpi_comm_world,ierr)
c
      if(iocean .ne. 1) then
         do iz=1,nnz
            if(zi .ge. z(iz) .and.
     +         zi .lt. z(iz+1)) itop = iz
         enddo
      else
         do iz=1,nnz
            if(zi .le. z(iz) .and.
     +         zi .gt. z(iz+1)) itop = iz
         enddo
      endif
c
c ---------- get zi/L
c
      hol = zi/amonin
c
c     if(l_root) write(6,7001) myid,zi,itop
 7001 format(' 7001 in pbltop myid = ',i4,' zi = ',e15.6,
     +       ' itop = ',i3)
c
      return
      end
      subroutine get_zi(gradmax,gradout,len,itype)
c
      use pars
      real gradmax(*), gradout(*)
c
c     write(nprt,2001) myid, len
c2001 format(' 2001 in get_zi myid = ',i4,' len = ',i8)
c     write(nprt,2002) (i,gradmax(i),gradmax(i+1),i=1,len,2)
c2002 format(' i ',5x,' grad ',5x,' location ',/,
c    +      (i5,2e15.6))
c
      do i=1,len,2
         if(gradmax(i) .gt. gradout(i)) then
              gradout(i)   = gradmax(i)
              gradout(i+1) = gradmax(i+1)
         endif
      enddo
c
      return
      end
      subroutine print(lu,it,iz_strt,iz_end)
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      write(lu,4000)
 4000 format(30X,' --- SOLUTION ---')
      write(lu,4100) it,time,dt,zi,tsfcc(1),uwsfc,vwsfc,wtsfc(1),
     +              zol, hol, amonin, ucfl, vcfl, wcfl, utau,
     +              t_ref, drg_x, drg_y, p_iter, p_err
 4100 format(' IT=',I7,5x,'TIME (s) = ',e15.8,',  DT(s) = ',e15.6,/,
     +       10x,'ZTOP = ',e15.6,
     +       ',  TSFC = ',e15.6,
     +       ',  UW = ',e15.6,',  VW = ',e15.6,/,10x,
     +       'WT = ',e15.6,',  Z1/L =',e15.6,
     +       ',  ZI/L = ',e15.6,',  MO-L = ',e15.6,/,10x,
     +       'U_cfl = ',e15.6,
     +       ',  V_cfl = ',e15.6,',  W_cfl = ',e15.6,
     +       ',  UTAU = ',e15.6,/,10x,
     +       'Theta Ref = ',e15.6,/,10x,'P-Drag_x  = ',e15.6,
     +       ',   P-Drag_y = ',e15.6,/,10x,
     +       'P iterations = ',e15.6,',   P error = ',e15.6)
      write(lu,4200)
 4200 format(//,20x,'--------- HORIZONTAL MEAN VALUES ---------- ',
     +       //,2x,'IZ',4x,'T_MEAN',7x,
     + 'DIVG',8X,'LE_KE',6X,'SGS_KE',7X,'LE_WT',6X,
     + 'SGS_WT',7X,'SHRZ',8X,'BUOY')
      do 19999 iz=iz_end,iz_strt,-1
         write(lu,4300)iz,txym(iz,1)-t_ref,divz(iz),
     +              englez(iz),eavg(iz),wtle(iz,1),
     +              wtsb(iz,1),shrz(iz),buyz(iz)
 4300    format(1X,I3,e12.4,7e12.4)
19999 continue
      write(lu,4400)tsfcc(1),wtsfc(1)
 4400 format('  SURFACE VALUES: TXYM = ',E15.6,'  WTSB = ',E15.6)
      if(iocean .eq. 1) then
         write(lu,4500) stokess,udrift,vdrift
 4500    format(/,' STOKESS = ',e12.4,' UDRIFT = ',e12.4,
     +          ' VDRIFT = ',e12.4)
      endif
      write(lu,4600) (iz,uxym(iz)+ugal,vxym(iz),uwle(iz),
     +       uwsb(iz),vwle(iz),vwsb(iz),pm(iz),iz=iz_strt,iz_end)
 4600 format(//,' IZ',5x,' UXYM + UGAL',8x,' VXYM',10x,' UWLE',11x,
     +          ' UWSB',11x,' VWLE',11x,' VWSB',11x,' P MEAN'
     +       ,/,(1x,i4,7(3x,e13.6)))
      if(ivis .eq. 1) then
         write(lu,4800) xksurf, nmatch, viscon, vise
 4800    format(//,' XKSURF = ',e15.6,' NMATCH = ',i4,/,
     +             ' VISCON = ',e15.6,' VISE = ',e15.6)
         write(lu,4700) (iz,dfac(iz),iz=iz_strt,iz_end)
 4700    format(//,'   IZ',5x,'  DFAC',/,(1x,i4,3x,e15.6))
      endif
c
c --------------- output additional scalars
c
      if(nscl .eq. 2) then
      write(lu,5005)tsfcc(2),wtsfc(2)
 5005 format(/,'  SURFACE VALUE: TXYM(2) =',e15.6,' WTSFC(2) = ',e15.6)
      write(lu,5100) (iz,txym(iz,2),wtle(iz,2),
     +              wtsb(iz,2),iz=iz_strt,iz_end)
 5100 format(//,' IZ',5x,' SCALAR-1 MEAN',8x,' WS1LE',10x,
     +          ' WS1SB',10x
     +       ,/,(1x,i4,3(3x,e13.6)))
c     else if (nscl .eq. 3) then
c     write(lu,5205)tsfcc(2),wtsfc(2),tsfcc(3),wtsfc(3)
c5205 format(/,'  SURFACE VALUE: TXYM(2) =',e15.6,' WTSFC(2) = ',e15.6,
c    +       /,'  SURFACE VALUE: TXYM(3) =',e15.6,' WTSFC(3) = ',e15.6)
c     write(lu,5200) (iz,txym(iz,2),txym(iz,3),wtle(iz,2),
c    +    wtsb(iz,2),wtle(iz,3),wtsb(iz,3),iz=iz_strt,iz_end)
c5200 format(//,' IZ',5x,' SCALAR-1 MEAN',8x,' SCALAR-2 MEAN',10x,
c    +          ' WS1LE',10x,' WS1SB',10x,' WS2LE',10x,' WS1SB'
c    +       ,/,(1x,i4,6(3x,e13.6)))
      endif
c
      return
      end
      subroutine xy_stats
c
c ------------ get statistics.
c              fluxes that appear in transformed equations
c              are based on (UF,VF,WF-Z_t) located at cell centers
c              and upper cell face.
c
c              f_avg = [\int \int f * fac dx dy]/area
c              fac   = [sqrt(zeta_x^2 + zeta_y^2 + zeta_z^2]/J
c
      use pars
      use fields
      use con_data
      use con_stats
      use particles
c
c ------- indices for indexing array stat(.,.)
c         js = number of non-scalar stats
c         ns = number of scalar stats
c
      parameter(js = 37, ns = 6, nstat = js + ns*nscl)
      real stat(1:nnz,nstat), buf(4)
c
c -------- stat(.,1) = u*u = ups
c          stat(.,2) = v*v = vps
c          stat(.,3) = w*w = wps
c          stat(.,4) = w**3 = wcube
c          stat(.,5) = w**4 = wfour
c          stat(.,6) = resolved tke at zu = englez
c          stat(.,7) = sgs e at zu = engsbz
c          stat(.,8) = sgs e at zu = eavg
c          stat(.,9) = resolved uw_f at zw = uwle
c          stat(.,10) = resolved vw_f at zw = vwle
c          stat(.,11) = resolved uw at zu = uwcle
c          stat(.,12) = resolved vw at zu = vwcle
c          stat(.,13) = form drag in x direction at zu
c          stat(.,14) = form drag in y direction at zu
c          stat(.,15) = x form drag at zu correlated with surface [p(z)*z_xi(1)]
c          stat(.,16) = y form drag at zu correlated with surface [p(z)*z_eta(1)]
c          stat(.,17) = resolved w*w_f flux at zw = wwle
c          stat(.,18) = particle mass in each horizontal layer
c          stat(.,19) = vpsum(1) (mean horizontal particle velocity)
c          stat(.,20) = vpsum(2) (mean spanwise particle velocity)
c          stat(.,21) = vpsum(3) (mean vertical particle velocity)
c          stat(.,22) = vpsqrsum(1)
c          stat(.,23) = vpsqrsum(2)
c          stat(.,24) = vpsqrsum(3)
c          stat(.,25) = upwp = up'*wp'
c          stat(.,26) = Tpsum
c          stat(.,27) = Tpsqrsum
c          stat(.,28) = wpTpsum
c          stat(.,29) = Tfsum
c          stat(.,30) = qfsum
c          stat(.,31) = radsum
c          stat(.,32) = rad2sum
c          stat(.,33) = qstarsum
c          stat(.,34) = zconc_U
c          stat(.,35) = pflux
c          stat(.,36) = pfluxup
c          stat(.,37) = pfluxdn
c          stat(.,m1) = resolved scalar flux w_ft at zw = wtle
c          stat(.,m2) = resolved scalar flux u_ft at zu = utle
c          stat(.,m3) = resolved scalar flux v_ft at zu = vtle
c          stat(.,m4) = scalar t*t at zu = tps
c          stat(.,m5) = scalar t*t*t at zu = tcube
c          stat(.,m6) = resolved scalar flux wct at zu = wctle
c
c --------- use a trick with mpi reduce over all z to get averages
c           by setting stat array = 0 for all z on each process
c
      do i=1,nstat
      do k=1,nnz
         stat(k,i) = 0.0
      enddo
      enddo
c
c -------- indices for scalars
c
      m1 = js
      m2 = js + nscl
      m3 = js + 2*nscl
      m4 = js + 3*nscl
      m5 = js + 4*nscl
      m6 = js + 5*nscl
c
      sgn = 1.0
      if(iocean .eq. 1 .and. iupwnd .eq. 1) sgn = -1.0
c
      do k=izs,ize
c
      kp2 = k + 2
      kp1 = k + 1
      km1 = k - 1
c
c ------- simple average along coordinate lines
c
      do j=iys,iye
      do i=1,nnx
         stat(k,1) = stat(k,1) + (u(i,j,k) - uxym(k))**2
         stat(k,2) = stat(k,2) + (v(i,j,k) - vxym(k))**2
         stat(k,3) = stat(k,3) + (w(i,j,k) - wxym(k))**2
         stat(k,4) = stat(k,4) + (w(i,j,k) - wxym(k))**3
         stat(k,5) = stat(k,5) + (w(i,j,k) - wxym(k))**4
         stat(k,6) = stat(k,6) +
     +                ((w(i,j,k)- wxym(k))**2 +
     +                 (u(i,j,k)- uxym(k))**2 +
     +                 (v(i,j,k)- vxym(k))**2)*0.5
c
         stat(k,7) = stat(k,7) + e(i,j,k)
         stat(k,8) = stat(k,8) + e(i,j,k)
         stat(k,9) = stat(k,9) + 0.5*(wf(i,j,k) - gm(i,j,k,5))*
     +               (u(i,j,k)-uxym(k) + u(i,j,kp1) - uxym(kp1))
         stat(k,10) = stat(k,10) + 0.5*(wf(i,j,k) - gm(i,j,k,5))*
     +               (v(i,j,k)-vxym(k) + v(i,j,kp1) - vxym(kp1))
         stat(k,11) = stat(k,11) + (w(i,j,k)-wxym(k))*
     +                           (u(i,j,k)-uxym(k))
         stat(k,12) = stat(k,12) + (w(i,j,k)-wxym(k))*
     +                             (v(i,j,k)-vxym(k))
c
c --------- get drag [p(z)*z_xi(z), p(z)*z_eta(z)]
c
         stat(k,13) = stat(k,13) + p(i,j,k)*gm(i,j,k,1)*gm(i,j,k,2)
         stat(k,14) = stat(k,14) + p(i,j,k)*gm(i,j,k,1)*gm(i,j,k,3)
c
c --------- get drag [p(z)*z_xi(1), p(z)*z_eta(1)]
c
         stat(k,15) = stat(k,15) + p(i,j,k)*gms(i,j,1)*gms(i,j,2)
         stat(k,16) = stat(k,16) + p(i,j,k)*gms(i,j,1)*gms(i,j,3)
c
         stat(k,17) = stat(k,17) + 0.5*(wf(i,j,k) - gm(i,j,k,5))*
     +               (w(i,j,k)-wxym(k) + w(i,j,kp1) - wxym(kp1))

         stat(k,18) = stat(k,18) + partcount(i,j,k)
         stat(k,19) = stat(k,19) + vpsum(i,j,k,1)
         stat(k,20) = stat(k,20) + vpsum(i,j,k,2)
         stat(k,21) = stat(k,21) + vpsum(i,j,k,3)
         stat(k,22) = stat(k,22) + vpsqrsum(i,j,k,1)
         stat(k,23) = stat(k,23) + vpsqrsum(i,j,k,2)
         stat(k,24) = stat(k,24) + vpsqrsum(i,j,k,3)
         stat(k,25) = stat(k,25) + upwp(i,j,k)
         stat(k,26) = stat(k,26) + Tpsum(i,j,k)
         stat(k,27) = stat(k,27) + Tpsqrsum(i,j,k)
         stat(k,28) = stat(k,28) + wpTpsum(i,j,k)
         stat(k,29) = stat(k,29) + Tfsum(i,j,k)
         stat(k,30) = stat(k,30) + qfsum(i,j,k)
         stat(k,31) = stat(k,31) + radsum(i,j,k)
         stat(k,32) = stat(k,32) + rad2sum(i,j,k)
         stat(k,33) = stat(k,33) + qstarsum(i,j,k)
         stat(k,34) = stat(k,34) + partcountU(i,j,k)
         stat(k,35) = stat(k,35) + partflux(i,j,k)
         stat(k,36) = stat(k,36) + partfluxup(i,j,k)
         stat(k,37) = stat(k,37) + partfluxdn(i,j,k)
      enddo
      enddo
c
c ------------ get scalar resolved fluxes and variances
c
      do l=1,nscl
         if(iupwnd .ne. 1 .or. k .eq. nnz) then
            do j=iys,iye
            do i=1,nnx
               stat(k,m1+l)=stat(k,m1+l) +
     +               0.5*(wf(i,j,k) - gm(i,j,k,5))*
     +          (t(i,j,l,k)-txym(k,l) + t(i,j,l,kp1) - txym(kp1,l))
            enddo
            enddo
         else
c
c ------------------- monotone fluxes
c
           do j=iys,iye
           do i=1,nnx
              wflx = wf(i,j,k) - gm(i,j,k,5)
              stat(k,m1+l) = stat(k,m1+l) +
     +    amax1(sgn*wflx,0.)*(t(i,j,l,k) +
     + rlim(t(i,j,l,kp1),t(i,j,l,k),t(i,j,l,km1))) +
     +    amin1(sgn*wflx,0.)*(t(i,j,l,kp1) +
     + rlim(t(i,j,l,k),t(i,j,l,kp1),t(i,j,l,kp2)))
           enddo
           enddo
         endif
         stat(k,m1+l)= sgn*stat(k,m1+l)
c
           do j=iys,iye
           do i=1,nnx
              stat(k,m6+l) = stat(k,m6+l) +
     +                       (w(i,j,k) - wxym(k))*
     +                       (t(i,j,l,k)-txym(k,l))
           enddo
           enddo
c
c ------------ get horizontal scalar resolved fluxes
c
         do j=iys,iye
         do i=1,nnx
            stat(k,m2+l) = stat(k,m2+l)+
     +               uf(i,j,k)*(t(i,j,l,k)-txym(k,l))
            stat(k,m3+l) = stat(k,m3+l)+
     +               vf(i,j,k)*(t(i,j,l,k)-txym(k,l))
         enddo
         enddo
c
c ------------------- scalar variances & higher moments
c
         do j=iys,iye
         do i=1,nnx
            stat(k,m4+l) = stat(k,m4+l) +
     +                (t(i,j,l,k) - txym(k,l))**2
            stat(k,m5+l) = stat(k,m5+l) +
     +                (t(i,j,l,k) - txym(k,l))**3
         enddo
         enddo
c
c ------ end scalar loop
c
      enddo
c
c ------ end z loop
c
      enddo
c
c -------- add partial sums and send it to all
c
      call mpi_sum_z(stat(1,1),i_root,myid,nnz*nstat,1)
c
c ------ fill arrays for printout and constant file
c
      do k=1,nnz
c
      ups(k)      =   stat(k,1)*fnxy
      vps(k)      =   stat(k,2)*fnxy
      wps(k)      =   stat(k,3)*fnxy
      wcube(k)    =   stat(k,4)*fnxy
      wfour(k)    =   stat(k,5)*fnxy
      englez(k)   =   stat(k,6)*fnxy
      engsbz(k)   =   stat(k,7)*fnxy
      eavg(k)     =   stat(k,8)*fnxy
      uwle(k)     =   stat(k,9)*fnxy
      vwle(k)     =   stat(k,10)*fnxy
      uwcle(k)    =   stat(k,11)*fnxy
      vwcle(k)    =   stat(k,12)*fnxy
      uw_tot(k)   =   uwle(k) + uwsb(k)
      vw_tot(k)   =   vwle(k) + vwsb(k)
      drg_z_x(k)  =   stat(k,13)*fnxy
      drg_z_y(k)  =   stat(k,14)*fnxy
      drg_zs_x(k) =  stat(k,15)*fnxy
      drg_zs_y(k) =  stat(k,16)*fnxy
      wwle(k)     =  stat(k,17)*fnxy
      zconc(k)    = stat(k,18)/xl/yl/dzw(k)
      zconc_U(k)  = stat(k,34)/xl/yl/dz_U
      vp1mean(k)  = stat(k,19)/stat(k,34)
      vp2mean(k)  = stat(k,20)/stat(k,34)
      vp3mean(k)  = stat(k,21)/stat(k,34)
      vp1msqr(k)  = sqrt(stat(k,22)/stat(k,34) - vp1mean(k)**2)
      vp2msqr(k)  = sqrt(stat(k,23)/stat(k,34) - vp2mean(k)**2)
      vp3msqr(k)  = sqrt(stat(k,24)/stat(k,34) - vp3mean(k)**2)
      upwpm(k)    = stat(k,25)/stat(k,34)-(vp1mean(k)*vp3mean(k))
      Tpmean(k)   = stat(k,26)/stat(k,34)
      Tpmsqr(k)   = sqrt(stat(k,27)/stat(k,34)-Tpmean(k)**2)
      wpTpm(k)    = stat(k,28)/stat(k,34) - (Tpmean(k)*vp3mean(k))
      Tfmean(k)   = stat(k,29)/stat(k,34)
      qfmean(k)   = stat(k,30)/stat(k,34)
      radmean(k)   = stat(k,31)/stat(k,34)
      rad2mean(k)   = stat(k,32)/stat(k,34)
      qstarmean(k)   = stat(k,33)/stat(k,34)
      pflux(k)   = stat(k,35)
      pfluxup(k)   = stat(k,36)
      pfluxdn(k)   = stat(k,37)

c
c ------------ get scalar resolved fluxes and variances
c
      do l=1,nscl
         wtle(k,l)   = stat(k,m1+l)*fnxy
         utle(k,l)   = stat(k,m2+l)*fnxy
         vtle(k,l)   = stat(k,m3+l)*fnxy
         tps(k,l)    = stat(k,m4+l)*fnxy
         tcube(k,l)  = stat(k,m5+l)*fnxy
         wctle(k,l)  = stat(k,m6+l)*fnxy
         wt_tot(k,l) = wtle(k,l) + wtsb(k,l)
      enddo
      enddo
c
c ---------- drag at first zu grid point
c
      drg_x = drg_z_x(1)
      drg_y = drg_z_y(1)
c
      return
      end
      subroutine tke_budget
c
c -------- get terms in resolved scale tke budget
c          as in gabls writeup at w-points
c          at istage = 1.
c          t_diss and tr_tau are in comp1, pxym is in get_means
c
      use pars
      use fields
      use con_data
      use con_stats
c
      real stat(1:nnz,2)
c
c -------- stat(.,1) = tke transport  = wq
c          stat(.,2) = pressure transport  = wp
c
      do k=1,nnz
         stat(k,1) = 0.0
         stat(k,2) = 0.0
      enddo
c
c --------------- get transport terms as vertical arrays
c
      do k=izs,ize
c
         km1 = k - 1
         do j=iys,iye
         do i=1,nnx
c
c --------- get estimate of turbulent transport term
c
            ufluc   = u(i,j,k) - uxym(k)
            vfluc   = v(i,j,k) - vxym(k)
            wfluc   = w(i,j,k) - wxym(k)
            stat(k,1)  = stat(k,1) + 0.5*wfluc*
     +             (ufluc**2 + vfluc**2 + wfluc**2)
c
c --------- get estimate of pressure transport term
c
            pfluc     = pm(k) + p(i,j,k) - c23*e(i,j,k)
            stat(k,2) = stat(k,2) + pfluc*wfluc
         enddo
         enddo
         stat(k,1) = stat(k,1)*fnxy
         stat(k,2) = stat(k,2)*fnxy
      enddo
      call mpi_sum_z(stat(1,1),i_root,myid,nnz*2,1)
c
c ------ we have all terms on all processors for all z, add them up
c        treat tr_tau at bottom special, tr_tau = 0
c
      tr_tau(0) = 0.0
      do k=1,nnz
c
         kp1 = k + 1
         km1 = k - 1
         if(k .eq. nnz) then
            t_tau(k) = 0.0
            t_wp(k)  = 0.0
            t_wq(k)  = 0.0
         else
            t_tau_u   = 0.5*(tr_tau(kp1) + tr_tau(k))
            t_tau_l   = 0.5*(tr_tau(km1) + tr_tau(k))
            t_tau(k) = -(t_tau_u - t_tau_l)*dzu_i(kp1)
            t_wq(k)  = -(stat(kp1,1) - stat(k,1))*dzu_i(kp1)
            t_wp(k)  = -(stat(kp1,2) - stat(k,2))*dzu_i(kp1)
         endif
         dudz = (uxym(kp1) - uxym(k))*dzu_i(kp1)
         dvdz = (vxym(kp1) - vxym(k))*dzu_i(kp1)
c
c ------------- gather all the budget terms
c
         t_tran(k)  = t_wq(k) + t_wp(k) + t_tau(k)
         t_rprod(k) = -(dudz*uwle(k) + dvdz*vwle(k))
         t_sprod(k) =  (dudz*uwsb(k) + dvdz*vwsb(k))
         t_buoy(k)  =  batag*wtle(k,1)
c
      enddo
c
      return
      end
      subroutine write_his(iloc)
c
c ----- write history file with global parameters
c       write tsfcc specially to preserve digits!
c
      use pars
      use fields
      use con_data
      use con_stats
      use particles
c
      real save_theta(0:nnz)
c
      divgmax = 0.0
      do iz=1,nnz
         divgmax = amax1(divgmax, divz(iz))
      enddo
c
      ziavg = zi
      holtop = hol
      wt_min = wtsb(iloc,1)
      wt_le  = wtle(iloc,1)
      krec = krec + 1
      mid = nnz/4
c
      write(nhis1,6000) time,dt,utau,ziavg,amonin,holtop,
     +         (tsfcc(1)-t_ref),uusfc,uvsfc,uwsfc,vvsfc,vwsfc,wwsfc,
     +         divgmax, wt_min, wt_le,
     +         ucfl, vcfl, wcfl, wtsfc(1),
     +         ups(mid),vps(mid),wps(mid),tps(mid,1),
     +         uwle(mid),uwsb(mid),uw_tot(mid),
     +         vwle(mid),vwsb(mid),vw_tot(mid),
     +         uwcle(mid),vwcle(mid),
     +         wtle(mid,1),wtsb(mid,1),wctle(mid,1),wt_tot(mid,1),
     +         englez(mid),eavg(mid), wabs, float(ip_max),
     +         p_iter, p_err, drg_x, drg_y,
     +         Rep_avg,float(his_reintro)
c    +         tps(mid,2), tps(mid,3),
c    +         wtle(mid,2),wtsb(mid,2),wt_tot(mid,2),
c    +         wtle(mid,3),wtsb(mid,3),wt_tot(mid,3)
 6000 format(46e17.8)
c
c ------ to increase resolution of theta
c        subtract reference before saving history file
c
      do k=1,nnz
         save_theta(k) = txym(k,1)
         txym(k,1)     = txym(k,1) - t_ref
      enddo
c
c -------------- write profile information
c
      call write_prof(nhisp,krec,isize,c_s%wwsb)
c
      do k=1,nnz
         txym(k,1) = save_theta(k)
      enddo
c
      return
      end
      subroutine write_prof(nhisp,krec,num,f)
      real f(num)
      real*4 f32(num)
c
c -------------- build special 32 bit arrays for profiles
c
      do i=1,num
         f32(i) = f(i)
      enddo
c
      write(nhisp,err=999,rec=krec) (f32(i),i=1,num)
c
      return
c --------------- errors
  999 continue
      write(6,9000) num,krec
 9000 format(' 9000, trouble in ',
     +       'SR. save_prof cannot write profile data ',/,
     +       ' num = ',i8, 'krec = ',i6)
      stop
      end
      subroutine close_his
c
c ---- close history files
c
      use pars
      logical there
c
c ---- root closes and checks the files
c
      close(nhis1)
      close(nhisp)
      inquire(file=path_sav_h,exist=there)
      if(.not.there) then
         write(6,8000) path_sav_h
         call mpi_finalize(ierr)
         stop
      endif
      inquire(file=path_sav_hp,exist=there)
      if(.not.there) then
         write(6,8100) path_sav_hp
         call mpi_finalize(ierr)
         stop
      endif
      write(6,7000) path_sav_h
      write(6,7100) path_sav_hp
c
      return
c -------------------- process write errors
 7000 format(' HISTORY DATA IS WRITTEN IN FILE  ',a80)
 7100 format(' PROFILE HISTORY DATA IS WRITTEN IN FILE  ',a80)
 8000 format(' SR. SAVE_HIS: Truble history file not in path =',a80)
 8100 format(' SR. SAVE_HIS: Trouble profile history file',
     +       ' not in path =',a80)
      end
      subroutine dealias_grid_speed
c
c --------- wave cutoff filter using 2d fft
c           for boundary condition values
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real wve(nny,jxs:jxe), gm_work(nnx+2,iys:iye)
c
c --------- sharp spectral cutoff, specific to current 2dfft
c
      ix_cut   = 2*int(float(nnx)/3.) + 3
      iy_cut_l = int(float(nny)/3.) + 2
      iy_cut_u = nnyp2 - iy_cut_l
c
      if(iss == 0) then
c
      do j=iys,iye
      do i=1,nnx
c        gm_work(i,j) = gm(i,j,izs-1,5)
         gm_work(i,j) = gm2(i,j,izs-1,5)
      enddo
      enddo
      call fft2d_mpi(gm_work(1,iys),wve(1,jxs),
     +           trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           1,1,myid,ncpu_s,numprocs,-2)
      call sharp_bndy(wve,ix_cut,iy_cut_l,iy_cut_u)
      call fft2d_mpi(gm_work(1,iys),wve(1,jxs),
     +           trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           1,1,myid,ncpu_s,numprocs,2)
      do j=iys,iye
      do i=1,nnx
c        gm(i,j,izs-1,5)  = gm_work(i,j)
         gm2(i,j,izs-1,5) = gm_work(i,j)
      enddo
      enddo
c
      endif
c
      if(ise == numprocs-1) then
c
      do j=iys,iye
      do i=1,nnx
c        gm_work(i,j) = gm(i,j,ize,5)
         gm_work(i,j) = gm2(i,j,ize,5)
      enddo
      enddo
      call fft2d_mpi(gm_work(1,iys),wve(1,jxs),
     +           trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           1,1,myid,ncpu_s,numprocs,-2)
      call sharp_bndy(wve,ix_cut,iy_cut_l,iy_cut_u)
      call fft2d_mpi(gm_work(1,iys),wve(1,jxs),
     +           trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           1,1,myid,ncpu_s,numprocs,2)
      do j=iys,iye
      do i=1,nnx
c        gm(i,j,ize,5)  = gm_work(i,j)
         gm2(i,j,ize,5) = gm_work(i,j)
      enddo
      enddo
c
      endif
c
      return
      end
      subroutine dealias
c
c --------- wave cutoff filter using 2d fft
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      real wve(nny,jxs:jxe,izs:ize)
      real wves(nnxp2,iys:iye,izs:ize)
c
c --------- sharp spectral cutoff, specific to current 2dfft
c
      ix_cut   = 2*int(float(nnx)/3.) + 3
      iy_cut_l = int(float(nny)/3.) + 2
c
c     ix_cut   = 2*int(float(nnx_c)/3.) + 3
c     iy_cut_l = int(float(nny_c)/3.) + 2
c
      iy_cut_u = nnyp2 - iy_cut_l
c
c ---------- u-equation
c
      call fft2d_mpi(u(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)
      call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
      call fft2d_mpi(u(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,2)
c
c ---------- v-equation
c
      call fft2d_mpi(v(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)
      call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
      call fft2d_mpi(v(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,2)
c
c ---------- w-equation
c
      call fft2d_mpi(w(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)
      call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
      call fft2d_mpi(w(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,2)
c
c ---------- e-equation
c
      call fft2d_mpi(e(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)
      call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
      call fft2d_mpi(e(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,2)
c
c ---------- uf-equation
c
      call fft2d_mpi(uf(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)
      call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
      call fft2d_mpi(uf(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,2)
c
c ---------- vf-equation
c
      call fft2d_mpi(vf(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)
      call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
      call fft2d_mpi(vf(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,2)
c
c ---------- wf-equation
c
      call fft2d_mpi(wf(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)
      call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
      call fft2d_mpi(wf(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +           nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,2)
c
c ------------- scalars, not stored in correct order
c
      do iscl=1,nscl
         do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            wves(ix,iy,iz) = t(ix,iy,iscl,iz)
         enddo
         enddo
         enddo
         call fft2d_mpi(wves(1,iys,izs),wve(1,jxs,izs),trigx(1,1),
     +           trigc,nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,-2)
         call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
         call fft2d_mpi(wves(1,iys,izs),wve(1,jxs,izs),trigx(1,1),
     +           trigc,nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +           izs,ize,myid,ncpu_s,numprocs,2)
         do iz=izs,ize
         do iy=iys,iye
         do ix=1,nnx
            t(ix,iy,iscl,iz) = wves(ix,iy,iz)
         enddo
         enddo
         enddo
      enddo
c
      return
      end
      subroutine sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
c
c --------- sharp cutoff filter for field wve stored
c           in 2d-fft order
c
      use pars
      real wve(nny,jxs:jxe,izs:ize)
c
      do iz=izs,ize
         do ix=jxs,jxe
         do iy=iy_cut_l,iy_cut_u
            wve(iy,ix,iz) = 0.0
         enddo
         enddo
      enddo
c
      if(jxe .lt. ix_cut) go to 999
c
         do iz=izs,ize
            do ix=max(jxs,ix_cut),jxe
            do iy=1,nny
               wve(iy,ix,iz) = 0.0
            enddo
            enddo
         enddo
c
  999 continue
c
      return
      end
      subroutine sharp_bndy(wve,ix_cut,iy_cut_l,iy_cut_u)
c
c --------- sharp cutoff filter for boundary array
c           in 2d-fft order
c
      use pars
      real wve(nny,jxs:jxe)
c
      do ix=jxs,jxe
      do iy=iy_cut_l,iy_cut_u
         wve(iy,ix) = 0.0
      enddo
      enddo
c
      if(jxe .lt. ix_cut) go to 999
c
         do ix=max(jxs,ix_cut),jxe
         do iy=1,nny
            wve(iy,ix) = 0.0
         enddo
         enddo
c
  999 continue
c
      return
      end
      subroutine gridd
c
c ----------- allocate space and pass arrays using modules
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
      use particles
c
c ------------- establish association between pointers
c               and data structures
c
      call fill_cc
      call fill_cs
c
      if (myid==0) write(6,5001) isize
 5001 format(' size of stats array = ',i8)
c
c ---------------- debug for arrays
c
      big = -99.0e+300
c
c ---------------- setup grid
c
      nnx = nxg1
      nny = nyg1
      nnz = nzg1
c     izs = 1
c     ize = nnz
c
c
c ----------- make sure problem and cpu's match
c
      maxp   = numprocs-1
      ncpu_z = numprocs/ncpu_s
      if(mod(numprocs,ncpu_s) .ne. 0 .or.
     +   ncpu_z .gt. nnz) then
         go to 999
      endif
      if(l_root) write(6, 1100) ncpu_s, ncpu_z, numprocs,
     +                          maxp
      write(nprt,1100) ncpu_s, ncpu_z, numprocs, maxp
 1100 format(' Number of x-y slab cpus = ',i5,/,
     +       ' Number of z-level cpus  = ',i5,/,
     +       ' Total number of cpus    = ',i5,/,
     +       ' Max-p for index arrays  = ',i5)
c
c ---------------- allocate arrays for (i,j,k)-indexing on
c                  each processor (see set_range)
c
      allocate(ix_s(0:maxp), ix_e(0:maxp),
     +         jx_s(0:maxp), jx_e(0:maxp),
     +         kx_s(0:maxp), kx_e(0:maxp),
     +         mx_s(0:maxp), mx_e(0:maxp),
     +         iy_s(0:maxp), iy_e(0:maxp),
     +         jy_s(0:maxp), jy_e(0:maxp),
     +         is_s(0:maxp), is_e(0:maxp),
     +         iz_s(0:maxp), iz_e(0:maxp),
     +         i2x_s(0:maxp), i2x_e(0:maxp))
c
c ---------------- setup array sizes and variable dimensions
c
      nxy   = nnx*nny
      ncx   = nnx/2 + 1
      ncy   = nny/2 + 1
      nnxp1 = nnx + 1
      nnyp1 = nny + 1
      nnxp2 = nnx + 2
      nnyp2 = nny + 2
      nnzp1 = nnz + 1
      nnzm1 = nnz - 1
      ivis = ivis0
      fnxy  = 1.0/float(nnx*nny)
c
      write(nprt,7001) nnx,nny,nnz
 7001 format(' 7001 gridd nnx = ',i4,' nny = ',i4,' nnz = ',i4)
c
      call set_range
c
      num_y = iye + 1 - iys
c
c ------------- allocate solution arrays
c               account for nnxp2 for fields but not in rhs
c               and possible monotone for scalars
c
      allocate(u(nnxp2,iys:iye,izs-1:ize+1),
     +         v(nnxp2,iys:iye,izs-1:ize+1),
     +         w(nnxp2,iys:iye,izs-1:ize+1),
     +         t(nnxp2,iys:iye,nscl,izs-2:ize+2),
     +         e(nnxp2,iys:iye,izs-1:ize+1),
     +         r1(nnx,iys:iye,izs-1:ize+1),
     +         r2(nnx,iys:iye,izs-1:ize+1),
     +         r3(nnx,iys:iye,izs-1:ize+1),
     +         r4(nnx,iys:iye,nscl,izs-1:ize+1),
     +         r5(nnx,iys:iye,izs-1:ize+1),
     +         rj(nnx,iys:iye,izs-1:ize+1),
     +         uf(nnxp2,iys:iye,izs-1:ize+1),
     +         vf(nnxp2,iys:iye,izs-1:ize+1),
     +         wf(nnxp2,iys:iye,izs-1:ize+1),
     +         ufi(nnxp2,iys:iye,izs-1:ize+1),
     +         vfi(nnxp2,iys:iye,izs-1:ize+1),
     +         wfi(nnxp2,iys:iye,izs-1:ize+1))
c
c ------------- allocate space for grid and metrics
c
      allocate(xgrid(nnx),
     +         ygrid(nny),
     +         zgrid_w(nnx,iys:iye,izs-1:ize+1,3),
     +         zgrid_u(nnx,iys:iye,izs-1:ize+1),
     +         gm(nnx,iys:iye,izs-1:ize+1,5),
     +         gm2(nnx,iys:iye,izs-1:ize+1,5),
     +         gm3(nnx,iys:iye,izs-1:ize+1,5),
     +         gms(nnx,iys:iye,3),
     +         bndy(nnxp2,iys:iye),
     +         hdot(nnxp2,iys:iye),
     +         u_orbit(nnxp2,iys:iye),
     +         v_orbit(nnxp2,iys:iye),
     +         w_orbit(nnxp2,iys:iye),
     +         bndy_t(nnx,iys:iye,3),
     +         hdot_t(nnx,iys:iye,3),
     +         z1b(nnx,iys:iye),
     +         z1b_2(nnx,iys:iye),
     +         z1b_3(nnx,iys:iye))

c
c ------------- allocate extended arrays for interpolation of
c               particle/spray location
c
      if (ispray==1) then
      allocate(uext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),
     +         vext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),
     +         wext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),
     +         Text(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),
     +         T2ext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3))
      !Transposed velocities to do the uf interpolation:
      allocate(u_t(0:nnz+1,iys:iye,mxs:mxe),
     +         v_t(0:nnz+1,iys:iye,mxs:mxe),
     +         w_t(0:nnz+1,iys:iye,mxs:mxe),
     +         T_t(0:nnz+1,iys:iye,mxs:mxe),
     +         T2_t(0:nnz+1,iys:iye,mxs:mxe))
      end if
      !Keep track of particle counts at each node (its surrounding volume)
      !NOTE: allocate even if ispray == 0, since it's in xy_stats
      allocate(partcount(nnx,iys:iye,izs-1:ize+1))
      allocate(partcount_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(partcountU(nnx,iys:iye,izs-1:ize+1))
      allocate(partcountU_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(partflux(nnx,iys:iye,izs-1:ize+1))
      allocate(partflux_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(partfluxup(nnx,iys:iye,izs-1:ize+1))
      allocate(partfluxup_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(partfluxdn(nnx,iys:iye,izs-1:ize+1))
      allocate(partfluxdn_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(upwp_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(upwp(nnx,iys:iye,izs-1:ize+1))
      allocate(vpsum(nnx,iys:iye,izs-1:ize+1,1:3))
      allocate(vpsum_t(0:nnz+1,iys:iye,mxs:mxe,1:3))
      allocate(vpsqrsum(nnx,iys:iye,izs-1:ize+1,1:3))
      allocate(vpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe,1:3))
      allocate(Tpsum(nnx,iys:iye,izs-1:ize+1))
      allocate(Tpsum_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(Tpsqrsum(nnx,iys:iye,izs-1:ize+1))
      allocate(Tpsqrsum_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(wpTpsum(nnx,iys:iye,izs-1:ize+1))
      allocate(wpTpsum_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(partsrc(nnx,iys:iye,izs-1:ize+1,1:3))
      allocate(partsrc_t(0:nnz+1,iys:iye+1,mxs:mxe+1,1:3))
      allocate(partTsrc(nnx,iys:iye,izs-1:ize+1))
      allocate(partTsrc_t(0:nnz+1,iys:iye+1,mxs:mxe+1))
      allocate(Tfsum(nnx,iys:iye,izs-1:ize+1))
      allocate(Tfsum_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(qfsum(nnx,iys:iye,izs-1:ize+1))
      allocate(qfsum_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(radsum(nnx,iys:iye,izs-1:ize+1))
      allocate(radsum_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(rad2sum(nnx,iys:iye,izs-1:ize+1))
      allocate(rad2sum_t(0:nnz+1,iys:iye,mxs:mxe))
      allocate(qstarsum(nnx,iys:iye,izs-1:ize+1))
      allocate(qstarsum_t(0:nnz+1,iys:iye,mxs:mxe))

c
c ------------- allocate space for boundary condition arrays
c               on top and bottom of domain
c
      allocate(ubc(nnx,iys:iye,2),
     +         vbc(nnx,iys:iye,2),
     +         wbc(nnx,iys:iye,2),
     +         tbc(nnx,iys:iye,nscl,2),
     +         ebc(nnx,iys:iye,2),
     +         pbc(nnx,iys:iye,2),
     +         pbc2(nnx,iys:iye,2),
     +         ufbc(nnx,iys:iye,2),
     +         vfbc(nnx,iys:iye,2),
     +         wfbc(nnx,iys:iye,4))
c
c ------------ allocate space for wind and surface arrays
c
      allocate(wind(nnx,iys:iye),
     +         tau13m(nnx,iys:iye),
     +         tau23m(nnx,iys:iye),
     +         taut3m(nnx,iys:iye,nscl),
     +         t_grnd(nnx,iys:iye,nscl))
c
c ------------------- allocate space for derivative arrays
c
      allocate(ux(nnx,iys:iye,izs-1:ize+1),
     +         uy(nnx,iys:iye,izs-1:ize+1),
     +         vx(nnx,iys:iye,izs-1:ize+1),
     +         vy(nnx,iys:iye,izs-1:ize+1),
     +         wx(nnx,iys:iye,izs-1:ize+1),
     +         wy(nnx,iys:iye,izs-1:ize+1))
c
c ------------- allocate space for pressure, pressure bcs
c
      allocate(p(nnxp2,iys:iye,izs-1:ize+1),
     +         p_old(nnxp2,iys:iye,izs-1:ize+1),
     +         ptop(nnxp2,iys:iye,2))
c
c ------------- allocate space for viscosity and diffusivity
c
      allocate(vis_m(nnx,iys:iye,izs-1:ize+1),
     +         vis_s(nnx,iys:iye,izs-1:ize+1))
c
c ------------- allocate space for tau bc
c
      allocate(t11bc(nnx,iys:iye),
     +         t12bc(nnx,iys:iye),
     +         t13bc(nnx,iys:iye),
     +         t22bc(nnx,iys:iye),
     +         t23bc(nnx,iys:iye),
     +         t33bc(nnx,iys:iye),
     +         ustar(nnx,iys:iye),
     +         amol(nnx,iys:iye),
     +         t1tbc(nnx,iys:iye,nscl),
     +         t2tbc(nnx,iys:iye,nscl),
     +         t3tbc(nnx,iys:iye,nscl))
c
c ------------- allocate space for fft trig factors
c
      nq_trig = max(nnx,nny)
      allocate(trigx(2*nq_trig+15,2),
     +         trigc(4*nq_trig+15), trigcx(4*nq_trig+15))
      return
  999 continue
c
      if(l_root) write(6,1000) numprocs, ncpu_s, mmz
      write(nprt,1000) numprocs, ncpu_s, nnz
 1000 format(' Gridd Trouble number of processors and grid',
     +          ' partitioning do not match!',/,
     +          ' Total num of cpus   = ',i5,
     +          ' Num cpu on x-y slab = ',i5,/,
     +          ' Num of z-levels     = ',i5)
      call mpi_finalize(ierr)
      end
      subroutine restart
c
c ----------- get restart file from local directory
c
      use pars
      use fields
      use con_data
      use con_stats
      character*80 path_res_c
      logical there
c
c --------------------- check if file is there
c
      inquire(file=path_res,exist=there)
      if(there) then
         if(l_root) write(6,6001) path_res
      else
         if(l_root) write(6,6005) path_res
         stop
      endif
c
c ------------------ get constant file
c
      iloc = index(path_res,' ')
      path_res_c = path_res(1:iloc-1)//'.con'
      inquire(file=path_res_c,exist=there)
      if(there) then
         if(l_root) write(6,6002) path_res_c
      else
         if(l_root) write(6,6006) path_res_c
         stop
      endif
      open(nvelc,err=200,file=path_res_c,form='unformatted',
     +        status='old')
c
      call read_res
c
      return
c ---------------------------- process errors
  100 continue
      write(6,9000) path_res, nvel
      call mpi_finalize(ierr)
      stop
c -----------------------
  200 continue
      write(6,9001) path_res_c, nvelc
      call mpi_finalize(ierr)
      stop
c -----------------------
 6001 format(' SR. RESTART: FILE READ = ',A80)
 6002 format(' SR. RESTART: CONSTANT FILE READ = ',A80)
 6005 format(' 6005, SR. RESTART: cannot find restart file = ',a80)
 6006 format(' 6005, SR. RESTART: cannot find constant file = ',a80)
 9000 format(' 9000, SR. RESTART: cannot open file =',a80,/,
     +       ' to unit number = ',i2)
 9001 format(' 9001, SR. RESTART: cannot open file =',a80,/,
     +       ' to unit number = ',i2)
      end
      subroutine read_res
c
c -------------- read restart file including constant file
c                changed for iys:iye
c
      use pars
      use fields
      use con_data
      use con_stats
#if defined(SWAP)
      use module_byteswap
#endif
      include 'mpif.h'
c
      integer status(mpi_status_size), ierr
      integer(kind=mpi_offset_kind) :: offset, disp
      integer(kind=k8)              :: nsize, nsize2
      real, allocatable, dimension(:,:,:) :: temp
      allocate(temp(nvar+4,nnx,iys:iye))
c
c ---- open file
c
      call mpi_file_open(mpi_comm_world, path_res,
     +                   mpi_mode_create+mpi_mode_rdwr,
     +                   mpi_info_null, nvel, ierr)
c
c ---- set file view
c
      disp = 0
      call mpi_file_set_view(nvel,disp,mpi_real8,mpi_real8,
     +                      'native',mpi_info_null,ierr)
c
c ------------ read 3d fields
c
      nsize  = int(nvar+4,k8)*nnx*nny
      nsize2 = int(nvar+4,k8)*nnx*(iys-1)
      n_read = (nvar+4)*nnx*(iye+1-iys)
c
      do k=izs,ize
         offset = int((k-1),k8)*nsize + nsize2
         call mpi_file_read_at_all(nvel,offset,temp,n_read,
     +                              mpi_real8,status,ierr)
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
         enddo
         enddo
         do is = 1,nscl
            do j = iys,iye
            do i = 1,nnx
               t(i,j,is,k) = temp(6+is,i,j)
            enddo
            enddo
         enddo
c
      enddo
c
c ---- close file
c
      call mpi_file_close(nvel, ierr)
c
      deallocate(temp)
c
c ------------ every mpi process reads constant file
c
      rewind(nvelc)
      read(nvelc,err=9993) c_c, c_s, case
      close(nvelc)
c
      if(l_root) write(6,4001) case
 4001 format(' 4001, SR. RESTART: case from restart = ',a3)
c
c ----- special restart conditions -------------------------------------
c
c -------- set case name to case input
c
      case   = case_inp
      if(l_root) write(6,4002) case_inp, utau, utausv
 4002 format(' 4002, SR. RESTART:',/,
     +       ' files will be saved with case name = ',a3,/,
     +       ' utau = ',e15.6,' utausv = ',e15.6)
c
c ------------------- if new vis model set match point for
c                     outer grid
      nmatch = 48
      utau = utausv
c
c -------- hand coded changes to restart if needed
c
        qstars = 0.000
        wtsfcs = 0.000
        zos    = 250.0e-06
c
c
c ------ reset qstar and wtsfc for no heat flux
c
c              qstar(1) = qstars
c              wtsfc(1) = wtsfcs
c              qstar(2) = qstars
c              wtsfc(2) = wtsfcs
c              zo       = zos
c ------ redefine case id to input value
c              case = cases
c
c
c --------- flaw in t_zero ... fix in later runs
c
      if(iti .lt. i_move_step) then
         t_zero = time
      endif
c
      if(l_root) write(6,4012) time, t_zero, dt, dt1, dt_new
      if(l_root) write(6,4013) qstar(1) , nmatch, case, zo
c
      call get_dz
c
      return
c ------------------------  process errors from read
c9991 continue
c     write(6,6000) nvel,iz
c6000 format(' SR. READ_RES: hit end of file on unit number = ',i2,/,
c    +       '               at iz = ',i4)
c     call mpi_finalize(ierr)
c     stop
c ---------------------
 9992 continue
      write(6,6100) nvel,iz
 6100 format(' SR. READ_RES: error reading file on unit number = ',i2,/,
     +       '               at iz = ',i4)
      call mpi_finalize(ierr)
      stop
c ---------------------
 9993 continue
      write(6,6200) nvelc
 6200 format(' SR. READ_RES:',/,
     +       '    error reading constant file on unit number = ',i2)
      call mpi_finalize(ierr)
      stop
c ---------------------
 4012 format(' SR. RESTART: ',/,
     +       ' T = ',e15.6,5x,' T_zero = ',e15.6,/,
     +       ' dt = ',e15.6,' dt1 = ',e15.6,' dt_new = ',e15.6)
 4013 format(' after restart qstar = ',e15.6,' nmatch = ',i5,
     +       ' case = ',a3,' zo = ',e15.6)
      end
      subroutine save_viz(it)
c
c --------------- save multiple (x-y), (x-z), (y-z), planes of data .
c                 modify recl in all open statements for more or less
c                 variables.
c                 Constant - x, implies yz planes
c                 Constant - y, implies xz planes
c                 Constant - z, implies xy planes
c
c ------------- routine uses send/recv to get information in y-z planes
c
      use pars
      use fields
      use con_data
      use con_stats
      use fftwk
#if defined(SWAP)
      use module_byteswap
#endif
      include 'mpif.h'
c
c ------ for xz and yz planes output z grid also
c
      parameter(nvar_o = 6, nvar_o1 = nvar_o + 1)
c
      integer ix_pick(maxnx),  iy_pick(maxny),  iz_pick(maxnz),
     +        ix_order(maxnx), iy_order(maxny), iz_order(maxnz)
c
      integer istatus(mpi_status_size), ierr
      integer(kind=mpi_offset_kind) :: offset, disp
      integer(kind=k8)              :: nsize, nsize2
c
      real(kind=4), dimension(nvar_o1,nny,izs:ize) :: temp_x
      real(kind=4), dimension(nvar_o1,nnx,izs:ize) :: temp_y
      real(kind=4), dimension(nvar_o,nnx,iys:iye)  :: temp_z
      real, dimension(nvar_o1,iys:iye,izs:ize)     :: buf_send
      real(kind=8), dimension(nnx,0:nnz)           :: z_xz
      real(kind=8), dimension(nny,0:nnz)           :: z_yz
c
c ------------- don't touch
c
      data iviz_x,  iviz_y,  iviz_z  /0, 0, 0/
      data ionce_x, ionce_y, ionce_z, istuff /0, 0, 0, 0/
      data ix_pick, iy_pick, iz_pick /maxnx*0, maxny*0, maxnz*0/
      data ix_order, iy_order, iz_order /maxnx*0, maxny*0, maxnz*0/
      save iviz_x,  iviz_y,  iviz_z,
     +     ix_pick, iy_pick, iz_pick,
     +     ix_order, iy_order, iz_order,
     +     ionce_x, ionce_y, ionce_z, istuff,
     +     npln_x, npln_y, npln_z
c
c
c ----------- turn on z levels to save. Customize for your own use.
c             Set iz_pick(iz) = iz, ix_pick(ix) = ix, iy_pick(iy) = iy
c             Data is round-robin alternated in the data file for more than
c             1 plane for any particular view.
c
      iz_pick(3)   = 3
      iz_pick(6)   = 6
      iz_pick(10)  = 10
      iz_pick(20)  = 20
      iz_pick(40)  = 40
      iz_pick(85)  = 85
      iz_pick(100) = 100
c
c -------------- pick an x-z plane of data (can add more)
c
      iy_pick(1)         = 1
      iy_pick(nny/4)     = nny/4
      iy_pick(nny/2)     = nny/2
      iy_pick((3*nny)/4) = (3*nny)/4
c
c -------------- pick a y-z plane of data (can add more)
c
      ix_pick(nnx/4) = nnx/4
      ix_pick(nnx/2) = nnx/2
      ix_pick(nnx)   = nnx
c
c ------ find total number of z's turned on and open file once
c
      if(ionce_z .eq. 0) then
         npln_z = 0
         do k=1,nnz
            if(iz_pick(k) .eq. k) then
               npln_z = npln_z + 1
               iz_order(k) = npln_z
            endif
         enddo
         ionce_z = 1
         iviz_z =  -npln_z
         if(npln_z .ne. 0) then
            call mpi_file_open(mpi_comm_world, path_viz_xy,
     +                         mpi_mode_create+mpi_mode_rdwr,
     +                         mpi_info_null, nviz_z, ierr)
            disp = 0
            call mpi_file_set_view(nviz_z,disp,mpi_real4,mpi_real4,
     +                            'native',mpi_info_null,ierr)
         endif
      endif
c
c ------ find total number of y's turned on and open file once
c
      if (ionce_y .eq. 0) then
         npln_y = 0
         do j=1,nny
            if(iy_pick(j) .eq. j) then
               npln_y = npln_y + 1
               iy_order(j) = npln_y
            endif
         enddo
         ionce_y = 1
         iviz_y  = -npln_y
         if(npln_y .ne. 0) then
            call mpi_file_open(mpi_comm_world, path_viz_xz,
     +                         mpi_mode_create+mpi_mode_rdwr,
     +                         mpi_info_null, nviz_y, ierr)
            disp = 0
            call mpi_file_set_view(nviz_y,disp,mpi_real4,mpi_real4,
     +                            'native',mpi_info_null,ierr)
         endif
      endif
c
c ------ find total number of x's turned on and open file once
c
      if(ionce_x .eq. 0) then
         npln_x = 0
         do i=1,nnx
            if(ix_pick(i) .eq. i) then
               npln_x = npln_x + 1
               ix_order(i) = npln_x
            endif
         enddo
         ionce_x = 1
         iviz_x  = -npln_x
         if(npln_x .ne. 0) then
            call mpi_file_open(mpi_comm_world, path_viz_yz,
     +                         mpi_mode_create+mpi_mode_rdwr,
     +                         mpi_info_null, nviz_x, ierr)
            disp = 0
            call mpi_file_set_view(nviz_x,disp,mpi_real4,mpi_real4,
     +                            'native',mpi_info_null,ierr)
         endif
      endif
c
      if(istuff .eq. 0 .and. l_root) then
         open(nviz_s,file=path_stuf)
         istuff = 1
      endif
c
c --------- write data, subtract t_ref to increase
c           resolution on 32 bit machines
c
c ---------- xy planes of data
c
      iviz_z  = iviz_z + npln_z
      nsize   = int(nvar_o,k8)*nnx*nny
      nsize2  = int(nvar_o,k8)*nnx*(iys-1)
      n_write = nvar_o*nnx*(iye+1-iys)
      do k=izs,ize
         if(iz_pick(k) .eq. k) then
            km1 = k - 1
            do j=iys,iye
            do i=1,nnx
               temp_z(1,i,j) = u(i,j,k)
               temp_z(2,i,j) = v(i,j,k)
               temp_z(3,i,j) = w(i,j,k)
               temp_z(4,i,j) = (t(i,j,1,k) - t_ref)
c
c ---------- get the pressure field assumes flux form for advection
c
               temp_z(5,i,j) = p(i,j,k) - c23*e(i,j,k) - pxym(k)
c
c ---------- get the flux velocity minus the grid speed
c
               temp_z(6,i,j) = wf(i,j,k) - gm(i,j,k,5)
            enddo
            enddo
#if defined(SWAP)
            call byteswap(temp_z)
#endif
            offset = int((iviz_z + iz_order(k) - 1),k8)*nsize + nsize2
            call mpi_file_write_at(nviz_z,offset,temp_z,n_write,
     +                             mpi_real4,istatus,ierr)
            if (ierr .ne. 0) go to 9991
         endif
      enddo
c
c ---------- xz planes of data
c
      iviz_y = iviz_y + npln_y
      nsize  = int(nvar_o1,k8)*nnx*nnz
      nsize2 = int(nvar_o1,k8)*nnx*(izs-1)
      nwrite = (nvar_o1)*nnx*(ize+1-izs)
      do j=iys,iye
         if(iy_pick(j) .eq. j) then
            do k=izs,ize
            km1 = k - 1
            do i=1,nnx
               temp_y(1,i,k) = u(i,j,k)
               temp_y(2,i,k) = v(i,j,k)
               temp_y(3,i,k) = w(i,j,k)
               temp_y(4,i,k) = (t(i,j,1,k) - t_ref)
c
c ---------- get the pressure field
c
               temp_y(5,i,k) =  p(i,j,k) - pxym(k) - c23*e(i,j,k)
c
c ---------- get the flux velocity minus the grid speed
c
               temp_y(6,i,k) = wf(i,j,k) - gm(i,j,k,5)
c
c ---------- the shape of the gridlines
c
               temp_y(7,i,k) =  0.5*(zgrid_w(i,j,km1,1) +
     +                               zgrid_w(i,j,k,1))
            enddo
            enddo
#if defined(SWAP)
            call byteswap(temp_y)
#endif
            offset = int((iviz_y + iy_order(j) - 1),k8)*nsize + nsize2
            call mpi_file_write_at(nviz_y,offset,temp_y,nwrite,
     +                                mpi_real4,istatus,ierr)
            if (ierr .ne. 0) goto 9992
         endif
      enddo
c
c ---------- yz planes that cut across all processors
c            just have root node on that slab write data
c
      iviz_x  = iviz_x + npln_x
      n_write = nvar_o1*nny*(ize+1-izs)
      nsize   = int(nvar_o1,k8)*nny*nnz
      nsize2  = int(nvar_o1,k8)*nny*(izs-1)
      n_send  = nvar_o1*(ize+1-izs)*(iye+1-iys)
      do i=1,nnx
         if(ix_pick(i) .eq. i) then
c
c ----------- build send buffer
c
            do k=izs,ize
            km1 = k - 1
            do j=iys,iye
               buf_send(1,j,k) = u(i,j,k)
               buf_send(2,j,k) = v(i,j,k)
               buf_send(3,j,k) = w(i,j,k)
               buf_send(4,j,k) = (t(i,j,1,k) - t_ref)
c
c ---------- get the pressure field
c
               buf_send(5,j,k) = p(i,j,k) - pxym(k) - c23*e(i,j,k)
c
c ---------- get the flux velocity minus the grid speed
c
               buf_send(6,j,k) = wf(i,j,k) - gm(i,j,k,5)
c
c ---------- shape of the gridlines
c
               buf_send(7,j,k) = 0.5*(zgrid_w(i,j,km1,1) +
     +                                zgrid_w(i,j,k,1))
            enddo
            enddo
            if(myid .ne. iss) then
              call mpi_send(buf_send(1,iys,izs),n_send,
     +                      mpi_real8,iss,1,
     +                      mpi_comm_world,ierr)
            else
              do k=izs,ize
              do j=iys,iye
              do ii=1,nvar_o1
                 temp_x(ii,j,k) = buf_send(ii,j,k)
              enddo
              enddo
              enddo
              do l=iss+1,ise
                 call recv_yz_var(temp_x,nvar_o1,nny,
     +                            iy_s(l),iy_e(l),izs,ize,l)
              enddo
#if defined(SWAP)
              call byteswap(temp_x)
#endif
              offset = int((iviz_x + ix_order(i) - 1),k8)*nsize + nsize2
              call mpi_file_write_at(nviz_x,offset,temp_x,n_write,
     +                          mpi_real4,istatus,ierr)
              if (ierr .ne. 0) goto 9993
            endif
         endif
      enddo
c
c ------------- ascii file with facts in it that goes
c               with visualization
c
      if(l_root) then
         write(nviz_s,5000) time, amonin, zi, utau
 5000    format(4e20.8)
      endif
c
c ---- last time step or hit mtape close the files
c
      if (it .eq. itmax .or. mtape) then
         if(npln_z .ne. 0) then
            call mpi_file_close(nviz_z, ierr)
            ionce_z = 0
         endif
         if(npln_y .ne. 0) then
            call mpi_file_close(nviz_y, ierr)
            ionce_y = 0
         endif
         if(npln_x .ne. 0) then
            call mpi_file_close(nviz_x, ierr)
            ionce_x = 0
         endif
         if(l_root) then
            close(nviz_s)
            istuff = 0
         endif
      endif
c
      return
c --------------------------  errors in writing viz file
 9991 continue
      write(6,6000) nviz_z, iz
 6000 format(' SR. SAVE_VIS:',/,
     +       '    trouble cannot write xy viz file on unit = ',i2,/,
     +       '             at iz = ',i4)
      call mpi_finalize(ierr)
      stop
c --------------------------  errors in writing viz file
 9992 continue
      write(6,6100) nviz_y, iz, iviz_y
 6100 format(' SR. SAVE_VIS:',/,
     +       '    trouble cannot write xz viz file on unit = ',i2,/,
     +       '             at iz = ',i4,/,
     +       '            iviz_y = ',i8)
c --------------------------  errors in writing viz file
 9993 continue
      write(6,6200) nviz_x, iz, iviz_x
 6200 format(' SR. SAVE_VIS:',/,
     +       '    trouble cannot write yz viz file on unit = ',i2,/,
     +       '             at iz = ',i4,/,
     +       '            iviz_x = ',i8)
      call mpi_finalize(ierr)
      stop
      end
      subroutine save_bndy(it)
c
c --------------- save boundary shape and boundary speed
c                 as function of time
c                 at same time step as viz files
c
      use pars
      use fields
      use con_data
      use con_stats
      use fftwk
#if defined(SWAP)
      use module_byteswap
#endif
      include 'mpif.h'
c
      parameter(nvar_o = 2)
c
      integer istatus(mpi_status_size), ierr
      integer(kind=mpi_offset_kind) :: offset, disp
      integer(kind=k8)              :: nsize, nsize2
c
      real(kind=4), dimension(2,nnx,iys:iye) :: temp_z
c
c ------------- don't touch
c
      data ionce_z, iviz_z /0, 0/
      save ionce_z, iviz_z
c
c ------ open file once
c
c     write(nprt,4341) nbndy, path_bndy
 4341 format(' 4341 nbndy = ',i5,' path_bndy = ',a80)
c
      if(ionce_z .eq. 0) then
        ionce_z = 1
        call mpi_file_open(mpi_comm_world, path_bndy,
     +                     mpi_mode_create+mpi_mode_rdwr,
     +                     mpi_info_null, nbndy, ierr)
        disp = 0
        call mpi_file_set_view(nbndy,disp,mpi_real4,mpi_real4,
     +                        'native',mpi_info_null,ierr)
      endif
c
c ---------- xy plane of data
c
      iviz_z  = iviz_z + 1
      nsize   = int(nvar_o,k8)*nnx*nny
      nsize2  = int(nvar_o,k8)*nnx*(iys-1)
      n_write = nvar_o*nnx*(iye+1-iys)
c
      if(iss .ne. 0) go to 999
c
      k = 0
c
      do j=iys,iye
      do i=1,nnx
         temp_z(1,i,j) = bndy_t(i,j,1)
         temp_z(2,i,j) = gm(i,j,k,5)
      enddo
      enddo
#if defined(SWAP)
      call byteswap(temp_z)
#endif
      offset = int((iviz_z - 1),k8)*nsize + nsize2
      call mpi_file_write_at(nbndy,offset,temp_z,n_write,
     +                       mpi_real4,istatus,ierr)
      if (ierr .ne. 0) go to 9991
c
  999 continue
c
c ---- last time step or hit mtape close the files
c
      if(it .eq. itmax .or. mtape) then
         call mpi_file_close(nbndy, ierr)
         ionce_z = 0
         iviz_z  = 0
      endif
c
c
      return
c --------------------------  errors in writing bndy file
 9991 continue
      write(6,6000) nbndy
 6000 format(' SR. SAVE_BNDY:',/,
     +       '    trouble cannot write boundary file on unit = ',i2)
      call mpi_finalize(ierr)
      stop
      end
      subroutine recv_yz_var(temp_x,nvar,nny,iys,iye,izs,ize,ir)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      real buf(nvar,iys:iye,izs:ize)
      real(kind=4), dimension(nvar,nny,izs:ize) :: temp_x
c
      num = nvar*(ize+1-izs)*(iye+1-iys)
      call mpi_recv(buf(1,iys,izs),num,mpi_real8,ir,1,
     +             mpi_comm_world,istatus,ierr)
      do k=izs,ize
      do j=iys,iye
      do ii=1,nvar
         temp_x(ii,j,k) = buf(ii,j,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine save_v(it)
c
c --------------- save 3d fields
c
      use pars
      use fields
#if defined(SWAP)
      use module_byteswap
#endif
      include 'mpif.h'
      logical there
c
      integer status(mpi_status_size), ierr
      integer(kind=mpi_offset_kind) :: offset, disp
      integer(kind=k8)                 nsize, nsize2
c
      real, allocatable, dimension(:,:,:) :: temp
      allocate(temp(nvar+4,nnx,iys:iye))
c
c ---- open file
c
      call mpi_file_open(mpi_comm_world, path_sav_v,
     +                   mpi_mode_create+mpi_mode_rdwr,
     +                   mpi_info_null, nvel, ierr)
c
c ---- set file view
c
      disp = 0
      call mpi_file_set_view(nvel,disp,mpi_real8,mpi_real8,
     +                      'native',mpi_info_null,ierr)
c
c ---- write data
c
      nsize   = int(nvar+4,k8)*nnx*nny
      nsize2  = int(nvar+4,k8)*nnx*(iys-1)
      n_write = (nvar+4)*nnx*(iye+1-iys)
c
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
         enddo
         enddo
         do is = 1,nscl
            do j = iys,iye
            do i = 1,nnx
               temp(6+is,i,j) = t(i,j,is,k)
            enddo
            enddo
         enddo
#if defined(SWAP)
      call byteswap(temp)
#endif
c
         offset = int((k-1),k8)*nsize + nsize2
c        call mpi_file_write_at_all(nvel,offset,temp,n_write,
         call mpi_file_write_at(nvel,offset,temp,n_write,
     +                              mpi_real8,status,ierr)
         if (ierr /= 0) goto 9991
c
      enddo
c
c ---- close file
c
      call mpi_file_close(nvel, ierr)
c
c ---- check file
c
      if (l_root) then
         inquire(file=path_sav_v,exist=there)
         if(.not.there) then
            write(6,8000) nvel,myid
            call mpi_finalize(ierr)
            stop
         endif
         write(6,7000) it,path_sav_v
      endif
c
      deallocate(temp)
c
      return
c --------------------------  errors in writing restart file
 9991 continue
      write(6,6000) nvel, iz
 6000 format(' SR. SAVE_V:',/,
     +       '    trouble cannot write restart file on unit = ',i2,/,
     +       '             at iz = ',i4)
      call mpi_finalize(ierr)
      stop
c --------------------
 7000 format(' **** DATA SET AT IT = ',I6,/,
     +       '      VELOCITY DATA IS WRITTEN IN FILE  ',a80)
 8000 format(' in SAVE_V: trouble writing file ',i5,'  myid = ',i5,
     +       ' at iz = ',i5)
      end
      subroutine save_c(it)
c
c --------------- root process writes constant file
c                 sequential fortan binary
c
      use pars
      use fields
      use fftwk
      use con_data
      use con_stats
c
      logical there
      character options*8, passwd*1
c
c ---- open file
c
      open(nvelc,err=9992,file=path_sav_c,form='unformatted',
     +                status='unknown')
      write(nvelc,err=9992) c_c, c_s, case
      close(nvelc)
c
        inquire(file=path_sav_c,exist=there)
        if(.not.there) then
           write(6,8001) path_sav_c
           call mpi_finalize(ierr)
           stop
        endif
c -----------------------------  output ok message
      write(6,7001) path_sav_c
c
      return
c --------------------------  errors in writing constant file
 9992 continue
      write(6,6100) nvelc
 6100 format(' SR. SAVE_V:',/,
     +  '    trouble cannot open/write constant file on unit = ',i2)
      call mpi_finalize(ierr)
      stop
c ---------------------
 7001 format('      CONSTANT DATA IS WRITTEN IN FILE  ',a80)
 8001 format(' SR. SAVE_C: Trouble constant file not in path =',a80)
      end
      subroutine save_p
c
c -------------- save pressure file
c
      use pars
      use fields
#if defined(SWAP)
      use module_byteswap
#endif
      include 'mpif.h'
      logical there
c
      integer status(mpi_status_size), ierr
      integer(kind=mpi_offset_kind) :: offset, disp
      integer(kind=k8)              :: nsize, nsize2
c
      real temp(nnx,iys:iye)
c
c ---- open file
c
      call mpi_file_open(mpi_comm_world, path_sav_p,
     +                   mpi_mode_create+mpi_mode_rdwr,
     +                   mpi_info_null, npre, ierr)
c
c ---- set file view
c
      disp = 0
      call mpi_file_set_view(npre,disp,mpi_real8,mpi_real8,
     +                      'native',mpi_info_null,ierr)
c
c ---- write data
c
      nsize   = int(nnx,k8)*nny
      nsize2  = int(nnx,k8)*(iys -1)
      n_write = nnx*(iye+1-iys)
      do k=izs,ize
         do j=iys,iye
         do i=1,nnx
            temp(i,j) = p(i,j,k)
         enddo
         enddo
#if defined(SWAP)
      call byteswap(temp)
#endif
         offset = int((k-1),k8)*nsize + nsize2
         call mpi_file_write_at(npre,offset,temp,n_write,
     +                              mpi_real8,status,ierr)
      enddo
c
c ---- close file
c
      call mpi_file_close(npre, ierr)
c
c ---- check file
c
      if (l_root) then
         inquire(file=path_sav_p,exist=there)
         if(.not.there) then
            write(6,8000) path_sav_p
            call mpi_finalize(ierr)
            stop
         endif
         write(6,7000) path_sav_p
      endif
c
      return
c -------------------- process write errors
 9991 continue
      write(6,6000) npre, iz
 6000 format(' SR. SAVE_P:',/,
     +       '    trouble cannot write pressure file on unit = ',i2,/,
     +       '             at iz = ',i4)
      call mpi_finalize(ierr)
      stop
c -----------------------
 7000 format('      PRESSURE DATA IS WRITTEN IN FILE  ',a80)
 8000 format(' SR. SAVE_P: Trouble pressure file not in path =',a80)
      end
      subroutine get_units
c
      use pars
c
c -------------- unit numbers for files
c
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
c
c ------------- unit number for standard print out
c               for each mpi task
c
      nprt = 1
c
c ------------- open unit for standard printout
c
!      path_prt = case_inp(1:3)//'.le.xxxxx.out'
!      write(path_prt(8:12),'(i5.5)') myid
!      open(nprt,file=path_prt,form='formatted')
c
      return
      end
      subroutine get_output_filenames
c
c ----------- build file names for velocity, pressure, and constants
c
      use pars
      include 'mpif.h'
      character cgrid*10, num*3
c
c --------------- build character strings for file name
c
      cgrid = '.le.'
      write(num,'(i3.3)') itn
      call blnk(path_sav_v)
      call blnk(path_sav_p)
      call blnk(path_sav_c)
      call blnk(path_sav_part)
      path_sav_v = trim(path_sav)//'/u'//
     +                 cgrid(1:4)//case(1:3)//num(1:3)
      path_sav_p = trim(path_sav)//'/p'//
     +                 cgrid(1:4)//case(1:3)//num(1:3)
      path_sav_c = trim(path_sav)//'/u'//
     +                 cgrid(1:4)//case(1:3)//num(1:3)//'.con'
      path_sav_part = trim(path_sav)//'/part'//
     +                 cgrid(1:4)//case(1:3)//num(1:3)
c     path_sav_g = trim(path_sav)//'/'//case(1:3)//'.grid'
c
      return
      end
      subroutine open_his(istep)
c
c ------------------- open history files by root
c                     isize determined in sr. fill_cs
c
      use pars
      include 'mpif.h'
      character cgrid*4, iblks*16
      logical there
c
c --------------- build character strings for ascii history file name
c
      cgrid = '.le.'
      call blnk(iblks)
      write(iblks(1:7),'(i7.7)') istep
      iblks(8:8) = '_'
      write(iblks(9:15),'(i7.7)') (istep + itape)
      iblnk = index(path_his,' ')
      call blnk(path_sav_h)
      path_sav_h = path_his(1:iblnk-1)//'/his'//
     +         cgrid(1:4)//case(1:3)//'.'//iblks(1:15)
c
c --------------- build character strings for ieee profile history file
c                 set record counter for direct access file = 0
c
      krec = 0
      cgrid = '.le.'
      call blnk(iblks)
      write(iblks(1:7),'(i7.7)') istep
      iblks(8:8) = '_'
      write(iblks(9:15),'(i7.7)') (istep + itape)
      iblnk = index(path_his,' ')
      call blnk(path_sav_hp)
      path_sav_hp = path_his(1:iblnk-1)//'/his'//
     +         cgrid(1:4)//case(1:3)//'.'//iblks(1:15)//'.ieee'
c
c ----------------- save data in directory
c
      if(l_root) then

      close(nhis1)
      open(nhis1,err=3000,file=path_sav_h,form='formatted')
c
      close(nhisp)
      open(nhisp,err=4000,file=path_sav_hp,
     +        form='unformatted',access='direct',recl=isize*j_recl,
     +        status='unknown')
      endif
c
      return
c ------------------- process errors
 3000 continue
      write(6,6301) nhis1, path_sav_h
 6301 format(' 6301, SR. OPEN_HIS:',/,
     +       '    cannot open history1 file on unit = ',i2,/,
     +       '    path = ',a80)
      stop
c-------------------
 4000 continue
      write(6,6302) nhisp, path_sav_hp
 6302 format(' 6302, SR. OPEN_HIS:',/,
     +       '    cannot open history profile file on unit = ',i2,/,
     +       '    path = ',a80)
      stop
      end
      subroutine range(n1,n2,nprocs,irank,ista,iend)
c
c ---------- the ibm range finder to balance load
c
      iwork1 = (n2 - n1 + 1)/nprocs
      iwork2 = mod(n2 - n1 +1, nprocs)
      ista = irank*iwork1 + n1 + min(irank,iwork2)
      iend = ista + iwork1 - 1
      if(iwork2 .gt. irank) iend = iend + 1
c
      return
      end
      subroutine viz_output_filename(istep)
c
c ------------------- set visualization files,
c                     leaves files in scratch directory
c
      use pars
      include 'mpif.h'
      character iblks*16
c
c --------------- build character strings for file names
c                 with time step
c
      call blnk(iblks)
      iblks(1:1) = '.'
      write(iblks(2:8),'(i7.7)') istep
      iblks(9:9) = '_'
      write(iblks(10:16),'(i7.7)') (istep + itape)
c
      iloc = index(path_seed,' ')
      path_viz_xy = path_seed(1:iloc-1)
     +         //'/viz.'//case(1:3)//iblks(1:16)//'.xy.data'
c
      path_viz_xz = path_seed(1:iloc-1)
     +         //'/viz.'//case(1:3)//iblks(1:16)//'.xz.data'
c
      path_viz_yz = path_seed(1:iloc-1)
     +         //'/viz.'//case(1:3)//iblks(1:16)//'.yz.data'
c
      path_stuf = path_seed(1:iloc-1)
     +         //'/stuff.'//case(1:3)//iblks(1:16)//'.data'
c
      path_bndy = path_seed(1:iloc-1)
     +         //'/bndy.'//case(1:3)//iblks(1:16)//'.data'
c
c     if(l_root) then
c        write(6,8001) path_viz_xy
c8001    format(' 8001: viz file = ',a80)
c        write(6,8001) path_viz_xz
c        write(6,8001) path_viz_yz
c        write(6,8001) path_stuf
c        write(6,8001) path_seed
c     endif
c
      return
      end
      subroutine set_range
c
c ---- build special x,y,z-ranges. dimensioned for 0:numprocs-1
c      indexed with myid
c
c      [ix_s:ix_e] x-range for computing y-derivatives nx-pts/ncpu_s
c                  in xtoy and ytox tranposes
c
c      [i2x_s:2ix_e] x-range for computing complx y-derivatives 2*nx-pts/ncpu_s
c                  in xtoy and ytox tranposes. special for moving surfaces
c
c      [jx_s:jx_e] x-range for computing 2d fft (nx+2)-pts/ncpu_s
c                  must be even in each x-interval for complex fft in y
c
c      [kx_s:kx_e] x-range for pressure solver transpose (nx+2)-pts/ncpu_z
c                  nx+2 fourier coefficients for xtoz and ztox transposes
c
c      [mx_s:mx_e] x-range split across z cpus as nx-pts/ncpu_z
c                  for use in surface layer routines
c
c      [is_s:is_e] starting and ending processor id's for a
c                  particular z-level
c
c      [iy_s:iy_e] y-range for computing y-derivatives ny-pts/ncpu_s
c                  in xtoy and ytox tranposes
c
c      [jy_s:jy_e] y-range for use in xtoz and ztox transposes
c                  in pressure solution
c
c      [iz_s:iz_e] z-range for a particular vertical slab
c
c
      use pars
c
      ii = -1
      do nn=0,ncpu_z-1
         call range(1,nnx+2,ncpu_z,nn,lx_s,lx_e)
         call range(1,nnx,ncpu_z,nn,nx_s,nx_e)
         call range(1,nny,ncpu_z,nn,ly_s,ly_e)
         call range(1,nnz,ncpu_z,nn,mz_s,mz_e)
         do mm=0,ncpu_s-1
            call range(1,nny,ncpu_s,mm,ny_s,ny_e)
            call range(1,nnx,ncpu_s,mm,nxy_s,nxy_e)
            call range(1,ncx,ncpu_s,mm,l2x_s,l2x_e)
            ii       = ii + 1
c
            ix_s(ii)  = nxy_s
            ix_e(ii)  = nxy_e
            i2x_s(ii) = (nxy_s - 1)*2 + 1
            i2x_e(ii) = nxy_e*2
            jx_s(ii)  = (l2x_s - 1)*2 + 1
            jx_e(ii)  = l2x_e*2
            kx_s(ii)  = lx_s
            kx_e(ii)  = lx_e
            mx_s(ii)  = nx_s
            mx_e(ii)  = nx_e
c
            iy_s(ii)  = ny_s
            iy_e(ii)  = ny_e
            jy_s(ii)  = ly_s
            jy_e(ii)  = ly_e
c
            iz_s(ii)  = mz_s
            iz_e(ii)  = mz_e
c
            is_s(ii)  = (ii/ncpu_s)*ncpu_s
            is_e(ii)  = is_s(ii) + ncpu_s - 1
         enddo
      enddo
c
      iys  =  iy_s(myid)
      iye  =  iy_e(myid)
      jys  =  jy_s(myid)
      jye  =  jy_e(myid)
      ixs  =  ix_s(myid)
      ixe  =  ix_e(myid)
      jxs  =  jx_s(myid)
      jxe  =  jx_e(myid)
      kxs  =  kx_s(myid)
      kxe  =  kx_e(myid)
      mxs  =  mx_s(myid)
      mxe  =  mx_e(myid)
      izs  =  iz_s(myid)
      ize  =  iz_e(myid)
      i2xs =  i2x_s(myid)
      i2xe =  i2x_e(myid)
c
c ----------- get starting and  ending processor id's on each
c             vertical slab
c
      iss = is_s(myid)
      ise = is_e(myid)
c
c ------------ debug ranges
c
      if(l_debug) then
         write(nprt,1200) myid, (nn, ix_s(nn), ix_e(nn), jx_s(nn),
     +                     jx_e(nn), kx_s(nn), kx_e(nn),
     +                     i2x_s(nn), i2x_e(nn),
     +                     nn = 0,numprocs-1)
 1200    format(' myid =  ',i4,/,
     +       ' nn',5x,' ixs ',5x,' ixe ',5x,' jxs ',5x,' jxe '
     +       ,5x,' kxs ',5x,' kxe',5x,' i2xs ',5x,' i2xe',/,(9i6))
c
         write(nprt,1213) myid, (nn, iy_s(nn), iy_e(nn),
     +                   jy_s(nn), jy_e(nn),
     +                   iz_s(nn), iz_e(nn), is_s(nn), is_e(nn),
     +                   nn=0,numprocs-1)
 1213    format(' myid = ',i3,/,
     +       ' nn ',3x,' iys ',5x,' iye ',5x,
     +       ' jys ',5x,' jye ',5x,
     +       ' izs ',5x,' ize',5x,' iss ',5x,' ise ',/,
     +       (9i6))
      endif
c
      return
      end
      subroutine mpi_sum_xy(f,myid,iss,ise,nsend)
c
c --------- get horizontal x-y sum over a set of proccessors [iss:ise]
c           for vector f(i). f(i) is overwritten. skip if single processor
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
      real work(nsend,iss:ise), f(nsend)
c
      if(iss .eq. ise) go to 999
c
      do j=1,nsend
         work(j,myid) = f(j)
         f(j)         = 0.0
      enddo
      do i=iss,ise
         if(i .ne. myid) then
            call mpi_sendrecv(work(1,myid),nsend,mpi_real8,i,1,
     +               work(1,i),nsend,mpi_real8,i,1,
     +           mpi_comm_world,istatus,ierr)
         endif
      enddo
      do i=iss,ise
      do j=1,nsend
         f(j) = f(j) + work(j,i)
      enddo
      enddo
c
  999 continue
c
      return
      end
      subroutine mpi_sum_z(f,i_root,myid,nsend,iall)
c
c --------- get sums on root or all processors
c           for all z for vector f(i)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
      real recv_b(nsend), f(nsend)
c
c -------- just root gets the result
c
      if(iall .ne. 1) then
         call mpi_reduce(f(1),recv_b(1),nsend,mpi_real8,mpi_sum,i_root,
     +                  mpi_comm_world,ierr)
         if(myid .eq. i_root) then
            do i=1,nsend
               f(i) = recv_b(i)
            enddo
         endif
      else
c
c -------- everyone gets the result
c
         call mpi_allreduce(f(1),recv_b(1),nsend,mpi_real8,mpi_sum,
     +                  mpi_comm_world,ierr)
         do i=1,nsend
            f(i) = recv_b(i)
         enddo
      endif
c
      return
      end
      subroutine mpi_sum_z_s(f,i_root,myid,nsend,nscl,iall)
c
c --------- get sums on root or all processors
c           for all z for vector f(i,nscl)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
      real recv_b(nsend,nscl), f(nsend,nscl)
c
      if(iall .ne. 1) then
         call mpi_reduce(f(1,1),recv_b(1,1),nsend*nscl,mpi_real8,
     +        mpi_sum,i_root,mpi_comm_world,ierr)
         if(myid .eq. i_root) then
            do iscl=1,nscl
            do i=1,nsend
               f(i,iscl) = recv_b(i,iscl)
            enddo
            enddo
         endif
      else
         call mpi_allreduce(f(1,1),recv_b(1,1),nsend*nscl,mpi_real8,
     +        mpi_sum, mpi_comm_world,ierr)
         do iscl=1,nscl
         do i=1,nsend
            f(i,iscl) = recv_b(i,iscl)
         enddo
         enddo
      endif
c
      return
      end
      subroutine mpi_gath_root(fs,fr,iz_s,iz_e,izs,ize,nz,myid,np,ns)
c
c ---------- gather results on root processors
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
      integer iz_s(0:np-1), iz_e(0:np-1)
      real fs(izs:ize), fr(1:nz)
c
      if(np .eq. 1) go to 999
c
      irow_r = mod(myid,ns)
      if(myid .gt. ns) then
        call mpi_send(fs(izs),ize+1-izs,mpi_real8,irow_r,1,
     +       mpi_comm_world,ierr)
      else
        do l=irow_r+ns,np-1,ns
           ind = iz_s(l) + 1
           num = iz_e(l) + 1 - iz_s(l)
           call mpi_recv(fr(ind),num,mpi_real8,l,1,
     +       mpi_comm_world,istatus,ierr)
        enddo
      endif
c
  999 continue
c
      return
      end
      subroutine mpi_send_root(fs,num,myid,np,ns)
c
c ---------- send root results to other processors above it
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
      real fs(num)
c
      if(np .eq. 1) go to 999
c
      irow_r = mod(myid,ns)
      if(myid .ge. ns) then
        call mpi_recv(fs(1),num,mpi_real8,irow_r,1,
     +       mpi_comm_world,istatus,ierr)
      else
        do l=irow_r+ns,np-1,ns
           call mpi_send(fs(1),num,mpi_real8,l,1,
     +          mpi_comm_world,ierr)
        enddo
      endif
c
  999 continue
c
      return
      end
      subroutine xtoy_trans(f,g,nx,ny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,
     +           myid,ncpu_s,np)
c
c ------- transpose array  f(nx,iys:iye,iz1:iz2) ---> g(ny,ixs:ixe,iz1:iz2)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      real f(nx,iys:iye,iz1:iz2),
     +     g(ny,ixs:ixe,iz1:iz2)
      real ft(nx*(iye+1-iys)*(iz2 - iz1 + 1)),
     +     gt(ny*(ixe+1-ixs)*(iz2 - iz1 + 1))
      integer ix_s(0:np-1), ix_e(0:np-1),
     +        iy_s(0:np-1), iy_e(0:np-1)
c
      jk = (iye - iys + 1)*(iz2 - iz1 + 1)
      ik = (ixe - ixs + 1)*(iz2 - iz1 + 1)
c
c ----------- get cpus on slab for myid
c
      islab = myid/ncpu_s
      iss   = islab*ncpu_s
      ise   = iss + ncpu_s - 1
c
      do i=iss,ise
         nsend = (ix_e(i) - ix_s(i) + 1)*jk
         nrecv = (iy_e(i) - iy_s(i) + 1)*ik
         if(i .eq. myid) then
            call send_xtoy(f,gt(1),nx,ix_s(i),ix_e(i),
     +                  iy_s(myid),iy_e(myid),iz1,iz2)
         else
            call send_xtoy(f,ft(1),nx,ix_s(i),ix_e(i),
     +                  iy_s(myid),iy_e(myid),iz1,iz2)
            call mpi_sendrecv(
     +          ft(1),nsend,mpi_real8,i,1,
     +          gt(1),nrecv,mpi_real8,i,1,
     +          mpi_comm_world,istatus,ierr)
         endif
         call recv_xtoy(g,gt(1),ny,ix_s(myid),ix_e(myid),
     +                  iy_s(i),iy_e(i),iz1,iz2)
      enddo
c
      return
      end
      subroutine send_xtoy(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
c
c ------------- grab correct chunk of array to be sent
c
      real f(nx,iys:iye,izs:ize), ft(ixs:ixe,iys:iye,izs:ize)
c
      do k=izs,ize
      do j=iys,iye
      do i=ixs,ixe
         ft(i,j,k) = f(i,j,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine recv_xtoy(g,gt,ny,ixs,ixe,iys,iye,izs,ize)
      real g(ny,ixs:ixe,izs:ize), gt(ixs:ixe,iys:iye,izs:ize)
c
      do k=izs,ize
      do j=iys,iye
      do i=ixs,ixe
         g(j,i,k) = gt(i,j,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine ytox_trans(g,f,nx,ny,ixs,ixe,ix_s,ix_e,
     +           iys,iye,iy_s,iy_e,iz1,iz2,
     +           myid,ncpu_s,np)
c
c ------- transpose array g(ny,ixs:ixe,iz1:iz2) ---> f(nx,iys:iye,iz1:iz2)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      real f(nx,iys:iye,iz1:iz2),
     +     g(ny,ixs:ixe,iz1:iz2)
      real ft(nx*(iye+1-iys)*(iz2 - iz1 + 1)),
     +     gt(ny*(ixe+1-ixs)*(iz2 - iz1 + 1))
c
      integer ix_s(0:np-1), ix_e(0:np-1),
     +        iy_s(0:np-1), iy_e(0:np-1)
c
      jk = (iye - iys + 1)*(iz2 - iz1 + 1)
      ik = (ixe - ixs + 1)*(iz2 - iz1 + 1)
c
c ----------- get cpus on slab for myid
c
      islab = myid/ncpu_s
      iss   = islab*ncpu_s
      ise   = iss + ncpu_s - 1
      do i=iss,ise
         nsend = (iy_e(i) - iy_s(i) + 1)*ik
         nrecv = (ix_e(i) - ix_s(i) + 1)*jk
         if(i .eq. myid) then
            call send_ytox(g,ft(1),ny,ix_s(myid),ix_e(myid),
     +                  iy_s(i),iy_e(i),iz1,iz2)
         else
            call send_ytox(g,gt(1),ny,ix_s(myid),ix_e(myid),
     +                  iy_s(i),iy_e(i),iz1,iz2)
            call mpi_sendrecv(
     +          gt(1),nsend,mpi_real8,i,1,
     +          ft(1),nrecv,mpi_real8,i,1,
     +          mpi_comm_world,istatus,ierr)
         endif
         call recv_ytox(f,ft(1),nx,ix_s(i),ix_e(i),
     +                  iy_s(myid),iy_e(myid),iz1,iz2)
      enddo
c
      return
      end
      subroutine send_ytox(g,gt,ny,ixs,ixe,iys,iye,izs,ize)
c
c ------------- grab correct chunk of array to be sent
c
      real g(ny,ixs:ixe,izs:ize), gt(iys:iye,ixs:ixe,izs:ize)
c
      do k=izs,ize
      do i=ixs,ixe
      do j=iys,iye
         gt(j,i,k) = g(j,i,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine recv_ytox(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
      real f(nx,iys:iye,izs:ize), ft(iys:iye,ixs:ixe,izs:ize)
c
      do k=izs,ize
      do i=ixs,ixe
      do j=iys,iye
         f(i,j,k) = ft(j,i,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine xtoz_trans(f,g,nx,nz,ixs,ixe,ix_s,ix_e,
     +           iys,iye,izs,ize,iz_s,iz_e,
     +           myid,ncpu_s,numprocs)
c
c ------- transpose array  f(nx,iys:iye,izs-1:ize+1)
c                     ---> g(0:nz+1,iys:iye,ixs:ixe)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      real f(nx,iys:iye,izs-1:ize+1), g(0:nz+1,iys:iye,ixs:ixe)
      real ft(nx*(iye+1-iys)*(ize-izs+1)),
     +     gt(nz*(ixe+1-ixs)*(iye-iys+1))
      integer ix_s(0:numprocs-1), ix_e(0:numprocs-1),
     +        iz_s(0:numprocs-1), iz_e(0:numprocs-1)
c
      jk = (ize - izs + 1)*(iye - iys + 1)
      ij = (ixe - ixs + 1)*(iye - iys + 1)
c
c ----------- get starting location
c
      iss = myid - (myid/ncpu_s)*ncpu_s
c
      do i=iss,numprocs-1,ncpu_s
         nsend = (ix_e(i) - ix_s(i) + 1)*jk
         nrecv = (iz_e(i) - iz_s(i) + 1)*ij
         if(i .eq. myid) then
            call send_xtoz(f,gt(1),nx,ix_s(i),ix_e(i),
     +                  iys,iye,iz_s(myid),iz_e(myid))
         else
            call send_xtoz(f,ft(1),nx,ix_s(i),ix_e(i),
     +                  iys,iye,iz_s(myid),iz_e(myid))
            call mpi_sendrecv(
     +          ft(1),nsend,mpi_real8,i,1,
     +          gt(1),nrecv,mpi_real8,i,1,
     +          mpi_comm_world,istatus,ierr)
         endif
         call recv_xtoz(g,gt(1),nz,ix_s(myid),ix_e(myid),
     +                  iys,iye,iz_s(i),iz_e(i))
      enddo
c
      return
      end
      subroutine send_xtoz(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
c
c ------- grab correct chunk of array to be sent and skip ghost points
c
      real f(nx,iys:iye,izs-1:ize+1), ft(ixs:ixe,iys:iye,izs:ize)
c
      do k=izs,ize
      do j=iys,iye
      do i=ixs,ixe
         ft(i,j,k) = f(i,j,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine recv_xtoz(g,gt,nz,ixs,ixe,iys,iye,izs,ize)
      real g(0:nz+1,iys:iye,ixs:ixe), gt(ixs:ixe,iys:iye,izs:ize)
c
      do k=izs,ize
      do j=iys,iye
      do i=ixs,ixe
         g(k,j,i) = gt(i,j,k)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine ztox_trans(g,f,nx,nz,ixs,ixe,ix_s,ix_e,
     +           iys,iye,izs,ize,iz_s,iz_e,
     +           myid,ncpu_s,numprocs)
c
c ------- transpose array g(0:nz+1,iys:iye,ixs:ixe)
c                    ---> f(nx,iys:iye,izs-1:ize+1)
c
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      real f(nx,iys:iye,izs-1:ize+1), g(0:nz+1,iys:iye,ixs:ixe)
      real ft(nx*(iye+1-iys)*(ize-izs+3)),
     +     gt((nz+3)*(iye+1-iys)*(ixe-ixs+1))
c
      integer ix_s(0:numprocs-1), ix_e(0:numprocs-1),
     +        iz_s(0:numprocs-1), iz_e(0:numprocs-1)
c
      jk = (ize - izs + 3)*(iye - iys + 1)
      ij = (ixe - ixs + 1)*(iye - iys + 1)
c
c ------------- get starting location
c
      iss = myid - (myid/ncpu_s)*ncpu_s
c
      do i=iss,numprocs-1,ncpu_s
         nsend = (iz_e(i) - iz_s(i) + 3)*ij
         nrecv = (ix_e(i) - ix_s(i) + 1)*jk
         if(i .eq. myid) then
            call send_ztox(g,ft(1),nz,ix_s(myid),ix_e(myid),
     +                  iys,iye,iz_s(i),iz_e(i))
         else
            call send_ztox(g,gt(1),nz,ix_s(myid),ix_e(myid),
     +                  iys,iye,iz_s(i),iz_e(i))
            call mpi_sendrecv(
     +          gt(1),nsend,mpi_real8,i,1,
     +          ft(1),nrecv,mpi_real8,i,1,
     +          mpi_comm_world,istatus,ierr)
         endif
         call recv_ztox(f,ft(1),nx,ix_s(i),ix_e(i),
     +                  iys,iye,iz_s(myid),iz_e(myid))
      enddo
c
      return
      end
      subroutine send_ztox(g,gt,nz,ixs,ixe,iys,iye,izs,ize)
c
c ------------- grab correct chunk of array to be sent,
c               account for ghost points
c
      real g(0:nz+1,iys:iye,ixs:ixe), gt(izs-1:ize+1,iys:iye,ixs:ixe)
c
      do j=iys,iye
      do i=ixs,ixe
      do k=izs-1,ize+1
         gt(k,j,i) = g(k,j,i)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine recv_ztox(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
      real f(nx,iys:iye,izs-1:ize+1), ft(izs-1:ize+1,iys:iye,ixs:ixe)
c
      do i=ixs,ixe
      do j=iys,iye
      do k=izs-1,ize+1
         f(i,j,k) = ft(k,j,i)
      enddo
      enddo
      enddo
c
      return
      end
      subroutine exchange
c
c ------------- exchange ghost points with mpi,
c               nb and nt are the destination and
c               source nodes. Allows for 1z per cpu
c
      use pars
      use fields
c     use fftwk
      include 'mpif.h'
c
      real fs(nnx,iys:iye,(5+nscl)),fr(nnx,iys:iye,(5+nscl))
      integer istatus(mpi_status_size)
c
      nb = myid - ncpu_s
      nt = myid + ncpu_s
c
c ------------ account for endpoints
c
      if(iss .eq. 0) then
         nb = mpi_proc_null
      endif
      if(ise .eq. numprocs-1) then
         nt = mpi_proc_null
      endif
      nsend = nnx*(iye + 1 - iys)*(5+nscl)
      nrecv = nsend
c
c --------- send top of myid, receive bottom from myid - ncpu_s
c
      do iy=iys,iye
      do ix=1,nnx
         fs(ix,iy,1) = u(ix,iy,ize)
         fs(ix,iy,2) = v(ix,iy,ize)
         fs(ix,iy,3) = w(ix,iy,ize)
         fs(ix,iy,4) = e(ix,iy,ize)
         fs(ix,iy,5) = wf(ix,iy,ize)
      enddo
      enddo
      do iscl=1,nscl
         jloc = 5 + iscl
         do iy=iys,iye
         do ix=1,nnx
            fs(ix,iy,jloc) = t(ix,iy,iscl,ize)
         enddo
         enddo
      enddo
c
      call mpi_sendrecv(
     +     fs(1,iys,1),nsend,mpi_real8,nt,0,
     +     fr(1,iys,1),nrecv,mpi_real8,nb,0,
     +     mpi_comm_world,istatus,ierr)
c
      if(iss .ne. 0) then
         izm1 = izs-1
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,izm1) = fr(ix,iy,1)
            v(ix,iy,izm1) = fr(ix,iy,2)
            w(ix,iy,izm1) = fr(ix,iy,3)
            e(ix,iy,izm1) = fr(ix,iy,4)
            wf(ix,iy,izm1) = fr(ix,iy,5)
         enddo
         enddo
         do iscl=1,nscl
            jloc = 5 + iscl
            do iy=iys,iye
            do ix=1,nnx
               t(ix,iy,iscl,izm1) = fr(ix,iy,jloc)
            enddo
            enddo
         enddo
      endif
c
c -------- send bottom of myid, receive bottom from myid + ncpu_s
c
      do iy=iys,iye
      do ix=1,nnx
         fs(ix,iy,1) = u(ix,iy,izs)
         fs(ix,iy,2) = v(ix,iy,izs)
         fs(ix,iy,3) = w(ix,iy,izs)
         fs(ix,iy,4) = e(ix,iy,izs)
      enddo
      enddo
      do iscl=1,nscl
         jloc = 4 + iscl
         do iy=iys,iye
         do ix=1,nnx
            fs(ix,iy,jloc) = t(ix,iy,iscl,izs)
         enddo
         enddo
      enddo
      nsend = nnx*(iye + 1 - iys)*(4+nscl)
      nrecv = nsend
c
      call mpi_sendrecv(
     +     fs(1,iys,1),nsend,mpi_real8,nb,1,
     +     fr(1,iys,1),nrecv,mpi_real8,nt,1,
     +     mpi_comm_world,istatus,ierr)
c
      if(ise .ne. numprocs-1) then
         izp1 = ize+1
         do iy=iys,iye
         do ix=1,nnx
            u(ix,iy,izp1) = fr(ix,iy,1)
            v(ix,iy,izp1) = fr(ix,iy,2)
            w(ix,iy,izp1) = fr(ix,iy,3)
            e(ix,iy,izp1) = fr(ix,iy,4)
         enddo
         enddo
         do iscl=1,nscl
            jloc = 4 + iscl
            do iy=iys,iye
            do ix=1,nnx
               t(ix,iy,iscl,izp1) = fr(ix,iy,jloc)
            enddo
            enddo
         enddo
      endif
c
c --------------- send extra scalar points
c
      nsend = nnx*(iye + 1 - iys)*nscl
      nrecv = nsend
c
c -------------- send top of myid, receive bottom from myid - ncpu_s
c
      izm1 = ize-1
      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            fs(ix,iy,iscl) = t(ix,iy,iscl,izm1)
         enddo
         enddo
      enddo
c
      call mpi_sendrecv(
     +     fs(1,iys,1),nsend,mpi_real8,nt,0,
     +     fr(1,iys,1),nrecv,mpi_real8,nb,0,
     +     mpi_comm_world,istatus,ierr)
c
      if(iss .ne. 0) then
         izm2 = izs-2
         do iscl=1,nscl
            do iy=iys,iye
            do ix=1,nnx
               t(ix,iy,iscl,izm2) = fr(ix,iy,iscl)
            enddo
            enddo
         enddo
      endif
c
c -------------- send bottom of myid, receive bottom from myid + ncpu_s
c
      izp1 = izs+1
      do iscl=1,nscl
         do iy=iys,iye
         do ix=1,nnx
            fs(ix,iy,iscl) = t(ix,iy,iscl,izp1)
         enddo
         enddo
      enddo
c
      call mpi_sendrecv(
     +     fs(1,iys,1),nsend,mpi_real8,nb,1,
     +     fr(1,iys,1),nrecv,mpi_real8,nt,1,
     +     mpi_comm_world,istatus,ierr)
c
      if(ise .ne. numprocs-1) then
         izp2 = ize+2
         do iscl=1,nscl
            do iy=iys,iye
            do ix=1,nnx
               t(ix,iy,iscl,izp2) = fr(ix,iy,iscl)
            enddo
            enddo
         enddo
      endif
c
      return
      end
      subroutine bcast_pbc
c
c ---- send upper boundary conditions to other processors
c      for fft solution of pressure
c
      use pars
      use fields
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      if(numprocs .eq. 1) go to 999
c
      irow_r = mod(myid,ncpu_s)
      irow_t = is_s(numprocs-1) + irow_r
      num = nnx*(iye+1-iys)
c
c ----- check which row myid is in
c
      if(iss .ne. is_s(numprocs-1)) then
c
c ------ not in the top row, receive from top
c
        call mpi_recv(pbc(1,iys,1),num,mpi_real8,irow_t,1,
     +       mpi_comm_world,istatus,ierr)
      else
c
c ------ myid is in the top row, send to everyone below
c
        do l=irow_r,irow_t-ncpu_s,ncpu_s
           call mpi_send(pbc(1,iys,1),num,mpi_real8,l,1,
     +          mpi_comm_world,ierr)
        enddo
      endif
c
c --------- same thing for another variable
c
      if(iss .ne. is_s(numprocs-1)) then
c
c ------ not in the top row, receive from top
c
        call mpi_recv(pbc2(1,iys,1),num,mpi_real8,irow_t,1,
     +       mpi_comm_world,istatus,ierr)
      else
c
c ------ in the top row, send to everyone below
c
        do l=irow_r,irow_t-ncpu_s,ncpu_s
           call mpi_send(pbc2(1,iys,1),num,mpi_real8,l,1,
     +          mpi_comm_world,ierr)
        enddo
      endif
c
  999 continue
c
      return
      end
      subroutine bcast_surf
c
c ---- send bottom geometry to all other processors
c      for computation of pressure drag at all z
c
      use pars
      use fields
      include 'mpif.h'
      integer istatus(mpi_status_size)
c
      if(iss == 0) then
        k = 1
        do j=iys,iye
        do i=1,nnx
           gms(i,j,1) = gm(i,j,k,1)
           gms(i,j,2) = gm(i,j,k,2)
           gms(i,j,3) = gm(i,j,k,3)
        enddo
        enddo
      endif
c
      if(numprocs == 1) go to 999
c
      irow_r = mod(myid,ncpu_s)
      irow_t = is_s(numprocs-1) + irow_r
      num    = nnx*(iye+1-iys)
c
c ----- check which row myid is in
c
      if(iss /= is_s(0)) then
c
c ------ not in the bottom row, receive from bottom
c
        call mpi_recv(gms(1,iys,1),num,mpi_real8,irow_r,1,
     +       mpi_comm_world,istatus,ierr)
      else
c
c ------ myid is in the bottom row, send to everyone else
c
        do l=irow_r+ncpu_s,irow_t,ncpu_s
           call mpi_send(gms(1,iys,1),num,mpi_real8,l,1,
     +          mpi_comm_world,ierr)
        enddo
      endif
c
c --------- same thing for another variable
c
      if(iss /= is_s(0)) then
c
c ------ not in the bottom row, receive from bottom
c
        call mpi_recv(gms(1,iys,2),num,mpi_real8,irow_r,1,
     +       mpi_comm_world,istatus,ierr)
      else
c
c ------ myid is in the bottom row, send to everyone else
c
        do l=irow_r+ncpu_s,irow_t,ncpu_s
           call mpi_send(gms(1,iys,2),num,mpi_real8,l,1,
     +          mpi_comm_world,ierr)
        enddo
      endif
c
c --------- same thing for another variable
c
      if(iss /= is_s(0)) then
c
c ------ not in the bottom row, receive from bottom
c
        call mpi_recv(gms(1,iys,3),num,mpi_real8,irow_r,1,
     +       mpi_comm_world,istatus,ierr)
      else
c
c ------ myid is in the bottom row, send to everyone else
c
        do l=irow_r+ncpu_s,irow_t,ncpu_s
           call mpi_send(gms(1,iys,3),num,mpi_real8,l,1,
     +          mpi_comm_world,ierr)
        enddo
      endif
  999 continue
c
      return
      end


      function ran2(idum)
      integer :: idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real :: ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     +     IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     +     IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER :: idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/,iv/NTAB*0/,iy/0/

      if (idum .le. 0) then
          idum=max(-idum,1)
          idum2 = idum
          do j = NTAB+8,1,-1
             k=idum/IQ1
             idum=IA1*(idum-k*IQ1)-k*IR1
             if (idum .lt. 0) idum=idum+IM1
             if (j .le. NTAB) iv(j) = idum
          end do
          iy=iv(1)
      end if
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum .lt. 0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2 .lt. 0) idum2=idum2+IM2
      j = 1+iy/NDIV
      iy = iv(j) - idum2
      iv(j) = idum
      if (iy .lt. 1) iy = iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      end function ran2

      subroutine linInterp(xa,ya,n,x,y)
      integer n
      real x, xper, y, xa(n), ya(n)
      integer i1, i2
      xper = x
      if (xper .ge. xa(n)) then
        xper = xper-xa(n)
      else if (xper < 0) then
        xper = xa(n)+xper
      endif
      i1 = floor((n-1)*xper/xa(n)) + 1
      i2 = i1+1
      y = ya(i1) + (xper-xa(i1))*(ya(i2)-ya(i1))/(xa(i2)-xa(i1))
      return
      end
