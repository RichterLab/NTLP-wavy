
      subroutine get_max

c --------- routine computes max velocities as sweep through
c           the flux velocity field

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        include 'mpif.h'

        real u_send(6), u_recv(6)

        dx_i = 1.0/dx
        dy_i = 1.0/dy

        u_temp   = 0.0
        v_temp   = 0.0
        w_temp   = 0.0
        e_temp   = 0.0
        vis_temp = 0.0
        do iz=izs,ize

          u_xy = 0.0
          v_xy = 0.0
          w_xy = 0.0
          e_xy = 0.0
          g_xy = gm(1,iys,iz,1)
          do iy=iys,iye
            do ix=1,nnx
              u_xy = amax1(u_xy,abs(uf(ix,iy,iz)))
              v_xy = amax1(v_xy,abs(vf(ix,iy,iz)))
              w_xy = amax1(w_xy,abs(wf(ix,iy,iz)),abs(gm(ix,iy,iz,5)))
              e_xy = amax1(e_xy,e(ix,iy,iz))
              g_xy = amax1(g_xy,gm(ix,iy,iz,1))
            enddo
          enddo
          u_xy   = u_xy*dx_i
          v_xy   = v_xy*dy_i
          wsav   = w_xy
          w_xy   = w_xy/abs(dzw(iz))

c ------- get max allowable diffusivity for scalars

          vis_xy = dfac(iz)*3.0*ck*dsl_z(iz)*(g_xy**pow3)*
     +        sqrt(e_xy)/amin1(dx,dy,dzw(iz))**2

          u_temp   = amax1(u_xy,u_temp)
          v_temp   = amax1(v_xy,v_temp)
          w_temp   = amax1(w_xy,w_temp)
          e_temp   = amax1(e_xy,e_temp)
          vis_temp = amax1(vis_xy,vis_temp)

        enddo
        u_send(1) = u_temp
        u_send(2) = v_temp
        u_send(3) = w_temp
        u_send(4) = wsav
        u_send(5) = e_temp
        u_send(6) = vis_temp

        call mpi_allreduce(u_send,u_recv,6,mpi_real8,
     +      mpi_max,mpi_comm_world,ierror)

        umax   = u_recv(1)
        vmax   = u_recv(2)
        wmax   = u_recv(3)
        wabs   = u_recv(4)
        emax   = u_recv(5)
        vismax = u_recv(6)

        return
      end
      subroutine get_dt(it)

c ---------- routine computes max time step for given cfl number
c            from max's found previously

        use pars
        use con_data
        use con_stats

        ucfl = umax
        vcfl = vmax
        wcfl = wmax
        vel_max = amax1(ucfl, vcfl, wcfl)
        if(vel_max <= 0.0) then
          write(6,6000) ucfl, vcfl, wcfl, vel_max
 6000     format('6000, sr. get_dt bad news, umax = ',e15.6,/,
     +        ' vmax = ',e15.6,' wmax = ',e15.6,/,
     +        ' vel_max = ',e15.5,/,
     +        ' infinite time step !!!')
          stop
        endif

c ---------------- choose fixed or variable time step

        if(ifix_dt /= 0) then

c ------------- if used, change to fit your problem

          dt_new = 0.000370099
        else

c ------------------- new estimate of best time step
c                     from cfl constraint

          dt_new  = cfl/vel_max
          dt_new = amin1(dt_new, 5.0)
        endif

c ---------------- compare against viscous stability limit

        if(vismax*dt_new > 0.5) then
          dt_cfl = dt_new
          dt_new = 0.5/vismax
          if(l_root) then
            write(6,6200) dt_new, dt_cfl, vismax
 6200       format(' 6200 get_dt: cfl time step too large',/,
     +          '   viscous time step = ',e15.6,
     +          ' cfl time step = ',e15.6,' vismax = ',e15.6)
          endif
        endif

c -------- for safety if restart set timestep = saved timestep in
c          restart file

CHECK ... for moving grids?

        if(it == iti .and. iti /= 0) then
          dt_new = dt1
          dt1    = dt
          write(nprt,4001) it, iti, dt_new, dt1, dt
 4001     format(' 4001: it = ',i5,' iti = ',i5,
     +        ' dt_new = ',e15.6,' dt1 = ',e15.6,
     +        ' dt = ',e15.6)
        endif

        return
      end
      subroutine lterp(n,zary,zpt,i,ip1,ratio)

c ---- linear interpolation for zpt in zary, where zary is
c      monotonic increasing or decreasing function

        dimension zary(*)
        nm1 = n-1
        if(n<=1) then
          i = 1
          ip1 = 1
          ratio = 0.0
          go to 999
        endif
        if(zary(1) < zary(2)) go to 1
        go to 101
    1   continue

c **** monotonic increasing array

        if(zpt < zary(1)) then
          i = 1
          ip1 = 1
          ratio = 0.0
          go to 999
        elseif(zpt > zary(n)) then
          i = n
          ip1 = n
          ratio = 1.0
          go to 999
        endif
        do j=1,nm1
          if(zpt >= zary(j) .and.
     $        zpt <= zary(j+1)) then
            i = j
            ip1 = j+1
            ratio = (zpt - zary(i))/(zary(ip1) - zary(i))
            go to 999
          endif
        enddo

c **** decreasing array

  101   continue
        if(zpt > zary(1)) then
          i = 1
          ip1 = 1
          ratio = 0.0
          go to 999
        elseif(zpt < zary(n)) then
          i = n
          ip1 = n
          ratio = 1.0
          go to 999
        endif
        do j=1,nm1
          if(zpt <= zary(j) .and.
     $        zpt >= zary(j+1)) then
            i = j
            ip1 = j+1
            ratio = (zpt - zary(i))/(zary(ip1) - zary(I))
            go to 999
          endif
        enddo
  999   continue
        return
      end
      subroutine setcon

        use pars
        use fftwk
        use con_data
        use con_stats
        include 'mpif.h'
        external get_zi

c ----------------- get machine type, can read datadir also

        open(99,file='mach.file')
        read(99,9000) imach
 9000   format(i1)

        close(99)

c ----------- initialize fft

        call rffti(nnx,trigx(1,1))
        call rffti(nny,trigx(1,2))
        call cffti(nny,trigc(1))
        call cffti(nnx,trigcx(1))

c ----------- start step for history and viz files

        it_his_nxt = it_his
        it_viz_nxt = it_viz

c ---------------- set min value of e

        if(iocean == 1) then

          smal_e = 0.0
          smal_e = 1.0e-12
c        smal_e = 6.0e-03
        else
          smal_e = 1.0e-09
c        smal_e = 0.0
        endif

c ---------------------- set constants in eddy viscosity model

        ck       = 0.1
        ceps     = 0.93
        csmag    = sqrt(ck*sqrt(ck/ceps))
        stab_c   = 0.76
        pow3     = 1.0/3.0

c ----------------- set stability constant

        stabmin = 1.0e-12

c ---------------- minimum dsl length constant

        almin_c = 0.0001

c -------------------- set viscosity model parameters
c                      to usual values, overridden by iso
c                      dfac is in constant file for calculation
c                      of max eddy viscosity

        viscon      = 0.0
        xksurf      = 0.0
        nmatch      = -1
        myid_newvis = 0
        do iz=1,nnz
          dfac(iz) = 1.0
        enddo

c -------------- initialize grid restart flag

        igrdr = 1

c ---------- set constants for pressure iteration

        p_surf = 0.0

        if(iflat == 0) then
          ip_max    = 2
          ptol      = 0.0
          accel_max = 1.0
          accel_min = 1.0
          do i=1,ip_accel
            accel(i) = 1.0
          enddo
        else

c ------- with effectively flat just use 2 iterations

          ip_max  = 2
          ptol    = 1.0e-07
c        ptol    = 1.0e-20

          if(i3d_surf == 3 .or. i3d_surf == 4) then   !  strategy for breaking waves in tank
            ptol      = 1.0e-10
            ip_max    = 50
            ip_max    = 35
            accel_max = 1.0
            accel_min = 0.1
            ip_ramp   = 30

c ---------- for initialization with flat

          else                     !  strategy for pbl with surface waves or hills
            ptol      = 1.0e-10
            ip_max    = 30
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
 1017     format(' acceleration parameters',/,
     +        ' i ',5x,' accel',/,(i5,e15.6))
        endif

c -------------- set grid generation factor for type of simulation
c                in routines exp_grid, exp_grid2

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
        elseif(i3d_surf == 4) then
          a_exp_grd   = 11.65
          a_exp_grd2  = -1000.0 ! fix for tank simulations
          a_decay     = 3.0
        endif

c -------------- create mpi operation to find max and location
c                using local gradient method

        call mpi_op_create(get_zi,.true.,ziloc,ierror)

c ---------- define coefficients for 3-order runge-kutta
c            time stepping scheme, borrowed from Spalart,
c            Moser and Rogers, J. Comp. Physics 3/21/90
c            Note this is a simplier version since all terms
c            are lumped in the non-linear terms.
c            cfl number is for an entire runge-kutta step
c            in this case three stages. cfl = max(u)*dt/dx

c            see jcp article by berkooz for time steps

        zetas(1) = 0.0
        zetas(2) = -17.0/60.0
        zetas(3) = -5.0/12.0
        gama(1)  = 8.0/15.0
        gama(2)  = 5.0/12.0
        gama(3)  = 3.0/4.0

        etas(1)  = 0.0
        etas(2)  = 8.0/15.0
        etas(3)  = 2.0/3.0

c ----------- a full step, at the new time

        etas(4)  = 1.0

c --------- set stationary flag on a restart
c           see sr. set_sav for tests

        imesh = 0

        return

      end
      subroutine setup(it)

        use particles
        use pars
        use fields
        use fftwk
        use con_data
        use con_stats

        it = iti

c ------------ turn on new sgs model at a particular step

        if(it >= new_vis .and. ivis0 == 1) then
          ivis = 1
        else
          ivis = 0
        endif

        if(igrdr . eq. 3) then
          if(l_root) then
            write(6,6)iti,utau,tsfcc(1) ,qstar(1)
            write(6,510)
            write(6,520) wtsfc(1),utau,amonin,dtdzf(1),zody,zo
     +          ,cdbtm,ugcont
            call print(6,it,1,nnz)
          endif
          if(l_debug) then
            write(nprt,6)iti,utau,tsfcc(1) ,qstar(1)
            write(nprt,510)
            write(nprt,520) wtsfc(1),utau,amonin,dtdzf(1),zody,zo
     +          ,cdbtm,ugcont
            call print(nprt,it,izs,ize)
          endif
        endif
c     if(ifilt==1)call filter
        if(l_root) then
          write(6,1) nnx,nny,nnz,ismlt,ifilt,iti,itmax,
     +        iupwnd,ibuoy,itcut,
     +        dt,zo,tsfcc(1),isubs,ibrcl,
     +        method, iocean, ivis
        endif
        if(l_debug) then
          write(nprt,1) nnx,nny,nnz,ismlt,ifilt,iti,itmax,
     +        iupwnd,ibuoy,itcut,
     +        dt,zo,tsfcc(1),isubs,ibrcl,
     +        method, iocean, ivis
        endif

c -------------- boundary condition flags

      !ibcu = 2
        ibcu = -1
        ibcl = 0

c -------------------- wavenumbers, introduce a normalized
c                      set of wavenumbers to eliminate computation
c                      in derivatives , xderiv, yderiv

        do i=1,nnx
          xkn(i) = float(i-1)*pi2/xl
          if(i>ncx)xkn(i) = -float(nnx-i+1)*pi2/xl
        enddo
        fn = 1.0/float(nnx)
        do i=1,nnx
          xk(i)  = xkn(i)*fn
          xkk(i) = xkn(i)*xkn(i)*fn
        enddo
        do i=1,nny
          ykn(i) = float(i-1)*pi2/yl
          if(i>ncy)ykn(i) = -float(nny-i+1)*pi2/yl
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

c ----------- choose correct sign so gravity waves
c             propagate out of the domain

        sgn = 1.0
        if(iocean == 1) sgn = -1.0
        if(ibcu==1) then
          do iy=1,nny
            do ix=1,nnxp2
              if(xks(ix,iy) <= 0.) then
                wavexy(ix,iy) = 0.0
              else
                wavexy(ix,iy) = sgn*sqrt(xnn/xks(ix,iy))
              endif
            enddo
          enddo
        endif

c -------------------- set length scale for SGS model

        if(iz_space == 0) then

c ------------- uniform vertical spacing

          dx32 = dx*3./2.
          dy32 = dy*3./2.
          dsl  = (abs(dx32*dy32*dzw(1)))**(1./3.)
          dslg = dsl

c -------

          do iz=0,nnzp1
            dsl_z(iz) = dslg
          enddo
          if(l_root)  write(6,2000) dsl
          if(l_debug) write(nprt,2000) dsl
          if(l_root) write(6, 4033) (iz,dsl_z(iz),iz=0,nnz+1)

c ------------- variable vertical spacing

        else

c ----------- just estimate dsl for average spacing

          dx32 = dx*3./2.
          dy32 = dy*3./2.

          dsl_max = (abs(dx32*dy32*dzw(0)))**pow3
          do iz=0,nnzp1
            dsl_z(iz) = (abs(dx32*dy32*dzw(iz)))**pow3
            if(dsl_z(iz) > dsl_max) dsl_max = dsl_z(iz)
          enddo
          dsl  = dsl_max
          dslg = dsl
          if(l_root) write(6, 4033) (iz,dsl_z(iz),iz=0,nnz+1)
        endif

c --------------------- create dsl array for easy indexing in comp1

        if(l_root) write(6, 4033) (iz,dsl_z(iz),iz=0,nnz+1)

        gridr = 1.0
        sml_eg = smal_e*gridr

c ------------------- set stokes velocity for atmos/oceanic flow

        call stokesv

c --------- special for moving 3d surface

        if(i3d_surf == 1) then
          call surf_setup
        elseif(i3d_surf == 2) then
          call surf_setup_gabls
        elseif(i3d_surf == 3) then
          call surf_setup_tank
        elseif(i3d_surf == 4) then
          call surf_setup_real
        endif

c --------- can add a time factor so as to skip into any part of
c           the specified geostrophic arrays. time factor in seconds

        t_factor = 7200.0

c -------------- for wind tunnel simulations
c                account for box height

        t_ref      = 0.0

        dpdx_con_z = rhoa*ustar_tank**2/zl

c -------------------- specify cooling rate and initial
c                      temperature even for restarts
c                      from runs with DM

        c_rate   = 0.375/3600.0
        t_surf_i = 287.99999

c -------------------- do not look for zi below zi_min

        zi_min = 30.0
        if(iocean == 1) zi_min = -5.0
        iz_min = 1
        do iz=1,nnz-1
          if(zz(iz) < zi_min .and.
     +        zz(iz+1) >= zi_min) iz_min = iz
        enddo
        if(l_root) then
          write(6,9000) zi_min, iz_min
        endif

c -------- for curvy bottom set the distribution
c          of specified heat flux. customize for
c          particular use

        call set_lower_bc

 9998   continue
        return
c --------------------------- format statements
    6   format(///,' DATA FROM RESTART FILE AT STEP =',I5,
     +      ' U_* = ',e15.6,' TS = ',e15.6,' Q_* = ',e15.6,///)
  510   format(' RESTART ***** CASE WITH : ******',/)
  520   format(' WT = ',e12.4,', U_* = ',e12.4,', L = ',e12.4,
     +      ', DTDZ FREE = ',e12.4,', ZODY = ',e12.4,/,10x,
     +      '  ZO(BTM) = ',e12.4,', CDBTM = ',e12.4,
     +      ', UG = ',e12.4)
    1   format(10x,' NNX = ',i5,',  NNY = ',i5,
     +      ',  NNZ = ',i5,/,10x,' SFC SMLT = ',i1,
     +      ',  FILTER = ',i1,
     +      ',  ITI = ',i6,',  ITMAX = ',i6,/,10x,
     +      ' IUPWIND = ',i1,',  BUYNCY = ',i1,
     +      ',  ITCUT = ',i1,/,10x,
     +      ' DT = ',e15.6,',  ZO = ',e15.6,',  TS = ',e15.6,
     +      ',  SUBSD = ',i1,/,
     +      10x,' BRCLICITY = ',i1,',  METHOD = ',i1,',  IOCEAN = ',i1,
     +      ',  IVIS = ',i1)
 4033   format(' Eddy viscosity length scale for variable z meshes',/,
     +      ' iz ',5x, ' dsl ',/,(i3,5x,e15.6))
 2000   format(10x,' DSL = ',e15.6)
 2100   format(10x,' DSL coarse = ',e15.6)
 9000   format(' Search for zi above the height = ',e15.6,/,
     +      ' iz_min = ',i5)
      end
      subroutine set_lower_bc

c -------- for curvy lower boundary set
c          distribution of heat flux or other
c          scalar sources and sinks

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats

        do iscl = 1,nscl

          if(isfc(iscl) == 0) then ! given qstar compute surface temp

            do j=iys,iye
              do i=1,nnx
                t3tbc(i,j,iscl)  = qstar(iscl)
                t_grnd(i,j,iscl) = tsfcc(iscl)
              enddo
            enddo

          else ! given surface temp compute qstar

            do j=iys,iye
              do i=1,nnx
                t3tbc(i,j,iscl)  = -99999.0
                t_grnd(i,j,iscl) = tsfcc(iscl)
              enddo
            enddo

          endif

        enddo

        return
      end
      subroutine nblnk(word)
        parameter (nmax=304)
        character wordt*304, word*(*)
        nchar = len(word)
        if(nchar > nmax) then
          write(6,6000) nchar,nmax
 6000     format(' TROUBLE, IN SR. NBLNK : NCHAR = ',i6,
     +        ' EXCEEDS NMAX = ',i6)
          stop
        endif
        jj = 0
        do j=1,nchar
          if(word(j:j) /= ' ') then
            jj = jj + 1
            wordt(jj:jj) = word(j:j)
          endif
          word(j:j) = ' '
        enddo
        do j=1,jj
          word(j:j) = wordt(j:j)
        enddo

        return
      end
      subroutine blnk(word)
        character word*(*)
        nchar = len(word)
        do j=1,nchar
          word(j:j) = ' '
        enddo

        return
      end
      subroutine iso(it)

CHECK new eddy viscosity ivis0 = 0 for testing

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        include 'mpif.h'

        real sfk(1:nnz)

c ---- get isotropy factor and scale it to match at the matching
c      height. uses boundary conditions from lower and upper.

        do iz=1,nnz
c        dfac(iz) = 1.0
          dfac(iz) = 0.0
          sfk(iz)  = 0.0
        enddo
        do iz=izs,ize
          dfac(iz) = 1.0
        enddo

c ------ set nmatch equal to fraction of initial zi in sr. random

        nmatch = 48
        do i=0,numprocs-1,ncpu_s
          if(nmatch >= iz_s(i) .and.
     +        nmatch <= iz_e(i)) myid_newvis = i
        enddo

        do iz=izs,min(ize,nmatch)
          izp1 = iz + 1
          izm1 = iz - 1

c ---- get fluctuating strains

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
     +            v_mn(iz) - v_mn(izp1))*dzu_i(izp1) +
     +            wy(i,j,iz)))**2
              sfk(iz) = sfk(iz) + 2.0*(s11 + s22 + s33) +
     +            s12 + s13 + s23
            enddo
          enddo
          sfk(iz) = sfk(iz)*fnxy
        enddo
        call mpi_sum_z(sfk,i_root,myid,nnz,1)

        do iz=izs,min(ize,nmatch)
          izp1 = iz + 1
          izm1 = iz - 1

          sfk(iz) = sqrt(sfk(iz))
          smk = sqrt((u_mn(izp1)-u_mn(iz))**2 +
     +        (v_mn(izp1)-v_mn(iz))**2)*abs(dzu_i(izp1))
          if(sfk(iz) <= 0. .and. smk <= 0.) then
            dfac(iz) = 1.0
          else
            dfac(iz) = sfk(iz)/(sfk(iz) + smk)
          endif
c     if(l_root) write(6,6001) iz,sfk(iz),smk,dfac(iz)
 6001     format(' iz = ',i3,' sfk = ',e15.6,
     +        ' smk = ',e15.6,' dfac = ',e15.6)
        enddo

c ---- rescale ratio to give unity at match height
c      and if nested grid match value at upper boundary
c      of coarser grid

        if(myid == myid_newvis) then
          dfacm = dfac(nmatch)
        endif

        call mpi_bcast(dfacm,1,mpi_real8,
     +      myid_newvis,mpi_comm_world,ierr)

        do iz=izs,min(ize,nmatch)
          dfac(iz) = dfac(iz)/dfacm
          dfac(iz) = amax1(dfac(iz), 0.1)
          dfac(iz) = amin1(dfac(iz), 1.0)
        enddo

c --------- gather dfac on all processes for printing and use in tke_vis
c           use reduce and divide by number of slab cpus

        call mpi_sum_z(dfac,i_root,myid,nnz,1)
        fncpu_s = 1.0/float(ncpu_s)
        do iz=1,nnz
          dfac(iz) = dfac(iz)*fncpu_s
        enddo

 6000   format(' in sr. iso, nmatch = ',i3,/,
     +      ' ivis = ',i3,'iz',5x,'dfac',/,(i3,1x,e15.6))
c     write(nprt,3001) (iz,dfac(iz),iz=1,nnz)
 3001   format(' iz ',5x,' dfac ',/,(i5,e15.6))
        return
      end
      subroutine surfvis(it)

c FIX for transformed grids

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        include 'mpif.h'
        real xkvis(nnx,iys:iye), alwk(nnx,iys:iye)

        real send(3), buf(3)

        xksurf = 0.0
        viscon = 0.0
        vise   = 0.0

c ----------- only root process(es) compute

        if(iss == 0) then

c     ck = 0.1
c     csmag = 0.18
c     xkmax  = dzdz/dt/5.
          iz   = 1
          izm1 = iz - 1
          izp1 = iz + 1
c     xkmax  = dzu(izp1)*dzu(izp1)/(5.0*dt)
          dz_i = dzu_i(izp1)
          if(iocean == 1) then
            call sufto(it)
          else
            call suft(it)
          endif
          if(qstar(1) == 0.) then
            zeta = 0.0
          else
            zeta = abs(z(1))/amonin
          endif
          if(ismlt == 1) then
            call busngr(zeta,phim,phis,psim,psis)
          else
            call fzol(zeta,phim,phis,psim,psis)
          endif
          viscon = vk*abs(z(1))/(utau*phim)
          vise   = utau*vk*abs(z(1))/phim

c ---- get special value at z1 to match with surface layer

          uws = 0.0
          vws = 0.0
          do iy=iys,iye
            do ix=1,nnx
              uws = uws + 0.5*(u(ix,iy,iz)-u_mn(iz) +
     +            u(ix,iy,izp1) - u_mn(izp1))*(w(ix,iy,iz)-w_mn(iz))
              vws = vws + 0.5*(v(ix,iy,iz)-v_mn(iz) +
     +            v(ix,iy,izp1) - v_mn(izp1))*(w(ix,iy,iz)-w_mn(iz))
            enddo
          enddo
          uws = uws*fnxy
          vws = vws*fnxy

c ---- get average fluctuating eddy viscsoity

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
              xkvis(ix,iy)=ck*alwk(ix,iy)*sqrt(e(ix,iy,iz))*dfac(1)
            enddo
          enddo

c ---- get average viscosity

          xkavg = 0.0
          do iy=iys,iye
            do ix=1,nnx
              xkavg = xkavg + xkvis(ix,iy)
            enddo
          enddo
          xkavg = xkavg*fnxy

          buf(1) = uws
          buf(2) = vws
          buf(3) = xkavg
          call mpi_sum_xy(buf,myid,iss,ise,3)

          uws   = buf(1)
          vws   = buf(2)
          xkavg = buf(3)

          xkz1 = vise - sqrt(uws**2 + vws**2)*viscon
          xksurf =  xkz1 - xkavg
          xksurf = amax1(xksurf,0.0)
          xksurf = amin1(xksurf,vise)
 6000     format(' dfac = ',e12.4,' xkavg = ',e12.4,' xkz1 = ',e12.4,/,
     +        ' vise = ',e12.4,' xksurf = ',e12.4)

        endif

c ---------- broadcast values to other processes

        send(1) = xksurf
        send(2) = viscon
        send(3) = vise

        call mpi_bcast(send,3,mpi_real8,
     +      i_root,mpi_comm_world,ierr)

        xksurf = send(1)
        viscon = send(2)
        vise   = send(3)

        return
      end
      subroutine comp1(istep,it)

c ----- 3-order runge-kutta time stepping and monotone scalar fluxes in x,y,z.
c       designed to use mpi in x & y directions.
c       variables advanced are (u,v,w,e,theta)/J

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        include 'mpif.h'
        integer istatus(mpi_status_size)

        parameter(js = 7, ns = 3, nstat = js + ns*nscl)
        real stat(1:nnz,nstat)

c ------ temp arrays to hold rhs from step n-1 and
c        field variables from step n

        real urhs(nnx,iys:iye,izs:ize),
     +      vrhs(nnx,iys:iye,izs:ize),
     +      wrhs(nnx,iys:iye,izs:ize),
     +      erhs(nnx,iys:iye,izs:ize),
     +      trhs(nnx,iys:iye,nscl,izs:ize)

c ------------ build rhs with u/J from last time step
c              when the routine is entered (u,v,w,e,theta) are
c              free of (1/J)

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

        do k=izs,ize
          do l=1,nscl
            do j=iys,iye
              do i=1,nnx
                trhs(i,j,l,k)=t(i,j,l,k)*gm(i,j,k,1)+dtzeta*r4(i,j,l,k)
              enddo
            enddo
          enddo
        enddo

c --------- get viscosity and rhs of (e,u,v,w)-equations
c           at next step

        if (iDNS == 1) then
          call dns_vis
          call rhs_uvw_dns(istep)
        else
          call tke_vis(istep)
          call rhs_uvw(istep)
        endif

c -------- evaluate rhs of scalar equations

        do l=1,nscl
          if (iDNS == 1) then
            call rhs_scl_dns(istep,l)
          else
            call rhs_scl(istep,l)
          endif
        enddo

c ---------- gather stat sums on root processor
c            using mpi_reduction over all processors

        if(istep == 1) then

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

c -------- end if block

        endif

c ------- save old rhs in field variables for RK-advancement
c         add new rhs for computation of flux velocities
c         at this stage field variables are (u,v,w,e,theta)/J

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

c --------- generate rhs of contravariant flux velocities equations
c           at cell faces

        call rhs_flux_vel(istep)

        return
      end
      subroutine rhs_flux_vel(istep)

c ------ get rhs of intermediate velocity (uf,vf,wf)
c        by interpolating from cell center values.
c        extrapolate from interior to get surface value for wf.
c        assumes that rhs of the (u,v,w) equations are
c        multiplied by (1/J) from comp1, ie,
c        the right sides used in the interpolation here
c        are (u/J,v/J,w/J).

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        include 'mpif.h'
        real fs(nnx,iys:iye,3), fr(nnx,iys:iye,3)
        integer istatus(mpi_status_size)
        real wfi_sum(0:nnz)

c --------- Send updated (u,v,w) (from comp1)
c           to processor below the current myid to
c           get the contravariant flux velocities.

        nb = myid - ncpu_s
        nt = myid + ncpu_s
        if(iss == 0) then
          nb = mpi_proc_null
        endif
        if(ise == numprocs-1) then
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
     +      fs(1,iys,1),nsend,mpi_real8,nb,2,
     +      fr(1,iys,1),nrecv,mpi_real8,nt,2,
     +      mpi_comm_world,istatus,ierr)
        if(ise /= numprocs-1) then
          do j=iys,iye
            do i=1,nnx
              u(i,j,ize+1)  = fr(i,j,1)
              v(i,j,ize+1)  = fr(i,j,2)
              w(i,j,ize+1)  = fr(i,j,3)
            enddo
          enddo
        endif

c ---------- get mean of wfi

        do k=0,nnz
          wfi_sum(k) = 0.0
        enddo

        do k=izs,ize
          km1 = k - 1
          kp1 = k + 1
          if(k == nnz) then

c ------------ don't worry about upper bc for wf,
c              correct value is used in psolver

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
     +              u(i,j,k)*gm(i,j,k,2)  +
     +              v(i,j,k)*gm(i,j,k,3)
                wf_kp1     = w(i,j,kp1)*gm(i,j,kp1,4)  +
     +              u(i,j,kp1)*gm(i,j,kp1,2)  +
     +              v(i,j,kp1)*gm(i,j,kp1,3)
                wfi(i,j,k) = wf_k*weitm(k) + wf_kp1*weit(k)
              enddo
            enddo
          endif

          do j=iys,iye
            do i=1,nnx
              wfi_sum(k) = wfi_sum(k) + wfi(i,j,k)
            enddo
          enddo
          wfi_sum(k) = wfi_sum(k)*fnxy

c --------- linear extrapolation to get wfi surface value

          if(k == 1) then
            w1 = -dzw(k)/(dzw(k) + dzw(kp1))
            w2 = 1.0 - w1
            do j=iys,iye
              do i=1,nnx
                ufi(i,j,km1) = u(i,j,k)
                vfi(i,j,km1) = v(i,j,k)
                wf_k         = w(i,j,k)*gm(i,j,k,4)  +
     +              u(i,j,k)*gm(i,j,k,2)  +
     +              v(i,j,k)*gm(i,j,k,3)
                wf_kp1       = w(i,j,kp1)*gm(i,j,kp1,4)  +
     +              u(i,j,kp1)*gm(i,j,kp1,2)  +
     +              v(i,j,kp1)*gm(i,j,kp1,3)
                wfi(i,j,km1) = wf_k*w2 + wf_kp1*w1
                wfi_sum(km1) = wfi_sum(km1) + wfi(i,j,km1)
              enddo
            enddo
            wfi_sum(km1) = wfi_sum(km1)*fnxy
          endif

        enddo

        call mpi_sum_z(wfi_sum(0),i_root,myid,nnz+1,1)

c -------  take out mean here

        do k=izs,ize
          do j=iys,iye
            do i=1,nnx
              wfi(i,j,k) = wfi(i,j,k) - wfi_sum(k)
            enddo
          enddo
        enddo

        if(izs == 1) then
          k = 0
          do j=iys,iye
            do i=1,nnx
              wfi(i,j,k) = wfi(i,j,k) - wfi_sum(k)
            enddo
          enddo
        endif
        return
      end
      subroutine rhs_uvw(istep)

c ---------- get right hand sides of (u,v,w) equations
c            for pencil size (nnx, iys:iye, izs:ize).
c            rhs are (r1,r2,r3)/J, but the field
c            variables used to build rhs are (u,v,w,e,theta).

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats

        real fntd(nnx,iys:iye,izs:ize)
        real fnt1(nnx,iys:iye), fnt2(nnx,iys:iye), fnt3(nnx,iys:iye)
        real fnt4(nnx,iys:iye), fnt5(nnx,iys:iye)
        real tauij(nnx,iys:iye,6,2), sij(nnx,iys:iye,6)
        real pm_z(1:nnz)

        if(iskew == 0) then
          upfac = 0.5
        else
          upfac = 1.0
        endif

c --------- special first grid point arrays, to avoid repeating

        if(izs == 1) then
          k   = 1
          kp1 = k + 1
          km1 = k - 1
          do j=iys,iye
            do i=1,nnx
              fnt4(i,j) = weit(k)*gm(i,j,kp1,1) + weitm(k)*gm(i,j,k,1)
              fnt4(i,j) = 2.0*(weit(k)*vis_m(i,j,kp1) +
     +            weitm(k)*vis_m(i,j,k))*fnt4(i,j)
              fnt5(i,j)=weit(km1)*gm(i,j,k,1)+weitm(km1)*gm(i,j,km1,1)
            enddo
          enddo
        endif

        do k=izs,ize

          km1 = k - 1
          kp1 = k + 1

c ------------ geostrophic wind

c     u_geo = -fcor*vg(k)
c     v_geo =  fcor*(ug(k)-ugal)

c ----------- for wind tunnel type simulations, account for box height

          u_geo = dpdx_con_z
          v_geo = 0.0

          do j=iys,iye
            do i=1,nnx

c ----------- boussinesq approximation for buoyancy (no hydro static piece)

              w_buy = batag*(t(i,j,1,k) - t00)

c ------------ coriolis, vertical and horizontal components

              u_cor     =  fcor*v(i,j,k) - fcor_h*w(i,j,k)
              v_cor     = -fcor*u(i,j,k)
              w_cor     =  fcor_h*u(i,j,k)

c ------------ build rhs times 1/J

              r1(i,j,k) = (u_cor + u_geo)*gm(i,j,k,1)
              r2(i,j,k) = (v_cor + v_geo)*gm(i,j,k,1)
              r3(i,j,k) = (w_cor + w_buy)*gm(i,j,k,1)
            enddo
          enddo

        enddo ! end z-loop

c ---------- advection for skew symmetric form

        if(iskew == 0) then    ! skew symmetric form no grid speed

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

                u_avg_u = u(i,j,k)*weitm(k) + u(i,j,kp1)*weit(k)
                v_avg_u = v(i,j,k)*weitm(k) + v(i,j,kp1)*weit(k)
                w_avg_u = w(i,j,k)*weitm(k) + w(i,j,kp1)*weit(k)
                u_avg_l = u(i,j,k)*weit(km1) + u(i,j,km1)*weitm(km1)
                v_avg_l = v(i,j,k)*weit(km1) + v(i,j,km1)*weitm(km1)
                w_avg_l = w(i,j,k)*weit(km1) + w(i,j,km1)*weitm(km1)

                u_adv = 0.5*(uf(i,j,k)*ux(i,j,k) + vf(i,j,k)*uy(i,j,k) +
     +              (wf(i,j,k)*u_avg_u-wf(i,j,km1)*u_avg_l)*dzw_i(k))+
     +              0.25*(wf(i,j,k)*uz + wf(i,j,km1)*uzm)
                v_adv = 0.5*(uf(i,j,k)*vx(i,j,k) + vf(i,j,k)*vy(i,j,k) +
     +              (wf(i,j,k)*v_avg_u-wf(i,j,km1)*v_avg_l)*dzw_i(k))+
     +              0.25*(wf(i,j,k)*vz + wf(i,j,km1)*vzm)
                w_adv = 0.5*(uf(i,j,k)*wx(i,j,k) + vf(i,j,k)*wy(i,j,k) +
     +              (wf(i,j,k)*w_avg_u-wf(i,j,km1)*w_avg_l)*dzw_i(k))+
     +              0.25*(wf(i,j,k)*wz + wf(i,j,km1)*wzm)

c ------------ totals

                r1(i,j,k) = r1(i,j,k) - u_adv
                r2(i,j,k) = r2(i,j,k) - v_adv
                r3(i,j,k) = r3(i,j,k) - w_adv
              enddo
            enddo
          enddo ! end z-loop

        else   ! flux form with grid speed

          do k=izs,ize
            km1 = k - 1
            kp1 = k + 1
            do j=iys,iye
              do i=1,nnx

                u_avg_u = u(i,j,k)*weitm(k) + u(i,j,kp1)*weit(k)
                v_avg_u = v(i,j,k)*weitm(k) + v(i,j,kp1)*weit(k)
                w_avg_u = w(i,j,k)*weitm(k) + w(i,j,kp1)*weit(k)
                u_avg_l = u(i,j,k)*weit(km1) + u(i,j,km1)*weitm(km1)
                v_avg_l = v(i,j,k)*weit(km1) + v(i,j,km1)*weitm(km1)
                w_avg_l = w(i,j,k)*weit(km1) + w(i,j,km1)*weitm(km1)

                wflx_u  = wf(i,j,k) - gm(i,j,k,5)
                wflx_l  = wf(i,j,km1) - gm(i,j,km1,5)
                u_adv   = (wflx_u*u_avg_u - wflx_l*u_avg_l)*dzw_i(k)
                v_adv   = (wflx_u*v_avg_u - wflx_l*v_avg_l)*dzw_i(k)
                w_adv   = (wflx_u*w_avg_u - wflx_l*w_avg_l)*dzw_i(k)

c ------------ totals

                r1(i,j,k) = r1(i,j,k) - u_adv
                r2(i,j,k) = r2(i,j,k) - v_adv
                r3(i,j,k) = r3(i,j,k) - w_adv
              enddo
            enddo
          enddo ! end z-loop

        endif ! end of skew check

c ----------- SGS cell center x-contributions + resolved fluxes
c             advective terms can be written in skew symmetric
c             or flux form. treat first point special if at lower boundary,
c             average boundary and upper cell tau flux

        do k=izs,ize
          if(k == 1) then
            call strain_u(sij,k)
            do j=iys,iye
              do i=1,nnx
                tauij(i,j,1,2) = 0.5*(fnt4(i,j)*sij(i,j,1) -
     +              t11bc(i,j)*fnt5(i,j))
                tauij(i,j,2,2) = 0.5*(fnt4(i,j)*sij(i,j,2) -
     +              t12bc(i,j)*fnt5(i,j))
                tauij(i,j,3,2) = 0.5*(fnt4(i,j)*sij(i,j,3) -
     +              t13bc(i,j)*fnt5(i,j))
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

c ----------- SGS z-contributions

        if(izs == 1) then
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
     +              weitm(k)*vis_m(i,j,k))
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
     +            weitm(k)*gm(i,j,k,1)*gm(i,j,k,2)
              gm_xl     = weit(km1)*gm(i,j,k,1)*gm(i,j,k,2) +
     +            weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,2)
              gm_yu     = weit(k)*gm(i,j,kp1,1)*gm(i,j,kp1,3) +
     +            weitm(k)*gm(i,j,k,1)*gm(i,j,k,3)
              gm_yl     = weit(km1)*gm(i,j,k,1)*gm(i,j,k,3) +
     +            weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,3)
              tau_up    = tauij(i,j,1,2)*gm_xu + tauij(i,j,2,2)*gm_yu +
     +            tauij(i,j,3,2)
              tau_lo    = tauij(i,j,1,1)*gm_xl + tauij(i,j,2,1)*gm_yl +
     +            tauij(i,j,3,1)
              uw_sum    = uw_sum + tau_up
              r1(i,j,k) = r1(i,j,k) + (tau_up - tau_lo)*dzw_i(k)
              tau_up    = tauij(i,j,2,2)*gm_xu + tauij(i,j,4,2)*gm_yu +
     +            tauij(i,j,5,2)
              tau_lo    = tauij(i,j,2,1)*gm_xl + tauij(i,j,4,1)*gm_yl +
     +            tauij(i,j,5,1)
              vw_sum    = vw_sum + tau_up
              r2(i,j,k) = r2(i,j,k) + (tau_up - tau_lo)*dzw_i(k)
              tau_up    = tauij(i,j,3,2)*gm_xu + tauij(i,j,5,2)*gm_yu +
     +            tauij(i,j,6,2)
              tau_lo    = tauij(i,j,3,1)*gm_xl + tauij(i,j,5,1)*gm_yl +
     +            tauij(i,j,6,1)
              ww_sum    = ww_sum + tau_up
              r3(i,j,k) = r3(i,j,k) + (tau_up - tau_lo)*dzw_i(k)
            enddo
          enddo

c -------- save SGS fluxes for printout

          if(istep == 1) then
            tr_tau(k) = 0.0
            do j=iys,iye
              do i=1,nnx
                ufluc    = (u(i,j,kp1) - uxym(kp1))*weit(k) +
     +              (u(i,j,k) - uxym(k))*weitm(k)
                vfluc    = (v(i,j,kp1) - vxym(kp1))*weit(k) +
     +              (v(i,j,k) - vxym(k))*weitm(k)
                tr_tau(k) = tr_tau(k) -
     +              tauij(i,j,3,2)*ufluc - tauij(i,j,5,2)*vfluc
              enddo
            enddo
            uwsb(k)   = -uw_sum*fnxy
            vwsb(k)   = -vw_sum*fnxy
            wwsb(k)   = -ww_sum*fnxy
            tr_tau(k) = tr_tau(k)*fnxy
          endif

c ---------- end z loop

        enddo

c ---------- SGS fluxes tau_12, tau_22, tau_23 that depend on
c            y-derivatives at cell centers. special treatment
c            of first point at lower boundary. assumes strain
c            at k = 1 is at upper face, see strain_12,_22,_23

        do k=izs,ize
          call strain_12(sij(1,iys,1),k)
          if(k == 1) then
            do j=iys,iye
              do i=1,nnx
                tau_c       = 0.5*(fnt4(i,j)*sij(i,j,1) -
     +              t12bc(i,j)*fnt5(i,j))
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
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
        do k=izs,ize
          do j=iys,iye
            do i=1,nnx
              r1(i,j,k)   = r1(i,j,k) + fntd(i,j,k)
            enddo
          enddo
        enddo

        do k=izs,ize
          call strain_22(sij(1,iys,1),k)
          if(k == 1) then
            do j=iys,iye
              do i=1,nnx
                tau_c       = 0.5*(fnt4(i,j)*sij(i,j,1) -
     +              t22bc(i,j)*fnt5(i,j))
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
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
        do k=izs,ize
          do j=iys,iye
            do i=1,nnx
              r2(i,j,k)   = r2(i,j,k) + fntd(i,j,k)
            enddo
          enddo
        enddo

        do k=izs,ize
          call strain_23(sij(1,iys,1),k)
          if(k == 1) then
            do j=iys,iye
              do i=1,nnx
                tau_c       = 0.5*(fnt4(i,j)*sij(i,j,1) -
     +              t23bc(i,j)*fnt5(i,j))
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
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
        do k=izs,ize
          do j=iys,iye
            do i=1,nnx
              r3(i,j,k) = r3(i,j,k) + fntd(i,j,k)
            enddo
          enddo
        enddo

c --------- get average of r3 as estimate of mean pressure gradient
c           from turbulence and add it to right hand sides of all equations

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

        do k=izs,ize
          do j=iys,iye
            do i=1,nnx
              r3(i,j,k) = r3(i,j,k) - pm_z(k)
            enddo
          enddo
        enddo

c ------- integrate p-gradient to get values at wf points

        pm(0) = p_surf
        do k=1,nnz
          pm(k) = pm(k-1) + pm_z(k)*dzw(k)
        enddo

        if(no_pgrad == 1) then

c --------- add mean p-gradients to right hand sides
c           of u and v equations using non-conservative form
c           to avoid specifying the surface pressure

          do k=izs,ize
            do j=iys,iye
              do i=1,nnx
                r1(i,j,k) = r1(i,j,k) - pm_z(k)*gm(i,j,k,1)*gm(i,j,k,2)
                r2(i,j,k) = r2(i,j,k) - pm_z(k)*gm(i,j,k,1)*gm(i,j,k,3)
              enddo
            enddo
          enddo

        endif ! check for no mean dpdx in (u,v)

c --------- leave with (r1,r2,r3)/J

        return
      end
      subroutine rhs_uvw_dns(istep)

c ---------- get right hand sides of (u,v,w) equations
c            for pencil size (nnx, iys:iye, izs:ize).
c            rhs are (r1,r2,r3)/J, but the field
c            variables used to build rhs are (u,v,w,e,theta).

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats

        real fntd(nnx,iys:iye,izs:ize)
        real fnt1(nnx,iys:iye), fnt2(nnx,iys:iye), fnt3(nnx,iys:iye)
        real fnt4(nnx,iys:iye), fnt5(nnx,iys:iye)
        real tauij(nnx,iys:iye,6,2), sij(nnx,iys:iye,6)
        real pm_z(1:nnz)
        real sfc_flx(2)

        if(iskew == 0) then
          upfac = 0.5
        else
          upfac = 1.0
        endif
        sfc_flx = 0.0

c --------- special first grid point arrays, to avoid repeating

        if(izs == 1) then
          k   = 1
          kp1 = k + 1
          km1 = k - 1
          do j=iys,iye
            do i=1,nnx
              fnt4(i,j) = weit(k)*gm(i,j,kp1,1) + weitm(k)*gm(i,j,k,1)
              fnt4(i,j) = 2.0*(weit(k)*vis_m(i,j,kp1) +
     +            weitm(k)*vis_m(i,j,k))*fnt4(i,j)
              fnt5(i,j)=weit(km1)*gm(i,j,k,1)+weitm(km1)*gm(i,j,km1,1)
            enddo
          enddo
        endif

        do k=izs,ize

          km1 = k - 1
          kp1 = k + 1

c ------------ geostrophic wind

c     u_geo = -fcor*vg(k)
c     v_geo =  fcor*(ug(k)-ugal)

c ----------- for wind tunnel type simulations, account for box height

          u_geo = dpdx_con_z
          v_geo = 0.0

          do j=iys,iye
            do i=1,nnx

c ----------- boussinesq approximation for buoyancy (no hydro static piece)

              w_buy = batag*(t(i,j,1,k) - t00)

c ------------ coriolis, vertical and horizontal components

              u_cor     =  fcor*v(i,j,k) - fcor_h*w(i,j,k)
              v_cor     = -fcor*u(i,j,k)
              w_cor     =  fcor_h*u(i,j,k)

c ------------ build rhs times 1/J

              r1(i,j,k) = (u_cor + u_geo)*gm(i,j,k,1)
              r2(i,j,k) = (v_cor + v_geo)*gm(i,j,k,1)
              r3(i,j,k) = (w_cor + w_buy)*gm(i,j,k,1)
            enddo
          enddo

        enddo ! end z-loop

c ---------- advection for skew symmetric form

        if(iskew == 0) then    ! skew symmetric form no grid speed

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

                u_avg_u = u(i,j,k)*weitm(k) + u(i,j,kp1)*weit(k)
                v_avg_u = v(i,j,k)*weitm(k) + v(i,j,kp1)*weit(k)
                w_avg_u = w(i,j,k)*weitm(k) + w(i,j,kp1)*weit(k)
                u_avg_l = u(i,j,k)*weit(km1) + u(i,j,km1)*weitm(km1)
                v_avg_l = v(i,j,k)*weit(km1) + v(i,j,km1)*weitm(km1)
                w_avg_l = w(i,j,k)*weit(km1) + w(i,j,km1)*weitm(km1)

                u_adv = 0.5*(uf(i,j,k)*ux(i,j,k) + vf(i,j,k)*uy(i,j,k) +
     +              (wf(i,j,k)*u_avg_u-wf(i,j,km1)*u_avg_l)*dzw_i(k))+
     +              0.25*(wf(i,j,k)*uz + wf(i,j,km1)*uzm)
                v_adv = 0.5*(uf(i,j,k)*vx(i,j,k) + vf(i,j,k)*vy(i,j,k) +
     +              (wf(i,j,k)*v_avg_u-wf(i,j,km1)*v_avg_l)*dzw_i(k))+
     +              0.25*(wf(i,j,k)*vz + wf(i,j,km1)*vzm)
                w_adv = 0.5*(uf(i,j,k)*wx(i,j,k) + vf(i,j,k)*wy(i,j,k) +
     +              (wf(i,j,k)*w_avg_u-wf(i,j,km1)*w_avg_l)*dzw_i(k))+
     +              0.25*(wf(i,j,k)*wz + wf(i,j,km1)*wzm)

c ------------ totals

                r1(i,j,k) = r1(i,j,k) - u_adv
                r2(i,j,k) = r2(i,j,k) - v_adv
                r3(i,j,k) = r3(i,j,k) - w_adv
              enddo
            enddo
          enddo ! end z-loop

        else   ! flux form with grid speed

          do k=izs,ize
            km1 = k - 1
            kp1 = k + 1
            do j=iys,iye
              do i=1,nnx

                u_avg_u = u(i,j,k)*weitm(k) + u(i,j,kp1)*weit(k)
                v_avg_u = v(i,j,k)*weitm(k) + v(i,j,kp1)*weit(k)
                w_avg_u = w(i,j,k)*weitm(k) + w(i,j,kp1)*weit(k)
                u_avg_l = u(i,j,k)*weit(km1) + u(i,j,km1)*weitm(km1)
                v_avg_l = v(i,j,k)*weit(km1) + v(i,j,km1)*weitm(km1)
                w_avg_l = w(i,j,k)*weit(km1) + w(i,j,km1)*weitm(km1)

                wflx_u  = wf(i,j,k) - gm(i,j,k,5)
                wflx_l  = wf(i,j,km1) - gm(i,j,km1,5)
                u_adv   = (wflx_u*u_avg_u - wflx_l*u_avg_l)*dzw_i(k)
                v_adv   = (wflx_u*v_avg_u - wflx_l*v_avg_l)*dzw_i(k)
                w_adv   = (wflx_u*w_avg_u - wflx_l*w_avg_l)*dzw_i(k)

c ------------ totals

                r1(i,j,k) = r1(i,j,k) - u_adv
                r2(i,j,k) = r2(i,j,k) - v_adv
                r3(i,j,k) = r3(i,j,k) - w_adv
              enddo
            enddo
          enddo ! end z-loop

        endif ! end of skew check

c ----------- SGS cell center x-contributions + resolved fluxes
c             advective terms can be written in skew symmetric
c             or flux form. treat first point special if at lower boundary,
c             average boundary and upper cell tau flux

        do k=izs,ize
          call strain_c_dns(sij,k)
          do j=iys,iye
            do i=1,nnx
              tauij(i,j,1,2) = 2.0*vis_m(i,j,k)*sij(i,j,1)*gm(i,j,k,1)
              tauij(i,j,2,2) = 2.0*vis_m(i,j,k)*sij(i,j,2)*gm(i,j,k,1)
              tauij(i,j,3,2) = 2.0*vis_m(i,j,k)*sij(i,j,3)*gm(i,j,k,1)
            enddo
          enddo
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

c ----------- SGS z-contributions

        call strain_u(sij,izs-1)  !NOTE: strain_u does not need dns revision -- doesn't treat k=1 specially
        do n=1,6
          do j=iys,iye
            do i=1,nnx
              fac   = 2.0*(weit(izs-1)*vis_m(i,j,izs) +
     +            weitm(izs-1)*vis_m(i,j,izs-1))
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
     +              weitm(k)*vis_m(i,j,k))
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
     +            weitm(k)*gm(i,j,k,1)*gm(i,j,k,2)
              gm_xl     = weit(km1)*gm(i,j,k,1)*gm(i,j,k,2) +
     +            weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,2)
              gm_yu     = weit(k)*gm(i,j,kp1,1)*gm(i,j,kp1,3) +
     +            weitm(k)*gm(i,j,k,1)*gm(i,j,k,3)
              gm_yl     = weit(km1)*gm(i,j,k,1)*gm(i,j,k,3) +
     +            weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,3)
              tau_up    = tauij(i,j,1,2)*gm_xu + tauij(i,j,2,2)*gm_yu +
     +            tauij(i,j,3,2)
              tau_lo    = tauij(i,j,1,1)*gm_xl + tauij(i,j,2,1)*gm_yl +
     +            tauij(i,j,3,1)
              uw_sum    = uw_sum + tau_up
              r1(i,j,k) = r1(i,j,k) + (tau_up - tau_lo)*dzw_i(k)

              if (k==1) then
                sfc_flx(1) = sfc_flx(1) + tau_lo
              endif

              tau_up    = tauij(i,j,2,2)*gm_xu + tauij(i,j,4,2)*gm_yu +
     +            tauij(i,j,5,2)
              tau_lo    = tauij(i,j,2,1)*gm_xl + tauij(i,j,4,1)*gm_yl +
     +            tauij(i,j,5,1)
              vw_sum    = vw_sum + tau_up

              if (k==1) then
                sfc_flx(2) = sfc_flx(2) + tau_lo
              endif

              r2(i,j,k) = r2(i,j,k) + (tau_up - tau_lo)*dzw_i(k)
              tau_up    = tauij(i,j,3,2)*gm_xu + tauij(i,j,5,2)*gm_yu +
     +            tauij(i,j,6,2)
              tau_lo    = tauij(i,j,3,1)*gm_xl + tauij(i,j,5,1)*gm_yl +
     +            tauij(i,j,6,1)
              ww_sum    = ww_sum + tau_up
              r3(i,j,k) = r3(i,j,k) + (tau_up - tau_lo)*dzw_i(k)

            enddo
          enddo

c -------- save SGS fluxes for printout

          if(istep == 1) then
            tr_tau(k) = 0.0
            do j=iys,iye
              do i=1,nnx
                ufluc    = (u(i,j,kp1) - uxym(kp1))*weit(k) +
     +              (u(i,j,k) - uxym(k))*weitm(k)
                vfluc    = (v(i,j,kp1) - vxym(kp1))*weit(k) +
     +              (v(i,j,k) - vxym(k))*weitm(k)
                tr_tau(k) = tr_tau(k) -
     +              tauij(i,j,3,2)*ufluc - tauij(i,j,5,2)*vfluc
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
          endif

c ---------- end z loop

        enddo

c ---------- SGS fluxes tau_12, tau_22, tau_23 that depend on
c            y-derivatives at cell centers. special treatment
c            of first point at lower boundary. assumes strain
c            at k = 1 is at upper face, see strain_12,_22,_23

        do k=izs,ize
          call strain_12_dns(sij(1,iys,1),k)
          do j=iys,iye
            do i=1,nnx
              tau_c       = 2.0*gm(i,j,k,1)*vis_m(i,j,k)*sij(i,j,1)
              fntd(i,j,k) = tau_c - upfac*vf(i,j,k)*u(i,j,k)
            enddo
          enddo
        enddo
        call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
        do k=izs,ize
          do j=iys,iye
            do i=1,nnx
              r1(i,j,k)   = r1(i,j,k) + fntd(i,j,k)
            enddo
          enddo
        enddo

        do k=izs,ize
          call strain_22_dns(sij(1,iys,1),k)

          do j=iys,iye
            do i=1,nnx
              tau_c       = 2.0*gm(i,j,k,1)*vis_m(i,j,k)*sij(i,j,1)
              fntd(i,j,k) = tau_c - upfac*vf(i,j,k)*v(i,j,k)
            enddo
          enddo
        enddo
        call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
        do k=izs,ize
          do j=iys,iye
            do i=1,nnx
              r2(i,j,k)   = r2(i,j,k) + fntd(i,j,k)
            enddo
          enddo
        enddo

        do k=izs,ize
          call strain_23_dns(sij(1,iys,1),k)
          do j=iys,iye
            do i=1,nnx
              tau_c       = 2.0*gm(i,j,k,1)*vis_m(i,j,k)*sij(i,j,1)
              fntd(i,j,k) = tau_c - upfac*vf(i,j,k)*w(i,j,k)
            enddo
          enddo
        enddo
        call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
        do k=izs,ize
          do j=iys,iye
            do i=1,nnx
              r3(i,j,k) = r3(i,j,k) + fntd(i,j,k)
            enddo
          enddo
        enddo

c --------- get average of r3 as estimate of mean pressure gradient
c           from turbulence and add it to right hand sides of all equations

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

        do k=izs,ize
          do j=iys,iye
            do i=1,nnx
              r3(i,j,k) = r3(i,j,k) - pm_z(k)
            enddo
          enddo
        enddo

c ------- integrate p-gradient to get values at wf points

        pm(0) = p_surf
        do k=1,nnz
          pm(k) = pm(k-1) + pm_z(k)*dzw(k)
        enddo

        if(no_pgrad == 1) then

c --------- add mean p-gradients to right hand sides
c           of u and v equations using non-conservative form
c           to avoid specifying the surface pressure

          do k=izs,ize
            do j=iys,iye
              do i=1,nnx
                r1(i,j,k) = r1(i,j,k) - pm_z(k)*gm(i,j,k,1)*gm(i,j,k,2)
                r2(i,j,k) = r2(i,j,k) - pm_z(k)*gm(i,j,k,1)*gm(i,j,k,3)
              enddo
            enddo
          enddo

        endif ! check for no mean dpdx in (u,v)

c --------- leave with (r1,r2,r3)/J

        return
      end
      subroutine rhs_uvw_old(istep)

c ---------- get right hand sides of (u,v,w) equations
c            for pencil size (nnx, iys:iye, izs:ize).
c            rhs are (r1,r2,r3)/J, but the field
c            variables used to build rhs are (u,v,w,e,theta).

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats

        real fntd(nnx,iys:iye,izs:ize)
        real fnt1(nnx,iys:iye), fnt2(nnx,iys:iye), fnt3(nnx,iys:iye)
        real tauij(nnx,iys:iye,6,2), sij(nnx,iys:iye,6)
        real r3_sum(1:nnz)

        if(iskew == 0) then
          upfac = 0.5
        else
          upfac = 1.0
        endif

        do k=izs,ize

          km1 = k - 1
          kp1 = k + 1

c ------------ geostrophic wind

          u_geo = -fcor*vg(k)
          v_geo =  fcor*(ug(k)-ugal)

c ----------- for wind tunnel type simulations

c     u_geo = dpdx_con
c     v_geo = 0.0

          do j=iys,iye
            do i=1,nnx

c ------------ buoyancy (with hydrostatic part)
c              will correct at the bottom
c              of the routine

c        w_buy = batag*(t(i,j,1,k) - t_mn(k,1))
              w_buy = batag*t(i,j,1,k)

c ------------ coriolis, vertical and horizontal components

              u_cor     =  fcor*v(i,j,k) - fcor_h*w(i,j,k)
              v_cor     = -fcor*u(i,j,k)
              w_cor     =  fcor_h*u(i,j,k)

c ------------ build rhs times 1/J

              r1(i,j,k) = (u_cor + u_geo)*gm(i,j,k,1)
              r2(i,j,k) = (v_cor + v_geo)*gm(i,j,k,1)
              r3(i,j,k) = (w_cor + w_buy)*gm(i,j,k,1)
            enddo
          enddo

        enddo ! end z-loop

c ---------- advection for skew symmetric form, fix for moving grids

        if(iskew == 0) then    ! skew symmetric form no grid speed

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

                u_avg_u = u(i,j,k)*weitm(k) + u(i,j,kp1)*weit(k)
                v_avg_u = v(i,j,k)*weitm(k) + v(i,j,kp1)*weit(k)
                w_avg_u = w(i,j,k)*weitm(k) + w(i,j,kp1)*weit(k)
                u_avg_l = u(i,j,k)*weit(km1) + u(i,j,km1)*weitm(km1)
                v_avg_l = v(i,j,k)*weit(km1) + v(i,j,km1)*weitm(km1)
                w_avg_l = w(i,j,k)*weit(km1) + w(i,j,km1)*weitm(km1)

                u_adv = 0.5*(uf(i,j,k)*ux(i,j,k) + vf(i,j,k)*uy(i,j,k) +
     +              (wf(i,j,k)*u_avg_u-wf(i,j,km1)*u_avg_l)*dzw_i(k))+
     +              0.25*(wf(i,j,k)*uz + wf(i,j,km1)*uzm)
                v_adv = 0.5*(uf(i,j,k)*vx(i,j,k) + vf(i,j,k)*vy(i,j,k) +
     +              (wf(i,j,k)*v_avg_u-wf(i,j,km1)*v_avg_l)*dzw_i(k))+
     +              0.25*(wf(i,j,k)*vz + wf(i,j,km1)*vzm)
                w_adv = 0.5*(uf(i,j,k)*wx(i,j,k) + vf(i,j,k)*wy(i,j,k) +
     +              (wf(i,j,k)*w_avg_u-wf(i,j,km1)*w_avg_l)*dzw_i(k))+
     +              0.25*(wf(i,j,k)*wz + wf(i,j,km1)*wzm)

c ------------ totals

                r1(i,j,k) = r1(i,j,k) - u_adv
                r2(i,j,k) = r2(i,j,k) - v_adv
                r3(i,j,k) = r3(i,j,k) - w_adv
              enddo
            enddo
          enddo ! end z-loop

        else   ! flux form with grid speed

          do k=izs,ize
            km1 = k - 1
            kp1 = k + 1
            do j=iys,iye
              do i=1,nnx

                u_avg_u = u(i,j,k)*weitm(k) + u(i,j,kp1)*weit(k)
                v_avg_u = v(i,j,k)*weitm(k) + v(i,j,kp1)*weit(k)
                w_avg_u = w(i,j,k)*weitm(k) + w(i,j,kp1)*weit(k)
                u_avg_l = u(i,j,k)*weit(km1) + u(i,j,km1)*weitm(km1)
                v_avg_l = v(i,j,k)*weit(km1) + v(i,j,km1)*weitm(km1)
                w_avg_l = w(i,j,k)*weit(km1) + w(i,j,km1)*weitm(km1)

                wflx_u  = wf(i,j,k) - gm(i,j,k,5)
                wflx_l  = wf(i,j,km1) - gm(i,j,km1,5)
                u_adv   = (wflx_u*u_avg_u - wflx_l*u_avg_l)*dzw_i(k)
                v_adv   = (wflx_u*v_avg_u - wflx_l*v_avg_l)*dzw_i(k)
                w_adv   = (wflx_u*w_avg_u - wflx_l*w_avg_l)*dzw_i(k)

c ------------ totals

                r1(i,j,k) = r1(i,j,k) - u_adv
                r2(i,j,k) = r2(i,j,k) - v_adv
                r3(i,j,k) = r3(i,j,k) - w_adv
              enddo
            enddo
          enddo ! end z-loop

        endif ! end of skew check

c ----------- SGS cell center x-contributions + resolved fluxes
c             advective terms written in skew symmetric form

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

c ----------- SGS z-contributions

        if(izs == 1) then
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
     +              weitm(k)*vis_m(i,j,k))
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
     +            weitm(k)*gm(i,j,k,1)*gm(i,j,k,2)
              gm_xl     = weit(km1)*gm(i,j,k,1)*gm(i,j,k,2) +
     +            weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,2)
              gm_yu     = weit(k)*gm(i,j,kp1,1)*gm(i,j,kp1,3) +
     +            weitm(k)*gm(i,j,k,1)*gm(i,j,k,3)
              gm_yl     = weit(km1)*gm(i,j,k,1)*gm(i,j,k,3) +
     +            weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,3)
              tau_up    = tauij(i,j,1,2)*gm_xu + tauij(i,j,2,2)*gm_yu +
     +            tauij(i,j,3,2)
              tau_lo    = tauij(i,j,1,1)*gm_xl + tauij(i,j,2,1)*gm_yl +
     +            tauij(i,j,3,1)
              uw_sum    = uw_sum + tau_up
              r1(i,j,k) = r1(i,j,k) + (tau_up - tau_lo)*dzw_i(k)
              tau_up    = tauij(i,j,2,2)*gm_xu + tauij(i,j,4,2)*gm_yu +
     +            tauij(i,j,5,2)
              tau_lo    = tauij(i,j,2,1)*gm_xl + tauij(i,j,4,1)*gm_yl +
     +            tauij(i,j,5,1)
              vw_sum    = vw_sum + tau_up
              r2(i,j,k) = r2(i,j,k) + (tau_up - tau_lo)*dzw_i(k)
              tau_up    = tauij(i,j,3,2)*gm_xu + tauij(i,j,5,2)*gm_yu +
     +            tauij(i,j,6,2)
              tau_lo    = tauij(i,j,3,1)*gm_xl + tauij(i,j,5,1)*gm_yl +
     +            tauij(i,j,6,1)
              ww_sum    = ww_sum + tau_up
              r3(i,j,k) = r3(i,j,k) + (tau_up - tau_lo)*dzw_i(k)
            enddo
          enddo

c -------- save SGS fluxes for printout

          if(istep == 1) then
            tr_tau(k) = 0.0
            do j=iys,iye
              do i=1,nnx
                ufluc    = (u(i,j,kp1) - uxym(kp1))*weit(k) +
     +              (u(i,j,k) - uxym(k))*weitm(k)
                vfluc    = (v(i,j,kp1) - vxym(kp1))*weit(k) +
     +              (v(i,j,k) - vxym(k))*weitm(k)
                tr_tau(k) = tr_tau(k) -
     +              tauij(i,j,3,2)*ufluc - tauij(i,j,5,2)*vfluc
              enddo
            enddo
            uwsb(k)   = -uw_sum*fnxy
            vwsb(k)   = -vw_sum*fnxy
            wwsb(k)   = -ww_sum*fnxy
            tr_tau(k) = tr_tau(k)*fnxy
          endif

c ---------- end z loop

        enddo

c ---------- SGS fluxes tau_12, tau_22, tau_23 that depend on
c            y-derivatives at cell centers

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
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
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
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
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
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)

c --------- add a constant to r3 to account for hydrostatic pressure
c           and for periodic sidewalls. we have r3/J therefore the
c           constant to be added  = - <J*r3/J>/J

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

        call mpi_sum_z(r3_sum,i_root,myid,nnz,1)

c ------- make sure <r3/J> = 0 and set r3 = 0 at top

        do k=izs,ize
          if(k == nnz) then
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

c --------- leave with (r1,r2,r3)/J

        return
      end
      subroutine strain_u(sij_u,k)

c ---------- get strain matrix at the upper face at the k-level

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        real sij_u(nnx,iys:iye,6)

        km1  = k - 1
        kp1  = k + 1
        w1   = weit(k)
        w1m  = weitm(k)
        w2   = weit(km1)
        w2m  = weitm(km1)

c -------- at upper face

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

        return
      end
      subroutine strain_c_dns(sij_c,k)

c ---------- get strain matrix at the cell center
c ---------- for DNS just don't treak k = 1 any different

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        real sij_c(nnx,iys:iye,6)

        km1  = k - 1
        kp1  = k + 1
        w1   = weit(k)
        w1m  = weitm(k)
        w2   = weit(km1)
        w2m  = weitm(km1)

c -------- at cell center

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
     +          vx(i,j,k)  + vzet*gm_x)
            sij_c(i,j,3) = 0.5*(uzet*gm_z + wx(i,j,k)  + wzet*gm_x)
            sij_c(i,j,4) = vy(i,j,k) + vzet*gm_y
            sij_c(i,j,5) = 0.5*(vzet*gm_z + wy(i,j,k) + wzet*gm_y)
            sij_c(i,j,6) = wzet*gm_z
          enddo
        enddo

!      endif

        return
      end
      subroutine strain_c(sij_c,k)

c ---------- get strain matrix at the cell center
c            special if at first level off the wall

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        real sij_c(nnx,iys:iye,6)

        km1  = k - 1
        kp1  = k + 1
        w1   = weit(k)
        w1m  = weitm(k)
        w2   = weit(km1)
        w2m  = weitm(km1)

        if(k == 1) then

c ------- vertical gradients at upper face
Cnew      now computed at cell centers using surface bcs

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
     +            vx(i,j,k)  + vzet*gm_x)
              sij_c(i,j,3) = 0.5*(uzet*gm_z + wx(i,j,k)  + wzet*gm_x)
              sij_c(i,j,4) = vy(i,j,k) + vzet*gm_y
              sij_c(i,j,5) = 0.5*(vzet*gm_z + wy(i,j,k) + wzet*gm_y)
              sij_c(i,j,6) = wzet*gm_z
            enddo
          enddo

        else

c -------- at cell center

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
     +            vx(i,j,k)  + vzet*gm_x)
              sij_c(i,j,3) = 0.5*(uzet*gm_z + wx(i,j,k)  + wzet*gm_x)
              sij_c(i,j,4) = vy(i,j,k) + vzet*gm_y
              sij_c(i,j,5) = 0.5*(vzet*gm_z + wy(i,j,k) + wzet*gm_y)
              sij_c(i,j,6) = wzet*gm_z
            enddo
          enddo

        endif

        return
      end
      subroutine strain_12_dns(sij_y,k)

c ---------- get strain s12
c            for special evaluation of eta derivatives
c            using matrix transpose

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        real sij_y(nnx,iys:iye)

        km1  = k - 1
        kp1  = k + 1
        w1   = weit(k)
        w1m  = weitm(k)
        w2   = weit(km1)
        w2m  = weitm(km1)

c -------- at cell center

        do j=iys,iye
          do i=1,nnx
            uup   = w1*u(i,j,kp1) + w1m*u(i,j,k)
            ulo   = w2*u(i,j,k) + w2m*u(i,j,km1)
            uzet  = (uup - ulo)*dzw_i(k)
            vup   = w1*v(i,j,kp1) + w1m*v(i,j,k)
            vlo   = w2*v(i,j,k) + w2m*v(i,j,km1)
            vzet  = (vup - vlo)*dzw_i(k)
            sij_y(i,j) = 0.5*(uy(i,j,k) + uzet*gm(i,j,k,3) +
     +          vx(i,j,k) + vzet*gm(i,j,k,2))
          enddo
        enddo

!      endif

        return
      end
      subroutine strain_12(sij_y,k)

c ---------- get strain s12
c            for special evaluation of eta derivatives
c            using matrix transpose

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        real sij_y(nnx,iys:iye)

        km1  = k - 1
        kp1  = k + 1
        w1   = weit(k)
        w1m  = weitm(k)
        w2   = weit(km1)
        w2m  = weitm(km1)

        if(k == 1) then

c ------ at upper face

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

        else

c -------- at cell center

          do j=iys,iye
            do i=1,nnx
              uup   = w1*u(i,j,kp1) + w1m*u(i,j,k)
              ulo   = w2*u(i,j,k) + w2m*u(i,j,km1)
              uzet  = (uup - ulo)*dzw_i(k)
              vup   = w1*v(i,j,kp1) + w1m*v(i,j,k)
              vlo   = w2*v(i,j,k) + w2m*v(i,j,km1)
              vzet  = (vup - vlo)*dzw_i(k)
              sij_y(i,j) = 0.5*(uy(i,j,k) + uzet*gm(i,j,k,3) +
     +            vx(i,j,k) + vzet*gm(i,j,k,2))
            enddo
          enddo

        endif

        return
      end
      subroutine strain_22_dns(sij_y,k)

c ---------- get strain s22 at cell center
c            for special evaluation of eta derivatives
c            using matrix transpose

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        real sij_y(nnx,iys:iye)

        km1  = k - 1
        kp1  = k + 1
        w1   = weit(k)
        w1m  = weitm(k)
        w2   = weit(km1)
        w2m  = weitm(km1)

c -------- at cell center

        do j=iys,iye
          do i=1,nnx
            vup   = w1*v(i,j,kp1) + w1m*v(i,j,k)
            vlo   = w2*v(i,j,k) + w2m*v(i,j,km1)
            vzet  = (vup - vlo)*dzw_i(k)
            sij_y(i,j) = vy(i,j,k) + vzet*gm(i,j,k,3)
          enddo
        enddo

!      endif

        return
      end
      subroutine strain_22(sij_y,k)

c ---------- get strain s22 at cell center
c            for special evaluation of eta derivatives
c            using matrix transpose
        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        real sij_y(nnx,iys:iye)

        km1  = k - 1
        kp1  = k + 1
        w1   = weit(k)
        w1m  = weitm(k)
        w2   = weit(km1)
        w2m  = weitm(km1)

        if(k == 1) then

c -------at upper face

          do j=iys,iye
            do i=1,nnx
              vzet       = (v(i,j,kp1) - v(i,j,k))*dzu_i(kp1)
              veta       = w1*vy(i,j,kp1) + w1m*vy(i,j,k)
              gm_y       = w1*gm(i,j,kp1,3) + w1m*gm(i,j,k,3)
              sij_y(i,j) = veta + vzet*gm_y
            enddo
          enddo
        else

c -------- at cell center

          do j=iys,iye
            do i=1,nnx
              vup   = w1*v(i,j,kp1) + w1m*v(i,j,k)
              vlo   = w2*v(i,j,k) + w2m*v(i,j,km1)
              vzet  = (vup - vlo)*dzw_i(k)
              sij_y(i,j) = vy(i,j,k) + vzet*gm(i,j,k,3)
            enddo
          enddo

        endif

        return
      end
      subroutine strain_23_dns(sij_y,k)

c ---------- get strain s23 at cell center
c            for special evaluation of eta derivatives
c            using matrix transpose
        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        real sij_y(nnx,iys:iye)

        km1  = k - 1
        kp1  = k + 1
        w1   = weit(k)
        w1m  = weitm(k)
        w2   = weit(km1)
        w2m  = weitm(km1)

c -------- at cell center

        do j=iys,iye
          do i=1,nnx
            vup   = w1*v(i,j,kp1) + w1m*v(i,j,k)
            vlo   = w2*v(i,j,k) + w2m*v(i,j,km1)
            vzet  = (vup - vlo)*dzw_i(k)
            wup   = w1*w(i,j,kp1) + w1m*w(i,j,k)
            wlo   = w2*w(i,j,k) + w2m*w(i,j,km1)
            wzet  = (wup - wlo)*dzw_i(k)
            sij_y(i,j) = 0.5*(vzet*gm(i,j,k,4) +
     +          wy(i,j,k) + wzet*gm(i,j,k,3))
          enddo
        enddo

!      endif

        return
      end
      subroutine strain_23(sij_y,k)

c ---------- get strain s23 at cell center
c            for special evaluation of eta derivatives
c            using matrix transpose

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        real sij_y(nnx,iys:iye)

        km1  = k - 1
        kp1  = k + 1
        w1   = weit(k)
        w1m  = weitm(k)
        w2   = weit(km1)
        w2m  = weitm(km1)

        if(k == 1) then

c -------- at upper face

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

c -------- at cell center

          do j=iys,iye
            do i=1,nnx
              vup   = w1*v(i,j,kp1) + w1m*v(i,j,k)
              vlo   = w2*v(i,j,k) + w2m*v(i,j,km1)
              vzet  = (vup - vlo)*dzw_i(k)
              wup   = w1*w(i,j,kp1) + w1m*w(i,j,k)
              wlo   = w2*w(i,j,k) + w2m*w(i,j,km1)
              wzet  = (wup - wlo)*dzw_i(k)
              sij_y(i,j) = 0.5*(vzet*gm(i,j,k,4) +
     +            wy(i,j,k) + wzet*gm(i,j,k,3))
            enddo
          enddo

        endif

        return
      end
      subroutine rhs_scl_dns(istep,iscl)

c ------ get right hand side of scalar equation (iscl)
c        monotone scalar fluxes only in z
c        for pencil size (nnx, iys:iye, izs:ize)
c        care is taken so that if monotone is on then
c        conservative horizontal flux form is used!
c        uses flux conservation form for advection always

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats

        real fnt1(nnx,iys:iye), fntd(nnx,iys:iye,izs:ize)
        real tx(nnx,iys:iye,izs-1:ize+1), ty(nnx,iys:iye,izs-1:ize+1)
        real flux_u(nnx,iys:iye), flux_l(nnx,iys:iye)
        real tflx(nnx,iys:iye,3,2), scl_g(nnx,iys:iye,3)
        real sfc_flx

c --------- set sign for ocean simulations that use monotone

        sgn = 1.0
        if(iocean == 1) sgn = -1.0
        if(iupwnd == 1) then
          upfac = 1.0
        else
          upfac = 0.5
        endif
        sfc_flx = 0.0

c --------- get x & y derivatives for [izs-1:ize+1]

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
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)

c -------- SGS term tflx(.,.,1) + resolved flux at cell center

        do k=izs,ize
          kp1 = k + 1
          km1 = k - 1
          call grads_c_dns(scl_g,k,tx,ty,iscl)
          do j=iys,iye
            do i=1,nnx
              tflx(i,j,1,1) = vis_s(i,j,k)*scl_g(i,j,1)*gm(i,j,k,1)
              tflx(i,j,2,1) = vis_s(i,j,k)*scl_g(i,j,2)*gm(i,j,k,1)
              fnt1(i,j)     = tflx(i,j,1,1) -
     +            upfac*uf(i,j,k)*t(i,j,iscl,k)
              fntd(i,j,k)   = tflx(i,j,2,1) -
     +            upfac*vf(i,j,k)*t(i,j,iscl,k)
            enddo
          enddo
!         endif
          call xderivp(fnt1(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
          do j=iys,iye
            do i=1,nnx
              r4(i,j,iscl,k) = fnt1(i,j)
            enddo
          enddo

c -------- save SGS fluxes for printout

          if(istep == 1) then
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

c ---------- SGS fluxes tflx(.,.,2) that depend on
c            y-derivative at cell centers

        call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
        do k=izs,ize
          do j=iys,iye
            do i=1,nnx
              r4(i,j,iscl,k) = r4(i,j,iscl,k) + fntd(i,j,k)
            enddo
          enddo
        enddo

c ------- resolved advection for skew symmetric form
c         for no monotone

        if(iupwnd /= 1) then

          do k=izs,ize
            do j=iys,iye
              do i=1,nnx
                r4(i,j,iscl,k) = r4(i,j,iscl,k) -
     +              upfac*(uf(i,j,k)*tx(i,j,k) + vf(i,j,k)*ty(i,j,k))
              enddo
            enddo
          enddo
        endif ! end of iupwnd check for horizontal advection

c ----------- SGS z-contributions

        call grads_u(scl_g,izs-1,tx,ty,iscl)
        do n=1,3
          do j=iys,iye
            do i=1,nnx
              fac   = weit(izs-1)*vis_s(i,j,izs) +
     +            weitm(izs-1)*vis_s(i,j,izs-1)
              tflx(i,j,n,2) = fac*scl_g(i,j,n)
            enddo
          enddo
        enddo
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
     +            weitm(k)*gm(i,j,k,1)*gm(i,j,k,2)
              gm_xl      = weit(km1)*gm(i,j,k,1)*gm(i,j,k,2) +
     +            weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,2)
              gm_yu      = weit(k)*gm(i,j,kp1,1)*gm(i,j,kp1,3) +
     +            weitm(k)*gm(i,j,k,1)*gm(i,j,k,3)
              gm_yl      = weit(km1)*gm(i,j,k,1)*gm(i,j,k,3) +
     +            weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,3)
              tflx_up    = tflx(i,j,1,2)*gm_xu + tflx(i,j,2,2)*gm_yu +
     +            tflx(i,j,3,2)
              tflx_lo    = tflx(i,j,1,1)*gm_xl + tflx(i,j,2,1)*gm_yl +
     +            tflx(i,j,3,1)
              r4(i,j,iscl,k) = r4(i,j,iscl,k) +
     +            (tflx_up - tflx_lo)*dzw_i(k)
              if (k==1) then
                sfc_flx = sfc_flx + tflx_lo
              endif
            enddo
          enddo

c -------- save SGS fluxes for printout, gather sums on exit

          if(istep == 1) then
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
          endif

c ---------- end z loop

        enddo

c ---------- resolved vertical advection

        if(iupwnd /= 1) then

c --------- skew symmetric advective form for
c           vertical flux = 0.5(wf*dt/dz + d/dz(wf*t))
c           no grid speed

          do k=izs,ize
            km1 = k - 1
            kp1 = k + 1
            do j=iys,iye
              do i=1,nnx
                tzm     = (t(i,j,iscl,k) - t(i,j,iscl,km1))*dzu_i(k)
                tz      = (t(i,j,iscl,kp1) - t(i,j,iscl,k))*dzu_i(kp1)
                theta_u = weitm(k)*t(i,j,iscl,k) +
     +              weit(k)*t(i,j,iscl,kp1)
                theta_l = weit(km1)*t(i,j,iscl,k) +
     +              weitm(km1)*t(i,j,iscl,km1)
                t_adv   = 0.5*(wf(i,j,k)*theta_u -
     +              wf(i,j,km1)*theta_l)*dzw_i(k) +
     +              0.25*(wf(i,j,k)*tz + wf(i,j,km1)*tzm)
                r4(i,j,iscl,k) = r4(i,j,iscl,k) - t_adv
              enddo
            enddo
          enddo

        else

c ----------- z-direction special with monotone

          do k=izs,ize
            km2 = k - 2
            km1 = k - 1
            kp1 = k + 1
            kp2 = k + 2
            if(k == 1) then
              do j=iys,iye
                do i=1,nnx
                  wflx_l      = wf(i,j,km1) - gm(i,j,km1,5)
                  wflx_u      = wf(i,j,k)   - gm(i,j,k,5)
                  flux_l(i,j) = sgn*0.5*wflx_l*
     +                (t(i,j,iscl,km1)+t(i,j,iscl,k))
                  flux_u(i,j) =
     +                amax1(sgn*wflx_u,0.)*(t(i,j,iscl,k) +
     +                rlim(t(i,j,iscl,kp1),t(i,j,iscl,k),
     +                t(i,j,iscl,km1))) +
     +                amin1(sgn*wflx_u,0.)*(t(i,j,iscl,kp1) +
     +                rlim(t(i,j,iscl,k),t(i,j,iscl,kp1),
     +                t(i,j,iscl,kp2)))
                enddo
              enddo
            elseif(k == nnz) then
              do j=iys,iye
                do i=1,nnx
                  wflx_l      = wf(i,j,km1) - gm(i,j,km1,5)
                  wflx_u      = wf(i,j,k) - gm(i,j,k,5)
                  flux_u(i,j) = sgn*0.5*wflx_u*
     +                (t(i,j,iscl,kp1)+t(i,j,iscl,k))
                  flux_l(i,j) =
     +                amax1(sgn*wflx_l,0.)*(t(i,j,iscl,km1) +
     +                rlim(t(i,j,iscl,k),t(i,j,iscl,km1),
     +                t(i,j,iscl,km2))) +
     +                amin1(sgn*wflx_l,0.)*(t(i,j,iscl,k) +
     +                rlim(t(i,j,iscl,km1),t(i,j,iscl,k),
     +                t(i,j,iscl,kp1)))
                enddo
              enddo
            else
              do j=iys,iye
                do i=1,nnx
                  wflx_l      = wf(i,j,km1) - gm(i,j,km1,5)
                  wflx_u      = wf(i,j,k) - gm(i,j,k,5)
                  flux_u(i,j) =
     +                amax1(sgn*wflx_u,0.)*(t(i,j,iscl,k) +
     +                rlim(t(i,j,iscl,kp1),t(i,j,iscl,k),
     +                t(i,j,iscl,km1))) +
     +                amin1(sgn*wflx_u,0.)*(t(i,j,iscl,kp1) +
     +                rlim(t(i,j,iscl,k),t(i,j,iscl,kp1),
     +                t(i,j,iscl,kp2)))
                  flux_l(i,j) =
     +                amax1(sgn*wflx_l,0.)*(t(i,j,iscl,km1) +
     +                rlim(t(i,j,iscl,k),t(i,j,iscl,km1),
     +                t(i,j,iscl,km2))) +
     +                amin1(sgn*wflx_l,0.)*(t(i,j,iscl,k) +
     +                rlim(t(i,j,iscl,km1),t(i,j,iscl,k),
     +                t(i,j,iscl,kp1)))
                enddo
              enddo
            endif

c ---------- sum vertical monotone flux

            do j=iys,iye
              do i=1,nnx
                r4(i,j,iscl,k) = r4(i,j,iscl,k)
     +              - sgn*(flux_u(i,j) - flux_l(i,j))*dzw_i(k)
              enddo
            enddo

c ---------- end z loop

          enddo

c -------- end monotone if block

        endif

        return
      end
      subroutine rhs_scl(istep,iscl)

c ------ get right hand side of scalar equation (iscl)
c        monotone scalar fluxes only in z
c        for pencil size (nnx, iys:iye, izs:ize)
c        care is taken so that if monotone is on then
c        conservative horizontal flux form is used!
c        uses flux conservation form for advection always

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats

        real fnt1(nnx,iys:iye), fntd(nnx,iys:iye,izs:ize)
        real tx(nnx,iys:iye,izs-1:ize+1), ty(nnx,iys:iye,izs-1:ize+1)
        real flux_u(nnx,iys:iye), flux_l(nnx,iys:iye)
        real tflx(nnx,iys:iye,3,2), scl_g(nnx,iys:iye,3)

c --------- set sign for ocean simulations that use monotone

        sgn = 1.0
        if(iocean == 1) sgn = -1.0
        if(iupwnd == 1) then
          upfac = 1.0
        else
          upfac = 0.5
        endif

c --------- get x & y derivatives for [izs-1:ize+1]

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
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)

c -------- SGS term tflx(.,.,1) + resolved flux at cell center
c          special treatment of first grid point off the wall,
c          average upper cell and boundary value to get cell
c          center value.

        do k=izs,ize
          kp1 = k + 1
          km1 = k - 1
          if(k == 1) then
            call grads_u(scl_g,k,tx,ty,iscl)
            do j=iys,iye
              do i=1,nnx
                gm_u  = weit(k)*gm(i,j,kp1,1) + weitm(k)*gm(i,j,k,1)
                gm_l  = weit(km1)*gm(i,j,k,1) + weitm(km1)*gm(i,j,km1,1)
                fac   = weit(k)*vis_s(i,j,kp1) + weitm(k)*vis_s(i,j,k)
                tflx(i,j,1,1) = 0.5*(fac*scl_g(i,j,1)*gm_u -
     +              t1tbc(i,j,iscl)*gm_l)
                tflx(i,j,2,1) = 0.5*(fac*scl_g(i,j,2)*gm_u -
     +              t2tbc(i,j,iscl)*gm_l)
                fnt1(i,j)     = tflx(i,j,1,1) -
     +              upfac*uf(i,j,k)*t(i,j,iscl,k)
                fntd(i,j,k)   = tflx(i,j,2,1) -
     +              upfac*vf(i,j,k)*t(i,j,iscl,k)
              enddo
            enddo
          else
            call grads_c(scl_g,k,tx,ty,iscl)
            do j=iys,iye
              do i=1,nnx
                tflx(i,j,1,1) = vis_s(i,j,k)*scl_g(i,j,1)*gm(i,j,k,1)
                tflx(i,j,2,1) = vis_s(i,j,k)*scl_g(i,j,2)*gm(i,j,k,1)
                fnt1(i,j)     = tflx(i,j,1,1) -
     +              upfac*uf(i,j,k)*t(i,j,iscl,k)
                fntd(i,j,k)   = tflx(i,j,2,1) -
     +              upfac*vf(i,j,k)*t(i,j,iscl,k)
              enddo
            enddo
          endif
          call xderivp(fnt1(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
          do j=iys,iye
            do i=1,nnx
              r4(i,j,iscl,k) = fnt1(i,j)
            enddo
          enddo

c -------- save SGS fluxes for printout

          if(istep == 1) then
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

c ---------- SGS fluxes tflx(.,.,2) that depend on
c            y-derivative at cell centers

        call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
        do k=izs,ize
          do j=iys,iye
            do i=1,nnx
              r4(i,j,iscl,k) = r4(i,j,iscl,k) + fntd(i,j,k)
            enddo
          enddo
        enddo

c ------- resolved advection for skew symmetric form
c         for no monotone

        if(iupwnd /= 1) then

          do k=izs,ize
            do j=iys,iye
              do i=1,nnx
                r4(i,j,iscl,k) = r4(i,j,iscl,k) -
     +              upfac*(uf(i,j,k)*tx(i,j,k) + vf(i,j,k)*ty(i,j,k))
              enddo
            enddo
          enddo
        endif ! end of iupwnd check for horizontal advection

c ----------- SGS z-contributions

        if(izs == 1) then
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
     +            weitm(k)*gm(i,j,k,1)*gm(i,j,k,2)
              gm_xl      = weit(km1)*gm(i,j,k,1)*gm(i,j,k,2) +
     +            weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,2)
              gm_yu      = weit(k)*gm(i,j,kp1,1)*gm(i,j,kp1,3) +
     +            weitm(k)*gm(i,j,k,1)*gm(i,j,k,3)
              gm_yl      = weit(km1)*gm(i,j,k,1)*gm(i,j,k,3) +
     +            weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,3)
              tflx_up    = tflx(i,j,1,2)*gm_xu + tflx(i,j,2,2)*gm_yu +
     +            tflx(i,j,3,2)
              !A hack to enforce a constant flux of heat/moisture
              !through the domain
              if (iz==nnz) then
                tflx_up = wtsfc(iscl)
              endif
              tflx_lo    = tflx(i,j,1,1)*gm_xl + tflx(i,j,2,1)*gm_yl +
     +            tflx(i,j,3,1)
              r4(i,j,iscl,k) = r4(i,j,iscl,k) +
     +            (tflx_up - tflx_lo)*dzw_i(k)
            enddo
          enddo

c -------- save SGS fluxes for printout, gather sums on exit

          if(istep == 1) then
            wtsb(k,iscl) = 0.0
            do j=iys,iye
              do i=1,nnx
                wtsb(k,iscl) = wtsb(k,iscl) - tflx(i,j,3,2)
              enddo
            enddo
            wtsb(k,iscl) = wtsb(k,iscl)*fnxy
          endif

c ---------- end z loop

        enddo

c ---------- resolved vertical advection

        if(iupwnd /= 1) then

c --------- skew symmetric advective form for
c           vertical flux = 0.5(wf*dt/dz + d/dz(wf*t))
c           no grid speed

          do k=izs,ize
            km1 = k - 1
            kp1 = k + 1
            do j=iys,iye
              do i=1,nnx
                tzm     = (t(i,j,iscl,k) - t(i,j,iscl,km1))*dzu_i(k)
                tz      = (t(i,j,iscl,kp1) - t(i,j,iscl,k))*dzu_i(kp1)
                theta_u = weitm(k)*t(i,j,iscl,k) +
     +              weit(k)*t(i,j,iscl,kp1)
                theta_l = weit(km1)*t(i,j,iscl,k) +
     +              weitm(km1)*t(i,j,iscl,km1)
                t_adv   = 0.5*(wf(i,j,k)*theta_u -
     +              wf(i,j,km1)*theta_l)*dzw_i(k) +
     +              0.25*(wf(i,j,k)*tz + wf(i,j,km1)*tzm)
                r4(i,j,iscl,k) = r4(i,j,iscl,k) - t_adv
              enddo
            enddo
          enddo

        else

c ----------- z-direction special with monotone

          do k=izs,ize
            km2 = k - 2
            km1 = k - 1
            kp1 = k + 1
            kp2 = k + 2
            if(k == 1) then
              do j=iys,iye
                do i=1,nnx
                  wflx_l      = wf(i,j,km1) - gm(i,j,km1,5)
                  wflx_u      = wf(i,j,k)   - gm(i,j,k,5)
                  flux_l(i,j) = sgn*0.5*wflx_l*
     +                (t(i,j,iscl,km1)+t(i,j,iscl,k))
                  flux_u(i,j) =
     +                amax1(sgn*wflx_u,0.)*(t(i,j,iscl,k) +
     +                rlim(t(i,j,iscl,kp1),t(i,j,iscl,k),
     +                t(i,j,iscl,km1))) +
     +                amin1(sgn*wflx_u,0.)*(t(i,j,iscl,kp1) +
     +                rlim(t(i,j,iscl,k),t(i,j,iscl,kp1),
     +                t(i,j,iscl,kp2)))
                enddo
              enddo
            elseif(k == nnz) then
              do j=iys,iye
                do i=1,nnx
                  wflx_l      = wf(i,j,km1) - gm(i,j,km1,5)
                  wflx_u      = wf(i,j,k) - gm(i,j,k,5)
!                  flux_u(i,j) = sgn*0.5*wflx_u*
!     +                (t(i,j,iscl,kp1)+t(i,j,iscl,k))
                  flux_u(i,j) = wtsfc(iscl)
                  flux_l(i,j) =
     +                amax1(sgn*wflx_l,0.)*(t(i,j,iscl,km1) +
     +                rlim(t(i,j,iscl,k),t(i,j,iscl,km1),
     +                t(i,j,iscl,km2))) +
     +                amin1(sgn*wflx_l,0.)*(t(i,j,iscl,k) +
     +                rlim(t(i,j,iscl,km1),t(i,j,iscl,k),
     +                t(i,j,iscl,kp1)))
                enddo
              enddo
            else
              do j=iys,iye
                do i=1,nnx
                  wflx_l      = wf(i,j,km1) - gm(i,j,km1,5)
                  wflx_u      = wf(i,j,k) - gm(i,j,k,5)
                  flux_u(i,j) =
     +                amax1(sgn*wflx_u,0.)*(t(i,j,iscl,k) +
     +                rlim(t(i,j,iscl,kp1),t(i,j,iscl,k),
     +                t(i,j,iscl,km1))) +
     +                amin1(sgn*wflx_u,0.)*(t(i,j,iscl,kp1) +
     +                rlim(t(i,j,iscl,k),t(i,j,iscl,kp1),
     +                t(i,j,iscl,kp2)))
                  flux_l(i,j) =
     +                amax1(sgn*wflx_l,0.)*(t(i,j,iscl,km1) +
     +                rlim(t(i,j,iscl,k),t(i,j,iscl,km1),
     +                t(i,j,iscl,km2))) +
     +                amin1(sgn*wflx_l,0.)*(t(i,j,iscl,k) +
     +                rlim(t(i,j,iscl,km1),t(i,j,iscl,k),
     +                t(i,j,iscl,kp1)))
                enddo
              enddo
            endif

c ---------- sum vertical monotone flux

            do j=iys,iye
              do i=1,nnx
                r4(i,j,iscl,k) = r4(i,j,iscl,k)
     +              - sgn*(flux_u(i,j) - flux_l(i,j))*dzw_i(k)
              enddo
            enddo

c ---------- end z loop

          enddo

c -------- end monotone if block

        endif

        return
      end
      subroutine rhs_scl_old(istep,iscl)

c ------ get right hand side of scalar equation (iscl)
c        monotone scalar fluxes only in z
c        for pencil size (nnx, iys:iye, izs:ize)
c        care is taken so that if monotone is on then
c        conservative horizontal flux form is used!
c        uses flux conservation form for advection always

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats

        real fnt1(nnx,iys:iye), fntd(nnx,iys:iye,izs:ize)
        real tx(nnx,iys:iye,izs-1:ize+1), ty(nnx,iys:iye,izs-1:ize+1)
        real flux_u(nnx,iys:iye), flux_l(nnx,iys:iye)
        real tflx(nnx,iys:iye,3,2), scl_g(nnx,iys:iye,3)

c --------- set sign for ocean simulations that use monotone

        sgn = 1.0
        if(iocean == 1) sgn = -1.0
        if(iupwnd == 1) then
          upfac = 1.0
        else
          upfac = 0.5
        endif

c --------- get x & y derivatives for [izs-1:ize+1]

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
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)

c -------- SGS term tflx(.,.,1) + resolved flux at cell center

        do k=izs,ize
          call grads_c(scl_g,k,tx,ty,iscl)
          do j=iys,iye
            do i=1,nnx
              tflx(i,j,1,1) = vis_s(i,j,k)*scl_g(i,j,1)
              tflx(i,j,2,1) = vis_s(i,j,k)*scl_g(i,j,2)
              fnt1(i,j)     = gm(i,j,k,1)*tflx(i,j,1,1) -
     +            upfac*uf(i,j,k)*t(i,j,iscl,k)
              fntd(i,j,k)   = gm(i,j,k,1)*tflx(i,j,2,1) -
     +            upfac*vf(i,j,k)*t(i,j,iscl,k)
            enddo
          enddo
          call xderivp(fnt1(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
          do j=iys,iye
            do i=1,nnx
              r4(i,j,iscl,k) = fnt1(i,j)
            enddo
          enddo

c -------- save SGS fluxes for printout

          if(istep == 1) then
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

c ---------- SGS fluxes tflx(.,.,2) that depend on
c            y-derivative at cell centers

        call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
        do k=izs,ize
          do j=iys,iye
            do i=1,nnx
              r4(i,j,iscl,k) = r4(i,j,iscl,k) + fntd(i,j,k)
            enddo
          enddo
        enddo

c ------- resolved advection for skew symmetric form
c         for no monotone

        if(iupwnd /= 1) then

          do k=izs,ize
            do j=iys,iye
              do i=1,nnx
                r4(i,j,iscl,k) = r4(i,j,iscl,k) -
     +              upfac*(uf(i,j,k)*tx(i,j,k) + vf(i,j,k)*ty(i,j,k))
              enddo
            enddo
          enddo
        endif

c ----------- SGS z-contributions

        if(izs == 1) then
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
     +            weitm(k)*gm(i,j,k,1)*gm(i,j,k,2)
              gm_xl      = weit(km1)*gm(i,j,k,1)*gm(i,j,k,2) +
     +            weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,2)
              gm_yu      = weit(k)*gm(i,j,kp1,1)*gm(i,j,kp1,3) +
     +            weitm(k)*gm(i,j,k,1)*gm(i,j,k,3)
              gm_yl      = weit(km1)*gm(i,j,k,1)*gm(i,j,k,3) +
     +            weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,3)
              tflx_up    = tflx(i,j,1,2)*gm_xu + tflx(i,j,2,2)*gm_yu +
     +            tflx(i,j,3,2)
              tflx_lo    = tflx(i,j,1,1)*gm_xl + tflx(i,j,2,1)*gm_yl +
     +            tflx(i,j,3,1)
              r4(i,j,iscl,k) = r4(i,j,iscl,k) +
     +            (tflx_up - tflx_lo)*dzw_i(k)
            enddo
          enddo

c -------- save SGS fluxes for printout, gather sums on exit

          if(istep == 1) then
            wtsb(k,iscl) = 0.0
            do j=iys,iye
              do i=1,nnx
                wtsb(k,iscl) = wtsb(k,iscl) - tflx(i,j,3,2)
              enddo
            enddo
            wtsb(k,iscl) = wtsb(k,iscl)*fnxy
          endif

c ---------- end z loop

        enddo

c ---------- resolved vertical advection

        if(iupwnd /= 1) then

c --------- skew symmetric advective form for
c           vertical flux = 0.5(wf*dt/dz + d/dz(wf*t))
c           no grid speed

          do k=izs,ize
            km1 = k - 1
            kp1 = k + 1
            do j=iys,iye
              do i=1,nnx
                tzm     = (t(i,j,iscl,k) - t(i,j,iscl,km1))*dzu_i(k)
                tz      = (t(i,j,iscl,kp1) - t(i,j,iscl,k))*dzu_i(kp1)
                theta_u = weitm(k)*t(i,j,iscl,k) +
     +              weit(k)*t(i,j,iscl,kp1)
                theta_l = weit(km1)*t(i,j,iscl,k) +
     +              weitm(km1)*t(i,j,iscl,km1)
                t_adv   = 0.5*(wf(i,j,k)*theta_u -
     +              wf(i,j,km1)*theta_l)*dzw_i(k) +
     +              0.25*(wf(i,j,k)*tz + wf(i,j,km1)*tzm)
                r4(i,j,iscl,k) = r4(i,j,iscl,k) - t_adv
              enddo
            enddo
          enddo

        else

c ----------- z-direction special with monotone

          do k=izs,ize
            km2 = k - 2
            km1 = k - 1
            kp1 = k + 1
            kp2 = k + 2
            if(k == 1) then
              do j=iys,iye
                do i=1,nnx
                  wflx_l      = wf(i,j,km1) - gm(i,j,km1,5)
                  wflx_u      = wf(i,j,k)   - gm(i,j,k,5)
                  flux_l(i,j) = sgn*0.5*wflx_l*
     +                (t(i,j,iscl,km1)+t(i,j,iscl,k))
                  flux_u(i,j) =
     +                amax1(sgn*wflx_u,0.)*(t(i,j,iscl,k) +
     +                rlim(t(i,j,iscl,kp1),t(i,j,iscl,k),
     +                t(i,j,iscl,km1))) +
     +                amin1(sgn*wflx_u,0.)*(t(i,j,iscl,kp1) +
     +                rlim(t(i,j,iscl,k),t(i,j,iscl,kp1),
     +                t(i,j,iscl,kp2)))
                enddo
              enddo
            elseif(k == nnz) then
              do j=iys,iye
                do i=1,nnx
                  wflx_l      = wf(i,j,km1) - gm(i,j,km1,5)
                  wflx_u      = wf(i,j,k) - gm(i,j,k,5)
                  flux_u(i,j) = sgn*0.5*wflx_u*
     +                (t(i,j,iscl,kp1)+t(i,j,iscl,k))
                  flux_l(i,j) =
     +                amax1(sgn*wflx_l,0.)*(t(i,j,iscl,km1) +
     +                rlim(t(i,j,iscl,k),t(i,j,iscl,km1),
     +                t(i,j,iscl,km2))) +
     +                amin1(sgn*wflx_l,0.)*(t(i,j,iscl,k) +
     +                rlim(t(i,j,iscl,km1),t(i,j,iscl,k),
     +                t(i,j,iscl,kp1)))
                enddo
              enddo
            else
              do j=iys,iye
                do i=1,nnx
                  wflx_l      = wf(i,j,km1) - gm(i,j,km1,5)
                  wflx_u      = wf(i,j,k) - gm(i,j,k,5)
                  flux_u(i,j) =
     +                amax1(sgn*wflx_u,0.)*(t(i,j,iscl,k) +
     +                rlim(t(i,j,iscl,kp1),t(i,j,iscl,k),
     +                t(i,j,iscl,km1))) +
     +                amin1(sgn*wflx_u,0.)*(t(i,j,iscl,kp1) +
     +                rlim(t(i,j,iscl,k),t(i,j,iscl,kp1),
     +                t(i,j,iscl,kp2)))
                  flux_l(i,j) =
     +                amax1(sgn*wflx_l,0.)*(t(i,j,iscl,km1) +
     +                rlim(t(i,j,iscl,k),t(i,j,iscl,km1),
     +                t(i,j,iscl,km2))) +
     +                amin1(sgn*wflx_l,0.)*(t(i,j,iscl,k) +
     +                rlim(t(i,j,iscl,km1),t(i,j,iscl,k),
     +                t(i,j,iscl,kp1)))
                enddo
              enddo
            endif

c ---------- sum vertical monotone flux

            do j=iys,iye
              do i=1,nnx
                r4(i,j,iscl,k) = r4(i,j,iscl,k)
     +              - sgn*(flux_u(i,j) - flux_l(i,j))*dzw_i(k)
              enddo
            enddo

c ---------- end z loop

          enddo

c -------- end monotone if block

        endif

        return
      end
      subroutine grads_u(scl_u,k,tx,ty,iscl)

c ---------- get scalar gradient at upper face at the k-level

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        real scl_u(nnx,iys:iye,3)
        real tx(nnx,iys:iye,izs-1:ize+1), ty(nnx,iys:iye,izs-1:ize+1)

        km1  = k - 1
        kp1  = k + 1
        w1   = weit(k)
        w1m  = weitm(k)
        w2   = weit(km1)
        w2m  = weitm(km1)

c -------- at upper face

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

        return
      end
      subroutine grads_c_dns(scl_c,k,tx,ty,iscl)

c ---------- get scalar gradient at the cell center

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        real scl_c(nnx,iys:iye,3)
        real tx(nnx,iys:iye,izs-1:ize+1), ty(nnx,iys:iye,izs-1:ize+1)

        km1  = k - 1
        kp1  = k + 1
        w1   = weit(k)
        w1m  = weitm(k)
        w2   = weit(km1)
        w2m  = weitm(km1)

        if(k == 1) then

          do j=iys,iye
            do i=1,nnx
              tzet  = (t(i,j,iscl,kp1) - t(i,j,iscl,k))*dzu_i(kp1)
              scl_c(i,j,1) = tx(i,j,k) + tzet*gm(i,j,k,2)
              scl_c(i,j,2) = ty(i,j,k) + tzet*gm(i,j,k,3)
              scl_c(i,j,3) = tzet*gm(i,j,k,4)
            enddo
          enddo
        else

c -------- at cell center

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

        endif

        return
      end
      subroutine grads_c(scl_c,k,tx,ty,iscl)

c ---------- get scalar gradient at the cell center

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        real scl_c(nnx,iys:iye,3)
        real tx(nnx,iys:iye,izs-1:ize+1), ty(nnx,iys:iye,izs-1:ize+1)

        km1  = k - 1
        kp1  = k + 1
        w1   = weit(k)
        w1m  = weitm(k)
        w2   = weit(km1)
        w2m  = weitm(km1)

        if(k == 1) then

          do j=iys,iye
            do i=1,nnx
              tzet  = (t(i,j,iscl,kp1) - t(i,j,iscl,k))*dzu_i(kp1)
              scl_c(i,j,1) = tx(i,j,k) + tzet*gm(i,j,k,2)
              scl_c(i,j,2) = ty(i,j,k) + tzet*gm(i,j,k,3)
              scl_c(i,j,3) = tzet*gm(i,j,k,4)
            enddo
          enddo
        else

c -------- at cell center

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

        endif

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

c ------ get viscosity using deardorff tke model with
c        stability correction. fixes for surface layer.
c        get rhs of e-equation/J. field variables are
c        (u,v,w,e,theta) all at cell center. get viscosity
c        at cell center also.
c        "length scale" based on the volume in physical space

c        volume = gm(.,.,.,1)*d_zeta*d_eta*d_xi
c        (fletcher, vol2 1988, p. 51)

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats

        real ex(nnx,iys:iye,izs-1:ize+1), ey(nnx,iys:iye,izs-1:ize+1)
        real fnt1(nnx,iys:iye), fntd(nnx,iys:iye,izs:ize)
        real dissp(nnx,iys:iye), alk(nnx,iys:iye), sij(nnx,iys:iye,6)
        real eflx(nnx,iys:iye,3,2), e_g(nnx,iys:iye,3)
        real rj3(nnx,iys:iye)

c -------- see par statments for setting k_stab

        do k=max(izs-1,1),min(nnz,ize+1)

          kp1  = k + 1
          km1  = k - 1
          w1   = weit(k)
          w1m  = weitm(k)
          w2   = weit(km1)
          w2m  = weitm(km1)
          dslk = dsl_z(k)

          do j=iys,iye
            do i=1,nnx
              rj3(i,j) = gm(i,j,k,1)**pow3
            enddo
          enddo

          if(k > 0) dslk  = amin1(dsl_z(k),vk*abs(z(k))/csmag)
          almin = almin_c*dsl_z(k)
          if(k == 0 .or. k == nnz+1) then
            dfack = 1.0
          else
            dfack = dfac(k)
          endif
c     if(ivis == 1 .and. k .le. nmatch) then

          if(k <= k_stab) then

c --------------- no stability corrected length scales

            do j=iys,iye
              do i=1,nnx
                alk(i,j) = dslk*rj3(i,j)
              enddo
            enddo
          else
            do j=iys,iye
              do i=1,nnx
                alk(i,j) = dslk*rj3(i,j)
                tup  = w1*t(i,j,1,kp1) + w1m*t(i,j,1,k)
                tlo  = w2*t(i,j,1,k) + w2m*t(i,j,1,km1)
                stab = gm(i,j,k,4)*batag*(tup - tlo)*dzw_i(k)
                if(stab > stabmin) then
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
     +            vis_m(i,j,k)
            enddo
          enddo

c -------------- special case at bottom and top

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

          if(k >= izs .and. k <= ize) then

c -------- to avoid storing length scale get dissipation
c          and shear production with cell centered strains
c          special at lower boundary

            call strain_c(sij,k)
            if(k == 1) then
              do j=iys,iye
                do i=1,nnx
                  fac       = 2.0*(weit(k)*vis_m(i,j,k+1) +
     +                weitm(k)*vis_m(i,j,k))
                  tau11_s11=0.5*(fac*sij(i,j,1)-t11bc(i,j))*sij(i,j,1)
                  tau12_s12=0.5*(fac*sij(i,j,2)-t12bc(i,j))*sij(i,j,2)
                  tau13_s13=0.5*(fac*sij(i,j,3)-t13bc(i,j))*sij(i,j,3)
                  tau22_s22=0.5*(fac*sij(i,j,4)-t22bc(i,j))*sij(i,j,4)
                  tau23_s23=0.5*(fac*sij(i,j,5)-t23bc(i,j))*sij(i,j,5)
                  tau33_s33=0.5*(fac*sij(i,j,6)-t33bc(i,j))*sij(i,j,6)
                  fnt1(i,j) = tau11_s11 + tau22_s22 + tau33_s33 +
     +                2.0*(tau12_s12 + tau13_s13 + tau23_s23)
                enddo
              enddo
            else
              do j=iys,iye
                do i=1,nnx
                  fnt1(i,j)  = 2.0*vis_m(i,j,k)*(
     +                sij(i,j,1)**2 + sij(i,j,4)**2 + sij(i,j,6)**2 +
     +                2.0*(sij(i,j,2)**2+sij(i,j,3)**2+sij(i,j,5)**2))
                enddo
              enddo
            endif

            do j=iys,iye
              do i=1,nnx
                dissp(i,j) = (0.19 + 0.74*alk(i,j)/(dslk*rj3(i,j)))*
     +              e(i,j,k)*sqrt(e(i,j,k))/alk(i,j)
                r5(i,j,k) = (fnt1(i,j) - dissp(i,j))*gm(i,j,k,1)
              enddo
            enddo

            if(istep == 1) then
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

          endif   ! endif for k = [izs,ize]

c -------------- end z loop

        enddo

c --------- get x & y derivatives for [izs-1:ize+1]

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
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)

c -------- diffusion term eflx(.,.,1) + resolved flux at cell center

        do k=izs,ize
          call grd_e_c(e_g,k,ex,ey)
          do j=iys,iye
            do i=1,nnx
              eflx(i,j,1,1) = 2.0*vis_m(i,j,k)*e_g(i,j,1)
              eflx(i,j,2,1) = 2.0*vis_m(i,j,k)*e_g(i,j,2)
              fnt1(i,j)     = gm(i,j,k,1)*eflx(i,j,1,1) -
     +            uf(i,j,k)*e(i,j,k)
              fntd(i,j,k)   = gm(i,j,k,1)*eflx(i,j,2,1) -
     +            vf(i,j,k)*e(i,j,k)
            enddo
          enddo
          call xderivp(fnt1(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
          do j=iys,iye
            do i=1,nnx
              r5(i,j,k) = r5(i,j,k) + fnt1(i,j)
            enddo
          enddo
        enddo

c -------- diffusion term eflx(.,.,2) + resolved flux at cell center

        call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
        do k=izs,ize
          do j=iys,iye
            do i=1,nnx
              r5(i,j,k) = r5(i,j,k) + fntd(i,j,k)
            enddo
          enddo
        enddo

c --------- z-diffusion term + resolved vertical flux
c           no divergence of flux at first gridpoint

        k   = izs
        km1 = izs - 1
        kp1 = izs + 1
        if(izs == 1) then
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
     +            weitm(k)*vis_m(i,j,k)
              eflx(i,j,n,2) = 2.0*fac*e_g(i,j,n)
            enddo
          enddo
        enddo

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
     +            weitm(k)*gm(i,j,k,1)*gm(i,j,k,2)
              gm_xl     = weit(km1)*gm(i,j,k,1)*gm(i,j,k,2) +
     +            weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,2)
              gm_yu     = weit(k)*gm(i,j,kp1,1)*gm(i,j,kp1,3) +
     +            weitm(k)*gm(i,j,k,1)*gm(i,j,k,3)
              gm_yl     = weit(km1)*gm(i,j,k,1)*gm(i,j,k,3) +
     +            weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,3)
              eflx_up   = eflx(i,j,1,2)*gm_xu + eflx(i,j,2,2)*gm_yu +
     +            eflx(i,j,3,2)
              eflx_lo   = eflx(i,j,1,1)*gm_xl + eflx(i,j,2,1)*gm_yl +
     +            eflx(i,j,3,1)
              r5(i,j,k) = r5(i,j,k) + (eflx_up - eflx_lo)*dzw_i(k)
            enddo
          enddo

c ---------- resolved vertical advection, flux form
c            plus grid speed

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

c --------- get buoyancy term/J in sgs e equation

          if(k == 1) then
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
                r5(i,j,k)=r5(i,j,k)-vis_s(i,j,k)*tzet*batag*gm(i,j,k,1)
              enddo
            enddo
          endif

c -------- end z loop

        enddo

        return
      end
      subroutine tke_vis_old(istep)

c ------ get viscosity using deardorff tke model with
c        stability correction. fixes for surface layer.
c        get rhs of e-equation/J. field variables are
c        (u,v,w,e,theta) all at cell center. get viscosity
c        at cell center also.
c        "length scale" based on the volume in physical space

c        volume = gm(.,.,.,1)*d_zeta*d_eta*d_xi
c        (fletcher, vol2 1988, p. 51)

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats

        real ex(nnx,iys:iye,izs-1:ize+1), ey(nnx,iys:iye,izs-1:ize+1)
        real fnt1(nnx,iys:iye), fntd(nnx,iys:iye,izs:ize)
        real dissp(nnx,iys:iye), alk(nnx,iys:iye), sij(nnx,iys:iye,6)
        real eflx(nnx,iys:iye,3,2), e_g(nnx,iys:iye,3)

        do k=max(izs-1,1),min(nnz,ize+1)

          kp1  = k + 1
          km1  = k - 1
          w1   = weit(k)
          w1m  = weitm(k)
          w2   = weit(km1)
          w2m  = weitm(km1)
          dslk = dsl_z(k)
          if(k > 0) dslk  = amin1(dsl_z(k),vk*abs(z(k))/csmag)
          almin = almin_c*dsl_z(k)
          if(k == 0 .or. k == nnz+1) then
            dfack = 1.0
          else
            dfack = dfac(k)
          endif

          if(k <= ((3*nnz)/6)) then

c --------------- no stability corrected length scales

            do j=iys,iye
              do i=1,nnx
                alk(i,j) = dslk*gm(i,j,k,1)
              enddo
            enddo
          else
            do j=iys,iye
              do i=1,nnx
                alk(i,j) = dslk*gm(i,j,k,1)
                tup  = w1*t(i,j,1,kp1) + w1m*t(i,j,1,k)
                tlo  = w2*t(i,j,1,k) + w2m*t(i,j,1,km1)
                stab = gm(i,j,k,4)*batag*(tup - tlo)*dzw_i(k)
                if(stab > stabmin) then
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
     +            vis_m(i,j,k)
            enddo
          enddo

c -------------- special case at bottom and top

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

          if(k >= izs .and. k <= ize) then

c -------- to avoid storing length scale get dissipation
c          and shear production with cell centered strains

            call strain_c(sij,k)
            do j=iys,iye
              do i=1,nnx
                dissp(i,j) = (0.19 + 0.74*alk(i,j)/(dslk*gm(i,j,k,1)))*
     +              e(i,j,k)*sqrt(e(i,j,k))/alk(i,j)
                fnt1(i,j)  = 2.0*vis_m(i,j,k)*(
     +              sij(i,j,1)**2 + sij(i,j,4)**2 + sij(i,j,6)**2 +
     +              2.0*(sij(i,j,2)**2 + sij(i,j,3)**2 + sij(i,j,5)**2))
                r5(i,j,k) = (fnt1(i,j) - dissp(i,j))*gm(i,j,k,1)
              enddo
            enddo

            if(istep == 1) then
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

          endif   ! endif for k = [izs,ize]

c -------------- end z loop

        enddo

c --------- get x & y derivatives for [izs-1:ize+1]

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
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)

c -------- diffusion term eflx(.,.,1) + resolved flux at cell center

        do k=izs,ize
          call grd_e_c(e_g,k,ex,ey)
          do j=iys,iye
            do i=1,nnx
              eflx(i,j,1,1) = 2.0*vis_m(i,j,k)*e_g(i,j,1)
              eflx(i,j,2,1) = 2.0*vis_m(i,j,k)*e_g(i,j,2)
              fnt1(i,j)     = gm(i,j,k,1)*eflx(i,j,1,1) -
     +            uf(i,j,k)*e(i,j,k)
              fntd(i,j,k)   = gm(i,j,k,1)*eflx(i,j,2,1) -
     +            vf(i,j,k)*e(i,j,k)
            enddo
          enddo
          call xderivp(fnt1(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
          do j=iys,iye
            do i=1,nnx
              r5(i,j,k) = r5(i,j,k) + fnt1(i,j)
            enddo
          enddo
        enddo

c -------- diffusion term eflx(.,.,2) + resolved flux at cell center

        call yd_mpi(fntd(1,iys,izs),trigx(1,2),yk(1),
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
        do k=izs,ize
          do j=iys,iye
            do i=1,nnx
              r5(i,j,k) = r5(i,j,k) + fntd(i,j,k)
            enddo
          enddo
        enddo

c --------- z-diffusion term + resolved vertical flux
c           no divergence of flux at first gridpoint

        k   = izs
        km1 = izs - 1
        kp1 = izs + 1
        if(izs == 1) then
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
     +            weitm(k)*vis_m(i,j,k)
              eflx(i,j,n,2) = 2.0*fac*e_g(i,j,n)
            enddo
          enddo
        enddo

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
     +            weitm(k)*gm(i,j,k,1)*gm(i,j,k,2)
              gm_xl     = weit(km1)*gm(i,j,k,1)*gm(i,j,k,2) +
     +            weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,2)
              gm_yu     = weit(k)*gm(i,j,kp1,1)*gm(i,j,kp1,3) +
     +            weitm(k)*gm(i,j,k,1)*gm(i,j,k,3)
              gm_yl     = weit(km1)*gm(i,j,k,1)*gm(i,j,k,3) +
     +            weitm(km1)*gm(i,j,km1,1)*gm(i,j,km1,3)
              eflx_up   = eflx(i,j,1,2)*gm_xu + eflx(i,j,2,2)*gm_yu +
     +            eflx(i,j,3,2)
              eflx_lo   = eflx(i,j,1,1)*gm_xl + eflx(i,j,2,1)*gm_yl +
     +            eflx(i,j,3,1)
              r5(i,j,k) = r5(i,j,k) + (eflx_up - eflx_lo)*dzw_i(k)
            enddo
          enddo

c ---------- resolved vertical advection, flux form
c            plus grid speed

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

c --------- get buoyancy term/J in sgs e equation

          if(k == 1) then
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
                r5(i,j,k)=r5(i,j,k)-vis_s(i,j,k)*tzet*batag*gm(i,j,k,1)
              enddo
            enddo
          endif

c -------- end z loop

        enddo

        return
      end
      subroutine grd_e_u(scl_u,k,ex,ey)

c ---------- get scalar gradient at upper face at the k-level

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        real scl_u(nnx,iys:iye,3)
        real ex(nnx,iys:iye,izs-1:ize+1), ey(nnx,iys:iye,izs-1:ize+1)

        km1  = k - 1
        kp1  = k + 1
        w1   = weit(k)
        w1m  = weitm(k)
        w2   = weit(km1)
        w2m  = weitm(km1)

c -------- at upper face

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

        return
      end
      subroutine grd_e_c(scl_c,k,ex,ey)

c ---------- get scalar gradient at the cell center

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        real scl_c(nnx,iys:iye,3)
        real ex(nnx,iys:iye,izs-1:ize+1), ey(nnx,iys:iye,izs-1:ize+1)

        km1  = k - 1
        kp1  = k + 1
        w1   = weit(k)
        w1m  = weitm(k)
        w2   = weit(km1)
        w2m  = weitm(km1)

        if(k == 1) then

          do j=iys,iye
            do i=1,nnx
              ezet  = (e(i,j,kp1) - e(i,j,k))*dzu_i(kp1)
              scl_c(i,j,1) = ex(i,j,k) + ezet*gm(i,j,k,2)
              scl_c(i,j,2) = ey(i,j,k) + ezet*gm(i,j,k,3)
              scl_c(i,j,3) = ezet*gm(i,j,k,4)
            enddo
          enddo
        else

c -------- at cell center

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

        endif

        return
      end
      subroutine lower(it)

c ------ setup lower boundary condition for entire plane at (k = 1)
c        using either businger or large formulas with wind.
c        index f(.,.,2)  indicates lower.

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        real sfc_flx(2+nscl)

        k    = 1
        km1  = k - 1
        dz_i = dzu_i(k)

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

        if(iocean == 1) then ! water surface layer
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

        else  !   land surface layer

          call suft(it)
          fac = -utau**2/(windm*sqrt(u1xy**2 + v1xy**2))
          do j=iys,iye
            do i=1,nnx
              tau13m(i,j) = fac*(windm*(u(i,j,k)+ugal-u1xy)+
     +            wind(i,j)*u1xy)
              tau23m(i,j) = fac*(windm*(v(i,j,k)-v1xy)+
     +            wind(i,j)*v1xy)
            enddo
          enddo
          do iscl=1,nscl
            dnom3=t10xy(iscl)*windm
            if(dnom3 /= 0.) then
              dnom_i = 1.0/dnom3
              do j=iys,iye
                do i=1,nnx
                  taut3m(i,j,iscl) = aut3m(iscl)*
     +                (windm*(t(i,j,iscl,k) - t1xy(iscl))+
     +                wind(i,j)*(t1xy(iscl) - tsfcc(iscl)))*dnom_i
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

        endif ! end of surface layer type

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

c -------- partial sums of surface fluxes and mean scalar

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
 2345   format(' in lower 2345 uwsfc = ',e15.6,' vwsfc = ',e15.6,
     +      ' wtsfc = ',e15.6,' tsfcc = ',e15.6)

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

c ------------ initialize u, v, w, t derivatives are
c              computed in get_derv at km1
c              no need to call derivatives here since
c              wbc = 0, change for more general lower bc

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

        return
      end
      subroutine lower_free(it)

c --------------- setup lower boundary condition for free
c                 convection where each processor applies
c                 log-law at several (ix,iy) for k = 1.
c                 index f(.,.,2)  indicates lower

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        include 'mpif.h'

        real u_level1(nnx,iys:iye,2+nscl), buf(2+2*nscl)
        real sbuf(2+2*nscl,mxs:mxe,iys:iye)
        real rbuf((2+2*nscl)*nnx*(iye+1-iys))
        integer istatus(mpi_status_size)

c -------------- broadcast level 1 data everywhere

        if(iss == 0) then
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

c ------ send all of root data to other processors

        call mpi_send_root(u_level1(1,iys,1),
     +      num,myid,numprocs,ncpu_s)

c --------- every task gets their own fluxes and surface scalars

        call suft2(u_level1,it)

c --------- send surface scalars and momentum fluxes
c           back to root(s)

        if(numprocs == 1) go to 999

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

        irow_r = mod(myid,ncpu_s)
        if(myid >= ncpu_s) then
          num = (2+2*nscl)*(mxe+1-mxs)*(iye+1-iys)
          call mpi_send(sbuf(1,mxs,iys),num,mpi_real8,irow_r,1,
     +        mpi_comm_world,ierr)
        else
          do l=irow_r+ncpu_s,numprocs-1,ncpu_s
            num = (2+2*nscl)*(mx_e(l)+1-mx_s(l))*(iye+1-iys)
            call mpi_recv(rbuf(1),num,mpi_real8,l,1,
     +          mpi_comm_world,istatus,ierr)
c          call f_suft2(rbuf,maxnx,maxny,mx_s(l),mx_e(l),iys,iye,nscl,
            call f_suft2(rbuf,nnx,mx_s(l),mx_e(l),iys,iye,nscl,
     +          tau13m,tau23m,taut3m,t_grnd)
          enddo
        endif

  999   continue

c ------------ only for root row = 0
c              get sums of surface conditions
c              and set surface boundary conditions

        if(iss == 0) then

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

          call mpi_sum_xy(buf,myid,iss,ise,2+2*nscl)
          uwsfc = buf(1)*fnxy
          vwsfc = buf(2)*fnxy
          do iscl=1,nscl
            wtsfc(iscl) = buf(2+iscl)*fnxy
            tsfcc(iscl) = buf(2+nscl+iscl)*fnxy
          enddo

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

          iz   = 1
          izm1 = iz - 1
          dz_i = dzu_i(iz)

          do iy=iys,iye
            do ix=1,nnx
              ebc(ix,iy,2)=amax1(e(ix,iy,iz),sml_eg)
              wbc(ix,iy,2)= 0.0
              pbc(ix,iy,2) = 0.0
              pbc2(ix,iy,2) = 0.0

              wfbc(i,j,2)  = 0.0
              ufbc(i,j,2)  = 0.0
              vfbc(i,j,2)  = 0.0
            enddo
          enddo

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

c ------------ initialize u, v, w, t and derivatives at izm1

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

          do iscl=1,nscl
            do iy=iys,iye
              do ix=1,nnx
                t(ix,iy,iscl,izm1) = tbc(ix,iy,iscl,2)
              enddo
            enddo
          enddo

c ----- end of if block for root row

        endif

 7999   continue

        return
      end
      subroutine f_suft2(rbuf,nnx,mxs,mxe,iys,iye,nscl,
     +      tau13m,tau23m,taut3m,t_grnd)

c ------ fill surface arrays on root processors

        real rbuf(2+2*nscl,mxs:mxe,iys:iye)
        real tau13m(nnx,iys:iye), tau23m(nnx,iys:iye),
     +      taut3m(nnx,iys:iye,nscl), t_grnd(nnx,iys:iye,nscl)

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

        return
      end
      subroutine lower_curvy(it,istage)

c --------- setup lower boundary condition for case
c           with curvy lower boundary. each processor
c           applies log-law at several (i,j) for k = 1

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        include 'mpif.h'

c ------- set the number of variables to pass and
c         get back

        parameter (mvar= 9+3*nscl, jvar = 8+2*nscl)
        real u_level1(nnx,iys:iye,mvar), buf(jvar)
        real sbuf(jvar,mxs:mxe,iys:iye)
        real rbuf(jvar*nnx*(iye+1-iys))
        integer istatus(mpi_status_size)

        t_bndy = t_stage - t_zero

        if(i3d_surf == 0) then
          call l_bndy(t_bndy)
        elseif(i3d_surf == 1) then
          call get_surf(t_bndy)
        elseif(i3d_surf == 2) then
          call get_surf_gabls(t_bndy)
        elseif(i3d_surf == 3) then
          call get_surf_tank(t_bndy)
        elseif(i3d_surf == 4) then
          call get_surf_real(t_bndy)
        endif

c     call dealias_grid_speed

c ----- broadcast level 1 data to processors above root

        k   = 1
        km1 = k - 1

        if(iss == 0) then
          do j=iys,iye
            do i=1,nnx

c -------- subtract off surface movement when applying the drag law
c          plus the orbital velocity of the wave.
c          w at the boundary is close to grid speed.
c          Note for iz_space = 4 grid speed is constant across
c          the first grid cell and thus z_t at w point = z_t at
c          the boundary. Can improve the surface estimate FIX
c          and couple it to the underlying water using the
c          fuller definition

              u_level1(i,j,1) = u(i,j,k) - u_orbit(i,j)
              u_level1(i,j,2) = v(i,j,k) - v_orbit(i,j)
              u_level1(i,j,3) = w(i,j,k) - w_orbit(i,j)
              u_level1(i,j,4) = gm(i,j,k,1)
              u_level1(i,j,5) = gm(i,j,k,2)
              u_level1(i,j,6) = gm(i,j,k,3)
              u_level1(i,j,7) = u_orbit(i,j)
              u_level1(i,j,8) = v_orbit(i,j)
              u_level1(i,j,9) = gm(i,j,km1,5)
            enddo
          enddo
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

        num = nnx*(iye + 1 - iys)*mvar

c ------ send all of root data to other processors

        call mpi_send_root(u_level1(1,iys,1),
     +      num,myid,numprocs,ncpu_s)

c --------- every task gets their own fluxes and surface scalars

        do j=iys,iye
          do i=1,nnx
            ustar(i,j)    = -999999.99
c        t_grnd(i,j,1) = -999999.99
          enddo
        enddo

        call drglaw_l(u_level1,mvar,it)

c --------- send surface scalars and momentum fluxes
c           back to root(s)

        if(numprocs == 1) go to 999

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

        irow_r = mod(myid,ncpu_s)
        if(myid >= ncpu_s) then
          num = jvar*(mxe+1-mxs)*(iye+1-iys)
          call mpi_send(sbuf(1,mxs,iys),num,mpi_real8,irow_r,1,
     +        mpi_comm_world,ierr)
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

  999   continue

        if(iss == 0) then

c ------------ only for root row = 0
c              get sums of surface conditions
c              and set surface boundary conditions

          do i=1,jvar
            buf(i) = 0.0
          enddo

c ------- add SGS pieces that appear on rhs of [u,v,w] equations

          do j=iys,iye
            do i=1,nnx
              fac1   = gm(i,j,k,1)*gm(i,j,k,2)
              fac2   = gm(i,j,k,1)*gm(i,j,k,3)
              buf(1) = buf(1) + t11bc(i,j)
              buf(2) = buf(2) + t12bc(i,j)
              buf(3) = buf(3) + t13bc(i,j) +
     +            fac1*t11bc(i,j) + fac2*t12bc(i,j)
              buf(4) = buf(4) + t22bc(i,j)
              buf(5) = buf(5) + t23bc(i,j) +
     +            fac1*t12bc(i,j) + fac2*t22bc(i,j)
              buf(6) = buf(6) + t33bc(i,j) +
     +            fac1*t13bc(i,j) + fac2*t23bc(i,j)
              buf(7) = buf(7) + ustar(i,j)**2
              buf(8) = buf(8) + amol(i,j)
            enddo
          enddo
          do iscl=1,nscl
            do j=iys,iye
              do i=1,nnx
                buf(8+iscl)      = buf(8+iscl) + t3tbc(i,j,iscl)
                buf(8+nscl+iscl) = buf(8+nscl+iscl) +
     +              t_grnd(i,j,iscl)
              enddo
            enddo
          enddo

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

          dz_i = dzu_i(k)

          do j=iys,iye
            do i=1,nnx
            !ebc(i,j,2)  = amax1(e(i,j,k),sml_eg)
              ebc(i,j,2) = 0.0
              pbc(i,j,2)  = 0.0
              pbc2(i,j,2) = 0.0

              ubc(i,j,2)  = 2.0*u_orbit(i,j) - u(i,j,k)
              vbc(i,j,2)  = 2.0*v_orbit(i,j) - v(i,j,k)

Cnew ------ bc for vertical velocity

              wbc(i,j,2)  = 2.0*gm(i,j,km1,5) - w(i,j,k)

              ufbc(i,j,2)  = -uf(i,j,k)
              vfbc(i,j,2)  = -vf(i,j,k)

c ------- no flow through the surface at t = t_stage

              wfbc(i,j,2)  = gm(i,j,km1,5)
            enddo
          enddo

          do iscl=1,nscl
            do j=iys,iye
              do i=1,nnx
                tbc(i,j,iscl,2) = 2.0*t_grnd(i,j,iscl) - t(i,j,iscl,k)
              enddo
            enddo
          enddo

c ------------ initialize u, v, w, t and derivatives at km1

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

              wf(i,j,km1) = wfbc(i,j,2)
              uf(i,j,km1) = ufbc(i,j,2)
              vf(i,j,km1) = vfbc(i,j,2)
            enddo
          enddo

          do iscl=1,nscl
            do j=iys,iye
              do i=1,nnx
                t(i,j,iscl,km1) = tbc(i,j,iscl,2)
              enddo
            enddo
          enddo

        endif ! end of if block for root row

        call mpi_bcast(wtsfc,2,mpi_real8,0,mpi_comm_world,ierr)

        return
      end
      subroutine lower_curvy_dns(it,istage)

c --------- setup lower boundary condition for DNS
c           with curvy lower boundary.
c           This is no-slip for a moving wavy surface

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        include 'mpif.h'

c ------- set the number of variables to pass and
c         get back

        parameter (mvar= 9+3*nscl, jvar = 8+2*nscl)
        real u_level1(nnx,iys:iye,mvar), buf(jvar)
        real sbuf(jvar,mxs:mxe,iys:iye)
        real rbuf(jvar*nnx*(iye+1-iys))
        integer istatus(mpi_status_size)

        t_bndy = t_stage - t_zero

        if(i3d_surf == 0) then
          call l_bndy(t_bndy)
        elseif(i3d_surf == 1) then
          call get_surf(t_bndy)
        elseif(i3d_surf == 2) then
          call get_surf_gabls(t_bndy)
        elseif(i3d_surf == 3) then
          call get_surf_tank(t_bndy)
        elseif(i3d_surf == 4) then
          call get_surf_real(t_bndy)
        endif

c ----- broadcast level 1 data to processors above root

        k   = 1
        km1 = k - 1

        if(iss == 0) then

c ------------ only for root row = 0
c              get sums of surface conditions
c              and set surface boundary conditions

          do i=1,jvar
            buf(i) = 0.0
          enddo

          dz_i = dzu_i(k)

          do j=iys,iye
            do i=1,nnx
              ebc(i,j,2)  = amax1(e(i,j,k),sml_eg)
              pbc(i,j,2)  = 0.0
              pbc2(i,j,2) = 0.0

              ubc(i,j,2)  = 2.0*u_orbit(i,j) - u(i,j,k)
              vbc(i,j,2)  = 2.0*v_orbit(i,j) - v(i,j,k)

Cnew ------ bc for vertical velocity

              wbc(i,j,2)  = 2.0*gm(i,j,km1,5) - w(i,j,k)
c           wbc(i,j,2)  = -w(i,j,k)

              ufbc(i,j,2)  = -uf(i,j,k)
              vfbc(i,j,2)  = -vf(i,j,k)

c ------- no flow through the surface at t = t_stage

              wfbc(i,j,2)  = gm(i,j,km1,5)
            enddo
          enddo

          do iscl=1,nscl
            do j=iys,iye
              do i=1,nnx
                if (isfc(iscl)==1) then
                  tbc(i,j,iscl,2) = 2.0*Tbot(iscl) - t(i,j,iscl,k)
                endif
                if (isfc(iscl)==0) then
                  tbc(i,j,iscl,2) = t(i,j,iscl,k) +
     +                dzu(0)*wtsfc(1)/vis_s(i,j,k)
                endif
               !tbc(i,j,iscl,2) = 2.0*t_grnd(i,j,iscl) - t(i,j,iscl,k)
              enddo
            enddo
          enddo

c ------------ initialize u, v, w, t and derivatives at km1

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

              wf(i,j,km1) = wfbc(i,j,2)
              uf(i,j,km1) = ufbc(i,j,2)
              vf(i,j,km1) = vfbc(i,j,2)
            enddo
          enddo

          do iscl=1,nscl
            do j=iys,iye
              do i=1,nnx
                t(i,j,iscl,km1) = tbc(i,j,iscl,2)
              enddo
            enddo
          enddo

        endif ! end of if block for root row

        return
      end
      subroutine f_suft_curvy(rbuf,nnx,mxs,mxe,iys,iye,nscl,jvar,
     +      t11bc,t12bc,t13bc,t22bc,t23bc,t33bc,ustar,t3tbc,t_grnd,
     +      amol)

c ------ fill surface arrays on root processors

        real rbuf(jvar,mxs:mxe,iys:iye)
        real t11bc(nnx,iys:iye), t12bc(nnx,iys:iye),
     +      t13bc(nnx,iys:iye), t22bc(nnx,iys:iye),
     +      t23bc(nnx,iys:iye), t33bc(nnx,iys:iye),
     +      ustar(nnx,iys:iye), amol(nnx,iys:iye),
     +      t3tbc(nnx,iys:iye,nscl), t_grnd(nnx,iys:iye,nscl)

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

        return
      end
      subroutine drglaw_l(u_level1,mvar,it)

c --------- apply drag law across a set of gridpoints
c           with curvy lower boundary. each processor
c           applies log-law at several (i,j) for k = 1

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

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats

        real u_level1(nnx,iys:iye,mvar)
        real aij(3,3), us(3), usn(3)
        real tlv1(nscl), qflux(nscl), tsurf(nscl)

c ------- get matrix of direction cosines for surface coordinate system
c         and compute surface fluxes for a set of gridpoints [mxs:mxe,iys:iye]

        do j=iys,iye
          do i=mxs,mxe

c -------- vector aligned with \xi (t1)

            e1   = 1.0
            e3   = -u_level1(i,j,4)*u_level1(i,j,5)
            emag = sqrt(e1*e1 + e3*e3)
            e1   = e1/emag
            e3   = e3/emag

c -------- vector aligned with \eta (t2)

            f2   = 1.0
            f3   = -u_level1(i,j,4)*u_level1(i,j,6)
            fmag = sqrt(f2*f2 + f3*f3)
            f2   = f2/fmag
            f3   = f3/fmag

c -------- surface normal n

            g1   = -e3*f2
            g2   = -e1*f3
            g3   = e1*f2
            gmag = sqrt(g1*g1 + g2*g2 + g3*g3)
            g1   = g1/gmag
            g2   = g2/gmag
            g3   = g3/gmag

c -------- surface wind (2 vectors, us_1 and us_2)

            us1a      = u_level1(i,j,1)*e1 + u_level1(i,j,3)*e3
            us2a      = u_level1(i,j,2)*f2 + u_level1(i,j,3)*f3
            us(1)     = us1a*e1
            us(2)     = us2a*f2
            us(3)     = us1a*e3 + us2a*f3
            wind(i,j) = sqrt(us(1)**2 + us(2)**2 + us(3)**2)
            usn(1)    = us(1)/wind(i,j)
            usn(2)    = us(2)/wind(i,j)
            usn(3)    = us(3)/wind(i,j)

c -------- third dimension r = n cross us

            rb1   = g2*usn(3) - g3*usn(2)
            rb2   = g3*usn(1) - g1*usn(3)
            rb3   = g1*usn(2) - g2*usn(1)
            rbmag = sqrt(rb1**2 + rb2**2 + rb3**2)
            rb1   = rb1/rbmag
            rb2   = rb2/rbmag
            rb3   = rb3/rbmag

c -------- matrix of direction cosines

            aij(1,1) = usn(1)
            aij(1,2) = usn(2)
            aij(1,3) = usn(3)
            aij(2,1) = rb1
            aij(2,2) = rb2
            aij(2,3) = rb3
            aij(3,1) = g1
            aij(3,2) = g2
            aij(3,3) = g3

c --------- set simple scalar arrays
c           fix here if specified surface value

            do l=1,nscl
              qflux(l) = u_level1(i,j,9+l)
              tlv1(l)  = u_level1(i,j,9+nscl+l)
              tsurf(l) = u_level1(i,j,9+2*nscl+l)
            enddo

            zody = alog(z1b(i,j)/zo + 1.0)
            call mo(ibuoy,isfc,ismlt,vk,batag,zody,z1b(i,j),
     +          wind(i,j), ustar(i,j), nscl,
     +          tlv1(1), tsurf(1), qflux(1),
     +          zol, cd, ch, amol(i,j),RHB)

c ------------ build the momentum fluxes

            utau2      = -ustar(i,j)**2
            t11bc(i,j) = utau2*2.0*aij(1,1)*aij(3,1)
            t12bc(i,j) = utau2*(aij(1,1)*aij(3,2) +
     +          aij(3,1)*aij(1,2))
            t13bc(i,j) = utau2*(aij(1,1)*aij(3,3) +
     +          aij(1,3)*aij(3,1))
            t22bc(i,j) = utau2*2.0*aij(1,2)*aij(3,2)
            t23bc(i,j) = utau2*(aij(1,2)*aij(3,3) +
     +          aij(1,3)*aij(3,2))
            t33bc(i,j) = utau2*2.0*aij(1,3)*aij(3,3)

            do iscl=1,nscl
              t1tbc(i,j,iscl)  = 0.0
              t2tbc(i,j,iscl)  = 0.0
              t3tbc(i,j,iscl)  = qflux(iscl)
              t_grnd(i,j,iscl) = tsurf(iscl)
            enddo
          enddo
        enddo

        return
      end
      subroutine mo(ibuoy,isfc,ismlt,vk,batag,zody,z1,
     +      wind,utau,nscl,t1xy,tsfcc,qstar,
     +      zeta,cd,ch,amonin,RHB)
        use particles

c ---------- routine to find surface fluxes and values
c            using MO theory. Iterate for z/L

c            Gridpoints [mxs:mxe,iys:iye].

c    [assume:  z1 is constant, zody is constant, minimum wind ]

c            isfc = 1, then compute tstar and qstar given tsfcc
c                 = 0,   "     "    tsfcc given qstar

c    can add more scalars using definition below

        parameter (iter_mo = 30, zeta_min = -6.0, zeta_max = 3.0)
        parameter (wind_min = 0.00001)
        real t1xy(nscl), tsfcc(nscl), qstar(nscl), RHB
        integer isfc(nscl)

        if(ismlt == 1) then
          vk74in = 0.74/vk
        else
          vk74in = 1.0/vk
        endif
        batagk = batag*vk

        windm = amax1(wind_min,wind)

c ---------- limits for zeta

        zeta_mn = zeta_min
        zeta_mx = zeta_max
        if(isfc(1) == 0) then
          f_con = z1*batagk*qstar(1)/((windm*vk)**3)
        else
          d_theta = vk74in*(tsfcc(1) - t1xy(1))
          f_con   = z1*batagk*d_theta/((windm*vk)**2)
        endif

c --------- iteration for zeta

        do iter=1,iter_mo
          zeta_a = 0.5*(zeta_mn + zeta_mx)
          if(ismlt == 1) then
            call busngr(zeta_a,phim,phis,psim,psis)
          else
            call fzol(zeta_a,phim,phis,psim,psis)
          endif
          u_fac = (zody - psim)
          if(isfc(1) == 0) then
            f_new =  zeta_a + f_con*u_fac**3
          else
            t_fac = 1.0/amax1((zody - psis),0.01)
            f_new =  zeta_a + f_con*u_fac*u_fac*t_fac
          endif
          if(f_new < 0.0) then
            zeta_mn = zeta_a
          else
            zeta_mx = zeta_a
          endif
        enddo

c --------- check if neutral surface layer

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
        if(dnom_ch /= 0.) then
          ch = qstar(1)/dnom_ch
        endif

      !Do the lower BC for RH:

        if (isfc(2) == 0) then  !Set surface flux, solve for tsfcc

          thstar = -qstar(2)/utau
          tsfcc(2) = t1xy(2) - thstar*dnom
        elseif (isfc(2) == 1) then !Set tsfcc, solve for qstar

          tsfcc(2) = RHB/100.0*Mw/Ru/tsfcc(1)*mod_magnus(tsfcc(1))/rhoa
          thstar = (t1xy(2) - tsfcc(2))/dnom
          qstar(2) = -utau*thstar
        endif

c ------ example of other scalars (specified surface value and source)

c     qstar(2) = qstar(1)*(t1xy(2) - tsfcc(2))/(t1xy(1) - tsfcc(1))

c     tsfcc(2) = t1xy(2) - qstar(2)*(t1xy(1) - tsfcc(1))/qstar(1)

        return
      end
      subroutine upper

c ---- set boundary condition on upper boundary iz=nnz
c      option for special radiation boundary condition
c                 index f(.,.,1)  indicates upper.

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        include 'mpif.h'
        integer istatus(mpi_status_size)

        iz   = nnz
        izm1 = iz - 1
        izm2 = iz - 2
        izp1 = iz + 1
        izp2 = iz + 2

c     call dealias_grid_speed

        if(ibcu == -1) then ! channel flow type simulation

c --------- boundary conditions are gradient conditions

c            dudzbc = 0.0
c            dvdzbc = 0.0
c            dtdzbc = 0.0
c            wbc    = 0.0
c            dedzbc = 0.0

          do iy=iys,iye
            do ix=1,nnx

Cnew ------ bc for vertical velocity

              wbc(ix,iy,1) = 2.0*gm(ix,iy,iz,5) - w(ix,iy,iz)
              ebc(ix,iy,1) = e(ix,iy,iz)
              ubc(ix,iy,1) = u(ix,iy,iz)
              vbc(ix,iy,1) = v(ix,iy,iz)
              pbc(ix,iy,1) = 0.0
              pbc2(ix,iy,1)= 0.0
              ufbc(ix,iy,1) = uf(ix,iy,iz)
              vfbc(ix,iy,1) = vf(ix,iy,iz)

c --------- upper boundary = grid speed

              wfbc(ix,iy,1) = gm(ix,iy,iz,5)
            enddo
          enddo
          do iscl=1,nscl

c ---------- average scalar gradient

            dtdzf(iscl) = 0.0
          enddo

          call mpi_sum_xy(dtdzf,myid,iss,ise,nscl)

          do iscl=1,nscl
            do iy=iys,iye
              do ix=1,nnx
                tbc(ix,iy,iscl,1) = t(ix,iy,iscl,iz) +
     +              dtdzf(iscl)*dzu(nnzp1)
              enddo
            enddo
          enddo

        elseif(ibcu == 0) then

c --------- boundary conditions are gradient conditions

          do iy=iys,iye
            do ix=1,nnx
              wbc(ix,iy,1) = 2.0*gm(ix,iy,iz,5) - w(ix,iy,iz)
              ebc(ix,iy,1) = 0.0
              ubc(ix,iy,1) = u(ix,iy,iz)
              vbc(ix,iy,1) = v(ix,iy,iz)
              pbc(ix,iy,1) = 0.0
              pbc2(ix,iy,1)= 0.0
              ufbc(ix,iy,1) = uf(ix,iy,iz)
              vfbc(ix,iy,1) = vf(ix,iy,iz)

c --------- upper boundary = grid speed

              wfbc(ix,iy,1) = gm(ix,iy,iz,5)
            enddo
          enddo
          do iscl=1,nscl

c ---------- first get average scalar gradient

            dtdzf(iscl) = 0.0
            do iy=iys,iye
              do ix=1,nnx
                dtdzf(iscl) = dtdzf(iscl) + (t(ix,iy,iscl,nnz) -
     +              t(ix,iy,iscl,nnz-1))*dzu_i(nnz)
              enddo
            enddo
            dtdzf(iscl) = dtdzf(iscl)*fnxy
          enddo

          call mpi_sum_xy(dtdzf,myid,iss,ise,nscl)

          do iscl=1,nscl
            do iy=iys,iye
              do ix=1,nnx
                tbc(ix,iy,iscl,1) = t(ix,iy,iscl,iz) +
     +              dtdzf(iscl)*dzu(nnzp1)
              enddo
            enddo
          enddo
        elseif(ibcu == 2) then

c --------- boundary conditions are Sullivan et al. (2000) Couette flow

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

c --------- upper boundary = grid speed

              wfbc(ix,iy,1) = gm(ix,iy,iz,5)
            enddo
          enddo
          do iscl=1,nscl

c ---------- first get average scalar gradient

            dtdzf(iscl) = 0.0
            do iy=iys,iye
              do ix=1,nnx
                dtdzf(iscl) = dtdzf(iscl) + (t(ix,iy,iscl,nnz) -
     +              t(ix,iy,iscl,nnz-1))*dzu_i(nnz)
              enddo
            enddo
            dtdzf(iscl) = dtdzf(iscl)*fnxy
          enddo

          call mpi_sum_xy(dtdzf,myid,iss,ise,nscl)

          do iscl=1,nscl
            do iy=iys,iye
              do ix=1,nnx
                tbc(ix,iy,iscl,1) = t(ix,iy,iscl,iz) +
     +              dtdzf(iscl)*dzu(nnzp1)
              enddo
            enddo
          enddo
        elseif(ibcu == 1) then

c ------------- special if iradup boundary condition
c               get estimate of w from continuity and
c               linearized relation for pressure

          xmeanp = 0.0
          grad_ug = ug(nnz) - ug((nnz-1))
          do iy=iys,iye
            do ix=1,nnx
              wfbc(ix,iy,1) = wf(ix,iy,izm1)-
     +            (ux(ix,iy,iz) + vy(ix,iy,iz))*dzw(iz)
              pbc(ix,iy,1) = .5*(wf(ix,iy,izm1) + wfbc(ix,iy,1))
              ebc(ix,iy,1) = 0.0
              ubc(ix,iy,1) = u(ix,iy,iz) + grad_ug
              vbc(ix,iy,1) = v(ix,iy,iz)
              pbc2(ix,iy,1)=0.5*(u(ix,iy,iz)**2 + v(ix,iy,iz)**2) +
     +            0.25*(wf(ix,iy,izm1)**2 + wfbc(ix,iy,1)**2)
              xmeanp = xmeanp + pbc2(ix,iy,1)

c -------- FIX upper radiation boundary conditions

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
     +              dtdzf(iscl)*dzu(nnzp1)
              enddo
            enddo
          enddo
          xmeanp = xmeanp*fnxy
          do iy=iys,iye
            do ix=1,nnx
              pbc2(ix,iy,1) = pbc2(ix,iy,1) - xmeanp
            enddo
          enddo

c ---------- end if block

        endif

        do iy=iys,iye
          do ix=1,nnx
            w(ix,iy,izp1)   = wbc(ix,iy,1)
            e(ix,iy,izp1)   = ebc(ix,iy,1)

            wf(ix,iy,iz)   = wfbc(ix,iy,1)
            uf(ix,iy,izp1) = ufbc(ix,iy,1)
            vf(ix,iy,izp1) = vfbc(ix,iy,1)

            u(ix,iy,izp1) = ubc(ix,iy,1)
            v(ix,iy,izp1) = vbc(ix,iy,1)
            r3(ix,iy,izp1)= 0.0
            r5(ix,iy,izp1)= 0.0

c ---------- set derivatives at top of box (wx,wy not needed)
c            ux,uy,vx,vy are used in e production, but neglect
c            at top of box becuase of bc

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

        return
      end
      subroutine upper_dns

c ---- set boundary condition on upper boundary iz=nnz
c ---- for DNS conditions (no-slip, molecular diffusion)
c ---- assumes the grid transformation has completely decayed

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        include 'mpif.h'
        integer istatus(mpi_status_size)

        iz   = nnz
        izm1 = iz - 1
        izm2 = iz - 2
        izp1 = iz + 1
        izp2 = iz + 2

c --------- boundary conditions are Sullivan et al. (2000) Couette flow

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

c --------- upper boundary = grid speed

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
              endif
              if (isfc(iscl)==0) then
                tbc(ix,iy,iscl,1) = t(ix,iy,iscl,iz) +
     +              dzu(izp1)*wtsfc(iscl)/vis_s(ix,iy,izp1)
              endif
            enddo
          enddo
        enddo

        do iy=iys,iye
          do ix=1,nnx
            w(ix,iy,izp1)   = wbc(ix,iy,1)
            e(ix,iy,izp1)   = ebc(ix,iy,1)

            wf(ix,iy,iz)   = wfbc(ix,iy,1)
            uf(ix,iy,izp1) = ufbc(ix,iy,1)
            vf(ix,iy,izp1) = vfbc(ix,iy,1)

            u(ix,iy,izp1) = ubc(ix,iy,1)
            v(ix,iy,izp1) = vbc(ix,iy,1)
c ------------- note w and e nnz+1 values are not needed
c        w(ix,iy,izp1) = wbc(ix,iy,1)
c        e(ix,iy,izp1) = ebc(ix,iy,1)
            r3(ix,iy,izp1)= 0.0
            r5(ix,iy,izp1)= 0.0

c ---------- set derivatives at top of box (wx,wy not needed)
c            ux,uy,vx,vy are used in e production, but neglect
c            at top of box becuase of bc

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

        return
      end
      subroutine psolver(istage)

c -------- solver poisson equation for the pressure

c          iflat == 0 ; no iteration flat bottom
c          iflat /= 0 ; iterate using diagonal preconditioner
c                       set max number of iterations ip_max
c                       and convergence tolerance p_tol
c          solve for fluctuating pressure then mean value.
c          array p(.,.,.) is p' array pm(.) is the mean value

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats

        call get_avg_j

        do k=izs-1,ize+1
          do j=iys,iye
            do i=1,nnx
              p(i,j,k) = 0.0
            enddo
          enddo
        enddo

c -------- special sequence for flat bottom

        if(iflat == 0) then
          ip_iter = 0
          call build_rhs(ip_iter)
          call div_flux_vel
          call pressure
          call build_rhs(ip_iter)
        else

c ------- iteration sequence

          ip_iter = 0
    1     continue
          ip_iter = ip_iter + 1
          call build_rhs(ip_iter)
          call div_flux_vel
          call precon_p(ip_iter)
          call pressure
          call update_p(istage,ip_iter,ip_con)
          if(ip_iter < ip_max .and. ip_con == 0) go to 1

c -------- use last pressure iteration in flux velocity

          call build_rhs(ip_iter)
        endif

        return
      end
      subroutine build_rhs(ip_iter)

c ---------- build the rhs of the poisson equation
c            and update the flux velocities
c            use new time metrics. use wf equation to
c            get dp/dz at k = 0

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

c ------------ Send p to processor above the current myid.
c              to get p(.,.,izs-1)

        nb = myid - ncpu_s
        nt = myid + ncpu_s
        if(iss == 0) then
          nb = mpi_proc_null
        endif
        if(ise == numprocs-1) then
          nt = mpi_proc_null
        endif
        nsend = nnx*(iye + 1 - iys)
        nrecv = nsend
        do j=iys,iye
          do i=1,nnx
            fs(i,j) = p(i,j,ize)
          enddo
        enddo

        call mpi_sendrecv(
     +      fs(1,iys),nsend,mpi_real8,nt,2,
     +      fr(1,iys),nrecv,mpi_real8,nb,2,
     +      mpi_comm_world,istatus,ierr)
        if(iss /= 0) then
          do j=iys,iye
            do i=1,nnx
              p(i,j,izs-1) = fr(i,j)
            enddo
          enddo
        endif

        if(iss == 0) then

c -------- flux velocity is set = grid speed at t = t_stage2
c          updated boundary speed at the next time step
c          picked to match future grid speed

          do j=iys,iye
            do i=1,nnx
              wf(i,j,izs-1) = wfbc(i,j,3)
            enddo
          enddo
        endif

c --------- get bottom ghost point for fluctuating pressure

        if(izs == 1) then

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
     +          nnx,nny,ixs,ixe,ix_s,ix_e,
     +          iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
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

          elseif(lpbc == 1) then

            w1 = -2.0*dzw(izs)/(dzw(izs) + dzw(izs+1))
            w2 = 1.0 - w1
            do j=iys,iye
              do i=1,nnx
                p(i,j,izs-1) = p(i,j,izs)*w2 + w1*p(i,j,izs+1)
              enddo
            enddo

          else

            do j=iys,iye
              do i=1,nnx
                p(i,j,izs-1) = p(i,j,izs)
              enddo
            enddo
          endif ! end of check for lpbc

        endif ! end of check for izs

        if(ize == nnz) then
          do j=iys,iye
            do i=1,nnx
              p(i,j,ize+1) = p(i,j,ize)
            enddo
          enddo
        endif

c --------- save the last pressure iteration

        do k=izs-1,ize+1
          do j=iys,iye
            do i=1,nnx
              p_old(i,j,k) = p(i,j,k)
            enddo
          enddo
        enddo

        do k=izs,ize+1
          do j=iys,iye
            do i=1,nnx
              fnt1(i,j,k) = p_old(i,j,k)
            enddo
          enddo
          call xderivp(fnt1(1,iys,k),trigx(1,1),xk(1),nnx,iys,iye)
        enddo

        do k=izs,ize
          kp1   = k + 1
          km1   = k - 1
          do j=iys,iye
            do i=1,nnx
              pu_zeta   = weit(k)*p_old(i,j,kp1) + weitm(k)*p_old(i,j,k)
              pl_zeta   = weit(km1)*p_old(i,j,k) +
     +            weitm(km1)*p_old(i,j,km1)
              dpdz_c    = (pu_zeta - pl_zeta)*dzw_i(k)
              uf(i,j,k) = ufi(i,j,k) - dtgama*dpdz_c*
     +            gm2(i,j,k,2)*gm2(i,j,k,1)
              vf(i,j,k) = vfi(i,j,k) - dtgama*dpdz_c*
     +            gm2(i,j,k,3)*gm2(i,j,k,1)
              dpdz      = (p_old(i,j,kp1) - p_old(i,j,k))*dzu_i(kp1)

              facu      = gm2(i,j,kp1,1)*(gm2(i,j,kp1,2)**2 +
     +            gm2(i,j,kp1,3)**2 + gm2(i,j,kp1,4)**2)
              facl      = gm2(i,j,k,1)*(gm2(i,j,k,2)**2 +
     +            gm2(i,j,k,3)**2 + gm2(i,j,k,4)**2)
              wf(i,j,k) = wfi(i,j,k) - dtgama*(weit(k)*facu +
     +            weitm(k)*facl)*dpdz
            enddo
          enddo
          do j=iys,iye
            do i=1,nnx
              uf(i,j,k) = uf(i,j,k) - dtgama*fnt1(i,j,k)*gm2(i,j,k,1)
              facu      = gm2(i,j,kp1,1)*gm2(i,j,kp1,2)*fnt1(i,j,kp1)
              facl      = gm2(i,j,k,1)*gm2(i,j,k,2)*fnt1(i,j,k)
              wf(i,j,k) = wf(i,j,k) - dtgama*
     +            (weit(k)*facu + weitm(k)*facl)
            enddo
          enddo
        enddo

c -------- get the y contributions

        do k=izs,ize+1
          do j=iys,iye
            do i=1,nnx
              fnt1(i,j,k) = p_old(i,j,k)
            enddo
          enddo
        enddo
        call yd_mpi(fnt1(1,iys,izs),trigx(1,2),yk(1),
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs,ize+1,myid,ncpu_s,numprocs)

        do k=izs,ize
          kp1 = k + 1
          km1 = k - 1
          do j=iys,iye
            do i=1,nnx
              vf(i,j,k) = vf(i,j,k) - dtgama*fnt1(i,j,k)*gm2(i,j,k,1)
              facu      = gm2(i,j,kp1,1)*gm2(i,j,kp1,3)*fnt1(i,j,kp1)
              facl      = gm2(i,j,k,1)*gm2(i,j,k,3)*fnt1(i,j,k)
              wf(i,j,k) = wf(i,j,k) - dtgama*
     +            (weit(k)*facu + weitm(k)*facl)
            enddo
          enddo
        enddo

        if(ise == numprocs-1) then

c -------- upper boundary speed at t = t_stage2

          do j=iys,iye
            do i=1,nnx
              wf(i,j,ize) = wfbc(i,j,4)
            enddo
          enddo
        endif

c ------------ make sure continuity is satisfied

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

c ---------- gather partial sums

        call mpi_sum_z(wf_sum(1),i_root,myid,nnz,1)

        do k=izs,ize
          do j=iys,iye
            do i=1,nnx
              wf(i,j,k) = wf(i,j,k) - wf_sum(k)
            enddo
          enddo
        enddo

        return
      end
      subroutine div_flux_vel

c ----------- get the divergence of the flux velocity
c             and impose the boundary conditions on wf
c             new time metrics and weight for pressure

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        include 'mpif.h'
        integer istatus(mpi_status_size)

        real fnt1(nnx,iys:iye,izs:ize)
        real fs(nnx,iys:iye), fr(nnx,iys:iye)
        real sum_div(1:nnz)

        gami = 1.0/dtgama

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
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)

c ------------ Send updated wf (from psolver)
c              to processor above the current myid.

        nb = myid - ncpu_s
        nt = myid + ncpu_s
        if(iss == 0) then
          nb = mpi_proc_null
        endif
        if(ise == numprocs-1) then
          nt = mpi_proc_null
        endif
        nsend = nnx*(iye + 1 - iys)
        nrecv = nsend
        do j=iys,iye
          do i=1,nnx
            fs(i,j) = wf(i,j,ize)
          enddo
        enddo

        call mpi_sendrecv(
     +      fs(1,iys),nsend,mpi_real8,nt,2,
     +      fr(1,iys),nrecv,mpi_real8,nb,2,
     +      mpi_comm_world,istatus,ierr)
        if(iss /= 0) then
          do j=iys,iye
            do i=1,nnx
              wf(i,j,izs-1) = fr(i,j)
            enddo
          enddo
        endif

        do k=izs,ize
          do j=iys,iye
            do i=1,nnx
              div_wf   = (wf(i,j,k) - wf(i,j,k-1))*dzw_i(k)
              p(i,j,k) = (p(i,j,k) + fnt1(i,j,k) + div_wf)*gami
            enddo
          enddo
        enddo

c debug -------

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

c ---------- gather partial sums for rj_dif

        call mpi_sum_z(sum_div,i_root,myid,nnz,1)

        return
      end
      subroutine get_avg_j

c ------- compute the average of (1/J,J) for use in preconditioner
c         metrics are at the new time

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats

c --------- 1/J values stored at cell centers

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

c ---------- gather partial sums for avg_j

        call mpi_sum_z(avg_j(1),i_root,myid,nnz,1)

c ---------- J values stored at cell faces

        avg_j(0)     = avg_j(1)
        avg_j(nnz+1) = avg_j(nnz)
        do k=1,nnz
          fac        = avg_j(k+1)*weit(k) + avg_j(k)*weitm(k)
          avg_j_i(k) = 1.0/fac
        enddo
        avg_j_i(0)     = 1.0/avg_j(1)
        avg_j_i(nnz+1) = 1.0/avg_j(nnz)

        return
      end
      subroutine precon_p(ip_iter)

c -------- evaluate the preconditioner at the last iteration
c          for p, p_old has last iteration
c          new weight for pressure

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats

        real fnt1(nnx,iys:iye,izs:ize)

        gami   = 1.0/dtgama
        relax  = accel(ip_iter)
        relax1 = 1.0 - relax
        relax1 = 1.0

c ------- get x contribution

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

c ------- get y contribution

        do k=izs,ize
          do j=iys,iye
            do i=1,nnx
              fnt1(i,j,k) = p_old(i,j,k)
            enddo
          enddo
        enddo
        call ydd_mpi(fnt1(1,iys,izs),trigx(1,2),ykk(1),
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
        do k=izs,ize
          do j=iys,iye
            do i=1,nnx
              p(i,j,k) = p(i,j,k) + fnt1(i,j,k)*avg_j(k)*relax1
            enddo
          enddo
        enddo

c ------- get z contribution
c         always use dp/dz = 0 for pre-conditioner
c         at izs=1

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
              dpdz_u=facu_i*(p_old(i,j,kp1)-p_old(i,j,k))*dzu_i(kp1)
              p(i,j,k) = p(i,j,k) + dpdz_u*dzw_i(k)*relax1
            enddo
          enddo
        endif

        return
      end
      subroutine xdd(ax,trigx,xkk,nnx,iys,iye)

c -------- get multiple x second derivatives using fftpack routines
c          use fftpack storage a0, (a1,b1), (a2,b2),...,an
c          assumes that first wavenumber xk(1) = 0.0

c          assumes that wavenumbers are normalized by number of points

        real xkk(nnx), trigx(2*nnx+15), ax(nnx,iys:iye)

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

        return
      end
      subroutine ydd_mpi(ay,trigy,ykk,
     +      nx,ny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)

c -------- get multiple y second derivatives using fftpack routines and mpi
c          use fftpack storage a_0, (a1,b1), (a2,b2), ...,
c          assumes that first wavenumber yk(1) = 0.0
c          wavenumbers are normalized by number of points, ny

        real ykk(ny), trigy(2*ny+15), ay(nx,iys:iye,iz1:iz2)
        real ayt(ny,ixs:ixe,iz1:iz2)

        integer ix_s(0:np-1), ix_e(0:np-1),
     +      iy_s(0:np-1), iy_e(0:np-1)

        call xtoy_trans(ay,ayt,nx,ny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)

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
     +      iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)

        return
      end
      subroutine update_p(istage,ip_iter,ip_con)

c -------- check convergence and update p guess

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        include 'mpif.h'
        integer istatus(mpi_status_size)
        data small_p /1.0e-10/
        save small_p

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

        call mpi_allreduce(ptest_loc,ptest,1,mpi_real8,
     +      mpi_max,mpi_comm_world,ierror)

        ip_con = 0
        if(ptest < ptol) ip_con = 1
        if(istage == 1) then
          p_iter = float(ip_iter)
          p_err  = ptest
        endif

        return
      end
      subroutine comp2

c ------- add p gradients to rhs of (u,v,w). Use already defined p
c         at izs-1 and ize+1 to get (u,v,w). see build_rhs

c         multiplied by gm2 to get (u,v,w,e,theta) at the new time
c         pressure gradient is at new time and is weighted by dtgama

c         included mean pressure in vertical gradients in sr. rhs_uvw

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats

        real fnt1(nnx,iys:iye,izs:ize)

c ---------  ghost point pressure comes from build_rhs

        do k=izs,ize
          km1 = k - 1
          kp1 = k + 1
          do j=iys,iye
            do i=1,nnx
              fnt1(i,j,k) = p(i,j,k)*gm2(i,j,k,1)
            enddo
          enddo
          call xderivp(fnt1(1,iys,k),trigx(1,1),xk(1),nnx,iys,iye)

          do j=iys,iye
            do i=1,nnx
              p_kp1     = p(i,j,kp1)
              p_k       = p(i,j,k)
              p_km1     = p(i,j,km1)

c ------- u equation

              pu_zeta   = weit(k)*p_kp1*gm2(i,j,kp1,2)*
     +            gm2(i,j,kp1,1) +
     +            weitm(k)*p_k*gm2(i,j,k,2)*gm2(i,j,k,1)
              pl_zeta   = weit(km1)*p_k*gm2(i,j,k,2)*gm2(i,j,k,1) +
     +            weitm(km1)*p_km1*gm2(i,j,km1,2)*
     +            gm2(i,j,km1,1)

              dpdz_c    = (pu_zeta - pl_zeta)*dzw_i(k)

c --------  add pressure gradient to rhs

              r1(i,j,k) = r1(i,j,k) - (fnt1(i,j,k) + dpdz_c)
              u(i,j,k)  = (u(i,j,k) - dtgama*(fnt1(i,j,k) + dpdz_c))*
     +            gm2(i,j,k,4)

c ------- w equation (use zeta_z/J = 1)

              pu_zeta   = weit(k)*p_kp1 + weitm(k)*p_k
              pl_zeta   = weit(km1)*p_k + weitm(km1)*p_km1
              dpdz_c    = (pu_zeta - pl_zeta)*dzw_i(k)

c -------- add pressure gradient to rhs

              r3(i,j,k) = r3(i,j,k) - dpdz_c
              w(i,j,k)  = (w(i,j,k) - dtgama*dpdz_c)*gm2(i,j,k,4)

c ------- e equation correct for Jacobian, result is e (not e/J)

              e(i,j,k) = e(i,j,k)*gm2(i,j,k,4)
            enddo
          enddo
        enddo

c ------ v equation

        do k=izs,ize
          do j=iys,iye
            do i=1,nnx
              fnt1(i,j,k) = p(i,j,k)*gm2(i,j,k,1)
            enddo
          enddo
        enddo
        call yd_mpi(fnt1(1,iys,izs),trigx(1,2),yk(1),
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
        do k=izs,ize
          kp1   = k + 1
          km1   = k - 1
          do j=iys,iye
            do i=1,nnx
              p_kp1     = p(i,j,kp1)
              p_k       = p(i,j,k)
              p_km1     = p(i,j,km1)
              pu_zeta   = weit(k)*p_kp1*gm2(i,j,kp1,3)*
     +            gm2(i,j,kp1,1) +
     +            weitm(k)*p_k*gm2(i,j,k,3)*
     +            gm2(i,j,k,1)
              pl_zeta   = weit(km1)*p_k*gm2(i,j,k,3)*
     +            gm2(i,j,k,1) +
     +            weitm(km1)*p_km1*gm2(i,j,km1,3)*
     +            gm2(i,j,km1,1)
              dpdz_c    = (pu_zeta - pl_zeta)*dzw_i(k)

c --------  add pressure gradient to rhs now

              r2(i,j,k) = r2(i,j,k) - (fnt1(i,j,k) + dpdz_c)
              v(i,j,k)  = (v(i,j,k) - dtgama*(fnt1(i,j,k) + dpdz_c))*
     +            gm2(i,j,k,4)
            enddo
          enddo
        enddo

c ------------ scalar equations correct for J, result is theta (not theta/J)

        do k=izs,ize
          do l=1,nscl
            do j=iys,iye
              do i=1,nnx
                t(i,j,l,k) = t(i,j,l,k)*gm2(i,j,k,4)
              enddo
            enddo
          enddo
        enddo

        return
      end
      subroutine pressure

c -------- solve for fluctuating pressure using a matrix transpose
c          across mpi tasks and tridiagonal solver.
c          The transposed array
c          is dimensioned (0:nnz+1). Values
c          (0 & nnz+1) are not needed but are useful in the
c          matrix transpose when we return (see send_ztox).
c          On exit p is defined at all [izs-1:ize+1].

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

c ------------ Fourier analyze the right hand side
c              at all iz = izs,ize. results are in pfft

        call fft2d_mpi(p(1,iys,izs),pfft(1,jxs,izs),trigx(1,1),trigc,
     +      nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +      izs,ize,myid,ncpu_s,numprocs,-2)

c ------------ Fourier analyze the radiation bc arrays

        if(ibcu == 1) then
          call fft2d_mpi(ptop(1,iys,1),ptopfft(1,jxs,1),
     +        trigx(1,1),trigc,
     +        nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +        1,2,myid,ncpu_s,numprocs,-2)
        endif

c ---------- transpose first and last index of array
c            the order of pfft is (y,x,z)

        call xtoz_trans(pfft,pt,nny,nnz,jys,jye,jy_s,jy_e,
     +      jxs,jxe,izs,ize,iz_s,iz_e,myid,ncpu_s,
     +      numprocs)
        call solve_trid(pt, ptopfft)

c ------------- transpose back

        call ztox_trans(pt,pfft,nny,nnz,jys,jye,jy_s,jy_e,
     +      jxs,jxe,izs,ize,iz_s,iz_e,myid,ncpu_s,
     +      numprocs)

        iz_ee = ize+1
        if(ise == numprocs-1) then
          iz_ee = ize
        endif

c --------- inverse fft at all iz=izs,iz_ee to get p
c           see z indices

        call fft2d_mpi(p(1,iys,izs),pfft(1,jxs,izs),trigx(1,1),trigc,
     +      nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +      izs,iz_ee,myid,ncpu_s,numprocs,2)

c -------- partial sums for pressure

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

        do iz=izs,iz_ee
c        psum(iz) = -psum(iz) + engz(iz) + c23*engsbz(iz)
          do iy=iys,iye
            do ix=1,nnx
              p(ix,iy,iz) = p(ix,iy,iz) - psum(iz)
            enddo
          enddo
        enddo

        return
      end
      subroutine solve_trid(pt, ptop)

c --------- tridiagonal solver. odd order for ptop, ptop2
c           because of 2d-fft

        use pars
        use con_data
        use con_stats

        real ptop(nny,jxs:jxe,1:2)
        real pt(0:nnz+1,jxs:jxe,jys:jye)
        real aa(nnz,jxs:jxe),bb(nnz,jxs:jxe),
     +      dd(nnz,jxs:jxe),rh(nnz,jxs:jxe)
        real fac_u(nnz), fac_l(nnz), fac_a(nnz)

        do iz=1,nnz
          fac_u(iz) = avg_j_i(iz)/(dzw(iz)*dzu(iz+1))
          fac_l(iz) = avg_j_i(iz-1)/(dzw(iz)*dzu(iz))
          fac_a(iz) = fac_l(iz) + fac_u(iz)
        enddo

        do kp=jys,jye
          do lp=jxs,jxe
            do iz=2,nnz-1
              bb(iz,lp)  = fac_l(iz)
              aa(iz,lp)  = fac_u(iz)
              dd(iz,lp)  = -xks(lp,kp)*avg_j(iz) - fac_a(iz)
              rh(iz,lp)  = pt(iz,lp,kp)
            enddo
          enddo

c --------------- lower boundary, fill exterior pressure (not used)

          do lp=jxs,jxe
            bb(1,lp)  = 1.0
            aa(1,lp)  = fac_u(1)
            dd(1,lp)  = -xks(lp,kp)*avg_j(1) - fac_u(1)
            rh(1,lp)  = pt(1,lp,kp)
            pt(0,lp,kp) = 0.0
          enddo

c --------------- upper boundary, fill exterior pressure (not used)

          if(ibcu == 1) then
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

c ---------------- special situation for zeroth mode
c                  makes mean pressure = 0

          if(kp == 1 .and. jxs == 1) then
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

c --------------- solve system

          call tridv(bb,dd,aa,rh,nnz,jxs,jxe)
          do lp=jxs,jxe
            do iz=1,nnz
              pt(iz,lp,kp) = rh(iz,lp)
            enddo
          enddo
        enddo

        return
      end
      subroutine solve_mean_p(wf_sum,sum_fz)

c --------- tridiagonal solver for mean pressure from wf equation
c           set ghost values of pm also

        use pars
        use con_data
        use con_stats

        real wf_sum(0:nnz), sum_fz(nnz)
        real aa(nnz),bb(nnz),dd(nnz),rh(nnz)

        gami = 1.0/dtgama
c     tiny = 1.0e-20
        tiny = 0.0

        do k=2,nnz-1
          fac_l  = sum_fz(k-1)*dzu_i(k)
          fac_u  = sum_fz(k)*dzu_i(k+1)
          bb(k)  = fac_l
          aa(k)  = fac_u
          dd(k)  = -(fac_u + fac_l + tiny)
          rh(k)  = gami*(wf_sum(k) - wf_sum(k-1))
        enddo

c --------------- lower boundary, assume average wf at boundary = 0

        k       = 1
        fac_u   = sum_fz(k)*dzu_i(k+1)
        bb(k)   = 0.0
        aa(k)   = fac_u
        dd(k)   = -(fac_u + tiny)
        rh(k)   = gami*(wf_sum(k) - wf_sum(k-1))

c --------------- upper boundary, average wf at boundary = 0

        k       = nnz
        fac_l   = sum_fz(k-1)*dzu_i(k)
        bb(k)   = fac_l
        aa(k)   = 0.0
        dd(k)   = -(fac_l + tiny)
        rh(k)   = gami*(wf_sum(k) - wf_sum(k-1))

c --------------- solve system

        call trids(bb,dd,aa,rh,nnz)
        do k=1,nnz
          pm(k) = rh(k)
        enddo

        if(lpbc == 1) then
          k       = 1
          w1      = -2.0*dzw(k)/(dzw(k) + dzw(k+1))
          w2      = 1.0 - w1
          pm(k-1) = pm(k)*w2 + w1*pm(k+1)
        elseif(lpbc == 0) then
          k       = 1
          pm(k-1) = pm(k)
        endif

        pm(nnz+1) = pm(nnz)

        return
      end
      subroutine solve_mean_p_old(wf_sum,sum_fz)

c --------- tridiagonal solver for mean pressure from wf equation
c           set ghost values of pm also

        use pars
        use con_data
        use con_stats

        real wf_sum(nnz), sum_fz(nnz)
        real aa(nnz),bb(nnz),dd(nnz),rh(nnz)

        gami = 1.0/dtgama
c     tiny = 1.0e-20
        tiny = 0.0

        do k=2,nnz-1
          fac_l  = sum_fz(k-1)/(dzw(k)*dzu(k))
          fac_u  = sum_fz(k)/(dzw(k)*dzu(k+1))
          bb(k)  = fac_l
          aa(k)  = fac_u
          dd(k)  = -(fac_u + fac_l + tiny)
          rh(k)  = gami*(wf_sum(k) - wf_sum(k-1))*dzw_i(k)
        enddo

c --------------- lower boundary, assume average wf at boundary = 0

        k       = 1
        fac_u   = sum_fz(k)/(dzw(k)*dzu(k+1))
        bb(k)   = 0.0
        aa(k)   = fac_u
        dd(k)   = -(fac_u + tiny)
        rh(k)   = gami*wf_sum(k)*dzw_i(k)

c --------------- upper boundary, average wf at boundary = 0

        k       = nnz
        fac_l   = sum_fz(k-1)/(dzw(k)*dzu(k))
        bb(k)   = fac_l
        aa(k)   = 0.0
        dd(k)   = -(fac_l + tiny)
        rh(k)   = -gami*wf_sum(k-1)*dzw_i(k)

c --------------- solve system

        call trids(bb,dd,aa,rh,nnz)
        do k=1,nnz
          pm(k) = rh(k)
        enddo

        if(lpbc == 1) then
          k       = 1
          w1      = -2.0*dzw(k)/(dzw(k) + dzw(k+1))
          w2      = 1.0 - w1
          pm(k-1) = pm(k)*w2 + w1*pm(k+1)
        elseif(lpbc == 0) then
          k       = 1
          pm(k-1) = pm(k)
        endif

        pm(nnz+1) = pm(nnz)

        return
      end
      subroutine trids(b,d,a,r,n)

c --- tridiagonal matrix solver with single vector

c --- input:   n   size of a,b,d and r
c              b   below diagonal elements (b(1) not used)
c              d   diagonal elements
c              a   above diagonal elements (a(n) not used)
c              r   right hand side

c --- output:  r   solution vector

        real b(n), d(n), a(n), r(n)

        if(n <= 1 ) then
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
  999   continue

        return
      end
      subroutine tridv(b,d,a,r,n,j1,j2)

c --- tridiagonal matrix solver with multiple vectors
c     (note j and i loops are reversed from cray version)

c --- input:   n   size of a,b,d and r
c              b   below diagonal elements (b(1) not used)
c              d   diagonal elements
c              a   above diagonal elements (a(n) not used)
c              r   right hand side
c              j1:j2  range of input vectors

c --- output:  r   solution vector

        real b(n,j1:j2), d(n,j1:j2), a(n,j1:j2), r(n,j1:j2)

        if(n <= 1 ) then
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
  999   continue

        return
      end
      subroutine get_derv

c ------- get ux,uy,vx,vy at all z for this node
c         using parallel fft. can be improved (?)
c         by using exchange to send derivatives

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats

        iz_ss = izs-1
        iz_ee = ize+1
        if(iss == 0) then
          iz_ss = izs
        endif
        if(ise == numprocs-1) then
          iz_ee = ize
        endif

c ------- make sure <w> = 0

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
     +        nnx,iys,iye)
          call xderivp(vx(1,iys,iz),trigx(1,1),xk(1),
     +        nnx,iys,iye)
          call xderivp(wx(1,iys,iz),trigx(1,1),xk(1),
     +        nnx,iys,iye)
        enddo

c ---------- get y derivatives for (u,v,w)

        call yd_mpi(uy(1,iys,izs-1),trigx(1,2),yk(1),
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
        call yd_mpi(vy(1,iys,izs-1),trigx(1,2),yk(1),
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
        call yd_mpi(wy(1,iys,izs-1),trigx(1,2),yk(1),
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)

        return
      end
      subroutine chk_div

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        real fnt1(nnx,iys:iye,izs:ize+1), fnt2(nnx,iys:iye)
        real div_c(nnz)

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
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
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
     +            (fnt2(i,j) + fnt1(i,j,k) +
     +            (wf(i,j,k) - wf(i,j,km1))*dzw_i(k))**2
            enddo
          enddo
          div_c(k) = div_c(k)*fnxy
        enddo
        call mpi_sum_z(div_c(1),i_root,myid,nnz,1)

        write(nprt,9010) (k, div_c(k),k=izs,ize)
 9010   format(' 9010 debug chk div',/,
     +      ' k ',5x,' div ',/,(i5,e20.10))

        return
      end
      subroutine get_means(istage)

c ------------ get means for all variables
c              for use in iso, surfvis, comp1, tke_stats.
c              simple average along coordinate lines

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        real fnt1(nnx,iys:iye,izs:ize+1), fnt2(nnx,iys:iye)

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
        if(ize == nnz) iz_ee = nnzp1
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

c -------- set e to minimum value

        do iz=izs-1,ize+1
          do iy=iys,iye
            do ix=1,nnx
              e(ix,iy,iz) = amax1(e(ix,iy,iz ),sml_eg)
            enddo
          enddo
        enddo

c ------- get terms which contribute to mean pressure
c         careful with the sum, get the mean p_star pressure
c         CHECK if rotational form of the equations correct pxym see flat code

        do iz=izs,ize
          izm1 = iz - 1
          do iy=iys,iye
            do ix=1,nnx
              q_temp     =  0.5*(u(ix,iy,iz)**2 +
     +            v(ix,iy,iz)**2 +
     +            w(ix,iy,iz)**2)
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

c ------------ save means and divergence for printout and compmn
c              all cpus have means over all z

        if(istage == 1) then
          do k=izs,ize
            do j=iys,iye
              do i=1,nnx
                fnt1(i,j,k) = vf(i,j,k)
              enddo
            enddo
          enddo
          call yd_mpi(fnt1(1,iys,izs),trigx(1,2),yk(1),
     +        nnx,nny,ixs,ixe,ix_s,ix_e,
     +        iys,iye,iy_s,iy_e,izs,ize,myid,ncpu_s,numprocs)
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
     +              (fnt2(i,j) + fnt1(i,j,k) +
     +              (wf(i,j,k) - wf(i,j,km1))*dzw_i(k))**2
              enddo
            enddo
            divz(k) = divz(k)*fnxy
          enddo
          call mpi_sum_z(divz(1),i_root,myid,nnz,1)

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

        return
      end
      subroutine spline(x,y,n,yp1,ypn,y2)
        integer n, nmax
        real yp1, ypn, x(n), y(n), y2(n)
        parameter (nmax=5000)
        integer i, k
        real p, qn, sig, un, u(nmax)
        if(yp1 > .99e30) then
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
     +        /(x(i) - x(i-1)))/(x(i+1) - x(i-1)) - sig*u(i-1))/p
        enddo
        if(ypn > .99e+30) then
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

        return
      end
      subroutine splint(xa,ya,y2a,n,x,y)
        integer n
        real x,y, xa(n), y2a(n), ya(n)
        integer k,khi,klo
        real a,b,h
        klo = 1
        khi = n
    1   continue
        if(khi - klo > 1) then
          k = (khi + klo)/2
          if(xa(k) > x) then
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
     +      ((a**3 - a)*y2a(klo) + (b**3 - b)*y2a(khi))*(h**2)/6.0

        return
      end
      subroutine hx_deriv(ax_s,trigx,xk,nnx,ax)

c -------- get an x derivative saving the input. Use fftpack routines
c          with fftpack storage a0, (a1,b1), (a2,b2),...,an
c          assumes that first wavenumber xk(1) = 0.0

c          assumes that wavenumbers are normalized by number of points

        real xk(nnx), trigx(2*nnx+15), ax_s(nnx), ax(nnx)

        fn = 1.0/float(nnx)
        do i=1,nnx
          ax(i) = ax_s(i)*fn
        enddo

        call rfftf(nnx,ax(1),trigx)

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

        return
      end
      subroutine gal_interp(f,trigx,xk,tb,fi)

c ---------- propagate array f forward in space by distance dx
c            input is fourier coefficients in fftpack order

        use pars
        use tank_pars

        real f(nnx), trigx(15+4*nnx), fi(nnx), xk(nnx), c_wav

        ncx = nnx/2 + 1

c --------- generate new coefficients

        cc    = 1.0
        ss    = 0.0
        fi(1) = f(1)*cc
        kk    = 1
        do i=2,nnx-1,2
          kk = kk + 1
          if(i3d_surf == 4) then
            c_wav  = sqrt(grav/xk(i))
            if(i>ncx) c_wav = sqrt(-grav/xk(i))
          else
            c_wav = c_speed
          endif
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

c --------- transform back to get function at correct x points

        call rfftb(nnx,fi(1),trigx)

        return
      end
      subroutine hdot_interp(f,trigx,xk,tb,fi)

c ----------- get hdot assuming packet moves with speed cs
        use pars
        use tank_pars

        real f(nnx), trigx(15+4*nnx), fi(nnx), xk(nnx), c_wav

        ncx = nnx/2 + 1

c --------- generate modified fourier coefficients

        fi(1) = 0.0
        kk    = 1
        do i=2,nnx-1,2
          kk = kk + 1
          if(i3d_surf == 4) then
            c_wav  = sqrt(grav/xk(i))
            if(i>ncx) c_wav = sqrt(-grav/xk(i))
          else
            c_wav = c_speed
          endif
          cc = cos(xk(kk)*-c_wav*tb)
          ss = sin(xk(kk)*-c_wav*tb)
          aa = f(i)
          bb = f(i+1)
          uk = xk(kk)*c_wav
          fi(i)     = uk*(aa*ss + bb*cc)
          fi(i+1)   = uk*(bb*ss - aa*cc)
        enddo
        fi(nnx) = 0.0

c --------- transform back to get function at correct x points

        call rfftb(nnx,fi(1),trigx)

        return
      end
      subroutine xderivp(ax,trigx,xk,nnx,iys,iye)

c -------- get multiple x derivatives using fftpack routines
c          use fftpack storage a0, (a1,b1), (a2,b2),...,an
c          assumes that first wavenumber xk(1) = 0.0

c          assumes that wavenumbers are normalized by number of points

        real xk(nnx), trigx(2*nnx+15), ax(nnx,iys:iye)

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

        return
      end
      subroutine fft2d_mpi(ax,at,trigx,trigc,nx,ny,
     +      jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +      iz1,iz2,myid,ncpu,np,isgn)

c -------- get 2d fft using fftpack routines and parallel mpi
c          use fftpack storage a0, (a1,b1), (a2,b2),...,

c         isgn = -1 do forward transform, get coefficients
c                   incoming array is ax(nx+2,iys:iye,iz1:iz2)
c                   outgoing array is ax(nx+2,iys:iye,iz1:iz2)

c         isgn = -2 do forward transform, get coefficients
c                   incoming array is ax(nx+2,iys:iye,iz1:iz2)
c                   outgoing array is at(ny,jxs:jxe,iz1:iz2)

c         isgn =  1 do inverse transform, move to physical space
c                   incoming array is ax(nx+2,iys:iye,iz1:iz2)
c                   outgoing array is ax(nx+2,iys:iye,iz1:iz2)

c         isgn =  2 do inverse transform, move to physical space
c                   incoming array is at(ny,jxs:jxe,iz1:iz2)
c                   outgoing array is ax(nx+2,iys:iye,iz1:iz2)

        real ax(nx+2,iys:iye,iz1:iz2), at(ny,jxs:jxe,iz1:iz2),
     +      trigx(2*nx+15), trigc(4*ny+15),
     +      a2d(2,ny), a_wrk(nx)
        integer jx_s(0:np-1), jx_e(0:np-1),
     +      iy_s(0:np-1), iy_e(0:np-1)

        nxp2 = nx + 2
        if(isgn < 0) then
          fn   = 1.0/(float(nx)*float(ny))

c ------ 1d fft in x over [iys,iye] for all z

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

c ------ 1d fft in y over [jxs,jxe] for all z

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

c ---- decide whether to transpose back or leave as is

          if(isgn == -1) then
            call ytox_trans(at,ax,nxp2,ny,jxs,jxe,jx_s,jx_e,
     +          iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
          endif

        else

c ---- decide whether to first transpose or leave as is

          if(isgn == 1) then
            call xtoy_trans(ax,at,nxp2,ny,jxs,jxe,jx_s,jx_e,
     +          iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)
          endif

c ------ 1d fft in y over [jxs,jxe] for all z

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

c ------  1d fft in x over [iys,iye] for all z

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

        return
      end
      subroutine fft2d_cmplx_mpi(ax,at,trigx,trigy,nx,nx2,ny,
     +      i2xs,i2xe,i2x_s,i2x_e,iys,iye,iy_s,iy_e,
     +      iz1,iz2,myid,ncpu,np)

c -------- get 2d inverse fft using complex fftpack routines
c          move from fourier space to physical space. specifically
c          designed for flow over waves with full complex wave spectrum

        real trigx(4*nx+15), trigy(4*ny+15), ax(nx2,iys:iye,iz1:iz2)
        real a2d(2,ny), at(ny,i2xs:i2xe,iz1:iz2)
        integer i2x_s(0:np-1), i2x_e(0:np-1), iy_s(0:np-1), iy_e(0:np-1)

c --------- transform in x direction

        do k=iz1,iz2
          do j=iys,iye
            call cfftb(nx,ax(1,j,k),trigx(1))
          enddo
        enddo

c --------- transform in y direction

        call xtoy_trans(ax,at,nx2,ny,i2xs,i2xe,i2x_s,i2x_e,
     +      iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)

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

        call ytox_trans(at,ax,nx2,ny,i2xs,i2xe,i2x_s,i2x_e,
     +      iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)

        return
      end
      subroutine yderiv(ay,trigy,yk,nnx,nny)

c -------- get multiple y derivatives using fftpack routines
c          use fftpack storage a_0, (a1,b1), (a2,b2), ...,
c          assumes that first wavenumber yk(1) = 0.0

c          assumes that wavenumbers are normalized by number of points

        real yk(nny), trigy(2*nny+15), ay(nnx,nny)
        real a_trans(nny)

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

        return
      end
      subroutine yd_mpi(ay,trigy,yk,
     +      nx,ny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)

c -------- get multiple y derivatives using fftpack routines and mpi
c          use fftpack storage a_0, (a1,b1), (a2,b2), ...,
c          assumes that first wavenumber yk(1) = 0.0
c          wavenumbers are normalized by number of points, ny

        real yk(ny), trigy(2*ny+15), ay(nx,iys:iye,iz1:iz2)
        real ayt(ny,ixs:ixe,iz1:iz2)

        integer ix_s(0:np-1), ix_e(0:np-1),
     +      iy_s(0:np-1), iy_e(0:np-1)

        call xtoy_trans(ay,ayt,nx,ny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)

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
     +      iys,iye,iy_s,iy_e,iz1,iz2,myid,ncpu,np)

        return
      end
      function rlim(d1,d2,d3)

c ------------- Cees's kappa=1/3 scheme
        r = (d1-d2+1.e-100)/(d2-d3-1.e-100)
        rlim = (d2-d3)*amax1(0.,amin1(r,amin1(1./6.+1./3.*r,1.)))
c ------------- Cees's kappa=-1 scheme
c     r = (d1-d2+1.e-100)/(d2-d3-1.e-100)
c     rlim = (d2-d3)*amin1(abs(r),0.5)
c ------------- first order upwind
c     rlim = 0.0
c ------------- QUICK scheme
c     rlim = -0.25*d2 - 0.125*d3 + 0.375*d1

        return
      end

      function ran1(idum)

c ----------- stolen from numerical recipes,p. 271

        integer idum, ia, im, iq, ir, ntab, ndiv
        real ran1, am, eps, rnmx
        parameter (ia=16807,im=2147483647,am=1.0/im,iq=127773,ir=2836.0,
     +      ntab=32,ndiv=1+(im-1)/ntab,eps=1.2e-07,rnmx=1.0-eps)
        integer j, k, iv(ntab), iy
        save iv, iy
        data iv /ntab*0/, iy /0/
        if(idum <= 0 .or. iy == 0) then
          idum = max(-idum,1)
          do j=ntab+8,1,-1
            k = idum/iq
            idum = ia*(idum - k*iq) - ir*k
            if(idum < 0) idum = idum + im
            if(j <= ntab) iv(j) = idum
          enddo
          iy = iv(1)
        endif
        k     = idum/iq
        idum  = ia*(idum - k*iq) - ir*k
        if(idum < 0) idum = idum + im
        j     = 1 + iy/ndiv
        iy    = iv(j)
        iv(j) = idum
        ran1  = min(am*iy, rnmx)

        return
      end

      function ranf()
        data inc /1/
        save inc, ix, ia, m, fm
        if(inc==1) then
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

c ----------- get stokes drift velocity for assumed wavelength stokesw
c             and wave amplitude stokesa. Changed sign for z.

        use pars
        use con_data
        use con_stats
        include 'mpif.h'

        if(iocean == 1) then

c ----------- compute stokes velocity for ocean pbls

c         stokesw = pi2/20.0
          stokesw = pi2/76.5
c         ak      = 0.04
          ak      = 0.00
c         stokesa = 1.0
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

        else

c ----------------- set stokes velocity = 0 for atmos. pbls

          do iz=1,nnzp1
            stokes(iz) = 0.0
          enddo
          stokess = 0.0
          udrift = 0.0
          vdrift = 0.0
        endif

        return
      end
      subroutine busngr(zeta,phim,phis,psim,psis)

c ---- Businger's version of similarity theory

        data pih /1.57079633/
        save pih

        if(zeta < 0.) then
          x=(1.0 - 15.0*zeta)**0.25
          phim = 1.0/x
          psim = 2.0*alog((1.0+x)/2.0) + alog((1.0+x*x)/2.0) -
     +        2.0*atan(x)+pih
          if(psim>2.0)psim=2.0
          y = sqrt(1.0-9.0*zeta)
          phis = 0.74/y
          psis = alog((1.0+y)/2.0)*2.0
        elseif(zeta > 0) then
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

        psimu(Y)  = 1.571 + 2.0*(alog(0.5*(1.0 + Y)) - atan(Y)) +
     +      alog(0.5 + 0.5*Y**2)
        psisu(Y)  = 2.0*alog(0.5 + 0.5*Y)
        psicu(Y,G)= (1.0 - G)*alog(abs(Y - 1.0))
     +      + 0.5*(G + 2.0)*alog(abs(Y**2 + Y + 1.0))
     +      - (2.0*G + 1.0) / sqrt(3.0) *
     +      atan((Y + 0.5)*2.0/sqrt(3.0))
        Xm(zol)   = (1.0 - 16.0*zol)**0.25
        Xs(zol)   = sqrt(1.0 - 16.0*zol)
        Xc(zol,f) =  abs(1.0 - f*zol)**(4.0/3.0)/(1.0 - f*zol)

        if(zeta>=0.0)       then
c                                          STABLE
          if(zeta<=1.0) then
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
          if(zeta>=zetam) then
            phim = 1.0 / Xm(zeta)
            psim = psimu(Xm(zeta))
          else
c                            use convective limit for momentum
            X = (1.0 - b3/a3 * zeta)**(1.0/3.0)

            fm = a3**(-1.0/3.0)
            phim = fm / Xc(zeta,b3/a3)
            psim = psimu(Xm(zetam))
     *      + psicu(Xc(zeta,b3/a3),fm)
     *      - psicu(Xc(zetam,b3/a3),fm)
          endif

c                                         UNSTABLE scalars
          if(zeta>=zetas) then
            phis = 1.0/Xs(zeta)
            psis = psisu(Xs(zeta))
          else
c                              use convective limit for scalars
            fs =   abs(a4)**(-1.0/3.0)*abs(a4)/a4
            phis = (a4 - b4*zeta)**(-1.0/3.0)
            psis = psisu(Xs(zetas))
     *      + psicu(Xc(zeta,b4/a4),fs)
     *      - psicu(Xc(zetas,b4/a4),fs)
          endif

        endif
        return
      end
      subroutine suft(it)

c ---------- iterate for zeta = z/L using bisection method
c            either businger or large functions can be specified

c            isfc = 0, specified surface heat flux
c                 = 1, specified surface temperature

        use pars
        use fields
        use con_data
        use con_stats
        real buf(3+nscl)
        include 'mpif.h'
        integer istatus(mpi_status_size)

        parameter (iter_mo = 30, zeta_min = -6.0, zeta_max = 3.0)

c ---------- limiting value for wind

        ufree = 0.07*(abs(batag*qstar(1)*dzw(1)))**(1./3.)

c ---- save old utau

        utausv = utau
        utau2  = utau*utau

        iz   = 1
        izp1 = iz + 1
        izm1 = iz - 1

        buf(1)  = 0.0
        buf(2)  = 0.0
        buf(3)  = 0.0
        tol     = 0.01
        do iy=iys,iye
          do ix=1,nnx
            buf(1) = buf(1) + u(ix,iy,iz)
            buf(2) = buf(2) + v(ix,iy,iz)
            wind(ix,iy) = sqrt((u(ix,iy,iz)+ugal)**2
     +          +v(ix,iy,iz)*v(ix,iy,iz))
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

c -------- get x-y slab sums

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

c ---------- limits for zeta

        zeta_mn = zeta_min
        zeta_mx = zeta_max
        if(isfc(1) == 0) then
          f_con = z1*batag*vk*qstar(1)/((windm*vk)**3)
        else
          d_theta = vk74in*(tsfcc(1) - t1xy(1))
          f_con   = z1*batag*vk*d_theta/((windm*vk)**2)
        endif

c --------- iteration for zeta

        do iter=1,iter_mo
          zeta_a = 0.5*(zeta_mn + zeta_mx)
          if(ismlt == 1) then
            call busngr(zeta_a,phim,phis,psim,psis)
          else
            call fzol(zeta_a,phim,phis,psim,psis)
          endif
          u_fac = (zody - psim)
          if(isfc(1) == 0) then
            f_new =  zeta_a + f_con*u_fac**3
          else
            t_fac = 1.0/(zody - psis)
            f_new =  zeta_a + f_con*u_fac*u_fac*t_fac
          endif
          if(f_new < 0.0) then
            zeta_mn = zeta_a
          else
            zeta_mx = zeta_a
          endif

        enddo

c --------- check if neutral surface layer

        if (ibuoy==0 .or. qstar(1) == 0.) then
          amonin    = 1000.
          zeta      = 0.0
          utau      = windm*vk/zody
          thstar(1) = 0.0
          t10xy(1)  = 0.0
          tsfcc(1)  = t1xy(1)
        else
          utau = windm*vk/(zody-psim)
          dnom = (zody-psis)*vk74in
          if(isfc(1) == 0) then
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

        if (utau>10.0) then
          write(6,9000)
          write(6,9200) utau,windm
          go to 9999
        endif
        if (t10xy(1)>0. .and. qstar(1) > 0.) then
          write(6,9000)
          write(6,9300) u1xy,v1xy,t1xy(1),
     +        tsfcc(1),amonin,utau,it
          go to 9999
        endif

      !Second scalar:
        if(isfc(2) == 0) then
          thstar(2) = -qstar(2)/utau
          tsfcc(2)  = t1xy(2)-thstar(2)*dnom
          t10xy(2)  = thstar(2)*dnom
        else
          thstar(2) = (t1xy(2) - tsfcc(2))/dnom
          t10xy(2)  = thstar(2)*dnom
          qstar(2)  = -utau*thstar(2)
        endif

c ---------- examples of two other scalars

c     c
c     c **** get flux of b scalar, specified surface value
c     c
c           dnom      = (zody-psis)*vk74in
c           thstar(2) = (t1xy(2)-tsfcc(2))/dnom
c           qstar(2)  = -thstar(2)*utau
c           t10xy(2)  = thstar(2)*dnom
c           aut3m(2)  =  qstar(2)

c **** get surface value of c scalar, specified surface flux

c     dnom      = (zody-psis)*vk74in
c     thstar(2) = -qstar(2)/utau
c     tsfcc(2)  = t1xy(2) - dnom*thstar(2)
c     t10xy(2)  = thstar(2)*dnom
c     aut3m(2)  = qstar(2)

        zol = zeta
        hol = zol*zi/z1

c ---- note roundoff problem in angles if close to multiples of pi

        tep = u1xy/vsfc
        if(tep>1.)  tep = 1.0
        if(tep<-1.) tep = -1.0
        thta      = acos(tep)
        utau2     = utau*utau
        au13m     = -utau2*cos(thta)
        au23m     = -utau2*sin(thta)*sign(1.,v1xy)
        aut3m(1)  =  qstar(1)

        return

c -------- iteration did not converge

 9999   continue
 9000   format(' Trouble in SR. suft')
 9200   format(' Stop because utau = ',e15.6,' windm = ',e15.6)
 9300   format(' ** CHECK SFC U = ',e15.6,' V=',e15.6,' T,TS = ',2e15.6,
     +      ' L =',e15.6,' U_* = ',e15.6,' AT IT = ',i5)
        call mpi_abort(mpi_comm_world,icode,ierr)
      end
      subroutine sufto(it)

        use pars
        use fields
        use con_data
        use con_stats
        real buf(3+nscl)

c ------- version of similarity theory adpated for ocean flows
c      option to use businger or large version of similarity theory

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
     +          +v(ix,iy,iz)*v(ix,iy,iz))
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

c -------- get x-y slab sums

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

        t10xy(1)=-qstar(1)/utau*zody*vk74in

c ---- check for temperature boundary condition

        if(isfc(1) == 0 ) then
          tsfcc(1)=t1xy(1)-t10xy(1)
        endif

c ----------- input surface wind stress (tau = 0.0184n/m*m)
c             density rho = 1000kg/m^3

        utau = 7.00e-03

c **** save old utau
        utausv = utau
        utau2  = utau*utau
        if (ibuoy==0 .or. qstar(1) == 0.) then
          amonin    = 1000.
          zeta      = 0.
          thstar(1) = 0.0
          t10xy(1)  = 0.0
        else
          amonin = -utau2*utau/(batagk*qstar(1))
          zeta   = z1_a/amonin
        endif
        if (t10xy(1)<0. .and. qstar(1) < 0.) then
          write(6,1234)u1xy,v1xy,t1xy(1),tsfcc(1),amonin,utau,it
 1234     format(' ** check sfc u=',e12.3,' v=',e12.3,' t,ts=',2f10.3,
     +        ' l=',e12.3,' u*=',e12.3,' at it=',i5)
          go to 9999
        endif

c -------- for stable,neutral and unstable pbl get drift velocity

        if(ismlt == 1) then
          call busngr(zeta,phim,phis,psim,psis)
        else
          call fzol(zeta,phim,phis,psim,psis)
        endif
        udrift = windm + stokes(1) - stokess + utau*(zody-psim)*vkin
        vdrift = 0.0
        dnom      = (zody-psis)*vk74in
        if (isfc(1)==1) then
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

c ---------- examples of two other scalars

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

c **** note roundoff problem in angles are close to multiples of pi
        utau2 = utau*utau
        au13m = utau2
        au23m = 0.0
        aut3m(1)= qstar(1)

        return

c --------- trouble in sl routine

 9999   continue

        write(nprt,9000)
 9000   format(' Trouble in SR. sufto')
        call mpi_finalize(ierr)
        stop
      end
      subroutine suft2(u_level1,it)

        use pars
        use fields
        use con_data
        use con_stats

        real u_level1(nnx,iys:iye,2+nscl)

        tol = 0.01
        ufree=0.07*(abs(batag*qstar(1)*dzw(1)))**(1./3.)
        zeta_mn = -6.0
        zeta_mn_i = 1.0/zeta_mn
        iz   = 1

        do iy=iys,iye
          do ix=mxs,mxe

c ----------------- first guess for utau

            utau = .001

            t10xy(1) = -qstar(1)/utau*zody*vk74in
            tsfcc(1) = u_level1(ix,iy,3) - t10xy(1)
            vsfc2    = u_level1(ix,iy,1)**2 + u_level1(ix,iy,2)**2
            vsfc     = sqrt(vsfc2)
            windm    = ufree+vsfc
            utausv   = utau
            utau2    = utau*utau
            amonin   = -utau2*utau/(batagk*qstar(1))
            if(amonin==0.) then
              write(6,5050) ix,iy,it,utau,amonin
 5050         format(' 5050, sr. suft2, trouble at ',/,
     +            ' ix = ',i6,'iy = ',i6,' it = ',i6,' utau = ',e15.6,
     +            ' amonin = ',e15.6)
              stop
            endif

c ---- for unstable, free convection pbl

            iter = 0
 100        continue

c ----------------- limit the min (-l/z) change to accmmodate stable flow

            zeta_i = amin1(amonin/z1,zeta_mn_i)
            zeta_a = 1.0/zeta_i

            if(ismlt == 1) then
              call busngr(zeta_a,phim,phis,psim,psis)
            else
              call fzol(zeta_a,phim,phis,psim,psis)
            endif
            utau     = windm*vk/(zody-psim)
            thstar(1)=-qstar(1)/utau
            amonold  = amonin
            amonin   = utau*utau/(batagk*thstar(1))
            diff     = abs(amonin - amonold)
            iter = iter+1
            if(iter>10)go to 1000
            if(diff>abs(tol*amonin)) go to 100
 1000       continue

 2000       continue

            if (utau>10.) then
              write(6,232)utau,windm
  232         format(' stop because utau=',e15.6,' windm=',e15.6)
              stop 9999
            endif
            t10xy(1) = -qstar(1)/utau*vk74in*(zody-psis)
            t_grnd(ix,iy,1) = u_level1(ix,iy,3) - t10xy(1)

            zol = zeta_a
            hol = zol*zi/z1
            tep = u_level1(ix,iy,1)/windm
            if(tep>1.)  tep = 1.0
            if(tep<-1.) tep = -1.0
            thta  = acos(tep)
            utau2 = utau*utau

            tau13m(ix,iy)   = -utau2*cos(thta)
            tau23m(ix,iy)=-utau2*sin(thta)*sign(1.,u_level1(ix,iy,2))
            taut3m(ix,iy,1) = qstar(1)

c ------- end of x-y loops

          enddo
        enddo

        return
      end
      subroutine init

        use pars
        use particles
        use fields
        use con_data
        use con_stats

        pi   = 4.0*atan(1.0)
        pi2  = 2.0*pi
        bfac = 1.0
        if(ibuoy==0) bfac = 0.

c -------------------- case specific data

        if(iocean == 1) then
          t00b    = 5000.0
          cp      = 4.20e03
          gcp     = grav/cp
          batag   = bfac*grav/t00b
          fcor_h  = 0.0
          vgcont  = 0.
          wtsfc(1)=1.190476e-06
          qstar(1)=wtsfc(1)
          dtdzf(1)=0.2548
          dtjump  = 0.
          divgls  = 0.
          zo      = 0.0001
          zi      = -5.
          izi     = 55
          xl      = 50.
          yl      = 50.
          zl      = -20.

c ---------- if stretched grid specify location of first point

        else
          cp      = 1.e3
          gcp     = grav/cp
          batag   = bfac*grav/t00
          fcor_h  = 0.0

c ----------- wind tunnel simulation

          wtsfc(1)=qstar(1)
          wtsfc(2)=qstar(2)
          dtdzf(1)=0.000
          dtdzf(2)=0.000
          dtjump  = 0.0
          divgls  = 0.0
          izi     = 102

c ------------ check type of grid generation, adjust box size
c              so that you get the desired size in physical space

          if(iz_space == 4) then
            zl = zl + zw1
          endif
        endif

        time     = 0.0
        t_stage  = 0.0
        t_stage2 = 0.0
        t_zero   = 0.0

c ---------- set mesh stationary flag

        imesh = 0

c ---------- outermost coarse grid  indicies are bounds of grid

        izlow = 1
        izup  = nnz
        dz    = zl/nnz
        dzg   = abs(dz)
        if(l_root) write(6,4040) zl,nnz,dzg

c --------------- generate z grids for particular mesh from
c                 iz = 0,1,...,nnz+1; this allows indexing
c                 to array elements z(0), etc.

        zwstrt = 0.0

c ------------ if uniform vertical spacing then

        if(iz_space == 0) then

c ------------ build z grid for w points

          do iz=0,nnz+1
            z(iz) = dz*float(iz) + zwstrt
          enddo
        else
          z_gues = 1.1
        !call vgrid(zw1,zi,zl,z_gues,nnz,z(0),l_root,l_debug)
          call vgrid_channel(zw1,zi,zl,nnz,z(0),l_root,l_debug)
        endif

        call get_dz

        if(l_root) then
          write(6,8002) zwstrt
          write(6,8003) (iz,z(iz),zz(iz),iz=0,nnz+1)
        endif

        nnzm1 = nnz-1
        dx    = xl/nnx
        dy    = yl/nny
        fnxy  = 1./float(nxy)
        dzdz  = dzw(1)*dzw(1)
        z1    = zz(1)

        c23  = 2.0/3.0
        dsl  = (dx*1.5*dy*1.5*abs(dzw(1)))**(1./3.)
        dslg = dsl
        cs   = 0.2

        vk     = 0.4
        batagk = batag*vk
        vkin   = 1./vk
        ttmean = 0.
        zody   = alog(abs(z1/zo))
        write(nprt, 9901) z1,zo,zody
 9901   format(' 9901 z1 = ',e15.6,' zo = ',e15.6,/,
     +      ' zody = ',e15.6)
        zodyin = 1./zody
        wstar  = abs(batag*zi*wtsfc(1))**(1./3.)
        if(ismlt == 1) then

c ---- set constants for businger similarity functions

          vk74   = vk*0.74
          vk74in = 0.74/vk
          zody74 = zody*0.74
        else

c ---- set constants for large similarity functions

          vk74    = vk
          vk74in  = 1.0/vk
          zody74  = zody
        endif
        ugal   = 0.0
c     ugal   = ugcont*0.5
c     ugcont = ugcont - ugal
        cdbtm  = vk*vk/zody/zody
        if(iocean == 1) then
c ----------- set surface friction velocity here and in sr. sufto
c        utau = 4.29e-03
          utau = 7.00e-03
        else
          ufree = 0.07*(abs(batag*qstar(1)*dzw(1)))**(1./3.)

c ---- note : new estimate for utau !!!

          utau  = vk*(ufree+ugcont)/zody
c        utau  = vk*(ufree)/zody
        endif
        utau2    = utau*utau
        if(ibuoy == 0 .or. qstar(1) == 0.) then
          amonin = 1000.0
        else
          amonin = -utau2*utau/(batagk*qstar(1))
        endif
        hol   = abs(zi)/amonin
        zol   = abs(z1)/amonin
        uwsfc = -utau*utau
        vwsfc = -utau*utau

c ------- make sure tsfcc is gt than t00 for both isfc=0 or 1

        if(l_root) then
          write(6,80)
          write(6,2)wtsfc(1),utau,amonin,dtdzf(1),zody,zo
     +        ,cdbtm,ugcont
        endif

        if(l_debug) then
          write(nprt,80)
          write(nprt,2)wtsfc(1),utau,amonin,dtdzf(1),zody,zo
     +        ,cdbtm,ugcont
        endif

        return
c ------------------------
   2    format(10x,' WT =',e12.4,',  U* =',e12.4,',  L =',e12.4,/,
     +      10x,' DTDZ FREE =',e12.4,',  ZODY=',e12.4,/,10x,
     +      ' ZO(BTM) =',e12.4,',  CDBTM=',e12.4,
     +      ',  UG = ',e12.4)
  80    format(///,' ***** SCRATCH RUN ***** ',//)
 4040   format(' zl = ',e15.6,' nnz = ',i5,' dzg = ',e15.6)
 4043   format(' znest = ',e15.6,' nnz = ',i5,' dzg = ',e15.6)
 8002   format(' zwstrt = ',e12.4)
 8003   format(' iz ',5x,' zw',5x,' zu ',5x,/,(i3,2e12.4))
      end
      subroutine cgrid(zl,nnz,z)

c --------- uniform vertical grid

        real z(0:nnz+1)

        ddz = zl/float(nnz)
        do k=1,nnz
          z(k) = float(k)*ddz
        enddo
        z(nnz) = zl
        z(0)   = 0.0
        z(nnz+1) = z(nnz) + (z(nnz) - z(nnz-1))

        return
      end
      subroutine surf_setup_gabls

c ------- set basic parameters of gabls surface and evaluate pieces that are
c         independent of time

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        use wavy_pars

c --------- allocate surface arrays

        allocate(eta0(2,nnx,nny))

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

        return
      end
      subroutine get_surf_gabls(t_bndy)

c ------- generate new gabls surface at t = t_bndy. smoothly add or subtract
c         bumps by varying the amplitude

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        use wavy_pars

        t1_st  = 100.0
        t1_en  = 400.0
        b_rate = 1.0/(t1_en - t1_st)

        if(t_bndy < t1_st) then
          ampl_fac  = 1.0
          ampl_rate = 0.0
        elseif(t_bndy >= t1_st .and. t_bndy < t1_en) then
          ampl_fac  = amax1((1.0 - b_rate*(t_bndy - t1_st)),0.0)
          ampl_rate = -b_rate
        elseif(t_bndy >= t1_en) then
          ampl_fac  = 0.0
          ampl_rate = 0.0
        endif

        do j=iys,iye
          do i=1,nnx
            bndy(i,j)    = ampl_fac*eta0(1,i,j)
            hdot(i,j)    = ampl_rate*eta0(1,i,j)
            w_orbit(i,j) = hdot(i,j)
          enddo
        enddo

        return
      end
      subroutine surf_setup

c ------- set basic parameters of surface and evaluate pieces that are
c         independent of time

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        use wavy_pars

        real amp(nnx,nny), ran_amp(2,nnx,nny)

c --------- allocate surface arrays

        allocate(omega(nnx,nny),
     +      eta0(2,nnx,nny),
     +      s_wrk(2,nnx,iys:iye),
     +      v_wrk(2,nnx,iys:iye),
     +      etax(nnxp2,iys:iye),
     +      etay(nnxp2,iys:iye),
     +      eta_wrk(nny,i2xs:i2xe))

        u10    = 15.0
        age    = 1.2
        angle  = atan(0.0)
        wind_x = cos(angle)
        wind_y = sin(angle)
        dkx    = pi2/xl
        dky    = pi2/yl

c ------------- generate an array of random phases or amplitudes
c               and radial frequencies

        do j=1,nny
          do i=1,nnx
            omega(i,j) = sqrt(sqrt(xkn(i)**2 + ykn(j)**2)*grav)
          enddo
        enddo
        idum = -1
        if(i_ranp == 0) then
          do j=1,nny
            do i=1,nnx
              aa_phase       = pi2*(ran1(idum) - 0.5)
              ran_amp(1,i,j) = cos(aa_phase)
              ran_amp(2,i,j) = sin(aa_phase)
            enddo
          enddo
        else

c ------- normalize amplitudes so that their squared sum = 1

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
 1345     format(' 1345 get surf:',/,
     +        '      sum1 = ',e15.6,' var1 = ',e15.6,/,
     +        '      sum2 = ',e15.6,' var2 = ',e15.6)
        endif

        sum1  = 0.0
        sum2  = 0.0
        p_max = 0.0
        p_min = 0.0
        do j=1,nny
          do i=1,nnx
            sum1 = sum1 + ran_amp(1,i,j)
            sum2 = sum2 + ran_amp(2,i,j)
            if(ran_amp(1,i,j) > p_max) p_max = ran_amp(1,i,j)
            if(ran_amp(1,i,j) < p_min) p_max = ran_amp(1,i,j)
          enddo
        enddo
        sum1 = sum1*fnxy
        sum2 = sum2*fnxy
        write(nprt,5000) p_max, p_min, sum1, sum2
 5000   format(' phase max = ',e15.6,/,
     +      ' phase min = ',e15.6,/,
     +      ' sum r     = ',e15.6,/,
     +      ' sum c     = ',e15.6)

c --------- first construct the pieces independent of t

        fac_a = sqrt(2.0)
        if(i_pm == 1) then

c ---------- pierson-moskowitz spectrum

          call pm64(dkx,dky,xkn,ykn,nnx,nny,u10,amp)

        else

c ---------- donelan spectrum

          call donelan(dkx,dky,xkn,ykn,nnx,nny,u10,age,angle,amp)

        endif

c ------- use a gaussian filter to smooth the wave shape

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

        write(6,8001) sum_spec
 8001   format(' 8001 sum_spec = ',e15.6)

c --------- build coefficients at t = 0

        do j=1,nny
          eta0(1,ncx,j) = 0.0
          eta0(2,ncx,j) = 0.0
        enddo
        do i=1,nnx
          eta0(1,i,ncy) = 0.0
          eta0(2,i,ncy) = 0.0
        enddo

c --------- make surface smooth

        m_cut_x = nnx/3
        m_cut_y = nny/3

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

        return
      end
      subroutine surf_setup_tank

c ------- basic parameters of wavy surface from lab experiments

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        use wavy_pars
        use tank_pars

        real xx(nx_tank), zz_tank(nx_tank), us(nx_tank), ws(nx_tank),
     +      y2(nx_tank), u2(nx_tank), w2(nx_tank), xles(nnx)

c --------- allocate surface arrays

        allocate(z_tank(nnx), u_tank(nnx), w_tank(nnx),
     +      z_tank_f(nnx), u_tank_f(nnx), w_tank_f(nnx),
     +      z_tank_i(nnx), u_tank_i(nnx), w_tank_i(nnx),
     +      hdot_i(nnx), hx_i(nnx), hx2_i(nnx))

        lut    = 89
        close(lut)
        wave_k  = pi2/wave_l

c---------- count number of points in wave file
        open(lut,file=wave_file)
        npts_t = -4 ! don't count top 4 comment lines
        do
          read (lut,*, end=10) ! break at end of file
          npts_t = npts_t + 1
        enddo
10      close (lut)
        open (lut, file = wave_file) !re-open to read

c--------- read in wave file
        do i=1,4
          read(lut,'(a1)') idum
        enddo
        do i=1,npts_t
          read(lut,*) xx(i), zz_tank(i), us(i), ws(i)
        enddo
        if (myid==0) write(*,*) 'Finished reading tank data'

        write(nprt,1001) (i,xx(i), zz_tank(i),us(i),ws(i),
     +      i=1,nnx)
 1001   format(' 1001: tank ',/,
     +      ' i ',5x,' xx ',5x,' zz ',5x,' us ',5x,
     +      ' ws ',/,(i5,4e15.6))

c ---------- fit a spline to the data

        yp1 = 2.0e+30
        ypn = yp1
        call spline(xx,zz_tank,npts_t,yp1,ypn,y2)
        call spline(xx,us,npts_t,yp1,ypn,u2)
        call spline(xx,ws,npts_t,yp1,ypn,w2)

c ---------- interpolate into the spline function at les points

        ddx = xl/float(nnx)
        do i=1,nnx
          xles(i) = float(i-1)*ddx
          call splint(xx,zz_tank,y2,npts_t,xles(i),z_tank(i))
          call splint(xx,us,u2,npts_t,xles(i),u_tank(i))
          call splint(xx,ws,w2,npts_t,xles(i),w_tank(i))
        enddo

        sum = 0.0
        do i=1,nnx
          sum = sum + z_tank(i)
        enddo
        sum = sum/float(nnx)
        do i=1,nnx
          z_tank(i) = z_tank(i) - sum
        enddo

c --------- force mean w = 0 from the data

        sum_w = 0.0
        do i=1,nnx
          sum_w = sum_w + w_tank(i)
        enddo
        sum_w = sum_w/float(nnx)
        do i=1,nnx
          w_tank(i) = w_tank(i) - sum_w
        enddo

        write(nprt,4001) sum, sum_w
 4001   format(' Mean wave height = ',e15.6,' sum w = ',e15.6)

c --------- get fourier coefficients at t = 0
c           to use in building shape at t > 0

        fnx = 1.0/float(nnx)
        do i=1,nnx
          z_tank_f(i) = z_tank(i)*fnx
          u_tank_f(i) = u_tank(i)*fnx
          w_tank_f(i) = w_tank(i)*fnx
        enddo
        call rfftf(nnx,z_tank_f(1),trigx(1,1))
        call rfftf(nnx,u_tank_f(1),trigx(1,1))
        call rfftf(nnx,w_tank_f(1),trigx(1,1))

c ------ dealias the inputs to match the grid resolution
c        careful with the indices, match fftpack order

        ix_cut   = 2*int(float(nnx)/3.) + 2
c     ix_cut   = 44

        do i=ix_cut,nnx
          z_tank_f(i) = 0.0
          u_tank_f(i) = 0.0
          w_tank_f(i) = 0.0
        enddo

        return
      end
      subroutine surf_setup_real

c ------- basic parameters of wavy surface from lab experiments

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        use wavy_pars
        use tank_pars

        real xx(nx_tank), zz_tank(nx_tank), us(nx_tank), ws(nx_tank),
     +      y2(nx_tank), u2(nx_tank), w2(nx_tank), xles(nnx)

c --------- allocate surface arrays

        allocate(z_tank(nnx), u_tank(nnx), w_tank(nnx),
     +      z_tank_f(nnx), u_tank_f(nnx), w_tank_f(nnx),
     +      z_tank_i(nnx), u_tank_i(nnx), w_tank_i(nnx),
     +      hdot_i(nnx), hx_i(nnx), hx2_i(nnx))

        lut    = 89
        close(lut)

c---------- count number of points in wave file
        open(lut,file=wave_file)
        npts_t = -4 ! don't count top 4 comment lines
        do
          read (lut,*, end=10) ! break at end of file
          npts_t = npts_t + 1
        enddo
10      close (lut)
        open (lut, file = wave_file) !re-open to read
       
c--------- read in wave file 
        do i=1,4
          read(lut,'(a1)') idum
        enddo
        do i=1,npts_t
          read(lut,*) xx(i), zz_tank(i), us(i), ws(i)
        enddo
        if (myid==0) write(*,*) 'Finished reading tank data'

        write(nprt,1001) (i,xx(i), zz_tank(i),us(i),ws(i),
     +      i=1,nnx)
 1001   format(' 1001: tank ',/,
     +      ' i ',5x,' xx ',5x,' zz ',5x,' us ',5x,
     +      ' ws ',/,(i5,4e15.6))

c ---------- fit a spline to the data

        yp1 = 2.0e+30
        ypn = yp1
        call spline(xx,zz_tank,npts_t,yp1,ypn,y2)

c ---------- interpolate into the spline function at les points

        ddx = xl/float(nnx)
        do i=1,nnx
          xles(i) = float(i-1)*ddx
          call splint(xx,zz_tank,y2,npts_t,xles(i),z_tank(i))
        enddo

        sum = 0.0
        do i=1,nnx
          sum = sum + z_tank(i)
        enddo
        sum = sum/float(nnx)
        do i=1,nnx
          z_tank(i) = z_tank(i) - sum
        enddo

c --------- zero orbitals for calculation

        do i=1,nnx
          u_tank(i) = 0.0
          w_tank(i) = 0.0
        enddo

        write(nprt,4001) sum, sum_w
 4001   format(' Mean wave height = ',e15.6,' sum w = ',e15.6)

c --------- get fourier coefficients at t = 0
c           to use in building shape at t > 0

        fnx = 1.0/float(nnx)
        do i=1,nnx
          z_tank_f(i) = z_tank(i)*fnx
        enddo
        call rfftf(nnx,z_tank_f(1),trigx(1,1))

c ------ dealias the inputs to match the grid resolution
c        careful with the indices, match fftpack order

        ix_cut   = 2*int(float(nnx)/3.) + 2
c     ix_cut   = 44

        do i=ix_cut,nnx
          z_tank_f(i) = 0.0
        enddo

        return
      end
      subroutine get_surf(t_bndy)

c ------- generate new surface at t = t_bndy
c         assuming its the real part of the complex fft
c         allow for time varying amplitude to smoothly add wave
c         spectrum from a restart

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        use wavy_pars
        real bndyx(nnx,iys:iye), bndyy(nnx,iys:iye,izs:izs)

c ------------ for doing restart with gradually growing bumps
c              specify initial and end times to blend over

        t1_st  = 0.59781897E+04 + 0.5 ! time for volume ja1062
        t1_en  = t1_st + 400.0
        b_rate = 1.0/(t1_en - t1_st)

c --------- for debug turn on waves right away

        if(t_bndy < t1_st) then
          ampl_fac  = 0.0
          ampl_rate = 0.0
        elseif(t_bndy >= t1_st .and. t_bndy < t1_en) then
          ampl_fac  = amax1((b_rate*(t_bndy - t1_st)),0.0)
          ampl_rate = b_rate
        elseif(t_bndy >= t1_en) then
          ampl_fac  = 1.0
          ampl_rate = 0.0
        endif

c -------------- make surface flat, for initialization

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

c --------- get the u_orbital velocity of the wave

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
        if(iys == 1) then
          v_wrk(1,i,iys) = 0.0
          v_wrk(2,i,iys) = 0.0
        endif
        call fft2d_cmplx_mpi(v_wrk(1,1,iys),eta_wrk(1,i2xs),trigcx,
     +      trigc,nnx,nx2,nny,i2xs,i2xe,i2x_s,i2x_e,iys,iye,iy_s,iy_e,1,
     +      1,myid,ncpu_s,numprocs)
        do j=iys,iye
          do i=1,nnx
            u_orbit(i,j) = v_wrk(1,i,j)
          enddo
        enddo

c --------- get the v_orbital velocity of the wave

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
        if(iys == 1) then
          v_wrk(1,i,iys) = 0.0
          v_wrk(2,i,iys) = 0.0
        endif
        call fft2d_cmplx_mpi(v_wrk(1,1,iys),eta_wrk(1,i2xs),trigcx,
     +      trigc,nnx,nx2,nny,i2xs,i2xe,i2x_s,i2x_e,iys,iye,iy_s,iy_e,1,
     +      1,myid,ncpu_s,numprocs)
        do j=iys,iye
          do i=1,nnx
            v_orbit(i,j) = v_wrk(1,i,j)
          enddo
        enddo

c ------------- get hdot

        do j=iys,iye
          do i=2,nnx
            v_wrk(1,i,j) = omega(i,j)*s_wrk(2,i,j)*ampl_fac +
     +          ampl_rate*s_wrk(1,i,j)
            v_wrk(2,i,j) = -omega(i,j)*s_wrk(1,i,j)*ampl_fac +
     +          ampl_rate*s_wrk(2,i,j)
          enddo
        enddo
        i  = 1
        do j=max(2,iys),iye
          v_wrk(1,i,j) = omega(i,j)*s_wrk(2,i,j)*ampl_fac +
     +        ampl_rate*s_wrk(1,i,j)
          v_wrk(2,i,j) = -omega(i,j)*s_wrk(1,i,j)*ampl_fac +
     +        ampl_rate*s_wrk(2,i,j)
        enddo
        if(iys == 1) then
          v_wrk(1,i,iys) = 0.0
          v_wrk(2,i,iys) = 0.0
        endif
        call fft2d_cmplx_mpi(v_wrk(1,1,iys),eta_wrk(1,i2xs),trigcx,
     +      trigc,nnx,nx2,nny,i2xs,i2xe,i2x_s,i2x_e,iys,iye,iy_s,iy_e,1,
     +      1,myid,ncpu_s,numprocs)
        do j=iys,iye
          do i=1,nnx
            hdot(i,j) = v_wrk(1,i,j)
          enddo
        enddo

c -------- get the surface

        call fft2d_cmplx_mpi(s_wrk(1,1,iys),eta_wrk(1,i2xs),trigcx,
     +      trigc,nnx,nx2,nny,i2xs,i2xe,i2x_s,i2x_e,iys,iye,iy_s,iy_e,1,
     +      1,myid,ncpu_s,numprocs)
        do j=iys,iye
          do i=1,nnx
            bndy(i,j)          = s_wrk(1,i,j)*ampl_fac
            bndyx(i,j)         = bndy(i,j)
            bndyy(i,j,izs:izs) = bndy(i,j)
          enddo
        enddo

c -------- get variance of surface

        sumb = 0.0
        do j=iys,iye
          do i=1,nnx
            sumb = sumb + bndy(i,j)
          enddo
        enddo
        sumb = sumb*fnxy
        call mpi_sum_xy(sumb,myid,iss,ise,1)

        varb = 0.0
        do j=iys,iye
          do i=1,nnx
            varb = varb + (bndy(i,j) - sumb)**2
          enddo
        enddo
        varb = varb*fnxy
        call mpi_sum_xy(varb,myid,iss,ise,1)

c --------  w_o = hdot +  u_o*h_x + v_o*h_y

        call xderivp(bndyx(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
        call yd_mpi(bndyy(1,iys,izs),trigx(1,2),yk(1),
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
        do j=iys,iye
          do i=1,nnx
            w_orbit(i,j) = hdot(i,j) + u_orbit(i,j)*bndyx(i,j) +
     +          v_orbit(i,j)*bndyy(i,j,izs)
          enddo
        enddo

c ----------- make sure mean of w_orbit = 0

        sum_w = 0.0
        do j=iys,iye
          do i=1,nnx
            sum_w = sum_w + w_orbit(i,j)
          enddo
        enddo
        sum_w = sum_w*fnxy
        call mpi_sum_xy(sum_w,myid,iss,ise,1)

        do j=iys,iye
          do i=1,nnx
            w_orbit(i,j) = w_orbit(i,j) - sum_w
          enddo
        enddo

        return
      end
      subroutine get_surf_tank(t_bndy)

c ------- generate new surface at t = t_bndy
c         for the laboratory wave. Build w_o and hdot so
c         that they satisfy the inputs for u_o and h
c         with fix sr. hx_deriv

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        use wavy_pars
        use tank_pars

c ------------ for doing restart with gradually growing bumps
c              specify initial and end times to blend over
c              more aggressive strategy to trigger turbulence

c     t1_st  = 0.85340788E+01 + 0.01 ! time from volume ub1020
        t1_st  =  0.0 ! time relative to t_zero set by move_step
c     t1_st  = 100000.0 ! to make surface always flat for initialization

        t1_en  = t1_st + 0.25
        b_rate = 1.0/(t1_en - t1_st)

c --------- for debug turn on waves right away

        if(t_bndy < t1_st) then
          ampl_fac  = 0.0
          ampl_rate = 0.0
        elseif(t_bndy >= t1_st .and. t_bndy < t1_en) then
          ampl_fac  = amax1((b_rate*(t_bndy - t1_st)),0.0)
          ampl_rate = b_rate
        elseif(t_bndy >= t1_en) then
          ampl_fac  = 1.0
          ampl_rate = 0.0
        endif

c ------------ get (u,w) components of surface velocity and wave
c              height. Make h, hdot, uo, wo match

        call gal_interp (u_tank_f, trigx(1,1),xkn, t_bndy, u_tank_i)
        call gal_interp (z_tank_f, trigx(1,1),xkn, t_bndy, z_tank_i)
        call hdot_interp(z_tank_f, trigx(1,1),xkn, t_bndy, hdot_i)
        call hx_deriv   (z_tank_i, trigx(1,1),xkn, nnx, hx_i)
        call hx_deriv   (hx_i,     trigx(1,1),xkn, nnx, hx2_i)

        sum   = 0.0
        sum_w = 0.0
        do i=1,nnx
          sum = sum + z_tank_i(i)
        enddo
        sum = sum/float(nnx)

        ampl_fac = 1.0
        do j=iys,iye
          do i=1,nnx
            u_orbit(i,j) = u_tank_i(i)*ampl_fac
            v_orbit(i,j) = 0.0
            bndy(i,j)    = z_tank_i(i)*ampl_fac
          enddo
        enddo

c -------- get variance of surface

        sumb = 0.0
        do j=iys,iye
          do i=1,nnx
            sumb = sumb + bndy(i,j)
          enddo
        enddo
        sumb = sumb*fnxy
        call mpi_sum_xy(sumb,myid,iss,ise,1)

        varb = 0.0
        do j=iys,iye
          do i=1,nnx
            varb = varb + (bndy(i,j) - sumb)**2
          enddo
        enddo
        varb = varb*fnxy
        call mpi_sum_xy(varb,myid,iss,ise,1)

c -------- use w_o = hdot + u_o*h_x + v_o*h_y

        do j=iys,iye
          do i=1,nnx
            hdot(i,j)    = hdot_i(i)*ampl_fac
            w_orbit(i,j) = hdot(i,j) + u_orbit(i,j)*hx_i(i)*ampl_fac
          enddo
        enddo

c ----------- make sure mean of w_orbit = 0

        sum_w = 0.0
        do j=iys,iye
          do i=1,nnx
            sum_w = sum_w + w_orbit(i,j)
          enddo
        enddo
        sum_w = sum_w*fnxy
        call mpi_sum_xy(sum_w,myid,iss,ise,1)

        do j=iys,iye
          do i=1,nnx
            w_orbit(i,j)    = w_orbit(i,j) - sum_w
          enddo
        enddo

        ! if (myid==0) then
        !   write(*,*) "z_tank_i u_orbit w_orbit"
        !   do i=1,nnx
        !     write(*,*) z_tank_i(i), u_orbit(i,iys), w_orbit(i,iys)
        !   end do
        ! end if

        return
      end
      subroutine get_surf_real(t_bndy)

c ------- generate new surface at t = t_bndy
c         for the laboratory wave. Build w_o and hdot so
c         that they satisfy the inputs for u_o and h
c         with fix sr. hx_deriv

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        use wavy_pars
        use tank_pars

        real gam, k_main, f_main

c ------------ for doing restart with gradually growing bumps
c              specify initial and end times to blend over
c              more aggressive strategy to trigger turbulence

c     t1_st  = 0.85340788E+01 + 0.01 ! time from volume ub1020
        t1_st  =  0.0 ! time relative to t_zero set by move_step
c     t1_st  = 100000.0 ! to make surface always flat for initialization

        t1_en  = t1_st + 0.25
        b_rate = 1.0/(t1_en - t1_st)

c --------- for debug turn on waves right away

        if(t_bndy < t1_st) then
          ampl_fac  = 0.0
          ampl_rate = 0.0
        elseif(t_bndy >= t1_st .and. t_bndy < t1_en) then
          ampl_fac  = amax1((b_rate*(t_bndy - t1_st)),0.0)
          ampl_rate = b_rate
        elseif(t_bndy >= t1_en) then
          ampl_fac  = 1.0
          ampl_rate = 0.0
        endif
c
c ------------ get (u,w) components of surface velocity and wave
c              height. Make h, hdot, uo, wo match

        call gal_interp (z_tank_f, trigx(1,1),xkn, t_bndy, z_tank_i)
        call hdot_interp(z_tank_f, trigx(1,1),xkn, t_bndy, hdot_i)
        call hx_deriv   (z_tank_i, trigx(1,1),xkn, nnx, hx_i)
        call hx_deriv   (hx_i,     trigx(1,1),xkn, nnx, hx2_i)

        sum   = 0.0
        sum_w = 0.0
        do i=1,nnx
          sum = sum + z_tank_i(i)
        enddo
        sum = sum/float(nnx)

        ampl_fac = 1.0

c
c ------------- Calculate u_orbit with semilinear eq 7 from Grue 2008
c       ! note: using u_tank_[fi] as calculation vars
        fnx = 1.0/float(nnx)
        do i=1,nnx
          u_tank_i(i) = hdot_i(i)*fnx  ! pre-normalize as rfftf does not
        enddo
        call rfftf(nnx,u_tank_i,trigx)
        j = 1
        do i=2, nnx-1, 2
          j = j + 1
          u_tank_i(i)   = u_tank_i(i)   / xkn(j)
          u_tank_i(i+1) = u_tank_i(i+1) / xkn(j)
        enddo
        call rfftb(nnx,u_tank_i,trigx)
        call hx_deriv(u_tank_i, trigx, xkn, nnx, u_tank_f)

        do j=iys,iye
          do i=1,nnx
            u_orbit(i,j) = u_tank_f(i)
     +          - (hdot_i(i) + hx_i(i)*u_tank_f(i)) / (1 + hx_i(i)**2)
            v_orbit(i,j) = 0.0
            bndy(i,j)    = z_tank_i(i)*ampl_fac
          enddo
        enddo

c -------- use w_o = hdot + u_o*h_x + v_o*h_y

        do j=iys,iye
          do i=1,nnx
            hdot(i,j)    = hdot_i(i)*ampl_fac
            w_orbit(i,j) = hdot(i,j) + u_orbit(i,j)*hx_i(i)*ampl_fac
          enddo
        enddo

c ----------- make sure mean of w_orbit = 0

        sum_w = 0.0
        do j=iys,iye
          do i=1,nnx
            sum_w = sum_w + w_orbit(i,j)
          enddo
        enddo
        sum_w = sum_w*fnxy
        call mpi_sum_xy(sum_w,myid,iss,ise,1)

        do j=iys,iye
          do i=1,nnx
            w_orbit(i,j)    = w_orbit(i,j) - sum_w
          enddo
        enddo

c ----------- add drift velocity to u_orbit

        gam = 0.5           ! set gamma here
        k_main = 2*pi/2.35  ! set dominant wavelength
        f_main = sqrt(grav*k_main)
        do j=iys,iye
          do i = 1, nnx
            u_orbit(i,j) = u_orbit(i,j)
     +          + gam * (1 + cos(k_main*bndy(i,j) - f_main*t_bndy))
          enddo
        enddo

c -------- get variance of surface

        sumb = 0.0
        do j=iys,iye
          do i=1,nnx
            sumb = sumb + bndy(i,j)
          enddo
        enddo
        sumb = sumb*fnxy
        call mpi_sum_xy(sumb,myid,iss,ise,1)

        varb = 0.0
        do j=iys,iye
          do i=1,nnx
            varb = varb + (bndy(i,j) - sumb)**2
          enddo
        enddo
        varb = varb*fnxy
        call mpi_sum_xy(varb,myid,iss,ise,1)

        return
      end
      subroutine pm64(dkx,dky,xk,yk,nx,ny,u10,amp)

c --------- evaluate PM wave spectrum in terms of wavenumber

        real, parameter ::
     +      grav  = 9.81, f_pm4 = 0.13, a_pm4 = 0.0081, b_pm4 = 1.25
        real kmag
        real amp(nx,ny), xk(nx), yk(ny)

        pi2 = 8.0*atan(1.0)

        sum_dir  = 0.0
        do j=1,ny
          do i=1,nx

c ----------- convert wave number into frequency

            kmag = sqrt(xk(i)**2 + yk(j)**2)
            f    = sqrt(grav*kmag)/pi2

            if(f <= 0.0) then
              amp(i,j) = 0.0
            else
              t1        = a_pm4*grav*grav/(f*(pi2*f)**4)
              fpm       = grav*f_pm4/u10
              t4        = b_pm4*(fpm/f)**4
              amp(i,j)  = t1*exp(-t4)

c ----------- convert into wavenumber

              amp(i,j)  = amp(i,j)*0.5*sqrt(grav/kmag)/(kmag*pi2)

c ----------- simple directional spectrum

              if(xk(i) >= 0.0 .and. yk(j) >= 0.0) then
                angle = acos(xk(i)/kmag)
              elseif(xk(i) >= 0.0 .and. yk(j) <= 0.0) then
                angle =  -acos(xk(i)/kmag)
              elseif(xk(i) <= 0.0 .and. yk(j) >= 0.0) then
                angle = pi2*0.5 - acos(abs(xk(i))/kmag)
              elseif(xk(i) <= 0.0 .and. yk(j) <= 0.0) then
                angle = -pi2*0.5 + acos(abs(xk(i))/kmag)
              endif
              if(angle < -pi2/4.0 .or. angle > pi2/4.0) then
                fac = 0.0
              else
                fac1     = cos(angle)
                fac      = fac1**8
              endif
              amp(i,j) = amp(i,j)*fac
              sum_dir  = sum_dir + fac*dkx*dky
            endif

          enddo
        enddo

        sum_dir_i = 1.0/sum_dir
        sum_spec  = 0.0
        do j=1,ny
          do i=1,nx
            amp(i,j) = amp(i,j)*sum_dir_i
            sum_spec = sum_spec + amp(i,j)*dkx*dky
          enddo
        enddo

        write(6,3000) sum_spec, sum_dir
 3000   format(' 3000: PM 64 sum_spec = ',e15.6,' sum_dir = ',e15.6)

        return
      end
      subroutine donelan(dkx,dky,xk,yk,nx,ny,u10,age,dir_mean,amp)

c -------- get amplitude of wavenumber spectrum according to donelan
c          prescription in Komen etal, p. 187 that includes
c          directional dependence. The routine converts the
c          frequency-directional spectrum F(omega,phi)
c          into a wavenumber-directional spectrum G(k,phi), see notes.

c     Inputs are :

c     (xk,yk) --- (kx,ky) wavenumbers
c     u10     --- wind speed at 10 m
c     age     --- wave age = cp/u10
c     dir_mean -- mean direction of wave field in radians from x-axis

c     Output :

c     amp --- amplitude of spectrum at wavenumber (xk,yk)

        real xk(nx), yk(ny), amp(nx,ny)

        data ionce /0/
        save grav, pi2, ionce
        real kmag

        if(ionce == 0) then
          grav  = 9.81
          pi2   = 8.0*atan(1.0)
          ionce = 1
        endif

c     write(6,5601) u10, age, dir_mean
 5601   format(' 5601: u10 = ',e15.6,' age = ',e15.6,/,
     +      ' dir mean = ',e15.6)

        cp      = age*u10
        omega_p = grav/cp
        age_i   = 1.0/age

        aa      = amin1(age_i,5.0)
        aa      = amax1(aa, 1.0)
        sigma_d = 0.08*(1.0 + 4.0/(aa**3))

        aa      = amin1(age_i,5.0)
        aa      = amax1(aa, 0.83)
        alpha_d = 0.006*(aa**0.55)

        if(age_i <= 1.0) then
          gamma_d = 1.7
        else
          aa = amin1(age_i,5.0)
          gamma_d = 1.7 + 6.0*alog10(age_i)
        endif

        sum_spec = 0.0
        sum_dir  = 0.0

        do j=1,ny
          do i=1,nx

            kmag    = sqrt(xk(i)**2 + yk(j)**2)
            omega   = sqrt(grav*kmag)
            if(kmag == 0.0) then
              amp(i,j) = 0.0
              go to 999
            endif

            if(xk(i) >= 0.0 .and. yk(j) >= 0.0) then
              angle = acos(xk(i)/kmag)
            elseif(xk(i) >= 0.0 .and. yk(j) <= 0.0) then
              angle =  -acos(xk(i)/kmag)
            elseif(xk(i) <= 0.0 .and. yk(j) >= 0.0) then
              angle = pi2*0.5 - acos(abs(xk(i))/kmag)
            elseif(xk(i) <= 0.0 .and. yk(j) <= 0.0) then
              angle = -pi2*0.5 + acos(abs(xk(i))/kmag)
            endif

            gamma = exp(-(omega - omega_p)**2/
     +          (2.0*(sigma_d**2)*omega_p**2))
            fac   = exp(-(omega_p/omega)**4)
            fac1  = gamma_d**gamma
            spec  = alpha_d*(grav**2)*(omega/omega_p)*
     +          fac*fac1/(omega**5)

c --------- directional spreading

            ratio = omega/omega_p
            if(ratio > 0.56 .and. ratio < 0.95) then
              beta = 2.61*ratio**1.3
            elseif(ratio >= 0.95 .and. ratio < 1.6) then
              beta = 2.28/ratio**1.3
            else
              beta = 1.24
            endif

            fac      = (angle - dir_mean)*beta
            sech     = 2.0/(exp(fac) + exp(-fac))
            spec_dir = beta*0.5*sech**2

            amp(i,j)  = spec*spec_dir

c --------- weight amp by the conversion factor to put it
c           into wavenumber space ... see notes

            amp(i,j) = amp(i,j)*0.5*sqrt(grav/kmag)/kmag

c --------- get integrals in wavenumber space

            sum_spec  = sum_spec + amp(i,j)*dkx*dky
            sum_dir   = sum_dir + spec_dir*dkx*dky

  999       continue

          enddo
        enddo

        write(6,6000) sum_spec, sum_dir
 6000   format(' 6000: in Donelan ',/,
     +      ' Sum spectrum = ',e15.6,/,
     +      ' Sum directional spectrum = ',e15.6)

c ---------- rescale the directional
c            spectrum to make sure its integral = 1

        sum_dir_i = 1.0/sum_dir

        return
      end

      real function gasdev(idum)
        integer idum, iset
        real fac, gset, rsq, v1, v2, ran1
        save iset, gset
        data iset /0/
        if(iset == 0) then
    1     continue
          v1  = 2.0*ran1(idum) - 1.0
          v2  = 2.0*ran1(idum) - 1.0
          rsq = v1**2 + v2**2
          if(rsq >= 1.0 .or. rsq == 0.0) go to 1
          fac    = sqrt(-2.0*log(rsq)/rsq)
          gset   = v1*fac
          gasdev = v2*fac
          iset   = 1
        else
          gasdev = gset
          iset   = 0
        endif

        return
      end
      subroutine l_bndy(t_bndy)

c ---------- build the shape of the lower boundary
c            and compute its speed of movement
c            at time t_bndy

        use pars
        use fields
        use con_data
        use con_stats

c -------- pick boundary shape

        ikase = 1 ! wavy bottom with 2 waves
        ikase = 4 ! 3-D hill
        ikase = 2 ! monochromatic wave
        ikase = 6 ! propagating monochromatic wave plus cross swell
        ikase = 7 ! stationary monochromatic wave , ak = 0.5
        ikase = 3 ! 2-D hill
        ikase = 5 ! propagating monochromatic wave

        if(ikase == 1) then ! multiple waves

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
     +            ampl2*sin(fac2*xgrid(i))
              hdot(i,j) = 0.0
              u_orbit(i,j) = 0.0
              v_orbit(i,j) = 0.0
            enddo
          enddo

        elseif(ikase == 2) then ! monochromatic wave

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

        elseif(ikase == 3) then ! 2-D hill

          a      = 25.0
          xcntr  = 0.5*xl
          hill_l = 67.0

          a      = 25.0
          hill_l = 67.0

c --------- flat hill for generating wt flow

c        a      = 0.0

          do j=iys,iye
            do i=1,nnx
              if(xgrid(i) >= (xcntr - 2.0*hill_l) .and.
     +            xgrid(i) <= (xcntr + 2.0*hill_l)) then
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

        elseif(ikase == 4) then ! 3-D hill

          xo    = xl/2.0
          yo    = yl/2.0
          ampl  = 60.000
          rlen  = 1000.0
          do j=iys,iye
            do i=1,nnx
              rad = sqrt((xgrid(i)-xo)**2 + (ygrid(j)-yo)**2)
              if(rad <= rlen*0.5) then
                bndy(i,j) = ampl*(1.0 + cos(pi2*rad/rlen))
              else
                bndy(i,j) = 0.0
              endif
              hdot(i,j) = 0.0
              u_orbit(i,j) = 0.0
              v_orbit(i,j) = 0.0
            enddo
          enddo

        elseif(ikase == 5) then ! propagating monochromatic wave

          wave_l  = xl/12.0
          wave_k  = pi2/wave_l
          ak      = 0.1
          ampl    = ak*wave_l/pi2
          c_speed = sqrt(grav/wave_k)

c ---------- get Cartesian orbital velocities

          f_rad = sqrt(grav*sqrt(wave_k**2))
          con_u = grav*wave_k/f_rad
          do j=iys,iye
            do i=1,nnx
              bndy(i,j) = ampl*cos(wave_k*(xgrid(i) - c_speed*t_bndy))
              hdot(i,j) = ampl*wave_k*c_speed*
     +            sin(wave_k*(xgrid(i) - c_speed*t_bndy))
              u_orbit(i,j) = bndy(i,j)*con_u
              v_orbit(i,j) = 0.0
            enddo
          enddo

        elseif(ikase == 6) then ! propagating monochromatic wave
                                 ! long wave moving vertically

          wave_l  = xl/12.0
          wave_k  = pi2/wave_l
          ak      = 0.1
          ampl    = ak*wave_l/pi2
          c_speed = sqrt(grav/wave_k)

          wave_l2 = xl/6.0
          wave_k2 = pi2/wave_l2
          ak2     = 0.025
          ampl2   = ak2*wave_l2/pi2
          c_spd2  = sqrt(grav/wave_k2)
          do j=iys,iye
            do i=1,nnx
              bndy(i,j) = ampl*cos(wave_k*(xgrid(i) - c_speed*t_bndy))
              hdot(i,j) = ampl*wave_k*c_speed*
     +            sin(wave_k*(xgrid(i) - c_speed*t_bndy))
              bndy2     = ampl2*cos(wave_k2*(ygrid(j) - c_spd2*t_bndy))
              bndy(i,j) = bndy2 + bndy(i,j)
            enddo
          enddo

c ---------- get Cartesian orbital velocities

          f_rad = sqrt(grav*sqrt(wave_k**2 + wave_k2**2))
          con_u = grav*wave_k/f_rad
          con_v = grav*wave_k2/f_rad
          do j=iys,iye
            do i=1,nnx
              u_orbit(i,j) = bndy(i,j)*con_u
              v_orbit(i,j) = bndy(i,j)*con_v
            enddo
          enddo
        elseif(ikase == 7) then ! stationary monochromatic wave
                                 ! long wave moving vertically

          wave_l  = xl/2.0
          wave_k  = pi2/wave_l
          ak      = 0.5

c -------- make it flat

          ak      = 0.0
          ampl    = ak*wave_l/pi2
          c_speed = sqrt(grav/wave_k)
          c_speed = 0.0

          do j=iys,iye
            do i=1,nnx
              bndy(i,j) = ampl*cos(wave_k*(xgrid(i)))
              hdot(i,j) = 0.0
            enddo
          enddo

c ---------- get Cartesian orbital velocities

          f_rad = sqrt(grav*sqrt(wave_k**2))
          con_u = grav*wave_k/f_rad
          con_u = 0.0
          do j=iys,iye
            do i=1,nnx
              u_orbit(i,j) = 0.0
              v_orbit(i,j) = 0.0
            enddo
          enddo
        endif
        return
      end
      subroutine mesh

c ---------- get the mesh of the terrain following or flat grid
c            metrics stored at cell centers
c            should be made more general to avoid grid generation
c            by every processor

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

        integer(kind=mpi_offset_kind) :: offset, disp
        integer(kind=k8)                 nsize, nsize2
        parameter(nvar_g = 4)

        real ztemp(0:nnz+2), fnt1(nnx,iys:iye,izs-1:ize+1)
        real, allocatable, dimension(:,:,:) :: temp
        allocate(temp(nvar_g,nnx,iys:iye))

c ------------ metric definitions

c     gm(.,.,.,1)  = z_zeta = 1/J
c     gm(.,.,.,2)  = zeta_x = - z_xi*J
c     gm(.,.,.,3)  = zeta_y = - z_eta*J
c     gm(.,.,.,4)  = zeta_z = 1/z_zeta = J
c     gm(.,.,.,5)  = z_t = grid speed; zeta_t = - z_t*J

        dx = xl/float(nnx)
        do i=1,nnx
          xgrid(i) = float(i-1)*dx
        enddo
        dy = yl/float(nny)
        do j=1,nny
          ygrid(j) = float(j-1)*dy
        enddo

        if(iflat == 0) then
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

          if(l_debug) then
            ix = nnx/2
            iy = iys
            write(nprt,1020) (k,xgrid(ix),zgrid_u(ix,iy,k),
     +          zgrid_w(ix,iy,k,1), k=izs-1,ize+1)
          endif

        else
          t_bndy = amax1(0.,time - t_zero)

          write(nprt,2001) t_bndy, time, t_zero
 2001     format(' 2001: mesh ',/,
     +        ' t_bndy = ',e15.6,' time = ',e15.6,
     +        ' t_zero = ',e15.6)

          if(i3d_surf == 0) then
            call l_bndy(t_bndy)
          elseif(i3d_surf == 1) then
            call get_surf(t_bndy)
          elseif(i3d_surf == 2) then
            call get_surf_gabls(t_bndy)
          elseif(i3d_surf == 3) then
            call get_surf_tank(t_bndy)
          elseif(i3d_surf == 4) then
            call get_surf_real(t_bndy)
          endif

c -------- initial guess for stretching factor
          z_gues = 1.1

          do j=iys,iye
            do i=1,nnx
              if(iz_space == 1) then
                call vgrid(zw1,zi,zl-bndy(i,j),z_gues,nnz,
     +              ztemp(0),l_root,l_debug)
              elseif(iz_space == 2) then
                z_first = zw1*(1.0 - bndy(i,j)/zl)
                call vgrid(z_first,zi,zl-bndy(i,j),z_gues,nnz,
     +              ztemp(0),l_root,l_debug)
              elseif(iz_space == 3) then
                call exp_grid(bndy(i,j),zl,z(0),ztemp(0),nnz,a_exp_grd)
              elseif(iz_space == 4) then
                call exp_grid2(bndy(i,j),zl,z(0),ztemp(0),nnz,
     +              a_exp_grd2)
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
     +              zgrid_u(i,j,izs)
              else
                zgrid_u(i,j,izs-1) = bndy(i,j) +
     +              0.5*(ztemp(izs-1) + ztemp(izs-2))
              endif

c --------- distance off the surface (every processor)

              z1b(i,j)      = 0.5*(ztemp(1) + ztemp(0))
              bndy_t(i,j,1) = bndy(i,j)
              hdot_t(i,j,1) = hdot(i,j)
            enddo
          enddo

          if(l_debug) then
            ix = nnx/2
            iy = iys
            write(nprt,1020) (k,xgrid(ix),zgrid_u(ix,iy,k),
     +          zgrid_w(ix,iy,k,1),
     +          k=izs-1,ize+1)
          endif

c ------ generate metric components at cell centers
c        set metric components outside domain = to boundary values

          do k=izs,ize+1
            do j=iys,iye
              do i=1,nnx
                gm(i,j,k,1) = (zgrid_w(i,j,k,1) - zgrid_w(i,j,k-1,1))*
     +              dzw_i(k)
              enddo
            enddo
          enddo

c ---------- pass the vertical gradient at ize
c            to the processor above at level izs-1

          nb = myid - ncpu_s
          nt = myid + ncpu_s
          if(iss == 0) then
            nb = mpi_proc_null
          endif
          if(ise == numprocs-1) then
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
     +        fs(1,iys),nsend,mpi_real8,nt,2,
     +        fr(1,iys),nrecv,mpi_real8,nb,2,
     +        mpi_comm_world,istatus,ierr)
          if(iss /= 0) then
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
          if(ize == nnz) then
            k = ize + 1
            j = (iys + iye )/2
            write(nprt,3131) (i,zgrid_u(i,j,k),gm(i,j,k,4),
     +          gm(i,j,k,2),i=1,nnx)
 3131       format('3131 mesh',/,
     +          ' i ',5x,' zu ',5x,' gm 4 ',5x,' gm 2',/,
     +          (i5,3e15.6))
          endif

          do k=izs-1,ize+1
            do j=iys,iye
              do i=1,nnx
                fnt1(i,j,k) = zgrid_u(i,j,k)
              enddo
            enddo
          enddo
          call yd_mpi(fnt1(1,iys,izs-1),trigx(1,2),yk(1),
     +        nnx,nny,ixs,ixe,ix_s,ix_e,
     +        iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
          do k=izs-1,ize+1
            do j=iys,iye
              do i=1,nnx
                gm(i,j,k,3) = -fnt1(i,j,k)*gm(i,j,k,4)
              enddo
            enddo
          enddo

        endif
        if(l_debug) then
          ix = nnx/2
          iy = iys
          write(nprt,1010) (k,gm(ix,iy,k,1),gm(ix,iy,k,2),
     +        gm(ix,iy,k,3),gm(ix,iy,k,4),
     +        gm(ix,iy,k,5),k=izs-1,ize+1)
 1010     format(' ind ',5x,' 1/J',5x,' -z_xi*J',5x,' -z_eta*J',5x,
     +        ' J',5x,' -z_t*J',/,
     +        (i5,5e15.6))
 1020     format(' k',5x,' xgrid',5x,' zgrid_u',5x,'zgrid_w',/,
     +        (i5,3e15.6))
          k  = izs-1
          iy = iys
          write(nprt,1010) (ix,gm(ix,iy,k,1),gm(ix,iy,k,2),
     +        gm(ix,iy,k,3),gm(ix,iy,k,4),
     +        gm(ix,iy,k,5),ix=1,nnx)
          k  = izs-1
          ix = nnx/2
          write(nprt,1010) (iy,gm(ix,iy,k,1),gm(ix,iy,k,2),
     +        gm(ix,iy,k,3),gm(ix,iy,k,4),
     +        gm(ix,iy,k,5),iy=iys,iye)
        endif

c ----------- save the grid file

        path_sav_g = trim(path_sav)//'/'//case(1:3)//'.grid'

        call mpi_file_open(mpi_comm_world, path_sav_g,
     +      mpi_mode_create+mpi_mode_rdwr,
     +      mpi_info_null, ngrd, ierr)

c ---- set file view

        disp = 0
        call mpi_file_set_view(ngrd,disp,mpi_real8,mpi_real8,
     +      'native',mpi_info_null,ierr)

c ---- write data

        nsize   = int(nvar_g,k8)*nnx*nny
        nsize2  = int(nvar_g,k8)*nnx*(iys-1)
        n_write = nvar_g*nnx*(iye+1-iys)

        if(ize == nnz) then
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

          offset = int((k-1),k8)*nsize + nsize2

c -------- mpi i/o with non-uniform blocks of data use _at

          call mpi_file_write_at(ngrd,offset,temp,n_write,
     +        mpi_real8,status,ierr)
          if (ierr /= 0) goto 9991

        enddo

c ---- close file

        call mpi_file_close(ngrd, ierr)

c ---- check file

        if (l_root) then
          inquire(file=path_sav_g,exist=there)
          if(.not.there) then
            write(6,8000) ngrd,myid
            call mpi_abort(mpi_comm_world,icode,ierr)
          endif
          write(6,7000) path_sav_g
        endif

        deallocate(temp)

c ----------- store metrics for possible time variation

        do l=1,5
          do k=izs-1,ize+1
            do j=iys,iye
              do i=1,nnx
                gm2(i,j,k,l) = gm(i,j,k,l)
              enddo
            enddo
          enddo
        enddo

        do k=izs-1,ize+1
          do j=iys,iye
            do i=1,nnx
              zgrid_w(i,j,k,2) = zgrid_w(i,j,k,1)
            enddo
          enddo
        enddo

        do j=iys,iye
          do i=1,nnx
            z1b_2(i,j)    = z1b(i,j)
            bndy_t(i,j,2) = bndy_t(i,j,1)
            hdot_t(i,j,2) = hdot_t(i,j,1)
          enddo
        enddo

        return
c --------------------------  errors in writing grid file
 9991   continue
        write(6,6000) ngrd, k
 6000   format(' SR. MESH:',/,
     +      '    trouble cannot write grid file on unit = ',i2,/,
     +      '             at iz = ',i4)
        call mpi_abort(mpi_comm_world,icode,ierr)
c --------------------
 7000   format(' **** GRID IS WRITTEN IN FILE  ',a80)
 8000   format(' in MESH: trouble writing file ',i5,'  myid = ',i5)
      end
      subroutine new_mesh(it,istage)

c ---------- get a new mesh at time step t_stage3

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        include 'mpif.h'
        logical there
        real fs(nnx,iys:iye), fr(nnx,iys:iye)
        integer istatus(mpi_status_size)

        real ztemp(0:nnz+2), fnt1(nnx,iys:iye,izs-1:ize+1)

c ------------ metric definitions

c     gm2(.,.,.,1)  = z_zeta = 1/J
c     gm2(.,.,.,2)  = zeta_x = - z_xi*J
c     gm2(.,.,.,3)  = zeta_y = - z_eta*J
c     gm2(.,.,.,4)  = zeta_z = 1/z_zeta = J
c     gm2(.,.,.,5)  = z_t = grid speed; zeta_t = - z_t*J

c --------- update the metrics [note gm2(.,.,.,5) = 0]

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

        do k=izs-1,ize+1
          do j=iys,iye
            do i=1,nnx
              zgrid_w(i,j,k,1) = zgrid_w(i,j,k,2)
              zgrid_w(i,j,k,2) = zgrid_w(i,j,k,3)
            enddo
          enddo
        enddo

c --------- distance off the surface

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

c --------- get the new boundary shape
c           at t = t_stage3

        t_bndy = t_stage3 - t_zero

        if(i3d_surf == 0) then
          call l_bndy(t_bndy)
        elseif(i3d_surf == 1) then
          call get_surf(t_bndy)
        elseif(i3d_surf == 2) then
          call get_surf_gabls(t_bndy)
        elseif(i3d_surf == 3) then
          call get_surf_tank(t_bndy)
        elseif(i3d_surf == 4) then
          call get_surf_real(t_bndy)
        endif

c -------- initial guess for stretching factor

        z_gues = 1.1

        do j=iys,iye
          do i=1,nnx
            if(iz_space == 1) then
              call vgrid(zw1,zi,zl-bndy(i,j),z_gues,nnz,
     +            ztemp(0),l_root,l_debug)
            elseif(iz_space == 2) then
              z_first = zw1*(1.0 - bndy(i,j)/zl)
              call vgrid(z_first,zi,zl-bndy(i,j),z_gues,nnz,
     +            ztemp(0),l_root,l_debug)
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
     +            zgrid_u(i,j,izs)
            else
              zgrid_u(i,j,izs-1) = bndy(i,j) +
     +            0.5*(ztemp(izs-1) + ztemp(izs-2))
            endif

c --------- distance off the surface (every processor)

            z1b_3(i,j)    = 0.5*(ztemp(1) + ztemp(0))
            bndy_t(i,j,3) = bndy(i,j)
            hdot_t(i,j,3) = hdot(i,j)
          enddo
        enddo

c ------ generate metric components at cell centers
c        set metric components outside domain = to boundary values

        do k=izs,ize+1
          do j=iys,iye
            do i=1,nnx
              gm3(i,j,k,1) = (zgrid_w(i,j,k,3) - zgrid_w(i,j,k-1,3))*
     +            dzw_i(k)
            enddo
          enddo
        enddo

c ---------- pass the 1/J at ize
c            to the processor above at level izs-1

        nb = myid - ncpu_s
        nt = myid + ncpu_s
        if(iss == 0) then
          nb = mpi_proc_null
        endif
        if(ise == numprocs-1) then
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
     +      fs(1,iys),nsend,mpi_real8,nt,2,
     +      fr(1,iys),nrecv,mpi_real8,nb,2,
     +      mpi_comm_world,istatus,ierr)
        if(iss /= 0) then
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

c ----- compute grid speed gm2(.,.,.,5) in sr. grid_speed

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

        do k=izs-1,ize+1
          do j=iys,iye
            do i=1,nnx
              fnt1(i,j,k) = zgrid_u(i,j,k)
            enddo
          enddo
        enddo
        call yd_mpi(fnt1(1,iys,izs-1),trigx(1,2),yk(1),
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
        do k=izs-1,ize+1
          do j=iys,iye
            do i=1,nnx
              gm3(i,j,k,3) = -fnt1(i,j,k)*gm3(i,j,k,4)
            enddo
          enddo
        enddo

        return
      end
      subroutine start_mesh(it,istage)

c ---------- case with moving mesh
c            initialize results at time step t_stage2

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        include 'mpif.h'
        logical there
        real fs(nnx,iys:iye), fr(nnx,iys:iye)
        integer istatus(mpi_status_size)

        real ztemp(0:nnz+2), fnt1(nnx,iys:iye,izs-1:ize+1)

c ------------ metric definitions

c     gm2(.,.,.,1)  = z_zeta = 1/J
c     gm2(.,.,.,2)  = zeta_x = - z_xi*J
c     gm2(.,.,.,3)  = zeta_y = - z_eta*J
c     gm2(.,.,.,4)  = zeta_z = 1/z_zeta = J
c     gm2(.,.,.,5)  = z_t = grid speed; zeta_t = - z_t*J

c --------- update the metrics [note gm2(.,.,.,5) = 0]

        do l=1,5
          do k=izs-1,ize+1
            do j=iys,iye
              do i=1,nnx
                gm2(i,j,k,l)  = gm(i,j,k,l)
              enddo
            enddo
          enddo
        enddo

        do j=iys,iye
          do i=1,nnx
            z1b_2(i,j)    = z1b(i,j)
            bndy_t(i,j,2) = bndy_t(i,j,1)
            hdot_t(i,j,2) = hdot_t(i,j,1)
          enddo
        enddo

c --------- get the new boundary shape
c           at t = t_stage2

        t_bndy = t_stage2 - t_zero

        if(i3d_surf == 0) then
          call l_bndy(t_bndy)
        elseif(i3d_surf == 1) then
          call get_surf(t_bndy)
        elseif(i3d_surf == 2) then
          call get_surf_gabls(t_bndy)
        elseif(i3d_surf == 3) then
          call get_surf_tank(t_bndy)
        elseif(i3d_surf == 4) then
          call get_surf_real(t_bndy)
        endif

c -------- initial guess for stretching factor

        z_gues = 1.1

        do j=iys,iye
          do i=1,nnx
            if(iz_space == 1) then
              call vgrid(zw1,zi,zl-bndy(i,j),z_gues,nnz,
     +            ztemp(0),l_root,l_debug)
            elseif(iz_space == 2) then
              z_first = zw1*(1.0 - bndy(i,j)/zl)
              call vgrid(z_first,zi,zl-bndy(i,j),z_gues,nnz,
     +            ztemp(0),l_root,l_debug)
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
     +            zgrid_u(i,j,izs)
            else
              zgrid_u(i,j,izs-1) = bndy(i,j) +
     +            0.5*(ztemp(izs-1) + ztemp(izs-2))
            endif

c --------- distance off the surface (every processor)

            z1b_3(i,j)    =  0.5*(ztemp(1) + ztemp(0))
            bndy_t(i,j,3) =  bndy(i,j)
            hdot_t(i,j,3) =  hdot(i,j)
          enddo
        enddo

c ------ generate metric components at cell centers
c        set metric components outside domain = to boundary values

        do k=izs,ize+1
          do j=iys,iye
            do i=1,nnx
              gm3(i,j,k,1) = (zgrid_w(i,j,k,3) - zgrid_w(i,j,k-1,3))*
     +            dzw_i(k)
            enddo
          enddo
        enddo

c ---------- pass the 1/J at ize
c            to the processor above at level izs-1

        nb = myid - ncpu_s
        nt = myid + ncpu_s
        if(iss == 0) then
          nb = mpi_proc_null
        endif
        if(ise == numprocs-1) then
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
     +      fs(1,iys),nsend,mpi_real8,nt,2,
     +      fr(1,iys),nrecv,mpi_real8,nb,2,
     +      mpi_comm_world,istatus,ierr)
        if(iss /= 0) then
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

c ----- compute grid speed gm2(.,.,.,5) in sr. grid_speed

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

        do k=izs-1,ize+1
          do j=iys,iye
            do i=1,nnx
              fnt1(i,j,k) = zgrid_u(i,j,k)
            enddo
          enddo
        enddo
        call yd_mpi(fnt1(1,iys,izs-1),trigx(1,2),yk(1),
     +      nnx,nny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,izs-1,ize+1,myid,ncpu_s,numprocs)
        do k=izs-1,ize+1
          do j=iys,iye
            do i=1,nnx
              gm3(i,j,k,3) = -fnt1(i,j,k)*gm3(i,j,k,4)
            enddo
          enddo
        enddo

        return
      end
      subroutine grid_speed(it,istage)

c ---------- get grid speeds for the mesh at t = t_stage

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        real fs(nnx,iys:iye), fr(nnx,iys:iye), rh(nnz)
        real gt(0:nnz+1,iys:iye,mxs:mxe)
        include 'mpif.h'
        integer istatus(mpi_status_size)

        if(istage <= 2) then
          dtgama_2 = dt*gama(istage+1)
          dtzeta_2 = dt*zetas(istage+1)
        else
          dtgama_2 = dt1*gama(istage-2)
          dtzeta_2 = dt1*zetas(istage-2)
        endif

c ------- get the grid speed at t = t_stage

        gami = 1.0/dtgama_2

        do k=izs,ize
          do j=iys,iye
            do i=1,nnx
              gm2(i,j,k,5) = dzw(k)*gami*(gm3(i,j,k,1) - gm2(i,j,k,1) -
     +            dtzeta_2*rj(i,j,k))
            enddo
          enddo
        enddo
        call xtoz_trans(gm2(1,iys,izs-1,5),gt,nnx,nnz,mxs,mxe,mx_s,mx_e,
     +      iys,iye,izs,ize,iz_s,iz_e,myid,ncpu_s,
     +      numprocs)

c ---------- integrate from bottom up to find grid speeds
c            using surface value as boundary condition
c            hdot comes from boundary routine

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

c --------- transpose back defining ghost point values [izs-1,ize+1]
c           surface value [izs-1] /= value at [izs]
c           with non-uniform spacing of first gridline.

        call ztox_trans(gt,gm2(1,iys,izs-1,5),nnx,nnz,mxs,mxe,mx_s,mx_e,
     +      iys,iye,izs,ize,iz_s,iz_e,myid,ncpu_s,
     +      numprocs)

c ------- get the new rhs of the 1/J equation

        do k=izs,ize
          do j=iys,iye
            do i=1,nnx
              rj(i,j,k) = (gm2(i,j,k,5) - gm2(i,j,k-1,5))*dzw_i(k)
            enddo
          enddo
        enddo

c ---------- determine boundary speeds at t = t_stage
c            for use in psolver

        if(iss == 0) then
          k = 0
          do j=iys,iye
            do i=1,nnx
              wfbc(i,j,3) = gm2(i,j,k,5)
            enddo
          enddo
        endif

        if(ise == numprocs-1) then
          k = nnz
          do j=iys,iye
            do i=1,nnx
              wfbc(i,j,4) = gm2(i,j,k,5)
            enddo
          enddo
        endif

        return
      end
      subroutine algeb_grid(h,zl,zeta,zw,nz,alp)

c ----------- algebraic mapping with grid decay.
c             nicely maps boundary shape as function of(x,y)
c             on to computational coordinates. Use a_decay (alp)
c             to control how fast bumps decay with height.
c             Large values cause non-uniform spacing of
c             first grid point. see hres test code

        real zeta(0:nz+1), zw(0:nz+1)

        do i=0,nz
          zw(i) = zeta(i) + h*(1.0 - zeta(i)/zl)**alp
        enddo
        zw(nz+1) = zw(nz) + (zw(nz) - zw(nz-1))

c --------- fix to match rest of code

        do i=0,nz+1
          zw(i) = zw(i) - h
        enddo

        return
      end
      subroutine exp_grid(h,zl,zeta,zw,nz,alp)

c ----------- simple decaying exponential grid

        real zeta(0:nz+1), zw(0:nz+1)

c --------- choose alp to match physical problem

c     alp = 0.017
c     alp = 0.020

c ----------- tank simulations

c     alp = 66.67
c     alp = 50.00

        do i=0,nz
          zw(i) = zeta(i) + h*exp(-alp*zeta(i))
        enddo
c     zw(nz)   = zeta(nz)
        zw(nz+1) = zw(nz) + (zw(nz) - zw(nz-1))

c --------- fix to match rest of code

        do i=0,nz+1
          zw(i) = zw(i) - h
        enddo

        return
      end
      subroutine exp_grid2(h,zl,zeta,zw,nz,alp)

c ----------- simple decaying exponential grid
c             but with first grid z(1) = zeta(1)
c             works best to pick box height
c             zl(code) = zl(desired) + zeta(1)
c             otherwise top of box in physical space
c             will be shorter by first grid point off the
c             surface

        real zeta(0:nz+1), zw(0:nz+1)

c     alp = 0.020
c     alp = 0.030
c     alp = 0.017
c     alp = 0.00005
c     alp = 0.00010

        zw(0) = h
        do i=1,nz
          d_zeta = zeta(i) - zeta(1)
          zw(i) = d_zeta + (h + zeta(1))*exp(-alp*d_zeta**2)
        enddo
c     zw(nz)   = zeta(nz)
        zw(nz+1) = zw(nz) + (zw(nz) - zw(nz-1))

c --------- fix to match rest of code

        do i=0,nz+1
          zw(i) = zw(i) - h
        enddo

        return
      end
      subroutine vgrid(z1,zi,zl,z_gues,nnz,z,l_root,ldebug)

c ------ iterate sufficiently so that z(nnz) = zl (use r8)

        real z(0:nnz+1)
        logical l_root, l_debug
        data ig_max /500/
        save ig_max

c ----------------- build grid up to zi first

        z_frst = z1
        z_cntr = zl
        n_pbl  = nnz
        z_fac1 = z_cntr/z_frst
        z_fac2 = 1.0/float(nnz)
        z_fac  = z_gues
        knt    = 0
        tol    = 1.0e-10
   10   continue
        knt = knt + 1
        z_facn = (z_fac1*(z_fac - 1.0) + 1.0)**z_fac2
        test   = abs(1.0 - z_facn/z_fac)
        if(knt > ig_max) then
          if(l_root) write(6,9000) z_fac, z_facn, knt
 9000     format(' Cannot find stretching factor',/,
     +        ' z_fac = ',e15.6,' z_facn = ',e15.6,' knt = ',i3)
          stop
        endif
        z_fac = z_facn
        if(test > tol) go to 10
        z(1) = z_frst
        do iz=2,nnz
          z(iz) = z_frst*(z_fac**(float(iz)) - 1.0)/(z_fac - 1.0)
        enddo
        z(nnz)   = zl
        z(0)     = 0.0
        z(nnz+1) = z(nnz) + (z(nnz) - z(nnz-1))

        return
      end
      subroutine get_dz

c --------------- compute spacing for given vertical
c                 point distribution

        use pars
        use fields
        use con_data
        use con_stats
        include 'mpif.h'

        do iz=1,nnz+1
          dzw(iz) = z(iz) - z(iz-1)
        enddo
        dzw(0)     = dzw(1)
        dzw(nnz+2) = dzw(nnz+1)
        do iz=0,nnz+2
          dzw_i(iz) = 1.0/dzw(iz)
        enddo

c ------------ build weights for vertical interpolation
c              of a centered variable to a wf-level

        do iz=0,nnz+1
          weit(iz)  = dzw(iz)/(dzw(iz) + dzw(iz+1))
          weitm(iz) = 1.0 - weit(iz)
        enddo

c ------------ build z grid for u points

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

        return
      end
      subroutine vgrid_channel(z1,zi,zl,nnz,z,l_root,ldebug)

        real z(0:nnz+1)
        integer :: zidx
        logical l_root, l_debug

c ----------------- build grid up to zi first

        z_frst = z1
        z_cntr = zi*0.5
        n_pbl  = nnz
c     n_pbl  = (5*nnz)/8
        z_fac1 = z_cntr/z_frst
        z_fac2 = 1.0/float(n_pbl/2)
        z_fac  = 1.1
        knt = 0
        tol = 0.00001
   10   continue
        knt = knt + 1
        z_facn = (z_fac1*(z_fac - 1.0) + 1.0)**z_fac2
        test   = abs(1.0 - z_facn/z_fac)
        if(knt > 50) then
          if(l_root) write(6,9000) z_fac, z_facn, knt
 9000     format(' Cannot find stretching factor',/,
     +        ' z_fac = ',e15.6,' z_facn = ',e15.6,' knt = ',i3)
          stop
        endif
        z_fac = z_facn
        if(test > tol) go to 10
        if(l_root) write(6,9100) z_fac, z_cntr, z1, knt
 9100   format(' Stretching factor = ',e15.6,/,
     +      ' Match point       = ',e15.6,/,
     +      ' First z           = ',e15.6,/,
     +      ' Number of iters   = ',i4)
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

        if(l_root) write(6,5300) n_pbl
 5300   format(' n_pbl = ',i4)

c -------------- build grid from zi on up
!     For the channel, zi represents the channel centerline
!     Want the mesh to be a mirror image across this:

        z(nnz+1) = z(nnz) + (z(nnz) - z(nnz-1))

        return
      end subroutine vgrid_channel
      subroutine random

c ----------- geostrophic winds designed for comparison case

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        use particles
        real psi(nnx,iys:iye), psix(nnx,iys:iye),
     +      psiy(nnx,iys:iye,izs:izs), uxx(nnx,iys:iye),
     +      vyy(nnx,iys:iye,izs:izs)

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

c ------------ note set nmatch in sr. iso so that
c              it is compatible with conditions here

        do iz=1,nnz
c        ug(iz)   = ugcont*(zz(iz)/zl)
          ug(iz)   = ugcont
          vg(iz)   = vgcont
          divz(iz) = 0.0
        enddo

        izi = 92
        zi  = z(izi)

c -------- alternate strategy for specifying theta

        z_lower = 300.0
        t_lower = 290.0

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

              uf(ix,iy,iz)  = u(ix,iy,iz)
              vf(ix,iy,iz)  = v(ix,iy,iz)
              wf(ix,iy,iz)  = 0.0
            enddo
          enddo
        enddo

c ------------- set initial random field to be
c               divergence free

        idum = -1 - myid
        do iz=izs,ize

c ----------- ampv and ampt are max amplitudes of random
c             velocity and temperature fields
c             make sure ampv is set if free convection so
c             that we have motions at first time step

          ampv = 0.0
          ampt = 0.10

c ------- simple random field scaled between -0.5 and 0.5

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
     +        nnx,nny,ixs,ixe,ix_s,ix_e,
     +        iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
          vmaxx = 0.0
          do iy=iys,iye
            do ix=1,nnx
              vmag = sqrt(psix(ix,iy)**2 + psiy(ix,iy,izs)**2)
              if(vmag > vmaxx) vmaxx = vmag
            enddo
          enddo
          facv = ampv/vmaxx

          if (z(iz) <= 50.0) then
            do iy=iys,iye
              do ix=1,nnx
                u(ix,iy,iz)   = u(ix,iy,iz) - psiy(ix,iy,izs)*facv
                v(ix,iy,iz)   = v(ix,iy,iz) + psix(ix,iy)*facv
                t(ix,iy,1,iz) = t(ix,iy,1,iz) + psi(ix,iy)*ampt
              enddo
            enddo
          endif

          if(z(iz) <= 250.0) then
            do iy=iys,iye
              do ix=1,nnx
                e(ix,iy,iz) = 0.4*(1.0 - z(iz)/250.0)**3
              enddo
            enddo
          endif

c ---------- check divergence of initial field

          do iy=iys,iye
            do ix=1,nnx
              uxx(ix,iy)     = uf(ix,iy,iz)
              vyy(ix,iy,izs) = vf(ix,iy,iz)
            enddo
          enddo
          call xderivp(uxx(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
          call yd_mpi(vyy(1,iys,izs),trigx(1,2),yk(1),
     +        nnx,nny,ixs,ixe,ix_s,ix_e,
     +        iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
          do iy=iys,iye
            do ix=1,nnx
              divz(iz) = divz(iz) + (uxx(ix,iy) + vyy(ix,iy,izs))**2
            enddo
          enddo
          divz(iz) = divz(iz)*fnxy

c -------- end z loop

        enddo

        call mpi_sum_z(divz(1),i_root,myid,nnz,1)

        write(nprt,6000)
 6000   format(' check of divergence for initial state',/,
     +      ' iz ',5x,' divergence')
        write(nprt,6100) (iz,divz(iz),iz=izs,ize)
 6100   format(i5,e15.6)

        return
      end
      subroutine random_wt

c ----------- random conditions for a wind tunnel flow
c             made dimensionless by u_* and wavelength

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        real psi(nnx,iys:iye), psix(nnx,iys:iye),
     +      psiy(nnx,iys:iye,izs:izs), uxx(nnx,iys:iye),
     +      vyy(nnx,iys:iye,izs:izs)

c ------------ note set nmatch in sr. iso so that
c              it is compatible with conditions here

        do iz=1,nnz
c        ug(iz)   = ugcont*(zz(iz)/zl)
          ug(iz)   = ugcont
          vg(iz)   = vgcont
          divz(iz) = 0.0
        enddo

        izi = 100
        zi  = z(izi)

        ustar_gues = 1.0
        ustar_vk   = ustar_gues/vk
        do iz=izs,ize
          do iy=iys,iye
            do ix=1,nnx
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
              uf(ix,iy,iz)  = u(ix,iy,iz)
              vf(ix,iy,iz)  = v(ix,iy,iz)
              wf(ix,iy,iz)  = 0.0
            enddo
          enddo
        enddo

c ------------- set initial random field to be
c               divergence free

        idum = -1 - myid
        do iz=izs,ize

c ----------- ampv and ampt are max amplitudes of random
c             velocity and temperature fields
c             make sure ampv is set if free convection so
c             that we have motions at first time step
c             made dimensionless by ustar

          ampv = 0.01
          ampt = 0.01

c ------- simple random field scaled between -0.5 and 0.5

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
     +        nnx,nny,ixs,ixe,ix_s,ix_e,
     +        iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
          vmaxx = 0.0
          do iy=iys,iye
            do ix=1,nnx
              vmag = sqrt(psix(ix,iy)**2 + psiy(ix,iy,izs)**2)
              if(vmag > vmaxx) vmaxx = vmag
            enddo
          enddo
          facv = ampv/vmaxx

          if (z(iz) <= 0.2) then
            do iy=iys,iye
              do ix=1,nnx
                u(ix,iy,iz)   = u(ix,iy,iz) - psiy(ix,iy,izs)*facv
                v(ix,iy,iz)   = v(ix,iy,iz) + psix(ix,iy)*facv
                t(ix,iy,1,iz) = t(ix,iy,1,iz) + psi(ix,iy)*ampt
              enddo
            enddo
          endif

          if(z(iz) <= 0.2) then
            do iy=iys,iye
              do ix=1,nnx
c           e(ix,iy,iz) = 0.4*(1.0 - z(iz)/zl)**3
                e(ix,iy,iz) = 0.05*(1.0 - z(iz)/zl)**3
              enddo
            enddo
          endif

c ---------- check divergence of initial field

          do iy=iys,iye
            do ix=1,nnx
              uxx(ix,iy)     = uf(ix,iy,iz)
              vyy(ix,iy,izs) = vf(ix,iy,iz)
            enddo
          enddo
          call xderivp(uxx(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
          call yd_mpi(vyy(1,iys,izs),trigx(1,2),yk(1),
     +        nnx,nny,ixs,ixe,ix_s,ix_e,
     +        iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
          do iy=iys,iye
            do ix=1,nnx
              divz(iz) = divz(iz) + (uxx(ix,iy) + vyy(ix,iy,izs))**2
            enddo
          enddo
          divz(iz) = divz(iz)*fnxy

c -------- end z loop

        enddo

        call mpi_sum_z(divz(1),i_root,myid,nnz,1)

        write(nprt,6000)
 6000   format(' check of divergence for initial state',/,
     +      ' iz ',5x,' divergence')
        write(nprt,6100) (iz,divz(iz),iz=izs,ize)
 6100   format(i5,e15.6)

        return
      end
      subroutine random_f

c ---------- example of using given (sparse) initial
c            sounding profiles (FIX for ncpu_s).

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        real psi(nnx,iys:iye), psix(nnx,iys:iye),
     +      psiy(nnx,iys:iye,izs:izs), uxx(nnx,iys:iye),
     +      vyy(nnx,iys:iye,izs:izs)

        parameter (nt=12, nz=11)
        real zg(nz), u_i(nz,nt), v_i(nz,nt), theta_i(nz,nt)
        real ui_temp(nz), vi_temp(nz), ti_temp(nz)
        real time_g(nt)

        data time_g /
     +      0.0000E+00,0.3600E+04,0.7200E+04,0.1080E+05,0.1440E+05,
     +      0.1800E+05,0.2160E+05,0.2520E+05,0.2880E+05,0.3240E+05,
     +      0.3600E+05,  0.3960E+05
     +      /
        data zg /
     +      0.1000E+02,0.3000E+02,0.5500E+02,0.9000E+02,0.1400E+03,
     +      0.2150E+03,0.3300E+03,0.5000E+03,0.7500E+03,0.1100E+04,
     +      0.1600E+04
     +      /
        data u_i /
     +      -0.1510E+01,-0.1560E+01,-0.1580E+01,-0.1580E+01,-0.1560E+01,
     +      -0.1530E+01,-0.1510E+01,-0.9000E+00,-0.1390E+01,-0.1220E+01,
     +      -0.5100E+00,
     +      -0.1090E+01,-0.1110E+01,-0.1120E+01,-0.1120E+01,-0.1030E+01,
     +      -0.9900E+00,-0.9500E+00,-0.6200E+00,-0.1230E+01,-0.9400E+00,
     +      0.2800E+00,
     +      -0.9100E+00,-0.9200E+00,-0.9100E+00,-0.9000E+00,-0.8800E+00,
     +      -0.8400E+00,-0.8000E+00,-0.6500E+00,-0.1510E+01,-0.1070E+01,
     +      0.2400E+00,
     +      -0.8900E+00,-0.8900E+00,-0.8900E+00,-0.8800E+00,-0.8700E+00,
     +      -0.8500E+00,-0.8100E+00,-0.7000E+00,-0.1830E+01,-0.8400E+00,
     +      0.3500E+00,
     +      -0.1250E+01,-0.1260E+01,-0.1260E+01,-0.1250E+01,-0.1240E+01,
     +      -0.1220E+01,-0.1160E+01,-0.8800E+00,-0.1980E+01,-0.1900E+00,
     +      0.7500E+00,
     +      -0.1800E+01,-0.1810E+01,-0.1820E+01,-0.1820E+01,-0.1800E+01,
     +      -0.1780E+01,-0.1710E+01,-0.1150E+01,-0.1960E+01,0.3900E+00,
     +      0.9200E+00,
     +      -0.2110E+01,-0.2130E+01,-0.2140E+01,-0.2140E+01,-0.2130E+01,
     +      -0.2110E+01,-0.2050E+01,-0.9300E+00,-0.1400E+01,0.8800E+00,
     +      0.9600E+00,
     +      -0.2250E+01,-0.2280E+01,-0.2290E+01,-0.2300E+01,-0.2290E+01,
     +      -0.2260E+01,-0.2070E+01,-0.4000E-01,-0.1600E+00,0.1440E+01,
     +      0.1190E+01,
     +      -0.2160E+01,-0.2200E+01,-0.2220E+01,-0.2220E+01,-0.2220E+01,
     +      -0.2190E+01,-0.1610E+01,0.1470E+01,0.1420E+01,0.2050E+01,
     +      0.1610E+01,
     +      -0.2230E+01,-0.2270E+01,-0.2290E+01,-0.2300E+01,-0.2300E+01,
     +      -0.2260E+01,-0.1350E+01,0.2480E+01,0.2380E+01,0.2320E+01,
     +      0.1740E+01,
     +      -0.1890E+01,-0.1930E+01,-0.1950E+01,-0.1950E+01,-0.1940E+01,
     +      -0.1890E+01,-0.1120E+01,0.3010E+01,0.3030E+01,0.2800E+01,
     +      0.2000E+01,
     +      -0.1210E+01,-0.1230E+01,-0.1240E+01,-0.1230E+01,-0.1210E+01,
     +      -0.1140E+01,-0.4600E+00,0.3320E+01,0.3510E+01,0.3420E+01,
     +      0.2340E+01
     +      /
        data v_i /
     +      0.4800E+00,0.5100E+00,0.5300E+00,0.5700E+00,0.6900E+00,
     +      0.7300E+00,0.7600E+00,0.1410E+01,-0.4200E+00,-0.3060E+01,
     +      -0.3500E+01,
     +      0.7800E+00,0.8100E+00,0.8400E+00,0.8900E+00,0.1060E+01,
     +      0.1110E+01,0.1130E+01,0.1190E+01,-0.1040E+01,-0.2900E+01,
     +      -0.3440E+01,
     +      0.3000E+00,0.3200E+00,0.3400E+00,0.3800E+00,0.4800E+00,
     +      0.5300E+00,0.5800E+00,0.5300E+00,-0.1330E+01,-0.2040E+01,
     +      -0.2830E+01,
     +      -0.2700E+00,-0.2600E+00,-0.2400E+00,-0.2200E+00,-0.1800E+00,
     +      -0.1300E+00,-0.5000E-01,0.1000E+00,-0.1170E+01,-0.1100E+01,
     +      -0.2370E+01,
     +      -0.5500E+00,-0.5400E+00,-0.5300E+00,-0.5100E+00,-0.4800E+00,
     +      -0.4100E+00,-0.2600E+00,0.1700E+00,-0.4200E+00,-0.2200E+00,
     +      -0.2080E+01,
     +      -0.2700E+00,-0.2600E+00,-0.2500E+00,-0.2400E+00,-0.2100E+00,
     +      -0.1600E+00,-0.1000E-01,0.8500E+00,0.9700E+00,0.3500E+00,
     +      -0.2250E+01,
     +      0.5300E+00,0.5400E+00,0.5600E+00,0.5700E+00,0.6000E+00,
     +      0.6500E+00,0.7600E+00,0.1960E+01,0.2280E+01,0.3600E+00,
     +      -0.2590E+01,
     +      0.1590E+01,0.1630E+01,0.1650E+01,0.1680E+01,0.1720E+01,
     +      0.1780E+01,0.2010E+01,0.3260E+01,0.3110E+01,0.1600E+00,
     +      -0.2580E+01,
     +      0.2560E+01,0.2620E+01,0.2660E+01,0.2690E+01,0.2740E+01,
     +      0.2830E+01,0.3400E+01,0.4030E+01,0.3030E+01,-0.7000E-01,
     +      -0.2320E+01,
     +      0.3500E+01,0.3600E+01,0.3650E+01,0.3700E+01,0.3750E+01,
     +      0.3860E+01,0.4580E+01,0.4100E+01,0.2450E+01,0.6000E-01,
     +      -0.1770E+01,
     +      0.4500E+01,0.4640E+01,0.4700E+01,0.4760E+01,0.4830E+01,
     +      0.4930E+01,0.5420E+01,0.3960E+01,0.2000E+01,0.5000E+00,
     +      -0.1150E+01,
     +      0.5290E+01,0.5470E+01,0.5550E+01,0.5620E+01,0.5690E+01,
     +      0.5790E+01,0.6070E+01,0.4000E+01,0.1910E+01,0.9700E+00,
     +      -0.5600E+00
     +      /
        data theta_i /
     +      0.2936E+03,0.2936E+03,0.2937E+03,0.2937E+03,0.2938E+03,
     +      0.2942E+03,0.2948E+03,0.2980E+03,0.3027E+03,0.3092E+03,
     +      0.3186E+03,
     +      0.2937E+03,0.2937E+03,0.2937E+03,0.2938E+03,0.2939E+03,
     +      0.2942E+03,0.2946E+03,0.2978E+03,0.3024E+03,0.3090E+03,
     +      0.3184E+03,
     +      0.2939E+03,0.2939E+03,0.2939E+03,0.2939E+03,0.2939E+03,
     +      0.2941E+03,0.2944E+03,0.2976E+03,0.3023E+03,0.3089E+03,
     +      0.3182E+03,
     +      0.2940E+03,0.2939E+03,0.2939E+03,0.2939E+03,0.2939E+03,
     +      0.2941E+03,0.2943E+03,0.2975E+03,0.3022E+03,0.3087E+03,
     +      0.3181E+03,
     +      0.2940E+03,0.2940E+03,0.2939E+03,0.2939E+03,0.2939E+03,
     +      0.2940E+03,0.2942E+03,0.2974E+03,0.3021E+03,0.3086E+03,
     +      0.3180E+03,
     +      0.2941E+03,0.2940E+03,0.2940E+03,0.2940E+03,0.2941E+03,
     +      0.2941E+03,0.2941E+03,0.2973E+03,0.3019E+03,0.3085E+03,
     +      0.3179E+03,
     +      0.2942E+03,0.2942E+03,0.2942E+03,0.2942E+03,0.2941E+03,
     +      0.2941E+03,0.2941E+03,0.2973E+03,0.3020E+03,0.3086E+03,
     +      0.3179E+03,
     +      0.2943E+03,0.2943E+03,0.2943E+03,0.2943E+03,0.2943E+03,
     +      0.2943E+03,0.2943E+03,0.2975E+03,0.3022E+03,0.3087E+03,
     +      0.3181E+03,
     +      0.2945E+03,0.2945E+03,0.2945E+03,0.2945E+03,0.2945E+03,
     +      0.2944E+03,0.2946E+03,0.2978E+03,0.3025E+03,0.3090E+03,
     +      0.3184E+03,
     +      0.2947E+03,0.2947E+03,0.2947E+03,0.2947E+03,0.2946E+03,
     +      0.2946E+03,0.2949E+03,0.2980E+03,0.3027E+03,0.3093E+03,
     +      0.3187E+03,
     +      0.2949E+03,0.2949E+03,0.2949E+03,0.2948E+03,0.2948E+03,
     +      0.2948E+03,0.2950E+03,0.2982E+03,0.3028E+03,0.3094E+03,
     +      0.3188E+03,
     +      0.2950E+03,0.2950E+03,0.2950E+03,0.2950E+03,0.2950E+03,
     +      0.2950E+03,0.2950E+03,0.2982E+03,0.3029E+03,0.3095E+03,
     +      0.3188E+03
     +      /

        save time_g, zg, u_i, v_i, theta_i

c --------- find time location of initial profiles

        call lterp(nt,time_g,t_factor,jt,jtp1,t_weit)

        do iz=1,nz
          ui_temp(iz) = u_i(iz,jt)*(1.0 - t_weit) +
     +        u_i(iz,jtp1)*t_weit
          vi_temp(iz) = v_i(iz,jt)*(1.0 - t_weit) +
     +        v_i(iz,jtp1)*t_weit
          ti_temp(iz) = theta_i(iz,jt)*(1.0 - t_weit) +
     +        theta_i(iz,jtp1)*t_weit
        enddo

c ----------- interpolate vertically

        do iz=izs,ize
          call lterp(nz,zg,zz(iz),kk,kkp1,v_weit)
          u_temp = ui_temp(kk)*(1.0 - v_weit) +
     +        ui_temp(kkp1)*v_weit
          v_temp = vi_temp(kk)*(1.0 - v_weit) +
     +        vi_temp(kkp1)*v_weit
          theta_temp = ti_temp(kk)*(1.0 - v_weit) +
     +        ti_temp(kkp1)*v_weit

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

c ------------- set initial random field to be
c               divergence free

        idum = -1
        do iz=izs,ize
          if (iz<=8) then

c ----------- ampv and ampt are max amplitudes of random
c             velocity and temperature fields

            ampv = 0.5
            ampt = 0.1

c ------- simple random field scaled between 0 and 1

            do iy=iys,iye
              do ix=1,nnx
                psi(ix,iy) = ran1(idum)
              enddo
            enddo

            do iy=iys,iye
              do ix=1,nnx
                psix(ix,iy) = psi(ix,iy)
                psiy(ix,iy,izs) = psi(ix,iy)
              enddo
            enddo
            call xderivp(psix(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
            call yd_mpi(psiy(1,iys,izs),trigx(1,2),yk(1),
     +          nnx,nny,ixs,ixe,ix_s,ix_e,
     +          iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)

            vmaxx = 0.0
            do iy=iys,iye
              do ix=1,nnx
                vmag = sqrt(psix(ix,iy)**2 + psiy(ix,iy,izs)**2)
                if(vmag > vmaxx) vmaxx = vmag
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

c ---------- check divergence of initial field

          do iy=iys,iye
            do ix=1,nnx
              uxx(ix,iy)     = u(ix,iy,iz)
              vyy(ix,iy,izs) = v(ix,iy,iz)
            enddo
          enddo
          call xderivp(uxx(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
          call yd_mpi(vyy(1,iys,izs),trigx(1,2),yk(1),
     +        nnx,nny,ixs,ixe,ix_s,ix_e,
     +        iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
          do iy=iys,iye
            do ix=1,nnx
              divz(iz) = divz(iz) + (uxx(ix,iy) + vyy(ix,iy,izs))**2
            enddo
          enddo
          divz(iz) = divz(iz)*fnxy

c -------- end z loop

        enddo

        call mpi_sum_z(divz(1),i_root,myid,nnz,1)

        write(nprt,6000)
 6000   format(' check of divergence for initial state',/,
     +      ' iz ',5x,' divergence')
        write(nprt,6100) (iz,divz(iz),iz=izs,ize)
 6100   format(i5,e15.6)

        return
      end
      subroutine randoc

c -------- random initial conditions for an
c          ocean simulation

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        real psi(nnx,iys:iye), psix(nnx,iys:iye),
     +      psiy(nnx,iys:iye,izs:izs), uxx(nnx,iys:iye),
     +      vyy(nnx,iys:iye,izs:izs)

        izi=(5*nnz)/20
        zi=z(izi)
        tmixed = 283.0
        do iz=izs,ize
          if (iz<=izi) then
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
          if (iz>izi) then
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

c ------------- set initial random field to be
c               divergence free

        idum = -1
        do iz=izs,ize
          if (iz<=4) then

c ----------- ampv and ampt are max amplitudes of random
c             velocity and temperature fields

            ampv = 0.01
c        ampt = 0.00
            ampt = 0.0001

c ------- simple random field scaled between 0 and 1

            do iy=iys,iye
              do ix=1,nnx
                psi(ix,iy) = ran1(idum)
              enddo
            enddo

            do iy=iys,iye
              do ix=1,nnx
                psix(ix,iy) = psi(ix,iy)
                psiy(ix,iy,izs) = psi(ix,iy)
              enddo
            enddo
            call xderivp(psix(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
            call yd_mpi(psiy(1,iys,izs),trigx(1,2),yk(1),
     +          nnx,nny,ixs,ixe,ix_s,ix_e,
     +          iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)

            vmaxx = 0.0
            do iy=iys,iye
              do ix=1,nnx
                vmag = sqrt(psix(ix,iy)**2 + psiy(ix,iy,izs)**2)
                if(vmag > vmaxx) vmaxx = vmag
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

c ---------- check divergence of initial field

          do iy=iys,iye
            do ix=1,nnx
              uxx(ix,iy) = u(ix,iy,iz)
              vyy(ix,iy,izs) = v(ix,iy,iz)
            enddo
          enddo
          call xderivp(uxx(1,iys),trigx(1,1),xk(1),nnx,iys,iye)
          call yd_mpi(vyy(1,iys,izs),trigx(1,2),yk(1),
     +        nnx,nny,ixs,ixe,ix_s,ix_e,
     +        iys,iye,iy_s,iy_e,izs,izs,myid,ncpu_s,numprocs)
          do iy=iys,iye
            do ix=1,nnx
              divz(iz) = divz(iz) + (uxx(ix,iy) + vyy(ix,iy,izs))**2
            enddo
          enddo
          divz(iz) = divz(iz)*fnxy

c -------- end z loop

        enddo

        call mpi_sum_z(divz(1),i_root,myid,nnz,1)

        write(nprt,6000)
 6000   format(' check of divergence for initial state',/,
     +      ' iz ',5x,' divergence')
        write(nprt,6100) (iz,divz(iz),iz=izs,ize)
 6100   format(i5,e15.6)

        do iz=izs,ize
          ug(iz)=ugcont
          vg(iz)=vgcont
        enddo

        return
      end
      subroutine forcing

c ----------- update surface temperature based on a
c             constant cooling rate

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        include 'mpif.h'

        t_surface = t_surf_i - c_rate*time

        do j=iys,iye
          do i=1,nnx
            t_grnd(i,j,1) = t_surface
          enddo
        enddo

        return
      end
      subroutine set_sav(it,istart)

c --------- time stepping layout

c ----------------------------------------------------------------------
c      it-1          it             it+1
c       0 ----------> 1               2
c       *             *               *
c             dt             dt1             dt_new
c     time
c ----------------------------------------------------------------------

c                it-1 = 1 ---------> it=2           3
c                     *               *             *
c                         dt = dt1       dt1 = dt_new
c                                                update dt_new from it = 1
c                    time
c ----------------------------------------------------------------------

c --------- dt  = current increment between [it-1,it]
c           dt1 = next increment between [it,it+1]
c           dt_new from flow field at end of current it+1 step

        use pars
        use fields
        use con_data
        use con_stats
        use netcdf_io

        data ionce /0/
        save ionce

        if(it /= istart) then

c ---------- increment time if not first time through

          time=time+dt
        endif

        if(istart == 0) then
          dt1 = dt_new
        endif

        it=it+1

        dt    = dt1
        dt1   = dt_new
        mnout = (mod(it,imean)==0).or. (it==1)
        mtape = (mod(it,itape)==0)
        micut = (mod(it,itcut)==0)
        mviz  = (mod(it,iviz)==0)

c ---------- decide frequency history fields are saved

        if(ihst < 0) then
          mhis = .false.
        else
          mhis = (mod(it,ihst)==0 .and. it >= it_his)
        endif
        mtrans = (mod(it+1,ihst)==0)

c ---------- decide frequency viz fields are saved

        if(iviz < 0) then
          msave_v = .false.
        else
          msave_v = (mod(it,iviz)==0 .and. it >= it_viz)
        endif

c ---------- decide whether velocity fields are saved

        msave = .false.
        if(it >= itstr .and. mtape) then
          msave = .true.
          call get_output_filenames
          itn=itn+1
        endif

c --------- decide whether history files are to be saved

        if((ihst > 0) .and. (it >= it_his_nxt)) then
          if (inetcdf == 1) then
            call open_his_netcdf
          else
            call open_his(it)
          endif
          it_his_nxt = it_his_nxt + itape
        endif

c --------- decide whether viz files are to be saved

        if((iviz > 0) .and. (it >= it_viz_nxt)) then
          call viz_output_filename(it)
          it_viz_nxt = it_viz_nxt + itape
        endif

c ---------- if first time thru and moving mesh
c            setup the grid and set t_zero = reference time

        if(move_step < 0 .or.
     +      it < move_step) then
          imesh  = 0
          t_zero = time
          go to 999
        endif

        if(it == move_step .and. imesh == 0) then

c -------- first time thru set W = grid speed everywhere

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

        elseif(it > move_step .and. imesh == 0) then

c ---------- restart with moving grid

          imesh    = 1
          t_stage  = time + dt*etas(1)
          t_stage2 = time + dt*etas(2)
          call start_mesh(it,1)
          call grid_speed(it,0)
        endif

  999   continue

        return
      end
      subroutine pbltop(itop)

c ---------- get estimate of pbl top

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

c ------------ if method uses average statistics then only root
c              process need find zi

        use pars
        use fields
        use con_data
        use con_stats
        real trun(maxnz)
        include 'mpif.h'
        real gradloc(2,nnx,nny), gradmax(2,nnx,nny)
        external get_zi

        if(method <= 2 .and. l_root) then

          sgn = 1.0
          if(iocean == 1) sgn = 1.0
          if (method <= 0 .or. method > 2) then
            itop=1
            wttot=wtle(1,1)+wtsb(1,1)
            wtmin=wttot*sgn
            do iz=2,nnz
              wttot=(wtle(iz,1)+wtsb(iz,1))*sgn
              if (wttot<=wtmin) then
                itop=iz
                wtmin=wttot
              endif
            enddo
            zi=z(itop)
          elseif (method == 1) then
            itop = 1
            crit = 0.05
            uwsf = utau*utau
            do iz=1,nnzm1
              uwtot = (uwle(iz) + uwsb(iz))**2 +
     $            (vwle(iz) + vwsb(iz))**2
              uwtot = sqrt(uwtot)
              if(uwtot/uwsf > crit) then
                itop=iz
              endif
            enddo
            zi=z(itop)
          elseif (method == 2) then
            trun(1) = txym(1,1)
            do iz=2,nnz
              weight = z(iz-1)/z(iz)
              trun(iz) = trun(iz-1)*weight + (1.0-weight)*txym(iz,1)
            enddo
            itop = 1
            tcrit = 0.25
            if(iocean == 1) tcrit = 0.1
            do iz=2,nnz
              if(txym(iz,1) > (trun(iz) + tcrit)) then
                itop = iz
                go to 320
              endif
            enddo
  320       continue
            zi=z(itop)
          endif
          do iy=1,nny
            do ix=1,nnx
              gradmax(2,ix,iy) = zi
            enddo
          enddo

c ----------- use gradient method, every process computes

        elseif(method == 3) then

c ------- similar to zeroing the stat array in sr. mean_stat

          do iy=1,nny
            do ix=1,nnx
              gradloc(1,ix,iy) = 0.0
              gradloc(2,ix,iy) = z(iz_min)
            enddo
          enddo

c ------------- now all z in this process

          if(iz_min <= ize) then
            do iz=max(izs,iz_min),ize
              izp1 = iz + 1
              do iy=iys,iye
                do ix=1,nnx
                  grad = (t(ix,iy,1,izp1) - t(ix,iy,1,iz))*dzu_i(izp1)
                  if(grad > gradloc(1,ix,iy)) then
                    gradloc(1,ix,iy) = grad
                    gradloc(2,ix,iy) = z(iz)
                  endif
                enddo
              enddo
            enddo
          endif

c ----------- alternate version using already defined function in mpi
c             passes 2 real8 variables

          call mpi_reduce(gradloc,gradmax,nnx*nny,mpi_2double_precision,
     +        mpi_maxloc,i_root,mpi_comm_world,ierror)

c ------------ get average on root process

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
              if(zi > z(k) .and. zi <= z(k+1)) then
                itop = k
              endif
            enddo
          endif

        endif

c -------- send average zi everywhere

        call mpi_bcast(zi,1,mpi_real8,
     +      i_root,mpi_comm_world,ierr)

        if(iocean /= 1) then
          do iz=1,nnz
            if(zi >= z(iz) .and.
     +          zi < z(iz+1)) itop = iz
          enddo
        else
          do iz=1,nnz
            if(zi <= z(iz) .and.
     +          zi > z(iz+1)) itop = iz
          enddo
        endif

c ---------- get zi/L

        hol = zi/amonin

c     if(l_root) write(6,7001) myid,zi,itop
 7001   format(' 7001 in pbltop myid = ',i4,' zi = ',e15.6,
     +      ' itop = ',i3)

        return
      end
      subroutine get_zi(gradmax,gradout,len,itype)

        use pars
        real gradmax(*), gradout(*)

        do i=1,len,2
          if(gradmax(i) > gradout(i)) then
            gradout(i)   = gradmax(i)
            gradout(i+1) = gradmax(i+1)
          endif
        enddo

        return
      end
      subroutine print(lu,it,iz_strt,iz_end)

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats

        write(lu,4000)
 4000   format(30X,' --- SOLUTION ---')
        write(lu,4100) it,time,dt,zi,tsfcc(1),uwsfc,vwsfc,wtsfc(1),
     +      zol, hol, amonin, ucfl, vcfl, wcfl, utau,
     +      t_ref, drg_x, drg_y, p_iter, p_err
 4100   format(' IT=',I7,5x,'TIME (s) = ',e15.8,',  DT(s) = ',e15.6,/,
     +      10x,'ZTOP = ',e15.6,
     +      ',  TSFC = ',e15.6,
     +      ',  UW = ',e15.6,',  VW = ',e15.6,/,10x,
     +      'WT = ',e15.6,',  Z1/L =',e15.6,
     +      ',  ZI/L = ',e15.6,',  MO-L = ',e15.6,/,10x,
     +      'U_cfl = ',e15.6,
     +      ',  V_cfl = ',e15.6,',  W_cfl = ',e15.6,
     +      ',  UTAU = ',e15.6,/,10x,
     +      'Theta Ref = ',e15.6,/,10x,'P-Drag_x  = ',e15.6,
     +      ',   P-Drag_y = ',e15.6,/,10x,
     +      'P iterations = ',e15.6,',   P error = ',e15.6)
        write(lu,4200)
 4200   format(//,20x,'--------- HORIZONTAL MEAN VALUES ---------- ',
     +      //,2x,'IZ',4x,'T_MEAN',7x,
     +      'DIVG',8X,'LE_KE',6X,'SGS_KE',7X,'LE_WT',6X,
     +      'SGS_WT',7X,'SHRZ',8X,'BUOY')
        do 19999 iz=iz_end,iz_strt,-1
          write(lu,4300)iz,txym(iz,1)-t_ref,divz(iz),
     +        englez(iz),eavg(iz),wtle(iz,1),
     +        wtsb(iz,1),shrz(iz),buyz(iz)
 4300     format(1X,I3,e12.4,7e12.4)
19999   continue
        write(lu,4400)tsfcc(1),wtsfc(1)
 4400   format('  SURFACE VALUES: TXYM = ',E15.6,'  WTSB = ',E15.6)
        if(iocean == 1) then
          write(lu,4500) stokess,udrift,vdrift
 4500     format(/,' STOKESS = ',e12.4,' UDRIFT = ',e12.4,
     +        ' VDRIFT = ',e12.4)
        endif
        write(lu,4600) (iz,uxym(iz)+ugal,vxym(iz),uwle(iz),
     +      uwsb(iz),vwle(iz),vwsb(iz),pm(iz),iz=iz_strt,iz_end)
 4600   format(//,' IZ',5x,' UXYM + UGAL',8x,' VXYM',10x,' UWLE',11x,
     +      ' UWSB',11x,' VWLE',11x,' VWSB',11x,' P MEAN'
     +      ,/,(1x,i4,7(3x,e13.6)))
        if(ivis == 1) then
          write(lu,4800) xksurf, nmatch, viscon, vise
 4800     format(//,' XKSURF = ',e15.6,' NMATCH = ',i4,/,
     +        ' VISCON = ',e15.6,' VISE = ',e15.6)
          write(lu,4700) (iz,dfac(iz),iz=iz_strt,iz_end)
 4700     format(//,'   IZ',5x,'  DFAC',/,(1x,i4,3x,e15.6))
        endif

c --------------- output additional scalars

        if(nscl == 2) then
          write(lu,5005)tsfcc(2),wtsfc(2)
 5005     format(/,'  SURFACE VALUE: TXYM(2) =',e15.6,' WTSFC(2) = ',
     +        e15.6)
          write(lu,5100) (iz,txym(iz,2),wtle(iz,2),
     +        wtsb(iz,2),iz=iz_strt,iz_end)
 5100     format(//,' IZ',5x,' SCALAR-1 MEAN',8x,' WS1LE',10x,
     +        ' WS1SB',10x
     +        ,/,(1x,i4,3(3x,e13.6)))
        endif

        return
      end
      subroutine xy_stats

c ------------ get statistics.
c              fluxes that appear in transformed equations
c              are based on (UF,VF,WF-Z_t) located at cell centers
c              and upper cell face.

c              f_avg = [\int \int f * fac dx dy]/area
c              fac   = [sqrt(zeta_x^2 + zeta_y^2 + zeta_z^2]/J

        use pars
        use fields
        use con_data
        use con_stats
        use particles

c ------- indices for indexing array stat(.,.)
c         js = number of non-scalar stats
c         ns = number of scalar stats

        parameter(js = 38, ns = 6, nstat = js + ns*nscl)
        real stat(1:nnz,nstat), buf(4), RHtmp, mod_magnus

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
c          stat(.,38) = mean RH
c          stat(.,m1) = resolved scalar flux w_ft at zw = wtle
c          stat(.,m2) = resolved scalar flux u_ft at zu = utle
c          stat(.,m3) = resolved scalar flux v_ft at zu = vtle
c          stat(.,m4) = scalar t*t at zu = tps
c          stat(.,m5) = scalar t*t*t at zu = tcube
c          stat(.,m6) = resolved scalar flux wct at zu = wctle

c --------- use a trick with mpi reduce over all z to get averages
c           by setting stat array = 0 for all z on each process

        do i=1,nstat
          do k=1,nnz
            stat(k,i) = 0.0
          enddo
        enddo

c -------- indices for scalars

        m1 = js
        m2 = js + nscl
        m3 = js + 2*nscl
        m4 = js + 3*nscl
        m5 = js + 4*nscl
        m6 = js + 5*nscl

        sgn = 1.0
        if(iocean == 1 .and. iupwnd == 1) sgn = -1.0

        do k=izs,ize

          kp2 = k + 2
          kp1 = k + 1
          km1 = k - 1

c ------- simple average along coordinate lines

          do j=iys,iye
            do i=1,nnx
              stat(k,1) = stat(k,1) + (u(i,j,k) - uxym(k))**2
              stat(k,2) = stat(k,2) + (v(i,j,k) - vxym(k))**2
              stat(k,3) = stat(k,3) + (w(i,j,k) - wxym(k))**2
              stat(k,4) = stat(k,4) + (w(i,j,k) - wxym(k))**3
              stat(k,5) = stat(k,5) + (w(i,j,k) - wxym(k))**4
              stat(k,6) = stat(k,6) +
     +            ((w(i,j,k)- wxym(k))**2 +
     +            (u(i,j,k)- uxym(k))**2 +
     +            (v(i,j,k)- vxym(k))**2)*0.5

              stat(k,7) = stat(k,7) + e(i,j,k)
              stat(k,8) = stat(k,8) + e(i,j,k)
              stat(k,9) = stat(k,9) + 0.5*(wf(i,j,k) - gm(i,j,k,5))*
     +            (u(i,j,k)-uxym(k) + u(i,j,kp1) - uxym(kp1))
              stat(k,10) = stat(k,10) + 0.5*(wf(i,j,k) - gm(i,j,k,5))*
     +            (v(i,j,k)-vxym(k) + v(i,j,kp1) - vxym(kp1))
              stat(k,11) = stat(k,11) + (w(i,j,k)-wxym(k))*
     +            (u(i,j,k)-uxym(k))
              stat(k,12) = stat(k,12) + (w(i,j,k)-wxym(k))*
     +            (v(i,j,k)-vxym(k))

c --------- get drag [p(z)*z_xi(z), p(z)*z_eta(z)]

              stat(k,13) = stat(k,13) + p(i,j,k)*gm(i,j,k,1)*gm(i,j,k,2)
              stat(k,14) = stat(k,14) + p(i,j,k)*gm(i,j,k,1)*gm(i,j,k,3)

c --------- get drag [p(z)*z_xi(1), p(z)*z_eta(1)]

              stat(k,15) = stat(k,15) + p(i,j,k)*gms(i,j,1)*gms(i,j,2)
              stat(k,16) = stat(k,16) + p(i,j,k)*gms(i,j,1)*gms(i,j,3)

              stat(k,17) = stat(k,17) + 0.5*(wf(i,j,k) - gm(i,j,k,5))*
     +            (w(i,j,k)-wxym(k) + w(i,j,kp1) - wxym(kp1))

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

              RHtmp = t(i,j,1,k)*t(i,j,2,k)*
     +            Ru/Mw/mod_magnus(t(i,j,1,k))*rhoa*100.0

              stat(k,38) = stat(k,38) + RHtmp
            enddo
          enddo

c ------------ get scalar resolved fluxes and variances

          do l=1,nscl
            if(iupwnd /= 1 .or. k == nnz) then
              do j=iys,iye
                do i=1,nnx
                  stat(k,m1+l)=stat(k,m1+l) +
     +                0.5*(wf(i,j,k) - gm(i,j,k,5))*
     +                (t(i,j,l,k)-txym(k,l)+t(i,j,l,kp1)-txym(kp1,l))
                enddo
              enddo
            else

c ------------------- monotone fluxes

              do j=iys,iye
                do i=1,nnx
                  wflx = wf(i,j,k) - gm(i,j,k,5)
                  stat(k,m1+l) = stat(k,m1+l) +
     +                amax1(sgn*wflx,0.)*(t(i,j,l,k) +
     +                rlim(t(i,j,l,kp1),t(i,j,l,k),t(i,j,l,km1))) +
     +                amin1(sgn*wflx,0.)*(t(i,j,l,kp1) +
     +                rlim(t(i,j,l,k),t(i,j,l,kp1),t(i,j,l,kp2)))
                enddo
              enddo
            endif
            stat(k,m1+l)= sgn*stat(k,m1+l)

            do j=iys,iye
              do i=1,nnx
                stat(k,m6+l) = stat(k,m6+l) +
     +              (w(i,j,k) - wxym(k))*
     +              (t(i,j,l,k)-txym(k,l))
              enddo
            enddo

c ------------ get horizontal scalar resolved fluxes

            do j=iys,iye
              do i=1,nnx
                stat(k,m2+l) = stat(k,m2+l)+
     +              uf(i,j,k)*(t(i,j,l,k)-txym(k,l))
                stat(k,m3+l) = stat(k,m3+l)+
     +              vf(i,j,k)*(t(i,j,l,k)-txym(k,l))
              enddo
            enddo

c ------------------- scalar variances & higher moments

            do j=iys,iye
              do i=1,nnx
                stat(k,m4+l) = stat(k,m4+l) +
     +              (t(i,j,l,k) - txym(k,l))**2
                stat(k,m5+l) = stat(k,m5+l) +
     +              (t(i,j,l,k) - txym(k,l))**3
              enddo
            enddo

c ------ end scalar loop

          enddo

c ------ end z loop

        enddo

c -------- add partial sums and send it to all

        call mpi_sum_z(stat(1,1),i_root,myid,nnz*nstat,1)

c ------ fill arrays for printout and constant file

        do k=1,nnz

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
          RHxym(k)    = stat(k,38)*fnxy

c ------------ get scalar resolved fluxes and variances

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

c ---------- drag at first zu grid point

        drg_x = drg_z_x(1)
        drg_y = drg_z_y(1)

        return
      end
      subroutine tke_budget

c -------- get terms in resolved scale tke budget
c          as in gabls writeup at w-points
c          at istage = 1.
c          t_diss and tr_tau are in comp1, pxym is in get_means

        use pars
        use fields
        use con_data
        use con_stats

        real stat(1:nnz,2)

c -------- stat(.,1) = tke transport  = wq
c          stat(.,2) = pressure transport  = wp

        do k=1,nnz
          stat(k,1) = 0.0
          stat(k,2) = 0.0
        enddo

c --------------- get transport terms as vertical arrays

        do k=izs,ize

          km1 = k - 1
          do j=iys,iye
            do i=1,nnx

c --------- get estimate of turbulent transport term

              ufluc   = u(i,j,k) - uxym(k)
              vfluc   = v(i,j,k) - vxym(k)
              wfluc   = w(i,j,k) - wxym(k)
              stat(k,1)  = stat(k,1) + 0.5*wfluc*
     +            (ufluc**2 + vfluc**2 + wfluc**2)

c --------- get estimate of pressure transport term

              pfluc     = pm(k) + p(i,j,k) - c23*e(i,j,k)
              stat(k,2) = stat(k,2) + pfluc*wfluc
            enddo
          enddo
          stat(k,1) = stat(k,1)*fnxy
          stat(k,2) = stat(k,2)*fnxy
        enddo
        call mpi_sum_z(stat(1,1),i_root,myid,nnz*2,1)

c ------ we have all terms on all processors for all z, add them up
c        treat tr_tau at bottom special, tr_tau = 0

        tr_tau(0) = 0.0
        do k=1,nnz

          kp1 = k + 1
          km1 = k - 1
          if(k == nnz) then
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

c ------------- gather all the budget terms

          t_tran(k)  = t_wq(k) + t_wp(k) + t_tau(k)
          t_rprod(k) = -(dudz*uwle(k) + dvdz*vwle(k))
          t_sprod(k) =  (dudz*uwsb(k) + dvdz*vwsb(k))
          t_buoy(k)  =  batag*wtle(k,1)

        enddo

        return
      end
      subroutine write_his(iloc)

c ----- write history file with global parameters
c       write tsfcc specially to preserve digits!

        use pars
        use fields
        use con_data
        use con_stats
        use particles

        real save_theta(0:nnz)

        divgmax = 0.0
        do iz=1,nnz
          divgmax = amax1(divgmax, divz(iz))
        enddo

        ziavg = zi
        holtop = hol
        wt_min = wtsb(iloc,1)
        wt_le  = wtle(iloc,1)
        krec = krec + 1
        mid = nnz/4

        write(nhis1,6000) time,dt,utau,ziavg,amonin,holtop,
     +      (tsfcc(1)-t_ref),uusfc,uvsfc,uwsfc,vvsfc,vwsfc,wwsfc,
     +      divgmax, wt_min, wt_le,
     +      ucfl, vcfl, wcfl, wtsfc(1),
     +      ups(mid),vps(mid),wps(mid),tps(mid,1),
     +      uwle(mid),uwsb(mid),uw_tot(mid),
     +      vwle(mid),vwsb(mid),vw_tot(mid),
     +      uwcle(mid),vwcle(mid),
     +      wtle(mid,1),wtsb(mid,1),wctle(mid,1),wt_tot(mid,1),
     +      englez(mid),eavg(mid), wabs, float(ip_max),
     +      p_iter, p_err, drg_x, drg_y,
     +      Rep_avg,float(his_reintro)
 6000   format(46e17.8)

c ------ to increase resolution of theta
c        subtract reference before saving history file

        do k=1,nnz
          save_theta(k) = txym(k,1)
          txym(k,1)     = txym(k,1) - t_ref
        enddo

c -------------- write profile information

        call write_prof(nhisp,krec,isize,c_s%wwsb)

        do k=1,nnz
          txym(k,1) = save_theta(k)
        enddo

        return
      end
      subroutine write_prof(nhisp,krec,num,f)
        real f(num)
        real*4 f32(num)

c -------------- build special 32 bit arrays for profiles

        do i=1,num
          f32(i) = f(i)
        enddo

        write(nhisp,err=999,rec=krec) (f32(i),i=1,num)

        return
c --------------- errors
  999   continue
        write(6,9000) num,krec
 9000   format(' 9000, trouble in ',
     +      'SR. save_prof cannot write profile data ',/,
     +      ' num = ',i8, 'krec = ',i6)
        stop
      end
      subroutine close_his

c ---- close history files

        use pars
        logical there

c ---- root closes and checks the files

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

        return
c -------------------- process write errors
 7000   format(' HISTORY DATA IS WRITTEN IN FILE  ',a80)
 7100   format(' PROFILE HISTORY DATA IS WRITTEN IN FILE  ',a80)
 8000   format(' SR. SAVE_HIS: Truble history file not in path =',a80)
 8100   format(' SR. SAVE_HIS: Trouble profile history file',
     +      ' not in path =',a80)
      end
      subroutine dealias_grid_speed

c --------- wave cutoff filter using 2d fft
c           for boundary condition values

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        real wve(nny,jxs:jxe), gm_work(nnx+2,iys:iye)

c --------- sharp spectral cutoff, specific to current 2dfft

        ix_cut   = 2*int(float(nnx)/3.) + 3
        iy_cut_l = int(float(nny)/3.) + 2
        iy_cut_u = nnyp2 - iy_cut_l

        if(iss == 0) then

          do j=iys,iye
            do i=1,nnx
c        gm_work(i,j) = gm(i,j,izs-1,5)
              gm_work(i,j) = gm2(i,j,izs-1,5)
            enddo
          enddo
          call fft2d_mpi(gm_work(1,iys),wve(1,jxs),
     +        trigx(1,1),trigc,
     +        nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +        1,1,myid,ncpu_s,numprocs,-2)
          call sharp_bndy(wve,ix_cut,iy_cut_l,iy_cut_u)
          call fft2d_mpi(gm_work(1,iys),wve(1,jxs),
     +        trigx(1,1),trigc,
     +        nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +        1,1,myid,ncpu_s,numprocs,2)
          do j=iys,iye
            do i=1,nnx
c        gm(i,j,izs-1,5)  = gm_work(i,j)
              gm2(i,j,izs-1,5) = gm_work(i,j)
            enddo
          enddo

        endif

        if(ise == numprocs-1) then

          do j=iys,iye
            do i=1,nnx
c        gm_work(i,j) = gm(i,j,ize,5)
              gm_work(i,j) = gm2(i,j,ize,5)
            enddo
          enddo
          call fft2d_mpi(gm_work(1,iys),wve(1,jxs),
     +        trigx(1,1),trigc,
     +        nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +        1,1,myid,ncpu_s,numprocs,-2)
          call sharp_bndy(wve,ix_cut,iy_cut_l,iy_cut_u)
          call fft2d_mpi(gm_work(1,iys),wve(1,jxs),
     +        trigx(1,1),trigc,
     +        nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +        1,1,myid,ncpu_s,numprocs,2)
          do j=iys,iye
            do i=1,nnx
c        gm(i,j,ize,5)  = gm_work(i,j)
              gm2(i,j,ize,5) = gm_work(i,j)
            enddo
          enddo

        endif

        return
      end
      subroutine dealias

c --------- wave cutoff filter using 2d fft

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        real wve(nny,jxs:jxe,izs:ize)
        real wves(nnxp2,iys:iye,izs:ize)

c --------- sharp spectral cutoff, specific to current 2dfft

        ix_cut   = 2*int(float(nnx)/3.) + 3
        iy_cut_l = int(float(nny)/3.) + 2
        iy_cut_u = nnyp2 - iy_cut_l

c ---------- u-equation

        call fft2d_mpi(u(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +      nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +      izs,ize,myid,ncpu_s,numprocs,-2)
        call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
        call fft2d_mpi(u(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +      nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +      izs,ize,myid,ncpu_s,numprocs,2)

c ---------- v-equation

        call fft2d_mpi(v(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +      nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +      izs,ize,myid,ncpu_s,numprocs,-2)
        call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
        call fft2d_mpi(v(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +      nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +      izs,ize,myid,ncpu_s,numprocs,2)

c ---------- w-equation

        call fft2d_mpi(w(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +      nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +      izs,ize,myid,ncpu_s,numprocs,-2)
        call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
        call fft2d_mpi(w(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +      nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +      izs,ize,myid,ncpu_s,numprocs,2)

c ---------- e-equation

        call fft2d_mpi(e(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +      nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +      izs,ize,myid,ncpu_s,numprocs,-2)
        call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
        call fft2d_mpi(e(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +      nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +      izs,ize,myid,ncpu_s,numprocs,2)

c ---------- uf-equation

        call fft2d_mpi(uf(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +      nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +      izs,ize,myid,ncpu_s,numprocs,-2)
        call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
        call fft2d_mpi(uf(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +      nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +      izs,ize,myid,ncpu_s,numprocs,2)

c ---------- vf-equation

        call fft2d_mpi(vf(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +      nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +      izs,ize,myid,ncpu_s,numprocs,-2)
        call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
        call fft2d_mpi(vf(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +      nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +      izs,ize,myid,ncpu_s,numprocs,2)

c ---------- wf-equation

        call fft2d_mpi(wf(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +      nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +      izs,ize,myid,ncpu_s,numprocs,-2)
        call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
        call fft2d_mpi(wf(1,iys,izs),wve(1,jxs,izs),trigx(1,1),trigc,
     +      nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +      izs,ize,myid,ncpu_s,numprocs,2)

c ------------- scalars, not stored in correct order

        do iscl=1,nscl
          do iz=izs,ize
            do iy=iys,iye
              do ix=1,nnx
                wves(ix,iy,iz) = t(ix,iy,iscl,iz)
              enddo
            enddo
          enddo
          call fft2d_mpi(wves(1,iys,izs),wve(1,jxs,izs),trigx(1,1),
     +        trigc,nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +        izs,ize,myid,ncpu_s,numprocs,-2)
          call sharp(wve,ix_cut,iy_cut_l,iy_cut_u)
          call fft2d_mpi(wves(1,iys,izs),wve(1,jxs,izs),trigx(1,1),
     +        trigc,nnx,nny,jxs,jxe,jx_s,jx_e,iys,iye,iy_s,iy_e,
     +        izs,ize,myid,ncpu_s,numprocs,2)
          do iz=izs,ize
            do iy=iys,iye
              do ix=1,nnx
                t(ix,iy,iscl,iz) = wves(ix,iy,iz)
              enddo
            enddo
          enddo
        enddo

        return
      end
      subroutine sharp(wve,ix_cut,iy_cut_l,iy_cut_u)

c --------- sharp cutoff filter for field wve stored
c           in 2d-fft order

        use pars
        real wve(nny,jxs:jxe,izs:ize)

        do iz=izs,ize
          do ix=jxs,jxe
            do iy=iy_cut_l,iy_cut_u
              wve(iy,ix,iz) = 0.0
            enddo
          enddo
        enddo

        if(jxe < ix_cut) go to 999

        do iz=izs,ize
          do ix=max(jxs,ix_cut),jxe
            do iy=1,nny
              wve(iy,ix,iz) = 0.0
            enddo
          enddo
        enddo

  999   continue

        return
      end
      subroutine sharp_bndy(wve,ix_cut,iy_cut_l,iy_cut_u)

c --------- sharp cutoff filter for boundary array
c           in 2d-fft order

        use pars
        real wve(nny,jxs:jxe)

        do ix=jxs,jxe
          do iy=iy_cut_l,iy_cut_u
            wve(iy,ix) = 0.0
          enddo
        enddo

        if(jxe < ix_cut) go to 999

        do ix=max(jxs,ix_cut),jxe
          do iy=1,nny
            wve(iy,ix) = 0.0
          enddo
        enddo

  999   continue

        return
      end
      subroutine gridd

c ----------- allocate space and pass arrays using modules

        use pars
        use fields
        use fftwk
        use con_data
        use con_stats
        use particles

c ------------- establish association between pointers
c               and data structures

        call fill_cc
        call fill_cs

        if (myid==0) write(6,5001) isize
 5001   format(' size of stats array = ',i8)

c ---------------- debug for arrays

        big = -99.0e+300

c ---------------- setup grid

        nnx = nxg1
        nny = nyg1
        nnz = nzg1

c ----------- make sure problem and cpu's match

        maxp   = numprocs-1
        ncpu_z = numprocs/ncpu_s
        if(mod(numprocs,ncpu_s) /= 0 .or.
     +      ncpu_z > nnz) then
          go to 999
        endif
        if(l_root) write(6, 1100) ncpu_s, ncpu_z, numprocs, maxp
        write(nprt,1100) ncpu_s, ncpu_z, numprocs, maxp
 1100   format(' Number of x-y slab cpus = ',i5,/,
     +      ' Number of z-level cpus  = ',i5,/,
     +      ' Total number of cpus    = ',i5,/,
     +      ' Max-p for index arrays  = ',i5)

c ---------------- allocate arrays for (i,j,k)-indexing on
c                  each processor (see set_range)

        allocate(ix_s(0:maxp), ix_e(0:maxp),
     +      jx_s(0:maxp), jx_e(0:maxp),
     +      kx_s(0:maxp), kx_e(0:maxp),
     +      mx_s(0:maxp), mx_e(0:maxp),
     +      iy_s(0:maxp), iy_e(0:maxp),
     +      jy_s(0:maxp), jy_e(0:maxp),
     +      is_s(0:maxp), is_e(0:maxp),
     +      iz_s(0:maxp), iz_e(0:maxp),
     +      i2x_s(0:maxp), i2x_e(0:maxp))

c ---------------- setup array sizes and variable dimensions

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

        write(nprt,7001) nnx,nny,nnz
 7001   format(' 7001 gridd nnx = ',i4,' nny = ',i4,' nnz = ',i4)

        call set_range

        num_y = iye + 1 - iys

c ------------- allocate solution arrays
c               account for nnxp2 for fields but not in rhs
c               and possible monotone for scalars

        allocate(u(nnxp2,iys:iye,izs-1:ize+1),
     +      v(nnxp2,iys:iye,izs-1:ize+1),
     +      w(nnxp2,iys:iye,izs-1:ize+1),
     +      t(nnxp2,iys:iye,nscl,izs-2:ize+2),
     +      e(nnxp2,iys:iye,izs-1:ize+1),
     +      r1(nnx,iys:iye,izs-1:ize+1),
     +      r2(nnx,iys:iye,izs-1:ize+1),
     +      r3(nnx,iys:iye,izs-1:ize+1),
     +      r4(nnx,iys:iye,nscl,izs-1:ize+1),
     +      r5(nnx,iys:iye,izs-1:ize+1),
     +      rj(nnx,iys:iye,izs-1:ize+1),
     +      uf(nnxp2,iys:iye,izs-1:ize+1),
     +      vf(nnxp2,iys:iye,izs-1:ize+1),
     +      wf(nnxp2,iys:iye,izs-1:ize+1),
     +      ufi(nnxp2,iys:iye,izs-1:ize+1),
     +      vfi(nnxp2,iys:iye,izs-1:ize+1),
     +      wfi(nnxp2,iys:iye,izs-1:ize+1))

c ------------- allocate space for grid and metrics

        allocate(xgrid(nnx),
     +      ygrid(nny),
     +      zgrid_w(nnx,iys:iye,izs-1:ize+1,3),
     +      zgrid_u(nnx,iys:iye,izs-1:ize+1),
     +      gm(nnx,iys:iye,izs-1:ize+1,5),
     +      gm2(nnx,iys:iye,izs-1:ize+1,5),
     +      gm3(nnx,iys:iye,izs-1:ize+1,5),
     +      gms(nnx,iys:iye,3),
     +      bndy(nnxp2,iys:iye),
     +      hdot(nnxp2,iys:iye),
     +      u_orbit(nnxp2,iys:iye),
     +      v_orbit(nnxp2,iys:iye),
     +      w_orbit(nnxp2,iys:iye),
     +      bndy_t(nnx,iys:iye,3),
     +      hdot_t(nnx,iys:iye,3),
     +      z1b(nnx,iys:iye),
     +      z1b_2(nnx,iys:iye),
     +      z1b_3(nnx,iys:iye))

c ------------- allocate extended arrays for interpolation of
c               particle/spray location

        if (ispray > 0) then
          allocate(uext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),
     +        vext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),
     +        wext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),
     +        Text(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3),
     +        T2ext(0:nnz+1,iys-2:iye+3,mxs-2:mxe+3))
      !Transposed velocities to do the uf interpolation:
          allocate(u_t(0:nnz+1,iys:iye,mxs:mxe),
     +        v_t(0:nnz+1,iys:iye,mxs:mxe),
     +        w_t(0:nnz+1,iys:iye,mxs:mxe),
     +        T_t(0:nnz+1,iys:iye,mxs:mxe),
     +        T2_t(0:nnz+1,iys:iye,mxs:mxe))
        endif
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

c ------------- allocate space for boundary condition arrays
c               on top and bottom of domain

        allocate(ubc(nnx,iys:iye,2),
     +      vbc(nnx,iys:iye,2),
     +      wbc(nnx,iys:iye,2),
     +      tbc(nnx,iys:iye,nscl,2),
     +      ebc(nnx,iys:iye,2),
     +      pbc(nnx,iys:iye,2),
     +      pbc2(nnx,iys:iye,2),
     +      ufbc(nnx,iys:iye,2),
     +      vfbc(nnx,iys:iye,2),
     +      wfbc(nnx,iys:iye,4))

c ------------ allocate space for wind and surface arrays

        allocate(wind(nnx,iys:iye),
     +      tau13m(nnx,iys:iye),
     +      tau23m(nnx,iys:iye),
     +      taut3m(nnx,iys:iye,nscl),
     +      t_grnd(nnx,iys:iye,nscl))

c ------------------- allocate space for derivative arrays

        allocate(ux(nnx,iys:iye,izs-1:ize+1),
     +      uy(nnx,iys:iye,izs-1:ize+1),
     +      vx(nnx,iys:iye,izs-1:ize+1),
     +      vy(nnx,iys:iye,izs-1:ize+1),
     +      wx(nnx,iys:iye,izs-1:ize+1),
     +      wy(nnx,iys:iye,izs-1:ize+1))

c ------------- allocate space for pressure, pressure bcs

        allocate(p(nnxp2,iys:iye,izs-1:ize+1),
     +      p_old(nnxp2,iys:iye,izs-1:ize+1),
     +      ptop(nnxp2,iys:iye,2))

c ------------- allocate space for viscosity and diffusivity

        allocate(vis_m(nnx,iys:iye,izs-1:ize+1),
     +      vis_s(nnx,iys:iye,izs-1:ize+1))

c ------------- allocate space for tau bc

        allocate(t11bc(nnx,iys:iye),
     +      t12bc(nnx,iys:iye),
     +      t13bc(nnx,iys:iye),
     +      t22bc(nnx,iys:iye),
     +      t23bc(nnx,iys:iye),
     +      t33bc(nnx,iys:iye),
     +      ustar(nnx,iys:iye),
     +      amol(nnx,iys:iye),
     +      t1tbc(nnx,iys:iye,nscl),
     +      t2tbc(nnx,iys:iye,nscl),
     +      t3tbc(nnx,iys:iye,nscl))

c ------------- allocate space for fft trig factors

        nq_trig = max(nnx,nny)
        allocate(trigx(2*nq_trig+15,2),
     +      trigc(4*nq_trig+15), trigcx(4*nq_trig+15))
        return
  999   continue

        if(l_root) write(6,1000) numprocs, ncpu_s, mmz
        write(nprt,1000) numprocs, ncpu_s, nnz
 1000   format(' Gridd Trouble number of processors and grid',
     +      ' partitioning do not match!',/,
     +      ' Total num of cpus   = ',i5,
     +      ' Num cpu on x-y slab = ',i5,/,
     +      ' Num of z-levels     = ',i5)
        call mpi_finalize(ierr)
      end
      subroutine restart

c ----------- get restart file from local directory

        use pars
        use fields
        use con_data
        use con_stats
        character*80 path_res_c
        logical there

c --------------------- check if file is there

        inquire(file=path_res,exist=there)
        if(there) then
          if(l_root) write(6,6001) path_res
        else
          if(l_root) write(6,6005) path_res
          stop
        endif

c ------------------ get constant file

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
     +      status='old')

        call read_res

        return
c ---------------------------- process errors
  100   continue
        write(6,9000) path_res, nvel
        call mpi_finalize(ierr)
        stop
c -----------------------
  200   continue
        write(6,9001) path_res_c, nvelc
        call mpi_finalize(ierr)
        stop
c -----------------------
 6001   format(' SR. RESTART: FILE READ = ',A80)
 6002   format(' SR. RESTART: CONSTANT FILE READ = ',A80)
 6005   format(' 6005, SR. RESTART: cannot find restart file = ',a80)
 6006   format(' 6005, SR. RESTART: cannot find constant file = ',a80)
 9000   format(' 9000, SR. RESTART: cannot open file =',a80,/,
     +      ' to unit number = ',i2)
 9001   format(' 9001, SR. RESTART: cannot open file =',a80,/,
     +      ' to unit number = ',i2)
      end
      subroutine recv_yz_var(temp_x,nvar,nny,iys,iye,izs,ize,ir)

        include 'mpif.h'
        integer istatus(mpi_status_size)

        real buf(nvar,iys:iye,izs:ize)
        real(kind=4), dimension(nvar,nny,izs:ize) :: temp_x

        num = nvar*(ize+1-izs)*(iye+1-iys)
        call mpi_recv(buf(1,iys,izs),num,mpi_real8,ir,1,
     +      mpi_comm_world,istatus,ierr)
        do k=izs,ize
          do j=iys,iye
            do ii=1,nvar
              temp_x(ii,j,k) = buf(ii,j,k)
            enddo
          enddo
        enddo

        return
      end
      subroutine range(n1,n2,nprocs,irank,ista,iend)

c ---------- the ibm range finder to balance load

        iwork1 = (n2 - n1 + 1)/nprocs
        iwork2 = mod(n2 - n1 +1, nprocs)
        ista = irank*iwork1 + n1 + min(irank,iwork2)
        iend = ista + iwork1 - 1
        if(iwork2 > irank) iend = iend + 1

        return
      end
      subroutine set_range

c ---- build special x,y,z-ranges. dimensioned for 0:numprocs-1
c      indexed with myid

c      [ix_s:ix_e] x-range for computing y-derivatives nx-pts/ncpu_s
c                  in xtoy and ytox tranposes

c      [i2x_s:2ix_e] x-range for computing complx y-derivatives 2*nx-pts/ncpu_s
c                  in xtoy and ytox tranposes. special for moving surfaces

c      [jx_s:jx_e] x-range for computing 2d fft (nx+2)-pts/ncpu_s
c                  must be even in each x-interval for complex fft in y

c      [kx_s:kx_e] x-range for pressure solver transpose (nx+2)-pts/ncpu_z
c                  nx+2 fourier coefficients for xtoz and ztox transposes

c      [mx_s:mx_e] x-range split across z cpus as nx-pts/ncpu_z
c                  for use in surface layer routines

c      [is_s:is_e] starting and ending processor id's for a
c                  particular z-level

c      [iy_s:iy_e] y-range for computing y-derivatives ny-pts/ncpu_s
c                  in xtoy and ytox tranposes

c      [jy_s:jy_e] y-range for use in xtoz and ztox transposes
c                  in pressure solution

c      [iz_s:iz_e] z-range for a particular vertical slab

        use pars

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

            iy_s(ii)  = ny_s
            iy_e(ii)  = ny_e
            jy_s(ii)  = ly_s
            jy_e(ii)  = ly_e

            iz_s(ii)  = mz_s
            iz_e(ii)  = mz_e

            is_s(ii)  = (ii/ncpu_s)*ncpu_s
            is_e(ii)  = is_s(ii) + ncpu_s - 1
          enddo
        enddo

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

c ----------- get starting and  ending processor id's on each
c             vertical slab

        iss = is_s(myid)
        ise = is_e(myid)

c ------------ debug ranges

        if(l_debug) then
          write(nprt,1200) myid, (nn, ix_s(nn), ix_e(nn), jx_s(nn),
     +        jx_e(nn), kx_s(nn), kx_e(nn),
     +        i2x_s(nn), i2x_e(nn),
     +        nn = 0,numprocs-1)
 1200     format(' myid =  ',i4,/,
     +        ' nn',5x,' ixs ',5x,' ixe ',5x,' jxs ',5x,' jxe '
     +        ,5x,' kxs ',5x,' kxe',5x,' i2xs ',5x,' i2xe',/,(9i6))

          write(nprt,1213) myid, (nn, iy_s(nn), iy_e(nn),
     +        jy_s(nn), jy_e(nn),
     +        iz_s(nn), iz_e(nn), is_s(nn), is_e(nn),
     +        nn=0,numprocs-1)
 1213     format(' myid = ',i3,/,
     +        ' nn ',3x,' iys ',5x,' iye ',5x,
     +        ' jys ',5x,' jye ',5x,
     +        ' izs ',5x,' ize',5x,' iss ',5x,' ise ',/,
     +        (9i6))
        endif

        return
      end
      subroutine mpi_sum_xy(f,myid,iss,ise,nsend)

c --------- get horizontal x-y sum over a set of proccessors [iss:ise]
c           for vector f(i). f(i) is overwritten. skip if single processor

        include 'mpif.h'
        integer istatus(mpi_status_size)
        real work(nsend,iss:ise), f(nsend)

        if(iss == ise) go to 999

        do j=1,nsend
          work(j,myid) = f(j)
          f(j)         = 0.0
        enddo
        do i=iss,ise
          if(i /= myid) then
            call mpi_sendrecv(work(1,myid),nsend,mpi_real8,i,1,
     +          work(1,i),nsend,mpi_real8,i,1,
     +          mpi_comm_world,istatus,ierr)
          endif
        enddo
        do i=iss,ise
          do j=1,nsend
            f(j) = f(j) + work(j,i)
          enddo
        enddo

  999   continue

        return
      end
      subroutine mpi_sum_z(f,i_root,myid,nsend,iall)

c --------- get sums on root or all processors
c           for all z for vector f(i)

        include 'mpif.h'
        integer istatus(mpi_status_size)
        real recv_b(nsend), f(nsend)

c -------- just root gets the result

        if(iall /= 1) then
          call mpi_reduce(f(1),recv_b(1),nsend,mpi_real8,mpi_sum,i_root,
     +        mpi_comm_world,ierr)
          if(myid == i_root) then
            do i=1,nsend
              f(i) = recv_b(i)
            enddo
          endif
        else

c -------- everyone gets the result

          call mpi_allreduce(f(1),recv_b(1),nsend,mpi_real8,mpi_sum,
     +        mpi_comm_world,ierr)
          do i=1,nsend
            f(i) = recv_b(i)
          enddo
        endif

        return
      end
      subroutine mpi_sum_z_s(f,i_root,myid,nsend,nscl,iall)

c --------- get sums on root or all processors
c           for all z for vector f(i,nscl)

        include 'mpif.h'
        integer istatus(mpi_status_size)
        real recv_b(nsend,nscl), f(nsend,nscl)

        if(iall /= 1) then
          call mpi_reduce(f(1,1),recv_b(1,1),nsend*nscl,mpi_real8,
     +        mpi_sum,i_root,mpi_comm_world,ierr)
          if(myid == i_root) then
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

        return
      end
      subroutine mpi_gath_root(fs,fr,iz_s,iz_e,izs,ize,nz,myid,np,ns)

c ---------- gather results on root processors

        include 'mpif.h'
        integer istatus(mpi_status_size)
        integer iz_s(0:np-1), iz_e(0:np-1)
        real fs(izs:ize), fr(1:nz)

        if(np == 1) go to 999

        irow_r = mod(myid,ns)
        if(myid > ns) then
          call mpi_send(fs(izs),ize+1-izs,mpi_real8,irow_r,1,
     +        mpi_comm_world,ierr)
        else
          do l=irow_r+ns,np-1,ns
            ind = iz_s(l) + 1
            num = iz_e(l) + 1 - iz_s(l)
            call mpi_recv(fr(ind),num,mpi_real8,l,1,
     +          mpi_comm_world,istatus,ierr)
          enddo
        endif

  999   continue

        return
      end
      subroutine mpi_send_root(fs,num,myid,np,ns)

c ---------- send root results to other processors above it

        include 'mpif.h'
        integer istatus(mpi_status_size)
        real fs(num)

        if(np == 1) go to 999

        irow_r = mod(myid,ns)
        if(myid >= ns) then
          call mpi_recv(fs(1),num,mpi_real8,irow_r,1,
     +        mpi_comm_world,istatus,ierr)
        else
          do l=irow_r+ns,np-1,ns
            call mpi_send(fs(1),num,mpi_real8,l,1,
     +          mpi_comm_world,ierr)
          enddo
        endif

  999   continue

        return
      end
      subroutine xtoy_trans(f,g,nx,ny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,iz1,iz2,
     +      myid,ncpu_s,np)

c ------- transpose array  f(nx,iys:iye,iz1:iz2) ---> g(ny,ixs:ixe,iz1:iz2)

        include 'mpif.h'
        integer istatus(mpi_status_size)

        real f(nx,iys:iye,iz1:iz2),
     +      g(ny,ixs:ixe,iz1:iz2)
        real ft(nx*(iye+1-iys)*(iz2 - iz1 + 1)),
     +      gt(ny*(ixe+1-ixs)*(iz2 - iz1 + 1))
        integer ix_s(0:np-1), ix_e(0:np-1),
     +      iy_s(0:np-1), iy_e(0:np-1)

        jk = (iye - iys + 1)*(iz2 - iz1 + 1)
        ik = (ixe - ixs + 1)*(iz2 - iz1 + 1)

c ----------- get cpus on slab for myid

        islab = myid/ncpu_s
        iss   = islab*ncpu_s
        ise   = iss + ncpu_s - 1

        do i=iss,ise
          nsend = (ix_e(i) - ix_s(i) + 1)*jk
          nrecv = (iy_e(i) - iy_s(i) + 1)*ik
          if(i == myid) then
            call send_xtoy(f,gt(1),nx,ix_s(i),ix_e(i),
     +          iy_s(myid),iy_e(myid),iz1,iz2)
          else
            call send_xtoy(f,ft(1),nx,ix_s(i),ix_e(i),
     +          iy_s(myid),iy_e(myid),iz1,iz2)
            call mpi_sendrecv(
     +          ft(1),nsend,mpi_real8,i,1,
     +          gt(1),nrecv,mpi_real8,i,1,
     +          mpi_comm_world,istatus,ierr)
          endif
          call recv_xtoy(g,gt(1),ny,ix_s(myid),ix_e(myid),
     +        iy_s(i),iy_e(i),iz1,iz2)
        enddo

        return
      end
      subroutine send_xtoy(f,ft,nx,ixs,ixe,iys,iye,izs,ize)

c ------------- grab correct chunk of array to be sent

        real f(nx,iys:iye,izs:ize), ft(ixs:ixe,iys:iye,izs:ize)

        do k=izs,ize
          do j=iys,iye
            do i=ixs,ixe
              ft(i,j,k) = f(i,j,k)
            enddo
          enddo
        enddo

        return
      end
      subroutine recv_xtoy(g,gt,ny,ixs,ixe,iys,iye,izs,ize)
        real g(ny,ixs:ixe,izs:ize), gt(ixs:ixe,iys:iye,izs:ize)

        do k=izs,ize
          do j=iys,iye
            do i=ixs,ixe
              g(j,i,k) = gt(i,j,k)
            enddo
          enddo
        enddo

        return
      end
      subroutine ytox_trans(g,f,nx,ny,ixs,ixe,ix_s,ix_e,
     +      iys,iye,iy_s,iy_e,iz1,iz2,
     +      myid,ncpu_s,np)

c ------- transpose array g(ny,ixs:ixe,iz1:iz2) ---> f(nx,iys:iye,iz1:iz2)

        include 'mpif.h'
        integer istatus(mpi_status_size)

        real f(nx,iys:iye,iz1:iz2),
     +      g(ny,ixs:ixe,iz1:iz2)
        real ft(nx*(iye+1-iys)*(iz2 - iz1 + 1)),
     +      gt(ny*(ixe+1-ixs)*(iz2 - iz1 + 1))

        integer ix_s(0:np-1), ix_e(0:np-1),
     +      iy_s(0:np-1), iy_e(0:np-1)

        jk = (iye - iys + 1)*(iz2 - iz1 + 1)
        ik = (ixe - ixs + 1)*(iz2 - iz1 + 1)

c ----------- get cpus on slab for myid

        islab = myid/ncpu_s
        iss   = islab*ncpu_s
        ise   = iss + ncpu_s - 1
        do i=iss,ise
          nsend = (iy_e(i) - iy_s(i) + 1)*ik
          nrecv = (ix_e(i) - ix_s(i) + 1)*jk
          if(i == myid) then
            call send_ytox(g,ft(1),ny,ix_s(myid),ix_e(myid),
     +          iy_s(i),iy_e(i),iz1,iz2)
          else
            call send_ytox(g,gt(1),ny,ix_s(myid),ix_e(myid),
     +          iy_s(i),iy_e(i),iz1,iz2)
            call mpi_sendrecv(
     +          gt(1),nsend,mpi_real8,i,1,
     +          ft(1),nrecv,mpi_real8,i,1,
     +          mpi_comm_world,istatus,ierr)
          endif
          call recv_ytox(f,ft(1),nx,ix_s(i),ix_e(i),
     +        iy_s(myid),iy_e(myid),iz1,iz2)
        enddo

        return
      end
      subroutine send_ytox(g,gt,ny,ixs,ixe,iys,iye,izs,ize)

c ------------- grab correct chunk of array to be sent

        real g(ny,ixs:ixe,izs:ize), gt(iys:iye,ixs:ixe,izs:ize)

        do k=izs,ize
          do i=ixs,ixe
            do j=iys,iye
              gt(j,i,k) = g(j,i,k)
            enddo
          enddo
        enddo

        return
      end
      subroutine recv_ytox(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
        real f(nx,iys:iye,izs:ize), ft(iys:iye,ixs:ixe,izs:ize)

        do k=izs,ize
          do i=ixs,ixe
            do j=iys,iye
              f(i,j,k) = ft(j,i,k)
            enddo
          enddo
        enddo

        return
      end
      subroutine xtoz_trans(f,g,nx,nz,ixs,ixe,ix_s,ix_e,
     +      iys,iye,izs,ize,iz_s,iz_e,
     +      myid,ncpu_s,numprocs)

c ------- transpose array  f(nx,iys:iye,izs-1:ize+1)
c                     ---> g(0:nz+1,iys:iye,ixs:ixe)

        include 'mpif.h'
        integer istatus(mpi_status_size)

        real f(nx,iys:iye,izs-1:ize+1), g(0:nz+1,iys:iye,ixs:ixe)
        real ft(nx*(iye+1-iys)*(ize-izs+1)),
     +      gt(nz*(ixe+1-ixs)*(iye-iys+1))
        integer ix_s(0:numprocs-1), ix_e(0:numprocs-1),
     +      iz_s(0:numprocs-1), iz_e(0:numprocs-1)

        jk = (ize - izs + 1)*(iye - iys + 1)
        ij = (ixe - ixs + 1)*(iye - iys + 1)

c ----------- get starting location

        iss = myid - (myid/ncpu_s)*ncpu_s

        do i=iss,numprocs-1,ncpu_s
          nsend = (ix_e(i) - ix_s(i) + 1)*jk
          nrecv = (iz_e(i) - iz_s(i) + 1)*ij
          if(i == myid) then
            call send_xtoz(f,gt(1),nx,ix_s(i),ix_e(i),
     +          iys,iye,iz_s(myid),iz_e(myid))
          else
            call send_xtoz(f,ft(1),nx,ix_s(i),ix_e(i),
     +          iys,iye,iz_s(myid),iz_e(myid))
            call mpi_sendrecv(
     +          ft(1),nsend,mpi_real8,i,1,
     +          gt(1),nrecv,mpi_real8,i,1,
     +          mpi_comm_world,istatus,ierr)
          endif
          call recv_xtoz(g,gt(1),nz,ix_s(myid),ix_e(myid),
     +        iys,iye,iz_s(i),iz_e(i))
        enddo

        return
      end
      subroutine send_xtoz(f,ft,nx,ixs,ixe,iys,iye,izs,ize)

c ------- grab correct chunk of array to be sent and skip ghost points

        real f(nx,iys:iye,izs-1:ize+1), ft(ixs:ixe,iys:iye,izs:ize)

        do k=izs,ize
          do j=iys,iye
            do i=ixs,ixe
              ft(i,j,k) = f(i,j,k)
            enddo
          enddo
        enddo

        return
      end
      subroutine recv_xtoz(g,gt,nz,ixs,ixe,iys,iye,izs,ize)
        real g(0:nz+1,iys:iye,ixs:ixe), gt(ixs:ixe,iys:iye,izs:ize)

        do k=izs,ize
          do j=iys,iye
            do i=ixs,ixe
              g(k,j,i) = gt(i,j,k)
            enddo
          enddo
        enddo

        return
      end
      subroutine ztox_trans(g,f,nx,nz,ixs,ixe,ix_s,ix_e,
     +      iys,iye,izs,ize,iz_s,iz_e,
     +      myid,ncpu_s,numprocs)

c ------- transpose array g(0:nz+1,iys:iye,ixs:ixe)
c                    ---> f(nx,iys:iye,izs-1:ize+1)

        include 'mpif.h'
        integer istatus(mpi_status_size)

        real f(nx,iys:iye,izs-1:ize+1), g(0:nz+1,iys:iye,ixs:ixe)
        real ft(nx*(iye+1-iys)*(ize-izs+3)),
     +      gt((nz+3)*(iye+1-iys)*(ixe-ixs+1))

        integer ix_s(0:numprocs-1), ix_e(0:numprocs-1),
     +      iz_s(0:numprocs-1), iz_e(0:numprocs-1)

        jk = (ize - izs + 3)*(iye - iys + 1)
        ij = (ixe - ixs + 1)*(iye - iys + 1)

c ------------- get starting location

        iss = myid - (myid/ncpu_s)*ncpu_s

        do i=iss,numprocs-1,ncpu_s
          nsend = (iz_e(i) - iz_s(i) + 3)*ij
          nrecv = (ix_e(i) - ix_s(i) + 1)*jk
          if(i == myid) then
            call send_ztox(g,ft(1),nz,ix_s(myid),ix_e(myid),
     +          iys,iye,iz_s(i),iz_e(i))
          else
            call send_ztox(g,gt(1),nz,ix_s(myid),ix_e(myid),
     +          iys,iye,iz_s(i),iz_e(i))
            call mpi_sendrecv(
     +          gt(1),nsend,mpi_real8,i,1,
     +          ft(1),nrecv,mpi_real8,i,1,
     +          mpi_comm_world,istatus,ierr)
          endif
          call recv_ztox(f,ft(1),nx,ix_s(i),ix_e(i),
     +        iys,iye,iz_s(myid),iz_e(myid))
        enddo

        return
      end
      subroutine send_ztox(g,gt,nz,ixs,ixe,iys,iye,izs,ize)

c ------------- grab correct chunk of array to be sent,
c               account for ghost points

        real g(0:nz+1,iys:iye,ixs:ixe), gt(izs-1:ize+1,iys:iye,ixs:ixe)

        do j=iys,iye
          do i=ixs,ixe
            do k=izs-1,ize+1
              gt(k,j,i) = g(k,j,i)
            enddo
          enddo
        enddo

        return
      end
      subroutine recv_ztox(f,ft,nx,ixs,ixe,iys,iye,izs,ize)
        real f(nx,iys:iye,izs-1:ize+1), ft(izs-1:ize+1,iys:iye,ixs:ixe)

        do i=ixs,ixe
          do j=iys,iye
            do k=izs-1,ize+1
              f(i,j,k) = ft(k,j,i)
            enddo
          enddo
        enddo

        return
      end
      subroutine exchange

c ------------- exchange ghost points with mpi,
c               nb and nt are the destination and
c               source nodes. Allows for 1z per cpu

        use pars
        use fields
c     use fftwk
        include 'mpif.h'

        real fs(nnx,iys:iye,(5+nscl)),fr(nnx,iys:iye,(5+nscl))
        integer istatus(mpi_status_size)

        nb = myid - ncpu_s
        nt = myid + ncpu_s

c ------------ account for endpoints

        if(iss == 0) then
          nb = mpi_proc_null
        endif
        if(ise == numprocs-1) then
          nt = mpi_proc_null
        endif
        nsend = nnx*(iye + 1 - iys)*(5+nscl)
        nrecv = nsend

c --------- send top of myid, receive bottom from myid - ncpu_s

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

        call mpi_sendrecv(
     +      fs(1,iys,1),nsend,mpi_real8,nt,0,
     +      fr(1,iys,1),nrecv,mpi_real8,nb,0,
     +      mpi_comm_world,istatus,ierr)

        if(iss /= 0) then
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

c -------- send bottom of myid, receive bottom from myid + ncpu_s

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

        call mpi_sendrecv(
     +      fs(1,iys,1),nsend,mpi_real8,nb,1,
     +      fr(1,iys,1),nrecv,mpi_real8,nt,1,
     +      mpi_comm_world,istatus,ierr)

        if(ise /= numprocs-1) then
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

c --------------- send extra scalar points

        nsend = nnx*(iye + 1 - iys)*nscl
        nrecv = nsend

c -------------- send top of myid, receive bottom from myid - ncpu_s

        izm1 = ize-1
        do iscl=1,nscl
          do iy=iys,iye
            do ix=1,nnx
              fs(ix,iy,iscl) = t(ix,iy,iscl,izm1)
            enddo
          enddo
        enddo

        call mpi_sendrecv(
     +      fs(1,iys,1),nsend,mpi_real8,nt,0,
     +      fr(1,iys,1),nrecv,mpi_real8,nb,0,
     +      mpi_comm_world,istatus,ierr)

        if(iss /= 0) then
          izm2 = izs-2
          do iscl=1,nscl
            do iy=iys,iye
              do ix=1,nnx
                t(ix,iy,iscl,izm2) = fr(ix,iy,iscl)
              enddo
            enddo
          enddo
        endif

c -------------- send bottom of myid, receive bottom from myid + ncpu_s

        izp1 = izs+1
        do iscl=1,nscl
          do iy=iys,iye
            do ix=1,nnx
              fs(ix,iy,iscl) = t(ix,iy,iscl,izp1)
            enddo
          enddo
        enddo

        call mpi_sendrecv(
     +      fs(1,iys,1),nsend,mpi_real8,nb,1,
     +      fr(1,iys,1),nrecv,mpi_real8,nt,1,
     +      mpi_comm_world,istatus,ierr)

        if(ise /= numprocs-1) then
          izp2 = ize+2
          do iscl=1,nscl
            do iy=iys,iye
              do ix=1,nnx
                t(ix,iy,iscl,izp2) = fr(ix,iy,iscl)
              enddo
            enddo
          enddo
        endif

        return
      end
      subroutine bcast_pbc

c ---- send upper boundary conditions to other processors
c      for fft solution of pressure

        use pars
        use fields
        include 'mpif.h'
        integer istatus(mpi_status_size)

        if(numprocs == 1) go to 999

        irow_r = mod(myid,ncpu_s)
        irow_t = is_s(numprocs-1) + irow_r
        num = nnx*(iye+1-iys)

c ----- check which row myid is in

        if(iss /= is_s(numprocs-1)) then

c ------ not in the top row, receive from top

          call mpi_recv(pbc(1,iys,1),num,mpi_real8,irow_t,1,
     +        mpi_comm_world,istatus,ierr)
        else

c ------ myid is in the top row, send to everyone below

          do l=irow_r,irow_t-ncpu_s,ncpu_s
            call mpi_send(pbc(1,iys,1),num,mpi_real8,l,1,
     +          mpi_comm_world,ierr)
          enddo
        endif

c --------- same thing for another variable

        if(iss /= is_s(numprocs-1)) then

c ------ not in the top row, receive from top

          call mpi_recv(pbc2(1,iys,1),num,mpi_real8,irow_t,1,
     +        mpi_comm_world,istatus,ierr)
        else

c ------ in the top row, send to everyone below

          do l=irow_r,irow_t-ncpu_s,ncpu_s
            call mpi_send(pbc2(1,iys,1),num,mpi_real8,l,1,
     +          mpi_comm_world,ierr)
          enddo
        endif

  999   continue

        return
      end
      subroutine bcast_surf

c ---- send bottom geometry to all other processors
c      for computation of pressure drag at all z

        use pars
        use fields
        include 'mpif.h'
        integer istatus(mpi_status_size)

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

        if(numprocs == 1) go to 999

        irow_r = mod(myid,ncpu_s)
        irow_t = is_s(numprocs-1) + irow_r
        num    = nnx*(iye+1-iys)

c ----- check which row myid is in

        if(iss /= is_s(0)) then

c ------ not in the bottom row, receive from bottom

          call mpi_recv(gms(1,iys,1),num,mpi_real8,irow_r,1,
     +        mpi_comm_world,istatus,ierr)
        else

c ------ myid is in the bottom row, send to everyone else

          do l=irow_r+ncpu_s,irow_t,ncpu_s
            call mpi_send(gms(1,iys,1),num,mpi_real8,l,1,
     +          mpi_comm_world,ierr)
          enddo
        endif

c --------- same thing for another variable

        if(iss /= is_s(0)) then

c ------ not in the bottom row, receive from bottom

          call mpi_recv(gms(1,iys,2),num,mpi_real8,irow_r,1,
     +        mpi_comm_world,istatus,ierr)
        else

c ------ myid is in the bottom row, send to everyone else

          do l=irow_r+ncpu_s,irow_t,ncpu_s
            call mpi_send(gms(1,iys,2),num,mpi_real8,l,1,
     +          mpi_comm_world,ierr)
          enddo
        endif

c --------- same thing for another variable

        if(iss /= is_s(0)) then

c ------ not in the bottom row, receive from bottom

          call mpi_recv(gms(1,iys,3),num,mpi_real8,irow_r,1,
     +        mpi_comm_world,istatus,ierr)
        else

c ------ myid is in the bottom row, send to everyone else

          do l=irow_r+ncpu_s,irow_t,ncpu_s
            call mpi_send(gms(1,iys,3),num,mpi_real8,l,1,
     +          mpi_comm_world,ierr)
          enddo
        endif
  999   continue

        return
      end

      subroutine change_RH_bcs_to_q
      !This converts the RH provided in the params file to a value of q that code needs
        use pars
        use con_data
        use particles
        implicit none
        real :: mod_magnus

      !Assuming tsfcc(2),Ttop(2),Tbot(2) in params.in are giving RH:
        if (iDNS == 1) then

          RHT = Ttop(2)
          RHB = Tbot(2)

         !Convert RH given in input file into specific humidity
          Ttop(2) = RHT/100.0*Mw/Ru/Ttop(1)*mod_magnus(Ttop(1))/rhoa
          Tbot(2) = RHB/100.0*Mw/Ru/Tbot(1)*mod_magnus(Tbot(1))/rhoa

          if (myid==0) then
            write(*,*) 'SURFACE HUMIDITY: %RH, q = ',RHB,Tbot(2)
          endif

        else !(doing LES)

          RHB = tsfcc(2)

          tsfcc(2) = RHB/100.0*Mw/Ru/tsfcc(1)*mod_magnus(tsfcc(1))/rhoa

          if (myid==0) then
            write(*,*) 'SURFACE HUMIDITY: %RH, q = ',RHB,tsfcc(2)
          endif

        endif



      end subroutine change_RH_bcs_to_q

      function ran2(idum)
        integer :: idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
        real :: ran2,AM,EPS,RNMX
        PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     +      IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     +      IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
        INTEGER :: idum2,j,k,iv(NTAB),iy
        SAVE iv,iy,idum2
        DATA idum2/123456789/,iv/NTAB*0/,iy/0/

        if (idum <= 0) then
          idum=max(-idum,1)
          idum2 = idum
          do j = NTAB+8,1,-1
            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if (idum < 0) idum=idum+IM1
            if (j <= NTAB) iv(j) = idum
          enddo
          iy=iv(1)
        endif
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if (idum < 0) idum=idum+IM1
        k=idum2/IQ2
        idum2=IA2*(idum2-k*IQ2)-k*IR2
        if (idum2 < 0) idum2=idum2+IM2
        j = 1+iy/NDIV
        iy = iv(j) - idum2
        iv(j) = idum
        if (iy < 1) iy = iy+IMM1
        ran2=min(AM*iy,RNMX)
        return
      end function ran2

      subroutine linInterp(xa,ya,n,x,y)
        integer n
        real x, xper, y, xa(n), ya(n)
        integer i1, i2
        xper = x
      ! enforce periodic domain
        if (xper >= xa(n)) then
          xper = xper-xa(n)
        elseif (xper < 0) then
          xper = xa(n)+xper
        endif
      ! calc neighboring indices
        i1 = floor((n-1)*xper/xa(n)) + 1
        i2 = i1+1
        y = ya(i1) + (xper-xa(i1))*(ya(i2)-ya(i1))/(xa(i2)-xa(i1))
        return
      end

      subroutine linInterpSlope(xa,ya,n,x,dy_dx)
        integer n
        real x, xper, dy_dx, xa(n), ya(n)
        integer i1, i2
        xper = x
      ! enforce periodic domain
        if (xper >= xa(n)) then
          xper = xper-xa(n)
        elseif (xper < 0) then
          xper = xa(n)+xper
        endif
      ! calc neighboring indices
        i1 = floor((n-1)*xper/xa(n)) + 1
        i2 = i1+1

        dy_dx = (ya(i2)-ya(i1))/(xa(i2)-xa(i1))
        return
      end

      function mod_magnus(T)
        implicit none

      !Take in T in Kelvin and return saturation vapor pressure using function of Alduchov and Eskridge, 1996
        real,intent(in) :: T
        real :: mod_magnus

        mod_magnus=610.94*exp((17.6257*(T-273.15))/(243.04+(T-273.15)))


      end function mod_magnus
