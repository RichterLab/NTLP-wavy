
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
            zosdy = alog(z1b(i,j)/zos + 1.0)
            call mo(ibuoy,isfc,ismlt,vk,batag,zody,zosdy,z1b(i,j),
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
      subroutine mo(ibuoy,isfc,ismlt,vk,batag,zody,zosdy,z1,wind,
     +      utau,nscl,t1xy,tsfcc,qstar,zeta,cd,ch,amonin,RHB)
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
            t_fac = 1.0/amax1((zosdy - psis),0.01)
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
          dnom = zosdy*vk74in
        else
          dnom      = amax1((zosdy-psis),0.01)*vk74in
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
        integer values(8)

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

        if (use_seed==0) then
          call date_and_time(VALUES=values)
          idum = -(values(8)+values(7)+values(6))
        else
          idum = -use_seed
        endif

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
        real, allocatable :: xx(:), zz_tank(:), us(:), ws(:),
     +      y2(:), u2(:), w2(:), xles(:)

c---------- count number of points in wave file
        lut    = 89
        close(lut)
        open(lut,file=wave_file)
        npts_t = -4 ! don't count top 4 comment lines
        do
          read (lut,*, end=10) ! break at end of file
          npts_t = npts_t + 1
        enddo
10      close (lut)

c --------- allocate surface arrays
        allocate(z_tank(nnx), u_tank(nnx), w_tank(nnx),
     +      z_tank_f(nnx), u_tank_f(nnx), w_tank_f(nnx),
     +      z_tank_i(nnx), u_tank_i(nnx), w_tank_i(nnx),
     +      hdot_i(nnx), hx_i(nnx), hx2_i(nnx),
     +      xx(npts_t), zz_tank(npts_t), us(npts_t), ws(npts_t),
     +      y2(npts_t), u2(npts_t), w2(npts_t), xles(nnx))

c--------- read in wave file
        open (lut, file = wave_file) !re-open to read from start
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


        wave_k  = pi2/wave_l ! calc main wave number

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

        real, allocatable :: xx(:), zz_tank(:), us(:), ws(:),
     +      y2(:), u2(:), w2(:), xles(:)

c---------- count number of points in wave file
        lut    = 89
        close(lut)
        open(lut,file=wave_file)
        npts_t = -4 ! don't count top 4 comment lines
        do
          read (lut,*, end=10) ! break at end of file
          npts_t = npts_t + 1
        enddo
10      close (lut)

c --------- allocate surface arrays
        allocate(z_tank(nnx), u_tank(nnx), w_tank(nnx),
     +      z_tank_f(nnx), u_tank_f(nnx), w_tank_f(nnx),
     +      z_tank_i(nnx), u_tank_i(nnx), w_tank_i(nnx),
     +      hdot_i(nnx), hx_i(nnx), hx2_i(nnx),
     +      xx(npts_t), zz_tank(npts_t), us(npts_t), ws(npts_t),
     +      y2(npts_t), u2(npts_t), w2(npts_t), xles(nnx))

c--------- read in wave file
        open (lut, file = wave_file) !re-open to read
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

        real gam, f_main

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

        gam = 0.0   ! set gamma here (0.7*0.9) 
        wave_k = 2*pi/wave_l  ! set dominant wavelength wave_l in params.in
        f_main = sqrt(grav*wave_k)
        do j=iys,iye
          do i = 1, nnx
            u_orbit(i,j) = u_orbit(i,j)
     +          + gam * (1 + cos(wave_k*bndy(i,j) - f_main*t_bndy))
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