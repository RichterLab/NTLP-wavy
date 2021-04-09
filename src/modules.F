c ----------------------------------------------------------------------
      module pars
c ----------------------------------------------------------------------
      integer :: iti,itmax,imean,ihst,itape
      integer :: itstr,it_his,it_viz,iviz

      integer :: ismlt,ifree,iradup,
     + iupwnd,iskew,ibuoy,ifilt,
     + itcut,isubs,ibrcl,i_move_step,
     + iocean,method,idebug,iz_space,
     + ivis0,ifix_dt,new_vis,iflat,
     + i3d_surf,i_pm,lpbc,k_stab,
     + no_pgrad,ispray,iDNS,
     + ifields,ievap,i_tank,i_zeta,
     + inetcdf,i_pinertia

c
c -------- j_recl = record length in "bytes" for history file
c          k8     = kind parameter for integers in mpi_io routines
c
         integer, parameter :: j_recl=4, k8=8
c
c -------- SGI diamond
c
c        integer, parameter :: j_recl=1, k8=8
c
c -------- number of vars, size of problem in (x,y,z), max size, max processors
c
         integer, parameter :: nscl = 2, nvar = (4+nscl)
!         integer, parameter :: nxg1  = 256, nyg1  = 256, nzg1 = 128
!         integer, parameter :: maxnx = 256, maxny = 256, maxnz = 128
         integer, parameter :: nxg1  = 128, nyg1  = 128, nzg1 = 128
         integer, parameter :: maxnx = 128, maxny = 128, maxnz = 128
!         integer, parameter :: nxg1  = 64, nyg1  = 64, nzg1 = 64
!         integer, parameter :: maxnx = 64, maxny = 64, maxnz = 64

      !Boundary conditions and other parameters
      real :: Uo,Ttop(nscl),Tbot(nscl),qstar(nscl),tsfcc(nscl),RHbot
      real :: ugcont,vgcont,fcor,t00,grav,ustar_tank,waveslope
      real :: zi, zl, xl, yl, zw1, zo, dt_new
      integer :: isfc(nscl)

c
c ------------ leave following definitions as is
c
         integer, parameter :: maxnz1 = maxnz + 1,
     +                         maxnz2 = maxnz + 2,
     +                         maxnx2 = maxnx + 2,
     +                         maxny2 = maxny + 2
c ----------------------------------------------------------------------
      integer ::    nnx, nny, nnz, nxy, ncx, nnxp1, nnyp1, ncy,
     +              nnxp2, nnyp2, nnzp1, ivis, nnzm1, isize, krec,
     +              izs, ize, ixs, ixe, jxs, jxe, kxs, kxe,
     +              mxs, mxe, iss, ise, iys, iye, jys, jye, i2xs, i2xe,
     +              nnx_c, nny_c, nnz_c, imesh
c ----------------------------------------------------------------------
      character case*3
c ----------------------------------------------------------------------
      integer  ::   nvel, npre, nhis1, nprt,
     +              nhisp, nvelc,
     +              nviz_z, nviz_y,
     +              nviz_x, nviz_s, ngrd, nbndy,
     +              kfile, jfile, ibcl, ibcu,
     +              igrdr, imach, itn, it_his_nxt, it_viz_nxt
      logical ::    mnout, micut, mtape, mhis, msave,mtrans,
     +              l_root, l_debug, msave_v, mviz
c ----------------------------------------------------------------------
      real    ::    windm,u1xy,v1xy,t1xy(nscl),
     +              t10xy(nscl),au13m,au23m,aut3m(nscl),tsfcm(nscl),
     +              thstar(nscl), eavg(maxnz), tr_tau(0:maxnz),
     +              pxym(0:maxnz1), zi_min, a_exp_grd, a_exp_grd2,
     +              a_decay, p_surf
      integer ::    izi, iz_min, i_match
      real, allocatable ::
     +              wind(:,:), tau13m(:,:), tau23m(:,:),
     +              taut3m(:,:,:), t_grnd(:,:,:)
c ----------------------------------------------------------------------
      real ::       u_mn(0:maxnz1), v_mn(0:maxnz1),
     +              w_mn(0:maxnz1), t_mn(0:maxnz1,nscl)
c ----------------------------------------------------------------------
      real ::       dzw(0:maxnz2), dzu(0:maxnz2),
     +              dzw_i(0:maxnz2), dzu_i(0:maxnz2),
     +              weit(0:maxnz1), weitm(0:maxnz1)
c ----------------------------------------------------------------------
      real ::       t_factor, t_ref, c_rate, t_surf_i
c ----------------------------------------------------------------------
      real ::       dsl_z(0:maxnz1),
     +              xksurf, viscon, vise, almin_c, stabmin,
     +              ck, ceps, csmag, stab_c, vis_mean(0:maxnz), pow3
      integer ::    nmatch
c ----------------------------------------------------------------------
      real ::       zetas(3), gama(3), etas(4),
     +              umax,vmax,wmax, wabs, emax, vismax,
     +              cfl, tzero,
     +              ucfl, vcfl, wcfl,
     +              t_stage, t_stage2, t_stage3,
     +              dtalp
c ----------------------------------------------------------------------
      real ::       dpdx_con, dpdx_con_z
c ----------------------------------------------------------------------
      integer, parameter :: ip_accel=100
      integer ::    ip_max, ip_ramp
      real ::       ptol, p_iter, p_err, accel_min, accel_max,
     +              accel(ip_accel),
     +              avg_j(0:maxnz1), avg_j_i(0:maxnz1)
c ----------------------------------------------------------------------
      real, allocatable ::
     +              t11bc(:,:), t12bc(:,:), t13bc(:,:), t22bc(:,:),
     +              t23bc(:,:), t33bc(:,:), ustar(:,:), t1tbc(:,:,:),
     +              t2tbc(:,:,:), t3tbc(:,:,:), amol(:,:)
c ----------------------------------------------------------------------
      character*80  path_res, path_sav, path_his, path_prt,
     +              path_sav_v, path_sav_c, path_part,
     +              path_sav_p, path_sav_h, path_sav_part,
     +              path_sav_hp, path_sav_g, path_bndy,path_ran,
     +              path_histog, path_sav_hist
      character case_inp*3
      character*80 path_viz_xy, path_viz_xz, path_viz_yz, path_stuf,
     +             path_seed
c ----------------------------------------------------------------------
      integer ::    myid, numprocs, i_root,
     +              ziloc, myid_newvis, ncpu_s, ncpu_z, maxp
      integer, allocatable, dimension(:) ::
     +              ix_s, ix_e, jx_s, jx_e,
     +              kx_s, kx_e, mx_s, mx_e,
     +              iy_s, iy_e, jy_s, jy_e,
     +              is_s, is_e, iz_s, iz_e,
     +              i2x_s, i2x_e
      end module pars
c ----------------------------------------------------------------------
      module particles
      integer :: rproc,trproc,tproc,tlproc,lproc,blproc,bproc,brproc
      integer :: pr_r,pl_r,pt_r,pb_r,ptr_r,ptl_r,pbl_r,pbr_r
      integer :: pr_s,pl_s,pt_s,pb_s,ptr_s,ptl_s,pbl_s,pbr_s
      real :: ymin,ymax,zmin,zmax,xmax,xmin
      real, allocatable :: uext(:,:,:), vext(:,:,:), wext(:,:,:)
      real, allocatable :: u_t(:,:,:), v_t(:,:,:), w_t(:,:,:)
      real, allocatable :: Text(:,:,:),T_t(:,:,:)
      real, allocatable :: T2ext(:,:,:),T2_t(:,:,:)
      real, allocatable :: partTsrc(:,:,:),partTsrc_t(:,:,:)
      real, allocatable :: partcount_t(:,:,:),partsrc_t(:,:,:,:)
      real, allocatable :: vpsum_t(:,:,:,:),vpsqrsum_t(:,:,:,:)
      real, allocatable :: upwp_t(:,:,:),upwp(:,:,:)
      real, allocatable :: partcount(:,:,:),partsrc(:,:,:,:)
      real, allocatable :: partcountU_t(:,:,:),partcountU(:,:,:)
      real, allocatable :: partflux(:,:,:),partflux_t(:,:,:)
      real, allocatable :: partfluxup(:,:,:),partfluxup_t(:,:,:)
      real, allocatable :: partfluxdn(:,:,:),partfluxdn_t(:,:,:)
      real, allocatable :: vpsum(:,:,:,:),vpsqrsum(:,:,:,:)
      real, allocatable :: Tpsum(:,:,:),Tpsum_t(:,:,:)
      real, allocatable :: Tpsqrsum(:,:,:),Tpsqrsum_t(:,:,:)
      real, allocatable :: wpTpsum(:,:,:),wpTpsum_t(:,:,:)
      real, allocatable :: Tfsum(:,:,:),Tfsum_t(:,:,:)
      real, allocatable :: qfsum(:,:,:),qfsum_t(:,:,:)
      real, allocatable :: radsum(:,:,:),radsum_t(:,:,:)
      real, allocatable :: rad2sum(:,:,:),rad2sum_t(:,:,:)
      real, allocatable :: qstarsum(:,:,:),qstarsum_t(:,:,:)
      integer :: particletype,pad_diff
      integer :: numpart,tnumpart,ngidx
      integer :: num_reintro,tot_reintro,his_reintro
      integer :: iseed
      integer,parameter :: rbins=4
      real :: Rep_avg,part_grav(3),Vpmax,dz_U
      real :: vp_init(3),Tp_init,radius_init

      real :: rhoa,rhow,nuf
      real :: CpaCpp,Pra,Sc,Cpa,Cpp
      real :: Mw,Ru,Sal,Ms,Gam,Ion,Os

      integer, parameter :: histbins = 256,numhistvars=5
      integer :: hist_flag=0
      real :: hist_locations(histbins+2,histbins+2)
      real :: hist_restime(histbins+2,histbins+2)
      real :: hist_height(histbins+2,histbins+2)
      real :: hist_momentum(histbins+2,histbins+2)
      real :: bins(histbins+2,5)  !residence time, wave-relative origin, wave-relative death, max height


      !REMEMBER: IF ADDING ANYTHING, MUST UPDATE MPI DATATYPE!
      type :: particle
      integer :: pidx,procidx
      real :: vp(3),xp(3),xip(3),uf(3),xrhs(3),vrhs(3),Tp,Tprhs_s,Tf
      real :: Tprhs_L,radius,radrhs,qinf,qstar,m_s,xip_i(3),xp_0(3),t0
      real :: zpmax,xp_i(3),u0
      type(particle), pointer :: prev,next
      end type particle

      type(particle), pointer :: part,first_particle
      end module particles
c     ---------------------------------------------------------------------
      module wavy_pars
c
      integer, parameter ::
     +           i_ranp = 0
c
      real :: u10, angle, wind_x, wind_y, dkx, dky
c
c --------- data arrays
c
      real, allocatable ::
     +      omega(:,:), eta_wrk(:,:),
     +      etax(:,:), etay(:,:), eta0(:,:,:),
     +      s_wrk(:,:,:), v_wrk(:,:,:)
      end module wavy_pars
c ----------------------------------------------------------------------
      module tank_pars
c
      integer, parameter ::
     +           nx_tank = 2000
c
      real :: wave_l, c_speed, minwavez
c
c --------- data arrays
c
      real, allocatable ::
     +      z_tank(:), u_tank(:), w_tank(:),
     +      z_tank_f(:), u_tank_f(:), w_tank_f(:),
     +      z_tank_i(:), u_tank_i(:), w_tank_i(:),
     +      hdot_i(:), hx_i(:)
      end module tank_pars
c ----------------------------------------------------------------------
      module fields
        real, allocatable ::
     +      u(:,:,:), v(:,:,:), w(:,:,:), t(:,:,:,:), e(:,:,:),
     +      r1(:,:,:), r2(:,:,:), r3(:,:,:), r4(:,:,:,:), r5(:,:,:),
     +      uf(:,:,:), vf(:,:,:), wf(:,:,:),
     +      ufi(:,:,:), vfi(:,:,:), wfi(:,:,:), rj(:,:,:)
        real, allocatable ::
     +      ux(:,:,:), uy(:,:,:), vx(:,:,:), vy(:,:,:),
     +      wx(:,:,:), wy(:,:,:),
     +      p(:,:,:), ptop(:,:,:), p_old(:,:,:),
     +      vis_m(:,:,:), vis_s(:,:,:)
        real, allocatable ::
     +      ubc(:,:,:), vbc(:,:,:), wbc(:,:,:), tbc(:,:,:,:),
     +      ebc(:,:,:), pbc(:,:,:), pbc2(:,:,:),
     +      ufbc(:,:,:), vfbc(:,:,:), wfbc(:,:,:)
        real, allocatable ::
     +      xgrid(:), ygrid(:), zgrid_w(:,:,:,:), zgrid_u(:,:,:),
     +      gm(:,:,:,:), gm2(:,:,:,:), gm3(:,:,:,:), gms(:,:,:),
     +      bndy(:,:), hdot(:,:),
     +      w_orbit(:,:), u_orbit(:,:), v_orbit(:,:),
     +      bndy_t(:,:,:), hdot_t(:,:,:),
     +      z1b(:,:), z1b_2(:,:), z1b_3(:,:)
      end module fields
c ----------------------------------------------------------------------
      module fftwk
        real, allocatable :: trigx(:,:), trigc(:), trigcx(:)
      end module fftwk
c ----------------------------------------------------------------------
      module con_data
c ----------------------------------------------------------------------
        use pars, only : nscl
        type con_d
        sequence
        real ::  vk, vkin, vk74, vk74in,
     +           gcp, fcor_h, pi2,
     +           batagk, batag,
     +           cdbtm, dtdzf(nscl), dtjump, ugal, divgls,
     +           z1, utausv, dx, dy, dz,
     +           dt, dt1,
     +           fnxy, dzdz, dsl, c23, dtgama, dtzeta, xkmax,
     +           time, t_zero, zody, zody74,
     +           utau, wtsfc(nscl),
     +           uusfc, uvsfc, uwsfc, vvsfc, vwsfc, wwsfc, amonin,
     +           zol, hol, smal_e, sml_eg,drg_x,drg_y
        end type con_d
        type(con_d), target :: c_c
        real, pointer ::
     +           vk, vkin, vk74, vk74in,
     +           gcp, fcor_h, pi2,
     +           batagk, batag,
     +           cdbtm, dtdzf(:), dtjump, ugal, divgls,
     +           z1, utausv, dx, dy, dz,
     +           dt, dt1,
     +           fnxy, dzdz, dsl, c23, dtgama, dtzeta, xkmax,
     +           time, t_zero, zody, zody74,
     +           utau, wtsfc(:),
     +           uusfc, uvsfc, uwsfc, vvsfc, vwsfc, wwsfc, amonin,
     +           zol, hol, smal_e, sml_eg,drg_x,drg_y
      contains
         subroutine fill_cc
c
c --------------- pointer associations for constant variables
c
             vk     => c_c%vk
             vkin   => c_c%vkin
             vk74   => c_c%vk74
             vk74in => c_c%vk74in
             gcp    => c_c%gcp
             fcor_h => c_c%fcor_h
             pi2    => c_c%pi2
             batagk => c_c%batagk
             batag  => c_c%batag
             cdbtm  => c_c%cdbtm
             dtdzf  => c_c%dtdzf
             dtjump => c_c%dtjump
             ugal   => c_c%ugal
             divgls => c_c%divgls
             z1     => c_c%z1
             utausv => c_c%utausv
             dx     => c_c%dx
             dy     => c_c%dy
             dz     => c_c%dz
             dt     => c_c%dt
             dt1    => c_c%dt1
             fnxy   => c_c%fnxy
             dzdz   => c_c%dzdz
             dsl    => c_c%dsl
             c23    => c_c%c23
             dtgama => c_c%dtgama
             dtzeta => c_c%dtzeta
             xkmax  => c_c%xkmax
             time   => c_c%time
             t_zero => c_c%t_zero
             zody   => c_c%zody
             zody74 => c_c%zody74
             utau   => c_c%utau
             wtsfc  => c_c%wtsfc
             uusfc  => c_c%uusfc
             uvsfc  => c_c%uvsfc
             uwsfc  => c_c%uwsfc
             vvsfc  => c_c%vvsfc
             vwsfc  => c_c%vwsfc
             wwsfc  => c_c%wwsfc
             amonin => c_c%amonin
             zol    => c_c%zol
             hol    => c_c%hol
             smal_e => c_c%smal_e
             sml_eg => c_c%sml_eg
             drg_x  => c_c%drg_x
             drg_y  => c_c%drg_y
         return
         end subroutine fill_cc
      end module con_data
c ----------------------------------------------------------------------
      module con_stats
        use pars
        type con_s
        sequence
        real ::  wwsb(maxnz),engz(0:maxnz1),
     +           engsbz(0:maxnz1),
     +           englez(maxnz),uxym(0:maxnz1),
     +           vxym(0:maxnz1),wxym(0:maxnz1), pm(0:maxnz1),
     +           txym(0:maxnz1,nscl),divz(0:maxnz1),
     +           utle(maxnz,nscl), utsb(maxnz,nscl),
     +           vtle(maxnz,nscl), vtsb(maxnz,nscl),
     +           wtle(maxnz,nscl), wtsb(maxnz,nscl),
     +           wctle(maxnz,nscl), wt_tot(maxnz,nscl),
     +           z(0:maxnz1),zz(0:maxnz1),
     +           shrz(maxnz),buyz(maxnz),
     +           triz(maxnz),
     +           uwsb(maxnz),vwsb(maxnz),
     +           uwle(maxnz),vwle(maxnz),wwle(maxnz),
     +           uwcle(maxnz),vwcle(maxnz),
     +           uw_tot(maxnz),vw_tot(maxnz),
     +           drg_z_x(maxnz),drg_z_y(maxnz),
     +           drg_zs_x(maxnz),drg_zs_y(maxnz),
     +           wcube(maxnz), wfour(maxnz),
     +           tcube(maxnz,nscl),
     +           ups(maxnz), vps(maxnz),
     +           wps(maxnz), tps(maxnz,nscl),
     +           t_rprod(maxnz), t_wq(maxnz),
     +           t_wp(maxnz), t_tau(maxnz),
     +           t_tran(maxnz), t_buoy(maxnz),
     +           t_diss(maxnz), t_sprod(maxnz),
     +           zconc(maxnz),zconc_U(maxnz),
     +           vp1mean(maxnz),vp2mean(maxnz),vp3mean(maxnz),
     +           vp1msqr(maxnz),vp2msqr(maxnz),vp3msqr(maxnz),
     +           upwpm(maxnz),
     +           Tpmean(maxnz),Tpmsqr(maxnz),wpTpm(maxnz),
     +           Tfmean(maxnz),qfmean(maxnz),radmean(maxnz),
     +           rad2mean(maxnz),qstarmean(maxnz),pflux(maxnz),
     +           pfluxup(maxnz),pfluxdn(maxnz)
        real ::  xkn(maxnx),ykn(maxny),
     +           xk(maxnx),yk(maxny),
     +           xkk(maxnx),ykk(maxny),
     +           xks(maxnx2,maxny),wavexy(maxnx2,maxny)
        real ::  ug(maxnz),vg(maxnz),
     +           wls(maxnz),uls(maxnx)
        real ::  udrift,vdrift,
     +           stokesw,stokesa,
     +           stokess,stokes(maxnz1)
        real ::  dtg, dslg, dzg
        real ::  dfac(maxnz)
        end type con_s
        type(con_s), target :: c_s
        real, pointer ::
     +           wwsb(:), engz(:), engsbz(:),
     +           englez(:), uxym(:), vxym(:), wxym(:), pm(:),
     +           txym(:,:), divz(:), utle(:,:), utsb(:,:),
     +           vtle(:,:), vtsb(:,:), wtle(:,:), wtsb(:,:),
     +           wctle(:,:), wt_tot(:,:), z(:), zz(:), shrz(:), buyz(:),
     +           triz(:), uwsb(:), vwsb(:), uwle(:), vwle(:), wwle(:),
     +           uwcle(:), vwcle(:),
     +           uw_tot(:), vw_tot(:), drg_z_x(:), drg_z_y(:),
     +           drg_zs_x(:), drg_zs_y(:),
     +           wcube(:), wfour(:),
     +           tcube(:,:), ups(:), vps(:),
     +           wps(:), tps(:,:), t_rprod(:), t_wq(:),
     +           t_wp(:), t_tau(:), t_tran(:), t_buoy(:),
     +           t_diss(:), t_sprod(:),
     +           zconc(:),zconc_U(:),
     +           vp1mean(:),vp2mean(:),vp3mean(:),
     +           vp1msqr(:),vp2msqr(:),vp3msqr(:),
     +           upwpm(:),
     +           Tpmean(:),Tpmsqr(:),wpTpm(:),
     +           Tfmean(:),qfmean(:),radmean(:),
     +           rad2mean(:),qstarmean(:),pflux(:),
     +           pfluxup(:),pfluxdn(:)
        real, pointer ::
     +           xkn(:), ykn(:), xk(:), yk(:), xkk(:), ykk(:),
     +           xks(:,:), wavexy(:,:)
        real, pointer ::
     +           ug(:), vg(:), wls(:), uls(:)
        real, pointer ::
     +           udrift, vdrift, stokesw, stokesa,
     +           stokess, stokes(:)
        real, pointer ::
     +           dtg, dslg, dzg
        real, pointer ::
     +           dfac(:)
      contains
         subroutine fill_cs
c
c -------------- pointer association for stat arrays and get size
c                of stat arrays isize for history files
c
             isize = 0
             wwsb     => c_s%wwsb      ; isize = isize + size(wwsb)
             engz     => c_s%engz      ; isize = isize + size(engz)
             engsbz   => c_s%engsbz    ; isize = isize + size(engsbz)
             englez   => c_s%englez    ; isize = isize + size(englez)
             uxym     => c_s%uxym      ; isize = isize + size(uxym)
             vxym     => c_s%vxym      ; isize = isize + size(vxym)
             wxym     => c_s%wxym      ; isize = isize + size(wxym)
             pm       => c_s%pm        ; isize = isize + size(pm)
             txym     => c_s%txym      ; isize = isize + size(txym)
             divz     => c_s%divz      ; isize = isize + size(divz)
             utle     => c_s%utle      ; isize = isize + size(utle)
             utsb     => c_s%utsb      ; isize = isize + size(utsb)
             vtle     => c_s%vtle      ; isize = isize + size(vtle)
             vtsb     => c_s%vtsb      ; isize = isize + size(vtsb)
             wtle     => c_s%wtle      ; isize = isize + size(wtle)
             wctle    => c_s%wctle     ; isize = isize + size(wctle)
             wtsb     => c_s%wtsb      ; isize = isize + size(wtsb)
             wt_tot   => c_s%wt_tot    ; isize = isize + size(wt_tot)
             z        => c_s%z         ; isize = isize + size(z)
             zz       => c_s%zz        ; isize = isize + size(zz)
             shrz     => c_s%shrz      ; isize = isize + size(shrz)
             buyz     => c_s%buyz      ; isize = isize + size(buyz)
             triz     => c_s%triz      ; isize = isize + size(triz)
             uwsb     => c_s%uwsb      ; isize = isize + size(uwsb)
             vwsb     => c_s%vwsb      ; isize = isize + size(vwsb)
             uwle     => c_s%uwle      ; isize = isize + size(uwle)
             vwle     => c_s%vwle      ; isize = isize + size(vwle)
             wwle     => c_s%wwle      ; isize = isize + size(wwle)
             uwcle    => c_s%uwcle     ; isize = isize + size(uwcle)
             vwcle    => c_s%vwcle     ; isize = isize + size(vwcle)
             uw_tot   => c_s%uw_tot    ; isize = isize + size(uw_tot)
             vw_tot   => c_s%vw_tot    ; isize = isize + size(vw_tot)
             drg_z_x  => c_s%drg_z_x   ; isize = isize + size(drg_z_x)
             drg_z_y  => c_s%drg_z_y   ; isize = isize + size(drg_z_y)
             drg_zs_x => c_s%drg_zs_x  ; isize = isize + size(drg_zs_x)
             drg_zs_y => c_s%drg_zs_y  ; isize = isize + size(drg_zs_y)
             wcube    => c_s%wcube     ; isize = isize + size(wcube)
             wfour    => c_s%wfour     ; isize = isize + size(wfour)
             tcube    => c_s%tcube     ; isize = isize + size(tcube)
             ups      => c_s%ups       ; isize = isize + size(ups)
             vps      => c_s%vps       ; isize = isize + size(vps)
             wps      => c_s%wps       ; isize = isize + size(wps)
             tps      => c_s%tps       ; isize = isize + size(tps)
             t_rprod  => c_s%t_rprod   ; isize = isize + size(t_rprod)
             t_wq     => c_s%t_wq      ; isize = isize + size(t_wq)
             t_wp     => c_s%t_wp      ; isize = isize + size(t_wp)
             t_tau    => c_s%t_tau     ; isize = isize + size(t_tau)
             t_tran   => c_s%t_tran    ; isize = isize + size(t_tran)
             t_buoy   => c_s%t_buoy    ; isize = isize + size(t_buoy)
             t_diss   => c_s%t_diss    ; isize = isize + size(t_diss)
             t_sprod  => c_s%t_sprod   ; isize = isize + size(t_sprod)
             zconc   => c_s%zconc    ; isize = isize + size(zconc)
             zconc_U => c_s%zconc_U  ; isize = isize + size(zconc_U)
             vp1mean  => c_s%vp1mean ; isize = isize + size(vp1mean)
             vp2mean  => c_s%vp2mean ; isize = isize + size(vp2mean)
             vp3mean  => c_s%vp3mean ; isize = isize + size(vp3mean)
             vp1msqr  => c_s%vp1msqr ; isize = isize + size(vp1msqr)
             vp2msqr  => c_s%vp2msqr ; isize = isize + size(vp2msqr)
             vp3msqr  => c_s%vp3msqr ; isize = isize + size(vp3msqr)
             upwpm   => c_s%upwpm    ; isize = isize + size(upwpm)
             Tpmean  => c_s%Tpmean   ; isize = isize + size(Tpmean)
             Tpmsqr  => c_s%Tpmsqr   ; isize = isize + size(Tpmsqr)
             wpTpm   => c_s%wpTpm    ; isize = isize + size(wpTpm)
             Tfmean  => c_s%Tfmean   ; isize = isize + size(Tfmean)
             qfmean  => c_s%qfmean   ; isize = isize + size(qfmean)
             radmean  => c_s%radmean ; isize = isize + size(radmean)
             rad2mean => c_s%rad2mean; isize = isize + size(rad2mean)
             qstarmean => c_s%qstarmean; isize = isize + size(qstarmean)
             pflux => c_s%pflux; isize = isize + size(pflux)
             pfluxup => c_s%pfluxup; isize = isize + size(pfluxup)
             pfluxdn => c_s%pfluxdn; isize = isize + size(pfluxdn)
             xkn      => c_s%xkn
             ykn      => c_s%ykn
             xk       => c_s%xk
             yk       => c_s%yk
             xkk      => c_s%xkk
             ykk      => c_s%ykk
             xks      => c_s%xks
             wavexy   => c_s%wavexy
             ug       => c_s%ug
             vg       => c_s%vg
             wls      => c_s%wls
             uls      => c_s%uls
             udrift   => c_s%udrift
             vdrift   => c_s%vdrift
             stokesw  => c_s%stokesw
             stokesa  => c_s%stokesa
             stokess  => c_s%stokess
             stokes   => c_s%stokes
             dtg      => c_s%dtg
             dslg     => c_s%dslg
             dzg      => c_s%dzg
             dfac     => c_s%dfac
         return
         end subroutine fill_cs
      end module con_stats
c ----------------------------------------------------------------------
