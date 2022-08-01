module netcdf_io
implicit none

integer :: ncid
integer :: time_dimid,zu_vid,zu_dimid,s_dimid
integer :: time_vid,dt_vid
integer :: utau_vid,uwsfc_vid,vwsfc_vid,drg_x_vid,drg_y_vid
integer :: tnumpart_vid
integer :: zw_vid,zw_dimid
integer :: uxym_vid,vxym_vid,wxym_vid,txym_vid,RHxym_vid
integer :: ups_vid,vps_vid,wps_vid,tps_vid
integer :: wtle_vid,wtsb_vid
integer :: uwle_vid,uwsb_vid
integer :: vwle_vid,vwsb_vid
integer :: zconc_vid
integer :: vp1mean_vid,vp2mean_vid,vp3mean_vid
integer :: vp1msqr_vid,vp2msqr_vid,vp3msqr_vid
integer :: Tpmean_vid,Tpmsqr_vid
integer :: Tfmean_vid,qfmean_vid
integer :: radmean_vid,rad2mean_vid
integer :: dimids(1),dimids_zu(2),dimids_zw(2),dimids_zu_s(3),dimids_zw_s(3)
integer :: his_counter
character(len=80) :: path_netcdf_his

CONTAINS

subroutine netcdf_check(status)
      use netcdf
      implicit none
      integer, intent ( in) :: status

      if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
      end if

end subroutine netcdf_check

subroutine netcdf_init
      use netcdf
      use pars
      implicit none

      path_netcdf_his = trim(adjustl(path_his))//"history.nc"

      call netcdf_check( nf90_create(path_netcdf_his,nf90_clobber,ncid))

      call netcdf_check( nf90_def_dim(ncid, "time",NF90_UNLIMITED, time_dimid) )

      call netcdf_check( nf90_def_dim(ncid, "zu",nnz, zu_dimid) )
      call netcdf_check( nf90_def_dim(ncid, "zw",nnz+1, zw_dimid) )

      call netcdf_check( nf90_def_dim(ncid,"nscl",nscl,s_dimid) )

      dimids = (/ time_dimid /)
      dimids_zu = (/ zu_dimid, time_dimid/)
      dimids_zw = (/ zw_dimid, time_dimid/)
      dimids_zu_s = (/ zu_dimid, s_dimid, time_dimid/)
      dimids_zw_s = (/ zw_dimid, s_dimid, time_dimid/)

      call netcdf_check( nf90_def_var(ncid,"time",NF90_REAL,dimids,time_vid) )
      call netcdf_check( nf90_put_att(ncid,time_vid,"title","Simulation time") )

      call netcdf_check( nf90_def_var(ncid, "dt", NF90_REAL, dimids,dt_vid) )
      call netcdf_check( nf90_put_att(ncid,dt_vid,"title","Model time step") )

      call netcdf_check( nf90_def_var(ncid, "utau", NF90_REAL, dimids,utau_vid) )
      call netcdf_check( nf90_put_att(ncid,utau_vid,"title","Friction velocity ustar (no pressure)") )

      call netcdf_check( nf90_def_var(ncid, "uwsfc", NF90_REAL, dimids,uwsfc_vid) )
      call netcdf_check( nf90_put_att(ncid,uwsfc_vid,"title","Lower surface friction stress in x-dir (no pressure)") )

      call netcdf_check( nf90_def_var(ncid, "vwsfc", NF90_REAL, dimids,vwsfc_vid) )
      call netcdf_check( nf90_put_att(ncid,vwsfc_vid,"title","Lower surface friction stress in y-dir (no pressure)") )

      call netcdf_check( nf90_def_var(ncid, "drg_x", NF90_REAL, dimids,drg_x_vid) )
      call netcdf_check( nf90_put_att(ncid,drg_x_vid,"title","Form stress in x due to pressure") )

      call netcdf_check( nf90_def_var(ncid, "drg_y", NF90_REAL, dimids,drg_y_vid) )
      call netcdf_check( nf90_put_att(ncid,drg_y_vid,"title","Form stress in y due to pressure") )

      call netcdf_check( nf90_def_var(ncid, "tnumpart", NF90_REAL, dimids,tnumpart_vid) )
      call netcdf_check( nf90_put_att(ncid,tnumpart_vid,"title","Total number of particles") )

      call netcdf_check( nf90_def_var(ncid, "zu", NF90_REAL, dimids_zu,zu_vid) )
      call netcdf_check( nf90_put_att(ncid,zu_vid,"title","z levels at u-points") )

      call netcdf_check( nf90_def_var(ncid, "zw", NF90_REAL, dimids_zw,zw_vid) )
      call netcdf_check( nf90_put_att(ncid,zw_vid,"title","z levels at w-points") )

      call netcdf_check( nf90_def_var(ncid,"uxym",NF90_REAL, dimids_zu,uxym_vid) )
      call netcdf_check( nf90_put_att(ncid,uxym_vid,"title","Horiz. avg. u vel") )

      call netcdf_check( nf90_def_var(ncid,"vxym",NF90_REAL, dimids_zu,vxym_vid) )
      call netcdf_check( nf90_put_att(ncid,vxym_vid,"title","Horiz. avg. v vel") )

      call netcdf_check( nf90_def_var(ncid,"wxym",NF90_REAL, dimids_zw,wxym_vid) )
      call netcdf_check( nf90_put_att(ncid,wxym_vid,"title","Horiz. avg. w vel") )

      call netcdf_check( nf90_def_var(ncid,"txym",NF90_REAL, dimids_zu_s,txym_vid) )
      call netcdf_check( nf90_put_att(ncid,txym_vid,"title","Horiz. avg. scalars") )

      call netcdf_check( nf90_def_var(ncid,"RHxym",NF90_REAL, dimids_zu,RHxym_vid) )
      call netcdf_check( nf90_put_att(ncid,RHxym_vid,"title","Horiz. avg. relative humidity") )

      call netcdf_check( nf90_def_var(ncid,"ups",NF90_REAL, dimids_zu,ups_vid) )
      call netcdf_check( nf90_put_att(ncid,ups_vid,"title","Fluctuating velocity <u'^2>") )

      call netcdf_check( nf90_def_var(ncid,"vps",NF90_REAL, dimids_zu,vps_vid) )
      call netcdf_check( nf90_put_att(ncid,vps_vid,"title","Fluctuating velocity <v'^2>") )

      call netcdf_check( nf90_def_var(ncid,"wps",NF90_REAL, dimids_zw,wps_vid) )
      call netcdf_check( nf90_put_att(ncid,wps_vid,"title","Fluctuating velocity <w'^2>") )

      call netcdf_check( nf90_def_var(ncid,"tps",NF90_REAL, dimids_zu_s,tps_vid) )
      call netcdf_check( nf90_put_att(ncid,tps_vid,"title","Fluctuating scalars <t'^2>") )

      call netcdf_check( nf90_def_var(ncid,"uwle",NF90_REAL, dimids_zw,uwle_vid) )
      call netcdf_check( nf90_put_att(ncid,uwle_vid,"title","Resolved <u'w'>") )

      call netcdf_check( nf90_def_var(ncid,"uwsb",NF90_REAL, dimids_zw,uwsb_vid) )
      call netcdf_check( nf90_put_att(ncid,uwsb_vid,"title","Subgrid <u'w'>") )

      call netcdf_check( nf90_def_var(ncid,"vwle",NF90_REAL, dimids_zw,vwle_vid) )
      call netcdf_check( nf90_put_att(ncid,vwle_vid,"title","Resolved <v'w'>") )

      call netcdf_check( nf90_def_var(ncid,"vwsb",NF90_REAL, dimids_zw,vwsb_vid) )
      call netcdf_check( nf90_put_att(ncid,vwsb_vid,"title","Subgrid <v'w'>") )

      call netcdf_check( nf90_def_var(ncid,"wtle",NF90_REAL, dimids_zw_s,wtle_vid) )
      call netcdf_check( nf90_put_att(ncid,wtle_vid,"title","Resolved <w't'>") )

      call netcdf_check( nf90_def_var(ncid,"wtsb",NF90_REAL, dimids_zw_s,wtsb_vid) )
      call netcdf_check( nf90_put_att(ncid,wtsb_vid,"title","Subgrid <w't'>") )

      call netcdf_check( nf90_def_var(ncid,"zconc",NF90_REAL, dimids_zu,zconc_vid) )
      call netcdf_check( nf90_put_att(ncid,zconc_vid,"title","Particle concentration") )

      call netcdf_check( nf90_def_var(ncid,"vp1mean",NF90_REAL, dimids_zu,vp1mean_vid) )
      call netcdf_check( nf90_put_att(ncid,vp1mean_vid,"title","Horiz. avg. particle velocity u") )

      call netcdf_check( nf90_def_var(ncid,"vp2mean",NF90_REAL, dimids_zu,vp2mean_vid) )
      call netcdf_check( nf90_put_att(ncid,vp2mean_vid,"title","Horiz. avg. particle velocity v") )

      call netcdf_check( nf90_def_var(ncid,"vp3mean",NF90_REAL, dimids_zu,vp3mean_vid) )
      call netcdf_check( nf90_put_att(ncid,vp3mean_vid,"title","Horiz. avg. particle velocity w") )

      call netcdf_check( nf90_def_var(ncid,"vp1msqr",NF90_REAL, dimids_zu,vp1msqr_vid) )
      call netcdf_check( nf90_put_att(ncid,vp1msqr_vid,"title","Mean squared particle velocity <u^2>") )

      call netcdf_check( nf90_def_var(ncid,"vp2msqr",NF90_REAL, dimids_zu,vp2msqr_vid) )
      call netcdf_check( nf90_put_att(ncid,vp2msqr_vid,"title","Mean squared particle velocity <v^2>") )

      call netcdf_check( nf90_def_var(ncid,"vp3msqr",NF90_REAL, dimids_zu,vp3msqr_vid) )
      call netcdf_check( nf90_put_att(ncid,vp3msqr_vid,"title","Mean squared particle velocity <w^2>") )

      call netcdf_check( nf90_def_var(ncid,"Tpmean",NF90_REAL, dimids_zu,Tpmean_vid) )
      call netcdf_check( nf90_put_att(ncid,Tpmean_vid,"title","Horiz. avg. particle temp") )

      call netcdf_check( nf90_def_var(ncid,"Tpmsqr",NF90_REAL, dimids_zu,Tpmsqr_vid) )
      call netcdf_check( nf90_put_att(ncid,Tpmsqr_vid,"title","Mean squared particle temp <Tp^2>") )

      call netcdf_check( nf90_def_var(ncid,"Tfmean",NF90_REAL, dimids_zu,Tfmean_vid) )
      call netcdf_check( nf90_put_att(ncid,Tfmean_vid,"title","Horiz. avg. fluid temp at particle") )

      call netcdf_check( nf90_def_var(ncid,"qfmean",NF90_REAL, dimids_zu,qfmean_vid) )
      call netcdf_check( nf90_put_att(ncid,qfmean_vid,"title","Horiz. avg. fluid qv at particle") )

      call netcdf_check( nf90_def_var(ncid,"radmean",NF90_REAL, dimids_zu,radmean_vid) )
      call netcdf_check( nf90_put_att(ncid,radmean_vid,"title","Horiz. avg. particle radius") )

      call netcdf_check( nf90_def_var(ncid,"rad2mean",NF90_REAL, dimids_zu,rad2mean_vid) )
      call netcdf_check( nf90_put_att(ncid,rad2mean_vid,"title","Mean squared particle radius <rp^2>") )

      call netcdf_check( nf90_enddef(ncid) )

      his_counter = 1

end subroutine netcdf_init

subroutine write_his_netcdf
      use netcdf
      use pars
      use fields
      use con_data
      use con_stats
      use particles

      implicit none
      real :: tmp(0:nnz),tmp_s(0:nnz,nscl)

      call netcdf_check( nf90_put_var(ncid, time_vid, real(time),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, dt_vid, real(dt),start=(/his_counter/)) )

      call netcdf_check( nf90_put_var(ncid, utau_vid, real(utau),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, uwsfc_vid, real(uwsfc),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, vwsfc_vid, real(vwsfc),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, drg_x_vid, real(drg_x),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, drg_y_vid, real(drg_y),start=(/his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, tnumpart_vid, real(tnumpart),start=(/his_counter/)) )

      call netcdf_check( nf90_put_var(ncid, zu_vid, real(zz(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid, zw_vid, real(z(0:nnz)),start=(/1, his_counter/)) )

      call netcdf_check( nf90_put_var(ncid,uxym_vid,real(uxym(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,vxym_vid,real(vxym(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,wxym_vid,real(wxym(0:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,txym_vid,real(txym(1:nnz,1:nscl)),start=(/1,1,his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,RHxym_vid,real(RHxym(1:nnz)),start=(/1, his_counter/)) )

      call netcdf_check( nf90_put_var(ncid,ups_vid,real(ups(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,vps_vid,real(vps(1:nnz)),start=(/1, his_counter/)) )

      tmp(0) = 0.0
      tmp(1:nnz) = wps(1:nnz)
      call netcdf_check( nf90_put_var(ncid,wps_vid,real(tmp),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,tps_vid,real(tps(1:nnz,1:nscl)),start=(/1,1,his_counter/)) )

      tmp(0) = 0
      tmp(1:nnz) = uwle(1:nnz)
      call netcdf_check( nf90_put_var(ncid,uwle_vid,real(tmp),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,uwsb_vid,real(uwsb(nnz)),start=(/1, his_counter/)) )

      tmp(0) = 0
      tmp(1:nnz) = vwle(1:nnz)
      call netcdf_check( nf90_put_var(ncid,vwle_vid,real(tmp),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,vwsb_vid,real(vwsb(nnz)),start=(/1, his_counter/)) )

      tmp_s(0,1:nscl) = 0
      tmp_s(1:nnz,1:nscl) = wtle(1:nnz,1:nscl)
      call netcdf_check( nf90_put_var(ncid,wtle_vid,real(tmp_s),start=(/1,1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,wtsb_vid,real(wtsb(nnz,1:nscl)),start=(/1,1, his_counter/)) )

      call netcdf_check( nf90_put_var(ncid,zconc_vid,real(zconc(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,vp1mean_vid,real(vp1mean(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,vp2mean_vid,real(vp2mean(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,vp3mean_vid,real(vp3mean(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,vp1msqr_vid,real(vp1msqr(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,vp2msqr_vid,real(vp2msqr(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,vp3msqr_vid,real(vp3msqr(1:nnz)),start=(/1, his_counter/)) )

      call netcdf_check( nf90_put_var(ncid,Tpmean_vid,real(Tpmean(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,Tpmsqr_vid,real(Tpmsqr(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,Tfmean_vid,real(Tfmean(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,qfmean_vid,real(qfmean(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,radmean_vid,real(radmean(1:nnz)),start=(/1, his_counter/)) )
      call netcdf_check( nf90_put_var(ncid,rad2mean_vid,real(rad2mean(1:nnz)),start=(/1, his_counter/)) )

      his_counter = his_counter + 1

end subroutine write_his_netcdf
subroutine close_his_netcdf
      use netcdf
      use pars
      implicit none

      call netcdf_check( nf90_close(ncid) )

end subroutine close_his_netcdf
subroutine open_his_netcdf
      use netcdf
      use pars
      implicit none

      call netcdf_check( nf90_open(path_netcdf_his,NF90_WRITE,ncid) )

end subroutine open_his_netcdf
subroutine netcdf_res
!https://www.unidata.ucar.edu/software/netcdf/docs-fortran/f90-use-of-the-netcdf-library.html#f90-writing-data-in-an-existing-netcdf-dataset
      use netcdf
      use pars
      implicit none

      integer :: dimids(1),dimids_zu(2),dimids_zw(2),dimids_zu_s(3),dimids_zw_s(3)
      integer :: Ntime

      path_netcdf_his = trim(adjustl(path_seed))//"history.nc"

      call netcdf_check( nf90_open(path_netcdf_his,NF90_WRITE,ncid) )

      call netcdf_check (nf90_inq_dimid(ncid,'time',time_dimid) )
      call netcdf_check (nf90_inq_dimid(ncid,'zu',zu_dimid) )
      call netcdf_check (nf90_inq_dimid(ncid,'zw',zw_dimid) )
      call netcdf_check (nf90_inq_dimid(ncid,'nscl',s_dimid) )

      dimids = (/ time_dimid /)
      dimids_zu = (/ zu_dimid, time_dimid/)
      dimids_zw = (/ zw_dimid, time_dimid/)
      dimids_zu_s = (/ zu_dimid, s_dimid, time_dimid/)
      dimids_zw_s = (/ zw_dimid, s_dimid, time_dimid/)

      !Get the length of the unlimited time dimension
      call netcdf_check( nf90_inquire_dimension(ncid=ncid,dimid=time_dimid,len=Ntime) )
      his_counter = Ntime + 1

      if (myid==0) write(*,*) 'Restart his_counter for netCDF = ',his_counter

!!! Single quantities
      call netcdf_check( nf90_inq_varid(ncid,"time",time_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"dt",dt_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"utau",utau_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"uwsfc",uwsfc_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"vwsfc",vwsfc_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"drg_x",drg_x_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"drg_y",drg_y_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"tnumpart",tnumpart_vid) )
!!! Profiles
      call netcdf_check( nf90_inq_varid(ncid,"zu",zu_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"zw",zw_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"uxym",uxym_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"vxym",vxym_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"wxym",wxym_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"txym",txym_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"RHxym",RHxym_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"ups",ups_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"vps",vps_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"wps",wps_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"tps",tps_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"uwle",uwle_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"uwsb",uwsb_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"vwle",vwle_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"vwsb",vwsb_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"wtle",wtle_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"wtsb",wtsb_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"zconc",zconc_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"vp1mean",vp1mean_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"vp2mean",vp2mean_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"vp3mean",vp3mean_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"vp1msqr",vp1msqr_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"vp2msqr",vp2msqr_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"vp3msqr",vp3msqr_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"Tpmean",Tpmean_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"Tpmsqr",Tpmsqr_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"Tfmean",Tfmean_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"qfmean",qfmean_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"radmean",radmean_vid) )
      call netcdf_check( nf90_inq_varid(ncid,"rad2mean",rad2mean_vid) )

end subroutine netcdf_res

end module netcdf_io
