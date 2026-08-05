module mo_gwd_cpp_interface
  use, intrinsic :: iso_c_binding
  implicit none

  interface
    subroutine c_ugwpv1_gsldrag_run(columns, layers, dtp, &
        ugrs, vgrs, tgrs, q1, prsl, prsi, prslk, phil, phii, del, &
        do_gsl_drag_ls_bl, do_gsl_drag_ss, do_gsl_drag_tofd, do_ngw_ec, &
        hprime, oc, theta, sigma, gamma, elvmax, clx, oa4, varss, &
        dx, xlat, area, &
        dudt_ogw, dvdt_ogw, dudt_ngw, dvdt_ngw, dtdt_ngw, &
        dudt_ofd, dvdt_ofd, tau_ogw, tau_ngw) &
        bind(C, name="c_ugwpv1_gsldrag_run")
      import :: c_size_t, c_double, c_int
      integer(c_size_t), value :: columns
      integer(c_size_t), value :: layers
      real(c_double), value    :: dtp
      real(c_double), intent(in)    :: ugrs(*)
      real(c_double), intent(in)    :: vgrs(*)
      real(c_double), intent(in)    :: tgrs(*)
      real(c_double), intent(in)    :: q1(*)
      real(c_double), intent(in)    :: prsl(*)
      real(c_double), intent(in)    :: prsi(*)
      real(c_double), intent(in)    :: prslk(*)
      real(c_double), intent(in)    :: phil(*)
      real(c_double), intent(in)    :: phii(*)
      real(c_double), intent(in)    :: del(*)
      integer(c_int), value         :: do_gsl_drag_ls_bl
      integer(c_int), value         :: do_gsl_drag_ss
      integer(c_int), value         :: do_gsl_drag_tofd
      integer(c_int), value         :: do_ngw_ec
      real(c_double), intent(in)    :: hprime(*)
      real(c_double), intent(in)    :: oc(*)
      real(c_double), intent(in)    :: theta(*)
      real(c_double), intent(in)    :: sigma(*)
      real(c_double), intent(in)    :: gamma(*)
      real(c_double), intent(in)    :: elvmax(*)
      real(c_double), intent(in)    :: clx(*)
      real(c_double), intent(in)    :: oa4(*)
      real(c_double), intent(in)    :: varss(*)
      real(c_double), intent(in)    :: dx(*)
      real(c_double), intent(in)    :: xlat(*)
      real(c_double), intent(in)    :: area(*)
      real(c_double), intent(out)   :: dudt_ogw(*)
      real(c_double), intent(out)   :: dvdt_ogw(*)
      real(c_double), intent(out)   :: dudt_ngw(*)
      real(c_double), intent(out)   :: dvdt_ngw(*)
      real(c_double), intent(out)   :: dtdt_ngw(*)
      real(c_double), intent(out)   :: dudt_ofd(*)
      real(c_double), intent(out)   :: dvdt_ofd(*)
      real(c_double), intent(out)   :: tau_ogw(*)
      real(c_double), intent(out)   :: tau_ngw(*)
    end subroutine c_ugwpv1_gsldrag_run

    subroutine c_ugwpv1_gsldrag_post_run(im, levs, ldiag_ugwp, dtf, &
        zobl, zlwb, zogw, tau_ogw, tau_ngw, du_ofdcol, du_oblcol, &
        tot_mtb, tot_ogw, tot_tofd, tot_ngw, &
        tot_zmtb, tot_zlwb, tot_zogw, &
        dudt_gw, dvdt_gw, dudt_obl, dudt_ofd, dudt_ogw, &
        du3dt_mtb, du3dt_tms, du3dt_ogw, &
        du3dt_ngw, dv3dt_ngw) &
        bind(C, name="c_ugwpv1_gsldrag_post_run")
      import :: c_size_t, c_double, c_int
      integer(c_size_t), value :: im
      integer(c_size_t), value :: levs
      integer(c_int), value    :: ldiag_ugwp
      real(c_double), value    :: dtf
      real(c_double), intent(in)    :: zobl(*)
      real(c_double), intent(in)    :: zlwb(*)
      real(c_double), intent(in)    :: zogw(*)
      real(c_double), intent(in)    :: tau_ogw(*)
      real(c_double), intent(in)    :: tau_ngw(*)
      real(c_double), intent(in)    :: du_ofdcol(*)
      real(c_double), intent(in)    :: du_oblcol(*)
      real(c_double), intent(out)   :: tot_mtb(*)
      real(c_double), intent(out)   :: tot_ogw(*)
      real(c_double), intent(out)   :: tot_tofd(*)
      real(c_double), intent(out)   :: tot_ngw(*)
      real(c_double), intent(out)   :: tot_zmtb(*)
      real(c_double), intent(out)   :: tot_zlwb(*)
      real(c_double), intent(out)   :: tot_zogw(*)
      real(c_double), intent(in)    :: dudt_gw(*)
      real(c_double), intent(in)    :: dvdt_gw(*)
      real(c_double), intent(in)    :: dudt_obl(*)
      real(c_double), intent(in)    :: dudt_ofd(*)
      real(c_double), intent(in)    :: dudt_ogw(*)
      real(c_double), intent(out)   :: du3dt_mtb(*)
      real(c_double), intent(out)   :: du3dt_tms(*)
      real(c_double), intent(out)   :: du3dt_ogw(*)
      real(c_double), intent(out)   :: du3dt_ngw(*)
      real(c_double), intent(out)   :: dv3dt_ngw(*)
    end subroutine c_ugwpv1_gsldrag_post_run
  end interface

end module mo_gwd_cpp_interface
