import sys

def apply_fix_sascnvn():
    filename = 'physics/CONV/SAS/sascnvn.F'
    with open(filename, 'r') as f:
        lines = f.readlines()

    out = []
    i = 0
    done_use = False
    done_vars = False
    done_init = False
    done_mass_calc = False
    done_fb = False
    done_end = False

    while i < len(lines):
        line = lines[i]

        # 1. Add USE
        if not done_use and 'use funcphys , only : fpvs' in line:
            out.append(line)
            out.append('      use tracer_transport_mod, only : calc_column_mass,                &\n')
            out.append('     &    transport_tracer_flux_form, transport_tracer_hole_filling,    &\n')
            out.append('     &    transport_tracer_advective, apply_mass_fixer\n')
            done_use = True
            i += 1
            continue

        # 2. Local variables
        if not done_vars and 'real(kind=kind_phys) ps(im),      del(im,km), prsl(im,km)' in line:
            out.append(line)
            out.append('      real(kind=kind_phys) mass_init_q1(im), mass_init_qlc(im),         &\n')
            out.append('     &                     mass_init_qli(im), zero_sink(im),            &\n')
            out.append('     &                     mnet(im,km)\n')
            out.append('      integer :: transport_opt = 1\n')
            out.append('      logical :: use_mass_fixer = .true.\n')
            done_vars = True
            i += 1
            continue

        # 3. zero_sink init
        if not done_init and 'cnvflg(i) = .true.' in line:
            out.append(line)
            out.append('        zero_sink(i) = 0.0\n')
            done_init = True
            i += 1
            continue

        # 4. Initial mass calc
        if not done_mass_calc and 'del  = delp  * 0.001' in line:
            out.append(line)
            out.append('!************************************************************************\n')
            out.append('      call calc_column_mass(im, km, cnvflg, q1, delp, mass_init_q1)\n')
            out.append('      call calc_column_mass(im, km, cnvflg, qlc, delp, mass_init_qlc)\n')
            out.append('      call calc_column_mass(im, km, cnvflg, qli, delp, mass_init_qli)\n')
            out.append('!************************************************************************\n')
            done_mass_calc = True
            i += 1
            continue

        # 5. Feedback section
        if not done_fb and 'delubar(i) = 0.' in line:
            out.append(line)
            i += 1
            while 'qcond(i) = 0.' not in lines[i]:
                out.append(lines[i])
                i += 1
            out.append(lines[i]) # qcond(i) = 0.
            i += 1
            out.append(lines[i]) # enddo
            i += 1

            # Now insert mnet and transport
            out.append('      do i = 1, im\n')
            out.append('        if (cnvflg(i)) then\n')
            out.append('          do k = 1, km - 1\n')
            out.append('            mnet(i,k) = eta(i,k) * xmb(i) -                             &\n')
            out.append('     &                  edto(i) * etad(i,k) * xmb(i)\n')
            out.append('          enddo\n')
            out.append('          mnet(i,km) = 0.0\n')
            out.append('        endif\n')
            out.append('      enddo\n')
            out.append('      select case(transport_opt)\n')
            out.append('        case(1)\n')
            out.append('          call transport_tracer_flux_form(im, km, dt2, cnvflg, delp,     &\n')
            out.append('     &                                    mnet, q1)\n')
            out.append('          call transport_tracer_flux_form(im, km, dt2, cnvflg, delp,     &\n')
            out.append('     &                                    mnet, qlc)\n')
            out.append('          call transport_tracer_flux_form(im, km, dt2, cnvflg, delp,     &\n')
            out.append('     &                                    mnet, qli)\n')
            out.append('        case(2)\n')
            out.append('          call transport_tracer_hole_filling(im, km, dt2, cnvflg, delp,  &\n')
            out.append('     &                                       mnet, q1)\n')
            out.append('          call transport_tracer_hole_filling(im, km, dt2, cnvflg, delp,  &\n')
            out.append('     &                                       mnet, qlc)\n')
            out.append('          call transport_tracer_hole_filling(im, km, dt2, cnvflg, delp,  &\n')
            out.append('     &                                       mnet, qli)\n')
            out.append('        case(3)\n')
            out.append('          call transport_tracer_advective(im, km, dt2, cnvflg, delp,     &\n')
            out.append('     &                                    mnet, q1)\n')
            out.append('          call transport_tracer_advective(im, km, dt2, cnvflg, delp,     &\n')
            out.append('     &                                    mnet, qlc)\n')
            out.append('          call transport_tracer_advective(im, km, dt2, cnvflg, delp,     &\n')
            out.append('     &                                    mnet, qli)\n')
            out.append('      end select\n')
            done_fb = True
            continue

        # 6. Comment out old updates
        if 'q1(i,k) = q1(i,k) + dellaq(i,k) * xmb(i) * dt2' in line:
            out.append('! ' + line)
            i += 1
            continue
        if 'q1(i,k) = q1(i,k) + qevap(i)' in line:
            out.append('! ' + line)
            i += 1
            continue
        if 't1(i,k) = to(i,k)' in line and 'reset the updated state variables' in lines[i-3]:
            out.append('! ' + line)
            i += 1
            continue
        if 'q1(i,k) = qo(i,k)' in line and 'reset the updated state variables' in lines[i-4]:
            out.append('! ' + line)
            i += 1
            continue
        if 'u1(i,k) = uo(i,k)' in line and 'reset the updated state variables' in lines[i-5]:
            out.append('! ' + line)
            i += 1
            continue
        if 'v1(i,k) = vo(i,k)' in line and 'reset the updated state variables' in lines[i-6]:
            out.append('! ' + line)
            i += 1
            continue

        # 7. Mass fixer before return
        if not done_end and line.strip() == 'return' and i > len(lines) - 50:
             out.append('      if (use_mass_fixer) then\n')
             out.append('        call apply_mass_fixer(im, km, cnvflg, delp, rn,                 &\n')
             out.append('     &                        mass_init_q1, q1)\n')
             out.append('        call apply_mass_fixer(im, km, cnvflg, delp, zero_sink,          &\n')
             out.append('     &                        mass_init_qlc, qlc)\n')
             out.append('        call apply_mass_fixer(im, km, cnvflg, delp, zero_sink,          &\n')
             out.append('     &                        mass_init_qli, qli)\n')
             out.append('      endif\n')
             out.append(line)
             done_end = True
             i += 1
             continue

        out.append(line)
        i += 1

    with open(filename, 'w') as f:
        f.writelines(out)

def apply_fix_shalcnv():
    filename = 'physics/CONV/SAS/shalcnv.F'
    with open(filename, 'r') as f:
        lines = f.readlines()

    out = []
    i = 0
    done_use = False
    done_vars = False
    done_init = False
    done_mass_calc = False
    done_fb = False
    done_end = False

    while i < len(lines):
        line = lines[i]

        if not done_use and 'use funcphys , only : fpvs' in line:
            out.append(line)
            out.append('      use tracer_transport_mod, only : calc_column_mass,                &\n')
            out.append('     &    transport_tracer_flux_form, transport_tracer_hole_filling,    &\n')
            out.append('     &    transport_tracer_advective, apply_mass_fixer\n')
            done_use = True
            i += 1
            continue

        if not done_vars and 'real(kind=kind_phys) cincr, cincrmax, cincrmin' in line:
            out.append('      real(kind=kind_phys) mass_init_q1(im), mass_init_qlc(im),         &\n')
            out.append('     &                     mass_init_qli(im), zero_sink(im),            &\n')
            out.append('     &                     mnet(im,km)\n')
            out.append('      integer :: transport_opt = 1\n')
            out.append('      logical :: use_mass_fixer = .true.\n')
            out.append(line)
            done_vars = True
            i += 1
            continue

        if not done_init and 'cnvflg(i) = .true.' in line:
            out.append(line)
            out.append('        zero_sink(i) = 0.0\n')
            done_init = True
            i += 1
            continue

        if not done_mass_calc and 'del  = delp  * 0.001' in line:
            out.append(line)
            out.append('!************************************************************************\n')
            out.append('      call calc_column_mass(im, km, cnvflg, q1, delp, mass_init_q1)\n')
            out.append('      call calc_column_mass(im, km, cnvflg, qlc, delp, mass_init_qlc)\n')
            out.append('      call calc_column_mass(im, km, cnvflg, qli, delp, mass_init_qli)\n')
            out.append('!************************************************************************\n')
            done_mass_calc = True
            i += 1
            continue

        if not done_fb and 'delubar(i) = 0.' in line:
            out.append(line)
            i += 1
            while 'qcond(i) = 0.' not in lines[i]:
                out.append(lines[i])
                i += 1
            out.append(lines[i]) # qcond(i) = 0.
            i += 1
            out.append(lines[i]) # enddo
            i += 1

            out.append('      do i = 1, im\n')
            out.append('        if (cnvflg(i)) then\n')
            out.append('          do k = 1, km - 1\n')
            out.append('            mnet(i,k) = eta(i,k) * xmb(i)\n')
            out.append('          enddo\n')
            out.append('          mnet(i,km) = 0.0\n')
            out.append('        endif\n')
            out.append('      enddo\n')
            out.append('      select case(transport_opt)\n')
            out.append('        case(1)\n')
            out.append('          call transport_tracer_flux_form(im, km, dt2, cnvflg, delp,     &\n')
            out.append('     &                                    mnet, q1)\n')
            out.append('          call transport_tracer_flux_form(im, km, dt2, cnvflg, delp,     &\n')
            out.append('     &                                    mnet, qlc)\n')
            out.append('          call transport_tracer_flux_form(im, km, dt2, cnvflg, delp,     &\n')
            out.append('     &                                    mnet, qli)\n')
            out.append('        case(2)\n')
            out.append('          call transport_tracer_hole_filling(im, km, dt2, cnvflg, delp,  &\n')
            out.append('     &                                       mnet, q1)\n')
            out.append('          call transport_tracer_hole_filling(im, km, dt2, cnvflg, delp,  &\n')
            out.append('     &                                       mnet, qlc)\n')
            out.append('          call transport_tracer_hole_filling(im, km, dt2, cnvflg, delp,  &\n')
            out.append('     &                                       mnet, qli)\n')
            out.append('        case(3)\n')
            out.append('          call transport_tracer_advective(im, km, dt2, cnvflg, delp,     &\n')
            out.append('     &                                    mnet, q1)\n')
            out.append('          call transport_tracer_advective(im, km, dt2, cnvflg, delp,     &\n')
            out.append('     &                                    mnet, qlc)\n')
            out.append('          call transport_tracer_advective(im, km, dt2, cnvflg, delp,     &\n')
            out.append('     &                                    mnet, qli)\n')
            out.append('      end select\n')
            done_fb = True
            continue

        if 'q1(i,k) = q1(i,k) + dellaq(i,k) * xmb(i) * dt2' in line:
            out.append('! ' + line)
            i += 1
            continue
        if 'q1(i,k) = q1(i,k) + qevap(i)' in line:
            out.append('! ' + line)
            i += 1
            continue

        if not done_end and line.strip() == 'return' and i > len(lines) - 50:
            out.append('      if (use_mass_fixer) then\n')
            out.append('        call apply_mass_fixer(im, km, cnvflg, delp, zero_sink,          &\n')
            out.append('     &                        mass_init_q1, q1)\n')
            out.append('        call apply_mass_fixer(im, km, cnvflg, delp, zero_sink,          &\n')
            out.append('     &                        mass_init_qlc, qlc)\n')
            out.append('        call apply_mass_fixer(im, km, cnvflg, delp, zero_sink,          &\n')
            out.append('     &                        mass_init_qli, qli)\n')
            out.append('      endif\n')
            out.append(line)
            done_end = True
            i += 1
            continue

        out.append(line)
        i += 1

    with open(filename, 'w') as f:
        f.writelines(out)

apply_fix_sascnvn()
apply_fix_shalcnv()
