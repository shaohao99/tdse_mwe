!
!... All input parameters are documented in file Initialize.TDSE.f90.

program prop_TDSE
use PropParams_TDSE
use Initialize_TDSE
use TDSE

implicit none

include 'mpif.h'
integer status(mpi_status_size)
!... other mpi integer variables declared in PropParams.TDSE.f90

integer :: error = 0
integer n, nrpt, nstart(nnr), nstop(nnr)
integer i, iw, it, ir
integer iz1, iz2, izstop2
real tstart,tstop,lastl,totlost
real zz,znext,dens0,e_pt
real ph_new(nnt), ph_old(nnt)
real simpson,time_av,spec_av
real specfun(nnr),emax(nnr)
real :: dens_z(nnz) = 0., pt(nnt) = 0.
real :: prob2(nnt) = 1., grnd(nnt) = 1.
double complex :: x_ft(nnt) = zero, e_init(nnt) = zero
double complex :: x_ft1(nnt) = zero, e_init1(nnt) = zero
double complex :: ew(nnw,nnr) = zero, et(nnt,nnr) = zero
double complex :: gt(nnt) = zero
double complex :: gh_tot(nnw,nnr) = zero, gh(nnw,nnr) = zero
double complex :: g_tot(nnw,nnr) = zero, g(nnw,nnr) = zero
double complex :: e(nnw,0:nnr+1) = zero
double complex :: eold(nnw,0:nnr+1) = zero
double complex :: eh(nnw,0:nnr+1) = zero
double complex :: ehold(nnw,0:nnr+1) = zero
character*24 FMT1
character*27 FMT2

call mpi_init(ierr)
call mpi_comm_size(mpi_comm_world,nprocs,ierr)
call mpi_comm_rank(mpi_comm_world,myid,ierr)

!... useful constants
pi = 4.*atan(1._8)   ! for propagation
pie = 4.*atan(1._8)  ! for TDSE

!... set initial conditions
if (myid == 0) then
   call init_params(error)         ! for propagation
   call TDSE_inputs(atom,error)    ! for TDSE
endif
call mpi_bcast(error,1,mpi_integer,0,mpi_comm_world,ierr)
if (error == 1) then
   call mpi_finalize(ierr)
   stop
endif

!... broadcast all input variables
call mpi_bcast(zl,1,mpi_double_precision,0,mpi_comm_world,ierr)
call mpi_bcast(nz,1,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(rm,1,mpi_double_precision,0,mpi_comm_world,ierr)
call mpi_bcast(nr,1,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(ngob,1,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(tau,1,mpi_double_precision,0,mpi_comm_world,ierr)
call mpi_bcast(npulse,1,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(nmult,1,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(phiabs,1,mpi_double_precision,0,mpi_comm_world,ierr)
call mpi_bcast(i0,1,mpi_double_precision,0,mpi_comm_world,ierr)
call mpi_bcast(lambda,1,mpi_double_precision,0,mpi_comm_world,ierr)
call mpi_bcast(bc,1,mpi_double_precision,0,mpi_comm_world,ierr)
call mpi_bcast(b,1,mpi_double_precision,0,mpi_comm_world,ierr)
call mpi_bcast(zfoc,1,mpi_double_precision,0,mpi_comm_world,ierr)
call mpi_bcast(i1,1,mpi_double_precision,0,mpi_comm_world,ierr)
call mpi_bcast(wmult,1,mpi_double_precision,0,mpi_comm_world,ierr)
call mpi_bcast(tau1,1,mpi_double_precision,0,mpi_comm_world,ierr)
call mpi_bcast(phi_t,1,mpi_double_precision,0,mpi_comm_world,ierr)
call mpi_bcast(bc1,1,mpi_double_precision,0,mpi_comm_world,ierr)
call mpi_bcast(b1,1,mpi_double_precision,0,mpi_comm_world,ierr)
call mpi_bcast(zfoc1,1,mpi_double_precision,0,mpi_comm_world,ierr)
call mpi_bcast(ip,1,mpi_double_precision,0,mpi_comm_world,ierr)
call mpi_bcast(density,1,mpi_double_precision,0,mpi_comm_world,ierr)
call mpi_bcast(nz1,1,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(sa_step,1,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(atom,2,mpi_character,0,mpi_comm_world,ierr)
call mpi_bcast(nrSA,1,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(nrgob,1,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(nl,1,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(nr_ion,1,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(l0,1,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(m0,1,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(tratio,1,mpi_integer,0,mpi_comm_world,ierr)
call mpi_bcast(psi0(1:nrSA),nrSA,mpi_double_precision,0,mpi_comm_world,ierr)
call mpi_bcast(pot,3*kr,mpi_double_precision,0,mpi_comm_world,ierr)
call mpi_bcast(drx,1,mpi_double_precision,0,mpi_comm_world,ierr)

call mpi_bcast(i2,1,mpi_double_precision,0,mpi_comm_world,ierr)  ! shaohao 2011.9
call mpi_bcast(tau2,1,mpi_double_precision,0,mpi_comm_world,ierr)  
call mpi_bcast(phi_t2,1,mpi_double_precision,0,mpi_comm_world,ierr)  
call mpi_bcast(b2,1,mpi_double_precision,0,mpi_comm_world,ierr)  
call mpi_bcast(zfoc2,1,mpi_double_precision,0,mpi_comm_world,ierr)  
call mpi_bcast(ohhg,4,mpi_double_precision,0,mpi_comm_world,ierr)  
call mpi_bcast(ep_hhg,4,mpi_double_precision,0,mpi_comm_world,ierr)  
call mpi_bcast(phi_hhg,4,mpi_double_precision,0,mpi_comm_world,ierr)  

!... n is actual number of time and frequency points for propagation
call init_grids(n,error)
if (error == 1) then
   call mpi_finalize(ierr)
   stop
endif
call init_pulse(n,tstart,tstop,x_ft,x_ft1,error)
if (error == 1) then
   call mpi_finalize(ierr)
   stop
endif
!... nstart and nstop are integration limits
nstart(1:nr) = NINT((tstart-tt(1))/dt) + 1
nstop(1:nr) = n-NINT((tt(n)-tstop)/dt) - 1

nsteps = tratio*n   ! for TDSE
!... initialize arrays for TDSE calculations
call TDSE_setup

! cdt = i*dt/2 in atomic units for TDSE
cdt = (ci/2.)*(dt/tratio)*(lambda/c)/2.4188956e-17

!... Define initial gas profile dens_z
if (nz1 <= 0) then    ! flat density profile
   forall (iz1=1:nz) dens_z(iz1) = density
else                  ! trapezoid density profile
   !... Note inserting -0. is used to convert to type real arithmetic
   forall (iz1=1:nz1) dens_z(iz1) = ((iz1-0.)/nz1)*density
   dens_z(nz1+1:nz-nz1) = density
   forall (iz1=nz-nz1+1:nz) dens_z(iz1) = ((nz-iz1+1.)/nz1)*density
endif

!... initialize FFT's
call zffti(n,wsave)

!... Calculate spectrum of initial pulse
call FFTf(n,x_ft,e_init)
call FFTf(n,x_ft1,e_init1)  ! second color field

!... Define laser beam in space. Result is returned in ew
zz = z(1)
call init_beam(n,zz,e_init,e_init1,ew)

!... for propagation (technically only the master processor needs this)
do ir=1,nr
   e(1:n,ir) = sqrt(rr(ir))*ew(1:n,ir)
enddo

!... Back-transform E(w,r) to time-domain E(t,r)
do ir=1,nr
   x_ft(1:n) = ew(1:n,ir)
   call FFTb(n,x_ft,et(:,ir))
enddo

!... Master processor writes initial files
if (myid == 0) then
   !... write initial spectrum
   open(9,file='spec.init',status='unknown')
   do iw=1,n
      specfun(1:nr) = abs(e(iw,1:nr))**2
      spec_av = simpson(specfun,1,nr,dr)
      write(9,'(4e15.6)') w(iw)/2./pi,spec_av,(abs(ew(iw,1)))**2, &
                atan2(imag(ew(iw,1)),real(ew(iw,1)))
   enddo
   close(unit=9)

   !... Write field at entrance to medium
   !... First write headers to indicate which harmonics are included
   open(14,file='zprofile.harm',access='append',status='unknown')
     write(FMT1,'(A8,I2,A14)') '(A1,21x,', nharm, '(A1,I2.2,12x))'
     write(14,FMT1) 'z', ('H',NINT(w(H_out(i))/2./pi), i=1,nharm)
     close(14,status='keep')
   call out_zprof(zz,nr,n,nharm,H_out,e,eh)

   open(15,file='rprofiles.harm',access='append',status='unknown')
     write(FMT2,'(A8,I2,A17)') '(A1,21x,', nharm, '(A1,I2.2,12x),A2)'
     write(15,FMT2) 'r', ('H',NINT(w(H_out(i))/2./pi), i=1,nharm), 'z '
     close(15,status='keep')
   call out_rprof(zz,nr,n,nharm,H_out,e,eh)

   !... Write radially integrated time profile
   open(92,file='time.init.integr',status='unknown')
   do it=1,n
      specfun(1:nr) = rr(1:nr)*abs(et(it,1:nr))**2
      spec_av = simpson(specfun,1,nr,dr)
      write(92,'(2e15.6)') tt(it),spec_av
   enddo
   close(unit=92)

   !... Time-profile of fundamental on axis, evolution with z
!   open(95,file='time.fund.z',access='append',status='unknown')
   open(95,file='time.fund.z',status='unknown')
   write(95,*) lambda,n,tmin,tmax
   write(95,*) iw0,dt,nz,zl

endif    ! end Master writing files

!.. for perturbation theory; I_pt is an input known only to the master
if (myid == 0) e_pt = 27.44925e-9*sqrt(I_pt)

!..................................................................
!
!........ Loop on z - march through the medium ..............
do iz1 = 1,nz-1,sa_step
   if (myid == 0) write(6,*) 'iz1:', iz1
   zz = z(iz1)
   znext = z(iz1+1)
   dens0 = dens_z(iz1)

   !... Find r_pt, beyond which the peak intensity is less than the user
   !... input I_pt W/cm2 - we use perturbation theory for r > r_pt
   if (myid == 0) then
      emax = maxval(abs(real(et)), 1)    ! maximum electric field at each r
      nrpt = nr
      do ir = nr,1,-1
         if (emax(ir) >= e_pt) then
            nrpt = ir
            exit
         endif
      end do
      write(*,*) 'nrpt,emax(nrpt)',nrpt,emax(nrpt)
   endif
   call mpi_bcast(nrpt,1,mpi_integer,0,mpi_comm_world,ierr)
   call mpi_bcast(ew,nnt*nnr,mpi_double_complex,0,mpi_comm_world,ierr)

   if (myid /= 0) then
      do ir=myid,nrpt-1,nprocs-1
         ! find integration limits
         call find_limits(ew(1:n,ir),ir,n,nstart(ir),nstop(ir))
         ! interpolate for increased time resolution, et returned in eTDSE
         call interpolate(ew(1:n,ir),n,tt(1),dw)
         ! calculate single-atom response using TDSE
         call run_TDSE(pt,prob2,grnd,gt,n,nstart(ir),nstop(ir),lastl,totlost)
         ! Calculate frequency dependent source terms for march routine
         call calc_source_TDSE(n,nstart(ir),nstop(ir),ir,pt,gt,dens0,&
                               prob2,gh,g)

         ! Communication to master
         call mpi_send(gh(1:n,ir),n,mpi_double_complex,0,ir,mpi_comm_world,ierr)
         call mpi_send(g(1:n,ir),n,mpi_double_complex,0,ir+nrpt,mpi_comm_world,ierr)
 
         if (ir == 1) then
            ! ... TDSE outputs for diagnostic, evolution with z
            open(100,file='TDSE.out.z',access='append',status='unknown')
            write(100,'(5e15.6)') zz,totlost,lastl,1.-prob2(n),1.-grnd(n)
            close(100)
         endif
      end do
   endif

   if (myid == 0) then
      ir = nrpt
      call find_limits(ew(1:n,ir),ir,n,nstart(ir),nstop(ir))
      call interpolate(ew(1:n,ir),n,tt(1),dw)
      call run_TDSE(pt,prob2,grnd,gt,n,nstart(ir),nstop(ir),lastl,totlost)
      call calc_source_TDSE(n,nstart(ir),nstop(ir),ir,pt,gt,dens0,&
                            prob2,gh,g)

      do ir=1,nrpt-1
         ! Communication from slaves
         call mpi_recv(gh(1:n,ir),n,mpi_double_complex,mpi_any_source,ir, &
                       mpi_comm_world,status,ierr)
         call mpi_recv(g(1:n,ir),n,mpi_double_complex,mpi_any_source,ir+nrpt, &
                       mpi_comm_world,status,ierr)
      end do

      do ir = nrpt+1,nr
!        Set plasma term to zero, use perturbation theory for harmonic term
         g(1:n,ir) = zero
         gh(1:n,ir)=gh(1:n,nrpt) &
                    *(abs(ew(iw0,ir))/abs(ew(iw0,nrpt)))**(abs(w(1:n))/2./pi)
      end do
      gh_tot = gh + 1.e9*g
      g_tot = g + 1.e-9*gh
!      g_tot = g   ! for no absorption/not self-consistent

      ! Use E(w,r) from before as E_old(w,r); these are array assignments
      eold = e
      ehold = eh

      ! March one step, find E(w,r) at z+dz using E_old(w,r) and g_tot(w,r)
      call march(nr,ngob,n,eold,g_tot,e)
      ! March the harmonic field to z+dz, using Eh_old(w,r) and gh_tot(w,r)
      ! This is only necessary for no absorption case
      call march(nr,ngob,n,ehold,gh_tot,eh)

      call out_zprof(znext,nr,n,nharm,H_out,e,eh)

      ! march some steps without re-calculating sa_response or ionization
      ! protect against reaching end of medium
      izstop2 = min(iz1+sa_step-1,nz-1)
      do iz2 = iz1+1,izstop2,1
!         write(6,*) iz2, 'iz2 loop'
         zz = z(iz2)
         znext = z(iz2+1)
         dens0 = dens_z(iz2)

         ! Scale the atomic response by the new density 
         ! Incorporate the change in the geometrical phase in harmonic response
         do ir=1,nr
            ph_new(1:n)=atan2(imag(e(1:n,ir)),real(e(1:n,ir)))
            ph_old(1:n)=atan2(imag(eold(1:n,ir)),real(eold(1:n,ir)))
            g(1:n,ir)=dens_z(iz2)/dens_z(iz2-1)*g(1:n,ir)
            gh(1:n,ir)=dens_z(iz2)/dens_z(iz2-1)*gh(1:n,ir) &
                       *exp(ci*(ph_new(1:n)-ph_old(1:n)))
         enddo
         gh_tot = gh + 1.e9*g
         g_tot = g + 1.e-9*gh
!         g_tot = g   ! for no absorption/not self-consistent

         ! Use E(w,r) from before as E_old(w,r)
         eold = e
         ehold = eh

         ! March one step, find E(w,r) at z+dz, using E_old(w,r) and G(w,r)
         call march(nr,ngob,n,eold,g_tot,e)
         ! March the harmonic field to z+dz, using Eh_old(w,r) and Gh(w,r)
         ! This is only necessary for no absorption case
         call march(nr,ngob,n,ehold,gh_tot,eh)

         call out_zprof(znext,nr,n,nharm,H_out,e,eh)
      end do  ! end iz2 loop

      ! Back-transform E(w,r) to E(t,r), needed at top of z loop
      do ir=1,nr
         ew(1:n,ir) = e(1:n,ir)/sqrt(rr(ir))
         x_ft(1:n) = ew(1:n,ir)
         call FFTb(n,x_ft,et(:,ir))
      enddo

      ! Master processor writes fields on axis (these are in frequency domain)
      call out_rprof(znext,nr,n,nharm,H_out,e,eh)
      call out_zspec(znext,nr,n,e,eh)

      ! Write on-axis time-profile, evolution with z
      do it=1,n
         write(95,'(5e15.6)') tt(it),abs(et(it,1))**2, &
                 real(et(it,1)),imag(et(it,1)),znext
      enddo
      write(95,*)
   endif ! for master
   ! Repeat (go to next z)
enddo

!... End of medium 

!... Master processor writes output files
if (myid == 0) then
   ! Write E(w,r) at exit of medium
   open(82,file='ew.exit',status='unknown',form='unformatted')
   write(82) tmin,tmax,dt,iw0
   write(82) nr,rm
   do ir=1,nr
      do iw=1,n/2
         write(82) ew(iw,ir)
      enddo
   enddo
   close(unit=82)

   ! Write Eh(w,r) at exit of medium
   open(83,file='ewharm.exit',status='unknown',form='unformatted')
   write(83) tmin,tmax,dt,iw0
   write(83) nr,rm
   do ir=1,nr
      do iw=1,n/2
         write(83) eh(iw,ir)/sqrt(rr(ir))
      enddo
   enddo
   close(unit=83)

   ! Write time profile (integrated and on axis) at exit of medium
   open(81,file='time.final',status='unknown')
   do it=1,n
      specfun(1:nr) = abs(et(it,1:nr))**2*rr(1:nr)
      time_av = simpson(specfun,1,nr,dr)
 	  write(81,'(4e15.6)') tt(it),time_av,abs(et(it,1))**2, &
 	                       atan2(imag(et(it,1)),real(et(it,1)))
   enddo
   close(unit=81)

   ! Write spectrum (integrated and on axis) at exit of medium
   ! Choose radii between 1 and last nrpt to write out
   ir = INT((nrpt-1.)/3.)
   open(91,file='spec.final',status='unknown')
   do iw=n/2+2,n
      specfun(1:nr) = abs(e(iw,1:nr))**2
      spec_av = simpson(specfun,1,nr,dr)
 	  write(91,'(7e15.6)') w(iw)/2./pi,spec_av,abs(e(iw,1))**2/rr(1), &
              atan2(imag(e(iw,1)),real(e(iw,1))), &
              abs(e(iw,ir))**2/rr(ir),abs(e(iw,2*ir))**2/rr(2*ir), &
              abs(e(iw,3*ir))**2/rr(3*ir)
   enddo
   close(unit=91)

   open(92,file='spec.harm.final',status='unknown')
   do iw=n/2+2,n
      specfun(1:nr) = abs(eh(iw,1:nr))**2
      spec_av = simpson(specfun,1,nr,dr)
      write(92,'(7e15.6)') w(iw)/2./pi,spec_av,abs(eh(iw,1))**2/rr(1), &
               atan2(imag(eh(iw,1)),real(eh(iw,1))), &
               abs(eh(iw,ir))**2/rr(ir),abs(eh(iw,2*ir))**2/rr(2*ir), &
               abs(eh(iw,3*ir))**2/rr(3*ir)
   enddo
   close(unit=92)
   write(*,*) 'Master finished writing files'
endif

call mpi_finalize(ierr)
stop
end program prop_TDSE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calc_source_TDSE(n,nstart,nstop,ir,pt,gt,dens_at,&
                            prob2,gh,g)
!... calculate source terms for march routine
   use PropParams_TDSE
   implicit none
   integer, intent(in) :: n, ir, nstart, nstop
   integer it, iw
   real dens_at, nfilt, a2
   real pt(nnt), prob2(nnt), filt(nnt)
   double complex gt(nnt), gh(nnw,nnr), g(nnw,nnr)
   double complex x_ft(nnt)

!... Transform acceleration a(t) to frequency domain, multiply with sqrt(r)
!... units conversion - 8.47842e-28 is conversion to C*cm: e*a0*100
!... extra factor of 2 comes from double occupancy of p orbital
!... dens_at is scaling by density of neutral atoms
!... For harmonic source term, multiply dipole with -mu_0*w^2
!... Since a(w) = -d(w)*w^2, multiply a(w) by mu_0*units conversion, gives
!... factor of 1.25664e-8(H/cm) * (4.13414e16)**2
!... The non-linear polarization is returned in units of C/cm2

   ! Hann window
   filt = 0._8
   nfilt = nstop-nstart+2.
   forall (it=nstart:nstop) filt(it) = 0.5*(1.-cos(2.*pi*(it-nstart+1.)/nfilt))

   a2 = 2.*8.47842e-28*dens_at*1.25664e-8
   x_ft(1:n) = a2*pt(1:n)*filt(1:n)
   call FFTf(n,x_ft,gh(:,ir))
   gh(1:n,ir) = gh(1:n,ir)*sqrt(rr(ir))*(4.13414e16)**2

!... Define G(t,r) = e(t)*dens_e*e^2/(eps0*m_e*c^2)
!... where dens_e = dens_at*(1.-(prob2(k))**2)
!... gt comes in containing the electric field in atomic units
   gt(1:n) = 3.54114e-12*(gt(1:n)*5.14225)*dens_at*(1.-(prob2(1:n))**2)

!... Transform G(t,r) to frequency domain, G(w,r)
!... Multiply with sqrt(r)
   call FFTf(n,gt,g(:,ir))
   g(1:n,ir) = g(1:n,ir) * sqrt(rr(ir))
end subroutine calc_source_TDSE



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real function simpson(f,i1,i2,del)
    implicit none
    real f(*), del, t
    integer i1, i2, i3, i
    i3=i2-2
    t=0.0
    do i=i1,i3,2
       t=t+del/3.*(f(i)+4.*f(i+1)+f(i+2))
    enddo
    simpson=t
    return
end
