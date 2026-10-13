#if defined(ROW_LAND)
#define SEA_P .true.
#define SEA_U .true.
#define SEA_V .true.
#define SEA_Q .true.
#elif defined(ROW_ALLSEA)
#define SEA_P allip(j).or.ip(i,j).ne.0
#define SEA_U alliu(j).or.iu(i,j).ne.0
#define SEA_V alliv(j).or.iv(i,j).ne.0
#define SEA_Q alliq(j).or.iq(i,j).ne.0
#else
#define SEA_P ip(i,j).ne.0
#define SEA_U iu(i,j).ne.0
#define SEA_V iv(i,j).ne.0
#define SEA_Q iq(i,j).ne.0
#endif
!  Scientific results obtained based on this parameterization scheme and code must cite the following references.
!
!    Zhang, J., Zhang, Z., & Qiu, B. (2023).  
!    Parameterizing Submesoscale Vertical Buoyancy Flux by Simultaneously Considering Baroclinic Instability and 
!    Strain‐Induced Frontogenesis. Geophysical Research Letters, 50(8), e2022GL102292.
!*******************************************************************************************************************
      module mod_z23vbf
      use mod_xc    ! HYCOM communication interface
      use mod_pipe  ! HYCOM debugging interface
      implicit none   
!
! --- module for tsadvc and related routines
!
      private !! default is private
      public  :: z23_vbf  
      integer, parameter, dimension (0:4) :: &
      mbdy_advtyp = (/ 2,     & !PCM
                         5,     & !MPDATA
                         5,     & !FCT2
                         0,     & !N/A
                         5 /)     !FCT4


#if defined(RELO)
      real, save, allocatable, dimension(:,:) :: &
#else
      real,    save, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: &
#endif
      dUdx, dUde, dVdx, dVde, tmp_psi
 

#if defined(RELO)
      real, save, allocatable, dimension(:,:,:) :: &
#else
      real,    save, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kdm) :: &
#endif
      dTdx_z, dTde_z, dSdx_z, dSde_z, dBdx_z, dBde_z, alpha_rho, dp_m, z_p, VBF_para, VHF_para
      contains 
      subroutine z23_vbf (m, n) 
        use mod_cb_arrays  ! HYCOM saved arrays
        implicit none
        integer m,n    
        
	real*8 :: C_sub = 8.0
        real*8 :: B_sub = 4.0
        real*8 :: f_min = 7.0e-6 ! limit smallest f at near 2-3°N
        real*8 :: gorho0 = 9.806/1000.0 ! g/rho0 (approximate) rho0取1025是因为ROMS模式中取的该值
        real*8 :: cff,cff1,S_sub,twc_tmp,swc_tmp,dBdx,dBde,dSdx,dSde, &
                  dTdx,dTde,Hz_sum,alpha_bar,mld_m,t_bar,s_bar,twc_hbar, &
                  swc_hbar,z_H,dz,uz,twc_nodz,swc_nodz
        integer i,j,k,mbdy_a,margin,ktr,htr

# include "stmt_fns.h"

#if defined(RELO)
        if (.not.allocated(dUdx)) then
          allocate( &
                 dUdx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                 dUde(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                 dVdx(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
                 dVde(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), &
              tmp_psi(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ) 
          call mem_stat_add( 5*(idm+2*nbdy)*(jdm+2*nbdy) )
          allocate( &
              dTdx_z(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kdm), &   
              dTde_z(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kdm), &  
              dSdx_z(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kdm), &   
              dSde_z(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kdm), &   
              dBdx_z(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kdm), &  
              dBde_z(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kdm), &
                dp_m(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kdm), &
                 z_p(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kdm), &
            VBF_para(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kdm), &
            VHF_para(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kdm), &
           alpha_rho(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,1:kdm) )    
          call mem_stat_add( 11*(idm+2*nbdy)*(jdm+2*nbdy)*kdm ) 
                 dUdx = 0.0
                 dUde = 0.0
                 dVdx = 0.0
                 dVde = 0.0
              tmp_psi = 0.0
                 dp_m = 0.0
                  z_p = 0.0
             VBF_para = 0.0
             VHF_para = 0.0
               dTdx_z = 0.0
               dTde_z = 0.0
               dSdx_z = 0.0
               dSde_z = 0.0
               dBdx_z = 0.0
               dBde_z = 0.0
            alpha_rho = 0.0
       endif
#endif	                   
! Calculate Mesoscale Strain Rate
        mbdy_a = mbdy_advtyp(1)  ! = 5
        margin = mbdy_a - 1

        do k=1,kk
          do j=1-margin,jj+margin
            do i=1-margin,ii+margin
              dUdx(i,j)=0.0
              dUde(i,j)=0.0
              dVdx(i,j)=0.0
              dVde(i,j)=0.0
              tmp_psi(i,j)=0.0
              if (SEA_Q) then
                dUde(i,j) = 1/scuy(i,j)*u(i,j,k,m)-1/scuy(i,j-1)*u(i,j-1,k,m)
                dVdx(i,j) = 1/scvx(i,j)*v(i,j,k,m)-1/scvx(i-1,j)*v(i-1,j,k,m)
              endif
              if (SEA_P) then
                dUdx(i,j) = 1/scpx(i,j)*(u(i+1,j,k,m)-u(i,j,k,m))
                dVde(i,j) = 1/scpy(i,j)*(v(i,j+1,k,m)-v(i,j,k,m))
              endif
              if (.not.( ip(i,j).ne.0 .and. ip(i+1,j).ne.0 .and. &
                 ip(i-1,j).ne.0 .and. ip(i,j+1).ne.0 .and. &
                 ip(i,j-1).ne.0 )) then
                dUdx(i,j)=0.0
                dUde(i,j)=0.0
                dVdx(i,j)=0.0
                dVde(i,j)=0.0
              endif
              tmp_psi(i,j) = (scq2i(i,j)*(dVdx(i,j) + dUde(i,j)))**2
            enddo
          enddo
          do j=1-margin,jj+margin
            do i=1-margin,ii+margin
              alpha_rho(i,j,k) = sqrt(0.25 *                              &
             &               (tmp_psi(i  ,j) + tmp_psi(i  ,j+1) +        &
             &               tmp_psi(i+1,j) + tmp_psi(i+1,j+1)) +        &
             &               (dUdx(i,j) - dVde(i,j))**2)
            enddo
          enddo
        enddo  
   

!$OMP PARALLEL DO PRIVATE(j,i,k,cff,cff1,dTdx,dSdx,dBdx,dTde,dSde,dBde) &  
!$OMP          SCHEDULE(STATIC,jblk)
        do k=1,kk
          do j=1-margin,jj+margin
            do i=1-margin,ii+margin
              dTdx_z(i,j,k)=0.0
              dTde_z(i,j,k)=0.0
              dSdx_z(i,j,k)=0.0
              dSde_z(i,j,k)=0.0
              dBdx_z(i,j,k)=0.0
              dBde_z(i,j,k)=0.0
              if (SEA_P) then
                cff =0.5*(1/scpx(i-1,j)+1/scpx(i  ,j))
                cff1=0.5*(1/scpx(i  ,j)+1/scpx(i+1,j))
                dTdx=0.5 *                                      &
                &       (cff*(temp(i  ,j,k,m)-temp(i-1,j,k,m)) +  &
                &       cff1*(temp(i+1,j,k,m)-temp(i  ,j,k,m)))
          
                dSdx=0.5 *                                      &
                &       (cff*(saln(i  ,j,k,m)-saln(i-1,j,k,m)) +  &
                &       cff1*(saln(i+1,j,k,m)-saln(i  ,j,k,m)))  
                dBdx=-0.5 * gorho0 *                                 &
                &         (cff*(th3d(i  ,j,k,m)-th3d(i-1,j,k,m)) +            &   
                &         cff1*(th3d(i+1,j,k,m)-th3d(i  ,j,k,m)))    
                cff =0.5*(1/scpy(i,j-1)+1/scpy(i,j  ))
                cff1=0.5*(1/scpy(i,j  )+1/scpy(i,j+1))
                dTde=0.5 *                                      &
                &       (cff*(temp(i,j  ,k,m)-temp(i,j-1,k,m)) +  &
                &       cff1*(temp(i,j+1,k,m)-temp(i,j  ,k,m)))  
                dSde=0.5 *                                      &
                &       (cff*(saln(i,j  ,k,m)-saln(i,j-1,k,m)) +  &
                &       cff1*(saln(i,j+1,k,m)-saln(i,j  ,k,m)))
                dBde=-0.5 * gorho0 *                                 &
                &       (cff*(th3d(i,j  ,k,m)-th3d(i,j-1,k,m))+               &
                &       cff1*(th3d(i,j+1,k,m)-th3d(i,j  ,k,m)))
              endif
              if (.not.( ip(i,j).ne.0 .and. ip(i+1,j).ne.0 .and. &
                 ip(i-1,j).ne.0 .and. ip(i,j+1).ne.0 .and. &
                 ip(i,j-1).ne.0 )) then
                dBdx=0.0
                dBde=0.0
                dTdx=0.0
                dTde=0.0
                dSdx=0.0
                dSde=0.0
              endif
              dTdx_z(i,j,k) = dTdx
              dTde_z(i,j,k) = dTde
              dSdx_z(i,j,k) = dSdx
              dSde_z(i,j,k) = dSde
              dBdx_z(i,j,k) = dBdx
              dBde_z(i,j,k) = dBde
            enddo
          enddo
        enddo
!$OMP END PARALLEL DO
!---------------------------------------------------------------------
!                                                                     
! Zhang23 Paratermeterlization of Submesoscale Vertical Bouyancy Flux                              
!                                                                     
!---------------------------------------------------------------------
        do j=1-margin,jj+margin
          do i=1-margin,ii+margin
            dp_m(i,j,1) = 0.0
              do k = 1,kk
                dp_m(i,j,k+1)=dp_m(i,j,k)+dp(i,j,k,m)/9806.0;   
                z_p(i,j,k) = 0.5 * ( dp_m(i,j,k) + dp_m(i,j,k+1) )
              enddo
          enddo
        enddo

!$OMP PARALLEL DO PRIVATE(j,i,k,S_sub,Hz_sum,alpha_bar,t_bar,s_bar, &
!$OMP                     twc_hbar,swc_hbar,mld_m,twc_tmp,swc_tmp, &
!$OMP                     z_H,dz,uz,twc_nodz,swc_nodz) &
!$OMP          SCHEDULE(STATIC,jblk) 
        do j=1-margin,jj+margin
          do i=1-margin,ii+margin
            if (SEA_P) then
              S_sub = 0.5*((scpx(i, j)) + (scpy(i, j)))  !!!deltx+delty
              Hz_sum=0.0
              alpha_bar=0.0
              t_bar=0.0
              s_bar=0.0
              twc_hbar=0.0
              swc_hbar=0.0
              mld_m=0.0
              do k = 1, kk  
                if (th3d(i,j,k,m) - th3d(i,j,1,m) .ge. 0.03) then
                  mld_m = z_p(i,j,k)+(z_p(i,j,k+1)-z_p(i,j,k)) /  &  
                  &         (th3d(i,j,k+1,m)-th3d(i,j,k,m)) *                 &
                  &         (th3d(i,j,1,m)+0.03-th3d(i,j,k,m))
                  exit 
                else
                  Hz_sum=Hz_sum + dp(i,j,k,m)/9806.0
                  alpha_bar=alpha_bar+dp(i,j,k,m)/9806.0*alpha_rho(i,j,k)
                  t_bar = t_bar +                 &
                    &               (dTde_z(i,j,k)*dBde_z(i,j,k) +   &
                    &               dTdx_z(i,j,k)*dBdx_z(i,j,k)) *   &
                    &               dp(i,j,k,m)/9806.0
                  s_bar = s_bar +                 &
                    &               (dSde_z(i,j,k)*dBde_z(i,j,k) +   &
                    &               dSdx_z(i,j,k)*dBdx_z(i,j,k)) *   &
                    &               dp(i,j,k,m)/9806.0
                endif
              enddo
              t_bar = t_bar/Hz_sum
              s_bar = s_bar/Hz_sum
              alpha_bar = alpha_bar/Hz_sum
              if ( (mld_m .gt. 0) .and. (dp_m(i,j,kk+1) .ge. 50) ) then
                twc_hbar = C_sub * B_sub * S_sub *             &
                &                 mld_m**2* t_bar *       &
                &                 alpha_bar**2/max(f_min, abs(corio(i,j)))**2
                swc_hbar = C_sub * B_sub * S_sub *             &
                &                 mld_m**2* s_bar *       &
                &                 alpha_bar**2/max(f_min, abs(corio(i,j)))**2
              endif
              if(twc_hbar .ge. 0) then
                twc_hbar = min(5e-4, twc_hbar)
              else
                twc_hbar = max(-5e-4, twc_hbar)
              endif
              if(swc_hbar .ge. 0) then
                swc_hbar = min(5e-4, swc_hbar)
              else
                swc_hbar = max(-5e-4, swc_hbar)
              endif
              do k = 1,kk
                twc_tmp=0.0
                swc_tmp=0.0
                twc_nodz=0.0
                swc_nodz=0.0
                if (mld_m .gt. 0) then
                  z_H = (z_p(i,j,k)/mld_m)**2
                  dz=1.0/mld_m*(64.0/27.0*z_H-4.0/3.0) *       &
                  &     exp(-8.0/9.0*z_H+1.0/2.0)
                  uz=4.0/3.0*sqrt(z_H)*exp(-8.0/9.0*z_H+1.0/2.0)
                endif
!
! Calculate Submesoscale Vertical Bouyancy Flux at Each Level
                twc_tmp = dz*twc_hbar
                swc_tmp = dz*swc_hbar
                twc_nodz= uz*twc_hbar
                swc_nodz= uz*swc_hbar
                temp(i,j,k,n) = temp(i,j,k,n) - delt1 * twc_tmp
                saln(i,j,k,n) = saln(i,j,k,n) - delt1 * swc_tmp
                th3d(i,j,k,n) = sig(temp(i,j,k,n),saln(i,j,k,n))-thbase
                VBF_para(i,j,k)=9.806*(1.7e-4*twc_nodz - 7.6e-4*swc_nodz)
                VHF_para(i,j,k)=1030.0*3900.0*twc_nodz
              enddo
            endif
          enddo
        enddo
!$OMP END PARALLEL DO
        ktr = itracr(804)
        if     (ktr.ne.0) then
!$OMP   PARALLEL DO PRIVATE(j,i,k) &
!$OMP            SCHEDULE(STATIC,jblk)
          do j=1-margin,jj+margin
            do i=1-margin,ii+margin
              do k = 1, kk 
                tracer(i,j,k,n,ktr)=VBF_para(i,j,k)
              enddo
            enddo
          enddo
!$OMP   END PARALLEL DO
        endif !itracr(804)
!!!!!!!
        htr = itracr(805)
        if     (htr.ne.0) then
!$OMP   PARALLEL DO PRIVATE(j,i,k) &
!$OMP            SCHEDULE(STATIC,jblk)
          do j=1-margin,jj+margin
            do i=1-margin,ii+margin
              do k = 1, kk 
                tracer(i,j,k,n,htr)=VHF_para(i,j,k)
              enddo
            enddo
          enddo
!$OMP   END PARALLEL DO
        endif !itracr(805)		
      return
      end subroutine z23_vbf
    end module mod_z23vbf

