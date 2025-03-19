! #######################   Main Subroutines   #######################
Subroutine Modes(ng,nk,nb,nx,ny,nz,nmat,order,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,chi,&
                 nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,x_smax,x_smin,xdiv,ydiv,zdiv,&
                 y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_bott,z_top,ntem,stab,&
				 cf,tin,cftem,cmtem,ccden,bpos,nfpin,cmflow,rf,tg,tc,ppitch,&
				 csiga,csigtr,csigf,cnuf,nstep,xsize,ysize,bmap,pra,&
				 fsiga,fsigtr,fsigf,fnuf,msiga,msigtr,msigf,mnuf,lsiga,lsigtr,lsigf,lnuf,&
				 dsiga,dsigtr,dsigf,dnuf,bcon,rbcon,rftem,rmtem,rcden,csigs,fsigs,msigs,&
				 lsigs,dsigs,powr,ppow,pos0,ssize,mode,&
				 fbpos,tmove,bspeed,ibeta,lamb,velo,ttot,tstep1,tstep2,tdiv)
Implicit None
Integer, Parameter :: dp = Selected_Real_Kind(10, 15)
! #######################   User Input Paramaters   #######################
Integer, Intent(in) :: ng,pra                            ! Number of Groups
Integer, Intent(in) :: nk                             ! Number of Nodes
Integer, Intent(in) :: order                          ! Order Polynomial NEM
Integer, Intent(in) :: nx, ny, nz                     ! Number of Assemblies in x, y & z Directions
Integer, Intent(in) :: nxx, nyy, nzz                  ! Number of Nodes in x, y, & z Directions
Integer, Intent(in) :: nmat                           ! Number of Materials
Integer, Dimension(nx), Intent(in) :: xdiv            ! Assembly Division
Integer, Dimension(ny), Intent(in) :: ydiv            ! Assembly Division
Integer, Dimension(nz), Intent(in) :: zdiv            ! Assembly Division
Integer, Dimension(nk), Intent(in) :: mat             ! Materials
Real(dp), Dimension(nx), Intent(in) :: xsize           ! Assembly Size & Division
Real(dp), Dimension(ny), Intent(in) :: ysize           ! Assembly Size & Division
Integer, Dimension(nxx), Intent(in) :: x_smax, x_smin ! imax and imin along x Direction for Staggered Nodes
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin ! imax and imin along y Direction for Staggered Nodes
Integer, Intent(in) :: x_east, x_west, y_north        ! Boundary Conditions
Integer, Intent(in) :: y_south, z_top, z_bott         ! Boundary Conditions
Real(dp), Dimension(nxx), Intent(in) :: delx          ! Delta_x Volume in cm3
Real(dp), Dimension(nyy), Intent(in) :: dely          ! Delta_y Volume in cm3
Real(dp), Dimension(nzz), Intent(in) :: delz          ! Delta_z Volume in cm3
Real(dp), Dimension(nk), Intent(in) :: delv           ! Nod's Volume in cm3
Real(dp), Dimension(nk), Intent(in) :: ix, iy, iz     ! x-Position, y-Position and z-Position of the Node
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz   ! xyz-Position of the Node
Character(Len=100), intent(in) :: mode
! Fixed Source Paramaters
! Real(dp), Dimension(nk,ng), Intent(in) ::  eSrc         ! External Source
! Integer, Intent(in) :: nS                             ! Sources Number
! Real(dp), Dimension(nS), Intent(in) :: dS             ! Source Density
! Real(dp), Dimension(nS,ng), Intent(in) :: spS         ! Source Spectrum
! Integer, Intent(in) :: npor                           ! Number of Radial Fixed Sources Position
! Integer, Intent(in) :: npox                           ! Number of Axial Fixed Sources Position
! Integer, Dimension(nS,npox), Intent(in) :: zpos       ! Axial Position of Extra Sources
! Integer, Dimension(nS,npor), Intent(in) :: xpos, ypos ! Radial Position of Extra Sources
! CXs Assigned to Materials
Real(dp), Dimension(nmat,ng), Intent(INOUT) :: chi       ! Neutron Fission Spectrum
Real(dp), Dimension(nmat,ng), Intent(INOUT) :: x_siga    ! Absorption Macroscopic CX
Real(dp), Dimension(nmat,ng), Intent(INOUT) :: x_sigtr   ! Transport Macroscopic CX
Real(dp), Dimension(nmat,ng), Intent(INOUT) :: x_sigf    ! Fission Macroscopic CX
Real(dp), Dimension(nmat,ng), Intent(INOUT) :: x_nu_sigf ! nu * Fission Macroscopic CX
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: x_sigs ! Scattering Macroscopic CX
! ##########################################################################

! Transient parameters
! Crod changes
Integer, Intent(in) :: nb                             ! Number of CR banks
REAL(DP), Intent(in) :: nstep		! Number of steps
REAL(DP), INTENT(IN) :: pos0, ssize		! Zero step position and step size
INTEGER, PARAMETER :: nf = 6			! Number of delaye dneutron precusor family
Real(dp), Dimension(nb), intent(in) :: bpos
INTEGER, DIMENSION(nx, ny), intent(in) :: bmap       ! Radial control rod bank map (assembly wise)                                      ! WHEN SECOND TIME STEP APPLY
INTEGER,  INTENT(IN) :: ntem    ! Number of temperature in steam table
REAL(DP), DIMENSION(ntem,pra), intent(in) :: stab  ! Steam table matrix
!ROD EJECTION DATA
REAL(DP), DIMENSION(nb), intent(in) :: fbpos    ! Final CR bank position
REAL(DP), DIMENSION(nb), intent(in) :: tmove    ! Time when CR bank starts moving
REAL(DP), DIMENSION(nb), intent(in) :: bspeed   ! CR bank movement speed
REAL(DP), DIMENSION(nf), intent(in) :: ibeta, lamb	! beta (delayed neutron fraction) and precusor decay constant
REAL(DP), DIMENSION(ng), intent(in) :: velo         ! Neutron velocity
REAL(DP), intent(in) :: ttot						! TOTAL SIMULATION TIME
REAL(DP), intent(in) :: tstep1						! FIRST TIME STEP
REAL(DP), intent(in) :: tstep2						! SECOND TIME STEP
REAL(DP), intent(in) :: tdiv 

! Thermal-hydraulics parameters
REAL(DP), intent(in) :: cf		! heat fraction deposited into coolant
REAL(DP), INTENT(IN) :: tin		! coolant inlet temperature (kelvin)
Real(dp), intent(in) :: cftem, cmtem, ccden
INTEGER, intent(in) :: nfpin	! Number of fuel pin and guide tubes
REAL(DP), intent(in) :: cmflow
REAL(DP), intent(in) :: rf, tg, tc, ppitch	! Fuel meat radius, gap thickness, clad thickness, and pin picth (m)
INTEGER, PARAMETER :: nm = 10      ! Fuel meat divided into 10 mesh
INTEGER, PARAMETER :: nt = nm + 2      ! two more mesh for gap and clad
REAL(DP), DIMENSION(nt) :: rdel		! mesh delta
REAL(DP), DIMENSION(nt) :: rpos		! mesh position

Real(dp), Dimension(nmat,ng), Intent(INOUT) :: csiga, csigtr, csigf, cnuf,&
fsiga, fsigtr, fsigf, fnuf, msiga, msigtr, msigf, mnuf, &
lsiga, lsigtr, lsigf, lnuf, dsiga, dsigtr, dsigf, dnuf
Real(dp), Intent(in) :: bcon, rbcon ,rftem, rmtem, rcden
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: csigs, fsigs, msigs, lsigs, dsigs
REAL(DP), Intent(in):: powr, ppow

!##########################################################################

Select Case(mode)
    ! Case('Fixed Source')
        ! Call FixedSrc(ng,nk,nx,ny,nz,nmat,order,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,chi,&
                      ! nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,x_smax,x_smin,y_smax,y_smin,xdiv,ydiv,zdiv,&
                      ! xyz,x_east,x_west,y_north,y_south,z_top,z_bott,mode,eSrc)
 Case('Adjoint')
     Call Adjoint(ng,nk,nx,ny,nz,nmat,order,mat,x_sigs,x_siga,x_sigtr,x_sigf,&
     		  x_nu_sigf,chi,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,&
     		  x_smax,x_smin,xdiv,ydiv,zdiv,y_smax,y_smin,xyz,x_east,x_west,&
     		  y_north,y_south,z_top,z_bott,mode)
 Case('Forward')
     Call Forward(ng,nk,nx,ny,nz,nmat,order,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,chi,&
                  nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,x_smax,x_smin,xdiv,ydiv,zdiv,&
                  y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott,mode)
 Case('CBC')
		Call cbsearch(ng,nk,nb,nx,ny,nz,nmat,order,mat,x_sigs,x_siga,&
			      x_sigtr,x_sigf,x_nu_sigf,chi,nxx,nyy,nzz,ix,iy,iz,&
			      delx,dely,delz,delv,x_smax,x_smin,xdiv,ydiv,zdiv,&
			      y_smax,y_smin,xyz,x_east,x_west,y_north,&
			      y_south,z_top,z_bott,mode,rbcon,csiga,csigtr,csigf,&
			      cnuf,csigs,cftem,rftem,fsiga,fsigtr,fsigf,fnuf,fsigs,&
			      msiga,msigtr,msigf,mnuf,msigs,rmtem,cmtem,lsiga,lsigtr,&
			      lsigf,lnuf,lsigs,rcden,ccden,dsiga,dsigtr,dsigf,dnuf,&
			      dsigs,pos0,ssize,bpos,bmap,nstep)!    To search critical boron concentration
	Case('CBCTH')
		Call cbsearchth(ng,nk,nb,nx,ny,nz,nmat,order,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,chi,&
					   nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,x_smax,x_smin,xdiv,ydiv,zdiv,xsize,ysize,&
					   y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott,mode,rbcon,&
					   csiga,csigtr,csigf,cnuf,csigs,cftem,rftem,fsiga,fsigtr,fsigf,fnuf,fsigs,&
					   msiga,msigtr,msigf,mnuf,msigs,rmtem,cmtem,&
					   lsiga,lsigtr,lsigf,lnuf,lsigs,rcden,ccden,dsiga,dsigtr,dsigf,dnuf,dsigs,&
					   pos0,ssize,bpos,bmap,nstep,&
					   ntem,stab,cf,tin,nfpin,cmflow,rf,tg,tc,ppitch,pra,powr,ppow)
   Case('RODEJCT')
		!Call rod_eject(ng,nk,nb,nx,ny,nz,nmat,order,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,chi,&
		!			   nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,x_smax,x_smin,xdiv,ydiv,zdiv,xsize,ysize,&
		!			   y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott,mode,&
		!			   dsiga,dsigtr,dsigf,dnuf,dsigs,pos0,ssize,bpos,bmap,nstep,fbpos,tmove,bspeed,ibeta,&
		!			   lamb,velo,ttot,tstep1,tstep2,tdiv)
		Call rod_eject(ng,nk,nb,nx,ny,nz,nmat,order,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,chi,&
			       nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,x_smax,x_smin,xdiv,ydiv,zdiv,xsize,ysize,&
			       y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott,mode,&
			       dsiga,dsigtr,dsigf,dnuf,dsigs,pos0,ssize,bpos,bmap,nstep,fbpos,tmove,bspeed,ibeta,&
			       lamb,velo,ttot,tstep1,tstep2,tdiv)
  Case('THRODEJCT')
		Call throd_eject(ng,nk,nb,nx,ny,nz,nmat,order,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,chi,&
						 nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,x_smax,x_smin,xdiv,ydiv,zdiv,xsize,ysize,&
						 y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott,mode,&
						 dsiga,dsigtr,dsigf,dnuf,pos0,ssize,bpos,bmap,nstep,fbpos,tmove,bspeed,ibeta,&
						 lamb,velo,ttot,tstep1,tstep2,tdiv,&
						 ntem,stab,pra,cf,tin,cftem,cmtem,ccden,nfpin,cmflow,rf,tg,tc,ppitch,&
						 csiga,csigtr,csigf,cnuf,fsiga,fsigtr,fsigf,fnuf,msiga,msigtr,msigf,mnuf,&
						 lsiga,lsigtr,lsigf,lnuf,bcon,rbcon,rftem,rmtem,rcden,csigs,&
						 fsigs,msigs,lsigs,dsigs,powr,ppow)
   !Case('FTEMP')
	!	Call ftemtest(ng,nk,nb,nx,ny,nz,nmat,order,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,chi,&
	!					  nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,x_smax,x_smin,xdiv,ydiv,zdiv,xsize,ysize,&
	!				      y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott,mode,&
	!				      bpos,ibeta,lamb,velo,ttot,tstep1,tstep2,tdiv,fsigs,fsiga,fsigtr,fsigf,fnuf,rftem,cftem,&
	!					  dsiga,dsigtr,dsigf,dnuf,dsigs,bmap,pos0,ssize,nstep,fbpos,bspeed,tmove,& 
	!					  cmtem,ccden,&
	!					  csiga,csigtr,csigf,cnuf,msiga,msigtr,msigf,mnuf,&
	!					  lsiga,lsigtr,lsigf,lnuf,bcon,rbcon,rmtem,rcden,csigs,&
	!					  msigs,lsigs) !    To perform rod ejection simulation
End Select

End Subroutine

    Subroutine Forward(ng,nk,nx,ny,nz,nmat,order,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,chi,&
                       nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,x_smax,x_smin,xdiv,ydiv,zdiv,&
                       y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott,mode) ! To Solve Forward (Normal) Problems
Implicit None
Integer, Parameter :: dp = Selected_Real_Kind(10,15)
! #######################   User Input Paramaters   #######################
Integer, Intent(in) :: ng, nk, order, nx, ny, nz, nmat, nxx, nyy, nzz, x_east, x_west, y_north, y_south, z_top, z_bott
Integer, Dimension(nxx), Intent(in) :: x_smax, x_smin
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin
Integer, Dimension(nk), Intent(in) :: mat
Integer, Dimension(nx), Intent(in) :: xdiv
Integer, Dimension(ny), Intent(in) :: ydiv
Integer, Dimension(nz), Intent(in) :: zdiv
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk), Intent(in) :: delv, ix, iy, iz
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz
Character(Len=100), intent(in) :: mode
! CXs Assigned to Materials
Real(dp), Dimension(nmat,ng), Intent(in) :: chi, x_siga, x_sigtr, x_sigf, x_nu_sigf
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: x_sigs

!##########################################################################
! Local Variables
Real(dp) :: Keff   ! Effectif Multiplication Factor
Real(dp) :: fer, fser	! Flux and Fission Source Error in BCSEARCH calcs.
Real(dp), Dimension(nk,ng) :: f0, fx1, fy1, fz1, fx2, fy2, fz2      ! Flux and Flux Moments
Real(dp), Dimension(nk,ng,6) :: jo   ! Nodals' Outgoing Currents (X+,X-,Y+,Y-,Z+,Z-)
Real(dp), Dimension(nk,ng,6) :: ji   ! Nodals' Ingoing Currents  (X+,X-,Y+,Y-,Z+,Z-)
Real(dp), Dimension(nk,ng,3) :: L0   ! Zeroth Transverse Leakages (Lx,Ly,Lz)
Real(dp), Dimension(nk,ng,6) :: al   ! Assembly Discontinuity Factor
Real(dp), Dimension(nk,ng,6,6) :: R4, R2, P2    ! Response Matrix
Real(dp), Dimension(nk,ng,6,7) :: P4            ! Response Matrix
! Real(dp), Dimension(nk,ng,7) :: Q      ! Nodal's Source and Source Moments (0,x1,y1,z1,x2,y2,z2)
! CXs Assigned to Nodes
Real(dp), Dimension(nk,ng,ng) :: sigs  ! Scattering Macroscopic CX
Real(dp), Dimension(nk,ng) :: siga     ! Absorption Macroscopic CX
Real(dp), Dimension(nk,ng) :: sigtr    ! Transport Macroscopic CX
Real(dp), Dimension(nk,ng) :: sigf     ! Fission Macroscopic CX
Real(dp), Dimension(nk,ng) :: nu_sigf  ! nu * Fission Macroscopic CX
Real(dp), Dimension(nk,ng) :: D        ! Diffusion Coefficient
Real(dp), Dimension(nk,ng) :: sigr     ! Scattering Macroscopic CX
Real(dp), Dimension(nk) :: power
Real(kind=4) :: st, fn
Integer, Parameter :: Ounit = 100      ! Output

Open (Unit=ounit, File='app/Output/NEM.out')
Write(ounit,*)
Write(ounit,*)
Write(ounit,*) ' ==============================================' &
// '=================================================='
Write(ounit,*) &
'                                       Calculation Results'
Write(ounit,*) ' ==============================================' &
// '=================================================='
Write(ounit,*)
Write(ounit,*) '  Iteration    Keff     FissSrc Error   Inner Error'
Write(ounit,*) ' ----------------------------------------------------'
Write(*,*)
Write(*,*)
Write(*,*) ' ==============================================' &
// '===================='
Write(*,*) &
'                         Calculation Results'
Write(*,*) ' ==============================================' &
// '===================='
Write(*,*)
Write(*,*) '  Iteration    Keff     FissSrc Error   Inner Error'
Write(*,*) ' ---------------------------------------------------------'

Call CPU_TIME(st)

! #######################   Initialiaze Guess Values   #######################
    Call Init(keff,jo,ji,L0,al,f0,fx1,fy1,fz1,fx2,fy2,fz2,ng,nk)

! #######################   Update XSec   #######################
    Call XS_updt(sigs,siga,sigtr,sigf,nu_sigf,D,sigr,ng,nk,nmat,&
                 mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf)

! #######################   Update Diffusion Coefficient and Removal XS   #######################
    Call D_sigr(D,sigr,siga,sigs,sigtr,ng,nk)

! #######################   Calculate Nodal Coupling Matrices   #######################
    Call Nodal_Coup4(P4,R4,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
    Call Nodal_Coup2(P2,R2,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)

! #######################   Call Outer Iteration   #######################
    Call Outer(Keff,fer,fser,f0,fx1,fy1,fz1,fx2,fy2,fz2,ng,nk,order,R2,P2,R4,P4,nmat,mat,al,jo,ji,L0,&
               D,sigr,chi,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,xyz,&
               x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott)

! #######################   Flux & Power Distribution Output   #######################
    Call PowDis(power,f0,sigf,delv,ng,nk,mode)
! Radial Power Distribution
    Call AsmPow(nk,nx,ny,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,xdiv,ydiv,y_smax,y_smin,power)
! Axial Power Distribution
    Call AxiPow(nk,nz,nxx,nyy,nzz,ix,iy,iz,xyz,delz,delv,zdiv,y_smax,y_smin,power)
! Radial Flux Distribution
    Call AsmFlux(ng,nk,nx,ny,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,xdiv,ydiv,y_smax,y_smin,f0,1.e0_DP)
! Axial Flux Distribution
    ! Call AxiFlux(ng,nk,nz,nxx,nyy,nzz,ix,iy,iz,xyz,delz,delv,zdiv,y_smax,y_smin,f0) 

Call CPU_TIME(fn)
Write(ounit,*)
Write(*,*)
Write(ounit,*) "Total Time : ", fn-st, "Seconds"
Write(*,*) "Total Time : ", fn-st, "Seconds"

End Subroutine

    Subroutine Adjoint(ng,nk,nx,ny,nz,nmat,order,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,chi,&
                       nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,x_smax,x_smin,xdiv,ydiv,zdiv,&
                       y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott,mode) ! To Solve Adjoint Problems
Implicit None
Integer, Parameter :: dp = Selected_Real_Kind(10,15)
! #######################   User Input Paramaters   #######################
Integer, Intent(in) :: ng, nk, order, nx, ny, nz, nmat, nxx, nyy, nzz, x_east, x_west, y_north, y_south, z_top, z_bott
Integer, Dimension(nxx), Intent(in) :: x_smax, x_smin
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin
Integer, Dimension(nk), Intent(in) :: mat
Integer, Dimension(nx), Intent(in) :: xdiv
Integer, Dimension(ny), Intent(in) :: ydiv
Integer, Dimension(nz), Intent(in) :: zdiv
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk), Intent(in) :: delv, ix, iy, iz
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz
Character(Len=100), intent(in) :: mode
! CXs Assigned to Materials
Real(dp), Dimension(nmat,ng), Intent(in) :: chi, x_siga, x_sigtr, x_sigf, x_nu_sigf
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: x_sigs
!##########################################################################
! Local Variables
Real(dp) :: Keff   ! Effectif Multiplication Factor
Real(dp), Dimension(nk,ng) :: f0, fx1, fy1, fz1, fx2, fy2, fz2      ! Flux and Flux Moments
Real(dp), Dimension(nk,ng,6) :: jo   ! Nodals' Outgoing Currents (X+,X-,Y+,Y-,Z+,Z-)
Real(dp), Dimension(nk,ng,6) :: ji   ! Nodals' Ingoing Currents  (X+,X-,Y+,Y-,Z+,Z-)
Real(dp), Dimension(nk,ng,3) :: L0   ! Zeroth Transverse Leakages (Lx,Ly,Lz)
Real(dp), Dimension(nk,ng,6) :: al   ! Assembly Discontinuity Factor
Real(dp), Dimension(nk,ng,6,6) :: R4, R2, P2    ! Response Matrix
Real(dp), Dimension(nk,ng,6,7) :: P4            ! Response Matrix
! Real(dp), Dimension(nk,ng,7) :: Q      ! Nodal's Source and Source Moments (0,x1,y1,z1,x2,y2,z2)
! CXs Assigned to Nodes
Real(dp), Dimension(nk,ng,ng) :: sigs  ! Scattering Macroscopic CX
Real(dp), Dimension(nk,ng) :: siga     ! Absorption Macroscopic CX
Real(dp), Dimension(nk,ng) :: sigtr    ! Transport Macroscopic CX
Real(dp), Dimension(nk,ng) :: sigf     ! Fission Macroscopic CX
Real(dp), Dimension(nk,ng) :: nu_sigf  ! nu * Fission Macroscopic CX
Real(dp), Dimension(nk,ng) :: D        ! Diffusion Coefficient
Real(dp), Dimension(nk,ng) :: sigr     ! Scattering Macroscopic CX
Real(dp), Dimension(nk) :: power
Real(dp) :: st, fn
Integer, Parameter :: Ounit = 100      ! Output

Open (Unit=ounit, File='app/Output/NEM.out')
Write(ounit,*)
Write(ounit,*)
Write(ounit,*) ' ==============================================' &
// '=================================================='
Write(ounit,*) &
'                                       Calculation Results'
Write(ounit,*) ' ==============================================' &
// '=================================================='
Write(ounit,*)
Write(ounit,*) '  Iteration    Keff     FissSrc Error   Inner Error'
Write(ounit,*) ' ----------------------------------------------------'
Write(*,*)
Write(*,*)
Write(*,*) ' ==============================================' &
// '===================='
Write(*,*) &
'                         Calculation Results'
Write(*,*) ' ==============================================' &
// '===================='
Write(*,*)
Write(*,*) '  Iteration    Keff     FissSrc Error   Inner Error'
Write(*,*) ' ---------------------------------------------------------'

Call CPU_TIME(st)

! #######################   Initialiaze Guess Values   #######################
    Call Init(keff,jo,ji,L0,al,f0,fx1,fy1,fz1,fx2,fy2,fz2,ng,nk)

! #######################   Update XSec   #######################
    Call XS_updt(sigs,siga,sigtr,sigf,nu_sigf,D,sigr,ng,nk,nmat,&
                 mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf)

! #######################   Update Diffusion Coefficient and Removal XS   #######################
    Call D_sigr(D,sigr,siga,sigs,sigtr,ng,nk)

! #######################   Calculate Nodal Coupling Matrices   #######################
    Call Nodal_Coup4(P4,R4,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
    Call Nodal_Coup2(P2,R2,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)

! #######################   Call Outer Iteration   #######################
    Call OuterAdj(Keff,f0,fx1,fy1,fz1,fx2,fy2,fz2,ng,nk,order,R2,P2,R4,P4,nmat,chi,mat,al,jo,ji,L0,&
                  D,sigr,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,&
                  x_smax,x_smin,y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott)

! #######################   Flux & Power Distribution Output   #######################
    Call PowDis(power,f0,sigf,delv,ng,nk,mode)
! Radial Power Distribution
    Call AsmPow(nk,nx,ny,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,xdiv,ydiv,y_smax,y_smin,power)
! Axial Power Distribution
    Call AxiPow(nk,nz,nxx,nyy,nzz,ix,iy,iz,xyz,delz,delv,zdiv,y_smax,y_smin,power)
! Radial Flux Distribution
    Call AsmFlux(ng,nk,nx,ny,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,xdiv,ydiv,y_smax,y_smin,f0,1.e0_DP)
Call CPU_TIME(fn)
Write(ounit,*) "Total Time : ", fn-st, " Seconds"

End Subroutine

    Subroutine FixedSrc(ng,nk,nx,ny,nz,nmat,order,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,chi,&
                        nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,x_smax,x_smin,y_smax,y_smin,xdiv,ydiv,zdiv,&
                        xyz,x_east,x_west,y_north,y_south,z_top,z_bott,mode,eSrc) ! To Solve Fixed Sources Problems
Implicit None
Integer, Parameter :: dp = Selected_Real_Kind(10,15)
! #######################   User Input Paramaters   #######################
Integer, Intent(in) :: ng, nk, order, nx, ny, nz, nmat, nxx, nyy, nzz, x_east, x_west, y_north, y_south, z_top, z_bott
Integer, Dimension(nxx), Intent(in) :: x_smax, x_smin
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin
Integer, Dimension(nk), Intent(in) :: mat
Integer, Dimension(nx), Intent(in) :: xdiv
Integer, Dimension(ny), Intent(in) :: ydiv
Integer, Dimension(nz), Intent(in) :: zdiv
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk), Intent(in) :: delv, ix, iy, iz
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz
Real(dp), Dimension(nk,ng), Intent(in) :: eSrc
! Integer, Intent(in) :: nS
! Real(dp), Dimension(nS), Intent(in) :: dS
! Real(dp), Dimension(nS,ng), Intent(in) :: spS
! Integer, Dimension(nS,npox), Intent(in) :: zpos
! Integer, Dimension(nS,npor), Intent(in) :: xpos, ypos
Character(Len=100), intent(in) :: mode
! CXs Assigned to Materials
Real(dp), Dimension(nmat,ng), Intent(in) :: chi, x_siga, x_sigtr, x_sigf, x_nu_sigf
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: x_sigs
!##########################################################################
! Local Variables
Real(dp) :: Keff   ! Effectif Multiplication Factor
Real(dp), Dimension(nk,ng) :: f0, fx1, fy1, fz1, fx2, fy2, fz2      ! Flux and Flux Moments
Real(dp), Dimension(nk,ng,6) :: jo   ! Nodals' Outgoing Currents (X+,X-,Y+,Y-,Z+,Z-)
Real(dp), Dimension(nk,ng,6) :: ji   ! Nodals' Ingoing Currents  (X+,X-,Y+,Y-,Z+,Z-)
Real(dp), Dimension(nk,ng,3) :: L0   ! Zeroth Transverse Leakages (Lx,Ly,Lz)
Real(dp), Dimension(nk,ng,6) :: al   ! Assembly Discontinuity Factor
Real(dp), Dimension(nk,ng,6,6) :: R4, R2, P2    ! Response Matrix
Real(dp), Dimension(nk,ng,6,7) :: P4            ! Response Matrix
! Real(dp), Dimension(nk,ng,7) :: Q      ! Nodal's Source and Source Moments (0,x1,y1,z1,x2,y2,z2)
! CXs Assigned to Nodes
Real(dp), Dimension(nk,ng,ng) :: sigs  ! Scattering Macroscopic CX
Real(dp), Dimension(nk,ng) :: siga     ! Absorption Macroscopic CX
Real(dp), Dimension(nk,ng) :: sigtr    ! Transport Macroscopic CX
Real(dp), Dimension(nk,ng) :: sigf     ! Fission Macroscopic CX
Real(dp), Dimension(nk,ng) :: nu_sigf  ! nu * Fission Macroscopic CX
Real(dp), Dimension(nk,ng) :: D        ! Diffusion Coefficient
Real(dp), Dimension(nk,ng) :: sigr     ! Scattering Macroscopic CX
Real(dp), Dimension(nk) :: power
Real(dp) :: st, fn
Integer, Parameter :: Ounit = 100      ! Output

Open (Unit=ounit, File='app/Output/NEM.out')

! #######################   Print Extra Sources   #######################
    ! Call ExSrc(eSrc,nk,ng,nx,ny,nz,nxx,nyy,nzz,nS,dS,spS,&
               ! xyz,xdiv,ydiv,zdiv,xpos,ypos,zpos,npox,npor)

Write(ounit,*)
Write(ounit,*)
Write(ounit,*) ' ==============================================' &
// '=================================================='
Write(ounit,*) &
'                                       Calculation Results'
Write(ounit,*) ' ==============================================' &
// '=================================================='
Write(ounit,*)
Write(ounit,*) '  Iteration    Keff     FissSrc Error   Inner Error'
Write(ounit,*) ' ----------------------------------------------------'
Write(*,*)
Write(*,*)
Write(*,*) ' ==============================================' &
// '===================='
Write(*,*) &
'                         Calculation Results'
Write(*,*) ' ==============================================' &
// '===================='
Write(*,*)
Write(*,*) '  Iteration    Keff     FissSrc Error   Inner Error'
Write(*,*) ' ---------------------------------------------------------'

Call CPU_TIME(st)

! #######################   Initialiaze Guess Values   #######################
    Call Init(keff,jo,ji,L0,al,f0,fx1,fy1,fz1,fx2,fy2,fz2,ng,nk)

! #######################   Update XSec   #######################
    Call XS_updt(sigs,siga,sigtr,sigf,nu_sigf,D,sigr,ng,nk,nmat,&
                 mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf)

! #######################   Update Diffusion Coefficient and Removal XS   #######################
    Call D_sigr(D,sigr,siga,sigs,sigtr,ng,nk)

! #######################   Calculate Nodal Coupling Matrices   #######################
    Call Nodal_Coup4(P4,R4,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
    Call Nodal_Coup2(P2,R2,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)

! #######################   Call Outer Iteration   #######################
    Call OuterFx(f0,fx1,fy1,fz1,fx2,fy2,fz2,ng,nk,order,R2,P2,R4,P4,nmat,chi,mat,al,jo,ji,L0,&
                 D,sigr,sigs,siga,nu_sigf,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,eSrc,&
                 x_smax,x_smin,y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott)

! #######################   Flux & Power Distribution Output   #######################
    Call PowDis(power,f0,sigf,delv,ng,nk,mode)
! Radial Power Distribution
    Call AsmPow(nk,nx,ny,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,xdiv,ydiv,y_smax,y_smin,power)
! Axial Power Distribution
    Call AxiPow(nk,nz,nxx,nyy,nzz,ix,iy,iz,xyz,delz,delv,zdiv,y_smax,y_smin,power)
! Radial Flux Distribution
    Call AsmFlux(ng,nk,nx,ny,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,xdiv,ydiv,y_smax,y_smin,f0,1.e0_DP)
    
Call CPU_TIME(fn)
Write(ounit,*) "Total Time : ", fn-st, " Seconds"

End Subroutine

	SUBROUTINE cbsearch(ng,nk,nb,nx,ny,nz,nmat,order,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,chi,&
					    nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,x_smax,x_smin,xdiv,ydiv,zdiv,&
					    y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott,mode,rbcon,&
						csiga,csigtr,csigf,cnuf,csigs,cftem,rftem,fsiga,fsigtr,fsigf,fnuf,fsigs,&
						msiga,msigtr,msigf,mnuf,msigs,rmtem,cmtem,&
						lsiga,lsigtr,lsigf,lnuf,lsigs,rcden,ccden,dsiga,dsigtr,dsigf,dnuf,dsigs,&
						pos0,ssize,bpos,bmap,nstep)!    To search critical boron concentration

Implicit None
Integer, Parameter :: dp = Selected_Real_Kind(10,15)
! #######################   User Input Paramaters   #######################
Integer, Intent(in) :: ng, nk,nb, order, nx, ny, nz, nmat, nxx, nyy, nzz, x_east, x_west, y_north, y_south, z_top, z_bott
Integer, Dimension(nxx), Intent(in) :: x_smax, x_smin
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin
Integer, Dimension(nk), Intent(in) :: mat
Integer, Dimension(nx), Intent(in) :: xdiv
Integer, Dimension(ny), Intent(in) :: ydiv
Integer, Dimension(nz), Intent(in) :: zdiv
!Real(dp), Dimension(nx), Intent(in) :: xsize           ! Assembly Size & Division
!Real(dp), Dimension(ny), Intent(in) :: ysize           ! Assembly Size & Division
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk), Intent(in) :: delv, ix, iy, iz
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz
Character(Len=100), intent(in) :: mode
Real(dp), Dimension(nmat,ng), Intent(in) :: chi, x_siga, x_sigtr, x_sigf, x_nu_sigf
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: x_sigs
Real(dp), Intent(in) :: rbcon
Real(dp), Dimension(nmat,ng), Intent(in) :: csiga, csigtr, csigf, cnuf,&
fsiga,fsigtr,fsigf,fnuf,msiga,msigtr,msigf,mnuf,lsiga,lsigtr,lsigf,lnuf,&
dsiga, dsigtr, dsigf, dnuf
REAL(DP), INTENT(IN) :: pos0, ssize			! Zero step position and step size
Real(dp), Dimension(nb), intent(in) :: bpos
INTEGER, DIMENSION(nx, ny), intent(in) :: bmap       ! Radial control rod bank map (assembly wise)
INTEGER, DIMENSION(nxx,nyy) :: fbmap	! Radial control rod bank map (node wise)
REAL(DP), Intent(in) :: nstep		! Number of steps
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: csigs,fsigs,msigs,lsigs,dsigs
! Local Variables
Real(dp) :: Keff   ! Effectif Multiplication Factor
Real(dp) :: fer, fser	! Flux and Fission Source Error in BCSEARCH calcs.
Real(dp), Dimension(nk,ng) :: f0, fx1, fy1, fz1, fx2, fy2, fz2      ! Flux and Flux Moments
Real(dp), Dimension(nk,ng,6) :: jo   ! Nodals' Outgoing Currents (X+,X-,Y+,Y-,Z+,Z-)
Real(dp), Dimension(nk,ng,6) :: ji   ! Nodals' Ingoing Currents  (X+,X-,Y+,Y-,Z+,Z-)
Real(dp), Dimension(nk,ng,3) :: L0   ! Zeroth Transverse Leakages (Lx,Ly,Lz)
Real(dp), Dimension(nk,ng,6) :: al   ! Assembly Discontinuity Factor
Real(dp), Dimension(nk,ng,6,6) :: R4, R2, P2    ! Response Matrix
Real(dp), Dimension(nk,ng,6,7) :: P4            ! Response Matrix
! Real(dp), Dimension(nk,ng,7) :: Q      ! Nodal's Source and Source Moments (0,x1,y1,z1,x2,y2,z2)
! CXs Assigned to Nodes
Real(dp), Dimension(nk,ng,ng) :: sigs  ! Scattering Macroscopic CX
Real(dp), Dimension(nk,ng) :: siga     ! Absorption Macroscopic CX
Real(dp), Dimension(nk,ng) :: sigtr    ! Transport Macroscopic CX
Real(dp), Dimension(nk,ng) :: sigf     ! Fission Macroscopic CX
Real(dp), Dimension(nk,ng) :: nu_sigf  ! nu * Fission Macroscopic CX
Real(dp), Dimension(nk,ng) :: D        ! Diffusion Coefficient
Real(dp), Dimension(nk,ng) :: sigr     ! Scattering Macroscopic CX
REAL(DP), DIMENSION(nk):: power	! nodes power (watt)
Real(dp) :: ferc = 1.e-5	! Flux Error Criteria
Real(dp) :: fserc = 1.e-5	! Fission Source Error Criteria
Integer, Parameter :: Ounit = 100	! Output
REAL(DP)  :: bc1, bc2    ! Boron Concentration
REAL(DP) :: ke1, ke2
INTEGER :: n, j, i
INTEGER :: xtot, ytot, ly, lx

REAL(DP) :: bcon       ! Boron concentration in ppm
Real(dp), Dimension(nk) :: ftem, mtem,cden
Real(dp), Intent(in) :: rftem,cftem,rmtem,cmtem,rcden,ccden
Real(dp) :: st, fn

! File Output
WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) ' ==============================================' &
            // '=================================================='
WRITE(ounit,*) &
               '                       CRITICAL BORON CONCENTRATION SEARCH'
WRITE(ounit,*) ' ==============================================' &
            // '=================================================='
WRITE(ounit,*)
WRITE(ounit,*) '  Itr  Boron Concentration          K-EFF    FLUX REL. ERROR' &
               //'   FISS. SOURCE REL. ERROR    DOPPLER ERROR'
WRITE(ounit,*) ' -----------------------------------------------------------' &
              // '-------------------------------------------'

! Terminal Output
WRITE(*,*)
WRITE(*,*)
WRITE(*,*) ' ==============================================' &
            // '=========='
WRITE(*,*) &
               '           CRITICAL BORON CONCENTRATION SEARCH'
WRITE(*,*) ' ==============================================' &
            // '=========='
WRITE(*,*)
WRITE(*,*) '  Itr  Boron Concentration          K-EFF    '
WRITE(*,*) ' --------------------------------------------------------'

Call Init(keff,jo,ji,L0,al,f0,fx1,fy1,fz1,fx2,fy2,fz2,ng,nk)
!!! Convert assembly wise CR bank map to node wise CR bank map

Call Output_crod (fbmap,ng,nmat,nb,nx,ny,nxx,nyy,nzz,xdiv,ydiv,delz,nstep,bpos,&
				  pos0,ssize,dsiga,dsigtr,dsigf,dnuf,dsigs,bmap)
Call Output_ftem (ftem,ng,nk,nmat,cftem,rftem,fsigs,&
					 fsiga,fsigtr,fsigf,fnuf)
Call Output_mtem (mtem,ng,nk,nmat,cmtem,rmtem,msigs,&
					 msiga,msigtr,msigf,mnuf)
Call Output_cden (cden,ng,nk,nmat,ccden,rcden,lsigs,&
					 lsiga,lsigtr,lsigf,lnuf)

bcon = rbcon

Call XS_updtcfmd(sigs,siga,sigtr,sigf,nu_sigf,D,sigr,ng,nk,nmat,mat,x_sigs,x_siga,x_sigtr,x_sigf,&
				x_nu_sigf,csiga,csigtr,csigf,cnuf,fsiga,fsigtr,fsigf,fnuf,msiga,msigtr,msigf,mnuf,&
				lsiga,lsigtr,lsigf,lnuf,dsiga,dsigtr,dsigf,dnuf,csigs,fsigs,msigs,lsigs,dsigs,bcon,&
				rbcon,rftem,rmtem,rcden,pos0,ssize,f0,fz1,fz2,jo,ji,nxx,nyy,nzz,&
				delz,xyz,nb,ftem,mtem,cden,bpos,fbmap)
CALL Nodal_Coup4(P4,R4,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
Call Nodal_Coup2(P2,R2,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
CALL Outer(Keff,fer,fser,f0,fx1,fy1,fz1,fx2,fy2,fz2,ng,nk,order,R2,P2,R4,P4,nmat,mat,al,jo,ji,L0,&
		   D,sigr,chi,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,xyz,&
		   x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott)
bc1 = bcon
ke1 = Keff

WRITE(ounit,'(I5, F15.2, F23.5, ES16.5, ES21.5, ES22.5)') 1, bc1, Ke1, fser, fer
WRITE(*,'(I5, F15.2, F23.5)') 1, bc1, Ke1

bcon = bcon + (Keff - 1.) * bcon   ! Guess next critical boron concentration

Call XS_updtcfmd(sigs,siga,sigtr,sigf,nu_sigf,D,sigr,ng,nk,nmat,mat,x_sigs,x_siga,x_sigtr,x_sigf,&
				x_nu_sigf,csiga,csigtr,csigf,cnuf,fsiga,fsigtr,fsigf,fnuf,msiga,msigtr,msigf,mnuf,&
				lsiga,lsigtr,lsigf,lnuf,dsiga,dsigtr,dsigf,dnuf,csigs,fsigs,msigs,lsigs,dsigs,bcon,&
				rbcon,rftem,rmtem,rcden,pos0,ssize,f0,fz1,fz2,jo,ji,nxx,nyy,nzz,&
				delz,xyz,nb,ftem,mtem,cden,bpos,fbmap)
CALL Nodal_Coup4(P4,R4,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
Call Nodal_Coup2(P2,R2,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
CALL Outer(Keff,fer,fser,f0,fx1,fy1,fz1,fx2,fy2,fz2,ng,nk,order,R2,P2,R4,P4,nmat,mat,al,jo,ji,L0,&
		   D,sigr,chi,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,xyz,&
		   x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott)
bc2 = bcon
ke2 = Keff

WRITE(ounit,'(I5, F15.2, F23.5, ES16.5, ES21.5, ES22.5)') 2, bc2, Ke2, fser, fer
WRITE(*,'(I5, F15.2, F23.5)') 2, bc2, Ke2

n = 3
DO
  bcon = bc2 + (1._DP - ke2) / (ke1 - ke2) * (bc1 - bc2)

  Call XS_updtcfmd(sigs,siga,sigtr,sigf,nu_sigf,D,sigr,ng,nk,nmat,mat,x_sigs,x_siga,x_sigtr,x_sigf,&
				x_nu_sigf,csiga,csigtr,csigf,cnuf,fsiga,fsigtr,fsigf,fnuf,msiga,msigtr,msigf,mnuf,&
				lsiga,lsigtr,lsigf,lnuf,dsiga,dsigtr,dsigf,dnuf,csigs,fsigs,msigs,lsigs,dsigs,bcon,&
				rbcon,rftem,rmtem,rcden,pos0,ssize,f0,fz1,fz2,jo,ji,nxx,nyy,nzz,&
				delz,xyz,nb,ftem,mtem,cden,bpos,fbmap)
  CALL Nodal_Coup4(P4,R4,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
  Call Nodal_Coup2(P2,R2,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
  CALL Outer(Keff,fer,fser,f0,fx1,fy1,fz1,fx2,fy2,fz2,ng,nk,order,R2,P2,R4,P4,nmat,mat,al,jo,ji,L0,&
			 D,sigr,chi,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,xyz,&
			 x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott)
  bc1 = bc2
  bc2 = bcon
  ke1 = ke2
  ke2 = Keff
  WRITE(ounit,'(I5, F15.2, F23.5, ES16.5, ES21.5, ES22.5)') n, bcon, Keff, fser, fer
  WRITE(*,'(I5, F15.2, F23.5)') n, bcon, Keff
    IF ((ABS(Keff - 1._DP) < 1.e-5_DP) .AND. (fser < 1.e-5_DP) .AND. (fer < 1.e-5_DP)) EXIT
    n = n + 1
    IF (bcon > 3000.) THEN
        WRITE(ounit,*) '  CRITICAL BORON CONCENTRATION EXCEEDS THE LIMIT(3000 ppm)'
        WRITE(ounit,*) '  OPENNODE IS STOPPING'
        WRITE(*,*) '  CRITICAL BORON CONCENTRATION EXCEEDS THE LIMIT(3000 ppm)'
        STOP
    END IF
    IF (bcon < 0.) THEN
        WRITE(ounit,*) '  CRITICAL BORON CONCENTRATION IS NOT FOUND (LESS THAN ZERO)'
        WRITE(ounit,*) '  OPENNODE IS STOPPING'
        WRITE(*,*) '  CRITICAL BORON CONCENTRATION IS NOT FOUND (LESS THAN ZERO)'
        STOP
    END IF
    IF (n == 20) THEN
        WRITE(ounit,*) '  MAXIMUM ITERATION FOR CRITICAL BORON SEARCH IS REACHING MAXIMUM'
        WRITE(ounit,*) '  OPENNODE IS STOPPING'
        WRITE(*,*) '  MAXIMUM ITERATION FOR CRITICAL BORON SEARCH IS REACHING MAXIMUM'
        STOP
    END IF
END DO

	CALL PowDis(power,f0,sigf,delv,ng,nk,mode)

    Call AsmPow(nk,nx,ny,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,xdiv,ydiv,y_smax,y_smin,power)

    Call AxiPow(nk,nz,nxx,nyy,nzz,ix,iy,iz,xyz,delz,delv,zdiv,y_smax,y_smin,power)

    Call AsmFlux(ng,nk,nx,ny,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,xdiv,ydiv,y_smax,y_smin,f0,1._DP)

Call CPU_TIME(fn)
Write(ounit,*)
Write(*,*)
Write(ounit,*) "Total Time : ", fn-st, "Seconds"
Write(*,*) "Total Time : ", fn-st, "Seconds"

END SUBROUTINE

	SUBROUTINE cbsearchth(ng,nk,nb,nx,ny,nz,nmat,order,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,chi,&
						  nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,x_smax,x_smin,xdiv,ydiv,zdiv,xsize,ysize,&
					      y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott,mode,rbcon,&
 					      csiga,csigtr,csigf,cnuf,csigs,cftem,rftem,fsiga,fsigtr,fsigf,fnuf,fsigs,&
					      msiga,msigtr,msigf,mnuf,msigs,rmtem,cmtem,&
					      lsiga,lsigtr,lsigf,lnuf,lsigs,rcden,ccden,dsiga,dsigtr,dsigf,dnuf,dsigs,&
				   	      pos0,ssize,bpos,bmap,nstep,&
					      ntem,stab,cf,tin,nfpin,cmflow,rf,tg,tc,ppitch,pra,powr,ppow)!    To search critical boron concentration with thermal feedback
						
Implicit None
Integer, Parameter :: dp = Selected_Real_Kind(10,15)
! #######################   User Input Paramaters   #######################
Integer, Intent(in) :: ng, nk, pra,nb, order, nx, ny, nz, nmat, nxx, nyy, nzz, x_east, x_west, y_north, y_south, z_top, z_bott
Integer, Dimension(nxx), Intent(in) :: x_smax, x_smin
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin
Integer, Dimension(nk), Intent(in) :: mat
Integer, Dimension(nx), Intent(in) :: xdiv
Integer, Dimension(ny), Intent(in) :: ydiv
Integer, Dimension(nz), Intent(in) :: zdiv
Real(dp), Dimension(nx), Intent(in) :: xsize           ! Assembly Size & Division
Real(dp), Dimension(ny), Intent(in) :: ysize           ! Assembly Size & Division
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk), Intent(in) :: delv, ix, iy, iz
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz
Character(Len=100), intent(in) :: mode
INTEGER, DIMENSION(nx, ny), intent(in) :: bmap       ! Radial control rod bank map (assembly wise)

! Thermal-hydraulics parameters
INTEGER,  INTENT(IN) :: ntem    ! Number of temperature in steam table
REAL(DP), DIMENSION(ntem,pra), intent(in) :: stab  ! Steam table matrix
REAL(DP), intent(in) :: cf		! heat fraction deposited into coolant
REAL(DP), INTENT(IN) :: tin		! coolant inlet temperature (kelvin)
Real(dp), intent(in) :: cftem, cmtem, ccden
Real(dp), Dimension(nb), intent(in) :: bpos
INTEGER, intent(in) :: nfpin		! Number of fuel pin and guide tubes
REAL(DP), intent(in) :: cmflow
REAL(DP), intent(in) :: rf, tg, tc, ppitch		! Fuel meat radius, gap thickness, clad thickness, and pin picth (m)
INTEGER, PARAMETER :: nm = 10      ! Fuel meat divided into 10 mesh
INTEGER, PARAMETER :: nt = nm + 2      ! two more mesh for gap and clad
REAL(DP), DIMENSION(nt) :: rpos	! mesh position
REAL(DP), DIMENSION(nk,nt+1):: tfm	! Fuel pin mesh temperature for each nodes
! CXs Assigned to Materials
Real(dp), Dimension(nmat,ng), Intent(in) :: chi, x_siga, x_sigtr, x_sigf, x_nu_sigf
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: x_sigs
Real(dp), Dimension(nmat,ng), Intent(in) :: csiga, csigtr, csigf, cnuf,&
fsiga, fsigtr, fsigf, fnuf, msiga, msigtr, msigf, mnuf, &
lsiga, lsigtr, lsigf, lnuf, dsiga, dsigtr, dsigf, dnuf
Real(dp), Intent(in) :: rbcon ,rftem, rmtem, rcden
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: csigs, fsigs, msigs, lsigs, dsigs
REAL(DP), Intent(in):: powr, ppow
REAL(DP), DIMENSION(nxx, nyy) :: node_nf	! Number of fuel pin per node
REAL(DP), INTENT(IN) :: pos0, ssize			! Zero step position and step size
Real(dp), Dimension(nk) :: ftem, mtem, cden
INTEGER, DIMENSION(nxx,nyy):: fbmap	! Radial control rod bank map (node wise)
REAL(DP), Intent(in) :: nstep		! Number of steps
Real(dp), Dimension(nk) :: fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2     ! Fission Source & Fission Source Moments
Real(dp) :: st, fn

!##########################################################################
! Local Variables
Real(dp), Dimension(nk,ng) :: f0     ! Flux
! CXs Assigned to Nodes
Real(dp), Dimension(nk,ng) :: sigf	! Fission Macroscopic CX
Real(dp), Dimension(nk,ng) :: nu_sigf	! nu * Fission Macroscopic CX
REAL(DP), DIMENSION(nk):: power	! nodes power (watt)
REAL(DP) :: th_err		! Doppler error
Real(dp) :: ferc = 1.e-5	! Flux Error Criteria
Real(dp) :: fserc = 1.e-5	! Fission Source Error Criteria
Integer, Parameter :: Ounit = 100	! Output
REAL(DP)  :: bc1, bc2    ! Boron Concentration
REAL(DP) :: ke1, ke2
INTEGER :: n, no = 10
REAL(DP) :: tf, tm, mtm, mtf
Real(dp) :: Keff	! Effectif Multiplication Factor
Real(dp) :: fer, fser	! Flux and Fission Source Error in BCSEARCH calcs.
REAL(DP) :: bcon       ! Boron concentration in ppm
Real(dp), Dimension(nk,ng,6) :: jo   ! Nodals' Outgoing Currents (X+,X-,Y+,Y-,Z+,Z-)
Real(dp), Dimension(nk,ng,6) :: ji   ! Nodals' Ingoing Currents  (X+,X-,Y+,Y-,Z+,Z-)
Real(dp), Dimension(nk,ng,3) :: L0   ! Zeroth Transverse Leakages (Lx,Ly,Lz)
Real(dp), Dimension(nk,ng,6) :: al   ! Assembly Discontinuity Factor
Real(dp), Dimension(nk,ng) :: fx1, fy1, fz1, fx2, fy2, fz2      ! Flux and Flux Moments
REAL(DP), DIMENSION(nk) :: heatf	! Heat flux (W/m2)
REAL(DP) :: dum
REAL(DP), DIMENSION(nt) :: rdel	! mesh delta
REAL(DP) :: rg, rc	! Outer radius of gap and cladding
INTEGER :: i
REAL(DP) :: dia	! Pi diameter, Hydraulic diameter (m) and sub-channel area (m2)
REAL(DP), PARAMETER :: pi = 3.14159265
REAL(DP) :: cflow
REAL(DP) :: dh, farea	! Pi diameter, Hydraulic diameter (m) and sub-channel area (m2)
REAL(DP), DIMENSION(nk) :: ent	! Coolant Enthalpy (J/Kg)

! Guess fuel and moderator temperature
! Used for radial fuel temperature distribution
tfm = 1200.
! Initial heat-flux rate
heatf = 0.

dum = rf / REAL(nm)
DO i = 1, nm
    rdel(i) = dum
END DO

! Fuel pin mesh size
rdel(nm+1) = tg
rdel(nm+2) = tc
! Fuel pin mesh position
rpos(1) = 0.5 * rdel(1)
DO i = 2, nt
    rpos(i) = rpos(i-1) + 0.5 * (rdel(i-1) + rdel(i))
END DO
!write(*,*) 'rpos=', rpos(:)
! Calculate outer radius of gap and cladding
rg = rf + tg
rc = rf + tg + tc
! Calculate fuel pin diameter
dia = 2. * rc
! Calculate hydraulic diameter
dh = dia * ((4./pi) * (ppitch/dia)**2 - 1.)
! Calculate sub-channel mass flow rate
cflow = cmflow / REAL(nfpin)
! Calculate sub-channel area
farea = ppitch**2 - 0.25*pi*dia**2

!WRITE(*,*) 'stab= ', stab
!WRITE(*,*) 'tin= ', tin
!WRITE(*,*) 'nfpin= ', nfpin
!WRITE(*,*) 'cmflow= ', cmflow
!WRITE(*,*) 'rf= ', rf
!WRITE(*,*) 'tg= ', tg
!WRITE(*,*) 'tc= ', tc
!WRITE(*,*) 'ppitch= ', ppitch
!WRITE(*,*) 'pra= ', pra
!WRITE(*,*) 'powr= ', powr
!WRITE(*,*) 'ppow= ', ppow


Call Output_cbcs (ng,nmat,rbcon,csigtr,csiga,cnuf,csigf,csigs)
Call Output_crod (fbmap,ng,nmat,nb,nx,ny,nxx,nyy,nzz,xdiv,ydiv,delz,nstep,bpos,&
				  pos0,ssize,dsiga,dsigtr,dsigf,dnuf,dsigs,bmap)
Call Output_ther (node_nf,ng,nk,nx,ny,nmat,nxx,nyy,&
				  y_smax,y_smin,xdiv,xsize,ydiv,ysize,powr,ppow,tin,&
				  rf,tg,tc,ppitch,cf,pra,&
				  cmflow,nfpin,ntem,stab)
Call Output_ftem (ftem,ng,nk,nmat,cftem,rftem,fsigs,&
					 fsiga,fsigtr,fsigf,fnuf)
Call Output_mtem (mtem,ng,nk,nmat,cmtem,rmtem,msigs,&
					 msiga,msigtr,msigf,mnuf)
Call Output_cden (cden,ng,nk,nmat,ccden,rcden,lsigs,&
					 lsiga,lsigtr,lsigf,lnuf)

! File Output
WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) ' ==============================================' &
            // '=================================================='
WRITE(ounit,*) &
               '                       CRITICAL BORON CONCENTRATION SEARCH'
WRITE(ounit,*) ' ==============================================' &
            // '=================================================='
WRITE(ounit,*)
WRITE(ounit,*) '  Itr  Boron Concentration          K-EFF    FLUX REL. ERROR' &
               //'   FISS. SOURCE REL. ERROR    DOPPLER ERROR'
WRITE(ounit,*) ' -----------------------------------------------------------' &
              // '-------------------------------------------'

! Terminal Output
WRITE(*,*)
WRITE(*,*)
WRITE(*,*) ' ==============================================' &
            // '=========='
WRITE(*,*) &
               '           CRITICAL BORON CONCENTRATION SEARCH'
WRITE(*,*) ' ==============================================' &
            // '=========='
WRITE(*,*)
WRITE(*,*) '  Itr  Boron Concentration          K-EFF    '
WRITE(*,*) ' --------------------------------------------------------'

Call Init(keff,jo,ji,L0,al,f0,fx1,fy1,fz1,fx2,fy2,fz2,ng,nk)
bcon = rbcon

CALL th_iter(th_err,tfm,Keff,fer,fser,f0,sigf,nu_sigf,ent,ng,nk,nx,ny,nz,nmat,order,mat,x_sigs,x_siga,&
			 x_sigtr,x_sigf,x_nu_sigf,chi,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,&
			 x_smax,x_smin,xdiv,ydiv,zdiv,y_smax,y_smin,xyz,x_east,x_west,y_north,&
			 y_south,z_top,z_bott,mode,ntem,stab,cf,tin,ftem,mtem,cden,bcon,bpos,&
			 rf,tg,tc,ppitch,pra,jo,ji,L0,al,fx1,fx2,fy1,fy2,fz1,fz2,heatf,dh,farea,&
			 csiga,csigtr,csigf,cnuf,fsiga,fsigtr,fsigf,fnuf,msiga,msigtr,msigf,&
			 mnuf,lsiga,lsigtr,lsigf,lnuf,dsiga,dsigtr,dsigf,dnuf,rbcon,rftem,rmtem,&
			 rcden,csigs,fsigs,msigs,lsigs,dsigs,powr,ppow,node_nf,pos0,ssize,nb,fbmap,&
			 cflow,dia,rg,rc,rpos,rdel,fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2)    ! Start thermal hydarulic iteration with current paramters
bc1 = bcon
ke1 = Keff

WRITE(ounit,'(I5, F15.2, F23.5, ES16.5, ES21.5, ES22.5)') 1, bc1, Ke1, fser, fer, th_err
WRITE(*,'(I5, F15.2, F23.5)') 1, bc1, Ke1

bcon = bcon + (Keff - 1.) * bcon   ! Guess next critical boron concentration
CALL th_iter(th_err,tfm,Keff,fer,fser,f0,sigf,nu_sigf,ent,ng,nk,nx,ny,nz,nmat,order,mat,x_sigs,x_siga,&
			 x_sigtr,x_sigf,x_nu_sigf,chi,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,&
			 x_smax,x_smin,xdiv,ydiv,zdiv,y_smax,y_smin,xyz,x_east,x_west,y_north,&
			 y_south,z_top,z_bott,mode,ntem,stab,cf,tin,ftem,mtem,cden,bcon,bpos,&
			 rf,tg,tc,ppitch,pra,jo,ji,L0,al,fx1,fx2,fy1,fy2,fz1,fz2,heatf,dh,farea,&
			 csiga,csigtr,csigf,cnuf,fsiga,fsigtr,fsigf,fnuf,msiga,msigtr,msigf,&
			 mnuf,lsiga,lsigtr,lsigf,lnuf,dsiga,dsigtr,dsigf,dnuf,rbcon,rftem,rmtem,&
			 rcden,csigs,fsigs,msigs,lsigs,dsigs,powr,ppow,node_nf,pos0,ssize,nb,fbmap,&
			 cflow,dia,rg,rc,rpos,rdel,fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2)                 ! Perform second thermal hydarulic iteration with updated parameters
bc2 = bcon
ke2 = Keff

WRITE(ounit,'(I5, F15.2, F23.5, ES16.5, ES21.5, ES22.5)') 2, bc2, Ke2, fser, fer, th_err
WRITE(*,'(I5, F15.2, F23.5)') 2, bc2, Ke2

n = 3
DO 
    bcon = bc2 + (1._DP - ke2) / (ke1 - ke2) * (bc1 - bc2)
    CALL th_iter(th_err,tfm,Keff,fer,fser,f0,sigf,nu_sigf,ent,ng,nk,nx,ny,nz,nmat,order,mat,x_sigs,x_siga,&
				x_sigtr,x_sigf,x_nu_sigf,chi,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,&
				x_smax,x_smin,xdiv,ydiv,zdiv,y_smax,y_smin,xyz,x_east,x_west,y_north,&
				y_south,z_top,z_bott,mode,ntem,stab,cf,tin,ftem,mtem,cden,bcon,bpos,&
				rf,tg,tc,ppitch,pra,jo,ji,L0,al,fx1,fx2,fy1,fy2,fz1,fz2,heatf,dh,farea,&
				csiga,csigtr,csigf,cnuf,fsiga,fsigtr,fsigf,fnuf,msiga,msigtr,msigf,&
				mnuf,lsiga,lsigtr,lsigf,lnuf,dsiga,dsigtr,dsigf,dnuf,rbcon,rftem,rmtem,&
				rcden,csigs,fsigs,msigs,lsigs,dsigs,powr,ppow,node_nf,pos0,ssize,nb,fbmap,&
				cflow,dia,rg,rc,rpos,rdel,fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2)  
    bc1 = bc2
    bc2 = bcon
    ke1 = ke2
    ke2 = Keff
    WRITE(ounit,'(I5, F15.2, F23.5, ES16.5, ES21.5, ES22.5)') n, bcon, Keff, fser, fer, th_err
    WRITE(*,'(I5, F15.2, F23.5)') n, bcon, Keff
    IF ((ABS(Keff - 1._DP) < 1.e-5_DP) .AND. (fser < fserc) .AND. (fer < ferc)) EXIT
    n = n + 1
    IF (bcon > 3000.) THEN
        WRITE(ounit,*) '  CRITICAL BORON CONCENTRATION EXCEEDS THE LIMIT(3000 ppm)'
        WRITE(ounit,*) '  OPENNODE IS STOPPING'
        WRITE(*,*) '  CRITICAL BORON CONCENTRATION EXCEEDS THE LIMIT(3000 ppm)'
        STOP
    END IF
    IF (bcon < 0.) THEN
        WRITE(ounit,*) '  CRITICAL BORON CONCENTRATION IS NOT FOUND (LESS THAN ZERO)'
        WRITE(ounit,*) '  OPENNODE IS STOPPING'
        WRITE(*,*) '  CRITICAL BORON CONCENTRATION IS NOT FOUND (LESS THAN ZERO)'
        STOP
    END IF
    IF (n == 30) THEN
        WRITE(ounit,*) '  MAXIMUM ITERATION FOR CRITICAL BORON SEARCH IS REACHING MAXIMUM'
        WRITE(ounit,*) '  OPENNODE IS STOPPING'
        WRITE(*,*) '  MAXIMUM ITERATION FOR CRITICAL BORON SEARCH IS REACHING MAXIMUM'
        exit
    END IF
END DO

	CALL PowDis(power,f0,sigf,delv,ng,nk,mode)

    Call AsmPow(nk,nx,ny,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,xdiv,ydiv,y_smax,y_smin,power)

    Call AxiPow(nk,nz,nxx,nyy,nzz,ix,iy,iz,xyz,delz,delv,zdiv,y_smax,y_smin,power)

    Call AsmFlux(ng,nk,nx,ny,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,xdiv,ydiv,y_smax,y_smin,f0,1._DP)

CALL par_ave_f(tf,ftem,nk,ng,nu_sigf,delv)
CALL par_ave(tm,mtem,nk,delv)
Call par_max(mtf,tfm(:,1),nk)
CALL par_max(mtm,mtem,nk)
CALL getfq(nk,nxx,nyy,nzz,ix,iy,iz,delz,node_nf,power)

! Write Output
WRITE(ounit,*)
WRITE(ounit, 5001) tf, tf-273.15
WRITE(ounit, 5002)  mtf, mtf-273.15
WRITE(ounit, 5003) tm, tm-273.15
WRITE(ounit, 5004) mtm, mtm-273.15

5001 FORMAT(2X, 'AVERAGE DOPPLER TEMPERATURE     : ', F7.1, ' K (', F7.1, ' C)')
5002 FORMAT(2X, 'MAX FUEL CENTERLINE TEMPERATURE : ', F7.1, ' K (', F7.1, ' C)')
5003 FORMAT(2X, 'AVERAGE MODERATOR TEMPERATURE   : ', F7.1, ' K (', F7.1, ' C)')
5004 FORMAT(2X, 'MAXIMUM MODERATOR TEMPERATURE   : ', F7.1, ' K (', F7.1, ' C)')

Call CPU_TIME(fn)
Write(ounit,*)
Write(*,*)
Write(ounit,*) "Total Time : ", fn-st, "Seconds"
Write(*,*) "Total Time : ", fn-st, "Seconds"

END SUBROUTINE

	SUBROUTINE rod_eject(ng,nk,nb,nx,ny,nz,nmat,order,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,chi,&
			     nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,x_smax,x_smin,xdiv,ydiv,zdiv,xsize,ysize,&
			     y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott,mode,&
			     dsiga,dsigtr,dsigf,dnuf,dsigs,pos0,ssize,bpos,bmap,nstep,fbpos,tmove,bspeed,ibeta,&
			     lamb,velo,ttot,tstep1,tstep2,tdiv)
!    To perform rod ejection simulation

Implicit None
Integer, Parameter :: dp = Selected_Real_Kind(10,15)
! #######################   User Input Paramaters   #######################
Integer, Intent(in) :: ng, nk, nb, order, nx, ny, nz, nmat, nxx, nyy, nzz, x_east, x_west, y_north, y_south, z_top, z_bott
Integer, Dimension(nxx), Intent(in) :: x_smax, x_smin
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin
Integer, Dimension(nk), Intent(in) :: mat
Integer, Dimension(nx), Intent(in) :: xdiv
Integer, Dimension(ny), Intent(in) :: ydiv
Integer, Dimension(nz), Intent(in) :: zdiv
Real(dp), Dimension(nx), Intent(in) :: xsize           ! Assembly Size & Division
Real(dp), Dimension(ny), Intent(in) :: ysize           ! Assembly Size & Division
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk), Intent(in) :: delv, ix, iy, iz
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz
Character(Len=100), intent(in) :: mode
Real(dp), Dimension(nmat,ng), Intent(inout) :: chi, x_siga, x_sigtr, x_sigf, x_nu_sigf
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: x_sigs
Real(dp), Dimension(nmat,ng), Intent(inout) :: dsiga, dsigtr, dsigf, dnuf
REAL(DP), INTENT(IN) :: pos0, ssize			! Zero step position and step size
INTEGER, DIMENSION(nx, ny), intent(in) :: bmap       ! Radial control rod bank map (assembly wise)
REAL(DP), Intent(in) :: nstep		! Number of steps
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: dsigs

REAL(dp), Dimension(nk) :: ftem, mtem, cden
INTEGER, DIMENSION(nxx,nyy) :: fbmap	! Radial control rod bank map (node wise)
! Crod changes
REAL(DP), DIMENSION(nb), intent(in) :: bpos  ! CR bank position
REAL(DP), DIMENSION(nb), intent(in) :: fbpos    ! Final CR bank position
REAL(DP), DIMENSION(nb), intent(in) :: tmove    ! Time when CR bank starts moving
REAL(DP), DIMENSION(nb), intent(in) :: bspeed   ! CR bank movement speed
INTEGER, DIMENSION(nb) :: mdir  ! To indicate CR movement direction (0=do not move, 1=down, 2 = up)
! Transient parameters
INTEGER, PARAMETER :: nf = 6			! Number of delaye dneutron precusor family
REAL(DP), DIMENSION(nf), intent(in) :: ibeta, lamb                 ! beta (delayed neutron fraction) and precusor decay constant
REAL(DP) :: tbeta                                      ! total beta
REAL(DP), DIMENSION(ng), intent(in) :: velo            ! Neutron velocity
REAL(DP), intent(in) :: ttot                                       ! TOTAL SIMULATION TIME
REAL(DP), intent(in) :: tstep1                                     ! FIRST TIME STEP
REAL(DP), intent(in) :: tstep2                                     ! SECOND TIME STEP
REAL(DP), intent(in) :: tdiv                                       ! WHEN SECOND TIME STEP APPLY
REAL(DP), DIMENSION(nk,ng) :: omeg         ! Exponential transformation constant
LOGICAL :: tranw = .FALSE.                        ! To activate unconverged  outer iteration warning

Real(dp), Dimension(nk) :: fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2     ! Fission Source & Fission Source Moments
REAL(DP), DIMENSION(nk,nf) :: c0, cx1, cy1, cz1, cx2, cy2, cz2  ! neutron precusor density
REAL(DP), DIMENSION(nk,ng) :: ft, ftx1, fty1, ftz1, ftx2, fty2, ftz2  ! Parameters at previous time step
! Local Variables
Real(dp) :: Keff   ! Effectif Multiplication Factor
Real(dp) :: fer, fser	! Flux and Fission Source Error in BCSEARCH calcs.
Real(dp), Dimension(nk,ng) :: f0, fx1, fy1, fz1, fx2, fy2, fz2      ! Flux and Flux Moments
Real(dp), Dimension(nk,ng,6) :: jo   ! Nodals' Outgoing Currents (X+,X-,Y+,Y-,Z+,Z-)
Real(dp), Dimension(nk,ng,6) :: ji   ! Nodals' Ingoing Currents  (X+,X-,Y+,Y-,Z+,Z-)
Real(dp), Dimension(nk,ng,3) :: L0   ! Zeroth Transverse Leakages (Lx,Ly,Lz)
Real(dp), Dimension(nk,ng,6) :: al   ! Assembly Discontinuity Factor
Real(dp), Dimension(nk,ng,6,6) :: R4, R2, P2    ! Response Matrix
Real(dp), Dimension(nk,ng,6,7) :: P4            ! Response Matrix
 Real(dp), Dimension(nk,ng,7) :: Q      ! Nodal's Source and Source Moments (0,x1,y1,z1,x2,y2,z2)
! CXs Assigned to Nodes
Real(dp), Dimension(nk,ng,ng) :: sigs  ! Scattering Macroscopic CX
Real(dp), Dimension(nk,ng) :: siga     ! Absorption Macroscopic CX
Real(dp), Dimension(nk,ng) :: sigtr    ! Transport Macroscopic CX
Real(dp), Dimension(nk,ng) :: sigf     ! Fission Macroscopic CX
Real(dp), Dimension(nk,ng) :: nu_sigf  ! nu * Fission Macroscopic CX
Real(dp), Dimension(nk,ng) :: D        ! Diffusion Coefficient
Real(dp), Dimension(nk,ng) :: sigr     ! Scattering Macroscopic CX
REAL(DP), DIMENSION(nk):: power	! nodes power (watt)
Real(dp) :: ferc = 1.e-5	! Flux Error Criteria
Real(dp) :: fserc = 1.e-5	! Fission Source Error Criteria
Integer, Parameter :: Ounit = 100	! Output
REAL(DP), DIMENSION(nk, ng) :: af		! adjoint flux
REAL(DP), DIMENSION(nk, ng) :: sigrp	! Temporary sigr
Real(dp), Dimension(nk) :: errn, erro
real(dp), dimension(nb) :: bpos_local        ! Local copy of bpos
Real(dp) :: st, fn

REAL(DP) :: rho
REAL(DP) :: t1, t2
REAL(DP) :: tpow1, tpow2
INTEGER :: n, i, j, g, imax, step
LOGICAL :: maxi   ! Maximum Outer Iteration Reached?

Call Init(keff,jo,ji,L0,al,f0,fx1,fy1,fz1,fx2,fy2,fz2,ng,nk)
Call Output_crod (fbmap,ng,nmat,nb,nx,ny,nxx,nyy,nzz,xdiv,ydiv,delz,nstep,bpos,&
				  pos0,ssize,dsiga,dsigtr,dsigf,dnuf,dsigs,bmap)
Call Output_ejct (mdir,ng,nk,nb,nmat,bpos,fbpos,tmove,bspeed,&
					ibeta,lamb,velo,ttot,tstep1,tstep2,tdiv)
					
! Update xsec
CALL XS_updtrod(sigs,siga,sigtr,sigf,nu_sigf,D,sigr,ng,nk,nmat,mat,x_sigs,x_siga,x_sigtr,x_sigf,&
				x_nu_sigf,dsiga,dsigtr,dsigf,dnuf,dsigs,pos0,ssize,f0,fz1,fz2,jo,ji,nxx,nyy,nzz,&
				delz,xyz,nb,bpos,fbmap)
! Calculate forward flux at t=0 and check if keff=1
CALL Nodal_Coup4(P4,R4,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
Call Nodal_Coup2(P2,R2,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
CALL Outerfs(Keff,fer,fser,f0,fx1,fy1,fz1,fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,ng,nk,&
			order,R2,P2,R4,P4,nmat,mat,al,jo,ji,L0,D,sigr,chi,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,&
			delx,dely,delz,delv,xyz,x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott)

!Write(*,*) 'before CALL KNE1=' , 'keff=', Keff
!DO n = 1, nk
!	if (n>50 .and. n<65) then
!		write(*,*) 'f0=', f0(n,:)
!		write(*,*) 'fs0=', fs0(n)
!	endif
!enddo
! If K-EFF NOT EQUAL TO 1.0
IF (ABS(Keff - 1._DP) > 1.e-5_DP) Then 
CALL KNE1(ng,nk,nmat,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,dsiga,dsigtr,dsigf,dnuf,dsigs,&
			pos0,ssize,nxx,nyy,nzz,delz,xyz,nb,bpos,fbmap,Keff,f0,fx1,fy1,fz1,fx2,fy2,fz2,&
			fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,order,R2,P2,R4,P4,al,jo,ji,L0,chi,ix,iy,iz,delx,dely,&
			delv,x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott)
Endif
! Update xsec
CALL XS_updtrod(sigs,siga,sigtr,sigf,nu_sigf,D,sigr,ng,nk,nmat,mat,x_sigs,x_siga,x_sigtr,x_sigf,&
				x_nu_sigf,dsiga,dsigtr,dsigf,dnuf,dsigs,pos0,ssize,f0,fz1,fz2,jo,ji,nxx,nyy,nzz,&
				delz,xyz,nb,bpos,fbmap)

!Write(*,*) 'before CALL OuterAdjfs=' , 'keff=', Keff
!DO n = 1, nk
!	if (n>50 .and. n<65) then
!		write(*,*) 'n=', n, 'f0=', f0(n,:)
!		write(*,*) 'n=', n, 'fs0=', fs0(n)
!		write(*,*) 'n=', n, 'nu_sigf=', nu_sigf(n,:)
!	endif
!enddo
! Calculate Adjoint flux
!Write(*,*) 'CALL OuterAdjfs=' 
CALL Nodal_coup4(P4,R4,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
Call Nodal_Coup2(P2,R2,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
CALL OuterAdjfs(Keff,f0,fx1,fy1,fz1,fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,ng,nk,order,&
				R2,P2,R4,P4,nmat,chi,mat,al,jo,ji,L0,D,sigr,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,&
				delx,dely,delz,delv,x_smax,x_smin,y_smax,y_smin,xyz,x_east,x_west,y_north,&
				y_south,z_top,z_bott)
!CALL Outerfs(Keff,fer,fser,f0,fx1,fy1,fz1,fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,ng,nk,&
!			order,R2,P2,R4,P4,nmat,mat,al,jo,ji,L0,D,sigr,chi,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,&
!			delx,dely,delz,delv,xyz,x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott)
!
!Write(*,*) 'after CALL OuterAdjfs=' , 'keff=', Keff
!DO n = 1, nk
!	if (n>50 .and. n<65) then
!		write(*,*) 'n=', n, 'f0=', f0(n,:)
!		write(*,*) 'n=', n, 'fs0=', fs0(n)
!		write(*,*) 'n=', n, 'nu_sigf=', nu_sigf(n,:)
!	endif
!enddo

af = f0   ! Save adjoint flux to af
  DO n = 1, nk
		if (n>40 .and. n>43) then
		!	write(*,*) 'f0=', f0(n,:)
		!	write(*,*) 'fs0=', fs0(n)
	!		write(*,*) 'af=', af(n,:)
		endif
	enddo
! ReCalculate forward flux
!Write(*,*) 'CALL Outer=' 
CALL Outerfs(Keff,fer,fser,f0,fx1,fy1,fz1,fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,ng,nk,&
			order,R2,P2,R4,P4,nmat,mat,al,jo,ji,L0,D,sigr,chi,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,&
			delx,dely,delz,delv,xyz,x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott)
!Write(*,*) 'CALL Outerfs=' 
  DO n = 1, nk
		if (n>40 .and. n>43) then
	!		write(*,*) 'f0=', f0(n,:)
	!		write(*,*) 'fs0=', fs0(n)
	!		write(*,*) 'af=', af(n,:)
		endif
	enddo
! Calculate Initial precursor density
CALL iPden(c0,cx1,cy1,cz1,cx2,cy2,cz2,nk,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,iBeta,lamb)

! Calculate initial power
CALL PowTot (tpow1,nk,ng,f0,sigf,delv)
! Total beta
tbeta = 0.
DO j = 1, nf
  tbeta = tbeta + iBeta(j)
END DO

! Calculate reactivity
CALL react(rho,nk,ng,nmat,nxx,nyy,nzz,mat,sigs,chi,f0,fs0,delx,dely,delz,delv,ix,iy,iz,&
		   L0,af,sigr)
		   
!WRITE(*, *) " RHO1=", rho

! File output
WRITE(ounit, *)
WRITE(ounit, *) " TRANSIENT RESULTS :"
WRITE(ounit, *)
WRITE(ounit, *) " Step  Time(s)  React.($)   Rel. Power   CR Bank Pos. (1-end)"
WRITE(ounit, *) "--------------------------------------------------------------"
WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 12F9.2)') 0, 0., 0., &
1.0, (bpos(n), n = 1, nb)

! Terminal output
WRITE(*,*)
WRITE(*,*)
WRITE(*, *) " TRANSIENT RESULTS :"
WRITE(*, *)
WRITE(*, *) " Step  Time(s)  React.($)   Rel. Power"
WRITE(*, *) "-------------------------------------------"
WRITE(*,'(I4, F10.3, F10.4, ES15.4)') 0, 0., 0., 1.0

! Start transient calculation
step = 0
t2 = 0.
imax = NINT(tdiv/tstep1)

!Do g= 1, ng
!    Do n = 1, nk
!		ft(n,g)  = 1._dp
!    End Do
!End Do

! Initialize the local copy of bpos
bpos_local = bpos

! First Time Step
DO i = 1, imax
!write(*,*) '***********************************'
!write(*,*) 'Step=', i
!write(*,*) '***********************************'
    step = step + 1
    t1 = t2
    t2 = REAL(i)*tstep1

    IF (i > 1) THEN
       omeg = LOG(f0 / ft) / tstep1
    ELSE
       omeg = 0.
    END IF

    ! Rod bank changes
    DO n = 1, nb
        IF (mdir(n) == 1) THEN   ! If CR moving down
            IF (t2-tmove(n) > 1.e-5_DP .AND. fbpos(n)-bpos_local(n) < 1.e-5_DP) THEN
                bpos_local(n) = bpos_local(n) - tstep1 *  bspeed(n)
                IF (bpos_local(n) < fbpos(n)) bpos_local(n) = fbpos(n)  ! If bpos exceed, set to fbpos
            END IF
        ELSE IF (mdir(n) == 2) THEN ! If CR moving up
            IF (t2-tmove(n) > 1.e-5_DP .AND. fbpos(n)-bpos_local(n) > 1.e-5_DP) THEN
                bpos_local(n) = bpos_local(n) + tstep1 *  bspeed(n)
                IF (bpos_local(n) > fbpos(n)) bpos_local(n) = fbpos(n)  ! If bpos exceed, set to fbpos
            END IF
        ELSE
            CONTINUE
        END IF
     END DO

			!	do n=1, nb
			!			write(*,*) '***********************************'
			!			write(*,*) '***********************************'
			!				write(*,*) 'AFTER Rod bank changes'
			!				write(*,*) 'Rod bank=', n, 'bpos=', bpos(n)
			!			write(*,*) '***********************************'
			!	enddo

    ! Calculate xsec after pertubation
    CALL XS_updtrod(sigs,siga,sigtr,sigf,nu_sigf,D,sigr,ng,nk,nmat,mat,x_sigs,x_siga,x_sigtr,x_sigf,&
					x_nu_sigf,dsiga,dsigtr,dsigf,dnuf,dsigs,pos0,ssize,f0,fz1,fz2,jo,ji,nxx,nyy,nzz,&
					delz,xyz,nb,bpos_local,fbmap)

!	write(*,*) 'AFTER call XS_updtcfmd AFTER Rod bank changes'
!DO n = 1, nk
!	if (n>40 .and. n<43) then
!		write(*,*) 'time step', i, 'sigs=', sigs(n,:,:)
!		write(*,*) 'time step', i, 'siga=', siga(n,:)
!		write(*,*) 'time step', i, 'sigtr=', sigtr(n,:)
!		write(*,*) 'time step', i, 'D=', D(n,:)
!		write(*,*) 'time step', i, 'sigr=', sigr(n,:)
!		write(*,*) 'time step', i, 'sigf=', sigf(n,:)
!		write(*,*) 'time step', i, 'nu_sigf=', nu_sigf(n,:)
!	endif
!enddo

    ! Modify removal xsec
    sigrp = sigr    ! Save sigr to sigrp
    DO g = 1, ng
       DO n = 1, nk
          sigr(n,g) = sigr(n,g) + 1. / (velo(g) * tstep1) + omeg(n,g) / velo(g)
       END DO
    END DO

!write(*,*) 'AFTER Rod bank changes'
do g = 1, ng
		DO n = 1, nk
			if (n>1 .and. n<70) then
			!	write(*,*) 'group=', g, 'node', n, 'sigs=', sigs(n,g,g)
			!	write(*,*) 'group=', g, 'node', n, 'siga=', siga(n,g)
			!	write(*,*) 'group=', g, 'node', n, 'sigtr=', sigtr(n,g)
			!	write(*,*) 'group=', g, 'node', n, 'D=', D(n,g)
			!	write(*,*) 'group=', g, 'node', n, 'sigr=', sigr(n,g)
			!	write(*,*) 'group=', g, 'node', n, 'sigf=', sigf(n,g)
			!	write(*,*) 'group=', g, 'node', n, 'nu_sigf=', nu_sigf(n,g)
			endif
		enddo
enddo
		
    ! Save the previous fluxes
    ft = f0
    ftx1 = fx1; fty1 = fy1; ftz1 = fz1
    ftx2 = fx2; fty2 = fy2; ftz2 = fz2
!DO n = 1, nk
!	if (n==50) then
!		write(*,*) 'time step', i, 'ft=', ft(n,:)
!		write(*,*) 'time step', i, 'ftx1=', ftx1(n,:)
!	endif
!enddo
    ! Transient calculation
    CALL Nodal_Coup4(P4,R4,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
	Call Nodal_Coup2(P2,R2,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)

DO n = 1, nk
		if (n>40 .and. n>43) then
	!	WRITE(*, *) "First Time Step   BEFORE outertf"
	!	WRITE(*, *) "step", i, " nu_sigf=", nu_sigf(n,:)
	!	WRITE(*, *) "step", i, " sigr=", sigr(n,:)
	!	WRITE(*, *) "step", i, " sigs=", sigs(n,:,:)
	!	WRITE(*, *) "step", i, " chi=", chi(:,:)
	!	WRITE(*, *) "step", i, " D=", D(n,:)
	!	WRITE(*, *) "step", i, " fs0=", fs0(n)
	!	WRITE(*, *) "step", i, " fsx1=", fsx1(n)
	!	WRITE(*, *) "step", i, " fsy1=", fsy1(n)
	!	WRITE(*, *) "step", i, " fsz1=", fsz1(n)
	!	WRITE(*, *) "step", i, " fsx2=", fsx2(n)
	!	WRITE(*, *) "step", i, " fsy2=", fsy2(n)
	!	WRITE(*, *) "step", i, " fsz2=", fsz2(n)
	!	WRITE(*, *) "step", i, " f0=", f0(n,:)
	!	WRITE(*, *) "step", i, " fx1=", fx1(n,:)
	!	WRITE(*, *) "step", i, " fy1=", fy1(n,:)
	!	WRITE(*, *) "step", i, " fz1=", fz1(n,:)
	!	WRITE(*, *) "step", i, " fx2=", fx2(n,:)
	!	WRITE(*, *) "step", i, " fy2=", fy2(n,:)
	!	WRITE(*, *) "step", i, " fz2=", fz2(n,:)
	!	WRITE(*, *) "step", i, " c0=", c0(n,:)
	!	WRITE(*, *) "step", i, " cx1=", cx1(n,:)
	!	WRITE(*, *) "step", i, " cy1=", cy1(n,:)
	!	WRITE(*, *) "step", i, " cz1=", cz1(n,:)
	!	WRITE(*, *) "step", i, " cx2=", cx2(n,:)
	!	WRITE(*, *) "step", i, " cy2=", cy2(n,:)
	!	WRITE(*, *) "step", i, " cz2=", cz2(n,:)
	!	WRITE(*, *) "step", i, " ft=", ft(n,:)
	!	WRITE(*, *) "step", i, " ftx1=", ftx1(n,:)
	!	WRITE(*, *) "step", i, " fty1=", fty1(n,:)
	!	WRITE(*, *) "step", i, " ftz1=", ftz1(n,:)
	!	WRITE(*, *) "step", i, " ftx2=", ftx2(n,:)
	!	WRITE(*, *) "step", i, " fty2=", fty2(n,:)
	!	WRITE(*, *) "step", i, " ftz2=", ftz2(n,:)
	!	WRITE(*, *) "step", i, " maxi=", maxi
	!	WRITE(*, *) "step", i, " R2=", (R2(n,g,6,6), g = 1, ng)
	!	WRITE(*, *) "step", i, " P2=", (P2(n,g,6,6), g = 1, ng)
	!	WRITE(*, *) "step", i, " R4=", (R4(n,g,6,6), g = 1, ng)
	!	WRITE(*, *) "step", i, " P4=", (P4(n,g,6,7), g = 1, ng)
	!	WRITE(*, *) "step", i, " mat=", mat(n)
	!	WRITE(*, *) "step", i, " al=", (al(n,g,6), g = 1, ng)
	!	WRITE(*, *) "step", i, " jo=", (jo(n,g,:), g = 1, ng)
	!	WRITE(*, *) "step", i, " ji=", (ji(n,g,:), g = 1, ng)
	!	WRITE(*, *) "step", i, " L0=", (L0(n,g,:), g = 1, ng)
	!	WRITE(*, *) "step", i, " delx=", delx(3)
	!	WRITE(*, *) "step", i, " dely=", dely(3)
	!	WRITE(*, *) "step", i, " delz=", delz(3)
	!	WRITE(*, *) "step", i, " delv=", delv(3)
	!	WRITE(*, *) "step", i, " xyz=", xyz(3,3,3)
	!	WRITE(*, *) "step", i, " iBeta=", iBeta(:)
	!	WRITE(*, *) "step", i, " lamb=", lamb(:)
	!	WRITE(*, *) "step", i, " tbeta=", tbeta
	!	WRITE(*, *) "step", i, " omeg=", omeg(n,:)
	!	WRITE(*, *) "step", i, " velo=", velo(:)
	!	WRITE(*, *) "step", i, " tstep1=", tstep1
	!	WRITE(*, *) "step", i, " errn=", errn(n)
	!	WRITE(*, *) "step", i, " erro=", erro(n)
	!	WRITE(*, *) "step", i, " fser=", fser
	!	WRITE(*, *) "step", i, " fer=", fer
	!	WRITE(*, *) "step", i, " Keff=", Keff
	endif
enddo

	!	CALL outertf (Keff,fer,fser,f0,fx1,fy1,fz1,fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,maxi,ng,nk,order,&
	!			  R2,P2,R4,P4,nmat,mat,al,jo,ji,L0,D,sigr,chi,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,&
	!			  delv,xyz,x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott,iBeta,lamb,&
	!			  tbeta,c0,cx1,cy1,cz1,cx2,cy2,cz2,ft,ftx1,fty1,ftz1,ftx2,fty2,ftz2,omeg,velo,tstep1,errn,erro)
!	WRITE(*, *) "Calling outertf"
!	WRITE(*, *) "Keo", Keo
!	WRITE(*, *) "fs0d", fs0d
!	WRITE(*, *) "f0d", f0d
!	WRITE(*, *) "fsd", fsd
do n=1, nk
		if (n>40 .and. n>43) then
		!	WRITE(*, *) "Q", (Q(n,g,:), g = 1, ng)
	endif
enddo

CALL outertf (Keff,fer,fser,f0,fx1,fy1,fz1,fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,maxi,ng,nk,order,&
		  R2,P2,R4,P4,nmat,mat,al,jo,ji,L0,D,sigr,chi,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,&
		  delv,xyz,x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott,iBeta,lamb,&
		  tbeta,c0,cx1,cy1,cz1,cx2,cy2,cz2,ft,ftx1,fty1,ftz1,ftx2,fty2,ftz2,omeg,velo,tstep1,errn,erro)

!CALL Outerfs(Keff,fer,fser,f0,fx1,fy1,fz1,fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,ng,nk,&
!			order,R2,P2,R4,P4,nmat,mat,al,jo,ji,L0,D,sigr,chi,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,&
!			Keo,fs0d,f0d,fsd,Q,delx,dely,delz,delv,xyz,x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott)
 DO n = 1, nk
	if (n>40 .and. n<43) then
		!WRITE(*, *) "First Time Step   AFTER outertf"
		!WRITE(*, *) "step", i, " fs0=", fs0(n)
		!WRITE(*, *) "step", i, " f0=", f0(n,:)
		!WRITE(*, *) "step", i, " af=", af(n,:)
		!WRITE(*, *) "step", i, " ft=", ft(n,:)
	endif
enddo
!Write(*,*) 'CALL outertf__1=' 
! DO n = 1, nk
!		if (n==50) then
!			write(*,*) 'f0=', f0(n,:)
!			write(*,*) 'fs0=', fs0(n)
!			write(*,*) 'af=', af(n,:)
!		endif
!	enddo

    ! Update precursor density
    CALL uPden(c0,cx1,cy1,cz1,cx2,cy2,cz2,nk,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,iBeta,lamb,tstep1)
! DO n = 1, nk
!	if (n==50) then
!		WRITE(*, *) "First Time Step   BEFORE PowTot"
!		WRITE(*, *) "step", i, " f0=", f0(n,:)
!		WRITE(*, *) "step", i, " sigf=", sigf(n,:)
!	endif
!enddo
    ! Calculate power
    CALL PowTot(tpow2,nk,ng,f0,sigf,delv)
! DO n = 1, nk
!		if (n==50) then
!			WRITE(*, *) "First Time Step   AFTER PowTot"
!			write(*,*) 'time step', i, 'tpow2=', tpow2
!		endif
!	enddo
! DO n = 1, nk
!	if (n==50) then
!		WRITE(*, *) "First Time Step   BEFORE react"
!		WRITE(*, *) "step", i, " sigrp=", sigrp(n,:)
!	endif
!enddo

!	WRITE(*, *) "First Time Step   BEFORE react"
! DO n = 1, nk
!	if (n==50) then
!		WRITE(*, *) "step", i, " fs0=", fs0(n)
!		WRITE(*, *) "step", i, " f0=", f0(n,:)
!		WRITE(*, *) "step", i, " af=", af(n,:)
!	endif
!enddo
    ! Calculate reactivity
    CALL react(rho,nk,ng,nmat,nxx,nyy,nzz,mat,sigs,chi,f0,fs0,delx,dely,delz,delv,ix,iy,iz,&
			   L0,af,sigrp)

	!WRITE(*, *) "First Time Step   AFTER react"
	!WRITE(*, *) "step", i, " RHO2=", rho
    WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 12F9.2)') step, t2, rho/tbeta, &
    tpow2/tpow1, (bpos_local(n), n = 1, nb)

    IF (maxi) THEN
        WRITE(*,'(I4, F10.3, F10.4, ES15.4, A35)') step, t2, rho/tbeta, &
        tpow2/tpow1, 'OUTER ITERATION DID NOT CONVERGE'
    ELSE
        WRITE(*,'(I4, F10.3, F10.4, ES15.4)') step, t2, rho/tbeta, &
        tpow2/tpow1
    END IF

    IF (maxi) tranw = .TRUE.

END DO



! Second Time Step
imax = NINT((ttot-tdiv)/tstep2)

DO i = 1, imax

    step = step + 1
    t1 = t2
    t2 = tdiv + REAL(i)*tstep2

    omeg = LOG(f0 / ft) / tstep2

    ! Rod bank changes
    DO n = 1, nb
        IF (mdir(n) == 1) THEN   ! If CR moving down
            IF (t2-tmove(n) > 1.e-5_DP .AND. fbpos(n)-bpos_local(n) < 1.e-5_DP) THEN
                bpos_local(n) = bpos_local(n) - tstep2 *  bspeed(n)
                IF (bpos_local(n) < fbpos(n)) bpos_local(n) = fbpos(n)  ! If bpos exceed, set to fbpos
            END IF
        ELSE IF (mdir(n) == 2) THEN ! If CR moving up
            IF (t2-tmove(n) > 1.e-5_DP .AND. fbpos(n)-bpos_local(n) > 1.e-5_DP) THEN
                bpos_local(n) = bpos_local(n) + tstep2 *  bspeed(n)
                IF (bpos_local(n) > fbpos(n)) bpos_local(n) = fbpos(n)  ! If bpos exceed, set to fbpos
            END IF
        ELSE
            CONTINUE
        END IF
     END DO

    ! Calculate xsec after perturbation
	CALL XS_updtrod(sigs,siga,sigtr,sigf,nu_sigf,D,sigr,ng,nk,nmat,mat,x_sigs,x_siga,x_sigtr,x_sigf,&
					x_nu_sigf,dsiga,dsigtr,dsigf,dnuf,dsigs,pos0,ssize,f0,fz1,fz2,jo,ji,nxx,nyy,nzz,&
					delz,xyz,nb,bpos_local,fbmap)
    ! Modify removal xsec
    sigrp = sigr    ! Save sigr to sigrp
    DO g = 1, ng
       DO n = 1, nk
          sigr(n,g) = sigr(n,g) + 1. / (velo(g) * tstep2) + omeg(n,g) / velo(g)
       END DO
    END DO

    ! Save the previous fluxes
    ft = f0
    ftx1 = fx1; fty1 = fy1; ftz1 = fz1
    ftx2 = fx2; fty2 = fy2; ftz2 = fz2

    ! Transient calculation
	CALL Nodal_Coup4(P4,R4,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
    Call Nodal_Coup2(P2,R2,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
 DO n = 1, nk
	if (n==50) then
	!	WRITE(*, *) "Second Time Step   Before outertf"
	!	WRITE(*, *) "step", i, " nu_sigf=", nu_sigf(n,:)
	!	WRITE(*, *) "step", i, " fs0=", fs0(n)
	!	WRITE(*, *) "step", i, " f0=", f0(n,:)
	!	WRITE(*, *) "step", i, " af=", af(n,:)
	endif
enddo
    CALL outertf (Keff,fer,fser,f0,fx1,fy1,fz1,fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,maxi,ng,nk,order,&
				  R2,P2,R4,P4,nmat,mat,al,jo,ji,L0,D,sigr,chi,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,&
				  delv,xyz,x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott,iBeta,lamb,&
				  tbeta,c0,cx1,cy1,cz1,cx2,cy2,cz2,ft,ftx1,fty1,ftz1,ftx2,fty2,ftz2,omeg,velo,tstep2,errn, erro)
 DO n = 1, nk
	if (n==50) then
	!	WRITE(*, *) "Second Time Step   AFTER outertf"
	!	WRITE(*, *) "step", i, " fs0=", fs0(n)
	!	WRITE(*, *) "step", i, " f0=", f0(n,:)
	!	WRITE(*, *) "step", i, " af=", af(n,:)
	endif
enddo
!Write(*,*) 'CALL outertf=' 
! DO n = 1, nk
!	if (n==50) then
!		write(*,*) 'f0=', f0(n,:)
!		write(*,*) 'fs0=', fs0(n)
!		write(*,*) 'af=', af(n,:)
!	endif
!enddo

				  ! Update precursor density
    CALL uPden(c0,cx1,cy1,cz1,cx2,cy2,cz2,nk,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,iBeta,lamb,tstep2)

    ! Calculate power
    CALL PowTot(tpow2,nk,ng,f0,sigf,delv)

    ! Calculate reactivity
    CALL react(rho,nk,ng,nmat,nxx,nyy,nzz,mat,sigs,chi,f0,fs0,delx,dely,delz,delv,ix,iy,iz,&
			   L0,af,sigrp)
! DO n = 1, nk
!	if (n==50) then
!		WRITE(*, *) "step", i, " fs0=", fs0(n)
!		WRITE(*, *) "step", i, " f0=", f0(n,:)
!		WRITE(*, *) "step", i, " af=", af(n,:)
!	endif
!enddo
	!	WRITE(*, *) "step", i, " RHO3=", rho

    WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 12F9.2)') step, t2, rho/tbeta, &
    tpow2/tpow1, (bpos_local(n), n = 1, nb)

    IF (maxi) THEN
        WRITE(*,'(I4, F10.3, F10.4, ES15.4, A35)') step, t2, rho/tbeta, &
        tpow2/tpow1, 'OUTER ITERATION DID NOT CONVERGE'
    ELSE
        WRITE(*,'(I4, F10.3, F10.4, ES15.4)') step, t2, rho/tbeta, &
        tpow2/tpow1
    END IF

    IF (maxi) tranw = .TRUE.

END DO

Call CPU_TIME(fn)
Write(ounit,*)
Write(*,*)
Write(ounit,*) "Total Time : ", fn-st, "Seconds"
Write(*,*) "Total Time : ", fn-st, "Seconds"

END SUBROUTINE

	SUBROUTINE throd_eject(ng,nk,nb,nx,ny,nz,nmat,order,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,chi,&
						   nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,x_smax,x_smin,xdiv,ydiv,zdiv,xsize,ysize,&
					       y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott,mode,&
					       dsiga,dsigtr,dsigf,dnuf,pos0,ssize,bpos,bmap,nstep,fbpos,tmove,bspeed,ibeta,&
						   lamb,velo,ttot,tstep1,tstep2,tdiv,&
						   ntem,stab,pra,cf,tin,cftem,cmtem,ccden,nfpin,cmflow,rf,tg,tc,ppitch,&
						   csiga,csigtr,csigf,cnuf,fsiga,fsigtr,fsigf,fnuf,msiga,msigtr,msigf,mnuf,&
						   lsiga,lsigtr,lsigf,lnuf,bcon,rbcon,rftem,rmtem,rcden,csigs,&
						   fsigs,msigs,lsigs,dsigs,powr,ppow)
!    To perform rod ejection simulation with TH feedbacks

Implicit None
Integer, Parameter :: dp = Selected_Real_Kind(10,15)
! #######################   User Input Paramaters   #######################
Integer, Intent(in) :: ng, nk, pra, nb, order, nx, ny, nz, nmat, nxx, nyy, nzz, x_east, x_west, y_north, y_south, z_top, z_bott
Integer, Dimension(nxx), Intent(in) :: x_smax, x_smin
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin
Integer, Dimension(nk), Intent(in) :: mat
Integer, Dimension(nx), Intent(in) :: xdiv
Integer, Dimension(ny), Intent(in) :: ydiv
Integer, Dimension(nz), Intent(in) :: zdiv
Real(dp), Dimension(nx), Intent(in) :: xsize           ! Assembly Size & Division
Real(dp), Dimension(ny), Intent(in) :: ysize           ! Assembly Size & Division
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk), Intent(in) :: delv, ix, iy, iz
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz
Character(Len=100), intent(in) :: mode

Real(dp), Dimension(nmat,ng), Intent(inout) :: chi, x_siga, x_sigtr, x_sigf, x_nu_sigf
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: x_sigs

REAL(DP), INTENT(IN) :: pos0, ssize			! Zero step position and step size

INTEGER, DIMENSION(nx, ny), intent(in) :: bmap       ! Radial control rod bank map (assembly wise)
REAL(DP), Intent(in) :: nstep		! Number of steps
! CXs Assigned to Materials
Real(dp), Dimension(nmat,ng), Intent(inout) :: csiga, csigtr, csigf, cnuf,&
fsiga, fsigtr, fsigf, fnuf, msiga, msigtr, msigf, mnuf, &
lsiga, lsigtr, lsigf, lnuf, dsiga, dsigtr, dsigf, dnuf
Real(dp), Intent(in) :: bcon, rbcon ,rftem, rmtem, rcden
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: csigs, fsigs, msigs, lsigs, dsigs
REAL(DP), Intent(in):: powr, ppow

! Crod changes
REAL(DP), DIMENSION(nb), intent(in) :: fbpos    ! Final CR bank position
REAL(DP), DIMENSION(nb), intent(in) :: tmove    ! Time when CR bank starts moving
REAL(DP), DIMENSION(nb), intent(in) :: bspeed   ! CR bank movement speed
REAL(DP), DIMENSION(nb), intent(in) :: bpos  ! CR bank position
INTEGER, DIMENSION(nb) :: mdir  ! To indicate CR movement direction (0=do not move, 1=down, 2 = up)
! Thermal-hydraulics parameters
INTEGER,  INTENT(IN) :: ntem    ! Number of temperature in steam table
REAL(DP), DIMENSION(ntem,pra), intent(in) :: stab  ! Steam table matrix
REAL(DP), intent(in) :: cf		! heat fraction deposited into coolant
REAL(DP), INTENT(IN) :: tin		! coolant inlet temperature (kelvin)
Real(dp), intent(in) :: cftem, cmtem, ccden
INTEGER, intent(in) :: nfpin		! Number of fuel pin and guide tubes
REAL(DP), intent(in) :: cmflow
REAL(DP), intent(in) :: rf, tg, tc, ppitch		! Fuel meat radius, gap thickness, clad thickness, and pin picth (m)
! Transient parameters 
INTEGER, PARAMETER :: nf = 6			! Number of delaye dneutron precusor family
REAL(DP), DIMENSION(nf), intent(in) :: ibeta, lamb                 ! beta (delayed neutron fraction) and precusor decay constant
REAL(DP) :: tbeta                                      ! total beta
REAL(DP), DIMENSION(ng), intent(in) :: velo            ! Neutron velocity
REAL(DP), intent(in) :: ttot                                       ! TOTAL SIMULATION TIME
REAL(DP), intent(in) :: tstep1                                     ! FIRST TIME STEP
REAL(DP), intent(in) :: tstep2                                     ! SECOND TIME STEP
REAL(DP), intent(in) :: tdiv                                       ! WHEN SECOND TIME STEP APPLY
										   !Frequency transformation constant
REAL(DP), DIMENSION(nk,ng) :: omeg         ! Exponential transformation constant
LOGICAL :: tranw = .FALSE.                        ! To activate unconverged  outer iteration warning
Real(dp), Dimension(nk) :: fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2     ! Fission Source & Fission Source Moments
REAL(DP), DIMENSION(nk,nf) :: c0, cx1, cy1, cz1, cx2, cy2, cz2  ! neutron precusor density
REAL(DP), DIMENSION(nk,ng) :: ft, ftx1, fty1, ftz1, ftx2, fty2, ftz2  ! Parameters at previous time step
REAL(dp), Dimension(nk) :: ftem, mtem, cden
INTEGER, PARAMETER :: nm = 10      ! Fuel meat divided into 10 mesh
INTEGER, PARAMETER :: nt = nm + 2      ! two more mesh for gap and clad
REAL(DP), DIMENSION(nt) :: rpos	! mesh position
REAL(DP), DIMENSION(nk,nt+1):: tfm	! Fuel pin mesh temperature for each nodes
INTEGER, DIMENSION(nxx,nyy) :: fbmap	! Radial control rod bank map (node wise)
! Local Variables
Real(dp) :: Keff   ! Effectif Multiplication Factor
Real(dp) :: fer, fser	! Flux and Fission Source Error in BCSEARCH calcs.
Real(dp), Dimension(nk,ng) :: f0, fx1, fy1, fz1, fx2, fy2, fz2      ! Flux and Flux Moments
Real(dp), Dimension(nk,ng,6) :: jo, ji, al  
Real(dp), Dimension(nk,ng,3) :: L0   ! Zeroth Transverse Leakages (Lx,Ly,Lz)
Real(dp), Dimension(nk,ng,6,6) :: R4, R2, P2    ! Response Matrix
Real(dp), Dimension(nk,ng,6,7) :: P4            ! Response Matrix
! Real(dp), Dimension(nk,ng,7) :: Q      ! Nodal's Source and Source Moments (0,x1,y1,z1,x2,y2,z2)
! CXs Assigned to Nodes
Real(dp), Dimension(nk,ng,ng) :: sigs  ! Scattering Macroscopic CX
Real(dp), Dimension(nk,ng) :: siga, sigtr, sigf, nu_sigf, D, sigr
REAL(DP), DIMENSION(nk):: power	! nodes power (watt)
Real(dp) :: ferc = 1.e-5	! Flux Error Criteria
Real(dp) :: fserc = 1.e-5	! Fission Source Error Criteria
Integer, Parameter :: Ounit = 100	! Output
REAL(DP), DIMENSION(nk, ng) :: af		! adjoint flux
REAL(DP), DIMENSION(nk, ng) :: sigrp	! Temporary sigr
REAL(DP) :: rho
REAL(DP) :: t1, t2
REAL(DP) :: tpow1, tpow2
INTEGER :: n, i, j, g, imax, step
LOGICAL :: maxi   ! Maximum Outer Iteration Reached?
REAL(DP), DIMENSION(nk) :: pline       ! Linear power density
REAL(DP) :: xppow
REAL(DP) :: tf, tm, mtf, mtm
REAL(DP), DIMENSION(nk) :: ent	! Coolant Enthalpy (J/Kg)

REAL(DP), DIMENSION(nxx, nyy) :: node_nf	! Number of fuel pin per node
Real(dp) :: st, fn

!##########################################################################
! Local Variables
REAL(DP) :: th_err		! Doppler error
INTEGER :: no = 10
REAL(DP), DIMENSION(nk) :: heatf	! Heat flux (W/m2)
REAL(DP) :: dum
REAL(DP), DIMENSION(nt) :: rdel	! mesh delta
REAL(DP) :: rg, rc	! Outer radius of gap and cladding
REAL(DP) :: dia	! Pi diameter, Hydraulic diameter (m) and sub-channel area (m2)
REAL(DP), PARAMETER :: pi = 3.14159265
REAL(DP) :: cflow
REAL(DP) :: dh, farea	! Pi diameter, Hydraulic diameter (m) and sub-channel area (m2)
Real(dp), Dimension(nk) :: errn, erro
real(dp), dimension(nb) :: bpos_local        ! Local copy of bpos

! Guess fuel and moderator temperature
! Used for radial fuel temperature distribution
tfm = 1200.
! Initial heat-flux rate
heatf = 0.

dum = rf / REAL(nm)
DO i = 1, nm
    rdel(i) = dum
END DO

! Fuel pin mesh size
rdel(nm+1) = tg
rdel(nm+2) = tc
! Fuel pin mesh position
rpos(1) = 0.5 * rdel(1)
DO i = 2, nt
    rpos(i) = rpos(i-1) + 0.5 * (rdel(i-1) + rdel(i))
END DO
!write(*,*) 'rpos=', rpos(:)
! Calculate outer radius of gap and cladding
rg = rf + tg
rc = rf + tg + tc
! Calculate fuel pin diameter
dia = 2. * rc
! Calculate hydraulic diameter
dh = dia * ((4./pi) * (ppitch/dia)**2 - 1.)
! Calculate sub-channel mass flow rate
cflow = cmflow / REAL(nfpin)
! Calculate sub-channel area
farea = ppitch**2 - 0.25*pi*dia**2

Call Init(keff,jo,ji,L0,al,f0,fx1,fy1,fz1,fx2,fy2,fz2,ng,nk)

Call Output_crod (fbmap,ng,nmat,nb,nx,ny,nxx,nyy,nzz,xdiv,ydiv,delz,nstep,bpos,&
				  pos0,ssize,dsiga,dsigtr,dsigf,dnuf,dsigs,bmap)
Call Output_ejct (mdir,ng,nk,nb,nmat,bpos,fbpos,tmove,bspeed,&
					ibeta,lamb,velo,ttot,tstep1,tstep2,tdiv)
Call Output_bcon (ng,nmat,bcon,rbcon,csigtr,csiga,cnuf,csigf,csigs)
Call Output_ther (node_nf,ng,nk,nx,ny,nmat,nxx,nyy,&
				  y_smax,y_smin,xdiv,xsize,ydiv,ysize,powr,ppow,tin,&
				  rf,tg,tc,ppitch,cf,pra,cmflow,nfpin,ntem,stab)
Call Output_ftem (ftem,ng,nk,nmat,cftem,rftem,fsigs,&
					 fsiga,fsigtr,fsigf,fnuf)
Call Output_mtem (mtem,ng,nk,nmat,cmtem,rmtem,msigs,&
					 msiga,msigtr,msigf,mnuf)
Call Output_cden (cden,ng,nk,nmat,ccden,rcden,lsigs,&
					 lsiga,lsigtr,lsigf,lnuf)

! Determine th paramters distribution
CALL th_iter(th_err,tfm,Keff,fer,fser,f0,sigf,nu_sigf,ent,ng,nk,nx,ny,nz,nmat,order,mat,x_sigs,x_siga,&
			 x_sigtr,x_sigf,x_nu_sigf,chi,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,&
			 x_smax,x_smin,xdiv,ydiv,zdiv,y_smax,y_smin,xyz,x_east,x_west,y_north,&
			 y_south,z_top,z_bott,mode,ntem,stab,cf,tin,ftem,mtem,cden,bcon,bpos,&
			 rf,tg,tc,ppitch,pra,jo,ji,L0,al,fx1,fx2,fy1,fy2,fz1,fz2,heatf,dh,farea,&
			 csiga,csigtr,csigf,cnuf,fsiga,fsigtr,fsigf,fnuf,msiga,msigtr,msigf,&
			 mnuf,lsiga,lsigtr,lsigf,lnuf,dsiga,dsigtr,dsigf,dnuf,rbcon,rftem,rmtem,&
			 rcden,csigs,fsigs,msigs,lsigs,dsigs,powr,ppow,node_nf,pos0,ssize,nb,fbmap,&
			 cflow,dia,rg,rc,rpos,rdel,fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2)
!Write(*,*) 'CALL th_iter=' 

		!	write(*,*) 'Keff=', Keff
 ! DO n = 1, nk
	!	if (n==50) then
		!	write(*,*) 'fbmap=', fbmap(5,5)
		!	write(*,*) 'csigf=', csigf(2,:)
		!	write(*,*) 'node_nf=', node_nf(5,5)
		!	write(*,*) 'ftem=', ftem(50)
		!	write(*,*) 'mtem=', mtem(50)
		!	write(*,*) 'cden=', cden(50)
		!	write(*,*) 'f0=', f0(n,:)
		!	write(*,*) 'tfm=', tfm(n,:)
	!	endif
	!enddo! If K-EFF NOT EQUAL TO 1.0

CALL XS_updtcfmd(sigs,siga,sigtr,sigf,nu_sigf,D,sigr,ng,nk,nmat,mat,x_sigs,x_siga,x_sigtr,x_sigf,&
			 x_nu_sigf,csiga,csigtr,csigf,cnuf,fsiga,fsigtr,fsigf,fnuf,msiga,msigtr,msigf,mnuf,&
			 lsiga,lsigtr,lsigf,lnuf,dsiga,dsigtr,dsigf,dnuf,csigs,fsigs,msigs,lsigs,dsigs,bcon,&
			 rbcon,rftem,rmtem,rcden,pos0,ssize,f0,fz1,fz2,jo,ji,nxx,nyy,nzz,&
			 delz,xyz,nb,ftem,mtem,cden,bpos,fbmap)
!    CALL XS_updtrod(sigs,siga,sigtr,sigf,nu_sigf,D,sigr,ng,nk,nmat,mat,x_sigs,x_siga,x_sigtr,x_sigf,&
!					x_nu_sigf,dsiga,dsigtr,dsigf,dnuf,dsigs,pos0,ssize,f0,fz1,fz2,jo,ji,nxx,nyy,nzz,&
!					delz,xyz,nb,bpos,fbmap)
! Update nodal couplings
CALL Nodal_coup4(P4,R4,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
Call Nodal_Coup2(P2,R2,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)

!Write(*,*) 'BEFORE CALL Outerfs=' 
DO n = 1, nk
	!if (n>40 .and. n<44) then
	!		write(*,*) 'n=', n, 'f0=', f0(n,:)
	!		write(*,*) 'n=', n, 'fx1=', fx1(n,:)
	!		write(*,*) 'n=', n, 'fs0=', fs0(n)
	!	endif
	enddo! If K-EFF NOT EQUAL TO 1.0

! Calculate forward flux
CALL Outerfs(Keff,fer,fser,f0,fx1,fy1,fz1,fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,ng,nk,&
			order,R2,P2,R4,P4,nmat,mat,al,jo,ji,L0,D,sigr,chi,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,&
			delx,dely,delz,delv,xyz,x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott)
!Write(*,*) 'AFTER CALL Outerfs=' 
!DO n = 1, nk
!	if (n>40 .and. n<44) then
!			write(*,*) 'n=', n, 'f0=', f0(n,:)
!			write(*,*) 'n=', n, 'fx1=', fx1(n,:)
!			write(*,*) 'n=', n, 'fs0=', fs0(n)
!		endif
!	enddo! If K-EFF NOT EQUAL TO 1.0

!Write(*,*) 'CALL Outerfs=' 
!DO n = 1, nk
!	if (n==50) then
!		write(*,*) 'f0=', f0(n,:)
!		write(*,*) 'fs0=', fs0(n)
!	endif
!enddo! If K-EFF NOT EQUAL TO 1.0
IF (ABS(Keff - 1._DP) > 1.e-5_DP) then
	!CALL KNE1(ng,nk,nmat,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,dsiga,dsigtr,dsigf,dnuf,dsigs,&
	!		pos0,ssize,nxx,nyy,nzz,delz,xyz,nb,bpos,fbmap,Keff,f0,fx1,fy1,fz1,fx2,fy2,fz2,&
	!		fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,order,R2,P2,R4,P4,al,jo,ji,L0,chi,ix,iy,iz,delx,dely,&
	!		delv,x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott)

	Call KNEth(ng,nk,nmat,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,dsiga,dsigtr,dsigf,dnuf,dsigs,&
			   pos0,ssize,nxx,nyy,nzz,delz,xyz,nb,bpos,fbmap,Keff,f0,fx1,fy1,fz1,fx2,fy2,fz2,&
			   fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,order,R2,P2,R4,P4,al,jo,ji,L0,chi,ix,iy,iz,delx,dely,&
			   delv,x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott,&
			   csiga,csigtr,csigf,cnuf,fsiga,fsigtr,fsigf,fnuf,&
			   msiga,msigtr,msigf,mnuf,lsiga,lsigtr,lsigf,lnuf,bcon,rbcon,rftem,rmtem,rcden,csigs,&
			   fsigs,msigs,lsigs,ftem,mtem,cden)
endif
!Write(*,*) 'CALL KNE1=' 
! DO n = 1, nk
!		if (n==50) then
!			write(*,*) 'f0=', f0(n,:)
!			write(*,*) 'fs0=', fs0(n)
!			write(*,*) 'Keff=', Keff
!			write(*,*) 'x_nu_sigf=', x_nu_sigf(5,:)
!			write(*,*) 'dnuf=', dnuf(5,:)
!		endif
!	enddo! If K-EFF NOT EQUAL TO 1.0


!Write(*,*) 'BEFORE CALL OuterAdjfs=' 	
!DO n = 1, nk
!		if (n==50) then
!			write(*,*) 'f0=', f0(n,:)
!			write(*,*) 'fs0=', fs0(n)
!			write(*,*) 'Keff=', Keff
!			write(*,*) 'nu_sigf=', nu_sigf(n,:)
!		endif
!	enddo! If K-EFF NOT EQUAL TO 1.0
!    CALL XS_updtrod(sigs,siga,sigtr,sigf,nu_sigf,D,sigr,ng,nk,nmat,mat,x_sigs,x_siga,x_sigtr,x_sigf,&
!					x_nu_sigf,dsiga,dsigtr,dsigf,dnuf,dsigs,pos0,ssize,f0,fz1,fz2,jo,ji,nxx,nyy,nzz,&
!					delz,xyz,nb,bpos,fbmap)

CALL XS_updtcfmd(sigs,siga,sigtr,sigf,nu_sigf,D,sigr,ng,nk,nmat,mat,x_sigs,x_siga,x_sigtr,x_sigf,&
			 x_nu_sigf,csiga,csigtr,csigf,cnuf,fsiga,fsigtr,fsigf,fnuf,msiga,msigtr,msigf,mnuf,&
			 lsiga,lsigtr,lsigf,lnuf,dsiga,dsigtr,dsigf,dnuf,csigs,fsigs,msigs,lsigs,dsigs,bcon,&
			 rbcon,rftem,rmtem,rcden,pos0,ssize,f0,fz1,fz2,jo,ji,nxx,nyy,nzz,&
			 delz,xyz,nb,ftem,mtem,cden,bpos,fbmap)
!Write(*,*) 'CALL OuterAdjfs=' 	
! Update nodal couplings
CALL Nodal_coup4(P4,R4,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
Call Nodal_Coup2(P2,R2,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
! Calculate Adjoint flux
CALL OuterAdjfs(Keff,f0,fx1,fy1,fz1,fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,ng,nk,order,&
				R2,P2,R4,P4,nmat,chi,mat,al,jo,ji,L0,D,sigr,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,&
				delx,dely,delz,delv,x_smax,x_smin,y_smax,y_smin,xyz,x_east,x_west,y_north,&
				y_south,z_top,z_bott)
af = f0   ! Save adjoint flux to af
! DO n = 1, nk
!		if (n==50) then
!			write(*,*) 'f0=', f0(n,:)
!			write(*,*) 'fs0=', fs0(n)
!			write(*,*) 'Keff=', Keff
!		endif
!	enddo! If K-EFF NOT EQUAL TO 1.0
	

! ReCalculate forward flux
CALL Outerfs(Keff,fer,fser,f0,fx1,fy1,fz1,fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,ng,nk,&
			order,R2,P2,R4,P4,nmat,mat,al,jo,ji,L0,D,sigr,chi,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,&
			delx,dely,delz,delv,xyz,x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott)

! Calculate power
CALL PowTot (tpow1,nk,ng,f0,sigf,delv)

! Calculate Initial precursor density
CALL iPden(c0,cx1,cy1,cz1,cx2,cy2,cz2,nk,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,iBeta,lamb)

! Total beta
tbeta = 0.
DO j = 1, nf
  tbeta = tbeta + iBeta(j)
END DO
!Write(*,*) 'BEFORE CALL react1='

! Calculate reactivity
CALL react(rho,nk,ng,nmat,nxx,nyy,nzz,mat,sigs,chi,f0,fs0,delx,dely,delz,delv,ix,iy,iz,&
		   L0,af,sigr)
!Write(*,*) 'AFTER CALL react1='
! ReCalculate forward flux
!DO n = 1, nk
!		if (n==50) then
!			write(*,*) 'f0=', f0(n,:)
!			write(*,*) 'fs0=', fs0(n)
!			write(*,*) 'sigr=', sigr(n,:)
!			write(*,*) 'rho=', rho
!		endif
!	enddo! If K-EFF NOT EQUAL TO 1.0

CALL par_ave_f(tf,ftem,nk,ng,nu_sigf,delv)
Call par_max(mtf,tfm(:,1),nk)
CALL par_ave(tm,mtem,nk,delv)
CALL par_max(mtm,mtem,nk)


! File output
WRITE(ounit, *)
WRITE(ounit, *) " TRANSIENT RESULTS :"
WRITE(ounit, *)
WRITE(ounit, *) " Step  Time(s)  React.($)   Rel. Power   Avg. Tm   Max. Tm   Avg. Tf   Max. Tf"
WRITE(ounit, *) "--------------------------------------------------------------------------------"
WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 4F10.2)') 0, 0., 0., &
ppow*0.01, tm-273.15, mtm-273.15, tf-273.15, mtf-273.15

! Terminal output
WRITE(*,*)
WRITE(*,*)
WRITE(*, *) " TRANSIENT RESULTS :"
WRITE(*, *)
WRITE(*, *) " Step  Time(s)  React.($)   Rel. Power"
WRITE(*, *) "-------------------------------------------"
WRITE(*,'(I4, F10.3, F10.4, ES15.4)') 0, 0., 0., ppow*0.01

! Start transient calculation
step = 0
t2 = 0.
imax = NINT(tdiv/tstep1)

! Initialize the local copy of bpos
bpos_local = bpos

! First Time Step
DO i = 1, imax

  step = step + 1
  t1 = t2
  t2 = REAL(i)*tstep1

  IF (i > 1) THEN
     omeg = LOG(f0 / ft) / tstep1
  ELSE
     omeg = 0.
  END IF


  ! Rod bank changes
  DO n = 1, nb
      IF (mdir(n) == 1) THEN   ! If CR moving down
          IF (t2-tmove(n) > 1.e-5_DP .AND. fbpos(n)-bpos_local(n) < 1.e-5_DP) THEN
              bpos_local(n) = bpos_local(n) - tstep1 *  bspeed(n)
              IF (bpos_local(n) < fbpos(n)) bpos_local(n) = fbpos(n)  ! If bpos exceed, set to fbpos
          END IF
      ELSE IF (mdir(n) == 2) THEN ! If CR moving up
          IF (t2-tmove(n) > 1.e-5_DP .AND. fbpos(n)-bpos_local(n) > 1.e-5_DP) THEN
              bpos_local(n) = bpos_local(n) + tstep1 *  bspeed(n)
              IF (bpos_local(n) > fbpos(n)) bpos_local(n) = fbpos(n)  ! If bpos exceed, set to fbpos
          END IF
      ELSE
          CONTINUE
      END IF
   END DO

  ! Calculate xsec after pertubation
   !CALL XS_updtrod(sigs,siga,sigtr,sigf,nu_sigf,D,sigr,ng,nk,nmat,mat,x_sigs,x_siga,x_sigtr,x_sigf,&
	!				x_nu_sigf,dsiga,dsigtr,dsigf,dnuf,dsigs,pos0,ssize,f0,fz1,fz2,jo,ji,nxx,nyy,nzz,&
	!				delz,xyz,nb,bpos,fbmap)
   CALL XS_updtcfmd(sigs,siga,sigtr,sigf,nu_sigf,D,sigr,ng,nk,nmat,mat,x_sigs,x_siga,x_sigtr,x_sigf,&
					 x_nu_sigf,csiga,csigtr,csigf,cnuf,fsiga,fsigtr,fsigf,fnuf,msiga,msigtr,msigf,mnuf,&
					 lsiga,lsigtr,lsigf,lnuf,dsiga,dsigtr,dsigf,dnuf,csigs,fsigs,msigs,lsigs,dsigs,bcon,&
					 rbcon,rftem,rmtem,rcden,pos0,ssize,f0,fz1,fz2,jo,ji,nxx,nyy,nzz,&
					 delz,xyz,nb,ftem,mtem,cden,bpos_local,fbmap)
  ! Modify removal xsec
  sigrp = sigr    ! Save sigr to sigrp
    DO g = 1, ng
       DO n = 1, nk
          sigr(n,g) = sigr(n,g) + 1. / (velo(g) * tstep1) + omeg(n,g) / velo(g)
       END DO
    END DO
    ! Save the previous fluxes
    ft = f0
    ftx1 = fx1; fty1 = fy1; ftz1 = fz1
    ftx2 = fx2; fty2 = fy2; ftz2 = fz2
 		!	write(*,*) 'BEFORE CALL outertf'
 DO n = 1, nk
	if (n>40 .and. n<44) then
		!	write(*,*) 'f0=', f0(n,:)
		!	write(*,*) 'fx1=', fx1(n,:)
		!	write(*,*) 'fs0=', fs0(n)
!			write(*,*) 'c0=', c0(n,:)
!			write(*,*) 'ft=', ft(n,:)
!!			write(*,*) 'tbeta=', tbeta
!!			write(*,*) 'omeg=', omeg(n,:)
!!			write(*,*) 'velo=', velo(:)
!!			write(*,*) 'lamb=', lamb(:)
!!			write(*,*) 'iBeta=', iBeta(:)
!!			write(*,*) 'tstep1=', tstep1
!			write(*,*) 'nu_sigf=', nu_sigf(n,:)
!			write(*,*) 'sigr=', sigr(n,:)
!			write(*,*) 'sigs=', sigs(n,:,:)
!!			write(*,*) 'P4=', P4(n,:,6,7)
		endif
	enddo
! If K-EFF NOT EQUAL TO 1.0
! DO g = 1, ng
!		if (n==50) then
!			write(*,*) 'ji=', P4(n,g,6,7)
!		endif
!enddo
! DO n = 1, nk
!		if (n==50) then
!			write(*,*) 'fs0=', fs0(n)
!		endif
!	enddo! If K-EFF NOT EQUAL TO 1.0
!maxi= .FALSE.
    ! Transient calculation
	CALL Nodal_Coup4(P4,R4,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
    Call Nodal_Coup2(P2,R2,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
DO n = 1, nk
		if (n>40 .and. n<43) then
	!	WRITE(*, *) "First Time Step   BEFORE outertf"
	!	WRITE(*, *) "step", i, 'node=', n,  " nu_sigf=", nu_sigf(n,:)
	!	WRITE(*, *) "step", i, 'node=', n,  " sigr=", sigr(n,:)
	!	WRITE(*, *) "step", i, 'node=', n,  " sigs=", sigs(n,:,:)
	!	WRITE(*, *) "step", i, 'node=', n,  " chi=", chi(:,:)
	!	WRITE(*, *) "step", i, 'node=', n,  " D=", D(n,:)
	!	WRITE(*, *) "step", i, 'node=', n,  " fs0=", fs0(n)
	!	WRITE(*, *) "step", i, 'node=', n,  " fsx1=", fsx1(n)
	!	WRITE(*, *) "step", i, 'node=', n,  " fsy1=", fsy1(n)
	!	WRITE(*, *) "step", i, 'node=', n,  " fsz1=", fsz1(n)
	!	WRITE(*, *) "step", i, 'node=', n,  " fsx2=", fsx2(n)
	!	WRITE(*, *) "step", i, 'node=', n,  " fsy2=", fsy2(n)
	!	WRITE(*, *) "step", i, 'node=', n,  " fsz2=", fsz2(n)
	!	WRITE(*, *) "step", i, 'node=', n,  " f0=", f0(n,:)
	!	WRITE(*, *) "step", i, 'node=', n,  " fx1=", fx1(n,:)
	!	WRITE(*, *) "step", i, 'node=', n,  " fy1=", fy1(n,:)
	!	WRITE(*, *) "step", i, 'node=', n,  " fz1=", fz1(n,:)
	!	WRITE(*, *) "step", i, 'node=', n,  " fx2=", fx2(n,:)
	!	WRITE(*, *) "step", i, 'node=', n,  " fy2=", fy2(n,:)
	!	WRITE(*, *) "step", i, 'node=', n,  " fz2=", fz2(n,:)
	!	WRITE(*, *) "step", i, 'node=', n,  " c0=", c0(n,:)
	!	WRITE(*, *) "step", i, 'node=', n,  " cx1=", cx1(n,:)
	!	WRITE(*, *) "step", i, 'node=', n,  " cy1=", cy1(n,:)
	!	WRITE(*, *) "step", i, 'node=', n,  " cz1=", cz1(n,:)
	!	WRITE(*, *) "step", i, 'node=', n,  " cx2=", cx2(n,:)
	!	WRITE(*, *) "step", i, 'node=', n,  " cy2=", cy2(n,:)
	!	WRITE(*, *) "step", i, 'node=', n,  " cz2=", cz2(n,:)
	!	WRITE(*, *) "step", i, 'node=', n,  " ft=", ft(n,:)
	!	WRITE(*, *) "step", i, 'node=', n,  " ftx1=", ftx1(n,:)
	!	WRITE(*, *) "step", i, 'node=', n,  " fty1=", fty1(n,:)
	!	WRITE(*, *) "step", i, 'node=', n,  " ftz1=", ftz1(n,:)
	!	WRITE(*, *) "step", i, 'node=', n,  " ftx2=", ftx2(n,:)
	!	WRITE(*, *) "step", i, 'node=', n,  " fty2=", fty2(n,:)
	!	WRITE(*, *) "step", i, 'node=', n,  " ftz2=", ftz2(n,:)
	!	WRITE(*, *) "step", i, 'node=', n,  " maxi=", maxi
	!	WRITE(*, *) "step", i, 'node=', n,  " R2=", (R2(n,g,6,6), g = 1, ng)
	!	WRITE(*, *) "step", i, 'node=', n,  " P2=", (P2(n,g,6,6), g = 1, ng)
	!	WRITE(*, *) "step", i, 'node=', n,  " R4=", (R4(n,g,6,6), g = 1, ng)
	!	WRITE(*, *) "step", i, 'node=', n,  " P4=", (P4(n,g,6,7), g = 1, ng)
	!	WRITE(*, *) "step", i, 'node=', n,  " mat=", mat(n)
	!	WRITE(*, *) "step", i, 'node=', n,  " al=", (al(n,g,6), g = 1, ng)
	!	WRITE(*, *) "step", i, 'node=', n,  " jo=", (jo(n,g,:), g = 1, ng)
	!	WRITE(*, *) "step", i, 'node=', n,  " ji=", (ji(n,g,:), g = 1, ng)
	!	WRITE(*, *) "step", i, 'node=', n,  " L0=", (L0(n,g,:), g = 1, ng)
	!	WRITE(*, *) "step", i, 'node=', n,  " delx=", delx(3)
	!	WRITE(*, *) "step", i, 'node=', n,  " dely=", dely(3)
	!	WRITE(*, *) "step", i, 'node=', n,  " delz=", delz(3)
	!	WRITE(*, *) "step", i, 'node=', n,  " delv=", delv(3)
	!	WRITE(*, *) "step", i, 'node=', n,  " xyz=", xyz(3,3,3)
	!	WRITE(*, *) "step", i, 'node=', n,  " iBeta=", iBeta(:)
	!	WRITE(*, *) "step", i, 'node=', n,  " lamb=", lamb(:)
	!	WRITE(*, *) "step", i, 'node=', n,  " tbeta=", tbeta
	!	WRITE(*, *) "step", i, 'node=', n,  " omeg=", omeg(n,:)
	!	WRITE(*, *) "step", i, 'node=', n,  " velo=", velo(:)
	!	WRITE(*, *) "step", i, 'node=', n,  " tstep1=", tstep1
	!	WRITE(*, *) "step", i, 'node=', n,  " errn=", errn(n)
	!	WRITE(*, *) "step", i, 'node=', n,  " erro=", erro(n)
	!	WRITE(*, *) "step", i, 'node=', n,  " fser=", fser
	!	WRITE(*, *) "step", i, 'node=', n,  " fer=", fer
	!	WRITE(*, *) "step", i, 'node=', n,  " Keff=", Keff
	endif
enddo
	CALL outertf (Keff,fer,fser,f0,fx1,fy1,fz1,fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,maxi,ng,nk,order,&
				  R2,P2,R4,P4,nmat,mat,al,jo,ji,L0,D,sigr,chi,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,&
				  delv,xyz,x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott,iBeta,lamb,&
				  tbeta,c0,cx1,cy1,cz1,cx2,cy2,cz2,ft,ftx1,fty1,ftz1,ftx2,fty2,ftz2,omeg,velo,tstep1,errn, erro)

!CALL Outerfs(Keff,fer,fser,f0,fx1,fy1,fz1,fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,ng,nk,&
!			order,R2,P2,R4,P4,nmat,mat,al,jo,ji,L0,D,sigr,chi,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,&
!			delx,dely,delz,delv,xyz,x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott)


DO n = 1, nk
		if (n>40 .and. n<43) then
		!WRITE(*, *) "First Time Step   AFTER outertf"
		!WRITE(*, *) "step", i, 'node=', n,  " fs0=", fs0(n)
		!WRITE(*, *) "step", i, 'node=', n,  " fsx1=", fsx1(n)
		!WRITE(*, *) "step", i, 'node=', n,  " fsy1=", fsy1(n)
		!WRITE(*, *) "step", i, 'node=', n,  " fsz1=", fsz1(n)
		!WRITE(*, *) "step", i, 'node=', n,  " fsx2=", fsx2(n)
		!WRITE(*, *) "step", i, 'node=', n,  " fsy2=", fsy2(n)
		!WRITE(*, *) "step", i, 'node=', n,  " fsz2=", fsz2(n)
		!WRITE(*, *) "step", i, 'node=', n,  " f0=", f0(n,:)
		!WRITE(*, *) "step", i, 'node=', n,  " fx1=", fx1(n,:)
		!WRITE(*, *) "step", i, 'node=', n,  " fy1=", fy1(n,:)
		!WRITE(*, *) "step", i, 'node=', n,  " fz1=", fz1(n,:)
		!WRITE(*, *) "step", i, 'node=', n,  " fx2=", fx2(n,:)
		!WRITE(*, *) "step", i, 'node=', n,  " fy2=", fy2(n,:)
		!WRITE(*, *) "step", i, 'node=', n,  " fz2=", fz2(n,:)
		!WRITE(*, *) "step", i, 'node=', n,  " maxi=", maxi
		!WRITE(*, *) "step", i, 'node=', n,  " errn=", errn(n)
		!WRITE(*, *) "step", i, 'node=', n,  " erro=", erro(n)
		!WRITE(*, *) "step", i, 'node=', n,  " fser=", fser
		!WRITE(*, *) "step", i, 'node=', n,  " fer=", fer
		!WRITE(*, *) "step", i, 'node=', n,  " Keff=", Keff
	endif
enddo
 		!	write(*,*) 'AFTER CALL outertf'
 DO n = 1, nk
	if (n>40 .and. n<44) then
		!	write(*,*) 'f0=', f0(n,:)
		!	write(*,*) 'fx1=', fx1(n,:)
		!	write(*,*) 'fs0=', fs0(n)
		endif
	enddo! If K-EFF NOT EQUAL TO 1.0

    ! Update precursor density
    CALL uPden(c0,cx1,cy1,cz1,cx2,cy2,cz2,nk,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,iBeta,lamb,tstep1)

    ! Calculate power
	CALL PowTot (tpow2,nk,ng,f0,sigf,delv)

    ! Calculate node power distribution
    Call PowDis(power,f0,sigf,delv,ng,nk,mode)

    ! Power change
    xppow = ppow * tpow2/tpow1 * 0.01
				   
    ! Calculate linear power density for each nodes (W/cm)
    DO n = 1, nk
       pline(n) = power(n) * powr * xppow &
       / (node_nf(ix(n),iy(n)) * delz(iz(n)))     ! Linear power density (W/cm)
    END DO
!WRITE(*,*) 'i=', i, 'pline=',pline(50)

    ! TH transient
   CALL th_trans(ftem,tfm,heatf,nk,nxx,nyy,nzz,pline,tstep1,ix,iy,iz,delz,ntem,stab,cf,&
				  tin,ent,rf,tg,tc,ppitch,pra,rpos,rdel,rg,rc,dia,cflow,dh,farea,mtem,cden)
!do n=1,nk
!			if (n == 30) then
!		!		WRITE(*,*) '   ftem=',ftem(n)
!			end if
!enddo
    ! Calculate reactivity
    CALL react(rho,nk,ng,nmat,nxx,nyy,nzz,mat,sigs,chi,f0,fs0,delx,dely,delz,delv,ix,iy,iz,&
			   L0,af,sigrp)
!Write(*,*) 'CALL react2='
! ReCalculate forward flux
! DO n = 1, nk
!		if (n==50) then
!			write(*,*) 'f0=', f0(n,:)
!			write(*,*) 'fs0=', fs0(n)
!			write(*,*) 'rho=', rho
!		endif
!	enddo! If K-EFF NOT EQUAL TO 1.0

	CALL par_ave_f(tf,ftem,nk,ng,nu_sigf,delv)
	Call par_max(mtf,tfm(:,1),nk)
	CALL par_ave(tm,mtem,nk,delv)
	CALL par_max(mtm,mtem,nk)

    WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 4F10.2)') step, t2, rho/tbeta, &
    xppow, tm-273.15, mtm-273.15, tf-273.15, mtf-273.15

    IF (maxi) THEN
        WRITE(*,'(I4, F10.3, F10.4, ES15.4, A35)') step, t2, rho/tbeta, &
        xppow, 'OUTER ITERATION DID NOT CONVERGE'
    ELSE
        WRITE(*,'(I4, F10.3, F10.4, ES15.4)') step, t2, rho/tbeta, &
        xppow
    END IF

    IF (maxi) tranw = .TRUE.

END DO

!Second Time Step
imax = NINT((ttot-tdiv)/tstep2)

DO i = 1, imax

  step = step + 1
  t1 = t2
  t2 = tdiv + REAL(i)*tstep2

  omeg = LOG(f0 / ft) / tstep2

  ! Rod bank changes
  DO n = 1, nb
      IF (mdir(n) == 1) THEN   ! If CR moving down
          IF (t2-tmove(n) > 1.e-5_DP .AND. fbpos(n)-bpos_local(n) < 1.e-5_DP) THEN
              bpos_local(n) = bpos_local(n) - tstep2 *  bspeed(n)
              IF (bpos_local(n) < fbpos(n)) bpos_local(n) = fbpos(n)  ! If bpos exceed, set to fbpos
          END IF
      ELSE IF (mdir(n) == 2) THEN ! If CR moving up
          IF (t2-tmove(n) > 1.e-5_DP .AND. fbpos(n)-bpos_local(n) > 1.e-5_DP) THEN
              bpos_local(n) = bpos_local(n) + tstep2 *  bspeed(n)
              IF (bpos_local(n) > fbpos(n)) bpos_local(n) = fbpos(n)  ! If bpos exceed, set to fbpos
          END IF
      ELSE
          CONTINUE
      END IF
   END DO

  ! Calculate xsec after pertubation
   CALL XS_updtcfmd(sigs,siga,sigtr,sigf,nu_sigf,D,sigr,ng,nk,nmat,mat,x_sigs,x_siga,x_sigtr,x_sigf,&
					 x_nu_sigf,csiga,csigtr,csigf,cnuf,fsiga,fsigtr,fsigf,fnuf,msiga,msigtr,msigf,mnuf,&
					 lsiga,lsigtr,lsigf,lnuf,dsiga,dsigtr,dsigf,dnuf,csigs,fsigs,msigs,lsigs,dsigs,bcon,&
					 rbcon,rftem,rmtem,rcden,pos0,ssize,f0,fz1,fz2,jo,ji,nxx,nyy,nzz,&
					 delz,xyz,nb,ftem,mtem,cden,bpos_local,fbmap)
					 
  ! Modify removal xsec
    sigrp = sigr    ! Save sigr to sigrp
    DO g = 1, ng
       DO n = 1, nk
          sigr(n,g) = sigr(n,g) + 1. / (velo(g) * tstep2) + omeg(n,g) / velo(g)
       END DO
    END DO

    ! Save the previous fluxes
    ft = f0
    ftx1 = fx1; fty1 = fy1; ftz1 = fz1
    ftx2 = fx2; fty2 = fy2; ftz2 = fz2

    ! Transient calculation
	CALL Nodal_Coup4(P4,R4,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
    Call Nodal_Coup2(P2,R2,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
    CALL outertf(Keff,fer,fser,f0,fx1,fy1,fz1,fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,maxi,ng,nk,order,&
				  R2,P2,R4,P4,nmat,mat,al,jo,ji,L0,D,sigr,chi,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,&
				  delv,xyz,x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott,iBeta,lamb,&
				  tbeta,c0,cx1,cy1,cz1,cx2,cy2,cz2,ft,ftx1,fty1,ftz1,ftx2,fty2,ftz2,omeg,velo,tstep2,errn, erro)


    ! Update precursor density
    CALL uPden(c0,cx1,cy1,cz1,cx2,cy2,cz2,nk,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,iBeta,lamb,tstep2)

    ! Calculate power
    CALL PowTot(tpow2,nk,ng,f0,sigf,delv)

    ! Calculate node power distribution
    Call PowDis(power,f0,sigf,delv,ng,nk,mode)


    ! Power change
    xppow = ppow * tpow2/tpow1 * 0.01
				   
    ! Calculate linear power density for each nodes (W/cm)
    DO n = 1, nk
       pline(n) = power(n) * powr * xppow &
       / (node_nf(ix(n),iy(n)) * delz(iz(n)))     ! Linear power density (W/cm)
    END DO

    ! TH transient
    CALL th_trans(ftem,tfm,heatf,nk,nxx,nyy,nzz,pline,tstep2,ix,iy,iz,delz,ntem,stab,cf,&
				  tin,ent,rf,tg,tc,ppitch,pra,rpos,rdel,rg,rc,dia,cflow,dh,farea,mtem,cden)
    ! Calculate reactivity
    CALL react(rho,nk,ng,nmat,nxx,nyy,nzz,mat,sigs,chi,f0,fs0,delx,dely,delz,delv,ix,iy,iz,&
			   L0,af,sigrp)
!Write(*,*) 'CALL react3='
! ReCalculate forward flux
!DO n = 1, nk
!		if (n==50) then
!			write(*,*) 'f0=', f0(n,:)
!			write(*,*) 'fs0=', fs0(n)
!			write(*,*) 'rho=', rho
!		endif
!	enddo! If K-EFF NOT EQUAL TO 1.0

	CALL par_ave_f(tf,ftem,nk,ng,nu_sigf,delv)
	Call par_max(mtf,tfm(:,1),nk)
	CALL par_ave(tm,mtem,nk,delv)
	CALL par_max(mtm,mtem,nk)

    WRITE(ounit,'(I4, F10.3, F10.4, ES15.4, 4F10.2)') step, t2, rho/tbeta, &
    xppow, tm-273.15, mtm-273.15, tf-273.15, mtf-273.15

    IF (maxi) THEN
        WRITE(*,'(I4, F10.3, F10.4, ES15.4, A35)') step, t2, rho/tbeta, &
        xppow, 'OUTER ITERATION DID NOT CONVERGE'
    ELSE
        WRITE(*,'(I4, F10.3, F10.4, ES15.4)') step, t2, rho/tbeta, &
        xppow
    END IF

    IF (maxi) tranw = .TRUE.

END DO

	CALL PowDis(power,f0,sigf,delv,ng,nk,mode)

    Call AsmPow(nk,nx,ny,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,xdiv,ydiv,y_smax,y_smin,power)

    Call AxiPow(nk,nz,nxx,nyy,nzz,ix,iy,iz,xyz,delz,delv,zdiv,y_smax,y_smin,power)

    Call AsmFlux(ng,nk,nx,ny,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,xdiv,ydiv,y_smax,y_smin,f0,1.e0_DP)

Call CPU_TIME(fn)
Write(ounit,*)
Write(*,*)
Write(ounit,*) "Total Time : ", fn-st, "Seconds"
Write(*,*) "Total Time : ", fn-st, "Seconds"

END SUBROUTINE 

		SUBROUTINE th_iter(th_err,tfm,Keff,fer,fser,f0,sigf,nu_sigf,ent,ng,nk,nx,ny,nz,nmat,order,mat,x_sigs,x_siga,&
				x_sigtr,x_sigf,x_nu_sigf,chi,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,&
				x_smax,x_smin,xdiv,ydiv,zdiv,y_smax,y_smin,xyz,x_east,x_west,y_north,&
				y_south,z_top,z_bott,mode,ntem,stab,cf,tin,ftem,mtem,cden,bcon,bpos,&
				rf,tg,tc,ppitch,pra,jo,ji,L0,al,fx1,fx2,fy1,fy2,fz1,fz2,heatf,dh,farea,&
				csiga,csigtr,csigf,cnuf,fsiga,fsigtr,fsigf,fnuf,msiga,msigtr,msigf,&
				mnuf,lsiga,lsigtr,lsigf,lnuf,dsiga,dsigtr,dsigf,dnuf,rbcon,rftem,rmtem,&
				rcden,csigs,fsigs,msigs,lsigs,dsigs,powr,ppow,node_nf,pos0,ssize,nb,fbmap,&
				cflow,dia,rg,rc,rpos,rdel,fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2)
!    To do thermal-hydraulics iteration
Implicit None
Integer, Parameter :: dp = Selected_Real_Kind(10,15)
! #######################   User Input Paramaters   #######################
Integer, Intent(in) :: ng, nk, pra,nb, order, nx, ny, nz, nmat, nxx, nyy, nzz, x_east, x_west, y_north, y_south, z_top, z_bott
Integer, Dimension(nxx), Intent(in) :: x_smax, x_smin
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin
Integer, Dimension(nk), Intent(in) :: mat
Integer, Dimension(nx), Intent(in) :: xdiv
Integer, Dimension(ny), Intent(in) :: ydiv
Integer, Dimension(nz), Intent(in) :: zdiv
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk), Intent(in) :: delv, ix, iy, iz
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz
Character(Len=100), intent(in) :: mode
INTEGER,  INTENT(IN) :: ntem    ! Number of temperature in steam table
REAL(DP), DIMENSION(ntem,pra), intent(in) :: stab  ! Steam table matrix
Real(dp), Dimension(nk,ng), intent(inout) :: fx1, fy1, fz1, fx2, fy2, fz2      ! Flux and Flux Moments
Real(dp), Dimension(nk), intent(inout) :: fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2     ! Fission Source & Fission Source Moments

! Thermal-hydraulics parameters
REAL(DP), intent(in) :: cf	! heat fraction deposited into coolant
REAL(DP), INTENT(IN) :: tin	! coolant inlet temperature (kelvin)
Real(dp), Dimension(nk), intent(inout) :: ftem, mtem, cden
Real(dp), Dimension(nb), intent(in) :: bpos
REAL(DP), intent(in) :: rf, tg, tc, ppitch	! Fuel meat radius, gap thickness, clad thickness, and pin picth (m)
INTEGER, PARAMETER :: nm = 10      ! Fuel meat divided into 10 mesh
INTEGER, PARAMETER :: nt = nm + 2      ! two more mesh for gap and clad
REAL(DP), DIMENSION(nt), intent(in) :: rdel	! mesh delta
REAL(DP), DIMENSION(nt), intent(in) :: rpos	! mesh position

! CXs Assigned to Materials
Real(dp), Dimension(nmat,ng), Intent(in) :: chi, x_siga, x_sigtr, x_sigf, x_nu_sigf
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: x_sigs
Real(dp), Dimension(nmat,ng), Intent(in) :: csiga, csigtr, csigf, cnuf,&
fsiga, fsigtr, fsigf, fnuf, msiga, msigtr, msigf, mnuf, &
lsiga, lsigtr, lsigf, lnuf, dsiga, dsigtr, dsigf, dnuf
Real(dp), Intent(in) :: rbcon ,rftem, rmtem, rcden
INTEGER, DIMENSION(nxx,nyy), intent(in) :: fbmap	! Radial control rod bank map (node wise)

Real(dp), Dimension(nmat,ng,ng), Intent(in) :: csigs, fsigs, msigs, lsigs, dsigs
REAL(DP), Intent(in):: powr, ppow
REAL(DP), DIMENSION(nxx, nyy), intent(in) :: node_nf       ! Number of fuel pin per node
REAL(DP), INTENT(IN) :: pos0, ssize			! Zero step position and step size
REAL(DP), intent(out) :: th_err                                     ! Doppler error
Real(dp), intent(inout) :: Keff   ! Effectif Multiplication Factor
Real(dp), intent(out) :: fer, fser                ! Flux and Fission Source Error in BCSEARCH calcs.
Real(dp), Dimension(nk,ng), intent(out) :: sigf     ! Fission Macroscopic CX
Real(dp), Dimension(nk,ng), intent(out)  :: nu_sigf  ! nu * Fission Macroscopic CX
Real(dp), Dimension(nk,ng), intent(inout) :: f0
REAL(DP), DIMENSION(nk,nt+1), intent(inout):: tfm	! Fuel pin mesh temperature for each nodes
REAL(DP), DIMENSION(nk), intent(inout) :: heatf	! Heat flux (W/m2)
REAL(DP), intent(in) :: dh, farea	! Pi diameter, Hydraulic diameter (m) and sub-channel area (m2)
REAL(DP), DIMENSION(nk), intent(out) :: ent	! Coolant Enthalpy (J/Kg)

!##########################################################################
! Local Variables
Real(dp), Dimension(nk,ng,6), intent(inout) :: jo   ! Nodals' Outgoing Currents (X+,X-,Y+,Y-,Z+,Z-)
Real(dp), Dimension(nk,ng,6), intent(inout) :: ji   ! Nodals' Ingoing Currents  (X+,X-,Y+,Y-,Z+,Z-)
Real(dp), Dimension(nk,ng,3), intent(inout) :: L0   ! Zeroth Transverse Leakages (Lx,Ly,Lz)
Real(dp), Dimension(nk,ng,6), intent(inout) :: al   ! Assembly Discontinuity Factor
Real(dp), Dimension(nk,ng,6,6) :: R4, R2, P2    ! Response Matrix
Real(dp), Dimension(nk,ng,6,7) :: P4            ! Response Matrix
! Real(dp), Dimension(nk,ng,7) :: Q      ! Nodal's Source and Source Moments (0,x1,y1,z1,x2,y2,z2)
! CXs Assigned to Nodes
Real(dp), Dimension(nk,ng,ng) :: sigs  ! Scattering Macroscopic CX
Real(dp), Dimension(nk,ng) :: siga     ! Absorption Macroscopic CX
Real(dp), Dimension(nk,ng) :: sigtr    ! Transport Macroscopic CX
Real(dp), Dimension(nk,ng) :: D        ! Diffusion Coefficient
Real(dp), Dimension(nk,ng) :: sigr     ! Scattering Macroscopic CX
REAL(DP), DIMENSION(nk) :: pline
REAL(DP), DIMENSION(nk) :: otem
Real(dp), Dimension(nk) :: power
INTEGER :: n, l, k
INTEGER :: th_niter = 20                           ! Maximum number of thermal-hydraulics iteration
!INTEGER :: th_niter = 20                           ! Maximum number of thermal-hydraulics iteration
Real(dp), intent(in) :: bcon
Integer, Parameter :: Ounit = 100	! Output
Integer, Parameter :: Oun = 101	! Output
REAL(DP) :: dum

REAL(DP), intent(in) :: dia	! Pi diameter, Hydraulic diameter (m) and sub-channel area (m2)
!REAL(DP), PARAMETER :: pi = 3.14159265
REAL(DP), intent(in) :: cflow
REAL(DP), intent(in) :: rg, rc	! Outer radius of gap and cladding
INTEGER :: i

!Call Init(keff,jo,ji,L0,al,f0,fx1,fy1,fz1,fx2,fy2,fz2,ng,nk)
Open (Unit=Oun, File='app/sh.txt')

	th_err = 1.
!tfm = 1200.
! Calculate fuel pin mesh delta and position


  DO l = 1, th_niter

      ! Save old fuel temp
      otem = ftem
	! if (l == 2) THEN 
	!	WRITE(*,*) 'itr=',l , 'BEFORE'
	!	do k=1,nk
	!		if (k>65 .AND. k<85) then
	!			WRITE(*,*) 'itr=',l , '   tfm=',tfm(k,1)
	!			WRITE(*,*) 'n=', k, '   sigf=',sigf(k,:)
	!		end if
	!	enddo
	! ! exit
	! end if 

	!			WRITE(*,*) 'BEFORE',  'iteration =',  l
    !  DO n = 1, nk
!	!  ftem(n) = 558.55
	!		IF (n>65 .AND. n<85) then 
	!			WRITE(*,*) 'n=', n, '   siga=',siga(n, :)
	!		end if
    !  END DO
	  
      ! Update XS
      CALL XS_updtcfmd(sigs,siga,sigtr,sigf,nu_sigf,D,sigr,ng,nk,nmat,mat,x_sigs,x_siga,x_sigtr,x_sigf,&
					 x_nu_sigf,csiga,csigtr,csigf,cnuf,fsiga,fsigtr,fsigf,fnuf,msiga,msigtr,msigf,mnuf,&
					 lsiga,lsigtr,lsigf,lnuf,dsiga,dsigtr,dsigf,dnuf,csigs,fsigs,msigs,lsigs,dsigs,bcon,&
					 rbcon,rftem,rmtem,rcden,pos0,ssize,f0,fz1,fz2,jo,ji,nxx,nyy,nzz,&
					 delz,xyz,nb,ftem,mtem,cden,bpos,fbmap)
					 

      ! Update nodal couplings
      CALL nodal_coup4(P4,R4,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
      Call Nodal_Coup2(P2,R2,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)


!Write(*,*) 'BEFORE CALL outer4th=' 
DO n = 1, nk
	!if (n>40 .and. n<44) then
	!		write(*,*) 'iter=', l, 'n=', n, 'f0=', f0(n,:)
	!		write(*,*) 'iter=', l, 'n=', n, 'fx1=', fx1(n,:)
	!		write(*,*) 'iter=', l, 'n=', n, 'fs0=', fs0(n)
	!	endif
	enddo! If K-EFF NOT EQUAL TO 1.0

      ! Perform outer inner iteration
      CALL outer4th(Keff,fer,fser,f0,fx1,fy1,fz1,fx2,fy2,fz2,20,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,&
					ng,nk,order,R2,P2,R4,P4,nmat,mat,al,jo,ji,L0,&
                    D,sigr,chi,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,xyz,&
                    x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott)
!CALL Outerfs(Keff,fer,fser,f0,fx1,fy1,fz1,fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,ng,nk,&
!			order,R2,P2,R4,P4,nmat,mat,al,jo,ji,L0,D,sigr,chi,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,&
!			delx,dely,delz,delv,xyz,x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott)

!Write(*,*) 'AFTER CALL outer4th=' 
DO n = 1, nk
	!if (n>40 .and. n<44) then
	!		write(*,*) 'iter=', l, 'n=', n, 'f0=', f0(n,:)
	!		write(*,*) 'iter=', l, 'n=', n, 'fx1=', fx1(n,:)
	!		write(*,*) 'iter=', l, 'n=', n, 'fs0=', fs0(n)
	!	endif
	enddo! If K-EFF NOT EQUAL TO 1.0


!WRITE(*,*) 'itr=',l , '   keff=',Keff
	!	if (l == 1) THEN
	!	WRITE(*,*) 'itr=',l , '   keff=',Keff
!DO n = 1, nk
!	IF (n>65 .AND. n<85) then 
!		WRITE(*,*) '   f0=',f0(n,:)
!	endif
!enddo
	!	end if 

      ! Calculate power density
      CALL PowDis(power,f0,sigf,delv,ng,nk,mode)
      ! Calculate linear power density for each nodes (W/cm)


	!			WRITE(*,*) 'BEFORE',  'iteration =',  l
    !  DO n = 1, nk
!	!  ftem(n) = 558.55
	!		IF (n>65 .AND. n<85) then 
	!			WRITE(*,*) 'n=', n, '   pline=',pline(n)
	!		end if
    !  END DO
    !  DO n = 1, nk
!	!  ftem(n) = 558.55
	!		IF (n>60 .AND. n<80) then 
	!			WRITE(*,*) 'n=', n, '   pline=',pline(n)
	!		end if
    !  END DO
	  
	  
      DO n = 1, nk
          pline(n) = power(n) * powr * ppow * 0.01 &
                   / (node_nf(ix(n),iy(n)) * delz(iz(n)))     ! Linear power density (W/cm)
	!!	if (n == 66) then
	!!		WRITE(*,*) 'itr=',l
	!!		write(*,*) 'pline=', pline(n)
	!!	end if
	!	!IF (n<150) then 
	!	!WRITE(*,*) 'n=', n, 'pline=',pline(n)
	!	!endif
  !
      END DO

	!			WRITE(*,*) 'AFTER',  'iteration =',  l
    !  DO n = 1, nk
!	!  ftem(n) = 558.55
	!		IF (n>65 .AND. n<85) then 
	!			WRITE(*,*) 'n=', n, '   pline=',pline(n)
	!		end if
    !  END DO
	  
	  
	  
      ! Update fuel, moderator temp. and coolant density
	!	WRITE(*,*) 'itr=',l
!WRITE(*,*) 'BEFORE TH UPDATE'
    !  DO n = 1, nk
		!	IF (n>900 .AND. n<950) then 
		!		WRITE(*,*) 'n:', n, 'tfm=',tfm(n,:)
		!		WRITE(*,*) 'n:', n, 'ftem=',ftem(n)
		!	end if
      !END DO
	  
	 ! DO n = 1, nk
	 ! Read(Oun,*) ftem(n)
	  !ftem(n) = 559.14999999999964
	 ! END DO

	  !DO n = 500, nk
	 ! Read(Oun,*) ftem(n)
	  !ftem(n) = 559.15000000000458
	  !END DO
      CALL th_upd(ftem,tfm,heatf,ent,mtem,cden,nk,nxx,nyy,nzz,pline,ix,iy,iz,delz,ntem,stab,cf,&
				  tin,rf,tg,tc,ppitch,pra,l,rpos,rdel,rg,rc,dia,cflow,dh,farea)

!Write(*,*) 'AFTER CALL TH UPDATE=' 
DO n = 1, nk
	!if (n>40 .and. n<44) then
	!		write(*,*) 'iter=', l, 'n=', n, 'ftem=', ftem(n)
	!		write(*,*) 'iter=', l, 'n=', n, 'mtem=', mtem(n)
	!		write(*,*) 'iter=', l, 'n=', n, 'cden=', cden(n)
	!		write(*,*) 'iter=', l, 'n=', n, 'tfm=', tfm(n,:)
	!		write(*,*) 'iter=', l, 'n=', n, 'heatf=', heatf(n)
	!		write(*,*) 'iter=', l, 'n=', n, 'ent=', ent(n)
	!		write(*,*) 'iter=', l, 'n=', n, 'pline=', pline(n)
	!	endif
	enddo! If K-EFF NOT EQUAL TO 1.0


!ftem = 558.55
!WRITE(*,*) 'AFTER TH UPDATE'
	!			WRITE(*,*) 'AFTER'
    !  DO n = 1, nk
!	!  ftem(n) = 558.55
	!		IF (n>65 .AND. n<85) then 
	!			WRITE(*,*) 'n=', n, '   ftem=',ftem(n)
	!		end if
    !  END DO
	!  
	  
			CALL AbsE(ftem, otem, th_err, nk)
      IF (th_err < 0.01) EXIT

  END DO
     IF (l >= 20) THEN
        WRITE(ounit,*) '  MAXIMUM TH ITERATION REACHED.'
        WRITE(ounit,*) '  CALCULATION MIGHT BE NOT CONVERGED OR CHANGE ITERATION CONTROL'
        STOP
     END IF

END SUBROUTINE

			SUBROUTINE th_upd(ftem,tfm,heatf,ent,mtem,cden,nk,nxx,nyy,nzz,xpline,ix,iy,iz,delz,ntem,stab,cf,&
				  tin,rf,tg,tc,ppitch,pra,l,rpos,rdel,rg,rc,dia,cflow,dh,farea)
				  !    To update thermal parameters
IMPLICIT NONE
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: nk, nxx, nyy, nzz, pra
REAL(DP), DIMENSION(nk), INTENT(IN) :: xpline    ! Linear Power Density (W/cm)
Real(dp), Dimension(nk),  Intent(in) :: ix, iy, iz
Real(dp), Dimension(nzz), Intent(in) :: delz
INTEGER,  INTENT(IN) :: ntem    ! Number of temperature in steam table
REAL(DP), DIMENSION(ntem,pra), intent(in) :: stab  ! Steam table matrix
! Thermal-hydraulics parameters
REAL(DP), intent(in) :: cf	! heat fraction deposited into coolant
REAL(DP), INTENT(IN) :: tin	! coolant inlet temperature (kelvin)
REAL(DP) :: Pr, kv, tcon ! Coolant Prandtl Number, Kinematic viscosity, and thermal conductivity
INTEGER, intent(in) ::  l	! Number of fuel pin and guide tubes
REAL(DP), intent(in) :: rf, tg, tc, ppitch	! Fuel meat radius, gap thickness, clad thickness, and pin picth (m)
INTEGER, PARAMETER :: nm = 10      ! Fuel meat divided into 10 mesh
INTEGER, PARAMETER :: nt = nm + 2      ! two more mesh for gap and clad
REAL(DP), DIMENSION(nk,nt+1), intent(inout) :: tfm	! Fuel pin mesh temperature for each nodes
REAL(DP), DIMENSION(nt), intent(in) :: rdel	! mesh delta
REAL(DP), DIMENSION(nt), intent(in) :: rpos	! mesh position
REAL(DP), DIMENSION(nk), intent(inout) :: heatf	! Heat flux (W/m2)
REAL(DP), DIMENSION(nk), intent(inout) :: ftem	! Fuel temperature in Kelvin for each nodes
Real(dp), Dimension(nk), intent(inout) :: mtem, cden
REAL(DP), intent(in) :: dh, farea	! Pi diameter, Hydraulic diameter (m) and sub-channel area (m2)

REAL(DP), DIMENSION(nk), intent(inout) :: ent	! Coolant Enthalpy (J/Kg)
REAL(DP), intent(in) :: dia	! Pi diameter, Hydraulic diameter (m) and sub-channel area (m2)
REAL(DP), PARAMETER :: pi = 3.14159265
REAL(DP), intent(in) :: cflow
REAL(DP), intent(in) :: rg, rc	! Outer radius of gap and cladding
!INTEGER :: nt = nm + 2	! Number Total mesh
INTEGER :: i, n
REAL(DP), DIMENSION(nt+1) :: a, b, c, d
REAL(DP) :: hs, Hg = 1.d4, kt, kt1, kt2
REAL(DP) :: alp = 0.7_DP, xa, xc
REAL(DP) :: pdens      ! power densisty  (W/m3)
REAL(DP) :: enti       ! Coolant inlet enthalpy
REAL(DP), DIMENSION(nxx, nyy) :: entm
REAL(DP) :: cpline     ! Coolant Linear power densisty (W/m)
!ent=858341.9
!ftem=891.19000000000005
!mtem=559.19000000000005 
!cden=0.74639999999999995
!tfm=1200.00
!heatf=0.0000


! Guess fuel and moderator temperature
! Used for radial fuel temperature distribution
!tfm = 1200.
!mtem = 0.0
CALL getent(tin, enti, ntem, stab,pra)
	!	WRITE(*,*) 'enti=',enti

DO n = 1, nk
   cpline = heatf(n) * pi * dia  + cf * xpline(n) * 100._DP	! Coolant Linear power densisty (W/m)
		!	if (n == 65) then
		!		WRITE(*,*)'Before Calling TridiaSolve'
		!		WRITE(*,*)'node=', n,  '   heatf=',heatf(n)
		!	end if
!IF (n<150) then 
!WRITE(*,*) 'n=', n, 'xpline=',xpline(n)
!endif

   IF (iz(n) == 1) THEN										! For most bootom channel
      ent(n) = enti + 0.5_DP * cpline * delz(iz(n)) * 0.01_DP / cflow	! Calculate coolant enthalpy
	  Call gettd(ent(n), mtem(n), cden(n), Pr, kv, tcon, ntem, stab,pra)  ! Get corresponding temp and density
      entm(ix(n),iy(n)) = 2._DP * ent(n) - enti                      ! Extrapolate enthalpy at node boundary
   ELSE
      ent(n) = entm(ix(n),iy(n)) + 0.5_DP * cpline * delz(iz(n)) * 0.01_DP / cflow
      CALL gettd(ent(n), mtem(n), cden(n), Pr, kv, tcon, ntem, stab,pra)
      entm(ix(n),iy(n)) = 2._DP * ent(n) - entm(ix(n),iy(n))
   END IF
   
		!	IF (n>65 .AND. n<85) then 
		!		WRITE(*,*) 'n=', n, '   mtem=',mtem(n)
		!		WRITE(*,*) 'n=', n, '   cden=',cden(n)
		!	endif 

	Call geths(hs,cden(n),Pr,kv,tcon,&
			   cflow,farea,dh)
   pdens = (1. - cf) * 100._DP * xpline(n) / (pi * rf**2)        ! Fuel pin Power Density (W/m3)

   ! Calculate tridiagonal matrix: a, b, c and source: d
   ! For nt=1 [FUEL CENTERLINE]

   Call getkf(kt1,tfm(n,1))			! Get thermal conductivity
   Call getkf(kt2,tfm(n,2))			! Get thermal conductivity

   kt  = 2._DP * kt1 * kt2 / (kt1 + kt2)
   xc  = kt * rpos(1) / rdel(1)
   b(1) =  xc
   c(1) = -xc
   d(1) = pdens * 0.5_DP * rpos(1)**2

   DO i = 2, nt-2
      kt1 = kt2
      Call getkf(kt2,tfm(n,i+1))
      kt  = 2._DP * kt1 * kt2 / (kt1 + kt2)
      xa = xc
      xc = kt * rpos(i) / rdel(i)
      a(i) = -xa
      b(i) =  xa + xc
      c(i) = -xc
      d(i) = pdens * 0.5_DP * (rpos(i)**2 - rpos(i-1)**2)
   END DO

   ! For nt-1 [FUEL-GAP INTERFACE]
   xa = xc
   xc = rg * Hg
   a(nt-1) = -xa
   b(nt-1) =  xa + xc
   c(nt-1) = -xc
   d(nt-1) = pdens * 0.5_DP * (rf**2 - rpos(nt-2)**2)

   ! For nt [GAP-CLADDING INTERFACE]
   Call getkc(kt1,tfm(n,nt))
   Call getkc(kt2,tfm(n,nt+1))

   kt  = 2._DP * kt1 * kt2 / (kt1 + kt2)     ! For cladding
   xa = xc
   xc = kt * rpos(nt) / rdel(nt)
   a(nt) = -xa
   b(nt) =  xa + xc
   c(nt) = -xc
   d(nt) = 0.

   ! For nt+1  [CLADDING-COOLANT INTERFACE]
   xa = xc
   a(nt+1) = -xa
   b(nt+1) =  xa + hs * rc
   d(nt+1) = rc * hs * mtem(n)

		!	if (n == 65) then
		!		WRITE(*,*)'Before Calling TridiaSolve'
		!		WRITE(*,*)'node=', n,  '   tfm=',tfm(n,:)
		!	end if

   ! Solve tridiagonal matrix
   CALL TridiaSolve(a,b,c,d,tfm(n, :),nk)


	  !if (l == 2) THEN 
	!	WRITE(*,*) 'itr=',l , '   ftem=',ftem
	!	do k=1,nk
		!	if (n == 66) then
		!		WRITE(*,*)'node=', n,  '   tfm=',tfm(n,:)
		!	end if
	!	enddo
	  ! exit
	 ! end if 
	
	
   ! Get lumped fuel temp
   ftem(n) = (1.-alp) * tfm(n,1) + alp * tfm(n,nt-1)

!			if (n == 30) then
!				WRITE(*,*) '   ftem=',ftem(n)
!			end if

	
   ! Calculate heat flux
   heatf(n) = hs * (tfm(n,nt+1) - mtem(n))
		!	if (n == 65) then
		!		WRITE(*,*)'After Calling TridiaSolve'
		!		WRITE(*,*)'node=', n,  '   heatf=',heatf(n)
		!	end if
!IF (n<70) then 
!write(*,*) 'n=', n, 'ftem=', heatf(n)
!endif
 !  			if (n == 30) then
	!			WRITE(*,*) '   heatf=',heatf(n)
	!		end if

END DO


!write(*,*) 'enti:', enti
!write(*,*) 'cpline:', cpline
!write(*,*) 'cflow:', cflow
!do n=1,nk
!	IF (n>70 .AND. n<250) then 
!		WRITE(*,*) '   heatf=',heatf(n)
!	endif
!enddo

!do n=1,nk
!	IF (n>70 .AND. n<350) then 
!		write(*,*) 'ent:', ent(n)
!	endif
!enddo
END SUBROUTINE

		SUBROUTINE th_trans(ftem,tfm,heatf,nk,nxx,nyy,nzz,xpline,h,ix,iy,iz,delz,ntem,stab,cf,&
						tin,ent,rf,tg,tc,ppitch,pra,rpos,rdel,rg,rc,dia,cflow,dh,farea,mtem,cden)
!    To perform fuel pin thermal transient

IMPLICIT NONE
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: nk, nxx, nyy, nzz, pra
REAL(DP), DIMENSION(nk), INTENT(IN) :: xpline    ! Linear Power Density (W/cm)
Real(dp), Dimension(nk),  Intent(in) :: ix, iy, iz
Real(dp), Dimension(nzz), Intent(in) :: delz
INTEGER,  INTENT(IN) :: ntem    ! Number of temperature in steam table
REAL(DP), DIMENSION(ntem,pra), intent(in) :: stab  ! Steam table matrix
! Thermal-hydraulics parameters
REAL(DP), intent(in) :: cf	! heat fraction deposited into coolant
REAL(DP), INTENT(IN) :: tin	! coolant inlet temperature (kelvin)
REAL(DP) :: Pr, kv, tcon ! Coolant Prandtl Number, Kinematic viscosity, and thermal conductivity
REAL(DP), intent(in) :: rf, tg, tc, ppitch	! Fuel meat radius, gap thickness, clad thickness, and pin picth (m)
INTEGER, PARAMETER :: nm = 10      ! Fuel meat divided into 10 mesh
INTEGER, PARAMETER :: nt = nm + 2      ! two more mesh for gap and clad
REAL(DP), DIMENSION(nk,nt+1), intent(out) :: tfm	! Fuel pin mesh temperature for each nodes
REAL(DP), DIMENSION(nt), intent(in) :: rdel	! mesh delta
REAL(DP), DIMENSION(nt), intent(in) :: rpos	! mesh position
REAL(DP), DIMENSION(nk), intent(inout) :: heatf	! Heat flux (W/m2)
REAL(DP), DIMENSION(nk), intent(out) :: ftem	! Fuel temperature in Kelvin for each nodes
Real(dp), Dimension(nk), intent(out) :: mtem, cden
REAL(DP), intent(in) :: dh, farea	! Pi diameter, Hydraulic diameter (m) and sub-channel area (m2)
REAL(DP), INTENT(IN) :: h			! Time step

REAL(DP), DIMENSION(nk), intent(inout) :: ent	! Coolant Enthalpy (J/Kg)
REAL(DP), intent(in) :: dia	! Pi diameter, Hydraulic diameter (m) and sub-channel area (m2)
REAL(DP), PARAMETER :: pi = 3.14159265
REAL(DP), intent(in) :: cflow
REAL(DP), intent(in) :: rg, rc	! Outer radius of gap and cladding
!INTEGER :: nt = nm + 2	! Number Total mesh
INTEGER :: i, n
REAL(DP), DIMENSION(nt+1) :: a, b, c, d
REAL(DP) :: hs, Hg = 1.d4, kt, kt1, kt2
REAL(DP) :: alpha = 0.7_DP
REAL(DP) :: xa, xc
REAL(DP) :: fdens = 10.412e3		! UO2 density (kg/m3)
REAL(DP) :: cdens = 6.6e3			! Cladding density (kg/m3)
REAL(DP) :: cp						! Specific heat capacity
REAL(DP) :: eps, eta
REAL(DP) :: pdens      				! power densisty  (W/m3)
REAL(DP) :: enti       				! Coolant inlet enthalpy
REAL(DP), DIMENSION(nxx, nyy) :: entm
REAL(DP) :: cpline     				! Coolant Linear power densisty (W/m)
REAL(DP), DIMENSION(nk) :: entp     ! previous enthalpy
REAL(DP) :: mdens, vol				! Coolant density and channel volume


CALL getent(tin, enti, ntem, stab,pra)
entp = ent
!write(*,*) 'tin=', tin
!write(*,*) 'enti=', enti

DO n = 1, nk
   mdens = cden(n) * 1000._DP                                    ! Coolant density (kg/m3)
   cpline = heatf(n) * pi * dia + cf * xpline(n) * 100._DP       ! Coolant Linear power densisty (W/m)
   vol   = farea * delz(iz(n)) * 0.01_DP                             ! channel node volume

   IF (iz(n) == 1) THEN                                              ! Calculate coolant enthalpy
       eps = mdens * vol / h
       ent(n) = (cpline * delz(iz(n)) * 0.01_DP + 2._DP * cflow * enti &
              + eps * entp(n)) / (eps + 2._DP * cflow)                             ! Calculate enthalpy
       	  Call gettd(ent(n), mtem(n), cden(n), Pr, kv, tcon, ntem, stab,pra)  ! Get corresponding temp and density        ! Get corresponding temp and density
       entm(ix(n),iy(n)) = 2._DP * ent(n) - enti
   ELSE
       eps = mdens * vol / h
       ent(n) = (cpline * delz(iz(n)) * 0.01_DP + 2._DP * cflow * entm(ix(n),iy(n)) &
              + eps * entp(n)) / (eps + 2._DP * cflow)
      CALL gettd(ent(n), mtem(n), cden(n), Pr, kv, tcon, ntem, stab,pra)
       entm(ix(n),iy(n)) = 2._DP * ent(n) - entm(ix(n),iy(n))
   END IF
if (n==50) then
	!write(*,*) 'mdens=', mdens
	!write(*,*) 'cpline=', cpline
	!write(*,*) 'delz=', delz(10)
	!write(*,*) 'cflow=', cflow
	!write(*,*) 'eps=', eps
	!write(*,*) 'entp=', entp(n)
	!write(*,*) '============'
	!write(*,*) 'ent=', ent(n)
	!write(*,*) 'mtem=', mtem(n)
	!write(*,*) 'cden=', cden(n)
endif
	Call geths(hs,cden(n),Pr,kv,tcon,&
			   cflow,farea,dh)
   pdens = (1. - cf) * 100._DP * xpline(n) / (pi * rf**2)                ! Fuel pin Power Density (W/m3)

   ! Calculate tridiagonal matrix: a, b, c and source: d
   ! For nt=1 [FUEL CENTERLINE]
   Call getkf(kt1,tfm(n,1))			! Get thermal conductivity
   Call getkf(kt2,tfm(n,2))			! Get thermal conductivity
   kt  = 2._DP * kt1 * kt2 / (kt1 + kt2)
   Call getcpf(cp,tfm(n,1))
   eta = fdens * cp * rpos(1)**2 / (2._DP * h)
   xc  = kt * rpos(1) / rdel(1)
   b(1) =  xc + eta
   c(1) = -xc
   d(1) = pdens * 0.5_DP * rpos(1)**2 + eta * tfm(n,1)

   DO i = 2, nt-2
       kt1 = kt2
	   Call getkf(kt2,tfm(n,i+1))			! Get thermal conductivity
       kt  = 2._DP * kt1 * kt2 / (kt1 + kt2)
	   Call getcpf(cp,tfm(n,i))
       eta = fdens * cp * (rpos(i)**2 - rpos(i-1)**2) / (2. * h)
       xa = xc
       xc = kt * rpos(i) / rdel(i)
       a(i) = -xa
       b(i) =  xa + xc + eta
       c(i) = -xc
       d(i) = pdens * 0.5_DP * (rpos(i)**2 - rpos(i-1)**2) + eta * tfm(n,i)
   END DO

   ! For nt-1 [FUEL-GAP INTERFACE]
   Call getcpf(cp,tfm(n,nt-1))
   eta = fdens * cp * (rf**2 - rpos(nt-2)**2) / (2. * h)
   xa = xc
   xc = rg * hg
   a(nt-1) = -xa
   b(nt-1) =  xa + xc + eta
   c(nt-1) = -xc
   d(nt-1) = pdens * 0.5_DP * (rf**2 - rpos(nt-2)**2) + eta * tfm(n,nt-1)

   ! For nt [GAP-CLADDING INTERFACE]
   Call getkc(kt1,tfm(n,nt))
   Call getkc(kt2,tfm(n,nt+1))
   kt  = 2._DP * kt1 * kt2 / (kt1 + kt2)     ! For cladding
   Call getcpc(cp,tfm(n,nt))
   eta = cdens * cp * (rpos(nt)**2 - rg**2) / (2. * h)
   xa = xc
   xc = kt * rpos(nt) / rdel(nt)
   a(nt) = -xa
   b(nt) =  xa + xc + eta
   c(nt) = -xc
   d(nt) = eta * tfm(n,nt)

   ! For nt+1  [CLADDING-COOLANT INTERFACE]
   Call getcpc(cp,tfm(n,nt+1))
   eta = cdens * cp * (rc**2 - rpos(nt)**2) / (2. * h)
   xa = xc
   xc = rc * hs
   a(nt+1) = -xa
   b(nt+1) =  xa + xc + eta
   d(nt+1) = rc * hs * mtem(n) + eta * tfm(n,nt+1)

   ! Solve tridiagonal matrix
   CALL TridiaSolve(a,b,c,d,tfm(n, :),nk)

   ! Get lumped fuel temp
   ftem(n) = (1.-alpha) * tfm(n, 1) + alpha * tfm(n, nt-1)

   ! Calculate heat flux
   heatf(n) = hs * (tfm(n, nt+1) - mtem(n))
END DO

END SUBROUTINE

        Subroutine Outer(Keff,fer,fser,f0,fx1,fy1,fz1,fx2,fy2,fz2,ng,nk,order,R2,P2,R4,P4,nmat,mat,al,jo,ji,L0,&
                         D,sigr,chi,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,xyz,&
                         x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott)    ! To Perform Outer Iterations
Implicit None
Integer, Parameter :: dp = Selected_Real_Kind(10,15), Ounit = 100
Integer, Intent(in) :: ng, nk, order, nmat, nxx, nyy, nzz, x_east, x_west, y_north, y_south, z_top, z_bott
Integer, Dimension(nk), Intent(in) :: mat
Integer, Dimension(nxx), Intent(in) :: x_smax, x_smin
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin 
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk), Intent(in) :: delv, ix, iy, iz
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz
Real(dp), Dimension(nk,ng,ng), Intent(in) :: sigs
Real(dp), Dimension(nk,ng), Intent(in) :: nu_sigf, D, sigr
Real(dp), Dimension(nmat,ng), Intent(in) :: chi
Real(dp), Dimension(nk,ng,6), Intent(inout):: al, jo, ji
Real(dp), Dimension(nk,ng,3), Intent(inout) :: L0
Real(dp), Dimension(nk,ng,6,6),Intent(in) :: R2, P2, R4
Real(dp), Dimension(nk,ng,6,7),Intent(in) :: P4
Real(dp), Dimension(nk,ng), intent(inout) :: f0, fx1, fy1, fz1, fx2, fy2, fz2
Real(dp), Intent(inout) :: Keff
! Local Variables
Integer :: nout = 500                ! Maximum Outer Iteration
Integer :: nac = 5                   ! Number of Outer Iteration Before Next Source Extrapolation
Real(dp) :: ferc = 1.e-5             ! Flux Error Criteria
Real(dp) :: fserc = 1.e-5            ! Fission Source Error Criteria
Real(dp), intent(out) :: fer, fser	! Flux and Fission Source Error in BCSEARCH calcs.
Real(dp), Dimension(nk) :: fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2     ! Fission Source & Fission Source Moments
Real(dp) :: Keo                      ! Old Multiplication Factor (Keff)
Real(dp), Dimension(nk) :: fs0d      ! Old Fission Source
Real(dp), Dimension(nk,ng) :: f0d    ! Old Flux
Real(dp) :: fsd                      ! Old Integrated Fission Sources
Real(dp) :: fs                       ! New Integrated Fission Sources
Real(dp) :: domiR, e1, e2
Integer :: g, o, npos, n
Real(dp), Dimension(nk,ng,7) :: Q
Real(dp), Dimension(nk) :: errn, erro

!#####################    Initialize Fission Source    ####################
    Call FSrc(fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,&
              f0,fx1,fy1,fz1,fx2,fy2,fz2,nu_sigf,ng,nk)

errn = 1._dp          ! Variable Fission Source Error
    Call Integrate(fs,fs0,delv,nk)
    Call Integrate(e1,errn,delv,nk)

!#######################    Start Outer Iteration   #######################
    Do o=1, nout
        f0d = f0          ! Save Old Fluxes
        fsd = fs          ! Save Old Integrated Fission Source
        fs0d  = fs0       ! Save Old Fission Source
        Keo = Keff        ! Save Old Multiplication Fctor
        erro = errn       ! Save Old Fission Source Error/Difference
 
        Do g = 1, ng
!#######################   Calculate Total Source   #######################
!				do n=1, nk
!					if (n==50) then
!							write(*,*) 'BEFORE CALL TSrc'
!							write(*,*) 'group=', g, 'f0=', f0(n,:)
!							write(*,*) 'group=', g, 'fx1=', fx1(n,:)
!							write(*,*) 'group=', g, 'fy1=', fy1(n,:)
!							write(*,*) 'group=', g, 'fz1=', fz1(n,:)
!							write(*,*) 'group=', g, 'fx2=', fx2(n,:)
!							write(*,*) 'group=', g, 'fy2=', fy2(n,:)
!							write(*,*) 'group=', g, 'fz2=', fz2(n,:)
!							write(*,*) 'group=', g, 'fs0=', fs0(n)
!							write(*,*) 'group=', g, 'fsx1=', fsx1(n)
!							write(*,*) 'group=', g, 'fsy1=', fsy1(n)
!							write(*,*) 'group=', g, 'fsz1=', fsz1(n)
!							write(*,*) 'group=', g, 'fsx2=', fsx2(n)
!							write(*,*) 'group=', g, 'fsy2=', fsy2(n)
!							write(*,*) 'group=', g, 'fsz2=', fsz2(n)
!							write(*,*) 'group=', g, 'Keff=', Keo
!					endif
!				enddo
            Call TSrc(Q,g,Keo,ng,nmat,nk,chi,mat,sigs,f0,fx1,fy1,fz1,&
                      fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2)
			!if (p==250) then
!				do n=1, nk
!					if (n==50) then
!							write(*,*) 'group=', g, 'source=', Q(n,g,:)
!							write(*,*) 'BEFORE CALL INNER'
!							write(*,*) 'group=', g, 'f0=', f0(n,:)
!					endif
!				enddo
			!endif
!#######################    Start Inner Iteration   #######################
            Call Inner(f0,fx1,fy1,fz1,fx2,fy2,fz2,g,ng,nk,order,jo,ji,L0,Q,R2,P2,R4,P4,al,D,sigr,nxx,nyy,nzz,ix,iy,iz,delx,&
                       dely,delz,x_smax,x_smin,y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott)
		!	if (p==250) then
			!	do n=1, nk
			!		if (n==50) then
			!				write(*,*) 'AFTER CALL INNER'
			!				write(*,*) 'group=', g, 'f0=', f0(n,:)
			!		endif
			!	enddo
			!endif
		End Do

!########################   Update Fission Source   #######################
        Call FSrc(fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,&
                  f0,fx1,fy1,fz1,fx2,fy2,fz2,nu_sigf,ng,nk)
        errn = fs0 - fs0d
        Call Integrate(e2,ABS(errn),delv,nk)    ! Variable Integrated Fission Source Error

!####################   Fission Source Extrapolation   ####################
        if (MOD(o,nac) == 0) Then
            domiR = e2 / e1
            npos = MAXLOC(ABS(erro),1)
            if (erro(npos) * errn(npos) < 0.0) domiR = -domiR
            fs0 = fs0 + domiR / (1._dp - domiR) * errn
            Write(ounit,*) '    ...Fission Source Extrapolated...'
       !     Write(*,*) '    ...Fission Source Extrapolated...'
        End if
        e1 = e2                         ! Save Integrated Fission Source Error

!############################    Update Keff    ###########################
        Call Integrate(fs,fs0,delv,nk)  ! Integrate Fission Source
        Keff = Keo * fs / fsd           ! Update Keff

!##################   Write Outer Iteration Evolution   ###################
        Call RelE(fs0,fs0d,fser,nk)     ! Search Maximum Point Wise Fission Source Relative Error
        Call RelEg(f0,f0d,fer,nk,ng)    ! Search Maximum Point Wise Flux Error
        if (.TRUE.) Then
            Write(ounit,'(3X,I5,2X,F13.6,2ES15.5)') o, Keff, fser, fer
          !  Write(*,'(3X,I5,2X,F13.6,2ES15.5)') o, Keff, fser, fer
        End if
        if ((fser < fserc) .AND. (fer < ferc)) Exit  ! if converge, exit.
    End Do

    if (o-1 == nout) Then
        Write(ounit,*)
        Write(ounit,*) '  Maximum Number of Outer Iteration is Reached.'
        Write(ounit,*) '  Iteration May Not Converge.'
        Write(ounit,*) '  Check Problem Specification or Change Iteration Control.'
       Write(ounit,*) '  Code is Stoping...'
        Write(*,*)
        Write(*,*) '  Maximum Number of Outer Iteration is Reached.'
        Write(*,*) '  Iteration May Not Converge.'
        Write(*,*) '  Check Problem Specification or Change Iteration Control.'
        Write(*,*) '  Code is Stoping...'
        Stop
    End if

    if (.TRUE.) Write(ounit,*)
    if (.TRUE.) Write(ounit,*) '***********************************************'
    if (.TRUE.) Write(ounit,   '(A40,F9.6)') 'Effective Factor Multiplication & 
    : Keff = ', Keff
    if (.TRUE.) Write(ounit,*) '***********************************************'
    if (.TRUE.) Write(*,*)
    if (.TRUE.) Write(*,*) '***********************************************'
    if (.TRUE.) Write(*, '(A40,F9.6)') 'Effective Factor Multiplication & 
    : Keff = ', Keff
    if (.TRUE.) Write(*,*) '***********************************************'

End Subroutine
        Subroutine Outerfs(Keff,fer,fser,f0,fx1,fy1,fz1,fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,ng,nk,&
						 order,R2,P2,R4,P4,nmat,mat,al,jo,ji,L0,D,sigr,chi,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,&
						 delx,dely,delz,delv,xyz,x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott)    ! To Perform Outer Iterations
Implicit None
Integer, Parameter :: dp = Selected_Real_Kind(10,15), Ounit = 100
Integer, Intent(in) :: ng, nk, order, nmat, nxx, nyy, nzz, x_east, x_west, y_north, y_south, z_top, z_bott
Integer, Dimension(nk), Intent(in) :: mat
Integer, Dimension(nxx), Intent(in) :: x_smax, x_smin
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin 
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk), Intent(in) :: delv, ix, iy, iz
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz
Real(dp), Dimension(nk,ng,ng), Intent(in) :: sigs
Real(dp), Dimension(nk,ng), Intent(in) :: nu_sigf, D, sigr
Real(dp), Dimension(nmat,ng), Intent(in) :: chi
Real(dp), Dimension(nk,ng,6), Intent(inout):: al, jo, ji
Real(dp), Dimension(nk,ng,3), Intent(inout) :: L0
Real(dp), Dimension(nk,ng,6,6),Intent(in) :: R2, P2, R4
Real(dp), Dimension(nk,ng,6,7),Intent(in) :: P4
Real(dp), Dimension(nk,ng), intent(inout) :: f0, fx1, fy1, fz1, fx2, fy2, fz2
Real(dp), Dimension(nk), intent(inout) :: fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2     ! Fission Source & Fission Source Moments
Real(dp), Intent(inout) :: Keff
! Local Variables
Integer :: nout = 500                ! Maximum Outer Iteration
Integer :: nac = 5                   ! Number of Outer Iteration Before Next Source Extrapolation
Real(dp) :: ferc = 1.e-5             ! Flux Error Criteria
Real(dp) :: fserc = 1.e-5            ! Fission Source Error Criteria
Real(dp), intent(out) :: fer, fser	! Flux and Fission Source Error in BCSEARCH calcs.
Real(dp) :: Keo                      ! Old Multiplication Factor (Keff)
Real(dp), Dimension(nk) :: fs0d      ! Old Fission Source
Real(dp), Dimension(nk,ng) :: f0d    ! Old Flux
Real(dp) :: fsd                      ! Old Integrated Fission Sources
Real(dp) :: fs                       ! New Integrated Fission Sources
Real(dp) :: domiR, e1, e2
Integer :: g, o, npos, n
Real(dp), Dimension(nk,ng,7) :: Q
Real(dp), Dimension(nk) :: errn, erro

!#####################    Initialize Fission Source    ####################
    Call FSrc(fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,&
              f0,fx1,fy1,fz1,fx2,fy2,fz2,nu_sigf,ng,nk)

errn = 1._dp          ! Variable Fission Source Error
    Call Integrate(fs,fs0,delv,nk)
!write (*,*) 'AFTER CALL Integrate1'
!	write(*,*)  'Keo=', Keo
!				do n=1, nk
!					if (n>50 .and. n<70) then
!						write(*,*)  'fs0=', fs0(n)
!					endif
!				enddo
!	write(*,*)  'fs=', fs
    Call Integrate(e1,errn,delv,nk)
!write (*,*) 'AFTER CALL IntegrateA'
!				do n=1, nk
!					if (n>50 .and. n<70) then
!						write(*,*)  'errn=', errn(n)
!					endif
!				enddo
!	write(*,*)  'e1=', e1
!#######################    Start Outer Iteration   #######################
    Do o=1, nout
        f0d = f0          ! Save Old Fluxes
        fsd = fs          ! Save Old Integrated Fission Source
        fs0d  = fs0       ! Save Old Fission Source
        Keo = Keff        ! Save Old Multiplication Fctor
        erro = errn       ! Save Old Fission Source Error/Difference
!	write(*,*) 'ITERATION=', o
        Do g = 1, ng
!#######################   Calculate Total Source   #######################
				do n=1, nk
					if (n>40 .and. n<43) then
					!	write(*,*) '***********************************'
					!	write(*,*) 'BEFORE CALL TSrc'
					!		write(*,*) 'group=', g, 'source=', Q(n,g,:)
					!		write(*,*) 'group=', g, 'f0=', f0(n,:)
					endif
				enddo
            Call TSrc(Q,g,Keo,ng,nmat,nk,chi,mat,sigs,f0,fx1,fy1,fz1,&
                      fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2)
!#######################    Start Inner Iteration   #######################
						!	write(*,*) 'BEFORE CALL INNER'
					!	write(*,*) '***********************************'
					!	write(*,*) 'BEFORE CALL INNER'
					!	write(*,*) '***********************************'
				do n=1, nk
				if (g==2) then
					if (n>1 .and. n<70) then
					!	write(*,*) 'group=', g, 'source=', Q(n,g,:)
						!	write(*,*) 'node=', n,'group=', g, 'node=', n,  'f0=', f0(n,g)
						!	write(*,*) 'node=', n,'group=', g,'node=', n,  'fx1=', fx1(n,g)
						!	write(*,*) 'node=', n,'group=', g,'node=', n,  'fy1=', fy1(n,g)
						!	write(*,*) 'node=', n,'group=', g,'node=', n,  'fz1=', fz1(n,g)
						!	write(*,*) 'node=', n,'group=', g,'node=', n,  'fx2=', fx2(n,g)
						!	write(*,*) 'node=', n,'group=', g,'node=', n,  'fy2=', fy2(n,g)
						!	write(*,*) 'node=', n,'group=', g,'node=', n,  'fz2=', fz2(n,g)
						!	write(*,*) 'node=', n,'group=', g,'node=', n,  'jo=', jo(n,g,:)
						!	write(*,*) 'node=', n,'group=', g,'node=', n,  'ji=', ji(n,g,:)
						!	write(*,*) 'node=', n,'group=', g,'node=', n,  'L0=', L0(n,g,:)
						!	write(*,*) 'node=', n,'group=', g,'node=', n,  'Q=', Q(n,g,:)
						!	write(*,*) 'node=', n,'group=', g,'node=', n,  'R4=', R4(n,g,6,6)
						!	write(*,*) 'node=', n,'group=', g,'node=', n,  'P4=', P4(n,g,6,7)
						!	write(*,*) 'node=', n,'group=', g,'node=', n,  'sigr=', sigr(n,g)
						!	write(*,*) 'node=', n,'group=', g,'node=', n,  'D=', D(n,g)
						!	write(*,*) 'node=', n,'group=', g,'node=', n,  'delx=', delx(3)
						!	write(*,*) 'node=', n,'group=', g,'node=', n,  'dely=', dely(3)
						!	write(*,*) 'node=', n,'group=', g,'node=', n,  'delz=', delz(3)
						!	write(*,*) 'node=', n,'group=', g,'node=', n,  'nxx=', nxx
						!	write(*,*) 'node=', n,'group=', g,'node=', n,  'nyy=', nyy
						!	write(*,*) 'node=', n,'group=', g,'node=', n,  'nzz=', nzz
						!	write(*,*) 'node=', n,'group=', g, 'fs0=', fs0(n)
						!	write(*,*) 'node=', n,'group=', g, 'fsx1=', fsx1(n)
						!	write(*,*) 'node=', n,'group=', g, 'fsy1=', fsy1(n)
						!	write(*,*) 'node=', n,'group=', g, 'fsz1=', fsz1(n)
						!	write(*,*) 'node=', n,'group=', g, 'fsx2=', fsx2(n)
						!	write(*,*) 'node=', n,'group=', g, 'fsy2=', fsy2(n)
						!	write(*,*) 'node=', n,'group=', g, 'fsz2=', fsz2(n)
						!	write(*,*) 'node=', n,'group=', g, 'Keff=', Keo
					endif
					endif
				enddo
            Call Inner(f0,fx1,fy1,fz1,fx2,fy2,fz2,g,ng,nk,order,jo,ji,L0,Q,R2,P2,R4,P4,al,D,sigr,nxx,nyy,nzz,ix,iy,iz,delx,&
                       dely,delz,x_smax,x_smin,y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott)
	!		if (p==250) then
				do n=1, nk
					if (n>40 .and. n<43) then
					!	write(*,*) '***********************************'
					!	write(*,*) '***********************************'
					!		write(*,*) 'AFTER CALL INNER'
					!		write(*,*) 'group=', g, 'node=', n,  'f0=', f0(n,g)
					!	write(*,*) '***********************************'
					endif
				enddo
	!		endif
				!		write(*,*) '***********************************'
					!	write(*,*) '***********************************'

 End Do

			!	do n=1, nk
			!		if (n>40 .and. n<45) then
			!			write(*,*) '***********************************'
			!				write(*,*) 'BEFORE CALL FSrc'
			!				write(*,*)  'f0=', f0(n,:)
			!			write(*,*) '***********************************'
			!		endif
			!	enddo

!########################   Update Fission Source   #######################
        Call FSrc(fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,&
                  f0,fx1,fy1,fz1,fx2,fy2,fz2,nu_sigf,ng,nk)
			!	do n=1, nk
			!		if (n>40 .and. n<45) then
			!			write(*,*) '***********************************'
			!				write(*,*) 'AFTER CALL FSrc'
			!				write(*,*)  'fs0=', fs0(n)
			!			write(*,*) '***********************************'
			!		endif
			!	enddo
        errn = fs0 - fs0d
        Call Integrate(e2,ABS(errn),delv,nk)    ! Variable Integrated Fission Source Error
!write (*,*) 'AFTER CALL IntegrateB'
!				do n=1, nk
!					if (n>50 .and. n<70) then
!						write(*,*)  'errn=', errn(n)
!					endif
!				enddo
!	write(*,*)  'e2=', e2

!####################   Fission Source Extrapolation   ####################
        if (MOD(o,nac) == 0) Then
            domiR = e2 / e1
            npos = MAXLOC(ABS(erro),1)
            if (erro(npos) * errn(npos) < 0.0) domiR = -domiR
            fs0 = fs0 + domiR / (1._dp - domiR) * errn
      !      Write(ounit,*) '    ...Fission Source Extrapolated...'
       !     Write(*,*) '    ...Fission Source Extrapolated...'
        End if
        e1 = e2                         ! Save Integrated Fission Source Error

!############################    Update Keff    ###########################
        Call Integrate(fs,fs0,delv,nk)  ! Integrate Fission Source

!write (*,*) '###########'
!write (*,*) 'AFTER CALL Integrate2'
!	write(*,*)  'Keo=', Keo
!				do n=1, nk
!					if (n>50 .and. n<70) then
!						write(*,*)  'fs0=', fs0(n)
!					endif
!				enddo
!	write(*,*)  'fs=', fs
!	write(*,*)  'fsd=', fsd
        Keff = Keo * fs / fsd           ! Update Keff

!##################   Write Outer Iteration Evolution   ###################
        Call RelE(fs0,fs0d,fser,nk)     ! Search Maximum Point Wise Fission Source Relative Error
!write (*,*) 'AFTER CALL RelE'
!				do n=1, nk
!					if (n==50) then
!						write(*,*)  'fs0=', fs0(n)
!						write(*,*)  'fs0d=', fs0d(n)
!					endif
!				enddo
!	write(*,*)  'fser=', fser
	
        Call RelEg(f0,f0d,fer,nk,ng)    ! Search Maximum Point Wise Flux Error
!write (*,*) 'AFTER CALL RelEg'
!				do n=1, nk
!					if (n==50) then
!						write(*,*)  'f0=', f0(n,:)
!						write(*,*)  'f0d=', f0d(n,:)
!					endif
!				enddo
!	write(*,*)  'fer=', fer

        if (.TRUE.) Then
            Write(ounit,'(3X,I5,2X,F13.6,2ES15.5)') o, Keff, fser, fer
         !   Write(*,'(3X,I5,2X,F13.6,2ES15.5)') o, Keff, fser, fer
        End if
        if ((fser < fserc) .AND. (fer < ferc)) Exit  ! if converge, exit.
    End Do

    if (o-1 == nout) Then
        Write(ounit,*)
        Write(ounit,*) '  Maximum Number of Outer Iteration is Reached.'
        Write(ounit,*) '  Iteration May Not Converge.'
        Write(ounit,*) '  Check Problem Specification or Change Iteration Control.'
       Write(ounit,*) '  Code is Stoping...'
  !      Write(*,*)
   !     Write(*,*) '  Maximum Number of Outer Iteration is Reached.'
    !    Write(*,*) '  Iteration May Not Converge.'
     !   Write(*,*) '  Check Problem Specification or Change Iteration Control.'
      !  Write(*,*) '  Code is Stoping...'
        Stop
    End if

    if (.TRUE.) Write(ounit,*)
    if (.TRUE.) Write(ounit,*) '***********************************************'
    if (.TRUE.) Write(ounit,   '(A40,F9.6)') 'Effective Factor Multiplication & 
    : Keff = ', Keff
    if (.TRUE.) Write(ounit,*) '***********************************************'
    if (.TRUE.) Write(*,*)
    !if (.TRUE.) Write(*,*) '***********************************************'
    !if (.TRUE.) Write(*,*) 'Forward Calculation'
    !if (.TRUE.) Write(*,*) '***********************************************'
    !if (.TRUE.) Write(*, '(A40,F9.6)') 'Effective Factor Multiplication & 
    !: Keff = ', Keff
    !if (.TRUE.) Write(*,*) '***********************************************'

End Subroutine

        Subroutine OuterAdj(Keff,f0,fx1,fy1,fz1,fx2,fy2,fz2,ng,nk,order,R2,P2,R4,P4,nmat,chi,mat,al,jo,ji,L0,&
                            D,sigr,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,&
                            x_smax,x_smin,y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott)    !    To Perform Adjoint Outer Iteration
Implicit None
Integer, Parameter :: dp = Selected_Real_Kind(10,15), Ounit = 100
Integer, Intent(in) :: ng, nk, order, nmat, nxx, nyy, nzz, x_east, x_west, y_north, y_south, z_top, z_bott
Integer, Dimension(nk), Intent(in) :: mat
Integer, Dimension(nxx), Intent(in) :: x_smax, x_smin
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin 
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk), Intent(in) :: delv, ix, iy, iz
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz
Real(dp), Dimension(nk,ng,ng), Intent(in) :: sigs
Real(dp), Dimension(nk,ng), Intent(in) :: nu_sigf, D, sigr
Real(dp), Dimension(nmat,ng), Intent(in) :: chi
Real(dp), Dimension(nk,ng,6), Intent(inout) :: al, jo, ji
Real(dp), Dimension(nk,ng,3), Intent(inout) :: L0
Real(dp), Dimension(nk,ng,6,6),Intent(in) :: R2, P2, R4
Real(dp), Dimension(nk,ng,6,7),Intent(in) :: P4
Real(dp), Dimension(nk,ng), intent(inout) :: f0, fx1, fy1, fz1, fx2, fy2, fz2
Real(dp), Intent(inout) :: Keff
! Local Variables
Integer :: nout = 500               ! Maximum Outer Iteration
Integer :: nac = 5                  ! Number of Outer Iteration Before Next Source EXTRAPOLATION
Real(dp) :: ferc = 1.e-5            ! Flux Error Criteria
Real(dp) :: fserc = 1.e-5           ! Fission Source Error Criteria
Real(dp) :: fer, fser               ! Flux and Fission Source Error in BCSEARCH calcs.
Real(dp), Dimension(nk) :: fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2     ! Fission Source and Fission Source Moments
Real(dp) :: Keo                     ! Old Multiplication Factor (Keff)
Real(dp), Dimension(nk) :: fs0d     ! Old Fission Source
Real(dp), Dimension(nk,ng) :: f0d   ! Old Flux
Real(dp) :: fsd                     ! Old Integrated Fission Sources
Real(dp) :: fs                      ! New Integrated Fission Sources
Real(dp) :: domiR, e1, e2
Integer :: g, o, npos
Real(dp), Dimension(nk,ng,7) :: Q
Real(dp), Dimension(nk) :: errn, erro

Open (Unit=ounit, File='app/Output/NEM.out')

!#####################    Initialize Fission Source    ####################
    Call FSrcAdj(fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,mat,&
                 chi,f0,fx1,fy1,fz1,fx2,fy2,fz2,ng,nk,nmat)
errn = 1._dp          ! Variable Fission Source Error
    Call Integrate(fs,fs0,delv,nk)
    Call Integrate(e1,errn,delv,nk)

!#######################    Start Outer Iteration   #######################
    Do o=1, nout
        f0d = f0          ! Save Old Fluxes
        fsd = fs          ! Save Old Integrated Fission Source
        fs0d  = fs0       ! Save Old Fission Source
        Keo = Keff        ! Save Old Multiplication Fctor
        erro = errn       ! Save Old Fission Source Error/Difference

        Do g = ng, 1, -1
!#######################   Calculate Total Source   #######################
            Call TSrcAdj(Q,Keo,sigs,nu_sigf,f0,fx1,fy1,fz1,fx2,fy2,fz2,&
                         fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,g,ng,nk)
!#######################    Start Inner Iteration   #######################
            Call Inner(f0,fx1,fy1,fz1,fx2,fy2,fz2,g,ng,nk,order,jo,ji,L0,Q,R2,P2,R4,P4,al,D,sigr,nxx,nyy,nzz,ix,iy,iz,delx,&
                       dely,delz,x_smax,x_smin,y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott)
		End Do

!########################   Update Fission Source   #######################
        Call FSrcAdj(fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,mat,&
                     chi,f0,fx1,fy1,fz1,fx2,fy2,fz2,ng,nk,nmat)
        errn = fs0 - fs0d
        Call Integrate(e2,ABS(errn),delv,nk)    ! Variable Integrated Fission Source Error

!####################   Fission Source Extrapolation   ####################
        if (MOD(o,nac) == 0) Then
            domiR = e2 / e1
            npos = MAXLOC(ABS(erro),1)
            if (erro(npos) * errn(npos) < 0.0) domiR = -domiR
            fs0 = fs0 + domiR / (1._dp - domiR) * errn
            Write(ounit,*) '    ...Fission Source Extrapolated...'
        End if
        e1 = e2                         ! Save Integrated Fission Source Error

!############################    Update Keff    ###########################
        Call Integrate(fs,fs0,delv,nk)  ! Integrate Fission Source
        Keff = Keo * fs / fsd

!##################   Write Outer Iteration Evolution   ###################
        Call RelE(fs0,fs0d,fser,nk)     ! Search Maximum Point Wise Fission Source Relative Error
        Call RelEg(f0,f0d,fer,nk,ng)    ! Search Maximum Point Wise Flux Error
        if (.TRUE.) Write(ounit,'(I5,F13.6,2ES15.5)') o, Keff, fser, fer
       ! if (.TRUE.) Write(*,'(I5,F13.6,2ES15.5)') o, Keff, fser, fer
        if ((fser < fserc) .AND. (fer < ferc)) Exit  ! if converge, exit.
    End Do

    if (o-1 == nout) Then
        Write(ounit,*)
        Write(ounit,*) '  Maximum Number of Outer Iteration is Reached.'
        Write(ounit,*) '  Iteration May Not Converge.'
        Write(ounit,*) '  Check Problem Specification or Change Iteration Control.'
        Write(ounit,*) '  Code is Stoping...'
        Write(*,*)
        Write(*,*) '  Maximum Number of Outer Iteration is Reached.'
        Write(*,*) '  Iteration May Not Converge.'
        Write(*,*) '  Check Problem Specification or Change Iteration Control.'
        Write(*,*) '  Code is Stoping...'
        Stop
    End if

    if (.TRUE.) Write(ounit,*)
    if (.TRUE.) Write(ounit,*) '***********************************************'
    if (.TRUE.) Write(ounit, '(A40,F9.6)') 'Effective Factor Multiplication & 
    : Keff = ', Keff
    if (.TRUE.) Write(ounit,*) '***********************************************'
    if (.TRUE.) Write(*,*)
    if (.TRUE.) Write(*,*) '***********************************************'
    if (.TRUE.) Write(*, '(A40,F9.6)') 'Effective Factor Multiplication & 
    : Keff = ', Keff
    if (.TRUE.) Write(*,*) '***********************************************'

End Subroutine
        Subroutine OuterAdjfs(Keff,f0,fx1,fy1,fz1,fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,ng,nk,order,&
							R2,P2,R4,P4,nmat,chi,mat,al,jo,ji,L0,D,sigr,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,&
							delx,dely,delz,delv,x_smax,x_smin,y_smax,y_smin,xyz,x_east,x_west,y_north,&
							y_south,z_top,z_bott)    !    To Perform Adjoint Outer Iteration
Implicit None
Integer, Parameter :: dp = Selected_Real_Kind(10,15), Ounit = 100
Integer, Intent(in) :: ng, nk, order, nmat, nxx, nyy, nzz, x_east, x_west, y_north, y_south, z_top, z_bott
Integer, Dimension(nk), Intent(in) :: mat
Integer, Dimension(nxx), Intent(in) :: x_smax, x_smin
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin 
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk), Intent(in) :: delv, ix, iy, iz
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz
Real(dp), Dimension(nk,ng,ng), Intent(in) :: sigs
Real(dp), Dimension(nk,ng), Intent(in) :: nu_sigf, D, sigr
Real(dp), Dimension(nmat,ng), Intent(in) :: chi
Real(dp), Dimension(nk,ng,6), Intent(inout) :: al, jo, ji
Real(dp), Dimension(nk,ng,3), Intent(inout) :: L0
Real(dp), Dimension(nk,ng,6,6),Intent(in) :: R2, P2, R4
Real(dp), Dimension(nk,ng,6,7),Intent(in) :: P4
Real(dp), Dimension(nk,ng), intent(inout) :: f0, fx1, fy1, fz1, fx2, fy2, fz2
Real(dp), Dimension(nk), intent(inout) :: fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2     ! Fission Source and Fission Source Moments
Real(dp), Intent(inout) :: Keff
! Local Variables
Integer :: nout = 500               ! Maximum Outer Iteration
Integer :: nac = 5                  ! Number of Outer Iteration Before Next Source EXTRAPOLATION
Real(dp) :: ferc = 1.e-5            ! Flux Error Criteria
Real(dp) :: fserc = 1.e-5           ! Fission Source Error Criteria
Real(dp) :: fer, fser               ! Flux and Fission Source Error in BCSEARCH calcs.
Real(dp) :: Keo                     ! Old Multiplication Factor (Keff)
Real(dp), Dimension(nk) :: fs0d     ! Old Fission Source
Real(dp), Dimension(nk,ng) :: f0d   ! Old Flux
Real(dp) :: fsd                     ! Old Integrated Fission Sources
Real(dp) :: fs                      ! New Integrated Fission Sources
Real(dp) :: domiR, e1, e2
Integer :: g, o, npos
Real(dp), Dimension(nk,ng,7) :: Q
Real(dp), Dimension(nk) :: errn, erro

Open (Unit=ounit, File='app/Output/NEM.out')

!#####################    Initialize Fission Source    ####################
    Call FSrcAdj(fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,mat,&
                 chi,f0,fx1,fy1,fz1,fx2,fy2,fz2,ng,nk,nmat)
errn = 1._dp          ! Variable Fission Source Error
    Call Integrate(fs,fs0,delv,nk)
    Call Integrate(e1,errn,delv,nk)

!#######################    Start Outer Iteration   #######################
    Do o=1, nout
        f0d = f0          ! Save Old Fluxes
        fsd = fs          ! Save Old Integrated Fission Source
        fs0d  = fs0       ! Save Old Fission Source
        Keo = Keff        ! Save Old Multiplication Fctor
        erro = errn       ! Save Old Fission Source Error/Difference

        Do g = ng, 1, -1
!#######################   Calculate Total Source   #######################
            Call TSrcAdj(Q,Keo,sigs,nu_sigf,f0,fx1,fy1,fz1,fx2,fy2,fz2,&
                         fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,g,ng,nk)
!#######################    Start Inner Iteration   #######################
            Call Inner(f0,fx1,fy1,fz1,fx2,fy2,fz2,g,ng,nk,order,jo,ji,L0,Q,R2,P2,R4,P4,al,D,sigr,nxx,nyy,nzz,ix,iy,iz,delx,&
                       dely,delz,x_smax,x_smin,y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott)
        End Do

!########################   Update Fission Source   #######################
        Call FSrcAdj(fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,mat,&
                     chi,f0,fx1,fy1,fz1,fx2,fy2,fz2,ng,nk,nmat)
        errn = fs0 - fs0d
        Call Integrate(e2,ABS(errn),delv,nk)    ! Variable Integrated Fission Source Error

!####################   Fission Source Extrapolation   ####################
        if (MOD(o,nac) == 0) Then
            domiR = e2 / e1
            npos = MAXLOC(ABS(erro),1)
            if (erro(npos) * errn(npos) < 0.0) domiR = -domiR
            fs0 = fs0 + domiR / (1._dp - domiR) * errn
            Write(ounit,*) '    ...Fission Source Extrapolated...'
        End if
        e1 = e2                         ! Save Integrated Fission Source Error

!############################    Update Keff    ###########################
        Call Integrate(fs,fs0,delv,nk)  ! Integrate Fission Source
        Keff = Keo * fs / fsd

!##################   Write Outer Iteration Evolution   ###################
        Call RelE(fs0,fs0d,fser,nk)     ! Search Maximum Point Wise Fission Source Relative Error
        Call RelEg(f0,f0d,fer,nk,ng)    ! Search Maximum Point Wise Flux Error
        if (.TRUE.) Write(ounit,'(I5,F13.6,2ES15.5)') o, Keff, fser, fer
     !   if (.TRUE.) Write(*,'(I5,F13.6,2ES15.5)') o, Keff, fser, fer
        if ((fser < fserc) .AND. (fer < ferc)) Exit  ! if converge, exit.
    End Do

    if (o-1 == nout) Then
        Write(ounit,*)
        Write(ounit,*) '  Maximum Number of Outer Iteration is Reached.'
        Write(ounit,*) '  Iteration May Not Converge.'
        Write(ounit,*) '  Check Problem Specification or Change Iteration Control.'
        Write(ounit,*) '  Code is Stoping...'
        Write(*,*)
        Write(*,*) '  Maximum Number of Outer Iteration is Reached.'
        Write(*,*) '  Iteration May Not Converge.'
        Write(*,*) '  Check Problem Specification or Change Iteration Control.'
        Write(*,*) '  Code is Stoping...'
        Stop
    End if

    if (.TRUE.) Write(ounit,*)
    if (.TRUE.) Write(ounit,*) '***********************************************'
    if (.TRUE.) Write(ounit, '(A40,F9.6)') 'Effective Factor Multiplication & 
    : Keff = ', Keff
    if (.TRUE.) Write(ounit,*) '***********************************************'
    !if (.TRUE.) Write(*,*)
    !if (.TRUE.) Write(*,*) '***********************************************'
    !if (.TRUE.) Write(*,*) 'Adjoint Calculation'
    !if (.TRUE.) Write(*,*) '***********************************************'
    !if (.TRUE.) Write(*, '(A40,F9.6)') 'Effective Factor Multiplication & 
    !: Keff = ', Keff
    !if (.TRUE.) Write(*,*) '***********************************************'

End Subroutine

        Subroutine OuterFx(f0,fx1,fy1,fz1,fx2,fy2,fz2,ng,nk,order,R2,P2,R4,P4,nmat,chi,mat,al,jo,ji,L0,&
                           D,sigr,sigs,siga,nu_sigf,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,eSrc,&
                           x_smax,x_smin,y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott)    !To perform fixed-source outer iteration
Implicit None
Integer, Parameter :: dp = Selected_Real_Kind(10,15), Ounit = 100
Integer, Intent(in) :: ng, nk, order, nmat, nxx, nyy, nzz, x_east, x_west, y_north, y_south, z_top, z_bott
Integer, Dimension(nk), Intent(in) :: mat
Integer, Dimension(nxx), Intent(in) :: x_smax, x_smin
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin 
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk), Intent(in) :: delv, ix, iy, iz
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz
Real(dp), Dimension(nk,ng,ng), Intent(in) :: sigs
Real(dp), Dimension(nk,ng), Intent(in) :: siga, nu_sigf, D, sigr, eSrc
Real(dp), Dimension(nmat,ng), Intent(in) :: chi
Real(dp), Dimension(nk,ng,6), Intent(inout):: al, jo, ji
Real(dp), Dimension(nk,ng,3), Intent(inout) :: L0
Real(dp), Dimension(nk,ng,6,6),Intent(in) :: R2, P2, R4
Real(dp), Dimension(nk,ng,6,7),Intent(in) :: P4
Real(dp), Dimension(nk,ng), intent(inout) :: f0, fx1, fy1, fz1, fx2, fy2, fz2
Real(dp) :: Keff
! Local Variables
Integer :: nout = 500               ! Maximum Outer Iteration
Integer :: nac = 5                  ! Number of Outer Iteration Before Next Source Extrapolation
Real(dp) :: ferc = 1.e-5            ! Flux Error Criteria
Real(dp) :: fserc = 1.e-5           ! Fission Source Error Criteria
Real(dp) :: fer, fser               ! Flux and Fission Source Error in BCSEARCH calcs.
Real(dp), Dimension(nk) :: fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2     ! Fission Source & Fission Source Moments
Real(dp), Dimension(nk) :: fs0d     ! Old Fission Source
Real(dp), Dimension(nk,ng) :: f0d   ! Old Flux
Real(dp) :: domiR, e1, e2
Integer :: g, o, npos
Real(dp), Dimension(nk,ng,7) :: Q
Real(dp), Dimension(nk) :: errn, erro

Open (Unit=ounit, File='app/Output/NEM.out')

!#####################    Initialize Fission Source    ####################
    Call FSrc(fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,&
              f0,fx1,fy1,fz1,fx2,fy2,fz2,nu_sigf,ng,nk)

errn = 1._dp          ! Variable Fission Source Error
    Call Integrate(e1,errn,delv,nk)

!#######################    Start Outer Iteration   #######################
    Do o=1, nout
        f0d = f0          ! Save Old Fluxes
        fs0d  = fs0       ! Save Old Fission Source
        erro = errn       ! Save Old Fission Source Error/Difference

        Do g = 1, ng
!#######################   Calculate Total Source   #######################
            Call TSrcFx(Q,g,ng,nmat,nk,chi,mat,sigs,f0,fx1,fy1,fz1,&
                        fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,eSrc)
!#######################    Start Inner Iteration   #######################
            Call Inner(f0,fx1,fy1,fz1,fx2,fy2,fz2,g,ng,nk,order,jo,ji,L0,Q,R2,P2,R4,P4,al,D,sigr,nxx,nyy,nzz,ix,iy,iz,delx,&
                       dely,delz,x_smax,x_smin,y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott)
		End Do

!########################   Update Fission Source   #######################
        Call FSrc(fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,&
                  f0,fx1,fy1,fz1,fx2,fy2,fz2,nu_sigf,ng,nk)

        errn = fs0 - fs0d
        Call Integrate(e2,ABS(errn),delv,nk)    ! Variable Integrated Fission Source Error

!####################   Fission Source Extrapolation   ####################
        if (Mod(o,nac) == 0) Then
            domiR = e2 / e1
            npos = MAXLOC(ABS(erro),1)
            if (erro(npos) * errn(npos) < 0.0) domiR = -domiR
            fs0 = fs0 + domiR / (1._dp - domiR) * errn
            Write(ounit,*) '    ...Fission Source Extrapolated...'
        End if
        e1 = e2                         ! Save Integrated Fission Source Error

!##################   Write Outer Iteration Evolution   ###################
        Call RelE(fs0,fs0d,fser,nk)     ! Search Maximum Point Wise Fission Source Relative Error
        Call RelEg(f0,f0d,fer,nk,ng)    ! Search Maximum Point Wise Flux Error
        Write(ounit,'(I5,2ES15.5)') o, fser, fer  ! Write Outer Iteration Evolution
      !  Write(*,'(I5,2ES15.5)') o, fser, fer  ! Write Outer Iteration Evolution
        if ((fser < fserc) .AND. (fer < ferc)) Exit  ! if converge, exit.
    End Do
    if (o-1 == nout) Then
        Write(ounit,*)
        Write(ounit,*) '  Maximum Number of Outer Iteration is Reached.'
        Write(ounit,*) '  Iteration May Not Converge.'
        Write(ounit,*) '  Check Problem Specification or Change Iteration Control.'
        Write(ounit,*) '  Code is Stoping...'
        Stop
    End if

!#########################     Calculate Keff     #########################
    Call MultF(Keff,ng,nk,siga,nu_sigf,f0,jo,ji,nxx,nyy,nzz,ix,iy,iz,&
               delx,dely,delz,x_smax,x_smin,y_smax,y_smin)

    if (o-1 == nout) Then
        Write(ounit,*)
        Write(ounit,*) '  Maximum Number of Outer Iteration is Reached.'
        Write(ounit,*) '  Iteration May Not Converge.'
        Write(ounit,*) '  Check Problem Specification or Change Iteration Control.'
        Write(ounit,*) '  Code is Stoping...'
        Write(*,*)
        Write(*,*) '  Maximum Number of Outer Iteration is Reached.'
        Write(*,*) '  Iteration May Not Converge.'
        Write(*,*) '  Check Problem Specification or Change Iteration Control.'
        Write(*,*) '  Code is Stoping...'
        Stop
    End if

    if (.TRUE.) Write(ounit,*)
    if (.TRUE.) Write(ounit,*) '***********************************************'
    if (.TRUE.) Write(ounit, '(A40,F9.6)') 'Effective Factor Multiplication & 
    : Keff = ', Keff
    if (.TRUE.) Write(ounit,*) '***********************************************'
    !if (.TRUE.) Write(*,*)
    !if (.TRUE.) Write(*,*) '***********************************************'
    !if (.TRUE.) Write(*,*) 'Fixed Source Calculation'
    !if (.TRUE.) Write(*,*) '***********************************************'
    !if (.TRUE.) Write(*, '(A40,F9.6)') 'Effective Factor Multiplication & 
    !: Keff = ', Keff
    !if (.TRUE.) Write(*,*) '***********************************************'

End Subroutine

		SUBROUTINE outer4th(Keff,fer,fser,f0,fx1,fy1,fz1,fx2,fy2,fz2,maxn,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,&
					ng,nk,order,R2,P2,R4,P4,nmat,mat,al,jo,ji,L0,&
                    D,sigr,chi,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,delv,xyz,&
                    x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott)!    To perform normal outer iteration when %THER card present
Implicit None
Integer, Parameter :: dp = Selected_Real_Kind(10,15), Ounit = 100
Integer, Intent(in) :: ng, nk, order, nmat, nxx, nyy, nzz, x_east, x_west, y_north, y_south, z_top, z_bott, maxn
Integer, Dimension(nk), Intent(in) :: mat
Integer, Dimension(nxx), Intent(in) :: x_smax, x_smin
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin 
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk), Intent(in) :: delv, ix, iy, iz
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz
Real(dp), Dimension(nk,ng,ng), Intent(in) :: sigs
Real(dp), Dimension(nk,ng), Intent(in) :: nu_sigf, D, sigr
Real(dp), Dimension(nmat,ng), Intent(in) :: chi
Real(dp), Dimension(nk,ng,6), Intent(inout):: al, jo, ji
Real(dp), Dimension(nk,ng,3), Intent(inout) :: L0
Real(dp), Dimension(nk,ng,6,6),Intent(in) :: R2, P2, R4
Real(dp), Dimension(nk,ng,6,7),Intent(in) :: P4
Real(dp), Dimension(nk,ng), intent(inout) :: f0, fx1, fy1, fz1, fx2, fy2, fz2
Real(dp), Dimension(nk), intent(inout) :: fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2     ! Fission Source and Fission Source Moments
Real(dp), Intent(inout) :: Keff
! Local Variables
Integer :: nac = 5                   ! Number of Outer Iteration Before Next Source Extrapolation
Real(dp) :: ferc = 1.e-5             ! Flux Error Criteria
Real(dp) :: fserc = 1.e-5            ! Fission Source Error Criteria
Real(dp), intent(out) :: fer, fser                ! Flux and Fission Source Error in BCSEARCH calcs.
Real(dp) :: Keo                      ! Old Multiplication Factor (Keff)
Real(dp), Dimension(nk) :: fs0d      ! Old Fission Source
Real(dp), Dimension(nk,ng) :: f0d    ! Old Flux
Real(dp) :: fsd                      ! Old Integrated Fission Sources
Real(dp) :: fs                       ! New Integrated Fission Sources
Real(dp) :: domiR, e1, e2
Integer :: g, o, npos
Real(dp), Dimension(nk,ng,7) :: Q
Real(dp), Dimension(nk) :: errn, erro
!-------------------

!#####################    Initialize Fission Source    ####################
    Call FSrc(fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,&
              f0,fx1,fy1,fz1,fx2,fy2,fz2,nu_sigf,ng,nk)

errn = 1._dp          ! Variable Fission Source Error
    Call Integrate(fs,fs0,delv,nk)
    Call Integrate(e1,errn,delv,nk)

!#######################    Start Outer Iteration   #######################
DO o=1, maxn
    f0d = f0          ! Save old fluxes
    fsd = fs            ! Save old integrated fission source
    fs0d  = fs0       ! Save old fission source
    Keo = Keff          ! Save old multiplication factor
    erro = errn       ! Save old fission source error/difference
    DO g = 1, ng
!#######################   Calculate Total Source   #######################
            Call TSrc(Q,g,Keo,ng,nmat,nk,chi,mat,sigs,f0,fx1,fy1,fz1,&
                      fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2)
!#######################    Start Inner Iteration   #######################
            Call Inner(f0,fx1,fy1,fz1,fx2,fy2,fz2,g,ng,nk,order,jo,ji,L0,Q,R2,P2,R4,P4,al,D,sigr,nxx,nyy,nzz,ix,iy,iz,delx,&
                       dely,delz,x_smax,x_smin,y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott)    
	END DO
!########################   Update Fission Source   #######################
        Call FSrc(fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,&
                  f0,fx1,fy1,fz1,fx2,fy2,fz2,nu_sigf,ng,nk)
    errn = fs0 - fs0d
        Call Integrate(e2,ABS(errn),delv,nk)    ! Variable Integrated Fission Source Error

!####################   Fission Source Extrapolation   ####################
    IF (MOD(o,nac) == 0) THEN   ! Fission source extrapolation
        domiR = e2 / e1
        npos = MAXLOC(ABS(erro),1)
        IF (erro(npos) * errn(npos) < 0.0) domiR = -domiR
        fs0 = fs0 + domiR / (1._DP - domiR) * errn
    END IF
    e1 = e2                       ! Save integrated fission source error
!############################    Update Keff    ###########################
        Call Integrate(fs,fs0,delv,nk)  ! Integrate Fission Source
        Keff = Keo * fs / fsd           ! Update Keff

!##################   Write Outer Iteration Evolution   ###################
        Call RelE(fs0,fs0d,fser,nk)     ! Search Maximum Point Wise Fission Source Relative Error
        Call RelEg(f0,f0d,fer,nk,ng)    ! Search Maximum Point Wise Flux Error
        if (.TRUE.) Write(ounit,'(I5,F13.6,2ES15.5)') o, Keff, fser, fer
      !  if (.TRUE.) Write(*,'(I5,F13.6,2ES15.5)') o, Keff, fser, fer
    IF ((fser < fserc) .AND. (fer < ferc)) EXIT  ! If converge, exit.
END DO
 !Write(*,*) '***********************************************'
 !   if (.TRUE.) Write(*,*) 'Thermal Calculation'
 !Write(*,*) '***********************************************'
 !WRITE(*, '(A36,F9.6)') 'MULTIPLICATION EFFECTIVE (K-EFF) = ', Keff
 !Write(*,*) '***********************************************'

END SUBROUTINE

		SUBROUTINE outertf (Keff,fer,fser,f0,fx1,fy1,fz1,fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,maxi,ng,nk,order,&
							R2,P2,R4,P4,nmat,mat,al,jo,ji,L0,D,sigr,chi,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,&
							delv,xyz,x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott,iBeta,lamb,&
							tbeta,c0,cx1,cy1,cz1,cx2,cy2,cz2,ft,ftx1,fty1,ftz1,ftx2,fty2,ftz2,omeg,velo,ht,&
							errn, erro)
							!    To perform outer iteration for transient with flux transformation
	
Implicit None
Integer, Parameter :: dp = Selected_Real_Kind(10,15), Ounit = 100	
Integer, Intent(in) :: ng, nk, order, nmat, nxx, nyy, nzz, x_east, x_west, y_north, y_south, z_top, z_bott
INTEGER, PARAMETER :: nf = 6                       ! Number of delaye dneutron precusor family
Integer, Dimension(nk), Intent(in) :: mat
Integer, Dimension(nxx), Intent(in) :: x_smax, x_smin
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin 
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk), Intent(in) :: delv, ix, iy, iz
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz
Real(dp), Dimension(nk,ng,ng), Intent(in) :: sigs
Real(dp), Dimension(nk,ng), Intent(in) :: nu_sigf, D, sigr
Real(dp), Dimension(nmat,ng), Intent(in) :: chi
Real(dp), Dimension(nk,ng,6), Intent(inout):: al, jo, ji
Real(dp), Dimension(nk,ng,3), Intent(inout) :: L0
Real(dp), Dimension(nk,ng,6,6),Intent(in) :: R2, P2, R4
Real(dp), Dimension(nk,ng,6,7),Intent(in) :: P4
Real(dp), Dimension(nk,ng), intent(inout) :: f0, fx1, fy1, fz1, fx2, fy2, fz2
Real(dp), Dimension(nk), intent(inout) :: fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2     ! Fission Source & Fission Source Moments
Real(dp), Intent(inout) :: Keff
REAL(DP), DIMENSION(nf), intent(in) :: iBeta, lamb	! beta (delayed neutron fraction) and precusor decay constant
REAL(DP), intent(in) :: tbeta	! total beta
REAL(DP), DIMENSION(nk,nf), intent(in) :: c0, cx1, cy1, cz1, cx2, cy2, cz2  ! neutron precusor density
REAL(DP), DIMENSION(nk,ng), intent(in) :: ft, ftx1, fty1, ftz1, ftx2, fty2, ftz2  ! Parameters at previous time step
REAL(DP), DIMENSION(nk,ng), intent(in) :: omeg	! Exponential transformation constant
REAL(DP), DIMENSION(ng), intent(in) :: velo            ! Neutron velocity
REAL(DP), INTENT(IN) :: ht
! Local Variables
Integer :: nout = 500                ! Maximum Outer Iteration
Integer :: nac = 5                   ! Number of Outer Iteration Before Next Source Extrapolation
Real(dp) :: ferc = 1.e-5             ! Flux Error Criteria
Real(dp) :: fserc = 1.e-5            ! Fission Source Error Criteria
Real(dp), intent(out) :: fer, fser	! Flux and Fission Source Error in BCSEARCH calcs.
Real(dp) :: Keo                      ! Old Multiplication Factor (Keff)
Real(dp), Dimension(nk) :: fs0d      ! Old Fission Source
Real(dp), Dimension(nk,ng) :: f0d    ! Old Flux
Real(dp) :: fsd                      ! Old Integrated Fission Sources
Real(dp) :: fs                       ! New Integrated Fission Sources
Real(dp) :: domiR, e1, e2
Integer :: g, o, npos, p, n, y
Real(dp), Dimension(nk,ng,7) :: Q
Real(dp), Dimension(nk), intent(inout) :: errn, erro

LOGICAL, INTENT(OUT) :: maxi

    !Call FSrc(fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,&
     !         f0,fx1,fy1,fz1,fx2,fy2,fz2,nu_sigf,ng,nk)

errn = 1._DP
Call Integrate(e1,errn,delv,nk)    

!Start outer iteration
DO p=1, nout
    f0d = f0          ! Save old fluxes
    fs0d  = fs0       ! Save old fission source
    erro = errn       ! Save old fission source error/difference 

!	write(*,*) 'ITERATION=', o
        Do g = 1, ng
!#######################   Calculate Total Source   #######################
				do n=1, nk
					if (n==50) then
						!write(*,*) '***********************************'
						!write(*,*) 'BEFORE CALL TSrcT'
						!	write(*,*) 'group=', g, 'source=', Q(n,g,:)
						!	write(*,*) 'group=', g, 'f0=', f0(n,:)
					endif
				enddo
        CALL TSrcT(Q,g,ng,nmat,nk,chi,mat,sigs,f0,fx1,fy1,fz1,&
				   fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,&
				   iBeta,lamb,tbeta,c0,cx1,cy1,cz1,cx2,cy2,cz2,&
				   ft,ftx1,fty1,ftz1,ftx2,fty2,ftz2,omeg,velo,ht)
!#######################    Start Inner Iteration   #######################
				do n=1, nk
					if (n==50) then
						!write(*,*) '***********************************'
						!write(*,*) 'AFTER CALL TSrcT'
						!	write(*,*) 'group=', g, 'source=', Q(n,g,:)
						!write(*,*) '***********************************'
						!	write(*,*) 'BEFORE CALL INNER'
						!	write(*,*) 'group=', g, 'f0=', f0(n,:)
					endif
				enddo
 
        !!!Inner Iteration
        Call Inner(f0,fx1,fy1,fz1,fx2,fy2,fz2,g,ng,nk,order,jo,ji,L0,Q,R2,P2,R4,P4,al,D,sigr,nxx,nyy,nzz,ix,iy,iz,delx,&
				   dely,delz,x_smax,x_smin,y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott)
				do n=1, nk
					if (n==50) then
					!	write(*,*) '***********************************'
					!		write(*,*) 'AFTER CALL INNER'
					!		write(*,*) 'group=', g, 'f0=', f0(n,:)
					!	write(*,*) '***********************************'
					endif
				enddo
					!	write(*,*) '***********************************'
					!	write(*,*) '***********************************'
    END DO
!########################   Update Fission Source   #######################
    Call FSrc(fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,&
              f0,fx1,fy1,fz1,fx2,fy2,fz2,nu_sigf,ng,nk)  !Update fission source
				do n=1, nk
					if (n==50) then
					!	write(*,*) '***********************************'
					!		write(*,*) 'AFTER CALL FSrc'
					!		write(*,*)  'fs0=', fs0(n)
					!	write(*,*) '***********************************'
					endif
				enddo
 
	errn = fs0 - fs0d
	Call Integrate(e2,ABS(errn),delv,nk)

    IF (MOD(p,nac) == 0) THEN   ! Fission source extrapolation
        domiR = e2 / e1
        npos = MAXLOC(ABS(erro),1)
        IF (erro(npos) * errn(npos) < 0.0) domiR = -domiR
        fs0 = fs0 + domiR / (1._DP - domiR) * errn
    END IF
    e1 = e2                       ! Save integrated fission source error

	Call RelE(fs0,fs0d,fser,nk)     ! Search Maximum Point Wise Fission Source Relative Error
	Call RelEg(f0,f0d,fer,nk,ng)    ! Search Maximum Point Wise Flux Error
		!Write(*,'(3X,I5,2X,F13.6,2ES15.5)') p, Keff, fser, fer
      !  if (.TRUE.) Write(ounit,'(I5,F13.6,2ES15.5)') o, Keff, fser, fer
      !  if (.TRUE.) Write(*,'(I5,F13.6,2ES15.5)') o, Keff, fser, fer

	if ((fser < fserc) .AND. (fer < ferc)) Exit  ! if converge, exit.
End Do
	!do n=1, nk
	!	do g=1, ng
	!		if (n==50) then
	!			WRITE(*, *) " source=", Q(n,g,:)
	!		endif
	!	enddo
	!enddo
IF (p==nout+1) THEN
    maxi = .TRUE.
ELSE
    maxi = .FALSE.
END IF
!Write(*,*) '***********************************************'
!WRITE(*, '(A36,F9.6)') 'MULTIPLICATION EFFECTIVE (K-EFF) = ', Ke
!Write(*,*) '***********************************************'

END SUBROUTINE

            Subroutine Inner(f0,fx1,fy1,fz1,fx2,fy2,fz2,g,ng,nk,order,jo,ji,L0,Q,R2,P2,R4,P4,al,D,sigr,nxx,nyy,nzz,ix,iy,iz,delx,&
                             dely,delz,x_smax,x_smin,y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott) ! To Perform Inner Iterations
Implicit None
Integer, Parameter :: dp = Selected_Real_Kind(10,15)
Integer, Intent(in) :: g, ng, nk, order, nxx, nyy, nzz, x_east, x_west, y_north, y_south, z_top, z_bott
Integer, Dimension(nxx), Intent(in) :: x_smax, x_smin
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin
Real(dp), Dimension(nk), Intent(in) :: ix, iy, iz
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz 
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz
Real(dp), Dimension(nk,ng), Intent(in) :: D, sigr
Real(dp), Dimension(nk,ng,7), Intent(in) :: Q
Real(dp), Dimension(nk,ng,6,6),Intent(in) :: R2, P2, R4
Real(dp), Dimension(nk,ng,6,7),Intent(in) :: P4
Real(dp), Dimension(nk,ng,6), Intent(in) :: al
Real(dp), Dimension(nk,ng,3), Intent(inout) :: L0
Real(dp), Dimension(nk,ng,6), Intent(inout) :: jo
Real(dp), Dimension(nk,ng,6) :: ji
Real(dp), Dimension(nk,ng), Intent(inout) :: f0, fx1, fy1, fz1, fx2, fy2, fz2
! Local Variables
Real(dp), Dimension(6) :: MulJin, MulMom
Real(dp), Dimension(nk,ng,7) :: L      ! Transverse Leakage Moments(0, Lx1, Ly1, Lz1, Lx2, Ly2, Lz2)
Integer :: nin = 5      ! Maximum Inner Iteration
Integer :: n, k

! Jout Nodals' outgoing currents+flux  (X+, X-, Y+, Y-, Z+, Z-)
! Jin  Nodals' ingoing currents+source (X+, X-, Y+, Y-, Z+, Z-)

! Start Inner Iteration
Do n = 1, nin
    Do k = 1, nk ! Loop on k for Jin & QTL

!#######################   Calculate Ingoing Partial Currents   #######################
        Call Jin(ji,g,k,jo,al,ng,nk,nxx,nyy,nzz,ix,iy,iz,x_smax,x_smin,&
                 y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott)

!#######################   Update Transverse Leakage Moments   #######################
        Call QTL(L(k,g,:),g,k,ng,nk,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,xyz,x_smax,x_smin,&
                 y_smax,y_smin,L0,x_east,x_west,y_north,y_south,z_top,z_bott)
   End Do
    Do k = 1, nk
        if      (order == 2) Then
            ! Call MatVec2(MulJin,MulMom,R,P,ji(k,g,:),Q(k,g,:),ng,nk,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,D,sigr)        Call Nodal_Coup2(P,R,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
            Call MatVec2(MulJin,MulMom,R2(k,g,:,:),P2(k,g,:,:),ji(k,g,:),Q(k,g,:))
        Else if (order == 4) Then
            ! Call MatVec4(MulJin,MulMom,R,P,ji(k,g,:),Q(k,g,:),L(k,g,:),ng,nk,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,D,sigr)        Call Nodal_Coup2(P,R,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
            Call MatVec4(MulJin,MulMom,R4(k,g,:,:),P4(k,g,:,:),ji(k,g,:),Q(k,g,:),L(k,g,:))
        End if
        ! Call MatVec(R,ji,MulJin,ng,nk,g,k,6,6,6)
        ! Call MatVec(P,Q,MulMom,ng,nk,g,k,6,7,7)

!#######################   Update Outgoing Partial Currents   #######################
       jo(k,g,:) = MulJin + MulMom     !nod(xyz(ix(1),iy(1),iz(3)),2)jo(3) : Jout(y+) 2grp (1x,1y,3z)'node
    ! Update Ingoing Partial Currents in Adjacent Nodes
        ! Call Jin(ji,g,k,jo,al,ng,nk,nxx,nyy,nzz,ix,iy,iz,x_smax,x_smin,&
                 ! y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott)    ! Update Zeroth Transverse Leakages

!#######################   Update Zeroth Transverse Leakages   #######################
	!	if (k==50) then
	!		Write(*,*) 'BEFORE CALL Lxyz' 	
	!		write(*,*) 'jo=', jo(k,g,:) !Jin
	!		write(*,*) 'ji=', ji(k,g,:) !Q-L
	!	endif

        Call Lxyz(L0,k,g,ng,nk,jo,ji)
	!	if (k==50) then
	!		Write(*,*) 'AFTER CALL Lxyz' 	
	!		write(*,*) 'L0=', L0(k,g,:)
	!	endif

	
!#######################   Update Flux and Flux Moments   #######################
        if      (order == 2) Then
            Call Flux2(f0,fx1,fy1,fz1,fx2,fy2,fz2,k,g,ng,nk,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,sigr,Q,L0)
        Else if (order == 4) Then
            Call Flux4(f0,fx1,fy1,fz1,fx2,fy2,fz2,k,g,ng,nk,nxx,nyy,nzz,&
                       D,sigr,delx,dely,delz,ix,iy,iz,jo,ji,L0,Q,L)
        End if
!		if (k==50) then
!			Write(*,*) 'AFTER CALL Flux4=' 	
!			write(*,*) 'f0=', f0(k,g)
!		endif
    End Do
    Do k = nk, 1, -1 ! Loop on k for Jin & QTL
    ! To Calculate Ingoing Partial Currents from Neighborhod Nodes    
        Call Jin(ji,g,k,jo,al,ng,nk,nxx,nyy,nzz,ix,iy,iz,x_smax,x_smin,&
                 y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott)
    ! Update Transverse Leakage Moments
        Call QTL(L(k,g,:),g,k,ng,nk,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,xyz,x_smax,x_smin,&
                 y_smax,y_smin,L0,x_east,x_west,y_north,y_south,z_top,z_bott)
    End Do
    Do k = nk, 1, -1
        if (order == 2) Then
            ! Call MatVec2(MulJin,MulMom,R,P,ji(k,g,:),Q(k,g,:),ng,nk,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,D,sigr)        Call Nodal_Coup2(P,R,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
            Call MatVec2(MulJin,MulMom,R2(k,g,:,:),P2(k,g,:,:),ji(k,g,:),Q(k,g,:))
        Else if (order == 4) Then
            ! Call MatVec4(MulJin,MulMom,R,P,ji(k,g,:),Q(k,g,:),L(k,g,:),ng,nk,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,D,sigr)        Call Nodal_Coup2(P,R,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz)
            Call MatVec4(MulJin,MulMom,R4(k,g,:,:),P4(k,g,:,:),ji(k,g,:),Q(k,g,:),L(k,g,:))
        End if
        ! Call MatVec(R,ji,MulJin,ng,nk,g,k,6,6,6)
        ! Call MatVec(P,Q,MulMom,ng,nk,g,k,6,7,7)
    ! Update Outgoing Partial Currents
       jo(k,g,:) = MulJin+MulMom     !nod(xyz(ix(1),iy(1),iz(3)),2)jo(3) : Jout(y+) 2grp (1x,1y,3z)'node
    ! Update Ingoing Partial Currents in Adjacent Nodes
        ! Call Jin(ji,g,k,jo,al,ng,nk,nxx,nyy,nzz,ix,iy,iz,x_smax,x_smin,&
                 ! y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott)    ! Update Zeroth Transverse Leakages
    ! Update Zeroth Transverse Leakages
        Call Lxyz(L0,k,g,ng,nk,jo,ji)
    ! Update Flux and Flux Moments
        if (order == 2) Then
            Call Flux2(f0,fx1,fy1,fz1,fx2,fy2,fz2,k,g,ng,nk,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,sigr,Q,L0)
        Else if (order == 4) Then
            Call Flux4(f0,fx1,fy1,fz1,fx2,fy2,fz2,k,g,ng,nk,nxx,nyy,nzz,&
                       D,sigr,delx,dely,delz,ix,iy,iz,jo,ji,L0,Q,L)
        End if
    End Do
End Do

End Subroutine

Subroutine Output(ng,np,nmat,nx,ny,nz,nxx,nyy,nzz,delx,dely,delz,xD,x_sigr,&
                  x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,chi,zdiv,node,zpln,&
                  x_east,x_west,y_north,y_south,z_top,z_bott,mode)
implicit none
Integer, Parameter :: dp = Selected_Real_Kind(10,15), Ounit = 100, gm = 36
Integer, Intent(in) :: ng, np, nmat, nx, ny, nz, nxx, nyy, nzz, x_east, x_west, y_north, y_south, z_top, z_bott
Integer, Dimension(nz), Intent(in) :: zpln, zdiv
Integer, Dimension(np,nxx,nyy), Intent(in) :: node
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: x_sigs
Real(dp), Dimension(nmat,ng), Intent(in) :: x_siga, x_sigtr, x_sigf, x_nu_sigf, chi, xD, x_sigr
! Local Variables
Character(Len=2), Dimension(nxx, nyy) :: mmap
Character(Len=100), intent(in) :: mode
Integer, Dimension(ng) :: group
Integer :: i, j, k, g, h, lz, ztot, ip, ipr, kp, xs, xf

Open (Unit=Ounit, File='app/Output/NEM.out')

! #######################           Calculation Mode            #######################
Select Case(mode)
    Case('CBC')
        Write(Ounit,*)    '                       =========================='
        Write(Ounit,*)    '  Calculation Mode :    Critical Boron Concentration Search Calculation '
        Write(Ounit,*)    '  ----------------     =========================='
        Write(*,*)        '                       =========================='
        Write(*,*)        '  Calculation Mode :    Critical Boron Concentration Search Calculation '
        Write(*,*)        '  -----------------    =========================='
    Case('CBCTH')
        Write(Ounit,*)    '                       =========================='
        Write(Ounit,*)    '  Calculation Mode :    Critical Boron Concentration Search Calculation with TH feedbacks '
        Write(Ounit,*)    '  ----------------     =========================='
        Write(*,*)        '                       =========================='
        Write(*,*)        '  Calculation Mode :    Critical Boron Concentration Search Calculation with TH feedbacks '
        Write(*,*)        '  -----------------    =========================='
    Case('RODEJCT')
        Write(Ounit,*)    '                       =========================='
        Write(Ounit,*)    '  Calculation Mode :    Rod Ejection Calculation '
        Write(Ounit,*)    '  ----------------     =========================='
        Write(*,*)        '                       =========================='
        Write(*,*)        '  Calculation Mode :    Rod Ejection Calculation '
        Write(*,*)        '  -----------------    =========================='
    Case('THRODEJCT')
        Write(Ounit,*)    '                       =========================='
        Write(Ounit,*)    '  Calculation Mode :    Rod Ejection Calculation with TH feedbacks '
        Write(Ounit,*)    '  ----------------     =========================='
        Write(*,*)        '                       =========================='
        Write(*,*)        '  Calculation Mode :    Rod Ejection Calculation with TH feedbacks '
        Write(*,*)        '  -----------------    =========================='
    Case('Fixed Source')
        Write(Ounit,*)    '                       =========================='
        Write(Ounit,*)    '  Calculation Mode :    Fixed Source Calculation '
        Write(Ounit,*)    '  ----------------     =========================='
        Write(*,*)        '                       =========================='
        Write(*,*)        '  Calculation Mode :    Fixed Source Calculation '
        Write(*,*)        '  -----------------    =========================='
    Case('Adjoint')
        Write(Ounit,*)    '                       ====================='
        Write(Ounit,*)    '  Calculation Mode :    Adjoint Calculation '
        Write(Ounit,*)    '  ----------------     ====================='
        Write(*,*)        '                       ====================='
        Write(*,*)        '  Calculation Mode :    Adjoint Calculation '
        Write(*,*)        '  ----------------     ====================='
    Case('Forward')
        Write(Ounit,*)    '                       ====================='
        Write(Ounit,*)    '  Calculation Mode :    Forward Calculation '
        Write(Ounit,*)    '  ----------------     ====================='
        Write(*,*)        '                       ====================='
        Write(*,*)        '  Calculation Mode :    Forward Calculation '
        Write(*,*)        '  ----------------     ====================='
End Select
Write(ounit,*)
Write(ounit,*) '           --------------------------------------------'
Write(*,*)
Write(*,*) '           --------------------------------------------'

Do g=1,ng
    group(g)=g
End Do

! #######################   Macroscopic Cross Sections Output   #######################

Write(Ounit,*) '           >>>> Writing Macroscopic Cross Sections <<<<'
Write(Ounit,*) '           --------------------------------------------'
!Write(*,*)     '           >>>> Writing Macroscopic Cross Sections <<<<'
!Write(*,*)     '           --------------------------------------------'
    Do k= 1, nmat
        Write(Ounit,1009) k
        Write(Ounit,1011)'Group', 'Transport', 'Diffusion', 'Absorption', &
        'Removal', 'Nu*Fiss', 'Kap*Fis','Fiss. Spctr'
    !    Write(*,1009) k
   !     Write(*,1011)'Group', 'Transport', 'Diffusion', 'Absorption', &
     !   'Removal', 'Nu*Fiss', 'Kap*Fis','Fiss. Spctr'
        Do g= 1, ng
            Write(Ounit,1010) g, x_sigtr(k,g), xD(k,g), x_siga(k,g), &
            x_sigr(k,g), x_nu_sigf (k,g), x_sigf(k,g), chi(k,g)
     !       Write(*,1010) g, x_sigtr(k,g), xD(k,g), x_siga(k,g), &
    !        x_sigr(k,g), x_nu_sigf (k,g), x_sigf(k,g), chi(k,g)
        End Do
        Write(Ounit,*)'        --Scattering Matrix--'
        Write(Ounit,'(10X, A5, 20I11)') "G/G'", (group(g), g=1,ng)
     !   Write(*,*)'        --Scattering Matrix--'
      !  Write(*,'(10X, A5, 20I11)') "G/G'", (group(g), g=1,ng)
        Do g= 1, ng
            Write(ounit,1015) g, (x_sigs(k,g,h), h=1,ng)
      !      Write(*,1015) g, (x_sigs(k,g,h), h=1,ng)
        End Do
    End Do
! #######################   Core Geometry Output   #######################
Write(Ounit,*)
Write(ounit,*) '           -------------------------------'
Write(ounit,*) '           >>>>>Wrting Core Geometry<<<<<'
Write(ounit,*) '           -------------------------------'
Write(ounit,*)' Number of Assembly in x, y and z Directions Respectively :'
Write(ounit,*) nx, ny, nz
Write(ounit,*)' Number of Nodes in x, y and z Directions Respectively :'
Write(ounit,*) nxx, nyy, nzz
Write(ounit,*)
Write(ounit,1016) 'x','x'
Write(ounit,'(2X,10F7.2)')(delx(i), i=1,nxx)
Write(ounit,1016) 'y','y'
Write(ounit,'(2X,10F7.2)')(dely(j), j=1,nyy)
Write(*,*)
Write(*,*) '           -------------------------------'
Write(*,*) '           >>>>>Wrting Core Geometry<<<<<'
Write(*,*) '           -------------------------------'
Write(*,*)' Number of Assembly in x, y and z Directions Respectively :'
Write(*,*) nx, ny, nz
Write(*,*)' Number of Nodes in x, y and z Directions Respectively :'
Write(*,*) nxx, nyy, nzz
Write(*,*)
Write(*,1016) 'x','x'
Write(*,'(2X,10F7.2)')(delx(i), i=1,nxx)
Write(*,1016) 'y','y'
Write(*,'(2X,10F7.2)')(dely(j), j=1,nyy)
Write(*,*)
    ip = nxx/gm
    ipr = MOD(nxx,gm) - 1
    Do k= 1,np
        Do j = 1, nyy
            Do i = 1, nxx
                if (node(zpln(k),i,j) == 0) Then
                    mmap(i,j) = '  '
                Else
                    Write (mmap(i,j),'(I2)') node(zpln(k),i,j)
                    mmap(i,j) = TRIM(ADJUSTL(mmap(i,j)))
                End if
            End Do
        End Do

        Write(ounit,1017) k
        Write(*,1017) k
        xs = 1; xf = gm
        Do kp = 1, ip
            Write(ounit,'(6X,100I3)') (i, i = xs, xf)
    !        Write(*,'(6X,100I3)') (i, i = xs, xf)
            Do j= nyy, 1, -1
                Write(ounit,'(2X,I4,1X,100A3)') j, (mmap(i,j), i=xs, xf)
       !         Write(*,'(2X,I4,1X,100A3)') j, (mmap(i,j), i=xs, xf)
            End Do
            xs = xs + gm
            xf = xf + gm
        End Do

        Write(ounit,'(6X,100I3)') (i, i = xs, xs+ipr)
 !       Write(*,'(6X,100I3)') (i, i = xs, xs+ipr)
        if (xs+ipr > xs) Then
            Do j= nyy, 1, -1
                Write(ounit,'(2X,I4,1X,100A3)') j, (mmap(i,j), i=xs, xs+ipr)
     !           Write(*,'(2X,I4,1X,100A3)') j, (mmap(i,j), i=xs, xs+ipr)
            End Do
       End if
    End Do
    
    Write(ounit,*)
    Write(ounit,1018)
    Write(ounit,*) '-------------------------------------'
    Write(ounit,*) '  Plane Number     Planar Region    Delta-z'
    Write(*,*)
    Write(*,1018)
    Write(*,*) '-------------------------------------'
    Write(*,*) '  Plane Number     Planar Region    Delta-z'
    ztot = nzz
    Do k= nz, 1, -1
       Do lz= 1, zdiv(k)
            If (ztot == nzz) Then
                Write(ounit,'(I9, A6, I13, F15.2)') ztot, ' (Top)', zpln(k), delz(ztot)
    !            Write(*,'(I9, A6, I13, F15.2)') ztot, ' (Top)', zpln(k), delz(ztot)
            Else if (ztot == 1) Then
                 Write(ounit,'(I9, A9, I10, F15.2)') ztot, ' (Bottom)', zpln(k), delz(ztot)
       !          Write(*,'(I9, A9, I10, F15.2)') ztot, ' (Bottom)', zpln(k), delz(ztot)
            Else
                Write(ounit,'(I9, I19, F15.2)') ztot, zpln(k), delz(ztot)
          !      Write(*,'(I9, I19, F15.2)') ztot, zpln(k), delz(ztot)
            End if
            ztot = ztot - 1
        End Do
    End Do

    Write(ounit,*)
    Write(ounit,*) '  Boundary Conditions'
    Write(ounit,*) '  -------------------'
    !Write(*,*)
 !   Write(*,*) '  Boundary Conditions'
  !  Write(*,*) '  -------------------'
    if (x_west == 0) Then
       Write(ounit,*)' X-Directed West   : Zero Flux'
     !   Write(*,*)' X-Directed West   : Zero Flux'
    Else if (x_west == 1) Then
        Write(ounit,*)' X-Directed West   : Zero Incoming Current'
  !      Write(*,*)' X-Directed West   : Zero Incoming Current'
    Else
        Write(ounit,*)' X-Directed West   : Reflective'
    !    Write(*,*)' X-Directed West   : Reflective'
    End if
!
    if (x_east == 0) Then
        Write(ounit,*)' X-Directed East   : Zero Flux'
   !     Write(*,*)' X-Directed East   : Zero Flux'
   Else if (x_east == 1) Then
        Write(ounit,*)' X-Directed East   : Zero Incoming Current'
    !    Write(*,*)' X-Directed East   : Zero Incoming Current'
    Else
        Write(ounit,*)' X-Directed East   : Reflective'
      !  Write(*,*)' X-Directed East   : Reflective'
    End if

   if (y_north == 0) Then
        Write(ounit,*)' Y-Directed North  : Zero Flux'
  !      Write(*,*)' Y-Directed North  : Zero Flux'
    Else if (y_north == 1) Then
        Write(ounit,*)' Y-Directed North  : Zero Incoming Current'
     !   Write(*,*)' Y-Directed North  : Zero Incoming Current'
    Else
        Write(ounit,*)' Y-Directed North  : Reflective'
     !   Write(*,*)' Y-Directed North  : Reflective'
    End if

    if (y_south == 0) Then
        Write(ounit,*)' Y-Directed South  : Zero Flux'
  !      Write(*,*)' Y-Directed South  : Zero Flux'
    Else if (y_south == 1) Then
        Write(ounit,*)' Y-Directed South  : Zero Incoming Current'
     !   Write(*,*)' Y-Directed South  : Zero Incoming Current'
   Else
       Write(ounit,*)' Y-Directed South  : Reflective'
  !      Write(*,*)' Y-Directed South  : Reflective'
   End if

    if (z_bott == 0) Then
        Write(ounit,*)' Z-Directed Bottom : Zero Flux'
      !  Write(*,*)' Z-Directed Bottom : Zero Flux'
    Else if (z_bott == 1) Then
       Write(ounit,*)' Z-Directed Bottom : Zero Incoming Current'
     !   Write(*,*)' Z-Directed Bottom : Zero Incoming Current'
    Else
      Write(ounit,*)' Z-Directed Bottom : Reflective'
     !   Write(*,*)' Z-Directed Bottom : Reflective'
  End if

    if (z_top == 0) Then
        Write(ounit,*)' Z-Directed Top    : Zero Flux'
   !     Write(*,*)' Z-Directed Top    : Zero Flux'
    Else if (z_top == 1) Then
        Write(ounit,*)' Z-Directed Top    : Zero Incoming Current'
      !  Write(*,*)' Z-Directed Top    : Zero Incoming Current'
    Else
        Write(ounit,*)' Z-Directed Top    : Reflective'
   !     Write(*,*)' Z-Directed Top    : Reflective'
    End if

 1031 Format(2X, 'Calculation Mode : ', A40)
1009 Format(5X, 'Material', I3)
1010 Format(2X, I6, F13.6, 3F12.6, 3F13.6)
1011 Format(2X, A7, A12, A13, A12, A11, 2A13, A15)
1015 Format(10X, I3, F16.6, 20F12.6)
1016 Format(2X,A,'-Directed Nodes Divison (Delta-',A,')')
1017 Format(2X, 'Planar Region : ', I2)
1018 Format(2X, 'Planar Region Assignment to Planes.')

Write(ounit,*)
Write(ounit,*) ' ...Core Geometry is Sucessfully Read...'
Write(ounit,*) '    ---------------------------------'
Write(*,*)
Write(*,*) ' ...Core Geometry is Sucessfully Read...'
Write(*,*) '    ---------------------------------'

end subroutine

SUBROUTINE Output_cbcs (ng,nmat,rbcon,csigtr,csiga,cnuf,csigf,csigs)!    To print boron concentration for bc search                 

IMPLICIT NONE
Integer, Parameter :: dp = Selected_Real_Kind(10,15), Ounit = 100
Integer, Intent(in) :: ng, nmat
Real(dp), Dimension(nmat,ng), Intent(in) :: csiga, csigtr, csigf, cnuf
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: csigs
Real(dp), Intent(in) :: rbcon 

INTEGER :: i, g, h
INTEGER, DIMENSION(ng) :: group

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '           >>>> READING BORON CONCENTRATION FOR BC SEARCH <<<<'
WRITE(ounit,*) '           --------------------------------------------------'

!! BCON PRINT OPTION
    WRITE(ounit,1422) rbcon

    WRITE(ounit,*)
    WRITE(ounit,*) ' MATERIAL CX CHANGES PER PPM BORON CHANGES : '
    DO i= 1, nmat
       WRITE(ounit,1429) i
        WRITE(ounit,1431)'GROUP', 'TRANSPORT', 'ABSORPTION', &
        'NU*FISS', 'FISSION'
        DO g= 1, ng
            WRITE(ounit,1430) g, csigtr(i,g), csiga(i,g), &
            cnuf(i,g), csigf(i,g)
            group(g) = g
        END DO
        WRITE(ounit,*)'  --SCATTERING MATRIX--'
        WRITE(ounit,'(4X, A5, 20I11)') "G/G'", (group(g), g=1,ng)
        DO g= 1, ng
            WRITE(ounit,1435)g, (csigs(i,g,h), h=1,ng)
        END DO
    END DO

1422 FORMAT(2X, 'BORON CONCENTRATION REFERENCE :', F8.2)
1429 FORMAT(4X, 'MATERIAL', I3)
1431 FORMAT(2X, A8, A12, A13, A10, A14)
1430 FORMAT(2X, I6, 4E15.6)
1435 FORMAT(4X, I3, E17.5, 20E13.5)


WRITE(ounit,*)
WRITE(ounit,*) ' ...Critical Boron Search card is sucessfully read...'

END SUBROUTINE
SUBROUTINE Output_bcon (ng,nmat,bcon,rbcon,csigtr,csiga,cnuf,csigf,csigs)!    To read boron concentration

IMPLICIT NONE
Integer, Parameter :: dp = Selected_Real_Kind(10,15), Ounit = 100
Integer, Intent(in) :: ng, nmat
Real(dp), Dimension(nmat,ng), Intent(in) :: csiga, csigtr, csigf, cnuf
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: csigs
Real(dp), Intent(in) :: rbcon, bcon

INTEGER :: i, g, h
INTEGER, DIMENSION(ng) :: group

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '           >>>>       READING BORON CONCENTRATION        <<<<'
WRITE(ounit,*) '           --------------------------------------------------'


    WRITE(ounit,1221) bcon
    WRITE(ounit,1222) rbcon

    WRITE(ounit,*)
    WRITE(ounit,*) ' MATERIAL CX CHANGES PER PPM BORON CHANGES : '
    DO i= 1, nmat
       WRITE(ounit,1229) i
        WRITE(ounit,1231)'GROUP', 'TRANSPORT', 'ABSORPTION', &
        'NU*FISS', 'FISSION'
        DO g= 1, ng
            WRITE(ounit,1230) g, csigtr(i,g), csiga(i,g), &
            cnuf(i,g), csigf(i,g)
            group(g) = g
        END DO
        WRITE(ounit,*)'  --SCATTERING MATRIX--'
        WRITE(ounit,'(4X, A5, 20I11)') "G/G'", (group(g), g=1,ng)
        DO g= 1, ng
            WRITE(ounit,1235)g, (csigs(i,g,h), h=1,ng)
        END DO
    END DO
!END IF

1221 FORMAT(2X, 'BORON CONCENTRATION SET       :', F8.2)
1222 FORMAT(2X, 'BORON CONCENTRATION REFERENCE :', F8.2)
1229 FORMAT(4X, 'MATERIAL', I3)
1231 FORMAT(2X, A8, A12, A13, A10, A14)
1230 FORMAT(2X, I6, F13.6, F12.6, 2E14.5)
1235 FORMAT(4X, I3, E17.5, 20E13.5)

WRITE(ounit,*)
WRITE(ounit,*) ' ...Boron Concentration card is sucessfully read...'

END SUBROUTINE
SUBROUTINE Output_crod (fbmap,ng,nmat,nb,nx,ny,nxx,nyy,nzz,xdiv,ydiv,delz,nstep,bpos,&
						pos0,ssize,dsiga,dsigtr,dsigf,dnuf,dsigs,bmap)!    To print control rod position

IMPLICIT NONE
Integer, Parameter :: dp = Selected_Real_Kind(10,15), Ounit = 100
Integer, Intent(in) :: ng, nb, nmat, nx, ny, nxx, nyy, nzz
Real(dp), Dimension(nzz), Intent(in) :: delz
REAL(DP), Intent(in) :: nstep		! Number of steps
Integer, Dimension(nx), Intent(in) :: xdiv
Integer, Dimension(ny), Intent(in) :: ydiv
Real(dp), Dimension(nb), intent(in) :: bpos
REAL(DP), INTENT(IN) :: pos0, ssize			! Zero step position and step size
Real(dp), Dimension(nmat,ng), Intent(in) :: dsiga, dsigtr, dsigf, dnuf
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: dsigs

INTEGER, DIMENSION(nxx,nyy), intent(out) :: fbmap	! Radial control rod bank map (node wise)
INTEGER, DIMENSION(ng) :: group
INTEGER, DIMENSION(nb) :: bank
INTEGER, DIMENSION(nx, ny), intent(in) :: bmap       ! Radial control rod bank map (assembly wise)
INTEGER :: i, j, g, h, k
Real(dp) :: coreh       			! Core Height
INTEGER :: xtot, ytot, ly, lx
! Rod cusping option
INTEGER :: cusp = 0


WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '           >>>> READING CONTROL RODS INSERTION <<<<'
WRITE(ounit,*) '           --------------------------------------------'

! Calculate Core Height
coreh = 0._dp
Do k = 1, nzz
    coreh = coreh + delz(k)
End Do


!!! Check Control Rod Bank POSITION
DO i = 1, nb
    IF (bpos(i) > REAL(nstep)) THEN
        WRITE(ounit,1999) 'ERROR: POSITION OF CONTROL ROD BANK ', i, ' IS ', bpos(i), ' WHICH IS HIGHER THAN NUMBER OF STEPS.'
        STOP
    END IF
    IF (bpos(i) < 0.) THEN
        WRITE(ounit,1999) 'ERROR: POSITION OF CONTROL ROD BANK ', i, ' IS ', bpos(i), ' WHICH IS LOWER THAN 0.'
        STOP
    END IF
    IF (coreh < bpos(i)*ssize) THEN
        WRITE(ounit,1998) 'ERROR: CORE HEIGHT ', coreh, ' IS LOWER THAN CONTROL ROD POSITION ', bpos(i)*ssize+pos0
        WRITE(ounit,*) ' BANK NUMBER ', i
        STOP
    END IF
END DO
1999 FORMAT (2X, A, I2, A, F5.1, A)
1998 FORMAT (2X, A, F6.2, A, F6.2)

DO j = ny, 1, -1
    DO i = 1, nx
        IF (bmap(i,j) > nb) THEN
            WRITE(ounit,*) '  ERROR: BANK NUMBER ON CR BANK MAP IS GREATER THAN NUMBER OF BANK'
            STOP
        END IF
    END DO
END DO


!! CROD PRINT OPTION
    WRITE(ounit,1201) nb
    WRITE(ounit,1216) NINT(nstep)
    WRITE(ounit,1202) pos0
    WRITE(ounit,1203) ssize
    IF (cusp == 0) THEN
       WRITE(ounit,*) ' CR CUSPING CORRECTION       : NOT ACTIVE'
    ELSE
       WRITE(ounit,*) ' CR CUSPING CORRECTION       : ACTIVE'
    END IF

    DO i = 1, nb
        bank(i) = i
    END DO
    WRITE(ounit,*) ' INITIAL CONTROL ROD BANK POSITION (STEPS) : '
    WRITE(ounit,*) ' (0 means fully inserted) '
    WRITE(ounit, 1204)(bank(i), i = 1, nb)
    WRITE(ounit, 1205)(bpos(i), i = 1, nb)

    WRITE(ounit,*)
    WRITE(ounit,*) ' CONTROL ROD BANK MAP : '
    DO j = ny, 1, -1
        WRITE(ounit,'(100I3)' ) (bmap(i,j), i = 1, nx)
    END DO

    WRITE(ounit,*)
    WRITE(ounit,*) ' MATERIAL CX INCREMENT OR DECREMENT DUE TO CR INSERTION : '
    DO i= 1, nmat
       WRITE(ounit,1209) i
        WRITE(ounit,1211)'GROUP', 'TRANSPORT', 'ABSORPTION', &
        'NU*FISS', 'FISSION'
        DO g= 1, ng
            WRITE(ounit,1210) g, dsigtr(i,g), dsiga(i,g), &
            dnuf(i,g), dsigf(i,g)
            group(g) = g
        END DO
        WRITE(ounit,*)'  --SCATTERING MATRIX--'
        WRITE(ounit,'(4X, A5, 20I9)') "G/G'", (group(g), g=1,ng)
        DO g= 1, ng
            WRITE(ounit,1215)g, (dsigs(i,g,h), h=1,ng)
        END DO
    END DO

1201 FORMAT(2X, 'NUMBER OF CONTROL ROD BANK  :', I3)
1216 FORMAT(2X, 'MAX. NUMBER OF STEPS        :', I4)
1202 FORMAT(2X, 'FULLY INSERTED POSITION (cm): ', F4.1, ' (FROM BOTTOM OF THE CORE)')
1203 FORMAT(2X, 'STEP SIZE (cm)              : ', F8.4)
1204 FORMAT(2X, 10(:, 2X, 'Bank ', I2))
1205 FORMAT(10(:, 2X, F7.1), /)
1209 FORMAT(4X, 'MATERIAL', I3)
1211 FORMAT(2X, A7, A12, A12, 2A13)
1210 FORMAT(2X, I6, F13.6, F12.6, 2F13.6)
1215 FORMAT(4X, I3, F14.6, 20F10.6)

!!! Convert assembly wise CR bank map to node wise CR bank map
ytot = 0
DO j= 1, ny
    DO ly= 1, ydiv(j)
        ytot = ytot+1
        xtot = 0
        DO i= 1, nx
            DO lx= 1, xdiv(i)
                 xtot = xtot+1
                 fbmap(xtot, ytot) = bmap(i,j)
            END DO
        END DO
    END DO
END DO

WRITE(ounit,*)
WRITE(ounit,*) ' ...Control Rods Insertion card is sucessfully read...'

END SUBROUTINE
SUBROUTINE Output_ther (node_nf,ng,nk,nx,ny,nmat,nxx,nyy,&
					    y_smax,y_smin,xdiv,xsize,ydiv,ysize,powr,ppow,tin,&
						rf,tg,tc,ppitch,cf,pra,&
						cmflow,nfpin,ntem,stab)!    To read thermalhydraulics parameters input
Implicit None                                                                                                
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15), Ounit = 100
!! CXs Assigned to Materials
Integer, Intent(in) :: ng, nk, nx, ny, nmat, nxx, nyy, pra
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin 
Integer, Dimension(nx), Intent(in) :: xdiv
Integer, Dimension(ny), Intent(in) :: ydiv
Real(dp), Dimension(nx), Intent(in) :: xsize           ! Assembly Size & Division
Real(dp), Dimension(ny), Intent(in) :: ysize           ! Assembly Size & Division

! Thermal-hydraulics parameters
REAL(DP), intent(in) :: powr	! Reactor power for given geometry (watt)
REAL(DP), intent(in) :: ppow	! Reactor percent power in percent
REAL(DP), intent(in) :: tin	! coolant inlet temperature (kelvin)
REAL(DP), intent(in) :: rf, tg, tc, ppitch	! Fuel meat radius, gap thickness, clad thickness, and pin picth (m)
REAL(DP), intent(in) :: cf	! heat fraction deposited into coolant
REAL(DP), DIMENSION(nxx, nyy), intent(out) :: node_nf       ! Number of fuel pin per node
INTEGER, PARAMETER :: nm = 10      ! Fuel meat divided into 10 mesh
INTEGER, PARAMETER :: nt = nm + 2      ! two more mesh for gap and clad
REAL(DP), DIMENSION(nt) :: rdel	! mesh delta
REAL(DP), DIMENSION(nt) :: rpos	! mesh position
REAL(DP), intent(in) :: cmflow
INTEGER, intent(in) :: nfpin	! Number of fuel pin and guide tubes
INTEGER,  INTENT(IN) :: ntem    ! Number of temperature in steam table
REAL(DP), DIMENSION(ntem,pra), intent(in) :: stab  ! Steam table matrix
INTEGER :: i, j
INTEGER :: ly, lx, ytot, xtot
REAL(DP), DIMENSION(nx,ny) :: area
REAL(DP) :: barea, div
REAL(DP) :: dum

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '           >>>>   READING THERMAL-HYDRAULIC DATA   <<<<'
WRITE(ounit,*) '           --------------------------------------------'

! Check gap and clad thickness
IF (tg > 0.25 * rf) THEN
    WRITE(ounit,*) '  ERROR: GAP THICKNESS IS TO LARGE (> 0.25*rf)'
    STOP
END IF
IF (tc > 0.25 * rf) THEN
    WRITE(ounit,*) '  ERROR: CLADDING THICKNESS IS TO LARGE (> 0.25*rf)'
    STOP
END IF

! Calculate total coolant mass flow rate and number of fuel pins per node
barea = 0.
DO j = 1, ny 
    DO i = 1, nx
        area(i,j) = xsize(i)*ysize(j)             ! assembly area
        IF (area(i,j) > barea) barea = area(i,j)  ! barea => largest assembly area for ref.
    END DO
END DO

node_nf = 0.

ytot = 0
DO j= 1, ny
    DO ly= 1, ydiv(j)
      ytot = ytot+1
        xtot = 0
        DO i= 1, nx
            DO lx= 1, xdiv(i)
                xtot = xtot+1
                IF ((xtot >= y_smin(ytot)) .AND. (xtot <= y_smax(ytot) )) THEN
                    div = REAL (ydiv(j) * xdiv(i))            ! Number of nodes in current assembly
                    node_nf(xtot,ytot) = area(i,j) * REAL(nfpin) / (barea * div)   ! Number of fuel pin for this node
                END IF
            END DO
        END DO
    END DO
END DO
!DO i= 1, nx
!					write(*,*) 'i=', i, 'xsize(i)=', xsize(i) 
!enddo
!DO j= 1, ny
!					write(*,*) 'j=', j, 'ysize(i)=', ysize(j) 
!enddo
!DO j= 1, ny
!DO i= 1, nx
	!				write(*,*) 'node_nf(i,j)=', node_nf(:,:) 
!enddo
!enddo

! Calculate fuel pin mesh delta and position
!nm = 10      ! Fuel meat divided into 10 mesh
!nt = nm + 2  ! two more mesh for gap and clad
!dum = rf / REAL(nm)

!DO i = 1, nm
 !   rdel(i) = dum
!END DO

! Fuel pin mesh size
!rdel(nm+1) = tg
!rdel(nm+2) = tc

! Fuel pin mesh position
!rpos(1) = 0.5 * rdel(1)
!DO i = 2, nt
!    rpos(i) = rpos(i-1) + 0.5 * (rdel(i-1) + rdel(i))
!END DO


!! THER PRINT OPTION
    WRITE(ounit,1309) ppow
    WRITE(ounit,1301) powr
    WRITE(ounit,1302) tin
    WRITE(ounit,1303) cmflow
    WRITE(ounit,1304) rf
    WRITE(ounit,1305) tg
    WRITE(ounit,1306) tc
    WRITE(ounit,1310) ppitch
    WRITE(ounit,1307) cf

WRITE(ounit,*)
WRITE(ounit,*) ' ...Thermal-hydraulic Card is sucessfully read...'

1309 FORMAT(2X, 'REACTOR PERCENT POWER (%)                : ', F12.5)
1301 FORMAT(2X, 'REACTOR POWER (Watt)                     : ', ES12.4)
1302 FORMAT(2X, 'COOLANT INLET TEMPERATURE (Kelvin)       : ', ES12.4)
1303 FORMAT(2X, 'FUEL ASSEMBLY MASS FLOW RATE (Kg/s)      : ', ES12.4)
1304 FORMAT(2X, 'FUEL MEAT RADIUS (m)                     : ', ES12.4)
1305 FORMAT(2X, 'GAP THICKNESS (m)                        : ', ES12.4)
1306 FORMAT(2X, 'CLAD THICKNESS (m)                       : ', ES12.4)
1310 FORMAT(2X, 'PIN PITCH(m)                             : ', ES12.4)
1307 FORMAT(2X, 'FRACTION OF HEAT DEPOSITED IN COOL.      : ', ES12.4)

END SUBROUTINE
SUBROUTINE Output_ftem (ftem,ng,nk,nmat,cftem,rftem,fsigs,&
					 fsiga,fsigtr,fsigf,fnuf)!    To read fuel temperature
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15), Ounit = 100
Integer, Intent(in) :: ng, nk, nmat
!! CXs Assigned to Materials
REAL(DP), intent(in) :: cftem, rftem
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: fsigs
Real(dp), Dimension(nk), Intent(out) :: ftem
Real(dp), Dimension(nmat,ng), Intent(in) :: fsiga, fsigtr, fsigf, fnuf
INTEGER :: i, g, h
INTEGER, DIMENSION(ng) :: group
INTEGER :: bther = 1

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '           >>>>      READING FUEL TEMPERATURE      <<<<'
WRITE(ounit,*) '           --------------------------------------------'

! ASSIGN CFTEM to FTEM
ftem = cftem

! Read CX changes fuel temperature change

!! FTEM PRINT OPTION

    IF (bther == 0) THEN
        WRITE(ounit,1241) cftem
    ELSE
        WRITE(ounit,1256) cftem
    END IF
    WRITE(ounit,1242) rftem

    WRITE(ounit,*)
    WRITE(ounit,*) ' MATERIAL CX CHANGES PER FUEL TEMPERATURE CHANGES : '
    DO i= 1, nmat
       WRITE(ounit,1249) i
        WRITE(ounit,1251)'GROUP', 'TRANSPORT', 'ABSORPTION', &
        'NU*FISS', 'FISSION'
        DO g= 1, ng
            WRITE(ounit,1250) g, fsigtr(i,g), fsiga(i,g), &
            fnuf(i,g), fsigf(i,g)
            group(g) = g
        END DO
        WRITE(ounit,*)'  --SCATTERING MATRIX--'
        WRITE(ounit,'(4X, A5, 20I11)') "G/G'", (group(g), g=1,ng)
        DO g= 1, ng
            WRITE(ounit,1255)g, (fsigs(i,g,h), h=1,ng)
        END DO
    END DO

1241 FORMAT(2X, 'AVERAGE FUEL TEMPERATURE   :', F6.2)
1242 FORMAT(2X, 'FUEL TEMPERATURE REFERENCE :', F6.2)
1249 FORMAT(4X, 'MATERIAL', I3)
1251 FORMAT(2X, A8, A12, A13, A10, A14)
1250 FORMAT(2X, I6, E14.5, E13.5, 2E14.5)
1255 FORMAT(4X, I3, E17.5, 20E13.5)
1256 FORMAT(2X, 'AVERAGE FUEL TEMPERATURE   :', F6.2, '  (NOT USED)')

WRITE(ounit,*)
WRITE(ounit,*) ' ...Fuel Temperature is card sucessfully read...'

END SUBROUTINE
SUBROUTINE Output_mtem (mtem,ng,nk,nmat,cmtem,rmtem,msigs,&
					 msiga,msigtr,msigf,mnuf)!    To read moderator temperature
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15), Ounit = 100
Integer, Intent(in) :: ng, nk, nmat
!! CXs Assigned to Materials
REAL(DP), intent(in) :: cmtem, rmtem
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: msigs
Real(dp), Dimension(nk), Intent(out) :: mtem
Real(dp), Dimension(nmat,ng), Intent(in) :: msiga, msigtr, msigf, mnuf
INTEGER :: i, g, h
INTEGER, DIMENSION(ng) :: group
INTEGER :: bther = 1

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '           >>>>   READING MODERATOR TEMPERATURE    <<<<'
WRITE(ounit,*) '           --------------------------------------------'

! ASSIGN CMTEM to MTEM
mtem = cmtem

! Read CX changes per moderator temperature change
!! MTEM PRINT OPTION

    IF (bther == 0) THEN
        WRITE(ounit,1261) cmtem
    ELSE
        WRITE(ounit,1276) cmtem
    END IF
    WRITE(ounit,1262) rmtem

    WRITE(ounit,*)
    WRITE(ounit,*) ' MATERIAL CX CHANGES PER MODERATOR TEMPERATURE CHANGES : '
    DO i= 1, nmat
       WRITE(ounit,1269) i
        WRITE(ounit,1271)'GROUP', 'TRANSPORT', 'ABSORPTION', &
        'NU*FISS', 'FISSION'
        DO g= 1, ng
            WRITE(ounit,1270) g, msigtr(i,g), msiga(i,g), &
            mnuf(i,g), msigf(i,g)
            group(g) = g
        END DO
        WRITE(ounit,*)'  --SCATTERING MATRIX--'
        WRITE(ounit,'(4X, A5, 20I11)') "G/G'", (group(g), g=1,ng)
        DO g= 1, ng
            WRITE(ounit,1275)g, (msigs(i,g,h), h=1,ng)
        END DO
    END DO

1261 FORMAT(2X, 'AVERAGE MODERATOR TEMPERATURE   :', F6.2)
1262 FORMAT(2X, 'MODERATOR TEMPERATURE REFERENCE :', F6.2)
1269 FORMAT(4X, 'MATERIAL', I3)
1271 FORMAT(2X, A8, A12, A13, A10, A14)
1270 FORMAT(2X, I6, E14.5, E13.5, 2E14.5)
1275 FORMAT(4X, I3, E17.5, 20E13.5)
1276 FORMAT(2X, 'AVERAGE MODERATOR TEMPERATURE   :', F6.2, '  (NOT USED)')



WRITE(ounit,*)
WRITE(ounit,*) ' ...Moderator Temperature Card is sucessfully read...'


END SUBROUTINE
SUBROUTINE Output_cden (cden,ng,nk,nmat,ccden,rcden,lsigs,&
					 lsiga,lsigtr,lsigf,lnuf)!    To read Coolant Density
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15), Ounit = 100
Integer, Intent(in) :: ng, nk, nmat
!! CXs Assigned to Materials
REAL(DP), intent(in) :: ccden, rcden
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: lsigs
Real(dp), Dimension(nk), Intent(out) :: cden
Real(dp), Dimension(nmat,ng), Intent(in) :: lsiga, lsigtr, lsigf, lnuf
INTEGER :: i, g, h
INTEGER, DIMENSION(ng) :: group
INTEGER :: bther = 1

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '           >>>>       READING COOLANT DENSITY      <<<<'
WRITE(ounit,*) '           --------------------------------------------'


!ASSIGN CCDEN TO CDEN
cden = ccden

! Read CX changes per Coolant Density change

!! CDEN PRINT OPTION
    IF (bther == 0) THEN
        WRITE(ounit,1361) ccden
    ELSE
        WRITE(ounit,1376) ccden
    END IF
    WRITE(ounit,1362) rcden

    WRITE(ounit,*)
    WRITE(ounit,*) ' MATERIAL CX CHANGES PER COOLANT DENSITY CHANGES : '
    DO i= 1, nmat
       WRITE(ounit,1369) i
        WRITE(ounit,1371)'GROUP', 'TRANSPORT', 'ABSORPTION', &
        'NU*FISS', 'FISSION'
        DO g= 1, ng
            WRITE(ounit,1370) g, lsigtr(i,g), lsiga(i,g), &
            lnuf(i,g), lsigf(i,g)
            group(g) = g
        END DO
        WRITE(ounit,*)'  --SCATTERING MATRIX--'
        WRITE(ounit,'(4X, A5, 20I11)') "G/G'", (group(g), g=1,ng)
        DO g= 1, ng
            WRITE(ounit,1375)g, (lsigs(i,g,h), h=1,ng)
        END DO
    END DO


1361 FORMAT(2X, 'AVERAGE COOLANT DENSITY   :', F8.4)
1362 FORMAT(2X, 'COOLANT DENSITY REFERENCE :', F8.4)
1369 FORMAT(4X, 'MATERIAL', I3)
1371 FORMAT(2X, A8, A12, A13, A10, A14)
1370 FORMAT(2X, I6, E14.5, E13.5, 2E14.5)
1375 FORMAT(4X, I3, E17.5, 20E13.5)
1376 FORMAT(2X, 'AVERAGE COOLANT DENSITY   :', F8.4, '  (USED AS GUESS)')



WRITE(ounit,*)
WRITE(ounit,*) ' ...Coolant Density Card is sucessfully read...'


END SUBROUTINE
SUBROUTINE Output_ejct (mdir,ng,nk,nb,nmat,bpos,fbpos,tmove,bspeed,&
					ibeta,lamb,velo,ttot,tstep1,tstep2,tdiv)!    To read rod ejection input

IMPLICIT NONE
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15), Ounit = 100
Integer, Intent(in) :: ng, nk, nmat, nb
! Crod changes
REAL(DP), DIMENSION(nb), intent(in) :: bpos  	! CR bank position
REAL(DP), DIMENSION(nb), intent(in) :: fbpos    ! Final CR bank position
REAL(DP), DIMENSION(nb), intent(in) :: tmove    ! Time when CR bank starts moving
REAL(DP), DIMENSION(nb), intent(in) :: bspeed   ! CR bank movement speed

! Transient parameters
INTEGER, PARAMETER :: nf = 6                       ! Number of delaye dneutron precusor family
REAL(DP), DIMENSION(nf), intent(in) :: ibeta, lamb                 ! beta (delayed neutron fraction) and precusor decay constant
REAL(DP), DIMENSION(ng), intent(in) :: velo            ! Neutron velocity
REAL(DP), intent(in) :: ttot				! TOTAL SIMULATION TIME
REAL(DP), intent(in) :: tstep1				! FIRST TIME STEP
REAL(DP), intent(in) :: tstep2				! SECOND TIME STEP
REAL(DP), intent(in) :: tdiv				! WHEN SECOND TIME STEP APPLY
LOGICAL :: tranw = .FALSE.					! To activate unconverged  outer iteration warning

INTEGER, DIMENSION(nb), intent(out)  :: mdir  	! To indicate CR movement direction (0=do not move, 1=down, 2 = up)


INTEGER :: ln   !Line number
INTEGER :: ios  ! IOSTAT status

INTEGER :: i, g
INTEGER :: popt
INTEGER, DIMENSION(nb) :: bank
CHARACTER(LEN=4) :: cnb         ! number of bank (character type)

WRITE(ounit,*)
WRITE(ounit,*)
WRITE(ounit,*) '           >>>>     READING ROD EJECTION DATA      <<<<'
WRITE(ounit,*) '           --------------------------------------------'

! Read Final CR bank position, time to start, and speed
DO i = 1, nb
    WRITE (cnb,'(I4)') nb
    cnb = TRIM(ADJUSTL(cnb))
    IF (ABS(fbpos(i)-bpos(i)) < 1.e-5_DP) THEN
        mdir(i) = 0
    ELSE IF (fbpos(i)-bpos(i) > 1.e-5_DP) THEN
        mdir(i) = 2
    ELSE
        mdir(i) = 1
    END IF
END DO

! Read time for CR to be ejected

! ttot must be bigger than tstep1 and tstep2
IF ((ttot < tstep1) .OR. (ttot < tstep2)) THEN
    WRITE(ounit,*) 'ERROR: TOTAL SIMULATION TIME SHALL BE GREATER THAN TIME STEPS'
    STOP
END IF

! tdiv must be bigger than tstep1
IF (tdiv < tstep1) THEN
    WRITE(ounit,*) 'ERROR: THE TIME WHEN SECOND TIME STEP STARTS SHALL BE GREATER THAN FIRST TIME STEP'
    STOP
END IF

! tdiv must be less than ttot
IF (tdiv > ttot) THEN
    WRITE(ounit,*) 'ERROR: THE TIME WHEN SECOND TIME STEP STARTS SHALL BE LESS THAN TOTAL TIME'
    STOP
END IF

! number of steps shall be less than 10,000
IF (NINT(tdiv/tstep1)+NINT((ttot-tdiv)/tstep2) > 10000) THEN
    WRITE(ounit,*) 'ERROR: NUMBER OF TOTAL TIME STEPS ARE MORE THAN 10,000'
    STOP
END IF

! Read beta (delayed neutron fraction)

! Read precusor decay constant

! Read neutron velocity


!! EJCT PRINT OPTION

    DO i = 1, nb
        bank(i) = i
    END DO
    WRITE(ounit, 1294)(bank(i), i = 1, nb)
    WRITE(ounit, 1295)(fbpos(i), i = 1, nb)
    WRITE(ounit, 1281)(tmove(i), i = 1, nb)
    WRITE(ounit, 1282)(bspeed(i), i = 1, nb)

    WRITE(ounit,*)
    WRITE(ounit,*) ' TIME PARAMETERS IN SECONDS : '
    WRITE(ounit,1297) ttot
    WRITE(ounit,1298) tstep1
    WRITE(ounit,1299) tstep2
    WRITE(ounit,1300) tdiv

    WRITE(ounit,*)
    WRITE(ounit,*) ' DELAYED NEUTRON FRACTION : '
    WRITE(ounit,'(100F11.5)') (iBeta(i), i = 1, nf)

    WRITE(ounit,*)
    WRITE(ounit,*) ' PRECUSOR DECAY CONSTANT (1/s): '
    WRITE(ounit,'(100F11.5)') (lamb(i), i = 1, nf)

    WRITE(ounit,*)
    WRITE(ounit,*) ' NEUTRON VELOCITY (cm/s) : '
    WRITE(ounit,'(100ES15.5)') (velo(g), g = 1, ng)

WRITE(ounit,*)
WRITE(ounit,*) ' ...Rod Ejection Card is sucessfully read...'

1294 FORMAT(25X, 99(:, 2X, 'Bank ', I2))
1295 FORMAT(2X, 'FINAL BANK POS. (STEP)', 99(:, 2X, F7.1), /)
1281 FORMAT(2X, 'STARTS MOVE (SECOND)  ', 99(:, 2X, F7.1), /)
1282 FORMAT(2X, 'SPEED (STEP/SECOND)   ', 99(:, 2X, F7.1), /)
1297 FORMAT(4X, 'TOTAL SIMULATION TIME         : ', F6.2)
1298 FORMAT(4X, 'FIRST TIME STEP               : ', F6.4)
1299 FORMAT(4X, 'SECOND TIME STEP              : ', F6.4)
1300 FORMAT(4X, 'WHEN SECOND TIME STEP APPLY?  : ', F6.2)

END SUBROUTINE

! ####################################################################

! #######################   Geometry Subroutines   #######################
Subroutine Nodxyz(nxx,nyy,nzz,nx,ny,nz,xdiv,ydiv,zdiv) ! To Calculate Number of Nodes in x,y & z Directions
Implicit None
Integer, Intent(in) :: nx, ny, nz
Integer, Dimension(nx), Intent(in) ::  xdiv            ! Assembly Size & Division
Integer, Dimension(ny), Intent(in) ::  ydiv            ! Assembly Size & Division
Integer, Dimension(nz), Intent(in) ::  zdiv            ! Assembly Size & Division

Integer, Intent(out) :: nxx, nyy, nzz
! Local Variables
Integer :: i, j, k

! x-direction
nxx=0
Do i= 1,nx
    nxx = nxx+xdiv(i)
End Do
! y-direction
nyy=0
Do j= 1,ny
    nyy = nyy+ydiv(j)
End Do
! z-direction
nzz=0
Do k= 1,nz
    nzz = nzz+zdiv(k)
End Do

End Subroutine

Subroutine Asmblg_delta(node,mnum,delx,dely,delz,nmat,np,nx,ny,nz,nxx,nyy,nzz,&
                        xsize,ysize,zsize,xdiv,ydiv,zdiv,zpln,assm) ! To assign Material into nodes & Calculate Delta x,y & z (node sizes)
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: np                              ! Number of planars
Integer, Intent(in) :: nmat, nx, ny, nz, nxx, nyy, nzz
Integer, Dimension(nx), Intent(in) ::  xdiv            ! Assembly Size & Division
Integer, Dimension(ny), Intent(in) ::  ydiv            ! Assembly Size & Division
Integer, Dimension(nz), Intent(in) ::  zdiv            ! Assembly Size & Division
Real(dp), Dimension(nx), Intent(in) :: xsize           ! Assembly Size & Division
Real(dp), Dimension(ny), Intent(in) :: ysize           ! Assembly Size & Division
Real(dp), Dimension(nz), Intent(in) :: zsize           ! Assembly Size & Division
Integer, Dimension(nz), Intent(in) :: zpln             ! Planar Assignment to Z Direction
Integer, Dimension(np,nx,ny), Intent(in) :: assm       ! Material Assignment into Assembly
Integer, Dimension(nxx,nyy,nzz), Intent(out) :: mnum
Integer, Dimension(np,nxx,nyy), Intent(out) :: node    ! Material Assignment into Nodes
Real(dp), Dimension(nxx), Intent(out) :: delx
Real(dp), Dimension(nyy), Intent(out) :: dely
Real(dp), Dimension(nzz), Intent(out) :: delz
! Local Variables
Integer :: i, j, k, lx, ly, lz, xtot, ytot, ztot
Real(dp) :: div
Integer, Parameter :: Ounit = 100   !output

Do k = 1, nz
    if (zpln(k) > np) Then
        Write(Ounit,'(2X,A15,I3,A35)') 'Error: Planar ', &
        zpln(k), ' is Greater than Number of Planar'
        Write(*,'(2X,A15,I3,A35)') 'Error: Planar ', &
        zpln(k), ' is Greater than Number of Planar'
        Stop
    End if
    if (zpln(k) < 1) Then
        Write(Ounit,'(2X,A)') 'Error: Planar Should Be at Least Equal 1'
        Write(*,'(2X,A)') 'Error: Planar Should Be at Least Equal 1'
        Stop
    End if
End Do

ztot = 0
Do k= 1, nz
    Do lz= 1, zdiv(k)
        ztot = ztot+1
        ytot = 0
        Do j= 1, ny
            Do ly= 1, ydiv(j)
                ytot = ytot+1
                xtot = 0
                Do i= 1, nx
                    Do lx= 1, xdiv(i)
                        xtot = xtot+1
                        mnum(xtot,ytot,ztot) = assm(zpln(k),j,i)
                        if (mnum(xtot,ytot,ztot) > nmat) Then
                            Write(ounit,'(2X,A17,I3,A37)') 'Error: Material ', &
                            mnum(xtot,ytot,ztot), ' is Greater than Number of Materials'
                            Write(*,'(2X,A17,I3,A37)') 'Error: Material ', &
                            mnum(xtot,ytot,ztot), ' is Greater than Number of Materials'
                            Stop
                        End if
                        if (mnum(xtot, ytot, ztot) < 0) Then
                            Write(ounit,'(2X,A)') 'Error: Negative Material Found'
                            Write(*,'(2X,A)') 'Error: Negative Material Found'
                            Stop
                        End if
                        node(zpln(k),xtot,ytot) = assm(zpln(k),j,i)
                    End Do
                End Do
            End Do
        End Do
    End Do
End Do

!Delta x
xtot=0
Do i= 1,nx
    div = xsize(i)/Real(xdiv(i))
    Do lx= 1, xdiv(i)
    xtot = xtot+1
    delx(xtot) = div
    End Do
End Do
!Delta y
ytot=0
Do j= 1,ny
    div = ysize(j)/Real(ydiv(j))
    Do ly= 1, ydiv(j)
    ytot = ytot+1
    dely(ytot) = div
    End Do
End Do
!Delta z
ztot=0
Do k= 1,nz
    div = zsize(k)/Real(zdiv(k))
    Do lz= 1, zdiv(k)
    ztot = ztot+1
    delz(ztot) = div
    End Do
End Do

End Subroutine

Subroutine Stagg(y_smax,y_smin,x_smax,x_smin,nxx,nyy,nzz,mnum) ! To Index non zero material for staggered mesh
Implicit None
Integer, Intent(in) :: nxx, nyy, nzz
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: mnum
Integer, Dimension(nxx), Intent(out) :: x_smax, x_smin ! imax and imin along x Direction for Staggered Nodes
Integer, Dimension(nyy), Intent(out) :: y_smax, y_smin ! imax and imin along y Direction for Staggered Nodes
! Local Variables
Integer :: i, j

! Along y Direction
Do j= 1, nyy
    y_smin(j) = nxx
    Do i = 1, nxx
        if (mnum(i,j,1) /= 0) Then
            y_smin(j) = i
            Exit
        End if
    End Do
End Do
Do j= 1, nyy
    y_smax(j) = 0
    Do i = nxx, 1, -1
        if (mnum(i,j,1) /= 0) Then
            y_smax(j) = i
            Exit
        End if
    End Do
End Do
! Along x Direction
Do i= 1, nxx
    x_smin(i) = nyy
    Do j = 1, nyy
        if (mnum(i,j,1) /= 0) Then
            x_smin(i) = j
            Exit
        End if
    End Do
End Do
Do i= 1, nxx
    x_smax(i) = 0
    Do j = nyy, 1, -1
        if (mnum(i,j,1) /= 0) Then
            x_smax(i) = j
            Exit
        End if
    End Do
End Do

End Subroutine

Subroutine Nod(nk,nyy,nzz,y_smax,y_smin) ! To Set Number of Nodes
Implicit None
Integer, Intent(in) :: nyy, nzz
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin
Integer, Intent(out):: nk
! Local Variables
Integer :: i, j, k

nk = 0
Do k = 1, nzz
    Do j = 1, nyy
        Do i = y_smin(j), y_smax(j)
            nk = nk + 1
        End Do
    End Do
End Do

End Subroutine

Subroutine PosiNod(ix,iy,iz,xyz,mat,y_smax,y_smin,nk,nxx,nyy,nzz,mnum) ! To Set ix,iy,iz & xyz
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: nk
Integer, Intent(in) :: nxx, nyy, nzz
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: mnum
Integer, Dimension(nxx,nyy,nzz), Intent(out) :: xyz
Real(dp), Dimension(nk), Intent(out) :: ix, iy, iz
Integer, Dimension(nk), Intent(out) :: mat
! Local Variables
Integer :: i, j, k, n

n = 0
Do k = 1, nzz      ! Don't change the order as it would affect th_upd and th_trans
    Do j = nyy, 1, -1
        Do i = y_smin(j), y_smax(j)
            n = n + 1
            ix(n) = i
            iy(n) = j
            iz(n) = k
            xyz(i,j,k) = n
            mat(n) = mnum(i,j,k)
        End Do
    End Do
End Do

End Subroutine

Subroutine Deltv(delv,nk,nxx,nyy,nzz,delx,dely,delz,ix,iy,iz) ! To Calculate Nodes' Volume
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: nk, nxx, nyy, nzz
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk), Intent(in) :: ix, iy, iz
Real(dp), Dimension(nk), Intent(out) :: delv
! Local Variables
Integer :: i

Do i = 1, nk
    delv(i) = delx(ix(i)) * dely(iy(i)) * delz(iz(i))
End Do

End Subroutine
! #########################################################################

Subroutine xD_xsigr(xD,x_sigr,x_siga,x_sigs,x_sigtr,ng,nk) ! To Update Diffusion Coefficient and Removal XS
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: ng, nk
Real(dp), Dimension(nk,ng,ng), Intent(in) :: x_sigs
Real(dp), Dimension(nk,ng), Intent(in) :: x_sigtr, x_siga
Real(dp), Dimension(nk,ng), Intent(out):: x_sigr, xD
! Local Variables
Integer :: k, g, h
Real(dp) :: som

Do k = 1, nk
    Do g = 1, ng
        xD(k,g) = 1./(3.*x_sigtr(k,g))
        som = 0.0
        Do h= 1, ng
            if (g/=h) som = som + x_sigs(k,g,h)
        End Do
        x_sigr(k,g) =  x_siga(k,g) + som
    End Do
End Do  

End Subroutine

Subroutine Init(keff,jo,ji,L0,al,f0,fx1,fy1,fz1,fx2,fy2,fz2,ng,nk) ! To provide Initial Guess Values for Keff, Jin, Jout, Flux & Flux Moments
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: ng, nk
Real(dp), Intent(out):: Keff
Real(dp), Dimension(nk,ng,6), Intent(out) :: jo, ji, al
Real(dp), Dimension(nk,ng,3), Intent(out) :: L0
Real(dp), Dimension(nk,ng), Intent(out) :: f0, fx1, fy1, fz1, fx2, fy2, fz2
! Local Variables
Integer :: g, k

Keff = 1._dp
Do g= 1, ng
    Do k = 1, nk
        jo(k,g,:) = 1._dp
        ji(k,g,:) = 1._dp
        al(k,g,:) = 0.0     ! Default alpha
        Call Lxyz(L0,k,g,ng,nk,jo,ji)
        f0(k,g)  = 1._dp
        fx1(k,g) = 1._dp
        fy1(k,g) = 1._dp
        fz1(k,g) = 1._dp
        fx2(k,g) = 1._dp
        fy2(k,g) = 1._dp
        fz2(k,g) = 1._dp
    End Do
End Do

End Subroutine

Subroutine XS_updt(sigs,siga,sigtr,sigf,nu_sigf,D,sigr,ng,nk,nmat,&
                   mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf) ! To Update Current XS to Base XS
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: ng, nk, nmat
!! CXs Assigned to Materials
Integer,  Dimension(nk), Intent(in) :: mat    ! Materials
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: x_sigs
Real(dp), Dimension(nmat,ng), Intent(in) :: x_siga, x_sigtr, x_sigf, x_nu_sigf
!! CXs Assigned to Nodes
Real(dp), Dimension(nk,ng,ng), Intent(out) :: sigs
Real(dp), Dimension(nk,ng), Intent(out) :: siga, sigtr, sigf, nu_sigf, sigr, D
Integer :: k

Do k = 1, nk
    sigs (k,1:,1:) = x_sigs (mat(k),1:,1:)
    siga (k,1:)    = x_siga (mat(k),1:)
    sigtr (k,1:)   = x_sigtr(mat(k),1:)
    sigf (k,1:)    = x_sigf (mat(k),1:)
    nu_sigf (k,1:) = x_nu_sigf (mat(k),1:)
End Do

    Call D_sigr(D,sigr,siga,sigs,sigtr,ng,nk)

End Subroutine
Subroutine XS_updtf(sigs,siga,sigtr,sigf,nu_sigf,D,sigr,ng,nk,nmat,mat,x_sigs,x_siga,x_sigtr,x_sigf,&
					 x_nu_sigf,fsiga,fsigtr,fsigf,fnuf,&
					 fsigs,rftem,f0,fz1,fz2,jo,ji,nxx,nyy,nzz,&
					 delz,xyz,xftem) ! To Update Current XS to Base XS
					 
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: ng, nk, nmat, nxx, nyy, nzz
!! CXs Assigned to Materials
Integer,  Dimension(nk), Intent(in) :: mat    ! Materials
Real(dp), Dimension(nzz), Intent(in) :: delz
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: x_sigs
Real(dp), Dimension(nmat,ng), Intent(in) :: x_siga, x_sigtr, x_sigf, x_nu_sigf
REAL(DP), DIMENSION(nk), INTENT(IN) :: xftem  ! Provided fuel temperature
Real(dp), Dimension(nmat,ng), Intent(in) :: fsiga, fsigtr, fsigf, fnuf
Real(dp), Dimension(nk,ng), Intent(in) :: f0, fz1, fz2
Real(dp), Dimension(nk,ng,6), Intent(in) :: jo, ji
Real(dp), Intent(in) :: rftem
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: fsigs
!! CXs Assigned to Nodes
Real(dp), Dimension(nk,ng,ng), Intent(out) :: sigs
Real(dp), Dimension(nk,ng), Intent(out) :: siga, sigtr, sigf, nu_sigf, sigr, D

Integer :: k, n

Do k = 1, nk
    sigs (k,1:,1:) = x_sigs (mat(k),1:,1:)
    siga (k,1:)    = x_siga (mat(k),1:)
    sigtr (k,1:)   = x_sigtr(mat(k),1:)
    sigf (k,1:)    = x_sigf (mat(k),1:)
    nu_sigf (k,1:) = x_nu_sigf (mat(k),1:)
End Do
 
	CALL ftem_updt(sigs,siga,sigtr,sigf,nu_sigf,sigr,ng,nk,nmat,&
				mat,fsigs,fsiga,fsigtr,fsigf,fnuf,xftem,rftem)
    Call D_sigr(D,sigr,siga,sigs,sigtr,ng,nk)
End Subroutine
Subroutine XS_updtrod(sigs,siga,sigtr,sigf,nu_sigf,D,sigr,ng,nk,nmat,mat,x_sigs,x_siga,x_sigtr,x_sigf,&
					 x_nu_sigf,dsiga,dsigtr,dsigf,dnuf,dsigs,pos0,ssize,f0,fz1,fz2,jo,ji,nxx,nyy,nzz,&
					 delz,xyz,nb,xbpos,fbmap) ! To Update Current XS to Base XS

Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: ng, nk, nmat, nxx, nyy, nzz, nb
!! CXs Assigned to Materials
Integer,  Dimension(nk), Intent(in) :: mat    ! Materials
Real(dp), Dimension(nzz), Intent(in) :: delz
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: x_sigs
Real(dp), Dimension(nmat,ng), Intent(in) :: x_siga, x_sigtr, x_sigf, x_nu_sigf

REAL(DP), DIMENSION(nb), INTENT(IN) :: xbpos  ! Provided control rod bank position
Real(dp), Dimension(nmat,ng), Intent(in) :: dsiga, dsigtr, dsigf, dnuf
Real(dp), Dimension(nk,ng), Intent(in) :: f0, fz1, fz2
Real(dp), Dimension(nk,ng,6), Intent(in) :: jo, ji
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: dsigs
REAL(DP), INTENT(IN) :: pos0, ssize			! Zero step position and step size
!! CXs Assigned to Nodes
Real(dp), Dimension(nk,ng,ng), Intent(out) :: sigs
Real(dp), Dimension(nk,ng), Intent(out) :: siga, sigtr, sigf, nu_sigf, sigr, D
INTEGER, DIMENSION(nxx,nyy), intent(in) :: fbmap                     ! Radial control rod bank map (node wise)

Integer :: k

Do k = 1, nk
    sigs (k,1:,1:) = x_sigs (mat(k),1:,1:)
    siga (k,1:)    = x_siga (mat(k),1:)
    sigtr (k,1:)   = x_sigtr(mat(k),1:)
    sigf (k,1:)    = x_sigf (mat(k),1:)
    nu_sigf (k,1:) = x_nu_sigf (mat(k),1:)
End Do
 
	CALL crod_updt(sigs,siga,sigtr,sigf,nu_sigf,sigr,ng,nk,nb,nmat,delz,xyz,&
				   nxx,nyy,nzz,mat,dsigs,dsiga,dsigtr,dsigf,dnuf,xbpos,pos0,ssize,&
				   f0,fz1,fz2,jo,ji,fbmap)
    Call D_sigr(D,sigr,siga,sigs,sigtr,ng,nk)

End Subroutine
Subroutine XS_updtcfmd(sigs,siga,sigtr,sigf,nu_sigf,D,sigr,ng,nk,nmat,mat,x_sigs,x_siga,x_sigtr,x_sigf,&
					 x_nu_sigf,csiga,csigtr,csigf,cnuf,fsiga,fsigtr,fsigf,fnuf,msiga,msigtr,msigf,mnuf,&
					 lsiga,lsigtr,lsigf,lnuf,dsiga,dsigtr,dsigf,dnuf,csigs,fsigs,msigs,lsigs,dsigs,xbcon,&
					 rbcon,rftem,rmtem,rcden,pos0,ssize,f0,fz1,fz2,jo,ji,nxx,nyy,nzz,&
					 delz,xyz,nb,xftem,xmtem,xcden,xbpos,fbmap) ! To Update Current XS to Base XS
					 
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: ng, nk, nmat, nxx, nyy, nzz, nb
!! CXs Assigned to Materials
Integer,  Dimension(nk), Intent(in) :: mat    ! Materials
Real(dp), Dimension(nzz), Intent(in) :: delz
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: x_sigs
Real(dp), Dimension(nmat,ng), Intent(in) :: x_siga, x_sigtr, x_sigf, x_nu_sigf
REAL(DP), DIMENSION(nk), INTENT(IN) :: xftem  ! Provided fuel temperature
REAL(DP), DIMENSION(nk), INTENT(IN) :: xmtem  ! Provided moderator temperature
REAL(DP), DIMENSION(nk), INTENT(IN) :: xcden  ! Provided coolant density
REAL(DP), DIMENSION(nb), INTENT(IN) :: xbpos  ! Provided control rod bank position
Real(dp), Dimension(nmat,ng), Intent(in) :: csiga, csigtr, csigf, cnuf,&
fsiga, fsigtr, fsigf, fnuf, msiga, msigtr, msigf, mnuf, &
lsiga, lsigtr, lsigf, lnuf, dsiga, dsigtr, dsigf, dnuf
Real(dp), Dimension(nk,ng), Intent(in) :: f0, fz1, fz2
Real(dp), Dimension(nk,ng,6), Intent(in) :: jo, ji
Real(dp), Intent(in) :: xbcon, rbcon ,rftem, rmtem, rcden
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: csigs, fsigs, msigs, lsigs, dsigs
REAL(DP), INTENT(IN) :: pos0, ssize			! Zero step position and step size
!! CXs Assigned to Nodes
Real(dp), Dimension(nk,ng,ng), Intent(out) :: sigs
Real(dp), Dimension(nk,ng), Intent(out) :: siga, sigtr, sigf, nu_sigf, sigr, D
INTEGER, DIMENSION(nxx,nyy), intent(in) :: fbmap                     ! Radial control rod bank map (node wise)

Integer :: k, n

Do k = 1, nk
    sigs (k,1:,1:) = x_sigs (mat(k),1:,1:)
    siga (k,1:)    = x_siga (mat(k),1:)
    sigtr (k,1:)   = x_sigtr(mat(k),1:)
    sigf (k,1:)    = x_sigf (mat(k),1:)
    nu_sigf (k,1:) = x_nu_sigf (mat(k),1:)
End Do
 
	CALL bcon_updt(sigs,siga,sigtr,sigf,nu_sigf,sigr,ng,nk,nmat,&
				mat,csigs,csiga,csigtr,csigf,cnuf,xbcon,rbcon)
	CALL ftem_updt(sigs,siga,sigtr,sigf,nu_sigf,sigr,ng,nk,nmat,&
				mat,fsigs,fsiga,fsigtr,fsigf,fnuf,xftem,rftem)
	CALL mtem_updt(sigs,siga,sigtr,sigf,nu_sigf,sigr,ng,nk,nmat,&
				mat,msigs,msiga,msigtr,msigf,mnuf,xmtem,rmtem)
	CALL cden_updt(sigs,siga,sigtr,sigf,nu_sigf,sigr,ng,nk,nmat,&
				mat,lsigs,lsiga,lsigtr,lsigf,lnuf,xcden,rcden)
	CALL crod_updt(sigs,siga,sigtr,sigf,nu_sigf,sigr,ng,nk,nb,nmat,delz,xyz,&
				nxx,nyy,nzz,mat,dsigs,dsiga,dsigtr,dsigf,dnuf,xbpos,pos0,ssize,&
				f0,fz1,fz2,jo,ji,fbmap)
 	!write(*,*) 'AFTER call crod_updt'
DO n = 1, nk
	if (n>1 .and. n<70) then
	!	write(*,*)  '======================'
	!	write(*,*)  'node=', n 
	!	write(*,*)  '======================'
	!	write(*,*)  'sigs=', sigs(n,:,:)
	!	write(*,*)  'siga=', siga(n,:)
	!	write(*,*)  'sigtr=', sigtr(n,:)
	!	write(*,*)  'D=', D(n,:)
	!	write(*,*)  'sigr=', sigr(n,:)
	!	write(*,*)  'sigf=', sigf(n,:)
	!	write(*,*)  'nu_sigf=', nu_sigf(n,:)
	!	write(*,*)  'xbpos=', xbpos(:)
	!	write(*,*)  'dsigs=', dsigs(:,:,:)
	!	write(*,*)  'dsiga=', dsiga(:,:)
	!	write(*,*)  'dsigtr=', dsigtr(:,:)
	!	write(*,*)  'dsigf=', dsigf(:,:)
	!	write(*,*)  'dnuf=', dnuf(:,:)
	endif
enddo

    Call D_sigr(D,sigr,siga,sigs,sigtr,ng,nk)


End Subroutine
SUBROUTINE bcon_updt (sigs,siga,sigtr,sigf,nu_sigf,sigr,ng,nk,nmat,mat,&
                      csigs,csiga,csigtr,csigf,cnuf,bcon,rbcon)!    To update CX for given boron concentration
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: ng, nk, nmat
Integer,  Dimension(nk), Intent(in) :: mat    ! Materials
!! CXs Assigned to Materials
Real(dp), Intent(in) :: bcon, rbcon
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: csigs
Real(dp), Dimension(nmat,ng), Intent(in) :: csiga, csigtr, csigf, cnuf
!! CXs Assigned to Nodes
Real(dp), Dimension(nk,ng,ng), Intent(inout) :: sigs
Real(dp), Dimension(nk,ng), Intent(inout) :: siga, sigtr, sigf, nu_sigf, sigr
Integer :: i, g, h

DO i = 1, nk
    DO g = 1, ng
        sigtr(i,g) 		= sigtr(i,g) 	+ csigtr(mat(i),g) * (bcon - rbcon)
        siga(i,g)  		= siga(i,g)  	+ csiga(mat(i),g)  * (bcon - rbcon)
        nu_sigf(i,g)   	= nu_sigf(i,g)  + cnuf(mat(i),g)   * (bcon - rbcon)
        sigf(i,g)  		= sigf(i,g)  	+ csigf(mat(i),g)  * (bcon - rbcon)
        DO h = 1, ng
            sigs(i,g,h) = sigs(i,g,h) 	+ csigs(mat(i),g,h)* (bcon - rbcon)
        END DO
    END DO
END DO


END SUBROUTINE

SUBROUTINE ftem_updt (sigs,siga,sigtr,sigf,nu_sigf,sigr,ng,nk,nmat,&
                      mat,fsigs,fsiga,fsigtr,fsigf,fnuf,ftem,rftem) !    To update CX for given fuel temp
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: ng, nk, nmat
!! CXs Assigned to Materials
Integer,  Dimension(nk), Intent(in) :: mat    ! Materials
Real(dp), Intent(in) :: rftem
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: fsigs
Real(dp), Dimension(nk), Intent(in) :: ftem
Real(dp), Dimension(nmat,ng), Intent(in) :: fsiga, fsigtr, fsigf, fnuf
!! CXs Assigned to Nodes
Real(dp), Dimension(nk,ng,ng), Intent(inout) :: sigs
Real(dp), Dimension(nk,ng), Intent(inout) :: siga, sigtr, sigf, nu_sigf, sigr
Integer :: i, g, h

DO i = 1, nk
    DO g = 1, ng
        sigtr(i,g) 		= sigtr(i,g)	+ fsigtr(mat(i),g) * (SQRT(ftem(i)) - SQRT(rftem))
        siga(i,g)  		= siga(i,g)  	+ fsiga(mat(i),g)  * (SQRT(ftem(i)) - SQRT(rftem))
        nu_sigf(i,g)   	= nu_sigf(i,g)  + fnuf(mat(i),g)   * (SQRT(ftem(i)) - SQRT(rftem))
        sigf(i,g)  		= sigf(i,g)  	+ fsigf(mat(i),g)  * (SQRT(ftem(i)) - SQRT(rftem))
        DO h = 1, ng
           sigs(i,g,h) 	= sigs(i,g,h) 	+ fsigs(mat(i),g,h)* (SQRT(ftem(i)) - SQRT(rftem))
        END DO
      END DO
END DO
END SUBROUTINE
SUBROUTINE mtem_updt (sigs,siga,sigtr,sigf,nu_sigf,sigr,ng,nk,nmat,&
                      mat,msigs,msiga,msigtr,msigf,mnuf,mtem,rmtem)!    To update CX for given moderator temperature
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: ng, nk, nmat
!! CXs Assigned to Materials
Integer,  Dimension(nk), Intent(in) :: mat    ! Materials
Real(dp), Intent(in) :: rmtem
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: msigs
Real(dp), Dimension(nk), Intent(in) :: mtem
Real(dp), Dimension(nmat,ng), Intent(in) :: msiga, msigtr, msigf, mnuf
!! CXs Assigned to Nodes
Real(dp), Dimension(nk,ng,ng), Intent(out) :: sigs
Real(dp), Dimension(nk,ng), Intent(out) :: siga, sigtr, sigf, nu_sigf, sigr
Integer :: i, g, h

DO i = 1, nk
    DO g = 1, ng
        sigtr(i,g) 		= sigtr(i,g) 	+ msigtr(mat(i),g) * (mtem(i) - rmtem)
        siga(i,g)  		= siga(i,g)  	+ msiga(mat(i),g)  * (mtem(i) - rmtem)
        nu_sigf(i,g)   	= nu_sigf(i,g)  + mnuf(mat(i),g)   * (mtem(i) - rmtem)
        sigf(i,g)  		= sigf(i,g)  	+ msigf(mat(i),g)  * (mtem(i) - rmtem)
        DO h = 1, ng
            sigs(i,g,h) = sigs(i,g,h) 	+ msigs(mat(i),g,h)* (mtem(i) - rmtem)
        END DO
    END DO
END DO


END SUBROUTINE
SUBROUTINE cden_updt (sigs,siga,sigtr,sigf,nu_sigf,sigr,ng,nk,nmat,&
                      mat,lsigs,lsiga,lsigtr,lsigf,lnuf,cden,rcden)!    To update CX for given coolant density
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: ng, nk, nmat
!! CXs Assigned to Materials
Integer,  Dimension(nk), Intent(in) :: mat    ! Materials
Real(dp), Intent(in) :: rcden
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: lsigs
Real(dp), Dimension(nk), Intent(in) :: cden
Real(dp), Dimension(nmat,ng), Intent(in) :: lsiga, lsigtr, lsigf, lnuf
!! CXs Assigned to Nodes
Real(dp), Dimension(nk,ng,ng), Intent(inout) :: sigs
Real(dp), Dimension(nk,ng), Intent(inout) :: siga, sigtr, sigf, nu_sigf, sigr
Integer :: i, g, h


DO i = 1, nk
    DO g = 1, ng
        sigtr(i,g) = sigtr(i,g) 		+ lsigtr(mat(i),g) * (cden(i) - rcden)
        siga(i,g)  = siga(i,g)  		+ lsiga(mat(i),g)  * (cden(i) - rcden)
        nu_sigf(i,g)   = nu_sigf(i,g)   + lnuf(mat(i),g)   * (cden(i) - rcden)
        sigf(i,g)  = sigf(i,g)  		+ lsigf(mat(i),g)  * (cden(i) - rcden)
        DO h = 1, ng
            sigs(i,g,h) = sigs(i,g,h)	+ lsigs(mat(i),g,h)* (cden(i) - rcden)
        END DO
    END DO
END DO


END SUBROUTINE
SUBROUTINE crod_updt (sigs,siga,sigtr,sigf,nu_sigf,sigr,ng,nk,nb,nmat,delz,xyz,&
                      nxx,nyy,nzz,mat,dsigs,dsiga,dsigtr,dsigf,dnuf,bpos,pos0,ssize,&
					  f0,fz1,fz2,jo,ji,fbmap)! TO UPDATE AND CALCUALTE VOLUME WEIGHTED HOMOGENIZED CX FOR RODDED NODE
IMPLICIT NONE
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: ng, nk, nb, nxx, nyy, nzz, nmat
Real(dp), Dimension(nk,ng), Intent(inout) :: siga, sigtr, sigf, nu_sigf, sigr
Real(dp), Dimension(nk,ng,ng), Intent(inout) :: sigs
REAL(DP), DIMENSION(nmat,ng), Intent(in) :: dsigtr, dsiga, dnuf, dsigf   ! CX incerement or decrement due to CR insertion
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: dsigs
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz
Real(dp), Dimension(nk,ng), Intent(in) :: f0, fz1, fz2
Real(dp), Dimension(nk,ng,6), Intent(in) :: jo, ji
Integer, Dimension(nk), Intent(in) :: mat             
Real(dp), Dimension(nzz), Intent(in) :: delz
REAL(DP), DIMENSION(nb), INTENT(IN) :: bpos
REAL(DP), INTENT(IN) :: pos0, ssize			! Zero step position and step size

INTEGER ::i, j, k, g, h, m
REAL(DP) :: rodh, vfrac
REAL(DP) :: dum, coreh
! Rod cusping option
INTEGER :: cusp = 0
INTEGER :: n, n1, n2, nmax
REAL(DP) :: del1, del2, eta1, eta2
REAL(DP) :: sum1, sum2, sum3, sum4, sumx
REAL(DP), DIMENSION(ng) :: sum5
REAL(DP) :: a1, a2, a3, a4, x, tx, f1, f2
REAL(DP), DIMENSION(50000) :: f
INTEGER, DIMENSION(nxx,nyy), intent(in) :: fbmap                     ! Radial control rod bank map (node wise)

! Calculate core height
coreh = 0._DP
DO m = 1, nzz
    coreh = coreh + delz(m)
END DO

DO j = 1, nyy
  DO i = 1, nxx
     IF (fbmap(i,j) > 0) THEN
        !!!(rodh -> posistion the tip of the control rod the top of core)
         rodh = coreh - pos0  - bpos(fbmap(i,j))*ssize
         dum = 0._DP
         DO k = nzz, 1, -1
           ! For partially rodded node, get volume weighted homogenized CX (0 < vfrac < 1.0)
           IF (rodh >= dum .AND. rodh < dum+delz(k)) THEN
              eta1 = rodh - dum
              eta2 = delz(k) - rodh + dum

              IF (cusp == 0 .OR. eta1 < 1. .OR. eta2 < 1) THEN    ! IF ROD CUSPING NOT ACTIVE
                vfrac = (rodh - dum) / delz(k)
                 sigtr(xyz(i,j,k),:) = sigtr(xyz(i,j,k),:) + &
                                    vfrac * dsigtr(mat(xyz(i,j,k)),:)
                 siga(xyz(i,j,k),:)  = siga(xyz(i,j,k),:) + &
                                    vfrac * dsiga(mat(xyz(i,j,k)),:)
                 nu_sigf(xyz(i,j,k),:)   = nu_sigf(xyz(i,j,k),:) + &
                                    vfrac * dnuf(mat(xyz(i,j,k)),:)
                 sigf(xyz(i,j,k),:)  = sigf(xyz(i,j,k),:) + &
                                    vfrac * dsigf(mat(xyz(i,j,k)),:)
                 sigs(xyz(i,j,k),:,:)  = sigs(xyz(i,j,k),:,:) + &
                                      vfrac * dsigs(mat(xyz(i,j,k)),:,:)
              ELSE                    ! IF ROD CUSPING ACTIVE
                 n1 = CEILING(rodh - dum)        ! Number of mesh in rodded area
                 del1 = (rodh - dum) / REAL(n1)  ! mesh size in rodded area
                 n2 = CEILING(delz(k) - rodh + dum)  ! Number of mesh in non-rodded area
                 del2 = (delz(k) - rodh + dum) / REAL(n2)  ! mesh size in non-rodded area

                 nmax = n1 + n2                     ! Total number of mesh

                 ! Calculate vectors a, b, c, d
					
                 DO g = 1, ng
                    ! Determine the flux coefficients
                    a1 = 2. * (jo(xyz(i,j,k),g,5) + ji(xyz(i,j,k),g,5) &
                       - jo(xyz(i,j,k),g,6) - ji(xyz(i,j,k),g,6))
                    a2 = 2. * (jo(xyz(i,j,k),g,5) + ji(xyz(i,j,k),g,5) &
                       + jo(xyz(i,j,k),g,6) + ji(xyz(i,j,k),g,6)) &
                       - 2. * f0(xyz(i,j,k),g)
                    a3 = 10. * a1 - 120. * fz1(xyz(i,j,k),g)
                    a4 = 35. * a2 - 700. * fz2(xyz(i,j,k),g)

                    ! Calculate fluxes in rodded area
                    x = 0.5 * delz(k)
                    tx = x / delz(k)
                    f1 = f0(xyz(i,j,k),g) + a1 * tx + a2 * (3*tx**2-0.25) &
                         + a3 * (tx*(tx+0.5)*(tx-0.5)) &
                         + a4 * ((tx**2-0.05)*(tx+0.5)*(tx-0.5))
                    DO n = 1, n1
                       x = x - del1
                       tx = x / delz(k)
                       f2 = f0(xyz(i,j,k),g) + a1 * tx + a2 * (3*tx**2-0.25) &
                            + a3 * (tx*(tx+0.5)*(tx-0.5)) &
                            + a4 * ((tx**2-0.05)*(tx+0.5)*(tx-0.5))
                      f(n) = 0.5 * (f1 + f2)
                      f1 = f2
                   END DO
                   ! Calculate fluxes in non-rodded area
                   DO n = n1+1, nmax
                      x = x - del2
                      tx = x / delz(k)
                      f2 = f0(xyz(i,j,k),g) + a1 * tx + a2 * (3*tx**2-0.25) &
                           + a3 * (tx*(tx+0.5)*(tx-0.5)) &
                           + a4 * ((tx**2-0.05)*(tx+0.5)*(tx-0.5))
                      f(n) = 0.5 * (f1 + f2)
                      f1 = f2
                   END DO

                           ! Calculate homogenized CXs
                           sumx = 0.
                           sum1 = 0.; sum2 = 0.; sum3 = 0.; sum4 = 0.; sum5 = 0.
                           DO n = 1, n1
                              sumx = sumx + f(n) * del1
                              sum1 = sum1 + f(n) * (sigtr(xyz(i,j,k),g) &
                                   + dsigtr(mat(xyz(i,j,k)),g)) * del1
                              sum2 = sum2 + f(n) * (siga(xyz(i,j,k),g) &
                                   + dsiga(mat(xyz(i,j,k)),g)) * del1
                              sum3 = sum3 + f(n) * (nu_sigf(xyz(i,j,k),g) &
                                   + dnuf(mat(xyz(i,j,k)),g)) * del1
                              sum4 = sum4 + f(n) * (sigf(xyz(i,j,k),g) &
                                   + dsigf(mat(xyz(i,j,k)),g)) * del1
                              DO h = 1, ng
                                 sum5(h) = sum5(h) + f(n) * (sigs(xyz(i,j,k),g,h) &
                                      + dsigs(mat(xyz(i,j,k)),g,h)) * del1
                              END DO
                           END DO

                           DO n = n1+1, nmax
                              sumx = sumx + f(n) * del2
                              sum1 = sum1 + f(n) * sigtr(xyz(i,j,k),g) * del2
                              sum2 = sum2 + f(n) * siga(xyz(i,j,k),g) * del2
                              sum3 = sum3 + f(n) * nu_sigf(xyz(i,j,k),g) * del2
                              sum4 = sum4 + f(n) * sigf(xyz(i,j,k),g) * del2
                              DO h = 1, ng
                                 sum5(h) = sum5(h) + f(n) &
                                         * sigs(xyz(i,j,k),g,h) * del2
                              END DO
                           END DO
						
                           sigtr(xyz(i,j,k),g) 	   = sum1 / sumx
                           siga(xyz(i,j,k),g) 	   = sum2 / sumx
                           nu_sigf(xyz(i,j,k),g)   = sum3 / sumx
                           sigf(xyz(i,j,k),g)  	   = sum4 / sumx
                           DO h = 1, ng
                              sigs(xyz(i,j,k),g,h) = sum5(h) / sumx
                           END DO

                       END DO
                    END IF

                    EXIT
                END IF


     ! For fully rodded node, vfrac = 1.
                sigtr(xyz(i,j,k),:) = sigtr(xyz(i,j,k),:) + &
                                       dsigtr(mat(xyz(i,j,k)),:)
                siga(xyz(i,j,k),:)  = siga(xyz(i,j,k),:) + &
                                       dsiga(mat(xyz(i,j,k)),:)
                nu_sigf(xyz(i,j,k),:)   = nu_sigf(xyz(i,j,k),:) + &
                                       dnuf(mat(xyz(i,j,k)),:)
                sigf(xyz(i,j,k),:)  = sigf(xyz(i,j,k),:) + &
                                       dsigf(mat(xyz(i,j,k)),:)
                sigs(xyz(i,j,k),:,:)  = sigs(xyz(i,j,k),:,:) + &
                                       dsigs(mat(xyz(i,j,k)),:,:)

                dum = dum + delz(k)
            END DO
            ! if negative CX found, Surpress CX to zero  and calculate D and sigr
            DO k = nzz, 1, -1
                DO g = 1, ng
                    IF (sigtr(xyz(i,j,k),g) < 0.) THEN
                        sigtr(xyz(i,j,k),g) = 0.
                    END IF
                    IF (siga(xyz(i,j,k),g) < 0.) THEN
                        siga(xyz(i,j,k),g) = 0.
                    END IF
                    IF (nu_sigf(xyz(i,j,k),g) < 0.) THEN
                        nu_sigf(xyz(i,j,k),g) = 0.
                    END IF
                    IF (sigf(xyz(i,j,k),g) < 0.) THEN
                        sigf(xyz(i,j,k),g) = 0.
                    END IF
                    DO h = 1, ng
                        IF (sigs(xyz(i,j,k),g,h) < 0.) THEN
                            sigs(xyz(i,j,k),g,h) = 0.
                        END IF
                    END DO
                END DO

            END DO
        END IF
    END DO
END DO

		!	write(*,*) 'crod upd'
 DO n = 1, nk
		if (n>40 .and. n<45) then
		!	write(*,*) 'sigtr=', sigtr(n,:)
		!	write(*,*) 'siga=', siga(n,:)
		!	write(*,*) 'sigf=', sigf(n,:)
		!	write(*,*) 'sigs=', sigs(n,:,:)
		!	write(*,*) 'nu_sigf=', nu_sigf(n,:)
		endif
enddo
 DO n = 1, nmat
			!	write(*,*) 'node=', n, 'dnuf=', dnuf(n,:)
enddo


END SUBROUTINE
SUBROUTINE getent(t,ent,ntem,stab,pra)!    To get enthalpy for given coolant temp. from steam table

IMPLICIT NONE
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15), Ounit = 100
INTEGER,  INTENT(IN) :: ntem, pra    ! Number of temperature in steam table
REAL(DP), INTENT(IN) :: t
REAL(DP), INTENT(OUT) :: ent

REAL(DP), DIMENSION(ntem,pra), intent(in) :: stab  ! Steam table matrix
REAL(DP) :: t1, ent1
REAL(DP) :: t2, ent2
INTEGER :: i

IF ((t < 473.15) .OR. (t > 617.91)) THEN
    WRITE(ounit,*) '  Coolant temp. : ', t
    WRITE(ounit,*) '  ERROR : MODERATOR TEMP. IS OUT OF THE RANGE OF DATA IN THE STEAM TABLE'
    WRITE(ounit,*) '  CHECK INPUT COOLANT MASS FLOW RATE OR CORE POWER'
    WRITE(*,*) '  Coolant temp. : ', t
    WRITE(*,*) '  ERROR : MODERATOR TEMP. IS OUT OF THE RANGE OF DATA IN THE STEAM TABLE'
    WRITE(*,*) '  CHECK INPUT COOLANT MASS FLOW RATE OR CORE POWER'
    STOP
END IF

t2 = stab(1,1); ent2 = stab(1,3)
DO i = 2, ntem
    t1 = t2
    ent1 = ent2
    t2 = stab(i,1); ent2 = stab(i,3)
    IF ((t >= t1) .AND. (t <= t2)) THEN
        ent = ent1 + (t - t1) / (t2 - t1) * (ent2 - ent1)
        EXIT
    END IF
END DO


END SUBROUTINE
SUBROUTINE gettd(ent,t,rho,prx,kvx,tcx,ntem,stab,pra)!    To get enthalpy for given coolant temp. from steam table

IMPLICIT NONE
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15), Ounit = 100
INTEGER,  INTENT(IN) :: ntem,pra    ! Number of temperature in steam table
REAL(DP), DIMENSION(ntem,pra), intent(in) :: stab  ! Steam table matrix
REAL(DP), INTENT(IN) :: ent
REAL(DP), INTENT(OUT) :: t, rho, prx, kvx, tcx

REAL(DP) :: t1, rho1, ent1, kv1, pr1, tc1
REAL(DP) :: t2, rho2, ent2, kv2, pr2, tc2
REAL(DP) :: ratx
INTEGER :: i

IF ((ent < 858341.5) .OR. (ent > 1624307.1)) THEN
    WRITE(ounit,*) '  Enthalpy. : ', ent
    WRITE(ounit,*) '  ERROR : ENTHALPY IS OUT OF THE RANGE OF DATA IN THE STEAM TABLE'
    WRITE(ounit,*) '  CHECK INPUT COOLANT MASS FLOW RATE OR CORE POWER'
    WRITE(*,*) '  Enthalpy. : ', ent
    WRITE(*,*) '  ERROR : ENTHALPY IS OUT OF THE RANGE OF DATA IN THE STEAM TABLE'
    WRITE(*,*) '  CHECK INPUT COOLANT MASS FLOW RATE OR CORE POWER'
    STOP
END IF


t2 = stab(1,1); rho2 = stab(1,2); ent2 = stab(1,3)
pr2 = stab(1,4); kv2 = stab(1,5); tc2 = stab(1,6)
DO i = 2, ntem
    t1 = t2
    ent1 = ent2
    rho1 = rho2
    pr1 = pr2
    kv1 = kv2
    tc1 = tc2
    t2 = stab(i,1); rho2 = stab(i,2); ent2 = stab(i,3)
    pr2 = stab(i,4); kv2 = stab(i,5); tc2 = stab(i,6)
    IF ((ent >= ent1) .AND. (ent <= ent2)) THEN
        ratx = (ent - ent1) / (ent2 - ent1)
        t   = t1   + ratx * (t2 - t1)
        rho = rho1 + ratx * (rho2 - rho1)
        prx = pr1  + ratx * (pr2 - pr1)
        kvx = kv1  + ratx * (kv2 - kv1)
        tcx = tc1 + ratx * (tc2 - tc1)
        EXIT
    END IF
END DO

END SUBROUTINE
SUBROUTINE geths(hs,xden,tc,kv,Pr,&
				 cflow,farea,dh)!    To calculate heat transfer coef.

IMPLICIT NONE
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
REAL(DP), INTENT(IN) :: xden  ! coolant densisty
REAL(DP), INTENT(IN) :: tc  ! coolant thermal conductivity
REAL(DP), INTENT(IN) :: kv  ! kinematic viscosity
REAL(DP), INTENT(IN) :: Pr  ! Prandtl Number
REAL(DP), intent(out) :: hs
REAL(DP), intent(in) :: dh, farea	! Pi diameter, Hydraulic diameter (m) and sub-channel area (m2)
REAL(DP), intent(in) :: cflow
REAL(DP) :: cvelo, Nu, Re

cvelo = cflow / (farea * xden * 1000._DP)        ! Calculate flow velocity (m/s)
Re = cvelo * dh / (kv * 1.e-6_DP)                 ! Calculate Reynolds Number
Nu = 0.023_DP*(Pr**0.4_DP)*(Re**0.8_DP)                ! Calculate Nusselt Number
hs = (tc / dh) * Nu                        ! Calculate heat transfer coefficient

END SUBROUTINE
SUBROUTINE getkf(tkf,t)!    To calculate thermal conductivity of fuel

IMPLICIT NONE
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
REAL(DP), INTENT(IN) :: t
REAL(DP), INTENT(out) :: tkf

tkf = 1.05_DP + 2150.0_DP / (t - 73.15_DP)

END SUBROUTINE
SUBROUTINE getkc(tkc,t)!    To calculate thermal conductivity of cladding
IMPLICIT NONE
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
REAL(DP), INTENT(IN) :: t
REAL(DP), INTENT(out) :: tkc

tkc = 7.51_DP + 2.09e-2_DP*t - 1.45e-5_DP*t**2 + 7.67e-9_DP*t**3

END SUBROUTINE
SUBROUTINE getcpf(cpf,t)!    To calculate specific heat capacity of fuel
IMPLICIT NONE
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
REAL(DP), INTENT(IN) :: t
REAL(DP), INTENT(out) :: cpf

cpf = 162.3_DP + 0.3038_DP*t - 2.391e-4_DP*t**2 + 6.404e-8_DP*t**3

END SUBROUTINE
SUBROUTINE getcpc(cpc,t)!    To calculate specific heat capacity of cladding
IMPLICIT NONE
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
REAL(DP), INTENT(IN) :: t
REAL(DP), INTENT(out) :: cpc

cpc = 252.54_DP + 0.11474_DP*t

END SUBROUTINE
SUBROUTINE TridiaSolve(a,b,c,d,x,nk)!    To solve tridiagonal matrix

IMPLICIT NONE
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
!REAL(DP), DIMENSION(:), INTENT(INOUT) :: a, b, c, d
!REAL(DP), DIMENSION(:), INTENT(OUT) :: x
Integer, Intent(in) :: nk
INTEGER, PARAMETER :: nm = 10      ! Fuel meat divided into 10 mesh
INTEGER, PARAMETER :: nt = nm + 2      ! two more mesh for gap and clad
REAL(DP), DIMENSION(nt+1), INTENT(INOUT) :: a, b, c, d
REAL(DP), DIMENSION(nt+1), INTENT(OUT) :: x

INTEGER :: i, n

n = SIZE(d)

! Gauss Elimination
c(1) = c(1)/b(1)
d(1) = d(1)/b(1)
DO i = 2, n
    c(i) = c(i) / (b(i) - a(i) * c(i-1))
    d(i) = (d(i) - a(i) * d(i-1)) / (b(i) - a(i) * c(i-1))
END DO

! Back Substitution
x(n) = d(n)
DO i = n-1, 1, -1
    x(i) = d(i) - c(i) * x(i+1)
END DO

END SUBROUTINE
SUBROUTINE AbsE(newF, oldF, rel, nk)!    To calculate Max Relative error

IMPLICIT NONE
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: nk
Real(dp), Dimension(nk),Intent(in) :: newF, oldF
Real(dp), Intent(out) :: rel
!Local Variables
Real(dp) :: error
Integer :: n

rel = 0.

DO n= 1, nk
    IF (ABS(newF(n)) > 1.e-10_DP) THEN
        error = ABS(newF(n) - oldF(n))
        IF (error > rel) rel = error
    END IF
END DO

END SUBROUTINE


SUBROUTINE KNE1(ng,nk,nmat,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,dsiga,dsigtr,dsigf,dnuf,dsigs,&
				pos0,ssize,nxx,nyy,nzz,delz,xyz,nb,bpos,fbmap,Keff,f0,fx1,fy1,fz1,fx2,fy2,fz2,&
				fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,order,R2,P2,R4,P4,al,jo,ji,L0,chi,ix,iy,iz,delx,dely,&
				delv,x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott)
				!    To adjuts the Keff to 1.0 if it is not equal to 1.0

IMPLICIT NONE
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: ng, nk, nb, order, nmat, nxx, nyy, nzz, x_east, x_west, y_north, y_south, z_top, z_bott
Integer, Dimension(nk), Intent(in) :: mat
Integer, Dimension(nxx), Intent(in) :: x_smax, x_smin
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin 
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk), Intent(in) :: delv, ix, iy, iz
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz
Real(dp), Dimension(nmat,ng), Intent(in) :: chi
Real(dp), Dimension(nk,ng,6), Intent(inout):: al, jo, ji
Real(dp), Dimension(nk,ng,3), Intent(inout) :: L0
Real(dp), Dimension(nk,ng,6,6),Intent(in) :: R2, P2, R4
Real(dp), Dimension(nk,ng,6,7),Intent(in) :: P4
Real(dp), Dimension(nk,ng), intent(inout) :: f0, fx1, fy1, fz1, fx2, fy2, fz2
Real(dp), Dimension(nk), intent(inout) :: fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2     ! Fission Source & Fission Source Moments
Real(dp), Intent(inout) :: Keff
!! CXs Assigned to Materials
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: x_sigs
Real(dp), Dimension(nmat,ng), Intent(inout) :: x_nu_sigf
Real(dp), Dimension(nmat,ng), Intent(in) :: x_siga, x_sigtr, x_sigf
REAL(DP), DIMENSION(nb), INTENT(IN) :: bpos
Real(dp), Dimension(nmat,ng), Intent(in) :: dsiga, dsigtr, dsigf
Real(dp), Dimension(nmat,ng), Intent(inout) :: dnuf
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: dsigs
REAL(DP), INTENT(IN) :: pos0, ssize		! Zero step position and step size
Real(dp) :: fer, fser	! Flux and Fission Source Error in BCSEARCH calcs.
!! CXs Assigned to Nodes
Real(dp), Dimension(nk,ng,ng) :: sigs
Real(dp), Dimension(nk,ng) :: siga, sigtr, sigf, nu_sigf, sigr, D
INTEGER, DIMENSION(nxx,nyy), intent(in) :: fbmap                     ! Radial control rod bank map (node wise)
Integer, Parameter :: Ounit = 100	! Output
INTEGER :: i


WRITE(ounit, *)
WRITE(ounit, '(A46,F9.6)') '  INITIAL MULTIPLICATION EFFECTIVE (K-EFF) = ', Keff
WRITE(ounit, *) '  WARNING: THE STEADY STATE K-EFF IS NOT EQUAL TO 1.0'
WRITE(ounit, *) '  AND NOW IT IS FORCED TO 1.0 BY MODIFYING THE nu*sigf CROSS SECTIONS '
WRITE(ounit, *)
DO i = 1, 10
   x_nu_sigf = x_nu_sigf / Keff
   dnuf = dnuf / Keff
	!Write(*,*) 'i=',i, 'x_nu_sigf=',x_nu_sigf(:,:)
	!Write(*,*) 'i=',i, 'dnuf=',dnuf(:,:)

    CALL XS_updtrod(sigs,siga,sigtr,sigf,nu_sigf,D,sigr,ng,nk,nmat,mat,x_sigs,x_siga,x_sigtr,x_sigf,&
					x_nu_sigf,dsiga,dsigtr,dsigf,dnuf,dsigs,pos0,ssize,f0,fz1,fz2,jo,ji,nxx,nyy,nzz,&
					delz,xyz,nb,bpos,fbmap)

	CALL Outerfs(Keff,fer,fser,f0,fx1,fy1,fz1,fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,ng,nk,&
				order,R2,P2,R4,P4,nmat,mat,al,jo,ji,L0,D,sigr,chi,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,&
				delx,dely,delz,delv,xyz,x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott)
			
	!Write(*,*) 'iteration=', 'i=',i , 'keff=', keff
	
   IF (ABS(Keff-1._DP) < 1.e-5_DP) EXIT
END DO
IF (i == 10) STOP "K-EFF STILL NOT EQUAL TO ONE. OPENNODE IS STOPPING"

END SUBROUTINE
SUBROUTINE KNEth(ng,nk,nmat,mat,x_sigs,x_siga,x_sigtr,x_sigf,x_nu_sigf,dsiga,dsigtr,dsigf,dnuf,dsigs,&
				pos0,ssize,nxx,nyy,nzz,delz,xyz,nb,bpos,fbmap,Keff,f0,fx1,fy1,fz1,fx2,fy2,fz2,&
				fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,order,R2,P2,R4,P4,al,jo,ji,L0,chi,ix,iy,iz,delx,dely,&
				delv,x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott,&
				csiga,csigtr,csigf,cnuf,fsiga,fsigtr,fsigf,fnuf,&
				msiga,msigtr,msigf,mnuf,lsiga,lsigtr,lsigf,lnuf,bcon,rbcon,rftem,rmtem,rcden,csigs,&
				fsigs,msigs,lsigs,ftem,mtem,cden)
				!    To adjuts the Keff to 1.0 if it is not equal to 1.0

IMPLICIT NONE
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: ng, nk, nb, order, nmat, nxx, nyy, nzz, x_east, x_west, y_north, y_south, z_top, z_bott
Integer, Dimension(nk), Intent(in) :: mat
Integer, Dimension(nxx), Intent(in) :: x_smax, x_smin
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin 
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk), Intent(in) :: delv, ix, iy, iz
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz
Real(dp), Dimension(nmat,ng), Intent(in) :: chi
Real(dp), Dimension(nk,ng,6), Intent(inout):: al, jo, ji
Real(dp), Dimension(nk,ng,3), Intent(inout) :: L0
Real(dp), Dimension(nk,ng,6,6),Intent(in) :: R2, P2, R4
Real(dp), Dimension(nk,ng,6,7),Intent(in) :: P4
Real(dp), Dimension(nk,ng), intent(inout) :: f0, fx1, fy1, fz1, fx2, fy2, fz2
Real(dp), Dimension(nk), intent(inout) :: fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2     ! Fission Source & Fission Source Moments
Real(dp), Intent(inout) :: Keff
!! CXs Assigned to Materials
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: x_sigs
Real(dp), Dimension(nmat,ng), Intent(inout) :: x_nu_sigf
Real(dp), Dimension(nmat,ng), Intent(in) :: x_siga, x_sigtr, x_sigf
REAL(DP), DIMENSION(nb), INTENT(IN) :: bpos
Real(dp), Dimension(nmat,ng), Intent(inout) :: dnuf
REAL(DP), INTENT(IN) :: pos0, ssize			! Zero step position and step size
Real(dp) :: fer, fser	! Flux and Fission Source Error in BCSEARCH calcs.
INTEGER, DIMENSION(nxx,nyy), intent(in) :: fbmap                     ! Radial control rod bank map (node wise)
Real(dp), Dimension(nmat,ng), Intent(in) :: csiga, csigtr, csigf, cnuf,&
fsiga,fsigtr,fsigf,fnuf,msiga,msigtr,msigf,mnuf,lsiga,lsigtr,lsigf,lnuf,&
dsiga, dsigtr, dsigf
Real(dp), Dimension(nmat,ng,ng), Intent(in) :: csigs,fsigs,msigs,lsigs,dsigs
REAL(dp), Dimension(nk), intent(in) :: ftem, mtem, cden
Real(dp), Intent(in) :: bcon, rbcon ,rftem, rmtem, rcden
!! CXs Assigned to Nodes
Real(dp), Dimension(nk,ng,ng) :: sigs
Real(dp), Dimension(nk,ng) :: siga, sigtr, sigf, nu_sigf, sigr, D

Integer, Parameter :: Ounit = 100	! Output
INTEGER :: i

WRITE(ounit, *)
WRITE(ounit, '(A46,F9.6)') '  INITIAL MULTIPLICATION EFFECTIVE (K-EFF) = ', Keff
WRITE(ounit, *) '  WARNING: THE STEADY STATE K-EFF IS NOT EQUAL TO 1.0'
WRITE(ounit, *) '  AND NOW IT IS FORCED TO 1.0 BY MODIFYING THE nu*sigf CROSS SECTIONS '
WRITE(ounit, *)
DO i = 1, 10
   x_nu_sigf = x_nu_sigf / Keff
   dnuf = dnuf / Keff
	!Write(*,*) 'i=',i, 'x_nu_sigf=',x_nu_sigf(:,:)
	!Write(*,*) 'i=',i, 'dnuf=',dnuf(:,:)

	CALL XS_updtcfmd(sigs,siga,sigtr,sigf,nu_sigf,D,sigr,ng,nk,nmat,mat,x_sigs,x_siga,x_sigtr,x_sigf,&
			 x_nu_sigf,csiga,csigtr,csigf,cnuf,fsiga,fsigtr,fsigf,fnuf,msiga,msigtr,msigf,mnuf,&
			 lsiga,lsigtr,lsigf,lnuf,dsiga,dsigtr,dsigf,dnuf,csigs,fsigs,msigs,lsigs,dsigs,bcon,&
			 rbcon,rftem,rmtem,rcden,pos0,ssize,f0,fz1,fz2,jo,ji,nxx,nyy,nzz,&
			 delz,xyz,nb,ftem,mtem,cden,bpos,fbmap)
	CALL Outerfs(Keff,fer,fser,f0,fx1,fy1,fz1,fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,ng,nk,&
				order,R2,P2,R4,P4,nmat,mat,al,jo,ji,L0,D,sigr,chi,sigs,nu_sigf,nxx,nyy,nzz,ix,iy,iz,&
				delx,dely,delz,delv,xyz,x_smax,x_smin,y_smax,y_smin,x_east,x_west,y_north,y_south,z_top,z_bott)
			
	!Write(*,*) 'iteration=', 'i=',i , 'keff=', keff
	
   IF (ABS(Keff-1._DP) < 1.e-5_DP) EXIT
END DO
IF (i == 10) STOP "K-EFF STILL NOT EQUAL TO ONE. OPENNODE IS STOPPING"

END SUBROUTINE
SUBROUTINE iPden(c0,cx1,cy1,cz1,cx2,cy2,cz2,nk,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,iBeta,lamb)
!    Calculate Initial precursor density

IMPLICIT NONE
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
INTEGER, PARAMETER :: nf = 6	! Number of delayed dneutron precusor family
Integer, Intent(in) :: nk
REAL(DP), DIMENSION(nf), intent(in) :: iBeta, lamb	! beta (delayed neutron fraction) and precusor decay constant
Real(dp), Dimension(nk), intent(in) :: fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2	! Fission Source & Fission Source Moments
REAL(DP), DIMENSION(nk,nf), intent(out) :: c0, cx1, cy1, cz1, cx2, cy2, cz2  ! neutron precusor density

INTEGER :: n, j

DO n = 1, nk
   DO j = 1, nf
      c0(n,j) = iBeta(j) * fs0(n) / lamb(j)
      cx1(n,j) = iBeta(j) * fsx1(n) / lamb(j)
      cy1(n,j) = iBeta(j) * fsy1(n) / lamb(j)
      cz1(n,j) = iBeta(j) * fsz1(n) / lamb(j)
      cx2(n,j) = iBeta(j) * fsx2(n) / lamb(j)
      cy2(n,j) = iBeta(j) * fsy2(n) / lamb(j)
      cz2(n,j) = iBeta(j) * fsz2(n) / lamb(j)
   END DO
END DO


END SUBROUTINE
SUBROUTINE PowTot (tpow,nk,ng,fx,sigf,delv)!    To calculate power distribution

IMPLICIT NONE
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: nk, ng
Real(dp), Dimension(nk,ng), Intent(in) :: sigf
Real(dp), Dimension(nk), Intent(in) :: delv
REAL(DP), DIMENSION(nk,ng), INTENT(IN) :: fx
REAL(DP), INTENT(OUT) :: tpow

REAL(DP), DIMENSION(nk) :: p
INTEGER :: g, n

p = 0.0
DO g= 1, ng
    DO n= 1, nk
        p(n) = p(n) + fx(n,g) * sigf(n,g) * delv(n)
    END DO
END DO


tpow = 0.
DO n = 1, nk
    tpow = tpow + p(n)
END DO

END SUBROUTINE
SUBROUTINE react(rho,nk,ng,nmat,nxx,nyy,nzz,mat,sigs,chi,f0,fs0,delx,dely,delz,delv,ix,iy,iz,&
				 L0,af,sigrp)!    To calculate dynamic reactivity

IMPLICIT NONE
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: nk, ng, nmat, nxx, nyy, nzz
Integer, Dimension(nk), Intent(in) :: mat
Real(dp), Dimension(nk,ng,ng), intent(in) :: sigs
Real(dp), Dimension(nk,ng), intent(in) :: f0
Real(dp), Dimension(nmat,ng), Intent(in) :: chi
Real(dp), Dimension(nk), intent(in) :: fs0
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk), Intent(in) :: delv, ix, iy, iz
Real(dp), Dimension(nk,ng,3),Intent(in) :: L0
REAL(DP), DIMENSION(nk,ng), INTENT(IN) :: af
REAL(DP), DIMENSION(nk,ng), INTENT(IN) :: sigrp
REAL(DP), INTENT(OUT) :: rho

INTEGER :: n, g, h
REAL(DP), DIMENSION(nk) :: scg
REAL(DP) :: rem, lea, src, fde


!   DO n = 1, nk
!		if (n==50) then
!			write(*,*) 'sigs=', sigs(n,:,:)
!			write(*,*) 'f0=', f0(n,:)
!			write(*,*) 'fs0=', fs0(n)
!			write(*,*) 'af=', af(n,:)
!			write(*,*) 'sigrp=', sigrp(n,:)
!		endif
!	enddo
!write(*,*) 'src=', src

src = 0.; rem = 0.; lea = 0.; fde = 0.
DO g = 1, ng
   scg = 0.
   DO h = 1, ng
      DO n = 1, nk
         IF (g /= h) scg(n) = scg(n) + sigs(n,h,g) * f0(n,h)
      END DO
   END DO
   DO n = 1, nk
      src = src + af(n,g) * (scg(n) + chi(mat(n),g) * fs0(n)) * delv(n)
      rem = rem + af(n,g) * sigrp(n,g) * f0(n,g) * delv(n)
      lea = lea + af(n,g) * &
                      (L0(n,g,1) * dely(iy(n)) * delz(iz(n)) + &
                       L0(n,g,2) * delx(ix(n)) * delz(iz(n)) + &
                       L0(n,g,3) * delx(ix(n)) * dely(iy(n)))
      fde = fde + af(n,g) * chi(mat(n),g) * fs0(n) * delv(n)
    END DO
END DO

!   DO n = 1, nk
!		if (n==50) then
!			write(*,*) 'scg=', scg(n)
!			write(*,*) 'af=', af(n,:)
!			write(*,*) 'fs0=', fs0(n)
!		endif
!	enddo
!write(*,*) 'src=', src
!
!write(*,*) 'lea=', lea
!write(*,*) 'rem=', rem
!write(*,*) 'fde=', fde
rho = (src - lea - rem) / fde
!write(*,*) 'rho=', rho

END SUBROUTINE
SUBROUTINE uPden(c0,cx1,cy1,cz1,cx2,cy2,cz2,nk,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,iBeta,lamb,h)
!    To update precursor density

IMPLICIT NONE
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
INTEGER, PARAMETER :: nf = 6	! Number of delayed dneutron precusor family
Integer, Intent(in) :: nk
REAL(DP), DIMENSION(nf), intent(in) :: iBeta, lamb	! beta (delayed neutron fraction) and precusor decay constant
Real(dp), Dimension(nk), intent(in) :: fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2	! Fission Source & Fission Source Moments
REAL(DP), DIMENSION(nk,nf), intent(out) :: c0, cx1, cy1, cz1, cx2, cy2, cz2  ! neutron precusor density

REAL(DP), INTENT(IN) :: h
REAL(DP) :: lat
INTEGER :: n, j

DO n = 1, nk
   DO j = 1, nf
      lat = (1. + lamb(j) * h)
      c0(n,j) = (c0(n,j) + iBeta(j) * h * fs0(n)) / lat
      cx1(n,j) = (cx1(n,j) + iBeta(j) * h * fsx1(n)) / lat
      cy1(n,j) = (cy1(n,j) + iBeta(j) * h * fsy1(n)) / lat
      cz1(n,j) = (cz1(n,j) + iBeta(j) * h * fsz1(n)) / lat
      cx2(n,j) = (cx2(n,j) + iBeta(j) * h * fsx2(n)) / lat
      cy2(n,j) = (cy2(n,j) + iBeta(j) * h * fsy2(n)) / lat
      cz2(n,j) = (cz2(n,j) + iBeta(j) * h * fsz2(n)) / lat
   END DO
END DO


END SUBROUTINE





Subroutine D_sigr(D,sigr,siga,sigs,sigtr,ng,nk) ! To Update Diffusion Coefficient and Removal XS
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: ng, nk
Real(dp), Dimension(nk,ng,ng), Intent(in) :: sigs
Real(dp), Dimension(nk,ng), Intent(in) :: sigtr, siga
Real(dp), Dimension(nk,ng), Intent(out):: sigr, D
! Local Variables
Integer :: k, g, h
Real(dp) :: som

Do k = 1, nk
    Do g = 1, ng
        D(k,g) = 1./(3.*sigtr(k,g))
        som = 0.
        Do h= 1, ng
            if (g/=h) som = som + sigs(k,g,h)
        End Do
        sigr(k,g) = siga(k,g) + som
    End Do
End Do

End Subroutine

Subroutine Nodal_Coup4(P4,R4,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz) ! To Calculate Response Matrices R & P
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: ng, nk, nxx, nyy, nzz
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk),  Intent(in) :: ix, iy, iz
Real(dp), Dimension(nk,ng), Intent(in) :: D, sigr
Real(dp), Dimension(nk,ng,6,6),Intent(out) :: R4
Real(dp), Dimension(nk,ng,6,7),Intent(out) :: P4
! Local Variables
Real(dp), Dimension(6,6) :: A, B
Real(dp), Dimension(6,7) :: C
Integer :: k, g
Real(dp) :: dx, dy, dz, lx, ly, lz
Real(dp) :: ax, ay, az, ax1, ay1, az1, xy, xz, yx
Real(dp) :: bx, by, bz, bx1, by1, bz1, yz, zx, zy

Do g= 1, ng
    Do k = 1, nk
        dx = D(k,g) / delx(ix(k))
        dy = D(k,g) / dely(iy(k))
        dz = D(k,g) / delz(iz(k))
        lx = 1._dp / sigr(k,g) / delx(ix(k))
        ly = 1._dp / sigr(k,g) / dely(iy(k))
        lz = 1._dp / sigr(k,g) / delz(iz(k))
    ! Matrix A
        ax = 1._dp+32.0*dx+120.0*dx*lx+960.0*dx*dx*lx+840.0*dx*dx*lx*lx
        ay = 1._dp+32.0*dy+120.0*dy*ly+960.0*dy*dy*ly+840.0*dy*dy*ly*ly
        az = 1._dp+32.0*dz+120.0*dz*lz+960.0*dz*dz*lz+840.0*dz*dz*lz*lz
         A(1,1) = ax; A(2,2) = ax
         A(3,3) = ay; A(4,4) = ay
         A(5,5) = az; A(6,6) = az
        ax1 = 8.0*dx+60.0*dx*lx+720.0*dx*dx*lx+840.0*dx*dx*lx*lx
        ay1 = 8.0*dy+60.0*dy*ly+720.0*dy*dy*ly+840.0*dy*dy*ly*ly
        az1 = 8.0*dz+60.0*dz*lz+720.0*dz*dz*lz+840.0*dz*dz*lz*lz
         A(1,2) = ax1; A(2,1) = ax1
         A(3,4) = ay1; A(4,3) = ay1
         A(5,6) = az1; A(6,5) = az1
        xy = 20.*dx*ly+840.0*dx*dx*lx*ly
        xz = 20.*dx*lz+840.0*dx*dx*lx*lz
        yx = 20.*dy*lx+840.0*dy*dy*ly*lx
        yz = 20.*dy*lz+840.0*dy*dy*ly*lz
        zx = 20.*dz*lx+840.0*dz*dz*lz*lx
        zy = 20.*dz*ly+840.0*dz*dz*lz*ly
         A(1,3) = xy; A(1,4) = xy
         A(2,3) = xy; A(2,4) = xy
         A(1,5) = xz; A(1,6) = xz
         A(2,5) = xz; A(2,6) = xz

         A(3,1) = yx; A(3,2) = yx
         A(4,1) = yx; A(4,2) = yx
         A(3,5) = yz; A(3,6) = yz
         A(4,5) = yz; A(4,6) = yz

         A(5,1) = zx; A(5,2) = zx
         A(6,1) = zx; A(6,2) = zx
         A(5,3) = zy; A(5,4) = zy
         A(6,3) = zy; A(6,4) = zy
    ! Matrix B
         B = A
        bx = 1._dp-32.0*dx+120.0*dx*lx-960.0*dx*dx*lx+840.0*dx*dx*lx*lx
        by = 1._dp-32.0*dy+120.0*dy*ly-960.0*dy*dy*ly+840.0*dy*dy*ly*ly
        bz = 1._dp-32.0*dz+120.0*dz*lz-960.0*dz*dz*lz+840.0*dz*dz*lz*lz
        bx1 = -8.0*dx+60.0*dx*lx-720.0*dx*dx*lx+840.0*dx*dx*lx*lx
        by1 = -8.0*dy+60.0*dy*ly-720.0*dy*dy*ly+840.0*dy*dy*ly*ly
        bz1 = -8.0*dz+60.0*dz*lz-720.0*dz*dz*lz+840.0*dz*dz*lz*lz
         B(1,1) = bx; B(2,2) = bx    !Replace..
         B(3,3) = by; B(4,4) = by
         B(5,5) = bz; B(6,6) = bz

         B(1,2) = bx1; B(2,1) = bx1
         B(3,4) = by1; B(4,3) = by1
         B(5,6) = bz1; B(6,5) = bz1
    ! Matrix C
         C = 0.0
        ax = 20.*dx*lx*delx(ix(k))+840.0*dx*dx*lx*lx*delx(ix(k))
        ay = 20.*dy*ly*dely(iy(k))+840.0*dy*dy*ly*ly*dely(iy(k))
        az = 20.*dz*lz*delz(iz(k))+840.0*dz*dz*lz*lz*delz(iz(k))
         C(1,1) = ax; C(2,1) = ax
         C(3,1) = ay; C(4,1) = ay
         C(5,1) = az; C(6,1) = az
        ax1 = 60.0*dx*lx*delx(ix(k))
        ay1 = 60.0*dy*ly*dely(iy(k))
        az1 = 60.0*dz*lz*delz(iz(k))
         C(1,2) =  ax1; C(2,2) = -ax1 
         C(3,3) =  ay1; C(4,3) = -ay1
         C(5,4) =  az1; C(6,4) = -az1
        ax1 = 140.0*dx*lx*delx(ix(k))
        ay1 = 140.0*dy*ly*dely(iy(k))
        az1 = 140.0*dz*lz*delz(iz(k))
         C(1,5) = ax1; C(2,5) = ax1
         C(3,6) = ay1; C(4,6) = ay1
         C(5,7) = az1; C(6,7) = az1

        Call Inverse(A,g,k,nk,ix,iy,iz)

        R4(k,g,:,:) = MATMUL(A,B)
        P4(k,g,:,:) = MATMUL(A,C)
    End Do
End Do

End Subroutine

Subroutine Nodal_Coup2(P2,R2,delx,dely,delz,D,sigr,ix,iy,iz,ng,nk,nxx,nyy,nzz) ! To Calculate Response Matrices R & P
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: ng, nk, nxx, nyy, nzz
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk),  Intent(in) :: ix, iy, iz
Real(dp), Dimension(nk,ng), Intent(in) :: D, sigr
Real(dp), Dimension(nk,ng,6,6),Intent(out) :: R2, P2
! Local Variables
Real(dp), Dimension(6,6) :: A, B, C
Integer :: k, g
Real(dp) :: dx, dy, dz, lx, ly, lz
Real(dp) :: ax, ay, az, ax1, ay1, az1, xy, xz, yx
Real(dp) :: bx, by, bz, bx1, by1, bz1, yz, zx, zy

Do g= 1, ng
    Do k = 1, nk
        dx = D(k,g) / delx(ix(k))
        dy = D(k,g) / dely(iy(k))
        dz = D(k,g) / delz(iz(k))
        lx = 1._dp / sigr(k,g) / delx(ix(k))
        ly = 1._dp / sigr(k,g) / dely(iy(k))
        lz = 1._dp / sigr(k,g) / delz(iz(k))
    ! Matrix A
        ax = 1._dp+8.0*dx+6.0*dx*lx
        ay = 1._dp+8.0*dy+6.0*dy*ly
        az = 1._dp+8.0*dz+6.0*dz*lz
         A(1,1) = ax; A(2,2) = ax
         A(3,3) = ay; A(4,4) = ay
         A(5,5) = az; A(6,6) = az
        ax1 = 4.0*dx+6.0*dx*lx
        ay1 = 4.0*dy+6.0*dy*ly
        az1 = 4.0*dz+6.0*dz*lz
         A(1,2) = ax1; A(2,1) = ax1
         A(3,4) = ay1; A(4,3) = ay1
         A(5,6) = az1; A(6,5) = az1
        xy = 6.*dx*ly
        xz = 6.*dx*lz
        yx = 6.*dy*lx
        yz = 6.*dy*lz
        zx = 6.*dz*lx
        zy = 6.*dz*ly
         A(1,3) = xy; A(1,4) = xy
         A(2,3) = xy; A(2,4) = xy
         A(1,5) = xz; A(1,6) = xz
         A(2,5) = xz; A(2,6) = xz

         A(3,1) = yx; A(3,2) = yx
         A(4,1) = yx; A(4,2) = yx
         A(3,5) = yz; A(3,6) = yz
         A(4,5) = yz; A(4,6) = yz

         A(5,1) = zx; A(5,2) = zx
         A(6,1) = zx; A(6,2) = zx
         A(5,3) = zy; A(5,4) = zy
         A(6,3) = zy; A(6,4) = zy
    ! Matrix B
         B = A
        bx = 1._dp-8.0*dx+6.0*dx*lx
        by = 1._dp-8.0*dy+6.0*dy*ly
        bz = 1._dp-8.0*dz+6.0*dz*lz
        bx1 = -4.0*dx+6.0*dx*lx
        by1 = -4.0*dy+6.0*dy*ly
        bz1 = -4.0*dz+6.0*dz*lz
         B(1,1) = bx; B(2,2) = bx    !Replace..
         B(3,3) = by; B(4,4) = by
         B(5,5) = bz; B(6,6) = bz

         B(1,2) = bx1; B(2,1) = bx1
         B(3,4) = by1; B(4,3) = by1
         B(5,6) = bz1; B(6,5) = bz1
    ! Matrix C
         C = 0.0
        ax1 = 6.0*dx*lx*delx(ix(k))
        ay1 = 6.0*dy*ly*dely(iy(k))
        az1 = 6.0*dz*lz*delz(iz(k))
         C(1,1) =  ax1; C(2,1) = ax1 
         C(3,1) =  ay1; C(4,1) = ay1
         C(5,1) =  az1; C(6,1) = az1

        Call Inverse(A,g,k,nk,ix,iy,iz)
        R2(k,g,:,:) = MATMUL(A,B)
        P2(k,g,:,:) = MATMUL(A,C)

    End Do
End Do

End Subroutine

    Subroutine Inverse(Mtx,g,n,nk,ix,iy,iz) ! To Perform Matrix Inverse by LU Decomposition 
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: nk, g, n
Real(dp), Dimension(nk), Intent(in) :: ix, iy, iz  ! x-position, y-position and z-position of the Node
Real(dp), Dimension(6,6), Intent(inout) :: Mtx
! Local Variables
Integer, Parameter :: Ounit = 100   !output
Real(dp), Dimension(6,6) :: L, U, imat, pmat
Real(dp), Dimension(6) :: y
Real(dp) :: piv, isum
Integer :: i, j, k

pmat = Mtx
U = Mtx
L = 0.0
! Start Matrix Decomposition
Do i= 1, 6
    if (ABS(Mtx(i,i)) < 10e-3) Then
      Write(ounit,*) 'Error in Matrix Decomposition: Diagonal Elements Close to Zero'
      Write(ounit,2001) g, ix(n), iy(n), iz(n)
      Write(*,*) 'Error in Matrix Decomposition: Diagonal Elements Close to Zero'
      Write(*,2001) g, ix(n), iy(n), iz(n)
      Stop
    End if
    L(i,i) = 1._dp
    Do j= i+1, 6
        piv = U(j,i)/U(i,i)
        L(j,i) = piv
        Do k= i, 6
            U(j,k) = U(j,k) - piv*U(i,k)
        End Do
        U(j,i) = 0.0
    End Do
End Do
! Check matrix decomposition
Do i = 1,6
    Do j = 1,6
        isum = 0.0
        Do k = 1,6
            isum = isum+L(i,k)*U(k,j)
        End Do
        if (ABS(Mtx(i,j)-isum)/ABS(Mtx(i,j)) > 1.e-3_dp) Then
            Write(ounit,*) 'Error in Matrix Decomposition: Decomposition Failed'
            Write(ounit,2001) g, ix(n), iy(n), iz(n)
            Write(*,*) 'Error in Matrix Decomposition: Decomposition Failed'
            Write(*,2001) g, ix(n), iy(n), iz(n)
            Stop
        End if
    End Do
End Do
!Initialiaze Identity matrix
imat = 0.0
Do i= 1, 6
    imat(i,i) = 1._dp
End Do
! Calculate matrix inverse
! Ref: https://www.gamedev.net/resources/_/technical/math-and-physics/matrix-inversion-using-lu-decomposition-r3637
Do j=1,6   ! For each column
    !Solve y in Ly = b (Forward substitution)
    y(1) = imat(1,j)
    Do i=2,6
        isum = 0.0
        Do k =1, i-1
            isum = isum + L(i,k)*y(k)
        End Do
        y(i) = imat(i,j)-isum
    End Do
    ! Solve x in Ux=y(Backward substitution) and store inverse matrix to input matrix 'mat'
    Mtx(6,j) = y(6)/U(6,6)
    Do i = 5,1,-1
        isum = 0.0
        Do k =i+1,6
            isum = isum + U(i,k)*Mtx(k,j)
        End Do
        Mtx(i,j) = (y(i)-isum) / U(i,i)
    End Do
End Do
!Check Matrix Inverse
Do i = 1,6
    Do j = 1,6
        isum = 0.0
        Do k = 1,6
            isum = isum+pmat(i,k)*Mtx(k,j)
        End Do
        if (ABS(imat(i,j)-isum) > 1.e-4_dp) Then
            Write(ounit,*) 'Error in Matrix Inversion'
            Write(ounit,2001) g, ix(n), iy(n), iz(n)
            Write(*,*) 'Error in Matrix Inversion'
            Write(*,2001) g, ix(n), iy(n), iz(n)
            Stop
        End if
    End Do
End Do

2001 Format(2X, 'Group = ', I2, ', I = ', I2, ', J = ', I2, ', K = ', I2)

End Subroutine

Subroutine FSrc(fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,&
                f0,fx1,fy1,fz1,fx2,fy2,fz2,nu_sigf,ng,nk) ! To Calculate Fission Source & Fission Source Moments
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: ng, nk
Real(dp), Dimension(nk,ng), Intent(in) :: nu_sigf, f0, fx1, fy1, fz1, fx2, fy2, fz2
Real(dp), Dimension(nk), Intent(out) :: fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2
! Local Variables
Integer :: g, k

fs0 = 0.; fsx1 = 0.; fsy1 = 0.; fsz1 = 0.; fsx2 = 0.; fsy2 = 0.; fsz2 = 0.
Do g = 1, ng
    Do k = 1, nk
       fs0(k)  = fs0(k)  + nu_sigf(k,g) * f0 (k,g)
       fsx1(k) = fsx1(k) + nu_sigf(k,g) * fx1(k,g) 
       fsy1(k) = fsy1(k) + nu_sigf(k,g) * fy1(k,g) 
       fsz1(k) = fsz1(k) + nu_sigf(k,g) * fz1(k,g)
       fsx2(k) = fsx2(k) + nu_sigf(k,g) * fx2(k,g) 
       fsy2(k) = fsy2(k) + nu_sigf(k,g) * fy2(k,g) 
       fsz2(k) = fsz2(k) + nu_sigf(k,g) * fz2(k,g)
    End Do 
End Do

End Subroutine

Subroutine FSrcAdj(fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,mat,chi,&
                   f0,fx1,fy1,fz1,fx2,fy2,fz2,ng,nk,nmat)    ! To Calculate Adjoint Fission Source & Fission Source Moments
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: ng, nk, nmat
Integer, Dimension(nk), Intent(in) :: mat
Real(dp), Dimension(nmat,ng), Intent(in) :: chi
Real(dp), Dimension(nk,ng), Intent(in) :: f0, fx1, fy1, fz1, fx2, fy2, fz2
Real(dp), Dimension(nk), Intent(out) :: fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2
! Local Variables
Integer :: g, k

fs0 = 0.; fsx1 = 0.; fsy1 = 0.; fsz1 = 0.; fsx2 = 0.; fsy2 = 0.; fsz2 = 0.
Do g = 1, ng
    Do k = 1, nk
      fs0(k)  = fs0(k)  + f0 (k,g) * chi(mat(k),g)
      fsx1(k) = fsx1(k) + fx1(k,g) * chi(mat(k),g)
      fsy1(k) = fsy1(k) + fy1(k,g) * chi(mat(k),g)
      fsz1(k) = fsz1(k) + fz1(k,g) * chi(mat(k),g)
      fsx2(k) = fsx2(k) + fx2(k,g) * chi(mat(k),g)
      fsy2(k) = fsy2(k) + fy2(k,g) * chi(mat(k),g)
      fsz2(k) = fsz2(k) + fz2(k,g) * chi(mat(k),g)
    End Do
End Do

End Subroutine 

Subroutine Integrate(ss,s,delv,nk) ! To Perform Fs Volume Integration
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: nk
Real(dp), Dimension(nk), Intent(in) :: delv, s
Real(dp), Intent(out) :: ss
! Local Variables
Integer:: k

ss = 0.
Do k = 1, nk
    ss = ss + delv(k) * s(k)
End Do

End Subroutine

Subroutine TSrc(Q,g,Keff,ng,nmat,nk,chi,mat,sigs,f0,fx1,fy1,fz1,&
                        fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2) ! To Update Total Source
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: g, ng, nk, nmat
Integer,  Dimension(nk), Intent(in) :: mat
Real(dp), Intent(in) :: Keff
Real(dp), Dimension(nmat,ng), Intent(in) :: chi
Real(dp), Dimension(nk,ng,ng), Intent(in) :: sigs
Real(dp), Dimension(nk,ng), Intent(in) :: f0, fx1, fy1, fz1, fx2, fy2, fz2
Real(dp), Dimension(nk), Intent(in) :: fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2
Real(dp), Dimension(nk,ng,7),Intent(out) :: Q
! Local Variables
Real(dp), Dimension(nk) :: s0, sx1, sy1, sz1, sx2, sy2, sz2     ! Scattering Source and Scattering Source Moments
Integer :: h, k

s0 = 0.; sx1 = 0.; sy1 = 0.; sz1 = 0.; sx2 = 0.; sy2 = 0.; sz2 = 0.
Do h = 1, ng
    Do k = 1, nk
        if (g /= h) Then
            s0(k)  = s0(k)  + sigs(k,h,g) * f0(k,h)
            sx1(k) = sx1(k) + sigs(k,h,g) * fx1(k,h)
            sy1(k) = sy1(k) + sigs(k,h,g) * fy1(k,h)
            sz1(k) = sz1(k) + sigs(k,h,g) * fz1(k,h)
            sx2(k) = sx2(k) + sigs(k,h,g) * fx2(k,h)
            sy2(k) = sy2(k) + sigs(k,h,g) * fy2(k,h)
            sz2(k) = sz2(k) + sigs(k,h,g) * fz2(k,h)
        End if
    End Do
End Do
Do k = 1, nk
    Q(k,g,1) = chi(mat(k),g) * fs0(k)/Keff  + s0(k)
    Q(k,g,2) = chi(mat(k),g) * fsx1(k)/Keff + sx1(k)
    Q(k,g,3) = chi(mat(k),g) * fsy1(k)/Keff + sy1(k)
    Q(k,g,4) = chi(mat(k),g) * fsz1(k)/Keff + sz1(k)
    Q(k,g,5) = chi(mat(k),g) * fsx2(k)/Keff + sx2(k)
    Q(k,g,6) = chi(mat(k),g) * fsy2(k)/Keff + sy2(k)
    Q(k,g,7) = chi(mat(k),g) * fsz2(k)/Keff + sz2(k)
End Do

End Subroutine

Subroutine TSrcAdj(Q,Keff,sigs,nu_sigf,f0,fx1,fy1,fz1,fx2,fy2,fz2,&
                   fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,g,ng,nk)      ! To Update Adjoint Total Source
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: g, ng, nk
Real(dp), Intent(in) :: Keff
Real(dp), Dimension(nk,ng,ng), Intent(in) :: sigs
Real(dp), Dimension(nk,ng), Intent(in) :: nu_sigf, f0, fx1, fy1, fz1, fx2, fy2, fz2
Real(dp), Dimension(nk), Intent(in) :: fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2
Real(dp), Dimension(nk,ng,7),Intent(out) :: Q
! Local Variables
Real(dp), Dimension(nk) :: s0, sx1, sy1, sz1, sx2, sy2, sz2     ! Scattering Source and Scattering Source Moments
Integer :: h, k

s0 = 0.; sx1 = 0.; sy1 = 0.; sz1 = 0.; sx2 = 0.; sy2 = 0.; sz2 = 0.
Do h = 1, ng
    Do k = 1, nk
        if (g /= h) Then
            s0(k)  = s0(k)  + sigs(k,g,h) * f0(k,h)
            sx1(k) = sx1(k) + sigs(k,g,h) * fx1(k,h)
            sy1(k) = sy1(k) + sigs(k,g,h) * fy1(k,h)
            sz1(k) = sz1(k) + sigs(k,g,h) * fz1(k,h)
            sx2(k) = sx2(k) + sigs(k,g,h) * fx2(k,h)
            sy2(k) = sy2(k) + sigs(k,g,h) * fy2(k,h)
            sz2(k) = sz2(k) + sigs(k,g,h) * fz2(k,h)
        End if
    End Do
End Do
Do k = 1, nk
    Q(k,g,1) = nu_sigf(k,g) * fs0(k)/Keff  + s0(k)
    Q(k,g,2) = nu_sigf(k,g) * fsx1(k)/Keff + sx1(k)
    Q(k,g,3) = nu_sigf(k,g) * fsy1(k)/Keff + sy1(k)
    Q(k,g,4) = nu_sigf(k,g) * fsz1(k)/Keff + sz1(k)
    Q(k,g,5) = nu_sigf(k,g) * fsx2(k)/Keff + sx2(k)
    Q(k,g,6) = nu_sigf(k,g) * fsy2(k)/Keff + sy2(k)
    Q(k,g,7) = nu_sigf(k,g) * fsz2(k)/Keff + sz2(k)
End Do

End Subroutine

Subroutine TSrcFx(Q,g,ng,nmat,nk,chi,mat,sigs,f0,fx1,fy1,fz1,&
                  fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,eSrc) ! To Update Total Source
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: g, ng, nk, nmat
Integer,  Dimension(nk), Intent(in) :: mat
Real(dp), Dimension(nmat,ng), Intent(in) :: chi
Real(dp), Dimension(nk,ng), Intent(in) :: eSrc
Real(dp), Dimension(nk,ng,ng), Intent(in) :: sigs
Real(dp), Dimension(nk,ng), Intent(in) :: f0, fx1, fy1, fz1, fx2, fy2, fz2
Real(dp), Dimension(nk), Intent(in) :: fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2
Real(dp), Dimension(nk,ng,7),Intent(out) :: Q
! Local Variables
Real(dp), Dimension(nk) :: s0, sx1, sy1, sz1, sx2, sy2, sz2     ! Scattering Source and Scattering Source Moments
Integer :: h, k

s0 = 0.; sx1 = 0.; sy1 = 0.; sz1 = 0.; sx2 = 0.; sy2 = 0.; sz2 = 0.
Do h = 1, ng
    Do k = 1, nk
        if (g /= h) Then
            s0(k)  = s0(k)  + sigs(k,h,g) * f0(k,h)
            sx1(k) = sx1(k) + sigs(k,h,g) * fx1(k,h)
            sy1(k) = sy1(k) + sigs(k,h,g) * fy1(k,h)
            sz1(k) = sz1(k) + sigs(k,h,g) * fz1(k,h)
            sx2(k) = sx2(k) + sigs(k,h,g) * fx2(k,h)
            sy2(k) = sy2(k) + sigs(k,h,g) * fy2(k,h)
            sz2(k) = sz2(k) + sigs(k,h,g) * fz2(k,h)
        End if
    End Do
End Do
Do k = 1, nk
    Q(k,g,1) = chi(mat(k),g) * fs0(k)  + s0(k) + eSrc(k,g)
    Q(k,g,2) = chi(mat(k),g) * fsx1(k) + sx1(k)
    Q(k,g,3) = chi(mat(k),g) * fsy1(k) + sy1(k)
    Q(k,g,4) = chi(mat(k),g) * fsz1(k) + sz1(k)
    Q(k,g,5) = chi(mat(k),g) * fsx2(k) + sx2(k)
    Q(k,g,6) = chi(mat(k),g) * fsy2(k) + sy2(k)
    Q(k,g,7) = chi(mat(k),g) * fsz2(k) + sz2(k)
End Do

End Subroutine

SUBROUTINE TSrcT(Q,g,ng,nmat,nk,chi,mat,sigs,f0,fx1,fy1,fz1,&
				 fx2,fy2,fz2,fs0,fsx1,fsy1,fsz1,fsx2,fsy2,fsz2,&
				 iBeta,lamb,tbeta,c0,cx1,cy1,cz1,cx2,cy2,cz2,&
				 ft,ftx1,fty1,ftz1,ftx2,fty2,ftz2,omeg,velo,ht)!   To update total source for transient calcs. with exponetial transformation

Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: g, ng, nk, nmat
INTEGER, PARAMETER :: nf = 6                       ! Number of delaye dneutron precusor family
Integer,  Dimension(nk), Intent(in) :: mat
Real(dp), Dimension(nmat,ng), Intent(in) :: chi
Real(dp), Dimension(nk,ng,ng), Intent(in) :: sigs
Real(dp), Dimension(nk,ng), Intent(in) :: f0, fx1, fy1, fz1, fx2, fy2, fz2
Real(dp), Dimension(nk), Intent(in) :: fs0, fsx1, fsy1, fsz1, fsx2, fsy2, fsz2
Real(dp), Dimension(nk,ng,7),Intent(out) :: Q
REAL(DP), DIMENSION(nf), intent(in) :: iBeta, lamb	! beta (delayed neutron fraction) and precusor decay constant
REAL(DP), intent(in) :: tbeta	! total beta
REAL(DP), DIMENSION(nk,nf), intent(in) :: c0, cx1, cy1, cz1, cx2, cy2, cz2  ! neutron precusor density
REAL(DP), DIMENSION(nk,ng), intent(in) :: ft, ftx1, fty1, ftz1, ftx2, fty2, ftz2  ! Parameters at previous time step
REAL(DP), DIMENSION(nk,ng), intent(in) :: omeg	! Exponential transformation constant
REAL(DP), DIMENSION(ng), intent(in) :: velo            ! Neutron velocity
REAL(DP), INTENT(IN) :: ht

! Local Variables
Real(dp), Dimension(nk) :: s0, sx1, sy1, sz1, sx2, sy2, sz2     ! Scattering Source and Scattering Source Moments
REAL(DP) :: dt, dtx1, dty1, dtz1, dtx2, dty2, dtz2, lat, dfis
INTEGER :: n, i, h

s0 = 0.; sx1 = 0.; sy1 = 0.; sz1 = 0.
sx2 = 0.; sy2 = 0.; sz2 = 0.

DO h = 1, ng
    DO n = 1, nk
        IF (g /= h) THEN
            s0(n)   = s0(n) + sigs(n,h,g) * f0(n,h)
            sx1(n) = sx1(n) + sigs(n,h,g) * fx1(n,h)
            sy1(n) = sy1(n) + sigs(n,h,g) * fy1(n,h)
            sz1(n) = sz1(n) + sigs(n,h,g) * fz1(n,h)
            sx2(n) = sx2(n) + sigs(n,h,g) * fx2(n,h)
            sy2(n) = sy2(n) + sigs(n,h,g) * fy2(n,h)
            sz2(n) = sz2(n) + sigs(n,h,g) * fz2(n,h)
        END IF
    END DO
END DO

DO n = 1, nk
     dt = 0.; dtx1 = 0.; dty1 = 0.; dtz1 = 0.; dtx2 = 0.; dty2 = 0.; dtz2 = 0.
     dfis = 0.
     DO i = 1, nf
        lat = 1. + lamb(i) * ht
        dt = dt  + lamb(i) * c0(n,i) / lat
        dtx1 = dtx1 + lamb(i) * cx1(n,i) / lat
        dty1 = dty1 + lamb(i) * cy1(n,i) / lat
        dtz1 = dtz1 + lamb(i) * cz1(n,i) / lat
        dtx2 = dtx2 + lamb(i) * cx2(n,i) / lat
        dty2 = dty2 + lamb(i) * cy2(n,i) / lat
        dtz2 = dtz2 + lamb(i) * cz2(n,i) / lat
        dfis = dfis + chi(mat(n),g) * iBeta(i) * lamb(i) * ht / lat
    END DO
    Q(n,g,1) = ((1. - tbeta) * chi(mat(n),g) + dfis) * fs0(n)  &
    + s0(n) + chi(mat(n),g) * dt + ft(n,g)  * EXP(omeg(n,g) * ht) / (velo(g) * ht)
    Q(n,g,2) = ((1. - tbeta) * chi(mat(n),g) + dfis) * fsx1(n)  &
    + sx1(n) + chi(mat(n),g) * dtx1 + ftx1(n,g)  * EXP(omeg(n,g) * ht) / (velo(g) * ht)
    Q(n,g,3) = ((1. - tbeta) * chi(mat(n),g) + dfis) * fsy1(n)  &
    + sy1(n) + chi(mat(n),g) * dty1 + fty1(n,g)  * EXP(omeg(n,g) * ht) / (velo(g) * ht)
    Q(n,g,4) = ((1. - tbeta) * chi(mat(n),g) + dfis) * fsz1(n)  &
    + sz1(n) + chi(mat(n),g) * dtz1 + ftz1(n,g)  * EXP(omeg(n,g) * ht) / (velo(g) * ht)
    Q(n,g,5) = ((1. - tbeta) * chi(mat(n),g) + dfis) * fsx2(n)  &
    + sx2(n) + chi(mat(n),g) * dtx2 + ftx2(n,g)  * EXP(omeg(n,g) * ht) / (velo(g) * ht)
    Q(n,g,6) = ((1. - tbeta) * chi(mat(n),g) + dfis) * fsy2(n)  &
    + sy2(n) + chi(mat(n),g) * dty2 + fty2(n,g)  * EXP(omeg(n,g) * ht) / (velo(g) * ht)
    Q(n,g,7) = ((1. - tbeta) * chi(mat(n),g) + dfis) * fsz2(n)  &
    + sz2(n) + chi(mat(n),g) * dtz2 + ftz2(n,g)  * EXP(omeg(n,g) * ht) / (velo(g) * ht)
END DO

END SUBROUTINE

Subroutine Jin(ji,g,k,jo,al,ng,nk,nxx,nyy,nzz,ix,iy,iz,x_smax,x_smin,&
                           y_smax,y_smin,xyz,x_east,x_west,y_north,y_south,z_top,z_bott) ! To Calculate Ingoing Partial Currents from Neighborhod Nodes
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: g, k, ng, nk, nxx, nyy, nzz, x_east, x_west, y_north, y_south, z_top, z_bott
Integer, Dimension(nxx), Intent(in) :: x_smax, x_smin
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin
Real(dp), Dimension(nk), Intent(in) :: ix, iy, iz
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz
Real(dp), Dimension(nk,ng,6), Intent(in) :: al, jo
Real(dp), Dimension(nk,ng,6), Intent(out):: ji

! g : Group, k : Nod, (Nmbr): Direction of Axes (+) or (-)         
    if (ix(k) == y_smax(iy(k))) Then    ! East (X+) BC
        Call BCond(ji,x_east,k,g,1,jo,ng,nk)
    Else
        ji(k,g,1) = (jo(xyz(ix(k)+1,iy(k),iz(k)),g,2) + &  
                    al(k,g,1)*jo(k,g,1)) / (1._dp-al(k,g,1))
    End if
    if (ix(k) == y_smin(iy(k))) Then    ! West (X-) BC
        Call BCond(ji,x_west,k,g,2,jo,ng,nk)
    Else
        ji(k,g,2) = (jo(xyz(ix(k)-1,iy(k),iz(k)),g,1) + &  
                    al(k,g,2)*jo(k,g,2)) / (1._dp-al(k,g,2))
    End if
    if (iy(k) == x_smax(ix(k))) Then    ! North (Y+) BC
        Call BCond(ji,y_north,k,g,3,jo,ng,nk)
    Else
        ji(k,g,3) = (jo(xyz(ix(k),iy(k)+1,iz(k)),g,4) + &  
                    al(k,g,3)*jo(k,g,3)) / (1._dp-al(k,g,3))
    End if
    if (iy(k) == x_smin(ix(k))) Then    ! South (Y-) BC
        Call BCond(ji,y_south,k,g,4,jo,ng,nk)
    Else
        ji(k,g,4) = (jo(xyz(ix(k),iy(k)-1,iz(k)),g,3) + &  
                    al(k,g,4)*jo(k,g,4)) / (1._dp-al(k,g,4))
    End if
    If (iz(k) == nzz) Then              ! Top (Z+) BC
        Call BCond(ji,z_top,k,g,5,jo,ng,nk)
    Else
        ji(k,g,5) = (jo(xyz(ix(k),iy(k),iz(k)+1),g,6) + &
                    al(k,g,5)*jo(k,g,5)) / (1._dp-al(k,g,5))
    End if
    If (iz(k) == 1) Then                ! Bottom (Z-)BC
        Call BCond(ji,z_bott,k,g,6,jo,ng,nk)
    Else
        ji(k,g,6) = (jo(xyz(ix(k),iy(k),iz(k)-1),g,5) + &
                    al(k,g,6)*jo(k,g,6)) / (1._dp-al(k,g,6))
    End if

End Subroutine

    Subroutine BCond(ji,bc,k,g,side,jo,ng,nk) ! To Provide proper Boundary Conditions
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: g, k, ng, nk, bc, side
Real(dp), Dimension(nk,ng,6),Intent(in) :: jo
Real(dp), Dimension(nk,ng,6),Intent(out):: ji

    if (bc == 0) Then
        ji(k,g,side) = -jo(k,g,side)
    Else if (bc == 1) Then                    ! Vacuum Boundary
        ji(k,g,side) = 0.0
    Else                                      ! Reflective Boundary
        ji(k,g,side) = jo(k,g,side)
    End if

End Subroutine

Subroutine QTL(L,g,k,ng,nk,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,xyz,x_smax,x_smin,&
                           y_smax,y_smin,L0,x_east,x_west,y_north,y_south,z_top,z_bott) ! To Calculate Transverse Leakage Moments using Quadratic Transverse Leakage Approximation
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: g, k, ng, nk, nxx, nyy, nzz, x_east, x_west, y_north, y_south, z_top, z_bott
Integer, Dimension(nxx), Intent(in) :: x_smax, x_smin
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin
Real(dp), Dimension(nk), Intent(in) :: ix, iy, iz
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz
Real(dp), Dimension(nk,ng,3), Intent(in) :: L0
Real(dp), Dimension(7), Intent(out) :: L     ! Transverse Leakage Moments(0, Lx1, Ly1, Lz1, Lx2, Ly2, Lz2)
!Local Variables
Real(dp) :: xm, xp, ym, yp, zm, zp
Real(dp) :: p1m, p2m, p1p, p2p, mm, pp, p1, p2, p3
Real(dp) :: p1xy, p2xy, p1xz, p2xz
Real(dp) :: p1yx, p2yx, p1yz, p2yz
Real(dp) :: p1zx, p2zx, p1zy, p2zy
         
! Set Paramaters for X_Direction Transverse Leakage
if (ix(k) == y_smax(iy(k))) Then                ! X+ Plan
    if (x_east == 0 .OR. x_east == 1) Then         ! Vacuum Boundary Conditions
        xm = delx(ix(k)-1)/delx(ix(k))
        p1m = xm+1._dp
        
        p1xy = 2.*( L0(xyz(ix(k),iy(k),iz(k)),g,2)   &
                  - L0(xyz(ix(k)-1,iy(k),iz(k)),g,2) &
                  ) / p1m
        p2xy = 0.0
        p1xz = 2.*( L0(xyz(ix(k),iy(k),iz(k)),g,3)   &
                  - L0(xyz(ix(k)-1,iy(k),iz(k)),g,3) &
                  ) / p1m
        p2xz = 0.0
    Else                                           ! Reflective Boundary Conditions
        xm = delx(ix(k)-1)/delx(ix(k))
        xp = 1._dp
        p1m = xm+1._dp; p2m = 2.*xm+1._dp; mm = p1m*p2m; p2 = xm+xp+2.;
        p1p = xp+1._dp; p2p = 2.*xp+1._dp; pp = p1p*p2p; p1 = pp-mm;
        p3 = p1m*p1p*(xm+xp+1._dp)
               
        P1xy = ( mm * L0(xyz(ix(k),iy(k),iz(k)),g,2)      &
               - pp * L0(xyz(ix(k)-1,iy(k),iz(k)),g,2)    &
               + p1 * L0(k,g,2)                           &
               ) / p3
        P2xy = ( p1m * L0(xyz(ix(k),iy(k),iz(k)),g,2)     &
               + p1p * L0(xyz(ix(k)-1,iy(k),iz(k)),g,2)   &
               - p2  * L0(xyz(ix(k),iy(k),iz(k)),g,2)     &
               ) / p3
        P1xz = ( mm * L0(xyz(ix(k),iy(k),iz(k)),g,3)      &
               - pp * L0(xyz(ix(k)-1,iy(k),iz(k)),g,3)    &
               + p1 * L0(k,g,3)                           &
               ) / p3
        P2xz = ( p1m * L0(xyz(ix(k),iy(k),iz(k)),g,3)     &
               + p1p * L0(xyz(ix(k)-1,iy(k),iz(k)),g,3)   &
               - p2  * L0(xyz(ix(k),iy(k),iz(k)),g,3)     &
               ) / p3
    End if
Else if (ix(k) == y_smin(iy(k))) Then           ! X- Plan                       
    if (x_west == 0 .OR. x_west == 1) Then         ! Vacuum Boundary Conditions
        xp = delx(ix(k)+1)/delx(ix(k))
        p1p = xp+1._dp
        
        p1xy = 2.*( L0(xyz(ix(k)+1,iy(k),iz(k)),g,2) &
                  - L0(xyz(ix(k),iy(k),iz(k)),g,2)   &
                  ) / p1p
        p2xy = 0.0
        p1xz = 2.*( L0(xyz(ix(k)+1,iy(k),iz(k)),g,3) &
                  - L0(xyz(ix(k),iy(k),iz(k)),g,3)   &
                  ) / p1p
        p2xz = 0.0
    Else                                           ! Reflective Boundary Conditions
        xm = 1._dp
        xp = delx(ix(k)+1)/delx(ix(k))
        p1m = xm+1._dp; p2m = 2.*xm+1._dp; mm = p1m*p2m; p2 = xm+xp+2.;
        p1p = xp+1._dp; p2p = 2.*xp+1._dp; pp = p1p*p2p; p1 = pp-mm;
        p3 = p1m*p1p*(xm+xp+1._dp)
         
        P1xy = ( mm * L0(xyz(ix(k)+1,iy(k),iz(k)),g,2)   &
               - pp * L0(xyz(ix(k),iy(k),iz(k)),g,2)     &
               + p1 * L0(k,g,2)                          &
               ) / p3
        P2xy = ( p1m * L0(xyz(ix(k)+1,iy(k),iz(k)),g,2)  &
               + p1p * L0(xyz(ix(k),iy(k),iz(k)),g,2)    &
               - p2  * L0(xyz(ix(k),iy(k),iz(k)),g,2)    &
               ) / p3
        P1xz = ( mm * L0(xyz(ix(k)+1,iy(k),iz(k)),g,3)   &
               - pp * L0(xyz(ix(k),iy(k),iz(k)),g,3)     &
               + p1 * L0(k,g,3)                          &
               ) / p3
        P2xz = ( p1m * L0(xyz(ix(k)+1,iy(k),iz(k)),g,3)  &
               + p1p * L0(xyz(ix(k),iy(k),iz(k)),g,3)    &
               - p2  * L0(xyz(ix(k),iy(k),iz(k)),g,3)    &
               ) / p3
    End if
Else
    xm = delx(ix(k)-1)/delx(ix(k))
    xp = delx(ix(k)+1)/delx(ix(k))     
    p1m = xm+1._dp; p2m = 2.*xm+1._dp; mm = p1m*p2m; p2 = xm+xp+2.    
    p1p = xp+1._dp; p2p = 2.*xp+1._dp; pp = p1p*p2p; p1 = pp-mm;
    p3 = p1m*p1p*(xm+xp+1._dp)
    
    P1xy = ( mm * L0(xyz(ix(k)+1,iy(k),iz(k)),g,2)   &
           - pp * L0(xyz(ix(k)-1,iy(k),iz(k)),g,2)   &
           + p1 * L0(k,g,2)                          &
           ) / p3
    P2xy = ( p1m * L0(xyz(ix(k)+1,iy(k),iz(k)),g,2)  &
           + p1p * L0(xyz(ix(k)-1,iy(k),iz(k)),g,2)  &
           - p2  * L0(xyz(ix(k),iy(k),iz(k)),g,2)    &
           ) / p3
    P1xz = ( mm * L0(xyz(ix(k)+1,iy(k),iz(k)),g,3)   &
           - pp * L0(xyz(ix(k)-1,iy(k),iz(k)),g,3)   &
           + p1 * L0(k,g,3)                          &
           ) / p3
    P2xz = ( p1m * L0(xyz(ix(k)+1,iy(k),iz(k)),g,3)  &
           + p1p * L0(xyz(ix(k)-1,iy(k),iz(k)),g,3)  &
           - p2  * L0(xyz(ix(k),iy(k),iz(k)),g,3)    &
           ) / p3
End if

! Set Paramaters for Y_Direction Transverse Leakage
if (iy(k) == x_smax(ix(k))) Then                ! Y+ Plan
    if (y_north == 0 .OR. y_north == 1) Then       ! Vacuum Boundary Conditions
        ym = dely(iy(k)-1)/dely(iy(k))
        p1m = ym+1._dp
        
        p1yx = 2.*( L0(xyz(ix(k),iy(k),iz(k)),g,1)   &
                  - L0(xyz(ix(k),iy(k)-1,iz(k)),g,1) &
                  ) / p1m
        p2yx = 0.0
        p1yz = 2.*( L0(xyz(ix(k),iy(k),iz(k)),g,3)   &
                  - L0(xyz(ix(k),iy(k)-1,iz(k)),g,3) &
                  ) / p1m
        p2yz = 0.0 
     Else                                          ! Reflective Boundary Conditions
        yp = 1._dp
        ym = dely(iy(k)-1)/dely(iy(k))
        p1m = ym+1._dp; p2m = 2.*ym+1._dp; mm = p1m*p2m; p2 = ym+yp+2.;
        p1p = yp+1._dp; p2p = 2.*yp+1._dp; pp = p1p*p2p; p1 = pp-mm;
        p3 = p1m*p1p*(ym+yp+1._dp)
               
        P1yx = ( mm * L0(xyz(ix(k),iy(k),iz(k)),g,1)      &
               - pp * L0(xyz(ix(k),iy(k)-1,iz(k)),g,1)    &
               + p1 * L0(k,g,1)                           &
               ) / p3
        P2yx = ( p1m * L0(xyz(ix(k),iy(k),iz(k)),g,1)     &
               + p1p * L0(xyz(ix(k),iy(k)-1,iz(k)),g,1)   &
               - p2  * L0(xyz(ix(k),iy(k),iz(k)),g,1)     &
               ) / p3
        P1yz = ( mm * L0(xyz(ix(k),iy(k),iz(k)),g,3)      &
               - pp * L0(xyz(ix(k),iy(k)-1,iz(k)),g,3)    &
               + p1 * L0(k,g,3)                           &
               ) / p3
        P2yz = ( p1m * L0(xyz(ix(k),iy(k),iz(k)),g,3)     &
               + p1p * L0(xyz(ix(k),iy(k)-1,iz(k)),g,3)   &
               - p2  * l0(xyz(ix(k),iy(k),iz(k)),g,3)     &
               ) / p3
    End If
Else if (iy(k) == x_smin(ix(k))) Then           ! Y- Plan
    if (y_south == 0 .OR. y_south == 1) Then       ! Vacuum Boundary Conditions
        yp = dely(iy(k)+1)/dely(iy(k))
        p1p = yp+1._dp
        
        p1yx = 2.*( L0(xyz(ix(k),iy(k)+1,iz(k)),g,1) &
                  - L0(xyz(ix(k),iy(k),iz(k)),g,1)   &
                  ) / p1p
        p2yx = 0.0
        p1yz = 2.*( L0(xyz(ix(k),iy(k)+1,iz(k)),g,3) &
                  - L0(xyz(ix(k),iy(k),iz(k)),g,3)   &
                  ) / p1p
        p2yz = 0.0 
     Else                                          ! Reflective Boundary Conditions
        ym = 1._dp
        yp = dely(iy(k)+1)/dely(iy(k))
        p1m = ym+1._dp; p2m = 2.*ym+1._dp; mm = p1m*p2m; p2 = ym+yp+2.;
        p1p = yp+1._dp; p2p = 2.*yp+1._dp; pp = p1p*p2p; p1 = pp-mm;
        p3 = p1m*p1p*(ym+yp+1._dp)
               
        P1yx = ( mm * L0(xyz(ix(k),iy(k)+1,iz(k)),g,1)      &
               - pp * L0(xyz(ix(k),iy(k),iz(k)),g,1)        &
               + p1 * L0(k,g,1)                             &
               ) / p3
        P2yx = ( p1m * L0(xyz(ix(k),iy(k)+1,iz(k)),g,1)     &
               + p1p * L0(xyz(ix(k),iy(k),iz(k)),g,1)       &
               - p2  * L0(xyz(ix(k),iy(k),iz(k)),g,1)       &
               ) / p3
        P1yz = ( mm * L0(xyz(ix(k),iy(k)+1,iz(k)),g,3)      &
               - pp * L0(xyz(ix(k),iy(k),iz(k)),g,3)        &
               + p1 * L0(k,g,3)                             &
               ) / p3
        P2yz = ( p1m * L0(xyz(ix(k),iy(k)+1,iz(k)),g,3)     &
               + p1p * L0(xyz(ix(k),iy(k),iz(k)),g,3)       &
               - p2  * L0(xyz(ix(k),iy(k),iz(k)),g,3)       &
               ) / p3
    End if
Else
    ym = dely(iy(k)-1)/dely(iy(k))
    yp = dely(iy(k)+1)/dely(iy(k))     
    p1m = ym+1._dp; p2m = 2.*ym+1._dp; mm = p1m*p2m; p2 = ym+yp+2.    
    p1p = yp+1._dp; p2p = 2.*yp+1._dp; pp = p1p*p2p; p1 = pp-mm;
    p3 = p1m*p1p*(ym+yp+1._dp)
    
    P1yx = ( mm * L0(xyz(ix(k),iy(k)+1,iz(k)),g,1)   &
           - pp * L0(xyz(ix(k),iy(k)-1,iz(k)),g,1)   &
           + p1 * L0(k,g,1)                          &
           ) / p3
    P2yx = ( p1m * L0(xyz(ix(k),iy(k)+1,iz(k)),g,1)  &
           + p1p * L0(xyz(ix(k),iy(k)-1,iz(k)),g,1)  &
           - p2  * L0(xyz(ix(k),iy(k),iz(k)),g,1)    &
           ) / p3
    P1yz = ( mm * L0(xyz(ix(k),iy(k)+1,iz(k)),g,3)   &
           - pp * L0(xyz(ix(k),iy(k)-1,iz(k)),g,3)   &
           + p1 * L0(k,g,3)                          &
           ) / p3
    P2yz = ( p1m * L0(xyz(ix(k),iy(k)+1,iz(k)),g,3)  &
           + p1p * L0(xyz(ix(k),iy(k)-1,iz(k)),g,3)  &
           - p2  * L0(xyz(ix(k),iy(k),iz(k)),g,3)    &
           ) / p3
End if

! Set Paramaters for Z_Direction Transverse Leakage
if (iz(k) == 1) Then                            ! Z- Plan
    if (z_bott == 0 .OR. z_bott == 1) Then         ! Vacuum Boundary Conditions
        zp = delz(iz(k)+1)/delz(iz(k))
        p1p = zp+1._dp
        
        p1zx = 2.*( L0(xyz(ix(k),iy(k),iz(k)+1),g,1) &
                  - L0(xyz(ix(k),iy(k),iz(k)),g,1)   &
                  ) / p1p
        p2zx = 0.0
        p1zy = 2.*( L0(xyz(ix(k),iy(k),iz(k)+1),g,2) &
                  - L0(xyz(ix(k),iy(k),iz(k)),g,2)   &
                  ) / p1p
        p2zy = 0.0
    Else                                           ! Reflective Boundary Conditions
        zm = 1._dp
        zp = delz(iz(k)+1)/delz(iz(k))
        p1m = zm+1._dp; p2m = 2.*zm+1._dp; mm = p1m*p2m; p2 = zm+zp+2.
        p1p = zp+1._dp; p2p = 2.*zp+1._dp; pp = p1p*p2p; p1 = pp-mm
        p3 = p1m*p1p*(zm+zp+1._dp)
        
        P1zx = ( mm * L0(xyz(ix(k),iy(k),iz(k)+1),g,1)   &
               - pp * L0(xyz(ix(k),iy(k),iz(k)),g,1)     &
               + p1 * L0(k,g,1)                          &
               ) / p3
        P2zx = ( p1m * L0(xyz(ix(k),iy(k),iz(k)+1),g,1)  &
               + p1p * L0(xyz(ix(k),iy(k),iz(k)),g,1)    &
               - p2  * L0(xyz(ix(k),iy(k),iz(k)),g,1)    &
               ) / p3
        P1zy = ( mm * L0(xyz(ix(k),iy(k),iz(k)+1),g,2)   &
               - pp * L0(xyz(ix(k),iy(k),iz(k)),g,2)     &
               + p1 * L0(k,g,2)                          &
               ) / p3
        P2zy = ( p1m * L0(xyz(ix(k),iy(k),iz(k)+1),g,2)  &
               + p1p * L0(xyz(ix(k),iy(k),iz(k)),g,2)    &
               - p2  * L0(xyz(ix(k),iy(k),iz(k)),g,2)    &
               ) / p3
    End if
Else if (iz(k) == nzz) Then                     ! Z+ Plan
    if (z_top == 0 .OR. z_top == 1) Then           ! Vacuum Boundary Conditions
        zm = delz(iz(k)-1)/delz(iz(k))
        p1m = zm+1._dp
        
        p1zx = 2.*( L0(xyz(ix(k),iy(k),iz(k)),g,1)   &
                  - L0(xyz(ix(k),iy(k),iz(k)-1),g,1) &
                  ) / p1m
        p2zx = 0.0
        p1zy = 2.*( L0(xyz(ix(k),iy(k),iz(k)),g,2)   &
                  - L0(xyz(ix(k),iy(k),iz(k)-1),g,2) &
                  ) / p1m
        p2zy = 0.0
    Else                                           ! Reflective Boundary Conditions
        zp = 1._dp
        zm = delz(iz(k)-1)/delz(iz(k))
        p1m = zm+1._dp; p2m = 2.*zm+1._dp; mm = p1m*p2m; p2 = zm+zp+2.
        p1p = zp+1._dp; p2p = 2.*zp+1._dp; pp = p1p*p2p; p1 = pp-mm
        p3 = p1m*p1p*(zm+zp+1._dp)
        
        P1zx = ( mm * L0(xyz(ix(k),iy(k),iz(k)),g,1)       &
               - pp * L0(xyz(ix(k),iy(k),iz(k)-1),g,1)     &
               + p1 * L0(k,g,1)                            &
               ) / p3
        P2zx = ( p1m * L0(xyz(ix(k),iy(k),iz(k)),g,1)      &
               + p1p * L0(xyz(ix(k),iy(k),iz(k)-1),g,1)    &
               - p2  * L0(xyz(ix(k),iy(k),iz(k)),g,1)      &
               ) / p3
        P1zy = ( mm * L0(xyz(ix(k),iy(k),iz(k)),g,2)       &
               - pp * L0(xyz(ix(k),iy(k),iz(k)-1),g,2)     &
               + p1 * L0(k,g,2)                            &
               ) / p3
        P2zy = ( p1m * L0(xyz(ix(k),iy(k),iz(k)),g,2)      &
               + p1p * L0(xyz(ix(k),iy(k),iz(k)-1),g,2)    &
               - p2  * L0(xyz(ix(k),iy(k),iz(k)),g,2)      &
               ) / p3
    End if
Else
    zm = delz(iz(k)-1)/delz(iz(k))
    zp = delz(iz(k)+1)/delz(iz(k))     
    p1m = zm+1._dp; p2m = 2.*zm+1._dp; mm = p1m*p2m; p2 = zm+zp+2.    
    p1p = zp+1._dp; p2p = 2.*zp+1._dp; pp = p1p*p2p; p1 = pp-mm;
    p3 = p1m*p1p*(zm+zp+1._dp)
    
    P1zx = ( mm * L0(xyz(ix(k),iy(k),iz(k)+1),g,1)   &
           - pp * L0(xyz(ix(k),iy(k),iz(k)-1),g,1)   &
           + p1 * L0(k,g,1)                          &
           ) / p3
    P2zx = ( p1m * L0(xyz(ix(k),iy(k),iz(k)+1),g,1)  &
           + p1p * L0(xyz(ix(k),iy(k),iz(k)-1),g,1)  &
           - p2  * L0(xyz(ix(k),iy(k),iz(k)),g,1)    &
           ) / p3
    P1zy = ( mm * L0(xyz(ix(k),iy(k),iz(k)+1),g,2)   &
           - pp * L0(xyz(ix(k),iy(k),iz(k)-1),g,2)   &
           + p1 * L0(k,g,2)                          &
           ) / p3
    P2zy = ( p1m * L0(xyz(ix(k),iy(k),iz(k)+1),g,2)  &
           + p1p * L0(xyz(ix(k),iy(k),iz(k)-1),g,2)  &
           - p2  * L0(xyz(ix(k),iy(k),iz(k)),g,2)    &
           ) / p3
End if

! Set Transverse Leakage Moments
L(1) = 0.0
L(2) = ( p1xy/dely(iy(k))+p1xz/delz(iz(k)) ) / 12. ! Lx1
L(3) = ( p1yx/delx(ix(k))+p1yz/delz(iz(k)) ) / 12. ! Ly1
L(4) = ( p1zx/delx(ix(k))+p1zy/dely(iy(k)) ) / 12. ! Lz1
L(5) = ( p2xy/dely(iy(k))+p2xz/delz(iz(k)) ) / 20. ! Lx2
L(6) = ( p2yx/delx(ix(k))+p2yz/delz(iz(k)) ) / 20. ! Ly2
L(7) = ( p2zx/delx(ix(k))+p2zy/dely(iy(k)) ) / 20. ! Lz2

End Subroutine

Subroutine MatVec2(MulJin,MulMom,R2,P2,ji,Q) ! To Perform Matrix vector Multiplication
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Real(dp), Dimension(6,6),Intent(in) :: P2, R2
Real(dp), Dimension(6),Intent(in) :: ji
Real(dp), Dimension(7),Intent(in) :: Q
Real(dp), Dimension(6), Intent(out) :: MulJin, MulMom    
!Local Variables
Integer :: i, j, m, n, v, w
Real(dp) :: isum6,isum7
Real(dp), Dimension(6) :: Q_L

Q_L(:) = Q(1)

m = Size(R2,1)
v = Size(ji,1)
Do i= 1, m
    isum6 = 0.
    Do j = 1, v
        isum6 = isum6 + R2(i,j) * ji(j)
    End Do
    MulJin(i) = isum6
End Do

n = Size(P2,1)
w = Size(Q_L,1)
Do i= 1, n
    isum7 = 0.
    Do j = 1, w
        isum7 = isum7 + P2(i,j) * Q_L(j)
    End Do
    MulMom(i) = isum7
End Do

End Subroutine    

Subroutine MatVec4(MulJin,MulMom,R4,P4,ji,Q,L) ! To Perform Matrix vector Multiplication
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Real(dp), Dimension(6,6),Intent(in) :: R4
Real(dp), Dimension(6,7),Intent(in) :: P4
Real(dp), Dimension(6),Intent(in) :: ji
Real(dp), Dimension(7),Intent(in) :: Q, L
Real(dp), Dimension(6), Intent(out) :: MulJin, MulMom    
!Local Variables
Integer :: i, j, m, n, v, w
Real(dp) :: isum6,isum7
Real(dp), Dimension(7) :: Q_L

Q_L(:) = Q(:) - L(:) 

m = Size(R4,1)
v = Size(ji,1)
Do i= 1, m
    isum6 = 0.
    Do j = 1, v
        isum6 = isum6 + R4(i,j) * ji(j)
    End Do
    MulJin(i) = isum6
End Do

n = Size(P4,1)
w = Size(Q_L,1)
Do i= 1, n
    isum7 = 0.
    Do j = 1, w
        isum7 = isum7 + P4(i,j) * Q_L(j)
    End Do
    MulMom(i) = isum7
End Do

End Subroutine    

Subroutine Lxyz(L0,k,g,ng,nk,jo,ji) ! To Update Transverse Leakages for g Group and n Nod
            ! Flat Leakage Approximation   L = J(+) - J(-)
            ! J(+) = Jout - Jin  &  -J(-) = Jout - Jin  (Net Currents)
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: g, k, ng, nk
Real(dp), Dimension(nk,ng,6),Intent(in) :: jo, ji
Real(dp), Dimension(nk,ng,3),Intent(out):: L0

! 1:x+, 2:x-, 3:y+, 4:y-, 5:z+, 6:z-
L0(k,g,1) = jo(k,g,1) - ji(k,g,1) & 
          - ji(k,g,2) + jo(k,g,2)  ! Lx
L0(k,g,2) = jo(k,g,3) - ji(k,g,3) &
          - ji(k,g,4) + jo(k,g,4)  ! Ly
L0(k,g,3) = jo(k,g,5) - ji(k,g,5) &
          - ji(k,g,6) + jo(k,g,6)  ! Lz

End Subroutine

Subroutine Flux4(f0,fx1,fy1,fz1,fx2,fy2,fz2,k,g,ng,nk,nxx,nyy,nzz,&
                             D,sigr,delx,dely,delz,ix,iy,iz,jo,ji,L0,Q,L) ! To Update Nod Averaged Flux & Flux Moments
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: g, k, ng, nk, nxx, nyy, nzz
Real(dp), Dimension(nk,ng), Intent(in) :: D, sigr
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk), Intent(in) :: ix, iy, iz
Real(dp), Dimension(nk,ng,6),Intent(in) :: jo, ji
Real(dp), Dimension(nk,ng,7),Intent(in) :: Q
Real(dp), Dimension(nk,ng,3),Intent(in) :: L0
Real(dp), Dimension(nk,ng,7), Intent(in) :: L
Real(dp), Dimension(nk,ng), Intent(out) :: f0, fx1, fy1, fz1, fx2, fy2, fz2
! Local Variables
Real(dp) :: Tx, Ty, Tz

! Update Zeroth Flux (Nod Averaged Flux)
    Call Flux2(f0,fx1,fy1,fz1,fx2,fy2,fz2,k,g,ng,nk,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,sigr,Q,L0)
if (f0(k,g) < 0.) f0(k,g) = 0.

! Set Parameters Tx, Ty & Tz   [ T = (J+) + (J-) ]
Tx = jo(k,g,1) - ji(k,g,1) - jo(k,g,2) + ji(k,g,2)
Ty = jo(k,g,3) - ji(k,g,3) - jo(k,g,4) + ji(k,g,4)
Tz = jo(k,g,5) - ji(k,g,5) - jo(k,g,6) + ji(k,g,6)
   
! Calculate Flux moments
fx1(k,g) = ( Q(k,g,2) - L(k,g,2) - 0.5*Tx/delx(ix(k))             &
           - 2.*(D(k,g)/delx(ix(k))**2)                       &
           *(jo(k,g,1) + ji(k,g,1) - jo(k,g,2) - ji(k,g,2)) ) &
           / sigr(k,g)

fy1(k,g) = ( Q(k,g,3) - L(k,g,3) - 0.5*Ty/dely(iy(k))             &
           - 2.*(D(k,g)/dely(iy(k))**2)                       &
           *(jo(k,g,3) + ji(k,g,3) - jo(k,g,4) - ji(k,g,4)) ) &
           / sigr(k,g)

fz1(k,g) = ( Q(k,g,4) - L(k,g,4) - 0.5*Tz/delz(iz(k))             &
           - 2.*(D(k,g)/delz(iz(k))**2)                       &
           *(jo(k,g,5) + ji(k,g,5) - jo(k,g,6) - ji(k,g,6)) ) &
           / sigr(k,g)
	!	if (k==50) then
	!		Write(*,*) 'BEFORE CALL fx2' 	
	!		write(*,*) 'L4=', L(k,g,4) !Jin
	!		write(*,*) 'Q4=', Q(k,g,4) !Jin
	!		write(*,*) 'L5=', L(k,g,5) !Jin
	!		write(*,*) 'Q5=', Q(k,g,5) !Jin
	!	endif
fx2(k,g) = ( Q(k,g,5) - L(k,g,5) - 0.5*L0(k,g,1)/delx(ix(k))      &
           - 6.*D(k,g)/delx(ix(k))**2                       &
           * (jo(k,g,1) + ji(k,g,1) + jo(k,g,2) + ji(k,g,2)   &
           - f0(k,g)) ) / sigr(k,g)
	!	if (k==50) then
	!		Write(*,*) 'AFTER CALL fx2' 	
	!		write(*,*) 'fx2=', fx2(k,:) !Jin
	!	endif
fy2(k,g) = ( Q(k,g,6) - L(k,g,6) - 0.5*L0(k,g,2)/dely(iy(k))      &
           - 6.*D(k,g)/dely(iy(k))**2                       &
           * (jo(k,g,3) + ji(k,g,3) + jo(k,g,4) + ji(k,g,4)   &
           - f0(k,g)) ) / sigr(k,g)

fz2(k,g) = ( Q(k,g,7) - L(k,g,7) - 0.5*L0(k,g,3)/delz(iz(k))      &
           - 6.*D(k,g)/delz(iz(k))**2                       &
           * (jo(k,g,5) + ji(k,g,5) + jo(k,g,6) + ji(k,g,6)   &
           - f0(k,g)) ) / sigr(k,g)

End Subroutine

Subroutine Flux2(f0,fx1,fy1,fz1,fx2,fy2,fz2,k,g,ng,nk,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,sigr,Q,L0) ! To Update Nod Averaged Flux 
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: g, k, ng, nk, nxx, nyy, nzz
Real(dp), Dimension(nk,ng), Intent(in) :: sigr
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk), Intent(in) :: ix, iy, iz
Real(dp), Dimension(nk,ng,7),Intent(in) :: Q
Real(dp), Dimension(nk,ng,3),Intent(in) :: L0
Real(dp), Dimension(nk,ng), Intent(out) :: f0, fx1, fy1, fz1, fx2, fy2, fz2

! Calculate Zeroth Flux
f0(k,g)  = ( Q(k,g,1)             &
         - L0(k,g,1)/delx(ix(k))   &
         - L0(k,g,2)/dely(iy(k))   &
         - L0(k,g,3)/delz(iz(k)) ) &
         / sigr(k,g)
fx1(k,g) = 0.0
fy1(k,g) = 0.0
fz1(k,g) = 0.0
fx2(k,g) = 0.0
fy2(k,g) = 0.0
fz2(k,g) = 0.0

End Subroutine

Subroutine RelE(newF,oldF,rel,nk) ! To Calculate Max Relative Error for FS
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: nk
Real(dp), Dimension(nk),Intent(in) :: newF, oldF
Real(dp), Intent(out) :: rel
!Local Variables
Real(dp) :: error
Integer :: k

rel = 0.
Do k= 1, nk
    if (ABS(newF(k)) > 1.e-10_dp) Then
        error = ABS(newF(k) - oldF(k)) / ABS(newF(k))
        if (error > rel) rel = error
    End if
End Do

End Subroutine

Subroutine RelEg(newF,oldF,rel,nk,ng) ! To Calculate Max Relative error for Flux
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: ng, nk     ! Number of Groups
Real(dp), Dimension(nk,ng),Intent(in) :: newF, oldF
Real(dp), Intent(out) :: rel
!Local Variables
Real(dp) :: error
Integer :: g, k

rel = 0.
Do k= 1, nk
    Do g= 1, ng
        if (ABS(newF(k,g)) > 1.d-10) Then
            error = ABS(newF(k,g) - oldF(k,g)) / ABS(newF(k,g))
            if (error > rel) rel = error
        End if
    End Do
End Do

End Subroutine

Subroutine MultF(Keff,ng,nk,siga,nu_sigf,f0,jo,ji,nxx,nyy,nzz,ix,iy,iz,&
                 delx,dely,delz,x_smax,x_smin,y_smax,y_smin)     ! To Calculate Keff for Fixed Source Problem
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: ng, nk, nxx, nyy, nzz
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk),  Intent(in) :: ix, iy, iz
Integer, Dimension(nxx), Intent(in) :: x_smax, x_smin
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin 
Real(dp), Dimension(nk,ng), Intent(in) :: siga, nu_sigf
Real(dp), Dimension(nk,ng), Intent(in) :: f0
Real(dp), Dimension(nk,ng,6),Intent(in) :: jo, ji
Real(dp), Intent(out) :: Keff
! Local Variables
Integer :: g, k
Real(dp) :: leak, absp, fiss
!! Jot Nodals' outgoing currents  (X+, X-, Y+, Y-, Z+, Z-)
!! Jin Nodals' ingoing currents   (X+, X-, Y+, Y-, Z+, Z-)
leak = 0.0
absp = 0.0
fiss = 0.0
Do g = 1, ng
    Do k = 1, nk
        !! Get leakages
        if (ix(k) == y_smax(iy(k))) leak = leak+(jo(k,g,1)-ji(k,g,1))*dely(iy(k))*delz(iz(k))
        if (ix(k) == y_smin(iy(k))) leak = leak+(jo(k,g,2)-ji(k,g,2))*dely(iy(k))*delz(iz(k))
        if (iy(k) == x_smax(ix(k))) leak = leak+(jo(k,g,3)-ji(k,g,3))*delx(ix(k))*delz(iz(k))
        if (iy(k) == x_smin(ix(k))) leak = leak+(jo(k,g,4)-ji(k,g,4))*delx(ix(k))*delz(iz(k))
        if (iz(k) == nzz) leak = leak+(jo(k,g,5)-ji(k,g,5))*delx(ix(k))*dely(iy(k))
        if (iz(k) == 1)   leak = leak+(jo(k,g,6)-ji(k,g,6))*delx(ix(k))*dely(iy(k))
        
        absp = absp + siga(k,g) * f0(k,g) * delx(ix(k)) * dely(iy(k)) * delz(iz(k))
        fiss = fiss + nu_sigf(k,g) * f0(k,g) * delx(ix(k)) * dely(iy(k)) * delz(iz(k))

    End Do
End Do

    Keff = fiss / (leak + absp)

End Subroutine

Subroutine PowDis(p,f0,sigf,delv,ng,nk,mode) ! To Calculate Power Distribution
Implicit None
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: ng, nk
Real(dp), Dimension(nk), Intent(in) :: delv
Real(dp), Dimension(nk,ng), Intent(in) :: sigf, f0
Real(dp), Dimension(nk), Intent(out) :: p        ! Power
Character(Len=100), intent(in) :: mode
! Local Variables
Integer :: g, k
Integer, Parameter :: Ounit = 100   !output
Real(dp) :: pow, tpow

Open (Unit=ounit, File='app/Output/NEM.out')
p = 0.0
Do g= 1, ng
    Do k= 1, nk
        pow = sigf(k,g)*f0(k,g)*delv(k)
        if (pow < 0.) pow = 0.
        p(k) = p(k)+pow
    End Do
End Do

! Normalize to 1._dp
tpow = 0.
Do k = 1, nk
    tpow = tpow + p(k)
End Do

if (tpow <= 0 .AND. mode /= 'Fixed Source') Then
   Write(ounit, *) '   ERROR: TOTAL NODES POWER IS ZERO OR LESS'
   Write(ounit, *) '   STOP IN SUBROUTINE POWDIS'
   Write(*, *) '   ERROR: TOTAL NODES POWER IS ZERO OR LESS'
   Write(*, *) '   STOP IN SUBROUTINE POWDIS'
   Stop
End if

Do k = 1, nk
    p(k) = p(k) / tpow
	!	if (k < 300) then
	!		WRITE(*,*) 'itr=',l
	!		write(*,*)'node=', k,  'power=', p(k)
	!	end if
		!IF (n<150) then 
		!WRITE(*,*) 'n=', n, 'pline=',pline(n)
		!endif
End Do

End Subroutine
Subroutine AsmPow(nk, nx, ny, nxx, nyy, nzz, ix, iy, iz, delx, dely, delz, xdiv, ydiv, y_smax, y_smin, p) 
! To Print Radially Averaged Power Distribution
implicit none
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: nk, nx, ny, nxx, nyy, nzz
Integer, Dimension(nx), Intent(in) :: xdiv
Integer, Dimension(ny), Intent(in) :: ydiv
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk), Intent(in) :: ix, iy, iz, p

! Local Variables
Real*4, Dimension(100) :: sdelx
Real*4, Dimension(100) :: sdely
Real*4 :: sx, sy
Real(dp) :: som, vsom, totp, fmax
Real(dp), Dimension(nx,ny) :: fasm
Real(dp), Dimension(nxx,nyy) :: fnode
Real(dp), Dimension(nxx,nyy,nzz) :: fx
Integer, Parameter :: xm = 12, Ounit = 100
Integer :: i, j, k, n, lx, ly, ip, ipr, xs, xf, ys, yf, nfuel, xmax, ymax, m
Character(Len=6), Dimension(nx,ny) :: cpow

fx = 0._dp
Do n = 1, nk
    fx(ix(n), iy(n), iz(n)) = p(n)
End Do

! Calculate Axially Averaged Node-Wise Distribution
fnode = 0._dp
Do j = 1, nyy
    Do i = y_smin(j), y_smax(j)
        som = 0._dp
        vsom = 0._dp
        Do k = 1, nzz
            som = som + fx(i, j, k) * delz(k)
            vsom = vsom + delz(k)
        End Do
        fnode(i, j) = som / vsom
    End Do
End Do

! Calculate Assembly Power
nfuel = 0
totp = 0._dp
ys = 1
yf = 0
Do j = 1, ny
    yf = yf + ydiv(j)
    xf = 0
    xs = 1
    Do i = 1, nx
        xf = xf + xdiv(i)
        som = 0._dp
        vsom = 0._dp
        Do ly = ys, yf
            Do lx = xs, xf
                som = som + fnode(lx, ly) * delx(lx) * dely(ly)
                vsom = vsom + delx(lx) * dely(ly)
            End Do
        End Do
        fasm(i, j) = som / vsom
        xs = xs + xdiv(i)
        if (fasm(i, j) > 0._dp) nfuel = nfuel + 1
        if (fasm(i, j) > 0._dp) totp = totp + fasm(i, j)
    End Do
    ys = ys + ydiv(j)
End Do

! Normalize Assembly Power to 1.0
xmax = 1; ymax = 1
fmax = 0._dp
Do j = 1, ny
    Do i = 1, nx
        if (totp > 0.) fasm(i, j) = Real(nfuel) / totp * fasm(i, j)
        if (fasm(i, j) > fmax) Then     ! Get max position
            xmax = i
            ymax = j
            fmax = fasm(i, j)
        End if
    End Do
End Do

! Print Assembly Power Distribution
Write(ounit, *)
Write(ounit, *)
Write(ounit, *) '    Radial Power Distribution'
Write(ounit, *) '  =============================='
Open(unit = 10, File = 'app/Output/RadialPower.out')
ip = nx / xm
ipr = MOD(nx, xm) - 1
xs = 1; xf = xm
sx = 0.0
sy = 0.0

Do m = 1, xm
    sdelx(m) = sx + delx(m) * xdiv(m)
    sx = sdelx(m)
End Do
Do m = xm, 1, -1
    sdely(m) = sy + dely(m) * ydiv(m)
    sy = sdely(m)
End Do

Do k = 1, ip
    Write(ounit, '(4X,100I8)') (i, i = xs, xf)
    Write(10, '(F9.3,100F9.3)') 0.000, (sdelx(i), i = xs, xf)
    Do j = ny, 1, -1
        Write(ounit, '(2X,I4,100A8)') j, (cpow(i, j), i = xs, xf)
        Write(10, '(2X,F9.3,2X,100A8)') sdely(j), (cpow(i, j), i = xs, xf)
    End Do
    Write(ounit, *)
    xs = xs + xm
    xf = xf + xm
End Do

Write(ounit, '(4X,100I8)') (i, i = xs, xs + ipr)
Write(10, '(F9.3,100F9.3)') 0.000, (sdelx(i), i = xs, xs + ipr)

if (xs + ipr > xs) Then
    Do j = ny, 1, -1
        Write(ounit, '(2X,I4,100A8)') j, (cpow(i, j), i = xs, xs + ipr)
        Write(10, '(2X,F9.3,2X,100A8)') sdely(j), (cpow(i, j), i = xs, xs + ipr)
    End Do
End if
close(10)

! Move the file after it's written
call system('mv app/Output/RadialPower.out app/Backup/RadialPower.out')  ! Adjust for your directory structure

Write(ounit, *)
Write(ounit, *) '  MAX POS.       Maximum Value'
Write(ounit, 1101) ymax, xmax, fasm(xmax, ymax)

1101 Format(2X, '(', I3, ',', I3, ')', F15.3)
End Subroutine

Subroutine AxiPow(nk,nz,nxx,nyy,nzz,ix,iy,iz,xyz,delz,delv,zdiv,y_smax,y_smin,p)    !   To print Axially Averaged Assembly-Wise Power Distribution
implicit none
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: nk, nz, nxx, nyy, nzz
Integer, Dimension(nz), Intent(in) :: zdiv
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin
Real(dp), Dimension(nk), Intent(in) :: delv
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk), Intent(in) :: ix, iy, iz, p
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz
! Local Variables
Real(dp) :: coreh       ! Core Height
Real(dp) :: som, vsom, totp, fmax
Real(dp), Dimension(nz) :: faxi
Real(dp), Dimension(nxx,nyy,nzz) :: fx
Integer, Parameter :: Ounit = 100
Integer :: i, j, k, n, ztot, lz, nfuel, amax

fx = 0._dp
Do n = 1, nk
    fx(ix(n),iy(n),iz(n)) = p(n)
End Do

! Calculate Axial Power
nfuel = 0
totp  = 0._dp
ztot = 0
Do k= 1, nz
    som = 0._dp
    vsom = 0._dp
    Do lz= 1, zdiv(k)
        ztot = ztot + 1
        Do j = 1, nyy
            Do i = y_smin(j), y_smax(j)
                som = som + fx(i,j,ztot)
                vsom = vsom + delv(xyz(i,j,ztot))
            End Do
        End Do
    End Do
    faxi(k) = som/vsom
    if (faxi(k) > 0._dp) nfuel = nfuel + 1
    if (faxi(k) > 0._dp) totp  = totp + faxi(k)
End Do

! Normalize Axial Power to 1.0
fmax = 0._dp
amax = 1
Do k = 1, nz
    faxi(k) = Real(nfuel) / totp * faxi(k)
    if (faxi(k) > fmax) Then
        amax = k   ! Get Max Position
        fmax = faxi(k)
    End if
End Do

! Calculate Core Height
coreh = 0._dp
Do k = 1, nzz
    coreh = coreh + delz(k)
End Do

! Print Axial power distribution
Write(ounit,*)
Write(ounit,*)
Write(ounit,*) '    Axial Power Density Distribution'
Write(ounit,*) '  ===================================='
Write(ounit,*)
Write(ounit,*) '    Plane Number        Power      Height'
Write(ounit,*) '   -----------------------------------------'
!Write(*,*)
!Write(*,*)
!Write(*,*) '    Axial Power Density Distribution'
!Write(*,*) '  ===================================='
!Write(*,*)
!Write(*,*) '    Plane Number        Power      Height'
!Write(*,*) '   -----------------------------------------'
Open (unit=20, File='app/Output/AxialPower.out')
som = 0.
ztot = nzz
        Write(20,*)  0, 100

Do k= nz, 1, -1
    if (k == nz) Then
        Write(ounit,'(2X,I8,A7,F13.3, F12.2)') k, ' (Top)', faxi(k), coreh-som
   !     Write(*,'(2X,I8,A7,F13.3, F12.2)') k, ' (Top)', faxi(k), coreh-som
        Write(20,'(2X, F12.2, F10.3)')  coreh-som, faxi(k)
    Else if (k == 1) Then
        Write(ounit,'(2X,I8,A10,F10.3, F12.2)') k, ' (Bottom)', faxi(k), coreh-som
    !    Write(*,'(2X,I8,A10,F10.3, F12.2)') k, ' (Bottom)', faxi(k), coreh-som
        Write(20,'(2X, F12.2, F10.3)')  coreh-som, faxi(k)
    Else
        Write(ounit,'(2X,I8,F20.3, F12.2)') k, faxi(k), coreh-som
     !   Write(*,'(2X,I8,F20.3, F12.2)') k, faxi(k), coreh-som
        Write(20,'(2X, F12.2, F10.3)')  coreh-som, faxi(k)
    End if
    Do lz = 1, zdiv(k)
        som = som + delz(ztot)
        ztot = ztot - 1
    End Do
End Do
close(20)
Write(ounit,*)
Write(ounit,*) '  MAX POS.       Maximum Value'
Write(ounit,1102)  amax, faxi(amax)
!Write(*,*)
!Write(*,*) '  MAX POS.       Maximum Value'
!Write(*,1102)  amax, faxi(amax)

1102 Format(4X, '(' , I3, ')', F18.3)

End Subroutine
Subroutine AxiFlux(ng,nk,nz,nxx,nyy,nzz,ix,iy,iz,xyz,delz,delv,zdiv,y_smax,y_smin,p)    !   To print Axially Averaged Assembly-Wise Power Distribution
implicit none
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: ng, nk, nz, nxx, nyy, nzz
Integer, Dimension(nz), Intent(in) :: zdiv
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin
Real(dp), Dimension(nk), Intent(in) :: delv
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk), Intent(in) :: ix, iy, iz
Real(dp), Dimension(nk,ng), Intent(in) :: p
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz
! Local Variables
Real(dp) :: coreh       ! Core Height
Real(dp) :: som, vsom, totp, fmax
Real(dp), Dimension(nz,ng) :: faxi
Real(dp), Dimension(ng) :: amax
Real(dp), Dimension(nxx,nyy,nzz,ng) :: fx
Integer, Parameter :: Ounit = 100
Integer :: i, j, k, n, f, g, ztot, lz, nfuel
Character(Len=20), Dimension(ng) :: Filename ! Name of the Complete File

fx = 0._dp
Do g = 1, ng
    Do n = 1, nk
        fx(ix(n),iy(n),iz(n),g) = p(n,g)
    End Do
End Do

! Calculate Axial Power
nfuel = 0
totp  = 0._dp
ztot = 0
Do g = 1, ng
    Do k= 1, nz
        som = 0._dp
        vsom = 0._dp
        Do lz= 1, zdiv(k)
            ztot = ztot + 1
            Do j = 1, nyy
                Do i = y_smin(j), y_smax(j)
                    som = som + fx(i,j,ztot,g)
                    vsom = vsom + delv(xyz(i,j,ztot))
                End Do
            End Do
        End Do
        faxi(k,g) = som/vsom
        if (faxi(k,g) > 0._dp) nfuel = nfuel + 1
        if (faxi(k,g) > 0._dp) totp  = totp + faxi(k,g)
    End Do
End Do

! Normalize Axial Power to 1.0
fmax = 0._dp
amax = 1
Do g = 1, ng
    Do k = 1, nz
        faxi(k,g) = Real(nfuel) / totp * faxi(k,g)
        if (faxi(k,g) > fmax) Then
            amax(g) = k   ! Get Max Position
            fmax = faxi(k,g)
        End if
    End Do
End Do

! Calculate Core Height
coreh = 0._dp
Do k = 1, nzz
    coreh = coreh + delz(k)
End Do

Do f=1,ng ! Loop through all Files
    Write(Filename(f),'(A17,I2)') 'app\Output\AxiFluxGr',f
End Do

Do g = 1, ng
    Write(ounit,'(A,I3)') '    Group : ', g
    Write(*,'(A,I3)')     '    Group : ', g
    ! Print Axial Flux Distribution
    Write(ounit,*)
    Write(ounit,*)
    Write(ounit,*) '    Axial Flux Density Distribution'
    Write(ounit,*) '  ===================================='
    Write(ounit,*)
    Write(ounit,*) '    Plane Number        Flux      Height'
    Write(ounit,*) '   -----------------------------------------'
  !  Write(*,*)
  !  Write(*,*)
  !  Write(*,*) '    Axial Flux Density Distribution'
  !  Write(*,*) '  ===================================='
  !  Write(*,*)
  !  Write(*,*) '    Plane Number        Flux      Height'
  !  Write(*,*) '   -----------------------------------------'
    Open (g,File=Filename(g))
    som = 0.
    ztot = nzz
    Write(Filename(g),*)  0, 100
    Do k= nz, 1, -1
        if (k == nz) Then
            Write(ounit,'(2X,I8,A7,F13.3, F12.2)') k, ' (Top)', faxi(k,g), coreh-som
        !    Write(*,'(2X,I8,A7,F13.3, F12.2)') k, ' (Top)', faxi(k,g), coreh-som
            Write(Filename(g),'(2X, F12.2, F10.3)')  coreh-som, faxi(k,g)
        Else if (k == 1) Then
            Write(ounit,'(2X,I8,A10,F10.3, F12.2)') k, ' (Bottom)', faxi(k,g), coreh-som
        !    Write(*,'(2X,I8,A10,F10.3, F12.2)') k, ' (Bottom)', faxi(k,g), coreh-som
            Write(Filename(g),'(2X, F12.2, F10.3)')  coreh-som, faxi(k,g)
        Else
            Write(ounit,'(2X,I8,F20.3, F12.2)') k, faxi(k,g), coreh-som
        !    Write(*,'(2X,I8,F20.3, F12.2)') k, faxi(k,g), coreh-som
            Write(Filename(g),'(2X, F12.2, F10.3)')  coreh-som, faxi(k,g)
        End if
        Do lz = 1, zdiv(k)
            som = som + delz(ztot)
            ztot = ztot - 1
        End Do
    End Do
End Do

close(20)
Do g = 1, ng
    Write(ounit,'(A,I3)') '    Group : ', g
 !   Write(*,'(A,I3)')     '    Group : ', g
    Write(ounit,*)
    Write(ounit,*) '  MAX POS.       Maximum Value'
    Write(ounit,1102)  amax(g), faxi(amax(g),g)
  !  Write(*,*)
  !  Write(*,*) '  MAX POS.       Maximum Value'
  !  Write(*,1102)  amax(g), faxi(amax(g),g)
End Do

1102 Format(4X, '(' , I3, ')', F18.3)

End Subroutine
Subroutine AsmFlux(ng,nk,nx,ny,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,xdiv,ydiv,y_smax,y_smin,p,norm)   !   To Print Radially Averaged Flux Distribution
implicit none
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: ng, nk, nx, ny, nxx, nyy, nzz
Integer, Dimension(nx), Intent(in) :: xdiv
Integer, Dimension(ny), Intent(in) :: ydiv
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk), Intent(in) :: ix, iy, iz
Real(dp), Dimension(nk,ng), Intent(in) :: p
Real(dp), Intent(in) :: norm
! Local Variables
Real(dp) :: som, vsom
Real*4, Dimension(100) :: sdelx
Real*4, Dimension(100) :: sdely
Real*4 :: sx, sy
Real(dp), Dimension(ng) :: totf
Real(dp), Dimension(nx,ny,ng) :: fasm
Real(dp), Dimension(nxx,nyy,ng) :: fnode
Real(dp), Dimension(nxx,nyy,nzz,ng) :: fx
Integer, Parameter :: xm = 12, Ounit = 100
Integer :: i, j, k, n, g, lx, ly, ip, ipr, xs, xf, ys, yf, negF, m, f
Character(Len=10), Dimension(nx,ny) :: cflx
Character(Len=20), Dimension(ng) :: Filename ! Name of the Complete File

fx = 0._dp
Do g = 1, ng
    Do n = 1, nk
        fx(ix(n),iy(n),iz(n),g) = p(n,g)
    End Do
End Do

!Calculate Axially Averaged Node-Wise Distribution
fnode = 0._dp
Do g = 1, ng
    Do j = 1, nyy
        Do i = y_smin(j), y_smax(j)
            som = 0._dp
            vsom = 0._dp
            Do k = 1, nzz
                som = som + fx(i,j,k,g)*delx(i)*dely(j)*delz(k)
                vsom = vsom + delx(i)*dely(j)*delz(k)
            End Do
            fnode(i,j,g)= som/vsom
        End Do
    End Do
End Do

!Calculate Radial Flux (Assembly Wise)
negF = 0
Do g = 1, ng
    totf(g)  = 0._dp
    ys = 1
    yf = 0
    Do j= 1, ny
        yf = yf + ydiv(j)
        xf = 0
        xs = 1
        Do i= 1, nx
            xf = xf + xdiv(i)
            som = 0._dp
            vsom = 0._dp
            Do ly= ys, yf
                Do lx= xs, xf
                    som = som + fnode(lx,ly,g)*delx(lx)*dely(ly)
                    vsom = vsom + delx(lx)*dely(ly)
                End Do
            End Do
            fasm(i,j,g) = som / vsom
            xs = xs + xdiv(i)
            if (fasm(i,j,g) > 0._dp) totf(g) = totf(g) + fasm(i,j,g)
            if (fasm(i,j,g) < 0._dp) negF = 1   ! Check if there is Negative Flux
        End Do
        ys = ys + ydiv(j)
    End Do
End Do

! Normalize Flux to Norm
Do g = 1, ng
    Do j = 1, ny
        Do i = 1, nx
            fasm(i,j,g) = norm / totf(g) * fasm(i,j,g) * norm
        End Do
    End Do
End Do

! Print Assembly Flux Distribution

Write(ounit,*)
Write(*,*)
if (negf > 0) Then
    Write(ounit,*) '    ....Warning: Negative Flux Encountered....'
    Write(*,*)     '    ....Warning: Negative Flux Encountered....'
End if

Write(ounit,*) '    Radial Flux Distribution'
Write(ounit,*) '  =============================='
!Write(*,*)     '    Radial Flux Distribution'
!Write(*,*)     '  =============================='

ip = nx/xm
ipr = Mod(nx,xm) - 1
sx = 0.0
sy = 0.0

Do m = 1, xm
    sdelx(m) = sx+delx(m)*xdiv(m)
    sx = sdelx(m)
End Do
Do m = xm, 1, -1
    sdely(m) = sy+dely(m)*ydiv(m)
    sy = sdely(m)
End Do

Do f=1,ng ! Loop through all Files
    Write(Filename(f),'(A17,I2)') 'app/Output/FluxGr',f
End Do

Do g = 1, ng
    Write(ounit,'(A,I3)') '    Group : ', g
   ! Write(*,'(A,I3)')     '    Group : ', g
    ! If Flux = 0 ==> Blank Spaces
    Do j = 1, ny
        Do i = 1, nx
            ! if ((fasm(i,j,g) - 0.) < 1.e-5_dp) Then
                ! cflx(i,j) = '         '
            ! Else
                Write (cflx(i,j),'(ES10.3)') fasm(i,j,g)
                cflx(i,j) = TRIM(ADJUSTL(cflx(i,j)))
            ! End if
        End Do
    End Do

    xs = 1; xf = xm
    ! Do f=1,ng
        Open (g,File=Filename(g))
        Do k = 1, ip
            Write(ounit,'(2X,100I11)') (i, i = xs, xf)
    !        Write(*,'(2X,100I11)') (i, i = xs, xf)
            Write(g,'(F9.3,100F9.3)') 0.000, ( sdelx(i) , i = xs, xf)

            Do j= ny, 1, -1
                Write(ounit,'(2X,I4,2X,100A11)') j, (cflx(i,j), i=xs, xf)
       !         Write(*,'(2X,I4,2X,100A11)') j, (cflx(i,j), i=xs, xf)
                Write(g,'(2X,F9.3,2X,100A11)') sdely(j), (cflx(i,j), i=xs, xf)
            End Do
            Write(ounit,*)
            Write(*,*)
            xs = xs + xm
            xf = xf + xm
        End Do

    Write(ounit,'(3X,100I11)') (i, i = xs, xs+ipr)
   ! Write(*,'(3X,100I11)') (i, i = xs, xs+ipr)
    Write(g,'(F9.3,100F9.3)') 0.000, ( sdelx(i), i = xs, xs+ipr)
    if (xs+ipr > xs) Then
        Do j= ny, 1, -1
            Write(ounit,'(2X,I4,2X,100A11)') j, (cflx(i,j), i=xs, xs+ipr)
      !      Write(*,'(2X,I4,2X,100A11)') j, (cflx(i,j), i=xs, xs+ipr)
            Write(g,'(2X,F9.3,2X,100A11)') sdely(j), (cflx(i,j), i=xs, xs+ipr)
        End Do
    End if
    Write(ounit,*)
 !   Write(*,*)
    Close(f)
End Do

1101 Format(2X, '(' , I3, ',', I3,')', F15.3)

End Subroutine
Subroutine AsmFlux_FxS(ng,nk,nx,ny,nxx,nyy,nzz,ix,iy,iz,delx,dely,delz,xdiv,ydiv,y_smax,y_smin,p)   !   To Print Radially Averaged Flux Distribution for Fixed Sources
implicit none
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: ng, nk, nx, ny, nxx, nyy, nzz
Integer, Dimension(nx), Intent(in) :: xdiv
Integer, Dimension(ny), Intent(in) :: ydiv
Integer, Dimension(nyy), Intent(in) :: y_smax, y_smin
Real(dp), Dimension(nxx), Intent(in) :: delx
Real(dp), Dimension(nyy), Intent(in) :: dely
Real(dp), Dimension(nzz), Intent(in) :: delz
Real(dp), Dimension(nk), Intent(in) :: ix, iy, iz
Real(dp), Dimension(nk,ng), Intent(in) :: p
! Local Variables
Real(dp) :: som, vsom
Real(dp), Dimension(ng) :: totf
Real(dp), Dimension(nx,ny,ng) :: fasm
Real(dp), Dimension(nxx,nyy,ng) :: fnode
Real(dp), Dimension(nxx,nyy,nzz,ng) :: fx
Integer, Parameter :: xm = 12, Ounit = 100
Integer :: i, j, k, n, g, lx, ly, ip, ipr, xs, xf, ys, yf, negF 
Character(Len=10), Dimension(nx,ny) :: cflx

fx = 0._dp
Do g = 1, ng
    Do n = 1, nk
        fx(ix(n),iy(n),iz(n),g) = p(n,g)
    End Do
End Do

!Calculate Axially Averaged Node-Wise Distribution
fnode = 0._dp
Do g = 1, ng
    Do j = 1, nyy
        Do i = y_smin(j), y_smax(j)
            som = 0._dp
            vsom = 0._dp
            Do k = 1, nzz
                som = som + fx(i,j,k,g)*delx(i)*dely(j)*delz(k)
                vsom = vsom + delx(i)*dely(j)*delz(k)
            End Do
            fnode(i,j,g)= som/vsom
        End Do
    End Do
End Do

!Calculate Radial Flux (Assembly Wise)
negF = 0
Do g = 1, ng
    totf(g)  = 0._dp
    ys = 1
    yf = 0
    Do j= 1, ny
        yf = yf + ydiv(j)
        xf = 0
        xs = 1
        Do i= 1, nx
            xf = xf + xdiv(i)
            som = 0._dp
            vsom = 0._dp
            Do ly= ys, yf
                Do lx= xs, xf
                    som = som + fnode(lx,ly,g)*delx(lx)*dely(ly)
                    vsom = vsom + delx(lx)*dely(ly)
                End Do
            End Do
            fasm(i,j,g) = som / vsom
            xs = xs + xdiv(i)
            if (fasm(i,j,g) > 0._dp) totf(g) = totf(g) + fasm(i,j,g)
            if (fasm(i,j,g) < 0._dp) negF = 1   ! Check if there is Negative Flux
        End Do
        ys = ys + ydiv(j)
    End Do
End Do

! Print Assembly Flux Distribution
Write(ounit,*)
if (negf > 0) Write(ounit,*) '    ....Warning: Negative Flux Encountered....'
Write(ounit,*) '    Radial Flux Distribution'
Write(ounit,*) '  =============================='

ip = nx/xm
ipr = Mod(nx,xm) - 1
Do g = 1, ng
    Write(ounit,'(A,I3)') '    Group : ', g
    ! If Flux = 0 ==> Blank Spaces
    Do j = 1, ny
        Do i = 1, nx
            if ((fasm(i,j,g) - 0.) < 1.e-5_dp) THEN
                cflx(i,j) = '         '
            Else
                Write (cflx(i,j),'(ES10.3)') fasm(i,j,g)
                cflx(i,j) = TRIM(ADJUSTL(cflx(i,j)))
            End if
        End Do
    End Do
    
    xs = 1; xf = xm
    Do k = 1, ip
        Write(ounit,'(2X,100I11)') (i, i = xs, xf)
        Do j= ny, 1, -1
            Write(ounit,'(2X,I4,2X,100A11)') j, (cflx(i,j), i=xs, xf)
        End Do
        Write(ounit,*)
        xs = xs + xm
        xf = xf + xm
    End Do

    Write(ounit,'(3X,100I11)') (i, i = xs, xs+ipr)
    if (xs+ipr > xs) Then
        Do j= ny, 1, -1
            Write(ounit,'(2X,I4,2X,100A11)') j, (cflx(i,j), i=xs, xs+ipr)
        End Do
    End if
    Write(ounit,*)
End Do

1101 Format(2X, '(' , I3, ',', I3,')', F15.3)

End Subroutine

SUBROUTINE par_ave_f(ave,par,nk,ng,nu_sigf,delv)!    To calculate average fuel temp (only for active core)

IMPLICIT NONE
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: nk, ng
REAL(DP), DIMENSION(nk,ng), INTENT(IN) :: nu_sigf            ! nu* fission macroscopic cx
Real(dp), Dimension(nk), Intent(in) :: delv
REAL(DP), DIMENSION(nk), INTENT(IN) :: par
REAL(DP), INTENT(OUT) :: ave
REAL(DP) :: dum, dum2
INTEGER :: n

dum = 0.; dum2 = 0.
DO n = 1, nk
   IF (nu_sigf(n,ng) > 0.) THEN
      dum = dum + par(n) * delv(n)
      dum2 = dum2 + delv(n)
   END IF
END DO

ave = dum / dum2

END SUBROUTINE 
SUBROUTINE par_ave(ave,par,nk,delv)!    To calculate average moderator temp (only for radially active core)
IMPLICIT NONE
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: nk
Real(dp), Dimension(nk), Intent(in) :: delv
REAL(DP), DIMENSION(nk), INTENT(IN) :: par
REAL(DP), INTENT(OUT) :: ave
REAL(DP) :: dum, dum2
INTEGER :: n

dum = 0.; dum2 = 0.
DO n = 1, nk
   dum = dum + par(n) * delv(n)
   dum2 = dum2 + delv(n)
END DO

ave = dum / dum2

END SUBROUTINE par_ave
SUBROUTINE par_max(pmax,par,nk)!    To calculate maximum fuel tem, coolant tem and density
IMPLICIT NONE
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: nk
REAL(DP), DIMENSION(nk), INTENT(IN) :: par
REAL(DP), INTENT(OUT) :: pmax
INTEGER :: n

pmax = 0.
DO n = 1, nk
   IF (par(n) > pmax) pmax = par(n)
END DO

END SUBROUTINE
SUBROUTINE  GetFq(nk,nxx,nyy,nzz,ix,iy,iz,delz,node_nf,fn)!    To get Heat Flux Hot Channel Factor
!    The maximum local linear power density in the core divided by the core average fuel rod linear power density.

IMPLICIT NONE
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15), Ounit = 100
Integer, Intent(in) :: nk, nxx, nyy, nzz
Real(dp), Dimension(nk), Intent(in) :: ix, iy, iz
Real(dp), Dimension(nzz), Intent(in) :: delz
REAL(DP), DIMENSION(nk), INTENT(IN) :: fn           ! Relative Power
REAL(DP), DIMENSION(nxx, nyy), intent(in) :: node_nf       ! Number of fuel pin per node

INTEGER :: n
REAL(DP), DIMENSION(nk) :: locp, xf
REAL(DP) :: totp, npmax, tleng, pave


totp = 0.; tleng = 0.
DO n = 1, nk
    IF (fn(n) > 0.) THEN
        xf(n) = fn(n) / node_nf(ix(n),iy(n))
        totp = totp + xf(n)
        tleng = tleng + delz(iz(n))
        locp(n) = xf(n) / delz(iz(n))
    END IF
END DO

pave = totp / tleng

npmax = 0.
DO n = 1, nk
    IF (fn(n) > 0.) THEN
       locp(n) = locp(n) / pave
       IF (locp(n) > npmax) npmax = locp(n)
    END IF
END DO


WRITE(ounit,*)
WRITE(ounit, 4001) npmax

4001 FORMAT (2X, ' HEAT FLUX HOT CHANNEL FACTOR :', F7.3)


END SUBROUTINE


Subroutine ExSrc(eSrc,nk,ng,nx,ny,nz,nxx,nyy,nzz,nS,dS,spS,&
                 xyz,xdiv,ydiv,zdiv,xpos,ypos,zpos,npox,npor) !   To Read Extra Sources if any
implicit none
Integer, Parameter :: dp = SELECTED_REAL_KIND(10, 15)
Integer, Intent(in) :: nk, ng, nx, ny, nz, nxx, nyy, nzz, npox, npor
Integer, Dimension(nx), Intent(in) :: xdiv
Integer, Dimension(ny), Intent(in) :: ydiv
Integer, Dimension(nz), Intent(in) :: zdiv
Integer, Dimension(nxx,nyy,nzz), Intent(in) :: xyz
Integer, Intent(in) :: nS
Real(dp), Dimension(nS), Intent(in) :: dS
Real(dp), Dimension(nS,ng), Intent(in) :: spS
Integer, Dimension(nS,npox), Intent(in) :: zpos
Integer, Dimension(nS,npor), Intent(in) :: xpos, ypos
Real(dp), Dimension(nk,ng), Intent(out) :: eSrc
! Local Variables
Character(Len=1), Dimension(nx,ny) :: posS         ! Source position
! Character(Len=2), Dimension(nxx, nyy) :: mmap
Integer, Parameter :: xm = 36, Ounit = 100
Integer :: i, j, k, g, h, n, xt, yt, zt, it, jt, kt, r, z
Real(dp) :: som

Open (Unit=Ounit, File='app/Output/NEM.out') 
Write(Ounit,*)
Write(Ounit,*)
Write(Ounit,*) '           >>>>> Reading Extra Sources <<<<<'
Write(Ounit,*) '           -------------------------------'
eSrc = 0._dp
Do n = 1, nS
    if (dS(n) <= 0.0) Then
        Write(ounit,*) '  Error: Source Density Shall Be Greater Than Zero'
        Stop
    End if

    ! Is Total Spectrum = 1.0?
    som = 0._dp
    Do g = 1, ng
        som = som + spS(n,g)
    End Do
    ! Check Total Spectrum
    if (ABS(som - 1._dp) > 1.e-5_dp) Then
        Write(ounit,*) 'Total Source Spectrum Is Not Equal To 1.0'
        Stop
    End if

    ! Write Output
    Write(Ounit,'(A12,I3)') '     Source ', n
    Write(Ounit,*)         '-----------------'
    Write(Ounit,'(A20,ES10.3, A11)') '  Source Density  : ', dS(n), '  n/(m^3*s)'
    Write(Ounit,'(A19,100F6.2)') '  Source Spectrum : ', (spS(n,g), g = 1, ng)
    Write(Ounit,*) ' Source Position '
    ! Read Source Position
! v = size(zpos,1)
! w = size(xpos,1)
    Do z = 1, npox
            if (zpos(n,z) < 1) Exit
            if (zpos(n,z) > nz) Then
                Write(Ounit,* ) '  Error: Wrong Extra Sources Position (ZPOS)'
                Write(Ounit, 2033) zpos(n,z)
                Stop
            End if
            posS = '0'
            Do r = 1, npor
                if (xpos(n,r) < 1 .OR. ypos(n,r) < 1) Exit
                if (xpos(n,r) > nx) Then
                    Write(ounit,* ) '  Error: Wrong Extra Sources Position (XPos)'
                    Write(ounit, 2033) xpos(n,r), ypos(n,r)
                    Stop
                End if
                if (ypos(n,r) > ny) Then
                    Write(ounit,* ) '  Error: Wrong Extra Sources Position (YPos)'
                    Write(ounit, 2033) xpos(n,r), ypos(n,r)
                    Stop
                End if
                posS(xpos(n,r),ypos(n,r)) = 'X'
                zt = 0
                kt = 1
                Do k = 1, zpos(n,z)
                    if (k > 1) kt = zt + 1
                    zt = zt + zdiv(k)
                End Do
                yt = 0
                jt = 1
                Do j = 1, ypos(n,r)
                    if (j > 1) jt = yt + 1
                    yt = yt + ydiv(j)
                End Do
                xt = 0
                it = 1
                Do i = 1, xpos(n,r)
                    if (i > 1) it = xt + 1
                    xt = xt + xdiv(i)
                End Do
                Do k = kt, zt
                    Do j = jt, yt
                        Do i = it, xt
                            Do h = 1, ng
                                eSrc(xyz(i,j,k),h) = eSrc(xyz(i,j,k),h) + dS(n) * spS(n,h)
                            End Do
                        End Do
                    End Do
                End Do
            End Do
            Write(Ounit,'(A18,I3)') '   Plane Number : ', zpos(n,z)
            Write(Ounit,'(7X,100I3)') (i, i = 1, nx)
            Do j = ny, 1, -1
                Write(Ounit,'(4X,I3, 100A3 )') j, (posS(i,j), i=1, nx)
            End Do
            Write(Ounit,*)
        End Do
    End Do

2033 Format(2X, I3, I3)

End Subroutine

Subroutine timestamp()
!      ------------------------------------------------------------------------
!      TIMESTAMP prints the current YMDHMS date as a time stamp.
!      Example:
!      31 May 2001   9:45:54.872 AM
!      Licensing:
!      This code is distributed under the GNU LGPL license.
!      Modified:
!      18 May 2013
!      Author:
!      John Burkardt
!      Parameters:
!      None
!      ------------------------------------------------------------------------
implicit none
Integer, Parameter ::  Ounit = 100
Character(len=8) :: ampm
Integer(kind=4) :: d
Integer(kind=4) :: h
Integer(kind=4) :: m
Integer(kind=4) :: mm
Character(len=9), Parameter, Dimension(12) :: month = (/ &
       'January  ', 'February ', 'March    ', 'April    ', &
       'May      ', 'June     ', 'July     ', 'August   ', &
       'September', 'October  ', 'November ', 'December ' /)
Integer(kind=4) :: n
Integer(kind=4) :: s
Integer(kind=4) :: values(8)
Integer(kind=4) :: y

    call date_and_time (values = values)
    y = values(1)
    m = values(2)
    d = values(3)
    h = values(5)
    n = values(6)
    s = values(7)
    mm = values(8)
    if (h < 12) Then
       ampm = 'AM'
    Else if (h == 12) Then
        if (n == 0 .and. s == 0) Then
            ampm = 'Noon'
        Else
            ampm = 'PM'
        End if
    Else
        h = h - 12
        if ( h < 12 ) Then
            ampm = 'PM'
        Else if ( h == 12 ) Then
            if (n == 0 .And. s == 0) Then
                ampm = 'Midnight'
            else
                ampm = 'AM'
            End if
        End if
    End if

    Write(*,*)
    write (*,'(2X,i6,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)') &
    d, trim (month(m)), y, h, ':', n, ':', s, '.', mm, trim (ampm)
    Return
    Write(*,*)
    Write(*,*)'             _______________________             '  

    Write(Ounit,*)
    write (Ounit,'(2X,i6,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)') &
    d, trim (month(m)), y, h, ':', n, ':', s, '.', mm, trim (ampm)
    Return
    Write(Ounit,*)
    Write(Ounit,*)'             _______________________             '  
End subroutine

Subroutine title1() 
implicit none
Integer, Parameter ::  Ounit = 100
Open (Unit=Ounit, File='app/Output/NEM.out') 
    Write(*,*) ''
    ! Write(*,*) 'CarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEM' &
    ! 'CarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEM'
    ! Write(*,*) ''
    ! Write(*,*) '     CCCCCCCCCC                                            NN           NN      EEEEEEEEE     MM             MM    '
    ! Write(*,*) '    CCCCCCCCCCCCC                                         NNNN         NNN     EEEEEEEEE     MMMM           MMMM   '
    ! Write(*,*) '   CCCC      CCCCC                                       NNNNNN       NNNN    EEEE          MMMMMM         MMMMMM  '
    ! Write(*,*) '  CCCC        CCCC                                       NNNNNNN      NNNN    EEEE          MMMMMMM       MMMMMMM  '
    ! Write(*,*) ' CCCC           CC                                       NNNN  NN     NNNN    EEEE          MMMM  MMM   MMM  MMMM  '
    ! Write(*,*) ' CCCC                                                    NNNN  NNN    NNNN    EEEE          MMMM   MMMMMMM   MMMM  '
    ! Write(*,*) ' CCCC                    AAAAAAAAAA         RRRRRRRRR    NNNN   NNN   NNNN    EEEEEEEEE     MMMM    MMMMM    MMMM  '
    ! Write(*,*) ' CCCC                   AAAAAAAAAAAA       RRRRRRRR      NNNN   NNN   NNNN    EEEEEEEEEE    MMMM     MMM     MMMM  '
    ! Write(*,*) ' CCCC                  AAAA      AAAA     RRR            NNNN    NNN  NNNN    EEEEEEEEE     MMMM             MMMM  '
    ! Write(*,*) ' CCCC           CC    AAAA       AAAA    RRRR            NNNN     NN  NNNN    EEEE          MMMM             MMMM  '
    ! Write(*,*) '  CCCC        CCCC    AAAA       AAAA    RRRR            NNNN      NNNNNNN    EEEE          MMMM             MMMM  '
    ! Write(*,*) '   CCCC      CCCCC     AAAA      AAAA    RRRR            NNNN       NNNNNN    EEEE          MMMM             MMMM  '
    ! Write(*,*) '    CCCCCCCCCCCCC       AAAAAAAAAAAAA    RRRR            NNN         NNNN      EEEEEEEEE    MMM               MMM  '
    ! Write(*,*) '     CCCCCCCCCC          AAAAAAAAAAAA    RRRR            NN           NN        EEEEEEEEE   MM                 MM  '
    ! Write(*,*)
    ! Write(*,*) 'CarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEM' &
    ! 'CarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEM'




    Write(*,*) '                                                                               '
    Write(*,*) '     /\               /\                                                       '
    Write(*,*) '    //\\             //\\                                                      '
    Write(*,*) '   //**\\           //**\\                                                     '
    Write(*,*) '  //****\\         //****\\                                                    '
    Write(*,*) ' //******\\       //******\\   ____                   _   _           _        '  
    Write(*,*) ' ---------   /*\   ---------  / __ \                 | \ | |         | |       '
    Write(*,*) '            (***)            | |  | |_ __   ___ _ __ |  \| | ___   __| | ___   '
    Write(*,*) '             \*/             | |  | | "_ \ / _ \ "_ \| . ` |/ _"\ / _" |/ _ \  '
    Write(*,*) '                             | |__| | |_) |  __/ | | | |\  | |_) | (_| |  __/  '
    Write(*,*) '             //\\             \____/| .__/ \___)_| |_|_| \_|\___/ \___/ \___|  '
    Write(*,*) '            //**\\                  | |                                        '
    Write(*,*) '           //****\\                 |_|                                        '
    Write(*,*) '          //******\\                                                           '
    Write(*,*) '          ----------                                                           '
    Write(*,*) '                                                                               '
    Write(*,*) '____________________________________________________________'
    Write(*,*) '         | The OpenNode, Nodal Methods for Neutron Diffusion'       
    Write(*,*) 'Version  | Version Number: 1.1                              '     
    Write(*,*) 'Copyright| 2021 Radiation & Nuclear System Laboratory       '
    Write(*,*) '         | University Abdelmalek Essaadi                    '
    Write(*,*) '         | Faculty of Sciences, Tetouan, Morocco            '
    Write(*,*) 'Source   | FORTRAN90 version                                ' 
    Write(*,*) 'GUI      | PyQt5                                            ' 
    Write(*,*) 'Method   | The Nodal Expansion Method (NEM)                 '  
    Write(*,*) 'Dimension| Three Dimensions (3D)                            '
    Write(*,*) 'Geometry | Cartesian                                        ' 
    Write(*,*) '____________________________________________________________'
    ! Write(*,*) ''
    ! Write(*,*) '     CCCCCCCCCC                                          &
       ! NN           NN      EEEEEEEEE     MM             MM    '
    ! Write(*,*) '    CCCCCCCCCCCCC                                       &
     ! NNNN         NNN     EEEEEEEEE     MMMM           MMMM   '
    ! Write(*,*) '   CCCC      CCCCC                                     &
    ! NNNNNN       NNNN    EEEE          MMMMMM         MMMMMM  '
    ! Write(*,*) '  CCCC        CCCC                                     &
    ! NNNNNNN      NNNN    EEEE          MMMMMMM       MMMMMMM  '
    ! Write(*,*) ' CCCC           CC                                     &
    ! NNNN  NN     NNNN    EEEE          MMMM  MMM   MMM  MMMM  '
    ! Write(*,*) ' CCCC                                                  &
    ! NNNN  NNN    NNNN    EEEE          MMMM   MMMMMMM   MMMM  '
    ! Write(*,*) ' CCCC                    AAAAAAAAAA         RRRRRRRRR  &
    ! NNNN   NNN   NNNN    EEEEEEEEE     MMMM    MMMMM    MMMM  '
    ! Write(*,*) ' CCCC                   AAAAAAAAAAAA       RRRRRRRR    &
    ! NNNN   NNN   NNNN    EEEEEEEEEE    MMMM     MMM     MMMM  '
    ! Write(*,*) ' CCCC                  AAAA      AAAA     RRR          &
    ! NNNN    NNN  NNNN    EEEEEEEEE     MMMM             MMMM  '
    ! Write(*,*) ' CCCC           CC    AAAA       AAAA    RRRR          &
    ! NNNN     NN  NNNN    EEEE          MMMM             MMMM  '
    ! Write(*,*) '  CCCC        CCCC    AAAA       AAAA    RRRR          &
    ! NNNN      NNNNNNN    EEEE          MMMM             MMMM  '
    ! Write(*,*) '   CCCC      CCCCC     AAAA      AAAA    RRRR          &
    ! NNNN       NNNNNN    EEEE          MMMM             MMMM  '
    ! Write(*,*) '    CCCCCCCCCCCCC       AAAAAAAAAAAAA    RRRR          &
    ! NNN         NNNN      EEEEEEEEE    MMM               MMM  '
    ! Write(*,*) '     CCCCCCCCCC          AAAAAAAAAAAA    RRRR          &
    ! NN           NN        EEEEEEEEE   MM                 MM  '
    Write(ounit,*) '                                                                               '
    Write(ounit,*) '     /\               /\                                                       '
    Write(ounit,*) '    //\\             //\\                                                      '
    Write(ounit,*) '   //**\\           //**\\                                                     '
    Write(ounit,*) '  //****\\         //****\\                                                    '
    Write(ounit,*) ' //******\\       //******\\   ____                   _   _           _        '  
    Write(ounit,*) ' ---------   /*\   ---------  / __ \                 | \ | |         | |       '
    Write(ounit,*) '            (***)            | |  | |_ __   ___ _ __ |  \| | ___   __| | ___   '
    Write(ounit,*) '             \*/             | |  | | "_ \ / _ \ "_ \| . ` |/ _"\ / _" |/ _ \  '
    Write(ounit,*) '                             | |__| | |_) |  __/ | | | |\  | |_) | (_| |  __/  '
    Write(ounit,*) '             //\\             \____/| .__/ \___)_| |_|_| \_|\___/ \___/ \___|  '
    Write(ounit,*) '            //**\\                  | |                                        '
    Write(ounit,*) '           //****\\                 |_|                                        '
    Write(ounit,*) '          //******\\                                                           '
    Write(ounit,*) '          ----------                                                           '
    Write(ounit,*) '                                                                               '
    Write(ounit,*) '____________________________________________________________'
    Write(ounit,*) '         | The OpenNode, Nodal Methods for Neutron Diffusion'       
    Write(ounit,*) 'Version  | Version Number: 1.1                              '     
    Write(ounit,*) 'Copyright| 2021 Radiation & Nuclear System Laboratory       '
    Write(ounit,*) '         | University Abdelmalk Essaadi                     '
    Write(ounit,*) '         | Faculty of Sciences, Tetouan, Morocco            '
    Write(ounit,*) 'Source   | FORTRAN90 version                                ' 
    Write(ounit,*) 'GUI      | PyQt5                                            ' 
    Write(ounit,*) 'Method   | The Nodal Expansion Method (NEM)                 '  
    Write(ounit,*) 'Dimension| Three Dimensions (3D)                            '
    Write(ounit,*) 'Geometry | Cartesian                                        ' 
    Write(ounit,*) '____________________________________________________________'
    ! Write(ounit,*) ''
    ! Write(ounit,*) 'CarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEM' &
    ! 'CarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEM'
    ! Write(ounit,*) ''    
    ! Write(ounit,*) '     CCCCCCCCCC                                          &
       ! NN           NN      EEEEEEEEE     MM             MM    '
    ! Write(ounit,*) '    CCCCCCCCCCCCC                                       &
     ! NNNN         NNN     EEEEEEEEE     MMMM           MMMM   '
    ! Write(ounit,*) '   CCCC      CCCCC                                     &
    ! NNNNNN       NNNN    EEEE          MMMMMM         MMMMMM  '
    ! Write(ounit,*) '  CCCC        CCCC                                     &
    ! NNNNNNN      NNNN    EEEE          MMMMMMM       MMMMMMM  '
    ! Write(ounit,*) ' CCCC           CC                                     &
    ! NNNN  NN     NNNN    EEEE          MMMM  MMM   MMM  MMMM  '
    ! Write(ounit,*) ' CCCC                                                  &
    ! NNNN  NNN    NNNN    EEEE          MMMM   MMMMMMM   MMMM  '
    ! Write(ounit,*) ' CCCC                    AAAAAAAAAA         RRRRRRRRR  &
    ! NNNN   NNN   NNNN    EEEEEEEEE     MMMM    MMMMM    MMMM  '
    ! Write(ounit,*) ' CCCC                   AAAAAAAAAAAA       RRRRRRRR    &
    ! NNNN   NNN   NNNN    EEEEEEEEEE    MMMM     MMM     MMMM  '
    ! Write(ounit,*) ' CCCC                  AAAA      AAAA     RRR          &
    ! NNNN    NNN  NNNN    EEEEEEEEE     MMMM             MMMM  '
    ! Write(ounit,*) ' CCCC           CC    AAAA       AAAA    RRRR          &
    ! NNNN     NN  NNNN    EEEE          MMMM             MMMM  '
    ! Write(ounit,*) '  CCCC        CCCC    AAAA       AAAA    RRRR          &
    ! NNNN      NNNNNNN    EEEE          MMMM             MMMM  '
    ! Write(ounit,*) '   CCCC      CCCCC     AAAA      AAAA    RRRR          &
    ! NNNN       NNNNNN    EEEE          MMMM             MMMM  '
    ! Write(ounit,*) '    CCCCCCCCCCCCC       AAAAAAAAAAAAA    RRRR          &
    ! NNN         NNNN      EEEEEEEEE    MMM               MMM  '
    ! Write(ounit,*) '     CCCCCCCCCC          AAAAAAAAAAAA    RRRR          &
    ! NN           NN        EEEEEEEEE   MM                 MM  '
    ! Write(ounit,*)
    ! Write(ounit,*) 'CarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEM' &
    ! 'CarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEMCarNEM'

End Subroutine

Subroutine title2()
       Write(*,*)''
       Write(*,*)'                         Finished                         '                             
       Write(*,*)'                         --------                         '  
End Subroutine
