!-*- F90 -*------------------------------------------------------------
! 
!  module: out_calc / meta3
!
!  output support functions and subroutines 
!
!----------------------------------------------------------------------


!======================================================================
!

module out_calc

  use constant
  use buflist

  implicit none
  private
  save

  ! --- Module Identifier

  character(len=20), private, parameter :: modname = 'OUT'

  ! --- Public Methods

  public :: PaScalar, PaVector, PaBuffer 

  ! --- Public Data


contains

!----------------------------------------------------------------------


!----------------------------------------------------------------------

    subroutine PaScalar(pa,fc,vc)

      integer :: pa
      M4_FTYPE :: fc
      real(kind=8) :: vc

      M4_IFELSE_CF({vc = real(fc)},{vc = fc})
      if ( pa .eq. 1 ) then
         M4_IFELSE_CF({
         vc = abs(fc)
         })
      else
         if ( pa .eq. 2 ) then
            M4_IFELSE_CF({
            vc = atan2(aimag(fc),real(fc))
            },{
            vc = 0.0 
            })
         endif
      endif

    end subroutine PaScalar

!----------------------------------------------------------------------

    subroutine PaVector(pa,fx,fy,fz,vx,vy,vz)

      integer :: pa
      M4_FTYPE :: fx,fy,fz
      real(kind=8) :: vx,vy,vz

      vx = real(fx)
      vy = real(fy)
      vz = real(fz)
      if ( pa .eq. 1 ) then
         M4_IFELSE_CF({
         vx = abs(fx)
         vy = abs(fy)
         vz = abs(fz)
         })
      else
         if ( pa .eq. 2 ) then
            M4_IFELSE_CF({
            vx = atan2(aimag(fx),real(fx))
            vy = atan2(aimag(fy),real(fy))
            vz = atan2(aimag(fz),real(fz))
            },{
            vx = 0.0 
            vy = 0.0 
            vz = 0.0 
            })
         endif
      endif
      
    end subroutine PaVector

!----------------------------------------------------------------------

    subroutine PaBuffer(pa, buf, p, val)

      integer :: pa
      type (T_BUF) :: buf
      integer :: p, i
      real(kind=8) :: val(:)
      
      val = buf%data(p,:)
      if ( pa .eq. 1 ) then
         M4_IFELSE_CF({
           do i = 1, buf%numslot
              val(i) = abs(buf%cdata(p,i))
           end do
         },{
           do i = 1, buf%numslot
              val(i) = abs(buf%data(p,i))
           end do
         })
      else
         if ( pa .eq. 2 ) then
            M4_IFELSE_CF({
             do i = 1, buf%numslot
                val(i) = atan2(aimag(buf%cdata(p,i)),real(buf%cdata(p,i)))
             end do
            },{
               val = 0.0
            })
         endif
      endif
      
    end subroutine PaBuffer

  
!----------------------------------------------------------------------

end module out_calc

!
! Authors:  J.Hamm 
! Modified: 6/1/2008
!
!======================================================================
