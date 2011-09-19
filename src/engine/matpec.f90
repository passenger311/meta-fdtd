!-*- F90 -*------------------------------------------------------------
!
!  module: matpec / meta
!
!  Pec material module.
!
!  subs:
!
!    InitializeMatPec
!    FinalizeMatPec
!    ReadMatPecObj
!    StepEMatPec
!    StepHMatPec
!    SumJEMatPec
!    SumKHMatPec
!
!----------------------------------------------------------------------


! =====================================================================
!
! The MatPec module sets the electric field components to 0. within
! a selected region.
!

module matpec

  use constant
  use parse
  use reglist
  use outlist
  use grid
  use fdtd

  implicit none
  private
  save

  M4_MATHEAD_DECL({MATPEC},100,{

  ! no parameters

  })

contains

!----------------------------------------------------------------------

  subroutine ReadMatPecObj(funit,lcount)

    M4_MODREAD_DECL({MATPEC}, funit,lcount,mat,reg,out)

    M4_WRITE_DBG(". enter ReadMatPecObj")
    
    M4_MODREAD_EXPR({MATPEC},funit,lcount,mat,reg,3,out,{ 

    ! no parameters to read

    })

    call CompressValRegObj(reg) ! compress filling factors

    M4_WRITE_DBG(". exit ReadMatPecObj")

  end subroutine ReadMatPecObj

!----------------------------------------------------------------------

  subroutine InitializeMatPec

    M4_WRITE_DBG(". enter InitializeMatPec")
    M4_WRITE_DBG(". exit InitializeMatPec")

  end subroutine InitializeMatPec

!----------------------------------------------------------------------

  subroutine FinalizeMatPec

    M4_WRITE_DBG(". enter FinalizeMatPec")
    M4_WRITE_DBG(". exit FinalizeMatPec")

  end subroutine FinalizeMatPec

!----------------------------------------------------------------------

  subroutine StepHMatPec(ncyc)

    integer :: ncyc

    return
  
  end subroutine StepHMatPec


!----------------------------------------------------------------------

  subroutine StepEMatPec(ncyc)

    integer :: ncyc, m, n
    M4_MODLOOP_DECL({MATPEC},mat)
    M4_REGLOOP_DECL(reg,p,i,j,k,w(3))

    M4_MODLOOP_EXPR({MATPEC},mat,{

       ! this loops over all mat structures, setting mat

    M4_MODOBJ_GETREG(mat,reg)

       M4_REGLOOP_EXPR(reg,p,i,j,k,w,{  

       if ( w(1) .ne. 0. ) Ex(i,j,k) = 0.
       if ( w(2) .ne. 0. ) Ey(i,j,k) = 0.
       if ( w(3) .ne. 0. ) Ez(i,j,k) = 0.

       })      

    })

  end subroutine StepEMatPec

!----------------------------------------------------------------------

  subroutine SumJEMatPec(mask, ncyc, sum, idx, mode)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    logical :: mode
    integer :: ncyc, idx
    real(kind=8) :: sum(MAXEBALCH)
 
  end subroutine SumJEMatPec

!----------------------------------------------------------------------

  subroutine SumKHMatPec(mask, ncyc, sum, idx, mode)

    logical, dimension(IMIN:IMAX,JMIN:JMAX,KMIN:KMAX) :: mask
    integer :: ncyc, idx
    logical :: mode
    real(kind=8) :: sum(MAXEBALCH)

  end subroutine SumKHMatPec
 
!----------------------------------------------------------------------

  subroutine DisplayMatPecObj(mat)

    type(T_MATPEC) :: mat
 
    M4_WRITE_INFO({"#",TRIM(i2str(mat%idx)),&
    	"pec"
    })
    call DisplayRegObj(regobj(mat%regidx))
    	
  end subroutine DisplayMatPecObj

!----------------------------------------------------------------------

   subroutine EchoMatPecObj(mat)

    type(T_MATPEC) :: mat

    M4_WRITE_INFO({"--- matpec # ",&
         TRIM(i2str(mat%idx))," ", TRIM(mat%type)})
    M4_WRITE_INFO({"defined over:"})
    call EchoRegObj(regobj(mat%regidx))
    

  end subroutine EchoMatPecObj

  
!----------------------------------------------------------------------

end module matpec

! Authors:  J.Hamm 
! Modified: 15/02/2008
! Changed: 7/07/2011 S.Wuestner
!
! =====================================================================


