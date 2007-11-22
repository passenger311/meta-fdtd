!----------------------------------------------------------------------
!
!  module: config / max3d
!
!  load configuration files and inject data to other modules.
!
!----------------------------------------------------------------------


module config

  use grid
  use upml
  implicit none
  save

  character(len=255) :: file

contains

  subroutine ReadConfig(insfx)

    implicit none

    character(len=*) :: insfx

    call ReadGrid(insfx)
!    call ReadOutput
!    call ReadTFSF

  end subroutine ReadConfig


! ---- ReadGrid

  subroutine ReadGrid(sfx)

    implicit none

    character(len=*) :: sfx

    integer :: err, i


    file = cat2(pfxgrid,sfx)

    open(UNITTMP,FILE=file,STATUS='unknown')
    read(UNITTMP,*) PARTITIONS
    read(UNITTMP,*) NCYCMAX
    read(UNITTMP,*) DT
    read(UNITTMP,*) IBEG, IEND    ! from ... to ranges
    read(UNITTMP,*) JBEG, JEND
    read(UNITTMP,*) KBEG, KEND
    read(UNITTMP,*) pmlpart
    read(UNITTMP,*) (planepml(i),i=1, 6)
    read(UNITTMP,*) PMLMAX
    read(UNITTMP,*) PotPml
    read(UNITTMP,*) SigmaMax
    read(UNITTMP,*) KappaMax 

    close(UNITTMP) 

    write(6,*) "PotPml = ", PotPml

    if ( PARTITIONS .ne. numproc .and. mpi_started .ne. 0 ) then
       write(6,*) "config: number of read in parititons does not match mpi numproc"
       stop
    end if

    ! Inner ranges are [IBEG,IEND] etc., outer ranges [IMIN,IMAX]

    IMIN = IBEG-1          
    IMAX = IEND+1
    JMIN = JBEG-1
    JMAX = JEND+1
    KMIN = KBEG-1
    KMAX = KEND+1

   ! Ranges for PML sheaths [ISIG,IEIG]

    ISIG=IBEG
    IEIG=IMAX
    JSIG=IBEG
    JEIG=JMAX
    KSIG=KBEG
    KEIG=KMAX

    if(planepml(1) .eq. 1) ISIG=IBEG+PMLMAX
    if(planepml(2) .eq. 1) IEIG=IMAX-PMLMAX
    if(planepml(3) .eq. 1) JSIG=JBEG+PMLMAX
    if(planepml(4) .eq. 1) JEIG=JMAX-PMLMAX
    if(planepml(5) .eq. 1) KSIG=KBEG+PMLMAX
    if(planepml(6) .eq. 1) KEIG=KMAX-PMLMAX

    if((IEIG.le.IBEG) .or. (JEIG.le.JBEG) .or. (KEIG.le.KBEG)) then
       write(STDERR,*) "Error with IEIG, JEIG or KEIG in InitPml()"
       stop
    endif
   
  end subroutine ReadGrid


end module config







