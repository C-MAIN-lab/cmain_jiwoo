!
!        Program read_data
!        Implicit none

        INTEGER :: atom_count, atom_type, check, k, i
        REAL*8   :: xlo, xhi, ylo, yhi, zlo, zhi
        REAL,ALLOCATABLE,DIMENSION(:,:) :: ATOM
        character(len=100) :: path, filename, aa


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call getcwd(path)        
        filename = trim(path) // "/" // trim("C_data")
        print *, filename
        open(1, file=filename, status='old', iostat=check)
        if (check == 0) then
                Go to 22
        else
                print *, 'No File'
                STOP
        endif

!       OPEN(10, file='asd')

22      READ(1,*)
        READ(1,*) atom_count
        READ(1,*) atom_type
        READ(1,*) xlo, xhi
        READ(1,*) ylo, yhi
        READ(1,*) zlo, zhi
        
        ALLOCATE (ATOM(atom_count,7))
!        Print *, atom_count, atom_type, xlo, xhi

10      READ(1,*) aa
        if (aa(1:5) == 'Atoms') then
                GO To 11
        endif
        GO To 10
11      READ(1,*)
        k = 0
55      k = k + 1
        READ(1,*,END=888) ATOM(k,1), &
          ATOM(k,2), ATOM(k,3), ATOM(k,4), ATOM(k,5), ATOM(k,6),&
                 ATOM(k,7)
!        Print *, ATOM(k,1), &
!          ATOM(k,2), ATOM(k,3), ATOM(k,4), ATOM(k,5), ATOM(k,6),&
!                 ATOM(k,7)
        GO To 55

888     close(1)
!        end program read_data
