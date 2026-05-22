
      Program MAKE
      Implicit none

      INTEGER   :: j, C_count, H_count, O_count, Bond_count, max_bond,&
                        A_c, max_Angle, D_c, max_Dihe
      INTEGER,ALLOCATABLE :: C(:), H(:), O(:)
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: Bond, Angle, Dihe
      REAL      :: ax, ay, az, dx, dy, dz, xlength, ylength, zlength
      REAL,ALLOCATABLE :: x(:), y(:), z(:)
      REAL,DIMENSION(:,:),ALLOCATABLE :: d


      INCLUDE 'readdata.f90'
!---------------------------------------

        OPEN(14, file='Bonds.txt')
        OPEN(15, file='Angles.txt')



        ALLOCATE (x(atom_count), y(atom_count), z(atom_count))
        ALLOCATE (d(atom_count,atom_count))

        xlength = xhi - xlo
        ylength = yhi - ylo
        zlength = zhi - zlo

        Do i=1, atom_count
                        x(i) = ATOM(i,5)
                        y(i) = ATOM(i,6)
                        z(i) = ATOM(i,7)
        End Do
        Do i=1, atom_count
                Do j=1, atom_count
                                ax = ABS(x(i)-x(j))
                                ay = ABS(y(i)-y(j))
                                az = ABS(z(i)-z(j))
                                if ( ax > (xlength/2)) then
                                        dx = xlength - ax
                                else
                                        dx = ax
                                end if
                                if ( ay > (ylength/2)) then
                                        dy = ylength - ay
                                else
                                        dy = ay
                                end if
                                if ( az > (zlength/2)) then
                                        dz = zlength - az
                                else
                                        dz = az
                                end if
                                d(i,j)  = sqrt(dx**2 + dy**2 + dz**2)
                End Do
        End Do
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!----------------------------------------------------------------
!----------------------------------------------------------------
      
        C_count = 0
        Do i = 1, atom_count
                if (ATOM(i,3) == 3 .OR. ATOM(i,3) == 1) then
                        C_count = C_count + 1
                end if
        End Do

        ALLOCATE (O(O_count), H(H_count), C(C_count))

        O_count = 0
        H_count = 0
        C_count= 0
        Do i = 1, atom_count
                if (ATOM(i,3) == 3 .OR. ATOM(i,3) == 1) then
                        C_count = C_count + 1
                        C(C_count)= ATOM(i,1)
                end if
        End Do
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!----------------------------------------------------------------
!-----------------------------------------------------------------
!----------------------------------------------------------------
!----------------------------------------------------------------

        max_bond = atom_count*(atom_count-1)

        ALLOCATE (Bond(max_bond,4))
        Bond_count = 0


        Do i = 1, C_count-1
                Do j = i+1, C_count
                        if ( d(C(i),C(j)) < 1.5 ) then
                                Bond_count = Bond_count + 1
                                Bond(Bond_count,1) = Bond_count
                                Bond(Bond_count,2) = 2
                                if (C(i) < C(j)) then
                                        Bond(Bond_count,3) = C(i)
                                        Bond(Bond_count,4) = C(j)
                                else
                                        Bond(Bond_count,3) = C(j)
                                        Bond(Bond_count,4) = C(i)
                                end if
                        endif
                End Do
        End Do

!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!----------------------------------------------------------------
!----------------------------------------------------------------
!----------------------------------------------------------------
        i = 0
        print *, 
        print *, "Bonds"
        print *, 


15      i = i + 1
        if (Bond_count < i) then
                Go To 16
        End if
        write(14,*) Bond(i,1),Bond(i,2),Bond(i,3),Bond(i,4)
        print *, Bond(i,1),Bond(i,2),Bond(i,3),Bond(i,4)
        Go To 15
        
16      max_Angle = Bond_count * Bond_count
        ALLOCATE ( Angle(max_Angle,5) )
        
        A_c = 0

        Do i = 1, Bond_count - 1
                Do j = i+1, Bond_count
                        if (Bond(i,3) == Bond(j,3)) then
                                A_c = A_c + 1
                                Angle(A_c,1) = A_c
                                Angle(A_c,2) = 9
                                Angle(A_c,3) = Bond(i,4)
                                Angle(A_c,4) = Bond(i,3)
                                Angle(A_c,5) = Bond(j,4)
                        elseif (Bond(i,3) == Bond(j,4)) then
                                A_c = A_c + 1
                                Angle(A_c,1) = A_c
                                Angle(A_c,2) = 9
                                Angle(A_c,3) = Bond(i,4)
                                Angle(A_c,4) = Bond(i,3)
                                Angle(A_c,5) = Bond(j,3)
                        elseif (Bond(i,4) == Bond(j,3)) then
                                A_c = A_c + 1
                                Angle(A_c,1) = A_c
                                Angle(A_c,2) = 9
                                Angle(A_c,3) = Bond(i,3)
                                Angle(A_c,4) = Bond(i,4)
                                Angle(A_c,5) = Bond(j,4)
                        elseif (Bond(i,4) == Bond(j,4)) then
                                A_c = A_c + 1
                                Angle(A_c,1) = A_c
                                Angle(A_c,2) = 9
                                Angle(A_c,3) = Bond(i,3)
                                Angle(A_c,4) = Bond(i,4)
                                Angle(A_c,5) = Bond(j,3)
                        endif
                End Do
        End Do

        i = 0

        print *, 
        print *, "Angles"
        print *, 




25      i = i + 1
        if (A_c < i) then
                Go To 26
        End if

        write(15,*) Angle(i,1),Angle(i,2),Angle(i,3),Angle(i,4),Angle(i,5)

!       print *, Angle(i,1),Angle(i,2),Angle(i,3),Angle(i,4),Angle(i,5)
        Go To 25

26      max_Dihe = A_c*(A_c-1)
        ALLOCATE (Dihe(max_Dihe,6))
        D_c = 0
        Do  i = 1, A_c-1
                Do j = i+1, A_c
                        if (Angle(i,3) == Angle(j,4) .AND. &
                                        Angle(i,4) == Angle(j,5)) then
                                D_c = D_c + 1
                                Dihe(D_c,1) = D_c
                                Dihe(D_c,2) = 9
                                Dihe(D_c,3) = Angle(j,3)
                                Dihe(D_c,4) = Angle(i,3)
                                Dihe(D_c,5) = Angle(i,4)
                                Dihe(D_c,6) = Angle(i,5)
                        else if (Angle(i,3) == Angle(j,4) .AND. &
                                        Angle(i,4) == Angle(j,3)) then
                                D_c = D_c + 1
                                Dihe(D_c,1) = D_c
                                Dihe(D_c,2) = 9
                                Dihe(D_c,3) = Angle(i,5)
                                Dihe(D_c,4) = Angle(j,3)
                                Dihe(D_c,5) = Angle(j,4)
                                Dihe(D_c,6) = Angle(j,5)
                        else if (Angle(i,4) == Angle(j,3) .AND. &
                                        Angle(i,5) == Angle(j,4)) then
                                D_c = D_c + 1
                                Dihe(D_c,1) = D_c
                                Dihe(D_c,2) = 9
                                Dihe(D_c,3) = Angle(i,3)
                                Dihe(D_c,4) = Angle(i,4)
                                Dihe(D_c,5) = Angle(i,5)
                                Dihe(D_c,6) = Angle(j,5)
                        else if (Angle(i,4) == Angle(j,5) .AND. &
                                        Angle(i,5) == Angle(j,4)) then
                                D_c = D_c + 1
                                Dihe(D_c,1) = D_c
                                Dihe(D_c,2) = 9
                                Dihe(D_c,3) = Angle(i,3)
                                Dihe(D_c,4) = Angle(i,4)
                                Dihe(D_c,5) = Angle(i,5)
                                Dihe(D_c,6) = Angle(j,3)
                        End if
                End Do
        End Do

        i = 0

        print *, 
        print *, "Dihedrals"
        print *, 


35      i = i + 1
        if (D_c < i) then
                Go To 36
        End if
!        print *, &
!        Dihe(i,1),Dihe(i,2),Dihe(i,3),Dihe(i,4),Dihe(i,5),Dihe(i,6)
        Go To 35

        



36      End Program MAKE

