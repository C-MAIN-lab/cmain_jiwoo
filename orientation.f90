program new_step

      Implicit none

      INTEGER ::  i, j,  &
                        atom_count, check, k, l,filter_count
      REAL,DIMENSION(:,:),ALLOCATABLE :: ATOM, at
      INTEGER,ALLOCATABLE :: SORT(:)
      REAL(16) :: xlength, ylength, zlength
      character(len=100) :: str, aa
      character(len=100) :: path, filename

      REAL :: save_syslength(3)
      INTEGER :: step, ct, cct, ctt, cttt, ctttt, countt, counttt,&
        st, seq, Max_step
      INTEGER,ALLOCATABLE :: Car_1(:), Nit_A(:), Nit_T (:), Car_33(:)
      Double precision,DIMENSION(:,:),ALLOCATABLE ::d,N_ACN_D,d1,d2,d3,d5
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: ACN
      Double precision,ALLOCATABLE :: x(:), y(:), z(:),&
                                        MIN_d(:)
      Double precision :: dx, dy, dz, ax, ay, az
      LOGICAL :: found
      REAL :: Angle,&
                x1, x2, x3, y1, y2, y3, z1, z2, z3, &
                ab(3), bc(3), an, vec, a, b, c, &
                x4, y4, z4, pi,&
                start_time, end_time,elapsed_time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call getcwd(path)
        print *, path
        filename = trim(path) // "/" // trim("dump")
        print *, filename
        open(1, file=filename, status='old', iostat=check)
        if (check /= 0) then
                print *, 'No File'
                STOP
        endif !check /=0
        step = 0
        REWIND(1)

        DO 
        READ(1,'(A)',IOSTAT=check)str
        IF (check /= 0 ) exit

                if (INDEX(str,'NUMBER OF ATOMS') > 0 ) then
                        READ(1,*) atom_count
                        cycle
                end if !NUMBER
        print *, atom_count
                if (INDEX(str, 'BOUNDS pp') > 0) then
                        READ(1,*) save_syslength(1:2)
                        xlength = save_syslength(2) - save_syslength(1)
                        READ(1,*) save_syslength(1:2)
                        ylength = save_syslength(2) - save_syslength(1)
                        READ(1,*) save_syslength(1:2)
                        zlength = save_syslength(2) - save_syslength(1)

                print *, xlength,ylength,zlength
                cycle
                end if !BOUNDS
                if (INDEX(str, 'id') > 0) then
                        step = step + 1
                ALLOCATE(ATOM(atom_count,5))
                ALLOCATE(at(atom_count, 5))
                ALLOCATE(SORT(atom_count))


           DO i = 1, atom_count
                READ(1,*) ATOM(i,1), &
          ATOM(i,2), ATOM(i,3), ATOM(i,4), ATOM(i,5)
           ENDDO !ATOM(i,5)
                DO l = 1, atom_count
                        SORT(l) = ATOM(l,1)
                        at(SORT(l),:) = ATOM(l,:)
                End Do ! sort
                Do k =1, atom_count
                at(k,3) = at(k,3)*xlength
                at(k,4) = at(k,4)*ylength
                at(k,5) = at(k,5)*zlength
                End Do  ! at * length
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        pi = 3.14159265358979323846
        call cpu_time(start_time)
        print *, atom_count

        Open(10, file='small.dat')
        Open(11, file='middle.dat')

        ALLOCATE(x(atom_count),y(atom_count),z(atom_count))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ct = 0
        ctt = 0        
        cttt= 0
        ctttt= 0
        Do i = 1, atom_count
      if (at(i,3) >= 83 .and. at(i,3) <= 107 .and.& 
         at(i,5) >= 43.0 .and. at(i,5) <= 60.0) then
            
                if (at(i,2) == 7) then
                        ctttt = ctttt + 1
                elseif (at(i,2) == 1) then
                        cttt = cttt + 1
                elseif (at(i,2) == 5) then
                        ctt = ctt + 1
                elseif (at(i,2) == 6) then
                        ct = ct + 1
                end if
        end if !distance limit
!        print *, cttt
        End do ! ctttt down
        
        ALLOCATE (Car_1(ct), Nit_A(ctt), Nit_T(cttt), Car_33(ctttt))
        ALLOCATE (ACN(ct,3))
        ALLOCATE (N_ACN_D(cttt,ct))
        ALLOCATE (d1(ct,ctttt))
        ALLOCATE (d2(ct,ctt))
        ALLOCATE (d3(ctt,ctttt))
        ALLOCATE (d5(cttt,ct))
        print *, ct,ctt,cttt,ctttt

        ct = 0
        ctt = 0
        cttt= 0
        ctttt= 0

        Do i = 1, atom_count
        if (at(i,3) >= 83 .and. at(i,3) <= 107 .and. &
         at(i,5) >= 43.0 .and. at(i,5) <= 60.0) then
                if (at(i,2) == 7) then
                        ctttt = ctttt +1
                        Car_33(ctttt) = at(i,1)
                end if
                if (at(i,2) == 1) then
                        cttt = cttt + 1
                        Nit_T(cttt)= at(i,1)
                end if
                if (at(i,2) == 5) then
                        ctt = ctt + 1
                        Nit_A(ctt)= at(i,1)
                end if
                if (at(i,2) == 6) then
                        ct = ct + 1
                        Car_1(ct)= at(i,1)
                end if
        endif !distance limit
        End Do ! ctttt down
        print *, ct,ctt,cttt,ctttt
        print *, 'Starting angle and distance calculations for step ', step
!----------------------------------------------------
!----------------START STEP LOOP---------------------
!----------------------------------------------------
               
                 Do i=1, ct 
                        Do j=1, ctttt
                                dx = at(Car_33(j), 3) - at(Car_1(i), 3)
                                dx = dx - nint(dx/xlength)*xlength
                                dy = at(Car_33(j), 4) - at(Car_1(i), 4)
                                dy = dy - nint(dy/ylength)*ylength
                                dz = at(Car_33(j), 5) - at(Car_1(i), 5)
                                dz = dz - nint(dz/zlength)*zlength

                                d1(i,j) = sqrt(dx**2 + dy**2 + dz**2)
                             EndDo !j
                      EndDo !i

                      DO i=1,ct
                        Do k=1,ctt
                                        dx = at(Nit_A(k), 3) -at(Car_1(i), 3)
                                        dx = dx -nint(dx/xlength)*xlength
                                        dy = at(Nit_A(k), 4) -at(Car_1(i), 4)
                                        dy = dy -nint(dy/ylength)*ylength
                                        dz = at(Nit_A(k), 5) -at(Car_1(i), 5)
                                        dz = dz -nint(dz/zlength)*zlength                     
                                        d2(i,k) = sqrt(dx**2 + dy**2 + dz**2)
                             EndDo !k
                           EndDo !i

                           Do k=1,ctt
                              Do j=1,ctttt
                                        dx = at(Nit_A(k), 3) -at(Car_33(j), 3)
                                        dx = dx -nint(dx/xlength)*xlength
                                        dy = at(Nit_A(k), 4) -at(Car_33(j), 4)
                                        dy = dy -nint(dy/ylength)*ylength
                                        dz = at(Nit_A(k), 5) -at(Car_33(j), 5)
                                        dz = dz-nint(dz/zlength)*zlength
                                        d3(k,j) = sqrt(dx**2 + dy**2 + dz**2)
                                EndDo !j
                                EndDo !k
                
                                print *, 'Finished distance'
               

                     
                DO i = 1, ct
                    ACN(i,1) = Car_1(i)
                    ACN(i,2) = 0
                    ACN(i,3) = 0
                    found = .false.
                    
                    DO j = 1, ctttt
                        if (found) exit
                        Do k = 1, ctt
                            if (d1(i,j) <1.46 .AND.d2(i,k) <1.16 &
                                    .AND. d3(k,j) <2.62) then
                                    
                                    ACN(i,2) = Car_33(j)
                                    ACN(i,3) = Nit_A(k)
                                    print *,d1(i,j),d2(i,k),d3(k,j), &
                                    Car_1(i),Car_33(j),Nit_A(k)
                                    found = .true.
                                    exit
                            end if ! distance ACN limit
                        End Do  !k
                    End Do  !j
                End Do !i
                
                print *, 'Total ACN found = ', ct
                
                DO i=1, cttt
                   Do j=1,ct
                        dx = at(Nit_T(i), 3) - at(ACN(j,1),3)
                        dx = dx - nint(dx/xlength)*xlength
                        dy = at(Nit_T(i), 4) - at(ACN(j,1),4)
                        dy = dy - nint(dy/ylength)*ylength
                        dz = at(Nit_T(i), 5) - at(ACN(j,1),5)
                        dz = dz - nint(dz/zlength)*zlength

                                d5(i,j) = sqrt(dx**2 + dy**2 + dz**2)
                        ENDdo
                        ENDdo

!---------------------------------------------------------------
                countt = 0
                counttt = 0
                Do i = 1, cttt
                !--------------Nit_T_Coordnation-----------
                        a = (at(Nit_T(i),3) - xlength/2)
                        b = (at(Nit_T(i),4) - ylength/2)
                        c = (at(Nit_T(i),5) - zlength/2)
                        
                        x1 = at(Nit_T(i),3) - a
                        y1 = at(Nit_T(i),4) - b
                        z1 = at(Nit_T(i),5) - c
                        Do j = 1, ct
                            if (ACN(j,2) == 0 .or. ACN(j,3) == 0) cycle
                       ! if ( d(Nit_T(i),Car_1(j)) <= 50 ) then
!-----------------------------------ACN_Car_1_coordination-----------------
                                        x2 = at(ACN(j,1),3) - a
                                                if ( x2 < 0 ) then
                                                        x2 = x2+xlength
                                                else if (x2 > xlength) then
                                                        x2 = x2-xlength
                                                end if
                                        y2 = at(ACN(j,1),4) - b
                                                if ( y2 < 0 ) then
                                                        y2 = y2+ylength
                                                else if (y2 > ylength) then
                                                        y2 = y2-ylength
                                                end if
                                        z2 = at(ACN(j,1),5) - c
                                                if ( z2 < 0 ) then
                                                        z2 = z2+zlength
                                                else if (z2 > zlength) then
                                                        z2 = z2-zlength
                                                end if
!-----------------------------------ACN_Nit_A_coordination-----------------
                                        x3 = at(ACN(j,3),3) - a
                                                if ( x3 < 0 ) then
                                                        x3 = x3+xlength
                                                else if (x3 > xlength) then
                                                        x3 = x3-xlength
                                                end if
                                        y3 = at(ACN(j,3),4) - b
                                                if ( y3 < 0 ) then
                                                        y3 = y3+ylength
                                                else if (y3 > ylength) then
                                                        y3 = y3-ylength
                                                end if
                                        z3 = at(ACN(j,3),5) - c
                                                if ( z3 < 0 ) then
                                                        z3 = z3+zlength
                                                else if (z3 > zlength) then
                                                        z3 = z3-zlength
                                                end if

!-----------------------------------ACN_Car_33_coordination-----------------
                                        x4 = at(ACN(j,2),3) - a
                                                if ( x4 < 0 ) then
                                                        x4 = x4+xlength
                                                else if (x4 > xlength) then
                                                        x4 = x4-xlength
                                                end if
                                        y4 = at(ACN(j,2),4) - b
                                                if ( y4 < 0 ) then
                                                        y4 = y4+ylength
                                                else if (y4 > ylength) then
                                                        y4 = y4-ylength
                                                end if
                                        z4 = at(ACN(j,2),5) - c
                                                if ( z4 < 0 ) then
                                                        z4 = z4+zlength
                                                else if (z4 > zlength) then
                                                        z4 = z4-zlength
                                                end if

                                        ab(1) = x2-x1
                                        ab(2) = y2-y1
                                        ab(3) = z2-z1
                                        bc(1) = x4-x3
                                        bc(2) = y4-y3
                                        bc(3) = z4-z3
                                        vec = &
                        (ab(1)*bc(1)+ab(2)*bc(2)+ab(3)*bc(3))/ &
                        (sqrt(ab(1)**2+ab(2)**2+ab(3)**2)*sqrt(bc(1)**2+bc(2)**2+bc(3)**2))
                                        an = vec
                                        Angle = ACOS(an) *(180/pi)
                                
                                if ( d5(i,j) <= 8.5) then
                                        write(10,*) Angle
                                elseif ( d5(i,j) <12) then
                                        write(11,*) Angle
                                End if  ! write                                  
                        End do !j
                End do !i
                print *, step


                Deallocate(ATOM,at,SORT,x,y,z,Car_1,Nit_A,Nit_T,Car_33,ACN,N_ACN_D,d1,d2,d3,d5)
                cycle
                End if
                End Do !Do

                close(1)
                print *, step
!--------------------------------------------------------------
!------------------END STEP LOOP-------------------------------
!--------------------------------------------------------------

        call cpu_time(end_time)
        elapsed_time = end_time - start_time
        print *, elapsed_time







      END Program new_step
