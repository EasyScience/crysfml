!!----
!!---- Menu: 3
!!---- Atoms Calculations
!!----
!!
 Module Menu_3
    !---- Use File ----!
    use Menu_0
    use CFML_Crystallographic_Symmetry
    use CFML_Atom_TypeDef
    use CFML_IO_Messages,      only: Wait_Message
    use CFML_string_utilities, only: getnum
    use CFML_Math_General,     only: cosd, acosd, sind, asind
    use CFML_Crystal_Metrics,  only: Crystal_Cell_Type, Write_Crystal_Cell, Cart_Vector
    use CFML_IO_Formats,       only: Readn_set_Xtal_Structure,err_form_mess,err_form,file_list_type
    use Menu_1,                only: Get_Wyckoff
    !---- Variables ----!
    implicit none

    logical                  :: structure_read=.false.
    type (space_group_type)  :: SpG
    type (Atom_list_Type)    :: A
    type (Crystal_Cell_Type) :: Cell
    real(kind=cp), parameter :: eps=0.00001
    real(kind=cp)            :: qcellp=0.0 , qcellm=0.0
    logical                  :: neutral=.true.

 Contains

    !!----
    !!---- Subroutine Menu_Princ3
    !!----
    !!
    Subroutine Menu_Princ3()
       !---- Local Variables ----!
       character(len=2)  :: car

       do
          call execute_command_line(clear_string)

          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "     GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "     Atomistic Calculations "
          write(unit=*,fmt="(a)") " =============================="
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " [0] Back..."
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " [1] Multiplicity of an atom position"
          write(unit=*,fmt="(a)") " [2] Read a CFL or CIF file to load a crystal structure"
          write(unit=*,fmt="(a)") " [3] Show current structure information"
          write(unit=*,fmt="(a)") " [4] Calculate Ionic Dipolar moment & polarisation of a symmetrized single unit cell"
          write(unit=*,fmt="(a)") " [5] Calculate Ionic polarisation with respect to a non-polar structure"
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)",advance="no") " OPTION: "
          read(*,'(a)') car
          if (len_trim(car) == 0) exit
          car=adjustl(car)

          select case (car)
             case ('0 ')
                exit

             case ('1 ')
                call Menu_ATOM_1()   !Multiplicity

             case ('2 ')
                call Menu_ATOM_2()   !Reading a CFL or CIF file

             case ('3 ')
                call Menu_ATOM_3()   !Showing the crystal structure

             case ('4 ')
                call Menu_ATOM_4()   !Calculating polarisation if charges are given

             case ('5 ')
                call Menu_ATOM_5()   !Calculating polarisation if charges are given

          end select
       end do

    End Subroutine Menu_Princ3

    !!----
    !!---- Subroutine Menu_ATOM_1
    !!----
    !!
    Subroutine Menu_Atom_1()

       !---- Local Variables ----!
       character(len=30)          :: line
       integer                    :: iv, mlt
       integer, dimension(3)      :: ivet
       real(kind=cp)              :: occ
       real(kind=cp), dimension(3):: vet
       type (Space_Group_type)    :: grp_espacial

       do
          call execute_command_line(clear_string)
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "     GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "     Space Groups Information "
          write(unit=*,fmt="(a)") " ================================"
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)",advance="no") " Space Group (HM/Hall/Num): "

          read(unit=*,fmt='(a)') line
          if (len_trim(line)==0) exit
          line=adjustl(line)
          call set_spacegroup(line,grp_espacial)
          if(Err_Symm) then
            write(unit=*,fmt="(a)") trim(ERR_Symm_Mess)
            call Wait_Message(" => Press <enter> to continue ...")
            cycle
          end if
          do
             call execute_command_line(clear_string)
             write(unit=*,fmt="(a)") " "
             write(unit=*,fmt="(a)") "       GENERAL CRYSTALLOGRAPHY CALCULATOR "
             write(unit=*,fmt="(a)") " "
             write(unit=*,fmt="(a)") "     Multiplicity and Occupancy of Position "
             write(unit=*,fmt="(a)") "   ==========================================="
             write(unit=*,fmt="(a)") " "
             write(unit=*,fmt="(a)") " "
             write(unit=*,fmt="(a)",advance="no") " Position: "

             read(*,'(a)') line
             if (len_trim(line)==0) exit
             write(unit=i_out,fmt="(a)") " "
             write(unit=i_out,fmt="(a)") "     Multiplicity and Occupancy of Position "
             write(unit=i_out,fmt="(a)") "   ==========================================="
             write(unit=i_out,fmt="(a)") " "
             line=adjustl(line)
             vet=0.0
             call getnum(line,vet,ivet,iv)
             if (iv == 3) then
                mlt=Get_Multip_Pos(vet,grp_espacial)
                occ=real(mlt,kind=cp)/real(grp_espacial%Multip,kind=cp)
                write(unit=*,fmt=*) " "
                write(unit=*,fmt='(a,i4,a,f7.4)') " Multiplicity: ",mlt, "     Occupancy(SHELX/FullProf) proportional to: ",occ
                call Wait_Message(" => Press <enter> to continue ...")
             end if
          end do
       end do

    End Subroutine Menu_Atom_1
    !!----
    !!---- Subroutine Menu_ATOM_2  Reading a CFL or CIF file
    !!----
    !!
    Subroutine Menu_Atom_2()
       !---- Local Variables ----!
       character(len=256)    :: line
       integer               :: i
       logical               :: esta

       do
          call execute_command_line(clear_string)
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "     GENERAL CRYSTALLOGRAPHY CALCULATOR "
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)") "     Atomistic Calculations "
          write(unit=*,fmt="(a)") " ================================"
          write(unit=*,fmt="(a)") " Reading a CFL or CIF file ..."
          write(unit=*,fmt="(a)") " "
          write(unit=*,fmt="(a)",advance="no") " => Enter the full name of the file: "

          read(unit=*,fmt='(a)') line
          if (len_trim(line)==0) exit
          write(unit=i_out,fmt="(a)") " "
          write(unit=i_out,fmt="(a)") "     Atomistic Calculations "
          write(unit=i_out,fmt="(a)") " ================================"
          write(unit=i_out,fmt="(a)") " Reading a CFL or CIF file ..."
          write(unit=i_out,fmt="(a)") " "

          line=adjustl(line)
          i=index(line,".cif")
          if(i /= 0) then
             inquire(file=trim(line),exist=esta)
             if (esta) then
                call Readn_set_Xtal_Structure(trim(line),Cell,SpG,A,Mode="CIF")
             else
                write(unit=*,fmt="(a)") " => File: "//trim(line)//"  does not exist!"
                call Wait_Message(" => Press <enter> to continue ...")
                cycle
             end if
          else
            i=index(line,".cfl")
            if(i /= 0) then
               inquire(file=trim(line),exist=esta)
               if (esta) then
                  call Readn_set_Xtal_Structure(trim(line),Cell,SpG,A,Mode="CFL")
               else
                  write(unit=*,fmt="(a)") " => File: "//trim(line)//"  does not exist!"
                  call Wait_Message(" => Press <enter> to continue ...")
                  cycle
               end if
            else
               write(unit=*,fmt="(a)") " => Illegal File name: "//trim(line)
               call Wait_Message(" => Press <enter> to continue ...")
               cycle
            end if
          end if
          if (err_form) then
             write(unit=*,fmt="(a)") trim(err_form_mess)
          else
             write(unit=*,fmt="(a)") " => File: "//trim(line)//" successfully read!"
             structure_read=.true.
             qcellp=0.0
             qcellm=0.0
             do i=1,A%natoms
               if(A%atom(i)%charge > 0.01)  then
                 qcellp=qcellp+A%atom(i)%Mult*A%atom(i)%charge
               else
                 qcellm=qcellm+A%atom(i)%Mult*abs(A%atom(i)%charge)
               end if
             end do
             if(abs(qcellp-qcellm) > eps) neutral=.false.
             call Wait_Message(" => Press <enter> to continue ...")
             exit
          end if

       end do

    End Subroutine Menu_Atom_2

    !!----
    !!---- Subroutine Menu_ATOM_3  Showing Cell, Space Group and Atom positions
    !!----
    !!
    Subroutine Menu_Atom_3()

       call execute_command_line(clear_string)
       write(unit=*,fmt="(a)") " "
       write(unit=*,fmt="(a)") "     GENERAL CRYSTALLOGRAPHY CALCULATOR "
       write(unit=*,fmt="(a)") " "
       write(unit=*,fmt="(a)") "     Atomistic Calculations "
       write(unit=*,fmt="(a)") " ================================"
       write(unit=*,fmt="(a)") " Showing Cell, Space Group and Atom positions ..."
       write(unit=*,fmt="(a)") " "
       if(.not. structure_read) then
         write(unit=*,fmt="(/a/)") " => A CIF or CFL file must read before using this option! "
         call Wait_Message(" => Press <enter> to continue ...")
         return
       end if
       write(unit=i_out,fmt="(a)") " "
       write(unit=i_out,fmt="(a)") "     Atomistic Calculations "
       write(unit=i_out,fmt="(a)") " ================================"
       write(unit=i_out,fmt="(a)") " Showing Cell, Space Group and Atom positions ..."
       write(unit=i_out,fmt="(a)") " "
       write(unit=i_out,fmt="(/a/)") " => UNIT CELL information: "
       write(unit=*,fmt="(/a/)") " => UNIT CELL information: "

       call Write_Crystal_Cell(Cell)
       call Write_Crystal_Cell(Cell,i_out)

       call Wait_Message(" => Press <enter> to continue ...")

       call execute_command_line(clear_string)
       write(unit=*,fmt="(/a/)") " => SPACE GROUP information: "
       call Write_SpaceGroup(SpG,Full=.true.)
       call Write_SpaceGroup(SpG,i_out,Full=.true.)
       call Wait_Message(" => Press <enter> to continue ...")

       call execute_command_line(clear_string)
       write(unit=*,fmt="(/a/)") " => ATOMS information: "
       call Write_Atom_List(A,level=1)
       call Write_Atom_List(A,level=1,lun=i_out)
       write(unit=*,fmt="(/a,f12.4)") " => Total positive charge per unit cell: ",qcellp
       write(unit=*,fmt="(a,f12.4/)") " => Total negative charge per unit cell: ",-qcellm
       write(unit=i_out,fmt="(/a,f12.4)") " => Total positive charge per unit cell: ",qcellp
       write(unit=i_out,fmt="(a,f12.4/)") " => Total negative charge per unit cell: ",-qcellm
       call Wait_Message(" => Press <enter> to continue ...")

    End Subroutine Menu_Atom_3

    Subroutine Menu_Atom_4()
       !---- Local Variables ----!
       integer                        :: i,j,Mult,m,lam,lbm,lcm,np,nm
       real(kind=cp)                  :: q, qp,pol,ang,ncell,dist,polc
       real(kind=cp), dimension(3)    :: pos,cpos,r_frac, r_pol,r_plus,r_minus,dir_pol
       real(kind=cp), dimension(3,384):: orb !increased to take into account the surface atoms
       real(kind=cp), dimension(  384):: qat !New charges
       logical                        :: calc_possible=.true.

       call execute_command_line(clear_string)
       write(unit=*,fmt="(a)") " "
       write(unit=*,fmt="(a)") "               GENERAL CRYSTALLOGRAPHY CALCULATOR "
       write(unit=*,fmt="(a)") " "
       write(unit=*,fmt="(a)") "                     Atomistic Calculations "
       write(unit=*,fmt="(a)") " ======================================================================="
       write(unit=*,fmt="(a)") " Calculating the ionic polarisation of a single symmetrized unit cell..."
       write(unit=*,fmt="(a)") " "
       do i=1,A%natoms
         if(abs(A%atom(i)%charge) <= 0.001)  then
           calc_possible=.false.
           exit
         end if
       end do
       if( .not. calc_possible) then
         write(unit=*,fmt="(a)") " => Calculation of P impossible. No charges have been provided! "
         call Wait_Message(" => Press <enter> to continue ...")
         return
       end if
       !if(SpG%Centred /= 1) then
       !  write(unit=*,fmt="(a)") " => Polarisation = 0.0 Centrosymmetric Crystal! "
       !  call Wait_Message(" => Press <enter> to continue ...")
       !  return
       !end if
      ! write(unit=*,fmt="(a)",advance="no")  " => Enter supercell extend (la,lb,lcm): "
      ! read(unit=*,fmt=*) lam,lbm,lcm
       lam=0; lbm=0; lcm=0
       ncell=(2*lam+1)*(2*lbm+1)*(2*lcm+1)
       !write(unit=*,fmt="(a,i10,a)")  " => Calculation for: ",nint(ncell)," unit cells"
       write(unit=i_out,fmt="(a)") " "
       write(unit=i_out,fmt="(a)") "                     Atomistic Calculations "
       write(unit=i_out,fmt="(a)") " ======================================================================="
       write(unit=i_out,fmt="(a)") " Calculating the ionic polarisation of a single symmetrized unit cell..."
       write(unit=i_out,fmt="(a)") " "
       r_plus=0.0; r_minus=0.0 ; np=0; nm=0
       do i=1,A%natoms
          pos=A%atom(i)%x
          q=A%atom(i)%charge
          call Get_Orbit(pos,SpG,Mult,orb) !here the orbit has always positive coordinates
          !Treatment of atoms at origin, edges and faces
          m=0; qat=0.0
          do j=1,Mult
            !Testing 0,0,0
            if(abs(orb(1,j)) < eps .and. abs(orb(2,j)) < eps  .and. abs(orb(3,j)) < eps  ) Then !atom at the origin
              qp=0.125*q
              qat(j)=qp
              m=m+1
              orb(:,Mult+m) = [1.0_cp,0.0_cp,0.0_cp] ; qat(Mult+m)=qp
              m=m+1
              orb(:,Mult+m) = [0.0_cp,1.0_cp,0.0_cp] ; qat(Mult+m)=qp
              m=m+1
              orb(:,Mult+m) = [0.0_cp,0.0_cp,1.0_cp] ; qat(Mult+m)=qp
              m=m+1
              orb(:,Mult+m) = [1.0_cp,1.0_cp,0.0_cp] ; qat(Mult+m)=qp
              m=m+1
              orb(:,Mult+m) = [1.0_cp,0.0_cp,1.0_cp] ; qat(Mult+m)=qp
              m=m+1
              orb(:,Mult+m) = [0.0_cp,1.0_cp,1.0_cp] ; qat(Mult+m)=qp
              m=m+1
              orb(:,Mult+m) = [1.0_cp,1.0_cp,1.0_cp] ; qat(Mult+m)=qp
            else if(abs(orb(1,j)) < eps .and. abs(orb(2,j)) < eps  ) Then !atom in z-edge
              qp=0.25*q
              qat(j)=qp
              m=m+1
              orb(:,Mult+m) = orb(:,j) + [1.0_cp,0.0_cp,0.0_cp]  ; qat(Mult+m)=qp
              m=m+1
              orb(:,Mult+m) = orb(:,j) + [0.0_cp,1.0_cp,0.0_cp]  ; qat(Mult+m)=qp
              m=m+1
              orb(:,Mult+m) = orb(:,j) + [1.0_cp,1.0_cp,0.0_cp]  ; qat(Mult+m)=qp
            else if(abs(orb(1,j)) < eps .and. abs(orb(3,j)) < eps  ) Then !atom in y-edge
              qp=0.25*q
              qat(j)=qp
              m=m+1
              orb(:,Mult+m) = orb(:,j) + [1.0_cp,0.0_cp,0.0_cp]  ; qat(Mult+m)=qp
              m=m+1
              orb(:,Mult+m) = orb(:,j) + [0.0_cp,0.0_cp,1.0_cp]  ; qat(Mult+m)=qp
              m=m+1
              orb(:,Mult+m) = orb(:,j) + [1.0_cp,0.0_cp,1.0_cp]  ; qat(Mult+m)=qp
            else if(abs(orb(2,j)) < eps .and. abs(orb(3,j)) < eps ) Then !atom in x-edge
              qp=0.25*q
              qat(j)=qp
              m=m+1
              orb(:,Mult+m) = orb(:,j) + [0.0_cp,1.0_cp,0.0_cp]  ; qat(Mult+m)=qp
              m=m+1
              orb(:,Mult+m) = orb(:,j) + [0.0_cp,0.0_cp,1.0_cp]  ; qat(Mult+m)=qp
              m=m+1
              orb(:,Mult+m) = orb(:,j) + [0.0_cp,1.0_cp,1.0_cp]  ; qat(Mult+m)=qp
            else if(abs(orb(1,j)) < eps ) Then !atom in yz-plane
              qp=0.5*q
              qat(j)=qp
              m=m+1
              orb(:,Mult+m) = orb(:,j) + [1.0_cp,0.0_cp,0.0_cp]  ; qat(Mult+m)=qp
            else if(abs(orb(2,j)) < eps ) Then !atom in xz-plane
              qp=0.5*q
              qat(j)=qp
              m=m+1
              orb(:,Mult+m) = orb(:,j) + [0.0_cp,1.0_cp,0.0_cp]  ; qat(Mult+m)=qp
            else if(abs(orb(3,j)) < eps ) Then !atom in xy-plane
              qp=0.5*q
              qat(j)=qp
              m=m+1
              orb(:,Mult+m) = orb(:,j) + [0.0_cp,0.0_cp,1.0_cp]  ; qat(Mult+m)=qp
            else
              qat(j)=q
            end if
          end do
          Mult=Mult+m
          write(unit=*,fmt="(a,a15,f8.2,a,i3)") " => Atom: ",trim(A%Atom(i)%Lab)//",  Charge: ", A%Atom(i)%Charge, &
                                                "     Modified Multiplicity:",mult
          write(unit=i_out,fmt="(a,a15,f8.2,a,i3)") " => Atom: ",trim(A%Atom(i)%Lab)//",  Charge: ", A%Atom(i)%Charge, &
                                                "     Modified Multiplicity:",mult
          cpos=0.0
          do j=1,Mult
            write(unit=*,fmt="(a,3f9.5,a,f9.5)") " => (x,y,z): ",orb(:,j),"  Qeff=",qat(j)
            write(unit=i_out,fmt="(a,3f9.5,a,f9.5)") " => (x,y,z): ",orb(:,j),"  Qeff=",qat(j)
            !do la=-lam,lam
            !  do lb=-lbm,lbm
            !    do lc=-lcm,lcm
                  pos=orb(:,j)-[0.5_cp,0.5_cp,0.5_cp]   !+[la,lb,lc]
                  r_frac=r_frac+ pos*qat(j)
                  cpos=cpos+pos
                  if(qat(j) > eps) then
                    np=np+1
                    r_plus=r_plus+orb(:,j)
                  else
                    nm=nm+1
                    r_minus=r_minus+orb(:,j)
                  end if
            !    end do
            !  end do
            !end do
          end do
         cpos=cpos/(ncell * real(mult,kind=cp))
         write(unit=*,fmt="(a,3f10.5,a)") "    Average of vector positions w.r.to (1/2,1/2,1/2): (",cpos," )"
         write(unit=i_out,fmt="(a,3f10.5,a)") "    Average of vector positions w.r.to (1/2,1/2,1/2): (",cpos," )"
       end do
       r_plus=r_plus/np
       r_minus=r_minus/nm
       !r_frac=r_frac/real(ncell)
       r_frac=(r_plus-r_minus) * qcellp
       r_pol=Cart_Vector("D",r_frac,Cell)
       where(abs(r_pol) <= 0.0001) r_pol=0.0
       pol=sqrt(dot_product(r_pol,r_pol))
       dir_pol=[0.0,0.0,1.0]
       if(pol > eps) dir_pol=r_pol/pol
       cpos=Cart_Vector("D",r_plus-r_minus,Cell)
       dist=sqrt(dot_product(cpos,cpos))
       polc=1602.176565*pol/Cell%CellVol
       write(unit=*,fmt="(/,a,3f8.4,a)")" => Geometric centre of Q+ (fract. coordinates)       :(",r_plus," )"
       write(unit=*,fmt="(a,3f8.4,a)")  " => Geometric centre of Q- (fract. coordinates)       :(",r_minus," )"
       write(unit=*,fmt="(a,f8.4)")     " => Distance between  Q+  and  Q- (Angstrom)          : ",dist
       write(unit=*,fmt="(a,g11.4)")    " => Ionic Dipolar Moment/UnitCell  (Coulomb.Metre)    :  ",1.602176565e-29*pol
       write(unit=*,fmt="(a,f8.4)")     " => Ionic Dipolar Moment/UnitCell  (electron.Angstrom): ",pol
       write(unit=*,fmt="(a,3f8.4,a)")  " => Cartesian Dipolar Moment vector(electron.Angstrom):(",r_pol," )"
       write(unit=*,fmt="(a,f8.4)")     " => Ionic Polarization/UnitVolume  (uCoulomb/cm^2)    : ",polc
       write(unit=*,fmt="(a,3f8.4,a)")  " => Cartesian Polarization vector  (uCoulomb/cm^2)    :(",polc*dir_pol," )"
       write(unit=i_out,fmt="(/,a,3f8.4,a)")" => Geometric centre of Q+ (fract. coordinates)       :(",r_plus," )"
       write(unit=i_out,fmt="(a,3f8.4,a)")  " => Geometric centre of Q- (fract. coordinates)       :(",r_minus," )"
       write(unit=i_out,fmt="(a,f8.4)")     " => Distance between  Q+  and  Q- (Angstrom)          : ",dist
       write(unit=i_out,fmt="(a,g11.4)")    " => Ionic Dipolar Moment/UnitCell  (Coulomb.Metre)    :  ",1.602176565e-29*pol
       write(unit=i_out,fmt="(a,f8.4)")     " => Ionic Dipolar Moment/UnitCell  (electron.Angstrom): ",pol
       write(unit=i_out,fmt="(a,3f8.4,a)")  " => Cartesian Dipolar Moment vector(electron.Angstrom):(",r_pol," )"
       write(unit=i_out,fmt="(a,f8.4)")     " => Ionic Polarization/UnitVolume  (uCoulomb/cm^2)    : ",polc
       write(unit=i_out,fmt="(a,3f8.4,a)")  " => Cartesian Polarization vector  (uCoulomb/cm^2)    :(",polc*dir_pol," )"
       if(pol > eps) then
         cpos=Cart_Vector("D",[1.0_cp,0.0_cp,0.0_cp],Cell)
         ang=Angle_Vect(cpos, r_pol)
         write(unit=*,fmt="(a,f7.1,a)")      " => Angle of Polarisation vector with a-axis          :",ang," degrees"
         write(unit=i_out,fmt="(a,f7.1,a)")  " => Angle of Polarisation vector with a-axis          :",ang," degrees"
         cpos=Cart_Vector("D",[0.0_cp,1.0_cp,0.0_cp],Cell)
         ang=Angle_Vect(cpos, r_pol)
         write(unit=*,fmt="(a,f7.1,a)")      " => Angle of Polarisation vector with b-axis          :",ang," degrees"
         write(unit=i_out,fmt="(a,f7.1,a)")  " => Angle of Polarisation vector with b-axis          :",ang," degrees"
         cpos=Cart_Vector("D",[0.0_cp,0.0_cp,1.0_cp],Cell)
         ang=Angle_Vect(cpos, r_pol)
         write(unit=*,fmt="(a,f7.1,a)")      " => Angle of Polarisation vector with c-axis          :",ang," degrees"
         write(unit=i_out,fmt="(a,f7.1,a)")  " => Angle of Polarisation vector with c-axis          :",ang," degrees"
       end if
       call Wait_Message(" => Press <enter> to continue ...")

    End Subroutine Menu_Atom_4

    Subroutine Menu_Atom_5()
       !---- Local Variables ----!
       integer                        :: i,j,Mult,m,lam,lbm,lcm,np,nm
       real(kind=cp)                  :: q, qp,pol,ang,ncell,dist,polc
       real(kind=cp), dimension(3)    :: pos,cpos,r_frac, r_pol,r_plus,r_minus,dir_pol
       real(kind=cp), dimension(3,384):: orb !increased to take into account the surface atoms
       real(kind=cp), dimension(  384):: qat !New charges
       logical                        :: calc_possible=.true.

       call execute_command_line(clear_string)
       write(unit=*,fmt="(a)") " "
       write(unit=*,fmt="(a)") "               GENERAL CRYSTALLOGRAPHY CALCULATOR "
       write(unit=*,fmt="(a)") " "
       write(unit=*,fmt="(a)") "                     Atomistic Calculations "
       write(unit=*,fmt="(a)") " ======================================================================="
       write(unit=*,fmt="(a)") " Calculating Ionic polarisation with respect to a non-polar structure..."
       write(unit=*,fmt="(a)") " The non-polar structure should be given in the same spacegroup as the"
       write(unit=*,fmt="(a)") " polar structure in order to calculate the displacement parameters in a"
       write(unit=*,fmt="(a)") " straightforward way by calculating difference of coordinates."
       write(unit=*,fmt="(a)") "  "
       write(unit=*,fmt="(a)",advance="no") " => Enter the name of the CFL file containing the non polar structure: "
       read(unit=*,fmt="(a)")


       do i=1,A%natoms
         if(abs(A%atom(i)%charge) <= 0.001)  then
           calc_possible=.false.
           exit
         end if
       end do
       if( .not. calc_possible) then
         write(unit=*,fmt="(a)") " => Calculation of P impossible. No charges have been provided! "
         call Wait_Message(" => Press <enter> to continue ...")
         return
       end if
       !if(SpG%Centred /= 1) then
       !  write(unit=*,fmt="(a)") " => Polarisation = 0.0 Centrosymmetric Crystal! "
       !  call Wait_Message(" => Press <enter> to continue ...")
       !  return
       !end if
      ! write(unit=*,fmt="(a)",advance="no")  " => Enter supercell extend (la,lb,lcm): "
      ! read(unit=*,fmt=*) lam,lbm,lcm
       lam=0; lbm=0; lcm=0
       ncell=(2*lam+1)*(2*lbm+1)*(2*lcm+1)
       !write(unit=*,fmt="(a,i10,a)")  " => Calculation for: ",nint(ncell)," unit cells"
       write(unit=i_out,fmt="(a)") " "
       write(unit=i_out,fmt="(a)") "                     Atomistic Calculations "
       write(unit=i_out,fmt="(a)") " ======================================================================="
       write(unit=i_out,fmt="(a)") " Calculating the ionic polarisation of a single symmetrized unit cell..."
       write(unit=i_out,fmt="(a)") " "
       r_plus=0.0; r_minus=0.0 ; np=0; nm=0
       do i=1,A%natoms
          pos=A%atom(i)%x
          q=A%atom(i)%charge
          call Get_Orbit(pos,SpG,Mult,orb) !here the orbit has always positive coordinates
          !Treatment of atoms at origin, edges and faces
          m=0; qat=0.0
          do j=1,Mult
            !Testing 0,0,0
            if(abs(orb(1,j)) < eps .and. abs(orb(2,j)) < eps  .and. abs(orb(3,j)) < eps  ) Then !atom at the origin
              qp=0.125*q
              qat(j)=qp
              m=m+1
              orb(:,Mult+m) = [1.0_cp,0.0_cp,0.0_cp] ; qat(Mult+m)=qp
              m=m+1
              orb(:,Mult+m) = [0.0_cp,1.0_cp,0.0_cp] ; qat(Mult+m)=qp
              m=m+1
              orb(:,Mult+m) = [0.0_cp,0.0_cp,1.0_cp] ; qat(Mult+m)=qp
              m=m+1
              orb(:,Mult+m) = [1.0_cp,1.0_cp,0.0_cp] ; qat(Mult+m)=qp
              m=m+1
              orb(:,Mult+m) = [1.0_cp,0.0_cp,1.0_cp] ; qat(Mult+m)=qp
              m=m+1
              orb(:,Mult+m) = [0.0_cp,1.0_cp,1.0_cp] ; qat(Mult+m)=qp
              m=m+1
              orb(:,Mult+m) = [1.0_cp,1.0_cp,1.0_cp] ; qat(Mult+m)=qp
            else if(abs(orb(1,j)) < eps .and. abs(orb(2,j)) < eps  ) Then !atom in z-edge
              qp=0.25*q
              qat(j)=qp
              m=m+1
              orb(:,Mult+m) = orb(:,j) + [1.0_cp,0.0_cp,0.0_cp]  ; qat(Mult+m)=qp
              m=m+1
              orb(:,Mult+m) = orb(:,j) + [0.0_cp,1.0_cp,0.0_cp]  ; qat(Mult+m)=qp
              m=m+1
              orb(:,Mult+m) = orb(:,j) + [1.0_cp,1.0_cp,0.0_cp]  ; qat(Mult+m)=qp
            else if(abs(orb(1,j)) < eps .and. abs(orb(3,j)) < eps  ) Then !atom in y-edge
              qp=0.25*q
              qat(j)=qp
              m=m+1
              orb(:,Mult+m) = orb(:,j) + [1.0_cp,0.0_cp,0.0_cp]  ; qat(Mult+m)=qp
              m=m+1
              orb(:,Mult+m) = orb(:,j) + [0.0_cp,0.0_cp,1.0_cp]  ; qat(Mult+m)=qp
              m=m+1
              orb(:,Mult+m) = orb(:,j) + [1.0_cp,0.0_cp,1.0_cp]  ; qat(Mult+m)=qp
            else if(abs(orb(2,j)) < eps .and. abs(orb(3,j)) < eps ) Then !atom in x-edge
              qp=0.25*q
              qat(j)=qp
              m=m+1
              orb(:,Mult+m) = orb(:,j) + [0.0_cp,1.0_cp,0.0_cp]  ; qat(Mult+m)=qp
              m=m+1
              orb(:,Mult+m) = orb(:,j) + [0.0_cp,0.0_cp,1.0_cp]  ; qat(Mult+m)=qp
              m=m+1
              orb(:,Mult+m) = orb(:,j) + [0.0_cp,1.0_cp,1.0_cp]  ; qat(Mult+m)=qp
            else if(abs(orb(1,j)) < eps ) Then !atom in yz-plane
              qp=0.5*q
              qat(j)=qp
              m=m+1
              orb(:,Mult+m) = orb(:,j) + [1.0_cp,0.0_cp,0.0_cp]  ; qat(Mult+m)=qp
            else if(abs(orb(2,j)) < eps ) Then !atom in xz-plane
              qp=0.5*q
              qat(j)=qp
              m=m+1
              orb(:,Mult+m) = orb(:,j) + [0.0_cp,1.0_cp,0.0_cp]  ; qat(Mult+m)=qp
            else if(abs(orb(3,j)) < eps ) Then !atom in xy-plane
              qp=0.5*q
              qat(j)=qp
              m=m+1
              orb(:,Mult+m) = orb(:,j) + [0.0_cp,0.0_cp,1.0_cp]  ; qat(Mult+m)=qp
            else
              qat(j)=q
            end if
          end do
          Mult=Mult+m
          write(unit=*,fmt="(a,a15,f8.2,a,i3)") " => Atom: ",trim(A%Atom(i)%Lab)//",  Charge: ", A%Atom(i)%Charge, &
                                                "     Modified Multiplicity:",mult
          write(unit=i_out,fmt="(a,a15,f8.2,a,i3)") " => Atom: ",trim(A%Atom(i)%Lab)//",  Charge: ", A%Atom(i)%Charge, &
                                                "     Modified Multiplicity:",mult
          cpos=0.0
          do j=1,Mult
            write(unit=*,fmt="(a,3f9.5,a,f9.5)") " => (x,y,z): ",orb(:,j),"  Qeff=",qat(j)
            write(unit=i_out,fmt="(a,3f9.5,a,f9.5)") " => (x,y,z): ",orb(:,j),"  Qeff=",qat(j)
            !do la=-lam,lam
            !  do lb=-lbm,lbm
            !    do lc=-lcm,lcm
                  pos=orb(:,j)-[0.5_cp,0.5_cp,0.5_cp]   !+[la,lb,lc]
                  r_frac=r_frac+ pos*qat(j)
                  cpos=cpos+pos
                  if(qat(j) > eps) then
                    np=np+1
                    r_plus=r_plus+orb(:,j)
                  else
                    nm=nm+1
                    r_minus=r_minus+orb(:,j)
                  end if
            !    end do
            !  end do
            !end do
          end do
         cpos=cpos/(ncell * real(mult,kind=cp))
         write(unit=*,fmt="(a,3f10.5,a)") "    Average of vector positions w.r.to (1/2,1/2,1/2): (",cpos," )"
         write(unit=i_out,fmt="(a,3f10.5,a)") "    Average of vector positions w.r.to (1/2,1/2,1/2): (",cpos," )"
       end do
       r_plus=r_plus/np
       r_minus=r_minus/nm
       !r_frac=r_frac/real(ncell)
       r_frac=(r_plus-r_minus) * qcellp
       r_pol=Cart_Vector("D",r_frac,Cell)
       where(abs(r_pol) <= 0.0001) r_pol=0.0
       pol=sqrt(dot_product(r_pol,r_pol))
       dir_pol=[0.0,0.0,1.0]
       if(pol > eps) dir_pol=r_pol/pol
       cpos=Cart_Vector("D",r_plus-r_minus,Cell)
       dist=sqrt(dot_product(cpos,cpos))
       polc=1602.176565*pol/Cell%CellVol
       write(unit=*,fmt="(/,a,3f8.4,a)")" => Geometric centre of Q+ (fract. coordinates)       :(",r_plus," )"
       write(unit=*,fmt="(a,3f8.4,a)")  " => Geometric centre of Q- (fract. coordinates)       :(",r_minus," )"
       write(unit=*,fmt="(a,f8.4)")     " => Distance between  Q+  and  Q- (Angstrom)          : ",dist
       write(unit=*,fmt="(a,g11.4)")    " => Ionic Dipolar Moment/UnitCell  (Coulomb.Metre)    :  ",1.602176565e-29*pol
       write(unit=*,fmt="(a,f8.4)")     " => Ionic Dipolar Moment/UnitCell  (electron.Angstrom): ",pol
       write(unit=*,fmt="(a,3f8.4,a)")  " => Cartesian Dipolar Moment vector(electron.Angstrom):(",r_pol," )"
       write(unit=*,fmt="(a,f8.4)")     " => Ionic Polarization/UnitVolume  (uCoulomb/cm^2)    : ",polc
       write(unit=*,fmt="(a,3f8.4,a)")  " => Cartesian Polarization vector  (uCoulomb/cm^2)    :(",polc*dir_pol," )"
       write(unit=i_out,fmt="(/,a,3f8.4,a)")" => Geometric centre of Q+ (fract. coordinates)       :(",r_plus," )"
       write(unit=i_out,fmt="(a,3f8.4,a)")  " => Geometric centre of Q- (fract. coordinates)       :(",r_minus," )"
       write(unit=i_out,fmt="(a,f8.4)")     " => Distance between  Q+  and  Q- (Angstrom)          : ",dist
       write(unit=i_out,fmt="(a,g11.4)")    " => Ionic Dipolar Moment/UnitCell  (Coulomb.Metre)    :  ",1.602176565e-29*pol
       write(unit=i_out,fmt="(a,f8.4)")     " => Ionic Dipolar Moment/UnitCell  (electron.Angstrom): ",pol
       write(unit=i_out,fmt="(a,3f8.4,a)")  " => Cartesian Dipolar Moment vector(electron.Angstrom):(",r_pol," )"
       write(unit=i_out,fmt="(a,f8.4)")     " => Ionic Polarization/UnitVolume  (uCoulomb/cm^2)    : ",polc
       write(unit=i_out,fmt="(a,3f8.4,a)")  " => Cartesian Polarization vector  (uCoulomb/cm^2)    :(",polc*dir_pol," )"
       if(pol > eps) then
         cpos=Cart_Vector("D",[1.0_cp,0.0_cp,0.0_cp],Cell)
         ang=Angle_Vect(cpos, r_pol)
         write(unit=*,fmt="(a,f7.1,a)")      " => Angle of Polarisation vector with a-axis          :",ang," degrees"
         write(unit=i_out,fmt="(a,f7.1,a)")  " => Angle of Polarisation vector with a-axis          :",ang," degrees"
         cpos=Cart_Vector("D",[0.0_cp,1.0_cp,0.0_cp],Cell)
         ang=Angle_Vect(cpos, r_pol)
         write(unit=*,fmt="(a,f7.1,a)")      " => Angle of Polarisation vector with b-axis          :",ang," degrees"
         write(unit=i_out,fmt="(a,f7.1,a)")  " => Angle of Polarisation vector with b-axis          :",ang," degrees"
         cpos=Cart_Vector("D",[0.0_cp,0.0_cp,1.0_cp],Cell)
         ang=Angle_Vect(cpos, r_pol)
         write(unit=*,fmt="(a,f7.1,a)")      " => Angle of Polarisation vector with c-axis          :",ang," degrees"
         write(unit=i_out,fmt="(a,f7.1,a)")  " => Angle of Polarisation vector with c-axis          :",ang," degrees"
       end if
       call Wait_Message(" => Press <enter> to continue ...")

    End Subroutine Menu_Atom_5

    Function Angle_vect(u,v) Result(angle)
      real(kind=cp), dimension(:), intent(in) :: u,v
      real(kind=cp) :: angle
      real(kind=cp) :: mu, mv
      mu=sqrt(dot_Product(u,u))
      mv=sqrt(dot_Product(v,v))
      if(mu < eps .or. mv < eps) then
          write(unit=*,fmt="(a)") " -> One of the directions is [0 0 0] ... retry!"
          call Wait_Message(" => Press <enter> to continue ...")
          return
      end if
      angle=dot_Product(u,v)/mu/mv
      if(angle > 1.0) angle=1.0
      if(angle < -1.0) angle=-1.0
      angle=acosd(angle)
      return
    End Function Angle_vect

end module Menu_3
