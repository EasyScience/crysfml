!!---------------------------------------------------------------------------------
!!--- Program SIMBO
!!--- Purpose: Generate files with neighbouring information around magnetic atoms
!!---          for simulation purposes.
!!---          Programs using output files from SIMBO: EnerMag, MC_MAG
!!---          Adapted from the old program SIMBO by P. Lacorre and J.Rodriguez-Carvajal
!!---          which was written from SIMilar and BOndstr (JRC) in January 1995.
!!--- (C) CopyLeft JRC-LLB, August 1999.
!!--- Author: Juan Rodriguez-Carvajal (Laboratoire Leon Brillouin, CEA-CNRS)
!!---         Transformed to F-language in July 2002.
!!---------------------------------------------------------------------------------
  Module Simbo_mod

    Use CFML_Math_general,              only: sort, negligible
    Use CFML_Crystallographic_Symmetry, only: Space_Group_Type, nlat
    Use CFML_String_Utilities,          only: Pack_String, Frac_Trans_1Dig
    use CFML_Geometry_Calc,             only: angle_uv, angle_dihedral,P1_dist, &
                                              Err_Geom,Err_Geom_Mess
    Use CFML_Atom_TypeDef,              only: Atoms_Cell_Type,  equiv_atm
    use CFML_crystal_metrics,           only: Crystal_Cell_Type
    Use Super_Exchange

    implicit none
    private

    public:: Exchange_Paths, construct_jxch
    real, parameter, private :: epsi=0.001

    !Transfomation matrices to put coordinates in a primitive cell
    !according to the equation Xp=M_Z Xz, where X are column matrices
    !of coordinate triplets.

    real, dimension(3,3), parameter,public :: M_A=reshape((/ 1.0, 0.0, 0.0, &
                                                             0.0, 1.0, 1.0, &
                                                             0.0,-1.0, 1.0 /),(/3,3/))
    real, dimension(3,3), parameter,public :: M_B=reshape((/ 1.0, 0.0,-1.0, &
                                                             0.0, 1.0, 0.0, &
                                                             1.0, 0.0, 1.0 /),(/3,3/))
    real, dimension(3,3), parameter,public :: M_C=reshape((/ 1.0,-1.0, 0.0, &
                                                             1.0, 1.0, 0.0, &
                                                             0.0, 0.0, 1.0 /),(/3,3/))
    real, dimension(3,3), parameter,public :: M_I=reshape((/ 1.0, 0.0, 0.0, &
                                                             0.0, 1.0, 0.0, &
                                                            -1.0,-1.0, 2.0 /),(/3,3/))
    real, dimension(3,3), parameter,public :: M_R=reshape((/ 1.0, 0.0,-1.0, &
                                                             1.0, 1.0, 1.0, &
                                                            -1.0,-1.0, 2.0 /),(/3,3/))
    real, dimension(3,3), parameter,public :: M_F=reshape((/ 1.0, 1.0,-1.0, &
                                                            -1.0, 1.0, 1.0, &
                                                             1.0,-1.0, 1.0 /),(/3,3/))
    real, dimension(3,3),public :: M_coor, M_basis

    contains

   Subroutine Exchange_Paths(lun,iprin,dmax,dbond,angm,angn,directex1,directex2,Cell,SpG,Ac,spaths)
     integer,                  intent(in)       :: lun
     logical,                  intent(in)       :: iprin
     real,                     intent(in)       :: dmax
     real,                     intent(in)       :: dbond,angm,angn,directex1,directex2
     type (Crystal_cell_type), intent(in)       :: Cell
     type (Space_Group_type),  intent(in)       :: SpG
     type (Atoms_Cell_type),   intent(in out)   :: Ac
     type (SE_Connection), dimension (:,:),intent(in out) :: spaths
     !-- Local Variables --!
     integer, dimension(Ac%nat) :: ind_mag
     integer                    :: i,j,k,ki,kk,ji,jk,nsij,nssij
     integer                    :: n_mag,im,km
     real, dimension(3)         :: vm,vmp,va,vap,vmmp,vma,vmpap,vaap
     real                       :: d2, ang, ang1,ang2, ang3, dis,dir1,dir2
     character (len=60)         :: tangl
     character (len=40)         :: translat

     d2= dbond*dbond   !Square of the distance cation-anion for s-ex and ss-ex paths
                       !Square of the distance  anion-anion for s-ex and ss-ex paths
     dir1=directex1*directex1   !Square of the distance for minimal direct exchange
     dir2=directex2*directex2   !Square of the distance for maximal direct exchange

     !Construct full connectivity of atoms in the cell without printing
     !The structure Ac(Atoms_Cell_type) contains all atoms + coordinations, etc...
     call P1_dist( dmax, Cell, SpG, Ac)
     !call P1_dist( dmax, Cell, SpG, Ac,lun)


     if(Err_Geom) then
       write(unit=*,fmt="(/,a)") " => WARNING"
       write(unit=*,fmt="(a)") trim(ERR_Geom_Mess)
       write(unit=*,fmt="(a,/)") " => The calculation of exchange paths may be wrong!"
     end if
     !Print the full structure Ac
     !call print_AC(Ac,directex1,directex2,lun)

     !Inititalize spaths
     spaths(:,:)%nd=0
     spaths(:,:)%ns=0
     spaths(:,:)%nss=0


     do i=1,num_de
       spaths(:,:)%DE(i)%nam(1:)=" "
       do j=1,3
        do k=1,2
         spaths(:,:)%DE(i)%coord(j,k)=0.0
        end do
       end do
     end do

     do i=1,num_se
      spaths(:,:)%SE(i)%nam(1:)=" "
      do j=1,4
       spaths(:,:)%SE(i)%geom(j)=0.0
      end do
       do j=1,3
        do k=1,3
         spaths(:,:)%SE(i)%coord(j,k)=0.0
        end do
       end do
     end do

     do i=1,num_sse
      spaths(:,:)%SSE(i)%nam(1:)=" "
      do j=1,7
       spaths(:,:)%SSE(i)%geom(j)=0.0
      end do
       do j=1,3
        do k=1,4
         spaths(:,:)%SSE(i)%coord(j,k)=0.0
        end do
       end do
     end do

     if(iprin) then
       write(unit=lun,fmt="(/,/,a)") "   ANALYSIS OF EXCHANGE PATHS"
       write(unit=lun,fmt="( a,/)") "   =========================="
     end if
     write(unit=lun,fmt="(a,f9.3,a)")  "   Maximum distance between magnetic atoms:",dmax,  " angstroms"
     write(unit=lun,fmt="(a,f9.3,a)")  "   Minimum distance for direct exchange   :",directex1," angstrom"
     write(unit=lun,fmt="(a,f9.3,a)")  "   Maximum distance for direct exchange   :",directex2," angstrom"
     write(unit=lun,fmt="(a,f9.3,a)")  "   Maximum distance between anions        :",dbond, &
                                       " angstroms (also for cation-anion bonds)"
     write(unit=lun,fmt="(a,f9.3,a )") "   Maximum angle between M-M' and M-A     :",angm , " degrees"
     write(unit=lun,fmt="(a,f9.3,a,/)")"   Maximum angle between M-M' and M-A     :",angn , " degrees"

     !     Determine the magnetic atoms in the list Ac and number them
     n_mag=0
     do i=1,Ac%nat
        if(Ac%moment(i) < 0.01) cycle
        n_mag=n_mag+1
        ind_mag(n_mag)=i
     end do
     !
     !Working now with the structure of Ac
     !to determine the exchange paths
     !
     nsij=0    !Number of super-exchange paths
     nssij=0   !Number of super-super-exchange paths
     do i=1,Ac%nat
       if(Ac%moment(i) < 0.01) cycle              !Select magnetic atom i (M)
       !Determine the number of the magnetic atom
       do im=1,n_mag
         if(ind_mag(im) == i) exit !im gets the correct value
       end do
       vm=MATMUL(Cell%Cr_Orth_cel,Ac%xyz(:,i))   !cartesian coordinates of atom i (M)

       do j=1,Ac%neighb(i)                       !loop over the neighbours of atom i
        k=Ac%neighb_atom(j,i)                    !k is the index in Ac of the j neighbour of i
        if(Ac%moment(k) < 0.01) cycle            !Select magnetic atom k (M') connected to i (M)
        !now k is a magnetic atom
        !Determine the number of the magnetic atom
        do km=1,n_mag
         if(ind_mag(km) == k) exit !km gets the correct value
        end do

        !We will construct the matrix element (im,km) with all superexchange paths
        vmp=MATMUL(Cell%Cr_Orth_cel,Ac%xyz(:,k)+Ac%trans(:,j,i))   !cartesian coordinates of atom j (M')
        vmmp=vmp-vm                    !interatomic vector between magnetic atoms in cartesian components
        call Frac_Trans_1Dig(Ac%trans(:,j,i),translat)
        !translation of the atom "k" w.r.t to that situated within the reference cell
        translat=Pack_String(translat)


        !look for non-magnetic atoms connected to i and k
         do ji=1,Ac%neighb(i)
           ki=Ac%neighb_atom(ji,i)    !non-magnetic atom (anion A) connected to i (M)
           if(Ac%moment(ki)> 0.01 .or. Ac%charge(ki) > 0.0 ) cycle    !discard cations
           va=MATMUL(Cell%Cr_Orth_cel,Ac%xyz(:,ki)+Ac%trans(:,ji,i))  !Cartesian coordinates of anion A
           vma=va-vm     !interatomic vector MA
                         !                              A--->A'
                         !discard paths with atoms: vma |    |
                         !and |vma|>dbond               M--->M'
                         !                               vmmp
           if(Angle_uv(vmmp,vma) > angm ) then
             !  write(unit=*,fmt=*) "  Angle A-M-M': ", Angle_uv(vmmp,vma), "  should be below:", angm
             cycle
           end if
           if(dot_PRODUCT(vma,vma) > d2) cycle   !d2 is the maximum allowed square of d(M-A)


          !In the analysis of the atoms connected to k one has to take into account that
          !the construction of neighbouring atoms has been made using the atoms inside
          !the unit cell, so a translation Ac%trans(:,j,i) (translation of connected k-atom to i
          !w.r.t. the original) must be added.
          !
          do jk=1,Ac%neighb(k)
            kk=Ac%neighb_atom(jk,k)   !non-magnetic atom (anion A') connected to k (M')
            if(Ac%moment(kk)> 0.01 .or. Ac%charge(kk) > 0.0 ) cycle   !discard cations
            vap=MATMUL(Cell%Cr_Orth_cel,Ac%xyz(:,kk)+Ac%trans(:,jk,k)+Ac%trans(:,j,i))
            vmpap=vap-vmp !interatomic vector M'A'
                          !                             A<---A'
                          !discard paths with atoms:    |    | vmpap
                          !and |vmpap|>dbond            M<---M'
                          !                              -vmmp
            if(Angle_uv(-vmmp,vmpap) > angm ) then
             !write(unit=*,fmt=*) " Angle A'-M'-M: ", Angle_uv(-vmmp,vmpap), " should be below:",angm
             cycle
            end if
            if(dot_PRODUCT(vmpap,vmpap) > d2) cycle

             vaap=vap-va

            if(negligible(SUM( abs(vaap(:)) )) ) then  !Eventual super-exchange path
                                        !                             A=A'
                                        !Super-Exchange path         /   \
                                        !                           M-----M'
              nsij=nsij+1               !                             vmmp
              tangl=" "
              write(unit=tangl,fmt="(5a)") trim(Ac%noms(i)),"-",trim(Ac%noms(ki)),"-",trim(Ac%noms(k))

              ang=Angle_uv(vma,vmpap)
              !construct spaths
              spaths(im,km)%ns=spaths(im,km)%ns+1
              if(spaths(im,km)%ns > num_se) then
                write(unit=*,fmt="(4a)") &
                " => WARNING!: too many Super-Exchange paths between atoms:", &
                Ac%noms(i)," and ",Ac%noms(k)
                cycle
              end if
              spaths(im,km)%SE(spaths(im,km)%ns)%nam1=trim(Ac%noms(i))
              spaths(im,km)%SE(spaths(im,km)%ns)%nam2=trim(Ac%noms(ki))
              spaths(im,km)%SE(spaths(im,km)%ns)%nam3=trim(Ac%noms(k))
              spaths(im,km)%SE(spaths(im,km)%ns)%nam=trim(tangl)//trim(translat)
              spaths(im,km)%SE(spaths(im,km)%ns)%geom(1)=Ac%distance(ji,i)
              spaths(im,km)%SE(spaths(im,km)%ns)%geom(2)=Ac%distance(jk,k)
              spaths(im,km)%SE(spaths(im,km)%ns)%geom(3)=ang
              spaths(im,km)%SE(spaths(im,km)%ns)%geom(4)=Ac%distance(j,i)
              spaths(im,km)%SE(spaths(im,km)%ns)%coord(:,1)=Ac%xyz(:,i)
              spaths(im,km)%SE(spaths(im,km)%ns)%coord(:,2)=Ac%xyz(:,ki)+Ac%trans(:,ji,i)
              spaths(im,km)%SE(spaths(im,km)%ns)%coord(:,3)=Ac%xyz(:,k) +Ac%trans(:,j,i)
              spaths(im,km)%SE(spaths(im,km)%ns)%carte(:,1)=vm(:)
              spaths(im,km)%SE(spaths(im,km)%ns)%carte(:,2)=va(:)
              spaths(im,km)%SE(spaths(im,km)%ns)%carte(:,3)=vmp(:)

            else

                   !Possible super-super-exchange path
                   !                          A'
                   !                        / |
                   !                  vaap /  |
                   !discard paths with:   /)a |
                   !  a > 89             A--  |
                   !and |vaap|>dbond     |    | vmpap
                   !                     M--->M'
                   !                      vmmp

             if(Angle_uv(vmmp,vaap) > angm) then
                !write(unit=*,fmt=*) " Angle AA'^MM': ", Angle_uv(vmmp,vaap), "  should be below:", angm
                cycle
             end if
             dis=dot_PRODUCT(vaap,vaap)
             if(dis > d2) cycle

             ang1=Angle_uv(-vma,vaap)
             ang2=Angle_uv(vaap,vmpap)             !

             if(ang1 < angn .or. ang2 < angn) then
                !write(unit=*,fmt=*) " Angles MAA' and M'A'A: ", ang1,ang2, "  must be >", angn
                cycle
             end if
                    !a super-exchange path exist
                    !
                    !
                    !discard paths with:
                    ! a1 <90 or a2 <90      A-------A'
                    !and |vaap|>dbond      / a1   a2 \ vmpap
                    !                     M----------->M'
                    !                       vmmp

             nssij=nssij+1
             dis=sqrt(dis)
             ang3=angle_dihedral(vma,vaap,-vmpap)

             tangl=" "
             write(unit=tangl,fmt="(7a)") trim(Ac%noms(i)),"-",trim(Ac%noms(ki)),"-",&
                                 trim(Ac%noms(kk)),"-",trim(Ac%noms(k))

              !construct spaths

              spaths(im,km)%nss=spaths(im,km)%nss+1
              if(spaths(im,km)%nss > num_sse) then
                write(unit=*,fmt="(4a)") &
                " => WARNING!: too many Super-Super-Exchange paths between atoms:", &
                Ac%noms(i)," and ",Ac%noms(k)
                cycle
              end if
              spaths(im,km)%SSE(spaths(im,km)%nss)%nam1=trim(Ac%noms(i))
              spaths(im,km)%SSE(spaths(im,km)%nss)%nam2=trim(Ac%noms(ki))
              spaths(im,km)%SSE(spaths(im,km)%nss)%nam3=trim(Ac%noms(kk))
              spaths(im,km)%SSE(spaths(im,km)%nss)%nam4=trim(Ac%noms(k))
              spaths(im,km)%SSE(spaths(im,km)%nss)%nam=trim(tangl)//trim(translat)
              spaths(im,km)%SSE(spaths(im,km)%nss)%geom(1)=Ac%distance(ji,i)
              spaths(im,km)%SSE(spaths(im,km)%nss)%geom(2)=dis
              spaths(im,km)%SSE(spaths(im,km)%nss)%geom(3)=Ac%distance(jk,k)
              spaths(im,km)%SSE(spaths(im,km)%nss)%geom(4)=ang1
              spaths(im,km)%SSE(spaths(im,km)%nss)%geom(5)=ang2
              spaths(im,km)%SSE(spaths(im,km)%nss)%geom(6)=ang3
              spaths(im,km)%SSE(spaths(im,km)%nss)%geom(7)=Ac%distance(j,i)
              spaths(im,km)%SSE(spaths(im,km)%nss)%coord(:,1)=Ac%xyz(:,i)
              spaths(im,km)%SSE(spaths(im,km)%nss)%coord(:,2)=Ac%xyz(:,ki)+Ac%trans(:,ji,i)
              spaths(im,km)%SSE(spaths(im,km)%nss)%coord(:,3)=Ac%xyz(:,kk)+Ac%trans(:,jk,k)+Ac%trans(:,j,i)
              spaths(im,km)%SSE(spaths(im,km)%nss)%coord(:,4)=Ac%xyz(:,k)+Ac%trans(:,j,i)
              spaths(im,km)%SSE(spaths(im,km)%nss)%carte(:,1)=vm(:)
              spaths(im,km)%SSE(spaths(im,km)%nss)%carte(:,2)=va(:)
              spaths(im,km)%SSE(spaths(im,km)%nss)%carte(:,3)=vap(:)
              spaths(im,km)%SSE(spaths(im,km)%nss)%carte(:,4)=vmp(:)

            end if

          end do  !jk
         end do  !ji

         ! Test for direct exchange
         dis=dot_PRODUCT(vmmp,vmmp)
         if(dis <= dir2 .and. dis > dir1) then
            spaths(im,km)%nd=spaths(im,km)%nd+1
            tangl=" "
            write(unit=tangl,fmt="(3a)") trim(Ac%noms(i)),"-",trim(Ac%noms(k))
            spaths(im,km)%DE(spaths(im,km)%nd)%nam1=trim(Ac%noms(i))
            spaths(im,km)%DE(spaths(im,km)%nd)%nam2=trim(Ac%noms(k))
            spaths(im,km)%DE(spaths(im,km)%nd)%nam=trim(tangl)//trim(translat)
            spaths(im,km)%DE(spaths(im,km)%nd)%coord(:,1)=Ac%xyz(:,i)
            spaths(im,km)%DE(spaths(im,km)%nd)%coord(:,2)=Ac%xyz(:,k)+Ac%trans(:,j,i)
            spaths(im,km)%DE(spaths(im,km)%nd)%carte(:,1)=vm(:)
            spaths(im,km)%DE(spaths(im,km)%nd)%carte(:,2)=vmp(:)
            spaths(im,km)%DE(spaths(im,km)%nd)%dist=sqrt(dis)
         end if

       end do   !j magnetic neighbours of i
     end do !i

     if(iprin) then
      do im=1,n_mag
       do km=1,n_mag
        write(unit=lun,fmt="(/,/,a)")       "   ------------------------------------------------------"
        write(unit=lun,fmt="(a,i2,a,i2,a)") "   Exchange paths Contributing to matrix element: (",im,",",km,")"
        write(unit=lun,fmt="(a,/)")         "   ------------------------------------------------------"
        do nsij=1,spaths(im,km)%ns
          write(unit=lun,fmt="(/,a)")   "   --------------------------------------------------"
          write(unit=lun,fmt="(4a)")    "   Super-Exchange paths between atoms:",trim(spaths(im,km)%SE(nsij)%nam1),&
                                        " and ",trim(spaths(im,km)%SE(nsij)%nam3)
          write(unit=lun,fmt="(a,/)")   "   --------------------------------------------------"

          write(unit=lun,fmt="(/,/,a,i5,a,a,a,f8.4/)") "   Super-Exchange Path (M-A-M'):",nsij, &
                                "  ",trim(spaths(im,km)%SE(nsij)%nam), "  Distance = ", spaths(im,km)%SE(nsij)%geom(4)
          write(unit=lun,fmt="(a,f8.4,/,a,f8.4,/,a,f8.2)")  &
          "    Distance    d1(M-A ):",spaths(im,km)%SE(nsij)%geom(1), &
          "    Distance    d2(A-M'):",spaths(im,km)%SE(nsij)%geom(2), &
          "    Super-Exchange angle:",spaths(im,km)%SE(nsij)%geom(3)
          write(unit=lun,fmt="(a)") &
          "                         x      y     z          xc      yc      zc"
          write(unit=lun,fmt="(3(a,a,a,3f7.4,a,3f8.4,/))")  &
          "          ",spaths(im,km)%SE(nsij)%nam1, " ",spaths(im,km)%SE(nsij)%coord(:,1),"   ",spaths(im,km)%SE(nsij)%carte(:,1),&
          "          ",spaths(im,km)%SE(nsij)%nam2, " ",spaths(im,km)%SE(nsij)%coord(:,2),"   ",spaths(im,km)%SE(nsij)%carte(:,2),&
          "          ",spaths(im,km)%SE(nsij)%nam3, " ",spaths(im,km)%SE(nsij)%coord(:,3),"   ",spaths(im,km)%SE(nsij)%carte(:,3)

        end do !nsij

        do nsij=1,spaths(im,km)%nss
          write(unit=lun,fmt="(/,a)")   "   ---------------------------------------------------------"
          write(unit=lun,fmt="(4a)")    "   Super-Super-Exchange paths between atoms:",trim(spaths(im,km)%SSE(nsij)%nam1),&
                                        " and ",trim(spaths(im,km)%SSE(nsij)%nam4)
          write(unit=lun,fmt="(a,/)")   "   ---------------------------------------------------------"

          write(unit=lun,fmt="(/,/,a,i5,a,a,a,f8.4/)") "   Super-Super-Exchange Path (M-A-A'-M'):",nsij, &
                                "  ",trim(spaths(im,km)%SSE(nsij)%nam), "  Distance = ", spaths(im,km)%SSE(nsij)%geom(7)
          write(unit=lun,fmt="(3(a,f8.4,/),3(a,f8.2,/))")  &
          "    Distance         d1(M -A ):",spaths(im,km)%SSE(nsij)%geom(1), &
          "    Distance         d2(A -A'):",spaths(im,km)%SSE(nsij)%geom(2), &
          "    Distance         d2(A'-M'):",spaths(im,km)%SSE(nsij)%geom(3), &
          "    Bond-Angle       (M-A -A'):",spaths(im,km)%SSE(nsij)%geom(4), &
          "    Bond-Angle       (A-A'-M'):",spaths(im,km)%SSE(nsij)%geom(5), &
          "    Dihedral-angle (M-A-A'-M'):",spaths(im,km)%SSE(nsij)%geom(6)
          write(unit=lun,fmt="(a)") &
          "                         x      y     z          xc      yc      zc"
          write(unit=lun,fmt="(4(tr10,a,a,3f7.4,a,3f8.4,/))")  &
          spaths(im,km)%SSE(nsij)%nam1," ",spaths(im,km)%SSE(nsij)%coord(:,1),"   ",spaths(im,km)%SSE(nsij)%carte(:,1),&
          spaths(im,km)%SSE(nsij)%nam2," ",spaths(im,km)%SSE(nsij)%coord(:,2),"   ",spaths(im,km)%SSE(nsij)%carte(:,2),&
          spaths(im,km)%SSE(nsij)%nam3," ",spaths(im,km)%SSE(nsij)%coord(:,3),"   ",spaths(im,km)%SSE(nsij)%carte(:,3),&
          spaths(im,km)%SSE(nsij)%nam4," ",spaths(im,km)%SSE(nsij)%coord(:,4),"   ",spaths(im,km)%SSE(nsij)%carte(:,4)

        end do !nsij

        do nsij=1,spaths(im,km)%nd
        !if(spaths(im,km)%nd ==1) then
          write(unit=lun,fmt="(/,a)")   "   ---------------------------------------------------------"
          write(unit=lun,fmt="(4a)")    "   Direct-Exchange between atoms:",trim(spaths(im,km)%DE(nsij)%nam1), &
                                        " and ",trim(spaths(im,km)%DE(nsij)%nam2)
          write(unit=lun,fmt="(a,/)")   "   ---------------------------------------------------------"

          write(unit=lun,fmt="(/,/,a,a,a,a,f8.4/)") "   Direct-Exchange :", &
                                "  ",trim(spaths(im,km)%DE(nsij)%nam), "  Distance = ", spaths(im,km)%DE(nsij)%dist
          write(unit=lun,fmt="(a)") &
          "                         x      y     z          xc      yc      zc"
          write(unit=lun,fmt="(2(tr10,a,a,3f7.4,a,3f8.4,/))")  &
          spaths(im,km)%DE(nsij)%nam1," ",spaths(im,km)%DE(nsij)%coord(:,1),"   ",spaths(im,km)%DE(nsij)%carte(:,1),&
          spaths(im,km)%DE(nsij)%nam2," ",spaths(im,km)%DE(nsij)%coord(:,2),"   ",spaths(im,km)%DE(nsij)%carte(:,2)

        !end if !nd
        end do !nsij

       end do  !km
      end do   !im

     end if  !iprin

     return
   End Subroutine Exchange_Paths

   Subroutine Print_AC(Ac,d1,d2,lun)
     type (Atoms_Cell_type),   intent(in)   :: Ac
     real,                     intent(in)   :: d1,d2
     integer,                  intent(in)   :: lun
     integer :: i,j,k
     real    :: dd
     real, dimension(3) :: tn
     write(unit=lun,fmt="(a)") "  ------------------------------------------------------- "
     write(unit=lun,fmt="(a)") "  Connectivity of all magnetic atoms with direct Exchange "
     write(unit=lun,fmt="(a)") "  ------------------------------------------------------- "
     write(unit=lun,fmt="(a)") "    "
     do i=1,Ac%nat
        if(Ac%moment(i) < 0.01) cycle
        write(unit=lun,fmt="(/,a,i6)") " => Neighbours of atom "//trim(Ac%noms(i))//"      Total number: ",Ac%neighb(i)
        do j=1,Ac%neighb(i)
            k=Ac%neighb_atom(j,i)
            if(Ac%moment(k) < 0.01) cycle
            dd=Ac%distance(j,i)
            if(dd < d1 .or. dd > d2) cycle
            tn=Ac%trans(:,j,i)
            write(unit=lun,fmt="(a,f10.5,a,3f8.3,a)") " => Atom "//trim(Ac%noms(k))//" at distance: ",dd, "  Tn = (",tn,")"
        end do
     end do
     write(unit=lun,fmt="(/,a)") " "
     return
   End Subroutine Print_AC

   Subroutine construct_jxch(lun,n_mag,spaths,Acm,kf)
    integer,                            intent(in) :: lun,n_mag
    type (SE_Connection),dimension(:,:),intent(in) :: spaths
    type(Atoms_Cell_Type),              intent(in) :: Acm
    real, optional,                     intent(in) :: kf
    !
    !  Subroutine to determine and write the different exchange interactions
    !  given the list of exchange paths
    !
    integer, parameter                  :: max_jx=96
    integer                             :: i,j,k,im,km,j1,j2,n,L,nj,nt,tot_neigh
    integer, dimension(n_mag,n_mag)     :: nterms
    integer,          dimension(max_jx) :: p
    character(len=25),dimension(max_jx) :: trans
    !character(len=25),dimension(max_jx) :: Nam_12
    character(len=180)                  :: text
    character(len=60)                   :: expo
    character(len=25)                   :: transla
    real, dimension(3)                  :: vect,pos
    real :: dmin
    type (Exchange_interaction), dimension(n_mag,n_mag,max_jx) :: jota
    type (Exchange_interaction), dimension(max_jx)             :: jxch
    logical :: newt

    nterms=0
    nj=0
    do k=1,max_jx
      call init_exchange_interaction(jxch(k))
      do j=1,n_mag
         do i=1,n_mag
          call init_exchange_interaction(jota(i,j,k))
         end do
      end do
    end do

    do im=1,n_mag
      do km=1,n_mag
        n=0
        do i=1,spaths(im,km)%ns
           j1=index(spaths(im,km)%SE(i)%nam,"(")
           j2=index(spaths(im,km)%SE(i)%nam,")")
           n=n+1
           trans(n)=spaths(im,km)%SE(i)%nam(j1:j2)
        end do

        do i=1,spaths(im,km)%nss
           j1=index(spaths(im,km)%SSE(i)%nam,"(")
           j2=index(spaths(im,km)%SSE(i)%nam,")")
           n=n+1
           trans(n)=spaths(im,km)%SSE(i)%nam(j1:j2)
        end do

        do i=1,spaths(im,km)%nd
        !if(spaths(im,km)%nd ==1) then
           j1=index(spaths(im,km)%DE(i)%nam,"(")
           j2=index(spaths(im,km)%DE(i)%nam,")")
           n=n+1
           trans(n)=spaths(im,km)%DE(i)%nam(j1:j2)
        !end if
        end do

        L=0
        !Loop for making a pointer to paths for different terms in the element im,km
        if(n /= 0) then
          p(:)= 0
          p(1)=1
          L=1
          do i=2,n
            newt=.true.
            do j=i-1,1,-1
              if(trans(j) == trans(i)) then
                p(i)=p(j)
                newt=.false.
                exit
              end if
            end do
            if(newt) then
              L=L+1
              p(i)=L
            end if
          end do
        end if
        nterms(im,km)=L

        do nt=1,nterms(im,km)

          if(spaths(im,km)%ns /= 0) then
           j1=index(spaths(im,km)%SE(1)%nam1,"_")
           jota(im,km,nt)%nam1=spaths(im,km)%SE(1)%nam1(1:j1-1)
           j1=index(spaths(im,km)%SE(1)%nam3,"_")
           jota(im,km,nt)%nam2=spaths(im,km)%SE(1)%nam3(1:j1-1)

          else if(spaths(im,km)%nss /= 0) then
           j1=index(spaths(im,km)%SSE(1)%nam1,"_")
           jota(im,km,nt)%nam1=spaths(im,km)%SSE(1)%nam1(1:j1-1)
           j1=index(spaths(im,km)%SSE(1)%nam4,"_")
           jota(im,km,nt)%nam2=spaths(im,km)%SSE(1)%nam4(1:j1-1)

          else if(spaths(im,km)%nd /= 0) then
           j1=index(spaths(im,km)%DE(1)%nam1,"_")
           jota(im,km,nt)%nam1=spaths(im,km)%DE(1)%nam1(1:j1-1)
           j1=index(spaths(im,km)%DE(1)%nam2,"_")
           jota(im,km,nt)%nam2=spaths(im,km)%DE(1)%nam2(1:j1-1)
           jota(im,km,nt)%nde=1

          end if

           n=0
           do i=1,spaths(im,km)%ns
              n=n+1
              if(nt == p(n)) then
                jota(im,km,nt)%ns=jota(im,km,nt)%ns+1
                L=jota(im,km,nt)%ns
                jota(im,km,nt)%se_geo(1:3,L) = spaths(im,km)%SE(i)%geom(1:3)
                jota(im,km,nt)%dist          = spaths(im,km)%SE(i)%geom(4)
                jota(im,km,nt)%s_nam(L)      = spaths(im,km)%SE(i)%nam
              end if
           end do

           do i=1,spaths(im,km)%nss
              n=n+1
              if(nt == p(n)) then
                jota(im,km,nt)%nss=jota(im,km,nt)%nss+1
                L=jota(im,km,nt)%nss
                jota(im,km,nt)%sse_geo(1:6,L) = spaths(im,km)%SSE(i)%geom(1:6)
                jota(im,km,nt)%dist           = spaths(im,km)%SSE(i)%geom(7)
                jota(im,km,nt)%ss_nam(L)      = spaths(im,km)%SSE(i)%nam
              end if
           end do

           do i=1,spaths(im,km)%nd
           !if(spaths(im,km)%nd ==1) then
              n=n+1
              if(nt == p(n)) then
                jota(im,km,nt)%dist = spaths(im,km)%DE(i)%dist
                jota(im,km,nt)%de_nam = spaths(im,km)%DE(i)%nam
                jota(im,km,nt)%nde=1
              end if
           !end if
           end do
        end do

      end do !km
    end do !im

    dmin=1000.0
    do im=1,n_mag
      do km=1,n_mag
        do nt=1,nterms(im,km)
          if(jota(im,km,nt)%dist < dmin) then
            dmin= jota(im,km,nt)%dist
            i=im
            k=km
            n=nt
          end if
        end do
      end do
    end do

    jxch(1)=jota(i,k,n)  !i,j,k is selected from the minimum distance in the previous loop
    jxch(1)%J= "J1"
    jxch(1)%valj=-10.0
    jota(i,k,n)%valj=-10.0
    nj=1
    do im=1,n_mag
      do km=1,n_mag
        do nt=1,nterms(im,km)
          newt=.true.
          do i=1,nj
            if(Equiv_jotas(jota(im,km,nt),jxch(i)) ) then
              jota(im,km,nt)%J=jxch(i)%J
              jota(im,km,nt)%valj=jxch(i)%valj
              !Equalisation of the two types THIS IS WRONG!!!!!
              !the components of "jxch" are not good to be trasferred
              !jota(im,km,nt) = jxch(i)
              newt=.false.
              exit
            end if
          end do
          if(newt) then
            nj=nj+1
            jxch(nj)=jota(im,km,nt)
            if(nj <10) then
              write(unit=jxch(nj)%J,fmt="(a,i1)") "J",nj
            else
              write(unit=jxch(nj)%J,fmt="(a,i2)") "J",nj
            end if
            if(present(kf) .and. (jxch(nj)%ns + jxch(nj)%nss == 0) ) then
               jxch(nj)%valj=rkky(jxch(nj)%dist,kf)
            else
               jxch(nj)%valj=-10.0*(dmin/jxch(nj)%dist)**12
            end if
            jota(im,km,nt)%valj=jxch(nj)%valj
            jota(im,km,nt)%J=jxch(nj)%J
          end if
        end do
      end do
    end do


    write(unit=lun,fmt="(/,/,a)") "  ============================================================"
    write(unit=lun,fmt="(a)")     "  ======>  LIST OF INDEPENDENT EXCHANGE INTERACTIONS  <======="
    write(unit=lun,fmt="(a,/,/)") "  ============================================================"
    write(unit=lun,fmt="(a)")     " (Warning:  direct exchange interactions may not be independent)"

    do i=1,nj
        call write_exchange_interaction(lun,jxch(i))
    end do

    do im=1,n_mag
      do km=1,n_mag
        write(unit=lun,fmt="(/,/,a)")       "  ==========================================================="
        write(unit=lun,fmt="(a,i2,a,i2,a)") "  ====>  TERMS OF THE ELEMENT [",im,",",km,"] of the J(k)-MATRIX"
        write(unit=lun,fmt="(a,/,/)")       "  ==========================================================="
        do nt=1,nterms(im,km)
           write(unit=lun,fmt="(/,a,i3)")   " => Term number: ",nt
           call write_exchange_interaction(lun,jota(im,km,nt))
        end do
      end do
    end do


    write(unit=lun,fmt="(/,/,a)") " -----------------------------"
    write(unit=lun,fmt="(a)"  )   " Effective Neighbouring matrix"
    write(unit=lun,fmt="(a)")     " -----------------------------"
    write(unit=lun,fmt="(a)")
    write(unit=lun,fmt="(a,24i4)")"     ",(i,i=1,n_mag)
    write(unit=lun,fmt="(a)")
    DO i=1,n_mag
     write(unit=lun,fmt="(i3,a,24i4)")i,"  ",(nterms(i,k),k=1,n_mag)
    END DO
    write(unit=lun,fmt="(a,/,/)")
    write(unit=lun,fmt="(/,a,44a1,/,a,/,a,44a1,/)") " ", ("-",j=1,44),  &
       " => Terms of the exchange interaction matrix:"," ",("-",j=1,44)

    DO im=1,n_mag

      tot_neigh=sum(nterms(im,:))
      write(unit=3,fmt="(a)")"!Site Neighb   Dsing_Anis      Dir                Name        x         y         z"
      pos=Matmul(M_coor,Acm%xyz(:,im))
      write(unit=3,fmt="(i4,i5,i10,Tr6,3i4,a,a6,3f10.5)")im,tot_neigh,0,0,0,0,"          :: ",Acm%noms(im), pos
      write(unit=3,fmt="(a)")"!    Nav   Av  Bv  Cv        J"

      DO km=1,n_mag

        write(unit=4,fmt="(2i4,i5)") im,km,nterms(im,km)
        write(unit=lun,fmt="(/,a,3(i3,a),/)") " => J(",im,",",km,")[K]   (",nterms(im,km), " terms)"

        Do nt=1,nterms(im,km)
            transla=" "
          if(jota(im,km,nt)%ns /= 0) then
            i=index(jota(im,km,nt)%s_nam(1),"(")
            j=index(jota(im,km,nt)%s_nam(1),")")
            transla=jota(im,km,nt)%s_nam(1)(i:j)
          else if(jota(im,km,nt)%nss /= 0) then
            i=index(jota(im,km,nt)%ss_nam(1),"(")
            j=index(jota(im,km,nt)%ss_nam(1),")")
            transla=jota(im,km,nt)%ss_nam(1)(i:j)
          else
            i=index(jota(im,km,nt)%de_nam,"(")
            j=index(jota(im,km,nt)%de_nam,")")
            transla=jota(im,km,nt)%de_nam(i:j)
          end if

          text=jota(im,km,nt)%J//" : "//trim(jota(im,km,nt)%s_nam(1))//trim(jota(im,km,nt)%ss_nam(1))
          if(len_trim(text) == 0) text=trim(jota(im,km,nt)%nam1)//trim(jota(im,km,nt)%nam2)//" <-Direct-Exchange"
          !write(*,"(a,2i3,a,i3,a)") " => Element: (",im,km,") -> Term: ",nt,"  Trans:"//transla
          call Get_Expo(transla,expo)
          call Get_vect(transla,vect)
         !  read(1,"(5x,3f9.5,f10.3,3x,2a)")  &  !Enermag
         !      (trans(m,i,j,nt),m=1,3),exch,jota,name_jota
          write(unit=4,fmt="(a,3f9.5,f10.3,a,a,f8.4)") "     ",vect,&
                        jota(im,km,nt)%valj,"   -> ",trim(text)//" --> dist=", jota(im,km,nt)%dist
          write(unit=lun,fmt="(a,a,a,f8.4,a,a,a)")  &
                  "    Rn=", transla," dist=",jota(im,km,nt)%dist," --> ",jota(im,km,nt)%J,expo
          vect=Matmul(M_coor,vect)
          write(unit=3,fmt="(a,i3,a,3i4,a,f10.4,a)") "     ", &
                      km," ",nint(vect),"  ",jota(im,km,nt)%valj,"    "//jota(im,km,nt)%J
        end do

      END DO
    END DO
   ! Write file *.mcm for MCMAG

      !IF(nlat == 1) THEN
      !  open(unit=3,file=outfil(1:ln)//".mcm",status="replace",action="write",position="rewind")
      !  write(unit=3,fmt="(a)")title
      !  write(unit=3,fmt="(a)")" File created by program SIMBO"
      !  write(unit=3,fmt="(a,i4,a,i3)")" ",-Acm%nat,"   0 ",Acm%nat
      !  DO i=1,Acm%nat
      !    write(unit=3,fmt="(a,i3,a,i3,a,a6)")" ",i," ",Acm%neighb(i),"   ",Acm%noms(i)
      !    DO k=1,Acm%neighb(i)
      !       j=indg(i,k)
      !       n=jin(i,j)
      !       if(n == 0) cycle
      !      ! n=jin(i,k)
      !      write(unit=3,fmt="(a,i3,a,3i4,a,f9.4)") "     ", &
      !             Acm%neighb_atom(i,j)," ",INT(Acm%trans(:,i,j)),"  ",valj(n)
      !    END DO
      !  END DO
      !  write(unit=3,fmt="(a,i3,a)")"  1 ",Acm%nat,"  1.00"
      !  write(unit=3,fmt="(a,6f11.5)")" ", cell%cell,cell%ang
      !END IF

! END Write file *.mcm for MCMAG

    return
   End Subroutine construct_jxch

  End Module Simbo_mod

  Program Simbo
   Use CFML_Crystallographic_Symmetry, only: Set_SpaceGroup, Write_SpaceGroup, Space_Group_Type
   Use CFML_String_Utilities,          only: l_case, number_lines, reading_lines,u_case
   use CFML_Math_general,              only: negligible
   use CFML_Math_3D,                   only: set_eps,invert_a
   use CFML_Geometry_Calc,             only: Allocate_Coordination_Type, calc_dist_angle, &
                                             Err_Geom,ERR_Geom_Mess
   Use CFML_Crystal_Metrics
   Use CFML_IO_Formats
   Use CFML_Atom_TypeDef
   Use Simbo_mod
   Use Super_Exchange

   Implicit None

   character(len=120), allocatable, dimension(:) :: file_dat
   integer, parameter :: max_magt=96
   type (Crystal_Cell_Type) :: Cell, Celln
   type (Space_Group_Type)  :: Spg, gP1
   type (Atom_list_Type)    :: A       !Original list of atoms in the asymmetric unit
   type (Atom_list_Type)    :: Ap      !List of atoms inside a primitive cell
   type (Atoms_Cell_Type)   :: Acm, Ac !Magnetic atoms and all atoms inside a primitive cell
   type (SE_Connection), dimension (max_magt,max_magt) :: spaths     !a maximun of max_magt magnetic atoms in the cell
   type (Job_Info_type)     :: Job_Info
   character(len=1)    :: ans
   character(len=20)   :: sp1
   character(len=80)   :: title
   character(len=256)  :: infil,outfil,texto
   integer, parameter  :: lun1=1,lun2=6,lun=2
   integer :: i, j, numops, ln, nauas, nmag, lr, max_coord, L
   integer, dimension(:), allocatable :: ptr
   integer, dimension(10) :: nif
   real, dimension(10)  :: mom
   character(len=4), dimension(10) :: scf=" "
   real  :: dmax=6.0, & !Maximum distance for distance calculations
            dangl=0.0,& !Maximun distance for angle calculations
            dbond=3.0,& !Maximun distance between anions involved in a SS-exchange path
            angm=89.5,& !Maximum angle between M-M' and M(M')-A(A') or M-M' and AA'
            angn=89.5,& !Minimun angle between M-A  and AA' to consider a possible SSE path
            directex1=0.0, directex2=0.0,&  !Minimum and Maximum distances to consider a possible direct exchange path
            kF=0.0     ! Fermi wave-vector
   real  :: seconds, End_time, start_time, rminutes, hours
   logical :: iprin, mag=.false.
   integer :: narg

   !---- Arguments on the command line ----!
   lr=0
   ln=0
   narg=command_argument_count()
   if(narg /= 0) then
     call get_command_argument(1,infil)
     i=index(infil,".cfl")
     if(i /= 0) infil=infil(1:i-1)
     ln=len_trim(infil)
     outfil=infil
     lr=ln
   end if
   if(narg > 1) then
     call get_command_argument(2,outfil)
     lr=len_trim(outfil)
   end if

   write(unit=*,fmt="(/,/,6(a,/))")                                       &
        "                  ------ PROGRAM SIMBO ------"                 , &
        "                 ---- Version 2.0 Oct-2003----"                , &
        "    ********************************************************"  , &
        "    * Generates neighboring files for magnetic simulations *"  , &
        "    ********************************************************"  , &
        "                      (JRC- October 2003 )"
   write(unit=*,fmt=*) " "

   do
      if(lr == 0) then
            write(unit=*,fmt="(a)",advance="no") " => Code of the file xx.cfl (give xx): "
            read(unit=*,fmt="(a)") infil
            write(unit=*,fmt="(a,a,a)",advance="no") &
                    " => Code of the output files (.nei,.exc) ( <cr>= ",trim(infil),") :"
            read(unit=*,fmt="(a)") outfil
            lr=len_trim(outfil)
            if(lr == 0) outfil=infil
      end if

      CALL CPU_TIME(start_time)
      call Readn_set_Xtal_Structure(trim(infil)//".cfl",Cell,SpG,A,Iphase=1,Job_Info=Job_Info)
      If(err_form) then
         write(unit=*,fmt="(a)") trim(ERR_Form_Mess)
         stop
      else
         exit
      end if
   end do

   call set_eps(0.001)  !In mathematical comparisons -> math_mod

   open(unit=lun,file=trim(outfil)//".nei",status="replace",action="write",position="rewind")
   !Writing titles and content of the imput file
     write(unit=lun,fmt="(/,/,6(a,/))")                                     &
          "                  ------ PROGRAM SIMBO ------"                 , &
          "                 ---- Version 2.0 Oct-2003----"                , &
          "    ********************************************************"  , &
          "    * Generates neighboring files for magnetic simulations *"  , &
          "    ********************************************************"  , &
          "                      (JRC- October 2003 )"

   write(unit=lun,fmt="(/,a,a,/)")" => Information deduced from the content of the input file: ",infil(1:ln)//".cfl"
   title=Job_Info%title
   write(unit=lun,fmt="(a,a,/)") " => TITLE: ",trim(title)
   call Write_Crystal_Cell(Cell,lun)
   call Write_SpaceGroup(SpG,lun,.true.)
   call Write_Atom_List(A,Level=1,Lun=lun,Cell=Cell)


   nauas=A%natoms
   mag=.false.
   write(unit=lun,fmt="(/,/,a,/)") " => List of detected magnetic atoms:"
   write(unit=lun,fmt="(a)") " Label  Type     x       y       z      occ     Biso  moment  Charge"
   do i=1,A%natoms
      write(unit=texto,fmt="(a,a4,a,a5,5f8.4,2f8.3)")" ", A%atom(i)%lab,"  ",&
      A%atom(i)%ChemSymb,A%atom(i)%x,A%atom(i)%occ,A%atom(i)%Biso,A%atom(i)%moment,A%atom(i)%charge
      if(A%atom(i)%moment > 0.01 ) then
        mag=.true.
        j=i+1
      else
        exit
      end if
      write(unit=lun,fmt="(a)") trim(texto)
   end do
   if(mag) then
      do i=j,A%natoms
       if(A%atom(i)%moment > 0.01 ) then
        write(unit=*,fmt="(a)") " => The magnetic atoms are not on the top of the list"
        write(unit=*,fmt="(a)") " => Please reorder the input file putting on top magnetic atoms"
        stop
       end if
      end do
   else
      write(unit=*,fmt="(a)") " => There is no magnetic atoms in the input list or"
      write(unit=*,fmt="(a)") "    the magnetic atoms are not on the top of the list"
      write(unit=*,fmt="(a)") " => Please reorder the input file putting on top magnetic atoms"
      stop
   end if

   write(unit=*,fmt="(a,f9.3)")" => Maximum bond-distance (Dmax)                       : ",dmax
   write(unit=*,fmt="(a,f9.3)")" => Minimum distance for direct exchange (Direct1)     : ",directex1
   write(unit=*,fmt="(a,f9.3)")" => Maximum distance for direct exchange (Direct2)     : ",directex2
   write(unit=*,fmt="(a,f9.3)")" => Maximum distance for angle calculation(Dangl)      : ",dangl
   write(unit=*,fmt="(a     )")" => Maximum distance for anion-anion bond (Dbond)        "
   write(unit=*,fmt="(a,f9.3)")"              (also for cation-anion bonds)            : ",dbond
   write(unit=*,fmt="(a,f9.3)")" => Maximum angle M-M'^M(M')-A(A') for S-E paths (Angm): ",angm
   write(unit=*,fmt="(a,f9.3)")" => Minimum angle M-A-A'/M'-A'-A   for S-E paths (Angn): ",angn
   write(unit=*,fmt="(a,/)")   "                  (if Dangl=0 no angles are calculated)"
   write(unit=*,fmt="(a  )",advance="no")   "    any change ? (Y/N=<cr>): "

   read(unit=*,fmt="(a)") ans
   if(ans == "Y" .OR. ans == "y") then
     write(unit=*,fmt="(a)",advance="no") " => Give new values for Dmax, Direct1,Direct2, Dangl, Dbond, Angm and Angn : "
     read(unit=*,fmt=*) dmax, directex1,directex2, dangl, dbond, angm, angn
     if(dbond < 2.0 ) then
       write(unit=*,fmt="(a)") " => Warning ! Dbond cannot be lower than 2, Dbond=2.0! "
       dbond=2.0
     end if
     if(directex1 > 0.001) then
       write(unit=*,fmt="(a)",advance="no") " => Give the values for kF (Fermi wavevector) : "
       read(unit=*,fmt=*) kF
     end if
   end if

   write(unit=lun,fmt="(/,/,/,a,f8.4,a)")" => Maximum distances: ",dmax, " for distances calculation"
   write(unit=lun,fmt="(      a,f8.4,a)")"                       ",dangl," for angles calculation"
   write(unit=lun,fmt="(      a,f8.4,a)")"                       ",dbond," for anion-anion/cation-anion bond-distances"

   write(unit=lun,fmt="(a,f8.3)")" => Minimum distance for direct exchange: ", directex1
   write(unit=lun,fmt="(a,f8.3)")" => Maximum distance for direct exchange: ", directex2
   write(unit=lun,fmt="(a,f8.2)")   &
   " => Maximum angle M-M'^M(M')-A(A') for Super(Super)-Exchange paths (Angm):",angm
   write(unit=lun,fmt="(a,f8.2,/,/)")   &
   " => Minimum angle M-A-A'/M'-A'-A   for Super-Super -Exchange paths (Angn):",angn


   ! Generates all atoms in a primitive cell to make calculations in P1
   !-------------------------------------------------------------------

   numops=SpG%numops
   if(SpG%Centred == 2) numops=2*numops
   Call Allocate_Atoms_Cell(nauas,numops,dmax,Ac)
   call Allocate_Coordination_Type(Ac%nat,1,dmax,max_coord)
   write(unit=*,fmt="(a,i6)") " => Calculated maximum coordination: ",max_coord
   call Multi(lun,.true.,.false.,SpG,A,Ac)  !construct Ac,  Multi(Lun,Iprin,Conven,Spg,A,Ac)
   Call deAllocate_atom_list(A)

   call Allocate_Atom_List(Ac%nat,Ap)       !allocate space for Ap
   call Atoms_Cell_to_List(Ac,Ap)
   sp1=SpG%SPG_Symb(1:1)//" 1"
   call Set_SpaceGroup(sp1,gP1)     !construct space group P1/A1/B1/C1/I1/R1/F1
   iprin=.false.
   write(unit=*,fmt="(a)",advance="no")" => List all distances & angles (y/n)? (def=n): "
   read(unit=*,fmt="(a)") ans
   if(ans == "y" .or. ans == "Y") iprin=.true.
   if(iprin) then
      call Calc_Dist_Angle(Dmax, Dangl, Cell, gP1, Ap, Lun)
   else
      call Calc_Dist_Angle(Dmax, Dangl, Cell, gP1, Ap)
   end if
   if(Err_Geom) then
     write(unit=*,fmt="(/,a)") " => WARNING"
     write(unit=*,fmt="(a)") trim(ERR_Geom_Mess)
     write(unit=*,fmt="(a,/)") " => The calculation of exchange paths may be wrong!"
   end if

   ! Extraction of a part (magnetic atoms) of the object Ap by copying into Acm
   if(.not. allocated(ptr)) allocate (ptr(Ap%natoms))
   nmag=0
   !    call Write_Atom_List(Ap)
   do i=1,Ap%natoms
    if(Ap%atom(i)%moment > 0.01) then   !Select magnetic atoms
      nmag=nmag+1
      ptr(nmag)=i
    end if
   end do
   if(nmag == 0) then
     write(unit=*,fmt="(a)") " => No magnetic atoms found!"
     stop
   end if

   ! Look here for Super-Exchange and Super-Super-Exchange paths
   iprin=.false.
   write(unit=*,fmt="(a)",advance="no")" => List detailed exchange paths (y/n)? (def=n): "
   read(unit=*,fmt="(a)") ans
   if(ans == "y" .or. ans == "Y") iprin=.true.
   if(.not. negligible(dbond)) &
       call Exchange_Paths(lun,iprin,dmax,dbond,angm,angn,directex1,directex2,Cell,gP1,Ac,spaths)
   Call deAllocate_Atoms_Cell(Ac)    !From Ac we have conserved only "spaths"
   Call Allocate_Atoms_Cell(nmag, 1, dmax, Acm)

   Acm%nat=nmag
   do i=1,nmag
     Acm%noms(i)  =Ap%atom(ptr(i))%lab
     Acm%xyz (:,i)=Ap%atom(ptr(i))%x(:)
     Acm%moment(i)=Ap%atom(ptr(i))%moment
     Acm%charge(i)=Ap%atom(ptr(i))%charge
   end do
   nif(1)=1
   mom(1)=Acm%moment(1)
   scf(1)="M"//u_case(Ap%atom(ptr(1))%ChemSymb)
   write(scf(1)(4:4),fmt="(i1)") nint(Ap%atom(ptr(1))%charge)
   j=0; L=1
   do i=1,nmag
     j=j+1
     if(abs(Acm%moment(i)-mom(L)) > 0.001) then
        L=L+1
        mom(L)=Acm%moment(i)
        scf(L)="M"//u_case(Ap%atom(ptr(i))%ChemSymb)
        write(scf(L)(4:4),fmt="(i1)") nint(Ap%atom(ptr(i))%charge)
        nif(L)=i
     end if
   end do
   nif(L+1)=nmag+1
   Call deAllocate_Atom_List(Ap)   !Ap is no more needed
   if(allocated(ptr)) deallocate(ptr)

   !  Write terms of the Fourier transform of the exchange interactions
   !  Write also the information to the file *.exc(unit=4) and MCMAG

   open(unit=4,file=trim(outfil)//".exc",status="replace",action="write",position="rewind")
   write(unit=4,fmt="(a)") title
   write(unit=4,fmt="(i4,f8.4)")Acm%nat,dmax
   write(unit=4,fmt="(a)") SpG%SPG_Symb
   write(unit=4,fmt="(3f8.4,3f8.3)") cell%cell,cell%ang
   DO i=1,Acm%nat
     write(unit=4,fmt="(a6,a,4f9.5)")Acm%noms(i)," ",Acm%xyz(:,i),Acm%moment(i)
   END DO
   ! Write file *.mcm for MCMAG
   open(unit=3,file=trim(outfil)//".mcm",status="replace",action="write",position="rewind")
   write(unit=3,fmt="(a)")  title
   write(unit=3,fmt="(a)")" Template File created by program SIMBO for MCMAG (Isotropic Interactions)"
   write(unit=3,fmt="(a)")"!The file should be modified to adapt it to the user needs."
   write(unit=3,fmt="(a)")"!NA(sites)  JCod    Z"
   write(unit=3,fmt="(3i7)")-Acm%nat,0,Acm%nat
   !Calculation for McMag
   Select Case(SpG%SPG_Symb(1:1))
      case("A")
        M_coor=M_A
      case("B")
        M_coor=M_B
      case("C")
        M_coor=M_C
      case("I")
        M_coor=M_I
      case("R")
        M_coor=M_R
      case("F")
        M_coor=M_F
      case default
        M_coor=reshape((/ 1.0, 0.0, 0.0, &
                          0.0, 1.0, 0.0, &
                          0.0, 0.0, 1.0 /),(/3,3/))
   End Select
   M_basis=transpose(invert_a(M_coor))
   call Change_Setting_Cell(Cell,M_basis,Celln)
   !End of calculation for McMag

   Call construct_jxch(lun,nmag,spaths,Acm)

   write(unit=3,fmt="(a)")  "!   Ni    Nf   Spin     ScattFact"
   do i=1,L
      write(unit=3,fmt="(2i6,f8.4,tr6,a)") nif(i),nif(i+1)-1,mom(i),scf(i)
   end do
   write(unit=3,fmt="(a)") "!    Primitive unit cell  "
   write(unit=3,fmt="(a)") "!     a          b          c        alpha       beta      gamma   "

   write(unit=3,fmt="(6f11.5)")celln%cell,celln%ang
   write(unit=3,fmt="(a)") "!"
   write(unit=3,fmt="(a)") "!  The conditions below should be adapted to the problem by the user "
   write(unit=3,fmt="(a)") "!"
   write(unit=3,fmt="(a)") "SpinModel    Heisenberg "
   write(unit=3,fmt="(a)") "  "
   write(unit=3,fmt="(a)") "Title  Simulation of classical Spins:"//trim(title)
   write(unit=3,fmt="(a)") "  "
   write(unit=3,fmt="(a)") "!  Simulation box "
   write(unit=3,fmt="(a)") "Ncells    3 3 3   "
   write(unit=3,fmt="(a)") "  "
   write(unit=3,fmt="(a)") "!  Initial configuration (R,I) "
   write(unit=3,fmt="(a)") "InitConf  R "
   write(unit=3,fmt="(a)") "  "
   write(unit=3,fmt="(a)") "! boundary conditions (Free,Periodic,Mixed)"
   write(unit=3,fmt="(a)") "Boundary  Periodic"
   write(unit=3,fmt="(a)") "  "
   write(unit=3,fmt="(a)") "! Scaling (sample,cell,site,mole)"
   write(unit=3,fmt="(a)") "Scale     cell"
   write(unit=3,fmt="(a)") "  "
   write(unit=3,fmt="(a)") "!  Sites for output during simulation"
   write(unit=3,fmt="(a)") "Sites   1 2 3"
   write(unit=3,fmt="(a)") "  "
   write(unit=3,fmt="(a)") "!         T_ini   Coef  T_final"
   write(unit=3,fmt="(a)") "schedule    500   0.95    0.5"
   write(unit=3,fmt="(a)") "  "
   write(unit=3,fmt="(a)") "!  Magnetic Field"
   write(unit=3,fmt="(a)") "hfield    0  0  0  1"
   write(unit=3,fmt="(a)") "  "
   write(unit=3,fmt="(a)") "!  Number of MC cycles and thermalization"
   write(unit=3,fmt="(a)") "mcyc   5000  500"
   write(unit=3,fmt="(a)") "  "
   write(unit=3,fmt="(a)") "print  E"
   write(unit=3,fmt="(a)") "averages"
   write(unit=3,fmt="(a)") "cryst   1  1 0 0"

   close(unit=3)

   Call deAllocate_Atoms_Cell(Acm)

   if(allocated(file_dat)) deallocate(file_dat)

   CALL CPU_TIME(End_time)
   seconds=End_time-start_time
   rminutes=seconds/60.0
   hours=rminutes/60.0
   write(unit=*,fmt="(/,a,/  )") "                  Total CPU-Time "
   write(unit=*,fmt="(a,f14.2)") "             CPU-seconds: ",seconds
   write(unit=*,fmt="(a,f14.2)") "             CPU-minutes: ",rminutes
   write(unit=*,fmt="(a,f14.2)") "             CPU-hours  : ",hours
   write(unit=*,fmt="(/,a)")" => Results in files:"
   write(unit=*,fmt="(2a)") "                     ",outfil(1:lr)//".nei"
   write(unit=*,fmt="(2a)") "                     ",outfil(1:lr)//".exc -> input for ENERMAG"
   write(unit=*,fmt="(2a)") "                     ",outfil(1:lr)//".mcm -> input for MCMAG"

   stop

  End Program Simbo
