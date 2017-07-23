module aenettinker

  !--------------------------------------------------------------------!
  !                      aenet - Tinker Interface                      !
  !                                                                    !
  ! Use the atomic energy network (aenet) library for machine-learning !
  ! potentials with Tinker.                                            !
  !                                                                    !
  ! To compile this module and link it with the Tinker lib, the module !
  ! file 'aenet.mod' from the aenet distribution needs to be present.  !
  !                                                                    !
  ! To work, the following keywords are required in the Tinker input   !
  ! file:                                                              !
  !       VDW-LIST         The vdW neighbor list is used to keep track !
  !                        of atomic environments.                     !
  !       VDW-CUTOFF <R>   To set the cutoff distance of the potential !
  !       EXTRATERM only   To exclusively use the aenet energy term    !
  !--------------------------------------------------------------------!

  use aenet,  only: aenet_init,                     &
                    aenet_final,                    &
                    aenet_load_potential,           &
                    aenet_atomic_energy,            &
                    aenet_atomic_energy_and_forces, &
                    aenet_print_info,               &
                    aenet_nnb_max,                  &
                    aenet_nbl_init,                 &
                    aenet_nbl_final,                &
                    aenet_nbl_neighbors,            &
                    AENET_OK, AENET_TRUE

  use atoms,  only: tinker_n => n,                  &
                    tinker_type => type,            &
                    tinker_x => x,                  &
                    tinker_y => y,                  &
                    tinker_z => z
  use atomid, only: tinker_name => name
  use boxes,  only: tinker_lvec => lvec
  use deriv,  only: tinker_dex => dex
  use energi, only: tinker_ex => ex
  use sizes,  only: tinker_maxtyp => maxtyp
  use units,  only: evolt, hartree

  implicit none
  private
  save

  public :: aenet_tinker_init, aenet_eval_energy, aenet_eval_energy_and_forces

  double precision, parameter, private :: eV2kcal = hartree/evolt

  logical,                                       private :: isInit = .false.
  integer,                                       private :: ntypes
  character(len=3), dimension(tinker_maxtyp),    private :: atom_types
  integer,          dimension(tinker_maxtyp),    private :: t2a_type_ID

  integer,          dimension(:),   allocatable, private :: at_type_ID
  double precision, dimension(:,:), allocatable, private :: at_coo_cart
  double precision, dimension(:,:), allocatable, private :: at_for_cart
  double precision, dimension(:,:), allocatable, private :: at_for_local

  double precision, dimension(:,:), allocatable, private :: at_nbcoo
  double precision, dimension(:),   allocatable, private :: at_nbdist
  integer,          dimension(:),   allocatable, private :: at_nbtype
  integer,          dimension(:),   allocatable, private :: at_nblist
contains

  !--------------------------------------------------------------------!
  !                      Initialize the module                         !
  ! The init subroutine can be called arbitrary times.  It will just   !
  ! return without doing anything after the first call.                !
  !                                                                    !
  ! For each atom type T an ANN potential file 'T.ann' is expected in  !
  ! the working directory.                                             !
  !--------------------------------------------------------------------!

  subroutine aenet_tinker_init()

    implicit none

    integer :: stat

    character(len=7) :: pot_file
    integer :: iatom, itype, i

    if (isInit) return

    ! determine number and name of atom types
    ntypes = 0
    do iatom = 1, tinker_n
       if (.not. in(tinker_name(iatom), atom_types(1:ntypes))) then
          ntypes = ntypes + 1
          atom_types(ntypes) = tinker_name(iatom)
       end if
    end do

    ! index for aenet type ID from Tinker type ID
    do itype = 1, ntypes
       iatom = 1
       do while (tinker_name(iatom) /= atom_types(itype))
          iatom = iatom + 1
       end do
       t2a_type_ID(tinker_type(iatom)) = itype
    end do

    call aenet_init(atom_types(1:ntypes), stat)
    if (stat /= AENET_OK) then
       write(0, *) 'Error: aenet module could not be initialized.'
       stop
    end if

    ! load ANN potentials of all atom types
    do itype = 1, ntypes
       pot_file = ""
       pot_file = trim(adjustl(atom_types(itype))) // '.ann'
       call aenet_load_potential(itype, trim(pot_file), stat)
       if (stat /= AENET_OK) then
          write(0, '(1x, A, I3)') &
               "Error: could not load ANN potential for type: ", itype
          write(0, '(1x, "       File: ", A)') trim(pot_file)
          stop
       end if
    end do

    call aenet_print_info()

    ! allocate memory
    allocate(at_coo_cart(3, tinker_n), &
             at_for_cart(3, tinker_n), &
             at_for_local(3, tinker_n), &
             at_type_ID(tinker_n), &
             at_nbcoo(3, aenet_nnb_max),  &
             at_nbdist(aenet_nnb_max), &
             at_nbtype(aenet_nnb_max), &
             at_nblist(aenet_nnb_max))

    ! atom type ID of each atom
    do iatom = 1, tinker_n
       at_type_ID(iatom) = t2a_type_ID(tinker_type(iatom))
    end do

    isInit = .true.

  end subroutine aenet_tinker_init

  !--------------------------------------------------------------------!

  subroutine aenet_tinker_final()

    implicit none

    integer :: stat

    if (.not. isInit) return

    deallocate(at_coo_cart, at_for_cart, at_for_local, at_type_ID, &
               at_nbcoo, at_nbdist, at_nbtype, at_nblist)

    call aenet_final(stat)
    if (stat /= AENET_OK) then
       write(0,*) 'Error: Unable to finalize aenet module.'
       return
    end if

  end subroutine aenet_tinker_final

  !--------------------------------------------------------------------!
  !                          Energy evaluation                         !
  !--------------------------------------------------------------------!

  subroutine aenet_eval_energy()

    implicit none

    double precision                      :: E, E_i
    double precision, dimension(3,3)      :: latt_vec
    double precision, dimension(3)        :: coo_i
    integer                               :: type_i
    integer                               :: iat1, iat2
    integer                               :: i, j, nnb
    integer                               :: stat

    tinker_ex = 0.0d0

    ! updated lattice vectors
    ! TODO: only do this in variable-cell simulations
    latt_vec(1:3, 1) = [tinker_lvec(1,1), tinker_lvec(1,2), tinker_lvec(1,3)]
    latt_vec(1:3, 2) = [tinker_lvec(2,1), tinker_lvec(2,2), tinker_lvec(2,3)]
    latt_vec(1:3, 3) = [tinker_lvec(3,1), tinker_lvec(3,2), tinker_lvec(3,3)]

    ! updated coordinates array
    do i=1, tinker_n
       at_coo_cart(1:3, i) = [tinker_x(i), tinker_y(i), tinker_z(i)]
    end do

    ! construct neighbor list
    call aenet_nbl_init(&
         latt_vec, tinker_n, at_type_ID, at_coo_cart, AENET_TRUE, AENET_TRUE)

    E = 0.0d0
!$OMP PARALLEL default(shared) &
!$OMP private(iat1, coo_i, type_i, nnb, at_nbcoo, at_nbdist, at_nblist, &
!$OMP         at_nbtype, stat, E_i)

    !$OMP DO schedule(guided) reduction(+:E)
    do iat1 = 1, tinker_n
       coo_i(1:3) = at_coo_cart(1:3, iat1)
       type_i = t2a_type_ID(tinker_type(iat1))
       nnb = aenet_nnb_max
       call aenet_nbl_neighbors(&
            iat1, nnb, at_nbcoo, at_nbdist, at_nblist, at_nbtype)
       call aenet_atomic_energy(&
            coo_i, type_i, nnb, at_nbcoo, at_nbtype, E_i, stat)
       E = E + E_i
       if (stat /= AENET_OK) then
          write(0, *) 'Error: evaluation of atomic energy failed.'
          stop
       end if
    end do
    !$OMP END DO

!$OMP END PARALLEL

    call aenet_nbl_final()

    tinker_ex = E*eV2kcal

  end subroutine aenet_eval_energy

  !--------------------------------------------------------------------!

  subroutine aenet_eval_energy_and_forces()

    implicit none

    ! double precision, dimension(tinker_n) :: E_i
    double precision                      :: E, E_i
    double precision, dimension(3,3)      :: latt_vec
    double precision, dimension(3)        :: coo_i
    integer                               :: type_i
    integer                               :: iat1, iat2
    integer                               :: i, j, nnb
    integer                               :: stat

    tinker_ex = 0.0d0
    tinker_dex = 0.0d0

    ! updated lattice vectors
    ! TODO: only do this in variable-cell simulations
    latt_vec(1:3, 1) = [tinker_lvec(1,1), tinker_lvec(1,2), tinker_lvec(1,3)]
    latt_vec(1:3, 2) = [tinker_lvec(2,1), tinker_lvec(2,2), tinker_lvec(2,3)]
    latt_vec(1:3, 3) = [tinker_lvec(3,1), tinker_lvec(3,2), tinker_lvec(3,3)]

    ! updated coordinates array
    do i = 1, tinker_n
       at_coo_cart(1:3, i) = [tinker_x(i), tinker_y(i), tinker_z(i)]
    end do

    ! construct neighbor list
    call aenet_nbl_init(&
         latt_vec, tinker_n, at_type_ID, at_coo_cart, AENET_TRUE, AENET_TRUE)

    E = 0.0d0
    at_for_cart(:,:) = 0.0d0
!$OMP PARALLEL default(shared) &
!$OMP private(iat1, coo_i, type_i, nnb, at_nbcoo, at_nbdist, at_nblist, &
!$OMP         at_nbtype, stat, at_for_local, E_i)
    at_for_local = 0.0d0
    !$OMP DO schedule(guided) reduction(+:E)
    do iat1 = 1, tinker_n
       coo_i(1:3) = at_coo_cart(1:3, iat1)
       type_i = t2a_type_ID(tinker_type(iat1))
       nnb = aenet_nnb_max
       call aenet_nbl_neighbors(&
            iat1, nnb, at_nbcoo, at_nbdist, at_nblist, at_nbtype)
       call aenet_atomic_energy_and_forces( &
            coo_i, type_i, iat1, nnb, at_nbcoo, at_nbtype, at_nblist, &
            tinker_n, E_i, at_for_local, stat)
       E = E + E_i
       if (stat /= AENET_OK) then
          write(0, *) 'Error: evaluation of atomic energy and forces failed.'
          stop
       end if
    end do
    !$OMP END DO
    !$OMP CRITICAL
    at_for_cart(:,:) = at_for_cart(:,:) + at_for_local(:,:)
    !$OMP END CRITICAL
!$OMP END PARALLEL

    call aenet_nbl_final()

    tinker_ex = E*eV2kcal
    tinker_dex = -at_for_cart*eV2kcal

  end subroutine aenet_eval_energy_and_forces

  !--------------------------------------------------------------------!
  !                        auxiliary functions                         !
  !--------------------------------------------------------------------!

  function in(record, list) result(in_list)

    implicit none

    character(len=*),               intent(in) :: record
    character(len=*), dimension(:), intent(in) :: list
    logical                                    :: in_list

    integer :: i

    in_list = .false.
    do i = 1, size(list)
       if (trim(record) == trim(list(i))) then
          in_list = .true.
          return
       end if
    end do

  end function in

end module aenettinker
