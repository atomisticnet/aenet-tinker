module aenettinker

  !--------------------------------------------------------------------!
  !                      aenet - Tinker Interface                      !
  !                                                                    !
  ! Use the atomic energy network (aenet) library for machine-learning !
  ! potentials with Tinker.                                            !
  !                                                                    !
  ! To make use of this interface, the accompanying files 'extra.f'    !
  ! and 'extra1.f' need to be copied to Tinker's 'source' directory.   !
  !                                                                    !
  ! To work, the following keywords are required in the Tinker input   !
  ! file:                                                              !
  !                                                                    !
  !       EXTRATERM only   To exclusively use the aenet energy term    !
  !                                                                    !
  ! The number of threads used for OpenMP parallelization can be set   !
  ! using the 'OPENMP-THREADS' keyword.                                !
  !                                                                    !
  ! NOTE: Currently, only the GNU Fortran compiler is compatible with  !
  !       the OpenMP parallelized section.  Intel Fortran compiler 15  !
  !       gives WRONG results!
  !--------------------------------------------------------------------!
  ! 2017-07-23 Nongnuch Artrith and Alexander Urban                    !
  !--------------------------------------------------------------------!

  use, intrinsic :: iso_c_binding, only: c_bool

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
                    AENET_OK, AENET_TRUE, AENET_FALSE

  use atoms,  only: tinker_n => n,                  &
                    tinker_type => type,            &
                    tinker_x => x,                  &
                    tinker_y => y,                  &
                    tinker_z => z
  use atomid, only: tinker_name => name
  use bound,  only: tinker_use_bounds => use_bounds
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

  logical(kind=c_bool),                          private :: at_pbc

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

    ! periodic boundary conditions?
    if (tinker_use_bounds) then
       at_pbc = AENET_TRUE
    else
       at_pbc = AENET_FALSE
    end if

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

    ! updated coordinates array
    do i=1, tinker_n
       at_coo_cart(1:3, i) = [tinker_x(i), tinker_y(i), tinker_z(i)]
    end do

    if (at_pbc) then
       ! updated lattice vectors
       ! TODO: only do this in variable-cell simulations
       latt_vec(1:3,1) = [tinker_lvec(1,1), tinker_lvec(1,2), tinker_lvec(1,3)]
       latt_vec(1:3,2) = [tinker_lvec(2,1), tinker_lvec(2,2), tinker_lvec(2,3)]
       latt_vec(1:3,3) = [tinker_lvec(3,1), tinker_lvec(3,2), tinker_lvec(3,3)]
    else
       ! Determine bounding box of isolated system. Also shift all atoms
       ! to the bounding box.
       call bounding_box(at_coo_cart, latt_vec)
    end if

    ! construct neighbor list
    call aenet_nbl_init(&
         latt_vec, tinker_n, at_type_ID, at_coo_cart, AENET_TRUE, at_pbc)

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

    ! updated coordinates array
    do i = 1, tinker_n
       at_coo_cart(1:3, i) = [tinker_x(i), tinker_y(i), tinker_z(i)]
    end do

    if (at_pbc) then
       ! updated lattice vectors
       ! TODO: only do this in variable-cell simulations
       latt_vec(1:3,1) = [tinker_lvec(1,1), tinker_lvec(1,2), tinker_lvec(1,3)]
       latt_vec(1:3,2) = [tinker_lvec(2,1), tinker_lvec(2,2), tinker_lvec(2,3)]
       latt_vec(1:3,3) = [tinker_lvec(3,1), tinker_lvec(3,2), tinker_lvec(3,3)]
    else
       ! Determine bounding box of isolated system. Also shift all atoms
       ! to the bounding box.
       call bounding_box(at_coo_cart, latt_vec)
    end if

    ! construct neighbor list
    call aenet_nbl_init(&
         latt_vec, tinker_n, at_type_ID, at_coo_cart, AENET_TRUE, at_pbc)

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

! DEBUG
!   write(*,*) '**** aenet **** Energy (eV) = ', E
! END DEBUG

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

  !--------------------------------------------------------------------!
  !             bounding box for non-periodic simulations              !
  !                                                                    !
  ! This routine is essentially copied from aenet's geometry.f90       !
  ! Will consider making it available through the aenet API in later   !
  ! versions.                                                          !
  !--------------------------------------------------------------------!

  subroutine bounding_box(cooCart, box)

    implicit none

    double precision, dimension(:,:), intent(inout) :: cooCart
    double precision, dimension(3,3), intent(out)   :: box

    double precision, dimension(3) :: shift

    double precision :: x_min, x_max
    double precision :: y_min, y_max
    double precision :: z_min, z_max

    integer :: iat

    ! To avoid numerical problems we make the bounding box
    ! slightly larger than necessary and also consider this in
    ! the shift of the coordinates. This should guarantee that all
    ! scaled coordinates are in [0,1).
    double precision, parameter :: EPS = 2.0d-6

    x_min = minval(cooCart(1,:))
    x_max = maxval(cooCart(1,:))
    y_min = minval(cooCart(2,:))
    y_max = maxval(cooCart(2,:))
    z_min = minval(cooCart(3,:))
    z_max = maxval(cooCart(3,:))

    ! orthogonal bounding box
    box(:,1) = (/ x_max - x_min + EPS, 0.0d0, 0.0d0  /)
    box(:,2) = (/ 0.0d0, y_max - y_min + EPS, 0.0d0  /)
    box(:,3) = (/ 0.0d0, 0.0d0, z_max - z_min + EPS  /)

    ! origin of the bounding box
    shift(1) = x_min - 0.5d0*EPS
    shift(2) = y_min - 0.5d0*EPS
    shift(3) = z_min - 0.5d0*EPS

    ! shift coordinates to bounding box
    do iat = 1, size(cooCart(1,:))
       cooCart(1:3,iat) = cooCart(1:3,iat) - shift(1:3)
    end do

  end subroutine bounding_box

end module aenettinker
