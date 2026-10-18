module irmsd_refine
  !***************************************************************
  !* Bond-consistent refinement of an iRMSD atom assignment.
  !* Terminal (degree-1) atoms are re-assigned only among the
  !* same-rank terminal neighbours of their matched parent, in
  !* alternation with a Kabsch alignment of all atoms, until the
  !* assignment no longer changes.
  !***************************************************************
  use,intrinsic :: iso_fortran_env,only:int64
  use crest_parameters
  use strucrd
  use irmsd_module,only:rmsd
  implicit none
  private

  public :: refine_terminal_assignment

  integer,parameter :: maxiter = 50
  !> number of starting assignments for the alternating loop
  integer,parameter :: nstart = 4
  !> largest sibling group solved by full enumeration (8! permutations)
  integer,parameter :: maxbrute = 8

contains

  subroutine refine_terminal_assignment(ref,mol,rank,niter,io)
    !*******************************************************************
    !* Re-assign the terminal atoms of mol among the same-rank terminal
    !* neighbours of their matched parent ("sibling groups"), alternating
    !* with a Kabsch alignment on all atoms until the assignment is
    !* stable. Kabsch minimizes the RMSD for a fixed assignment and the
    !* assignment step minimizes it for a fixed frame, so each run ends
    !* in a local minimum. The loop is started from the incoming mapping,
    !* from all groups reversed, and from seeded random assignments; the
    !* lowest result is kept. Non-terminal atoms keep their mapping.
    !*
    !* ref   : reference structure
    !* mol   : structure in matched order (mol atom i <-> ref atom i);
    !*         reordered and aligned onto ref on output
    !* rank  : canonical rank of each matched atom pair
    !* niter : (out) total number of assignment updates over all starts
    !* io    : (out) 0 on success; 1 if the terminal atoms around some
    !*         parent differ between the two structures, in which case
    !*         mol is returned unchanged
    !*******************************************************************
    implicit none
    type(coord),intent(in) :: ref
    type(coord),intent(inout) :: mol
    integer,intent(in) :: rank(:)
    integer,intent(out) :: niter
    integer,intent(out) :: io

    type(coord) :: work
    integer,allocatable :: gstart(:),gref(:),gmol(:)
    integer,allocatable :: perm(:),bestperm(:),used(:)
    real(wp) :: val,bestval,U(3,3)
    integer :: nat,ngroup,istart,g,n,j,iter
    integer(int64) :: seed

    io = 0
    niter = 0
    nat = ref%nat

    call sibling_groups(ref,mol,rank,ngroup,gstart,gref,gmol,io)
    if (io /= 0) return
    if (ngroup == 0) return

    allocate (perm(nat),bestperm(nat),used(nat))
    work = mol
    bestval = huge(bestval)
    seed = 1

    STARTS: do istart = 1,nstart
!>-- starting assignment
      perm(:) = [(j,j=1,nat)]
      if (istart > 1) then
        do g = 1,ngroup
          n = gstart(g+1)-gstart(g)
          call start_order(istart,n,gmol(gstart(g):gstart(g+1)-1), &
          &                perm,gref(gstart(g):gstart(g+1)-1),seed)
        end do
      end if

!>-- the groups must combine into a permutation of all atoms
      used(:) = 0
      do j = 1,nat
        used(perm(j)) = used(perm(j))+1
      end do
      if (any(used /= 1)) then
        io = 1
        return
      end if

      call alternate(ref,mol,work,ngroup,gstart,gref,gmol,perm,iter)
      niter = niter+iter

      work%xyz(:,:) = mol%xyz(:,perm)
      val = rmsd(ref,work)
      if (val < bestval) then
        bestval = val
        bestperm(:) = perm(:)
      end if
    end do STARTS

!>-- apply the best assignment and alignment
    mol%xyz(:,:) = mol%xyz(:,bestperm)
    mol%at(:) = mol%at(bestperm)
    val = rmsd(ref,mol,rotmat=U)
    mol%xyz = matmul(U,mol%xyz)
  end subroutine refine_terminal_assignment

!========================================================================================!

  subroutine sibling_groups(ref,mol,rank,ngroup,gstart,gref,gmol,io)
    !*****************************************************************
    !* Collect the sibling groups: for every non-terminal atom p, the
    !* terminal neighbours of a given rank in ref (gref) and the
    !* terminal neighbours of the same rank of mol atom p (gmol).
    !* Group g occupies positions gstart(g):gstart(g+1)-1.
    !* io = 1 if a group has a different size in the two structures.
    !*****************************************************************
    implicit none
    type(coord),intent(in) :: ref,mol
    integer,intent(in) :: rank(:)
    integer,intent(out) :: ngroup
    integer,allocatable,intent(out) :: gstart(:),gref(:),gmol(:)
    integer,intent(out) :: io
    integer,allocatable :: aref(:,:),amol(:,:),degref(:),degmol(:)
    integer :: nat,p,t,j,n,k

    io = 0
    nat = ref%nat
    call ref%cn_to_adjacency(aref)
    call mol%cn_to_adjacency(amol)
    degref = sum(aref,dim=1)
    degmol = sum(amol,dim=1)

    allocate (gstart(nat+1),gref(nat),gmol(nat))
    ngroup = 0
    gstart(1) = 1
    do p = 1,nat
      if (degref(p) < 2) cycle
      do t = 1,nat
        if (.not.is_child(aref,degref,p,t)) cycle
        !> handle each rank once per parent, at its first child
        if (any([(is_child(aref,degref,p,j).and.rank(j) == rank(t),j=1,t-1)])) cycle
        n = gstart(ngroup+1)-1
        k = n
        do j = 1,nat
          if (rank(j) /= rank(t)) cycle
          if (is_child(aref,degref,p,j)) then
            n = n+1
            gref(n) = j
          end if
          if (is_child(amol,degmol,p,j)) then
            k = k+1
            if (k > nat) then
              io = 1
              return
            end if
            gmol(k) = j
          end if
        end do
        if (k /= n) then
          io = 1
          return
        end if
        ngroup = ngroup+1
        gstart(ngroup+1) = n+1
      end do
    end do
  end subroutine sibling_groups

!========================================================================================!

  subroutine start_order(istart,n,mm,perm,rr,seed)
    !***********************************************************
    !* Starting assignment of one sibling group: reversed order
    !* for istart = 2, a seeded random permutation otherwise.
    !***********************************************************
    implicit none
    integer,intent(in) :: istart,n
    integer,intent(in) :: mm(n),rr(n)
    integer,intent(inout) :: perm(:)
    integer(int64),intent(inout) :: seed
    integer :: s(n),a,b,tmp

    s(:) = [(a,a=1,n)]
    if (istart == 2) then
      s(:) = s(n:1:-1)
    else
      do a = n,2,-1
        seed = mod(1103515245_int64*seed+12345_int64,2147483647_int64)
        b = int(mod(seed,int(a,int64)))+1
        tmp = s(a); s(a) = s(b); s(b) = tmp
      end do
    end if
    do a = 1,n
      perm(rr(a)) = mm(s(a))
    end do
  end subroutine start_order

!========================================================================================!

  subroutine alternate(ref,mol,work,ngroup,gstart,gref,gmol,perm,niter)
    !***********************************************************
    !* Alternate Kabsch alignment and sibling-group assignment
    !* from the starting assignment perm until it is stable.
    !***********************************************************
    implicit none
    type(coord),intent(in) :: ref,mol
    type(coord),intent(inout) :: work
    integer,intent(in) :: ngroup
    integer,intent(in) :: gstart(:),gref(:),gmol(:)
    integer,intent(inout) :: perm(:)
    integer,intent(out) :: niter
    integer,allocatable :: newperm(:),sigma(:)
    real(wp),allocatable :: aligned(:,:),cost(:,:)
    real(wp) :: U(3,3),cm(3),cr(3),dum
    integer :: nat,g,i0,n,a,b,j,iter

    nat = ref%nat
    niter = 0
    allocate (newperm(nat),sigma(nat),aligned(3,nat))
    cm(:) = sum(mol%xyz,dim=2)/real(nat,wp)
    cr(:) = sum(ref%xyz,dim=2)/real(nat,wp)

    do iter = 1,maxiter
      work%xyz(:,:) = mol%xyz(:,perm)
      dum = rmsd(ref,work,rotmat=U)
      do j = 1,nat
        aligned(:,j) = matmul(U,mol%xyz(:,j)-cm)+cr
      end do

      newperm(:) = perm(:)
      do g = 1,ngroup
        i0 = gstart(g)-1
        n = gstart(g+1)-gstart(g)
        allocate (cost(n,n))
        do a = 1,n
          do b = 1,n
            cost(a,b) = sum((ref%xyz(:,gref(i0+a))-aligned(:,gmol(i0+b)))**2)
          end do
        end do
        call best_assignment(n,cost,sigma)
        deallocate (cost)
        do a = 1,n
          newperm(gref(i0+a)) = gmol(i0+sigma(a))
        end do
      end do

      if (all(newperm == perm)) exit
      perm(:) = newperm(:)
      niter = niter+1
    end do
  end subroutine alternate

!========================================================================================!

  pure logical function is_child(A,deg,p,j)
    !*************************************************
    !* True if j is a terminal (degree-1) neighbour
    !* of p in the graph with adjacency A.
    !*************************************************
    implicit none
    integer,intent(in) :: A(:,:),deg(:)
    integer,intent(in) :: p,j
    is_child = (A(j,p) /= 0).and.(deg(j) == 1)
  end function is_child

!========================================================================================!

  subroutine best_assignment(n,cost,sigma)
    !***********************************************************
    !* Assignment a -> sigma(a) minimizing sum cost(a,sigma(a)).
    !* Full enumeration (Heap's algorithm) up to maxbrute atoms,
    !* greedy beyond that; sibling groups are almost always
    !* 2 or 3 atoms (CH2, CH3, CF3).
    !***********************************************************
    implicit none
    integer,intent(in) :: n
    real(wp),intent(in) :: cost(n,n)
    integer,intent(out) :: sigma(:)
    integer :: s(n),c(n),i,a,b,tmp
    real(wp) :: val,best
    logical :: taken(n)

    sigma(1:n) = [(i,i=1,n)]
    if (n < 2) return

    if (n > maxbrute) then
      taken(:) = .false.
      do a = 1,n
        best = huge(best)
        do b = 1,n
          if (.not.taken(b).and.cost(a,b) < best) then
            best = cost(a,b)
            sigma(a) = b
          end if
        end do
        taken(sigma(a)) = .true.
      end do
      return
    end if

    s(:) = [(i,i=1,n)]
    c(:) = 1
    best = sum([(cost(a,s(a)),a=1,n)])
    i = 1
    do while (i <= n)
      if (c(i) < i) then
        if (mod(i,2) == 1) then
          tmp = s(1); s(1) = s(i); s(i) = tmp
        else
          tmp = s(c(i)); s(c(i)) = s(i); s(i) = tmp
        end if
        val = sum([(cost(a,s(a)),a=1,n)])
        if (val < best) then
          best = val
          sigma(1:n) = s(:)
        end if
        c(i) = c(i)+1
        i = 1
      else
        c(i) = 1
        i = i+1
      end if
    end do
  end subroutine best_assignment

end module irmsd_refine
