module rng_c
  use v3d_func_rep
  implicit none
  private
  public :: rng_init_seed, rng_get_rotmat, rng_get_vec_ss, rng_get_double

contains
  subroutine rng_init_seed(rand_flag)
    !~ Initializes the rng depending on rand_flag with
    !~ either a seed based on time and pid or the seed
    !~ used last or the machine default seed
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, dt(8), pid, t(2), s, ierr
    integer(8) :: count, tms
    logical :: rand_flag

    call random_seed(size=n)
    allocate(seed(n))

    if (rand_flag .eqv. .true.) then
     ! XOR:ing the current time and pid. The PID is
     ! useful in case one launches multiple instances of the same
     ! program in parallel.
      call system_clock(count)
      if (count /= 0) then
        t = transfer(count, t)
      else
        call date_and_time(values=dt)
        tms = (dt(1) - 1980) * 365_8 * 24 * 60 * 60 * 1000 &
            + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
            + dt(3) * 24 * 60 * 60 * 1000 &
            - dt(4) * 60 * 1000 + dt(5) * 60 * 60 * 1000 &
            + dt(6) * 60 * 1000 + dt(7) * 1000 &
            + dt(8)
        t = transfer(tms, t)
      end if
      s = ieor(t(1), t(2))
      pid = getpid() + 1099279 ! Add a prime
      s = ieor(s, pid)
      if (n >= 3) then
        seed(1) = t(1) + 36269
        seed(2) = t(2) + 72551
        seed(3) = pid
        if (n > 3) then
          seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
        end if
      else
        seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
      end if
    else
      open(404,file='.rng_last_seed',form='unformatted',status='old',iostat=ierr)
      if (ierr == 0) then
        read(404) seed
        close(404)
      else
        call random_seed(get=seed)
      end if
    end if
    open(405,file='.rng_last_seed',form='unformatted',status='replace')
    write(405) seed
    close(405)
    call random_seed(put=seed)
    deallocate(seed)
  end subroutine
  
  function rng_get_rotmat() result(mat)
    !~ Creates a random rotation matrix
    double precision, dimension(4) :: x
    double precision, dimension(3,3) :: mat
    
    x = get_po3sph(get_ptoncir(),get_ptoncir())
    
    mat(1,1) = 1-2*(x(2)**2+x(3)**2)
    mat(2,2) = 1-2*(x(1)**2+x(3)**2)
    mat(3,3) = 1-2*(x(1)**2+x(2)**2)
    mat(1,2) = 2*(x(1)*x(2)-x(3)*x(4))
    mat(1,3) = 2*(x(1)*x(3)+x(2)*x(4))
    mat(2,3) = 2*(x(2)*x(3)-x(1)*x(4))
    mat(2,1) = 2*(x(1)*x(2)+x(3)*x(4))
    mat(3,1) = 2*(x(1)*x(3)-x(2)*x(4))
    mat(3,2) = 2*(x(2)*x(3)+x(1)*x(4))
  end function
  
  function get_po3sph(p1,p2) result(x)
    !~ Used by rng_get_rotmat
    !~ Calculates a point on the S3 using two given points
    !~ inside the unit circle
    !~ Parameters:
    !~ p1, p2 ... points on the unit circle, x ... point on S3
    double precision, dimension(2) :: p1, p2
    double precision, dimension(4) :: x
    double precision :: c
    
    c=dsqrt((1-p1(1)**2-p1(2)**2)/(p2(1)**2+p2(2)**2))
    
    x(1:2) = p1
    x(3) = p2(1)*c
    x(4) = p2(2)*c
  end function
  
  function get_ptoncir() result(p)
    !~ Used by rng_get_rotmat
    !~ Calculates a random point inside the unit circle
    double precision, dimension(2) :: p
    
    do
      call random_number(p(1))
      call random_number(p(2))
      p = p*2-1.d0
      if (p(1)**2+p(2)**2 .lt. 1.d0) exit
    end do
  end function
  
  function rng_get_vec_ss(minr,maxr) result(vec)
    !~ Calculates a random point on a spherical shell defined
    !~ by the upper and lower radius
    double precision :: minr, maxr, minb, nrm
    double precision, dimension(3) :: vec
    integer :: i
    
    if (abs(maxr/minr-1.d0) .lt. 1.d-1) then
      vec(:) = rng_get_double(minr,maxr)*get_lin_map(rng_get_rotmat(),(/1.d0,0.d0,0.d0/))
    else
      minb = minr/dsqrt(3.d0)
      do
        do i=1,3
          vec(i) = rng_get_double(-maxr,maxr)
        end do
        nrm = get_norm(vec)
        if (nrm .ge. minr .and. nrm .le. maxr) exit
      end do
    end if
  end function
  
  function rng_get_double(minx,maxx) result(x)
    double precision :: minx, maxx, x, tmp
    
    call random_number(tmp)
    x = minx + (maxx - minx)*tmp
  
  end function
  
end module
