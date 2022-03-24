program test_bobyqa

    use kind_module
    use bobyqa_module
    
    implicit none

!*****************************************************************************************
!>
!  Test problem for [[bobyqa]], the objective function being the sum of
!  the reciprocals of all pairwise distances between the points P_I,
!  I=1,2,...,M in two dimensions, where M=N/2 and where the components
!  of P_I are X(2*I-1) and X(2*I). Thus each vector X of N variables
!  defines the M points P_I. The initial X gives equally spaced points
!  on a circle. Four different choices of the pairs (N,NPT) are tried,
!  namely (10,16), (10,21), (20,26) and (20,41). Convergence to a local
!  minimum that is not global occurs in both the N=10 cases. The details
!  of the results are highly sensitive to computer rounding errors. The
!  choice IPRINT=2 provides the current X and optimal F so far whenever
!  RHO is reduced. The bound constraints of the problem require every
!  component of X to be in the interval [-1,1].

    ! subroutine bobyqa_test()

     !   implicit none
        
        real(wp),dimension(100) :: x, xl, xu        
        integer :: i,j,m,n,jcase,npt
        real(wp) :: temp
 
        real(wp),parameter :: twopi  = 8.0_wp * atan (1.0_wp)
        real(wp),parameter :: bdl    = - 1.0_wp
        real(wp),parameter :: bdu    = 1.0_wp
        integer,parameter  :: iprint = 2
        integer,parameter  :: maxfun = 500000
        real(wp),parameter :: rhobeg = 1.0e-1_wp
        real(wp),parameter :: rhoend = 1.0e-6_wp
        
        m = 5
        do
            n = 2 * m
            do i = 1, n
                xl (i) = bdl
                xu (i) = bdu
            end do
            do jcase = 1, 2
                npt = n + 6
                if (jcase == 2) npt = 2 * n + 1
                print 30, m, n, npt
30              format (/ / 5 x, '2D output with M =', i4, ',  N =', i4, '  and  NPT =', i4)
                do j = 1, m
                    temp = real (j, wp) * twopi / real (m, wp)
                    x (2*j-1) = cos (temp)
                    x (2*j) = sin (temp)
                end do
                call bobyqa (n, npt, x, xl, xu, rhobeg, rhoend, iprint, maxfun, calfun)
            end do
            m = m + m
            if (m > 10) exit
        end do
 
    contains
 
        subroutine calfun (n, x, f)
        
            implicit none
            
            integer,intent(in)               :: n
            real(wp),dimension(:),intent(in) :: x
            real(wp),intent(out)             :: f
            
            integer :: i,j
            real(wp) :: temp
            
            f = 0.0_wp
            do i = 4, n, 2
                do j = 2, i - 2, 2
                    temp = (x(i-1)-x(j-1)) ** 2 + (x(i)-x(j)) ** 2
                    temp = max (temp, 1.0e-6_wp)
                    f = f + 1.0_wp / sqrt (temp)
                end do
            end do
            
        end subroutine calfun
 
    end subroutine bobyqa_test



end program test_bobyqa