          implicit real*8 (a-h,o-z)
          
          dimension x(100),a(100),b(100),c(100),r(100),v(100),
     &              AA(100),BB(100),rb(100),u(100),th_u(100)
          
          Pi=dacos(-1.d0)
          
          nx=100
          dx=0.01

          dx2=dx**2 

          do 100 i=1,nx 
          x(i)=(i*dx) 
          th_u(i)=(x(i))**4      !Theoretical value
100       continue 
         
          do 700 i=1,nx
 
          rb(i)=12.*(x(i))**2
          AA(i)=0.d0
          BB(i)=0.d0   
     
          a(i)=1.-(dx/2.)*AA(i)
          b(i)=dx2*BB(i)-2.
          c(i)=1.+(dx/2.)*AA(i)

700       continue
          
          do 20 i=2,nx-1
          r(i)=dx2*rb(i)
20        enddo   
        
          r(1)=dx2*rb(1)-a(1)*(0.d0)
          r(100)=dx2*rb(100)-c(100)*(1.d0)
          
          call solve_tridiag(nx,a,b,c,r,v)
        
          do 800 i=1,nx
          u(i)=v(i)   
          write(*,*)x(i),th_u(i),u(i)
800       continue
        
          end
c--------------------------------------------------------------------    
          
          subroutine tridag(n,a,b,c,r,u)
          implicit real*8 (a-h,o-z)

          PARAMETER (NMAX=500)
    
          dimension a(n),b(n),c(n),r(n),u(n),gam(NMAX)

          if(b(1).eq.0.)pause 'tridag: rewrite equations'
         
          bet=b(1)
          u(1)=r(1)/bet
          
          do 11 j=2,n
          gam(j-1)=c(j-1)/bet
          bet=b(j)-a(j)*gam(j-1)
          if(bet.eq.0.)pause 'tridag failed'
          u(j)=(r(j)-a(j)*u(j-1))/bet
11        enddo
 
          do 12 j=n-1,1,-1
          u(j)=u(j)-gam(j+1)*u(j+1)
12        enddo
          
          return
          END
c-------------------------------------------------------------------
          
           subroutine solve_tridiag(n,a,b,c,d,x)
           implicit none

!	 a - sub-diagonal (means it is the diagonal below the main diagonal)
!	 b - the main diagonal
!	 c - sup-diagonal (means it is the diagonal above the main diagonal)
!	 d - right part
!	 x - the answer
!	 n - number of equations

        integer,intent(in) :: n
        real(8),dimension(n),intent(in) :: a,b,c,d
        real(8),dimension(n),intent(out) :: x
        real(8),dimension(n) :: cp,dp
        real(8) :: m
        integer i

! initialize c-prime and d-prime
        cp(1) = c(1)/b(1)
        dp(1) = d(1)/b(1)
! solve for vectors c-prime and d-prime
         do i = 2,n
           m = b(i)-cp(i-1)*a(i)
           cp(i) = c(i)/m
           dp(i) = (d(i)-dp(i-1)*a(i))/m
         enddo
! initialize x
         x(n) = dp(n)
! solve for x from the vectors c-prime and d-prime
        do i = n-1, 1, -1
          x(i) = dp(i)-cp(i)*x(i+1)
        end do

        return
        end 
c---------------------------------------------------------------------
