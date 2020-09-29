          implicit real*8 (a-h,o-z)
          
          dimension x(0:101),xi(0:102),u(0:101), 
     &              q1(-1:102),q1new(1:100),
     &              q2(-1:102),q2new(1:100),
     &              q3(-1:102),q3new(1:100), 
     &              f1(101),f2(101),f3(101),
     &              Pr(0:101)
          
          x0=0  
          x1=1
          nx=100
          dx=0.01
          gama=3.   
          fac=0.02
c-----------------------------------------------------------      
c---------SETTING  THE  GRID--------------------------------
    
          do 100 i=0,nx+1 
          x(i)=(i-1.)*dx 
100       continue

          do 200 i=1,nx+1
          xi(i)=0.5*(x(i)+x(i-1))
200       continue

          xi(0)=x(0)-0.5*(x(1)-x(0))
          xi(nx+2)=x(nx+1)+0.5*(x(nx+1)-x(nx)) 
c-----------------------------------------------------------
c--------ASSIGNING INITIAL CONDITION------------------------
          
          dg=0.1*(x1-x0)
          xmid=0.5*(x0+x1)

          do 300 i=1,nx
          
          q1(i)=dexp(-((x(i)-xmid)**2)/(2.*dg**2))
          q2(i)=0.5*q1(i)
          q3(i)=0.5*q1(i)*fac+0.5*q2(i)**2/q1(i)  

c          if(x(i).le.0.5)then
c          q1(i)=1.0d0
c          q2(i)=q1(i)*(-1.0)
c          q3(i)=0.5*q1(i)*(-1.0)**2+(0.4)/(gama-1.)
c          else
c          q1(i)=1.0
c          q2(i)=q1(i)*1.0
c          q3(i)=0.5*q1(i)*(1.0)**2+(0.4)/(gama-1.)
c          endif
c          write(*,*)x(i),q1(i)

300       continue   
c-----------------------------------------------------------
c---------ADVECTION-----------------------------------------
          
          dt=0.003d0
          tend=0.65
    
          do 400 t=0.0d0,tend,dt

c----------------------------------------------------------

          q1(-1)=q1(nx-1)     ! Periodic Boundary Conditions
          q1(0)=q1(nx)        ! Periodic Boundary Conditions
          q1(nx+1)=q1(1)
          q1(nx+2)=q1(2)
          
          q2(-1)=q2(nx-1)     ! Periodic Boundary Conditions
          q2(0)=q2(nx)        ! Periodic Boundary Conditions
          q2(nx+1)=q2(1)
          q2(nx+2)=q2(2)
          
          q3(-1)=q3(nx-1)     ! Periodic Boundary Conditions
          q3(0)=q3(nx)        ! Periodic Boundary Conditions
          q3(nx+1)=q3(1)
          q3(nx+2)=q3(2)
c----------------------------------------------------------
          do 444 i=0,nx+1

c-----------SETTING  VELOCITY & PRESSURE-------------------
          u(i)=fac*t*(x(i)-xmid)/(fac*(t**2)+dg**2)     
c          u(i)=(q2(i)/q1(i)+q2(i-1)/q1(i-1))/2. 
          Pr(i)=(gama-1.)*(q3(i)-0.5*q2(i)**2/q1(i))
444       continue
c---------------------------------------------------------         
         do 555 i=1,nx+1    

         if(u(i).ge.0)then
         f1(i)=q1(i-1)*u(i)+(u(i)/4.)*(1.-u(i)*dt/dx)*(q1(i)-q1(i-2))
         f2(i)=q2(i-1)*u(i)+(u(i)/4.)*(1.-u(i)*dt/dx)*(q2(i)-q2(i-2))
         f3(i)=q3(i-1)*u(i)+(u(i)/4.)*(1.-u(i)*dt/dx)*(q3(i)-q3(i-2)) 
         else
         f1(i)=q1(i)*u(i)-(u(i)/4.)*(1.+u(i)*dt/dx)*(q1(i+1)-q1(i-1))
         f2(i)=q2(i)*u(i)-(u(i)/4.)*(1.+u(i)*dt/dx)*(q2(i+1)-q2(i-1))
         f3(i)=q3(i)*u(i)-(u(i)/4.)*(1.+u(i)*dt/dx)*(q3(i+1)-q3(i-1)) 
         endif

555       continue   

          do 600 i=1,nx       
            
          q1new(i)=q1(i)-dt*(f1(i+1)-f1(i))/(xi(i+1)-xi(i)) 
          q2new(i)=q2(i)-dt*(f2(i+1)-f2(i))/(xi(i+1)-xi(i))
     &            -(dt/(2.*dx))*(Pr(i+1)-Pr(i-1))
          q3new(i)=q3(i)-(dt/(xi(i+1)-xi(i)))*(f3(i+1)-f3(i))
     &            -(dt/(2.*dx))*(Pr(i+1)*u(i+1)-Pr(i-1)*u(i-1))
          
          q1(i)=q1new(i)
          q2(i)=q2new(i)
          q3(i)=q3new(i)          

600       continue 
c          write(*,*)t,(q3(50)-0.5*q1(50)*(u(50)**2))/(0.5*q1(50))
400       continue      

          do i=1,nx
          write(*,*)x(i),u(i)    !q1(i),q3(i)    !(q3(i)-0.5*q1(i)*u(i)**2)/q1(i)
          enddo

          end
c--------------------------------------------------------------------    
c--------------------------------------------------------------------
