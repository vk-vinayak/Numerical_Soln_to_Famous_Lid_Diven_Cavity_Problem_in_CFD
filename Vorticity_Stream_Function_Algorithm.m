%This code solves famous lid driven cavity problem in CFD with stream
%function-vorticity method. It uses finite difference method for solving
%transport and poisson's equation. Size of the cavity is 0.1 x 0.1m and
%velocity of the lid is 0.1 m/s with characteristic Re=100. The code asks
%for inputs like x & y grid (128) and time (20 sec). The code took around
%40 mins to run on the machine with 8gb RAM and intel-i7 processor. The
%code outputs drag force on the lid, compares our u-velocity along the vertical
%line through geometric center of the cavity with Ghia paper. The code also
%outputs vorticity, stream function, u-velocity & y-velocity contour.
function a10_2
%geometry and input constants
L=0.1;
H=0.1;
U=0.1;
rho=1;
gamma=10^-4;
x_grid=input('Enter x_grid');
y_grid=input('Enter y_grid');
T=input('Time till which you want analysis');
delta_x=L/x_grid;
delta_y=H/y_grid;
X=x_grid+1;
Y=y_grid+1;
%delta_t is calculated on paper such that for 128 by 128 grids and u=0.1
%stability condition is satisfied
delta_t=2*10^-4;
ww=1;
%Initial condition
w=zeros(Y,X);
si=zeros(Y,X);
u=zeros(Y,X);
u(end,:)=U;
v=zeros(Y,X);
%Boundary Condition
%si
si(1,:)=0;
si(end,:)=0;
si(:,1)=0;
si(:,end)=0;
t=0;
%stedy state is reached at high time;here20
T_max=(T+1-delta_t)./delta_t;
%calculation of w,si,u,v
for k=1:T_max
    t=t+delta_t;
for i=2:X-1
    for j=2:Y-1
%w boundary condition for w
w(1,:)=-si(2,:)*2/(delta_y.^2);
w(end,:)=-si(end-1,:)*2/(delta_y.^2)-2*U/delta_y;
w(:,1)=-si(:,2)*2/((delta_x)^2);
w(:,end)=-si(:,end-1)*2/((delta_x)^2);
%constants
        c1(j,i)=u(j,i)*delta_t/delta_x;
        c2(j,i)=v(j,i)*delta_t/delta_y;
        D1=gamma*delta_t/(delta_x^2);
        D2=gamma*delta_t/(delta_y^2);
        criteria1(j,i)=c1(j,i)^2/(2*D1);
        criteria2(j,i)=c2(j,i)^2/(2*D1);
        if(criteria1(j,i)<=1 && criteria2(j,i)<=1)
        %w calculation by ftcs(2-D Transport Equation)
      w(j,i)=w(j,i)+ww*(-1.*c1(j,i)*(w(j,i+1)-w(j,i-1))*0.5-c2(j,i)*(w(j+1,i)-w(j-1,i))*0.5+D1*(w(j,i+1)-2*w(j,i)+w(j,i-1))+D2*(w(j+1,i)-2*w(j,i)+w(j-1,i)));
      %si calculations (poissons equation)
      si(j,i)=((delta_y^2)*(si(j,i+1)+si(j,i-1))+(delta_x^2)*(si(j+1,i)+si(j-1,i))+w(j,i)*(delta_x^2)*(delta_y^2))/(2*delta_x^2+2*delta_y^2);
      %u,v calculation by central difference scheme
      u(j,i)=(si(j+1,i)-si(j-1,i))/(2*delta_y);
      v(j,i)=-(si(j,i+1)-si(j,i-1))/(2*delta_x);
        else
            fprintf('Criteria is not satisfied');
            criteria1(j,i)
            criteria2(j,i)
            
        end
        
    end
    
end

end
Drag_force=0;
for i=1:X
Drag_force=Drag_force+delta_x*gamma*rho*(u(end,i)-u(end-1,i))/(delta_y);
end
 Drag_force   
    

contourf(w)
title('Vorticity Contour');
figure

contourf(si)
title('Stream Function Contour');
figure

contourf(u)
title('u-velocity contour');
figure

contourf(v)
title('v-velocity contour');

end