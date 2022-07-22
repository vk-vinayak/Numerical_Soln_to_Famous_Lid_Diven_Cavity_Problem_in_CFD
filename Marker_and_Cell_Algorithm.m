%This code solves famous lid driven cavity problem in CFD with Marker & Cell Algorithm. 
%It uses finite difference method for solving
%transport and poisson's equation. Size of the cavity is 0.1 x 0.1m and
%velocity of the lid is 0.1 m/s with characteristic Re=100.  The code 
%outputs vorticity, stream function, u-velocity & y-velocity contour with 10 x 10 grid size.
function a11_new_1
%Basic Geometry definition
%No. of nodes
n=10;
nn=n+2;
l=0.1;
U=0.1; 
nu=10^-4;
dt=0.00005;
dx=l/n;
dy=l/n;
kk=0;
x=zeros(nn,nn);
y=zeros(nn,nn);
x(1,1)=dx/2;
y(1,1)=dy/2;
for i=2:nn
    x(i,:)=x(i-1,:)+dx;
end
for j=2:nn
    y(:,j)=y(:,j-1)+dy;
end

%initial conditions
u=zeros(nn,nn);
v=zeros(nn,nn);
p=ones(nn,nn);
R=zeros(nn,nn);
D=zeros(nn,nn);
Q=zeros(nn,nn);
uc=zeros(nn,nn);
vc=zeros(nn,nn);
psi=zeros(nn,nn);

%boundary conditions
    u(1,:)=0;
    v(1,:)=-v(2,:);
    p(1,:)=p(2,:)-2*nu.*u(2,:)/dx;
    u(nn-1,:)=0;
    u(nn,:)=u(nn-2,:);
    v(nn,:)=-v(nn-1,:);
    p(nn,:)=p(nn-1,:)+2*nu.*u(nn-1,:)/dx; 
    u(:,1)=-u(:,2);
    v(:,1)=0;
    p(:,1)=p(:,2)-2*nu.*v(:,2)/dy;
    u(:,nn)=2*U-u(:,nn-1);
    v(:,nn)=v(:,nn-2);
    v(:,nn-1)=0;
    p(:,nn)=p(:,nn-1)+2*nu*v(:,nn-1)/dy;
u0=zeros(nn);
v0=zeros(nn);

err=1;
time=0;
%time till which we want to do analysis
final_time=0.05;
%Loop for every time step
for time=0:dt:final_time
    %Loop to converge continuity
    for iter=1:25
        %Defining uc at requirednodes
        for i=2:nn
            for j=2:nn
                uc(i,j) = (u(i-1,j)+u(i,j))/2;
                vc(i,j) = (v(i,j-1)+v(i,j))/2;
            end
        end
        
        for j=1:nn
            u0(j) = u(2,j);
            uc(1,j) = (u0(j)+u(1,j))/2;
        end
        for i=1:nn
            v0(i) = v(i,2);
            vc(i,1) = (v0(i)+v(i,2))/2;
        end
        
        u_old=u;
        v_old=v;
        uc_old=uc;
        vc_old=vc;
        p_old=p;
        err_p=1;
        %Loop for pressure calculations
        for iter_p=1:18
            n=0;
            
            %Boundary condition for D
            for j=1:nn
                D(1,:)=0;
                D(nn,:)=0;
            end
            for i=1:nn
                D(:,1)=0;
                D(:,nn)=0;
            end
            %Calculation of Q at every node
            for i=2:nn-1
                for j=2:nn-1
                    n=n+1;
                    D(i,j)=(u_old(i,j)-u_old(i-1,j))/dx+(v_old(i,j)-v_old(i,j-1))/dy;
                    Q_term1 = (uc_old(i+1,j)^2-2*uc_old(i,j)^2+uc_old(i-1,j)^2)/(dx^2);
                    Q_term2 = (vc_old(i,j+1)^2-2*vc_old(i,j)^2+vc_old(i,j-1)^2)/(dy^2);
                    Q_term3 = 2*((u_old(i,j)+u_old(i,j+1))/2*(v_old(i,j)+v_old(i+1,j))/2 ...
                        +(u_old(i-1,j)+u_old(i-1,j-1))/2*(v_old(i,j-1)+v_old(i-1,j-1))/2 ...
                        +(u_old(i,j)+u_old(i,j-1))/2*(v_old(i,j-1)+v_old(i+1,j-1))/2 ...
                        +(u_old(i-1,j)+u_old(i-1,j+1))/2*(v_old(i,j)+v_old(i-1,j))/2)/(dx*dy);
                    
                    Q(i,j) = Q_term1+Q_term2+Q_term3;
                    
                end
            end
            %Calculation of R, Pressure at every node
            for i=2:nn-1
                for j=2:nn-1
                    R(i,j) = -Q(i,j) - (-D(i,j))/dt ...
                        +nu*((D(i+1,j)^2-2*D(i,j)^2-D(i-1,j)^2)/(dx^2) ...
                        +(D(i,j+1)^2-2*D(i,j)^2-D(i,j-1)^2)/(dy^2));
                    %Poisson equation for pressure
                    p(i,j) = (-R(i,j)+(p_old(i-1,j)+p_old(i+1,j))/(dx^2)+(p_old(i,j-1)+p_old(i,j+1))/(dy^2))/(2/dx^2+2/dy^2);
                end
            end
            %Extra elements(outside) for pressure
            %Top and bottom wall
            for j=1:nn
                p(1,j)=p(2,j)-2*nu*u(2,j)/dx;
                p(nn,j)=p(nn-1,j)+2*nu*u(nn-1,j)/dx;
            end
            %side walls
            for i=1:nn
                p(i,1)=p(i,2)-2*nu*v(i,2)/dy;
                p(i,nn)=p(i,nn-1)+2*nu*v(i,nn-1)/dy;
            end
            err_p = max(max(abs(p-p_old)));
        end
        
        %calculation of u
        
        for i=2:nn-2
            for j=2:nn-1
                
                pr_term = (p_old(i,j)-p_old(i+1,j))/dx;
                diff = nu*((u_old(i+1,j)-2*u_old(i,j)+u_old(i-1,j))/(dx*dx) + (u_old(i,j+1)-2*u_old(i,j)+u_old(i,j-1))/(dy*dy));
                conv_u = (uc_old(i,j)^2-uc_old(i+1,j)^2)/dx;
                conv_v = ((u_old(i,j)+u_old(i,j-1))/2*(v_old(i,j-1)+v_old(i+1,j-1))/2 ...
                    -(u_old(i,j)+u_old(i,j+1))/2*(v_old(i,j)+v_old(i+1,j))/2)/dx;
                u(i,j) = u_old(i,j) + dt*(pr_term+diff+conv_u+conv_v);
            end
        end
        %calculation of v
        for i=2:nn-1
            for j=2:nn-2 
                pr_term = (p_old(i,j)-p_old(i,j+1))/dy;
                diff = nu*((v_old(i+1,j)-2*v_old(i,j)+v_old(i-1,j))/(dx*dx) + (v_old(i,j+1)-2*v_old(i,j)+v_old(i,j-1))/(dy*dy));
                conv_v = (vc_old(i,j)^2-vc_old(i,j+1)^2)/dy;
                conv_u = ((u_old(i-1,j)+u_old(i-1,j+1))/2*(v_old(i,j)+v_old(i-1,j))/2 ...
                    -(u_old(i,j)+u_old(i,j+1))/2*(v_old(i,j)+v_old(i+1,j))/2)/dx;
                v(i,j) = v_old(i,j) + dt*(pr_term+diff+conv_u+conv_v);
                
            end
        end
        
        %updating boundary conditions
        
        %left face
        for j=1:nn
            u(1,:)=0;
            v(1,:)=-v(2,:);
            u(nn-1,:)=0;
            u(nn,:)=u(nn-2,:);
            v(nn,:)=-v(nn-1,:);
        end
        for i=1:nn
            u(:,1)=-u(:,2);
            v(:,1)=0;
            u(:,nn)=2*U-u(:,nn-1);
            v(:,nn)=v(:,nn-2);
            v(:,nn-1)=0;
        end
        for i=2:nn-1
            for j=2:nn-1
                D(i,j) = (u(i,j)-u(i-1,j))/dx+(v(i,j)-v(i,j-1))/dy;
            end
        end
        err = max(max(abs(D)));
        kk=kk+1;
        g(kk)=err;
    end
    
end
%All calculation at every time step and every node ended for p,u,v
for i=1:nn
    for j=1:nn
        uu(i,j) = u(i,j);
        vv(i,j) = v(i,j);
    end
end
%psi calculations
psi(1,:)=-psi(2,:);
psi(nn,:)=-psi(nn-1,:);
psi(:,1)=-psi(:,2);
psi(:,nn)=-psi(:,nn-1);
for i=2:nn-1
    for j=2:nn-2
        psi(i,j+1)=u(i,j)*dt+psi(j);
    end
end
%% plotting Graphs

figure
contourf(x,y,psi,20)
title('Streamline contours');
xlabel('x');
ylabel('y');

figure
contourf(x,y,u,20)
title('U-velocity contours');
xlabel('x');
ylabel('y');

figure
contourf(x,y,v,20)
title('V-velocity contours');
xlabel('x');
ylabel('y');

end