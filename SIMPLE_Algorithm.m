%This code solves famous lid driven cavity problem in CFD with SIMPLE
%Algorithm.It outputs pressure, u-velocity, v-velocity contours. It
%compares u-velocity and v-velocity plot along the central vertical line of the cavity with Ghia's paper
function a12_new
rho=1;
nu=.0001;
un=0.1; 
mu=rho*nu;
l=0.1;
w=0.1;
Re=un*l/nu;
Dh=0.001;
Dl=0.001;
dl=Dh;
dh=Dl;
sx=(l/dl)+2;
sy=(w/dh)+2;

  
%% Initializing pressure and velocity

%initial pressure 
p=zeros(sx,sy);
            
% u velocity boundary conditions and initial condition.
u=zeros(sy,sx-1);
u(1:sy,1)=0;              %left wall boundary
u(1:sy,sx-1)=0;         %right wall boundary condition
u(1,2:sx-1)=0;            %bottom boundary
u(sy,1:sx-1)=un;         %top boundary
uold=u;
ustar=u;
            
% v velocity boundary conditions and initial condition.
v=zeros(sy-1,sx);
vstar=v;
vold=v;
            
%time step
Dt_a=Dh/un;
Dt_d=2*nu*(1/Dh^2+1/Dl^2);
Dt_d=1/Dt_d;
dt=min(Dt_d,Dt_a);
dt=0.5*dt;
itr=0;

%% SIMPLE method

% main while
  epsi=1;
 
  while epsi>1e-3
 
   
%% U control volume fluxes
  
 %fluxes in x drection
    for i=1:sx-2
       for j=2:sy-1
           mxu(j,i)=rho*(u(j,i)+u(j,i+1))/2;
           dxu(j,i)=mu*(u(j,i+1)-u(j,i))/Dh ;
           pxu=mxu(j,i)/dxu(j,i);
           axu(j,i)=max(mxu(j,i),0)*(dxu(j,i)*(max(0,(1-0.1*abs(pxu))^5))+max(mxu(j,i),0)-(dxu(j,i))^5)*u(j,i)+min(mxu(j,i),0)*(dxu(j,i)*(max(0,(1-0.1*abs(pxu))^5))+min(mxu(j,i),0)-(dxu(j,i))^5)*u(j,i+1);
           if Dh==0.005
           axu(j,i)=max(mxu(j,i),0)*u(j,i)+min(mxu(j,i),0)*u(j,i+1);
           end
       end
    end
    
   %fluxes in y drection   
      for i=2:sx-2
       for j=1:sy-1
           myu(j,i)=rho*(v(j,i)+v(j,i+1))/2;
           ayu(j,i)=max(myu(j,i),0)*u(j,i)+min(myu(j,i),0)*u(j+1,i);
           dyu(j,i)=mu*(u(j+1,i)-u(j,i))/Dl ;
           pyu=myu(j,i)/dyu(j,i);
           ayu(j,i)=max(myu(j,i),0)*(dyu(j,i)*(max(0,(1-0.1*abs(pyu))^5))+max(myu(j,i),0)-(dyu(j,i))^5)*u(j,i)+min(myu(j,i),0)*(dyu(j,i)*(max(0,(1-0.1*abs(pyu))^5))+min(myu(j,i),0)-(dyu(j,i))^5)*u(j,i+1);
           if Dh==0.005
           ayu(j,i)=max(myu(j,i),0)*u(j,i)+min(myu(j,i),0)*u(j+1,i);
           end
       end
      end

   
%% U control volume Predictor Step

   for  i=2:sx-2
       for j=2:sy-1
           Add(j,i)=(axu(j,i)-axu(j,i-1))*Dl+(ayu(j,i)-ayu(j-1,i))*dl;
           Diff(j,i)=(dxu(j,i)-dxu(j,i-1))*Dl+(dyu(j,i)-dyu(j-1,i))*dl;
           Sourc(j,i)=(p(j,i)-p(j,i+1))*Dl;
           ustar(j,i)=uold(j,i)+(dt/(rho*Dl*dl))*(Diff(j,i)-Add(j,i)+Sourc(j,i));
       end
   end 
 
%% V control volume fluxes

  %fluxes in x direction.
  for i=1:sx-1
      for j=2:sy-2
          mxv(j,i)=rho*(u(j,i)+u(j+1,i))/2;
          axv(j,i)=max(mxv(j,i),0)*v(j,i)+min(mxv(j,i),0)*v(j,i+1);
          dxv(j,i)=mu*(v(j,i+1)-v(j,i))/dl;
          pxv=mxv(j,i)/dxv(j,i);
          axv(j,i)=min(mxv(j,i),0)*(dxv(j,i)*(max(0,(1-0.1*abs(pxv))^5))+max(mxv(j,i),0)-(dxv(j,i))^5)*v(j,i+1)+max(mxv(j,i),0)*(dxv(j,i)*(max(0,(1-0.1*abs(pxv))^5))+min(mxv(j,i),0)-(dxv(j,i))^5)*v(j,i);
          if Dh==0.005
          axv(j,i)=max(mxv(j,i),0)*v(j,i)+min(mxv(j,i),0)*v(j,i+1);
          end
      end
  end
  
  %fluxes in y direction.
  for i=2:sx-1
      for j=1:sy-2
          myv(j,i)=rho*(v(j+1,i)+v(j,i))/2;
          ayv(j,i)=max(myv(j,i),0)*v(j,i)+min(myv(j,i),0)*v(j+1,i);
          dyv(j,i)=mu*(v(j+1,i)-v(j,i))/Dl;
          pyv=myv(j,i)/dyv(j,i);
          ayv(j,i)=min(myv(j,i),0)*(dyv(j,i)*(max(0,(1-0.1*abs(pyv))^5))+max(myv(j,i),0)-(dyv(j,i))^5)*v(j,i+1)+max(myv(j,i),0)*(dyv(j,i)*(max(0,(1-0.1*abs(pyv))^5))+min(myv(j,i),0)-dyv(j,i))*v(j,i);
          if Dh==0.005
          ayv(j,i)=max(myv(j,i),0)*v(j,i)+min(myv(j,i),0)*v(j+1,i);
          end
      end
  end
  
%% V control volume Predictor Step

   for  i=2:sx-1
       for j=2:sy-2
           Add(j,i)=(axv(j,i)-axv(j,i-1))*dh+(ayv(j,i)-ayv(j-1,i))*Dh;
           Diff(j,i)=(dxv(j,i)-dxv(j,i-1))*dh+(dyv(j,i)-dyv(j-1,i))*Dh;
           Sourc(j,i)=(p(j,i)-p(j+1,i))*Dh;
           vstar(j,i)=vold(j,i)+(dt/(rho*Dh*dh))*(Diff(j,i)-Add(j,i)+Sourc(j,i));
       end
   end 
   
 
%% SIMPLE method loop

 pdash(1:sy,1:sx)=0;
 error=1;
 
   while error>1e-3
    
        pdash(1:sy,1)=pdash(1:sy,2);         %left boundary contion
        pdash(1:sy,sx)=pdash(1:sy,sx-1); %right boundary condition
        pdash(1,1:sx)=pdash(2,1:sx);         %bottom boundary
        pdash(sy,1:sx)=pdash(sy-1,1:sx); %top boundary condition
        
       for i=2:sx-1
           for j=2:sy-1
                 b(j,i)=(ustar(j,i)-ustar(j,i-1))*rho*Dl+(vstar(j,i)-vstar(j-1,i))*rho*Dh;
           end
       end
       
 %% point by point gauss seidel method
       
        aE=rho*dt*Dl/Dh;
        aW=aE;
        aN=rho*dt*Dh/Dl;
        aS=aN;
        aP=aW+aE+aS+aN;
        
       for i=2:sx-1
           for j=2:sy-1
               pdash(j,i)=aE*pdash(j,i+1)+aW*pdash(j,i-1)+aN*pdash(j+1,i)+aS*pdash(j-1,i)-b(j,i);
               pdash(j,i)=pdash(j,i)/aP;
               p(j,i)=p(j,i)+pdash(j,i);
           end
       end
       for i=2:sx-1
           for j=2:sy-1
               if i<(sx-1)
                   ustar(j,i)=ustar(j,i)+(dt/(rho*dl*Dl))*(pdash(j,i)-pdash(j,i+1))*Dl;
               end
               if j<(sy-1)
                   vstar(j,i)=vstar(j,i)+(dt/(rho*Dh*dh))*(pdash(j,i)-pdash(j+1,i))*Dh;
               end
           end
       end
       
    error=max(max(abs(b)));
   end
    
    espiu=max(max(abs(ustar-uold)))/dt;
    epsiv=max(max(abs(vstar-vold)))/dt;
    
    epsi=max(espiu,epsiv);
    
    vold=vstar;
    v=vstar;
    uold=ustar;
    u=ustar;
    
    itr=itr+1;
%     plot(itr,epsi,'ro');
%     hold on;
    
  end
  
  
  %% Output
  
  
  %plotting of Grid points
  x(1)=0;
  x(sx-1)=l;
  
  for i=2:sx-2
  x(i)=x(i-1)+Dh;
  end
  
  y=x;
  
  %%plotting control surfaces
  xc(1)=x(1);xc(sx)=x(sx-1);
  for i=2:sx-1
      xc(i)=(x(i)+x(i-1))/2;
  end
  
  yc=xc;
  
  [xx yy]=meshgrid(xc,yc);
   
  % vector plot
  U=zeros(sy,sx);
  V=zeros(sy,sx);
  U(sy,1:sx)=un;
     for j=2:sy-1
         for i=2:sx-1
             U(j,i)=(u(j,i)+u(j,i-1))/2;
         end
     end
     for j=2:sy-1
         for i=2:sx-1
             V(j,i)=(v(j,i)+v(j-1,i))/2;
         end
     end
     

    


    figure
    contour(xc,yc,p,15);
    xlim([0 l])
    ylim([0 w])
    xlabel('X cordinates');
    ylabel('Y cordinates');
    title('pressure contours')
    
    figure
    colormap jet
    contour(xc,yc,U,15);
    xlim([0 l])
    ylim([0 w])
    xlabel('X cordinates');
    ylabel('Y cordinates');
    title('U-velocity contour');
    
    
    %% calculation for midsection velocities
    
    c=mod(sx,2);
    if c==0
        n=sx/2;
        avgu=(U(:,n)+U(:,n+1))/2;
        avgv=(V(n,:)+V(n+1,:))/2;
    else
        n=(sx+1)/2;        
        avgu=U(:,n);        
        avgv=V(n,:);
    end
    
    %% Comparison with ghia et al
    
    % Note: Ghia data at Re 100 is calculated. 
    
    figure
    plot(avgu/un,yc/(w),'b-o','MarkerFaceColor','r','MarkerSize',3)
    ylim([0 1])
    xlabel('U');
    ylabel('Y cordinates');
    title('U-velocity plot');
    hold on
    u_ghia=[1.0000 0.65928 0.78871 0.73722 0.68717 0.23151 0.00332 -0.136641 -0.20581 -0.2109 -0.15662 -0.1015 -0.063434 -0.04775 -0.04192 -0.03717 0.00000];
    y_ghia=[1 0.9766 0.9688 0.9609 0.9531 0.8516 0.7344 0.6172 0.5 0.4531 0.2813 0.1719 0.1016 0.0703 0.0625 0.0547 0.0000];
   
        plot(u_ghia,y_ghia,'r-*','MarkerFaceColor','r');
        legend('Numerical data','ghia data');
   
    

    figure
    plot(xc/(l),avgv/un,'b-o','MarkerFaceColor','r','MarkerSize',3);
    xlim([0 1])
    xlabel('X cordinates');
    ylabel('V-velocity');
    title('V-velocity plot');
    hold on;
    x_ghia=[1 0.9688 0.9609 0.9531 0.9453 0.9063 0.8594 0.8047 0.5 0.2344 0.2266 0.1563 0.0938 0.0781 0.0703 0.0625 0];
    v_ghia= [0 -0.05906 -0.07391 -0.08864 -0.10313 -0.16914 -0.22445 -0.24533 0.05454 0.17527 0.17507 0.16077 0.12317 0.1089 0.100091 0.09233 0];
    
  
        plot(x_ghia,v_ghia,'r-d','MarkerFaceColor','r');
        legend('Numerical data','ghia data');
   
  
   
    

hold off
end


