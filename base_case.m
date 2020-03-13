function base_case()
  clear all
  addpath matlab_script/
  input_p = 'fringe_m30.f00008';
  Re = 5.33333e5;
  nelx = 165;
  nely = 40;
  [data_p,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek(input_p);
  [xx,yy,vx,vy,p,t] = reshapenek(data_p,nelx,nely);

  [ny,nx] = size(xx);
  xa = xx(1,:);
  pa = p(1,:);

  dx = xx(1,2:end) - xx(1,1:end-1);
  dy = yy(1,2:end) - yy(1,1:end-1);
  N = zeros(ny,nx);
  cos_a = zeros(ny,nx);
  sin_a = zeros(ny,nx);
  for i=1:nx
    dx = xx(:,i) - xx(1,i);
    dy = yy(:,i) - yy(1,i);
    N(:,i) = sqrt(dx.^2+dy.^2);
    cos_a(:,i) = dy(2)/N(2,i);
    sin_a(:,i) = -dx(2)/N(2,i);
  end
  U0 = 1;
  pinf = 0;
  pref = U0^2/2+pinf;
                                %pref = max(pa);
  Uv = sqrt((pref-pa)*2);
  cp = (pa-pinf)/(pref - pinf);

  ut = vx.*cos_a + vy.*sin_a;
  dUTdn = (ut(2,:) - ut(1,:))./N(2,:);

  dth = zeros(nx,1);
  y99 = zeros(nx,1);
  u99 = zeros(nx,1);

  for i=1:nx
    for j=1:ny
      if ut(j,i)>=0.99*abs(Uv(i))
        u99(i) = ut(j,i);
        y99(i) = N(j,i);
        integrand = 1 - abs(ut(1:j,i))/Uv(i);
        if j==1
          dth(i) = 0;
        else
          dth(i) = trapz(N(1:j,i),integrand);
        end
        break
      end
    end
  end

  [var, ind] = min(xx(1,:));
  dx = xx(1,ind+1:end) - xx(1,ind:end-1);
  dy = yy(1,ind+1:end) - yy(1,ind:end-1);
  ss = zeros(nx,1);
  ss(ind+1:end) = cumsum(sqrt(dx.^2 + dy.^2));
  Re_s = Re*ss;

  figure(2000)
  plot(Re_s(ind:end),p(1,ind:end))
  xlabel('$Re_s$','Interpreter','latex')
  ylabel('P')

  figure(2001)
  contourf(xx,yy,ut)
  xlabel('x')
  ylabel('y')
  colorbar()

  figure(2002)
  plot(Re_s(ind:end), Uv(ind:end))
  xlabel('$Re_s$','Interpreter','latex')
  ylabel("$U_{99}$",'Interpreter','latex')

  figure(2003)
  plot(Re_s(ind:end),dth(ind:end))
  ylabel('$\delta^*$','Interpreter','latex')
  xlabel('$Re_s$','Interpreter','latex')

  figure(2004)
  plot(Re_s(ind:end),y99(ind:end))
  xlabel('$Re_s$','Interpreter','latex')
  ylabel('$\delta$','Interpreter','latex')

  figure(2005)
  plot(Re_s(ind:end),dUTdn(ind:end))
  xlabel('$Re_s$','Interpreter','latex')
  ylabel('$\partial U_t \partial n$','Interpreter','latex')
end
