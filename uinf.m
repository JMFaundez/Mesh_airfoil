function [u99,y99,dth] = uinf(x)

    addpath matlab_script/
    input_p = 'fringe_m30.f00008';
    nelx = 165;
    nely = 40;
    [data_p,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek(input_p);
    [xx,yy,vx,vy,p,t] = reshapenek(data_p,nelx,nely);

    xa = xx(1,:);
    [x0, in] = min(abs(xa));
    in = 550;
    xa = xa(in:end);
    pa = p(1,in:end);
    [ny,nx] = size(xx);
    alpha = zeros(ny,nx);
    dx = xx(10,:) - xx(1,:);
    dy = yy(10,:) - yy(1,:);
    L = sqrt(dx.^2 + dy.^2);
    a = acos(dy./L);
    for i=1:nx
        alpha(:,i) = a(i);
    end
    U0 = 1;
    pinf = 0;
    pref = U0^2/2+pinf;

    Uv = sqrt((pref-pa)*2);
    cp = (pa-pinf)/(pref - pinf);
    
    ut = vx.*cos(alpha) + vy.*sin(alpha);
    %figure()
    %contourf(xx,yy,ut)
    %colorbar()
    [val, ind] = min(abs(xa -x));
    Uinf = Uv(ind);
    [val, ind] = min(abs(x-xx(1,:)));
    for j=1:ny
        if ut(j,ind)>=0.99*Uinf
            u99 = ut(j,ind);
            y99 = sqrt((yy(j,ind)-yy(1,ind))^2+(xx(j,ind)-xx(1,ind))^2);
            break
        end
    end
    
    integrand = 1 - ut(1:j,ind)/Uinf;
    dth = trapz(yy(1:j,ind),integrand);
    
    %figure()
    %hold on
    %plot(xx(1,:),p(1,:))
    %plot(xa,pa)


    %figure()
    %plot(xa,cp)

    %figure()
    %plot(xa, Uinf)
end

