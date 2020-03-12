function Nval = nLen(xa,ya)

    addpath matlab_script/

    input_p = 'fringe_m30.f00008';
    nelx = 165;
    nely = 40;
    [data_p,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek(input_p);
    [xx,yy,vx,vy,p,t] = reshapenek(data_p,nelx,nely);
    [ny, nx] = size(xx);
    spr = zeros(nx,1);
    dx = xx(1,2:end) - xx(1,1:end-1);
    dy = yy(1,2:end) - yy(1,1:end-1);
    spr(2:end) = cumsum(sqrt(dx.^2+dy.^2));
    N = zeros(ny,nx);
    for i=1:nx
        dx = xx(:,i) - xx(1,i);
        dy = yy(:,i) - yy(1,i);
        N(:,i) = sqrt(dx.^2+dy.^2);
    end
    nxa = length(xa);
    Nval = zeros(nxa,1);
    
    x = xx(1,:);
    y = yy(1,:);
    cond = y>=0;
    x = x(cond);
    y = y(cond);
    spr = spr(cond);
    xprfun = csapi(x,spr);
    xaa = xa(ya>=0);
    sa = fnval(xprfun,xaa);
    %method = 'spline';
    %sa = interp1(x,spr,xa, method)
    Nff = N(end,:);
    Nff = Nff(cond);
    nfn = csapi(spr, Nff);
    Nval(nxa-length(xaa)+1:nxa) = fnval(nfn,sa);
    %Nval = interp1(spr,N(end,:),sa,method);
    
    x = xx(1,:);
    y = yy(1,:);
    cond = y<0;
    x = x(cond);
    y = y(cond);
    spr = spr(cond);
    xprfun = csapi(x,spr);
    xaa = xa(ya<0);
    sa = fnval(xprfun,xaa);
    %method = 'spline';
    %sa = interp1(x,spr,xa, method)
    Nff = N(end,:);
    Nff = Nff(cond);
    nfn = csapi(spr, Nff);
    Nval(1:length(xaa)) = fnval(nfn,sa);
    %Nval = interp1(spr,N(end,:),sa,method);
   
    % figure()
    % hold on
    % plot(xx(1,:),N(end,:),'k.')
    % plot(S(1,:),N(end,:),'b.')
end