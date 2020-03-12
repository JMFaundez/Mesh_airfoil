function [dUTdn, x] = dUdn()

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
    
    x = xx(1,:);
    y = yy(1,:);
    
    UT = vx.*cos_a + vy.*sin_a;
    dUTdn = (UT(2,:) - UT(1,:))./N(2,:);
    
%     figure()
%     contourf(xx,yy,UT)
%     colorbar()
%     
%     figure()
%     plot(x,dUTdn(1,:))
%     
%     figure()
%     contourf(xx,yy,cos_a)
%     colorbar()
%     title('cosine')
%     
%     figure()
%     contourf(xx,yy,sin_a)
%     colorbar()
%     title('sin')
%     
%     figure()
%     hold on
%     plot(x,cos_a(1,:))
%     plot(y,sin_a(1,:))
%     legend('cosine','sine')
end