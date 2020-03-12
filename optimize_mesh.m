clear all
close all

% Orignal mesh nelx=165, rexp=0.55
nelx = 200;
rexp = 0.55;

to_gll = 0.17;
Re = 5.3333e5;
nely = 40;
obj1 = 0.0012;

% Get the reference dUTdn from old simulation
[dUTdn, xr] = dUdn();
[val , ind] = min(xr);
Re_tauref = sqrt(dUTdn(ind:end)*Re);
xr = xr(ind:end);


% get and unroll data for new mesh
data = mesh_values(nelx,nely,rexp);

xp = data.xpr;
yp = data.ypr;
sp = data.spr;
xBC = data.xBC;
yBC = data.yBC;

xe = (xp(2:end) + xp(1:end-1))/2;
ye = (yp(2:end) + yp(1:end-1))/2;
xeBC = (xBC(2:end) + xBC(1:end-1))/2;
yeBC = (yBC(2:end) + yBC(1:end-1))/2;

ds = sp(2:end)-sp(1:end-1);
Re_tau = interp1(xr,Re_tauref, xe);

% Compute the length of the elements at the boundary
dx = xBC(2:end)-xBC(1:end-1);
dy = yBC(2:end)-yBC(1:end-1);
sBC = zeros(nelx,1);
sBC = sqrt(dx.^2 + dy.^2);

% find the right index to get the minimum
[val, ind] = min(abs(xBC.*sign(yBC) + 0.08));
SBC_max = max(sBC(1:ind))

ds_plus = Re_tau.*ds*to_gll;


figure()
plot(xe,ds_plus)
xlabel('Chord')
ylabel("$\Delta s^+$", 'Interpreter','latex')
title('$\Delta s^+$ along the chord','Interpreter','latex')

figure()
cond = xeBC<-0.07;
hold on
plot(yeBC(cond), sBC(cond))
plot(yeBC(cond), obj1*ones(length(yeBC(cond)),1))
xlabel('y')
ylabel("$\Delta s$", 'Interpreter','latex')
title('$\Delta s$ at FST BC', 'Interpreter','latex')
legend('Mesh values','Needed')
hold off

figure()
ind = nelx;
hold on
scatter(xeBC(1:ind), yeBC(1:ind), [],sBC(1:ind), 'filled')
scatter(xe(1:ind),ye(1:ind),[],ds(1:ind),'filled')
colorbar()
axis('equal')
xlabel('x')
ylabel('y')
title('$\Delta S$','Interpreter','latex')