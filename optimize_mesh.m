function optimize_mesh(m)
  %clear all
  %close all
  
  
  % Some information
  to_gll = 0.17;
  Re = 5.3333e5;
  dstar = 4.2e-4; % corresponding to Res=30300
  obj1 = 2*pi*dstar/3
  
  dis_th = 1.536e-3; % at chord=0.4
  bl_th = 4.048e-3; %at chord=0.4
  dis_th2 = .372e-3; % at chord=0.04
  bl_th2 = 1.059e-3; % at chord=0.04
  
  
  % Orignal mesh nelx=165; rexp=0.55
  % lin =0.3; doutup = 0.1; doutlw=0.08
  % xinup = xinlw = -.05 ;
  
  % Parameters to tune
  nelx = [200, 165, 165, 200, 200, 200];
  nely = [40, 40, 40, 34, 34, 34];
  rexp = [0.55, 0.55, 0.55, 0.55, 0.35, 0.35];
  xprcutup = [.5, 0.4, 0.4, 0.4, 0.4, 0.4];
  xprcutlw = [.04, .04, .04, .04, .04, .04];
  lin = [.3, .3, .3, .06, 0.06, .06];
  xinup = [-.05, -.05, -.05, -.01, -.01, -.01];
  xinlw = [-.05, -.05, -.05, -.01, -.01, -.01];
  doutup = [.1, .1 , .1, 15*dis_th, 20*dis_th, 30*dis_th];
  doutlw = [.08, .08, .08, 10*dis_th, 15*dis_th, 20*dis_th ];
  
  span = [0.03, 0.03, 0.03, 0.02, 0.02, 0.02];
  nelz = [18, 18, 18, 12, 12, 12];
  
  
  
  % Get the reference dUTdn from base flow sim
  [dUTdn, xr] = dUdn();
  [val , ind] = min(xr);
  Re_tauref = sqrt(dUTdn(ind:end)*Re);
  xr = xr(ind:end);
  dz = span(m)/nelz(m)*to_gll;
 
  % get and unroll data for new mesh
  data = mesh_values(nelx(m),nely(m),rexp(m), xprcutup(m), xprcutlw(m),lin(m), xinup(m), xinlw(m), doutup(m), doutlw(m));
  
  xp = data.xpr;
  yp = data.ypr;
  x2 = data.x2;
  y2 = data.y2;
  sp = data.spr;
  xBC = data.xBC;
  yBC = data.yBC;
  
  
  N_total = nelx(m)*nely(m)*nelz(m)
  lambda_start = dstar/0.23*2*pi
  
  
  % position at the mid point of elemnts
  xe = (xp(2:end) + xp(1:end-1))/2;
  ye = (yp(2:end) + yp(1:end-1))/2;
  xeBC = (xBC(2:end) + xBC(1:end-1))/2;
  yeBC = (yBC(2:end) + yBC(1:end-1))/2;
  
  dn = sqrt((x2-xp).^2 + (y2-yp).^2);
  ds = sp(2:end)-sp(1:end-1);
  
  % Compute the length of the elements at the boundary
  dx = xBC(2:end)-xBC(1:end-1);
  dy = yBC(2:end)-yBC(1:end-1);
  sBC = zeros(nelx(m),1);
  sBC = sqrt(dx.^2 + dy.^2);
  
  % find the right index to get the minimum
  [val, ind] = min(abs(xBC.*sign(yBC)));
  SBC_max = max(sBC(1:ind))
  
  Re_tau = interp1(xr,Re_tauref, xe);
  ds_plus = Re_tau.*ds*to_gll;
  
  Re_tau = interp1(xr,Re_tauref, xp);
  dn_plus = Re_tau.*dn;
  
  figure(9900)
  hold on
  plot(xp,dz*Re_tau)
  xlabel('Chord')
  ylabel("$\Delta z^+_{GLL}$", 'Interpreter','latex')
  title('$\Delta z^+_{GLL}$ along the chord','Interpreter','latex')
  
  figure(1000)
  hold on
  plot(xe,ds_plus)
  xlabel('Chord')
  ylabel("$\Delta s^+_{GLL}$", 'Interpreter','latex')
  title('$\Delta s^+_{GLL}$ along the chord','Interpreter','latex')
  
  figure(1001)
  hold on
  plot(xp,dn_plus)
  xlabel('Chord')
  ylabel("$\Delta n^+_{el}$", 'Interpreter','latex')
  title('$\Delta n^+_{el}$ along the chord','Interpreter','latex')
  
  figure(1002)
  cond = xeBC<xinup(m);
  hold on
  plot(yeBC(cond), obj1./(sBC(cond)*to_gll))
  %plot(yeBC(cond), obj1*ones(length(yeBC(cond)),1))
  xlabel('y')
  ylabel("$\lambda_{end}/\Delta s_{gll}$", 'Interpreter','latex')
  title('$\lambda_{end}/\Delta s_{gll}$ at FST BC', 'Interpreter','latex')
  %legend('Mesh values','Required')
  hold off
  
  figure(1003)
  ind = nelx(m);
  hold on
  scatter(xeBC(1:ind), yeBC(1:ind), [],sBC(1:ind), 'filled')
  scatter(xe(1:ind),ye(1:ind),[],ds(1:ind),'filled')
  colorbar()
  axis('equal')
  xlabel('x')
  ylabel('y')
  title('$\Delta S_{el}$','Interpreter','latex')
end
