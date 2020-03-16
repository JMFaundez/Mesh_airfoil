function d = mesh_values(nelx,nely, rexp, xprcutup, xprcutlw, lin, xinup, xinlw, doutup, doutlw)

	% - boundary/initial condition FLUENT data
	dxbox = 2e-3; xbox = [-.2 .8]; nxbox = ceil(diff(xbox)/dxbox) + 1;
	dybox = 2e-3; ybox = [-.2 .2]; nybox = ceil(diff(ybox)/dybox) + 1;
	
	bcsource = 'FLUENT';
	bcdata = 'naca0008_fluent-data.mat';
	
	% -- profile
	%xprcutup = .4;.3;
	%xprcutlw = .04;
	
	% -- inflow
	%lin = .3;.25;.5;
	%xinup = -.05;-.00;-.05;-.00;
	%xinlw = -.05;-.00;-.05;-.00;
	
	% -- outflows
	%doutup = .10;
	%doutlw = .08;
	
	% 1/-- fringe
	%wdfrup = .050; rsfrup = .015;
	%wdfrlw = .030; rsfrlw = .010;

    %% Load original profile coordinates
    iaf.designation='0008';
    iaf.n = 1400;
    iaf.HalfCosineSpacing=1;
    iaf.wantFile=0;
    iaf.datFilePath='./'; % Current folder
    iaf.is_finiteTE=0;
    data = naca4gen(iaf); xpro =data.x'; ypro = flip(data.z)';
    clear data
    % tangential coordinate
    dx = xpro(2:end) - xpro(1:end-1);
    dy = ypro(2:end) - ypro(1:end-1);

    ds = sqrt(dx.^2+dy.^2);
    spro = cumsum([0 ds]); [~,ispro0] = min(abs(xpro)); spro = spro - spro(ispro0);

    % profile spline
    xprfun = csapi(spro,xpro);
    yprfun = csapi(spro,ypro);


    % curvature radius (spline differentiation)
    d1x = fnval(fnder(xprfun,1),spro); d2x = fnval(fnder(xprfun,2),spro);
    d1y = fnval(fnder(yprfun,1),spro); d2y = fnval(fnder(yprfun,2),spro);

    rpro = abs(sqrt((d1x.^2 + d1y.^2).^3)./ ...
                  (d1x.*d2y - d2x.*d1y)   );

    rprfun = csapi(spro,rpro);

    % tangential and normal directions
    tpro = [  d1x ;
              d1y ]; tpro = tpro ./ ([1;1]*(sqrt(tpro(1,:).^2 + tpro(2,:).^2)));
    npro = [ -d1y ;
              d1x ]; npro = npro ./ ([1;1]*(sqrt(npro(1,:).^2 + npro(2,:).^2)));
    tprfun{1} = csapi(spro,tpro(1,:)); tprfun{2} = csapi(spro,tpro(2,:));
    nprfun{1} = csapi(spro,npro(1,:)); nprfun{2} = csapi(spro,npro(2,:));

    %% Cut and re-interpolate profile
    % cut profile
    icb = find(xpro < xprcutlw,1,'first');
    ice = find(xpro < xprcutup,1,'last'); scut = spro([icb,ice]);

    scut(1) = fzero(@(s) fnval(xprfun,s)-xprcutlw,scut(1));
    scut(2) = fzero(@(s) fnval(xprfun,s)-xprcutup,scut(2));

    %scut_new = fzero(@(s) fnval(xprfun,s)-0.4,0.4)
    %Nn = nelx*(scut_new-scut(1))/(scut(2)-scut(1))

    % new spacing: weight function (wpr) for ds based on the profile curvature (1/rpr)
    wprfun = fnint(csapi(spro,abs(1./rpro).^(rexp)));
    wpro = fnval(wprfun,spro); wprinv = csapi(wpro,spro);

    wpr = linspace(fnval(wprfun,scut(1)),fnval(wprfun,scut(2)),nelx+1);

    % interpolation
    spr = fnval(wprinv,wpr);
    xpr = fnval(xprfun,spr);
    ypr = fnval(yprfun,spr);
    tpr =[fnval(tprfun{1},spr);
        fnval(tprfun{2},spr)];
    npr =[fnval(nprfun{1},spr);
        fnval(nprfun{2},spr)];
    rpr = fnval(rprfun,spr);
  	

	
	%% Compute upper/lower boundary from FLUENT/Morino data (streamlines)
	xic = linspace(xbox(1),xbox(2),nxbox);
	yic = linspace(ybox(1),ybox(2),nybox);
	
	[xxic,yyic] = meshgrid(xic,yic);
	
	% interpolate data on a rectangular box (xf,yf)
	bc = load(bcdata);
	
	switch bcsource
	    case 'FLUENT'
	        uuic = griddata(bc.xx,bc.yy,bc.uu,xxic,yyic,'cubic');
	        vvic = griddata(bc.xx,bc.yy,bc.vv,xxic,yyic,'cubic');
	        ppic = griddata(bc.xx,bc.yy,bc.pp,xxic,yyic,'cubic');
	end
	
	% streamline (upper boundaries)
	x0up = xpr(end) + doutup*npr(1,end);
	y0up = ypr(end) + doutup*npr(2,end);
	line = stream2(xic,yic,-uuic,-vvic,x0up,y0up); xupo = line{1}(:,1)'; yupo = line{1}(:,2)';
	
	% - tangential coordinate
	dx = xupo(2:end) - xupo(1:end-1);
	dy = yupo(2:end) - yupo(1:end-1);
	ds = sqrt(dx.^2+dy.^2); supo = cumsum([0 ds]);
	
	xupfun = csapi(supo,xupo); xupinv = csapi(xupo,supo);
	yupfun = csapi(supo,yupo);
	
	
	% streamline (lower boundaries)
	x0lw = xpr(1) + doutlw*npr(1,1);
	y0lw = ypr(1) + doutlw*npr(2,1);
	line = stream2(xic,yic,-uuic,-vvic,x0lw,y0lw); xlwo = line{1}(:,1)'; ylwo = line{1}(:,2)';
	
	% - tangential coordinate
	dx = xlwo(2:end) - xlwo(1:end-1);
	dy = ylwo(2:end) - ylwo(1:end-1);
	ds = sqrt(dx.^2+dy.^2); slwo = cumsum([0 ds]);
	
	xlwfun = csapi(slwo,xlwo); xlwinv = csapi(xlwo,slwo);
	ylwfun = csapi(slwo,ylwo);
	
	
	%% Connect and merge streamlines boundaries (freestream)
	
	% cut upper boundary (by interpolation)
	xup = linspace(xinup,xupo(1),ceil((xupo(1)-xinup)/dxbox)+1);
	sup = fnval(xupinv,xup);
	yup = fnval(yupfun,sup);
	iup = 3 + 0*xup;
	
	
	% cut lower boundary (by interpolation)
	xlw = linspace(xinlw,xlwo(1),ceil((xlwo(1)-xinlw)/dxbox)+1);
	slw = fnval(xlwinv,xlw);
	ylw = fnval(ylwfun,slw);
	ilw = 5 + 0*xlw;
	
	% flip lower boundary for consistency
	xlw = fliplr(xlw); ylw = fliplr(ylw); slw = fliplr(slw); slw = slw - slw(1);
	
	
	% inflow (spline from upper to lower streamline)
	xx = [         xlw(end-2:end) xup(1:3)           ];
	yy = [         ylw(end-2:end) yup(1:3)           ];
	ss = [slw(end-2:end)-slw(end) sup(1)-sup(1:3)+lin];
	
	sin = linspace(0,lin,ceil(lin/dxbox)+1);
	xin = interp1(ss,xx,sin,'spline');
	yin = interp1(ss,yy,sin,'spline');
	iin = 4 + 0*xin;
	
	
	% merge streamlines and inflow (freestream)
	xfso = [ xlw(1:end-1) xin(1:end) xup(2:end)];
	yfso = [ ylw(1:end-1) yin(1:end) yup(2:end)];
	ifso = [ ilw(1:end-1) iin(1:end) iup(2:end)];
	
	% - tangential coordinate
	dx = xfso(2:end) - xfso(1:end-1);
	dy = yfso(2:end) - yfso(1:end-1);
	ds = sqrt(dx.^2+dy.^2); sfso = cumsum([0 ds]);
	
	% freestream spline
	xfsfun = csapi(sfso,xfso);
	yfsfun = csapi(sfso,yfso); 
	
	% curvature radius (spline differentiation)
	d1x = fnval(fnder(xfsfun,1),sfso); d2x = fnval(fnder(xfsfun,2),sfso);
	d1y = fnval(fnder(yfsfun,1),sfso); d2y = fnval(fnder(yfsfun,2),sfso);
	
	rfso = abs(sqrt((d1x.^2 + d1y.^2).^3)./ ...
	              (d1x.*d2y - d2x.*d1y)   );
	
	rfsfun = csapi(spro,rpro);
	
	% tangential and normal directions
	tfso = -[  d1x ;
	           d1y ]; tfso = tfso ./ ([1;1]*(sqrt(tfso(1,:).^2 + tfso(2,:).^2)));
	nfso = -[ -d1y ;
	           d1x ]; nfso = nfso ./ ([1;1]*(sqrt(nfso(1,:).^2 + nfso(2,:).^2)));
	
	tfsfun{1} = csapi(sfso,tfso(1,:)); tfsfun{2} = csapi(sfso,tfso(2,:));
	nfsfun{1} = csapi(sfso,nfso(1,:)); nfsfun{2} = csapi(sfso,nfso(2,:));
	
	
	%% Map the free-stream boundary via profile normals
	% define tangential distance function
	fun = @(s,t,x0) t' * ([fnval(xfsfun,s);
	                       fnval(yfsfun,s)] - x0);
	
	% map (by finding intersections)
	sfs = zeros(size(spr));
	for i = 2:length(spr)
	    sfs(i) = fzero(@(s) fun(s,tpr(:,i),[xpr(i);
	                                        ypr(i)]),sfs(i-1));
	end
	xfs = fnval(xfsfun,sfs);
	yfs = fnval(yfsfun,sfs);
	
	tfs =[fnval(tfsfun{1},sfs);
	      fnval(tfsfun{2},sfs)];
	nfs =[fnval(nfsfun{1},sfs);
	      fnval(nfsfun{2},sfs)];
	rfs = fnval(rfsfun,sfs);
	
	ifs = interp1(sfso,ifso,sfs,'linear'); ifs = floor(ifs);
	
	map = (1-cos(linspace(0,1,nely+1)*pi/2));
	
	xxgr = zeros(nely+1,nelx+1);
	yygr = zeros(nely+1,nelx+1);
	
	%% curvature radius (spline differentiation)
	d1x = fnval(fnder(xprfun,1),spro); d2x = fnval(fnder(xprfun,2),spro);
	d1x = fnval(fnder(xprfun,1),spro); d2x = fnval(fnder(xprfun,2),spro);
	
	for i = 1:nelx+1
	    % wall-normal direction
	    xxgr(:,i) = interp1([0 1],[xpr(i) xfs(i)],map,'linear');
	    yygr(:,i) = interp1([0 1],[ypr(i) yfs(i)],map,'linear');
	
	end
	
	%% Grid-quality check
	xxel = zeros(nely,nelx); yyel = zeros(nely,nelx);
	dtel = zeros(nely,nelx); dnel = zeros(nely,nelx);
	arel = zeros(nely,nelx); skel = zeros(nely,nelx);
	
	for i = 1:nelx
	    for j = 1:nely
	        v1 = [xxgr(j  ,i  ); yygr(j  ,i  )];
	        v2 = [xxgr(j  ,i+1); yygr(j  ,i+1)];
	        v3 = [xxgr(j+1,i+1); yygr(j+1,i+1)];
	        v4 = [xxgr(j+1,i  ); yygr(j+1,i  )];
	
	        l1 = v2-v1; ll1 = norm(l1);
	        l2 = v3-v2; ll2 = norm(l2);
	        l3 = v4-v3; ll3 = norm(l3);
	        l4 = v1-v4; ll4 = norm(l4);
	
	        th1 = acos((l1(1)*l4(1) + l1(2)*l4(2))/(ll1*ll4));
	        th2 = acos((l2(1)*l1(1) + l2(2)*l1(2))/(ll2*ll1));
	        th3 = acos((l3(1)*l2(1) + l3(2)*l2(2))/(ll3*ll2));
	        th4 = acos((l4(1)*l3(1) + l4(2)*l3(2))/(ll4*ll3));
	
	        d1 = v3 - v1;
	        d2 = v4 - v2;
	
	% element grid (element center of mass)
	        xxel(j,i) = (v1(1) + v2(1) + v3(1) + v4(1))/4;
	        yyel(j,i) = (v1(2) + v2(2) + v3(2) + v4(2))/4;
	
	% aspect ratio
	        arel(j,i) = max([ll1 ll2 ll3 ll4])/min([ll1 ll2 ll3 ll4]);
	
	% skewness
	        skel(j,i) = max(abs([th1 th2 th3 th4] - pi/2))/(pi/2);
	
	% resolution
	        dnel(j,i) = max([ll2 ll4]);
	        dtel(j,i) = max([ll1 ll3]);
	
	
	    end
	end
	
	figure()
	mesh(xxgr, yygr, 0*xxgr)
	axis('equal')
	view(2)
  	
  	d.spr = spr;
  	d.xpr = xpr;
  	d.ypr = ypr;
  	d.x2 = xxgr(2,:);
  	d.y2 = yygr(2,:);
  	d.tpr = tpr;
  	d.xBC = xxgr(end,:);
  	d.yBC = yygr(end,:);
end



