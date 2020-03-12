%% compute flow quantities

%Compute flow quantities
for i=1:length(top)
    top(i).tauw=mu*top(i).dUdy(1);
    top(i).ut=sqrt(abs(top(i).tauw)/rho);
    top(i).Up=top(i).U/top(i).ut;
    top(i).Vp=top(i).V/top(i).ut;
    top(i).yp=top(i).yn*top(i).ut/nu;
    top(i).dUpdyp=nu/top(i).ut^2*top(i).dUdy;
end

for i=1:length(bottom)
    bottom(i).tauw=mu*bottom(i).dUdy(1);
    bottom(i).ut=sqrt(abs(bottom(i).tauw)/rho);
    bottom(i).Up=bottom(i).U/bottom(i).ut;
    bottom(i).Vp=bottom(i).V/bottom(i).ut;
    bottom(i).yp=bottom(i).yn*bottom(i).ut/nu;
    bottom(i).dUpdyp=nu/bottom(i).ut^2*bottom(i).dUdy;
end

%Compute RMSs
for i=1:length(top)
    top(i).uup=top(i).uu/top(i).ut^2;
    top(i).vvp=top(i).vv/top(i).ut^2;
    top(i).wwp=top(i).ww/top(i).ut^2;
    top(i).uvp=top(i).uv/top(i).ut^2;
end

for i=1:length(bottom)
    bottom(i).uup=bottom(i).uu/bottom(i).ut^2;
    bottom(i).vvp=bottom(i).vv/bottom(i).ut^2;
    bottom(i).wwp=bottom(i).ww/bottom(i).ut^2;
    bottom(i).uvp=bottom(i).uv/bottom(i).ut^2;
end

%Use diagnostic plot to compute BL edge and integral quantities
for i=30:79
    [top(i).deltap99,top(i).Uinfp,top(i).i99,top(i).H,top(i).deltasp, ...
        top(i).thetap]=diagnostic(top(i).yp,top(i).Up,sqrt(top(i).uup));
    top(i).delta99=top(i).deltap99*nu/top(i).ut;
    top(i).deltas=top(i).deltasp*nu/top(i).ut;
    top(i).theta=top(i).thetap*nu/top(i).ut;
    top(i).Uinf=top(i).Uinfp*top(i).ut;
    top(i).Ret=top(i).deltap99;
    top(i).Reds=top(i).Uinf*top(i).deltas/nu;
    top(i).Reth=top(i).Uinf*top(i).theta/nu;
    
end

for i=29:79
    [bottom(i).deltap99,bottom(i).Uinfp,bottom(i).i99,bottom(i).H,bottom(i).deltasp, ...
        bottom(i).thetap]=diagnostic(bottom(i).yp,bottom(i).Up,sqrt(bottom(i).uup));
    bottom(i).delta99=bottom(i).deltap99*nu/bottom(i).ut;
    bottom(i).deltas=bottom(i).deltasp*nu/bottom(i).ut;
    bottom(i).theta=bottom(i).thetap*nu/bottom(i).ut;
    bottom(i).Uinf=bottom(i).Uinfp*bottom(i).ut;
    bottom(i).Ret=bottom(i).deltap99;
    bottom(i).Reds=bottom(i).Uinf*bottom(i).deltas/nu;
    bottom(i).Reth=bottom(i).Uinf*bottom(i).theta/nu;
    
end


%Compute flow quantities
for i=1:length(top)
    top(i).Ue=top(i).U(top(i).i99);
    top(i).Ve=top(i).V(top(i).i99);
end

for i=1:length(bottom)
    bottom(i).Ue=bottom(i).U(bottom(i).i99);    
    bottom(i).Ve=bottom(i).V(bottom(i).i99);    
end


%Assign pressure and pressure derivative
for j=1:2*np-1
    if j<=np
        for i=1:ln
            top(j).P(i,1)=P(i,j);
            top(j).dPdx(i,1)=D7(i,j);
        end
    else
        for i=1:ln
            bottom(j-np+1).P(i,1)=P(i,j);
            bottom(j-np+1).dPdx(i,1)=D7(i,j);
        end
    end
end

for j=1
    for i=1:ln
        bottom(j).P(i,1)=P(i,j);
        bottom(j).dPdx(i,1)=D7(i,j);
    end
end

%Assign RMS of pressure 
for j=1:2*np-1
    if j<=np
        for i=1:ln
            top(j).pp(i,1)=pp(i,j);
            top(j).prms(i,1)=prms(i,j);
        end
    else
        for i=1:ln
            bottom(j-np+1).pp(i,1)=pp(i,j);
            bottom(j-np+1).prms(i,1)=prms(i,j);
        end
    end
end

for j=1
    for i=1:ln
        bottom(j).pp(i,1)=pp(i,j);
        bottom(j).prms(i,1)=prms(i,j);
    end
end

%Compute pressure gradient parameter
for i=30:79
    top(i).Delta=top(i).Uinfp*top(i).deltas;
    top(i).dUinfdx=top(i).dUdx(top(i).i99);
    top(i).beta=top(i).deltas/top(i).tauw*top(i).dPdx(top(i).i99);
end

for i=29:79
    bottom(i).Delta=bottom(i).Uinfp*bottom(i).deltas;
    bottom(i).dUinfdx=bottom(i).dUdx(bottom(i).i99);
    bottom(i).beta=bottom(i).deltas/bottom(i).tauw*bottom(i).dPdx(bottom(i).i99);
end

%compute cf:
for i=30:79
    if top(i).ut>=0
        top(i).Cf=2*(top(i).ut/top(i).Uinf)^2;
    else
        top(i).Cf=-2*(top(i).ut/top(i).Uinf)^2;
    end
end

for i=30:79
    bottom(i).Cf=2*(bottom(i).ut/bottom(i).Uinf)^2;
end


%compute Cp
pinf = max(max(P))-0.5;
rhoinf = rho;
Uinf = 1;

for i = 30:79
    top(i).Cp = 2*(top(i).P(1)-pinf)/(rhoinf*Uinf^2);
    bottom(i).Cp = 2*(bottom(i).P(1)-pinf)/(rhoinf*Uinf^2);
end
