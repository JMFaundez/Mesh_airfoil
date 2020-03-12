function [ int ] = surfintegral( x,y,f )
%surfintegral Summary of this function goes here


deltax=diff(x);
deltay=diff(y);
f2=(f(1:end-1)+f(2:end))/2.0;

npts=length(deltax);
int=0;

for ipts=1:npts
    deltal=sqrt(deltax(ipts)^2+deltay(ipts)^2);
    int=int+deltal*f2(ipts);
end

end

