function [data,EL,header,tag,N,nel,nfields,status,t] = readnek_reshape(endian,fname)

% This function reads binary data from the nek5000 output files

% open the file
[fid,message] = fopen(fname,'r',['ieee-' endian]);
if fid == -1, disp(message), return, end

% read the header
header = fread(fid,132,'*char')';
%%%%
std = header(6);
if str2num(std)==4
    precision='*float32';
else
    precision='*float64';
end
%%%%
N   = str2double(header(8:9));
%%%%
nz  = header(14:15);
if str2num(nz)==1
    dim =2;
else
    dim=3;
end
%%%%
nel = str2double(header(20:30));
%%%%
t = str2double(header(39:60));
%%%%
tag   = fread(fid,1,'*float32');
%%%% elements indices (not sorted) 
lglel = fread(fid,nel,'*int32');
%%%%
fields=header(84:87);
var=zeros(1,4);
if(~isempty(strfind(fields,'X')))
    var(1)=1;
end
if(~isempty(strfind(fields,'U')))
    var(2)=1;
end
if(~isempty(strfind(fields,'P')))
    var(3)=1;
end
(~isempty(strfind(fields,'T')))
if(~isempty(strfind(fields,'T')))
    var(4)=1;
end
nfields = var(1)*dim+var(2)*dim+var(3)+var(4);
%%%%
% tic
% %  read the data
% for ivar=1:sum(var(1:2))
%     for iel=1:nel
%         data(lglel(iel),:,(ivar-1)*dim+1:ivar*dim)=reshape(fread(fid,dim*(N^dim),precision),N^dim,dim);
%     end
% end
% toc
if dim ==3
    data = cell(N,N,N,nfields);
    if var(1)==1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% read coordinate points (x,y,z)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for iel=1:nel
            dum = fread(fid,dim*(N^dim),precision);
            EL(lglel(iel)).GLL          = reshape(dum,N^dim,dim);
            data{lglel(iel)}(:,:,:,1:3) = reshape(dum,N,N,N,dim);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% read velocities (u,v,w)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iel=1:nel
        dum = fread(fid,dim*(N^dim),precision);
        EL(lglel(iel)).VEL          = reshape(dum,N^dim,dim);
        data{lglel(iel)}(:,:,:,4:6) = reshape(dum,N,N,N,dim);
    end
    if var(3)==1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% read pressure p
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for iel=1:nel
            dum = fread(fid,N^dim,precision);
            EL(lglel(iel)).P          = reshape(dum,N^dim,1);
            data{lglel(iel)}(:,:,:,7) = reshape(dum,N,N,N,1);
        end
    end
elseif dim==2
    data = cell(N,N,nfields);
    if var(1)==1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% read coordinate points (x,y,z)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for iel=1:nel
            dum = fread(fid,dim*(N^dim),precision);
            EL(lglel(iel)).GLL        = reshape(dum,N^dim,dim);
            %data{lglel(iel)}(:,:,1:2) =reshape(dum,N,N,dim);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% read velocities (u,v,w)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iel=1:nel
        dum = fread(fid,dim*(N^dim),precision);
        EL(lglel(iel)).VEL        = reshape(dum,N^dim,dim);
        %data{lglel(iel)}(:,:,3:4) = reshape(dum,N,N,dim);
    end
    if var(3)==1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% read pressure p
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for iel=1:nel
            dum = fread(fid,N^dim,precision);
            EL(lglel(iel)).P        = reshape(dum,N^dim,1);
            %data{lglel(iel)}(:,:,5) = reshape(dum,N,N,1);
        end
    end
    if var(4)==1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% read pressure T
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for iel=1:nel
            dum = fread(fid,N^dim,precision);
            EL(lglel(iel)).T        = reshape(dum,N^dim,1);
            %data{lglel(iel)}(:,:,5) = reshape(dum,N,N,1);
        end
    end
    
end

status = fclose(fid);

return
