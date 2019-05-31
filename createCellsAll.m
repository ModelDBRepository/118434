
warning off MATLAB:griddata:DuplicateDataPoints

ti = 0:0.01:1; 
[X,Y] = meshgrid(ti,ti);

HPC=[];
wb=waitbar(0,'Please wait... Creating Place Fields.');

int=1:M;

for i=1:N
    HPC(:,:,i)=griddata(location(int,1)/envSize, location(int,2)/envSize, rate(int,i), X, Y);
    waitbar(i/N,wb);
end;
close(wb);

HPC(find(isnan(HPC)==1))=0;