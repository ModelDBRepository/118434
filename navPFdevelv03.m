

clear all
% close all
clc


fd=fopen('rates.dat','w+');

t0 = clock;

warning off MATLAB:griddata:DuplicateDataPoints

% BEGIN Pradiniai nustatymai
n=500;
m=1;
N=n*m;             %%% Dimensions of output layer
M=5000;            %%% Size of input data
envSize=10000;     %%% Size of environment
nNV=1;             %%% Number of features we do not see and are predicted by previous
d=8;               %%% Number of features
noiseAmpl=0.03;    %%% Noise amplitude
alfa=0;            %%% View angle at starting position
ceateHPC=0;        %%% 0 - No cells; 1 - All cells; 2 - NonZero Cells 

miu=0.01;
ff(1:N)=0;
LR(1:N)=1;

sigma=0.07;

spe_1b=1;          %%% greiciai x kryptimi
spe_2b=1;          %%% greiciai x kryptimi
deltat=100;         %%% nuo 1 iki 100
i1=1000;           %%% pradine x koordinate
i2=1000;           %%% pradine y koordinate

oldX=i1;
oldY=i2;
% END Pradiniai nustatymai

% BEGIN Initializations of parameters for odors

nOS=4;
odrSize=100;
g=0.1;
s=15;
a=0.01;

[Xodr,Yodr] = meshgrid(1:1:odrSize, 1:1:odrSize);

OS(1,:)=[odrSize*a, odrSize*a]+round(randn(1,2)*odrSize*g); 
OS(2,:)=[odrSize*a, odrSize*(1-a)]+round(randn(1,2)*odrSize*g); 
OS(3,:)=[odrSize*(1-a), odrSize*(1-a)]+round(randn(1,2)*odrSize*g); 
OS(4,:)=[odrSize*(1-a), odrSize*a]+round(randn(1,2)*odrSize*g); 
for kOS=1:nOS
    for i=1:length(Xodr)
        for j=1:length(Yodr)
            sx=s+i+5*sin(0.1*i);
            sy=s+j-5*sin(0.1*j);
            ODR(i,j,kOS)=exp(-((Xodr(i,j)-OS(kOS,1))^2/(2*sx^2) + (Yodr(i,j)-OS(kOS,2))^2/(2*sy^2)));
        end;
    end;
end;
for kOS=1:nOS
    ODR(:,:,kOS)=ODR(:,:,kOS)/max(max(ODR(:,:,kOS)));
end;

% END Initializations of parameters for odors


w0=(1+exp((rand(N,8)-0.5)/(2*0.20^2))).^-1;
for di=1:d
    w0(:,di)=w0(:,di)-min(w0(:,di));
    w0(:,di)=w0(:,di)/max(w0(:,di));
end;

% save 'weights.dat' -ASCII w0
% w0=load('weights.dat');

w=w0;
% figure, axis([0 1 0 1]), hold on, plot(w(:,1),w(:,2),'r*')
% figure, hist(w0)

nDivData=0;

wb=waitbar(0,'Please wait...');
for j=1:M
      
    location(j,:)=[i1 i2];      
    data(j,1:2)=location(j,1:2)/envSize;
    
    if d>=4
        data(j,3:4)=1-data(j,1:2);
        data(j,5:8)=zeros(1,4);
    end;         
       
    if d==8     
        I1=round(i1/100);
        I2=round(i2/100);
        if I1==0
            I1=1;
        end;
        if I2==0
            I2=1;
        end;
        data(j,5:8)=ODR(I1,I2,1:4);
    end;
    

    noise1=[rand(1,4)-0.5]*2*noiseAmpl;
    noise2=[rand(1,4)-0.5]*2*noiseAmpl;
    
    if d>=4
        data(j,1:4)=data(j,1:4)+data(j,1:4).*noise1;
    end;
    
    if d==8
        data(j,5:8)=data(j,5:8)+[1-data(j,5:8)].*noise2;
    end;
        
    for di=1:8
        if data(j,di)<0
            data(j,di)=0;
        elseif data(j,di)>1
            data(j,di)=1;
        end;
    end;
    
    newX=i1;
    newY=i2;
    if (newX-oldX)==0 & (newY-oldY)==0
    else
        if newY-oldY>0
            alfa=acos((newX-oldX)/sqrt((newX-oldX)^2+(newY-oldY)^2))*(180/pi);
        else
            alfa=360-acos((newX-oldX)/sqrt((newX-oldX)^2+(newY-oldY)^2))*(180/pi);
        end;
    end;

    oldX=newX;
    oldY=newY;
    
    if nNV>0 & d>=4
        R=[];
        if alfa<20 | alfa>340
                R=[1];
        elseif alfa>=20 & alfa<=70
            if newX<envSize*0.25 & newY<envSize*0.25
                R=[];
            else
                R=[1 2];
            end;
        elseif alfa>70 & alfa<110
                R=[2];
        elseif alfa>=110 & alfa<=160
            if newX>envSize*0.75 & newY<envSize*0.25
                R=[];
            else
                R=[2 3];
            end;
        elseif alfa>160 & alfa<200
                R=[3];
        elseif alfa>=200 & alfa<=250
            if newX>envSize*0.75 & newY>envSize*0.75
                R=[];
            else
                R=[3 4];
            end;
        elseif alfa>250 & alfa<290
                R=[4];
        elseif alfa>=290 & alfa<=340
             if newX<envSize*0.25 & newY>envSize*0.75
                R=[];
            else
                R=[4 1];
            end;
        end;
        
        if j > 1
            data(j,R)=data(j-1,R);
        else
            data(j,R)=0;
        end;
    end;
    
    if d==4
        data(j,5:8)=data(j,1:4);
    end;
    
    
    v=[];
    rate=[];
    for i=1:N
        v(i)=(sum([data(j,:)-w(i,:)].^2))/8;
        rate(i)=exp(-(v(i)^2)/(2*sigma^2));
    end;
    fprintf(fd,'\t%f',rate);
    fprintf(fd,'\n');
    
    [vm vmi]=min(v);
    w(vmi,:)=w(vmi,:)+miu*[data(j,:)-w(vmi,:)];
      
    
    spe_1=spe_1b+round(rand*6-3); 
    spe_2=spe_2b+round(rand*6-3);
    
    if spe_1==0 
        if rand<=0.5
            spe_1=-1;
        else
            spe_1=1;
        end;
    end;
    if spe_2==0 
        if rand<=0.5
            spe_2=-1;
        else
            spe_2=1;
        end;
    end;
    
    i1=i1+deltat*spe_1;
    i2=i2+deltat*spe_2;

    % correction of borders
    if i1>envSize
         i1=envSize;
         spe_1b=-1;
    end;
    if i2>envSize
        i2=envSize;
        spe_2b=-1;
    end;
    if i1<0
         i1=0;
         spe_1b=1;
    end;
    if i2<0;
        i2=0;
        spe_2b=1;
    end; 
     
    waitbar(j/M,wb);
end;
close(wb);

fclose(fd);

% nDivData
% MeanNoActiveCells=length(find(cellNumber>0))/M

% figure
% plot(location(:,1),location(:,2))
% axis([0 envSize 0 envSize])

rate=[];
rate=load('rates.dat');

elapseTime=etime(clock,t0)/60

createCellsAll
plotPFM(1,10,10,X,Y,HPC,0.01)
plotPFM(101,10,10,X,Y,HPC,0.01)
plotPFM(201,10,10,X,Y,HPC,0.01)
plotPFM(301,10,10,X,Y,HPC,0.01)
plotPFM(401,10,10,X,Y,HPC,0.01)