clear all
close all
clc

warning off MATLAB:griddata:DuplicateDataPoints

% BEGIN Initializations of parameters
n=500;
m=1;
N=n*m;              %%% Dimensions of output layer (number of place cells)
NofSteps=200;      %%% Number of steps
envSize=10000;      %%% Size of environment
noiseAmpl=0.03;     %%% Amplitude of noise
nNV=1;              %%% Number of features we do not see and are predicted by previous
d=8;                %%% Number of features
NofRuns=30;         %%% Number of runs

ceateHPC=0;         %%% 0 - No cells; 1 - All cells; 2 - NonZero Cells 
showFig=1;

miu=0.01;
sigma=0.07;
threshold=0.5;

deltat=500; 
rand_ampl=deltat*0.20;
border=deltat+rand_ampl;

explorP=0.0;
maxGain=1;
alfa=0.7;
gama=0.7;
weightGain=0.0;

% w=rand(N,d); 
% w0=w;
% weightsGenerator

w=(1+exp((rand(N,8)-0.5)/(2*0.20^2))).^-1;
for di=1:8
    w(:,di)=w(:,di)-min(w(:,di));
    w(:,di)=w(:,di)/max(w(:,di));
end;

totalT=0;
% END Initializations of parameters

% BEGIN Initializations of parameters for odors
nOS=4;
odrSize=100;
g=0.1;
s=15;
a=0.01;

OS(1,:)=[odrSize*a, odrSize*a]+round(randn(1,2)*odrSize*g); 
OS(2,:)=[odrSize*a, odrSize*(1-a)]+round(randn(1,2)*odrSize*g); 
OS(3,:)=[odrSize*(1-a), odrSize*(1-a)]+round(randn(1,2)*odrSize*g); 
OS(4,:)=[odrSize*(1-a), odrSize*a]+round(randn(1,2)*odrSize*g); 

[Xodr,Yodr] = meshgrid(1:1:odrSize, 1:1:odrSize);

for kOS=1:nOS
   for i=1:length(Xodr)
       for j=1:length(Yodr)
           sx=s+i+5*sin(0.1*j);
           sy=s+j+5*sin(0.1*i);
           ODR(i,j,kOS)=exp(-((Xodr(i,j)-OS(kOS,1))^2/(2*sx^2) + (Yodr(i,j)-OS(kOS,2))^2/(2*sy^2)));
       end;
   end;
end;
for kOS=1:nOS
   ODR(:,:,kOS)=ODR(:,:,kOS)/max(max(ODR(:,:,kOS)));
end;
% END Initializations of parameters for odors

%% BEGIN Initializations of parameters for Q-learning

searchTime(1:NofRuns)=-1;
r=0;

weightN=weightGain*rand(1,N);
weightS=weightGain*rand(1,N);
weightE=weightGain*rand(1,N);
weightW=weightGain*rand(1,N);
weights_allN=weightN;
weights_allS=weightS;
weights_allE=weightE;
weights_allW=weightW;
weightNE=weightGain*rand(1,N);
weightNW=weightGain*rand(1,N);
weightSE=weightGain*rand(1,N);
weightSW=weightGain*rand(1,N);
weights_allNE=weightNE;
weights_allNW=weightNW;
weights_allSE=weightSE;
weights_allSW=weightSW;
%% END Initializations of weights for Q-learning

radius=6;
urine=zeros(100,100);

fxl=3000;		%%% location of food
fxr=5000;
fyd=6500;
fyu=8500;

if showFig==1
    figure('NumberTitle','off','Name', 'Latencies', 'Position', [232 276 560 410])
    uimenu('Label','Restart','Callback','navPFQLUv07')
    subplot(2,2,3)
    axis([0 100 0 100])
    hold on
    plot([fxl fxl fxr fxr fxl]/100,[fyd fyu fyu fyd fyd]/100,'r-')
end;

for ri=1:NofRuns
    M=NofSteps;         %%% Size of input data
    
    r=0;
    food=0;

    i1=1000;            %%% pradine x koordinate
    i2=1000;            %%% pradine y koordinate
    i1s=i1;
    i2s=i2;
    
    oldX=i1;
    oldY=i2;
    angle=0;

    t=0;
    j=0;
    nDivData=0;
    location=[];
    data=[];
    cellNumber=[];
    rate=[];
    
    kuri=[];
    kuri_pries=[];

    NS=rand-0.5;
    EW=rand-0.5;
    NE_SW=rand-0.5;
    NW_SE=rand-0.5;
    
    new_dir=0;
    old_dir=0;
    
    wb=waitbar(0,'Please wait... Creating SOM.','Position',[20 60 300 50]);
    while (t<deltat*NofSteps)&(food==0)
              
        totalT=totalT+1;
        j=j+1;
        
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
                angle=acos((newX-oldX)/sqrt((newX-oldX)^2+(newY-oldY)^2))*(180/pi);
            else
                angle=360-acos((newX-oldX)/sqrt((newX-oldX)^2+(newY-oldY)^2))*(180/pi);
            end;
        end;

        oldX=newX;
        oldY=newY;
        
        if nNV>0 & d>=4
            R=[];
            if angle<20 | angle>340
                R=[1];
            elseif angle>=20 & angle<=70
                if newX<envSize*0.25 & newY<envSize*0.25
                    R=[];
                else
                    R=[1 2];
                end;
            elseif angle>70 & angle<110
                R=[2];
            elseif angle>=110 & angle<=160
                if newX>envSize*0.75 & newY<envSize*0.25
                    R=[];
                else
                    R=[2 3];
                end;
            elseif angle>160 & angle<200
                R=[3];
            elseif angle>=200 & angle<=250
                if newX>envSize*0.75 & newY>envSize*0.75
                    R=[];
                else
                    R=[3 4];
                end;
            elseif angle>250 & angle<290
                R=[4];
            elseif angle>=290 & angle<=340
                if newX<envSize*0.25 & newY>envSize*0.75
                    R=[];
                else
                    R=[4 1];
                end;
            end;
            
            if j>1
                data(j,R)=data(j-1,R);
            else
                data(j,R)=0;
            end;
        end;
        
        if d==4
            data(j,5:8)=data(j,1:4);
        end;  
        
        %%%BEGIN Looking for smell 
        ui1=round(i1/100);
        ui2=round(i2/100);
        if ui1<1
            ui1=1;
        end;
        if ui2<1
            ui2=1;
        end;
        
        xs=ui1-radius;
        if xs<1
            xs=1;
        end;
        xe=ui1+radius;
        if xe>100
            xe=100;
        end;
        ys=ui2-radius;
        if ys<1
            ys=1;
        end;
        ye=ui2+radius;
        if ye>100
            ye=100;
        end;
        
        Smell=0;
        smell_pos=[];
        for xi=xs:xe
            for yi=ys:ye
                if urine(xi,yi)>0
                    Smell=1;
                    smell_pos=[smell_pos; xi yi urine(xi,yi)];
                    smell_pos=sortrows(smell_pos,3);
                end;
            end;
        end;   
        %%%END Looking for smell   
        
        kuri_pries=[];
        kuri_pries=kuri;
        kuri=[]; 
        
        v=[];
        for i=1:N

            v(i)=(sum([data(j,:)-w(i,:)].^2))/8;
            
            rate(j,i)=exp(-(v(i)^2)/(2*sigma^2));
            if rate(j,i)>=threshold
                kuri=[kuri, i];
            end; 
        end;
        
        [vm vmi]=min(v);
        w(vmi,:)=w(vmi,:)+miu*[data(j,:)-w(vmi,:)];     

       if length(kuri)<1
           disp('Cells do no fire...');
       end;
%         disp(length(kuri))
               
        %%% BEGIN Food finding
        if (i1>=fxl)&(i1<=fxr)&(i2>=fyd)&(i2<=fyu)
            food=1;
            r=1;
        end;
        %%% END Food finding 
        
       %%% BEGIN Motor learning
        nkiekis=length(kuri);
        if nkiekis~=0
            maksimumas=max([mean(weightN(kuri)) mean(weightS(kuri)) mean(weightE(kuri)) mean(weightW(kuri)) mean(weightNE(kuri)) mean(weightNW(kuri)) mean(weightSE(kuri)) mean(weightSW(kuri))]);
            vidQ=mean([mean(weightN(kuri)) mean(weightS(kuri)) mean(weightE(kuri)) mean(weightW(kuri)) mean(weightNE(kuri)) mean(weightNW(kuri)) mean(weightSE(kuri)) mean(weightSW(kuri))]);
        else
            maksimumas=0;
            vidQ=0;
        end;
        maksimumas=maxGain*maksimumas; 
        % disp(maksimumas)
             
        skiekis=length(kuri_pries);
        pos_dir=[];
        pos_dir=randperm(8);
        
        if i2>i2s & abs(i2-i2s)>rand_ampl & abs(i1-i1s)<=rand_ampl   
            weightN(kuri_pries)=weightN(kuri_pries)+alfa*(r+gama*maksimumas-weightN(kuri_pries));  
            NS=1;
            EW=0;
            NE_SW=0;
            NW_SE=0;
            old_dir=1;
            pos_dir(2)=0;
        elseif i2<i2s & abs(i2-i2s)>rand_ampl & abs(i1-i1s)<=rand_ampl  
            weightS(kuri_pries)=weightS(kuri_pries)+alfa*(r+gama*maksimumas-weightS(kuri_pries));
            NS=-1;
            EW=0;
            NE_SW=0;
            NW_SE=0;
            old_dir=-1;
            pos_dir(1)=0;
        elseif i1>i1s & abs(i1-i1s)>rand_ampl & abs(i2-i2s)<=rand_ampl
            weightE(kuri_pries)=weightE(kuri_pries)+alfa*(r+gama*maksimumas-weightE(kuri_pries)); 
            NS=0;
            EW=1;
            NE_SW=0;
            NW_SE=0;
            old_dir=2;
            pos_dir(4)=0;
        elseif i1<i1s & abs(i1-i1s)>rand_ampl & abs(i2-i2s)<=rand_ampl
            weightW(kuri_pries)=weightW(kuri_pries)+alfa*(r+gama*maksimumas-weightW(kuri_pries));  
            NS=0;
            EW=-1;
            NE_SW=0;
            NW_SE=0;
            old_dir=-2;
            pos_dir(3)=0;
        elseif (i1>i1s) & (i2>i2s) & abs(i1-i1s)>rand_ampl & abs(i2-i2s)>rand_ampl
            weightNE(kuri_pries)=weightNE(kuri_pries)+alfa*(r+gama*maksimumas-weightNE(kuri_pries));  
            NS=0;
            EW=0;
            NE_SW=1;
            NW_SE=0;
            old_dir=3;
            pos_dir(6)=0;
        elseif (i1<i1s) & (i2<i2s) & abs(i1-i1s)>rand_ampl & abs(i2-i2s)>rand_ampl
            weightSW(kuri_pries)=weightSW(kuri_pries)+alfa*(r+gama*maksimumas-weightSW(kuri_pries));
            NS=0;
            EW=0;
            NE_SW=-1;
            NW_SE=0;
            old_dir=-3;
            pos_dir(5)=0;
        elseif (i1<i1s) & (i2>i2s) & abs(i1-i1s)>rand_ampl & abs(i2-i2s)>rand_ampl 
            weightNW(kuri_pries)=weightNW(kuri_pries)+alfa*(r+gama*maksimumas-weightNW(kuri_pries)); 
            NS=0;
            EW=0;
            NE_SW=0;
            NW_SE=1;
            old_dir=4;
            pos_dir(8)=0;
        elseif (i1>i1s) & (i2<i2s) & abs(i1-i1s)>rand_ampl & abs(i2-i2s)>rand_ampl
            weightSE(kuri_pries)=weightSE(kuri_pries)+alfa*(r+gama*maksimumas-weightSE(kuri_pries)); 
            NS=0;
            EW=0;
            NE_SW=0;
            NW_SE=-1;
            old_dir=-4;
            pos_dir(7)=0;
        else
            if j>1
                disp('Could not specify direction...');
            end;
        end;

        
        if nkiekis~=0
            maxMax=max([max(weightN(kuri)) max(weightS(kuri)) max(weightE(kuri))...
                            max(weightW(kuri)) max(weightNE(kuri)) max(weightNW(kuri))...
                            max(weightSE(kuri)) max(weightSW(kuri))]);
        else
            maxMax=0;
        end; 
       
        if r==1
            if urine(ui1,ui2)<1
                urine(ui1,ui2)=urine(ui1,ui2)+0.005;
            end;
        end;
        
        if Smell==1 & maksimumas/vidQ>1.50
            if urine(ui1,ui2)<1
                urine(ui1,ui2)=urine(ui1,ui2)+0.005;
            end;
        end; 
        
        if showFig==1
            figure(findobj('Name','Latencies'))
            subplot(2,2,3)
            hold on
            if urine(ui1,ui2)>0
                ms=10+urine(ui1,ui2)*100;
                if ms>25
                    ms=25;
                end;
                plot(ui1,ui2,'g.','Markersize',ms)
            end;
        end;
        
% saving weights
%         weights_allN=[weights_allN;weightN];
%         weights_allS=[weights_allS;weightS];
%         weights_allE=[weights_allE;weightE];
%         weights_allW=[weights_allW;weightW];
%         weights_allNE=[weights_allNE;weightNE];
%         weights_allNW=[weights_allNW;weightNW];
%         weights_allSE=[weights_allSE;weightSE];
%         weights_allSW=[weights_allSW;weightSW];

        %%% END Motor learning
        
        if r==0
            %%%BEGIN New path
            if Smell==1  
                
                mI=size(smell_pos,1);
                i1n=smell_pos(mI,1)*100;
                i2n=smell_pos(mI,2)*100; 
                
                if i1n==i1 & i2n==i2 & mI>1
                    i1n=smell_pos(mI-1,1)*100;
                    i2n=smell_pos(mI-1,2)*100;     
                end;
                
                if abs(i1n-i1)<=rand_ampl | abs(i2n-i2)<=rand_ampl
                    if i2n>i2 & abs(i2n-i2)>rand_ampl    
                        new_dir=1;
                    elseif i2n<i2 & abs(i2n-i2)>rand_ampl
                        new_dir=-1;
                    elseif i1n>i1 & abs(i1n-i1)>rand_ampl
                        new_dir=2;
                    elseif i1n<i1 & abs(i1n-i1)>rand_ampl
                        new_dir=-2;
                    end;
                else
                    if (i1n>i1) & (i2n>i2)
                        new_dir=3;
                    elseif (i1n<i1) & (i2n<i2)
                        new_dir=-3;
                    elseif (i1n<i1) & (i2n>i2) 
                        new_dir=4;
                    elseif (i1n>i1) & (i2n<i2)
                        new_dir=-4;
                    end;
                end;
                
            else

                dwv=0;
                dwh=0;
                dwne_sw=0;
                dwnw_se=0;
                
                dwv=sum([weightN(kuri)-weightS(kuri)]);
                dwh=sum([weightE(kuri)-weightW(kuri)]);
                dwne_sw=sum([weightNE(kuri)-weightSW(kuri)]);
                dwnw_se=sum([weightNW(kuri)-weightSE(kuri)]);
                
                motion(totalT)=3;
                if dwv==0 & dwh==0 & dwne_sw==0 & dwnw_se==0  
                    if rand<0.25
                        dwv=rand-0.5;
                        dwh=rand-0.5;
                        dwne_sw=rand-0.5;
                        dwnw_se=rand-0.5;
                    else             
                        dwv=NS;
                        dwh=EW;
                        dwne_sw=NE_SW;
                        dwnw_se=NW_SE;
                    end;
                    motion(totalT)=0;
                else      
                    if rand<0.00 & j>1  
                        dwv=NS;
                        dwh=EW;
                        dwne_sw=NE_SW;
                        dwnw_se=NW_SE;
                        motion(totalT)=1;
                    end;        
                end;
                
                if rand<0.10 %explorP
                    dwv=rand-0.5;
                    dwh=rand-0.5;
                    dwne_sw=rand-0.5;
                    dwnw_se=rand-0.5;
                    motion(totalT)=2;
                end; 

                if abs(dwv)>abs(dwh) & abs(dwv)>abs(dwne_sw) & abs(dwv)>abs(dwnw_se)
                    i1n=i1+(rand-0.5)*2*rand_ampl;
                    i2n=i2+(deltat+(rand-0.5)*2*rand_ampl)*sign(dwv);
                    if sign(dwv)>0
                        new_dir=1;
                    else
                        new_dir=-1;
                    end;
                elseif abs(dwh)>abs(dwv) & abs(dwh)>abs(dwne_sw) & abs(dwh)>abs(dwnw_se)
                    i1n=i1+(deltat+(rand-0.5)*2*rand_ampl)*sign(dwh); 
                    i2n=i2+(rand-0.5)*2*rand_ampl;
                    if sign(dwh)>0
                        new_dir=2;
                    else
                        new_dir=-2;
                    end; 
                elseif abs(dwne_sw)>abs(dwh) & abs(dwne_sw)>abs(dwv) & abs(dwne_sw)>abs(dwnw_se)
                    i1n=i1+(deltat/sqrt(2)+(rand-0.5)*2*rand_ampl)*sign(dwne_sw);  
                    i2n=i2+(deltat/sqrt(2)+(rand-0.5)*2*rand_ampl)*sign(dwne_sw);
                    if sign(dwne_sw)>0
                        new_dir=3;
                    else
                        new_dir=-3;
                    end;
                elseif abs(dwnw_se)>abs(dwh) & abs(dwnw_se)>abs(dwv) & abs(dwnw_se)>abs(dwne_sw) 
                    i1n=i1-(deltat/sqrt(2)+(rand-0.5)*2*rand_ampl)*sign(dwnw_se); 
                    i2n=i2+(deltat/sqrt(2)+(rand-0.5)*2*rand_ampl)*sign(dwnw_se);
                    if sign(dwnw_se)>0
                        new_dir=4;
                    else
                        new_dir=-4;
                    end;
                end;  
            end;
            
            if old_dir==new_dir*(-1)
                if abs(old_dir)==1
                    i1n=i1+(rand-0.5)*2*rand_ampl;
                    i2n=i2+(deltat+(rand-0.5)*2*rand_ampl)*sign(old_dir);
                elseif abs(old_dir)==2
                    i1n=i1+(deltat+(rand-0.5)*2*rand_ampl)*sign(old_dir); 
                    i2n=i2+(rand-0.5)*2*rand_ampl;
                elseif abs(old_dir)==3
                    i1n=i1+(deltat/sqrt(2)+(rand-0.5)*2*rand_ampl)*sign(old_dir);  
                    i2n=i2+(deltat/sqrt(2)+(rand-0.5)*2*rand_ampl)*sign(old_dir);
                elseif abs(old_dir)==4
                    i1n=i1-(deltat/sqrt(2)+(rand-0.5)*2*rand_ampl)*sign(old_dir);  
                    i2n=i2+(deltat/sqrt(2)+(rand-0.5)*2*rand_ampl)*sign(old_dir);
                end;   
            end;
            
            out_of_border=0;
            % Correction of borders
            if i1n<=0 | i1n>=envSize | i2n<=0 | i2n>=envSize
                out_of_border=1;
                %disp(['out of border ' num2str(totalT)])
                if i1<=border & i2<=border
                    pos_dir([8,2,6,4,7])=0;
                elseif i1<=border & i2>border & i2<envSize-border 
                    pos_dir([6,4,7])=0;
                elseif i1<=border & i2>=envSize-border 
                    pos_dir([6,4,7,1,5])=0;
                elseif i1>border & i1<envSize-border & i2>=envSize-border 
                    pos_dir([7,1,5])=0;
                elseif i1>=envSize-border & i2>=envSize-border 
                    pos_dir([7,1,5,3,8])=0;
                elseif i1>=envSize-border & i2>border & i2<envSize-border 
                    pos_dir([5,3,8])=0;
                elseif i1>=envSize-border & i2<=border
                    pos_dir([5,3,8,2,6])=0;
                elseif i1>border & i1<envSize-border & i2<=border
                    pos_dir([8,2,6])=0;
                end;
            end;
            
            if out_of_border==1
                [maxV, maxI]=max(pos_dir);
                if maxI==1
                    i1n=i1;
                    i2n=i2+(deltat+(rand-0.5)*2*rand_ampl)*(1);
                elseif maxI==2
                    i1n=i1;
                    i2n=i2+(deltat+(rand-0.5)*2*rand_ampl)*(-1);
                elseif maxI==3
                    i1n=i1+(deltat+(rand-0.5)*2*rand_ampl)*(1);
                    i2n=i2;
                elseif maxI==4
                    i1n=i1+(deltat+(rand-0.5)*2*rand_ampl)*(-1);
                    i2n=i2;
                elseif maxI==5
                    i1n=i1+(deltat/sqrt(2)+(rand-0.5)*2*rand_ampl)*(1);  
                    i2n=i2+(deltat/sqrt(2)+(rand-0.5)*2*rand_ampl)*(1);
                elseif maxI==6
                    i1n=i1+(deltat/sqrt(2)+(rand-0.5)*2*rand_ampl)*(-1);  
                    i2n=i2+(deltat/sqrt(2)+(rand-0.5)*2*rand_ampl)*(-1);
                elseif maxI==7
                    i1n=i1-(deltat/sqrt(2)+(rand-0.5)*2*rand_ampl)*(1);
                    i2n=i2+(deltat/sqrt(2)+(rand-0.5)*2*rand_ampl)*(1);
                elseif maxI==8
                    i1n=i1-(deltat/sqrt(2)+(rand-0.5)*2*rand_ampl)*(-1);
                    i2n=i2+(deltat/sqrt(2)+(rand-0.5)*2*rand_ampl)*(-1);
                end;
            end;
            
            i1s=i1;
            i2s=i2;
            i1=i1n;
            i2=i2n;
            %%% END New path
        end;
        
        t=t+deltat;
        
        waitbar(j/M,wb);
    end;
    close(wb);
    
    M=j;
%      nDivData
%      MeanNoActiveCells=length(find(cellNumber~=0))/M
    
    if nDivData==0 
        searchTime(ri)=M;
    else
        break;
    end;
    
    if showFig==1
        figure(findobj('Name','Latencies'))
        subplot(2,2,1)
        cla
        hold on
        plot([fxl fxl fxr fxr fxl],[fyd fyu fyu fyd fyd],'r-')
        plot(location(:,1),location(:,2))
        axis([0 envSize 0 envSize])
        subplot(2,2,2)
        hold on
        stem(ri,searchTime(ri))
        set(gca,'xlim',[1 NofRuns])
    end;
    
    if ceateHPC==1
        createCellsAll;
    elseif ceateHPC==2
        createCells;
    end;

end;

if showFig==1
    figure(findobj('Name','Latencies'))
    subplot(2,2,4)
    mesh(urine')
end;

