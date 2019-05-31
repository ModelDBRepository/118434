clear all
close all
clc

warning off MATLAB:griddata:DuplicateDataPoints

% BEGIN Initializations of parameters

NofSteps=200;      %%% Number of steps
envSize=10000;      %%% Size of environment
NofRuns=60;         %%% Number of runs

deltat=500; 
rand_ampl=deltat*0.20;
border=deltat+rand_ampl;

radius=6;
urine=zeros(100,100);

explorP=0.0;

totalT=0;

urine=zeros(100,100);

searchTime(1:NofRuns)=-1;

fxl=3000;		%%% location of food
fxr=5000;
fyd=6500;
fyu=8500;

showFig=1;
% END Initializations of parameters

if showFig==1
    figure('NumberTitle','off','Name', 'Latencies', 'Position', [232 237 560 441])
    uimenu('Label','Restart','Callback','navUrineBasedv04')
    subplot(2,2,3)
    axis([0 100 0 100])
    hold on
    plot([fxl fxl fxr fxr fxl]/100,[fyd fyu fyu fyd fyd]/100,'r-')
end;

for ri=1:NofRuns
    
    i1=1000;            %%% pradine x koordinate
    i2=1000;            %%% pradine y koordinate
    i1s=i1;
    i2s=i2;
    
    M=NofSteps;         %%% Number of steps
    
    food=0;

    t=0;
    j=0;
    location=[];
    
    NS=rand-0.5;
    EW=rand-0.5;
    NE_SW=rand-0.5;
    NW_SE=rand-0.5;
    
    new_dir=0;
    old_dir=0;
    
    wb=waitbar(0,'Please wait...','Position',[20 60 300 50]);
    while (t<deltat*NofSteps)&(food==0)
    
        totalT=totalT+1;
        j=j+1;
        
        location(j,:)=[i1 i2];
        ui1=round(i1/100);
        ui2=round(i2/100);
        if ui1<1
            ui1=1;
        end;
        if ui2<1
            ui2=1;
        end;
        
        %%%BEGIN Looking for smell 
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
        
        %%% BEGIN Food finding
        if (i1>=fxl)&(i1<=fxr)&(i2>=fyd)&(i2<=fyu)
            food=1;
            r=1;
            
            if urine(ui1,ui2)<1
                urine(ui1,ui2)=urine(ui1,ui2)+0.005;
            end;
        end;
        %%% END Food finding 

         if Smell==1
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
         
        pos_dir=[];
        pos_dir=randperm(8);
        
        if i2>i2s & abs(i2-i2s)>rand_ampl & abs(i1-i1s)<=rand_ampl   
            NS=1;
            EW=0;
            NE_SW=0;
            NW_SE=0;
            old_dir=1;
            pos_dir(2)=0;
        elseif i2<i2s & abs(i2-i2s)>rand_ampl & abs(i1-i1s)<=rand_ampl  
            NS=-1;
            EW=0;
            NE_SW=0;
            NW_SE=0;
            old_dir=-1;
            pos_dir(1)=0;
        elseif i1>i1s & abs(i1-i1s)>rand_ampl & abs(i2-i2s)<=rand_ampl
            NS=0;
            EW=1;
            NE_SW=0;
            NW_SE=0;
            old_dir=2;
            pos_dir(4)=0;
        elseif i1<i1s & abs(i1-i1s)>rand_ampl & abs(i2-i2s)<=rand_ampl
            NS=0;
            EW=-1;
            NE_SW=0;
            NW_SE=0;
            old_dir=-2;
            pos_dir(3)=0;
        elseif (i1>i1s) & (i2>i2s) & abs(i1-i1s)>rand_ampl & abs(i2-i2s)>rand_ampl
            NS=0;
            EW=0;
            NE_SW=1;
            NW_SE=0;
            old_dir=3;
            pos_dir(6)=0;
        elseif (i1<i1s) & (i2<i2s) & abs(i1-i1s)>rand_ampl & abs(i2-i2s)>rand_ampl
            NS=0;
            EW=0;
            NE_SW=-1;
            NW_SE=0;
            old_dir=-3;
            pos_dir(5)=0;
        elseif (i1<i1s) & (i2>i2s) & abs(i1-i1s)>rand_ampl & abs(i2-i2s)>rand_ampl 
            NS=0;
            EW=0;
            NE_SW=0;
            NW_SE=1;
            old_dir=4;
            pos_dir(8)=0;
        elseif (i1>i1s) & (i2<i2s) & abs(i1-i1s)>rand_ampl & abs(i2-i2s)>rand_ampl
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
            % Koordinates skaiciavimas 
            
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
%             disp(['out of border ' num2str(totalT)])
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
%             disp(pos_dir)
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
        
        t=t+deltat;
        
        waitbar(j/M,wb);
    end;
    close(wb);
    
    M=j;
    
    searchTime(ri)=M;
    
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
    
% pause(0.5)
end;

if showFig==1
    figure(findobj('Name','Latencies'))
    subplot(2,2,4)
    mesh(urine')
end;