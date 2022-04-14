
%Author: Alka Potdar
% Run using MATLAB R2009b

%%Predicting protein production with and without translational silencing by
%%GAIT system (G) in the presence of mini-EPRS (M)
% we predict the protein translation from 0 to 24 hours (which is time
% after IFN-gamma (IFN) treatment.

%%%%%%SPECIES concentration%%%%%%%
% protein (y2)
% GAIT complex (G) and mini EPRS (M)
% m-rna (y1) free form not bound to anything (G or M)
% m-rna in the active form is either y1 or bound to M (My1)
% Inactive m-rna (Gy1)
% IFN -gamma (ifn)

%%%%%%03/31/11
%%%%%%incorporating negative feeback to GATE element by DAPK-ZIPK proteins
%%%%%%\resulting in degradation of GAIT element after 32 hrs


%%%%%M-RNA SATURATED AROUND 24, no half life in the model

%%%%%%note there are 2 ways to model transcriptional surge due ifn
%%%%%%bolus%%%%
%%%%%041111: note there is a cut-off value of ifn0 that is required for
%%%%%mrna surge as well as GAIT production (which is 20 uM)

%%%%%%%%
%%%% values of forward and reverse rate constants (from biocore binding
%%%% expt)

%%%%%here we assume that the RT-PCR and protein synthesis rates are in nM
%%%%%range
% Ka_eprs=2.65*10^5  (1/Ms) %%%convert to (nM.h)-1 = approx 0.75!!!
% Kr_eprs=3.07*10^-3 (1/s)
% Ka_eprsn1=2.03*10^5 (1/Ms)
% Kr_eprsn1=5.13*10^-5 (1/s)
%%%%%%%%%%



%%%%%%Performing parameter optimization using lsqcurvefit using the
%%%%%%experimnetally know total mrna and protein synthesis rate profiles

%%%%%%parameters involved: G,M,kd2

function gait_v10_lsq_try
clear
clc

global ifnogam ifn0 ydata1 ydata2 final_mrna final_dpdt



%least square fit for 24/8 hour data
%ydata1=importdata('mrna_24/8');


%%%%mrna for 8 hour %%%%%%

% ydata1=[1.0000
% 0.83000
% 0.68000
% 1.3600
% 1.9100
% 2.2700
% 2.0600
% 2.7500
% 2.3800
% 3.0700
% 3.3400
% 4.1100
% ];

%%%%mrna for 24 hour %%%%%%

ydata1=[1
1.44
0.97
1.75
2.33
3.36
4.5
4.23
6.41
7.67
7.84
8.06];




%ydata2=importdata('dprot for 8h');
% ydata2=[1.0000
% 5.5283
% 10.048
% 15.722
% 18.601
% 26.200
% 31.369
% 36.775
% 31.041
% 37.853
% 39.185
% 37.529
% ];



%ydata2 for 24 hour%%%%

ydata2=[1
1.492877412
1.404164144
2.778810537
5.264814697
3.58699764
4.948568599
3.843601494
3.795281407
4.113178139
5.010314971
5.697637282];

ydata1=transpose(ydata1);
ydata2=transpose(ydata2);

%ifnogam=[0;10;20;30;50;100;200;500];
%xdata=transpose(ifnogam);

ifnogam=0.0001:25:501;

%options=optimset('TOLX',1e-5);

%[x1,RESNORM,RESIDUAL,EXITFLAG,OUTPUT]=lsqcurvefit(@rna,[500 1 0.2],xdata,ydata1,[100 0.01 0.01],[2000 20 10],options)


%[x,RESNORM,RESIDUAL,EXITFLAG,OUTPUT]=lsqcurvefit(@before,[G0 M0 kd2 k3 scale],xdata,ydata,[0.000001 0.000001],[1000 1000])
%[x2,RESNORM,RESIDUAL,EXITFLAG,OUTPUT]=lsqcurvefit(@prot,[500 1 0.2],xdata,ydata2,[100 0.01 0.01],[8000 20 10],options)

%[x,RESNORM,RESIDUAL,EXITFLAG,OUTPUT]=lsqcurvefit(@before,[k3],xdata,ydata,[0.000001 0.000001],[1000 1000])

 
%function mrna=rna(x1,xdata)
%global ifnogam ifn0 ydata1 final_mrna
%mrna=zeros(1,8);
%dpdt=zeros(1,8);
for i=1:21
   
ifn0=ifnogam(i);
        % time after IFN treatment in the experiments, tout in hours
        t0=0.0;
        
        % change tf according to time point either 8 or 24
        tf=24.0;
        step=1.0;
        tout=t0:1:tf;

       
        % initial concentration
        y1(1)=1;
        y2(1)=0;
        My1(1)=0;
        Gy1(1)=0;


        u(1)=y1(1);
        u(2)=Gy1(1);
        u(3)=My1(1);
        u(4)=y2(1);

        u0=[u(1) u(2) u(3) u(4)];


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [t,u] = ode15s(@extra,tout,u0);



        tnew=t0:1:tf-1;


        %rate of protein change
        new=diff(u(:,4));

        %total mrna
        u(:,1)+u(:,3)+u(:,2);


        %u(8,1)+u(8,2)+u(8,3)
        %new(8)
        total_mrna=u(tf,1)+u(tf,2)+u(tf,3);
        proteinsyn=new(tf);
        
        
        if  ifn0==0.0
            control=new(tf)
        end
        
        mrna(i)=u(tf,1)+u(tf,2)+u(tf,3);
        dpdt(i)=new(tf);
  

% figure(1)
% plot(t,u(:,1)+u(:,3)+u(:,2));
% xlabel('time after ifn bolus (hours)')
% ylabel('concentration (AU)')
% legend('mrna-unbound','Gy1-Gait bound','My1, mini-eprs bound','total mrna')
%   
% hold on;
%         
% figure(2)
% plot(tnew,diff(u(:,4)),'b');
% xlabel('time after ifn bolus (hours)');
% ylabel('rate of protein synthesis,dP/dt (AU/hour)')
% %title('red:no GAIT blue:GAIT, no M-EPRS green: GAIT, with M-EPRS')
% hold on;

end
transpose(ifnogam)
final_mrna=mrna
final_dpdt=dpdt



figure(3)
plot(final_mrna,final_dpdt,'b*-',ydata1,ydata2,'r*')
xlabel('total mrna')
ylabel('protein synthesis')
legend('predicted','experimental');
hold on;

transpose(final_dpdt)
transpose(mrna)


function ut=extra(t,u)

            global ifn0
            % rate constants

            k1= 0.75; % forward rate constant for y1 + G --> Gy1 IN (conc.time)-1
            kd1=10.6; % dissociation constant for GAIT complex protein in (conc)-1 spec in (nM)-1
            %k2= 0.4;%forward rate constant for y1 + M --> My1 IN (conc.time)-1
            k2=0.75; % assuming it same as that of full length EPRS
            kd2=0.2; % dissociation constant for mini-EPRS  protein in (conc)-1 spec in (nM)-1
            k3=2; % rate of translation of mrna to protein, (time around 30 minutes)so 2 per hour 
            %v3=0.1; %maximal rate for proten synthesis
            k=0.08;% rate constant for transcriptional surge by ifn
            scale=11.4;%scaling for mrna and ifn0 calibration
            scale2=5.0;
            
            %tmrna=8;%in hrs 
            %y1_ss=800; %units/ml
            
            y1_ss=800; %units/ml


            % concentration of of GAIT comples (G) for translational inhibition is 25%
            % of intial mrna (1000) using info that M is 5% eprs

   
            %tau=20/60; % experimentally known half-life time of ifn-gamma in hours
            %ifn=ifn0*exp(-t/tau);

            cutoff=0.0;

            if t<=32
                if ifn0 >cutoff %cut-off ifn0 for gait activation
                G=y1_ss*heaviside(t-14)-u(2);% y1ss is the maximum value of mrna at ss after ifn stimulation
                % total GAIT protein (bound (u(2)) + unbound remains constant

                %G=0;
                else
                G=0.0;
                end
            else
                G=0.0;
                %G=y1_ss*exp(-t*10);
            end


            if ifn0 > cutoff %cut-off ifno for mrna transcriptional surge
                mrna_max=ifn0*scale/(200+ifn0); % hyperbolic relation between mrna_max(t=24) versus ifno
            %elseif or(ifn0>0,ifn0<20)
             %   k_base=0.025;
              %  mrna_max=k_base*ifn0;
            else
                 mrna_max=0;
            end


            % total mini-eprs remains constant, M_initial=2
            initial_m=0;

            M=initial_m-u(3);
            %M=0;

            if t<=32
            ut(1,1)=k*mrna_max*exp(-t*k)+ k1*(-G*u(1) + kd1*u(2)) + k2*(-M*u(1) + kd2*u(3)); 
            ut(2,1)=k1*(u(1)*G - kd1*u(2)); %inactive m-rna
            else
            ut(1,1)=k*mrna_max*exp(-t*k)+ k2*(-M*u(1) + kd2*u(3)); 
            ut(2,1)=0; %inactive m-rna
            end
            %ut(2,1)=k1*(u(1)*G - kd1*u(2)); %inactive m-rna
            ut(3,1)=k2*(M*u(1) - kd2*u(3)); % My1 active bound mrna
            
            if ifn0 > cutoff
            ut(4,1)=scale2*k3*(u(1)+u(3));
            else
            ut(4,1)=k3*(u(1)+u(3));
            end
            
            
            





