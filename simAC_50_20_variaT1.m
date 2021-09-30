clear;

%vector_variation=[50,100,150,200,250,300,350,400,450,500]/(30*10);
vector_variation=[0.5,1.0,1.5,2.0,2.5,3.0];
%vector_variation=0.5;



for loop_variable=vector_variation
fprintf('Sim for param: %f \n',loop_variable);

rng(740);       %random number generation seed






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CONSTANTS:
%Constants:
config.RAN_SLICING_AT_SPECTRUM_PLANNING=0;
config.RAN_SLICING_AT_PS=1;
config.RAN_SLICING_AT_AC=2;

%Slicing algorithm:
config.RAN_slicing_algorithm=config.RAN_SLICING_AT_AC;
%config.RAN_slicing_algorithm=config.RAN_SLICING_AT_PS;
%config.RAN_slicing_algorithm=config.RAN_SLICING_AT_SPECTRUM_PLANNING;



%AC algorithms:
config.NO_SLICING=0;  %The AC only accounts for the 1st check (global check).
config.SLICING_NO_DELTA=1;  %The AC does not account for the delta
config.SLICING_DELTA=2;     %The AC account for delta parameters

config.GBR=1;
config.NonGBR=0;

%config.AC_algorithm=config.NO_SLICING;
config.AC_algorithm=config.SLICING_NO_DELTA;  
%config.AC_algorithm=config.SLICING_DELTA;

%General scenario parameters
config.ISD=200;  %m  IntersiteDistance
config.cell_R=(config.ISD/2)*2/sqrt(3);   %Cell Radius
config.num_cells=1;
config.num_tenants=2;
config.time_step=1;     %Duration of the simulation time step in s.
config.simulation_duration=20000.0;    %Simulation duration in s.


config.activate_arrival_log=1;   %Set it to 1 to store the log with each arrival (very long output file with large number of cells!!!)
config.activate_list_UEs_log=1;  %Set it to 1 to store the log with all the UEs in all the time steps (very long output file).
                                 %If not activated, there are some
                                 %statistics (percentiles, stats per UE,
                                 %etc.) that will not be created.

config.num_RBs=275;   %Number of RBs (default for all the cells)
config.B_RB=30*12; %Bandwidth of one RB in kHz (12 subcarriers with Deltaf={15,30,60,120,240,480} kHz
config.Ptot=41; %dBm. Total power per LTE carrier (default for all the cells)
config.P_RB=config.Ptot-10.0*log10(config.num_RBs);  % dBm.  Power per RB.

config.antenna_gain=5.0;

config.noise_figure_UE=9;   %dB.
config.Pnoise_RB=-174+config.noise_figure_UE+10*log10(config.B_RB*1E3);   %Noise power per RB.
config.Pnoise_RB=power(10,0.1*config.Pnoise_RB);  %mW.  Noise power per RB.

%Propagation model parameters
config.prop_model_params.height_BS=10;   %meters
config.prop_model_params.height_UE=1.5;
config.prop_model_params.f=3.6;  %GHz
config.prop_model_params.d_BP=4*(config.prop_model_params.height_BS-1)*(config.prop_model_params.height_UE-1)*config.prop_model_params.f*1E9/3E8;
config.prop_model_params.dmin=10;     %minimum distance between a UE and a cell.
config.prop_model_params.sigma_LOS=3;     %dB
config.prop_model_params.sigma_NLOS=4;    %dB

%Spectral Efficiency computation parameters
config.spec_eff_params.SINRmin=-10.0;  % dB
config.spec_eff_params.SINRmin=power(10,0.1*config.spec_eff_params.SINRmin); %linear units
config.spec_eff_params.alfa=0.6;
config.spec_eff_params.Smax=8.8;  % b/s/Hz
config.spec_eff_params.SINRmax=power(2,config.spec_eff_params.Smax/config.spec_eff_params.alfa)-1;  %linear units


%Service profiles:
config.num_services=8;
config.services.name{1}='PREMIUM_Video';
config.services.type(1)=config.GBR;
config.services.priority(1)=40;
config.services.ARP(1)=2;
config.services.GBR_UL(1)=10000;
config.services.GBR_DL(1)=10000;
config.services.AMBR_UL(1)=0;
config.services.AMBR_DL(1)=0;
config.services.duration(1)=120; %seconds (average duration)
config.services.activity(1)=1; %Prob. that the service is active

config.services.name{2}='PREMIUM_Data';
config.services.type(2)=config.NonGBR;
config.services.priority(2)=60;
config.services.ARP(2)=2;
config.services.GBR_UL(2)=0;
config.services.GBR_DL(2)=0;
config.services.AMBR_UL(2)=300000;
config.services.AMBR_DL(2)=300000;
config.services.duration(2)=120; %seconds (average duration)
config.services.activity(2)=0.2; %Prob. that the service is active

config.services.name{3}='BASIC_Video';
config.services.type(3)=config.GBR;
config.services.priority(3)=40;
config.services.ARP(3)=3;
config.services.GBR_UL(3)=1000;
config.services.GBR_DL(3)=1000;
config.services.AMBR_UL(3)=0;
config.services.AMBR_DL(3)=0;
config.services.duration(3)=120; %seconds (average duration)
config.services.activity(3)=1; %Prob. that the service is active

config.services.name{4}='BASIC_Data';
config.services.type(4)=config.NonGBR;
config.services.priority(4)=80;
config.services.ARP(4)=3;
config.services.GBR_UL(4)=0;
config.services.GBR_DL(4)=0;
config.services.AMBR_UL(4)=100000;
config.services.AMBR_DL(4)=100000;
config.services.duration(4)=120; %seconds (average duration)
config.services.activity(4)=0.2; %Prob. that the service is active

config.services.name{5}='Corporative_BUS';
config.services.type(5)=config.NonGBR;
config.services.priority(5)=80;
config.services.ARP(5)=2;
config.services.GBR_UL(5)=0;
config.services.GBR_DL(5)=0;
config.services.AMBR_UL(5)=792000;
config.services.AMBR_DL(5)=792000;
config.services.duration(5)=120; %seconds (average duration)
config.services.activity(5)=0.8; %Prob. that the service is active

config.services.name{6}='MC_Video';
config.services.type(6)=config.GBR;
config.services.priority(6)=40;
config.services.ARP(6)=2;
config.services.GBR_UL(6)=2000;
config.services.GBR_DL(6)=2000;
config.services.AMBR_UL(6)=0;
config.services.AMBR_DL(6)=0;
config.services.duration(6)=120; %seconds (average duration)
config.services.activity(6)=1; %Prob. that the service is active

config.services.name{7}='MC_PTT';
config.services.type(7)=config.GBR;
config.services.priority(7)=7;
config.services.ARP(7)=1;
config.services.GBR_UL(7)=10;
config.services.GBR_DL(7)=10;
config.services.AMBR_UL(7)=0;
config.services.AMBR_DL(7)=0;
config.services.duration(7)=120; %seconds (average duration)
config.services.activity(7)=1; %Prob. that the service is active

config.services.name{8}='MC_Data';
config.services.type(8)=config.NonGBR;
config.services.priority(8)=55;
config.services.ARP(8)=3;
config.services.GBR_UL(8)=0;
config.services.GBR_DL(8)=0;
config.services.AMBR_UL(8)=100000;
config.services.AMBR_DL(8)=100000;
config.services.duration(8)=120; %seconds (average duration)
config.services.activity(8)=0.2; %Prob. that the service is active

%Parameter for the scheduling of NonGBR UEs
config.exponent_priority_nonGBR=1;  %Exponent of the resource sharing. (1/Priority)^exponent.

%Traffic parameters (default values: can be changed later on)
config.traffic_params.lambda=1.0;   %Sessions/sec. Session generation rate (default value)

%Admission parameters (default parameters: can be changed later on)
config.admission_params.alfa_th=1.0; %Admission threshold (default)
config.admission_params.beta=1.0; 
config.admission_params.gamma=1.0;
config.admission_params.Cextra_min=0.0;

config.time_window_utilisation_averaging=30.0;    %s Time window for averaging the RB utilisation
config.time_window_utilisation_averaging_samples=config.time_window_utilisation_averaging/config.time_step;

config.time_window_delta_averaging=300.0;    %s Time window for averaging the delta parameters of the AC.
config.time_window_delta_averaging_samples=config.time_window_delta_averaging/config.time_step;

config.time_window_bit_rate_averaging=30.0;     %s Time window for estimating the bit rate per RB
config.time_window_bit_rate_averaging_samples=config.time_window_bit_rate_averaging/config.time_step;


%traffic spatial parameters (gaussian distribution centered at (x,y) with
%deviation stdev (m)
%general params:
config.spatial_distrib_radius=200;  %The hotspot of each tenant will be distributed in a circle with radius 200m and with different angles.
config.spatial_distrib_stdev_general=300;
config.spatial_distrib_stdev_reference=300;  %This is for the reference computation of the total traffic with the hotspot at the center

%Tenant s: X=spatial_distrib_center(s,1); Y=spatial_distrib_center(s,2); spatial_distrib_stdev(s)
for s=1:config.num_tenants
    config.spatial_distrib_center(s,1)=config.spatial_distrib_radius*cos(2*pi*(s-1)/config.num_tenants);
    config.spatial_distrib_center(s,2)=config.spatial_distrib_radius*sin(2*pi*(s-1)/config.num_tenants);
    config.spatial_distrib_stdev(s)=config.spatial_distrib_stdev_general;
    config.spatial_distrib_stdev_ref(s)=config.spatial_distrib_stdev_reference;
end;    
    


%Capacity share parameters per tenant:
%config.C(1)=200/356.4;
config.C(1)=0.5;
%config.C(2)=50/356.4;
config.C(2)=0.2;
%config.C(3)=0.15;
%config.C(4)=0.15;




config.C_avg_multi_cell=zeros(config.num_tenants,1);
config.C_avg_multi_cell_samples=zeros(1+config.simulation_duration/config.time_step,config.num_tenants);
config.aggregate_avg_Rb_multi_cell=zeros(config.num_tenants,1);  %Measures the aggregate avg Rb of each tenant in the whole scenario

config.C_avg_multi_cell_GBR=zeros(config.num_tenants,1);
config.C_avg_multi_cell_GBR_samples=zeros(1+config.simulation_duration/config.time_step,config.num_tenants);
config.aggregate_avg_Rb_multi_cell_GBR=zeros(config.num_tenants,1);  %Measures the aggregate avg Rb of each tenant in the whole scenario

stats.num_adm_above_global_SAGBR=zeros(1,config.num_tenants);  %Measures the number of admissions above the global SAGBR of the whole scenario
stats.num_rej_below_global_SAGBR=zeros(1,config.num_tenants);  %Measures the number of rejections below the global SAGBR of the whole scenario

config.Cell_Capacity_theoretical=config.num_RBs*config.B_RB*config.spec_eff_params.Smax; %In kb/s
config.correction_factor=0.7757;   %Parameter Theta of the algorithm  (Factor equal to the ratio Effective Capacity/TotalCapacity, where Effective capacity reflects the maximum offered load for a max. blocking probability while TotalCapacity reflects the max capacity of a cell based on the amount of RBs).
config.Cell_Capacity_effective=config.Cell_Capacity_theoretical*config.correction_factor;     %Empirical correction (for accounting blocking 2%)


%SAGBR for the total scenario
config.SAGBR(1)=config.num_cells*config.C(1)*config.Cell_Capacity_effective;
config.SAGBR(2)=config.num_cells*config.C(2)*config.Cell_Capacity_effective;
%config.SAGBR(3)=config.num_cells*config.C(3)*config.Cell_Capacity_effective;
%config.SAGBR(4)=config.num_cells*config.C(4)*config.Cell_Capacity_effective;

%SAGBR "per cell"
config.cellSAGBR(1)=config.C(1)*config.Cell_Capacity_effective;
config.cellSAGBR(2)=config.C(2)*config.Cell_Capacity_effective;
%config.cellSAGBR(3)=config.C(3)*config.Cell_Capacity_effective;
%config.cellSAGBR(4)=config.C(4)*config.Cell_Capacity_effective;


%Distribute and initialize the cells:
if config.num_cells>19
    fprintf('ERROR (num_cells exceeds 19: not supported)!!!!!\n');
end;

for n=1:config.num_cells
    BS(n)=base;   
    BS(n).id=n;
    BS(n).init_BS(config);
end;

%If we want to have different traffic parameters per cell/tenant, specify
%them here (otherwise the traffic parameters are set to the default
%values).

%APPLY THE SPATIAL DISTRIBUTION TO DETERMINE THE TRAFFIC AT EACH CELL:
%Generation rates at the central point of the hot spot for each tenant:
lambda_ini(1)=0.27;
lambda_ini(2)=0.2;
%lambda_ini(3)=0.1;
%lambda_ini(4)=0.1;

ref_center_traffic=zeros(1,config.num_tenants);
ref_current_traffic=zeros(1,config.num_tenants);
for s=1:config.num_tenants
    ref_center_traffic(s)=0;   %Total generation rate for all the cells if the hotspot was centered at the origin
    ref_current_traffic(s)=0;  %Total generation rate for all the cells with the current hotspot position.
    for n=1:config.num_cells
        %Spatial distribution:
        d=sqrt(power(BS(n).pos(1)-config.spatial_distrib_center(s,1),2)+power(BS(n).pos(2)-config.spatial_distrib_center(s,2),2));
        dcenter_ref=sqrt(power(BS(n).pos(1),2)+power(BS(n).pos(2),2));
        
        BS(n).lambda_ini(s)=exp(-0.5*power(d/config.spatial_distrib_stdev(s),2));
        
        ref_center_traffic(s)=ref_center_traffic(s)+exp(-0.5*power(dcenter_ref/config.spatial_distrib_stdev_ref(s),2));
        ref_current_traffic(s)=ref_current_traffic(s)+BS(n).lambda_ini(s);
        
    end;
end;

%Now, let's determine the total generation rate by weighting the
%"lambda_ini" with the correction factor to ensure that the total traffic is the
%same as if the hotspot was centered at the origin
for s=1:config.num_tenants
    for n=1:config.num_cells
        %Spatial distribution:
        BS(n).lambda_ini(s)=BS(n).lambda_ini(s)*lambda_ini(s)*ref_center_traffic(s)/ref_current_traffic(s);
        
        %Time variation
        BS(n).lambda(s)=max(BS(n).lambda_ini(s)+BS(n).delta_lambda(s)*cos(2*pi*BS(n).freq_traffic_period(s)*(0-BS(n).time_shift(s))),1E-8); 
    end;
end;


%If we want to disable the traffic generation based on hotspots, specify
%the traffic per cell/tenant here:

BS(1).lambda_ini(1)=loop_variable;
BS(1).lambda_ini(2)=2.0;

BS(1).lambda(1)=max(BS(1).lambda_ini(1)+BS(1).delta_lambda(1)*cos(2*pi*BS(1).freq_traffic_period(1)*(0-BS(1).time_shift(1))),1E-8);
BS(1).lambda(2)=max(BS(1).lambda_ini(2)+BS(1).delta_lambda(2)*cos(2*pi*BS(1).freq_traffic_period(2)*(0-BS(1).time_shift(2))),1E-8);


%Specify the service mix for each cell and tenant(slice):
BS(1).service_mix(1,:)=[0.1,0.2,0.3,0.4,0,0,0,0]; %Tenant 1 (eMBB)
BS(1).service_mix(2,:)=[0,0,0,0,0,0.1,0.5,0.4]; %Tenant 2 (PS)



%If we want to have different frequencies per cell, specify them here.
%Otherwise, by default all the cells use the same frequency
BS(1).freq_index=1;
%BS(2).freq_index=1;
%BS(3).freq_index=3;

%PARAMETERS FOR THE SLICING AT PS AND SLICING AT SPECTRUM PLANNING
%Number of RBs per tenant in case of RAN slicing at PS (if not initialised
%here the default value is config.num_RB*C(s)) or at Spectrum Planning 
BS(1).num_RBs_per_tenant(1)=config.num_RBs*0.7;  %To be used in case of RAN slicing at PS or at Sp.Planning
BS(1).num_RBs_per_tenant(2)=config.num_RBs*0.3;
BS(1).alfa_th_GBR(1)=5/7;  %This factor is used to control the amount of RBs to be assigned to GBR traffic
BS(1).alfa_th_GBR(2)=2/3;



%After having modified cell-specific parameters, initialize the radio
%parameters of each cell and schedule the next session arrival rates:

for n=1:config.num_cells
    BS(n).init_radio_and_next_arrivals(config);
end;




%Main simulation
t_index=0;

%preallocate matrix for the logs:
if (config.activate_arrival_log==1) 
    arrival_log=zeros(1E6,600);
end;
 
num_entries=0;  %index to register the simulation_log

num_entries_UEs=0;  %index to register the list of UEs
%preallocate matrix for list of UEs and for the logs:
if config.activate_list_UEs_log==1
    list_UEs=zeros(1E8,11);
end;

for time=0:config.time_step:config.simulation_duration
    if mod(time,10)==0
        fprintf('Simulating time: %f \n',time);
    end;
    
    %if t_index==416
    %    fprintf('Ara \n',time);
    %end;
    
  
    changing_conditions=0;  %To identify if there is some change (arrival/end) in the time step.
    
    t_index=t_index+1;
    
    %Check session finalisations:
    for n=1:config.num_cells
        for s=1:config.num_tenants
            if BS(n).numUEs(s)>0
                end_process=0;
                i=1;
            else
                end_process=1;
            end;
            while ~end_process
                if BS(n).UElist{s}(i).end_session_time<=time
                    %Remove UE i from the list.
                    BS(n).UElist{s}=horzcat(BS(n).UElist{s}(1:i-1),BS(n).UElist{s}(i+1:BS(n).numUEs(s)));
                    BS(n).numUEs(s)=BS(n).numUEs(s)-1;   
                    
                    
                    changing_conditions=1;
                    
                    
                    %Note: the next UE to check is still the index i!!!
                    %(because we have shifted the UEs in the array)
                    %Then, we only increase i when we do not remove the UE.
                else
                    i=i+1;
                end;
                if i>BS(n).numUEs(s)
                    end_process=1;
                end;
            end;
        end;
    end;
    
    %Check session starts:
    for n=1:config.num_cells
        for s=1:config.num_tenants
            %Compute lambda
            BS(n).lambda(s)=max(BS(n).lambda_ini(s)+BS(n).delta_lambda(s)*cos(2*pi*BS(n).freq_traffic_period(s)*(time-BS(n).time_shift(s))),1E-8);
            
            
            %if BS(n).time_next_session_arrival(s)<=time
            while BS(n).time_next_session_arrival(s)<=time  %By putting the while, we allow multiple arrivals in a time step.
                %New session
                %First, execute the admission process (assuming it is independent of the
                %UE position):
                
                
                %Determine the service associated to this UE based on the
                %service_mix:
                aux_prob=rand();
                service=1;
                sum_aux=BS(n).service_mix(s,service);
                while (aux_prob>sum_aux)
                    service=service+1;
                    sum_aux=sum_aux+BS(n).service_mix(s,service);
                end;
                %Selected service is stored in variable "service".
                
                BS(n).num_session_attempts(s)=BS(n).num_session_attempts(s)+1;
                BS(n).num_session_attempts_per_service_per_tenant(s,service)=BS(n).num_session_attempts_per_service_per_tenant(s,service)+1;
                
                
                %Compute duration (even if it is not admitted later on)
                %duration_session=(-BS(n).duration(s))*log(1-rand());
                duration_session=(-config.services.duration(service))*log(1-rand());
                
                %BS(n).offered_load(s)=BS(n).offered_load(s)+duration_session*BS(n).Rbreq(s);
                BS(n).offered_load(s)=BS(n).offered_load(s)+duration_session*config.services.GBR_DL(service); %Offered load only accounts for GBR
                
                if config.services.type(service)==config.GBR
                    %Only execute admission for GBR services.
                    %admit=BS(n).admission(s,BS(n).Rbreq(s),config);
                    admit=BS(n).admission(s,config.services.GBR_DL(service),config,service);
                else
                    admit=1;
                end;
                
                if admit
                    %Generate the new UE of tenant s
                    BS(n).numUEs(s)=BS(n).numUEs(s)+1;
                    BS(n).UElist{s}(BS(n).numUEs(s))=UE;
                    BS(n).UElist{s}(BS(n).numUEs(s)).init_UE(config,BS,n,s);
                    BS(n).UElist{s}(BS(n).numUEs(s)).service=service;
                    BS(n).UElist{s}(BS(n).numUEs(s)).Rbreq=config.services.GBR_DL(service);
                    BS(n).UElist{s}(BS(n).numUEs(s)).activity_state=1; 
                    
                    %Compute the end session time for this UE:
                    BS(n).UElist{s}(BS(n).numUEs(s)).end_session_time=time+duration_session;
                                        
                    %Compute statistics:
                    %if (BS(n).avg_bit_rate_assigned_per_tenant(s)+BS(n).Rbreq(s))>config.cellSAGBR(s)
                    if (BS(n).avg_bit_rate_assigned_per_tenant(s)+config.services.GBR_DL(service))>config.cellSAGBR(s)
                        %Measure of SAGBR at cell level
                        BS(n).num_adm_above_SAGBR(s)=BS(n).num_adm_above_SAGBR(s)+1;
                    end;   
                    
                    %if (config.aggregate_avg_Rb_multi_cell(s)+BS(n).Rbreq(s))>config.SAGBR(s)
                    if (config.aggregate_avg_Rb_multi_cell(s)+config.services.GBR_DL(service))>config.SAGBR(s)
                        %Measure of SAGBR at the whole scenario
                        stats.num_adm_above_global_SAGBR(s)=stats.num_adm_above_global_SAGBR(s)+1;
                    end; 
                    
                    changing_conditions=1;
                
                else
                    %Count a blocking:
                    BS(n).num_blocks(s)=BS(n).num_blocks(s)+1;
                    BS(n).num_blocks_per_service_per_tenant(s,service)=BS(n).num_blocks_per_service_per_tenant(s,service)+1;
                     %Compute statistics: (for SAGBR we consider only the
                     %GBR bit rate
                    %if (BS(n).avg_bit_rate_assigned_per_tenant(s)+BS(n).Rbreq(s))<config.cellSAGBR(s)
                    if (sum(BS(n).avg_bit_rate_assigned_per_service_per_tenant(s,config.services.type(:)==config.GBR))+config.services.GBR_DL(service))<config.cellSAGBR(s)
                        %Measure of SAGBR at cell level
                        BS(n).num_rej_below_SAGBR(s)=BS(n).num_rej_below_SAGBR(s)+1;
                    end;
                    
                    %if (config.aggregate_avg_Rb_multi_cell(s)+BS(n).Rbreq(s))<config.SAGBR(s)
                    if (config.aggregate_avg_Rb_multi_cell_GBR(s)+config.services.GBR_DL(service))<config.SAGBR(s)
                        %Measure of SAGBR at the whole scenario
                        stats.num_rej_below_global_SAGBR(s)=stats.num_rej_below_global_SAGBR(s)+1;
                    end; 
                    
                end;
                %Schedule the arrival of next session for the tenant s:
                
                BS(n).time_next_session_arrival(s)=BS(n).time_next_session_arrival(s)+(-1/BS(n).lambda(s))*log(1-rand());
                %NOTE: We schedule not after "time" but after
                %"time_next_session_arrival" allowing multiple calls in a
                %time step.
                
                
                %Register the arrival and the system state
                
                if (config.activate_arrival_log==1) 
                %LOG FORMAT
                %[Time,Tenant,Cell,Duration,Rbreq,Admit,Bitrate per tenant
                %(inst), Bitrate per tenant(avg), NumRBper tenant (inst),
                %NumRB per tenant (avg), DeltaC, DeltaCext, DeltaCbal, CongestionStatus]
                num_entries=num_entries+1;
                arrival_log(num_entries,1:7)=[time,s,n,duration_session,service,config.services.GBR_DL(service),admit];
                index_log=8;
                for naux=1:config.num_cells
                    arrival_log(num_entries,index_log:(index_log+config.num_tenants-1))=BS(n).bit_rate_assigned_per_tenant;
                    arrival_log(num_entries,(index_log+config.num_tenants):(index_log+2*config.num_tenants-1))=BS(n).avg_bit_rate_assigned_per_tenant;
                    arrival_log(num_entries,(index_log+2*config.num_tenants):(index_log+3*config.num_tenants-1))=BS(n).num_assigned_RB_per_tenant;
                    arrival_log(num_entries,(index_log+3*config.num_tenants):(index_log+4*config.num_tenants-1))=BS(n).avg_num_RB_per_tenant;
                    arrival_log(num_entries,(index_log+4*config.num_tenants):(index_log+5*config.num_tenants-1))=BS(n).DeltaC;
                    arrival_log(num_entries,(index_log+5*config.num_tenants):(index_log+6*config.num_tenants-1))=BS(n).DeltaCext;
                    arrival_log(num_entries,(index_log+6*config.num_tenants):(index_log+7*config.num_tenants-1))=BS(n).DeltaCbal;
                    arrival_log(num_entries,(index_log+7*config.num_tenants):(index_log+7*config.num_tenants))=BS(n).congestion_status;
                    index_log=index_log+7*config.num_tenants+1;
                end;
                end;
                
            end;
        end;
    end;
    
    %Check activity of the UEs:
    for n=1:config.num_cells
        for s=1:config.num_tenants
            for i=1:BS(n).numUEs(s)
                previous_activity=BS(n).UElist{s}(i).activity_state;
                if rand()<config.services.activity(BS(n).UElist{s}(i).service)
                    BS(n).UElist{s}(i).activity_state=1;
                else
                    BS(n).UElist{s}(i).activity_state=0;
                end;
                if previous_activity ~= BS(n).UElist{s}(i).activity_state 
                    changing_conditions=1;
                end;
            end;    
        end;
    end;
    
    %Assess performance for the current time step.
    %ONLY IF SOME CHANGE HAS OCCURRED (NEW UEs OR UES ending).
    if changing_conditions
        compute_occupation(config,BS);
    end;
    
    %1) Update the averages of 
    %BS(n).avg_num_RB_per_tenant   %Average number of RBs used by each tenant.
    %BS(n).avg_RB_utilisation_per_tenant    %Average real capacity share (RB_utilisation) per tenant.
    %BS(n).Rb_estimate_per_RB             %Estimate of bit rate per RB (in kb/s) achieved in the cell
    %BS(n).avg_bit_rate_assigned_per_tenant
    
    for n=1:config.num_cells
       BS(n).total_assigned_RB_sample(t_index)=sum(BS(n).num_assigned_RB_per_tenant);
       BS(n).total_assigned_RB_per_service_sample(t_index,:)=sum(BS(n).num_assigned_RB_per_service_per_tenant);
       
        
       for s=1:config.num_tenants
           BS(n).num_assigned_RB_per_tenant_sample(t_index,s)=BS(n).num_assigned_RB_per_tenant(s);
           BS(n).num_assigned_RB_per_service_per_tenant_sample(t_index,s,:)=BS(n).num_assigned_RB_per_service_per_tenant(s,:);
           
           BS(n).bit_rate_assigned_per_tenant_sample(t_index,s)=BS(n).bit_rate_assigned_per_tenant(s);
           BS(n).bit_rate_assigned_per_service_per_tenant_sample(t_index,s,:)=BS(n).bit_rate_assigned_per_service_per_tenant(s,:);
            
           if t_index<=config.time_window_utilisation_averaging_samples
               BS(n).avg_num_RB_per_tenant(s)=mean(BS(n).num_assigned_RB_per_tenant_sample(1:t_index,s));
               BS(n).avg_bit_rate_assigned_per_tenant(s)=mean(BS(n).bit_rate_assigned_per_tenant_sample(1:t_index,s));
               
               BS(n).avg_num_RB_per_service_per_tenant(s,:)=mean(BS(n).num_assigned_RB_per_service_per_tenant_sample(1:t_index,s,:),1);
               BS(n).avg_bit_rate_assigned_per_service_per_tenant(s,:)=mean(BS(n).bit_rate_assigned_per_service_per_tenant_sample(1:t_index,s,:),1);
           else
               BS(n).avg_num_RB_per_tenant(s)=mean(BS(n).num_assigned_RB_per_tenant_sample((t_index-config.time_window_utilisation_averaging_samples):t_index,s));
               BS(n).avg_bit_rate_assigned_per_tenant(s)=mean(BS(n).bit_rate_assigned_per_tenant_sample((t_index-config.time_window_utilisation_averaging_samples):t_index,s));
               
               BS(n).avg_num_RB_per_service_per_tenant(s,:)=mean(BS(n).num_assigned_RB_per_service_per_tenant_sample((t_index-config.time_window_utilisation_averaging_samples):t_index,s,:),1);
               BS(n).avg_bit_rate_assigned_per_service_per_tenant(s,:)=mean(BS(n).bit_rate_assigned_per_service_per_tenant_sample((t_index-config.time_window_utilisation_averaging_samples):t_index,s,:),1);
           end; 
           BS(n).avg_num_RB_per_tenant_sample(t_index,s)=BS(n).avg_num_RB_per_tenant(s);
           BS(n).avg_RB_utilisation_per_tenant(s)=BS(n).avg_num_RB_per_tenant(s)/BS(n).num_RBs;
           BS(n).avg_RB_utilisation_per_tenant_sample(t_index,s)=BS(n).avg_RB_utilisation_per_tenant(s);
           
           BS(n).avg_bit_rate_assigned_per_tenant_sample(t_index,s)=BS(n).avg_bit_rate_assigned_per_tenant(s);
           
           BS(n).avg_num_RB_per_service_per_tenant_sample(t_index,s,:)=BS(n).avg_num_RB_per_service_per_tenant(s,:);
           BS(n).avg_RB_utilisation_per_service_per_tenant(s,:)=BS(n).avg_num_RB_per_service_per_tenant(s,:)/BS(n).num_RBs;
           BS(n).avg_RB_utilisation_per_service_per_tenant_sample(t_index,s,:)=BS(n).avg_RB_utilisation_per_service_per_tenant(s,:);
           
           BS(n).avg_bit_rate_assigned_per_service_per_tenant_sample(t_index,s,:)=BS(n).avg_bit_rate_assigned_per_service_per_tenant(s,:);
           
           
           %Compute data volume
           for i=1:BS(n).numUEs(s)
               BS(n).data_volume_per_tenant(s)=BS(n).data_volume_per_tenant(s)+BS(n).UElist{s}(i).assigned_bit_rate*config.time_step;
               serv=BS(n).UElist{s}(i).service;
               
               %Statistics at UE level:
               BS(n).sum_SINR(s)=BS(n).sum_SINR(s)+10*log10(BS(n).UElist{s}(i).SINR); %Average of the SINR in dB
               BS(n).sum_required_RBs(s)=BS(n).sum_required_RBs(s)+BS(n).UElist{s}(i).required_RBs;
               BS(n).sum_required_RBs_per_service_per_tenant(s,serv)=BS(n).sum_required_RBs_per_service_per_tenant(s,serv)+BS(n).UElist{s}(i).required_RBs;
               if (BS(n).UElist{s}(i).assigned_RBs<BS(n).UElist{s}(i).required_RBs)
                   BS(n).num_degraded_UEs(s)=BS(n).num_degraded_UEs(s)+1;
                   BS(n).num_degraded_UEs_per_service_per_tenant(s,serv)=BS(n).num_degraded_UEs_per_service_per_tenant(s,serv)+1;
               end;
               BS(n).sum_Rb_degradation(s)= BS(n).sum_Rb_degradation(s)+100*(BS(n).UElist{s}(i).Rbreq-BS(n).UElist{s}(i).assigned_bit_rate)/BS(n).UElist{s}(i).Rbreq;
               BS(n).num_samples_UEs(s)=BS(n).num_samples_UEs(s)+1; 
               BS(n).sum_Rb_degradation_per_service_per_tenant(s,serv)=BS(n).sum_Rb_degradation_per_service_per_tenant(s,serv)+100*(BS(n).UElist{s}(i).Rbreq-BS(n).UElist{s}(i).assigned_bit_rate)/BS(n).UElist{s}(i).Rbreq;
               BS(n).num_samples_UEs_per_service_per_tenant(s,serv)=BS(n).num_samples_UEs_per_service_per_tenant(s,serv)+1; 
               
               
               
               if config.activate_list_UEs_log==1
                    num_entries_UEs=num_entries_UEs+1;
                    list_UEs(num_entries_UEs,1)=time;
                    list_UEs(num_entries_UEs,2)=n;
                    list_UEs(num_entries_UEs,3)=s;
                    list_UEs(num_entries_UEs,4)=BS(n).UElist{s}(i).distance;
                    list_UEs(num_entries_UEs,5)=10*log10(BS(n).UElist{s}(i).SINR);
                    list_UEs(num_entries_UEs,6)=BS(n).UElist{s}(i).required_RBs;
                    list_UEs(num_entries_UEs,7)=BS(n).UElist{s}(i).assigned_RBs;
                    if (BS(n).UElist{s}(i).assigned_RBs<BS(n).UElist{s}(i).required_RBs)
                        list_UEs(num_entries_UEs,8)=1;
                    else
                        list_UEs(num_entries_UEs,8)=0;
                    end;
                    list_UEs(num_entries_UEs,9)=BS(n).UElist{s}(i).assigned_bit_rate;
                    list_UEs(num_entries_UEs,10)=BS(n).UElist{s}(i).service;
                    list_UEs(num_entries_UEs,11)=BS(n).UElist{s}(i).activity_state;
               end;
           end;
           
           %Check congestion status (per tenant)
           if BS(n).congestion_status_per_tenant(s)
                BS(n).num_congested_samples_per_tenant(s)=BS(n).num_congested_samples_per_tenant(s)+1;
           end;
       end;
       
       
       %Estimate of bit rate per RB:
       %Note: It accounts for the total bit rate (for all services)
       BS(n).total_assigned_bit_rate_sample(t_index)=sum(BS(n).bit_rate_assigned_per_tenant);
       BS(n).total_assigned_bit_rate_per_service_sample(t_index,:)=sum(BS(n).bit_rate_assigned_per_service_per_tenant);
       %BS(n).bit_rate_assigned_per_tenant_sample(t_index,:)=BS(n).bit_rate_assigned_per_tenant(:);
       
       if t_index<=config.time_window_bit_rate_averaging_samples
           aggregated_Rb=sum(BS(n).total_assigned_bit_rate_sample(1:t_index));
           aggregated_RBs=sum(BS(n).total_assigned_RB_sample(1:t_index));
       else
           aggregated_Rb=sum(BS(n).total_assigned_bit_rate_sample((t_index-config.time_window_bit_rate_averaging_samples):t_index));
           aggregated_RBs=sum(BS(n).total_assigned_RB_sample((t_index-config.time_window_bit_rate_averaging_samples):t_index));
       end;
       
       if aggregated_Rb>0
           %Update the new average of Bit rate
           BS(n).Rb_estimate_per_RB=aggregated_Rb/aggregated_RBs;
       end;
       %Note: if aggregated_Rb=0, meaning that no bit rate has been
       %obtained in the last period, the value of Rb_estimate_per_RB is
       %kept to the initial value.
       BS(n).Rb_estimate_per_RB_sample(t_index)=BS(n).Rb_estimate_per_RB;
       
       %Check congestion status:
       if BS(n).congestion_status
           BS(n).num_congested_samples=BS(n).num_congested_samples+1;
       end;
       
       
    end;
    
    
    %2) Compute the average RB utilisation of each tenant at multi-cell level
    %and the aggregate bit rate of each tenant at multi-cell level
    for s=1:config.num_tenants
        config.C_avg_multi_cell(s)=0;
        config.aggregate_avg_Rb_multi_cell(s)=0;
        
        config.C_avg_multi_cell_GBR(s)=0;
        config.aggregate_avg_Rb_multi_cell_GBR(s)=0;
        
        for n=1:config.num_cells
            config.C_avg_multi_cell(s)=config.C_avg_multi_cell(s)+BS(n).avg_RB_utilisation_per_tenant(s);
            config.aggregate_avg_Rb_multi_cell(s)=config.aggregate_avg_Rb_multi_cell(s)+BS(n).avg_bit_rate_assigned_per_tenant(s);
            
            config.C_avg_multi_cell_GBR(s)=config.C_avg_multi_cell_GBR(s)+sum(BS(n).avg_RB_utilisation_per_service_per_tenant(s,config.services.type(:)==config.GBR));
            config.aggregate_avg_Rb_multi_cell_GBR(s)=config.aggregate_avg_Rb_multi_cell_GBR(s)+sum(BS(n).avg_bit_rate_assigned_per_service_per_tenant(s,config.services.type(:)==config.GBR));
            
        end;
        config.C_avg_multi_cell(s)=config.C_avg_multi_cell(s)/config.num_cells;
        config.C_avg_multi_cell_samples(t_index,s)=config.C_avg_multi_cell(s);
        
        config.C_avg_multi_cell_GBR(s)=config.C_avg_multi_cell_GBR(s)/config.num_cells;
        config.C_avg_multi_cell_GBR_samples(t_index,s)=config.C_avg_multi_cell_GBR(s);
        
    end;

    
    %3) Compute and update the delta_C parameters
    %Computation of deltaCext:
    for n=1:config.num_cells
       for s=1:config.num_tenants
            BS(n).DeltaCext_sample(t_index,s)=0;
            for saux=1:config.num_tenants
               if  saux~=s
                    %BS(n).DeltaCext_sample(t_index,s)=BS(n).DeltaCext_sample(t_index,s)+config.C(saux)*config.correction_factor-BS(n).avg_RB_utilisation_per_tenant(saux);
                    
                    %We only account for the avg RB utilisation of GBR
                    %services!!
                    BS(n).DeltaCext_sample(t_index,s)=BS(n).DeltaCext_sample(t_index,s)+config.C(saux)*config.correction_factor-sum(BS(n).avg_RB_utilisation_per_service_per_tenant(saux,config.services.type(:)==config.GBR));
                    
               end; 
            end;
            BS(n).DeltaCext_sample(t_index,s)=max(BS(n).DeltaCext_sample(t_index,s),0);
            
            %Average:
            if t_index<=config.time_window_delta_averaging_samples
               BS(n).DeltaCext(s)=mean(BS(n).DeltaCext_sample(1:t_index,s));
            else
               BS(n).DeltaCext(s)=mean(BS(n).DeltaCext_sample((t_index-config.time_window_delta_averaging_samples):t_index,s));
            end;
            
            BS(n).DeltaCext_avg_sample(t_index,s)=BS(n).DeltaCext(s);
            
       end;
    end;
    
    %Computation of deltaCbal and DeltaC total:
    for s=1:config.num_tenants
        for n=1:config.num_cells
            BS(n).DeltaCbal_sample(t_index,s)=(config.num_cells-1)*config.C(s); 
            for naux=1:config.num_cells
                if naux~=n
                    %Only account for GBR services:
                    %BS(n).DeltaCbal_sample(t_index,s)=BS(n).DeltaCbal_sample(t_index,s)-BS(naux).avg_RB_utilisation_per_tenant(s);
                    BS(n).DeltaCbal_sample(t_index,s)=BS(n).DeltaCbal_sample(t_index,s)-sum(BS(naux).avg_RB_utilisation_per_service_per_tenant(s,config.services.type(:)==config.GBR));
                end;
            end;
            
            %Average:
            if t_index<=config.time_window_delta_averaging_samples
               BS(n).DeltaCbal(s)=mean(BS(n).DeltaCbal_sample(1:t_index,s));
            else
               BS(n).DeltaCbal(s)=mean(BS(n).DeltaCbal_sample((t_index-config.time_window_delta_averaging_samples):t_index,s));
            end;
            BS(n).DeltaCbal_avg_sample(t_index,s)=BS(n).DeltaCbal(s);
            
            %Computation of DeltaC total:
            
            if BS(n).DeltaCext_sample(t_index,s)>BS(n).Cextra_min
                BS(n).DeltaC_sample(t_index,s)=BS(n).beta*BS(n).DeltaCext_sample(t_index,s);
            else
                BS(n).DeltaC_sample(t_index,s)=BS(n).gamma*BS(n).DeltaCbal_sample(t_index,s);
            end;
            
            %Average:
            if t_index<=config.time_window_delta_averaging_samples
               BS(n).DeltaC(s)=mean(BS(n).DeltaC_sample(1:t_index,s));
            else
               BS(n).DeltaC(s)=mean(BS(n).DeltaC_sample((t_index-config.time_window_delta_averaging_samples):t_index,s));
            end;
            BS(n).DeltaC_avg_sample(t_index,s)=BS(n).DeltaC(s);
            
        end;
    end;
    
   
end;

%Reduce the size of the variable list_UEs to include only the actual
%entries
if config.activate_list_UEs_log==1
    list_UEs=list_UEs(1:num_entries_UEs,:);
end;

% Measure final statistics and plots.



%Average along the whole simulation:    
stats.avg_Cshare_multicell=mean(config.C_avg_multi_cell_samples);
stats.avg_Cshare_multicell_GBR=mean(config.C_avg_multi_cell_GBR_samples);

stats.num_blocks_total=zeros(1,config.num_tenants);
stats.num_session_attempts_total=zeros(1,config.num_tenants);

stats.num_blocks_total_per_service_per_tenant=zeros(config.num_tenants,config.num_services);
stats.num_session_attempts_total_per_service_per_tenant=zeros(config.num_tenants,config.num_services);

stats.offered_load_total=zeros(1,config.num_tenants);
stats.congestion_prob_per_cell=zeros(config.num_cells,1);

stats.num_adm_above_SAGBR_total=zeros(1,config.num_tenants); %For the cellSAGBR
stats.num_rej_below_SAGBR_total=zeros(1,config.num_tenants); %For the cellSAGBR 


for n=1:config.num_cells
    stats.avg_Cshare_per_cell(n,:)=mean(BS(n).avg_RB_utilisation_per_tenant_sample);
    stats.avg_Cshare_per_cell_GBR(n,:)=mean(sum(BS(n).avg_RB_utilisation_per_service_per_tenant_sample(:,:,config.services.type(:)==config.GBR),3));
    
    stats.avg_Cbal(n,:)=mean(BS(n).DeltaCbal_avg_sample);
    stats.avg_Cext(n,:)=mean(BS(n).DeltaCext_avg_sample);
    stats.avg_Ctot(n,:)=mean(BS(n).DeltaC_avg_sample);
    stats.avg_bit_rate_per_tenant_per_cell(n,:)=mean(BS(n).bit_rate_assigned_per_tenant_sample);
    stats.avg_bit_rate_per_tenant_per_cell_GBR(n,:)=mean(sum(BS(n).bit_rate_assigned_per_service_per_tenant_sample(:,:,config.services.type(:)==config.GBR),3));
    
    stats.avg_bit_rate_per_service_per_tenant_per_cell(n,:,:)=mean(BS(n).bit_rate_assigned_per_service_per_tenant_sample);
    
    
    for s=1:config.num_tenants
        stats.blocking_prob_per_cell(n,s)=BS(n).num_blocks(s)/BS(n).num_session_attempts(s);
        stats.num_blocks_total(s)=stats.num_blocks_total(s)+BS(n).num_blocks(s);
        stats.num_session_attempts_total(s)=stats.num_session_attempts_total(s)+BS(n).num_session_attempts(s);
        stats.offered_load_total(s)=stats.offered_load_total(s)+BS(n).offered_load(s);  %It's only GBR load.
        
        stats.blocking_prob_per_cell_per_service(n,s,:)=BS(n).num_blocks_per_service_per_tenant(s,:)./BS(n).num_session_attempts_per_service_per_tenant(s,:);
        stats.num_blocks_total_per_service_per_tenant(s,:)=stats.num_blocks_total_per_service_per_tenant(s,:)+BS(n).num_blocks_per_service_per_tenant(s,:);
        stats.num_session_attempts_total_per_service_per_tenant(s,:)=stats.num_session_attempts_total_per_service_per_tenant(s,:)+BS(n).num_session_attempts_per_service_per_tenant(s,:);
        
        stats.adm_prob_above_SAGBR_per_cell(n,s)=BS(n).num_adm_above_SAGBR(s)/BS(n).num_session_attempts(s);
        stats.rej_prob_below_SAGBR_per_cell(n,s)=BS(n).num_rej_below_SAGBR(s)/BS(n).num_session_attempts(s);
        
        stats.num_adm_above_SAGBR_total(s)=stats.num_adm_above_SAGBR_total(s)+BS(n).num_adm_above_SAGBR(s);
        stats.num_rej_below_SAGBR_total(s)=stats.num_rej_below_SAGBR_total(s)+BS(n).num_rej_below_SAGBR(s);
        
        stats.data_volume_per_cell(n,s)=BS(n).data_volume_per_tenant(s)/(8*1E6);  %Measured in GByte
        
        stats.congestion_prob_per_cell_per_tenant(n,s)=BS(n).num_congested_samples_per_tenant(s)/t_index;
        
        stats.avg_SINR(n,s)=BS(n).sum_SINR(s)/BS(n).num_samples_UEs(s);
        stats.avg_required_RBs(n,s)=BS(n).sum_required_RBs(s)/BS(n).num_samples_UEs(s);
        stats.prob_degraded_UE(n,s)=BS(n).num_degraded_UEs(s)/BS(n).num_samples_UEs(s);
        stats.avg_perc_Rb_degradation(n,s)=BS(n).sum_Rb_degradation(s)/BS(n).num_samples_UEs(s);
        
        stats.avg_required_RBs_per_service_per_tenant(n,s,:)=BS(n).sum_required_RBs_per_service_per_tenant(s,:)./BS(n).num_samples_UEs_per_service_per_tenant(s,:);
        stats.prob_degraded_UE_per_service_per_tenant(n,s,:)=BS(n).num_degraded_UEs_per_service_per_tenant(s,:)./BS(n).num_samples_UEs_per_service_per_tenant(s,:);
        stats.avg_perc_Rb_degradation_per_service_per_tenant(n,s,:)=BS(n).sum_Rb_degradation_per_service_per_tenant(s,:)./BS(n).num_samples_UEs_per_service_per_tenant(s,:);
        
        
    end;
    stats.avg_session_rate_per_cell(n,:)=BS(n).num_session_attempts/config.simulation_duration;
    stats.avg_session_rate_per_cell_per_service_per_tenant(n,:,:)=BS(n).num_session_attempts_per_service_per_tenant/config.simulation_duration;
    stats.offered_load_per_cell(n,:)=BS(n).offered_load/config.simulation_duration; %Measured in kb/s
    stats.congestion_prob_per_cell(n)=BS(n).num_congested_samples/t_index;
    stats.avg_RB_occupation(n)=mean(BS(n).total_assigned_RB_sample);
    stats.avg_RB_occupation_GBR(n)=mean(sum(BS(n).total_assigned_RB_per_service_sample(:,config.services.type(:)==config.GBR),2));
    
    stats.avg_RB_occupation_per_service(n,:)=mean(BS(n).total_assigned_RB_per_service_sample);
    stats.avg_RB_occupation_per_service_per_tenant(n,:,:)=mean(BS(n).num_assigned_RB_per_service_per_tenant_sample);
    
    stats.avg_bit_rate_per_RB(n)=mean(BS(n).Rb_estimate_per_RB_sample);
    stats.avg_sp_efficiency(n)=stats.avg_bit_rate_per_RB(n)/config.B_RB;
end;
stats.blocking_prob_total=stats.num_blocks_total./stats.num_session_attempts_total;
stats.avg_session_rate_total=stats.num_session_attempts_total/config.simulation_duration;

stats.blocking_prob_total_per_service_per_tenant=stats.num_blocks_total_per_service_per_tenant./stats.num_session_attempts_total_per_service_per_tenant;
stats.avg_session_rate_total_per_service_per_tenant=stats.num_session_attempts_total_per_service_per_tenant/config.simulation_duration;

stats.offered_load_total=stats.offered_load_total/config.simulation_duration;

stats.avg_bit_rate_per_tenant_total=sum(stats.avg_bit_rate_per_tenant_per_cell,1);
stats.avg_bit_rate_per_tenant_total_GBR=sum(stats.avg_bit_rate_per_tenant_per_cell_GBR,1);

stats.avg_bit_rate_per_service_per_tenant_total(:,:)=sum(stats.avg_bit_rate_per_service_per_tenant_per_cell,1);


stats.data_volume_per_tenant=sum(stats.data_volume_per_cell,1);

stats.adm_prob_above_SAGBR_total=stats.num_adm_above_SAGBR_total./stats.num_session_attempts_total; %For the cellSAGBR
stats.rej_prob_below_SAGBR_total=stats.num_rej_below_SAGBR_total./stats.num_session_attempts_total; %For the cellSAGBR

stats.adm_prob_above_GLOBAL_SAGBR=stats.num_adm_above_global_SAGBR./stats.num_session_attempts_total; %For the global SAGBR
stats.rej_prob_below_GLOBAL_SAGBR=stats.num_rej_below_global_SAGBR./stats.num_session_attempts_total; %For the global SAGBR



stats.avg_required_RBs_total_per_service_per_tenant=zeros(config.num_tenants,config.num_services);
stats.prob_degraded_UE_total_per_service_per_tenant=zeros(config.num_tenants,config.num_services);
stats.avg_perc_Rb_degradation_total_per_service_per_tenant=zeros(config.num_tenants,config.num_services);
for s=1:config.num_tenants
    stats.avg_SINR_total(s)=0;
    stats.avg_required_RBs_total(s)=0;
    stats.prob_degraded_UE_total(s)=0;
    stats.avg_perc_Rb_degradation_total(s)=0;
    num_samples_aux=0;
    num_samples_per_service_aux=zeros(1,config.num_services);
    
    for n=1:config.num_cells
        stats.avg_SINR_total(s)=stats.avg_SINR_total(s)+BS(n).sum_SINR(s);
        stats.avg_required_RBs_total(s)=stats.avg_required_RBs_total(s)+BS(n).sum_required_RBs(s);
        stats.prob_degraded_UE_total(s)=stats.prob_degraded_UE_total(s)+BS(n).num_degraded_UEs(s);
        stats.avg_perc_Rb_degradation_total(s)=stats.avg_perc_Rb_degradation_total(s)+BS(n).sum_Rb_degradation(s);
        num_samples_aux=num_samples_aux+BS(n).num_samples_UEs(s);
        
        for serv=1:config.num_services
            stats.avg_required_RBs_total_per_service_per_tenant(s,serv)=stats.avg_required_RBs_total_per_service_per_tenant(s,serv)+BS(n).sum_required_RBs_per_service_per_tenant(s,serv);
            stats.prob_degraded_UE_total_per_service_per_tenant(s,serv)=stats.prob_degraded_UE_total_per_service_per_tenant(s,serv)+BS(n).num_degraded_UEs_per_service_per_tenant(s,serv);
            stats.avg_perc_Rb_degradation_total_per_service_per_tenant(s,serv)=stats.avg_perc_Rb_degradation_total_per_service_per_tenant(s,serv)+BS(n).sum_Rb_degradation_per_service_per_tenant(s,serv);
            num_samples_per_service_aux(serv)=num_samples_per_service_aux(serv)+BS(n).num_samples_UEs_per_service_per_tenant(s,serv);
        end;
        
    end;
    stats.avg_SINR_total(s)=stats.avg_SINR_total(s)/num_samples_aux;
    stats.avg_required_RBs_total(s)=stats.avg_required_RBs_total(s)/num_samples_aux;
    stats.prob_degraded_UE_total(s)=stats.prob_degraded_UE_total(s)/num_samples_aux;
    stats.avg_perc_Rb_degradation_total(s)=stats.avg_perc_Rb_degradation_total(s)/num_samples_aux;
    
    stats.avg_required_RBs_total_per_service_per_tenant(s,:)=stats.avg_required_RBs_total_per_service_per_tenant(s,:)./num_samples_per_service_aux(1,:);
    stats.prob_degraded_UE_total_per_service_per_tenant(s,:)=stats.prob_degraded_UE_total_per_service_per_tenant(s,:)./num_samples_per_service_aux(1,:);
    stats.avg_perc_Rb_degradation_total_per_service_per_tenant(s,:)=stats.avg_perc_Rb_degradation_total_per_service_per_tenant(s,:)./num_samples_per_service_aux(1,:);
    
    
    if config.activate_list_UEs_log==1
        aux_matrix=list_UEs(find(list_UEs(:,3)==s),5);   %SINR of the UEs of tenant s (in any cell)
        stats.perc5_SINR_total(s)=prctile(aux_matrix,5);
        stats.perc95_SINR_total(s)=prctile(aux_matrix,95);
        
        aux_matrix=list_UEs(find(list_UEs(:,3)==s),9);   %Assigned Rb of the UEs of tenant s (in any cell)
        stats.avg_Rb_per_user_total(s)=mean(aux_matrix);
        stats.perc5_Rb_per_user_total(s)=prctile(aux_matrix,5);
        stats.perc95_Rb_per_user_total(s)=prctile(aux_matrix,95);
    
        for serv=1:config.num_services
            for n=1:config.num_cells
                aux_matrix=list_UEs(find((list_UEs(:,3)==s) & (list_UEs(:,10)==serv) & (list_UEs(:,2)==n) & (list_UEs(:,11)==1)),9);
                stats.avg_Rb_per_user_total_per_cell_per_service_per_tenant(n,s,serv)=mean(aux_matrix);
                stats.perc5_Rb_per_user_total_per_cell_per_service_per_tenant(n,s,serv)=prctile(aux_matrix,5);
                stats.perc95_Rb_per_user_total_per_cell_per_service_per_tenant(n,s,serv)=prctile(aux_matrix,95);
            end;
            aux_matrix=list_UEs(find((list_UEs(:,3)==s) & (list_UEs(:,10)==serv) & (list_UEs(:,11)==1)),9);
            stats.avg_Rb_per_user_total_per_service_per_tenant(s,serv)=mean(aux_matrix);
            stats.perc5_Rb_per_user_total_per_service_per_tenant(s,serv)=prctile(aux_matrix,5);
            stats.perc95_Rb_per_user_total_per_service_per_tenant(s,serv)=prctile(aux_matrix,95);
        end;
    end;
end;


   
%generate_plots(config, BS);
   

if (config.activate_arrival_log==1) 
    arrival_log=arrival_log(1:num_entries,1:index_log-1); %To reduce the size of the arrival_log
end;

name_output_file=['simAC_50_20_T1_',num2str(loop_variable),'_T2_2.mat'];
save(name_output_file,'-v7.3');  %The option '-v7.3' is needed to store very large variables such as list_UEs
%save('prova1.mat');

end;
