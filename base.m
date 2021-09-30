classdef base < handle
    % Class that represents an LTE cell. 
    properties
        id   %eNodeB ID. It is also the index in the eNodeB site array in the simulator
        pos  %Position in meters (x,y)
        lambda   %Total RAB(or session) generation rate. Vector wher every component corresponds to a tenant
        %duration  %Average duration of each session. Vector wher every component corresponds to a tenant
        %Rbreq     %Required bit rate of each session. Vector wher every component corresponds to a tenant
        UElist     %List of UEs. Each component is a cell. The cell contains the UEs of each tenant.
        numUEs
    
        num_RBs
        Ptot
        antenna_gain
        cell_R
        PIRE_RB
        P_RB
        
        freq_index                 %Number of the frequency channel
        
        alfa_th
        alfa_th_GBR                %Parameter in the AC applied to GBR resources with RAN slicing at PS or at PS (is to consider only a fraction of resources for GBR, leaving the rest for Non-GBR).
        beta
        gamma
        
        C       %Nominal capacity share per tenant
        
        DeltaC                        %Current value of DeltaC per tenant (average over the delta window)
        DeltaC_sample                 %Matrix with the instantaneous values of DeltaC per tenant and per time step
        DeltaC_avg_sample             %Matrix with the averaged values of DeltaCext per tenant and per time step
        
        DeltaCext                     %Current value of DeltaCext per tenant (average over the delta window)
        DeltaCext_sample              %Matrix with the instantaneous values of DeltaCext per tenant and per time step  
        DeltaCext_avg_sample          %Matrix with the averaged values of DeltaCext per tenant and per time step 
        
        Cextra_min                    %Minimum value to switch between DeltaCext and DeltaCbal
        
        DeltaCbal                     %Current value of DeltaCbal per tenant (average over the delta window)
        DeltaCbal_sample              %Matrix with the instantaneous values of DeltaCbal per tenant and per time step
        DeltaCbal_avg_sample          %Matrix with the averaged values of DeltaCbal per tenant and per time step 
        
        num_assigned_RB_per_tenant    %Number of assigned RB to each tenant in a time step
        num_assigned_RB_per_tenant_sample    %Matrix where each column is the samples of assigned RB per tenant in each time step
        total_assigned_RB_sample      %Vector with the samples of the total assigned RBs in each time step
        
        bit_rate_assigned_per_tenant  %Bit rate in the assigned RBs to each tenant in a time step
        bit_rate_assigned_per_tenant_sample %Matrix with the samples of the total assigned bit rate to each tenant in each time step
        total_assigned_bit_rate_sample    %Vector with the samples of the total assigned bit rate in each time step
        
        avg_bit_rate_assigned_per_tenant    %Average bit rate in the assigned RBs to each tenant
        avg_bit_rate_assigned_per_tenant_sample  %Matrix where each column is the average of assigned bit rate per tenant in each time step
        
        
      
        congestion_status             %Takes value 1 if the BS is congested in a time step.
        num_congested_samples         %Number of time steps with congestion
        
        congestion_status_per_tenant
        num_congested_samples_per_tenant
        
        
        avg_num_RB_per_tenant   %Average number of RBs used by each tenant.
        avg_num_RB_per_tenant_sample     %Matrix where each column is the average of assigned RB per tenant in each time step
        
        avg_RB_utilisation_per_tenant    %Average real capacity share (RB_utilisation) per tenant.
        avg_RB_utilisation_per_tenant_sample %Matrix where each column is the average of RB utilisation per tenant in each time step
        
        Rb_estimate_per_RB             %Estimate of bit rate per RB (in kb/s) achieved in the cell
        Rb_estimate_per_RB_sample      %Vector with the samples of the Rb estimate per RB in each time step
        
        %NEWPARAMS PER SERVICE
        num_assigned_RB_per_service_per_tenant    %Number of assigned RB to each tenant and service in a time step
        num_assigned_RB_per_service_per_tenant_sample    %3D matrix with the samples of assigned RB per tenant/service in each time step
        total_assigned_RB_per_service_sample      %Matrix where each row includes the samples of total assigned RBs in each time step
        
        bit_rate_assigned_per_service_per_tenant  %Bit rate in the assigned RBs to each tenant per service in a time step 
        bit_rate_assigned_per_service_per_tenant_sample %Matrix with the samples of the total assigned bit rate per service to each tenant in each time step
        total_assigned_bit_rate_per_service_sample    %Vector with the samples of the total assigned bit rate per service in each time step
        
        avg_bit_rate_assigned_per_service_per_tenant    %Average bit rate in the assigned RBs per service to each tenant
        avg_bit_rate_assigned_per_service_per_tenant_sample  %Matrix where each column is the average of assigned bit rate per service per tenant in each time step
        
        avg_num_RB_per_service_per_tenant   %Average number of RBs used by each tenant in each service.
        avg_num_RB_per_service_per_tenant_sample     %Matrix where each column is the average of assigned RB per tenant/service in each time step
        
        avg_RB_utilisation_per_service_per_tenant    %Average real capacity share (RB_utilisation) per tenant/service.
        avg_RB_utilisation_per_service_per_tenant_sample  %Samples of the Average real capacity share (RB_utilisation) per tenant/service.
        
        num_blocks_per_service_per_tenant           %Number of blocked sessions per service per tenant
        num_session_attempts_per_service_per_tenant   %Number of attempted sessions per service tenant
        
        %END NEW PARAMS
        
        num_RBs_per_tenant             %Number of RBs that the PS will allocate to each tenant (in case of RAN slicing at PS)
        
        
        time_next_session_arrival   
        
        num_blocks                     %Number of blocked sessions per tenant
        num_session_attempts           %Number of attempted sessions per tenant
        num_adm_above_SAGBR            %Number of admissions with avg bit rate above SAGBR (per tenant)
        num_rej_below_SAGBR            %Number of blockings with avg bit rate below SAGBR (per tenant)
        offered_load                   %Aggregation of the duration of the sessions multiplied by the bit rate.
        data_volume_per_tenant  %Total data volume (in kbits) transmitted per tenant along the simulation in the cell.
        
        sum_SINR                %Aggregate SINR (in dB) for all the UEs of a tenant
        sum_required_RBs        %Aggregate of the required RBs for all the UEs of a tenant
        num_degraded_UEs        %Number of degraded UEs (with less RBs than required) per tenant
        sum_Rb_degradation      %Percentage of Rb degradation with respect to Rbreq.
        num_samples_UEs         %Number of samples collected of the above statistics (per tenant)
        
        sum_required_RBs_per_service_per_tenant        %Aggregate of the required RBs for all the UEs of a tenant in a given service
        num_degraded_UEs_per_service_per_tenant        %Number of degraded UEs (with less RBs than required) per tenant
        sum_Rb_degradation_per_service_per_tenant      %Percentage of Rb degradation with respect to Rbreq.
        num_samples_UEs_per_service_per_tenant         %Number of samples collected of the above statistics (per tenant)
        
        
        %Parameters for the traffic variation per tenant.
        lambda_ini                      %Initial value of lambda
        delta_lambda                    %Increase of lambda
        freq_traffic_period             %Frequency of the periods 1/T where T is the period in s. (e.g. 1/86400 is 1 day) 
        time_shift                      %Time shift in (s) for the periodicity
        
        service_mix                    %vector with the fraction of session arrivals of each service
        
        
    end
    methods
        function init_BS(obj,config)
            n=obj.id;
            if n==1
                obj.pos=[0,0];
            elseif n<=7
                obj.pos=[config.ISD*cos(pi*(n-2)/3),config.ISD*sin(pi*(n-2)/3)];
            elseif n<=19
                if mod(n,2)==0
                    %Cell at 2ISD distance:
                    obj.pos=[2*config.ISD*cos(((n-8)/2)*pi/3),2*config.ISD*sin(((n-8)/2)*pi/3)];
                else
                    %Cell at 3R distance
                    obj.pos=[3*config.cell_R*cos((pi/6)+((n-9)/2)*pi/3),3*config.cell_R*sin((pi/6)+((n-9)/2)*pi/3)];
                end;
            end;
    
            for s=1:config.num_tenants
                %Assign initially the default values to the generation rate and
                %duration for all the tenants
                %obj.Rbreq(s)=config.traffic_params.Rbreq;
                %obj.duration(s)=config.traffic_params.duration;
                obj.lambda(s)=config.traffic_params.lambda;
                
                %Traffic variation parameters initially set to 0 (no variation).
                obj.lambda_ini(s)=config.traffic_params.lambda;
                obj.delta_lambda(s)=0;
                obj.freq_traffic_period(s)=0;
                obj.time_shift(s)=0;
                
            end;
    
            %Initialize the list of UEs of each tenant.
            obj.UElist=cell(config.num_tenants,1);
            obj.numUEs=zeros(config.num_tenants,1);
    
            %Initialize the number of RBs per BS (by default the total number)
            obj.Ptot=config.Ptot;
            obj.num_RBs=config.num_RBs;
            obj.antenna_gain=config.antenna_gain;
            obj.cell_R=config.cell_R;
            

            %Initialize admission parameters to default value.
            obj.alfa_th=config.admission_params.alfa_th;
            obj.alfa_th_GBR=ones(config.num_tenants,1);
            obj.beta=config.admission_params.beta;
            obj.gamma=config.admission_params.gamma;
    
            obj.C=config.C;   %Theoretical capacity share per tenant as specified in the configuration.
            
            for s=1:config.num_tenants
                obj.num_RBs_per_tenant(s)=obj.num_RBs*obj.C(s);
            end;
            
            
            
            sim_duration_samples=1+config.simulation_duration/config.time_step;
            
            obj.DeltaC=zeros(config.num_tenants,1);
            obj.DeltaC_sample=zeros(sim_duration_samples,config.num_tenants);
            obj.DeltaC_avg_sample=zeros(sim_duration_samples,config.num_tenants);
            
            obj.DeltaCext=zeros(config.num_tenants,1);
            obj.DeltaCext_sample=zeros(sim_duration_samples,config.num_tenants);
            obj.DeltaCext_avg_sample=zeros(sim_duration_samples,config.num_tenants);
            
            obj.Cextra_min=config.admission_params.Cextra_min;
            
            obj.DeltaCbal=zeros(config.num_tenants,1);
            obj.DeltaCbal_sample=zeros(sim_duration_samples,config.num_tenants);
            obj.DeltaCbal_avg_sample=zeros(sim_duration_samples,config.num_tenants);
            
            obj.congestion_status=0;
            obj.num_congested_samples=0;
            obj.congestion_status_per_tenant=zeros(config.num_tenants,1);
            obj.num_congested_samples_per_tenant=zeros(config.num_tenants,1);
            
            
            obj.num_assigned_RB_per_tenant=zeros(config.num_tenants,1);
            obj.avg_num_RB_per_tenant=zeros(config.num_tenants,1);  %Average number of RBs used by each tenant.
            obj.avg_RB_utilisation_per_tenant=zeros(config.num_tenants,1);    %Average real capacity share (RB_utilisation) per tenant.
            
            obj.bit_rate_assigned_per_tenant=zeros(config.num_tenants,1);
            obj.bit_rate_assigned_per_tenant_sample=zeros(sim_duration_samples,config.num_tenants);
            obj.total_assigned_bit_rate_sample=zeros(sim_duration_samples,1);
            
            obj.avg_bit_rate_assigned_per_tenant=zeros(config.num_tenants,1);
            obj.avg_bit_rate_assigned_per_tenant_sample=zeros(sim_duration_samples,config.num_tenants);
            
            
            obj.num_assigned_RB_per_tenant_sample=zeros(sim_duration_samples,config.num_tenants);    
            obj.total_assigned_RB_sample=zeros(sim_duration_samples,1);
            obj.total_assigned_bit_rate_sample=zeros(sim_duration_samples,1);
            obj.avg_num_RB_per_tenant_sample=zeros(sim_duration_samples,config.num_tenants);
            obj.avg_RB_utilisation_per_tenant_sample=zeros(sim_duration_samples,config.num_tenants);
        
            obj.Rb_estimate_per_RB_sample=zeros(sim_duration_samples,1);
            
            obj.num_blocks=zeros(config.num_tenants,1);
            obj.num_session_attempts=zeros(config.num_tenants,1);
            obj.num_adm_above_SAGBR=zeros(config.num_tenants,1);
            obj.num_rej_below_SAGBR=zeros(config.num_tenants,1);
            obj.offered_load=zeros(config.num_tenants,1);
            
            obj.data_volume_per_tenant=zeros(config.num_tenants,1);
            
            obj.freq_index=1;       %By default all the cells with the same freq. If we want a different configuration, it has to be set manually.
            
            %NEWPARAMS PER SERVICE
            obj.service_mix=zeros(config.num_tenants,config.num_services);
            
            obj.num_assigned_RB_per_service_per_tenant=zeros(config.num_tenants,config.num_services);
            obj.num_assigned_RB_per_service_per_tenant_sample=zeros(sim_duration_samples,config.num_tenants,config.num_services);
            obj.total_assigned_RB_per_service_sample=zeros(sim_duration_samples,config.num_services);
            
            obj.bit_rate_assigned_per_service_per_tenant=zeros(config.num_tenants,config.num_services);  
            obj.bit_rate_assigned_per_service_per_tenant_sample=zeros(sim_duration_samples,config.num_tenants,config.num_services);
            obj.total_assigned_bit_rate_per_service_sample=zeros(sim_duration_samples,config.num_services);
        
            obj.avg_bit_rate_assigned_per_service_per_tenant=zeros(config.num_tenants,config.num_services);
            obj.avg_bit_rate_assigned_per_service_per_tenant_sample=zeros(sim_duration_samples,config.num_tenants,config.num_services);  
            
            obj.avg_num_RB_per_service_per_tenant=zeros(config.num_tenants,config.num_services);
            obj.avg_num_RB_per_service_per_tenant_sample=zeros(sim_duration_samples,config.num_tenants,config.num_services);
            
            obj.avg_RB_utilisation_per_service_per_tenant=zeros(config.num_tenants,config.num_services);
            obj.avg_RB_utilisation_per_service_per_tenant_sample=zeros(sim_duration_samples,config.num_tenants,config.num_services);  
            
            obj.num_blocks_per_service_per_tenant=zeros(config.num_tenants,config.num_services);
            obj.num_session_attempts_per_service_per_tenant=zeros(config.num_tenants,config.num_services);
            %END NEW PARAMS
            
            obj.sum_SINR=zeros(config.num_tenants,1);
            obj.sum_required_RBs=zeros(config.num_tenants,1);
            obj.num_degraded_UEs=zeros(config.num_tenants,1);
            obj.sum_Rb_degradation=zeros(config.num_tenants,1);
            obj.num_samples_UEs=zeros(config.num_tenants,1);
            
            obj.sum_required_RBs_per_service_per_tenant=zeros(config.num_tenants,config.num_services);        
            obj.num_degraded_UEs_per_service_per_tenant=zeros(config.num_tenants,config.num_services);        
            obj.sum_Rb_degradation_per_service_per_tenant=zeros(config.num_tenants,config.num_services);       
            obj.num_samples_UEs_per_service_per_tenant=zeros(config.num_tenants,config.num_services);          
        
            
        end;
        
        function init_radio_and_next_arrivals(obj,config)
            
            obj.time_next_session_arrival=-(1./obj.lambda).*log(1-rand(1,config.num_tenants));
            
            %Radio (power per RB and PIRE per RB)
            obj.P_RB=obj.Ptot-10*log10(obj.num_RBs); 
            obj.PIRE_RB=obj.P_RB+obj.antenna_gain;     %dBm
            
            %Compute a reference SNR to estimate the RBs, based on e.g. a
            %fraction of the cell radius without interference
            L=prop_model(config, 0.8*config.cell_R,1); %Loss without shadowing
            SNR_ref=obj.PIRE_RB-L-10*log10(config.Pnoise_RB); %dB
            SNR_ref=power(10,0.1*SNR_ref);
            S=SpEff(config,SNR_ref);
            obj.Rb_estimate_per_RB=S*config.B_RB;  %in kb/s
            
        end;
        
        function admission_result=admission(obj,tenant,Rbreq,config,service)
            
            delta_rho=Rbreq/obj.Rb_estimate_per_RB; 
            arp=config.services.ARP(service);
            
            switch (config.RAN_slicing_algorithm) 
                case config.RAN_SLICING_AT_SPECTRUM_PLANNING
                    %There is only one check at tenant level, accounting
                    %for the RBs of the tenant.
                    
                    %Check: Usage at tenant level
                    if (sum(obj.avg_num_RB_per_service_per_tenant(tenant,(config.services.type(:)==config.GBR)&(config.services.ARP(:)<=arp)))+ delta_rho)> obj.num_RBs_per_tenant(tenant)*obj.alfa_th_GBR(tenant)*obj.alfa_th
                        %cannot be admitted
                        admission_result=0;
                    else
                        %admission OK
                        admission_result=1;
                    end;
                case config.RAN_SLICING_AT_PS
                    %There is only one check at tenant level, accounting
                    %for the RBs of the tenant.
                    
                    %Check: Usage at tenant level
                    if (sum(obj.avg_num_RB_per_service_per_tenant(tenant,(config.services.type(:)==config.GBR)&(config.services.ARP(:)<=arp)))+ delta_rho)> obj.num_RBs_per_tenant(tenant)*obj.alfa_th_GBR(tenant)*obj.alfa_th
                        %cannot be admitted
                        admission_result=0;
                    else
                        %admission OK
                        admission_result=1;
                    end;
                case config.RAN_SLICING_AT_AC
                    switch (config.AC_algorithm)
                        case (config.NO_SLICING)
                            %Check only the aggregate usage of GBR Rbs among all tenants
                            %rho_aggr=sum(obj.avg_num_RB_per_tenant);
                            rho_aggr=sum(sum(obj.avg_num_RB_per_service_per_tenant(:,(config.services.type(:)==config.GBR)&(config.services.ARP(:)<=arp))));
                    
                            if (rho_aggr+delta_rho) >  obj.num_RBs*obj.alfa_th
                                %cannot be admitted
                                admission_result=0;
                            else
                                %admission OK
                                admission_result=1;
                            end;
                    
                        case (config.SLICING_NO_DELTA)
                            %1st check: Aggregate usage among all tenants:
                            %rho_aggr=sum(obj.avg_num_RB_per_tenant);
                            rho_aggr=sum(sum(obj.avg_num_RB_per_service_per_tenant(:,(config.services.type(:)==config.GBR)&(config.services.ARP(:)<=arp))));
                            if (rho_aggr+delta_rho) >  obj.num_RBs*obj.alfa_th
                                %cannot be admitted
                                admission_result=0;
                            else
                                %2nd_check: Usage at tenant level
                                %if (obj.avg_num_RB_per_tenant(tenant)+ delta_rho) > obj.num_RBs*obj.alfa_th*obj.C(tenant)
                                if (sum(obj.avg_num_RB_per_service_per_tenant(tenant,(config.services.type(:)==config.GBR)&(config.services.ARP(:)<=arp)))+ delta_rho) > obj.num_RBs*obj.alfa_th*obj.C(tenant)
                                    %cannot be admitted
                                    admission_result=0;
                                else
                                    %admission OK
                                    admission_result=1;
                                end;
                            end;
                
                        case (config.SLICING_DELTA)    
                            %1st check: Aggregate usage among all tenants:
                            %rho_aggr=sum(obj.avg_num_RB_per_tenant);
                            rho_aggr=sum(sum(obj.avg_num_RB_per_service_per_tenant(:,(config.services.type(:)==config.GBR)&(config.services.ARP(:)<=arp))));
                            if (rho_aggr+delta_rho) >  obj.num_RBs*obj.alfa_th
                                %cannot be admitted
                                admission_result=0;
                            else
                                %2nd_check: Usage at tenant level
                                %if (obj.avg_num_RB_per_tenant(tenant)+ delta_rho) > obj.num_RBs*obj.alfa_th*(obj.C(tenant)+obj.DeltaC(tenant))
                                if (sum(obj.avg_num_RB_per_service_per_tenant(tenant,(config.services.type(:)==config.GBR)&(config.services.ARP(:)<=arp)))+ delta_rho) > obj.num_RBs*obj.alfa_th*(obj.C(tenant)+obj.DeltaC(tenant))
                                    %cannot be admitted
                                    admission_result=0;
                                else
                                    %admission OK
                                    admission_result=1;
                                end;
                            end;
            
                    end;
            end;
        end;
        

       
    end
end