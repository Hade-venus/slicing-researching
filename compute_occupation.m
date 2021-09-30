function compute_occupation( config, BS )
%COMPUTE_OCCUPATION Summary of this function goes here
%   Detailed explanation goes here

    %This function computes the total occupation of the different BSs at
    %the current time step.
    
    
    
    for n=1:config.num_cells
        num_assigned_RB_per_cell_previous=sum(BS(n).num_assigned_RB_per_tenant);
    end;
    %num_assigned_RB_per_cell_previous=zeros(config.num_cells,1);
    num_assigned_RB_per_cell_current=zeros(config.num_cells,1);
    
    
    switch (config.RAN_slicing_algorithm)
    
        case config.RAN_SLICING_AT_AC
            %Slicing at AC: 
            %It means that there is only a L3 descriptor specifying the %
            %of PRBs for each slice/tenant.
            %However, the PS does not make distinctions between
            %slices/tenants.
            %Characteristics:
            %1) Intercell interference: Received interference can come from
            %any RB of the neighbor cell, assigned to any tenant. Then, average interference
            %is weighted by: (sum(BS(j).num_assigned_RB_per_tenant)/BS(j).num_RBs)
            
            %2) Resource assignment: There are no limits in the amount of RBs that can be assigned on a per tenant basis. The only limit is that the total number of RBs assigned to all tenants should not exceed 
            %the limits of the cell. If this limit is exceeded,
            %there is congestion for the cell, so all the UEs are
            %degraded.
    
    
    
            end_process=0;
            num_iter=1;
    
    
            while ~end_process
                %Step 1: Compute the SINR seen by each UE of each BS in one RB, and
                %from this the number of required RBs by each UE.
    
                for n=1:config.num_cells
                    for s=1:config.num_tenants
                        for i=1:BS(n).numUEs(s)
                            Rx_level=BS(n).PIRE_RB-BS(n).UElist{s}(i).PL(n); %dBm
                            Rx_level=power(10,0.1*Rx_level); %mW
                
                            %Interference:
                            Interf=config.Pnoise_RB;  %mW
                            for j=1:config.num_cells
                                if (j~=n) && (BS(j).freq_index==BS(n).freq_index)
                                    %We only account for interference if the two cells
                                    %use the same frequency.
                        
                                    %the generated interference is weighted by the
                                    %total number of assigned_RBs
                                    Interf=Interf+power(10,0.1*(BS(j).PIRE_RB-BS(n).UElist{s}(i).PL(j)))*(sum(BS(j).num_assigned_RB_per_tenant)/BS(j).num_RBs);
                               end;
                            end;
                
                            BS(n).UElist{s}(i).SINR=Rx_level/Interf;
                            BS(n).UElist{s}(i).S=SpEff(config, BS(n).UElist{s}(i).SINR);
                
                            BS(n).UElist{s}(i).assigned_RBs=0;  %Initially, no assigned RBs are considered.
                            BS(n).UElist{s}(i).assigned_bit_rate=0;
                
                            %Account required_RBs only for GBR and active UEs:
                            if (config.services.type(BS(n).UElist{s}(i).service)==config.GBR) && (BS(n).UElist{s}(i).activity_state==1)
                                BS(n).UElist{s}(i).required_RBs=BS(n).UElist{s}(i).Rbreq/(config.B_RB*BS(n).UElist{s}(i).S);
                            else
                                BS(n).UElist{s}(i).required_RBs=0;
                            end;
                        end;
                    end;
                end;
    
                %Step 2: Determine the number of assigned RBs
                for n=1:config.num_cells
                    BS(n).num_assigned_RB_per_tenant=zeros(config.num_tenants,1);
                    BS(n).num_assigned_RB_per_service_per_tenant=zeros(config.num_tenants,config.num_services);
        
                    BS(n).congestion_status=0;
                    arp=min(config.services.ARP);
        
                    %Assignment of GBR services in increasing order of ARP (when
                    %reaching congestion the process ends).
                    while (arp<=max(config.services.ARP)) && (BS(n).congestion_status==0)
        
                        for s=1:config.num_tenants
                            for i=1:BS(n).numUEs(s)
                                if (config.services.type(BS(n).UElist{s}(i).service)==config.GBR) && (config.services.ARP(BS(n).UElist{s}(i).service)==arp) && (BS(n).UElist{s}(i).activity_state==1)
                                    BS(n).num_assigned_RB_per_tenant(s)=BS(n).num_assigned_RB_per_tenant(s)+BS(n).UElist{s}(i).required_RBs;
                                    BS(n).num_assigned_RB_per_service_per_tenant(s,BS(n).UElist{s}(i).service)=BS(n).num_assigned_RB_per_service_per_tenant(s,BS(n).UElist{s}(i).service)+BS(n).UElist{s}(i).required_RBs;
                                    BS(n).UElist{s}(i).assigned_RBs=BS(n).UElist{s}(i).required_RBs;
                                end;
                            end;
                        end;
            
                        %Check congestion status after having assigned RBs to the users
                        %with ARP=arp.
                        if sum(BS(n).num_assigned_RB_per_tenant)>BS(n).num_RBs
                            %Congestion!!!
                            %Criterion: we split the excess proportionally between all the
                            %UEs with ARP=arp
                            BS(n).congestion_status=1;
                
                            %Reference: Number of assignedRBs of the GBR UEs with
                            %ARP=arp
                            num_assigned_RBs_GBR_arp=sum(sum(BS(n).num_assigned_RB_per_service_per_tenant(:,(config.services.type(:)==config.GBR)&(config.services.ARP(:)==arp))));
                
                            excess_fraction=(sum(BS(n).num_assigned_RB_per_tenant)-BS(n).num_RBs)/num_assigned_RBs_GBR_arp;
                            %Meaning: if excess_fraction = 0.1, it means that the excess is
                            %the 10% of all the requirements of the GBR UEs with ARP, so we remove to each UE this
                            %fraction
                            for s=1:config.num_tenants
                                for i=1:BS(n).numUEs(s)
                                    if (config.services.type(BS(n).UElist{s}(i).service)==config.GBR) && (config.services.ARP(BS(n).UElist{s}(i).service)==arp) && (BS(n).UElist{s}(i).activity_state==1)
                                        %Modify new assignment to the UE:
                                        BS(n).UElist{s}(i).assigned_RBs=BS(n).UElist{s}(i).required_RBs*(1-excess_fraction); 
    
                                        %Remove the excess assignment from the total:
                                        BS(n).num_assigned_RB_per_tenant(s)=BS(n).num_assigned_RB_per_tenant(s)-BS(n).UElist{s}(i).required_RBs*excess_fraction;
                                        BS(n).num_assigned_RB_per_service_per_tenant(s,BS(n).UElist{s}(i).service)=BS(n).num_assigned_RB_per_service_per_tenant(s,BS(n).UElist{s}(i).service)-BS(n).UElist{s}(i).required_RBs*excess_fraction;
        
                                    end;
                                end;                    
                            end;
                        else
                            BS(n).congestion_status=0;
                        end;
                        arp=arp+1;
                    end;
        
                    if (BS(n).congestion_status==0)
                        %Assignment of non-GBR:
                        num_available_RBs=BS(n).num_RBs-sum(BS(n).num_assigned_RB_per_tenant);
            
                        %1.- Compute the denominator of the resource sharing. It is the
                        %aggregate of (1/priority)^exponent for all the UEs Non-GBR and
                        %active
                        sum_denominator=0;
                        for s=1:config.num_tenants
                            for i=1:BS(n).numUEs(s)
                                if (config.services.type(BS(n).UElist{s}(i).service)==config.NonGBR) && (BS(n).UElist{s}(i).activity_state==1)
                                    sum_denominator=sum_denominator+power(1/config.services.priority(BS(n).UElist{s}(i).service),config.exponent_priority_nonGBR);
                                end;
                            end;
                        end;
                        %2.- Compute the amount of assigned RBs: 
                        if (sum_denominator>0)
                            for s=1:config.num_tenants
                                for i=1:BS(n).numUEs(s)
                                    if (config.services.type(BS(n).UElist{s}(i).service)==config.NonGBR) && (BS(n).UElist{s}(i).activity_state==1)
                                        BS(n).UElist{s}(i).assigned_RBs=num_available_RBs*power(1/config.services.priority(BS(n).UElist{s}(i).service),config.exponent_priority_nonGBR)/sum_denominator;

                                        BS(n).num_assigned_RB_per_tenant(s)=BS(n).num_assigned_RB_per_tenant(s)+BS(n).UElist{s}(i).assigned_RBs;
                                        BS(n).num_assigned_RB_per_service_per_tenant(s,BS(n).UElist{s}(i).service)=BS(n).num_assigned_RB_per_service_per_tenant(s,BS(n).UElist{s}(i).service)+BS(n).UElist{s}(i).assigned_RBs;
                                    end;
                                end;
                            end;
                        end;
                    end;
        
                    num_assigned_RB_per_cell_current(n)=sum(BS(n).num_assigned_RB_per_tenant);
        
                    if BS(n).congestion_status==1
                        BS(n).congestion_status_per_tenant=ones(config.num_tenants,1);
                    else 
                        BS(n).congestion_status_per_tenant=zeros(config.num_tenants,1);
                    end;
        
        
                 end;
     
                %Step 3: Check finalisation condition:
     
                error=abs(num_assigned_RB_per_cell_current-num_assigned_RB_per_cell_previous);
     
                error_total=max(error);
                if (error_total<0.01)||(num_iter==100)
                    end_process=1;
         
                    %Compute the final bit rate assigned to each tenant.
                    for n=1:config.num_cells
                        BS(n).bit_rate_assigned_per_tenant=zeros(config.num_tenants,1);
                        BS(n).bit_rate_assigned_per_service_per_tenant=zeros(config.num_tenants,config.num_services);
                        for s=1:config.num_tenants
                            for i=1:BS(n).numUEs(s)
                                BS(n).UElist{s}(i).assigned_bit_rate=BS(n).UElist{s}(i).S*BS(n).UElist{s}(i).assigned_RBs*config.B_RB;
                                BS(n).bit_rate_assigned_per_tenant(s)=BS(n).bit_rate_assigned_per_tenant(s)+BS(n).UElist{s}(i).assigned_bit_rate;
        
                                BS(n).bit_rate_assigned_per_service_per_tenant(s,BS(n).UElist{s}(i).service)=BS(n).bit_rate_assigned_per_service_per_tenant(s,BS(n).UElist{s}(i).service)+BS(n).UElist{s}(i).assigned_bit_rate;
                            end;
                        end;
                    end;
                    if (num_iter==100) && (error_total>0.01)
                        fprintf('Warning: occupation not converged after %d iterations... error: %f\n',num_iter,error_total);
                    end;
                else
                    num_assigned_RB_per_cell_previous=num_assigned_RB_per_cell_current;
                    num_iter=num_iter+1;
                end;
     
            end;  
        case config.RAN_SLICING_AT_PS
            %Slicing at PS
            %In this case, there is an L2 descriptor that specifies the number of PRBs that the PS can allocate for each slice/tenant.
            %This is included in parameter: BS(j).num_RBs_per_tenant
            %If one slice/tenant does not use all the PRBs, they can be used by
            %the other slices/tenants

            %Characteristics:
            %1) Intercell interference: Received interference can come from
            %any RB of the neighbor cell, assigned to any tenant. Then, average interference
            %is weighted by: (sum(BS(j).num_assigned_RB_per_tenant)/BS(j).num_RBs)
            
            %2) Resource assignment: Congestion occurs if the total required RBs exceed the available RBs. In that case, congestion is only experienced
            %by the tenants whose required RBs exceed the value of
            %BS(j).num_RBs_per_tenant(s). If not, the "extra RBs" left
            %for this tenant can be used by the other tenants.
            end_process=0;
            num_iter=1;
    
            while ~end_process
                
                %Step 1: Compute the SINR seen by each UE of each BS in one RB, and
                %from this the number of required RBs by each UE.
                
                for n=1:config.num_cells
                    for s=1:config.num_tenants
                        for i=1:BS(n).numUEs(s)
                            Rx_level=BS(n).PIRE_RB-BS(n).UElist{s}(i).PL(n); %dBm
                            Rx_level=power(10,0.1*Rx_level); %mW
                
                            %Interference:
                            Interf=config.Pnoise_RB;  %mW
                            for j=1:config.num_cells
                                if (j~=n) && (BS(j).freq_index==BS(n).freq_index)
                                    %We only account for interference if the two cells
                                    %use the same frequency.
                        
                                    %the generated interference is weighted by the
                                    %total number of assigned_RBs 
                                    Interf=Interf+power(10,0.1*(BS(j).PIRE_RB-BS(n).UElist{s}(i).PL(j)))*(sum(BS(j).num_assigned_RB_per_tenant)/BS(j).num_RBs);
                                end;
                            end;
                
                            BS(n).UElist{s}(i).SINR=Rx_level/Interf;
                            BS(n).UElist{s}(i).S=SpEff(config, BS(n).UElist{s}(i).SINR);
                        
                            BS(n).UElist{s}(i).assigned_RBs=0;  %Initially, no assigned RBs are considered.
                            BS(n).UElist{s}(i).assigned_bit_rate=0;
                        
                            %Account required_RBs only for GBR and active UEs:
                            if (config.services.type(BS(n).UElist{s}(i).service)==config.GBR) && (BS(n).UElist{s}(i).activity_state==1)
                                BS(n).UElist{s}(i).required_RBs=BS(n).UElist{s}(i).Rbreq/(config.B_RB*BS(n).UElist{s}(i).S);
                            else
                                BS(n).UElist{s}(i).required_RBs=0;
                            end;
                        end;
                   end;
                end;
                
                %Step 2: Determine the number of assigned RBs, following
                %the same process as in the "Slicing at AC" but considering
                %each tenant s separately from the rest, and with the limit
                %of BS(n).num_RB_per_tenant(s).
                %The extra RBs not used by a tenant are distributed among
                %the other tenants based on the ARPs of their services.
                for n=1:config.num_cells
                    BS(n).num_assigned_RB_per_tenant=zeros(config.num_tenants,1);
                    BS(n).num_assigned_RB_per_service_per_tenant=zeros(config.num_tenants,config.num_services);
                    BS(n).congestion_status=0;
                    BS(n).congestion_status_per_tenant=zeros(config.num_tenants,1);
                    
                    extra_RBs=0;
                    
                    for s=1:config.num_tenants
                        %Get the minimum and maximum arp of the UEs of the
                        %tenant                        
                        arpmin=max(config.services.ARP); %Initial value
                        arpmax=0;
                        for i=1:BS(n).numUEs(s)
                            if config.services.ARP(BS(n).UElist{s}(i).service)<arpmin
                                arpmin=config.services.ARP(BS(n).UElist{s}(i).service);
                            end;
                            if config.services.ARP(BS(n).UElist{s}(i).service)>arpmax
                                arpmax=config.services.ARP(BS(n).UElist{s}(i).service);
                            end;
                        end;
                        arp=arpmin;
                        while (arp<=arpmax) && (BS(n).congestion_status_per_tenant(s)==0)
                            for i=1:BS(n).numUEs(s)
                                if (config.services.type(BS(n).UElist{s}(i).service)==config.GBR) && (config.services.ARP(BS(n).UElist{s}(i).service)==arp) && (BS(n).UElist{s}(i).activity_state==1)
                                    BS(n).num_assigned_RB_per_tenant(s)=BS(n).num_assigned_RB_per_tenant(s)+BS(n).UElist{s}(i).required_RBs;
                                    BS(n).num_assigned_RB_per_service_per_tenant(s,BS(n).UElist{s}(i).service)=BS(n).num_assigned_RB_per_service_per_tenant(s,BS(n).UElist{s}(i).service)+BS(n).UElist{s}(i).required_RBs;
                                    BS(n).UElist{s}(i).assigned_RBs=BS(n).UElist{s}(i).required_RBs;
                                end;
                            end;
                            %Check congestion status after having assigned
                            %RBs to the users of tenant s
                            %with ARP=arp.
                            if BS(n).num_assigned_RB_per_tenant(s)>BS(n).num_RBs_per_tenant(s)
                                %Congestion!!!
                                %Criterion: we split the excess proportionally between all the
                                %UEs with ARP=arp
                                BS(n).congestion_status_per_tenant(s)=1;
                
                                %Reference: Number of assignedRBs of the GBR UEs with
                                %ARP=arp
                                num_assigned_RBs_GBR_arp=sum(BS(n).num_assigned_RB_per_service_per_tenant(s,(config.services.type(:)==config.GBR)&(config.services.ARP(:)==arp)));
                
                                excess_fraction=(BS(n).num_assigned_RB_per_tenant(s)-BS(n).num_RBs_per_tenant(s))/num_assigned_RBs_GBR_arp;
                                %Meaning: if excess_fraction = 0.1, it means that the excess is
                                %the 10% of all the requirements of the GBR UEs with ARP, so we remove to each UE this
                                %fraction
                            
                                for i=1:BS(n).numUEs(s)
                                    if (config.services.type(BS(n).UElist{s}(i).service)==config.GBR) && (config.services.ARP(BS(n).UElist{s}(i).service)==arp) && (BS(n).UElist{s}(i).activity_state==1)
                                        %Modify new assignment to the UE:
                                        BS(n).UElist{s}(i).assigned_RBs=BS(n).UElist{s}(i).required_RBs*(1-excess_fraction); 
    
                                        %Remove the excess assignment from the total:
                                        BS(n).num_assigned_RB_per_tenant(s)=BS(n).num_assigned_RB_per_tenant(s)-BS(n).UElist{s}(i).required_RBs*excess_fraction;
                                        BS(n).num_assigned_RB_per_service_per_tenant(s,BS(n).UElist{s}(i).service)=BS(n).num_assigned_RB_per_service_per_tenant(s,BS(n).UElist{s}(i).service)-BS(n).UElist{s}(i).required_RBs*excess_fraction;
        
                                    end;
                                end;                    
                            
                            else
                                BS(n).congestion_status_per_tenant(s)=0;
                            end;
                            arp=arp+1;
                        end;    
                        if (BS(n).congestion_status_per_tenant(s)==0)
                            %Assignment of non-GBR:
                            num_available_RBs=BS(n).num_RBs_per_tenant(s)-BS(n).num_assigned_RB_per_tenant(s);
            
                            %1.- Compute the denominator of the resource sharing. It is the
                            %aggregate of (1/priority)^exponent for all the UEs Non-GBR and
                            %active
                            sum_denominator=0;

                            for i=1:BS(n).numUEs(s)
                                if (config.services.type(BS(n).UElist{s}(i).service)==config.NonGBR) && (BS(n).UElist{s}(i).activity_state==1)
                                    sum_denominator=sum_denominator+power(1/config.services.priority(BS(n).UElist{s}(i).service),config.exponent_priority_nonGBR);
                                end;
                            end;

                            %2.- Compute the amount of assigned RBs: 
                            if (sum_denominator>0)
                                for i=1:BS(n).numUEs(s)
                                    if (config.services.type(BS(n).UElist{s}(i).service)==config.NonGBR) && (BS(n).UElist{s}(i).activity_state==1)
                                        BS(n).UElist{s}(i).assigned_RBs=num_available_RBs*power(1/config.services.priority(BS(n).UElist{s}(i).service),config.exponent_priority_nonGBR)/sum_denominator;

                                        BS(n).num_assigned_RB_per_tenant(s)=BS(n).num_assigned_RB_per_tenant(s)+BS(n).UElist{s}(i).assigned_RBs;
                                        BS(n).num_assigned_RB_per_service_per_tenant(s,BS(n).UElist{s}(i).service)=BS(n).num_assigned_RB_per_service_per_tenant(s,BS(n).UElist{s}(i).service)+BS(n).UElist{s}(i).assigned_RBs;
                                    end;
                                end;
                            end;
                        end;
                        
                        %Update the extra_RBs not used by tenant s
                        extra_RBs=extra_RBs+max(BS(n).num_RBs_per_tenant(s)-BS(n).num_assigned_RB_per_tenant(s),0);
                        
                    end;
                    
                    %Distribute the extra RBs among the different tenants,
                    %in accordance with the ARP priorities (similar process
                    %as in the "slicing at AC":
                    
                    if (extra_RBs>1E-8)
                        %This condition is to avoid small differences associated
                        %to decimals in the computation (e.g. extra_RBs in
                        %the order of 1E-15 which actually mean that the
                        %extra_RBs are 0).
                        num_available_RBs=extra_RBs;
                    else
                        num_available_RBs=0;
                    end;
                    
                    arp=min(config.services.ARP);
        
                    %Assignment of GBR services in increasing order of ARP (when
                    %reaching congestion the process ends).
                    while (arp<=max(config.services.ARP)) && (num_available_RBs>0)
                        assigned_extra_RBs=0;    
                        for s=1:config.num_tenants
                            for i=1:BS(n).numUEs(s)
                                %Select the GBR users with ARP=arp,
                                %activity_state=1 and with less assigned
                                %RBs than required.
                                if (config.services.type(BS(n).UElist{s}(i).service)==config.GBR) && (config.services.ARP(BS(n).UElist{s}(i).service)==arp) && (BS(n).UElist{s}(i).activity_state==1) && (BS(n).UElist{s}(i).assigned_RBs<BS(n).UElist{s}(i).required_RBs)
                                    %In principle we assign all the RBs
                                    %that the UE needs until getting
                                    %"required_RBs". If we exceed the
                                    %remaining RBs, we will re-distribute
                                    %among all.
                                    assigned_extra_RBs=assigned_extra_RBs+BS(n).UElist{s}(i).required_RBs-BS(n).UElist{s}(i).assigned_RBs;
                                    BS(n).num_assigned_RB_per_tenant(s)=BS(n).num_assigned_RB_per_tenant(s)+BS(n).UElist{s}(i).required_RBs-BS(n).UElist{s}(i).assigned_RBs;
                                    BS(n).num_assigned_RB_per_service_per_tenant(s,BS(n).UElist{s}(i).service)=BS(n).num_assigned_RB_per_service_per_tenant(s,BS(n).UElist{s}(i).service)+BS(n).UElist{s}(i).required_RBs-BS(n).UElist{s}(i).assigned_RBs;
                                    BS(n).UElist{s}(i).assigned_RBs=BS(n).UElist{s}(i).required_RBs;
                                end;
                            end;
                        end;
            
                        %Check congestion status after having assigned the extra RBs to the users
                        %with ARP=arp.
                        if assigned_extra_RBs>num_available_RBs
                            %Criterion: we split the excess proportionally between all the
                            %UEs with ARP=arp
                            
                            %Reference: Number of assignedRBs of the GBR UEs with
                            %ARP=arp
                            num_assigned_RBs_GBR_arp=sum(sum(BS(n).num_assigned_RB_per_service_per_tenant(:,(config.services.type(:)==config.GBR)&(config.services.ARP(:)==arp))));
                
                            excess_fraction=(assigned_extra_RBs-num_available_RBs)/num_assigned_RBs_GBR_arp;
                            %Meaning: if excess_fraction = 0.1, it means that the excess is
                            %the 10% of all the requirements of the GBR UEs with ARP, so we remove to each UE this
                            %fraction
                            for s=1:config.num_tenants
                                for i=1:BS(n).numUEs(s)
                                    if (config.services.type(BS(n).UElist{s}(i).service)==config.GBR) && (config.services.ARP(BS(n).UElist{s}(i).service)==arp) && (BS(n).UElist{s}(i).activity_state==1)
                                        %Modify new assignment to the UE:
                                        BS(n).UElist{s}(i).assigned_RBs=BS(n).UElist{s}(i).required_RBs*(1-excess_fraction); 
    
                                        %Remove the excess assignment from the total:
                                        BS(n).num_assigned_RB_per_tenant(s)=BS(n).num_assigned_RB_per_tenant(s)-BS(n).UElist{s}(i).required_RBs*excess_fraction;
                                        BS(n).num_assigned_RB_per_service_per_tenant(s,BS(n).UElist{s}(i).service)=BS(n).num_assigned_RB_per_service_per_tenant(s,BS(n).UElist{s}(i).service)-BS(n).UElist{s}(i).required_RBs*excess_fraction;
        
                                    end;
                                end;                    
                            end;
                        end;
                        arp=arp+1;
                        num_available_RBs=max(num_available_RBs-assigned_extra_RBs,0);
                    end;
        
                    if num_available_RBs>0
                        %Assignment of non-GBR:
                       
                        %1.- Compute the denominator of the resource sharing. It is the
                        %aggregate of (1/priority)^exponent for all the UEs Non-GBR and
                        %active
                        sum_denominator=0;
                        for s=1:config.num_tenants
                            for i=1:BS(n).numUEs(s)
                                if (config.services.type(BS(n).UElist{s}(i).service)==config.NonGBR) && (BS(n).UElist{s}(i).activity_state==1)
                                    sum_denominator=sum_denominator+power(1/config.services.priority(BS(n).UElist{s}(i).service),config.exponent_priority_nonGBR);
                                end;
                            end;
                        end;
                        %2.- Compute the amount of assigned RBs: 
                        if (sum_denominator>0)
                            for s=1:config.num_tenants
                                for i=1:BS(n).numUEs(s)
                                    if (config.services.type(BS(n).UElist{s}(i).service)==config.NonGBR) && (BS(n).UElist{s}(i).activity_state==1)
                                        extra_assignment=num_available_RBs*power(1/config.services.priority(BS(n).UElist{s}(i).service),config.exponent_priority_nonGBR)/sum_denominator;
                                        BS(n).UElist{s}(i).assigned_RBs=BS(n).UElist{s}(i).assigned_RBs+extra_assignment;

                                        BS(n).num_assigned_RB_per_tenant(s)=BS(n).num_assigned_RB_per_tenant(s)+extra_assignment;
                                        BS(n).num_assigned_RB_per_service_per_tenant(s,BS(n).UElist{s}(i).service)=BS(n).num_assigned_RB_per_service_per_tenant(s,BS(n).UElist{s}(i).service)+extra_assignment;
                                    end;
                                end;
                            end;
                        end;
                    end;
                    
                    num_assigned_RB_per_cell_current(n)=sum(BS(n).num_assigned_RB_per_tenant);
                    
                end;
                
                %Step 3: Check finalisation condition:
     
                error=abs(num_assigned_RB_per_cell_current-num_assigned_RB_per_cell_previous);
     
                error_total=max(error);
                if (error_total<0.01)||(num_iter==100)
                    end_process=1;
         
                    %Compute the final bit rate assigned to each tenant.
                    for n=1:config.num_cells
                        BS(n).bit_rate_assigned_per_tenant=zeros(config.num_tenants,1);
                        BS(n).bit_rate_assigned_per_service_per_tenant=zeros(config.num_tenants,config.num_services);
                        
                        BS(n).congestion_status=0;
                        BS(n).congestion_status_per_tenant=zeros(config.num_tenants,1);
                        for s=1:config.num_tenants
                            for i=1:BS(n).numUEs(s)
                                BS(n).UElist{s}(i).assigned_bit_rate=BS(n).UElist{s}(i).S*BS(n).UElist{s}(i).assigned_RBs*config.B_RB;
                                BS(n).bit_rate_assigned_per_tenant(s)=BS(n).bit_rate_assigned_per_tenant(s)+BS(n).UElist{s}(i).assigned_bit_rate;
        
                                BS(n).bit_rate_assigned_per_service_per_tenant(s,BS(n).UElist{s}(i).service)=BS(n).bit_rate_assigned_per_service_per_tenant(s,BS(n).UElist{s}(i).service)+BS(n).UElist{s}(i).assigned_bit_rate;
                                
                                %Note: since the congestion status may have
                                %changed depending on the assignment of
                                %extra RBs, we recompute it here,
                                %considering that there will be congestion
                                %if some UEs get a bit rate below their
                                %requirement.
                                
                                if (BS(n).UElist{s}(i).Rbreq-BS(n).UElist{s}(i).assigned_bit_rate)>1E-8
                                    %Note: To avoid errors due to decimals,
                                    %we put a threshold 1E-8 in the
                                    %difference between Rbreq and
                                    %assigned_bit_rate to declare
                                    %congestion
                                    BS(n).congestion_status_per_tenant(s)=1;
                                    BS(n).congestion_status=1;  %If one tenant is congested, we assume congestion in the whole cell.
                                end;
                            end;
                        end;
                    end;
                    if (num_iter==100) && (error_total>0.01)
                        fprintf('Warning: occupation not converged after %d iterations... error: %f\n',num_iter,error_total);
                    end;
                else
                    num_assigned_RB_per_cell_previous=num_assigned_RB_per_cell_current;
                    num_iter=num_iter+1;
                end;                
    
            end;
        
        case config.RAN_SLICING_AT_SPECTRUM_PLANNING
            %Slicing at Spectrum Planning
            %It means that there is an L1 descriptor specifying orthogonal
            %RBs for each slice/tenant. These RBs are the same in all the
            %cells.
             %Characteristics:
            %1) Intercell interference: only the RBs of the specific tenant in
            %the neighbour cell generate interference. Then, average interference
            %is weighted by: (BS(j).num_assigned_RB_per_tenant(s)/BS(j).num_RBs_per_tenant(s))
            
            %2) Resource assignment: Each tenant can only get up to
            %BS(j).num_RBs_per_tenant(s)  RBs. If this value is exceeded,
            %there is congestion for the specific tenant, so its UEs are
            %degraded.
            
            end_process=0;
            num_iter=1;
    
            while ~end_process
                %Step 1: Compute the SINR seen by each UE of each BS in one RB, and
                %from this the number of required RBs by each UE.
    
                for n=1:config.num_cells
                    for s=1:config.num_tenants
                        for i=1:BS(n).numUEs(s)
                            Rx_level=BS(n).PIRE_RB-BS(n).UElist{s}(i).PL(n); %dBm
                            Rx_level=power(10,0.1*Rx_level); %mW
                
                            %Interference:
                            Interf=config.Pnoise_RB;  %mW
                            for j=1:config.num_cells
                                if (j~=n) && (BS(j).freq_index==BS(n).freq_index)
                                    %We only account for interference if the two cells
                                    %use the same frequency.
                        
                                    %the generated interference is weighted by the
                                    %total number of assigned_RBs
                                    Interf=Interf+power(10,0.1*(BS(j).PIRE_RB-BS(n).UElist{s}(i).PL(j)))*(BS(j).num_assigned_RB_per_tenant(s)/BS(j).num_RBs_per_tenant(s));
                               end;
                            end;
                
                            BS(n).UElist{s}(i).SINR=Rx_level/Interf;
                            BS(n).UElist{s}(i).S=SpEff(config, BS(n).UElist{s}(i).SINR);
                
                            BS(n).UElist{s}(i).assigned_RBs=0;  %Initially, no assigned RBs are considered.
                            BS(n).UElist{s}(i).assigned_bit_rate=0;
                
                            %Account required_RBs only for GBR and active UEs:
                            if (config.services.type(BS(n).UElist{s}(i).service)==config.GBR) && (BS(n).UElist{s}(i).activity_state==1)
                                BS(n).UElist{s}(i).required_RBs=BS(n).UElist{s}(i).Rbreq/(config.B_RB*BS(n).UElist{s}(i).S);
                            else
                                BS(n).UElist{s}(i).required_RBs=0;
                            end;
                        end;
                    end;
                end;
                
                %Step 2: Determine the number of assigned RBs. It is
                %equivalent to the case "Slicing at PS" but without
                %considering reallocation of extra RBs.
                for n=1:config.num_cells
                    BS(n).num_assigned_RB_per_tenant=zeros(config.num_tenants,1);
                    BS(n).num_assigned_RB_per_service_per_tenant=zeros(config.num_tenants,config.num_services);
                    BS(n).congestion_status=0;
                    BS(n).congestion_status_per_tenant=zeros(config.num_tenants,1);
                    
                    for s=1:config.num_tenants
                        %Get the minimum and maximum arp of the UEs of the
                        %tenant                        
                        arpmin=max(config.services.ARP); %Initial value
                        arpmax=0;
                        for i=1:BS(n).numUEs(s)
                            if config.services.ARP(BS(n).UElist{s}(i).service)<arpmin
                                arpmin=config.services.ARP(BS(n).UElist{s}(i).service);
                            end;
                            if config.services.ARP(BS(n).UElist{s}(i).service)>arpmax
                                arpmax=config.services.ARP(BS(n).UElist{s}(i).service);
                            end;
                        end;
                        arp=arpmin;
                        while (arp<=arpmax) && (BS(n).congestion_status_per_tenant(s)==0)
                            for i=1:BS(n).numUEs(s)
                                if (config.services.type(BS(n).UElist{s}(i).service)==config.GBR) && (config.services.ARP(BS(n).UElist{s}(i).service)==arp) && (BS(n).UElist{s}(i).activity_state==1)
                                    BS(n).num_assigned_RB_per_tenant(s)=BS(n).num_assigned_RB_per_tenant(s)+BS(n).UElist{s}(i).required_RBs;
                                    BS(n).num_assigned_RB_per_service_per_tenant(s,BS(n).UElist{s}(i).service)=BS(n).num_assigned_RB_per_service_per_tenant(s,BS(n).UElist{s}(i).service)+BS(n).UElist{s}(i).required_RBs;
                                    BS(n).UElist{s}(i).assigned_RBs=BS(n).UElist{s}(i).required_RBs;
                                end;
                            end;
                            %Check congestion status after having assigned
                            %RBs to the users of tenant s
                            %with ARP=arp.
                            if BS(n).num_assigned_RB_per_tenant(s)>BS(n).num_RBs_per_tenant(s)
                                %Congestion!!!
                                %Criterion: we split the excess proportionally between all the
                                %UEs with ARP=arp
                                BS(n).congestion_status_per_tenant(s)=1;
                
                                %Reference: Number of assignedRBs of the GBR UEs with
                                %ARP=arp
                                num_assigned_RBs_GBR_arp=sum(BS(n).num_assigned_RB_per_service_per_tenant(s,(config.services.type(:)==config.GBR)&(config.services.ARP(:)==arp)));
                
                                excess_fraction=(BS(n).num_assigned_RB_per_tenant(s)-BS(n).num_RBs_per_tenant(s))/num_assigned_RBs_GBR_arp;
                                %Meaning: if excess_fraction = 0.1, it means that the excess is
                                %the 10% of all the requirements of the GBR UEs with ARP, so we remove to each UE this
                                %fraction
                            
                                for i=1:BS(n).numUEs(s)
                                    if (config.services.type(BS(n).UElist{s}(i).service)==config.GBR) && (config.services.ARP(BS(n).UElist{s}(i).service)==arp) && (BS(n).UElist{s}(i).activity_state==1)
                                        %Modify new assignment to the UE:
                                        BS(n).UElist{s}(i).assigned_RBs=BS(n).UElist{s}(i).required_RBs*(1-excess_fraction); 
    
                                        %Remove the excess assignment from the total:
                                        BS(n).num_assigned_RB_per_tenant(s)=BS(n).num_assigned_RB_per_tenant(s)-BS(n).UElist{s}(i).required_RBs*excess_fraction;
                                        BS(n).num_assigned_RB_per_service_per_tenant(s,BS(n).UElist{s}(i).service)=BS(n).num_assigned_RB_per_service_per_tenant(s,BS(n).UElist{s}(i).service)-BS(n).UElist{s}(i).required_RBs*excess_fraction;
        
                                    end;
                                end;                    
                            
                            else
                                BS(n).congestion_status_per_tenant(s)=0;
                            end;
                            arp=arp+1;
                        end;    
                        if (BS(n).congestion_status_per_tenant(s)==0)
                            %Assignment of non-GBR:
                            num_available_RBs=BS(n).num_RBs_per_tenant(s)-BS(n).num_assigned_RB_per_tenant(s);
            
                            %1.- Compute the denominator of the resource sharing. It is the
                            %aggregate of (1/priority)^exponent for all the UEs Non-GBR and
                            %active
                            sum_denominator=0;

                            for i=1:BS(n).numUEs(s)
                                if (config.services.type(BS(n).UElist{s}(i).service)==config.NonGBR) && (BS(n).UElist{s}(i).activity_state==1)
                                    sum_denominator=sum_denominator+power(1/config.services.priority(BS(n).UElist{s}(i).service),config.exponent_priority_nonGBR);
                                end;
                            end;

                            %2.- Compute the amount of assigned RBs: 
                            if (sum_denominator>0)
                                for i=1:BS(n).numUEs(s)
                                    if (config.services.type(BS(n).UElist{s}(i).service)==config.NonGBR) && (BS(n).UElist{s}(i).activity_state==1)
                                        BS(n).UElist{s}(i).assigned_RBs=num_available_RBs*power(1/config.services.priority(BS(n).UElist{s}(i).service),config.exponent_priority_nonGBR)/sum_denominator;

                                        BS(n).num_assigned_RB_per_tenant(s)=BS(n).num_assigned_RB_per_tenant(s)+BS(n).UElist{s}(i).assigned_RBs;
                                        BS(n).num_assigned_RB_per_service_per_tenant(s,BS(n).UElist{s}(i).service)=BS(n).num_assigned_RB_per_service_per_tenant(s,BS(n).UElist{s}(i).service)+BS(n).UElist{s}(i).assigned_RBs;
                                    end;
                                end;
                            end;
                        end;
                    end;
                    num_assigned_RB_per_cell_current(n)=sum(BS(n).num_assigned_RB_per_tenant);
                    %If there is congestion for at least one tenant in the
                    %cell, it is assumed the cell is congested.
                    if sum(BS(n).congestion_status_per_tenant)>0  
                        BS(n).congestion_status=1;
                    else
                        BS(n).congestion_status=0;
                    end;
                end;
                
                %Step 3: Check finalisation condition:
     
                error=abs(num_assigned_RB_per_cell_current-num_assigned_RB_per_cell_previous);
     
                error_total=max(error);
                if (error_total<0.01)||(num_iter==100)
                    end_process=1;
         
                    %Compute the final bit rate assigned to each tenant.
                    for n=1:config.num_cells
                        BS(n).bit_rate_assigned_per_tenant=zeros(config.num_tenants,1);
                        BS(n).bit_rate_assigned_per_service_per_tenant=zeros(config.num_tenants,config.num_services);
                        
                       for s=1:config.num_tenants
                            for i=1:BS(n).numUEs(s)
                                BS(n).UElist{s}(i).assigned_bit_rate=BS(n).UElist{s}(i).S*BS(n).UElist{s}(i).assigned_RBs*config.B_RB;
                                BS(n).bit_rate_assigned_per_tenant(s)=BS(n).bit_rate_assigned_per_tenant(s)+BS(n).UElist{s}(i).assigned_bit_rate;
        
                                BS(n).bit_rate_assigned_per_service_per_tenant(s,BS(n).UElist{s}(i).service)=BS(n).bit_rate_assigned_per_service_per_tenant(s,BS(n).UElist{s}(i).service)+BS(n).UElist{s}(i).assigned_bit_rate;
                            end;
                        end;
                    end;
                    if (num_iter==100) && (error_total>0.01)
                        fprintf('Warning: occupation not converged after %d iterations... error: %f\n',num_iter,error_total);
                    end;
                else
                    num_assigned_RB_per_cell_previous=num_assigned_RB_per_cell_current;
                    num_iter=num_iter+1;
                end;                
                
            end;    
            
    end;
end

