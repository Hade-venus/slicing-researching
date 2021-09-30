num=num+1;
i=1;

%Offered load total per tenant (Mb/s)
results(num,i)=stats.offered_load_total(1)/1000; i=i+1;
results(num,i)=stats.offered_load_total(2)/1000; i=i+1;

%Average number of RBs per service (average for all cells)
results(num,i)=stats.avg_RB_occupation_per_service(1); i=i+1;
results(num,i)=stats.avg_RB_occupation_per_service(2); i=i+1;
results(num,i)=stats.avg_RB_occupation_per_service(3); i=i+1;
results(num,i)=stats.avg_RB_occupation_per_service(4); i=i+1;
results(num,i)=stats.avg_RB_occupation_per_service(6); i=i+1;
results(num,i)=stats.avg_RB_occupation_per_service(7); i=i+1;
results(num,i)=stats.avg_RB_occupation_per_service(8); i=i+1;
results(num,i)=stats.avg_RB_occupation; i=i+1;
results(num,i)=stats.avg_RB_occupation_GBR; i=i+1;

%Average bit rate (throughput) per tenant in Mb/s (aggregate for all cells)
results(num,i)=stats.avg_bit_rate_per_tenant_total(1)/1000; i=i+1;
results(num,i)=stats.avg_bit_rate_per_tenant_total(2)/1000;  i=i+1;
results(num,i)=(stats.avg_bit_rate_per_tenant_total(1)+stats.avg_bit_rate_per_tenant_total(2))/1000; i=i+1;
results(num,i)=stats.avg_bit_rate_per_tenant_total_GBR(1)/1000; i=i+1;
results(num,i)=stats.avg_bit_rate_per_tenant_total_GBR(2)/1000;  i=i+1;
results(num,i)=(stats.avg_bit_rate_per_tenant_total_GBR(1)+stats.avg_bit_rate_per_tenant_total_GBR(2))/1000; i=i+1;

%Average bit rate (throughput) per service in Mb/s (aggregate for all
%cells)
results(num,i)=stats.avg_bit_rate_per_service_per_tenant_total(1,1)/1000; i=i+1;
results(num,i)=stats.avg_bit_rate_per_service_per_tenant_total(1,2)/1000; i=i+1;
results(num,i)=stats.avg_bit_rate_per_service_per_tenant_total(1,3)/1000; i=i+1;
results(num,i)=stats.avg_bit_rate_per_service_per_tenant_total(1,4)/1000; i=i+1;
results(num,i)=stats.avg_bit_rate_per_service_per_tenant_total(2,6)/1000; i=i+1;
results(num,i)=stats.avg_bit_rate_per_service_per_tenant_total(2,7)/1000; i=i+1;
results(num,i)=stats.avg_bit_rate_per_service_per_tenant_total(2,8)/1000; i=i+1;

%Average bit rate per user per service in Mb/s (average for all the UEs in
%all the cells)
results(num,i)=stats.avg_Rb_per_user_total_per_service_per_tenant(1,1)/1000; i=i+1;
results(num,i)=stats.avg_Rb_per_user_total_per_service_per_tenant(1,2)/1000; i=i+1;
results(num,i)=stats.avg_Rb_per_user_total_per_service_per_tenant(1,3)/1000; i=i+1;
results(num,i)=stats.avg_Rb_per_user_total_per_service_per_tenant(1,4)/1000; i=i+1;
results(num,i)=stats.avg_Rb_per_user_total_per_service_per_tenant(2,6)/1000; i=i+1;
results(num,i)=stats.avg_Rb_per_user_total_per_service_per_tenant(2,7)/1000; i=i+1;
results(num,i)=stats.avg_Rb_per_user_total_per_service_per_tenant(2,8)/1000; i=i+1;

%Probability of degradation for GBR services in % (for all the UEs in all
%the cells)
results(num,i)=stats.prob_degraded_UE_total_per_service_per_tenant(1,1)*100; i=i+1;
results(num,i)=stats.prob_degraded_UE_total_per_service_per_tenant(1,3)*100; i=i+1;
results(num,i)=stats.prob_degraded_UE_total_per_service_per_tenant(2,6)*100; i=i+1;
results(num,i)=stats.prob_degraded_UE_total_per_service_per_tenant(2,7)*100; i=i+1;

%Percentage of degradation in the Rb of GBR services (average for all the
%UEs in all the cells)
results(num,i)=stats.avg_perc_Rb_degradation_total_per_service_per_tenant(1,1); i=i+1;
results(num,i)=stats.avg_perc_Rb_degradation_total_per_service_per_tenant(1,3); i=i+1;
results(num,i)=stats.avg_perc_Rb_degradation_total_per_service_per_tenant(2,6); i=i+1;
results(num,i)=stats.avg_perc_Rb_degradation_total_per_service_per_tenant(2,7); i=i+1;

%Percentile 5 of the bit rate per user for NonGBR services in Mb/s (for all
%UEs in all cells)
results(num,i)=stats.perc5_Rb_per_user_total_per_service_per_tenant(1,2)/1000; i=i+1;
results(num,i)=stats.perc5_Rb_per_user_total_per_service_per_tenant(1,4)/1000; i=i+1;
results(num,i)=stats.perc5_Rb_per_user_total_per_service_per_tenant(2,8)/1000; i=i+1;

%Percentile 95 of the bit rate per user for NonGBR services in Mb/s (for all
%UEs in all cells)
results(num,i)=stats.perc95_Rb_per_user_total_per_service_per_tenant(1,2)/1000; i=i+1;
results(num,i)=stats.perc95_Rb_per_user_total_per_service_per_tenant(1,4)/1000; i=i+1;
results(num,i)=stats.perc95_Rb_per_user_total_per_service_per_tenant(2,8)/1000; i=i+1;

%Blocking probability per tenant in all the cells (%)
results(num,i)=stats.blocking_prob_total(1)*100; i=i+1;
results(num,i)=stats.blocking_prob_total(2)*100; i=i+1;

%Blocking probability for GBR services in all the cells (%)
results(num,i)=stats.blocking_prob_total_per_service_per_tenant(1,1)*100; i=i+1;
results(num,i)=stats.blocking_prob_total_per_service_per_tenant(1,3)*100; i=i+1;
results(num,i)=stats.blocking_prob_total_per_service_per_tenant(2,6)*100; i=i+1;
results(num,i)=stats.blocking_prob_total_per_service_per_tenant(2,7)*100; i=i+1;

%Congestion probability (%) in cell 1. Total and per tenant.
results(num,i)=stats.congestion_prob_per_cell*100; i=i+1;
results(num,i)=stats.congestion_prob_per_cell_per_tenant(1)*100; i=i+1;
results(num,i)=stats.congestion_prob_per_cell_per_tenant(2)*100; i=i+1;

%Average bit rate (Mb/s) per RB in cell 1. Average spectral efficiency in
%cell 1.
results(num,i)=stats.avg_bit_rate_per_RB/1000; i=i+1;
results(num,i)=stats.avg_sp_efficiency; i=i+1;



