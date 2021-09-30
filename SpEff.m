function [ S ] = SpEff( config, SINR )

    %Computes the spectral efficiency S in (b/s/Hz) as a function of the SINR in
    %linear.
    %It is based on the model shown in annex A of 3GPP TR 36.942
    
   
    
    
     if SINR<=config.spec_eff_params.SINRmin
        S=0.0;
     elseif SINR<config.spec_eff_params.SINRmax
        S=config.spec_eff_params.alfa*log2(1+SINR);
     else S=config.spec_eff_params.Smax;
     end;

end

