function [ Loss ] = prop_model( config, d, no_shadow )
    %PROP_MODEL According to the Urban Micro model (UMi) assuming hexagonal
    %layout
    %config is a struct than contains different parameters of the scenario
    % d is distance in m
    % Loss is the total loss (PL + shadowing) in dB
    % no_shadow: if 1, computation will be without shadowing (used just as
    % a reference value). If 0, computation will be with shadowing (normal
    % operation case)
    
    
    if (d<config.prop_model_params.dmin)
        d=config.prop_model_params.dmin;
    end;
    
    LOS_prob=min(18/d,1)*(1-exp(-d/36))+exp(-d/36);
    
    if (rand()<=LOS_prob)
        %LOS model
        if (d<config.prop_model_params.d_BP)
            Loss=22*log10(d)+28+20*log10(config.prop_model_params.f);
        else
            Loss=40*log10(d)+7.8-18*log10(config.prop_model_params.height_BS-1)-18*log10(config.prop_model_params.height_UE-1)+2*log10(config.prop_model_params.f);
        end;
        if ~no_shadow 
            Loss=Loss+config.prop_model_params.sigma_LOS*randn;
        end;
    else
        %NLOS Model
        Loss=36.7*log10(d)+22.7+26*log10(config.prop_model_params.f);
        if ~no_shadow
            Loss=Loss+config.prop_model_params.sigma_NLOS*randn;
        end;
    end;
  
end

