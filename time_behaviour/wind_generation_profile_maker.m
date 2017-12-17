%config
T = 10; %seconds
T_import = 10*60; %seconds

wind_data = load('time_behaviour\wind_generation_data_sotavento_galicia.mat');

production_normalised = wind_data.EnergykWh/max(wind_data.EnergykWh);

time_step_factor = T_import/T;

production = kron(production_normalised, ones(time_step_factor,1));

p_ref_lower_limit_base = [[   0;    0; -0.2; -0.56]*ones(1,400), [   0;    0; -0.2  ; -0.56]*ones(1,length(production))];

p_ref_upper_limit_base = [[ 1.0;  0.6; -0.2; -0.56]*ones(1,400), [ 1.0; 0; -0.2; -0.56]*ones(1,length(production))];

p_ref_upper_limit_base(2,401:end) = production;
p_ref_lower_limit_base(2,401:end) = production;

save('p_ref_limits_wind.mat','p_ref_upper_limit_base','p_ref_lower_limit_base');