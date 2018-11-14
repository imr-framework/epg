function theta_list = get_goldenangle(Ns)

ga = 111.25; % Winklemann S. et al., IEEE TMI 2007
theta_list = zeros(1, Ns);
theta_list(1) = ga;
for t =2: Ns
    theta_list(t) = mod((theta_list(t-1)+ ga), 360); 
end

theta_list = sort(theta_list);
%% Can also include small angle golden angle