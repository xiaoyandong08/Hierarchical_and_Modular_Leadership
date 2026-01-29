function [all_module_bird_2times,module_bird_2times] = Generate_2times_Modules(One_LFnet,all_module_bird)

all_module_bird_2times = [];
for i = 1 : size(all_module_bird,1)
    if length(all_module_bird{i})>5
        [module_bird_2times{i},module_net_2times{i}] = Generate_Modularity(One_LFnet(all_module_bird{i},all_module_bird{i}));
    else
        module_bird_2times{i} = {[1:1:length(all_module_bird{i})]'};
        module_net_2times{i} = 0;
    end
    for j = 1 : size(module_bird_2times{i},1)
        module_bird_2times{i}{j,2} = all_module_bird{i}(module_bird_2times{i}{j,1});
        all_module_bird_2times = [all_module_bird_2times;module_bird_2times{i}(j,2)];
    end
end
end