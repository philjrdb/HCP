%% HCP Matlab-SPSS unravel
% Puts HCPsum and questionnaire variables 

n_blocks = 6;
pun_trials = 3:6;

RT_lowcutoff = 1; % valid RT lower cut-off in secs (any RT/#qs < low cutoff = invalid RT)
RT_highcutoff = 30; % valid RT upper cut-off in secs (any RT/#qs > high cutoff = invalid RT)

ques_names = {'dass_anx' 'dass_dep' 'bis_sr' 'bis_imp' 'aor_rno' 'aor_rpo'...
  'bisbas_drv' 'bisbas_rew' 'bisbas_fun' 'bisbas_bis'...
  'ipip_e' 'ipip_a' 'ipip_c' 'ipip_n' 'ipip_i'}; % questionnaire names

%% Prep variables
clear SPSS_labels SPSS_mat

start_col = 13;

ag_leng = length(HCP_aggr); 

f_names = fieldnames(HCP_aggr);
q_vars = union(f_names(contains(f_names,'ques')),ques_names);
q_leng = length(q_vars);

vars = fieldnames(HCPsum);
v_leng = length(vars);

SPSS_labels{start_col+v_leng*n_blocks+q_leng} = [];
SPSS_mat = NaN(ag_leng,start_col+v_leng*n_blocks+q_leng);
%cur_col = 0;

%% Demographics/parameters
SPSS_labels{1} = 'Idx';
SPSS_mat(:,1) = [1:ag_leng]';

SPSS_labels{2} = 'json';
u_list = unique({HCP_aggr(:).json},'stable');
for u = 1:length(u_list)
  SPSS_mat(ismember({HCP_aggr(:).json},u_list{u}),2) = u;
end

SPSS_labels{3} = 'Sampling';
u_list = unique({HCP_aggr(:).Sampling},'stable');
for u = 1:length(u_list)
  SPSS_mat(ismember({HCP_aggr(:).Sampling},u_list{u}),3) = u;
end

SPSS_labels{4} = 'Gender'; 
SPSS_mat(ismember({HCP_aggr(:).Gender},'male'),4) = 1;
SPSS_mat(ismember({HCP_aggr(:).Gender},'female'),4) = 2;
SPSS_mat(ismember({HCP_aggr(:).Gender},'other'),4) = 3;

SPSS_labels{5} = 'Age';
SPSS_mat(:,5) = vertcat(HCP_aggr(:).Age);

SPSS_labels{6} = 'Language'; % 1 = includes english, 2 = other
idx = cellfun(@(x) any(x),(regexpi({HCP_aggr(:).Language},'english')));
SPSS_mat(idx,6) = 1;
SPSS_mat(~idx,6) = 2;

SPSS_labels{7} = 'PunPlanet'; % 1 = left, 2 = right
SPSS_mat(ismember({HCP_aggr(:).PunPlanet},'left'),7) = 1;
SPSS_mat(ismember({HCP_aggr(:).PunPlanet},'right'),7) = 2;

SPSS_labels{8} = 'PunShip'; % 1 = TypeI, 2 = TypeII
SPSS_mat(ismember({HCP_aggr(:).PunShip},'TypeI'),8) = 1;
SPSS_mat(ismember({HCP_aggr(:).PunShip},'TypeII'),8) = 2;

SPSS_labels{9} = 'Catch'; % 1 = failed catch questions
SPSS_mat(:,9) = any(horzcat(HCP_aggr(:).Catch),1);

SPSS_labels{10} = 'MinValInfRT';
SPSS_mat(:,10) = cellfun(@(x) min(x,[],'all','omitnan'),{HCP_aggr(:).ValInf_AvgRT});

SPSS_labels{11} = 'MaxValInfRT';
SPSS_mat(:,11) = cellfun(@(x) max(x,[],'all','omitnan'),{HCP_aggr(:).ValInf_AvgRT});

SPSS_labels{12} = 'Exclude';
SPSS_mat(:,12) = any(horzcat(SPSS_mat(:,9),SPSS_mat(:,10)<RT_lowcutoff,SPSS_mat(:,11)>RT_highcutoff),2);

try
  groups = unique({HCP_aggr(:).Group});
  n_groups = length(groups);

  SPSS_labels{13} = ['Group' strjoin(groups,'_')];
  for n = 1:n_groups
    SPSS_mat(ismember(vertcat({HCP_aggr(:).Group})',groups(n)),13) = n;
  end
catch
  warning('No groups');
end

cur_col = start_col; % # previous SPSS columns

%% HCPsum vars
for v = 1:v_leng
  if size(HCPsum.(vars{v}),2) == n_blocks 
    SPSS_mat(:,cur_col+1:cur_col+n_blocks) = HCPsum.(vars{v});

    for p = 1:n_blocks
      SPSS_labels{cur_col+p} = [vars{v} '_B' num2str(p)];
    end
    
    cur_col = cur_col+n_blocks;
  elseif size(HCPsum.(vars{v}),2) == 1 && size(HCPsum.(vars{v}),1) == ag_leng  
    SPSS_mat(:,cur_col+1) = HCPsum.(vars{v});
    
    SPSS_labels{cur_col+1} = vars{v};
    cur_col = cur_col+1;
  end
end

%% Questionnaire vars
for v = 1:q_leng
	var_size = size(HCP_aggr(1).(q_vars{v}));

  if var_size(1)>var_size(2)
    SPSS_mat(:,cur_col+1:cur_col+var_size(1)) = horzcat(HCP_aggr(:).(q_vars{v}))';

  	SPSS_labels{cur_col+var_size(1)} = []; 
    
    for s = 1:var_size(1)
      SPSS_labels{cur_col+s} = [q_vars{v} '_' num2str(s)];
    end

    cur_col = cur_col+var_size(1);

  elseif var_size(2)>var_size(1)
    SPSS_mat(:,cur_col+1:cur_col+var_size(2)) = vertcat(HCP_aggr(:).(q_vars{v}));

    SPSS_labels{cur_col+var_size(2)} = [];
    
    for s = 1:var_size(2)
      SPSS_labels{cur_col+s} = [q_vars{v} '_' num2str(s)];
    end

    cur_col = cur_col+var_size(2);
    
  elseif var_size(2) == 1
    SPSS_mat(:,cur_col+1) = vertcat(HCP_aggr(:).(q_vars{v}));

    SPSS_labels{cur_col+1} = q_vars{v};

    cur_col = cur_col+1;
  end
end

any_idx = any(SPSS_mat(:,start_col+1:end),1);
SPSS_labels = SPSS_labels(:,[true(1,start_col),any_idx]);
SPSS_mat = SPSS_mat(:,[true(1,start_col),any_idx]);
SPSS_mat(isnan(SPSS_mat)) = -999;

T = array2table(SPSS_mat,'VariableNames',SPSS_labels);
writetable(T,[work_folder '\SPSSdata.csv']);

clearvars -except HCP* work_folder *keep SPSS*;

% %% Planet-Attack chain inference estimate (now computed in HCPsum)
% labels = SPSS_labels(cellfun(@(x) ~isempty(x),SPSS_labels));
% for b = 3:5
%   SPSS_labels{cur_col+b} = ['PunP_AtkInfEst_B' num2str(b)];
%   SPSS_mat(:,cur_col+b) = (SPSS_mat(:,...
%     ismember(labels,{['PunPlanet_PunShipInf_B' num2str(b)]})).*SPSS_mat(:,...
%     ismember(labels,{['PunShip_AttackInf_B' num2str(b)]}))/100)+(SPSS_mat(:,...
%     ismember(labels,{['PunPlanet_UnpShipInf_B' num2str(b)]})).*SPSS_mat(:,...
%     ismember(labels,{['UnpShip_AttackInf_B' num2str(b)]}))/100);
% end
% cur_col = cur_col+5;
% for b = 3:5
%   SPSS_labels{cur_col+b} = ['UnpP_AtkInfEst_B' num2str(b)];
%   SPSS_mat(:,cur_col+b) = (SPSS_mat(:,...
%     ismember(labels,['UnpPlanet_PunShipInf_B' num2str(b)])).*SPSS_mat(:,...
%     ismember(labels,['PunShip_AttackInf_B' num2str(b)]))/100)+(SPSS_mat(:,...
%     ismember(labels,['UnpPlanet_UnpShipInf_B' num2str(b)])).*SPSS_mat(:,...
%     ismember(labels,['UnpShip_AttackInf_B' num2str(b)]))/100);
% end
% cur_col = cur_col+5;