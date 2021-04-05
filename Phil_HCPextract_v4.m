%% HumanCondPun extraction
% Take json files within "Raw" folder (work_folder\Raw\jatos###.txt) - preload Raw folder with json files
% and extracts into HCP_aggr variable in workspace

redo_extract = 1; % if 1, will completely redo extraction
pun_trials = 3:5;

end_trial_wait = 1000;
block_duration = 180*1000+end_trial_wait;

HCP_parameters.signal_time = 2000;
HCP_parameters.show_ship_delay = 1500;
HCP_parameters.ship_duration = 6000;
HCP_parameters.ship_attack_time = 6000;
HCP_parameters.feedback_duration = 2500;
HCP_parameters.shield_charging_time = 3000;
HCP_parameters.shield_cost = -50; % point loss (i.e. -#)

HCP_parameters.BinLabels = {'ITI'; 'PunS_PreShield'; 'PunS_ShieldAv'; ...
        'PunS_Shielded'; 'PunS_NoShield'; 'PunS_ShieldedOutcome'; 'PunS_UnshieldedOutcome'; ...
        'UnpS_PreShield'; 'UnpS_ShieldAv'; 'UnpS_Shielded'; 'UnpS_NoShield'; ...
        'UnpS_ShieldedOutcome'; 'UnpS_UnshieldedOutcome'}; 
      % bins for HCP_windowbins.m 

phase_names = {'phase1' 'Phase2'}; % {'phase1' 'Phase2'} 2nd entry = pun block
ques_names = {'ques_dass' 'ques_bis' 'ques_aor' 'ques_bisbas' 'ques_ipip'}; % questionnaire names
ques_checks = {[] [] [8;0] [10;3] []}; % [index;correct answer] for catch questions
RT_check_grid = [3 3 6 6 6; 2 2 8 8 8; 2 2 8 8 8; 0 0 2 2 2; 0 0 2 2 2]; % # post-block val/inf questions
RT_lowcutoff = 1; % valid RT lower cut-off in secs (any RT/#qs < low cutoff = invalid RT)
RT_highcutoff = 30; % valid RT upper cut-off in secs (any RT/#qs > high cutoff = invalid RT)

work_folder = ...
   ['M:\Gavan McNally'...
   's Lab\Staff folders\Philip\Post-sub\Experiments\Collaborations\Human Conditioned Punishment\Exp 2'];

ev_ylabel = {'PunShip' 'UnpShip' 'Shield' 'PunClick' 'UnpClick' 'ShieldClick' 'OtherClick' '+100' 'Loss'};
ev_ytick = [12 11 10 8 7 6 5 2 1];
ev_binning_tolerance = 100; %ms

%% House-keeping
clc
close all 
fclose('all');
start_tic_keep = tic;

n_bins = length(HCP_parameters.BinLabels);
[~,ytick_idx] = sort(ev_ytick);

% Make folders, save scripts/functions to Scripts folder
mkdir(work_folder, 'Figs');
mkdir([work_folder '\Figs'],'Click location');
mkdir([work_folder '\Figs'],'Summary');
mkdir([work_folder '\Figs'],'Timecourse');
mkdir([work_folder '\Scripts'],date);
copy_file = {'Phil_HCPextract_v4.m' 'AvBbinning.m' 'bin_rate2.m' 'HCP_sumfig2.m' ...
  'HCP_windowbins.m' 'weighted_rate.m'};
for c = 1:length(copy_file)
   copyfile(copy_file{c},[work_folder '\Scripts\' date])
end

% Block structure fieldnames
cell_vars = {'NormClickLoc' 'PunClicks' 'UnpClicks' 'ShieldClicks' 'OtherClicks' 'PunShips' ...
  'UnpShips' 'PunShields' 'UnpShields' 'Attack' 'ShieldCost' 'PunRew' 'UnpRew' ...
  'BinPunRate' 'BinUnpRate' 'BinShieldRate' 'BinTime'};
d_vars = {'EndPoints' 'End' 'RewVal' 'PunPlanetVal' 'UnpPlanetVal' 'PunShipVal' 'UnpShipVal' 'AttackVal' ...
  'PunPlanet_RewInf' 'PunPlanet_RewInfconf' 'UnpPlanet_RewInf' 'UnpPlanet_RewInfconf' ...
  'PunPlanet_PunShipInf' 'PunPlanet_PunShipInfconf' 'PunPlanet_UnpShipInf' 'PunPlanet_UnpShipInfconf' ...
  'PunPlanet_AttackInf' 'PunPlanet_AttackInfconf' 'UnpPlanet_PunShipInf' 'UnpPlanet_PunShipInfconf' ...
  'UnpPlanet_UnpShipInf' 'UnpPlanet_UnpShipInfconf' 'UnpPlanet_AttackInf' 'UnpPlanet_AttackInfconf' ...
  'PunShip_AttackInf' 'PunShip_AttackInfconf' 'UnpShip_AttackInf' 'UnpShip_AttackInfconf' ...
  'ValCheckRT' 'InfCheck1RT' 'InfCheck2RT' 'InfCheckSh1RT' 'InfCheckSh2RT'};

diaryfile = [work_folder '\Scripts\' date '\HCP extract (' date ').txt'];
diary(diaryfile);
 
%% Searches work_folder for new data - load 
if redo_extract == 1
  fprintf('Redoing extraction. '); 
  clearvars HCP_aggr
  HCP_aggr(1).json = [];
  cur_row = 0;
else % Searches work_folder for new data - load 
  if exist('HCP_aggr','var') == 1
    fprintf('Aggregate already loaded. '); 
    cur_row = length(HCP_aggr); 
    if isempty(HCP_aggr(cur_row).json)
      error('HCP_aggr?!');
    end
  elseif exist([work_folder '\HCP.mat'],'file')
    load([work_folder '\HCP.mat']);
    fprintf('Previous aggregate loaded. ');
  else
    fprintf('Aggregate not loaded. '); 
    HCP_aggr(1).json = [];
    cur_row = 0;
  end
end

%Set block list to blocks not within HCP_meta
datasets = dir([work_folder '\Raw']);
datasets = datasets(endsWith({datasets.name},'.txt'));
dataset_new = {datasets.name}';
dataset_new = cellfun(@(x) x(1:end-4), dataset_new, 'UniformOutput', false);
if cur_row > 0
  dataset_new = setdiff(dataset_new, {HCP_aggr(:).json});
end

if ~isempty(dataset_new)
  fprintf('Add: \n');
  disp(dataset_new)

  %Preallocate extra cells (doesn't affect prev aggr cells)
  n_datasets = length(dataset_new);

  for d = 1:n_datasets
    data_name = dataset_new{d};
    fprintf(['Extracting ' data_name '. ']);

    json_text = importdata([work_folder '\Raw\' data_name '.txt']);

    n_j = length(json_text);
    fprintf([num2str(n_j) ' set(s) found\n']);

    if n_j > 0
      HCP_aggr(cur_row+n_j).raw{1,3} = [];
%       HCP_meta{cur_row+n_j,4} = [];

      for s = 1:n_j
        cur_row = cur_row+1;
        
        data = jsondecode(json_text{s});
        data_leng = length(data);
        HCP_aggr(cur_row).raw{data_leng,3} = [];
        HCP_aggr(cur_row).raw(:,3) = data;
        
        HCP_aggr(cur_row).json = data_name;
        HCP_aggr(cur_row).Set = s;
        HCP_aggr(cur_row).Sampling = HCP_aggr(cur_row).raw{1,3}.sample;

        for p = 1:data_leng
           if isfield(data{p},'phase')
              HCP_aggr(cur_row).raw{p,1} = data{p}.phase;
           end
           if isfield(data{p},'block_number')
              HCP_aggr(cur_row).raw{p,2} = data{p}.block_number;
           end 
        end
        
        fprintf('Extracted. ')
        toc(start_tic_keep);
      end
    end
  end
else
  fprintf('None new.\n');
end

fprintf('\n');

%% Process data into HCP_aggr structure
fprintf('Processing data...\n');
for r = 1:length(HCP_aggr)
  fprintf(['\n - aggr' num2str(r) '-\n']);
     
  % Demographics
  is_idx = find(cellfun(@(x) ~isempty(x),HCP_aggr(r).raw(:,1)));
  
  p = is_idx(ismember(HCP_aggr(r).raw(is_idx,1),'demographics'));
  [~,g_eIdx] = regexp(HCP_aggr(r).raw{p,3}.responses,'"gender":"');
  [a_sIdx,a_eIdx] = regexp(HCP_aggr(r).raw{p,3}.responses,'","age":"');
  [l_sIdx,l_eIdx] = regexp(HCP_aggr(r).raw{p,3}.responses,'","language":"');
  
  HCP_aggr(r).Gender = HCP_aggr(r).raw{p,3}.responses(g_eIdx+1:a_sIdx-1);
  HCP_aggr(r).Age = str2num(HCP_aggr(r).raw{p,3}.responses(a_eIdx+1:l_sIdx-1));
  if isempty(HCP_aggr(r).Age)
    HCP_aggr(r).Age = NaN;
  end
  HCP_aggr(r).Language = HCP_aggr(r).raw{p,3}.responses(l_eIdx+1:end-2);
  
  fprintf([HCP_aggr(r).Gender ',' num2str(HCP_aggr(r).Age) ',' HCP_aggr(r).Language '\n']);
  
  %% Prep variables 
  phase_idx = is_idx(ismember(HCP_aggr(r).raw(is_idx,1),phase_names));
  n_phases = length(phase_idx);
  
  punP = HCP_aggr(r).raw{phase_idx(1),3}.pun_planet_side;
  if punP == '1'
     unpP = '0';
     punP_n = 1;
     unpP_n = 0;
     HCP_aggr(r).PunPlanet = 'right';
     fprintf('PunPlanet = right |');
  else
     unpP = '1';
     punP_n = 0;
     unpP_n = 1;
     HCP_aggr(r).PunPlanet = 'left';
     fprintf('PunPlanet = left | ');
  end
  [~,idx] = regexp(HCP_aggr(r).raw{phase_idx(end),3}.pun_ship,'ship');
  punS = HCP_aggr(r).raw{phase_idx(end),3}.pun_ship(idx+1);
  if punS == '1'
    unpS = '2';
     HCP_aggr(r).PunShip = 'TypeI';
    fprintf('PunShip = Type I\n');
  elseif punS == '2'
    unpS = '1';
     HCP_aggr(r).PunShip = 'TypeII';
    fprintf('PunShip = Type II\n');
  else
    error('Ship types?');
  end

  click_loc{n_phases} = [];
  pun_clicks{n_phases} = [];
  unp_clicks{n_phases} = [];
  shield_clicks{n_phases} = [];
  other_clicks{n_phases} = [];
  pun_ships{n_phases} = [];
  unp_ships{n_phases} = [];
  
  % Prep saved variable structure
  for v = 1:length(cell_vars)
    HCP_aggr(r).Block.(cell_vars{v}){n_phases} = [];
  end
  for v = 1:length(d_vars)
    HCP_aggr(r).Block.(d_vars{v}) = NaN(1,n_phases);
  end
  
  click_fig = figure;
  normclick_fig = figure;
  t_fig = figure; hold on
  ylim([min(ev_ytick)*.5 max(ev_ytick)+min(ev_ytick)*.5]);
  yticks(ev_ytick(ytick_idx));
  yticklabels(ev_ylabel(ytick_idx));
  xlabel('Seconds');
  prev_times = 0;

  %% Loop through phases
  for p = 1:n_phases
    fprintf(['Block' num2str(p) '\n']);
    
    click_loc{p} = HCP_aggr(r).raw{phase_idx(p),3}.all_clicks.loc;
    if ~isempty(click_loc{p})
      view_size = HCP_aggr(r).raw{phase_idx(p),3}.viewport_size./2;
      if HCP_aggr(r).raw{phase_idx(p),3}.viewport_size(1)>HCP_aggr(r).raw{phase_idx(p),3}.viewport_size(2)
        HCP_aggr(r).Block.NormClickLoc{p} = click_loc{p};
          HCP_aggr(r).Block.NormClickLoc{p}(:,1) = (HCP_aggr(r).Block.NormClickLoc{p}(:,1)-view_size(1))./view_size(1);
          HCP_aggr(r).Block.NormClickLoc{p}(:,2) = -(HCP_aggr(r).Block.NormClickLoc{p}(:,2)-view_size(2))./(view_size(2));
      else
        HCP_aggr(r).Block.NormClickLoc{p} = fliplr(click_loc{p});
        HCP_aggr(r).Block.NormClickLoc{p}(:,1) = (HCP_aggr(r).Block.NormClickLoc{p}(:,1)-view_size(2))./view_size(2);
        HCP_aggr(r).Block.NormClickLoc{p}(:,2) = -(HCP_aggr(r).Block.NormClickLoc{p}(:,2)-view_size(1))./(view_size(1));
      end
    end % skip if no clicks
    
    clicks = HCP_aggr(r).raw{phase_idx(p),3}.all_clicks.element;
    if ~isempty(clicks)
      clicks(cellfun(@(x) isempty(x),HCP_aggr(r).raw{phase_idx(p),3}.all_clicks.element)) = {'none'};
    end
    
    %% Clicks (timestamps [per type])
    pun_clicks{p} = ismember(clicks,['planet-' punP]);
    HCP_aggr(r).Block.PunClicks{p} = ...
      HCP_aggr(r).raw{phase_idx(p),3}.all_clicks.timestamp(pun_clicks{p})/1000;
    unp_clicks{p} = ismember(clicks,['planet-' unpP]);
    HCP_aggr(r).Block.UnpClicks{p} = ...
      HCP_aggr(r).raw{phase_idx(p),3}.all_clicks.timestamp(unp_clicks{p})/1000;    
    shield_clicks{p} = ismember(clicks,'ship-shield-button');
    HCP_aggr(r).Block.ShieldClicks{p} = ...
      HCP_aggr(r).raw{phase_idx(p),3}.all_clicks.timestamp(shield_clicks{p})/1000;
    other_clicks{p} = ismember(clicks,'none');
    HCP_aggr(r).Block.OtherClicks{p} = ...
      HCP_aggr(r).raw{phase_idx(p),3}.all_clicks.timestamp(other_clicks{p})/1000;
    
    if ismember(p,pun_trials)
      %% Ships [pun v unp] (col1 = start, col2 = outcome (disappear = col2 + feedback_duration))
      pun_ships{p} = ismember(HCP_aggr(r).raw{phase_idx(p),3}.ships.type,punP);
      HCP_aggr(r).Block.PunShips{p} = zeros(sum(pun_ships{p}),2);
      HCP_aggr(r).Block.PunShips{p}(:,1) = ...
        HCP_aggr(r).raw{phase_idx(p),3}.ships.time_appear(pun_ships{p})/1000;
      HCP_aggr(r).Block.PunShips{p}(:,2) = ...
        (HCP_aggr(r).raw{phase_idx(p),3}.ships.time_outcome(pun_ships{p}))/1000;
      
      unp_ships{p} = ismember(HCP_aggr(r).raw{phase_idx(p),3}.ships.type,unpP);
      HCP_aggr(r).Block.UnpShips{p} = zeros(sum(unp_ships{p}),2);
      HCP_aggr(r).Block.UnpShips{p}(:,1) = ...
        HCP_aggr(r).raw{phase_idx(p),3}.ships.time_appear(unp_ships{p})/1000;
      HCP_aggr(r).Block.UnpShips{p}(:,2) = ...
        (HCP_aggr(r).raw{phase_idx(p),3}.ships.time_outcome(unp_ships{p}))/1000;

      %% Shields (col1 = start, col2 = rt activated, col3 = outcome (disappear = col3 + feedback_duration))
      idx = HCP_aggr(r).raw{phase_idx(p),3}.ships.shield_available(pun_ships{p});
      HCP_aggr(r).Block.PunShields{p} = NaN(sum(pun_ships{p}),3);
      HCP_aggr(r).Block.PunShields{p}(idx,1) = ...
        HCP_aggr(r).Block.PunShips{p}(idx,1)+HCP_parameters.shield_charging_time/1000;
      HCP_aggr(r).Block.PunShields{p}(idx,2) = ...
        (HCP_aggr(r).raw{phase_idx(p),3}.ships.rt_shield_activated(idx))/1000;
      HCP_aggr(r).Block.PunShields{p}(idx,3) = ...
        HCP_aggr(r).Block.PunShips{p}(idx,2);

      idx = HCP_aggr(r).raw{phase_idx(p),3}.ships.shield_available(unp_ships{p});
      HCP_aggr(r).Block.UnpShields{p} = NaN(sum(unp_ships{p}),3);
      HCP_aggr(r).Block.UnpShields{p}(idx,1) = ...
        HCP_aggr(r).Block.UnpShips{p}(idx,1)+HCP_parameters.shield_charging_time/1000;
      HCP_aggr(r).Block.UnpShields{p}(idx,2) = ...
        (HCP_aggr(r).raw{phase_idx(p),3}.ships.rt_shield_activated(idx))/1000;
      HCP_aggr(r).Block.UnpShields{p}(idx,3) = ...
        HCP_aggr(r).Block.UnpShips{p}(idx,2);
    end
    
    %% Point Outcomes (col1 = amount, col2 = time [per cause])
    HCP_aggr(r).Block.EndPoints(p) = HCP_aggr(r).raw{phase_idx(p),3}.points_total;
    
    loss_idx = HCP_aggr(r).raw{phase_idx(p),3}.all_outcomes.outcome < 0;  
    if any(loss_idx)
      % Attack
      idx = HCP_aggr(r).raw{phase_idx(p),3}.ships.outcome < 0;
      HCP_aggr(r).Block.Attack{p} = zeros(sum(idx),2);
      HCP_aggr(r).Block.Attack{p}(:,1) = HCP_aggr(r).raw{phase_idx(p),3}.ships.outcome(idx);
      HCP_aggr(r).Block.Attack{p}(:,2) = HCP_aggr(r).raw{phase_idx(p),3}.ships.time_outcome(idx);
      
      % Shield cost
      if HCP_parameters.shield_cost < 0
        points = HCP_aggr(r).raw{phase_idx(p),3}.all_outcomes.outcome(loss_idx);
        times = HCP_aggr(r).raw{phase_idx(p),3}.all_outcomes.time_outcome(loss_idx)/1000;
        idx = points == HCP_parameters.shield_cost;
        HCP_aggr(r).Block.ShieldCost{p} = zeros(sum(idx),2);
        HCP_aggr(r).Block.ShieldCost{p}(:,1) = points(idx);
        HCP_aggr(r).Block.ShieldCost{p}(:,2) = times(idx);
      else
        HCP_aggr(r).Block.ShieldCost{p} = [];
      end     
    else
      HCP_aggr(r).Block.ShieldCost{p} = [];
      HCP_aggr(r).Block.Attack{p} = [];
    end
    
    points = HCP_aggr(r).raw{phase_idx(p),3}.planets.outcome;
    times = HCP_aggr(r).raw{phase_idx(p),3}.planets.time_outcome/1000;
    
    [A_idx,B_idx] = AvBbinning(times,...
       HCP_aggr(r).Block.PunClicks{p},HCP_aggr(r).Block.UnpClicks{p},...
       HCP_parameters.signal_time/1000,ev_binning_tolerance/1000);
    
    HCP_aggr(r).Block.PunRew{p} = zeros(sum(A_idx),2);
    HCP_aggr(r).Block.PunRew{p}(:,1) = points(A_idx);
    HCP_aggr(r).Block.PunRew{p}(:,2) = times(A_idx);
    HCP_aggr(r).Block.UnpRew{p} = zeros(sum(B_idx),2);
    HCP_aggr(r).Block.UnpRew{p}(:,1) = points(B_idx);
    HCP_aggr(r).Block.UnpRew{p}(:,2) = times(B_idx);
    
%    time(1) = min([min(HCP_aggr(r).raw{phase_idx(p),3}.all_outcomes.time_outcome) ...
%      min(HCP_aggr(r).raw{phase_idx(p),3}.all_clicks.timestamp)])/1000;
    HCP_aggr(r).Block.End(p) = max([block_duration ...
      max(HCP_aggr(r).raw{phase_idx(p),3}.all_outcomes.time_outcome)+HCP_parameters.feedback_duration ...
      max(HCP_aggr(r).raw{phase_idx(p),3}.all_clicks.timestamp)])/1000;
    
    
    %% Value & inference checks
    HCP_aggr(r).Block.ValCheckRT(p) = HCP_aggr(r).raw{phase_idx(p)+1,3}.rt/1000;
    HCP_aggr(r).Block.InfCheck1RT(p) = HCP_aggr(r).raw{phase_idx(p)+2,3}.rt/1000;
    HCP_aggr(r).Block.InfCheck2RT(p) = HCP_aggr(r).raw{phase_idx(p)+3,3}.rt/1000;
    
    HCP_aggr(r).Block.RewVal(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+1,3}.val_1);
    
    % Planet values
    if punP == '0'
      HCP_aggr(r).Block.PunPlanetVal(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+1,3}.val_2);
      HCP_aggr(r).Block.UnpPlanetVal(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+1,3}.val_3);
      
      HCP_aggr(r).Block.PunPlanet_RewInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.rate_1);
      HCP_aggr(r).Block.PunPlanet_RewInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.conf_1);
      HCP_aggr(r).Block.UnpPlanet_RewInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.rate_1);   
      HCP_aggr(r).Block.UnpPlanet_RewInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.conf_1);   
    elseif punP == '1'
      HCP_aggr(r).Block.PunPlanetVal(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+1,3}.val_3);
      HCP_aggr(r).Block.UnpPlanetVal(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+1,3}.val_2);
      
      HCP_aggr(r).Block.PunPlanet_RewInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.rate_1);
      HCP_aggr(r).Block.PunPlanet_RewInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.conf_1);
      HCP_aggr(r).Block.UnpPlanet_RewInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.rate_1);   
      HCP_aggr(r).Block.UnpPlanet_RewInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.conf_1);  
    else
      error('Huh planet?')
    end

    if ismember(p,pun_trials) % if punished phase
      HCP_aggr(r).Block.InfCheckSh1RT(p) = HCP_aggr(r).raw{phase_idx(p)+4,3}.rt/1000;
      HCP_aggr(r).Block.InfCheckSh2RT(p) = HCP_aggr(r).raw{phase_idx(p)+5,3}.rt/1000;
      
      HCP_aggr(r).Block.AttackVal(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+1,3}.val_6);
      
      if punS == '1'
        HCP_aggr(r).Block.PunShipVal(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+1,3}.val_4);
        HCP_aggr(r).Block.UnpShipVal(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+1,3}.val_5);
        
        HCP_aggr(r).Block.PunShip_AttackInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+4,3}.rate_1);
        HCP_aggr(r).Block.PunShip_AttackInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+4,3}.conf_1);
        HCP_aggr(r).Block.UnpShip_AttackInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+5,3}.rate_1);
        HCP_aggr(r).Block.UnpShip_AttackInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+5,3}.conf_1);
        
        if punP == '0'
          HCP_aggr(r).Block.PunPlanet_PunShipInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.rate_2);
          HCP_aggr(r).Block.PunPlanet_PunShipInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.conf_2);
          HCP_aggr(r).Block.PunPlanet_UnpShipInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.rate_3);
          HCP_aggr(r).Block.PunPlanet_UnpShipInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.conf_3);
          HCP_aggr(r).Block.PunPlanet_AttackInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.rate_4);
          HCP_aggr(r).Block.PunPlanet_AttackInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.conf_4);
          
          HCP_aggr(r).Block.UnpPlanet_PunShipInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.rate_2);
          HCP_aggr(r).Block.UnpPlanet_PunShipInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.conf_2);
          HCP_aggr(r).Block.UnpPlanet_UnpShipInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.rate_3);
          HCP_aggr(r).Block.UnpPlanet_UnpShipInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.conf_3);
          HCP_aggr(r).Block.UnpPlanet_AttackInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.rate_4);
          HCP_aggr(r).Block.UnpPlanet_AttackInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.conf_4);
        elseif punP == '1'
          HCP_aggr(r).Block.PunPlanet_PunShipInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.rate_2);
          HCP_aggr(r).Block.PunPlanet_PunShipInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.conf_2);
          HCP_aggr(r).Block.PunPlanet_UnpShipInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.rate_3);
          HCP_aggr(r).Block.PunPlanet_UnpShipInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.conf_3);
          HCP_aggr(r).Block.PunPlanet_AttackInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.rate_4);
          HCP_aggr(r).Block.PunPlanet_AttackInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.conf_4);
          
          HCP_aggr(r).Block.UnpPlanet_PunShipInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.rate_2);
          HCP_aggr(r).Block.UnpPlanet_PunShipInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.conf_2);
          HCP_aggr(r).Block.UnpPlanet_UnpShipInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.rate_3);
          HCP_aggr(r).Block.UnpPlanet_UnpShipInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.conf_3);
          HCP_aggr(r).Block.UnpPlanet_AttackInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.rate_4);
          HCP_aggr(r).Block.UnpPlanet_AttackInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.conf_4);
        else
          error('Huh planet?')
        end
        
      elseif punS == '2'
        HCP_aggr(r).Block.PunShipVal(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+1,3}.val_5);
        HCP_aggr(r).Block.UnpShipVal(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+1,3}.val_4);
        
        HCP_aggr(r).Block.PunShip_AttackInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+5,3}.rate_1);
        HCP_aggr(r).Block.PunShip_AttackInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+5,3}.conf_1);
        HCP_aggr(r).Block.UnpShip_AttackInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+4,3}.rate_1);
        HCP_aggr(r).Block.UnpShip_AttackInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+4,3}.conf_1);
      
        if punP == '0'
          HCP_aggr(r).Block.PunPlanet_PunShipInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.rate_3);
          HCP_aggr(r).Block.PunPlanet_PunShipInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.conf_3);
          HCP_aggr(r).Block.PunPlanet_UnpShipInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.rate_2);
          HCP_aggr(r).Block.PunPlanet_UnpShipInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.conf_2);
          HCP_aggr(r).Block.PunPlanet_AttackInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.rate_4);
          HCP_aggr(r).Block.PunPlanet_AttackInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.conf_4);
          
          HCP_aggr(r).Block.UnpPlanet_PunShipInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.rate_3);
          HCP_aggr(r).Block.UnpPlanet_PunShipInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.conf_3);
          HCP_aggr(r).Block.UnpPlanet_UnpShipInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.rate_2);
          HCP_aggr(r).Block.UnpPlanet_UnpShipInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.conf_2);
          HCP_aggr(r).Block.UnpPlanet_AttackInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.rate_4);
          HCP_aggr(r).Block.UnpPlanet_AttackInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.conf_4);
        elseif punP == '1'        
          HCP_aggr(r).Block.PunPlanet_PunShipInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.rate_3);
          HCP_aggr(r).Block.PunPlanet_PunShipInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.conf_3);
          HCP_aggr(r).Block.PunPlanet_UnpShipInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.rate_2);
          HCP_aggr(r).Block.PunPlanet_UnpShipInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.conf_2);
          HCP_aggr(r).Block.PunPlanet_AttackInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.rate_4);
          HCP_aggr(r).Block.PunPlanet_AttackInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+3,3}.conf_4);
          
          HCP_aggr(r).Block.UnpPlanet_PunShipInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.rate_3);
          HCP_aggr(r).Block.UnpPlanet_PunShipInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.conf_3);
          HCP_aggr(r).Block.UnpPlanet_UnpShipInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.rate_2);
          HCP_aggr(r).Block.UnpPlanet_UnpShipInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.conf_2);
          HCP_aggr(r).Block.UnpPlanet_AttackInf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.rate_4);
          HCP_aggr(r).Block.UnpPlanet_AttackInfconf(p) = str2double(HCP_aggr(r).raw{phase_idx(p)+2,3}.conf_4);
        else
          error('Huh planet?')
        end
      else
        error('Huh ship?')
      end
      
      %% Click rate binning
      % Ship/shield bins (all else = ITI) 
      % {1-6} PunShip: PreShield, ShieldAvail, Shielded, NoShield, ShieldedOutcome, UnshieldedOutcome     
      % {7-12} UnpShip: PreShield, ShieldAvail, Shielded, NoShield, ShieldedOutcome, UnshieldedOutcome 
      windows = HCP_windowbins(HCP_aggr(r).Block.PunShips{p},HCP_aggr(r).Block.PunShields{p},...
        HCP_aggr(r).Block.UnpShips{p},HCP_aggr(r).Block.UnpShields{p},...
        HCP_parameters.shield_charging_time/1000,HCP_parameters.feedback_duration/1000);
           
      % Set BinTime defaults
      ITI_time = HCP_aggr(r).Block.End(p);
      times = zeros(length(windows),1);
      
      % Bin PunPlanet clicks
      if ~isempty(HCP_aggr(r).Block.PunClicks{p})
        [ITI_rate, ITI_time, rates, times] = bin_rate2(HCP_aggr(r).Block.PunClicks{p},...
          HCP_aggr(r).Block.End(p), windows);
        HCP_aggr(r).Block.BinPunRate{p} = [ITI_rate; vertcat(rates{:})];
      else
        HCP_aggr(r).Block.BinPunRate{p} = zeros(n_bins,1);
      end
      % Bin UnpPlanet clicks
      if ~isempty(HCP_aggr(r).Block.UnpClicks{p})
        [ITI_rate, ITI_time, rates, times] = bin_rate2(HCP_aggr(r).Block.UnpClicks{p},...
          HCP_aggr(r).Block.End(p), windows);
        HCP_aggr(r).Block.BinUnpRate{p} = [ITI_rate; vertcat(rates{:})];
      else
        HCP_aggr(r).Block.BinUnpRate{p} = zeros(n_bins,1);
      end      
      % Bin Shield button clicks
      if ~isempty(HCP_aggr(r).Block.ShieldClicks{p})
        [ITI_rate, ITI_time, rates, times] = bin_rate2(HCP_aggr(r).Block.ShieldClicks{p},...
          HCP_aggr(r).Block.End(p), windows);
        HCP_aggr(r).Block.BinShieldRate{p} = [ITI_rate; vertcat(rates{:})];
      else
        HCP_aggr(r).Block.BinShieldRate{p} = zeros(n_bins,1);
      end      
      
      HCP_aggr(r).Block.BinTime{p} = [ITI_time; times];
      
    else % not pun trial
      % PunRates
      HCP_aggr(r).Block.BinPunRate{p} = NaN(n_bins,1);
      HCP_aggr(r).Block.BinUnpRate{p} = NaN(n_bins,1);
      
      if ~isempty(HCP_aggr(r).Block.PunClicks{p})
        HCP_aggr(r).Block.BinPunRate{p}(1) = ...
          size(HCP_aggr(r).Block.PunClicks{p},1)/HCP_aggr(r).Block.End(p);
      else
        HCP_aggr(r).Block.BinPunRate{p}(1) = 0;
      end
      % UnpRates (only ITI)
      if ~isempty(HCP_aggr(r).Block.UnpClicks{p})
        HCP_aggr(r).Block.BinUnpRate{p}(1) = ...
          size(HCP_aggr(r).Block.UnpClicks{p},1)/HCP_aggr(r).Block.End(p);
      else
        HCP_aggr(r).Block.BinUnpRate{p}(1) = 0;
      end
      
      HCP_aggr(r).Block.BinTime{p} = [HCP_aggr(r).Block.End(p); zeros(12,1)];    
    end
    
    %% Plotting
    % Click location figure
    if ~isempty(click_loc{p})
      figure(click_fig); hold on
      scatter(click_loc{p}(:,1),-click_loc{p}(:,2),[],col_rep(p)); 

      figure(normclick_fig); hold on
      scatter(HCP_aggr(r).Block.NormClickLoc{p}(:,1),HCP_aggr(r).Block.NormClickLoc{p}(:,2),[],col_rep(p)); 
    end
    
    %% Plot intra-phase data
    figure(t_fig);
    subplot(2,1,1); hold on
    
    %Clicks 
    plot((HCP_aggr(r).Block.PunClicks{p}+prev_times),ones(sum(pun_clicks{p}),1)*ev_ytick(4),...
      '.','Color',col_rep(2),'LineStyle','none');
    plot((HCP_aggr(r).Block.UnpClicks{p}+prev_times),ones(sum(unp_clicks{p}),1)*ev_ytick(5),...
      '.','Color',col_rep(3),'LineStyle','none');
    plot((HCP_aggr(r).Block.ShieldClicks{p}+prev_times),ones(sum(shield_clicks{p}),1)*ev_ytick(6),...
      '.','Color',col_rep(1),'LineStyle','none');
    plot((HCP_aggr(r).Block.OtherClicks{p}+prev_times),ones(sum(other_clicks{p}),1)*ev_ytick(7),...
      'k.','LineStyle','none');
    
    %Rewards
    idx = HCP_aggr(r).Block.PunRew{p}(:,1)==100;
    plot((HCP_aggr(r).Block.PunRew{p}(idx,2)+prev_times),ones(sum(idx),1)*ev_ytick(8)-.1,...
      '+','Color',col_rep(2),'LineStyle','none');
    idx = HCP_aggr(r).Block.UnpRew{p}(:,1)==100;
    plot((HCP_aggr(r).Block.UnpRew{p}(idx,2)+prev_times),ones(sum(idx),1)*ev_ytick(8)+.1,...
      '+','Color',col_rep(3),'LineStyle','none');
    if ~isempty(HCP_aggr(r).Block.Attack{p})
      plot((HCP_aggr(r).Block.Attack{p}(:,2)+prev_times),...
        ones(size(HCP_aggr(r).Block.Attack{p},1),1)*ev_ytick(9),...
        'x','Color',col_rep(2),'LineStyle','none');
    end
    if ~isempty(HCP_aggr(r).Block.ShieldCost{p})
    plot((HCP_aggr(r).Block.ShieldCost{p}(:,2)+prev_times),...
      ones(size(HCP_aggr(r).Block.ShieldCost{p},1),1)*ev_ytick(9),...
      'k+','LineStyle','none');
    end
    
    if ismember(p,pun_trials)
    %Ships & shields
    for s = 1:size(HCP_aggr(r).Block.PunShips{p},1)
      plot([(HCP_aggr(r).Block.PunShips{p}(s,1)+prev_times) ...
        (HCP_aggr(r).Block.PunShips{p}(s,2)+prev_times)], ...
      	[ev_ytick(1) ev_ytick(1)],'Color',col_rep(2),'LineWidth',4);
      plot([(HCP_aggr(r).Block.PunShips{p}(s,2)+prev_times) ...
        (HCP_aggr(r).Block.PunShips{p}(s,2)+HCP_parameters.feedback_duration/1000 ...
        +prev_times)], [ev_ytick(1) ev_ytick(1)],'Color',col_rep(14),'LineWidth',4);        
    end 
    for s = 1:size(HCP_aggr(r).Block.UnpShips{p},1)
      plot([(HCP_aggr(r).Block.UnpShips{p}(s,1)+prev_times) ...
        (HCP_aggr(r).Block.UnpShips{p}(s,2)+prev_times)], ...
      	[ev_ytick(2) ev_ytick(2)],'Color',col_rep(3),'LineWidth',4);
      plot([(HCP_aggr(r).Block.UnpShips{p}(s,2)+prev_times) ...
        (HCP_aggr(r).Block.UnpShips{p}(s,2)+HCP_parameters.feedback_duration/1000 ...
        +prev_times)], [ev_ytick(2) ev_ytick(2)],'Color',col_rep(14),'LineWidth',4);        
    end 
    for s = 1:size(HCP_aggr(r).Block.PunShields{p},1)
      if ~isnan(HCP_aggr(r).Block.PunShields{p}(s,1))
        plot([(HCP_aggr(r).Block.PunShields{p}(s,1)+prev_times) ...
          (HCP_aggr(r).Block.PunShields{p}(s,3)+prev_times)], ...
          [ev_ytick(3) ev_ytick(3)],'Color',col_rep(1),'LineWidth',4);
        if ~isnan(HCP_aggr(r).Block.PunShields{p}(s,2))
          plot([(HCP_aggr(r).Block.PunShields{p}(s,1)+...
            HCP_aggr(r).Block.PunShields{p}(s,2)+prev_times) ...
            (HCP_aggr(r).Block.PunShields{p}(s,3)+prev_times)], ...
            [ev_ytick(3) ev_ytick(3)],'Color',col_rep(3),'LineWidth',4);
        end
      end
    end  
    for s = 1:size(HCP_aggr(r).Block.UnpShields{p},1)
      if ~isnan(HCP_aggr(r).Block.UnpShields{p}(s,1))
        plot([(HCP_aggr(r).Block.UnpShields{p}(s,1)+prev_times) ...
          (HCP_aggr(r).Block.UnpShields{p}(s,3)+prev_times)], ...
          [ev_ytick(3) ev_ytick(3)],'Color',col_rep(1),'LineWidth',4);
        if ~isnan(HCP_aggr(r).Block.UnpShields{p}(s,2))
          plot([(HCP_aggr(r).Block.UnpShields{p}(s,1)+HCP_aggr(r).Block.UnpShields{p}(s,2)+prev_times) ...
            (HCP_aggr(r).Block.UnpShields{p}(s,3)+prev_times)], ...
            [ev_ytick(3) ev_ytick(3)],'Color',col_rep(3),'LineWidth',4);
        end
      end
    end 
    end
    
    % Points
    subplot(2,1,2); hold on
    plot((HCP_aggr(r).raw{phase_idx(p),3}.all_outcomes.time_outcome/1000+prev_times),...
      HCP_aggr(r).raw{phase_idx(p),3}.all_outcomes.total,'d-','Color',col_rep(1));   

    prev_times = prev_times + HCP_aggr(r).Block.End(p) + 1;
  end
  
  %% Average RT for Val/Inf questions
  HCP_aggr(r).ValInf_AvgRT = vertcat(HCP_aggr(r).Block.ValCheckRT,...
    HCP_aggr(r).Block.InfCheck1RT, HCP_aggr(r).Block.InfCheck2RT,...
    HCP_aggr(r).Block.InfCheckSh1RT, HCP_aggr(r).Block.InfCheckSh2RT)./RT_check_grid;
  fprintf(['Val/Inf average RT per question (per screen): mean = ' ...
    num2str(mean(HCP_aggr(r).ValInf_AvgRT,'all','omitnan')) ...
    's (range = ' num2str(min(HCP_aggr(r).ValInf_AvgRT,[],'all','omitnan')) ...
    ' - ' num2str(max(HCP_aggr(r).ValInf_AvgRT,[],'all','omitnan')) 's)\n']);
  
  %% Questionnaires
  HCP_aggr(r).Catch = NaN(length(ques_names),1);

  for q = 1:length(ques_names)
    q_idx = is_idx(ismember(HCP_aggr(r).raw(is_idx,1),ques_names{q}));

    idx = regexp(HCP_aggr(r).raw{q_idx,3}.responses,':');
    HCP_aggr(r).(ques_names{q}) = HCP_aggr(r).raw{q_idx,3}.responses(idx+1) - '0';

    if ~isempty(ques_checks{q}) % Has catch question
      if HCP_aggr(r).(ques_names{q})(ques_checks{q}(1)) ==  ques_checks{q}(2)
        HCP_aggr(r).Catch(q) = false;
      else
        HCP_aggr(r).Catch(q) = true;
      end
      
      % Remove catch answer from questionnaire vector
      idx = 1:length(HCP_aggr(r).(ques_names{q}));
      idx = setdiff(idx,ques_checks{q}(1));
      HCP_aggr(r).(ques_names{q}) = HCP_aggr(r).(ques_names{q})(idx);
    end
  end

  %% Failed checks line
  if ~any(HCP_aggr(r).Catch) && ...
      ~any(HCP_aggr(r).ValInf_AvgRT<RT_lowcutoff,'all') && ~any(HCP_aggr(r).ValInf_AvgRT>RT_highcutoff,'all')
    cq_line = '';
    HCP_aggr(r).Exclude = false;
  elseif any(HCP_aggr(r).Catch) ...
      && any(any(HCP_aggr(r).ValInf_AvgRT<RT_lowcutoff,'all')|any(HCP_aggr(r).ValInf_AvgRT>RT_highcutoff,'all'))
    fprintf('FAILED BOTH CHECKS\n');
    cq_line = ' (Catch_RT_Fail)';
    HCP_aggr(r).Exclude = true;
  elseif any(HCP_aggr(r).Catch)
    fprintf('FAILED CATCH QUESTIONS\n');
    cq_line = ' (Catch_Fail)';
    HCP_aggr(r).Exclude = true;
  elseif any(HCP_aggr(r).ValInf_AvgRT<RT_lowcutoff,'all') || any(HCP_aggr(r).ValInf_AvgRT>RT_highcutoff,'all')
    fprintf('INVALID AVG Val/Inf RT VALUES\n');  
    cq_line = ' (RT_Fail)';
    HCP_aggr(r).Exclude = true;
  else
    error('Check check fail?');
  end
  
  %% 
  figure(t_fig);
  subplot(2,1,1); hold on
  prev_times = 0;
  for p = 1:n_phases
    prev_times = prev_times + HCP_aggr(r).Block.End(p) + 1;
    plot([prev_times-1 prev_times-1],ylim,'k--');
    plot([prev_times prev_times],ylim,'k--');
  end
  xlim([0 prev_times]);
  xlabel('Seconds');
  ylim([min(ev_ytick)*.5 max(ev_ytick)+min(ev_ytick)*.5]);
  yticks(ev_ytick(ytick_idx));
  yticklabels(ev_ylabel(ytick_idx));
  
  subplot(2,1,2); hold on
  prev_times = 0;
  for p = 1:n_phases
    ana_struct = ['Block' num2str(p)]; 
    prev_times = prev_times + HCP_aggr(r).Block.End(p) + 1;
    plot([prev_times-1 prev_times-1],ylim,'k--');
    plot([prev_times prev_times],ylim,'k--');
  end
  xlim([0 prev_times]);
  xlabel('Seconds');
  ylabel('Total Points');
    
  set(t_fig,'Position',get(0,'Screensize'));
  saveas(t_fig,[work_folder '\Figs\Timecourse\aggr' num2str(r) cq_line ...
  	' timecourse (' HCP_aggr(r).json ',set' num2str(HCP_aggr(r).Set) ').png']);
  close(t_fig);

  %% Click figures
  figure(click_fig); hold on   
  ylim([-HCP_aggr(r).raw{phase_idx(1),3}.viewport_size(2) 0]);
  xlim([0 HCP_aggr(r).raw{phase_idx(1),3}.viewport_size(1)]);
  legend({'Block1' 'Block2' 'Block3' 'Block4' 'Block5'});
  title(['Click locations (PunP = ' HCP_aggr(r).PunPlanet ')']);
  set(click_fig,'Position',get(0,'Screensize'));
  saveas(click_fig,[work_folder '\Figs\Click location\aggr' num2str(r) cq_line ...
  	' click location (' HCP_aggr(r).json ',set' num2str(HCP_aggr(r).Set) ').png']);
  close(click_fig);
  
  figure(normclick_fig); hold on   
  ylim([-1 1]);
  xlim([-1 1]);
  title(['Normalised click locations (PunP = ' HCP_aggr(r).PunPlanet ')']);
  legend({'Block1' 'Block2' 'Block3' 'Block4' 'Block5'});
  set(normclick_fig,'Position',get(0,'Screensize'));
  saveas(normclick_fig,[work_folder '\Figs\Click location\aggr' num2str(r) cq_line ...
  	' norm click location (' HCP_aggr(r).json ',set' num2str(HCP_aggr(r).Set) ').png']);
  close(normclick_fig);

  %% Summary plot    
  sum_fig = HCP_sumfig2(HCP_aggr(r).Block,pun_trials);
  
  % Save
  set(sum_fig,'Position',get(0,'Screensize'));
  saveas(figure(1), [work_folder '\Figs\Summary\aggr' num2str(r) cq_line ...
  	' summary (' HCP_aggr(r).json ',set' num2str(HCP_aggr(r).Set) ').png']);
  close all
end

%% Finish up
fprintf(['\nHCP extracted. N = ' num2str(length(HCP_aggr)) ' (-' ...
  num2str(sum(vertcat(HCP_aggr(:).Exclude))) ' with failed checks). ']);
toc(start_tic_keep);
fprintf('Saving... ');

clearvars -except HCP* work_folder *keep;
save([work_folder '\HCP.mat']);

fprintf(['\nSaved (' char(datetime) '). ']);
toc(start_tic_keep);
clearvars start_tic_keep;

diary OFF