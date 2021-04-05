%% HCP summary v2
% Run after HCPextract

n_blocks = 5;
unp_blocks = 1:2;
pun_blocks = 3:6;

bin_var = HCP_parameters.BinLabels;

map_sz = [64 40];
loc_lim = []; % color limits for click heatmap (normalized to # subjs+blocks)

%% Prep variables
clear HCPsum
close all

blockstruct = vertcat(HCP_aggr(:).Block);

n_subj = length(HCP_aggr);
exclude_idx = vertcat(HCP_aggr(:).Exclude);

ValInf_vars = fieldnames(HCP_aggr(1).Block);
ValInf_vars = ValInf_vars(endsWith(ValInf_vars,{'Val' 'Inf' 'Infconf'}));

HCPsum.nPunShldAv = NaN(n_subj,n_blocks);
HCPsum.nUnpShldAv = NaN(n_subj,n_blocks);
HCPsum.pcPunShldAv = NaN(n_subj,n_blocks);
HCPsum.pcUnpShldAv = NaN(n_subj,n_blocks);
HCPsum.PunShldRT = NaN(n_subj,n_blocks);
HCPsum.UnpShldRT = NaN(n_subj,n_blocks);
HCPsum.pcPunShldAv_PunB = NaN(n_subj,n_blocks);
HCPsum.pcUnpShldAv_PunB = NaN(n_subj,n_blocks);
HCPsum.PunShldRT_PunB = NaN(n_subj,n_blocks);
HCPsum.UnpShldRT_PunB = NaN(n_subj,n_blocks);

HCPsum.nPunShldTaken = NaN(n_subj,n_blocks);
HCPsum.nUnpShldTaken = NaN(n_subj,n_blocks);
HCPsum.pcPunShldTaken = NaN(n_subj,n_blocks);
HCPsum.pcUnpShldTaken = NaN(n_subj,n_blocks);
HCPsum.pcPunShldTaken_PunB = NaN(n_subj,1);
HCPsum.pcUnpShldTaken_PunB = NaN(n_subj,1);
HCPsum.Pref_pcShldTaken_PunB = NaN(n_subj,1);
HCPsum.PunS_suppr = NaN(n_subj,n_blocks);
HCPsum.UnpS_suppr = NaN(n_subj,n_blocks);
HCPsum.PunS_suppr_PunB = NaN(n_subj,1);
HCPsum.UnpS_suppr_PunB = NaN(n_subj,1);
HCPsum.PrefS_suppr_PunB = NaN(n_subj,1);

HCPsum.EndPoints = vertcat(blockstruct(:).EndPoints);

%% Fill HCPsum
n_var = length(ValInf_vars);
for v = 1:n_var
  HCPsum.(ValInf_vars{v}) = vertcat(blockstruct(:).(ValInf_vars{v}));
  
  % Pun Blocks
  HCPsum.([ValInf_vars{v} '_PunB']) = mean(HCPsum.(ValInf_vars{v})(:,pun_blocks),2);
  
  % Unp Blocks + PunB:UnpB ratio
  if any(HCPsum.(ValInf_vars{v})(:,unp_blocks),'all')
    HCPsum.([ValInf_vars{v} '_UnpB']) = mean(HCPsum.(ValInf_vars{v})(:,unp_blocks),2);
    HCPsum.([ValInf_vars{v} '_PunUnpBratio']) = HCPsum.([ValInf_vars{v} '_PunB'])./...
      (HCPsum.([ValInf_vars{v} '_PunB'])+HCPsum.([ValInf_vars{v} '_UnpB']));
  end
end

BinTime = cell2mat(vertcat(blockstruct(:).BinTime));
zero_idx = BinTime == 0;
BinPunRate = cell2mat(vertcat(blockstruct(:).BinPunRate));
BinPunRate(zero_idx) = 0;
BinPunClicks = BinPunRate.*BinTime;
BinUnpRate = cell2mat(vertcat(blockstruct(:).BinUnpRate));
BinUnpRate(zero_idx) = 0;
BinUnpClicks = BinUnpRate.*BinTime;
for r = 1:n_subj
  blockstruct(r).BinShieldRate{1} = zeros(13,1);
  blockstruct(r).BinShieldRate{2} = zeros(13,1);
end
BinShieldRate = cell2mat(vertcat(blockstruct(:).BinShieldRate));
BinShieldRate(zero_idx) = 0;
BinShieldClicks = BinShieldRate.*BinTime;

n_var = length(bin_var);
for v = 1:n_var
  HCPsum.([HCP_parameters.BinLabels{v} '_Time']) = BinTime(v:13:end,:);
  
  HCPsum.([HCP_parameters.BinLabels{v} '_PunClicks']) = BinPunClicks(v:13:end,:);
  HCPsum.([HCP_parameters.BinLabels{v} '_PunRate']) = BinPunRate(v:13:end,:);
  
  HCPsum.([HCP_parameters.BinLabels{v} '_UnpClicks']) = BinUnpClicks(v:13:end,:);
  HCPsum.([HCP_parameters.BinLabels{v} '_UnpRate']) = BinUnpRate(v:13:end,:);
  
  HCPsum.([HCP_parameters.BinLabels{v} '_ShieldClicks']) = BinShieldClicks(v:13:end,:);
  HCPsum.([HCP_parameters.BinLabels{v} '_ShieldRate']) = BinShieldRate(v:13:end,:);
  
  % Pun Blocks
  HCPsum.([HCP_parameters.BinLabels{v} '_Time_PunB']) = sum(BinTime(v:13:end,pun_blocks),2);
  
  HCPsum.([HCP_parameters.BinLabels{v} '_PunClicks_PunB']) = sum(BinPunClicks(v:13:end,pun_blocks),2);
  HCPsum.([HCP_parameters.BinLabels{v} '_PunRate_PunB']) = HCPsum.([HCP_parameters.BinLabels{v} '_PunClicks_PunB'])...
    ./HCPsum.([HCP_parameters.BinLabels{v} '_Time_PunB']);
  
  HCPsum.([HCP_parameters.BinLabels{v} '_UnpClicks_PunB']) = sum(BinUnpClicks(v:13:end,pun_blocks),2);
  HCPsum.([HCP_parameters.BinLabels{v} '_UnpRate_PunB']) = HCPsum.([HCP_parameters.BinLabels{v} '_UnpClicks_PunB'])...
    ./HCPsum.([HCP_parameters.BinLabels{v} '_Time_PunB']);
  
  HCPsum.([HCP_parameters.BinLabels{v} '_ShieldClicks_PunB']) = sum(BinShieldClicks(v:13:end,pun_blocks),2);
  HCPsum.([HCP_parameters.BinLabels{v} '_ShieldRate_PunB']) = HCPsum.([HCP_parameters.BinLabels{v} '_ShieldClicks_PunB'])...
    ./HCPsum.([HCP_parameters.BinLabels{v} '_Time_PunB']);

  if v == 1 % ITI for UnpB + PunB:UnpB ratio
    % Unp Blocks
    HCPsum.([HCP_parameters.BinLabels{v} '_Time_UnpB']) = sum(BinTime(v:13:end,unp_blocks),2);

    HCPsum.([HCP_parameters.BinLabels{v} '_PunClicks_UnpB']) = sum(BinPunClicks(v:13:end,unp_blocks),2);
    HCPsum.([HCP_parameters.BinLabels{v} '_PunRate_UnpB']) = HCPsum.([HCP_parameters.BinLabels{v} '_PunClicks_UnpB'])...
      ./HCPsum.([HCP_parameters.BinLabels{v} '_Time_UnpB']);

    HCPsum.([HCP_parameters.BinLabels{v} '_UnpClicks_UnpB']) = sum(BinUnpClicks(v:13:end,unp_blocks),2);
    HCPsum.([HCP_parameters.BinLabels{v} '_UnpRate_UnpB']) = HCPsum.([HCP_parameters.BinLabels{v} '_UnpClicks_UnpB'])...
      ./HCPsum.([HCP_parameters.BinLabels{v} '_Time_UnpB']);

    HCPsum.([HCP_parameters.BinLabels{v} '_ShieldClicks_UnpB']) = sum(BinShieldClicks(v:13:end,unp_blocks),2);
    HCPsum.([HCP_parameters.BinLabels{v} '_ShieldRate_UnpB']) = HCPsum.([HCP_parameters.BinLabels{v} '_ShieldClicks_UnpB'])...
      ./HCPsum.([HCP_parameters.BinLabels{v} '_Time_UnpB']);

    % Pun:Unp Block ratio
    HCPsum.([HCP_parameters.BinLabels{v} '_PunRate_PunUnpBratio']) = HCPsum.([HCP_parameters.BinLabels{v} '_PunRate_PunB'])...
      ./(HCPsum.([HCP_parameters.BinLabels{v} '_PunRate_PunB'])+HCPsum.([HCP_parameters.BinLabels{v} '_PunRate_UnpB']));

    HCPsum.([HCP_parameters.BinLabels{v} '_UnpRate_PunUnpBratio']) = HCPsum.([HCP_parameters.BinLabels{v} '_UnpRate_PunB'])...
      ./(HCPsum.([HCP_parameters.BinLabels{v} '_UnpRate_PunB'])+HCPsum.([HCP_parameters.BinLabels{v} '_UnpRate_UnpB']));
  else % Bin:ITI ratio   
    HCPsum.([HCP_parameters.BinLabels{v} '_ITIratio_Pun']) = HCPsum.([HCP_parameters.BinLabels{v} '_PunRate'])./...
      (HCPsum.([HCP_parameters.BinLabels{v} '_PunRate']) + HCPsum.ITI_PunRate);
    HCPsum.([HCP_parameters.BinLabels{v} '_ITIratio_Unp']) = HCPsum.([HCP_parameters.BinLabels{v} '_UnpRate'])./...
      (HCPsum.([HCP_parameters.BinLabels{v} '_UnpRate']) + HCPsum.ITI_UnpRate);
    HCPsum.([HCP_parameters.BinLabels{v} '_ITIratio_Both']) = ...
      (HCPsum.([HCP_parameters.BinLabels{v} '_PunRate'])+HCPsum.([HCP_parameters.BinLabels{v} '_UnpRate']))./...
      (HCPsum.([HCP_parameters.BinLabels{v} '_UnpRate'])+HCPsum.([HCP_parameters.BinLabels{v} '_PunRate'])...
      + HCPsum.ITI_PunRate+HCPsum.ITI_UnpRate);

    HCPsum.([HCP_parameters.BinLabels{v} '_ITIratio_Pun_PunB']) = HCPsum.([HCP_parameters.BinLabels{v} '_PunRate_PunB'])./...
      (HCPsum.([HCP_parameters.BinLabels{v} '_PunRate_PunB']) + HCPsum.ITI_PunRate_PunB);
    HCPsum.([HCP_parameters.BinLabels{v} '_ITIratio_Unp_PunB']) = HCPsum.([HCP_parameters.BinLabels{v} '_UnpRate_PunB'])./...
      (HCPsum.([HCP_parameters.BinLabels{v} '_UnpRate_PunB']) + HCPsum.ITI_UnpRate_PunB);
    HCPsum.([HCP_parameters.BinLabels{v} '_ITIratio_Both_PunB']) = ...
      (HCPsum.([HCP_parameters.BinLabels{v} '_PunRate_PunB'])+HCPsum.([HCP_parameters.BinLabels{v} '_UnpRate_PunB']))./...
      (HCPsum.([HCP_parameters.BinLabels{v} '_UnpRate_PunB'])+HCPsum.([HCP_parameters.BinLabels{v} '_PunRate_PunB'])...
      + HCPsum.ITI_PunRate_PunB + HCPsum.ITI_UnpRate_PunB);  
  end
  
  % UnpITI ratio
  HCPsum.([HCP_parameters.BinLabels{v} '_PunRate_PunUnpITIratio']) = HCPsum.([HCP_parameters.BinLabels{v} '_PunRate_PunB'])...
    ./(HCPsum.([HCP_parameters.BinLabels{v} '_PunRate_PunB'])+HCPsum.ITI_PunRate_UnpB);

  HCPsum.([HCP_parameters.BinLabels{v} '_UnpRate_PunUnpITIratio']) = HCPsum.([HCP_parameters.BinLabels{v} '_UnpRate_PunB'])...
    ./(HCPsum.([HCP_parameters.BinLabels{v} '_UnpRate_PunB'])+HCPsum.ITI_UnpRate_UnpB);
end

%% Unshielded rates + suppr
HCPsum.PunS_Unshld_Time = (BinTime(2:13:end,:)+BinTime(3:13:end,:)+BinTime(5:13:end,:));
HCPsum.PunS_Unshld_PunRate = ...
  (BinPunClicks(2:13:end,:)+BinPunClicks(3:13:end,:)+BinPunClicks(5:13:end,:))./HCPsum.PunS_Unshld_Time;
HCPsum.PunS_Unshld_UnpRate = ...
  (BinUnpClicks(2:13:end,:)+BinUnpClicks(3:13:end,:)+BinUnpClicks(5:13:end,:))./HCPsum.PunS_Unshld_Time;
HCPsum.PunS_Unshld_ShieldRate = ...
  (BinShieldClicks(2:13:end,:)+BinShieldClicks(3:13:end,:)+BinShieldClicks(5:13:end,:))./HCPsum.PunS_Unshld_Time;

HCPsum.PunS_Unshld_PunRate_PunB = ...
  sum((BinPunClicks(2:13:end,pun_blocks)+BinPunClicks(3:13:end,pun_blocks)+BinPunClicks(5:13:end,pun_blocks)),2)...
  ./sum(HCPsum.PunS_Unshld_Time(:,pun_blocks),2);
HCPsum.PunS_Unshld_UnpRate_PunB = ...
  sum((BinUnpClicks(2:13:end,pun_blocks)+BinUnpClicks(3:13:end,pun_blocks)+BinUnpClicks(5:13:end,pun_blocks)),2)...
  ./sum(HCPsum.PunS_Unshld_Time(:,pun_blocks),2);

HCPsum.UnpS_Unshld_Time = (BinTime(8:13:end,:)+BinTime(9:13:end,:)+BinTime(11:13:end,:));
HCPsum.UnpS_Unshld_PunRate = ...
  (BinPunClicks(8:13:end,:)+BinPunClicks(9:13:end,:)+BinPunClicks(11:13:end,:))./HCPsum.UnpS_Unshld_Time;
HCPsum.UnpS_Unshld_UnpRate = ...
  (BinUnpClicks(8:13:end,:)+BinUnpClicks(9:13:end,:)+BinUnpClicks(11:13:end,:))./HCPsum.UnpS_Unshld_Time;
HCPsum.UnpS_Unshld_ShieldRate = ...
  (BinShieldClicks(8:13:end,:)+BinShieldClicks(9:13:end,:)+BinShieldClicks(11:13:end,:))./HCPsum.UnpS_Unshld_Time;

HCPsum.UnpS_Unshld_PunRate_PunB = ...
  sum((BinPunClicks(8:13:end,pun_blocks)+BinPunClicks(9:13:end,pun_blocks)+BinPunClicks(11:13:end,pun_blocks)),2)...
  ./sum(HCPsum.UnpS_Unshld_Time(:,pun_blocks),2);
HCPsum.UnpS_Unshld_UnpRate_PunB = ...
  sum((BinUnpClicks(8:13:end,pun_blocks)+BinUnpClicks(9:13:end,pun_blocks)+BinUnpClicks(11:13:end,pun_blocks)),2)...
  ./sum(HCPsum.UnpS_Unshld_Time(:,pun_blocks),2);

HCPsum.PunS_suppr = (HCPsum.PunS_Unshld_PunRate+HCPsum.PunS_Unshld_UnpRate)./...
  ((HCPsum.PunS_Unshld_PunRate+HCPsum.PunS_Unshld_UnpRate)+(HCPsum.ITI_PunRate+HCPsum.ITI_UnpRate));
HCPsum.UnpS_suppr = (HCPsum.UnpS_Unshld_PunRate+HCPsum.UnpS_Unshld_UnpRate)./...
  ((HCPsum.UnpS_Unshld_PunRate+HCPsum.UnpS_Unshld_UnpRate)+(HCPsum.ITI_PunRate+HCPsum.ITI_UnpRate));
HCPsum.PunS_suppr_PunB = (HCPsum.PunS_Unshld_PunRate_PunB+HCPsum.PunS_Unshld_UnpRate_PunB)./...
  (HCPsum.PunS_Unshld_PunRate_PunB+HCPsum.PunS_Unshld_UnpRate_PunB+...
  HCPsum.ITI_PunRate_PunB+HCPsum.ITI_UnpRate_PunB);
HCPsum.UnpS_suppr_PunB = (HCPsum.UnpS_Unshld_PunRate_PunB+HCPsum.UnpS_Unshld_UnpRate_PunB)./...
  (HCPsum.UnpS_Unshld_PunRate_PunB+HCPsum.UnpS_Unshld_UnpRate_PunB+...
  HCPsum.ITI_PunRate_PunB+HCPsum.ITI_UnpRate_PunB);
HCPsum.PrefS_suppr_PunB = HCPsum.PunS_suppr_PunB./(HCPsum.PunS_suppr_PunB+HCPsum.UnpS_suppr_PunB);

%% ITI ratios
HCPsum.Pref = HCPsum.ITI_PunRate./(HCPsum.ITI_PunRate+HCPsum.ITI_UnpRate);
HCPsum.Pref_PunB = mean(HCPsum.Pref(:,pun_blocks),2);
HCPsum.Pref_UnpB = mean(HCPsum.Pref(:,unp_blocks),2);
HCPsum.Pref_PunUnpBratio = HCPsum.Pref_PunB./(HCPsum.Pref_UnpB+HCPsum.Pref_PunB);

HCPsum.PrefVal = HCPsum.PunPlanetVal./(HCPsum.PunPlanetVal+HCPsum.UnpPlanetVal);
HCPsum.PrefVal_PunB = mean(HCPsum.PrefVal(:,pun_blocks),2);
HCPsum.PrefVal_UnpB = mean(HCPsum.PrefVal(:,unp_blocks),2);
HCPsum.PrefVal_PunUnpBratio = HCPsum.PrefVal_PunB./(HCPsum.PrefVal_UnpB+HCPsum.PrefVal_PunB);

HCPsum.PrefRewInf = HCPsum.PunPlanet_RewInf./(HCPsum.PunPlanet_RewInf+HCPsum.UnpPlanet_RewInf);
HCPsum.PrefRewInf_PunB = mean(HCPsum.PrefRewInf(:,pun_blocks),2);
HCPsum.PrefRewInf_UnpB = mean(HCPsum.PrefRewInf(:,unp_blocks),2);
HCPsum.PrefRewInf_PunUnpBratio = HCPsum.PrefRewInf_PunB./(HCPsum.PrefRewInf_UnpB+HCPsum.PrefRewInf_PunB);

HCPsum.PrefRewInfconf = HCPsum.PunPlanet_RewInfconf./(HCPsum.PunPlanet_RewInfconf+HCPsum.UnpPlanet_RewInfconf);
HCPsum.PrefRewInfconf_PunB = mean(HCPsum.PrefRewInfconf(:,pun_blocks),2);
HCPsum.PrefRewInfconf_UnpB = mean(HCPsum.PrefRewInfconf(:,unp_blocks),2);
HCPsum.PrefRewInfconf_PunUnpBratio = HCPsum.PrefRewInfconf_PunB./(HCPsum.PrefRewInfconf_UnpB+HCPsum.PrefRewInfconf_PunB);

HCPsum.PrefPAtkInf = HCPsum.PunPlanet_AttackInf./(HCPsum.PunPlanet_AttackInf+HCPsum.UnpPlanet_AttackInf);
HCPsum.PrefPAtkInf_PunB = mean(HCPsum.PrefPAtkInf(:,pun_blocks),2);

HCPsum.PunTsuppr = HCPsum.ITI_PunRate./(HCPsum.ITI_PunRate+HCPsum.ITI_PunRate(:,2));
HCPsum.UnpTsuppr = HCPsum.ITI_UnpRate./(HCPsum.ITI_UnpRate+HCPsum.ITI_UnpRate(:,2));
HCPsum.PunValChange = HCPsum.PunPlanetVal./(HCPsum.PunPlanetVal+HCPsum.PunPlanetVal(:,2));
HCPsum.UnpValChange = HCPsum.UnpPlanetVal./(HCPsum.UnpPlanetVal+HCPsum.UnpPlanetVal(:,2));

%% Ships & Shields (n, percent, rt)
HCPsum.nPunShips = cellfun(@(x) size(x,1), vertcat(blockstruct.PunShields));
PunShields = cellfun(@(x) sum(~isnan(x),1), vertcat(blockstruct.PunShields),'UniformOutput',false);
idx = cellfun(@(x) ~isempty(x), PunShields);
HCPsum.nPunShldAv(idx) = cellfun(@(x) x(1), PunShields(idx));
HCPsum.nPunShldAv(~idx) = 0;
HCPsum.nPunShldTaken(idx) = cellfun(@(x) x(2), PunShields(idx));
HCPsum.nPunShldTaken(~idx) = 0;
HCPsum.pcPunShldAv = HCPsum.nPunShldAv./HCPsum.nPunShips;
HCPsum.pcPunShldTaken = HCPsum.nPunShldTaken./HCPsum.nPunShldAv;
PunShields = cellfun(@(x) mean(x,1,'omitnan'),vertcat(blockstruct.PunShields),'UniformOutput',false);
HCPsum.PunShldRT(idx) = cellfun(@(x) x(2), PunShields(idx));

HCPsum.pcPunShldAv_PunB = sum(HCPsum.nPunShldAv(:,pun_blocks),2)./sum(HCPsum.nPunShips(:,pun_blocks),2);
HCPsum.pcPunShldTaken_PunB = sum(HCPsum.nPunShldTaken(:,pun_blocks),2)./sum(HCPsum.nPunShldAv(:,pun_blocks),2);
HCPsum.PunShldRT_PunB = sum((HCPsum.PunShldRT(:,pun_blocks).*HCPsum.nPunShldTaken(:,pun_blocks)),2,'omitnan')...
  ./sum(HCPsum.nPunShldTaken(:,pun_blocks),2);

HCPsum.nUnpShips = cellfun(@(x) size(x,1), vertcat(blockstruct.UnpShields));
UnpShields = cellfun(@(x) sum(~isnan(x),1), vertcat(blockstruct.UnpShields),'UniformOutput',false);
idx = cellfun(@(x) ~isempty(x), UnpShields);
HCPsum.nUnpShldAv(idx) = cellfun(@(x) x(1), UnpShields(idx));
HCPsum.nUnpShldAv(~idx) = 0;
HCPsum.nUnpShldTaken(idx) = cellfun(@(x) x(2), UnpShields(idx));
HCPsum.nUnpShldTaken(~idx) = 0;
HCPsum.pcUnpShldAv = HCPsum.nUnpShldAv./HCPsum.nUnpShips;
HCPsum.pcUnpShldTaken = HCPsum.nUnpShldTaken./HCPsum.nUnpShldAv;
UnpShields = cellfun(@(x) mean(x,1,'omitnan'),vertcat(blockstruct.UnpShields),'UniformOutput',false);
HCPsum.UnpShldRT(idx) = cellfun(@(x) x(2), UnpShields(idx));

HCPsum.pcUnpShldAv_PunB = sum(HCPsum.nUnpShldAv(:,pun_blocks),2)./sum(HCPsum.nUnpShips(:,pun_blocks),2);
HCPsum.pcUnpShldTaken_PunB = sum(HCPsum.nUnpShldTaken(:,pun_blocks),2)./sum(HCPsum.nUnpShldAv(:,pun_blocks),2);
HCPsum.UnpShldRT_PunB = mean((HCPsum.UnpShldRT(:,pun_blocks).*HCPsum.nUnpShldTaken(:,pun_blocks)),2,'omitnan')...
  ./sum(HCPsum.nUnpShldTaken(:,pun_blocks),2);

%% Planet-Attack chain inference estimate
HCPsum.PunP_AtkInfEst = ...
  (HCPsum.PunPlanet_PunShipInf.*HCPsum.PunShip_AttackInf/100)+...
  (HCPsum.PunPlanet_UnpShipInf.*HCPsum.UnpShip_AttackInf/100);
HCPsum.UnpP_AtkInfEst = ...
  (HCPsum.UnpPlanet_PunShipInf.*HCPsum.PunShip_AttackInf/100)+...
  (HCPsum.UnpPlanet_UnpShipInf.*HCPsum.UnpShip_AttackInf/100);

HCPsum.PunP_PunS_AtkInfEst = HCPsum.PunPlanet_PunShipInf.*HCPsum.PunShip_AttackInf/100;
HCPsum.PunP_UnpS_AtkInfEst = HCPsum.PunPlanet_UnpShipInf.*HCPsum.UnpShip_AttackInf/100;
HCPsum.UnpP_UnpS_AtkInfEst = HCPsum.UnpPlanet_UnpShipInf.*HCPsum.UnpShip_AttackInf/100;
HCPsum.UnpP_PunS_AtkInfEst = HCPsum.UnpPlanet_PunShipInf.*HCPsum.PunShip_AttackInf/100;

%% Figures
%% Click locations
clicks = vertcat(blockstruct(:).NormClickLoc);
left_idx = ismember({HCP_aggr(:).PunPlanet},'left')';

avgUnpBleftP_clicks = vertcat(clicks{left_idx&~exclude_idx,unp_blocks});
avgUnpBleftP_clickmap = histcounts2(avgUnpBleftP_clicks(:,1),avgUnpBleftP_clicks(:,2)...
  ,-1:(2/map_sz(1)):1,-1:(2/map_sz(2)):1)/(length(unp_blocks)*sum(left_idx&~exclude_idx));
avgUnpBrightP_clicks = vertcat(clicks{~left_idx&~exclude_idx,unp_blocks});
avgUnpBrightP_clickmap = histcounts2(avgUnpBrightP_clicks(:,1),avgUnpBrightP_clicks(:,2)...
  ,-1:(2/map_sz(1)):1,-1:(2/map_sz(2)):1)/(length(unp_blocks)*sum(~left_idx&~exclude_idx));
avgPunBleftP_clicks = vertcat(clicks{left_idx&~exclude_idx,pun_blocks});
avgPunBleftP_clickmap = histcounts2(avgPunBleftP_clicks(:,1),avgPunBleftP_clicks(:,2)...
  ,-1:(2/map_sz(1)):1,-1:(2/map_sz(2)):1)/(length(pun_blocks)*sum(left_idx&~exclude_idx));
avgPunBrightP_clicks = vertcat(clicks{~left_idx&~exclude_idx,pun_blocks});
avgPunBrightP_clickmap = histcounts2(avgPunBrightP_clicks(:,1),avgPunBrightP_clicks(:,2)...
  ,-1:(2/map_sz(1)):1,-1:(2/map_sz(2)):1)/(length(pun_blocks)*sum(~left_idx&~exclude_idx));

if isempty(loc_lim)
  minofmax = min([max(avgUnpBleftP_clickmap,[],'all') max(avgUnpBrightP_clickmap,[],'all') ...
    max(avgPunBleftP_clickmap,[],'all') max(avgPunBrightP_clickmap,[],'all')]);
  loc_lim = [0 minofmax];
end

figure; 
subplot(2,2,1); hold on
imagesc(avgUnpBleftP_clickmap',loc_lim)
%hist3(vertcat(clicks{left_idx&~exclude_idx,unp_blocks}),'Ctrs',{-1:(2/80):1 -1:(2/50):1},'CdataMode','auto')
title('Pre-pun average click locations (PunPlanet=left)');

subplot(2,2,2); hold on
imagesc(avgPunBleftP_clickmap',loc_lim)
%hist3(vertcat(clicks{left_idx&~exclude_idx,pun_blocks}),'Ctrs',{-1:(2/80):1 -1:(2/50):1},'CdataMode','auto')
title('Pun phase average click locations (PunPlanet=left)');

subplot(2,2,3); hold on
imagesc(avgUnpBrightP_clickmap',loc_lim)
%hist3(vertcat(clicks{~left_idx&~exclude_idx,unp_blocks}),'Ctrs',{-1:(2/80):1 -1:(2/50):1},'CdataMode','auto')
title('Pre-pun average click locations (PunPlanet=right)');

subplot(2,2,4); hold on
imagesc(avgPunBrightP_clickmap',loc_lim)
%hist3(vertcat(clicks{~left_idx&~exclude_idx,pun_blocks}),'Ctrs',{-1:(2/80):1 -1:(2/50):1},'CdataMode','auto')
title('Pun phase average click locations (PunPlanet=right)');

for p = 1:4
  subplot(2,2,p); hold on
  axis tight off 
  colormap('hot');
  colorbar;
end

set(gcf,'Position',get(0,'Screensize'));
saveas(figure(1),[work_folder '\Figs\Click heatmap.png']);
close all

% Location change
LvRdiff_leftP_UnpB = avgUnpBleftP_clickmap(1:map_sz(1)/2,:)-flipud(avgUnpBleftP_clickmap(map_sz(1)/2+1:end,:));
LvRdiff_leftP_PunB = avgPunBleftP_clickmap(1:map_sz(1)/2,:)-flipud(avgPunBleftP_clickmap(map_sz(1)/2+1:end,:));
LvRdiff_rightP_UnpB = avgUnpBrightP_clickmap(1:map_sz(1)/2,:)-flipud(avgUnpBrightP_clickmap(map_sz(1)/2+1:end,:));
LvRdiff_rightP_PunB = avgPunBrightP_clickmap(1:map_sz(1)/2,:)-flipud(avgPunBrightP_clickmap(map_sz(1)/2+1:end,:));
minofmax = max([abs(min([max(LvRdiff_rightP_UnpB,[],'all') max(LvRdiff_rightP_PunB,[],'all')])) ...
	abs(max([min(LvRdiff_rightP_UnpB,[],'all') min(LvRdiff_rightP_PunB,[],'all')]))]);
loc_lim = [-minofmax minofmax];

figure; 
subplot(2,2,1); hold on
imagesc(vertcat(LvRdiff_leftP_UnpB,LvRdiff_leftP_PunB)',loc_lim)
title('Left-Right [UnpBlock / PunBlock] (PunPlanet=left)','interpreter','none');
plot([map_sz(1)/2+0.5 map_sz(1)/2+0.5],[0.5 map_sz(2)+0.5],'k--');
drawnow;

subplot(2,2,3); hold on
imagesc(vertcat(LvRdiff_rightP_UnpB,LvRdiff_rightP_PunB)',loc_lim)
title('Left-Right [UnpBlock / PunBlock] (PunPlanet=right)','interpreter','none');
plot([map_sz(1)/2+0.5 map_sz(1)/2+0.5],[0.5 map_sz(2)+0.5],'k--');
drawnow;

UnpPunDiff_leftP = avgPunBleftP_clickmap-avgUnpBleftP_clickmap;
UnpPunDiff_rightP = avgPunBrightP_clickmap-avgUnpBrightP_clickmap;
minofmax = max([abs(min([max(UnpPunDiff_leftP,[],'all') max(UnpPunDiff_rightP,[],'all')])) ...
	abs(max([min(UnpPunDiff_leftP,[],'all') min(UnpPunDiff_rightP,[],'all')]))]);
loc_lim = [-minofmax minofmax];

subplot(2,2,2); hold on
imagesc(UnpPunDiff_leftP',loc_lim)
title('Pun-UnpBlock click location (PunPlanet=left)','interpreter','none');

subplot(2,2,4); hold on
imagesc(UnpPunDiff_rightP',loc_lim)
title('Pun-UnpBlock click location (PunPlanet=right)','interpreter','none');
for p = 1:4
  subplot(2,2,p); hold on
  axis tight off 
  colormap(whitejet(100));
  colorbar;
end
set(gcf,'Position',get(0,'Screensize'));
saveas(figure(1),[work_folder '\Figs\Click change heatmap.png']);
close all


%% Planet clicks (per sess)
% ITI rates
figure; 
subplot(2,3,1); hold on
errorbar(mean(HCPsum.ITI_PunRate(~exclude_idx,:)*60,1),sem(HCPsum.ITI_PunRate(~exclude_idx,:)*60),...
  'k','CapSize',2);
errorbar(mean(HCPsum.ITI_UnpRate(~exclude_idx,:)*60,1),sem(HCPsum.ITI_UnpRate(~exclude_idx,:)*60),...
  'k','CapSize',2);
p1 = plot(mean(HCPsum.ITI_PunRate(~exclude_idx,:)*60,1),'Color',col_rep(2),'LineWidth',2);
p2 = plot(mean(HCPsum.ITI_UnpRate(~exclude_idx,:)*60,1),'Color',col_rep(3),'LineWidth',2);
xlabel('Block');
xlim([0.5 n_blocks+0.5]);
xticks(1:n_blocks);
ylims = ylim;
ylims(1) = 0;
ylim(ylims);
ylabel('Clicks/min','interpreter','none');
title('ITI planet click rates')
legend([p1 p2],{['PunPlanet (n=' num2str(sum(~exclude_idx)) ')'] 'UnpPlanet'},...
  'Location','best','interpreter','none');
legend('boxoff');

% T suppr
subplot(2,3,2); hold on
title('ITI ratio (relative to Blk2)');
incl_idx = ~exclude_idx & ~isnan(HCPsum.PunTsuppr(:,2)) & ~isnan(HCPsum.UnpTsuppr(:,2));
errorbar(mean(HCPsum.PunTsuppr(incl_idx,:),1),sem(HCPsum.PunTsuppr(incl_idx,:)),...
  'k.','CapSize',2);
errorbar(mean(HCPsum.UnpTsuppr(incl_idx,:),1),sem(HCPsum.UnpTsuppr(incl_idx,:)),...
  'k.','CapSize',2);
plot(HCPsum.PunTsuppr(incl_idx,:)','Color',col_rep(2)+((1-col_rep(2)).*.5),'LineWidth',0.5);
plot(HCPsum.UnpTsuppr(incl_idx,:)','Color',col_rep(3)+((1-col_rep(3)).*.5),'LineWidth',0.5);
p1 = plot(mean(HCPsum.PunTsuppr(incl_idx,:),1),'Color',col_rep(2),'LineWidth',2);
p2 = plot(mean(HCPsum.UnpTsuppr(incl_idx,:),1),'Color',col_rep(3),'LineWidth',2);
xlabel('Block');
xlim([0.5 n_blocks+0.5]);
xticks(1:n_blocks);
plot(xlim,[0.5 0.5],'k:');
ylim([0 1]);
legend([p1 p2],{['Pun T_suppr (n=' num2str(sum(incl_idx)) ')'] 'Unp T_suppr'},...
  'Location','best','interpreter','none');
legend('boxoff');

% Ratios
subplot(2,3,3); hold on
title('Suppression Ratios');
plot(HCPsum.Pref(~exclude_idx,:)','Color',col_rep(14),'LineWidth',0.5);
errorbar(mean(HCPsum.Pref(~exclude_idx,:),1),sem(HCPsum.Pref(~exclude_idx,:)),...
  'k.','CapSize',2);
p1 = plot(mean(HCPsum.Pref(~exclude_idx,:),1),'k','LineWidth',1.5);
exclude_idx2 = any(exclude_idx | isnan(HCPsum.PunS_suppr(:,pun_blocks)),2);
errorbar(mean(HCPsum.PunS_suppr(~exclude_idx2,:),1),sem(HCPsum.PunS_suppr(~exclude_idx2,:)),...
  'k.','CapSize',2);
p2 = plot(mean(HCPsum.PunS_suppr(~exclude_idx2,:),1),'d--','Color',col_rep(2),'LineWidth',2,...
  'MarkerFaceColor',col_rep(2));
exclude_idx2 = any(exclude_idx | isnan(HCPsum.UnpS_suppr(:,pun_blocks)),2);
errorbar(mean(HCPsum.UnpS_suppr(~exclude_idx2,:),1),sem(HCPsum.UnpS_suppr(~exclude_idx2,:)),...
  'k.','CapSize',2);
p3 = plot(mean(HCPsum.UnpS_suppr(~exclude_idx2,:),1),'o--','Color',col_rep(3),'LineWidth',2,...
  'MarkerFaceColor',col_rep(3));
errorbar(mean(HCPsum.pcPunShldTaken(~exclude_idx,:),1,'omitnan'),sem(HCPsum.pcPunShldTaken(~exclude_idx,:)),...
  'k.','CapSize',2);
p4 = plot(mean(HCPsum.pcPunShldTaken(~exclude_idx2,:),1,'omitnan'),'d-','Color',col_rep(1),'LineWidth',2,...
  'MarkerFaceColor',col_rep(2));
errorbar(mean(HCPsum.pcUnpShldTaken(~exclude_idx,:),1,'omitnan'),sem(HCPsum.pcUnpShldTaken(~exclude_idx,:)),...
  'k.','CapSize',2);
p5 = plot(mean(HCPsum.pcUnpShldTaken(~exclude_idx2,:),1,'omitnan'),'o-','Color',col_rep(1),'LineWidth',2,...
  'MarkerFaceColor',col_rep(3));
xlabel('Block');
xticks(1:n_blocks);
xlim([0.5 n_blocks+0.5]);
plot(xlim,[0.5 0.5],'k:');
ylim([0 1]);
ylabel('Ratio');
legend([p1 p2 p3 p4 p5],{'ITI PunP v UnpP' 'PunS Suppr' 'UnpS Suppr' 'PunS Shld %' 'UnpS Shld %'},...
  'Location','best','interpreter','none');
legend('boxoff')

% Values
subplot(2,3,4); hold on
title('Values');
errorbar(mean(HCPsum.PunPlanetVal(~exclude_idx,:),1),sem(HCPsum.PunPlanetVal(~exclude_idx,:)),...
  'k.','CapSize',2);
errorbar(mean(HCPsum.UnpPlanetVal(~exclude_idx,:),1),sem(HCPsum.UnpPlanetVal(~exclude_idx,:)),...
  'k.','CapSize',2);
errorbar(mean(HCPsum.PunShipVal(~exclude_idx,:),1),sem(HCPsum.PunShipVal(~exclude_idx,:)),...
  'k.','CapSize',2);
errorbar(mean(HCPsum.UnpShipVal(~exclude_idx,:),1),sem(HCPsum.UnpShipVal(~exclude_idx,:)),...
  'k.','CapSize',2);
errorbar(mean(HCPsum.RewVal(~exclude_idx,:),1),sem(HCPsum.RewVal(~exclude_idx,:)),...
  'k.','CapSize',2);
errorbar(mean(HCPsum.AttackVal(~exclude_idx,:),1),sem(HCPsum.AttackVal(~exclude_idx,:)),...
  'k.','CapSize',2);
p1 = plot(mean(HCPsum.PunPlanetVal(~exclude_idx,:),1),'Color',col_rep(2),'LineWidth',2);
p2 = plot(mean(HCPsum.UnpPlanetVal(~exclude_idx,:),1),'Color',col_rep(3),'LineWidth',2);
p3 = plot(mean(HCPsum.PunShipVal(~exclude_idx,:),1),'d--','Color',col_rep(2),'LineWidth',2,...
  'MarkerFaceColor',col_rep(2));
p4 = plot(mean(HCPsum.UnpShipVal(~exclude_idx,:),1),'o--','Color',col_rep(3),'LineWidth',2,...
  'MarkerFaceColor',col_rep(3));
p5 = plot(mean(HCPsum.RewVal(~exclude_idx,:),1),'k+:','LineWidth',1.5,'MarkerSize',10);
p6 = plot(mean(HCPsum.AttackVal(~exclude_idx,:),1),'x-','Color',col_rep(2),'LineWidth',1.5,...
  'MarkerFaceColor',col_rep(2),'MarkerSize',10);
xlabel('Block');
xlim([0.5 n_blocks+0.5]);
xticks(1:n_blocks);
ylim([0 100]);
ylabel('Value rating');
legend([p1 p2 p3 p4 p5 p6],{'PunPlanet' 'UnpPlanet' 'PunShip' 'UnpShip' 'Reward' 'Attack'},...
  'Location','best','interpreter','none');
legend('boxoff');

% Inferences
subplot(2,3,5); hold on
title('Inferences');
errorbar(mean(HCPsum.PunPlanet_RewInf(~exclude_idx,:),1),sem(HCPsum.PunPlanet_RewInf(~exclude_idx,:)),...
  'k.','CapSize',2);
errorbar(mean(HCPsum.UnpPlanet_RewInf(~exclude_idx,:),1),sem(HCPsum.UnpPlanet_RewInf(~exclude_idx,:)),...
  'k.','CapSize',2);
errorbar(mean(HCPsum.PunPlanet_PunShipInf(~exclude_idx,:),1),sem(HCPsum.PunPlanet_PunShipInf(~exclude_idx,:)),...
  'k.','CapSize',2);
errorbar(mean(HCPsum.UnpPlanet_UnpShipInf(~exclude_idx,:),1),sem(HCPsum.UnpPlanet_UnpShipInf(~exclude_idx,:)),...
  'k.','CapSize',2);
errorbar(mean(HCPsum.PunPlanet_UnpShipInf(~exclude_idx,:),1),sem(HCPsum.PunPlanet_UnpShipInf(~exclude_idx,:)),...
  'k.','CapSize',2);
errorbar(mean(HCPsum.UnpPlanet_PunShipInf(~exclude_idx,:),1),sem(HCPsum.UnpPlanet_PunShipInf(~exclude_idx,:)),...
  'k.','CapSize',2);
errorbar(mean(HCPsum.PunShip_AttackInf(~exclude_idx,:),1),sem(HCPsum.PunShip_AttackInf(~exclude_idx,:)),...
  'k.','CapSize',2);
errorbar(mean(HCPsum.UnpShip_AttackInf(~exclude_idx,:),1),sem(HCPsum.UnpShip_AttackInf(~exclude_idx,:)),...
  'k.','CapSize',2);
errorbar(mean(HCPsum.PunPlanet_AttackInf(~exclude_idx,:),1),sem(HCPsum.PunPlanet_AttackInf(~exclude_idx,:)),...
  'k.','CapSize',2);
errorbar(mean(HCPsum.UnpPlanet_AttackInf(~exclude_idx,:),1),sem(HCPsum.UnpPlanet_AttackInf(~exclude_idx,:)),...
  'k.','CapSize',2);
p1 = plot(mean(HCPsum.PunPlanet_RewInf(~exclude_idx,:),1),'+:','Color',col_rep(2),...
  'LineWidth',1.5,'MarkerSize',10);
p2 = plot(mean(HCPsum.UnpPlanet_RewInf(~exclude_idx,:),1),'+:','Color',col_rep(3),...
  'LineWidth',1.5,'MarkerSize',10);
p3 = plot(mean(HCPsum.PunPlanet_PunShipInf(~exclude_idx,:),1),'d--','Color',col_rep(2),...
  'LineWidth',2,'MarkerFaceColor',col_rep(2));
p9 = plot(mean(HCPsum.PunPlanet_UnpShipInf(~exclude_idx,:),1),'o--','Color',col_rep(2),...
  'LineWidth',2,'MarkerFaceColor',col_rep(3));
p4 = plot(mean(HCPsum.UnpPlanet_UnpShipInf(~exclude_idx,:),1),'o--','Color',col_rep(3),...
  'LineWidth',2,'MarkerFaceColor',col_rep(3));
p10 = plot(mean(HCPsum.UnpPlanet_PunShipInf(~exclude_idx,:),1),'d--','Color',col_rep(3),...
  'LineWidth',2,'MarkerFaceColor',col_rep(2));
p5 = plot(mean(HCPsum.PunShip_AttackInf(~exclude_idx,:),1),'x--','Color',col_rep(2),...
  'LineWidth',1.5,'MarkerSize',10);
p6 = plot(mean(HCPsum.UnpShip_AttackInf(~exclude_idx,:),1),'x--','Color',col_rep(3),...
  'LineWidth',1.5,'MarkerSize',10);
p7 = plot(mean(HCPsum.PunPlanet_AttackInf(~exclude_idx,:),1),'x:','Color',col_rep(2),...
  'LineWidth',1.5,'MarkerSize',10);
p8 = plot(mean(HCPsum.UnpPlanet_AttackInf(~exclude_idx,:),1),'x:','Color',col_rep(3),...
  'LineWidth',1.5,'MarkerSize',10);
xlabel('Block');
xlim([0.5 n_blocks+0.5]);
xticks(1:n_blocks);
ylim([0 100]);
ylabel('% likelihood');
legend([p1 p2 p3 p9 p4 p10 p5 p6 p7 p8],...
  {'PunP-Rew' 'UnpP-Rew' 'PunP-PunS' 'PunP-UnpS' 'UnpP-UnpS' 'UnpP-PunS' ...
  'PunS-Attack' 'UnpS-Attack' 'PunP-Attack' 'UnpP-Attack'},...
  'Location','best','interpreter','none');
legend('boxoff');

% Confidence
subplot(2,3,6); hold on
title('Confidence');
errorbar(mean(HCPsum.PunPlanet_RewInfconf(~exclude_idx,:),1),sem(HCPsum.PunPlanet_RewInfconf(~exclude_idx,:)),...
  'k.','CapSize',2);
errorbar(mean(HCPsum.UnpPlanet_RewInfconf(~exclude_idx,:),1),sem(HCPsum.UnpPlanet_RewInfconf(~exclude_idx,:)),...
  'k.','CapSize',2);
errorbar(mean(HCPsum.PunPlanet_PunShipInfconf(~exclude_idx,:),1),sem(HCPsum.PunPlanet_PunShipInfconf(~exclude_idx,:)),...
  'k.','CapSize',2);
errorbar(mean(HCPsum.UnpPlanet_UnpShipInfconf(~exclude_idx,:),1),sem(HCPsum.UnpPlanet_UnpShipInfconf(~exclude_idx,:)),...
  'k.','CapSize',2);
errorbar(mean(HCPsum.PunShip_AttackInfconf(~exclude_idx,:),1),sem(HCPsum.PunShip_AttackInfconf(~exclude_idx,:)),...
  'k.','CapSize',2);
errorbar(mean(HCPsum.UnpShip_AttackInfconf(~exclude_idx,:),1),sem(HCPsum.UnpShip_AttackInfconf(~exclude_idx,:)),...
  'k.','CapSize',2);
errorbar(mean(HCPsum.PunPlanet_AttackInfconf(~exclude_idx,:),1),sem(HCPsum.PunPlanet_AttackInfconf(~exclude_idx,:)),...
  'k.','CapSize',2);
errorbar(mean(HCPsum.UnpPlanet_AttackInfconf(~exclude_idx,:),1),sem(HCPsum.UnpPlanet_AttackInfconf(~exclude_idx,:)),...
  'k.','CapSize',2);
p1 = plot(mean(HCPsum.PunPlanet_RewInfconf(~exclude_idx,:),1),'+:','Color',col_rep(2),...
  'LineWidth',1.5,'MarkerSize',10);
p2 = plot(mean(HCPsum.UnpPlanet_RewInfconf(~exclude_idx,:),1),'+:','Color',col_rep(3),...
  'LineWidth',1.5,'MarkerSize',10);
p3 = plot(mean(HCPsum.PunPlanet_PunShipInfconf(~exclude_idx,:),1),'d--','Color',col_rep(2),...
  'LineWidth',2,'MarkerFaceColor',col_rep(2));
p4 = plot(mean(HCPsum.UnpPlanet_UnpShipInfconf(~exclude_idx,:),1),'o--','Color',col_rep(3),...
  'LineWidth',2,'MarkerFaceColor',col_rep(3));
p5 = plot(mean(HCPsum.PunShip_AttackInfconf(~exclude_idx,:),1),'x--','Color',col_rep(2),...
  'LineWidth',1.5,'MarkerSize',10);
p6 = plot(mean(HCPsum.UnpShip_AttackInfconf(~exclude_idx,:),1),'x--','Color',col_rep(3),...
  'LineWidth',1.5,'MarkerSize',10);
p7 = plot(mean(HCPsum.PunPlanet_AttackInfconf(~exclude_idx,:),1),'x:','Color',col_rep(2),...
  'LineWidth',1.5,'MarkerSize',10);
p8 = plot(mean(HCPsum.UnpPlanet_AttackInfconf(~exclude_idx,:),1),'x:','Color',col_rep(3),...
  'LineWidth',1.5,'MarkerSize',10);
xlabel('Block');
xlim([0.5 n_blocks+0.5]);
xticks(1:n_blocks);
ylim([0 100]);
ylabel('Confidence rating');
legend([p1 p2 p3 p4 p5 p6 p7 p8],...
  {'PunP-Rew' 'UnpP-Rew' 'PunP-PunS' 'UnpP-UnpS' 'PunS-Attack' 'UnpS-Attack' 'PunP-Attack' 'UnpP-Attack'},...
  'Location','best','interpreter','none');
legend('boxoff');

set(gcf,'Position',get(0,'Screensize'));
saveas(figure(1),[work_folder '\Figs\Aggregate summary.png']);
close all
% 
% % Planet value change
% subplot(2,3,5); hold on
% exclude_idx2 = exclude_idx | isnan(HCPsum.PunValChange(:,2))';
% errorbar(mean(HCPsum.PunValChange(~exclude_idx2,:),1),sem(HCPsum.PunValChange(~exclude_idx2,:)),...
%   'k.','CapSize',2);
% p1 = plot(mean(HCPsum.PunValChange(~exclude_idx2,:),1),'Color',col_rep(2),'LineWidth',2);
% exclude_idx2 = exclude_idx | isnan(HCPsum.UnpValChange(:,2))';
% errorbar(mean(HCPsum.UnpValChange(~exclude_idx2,:),1),sem(HCPsum.UnpValChange(~exclude_idx2,:)),...
%   'k.','CapSize',2);
% p2 = plot(mean(HCPsum.UnpValChange(~exclude_idx2,:),1),'Color',col_rep(3),'LineWidth',2);
% xlabel('Block');
% xticks(1:n_blocks);
% xlim([0.5 n_blocks+0.5]);
% plot(xlim,[0.5 0.5],'k:');
% ylim([0 1]);
% ylabel('Planet value ratio (relative to Blk2)');
% legend([p1 p2],{'PunP value change' 'UnpP value change'},'Location','best','interpreter','none');
% legend('boxoff');

clearvars -except HCP* work_folder *keep;