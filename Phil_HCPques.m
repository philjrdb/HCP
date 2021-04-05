%% HumanCondPun questionnaires

ques_names = {'dass_anx' 'dass_dep' 'bis_sr' 'bis_imp' 'aor_rno' 'aor_rpo'...
  'bisbas_drv' 'bisbas_rew' 'bisbas_fun' 'bisbas_bis'...
  'ipip_e' 'ipip_a' 'ipip_c' 'ipip_n' 'ipip_i'}; % questionnaire names

ques{1} = [1 0 1 0 1 1 0 0 1 0 0 1 1];
ques{2} = [0 1 0 1 0 0 1 1 0 1 1 0 0];
ques{3} = [1 0 0 1 1 1 0 0];
ques{4} = [0 1 1 0 0 0 1 1];
ques{5} = [1 0 1 0 1 0 1 0 1 0 1 0];
ques{6} = [0 1 0 1 0 -1 0 -1 0 -1 0 1];
ques{7} = [1 0 0 0 1 0 1 0 0 0 0 0 0];
ques{8} = [0 1 1 0 0 0 0 0 0 0 0 1 0];
ques{9} = [0 0 0 0 0 1 0 0 1 0 1 0 0];
ques{10} = [0 0 0 1 0 0 0 1 0 1 0 0 1];
ques{11} = [1 0 0 0 0 -1 0 0 0 0 1 0 0 0 0 -1 0 0 0 0];
ques{12} = [0 1 0 0 0 0 -1 0 0 0 0 1 0 0 0 0 -1 0 0 0];
ques{13} = [0 0 1 0 0 0 0 -1 0 0 0 0 1 0 0 0 0 -1 0 0];
ques{14} = [0 0 0 1 0 0 0 0 -1 0 0 0 0 1 0 0 0 0 -1 0];
ques{15} = [0 0 0 0 1 0 0 0 0 -1 0 0 0 0 -1 0 0 0 0 -1];

a_leng = length(HCP_aggr);

%%
for q = 1:length(ques_names)  
  n_idx = regexp(ques_names{q},'_');
  ques_field = ['ques_' ques_names{q}(1:n_idx-1)];
  q_idx = find(ques{q} ~= 0);
  q_scoring = ques{q}(q_idx);
  
  for r = 1:a_leng
    HCP_aggr(r).(ques_names{q}) = mean(HCP_aggr(r).(ques_field)(q_idx).*q_scoring);
  end
end

clearvars -except HCP* work_folder *keep;