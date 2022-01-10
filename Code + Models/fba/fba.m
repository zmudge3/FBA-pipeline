% input folder name
input_folder = 'Test';

% signify if gene knockdown screen
% if True, will take all settings (besides gene knockdown) from first file in input folder, and use these for all subsequent files
gene_knockdown_screen = false;

% signify if objective function screen
% if True, will take all settings (besides objective function) from first file in input folder, and use these for all subsequent files
objective_function_screen = false;

% number of parallel cores
parallel_cores = 1;

% look for already-completed input files
look_for_completed = false;
completed_output_folder = 'RAD_knockout_nadh';

% don't do already-completed files
if look_for_completed
    ds = {sprintf('results/%s',completed_output_folder)};
    overlap = {};
    for i = 1:length(ds)
        directory = dir(ds{i});
        directory = directory(3:end);
        for j = 1:length(directory)
            overlap{end+1} = directory(j).name;
        end
    end
else
    overlap = {};
end

directory = dir(sprintf('input/%s/*.xlsx',input_folder));
input_files = {};
for i = 1:length(directory)
    fn = sprintf('input/%s/%s',input_folder,directory(i).name);
    test = strsplit(fn,'/');
    if ~strcmp(extractBetween(test{end},1,2),'~$')
        temp = strsplit(test{end},'.xlsx');
        if ~any(strcmp(overlap,temp{1}))
            input_files{end+1} = fn;
        end
    end
end

% name output folder
if ~(exist(sprintf('results/%s',input_folder)) == 7)
    output_folder = input_folder;
else
    num = 1;
    while exist(sprintf('results/%s',sprintf('%s_%d',input_folder,num))) == 7
        num = num + 1;
    end
    output_folder = sprintf('%s_%d',input_folder,num);
end

% create output folder
mkdir(sprintf('results/%s',output_folder))

tic

% if gene knockdown screen
if gene_knockdown_screen

    % parse info from first input file
    input_file = input_files(1);

    % parse input file - general
    [type_of_analysis,pfba_fraction,fva_all_fraction,fva_select_fraction,fva_select_groups,fva_select_reactions,sampler_fraction,sampler_samples,sampler_steps,media_sensitivity_critical_fraction,media_sensitivity_critical_metabolite1,media_sensitivity_spectrum_fraction,media_sensitivity_spectrum_metabolite1,media_sensitivity_spectrum_metabolite2,media_sensitivity_spectrum_points] = parse_general(input_file{1});

    % parse input file - constraints - proliferation
    [constraints_proliferation,constraints_proliferation_type,constraints_proliferation_value] = parse_constraints_proliferation(input_file{1});

    % parse input file - constraints - value
    [constraints_value_type,constraints_value_value,constraints_value_reaction,constraints_value_count] = parse_constraints_value(input_file{1});

    % parse input file - constraints - fraction
    [constraints_fraction_type,constraints_fraction_value,constraints_fraction_reaction,constraints_fraction_count] = parse_constraints_fraction(input_file{1});
    
    % parse input file - constraints - max multi
    [constraints_maxmulti_weight,constraints_maxmulti_reaction,constraints_maxmulti_fraction,constraints_maxmulti_count] = parse_constraints_maxmulti(input_file{1});
    
    % parse input file - objective function
    [objective_function_direction,objective_function_weight,objective_function_reaction] = parse_objective(input_file{1});

    % parse input file - drugs
    [required_drugs] = parse_drugs(input_file{1});

    % parse input file - mutations
    [use_mutation] = parse_mutations(input_file{1});

    % parse input file - samples - TCGA
    [samples_all,samples_all_source,samples_all_source_protein] = parse_samples_tcga(input_file{1},constraints_proliferation,required_drugs,use_mutation);
    
    % parse input file - samples - CCLE
    [samples_all,samples_all_source,samples_all_source_protein] = parse_samples_ccle(input_file{1},constraints_proliferation,required_drugs,use_mutation,samples_all,samples_all_source,samples_all_source_protein);
    
    % parse input file - samples - GTEx
    [samples_all,samples_all_source,samples_all_source_protein] = parse_samples_gtex(input_file{1},constraints_proliferation,required_drugs,use_mutation,samples_all,samples_all_source,samples_all_source_protein);
    
    % parse input file - samples - other
    [samples_all,samples_all_source,samples_all_source_protein] = parse_samples_other(input_file{1},constraints_proliferation,required_drugs,use_mutation,samples_all,samples_all_source,samples_all_source_protein);
    
    % make sure at least one sample chosen
    if isempty(samples_all)
        error('ERROR - Samples - Must choose at least one sample')
    end

    % parse input file - media
    [media_choice,media_constraint_metabolite,media_constraint_uptake] = parse_media(input_file{1});

    % parse input file - concentrations
    [recalculate_thermodynamics,default_concentration_lb,default_concentration_ub,concentration_ranges_all,concentration_ranges_sample_folder,concentration_ranges_sample] = parse_concentrations(input_file{1});
    
    % parse input file - modules
    [catalase_version, module_5fu, module_cis, module_cpa, module_dox] = parse_modules(input_file{1});
    
    % iterate over input files
    %parpool(parallel_cores);
    %parfor a = 1:length(input_files)
    for a = 1:length(input_files)
        input_file = input_files(a);
    
        % parse input file - gene knockdown
        [knockdown_genes,knockdown_fractions] = parse_knockdown(input_file{1});
        
        % load model
        [model,concentration_mets,concentration_values,constraints_fraction_id] = load_model(constraints_proliferation,constraints_proliferation_type,constraints_proliferation_value,constraints_value_type,constraints_value_value,constraints_value_reaction,constraints_value_count,constraints_fraction_type,constraints_fraction_value,constraints_fraction_reaction,constraints_fraction_count,constraints_maxmulti_weight,constraints_maxmulti_reaction,constraints_maxmulti_count,objective_function_direction,objective_function_weight,objective_function_reaction);
        
        % load media
        [model,concentration_mets,concentration_values,exchange_rxn] = load_media(model,media_choice,media_constraint_metabolite,media_constraint_uptake,concentration_mets,concentration_values,module_5fu,module_cis,module_cpa,module_dox);

        % load deltag values
        [deltag_lower,deltag_upper] = load_deltag(model,recalculate_thermodynamics,default_concentration_lb,default_concentration_ub,concentration_ranges_all,concentration_ranges_sample_folder,concentration_ranges_sample,concentration_mets,concentration_values,samples_all,samples_all_source);

        % calculate Vmax for each sample
        [vmax_forward,vmax_reverse] = calculate_vmax(model,samples_all,samples_all_source,samples_all_source_protein,use_mutation,knockdown_genes,knockdown_fractions,gene_knockdown_screen);

        % run FBA for each sample
        run_fba(output_folder,input_file{1},model,type_of_analysis,samples_all,samples_all_source,samples_all_source_protein,vmax_forward,vmax_reverse,objective_function_direction,objective_function_weight,objective_function_reaction,exchange_rxn,deltag_lower,deltag_upper,constraints_fraction_count,constraints_fraction_type,constraints_fraction_reaction,constraints_fraction_value,constraints_fraction_id,constraints_maxmulti_fraction,constraints_maxmulti_count,pfba_fraction,fva_all_fraction,fva_select_fraction,fva_select_groups,fva_select_reactions,sampler_fraction,sampler_samples,sampler_steps,media_sensitivity_critical_fraction,media_sensitivity_critical_metabolite1,media_sensitivity_spectrum_fraction,media_sensitivity_spectrum_metabolite1,media_sensitivity_spectrum_metabolite2,media_sensitivity_spectrum_points,parallel_cores,knockdown_genes,knockdown_fractions,catalase_version,module_5fu,module_cis,module_cpa,module_dox);
    end

% if objective function screen
elseif objective_function_screen

    % parse info from first input file
    input_file = input_files(1);

    % parse input file - general
    [type_of_analysis,pfba_fraction,fva_all_fraction,fva_select_fraction,fva_select_groups,fva_select_reactions,sampler_fraction,sampler_samples,sampler_steps,media_sensitivity_critical_fraction,media_sensitivity_critical_metabolite1,media_sensitivity_spectrum_fraction,media_sensitivity_spectrum_metabolite1,media_sensitivity_spectrum_metabolite2,media_sensitivity_spectrum_points] = parse_general(input_file{1});

    % parse input file - constraints - proliferation
    [constraints_proliferation,constraints_proliferation_type,constraints_proliferation_value] = parse_constraints_proliferation(input_file{1});

    % parse input file - constraints - value
    [constraints_value_type,constraints_value_value,constraints_value_reaction,constraints_value_count] = parse_constraints_value(input_file{1});

    % parse input file - constraints - fraction
    [constraints_fraction_type,constraints_fraction_value,constraints_fraction_reaction,constraints_fraction_count] = parse_constraints_fraction(input_file{1});
    
    % parse input file - constraints - max multi
    [constraints_maxmulti_weight,constraints_maxmulti_reaction,constraints_maxmulti_fraction,constraints_maxmulti_count] = parse_constraints_maxmulti(input_file{1});

    % parse input file - drugs
    [required_drugs] = parse_drugs(input_file{1});

    % parse input file - mutations
    [use_mutation] = parse_mutations(input_file{1});

    % parse input file - samples - TCGA
    [samples_all,samples_all_source,samples_all_source_protein] = parse_samples_tcga(input_file{1},constraints_proliferation,required_drugs,use_mutation);
    
    % parse input file - samples - CCLE
    [samples_all,samples_all_source,samples_all_source_protein] = parse_samples_ccle(input_file{1},constraints_proliferation,required_drugs,use_mutation,samples_all,samples_all_source,samples_all_source_protein);
    
    % parse input file - samples - GTEx
    [samples_all,samples_all_source,samples_all_source_protein] = parse_samples_gtex(input_file{1},constraints_proliferation,required_drugs,use_mutation,samples_all,samples_all_source,samples_all_source_protein);
    
    % parse input file - samples - other
    [samples_all,samples_all_source,samples_all_source_protein] = parse_samples_other(input_file{1},constraints_proliferation,required_drugs,use_mutation,samples_all,samples_all_source,samples_all_source_protein);
    
    % make sure at least one sample chosen
    if isempty(samples_all)
        error('ERROR - Samples - Must choose at least one sample')
    end

    % parse input file - media
    [media_choice,media_constraint_metabolite,media_constraint_uptake] = parse_media(input_file{1});

    % parse input file - concentrations
    [recalculate_thermodynamics,default_concentration_lb,default_concentration_ub,concentration_ranges_all,concentration_ranges_sample_folder,concentration_ranges_sample] = parse_concentrations(input_file{1});
    
    % parse input file - modules
    [catalase_version, module_5fu, module_cis, module_cpa, module_dox] = parse_modules(input_file{1});
    
    % parse input file - gene knockdown
    [knockdown_genes,knockdown_fractions] = parse_knockdown(input_file{1});
    
    % iterate over input files
    %parpool(parallel_cores);
    %parfor a = 1:length(input_files)
    for a = 1:length(input_files)
        input_file = input_files(a);
    
        % parse input file - objective function
        [objective_function_direction,objective_function_weight,objective_function_reaction] = parse_objective(input_file{1});
        
        % load model
        [model,concentration_mets,concentration_values,constraints_fraction_id] = load_model(constraints_proliferation,constraints_proliferation_type,constraints_proliferation_value,constraints_value_type,constraints_value_value,constraints_value_reaction,constraints_value_count,constraints_fraction_type,constraints_fraction_value,constraints_fraction_reaction,constraints_fraction_count,constraints_maxmulti_weight,constraints_maxmulti_reaction,constraints_maxmulti_count,objective_function_direction,objective_function_weight,objective_function_reaction);
        
        % load media
        [model,concentration_mets,concentration_values,exchange_rxn] = load_media(model,media_choice,media_constraint_metabolite,media_constraint_uptake,concentration_mets,concentration_values,module_5fu,module_cis,module_cpa,module_dox);

        % load deltag values
        [deltag_lower,deltag_upper] = load_deltag(model,recalculate_thermodynamics,default_concentration_lb,default_concentration_ub,concentration_ranges_all,concentration_ranges_sample_folder,concentration_ranges_sample,concentration_mets,concentration_values,samples_all,samples_all_source);

        % calculate Vmax for each sample
        [vmax_forward,vmax_reverse] = calculate_vmax(model,samples_all,samples_all_source,samples_all_source_protein,use_mutation,knockdown_genes,knockdown_fractions,gene_knockdown_screen);

        % run FBA for each sample
        run_fba(output_folder,input_file{1},model,type_of_analysis,samples_all,samples_all_source,samples_all_source_protein,vmax_forward,vmax_reverse,objective_function_direction,objective_function_weight,objective_function_reaction,exchange_rxn,deltag_lower,deltag_upper,constraints_fraction_count,constraints_fraction_type,constraints_fraction_reaction,constraints_fraction_value,constraints_fraction_id,constraints_maxmulti_fraction,constraints_maxmulti_count,pfba_fraction,fva_all_fraction,fva_select_fraction,fva_select_groups,fva_select_reactions,sampler_fraction,sampler_samples,sampler_steps,media_sensitivity_critical_fraction,media_sensitivity_critical_metabolite1,media_sensitivity_spectrum_fraction,media_sensitivity_spectrum_metabolite1,media_sensitivity_spectrum_metabolite2,media_sensitivity_spectrum_points,parallel_cores,knockdown_genes,knockdown_fractions,catalase_version,module_5fu,module_cis,module_cpa,module_dox);
    end

% if not either
else
        
    % iterate over input files
    %parpool(parallel_cores);
    %parfor a = 1:length(input_files)
    for a = 1:length(input_files)
        input_file = input_files(a);

        % parse input file - general
        [type_of_analysis,pfba_fraction,fva_all_fraction,fva_select_fraction,fva_select_groups,fva_select_reactions,sampler_fraction,sampler_samples,sampler_steps,media_sensitivity_critical_fraction,media_sensitivity_critical_metabolite1,media_sensitivity_spectrum_fraction,media_sensitivity_spectrum_metabolite1,media_sensitivity_spectrum_metabolite2,media_sensitivity_spectrum_points] = parse_general(input_file{1});

        % parse input file - constraints - proliferation
        [constraints_proliferation,constraints_proliferation_type,constraints_proliferation_value] = parse_constraints_proliferation(input_file{1});

        % parse input file - constraints - value
        [constraints_value_type,constraints_value_value,constraints_value_reaction,constraints_value_count] = parse_constraints_value(input_file{1});

        % parse input file - constraints - fraction
        [constraints_fraction_type,constraints_fraction_value,constraints_fraction_reaction,constraints_fraction_count] = parse_constraints_fraction(input_file{1});
        
        % parse input file - constraints - max multi
        [constraints_maxmulti_weight,constraints_maxmulti_reaction,constraints_maxmulti_fraction,constraints_maxmulti_count] = parse_constraints_maxmulti(input_file{1});

        % parse input file - objective function
        [objective_function_direction,objective_function_weight,objective_function_reaction] = parse_objective(input_file{1});
        
        % parse input file - drugs
        [required_drugs] = parse_drugs(input_file{1});

        % parse input file - mutations
        [use_mutation] = parse_mutations(input_file{1});

        % parse input file - samples - TCGA
        [samples_all,samples_all_source,samples_all_source_protein] = parse_samples_tcga(input_file{1},constraints_proliferation,required_drugs,use_mutation);
        
        % parse input file - samples - CCLE
        [samples_all,samples_all_source,samples_all_source_protein] = parse_samples_ccle(input_file{1},constraints_proliferation,required_drugs,use_mutation,samples_all,samples_all_source,samples_all_source_protein);
        
        % parse input file - samples - GTEx
        [samples_all,samples_all_source,samples_all_source_protein] = parse_samples_gtex(input_file{1},constraints_proliferation,required_drugs,use_mutation,samples_all,samples_all_source,samples_all_source_protein);
         
        % parse input file - samples - other
        [samples_all,samples_all_source,samples_all_source_protein] = parse_samples_other(input_file{1},constraints_proliferation,required_drugs,use_mutation,samples_all,samples_all_source,samples_all_source_protein);
    
        % make sure at least one sample chosen
        if isempty(samples_all)
            error('ERROR - Samples - Must choose at least one sample')
        end

        % parse input file - media
        [media_choice,media_constraint_metabolite,media_constraint_uptake] = parse_media(input_file{1});

        % parse input file - concentrations
        [recalculate_thermodynamics,default_concentration_lb,default_concentration_ub,concentration_ranges_all,concentration_ranges_sample_folder,concentration_ranges_sample] = parse_concentrations(input_file{1});

        % parse input file - modules
        [catalase_version, module_5fu, module_cis, module_cpa, module_dox] = parse_modules(input_file{1});

        % parse input file - gene knockdown
        [knockdown_genes,knockdown_fractions] = parse_knockdown(input_file{1});

        % load model
        [model,concentration_mets,concentration_values,constraints_fraction_id] = load_model(constraints_proliferation,constraints_proliferation_type,constraints_proliferation_value,constraints_value_type,constraints_value_value,constraints_value_reaction,constraints_value_count,constraints_fraction_type,constraints_fraction_value,constraints_fraction_reaction,constraints_fraction_count,constraints_maxmulti_weight,constraints_maxmulti_reaction,constraints_maxmulti_count,objective_function_direction,objective_function_weight,objective_function_reaction);
        
        % load media
        [model,concentration_mets,concentration_values,exchange_rxn] = load_media(model,media_choice,media_constraint_metabolite,media_constraint_uptake,concentration_mets,concentration_values,module_5fu,module_cis,module_cpa,module_dox);

        % load deltag values
        [deltag_lower,deltag_upper] = load_deltag(model,recalculate_thermodynamics,default_concentration_lb,default_concentration_ub,concentration_ranges_all,concentration_ranges_sample_folder,concentration_ranges_sample,concentration_mets,concentration_values,samples_all,samples_all_source);

        % calculate Vmax for each sample
        [vmax_forward,vmax_reverse] = calculate_vmax(model,samples_all,samples_all_source,samples_all_source_protein,use_mutation,knockdown_genes,knockdown_fractions,gene_knockdown_screen);

        % run FBA for each sample
        run_fba(output_folder,input_file{1},model,type_of_analysis,samples_all,samples_all_source,samples_all_source_protein,vmax_forward,vmax_reverse,objective_function_direction,objective_function_weight,objective_function_reaction,exchange_rxn,deltag_lower,deltag_upper,constraints_fraction_count,constraints_fraction_type,constraints_fraction_reaction,constraints_fraction_value,constraints_fraction_id,constraints_maxmulti_fraction,constraints_maxmulti_count,pfba_fraction,fva_all_fraction,fva_select_fraction,fva_select_groups,fva_select_reactions,sampler_fraction,sampler_samples,sampler_steps,media_sensitivity_critical_fraction,media_sensitivity_critical_metabolite1,media_sensitivity_spectrum_fraction,media_sensitivity_spectrum_metabolite1,media_sensitivity_spectrum_metabolite2,media_sensitivity_spectrum_points,parallel_cores,knockdown_genes,knockdown_fractions,catalase_version,module_5fu,module_cis,module_cpa,module_dox);
    end    
end

toc


