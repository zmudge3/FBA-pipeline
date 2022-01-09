tic

source = 'sc-glio'

parallel_cores = 1;

load('../recon/recon3d_qflux.mat');

parsedGPR = GPRparser(model);

idh1_patients_H133Q = importdata('../idh1/H133Q.txt');
idh1_patients_A134D = importdata('../idh1/A134D.txt');
idh1_patients_R100Q = importdata('../idh1/R100Q.txt');
idh1_patients_R132H = importdata('../idh1/R132H.txt');
idh1_patients_R132C = importdata('../idh1/R132C.txt');
idh1_patients_R132G = importdata('../idh1/R132G.txt');
idh1_patients_R132W = importdata('../idh1/R132W.txt');
idh1_patients_R132A = importdata('../idh1/R132A.txt');
idh1_patients_R132Q = importdata('../idh1/R132Q.txt');
idh1_patients_R132K = importdata('../idh1/R132K.txt');
idh1_patients_R132N = importdata('../idh1/R132N.txt');
idh1_patients_all = [idh1_patients_H133Q;idh1_patients_A134D;idh1_patients_R100Q;idh1_patients_R132H;idh1_patients_R132C;idh1_patients_R132G;idh1_patients_R132W;idh1_patients_R132A;idh1_patients_R132Q;idh1_patients_R132K;idh1_patients_R132N];

% forward reaction
idh1_normal_H133Q = 45*60*60;
idh1_normal_R132A = 10.4*60*60;
idh1_normal_R132G = 9.3*60*60;
idh1_normal_R132Q = 9.2*60*60;
idh1_normal_R132K = 7.2*60*60;
idh1_normal_R100Q = 5.6*60*60;
idh1_normal_R132C = 4.4*60*60;
idh1_normal_R132H = 2.4*60*60;
idh1_normal_A134D = 2.3*60*60;
idh1_normal_R132W = 1.21*60*60;
idh1_normal_R132N = 0.047*60*60;

% neomorphic reaction
idh1_neomorphic_H133Q = 0*60*60;
idh1_neomorphic_A134D = 0*60*60;
idh1_neomorphic_R100Q = 0.34*60*60;
idh1_neomorphic_R132A = 0.37*60*60;
idh1_neomorphic_R132W = 0.54*60*60;
idh1_neomorphic_R132K = 0.57*60*60;
idh1_neomorphic_R132N = 0.79*60*60;
idh1_neomorphic_R132G = 1.59*60*60;
idh1_neomorphic_R132C = 1.60*60*60;
idh1_neomorphic_R132H = 4.2*60*60;
idh1_neomorphic_R132Q = 4.7*60*60;

mkdir(sprintf('%s',source))
mkdir(sprintf('%s/no_mutation',source))
mkdir(sprintf('%s/yes_mutation',source))

% get all data files
directory = dir(sprintf('../protein/output/%s/*.csv',source));

% iterate over data files
delete(gcp('nocreate'));
parpool(parallel_cores);
parfor a = 1:length(directory)
%for a = 1:length(directory)

%%% no mutation data

    % load kcat data
    % forward
    f = fopen('../kcat/kcat_forward.csv','r');
    kcat_forward = textscan(f,'%s %f %s','Delimiter',',','headerLines',1);
    fclose(f);

    % reverse
    f = fopen('../kcat/kcat_reverse.csv','r');
    kcat_reverse = textscan(f,'%s %f %s','Delimiter',',','headerLines',1);
    fclose(f);

    % load protein data
    f = fopen(sprintf('../protein/output/%s/%s',source,directory(a).name),'r');
    protein = textscan(f,'%s %f','Delimiter',',');
    protein{1} = strrep(protein{1},'"','');
    fclose(f);

    % load compartment gene fractions
    f = fopen('_data_/compartments/gene_fraction.tsv','r');
    gene_fraction = textscan(f,'%s %s %s %f','Delimiter','\t','headerLines',1);
    fclose(f);

    % convert ppm to mmol/gDW
    m = 50; % average mass of protein, in kDa
    fdw = 0.5; % fraction of dry weight that is protein
    Na = 6.0221409*10^23; % Avogadro constant
    gDa = 1.66054*10^-24; % number of grams per dalton

    protein{2} = protein{2}/1000000; % protein/million --> protein/protein
    protein{2} = protein{2} / Na; % protein/protein --> mol/protein
    protein{2} = protein{2} * 1000; % mol/protein --> mmol/protein
    protein{2} = protein{2} / (m*1000); % mmol/protein --> mmol/Daltons_protein
    protein{2} = protein{2} / gDa; % mmol/Daltons --> mmol/grams_protein
    protein{2} = protein{2} * fdw; % mmol/grams_protein --> mmol/gDW

    % get protein expression for each reaction with grRule and EC number
    protein_rxn = [];
    for m = 1:length(model.rxns)
        if (~isempty(model.rules{m})) && (~isempty(model.rxnECNumbers{m}))

            % get genes associated with reaction
            genes = find(model.rxnGeneMat(m,:));
            allfound = 1;
            for n = 1:length(genes)
                if ~any(strcmp(protein{1},model.geneSymbols{genes(n)}));
                    allfound = 0;
                    break
                end
            end
            if allfound == 1

                % create model structure for selectGeneFromGPR
                modelRxn = {};
                modelRxn.rxns = {model.rxns{m}};
                modelRxn.genes = {};
                for n = 1:length(genes)
                    modelRxn.genes{end+1} = model.genes{genes(n)};
                end
                modelRxn.genes = modelRxn.genes';

                % calculate reaction expression
                values = [];
                for n = 1:length(genes)
                    values(end+1) = protein{2}(strcmp(protein{1},model.geneSymbols{genes(n)})) * gene_fraction{4}(strcmp(gene_fraction{1},model.rxns{m}) & strcmp(gene_fraction{3},model.geneSymbols{genes(n)}));
                end
                % updated version of selectGeneFromGPR
                protein_rxn(end+1) = selectGeneFromGPR(modelRxn,modelRxn.genes,values',{parsedGPR{m}});

            else
               protein_rxn(end+1) = NaN; 
            end

        else
            protein_rxn(end+1) = NaN;
        end
    end

    % set vmax values
    vmax_forward_current = nan(length(model.rxns),1);
    vmax_reverse_current = nan(length(model.rxns),1);
    for m = 1:length(model.rxns)

        % forward vmax
        if model.ub(m) > 0
            if (~(isnan(protein_rxn(m)))) && (~(isnan(kcat_forward{2}(strcmp(kcat_forward{1},model.rxns{m})))))
                vmax_forward_current(m) = protein_rxn(m) * kcat_forward{2}(strcmp(kcat_forward{1},model.rxns{m}));
            end
        end

        % reverse vmax
        if model.lb(m) < 0
            if (~(isnan(protein_rxn(m)))) && (~(isnan(kcat_reverse{2}(strcmp(kcat_reverse{1},model.rxns{m})))))
                vmax_reverse_current(m) = protein_rxn(m) * kcat_reverse{2}(strcmp(kcat_reverse{1},model.rxns{m}));
            end
        end
    end

    % save results
    f = fopen(sprintf('%s/no_mutation/%s',source,directory(a).name),'w');
    fprintf(f,'REACTION,FORWARD [mmol/gDW/hr],REVERSE [mmol/gDW/hr]\n');
    for m = 1:length(model.rxns)
        fprintf(f,'%s,%.10f,%.10f\n',model.rxns{m},vmax_forward_current(m),vmax_reverse_current(m));
    end
    fclose(f);

%%% yes mutation data

    % if mutation data exists
    temp = strsplit(source,'_');
    prefix = temp{1};
    switch prefix
        case 'TCGA'
            lookup_name = 'TCGA';
        case 'NCI60'
            lookup_name = 'NCI60';
        otherwise
            lookup_name = source;
    end
    temp = strsplit(directory(a).name,'.');
    sample_name = temp{1};
    if (exist(sprintf('../mutation/%s/%s',lookup_name,directory(a).name),'file')) || (any(strcmp(idh1_patients_all,sample_name)))

        % load protein data
        f = fopen(sprintf('../protein/output/%s/%s',source,directory(a).name),'r');
        protein = textscan(f,'%s %f','Delimiter',',');
        protein{1} = strrep(protein{1},'"','');
        fclose(f);

        % load mutation data
        if (exist(sprintf('../mutation/%s/%s',lookup_name,directory(a).name),'file'))
            f = fopen(sprintf('../mutation/%s/%s',lookup_name,directory(a).name),'r');
            mutation = textscan(f,'%s %f','Delimiter',',','headerLines',1);
            fclose(f);

            % multiply protein expression data by envision scores
            for m = 1:length(mutation{1})
                if any(strcmp(protein{1},mutation{1}{m}))
                    protein{2}(find(strcmp(protein{1},mutation{1}{m}))) = protein{2}(find(strcmp(protein{1},mutation{1}{m}))) * mutation{2}(m);
                end
            end
        end

        % convert ppm to mmol/gDW
        m = 50; % average mass of protein, in kDa
        fdw = 0.5; % fraction of dry weight that is protein
        Na = 6.0221409*10^23; % Avogadro constant
        gDa = 1.66054*10^-24; % number of grams per dalton

        protein{2} = protein{2}/1000000; % protein/million --> protein/protein
        protein{2} = protein{2} / Na; % protein/protein --> mol/protein
        protein{2} = protein{2} * 1000; % mol/protein --> mmol/protein
        protein{2} = protein{2} / (m*1000); % mmol/protein --> mmol/Daltons_protein
        protein{2} = protein{2} / gDa; % mmol/Daltons --> mmol/grams_protein
        protein{2} = protein{2} * fdw; % mmol/grams_protein --> mmol/gDW

        % get protein expression for each reaction with grRule and EC number
        protein_rxn = [];
        for m = 1:length(model.rxns)
            if (~isempty(model.rules{m})) && (~isempty(model.rxnECNumbers{m}))

                % get genes associated with reaction
                genes = find(model.rxnGeneMat(m,:));
                allfound = 1;
                for n = 1:length(genes)
                    if ~any(strcmp(protein{1},model.geneSymbols{genes(n)}));
                        allfound = 0;
                        break
                    end
                end
                if allfound == 1

                    % create model structure for selectGeneFromGPR
                    modelRxn = {};
                    modelRxn.rxns = {model.rxns{m}};
                    modelRxn.genes = {};
                    for n = 1:length(genes)
                        modelRxn.genes{end+1} = model.genes{genes(n)};
                    end
                    modelRxn.genes = modelRxn.genes';

                    % calculate reaction expression
                    values = [];
                    for n = 1:length(genes)
                        values(end+1) = protein{2}(strcmp(protein{1},model.geneSymbols{genes(n)})) * gene_fraction{4}(strcmp(gene_fraction{1},model.rxns{m}) & strcmp(gene_fraction{3},model.geneSymbols{genes(n)}));
                    end
                    % updated version of selectGeneFromGPR
                    protein_rxn(end+1) = selectGeneFromGPR(modelRxn,modelRxn.genes,values',{parsedGPR{m}});

                else
                   protein_rxn(end+1) = NaN; 
                end

            else
                protein_rxn(end+1) = NaN;
            end
        end

        % IDH1 mutation kcat's
        if (any(strcmp(idh1_patients_all,sample_name)))
            if any(strcmp(idh1_patients_H133Q,sample_name))
                kcat_forward{2}(strcmp(kcat_forward{1},'ICDHy')) = idh1_normal_H133Q;
                kcat_forward{2}(strcmp(kcat_forward{1},'ICDHyp')) = idh1_normal_H133Q;
                kcat_forward{2}(strcmp(kcat_forward{1},'IDH1_R132')) = idh1_neomorphic_H133Q;
            elseif any(strcmp(idh1_patients_A134D,sample_name))
                kcat_forward{2}(strcmp(kcat_forward{1},'ICDHy')) = idh1_normal_A134D;
                kcat_forward{2}(strcmp(kcat_forward{1},'ICDHyp')) = idh1_normal_A134D;
                kcat_forward{2}(strcmp(kcat_forward{1},'IDH1_R132')) = idh1_neomorphic_A134D;
            elseif any(strcmp(idh1_patients_R100Q,sample_name))
                kcat_forward{2}(strcmp(kcat_forward{1},'ICDHy')) = idh1_normal_R100Q;
                kcat_forward{2}(strcmp(kcat_forward{1},'ICDHyp')) = idh1_normal_R100Q;
                kcat_forward{2}(strcmp(kcat_forward{1},'IDH1_R132')) = idh1_neomorphic_R100Q;
            elseif any(strcmp(idh1_patients_R132H,sample_name))
                kcat_forward{2}(strcmp(kcat_forward{1},'ICDHy')) = idh1_normal_R132H;
                kcat_forward{2}(strcmp(kcat_forward{1},'ICDHyp')) = idh1_normal_R132H;
                kcat_forward{2}(strcmp(kcat_forward{1},'IDH1_R132')) = idh1_neomorphic_R132H;
            elseif any(strcmp(idh1_patients_R132C,sample_name))
                kcat_forward{2}(strcmp(kcat_forward{1},'ICDHy')) = idh1_normal_R132C;
                kcat_forward{2}(strcmp(kcat_forward{1},'ICDHyp')) = idh1_normal_R132C;
                kcat_forward{2}(strcmp(kcat_forward{1},'IDH1_R132')) = idh1_neomorphic_R132C;
            elseif any(strcmp(idh1_patients_R132G,sample_name))
                kcat_forward{2}(strcmp(kcat_forward{1},'ICDHy')) = idh1_normal_R132G;
                kcat_forward{2}(strcmp(kcat_forward{1},'ICDHyp')) = idh1_normal_R132G;
                kcat_forward{2}(strcmp(kcat_forward{1},'IDH1_R132')) = idh1_neomorphic_R132G;
            elseif any(strcmp(idh1_patients_R132W,sample_name))
                kcat_forward{2}(strcmp(kcat_forward{1},'ICDHy')) = idh1_normal_R132W;
                kcat_forward{2}(strcmp(kcat_forward{1},'ICDHyp')) = idh1_normal_R132W;
                kcat_forward{2}(strcmp(kcat_forward{1},'IDH1_R132')) = idh1_neomorphic_R132W;
            elseif any(strcmp(idh1_patients_R132A,sample_name))
                kcat_forward{2}(strcmp(kcat_forward{1},'ICDHy')) = idh1_normal_R132A;
                kcat_forward{2}(strcmp(kcat_forward{1},'ICDHyp')) = idh1_normal_R132A;
                kcat_forward{2}(strcmp(kcat_forward{1},'IDH1_R132')) = idh1_neomorphic_R132A;
            elseif any(strcmp(idh1_patients_R132Q,sample_name))
                kcat_forward{2}(strcmp(kcat_forward{1},'ICDHy')) = idh1_normal_R132Q;
                kcat_forward{2}(strcmp(kcat_forward{1},'ICDHyp')) = idh1_normal_R132Q;
                kcat_forward{2}(strcmp(kcat_forward{1},'IDH1_R132')) = idh1_neomorphic_R132Q;
            elseif any(strcmp(idh1_patients_R132K,sample_name))
                kcat_forward{2}(strcmp(kcat_forward{1},'ICDHy')) = idh1_normal_R132K;
                kcat_forward{2}(strcmp(kcat_forward{1},'ICDHyp')) = idh1_normal_R132K;
                kcat_forward{2}(strcmp(kcat_forward{1},'IDH1_R132')) = idh1_neomorphic_R132K;
            elseif any(strcmp(idh1_patients_R132N,sample_name))
                kcat_forward{2}(strcmp(kcat_forward{1},'ICDHy')) = idh1_normal_R132N;
                kcat_forward{2}(strcmp(kcat_forward{1},'ICDHyp')) = idh1_normal_R132N;
                kcat_forward{2}(strcmp(kcat_forward{1},'IDH1_R132')) = idh1_neomorphic_R132N;
            end
        end

        % set vmax values
        vmax_forward_current = nan(length(model.rxns),1);
        vmax_reverse_current = nan(length(model.rxns),1);
        for m = 1:length(model.rxns)

            % forward vmax
            if model.ub(m) > 0
                if (~(isnan(protein_rxn(m)))) && (~(isnan(kcat_forward{2}(strcmp(kcat_forward{1},model.rxns{m})))))
                    vmax_forward_current(m) = protein_rxn(m) * kcat_forward{2}(strcmp(kcat_forward{1},model.rxns{m}));
                end
            end

            % reverse vmax
            if model.lb(m) < 0
                if (~(isnan(protein_rxn(m)))) && (~(isnan(kcat_reverse{2}(strcmp(kcat_reverse{1},model.rxns{m})))))
                    vmax_reverse_current(m) = protein_rxn(m) * kcat_reverse{2}(strcmp(kcat_reverse{1},model.rxns{m}));
                end
            end
        end

        % save results
        f = fopen(sprintf('%s/yes_mutation/%s',source,directory(a).name),'w');
        fprintf(f,'REACTION,FORWARD [mmol/gDW/hr],REVERSE [mmol/gDW/hr]\n');
        for m = 1:length(model.rxns)
            fprintf(f,'%s,%.10f,%.10f\n',model.rxns{m},vmax_forward_current(m),vmax_reverse_current(m));
        end
        fclose(f);
    end
end

toc
