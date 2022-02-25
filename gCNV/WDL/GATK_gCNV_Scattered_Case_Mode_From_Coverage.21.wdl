#
# A wrapper for the gCNV case workflow intended for lowering computing cost by making it feasible to use 
# preemptible cloud instances with low memory requirements. CPU, memory and disk requirements can be 
# lowered for GermlineCNVCaller and DetermineGermlineContigPloidy tasks. 
# 
#
# - Example invocation:
#
#       java -jar cromwell.jar run cnv_germline_case_scattered_workflow.wdl -i my_parameters.json
#
####################

import "https://api.firecloud.org/ga4gh/v1/tools/asmirnov-broad:GATK_gCNV_Case_Mode_From_Coverage/versions/18/plain-WDL/descriptor" as GermlineCNVCaseWorkflow

workflow CNVGermlineCaseScatteredWorkflow {

    ##################################
    #### required basic arguments ####
    ##################################
    File preprocessed_intervals
    Array[String]+ count_files
    Array[String]+ entity_ids
    File contig_ploidy_model_tar

	File file_gcnv_model_tars
    File file_calling_configs
    File file_denoising_configs
    File file_gcnvkernel_version
    File file_sharded_interval_lists
    
    Int num_intervals_per_scatter
    File ref_fasta_dict
    File ref_fasta_fai
    File ref_fasta
    String gatk_docker
    Int num_samples_per_scatter = 50

    ##################################
    #### optional basic arguments ####
    ##################################
    File? gatk4_jar_override
    Int? preemptible_attempts


    ######################################################################
    #### optional arguments for DetermineGermlineContigPloidyCaseMode ####
    ######################################################################
    Float? ploidy_mapping_error_rate
    Float? ploidy_sample_psi_scale
    Int? mem_gb_for_determine_germline_contig_ploidy
    Int? cpu_for_determine_germline_contig_ploidy
    Int? disc_for_determine_germline_contig_ploidy    

    ##########################################################
    #### optional arguments for GermlineCNVCallerCaseMode ####
    ##########################################################
    Float? gcnv_p_alt
    Float? gcnv_cnv_coherence_length
    Int? gcnv_max_copy_number
    Int? mem_gb_for_germline_cnv_caller
    Int? cpu_for_germline_cnv_caller
    Int? disc_for_germline_cnv_caller

    # optional arguments for germline CNV denoising model
    Float? gcnv_mapping_error_rate
    Float? gcnv_sample_psi_scale
    Float? gcnv_depth_correction_tau
    String? gcnv_copy_number_posterior_expectation_mode
    Int? gcnv_active_class_padding_hybrid_mode

    # optional arguments for Hybrid ADVI
    Float? gcnv_learning_rate
    Float? gcnv_adamax_beta_1
    Float? gcnv_adamax_beta_2
    Int? gcnv_log_emission_samples_per_round
    Float? gcnv_log_emission_sampling_median_rel_error
    Int? gcnv_log_emission_sampling_rounds
    Int? gcnv_max_advi_iter_first_epoch
    Int? gcnv_max_advi_iter_subsequent_epochs
    Int? gcnv_min_training_epochs
    Int? gcnv_max_training_epochs
    Float? gcnv_initial_temperature
    Int? gcnv_num_thermal_advi_iters
    Int? gcnv_convergence_snr_averaging_window
    Float? gcnv_convergence_snr_trigger_threshold
    Int? gcnv_convergence_snr_countdown_window
    Int? gcnv_max_calling_iters
    Float? gcnv_caller_update_convergence_threshold
    Float? gcnv_caller_internal_admixing_rate
    Float? gcnv_caller_external_admixing_rate
    Boolean? gcnv_disable_annealing

    ###################################################
    #### arguments for PostprocessGermlineCNVCalls ####
    ###################################################
    Int ref_copy_number_autosomal_contigs
    Array[String]? allosomal_contigs
    
    ###################################################
    call fof_usage_task as fof_usage_task_gcnv_model_tars{
    	input:
        	fof = file_gcnv_model_tars
    }

    call fof_usage_task as fof_usage_task_calling_configs{
    	input:
        	fof = file_calling_configs
    }

    call fof_usage_task as fof_usage_task_denoising_configs{
    	input:
        	fof = file_denoising_configs
    }

    call fof_usage_task as fof_usage_task_gcnvkernel_version{
    	input:
        	fof = file_gcnvkernel_version
    }

    call fof_usage_task as fof_usage_task_sharded_interval_lists{
    	input:
        	fof = file_sharded_interval_lists
    }
    Array[String] gcnv_model_tars = fof_usage_task_gcnv_model_tars.array_of_files
    Array[String] calling_configs = fof_usage_task_calling_configs.array_of_files
    Array[String] denoising_configs = fof_usage_task_denoising_configs.array_of_files
    Array[String] gcnvkernel_version = fof_usage_task_gcnvkernel_version.array_of_files
    Array[String] sharded_interval_lists = fof_usage_task_sharded_interval_lists.array_of_files
    
   ###################################################

    call SplitInputArray as SplitInputCountsList {
        input:
            input_array = count_files,
            num_inputs_in_scatter = num_samples_per_scatter,
            gatk_docker = gatk_docker
    }
    
    call SplitInputArray as SplitInputIdList {
        input:
            input_array = entity_ids,
            num_inputs_in_scatter = num_samples_per_scatter,
            gatk_docker = gatk_docker
    }

    Array[Array[File]] splitCountsInput = SplitInputCountsList.splitArray
    Array[Array[String]] splitIdsInput = SplitInputIdList.splitArray

    scatter (subarray_index in range(length(splitCountsInput))) {
        call GermlineCNVCaseWorkflow.CNVGermlineCaseWorkflowFromCoverage {
            input:
                preprocessed_intervals = preprocessed_intervals,
                count_files = splitCountsInput[subarray_index],
                entity_ids = splitIdsInput[subarray_index],
                contig_ploidy_model_tar = contig_ploidy_model_tar,
                gcnv_model_tars = gcnv_model_tars,
                calling_configs = calling_configs,
                denoising_configs = denoising_configs,
                gcnvkernel_version = gcnvkernel_version,
                sharded_interval_lists = sharded_interval_lists,
                num_intervals_per_scatter = num_intervals_per_scatter,
                ref_fasta_dict = ref_fasta_dict,
                ref_fasta_fai = ref_fasta_fai,
                ref_fasta = ref_fasta,
                gatk_docker = gatk_docker,
                gatk4_jar_override = gatk4_jar_override,
                preemptible_attempts = preemptible_attempts,
                ploidy_mapping_error_rate = ploidy_mapping_error_rate,
                ploidy_sample_psi_scale = ploidy_sample_psi_scale,
                mem_gb_for_determine_germline_contig_ploidy = mem_gb_for_determine_germline_contig_ploidy,
                cpu_for_determine_germline_contig_ploidy = cpu_for_determine_germline_contig_ploidy,
                disc_for_determine_germline_contig_ploidy = disc_for_determine_germline_contig_ploidy,
                gcnv_p_alt = gcnv_p_alt,
                gcnv_cnv_coherence_length = gcnv_cnv_coherence_length,
                gcnv_max_copy_number = gcnv_max_copy_number,
                mem_gb_for_germline_cnv_caller = mem_gb_for_germline_cnv_caller,
                cpu_for_germline_cnv_caller = cpu_for_germline_cnv_caller,
                disc_for_germline_cnv_caller = disc_for_germline_cnv_caller,
                gcnv_mapping_error_rate = gcnv_mapping_error_rate,
                gcnv_sample_psi_scale = gcnv_sample_psi_scale,
                gcnv_depth_correction_tau = gcnv_depth_correction_tau,
                gcnv_copy_number_posterior_expectation_mode = gcnv_copy_number_posterior_expectation_mode,
                gcnv_active_class_padding_hybrid_mode = gcnv_active_class_padding_hybrid_mode,
                gcnv_learning_rate = gcnv_learning_rate,
                gcnv_adamax_beta_1 = gcnv_adamax_beta_1,
                gcnv_adamax_beta_2 = gcnv_adamax_beta_2,
                gcnv_log_emission_samples_per_round = gcnv_log_emission_samples_per_round,
                gcnv_log_emission_sampling_median_rel_error = gcnv_log_emission_sampling_median_rel_error,
                gcnv_log_emission_sampling_rounds = gcnv_log_emission_sampling_rounds,
                gcnv_max_advi_iter_first_epoch = gcnv_max_advi_iter_first_epoch,
                gcnv_max_advi_iter_subsequent_epochs = gcnv_max_advi_iter_subsequent_epochs,
                gcnv_min_training_epochs = gcnv_min_training_epochs,
                gcnv_max_training_epochs = gcnv_max_training_epochs,
                gcnv_initial_temperature = gcnv_initial_temperature,
                gcnv_num_thermal_advi_iters = gcnv_num_thermal_advi_iters,
                gcnv_convergence_snr_averaging_window = gcnv_convergence_snr_averaging_window,
                gcnv_convergence_snr_trigger_threshold = gcnv_convergence_snr_trigger_threshold,
                gcnv_convergence_snr_countdown_window = gcnv_convergence_snr_countdown_window,
                gcnv_max_calling_iters = gcnv_max_calling_iters,
                gcnv_caller_update_convergence_threshold = gcnv_caller_update_convergence_threshold,
                gcnv_caller_internal_admixing_rate = gcnv_caller_internal_admixing_rate,
                gcnv_caller_external_admixing_rate = gcnv_caller_external_admixing_rate,
                gcnv_disable_annealing = gcnv_disable_annealing,
                ref_copy_number_autosomal_contigs = ref_copy_number_autosomal_contigs,
                allosomal_contigs = allosomal_contigs 
        }
    }

    output {
        Array[File] contig_ploidy_calls_tars = CNVGermlineCaseWorkflowFromCoverage.contig_ploidy_calls_tar
        # Array[Array[File]] gcnv_calls_tars = CNVGermlineCaseWorkflowFromCoverage.gcnv_call_tars
        # Array[Array[File]] gcnv_tracking_tars = CNVGermlineCaseWorkflowFromCoverage.gcnv_tracking_tars
        Array[Array[File]] genotyped_intervals_vcf = CNVGermlineCaseWorkflowFromCoverage.genotyped_intervals_vcf
        Array[Array[File]] genotyped_segments_vcf = CNVGermlineCaseWorkflowFromCoverage.genotyped_segments_vcf
        Array[Array[File]] denoised_copy_ratios = CNVGermlineCaseWorkflowFromCoverage.denoised_copy_ratios
        Array[Array[File]] read_counts = CNVGermlineCaseWorkflowFromCoverage.read_counts
    }
}

task SplitInputArray {
    Array[String] input_array
    Int num_inputs_in_scatter
    String gatk_docker

    Int machine_mem_mb = 4000
    Int disk_space_gb = 20
    Int cpu = 1
    Int preemptible_attempts = 5
    Boolean use_ssd = false

    File input_array_file = write_lines(input_array)

    command <<<
        python <<CODE
        import math
        with open("${input_array_file}", "r") as input_array_file:
            input_array = input_array_file.read().splitlines()
        num = ${num_inputs_in_scatter}
        values_to_write = [input_array[num*i:num*i+min(num, len(input_array)-num*i)] for i in range(int(math.ceil(len(input_array)/num)))]
        with open('input_array_split.tsv', 'w') as outfile:
            for i in range(len(values_to_write)):
                current_sub_array = values_to_write[i]
                for j in range(len(current_sub_array)):
                    outfile.write(current_sub_array[j] + "\t")
                outfile.write("\n")
        CODE
    >>>

    runtime {
        docker: "${gatk_docker}"
        bootDiskSizeGb: 12
        memory: machine_mem_mb + " MB"
        disks: "local-disk " + select_first([disk_space_gb, 150]) + if use_ssd then " SSD" else " HDD"
        cpu: select_first([cpu, 8])
        preemptible: select_first([preemptible_attempts, 2])
    }

    output {
        Array[Array[String]] splitArray = read_tsv("input_array_split.tsv")
    }
}

task fof_usage_task {
   File fof
   Array[File] my_files=read_lines(fof)
   command {
   #do stuff with arrays below
   #....
   }
   runtime {
       docker : "ubuntu:16.04"
       disks: "local-disk 50 HDD"
       memory: "2 GB"
   }
   output {
    Array[File] array_of_files = my_files
    }  
}