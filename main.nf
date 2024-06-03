#!/usr/bin/env nextflow

if(params.help) {
    usage = file("$baseDir/USAGE")
    cpu_count = Runtime.runtime.availableProcessors()

    bindings = ["use_isotropic":"$params.use_isotropic",
                "model_selection":"$params.model_selection",
                "register_processes":"$params.register_processes",
                "rbx_processes":"$params.rbx_processes",
                "single_dataset_size_GB":"$params.single_dataset_size_GB",
                "cpu_count":"$cpu_count"]

    engine = new groovy.text.SimpleTemplateEngine()
    template = engine.createTemplate(usage.text).make(bindings)
    print template.toString()
    return
}

log.info "SCIL MRDS pipeline"
log.info "=========================="
log.info ""
log.info "Start time: $workflow.start"
log.info ""

log.debug "[Command-line]"
log.debug "$workflow.commandLine"
log.debug ""

log.info "[Git Info]"
log.info "$workflow.repository - $workflow.revision [$workflow.commitId]"
log.info ""

log.info "Options"
log.info "======="
log.info ""
log.info "[Recobundles options]"
log.info "Use Isotropic Compartment: $params.use_isotropic"
log.info "Model Selector: $params.model_selection"
log.info ""
log.info ""

log.info "Input: $params.input"
root = file(params.input)
/* Watch out, files are ordered alphabetically in channel */
tractogram_for_todi = Channel
     .fromFilePairs("$root/**/{*tracking*.*,}",
                    size: -1,
                    maxDepth:1) {it.parent.name}

Channel
    .fromFilePairs("$params.input/**/*dwi.nii.gz",
        size: -1) { it.parent.name }
    .into{dwi_for_mrds; dwi_for_todi; dwi_for_modsel}

Channel
    .fromFilePairs("$params.input/**/*scheme",
        size: -1) { it.parent.name }
    .set{scheme_for_mrds}

Channel
    .fromFilePairs("$params.input/**/*mask.nii.gz",
        size: -1) { it.parent.name }
    .into{mask_for_mrds; mask_for_modsel}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

dwi_for_mrds
    .combine(scheme_for_mrds, by: 0)
    .combine(mask_for_mrds, by: 0)
    .set{dwi_scheme_mask_for_mrds}

process Fit_MRDS {
    input:
    set sid, file(dwi), file(scheme), file(mask) from dwi_scheme_mask_for_mrds

    output:
    set sid, "${sid}__MRDS_Diff_V1_COMP_SIZE.nii.gz",\
             "${sid}__MRDS_Diff_V1_EIGENVALUES.nii.gz",\
             "${sid}__MRDS_Diff_V1_ISOTROPIC.nii.gz",\
             "${sid}__MRDS_Diff_V1_NUM_COMP.nii.gz",\
             "${sid}__MRDS_Diff_V1_PDDs_CARTESIAN.nii.gz",\
             "${sid}__MRDS_Diff_V2_COMP_SIZE.nii.gz",\
             "${sid}__MRDS_Diff_V2_EIGENVALUES.nii.gz",\
             "${sid}__MRDS_Diff_V2_ISOTROPIC.nii.gz",\
             "${sid}__MRDS_Diff_V2_NUM_COMP.nii.gz",\
             "${sid}__MRDS_Diff_V2_PDDs_CARTESIAN.nii.gz",\
             "${sid}__MRDS_Diff_V3_COMP_SIZE.nii.gz",\
             "${sid}__MRDS_Diff_V3_EIGENVALUES.nii.gz",\
             "${sid}__MRDS_Diff_V3_ISOTROPIC.nii.gz",\
             "${sid}__MRDS_Diff_V3_NUM_COMP.nii.gz",\
             "${sid}__MRDS_Diff_V3_PDDs_CARTESIAN.nii.gz" into mrds_for_modsel

    file "${sid}__MRDS_Diff_${params.model_selection}_COMP_SIZE.nii.gz"
    file "${sid}__MRDS_Diff_${params.model_selection}_EIGENVALUES.nii.gz"
    file "${sid}__MRDS_Diff_${params.model_selection}_ISOTROPIC.nii.gz"
    file "${sid}__MRDS_Diff_${params.model_selection}_MSE.nii.gz"
    file "${sid}__MRDS_Diff_${params.model_selection}_NUM_COMP.nii.gz"
    file "${sid}__MRDS_Diff_${params.model_selection}_PDDs_CARTESIAN.nii.gz"
    file "${sid}__MRDS_Diff_V0_COMP_SIZE.nii.gz"
    file "${sid}__MRDS_Diff_V0_EIGENVALUES.nii.gz"
    file "${sid}__MRDS_Diff_V0_ISOTROPIC.nii.gz"
    file "${sid}__MRDS_Diff_V0_MSE.nii.gz"
    file "${sid}__MRDS_Diff_V0_NUM_COMP.nii.gz"
    file "${sid}__MRDS_Diff_V0_PDDs_CARTESIAN.nii.gz"
    file "${sid}__MRDS_Diff_V1_MSE.nii.gz"
    file "${sid}__MRDS_Diff_V2_MSE.nii.gz"
    file "${sid}__MRDS_Diff_V3_MSE.nii.gz"
    file "${sid}__MRDS_Equal_${params.model_selection}_COMP_SIZE.nii.gz"
    file "${sid}__MRDS_Equal_${params.model_selection}_EIGENVALUES.nii.gz"
    file "${sid}__MRDS_Equal_${params.model_selection}_ISOTROPIC.nii.gz"
    file "${sid}__MRDS_Equal_${params.model_selection}_MSE.nii.gz"
    file "${sid}__MRDS_Equal_${params.model_selection}_NUM_COMP.nii.gz"
    file "${sid}__MRDS_Equal_${params.model_selection}_PDDs_CARTESIAN.nii.gz"
    file "${sid}__MRDS_Equal_V0_COMP_SIZE.nii.gz"
    file "${sid}__MRDS_Equal_V0_EIGENVALUES.nii.gz"
    file "${sid}__MRDS_Equal_V0_ISOTROPIC.nii.gz"
    file "${sid}__MRDS_Equal_V0_MSE.nii.gz"
    file "${sid}__MRDS_Equal_V0_NUM_COMP.nii.gz"
    file "${sid}__MRDS_Equal_V0_PDDs_CARTESIAN.nii.gz"
    file "${sid}__MRDS_Equal_V1_COMP_SIZE.nii.gz"
    file "${sid}__MRDS_Equal_V1_EIGENVALUES.nii.gz"
    file "${sid}__MRDS_Equal_V1_ISOTROPIC.nii.gz"
    file "${sid}__MRDS_Equal_V1_MSE.nii.gz"
    file "${sid}__MRDS_Equal_V1_NUM_COMP.nii.gz"
    file "${sid}__MRDS_Equal_V1_PDDs_CARTESIAN.nii.gz"
    file "${sid}__MRDS_Equal_V2_COMP_SIZE.nii.gz"
    file "${sid}__MRDS_Equal_V2_EIGENVALUES.nii.gz"
    file "${sid}__MRDS_Equal_V2_ISOTROPIC.nii.gz"
    file "${sid}__MRDS_Equal_V2_MSE.nii.gz"
    file "${sid}__MRDS_Equal_V2_NUM_COMP.nii.gz"
    file "${sid}__MRDS_Equal_V2_PDDs_CARTESIAN.nii.gz"
    file "${sid}__MRDS_Equal_V3_COMP_SIZE.nii.gz"
    file "${sid}__MRDS_Equal_V3_EIGENVALUES.nii.gz"
    file "${sid}__MRDS_Equal_V3_ISOTROPIC.nii.gz"
    file "${sid}__MRDS_Equal_V3_MSE.nii.gz"
    file "${sid}__MRDS_Equal_V3_NUM_COMP.nii.gz"
    file "${sid}__MRDS_Equal_V3_PDDs_CARTESIAN.nii.gz"
    file "${sid}__MRDS_Fixed_${params.model_selection}_COMP_SIZE.nii.gz"
    file "${sid}__MRDS_Fixed_${params.model_selection}_EIGENVALUES.nii.gz"
    file "${sid}__MRDS_Fixed_${params.model_selection}_ISOTROPIC.nii.gz"
    file "${sid}__MRDS_Fixed_${params.model_selection}_MSE.nii.gz"
    file "${sid}__MRDS_Fixed_${params.model_selection}_NUM_COMP.nii.gz"
    file "${sid}__MRDS_Fixed_${params.model_selection}_PDDs_CARTESIAN.nii.gz"
    file "${sid}__MRDS_Fixed_V0_COMP_SIZE.nii.gz"
    file "${sid}__MRDS_Fixed_V0_EIGENVALUES.nii.gz"
    file "${sid}__MRDS_Fixed_V0_ISOTROPIC.nii.gz"
    file "${sid}__MRDS_Fixed_V0_MSE.nii.gz"
    file "${sid}__MRDS_Fixed_V0_NUM_COMP.nii.gz"
    file "${sid}__MRDS_Fixed_V0_PDDs_CARTESIAN.nii.gz"
    file "${sid}__MRDS_Fixed_V1_COMP_SIZE.nii.gz"
    file "${sid}__MRDS_Fixed_V1_EIGENVALUES.nii.gz"
    file "${sid}__MRDS_Fixed_V1_ISOTROPIC.nii.gz"
    file "${sid}__MRDS_Fixed_V1_MSE.nii.gz"
    file "${sid}__MRDS_Fixed_V1_NUM_COMP.nii.gz"
    file "${sid}__MRDS_Fixed_V1_PDDs_CARTESIAN.nii.gz"
    file "${sid}__MRDS_Fixed_V2_COMP_SIZE.nii.gz"
    file "${sid}__MRDS_Fixed_V2_EIGENVALUES.nii.gz"
    file "${sid}__MRDS_Fixed_V2_ISOTROPIC.nii.gz"
    file "${sid}__MRDS_Fixed_V2_MSE.nii.gz"
    file "${sid}__MRDS_Fixed_V2_NUM_COMP.nii.gz"
    file "${sid}__MRDS_Fixed_V2_PDDs_CARTESIAN.nii.gz"
    file "${sid}__MRDS_Fixed_V3_COMP_SIZE.nii.gz"
    file "${sid}__MRDS_Fixed_V3_EIGENVALUES.nii.gz"
    file "${sid}__MRDS_Fixed_V3_ISOTROPIC.nii.gz"
    file "${sid}__MRDS_Fixed_V3_MSE.nii.gz"
    file "${sid}__MRDS_Fixed_V3_NUM_COMP.nii.gz"
    file "${sid}__MRDS_Fixed_V3_PDDs_CARTESIAN.nii.gz"

    script:
    """
    mdtmrds ${dwi} ${scheme} ${sid}_ -correction 0 -response 0,0,0.003 -mask ${mask} -modsel ${params.model_selection.toLowerCase()} -each -intermediate -iso -mse -method diff
    """
}

dwi_for_todi
    .combine(tractogram_for_todi, by: 0)
    .set{dwi_tractogram_for_todi}

process Compute_TODI {
    input:
    set sid, file(dwi), file(tractogram) from dwi_tractogram_for_todi

    output:
    set sid, "${sid}__MRDS_Diff_${params.model_selection}_TOD_NUFO.nii.gz" into nufo_for_modsel
    file "${sid}__MRDS_Diff_${params.model_selection}_TOD_SH.nii.gz"

    script:
    """
    scil_mrds_compute_todi.py ${tractogram} ${dwi} ${sid}__MRDS_Diff_${params.model_selection}_TOD_SH.nii.gz ${sid}__MRDS_Diff_${params.model_selection}_TOD_NUFO.nii.gz
    """
}

nufo_for_modsel
    .combine(dwi_for_modsel, by: 0)
    .combine(mask_for_modsel, by: 0)
    .combine(mrds_for_modsel, by: 0)
    .set{dwi_nufo_mrds_for_modsel}

process Modsel_TODI {
    input:
    set sid, file(nufo), file(dwi), file(mask), file(n1_compsize), file(n1_eigen), file(n1_iso), file(n1_numcomp), file(n1_pdds),\
                                                file(n2_compsize), file(n2_eigen), file(n2_iso), file(n2_numcomp), file(n2_pdds),\
                                                file(n3_compsize), file(n3_eigen), file(n3_iso), file(n3_numcomp), file(n3_pdds) from dwi_nufo_mrds_for_modsel

    output:
    file "${sid}__MRDS_Diff_TODI_COMP_SIZE.nii.gz"
    file "${sid}__MRDS_Diff_TODI_EIGENVALUES.nii.gz"
    file "${sid}__MRDS_Diff_TODI_ISOTROPIC.nii.gz"
    file "${sid}__MRDS_Diff_TODI_NUM_COMP.nii.gz"
    file "${sid}__MRDS_Diff_TODI_PDDs_CARTESIAN.nii.gz"

    script:
    """
    scil_mrds_modsel_todi.py ${nufo} ${dwi} TODI --N1 ${n1_compsize} ${n1_eigen} ${n1_iso} ${n1_numcomp} ${n1_pdds} --N2 ${n2_compsize} ${n2_eigen} ${n2_iso} ${n2_numcomp} ${n2_pdds} --N3 ${n3_compsize} ${n3_eigen} ${n3_iso} ${n3_numcomp} ${n3_pdds} --prefix ${sid}_ --mask ${mask}
    """
}