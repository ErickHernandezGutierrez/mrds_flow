#!/usr/bin/env nextflow

if(params.help) {
    usage = file("$baseDir/USAGE")
    cpu_count = Runtime.runtime.availableProcessors()

    bindings = ["use_fslgrad":"$params.use_fslgrad",
                "use_isotropic":"$params.use_isotropic",
                "model_selection":"$params.model_selection",
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
log.info "[MRDS options]"
log.info "Use FSL Gradient: $params.use_fslgrad"
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
    .fromFilePairs("$params.input/**/*.bval",
                   size: -1) { it.parent.name }
    .set{in_bval}

Channel
    .fromFilePairs("$params.input/**/*.bvec",
                   size: -1) { it.parent.name }
    .set{in_bvec}

Channel
    .fromFilePairs("$params.input/**/*.b",
        size: -1) { it.parent.name }
    .set{scheme_for_mrds}

Channel
    .fromFilePairs("$params.input/**/*mask.nii.gz",
        size: -1) { it.parent.name }
    .into{mask_for_mrds; mask_for_modsel; mask_for_metrics}

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

in_bval
    .combine(in_bvec, by: 0)
    .set{bval_bvec_for_scheme}

process Convert_Scheme {
    input:
    set sid, file(bval), file(bvec) from bval_bvec_for_scheme

    output:
    set sid, "${sid}__mrtrix-scheme.b" into mrtrix_scheme_for_mrds

    when:
    params.use_fslgrad

    script:
    """
    scil_convert_gradients_fsl_to_mrtrix.py ${bval} ${bvec} ${sid}__mrtrix-scheme.b
    """
}

if (params.use_fslgrad) {
    dwi_for_mrds
        .combine(mrtrix_scheme_for_mrds, by: 0)
        .combine(mask_for_mrds, by: 0)
        .set{dwi_scheme_mask_for_mrds}
} else {
    dwi_for_mrds
        .combine(scheme_for_mrds, by: 0)
        .combine(mask_for_mrds, by: 0)
        .set{dwi_scheme_mask_for_mrds}
}

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
    file "${sid}__DTInolin_COMP_SIZE.nii.gz"
    file "${sid}__DTInolin_EIGENVALUES.nii.gz"
    file "${sid}__DTInolin_ISOTROPIC.nii.gz"
    file "${sid}__DTInolin_NUM_COMP.nii.gz"
    file "${sid}__DTInolin_PDDs_CARTESIAN.nii.gz"
    file "${sid}__DTInolin_ResponseAnisotropic.txt"
    file "${sid}__DTInolin_ResponseAnisotropicMask.nii.gz"
    file "${sid}__DTInolin_ResponseIsotropic.txt"
    file "${sid}__DTInolin_ResponseIsotropicMask.nii.gz"
    file "${sid}__DTInolin_Tensor.nii.gz"
    file "${sid}__MRDS_Diff_${params.model_selection}_COMP_SIZE.nii.gz"
    file "${sid}__MRDS_Diff_${params.model_selection}_EIGENVALUES.nii.gz"
    file "${sid}__MRDS_Diff_${params.model_selection}_ISOTROPIC.nii.gz"
    file "${sid}__MRDS_Diff_${params.model_selection}_MSE.nii.gz"
    file "${sid}__MRDS_Diff_${params.model_selection}_NUM_COMP.nii.gz"
    file "${sid}__MRDS_Diff_${params.model_selection}_PDDs_CARTESIAN.nii.gz"
    file "${sid}__MRDS_Diff_V0_COMP_SIZE.nii.gz", optinal: true
    file "${sid}__MRDS_Diff_V0_EIGENVALUES.nii.gz", optinal: true
    file "${sid}__MRDS_Diff_V0_ISOTROPIC.nii.gz", optinal: true
    file "${sid}__MRDS_Diff_V0_MSE.nii.gz", optinal: true
    file "${sid}__MRDS_Diff_V0_NUM_COMP.nii.gz", optinal: true
    file "${sid}__MRDS_Diff_V0_PDDs_CARTESIAN.nii.gz", optinal: true
    file "${sid}__MRDS_Diff_V1_MSE.nii.gz"
    file "${sid}__MRDS_Diff_V2_MSE.nii.gz"
    file "${sid}__MRDS_Diff_V3_MSE.nii.gz"
    file "${sid}__MRDS_Equal_${params.model_selection}_COMP_SIZE.nii.gz"
    file "${sid}__MRDS_Equal_${params.model_selection}_EIGENVALUES.nii.gz"
    file "${sid}__MRDS_Equal_${params.model_selection}_ISOTROPIC.nii.gz"
    file "${sid}__MRDS_Equal_${params.model_selection}_MSE.nii.gz"
    file "${sid}__MRDS_Equal_${params.model_selection}_NUM_COMP.nii.gz"
    file "${sid}__MRDS_Equal_${params.model_selection}_PDDs_CARTESIAN.nii.gz"
    file "${sid}__MRDS_Equal_V0_COMP_SIZE.nii.gz", optional: true
    file "${sid}__MRDS_Equal_V0_EIGENVALUES.nii.gz", optional: true
    file "${sid}__MRDS_Equal_V0_ISOTROPIC.nii.gz", optional: true
    file "${sid}__MRDS_Equal_V0_MSE.nii.gz", optional: true
    file "${sid}__MRDS_Equal_V0_NUM_COMP.nii.gz", optional: true
    file "${sid}__MRDS_Equal_V0_PDDs_CARTESIAN.nii.gz", optional: true
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
    file "${sid}__MRDS_Fixed_V0_COMP_SIZE.nii.gz", optional: true
    file "${sid}__MRDS_Fixed_V0_EIGENVALUES.nii.gz", optional: true
    file "${sid}__MRDS_Fixed_V0_ISOTROPIC.nii.gz", optional: true
    file "${sid}__MRDS_Fixed_V0_MSE.nii.gz", optional: true
    file "${sid}__MRDS_Fixed_V0_NUM_COMP.nii.gz", optional: true
    file "${sid}__MRDS_Fixed_V0_PDDs_CARTESIAN.nii.gz", optional: true
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
    scil_fit_mrds.py ${dwi} ${scheme} --mask ${mask} --modsel ${params.model_selection.toLowerCase()} --method Diff --prefix ${sid}_ ${params.use_isotropic ? '-iso' : ''}
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
    scil_compute_todi.py ${tractogram} --out_todi_sh ${sid}__MRDS_Diff_${params.model_selection}_TOD_SH.nii.gz --reference ${dwi} --sh_basis tournier07 -f
    scil_compute_fodf_metrics.py ${sid}__MRDS_Diff_${params.model_selection}_TOD_SH.nii.gz --not_all --nufo ${sid}__MRDS_Diff_${params.model_selection}_TOD_NUFO.nii.gz --sh_basis tournier07 --rt 0.2 -f
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
    set sid, "${sid}__MRDS_Diff_TODI_EIGENVALUES.nii.gz" into eigenvalues_for_metrics
    file "${sid}__MRDS_Diff_TODI_COMP_SIZE.nii.gz"
    file "${sid}__MRDS_Diff_TODI_ISOTROPIC.nii.gz"
    file "${sid}__MRDS_Diff_TODI_NUM_COMP.nii.gz"
    file "${sid}__MRDS_Diff_TODI_PDDs_CARTESIAN.nii.gz"

    script:
    """
    scil_mrds_modsel_todi.py ${nufo} ${dwi} TODI --N1 ${n1_compsize} ${n1_eigen} ${n1_iso} ${n1_numcomp} ${n1_pdds} --N2 ${n2_compsize} ${n2_eigen} ${n2_iso} ${n2_numcomp} ${n2_pdds} --N3 ${n3_compsize} ${n3_eigen} ${n3_iso} ${n3_numcomp} ${n3_pdds} --prefix ${sid}_ --mask ${mask}
    """
}

eigenvalues_for_metrics
    .combine(mask_for_metrics, by: 0)
    .set{eigenvalues_mask_for_metrics}

process MRDS_Metrics {
    input:
    set sid, eigenvalues, mask from eigenvalues_mask_for_metrics

    output:
    file "${sid}__MRDS_Diff_TODI_AD.nii.gz"
    file "${sid}__MRDS_Diff_TODI_RD.nii.gz"
    file "${sid}__MRDS_Diff_TODI_MD.nii.gz"
    file "${sid}__MRDS_Diff_TODI_FA.nii.gz"

    script:
    """
    scil_compute_mrds_metrics.py ${eigenvalues} ${sid}__MRDS_Diff_TODI_AD.nii.gz ${sid}__MRDS_Diff_TODI_RD.nii.gz ${sid}__MRDS_Diff_TODI_MD.nii.gz ${sid}__MRDS_Diff_TODI_FA.nii.gz
    """
}
