#!/usr/bin/env nextflow

nextflow.enable.dsl=2

log.info """\
         G A P - N F   P I P E L I N E
         ===================================
         outdir       : ${params.outdir}
         manifest     : ${params.manifest}
         samplesheet  : ${params.samplesheet}
         launchDir    : ${params.launchDir}
         recalculate  : ${params.recalculate}
         controls  : ${params.controls}
         """
         .stripIndent()


include { getBarcodes }           from './modules/modules'
include { gtcToFinalReport }      from './modules/modules'
include { mergeFinalReports }     from './modules/modules'
include { finalReportToOptical }  from './modules/modules'
include { OptiCall }              from './modules/modules'
include { OptiCallToGenSample }   from './modules/modules'

workflow {
  gtcToFinalReport( getBarcodes() )
  mergeFinalReports( gtcToFinalReport.out.finalreports.collect() )
  chr_list = finalReportToOptical( mergeFinalReports.out )
  OptiCall( chr_list.flatMap() )
  OptiCallToGenSample( OptiCall.out.calls,OptiCall.out.probs )
}
