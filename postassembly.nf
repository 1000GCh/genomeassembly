#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { POSTASSEMBLY } from './workflows/postassembly'

params.assembly                       = null
params.long_reads                     = null
params.hic                            = null
params.sample                         = 'sample'
params.outdir                         = 'results'
params.hic_aligner                    = 'minimap2'
params.hic_mapping_cram_chunk_size    = 10000
params.scaffolding_cool_bin_size      = 1000
params.busco_lineage                  = 'eukaryota_odb10'
params.busco_lineage_directory        = null
params.oatk_kmer_size                 = 1001
params.oatk_coverage_cutoff           = 50

workflow {
    if (!params.assembly || !params.long_reads || !params.hic) {
        error 'Required parameters: --assembly, --long_reads, and --hic'
    }

    def meta = [id: params.sample]
    def assembly = file(params.assembly, checkIfExists: true)
    def long_reads = file(params.long_reads, checkIfExists: true)
    def hic = file(params.hic, checkIfExists: true)

    if ([hic].flatten().size() != 2) {
        error "--hic must resolve to exactly two paired FASTQ files; found ${[hic].flatten().size()}"
    }

    POSTASSEMBLY(
        channel.of([meta, assembly]),
        channel.of([meta, long_reads]),
        channel.of([meta, hic]),
        params.hic_aligner,
        params.hic_mapping_cram_chunk_size,
        params.scaffolding_cool_bin_size,
        params.busco_lineage,
        params.busco_lineage_directory ? file(params.busco_lineage_directory, checkIfExists: true) : null,
        params.oatk_kmer_size,
        params.oatk_coverage_cutoff
    )

    ch_assembly_stats_files = POSTASSEMBLY.out.assembly_stats.map { _meta, files -> files }
    ch_busco_files = POSTASSEMBLY.out.busco.map { _meta, files -> files }
    ch_merqury_files = POSTASSEMBLY.out.merqury.map { _meta, files -> files }

    publish:
    scaffolding    = POSTASSEMBLY.out.scaffolding
    organelles     = POSTASSEMBLY.out.organelles
    assembly_stats = ch_assembly_stats_files
    busco          = ch_busco_files
    merqury        = ch_merqury_files
}

output {
    scaffolding {
        path { result ->
            result.output.scaffolding.values().flatten().findAll() >> 'scaffolding/'
        }
    }
    organelles {
        path { result ->
            result.output.oatk.values().flatten().findAll() >> 'organelles/'
        }
    }
    assembly_stats {
        path { files -> files >> 'qc/statistics/' }
    }
    busco {
        path { files -> files >> 'qc/busco/' }
    }
    merqury {
        path { files -> files >> 'qc/merqury/' }
    }
}
