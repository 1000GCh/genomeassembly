/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Scaffold a pre-assembled genome, assemble organelles, and assess genome quality
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SCAFFOLDING       } from '../subworkflows/local/scaffolding'
include { OATK_ASSEMBLY     } from '../subworkflows/local/oatk_assembly'
include { GENOME_STATISTICS } from '../subworkflows/sanger-tol/genome_statistics'
include { FASTK_FASTK       } from '../modules/nf-core/fastk/fastk'

workflow POSTASSEMBLY {
    take:
    ch_assembly                    // [meta, assembly FASTA]
    ch_long_reads                  // [meta, long-read FASTQ files]
    ch_hic                         // [meta, paired Hi-C FASTQ files]
    val_hic_aligner
    val_hic_mapping_cram_chunk_size
    val_scaffolding_cool_bin_size
    val_busco_lineage
    val_busco_lineage_directory
    val_oatk_kmer_size
    val_oatk_coverage_cutoff
    val_oatk_mito_hmm
    val_oatk_plastid_hmm

    main:
    ch_inputs = ch_assembly
        .join(ch_long_reads, by: 0)
        .join(ch_hic, by: 0)

    ch_assemblies = ch_inputs.map { meta, assembly, _reads, _hic ->
        def input_meta = meta + [id: "${meta.id}.input"]
        [input_meta, assembly, []]
    }

    ch_scaffolding_specs = ch_inputs.map { meta, _assembly, _reads, hic ->
        [
            id: "${meta.id}.scaffolded",
            stage: 'scaffolding',
            prevID: "${meta.id}.input",
            data: [hic: [reads: hic]],
            params: [:],
            tools: [:]
        ]
    }

    SCAFFOLDING(
        ch_scaffolding_specs,
        ch_assemblies,
        val_hic_aligner,
        val_hic_mapping_cram_chunk_size,
        val_scaffolding_cool_bin_size
    )

    ch_oatk_specs = ch_inputs.map { meta, _assembly, reads, _hic ->
        [
            id: "${meta.id}.organelles",
            stage: 'oatk',
            data: [long_read: [reads: reads]],
            params: [
                oatk_kmer_size: val_oatk_kmer_size,
                oatk_coverage_cutoff: val_oatk_coverage_cutoff,
                oatk_arguments: '',
                oatk_mito_hmm: val_oatk_mito_hmm,
                oatk_plastid_hmm: val_oatk_plastid_hmm
            ],
            tools: [:]
        ]
    }

    OATK_ASSEMBLY(ch_oatk_specs)

    ch_qc_reads = ch_long_reads
        .combine(SCAFFOLDING.out.scaffolded_assemblies)
        .filter { input_meta, _reads, scaffold_meta, _hap1, _hap2 ->
            scaffold_meta.id == "${input_meta.id}.scaffolded"
        }
        .map { _input_meta, reads, scaffold_meta, _hap1, _hap2 ->
            [scaffold_meta, reads]
        }

    FASTK_FASTK(ch_qc_reads)

    ch_fastk = FASTK_FASTK.out.hist
        .combine(FASTK_FASTK.out.ktab, by: 0)
        .map { meta, hist, ktab -> [meta, hist, ktab, [], []] }

    ch_busco_lineage = SCAFFOLDING.out.scaffolded_assemblies
        .map { meta, _hap1, _hap2 -> [meta, val_busco_lineage] }

    GENOME_STATISTICS(
        SCAFFOLDING.out.scaffolded_assemblies,
        ch_fastk,
        ch_busco_lineage,
        val_busco_lineage_directory
    )

    emit:
    scaffolds       = SCAFFOLDING.out.scaffolded_assemblies
    scaffolding     = SCAFFOLDING.out.scaffolding_output
    organelles      = OATK_ASSEMBLY.out.oatk_output
    assembly_stats  = GENOME_STATISTICS.out.stats
    busco           = GENOME_STATISTICS.out.busco
    merqury         = GENOME_STATISTICS.out.merqury
}
