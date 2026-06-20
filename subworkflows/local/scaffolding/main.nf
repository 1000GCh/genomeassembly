include { CRAM_MAP_ILLUMINA_HIC           } from '../../../subworkflows/sanger-tol/cram_map_illumina_hic'
include { BAM_STATS_SAMTOOLS              } from '../../../subworkflows/nf-core/bam_stats_samtools'
include { FASTA_BAM_SCAFFOLDING_YAHS      } from '../../../subworkflows/sanger-tol/fasta_bam_scaffolding_yahs'

include { TABIX_BGZIP as BGZIP_SCAFFOLDED } from '../../../modules/nf-core/tabix/bgzip'
include { MINIMAP2_ALIGN                  } from '../../../modules/nf-core/minimap2/align'

workflow SCAFFOLDING {
    take:
    ch_scaffolding_specs            // spec
    ch_assemblies                   // [meta, hap1, hap2]
    val_hic_aligner                 // "bwamem2" or "minimap2"
    val_hic_mapping_cram_chunk_size // int > 1
    val_cool_bin                    // int > 1

    main:
    //
    // Logic: join all the assemblies with the scaffolding specifications and
    // data, filter for those assemblies which are to be scaffolded, then
    // map out the data for the purging subworkflow.
    //
    ch_hic_mapping_inputs = ch_assemblies
        .combine(ch_scaffolding_specs)
        .filter { asm_meta, _asm1, _asm2, spec -> asm_meta.id == spec.prevID }
        .flatMap { _asm_meta, asm1, asm2, spec ->
            [asm1, asm2]
                .withIndex()
                .findAll { asm, _index -> asm }
                .collect { asm, index -> [spec + [_hap: "hap${index + 1}"], asm, spec.data.hic.reads] }
        }

    ch_hic_mapping_inputs_by_format = ch_hic_mapping_inputs.branch { meta, asm, reads ->
        def input_files = [reads].flatten()
        fastq: input_files.every { read -> read.name ==~ /.*\.(fastq|fq)(\.gz)?$/ }
        cram: input_files.every { read -> read.name.endsWith('.cram') }
        unsupported: true
    }

    ch_hic_mapping_inputs_by_format.unsupported.subscribe { _meta, _asm, reads ->
        error("Unsupported Hi-C input: ${reads}. Expected paired FASTQ or CRAM files.")
    }

    ch_cram_assemblies = ch_hic_mapping_inputs_by_format.cram.map { meta, asm, _reads -> [meta, asm] }
    ch_cram_reads = ch_hic_mapping_inputs_by_format.cram.map { meta, _asm, reads -> [meta, reads] }

    ch_fastq_assemblies = ch_hic_mapping_inputs_by_format.fastq.map { meta, asm, _reads -> [meta, asm] }
    ch_fastq_reads = ch_hic_mapping_inputs_by_format.fastq.map { meta, _asm, reads -> [meta, reads] }

    //
    // Subworkflow: Map Hi-C data to each assembly
    //
    CRAM_MAP_ILLUMINA_HIC(
        ch_cram_assemblies,
        ch_cram_reads,
        val_hic_aligner,
        val_hic_mapping_cram_chunk_size,
    )

    // Map paired Hi-C FASTQ directly to each assembly. The module sorts and indexes
    // the BAM before it is passed to the same statistics and YaHS path as CRAM data.
    MINIMAP2_ALIGN(
        ch_fastq_reads,
        ch_fastq_assemblies,
        true,
        'csi',
        false,
        false
    )

    ch_hic_assemblies = ch_cram_assemblies.mix(ch_fastq_assemblies)
    ch_hic_bam = CRAM_MAP_ILLUMINA_HIC.out.bam.mix(MINIMAP2_ALIGN.out.bam)
    ch_hic_bam_index = CRAM_MAP_ILLUMINA_HIC.out.bam_index
        .filter { _meta, idx -> idx.getExtension() == 'csi' }
        .mix(MINIMAP2_ALIGN.out.index)

    //
    // Subworkflow: Calculate stats for Hi-C mapping
    //
    ch_hic_mapping_stats_input = ch_hic_bam
        .combine(ch_hic_bam_index, by: 0)
        .combine(ch_hic_assemblies, by: 0)
        .multiMap { meta, bam, bai, asm ->
            bam: [ meta, bam, bai ]
            asm: [ meta, asm ]
        }

    BAM_STATS_SAMTOOLS(
        ch_hic_mapping_stats_input.bam,
        ch_hic_mapping_stats_input.asm
    )

    //
    // Subworkflow: scaffold assemblies using yahs and create contact maps
    //
    FASTA_BAM_SCAFFOLDING_YAHS(
        ch_hic_assemblies,
        ch_hic_bam,
        true,
        true,
        true,
        val_cool_bin
    )

    //
    // Module: bgzip all scaffolded assembly fasta
    //
    BGZIP_SCAFFOLDED(FASTA_BAM_SCAFFOLDING_YAHS.out.scaffolds_fasta)

    //
    // Logic: re-join pairs of assemblies from scaffolding to pass for genome statistics
    //
    ch_assemblies_scaffolded = FASTA_BAM_SCAFFOLDING_YAHS.out.scaffolds_fasta
        .map { meta, asm -> [meta - meta.subMap("_hap"), [meta._hap, asm]] }
        .groupTuple()
        .map { meta, hap_assemblies ->
            def asms = hap_assemblies.sort { a, b -> a[0] <=> b[0] }.collect { hap, asm -> asm }
            [meta, asms[0], asms.size() > 1 ? asms[1] : []]
        }

    //
    // Logic: combine all scaffolding outputs into a single map for ease of publishing
    //
    ch_scaffolding_output = BGZIP_SCAFFOLDED.out.output
        .join(ch_hic_bam, by: 0)
        .join(ch_hic_bam_index, by: 0)
        .join(BAM_STATS_SAMTOOLS.out.stats, by: 0)
        .join(BAM_STATS_SAMTOOLS.out.flagstat, by: 0)
        .join(BAM_STATS_SAMTOOLS.out.idxstats, by: 0)
        .join(FASTA_BAM_SCAFFOLDING_YAHS.out.scaffolds_agp, by: 0)
        .join(FASTA_BAM_SCAFFOLDING_YAHS.out.yahs_bin, by: 0)
        .join(FASTA_BAM_SCAFFOLDING_YAHS.out.yahs_inital, by: 0, remainder: true)
        .join(FASTA_BAM_SCAFFOLDING_YAHS.out.yahs_intermediate, by: 0, remainder: true)
        .join(FASTA_BAM_SCAFFOLDING_YAHS.out.yahs_log, by: 0)
        .join(FASTA_BAM_SCAFFOLDING_YAHS.out.pretext, by: 0, remainder: true)
        .join(FASTA_BAM_SCAFFOLDING_YAHS.out.pretext_png, by: 0, remainder: true)
        .join(FASTA_BAM_SCAFFOLDING_YAHS.out.cool, by: 0, remainder: true)
        .join(FASTA_BAM_SCAFFOLDING_YAHS.out.hic, by: 0, remainder: true)
        .map { spec, fasta, bam, bai, stats, flagstats, idxstats, agp, bin, initial, intermed, log, pretext, png, cool, hic ->
            return spec.subMap(["id", "stage", "data", "params", "tools"]) + [
                hap: spec._hap,
                output: [
                    scaffolding: [
                        fasta: fasta,
                        bam: bam,
                        bai: bai,
                        stats: stats,
                        flagstats: flagstats,
                        idxstats: idxstats,
                        yahs_agp: agp,
                        yahs_bin: bin,
                        yahs_initial: initial,
                        yahs_intermeriate: intermed,
                        yahs_log: log,
                        pretext: pretext,
                        pretext_png: png,
                        cool: cool,
                        hic: hic
                    ]
                ]
            ]
        }

    emit:
    scaffolded_assemblies = ch_assemblies_scaffolded
    scaffolding_output = ch_scaffolding_output
}
