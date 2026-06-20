# Post-assembly workflow

`postassembly.nf` scaffolds an existing nuclear assembly with Hi-C data, assembles
organellar genomes from the long reads with Oatk, and runs assembly statistics,
BUSCO, and MerquryFK QC on the scaffolded assembly.

## Inputs

- `--assembly`: assembled nuclear genome in FASTA format.
- `--long_reads`: long reads in FASTA/FASTQ format. Globs must be quoted.
- `--hic`: exactly two paired Hi-C FASTQ files. A quoted glob may be used.
- `--oatk_mito_hmm`: mitochondrial Oatk `.fam` profile. Optional when a plastid
  profile is supplied.
- `--oatk_plastid_hmm`: plastid Oatk `.fam` profile. Optional when a mitochondrial
  profile is supplied.

At least one Oatk profile is required. Each `.fam` file must be accompanied by
its `.h3f`, `.h3i`, `.h3m`, and `.h3p` files in the same directory. See the
[bundled database index](oatk-databases.md) for profiles under `dbs/v20230921`.

## Example

```bash
nextflow run postassembly.nf \
  --assembly assembly.fasta \
  --long_reads 'reads/*.fastq.gz' \
  --hic 'hic/*_R{1,2}.fastq.gz' \
  --oatk_mito_hmm genomeassembly/dbs/v20230921/insecta_mito.fam \
  --sample sample_name \
  --outdir results \
  -profile docker
```

The default BUSCO lineage is `eukaryota_odb10`. Override it with
`--busco_lineage`, or provide a local lineage root with
`--busco_lineage_directory`.

Hi-C FASTQ reads are aligned to the input assembly with minimap2's short-read
preset before YaHS scaffolding. The reusable scaffolding subworkflow continues
to accept CRAM input for the main specification-driven pipeline.

The workflow accepts a single (unphased) assembly. It also exposes the reusable
`POSTASSEMBLY` workflow in `workflows/postassembly.nf` for composition in another
DSL2 pipeline.
