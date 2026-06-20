# Oatk profile databases

The OatkDB release `v20230921` is available under `dbs/v20230921`. Pass the
selected `.fam` file to `--oatk_mito_hmm` or `--oatk_plastid_hmm`. The pipeline
automatically validates and stages the matching `.h3f`, `.h3i`, `.h3m`, and
`.h3p` files.

## Mitochondrial profiles

| Profile | File |
| --- | --- |
| Acrogymnospermae | `acrogymnospermae_mito.fam` |
| Actinopterygii | `actinopterygii_mito.fam` |
| Amphibia | `amphibia_mito.fam` |
| Anthocerotophyta | `anthocerotophyta_mito.fam` |
| Araneae | `araneae_mito.fam` |
| Aves | `aves_mito.fam` |
| Bryophyta | `bryophyta_mito.fam` |
| Crocodylia | `crocodylia_mito.fam` |
| Dikarya | `dikarya_mito.fam` |
| Embryophyta | `embryophyta_mito.fam` |
| Hymenoptera | `hymenoptera_mito.fam` |
| Insecta | `insecta_mito.fam` |
| Lepidosauria | `lepidosauria_mito.fam` |
| Lycopodiopsida | `lycopodiopsida_mito.fam` |
| Magnoliopsida | `magnoliopsida_mito.fam` |
| Mammalia | `mammalia_mito.fam` |
| Marchantiophyta | `marchantiophyta_mito.fam` |
| Mollusca | `mollusca_mito.fam` |
| Nematoda | `nematoda_mito.fam` |
| Polypodiopsida | `polypodiopsida_mito.fam` |
| Serpentes | `serpentes_mito.fam` |
| Testudines | `testudines_mito.fam` |

## Plastid profiles

| Profile | File |
| --- | --- |
| Acrogymnospermae | `acrogymnospermae_pltd.fam` |
| Anthocerotophyta | `anthocerotophyta_pltd.fam` |
| Bryophyta | `bryophyta_pltd.fam` |
| Embryophyta | `embryophyta_pltd.fam` |
| Lycopodiopsida | `lycopodiopsida_pltd.fam` |
| Magnoliopsida | `magnoliopsida_pltd.fam` |
| Marchantiophyta | `marchantiophyta_pltd.fam` |
| Polypodiopsida | `polypodiopsida_pltd.fam` |

Use the profile matching the sample's broad taxonomic group. Animal and fungal
samples generally use only a mitochondrial profile because plastids are absent.
