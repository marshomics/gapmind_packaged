#!/usr/bin/env nextflow
/*
 * GapMind (amino-acid biosynthesis + carbon catabolism) at scale.
 * Nextflow DSL2 + container (Singularity/Apptainer or Docker).
 *
 *   PREPARE_DB ---> GAPMIND_BATCH (per batch of genomes) ---> MERGE ---> PRESENCE
 *                                                                    \--> PLOTS
 *
 * Targets the prebuilt-database analysis path; the container carries the tools
 * and the patched PaperBLAST code, so there is no host-environment setup.
 */
nextflow.enable.dsl = 2

def helpMessage() {
    log.info """
    GapMind (Nextflow + Singularity)

    Usage:
      nextflow run . --input samplesheet.csv --outdir results \\
               -profile <local|sge|slurm>,<singularity|docker>

    Required:
      --input           CSV with a header and columns: sample,faa
                        (sample optional; derived from the filename if absent)

    Key options (see nextflow.config):
      --sets            'aa carbon'         which pathway sets to run
      --batch_size      500                 genomes per task
      --db_dir          null                use an existing db dir instead of downloading
      --pa_mode         probably|strict     presence/absence threshold
      --keep_cand       true                keep the large per-candidate table
      --knownsim        true                amino-acid known-gap comparison (needs usearch)
      --singularity_image  gapmind.sif      image for -profile singularity
    """.stripIndent()
}

// ---------------------------------------------------------------------------
process PREPARE_DB {
    tag 'prepare_db'
    label 'process_low'

    input:
    val sets

    output:
    path 'db', emit: db

    script:
    """
    prepare_db.sh "${sets}" db "${params.prebuilt_base}" ${params.search_tool}
    """

    stub:
    """
    for s in ${sets}; do
      mkdir -p db/path.\$s
      : > db/path.\$s/steps.db ; : > db/path.\$s/curated.db
      : > db/path.\$s/curated.faa ; : > db/path.\$s/curated.faa.dmnd
    done
    """
}

process GAPMIND_BATCH {
    tag "batch_${id}"
    label 'process_medium'

    input:
    tuple val(id), val(names), path(faas)
    path db

    output:
    path "batch_${id}", emit: results

    script:
    def org = [names, faas].transpose().collect { n, f -> "file:${f.name}:${n}" }.join('\n')
    """
    cat > batch.orgfile <<'ORGFILE'
${org}
ORGFILE
    run_batch.sh ${id} batch.orgfile ${db} "${params.sets}" ${task.cpus} ${params.search_tool} ${params.knownsim ? 1 : 0}
    """

    stub:
    """
    mkdir -p batch_${id}
    printf 'orgId\\tgdb\\tgid\\tgenomeName\\tnProteins\\n' > batch_${id}/orgs.org
    i=0
    for n in ${names.join(' ')}; do i=\$((i+1)); printf 'local__%05d\\tlocal\\t%05d\\t%s\\t3000\\n' \$i \$i "\$n" >> batch_${id}/orgs.org; done
    for s in ${params.sets}; do
      printf 'orgId\\tgdb\\tgid\\tpathway\\trule\\tnHi\\tnMed\\tnLo\\tscore\\texpandedPath\\tpath\\tpath2\\n' > batch_${id}/\$s.sum.rules
      printf 'orgId\\tgdb\\tgid\\tpathway\\tstep\\tonBestPath\\tscore\\tlocusId\\tsysName\\tscore2\\tlocusId2\\tsysName2\\n' > batch_${id}/\$s.sum.steps
    done
    """
}

process MERGE {
    label 'process_low'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path batchdirs
    val sets

    output:
    path 'merged', emit: merged

    script:
    """
    merge_tables.sh merged "${sets}" ${params.keep_cand ? 1 : 0} batch_*
    """

    stub:
    """
    mkdir -p merged
    printf 'orgId\\tgdb\\tgid\\tgenomeName\\tnProteins\\n' > merged/orgs.tsv
    for s in ${sets}; do : > merged/\$s.sum.rules ; : > merged/\$s.sum.steps ; done
    """
}

process PRESENCE {
    label 'process_low'
    publishDir "${params.outdir}/presence", mode: 'copy'

    input:
    path merged

    output:
    path '*.presence.tsv'
    path '*.confidence.tsv'
    path '*.pathways.tsv', optional: true

    script:
    """
    presence_absence.py --tables ${merged} --orgs ${merged}/orgs.tsv \\
        --sets "${params.sets}" --code-dir "\${GAPMIND_DIR:-}" --mode ${params.pa_mode} --out .
    """

    stub:
    """
    for s in ${params.sets}; do : > \${s}.presence.tsv ; : > \${s}.confidence.tsv ; done
    """
}

process PLOTS {
    label 'process_low'
    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path merged

    output:
    path '*.png'
    path '*.svg'
    path 'summary_stats.tsv'

    script:
    """
    make_plots.py --tables ${merged} --orgs ${merged}/orgs.tsv \\
        --sets "${params.sets}" --code-dir "\${GAPMIND_DIR:-}" --out .
    """

    stub:
    """
    : > summary_stats.tsv
    for s in ${params.sets}; do : > \${s}_pathway_prevalence.png ; : > \${s}_pathway_prevalence.svg ; done
    """
}

// ---------------------------------------------------------------------------
workflow {
    if (params.help) { helpMessage(); return }
    if (!params.input) { error "--input samplesheet.csv is required (columns: sample,faa)" }

    // Databases: use an existing dir, or download+format the prebuilt ones.
    if (params.db_dir) {
        ch_db = Channel.value(file(params.db_dir, checkIfExists: true))
    } else {
        PREPARE_DB(params.sets)
        ch_db = PREPARE_DB.out.db
    }

    // Samplesheet -> (name, faa) -> fixed-size batches, each tagged with an id.
    ch_batches = Channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def faa  = file(row.faa, checkIfExists: true)
            def name = (row.sample && row.sample.toString().trim()) ? row.sample.toString().trim() : faa.simpleName
            tuple(name.replaceAll(/[^A-Za-z0-9_.-]/, '_'), faa)
        }
        .collate(params.batch_size)
        .toList()
        .flatMap { batches ->
            batches.withIndex().collect { b, i ->
                tuple(String.format('%05d', i + 1), b.collect { it[0] }, b.collect { it[1] })
            }
        }

    GAPMIND_BATCH(ch_batches, ch_db)
    MERGE(GAPMIND_BATCH.out.results.collect(), params.sets)
    PRESENCE(MERGE.out.merged)
    PLOTS(MERGE.out.merged)
}
