// If a species is not provided, it bypasses the species verification
if (params.mlstSpecies{{ param_id }} == null){
   IN_expected_species_{{ pid }} = Channel.value("PASS")
} else {
    IN_expected_species_{{ pid }} = Channel.value(params.mlstSpecies{{ param_id }})
}

process mlst_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }
    // This process can only use a single CPU
    cpus 1

    input:
    set sample_id, file(assembly) from {{ input_channel }}
    val expected_species from IN_expected_species_{{ pid }}

    output:
    file '*.mlst.txt' into LOG_mlst_{{ pid }}
    set sample_id, file(assembly), file(".status") into MAIN_mlst_out_{{ pid }}
    {% with task_name="mlst" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    script:
    template "run_mlst.py"

}

process compile_mlst_{{ pid }} {

    publishDir "results/annotation/mlst_{{ pid }}/"

    input:
    file res from LOG_mlst_{{ pid }}.collect()

    output:
    file "mlst_report.tsv"

    script:
    """
    cat $res >> mlst_report.tsv
    """
}

{{ output_channel }} = Channel.create()
MAIN_mlst_out_{{ pid }}
    .filter{ it[2].text != "fail" }
    .map{ [it[0], it[1]] }
    .set{ {{output_channel}} }


{{ forks }}

