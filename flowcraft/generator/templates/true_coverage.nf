
process true_coverage_{{ pid }} {

    // Send POST request to platform
    {% include "post.txt" ignore missing %}

    tag { sample_id }

    input:
    set sample_id, file(fastq_pair) from {{ input_channel }}

    output:
    set sample_id, file(fastq_pair) into {{ output_channel }}
    {% with task_name="true_coverage" %}
    {%- include "compiler_channels.txt" ignore missing -%}
    {% endwith %}

    """
    {
        trueCoverage_rematch.py -f $fastq_pair --species $params.species \
        -i /NGStools/true_coverage/data --json
        if ls failing* 1> /dev/null 2>&1;
        then
            parse_true_coverage.py sample_*.json failing*.json
        else
            parse_true_coverage.py sample_*.json
        fi
        echo pass > .status
    } || {
        echo fail > .status
    }
    """

}

{{ forks }}
