#!/usr/bin/env python3

import os
import sys
import shutil
import logging
import argparse
import logging.config

from distutils.dir_util import copy_tree

from os.path import join, dirname

try:
    from __init__ import __version__, __build__
    from generator.engine import NextflowGenerator
    from generator.inspect import NextflowInspector
    from generator.report import FlowcraftReport
    from generator.process_collector import collect_process_map
    from generator.recipe import brew_innuendo, brew_recipe, list_recipes
    from generator.pipeline_parser import parse_pipeline, SanityError
    from generator.process_details import proc_collector, colored_print
    import generator.error_handling as eh
except ImportError as e:
    from flowcraft import __version__, __build__
    from flowcraft.generator.engine import NextflowGenerator
    from flowcraft.generator.inspect import NextflowInspector
    from flowcraft.generator.report import FlowcraftReport
    from flowcraft.generator.recipe import brew_innuendo, \
        brew_recipe, list_recipes
    from flowcraft.generator.pipeline_parser import parse_pipeline, \
        SanityError
    from flowcraft.generator.process_details import proc_collector, \
        colored_print
    import flowcraft.generator.error_handling as eh
    from flowcraft.generator.process_collector import collect_process_map

logger = logging.getLogger("main")


def get_args(args=None):

    parser = argparse.ArgumentParser(
        description="A Nextflow pipeline generator")

    subparsers = parser.add_subparsers(help="Select which mode to run",
                                       dest="main_op")

    # BUILD MODE
    build_parser = subparsers.add_parser("build",
                                         help="Build a nextflow pipeline")

    group_lists = build_parser.add_mutually_exclusive_group()

    build_parser.add_argument(
        "-t", "--tasks", type=str, dest="tasks",
        help="Space separated tasks of the pipeline")
    build_parser.add_argument(
        "-r", "--recipe", dest="recipe",
        help="Use one of the available recipes")
    build_parser.add_argument(
        "-o", dest="output_nf", help="Name of the pipeline file")
    build_parser.add_argument(
        "-n", dest="pipeline_name", default="flowcraft",
        help="Provide a name for your pipeline.")
    build_parser.add_argument(
        "--merge-params", dest="merge_params", action="store_true",
        help="Merges identical parameters from multiple components into the "
             "same one. Otherwise, the parameters will be separated and unique"
             " to each component.")
    build_parser.add_argument(
        "--pipeline-only", dest="pipeline_only", action="store_true",
        help="Write only the pipeline files and not the templates, bin, and"
             " lib folders.")
    build_parser.add_argument(
        "-nd", "--no-dependecy", dest="no_dep", action="store_false",
        help="Do not automatically add dependencies to the pipeline.")
    build_parser.add_argument(
        "-c", "--check-pipeline", dest="check_only", action="store_const",
        const=True, help="Check only the validity of the pipeline "
                         "string and exit.")
    group_lists.add_argument(
        "-L", "--component-list", action="store_const", dest="detailed_list",
        const=True, help="Print a detailed description for all the "
                         "currently available processes.")
    group_lists.add_argument(
        "-l", "--component-list-short", action="store_const", dest="short_list",
        const=True, help="Print a short list of the currently "
                         "available processes.")
    group_lists.add_argument(
        "--recipe-list", dest="recipe_list", action="store_const", const=True,
        help="Print a short list of the currently available recipes."
    )
    group_lists.add_argument(
        "--recipe-list-short", dest="recipe_list_short", action="store_const",
        const=True, help="Print a condensed list of the currently available "
                         "recipes"
    )
    build_parser.add_argument(
        "-cr", "--check-recipe", dest="check_recipe",
        action="store_const", const=True,
        help="Check tasks that the recipe contain and "
             "their flow. This option might be useful "
             "if a user wants to change some components "
             "of a given recipe, by using the -t option.")
    build_parser.add_argument(
        "--export-params", dest="export_params", action="store_const",
        const=True, help="Only export the parameters for the provided "
                         "components (via -t option) in JSON format to stdout. "
                         "No pipeline will be generated with this option."
    )
    build_parser.add_argument(
        "--export-directives", dest="export_directives", action="store_const",
        const=True, help="Only export the directives for the provided "
                         "components (via -t option) in JSON format to stdout. "
                         "No pipeline will be generated with this option."
    )
    build_parser.add_argument(
        "-ft", "--fetch-tags", dest="fetch_docker_tags",
        action="store_const", const=True, help="Allows to fetch all docker tags"
                                               " for the components listed with"
                                               " the -t flag."
    )

    # GENERAL OPTIONS
    parser.add_argument(
        "--debug", dest="debug", action="store_const", const=True,
        help="Set log to debug mode")
    parser.add_argument(
        "-v", "--version", dest="version", action="store_const", const=True,
        help="Show version and exit.")

    # INSPECT MODE
    inspect_parser = subparsers.add_parser("inspect",
                                           help="Inspect the progress of a "
                                                "pipeline execution")
    inspect_parser.add_argument(
        "-i", dest="trace_file", default="pipeline_stats.txt",
        help="Specify the nextflow trace file."
    )
    inspect_parser.add_argument(
        "-r", dest="refresh_rate", default=0.02,
        help="Set the refresh frequency for the continuous inspect functions"
    )
    inspect_parser.add_argument(
        "-m", "--mode", dest="mode", default="overview",
        choices=["overview", "broadcast"],
        help="Specify the inspection run mode."
    )
    inspect_parser.add_argument(
        "-u", "--url", dest="url", default="http://www.flowcraft.live:80/",
        help="Specify the URL to where the data should be broadcast"
    )
    inspect_parser.add_argument(
        "--pretty", dest="pretty", action="store_const", const=True,
        help="Pretty inspection mode that removes usual reporting processes."
    )

    # REPORT MODE
    reports_parser = subparsers.add_parser("report",
                                           help="Broadcast the report of "
                                                "a pipeline")
    reports_parser.add_argument(
        "-i", dest="report_file",
        default="pipeline_report/pipeline_report.json",
        help="Specify the path to the pipeline report JSON file."
    )
    reports_parser.add_argument(
        "-u", "--url", dest="url", default="http://www.flowcraft.live:80/",
        help="Specify the URL to where the data should be broadcast"
    )
    reports_parser.add_argument(
        "--trace-file", dest="trace_file", default="pipeline_stats.txt",
        help="Specify the nextflow trace file. Only applicable in combination "
             "with --watch option."
    )
    reports_parser.add_argument(
        "--log-file", dest="log_file", default=".nextflow.log",
        help="Specify the nextflow log file. Only applicable in combination "
             "with --watch option."
    )
    reports_parser.add_argument(
        "-w", "--watch", dest="watch",  action="store_const", const=True,
        help="Run the report in watch mode. This option will track the "
             "generation of reports during the execution of the pipeline, "
             "allowing for the visualization of the reports in real-time"
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    return parser.parse_args(args)


def validate_build_arguments(args):

    # Skip all checks when listing the processes
    if args.detailed_list or args.short_list:
        return

    # Skip all checks when exporting parameters AND providing at least one
    # component
    if args.export_params or args.export_directives or args.fetch_docker_tags:
        # Check if components provided
        if not args.tasks:
            logger.error(colored_print(
                "At least one component needs to be provided via the -t option"
                " when exporting parameters in JSON format"
            ))
            sys.exit(1)
        return

    # When none of the main run options is specified
    if not args.tasks and not args.recipe and not args.check_only \
            and not args.detailed_list and not args.short_list:
        logger.error(colored_print(
            "At least one of these options is required: -t, -r, -c, "
            "-l, -L", "red_bold"))
        sys.exit(1)

    # When the build mode is active via tasks or recipe, but no output file
    # option has been provided
    if (args.tasks or args.recipe) and not args.check_recipe \
            and not args.output_nf:
        logger.error(colored_print(
            "Please provide the path and name of the pipeline file using the"
            " -o option.", "red_bold"))
        sys.exit(1)

    if args.output_nf:
        if not os.path.basename(args.output_nf):
            logger.error(colored_print(
                "Output pipeline path '{}' missing a name (only the directory "
                "path was provided)".format(args.output_nf), "red_bold"))
            sys.exit(1)

        parsed_output_nf = (args.output_nf if args.output_nf.endswith(".nf")
                            else "{}.nf".format(args.output_nf.strip()))
        opath = parsed_output_nf
        if os.path.dirname(opath):
            parent_dir = os.path.dirname(opath)
            if not os.path.exists(parent_dir):
                logger.error(colored_print(
                    "The provided directory '{}' does not exist.".format(
                        parent_dir), "red_bold"))
                sys.exit(1)

        return parsed_output_nf


def copy_project(path):
    """

    Parameters
    ----------
    path

    Returns
    -------

    """

    # Get nextflow repo directory
    repo_dir = dirname(os.path.abspath(__file__))

    # Get target directory
    target_dir = dirname(path)

    # Copy templates
    copy_tree(join(repo_dir, "templates"), join(target_dir, "templates"))

    # Copy Helper scripts
    copy_tree(join(repo_dir, "lib"), join(target_dir, "lib"))

    # Copy resources dir
    copy_tree(join(repo_dir, "resources"), join(target_dir, "resources"))

    # Copy bin scripts
    copy_tree(join(repo_dir, "bin"), join(target_dir, "bin"))

    # Copy default config file
    shutil.copy(join(repo_dir, "nextflow.config"),
                join(target_dir, "nextflow.config"))

    # Copy static profiles file
    shutil.copy(join(repo_dir, "profiles.config"),
                join(target_dir, "profiles.config"))


def build(args):

    # Disable standard logging for stdout when the following modes are
    #  executed:
    if args.export_params or args.export_directives or args.fetch_docker_tags:
        logger.setLevel(logging.ERROR)

    if args.recipe_list_short:
        list_recipes()

    if args.recipe_list:
        list_recipes(full=True)

    welcome = [
        "========= F L O W C R A F T =========",
        "Build mode\n"
        "version: {}".format(__version__),
        "build: {}".format(__build__),
        "====================================="
    ]

    parsed_output_nf = validate_build_arguments(args)

    logger.info(colored_print("\n".join(welcome), "green_bold"))

    # If a recipe is specified, build pipeline based on the
    # appropriate recipe
    if args.recipe:
        if args.recipe == "innuendo":
            pipeline_string = brew_innuendo(args)
        else:
            # pipeline_string = available_recipes[args.recipe]
            pipeline_string = brew_recipe(args.recipe)
            if args.tasks:
                logger.warning(colored_print(
                    "-t parameter will be ignored for recipe: {}\n".format(
                        args.recipe), "yellow_bold")
                )

        if args.check_recipe:
            logger.info(colored_print("Pipeline string for recipe: {}"
                                      .format(args.recipe), "purple_bold"))
            logger.info(pipeline_string)
            sys.exit(0)
    else:
        pipeline_string = args.tasks

    process_map = collect_process_map()

    # used for lists print
    proc_collector(process_map, args, pipeline_string)

    try:
        logger.info(colored_print("Checking pipeline for errors..."))
        pipeline_list = parse_pipeline(pipeline_string)
    except SanityError as e:
        logger.error(colored_print(e.value, "red_bold"))
        sys.exit(1)
    logger.debug("Pipeline successfully parsed: {}".format(pipeline_list))

    # Exit if only the pipeline parser needs to be checked
    if args.check_only:
        sys.exit()

    nfg = NextflowGenerator(process_connections=pipeline_list,
                            nextflow_file=parsed_output_nf,
                            process_map=process_map,
                            pipeline_name=args.pipeline_name,
                            auto_dependency=args.no_dep,
                            merge_params=args.merge_params,
                            export_params=args.export_params)

    logger.info(colored_print("Building your awesome pipeline..."))

    if args.export_params:
        nfg.export_params()
        sys.exit(0)
    elif args.export_directives:
        nfg.export_directives()
        sys.exit(0)
    elif args.fetch_docker_tags:
        nfg.fetch_docker_tags()
        sys.exit(0)
    else:
        # building the actual pipeline nf file
        nfg.build()

    # copy template to cwd, to allow for immediate execution
    if not args.pipeline_only:
        copy_project(parsed_output_nf)

    logger.info(colored_print("DONE!", "green_bold"))


def inspect(args):

    try:
        nf_inspect = NextflowInspector(args.trace_file, args.refresh_rate,
                                       args.pretty, args.url)
        if args.mode == "overview":
            nf_inspect.display_overview()

        if args.mode == "broadcast":
            nf_inspect.broadcast_status()

    except eh.InspectionError as ie:
        logger.error(colored_print(ie.value, "red_bold"))
        sys.exit(1)

    except eh.LogError as le:
        logger.error(colored_print(le.value, "red_bold"))
        sys.exit(1)



def report(args):

    try:
        fc_report = FlowcraftReport(
            report_file=args.report_file,
            trace_file=args.trace_file,
            log_file=args.log_file,
            watch=args.watch,
            ip_addr=args.url)

        fc_report.broadcast_report()

    except eh.ReportError as re:
        logger.error(colored_print(re.value, "red_bold"))
        sys.exit(1)

    except eh.LogError as le:
        logger.error(colored_print(le.value, "red_bold"))
        sys.exit(1)


def main():

    args = get_args()

    if args.version:
        print(__version__)

    if args.debug:
        logger.setLevel(logging.DEBUG)

        # create formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s')

    else:
        logger.setLevel(logging.INFO)

        # create special formatter for info logs
        formatter = logging.Formatter('%(message)s')

    # create console handler and set level to debug
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG)

    # add formatter to ch
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    if args.main_op == "build":
        build(args)

    if args.main_op == "inspect":
        inspect(args)

    if args.main_op == "report":
        report(args)


if __name__ == '__main__':

    main()
