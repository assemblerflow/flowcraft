import logging
import sys

logger = logging.getLogger("main.{}".format(__name__))

COLORS = {
    "green_bold": "1;32m",
    "red_bold": "1;31m",
    "white_bold": "1;38m",
    "white_underline": "4;38m",
    "blue_bold": "1;36m",
    "purple_bold": "1;34m"
}


def colored_print(msg, color_label="white_bold"):
    """
    This function enables users to add a color to the print. It also enables
    to pass end_char to print allowing to print several strings in the same line
    in different prints.

    Parameters
    ----------
    color_string: str
        The color code to pass to the function, which enables color change as
        well as background color change.
    msg: str
        The actual text to be printed
    end_char: str
        The character in which each print should finish. By default it will be
        "\n".

    """

    col = COLORS[color_label]

    return "\x1b[{}{}\x1b[0m".format(col, msg)


def procs_dict_parser(procs_dict):
    """
    This function handles the dictionary of attributes of each Process class
    to print to stdout.

    Parameters
    ----------
    procs_dict: dict
        A dictionary with the class attributes used by the argument that prints
        the lists of processes, both for short_list and for detailed_list.


    """

    logger.info(colored_print(
        "\n===== L I S T   O F   P R O C E S S E S =====\n", "green_bold"))

    for template, dict_proc_info in procs_dict.items():
        template_str = "=> {}".format(template)
        logger.info(colored_print(template_str, "blue_bold"))

        for info in dict_proc_info:
            info_str = "{}:".format(info)

            if isinstance(dict_proc_info[info], list):
                if len(dict_proc_info[info]) == 0:
                    arg_msg = "None"
                else:
                    arg_msg = ", ".join(dict_proc_info[info])
            else:
                arg_msg = dict_proc_info[info]

            logger.info("   {} {}".format(
                colored_print(info_str, "white_underline"), arg_msg
            ))


def proc_collector(process_map, args, processes_list=None):
    """
    Function that collects all processes available and stores a dictionary of
    the required arguments of each process class to be passed to
    procs_dict_parser

    Parameters
    ----------
    process_map: dict
        The dictionary with the Processes currently available in flowcraft
        and their corresponding classes as values
    args: argparse.Namespace
        The arguments passed through argparser that will be access to check the
        type of list to be printed
    processes_list: list
        List with all the available processes of a recipe. In case no recipe
        is passed, the list should come empty.

    """

    arguments_list = []

    # prints a detailed list of the process class arguments
    if args.detailed_list:
        # list of attributes to be passed to proc_collector
        arguments_list += [
            "input_type",
            "output_type",
            "description",
            "dependencies",
            "conflicts"
        ]

    # prints a short list with each process and the corresponding description
    if args.short_list:
        arguments_list += [
            "description"
        ]

    if arguments_list:
        # dict to store only the required entries
        procs_dict = {}
        # loops between all process_map Processes
        for name, cls in process_map.items():

            if processes_list:
                # Skip process if process_list is provided and name not in
                # processes list
                if name not in processes_list:
                    continue

            # instantiates each Process class
            cls_inst = cls(template=name)
            d = {arg_key: vars(cls_inst)[arg_key] for arg_key in vars(cls_inst)
                 if arg_key in arguments_list}
            procs_dict[name] = d

        procs_dict_parser(procs_dict)

        sys.exit(0)
