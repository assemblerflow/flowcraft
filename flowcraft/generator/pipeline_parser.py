import os
import logging
import re
from difflib import SequenceMatcher

try:
    from generator.error_handling import SanityError
    from generator.process_details import colored_print
except ImportError:
    from flowcraft.generator.error_handling import SanityError
    from flowcraft.generator.process_details import colored_print

logger = logging.getLogger("main.{}".format(__name__))

# Set the tokens used for the main syntax
# Token signaling the start of a fork
FORK_TOKEN = "("
# Token separating different lanes from a fork
LANE_TOKEN = "|"
# Token that closes a fork
CLOSE_TOKEN = ")"


def guess_process(query_str, process_map):
    """
    Function to guess processes based on strings that are not available in
    process_map. If the string has typos and is somewhat similar (50%) to any
    process available in flowcraft it will print info to the terminal,
    suggesting the most similar processes available in flowcraft.

    Parameters
    ----------
    query_str: str
        The string of the process with potential typos
    process_map:
        The dictionary that contains all the available processes

    """

    save_list = []
    # loops between the processes available in process_map
    for process in process_map:
        similarity = SequenceMatcher(None, process, query_str)
        # checks if similarity between the process and the query string is
        # higher than 50%
        if similarity.ratio() > 0.5:
            save_list.append(process)

    # checks if any process is stored in save_list
    if save_list:
        logger.info(colored_print(
            "Maybe you meant:\n\t{}".format("\n\t".join(save_list)), "white"))

    logger.info(colored_print("Hint: check the available processes by using "
                              "the '-l' or '-L' flag.", "white"))


def remove_inner_forks(text):
    """Recursively removes nested brackets

    This function is used to remove nested brackets from fork strings using
    regular expressions

    Parameters
    ----------
    text: str
        The string that contains brackets with inner forks to be removed

    Returns
    -------
    text: str
        the string with only the processes that are not in inner forks, thus
        the processes that belong to a given fork.

    """

    n = 1  # run at least once for one level of fork
    # Then this loop assures that all brackets will get removed in a nested
    # structure
    while n:
        # this removes non-nested brackets
        text, n = re.subn(r'\([^()]*\)', '', text)

    return text

def empty_tasks(p_string):
    """
    Function to check if pipeline string is empty or has an empty string

    Parameters
    ----------
    p_string: str
         String with the definition of the pipeline, e.g.::
             'processA processB processC(ProcessD | ProcessE)'

    """
    if p_string.strip() == "":
        raise SanityError("'-t' parameter received an empty string or "
                          "an empty file.")


def brackets_but_no_lanes(p_string):
    """
    Function to check if a LANE_TOKEN is provided but no fork is initiated.
    Parameters
    ----------
    p_string: str
         String with the definition of the pipeline, e.g.::
             'processA processB processC(ProcessD | ProcessE)'

    """

    if "|" in p_string and "(" not in p_string:
        raise SanityError("No fork initiation character '(' was "
                          "provided but there is a fork lane separator "
                          "character '|'")


def brackets_insanity_check(p_string):
    """
    This function performs a check for different number of '(' and ')'
    characters, which indicates that some forks are poorly constructed.

    Parameters
    ----------
    p_string: str
         String with the definition of the pipeline, e.g.::
             'processA processB processC(ProcessD | ProcessE)'

    """

    if p_string.count(FORK_TOKEN) != p_string.count(CLOSE_TOKEN):
        # get the number of each type of bracket and state the one that has a
        # higher value
        dict_values = {
            FORK_TOKEN: p_string.count(FORK_TOKEN),
            CLOSE_TOKEN: p_string.count(CLOSE_TOKEN)
        }
        max_bracket = max(dict_values, key=dict_values.get)

        raise SanityError(
            "A different number of '(' and ')' was specified. There are "
            "{} extra '{}'. The number of '(' and ')'should be equal.".format(
                str(abs(p_string.count(FORK_TOKEN) - p_string.count(CLOSE_TOKEN))),
                max_bracket))


def lane_char_insanity_check(p_string):
    """
    This function performs a sanity check for multiple '|' character
    between two processes.

    Parameters
    ----------
    p_string: str
         String with the definition of the pipeline, e.g.::
             'processA processB processC(ProcessD | ProcessE)'

    """

    if LANE_TOKEN + LANE_TOKEN in p_string:
        raise SanityError("Duplicated fork separator character '|'.")


def final_char_insanity_check(p_string):
    """
    This function checks if lane token is the last element of the pipeline
    string.

    Parameters
    ----------
    p_string: str
         String with the definition of the pipeline, e.g.::
             'processA processB processC(ProcessD | ProcessE)'

    """

    # Check if last character of string is a LANE_TOKEN
    if p_string.endswith(LANE_TOKEN):
        raise SanityError("Fork separator character '|' cannot be the "
                          "last element of pipeline string")


def fork_procs_insanity_check(p_string):
    """
    This function checks if the pipeline string contains a process between
    the fork start token or end token and the separator (lane) token. Checks for
    the absence of processes in one of the branches of the fork ['|)' and '(|']
    and for the existence of a process before starting a fork (in an inner fork)
    ['|('].

    Parameters
    ----------
    p_string: str
         String with the definition of the pipeline, e.g.::
             'processA processB processC(ProcessD | ProcessE)'

    """

    # Check for the absence of processes in one of the branches of the fork
    # ['|)' and '(|'] and for the existence of a process before starting a fork
    # (in an inner fork) ['|('].
    if FORK_TOKEN + LANE_TOKEN in p_string or \
            LANE_TOKEN + CLOSE_TOKEN in p_string or \
            LANE_TOKEN + FORK_TOKEN in p_string:
        raise SanityError("There must be a process between the fork "
                          "start character '(' or end ')' and the separator of "
                          "processes character '|'")


def start_proc_insanity_check(p_string):
    """
    This function checks if there is a starting process after the beginning of
    each fork. It checks for duplicated start tokens ['(('].

    Parameters
    ----------
    p_string: str
         String with the definition of the pipeline, e.g.::
             'processA processB processC(ProcessD | ProcessE)'

    """

    if FORK_TOKEN + FORK_TOKEN in p_string:
        raise SanityError("There must be a starting process after the "
                          "fork before adding a new fork. E.g: proc1 ( proc2.1 "
                          "(proc3.1 | proc3.2) | proc 2.2 )")


def late_proc_insanity_check(p_string):
    """
    This function checks if there are processes after the close token. It
    searches for everything that isn't "|" or ")" after a ")" token.

    Parameters
    ----------
    p_string: str
         String with the definition of the pipeline, e.g.::
             'processA processB processC(ProcessD | ProcessE)'

    """

    if re.search('\{}[^|)]'.format(CLOSE_TOKEN), p_string):
        raise SanityError("After a fork it is not allowed to have any "
                          "alphanumeric value.")


def inner_fork_insanity_checks(pipeline_string):
    """
    This function performs two sanity checks in the pipeline string. The first
    check, assures that each fork contains a lane token '|', while the second
    check looks for duplicated processes within the same fork.

    Parameters
    ----------
    pipeline_string: str
         String with the definition of the pipeline, e.g.::
             'processA processB processC(ProcessD | ProcessE)'

    """

    # first lets get all forks to a list.
    list_of_forks = []  # stores forks
    left_indexes = []   # stores indexes of left brackets

    # iterate through the string looking for '(' and ')'.
    for pos, char in enumerate(pipeline_string):
        if char == FORK_TOKEN:
            # saves pos to left_indexes list
            left_indexes.append(pos)
        elif char == CLOSE_TOKEN and len(left_indexes) > 0:
            # saves fork to list_of_forks
            list_of_forks.append(pipeline_string[left_indexes[-1] + 1: pos])
            # removes last bracket from left_indexes list
            left_indexes = left_indexes[:-1]

    # sort list in descending order of number of forks
    list_of_forks.sort(key=lambda x: x.count(FORK_TOKEN), reverse=True)

    # Now, we can iterate through list_of_forks and check for errors in each
    # fork
    for fork in list_of_forks:
        # remove inner forks for these checks since each fork has its own entry
        # in list_of_forks. Note that each fork is now sorted in descending
        # order which enables to remove sequentially the string for the fork
        # potentially with more inner forks
        for subfork in list_of_forks:
            # checks if subfork is contained in fork and if they are different,
            # avoiding to remove itself
            if subfork in list_of_forks and subfork != fork:
                # removes inner forks. Note that string has no spaces
                fork_simplified = fork.replace("({})".format(subfork), "")
            else:
                fork_simplified = fork

        # Checks if there is no fork separator character '|' within each fork
        if not len(fork_simplified.split(LANE_TOKEN)) > 1:
            raise SanityError("One of the forks doesn't have '|' "
                              "separator between the processes to fork. This is"
                              " the prime suspect: '({})'".format(fork))


def insanity_checks(pipeline_str):
    """Wrapper that performs all sanity checks on the pipeline string

    Parameters
    ----------
    pipeline_str : str
        String with the pipeline definition
    """

    # Gets rid of all spaces in string
    p_string = pipeline_str.replace(" ", "").strip()

    # some of the check functions use the pipeline_str as the user provided but
    # the majority uses the parsed p_string.
    checks = [
        [p_string, [
            empty_tasks,
            brackets_but_no_lanes,
            brackets_insanity_check,
            lane_char_insanity_check,
            final_char_insanity_check,
            fork_procs_insanity_check,
            start_proc_insanity_check,
            late_proc_insanity_check
        ]],
        [pipeline_str, [
            inner_fork_insanity_checks
        ]]
    ]

    # executes sanity checks in pipeline string before parsing it.
    for param, func_list in checks:
        for func in func_list:
            func(param)


def parse_pipeline(pipeline_str):
    """Parses a pipeline string into a dictionary with the connections between
    process

    Parameters
    ----------
    pipeline_str : str
        String with the definition of the pipeline, e.g.::
            'processA processB processC(ProcessD | ProcessE)'

    Returns
    -------
    pipeline_links : list

    """

    if os.path.exists(pipeline_str):
        logger.debug("Found pipeline file: {}".format(pipeline_str))
        with open(pipeline_str) as fh:
            pipeline_str = "".join([x.strip() for x in fh.readlines()])

    logger.info(colored_print("Resulting pipeline string:\n"))
    logger.info(colored_print(pipeline_str + "\n"))

    # Perform pipeline insanity checks
    insanity_checks(pipeline_str)

    logger.debug("Parsing pipeline string: {}".format(pipeline_str))

    pipeline_links = []
    lane = 1

    # Get number of forks in the pipeline
    nforks = pipeline_str.count(FORK_TOKEN)
    logger.debug("Found {} fork(s)".format(nforks))

    # If there are no forks, connect the pipeline as purely linear
    if not nforks:
        logger.debug("Detected linear pipeline string : {}".format(
            pipeline_str))
        linear_pipeline = ["__init__"] + pipeline_str.split()
        pipeline_links.extend(linear_connection(linear_pipeline, lane))
        return pipeline_links

    for i in range(nforks):

        logger.debug("Processing fork {} in lane {}".format(i, lane))
        # Split the pipeline at each fork start position. fields[-1] will
        # hold the process after the fork. fields[-2] will hold the processes
        # before the fork.
        fields = pipeline_str.split(FORK_TOKEN, i + 1)

        # Get the processes before the fork. This may be empty when the
        # fork is at the beginning of the pipeline.
        previous_process = fields[-2].split(LANE_TOKEN)[-1].split()
        logger.debug("Previous processes string: {}".format(fields[-2]))
        logger.debug("Previous processes list: {}".format(previous_process))
        # Get lanes after the fork
        next_lanes = get_lanes(fields[-1])
        logger.debug("Next lanes object: {}".format(next_lanes))
        # Get the immediate targets of the fork
        fork_sink = [x[0] for x in next_lanes]
        logger.debug("The fork sinks into the processes: {}".format(fork_sink))

        # The first fork is a special case, where the processes before AND
        # after the fork (until the start of another fork) are added to
        # the ``pipeline_links`` variable. Otherwise, only the processes
        # after the fork will be added
        if i == 0:
            # If there are no previous process, the fork is at the beginning
            # of the pipeline string. In this case, inject the special
            # "init" process.
            if not previous_process:
                previous_process = ["__init__"]
                lane = 0
            else:
                previous_process = ["__init__"] + previous_process

            # Add the linear modules before the fork
            pipeline_links.extend(
                linear_connection(previous_process, lane))

        fork_source = previous_process[-1]
        logger.debug("Fork source is set to: {}".format(fork_source))
        fork_lane = get_source_lane(previous_process, pipeline_links)
        logger.debug("Fork lane is set to: {}".format(fork_lane))
        # Add the forking modules
        pipeline_links.extend(
            fork_connection(fork_source, fork_sink, fork_lane, lane))
        # Add the linear connections in the subsequent lanes
        pipeline_links.extend(
            linear_lane_connection(next_lanes, lane))

        lane += len(fork_sink)

    return pipeline_links


def get_source_lane(fork_process, pipeline_list):
    """Returns the lane of the last process that matches fork_process

    Parameters
    ----------
    fork_process : list
        List of processes before the fork.
    pipeline_list : list
        List with the pipeline connection dictionaries.

    Returns
    -------
    int
        Lane of the last process that matches fork_process
    """

    fork_source = fork_process[-1]
    fork_sig = [x for x in fork_process if x != "__init__"]

    for position, p in enumerate(pipeline_list[::-1]):

        if p["output"]["process"] == fork_source:

            lane = p["output"]["lane"]
            logger.debug("Possible source match found in position {} in lane"
                         " {}".format(position, lane))
            lane_sequence = [x["output"]["process"] for x in pipeline_list
                             if x["output"]["lane"] == lane]
            logger.debug("Testing lane sequence '{}' against fork signature"
                         " '{}'".format(lane_sequence, fork_sig))
            if lane_sequence == fork_sig:
                return p["output"]["lane"]

    return 0


def get_lanes(lanes_str):
    """From a raw pipeline string, get a list of lanes from the start
    of the current fork.

    When the pipeline is being parsed, it will be split at every fork
    position. The string at the right of the fork position will be provided
    to this function. It's job is to retrieve the lanes that result
    from that fork, ignoring any nested forks.

    Parameters
    ----------
    lanes_str : str
        Pipeline string after a fork split

    Returns
    -------
    lanes : list
        List of lists, with the list of processes for each lane

    """

    logger.debug("Parsing lanes from raw string: {}".format(lanes_str))

    # Temporarily stores the lanes string after removal of nested forks
    parsed_lanes = ""
    # Flag used to determined whether the cursor is inside or outside the
    # right fork
    infork = 0

    for i in lanes_str:

        # Nested fork started
        if i == FORK_TOKEN:
            infork += 1
        # Nested fork stopped
        if i == CLOSE_TOKEN:
            infork -= 1

        if infork < 0:
            break

        # Save only when in the right fork
        if infork == 0:
            # Ignore forking syntax tokens
            if i not in [FORK_TOKEN, CLOSE_TOKEN]:
                parsed_lanes += i

    return [x.split() for x in parsed_lanes.split(LANE_TOKEN)]


def linear_connection(plist, lane):
    """Connects a linear list of processes into a list of dictionaries

    Parameters
    ----------
    plist : list
        List with process names. This list should contain at least two entries.
    lane : int
        Corresponding lane of the processes

    Returns
    -------
    res : list
        List of dictionaries with the links between processes
    """

    logger.debug(
        "Establishing linear connection with processes: {}".format(plist))

    res = []
    previous = None

    for p in plist:
        # Skip first process
        if not previous:
            previous = p
            continue

        res.append({
            "input": {
                "process": previous,
                "lane": lane
            },
            "output": {
                "process": p,
                "lane": lane
            }
        })
        previous = p

    return res


def fork_connection(source, sink, source_lane, lane):
    """Makes the connection between a process and the first processes in the
    lanes to wich it forks.

    The ``lane`` argument should correspond to the lane of the source process.
    For each lane in ``sink``, the lane counter will increase.

    Parameters
    ----------
    source : str
        Name of the process that is forking
    sink : list
        List of the processes where the source will fork to. Each element
        corresponds to the start of a lane.
    source_lane : int
        Lane of the forking process
    lane : int
        Lane of the source process

    Returns
    -------
    res : list
        List of dictionaries with the links between processes
    """

    logger.debug("Establishing forking of source '{}' into processes"
                 " '{}'. Source lane set to '{}' and lane set to '{}'".format(
                    source, sink, source_lane, lane))

    res = []
    # Increase the lane counter for the first lane
    lane_counter = lane + 1

    for p in sink:
        res.append({
            "input": {
                "process": source,
                "lane": source_lane
            },
            "output": {
                "process": p,
                "lane": lane_counter
            }
        })
        lane_counter += 1

    return res


def linear_lane_connection(lane_list, lane):
    """

    Parameters
    ----------
    lane_list : list
        Each element should correspond to a list of processes for a given lane
    lane : int
        Lane counter before the fork start

    Returns
    -------
    res : list
        List of dictionaries with the links between processes
    """

    logger.debug(
        "Establishing linear connections for lanes: {}".format(lane_list))

    res = []
    # Increase the lane counter for the first lane
    lane += 1

    for l in lane_list:
        res.extend(linear_connection(l, lane))
        lane += 1

    return res
