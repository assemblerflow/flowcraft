import logging

logger = logging.getLogger("main.{}".format(__name__))

# Set the tokens used for the main syntax
# Token signaling the start of a fork
FORK_TOKEN = "("
# Token separating different lanes from a fork
LANE_TOKEN = "|"
# Token that closes a fork
CLOSE_TOKEN = ")"


class SanityError(Exception):
    """
    Class to raise a custom error for sanity checks
    """
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)

def insanity_check(pipeline_string):
    """
    A function that performs sanity checks in the pipeline string.
    It will check for:
    - Different number of '(' and ')' characters, which indicates that some
    forks are poorly constructed.
    - Duplicated '|' character between two processes.
    - The existence of processes between the fork starting and ending character
    and the process separator character.
    - The existence of a process between the fork starting character and the
    process separator character.
    - The presence of a process between two fork starting characters.
    - The presence of process separator character '|' within each fork.
    - Duplicated processes within a fork.

    Checks are performed in sequential order, as indicated above and will raise
    an exception `SanityError` each time a string violates any of these checks.
    One can test this function by using the -c option.

    Parameters
    ----------
    pipeline_string: str
        String with the definition of the pipeline, e.g.::
            'processA processB processC(ProcessD | ProcessE)'

    """

    # Gets rid of all spaces in string
    string = pipeline_string.replace(" ", "")

    # FIRST, CHECK if the number of brackets are equal, if not raise a warning
    if string.count("(") != string.count(")"):
        # get the number of each type of bracket and state the one that has a
        # higher value
        dict_values = {"(": string.count("("), ")": string.count(")")}
        max_bracket = max(dict_values, key=dict_values.get)

        raise SanityError("A different number of '(' and ')' was specified."
                          "There are {} extra '{}'. The number of '(' and ')'"
                          "should be equal.".format(
                              str(abs(string.count("(") - string.count(")"))),
                              max_bracket)
                          )

    # SECOND, CHECK for improper use of special characters.
    # this is in fact composed by three if statements.

    # Check for duplicated '|' character.
    if "||" in string:
        raise SanityError("Duplicated fork separator character '|'.")

    # Check for the absence of processes in one of the branches of the fork
    # ['|)' and '(|'] and for the existence of a process before starting a fork
    # (in an inner fork) ['|('].
    if "(|" in string or "|(" in string or "|)" in string:
        raise SanityError("There must be a process between the fork start "
                          "character '(' or end ')' and the separator of "
                          "processes character '|'")

    # Check for duplicated '(' character, which indicate that no process was
    # specified before forking twice.
    if "((" in string:
        raise SanityError("There must be a starting process after the fork "
                          "before adding a new fork. E.g: proc1 ( proc2.1 "
                          "(proc3.1 | proc3.2) | proc 2.2 )")

    # CHECKS INSIDE THE FORKS
    # first lets get all forks to a list.
    list_of_forks = []  # stores forks
    left_indexes = []   # stores indexes of left brackets

    # iterate through the string looking for '(' and ')'.
    for pos, char in enumerate(string):
        if char == "(":
            # saves pos to left_indexes list
            left_indexes.append(pos)
        elif char == ")" and len(left_indexes) > 0:
            # saves fork to list_of_forks
            list_of_forks.append(string[left_indexes[-1] + 1: pos])
            # removes last bracket from left_indexes list
            left_indexes = left_indexes[:-1]

    # sort list in descending order of number of forks
    list_of_forks.sort(key=lambda x: x.count("("), reverse=True)

    # Now, we can iterate through list_of_forks and check for errors in each fork
    for fork in list_of_forks:
        # remove inner forks for these checks since each fork has its own entry
        # in list_of_forks. Note that each fork is now sorted in descending order
        # which enables to remove sequentially the string for the fork potentially with
        # more inner forks
        for subfork in list_of_forks:
            # checks if subfork is contained in fork and if they are different,
            # avoiding to remove itself
            if subfork in list_of_forks and subfork != fork:
                # removes inner forks. Note that string has no spaces
                fork_simplified = fork.replace("({})".format(subfork), "")

        # Checks if there is no fork separator character '|' within each fork
        if not len(fork_simplified.split("|")) > 1:
            raise SanityError("One of the forks doesn't have '|' separator "
                              "between the processes to fork. This is the prime"
                              " suspect: '({})'".format(fork))

        # Check if there is a repeated process within a fork - linked with the above
        if len(fork_simplified.split("|")) != len(
                set(fork_simplified.split("|"))):
            raise SanityError("There are duplicated processes within a fork. "
                              "E.g.: proc1 (proc2.1 | proc2.1 | proc2.2)."
                              "This is the prime suspect: '({})'".format(fork))


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

    # executes sanity checks in pipeline string before parsing it.
    insanity_check(pipeline_str)

    pipeline_links = []
    lane = 1

    # Get number of forks in the pipeline
    nforks = pipeline_str.count(FORK_TOKEN)
    logger.debug("Found {} fork(s)".format(nforks))

    # If there are no forks, connect the pipeline as purely linear
    if not nforks:
        logger.debug("Detected linear pipeline string : {}".format(
            pipeline_str))
        pipeline_links.extend(linear_connection(pipeline_str.split(), lane))
        return pipeline_links

    for i in range(nforks):

        logger.debug("Processing fork {}".format(i))
        # Split the pipeline at each fork start position. fields[-1] will
        # hold the process after the fork. fields[-2] will hold the processes
        # before the fork.
        fields = pipeline_str.split(FORK_TOKEN, i + 1)

        # Get the processes before the fork. This may be empty when the
        # fork is at the beginning of the pipeline.
        previous_process = fields[-2].split()
        logger.debug("Previous processes: {}".format(fields[-2]))
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
                previous_process = ["init"]
                lane = 0

            # Add the linear modules before the fork
            pipeline_links.extend(
                linear_connection(previous_process, lane))

        fork_source = previous_process[-1]
        logger.debug("Fork source is set to: {}".format(fork_source))
        fork_lane = get_source_lane(fork_source, pipeline_links)
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
    fork_process : str
        Process name.
    pipeline_list : list
        List with the pipeline connection dictionaries.

    Returns
    -------
    int
        Lane of the last process that matches fork_process
    """

    for p in pipeline_list[::-1]:
        print(p)
        if p["output"]["process"] == fork_process:
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


def fork_connection(source, sink, fork_lane, lane):
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
    fork_lane : int
        Lane of the forking process
    lane : int
        Lane of the source process

    Returns
    -------
    res : list
        List of dictionaries with the links between processes
    """

    logger.debug("Establishing forking of source '{}' into processes"
                 " '{}'. Fork lane set to '{}' and lane set to "
                 "'{}'".format(source, sink, fork_lane, lane))

    res = []
    # Increase the lane counter for the first lane
    lane_counter = lane + 1

    for p in sink:
        res.append({
            "input": {
                "process": source,
                "lane": fork_lane
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
