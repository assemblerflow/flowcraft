try:
    from generator.process_details import colored_print
except ImportError:
    from flowcraft.generator.process_details import colored_print

from collections import OrderedDict
import sys
import logging

logger = logging.getLogger("main.{}".format(__name__))


class Recipe:

    def __init__(self):
        """Class to build automatic pipelines based on the processes provided.

        This class provides the methods to build the most eficient pipeline
        based on the processes provided. It automatic creates the
        flowcraft pipeline string based on the relationships between the
        possible processes.

        """

        self.count_forks = 0
        """
        int : counts the total possible number of forks
        """

        self.forks = []
        """
        list : a list with all the possible forks
        """

        self.pipeline_string = ""
        """
        str : the generated pipeline string
        """

        self.process_to_id = {}
        """
        dict: key value between the process name and its identifier
        """

        self.process_descriptions = {}

    @staticmethod
    def validate_pipeline(pipeline_string):
        """Validate pipeline string

        Validates the pipeline string by searching for forbidden characters

        Parameters
        ----------
        pipeline_string : str
            STring with the processes provided

        Returns
        -------

        """
        if "(" in pipeline_string or ")" in pipeline_string or "|" in \
                pipeline_string:
            logger.error(
                colored_print("Please provide a valid task list!", "red_bold")
            )
            return False

        return True

    def build_upstream(self, process_descriptions, task, all_tasks,
                       task_pipeline,
                       count_forks, total_tasks, forks):
        """Builds the upstream pipeline of the current process

        Checks for the upstream processes to the current process and
        adds them to the current pipeline fragment if they were provided in
        the process list.

        Parameters
        ----------
        process_descriptions : dict
            Information of processes input, output and if is forkable
        task : str
            Current process
        all_tasks : list
            A list of all provided processes
        task_pipeline : list
            Current pipeline fragment
        count_forks : int
            Current number of forks
        total_tasks : str
            All space separated processes
        forks : list
            Current forks
        Returns
        -------
        list : resulting pipeline fragment
        """
        if task in process_descriptions:
            if process_descriptions[task][1] is not None:
                if len(process_descriptions[task][1].split("|")) > 1:
                    local_forks = process_descriptions[task][1].split("|")

                    # Produces a new pipeline fragment for each forkable
                    #  process
                    for local_fork in local_forks:
                        if local_fork in total_tasks:
                            count_forks += 1
                            task_pipeline.insert(
                                0,
                                process_descriptions[task][1]
                            )
                            self.define_pipeline_string(
                                process_descriptions,
                                local_fork,
                                False,
                                True,
                                count_forks,
                                total_tasks,
                                forks
                            )

                    return task_pipeline
                else:
                    # Adds the process to the pipeline fragment in case it is
                    # provided in the task list
                    if process_descriptions[task][1] in total_tasks:
                        task_pipeline.insert(
                            0,
                            process_descriptions[task][1].split("|")[0]
                        )

                        # Proceeds building upstream until the input for a
                        # process is None
                        self.build_upstream(
                            process_descriptions,
                            process_descriptions[task][1].split("|")[0],
                            all_tasks,
                            task_pipeline,
                            count_forks,
                            total_tasks,
                            forks
                        )
                    else:
                        logger.error(
                            colored_print("{} not in provided protocols as "
                                          "input for {}".format(
                                process_descriptions[task][1], task), "red_bold"
                            )
                        )

                        sys.exit()

                    return task_pipeline
            else:
                return task_pipeline

    def build_downstream(self, process_descriptions, task, all_tasks,
                         task_pipeline,
                         count_forks, total_tasks, forks):
        """Builds the downstream pipeline of the current process

        Checks for the downstream processes to the current process and
        adds them to the current pipeline fragment.

        Parameters
        ----------
        process_descriptions : dict
            Information of processes input, output and if is forkable
        task : str
            Current process
        all_tasks : list
            A list of all provided processes
        task_pipeline : list
            Current pipeline fragment
        count_forks : int
            Current number of forks
        total_tasks : str
            All space separated processes
        forks : list
            Current forks
        Returns
        -------
        list : resulting pipeline fragment
        """

        if task in process_descriptions:
            if process_descriptions[task][2] is not None:
                if len(process_descriptions[task][2].split("|")) > 1:
                    local_forks = process_descriptions[task][2].split("|")

                    # Adds the process to the pipeline fragment downstream
                    # and defines a new pipeline fragment for each fork.
                    # Those will only look for downstream processes
                    for local_fork in local_forks:
                        if local_fork in total_tasks:
                            count_forks += 1
                            task_pipeline.append(process_descriptions[task][2])
                            self.define_pipeline_string(
                                process_descriptions,
                                local_fork,
                                False,
                                True,
                                count_forks,
                                total_tasks,
                                forks
                            )

                    return task_pipeline
                else:
                    if process_descriptions[task][2] in total_tasks:
                        task_pipeline.append(process_descriptions[task][2].split("|")[0])

                        # Proceeds building downstream until the output for a
                        # process is None
                        self.build_downstream(
                            process_descriptions,
                            process_descriptions[task][2].split("|")[0],
                            all_tasks,
                            task_pipeline,
                            count_forks,
                            total_tasks,
                            forks
                        )

                    return task_pipeline
            else:
                return task_pipeline

    def define_pipeline_string(self, process_descriptions, tasks,
                               check_upstream,
                               check_downstream, count_forks, total_tasks,
                               forks):
        """Builds the possible forks and connections between the provided
        processes

        This method loops through all the provided tasks and builds the
        upstream and downstream pipeline if required. It then returns all
        possible forks than need to be merged à posteriori`

        Parameters
        ----------
        process_descriptions : dict
            Information of processes input, output and if is forkable
        tasks : str
            Space separated processes
        check_upstream : bool
            If is to build the upstream pipeline of the current task
        check_downstream : bool
            If is to build the downstream pipeline of the current task
        count_forks : int
            Number of current forks
        total_tasks : str
            All space separated processes
        forks : list
            Current forks

        Returns
        -------
        list : List with all the possible pipeline forks
        """

        tasks_array = tasks.split()

        for task_unsplit in tasks_array:
            task = task_unsplit.split("=")[0]

            if task not in process_descriptions.keys():
                logger.error(
                    colored_print(
                        "{} not in the possible processes".format(task),
                        "red_bold"
                    )
                )

                sys.exit()
            else:
                process_split = task_unsplit.split("=")

                if len(process_split) > 1:
                    self.process_to_id[process_split[0]] = process_split[1]

            # Only uses the process if it is not already in the possible forks
            if not bool([x for x in forks if task in x]) and not bool([y for y in forks if process_descriptions[task][2] in y]):
                task_pipeline = []

                if task in process_descriptions:

                    if check_upstream:
                        task_pipeline = self.build_upstream(
                            process_descriptions,
                            task,
                            tasks_array,
                            task_pipeline,
                            count_forks,
                            total_tasks,
                            forks
                        )

                    task_pipeline.append(task)

                    if check_downstream:
                        task_pipeline = self.build_downstream(
                            process_descriptions,
                            task,
                            tasks_array,
                            task_pipeline,
                            count_forks,
                            total_tasks,
                            forks
                        )

                # Adds the pipeline fragment to the list of possible forks
                forks.append(list(OrderedDict.fromkeys(task_pipeline)))

            # Checks for task in fork. Case order of input processes is reversed
            elif bool([y for y in forks if process_descriptions[task][2] in y]):
                for fork in forks:
                    if task not in fork:
                        try:
                            dependent_index = fork.index(process_descriptions[task][2])
                            fork.insert(dependent_index, task)
                        except ValueError:
                            continue

        for i in range(0, len(forks)):
            for j in range(0, len(forks[i])):
                try:
                    if len(forks[i][j].split("|")) > 1:
                        forks[i][j] = forks[i][j].split("|")
                        tmp_fork = []
                        for s in forks[i][j]:
                            if s in total_tasks:
                                tmp_fork.append(s)

                        forks[i][j] = tmp_fork

                except AttributeError as e:
                    continue

        return forks

    def build_pipeline_string(self, forks):
        """Parses, filters and merge all possible pipeline forks into the
        final pipeline string

        This method checks for shared start and end sections between forks
        and merges them according to the shared processes::

            [[spades, ...], [skesa, ...], [...,[spades, skesa]]]
                -> [..., [[spades, ...], [skesa, ...]]]

        Then it defines the pipeline string by replacing the arrays levels
        to the flowcraft fork format::

            [..., [[spades, ...], [skesa, ...]]]
                -> ( ... ( spades ... | skesa ... ) )

        Parameters
        ----------
        forks : list
            List with all the possible pipeline forks.

        Returns
        -------
        str : String with the pipeline definition used as input for
        parse_pipeline
        """

        final_forks = []

        for i in range(0, len(forks)):
            needs_merge = [False, 0, 0, 0, 0, ""]
            is_merged = False
            for i2 in range(0, len(forks[i])):
                for j in range(i, len(forks)):
                    needs_merge[0] = False
                    for j2 in range(0, len(forks[j])):
                        try:
                            j2_fork = forks[j][j2].split("|")
                        except AttributeError:
                            j2_fork = forks[j][j2]

                        # Gets the indexes of the forks matrix that need to
                        # be merged
                        if forks[i][i2] in j2_fork and (i2 == 0 or j2 == 0) and i != j:
                            needs_merge[0] = True
                            needs_merge[1] = i
                            needs_merge[2] = i2
                            needs_merge[3] = j
                            needs_merge[4] = j2
                            needs_merge[5] = forks[i][i2]

                    if needs_merge[0]:
                        index_merge_point = forks[needs_merge[3]][-1].index(needs_merge[5])

                        # Merges the forks. If only one fork is possible,
                        # that fork is neglected and it merges into a single
                        # channel.
                        if needs_merge[2] == 0:
                            if len(forks[needs_merge[3]][-1]) < 2:
                                forks[needs_merge[3]] = forks[needs_merge[3]][:-1] + forks[needs_merge[1]][::]
                            else:
                                forks[needs_merge[3]][-1][index_merge_point] = forks[needs_merge[1]]

                        elif needs_merge[4] == 0:
                            if len(forks[needs_merge[3]][-1]) < 2:
                                forks[needs_merge[3]] = forks[needs_merge[3]][:-1] + forks[needs_merge[1]][::]
                            else:
                                forks[needs_merge[3]][-1][index_merge_point] = forks[needs_merge[1]]

                        is_merged = True

            # Adds forks that dont need merge to the final forks
            if needs_merge[0] is not None and not is_merged:
                if bool([nf for nf in forks[i] if "|" in nf]):
                    continue
                final_forks.append(forks[i])

        if len(final_forks) == 1:
            final_forks = str(final_forks[0])

        # parses the string array to the flowcraft nomenclature
        pipeline_string = " " + str(final_forks)\
            .replace("[[", "( ")\
            .replace("]]", " )")\
            .replace("]", " |")\
            .replace(", [", " ")\
            .replace("'", "")\
            .replace(",", "")\
            .replace("[", "")

        if pipeline_string[-1] == "|":
            pipeline_string = pipeline_string[:-1]

        to_search = " {} "
        to_replace = " {}={} "

        # Replace only names by names + process ids
        for key, val in self.process_to_id.items():
            # Case only one process in the pipeline
            pipeline_string = pipeline_string\
                .replace(to_search.format(key),
                         to_replace.format(key, val))

        return pipeline_string

    def run_auto_pipeline(self, tasks):
        """Main method to run the automatic pipeline creation

        This method aggregates the functions required to build the pipeline
        string that can be used as input for the workflow generator.

        Parameters
        ----------
        tasks : str
            A string with the space separated tasks to be included in the
            pipeline

        Returns
        -------
        str : String with the pipeline definition used as input for
        parse_pipeline
        """

        self.forks = self.define_pipeline_string(
            self.process_descriptions,
            tasks,
            True,
            True,
            self.count_forks,
            tasks,
            self.forks
        )

        self.pipeline_string = self.build_pipeline_string(self.forks)

        return self.pipeline_string

    # def get_process_info(self):
    #     return list(self.process_descriptions.keys())


class Innuendo(Recipe):
    """
    Recipe class for the INNUENDO Project. It has all the available in the
    platform for quick use of the processes in the scope of the project.
    """

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

        # The description of the processes
        # [forkable, input_process, output_process]
        self.process_descriptions = {
            "reads_download": [False, None,"integrity_coverage|seq_typing|patho_typing"],
            "patho_typing": [True, None, None],
            "seq_typing": [True, None, None],
            "integrity_coverage": [True, None, "fastqc_trimmomatic"],
            "fastqc_trimmomatic": [False, "integrity_coverage",
                                   "true_coverage"],
            "true_coverage": [False, "fastqc_trimmomatic",
                              "fastqc"],
            "fastqc": [False, "true_coverage", "check_coverage"],
            "check_coverage": [False, "fastqc", "spades"],
            "spades": [False, "fastqc_trimmomatic", "process_spades"],
            "process_spades": [False, "spades", "assembly_mapping"],
            "assembly_mapping": [False, "process_spades", "pilon"],
            "pilon": [False, "assembly_mapping", "mlst"],
            "mlst": [False, "pilon", "abricate|prokka|chewbbaca|sistr"],
            "sistr": [True, "mlst", None],
            "abricate": [True, "mlst", None],
            #"prokka": [True, "mlst", None],
            "chewbbaca": [True, "mlst", None]
        }


def brew_recipe(args, available_recipes):
    """Brews a given list of processes according to the recipe

    Parameters
    ----------
    args : argparse.Namespace
        The arguments passed through argparser that will be used to check the
        the recipe, tasks and brew the process

    Returns
    -------
    str
        The final pipeline string, ready for the engine.
    list
        List of process strings.
    """

    # Exit if recipe does not exist
    if args.recipe not in available_recipes:
        logger.error(
            colored_print("Please provide a recipe to use in automatic "
                          "mode.", "red_bold"))
        sys.exit(1)

    # Create recipe class instance
    automatic_pipeline = available_recipes[args.recipe]()

    if not args.tasks:
        input_processes = " ".join(
            automatic_pipeline.process_descriptions.keys())
    else:
        input_processes = args.tasks

    # Validate the provided pipeline processes
    validated = automatic_pipeline.validate_pipeline(input_processes)
    if not validated:
        sys.exit(1)
    # Get the final pipeline string
    pipeline_string = automatic_pipeline.run_auto_pipeline(input_processes)

    return pipeline_string


# A dictionary of quick recipes
available_recipes = {
    "innuendo": Innuendo,
    "plasmids": "integrity_coverage fastqc_trimmomatic (spades pilon "
              "(mash_dist | abricate) | mash_screen | mapping_patlas)",
    "plasmids_mapping": "integrity_coverage fastqc_trimmomatic mapping_patlas",
    "plasmids_assembly": "integrity_coverage fastqc_trimmomatic (spades pilon"
                         " mash_dist)",
    "plasmids_mash": "integrity_coverage fastqc_trimmomatic mash_screen",
    "den-im": "integrity_coverage fastqc_trimmomatic filter_poly remove_host bowtie retrieve_mapped viral_assembly assembly_mapping pilon split_assembly (dengue_typing | mafft raxml)",
}
