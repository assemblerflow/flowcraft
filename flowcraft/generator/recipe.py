try:
    from generator.process_details import colored_print
    import generator.error_handling as eh
    from generator import recipes
except ImportError:
    from flowcraft.generator.process_details import colored_print
    import flowcraft.generator.error_handling as eh
    from flowcraft.generator import recipes

from collections import OrderedDict
import sys
import json
import logging
import pkgutil

logger = logging.getLogger("main.{}".format(__name__))


class InnuendoRecipe:

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


class Innuendo(InnuendoRecipe):
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


def brew_innuendo(args):
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

    # Create recipe class instance
    automatic_pipeline = Innuendo()

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


class Recipe:

    def __init__(self):

        self.pipeline_str = None
        """
        str: The raw pipeline string, with no attribute or directives, except
        for number indicators for when there are duplicate components.
        
        e.g.: "fastqc trimmomatic spades"
        e.g.: "fastqc trimmomatic (spades#1 | spades#2)
        """

        self.directives = {}
        """
        dict: Dictionary with the parameters and directives for each component
        in the pipeline_str attribute. Missing components will be left with
        the default parameters and directives. 
        """

    def brew(self):

        if not hasattr(self, "name"):
            raise eh.RecipeError("Recipe class '{}' does not have a 'name' "
                                 "attribute set".format(self.__class__))

        if not self.pipeline_str:
            raise eh.RecipeError("Recipe with name '{}' does not have a "
                                 "pipeline_str attribute set".format(self.name))

        for component, vals in self.directives.items():

            params = vals.get("params", None)
            directives = vals.get("directives", None)

            # Check for component number symbol
            if "#" in component:
                _component = component.split("#")[0]
            else:
                _component = component

            component_str = self._get_component_str(_component, params,
                                                    directives)

            self.pipeline_str = self.pipeline_str.replace(component,
                                                          component_str)

        return self.pipeline_str

    @staticmethod
    def _get_component_str(component, params=None, directives=None):
        """ Generates a component string based on the provided parameters and
        directives

        Parameters
        ----------
        component : str
            Component name
        params : dict
            Dictionary with parameter information
        directives : dict
            Dictionary with directives information

        Returns
        -------
        str
            Component string with the parameters and directives, ready for
            parsing by flowcraft engine
        """

        final_directives = {}

        if directives:
            final_directives = directives

        if params:
            final_directives["params"] = params

        if final_directives:
            return "{}={}".format(
                component, json.dumps(final_directives, separators=(",", ":")))
        else:
            return component


def brew_recipe(recipe_name):
    """Returns a pipeline string from a recipe name.

    Parameters
    ----------
    recipe_name : str
        Name of the recipe. Must match the name attribute in one of the classes
        defined in :mod:`flowcraft.generator.recipes`

    Returns
    -------
    str
        Pipeline string ready for parsing and processing by flowcraft engine
    """

    # This will iterate over all modules included in the recipes subpackage
    # It will return the import class and the module name, algon with the
    # correct prefix
    prefix = "{}.".format(recipes.__name__)
    for importer, modname, _ in pkgutil.iter_modules(recipes.__path__, prefix):

        # Import the current module
        _module = importer.find_module(modname).load_module(modname)

        # Fetch all available classes in module
        _recipe_classes = [cls for cls in _module.__dict__.values() if
                           isinstance(cls, type)]

        # Iterate over each Recipe class, and check for a match with the
        # provided recipe name.
        for cls in _recipe_classes:
            # Create instance of class to allow fetching the name attribute
            recipe_cls = cls()
            if getattr(recipe_cls, "name", None) == recipe_name:
                return recipe_cls.brew()

    logger.error(
        colored_print("Recipe name '{}' does not exist.".format(recipe_name))
    )
    sys.exit(1)


def list_recipes(full=False):
    """Method that iterates over all available recipes and prints their
    information to the standard output

    Parameters
    ----------
    full : bool
        If true, it will provide the pipeline string along with the recipe name
    """

    logger.info(colored_print(
        "\n===== L I S T   O F   R E C I P E S =====\n",
        "green_bold"))

    # This will iterate over all modules included in the recipes subpackage
    # It will return the import class and the module name, algon with the
    # correct prefix
    prefix = "{}.".format(recipes.__name__)
    for importer, modname, _ in pkgutil.iter_modules(recipes.__path__, prefix):

        # Import the current module
        _module = importer.find_module(modname).load_module(modname)

        # Fetch all available classes in module
        _recipe_classes = [cls for cls in _module.__dict__.values() if
                           isinstance(cls, type)]

        # Iterate over each Recipe class, and check for a match with the
        # provided recipe name.
        for cls in _recipe_classes:

            recipe_cls = cls()

            if hasattr(recipe_cls, "name"):
                logger.info(colored_print("=> {}".format(recipe_cls.name), "blue_bold"))
                if full:
                    logger.info(colored_print("\t {}".format(recipe_cls.__doc__), "purple_bold"))
                    logger.info(colored_print("Pipeline string: {}\n".format(recipe_cls.pipeline_str), "yellow_bold"))

    sys.exit(0)
