import re
import pkgutil

try:
    from generator import components
except ImportError:
    from flowcraft.generator import components


def convert_camel_case(name):
    """Convers a CamelCase string into a snake_case one

    Parameters
    ----------
    name : str
        An arbitrary string that may be CamelCase

    Returns
    -------
    str
        The input string converted into snake_case

    """
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()


def collect_process_map():
    """Collects Process classes and return dict mapping templates to classes

    This function crawls through the components module and retrieves all
    classes that inherit from the Process class. Then, it converts the name
    of the classes (which should be CamelCase) to snake_case, which is used
    as the template name.

    Returns
    -------
    dict
        Dictionary mapping the template name (snake_case) to the corresponding
        process class.
    """

    process_map = {}

    prefix = "{}.".format(components.__name__)
    for importer, modname, _ in pkgutil.iter_modules(components.__path__,
                                                     prefix):

        _module = importer.find_module(modname).load_module(modname)

        _component_classes = [
            cls for cls in _module.__dict__.values() if
            isinstance(cls, type) and cls.__name__ != "Process"
        ]

        for cls in _component_classes:
            process_map[convert_camel_case(cls.__name__)] = cls

    return process_map
