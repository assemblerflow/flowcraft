import pytest
import pkgutil

from argparse import Namespace

from flowcraft.generator import error_handling as eh
from flowcraft.generator import recipes
from flowcraft.generator import recipe


def test_empty_recipe():

    r = recipe.Recipe()

    with pytest.raises(eh.RecipeError):
        r.brew()


def test_empty_pipeline_str():

    r = recipe.Recipe()

    r.name = "teste"

    with pytest.raises(eh.RecipeError):
        r.brew()


def test_basic_recipe():

    r = recipe.Recipe()

    r.name = "teste"
    r.pipeline_str = "teste"

    r.brew()

    assert r.pipeline_str == "teste"


def test_recipe_wdirectives():

    r = recipe.Recipe()

    r.name = "teste"
    r.pipeline_str = "componentA"
    r.directives = {
        "componentA": {
            "params": {
                "paramA": "val"
            },
            "directives": {
                "dirA": "val"
            }
        }
    }

    r.brew()

    assert '"params":{"paramA":"val"}' in r.pipeline_str and \
        '"dirA":"val"' in r.pipeline_str


def test_recipe_partial_directives():

    r = recipe.Recipe()

    r.name = "teste"
    r.pipeline_str = "componentA"
    r.directives = {
        "componentA": {
            "params": {
                "paramA": "val"
            },
        }
    }

    r.brew()

    assert '"params":{"paramA":"val"}' in r.pipeline_str


def test_recipe_partial_directives2():

    r = recipe.Recipe()

    r.name = "teste"
    r.pipeline_str = "componentA"
    r.directives = {
        "componentA": {
            "directives": {
                "dirA": "val"
            }
        }
    }

    r.brew()

    assert '"dirA":"val"' in r.pipeline_str


def test_component_str():

    r = recipe.Recipe()

    r.name = "teste"
    r.pipeline_str = "componentA"
    directives = {
        "dirA": "val"
    }

    res = r._get_component_str("componentA", directives=directives)

    assert '"dirA":"val"' in res


def test_component_str2():

    r = recipe.Recipe()

    r.name = "teste"
    r.pipeline_str = "componentA"
    directives = {
        "paramA": "val"
    }

    res = r._get_component_str("componentA", params=directives)
    print(res)

    assert '"params":{"paramA":"val"}' in res


def test_component_str3():

    r = recipe.Recipe()

    r.name = "teste"
    r.pipeline_str = "componentA"
    params = {
        "paramA": "val"
    }
    directives = {
        "dirA": "val"
    }

    res = r._get_component_str("componentA", params=params,
                               directives=directives)

    assert '"params":{"paramA":"val"}' in res and \
           '"dirA":"val"' in res


def test_brew_recipe():

    res = recipe.brew_recipe("innuca")

    assert res != ""


def test_bad_recipe_name():

    with pytest.raises(SystemExit):
        res = recipe.brew_recipe("bad_name")


def test_all_recipes():

    prefix = "{}.".format(recipes.__name__)
    for importer, modname, _ in pkgutil.iter_modules(recipes.__path__, prefix):

        _module = importer.find_module(modname).load_module(modname)

        _recipe_classes = [cls for cls in _module.__dict__.values() if
                           isinstance(cls, type)]

        for cls in _recipe_classes:
            cls()


def test_innuendo_recipe():

    args = Namespace(tasks=None)

    recipe.brew_innuendo(args)


def test_innuendo_partial_recipe():

    args = Namespace(tasks="integrity_coverage")

    recipe.brew_innuendo(args)


def test_list_recipes():

    with pytest.raises(SystemExit):
        recipe.list_recipes()

def test_list_recipes_full():

    with pytest.raises(SystemExit):
        recipe.list_recipes(True)