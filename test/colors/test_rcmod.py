import unittest
from unittest import TestCase
import numpy as np
import matplotlib as mpl
from distutils.version import LooseVersion
import nose
import matplotlib.pyplot as plt
import nose.tools as nt

from light.colors import rcmod


class RCParamTester(TestCase):

    def flatten_list(self, orig_list):

        iter_list = map(np.atleast_1d, orig_list)
        flat_list = [item for sublist in iter_list for item in sublist]
        return flat_list

    def assert_rc_params(self, params):

        for k, v in params.items():
            if k == "svg.embed_char_paths":
                # This param causes test issues and is deprecated anyway
                continue
            elif isinstance(v, np.ndarray):
                self.assertCountEqual(mpl.rcParams[k], v)
            else:
                self.assertEqual((k, mpl.rcParams[k]), (k, v))


class TestAxesStyle(RCParamTester):

    styles = ["white", "dark", "whitegrid", "darkgrid", "ticks"]

    def test_default_return(self):

        current = rcmod.axes_style()
        self.assert_rc_params(current)

    def test_key_usage(self):

        _style_keys = set(rcmod._style_keys)
        for style in self.styles:
            nt.assert_true(not set(rcmod.axes_style(style)) ^ _style_keys)

    def test_bad_style(self):

        with self.assertRaises(ValueError):
            rcmod.axes_style("i_am_not_a_style")

    def test_rc_override(self):

        rc = {"axes.facecolor": "blue", "foo.notaparam": "bar"}
        out = rcmod.axes_style("darkgrid", rc)
        self.assertEqual(out["axes.facecolor"], "blue")
        self.assertNotIn("foo.notaparam", out)

    def test_set_style(self):

        for style in self.styles:

            style_dict = rcmod.axes_style(style)
            rcmod.set_style(style)
            self.assert_rc_params(style_dict)

    def test_style_context_manager(self):

        rcmod.set_style("darkgrid")
        orig_params = rcmod.axes_style()
        with rcmod.axes_style("whitegrid"):
            context_params = rcmod.axes_style("whitegrid")
            self.assert_rc_params(context_params)
        self.assert_rc_params(orig_params)

    def test_style_context_independence(self):

        self.assertTrue(set(rcmod._style_keys) ^ set(rcmod._context_keys))

    def test_set_rc(self):

        rcmod.set(rc={"lines.linewidth": 4})
        self.assertEqual(mpl.rcParams["lines.linewidth"], 4)
        rcmod.set()

    def test_reset_defaults(self):

        # Changes to the rc parameters make this test hard to manage
        # on older versions of matplotlib, so we'll skip it
        if LooseVersion(mpl.__version__) < LooseVersion("1.3"):
            raise nose.SkipTest

        rcmod.reset_defaults()
        self.assert_rc_params(mpl.rcParamsDefault)
        rcmod.set()

    def test_reset_orig(self):

        # Changes to the rc parameters make this test hard to manage
        # on older versions of matplotlib, so we'll skip it
        if LooseVersion(mpl.__version__) < LooseVersion("1.3"):
            raise nose.SkipTest

        rcmod.reset_orig()
        self.assert_rc_params(mpl.rcParamsOrig)
        rcmod.set()


class TestPlottingContext(RCParamTester):

    contexts = ["paper", "notebook", "talk", "poster"]

    def test_default_return(self):

        current = rcmod.plotting_context()
        self.assert_rc_params(current)

    def test_key_usage(self):

        _context_keys = set(rcmod._context_keys)
        for context in self.contexts:
            missing = set(rcmod.plotting_context(context)) ^ _context_keys
            self.assertTrue(not missing)

    def test_bad_context(self):

        with self.assertRaises(ValueError):
            rcmod.plotting_context("i_am_not_a_context")

    def test_font_scale(self):

        notebook_ref = rcmod.plotting_context("notebook")
        notebook_big = rcmod.plotting_context("notebook", 2)

        font_keys = ["axes.labelsize", "axes.titlesize", "legend.fontsize",
                     "xtick.labelsize", "ytick.labelsize"]

        for k in font_keys:
            self.assertEqual(notebook_ref[k] * 2, notebook_big[k])

    def test_rc_override(self):

        key, val = "grid.linewidth", 5
        rc = {key: val, "foo": "bar"}
        out = rcmod.plotting_context("talk", rc=rc)
        self.assertEqual(out[key], val)
        self.assertNotIn("foo", out)

    def test_set_context(self):

        for context in self.contexts:

            context_dict = rcmod.plotting_context(context)
            rcmod.set_context(context)
            self.assert_rc_params(context_dict)

    def test_context_context_manager(self):

        rcmod.set_context("notebook")
        orig_params = rcmod.plotting_context()
        with rcmod.plotting_context("paper"):
            context_params = rcmod.plotting_context("paper")
            self.assert_rc_params(context_params)
        self.assert_rc_params(orig_params)


class TestFonts(TestCase):

    def test_set_font(self):

        rcmod.set(font="Verdana")

        _, ax = plt.subplots()
        ax.set_xlabel("foo")

        try:
            self.assertEqual(ax.xaxis.label.get_fontname(),
                             "Verdana")
        except AssertionError:
            if has_verdana():
                raise
            else:
                raise unittest.skip("Verdana font is not present")
        finally:
            rcmod.set()
            plt.close("all")

    def test_set_serif_font(self):

        rcmod.set(font="serif")

        _, ax = plt.subplots()
        ax.set_xlabel("foo")

        nt.assert_in(ax.xaxis.label.get_fontname(),
                     mpl.rcParams["font.serif"])

        rcmod.set()
        plt.close("all")

    def test_different_sans_serif(self):

        if LooseVersion(mpl.__version__) < LooseVersion("1.4"):
            raise nose.SkipTest

        rcmod.set()
        rcmod.set_style(rc={"font.sans-serif":
                            ["Verdana"]})

        _, ax = plt.subplots()
        ax.set_xlabel("foo")

        try:
            self.assertEqual(ax.xaxis.label.get_fontname(),
                             "Verdana")
        except AssertionError:
            if has_verdana():
                raise
            else:
                raise unittest.skip("Verdana font is not present")
        finally:
            rcmod.set()
            plt.close("all")


def has_verdana():
    """Helper to verify if Verdana font is present"""
    # This import is relatively lengthy, so to prevent its import for
    # testing other tests in this module not requiring this knowledge,
    # import font_manager here
    import matplotlib.font_manager as mplfm
    try:
        verdana_font = mplfm.findfont('Verdana', fallback_to_default=False)
    except:
        # if https://github.com/matplotlib/matplotlib/pull/3435
        # gets accepted
        return False
    # otherwise check if not matching the logic for a 'default' one
    try:
        unlikely_font = mplfm.findfont("very_unlikely_to_exist1234",
                                       fallback_to_default=False)
    except:
        # if matched verdana but not unlikely, Verdana must exist
        return True
    # otherwise -- if they match, must be the same default
    return verdana_font != unlikely_font
