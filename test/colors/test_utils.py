import warnings

from unittest import TestCase

import numpy as np
import matplotlib.pyplot as plt
from light.colors import utils, rcmod


a_norm = np.random.randn(100)


class TestPlottingUtilities(TestCase):
    """
    Tests for plotting utilities.
    """
    def test_pmf_hist_basics(self):
        """
        Test the function to return barplot args for pmf hist.
        """
        out = utils.pmf_hist(a_norm)
        self.assertEqual(len(out), 3)
        x, h, w = out
        self.assertEqual(len(x), len(h))

        # Test simple case
        a = np.arange(10)
        x, h, w = utils.pmf_hist(a, 10)
        self.assertTrue(np.all(h == h[0]))

    def test_pmf_hist_widths(self):
        """
        Test histogram width is correct.
        """
        x, h, w = utils.pmf_hist(a_norm)
        self.assertEqual(x[1] - x[0], w)

    def test_pmf_hist_normalization(self):
        """
        Test that output data behaves like a PMF.
        """
        x, h, w = utils.pmf_hist(a_norm)
        self.assertAlmostEqual(sum(h), 1)
        # self.assertAlmostEqual(h.max(), 1)

    def test_pmf_hist_bins(self):
        """
        Test bin specification.
        """
        x, h, w = utils.pmf_hist(a_norm, 20)
        self.assertEqual(len(x), 20)

    def test_desaturate(self):
        """
        Test color desaturation.
        """
        out1 = utils.desaturate("red", .5)
        self.assertEqual(out1, (.75, .25, .25))

        out2 = utils.desaturate("#00FF00", .5)
        self.assertEqual(out2, (.25, .75, .25))

        out3 = utils.desaturate((0, 0, 1), .5)
        self.assertEqual(out3, (.25, .25, .75))

        out4 = utils.desaturate("red", .5)
        self.assertEqual(out4, (.75, .25, .25))

    def test_saturate(self):
        """
        Test performance of saturation function.
        """
        out = utils.saturate((.75, .25, .25))
        self.assertEqual(out, (1, 0, 0))

    def test_iqr(self):
        """
        Test the IQR function.
        """
        a = np.arange(5)
        iqr = utils.iqr(a)
        self.assertEqual(iqr, 2)


class TestSpineUtils(TestCase):

    sides = ["left", "right", "bottom", "top"]
    outer_sides = ["top", "right"]
    inner_sides = ["left", "bottom"]

    offset = 10
    original_position = ("outward", 0)
    offset_position = ("outward", offset)

    def test_despine(self):
        f, ax = plt.subplots()
        for side in self.sides:
            self.assertTrue(ax.spines[side].get_visible())

        utils.despine()
        for side in self.outer_sides:
            self.assertTrue(~ax.spines[side].get_visible())
        for side in self.inner_sides:
            self.assertTrue(ax.spines[side].get_visible())

        utils.despine(**dict(zip(self.sides, [True] * 4)))
        for side in self.sides:
            self.assertTrue(~ax.spines[side].get_visible())

        plt.close("all")

    def test_despine_specific_axes(self):
        f, (ax1, ax2) = plt.subplots(2, 1)

        utils.despine(ax=ax2)

        for side in self.sides:
            self.assertTrue(ax1.spines[side].get_visible())

        for side in self.outer_sides:
            self.assertTrue(~ax2.spines[side].get_visible())
        for side in self.inner_sides:
            self.assertTrue(ax2.spines[side].get_visible())

        plt.close("all")

    def test_despine_with_offset(self):
        f, ax = plt.subplots()

        for side in self.sides:
            self.assertEqual(ax.spines[side].get_position(),
                             self.original_position)

        utils.despine(ax=ax, offset=self.offset)

        for side in self.sides:
            is_visible = ax.spines[side].get_visible()
            new_position = ax.spines[side].get_position()
            if is_visible:
                self.assertEqual(new_position, self.offset_position)
            else:
                self.assertEqual(new_position, self.original_position)

        plt.close("all")

    def test_despine_with_offset_specific_axes(self):
        f, (ax1, ax2) = plt.subplots(2, 1)

        utils.despine(offset=self.offset, ax=ax2)

        for side in self.sides:
            self.assertEqual(ax1.spines[side].get_position(),
                             self.original_position)
            if ax2.spines[side].get_visible():
                self.assertEqual(ax2.spines[side].get_position(),
                                 self.offset_position)
            else:
                self.assertEqual(ax2.spines[side].get_position(),
                                 self.original_position)
        plt.close("all")

    def test_despine_trim_spines(self):
        f, ax = plt.subplots()
        ax.plot([1, 2, 3], [1, 2, 3])
        ax.set_xlim(.75, 3.25)

        utils.despine(trim=True)
        for side in self.inner_sides:
            bounds = ax.spines[side].get_bounds()
            self.assertEqual(bounds, (1, 3))

        plt.close("all")

    def test_offset_spines_warns(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always", category=UserWarning)

            f, ax = plt.subplots()
            utils.offset_spines(offset=self.offset)
            self.assertTrue('deprecated' in str(w[0].message))
            self.assertTrue(issubclass(w[0].category, UserWarning))

        plt.close('all')

    def test_offset_spines(self):
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always", category=UserWarning)
            f, ax = plt.subplots()

            for side in self.sides:
                self.assertEqual(ax.spines[side].get_position(),
                                 self.original_position)

            utils.offset_spines(offset=self.offset)

            for side in self.sides:
                self.assertEqual(ax.spines[side].get_position(),
                                 self.offset_position)

        plt.close("all")

    def test_offset_spines_specific_axes(self):
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always", category=UserWarning)
            f, (ax1, ax2) = plt.subplots(2, 1)

            utils.offset_spines(offset=self.offset, ax=ax2)

            for side in self.sides:
                self.assertEqual(ax1.spines[side].get_position(),
                                 self.original_position)
                self.assertEqual(ax2.spines[side].get_position(),
                                 self.offset_position)
        plt.close("all")


def test_ticklabels_overlap():

    rcmod.set()
    f, ax = plt.subplots(figsize=(2, 2))
    f.tight_layout()  # This gets the Agg renderer working

    assert not utils.axis_ticklabels_overlap(ax.get_xticklabels())

    big_strings = "abcdefgh", "ijklmnop"
    ax.set_xlim(-.5, 1.5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(big_strings)

    assert utils.axis_ticklabels_overlap(ax.get_xticklabels())

    x, y = utils.axes_ticklabels_overlap(ax)
    assert x
    assert not y
