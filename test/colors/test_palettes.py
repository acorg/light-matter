from unittest import TestCase

import colorsys
import numpy as np
import matplotlib as mpl

from light.colors import palettes, utils, rcmod, husl
from light.colors.xkcd_rgb import xkcd_rgb


class TestColorPalettes(TestCase):

    def test_current_palette(self):

        pal = palettes.color_palette(["red", "blue", "green"], 3)
        rcmod.set_palette(pal, 3)
        self.assertEqual(pal, mpl.rcParams["axes.color_cycle"])
        rcmod.set()

    def test_palette_context(self):

        default_pal = palettes.color_palette()
        context_pal = palettes.color_palette("muted")

        with palettes.color_palette(context_pal):
            self.assertEqual(mpl.rcParams["axes.color_cycle"], context_pal)

        self.assertEqual(mpl.rcParams["axes.color_cycle"], default_pal)

    def test_big_palette_context(self):

        default_pal = palettes.color_palette()
        context_pal = palettes.color_palette("husl", 10)

        with palettes.color_palette(context_pal, 10):
            self.assertEqual(mpl.rcParams["axes.color_cycle"], context_pal)

        self.assertEqual(mpl.rcParams["axes.color_cycle"], default_pal)

    def test_seaborn_palettes(self):

        pals = "deep", "muted", "pastel", "bright", "dark", "colorblind"
        for name in pals:
            pal_out = palettes.color_palette(name)
            self.assertEqual(len(pal_out), 6)

    def test_hls_palette(self):

        hls_pal1 = palettes.hls_palette()
        hls_pal2 = palettes.color_palette("hls")
        self.assertEqual(hls_pal1, hls_pal2)

    def test_mpl_palette(self):

        mpl_pal1 = palettes.mpl_palette("Reds")
        mpl_pal2 = palettes.color_palette("Reds")
        self.assertEqual(mpl_pal1, mpl_pal2)

    def test_mpl_dark_palette(self):

        mpl_pal1 = palettes.mpl_palette("Blues_d")
        mpl_pal2 = palettes.color_palette("Blues_d")
        self.assertEqual(mpl_pal1, mpl_pal2)

    def test_bad_palette_name(self):

        with self.assertRaises(ValueError):
            palettes.color_palette("IAmNotAPalette")

    def test_terrible_palette_name(self):

        with self.assertRaises(ValueError):
            palettes.color_palette("jet")

    def test_bad_palette_colors(self):

        pal = ["red", "blue", "iamnotacolor"]
        with self.assertRaises(ValueError):
            palettes.color_palette(pal)

    def test_palette_desat(self):

        pal1 = palettes.husl_palette(6)
        pal1 = [utils.desaturate(c, .5) for c in pal1]
        pal2 = palettes.color_palette("husl", desat=.5)
        self.assertEqual(pal1, pal2)

    def test_palette_is_list_of_tuples(self):

        pal_in = np.array(["red", "blue", "green"])
        pal_out = palettes.color_palette(pal_in, 3)

        self.assertIsInstance(pal_out, list)
        self.assertIsInstance(pal_out[0], tuple)
        self.assertIsInstance(pal_out[0][0], float)
        self.assertEqual(len(pal_out[0]), 3)

    def test_palette_cycles(self):

        deep = palettes.color_palette("deep")
        double_deep = palettes.color_palette("deep", 12)
        self.assertEqual(double_deep, deep + deep)

    def test_cbrewer_qual(self):

        pal_short = palettes.mpl_palette("Set1", 4)
        pal_long = palettes.mpl_palette("Set1", 6)
        self.assertEqual(pal_short, pal_long[:4])

        pal_full = palettes.mpl_palette("Set2", 8)
        pal_long = palettes.mpl_palette("Set2", 10)
        self.assertEqual(pal_full, pal_long[:8])

    def test_mpl_reversal(self):

        pal_forward = palettes.mpl_palette("BuPu", 6)
        pal_reverse = palettes.mpl_palette("BuPu_r", 6)
        self.assertEqual(pal_forward, pal_reverse[::-1])

    def test_rgb_from_hls(self):

        color = .5, .8, .4
        rgb_got = palettes._color_to_rgb(color, "hls")
        rgb_want = colorsys.hls_to_rgb(*color)
        self.assertEqual(rgb_got, rgb_want)

    def test_rgb_from_husl(self):

        color = 120, 50, 40
        rgb_got = palettes._color_to_rgb(color, "husl")
        rgb_want = husl.husl_to_rgb(*color)
        self.assertEqual(rgb_got, rgb_want)

    def test_rgb_from_xkcd(self):

        color = "dull red"
        rgb_got = palettes._color_to_rgb(color, "xkcd")
        rgb_want = xkcd_rgb[color]
        self.assertEqual(rgb_got, rgb_want)

    def test_light_palette(self):

        pal_forward = palettes.light_palette("red")

        red = tuple(mpl.colors.colorConverter.to_rgba("red"))
        self.assertEqual(tuple(pal_forward[-1]), red)

        pal_cmap = palettes.light_palette("blue", as_cmap=True)
        self.assertIsInstance(pal_cmap, mpl.colors.LinearSegmentedColormap)

    def test_dark_palette(self):

        pal_forward = palettes.dark_palette("red")

        red = tuple(mpl.colors.colorConverter.to_rgba("red"))
        self.assertEqual(tuple(pal_forward[-1]), red)

        pal_cmap = palettes.dark_palette("blue", as_cmap=True)
        self.assertIsInstance(pal_cmap, mpl.colors.LinearSegmentedColormap)

    def test_blend_palette(self):

        colors = ["red", "yellow", "white"]
        pal_cmap = palettes.blend_palette(colors, as_cmap=True)
        self.assertIsInstance(pal_cmap, mpl.colors.LinearSegmentedColormap)

    def test_cubehelix_against_matplotlib(self):

        x = np.linspace(0, 1, 8)
        mpl_pal = mpl.cm.cubehelix(x)[:, :3].tolist()

        sns_pal = palettes.cubehelix_palette(8, start=0.5, rot=-1.5, hue=1,
                                             dark=0, light=1, reverse=True)
        self.assertEqual(sns_pal, mpl_pal)

    def test_cubehelix_n_colors(self):

        for n in [3, 5, 8]:
            pal = palettes.cubehelix_palette(n)
            self.assertEqual(len(pal), n)

    def test_cubehelix_reverse(self):

        pal_forward = palettes.cubehelix_palette()
        pal_reverse = palettes.cubehelix_palette(reverse=True)
        self.assertEqual(pal_forward, pal_reverse[::-1])

    def test_cubehelix_cmap(self):

        cmap = palettes.cubehelix_palette(as_cmap=True)
        self.assertIsInstance(cmap, mpl.colors.ListedColormap)

        x = np.linspace(0, 1, 6)
        cmap_rev = palettes.cubehelix_palette(as_cmap=True, reverse=True)
        x = np.linspace(0, 1, 6)
        pal_forward = cmap(x).tolist()
        pal_reverse = cmap_rev(x[::-1]).tolist()
        self.assertEqual(pal_forward, pal_reverse)

    def test_xkcd_palette(self):

        names = list(xkcd_rgb.keys())[10:15]
        colors = palettes.xkcd_palette(names)
        for name, color in zip(names, colors):
            as_hex = mpl.colors.rgb2hex(color)
            self.assertEqual(as_hex, xkcd_rgb[name])
