from unittest import TestCase

import matplotlib.pyplot as plt

from light.colors import miscplot as misc
from light.colors.colors import color_palette


class TestPalPlot(TestCase):
    """
    Test the function that visualizes a color palette.
    """
    def test_palplot_size(self):
        """
        The plot that is made must have the right dimensions.
        """
        pal4 = color_palette("husl", 4)
        misc.palplot(pal4)
        size4 = plt.gcf().get_size_inches()
        self.assertEqual(tuple(size4), (4, 1))

        pal5 = color_palette("husl", 5)
        misc.palplot(pal5)
        size5 = plt.gcf().get_size_inches()
        self.assertEqual(tuple(size5), (5, 1))

        palbig = color_palette("husl", 3)
        misc.palplot(palbig, 2)
        sizebig = plt.gcf().get_size_inches()
        self.assertEqual(tuple(sizebig), (6, 2))

        plt.close("all")
