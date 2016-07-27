from pymol.cgo import wire_text, COLOR, SPHERE  # get constants
from pymol.vfont import plain


def makeLegend(features):
    """
    A function to display a legend for the features and their respective
    colors.

    @param features: A C{set} of feature names.
    @return: A C{list} of objects to be plotted.
    """
    cgo = []

    colors = [[0.5, 1.0, 1.0], [1.0, 0.7, 0.2], [0.73, 0.55, 0.52],
              [1.0, 0.5, 0.5], [0.6, 0.6, 0.1], [0.1, 0.6, 0.6],
              [0.698, 0.13, 0.13], [0.25, 1.00, 0.75], [0.75, 0.75, 1.0],
              [0.75, 1.00, 0.25], [0.0, 0.5, 1.0], [1.0, 0.5, 0.0],
              [0.65, 0.9, 0.65], [1.0, 0.3, 0.3], [0.2, 1.0, 0.2],
              [1.0, 0.2, 0.2], [0.55, 0.25, 0.60], [0.85, 0.20, 0.50],
              [0.619607843, 0.388235294, 0.709803922],
              [0.341176471, 0.090196078, 0.560784314],
              [1.0, 0.819607843, 0.137254902], [0.2, 0.6, 0.2],
              [0.5, 0.5, 1.0], [0.3, 0.3, 1.0], [1.0, 0.5, 1.0],
              [0.819607843, 0.0, 0.309803922],
              [0.090196078, 0.329411765, 0.529411765],
              [1.0, 0.501960784, 0.0],
              [1.0, 0.0, 1.0], [1.0, 1.0, 0.0],
              [1.0, 1.0, 1.0]]

    axes = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]
    textY = 0.0
    sphereY = 1.0
    spheres = []
    for i, featureName in enumerate(features):
        wire_text(cgo, plain, [0.0, textY, 0.0], featureName[1], axes=axes)
        textY -= 3.0
        spheres.extend([COLOR] + colors[i] +
                       [SPHERE, -1.5, sphereY, 0.0, 1])
        sphereY -= 3.0
    cgo.extend(spheres)

    return cgo
