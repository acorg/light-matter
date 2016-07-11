from pymol.cgo import wire_text, COLOR, SPHERE  # get constants
from pymol.vfont import plain


def makeLegend():
    """
    A function to display a legend for the the features and their respective
    colors.

    @return: A C{list} of objects to be plotted.
    """
    cgo = []

    axes = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]]

    pos = [0.0, 0.0, 0.0]
    wire_text(cgo, plain, pos, 'AC AlphaHelix', axes)

    pos = [0.0, -3.0, 0.0]
    wire_text(cgo, plain, pos, 'AC AlphaHelix_3_10', axes=axes)

    pos = [0.0, -6.0, 0.0]
    wire_text(cgo, plain, pos, 'AC AlphaHelix_pi', axes=axes)

    pos = [0.0, -9.0, 0.0]
    wire_text(cgo, plain, pos, 'AC ExtendedStrand', axes)

    pos = [0.0, -12.0, 0.0]
    wire_text(cgo, plain, pos, 'AlphaHelix', axes=axes)

    pos = [0.0, -15.0, 0.0]
    wire_text(cgo, plain, pos, 'AlphaHelix_3_10', axes=axes)

    pos = [0.0, -18.0, 0.0]
    wire_text(cgo, plain, pos, 'AlphaHelix_pi', axes)

    pos = [0.0, -21.0, 0.0]
    wire_text(cgo, plain, pos, 'AminoAcidsLm', axes=axes)

    pos = [0.0, -24.0, 0.0]
    wire_text(cgo, plain, pos, 'BetaStrand', axes=axes)

    pos = [0.0, -27.0, 0.0]
    wire_text(cgo, plain, pos, 'BetaTurn', axes)

    pos = [0.0, -30.0, 0.0]
    wire_text(cgo, plain, pos, 'ClusterAlphaHelix', axes=axes)

    pos = [0.0, -33.0, 0.0]
    wire_text(cgo, plain, pos, 'GOR4AlphaHelix', axes=axes)

    pos = [0.0, -36.0, 0.0]
    wire_text(cgo, plain, pos, 'GOR4BetaStrand', axes)

    pos = [0.0, -39.0, 0.0]
    wire_text(cgo, plain, pos, 'GOR4Coil', axes=axes)

    pos = [0.0, -42.0, 0.0]
    wire_text(cgo, plain, pos, 'PDB AlphaHelix', axes=axes)

    pos = [0.0, -45.0, 0.0]
    wire_text(cgo, plain, pos, 'PDB AlphaHelix_3_10', axes=axes)

    pos = [0.0, -48.0, 0.0]
    wire_text(cgo, plain, pos, 'PDB AlphaHelix_pi', axes)

    pos = [0.0, -51.0, 0.0]
    wire_text(cgo, plain, pos, 'PDB ExtendedStrand', axes=axes)

    pos = [0.0, -54.0, 0.0]
    wire_text(cgo, plain, pos, 'Prosite', axes=axes)

    pos = [0.0, -57.0, 0.0]
    wire_text(cgo, plain, pos, 'RandomLandmark', axes)

    pos = [0.0, -60.0, 0.0]
    wire_text(cgo, plain, pos, 'THAlphaHelix', axes=axes)

    pos = [0.0, -63.0, 0.0]
    wire_text(cgo, plain, pos, 'AminoAcids', axes=axes)

    pos = [0.0, -66.0, 0.0]
    wire_text(cgo, plain, pos, 'IndividualPeaks', axes)

    pos = [0.0, -69.0, 0.0]
    wire_text(cgo, plain, pos, 'IndividualTroughs', axes=axes)

    pos = [0.0, -72.0, 0.0]
    wire_text(cgo, plain, pos, 'Peaks', axes=axes)

    pos = [0.0, -75.0, 0.0]
    wire_text(cgo, plain, pos, 'RandomTrigPoint', axes)

    pos = [0.0, -78.0, 0.0]
    wire_text(cgo, plain, pos, 'Troughs', axes=axes)

    pos = [0.0, -81.0, 0.0]
    wire_text(cgo, plain, pos, 'Volume', axes=axes)

    cgo.extend([COLOR, 0.5, 1.0, 1.0,
                SPHERE, -1.5, 1.0, 0.0, 1,
                COLOR, 1.0, 0.7, 0.2,
                SPHERE, -1.5, -2.0, 0.0, 1,
                COLOR, 0.73, 0.55, 0.52,
                SPHERE, -1.5, -5.0, 0.0, 1,
                COLOR, 0.6, 0.6, 0.1,
                SPHERE, -1.5, -8.0, 0.0, 1,
                COLOR, 1.0, 0.5, 0.5,
                SPHERE, -1.5, -11.0, 0.0, 1,
                COLOR, 0.1, 0.6, 0.6,
                SPHERE, -1.5, -14.0, 0.0, 1,
                COLOR, 0.698, 0.13, 0.13,
                SPHERE, -1.5, -17.0, 0.0, 1,
                COLOR, 0.25, 1.00, 0.75,
                SPHERE, -1.5, -20.0, 0.0, 1,
                COLOR, 0.75, 0.75, 1.0,
                SPHERE, -1.5, -23.0, 0.0, 1,
                COLOR, 0.75, 1.00, 0.25,
                SPHERE, -1.5, -26.0, 0.0, 1,
                COLOR, 0.0, 0.5, 1.0,
                SPHERE, -1.5, -29.0, 0.0, 1,
                COLOR, 1.0, 0.5, 0.0,
                SPHERE, -1.5, -32.0, 0.0, 1,
                COLOR, 0.65, 0.9, 0.65,
                SPHERE, -1.5, -35.0, 0.0, 1,
                COLOR, 1.0, 0.3, 0.3,
                SPHERE, -1.5, -38.0, 0.0, 1,
                COLOR, 0.2, 1.0, 0.2,
                SPHERE, -1.5, -41.0, 0.0, 1,
                COLOR, 1.0, 0.2, 0.2,
                SPHERE, -1.5, -44.0, 0.0, 1,
                COLOR, 0.55, 0.25, 0.60,
                SPHERE, -1.5, -47.0, 0.0, 1,
                COLOR, 0.85, 0.20, 0.50,
                SPHERE, -1.5, -50.0, 0.0, 1,
                COLOR, 0.619607843, 0.388235294, 0.709803922,
                SPHERE, -1.5, -53.0, 0.0, 1,
                COLOR, 0.341176471, 0.090196078, 0.560784314,
                SPHERE, -1.5, -56.0, 0.0, 1,
                COLOR, 1.000000000, 0.819607843, 0.137254902,
                SPHERE, -1.5, -59.0, 0.0, 1,
                COLOR, 0.2, 0.6, 0.2,
                SPHERE, -1.5, -62.0, 0.0, 1,
                COLOR, 0.5, 0.5, 1.0,
                SPHERE, -1.5, -65.0, 0.0, 1,
                COLOR, 0.3, 0.3, 1.0,
                SPHERE, -1.5, -68.0, 0.0, 1,
                COLOR, 1.0, 0.5, 1.0,
                SPHERE, -1.5, -71.0, 0.0, 1,
                COLOR, 0.819607843, 0.000000000, 0.309803922,
                SPHERE, -1.5, -74.0, 0.0, 1,
                COLOR, 0.090196078, 0.329411765, 0.529411765,
                SPHERE, -1.5, -77.0, 0.0, 1,
                COLOR, 1.000000000, 0.501960784, 0.000000000,
                SPHERE, -1.5, -80.0, 0.0, 1,
                COLOR, 1.0, 1.0, 1.0])

    return cgo
