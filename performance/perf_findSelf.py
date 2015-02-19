from unittest import TestCase

from dark.reads import AARead

from light.database import Database
from light.landmarks import (
    AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, BetaStrand, BetaTurn,
    ALL_LANDMARK_CLASSES)
from light.trig import (
    AminoAcids, Peaks, Troughs, IndividualPeaks, IndividualTroughs,
    ALL_TRIG_CLASSES)


class _FindSelfMixin(object):
    """
    Test if a database built from a sequence finds itself.
    """
    ID = 'gi|460111557|Oropouche virus TRVL-9760, glycoprotein, Segment M'
    SEQUENCE = ('DEDCLSKDIKITYQELHNCIGPKIMGDTCVSKSELYSDLFSKNLVTEYDKKYFEPDTVNDQ'
                'FNKIEFAQDAHRMILLERILYKTECEMLSLKKNSGPYNVAWRTYLKNHNIDLCSRHNYKMI'
                'CQCINAHSMCKNTDIDYNKEIETYYKSNAAAYRADLNTIMDTLKTAFRGLTKVLIENYIEK'
                'DDSDALKALFSNITDSVQDNYQMIGILKFASKLLDINLGRSTRSAHHSIMTNEIPKSNPFT'
                'DYSYSNVNIKECMSPESLKCFKKRDSTPHTNHLLCKIDNKYKAFDWPEIETIQKGQKLCLG'
                'DSHCNLEFTAITADRIMSLTNCYKESFTAQPADMQAGIKKCSADEIGECTTLEDKTWPIVF'
                'CSGKYYYSDSKEHAKDGSINNYCLTNKCSEQRFPIHENWFKKCNWDKTHKEFSTMRQINYN'
                'DITSYRKAIESEIGTDLMTHHYKPTKNLPHVVPRYHSIDVQGTESTEGIINGFIQNTIPAI'
                'SGLGVGYHLNFKGNQLFDIVIFVKKAVYKAQYQKAYTTGPSISINIEHNERCTGHCPEKIP'
                'AKEGWLTFSKEHTSSWGCEEYGCLAIDTGCLYGSCQDVIRPELDVYKKIGSEASLIEICIT'
                'LPHETYCNDMDILEPIIGDKLSASFQNTQTNQLPTLIAYKKGKIYTGQINDIGNTALQCGS'
                'IQVINGSTIGTGSPKFDYICHAMRRKDVIVRKCFNDNYQSCTRLEKRNDLIPYRKGDVIEV'
                'SKTGSNMGQMTFKIELGDINYKIFTKSVDLQMSGICAGCIDCAEGISCSINAEVPAETVCH'
                'CKTNCEDFINNIVFSPQIKNYNIKVHCKSKVEKITAHICGRDIDLQLTVKPYNQKIDLSQL'
                'DESNYIREEDLQCGTWLCKVQKEGIDIMFKGLFSGLGRYWTILIYSIIGVVIIVILVYILL'
                'PIGRLLKAFLIRHEIEYAMEQKIK')

    MAXDISTANCE = None
    MINDISTANCE = None
    LIMITPERLANDMARK = None

    def testFindSelf(self):
        """
        Does a sequence match itself using different landmark and trig point
        finders.
        """
        read = AARead(self.ID, self.SEQUENCE)
        database = Database(self.LANDMARKS, self.TRIG_POINTS,
                            maxDistance=self.MAXDISTANCE,
                            minDistance=self.MINDISTANCE,
                            limitPerLandmark=self.LIMITPERLANDMARK)
        database.addSubject(read)
        result = database.find(read)
        significant = list(result.significant())
        if 0 in significant:
            self.details = {
                'result': True,
                'score': result.analysis[0]['score'],
            }
        else:
            self.details = {
                'result': False,
            }


class TestFindSelfAllFinders(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using all finders.
    """
    LANDMARKS = ALL_LANDMARK_CLASSES
    TRIG_POINTS = ALL_TRIG_CLASSES


class TestFindSelfAHAllTrig(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using AlphaHelix
    landmarks and all trig point finders.
    """
    LANDMARKS = [AlphaHelix]
    TRIG_POINTS = ALL_TRIG_CLASSES


class TestFindSelfAH310AllTrig(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using AlphaHelix_3_10
    landmarks and all trig point finders.
    """
    LANDMARKS = [AlphaHelix_3_10]
    TRIG_POINTS = ALL_TRIG_CLASSES


class TestFindSelfAHpiAllTrig(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using AlphaHelix_pi
    landmarks and all trig point finders.
    """
    LANDMARKS = [AlphaHelix_pi]
    TRIG_POINTS = ALL_TRIG_CLASSES


class TestFindSelfBetaStrandAllTrig(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using BetaStrand
    landmarks and all trig point finders.
    """
    LANDMARKS = [BetaStrand]
    TRIG_POINTS = ALL_TRIG_CLASSES


class TestFindSelfBetaTurnAllTrig(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using BetaTurn
    landmarks and all trig point finders.
    """
    LANDMARKS = [BetaTurn]
    TRIG_POINTS = ALL_TRIG_CLASSES


class TestFindSelfAllLandmarkAA(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using all landmark
    finders and AminoAcid trig points.
    """
    LANDMARKS = ALL_LANDMARK_CLASSES
    TRIG_POINTS = [AminoAcids]


class TestFindSelfAllLandmarkPeaks(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using all landmark
    finders and Peaks trig points.
    """
    LANDMARKS = ALL_LANDMARK_CLASSES
    TRIG_POINTS = [Peaks]


class TestFindSelfAllLandmarkTroughs(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using all landmark
    finders and Troughs trig points.
    """
    LANDMARKS = ALL_LANDMARK_CLASSES
    TRIG_POINTS = [Troughs]


class TestFindSelfAllLandmarkIndividualPeaks(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using all landmark
    finders and IndividualPeaks trig points.
    """
    LANDMARKS = ALL_LANDMARK_CLASSES
    TRIG_POINTS = [IndividualPeaks]


class TestFindSelfAllLandmarkIndividualTroughs(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using all landmark
    finders and IndividualTroughs trig points.
    """
    LANDMARKS = ALL_LANDMARK_CLASSES
    TRIG_POINTS = [IndividualTroughs]


class TestFindSelfAllLandmarksAllTrigLowMinDistance(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using all finders and a
    low minDistance.
    """
    LANDMARKS = ALL_LANDMARK_CLASSES
    TRIG_POINTS = ALL_TRIG_CLASSES
    MINDISTANCE = 2


class TestFindSelfAllLandmarksAllTrigHighMinDistance(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using all finders and a
    high minDistance.
    """
    LANDMARKS = ALL_LANDMARK_CLASSES
    TRIG_POINTS = ALL_TRIG_CLASSES
    MINDISTANCE = 30


class TestFindSelfAllLandmarksAllTrigHighMaxDistance(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using all finders and a
    high maxDistance.
    """
    LANDMARKS = ALL_LANDMARK_CLASSES
    TRIG_POINTS = ALL_TRIG_CLASSES
    MAXDISTANCE = 300


class TestFindSelfAllLandmarksAllTrigLowMaxDistance(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using all finders and a
    low maxDistance.
    """
    LANDMARKS = ALL_LANDMARK_CLASSES
    TRIG_POINTS = ALL_TRIG_CLASSES
    MAXDISTANCE = 10


class TestFindSelfAllLandmarksAllTrigLowLimitPerLandmark(_FindSelfMixin,
                                                         TestCase):
    """
    A class for testing if a sequence matches itself using all finders and a
    low limitPerLandmark.
    """
    LANDMARKS = ALL_LANDMARK_CLASSES
    TRIG_POINTS = ALL_TRIG_CLASSES
    LIMITPERLANDMARK = 10


class TestFindSelfAllLandmarksAllTrigHighLimitPerLandmark(_FindSelfMixin,
                                                          TestCase):
    """
    A class for testing if a sequence matches itself using all finders and a
    high limitPerLandmark.
    """
    LANDMARKS = ALL_LANDMARK_CLASSES
    TRIG_POINTS = ALL_TRIG_CLASSES
    LIMITPERLANDMARK = 100
