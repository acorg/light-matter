from unittest import TestCase

from dark.reads import AARead

from light.database import Database
from light.landmarks import (
    AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, BetaStrand, BetaTurn)
from light.trig import (
    AminoAcids, Peaks, Troughs, IndividualPeaks, IndividualTroughs)


class _FindSelfMixin(object):
    """
    Test if a database build from a sequence finds itself.
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
    LANDMARKS = [AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, BetaStrand,
                 BetaTurn]
    TRIG_POINTS = [AminoAcids, Peaks, Troughs, IndividualPeaks,
                   IndividualTroughs]
    MAXDISTANCE = None
    MINDISTANCE = None
    LIMITPERLANDMARK = None


class TestFindSelfAHAllTrig(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using AlphaHelix
    landmarks and all trig point finders.
    """
    LANDMARKS = [AlphaHelix]
    TRIG_POINTS = [AminoAcids, Peaks, Troughs, IndividualPeaks,
                   IndividualTroughs]
    MAXDISTANCE = None
    MINDISTANCE = None
    LIMITPERLANDMARK = None


class TestFindSelfAH310AllTrig(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using AlphaHelix_3_10
    landmarks and all trig point finders.
    """
    LANDMARKS = [AlphaHelix_3_10]
    TRIG_POINTS = [AminoAcids, Peaks, Troughs, IndividualPeaks,
                   IndividualTroughs]
    MAXDISTANCE = None
    MINDISTANCE = None
    LIMITPERLANDMARK = None


class TestFindSelfAHpiAllTrig(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using AlphaHelix_pi
    landmarks and all trig point finders.
    """
    LANDMARKS = [AlphaHelix_pi]
    TRIG_POINTS = [AminoAcids, Peaks, Troughs, IndividualPeaks,
                   IndividualTroughs]
    MAXDISTANCE = None
    MINDISTANCE = None
    LIMITPERLANDMARK = None


class TestFindSelfBetaStrandAllTrig(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using BetaStrand
    landmarks and all trig point finders.
    """
    LANDMARKS = [BetaStrand]
    TRIG_POINTS = [AminoAcids, Peaks, Troughs, IndividualPeaks,
                   IndividualTroughs]
    MAXDISTANCE = None
    MINDISTANCE = None
    LIMITPERLANDMARK = None


class TestFindSelfBetaTurnAllTrig(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using BetaTurn
    landmarks and all trig point finders.
    """
    LANDMARKS = [BetaTurn]
    TRIG_POINTS = [AminoAcids, Peaks, Troughs, IndividualPeaks,
                   IndividualTroughs]
    MAXDISTANCE = None
    MINDISTANCE = None
    LIMITPERLANDMARK = None


class TestFindSelfAllLandmarkAA(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using all landmark
    finders and AminoAcid trig points.
    """
    LANDMARKS = [AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, BetaStrand,
                 BetaTurn]
    TRIG_POINTS = [AminoAcids]
    MAXDISTANCE = None
    MINDISTANCE = None
    LIMITPERLANDMARK = None


class TestFindSelfAllLandmarkPeaks(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using all landmark
    finders and Peaks trig points.
    """
    LANDMARKS = [AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, BetaStrand,
                 BetaTurn]
    TRIG_POINTS = [Peaks]
    MAXDISTANCE = None
    MINDISTANCE = None
    LIMITPERLANDMARK = None


class TestFindSelfAllLandmarkTroughs(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using all landmark
    finders and Troughs trig points.
    """
    LANDMARKS = [AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, BetaStrand,
                 BetaTurn]
    TRIG_POINTS = [Troughs]
    MAXDISTANCE = None
    MINDISTANCE = None
    LIMITPERLANDMARK = None


class TestFindSelfAllLandmarkIndividualPeaks(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using all landmark
    finders and IndividualPeaks trig points.
    """
    LANDMARKS = [AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, BetaStrand,
                 BetaTurn]
    TRIG_POINTS = [IndividualPeaks]
    MAXDISTANCE = None
    MINDISTANCE = None
    LIMITPERLANDMARK = None


class TestFindSelfAllLandmarkIndividualTroughs(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using all landmark
    finders and IndividualTroughs trig points.
    """
    LANDMARKS = [AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, BetaStrand,
                 BetaTurn]
    TRIG_POINTS = [IndividualTroughs]
    MAXDISTANCE = None
    MINDISTANCE = None
    LIMITPERLANDMARK = None


class TestFindSelfAllLandmarksAllTrigLowMinDistance(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using all finders and a
    low minDistance.
    """
    LANDMARKS = [AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, BetaStrand,
                 BetaTurn]
    TRIG_POINTS = [AminoAcids, Peaks, Troughs, IndividualPeaks,
                   IndividualTroughs]
    MAXDISTANCE = None
    MINDISTANCE = 2
    LIMITPERLANDMARK = None


class TestFindSelfAllLandmarksAllTrigHighMinDistance(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using all finders and a
    high minDistance.
    """
    LANDMARKS = [AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, BetaStrand,
                 BetaTurn]
    TRIG_POINTS = [AminoAcids, Peaks, Troughs, IndividualPeaks,
                   IndividualTroughs]
    MAXDISTANCE = None
    MINDISTANCE = 30
    LIMITPERLANDMARK = None


class TestFindSelfAllLandmarksAllTrigHighMaxDistance(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using all finders and a
    high maxDistance.
    """
    LANDMARKS = [AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, BetaStrand,
                 BetaTurn]
    TRIG_POINTS = [AminoAcids, Peaks, Troughs, IndividualPeaks,
                   IndividualTroughs]
    MAXDISTANCE = 300
    MINDISTANCE = None
    LIMITPERLANDMARK = None


class TestFindSelfAllLandmarksAllTrigLowMaxDistance(_FindSelfMixin, TestCase):
    """
    A class for testing if a sequence matches itself using all finders and a
    low maxDistance.
    """
    LANDMARKS = [AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, BetaStrand,
                 BetaTurn]
    TRIG_POINTS = [AminoAcids, Peaks, Troughs, IndividualPeaks,
                   IndividualTroughs]
    MAXDISTANCE = 10
    MINDISTANCE = None
    LIMITPERLANDMARK = None


class TestFindSelfAllLandmarksAllTrigLowLimitPerLandmark(_FindSelfMixin,
                                                         TestCase):
    """
    A class for testing if a sequence matches itself using all finders and a
    low limitPerLandmark.
    """
    LANDMARKS = [AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, BetaStrand,
                 BetaTurn]
    TRIG_POINTS = [AminoAcids, Peaks, Troughs, IndividualPeaks,
                   IndividualTroughs]
    MAXDISTANCE = None
    MINDISTANCE = None
    LIMITPERLANDMARK = 10


class TestFindSelfAllLandmarksAllTrigHighLimitPerLandmark(_FindSelfMixin,
                                                          TestCase):
    """
    A class for testing if a sequence matches itself using all finders and a
    high limitPerLandmark.
    """
    LANDMARKS = [AlphaHelix, AlphaHelix_3_10, AlphaHelix_pi, BetaStrand,
                 BetaTurn]
    TRIG_POINTS = [AminoAcids, Peaks, Troughs, IndividualPeaks,
                   IndividualTroughs]
    MAXDISTANCE = None
    MINDISTANCE = None
    LIMITPERLANDMARK = 100
