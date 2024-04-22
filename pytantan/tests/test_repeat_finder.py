import pickle
import unittest
import sys
import textwrap

import pytantan

def inline(x):
    return textwrap.dedent(x).replace("\n", "").strip()

TITIN_H16 = inline(
    """
    MTTQAPTFTQPLQSVVVLEGSTATFEAHISGFPVPEVSWFRDGQVISTST
    LPGVQISFSDGRAKLTIPAVTKANSGRYSLKATNGSGQATSTAELLVKAE
    TAPPNFVQRLQSMTVRQGSQVRLQVRVTGIPTPVVKFYRDGAEIQSSLDF
    QISQEGDLYSLLIAEAYPEDSGTYSVNATNSVGRATSTAELLVQGEEEVP
    AKKTKTIVSTAQISESRQTRIEKKIEAHFDARSIATVEMVIDGAAGQQLP
    HKTPHRIPPKPKSRSPTPPSIAAKAQLARQQSPSPIRHSPSPVRHVRAPT
    PSPVRSVSPAARISTSPIRSVRSPLLMRKTQASTVATGPEVPPPWKQEGY
    VASSSEAEMRETTLTTSTQIRTEERWEGRYGVQEQVTISGAAGAAASVSA
    SASYAAEAVATGAKEVKQDADKSAAVATVVAAVDMARVREPVISAVEQTA
    QRTTTTAVHIQPAQEQVRKEAEKTAVTKVVVAADKAKEQELKSRTKEVIT
    TKQEQMHVTHEQIRKETEKTFVPKVVISAAKAKEQETRISEEITKKQKQV
    TQEAIRQETEITAASMVVVATAKSTKLETVPGAQEETTTQQDQMHLSYEK
    IMKETRKTVVPKVIVATPKVKEQDLVSRGREGITTKREQVQITQEKMRKE
    AEKTALSTIAVATAKAKEQETILRTRETMATRQEQIQVTHGKVDVGKKAE
    AVATVVAAVDQARVREPREPGHLEESYAQQTTLEYGYKERISAAKVAEPP
    QRPASEPHVVPKAVKPRVIQAPSETHIKTTDQKGMHISSQIKKTTDLTTE
    """
)

class TestRepeatFinder(unittest.TestCase):

    def test_hard(self):
        matrix = pytantan.ScoreMatrix.dna()
        
        tantan = pytantan.RepeatFinder(matrix)
        self.assertEqual(
            tantan.mask_repeats("CATCATCATCATCATATCATATCATATCATATCATATCAT"),
            "CATcatcatcatcaTAtcatatcatatcatatcatatcat",
        )

        tantan = pytantan.RepeatFinder(matrix, repeat_start=0.5)
        self.assertEqual(
            tantan.mask_repeats("CATCATCATCATCATATCATATCATATCATATCATATCAT"),
            "CATcatcatcatcatatcatatcatatcatatcatatcat",
        )

    def titin_h16(self):
        matrix = pytantan.ScoreMatrix.protein()

        tantan = pytantan.RepeatFinder(matrix)
        self.assertEqual(
            tantan.mask_repeats(TITIN_H16, mask="x"),
            inline(
                """
                MTTQAPTFTQPLQSVVVLEGSTATFEAHISGFPVPEVSWFRDGQVISTST
                LPGVQISFSDGRAKLTIPAVTKANSGRYSLKATNGSGQATSTAELLVKAE
                TAPPNFVQRLQSMTVRQGSQVRLQVRVTGIPTPVVKFYRDGAEIQSSLDF
                QISQEGDLYSLLIAEAYPEDSGTYSVNATNSVGRATSTAELLVQGEEEVP
                AKKTKTIVSTAQISESRQTRIEKKIEAHFDARSIATVEMVIDGAAGQQLP
                HKTPHRIPPKPKSRSxxxPSIAAKAQLARQQSPSPIxxxxxxxxxxxxxx
                PSPVRSVSPAAxxxxxxxxxxxxxxLMRKTQASTVATGPEVPPPWKQEGY
                VASSSEAEMRETTLTTSTQIRTEERWEGRYGVQExxxxxxxxxxxxxxxx
                xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxDMARVREPVISAVEQTA
                QRTTTTAVHIQxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxISEEITxxxxxx
                xxxxxxxxxxxxxxxxxxxxxxxxxxxxxVPGAQEETTTQQxxxxxxxxx
                xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                xxxxxxxxxxxxxxxxPREPGHLEESYAQQTTLEYGYKERISAAKVAEPP
                QRPASEPHVVPKAVKPRVIQAPSETHIKTTDQKGMHISSQIKKTTDLTTE
                """
            )
        )

        tantan = pytantan.RepeatFinder(matrix)
        self.assertEqual(
            tantan.mask_repeats(TITIN_H16, mask="x"),
            inline(
                """
                MTTQAPTFTQPLQSVVVLEGSTATFEAHISGFPVPEVSWFRDGQVISTST
                LPGVQISFSDGRAKLTIPAVTKANSGRYSLKATNGSGQATSTAELLVKAE
                TAPPNFVQRLQSMTVRQGSQVRLQVRVTGIPTPVVKFYRDGAEIQSSLDF
                QISQEGDLYSLLIAEAYPEDSGTYSVNATNSVGRATSTAELLVQGEEEVP
                AKKTKTIVSTAQISESRQTRIEKKIEAHFDARSIATVEMVIDGAAGQQLP
                HKTPHRIPPKPKSRSxxxPSIAAKAQLARQQSPSPIxxxxxxxxxxxxxx
                PSPVRSVSPAAxxxxxxxxxxxxxxLMRKTQASTVATGPEVPPPWKQEGY
                VASSSEAEMRETTLTTSTQIRTEERWEGRYGVQExxxxxxxxxxxxxxxx
                xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxDMARVREPVISAVEQTA
                QRTTTTAVHIQxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxISEEITxxxxxx
                xxxxxxxxxxxxxxxxxxxxxxxxxxxxxVPGAQEETTTQQxxxxxxxxx
                xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
                xxxxxxxxxxxxxxxxPREPGHLEESYAQQTTLEYGYKERISAAKVAEPP
                QRPASEPHVVPKAVKPRVIQAPSETHIKTTDQKGMHISSQIKKTTDLTTE
                """
            )
        )
        self.assertEqual(
            tantan.mask_repeats(TITIN_H16, mask="X", threshold=0.9),
            inline(
                """
                MTTQAPTFTQPLQSVVVLEGSTATFEAHISGFPVPEVSWFRDGQVISTST
                LPGVQISFSDGRAKLTIPAVTKANSGRYSLKATNGSGQATSTAELLVKAE
                TAPPNFVQRLQSMTVRQGSQVRLQVRVTGIPTPVVKFYRDGAEIQSSLDF
                QISQEGDLYSLLIAEAYPEDSGTYSVNATNSVGRATSTAELLVQGEEEVP
                AKKTKTIVSTAQISESRQTRIEKKIEAHFDARSIATVEMVIDGAAGQQLP
                HKTPHRIPPKPKSRSPTPPSIAAKAQLARQQSPSPIRHXXXXXXXVRAPT
                PSPVRSVSPAARISTSPIRSVRSPLLMRKTQASTVATGPEVPPPWKQEGY
                VASSSEAEMRETTLTTSTQIRTEERWEGRYGVQEQVTISGAAGAXXXXXX
                XXXXXAEAVATGAKEVKQDADKSAAVATVVAAVDMARVREPVISAVEQTA
                QRTTTTAVHIQPAQEQVRKEAEKTAVTKVVVAADKAKEQELKSRTKEVIT
                TKQEQMHXXXXXXXXXXXXXXXXXXXXXXXXXXXXXTRISEEITKKQKXX
                XXXXXXXXXXXXXXXXXXXXXXXXXKLETVPGAQEETTTQQDQMHLSXXX
                XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                XXXXXXXXXXXARVREPREPGHLEESYAQQTTLEYGYKERISAAKVAEPP
                QRPASEPHVVPKAVKPRVIQAPSETHIKTTDQKGMHISSQIKKTTDLTTE
                """
            )
        )