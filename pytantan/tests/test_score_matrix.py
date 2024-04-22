import pickle
import unittest
import sys

import pytantan


class TestScoreMatrix(unittest.TestCase):

    def test_aa_blosum50(self):
        aa = pytantan.ScoreMatrix.protein("BLOSUM50")
        matrix = aa.matrix
        diagonal = [ matrix[i][i] for i in range(len(aa.columns)) ]
        self.assertEqual(len(diagonal), 24)
        self.assertEqual(diagonal, [5, 7, 7, 8, 13, 7, 6, 8, 10, 5, 5, 6, 7, 8, 10, 5, 5, 15, 8, 5, 5, 5 ,-1 ,1])

    def test_aa_blosum62(self):
        aa = pytantan.ScoreMatrix.protein("BLOSUM62")
        matrix = aa.matrix
        diagonal = [ matrix[i][i] for i in range(len(aa.columns)) ]
        self.assertEqual(len(diagonal), 24)
        self.assertEqual(diagonal, [4, 5, 6, 6, 9, 5, 5, 6, 8, 4, 4, 5, 5, 6, 7, 4, 5, 11, 7, 4, 4, 4, -1, 1])

    def test_aa_invalid_name(self):
        with self.assertRaises(ValueError):
            aa = pytantan.ScoreMatrix.protein("nonsensical")

    def test_matrix(self):
        aa = pytantan.ScoreMatrix.protein("BLOSUM50")
        matrix = aa.matrix
        columns = aa.columns
        self.assertEqual(len(columns), 24)
        self.assertEqual(len(matrix), 24)
        for row in matrix:
            self.assertEqual(len(row), 24)

    @unittest.skipUnless(sys.implementation.name == "cpython", "memoryview not supported")
    def test_memoryview(self):
        aa = pytantan.ScoreMatrix.protein("BLOSUM50")
        mem = memoryview(aa)
        self.assertEqual(mem.shape, (64, 64))
        self.assertEqual(mem[0, 0], 5) # A <-> A
        self.assertEqual(mem[3, 3], 6) # E <-> E
        
    def test_init_invalid_length(self):
        with self.assertRaises(ValueError):
            aa = pytantan.ScoreMatrix(
                pytantan.Alphabet("ATGC"),
                [
                    [0, 0, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 0, 0],
                ]
            )
        with self.assertRaises(ValueError):
            aa = pytantan.ScoreMatrix(
                pytantan.Alphabet("ATGC"),
                [
                    [0, 0, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 0, 0],
                    [0, 0, 0],
                ]
            )

    def test_eq(self):
        sm1 = pytantan.ScoreMatrix.protein("BLOSUM50")
        sm2 = pytantan.ScoreMatrix.protein("BLOSUM50")
        sm3 = pytantan.ScoreMatrix.protein("BLOSUM62")
        self.assertEqual(sm1, sm1)
        self.assertEqual(sm1, sm2)
        self.assertNotEqual(sm1, sm3)
        self.assertNotEqual(sm1, 12)

    def test_pickle(self):
        sm1 = pytantan.ScoreMatrix.protein()
        sm2 = pickle.loads(pickle.dumps(sm1))
        self.assertEqual(sm1.alphabet, sm2.alphabet)
        self.assertEqual(sm1.matrix, sm2.matrix)