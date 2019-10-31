from molvector import encode, decode, canonical_order, mutate, generate_random_smiles_orders

from rdkit.Chem import MolFromSmiles, MolToSmiles
import unittest
import functools

class TestMolVector(unittest.TestCase):
    def test_molvector(self):
        test = "NCCCCCOCC1OC(OCCc2c[nH]c3ccccc23)C(OCc2ccccc2)C(OCc2ccccc2)C1OCc1ccccc1"
        m = MolFromSmiles(test)
        smi = MolToSmiles(m, False)
        vectors = encode(m, functools.partial(generate_random_smiles_orders, N=100))
        for vector in vectors:
            smi2 = MolToSmiles(decode(vector))
            self.assertEqual(smi, smi2)

    def test_mutate(self):
        test3 = "NCCCCCOCC1OC(OCCc2c[nH]c3ccccc23)C(OCc2ccccc2)C(OCc2ccccc2)C1OCc1ccccc1"
        test4 = "NCCCCC(C(=O)NCCc1ccccc1)N1Cc2[nH]c3ccccc3c2CC(NC(=O)Cc2ccccc2)C1=O.O=C(O)C(F)(F)F"
        m = MolFromSmiles(test3)
        m2 = MolFromSmiles(test4)
    
        v = encode(m, canonical_order)[0]
        v2 = encode(m2, canonical_order)[0]
        from rdkit import rdBase
        rdBase.DisableLog("rdApp.*")
        good = 0
        for i in range(10000):
            r = mutate(v,v2)
            try:
                smi = MolToSmiles(decode(r))
                if MolFromSmiles(smi):
                    good += 1
            except:
                raise
        self.assertTrue(good > 100)
        
        

if __name__ == "__main__":
    unittest.main()
