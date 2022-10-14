import sys
sys.path.append('/Users/hoto7260/Hope_LIET/LIET2/liet')
import rnap_lib_data_proc as dp
import unittest



## A class in which we test the funcitons of our module
class ModuleTest(unittest.TestCase):
    
    # a function to test annotation in our module (must start with 'test_')
    
    def test_annotation(self):
        #make example result for string & integer inputs
        res1 = dp.annotation(chrom='chr1', start=2, stop='1', strand='+')
        expres1 = {'chrom': 'chr1', 'start': 2, 'stop': 1, 'strand': 1}
        #make example result for negative strand & variable inputs
        chrom = 'chr2'
        strand = '-'
        res2 = dp.annotation(chrom=chrom, start='2000', stop=2000, strand=strand)
        expres2 = {'chrom': 'chr2', 'start': 2000, 'stop': 2000, 'strand': -1}
        # Check that it produces a dictionary
        self.assertIsInstance(res1, dict, msg="test_annotation not a dict")
        # Check that it produces the correct outputs
        self.assertEqual(res1, expres1)
        self.assertEqual(res2, expres2)
    
    def test_annot_BED6_loader(self):
        #use sample annot bed_file (+, -, small & large #, X & Y chr, weird chr, & out of order for both sorting)
        annot_file_bed6 = "/Users/hoto7260/Hope_LIET/LIET2/liet/tests/Ex_files/samp_annot"
        res1 = dp.annot_BED6_loader(annot_file_bed6, pad5=1, pad3=3)
        expres1 =  ({'chr1': {(20190562, 20191401, 1): 'UBXN10'}, 'chr10': {(246767297, 246, 1): 'SCCPDH'}, 'chr9': {(15517161123412341234, 15719771123412341234, -1): 'SMIM12'}, 'chrWeirdValue': {(32470625, 32485415, 1): 'ZBTB8B'}, 'chrX': {(30713585, 30723551, -1): 'MATN1'}, 'chrY': {(28259848, 28280799, 1): 'SESN2'}}, {'UBXN10': (1, 3), 'SCCPDH': (1, 3), 'SESN2': (1, 3), 'MATN1': (1, 3), 'ZBTB8B':(1, 3), 'SMIM12': (1, 3)})
        
        # Check that it produces two dictionaries in a tuple
        self.assertIsInstance(res1, tuple, msg="test_annot_BED6_loader not right type")
        # Check that it produces the correct values
        self.assertEqual(res1, expres1, msg="test_annot_BED6_loader not right value")
        
#     def test_annot_loader(self):
#         # What happens if does not have chromosomes in human
#         annot_file_gtf = "/Users/hoto7260/Hope_LIET/LIET2/liet/tests/Ex_files/samp_annot_GTF"
#         res1 = dp.annot_loader(annot_file_gtf)
        
#     def test_prior_config(self)
        # priors is dict: {category: pname: pval
        # pname can be: method, optimizer, meanfield: boolean, tolerance: float, samples: float, pad: [float, float], cov_thresholds: [float, float], learning_rate: float, range_shift: float, antisense: float
#         priors = {}
#         res1 = dp.prior_config(priors, tss=1, tcs=2, frac_priors=False)
#         res2 = dp.prior_config(priors, tss=1, tcs=2, frac_priors=True)

    def test_bgline_cast(self):
        # make input & expected output
        bgline = 'chr	23	10000000000000000000	10'
        res1 = dp.bgline_cast(bgline)
        expres1 = ["chr",23, 10000000000000000000, 10]
        # Check that it produces a list
        self.assertIsInstance(res1, list)
        # Check that it produces the correct list
        self.assertEqual(res1, expres1)
        
    def test_bg2d(self):
        res1 = dp.bg2d(["chr",1, 3, 10])
        res2 = dp.bg2d(["chr", 100000234, 100000240, -145143])
        expres1 = {1:10, 2:10}
        expres2 = {100000234:145143, 100000235:145143, 100000236:145143, 100000237:145143, 100000238:145143, 100000239:145143}
        # Check empty dictionary if 0 counts
        self.assertEqual(dp.bg2d(["hello", 1, 2, 0]), {})
        # Check if correct output
        self.assertEqual(res1, expres1)
        # Check it takes a - count value & returns +, & works w/ large #s
        self.assertEqual(res2, expres2)
        
    def test_reads_d2l(self):
        res1 = dp.reads_d2l({1:10, 2:10, 62:4324523})
        expres1 = [10, 20, 268120426]
        # Check result is a list
        self.assertIsInstance(res1, list)
        #Check produces correct results
        self.assertEqual(res1, expres1)
    
#     def test_chrom_order_reader(self):
        
#     def test_bglist_check(self):
        
#     def test_add_bg_dict(self):
        
    def test_bgreads(self):
        
        
    def test_pad_calc(self):
        
    def test_bedgraph_loader(self):
        
    def test_cov_filter(self):
        
#     def test_gene_data(self):
        
#     def test_rng_shift(self):
        
#     def test_overlap_check2(self):
        
#     def test_overlap_check3(self):
        
#     def test_bg_iterator(self):
        
#     def test_bgreads_old(self):

        
        
if __name__ == '__main__':
    unittest.main()
    
