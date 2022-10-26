import argparse, os
import unittest
from unittest.mock import patch, call, ANY

from bin.get_candidate_ITR_reads_and_IS_contigs import *

class TestCandidateITRReadsAndISContigs(unittest.TestCase):
    TEST_SAM_FILE_IN_1 = 'tests/data/input/test_itr_in_1.sam.mapped.sorted'
    TEST_SAM_FILE_IN_2 = 'tests/data/input/test_itr_in_2.sam.mapped.sorted'
    TEST_FASTA1 = 'tests/data/input/reference_IR_1.fasta'
    TEST_FASTA2 = 'tests/data/input/reference_IR_2.fasta'
    TEST_TAB = 'tests/data/input/test_IR.tab'
    TEST_CONTIG = 'tests/data/input/test_scaffolds.fasta'
    TEST_OUTPUT_PREFIX = 'tests/data/output/test'

    size = get_array_size(TEST_FASTA1, TEST_FASTA2)

    MIN_IS_LEN = 500
    MAX_IS_LEN = 3000

    ir_positions_dict = {
        'NODE_10183_length_410_cov_29.9493': {
            'Seq4_nstart_SRS011144_nend_SRR062032.2555.1_f1': [(100, 126), True],
            'Seq1476536_nstart_SRS011144_nend_SRR062031.31005955.1_f1': [(643, 731), False]
        },
        'NODE_41901_length_110_cov_74.5455': {
            'Seq1477271_nstart_SRS011144_nend_SRR062031.31036872.2_f2': [(7, 83), False]
        },
        'NODE_823_length_1805_cov_1014.02': {
            'Seq134_nstart_SRS011144_nend_SRR062032.12043.1_f1': [(246, 271), True],
            'Seq687973_nstart_SRS011144_nend_SRR062032.29886260.1_f1': [(743, 771), True],
            'Seq1519667_nstart_SRS011144_nend_SRR062031.33018371.2_f2': [(1242, 1266), True]
        }
    }

    def test_get_ir_offset(self):
        actual = get_ir_offset(self.TEST_TAB, 2000000)
        self.assertEqual(actual[0][3], [48, 74])
        self.assertEqual(actual[1][1519666], [30, 56])

    def test_get_ir_positions(self):
        """
        Test contigs with a pair of ITRs
        """
        alloc_arrays = get_ir_offset(self.TEST_TAB, 2000000)
        actual = get_ir_positions([self.TEST_SAM_FILE_IN_1, self.TEST_SAM_FILE_IN_2], alloc_arrays)
        self.maxDiff = None
        self.assertEqual(actual, {
            'NODE_10183_length_410_cov_29.9493': {
                'Seq4_nstart_SRS011144_nend_SRR062032.2555.1_f1_LCoord_30_RCoord_60': [(100, 126), True]
            },
            'NODE_41901_length_110_cov_74.5455': {
                'Seq1477271_nstart_SRS011144_nend_SRR062031.31036872.2_f2_LCoord_30_RCoord_60': [(7, 83), False]
            },
            'NODE_823_length_1805_cov_1014.02': {
                'Seq134_nstart_SRS011144_nend_SRR062032.12043.1_f1_LCoord_30_RCoord_60': [(246, 271), True],
                'Seq687973_nstart_SRS011144_nend_SRR062032.29886260.1_f1_LCoord_30_RCoord_60': [(743, 771), True],
                'Seq1519667_nstart_SRS011144_nend_SRR062031.33018371.2_f2_LCoord_30_RCoord_60': [(251, 277), True]
            },
            'NODE_41829_length_121_cov_87.6667': {
                'Seq1476536_nstart_SRS011144_nend_SRR062031.31005955.1_f1_LCoord_30_RCoord_60': [(43, 131), False]
            }
        })


    def test_split_match_flag(self):
        match_flag = '97M2S'
        actual = split_match_flag(match_flag)

        self.assertEqual(actual, [0,97,2])

        match_flag = '3S80M'
        actual = split_match_flag(match_flag)

        self.assertEqual(actual, [3,80,0])

        match_flag = '3S80M2S'
        actual = split_match_flag(match_flag)

        self.assertEqual(actual, [3,80,2])

        match_flag = '80M'
        actual = split_match_flag(match_flag)

        self.assertEqual(actual, [0,80,0])

        match_flag = '1S25M1I2M1D71M'
        actual = split_match_flag(match_flag)
        self.assertEqual(actual, [1,98,0])

        match_flag = '9S69M1D20M3S'
        actual = split_match_flag(match_flag)
        self.assertEqual(actual, [9,90,3])

        match_flag = '90M1I4M6S'
        actual = split_match_flag(match_flag)
        self.assertEqual(actual, [0,93,6])

    def test_get_contigs_with_itrs(self):

        alloc_arrays = get_ir_offset(self.TEST_TAB, 2000000)
        alloc_arrays, contigs_w_ir = get_contigs_with_itrs(self.ir_positions_dict, alloc_arrays, 2000000, self.MIN_IS_LEN, self.MAX_IS_LEN, 'tests/data/output/test')

        self.assertEqual(alloc_arrays[0][3], 1)
        self.assertEqual(alloc_arrays[1][1476535], 1)
        self.assertEqual(alloc_arrays[0][133], 1)
        self.assertEqual(alloc_arrays[0][134], 0)
        self.assertEqual(alloc_arrays[0][687972], 0)
        self.assertEqual(alloc_arrays[0][1519666], 0)
        self.assertEqual(alloc_arrays[1][1519666], 1)
        self.assertEqual(contigs_w_ir, ['NODE_10183_length_410_cov_29.9493', 'NODE_823_length_1805_cov_1014.02'])

        tab_file = 'tests/data/output/test_contigs_reads_ir_position_info.tab'
        file = open(tab_file, "r")
        actual_tab = "".join(file.readlines())
        #os.remove(tab_file)
        self.assertEqual(actual_tab, """sample_id\tcontig\tread\tIR\tstart_position\tend_position\ntest\tNODE_10183_length_410_cov_29.9493\tSeq4_nstart_SRS011144_nend_SRR062032.2555.1_f1_LCoord_48_RCoord_74\tTrue\t100\t126\ntest\tNODE_823_length_1805_cov_1014.02\tSeq134_nstart_SRS011144_nend_SRR062032.12043.1_f1_LCoord_26_RCoord_51\tTrue\t246\t271\ntest\tNODE_823_length_1805_cov_1014.02\tSeq1519667_nstart_SRS011144_nend_SRR062031.33018371.2_f2_LCoord_30_RCoord_56\tTrue\t1242\t1266\n""")

    def test_write_fasta_files(self):
        """
        Test fasta file outputs
        """
        position_arrays = get_ir_offset(self.TEST_TAB, 2000000)
        ir_arrays, contigs_w_ir = get_contigs_with_itrs(self.ir_positions_dict, position_arrays, 2000000, self.MIN_IS_LEN, self.MAX_IS_LEN, 'tests/data/output/test')
        write_fasta_files(self.TEST_FASTA1, self.TEST_FASTA2, ir_arrays, position_arrays, self.TEST_OUTPUT_PREFIX)

        f1_name = self.TEST_OUTPUT_PREFIX + '_reads_with_candidate_itrs_1.fasta'
        f1 = open(f1_name, "r")
        actual = "".join(f1.readlines())
        os.remove(f1_name)
        self.assertEqual(actual, """>Seq4_nstart_SRS011144_nend_SRR062032.2555.1_f1_LCoord_48_RCoord_74\nACTGTACATAAAATATCTAACTACCAAAACTATTTATTTTATACCCGACATAAAATATCAAAGTACCCAAACTATATATAGTATACTGTACATCAAACA\n>Seq134_nstart_SRS011144_nend_SRR062032.12043.1_f1_LCoord_26_RCoord_51\nTTATATATCATAATACATATAAATGTGTATTATGTTATCTATAATTATATAATTTCATATATAAGATGTATAATATGTATTAT\n""")

        f2_name = self.TEST_OUTPUT_PREFIX + '_reads_with_candidate_itrs_2.fasta'
        f2 = open(f2_name, "r")
        actual = "".join(f2.readlines())
        os.remove(f2_name)
        self.assertEqual(actual, """>Seq1476536_nstart_SRS011144_nend_SRR062031.31005955.2_f2\nAAAGAAGAAGTCGCTAAGGGAACAGATCCAAAAAATCCAAACTCTAAACCAGAATCTAAACCAGAAGATCCATCAACTAAGGATACC\n>Seq1519667_nstart_SRS011144_nend_SRR062031.33018371.2_f2_LCoord_30_RCoord_56\nCATTATACATCTTATATATGAAATTTTACAATTATAGATAACATAATACACATTTGTATGTATTATGATATAAAACAAAGGTGCTA\n""")

    def test_write_contig_file(self):
        """
        Test contig file output
        """
        contigs = ['NODE_10183_length_410_cov_29.9493', 'NODE_823_length_1805_cov_1014.02']
        write_contig_file(contigs, self.TEST_CONTIG, self.TEST_OUTPUT_PREFIX)

        contig_file_name = self.TEST_OUTPUT_PREFIX + '_contigs_with_candidate_itrs.fa'
        file = open(contig_file_name, "r")
        actual = "".join(file.readlines())
        os.remove(contig_file_name)
        self.assertEqual(actual, """>NODE_823_length_1805_cov_1014.02\nGTTGCAGTGAGCAAAGATGGCGCCATTGTACCCCACCCTGTGCAACAGGACAAGACTCTG\nTCTAAAAAAAAATTATATATGATATATATTACATGTTATGTGCTATGCCTTATATGTAAC\nATATAACCTCATATATTTTATATGTCATAGTATATAATATTTATCCTGTGTCATCTTATA\nTATACTATATATATTACATAACATATATAATATATGATACATATTATACATCTTATATAT\nGAAATTATATAATTATAGATAACATAATACACATTTATATGTATTATGATATATAACAAA\nGGTACTATATCATATATGATATATAATATATATATTTGTAGTACATATTATAATTTTTAT\nATATTATAATTTATATAATGATAAAATTTTTTATAGAATATAAATAATTATATATGATTA\nTATAATTCTACATATTATAAATGAAAATATGTAGAATTACATGTTTATCTAGTTATATGT\nATAAAAATATAATTATATATATAGTTATATATATCACATAATAGACAACATACATATTAT\nATAATGTGTAATATATATTATAGAATATATGATACATGATATATATGATATTTAATATAA\nCATAATATAATCTATATTATAATATTATGTATGTTACTTATATTGGGTGATATGTAATAT\nATATTATGTAATATGAAATAATATAATATATATTATATTATGATATTTTATGTAAATTAT\nGTTATAAAAGTATATATAACATAATATATAGTTATATATAATATATTATATAGTTATATA\nTATTATATTACATATTTACATAAAATATTATACAGTTATATATAATTTCAAATAGTTATA\nTATAACAGAATATAAAACATATAGAATACATATCATAAAATATATATTGTATACCATATA\nTATTAGGTATCATATATTGTATGTTATATATCATATACTGTATATCATATATCATATATT\nATATGTTATATATATCACATTCCACATATGGTATATTATAAATTATACATAATATATTGC\nATTCTTTATTTTACATGTAATAAATTATACATTATAGGTAACATAGTATATATTGTATGT\nAACCGGTATATTCTATGTAAAATATATAATATATAACACATGTGACATGTAATTAATAAC\nATAGTCATATAATATATATTTAATATAGGTATTATACATACAACGTATTATATATAATAT\nATAATATATATTACATATGTTATAAAATATAATATATATTATATGTCGCATTTTATATAT\nACAATTATATATTATTATATGTTGTTTTATACTATATATAACATAGACATAATATATATA\nTTATGCATAATGGGTATTATATATAACAGATATATTATATATAATATATATAATGTATAA\nTTTATATTTAATACGTATATTAACCGTATAATATATTCTATATATGATATATAACATATA\nTTATATATGTCATAAATTATAATATATACTATATGTCATATTGTATAATTATAAAACTTA\nTATATTATCAATTATTATATTTTATACAATATATATAACATAAGATGCAGGATGTAAAAC\nAATATTATACATAGGTTATATATATTATATATATGATACTGTATATTATTATATACATAA\nTGTATACGTATGTTTGTGGGTGCCCTTTTTCCCATCTCATAACTTGTTTTAAGAAGCGCA\nGCCTAATAACGTGTGGGCTTGGGATTCAGTTCTTGAAACAAAACTCTGAGCCTTCAATGA\nCCTTTCGGTCTATGTAAAAGCACTCCTGTCTTCCTGGCAGCAGTTGGACCTCACAATGTG\nGATAG\n>NODE_10183_length_410_cov_29.9493\nCCAAACATTTATAATAAACTGTACATAAAATATCAAAGTACTCAAACTATATACTGTACA\nTAAAATATCTAAGTACCAAAACTATTTATTTTATACCCGACATAAAATATCAAAGTACCC\nAAACTATATATAGTATACTGTACATCAAATATCGAAGTACCCAAAGTATGTATTATATAC\nTGTACATAAAATATCCAACTACCGAAAGTATGTGTTATATACTGCACATAAAATATCAAA\nGTACCTAAACTATATATTATATACTGTATTTAAAATAACAAAGTATCCAAGGAATGTATT\nATGTACTGTACATAAAATGTCAAAGTACCCTAAGTATGTATTATCTATTGTACATAAAAC\nATCAAAGTACCCAAAGTATGCATTATATACTGTACATAAAATATCAAAGT\n""")


    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--contig_fasta', 'contig_path', '--sam_file1', 'sam_path1', '--sam_file2', 'sam_path2', '--fasta1', 'fasta1_path',
            '--fasta2', 'fasta2_path', '--tab_file', 'tab_file', '--min_is_len', '500', '--max_is_len', '3000', '--output_prefix', 'output_prefix'])
        self.assertEqual(actual,
                         argparse.Namespace(contig_fasta='contig_path', sam_file1='sam_path1', sam_file2='sam_path2',
                         fasta_file1='fasta1_path', fasta_file2='fasta2_path', tab_file='tab_file', min_is_len=500,
                         max_is_len=3000, output_prefix='output_prefix'))

    def test_arguments_short_options(self):
        actual = get_arguments().parse_args(
            ['-c', 'contig_path', '-s1', 'sam_path1', '-s2', 'sam_path2', '-f1', 'fasta1_path',
            '-f2', 'fasta2_path', '-t', 'tab_file', '-min', '500', '-max', '3000', '-o', 'output_prefix'])
        self.assertEqual(actual,
                         argparse.Namespace(contig_fasta='contig_path', sam_file1='sam_path1', sam_file2='sam_path2',
                         fasta_file1='fasta1_path', fasta_file2='fasta2_path', tab_file='tab_file', min_is_len=500,
                         max_is_len=3000, output_prefix='output_prefix'))

    @patch('bin.get_candidate_ITR_reads_and_IS_contigs.get_array_size')
    @patch('bin.get_candidate_ITR_reads_and_IS_contigs.get_ir_offset')
    @patch('bin.get_candidate_ITR_reads_and_IS_contigs.get_ir_positions')
    @patch('bin.get_candidate_ITR_reads_and_IS_contigs.get_contigs_with_itrs')
    @patch('bin.get_candidate_ITR_reads_and_IS_contigs.write_contig_file')
    @patch('bin.get_candidate_ITR_reads_and_IS_contigs.write_fasta_files')
    def test_main(self, mock_write_fasta_files, mock_write_contig_file, mock_get_contigs_with_itrs, mock_get_ir_positions, mock_get_ir_offset, mock_get_array_size):
        args = get_arguments().parse_args(
            ['--contig_fasta', 'contig_path', '--sam_file1', 'sam_path1', '--sam_file2', 'sam_path2', '--fasta1', 'fasta1_path',
            '--fasta2', 'fasta2_path', '--tab_file', 'tab_file', '--min_is_len', '500', '--max_is_len', '3000', '--output_prefix', 'output_prefix'])
        mock_get_array_size.return_value = 2000000
        mock_get_ir_offset.return_value = ['pos_array1', 'pos_array2']
        mock_get_contigs_with_itrs.return_value = ['read_array1', 'read_array2'], ['contig1', 'contig2']
        mock_get_ir_positions.return_value = {'contig1': {'read1': ['positions']}}

        main(args)

        mock_get_array_size.call_args_list = call('fasta1_path', 'fasta2_path')
        mock_get_ir_offset.call_args_list = call('tab_file', 2000000)
        mock_get_ir_positions.call_args_list = call([args.sam_file1, args.sam_file2], ['pos_array1', 'pos_array2'])
        mock_get_contigs_with_itrs.call_args_list = call({'contig1': {'read1': ['positions']}}, ['pos_array1', 'pos_array2'], 2000000, 500, 3000, args.output_prefix)
        mock_write_contig_file.call_args_list = call(['contig1', 'contig2'], 'contig_path', 'output_prefix')
        mock_write_fasta_files.call_args_list = call('fasta1_path', 'fasta2_path', ['read_array1', 'read_array2'], ['pos_array1', 'pos_array2'], 'output_prefix')
