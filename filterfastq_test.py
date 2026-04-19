import pytest
import os
from main import quality_filter
from main import filter_fastq
from main import is_dna
from main import gc_filter
from main import complement, reverse_complement, reverse
from main import file_to_dict
from main import gc_perc


class InputFormatTests:

    def test_empty_quality(self):
        with pytest.raises(ZeroDivisionError):
            quality_filter('', 30)


    def test_is_dna_works_with_nucleotide_sequences(self):
        assert is_dna([1,2,3,4,5]) is False


    def test_edge_case_for_gc_filter(self):
        assert gc_filter('GGGGGGGGGGG', (0, 10)) is False
        assert gc_filter('AAAAAAAAAAA', (100, 100)) is False
        assert gc_filter('C', (100, 100)) is True
        assert gc_filter('', (0, 0)) is True


class TestsForFileManipulations:

    def test_no_such_file_or_directory(self):
        with pytest.raises(FileNotFoundError):
            filter_fastq('epsteinfile.fastq', 100, 100, 0)

    def test_is_file_created(self, tmp_path):
        input = tmp_path/'input.fastq'

        input.write_text("""
                            @read1
                            ATCGATCGATCG
                            +
                            IIIIIIIIIIII
                            """)
        
        result = filter_fastq(
            str(input),
            gc_bounds=(0,100),
            length_bounds=(0,2**32),
            quality_threshold=0
        )

        assert os.path.exists(result)

    def test_file_to_dict_return_type(self, tmp_path):
        file = tmp_path/'test.fastq'
        file.write_text("""
        @read1
        ATCGATCGATCG
        +
        IIIIIIIIIIII
        """)

        result = file_to_dict(str(file))

        assert isinstance(result, dict)

        for name, (seq, quality) in result.items():
            assert isinstance(name, str)
            assert isinstance(seq, str)
            assert isinstance(quality, str)


class IntegrationTests:

    def test_integration_of_complement_and_reverse_functions(self):
        seq = 'ATGCATGCATGC'

        expected = reverse(complement(seq))
        result = reverse_complement(seq)

        assert result == expected

class PerformanceTests:

    def test_overload(self):
        huge_sequence = "AAATTTGGGCCC" * 10000000
        assert gc_perc(huge_sequence) < 100
