# Copyright (c) 2021 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
import argparse
import contextlib
import functools
import logging
from typing import Callable, Iterable, Iterator, List, Tuple

import dnaio

import xopen  # type: ignore

from ._filters import (
    AverageErrorRateFilter,
    DEFAULT_PHRED_SCORE_OFFSET,
    MaximumLengthFilter,
    MedianQualityFilter,
    MinimumLengthFilter,
    average_error_rate,
    qualmean,
    qualmedian,
)

__version__ = "1.0.0-dev"

__all__ = [
    "file_to_fastq_records",
    "fastq_records_to_file",
    "filter_fastq",
    "AverageErrorRateFilter",
    "MaximumLengthFilter",
    "MedianQualityFilter",
    "MinimumLengthFilter",
    "average_error_rate",
    "qualmean",
    "qualmedian",
    "DEFAULT_PHRED_SCORE_OFFSET"
]

DEFAULT_COMPRESSION_LEVEL = 2


def file_to_fastq_records(filepath: str) -> Iterator[dnaio.Sequence]:
    """Parse a FASTQ file into a generator of Sequence objects"""
    opener = functools.partial(xopen.xopen, threads=0)
    with dnaio.open(filepath, opener=opener) as record_h:  # type: ignore
        yield from record_h


def fastq_records_to_file(records: Iterable[dnaio.Sequence], filepath: str,
                          compression_level: int = DEFAULT_COMPRESSION_LEVEL):
    with xopen.xopen(filepath, mode='wb', threads=0,
                     compresslevel=compression_level) as output_h:
        for record in records:
            output_h.write(record.fastq_bytes())


def multiple_files_to_records(input_files: List[str], threads: int = 0
                              ) -> Iterator[Tuple[dnaio.SequenceRecord, ...]]:
    readers = [file_to_fastq_records(f) for f in input_files]
    iterators = [iter(reader) for reader in readers]

    # By differentiating between single, paired and multiple files we can
    # choose the fastest method for each situation.
    if len(iterators) == 1:
        yield from zip(iterators[0])  # This yields tuples of length 1.
    elif len(iterators) == 2:
        for record1, record2 in zip(*iterators):
            if not record1.is_mate(record2):
                raise dnaio.FastqFormatError(
                    f"Records are out of sync, names "
                    f"{record1}, f{record2} do not match.",
                    line=None)
            yield record1, record2
    else:
        for records in zip(*iterators):
            if not dnaio.records_are_mates(*records):
                raise dnaio.FastqFormatError(
                    f"Records are out of sync, names "
                    f"{', '.join(r.name for r in records)} do not match.",
                    line=None
                )
            yield records
    # Check if all iterators are exhausted.
    for iterator in iterators:
        try:
            _ = next(iterator)
            raise dnaio.FastqFormatError("Input files have an unequal number"
                                         " of FASTQ records.", line=None)
        except StopIteration:
            pass


def filter_fastq(input_files: List[str], output_files: List[str],
                 filters: List[Callable[[Tuple[dnaio.SequenceRecord, ...]], bool]],
                 compression_level: int = DEFAULT_COMPRESSION_LEVEL,
                 threads: int = 0):
    """
    Filter FASTQ input files with the filters in filters and write
    the results to the output file.

    :param filters: Functions that filter a tuple of dnaio.sequence records
    :param input_files: FASTQ input filenames. Compressed files are handled
    automatically.
    :param output_files: FASTQ output filenames. Compressed files are handled
    automatically.
    :param compression_level: Compression level for the output files (if
    applicable)
    """
    if len(input_files) != len(output_files):
        raise ValueError("Number of inputs and outputs should be equal.")
    filtered_fastq_records = multiple_files_to_records(input_files, threads=threads)
    for filter_func in filters:
        filtered_fastq_records = filter(filter_func, filtered_fastq_records)
    with contextlib.ExitStack() as output_stack:
        outputs = [output_stack.enter_context(
               xopen.xopen(output_file, threads=threads, mode="wb",
                       compresslevel=compression_level))
               for output_file in output_files]
        # Use faster methods for more common cases before falling back to
        # generic multiple files mode (which is slower).
        if len(outputs) == 1:
            output = outputs[0]
            for record, in filtered_fastq_records:
                output.write(record.fastq_bytes())
        elif len(outputs) == 2:
            output1 = outputs[0]
            output2 = outputs[1]
            for record1, record2 in filtered_fastq_records:
                output1.write(record1.fastq_bytes())
                output2.write(record2.fastq_bytes())
        elif len(outputs) == 3:
            output1 = outputs[0]
            output2 = outputs[1]
            output3 = outputs[2]
            for record1, record2, record3 in filtered_fastq_records:
                output1.write(record1.fastq_bytes())
                output2.write(record2.fastq_bytes())
                output3.write(record3.fastq_bytes())
        else:  # More than 3 files is quite uncommon.
            for records in filtered_fastq_records:
                for record, output in zip(records, outputs):
                    output.write(record.fastq_bytes())


def parse_threshold(value: str, converter=float):
    """Parse a threshold value that can be a single value or comma-separated list.
    
    Args:
        value: String like "30" or "30,25,20" for per-read thresholds
        converter: Function to convert each value (float or int)
    
    Returns:
        Single value or tuple of values
    """
    if ',' in value:
        return tuple(converter(v.strip()) for v in value.split(','))
    return converter(value)


def initiate_logger(verbose: int = 0, quiet: int = 0):
    log_level = logging.INFO - 10 * (verbose - quiet)
    logger = logging.getLogger("fastq-filter")
    logger.setLevel(log_level)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(log_level)
    formatter = logging.Formatter(
        "{asctime}:{levelname}:{name}: {message}",
        datefmt="%m/%d/%Y %I:%M:%S",
        style="{")
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.description = "Filter FASTQ files on various metrics."
    # Inputs: support both positional and -i/--input flags; combine them in main()
    parser.add_argument("positional_inputs",
                    help="Input FASTQ files. Compression format automatically detected. "
                        "Use '-' for stdin. Alternatively, repeat -i/--input for each file.",
                    nargs='*')
    parser.add_argument("-i", "--input",
                    dest="flag_inputs",
                    help="Input FASTQ file (repeat this flag once per file).",
                    action='append')
    # Outputs: require one -o per output to avoid greedy parsing
    parser.add_argument("-o", "--output",
                    help="Output FASTQ file (repeat this flag once per output). "
                        "Compression format is determined by extension. "
                        "Number of outputs must equal number of inputs. Default: stdout.",
                    action='append')
    parser.add_argument("-l", "--min-length", type=str,
                        help="The minimum length for a read. Use comma-separated "
                             "values for per-read thresholds (e.g., '50,60' for "
                             "R1 and R2).")
    parser.add_argument("-L", "--max-length", type=str,
                        help="The maximum length for a read. Use comma-separated "
                             "values for per-read thresholds (e.g., '150,150').")
    parser.add_argument("-e", "--average-error-rate", type=str,
                        help="The minimum average per base error rate. Use comma-separated "
                             "values for per-read thresholds (e.g., '0.01,0.02').")
    parser.add_argument("-q", "--mean-quality", type=str,
                        help="Average quality. Same as the "
                             "'--average-error-rate' option but specified "
                             "with a phred score. I.e '-q 30' is equivalent "
                             "to '-e 0.001'. Use comma-separated values for "
                             "per-read thresholds (e.g., '30,25').")
    parser.add_argument("-Q", "--median-quality", type=str,
                        help="The minimum median phred score. Use comma-separated "
                             "values for per-read thresholds (e.g., '30,25').")
    parser.add_argument("-c", "--compression-level", type=int,
                        default=DEFAULT_COMPRESSION_LEVEL,
                        help=f"Compression level for the output files. "
                             f"Relevant when output files have a .gz "
                             f"extension. Default: {DEFAULT_COMPRESSION_LEVEL}"
                        )
    parser.add_argument("--threads", type=int, default=0,
                        help="Number of threads to use for compression/decompression (passed to xopen). 0 means single-threaded.")
    parser.add_argument("--apply-filters-to", type=str,
                        help="Comma-separated 1-based indices of which reads to filter. "
                             "For example, '1' filters only R1, '1,2' filters both R1 and R2. "
                             "If not specified and per-read thresholds are given, all reads are filtered. "
                             "This allows filtering on one read while keeping pairs synchronized.")
    parser.add_argument("--verbose", action="count", default=0,
                        help="Report stats on individual filters.")
    parser.add_argument("--quiet", action="count", default=0,
                        help="Turn of logging output.")
    return parser


def main():
    args = argument_parser().parse_args()
    # Combine positional and flagged inputs
    inputs: List[str] = []
    if getattr(args, "positional_inputs", None):
        inputs.extend(args.positional_inputs)
    if getattr(args, "flag_inputs", None):
        inputs.extend(args.flag_inputs)
    output = args.output if args.output else ["-"]
    threads = getattr(args, "threads", 0)
    filters = []

    initiate_logger(args.verbose, args.quiet)
    log = logging.getLogger("fastq-filter")
    log.info(f"input files: {', '.join(inputs)}")
    log.info(f"output files: {', '.join(output)}")

    # Parse filter indices (1-based from CLI, converted to 0-based for C)
    filter_indices = None
    if args.apply_filters_to:
        indices_1based = [int(x.strip()) for x in args.apply_filters_to.split(',')]
        filter_indices = tuple(i - 1 for i in indices_1based)  # Convert to 0-based
        log.info(f"Filtering only reads at positions: {', '.join(map(str, indices_1based))}")
    
    # Filters are ordered from low cost to high cost.
    if args.min_length:
        threshold = parse_threshold(args.min_length, int)
        filters.append(MinimumLengthFilter(threshold, filter_indices=filter_indices))
    if args.max_length:
        threshold = parse_threshold(args.max_length, int)
        filters.append(MaximumLengthFilter(threshold, filter_indices=filter_indices))
    if args.average_error_rate:
        threshold = parse_threshold(args.average_error_rate, float)
        filters.append(AverageErrorRateFilter(threshold, filter_indices=filter_indices))
    if args.mean_quality:
        mean_qual = parse_threshold(args.mean_quality, int)
        if isinstance(mean_qual, tuple):
            error_rate = tuple(10 ** -(q / 10) for q in mean_qual)
        else:
            error_rate = 10 ** -(mean_qual / 10)
        filters.append(AverageErrorRateFilter(error_rate, filter_indices=filter_indices))
    if args.median_quality:
        threshold = parse_threshold(args.median_quality, int)
        filters.append(MedianQualityFilter(threshold, filter_indices=filter_indices))
    for filter in filters:
        threshold_str = filter.threshold if not isinstance(filter.threshold, tuple) else ','.join(map(str, filter.threshold))
        log.info(f"{filter.name}: {threshold_str}")
    if not filters:
        log.warning("No filters were applied. Was this intentional?")

    filter_fastq(filters=filters,
                 input_files=inputs,
                 output_files=output,
                 compression_level=args.compression_level,
                 threads=threads)

    if filters:
        total = filters[0].total
        passed = filters[-1].passed
        failed = total - passed
        log.info(f"processed {total} reads.")
        log.info(f"passed: {passed} ({passed * 100 / total :.2f}%)")
        log.info(f"failed: {failed} ({failed * 100 / total :.2f}%)")

    for filter in filters:
        log.debug(f"{filter.name}: "
                  f"{filter.total} processed, {filter.passed} passed "
                  f"({filter.passed * 100 / filter.total :.2f}%)")


if __name__ == "__main__":  # pragma: no cover
    main()
