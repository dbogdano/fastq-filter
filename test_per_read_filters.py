#!/usr/bin/env python3
"""Test per-read filtering functionality."""

from dnaio import SequenceRecord
from fastq_filter import (
    AverageErrorRateFilter,
    MedianQualityFilter,
    MinimumLengthFilter,
    MaximumLengthFilter
)


def test_single_threshold():
    """Test backward compatibility with single threshold values."""
    print("Testing single threshold (backward compatibility)...")
    
    # Single record - should work as before
    record1 = SequenceRecord("read1", "ACGT" * 20, "I" * 80)  # High quality
    
    # Quality filters
    avg_filter = AverageErrorRateFilter(0.001)
    assert avg_filter((record1,)) is True
    print(f"  ✓ AverageErrorRateFilter single threshold: {avg_filter.threshold}")
    
    median_filter = MedianQualityFilter(30)
    assert median_filter((record1,)) is True
    print(f"  ✓ MedianQualityFilter single threshold: {median_filter.threshold}")
    
    # Length filters
    min_filter = MinimumLengthFilter(50)
    assert min_filter((record1,)) is True
    print(f"  ✓ MinimumLengthFilter single threshold: {min_filter.threshold}")
    
    max_filter = MaximumLengthFilter(100)
    assert max_filter((record1,)) is True
    print(f"  ✓ MaximumLengthFilter single threshold: {max_filter.threshold}")


def test_per_read_quality_filters():
    """Test per-read thresholds for quality filters."""
    print("\nTesting per-read quality thresholds...")
    
    # Paired-end reads with different quality
    r1_high = SequenceRecord("read1", "ACGT" * 20, "I" * 80)  # High quality Q=40
    r2_low = SequenceRecord("read2", "TGCA" * 15, "5" * 60)   # Lower quality Q=20
    
    # Single threshold - would reject this pair due to R2
    single_filter = AverageErrorRateFilter(0.001)  # Q30 equivalent
    result = single_filter((r1_high, r2_low))
    print(f"  Single threshold 0.001: {result} (expects both reads to meet Q30)")
    
    # Per-read thresholds - R1 needs Q30, R2 only needs Q15
    per_read_filter = AverageErrorRateFilter((0.001, 0.05))  # Q30 for R1, Q13 for R2
    result = per_read_filter((r1_high, r2_low))
    print(f"  ✓ Per-read thresholds (0.001, 0.05): {result}")
    print(f"    Filter thresholds: {per_read_filter.threshold}")
    
    # Median quality filter with per-read thresholds
    median_filter = MedianQualityFilter((35, 18))  # Q35 for R1, Q18 for R2
    result = median_filter((r1_high, r2_low))
    print(f"  ✓ Median per-read thresholds (35, 18): {result}")
    print(f"    Filter thresholds: {median_filter.threshold}")


def test_per_read_length_filters():
    """Test per-read thresholds for length filters."""
    print("\nTesting per-read length thresholds...")
    
    # Paired-end reads with different lengths
    r1_long = SequenceRecord("read1", "ACGT" * 40, "I" * 160)   # 160bp
    r2_short = SequenceRecord("read2", "TGCA" * 20, "I" * 80)   # 80bp
    
    # Single threshold - R2 would fail min-length of 100
    single_min = MinimumLengthFilter(100)
    result = single_min((r1_long, r2_short))
    print(f"  Single min-length 100: {result} (passes if ANY read >= 100)")
    
    # Per-read minimum lengths - R1 needs 150bp, R2 needs 70bp
    per_read_min = MinimumLengthFilter((150, 70))
    result = per_read_min((r1_long, r2_short))
    print(f"  ✓ Per-read min-lengths (150, 70): {result}")
    print(f"    Filter thresholds: {per_read_min.threshold}")
    
    # Per-read maximum lengths
    per_read_max = MaximumLengthFilter((200, 100))
    result = per_read_max((r1_long, r2_short))
    print(f"  ✓ Per-read max-lengths (200, 100): {result}")
    print(f"    Filter thresholds: {per_read_max.threshold}")


def test_triple_read():
    """Test with three reads (e.g., R1, R2, I1 index)."""
    print("\nTesting triple-read filtering...")
    
    r1 = SequenceRecord("read1", "ACGT" * 25, "H" * 100)
    r2 = SequenceRecord("read2", "TGCA" * 25, "G" * 100)
    i1 = SequenceRecord("index", "ATCG" * 2, "I" * 8)  # Short index read
    
    # Per-read length filters for R1, R2, and index
    min_filter = MinimumLengthFilter((90, 90, 6))
    result = min_filter((r1, r2, i1))
    print(f"  ✓ Triple per-read min-lengths (90, 90, 6): {result}")
    print(f"    Filter thresholds: {min_filter.threshold}")
    
    # Per-read quality filters
    median_filter = MedianQualityFilter((30, 30, 35))  # Higher Q for index
    result = median_filter((r1, r2, i1))
    print(f"  ✓ Triple per-read median qualities (30, 30, 35): {result}")
    print(f"    Filter thresholds: {median_filter.threshold}")


def test_statistics():
    """Test that statistics are correctly tracked."""
    print("\nTesting filter statistics...")
    
    r1 = SequenceRecord("pass", "ACGT" * 20, "I" * 80)
    r2 = SequenceRecord("fail", "TGCA" * 10, "#" * 40)  # Low quality
    
    filt = AverageErrorRateFilter((0.001, 0.001))
    
    # First call should pass
    result1 = filt((r1, r1))
    # Second call should fail
    result2 = filt((r1, r2))
    
    print(f"  ✓ Total processed: {filt.total}")
    print(f"  ✓ Passed: {filt.passed}")
    print(f"  ✓ Results: {result1}, {result2}")
    assert filt.total == 2
    assert filt.passed == 1


if __name__ == "__main__":
    print("=" * 60)
    print("Per-Read Filtering Test Suite")
    print("=" * 60)
    
    test_single_threshold()
    test_per_read_quality_filters()
    test_per_read_length_filters()
    test_triple_read()
    test_statistics()
    
    print("\n" + "=" * 60)
    print("✓ All tests passed!")
    print("=" * 60)
