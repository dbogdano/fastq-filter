#!/usr/bin/env python3
"""Test filtering with three FASTQ files (R1, R2, R3)."""

import tempfile
import subprocess
from pathlib import Path
import sys

def create_triple_fastq_files(tmpdir):
    """Create three FASTQ files for testing."""
    r1 = Path(tmpdir) / "R1.fastq"
    r2 = Path(tmpdir) / "R2.fastq"
    r3 = Path(tmpdir) / "R3.fastq"
    
    # R1: 100bp, high quality
    r1.write_text(
        "@read1\n" + "A" * 100 + "\n+\n" + "I" * 100 + "\n" +
        "@read2\n" + "T" * 100 + "\n+\n" + "H" * 100 + "\n" +
        "@read3\n" + "G" * 100 + "\n+\n" + "G" * 100 + "\n"
    )
    
    # R2: 80bp, medium quality
    r2.write_text(
        "@read1\n" + "C" * 80 + "\n+\n" + "H" * 80 + "\n" +
        "@read2\n" + "G" * 80 + "\n+\n" + "G" * 80 + "\n" +
        "@read3\n" + "A" * 80 + "\n+\n" + "5" * 80 + "\n"  # Low quality
    )
    
    # R3: 10bp index, variable quality
    r3.write_text(
        "@read1\n" + "ATCGATCGAT" + "\n+\n" + "I" * 10 + "\n" +
        "@read2\n" + "GCTAGCTAGC" + "\n+\n" + "I" * 10 + "\n" +
        "@read3\n" + "TATATATATAT"[:10] + "\n+\n" + "5" * 10 + "\n"  # Low quality
    )
    
    return r1, r2, r3

def run_filter(python_exe, inputs, outputs, args):
    """Run fastq_filter."""
    cmd = [python_exe, "-m", "fastq_filter", *args]
    for out in outputs:
        cmd.extend(["-o", str(out)])
    cmd.extend([str(f) for f in inputs])
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result

def count_reads(filepath):
    """Count reads in FASTQ file."""
    if not filepath.exists():
        return 0
    return filepath.read_text().count("@read")

def main():
    tmpdir = tempfile.mkdtemp(prefix="triple_fastq_test_")
    print(f"Testing three-file filtering in: {tmpdir}\n")
    
    r1_in, r2_in, r3_in = create_triple_fastq_files(tmpdir)
    r1_out = Path(tmpdir) / "R1_out.fastq"
    r2_out = Path(tmpdir) / "R2_out.fastq"
    r3_out = Path(tmpdir) / "R3_out.fastq"
    
    python_exe = sys.executable
    
    print("=" * 70)
    print("INPUT DATA (3 files):")
    print("=" * 70)
    print("R1: 100bp, Qualities: read1=Q40, read2=Q39, read3=Q38")
    print("R2:  80bp, Qualities: read1=Q39, read2=Q38, read3=Q20")
    print("R3:  10bp, Qualities: read1=Q40, read2=Q40, read3=Q20")
    
    # Test 1: Single threshold
    print("\n" + "=" * 70)
    print("TEST 1: Single minimum length threshold")
    print("=" * 70)
    print("Command: -l 15")
    print("Expected: All 3 read sets pass (any read >= 15)")
    
    result = run_filter(
        python_exe,
        [r1_in, r2_in, r3_in],
        [r1_out, r2_out, r3_out],
        ["-l", "15"]
    )
    
    if result.returncode == 0:
        count = count_reads(r1_out)
        print(f"✓ Result: {count}/3 read sets passed")
        assert count == 3, f"Expected 3, got {count}"
    else:
        print(f"✗ Failed: {result.stderr}")
        return
    
    # Test 2: Per-read thresholds (3 values)
    print("\n" + "=" * 70)
    print("TEST 2: Per-read minimum lengths")
    print("=" * 70)
    print("Command: -l 95,75,8")
    print("Expected: All 3 pass (R1>=95✓, R2>=75✓, R3>=8✓)")
    
    result = run_filter(
        python_exe,
        [r1_in, r2_in, r3_in],
        [r1_out, r2_out, r3_out],
        ["-l", "95,75,8"]
    )
    
    if result.returncode == 0:
        count = count_reads(r1_out)
        print(f"✓ Result: {count}/3 read sets passed")
        assert count == 3, f"Expected 3, got {count}"
    else:
        print(f"✗ Failed: {result.stderr}")
        return
    
    # Test 3: Stricter per-read thresholds
    print("\n" + "=" * 70)
    print("TEST 3: Stricter per-read thresholds")
    print("=" * 70)
    print("Command: -l 95,75,8 -Q 39,38,35")
    print("Expected: 2 pass (read3 fails: R2 Q20<38 and R3 Q20<35)")
    
    result = run_filter(
        python_exe,
        [r1_in, r2_in, r3_in],
        [r1_out, r2_out, r3_out],
        ["-l", "95,75,8", "-Q", "39,38,35"]
    )
    
    if result.returncode == 0:
        count = count_reads(r1_out)
        print(f"✓ Result: {count}/3 read sets passed")
        assert count == 2, f"Expected 2, got {count}"
        
        # Show which reads passed
        if r1_out.exists():
            content = r1_out.read_text()
            passed_reads = [line[1:] for line in content.split('\n') if line.startswith('@')]
            print(f"  Passed reads: {', '.join(passed_reads)}")
    else:
        print(f"✗ Failed: {result.stderr}")
        return
    
    # Test 4: Per-read quality only
    print("\n" + "=" * 70)
    print("TEST 4: Per-read median quality thresholds")
    print("=" * 70)
    print("Command: -Q 38,30,30")
    print("Expected: 2 pass (read3 fails: R2 Q20<30)")
    
    result = run_filter(
        python_exe,
        [r1_in, r2_in, r3_in],
        [r1_out, r2_out, r3_out],
        ["-Q", "38,30,30"]
    )
    
    if result.returncode == 0:
        count = count_reads(r1_out)
        print(f"✓ Result: {count}/3 read sets passed")
        assert count == 2, f"Expected 2, got {count}"
    else:
        print(f"✗ Failed: {result.stderr}")
        return
    
    print("\n" + "=" * 70)
    print("✓ ALL TESTS PASSED - Three-file filtering works perfectly!")
    print("=" * 70)
    print(f"\nCleanup: rm -rf {tmpdir}")

if __name__ == "__main__":
    main()
