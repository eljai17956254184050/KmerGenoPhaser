#!/usr/bin/env python3
"""
mapping_counts_to_blocks.py
增强了列名处理、错误检查、日志
"""
import argparse
import os
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input_dir",      required=True,
                   help="map_kmers_to_genome.py 输出的目录 (含 *_mapping.tsv)")
    p.add_argument("--output_dir",     required=True,
                   help="输出目录 (block .txt 文件)")
    p.add_argument("--dominance_thr",  type=float, default=0.55,
                   help="主导物种阈值 (默认: 0.55)")
    p.add_argument("--min_counts",     type=int,   default=10,
                   help="最小 k-mer 计数 (默认: 10)")
    return p.parse_args()


def parse_tsv_header(header_line):
    """✅ FIXED: 灵活的表头解析"""
    header_line = header_line.lstrip("#").strip()
    parts = header_line.split("\t")
    if len(parts) < 3:
        raise ValueError(f"表头列数不足: {header_line}")
    species = [p.strip() for p in parts[2:]]
    if not species:
        raise ValueError(f"找不到物种列: {header_line}")
    return species


def call_dominant(counts_dict, dominance_thr, min_counts):
    """确定窗口的主导物种"""
    total = sum(counts_dict.values())
    if total < min_counts:
        return "LowInfo", total
    
    best_sp = max(counts_dict, key=counts_dict.get)
    best_frac = counts_dict[best_sp] / total
    
    if best_frac >= dominance_thr:
        return best_sp, total
    else:
        return "Mixed", total


def rle_merge(windows):
    """合并相邻的同物种窗口"""
    if not windows:
        return []
    
    merged = []
    cur_start, cur_end, cur_bl = windows[0]
    
    for start, end, bl in windows[1:]:
        if bl == cur_bl:
            cur_end = end
        else:
            merged.append((cur_start, cur_end, cur_bl))
            cur_start, cur_end, cur_bl = start, end, bl
    
    merged.append((cur_start, cur_end, cur_bl))
    return merged


def convert_file(tsv_path, output_dir, dominance_thr, min_counts):
    """✅ FIXED: 健壮的 TSV 解析和错误处理"""
    species = []
    windows = []
    windows_with_data = 0

    try:
        with open(tsv_path) as fh:
            for line_num, line in enumerate(fh, 1):
                line = line.strip()
                if not line:
                    continue
                
                # 解析表头
                if line.startswith("#") or (line_num == 1):
                    try:
                        species = parse_tsv_header(line)
                        continue
                    except ValueError as e:
                        print(f"  [WARN] {os.path.basename(tsv_path)}:{line_num} — {e}")
                        continue
                
                # 解析数据行
                parts = line.split("\t")
                if len(parts) < 2 + len(species):
                    continue
                
                try:
                    start = int(parts[0])
                    end = int(parts[1])
                    counts = {}
                    for i, sp in enumerate(species):
                        counts[sp] = int(parts[2 + i].strip())
                    
                    bl, total = call_dominant(counts, dominance_thr, min_counts)
                    windows.append((start, end, bl))
                    
                    if total > 0:
                        windows_with_data += 1
                
                except (ValueError, IndexError):
                    continue
    
    except IOError as e:
        print(f"  [ERROR] 无法读取 {tsv_path}: {e}")
        return os.path.basename(tsv_path).replace("_mapping.tsv", ""), 0, 0, 0

    # 验证
    if not windows:
        print(f"  [WARN] 未解析到有效窗口: {os.path.basename(tsv_path)}")
        return os.path.basename(tsv_path).replace("_mapping.tsv", ""), 0, 0, 0

    blocks = rle_merge(windows)
    chrom = os.path.basename(tsv_path).replace("_mapping.tsv", "")
    out_path = os.path.join(output_dir, f"{chrom}.txt")

    # 写输出
    try:
        with open(out_path, "w") as fh:
            fh.write("Start\tEnd\tBloodline\n")
            for start, end, bl in blocks:
                fh.write(f"{start}\t{end}\t{bl}\n")
    except IOError as e:
        print(f"  [ERROR] 写入 {out_path} 失败: {e}")
        return chrom, 0, len(windows), windows_with_data

    info = f"{len(blocks)} blocks from {len(windows)} windows"
    if windows_with_data < len(windows):
        info += f" ({windows_with_data} w/ k-mers)"
    
    print(f"  ✓ {chrom:<30} {info}")
    return chrom, len(blocks), len(windows), windows_with_data


def main():
    args = parse_args()
    
    if not os.path.isdir(args.input_dir):
        print(f"[ERROR] 输入目录不存在: {args.input_dir}")
        sys.exit(1)
    
    os.makedirs(args.output_dir, exist_ok=True)

    tsv_files = sorted([f for f in os.listdir(args.input_dir)
                       if f.endswith("_mapping.tsv")])
    
    if not tsv_files:
        print(f"[ERROR] 找不到 *_mapping.tsv 文件在: {args.input_dir}")
        print(f"        来源: map_kmers_to_genome.py")
        sys.exit(1)

    print(f"✓ 找到 {len(tsv_files)} 个 TSV 文件")
    print(f"  参数: dominance={args.dominance_thr}, min_counts={args.min_counts}\n")
    print("-" * 70)

    total_blocks = 0
    total_windows = 0
    total_windows_hit = 0

    for fname in tsv_files:
        chrom, n_blocks, n_windows, n_hit = convert_file(
            os.path.join(args.input_dir, fname),
            args.output_dir,
            args.dominance_thr,
            args.min_counts
        )
        total_blocks += n_blocks
        total_windows += n_windows
        total_windows_hit += n_hit

    print("-" * 70)
    print(f"\n✓ 完成!")
    print(f"  总blocks: {total_blocks:,}")
    print(f"  总windows: {total_windows:,}")
    print(f"  有数据windows: {total_windows_hit:,}")
    print(f"  输出: {args.output_dir}/\n")

    if total_blocks == 0:
        print("[ERROR] 未生成任何 block！检查输入数据")
        sys.exit(1)


if __name__ == "__main__":
    main()
