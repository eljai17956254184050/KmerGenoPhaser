#!/usr/bin/env python3
"""
map_kmers_to_genome.py
修复了canonical_kmer函数、列名处理、错误检查
"""
import os
import argparse
import time
import multiprocessing
import sys
from Bio import SeqIO


def get_canonical_kmer(kmer):
    """✅ FIXED: 正确的DNA反向互补和规范化"""
    if not kmer:
        return None
    
    kmer = kmer.upper()
    # 正确方式：先做互补，再反向
    complement_table = str.maketrans("ATCG", "TAGC")
    rev_complement = kmer.translate(complement_table)[::-1]
    # 返回字典序较小的那个（规范形式）
    return kmer if kmer <= rev_complement else rev_complement


def init_worker(db_ref):
    global global_kmer_db
    global_kmer_db = db_ref


def process_chromosome(args):
    chrom_name, chrom_seq, k, window_size, output_dir, species_list = args
    seq_len     = len(chrom_seq)
    num_windows = (seq_len + window_size - 1) // window_size
    counts      = {sp: [0] * num_windows for sp in species_list}
    
    # 统计
    window_hits = {sp: 0 for sp in species_list}
    for i in range(seq_len - k + 1):
        kmer = get_canonical_kmer(chrom_seq[i:i + k])
        if kmer and kmer in global_kmer_db:
            sp = global_kmer_db[kmer]
            window_idx = i // window_size
            if window_idx < num_windows:
                counts[sp][window_idx] += 1
                window_hits[sp] += 1

    # 输出
    out_p = os.path.join(output_dir, f"{chrom_name}_mapping.tsv")
    with open(out_p, "w") as f:
        f.write("#Start\tEnd\t" + "\t".join(species_list) + "\n")
        for w in range(num_windows):
            start = w * window_size
            end   = min(start + window_size, seq_len)
            f.write(f"{start}\t{end}\t"
                    + "\t".join(str(counts[sp][w]) for sp in species_list)
                    + "\n")
    
    total_hits = sum(window_hits.values())
    return chrom_name, total_hits, window_hits


def find_column_index(columns, possible_names):
    """✅ FIXED: 灵活的列名查找（不区分大小写）"""
    cols_lower = {h.lower().strip(): i for i, h in enumerate(columns)}
    for name in possible_names:
        if name.lower().strip() in cols_lower:
            return cols_lower[name.lower().strip()]
    return None


def main():
    import pandas as pd

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--merged_kmer_file", required=True,
                        help="从 Step 2 (equalize_and_sample.py) 输出的 k-mer 文件")
    parser.add_argument("--genome_file",      required=True,
                        help="目标基因组 FASTA")
    parser.add_argument("--output_dir",       required=True,
                        help="输出目录，存放 *_mapping.tsv")
    parser.add_argument("--species_list",     required=True,
                        help="物种列表（逗号分隔）")
    parser.add_argument("--k",                type=int, required=True,
                        help="K-mer 大小")
    parser.add_argument("--threads",          type=int, default=10,
                        help="线程数 (默认: 10)")
    parser.add_argument("--window_size",      type=int, default=100000,
                        help="窗口大小，单位bp (默认: 100000)")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # ✅ 第1步：加载 k-mer 库，加入充分的检查
    print("Loading k-mer database...")
    
    if not os.path.isfile(args.merged_kmer_file):
        print(f"[ERROR] K-mer 文件不存在: {args.merged_kmer_file}")
        print(f"        → 检查 Step 2 (equalize_and_sample.py) 是否成功完成")
        print(f"        → 检查文件路径是否正确")
        sys.exit(1)

    file_size = os.path.getsize(args.merged_kmer_file)
    if file_size < 100:
        print(f"[ERROR] K-mer 文件过小 ({file_size} bytes)，可能只有表头")
        print(f"        → 检查 Step 1.6 (merge) 的输出")
        sys.exit(1)

    try:
        df = pd.read_csv(args.merged_kmer_file, sep="\t")
    except Exception as e:
        print(f"[ERROR] 读取 k-mer 文件失败: {e}")
        sys.exit(1)

    # ✅ 灵活的列名检测
    print(f"  检测到列: {list(df.columns)}")
    kmer_col_idx = find_column_index(df.columns, ["Kmer", "K-mer", "kmer"])
    sp_col_idx = find_column_index(df.columns, ["Species", "species", "Sp"])

    if kmer_col_idx is None or sp_col_idx is None:
        print(f"[ERROR] 无法找到必需的列 (Kmer 或 Species)")
        print(f"        可用列: {list(df.columns)}")
        sys.exit(1)

    kmer_col = df.columns[kmer_col_idx]
    sp_col = df.columns[sp_col_idx]
    print(f"  使用列: Kmer='{kmer_col}', Species='{sp_col}'")

    local_db       = {}
    conflict_kmers = set()

    for _, row in df.iterrows():
        try:
            kmer = str(row[kmer_col]).strip().upper()
            sp = str(row[sp_col]).strip()
            canonical = get_canonical_kmer(kmer)
            if not canonical:
                continue
            
            if canonical in local_db and local_db[canonical] != sp:
                conflict_kmers.add(canonical)
            else:
                local_db[canonical] = sp
        except Exception:
            continue

    for c in conflict_kmers:
        if c in local_db:
            del local_db[c]

    print(f"  ✓ 加载了 {len(local_db):,} 个 k-mer")
    print(f"    (移除了 {len(conflict_kmers):,} 个冲突 k-mer)")

    if len(local_db) == 0:
        print("[ERROR] K-mer 库为空！")
        print("        可能原因:")
        print("          1) 所有 k-mer 都是非特异性的（冲突）")
        print("          2) K-mer 文件格式错误")
        sys.exit(1)

    # ✅ 第2步：解析基因组
    print(f"\nParsing genome...")
    species_list = [s.strip() for s in args.species_list.split(",")]
    tasks = []
    for rec in SeqIO.parse(args.genome_file, "fasta"):
        tasks.append((rec.id, str(rec.seq),
                      args.k, args.window_size,
                      args.output_dir, species_list))

    if not tasks:
        print("[ERROR] 基因组 FASTA 中没有序列!")
        sys.exit(1)

    print(f"Mapping {len(tasks)} chromosome(s) | window={args.window_size:,} bp | threads={args.threads}")
    t0 = time.time()

    # ✅ 第3步：并行处理
    with multiprocessing.Pool(processes=args.threads,
                              initializer=init_worker,
                              initargs=(local_db,)) as pool:
        results = pool.map(process_chromosome, tasks)

    elapsed = time.time() - t0

    # ✅ 第4步：总结报告
    total_hits = sum(r[1] for r in results)
    print(f"\n{'Chromosome':<25} {'Hits':>12} Species分布")
    print("-" * 70)
    for chrom_name, hits, sp_hits in results:
        dist = ", ".join([f"{sp}:{sp_hits.get(sp, 0)}" for sp in species_list])
        print(f"{chrom_name:<25} {hits:>12,} {dist}")
    print("-" * 70)
    print(f"{'TOTAL':<25} {total_hits:>12,}\n")
    print(f"✓ 完成 ({elapsed:.1f}s)")
    print(f"输出目录: {args.output_dir}/")
    print(f"下一步: 用 mapping_counts_to_blocks.py 转换 TSV → block .txt")


if __name__ == "__main__":
    main()
