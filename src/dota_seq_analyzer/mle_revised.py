#!/usr/bin/env python3
import math
from collections import Counter
from typing import List, Dict

def parse_and_analyze_perfect_corrected_revised(
    _16s_packets: List[Dict],
    p_match = 0.90,
    p_none = 0.09,
    p_error = 0.01,
    alpha_prior = 1.0,
    beta_prior = 9.0,
    min_confidence = 0.95,
    min_noise_reads = 2,
    noise_cutoff_ratio = 0.05,
):
    """
    Infer one corrected taxonomy path and contamination estimate per barcode.

    The algorithm works in three broad stages:
    1. Use maximum likelihood to choose the dominant taxon at each rank.
    2. Stop going deeper when the best rank assignment is not confident enough.
    3. Treat small off-path taxa as technical noise and larger off-path taxa as
       real contamination before computing a Bayesian contamination estimate.
    """
    
    # Taxonomy ranks are evaluated from broad to specific. Once confidence fails
    # at a rank, all deeper ranks are set to None.
    levels = ["R1", "P", "C", "O", "F", "G", "S"]

    reads_list = []
    for packet in _16s_packets:
        tax = packet["taxonomy"]
        tax.pop("classifiable")
        reads_list.append(tax)

    total_reads = len(reads_list)

    # ---------------- 1. Stepwise adaptive MLE taxonomy path inference ----------------
    final_path = {}
    last_valid_confidence = 1.0
    is_truncated = False

    for lvl_idx, lvl in enumerate(levels):
        if is_truncated:
            final_path[lvl] = None
            continue

        if lvl_idx == 0:
            valid_reads = reads_list
        else:
            prev_lvl = levels[lvl_idx - 1]
            prev_expected = final_path[prev_lvl]
            valid_reads = [r for r in reads_list if r.get(prev_lvl) == prev_expected]

        # 🛠️ 【微调点 1】：允许 None 成为合法候选者，不再用 is not None 过滤它
        # 这样当大多数 reads 走到某一层断掉变成 None 时，None 作为一个群体能抱团对抗噪声
        lvl_candidates = set(r[lvl] for r in valid_reads)

        # 如果彻底没有有效 reads 分支了，再截断
        if not lvl_candidates or (len(lvl_candidates) == 1 and None in lvl_candidates):
            is_truncated = True
            final_path[lvl] = None
            continue

        mle_scores = {}
        for cand in lvl_candidates:
            log_likelihood = 0.0
            for r in reads_list:
                path_conflict = False
                for p_idx in range(lvl_idx):
                    p_lvl = levels[p_idx]
                    if r.get(p_lvl) is not None and r.get(p_lvl) != final_path[p_lvl]:
                        path_conflict = True
                        break
                if path_conflict:
                    log_likelihood += math.log(p_error)
                    continue

                r_val = r.get(lvl)
                
                # 🛠️ 【微调点 2】：针对 None 候选者，设计合理的条件似然度打分
                if cand is None:
                    if r_val is None:
                        # 候选是 None，Read 也是 None -> 完美匹配，给 p_match
                        log_likelihood += math.log(p_match)
                    else:
                        # 候选是 None，Read 却有具体的菌名 -> 属于错配，给 p_error
                        log_likelihood += math.log(p_error)
                else:
                    # 候选是具体菌名（如 Vibrio），维持原来的传统打分逻辑不变
                    if r_val is None:
                        log_likelihood += math.log(p_none)
                    elif r_val == cand:
                        log_likelihood += math.log(p_match)
                    else:
                        log_likelihood += math.log(p_error)
                        
            mle_scores[cand] = log_likelihood

        # Convert log-likelihoods to relative probabilities with the
        # log-sum-exp trick, avoiding overflow or underflow on many reads.
        max_log = max(mle_scores.values())
        exp_scores = {
            cand: math.exp(score - max_log) for cand, score in mle_scores.items()
        }
        sum_exp = sum(exp_scores.values())

        best_cand = max(mle_scores, key=mle_scores.get)
        confidence = exp_scores[best_cand] / sum_exp

        # If confidence is high enough AND the winning candidate is not None,
        # accept this rank and continue. If None wins or confidence fails, truncate.
        if confidence >= min_confidence and best_cand is not None:
            final_path[lvl] = best_cand
            last_valid_confidence = confidence
        else:
            is_truncated = True
            final_path[lvl] = None

    # ---------------- 2. Separate technical noise (typing errors) from real contamination ----------------
    # Find the current safe truncation level.
    active_lvl_idx = 0
    for i, lvl in enumerate(levels):
        if final_path[lvl] is not None:
            active_lvl_idx = i

    target_lvl = levels[active_lvl_idx]
    target_bug_name = final_path[target_lvl]

    # Count each organism at the inferred terminal level.
    lvl_values = [
        r[target_lvl] for r in reads_list if r.get(target_lvl) is not None
    ]
    counts_summary = Counter(lvl_values)

    match_reads_count = 0
    real_contamination_count = 0
    technical_noise_count = 0

    for bug, count in counts_summary.items():
        if bug == target_bug_name:
            match_reads_count += count
        else:
            bug_ratio = count / total_reads
            if count < min_noise_reads or bug_ratio < noise_cutoff_ratio:
                technical_noise_count += count
            else:
                real_contamination_count += count

    # ---------------- 3. Calculate the corrected Bayesian posterior contamination rate ----------------
    corrected_match_count = match_reads_count + technical_noise_count

    alpha_post = alpha_prior + real_contamination_count
    beta_post = beta_prior + corrected_match_count
    bayesian_contamination_mean = alpha_post / (alpha_post + beta_post)

    taxa_result_str = " | ".join(
        f"{lvl} - {final_path[lvl] if final_path[lvl] else 'None'}"
        for lvl in levels
    )

    return total_reads, technical_noise_count, taxa_result_str, last_valid_confidence, bayesian_contamination_mean
