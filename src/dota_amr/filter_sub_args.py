import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import poisson
import warnings

# =======================================================

warnings.filterwarnings('ignore')

def exp_decay(x, a, b, c):
    """ y = a * e^(-bx) + c"""
    return a * np.exp(-b * x) + c

def run_sub_arg_denoising_pipeline(filepath, alpha, output_file):
    
    df = pd.read_csv(filepath, sep='\t')

    df["Family"] = df['Sub-ARG_Arbitrary_Name'].apply(lambda x: x.split("_")[0])
    #df['Family'] = df['Sub-ARG_Arbitrary_Name'].str.extract(r'^(.*?)_<\d+>$')
    print(df["Family"])
    
    # 2026-08-10: Infer the optional decay-control rate internally and fall back to family-wise decay.
    # Reason: primer panels and small read subsets may not contain enough meca observations.
    decay_control = df[df['Family'].str.casefold() == 'meca']['Cell_count'].sort_values(ascending=False).values
    p_err = decay_control[1] / decay_control[0] if len(decay_control) >= 2 else 0.0
    if len(decay_control) >= 2:
        print(f"[*] Estimated sequencing-noise rate: {p_err:.4%}")
    else:
        print("[*] Using family-wise decay filtering.")

    final_results = []
    
    # 3. 
    for family, group in df.groupby('Family'):
        group = group.sort_values(by='Cell_count', ascending=False)
        counts = group['Cell_count'].values
        args = group['Sub-ARG_Arbitrary_Name'].values
        
        if len(counts) == 0: continue
        
        true_peaks = [args[0]]
        true_counts = [counts[0]]
        
        if len(counts) > 1:
            # poisson cutoff
            expected_noise = counts[0] * p_err
            poisson_cutoff = poisson.ppf(1 - alpha, expected_noise) + 1
            
            # decay
            tail_idx = [i for i, c in enumerate(counts) if c < poisson_cutoff and i > 0]
            if len(tail_idx) < 4:
                tail_idx = list(range(max(1, len(counts)//2), len(counts)))
            
            popt = None
            if len(tail_idx) >= 3:
                x_tail = np.array(tail_idx)
                y_tail = counts[x_tail]
                try:
                    p0 = (max(y_tail), 0.5, min(y_tail))
                    popt, _ = curve_fit(exp_decay, x_tail, y_tail, p0=p0, maxfev=10000)
                    std_res = np.std(y_tail - exp_decay(x_tail, *popt))
                except: popt = None

            # 
            is_noise_tail = False
            for i in range(1, len(counts)):
                obs = counts[i]
                if not is_noise_tail:
                    # noise = failed_possion + decay
                    if poisson.sf(obs - 1, expected_noise) >= alpha:
                        is_noise_tail = True
                        continue
                    
                    is_peak = True
                    if popt is not None:
                        expected_y = exp_decay(i, *popt)
                        # 
                        if obs <= expected_y + max(3 * std_res, expected_y * 0.5):
                            is_peak = False
                            
                    if is_peak:
                        true_peaks.append(args[i])
                        true_counts.append(obs)
                    else:
                        is_noise_tail = True
                else:
                    break
        
        for a, c in zip(true_peaks, true_counts):
            final_results.append({'Sub-ARG_final': a, 'Cell_count': c})

    
    final_df = pd.DataFrame(final_results)
    final_df['Family'] = final_df['Sub-ARG_final'].str.extract(r'^(.*?)(?:_<\d+>)$')
    final_df = final_df.sort_values(by=['Family', 'Cell_count'], ascending=[True, False])

    final_df[['Sub-ARG_final', 'Cell_count']].to_csv(output_file, sep='\t', index=False)

    return final_df['Sub-ARG_final'].to_list()
