import sys
from collections import defaultdict

def process_file(input_file, output_filtered, output_unfiltered):
    posid_counts = defaultdict(lambda: {'plus_16_or_minus_0': 0, 'plus_0_or_minus_16': 0, 'total': 0})
    
    # ファイルの読み込みとカウント
    with open(input_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) != 4:
                continue
            posid, _, direction, value = parts
            
            posid_counts[posid]['total'] += 1
            
            if (direction == '+' and value == '16') or (direction == '-' and value == '0'):
                posid_counts[posid]['plus_16_or_minus_0'] += 1
            elif (direction == '+' and value == '0') or (direction == '-' and value == '16'):
                posid_counts[posid]['plus_0_or_minus_16'] += 1
    
    # 80%を超えるPOSIDをフィルタリング
    filtered = []
    unfiltered = []
    
    for posid, counts in posid_counts.items():
        total = counts['total']
        plus_16_or_minus_0_ratio = counts['plus_16_or_minus_0'] / total
        plus_0_or_minus_16_ratio = counts['plus_0_or_minus_16'] / total
        
        if plus_16_or_minus_0_ratio > 0.7:
            filtered.append(f"{posid}\t+")
        elif plus_0_or_minus_16_ratio > 0.7:
            filtered.append(f"{posid}\t-")
        else:
            unfiltered.append(posid)
    
    # 結果をファイルに書き込む
    with open(output_filtered, 'w') as f:
        f.write('\n'.join(filtered) + '\n')
    
    with open(output_unfiltered, 'w') as f:
        f.write('\n'.join(unfiltered) + '\n')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_file> <output_filtered> <output_unfiltered>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_filtered = sys.argv[2]
    output_unfiltered = sys.argv[3]
    
    process_file(input_file, output_filtered, output_unfiltered)
