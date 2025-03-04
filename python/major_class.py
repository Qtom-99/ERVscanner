import pandas as pd
import sys

def process_erv_data(input_file, output_file1, output_file2):
    # 入力ファイルをデータフレームとして読み込む
    df = pd.read_csv(input_file, sep='\t', header=None, names=['Position', 'Name', 'Class'])

    # ポジションごとにクラスのカウントを集計
    pos_class_count = df.groupby(['Position', 'Class']).size().reset_index(name='Count')

    # ポジションごとのクラスの割合を計算
    pos_total_count = df.groupby('Position').size().reset_index(name='Total')
    pos_class_count = pos_class_count.merge(pos_total_count, on='Position')
    pos_class_count['Ratio'] = pos_class_count['Count'] / pos_class_count['Total']

    # ポジションごとの最もメジャーなクラスを見つける
    major_class = pos_class_count.loc[pos_class_count.groupby('Position')['Count'].idxmax()]
    major_class = major_class[['Position', 'Class', 'Ratio']]

    # 出力ファイルに結果を書き込む
    pos_class_count.to_csv(output_file1, sep='\t', index=False)
    major_class.to_csv(output_file2, sep='\t', index=False, header=False)

if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file1 = sys.argv[2]
    output_file2 = sys.argv[3]
    
    process_erv_data(input_file, output_file1, output_file2)
