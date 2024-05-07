import os
import re
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests

# 用于提取似然值的正则表达式
likelihood_re = re.compile(r"lnL\(ntime:\s*\d+\s*np:\s*\d+\):\s*([-.\d]+)")

def extract_likelihood(file_path):
    """从mlc文件中提取最大似然值。"""
    with open(file_path, 'r') as file:
        for line in file:
            match = likelihood_re.search(line)
            if match:
                print(f"提取到{file_path}的似然值: {match.group(1)}")
                return float(match.group(1))
    print(f"未能在{file_path}中提取似然值。")
    return None

def calculate_LRT(null_likelihood, alt_likelihood):
    """计算似然比检验（LRT）的值。"""
    lrt_value = 2 * (alt_likelihood - null_likelihood)
    print(f"LRT值计算结果: {lrt_value}")
    return lrt_value

def apply_fdr_correction(p_values, alpha=0.01):
    """对P值应用FDR校正。"""
    _, p_adjusted, _, _ = multipletests(p_values, alpha=alpha, method='fdr_bh')
    return p_adjusted

def process_directory(directory_path):
    """处理指定目录，比较null和alternative模型的似然值，并进行LRT计算。"""
    files = os.listdir(directory_path)
    genes = set(filename.split('_')[0] for filename in files)
    
    lrt_values = []
    p_values = []
    gene_names = []
    for gene in genes:
        null_file = os.path.join(directory_path, f"{gene}_mlc00")
        alt_file = os.path.join(directory_path, f"{gene}_mlc20")
        print(f"正在处理基因{gene}...")
        if os.path.exists(null_file) and os.path.exists(alt_file):
            null_likelihood = extract_likelihood(null_file)
            alt_likelihood = extract_likelihood(alt_file)
            if null_likelihood is not None and alt_likelihood is not None:
                lrt_value = calculate_LRT(null_likelihood, alt_likelihood)
                p_value = chi2.sf(lrt_value, df=1)
                print(f"基因{gene}的P值: {p_value}")
                lrt_values.append(lrt_value)
                p_values.append(p_value)
                gene_names.append(gene)
        else:
            print(f"未找到基因{gene}的null或alternative模型文件。")
    
    return gene_names, lrt_values, p_values

# 主程序
directory_path = '/home/jiale/Bioinformatics/Project/06.Ka_Ks/20240228/Results/mlc'
output_file_path = '/home/jiale/Bioinformatics/Project/06.Ka_Ks/20240228/Results/significant_genes.csv'
gene_names, lrt_values, p_values = process_directory(directory_path)

if len(p_values) > 0:
    p_adjusted = apply_fdr_correction(p_values)
    # 保存显著的结果到文件
    with open(output_file_path, 'w') as outfile:
        outfile.write("Gene,LRT,P-value,Adjusted P-value,Significant\n")
        for gene, lrt, p_value, p_adj in zip(gene_names, lrt_values, p_values, p_adjusted):
            significant = "Yes" if p_adj < 0.05 else "No"
            outfile.write(f"{gene},{lrt},{p_value},{p_adj},{significant}\n")
    print("显著性结果已保存到文件。")
else:
    print("没有找到有效的P值进行FDR校正。")
