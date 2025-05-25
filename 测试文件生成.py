import os
import random

# 创建目录结构
os.makedirs('./diseaselist-deg2function', exist_ok=True)
os.makedirs('./druglist-pos1', exist_ok=True)
os.makedirs('./network', exist_ok=True)
os.makedirs('./output', exist_ok=True)


# 随机生成基因名
def generate_gene_names(n):
    return [f"Gene{chr(65 + i // 26)}{chr(65 + i % 26)}" for i in range(n)]


# 生成疾病文件
def generate_disease_file(filename, gene_names, up_ratio=0.5):
    num_genes = len(gene_names)
    up_genes = random.sample(gene_names, int(num_genes * up_ratio))
    with open(filename, 'w') as f:
        for gene in gene_names:
            logfc = round(random.uniform(0.5, 2), 2) if gene in up_genes else round(random.uniform(-2, -0.5), 2)
            pv = round(random.uniform(0.01, 0.05), 2)
            f.write(f"{gene}\t{logfc}\t{pv}\n")


# 生成药物文件
def generate_drug_file(filename, gene_names, pos_ratio=0.5):
    num_genes = len(gene_names)
    pos_genes = random.sample(gene_names, int(num_genes * pos_ratio))
    with open(filename, 'w') as f:
        for gene in gene_names:
            corr = round(random.uniform(0.5, 2), 2) if gene in pos_genes else round(random.uniform(-2, -0.5), 2)
            f.write(f"{gene}\t{corr}\n")


# 生成PPI网络文件
def generate_ppi_network(filename, gene_names, avg_degree=3):
    num_genes = len(gene_names)
    edges = set()
    while len(edges) < num_genes * avg_degree // 2:
        gene1, gene2 = random.sample(gene_names, 2)
        if gene1 != gene2:
            edges.add((gene1, gene2))
    with open(filename, 'w') as f:
        for gene1, gene2 in edges:
            score = random.randint(400, 1000)
            f.write(f"10090.{gene1}\t10090.{gene2}\t{score}\n")


# 生成基因符号文件
def generate_ensp2symbol_file(filename, gene_names):
    with open(filename, 'w') as f:
        for gene in gene_names:
            f.write(f"10090.{gene}\tENST{gene}\tENSG{gene}\t{gene}\n")


# 主函数生成所有文件
def main():
    num_genes = 50
    gene_names = generate_gene_names(num_genes)

    generate_disease_file('./diseaselist-deg2function/disease1.txt', gene_names)
    generate_drug_file('./druglist-pos1/drug1.txt', gene_names)
    generate_ppi_network('./network/testnet1.txt', gene_names)
    generate_ensp2symbol_file('./network/testensp2symbol.txt', gene_names)


if __name__ == "__main__":
    main()
