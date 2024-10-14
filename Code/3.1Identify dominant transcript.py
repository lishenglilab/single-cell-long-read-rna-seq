import multiprocessing as mp
def read_data(file_path):
    with open(file_path, 'r') as file:
        headers = file.readline().strip().split('\t')
        data = [line.strip().split('\t') for line in file]
    return headers, data

def process_gene(data):
    gene_rows, headers = data  
    gene_id = gene_rows[0][1] 
    results = []
    num_samples = len(gene_rows[0]) - 2 

    max_transcripts = [None] * num_samples
    for i in range(num_samples):
        max_usage = -1
        for row in gene_rows:
            current_usage = float(row[i + 2])  
            if current_usage > max_usage:
                max_usage = current_usage
                max_transcripts[i] = row[0] 

    for i in range(num_samples):
        results.append([gene_id, headers[i + 2], max_transcripts[i]])

    return results

def main():
    file_path = 'Result/transcript_usage.txt'
    headers, data = read_data(file_path)

    gene_dict = {}
    for row in data:
        gene_id = row[1]
        if gene_id not in gene_dict:
            gene_dict[gene_id] = []
        gene_dict[gene_id].append(row)

    with mp.Pool(mp.cpu_count() - 10) as pool:
        # 将 headers 作为第二个元素传递
        results = pool.map(process_gene, [(gene_rows, headers) for gene_rows in gene_dict.values()])

    flat_results = [item for sublist in results for item in sublist]

    with open('Result/dominant_transcripts_each_cell.txt', 'w') as f:
        f.write('Gene\tCellid\tDominantTranscript\n')
        for line in flat_results:
            f.write('\t'.join(line) + '\n')

if __name__ == '__main__':
    main()
