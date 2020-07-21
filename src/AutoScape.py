import networkx as nx
import sys
import csv

if __name__ == "__main__":
    

    arg_len = len(sys.argv)

    if arg_len < 3:
        print("ERROR: not enough arguments. Specify file and centrality measure")
        sys.exit()

    DATA_FILE = sys.argv[1]
    centrality = sys.argv[2]


    edges = []

    with open(DATA_FILE, "r") as dataFile:
        for row in dataFile:
            line = row.strip().split("\t")
            edges.append((line[2], line[3]))
            
    G = nx.Graph()
    G.add_edges_from(edges)

    nx.info(G)

    if centrality == "degree":

        avg_degree = s = sum(dict(G.degree()).values())/G.number_of_nodes()
        #print(avg_degree)



        deg = G.degree()
        #print(deg)

        deg_dict = {}
        for gene in deg:
            deg_dict[gene[0]] = gene[1]

        print(deg_dict)
        # Create graph usong only central genes (deg > avg_degree)


        central_genes = []

        for pair in deg:
            if pair[1] > avg_degree:
                central_genes.append(pair[0])
                
        print(len(central_genes))



        central_edges = []

        for gene1, gene2 in edges:
            if gene1 in central_genes and gene2 in central_genes:
                central_edges.append((gene1, gene2))
                
        print(central_edges)





        central_G = nx.Graph()

        central_G.add_edges_from(central_edges)

        nx.info(central_G)





        central_G.nodes()





        SAVE_PATH = DATA_FILE[:-4] + "_CENTRAL_GENES.txt"
        genes_list = list(central_G.nodes())

        with open(SAVE_PATH, "w+") as saveFile:
            saveFile.write("Gene\tDegree\n")
            for gene in genes_list:
                saveFile.write(gene + "\t" + str(deg_dict[gene]) + "\n")
        





        deg_list = sorted(list(G.degree()), reverse = True, key = lambda x:x[1])

        central_deg_list = []

        for gene in deg_list:
            if gene[0] in genes_list:
                central_deg_list.append(gene)
                
        sorted(central_deg_list, reverse = True, key = lambda x:x[1])

    else:
        print("Centralioty measrue not recognized: choose from [degree, betweenness, TO ADD]")
        sys.exit()
    # To add: take in input interaction file and centrality measure to use. Then, return a file with central genes: GENE_NAME DEG_ORIGINAL DEG_CENTRAL_NETWORK + avg degree in both networks (may be useful). To that for now, then, decide if consider edge weight.
